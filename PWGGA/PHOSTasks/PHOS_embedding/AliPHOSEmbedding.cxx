// Analysis task to perform embedding of the simulated AOD
// into real data (ESD)
// Output conntainer contain AOD with
//   0. Data before embedding
//   1. Data+with embedded signal
//   2. MC information from signal
// Data for embedding are obtained using standard Manager
// Signal is read from TChain with aodTree
//
// Authors: Dmitri Peressounko
// Date   : 28.05.2011

#include "TChain.h"
#include "TFile.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TH2.h"
#include "TGeoManager.h"
#include "TGeoGlobalMagField.h"

#include "AliGeomManager.h"
#include "AliGRPObject.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliCDBManager.h"
#include "AliAODHandler.h"
#include "AliMagF.h"

#include "AliPHOSReconstructor.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSDigit.h"
#include "AliPHOSGeometry.h"

#include "AliPHOSEmbedding.h"

ClassImp(AliPHOSEmbedding)


//________________________________________________________________________
AliPHOSEmbedding::AliPHOSEmbedding(const char *name) 
  : AliAnalysisTaskSE(name),
    fAODChain(0x0),
    fDigitsTree(0x0),
    fClustersTree(0x0),
    fTreeOut(0x0),
    fDigitsArr(0x0),
    fEmbeddedClusters(0x0),
    fEmbeddedCells(0x0),
    fCellsPHOS(0x0),
    fClusterizer(0x0),
    fPHOSReconstructor(0x0),
    fNSignal(0),
    fNCaloClustersOld(0),
    fInitialized(0)   
{
  // Constructor
  for(int i=0;i<5;i++)fOldPHOSCalibration[i]=0x0 ;

}

//________________________________________________________________________
void AliPHOSEmbedding::UserCreateOutputObjects()
{
 //prepare output containers
  
  //Create branch for output MC particles  
  TClonesArray *mcparticles = new TClonesArray("AliAODMCParticle",0);
  mcparticles->SetName(AliAODMCParticle::StdBranchName());
  AddAODBranch("TClonesArray", &mcparticles);

  //Branch with clusters after embedding
  fEmbeddedClusters = new TClonesArray("AliAODCaloCluster",0);
  fEmbeddedClusters->SetName("EmbeddedCaloClusters");
  AddAODBranch("TClonesArray", &fEmbeddedClusters);
  

  //Branch with cells after embedding
  fEmbeddedCells = new AliAODCaloCells();
  fEmbeddedCells->SetName("EmbeddedPHOScells");
  AddAODBranch("AliAODCaloCells", &fEmbeddedCells);

  
  fDigitsArr = new TClonesArray("AliPHOSDigit",3*56*64) ;

  AliAODHandler* handler = (AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
  if (handler) {
      fTreeOut = handler->GetTree();
  }
  else {
       AliWarning("No output AOD Event Handler connected.") ;
  }

}

//________________________________________________________________________
void AliPHOSEmbedding::UserExec(Option_t *) {
  // Main loop, called for each event
  // Perform embedding
  
  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!event) {
    Printf("ERROR: Could not retrieve event");
    PostData(0, fTreeOut);
    return;
  }

  AliCentrality *centrality = event->GetCentrality(); 
  Float_t acentrality=centrality->GetCentralityPercentileUnchecked("V0M");
  if( acentrality <= 0. || acentrality>80.){
    return;
  }

  
  //Read next AOD event
  //If necesary method checks if AOD event is good to embed
  //e.g. there are PHOS clusters etc.
  AliAODEvent * signal = GetNextSignalEvent() ;
  if (!signal) {
    Printf("ERROR: Could not retrieve signal");
    PostData(0, fTreeOut);
    return;
  }
  
  Init() ;

  //Remember PHOS digits before any modification
  if(fCellsPHOS)
    delete fCellsPHOS ;
  fCellsPHOS = new AliESDCaloCells(*(event->GetPHOSCells())) ;

  //First re-reconstruct existing digits to ensure identical reconstruction before and 
  //after embeding
  AliAODEvent * nullSignal=0x0 ;
  MakeEmbedding(event,nullSignal) ;

  //Convert ESD + MC information from signal to AOD
  ConvertESDtoAOD(event);
  
  //Now perform real embedding  
  //Embed signal clusters into event
  MakeEmbedding(event,signal) ;

  //Convert ESD to AOD
  ConvertEmbeddedClusters(event) ;
 
  //Fill MC information
  ConvertMCParticles(signal) ;

  delete signal ;

  PostData(0, fTreeOut);
}
//______________________________________________________________________________
void AliPHOSEmbedding::Init(){
  if(fInitialized)
    return ;
  
  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!event){
    AliError("Can not obtain InputEvent!") ;
    return ;
  }
  Int_t runNum = event->GetRunNumber();
  AliCDBManager::Instance()->SetRun(runNum);

  fPHOSReconstructor = new AliPHOSReconstructor() ;
  AliCDBPath path("PHOS","Calib","RecoParam");
  AliCDBEntry *entry=AliCDBManager::Instance()->Get(path.GetPath());
  if(!entry){
    AliError(Form("Can not get OCDB entry %s",path.GetPath().Data())) ;
    return ;
  }
  
  TObjArray* recoParamArray = (TObjArray*)entry->GetObject();
  AliPHOSRecoParam* recoParam = (AliPHOSRecoParam*)recoParamArray->At(2); 
  fPHOSReconstructor->SetRecoParam(recoParam) ;
  
  InitMF() ;  
  InitGeometry() ;

  fInitialized=kTRUE ;
}

//________________________________________________________________________
AliAODEvent* AliPHOSEmbedding::GetNextSignalEvent(){
  //Read AOD event from the chain
  
  if(fAODChain==0){
    AliError("No chain to read signal events") ;
    return 0 ;
  }
  
  AliAODEvent* aod = new AliAODEvent;
  aod->ReadFromTree(fAODChain);
  
  if(fNSignal>=fAODChain->GetEntries()){
    delete aod ;
    return 0 ;
  }
  fAODChain->GetEvent(fNSignal) ;
  fNSignal++ ;

  //check if event is read and there is at least 2 PHOS clusters
  //Otherwise fill histogram and read next event
  //TODO

  return aod ;
  
}

//________________________________________________________________________
void AliPHOSEmbedding::MakeEmbedding(AliESDEvent *event,AliAODEvent * signal){
  //Perform embedding of the signal to the event

  gROOT->cd(); //make sure that the digits and RecPoints Trees are memory resident
  if(fDigitsTree){
    delete fDigitsTree ;
  }
  fDigitsTree = new TTree("digitstree","digitstree");
  fDigitsTree->Branch("PHOS","TClonesArray", &fDigitsArr, 32000);
  
  if(fClustersTree){
    delete fClustersTree ;
  }
  fClustersTree = new TTree("clustertree","clustertree");

  //Remember number of Clusters before we added new ones...
  fNCaloClustersOld=event->GetNumberOfCaloClusters() ;


/*
  printf("---------------Before---------\n") ;
  Int_t n=event->GetNumberOfCaloClusters() ;
  for (Int_t iClust=0; iClust<event->GetNumberOfCaloClusters(); ++iClust) {
    AliESDCaloCluster * cluster = event->GetCaloCluster(iClust);
     if(!cluster->IsPHOS())
      continue;
   Float_t pos[3] = { 0.};
    cluster->GetPosition(pos);
    printf("Cluster(%d): E=%f, x=%f, z=%f, CPV=%f \n",iClust,cluster->E(),pos[0],pos[2],cluster->GetEmcCpvDistance()) ;
    UShort_t * index    = cluster->GetCellsAbsId() ;
    for(Int_t ic=0; ic < cluster->GetNCells(); ic++ )
    printf("Dig(%d)=%d, ",ic,index[ic]) ;
  printf("\n") ;
	
	
  }
  printf("---------------END Before---------\n") ;
*/

  
  //create digits
  MakeDigits(signal);	  
  
  //clusterize and make tracking
  fPHOSReconstructor->Reconstruct(fDigitsTree,fClustersTree) ;

  //Note that the cirrent ESDEvent will be modified!
  fPHOSReconstructor->FillESD(fDigitsTree, fClustersTree, event) ;


/*
   printf("---------------After---------\n") ;
  for (Int_t iClust=n; iClust<event->GetNumberOfCaloClusters(); ++iClust) {
    AliESDCaloCluster * cluster = event->GetCaloCluster(iClust);
    if(!cluster->IsPHOS())
      continue;
    Float_t pos[3] = { 0.};
    cluster->GetPosition(pos);
    printf("Cluster(%d): E=%f, x=%f, z=%f, CPV=%f, Label=%d \n",iClust-n,cluster->E(),pos[0],pos[2],cluster->GetEmcCpvDistance(),cluster->GetLabel()) ;
    UShort_t * index    = cluster->GetCellsAbsId() ;
    for(Int_t ic=0; ic < cluster->GetNCells(); ic++ )
      printf("Dig(%d)=%d, ",ic,index[ic]) ;
  printf("\n") ;
  }
  printf("---------------END After---------\n") ;
*/


}
//______________________________________________________________________________
void AliPHOSEmbedding::ConvertESDtoAOD(AliESDEvent* esd) 
{
  // ESD Filter analysis task executed for each event
  
  if(!esd)return;
    
  AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);
  AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillExtension(kTRUE);
      
  ConvertHeader(*esd);
  
  // Fetch Stack for debuggging if available 
  
  // Update the header
  Int_t nTracks    = esd->GetNumberOfTracks();
  Int_t nV0s      = esd->GetNumberOfV0s();
  Int_t nCascades = esd->GetNumberOfCascades();
  Int_t nKinks    = esd->GetNumberOfKinks();
  Int_t nVertices = nV0s + nCascades /*V0 wihtin cascade already counted*/+ nKinks + 1 /* = prim. vtx*/;
  Int_t nPileSPDVertices=1+esd->GetNumberOfPileupVerticesSPD(); // also SPD main vertex
  Int_t nPileTrkVertices=esd->GetNumberOfPileupVerticesTracks();
  nVertices+=nPileSPDVertices;
  nVertices+=nPileTrkVertices;
  Int_t nJets     = 0;
  Int_t nCaloClus = esd->GetNumberOfCaloClusters();
  Int_t nFmdClus  = 0;
  Int_t nPmdClus  = esd->GetNumberOfPmdTracks();
    
       
  AODEvent()->ResetStd(nTracks, nVertices, nV0s, nCascades, nJets, nCaloClus, nFmdClus, nPmdClus);
      
  ConvertPrimaryVertices(*esd);
  
  ConvertCaloClusters(*esd);
  
  ConvertEMCALCells(*esd);
  
  ConvertPHOSCells(*esd);
  
}
//______________________________________________________________________________
void AliPHOSEmbedding::ConvertPHOSCells(const AliESDEvent& esd)
{
  // fill PHOS cell info
  if (esd.GetPHOSCells()) { // protection against missing ESD information
    AliESDCaloCells &esdPHcells = *(esd.GetPHOSCells());
    Int_t nPHcell = esdPHcells.GetNumberOfCells() ;
    
    AliAODCaloCells &aodPHcells = *(AODEvent()->GetPHOSCells());
    aodPHcells.CreateContainer(nPHcell);
    aodPHcells.SetType(AliAODCaloCells::kPHOSCell);
    for (Int_t iCell = 0; iCell < nPHcell; iCell++) {      
      aodPHcells.SetCell(iCell,esdPHcells.GetCellNumber(iCell),esdPHcells.GetAmplitude(iCell));
    }
    aodPHcells.Sort();
  }
}

//______________________________________________________________________________
void AliPHOSEmbedding::ConvertEMCALCells(const AliESDEvent& esd)
{
  // fill EMCAL cell info
  if (esd.GetEMCALCells()) { // protection against missing ESD information
    AliESDCaloCells &esdEMcells = *(esd.GetEMCALCells());
    Int_t nEMcell = esdEMcells.GetNumberOfCells() ;
    
    AliAODCaloCells &aodEMcells = *(AODEvent()->GetEMCALCells());
    aodEMcells.CreateContainer(nEMcell);
    aodEMcells.SetType(AliAODCaloCells::kEMCALCell);
    for (Int_t iCell = 0; iCell < nEMcell; iCell++) {      
      aodEMcells.SetCell(iCell,esdEMcells.GetCellNumber(iCell),esdEMcells.GetAmplitude(iCell));
    }
    aodEMcells.Sort();
  }
}
//______________________________________________________________________________
void AliPHOSEmbedding::ConvertCaloClusters(const AliESDEvent& esd)
{
  
  // Access to the AOD container of clusters
  TClonesArray &caloClusters = *(AODEvent()->GetCaloClusters());
  Int_t jClusters(0);
  
  for (Int_t iClust=fNCaloClustersOld; iClust<esd.GetNumberOfCaloClusters(); ++iClust) {
    
    AliESDCaloCluster * cluster = esd.GetCaloCluster(iClust);
    
    Int_t  id        = cluster->GetID();
    Int_t  nLabel    = cluster->GetNLabels();
    Int_t *labels    = cluster->GetLabels();
    
    Float_t energy = cluster->E();
    Float_t posF[3] = { 0.};
    cluster->GetPosition(posF);
    
    AliAODCaloCluster *caloCluster = new(caloClusters[jClusters++]) AliAODCaloCluster(id,
                                                                                      nLabel,
                                                                                      labels,
                                                                                      energy,
                                                                                      posF,
                                                                                      NULL,
                                                                                      cluster->GetType(),0);
    
    Double_t dx=cluster->GetTrackDx() ;
    Double_t dz=cluster->GetTrackDz() ;
    Float_t cpv=99. ; //No track matched by default
    TArrayI * itracks = cluster->GetTracksMatched() ;  
    if(itracks->GetSize()>0){
      Int_t iTr = itracks->At(0);
      if(iTr>=0 && iTr<esd.GetNumberOfTracks()){
        AliESDtrack *track = esd.GetTrack(iTr) ;
        Double_t pt = track->Pt() ;
        Short_t charge = track->Charge() ;
        cpv=TestCPV(dx, dz, pt,charge) ;
      }
    }
    
    caloCluster->SetCaloCluster(cluster->GetDistanceToBadChannel(),
                                cluster->GetDispersion(),
                                cluster->GetM20(), cluster->GetM02(),
                                cpv,  
                                cluster->GetNExMax(),cluster->GetTOF()) ;
    
    caloCluster->SetPIDFromESD(cluster->GetPID());
    caloCluster->SetNCells(cluster->GetNCells());
    caloCluster->SetCellsAbsId(cluster->GetCellsAbsId());
    caloCluster->SetCellsAmplitudeFraction(cluster->GetCellsAmplitudeFraction());

    /*    
    TArrayI* matchedT = 	cluster->GetTracksMatched();
    if (fNumberOfTracks>0 && matchedT && cluster->GetTrackMatchedIndex() >= 0) {	
      for (Int_t im = 0; im < matchedT->GetSize(); im++) {
        Int_t iESDtrack = matchedT->At(im);;
        if (fAODTrackRefs->At(iESDtrack) != 0) {
          caloCluster->AddTrackMatched((AliAODTrack*)fAODTrackRefs->At(iESDtrack));
        }
      }
    }
    */
    
  } 
  caloClusters.Expand(jClusters); // resize TObjArray to 'remove' slots for pseudo clusters	 
}
//______________________________________________________________________________
void AliPHOSEmbedding::ConvertEmbeddedClusters(const AliESDEvent* esd)
{
   //Copy PHOS clusters and cells after embedding
  
  // Access to the AOD container of clusters
  Int_t jClusters(0);
  fEmbeddedClusters->Clear();
  fEmbeddedClusters->Expand(esd->GetNumberOfCaloClusters()-fNCaloClustersOld) ;
  for (Int_t iClust=fNCaloClustersOld; iClust<esd->GetNumberOfCaloClusters(); ++iClust) {
    
    AliESDCaloCluster * cluster = esd->GetCaloCluster(iClust);
    
    Int_t  id        = cluster->GetID();
    Int_t  nLabel    = cluster->GetNLabels();
    Int_t *labels    = cluster->GetLabels();
    
    Float_t energy = cluster->E();
    Float_t posF[3] = { 0.};
    cluster->GetPosition(posF);
    
    AliAODCaloCluster *caloCluster = new((*fEmbeddedClusters)[jClusters++]) AliAODCaloCluster(id,
                                                                                      nLabel,
                                                                                      labels,
                                                                                      energy,
                                                                                      posF,
                                                                                      NULL,
                                                                                      cluster->GetType(),0);

    
    Double_t dx=cluster->GetTrackDx() ;
    Double_t dz=cluster->GetTrackDz() ;
    Float_t cpv=99. ; //No track matched by default
    TArrayI * itracks = cluster->GetTracksMatched() ;  
    if(itracks->GetSize()>0){
      Int_t iTr = itracks->At(0);
      if(iTr>=0 && iTr<esd->GetNumberOfTracks()){
        AliESDtrack *track = esd->GetTrack(iTr) ;
        Double_t pt = track->Pt() ;
        Short_t charge = track->Charge() ;
        cpv=TestCPV(dx, dz, pt,charge) ;
      }
    }

    caloCluster->SetCaloCluster(cluster->GetDistanceToBadChannel(),
                                cluster->GetDispersion(),
                                cluster->GetM20(), cluster->GetM02(),
                                cpv,  
                                cluster->GetNExMax(),cluster->GetTOF()) ;
    
    caloCluster->SetPIDFromESD(cluster->GetPID());
    caloCluster->SetNCells(cluster->GetNCells());    
    caloCluster->SetCellsAbsId(cluster->GetCellsAbsId());
    caloCluster->SetCellsAmplitudeFraction(cluster->GetCellsAmplitudeFraction());            
  } 
  
  for(Int_t i=0;i<fEmbeddedClusters->GetEntriesFast();i++){
    AliAODCaloCluster *caloCluster =(AliAODCaloCluster *)fEmbeddedClusters->At(i) ;
    caloCluster->GetID() ;
    
  }
  
  //Now cells
  if (esd->GetPHOSCells()) { // protection against missing ESD information
    AliESDCaloCells &esdPHcells = *(esd->GetPHOSCells());
    Int_t nPHcell = esdPHcells.GetNumberOfCells() ;    
    fEmbeddedCells->CreateContainer(nPHcell);
    fEmbeddedCells->SetType(AliAODCaloCells::kPHOSCell);
    for (Int_t iCell = 0; iCell < nPHcell; iCell++) {      
      fEmbeddedCells->SetCell(iCell,esdPHcells.GetCellNumber(iCell),esdPHcells.GetAmplitude(iCell));
    }
    fEmbeddedCells->Sort();
  }

  
}
//______________________________________________________________________________
void AliPHOSEmbedding::ConvertMCParticles(const AliAODEvent* aod)
{
  //Copy MC branches to new AOD
  
  TClonesArray *mcArray = (TClonesArray*)aod->FindListObject(AliAODMCParticle::StdBranchName());
  TClonesArray * aodMcParticles = (TClonesArray*)AODEvent()->FindListObject(AliAODMCParticle::StdBranchName());
  for(Int_t i=0;i<mcArray->GetEntriesFast();i++){
    
    AliAODMCParticle* aodpart =  (AliAODMCParticle*) mcArray->At(i);
//    AliAODMCParticle * mcp =new ((*aodMcParticles)[i]) AliAODMCParticle(*aodpart);
    new ((*aodMcParticles)[i]) AliAODMCParticle(*aodpart);
    
  }

  /*  
  // MC header
  AliAODMCHeader * mcHeader = (AliAODMCHeader *)aod.FindListObject(AliAODMCHeader::StdBranchName());
  AliAODMCHeader * aodMcHeader = new AliAODMCHeader(*mcHeader);
  aodMcHeader->SetName(AliAODMCHeader::StdBranchName());
  AddAODBranch("AliAODMCHeader",&aodMcHeader);    
  */

}
//______________________________________________________________________________
void AliPHOSEmbedding::ConvertHeader(AliESDEvent & esd){
  
  AliAODHeader* header = AODEvent()->GetHeader();
  
  header->SetRunNumber(esd.GetRunNumber());
  header->SetOfflineTrigger(fInputHandler->IsEventSelected()); // propagate the decision of the physics selection

  TTree* tree = fInputHandler->GetTree();
  if (tree) {
    TFile* file = tree->GetCurrentFile();
    if (file) header->SetESDFileName(file->GetName());
  }
  
  header->SetBunchCrossNumber(esd.GetBunchCrossNumber());
  header->SetOrbitNumber(esd.GetOrbitNumber());
  header->SetPeriodNumber(esd.GetPeriodNumber());
  header->SetEventType(esd.GetEventType());
  
  header->SetEventNumberESDFile(esd.GetHeader()->GetEventNumberInFile());
  if(const_cast<AliESDEvent&>(esd).GetCentrality()){
    header->SetCentrality(const_cast<AliESDEvent&>(esd).GetCentrality());
  }
  else{
    header->SetCentrality(0);
  }
  if(const_cast<AliESDEvent&>(esd).GetEventplane()){
    header->SetEventplane(const_cast<AliESDEvent&>(esd).GetEventplane());
  }
  else{
    header->SetEventplane(0);
  }

  // Trigger
  header->SetFiredTriggerClasses(esd.GetFiredTriggerClasses());
  header->SetTriggerMask(esd.GetTriggerMask()); 
  header->SetTriggerCluster(esd.GetTriggerCluster());
  header->SetL0TriggerInputs(esd.GetHeader()->GetL0TriggerInputs());    
  header->SetL1TriggerInputs(esd.GetHeader()->GetL1TriggerInputs());    
  header->SetL2TriggerInputs(esd.GetHeader()->GetL2TriggerInputs());    
  
  header->SetMagneticField(esd.GetMagneticField());
  header->SetMuonMagFieldScale(esd.GetCurrentDip()/6000.);
  header->SetZDCN1Energy(esd.GetZDCN1Energy());
  header->SetZDCP1Energy(esd.GetZDCP1Energy());
  header->SetZDCN2Energy(esd.GetZDCN2Energy());
  header->SetZDCP2Energy(esd.GetZDCP2Energy());
  header->SetZDCEMEnergy(esd.GetZDCEMEnergy(0),esd.GetZDCEMEnergy(1));
  
  // ITS Cluster Multiplicty
  const AliMultiplicity *mult = esd.GetMultiplicity();
  for (Int_t ilay = 0; ilay < 6; ilay++) header->SetITSClusters(ilay, mult->GetNumberOfITSClusters(ilay));
  
  // TPC only Reference Multiplicty
  //  Int_t refMult  = fTPCaloneTrackCuts ? (Short_t)fTPCaloneTrackCuts->GetReferenceMultiplicity(&esd, kTRUE) : -1;
  header->SetTPConlyRefMultiplicity(-1);
  
  //
  Float_t diamxy[2]={esd.GetDiamondX(),esd.GetDiamondY()};
  Float_t diamcov[3]; 
  esd.GetDiamondCovXY(diamcov);
  header->SetDiamond(diamxy,diamcov);
  header->SetDiamondZ(esd.GetDiamondZ(),esd.GetSigma2DiamondZ());

  
  //PHOS matrixes
  char path[255] ;
  TGeoHMatrix * m ;
  for(Int_t mod=0; mod<5; mod++){
    snprintf(path,255,"/ALIC_1/PHOS_%d",mod+1) ; //In Geometry modules numbered 1,2,.,5
    if (gGeoManager->cd(path)){
      m = gGeoManager->GetCurrentMatrix() ;
      header->SetPHOSMatrix(new TGeoHMatrix(*m),mod) ;
    }
    else{
      header->SetPHOSMatrix(NULL,mod) ;
    }
  } 
  
  //EventPlane
  Double_t epQ=TPCrp(&esd);
  //Standard Eventplane setter is too complicated...
  //use durty hack
  header->SetZDCN1Energy(epQ) ;
  
  AliCentrality *centrality = esd.GetCentrality(); 
  Double_t acentrality=centrality->GetCentralityPercentileUnchecked("V0M");
  header->SetZDCN2Energy(acentrality) ;

  
}//______________________________________________________________________________
void AliPHOSEmbedding::ConvertPrimaryVertices(AliESDEvent const&esd){
  // Access to the AOD container of vertices

  // Access to the AOD container of vertices
  Int_t fNumberOfVertices = 0;
    
  Double_t pos[3] = { 0. };
  Double_t covVtx[6] = { 0. };

  // Add primary vertex. The primary tracks will be defined
  // after the loops on the composite objects (V0, cascades, kinks)
  const AliESDVertex *vtx = esd.GetPrimaryVertex();
  
  vtx->GetXYZ(pos); // position
  vtx->GetCovMatrix(covVtx); //covariance matrix
  
  TClonesArray * vertexes=AODEvent()->GetVertices();

  AliAODVertex* primaryVertex = new((*vertexes)[fNumberOfVertices++])
  AliAODVertex(pos, covVtx, vtx->GetChi2toNDF(), NULL, -1, AliAODVertex::kPrimary);
  primaryVertex->SetName(vtx->GetName());
  primaryVertex->SetTitle(vtx->GetTitle());
  
  TString vtitle = vtx->GetTitle();
  if (!vtitle.Contains("VertexerTracks")) 
    primaryVertex->SetNContributors(vtx->GetNContributors());
  
  if (fDebug > 0) primaryVertex->Print();  
  
  // Add SPD "main" vertex 
  const AliESDVertex *vtxS = esd.GetPrimaryVertexSPD();
  vtxS->GetXYZ(pos); // position
  vtxS->GetCovMatrix(covVtx); //covariance matrix
  AliAODVertex * mVSPD = new((*vertexes)[fNumberOfVertices++])
  AliAODVertex(pos, covVtx, vtxS->GetChi2toNDF(), NULL, -1, AliAODVertex::kMainSPD);
  mVSPD->SetName(vtxS->GetName());
  mVSPD->SetTitle(vtxS->GetTitle());
  mVSPD->SetNContributors(vtxS->GetNContributors()); 
  
  // Add SPD pileup vertices
  for(Int_t iV=0; iV<esd.GetNumberOfPileupVerticesSPD(); ++iV)
  {
    const AliESDVertex *vtxP = esd.GetPileupVertexSPD(iV);
    vtxP->GetXYZ(pos); // position
    vtxP->GetCovMatrix(covVtx); //covariance matrix
    AliAODVertex * pVSPD =  new((*vertexes)[fNumberOfVertices++])
    AliAODVertex(pos, covVtx, vtxP->GetChi2toNDF(), NULL, -1, AliAODVertex::kPileupSPD);
    pVSPD->SetName(vtxP->GetName());
    pVSPD->SetTitle(vtxP->GetTitle());
    pVSPD->SetNContributors(vtxP->GetNContributors()); 
  }
  
  // Add TRK pileup vertices
  for(Int_t iV=0; iV<esd.GetNumberOfPileupVerticesTracks(); ++iV)
  {
    const AliESDVertex *vtxP = esd.GetPileupVertexTracks(iV);
    vtxP->GetXYZ(pos); // position
    vtxP->GetCovMatrix(covVtx); //covariance matrix
    AliAODVertex * pVTRK = new((*vertexes)[fNumberOfVertices++])
    AliAODVertex(pos, covVtx, vtxP->GetChi2toNDF(), NULL, -1, AliAODVertex::kPileupTracks);
    pVTRK->SetName(vtxP->GetName());
    pVTRK->SetTitle(vtxP->GetTitle());
    pVTRK->SetNContributors(vtxP->GetNContributors());
  }

  
}
//__________________________________________________________________________________
void AliPHOSEmbedding::MakeDigits(AliAODEvent * signal){
  
    //-------------------------------------------------------------------------------------
  //Transform CaloCells into Digits
  //-------------------------------------------------------------------------------------
  
  fDigitsArr->Clear() ;
  fDigitsTree->Branch("PHOS","TClonesArray", &fDigitsArr, 32000);

  //First copy data digits
  Int_t ndigit=0 ;
  for (Short_t icell = 0; icell < fCellsPHOS->GetNumberOfCells(); icell++) {
    Short_t id=0;
    Double_t time=0., amp=0. ;
    if (fCellsPHOS->GetCell(icell, id, amp, time) != kTRUE)
      break;
        
    new((*fDigitsArr)[ndigit]) AliPHOSDigit(-1,id,float(amp),float(time),ndigit);
    ndigit++;
  }


  //Add Digits from Signal
  TClonesArray sdigits("AliPHOSDigit",0) ;
  Int_t isdigit=0 ;
  if(signal){
    AliAODCaloCells* cellsS = signal->GetPHOSCells();
//    Int_t cellLabels[cellsS->GetNumberOfCells()]={0} ;
//    Int_t cellSecondLabels[cellsS->GetNumberOfCells()]={0} ;
    Int_t cellLabels[1000]={0} ;       //1000 should be enough for simulated
    Int_t cellSecondLabels[1000]={0} ; //low-statistics event.
    for(Int_t i=0;i<cellsS->GetNumberOfCells();i++){
      cellLabels[i]=-1 ;
      cellSecondLabels[i]=-1;
   }
//  printf("===========Signal clusters==============\n") ;
    //------------------------------------------------------------------------------------
    //Ancestry information
    //------------------------------------------------------------------------------------
    sdigits.Expand(cellsS->GetNumberOfCells());
    for(Int_t i=0; i<signal->GetNumberOfCaloClusters(); i++) {    
      //cluster from embedded signal
      AliVCluster *clus = signal->GetCaloCluster(i);  
    
      if(!clus->IsPHOS())
        continue;
/*    
    printf("Signal clu(%d): E=%f \n",i,clus->E()) ;
    UShort_t * ind    = clus->GetCellsAbsId() ;
    for(Int_t ic=0; ic < clus->GetNCells(); ic++ )
      printf("Dig(%d)=%d, ",ic,ind[ic]) ;
  printf("\n") ;
*/   
    
      Int_t label = clus->GetLabel();
      Int_t label2 = -1 ;
      if (clus->GetNLabels()>=2) label2 = clus->GetLabelAt(1) ;
    
      UShort_t * index    = clus->GetCellsAbsId() ;
    
      for(Int_t ic=0; ic < clus->GetNCells(); ic++ ){
        for (Int_t icell = 0; icell < cellsS->GetNumberOfCells(); icell++){
	   Short_t cellNumber;
	   Double_t cellAmplitude=0., cellTime=0. ;
	   cellsS->GetCell(icell, cellNumber, cellAmplitude, cellTime) ;
	   if(cellNumber==index[ic]){
	      cellLabels[icell]=label;
              cellSecondLabels[icell]=label2;
	      break ;
	  }
        }
      }
    }
// printf("================End Signal==================\n") ;

    for (Int_t icell = 0; icell < cellsS->GetNumberOfCells(); icell++) {
      Short_t cellNumber;
      Double_t cellAmplitude=0., cellTime=0. ;
      if (cellsS->GetCell(icell, cellNumber, cellAmplitude, cellTime) != kTRUE)
        break;
      //Add only digits related to the cluster, no noisy digits...
      if(cellLabels[icell]==-1)
        continue ;
    
      new(sdigits[isdigit]) AliPHOSDigit(cellLabels[icell],cellNumber,float(cellAmplitude),float(cellTime),isdigit);    
      isdigit++;
    }
  }
  
  //If necessary, take into account difference between calibration used to produce ESD and
  //final calibration
  if(fOldPHOSCalibration[0]){ //there is a difference
    for(Int_t i=0; i<isdigit;i++){
      AliPHOSDigit * sdigit=static_cast<AliPHOSDigit*>(sdigits.At(i)) ;
      Int_t relId[4] ;
      AliPHOSGeometry::GetInstance()->AbsToRelNumbering(sdigit->GetId(),relId);
      Float_t calibESD=fOldPHOSCalibration[relId[0]-1]->GetBinContent(relId[2],relId[3]) ;
      Float_t calibNew=fPHOSReconstructor->Calibrate(1.,sdigit->GetId()) ;
      if(calibNew>0.)
        sdigit->SetEnergy(sdigit->GetEnergy()*calibESD/calibNew) ;
    }
  }
  
  //Merge digits  
  Int_t icurrent = 0 ; //index of the last used digit in underlying event
  fDigitsArr->Expand(ndigit+isdigit) ;
  for(Int_t i=0; i<isdigit;i++){
    AliPHOSDigit * sdigit=static_cast<AliPHOSDigit*>(sdigits.At(i)) ;
    Bool_t added=kFALSE ;
    for(Int_t id=icurrent;id<ndigit;id++){
      AliPHOSDigit * digit=static_cast<AliPHOSDigit*>(fDigitsArr->At(id)) ;
      if(sdigit->GetId() ==  digit->GetId() ){
        *digit += *sdigit ;  //add energies
	icurrent=id+1 ;
        added=kTRUE ;
	break ; //no more digits with same ID in the list
      }
      if(sdigit->GetId() < digit->GetId() ){
        icurrent=id ;
        break ; //no more digits with same ID in the list
      }
    }
    if(!added){
      new((*fDigitsArr)[ndigit]) AliPHOSDigit(*sdigit) ;
      ndigit++ ;
    }
  }
 
  //Change Amp back from Energy to ADC counts
  for(Int_t i=0; i<ndigit;i++){
     AliPHOSDigit * digit=static_cast<AliPHOSDigit*>(fDigitsArr->At(i)) ;
     Float_t calib = 0 ; 
     if(fOldPHOSCalibration[0]){ //Use old calibration
        Int_t relId[4] ;
	AliPHOSGeometry::GetInstance()->AbsToRelNumbering(digit->GetId(),relId);
        calib=fOldPHOSCalibration[relId[0]-1]->GetBinContent(relId[2],relId[3]) ;
     }
     else{
       calib=fPHOSReconstructor->Calibrate(1.,digit->GetId()) ;	
     }
     if(calib>0.)
       digit->SetEnergy(digit->GetEnergy()/calib) ;
  }  
  
   
  fDigitsArr->Sort() ;
  for (Int_t i = 0 ; i < ndigit ; i++) { 
    AliPHOSDigit *digit = static_cast<AliPHOSDigit*>( fDigitsArr->At(i) ) ; 
    digit->SetIndexInList(i) ;     
  }
  fDigitsTree->Fill() ;

}
//__________________________________________________________________________________
void AliPHOSEmbedding::SetOldCalibration(TH2F **calibH){ 
  //Copy histograms with calibration coeff
  gROOT->cd() ;
  char name[55] ;
  for(Int_t i=0;i<5;i++){
    fOldPHOSCalibration[i] = new TH2F(*((TH2F*)calibH[i])) ;
    snprintf(name,55,"OldPHOSCalibrationMod%d",i) ;
    fOldPHOSCalibration[i]->SetName(name) ;
  }
}  
//___________________________________________________________________________
Double_t AliPHOSEmbedding::TPCrp(const AliESDEvent * event){
  //calculate reaction plain in TPC
 
  Double_t rp=999 ; //!Reaction plain calculated with full TPC 
  
  Double_t cosF=0.,sinF=0.,nF=0.;
  for (Int_t i=0;i<event->GetNumberOfTracks();i++) {
    AliESDtrack *track = event->GetTrack(i) ;
    if(!SelectTrack(track))
      continue ;
    Double_t eta= track->Eta() ;
    if(TMath::Abs(eta)> 0.9)
      continue ;
    Double_t phi=track->Phi();
    Double_t cos2phi=TMath::Cos(2.*phi) ;
    Double_t sin2phi=TMath::Sin(2.*phi) ;
    cosF+=cos2phi ;
    sinF+=sin2phi ;
    nF+=1. ;
  }
  if(nF>0){
    rp=0.5*TMath::ATan2(sinF,cosF) ;
  }
  return rp ;
  
}
//___________________________________________________________________________
Bool_t AliPHOSEmbedding::SelectTrack(AliESDtrack * t){
  //estimate if this track can be used for the RP calculation
  Float_t pt=t->Pt();
  if(pt<0.15 || pt>5.)
    return kFALSE ;
//  if(!fESDtrackCuts->AcceptTrack(t))
//    return kFALSE ;
  if(t->GetTPCNcls()<70) 
    return kFALSE;
/*
  Float_t dcaXY=t->DCA();
  Float_t dcaZ=t->ZAtDCA() ;
  t->GetImpactParametersTPC(&dcaXY,&dcaZ);
  if (TMath::Abs(dcaZ)>3)
    return kFALSE;
  if (TMath::Abs(dcaXY)>3.) 
    return kFALSE;
*/  
  return kTRUE ;


}
//____________________________________________________________________________
void AliPHOSEmbedding::InitMF(){
  
  printf("............Init MF \n") ;
  //------------------------------------
  // Initialization of the Mag.Fiels from GRP entry
  // Copied from AliReconstruction::InitGRP()
  //------------------------------------
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  AliGRPObject * aGRPData=0 ;
  if (entry) {

    TMap* m = dynamic_cast<TMap*>(entry->GetObject());  // old GRP entry

    if (m) {
       AliInfo("Found a TMap in GRP/GRP/Data, converting it into an AliGRPObject");
       m->Print();
       aGRPData = new AliGRPObject();
       aGRPData->ReadValuesFromMap(m);
    }

    else {
       AliInfo("Found an AliGRPObject in GRP/GRP/Data, reading it");
       aGRPData = dynamic_cast<AliGRPObject*>(entry->GetObject());  // new GRP entry
       entry->SetOwner(0);
    }
  }

  if (!aGRPData) {
     AliError("No GRP entry found in OCDB!");
     return ;
  }
  //*** Dealing with the magnetic field map
    
  TString lhcState = aGRPData->GetLHCState();
  if (lhcState==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the LHC state ! Using UNKNOWN");
    lhcState = "UNKNOWN";
  }

  TString beamType = aGRPData->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
    beamType = "UNKNOWN";
  }

  Float_t beamEnergy = aGRPData->GetBeamEnergy();
  if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
    beamEnergy = 0;
  }

  TString runType = aGRPData->GetRunType();
  if (runType==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the run type ! Using UNKNOWN");
    runType = "UNKNOWN";
  }
  
  if ( TGeoGlobalMagField::Instance()->IsLocked() ) {
    if (TGeoGlobalMagField::Instance()->GetField()->TestBit(AliMagF::kOverrideGRP)) {
      AliInfo("PHOSEmbedding: MF is locked - GRP information will be ignored !");
      AliInfo("Running with the externally locked B field !");
    }
    else {
      AliInfo("Destroying existing B field instance!");
      delete TGeoGlobalMagField::Instance();
    }   
  }
  if ( !TGeoGlobalMagField::Instance()->IsLocked() ) {
    // Construct the field map out of the information retrieved from GRP.
    Bool_t ok = kTRUE;
    // L3
    Float_t l3Current = aGRPData->GetL3Current((AliGRPObject::Stats)0);
    if (l3Current == AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 current !");
      ok = kFALSE;
    }

    Char_t l3Polarity = aGRPData->GetL3Polarity();
    if (l3Polarity == AliGRPObject::GetInvalidChar()) {
      AliError("GRP/GRP/Data entry:  missing value for the L3 polarity !");
      ok = kFALSE;
    }

    // Dipole
    Float_t diCurrent = aGRPData->GetDipoleCurrent((AliGRPObject::Stats)0);
    if (diCurrent == AliGRPObject::GetInvalidFloat()) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole current !");
      ok = kFALSE;
    }

    Char_t diPolarity = aGRPData->GetDipolePolarity();
    if (diPolarity == AliGRPObject::GetInvalidChar()) {
      AliError("GRP/GRP/Data entry:  missing value for the dipole polarity !");
      ok = kFALSE;
    }

    // read special bits for the polarity convention and map type
    Int_t  polConvention = aGRPData->IsPolarityConventionLHC() ? AliMagF::kConvLHC : AliMagF::kConvDCS2008;
    Bool_t uniformB = aGRPData->IsUniformBMap();

    if (ok) {
      AliMagF* fld = AliMagF::CreateFieldMap(TMath::Abs(l3Current) * (l3Polarity ? -1:1),
                                             TMath::Abs(diCurrent) * (diPolarity ? -1:1),
                                             polConvention,uniformB,beamEnergy, beamType.Data());
      if (fld) {
        TGeoGlobalMagField::Instance()->SetField( fld );
        TGeoGlobalMagField::Instance()->Lock();
        AliInfo("Running with the B field constructed out of GRP !");
      }
      else AliFatal("Failed to create a B field map !");
    }
    else AliFatal("B field is neither set nor constructed from GRP ! Exitig...");
  }
  
}
//____________________________________________________________________________
void AliPHOSEmbedding::InitGeometry(){

  // Import ideal TGeo geometry and apply misalignment
  if (!gGeoManager) {
    AliGeomManager::LoadGeometry("geometry.root");
    if (!gGeoManager) {
      AliFatal("Can not load geometry");
    }
    if(!AliGeomManager::CheckSymNamesLUT("PHOS")) {
      AliFatal("CheckSymNamesLUT");
    }
  }
    

  TString detStr = "PHOS";
  TString loadAlObjsListOfDets = "PHOS";
        
  if(AliGeomManager::GetNalignable("GRP") != 0)
      loadAlObjsListOfDets.Prepend("GRP "); //add alignment objects for non-sensitive modules
  AliGeomManager::ApplyAlignObjsFromCDB(loadAlObjsListOfDets.Data());

  AliCDBManager::Instance()->UnloadFromCache("*/Align/*");
  AliCDBManager::Instance()->UnloadFromCache("GRP/Geometry/Data");
 
}
//____________________________________________________________________________
Float_t AliPHOSEmbedding::TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge){
  //Parameterization of LHC10h period
  //_true if neutral_
  
  Double_t z0=0.84 ; //Both charges
  Double_t sz=TMath::Min(2.75,2.*3.*3./(pt*pt+3.*3.)+1.);
  Double_t x0=TMath::Min(7.6, 27.4628*TMath::Exp(-(pt+13.3572)*(pt+13.3572)/2./6.29214/6.29214)+
                                    243.020*0.0526626*0.0526626/(pt*pt+0.0526626*0.0526626)) ;
  if(charge>0)
    x0=-x0 ;
  Double_t sx=TMath::Min(5.4,9.9*TMath::Exp(-pt/0.31)+0.5*0.4*0.4/((pt-1.2)*(pt-1.2)+0.4*0.4)+1.5);
  Double_t rz=(dz-z0)/sz ;
  Double_t rx=(dx-x0)/sx ;
  return TMath::Sqrt(rx*rx+rz*rz)/2. ;
}
