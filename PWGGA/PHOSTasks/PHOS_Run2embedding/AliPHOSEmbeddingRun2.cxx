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
// updated to support Run2, Daiki Sekihata


#include "TChain.h"
#include "TFile.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TObjArray.h"
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
#include "AliMultSelection.h"

#include "AliPHOSReconstructor.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSDigit.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSCalibData.h"

#include "AliOADBContainer.h"
#include "AliPHOSEmbeddingRun2.h"

ClassImp(AliPHOSEmbeddingRun2)


//________________________________________________________________________
AliPHOSEmbeddingRun2::AliPHOSEmbeddingRun2(const char *name) 
  : AliAnalysisTaskESDfilter(name),
    fAODChain1(0x0),
    fAODChain2(0x0),
    fAODChain3(0x0),
    fSignal1(0x0),
    fSignal2(0x0),
    fSignal3(0x0),
    fDigitsTree(0x0),
    fClustersTree(0x0),
    fTreeOut(0x0),
    fDigitsArr(0x0),
    fEmbeddedClusters1(0x0),
    fEmbeddedClusters2(0x0),
    fEmbeddedClusters3(0x0),
    fEmbeddedCells1(0x0),
    fEmbeddedCells2(0x0),
    fEmbeddedCells3(0x0),
    fCellsPHOS(0x0),
    fClusterizer(0x0),
    fPHOSReconstructor(0x0),
    fSignalCalibData(0x0),
    fPathPrivateOADB("alien://"),
    fSignalECorrection(1.),
    fNSignal(0),
  	fNSignalPi0(0),
	  fNSignalEta(0), 
	  fNSignalGamma(0),
    fNCaloClustersOld(0),
    fInitialized(0),
		fIsPi0Embedded(kFALSE),
		fIsEtaEmbedded(kFALSE),
		fIsGammaEmbedded(kFALSE)
{
  // Constructor
  for(int i=0;i<5;i++)fOldPHOSCalibration[i]=0x0 ;

}

//________________________________________________________________________
void AliPHOSEmbeddingRun2::UserCreateOutputObjects()
{
 //prepare output containers
  AliAnalysisTaskESDfilter::UserCreateOutputObjects() ;
  
  //Create branch for output MC particles  
  TClonesArray *mcparticles1 = new TClonesArray("AliAODMCParticle",0);
  mcparticles1->SetName(Form("%s_pi0",AliAODMCParticle::StdBranchName()));
  AddAODBranch("TClonesArray", &mcparticles1);

  //Branch with clusters after embedding  ---Pi0---
  fEmbeddedClusters1 = new TClonesArray("AliAODCaloCluster",0);
  fEmbeddedClusters1->SetName("EmbeddedPi0CaloClusters");
  AddAODBranch("TClonesArray", &fEmbeddedClusters1);
  
  //Branch with cells after embedding
  fEmbeddedCells1 = new AliAODCaloCells();
  fEmbeddedCells1->SetName("EmbeddedPi0PHOScells");
  AddAODBranch("AliAODCaloCells", &fEmbeddedCells1);

  //Create branch for output MC particles  ---Eta---
  TClonesArray *mcparticles2 = new TClonesArray("AliAODMCParticle",0);
  mcparticles2->SetName(Form("%s_eta",AliAODMCParticle::StdBranchName()));
  AddAODBranch("TClonesArray", &mcparticles2);

  //Branch with clusters after embedding
  fEmbeddedClusters2 = new TClonesArray("AliAODCaloCluster",0);
  fEmbeddedClusters2->SetName("EmbeddedEtaCaloClusters");
  AddAODBranch("TClonesArray", &fEmbeddedClusters2);
  
  //Branch with cells after embedding
  fEmbeddedCells2 = new AliAODCaloCells();
  fEmbeddedCells2->SetName("EmbeddedEtaPHOScells");
  AddAODBranch("AliAODCaloCells", &fEmbeddedCells2);

  //Create branch for output MC particles  ---Gamma---
  TClonesArray *mcparticles3 = new TClonesArray("AliAODMCParticle",0);
  mcparticles3->SetName(Form("%s_gamma",AliAODMCParticle::StdBranchName()));
  AddAODBranch("TClonesArray", &mcparticles3);

  //Branch with clusters after embedding
  fEmbeddedClusters3 = new TClonesArray("AliAODCaloCluster",0);
  fEmbeddedClusters3->SetName("EmbeddedGammaCaloClusters");
  AddAODBranch("TClonesArray", &fEmbeddedClusters3);
  
  //Branch with cells after embedding
  fEmbeddedCells3 = new AliAODCaloCells();
  fEmbeddedCells3->SetName("EmbeddedGammaPHOScells");
  AddAODBranch("AliAODCaloCells", &fEmbeddedCells3);

  
  fDigitsArr = new TClonesArray("AliPHOSDigit",4*56*64) ;

  AliAODHandler* handler = (AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
  if (handler) {
      fTreeOut = handler->GetTree();
  }
  else {
       AliWarning("No output AOD Event Handler connected.") ;
  }

}

//________________________________________________________________________
void AliPHOSEmbeddingRun2::UserExec(Option_t *) {
  // Main loop, called for each event
  // Perform embedding
  
  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());//for general information before embedding
  if (!event) {
    Printf("ERROR: Could not retrieve event");
    PostData(0, fTreeOut);
    return;
  }

  TString trigClasses = event->GetFiredTriggerClasses();
	//cout << "trigger classes = " << event->GetFiredTriggerClasses() << endl;
  AliInfo(Form("trig.class %s, period %d, bc %d, orbit %d", trigClasses.Data(),event->GetPeriodNumber(), event->GetBunchCrossNumber(),event->GetOrbitNumber()));

  Float_t acentrality=-1;
  AliCentrality *centrality = event->GetCentrality(); 

	//for Run2
	AliMultSelection *fMultSelection = (AliMultSelection*)event->FindListObject("MultSelection");
	if(!fMultSelection){
   AliWarning("AliMultSelection object not found!");
   return;
	}
	else{
		acentrality = fMultSelection->GetMultiplicityPercentile("V0M");
	}

  //Float_t acentrality=centrality->GetCentralityPercentile("V0M");
  if( acentrality <= 0. || acentrality>101.){
    return;
  }

	//cout << "acentrality = " << acentrality << endl;

  //Read next AOD event
  //If necesary method checks if AOD event is good to embed
  //e.g. there are PHOS clusters etc.
	//GetNextSignalEvent() ;
	if(fIsPi0Embedded) GetNextSignalEventPi0() ;
	if(fIsEtaEmbedded) GetNextSignalEventEta() ;
	if(fIsGammaEmbedded) GetNextSignalEventGamma() ;

//  if (!signal) {
//    Printf("ERROR: Could not retrieve signal");
//    PostData(0, fTreeOut);
//    return;
//  }
  
  Init() ;

  //Remember PHOS digits before any modification
  CopyRecalibrateDigits() ;

  AliInfo(Form("Ncluster in a real collision = %d",event->GetNumberOfCaloClusters()));
  AliInfo(Form("Ncells in real coolisions = %d.",fCellsPHOS->GetNumberOfCells()));

  //Convert ESD to AOD for real collisions
  AliAnalysisTaskESDfilter::UserExec("");

  AliInfo(Form("Ncluster after AliAnalysisTaskESDfilter::UserExec("") = %d",event->GetNumberOfCaloClusters()));

  AliInfo("Next, go to embedding! AliESDEvent will be modified.");

  //First re-reconstruct existing digits to ensure identical reconstruction before and after embeding
  AliAODEvent * nullSignal=0x0 ;
  MakeEmbedding(event,nullSignal);

  AliInfo(Form("Ncells in PHOS = %d.",fCellsPHOS->GetNumberOfCells()));
  AliInfo(Form("Ncluster after MakeEmbedding(event,nullSignal) = %d",event->GetNumberOfCaloClusters()));

	if(fIsPi0Embedded){
	  //Now perform real embedding  
	  //Embed signal clusters into event

	  MakeEmbedding(event,fSignal1);
    AliInfo(Form("Ncluster after MakeEmbedding for Pi0 = %d",event->GetNumberOfCaloClusters()));
	
	  //Convert ESD to AOD
	  ConvertEmbeddedClusters(event,1) ;
	 
	  //Fill MC information
	  ConvertMCParticles(fSignal1,1) ;
    AliInfo(Form("Ncells in PHOS = %d.",fCellsPHOS->GetNumberOfCells()));
    AliInfo(Form("Ncluster after embedding Pi0 = %d",event->GetNumberOfCaloClusters()));
	}


	if(fIsEtaEmbedded){
	  //Eta
	  MakeEmbedding(event,fSignal2) ;
    AliInfo(Form("Ncluster after MakeEmbedding for Eta = %d",event->GetNumberOfCaloClusters()));

	  //Convert ESD to AOD
	  ConvertEmbeddedClusters(event,2) ;

	  //Fill MC information
	  ConvertMCParticles(fSignal2,2) ;
    AliInfo(Form("Ncells in PHOS = %d.",fCellsPHOS->GetNumberOfCells()));
    AliInfo(Form("Ncluster after embedding Eta = %d",event->GetNumberOfCaloClusters()));
	}
	if(fIsGammaEmbedded){
	  //gamma
	  MakeEmbedding(event,fSignal3) ;
    AliInfo(Form("Ncluster after MakeEmbedding for Gamma = %d",event->GetNumberOfCaloClusters()));

	  //Convert ESD to AOD
	  ConvertEmbeddedClusters(event,3) ;

	  //Fill MC information
	  ConvertMCParticles(fSignal3,3) ;
    AliInfo(Form("Ncells in PHOS = %d.",fCellsPHOS->GetNumberOfCells()));
    AliInfo(Form("Ncluster after embedding Gamma = %d",event->GetNumberOfCaloClusters()));
 	} 
  
  AliInfo(Form("Ncluster after all embedding procedure = %d",event->GetNumberOfCaloClusters()));

  PostData(0, fTreeOut);
}
//______________________________________________________________________________
void AliPHOSEmbeddingRun2::CopyRecalibrateDigits(){
  //Recalibrate digits if there is better calibration ("new") 
  //exists than one used in reconstruction ("ESD")
  
  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
  if(fCellsPHOS)
    delete fCellsPHOS ;
  fCellsPHOS = new AliESDCaloCells() ;
  fCellsPHOS->CreateContainer(event->GetPHOSCells()->GetNumberOfCells());

  for (Short_t icell = 0; icell < event->GetPHOSCells()->GetNumberOfCells(); icell++) {
      Short_t id=0;
      Double_t time=0., amp=0. ;
      Int_t mclabel; 
      Double_t efrac =0. ;
      if (event->GetPHOSCells()->GetCell(icell, id, amp,time,mclabel,efrac) != kTRUE)
        break;      
      fCellsPHOS->SetCell(icell, id, amp, time,mclabel,efrac);     
   }




}
//______________________________________________________________________________
void AliPHOSEmbeddingRun2::Init(){
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

  //Get AODB calibration for this run for de-calibration of signal
  AliOADBContainer calibContainer("phosRecalibration");
  calibContainer.InitFromFile(fPathPrivateOADB,"phosRecalibration");
  //calibContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSCalibrations.root","phosRecalibration");
  TObjArray *recalib = (TObjArray*)calibContainer.GetObject(runNum,"PHOSRecalibration");
  if(!recalib){
    AliWarning(Form("Can not read calibrations for run %d, do not apply OADB de-calibration\n",runNum)) ;
  }
  else{
    const Int_t recoPass=1; 
    fSignalCalibData = (AliPHOSCalibData*)recalib->At(recoPass-1) ;
  }

  fInitialized=kTRUE ;
}
//________________________________________________________________________
void AliPHOSEmbeddingRun2::GetNextSignalEvent(){
  //Read signal AOD event from the chain
  
  if(fAODChain1==0 || fAODChain2==0 || fAODChain3==0){
    AliError(Form("No chain to read signal events: pi0=%p, eta=%p, gamma=%p",fAODChain1,fAODChain2,fAODChain3)) ;
    if(fSignal1) delete fSignal1 ; fSignal1=0x0 ;
    if(fSignal2) delete fSignal2 ; fSignal2=0x0 ;
    if(fSignal3) delete fSignal3 ; fSignal3=0x0 ;
    return  ;
  }
 
	//===== pi0 =====
	if(fIsPi0Embedded){
		if(fSignal1) delete fSignal1 ;
		fSignal1 = new AliAODEvent;
		fSignal1->ReadFromTree(fAODChain1);
		
		if(fNSignal>=fAODChain1->GetEntries()){
			delete fSignal1 ;
			fSignal1= 0x0 ;
		}
		else
			fAODChain1->GetEvent(fNSignal) ;
	} 
	//===== eta =====

	if(fIsEtaEmbedded){
		if(fSignal2) delete fSignal2 ;
		fSignal2 = new AliAODEvent;
		fSignal2->ReadFromTree(fAODChain2);
	
		if(fNSignal>=fAODChain1->GetEntries()){
			delete fSignal2 ;
			fSignal2= 0x0 ;
		}
		else
			fAODChain2->GetEvent(fNSignal) ;
	}
	//===== gamma =====
	if(fIsGammaEmbedded){
	  if(fSignal3) delete fSignal3 ;
	  fSignal3 = new AliAODEvent;
	  fSignal3->ReadFromTree(fAODChain3);
	  
	  if(fNSignal>=fAODChain1->GetEntries()){
	    delete fSignal3 ;
	    fSignal3= 0x0 ;
	  }
	  else
	    fAODChain3->GetEvent(fNSignal) ;  
 
	}
 
  fNSignal++ ;
  
}
//________________________________________________________________________
void AliPHOSEmbeddingRun2::GetNextSignalEventPi0(){
  //Read signal AOD event from the chain
  
  if(fAODChain1==0){
    AliError(Form("No chain to read signal events: pi0=%p",fAODChain1)) ;
    if(fSignal1) delete fSignal1 ; fSignal1=0x0 ;
    return  ;
  }
 
	//===== pi0 =====
 
  if(fSignal1) delete fSignal1 ;
  fSignal1 = new AliAODEvent;
  fSignal1->ReadFromTree(fAODChain1);
  
  if(fNSignalPi0>=fAODChain1->GetEntries()){
    delete fSignal1 ;
    fSignal1= 0x0 ;
  }
  else
    fAODChain1->GetEvent(fNSignalPi0) ;
 
  fNSignalPi0++ ;
  
}
//________________________________________________________________________
void AliPHOSEmbeddingRun2::GetNextSignalEventEta(){
  //Read signal AOD event from the chain
  
  if(fAODChain2==0){
    AliError(Form("No chain to read signal events: eta=%p",fAODChain2)) ;
    if(fSignal2) delete fSignal2 ; fSignal2=0x0 ;
    return  ;
  }
 
	//===== eta =====
  
  if(fSignal2) delete fSignal2 ;
  fSignal2 = new AliAODEvent;
  fSignal2->ReadFromTree(fAODChain2);
  
  if(fNSignalEta>=fAODChain2->GetEntries()){
    delete fSignal2 ;
    fSignal2= 0x0 ;
  }
  else
    fAODChain2->GetEvent(fNSignalEta) ;
 
  fNSignalEta++ ;
  
}
//________________________________________________________________________
void AliPHOSEmbeddingRun2::GetNextSignalEventGamma(){
  //Read signal AOD event from the chain
  
  if(fAODChain3==0){
    AliError(Form("No chain to read signal events: gamma=%p",fAODChain3)) ;
    if(fSignal3) delete fSignal3 ; fSignal3=0x0 ;
    return  ;
  }

	//===== gamma =====
  if(fSignal3) delete fSignal3 ;
  fSignal3 = new AliAODEvent;
  fSignal3->ReadFromTree(fAODChain3);
  
  if(fNSignalGamma>=fAODChain3->GetEntries()){
    delete fSignal3 ;
    fSignal3= 0x0 ;
  }
  else
    fAODChain3->GetEvent(fNSignalGamma) ;
 
  fNSignalGamma++ ;

}
//________________________________________________________________________
void AliPHOSEmbeddingRun2::MakeEmbedding(AliESDEvent *event,AliAODEvent * signal){
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

  Int_t nPHOSBefore=0;
  Int_t nCPVBefore=0;
  Int_t nEMCALBefore=0;

  for (Int_t iClust=0; iClust<event->GetNumberOfCaloClusters(); ++iClust) {
    AliESDCaloCluster * cluster = event->GetCaloCluster(iClust);
    //if(!cluster->IsPHOS()) continue;
    //if(cluster->GetType() != AliVCluster::kPHOSNeutral) continue;

    if(cluster->IsPHOS()){
      if(cluster->GetType() == AliVCluster::kPHOSNeutral) nPHOSBefore++;
      if(cluster->GetType() == AliVCluster::kPHOSCharged) nCPVBefore++;
    }
    else{
      nEMCALBefore++;
    }

//    if(cluster->GetType() == AliVCluster::kPHOSNeutral)
//      printf("MC label of PHOS cluster before embedding: ID = %d , energy = %f , label = %d\n",iClust,cluster->E(),cluster->GetLabel());

//    Float_t pos[3] = {0.};
//    cluster->GetPosition(pos);
//    printf("Cluster(%d): E=%f, x=%f, z=%f, CPV=%f \n",iClust,cluster->E(),pos[0],pos[2],cluster->GetEmcCpvDistance()) ;
//    UShort_t * index    = cluster->GetCellsAbsId() ;
//    for(Int_t ic=0; ic < cluster->GetNCells(); ic++ )
//      printf("Dig(%d)=%d, ",ic,index[ic]) ;
//    printf("\n") ;


  }

  AliInfo(Form("Before embedding: Nall = %d , nPHOS = %d , nCPV = %d , nEMCAL = %d.",fNCaloClustersOld, nPHOSBefore, nCPVBefore, nEMCALBefore)) ;

  if(signal){
    AliInfo(Form("MC Ncluster in sigle simulation = %d",signal->GetNumberOfCaloClusters()));
//    for(Int_t i=0;i<signal->GetNumberOfCaloClusters();i++){
//      AliAODCaloCluster *caloCluster = (AliAODCaloCluster*)signal->GetCaloCluster(i);
//      if(caloCluster->GetType() != AliVCluster::kPHOSNeutral) continue;
//      printf("MC label of PHOS cluster in simulation: ID = %d , energy = %f , label = %d\n",i,caloCluster->E(),caloCluster->GetLabel());
//    }

    TClonesArray *mcArray = (TClonesArray*)signal->FindListObject(AliAODMCParticle::StdBranchName());
    AliAODCaloCells *cellsEmb = (AliAODCaloCells*)signal->GetPHOSCells();
    Int_t Ncell = cellsEmb->GetNumberOfCells();
    AliInfo(Form("Before embedding: N PHOS cells in simulation = %d",Ncell));
  }



//  for(Int_t icell=0;icell<Ncell;icell++){
//    Int_t pdg=0;
//    if(cellsEmb->GetMCLabel(icell)>-1) pdg = ((AliAODMCParticle*)mcArray->At(cellsEmb->GetMCLabel(icell)))->GetPdgCode();
//    printf("cell(%d): AbsID = %d , Amp = %f , MC label = %d , pdg = %d\n",icell, cellsEmb->GetCellNumber(icell), cellsEmb->GetAmplitude(icell), cellsEmb->GetMCLabel(icell),pdg);
//  }


  //create digits
  MakeDigits(signal);	  
 
// for(Int_t i=0;i<fDigitsArr->GetEntriesFast();i++){  
//    AliPHOSDigit * d = (AliPHOSDigit *) fDigitsArr->At(i) ;
//    printf("  Digit(%d) = %d, E=%f \n",i,d->GetId(),d->GetEnergy()) ;
// }

  
  //clusterize and make tracking
  fPHOSReconstructor->Reconstruct(fDigitsTree,fClustersTree) ;

  //event->Clear(); 
  //Note that the current ESDEvent will be modified!
  fPHOSReconstructor->FillESD(fDigitsTree, fClustersTree, event) ;
  //this modified event should go to ConvertEmbeddedClusters

  Int_t nPHOSAfter=0;
  Int_t nCPVAfter=0;
  Int_t nEMCALAfter=0;

  //fNCaloClustersOld is the number of clusters just before this embddding step.
  AliInfo(Form("cluster loop to store embedded clusters + real clusters will run between %d-%d.",fNCaloClustersOld,event->GetNumberOfCaloClusters()));

  for (Int_t iClust=fNCaloClustersOld; iClust<event->GetNumberOfCaloClusters(); ++iClust) {
    AliESDCaloCluster * cluster = event->GetCaloCluster(iClust);
    //if(!cluster->IsPHOS()) continue;
    //if(cluster->GetType() != AliVCluster::kPHOSNeutral) continue;

    if(cluster->IsPHOS()){
      if(cluster->GetType() == AliVCluster::kPHOSNeutral) nPHOSAfter++;
      if(cluster->GetType() == AliVCluster::kPHOSCharged) nCPVAfter++;
    }
    else{
      nEMCALAfter++;
    }

//    if(cluster->GetType() == AliVCluster::kPHOSNeutral && signal){
//      if(cluster->GetLabel() > -1) printf("MC label of PHOS cluster after embedding: ID = %d , energy = %f , label = %d\n",iClust,cluster->E(),cluster->GetLabel());
//    }

//    Float_t pos[3] = { 0.};
//    cluster->GetPosition(pos);
//    printf("Cluster(%d): E=%f, x=%f, z=%f, CPV=%f, Label=%d \n",iClust-n,cluster->E(),pos[0],pos[2],cluster->GetEmcCpvDistance(),cluster->GetLabel()) ;
//    UShort_t * index    = cluster->GetCellsAbsId() ;
//    for(Int_t ic=0; ic < cluster->GetNCells(); ic++ )
//      printf("Dig(%d)=%d, ",ic,index[ic]) ;
//    printf("\n") ;

  }
  AliInfo(Form("After embedding: Nall = %d , nPHOSAfter = %d , nCPVAfter = %d , nEMCALAfter = %d.", event->GetNumberOfCaloClusters(), nPHOSAfter, nCPVAfter, nEMCALAfter));
  //AliInfo("At the end of MakeEmbedding(), Nall = N old cluster + N PHOS cluster after embedding + N CPV cluster after embedding + N EMCAL cluster after embedding.This means N old PHOS clusters are double counted. No problem!");

}
//______________________________________________________________________________
void AliPHOSEmbeddingRun2::ConvertEmbeddedClusters(const AliESDEvent* esd,Int_t what)
{
  //Copy PHOS clusters and cells after embedding

  // Access to the AOD container of clusters
  Int_t jClusters(0);
  TClonesArray * embeddedClusters =0x0 ;
  AliAODCaloCells * embeddedCells =0x0 ;
  if(what==1){ //Pi0
    embeddedClusters=fEmbeddedClusters1;
    embeddedCells=fEmbeddedCells1;
  }  
  if(what==2){ //Eta
    embeddedClusters=fEmbeddedClusters2;
    embeddedCells=fEmbeddedCells2;
  }  
  if(what==3){ //Gama
    embeddedClusters=fEmbeddedClusters3;
    embeddedCells=fEmbeddedCells3;
  }  

  embeddedClusters->Clear();
  embeddedClusters->Expand(esd->GetNumberOfCaloClusters()-fNCaloClustersOld) ;

  //fNCaloClustersOld is the number of calo clusters before embedding.
  //this cluster loop will be performed for real cluster + embedded cluster to store Emb_particle_CaloClusters. not only PHOS clusters, but also EMCAL/DCAL clusters will be stored
  for (Int_t iClust=fNCaloClustersOld; iClust<esd->GetNumberOfCaloClusters(); ++iClust) {

    AliESDCaloCluster * cluster = esd->GetCaloCluster(iClust);

    Int_t  id        = cluster->GetID();
    Int_t  nLabel    = cluster->GetNLabels();
    Int_t *labels    = cluster->GetLabels();

    Float_t energy = cluster->E();
    Float_t posF[3] = { 0.};
    cluster->GetPosition(posF);

    AliAODCaloCluster *caloCluster = new((*embeddedClusters)[jClusters++]) AliAODCaloCluster(id,
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
        cpv=TestCPVRun2(dx, dz, pt,charge) ;
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

  for(Int_t i=0;i<embeddedClusters->GetEntriesFast();i++){
    AliAODCaloCluster *caloCluster =(AliAODCaloCluster *)embeddedClusters->At(i) ;
    caloCluster->GetID() ;
  }

  //Now cells
  if (esd->GetPHOSCells()) { // protection against missing ESD information
    AliESDCaloCells &esdPHcells = *(esd->GetPHOSCells());
    Int_t nPHcell = esdPHcells.GetNumberOfCells() ;    
    embeddedCells->CreateContainer(nPHcell);
    embeddedCells->SetType(AliAODCaloCells::kPHOSCell);
    for (Int_t iCell = 0; iCell < nPHcell; iCell++) {      
      //embeddedCells->SetCell(iCell,esdPHcells.GetCellNumber(iCell),esdPHcells.GetAmplitude(iCell),esdPHcells.GetTime(iCell));
      embeddedCells->SetCell(iCell,esdPHcells.GetCellNumber(iCell),esdPHcells.GetAmplitude(iCell),esdPHcells.GetTime(iCell),esdPHcells.GetMCLabel(iCell),esdPHcells.GetEFraction(iCell),esdPHcells.GetHighGain(iCell));

      //cout << "iCell = " << iCell << " , CellNumber = " << esdPHcells.GetCellNumber(iCell) << " , esdPHcells.GetAmplitude(iCell) = " <<esdPHcells.GetAmplitude(iCell) << " , esdPHcells.GetTime(iCell) = " << esdPHcells.GetTime(iCell) 
      //<< " , esdPHcells.GetMCLabel(iCell) = " << esdPHcells.GetMCLabel(iCell) << " , esdPHcells.GetEFraction(iCell) = " << esdPHcells.GetEFraction(iCell) << " , esdPHcells.GetHighGain(iCell) = " << esdPHcells.GetHighGain(iCell) << endl;

    }
    embeddedCells->Sort();
  }

  AliInfo(Form("entry of embeddedClusters = %d in embedded array, after embedding.",embeddedClusters->GetEntriesFast()));

//  for(Int_t i=0;i<embeddedClusters->GetEntriesFast();i++){
//    AliAODCaloCluster *caloCluster =(AliAODCaloCluster *)embeddedClusters->At(i) ;
//    if(caloCluster->GetType() != AliVCluster::kPHOSNeutral) continue;
//    printf("MC label of PHOS cluster after embedding: ID = %d , energy = %f , label = %d\n",i,caloCluster->E(),caloCluster->GetLabel());
//  }
  



}
//______________________________________________________________________________
void AliPHOSEmbeddingRun2::ConvertMCParticles(const AliAODEvent* aod,Int_t what)
{
  //Copy MC branches to new AOD
  
  TClonesArray *mcArray = (TClonesArray*)aod->FindListObject(AliAODMCParticle::StdBranchName());
  TClonesArray * aodMcParticles = 0x0 ;
  if(what==1)
    aodMcParticles = (TClonesArray*)AODEvent()->FindListObject(Form("%s_pi0",AliAODMCParticle::StdBranchName()));
  if(what==2)
    aodMcParticles = (TClonesArray*)AODEvent()->FindListObject(Form("%s_eta",AliAODMCParticle::StdBranchName()));
  if(what==3)
    aodMcParticles = (TClonesArray*)AODEvent()->FindListObject(Form("%s_gamma",AliAODMCParticle::StdBranchName()));

  for(Int_t i=0;i<mcArray->GetEntriesFast();i++){
    AliAODMCParticle* aodpart =  (AliAODMCParticle*) mcArray->At(i);
    new ((*aodMcParticles)[i]) AliAODMCParticle(*aodpart);
    
  }
}
//__________________________________________________________________________________
void AliPHOSEmbeddingRun2::MakeDigits(AliAODEvent * signal){

  //-------------------------------------------------------------------------------------
  //Transform CaloCells into Digits which can be used for standard reconstruction
  //Add signal digits to the event
  //-------------------------------------------------------------------------------------

  fDigitsArr->Clear() ;
  fDigitsTree->Branch("PHOS","TClonesArray", &fDigitsArr, 32000);

  //First copy data digits
  Int_t ndigit=0 ;
  for (Short_t icell = 0; icell < fCellsPHOS->GetNumberOfCells(); icell++) {
    Short_t id=0;
    Double_t time=0., amp=0. ;
    Int_t mclabel;
    Double_t efrac =0. ;
    if (fCellsPHOS->GetCell(icell, id, amp, time,mclabel, efrac) != kTRUE)
      break;
    Int_t idLong=id ;
    if(id<0)idLong= -id+56*64*5 ; //CPV digits
    new((*fDigitsArr)[ndigit]) AliPHOSDigit(-1,idLong,float(amp),float(time),ndigit);
    ndigit++;
  }

  //Add Digits from Signal
  TClonesArray sdigits("AliPHOSDigit",0) ;
  Int_t isdigit=0 ;
  if(signal){
    AliAODCaloCells* cellsS = signal->GetPHOSCells();
    Int_t cellLabels[1000]={0} ;       //1000 should be enough for simulated
    Int_t cellSecondLabels[1000]={0} ; //low-statistics event.
    for(Int_t i=0;i<cellsS->GetNumberOfCells();i++){
      cellLabels[i]=-1 ;
      cellSecondLabels[i]=-1;
    }
    //------------------------------------------------------------------------------------
    //Ancestry information
    //Celect digits contributing to signal clusters and add primary information
    //(it is not stored in CaloCells)
    //------------------------------------------------------------------------------------
    sdigits.Expand(cellsS->GetNumberOfCells());
    for(Int_t i=0; i<signal->GetNumberOfCaloClusters(); i++) {    
      //cluster from embedded signal
      AliVCluster *clus = signal->GetCaloCluster(i);  

      if(!clus->IsPHOS())
        continue;

      Int_t label = clus->GetLabel();
      Int_t label2 = -1 ;
      if (clus->GetNLabels()>=2) label2 = clus->GetLabelAt(1) ;

      UShort_t * index    = clus->GetCellsAbsId() ;
      for(Int_t ic=0; ic < clus->GetNCells(); ic++ ){
        for (Int_t icell = 0; icell < cellsS->GetNumberOfCells(); icell++){
          Short_t cellNumber;
          Double_t cellAmplitude=0., cellTime=0. ;
          Int_t mclabel;
          Double_t efrac =0. ;
          cellsS->GetCell(icell, cellNumber, cellAmplitude, cellTime,mclabel,efrac) ;
          Int_t longCellNumber=cellNumber ;
          if(cellNumber<0)longCellNumber= -cellNumber+56*64*5 ; //CPV digits	   
          if(longCellNumber==index[ic]){
            cellLabels[icell]=label;
            cellSecondLabels[icell]=label2;
            break ;
          }
        }
      }
    }

    for (Int_t icell = 0; icell < cellsS->GetNumberOfCells(); icell++) {
      Short_t cellNumber;
      Double_t cellAmplitude=0., cellTime=0. ;
      Int_t mclabel;
      Double_t efrac =0. ;
      if (cellsS->GetCell(icell, cellNumber, cellAmplitude, cellTime,mclabel,efrac) != kTRUE)
        break;
      //Add only digits related to the cluster, no noisy digits...
      if(cellLabels[icell]==-1){
        continue ;
      }   

      cellAmplitude=DecalibrateSignal(cellAmplitude,cellNumber);//decalibration from Dmitri, Daiki added this line on 21.December.2016

      Int_t longCellNumber=cellNumber ;
      if(cellNumber<0)longCellNumber= -cellNumber+56*64*5 ; //CPV digits
      new(sdigits[isdigit]) AliPHOSDigit(cellLabels[icell],longCellNumber,float(cellAmplitude),float(cellTime),isdigit);    
      isdigit++;
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
  //Note that Reconstructor uses best ("new") calibration
  for(Int_t i=0; i<ndigit;i++){
    AliPHOSDigit * digit=static_cast<AliPHOSDigit*>(fDigitsArr->At(i)) ;
    Float_t calib =fPHOSReconstructor->Calibrate(1.,digit->GetId()) ;	
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
//____________________________________________________________________________
void AliPHOSEmbeddingRun2::InitMF(){
  
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
void AliPHOSEmbeddingRun2::InitGeometry(){

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
Float_t AliPHOSEmbeddingRun2::TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge){
  //Parameterization of LHC10h period
  //_true if neutral_
  
  Double_t meanX=0;
  Double_t meanZ=0.;
  Double_t sx=TMath::Min(5.4,2.59719e+02*TMath::Exp(-pt/1.02053e-01)+
              6.58365e-01*5.91917e-01*5.91917e-01/((pt-9.61306e-01)*(pt-9.61306e-01)+5.91917e-01*5.91917e-01)+1.59219);
  Double_t sz=TMath::Min(2.75,4.90341e+02*1.91456e-02*1.91456e-02/(pt*pt+1.91456e-02*1.91456e-02)+1.60) ;
  
  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
  Double_t mf = event->GetMagneticField(); //Positive for ++ and negative for --

  if(mf<0.){ //field --
    meanZ = -0.468318 ;
    if(charge>0)
      meanX=TMath::Min(7.3, 3.89994*1.20679*1.20679/(pt*pt+1.20679*1.20679)+0.249029+2.49088e+07*TMath::Exp(-pt*3.33650e+01)) ;
    else
      meanX=-TMath::Min(7.7,3.86040*0.912499*0.912499/(pt*pt+0.912499*0.912499)+1.23114+4.48277e+05*TMath::Exp(-pt*2.57070e+01)) ;
  }
  else{ //Field ++
    meanZ= -0.468318;
    if(charge>0)
      meanX=-TMath::Min(8.0,3.86040*1.31357*1.31357/(pt*pt+1.31357*1.31357)+0.880579+7.56199e+06*TMath::Exp(-pt*3.08451e+01)) ;
    else
      meanX= TMath::Min(6.85, 3.89994*1.16240*1.16240/(pt*pt+1.16240*1.16240)-0.120787+2.20275e+05*TMath::Exp(-pt*2.40913e+01)) ;     
  }

  Double_t rz=(dz-meanZ)/sz ;
  Double_t rx=(dx-meanX)/sx ;
  return TMath::Sqrt(rx*rx+rz*rz) ;
}
//____________________________________________________________________________
Float_t AliPHOSEmbeddingRun2::TestCPVRun2(Double_t dx, Double_t dz, Double_t pt, Int_t charge){
  //Parameterization of LHC15o period
  //_true if neutral_

    Double_t meanX=0.;
  Double_t meanZ=0.;

  Double_t sx = TMath::Min(5.2, 1.160 + 0.52 * TMath::Exp(-0.042 * pt*pt) + 5.1/TMath::Power(pt+0.62,3));
  Double_t sz = TMath::Min(3.3, 1.10  + 0.39 * TMath::Exp(-0.027 * pt*pt) + 0.70 /TMath::Power(pt+0.223,3));

  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
  Double_t mf = event->GetMagneticField(); //Positive for ++ and negative for --

  if(mf<0.){ //field --
    meanZ = 0.077;
    if(charge>0)
      meanX =  TMath::Min(5.8, 0.2 + 0.7 * TMath::Exp(-0.019 * pt*pt) + 34./TMath::Power(pt+1.39,3));
    else
      meanX = -TMath::Min(5.8, 0.1 + 0.7 * TMath::Exp(-0.014 * pt*pt) + 30./TMath::Power(pt+1.36,3));
  }
  else{ //Field ++
    meanZ= 0.077;
    if(charge>0)
      meanX = -TMath::Min(5.8, 0.3 + 0.7 * TMath::Exp(-0.012 * pt*pt) + 35./TMath::Power(pt+1.43,3));
    else
      meanX =  TMath::Min(5.8, 0.2 + 0.6 * TMath::Exp(-0.014 * pt*pt) + 28./TMath::Power(pt+1.27,3));
  }

  Double_t rz=(dz-meanZ)/sz ;
  Double_t rx=(dx-meanX)/sx ;
  return TMath::Sqrt(rx*rx+rz*rz) ;

  
}
//____________________________________________________________________________
Double_t AliPHOSEmbeddingRun2::DecalibrateSignal(Double_t cellAmplitude,Int_t cellNumber){
  //Apply de-calibration inverse to the calibration, stored in OADB

  AliInfo(Form("fSignalECorrection = %e.",fSignalECorrection));

  //Apply overall energy correction  
  cellAmplitude *= fSignalECorrection;
    
  if(!fSignalCalibData){
    AliInfo("fSignalCalibData is not applied.");
    return cellAmplitude;
  }

  Int_t relId[4];
  AliPHOSGeometry * phosgeom = AliPHOSGeometry::GetInstance();
  phosgeom->AbsToRelNumbering(cellNumber,relId);
  if(relId[1]!=0) //CPV
    return  cellAmplitude ;  
  Int_t   module = relId[0];
  Int_t   column = relId[3];
  Int_t   row    = relId[2];
  Double_t c = fSignalCalibData->GetADCchannelEmc(module,column,row);
  
  if(c>0)
    return  cellAmplitude/c ; 
  else
    return  cellAmplitude ; 
}
//____________________________________________________________________________

