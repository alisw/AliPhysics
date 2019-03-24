// $Id$
//
//
//
//

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDUtils.h"
#include "AliESDInputHandler.h"
#include "AliCentrality.h"
#include "AliESDpid.h"
#include "AliKFParticle.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "TParticle.h"
#include "AliAODMCParticle.h"
#include "AliDataFile.h"


#include "AliESDtrackCuts.h"
#include "AliESDv0.h"
#include "AliV0vertexer.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliAODEvent.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliOADBContainer.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "TLorentzVector.h"


#include "AliAnalysisTaskEMCALPhoton.h"
#include "TFile.h"


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEMCALPhoton)

//________________________________________________________________________
AliAnalysisTaskEMCALPhoton::AliAnalysisTaskEMCALPhoton() : 
  AliAnalysisTaskSE(), 
  fTrCuts(0),
  fPrTrCuts(0),
  fSelTracks(0),
  fSelPrimTracks(0),
  fTracks(0),
  fPhotConvArray(0),
  fMyClusts(0),
  fMyAltClusts(0),
  fMyCells(0),
  fMyTracks(0),
  fMyMcParts(0),
  fHeader(0x0),
  fOADBContainer(0),
  fCaloClusters(0),
  fCaloClustersNew(0),
  fAODMCParticles(0),
  fVCells(0),
  fGeom(0x0),
  fTimeResTOF(0),
  fMipResponseTPC(0),
  fGeoName("EMCAL_COMPLETEV1"),
  fPeriod("LHC11d"),
  fIsTrain(0),
  fIsMC(0),
  fDebug(0),
  fRedoV0(1),
  fIsGrid(0),
  fClusThresh(2.0),
  fClusterizer(0),
  fCaloClustersName("EmcalClusterv2"),
  fESD(0),
  fAOD(0),
  fVev(0),
  fMCEvent(0),
  fStack(0x0),
  fOutputList(0),
  fTree(0),
  fMyMcIndex(0),
  fNV0sBefAndAftRerun(0),
  fConversionVtxXY(0),
  fInvMassV0(0),
  fInvMassV0KF(0),
  fInvMassV0SS(0),
  fDedxPAll(0)
{
  // Default constructor.
  for(Int_t i = 0; i < 12;    i++)  fGeomMatrix[i] =  0;
}

//________________________________________________________________________
AliAnalysisTaskEMCALPhoton::AliAnalysisTaskEMCALPhoton(const char *name) : 
  AliAnalysisTaskSE(name), 
  fTrCuts(0),
  fPrTrCuts(0),
  fSelTracks(0),
  fSelPrimTracks(0),
  fTracks(0),
  fPhotConvArray(0),
  fMyClusts(0),
  fMyAltClusts(0),
  fMyCells(0),
  fMyTracks(0),
  fMyMcParts(0),
  fHeader(0),
  fOADBContainer(0),
  fCaloClusters(0),
  fCaloClustersNew(0),
  fAODMCParticles(0),
  fVCells(0),
  fGeom(0x0),
  fTimeResTOF(0),
  fMipResponseTPC(0),
  fGeoName("EMCAL_COMPLETEV1"),
  fPeriod("LHC11c"),
  fIsTrain(0),
  fIsMC(0),
  fDebug(0),
  fRedoV0(1),
  fIsGrid(0),
  fClusThresh(2.),
  fClusterizer(0),
  fCaloClustersName("EmcalClusterv2"),
  fESD(0),
  fAOD(0),
  fVev(0),
  fMCEvent(0),
  fStack(0x0),
  fOutputList(0),
  fTree(0),
  fMyMcIndex(0),
  fNV0sBefAndAftRerun(0),
  fConversionVtxXY(0),
  fInvMassV0(0),
  fInvMassV0KF(0),
  fInvMassV0SS(0),
  fDedxPAll(0)
{
  // Constructor
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPhoton::UserCreateOutputObjects()
{
  // Create histograms, called once.
  if(this->fDebug)
    printf("::UserCreateOutputObjects() starting\n");

  fSelTracks = new TObjArray();
  
  fSelPrimTracks = new TObjArray();
  
  if (TClass::GetClass("AliPhotonHeaderObj"))
      TClass::GetClass("AliPhotonHeaderObj")->IgnoreTObjectStreamer();
  fHeader = new AliPhotonHeaderObj;

  if (TClass::GetClass("AliPhotonConvObj"))
      TClass::GetClass("AliPhotonConvObj")->IgnoreTObjectStreamer();
  fPhotConvArray = new TClonesArray("AliPhotonConvObj");
  
  if (TClass::GetClass("AliPhotonClusterObj"))
      TClass::GetClass("AliPhotonClusterObj")->IgnoreTObjectStreamer();
  fMyClusts = new TClonesArray("AliPhotonClusterObj");
  
  if (TClass::GetClass("AliPhotonCellObj"))
      TClass::GetClass("AliPhotonCellObj")->IgnoreTObjectStreamer();
  fMyCells = new TClonesArray("AliPhotonCellObj");
  
  if (TClass::GetClass("AliPhotonTrackObj"))
      TClass::GetClass("AliPhotonTrackObj")->IgnoreTObjectStreamer();
  fMyTracks = new TClonesArray("AliPhotonTrackObj");

  if (fClusterizer || fIsTrain){
    if(fClusterizer)
      fCaloClustersName = fClusterizer->GetNewClusterArrayName();
    else {
      if(fPeriod.Contains("10h") || fPeriod.Contains("11h"))
	fCaloClustersName = "EmcalClustersv1";
      else
	fCaloClustersName = "EmcalClustersv2";
    }
    if (TClass::GetClass("AliPhotonClusterObj"))
	TClass::GetClass("AliPhotonClusterObj")->IgnoreTObjectStreamer();
    fMyAltClusts = new TClonesArray("AliPhotonClusterObj");
  }
  cout<<fCaloClustersName<<endl;
  if(fIsMC){
    if (TClass::GetClass("AliPhotonMcPartObj"))
        TClass::GetClass("AliPhotonMcPartObj")->IgnoreTObjectStreamer();
    fMyMcParts = new TClonesArray("AliPhotonMcPartObj");
  }
 
  fCaloClusters = new TClonesArray();
    
  fOutputList = new TList();
  fOutputList->SetOwner();// Container cleans up all histos (avoids leaks in merging) 

  if( !fTree){
    TFile *f = OpenFile(2); 
    TDirectory::TContext context(f);
    if (f) {
      f->SetCompressionLevel(2);
      fTree = new TTree("photon_ana_out", "StandaloneTree");
      fTree->SetDirectory(f);
      if (fIsTrain) {
        fTree->SetAutoFlush(-2*1024*1024);
        fTree->SetAutoSave(0);
      } else {
        fTree->SetAutoFlush(-32*1024*1024);
        fTree->SetAutoSave(0);
      }
    
    fTree->Branch("header", &fHeader, 16*1024, 99);    
    fTree->Branch("conversions", &fPhotConvArray, 8*16*1024, 99);
    fTree->Branch("clusters", &fMyClusts, 8*16*1024, 99);
    if(fClusterizer || fIsTrain)
      fTree->Branch(fCaloClustersName, &fMyAltClusts, 8*16*1024, 99);
    fTree->Branch("cells", &fMyCells, 8*16*1024, 99);
    fTree->Branch("IsoTracks", &fMyTracks, 8*16*1024, 99);
    if(fIsMC)
      fTree->Branch("mcparts", &fMyMcParts, 8*16*1024, 99);
    }
  }
  //if(fIsGrid)fOutputList->Add(fTree);
  fGeom = AliEMCALGeometry::GetInstance(fGeoName);
  fOADBContainer = new AliOADBContainer("AliEMCALgeo");
  fOADBContainer->InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALlocal2master.root").data(),"AliEMCALgeo");
  
  
  fNV0sBefAndAftRerun = new TH2F("hNV0sBefAndAftRerun","check if the number of v0s change with rerun;old v0 n;new v0 n",50,0.5,50.5,50,0.5,50.5);
  fOutputList->Add(fNV0sBefAndAftRerun);
  
  fConversionVtxXY = new TH2F("hConversionVtxXY","x and y of conversion vertex candidates;x;y",1000,-100,100,1000,-100,100);
  fOutputList->Add(fConversionVtxXY);
  
  fInvMassV0 = new TH1F("hInvMassV0","v0->GetEffMass();v0->GetEffMass();dN/dM",400,0,4);
  fOutputList->Add(fInvMassV0);
  
  fInvMassV0KF = new TH1F("hInvMassV0KF","Inv. mass calculated from AliKFParticle made of V0 tracks;mass_{TrTr};dN/dM",400,0,4);
  fOutputList->Add(fInvMassV0KF);
      
  fInvMassV0SS = new TH1F("hInvMassV0SS","Inv. mass (same sign) calculated from AliKFParticle made of V0 tracks;mass_{TrTr};dN/dM",400,0,4);
  fOutputList->Add(fInvMassV0SS);
      
  fDedxPAll = new TH2F("hDedxPAll","dE/dx vs p (all selected tracks);p (GeV/c);dE/dx (a.u.)",400,0,40, 150,0,150);
  fOutputList->Add(fDedxPAll);
  
  
  PostData(1, fOutputList);
  PostData(2, fTree);

  if(this->fDebug)
    printf("::UserCreateOutputObjects() DONE!\n");

}

//________________________________________________________________________
void AliAnalysisTaskEMCALPhoton::UserExec(Option_t *) 
{
  // User exec, called once per event.

  Bool_t isSelected = kTRUE;
  if(fPeriod.Contains("11")){
    if(fPeriod.Contains("11a"))
      isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMC1);
    if(fPeriod.Contains("11c") ||fPeriod.Contains("11d") )
      isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMC7);
    if(fPeriod.Contains("11h") )
      isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kAny);//kEMCEGA);

  }
  if(fIsMC){
    isSelected = kTRUE;  
    if(this->fDebug)
      printf("+++Message+++: MC input files.\n");
  }
  if(!isSelected){
    if(this->fDebug)
      printf("+++Message+++: Event didn't pass the selection\n");
    return;
  }
  if(this->fDebug){
    printf("::UserExec(): event accepted\n");
    if(fIsMC)
      printf("\t in MC mode\n");
  }

  TTree *tree = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetTree();
  TFile *inpfile = (TFile*)tree->GetCurrentFile();

  fVev = (AliVEvent*)InputEvent();
  if (!fVev) {
    printf("ERROR: event not available\n");
    return;
  }
  Int_t   runnumber = InputEvent()->GetRunNumber() ;
  if(this->fDebug)
    printf("run number = %d\n",runnumber);
  fESD = dynamic_cast<AliESDEvent*>(fVev);
  if(!fESD){
    fAOD = dynamic_cast<AliAODEvent*>(fVev);
    if(!fAOD){
      printf("ERROR: Invalid type of event!!!\n");
      return;
    }
    else if(this->fDebug)
      printf("AOD EVENT!\n");
  }
  
  AliVVertex *pv = (AliVVertex*)fVev->GetPrimaryVertex();
  Bool_t pvStatus = kTRUE;
  if(fESD){
    AliESDVertex *esdv = (AliESDVertex*)fESD->GetPrimaryVertex();
    pvStatus = esdv->GetStatus();
  }

  if(!pv) {
    printf("Error: no primary vertex found!\n");
    return;
  }
  if(!pvStatus && this->fDebug)
    printf("bad pv status\n");
  if(TMath::Abs(pv->GetZ())>15)
    return;
  if(this->fDebug)
    printf("+++Message+++: event passed the vertex cut\n");

  if(fESD)
    fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("Tracks"));
  if(fAOD)
    fTracks = dynamic_cast<TClonesArray*>(fAOD->GetTracks());

  if(!fTracks){
    AliError("Track array in event is NULL!");
    if(this->fDebug)
      printf("returning due to not finding tracks!\n");
    return;
  }

  const Int_t Ntracks = fTracks->GetEntriesFast();
  // Track loop to fill a pT spectrum
  for (Int_t iTracks = 0;  iTracks < Ntracks; ++iTracks) {
    AliVTrack *track = (AliVTrack*)fTracks->At(iTracks);
    if (!track)
      continue;
    AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(track);
    AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*>(track);
    if (esdTrack){
      if (fTrCuts && fTrCuts->IsSelected(track)){
	fSelTracks->Add(track);
	fDedxPAll->Fill(track->P(), track->GetTPCsignal());
      }
      if (fPrTrCuts && fPrTrCuts->IsSelected(track))
	fSelPrimTracks->Add(track);
    }
    else if(aodTrack)
      fSelPrimTracks->Add(track);
  }//track loop 

    

  fHeader->fInputFileName  = inpfile->GetName();
  fHeader->fRunNumber =  runnumber;
  fHeader->fTrClassMask    = fVev->GetHeader()->GetTriggerMask();
  fHeader->fTrCluster      = fVev->GetHeader()->GetTriggerCluster();
  AliCentrality *cent = InputEvent()->GetCentrality();
  fHeader->fV0Cent    = cent->GetCentralityPercentileUnchecked("V0M");
  fHeader->fCl1Cent   = cent->GetCentralityPercentileUnchecked("CL1");
  fHeader->fTrCent    = cent->GetCentralityPercentileUnchecked("TRK");

  TObjArray *matEMCAL=(TObjArray*)fOADBContainer->GetObject(runnumber,"EmcalMatrices");
  for(Int_t mod=0; mod < (fGeom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
    if(fGeoName=="EMCAL_FIRSTYEARV1" && mod>3)
      break;
    /*if(fESD)
      fGeom->SetMisalMatrix(fESD->GetEMCALMatrix(mod), mod);
      else*/
    // if(event->GetEMCALMatrix(mod))
    fGeomMatrix[mod] = (TGeoHMatrix*) matEMCAL->At(mod);
    fGeom->SetMisalMatrix(fGeomMatrix[mod] , mod);
  }
  Int_t trackMult = 0;
  if(fESD){
    AliESDtrackCuts *fTrackCuts = new AliESDtrackCuts();
    trackMult = fTrackCuts->GetReferenceMultiplicity(fESD);//kTrackletsITSTPC ,0.5); 
    if(this->fDebug)
      printf("ESD Track mult= %d\n",trackMult);
  }
  else if(fAOD){
    trackMult = Ntracks;
    if(this->fDebug)
      printf("AOD Track mult= %d\n",trackMult);
  }
  if(fESD){
    TList *l = fESD->GetList();
    fCaloClusters =  dynamic_cast<TClonesArray*>(l->FindObject("CaloClusters"));
    fVCells = fESD->GetEMCALCells();
    fHeader->fNCells = fVCells->GetNumberOfCells();
    if(this->fDebug)
      printf("ESD cluster mult= %d\n",fCaloClusters->GetEntriesFast());
  }
  else if(fAOD){
    fCaloClusters = dynamic_cast<TClonesArray*>(fAOD->GetCaloClusters());
    fVCells = fAOD->GetEMCALCells();
    fHeader->fNCells = fVCells->GetNumberOfCells();
    if(this->fDebug)
      printf("AOD cluster mult= %d\n",fCaloClusters->GetEntriesFast());
  }
    

  fHeader->fNClus = fCaloClusters->GetEntriesFast();
  fHeader->fTrackMult = trackMult;

  fMCEvent = MCEvent();
  if(fMCEvent){
    fStack = (AliStack*)fMCEvent->Stack();
    if(this->fStack)
      fHeader->fNMcParts = this->fStack->GetNtrack();
    else{
      printf("Stack not found\n");
      fHeader->fNMcParts = 0;
      fAODMCParticles = (TClonesArray*)fVev->FindListObject("mcparticles");
    }
    if(fAODMCParticles)
      fHeader->fNMcParts = fAODMCParticles->GetEntriesFast();
  }
  else{
    if(fIsMC){
      printf("ERROR: MC Event not available, returning...\n");
      return;
    }
  }

  
  FindConversions();
  if(this->fDebug)
    printf("FindConversions done\n");
  FillMyCells();
  if(this->fDebug)
    printf("FillMyCells done\n");
  FillMyClusters();
  if(this->fDebug)
    printf("FillMyClusters done\n");
  if(fCaloClustersNew)
    FillMyAltClusters();
  FillIsoTracks();
  if(fIsMC)
    GetMcParts();

  if(this->fDebug && fIsMC)
    printf("fMyMcParts nentries=%d",fMyMcParts->GetEntries());
  
  fTree->Fill();
  fTracks->Clear();
  fSelTracks->Clear();
  fSelPrimTracks->Clear();
  fPhotConvArray->Clear();
  fMyClusts->Clear();
  if(fMyAltClusts)
    fMyAltClusts->Clear();
  fMyCells->Clear();
  fMyTracks->Clear();
  if(fMyMcParts)
    fMyMcParts->Clear();
  fMyMcIndex = 0;
  fCaloClusters->Clear();
  if(fCaloClustersNew)
    fCaloClustersNew->Clear();
  if(fVCells)
    fVCells = 0;
  // Post output data.
  PostData(1, fOutputList);
  PostData(2, fTree);
}      

//________________________________________________________________________
void AliAnalysisTaskEMCALPhoton::FindConversions() //WARNING: not ready to use with AODs
{
  // Find conversion.
  if(!fESD)//not working with AODs yet
    return;
  if(this->fDebug)
    printf("::FindConversions() starting\n");
  if(!fTrCuts)
    return;
  Int_t iconv = 0;
  Int_t nV0Orig = 0;
  if(fESD)
    nV0Orig = fESD->GetNumberOfV0s();
  if(fAOD)
    nV0Orig = fAOD->GetNumberOfV0s();
  Int_t nV0s = nV0Orig;
  Int_t ntracks = fSelTracks->GetEntriesFast();
  if(fRedoV0 && !fAOD){
    fESD->ResetV0s();
    AliV0vertexer lV0vtxer;
    lV0vtxer.Tracks2V0vertices(fESD);
    nV0s = fESD->GetNumberOfV0s();
  }
  fNV0sBefAndAftRerun->Fill(nV0Orig, nV0s);
  for(Int_t iv0=0; iv0<nV0s; iv0++){
    AliESDv0 * v0 = fESD->GetV0(iv0);
    if(!v0)
      continue;
    Double_t vpos[3];
    fInvMassV0->Fill(v0->GetEffMass());
    v0->GetXYZ(vpos[0], vpos[1], vpos[2]);
    Int_t ipos = v0->GetPindex();
    Int_t ineg = v0->GetNindex();
    if(ipos<0 || ipos> ntracks)
      continue;
    if(ineg<0 || ineg> ntracks)
      continue;
    AliESDtrack *pos = static_cast<AliESDtrack*>(fSelTracks->At(ipos));
    AliESDtrack *neg = static_cast<AliESDtrack*>(fSelTracks->At(ineg));
    if(!pos || !neg)
      continue;
    /*if(pos->GetTPCsignal()<65 || neg->GetTPCsignal()<65)
      continue;*/
    const AliExternalTrackParam * paramPos = v0->GetParamP() ;
    const AliExternalTrackParam * paramNeg = v0->GetParamN() ;
    if(!paramPos || !paramNeg)
      continue;
    if(pos->GetSign() <0){//change tracks
      pos=neg ;
      neg=fESD->GetTrack(v0->GetPindex()) ;
      paramPos=paramNeg ;
      paramNeg=v0->GetParamP() ;
    }
    AliKFParticle negKF(*paramNeg,11);
    AliKFParticle posKF(*paramPos,-11);
    AliKFParticle photKF(negKF,posKF) ;
   
    if( neg->GetKinkIndex(0) > 0 || pos->GetKinkIndex(0) > 0) 
      continue ;
    
    if(pos->GetSign() == neg->GetSign()){ 
      fInvMassV0SS->Fill(photKF.GetMass());
      continue;
    }
    fInvMassV0KF->Fill(photKF.GetMass());
    TLorentzVector photLV(photKF.GetPx(), photKF.GetPy(),photKF.GetPz(), photKF.GetE()); 
    if(photLV.M()>0.05 || photLV.M()<0)
      continue;
    fConversionVtxXY->Fill(vpos[0], vpos[1]);//negKF.GetX(), negKF.GetY());
    Double_t convPhi = TMath::ATan2(photLV.Py(),photLV.Px());
    if(convPhi<0)
      convPhi+=TMath::TwoPi();
    TVector3 vecpos(vpos);
    Double_t v0Phi = 0;
    if(vecpos.Perp()>0)
      vecpos.Phi();
    if(v0Phi<0)
      v0Phi+=TMath::TwoPi();
    AliPhotonConvObj *myconv = static_cast<AliPhotonConvObj*>(fPhotConvArray->New(iconv++));
    myconv->fPt = photLV.Pt();
    myconv->fEta = photLV.Eta();
    myconv->fPhi = convPhi;
    myconv->fVR = vecpos.Perp();
    if(vecpos.Perp()>0)
      myconv->fVEta = vecpos.Eta();
    myconv->fVPhi = v0Phi;
    myconv->fMass = photLV.M();
    myconv->fMcLabel = -3; //WARNING: include the correct labeling
    //negative daughter
   myconv->fNegPt      =  negKF.GetPt();
   myconv->fNegEta     =  negKF.GetEta();
   Double_t trackPhi   =  negKF.GetPhi();
   if(trackPhi<0)
     trackPhi+=TMath::TwoPi();
   myconv->fNegPhi     =  trackPhi;
   myconv->fNegDedx    =  neg->GetTPCsignal();
   myconv->fNegMcLabel =  neg->GetLabel();
    //negative daughter
   myconv->fPosPt      =  posKF.GetPt();
   myconv->fPosEta     =  posKF.GetEta();
   trackPhi            =  posKF.GetPhi();
   if(trackPhi<0)
     trackPhi+=TMath::TwoPi();
   myconv->fPosPhi     =  trackPhi;
   myconv->fPosDedx    =  pos->GetTPCsignal();
   myconv->fPosMcLabel =  pos->GetLabel();
  }
  if(this->fDebug)
    printf("::FindConversions() returning...\n\n");
  return;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPhoton::FillMyCells() 
{
  // Fill cells.
  if(this->fDebug)
    printf("::FillMyCells() starting\n");

  if (!fVCells)
    return;
  Int_t ncells = fVCells->GetNumberOfCells();
  Int_t mcel = 0;//, maxcelid=-1;
  Double_t maxcellE = 0;//, maxcellEta=0, maxcellPhi=0;
  for(Int_t icell = 0; icell<ncells; icell++){
    Int_t absID = TMath::Abs(fVCells->GetCellNumber(icell));
    AliPhotonCellObj *mycell = static_cast<AliPhotonCellObj*>(fMyCells->New(mcel++));
    Float_t eta=-1, phi=-1;
    if(!fGeom){
      std::cout<<"+++fGeom not found! MyCells branch will not be filled for this event!+++\n";
      return;
    }
    if(!fGeom)
      return;
    /*if(!fIsMC)*/fGeom->EtaPhiFromIndex(absID,eta,phi);
    if(maxcellE<fVCells->GetCellAmplitude(absID)){
      maxcellE = fVCells->GetCellAmplitude(absID);
      /*maxcellEta = eta;
      maxcellPhi = phi;
      maxcelid = absID;*/
    }
    Float_t theta = 2*TMath::ATan(TMath::Exp(-eta));
    mycell->fAbsID = absID;
    mycell->fE = fVCells->GetCellAmplitude(absID);
    mycell->fEt = fVCells->GetCellAmplitude(absID)*TMath::Sin(theta);
    mycell->fEta = eta;
    mycell->fPhi = phi;
    mycell->fTime = fVCells->GetCellTime(absID);
  }
  if(this->fDebug)
    printf("::FillMyCells() returning...\n\n");
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPhoton::FillMyClusters() 
{
  // Fill clusters.
  if(this->fDebug)
    printf("::FillMyClusters() starting\n");

  if (!fCaloClusters){
    printf("CaloClusters is empty!\n");
    return;
  }
  Int_t nclus = fCaloClusters->GetEntries();
  if(0==nclus)
    printf("CaloClusters has ZERO entries\n");
  Int_t mcl = 0, maxcelid=-1;
  Double_t maxcellE=0, maxcellEtac=0,maxcellPhic=0;
  for(Int_t ic=0; ic < nclus; ic++){
    AliVCluster *clus = static_cast<AliVCluster*>(fCaloClusters->At(ic));
    if(!clus)
      continue;
    if(!clus->IsEMCAL())
      continue;
    if(clus->E() < fClusThresh)
      continue;
    if(fDebug)
      printf("cluster %d survived\n", ic);
    Float_t pos[3];
    clus->GetPosition(pos);
    TVector3 cpos(pos);
    TString cellsAbsID;
    AliPhotonClusterObj *myclus = static_cast<AliPhotonClusterObj*>(fMyClusts->New(mcl++));
    myclus->fE       = clus->E();
    myclus->fEt      = clus->E()*TMath::Sin(cpos.Theta());
    myclus->fR       = cpos.Perp();
    myclus->fEta     = cpos.Eta();
    myclus->fPhi     = cpos.Phi();
    if(cpos.Phi()<0){
      myclus->fPhi+=TMath::TwoPi();
    }
    myclus->fN       = clus->GetNCells();
    Short_t  id = -1;
    myclus->fEmax    = GetMaxCellEnergy( clus, id); 
    myclus->fIdmax   = id;
    if(maxcellE <  myclus->fEmax){
      maxcellE =  myclus->fEmax;
      maxcelid = id;
      maxcellEtac = cpos.Eta();
      maxcellPhic = cpos.Phi();
    }
    myclus->fTmax    = fVCells->GetCellTime(id);
    myclus->fEcross  = GetCrossEnergy( clus, id);
    myclus->fDisp    = clus->GetDispersion();
    myclus->fM20     = clus->GetM20();
    myclus->fM02     = clus->GetM02();
    myclus->fTrDEta  = clus->GetTrackDz();
    myclus->fTrDPhi  = clus->GetTrackDx();
    myclus->fTrIso01 = GetTrackIsolation( cpos.Eta(), cpos.Phi(), 0.1, 0.);
    myclus->fTrIso02 = GetTrackIsolation( cpos.Eta(), cpos.Phi(), 0.2, 0.);
    myclus->fTrIso03 = GetTrackIsolation( cpos.Eta(), cpos.Phi(), 0.3, 0.);
    myclus->fTrIso04 = GetTrackIsolation( cpos.Eta(), cpos.Phi(), 0.4, 0.);
    myclus->fTrPhiBand01 = GetPhiBandEt( cpos.Eta(), cpos.Phi(), 0.1, 0.);
    myclus->fTrPhiBand02 = GetPhiBandEt( cpos.Eta(), cpos.Phi(), 0.2, 0.);
    myclus->fTrPhiBand03 = GetPhiBandEt( cpos.Eta(), cpos.Phi(), 0.3, 0.);
    myclus->fTrPhiBand04 = GetPhiBandEt( cpos.Eta(), cpos.Phi(), 0.4, 0.);
    for(Int_t icell=0;icell<myclus->fN;icell++){
      Int_t absID = clus->GetCellAbsId(icell);
      cellsAbsID.Append(Form("%d",absID));
      cellsAbsID.Append(";");
    }
    myclus->fCellsAbsId = cellsAbsID;
    myclus->fMcLabel = clus->GetLabel(); 
    Int_t matchIndex = clus->GetTrackMatchedIndex();
    if(matchIndex<0 || matchIndex>fVev->GetNumberOfTracks()){
      myclus->fTrEp = -1;
      continue;
    }
    AliVTrack* track = static_cast<AliVTrack*>(fVev->GetTrack(matchIndex));
    if(!track)
      continue;
    if(fESD){
      if(!fPrTrCuts)
	continue;
      if(!fPrTrCuts->IsSelected(track))
	continue;
    }
    
    myclus->fTrEp = clus->E()/track->P();
    myclus->fTrDedx = track->GetTPCsignal();
  }
  if(this->fDebug){
    printf(" ---===+++ Max Cell among clusters: id=%d, E=%1.2f, eta-clus=%1.2f, phi-clus=%1.2f\n",maxcelid,maxcellE,maxcellEtac,maxcellPhic);
    printf("::FillMyClusters() returning...\n\n");
  }
  
}
//________________________________________________________________________
void AliAnalysisTaskEMCALPhoton::FillMyAltClusters() 
{
  // Fill clusters.
  if(this->fDebug)
    printf("::FillMyAltClusters() starting\n");

  if(!fCaloClustersNew)
    return;
  Int_t nclus = fCaloClustersNew->GetEntries();
  Int_t mcl = 0;
  for(Int_t ic=0; ic < nclus; ic++){
    AliVCluster *clus = static_cast<AliVCluster*>(fCaloClustersNew->At(ic));
    if(!clus)
      continue;
    if(!clus->IsEMCAL())
      continue;
    if(clus->E() < fClusThresh)
      continue;
    Float_t pos[3];
    clus->GetPosition(pos);
    TVector3 cpos(pos);
    TString cellsAbsID;
    AliPhotonClusterObj *myclus = static_cast<AliPhotonClusterObj*>(fMyAltClusts->New(mcl++));
    myclus->fE       = clus->E();
    myclus->fEt      = clus->E()*TMath::Sin(cpos.Theta());
    myclus->fR       = cpos.Perp();
    myclus->fEta     = cpos.Eta();
    myclus->fPhi     = cpos.Phi();
    if(cpos.Phi()<0){
      myclus->fPhi+=TMath::TwoPi();
    }
    myclus->fN       = clus->GetNCells();
    myclus->fDisp    = clus->GetDispersion();
    myclus->fM20     = clus->GetM20();
    myclus->fM02     = clus->GetM02();
    myclus->fTrDEta  = clus->GetTrackDz();
    myclus->fTrDPhi  = clus->GetTrackDx();
    myclus->fTrIso01 = GetTrackIsolation( cpos.Eta(), cpos.Phi(), 0.1, 0.);
    myclus->fTrIso02 = GetTrackIsolation( cpos.Eta(), cpos.Phi(), 0.2, 0.);
    myclus->fTrIso03 = GetTrackIsolation( cpos.Eta(), cpos.Phi(), 0.3, 0.);
    myclus->fTrIso04 = GetTrackIsolation( cpos.Eta(), cpos.Phi(), 0.4, 0.);
    for(Int_t icell=0;icell<myclus->fN;icell++){
      Int_t absID = clus->GetCellAbsId(icell);
      cellsAbsID.Append(Form("%d",absID));
      cellsAbsID.Append(";");
    }
    myclus->fCellsAbsId = cellsAbsID;
    myclus->fMcLabel = clus->GetLabel(); 
    Int_t matchIndex = clus->GetTrackMatchedIndex();
    if(matchIndex<0 || matchIndex>fVev->GetNumberOfTracks()){
      myclus->fTrEp = -1;
      continue;
    }
    AliVTrack* track = static_cast<AliVTrack*>(fVev->GetTrack(matchIndex));
    if(!track){
      myclus->fTrEp = -1;
      continue;
    }
    if(!fPrTrCuts){
      myclus->fTrEp = -1;
      continue;
    }
    if(!fPrTrCuts->IsSelected(track)){
      myclus->fTrEp = -1;
      continue;
    }
    myclus->fTrEp = clus->E()/track->P();
  }
  if(this->fDebug)
    printf("::FillMyAltClusters() returning...\n\n");
  
}
//________________________________________________________________________
void  AliAnalysisTaskEMCALPhoton::FillIsoTracks()
{
  // Fill high pt tracks.
  if(this->fDebug)
    printf("::FillIsoTracks() starting\n");

  if(!fSelPrimTracks)
    return;
  Int_t ntracks = fSelPrimTracks->GetEntries();
  Int_t imtr = 0;
  for(Int_t it=0;it<ntracks; it++){
    AliVTrack *track = static_cast<AliVTrack*>(fSelPrimTracks->At(it));
    if(!track)
      continue;
    /*if(track->Phi()<1.0 || track->Phi()>3.55)
      continue;*/
    if(TMath::Abs(track->Eta())>1.1)
      continue;
    /*if(track->Pt()<3)
      continue;*/
    AliPhotonTrackObj *mtr = static_cast<AliPhotonTrackObj*>(fMyTracks->New(imtr++));
    mtr->fPt = track->Pt();
    mtr->fEta = track->Eta();
    mtr->fPhi = track->Phi();
    mtr->fCharge = track->Charge();
    mtr->fDedx = track->GetTPCsignal();
    mtr->fMcLabel = track->GetLabel();
  }
  if(this->fDebug)
    printf("::FillIsoTracks() returning...\n\n");
}

//________________________________________________________________________
void  AliAnalysisTaskEMCALPhoton::GetMcParts()
{
   // Get MC particles.
  if(this->fDebug)
    printf("::GetMcParts() starting\n");

  if (!this->fStack && !fAODMCParticles)
    return;
  if(this->fDebug)
    printf("either stack or aodmcpaticles exists\n");
  const AliVVertex *evtVtx = 0;
  if(this->fStack)
     evtVtx = fMCEvent->GetPrimaryVertex();
  else
    printf("no such thing as mc vvtx\n");
  /*if (!evtVtx)
    return;*/
  if(this->fDebug)
    printf("mc vvertex ok\n");
  Int_t nTracks = 0;
  if(this->fStack)
    nTracks = this->fStack->GetNtrack();
  else if(fAODMCParticles)
    nTracks = fAODMCParticles->GetEntriesFast();
  TParticle *mcP = 0;
  AliAODMCParticle *amcP = 0;
  if(this->fDebug)
    printf("loop in the %d mc particles starting\n",nTracks);
  for (Int_t iTrack = 0; iTrack<nTracks; ++iTrack) {
    if(this->fStack)
      mcP = dynamic_cast<TParticle*>(this->fStack->Particle(iTrack));
    if(fAODMCParticles)
      amcP = dynamic_cast<AliAODMCParticle*>(fAODMCParticles->At(iTrack));

    // primary particle
    Double_t dR = 0;
    if(mcP){
      dR = TMath::Sqrt((mcP->Vx()-evtVtx->GetX())*(mcP->Vx()-evtVtx->GetX()) + 
		       (mcP->Vy()-evtVtx->GetY())*(mcP->Vy()-evtVtx->GetY()) +
		       (mcP->Vz()-evtVtx->GetZ())*(mcP->Vz()-evtVtx->GetZ()));
    }
    if((dR > 0.5))
      continue;
    
    
    // kinematic cuts
    Double_t pt = 0;
    Double_t eta = 0;
    Double_t phi  = 0;
    Int_t mother = -1;
    Int_t pdgcode = 0;
     if(mcP){ 
      pt = mcP->Pt() ;
      eta = mcP->Eta();
      phi = mcP->Phi();
      mother = mcP->GetMother(0);
      pdgcode = mcP->GetPdgCode();
     } else { 
       if(amcP){
	 pt = amcP->Pt();
	 eta = amcP->Eta();
	 phi = amcP->Phi();
	 mother = amcP->GetMother();
	 pdgcode = amcP->GetPdgCode();
       } else
       continue;
     }
    if (pt<0.5)
      continue;
    
    if (TMath::Abs(eta)>0.7)
      continue;

    if (phi<1.0||phi>3.3)
      continue;
    
    if (mother!=6 && mother!=7 )
      continue;


    // pion or eta meson or direct photon
    if(pdgcode == 111) {
    } else if(pdgcode == 221) {
    } else if(pdgcode == 22 ) {
    } else {
      continue;
    }

    FillMcPart( fMyMcIndex++, iTrack);
  }
  if(this->fDebug)
    printf("::GetMcParts() returning...\n\n");
}

//________________________________________________________________________
void  AliAnalysisTaskEMCALPhoton::FillMcPart(  Int_t itrack, Int_t label)
{
  // Fill MC particles.
  Int_t nTracks = 0;
  TParticle *mcP = 0;
  AliAODMCParticle *amcP= 0;
  if(this->fStack){
    nTracks = this->fStack->GetNtrack();
    mcP = dynamic_cast<TParticle*>(this->fStack->Particle(itrack));
  }
  else if(fAODMCParticles){
    nTracks = fAODMCParticles->GetEntriesFast();
    amcP = dynamic_cast<AliAODMCParticle*>(fAODMCParticles->At(itrack));
 }
  if(this->fDebug)
    printf("\t::FillMcParts() starting with label %d\n", itrack);
   TVector3 vmcv;
   if(mcP)
     vmcv.SetXYZ(mcP->Vx(),mcP->Vy(), mcP->Vz());
   else { 
     if(amcP)
       vmcv.SetXYZ(amcP->Xv(),amcP->Yv(), amcP->Zv());
     else
       return;
   }
 
  AliPhotonMcPartObj *mcp = static_cast<AliPhotonMcPartObj*>(fMyMcParts->New(itrack));
  mcp->fLabel = label ;
  if(mcP){
    mcp->fPdg = mcP->GetPdgCode() ;
    mcp->fPt = mcP->Pt() ;
    mcp->fEta = mcP->Eta() ;
    mcp->fPhi = mcP->Phi() ;
    mcp->fMother = mcP->GetMother(0) ;
    mcp->fFirstD = mcP->GetFirstDaughter() ;
    mcp->fLastD = mcP->GetLastDaughter() ;
    mcp->fStatus = mcP->GetStatusCode();
  }
  if(amcP){
    mcp->fPdg = amcP->GetPdgCode() ;
    mcp->fPt = amcP->Pt() ;
    mcp->fEta = amcP->Eta() ;
    mcp->fPhi = amcP->Phi() ;
    mcp->fMother = amcP->GetMother() ;
    mcp->fFirstD = amcP->GetDaughterLabel(0) ;
    mcp->fLastD = amcP->GetDaughterLabel(amcP->GetNDaughters()-1) ;
    mcp->fStatus = amcP->GetStatus();
  }
  mcp->fVR = vmcv.Perp();
  if(vmcv.Perp()>0){
    mcp->fVEta = vmcv.Eta() ;
    mcp->fVPhi = vmcv.Phi() ;
  }
  if(itrack == 8){
    mcp->fIso = AliAnalysisTaskEMCALPhoton::GetMcIsolation( itrack, 0.4 , 0.2);
    mcp->fIso3 = AliAnalysisTaskEMCALPhoton::GetMcIsolation( itrack, 0.3 , 0.2);
    }
  if(this->fDebug)
    printf("\t::FillMcParts(): label=%d, pdg=%d, pt=%1.1f, mom=%d, 1stD=%d,last=%d\n\t::FillMcParts() returning...\n\n", mcp->fLabel,mcp->fPdg,mcp->fPt,mcp->fMother,mcp->fFirstD,mcp->fLastD);
  for(Int_t id=mcp->fFirstD; id < mcp->fLastD; id++){
    if(id<=mcp->fMother)
      continue;
    if(id<0 || id>nTracks)
      continue;
    FillMcPart( fMyMcIndex++, id);
  }
  
}
//________________________________________________________________________                                                                                                                                   
Double_t AliAnalysisTaskEMCALPhoton::GetMcIsolation( Int_t itrack, Double_t radius, Double_t ptcut) const
{
  if(this->fDebug){
    printf("\t\t::GetMcIsolation() starting\n");
    //printf("\t\t   incoming particle: PDG = %d, itrack=%d;\n",mcP->GetPdgCode(),itrack);
  }
  if (!this->fStack && !this->fAODMCParticles && this->fDebug){
    printf("\t\t\tNo MC stack/array!\n");
    return -1;
  }
  if(itrack<6 || itrack>8)
    return -1;
  if(this->fDebug)
    printf("\t\t\tparticle of interest selected\n");
  TParticle *mcP = 0;
  AliAODMCParticle *amcP = 0;
  Int_t pdgcode = 0;
  Float_t eta = 0;
  Float_t phi = 0;
  Double_t sumpt=0;
  Float_t dR;
  Int_t nparts =  0;
  if(this->fStack){
    nparts = fStack->GetNtrack();
    mcP = dynamic_cast<TParticle*>(this->fStack->Particle(itrack));  
    eta = mcP->Eta();
    phi = mcP->Phi();
    pdgcode = mcP->GetPdgCode();
  }
  if(this->fAODMCParticles){
    nparts = fAODMCParticles->GetEntriesFast();
    amcP = dynamic_cast<AliAODMCParticle*>(this->fAODMCParticles->At(itrack));
    if(amcP){
      eta = amcP->Eta();
      phi = amcP->Phi();
      pdgcode = amcP->GetPdgCode();
    }
  }
  if(pdgcode!=22)
    return -1;
  TParticle *mcisop = 0;
  AliAODMCParticle *amcisop = 0;
  for(Int_t ip = 0; ip<nparts; ip++){
    if(ip==itrack)
      continue;
    if(this->fStack)
      mcisop = dynamic_cast<TParticle*>(this->fStack->Particle(ip));
    if(fAODMCParticles)
      amcisop = dynamic_cast<AliAODMCParticle*>(fAODMCParticles->At(ip));
    Int_t status = 0;
    Int_t mother = 0;
    Float_t pt = 0;
    Float_t isophi = -99;
    Float_t isoeta = -99;
    TVector3 vmcv;
    if(mcisop){
      status = mcisop->GetStatusCode();
      mother = mcisop->GetMother(0);
      pt = mcisop->Pt();
      isophi = mcisop->Phi();
      isoeta = mcisop->Eta();
      vmcv.SetXYZ(mcisop->Vx(),mcisop->Vy(), mcisop->Vz());
    }
    else {
      if(amcisop){
	status = amcisop->GetStatus();
	mother = amcisop->GetMother();
	pt = amcisop->Pt();
	isophi = amcisop->Phi();
	isoeta = amcisop->Eta();
	vmcv.SetXYZ(amcisop->Xv(),amcisop->Yv(), amcisop->Zv());
      }
      else
	continue;
    }
    if(status!=1)
      continue;
    if(mother == itrack)
      continue;
    if(pt<ptcut)
      continue;
    if(vmcv.Perp()>1)
      continue;
    dR = TMath::Sqrt((phi-isophi)*(phi-isophi)+(eta-isoeta)*(eta-isoeta));
    if(dR>radius)
      continue;
    sumpt += pt;
  }
  if(this->fDebug)
    printf("\t\t::GetMcIsolation() returning value %f ...\n\n",sumpt);
  return sumpt;
 }

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPhoton::GetTrackIsolation(Double_t cEta, Double_t cPhi, Double_t radius, Double_t pt) const
{
  // Compute isolation based on tracks.
  if(this->fDebug)
    printf("\t::GetTrackIsolation() starting\n");
   
  Double_t trkIsolation = 0;
  Double_t rad2 = radius*radius;
  Int_t ntrks = fSelPrimTracks->GetEntries();
  for(Int_t j = 0; j<ntrks; ++j) {
    AliVTrack *track = static_cast<AliVTrack*>(fSelPrimTracks->At(j));
    if (!track)
      continue;
    if (track->Pt()<pt)
      continue;
    
    Float_t eta = track->Eta();
    Float_t phi = track->Phi();
    Double_t phidiff = TVector2::Phi_mpi_pi(phi-cPhi);
    Double_t dist = (eta-cEta)*(eta-cEta)+phidiff*phidiff;
    if(dist>rad2)
      continue;
    trkIsolation += track->Pt();
  } 
  if(this->fDebug)
    printf("\t::GetTrackIsolation() returning\n\n");
  return trkIsolation;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPhoton::GetPhiBandEt(Double_t eta, Double_t phi, Double_t R, Double_t minpt) const
{
  // Get phi band.

  if(!fSelPrimTracks)
    return 0;
  Int_t ntracks = fSelPrimTracks->GetEntries();
  Double_t loweta = eta - R;
  Double_t upeta = eta + R;
  Double_t upphi = phi + R;
  Double_t lowphi = phi - R;
  Double_t et = 0;
  for(Int_t itr=0; itr<ntracks; itr++){
    AliVTrack *track = static_cast<AliVTrack*>(fSelPrimTracks->At(itr));
    if(!track)
      continue;
    if(track->Pt()<minpt)
      continue;
    if((track->Eta() < upeta) && (track->Eta() > loweta))
      continue;
    if((track->Phi() > upphi) || (track->Phi() < lowphi))
      continue;
    et+=track->Pt();
  }
  return et;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPhoton::GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax)
{
  // Calculate the energy of cross cells around the leading cell.
  if(!fVCells)
    return 0;

  if (!fGeom)
    return 0;

  Int_t iSupMod = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  Int_t iphis   = -1;
  Int_t ietas   = -1;

  Double_t crossEnergy = 0;

  fGeom->GetCellIndex(idmax,iSupMod,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphis,ietas);

  Int_t ncells = cluster->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Int_t cellAbsId = cluster->GetCellAbsId(i);
    fGeom->GetCellIndex(cellAbsId,iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);
    Int_t aphidiff = TMath::Abs(iphi-iphis);
    if (aphidiff>1)
      continue;
    Int_t aetadiff = TMath::Abs(ieta-ietas);
    if (aetadiff>1)
      continue;
    if ( (aphidiff==1 && aetadiff==0) ||
	(aphidiff==0 && aetadiff==1) ) {
      crossEnergy += fVCells->GetCellAmplitude(cellAbsId);
    }
  }

  return crossEnergy;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALPhoton ::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
{
  // Get maximum energy of attached cell.

  id = -1;
  if(!fVCells)
    return 0;

  Double_t maxe = 0;
  Int_t ncells = cluster->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Double_t e = fVCells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
    if (e>maxe) {
      maxe = e;
      id   = cluster->GetCellAbsId(i);
    }
  }
  return maxe;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPhoton::Terminate(Option_t *) 
{
  // Called once at the end of the query
/*  if(fIsGrid)
    return;*/
  if (fTree) {
    printf("***tree %s being saved***\n",fTree->GetName());
    TFile *f = OpenFile(2);
    TDirectory::TContext context(f);
    if (f) 
      fTree->Write();
  }
}
