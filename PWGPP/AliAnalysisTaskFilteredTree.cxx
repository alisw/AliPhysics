/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
   Class to process/filter reconstruction information from ESD, ESD friends, MC and provide them for later reprocessing
   Filtering schema - low pt part is downscaled - to have flat pt specatra of selected topologies (tracks and V0s)
   Downscaling schema is controlled by downscaling factors
   Usage: 
     1.) Filtering on Lego train
     2.) expert QA for tracking (resolution efficnecy)
     3.) pt reoslution studies using V0s 
     4.) dEdx calibration using V0s
     5.) pt resolution and dEdx studies using cosmic
     +
     6.) Info used for later raw data OFFLINE triggering  (highPt, V0, laser, cosmic, high dEdx)

   Exported trees (with full objects and dereived variables):
   1.) "highPt"     - filtered trees with esd tracks, derived variables(propagated tracks), optional MC info +optional space points
   2.) "V0s" -      - filtered trees with selected V0s (rough KF chi2 cut), KF particle and corresponding esd tracks + optionla space points
   3.) "Laser"      - dump laser tracks with space points if exests 
   4.) "CosmicTree" - cosmic track candidate (random or triggered) + esdTracks(up/down)+ optional points
   5.) "dEdx"       - tree with high dEdx tpc tracks
*/

#include "iostream"
#include "TSystem.h"
#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include "TChain.h"
#include "TTreeStream.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TList.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TRandom3.h"

#include "AliHeader.h"  
#include "AliGenEventHeader.h"  
#include "AliInputEventHandler.h"  
#include "AliStack.h"  
#include "AliTrackReference.h"  
#include "AliTrackPointArray.h"
#include "AliSysInfo.h"

#include "AliPhysicsSelection.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliTracker.h"
#include "AliVTrack.h"
#include "AliGeomManager.h"

#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "AliMultiplicity.h"

#include "AliESDtrackCuts.h"
#include "AliMCEventHandler.h"
#include "AliFilteredTreeEventCuts.h"
#include "AliFilteredTreeAcceptanceCuts.h"

#include "AliAnalysisTaskFilteredTree.h"
#include "AliKFParticle.h"
#include "AliESDv0.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "TVectorD.h"
#include "TStatToolkit.h"
using namespace std;

ClassImp(AliAnalysisTaskFilteredTree)

  //_____________________________________________________________________________
  AliAnalysisTaskFilteredTree::AliAnalysisTaskFilteredTree(const char *name) 
  : AliAnalysisTaskSE(name)
  , fESD(0)
  , fMC(0)
  , fESDfriend(0)
  , fOutput(0)
  , fPitList(0)
  , fUseMCInfo(kFALSE)
  , fUseESDfriends(kFALSE)
  , fReducePileUp(kTRUE)
  , fFillTree(kTRUE)
  , fFilteredTreeEventCuts(0)
  , fFilteredTreeAcceptanceCuts(0)
  , fFilteredTreeRecAcceptanceCuts(0)
  , fEsdTrackCuts(0)
  , fTrigger(AliTriggerAnalysis::kMB1) 
  , fAnalysisMode(kTPCAnalysisMode) 
  , fTreeSRedirector(0)
  , fCentralityEstimator(0)
  , fLowPtTrackDownscaligF(0)
  , fLowPtV0DownscaligF(0)
  , fFriendDownscaling(-3.)   
  , fProcessAll(kFALSE)
  , fProcessCosmics(kFALSE)
  , fProcessITSTPCmatchOut(kFALSE)  // swittch to process ITS/TPC standalone tracks
  , fHighPtTree(0)
  , fV0Tree(0)
  , fdEdxTree(0)
  , fLaserTree(0)
  , fMCEffTree(0)
  , fCosmicPairsTree(0)
  , fPtResPhiPtTPC(0)
  , fPtResPhiPtTPCc(0)
  , fPtResPhiPtTPCITS(0)
  , fPtResEtaPtTPC(0)
  , fPtResEtaPtTPCc(0)
  , fPtResEtaPtTPCITS(0)
  , fPtResCentPtTPC(0)
  , fPtResCentPtTPCc(0)
  , fPtResCentPtTPCITS(0)
  , fCurrentFileName("")
  , fDummyTrack(0)
{
  // Constructor

  // Define input and output slots here
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
  DefineOutput(5, TTree::Class());
  DefineOutput(6, TTree::Class());
  DefineOutput(7, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskFilteredTree::~AliAnalysisTaskFilteredTree()
{
  //
  // Destructor
  //
  delete fFilteredTreeEventCuts;
  delete fFilteredTreeAcceptanceCuts;
  delete fFilteredTreeRecAcceptanceCuts;
  delete fEsdTrackCuts;
}

//____________________________________________________________________________
Bool_t AliAnalysisTaskFilteredTree::Notify()
{
  //
  //
  
  static Int_t count = 0;
  count++;
  TTree *chain = (TChain*)GetInputData(0);
  if(chain){
    Printf("Processing %d. file: %s", count, chain->GetCurrentFile()->GetName());
  }  
  fCurrentFileName=chain->GetCurrentFile()->GetName();   // This is the way to get path infor processing on the Lego train
  //  
  if (fCurrentFileName.String().CountChar('/')<1){ 
    //  in case path to input is not absolute use other handle to guess chunk ID
    //  in case filtering done locally, use the Env varaible AliEn output path
    //    this will work only in case we are doing alien processing
    TString fns = gSystem->Getenv("OutputDir");
    if (!fns.IsNull() && !fns.EndsWith("/")) fns += "/";
    fns += chain->GetCurrentFile()->GetName();
    Printf("Processing %d. file: %s", count, fns.Data());
    fCurrentFileName = fns.Data();
  }

  return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskFilteredTree::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  //
  //get the output file to make sure the trees will be associated to it
  OpenFile(1);
  fTreeSRedirector = new TTreeSRedirector();

  //
  // Create trees
  fV0Tree = ((*fTreeSRedirector)<<"V0s").GetTree();
  fHighPtTree = ((*fTreeSRedirector)<<"highPt").GetTree();
  fdEdxTree = ((*fTreeSRedirector)<<"dEdx").GetTree();
  fLaserTree = ((*fTreeSRedirector)<<"Laser").GetTree();
  fMCEffTree = ((*fTreeSRedirector)<<"MCEffTree").GetTree();
  fCosmicPairsTree = ((*fTreeSRedirector)<<"CosmicPairs").GetTree();

  if (!fDummyTrack)  {
    fDummyTrack=new AliESDtrack();
  }

  // histogram booking

  Double_t minPt = 0.1; 
  Double_t maxPt = 100.; 
  Int_t nbinsPt = 30; 

  Double_t logminPt = TMath::Log10(minPt);
  Double_t logmaxPt = TMath::Log10(maxPt);
  Double_t binwidth = (logmaxPt-logminPt)/nbinsPt;
  Double_t *binsPt =  new Double_t[nbinsPt+1];
  binsPt[0] = minPt;
  for (Int_t i=1;i<=nbinsPt;i++) {
    binsPt[i] = minPt + TMath::Power(10,logminPt+i*binwidth);
  }

  // 1pT resol cov matrix bins
  Double_t min1PtRes = 0.; 
  Double_t max1PtRes = 0.3; 
  Int_t nbins1PtRes = 300; 
  Double_t bins1PtRes[301];
  for (Int_t i=0;i<=nbins1PtRes;i++) {
    bins1PtRes[i] = min1PtRes + i*(max1PtRes-min1PtRes)/nbins1PtRes;
  }

  // phi bins
  Double_t minPhi = 0.; 
  Double_t maxPhi = 6.5; 
  Int_t nbinsPhi = 100; 
  Double_t binsPhi[101];
  for (Int_t i=0;i<=nbinsPhi;i++) {
    binsPhi[i] = minPhi + i*(maxPhi-minPhi)/nbinsPhi;
  }

  // eta bins
  Double_t minEta = -1.;
  Double_t maxEta = 1.;
  Int_t nbinsEta = 20;
  Double_t binsEta[21];
  for (Int_t i=0;i<=nbinsEta;i++) {
    binsEta[i] = minEta + i*(maxEta-minEta)/nbinsEta;
  }

  // mult bins
  Double_t minCent = 0.;
  Double_t maxCent = 100;
  Int_t nbinsCent = 20;
  Double_t binsCent[101];
  for (Int_t i=0;i<=nbinsCent;i++) {
    binsCent[i] = minCent + i*(maxCent-minCent)/nbinsCent;
  }

  fPtResPhiPtTPC = new TH3D("fPtResPhiPtTPC","pt rel. resolution from cov. matrix TPC tracks",nbinsPt,binsPt,nbinsPhi,binsPhi,nbins1PtRes,bins1PtRes);
  fPtResPhiPtTPCc = new TH3D("fPtResPhiPtTPCc","pt rel. resolution from cov. matrix TPC constrained tracks",nbinsPt,binsPt,nbinsPhi,binsPhi,nbins1PtRes,bins1PtRes);
  fPtResPhiPtTPCITS = new TH3D("fPtResPhiPtTPCITS","pt rel. resolution from cov. matrix TPC+ITS tracks",nbinsPt,binsPt,nbinsPhi,binsPhi,nbins1PtRes,bins1PtRes);

  fPtResEtaPtTPC = new TH3D("fPtResEtaPtTPC","pt rel. resolution from cov. matrix TPC tracks",nbinsPt,binsPt,nbinsEta,binsEta,nbins1PtRes,bins1PtRes);
  fPtResEtaPtTPCc = new TH3D("fPtResEtaPtTPCc","pt rel. resolution from cov. matrix TPC constrained tracks",nbinsPt,binsPt,nbinsEta,binsEta,nbins1PtRes,bins1PtRes);
  fPtResEtaPtTPCITS = new TH3D("fPtResEtaPtTPCITS","pt rel. resolution from cov. matrix TPC+ITS tracks",nbinsPt,binsPt,nbinsEta,binsEta,nbins1PtRes,bins1PtRes);

  fPtResCentPtTPC = new TH3D("fPtResCentPtTPC","pt rel. resolution from cov. matrix TPC tracks",nbinsPt,binsPt,nbinsCent,binsCent,nbins1PtRes,bins1PtRes);
  fPtResCentPtTPCc = new TH3D("fPtResCentPtTPCc","pt rel. resolution from cov. matrix TPC constrained tracks",nbinsPt,binsPt,nbinsCent,binsCent,nbins1PtRes,bins1PtRes);
  fPtResCentPtTPCITS = new TH3D("fPtResCentPtTPCITS","pt rel. resolution from cov. matrix TPC+ITS tracks",nbinsPt,binsPt,nbinsCent,binsCent,nbins1PtRes,bins1PtRes);


  fOutput = new TList; 
  if(!fOutput) return;
  fOutput->SetOwner();

  fOutput->Add(fPtResPhiPtTPC);
  fOutput->Add(fPtResPhiPtTPCc);
  fOutput->Add(fPtResPhiPtTPCITS);
  fOutput->Add(fPtResEtaPtTPC);
  fOutput->Add(fPtResEtaPtTPCc);
  fOutput->Add(fPtResEtaPtTPCITS);
  fOutput->Add(fPtResCentPtTPC);
  fOutput->Add(fPtResCentPtTPCc);
  fOutput->Add(fPtResCentPtTPCITS);

  // post data to outputs

  PostData(1,fV0Tree);
  PostData(2,fHighPtTree);
  PostData(3,fdEdxTree);
  PostData(4,fLaserTree);
  PostData(5,fMCEffTree);
  PostData(6,fCosmicPairsTree);

  PostData(7,fOutput);
}

//_____________________________________________________________________________
void AliAnalysisTaskFilteredTree::UserExec(Option_t *) 
{
  //
  // Called for each event
  //

  // ESD event
  fESD = (AliESDEvent*) (InputEvent());
  if (!fESD) {
    Printf("ERROR: ESD event not available");
    return;
  }
  //if MC info available - use it.
  fMC = MCEvent();
  if (fMC){  
    // Bug fix 28.05.2016 - do not trust to presence of MC handler, check if the content is valid
    //                    - proper solution (autodetection of MC information) to be implemented 
    if (fMC->Stack()==NULL) {
      fMC=NULL; 
    }else{
      AliInfo("ToFix: MC stack not available. Prefered MCEvent() will return 0");
    }
  }

  if(fUseESDfriends) {
    //fESDfriend = dynamic_cast<AliESDfriend*>(fESD->FindListObject("AliESDfriend"));
    fESDfriend = ESDfriend();
    if(!fESDfriend) {
      Printf("ERROR: ESD friends not available");
    }
  }
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!inputHandler){
    return;
  }

  //if set, use the environment variables to set the downscaling factors
  //AliAnalysisTaskFilteredTree_fLowPtTrackDownscaligF
  //AliAnalysisTaskFilteredTree_fLowPtV0DownscaligF
  //AliAnalysisTaskFilteredTree_fFriendDownscaling
  TString env;
  env = gSystem->Getenv("AliAnalysisTaskFilteredTree_fLowPtTrackDownscaligF");
  if (!env.IsNull()){
    fLowPtTrackDownscaligF=env.Atof();
    AliInfo(Form("fLowPtTrackDownscaligF=%f",fLowPtTrackDownscaligF));
  }
  env = gSystem->Getenv("AliAnalysisTaskFilteredTree_fLowPtV0DownscaligF");
  if (!env.IsNull()){
    fLowPtV0DownscaligF=env.Atof();
    AliInfo(Form("fLowPtV0DownscaligF=%f",fLowPtTrackDownscaligF));
  }
  env = gSystem->Getenv("AliAnalysisTaskFilteredTree_fFriendDownscaling");
  if (!env.IsNull()){
    fFriendDownscaling=env.Atof();
    AliInfo(Form(" fFriendDownscaling=%f",fFriendDownscaling));
  }
  //
  //
  //
  if(fProcessAll) { 
    ProcessAll(fESD,fMC,fESDfriend); // all track stages and MC
  }
  else {
    Process(fESD,fMC,fESDfriend);    // only global and TPC tracks
  }
  //
  ProcessV0(fESD,fMC,fESDfriend);
  ProcessLaser(fESD,fMC,fESDfriend);
  ProcessdEdx(fESD,fMC,fESDfriend);
  if (fProcessCosmics) { ProcessCosmics(fESD,fESDfriend); }
  if(fMC) { ProcessMCEff(fESD,fMC,fESDfriend);}
  if (fProcessITSTPCmatchOut) ProcessITSTPCmatchOut(fESD, fESDfriend);
  printf("processed event %d\n", Int_t(Entry()));
}

//_____________________________________________________________________________
void AliAnalysisTaskFilteredTree::ProcessCosmics(AliESDEvent *const event, AliESDfriend* esdFriend)
{
  //
  // Find cosmic pairs (triggered or random) 
  //
  //
  AliESDVertex *vertexSPD =  (AliESDVertex *)event->GetPrimaryVertexSPD();
  AliESDVertex *vertexTPC =  (AliESDVertex *)event->GetPrimaryVertexTPC(); 
  const Double_t kMinPt=0.8;
  const Double_t kMinPtMax=0.8;
  const Double_t kMinNcl=50;
  const Double_t kMaxDelta[5]={2,600,0.02,0.02,0.1};
  Int_t ntracks=event->GetNumberOfTracks(); 
  UInt_t specie = event->GetEventSpecie();  // skip laser events
  if (specie==AliRecoParam::kCalib) return;
  Int_t ntracksFriend = esdFriend ? esdFriend->GetNumberOfTracks() : 0;


  for (Int_t itrack0=0;itrack0<ntracks;itrack0++) {
    AliESDtrack *track0 = event->GetTrack(itrack0);
    if (!track0) continue;
    if (!track0->IsOn(AliESDtrack::kTPCrefit)) continue;

    if (TMath::Abs(AliTracker::GetBz())>1 && track0->Pt() < kMinPt) continue;
    if (track0->Pt() < kMinPt) continue;
    if (track0->GetTPCncls() < kMinNcl) continue;
    if (TMath::Abs(track0->GetY())<kMaxDelta[0]) continue; 
    if (track0->GetKinkIndex(0)>0) continue;
    const Double_t * par0=track0->GetParameter(); //track param at rhe DCA
    //rm primaries
    //
    //track0->GetImpactParametersTPC(dcaTPC,covTPC);
    //if (TMath::Abs(dcaTPC[0])<kMaxDelta[0]) continue;
    //if (TMath::Abs(dcaTPC[1])<kMaxDelta[0]*2) continue;
    //    const AliExternalTrackParam * trackIn0 = track0->GetInnerParam();
    AliESDfriendTrack* friendTrack0=NULL;
    if (esdFriend &&!esdFriend->TestSkipBit()){
      if (itrack0<ntracksFriend){
	friendTrack0 = esdFriend->GetTrack(itrack0);
      } //this guy can be NULL
    }

    for (Int_t itrack1=itrack0+1;itrack1<ntracks;itrack1++) {
      AliESDtrack *track1 = event->GetTrack(itrack1);
      if (!track1) continue;  
      if (!track1->IsOn(AliESDtrack::kTPCrefit)) continue;
      if (track1->GetKinkIndex(0)>0) continue;
      if ((TMath::Abs(AliTracker::GetBz())>1) && (track1->Pt() < kMinPt)) continue;
      if (track1->Pt() < kMinPt) continue;
      if (track1->GetTPCncls()<kMinNcl) continue;
      if (TMath::Abs(AliTracker::GetBz())>1 && TMath::Max(track1->Pt(), track0->Pt())<kMinPtMax) continue;
      if (TMath::Abs(track1->GetY())<kMaxDelta[0]) continue;
      //track1->GetImpactParametersTPC(dcaTPC,covTPC);
      //      if (TMath::Abs(dcaTPC[0])<kMaxDelta[0]) continue;
      //if (TMath::Abs(dcaTPC[1])<kMaxDelta[0]*2) continue;
      //
      const Double_t* par1=track1->GetParameter(); //track param at rhe DCA
      //
      Bool_t isPair=kTRUE;
      for (Int_t ipar=0; ipar<5; ipar++){
        if (ipar==4&&TMath::Abs(AliTracker::GetBz())<1) continue; // 1/pt not defined for B field off
        if (TMath::Abs(TMath::Abs(par0[ipar])-TMath::Abs(par1[ipar]))>kMaxDelta[ipar]) isPair=kFALSE;
      }
      if (!isPair) continue;
      if (TMath::Abs(TMath::Abs(track0->GetAlpha()-track1->GetAlpha())-TMath::Pi())>kMaxDelta[2]) isPair=kFALSE;
      //delta with correct sign
      /*
	TCut cut0="abs(t1.fP[0]+t0.fP[0])<2"
	TCut cut3="abs(t1.fP[3]+t0.fP[3])<0.02"
	TCut cut4="abs(t1.fP[4]+t0.fP[4])<0.2"
      */
      if  (TMath::Abs(par0[0]+par1[0])>kMaxDelta[0]) isPair=kFALSE; //delta y   opposite sign
      if  (TMath::Abs(par0[3]+par1[3])>kMaxDelta[3]) isPair=kFALSE; //delta tgl opposite sign
      if  (TMath::Abs(AliTracker::GetBz())>1 && TMath::Abs(par0[4]+par1[4])>kMaxDelta[4]) isPair=kFALSE; //delta 1/pt opposite sign
      if (!isPair) continue;
      TString filename(AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile()->GetName());
      Int_t eventNumber = event->GetEventNumberInFile(); 
      //
      //               
      Int_t ntracksSPD = vertexSPD->GetNContributors();
      Int_t ntracksTPC = vertexTPC->GetNContributors();        
      Int_t runNumber     = event->GetRunNumber();        
      Int_t timeStamp    = event->GetTimeStamp();
      ULong64_t triggerMask = event->GetTriggerMask();
      Float_t magField    = event->GetMagneticField();
      TObjString triggerClass = event->GetFiredTriggerClasses().Data();

      // Global event id calculation using orbitID, bunchCrossingID and periodID
      ULong64_t orbitID      = (ULong64_t)event->GetOrbitNumber();
      ULong64_t bunchCrossID = (ULong64_t)event->GetBunchCrossNumber();
      ULong64_t periodID     = (ULong64_t)event->GetPeriodNumber();
      ULong64_t gid          = ((periodID << 36) | (orbitID << 12) | bunchCrossID); 
      

      AliESDfriendTrack* friendTrack1=NULL;
      if (esdFriend &&!esdFriend->TestSkipBit()){
	if (itrack1<ntracksFriend){
	  friendTrack1 = esdFriend->GetTrack(itrack1);
	} //this guy can be NULL
      }

      //
      AliESDfriendTrack *friendTrackStore0=friendTrack0;    // store friend track0 for later processing
      AliESDfriendTrack *friendTrackStore1=friendTrack1;    // store friend track1 for later processing
      if (fFriendDownscaling>=1){  // downscaling number of friend tracks
	if (gRandom->Rndm()>1./fFriendDownscaling){
	  friendTrackStore0 = 0;
	  friendTrackStore1 = 0;
	}
      }
      if (fFriendDownscaling<=0){
	if (((*fTreeSRedirector)<<"CosmicPairs").GetTree()){
	  TTree * tree = ((*fTreeSRedirector)<<"CosmicPairs").GetTree();
	  if (tree){
	    Double_t sizeAll=tree->GetZipBytes();
	    TBranch * br= tree->GetBranch("friendTrack0.fPoints");
	    Double_t sizeFriend=(br!=NULL)?br->GetZipBytes():0;
	    br= tree->GetBranch("friendTrack0.fCalibContainer");
	    if (br) sizeFriend+=br->GetZipBytes();
	    if (sizeFriend*TMath::Abs(fFriendDownscaling)>sizeAll) {
	      friendTrackStore0=0;
	      friendTrackStore1=0;
	    }
	  }
	}
      }
      if(!fFillTree) return;
      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"CosmicPairs"<<
        "gid="<<gid<<                         // global id of track
        "fileName.="<<&fCurrentFileName<<     // file name
        "runNumber="<<runNumber<<             // run number	    
        "evtTimeStamp="<<timeStamp<<          // time stamp of event
        "evtNumberInFile="<<eventNumber<<     // event number	    
        "trigger="<<triggerMask<<             // trigger mask
        "triggerClass="<<&triggerClass<<      // trigger class
        "Bz="<<magField<<                     // magnetic field
        //
        "multSPD="<<ntracksSPD<<              // event ultiplicity
        "multTPC="<<ntracksTPC<<              //  
        "vertSPD.="<<vertexSPD<<              // primary vertex -SPD
        "vertTPC.="<<vertexTPC<<              // primary vertex -TPC
        "t0.="<<track0<<                      // first half of comsic trak
        "t1.="<<track1<<                      // second half of cosmic track
        "friendTrack0.="<<friendTrackStore0<< // friend information first track  + points
        "friendTrack1.="<<friendTrackStore1<< // frined information first track  + points 
        "\n";      
    }
  }
}


//_____________________________________________________________________________
void AliAnalysisTaskFilteredTree::Process(AliESDEvent *const esdEvent, AliMCEvent * const mcEvent, AliESDfriend *const /*esdFriend*/)
{
  //
  // Select real events with high-pT tracks 
  //
  static Int_t downscaleCounter=0;
  // get selection cuts
  AliFilteredTreeEventCuts *evtCuts = GetEventCuts(); 
  AliFilteredTreeAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 

  if(!evtCuts || !accCuts  || !esdTrackCuts) {
    Printf("ERROR cuts not available");
    return;
  }

  // trigger selection
  Bool_t isEventTriggered = kTRUE;
  AliPhysicsSelection *physicsSelection = NULL;
  AliTriggerAnalysis* triggerAnalysis = NULL;

  // 
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  // trigger
  if(evtCuts->IsTriggerRequired())  
  {
    // always MB
    isEventTriggered = inputHandler->IsEventSelected() & AliVEvent::kMB;

    physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if(!physicsSelection) return;
    //SetPhysicsTriggerSelection(physicsSelection);

    if (isEventTriggered && (GetTrigger() == AliTriggerAnalysis::kV0AND)) {
      // set trigger (V0AND)
      triggerAnalysis = physicsSelection->GetTriggerAnalysis();
      if(!triggerAnalysis) return;
      isEventTriggered = triggerAnalysis->IsOfflineTriggerFired(esdEvent, GetTrigger());
    }
  }

  // centrality determination
  Float_t centralityF = -1;
  AliCentrality *esdCentrality = esdEvent->GetCentrality();
  centralityF = esdCentrality->GetCentralityPercentile(fCentralityEstimator.Data());

  // use MC information
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  TArrayF vtxMC(3);

  Int_t multMCTrueTracks = 0;
  if(mcEvent)
  {
    // get MC event header
    header = mcEvent->Header();
    if (!header) {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }
    // MC particle stack
    stack = mcEvent->Stack();
    if (!stack) {
      AliDebug(AliLog::kError, "Stack not available");
      return;
    }

    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);

    // multipliticy of all MC primary tracks
    // in Zv, pt and eta ranges)
    multMCTrueTracks = GetMCTrueTrackMult(mcEvent,evtCuts,accCuts);
  } // end bUseMC

  // get reconstructed vertex  
  //const AliESDVertex* vtxESD = 0; 
  AliESDVertex* vtxESD = 0; 
  if(GetAnalysisMode() == kTPCAnalysisMode) {
    vtxESD = (AliESDVertex*)esdEvent->GetPrimaryVertexTPC();
  }
  else if(GetAnalysisMode() == kTPCITSAnalysisMode) {
    vtxESD = (AliESDVertex*)esdEvent->GetPrimaryVertexTracks();
  }
  else {
    return;
  }

  AliESDVertex* vtxTPC = (AliESDVertex*)esdEvent->GetPrimaryVertexTPC();
  AliESDVertex* vtxSPD = (AliESDVertex*)esdEvent->GetPrimaryVertexSPD();

  if(!vtxESD) return;
  if(!vtxTPC) return;
  if(!vtxSPD) return;

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD); 
  //printf("isEventOK %d, isEventTriggered %d, status %d, vz %f \n",isEventOK, isEventTriggered, vtxESD->GetStatus(), vtxESD->GetZ());
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());
  Int_t ntracks = esdEvent->GetNumberOfTracks();

  // check event cuts
  if(isEventOK && isEventTriggered)
  {

    //
    // get IR information
    //
    AliESDHeader *esdHeader = 0;
    esdHeader = esdEvent->GetHeader();
    if(!esdHeader) return;
    //Int_t ir1 = esdHeader->GetTriggerIREntries(); //all ir-s
    //Int_t ir2 = esdHeader->GetTriggerIREntries(-1,1); // int2 set, 180 ms time interval

    // Use when Peter commit the changes in the header
    Int_t ir1 = 0;
    Int_t ir2 = 0;

    //
    //Double_t vert[3] = {0}; 
    //vert[0] = vtxESD->GetX();
    //vert[1] = vtxESD->GetY();
    //vert[2] = vtxESD->GetZ();
    Int_t mult = vtxESD->GetNContributors();
    Int_t multSPD = vtxSPD->GetNContributors();
    Int_t multTPC = vtxTPC->GetNContributors();

    Float_t bz = esdEvent->GetMagneticField();
    Int_t runNumber = esdEvent->GetRunNumber();
    Int_t evtTimeStamp = esdEvent->GetTimeStamp();
    Int_t evtNumberInFile = esdEvent->GetEventNumberInFile();
   
    // Global event id calculation using orbitID, bunchCrossingID and periodID 
    ULong64_t orbitID      = (ULong64_t)esdEvent->GetOrbitNumber();
    ULong64_t bunchCrossID = (ULong64_t)esdEvent->GetBunchCrossNumber();
    ULong64_t periodID     = (ULong64_t)esdEvent->GetPeriodNumber();
    ULong64_t gid          = ((periodID << 36) | (orbitID << 12) | bunchCrossID); 
    

    // high pT tracks
    for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++)
    {
      AliESDtrack *track = esdEvent->GetTrack(iTrack);
      if(!track) continue;
      if(track->Charge()==0) continue;
      if(!esdTrackCuts->AcceptTrack(track)) continue;
      if(!accCuts->AcceptTrack(track)) continue;

      // downscale low-pT tracks
      Double_t scalempt= TMath::Min(track->Pt(),10.);
      Double_t downscaleF = gRandom->Rndm();
      downscaleF *= fLowPtTrackDownscaligF;
      if( downscaleCounter>0 && TMath::Exp(2*scalempt)<downscaleF) continue;
      //printf("TMath::Exp(2*scalempt) %e, downscaleF %e \n",TMath::Exp(2*scalempt), downscaleF);

      AliExternalTrackParam * tpcInner = (AliExternalTrackParam *)(track->GetTPCInnerParam());
      if (!tpcInner) continue;
      // transform to the track reference frame 
      Bool_t isOK = kFALSE;
      isOK = tpcInner->Rotate(track->GetAlpha());
      isOK = tpcInner->PropagateTo(track->GetX(),esdEvent->GetMagneticField());
      if(!isOK) continue;

      // Dump to the tree 
      // vertex
      // TPC-ITS tracks
      //
      TObjString triggerClass = esdEvent->GetFiredTriggerClasses().Data();
      if(!fFillTree) return;
      if(!fTreeSRedirector) return;
      downscaleCounter++;
      (*fTreeSRedirector)<<"highPt"<<
        "gid="<<gid<<
        "fileName.="<<&fCurrentFileName<<            
        "runNumber="<<runNumber<<
        "evtTimeStamp="<<evtTimeStamp<<
        "evtNumberInFile="<<evtNumberInFile<<
        "triggerClass="<<&triggerClass<<      //  trigger
        "Bz="<<bz<<                           //  magnetic field
        "vtxESD.="<<vtxESD<<
        "ntracksESD="<<ntracks<<              // number of tracks in the ESD
        "IRtot="<<ir1<<                      // interaction record history info
        "IRint2="<<ir2<<
        "mult="<<mult<<                      // multiplicity of tracks pointing to the primary vertex
        "multSPD="<<multSPD<<                // multiplicity of tracks pointing to the SPD primary vertex
        "multTPC="<<multTPC<<                // multiplicity of tracks pointing to the TPC primary vertex
        "esdTrack.="<<track<<
        "centralityF="<<centralityF<<       
        "\n";
    }
  }

}


//_____________________________________________________________________________
void AliAnalysisTaskFilteredTree::ProcessLaser(AliESDEvent *const esdEvent, AliMCEvent * const /*mcEvent*/, AliESDfriend *const esdFriend)
{
  //
  // Process laser events -> dump tracks and clusters  to the special tree
  //
  const Double_t kMinPt = 5; 
  if(!fFillTree) return;
  if(!fTreeSRedirector) return;
  const AliESDHeader* esdHeader = esdEvent->GetHeader();
  if(esdHeader && esdHeader->GetEventSpecie()==AliRecoParam::kCalib) {
    Int_t countLaserTracks = 0;
    Int_t runNumber = esdEvent->GetRunNumber();
    Int_t evtTimeStamp = esdEvent->GetTimeStamp();
    Int_t evtNumberInFile = esdEvent->GetEventNumberInFile();
    Float_t bz = esdEvent->GetMagneticField();
    TObjString triggerClass = esdEvent->GetFiredTriggerClasses().Data();    
    // Global event id calculation using orbitID, bunchCrossingID and periodID
    ULong64_t orbitID      = (ULong64_t)esdEvent->GetOrbitNumber();
    ULong64_t bunchCrossID = (ULong64_t)esdEvent->GetBunchCrossNumber();
    ULong64_t periodID     = (ULong64_t)esdEvent->GetPeriodNumber();
    ULong64_t gid          = ((periodID << 36) | (orbitID << 12) | bunchCrossID); 
    for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++){
      AliESDtrack *track = esdEvent->GetTrack(iTrack);
      if(!track) continue;
      if (track->GetTPCInnerParam()==NULL) continue;
      if(track->GetTPCInnerParam()) countLaserTracks++;      
      AliESDfriendTrack* friendTrack=NULL;
      // suppress beam background and CE random reacks
      if (track->GetInnerParam()->Pt()<kMinPt) continue;
      Bool_t skipTrack=gRandom->Rndm()>1/(1+TMath::Abs(fFriendDownscaling));
      if (skipTrack) continue;
      if (esdFriend) {if (!esdFriend->TestSkipBit()) friendTrack = esdFriend->GetTrack(iTrack);} //this guy can be NULL      
      (*fTreeSRedirector)<<"Laser"<<
        "gid="<<gid<<                          // global identifier of event
        "fileName.="<<&fCurrentFileName<<              //
        "runNumber="<<runNumber<<
        "evtTimeStamp="<<evtTimeStamp<<
        "evtNumberInFile="<<evtNumberInFile<<
        "triggerClass="<<&triggerClass<<        //  trigger
        "Bz="<<bz<<                             //  magnetic field
        "multTPCtracks="<<countLaserTracks<<    //  multiplicity of tracks
	"track.="<<track<<                      //  track parameters
        "friendTrack.="<<friendTrack<<          //  friend track information
        "\n";
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskFilteredTree::ProcessAll(AliESDEvent *const esdEvent, AliMCEvent * const mcEvent, AliESDfriend *const esdFriend)
{
  //
  // Select real events with high-pT tracks
  // Calculate and stor in the output tree:
  //  TPC constrained tracks
  //  InnerParams constrained tracks
  //  TPC-ITS tracks
  //  ITSout-InnerParams tracks
  //  chi2 distance between TPC constrained and TPC-ITS tracks
  //  chi2 distance between TPC refitted constrained and TPC-ITS tracks
  //  chi2 distance between ITSout and InnerParams tracks
  //  MC information: 
  //   track references at ITSin, TPCin; InnerParam at first TPC track reference, 
  //   particle ID, mother ID, production mechanism ...
  // 
  // get selection cuts
  static Int_t downscaleCounter=0;
  AliFilteredTreeEventCuts *evtCuts = GetEventCuts(); 
  AliFilteredTreeAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 

  if(!evtCuts || !accCuts  || !esdTrackCuts) {
    AliDebug(AliLog::kError, "cuts not available");
    return;
  }

  // trigger selection
  Bool_t isEventTriggered = kTRUE;
  AliPhysicsSelection *physicsSelection = NULL;
  AliTriggerAnalysis* triggerAnalysis = NULL;

  // 
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  AliPIDResponse *pidResponse = inputHandler->GetPIDResponse();


  // trigger
  if(evtCuts->IsTriggerRequired())  
  {
    // always MB
    isEventTriggered = inputHandler->IsEventSelected() & AliVEvent::kMB;

    physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if(!physicsSelection) {AliInfo("no physics selection"); return;}
    //SetPhysicsTriggerSelection(physicsSelection);

    if (isEventTriggered && (GetTrigger() == AliTriggerAnalysis::kV0AND)) {
      // set trigger (V0AND)
      triggerAnalysis = physicsSelection->GetTriggerAnalysis();
      if(!triggerAnalysis) {AliInfo("no trigger analysis");return;}
      isEventTriggered = triggerAnalysis->IsOfflineTriggerFired(esdEvent, GetTrigger());
    }
  }

  // centrality determination
  Float_t centralityF = -1;
  AliCentrality *esdCentrality = esdEvent->GetCentrality();
  centralityF = esdCentrality->GetCentralityPercentile(fCentralityEstimator.Data());

  // use MC information
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  TArrayF vtxMC(3);
  Int_t mcStackSize=0;

  Int_t multMCTrueTracks = 0;
  if(mcEvent)
  {
    // get MC event header
    header = mcEvent->Header();
    if (!header) {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }
    // MC particle stack
    stack = mcEvent->Stack();
    if (!stack) {
      AliDebug(AliLog::kError, "Stack not available");
      return;
    }
    mcStackSize=stack->GetNtrack();

    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);

    // multipliticy of all MC primary tracks
    // in Zv, pt and eta ranges)
    multMCTrueTracks = GetMCTrueTrackMult(mcEvent,evtCuts,accCuts);

  } // end bUseMC

  // get reconstructed vertex  
  //const AliESDVertex* vtxESD = 0; 
  AliESDVertex* vtxESD = 0; 
  if(GetAnalysisMode() == kTPCAnalysisMode) {
    vtxESD = (AliESDVertex*)esdEvent->GetPrimaryVertexTPC();
  }
  else if(GetAnalysisMode() == kTPCITSAnalysisMode) {
    vtxESD = (AliESDVertex*)esdEvent->GetPrimaryVertexTracks();
  }
  else {
    AliInfo("no ESD vertex");
    return;
  }

  if(!vtxESD) return;

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD); 

  //
  // Vertex info comparison and track multiplicity
  //
  AliESDVertex *vertexSPD =  (AliESDVertex *)esdEvent->GetPrimaryVertexSPD();
  AliESDVertex *vertexTPC =  (AliESDVertex *)esdEvent->GetPrimaryVertexTPC(); 
  Int_t contSPD = vertexSPD->GetNContributors();
  Int_t contTPC = vertexTPC->GetNContributors();        
  TVectorD vertexPosTPC(3), vertexPosSPD(3);
  vertexSPD->GetXYZ(vertexPosSPD.GetMatrixArray());
  vertexTPC->GetXYZ(vertexPosTPC.GetMatrixArray());
  Int_t ntracksTPC=0;
  Int_t ntracksITS=0;
  for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++){
    AliESDtrack *track = esdEvent->GetTrack(iTrack);    
    if(!track) continue;
    if (track->IsOn(AliVTrack::kTPCrefit)) ntracksTPC++;
    if (track->IsOn(AliVTrack::kITSrefit)) ntracksITS++;
  }

  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK, isEventTriggered);
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());
  Int_t ntracks = esdEvent->GetNumberOfTracks();
  // Global event id calculation using orbitID, bunchCrossingID and periodID 
  ULong64_t orbitID      = (ULong64_t)esdEvent->GetOrbitNumber();
  ULong64_t bunchCrossID = (ULong64_t)esdEvent->GetBunchCrossNumber();
  ULong64_t periodID     = (ULong64_t)esdEvent->GetPeriodNumber();
  ULong64_t gid          = ((periodID << 36) | (orbitID << 12) | bunchCrossID); 
  TObjString triggerClass = esdEvent->GetFiredTriggerClasses().Data();
  Float_t bz = esdEvent->GetMagneticField();
  Int_t runNumber = esdEvent->GetRunNumber();
  Int_t evtTimeStamp = esdEvent->GetTimeStamp();
  Int_t evtNumberInFile = esdEvent->GetEventNumberInFile();
  Int_t mult = vtxESD->GetNContributors();
  (*fTreeSRedirector)<<"eventInfoTracks"<<
    "gid="<<gid<<
    "fileName.="<<&fCurrentFileName<<                // name of the chunk file (hopefully full)
    "runNumber="<<runNumber<<                             // runNumber
    "evtTimeStamp="<<evtTimeStamp<<           // time stamp of event (in seconds)
    "evtNumberInFile="<<evtNumberInFile<<     // event number
    "triggerClass="<<&triggerClass<<          // trigger class as a string
    "Bz="<<bz<<                               // solenoid magnetic field in the z direction (in kGaus)
    "mult="<<mult<<                           // multiplicity of tracks pointing to the primary vertex
    "ntracks="<<ntracks<<                     // number of the esd tracks (to take into account the pileup in the TPC)
    "isEventOK="<<isEventOK<<                 // flag - AliFilteredTreeEventCuts - track dumped only for selected events 
    "isEventTriggered="<<isEventTriggered<<   // flag - if tigger required - track dumped only for selected events 
    "\n";



  // check event cuts
  if(isEventOK && isEventTriggered)
  {
    //
    // get IR information
    //
    AliESDHeader *esdHeader = 0;
    esdHeader = esdEvent->GetHeader();
    if(!esdHeader) {AliInfo("no esdHeader");return;}
    //Int_t ir1 = esdHeader->GetTriggerIREntries(); //all ir-s
    //Int_t ir2 = esdHeader->GetTriggerIREntries(-1,1); // int2 set, 180 ms time interval
    //
    Int_t ir1 = 0;
    Int_t ir2 = 0;

    //
    Double_t vert[3] = {0}; 
    vert[0] = vtxESD->GetX();
    vert[1] = vtxESD->GetY();
    vert[2] = vtxESD->GetZ();
    Int_t mult = vtxESD->GetNContributors();
    Int_t numberOfTracks=esdEvent->GetNumberOfTracks();
    // high pT tracks
    for (Int_t iTrack = 0; iTrack < numberOfTracks; iTrack++)
    {
      AliESDtrack *track = esdEvent->GetTrack(iTrack);
      AliESDfriendTrack* friendTrack=NULL;
      Int_t numberOfFriendTracks=0;
      if (esdFriend) numberOfFriendTracks=esdFriend->GetNumberOfTracks();
      if (esdFriend && iTrack<numberOfFriendTracks) {if (!esdFriend->TestSkipBit()) friendTrack = esdFriend->GetTrack(iTrack);} //this guy can be NULL
      if(!track) continue;
      if(track->Charge()==0) continue;
      if(!esdTrackCuts->AcceptTrack(track)) continue;
      if(!accCuts->AcceptTrack(track)) continue;

      // downscale low-pT tracks
      Double_t scalempt= TMath::Min(track->Pt(),10.);
      Double_t downscaleF = gRandom->Rndm();
      downscaleF *= fLowPtTrackDownscaligF;
      if( downscaleCounter>0 && TMath::Exp(2*scalempt)<downscaleF) continue;
      //printf("TMath::Exp(2*scalempt) %e, downscaleF %e \n",TMath::Exp(2*scalempt), downscaleF);

      // Dump to the tree 
      // vertex
      // TPC constrained tracks
      // InnerParams constrained tracks
      // TPC-ITS tracks
      // ITSout-InnerParams tracks
      // chi2 distance between TPC constrained and TPC-ITS tracks
      // chi2 distance between TPC refitted constrained and TPC-ITS tracks
      // chi2 distance between ITSout and InnerParams tracks
      // MC information

      Double_t x[3]; track->GetXYZ(x);
      Double_t b[3]; AliTracker::GetBxByBz(x,b);

      //
      // Transform TPC inner params to track reference frame
      //
      Bool_t isOKtpcInner = kFALSE;
      AliExternalTrackParam * tpcInner = (AliExternalTrackParam *)(track->GetTPCInnerParam());
      if (tpcInner) {
        // transform to the track reference frame 
        isOKtpcInner = tpcInner->Rotate(track->GetAlpha());
        isOKtpcInner = tpcInner->PropagateTo(track->GetX(),esdEvent->GetMagneticField());
      }

      //
      // Constrain TPC inner to vertex
      // clone TPCinner has to be deleted
      //
      Bool_t isOKtpcInnerC = kFALSE;
      AliExternalTrackParam * tpcInnerC = new AliExternalTrackParam(*(track->GetTPCInnerParam()));
      if (tpcInnerC) {
        isOKtpcInnerC = ConstrainTPCInner(tpcInnerC,vtxESD,b);
        isOKtpcInnerC = tpcInnerC->Rotate(track->GetAlpha());
        isOKtpcInnerC = tpcInnerC->PropagateTo(track->GetX(),esdEvent->GetMagneticField());
      }

      //
      // Constrain TPC refitted tracks at inner TPC wall (InnerParams) to vertex  
      // Clone track InnerParams has to be deleted
      //
      Bool_t isOKtrackInnerC = kTRUE;
      AliExternalTrackParam * trackInnerC = NULL; 
      AliExternalTrackParam * trackInnerV =  new AliExternalTrackParam(*(track->GetInnerParam()));
      isOKtrackInnerC=AliTracker::PropagateTrackToBxByBz(trackInnerV,3,track->GetMass(),3,kFALSE);
      isOKtrackInnerC&= trackInnerV->Rotate(track->GetAlpha());
      isOKtrackInnerC&= trackInnerV->PropagateTo(track->GetX(),esdEvent->GetMagneticField());

      if (isOKtrackInnerC) {
	trackInnerC =  new AliExternalTrackParam(*trackInnerV);
        isOKtrackInnerC&= ConstrainTrackInner(trackInnerC,vtxESD,track->GetMass(),b);
        isOKtrackInnerC&= trackInnerC->Rotate(track->GetAlpha());
        isOKtrackInnerC&= trackInnerC->PropagateTo(track->GetX(),esdEvent->GetMagneticField());
      }       
      //
      // calculate chi2 between vi and vj vectors
      // with covi and covj covariance matrices
      // chi2ij = (vi-vj)^(T)*(covi+covj)^(-1)*(vi-vj)
      //
      TMatrixD deltaT(5,1), deltaTtrackC(5,1);
      TMatrixD delta(1,5),  deltatrackC(1,5);
      TMatrixD covarM(5,5), covarMtrackC(5,5);
      TMatrixD chi2(1,1);
      TMatrixD chi2trackC(1,1);

      if(isOKtpcInnerC && isOKtrackInnerC) 
      {
        for (Int_t ipar=0; ipar<5; ipar++) {
          deltaT(ipar,0)=tpcInnerC->GetParameter()[ipar]-track->GetParameter()[ipar];
          delta(0,ipar)=tpcInnerC->GetParameter()[ipar]-track->GetParameter()[ipar];

          deltaTtrackC(ipar,0)=trackInnerC->GetParameter()[ipar]-track->GetParameter()[ipar];
          deltatrackC(0,ipar)=trackInnerC->GetParameter()[ipar]-track->GetParameter()[ipar];

          for (Int_t jpar=0; jpar<5; jpar++) {
            Int_t index=track->GetIndex(ipar,jpar);
            covarM(ipar,jpar)=track->GetCovariance()[index]+tpcInnerC->GetCovariance()[index];
            covarMtrackC(ipar,jpar)=track->GetCovariance()[index]+trackInnerC->GetCovariance()[index];
          }
        }

        // chi2 distance TPC constrained and TPC+ITS
        TMatrixD covarMInv = covarM.Invert();
        TMatrixD mat2 = covarMInv*deltaT;
        chi2 = delta*mat2; 
        //chi2.Print();

        // chi2 distance TPC refitted constrained and TPC+ITS
        TMatrixD covarMInvtrackC = covarMtrackC.Invert();
        TMatrixD mat2trackC = covarMInvtrackC*deltaTtrackC;
        chi2trackC = deltatrackC*mat2trackC; 
        //chi2trackC.Print();
      }
      //
      // Find nearest combined and ITS standaalone tracks
      AliExternalTrackParam paramITS;     // nearest ITS track  -   chi2 distance at vertex
      AliExternalTrackParam paramITSC;    // nearest ITS track  -   to constrained track   chi2 distance at vertex
      AliExternalTrackParam paramComb;    // nearest comb. tack -   chi2 distance at inner wall
      Int_t indexNearestITS   = GetNearestTrack((trackInnerV!=NULL)? trackInnerV:track, iTrack, esdEvent,0,0,paramITS); 
      if (indexNearestITS<0)  indexNearestITS   = GetNearestTrack((trackInnerV!=NULL)? trackInnerV:track, iTrack, esdEvent,2,0,paramITS);
      Int_t indexNearestITSC  = GetNearestTrack((trackInnerC!=NULL)? trackInnerC:track, iTrack, esdEvent,0,0,paramITSC);  
      if (indexNearestITSC<0)  indexNearestITS   = GetNearestTrack((trackInnerC!=NULL)? trackInnerC:track, iTrack, esdEvent,2,0,paramITSC);
      Int_t indexNearestComb  = GetNearestTrack(track->GetInnerParam(), iTrack, esdEvent,1,1,paramComb);  
      
      //
      // Propagate ITSout to TPC inner wall 
      // and calculate chi2 distance to track (InnerParams)
      //
      const Double_t kTPCRadius=85; 
      const Double_t kStep=3; 

      // clone track InnerParams has to be deleted
      Bool_t isOKtrackInnerC2 = kFALSE;
      AliExternalTrackParam *trackInnerC2 = new AliExternalTrackParam(*(track->GetInnerParam()));
      if (trackInnerC2) {
        isOKtrackInnerC2 = AliTracker::PropagateTrackToBxByBz(trackInnerC2,kTPCRadius,track->GetMass(),kStep,kFALSE);
      }

      Bool_t isOKouterITSc = kFALSE;
      AliExternalTrackParam *outerITSc = NULL;
      TMatrixD chi2OuterITS(1,1);

      if(esdFriend && !esdFriend->TestSkipBit()) 
      {
        // propagate ITSout to TPC inner wall
        if(friendTrack) 
        {

          outerITSc = NULL;
          if (friendTrack->GetITSOut()) outerITSc = new AliExternalTrackParam(*(friendTrack->GetITSOut()));
          if(outerITSc) 
          {
            isOKouterITSc = AliTracker::PropagateTrackToBxByBz(outerITSc,kTPCRadius,track->GetMass(),kStep,kFALSE);
            isOKouterITSc = outerITSc->Rotate(trackInnerC2->GetAlpha());
            isOKouterITSc = outerITSc->PropagateTo(trackInnerC2->GetX(),esdEvent->GetMagneticField());

            //
            // calculate chi2 between outerITS and innerParams
            // cov without z-coordinate at the moment
            //
            TMatrixD deltaTouterITS(4,1);
            TMatrixD deltaouterITS(1,4);
            TMatrixD covarMouterITS(4,4);

            if(isOKtrackInnerC2 && isOKouterITSc) {
              Int_t kipar = 0;
              Int_t kjpar = 0;
              for (Int_t ipar=0; ipar<5; ipar++) {
                if(ipar!=1) {
                  deltaTouterITS(kipar,0)=outerITSc->GetParameter()[ipar]-trackInnerC2->GetParameter()[ipar];
                  deltaouterITS(0,kipar)=outerITSc->GetParameter()[ipar]-trackInnerC2->GetParameter()[ipar];
                }

                kjpar=0;
                for (Int_t jpar=0; jpar<5; jpar++) {
                  Int_t index=outerITSc->GetIndex(ipar,jpar);
                  if(ipar !=1 || jpar!=1) {
                    covarMouterITS(kipar,kjpar)=outerITSc->GetCovariance()[index]+trackInnerC2->GetCovariance()[index];
                  }
                  if(jpar!=1)  kjpar++;
                }
                if(ipar!=1) kipar++;
              }

              // chi2 distance ITSout and InnerParams
              TMatrixD covarMInvouterITS = covarMouterITS.Invert();
              TMatrixD mat2outerITS = covarMInvouterITS*deltaTouterITS;
              chi2OuterITS = deltaouterITS*mat2outerITS; 
              //chi2OuterITS.Print();
            } 
          }
        }
      }

      //
      // MC info
      //
      TParticle *particle=NULL, *particleTPC=NULL, *particleITS=NULL;
      TParticle *particleMother=NULL, *particleMotherTPC=NULL, *particleMotherITS=NULL;
      Int_t mech=-1, mechTPC=-1, mechITS=-1;
      Bool_t isPrim=kFALSE, isPrimTPC=kFALSE, isPrimITS=kFALSE;
      Bool_t isFromStrangess=kFALSE, isFromStrangessTPC=kFALSE, isFromStrangessITS=kFALSE;
      Bool_t isFromConversion=kFALSE, isFromConversionTPC=kFALSE, isFromConversionITS=kFALSE;
      Bool_t isFromMaterial=kFALSE, isFromMaterialTPC=kFALSE, isFromMaterialITS=kFALSE;

      AliTrackReference *refTPCIn = NULL;
      AliTrackReference *refTPCOut = NULL;
      AliTrackReference *refITS = NULL;
      AliTrackReference *refTRD = NULL;
      AliTrackReference *refTOF = NULL;
      AliTrackReference *refEMCAL = NULL;
      AliTrackReference *refPHOS = NULL;
      Int_t nrefTPC=0, nrefTRD=0, nrefTOF=0, nrefITS=0, nrefEMCAL=0, nrefPHOS=0;

      Bool_t isOKtrackInnerC3 = kFALSE;
      AliExternalTrackParam *trackInnerC3 = new AliExternalTrackParam(*(track->GetInnerParam()));
      if(mcEvent && stack) 
      {
        do //artificial loop (once) to make the continue statements jump out of the MC part
        {
          multMCTrueTracks = GetMCTrueTrackMult(mcEvent,evtCuts,accCuts);
          //
          // global track
          //
          Int_t label = TMath::Abs(track->GetLabel()); 
          if (label >= mcStackSize) continue;
          particle = stack->Particle(label);
          if (!particle) continue;
          if(particle && particle->GetPDG() && particle->GetPDG()->Charge()!=0.)
          {
            particleMother = GetMother(particle,stack);
            mech = particle->GetUniqueID();
            isPrim = stack->IsPhysicalPrimary(label);
            isFromStrangess  = IsFromStrangeness(label,stack);
            isFromConversion = IsFromConversion(label,stack);
            isFromMaterial   = IsFromMaterial(label,stack);
          }

          //
          // TPC track
          //
          Int_t labelTPC = TMath::Abs(track->GetTPCLabel()); 
          if (labelTPC >= mcStackSize) continue;
          particleTPC = stack->Particle(labelTPC);
          if (!particleTPC) continue;
          if(particleTPC && particleTPC->GetPDG() && particleTPC->GetPDG()->Charge()!=0.)
          {
            particleMotherTPC = GetMother(particleTPC,stack);
            mechTPC = particleTPC->GetUniqueID();
            isPrimTPC = stack->IsPhysicalPrimary(labelTPC);
            isFromStrangessTPC  = IsFromStrangeness(labelTPC,stack);
            isFromConversionTPC = IsFromConversion(labelTPC,stack);
            isFromMaterialTPC   = IsFromMaterial(labelTPC,stack);
          }

          //
          // store first track reference
          // for TPC track
          //
          TParticle *part=0;
          TClonesArray *trefs=0;
          Int_t status = mcEvent->GetParticleAndTR(TMath::Abs(labelTPC), part, trefs);

          if(status>0 && part && trefs && part->GetPDG() && part->GetPDG()->Charge()!=0.) 
          {
            Int_t nTrackRef = trefs->GetEntries();
            //printf("nTrackRef %d \n",nTrackRef);

            Int_t countITS = 0;
            for (Int_t iref = 0; iref < nTrackRef; iref++) 
            {
              AliTrackReference *ref = (AliTrackReference *)trefs->At(iref);

              // Ref. in the middle ITS 
              if(ref && ref->Label()==label && ref->DetectorId()==AliTrackReference::kITS)
              {
                nrefITS++;
                if(!refITS && countITS==2) {
                  refITS = ref;
                  //printf("refITS %p \n",refITS);
                }
                countITS++;
              }

              // TPC
              if(ref && ref->Label()==label  && ref->DetectorId()==AliTrackReference::kTPC)
              {
                nrefTPC++;
                refTPCOut=ref;
                if(!refTPCIn) {
                  refTPCIn = ref;
                  //printf("refTPCIn %p \n",refTPCIn);
                  //break;
                }
              }
              // TRD
              if(ref && ref->Label()==label && ref->DetectorId()==AliTrackReference::kTRD)
              {
                nrefTRD++;
                if(!refTRD) {
                  refTRD = ref;
                }
              }
              // TOF
              if(ref && ref->Label()==label  && ref->DetectorId()==AliTrackReference::kTOF)
              {
                nrefTOF++;
                if(!refTOF) {
                  refTOF = ref;
                }
              }
              // EMCAL
              if(ref && ref->Label()==label  && ref->DetectorId()==AliTrackReference::kEMCAL)
              {
                nrefEMCAL++;
                if(!refEMCAL) {
                  refEMCAL = ref;
                }
              }
              // PHOS
              //            if(ref && ref->Label()==label  && ref->DetectorId()==AliTrackReference::kPHOS)
              //             {
              // 	      nrefPHOS++;
              // 	      if(!refPHOS) {
              // 	        refPHOS = ref;
              // 	      }
              //             }
            }

            // transform inner params to TrackRef
            // reference frame
            if(refTPCIn && trackInnerC3) 
            {
              Double_t kRefPhi = TMath::ATan2(refTPCIn->Y(),refTPCIn->X());
              isOKtrackInnerC3 = trackInnerC3->Rotate(kRefPhi);
              isOKtrackInnerC3 = AliTracker::PropagateTrackToBxByBz(trackInnerC3,refTPCIn->R(),track->GetMass(),kStep,kFALSE);
            }
          }

          //
          // ITS track
          //
          Int_t labelITS = TMath::Abs(track->GetITSLabel()); 
          if (labelITS >= mcStackSize) continue;
          particleITS = stack->Particle(labelITS);
          if (!particleITS) continue;
          if(particleITS && particleITS->GetPDG() && particleITS->GetPDG()->Charge()!=0.)
          {
            particleMotherITS = GetMother(particleITS,stack);
            mechITS = particleITS->GetUniqueID();
            isPrimITS = stack->IsPhysicalPrimary(labelITS);
            isFromStrangessITS  = IsFromStrangeness(labelITS,stack);
            isFromConversionITS = IsFromConversion(labelITS,stack);
            isFromMaterialITS   = IsFromMaterial(labelITS,stack);
          }
        }
        while (0);
      }

        //
        Bool_t dumpToTree=kFALSE;

        if(isOKtpcInnerC  && isOKtrackInnerC) dumpToTree = kTRUE;
        //if(fUseESDfriends && isOKtrackInnerC2 && isOKouterITSc) dumpToTree = kTRUE;
        if(isOKtrackInnerC2 && isOKouterITSc) dumpToTree = kTRUE;
        if(mcEvent && isOKtrackInnerC3) dumpToTree = kTRUE;
        TObjString triggerClass = esdEvent->GetFiredTriggerClasses().Data();
        if (fReducePileUp){  
          //
          // 18.03 - Reduce pile-up chunks, done outside of the ESDTrackCuts for 2012/2013 data pile-up about 95 % of tracks
          // Done only in case no MC info 
          //
          Float_t dcaTPC[2];
          track->GetImpactParametersTPC(dcaTPC[0],dcaTPC[1]);
          Bool_t isRoughPrimary = TMath::Abs(dcaTPC[1])<10;
          Bool_t hasOuter=(track->IsOn(AliVTrack::kITSin))||(track->IsOn(AliVTrack::kTOFout))||(track->IsOn(AliVTrack::kTRDin));
          Bool_t keepPileUp=gRandom->Rndm()<0.05;
          if ( (!hasOuter) && (!isRoughPrimary) && (!keepPileUp)){
            dumpToTree=kFALSE;
          }
        }

        //init dummy objects
        static AliESDVertex dummyvtxESD;
        //if (!dummyvtxESD) 
        //{
        //  dummyvtxESD=new AliESDVertex();
        //  //dummyvtxESD->SetNContributors(1);
        //  //UShort_t pindices[1]; pindices[0]=0;
        //  //dummyvtxESD->SetIndices(1,pindices);
        //  //dummyvtxESD->SetNContributors(0);
        //}
        static AliExternalTrackParam dummyexternaltrackparam;
        //if (!dummyexternaltrackparam) dummyexternaltrackparam=new AliExternalTrackParam();
        static AliTrackReference dummytrackreference;
        //if (!dummytrackreference) dummytrackreference=new AliTrackReference();
        static TParticle dummyparticle;
        //if (!dummyparticle) dummyparticle=new TParticle();

        //assign the dummy objects if needed
        if (!track) {track=fDummyTrack;}
	AliESDfriendTrack *friendTrackStore=friendTrack;    // store friend track for later processing
	if (fFriendDownscaling>=1){  // downscaling number of friend tracks
	  friendTrackStore = (gRandom->Rndm()<1./fFriendDownscaling)? friendTrack:0;
	}
	if (fFriendDownscaling<=0){
	  if (((*fTreeSRedirector)<<"highPt").GetTree()){
	    TTree * tree = ((*fTreeSRedirector)<<"highPt").GetTree();
	    if (tree){
	      Double_t sizeAll=tree->GetZipBytes();
	      TBranch * br= tree->GetBranch("friendTrack.fPoints");
	      Double_t sizeFriend=(br!=NULL)?br->GetZipBytes():0;
	      br= tree->GetBranch("friendTrack.fCalibContainer");
	      if (br) sizeFriend+=br->GetZipBytes();
	      if (sizeFriend*TMath::Abs(fFriendDownscaling)>sizeAll) friendTrackStore=0;
	    }
	  }
	}


	//        if (!friendTrackStore && fFriendDownscaling<=1) {friendTrack=fDummyFriendTrack;}
        if (!vtxESD) {vtxESD=&dummyvtxESD;}
        if (mcEvent)
        {
          if (!refTPCIn) {refTPCIn=&dummytrackreference;}
          if (!refITS) {refITS=&dummytrackreference;}
          if (!particle) {particle=&dummyparticle;}
          if (!particleMother) {particleMother=&dummyparticle;}
          if (!particleTPC) {particleTPC=&dummyparticle;}
          if (!particleMotherTPC) {particleMotherTPC=&dummyparticle;}
          if (!particleITS) {particleITS=&dummyparticle;}
          if (!particleMotherITS) {particleMotherITS=&dummyparticle;}
        }

     
        // fill histograms
        FillHistograms(track, tpcInnerC, centralityF, (Double_t)chi2(0,0));
	TVectorD tofClInfo(5);                        // starting at 2014 - TOF infdo not part of the AliESDtrack
	tofClInfo[0]=track->GetTOFsignal();
	tofClInfo[1]=track->GetTOFsignalToT();
	tofClInfo[2]=track->GetTOFsignalRaw();
	tofClInfo[3]=track->GetTOFsignalDz();
	tofClInfo[4]=track->GetTOFsignalDx();

	//get the nSigma information; NB particle number ID in the vectors follow the convention of AliPID
        const Int_t nSpecies=AliPID::kSPECIES;
	TVectorD tpcNsigma(nSpecies); 
        TVectorD tofNsigma(nSpecies);
	TVectorD tpcPID(nSpecies); // bayes
        TVectorD tofPID(nSpecies);
	if(pidResponse){
          for (Int_t ispecie=0; ispecie<nSpecies; ++ispecie) {
            if (ispecie == Int_t(AliPID::kMuon)) continue;
            tpcNsigma[ispecie] = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType)ispecie);
            tofNsigma[ispecie] = pidResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType)ispecie);
	    //
          }	
	  pidResponse->ComputePIDProbability(AliPIDResponse::kTPC, track, nSpecies, tpcPID.GetMatrixArray());
	  pidResponse->ComputePIDProbability(AliPIDResponse::kTOF, track, nSpecies, tofPID.GetMatrixArray());	    
	}
        if(fTreeSRedirector && dumpToTree && fFillTree) {
	  downscaleCounter++;
          (*fTreeSRedirector)<<"highPt"<<
	    "downscaleCounter="<<downscaleCounter<<   
            "gid="<<gid<<
            "fileName.="<<&fCurrentFileName<<                // name of the chunk file (hopefully full)
            "runNumber="<<runNumber<<                // runNumber
            "evtTimeStamp="<<evtTimeStamp<<          // time stamp of event (in seconds)
            "evtNumberInFile="<<evtNumberInFile<<    // event number
            "triggerClass="<<&triggerClass<<         // trigger class as a string
            "Bz="<<bz<<                              // solenoid magnetic field in the z direction (in kGaus)
            "vtxESD.="<<vtxESD<<                    // vertexer ESD tracks (can be biased by TPC pileup tracks)
            "IRtot="<<ir1<<                         // interaction record (trigger) counters - coutner 1
            "IRint2="<<ir2<<                        // interaction record (trigger) coutners - counter 2
            "mult="<<mult<<                         // multiplicity of tracks pointing to the primary vertex
            "ntracks="<<ntracks<<                   // number of the esd tracks (to take into account the pileup in the TPC)
            //                                           important variables for the pile-up studies
            "contTPC="<< contTPC<<                    // number of contributors to the TPC primary vertex candidate
            "contSPD="<< contSPD<<                    // number of contributors to the SPD primary vertex candidate
            "vertexPosTPC.="<<&vertexPosTPC<<          // TPC vertex position
            "vertexPosSPD.="<<&vertexPosSPD<<          // SPD vertex position
            "ntracksTPC="<<ntracksTPC<<               // total number of the TPC tracks which were refitted
            "ntracksITS="<<ntracksITS<<               // total number of the ITS tracks which were refitted
            //
            "esdTrack.="<<track<<                  // esdTrack as used in the physical analysis
	    "tofClInfo.="<<&tofClInfo<<           // tof info
	    //            "friendTrack.="<<friendTrack<<      // esdFriendTrack associated to the esdTrack
	    "tofNsigma.="<<&tofNsigma<<
	    "tpcNsigma.="<<&tpcNsigma<<
	    "tofPID.="<<&tofPID<<                  // bayesian PID - without priors
	    "tpcPID.="<<&tpcPID<<                  // bayesian PID - without priors
	    
	    "friendTrack.="<<friendTrackStore<<      // esdFriendTrack associated to the esdTrack 
            "extTPCInnerC.="<<tpcInnerC<<          // TPC track from the first tracking iteration propagated and updated at vertex 
            "extInnerParamV.="<<trackInnerV<<      // TPC+TRD  inner param after refit  propagate to vertex 
            "extInnerParamC.="<<trackInnerC<<      // TPC+TRD  inner param after refit  propagate and updated at vertex 
            "extInnerParam.="<<trackInnerC2<<      // TPC+TRD  inner param after refit propagate to refernce TPC layer
            "extOuterITS.="<<outerITSc<<           // ITS outer track propagated to the TPC refernce radius
            "extInnerParamRef.="<<trackInnerC3<<   // TPC+TRD  inner param after refit propagated to the first TPC reference 
            "chi2TPCInnerC="<<chi2(0,0)<<           // chi2   of tracks ???
            "chi2InnerC="<<chi2trackC(0,0)<<        // chi2s  of tracks TPCinner to the combined
            "chi2OuterITS="<<chi2OuterITS(0,0)<<    // chi2s  of tracks TPC at inner wall to the ITSout
            "centralityF="<<centralityF;
	  // info for 2 track resolution studies and matching efficency studies 
	  //
	  (*fTreeSRedirector)<<"highPt"<<
	    "paramITS.="<<&paramITS<<                // nearest ITS track  -   chi2 distance at vertex
	    "paramITSC.="<<&paramITSC<<              // nearest ITS track  -  to constrained track   chi2 distance at vertex
	    "paramComb.="<<&paramComb<<              // nearest comb. tack -   chi2 distance at inner wall
	    "indexNearestITS="<<indexNearestITS<<    // index of  nearest ITS track
	    "indexNearestITSC="<<indexNearestITSC<<  // index of  nearest ITS track for constrained track
	    "indexNearestComb="<<indexNearestComb;   // index of  nearest track for constrained track

          if (mcEvent){
            static AliTrackReference refDummy;
            if (!refITS) refITS = &refDummy;
            if (!refTRD) refTRD = &refDummy;
            if (!refTOF) refTOF = &refDummy;
            if (!refEMCAL) refEMCAL = &refDummy;
            if (!refPHOS) refPHOS = &refDummy;
	    downscaleCounter++;
            (*fTreeSRedirector)<<"highPt"<<	
              "multMCTrueTracks="<<multMCTrueTracks<<   // mC track multiplicities
              "nrefITS="<<nrefITS<<              // number of track references in the ITS
              "nrefTPC="<<nrefTPC<<              // number of track references in the TPC
              "nrefTRD="<<nrefTRD<<              // number of track references in the TRD
              "nrefTOF="<<nrefTOF<<              // number of track references in the TOF
              "nrefEMCAL="<<nrefEMCAL<<              // number of track references in the TOF
              "nrefPHOS="<<nrefPHOS<<              // number of track references in the TOF
              "refTPCIn.="<<refTPCIn<<
              "refTPCOut.="<<refTPCOut<<
              "refITS.="<<refITS<<	    
              "refTRD.="<<refTRD<<	    
              "refTOF.="<<refTOF<<	    
              "refEMCAL.="<<refEMCAL<<	    
              "refPHOS.="<<refPHOS<<	    
              "particle.="<<particle<<
              "particleMother.="<<particleMother<<
              "mech="<<mech<<
              "isPrim="<<isPrim<<
              "isFromStrangess="<<isFromStrangess<<
              "isFromConversion="<<isFromConversion<<
              "isFromMaterial="<<isFromMaterial<<
              "particleTPC.="<<particleTPC<<
              "particleMotherTPC.="<<particleMotherTPC<<
              "mechTPC="<<mechTPC<<
              "isPrimTPC="<<isPrimTPC<<
              "isFromStrangessTPC="<<isFromStrangessTPC<<
              "isFromConversionTPC="<<isFromConversionTPC<<
              "isFromMaterialTPC="<<isFromMaterialTPC<<
              "particleITS.="<<particleITS<<
              "particleMotherITS.="<<particleMotherITS<<
              "mechITS="<<mechITS<<
              "isPrimITS="<<isPrimITS<<
              "isFromStrangessITS="<<isFromStrangessITS<<
              "isFromConversionITS="<<isFromConversionITS<<
              "isFromMaterialITS="<<isFromMaterialITS;
          }
          //finish writing the entry
          AliInfo("writing tree highPt");
          (*fTreeSRedirector)<<"highPt"<<"\n";
        }
        AliSysInfo::AddStamp("filteringTask",iTrack,numberOfTracks,numberOfFriendTracks,(friendTrackStore)?0:1);
        delete tpcInnerC;
        delete trackInnerC;
        delete trackInnerC2;
        delete outerITSc;
        delete trackInnerC3;
    }
  }
}


//_____________________________________________________________________________
void AliAnalysisTaskFilteredTree::ProcessMCEff(AliESDEvent *const esdEvent, AliMCEvent * const mcEvent, AliESDfriend *const /*esdFriend*/)
{
  //
  // Fill tree for efficiency studies MC only
  static Int_t downscaleCounter=0;
  AliInfo("we start!");
  if(!mcEvent) {
    AliDebug(AliLog::kError, "mcEvent not available");
    return;
  }

  // get selection cuts
  AliFilteredTreeEventCuts *evtCuts = GetEventCuts(); 
  AliFilteredTreeAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 

  if(!evtCuts || !accCuts  || !esdTrackCuts) {
    AliDebug(AliLog::kError, "cuts not available");
    return;
  }

  // trigger selection
  Bool_t isEventTriggered = kTRUE;
  AliPhysicsSelection *physicsSelection = NULL;
  AliTriggerAnalysis* triggerAnalysis = NULL;

  // 
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();


  // trigger
  if(evtCuts->IsTriggerRequired())  
  {
    // always MB
    isEventTriggered = inputHandler->IsEventSelected() & AliVEvent::kMB;
    physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if(!physicsSelection) return;

    if (isEventTriggered && (GetTrigger() == AliTriggerAnalysis::kV0AND)) {
      // set trigger (V0AND)
      triggerAnalysis = physicsSelection->GetTriggerAnalysis();
      if(!triggerAnalysis) return;
      isEventTriggered = triggerAnalysis->IsOfflineTriggerFired(esdEvent, GetTrigger());
    }
  }

  // centrality determination
  Float_t centralityF = -1;
  AliCentrality *esdCentrality = esdEvent->GetCentrality();
  centralityF = esdCentrality->GetCentralityPercentile(fCentralityEstimator.Data());

  // use MC information
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  Int_t mcStackSize=0;
  TArrayF vtxMC(3);

  Int_t multMCTrueTracks = 0;
  //
  if(!mcEvent) {
    AliDebug(AliLog::kError, "mcEvent not available");
    return;
  }
  // get MC event header
  header = mcEvent->Header();
  if (!header) {
    AliDebug(AliLog::kError, "Header not available");
    return;
  }
  // MC particle stack
  stack = mcEvent->Stack();
  if (!stack) {
    AliDebug(AliLog::kError, "Stack not available");
    return;
  }
  mcStackSize=stack->GetNtrack();

  // get MC vertex
  genHeader = header->GenEventHeader();
  if (!genHeader) {
    AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
    return;
  }
  genHeader->PrimaryVertex(vtxMC);

  // multipliticy of all MC primary tracks
  // in Zv, pt and eta ranges)
  multMCTrueTracks = GetMCTrueTrackMult(mcEvent,evtCuts,accCuts);


  // get reconstructed vertex  
  //const AliESDVertex* vtxESD = 0; 
  AliESDVertex* vtxESD = 0; 
  if(GetAnalysisMode() == kTPCAnalysisMode) {
    vtxESD = (AliESDVertex*)esdEvent->GetPrimaryVertexTPC();
  }
  else if(GetAnalysisMode() == kTPCITSAnalysisMode) {
    vtxESD = (AliESDVertex*)esdEvent->GetPrimaryVertexTracks();
  }
  else {
    return;
  }

  if(!vtxESD) return;   
  //
  // Vertex info comparison and track multiplicity
  //
  AliESDVertex *vertexSPD =  (AliESDVertex *)esdEvent->GetPrimaryVertexSPD();
  AliESDVertex *vertexTPC =  (AliESDVertex *)esdEvent->GetPrimaryVertexTPC(); 
  Int_t contSPD = vertexSPD->GetNContributors();
  Int_t contTPC = vertexTPC->GetNContributors();        
  TVectorD vertexPosTPC(3), vertexPosSPD(3);
  vertexSPD->GetXYZ(vertexPosSPD.GetMatrixArray());
  vertexTPC->GetXYZ(vertexPosTPC.GetMatrixArray());
  Int_t ntracksTPC=0;
  Int_t ntracksITS=0;
  for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++){
    AliESDtrack *track = esdEvent->GetTrack(iTrack);    
    if(!track) continue;
    if (track->IsOn(AliVTrack::kTPCrefit)) ntracksTPC++;
    if (track->IsOn(AliVTrack::kITSrefit)) ntracksITS++;
  }

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD); 
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK, isEventTriggered);
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());

  TObjString triggerClass = esdEvent->GetFiredTriggerClasses().Data();

  // check event cuts
  if(isEventOK && isEventTriggered)
  {
    //if(!stack) return;

    //
    // MC info
    //
    TParticle *particle=NULL;
    TParticle *particleMother=NULL;
    Int_t mech=-1;

    // reco event info
    Double_t vert[3] = {0}; 
    vert[0] = vtxESD->GetX();
    vert[1] = vtxESD->GetY();
    vert[2] = vtxESD->GetZ();
    Int_t mult = vtxESD->GetNContributors();
    Double_t bz = esdEvent->GetMagneticField();
    Double_t runNumber = esdEvent->GetRunNumber();
    Double_t evtTimeStamp = esdEvent->GetTimeStamp();
    Int_t evtNumberInFile = esdEvent->GetEventNumberInFile();
    // loop over MC stack
    for (Int_t iMc = 0; iMc < mcStackSize; ++iMc) 
    {
      particle = stack->Particle(iMc);
      if (!particle)
        continue;

      // only charged particles
      if(!particle->GetPDG()) continue;
      Double_t charge = particle->GetPDG()->Charge()/3.;
      if (TMath::Abs(charge) < 0.001)
        continue;

      // only primary particles
      Bool_t prim = stack->IsPhysicalPrimary(iMc);
      if(!prim) continue;

      // downscale low-pT particles
      Double_t scalempt= TMath::Min(particle->Pt(),10.);
      Double_t downscaleF = gRandom->Rndm();
      downscaleF *= fLowPtTrackDownscaligF;
      if (downscaleCounter>0 && TMath::Exp(2*scalempt)<downscaleF) continue;
      // is particle in acceptance
      if(!accCuts->AcceptTrack(particle)) continue;

      // check if particle reconstructed
      Bool_t isRec = kFALSE;
      Int_t trackIndex = -1;  
      Int_t trackLoopIndex = -1;  
      Int_t isESDtrackCut= 0;
      Int_t isAccCuts    = 0;
      Int_t nRec = 0;    // how many times reconstructed 
      Int_t nFakes = 0;  // how many times reconstructed as a fake track
      AliESDtrack *recTrack = NULL; 

      for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++)
      {
        AliESDtrack *track = esdEvent->GetTrack(iTrack);
        if(!track) continue;
        if(track->Charge()==0) continue;
        //
        Int_t label =  TMath::Abs(track->GetLabel());
        if (label >= mcStackSize) continue;
        if(label == iMc) {	  
          Bool_t isAcc=esdTrackCuts->AcceptTrack(track);
          if (isAcc) isESDtrackCut=1;
          if (accCuts->AcceptTrack(track)) isAccCuts=1;
          isRec = kTRUE;
          trackIndex = iTrack;

          if (recTrack){
            if (track->GetTPCncls()<recTrack->GetTPCncls()) continue; // in case looper tracks use longer track
            if (!isAcc) continue;
            trackLoopIndex = iTrack;
          }
          recTrack = esdEvent->GetTrack(trackIndex); 
          nRec++;
          if(track->GetLabel()<0) nFakes++;

          continue;
        }        
      }

      // Store information in the output tree
      if (trackLoopIndex>-1)  { 
        recTrack = esdEvent->GetTrack(trackLoopIndex); 
      } else if (trackIndex >-1) {
        recTrack = esdEvent->GetTrack(trackIndex); 
      } else {
        recTrack = fDummyTrack; 
      } 

      particleMother = GetMother(particle,stack);
      mech = particle->GetUniqueID();

      //MC particle track length
      Double_t tpcTrackLength = 0.;
      AliMCParticle *mcParticle = (AliMCParticle*) mcEvent->GetTrack(iMc);
      if(mcParticle) {
	Int_t counter=0;
        tpcTrackLength = mcParticle->GetTPCTrackLength(bz,0.05,counter,3.0);
      } 


      //
      if(fTreeSRedirector && fFillTree) {
	downscaleCounter++;
        (*fTreeSRedirector)<<"MCEffTree"<<
          "fileName.="<<&fCurrentFileName<<
          "triggerClass.="<<&triggerClass<<
          "runNumber="<<runNumber<<
          "evtTimeStamp="<<evtTimeStamp<<
          "evtNumberInFile="<<evtNumberInFile<<     // 
          "Bz="<<bz<<                               // magnetic field
          "vtxESD.="<<vtxESD<<                      // vertex info
          //
          "mult="<<mult<<                           // primary vertex 9whatewe found) multiplicity
          "multMCTrueTracks="<<multMCTrueTracks<<   // mC track multiplicities
          //                                           important variables for the pile-up studies
          "contTPC="<< contTPC<<                    // number of contributors to the TPC primary vertex candidate
          "contSPD="<< contSPD<<                    // number of contributors to the SPD primary vertex candidate
          "vertexPosTPC.="<<&vertexPosTPC<<          // TPC vertex position
          "vertexPosSPD.="<<&vertexPosSPD<<          // SPD vertex position
          "ntracksTPC="<<ntracksTPC<<               // total number of the TPC tracks which were refitted
          "ntracksITS="<<ntracksITS<<               // total number of the ITS tracks which were refitted
          //
          //
          "isAcc0="<<isESDtrackCut<<                // track accepted by ESD track cuts
          "isAcc1="<<isAccCuts<<                    // track accepted by acceptance cuts flag
          "esdTrack.="<<recTrack<<                  // reconstructed track (only the longest from the loopers)
          "isRec="<<isRec<<                         // track was reconstructed
          "tpcTrackLength="<<tpcTrackLength<<       // track length in the TPC r projection
          "particle.="<<particle<<                  // particle properties
          "particleMother.="<<particleMother<<      // particle mother
          "mech="<<mech<<                           // production mechanizm
          "nRec="<<nRec<<                           // how many times reconstruted
          "nFakes="<<nFakes<<                       // how many times reconstructed as a fake track
          "\n";
      }

      //if(trackIndex <0 && recTrack) delete recTrack; recTrack=0;
    }
  }

}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskFilteredTree::IsHighDeDxParticle(AliESDtrack * track) {
  //
  // check if particle is Z > 1 
  //
  if (track->GetTPCNcls() < 60) return kFALSE;
  Double_t mom = track->GetInnerParam()->GetP();
  if (mom < 0.2) return kFALSE; // protection against unexpected behavior of Aleph parameterization
  Float_t dca[2], bCov[3];
  track->GetImpactParameters(dca,bCov);
  //

  Double_t triggerDeDx = 4*AliExternalTrackParam::BetheBlochAleph((mom*2)/(0.938*3),1.0288,31.9806,5.04114e-11,2.13096,2.38541);

  if (track->GetTPCsignal() > triggerDeDx && track->GetTPCsignal()<1000 && TMath::Abs(dca[0])<3.) return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________
void AliAnalysisTaskFilteredTree::ProcessV0(AliESDEvent *const esdEvent, AliMCEvent * const mcEvent, AliESDfriend *const esdFriend)
{
  //
  // Select real events with V0 (K0s and Lambda and Gamma) high-pT candidates
  //
  static Int_t downscaleCounter=0;
  if(!esdEvent) {
    AliDebug(AliLog::kError, "esdEvent not available");
    return;
  }

  // get selection cuts
  AliFilteredTreeEventCuts *evtCuts = GetEventCuts(); 
  AliFilteredTreeAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 

  if(!evtCuts || !accCuts  || !esdTrackCuts) {
    AliDebug(AliLog::kError, "cuts not available");
    return;
  }

  // trigger selection
  Bool_t isEventTriggered = kTRUE;
  AliPhysicsSelection *physicsSelection = NULL;
  AliTriggerAnalysis* triggerAnalysis = NULL;

  // 
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  AliPIDResponse *pidResponse = inputHandler->GetPIDResponse();

  // trigger
  if(evtCuts->IsTriggerRequired())  
  {
    // always MB
    isEventTriggered = inputHandler->IsEventSelected() & AliVEvent::kMB;

    physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if(!physicsSelection) return;
    //SetPhysicsTriggerSelection(physicsSelection);

    if (isEventTriggered && (GetTrigger() == AliTriggerAnalysis::kV0AND)) {
      // set trigger (V0AND)
      triggerAnalysis = physicsSelection->GetTriggerAnalysis();
      if(!triggerAnalysis) return;
      isEventTriggered = triggerAnalysis->IsOfflineTriggerFired(esdEvent, GetTrigger());
    }
  }

  // centrality determination
  Float_t centralityF = -1;
  AliCentrality *esdCentrality = esdEvent->GetCentrality();
  centralityF = esdCentrality->GetCentralityPercentile(fCentralityEstimator.Data());


  // get reconstructed vertex  
  //const AliESDVertex* vtxESD = 0; 
  AliESDVertex* vtxESD = 0; 
  if(GetAnalysisMode() == kTPCAnalysisMode) {
    vtxESD = (AliESDVertex*)esdEvent->GetPrimaryVertexTPC();
  }
  else if(GetAnalysisMode() == kTPCITSAnalysisMode) {
    vtxESD = (AliESDVertex*)esdEvent->GetPrimaryVertexTracks();
  }
  else {
    return;
  }

  if(!vtxESD) return;

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD); 
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK, isEventTriggered);
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());
  Int_t ntracks = esdEvent->GetNumberOfTracks();
  // Global event id calculation using orbitID, bunchCrossingID and periodID 
  ULong64_t orbitID      = (ULong64_t)esdEvent->GetOrbitNumber();
  ULong64_t bunchCrossID = (ULong64_t)esdEvent->GetBunchCrossNumber();
  ULong64_t periodID     = (ULong64_t)esdEvent->GetPeriodNumber();
  ULong64_t gid          = ((periodID << 36) | (orbitID << 12) | bunchCrossID); 
  TObjString triggerClass = esdEvent->GetFiredTriggerClasses().Data();
  Float_t bz = esdEvent->GetMagneticField();
  Int_t run = esdEvent->GetRunNumber();
  Int_t time = esdEvent->GetTimeStamp();
  Int_t evtNumberInFile = esdEvent->GetEventNumberInFile();
  Int_t nV0s = esdEvent->GetNumberOfV0s();
  Int_t mult = vtxESD->GetNContributors();
  (*fTreeSRedirector)<<"eventInfoV0"<<
    "gid="<<gid<<
    "fileName.="<<&fCurrentFileName<<                // name of the chunk file (hopefully full)
    "run="<<run<<                             // runNumber
    "time="<<time<<                           // time stamp of event (in seconds)
    "evtNumberInFile="<<evtNumberInFile<<     // event number
    "triggerClass="<<&triggerClass<<          // trigger class as a string
    "Bz="<<bz<<                               // solenoid magnetic field in the z direction (in kGaus)
    "mult="<<mult<<                           // multiplicity of tracks pointing to the primary vertex
    "ntracks="<<ntracks<<                     // number of the esd tracks (to take into account the pileup in the TPC)
    "nV0s="<<nV0s<<                           // number of V0s
    "isEventOK="<<isEventOK<<                 // flag - AliFilteredTreeEventCuts - track dumped only for selected events 
    "isEventTriggered="<<isEventTriggered<<   // flag - if tigger required - track dumped only for selected events 
    "\n";




  // check event cuts
  if(isEventOK && isEventTriggered) {
    //
    // Dump the pt downscaled V0 into the tree
    // 
    Int_t ntracks = esdEvent->GetNumberOfTracks();
    Int_t evNr=esdEvent->GetEventNumberInFile();


    for (Int_t iv0=0; iv0<nV0s; iv0++){
 
      AliESDv0 * v0 = esdEvent->GetV0(iv0);
      if (!v0) continue;
      AliESDtrack * track0 = esdEvent->GetTrack(v0->GetIndex(0));
      AliESDtrack * track1 = esdEvent->GetTrack(v0->GetIndex(1));
      if (!track0) continue;
      if (!track1) continue;
      AliESDfriendTrack* friendTrack0=NULL;
      AliESDfriendTrack* friendTrack1=NULL;
      if (esdFriend)       {
        if (!esdFriend->TestSkipBit()){
	  Int_t ntracksFriend = esdFriend->GetNumberOfTracks();
	  if (v0->GetIndex(0)<ntracksFriend){
	    friendTrack0 = esdFriend->GetTrack(v0->GetIndex(0)); //this guy can be NULL
	  }
	  if (v0->GetIndex(1)<ntracksFriend){
	    friendTrack1 = esdFriend->GetTrack(v0->GetIndex(1)); //this guy can be NULL
	  }
        }
      }
      if (track0->GetSign()<0) {
        track1 = esdEvent->GetTrack(v0->GetIndex(0));
        track0 = esdEvent->GetTrack(v0->GetIndex(1));
      }

      //
      AliESDfriendTrack *friendTrackStore0=friendTrack0;    // store friend track0 for later processing
      AliESDfriendTrack *friendTrackStore1=friendTrack1;    // store friend track1 for later processing
      if (fFriendDownscaling>=1){  // downscaling number of friend tracks
	if (gRandom->Rndm()>1./fFriendDownscaling){
	  friendTrackStore0 = 0;
	  friendTrackStore1 = 0;
	}
      }
      if (fFriendDownscaling<=0){
	if (((*fTreeSRedirector)<<"V0s").GetTree()){
	  TTree * tree = ((*fTreeSRedirector)<<"V0s").GetTree();
	  if (tree){
	    Double_t sizeAll=tree->GetZipBytes();
	    TBranch * br= tree->GetBranch("friendTrack0.fPoints");
	    Double_t sizeFriend=(br!=NULL)?br->GetZipBytes():0;
	    br= tree->GetBranch("friendTrack0.fCalibContainer");
	    if (br) sizeFriend+=br->GetZipBytes();
	    if (sizeFriend*TMath::Abs(fFriendDownscaling)>sizeAll) {
	      friendTrackStore0=0;
	      friendTrackStore1=0;
	    }
	  }
	}
      }

      //
      Bool_t isDownscaled = IsV0Downscaled(v0);
      if (downscaleCounter>0 && isDownscaled) continue;
      AliKFParticle kfparticle; //
      Int_t type=GetKFParticle(v0,esdEvent,kfparticle);
      if (type==0) continue;   
      TObjString triggerClass = esdEvent->GetFiredTriggerClasses().Data();

      if(!fFillTree) return;
      if(!fTreeSRedirector) return;
      
      TVectorD tofClInfo0(5);                        // starting at 2014 - TOF infdo not part of the AliESDtrack
      TVectorD tofClInfo1(5);                        // starting at 2014 - TOF infdo not part of the AliESDtrack
      tofClInfo0[0]=track0->GetTOFsignal();
      tofClInfo0[1]=track0->GetTOFsignalToT();
      tofClInfo0[2]=track0->GetTOFsignalRaw();
      tofClInfo0[3]=track0->GetTOFsignalDz();
      tofClInfo0[4]=track0->GetTOFsignalDx();
      tofClInfo1[0]=track1->GetTOFsignal();
      tofClInfo1[1]=track1->GetTOFsignalToT();
      tofClInfo1[2]=track1->GetTOFsignalRaw();
      tofClInfo1[3]=track1->GetTOFsignalDz();
      tofClInfo1[4]=track1->GetTOFsignalDx();
 
      //get the nSigma information; NB particle number ID in the vectors follow the convention in AliPID
      const Int_t nSpecies=AliPID::kSPECIES;
      TVectorD tpcNsigma0(nSpecies); 
      TVectorD tofNsigma0(nSpecies);
      TVectorD tpcNsigma1(nSpecies); 
      TVectorD tofNsigma1(nSpecies);
      if(pidResponse){
        for (Int_t ispecie=0; ispecie<nSpecies; ++ispecie) {
          if (ispecie == Int_t(AliPID::kMuon)) continue;
          tpcNsigma0[ispecie] = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC, track0, (AliPID::EParticleType)ispecie);
          tofNsigma0[ispecie] = pidResponse->NumberOfSigmas(AliPIDResponse::kTOF, track0, (AliPID::EParticleType)ispecie);
          tpcNsigma1[ispecie] = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC, track1, (AliPID::EParticleType)ispecie);
          tofNsigma1[ispecie] = pidResponse->NumberOfSigmas(AliPIDResponse::kTOF, track1, (AliPID::EParticleType)ispecie);
        }
      }

      downscaleCounter++;
      (*fTreeSRedirector)<<"V0s"<<
        "gid="<<gid<<                         //  global id of event
        "isDownscaled="<<isDownscaled<<       //  
        "triggerClass="<<&triggerClass<<      //  trigger
        "Bz="<<bz<<                           //
        "fileName.="<<&fCurrentFileName<<     //  full path - file name with ESD
        "runNumber="<<run<<                   //
        "evtTimeStamp="<<time<<               //  time stamp of event in secons
        "evtNumberInFile="<<evNr<<            //  
        "type="<<type<<                       // type of V0-
        "ntracks="<<ntracks<<
        "v0.="<<v0<<
        "kf.="<<&kfparticle<<
        "track0.="<<track0<<                  // track
        "track1.="<<track1<<
	"tofClInfo0.="<<&tofClInfo0<<
	"tofClInfo1.="<<&tofClInfo1<<
      	"tofNsigma0.="<<&tofNsigma0<<
	"tofNsigma1.="<<&tofNsigma1<<
      	"tpcNsigma0.="<<&tpcNsigma0<<
	"tpcNsigma1.="<<&tpcNsigma1<<
        "friendTrack0.="<<friendTrackStore0<<
        "friendTrack1.="<<friendTrackStore1<<
        "centralityF="<<centralityF<<
        "\n";
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskFilteredTree::ProcessdEdx(AliESDEvent *const esdEvent, AliMCEvent * const mcEvent, AliESDfriend *const esdFriend)
{
  //
  // Select real events with large TPC dEdx signal
  //
  static Int_t downscaleCounter=0;
  if(!esdEvent) {
    AliDebug(AliLog::kError, "esdEvent not available");
    return;
  }

  // get selection cuts
  AliFilteredTreeEventCuts *evtCuts = GetEventCuts(); 
  AliFilteredTreeAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 

  if(!evtCuts || !accCuts  || !esdTrackCuts) {
    AliDebug(AliLog::kError, "cuts not available");
    return;
  }


  // trigger
  Bool_t isEventTriggered = kTRUE;
  AliPhysicsSelection *physicsSelection = NULL;
  AliTriggerAnalysis* triggerAnalysis = NULL;

  // 
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  AliPIDResponse *pidResponse = inputHandler->GetPIDResponse();

  if(evtCuts->IsTriggerRequired())  
  {
    // always MB
    isEventTriggered = inputHandler->IsEventSelected() & AliVEvent::kMB;

    physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if(!physicsSelection) return;

    if (isEventTriggered && (GetTrigger() == AliTriggerAnalysis::kV0AND)) {
      // set trigger (V0AND)
      triggerAnalysis = physicsSelection->GetTriggerAnalysis();
      if(!triggerAnalysis) return;
      isEventTriggered = triggerAnalysis->IsOfflineTriggerFired(esdEvent, GetTrigger());
    }
  }

  // get reconstructed vertex  
  AliESDVertex* vtxESD = 0; 
  if(GetAnalysisMode() == kTPCAnalysisMode) {
    vtxESD = (AliESDVertex*)esdEvent->GetPrimaryVertexTPC();
  }
  else if(GetAnalysisMode() == kTPCITSAnalysisMode) {
    vtxESD = (AliESDVertex*)esdEvent->GetPrimaryVertexTracks();
  }
  else {
    return;
  }
  if(!vtxESD) return;

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD); 
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK, isEventTriggered);
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());


  // check event cuts
  if(isEventOK && isEventTriggered)
  {
    Double_t vert[3] = {0}; 
    vert[0] = vtxESD->GetX();
    vert[1] = vtxESD->GetY();
    vert[2] = vtxESD->GetZ();
    Int_t mult = vtxESD->GetNContributors();
    Double_t bz = esdEvent->GetMagneticField();
    Double_t runNumber = esdEvent->GetRunNumber();
    Double_t evtTimeStamp = esdEvent->GetTimeStamp();
    Int_t evtNumberInFile = esdEvent->GetEventNumberInFile();
   
    // Global event id calculation using orbitID, bunchCrossingID and periodID 
    ULong64_t orbitID      = (ULong64_t)esdEvent->GetOrbitNumber();
    ULong64_t bunchCrossID = (ULong64_t)esdEvent->GetBunchCrossNumber();
    ULong64_t periodID     = (ULong64_t)esdEvent->GetPeriodNumber();
    ULong64_t gid          = ((periodID << 36) | (orbitID << 12) | bunchCrossID); 
    
    // large dEdx
    for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++)
    {
      AliESDtrack *track = esdEvent->GetTrack(iTrack);
      if(!track) continue;
      AliESDfriendTrack* friendTrack=NULL;
      if (esdFriend && !esdFriend->TestSkipBit()) {
	Int_t ntracksFriend = esdFriend->GetNumberOfTracks();
	if (iTrack<ntracksFriend){
	  friendTrack = esdFriend->GetTrack(iTrack);
	} //this guy can be NULL	
      }
      
      if(track->Charge()==0) continue;
      if(!esdTrackCuts->AcceptTrack(track)) continue;
      if(!accCuts->AcceptTrack(track)) continue;

      if(!IsHighDeDxParticle(track)) continue;
      TObjString triggerClass = esdEvent->GetFiredTriggerClasses().Data();

      if(!fFillTree) return;
      if(!fTreeSRedirector) return;


      //get the nSigma information; NB particle number ID in the vectors follow the convention of AliPID
      const Int_t nSpecies=AliPID::kSPECIES;
      TVectorD tpcNsigma(nSpecies); 
      TVectorD tofNsigma(nSpecies);
      if(pidResponse){
        for (Int_t ispecie=0; ispecie<nSpecies; ++ispecie) {
          if (ispecie == Int_t(AliPID::kMuon)) continue;
          tpcNsigma[ispecie] = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType)ispecie);
          tofNsigma[ispecie] = pidResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType)ispecie);
        }
      }
	
      downscaleCounter++;
      (*fTreeSRedirector)<<"dEdx"<<           // high dEdx tree
        "gid="<<gid<<                         // global id
        "fileName.="<<&fCurrentFileName<<     // file name
        "runNumber="<<runNumber<<
        "evtTimeStamp="<<evtTimeStamp<<
        "evtNumberInFile="<<evtNumberInFile<<
        "triggerClass="<<&triggerClass<<      //  trigger
        "Bz="<<bz<<
        "vtxESD.="<<vtxESD<<                  // 
        "mult="<<mult<<
        "esdTrack.="<<track<<
        "friendTrack.="<<friendTrack<<
        "tofNsigma.="<<&tofNsigma<<
        "tpcNsigma.="<<&tpcNsigma<<
        "\n";
    }
  }
}

//_____________________________________________________________________________
Int_t   AliAnalysisTaskFilteredTree::GetKFParticle(AliESDv0 *const v0, AliESDEvent * const event, AliKFParticle & kfparticle)
{
  //
  // Create KF particle in case the V0 fullfill selection criteria
  //
  // Selection criteria
  //  0. algorithm cut
  //  1. track cut
  //  3. chi2 cut
  //  4. rough mass cut
  //  5. Normalized pointing angle cut
  //
  const Double_t cutMass=0.2;
  const Double_t kSigmaDCACut=3;
  //
  // 0.) algo cut - accept only on the fly
  //
  if (v0->GetOnFlyStatus() ==kFALSE) return 0;     
  //
  // 1.) track cut
  // 
  AliESDtrack * track0 = event->GetTrack(v0->GetIndex(0));
  AliESDtrack * track1 = event->GetTrack(v0->GetIndex(1));
  /*
     TCut cutD="abs(track0.fD/sqrt(track0.fCdd))>2&&abs(track1.fD/sqrt(track1.fCdd))>2";
     TCut cutTheta="abs(track0.fP[3])<1&&abs(track1.fP[3])<1";
     TCut cutNcl="track0.GetTPCClusterInfo(2,1)>100&&track1.GetTPCClusterInfo(2,1)>100";
     */  
  if (TMath::Abs(track0->GetTgl())>1) return 0;
  if (TMath::Abs(track1->GetTgl())>1) return 0;
  if ((track0->GetTPCClusterInfo(2,1))<100) return 0;
  if ((track1->GetTPCClusterInfo(2,1))<100) return 0;
  Float_t pos0[2]={0}, cov0[3]={0};
  Float_t pos1[2]={0}, cov1[3]={0};
  track0->GetImpactParameters(pos0,cov0);
  track0->GetImpactParameters(pos1,cov1);
  //
  if (TMath::Abs(pos0[0])<kSigmaDCACut*TMath::Sqrt(cov0[0])) return 0;
  if (TMath::Abs(pos1[0])<kSigmaDCACut*TMath::Sqrt(cov1[0])) return 0;
  // 
  //
  // 3.) Chi2 cut
  //
  Double_t chi2KF = v0->GetKFInfo(2,2,2);
  if (chi2KF>25) return 0;
  //
  // 4.) Rough mass cut - 0.200 GeV
  //
  static Double_t masses[2]={-1};
  if (masses[0]<0){
    masses[0] = TDatabasePDG::Instance()->GetParticle("K_S0")->Mass();
    masses[1] = TDatabasePDG::Instance()->GetParticle("Lambda0")->Mass();
  }
  Double_t mass00=  v0->GetEffMass(0,0);
  Double_t mass22=  v0->GetEffMass(2,2);
  Double_t mass42=  v0->GetEffMass(4,2);
  Double_t mass24=  v0->GetEffMass(2,4);
  Bool_t massOK=kFALSE;
  Int_t type=0;
  Int_t ptype=0;
  Double_t dmass=1;
  Int_t p1=0, p2=0;
  if (TMath::Abs(mass00-0)<cutMass) {
    massOK=kTRUE; type+=1; 
    if (TMath::Abs(mass00-0)<dmass) {
      ptype=1;
      dmass=TMath::Abs(mass00-0);      
      p1=0; p2=0;
    } 
  }
  if (TMath::Abs(mass24-masses[1])<cutMass) {
    massOK=kTRUE; type+=2; 
    if (TMath::Abs(mass24-masses[1])<dmass){
      dmass = TMath::Abs(mass24-masses[1]);
      ptype=2;
      p1=2; p2=4;
    }
  }
  if (TMath::Abs(mass42-masses[1])<cutMass) {
    massOK=kTRUE; type+=4;
    if (TMath::Abs(mass42-masses[1])<dmass){
      dmass = TMath::Abs(mass42-masses[1]);
      ptype=4;
      p1=4; p2=2;
    }
  }
  if (TMath::Abs(mass22-masses[0])<cutMass) {
    massOK=kTRUE; type+=8;
    if (TMath::Abs(mass22-masses[0])<dmass){
      dmass = TMath::Abs(mass22-masses[0]);
      ptype=8;
      p1=2; p2=2;
    }
  }
  if (type==0) return 0;
  //
  const Int_t spdg[5]={kPositron,kMuonPlus,kPiPlus, kKPlus, kProton};
  const AliExternalTrackParam *paramP = v0->GetParamP();
  const AliExternalTrackParam *paramN = v0->GetParamN();
  if (paramP->GetSign()<0){
    paramP=v0->GetParamP();
    paramN=v0->GetParamN();
  }
  //Double_t *pparam1 = (Double_t*)paramP->GetParameter();
  //Double_t *pparam2 = (Double_t*)paramN->GetParameter();
  //
  AliKFParticle kfp1( *paramP, spdg[p1]  );
  AliKFParticle kfp2( *paramN, -1 *spdg[p2]  );
  AliKFParticle V0KF;
  (V0KF)+=kfp1;
  (V0KF)+=kfp2;
  kfparticle=V0KF;
  //
  // Pointing angle
  //
  Double_t  errPhi    = V0KF.GetErrPhi();
  Double_t  pointAngle= TMath::ACos(v0->GetV0CosineOfPointingAngle());
  if (pointAngle/errPhi>10) return 0;  
  //
  return ptype;  
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskFilteredTree::IsV0Downscaled(AliESDv0 *const v0)
{
  //
  // Downscale randomly low pt V0
  //
  //return kFALSE;
  Double_t maxPt= TMath::Max(v0->GetParamP()->Pt(), v0->GetParamN()->Pt());
  Double_t scalempt= TMath::Min(maxPt,10.);
  Double_t downscaleF = gRandom->Rndm();
  downscaleF *= fLowPtV0DownscaligF;
  //
  // Special treatment of the gamma conversion pt spectra is softer - 
  Double_t mass00=  v0->GetEffMass(0,0);
  const Double_t cutMass=0.2;
  if (TMath::Abs(mass00-0)<cutMass){
    downscaleF/=10.;  // 10 times smaller downscaling for the gamma concersion candidate
  }
  //printf("V0 TMath::Exp(2*scalempt) %e, downscaleF %e \n",TMath::Exp(2*scalempt), downscaleF);
  if (TMath::Exp(2*scalempt)<downscaleF) return kTRUE;
  return kFALSE;

  /*
     TH1F his1("his1","his1",100,0,10);
     TH1F his2("his2","his2",100,0,10);
     {for (Int_t i=0; i<10000; i++){
     Double_t rnd=gRandom->Exp(1);
     Bool_t isDownscaled =TMath::Exp(rnd)<100*gRandom->Rndm();
     his1->Fill(rnd); 
     if (!isDownscaled) his2->Fill(rnd); 
     }}

*/

}



//_____________________________________________________________________________
Bool_t AliAnalysisTaskFilteredTree::ConstrainTPCInner(AliExternalTrackParam *const tpcInnerC, const AliESDVertex* vtx, Double_t b[3])
{
  // Constrain TPC inner params constrained
  //
  if(!tpcInnerC) return kFALSE; 
  if(!vtx) return kFALSE;

  Double_t dz[2],cov[3];
  //AliESDVertex *vtx= (AliESDVertex *)esdEvent->GetPrimaryVertex();
  //if(!tpcInnerC->PropagateToDCA(vtx, esdEvent->GetMagneticField(), 3, dz, cov)) return kFALSE; 
  //if(!tpcInnerC->PropagateToDCA(vtx, Bz, 3, dz, cov)) return kFALSE; 
  if(!tpcInnerC->PropagateToDCABxByBz(vtx, b, 3, dz, cov)) return kFALSE; 


  Double_t covar[6]; vtx->GetCovMatrix(covar);
  Double_t p[2]={tpcInnerC->GetParameter()[0]-dz[0],tpcInnerC->GetParameter()[1]-dz[1]};
  Double_t c[3]={covar[2],0.,covar[5]};
  Double_t chi2C=tpcInnerC->GetPredictedChi2(p,c);
  if (chi2C>kVeryBig) return kFALSE; 

  if(!tpcInnerC->Update(p,c)) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskFilteredTree::ConstrainTrackInner(AliExternalTrackParam *const trackInnerC, const AliESDVertex* vtx, Double_t mass, Double_t b[3])
{
  // Constrain TPC inner params constrained
  //
  if(!trackInnerC) return kFALSE; 
  if(!vtx) return kFALSE;

  const Double_t kRadius  = 2.8; 
  const Double_t kMaxStep = 1.0; 

  Double_t dz[2],cov[3];

  //AliESDVertex *vtx= (AliESDVertex *)esdEvent->GetPrimaryVertex();
  //if(!trackInnerC->PropagateToDCA(vtx, esdEvent->GetMagneticField(), 3, dz, cov)) return kFALSE; 
  //if(!trackInnerC->PropagateToDCA(vtx, Bz, 3, dz, cov)) return kFALSE; 

  if(!AliTracker::PropagateTrackToBxByBz(trackInnerC,kRadius,mass,kMaxStep,kFALSE)) return kFALSE;
  if(!trackInnerC->PropagateToDCABxByBz(vtx, b, 3, dz, cov)) return kFALSE; 

  //
  Double_t covar[6]; vtx->GetCovMatrix(covar);
  Double_t p[2]={trackInnerC->GetParameter()[0]-dz[0],trackInnerC->GetParameter()[1]-dz[1]};
  Double_t c[3]={covar[2],0.,covar[5]};
  Double_t chi2C=trackInnerC->GetPredictedChi2(p,c);
  if (chi2C>kVeryBig) return kFALSE; 

  if(!trackInnerC->Update(p,c)) return kFALSE;

  return kTRUE;
}


//_____________________________________________________________________________
TParticle *AliAnalysisTaskFilteredTree::GetMother(TParticle *const particle, AliStack *const stack) 
{
  if(!particle) return NULL;
  if(!stack) return NULL;

  Int_t motherLabel = TMath::Abs(particle->GetMother(0));  
  TParticle* mother = NULL; 
  Int_t mcStackSize=stack->GetNtrack();
  if (motherLabel>=mcStackSize) return NULL;
  mother = stack->Particle(motherLabel); 

  return mother;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskFilteredTree::IsFromConversion(Int_t label, AliStack *const stack) 
{
  Bool_t isFromConversion = kFALSE;

  if(stack) {
    Int_t mcStackSize=stack->GetNtrack();
    if (label>=mcStackSize) return kFALSE;
    TParticle* particle = stack->Particle(label);
    if (!particle) return kFALSE;

    if(particle && particle->GetPDG() && particle->GetPDG()->Charge()!=0) 
    {
      Int_t mech = particle->GetUniqueID(); // production mechanism 
      Bool_t isPrim = stack->IsPhysicalPrimary(label);

      Int_t motherLabel = TMath::Abs(particle->GetMother(0));  
      if (motherLabel>=mcStackSize) return kFALSE;
      TParticle* mother = stack->Particle(motherLabel); 
      if(mother) {
        Int_t motherPdg = mother->GetPdgCode();

        if(!isPrim && mech==5 && motherPdg==kGamma) { 
          isFromConversion=kTRUE; 
        }
      }
    } 
  } 

  return isFromConversion;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskFilteredTree::IsFromMaterial(Int_t label, AliStack *const stack) 
{
  Bool_t isFromMaterial = kFALSE;

  if(stack) {
    Int_t mcStackSize=stack->GetNtrack();
    if (label>=mcStackSize) return kFALSE;
    TParticle* particle = stack->Particle(label);
    if (!particle) return kFALSE;

    if(particle && particle->GetPDG() && particle->GetPDG()->Charge()!=0) 
    {
      Int_t mech = particle->GetUniqueID(); // production mechanism 
      Bool_t isPrim = stack->IsPhysicalPrimary(label);

      Int_t motherLabel = TMath::Abs(particle->GetMother(0));  
      if (motherLabel>=mcStackSize) return kFALSE;
      TParticle* mother = stack->Particle(motherLabel); 
      if(mother) {
        if(!isPrim && mech==13) { 
          isFromMaterial=kTRUE; 
        }
      }
    } 
  } 

  return isFromMaterial;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskFilteredTree::IsFromStrangeness(Int_t label, AliStack *const stack) 
{
  Bool_t isFromStrangeness = kFALSE;

  if(stack) {
    Int_t mcStackSize=stack->GetNtrack();
    if (label>=mcStackSize) return kFALSE;
    TParticle* particle = stack->Particle(label);
    if (!particle) return kFALSE;

    if(particle && particle->GetPDG() && particle->GetPDG()->Charge()!=0) 
    {
      Int_t mech = particle->GetUniqueID(); // production mechanism 
      Bool_t isPrim = stack->IsPhysicalPrimary(label);

      Int_t motherLabel = TMath::Abs(particle->GetMother(0));  
      if (motherLabel>=mcStackSize) return kFALSE;
      TParticle* mother = stack->Particle(motherLabel); 
      if(mother) {
        Int_t motherPdg = mother->GetPdgCode();

        // K+-, lambda, antilambda, K0s decays
        if(!isPrim && mech==4 && 
            (TMath::Abs(motherPdg)==kKPlus || TMath::Abs(motherPdg)==kLambda0 || motherPdg==kK0Short))
        {
          isFromStrangeness = kTRUE;
        } 
      }
    } 
  } 

  return isFromStrangeness;
}


//_____________________________________________________________________________
void AliAnalysisTaskFilteredTree::FinishTaskOutput() 
{
  //
  // Called one at the end 
  // locally on working node
  //
  Bool_t deleteTrees=kTRUE;
  if ((AliAnalysisManager::GetAnalysisManager()))
  {
    if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() == 
        AliAnalysisManager::kProofAnalysis)
      deleteTrees=kFALSE;
  }
  if (deleteTrees) delete fTreeSRedirector;
  fTreeSRedirector=NULL;
}

//_____________________________________________________________________________
void AliAnalysisTaskFilteredTree::Terminate(Option_t *) 
{
  //
  // Called one at the end 
  //
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskFilteredTree::GetMCTrueTrackMult(AliMCEvent *const mcEvent, AliFilteredTreeEventCuts *const evtCuts, AliFilteredTreeAcceptanceCuts *const accCuts)
{
  //
  // calculate mc event true track multiplicity
  //
  if(!mcEvent) return 0;

  AliStack* stack = 0;
  Int_t mult = 0;

  // MC particle stack
  stack = mcEvent->Stack();
  if (!stack) return 0;

  //
  //printf("minZv %f, maxZv %f \n", evtCuts->GetMinZv(), evtCuts->GetMaxZv());
  //

  Bool_t isEventOK = evtCuts->AcceptMCEvent(mcEvent);
  if(!isEventOK) return 0; 

  Int_t nPart  = stack->GetNtrack();
  for (Int_t iMc = 0; iMc < nPart; ++iMc) 
  {
    TParticle* particle = stack->Particle(iMc);
    if (!particle)
      continue;

    // only charged particles
    if(!particle->GetPDG()) continue;
    Double_t charge = particle->GetPDG()->Charge()/3.;
    if (TMath::Abs(charge) < 0.001)
      continue;

    // physical primary
    Bool_t prim = stack->IsPhysicalPrimary(iMc);
    if(!prim) continue;

    // checked accepted without pt cut
    //if(accCuts->AcceptTrack(particle)) 
    if( particle->Eta() > accCuts->GetMinEta() && particle->Eta() < accCuts->GetMaxEta() ) 
    {
      mult++;
    }
  }

  return mult;  
}

//_____________________________________________________________________________
void AliAnalysisTaskFilteredTree::FillHistograms(AliESDtrack* const ptrack, AliExternalTrackParam* const ptpcInnerC, Double_t centralityF, Double_t chi2TPCInnerC) 
{
  //
  // Fill pT relative resolution histograms for 
  // TPC only, TPC only constrained to vertex and TPC+ITS tracking
  //
  if(!ptrack) return;    
  if(!ptpcInnerC) return;    

  const AliExternalTrackParam * innerParam = (AliExternalTrackParam *) ptrack->GetInnerParam();
  if(!innerParam) return;

  Float_t dxy, dz;
  ptrack->GetImpactParameters(dxy,dz);

  // TPC+ITS primary tracks 
  if( abs(ptrack->Eta())<0.8 && 
      ptrack->GetTPCClusterInfo(3,1)>120 && 
      ptrack->IsOn(0x40) && 
      ptrack->GetTPCclusters(0)>0.0 &&  
      ptrack->GetTPCnclsS()/ptrack->GetTPCclusters(0)<0.4 && 
      //abs(innerParam->GetX())>0.0 && 
      //abs(innerParam->GetY()/innerParam->GetX())<0.14 && 
      //abs(innerParam->GetTgl())<0.85 && 
      ptrack->IsOn(0x0004) && 
      ptrack->GetNcls(0)>0 &&
      ptrack->GetITSchi2()>0 && 
      sqrt(ptrack->GetITSchi2()/ptrack->GetNcls(0))<6 &&
      sqrt(chi2TPCInnerC)<6 &&
      (ptrack->HasPointOnITSLayer(0) || ptrack->HasPointOnITSLayer(1)) &&
      abs(dz)<2.0 && 
      abs(dxy)<(0.018+0.035*abs(ptrack->GetSigned1Pt())) )
  {
    fPtResPhiPtTPCITS->Fill(ptrack->Pt(),ptrack->Phi(),1./abs(ptrack->GetSigned1Pt())*TMath::Sqrt(ptrack->GetSigma1Pt2()));
    fPtResEtaPtTPCITS->Fill(ptrack->Pt(),ptrack->Eta(),1./abs(ptrack->GetSigned1Pt())*TMath::Sqrt(ptrack->GetSigma1Pt2()));
    fPtResCentPtTPCITS->Fill(ptrack->Pt(),centralityF,1./abs(ptrack->GetSigned1Pt())*TMath::Sqrt(ptrack->GetSigma1Pt2()));
  }

  // TPC primary tracks 
  // and TPC constrained primary tracks 

  AliExternalTrackParam *ptpcInner  = (AliExternalTrackParam *) ptrack->GetTPCInnerParam(); 
  if(!ptpcInner) return;


  Float_t dxyTPC, dzTPC;
  ptrack->GetImpactParametersTPC(dxyTPC,dzTPC);

  if( abs(ptrack->Eta())<0.8 && 
      ptrack->GetTPCClusterInfo(3,1)>120 && 
      ptrack->IsOn(0x40)&& 
      ptrack->GetTPCclusters(0)>0.0 &&  
      ptrack->GetTPCnclsS()/ptrack->GetTPCclusters(0)<0.4 && 
      //abs(innerParam->GetX())>0.0 && 
      //abs(innerParam->GetY()/innerParam->GetX())<0.14 && 
      //abs(innerParam->GetTgl())<0.85 && 
      abs(dzTPC)<3.2 && 
      abs(dxyTPC)<2.4 )
  {
    // TPC only
    fPtResPhiPtTPC->Fill(ptpcInner->Pt(),ptpcInner->Phi(),1./abs(ptpcInner->GetSigned1Pt())*TMath::Sqrt(ptpcInner->GetSigma1Pt2()));
    fPtResEtaPtTPC->Fill(ptpcInner->Pt(),ptpcInner->Eta(),1./abs(ptpcInner->GetSigned1Pt())*TMath::Sqrt(ptpcInner->GetSigma1Pt2()));
    fPtResCentPtTPC->Fill(ptpcInner->Pt(),centralityF,1./abs(ptpcInner->GetSigned1Pt())*TMath::Sqrt(ptpcInner->GetSigma1Pt2()));

    // TPC constrained to vertex 
    fPtResPhiPtTPCc->Fill(ptpcInnerC->Pt(),ptpcInnerC->Phi(),1./abs(ptpcInnerC->GetSigned1Pt())*TMath::Sqrt(ptpcInnerC->GetSigma1Pt2()));
    fPtResEtaPtTPCc->Fill(ptpcInnerC->Pt(),ptpcInnerC->Eta(),1./abs(ptpcInnerC->GetSigned1Pt())*TMath::Sqrt(ptpcInnerC->GetSigma1Pt2()));
    fPtResCentPtTPCc->Fill(ptpcInnerC->Pt(),centralityF,1./abs(ptpcInnerC->GetSigned1Pt())*TMath::Sqrt(ptpcInnerC->GetSigma1Pt2()));
  }
}


void AliAnalysisTaskFilteredTree::ProcessITSTPCmatchOut(AliESDEvent *const esdEvent, AliESDfriend *const esdFriend){
  //
  // Process ITS standalone tracks find match with closest TPC(or combined tracks) tracks 
  // marian.ivanov@cern.ch
  // 0.) Init variables
  // 1.) GetTrack parameters at TPC inner wall
  // 2.) Match closest TPC  track  (STANDALONE/global) - chi2 match criteria
  //
  // Logic to be used in reco:
  // 1.) Find matching ITSalone->TPCalone
  // 2.) if (!TPCalone.FindClose(TPCother))  TPCalone.Addopt(ITSalone)
  // 3.) ff ((ITSalone.FindClose(Global)==0) CreateGlobaltrack
  const Double_t radiusMatch=84.;    // redius to propagate
  //
  const Double_t dFastPhiCut=0.2;        // 6 sigma (200 MeV) fast angular cut
  const Double_t dFastThetaCut=0.12;     // 6 sigma (200 MeV) fast angular cut
  const Double_t dFastPosPhiCut=0.06;    // 6 sigma (200 MeV) fast angular cut
  const Double_t dFastZCut=6;            // 6 sigma (200 MeV) fast  z difference cut
  const Double_t dFastPtCut=2.;          // 6 sigma (200 MeV) fast 1/pt cut 
  const Double_t chi2Cut=100;            // chi2 matching cut
  //
  if (!esdFriend) return;  // not ITS standalone track
  if (esdFriend->TestSkipBit()) return; // friends tracks  not stored
  Int_t ntracks=esdEvent->GetNumberOfTracks();
  Float_t bz = esdEvent->GetMagneticField();
  //
  // 0.) Get parameters in reference radius TPC Inner wall
  //
  //
  TMatrixD vecPosR0(ntracks,6);   // possition and  momentum estimate at reference radius  
  TMatrixD vecMomR0(ntracks,6);   //
  TMatrixD vecPosR1(ntracks,6);   // possition and  momentum estimate at reference radius TPC track  
  TMatrixD vecMomR1(ntracks,6);   //
  Double_t xyz[3], pxyz[3];      //
  for (Int_t iTrack=0; iTrack<ntracks; iTrack++){
    AliESDtrack *track = esdEvent->GetTrack(iTrack);   
    if(!track) continue;
    if (track->GetInnerParam()){
      const AliExternalTrackParam *trackTPC=track->GetInnerParam();
      trackTPC->GetXYZAt(radiusMatch,bz,xyz);
      trackTPC->GetPxPyPzAt(radiusMatch,bz,pxyz);
      for (Int_t i=0; i<3; i++){
	vecPosR1(iTrack,i)=xyz[i];
	vecMomR1(iTrack,i)=pxyz[i];
      }
      vecPosR1(iTrack,3)= TMath::ATan2(xyz[1],xyz[0]);    // phi pos angle
      vecMomR1(iTrack,3)= TMath::ATan2(pxyz[1],pxyz[0]);  // phi mom angle
      vecMomR1(iTrack,4)= trackTPC->GetSigned1Pt();;    
      vecMomR1(iTrack,5)= trackTPC->GetTgl();;    
    }
    AliESDfriendTrack* friendTrack=esdFriend->GetTrack(iTrack);
    if(!friendTrack) continue;
    if (friendTrack->GetITSOut()){
      const AliExternalTrackParam *trackITS=friendTrack->GetITSOut();
      trackITS->GetXYZAt(radiusMatch,bz,xyz);
      trackITS->GetPxPyPzAt(radiusMatch,bz,pxyz);
      for (Int_t i=0; i<3; i++){
	vecPosR0(iTrack,i)=xyz[i];
	vecMomR0(iTrack,i)=pxyz[i];
      }
      vecPosR0(iTrack,3)= TMath::ATan2(xyz[1],xyz[0]);
      vecMomR0(iTrack,3)= TMath::ATan2(pxyz[1],pxyz[0]);
      vecMomR0(iTrack,4)= trackITS->GetSigned1Pt();;    
      vecMomR0(iTrack,5)= trackITS->GetTgl();;    
    }
  }
  //  
  // 1.) Find closest matching tracks, between the ITS standalone track 
  // and  the all other tracks
  //  a.) caltegory  - All
  //  b.) category   - without ITS 
  //      
  //
  Int_t ntracksPropagated=0;
  AliExternalTrackParam extTrackDummy;
  AliESDtrack           esdTrackDummy; 
  AliExternalTrackParam itsAtTPC;
  AliExternalTrackParam itsAtITSTPC;
  for (Int_t iTrack0=0; iTrack0<ntracks; iTrack0++){
    AliESDtrack *track0 = esdEvent->GetTrack(iTrack0);   
    if(!track0) continue;
    if (track0->IsOn(AliVTrack::kTPCin)) continue;
    AliESDfriendTrack* friendTrack0=esdFriend->GetTrack(iTrack0); 
    if (!friendTrack0) continue;
    //if (!track0->IsOn(AliVTrack::kITSpureSA)) continue;
    //if (!friendTrack0->GetITSOut()) continue;  // is there flag for ITS standalone?
    ntracksPropagated++;
    // 
    // 2.) find clostest TPCtrack
    //     a.) all tracks
    Double_t minChi2All=10000000;
    Double_t minChi2TPC=10000000;
    Double_t minChi2TPCITS=10000000;
    Int_t indexAll=-1;
    Int_t indexTPC=-1;
    Int_t indexTPCITS=-1;
    Int_t ncandidates0=0; // n candidates - rough cut
    Int_t ncandidates1=0; // n candidates - rough + chi2 cut
    itsAtTPC=*(friendTrack0->GetITSOut());
    itsAtITSTPC=*(friendTrack0->GetITSOut());
    for (Int_t iTrack1=0; iTrack1<ntracks; iTrack1++){
      AliESDtrack *track1 = esdEvent->GetTrack(iTrack1);   
      if(!track1) continue;
      if (!track1->IsOn(AliVTrack::kTPCin)) continue;
      // fast checks
      //
      if (TMath::Abs(vecPosR1(iTrack1,2)-vecPosR0(iTrack0,2))>dFastZCut) continue;
      if (TMath::Abs(vecPosR1(iTrack1,3)-vecPosR0(iTrack0,3))>dFastPosPhiCut) continue;
      if (TMath::Abs(vecMomR1(iTrack1,3)-vecMomR0(iTrack0,3))>dFastPhiCut) continue;
      if (TMath::Abs(vecMomR1(iTrack1,5)-vecMomR0(iTrack0,5))>dFastThetaCut) continue;
      if (TMath::Abs(vecMomR1(iTrack1,4)-vecMomR0(iTrack0,4))>dFastPtCut) continue;
      ncandidates0++;
      //
      const AliExternalTrackParam * param1= track1->GetInnerParam();
      if (!friendTrack0->GetITSOut()) continue;
      AliExternalTrackParam outerITS = *(friendTrack0->GetITSOut());
      if (!outerITS.Rotate(param1->GetAlpha())) continue;
      if (!outerITS.PropagateTo(param1->GetX(),bz)) continue; // assume track close to the TPC inner wall
      Double_t chi2 =  outerITS.GetPredictedChi2(param1);
      if (chi2>chi2Cut) continue;
      ncandidates1++;
      if (chi2<minChi2All){
	minChi2All=chi2;
	indexAll=iTrack1;
      }
      if (chi2<minChi2TPC && track1->IsOn(AliVTrack::kITSin)==0){
	minChi2TPC=chi2;
	indexTPC=iTrack1;
	itsAtTPC=outerITS;
      }
      if (chi2<minChi2TPCITS && track1->IsOn(AliVTrack::kITSin)){
	minChi2TPCITS=chi2;
	indexTPCITS=iTrack1;
	itsAtITSTPC=outerITS;
      }
    }
    //
    AliESDtrack * trackAll= (indexAll>=0)? esdEvent->GetTrack(indexAll):&esdTrackDummy;
    AliESDtrack * trackTPC= (indexTPC>=0)? esdEvent->GetTrack(indexTPC):&esdTrackDummy;
    AliESDtrack * trackTPCITS= (indexTPCITS>=0)? esdEvent->GetTrack(indexTPCITS):&esdTrackDummy;
    (*fTreeSRedirector)<<"itsTPC"<<
      "indexAll="<<indexAll<<          // index of closest track (chi2)
      "indexTPC="<<indexTPC<<          // index of closest TPCalone tracks
      "indexTPCITS="<<indexTPCITS<<    // index of closest cobined tracks
      "ncandidates0="<<ncandidates0<<  // number of candidates
      "ncandidates1="<<ncandidates1<<
      //
      "chi2All="<<minChi2All<<         // chi2 of closest  tracks
      "chi2TPC="<<minChi2TPC<<         
      "chi2TPCITS="<<minChi2TPCITS<<
      //
      "track0.="<<track0<<             // ITS standalone tracks
      "trackAll.="<<trackAll<<         // Closets other track
      "trackTPC.="<<trackTPC<<         // Closest TPC only track
      "trackTPCITS.="<<trackTPCITS<<   // closest combined track
      //
      "itsAtTPC.="<<&itsAtTPC<<        // ITS track parameters at the TPC alone track  frame
      "itsAtITSTPC.="<<&itsAtITSTPC<<  // ITS track parameters at the TPC combeined track  frame
      "\n"; 
  }
}

void AliAnalysisTaskFilteredTree::ProcessTrackMatch(AliESDEvent *const /*esdEvent*/, AliESDfriend *const /*esdFriend*/){
/*
    Use TPC standalone, ITS standalone and combined tracks to categorize/resp. recover track information.

    Track categories:
       -TPC+ITS
       -TPC only 
       -ITS only
    Topology categories:
       -Nice isolated tracks with full information 
       -Overlapped tracks - Refit and sign them
       -Multiple found (check overlap factor) - Merge and sign
       -Charge particle (kink) decays - Change of direction - Sign them) 
    Info:
       -Array  of indexes of closest tracks in each individual category
       -Chi2  distance to the closest tracks at reference radius of TPCin
       -Overlap factors  - fraction of shared clusters or missing  region
       -Chi2 distance at DCA
    Information matrix:   
       -matrix closest tracks from each categories
       -characterization - chi2, index,chi2,  overlap ratio
    
    Decissions:
       0.) Kink decay or catastophic multiple scaterring 
           (combining all track categories)
              - small chi2 at the DCA
              - significantly large deflection angle
              - Combinatorial algorithm - to decrease CPU time restriction of investigation to tracks with small DCA at  refernce radius

       1.) if (TPConly && !(TPC+ITS) && ITSonly match ) {
             if (!kink) TPCOnly.addoptITS
             if (kink) TPConly sign
           }

       2.) Overlap tracks - Refit with doUnfold
*/
}

Int_t   AliAnalysisTaskFilteredTree::GetNearestTrack(const AliExternalTrackParam * trackMatch, Int_t indexSkip, AliESDEvent*event, Int_t trackType, Int_t paramType, AliExternalTrackParam & paramNearest){
  //
  // Find track with closest chi2 distance  (assume all track ae propagated to the DCA)
  //   trackType = 0 - find closets ITS standalone
  //               1 - find closest track with TPC
  //               2 - closest track with ITS and TPC
  //   paramType = 0 - global track
  //               1 - track at inner wall of TPC
  //
  //          
  if (trackMatch==NULL){
    ::Error("AliAnalysisTaskFilteredTree::GetNearestTrack","invalid track pointer");
    return -1;
  }
  Int_t ntracks=event->GetNumberOfTracks();
  const Double_t ktglCut=0.1;
  const Double_t kqptCut=0.4;
  const Double_t kAlphaCut=0.2;
  //
  Double_t chi2Min=100000;
  Int_t indexMin=-1;
  for (Int_t itrack=0; itrack<ntracks; itrack++){
    if (itrack==indexSkip) continue;
    AliESDtrack *ptrack=event->GetTrack(itrack);
    if (ptrack==NULL) continue;
    if (trackType==0 && (ptrack->IsOn(0x1)==kFALSE || ptrack->IsOn(0x10)==kTRUE))  continue;     // looks for track without TPC information
    if (trackType==1 && (ptrack->IsOn(0x10)==kFALSE))   continue;                                // looks for tracks with   TPC information
    if (trackType==2 && (ptrack->IsOn(0x1)==kFALSE || ptrack->IsOn(0x10)==kFALSE)) continue;      // looks for tracks with   TPC+ITS information
    
    if (ptrack->GetKinkIndex(0)<0) continue;              // skip kink daughters
    const AliExternalTrackParam * track=0;                // 
    if (paramType==0) track=ptrack;                       // Global track         
    if (paramType==1) track=ptrack->GetInnerParam();      // TPC only track at inner wall of TPC
    if (track==NULL) {
      continue;
    }
    // first rough cuts
    // fP3 cut
    if (TMath::Abs((track->GetTgl()-trackMatch->GetTgl()))>ktglCut) continue; 
    // fP4 cut 
    if (TMath::Abs((track->GetSigned1Pt()-trackMatch->GetSigned1Pt()))>kqptCut) continue; 
    // fAlpha cut
    //Double_t alphaDist=TMath::Abs((track->GetAlpha()-trackMatch->GetAlpha()));
    Double_t alphaDist=TMath::Abs(TMath::ATan2(track->Py(),track->Px())-TMath::ATan2(trackMatch->Py(),trackMatch->Py()));
    if (alphaDist>TMath::Pi()) alphaDist-=TMath::TwoPi();
    if (alphaDist>kAlphaCut) continue;
    // calculate and extract track with smallest chi2 distance
    AliExternalTrackParam param(*track);
    if (param.Rotate(trackMatch->GetAlpha())==kFALSE) continue;
    if (param.PropagateTo(trackMatch->GetX(),trackMatch->GetBz())==kFALSE) continue;
    Double_t chi2=trackMatch->GetPredictedChi2(&param);
    if (chi2<chi2Min){
      indexMin=itrack;
      chi2Min=chi2;
      paramNearest=param;
    }
  }
  return indexMin;

}


void  AliAnalysisTaskFilteredTree::SetDefaultAliasesV0(TTree *tree){
  //
  // SetAliases and Metadata for the V0 trees
  //
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  Double_t massLambda = pdg->GetParticle("Lambda0")->Mass();
  Double_t massK0 = pdg->GetParticle("K0")->Mass();
  Double_t massPion = pdg->GetParticle("pi+")->Mass();
  Double_t massProton = pdg->GetParticle("proton")->Mass();
  const Double_t livetimeK0=2.684341668932;  // livetime in cm (surpisely missing info in PDG - see root forum)
  const Double_t livetimeLambda=7.8875395;  // livetime in cm (missing info in PDG - see root forum)
  //
  tree->SetAlias("massPion",Form("(%f+0)",massPion));
  tree->SetAlias("massProton",Form("(%f+0)",massProton));
  tree->SetAlias("massK0",Form("(%f+0)",massK0));
  tree->SetAlias("massLambda",Form("(%f+0)",massLambda));
  //
  tree->SetAlias("livetimeK0",TString::Format("%f",livetimeK0));
  tree->SetAlias("livetimeLambda",TString::Format("%f",livetimeLambda));
  tree->SetAlias("livetimeLikeK0",TString::Format("exp(-v0.fRr/(sqrt((v0.P()/%f)^2+1)*%f))",massK0, livetimeK0)); // 
  tree->SetAlias("livetimeLikeLambda",TString::Format("exp(-v0.fRr/(sqrt((v0.P()/%f)^2+1)*%f))",massLambda,livetimeLambda)); // 
  tree->SetAlias("livetimeLikeGamma","v0.fRr/80"); // 
  tree->SetAlias("livetimeLikeBkg","v0.fRr/80"); //   
  // delta of mass
  tree->SetAlias("K0Delta","(v0.GetEffMass(2,2)-massK0)");
  tree->SetAlias("LDelta","(v0.GetEffMass(4,2)-massLambda)");
  tree->SetAlias("ALDelta","(v0.GetEffMass(2,4)-massLambda)");
  tree->SetAlias("EDelta","(v0.GetEffMass(0,0))");
  // pull of the mass
  tree->SetAlias("K0Pull","(v0.GetEffMass(2,2)-massK0)/v0.GetKFInfo(2,2,1)");
  tree->SetAlias("LPull","(v0.GetEffMass(4,2)-massLambda)/v0.GetKFInfo(4,2,1)");
  tree->SetAlias("ALPull","(v0.GetEffMass(2,4)-massLambda)/v0.GetKFInfo(2,4,1)");
  tree->SetAlias("EPull","EDelta/v0.GetKFInfo(0,0,1)");
  // effective pull of the mass - (empirical values from fits)
  tree->SetAlias("K0PullEff","K0Delta/sqrt((3.63321e-03)**2+(5.68795e-04*v0.Pt())**2)");
  tree->SetAlias("LPullEff","LDelta/sqrt((1.5e-03)**2+(1.8e-04*v0.Pt())**2)");
  tree->SetAlias("ALPullEff","ALDelta/sqrt((1.5e-03)**2+(1.8e-04*v0.Pt())**2)");
  tree->SetAlias("EPullEff","v0.GetEffMass(0,0)/sqrt((5e-03)**2+(1.e-04*v0.Pt())**2)");
  // 
  tree->SetAlias("dEdx0DProton","AliMathBase::BetheBlochAleph(track0.fIp.P()/massProton)");
  tree->SetAlias("dEdx1DProton","AliMathBase::BetheBlochAleph(track1.fIp.P()/massProton)");
  tree->SetAlias("dEdx0DPion","AliMathBase::BetheBlochAleph(track0.fIp.P()/massPion)");
  tree->SetAlias("dEdx1DPion","AliMathBase::BetheBlochAleph(track1.fIp.P()/massPion)");
  // V0 - cuts - for PID 
  tree->SetAlias("cutDist","sqrt((track0.fIp.fP[0]-track1.fIp.fP[0])**2+(track0.fIp.fP[1]-track1.fIp.fP[1])**2)>3");
  tree->SetAlias("cutLong","track0.GetTPCClusterInfo(3,1,0)-5*abs(track0.fP[4])>130&&track1.GetTPCClusterInfo(3,1,0)>130-5*abs(track0.fP[4])");
  tree->SetAlias("cutPID","track0.fTPCsignal>0&&track1.fTPCsignal>0");
  tree->SetAlias("cutResol","sqrt(track0.fC[14]/track0.fP[4])<0.15&&sqrt(track1.fC[14]/track1.fP[4])<0.15");
  tree->SetAlias("cutV0","cutPID&&cutLong&&cutResol"); 
  //  
  tree->SetAlias("K0PullBkg","min(min(abs(LPull),abs(ALPull)),abs(EPull))+0");
  tree->SetAlias("LambdaPullBkg","min(min(abs(K0Pull),abs(ALPull)),abs(EPull)+0)");
  tree->SetAlias("ALambdaPullBkg","min(min(abs(K0Pull),abs(LPull)),abs(EPull)+0)");
  tree->SetAlias("EPullBkg","min(min(abs(K0Pull),abs(LPull)),abs(ALPull)+0)");
  //
  tree->SetAlias("K0Selected",      "abs(K0Pull)<3. &&abs(K0PullEff)<3.  && abs(LPull)>3  && abs(ALPull)>3 &&v0.PtArmV0()>0.11"); 
  tree->SetAlias("LambdaSelected",  "abs(LPull)<3.  &&abs(LPullEff)<3.   && abs(K0Pull)>3 && abs(EPull)>3  && abs(EDelta)>0.05");  
  tree->SetAlias("ALambdaSelected", "abs(ALPull)<3. &&abs(ALPullEff)<3   && abs(K0Pull)>3 && abs(EPull)>3  &&abs(EDelta)>0.05");
  tree->SetAlias("GammaSelected", "abs(EPull)<3     && abs(K0Pull)>3 && abs(LPull)>3 && abs(ALPull)>3");

  tree->SetAlias("K0Like0","exp(-K0Pull^2)*livetimeLikeK0");
  tree->SetAlias("LLike0","exp(-LPull^2)*livetimeLikeLambda");
  tree->SetAlias("ALLike0","exp(-ALPull^2)*livetimeLikeLambda");
  tree->SetAlias("ELike0","exp(-abs(EPull)*0.2)*livetimeLikeGamma");
  tree->SetAlias("BkgLike","0.000005*ntracks");  // backround coeefecint  to be fitted - depends on other cuts 
  //
  tree->SetAlias("V0Like","exp(-acos(v0.fPointAngle)*v0.fRr/0.36)*exp(-sqrt(kf.GetChi2())/0.5)");
  tree->SetAlias("ELike","(V0Like*ELike0)/(V0Like*(K0Like0+LLike0+ALLike0+ELike0)+BkgLike)");
  tree->SetAlias("K0Like","K0Like0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike)");
  tree->SetAlias("LLike","LLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike)");
  tree->SetAlias("ALLike","ALLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike)");
  //
  tree->SetAlias("K0PIDPull","(abs(track0.fTPCsignal/dEdx0DPion-50)+abs(track1.fTPCsignal/dEdx1DPion-50))/5.");
  tree->SetAlias("mpt","1/v0.Pt()");                 // 
  tree->SetAlias("tglV0","v0.Pz()/v0.Pt()");                 // 
  tree->SetAlias("alphaV0","atan2(v0.Py(),v0.Px()+0)");
  tree->SetAlias("dalphaV0","alphaV0-((int(36+9*(alphaV0/pi))-36)*pi/9.)");

}

void  AliAnalysisTaskFilteredTree::SetDefaultAliasesHighPt(TTree *tree){
  //
  // set shortcut aliases for some variables
  //
  tree->SetAlias("phiInner","atan2(esdTrack.fIp.Py(),esdTrack.fIp.Px()+0)");
  tree->SetAlias("secInner","9*(atan2(esdTrack.fIp.Py(),esdTrack.fIp.Px()+0)/pi)+18*(esdTrack.fIp.Py()<0)");
  tree->SetAlias("tgl","esdTrack.fP[3]");
  tree->SetAlias("alphaV","esdTrack.fAlpha");
  tree->SetAlias("qPt","esdTrack.fP[4]");
  tree->SetAlias("dalphaQ","sign(esdTrack.fP[4])*(esdTrack.fIp.fP[0]/esdTrack.fIp.fX)");
  TStatToolkit::AddMetadata(tree,"phiInner.Title","#phi_{TPCin}");
  TStatToolkit::AddMetadata(tree,"secInner.Title","sector_{TPCin}");
  TStatToolkit::AddMetadata(tree,"tgl.Title","#it{p_{z}}/#it{p}_{T}");
  TStatToolkit::AddMetadata(tree,"alphaV.Title","#phi_{vertex}");
  TStatToolkit::AddMetadata(tree,"qPt.Title","q/#it{p}_{T}");
  TStatToolkit::AddMetadata(tree,"phiInner.AxisTitle","#it{#phi}_{TPCin}");
  TStatToolkit::AddMetadata(tree,"secInner.AxisTitle","sector_{TPCin}");
  TStatToolkit::AddMetadata(tree,"tgl.AxisTitle","#it{p}_{z}/#it{p}_{t}");
  TStatToolkit::AddMetadata(tree,"alphaV.AxisTitle","#it{#phi}_{vertex}");
  TStatToolkit::AddMetadata(tree,"qPt.AxisTitle","q/#it{p}_{T} (1/GeV)");
  
  //
  tree->SetAlias("normChi2ITS","sqrt(esdTrack.fITSchi2/esdTrack.fITSncls)");
  tree->SetAlias("normChi2TPC","esdTrack.fTPCchi2/esdTrack.fTPCncls");
  tree->SetAlias("normChi2TRD","esdTrack.fTRDchi2/esdTrack.fTRDncls");
  tree->SetAlias("normDCAR","esdTrack.fdTPC/sqrt(1+esdTrack.fP[4]**2)");
  tree->SetAlias("normDCAZ","esdTrack.fzTPC/sqrt(1+esdTrack.fP[4]**2)");
  TStatToolkit::AddMetadata(tree,"normChi2ITS.Title","#sqrt{#chi2_{ITS}/N_{clITS}}");
  TStatToolkit::AddMetadata(tree,"normChi2TPC.Title","#chi2_{TPC}/N_{clTPC}");
  TStatToolkit::AddMetadata(tree,"normChi2ITS.AxisTitle","#sqrt{#chi2_{ITS}/N_{clITS}}");
  TStatToolkit::AddMetadata(tree,"normChi2TPC.AxisTitle","#chi2_{TPC}/N_{clTPC}");
  //
  tree->SetAlias("TPCASide","esdTrack.fIp.fP[1]>0");
  tree->SetAlias("TPCCSide","esdTrack.fIp.fP[1]<0");
  tree->SetAlias("TPCCross","esdTrack.fIp.fP[1]*esdTrack.fIp.fP[3]<0");
  tree->SetAlias("ITSOn","((esdTrack.fFlags&0x1)>0)");
  tree->SetAlias("TPCOn","((esdTrack.fFlags&0x10)>0)");
  tree->SetAlias("ITSRefit","((esdTrack.fFlags&0x4)>0)");
  tree->SetAlias("TPCRefit","((esdTrack.fFlags&0x40)>0)");
  tree->SetAlias("TOFOn","((esdTrack.fFlags&0x2000)>0)");
  tree->SetAlias("TRDOn","((esdTrack.fFlags&0x400)>0)");
  tree->SetAlias("ITSOn0","esdTrack.fITSncls>4&&esdTrack.HasPointOnITSLayer(0)&&esdTrack.HasPointOnITSLayer(1)");
  tree->SetAlias("ITSOn01","esdTrack.fITSncls>3&&(esdTrack.HasPointOnITSLayer(0)||esdTrack.HasPointOnITSLayer(1))");
  tree->SetAlias("nclCut","(esdTrack.GetTPCClusterInfo(3,1)+esdTrack.fTRDncls)>140-5*(abs(esdTrack.fP[4]))");
  tree->SetAlias("IsPrim4","sqrt((esdTrack.fD**2)/esdTrack.fCdd+(esdTrack.fZ**2)/esdTrack.fCzz)<4");
  tree->SetAlias("IsPrim4TPC","sqrt((esdTrack.fdTPC**2)/esdTrack.fCddTPC+(esdTrack.fzTPC**2)/esdTrack.fCzzTPC)<4");
  

  const char * chName[5]={"#it{r#phi}","#it{z}","sin(#phi)","tan(#it{#theta})", "q/#it{p}_{t}"};
  const char * chUnit[5]={"cm","cm","","", "(1/GeV)"};
  const char * refBranch=(tree->GetBranch("extInnerParamV."))? "extInnerParamV":"esdTrack.fTPCInner";

  for (Int_t iPar=0; iPar<5; iPar++){
    tree->SetAlias(TString::Format("covarP%dITS",iPar).Data(),TString::Format("sqrt(esdTrack.fC[%d]+0)",AliExternalTrackParam::GetIndex(iPar,iPar)).Data());    
    tree->SetAlias(TString::Format("covarP%d",iPar).Data(),TString::Format("sqrt(%s.fC[%d]+0)",refBranch,AliExternalTrackParam::GetIndex(iPar,iPar)).Data());
    tree->SetAlias(TString::Format("covarPC%d",iPar).Data(),TString::Format("sqrt(extInnerParamC.fC[%d]+0)",AliExternalTrackParam::GetIndex(iPar,iPar)).Data());
    tree->SetAlias(TString::Format("covarP%dNorm",iPar).Data(),TString::Format("sqrt(%s.fC[%d]+0)/sqrt(1+(1*esdTrack.fP[4])**2)/sqrt(1+(1*esdTrack.fP[4])**2)",refBranch,AliExternalTrackParam::GetIndex(iPar,iPar)).Data());
    tree->SetAlias(TString::Format("covarPC%dNorm",iPar).Data(),TString::Format("sqrt(extInnerParamC.fC[%d]+0)/sqrt(1+(5*esdTrack.fP[4])**2)",AliExternalTrackParam::GetIndex(iPar,iPar)).Data());

    tree->SetAlias(TString::Format("deltaP%d",iPar).Data(),TString::Format("(%s.fP[%d]-esdTrack.fCp.fP[%d])",refBranch,iPar,iPar).Data());
    tree->SetAlias(TString::Format("deltaPC%d",iPar).Data(),TString::Format("(extInnerParamC.fP[%d]-esdTrack.fCp.fP[%d])",iPar,iPar).Data());
    // rough  normalization of the residual sigma ~ sqrt(1+*k/pt)^2) 
    // histogramming pt normalized deltas enable us to use wider bins in q/pt and also less bins in deltas 
    tree->SetAlias(TString::Format("deltaP%dNorm",iPar).Data(),TString::Format("(%s.fP[%d]-esdTrack.fCp.fP[%d])/sqrt(1+(1*esdTrack.fP[4])**2)",refBranch,iPar,iPar).Data());
    tree->SetAlias(TString::Format("deltaPC%dNorm",iPar).Data(),TString::Format("(extInnerParamC.fP[%d]-esdTrack.fCp.fP[%d])/sqrt(1.+(5.*esdTrack.fP[4])**2)",iPar,iPar).Data());
    //
    tree->SetAlias(TString::Format("pullP%d",iPar).Data(),
		   TString::Format("(%s.fP[%d]-esdTrack.fCp.fP[%d])/sqrt(%s.fC[%d]+esdTrack.fCp.fC[%d])",refBranch,
				   iPar,iPar,refBranch,AliExternalTrackParam::GetIndex(iPar,iPar),AliExternalTrackParam::GetIndex(iPar,iPar)).Data());
    tree->SetAlias(TString::Format("pullPC%d",iPar).Data(),
		   TString::Format("(extInnerParamC.fP[%d]-esdTrack.fCp.fP[%d])/sqrt(extInnerParamC.fC[%d]+esdTrack.fCp.fC[%d])",
				   iPar,iPar,AliExternalTrackParam::GetIndex(iPar,iPar),AliExternalTrackParam::GetIndex(iPar,iPar)).Data());
    TStatToolkit::AddMetadata(tree,TString::Format("deltaP%d.AxisTitle",iPar).Data(),TString::Format("%s (%s)",chName[iPar], chUnit[iPar]).Data());
    TStatToolkit::AddMetadata(tree,TString::Format("deltaPC%d.AxisTitle",iPar).Data(),TString::Format("%s (%s)",chName[iPar], chUnit[iPar]).Data());
    TStatToolkit::AddMetadata(tree,TString::Format("pullP%d.AxisTitle",iPar).Data(),TString::Format("pull %s (unit)",chName[iPar]).Data());
    TStatToolkit::AddMetadata(tree,TString::Format("pullPC%d.AxisTitle",iPar).Data(),TString::Format("pull %s (unit)",chName[iPar]).Data());
    //
    TStatToolkit::AddMetadata(tree,TString::Format("deltaP%d.Title",iPar).Data(),TString::Format("%s",chName[iPar], chUnit[iPar]).Data());
    TStatToolkit::AddMetadata(tree,TString::Format("deltaPC%d.Title",iPar).Data(),TString::Format("%s",chName[iPar], chUnit[iPar]).Data());
    TStatToolkit::AddMetadata(tree,TString::Format("pullP%d.Title",iPar).Data(),TString::Format("pull %s",chName[iPar]).Data());
    TStatToolkit::AddMetadata(tree,TString::Format("pullPC%d.Title",iPar).Data(),TString::Format("pull %s",chName[iPar]).Data());
  }
  TStatToolkit::AddMetadata(tree, "mult.Title","N_{prim}");
  TStatToolkit::AddMetadata(tree, "mult.AxisTitle","N_{prim}");
  TStatToolkit::AddMetadata(tree, "ntracks.Title","N_{tr}");
  TStatToolkit::AddMetadata(tree, "ntracks.AxisTitle","N_{tr} (prim+sec+pile-up)");
} 
