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

#include "iostream"

#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include "TChain.h"
#include "TTreeStream.h"
#include "TTree.h"
#include "TH1F.h"
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
  , fProcessAll(kFALSE)
  , fProcessCosmics(kFALSE)
  , fHighPtTree(0)
  , fV0Tree(0)
  , fdEdxTree(0)
  , fLaserTree(0)
  , fMCEffTree(0)
  , fCosmicPairsTree(0)
{
  // Constructor

  // Define input and output slots here
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
  DefineOutput(5, TTree::Class());
  DefineOutput(6, TTree::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskFilteredTree::~AliAnalysisTaskFilteredTree()
{
  //the output trees not to be deleted in case of proof analysis
  Bool_t deleteTrees=kTRUE;
  //if ((AliAnalysisManager::GetAnalysisManager()))
  //{
  //  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() == 
  //           AliAnalysisManager::kProofAnalysis)
  //    deleteTrees=kFALSE;
  //}
  //if (deleteTrees) delete fTreeSRedirector;

  delete fFilteredTreeEventCuts;
  delete fFilteredTreeAcceptanceCuts;
  delete fFilteredTreeRecAcceptanceCuts;
  delete fEsdTrackCuts;
}

//____________________________________________________________________________
Bool_t AliAnalysisTaskFilteredTree::Notify()
{
  static Int_t count = 0;
  count++;
  TTree *chain = (TChain*)GetInputData(0);
  if(chain)
  {
    Printf("Processing %d. file: %s", count, chain->GetCurrentFile()->GetName());
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

  PostData(1,fV0Tree);
  PostData(2,fHighPtTree);
  PostData(3,fdEdxTree);
  PostData(4,fLaserTree);
  PostData(5,fMCEffTree);
  PostData(6,fCosmicPairsTree);
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

  //// MC event
  //if(fUseMCInfo) {
  //  fMC = MCEvent();
  //  if (!fMC) {
  //    Printf("ERROR: MC event not available");
  //    return;
  //  }
  //}
  
  //if MC info available - use it.
  fMC = MCEvent();

  if(fUseESDfriends) {
    fESDfriend = static_cast<AliESDfriend*>(fESD->FindListObject("AliESDfriend"));
      if(!fESDfriend) {
        Printf("ERROR: ESD friends not available");
    }
  }

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

  if (fProcessCosmics) { ProcessCosmics(fESD); }
  if(IsUseMCInfo()) { ProcessMCEff(fESD,fMC,fESDfriend);}
}

//_____________________________________________________________________________
void AliAnalysisTaskFilteredTree::ProcessCosmics(AliESDEvent *const event)
{
  //
  // Select real events with high-pT tracks 
  //
  if(!event) {
    AliDebug(AliLog::kError, "event not available");
    return;
  }

  // 
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!inputHandler)
  {
    Printf("ERROR: Could not receive input handler");
    return;
  }
   
  // get file name
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { 
    Printf("ERROR: Could not receive input chain");
    return;
  }
  TObjString fileName(chain->GetCurrentFile()->GetName());


    // check for cosmic pairs
    //
    // find cosmic pairs trigger by random trigger
    //
    //
    AliESDVertex *vertexSPD =  (AliESDVertex *)event->GetPrimaryVertexSPD();
    AliESDVertex *vertexTPC =  (AliESDVertex *)event->GetPrimaryVertexTPC(); 
    const Double_t kMinPt=0.8;
    const Double_t kMinPtMax=0.8;
    const Double_t kMinNcl=50;
    const Double_t kMaxDelta[5]={2,600,0.02,0.02,0.1};
    Int_t ntracks=event->GetNumberOfTracks(); 
    //  Float_t dcaTPC[2]={0,0};
    // Float_t covTPC[3]={0,0,0};

    UInt_t specie = event->GetEventSpecie();  // skip laser events
    if (specie==AliRecoParam::kCalib) return;
  


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
        //Bool_t hasFriend = kFALSE;
        //Bool_t hasITS=(track0->GetNcls(0)+track1->GetNcls(0)>4);
        //printf("DUMPHPTCosmic:%s|%f|%d|%d|%d\n",filename.Data(),(TMath::Min(track0->Pt(),track1->Pt())), eventNumber,hasFriend,hasITS);
        //      const AliExternalTrackParam * trackIn1 = track1->GetInnerParam();      
        //
        //               
        Int_t ntracksSPD = vertexSPD->GetNContributors();
        Int_t ntracksTPC = vertexTPC->GetNContributors();        
        Int_t runNumber     = event->GetRunNumber();        
        Int_t timeStamp    = event->GetTimeStamp();
        ULong64_t triggerMask = event->GetTriggerMask();
        Float_t magField    = event->GetMagneticField();
        TObjString triggerClass = event->GetFiredTriggerClasses().Data();
        
       //
      // Dump to the tree 
      // vertex
      // TPC-ITS tracks
      //

      //fCosmicPairsTree->Branch("fileName",&fileName,32000,0);
      //fCosmicPairsTree->Branch("runNumber",&runNumber,"runNumber/I");
      //fCosmicPairsTree->Branch("timeStamp",&timeStamp,"timeStamp/I");
      //fCosmicPairsTree->Branch("eventNumber",&eventNumber,"eventNumber/I");
      //fCosmicPairsTree->Branch("triggerMask",&triggerMask,32000,0);
      //fCosmicPairsTree->Branch("triggerClass",&triggerClass,32000,0);
      //fCosmicPairsTree->Branch("magField",&magField,"magField/F");
      //fCosmicPairsTree->Branch("ntracksSPD",&ntracksSPD,"ntracksSPD/I");
      //fCosmicPairsTree->Branch("ntracksTPC",&ntracksTPC,"ntracksTPC/I");
      //fCosmicPairsTree->Branch("vertexSPD",vertexSPD,32000,0);
      //fCosmicPairsTree->Branch("vertexTPC",vertexTPC,32000,0);
      //fCosmicPairsTree->Branch("track0",track0,32000,0);
      //fCosmicPairsTree->Branch("track1",track1,32000,0);
      //
      //fCosmicPairsTree->Fill();

      if(!fTreeSRedirector) return;
	  (*fTreeSRedirector)<<"CosmicPairs"<<
	    "fileName.="<<&fileName<<         // file name
	    "runNumber="<<runNumber<<              //  run number	    
	    "evtTimeStamp="<<timeStamp<<            //  time stamp of event
            "evtNumberInFile="<<eventNumber<<          //  event number	    
	    "trigger="<<triggerMask<<      //  trigger
 	    "triggerClass="<<&triggerClass<<      //  trigger
	    "Bz="<<magField<<             //  magnetic field
	    //
	    "multSPD="<<ntracksSPD<<
	    "multTPC="<<ntracksTPC<<
	    "vertSPD.="<<vertexSPD<<         //primary vertex -SPD
	    "vertTPC.="<<vertexTPC<<         //primary vertex -TPC
	    "t0.="<<track0<<              //track0
	    "t1.="<<track1<<              //track1
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
  if(!esdEvent) {
    AliDebug(AliLog::kError, "esdEvent not available");
    return;
  }

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
  if (!inputHandler)
  {
    Printf("ERROR: Could not receive input handler");
    return;
  }
   
  // get file name
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { 
    Printf("ERROR: Could not receive input chain");
    return;
  }
  TObjString fileName(chain->GetCurrentFile()->GetName());

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
  if(IsUseMCInfo())
  {
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

  if(!vtxESD) return;

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD); 
  //printf("isEventOK %d, isEventTriggered %d, status %d, vz %f \n",isEventOK, isEventTriggered, vtxESD->GetStatus(), vtxESD->GetZv());
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
    Double_t vert[3] = {0}; 
    vert[0] = vtxESD->GetXv();
    vert[1] = vtxESD->GetYv();
    vert[2] = vtxESD->GetZv();
    Int_t mult = vtxESD->GetNContributors();
    Float_t bz = esdEvent->GetMagneticField();
    Int_t runNumber = esdEvent->GetRunNumber();
    Int_t evtTimeStamp = esdEvent->GetTimeStamp();
    Int_t evtNumberInFile = esdEvent->GetEventNumberInFile();

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
      if(TMath::Exp(2*scalempt)<downscaleF) continue;
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

      //fHighPtTree->Branch("fileName",&fileName,32000,0);
      //fHighPtTree->Branch("runNumber",&runNumber,"runNumber/I");
      //fHighPtTree->Branch("evtTimeStamp",&evtTimeStamp,"evtTimeStamp/I");
      //fHighPtTree->Branch("evtNumberInFile",&evtNumberInFile,"evtNumberInFile/I");
      //fHighPtTree->Branch("triggerClass",&triggerClass,32000,0);
      //fHighPtTree->Branch("Bz",&bz,"Bz/F");
      //fHighPtTree->Branch("vtxESD",vtxESD,32000,0);
      //fHighPtTree->Branch("IRtot",&ir1,"IRtot/I");
      //fHighPtTree->Branch("IRint2",&ir2,"IRint2/I");
      //fHighPtTree->Branch("mult",&mult,"mult/I");
      //fHighPtTree->Branch("esdTrack",track,32000,0);
      //fHighPtTree->Branch("centralityF",&centralityF,"centralityF/F");

      //fHighPtTree->Fill();

      //Double_t vtxX=vtxESD->GetX();
      //Double_t vtxY=vtxESD->GetY();
      //Double_t vtxZ=vtxESD->GetZ();
      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"highPt"<<
        "fileName.="<<&fileName<<            
        "runNumber="<<runNumber<<
        "evtTimeStamp="<<evtTimeStamp<<
        "evtNumberInFile="<<evtNumberInFile<<
	"triggerClass="<<&triggerClass<<      //  trigger
        "Bz="<<bz<<                           //  magnetic field
        "vtxESD.="<<vtxESD<<
	"ntracksESD="<<ntracks<<              // number of tracks in the ESD
      //  "vtxESDx="<<vtxX<<
      //  "vtxESDy="<<vtxY<<
      //  "vtxESDz="<<vtxZ<<
	"IRtot="<<ir1<<                      // interaction record history info
	"IRint2="<<ir2<<
        "mult="<<mult<<                      // multiplicity of tracks pointing to the primary vertex
        "esdTrack.="<<track<<
        "centralityF="<<centralityF<<       
        "\n";
    }
  }
  
}


//_____________________________________________________________________________
void AliAnalysisTaskFilteredTree::ProcessLaser(AliESDEvent *const esdEvent, AliMCEvent * const /*mcEvent*/, AliESDfriend *const /*esdFriend*/)
{
  //
  // Process laser events
  //
  if(!esdEvent) {
    AliDebug(AliLog::kError, "esdEvent not available");
    return;
  }

  // get file name
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { 
    Printf("ERROR: Could not receive input chain");
    return;
  }
  TObjString fileName(chain->GetCurrentFile()->GetName());

  // laser events 
  const AliESDHeader* esdHeader = esdEvent->GetHeader();
  if(esdHeader && esdHeader->GetEventSpecie()==AliRecoParam::kCalib) 
  {
    Int_t countLaserTracks = 0;
    for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++)
    {
      AliESDtrack *track = esdEvent->GetTrack(iTrack);
      if(!track) continue;

      if(track->GetTPCInnerParam()) countLaserTracks++;
    }
       
    if(countLaserTracks > 100) 
    {      
      Int_t runNumber = esdEvent->GetRunNumber();
      Int_t evtTimeStamp = esdEvent->GetTimeStamp();
      Int_t evtNumberInFile = esdEvent->GetEventNumberInFile();
      Float_t bz = esdEvent->GetMagneticField();
      TObjString triggerClass = esdEvent->GetFiredTriggerClasses().Data();

      //fLaserTree->Branch("fileName",&fileName,32000,0);
      //fLaserTree->Branch("runNumber",&runNumber,"runNumber/I");
      //fLaserTree->Branch("evtTimeStamp",&evtTimeStamp,"evtTimeStamp/I");
      //fLaserTree->Branch("evtNumberInFile",&evtNumberInFile,"evtNumberInFile/I");
      //fLaserTree->Branch("triggerClass",&triggerClass,32000,0);
      //fLaserTree->Branch("Bz",&bz,"Bz/F");
      //fLaserTree->Branch("multTPCtracks",&countLaserTracks,"multTPCtracks/I");

      //fLaserTree->Fill();

      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"Laser"<<
        "fileName.="<<&fileName<<
        "runNumber="<<runNumber<<
        "evtTimeStamp="<<evtTimeStamp<<
        "evtNumberInFile="<<evtNumberInFile<<
	"triggerClass="<<&triggerClass<<      //  trigger
        "Bz="<<bz<<
        "multTPCtracks="<<countLaserTracks<<
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
  if (!inputHandler)
  {
    Printf("ERROR: Could not receive input handler");
    return;
  }
   
  // get file name
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { 
    Printf("ERROR: Could not receive input chain");
    return;
  }
  TObjString fileName(chain->GetCurrentFile()->GetName());

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
  if(IsUseMCInfo())
  {
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

  if(!vtxESD) return;

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD); 
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK, isEventTriggered);
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());


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
    //
    Int_t ir1 = 0;
    Int_t ir2 = 0;

    //
    Double_t vert[3] = {0}; 
    vert[0] = vtxESD->GetXv();
    vert[1] = vtxESD->GetYv();
    vert[2] = vtxESD->GetZv();
    Int_t mult = vtxESD->GetNContributors();
    Float_t bz = esdEvent->GetMagneticField();
    Int_t runNumber = esdEvent->GetRunNumber();
    Int_t evtTimeStamp = esdEvent->GetTimeStamp();
    Int_t evtNumberInFile = esdEvent->GetEventNumberInFile();

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
      if(TMath::Exp(2*scalempt)<downscaleF) continue;
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
      Bool_t isOKtrackInnerC = kFALSE;
      AliExternalTrackParam * trackInnerC =  new AliExternalTrackParam(*(track->GetInnerParam()));
      if (trackInnerC) {
        isOKtrackInnerC = ConstrainTrackInner(trackInnerC,vtxESD,track->GetMass(),b);
        isOKtrackInnerC = trackInnerC->Rotate(track->GetAlpha());
        isOKtrackInnerC = trackInnerC->PropagateTo(track->GetX(),esdEvent->GetMagneticField());
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

      if(esdFriend && esdFriend->TestSkipBit()==kFALSE) 
      {
        // propagate ITSout to TPC inner wall
        AliESDfriendTrack *friendTrack = esdFriend->GetTrack(iTrack);

        if(friendTrack) 
	{
          outerITSc = new AliExternalTrackParam(*friendTrack->GetITSOut());
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
      AliTrackReference *refITS = NULL;

      Bool_t isOKtrackInnerC3 = kFALSE;
      AliExternalTrackParam *trackInnerC3 = new AliExternalTrackParam(*(track->GetInnerParam()));

      if(IsUseMCInfo()) 
      {
        if(!stack) return;

        //
        // global track
	//
        Int_t label = TMath::Abs(track->GetLabel()); 
        particle = stack->Particle(label);
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
        particleTPC = stack->Particle(labelTPC);
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
        Int_t status = mcEvent->GetParticleAndTR(track->GetTPCLabel(), part, trefs);

	if(status>0 && part && trefs && part->GetPDG() && part->GetPDG()->Charge()!=0.) 
	{
	  Int_t nTrackRef = trefs->GetEntries();
	  //printf("nTrackRef %d \n",nTrackRef);

          Int_t countITS = 0;
	  for (Int_t iref = 0; iref < nTrackRef; iref++) 
          {
            AliTrackReference *ref = (AliTrackReference *)trefs->At(iref);

             // Ref. in the middle ITS 
            if(ref && ref->DetectorId()==AliTrackReference::kITS)
            {
	      if(!refITS && countITS==2) {
	        refITS = ref;
	        //printf("refITS %p \n",refITS);
	      }
	      countITS++;
            }

            // TPC
            if(ref && ref->DetectorId()==AliTrackReference::kTPC)
            {
	      if(!refTPCIn) {
	        refTPCIn = ref;
	        //printf("refTPCIn %p \n",refTPCIn);
	        //break;
	      }
            }
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
        particleITS = stack->Particle(labelITS);
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

      //
      Bool_t dumpToTree=kFALSE;
      
      if(isOKtpcInnerC  && isOKtrackInnerC) dumpToTree = kTRUE;
      if(fUseESDfriends && isOKtrackInnerC2 && isOKouterITSc) dumpToTree = kTRUE;
      if(fMC && isOKtrackInnerC3) dumpToTree = kTRUE;
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
      /////////////////
      //book keeping of created dummy objects (to avoid NULL in trees)
      Bool_t newvtxESD=kFALSE;
      Bool_t newtrack=kFALSE;
      Bool_t newtpcInnerC=kFALSE;
      Bool_t newtrackInnerC=kFALSE;
      Bool_t newtrackInnerC2=kFALSE;
      Bool_t newouterITSc=kFALSE;
      Bool_t newtrackInnerC3=kFALSE;
      Bool_t newrefTPCIn=kFALSE;
      Bool_t newrefITS=kFALSE;
      Bool_t newparticle=kFALSE;
      Bool_t newparticleMother=kFALSE;
      Bool_t newparticleTPC=kFALSE;
      Bool_t newparticleMotherTPC=kFALSE;
      Bool_t newparticleITS=kFALSE;
      Bool_t newparticleMotherITS=kFALSE;
      
      //check that the vertex is there and that it is OK, 
      //i.e. no null member arrays, otherwise a problem with merging
      //later on.
      //this is a very ugly hack!
      if (!vtxESD)
      {
        AliInfo("fixing the ESD vertex for streaming");
        vtxESD=new AliESDVertex(); 
        //vtxESD->SetNContributors(1);
        //UShort_t pindices[1]; pindices[0]=0;
        //vtxESD->SetIndices(1,pindices);
        //vtxESD->SetNContributors(0);
        newvtxESD=kTRUE;
      }
      //
      if (!track) {track=new AliESDtrack();newtrack=kTRUE;}
      if (!tpcInnerC) {tpcInnerC=new AliExternalTrackParam();newtpcInnerC=kTRUE;}
      if (!trackInnerC) {trackInnerC=new AliExternalTrackParam();newtrackInnerC=kTRUE;}
      if (!trackInnerC2) {trackInnerC2=new AliExternalTrackParam();newtrackInnerC2=kTRUE;}
      if (!outerITSc) {outerITSc=new AliExternalTrackParam();newouterITSc=kTRUE;}
      if (!trackInnerC3) {trackInnerC3=new AliExternalTrackParam();newtrackInnerC3=kTRUE;}
      if (fMC)
      {
        if (!refTPCIn) {refTPCIn=new AliTrackReference(); newrefTPCIn=kTRUE;}
        if (!refITS) {refITS=new AliTrackReference();newrefITS=kTRUE;}
        if (!particle) {particle=new TParticle(); newparticle=kTRUE;}
        if (!particleMother) {particleMother=new TParticle();newparticleMother=kTRUE;}
        if (!particleTPC) {particleTPC=new TParticle();newparticleTPC=kTRUE;}
        if (!particleMotherTPC) {particleMotherTPC=new TParticle();newparticleMotherTPC=kTRUE;}
        if (!particleITS) {particleITS=new TParticle();newparticleITS=kTRUE;}
        if (!particleMotherITS) {particleMotherITS=new TParticle();newparticleMotherITS=kTRUE;}
      }
      /////////////////////////

      //Double_t vtxX=vtxESD->GetX();
      //Double_t vtxY=vtxESD->GetY();
      //Double_t vtxZ=vtxESD->GetZ();

      AliESDVertex* pvtxESD = (AliESDVertex*)vtxESD->Clone();
      AliESDtrack* ptrack=(AliESDtrack*)track->Clone();
      AliExternalTrackParam* ptpcInnerC = (AliExternalTrackParam*)tpcInnerC->Clone();
      AliExternalTrackParam* ptrackInnerC = (AliExternalTrackParam*)trackInnerC->Clone();
      AliExternalTrackParam* ptrackInnerC2 = (AliExternalTrackParam*)trackInnerC2->Clone();
      AliExternalTrackParam* pouterITSc = (AliExternalTrackParam*)outerITSc->Clone();
      AliExternalTrackParam* ptrackInnerC3 = (AliExternalTrackParam*)trackInnerC3->Clone();
      Int_t ntracks = esdEvent->GetNumberOfTracks();
      
      if(fTreeSRedirector && dumpToTree) 
      {

        (*fTreeSRedirector)<<"highPt"<<
	  "fileName.="<<&fileName<<                // name of the chunk file (hopefully full)
	  "runNumber="<<runNumber<<                // runNumber
	  "evtTimeStamp="<<evtTimeStamp<<          // time stamp of event (in seconds)
	  "evtNumberInFile="<<evtNumberInFile<<    // event number
	  "triggerClass="<<&triggerClass<<         // trigger class as a string
	  "Bz="<<bz<<                              // solenoid magnetic field in the z direction (in kGaus)
	  "vtxESD.="<<pvtxESD<<                    // vertexer ESD tracks (can be biased by TPC pileup tracks)
	  //"vtxESDx="<<vtxX<<
	  //"vtxESDy="<<vtxY<<
	  //"vtxESDz="<<vtxZ<<
	  "IRtot="<<ir1<<                         // interaction record (trigger) counters - coutner 1
	  "IRint2="<<ir2<<                        // interaction record (trigger) coutners - counter 2
	  "mult="<<mult<<                         // multiplicity of tracks pointing to the primary vertex
	  "ntracks="<<ntracks<<                   // number of the esd tracks (to take into account the pileup in the TPC)
	  "esdTrack.="<<ptrack<<                  // esdTrack as used in the physical analysis
	  "extTPCInnerC.="<<ptpcInnerC<<          // ??? 
	  "extInnerParamC.="<<ptrackInnerC<<      // ???
	  "extInnerParam.="<<ptrackInnerC2<<      // ???
	  "extOuterITS.="<<pouterITSc<<           // ???
	  "extInnerParamRef.="<<ptrackInnerC3<<   // ???
	  "chi2TPCInnerC="<<chi2(0,0)<<           // chi2   of tracks ???
	  "chi2InnerC="<<chi2trackC(0,0)<<        // chi2s  of tracks TPCinner to the combined
	  "chi2OuterITS="<<chi2OuterITS(0,0)<<    // chi2s  of tracks TPC at inner wall to the ITSout
	  "centralityF="<<centralityF;
        if (fMC)
        {
          (*fTreeSRedirector)<<"highPt"<<
          "refTPCIn.="<<refTPCIn<<
          "refITS.="<<refITS<<
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
        (*fTreeSRedirector)<<"highPt"<<"\n";
      }

      delete pvtxESD;
      delete ptrack;
      delete ptpcInnerC;
      delete ptrackInnerC;
      delete ptrackInnerC2;
      delete pouterITSc;
      delete ptrackInnerC3;

      ////////////////////
      //delete the dummy objects we might have created.
      if (newvtxESD) {delete vtxESD; vtxESD=NULL;}
      if (newtrack) {delete track; track=NULL;}
      if (newtpcInnerC) {delete tpcInnerC; tpcInnerC=NULL;}
      if (newtrackInnerC) {delete trackInnerC; trackInnerC=NULL;}
      if (newtrackInnerC2) {delete trackInnerC2; trackInnerC2=NULL;}
      if (newouterITSc) {delete outerITSc; outerITSc=NULL;}
      if (newtrackInnerC3) {delete trackInnerC3; trackInnerC3=NULL;}
      if (newrefTPCIn) {delete refTPCIn; refTPCIn=NULL;}
      if (newrefITS) {delete refITS; refITS=NULL;}
      if (newparticle) {delete particle; particle=NULL;}
      if (newparticleMother) {delete particleMother; particleMother=NULL;}
      if (newparticleTPC) {delete particleTPC; particleTPC=NULL;}
      if (newparticleMotherTPC) {delete particleMotherTPC; particleMotherTPC=NULL;}
      if (newparticleITS) {delete particleITS; particleITS=NULL;}
      if (newparticleMotherITS) {delete particleMotherITS; particleMotherITS=NULL;}
      ///////////////

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
  AliInfo("we start!");

  if(!esdEvent) {
    AliDebug(AliLog::kError, "esdEvent not available");
    return;
  }

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
  if (!inputHandler)
  {
    Printf("ERROR: Could not receive input handler");
    return;
  }

  // get file name
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { 
    Printf("ERROR: Could not receive input chain");
    return;
  }
  TObjString fileName(chain->GetCurrentFile()->GetName());

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

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD); 
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK, isEventTriggered);
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());

  TObjString triggerClass = esdEvent->GetFiredTriggerClasses().Data();

  // check event cuts
  if(isEventOK && isEventTriggered)
  {
    if(!stack) return;

    //
    // MC info
    //
    TParticle *particle=NULL;
    TParticle *particleMother=NULL;
    Int_t mech=-1;

    // reco event info
    Double_t vert[3] = {0}; 
    vert[0] = vtxESD->GetXv();
    vert[1] = vtxESD->GetYv();
    vert[2] = vtxESD->GetZv();
    Int_t mult = vtxESD->GetNContributors();
    Double_t bz = esdEvent->GetMagneticField();
    Double_t runNumber = esdEvent->GetRunNumber();
    Double_t evtTimeStamp = esdEvent->GetTimeStamp();
    Int_t evtNumberInFile = esdEvent->GetEventNumberInFile();

    // loop over MC stack
    for (Int_t iMc = 0; iMc < stack->GetNtrack(); ++iMc) 
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
      if(TMath::Exp(2*scalempt)<downscaleF) continue;

      // is particle in acceptance
      if(!accCuts->AcceptTrack(particle)) continue;

      // check if particle reconstructed
      Bool_t isRec = kFALSE;
      Int_t  trackIndex = -1;
      for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++)
      {

        AliESDtrack *track = esdEvent->GetTrack(iTrack);
        if(!track) continue;
        if(track->Charge()==0) continue;
        if(esdTrackCuts->AcceptTrack(track) && accCuts->AcceptTrack(track)) 
        {
          Int_t label =  TMath::Abs(track->GetLabel());
          if(label == iMc) {
            isRec = kTRUE;
            trackIndex = iTrack;
            break;
          }
        } 
      }

      // Store information in the output tree
      AliESDtrack *recTrack = NULL; 
      if(trackIndex>-1)  { 
        recTrack = esdEvent->GetTrack(trackIndex); 
      } else {
        recTrack = new AliESDtrack(); 
      } 

      particleMother = GetMother(particle,stack);
      mech = particle->GetUniqueID();

      //MC particle track length
      Double_t tpcTrackLength = 0.;
      AliMCParticle *mcParticle = (AliMCParticle*) mcEvent->GetTrack(iMc);
      if(mcParticle) {
        Int_t counter;
        tpcTrackLength = mcParticle->GetTPCTrackLength(bz,0.05,counter,3.0);
      } 


      //
      if(fTreeSRedirector) {
        (*fTreeSRedirector)<<"MCEffTree"<<
          "fileName.="<<&fileName<<
          "triggerClass.="<<&triggerClass<<
          "runNumber="<<runNumber<<
          "evtTimeStamp="<<evtTimeStamp<<
          "evtNumberInFile="<<evtNumberInFile<<
          "Bz="<<bz<<
          "vtxESD.="<<vtxESD<<
          "mult="<<mult<<
          "esdTrack.="<<recTrack<<
          "isRec="<<isRec<<
          "tpcTrackLength="<<tpcTrackLength<<
          "particle.="<<particle<<
          "particleMother.="<<particleMother<<
          "mech="<<mech<<
          "\n";
      }

      if(trackIndex <0 && recTrack) delete recTrack; recTrack=0;
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
void AliAnalysisTaskFilteredTree::ProcessV0(AliESDEvent *const esdEvent, AliMCEvent * const mcEvent, AliESDfriend *const /*esdFriend*/)
{
  //
  // Select real events with V0 (K0s and Lambda and Gamma) high-pT candidates
  //
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
  if (!inputHandler)
  {
    Printf("ERROR: Could not receive input handler");
    return;
  }
   
  // get file name
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { 
    Printf("ERROR: Could not receive input chain");
    return;
  }
  TObjString fileName(chain->GetCurrentFile()->GetName());

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

  // check event cuts
  if(isEventOK && isEventTriggered) {
  //
  // Dump the pt downscaled V0 into the tree
  // 
  Int_t ntracks = esdEvent->GetNumberOfTracks();
  Int_t nV0s = esdEvent->GetNumberOfV0s();
  Int_t run = esdEvent->GetRunNumber();
  Int_t time= esdEvent->GetTimeStamp();
  Int_t evNr=esdEvent->GetEventNumberInFile();
  Double_t bz = esdEvent->GetMagneticField();


  for (Int_t iv0=0; iv0<nV0s; iv0++){
    AliESDv0 * v0 = esdEvent->GetV0(iv0);
    if (!v0) continue;
    AliESDtrack * track0 = esdEvent->GetTrack(v0->GetIndex(0));
    AliESDtrack * track1 = esdEvent->GetTrack(v0->GetIndex(1));
    if (!track0) continue;
    if (!track1) continue;
    if (track0->GetSign()<0) {
      track1 = esdEvent->GetTrack(v0->GetIndex(0));
      track0 = esdEvent->GetTrack(v0->GetIndex(1));
    }
    //
    Bool_t isDownscaled = IsV0Downscaled(v0);
    if (isDownscaled) continue;
    AliKFParticle kfparticle; //
    Int_t type=GetKFParticle(v0,esdEvent,kfparticle);
    if (type==0) continue;   
    TObjString triggerClass = esdEvent->GetFiredTriggerClasses().Data();

    if(!fTreeSRedirector) return;
    (*fTreeSRedirector)<<"V0s"<<
      "isDownscaled="<<isDownscaled<<
      "triggerClass="<<&triggerClass<<      //  trigger
      "Bz="<<bz<<
      "fileName.="<<&fileName<<
      "runNumber="<<run<<
      "evtTimeStamp="<<time<<
      "evtNumberInFile="<<evNr<<
      "type="<<type<<
      "ntracks="<<ntracks<<
      "v0.="<<v0<<
      "kf.="<<&kfparticle<<
      "track0.="<<track0<<
      "track1.="<<track1<<
      "centralityF="<<centralityF<<
      "\n";
  }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskFilteredTree::ProcessdEdx(AliESDEvent *const esdEvent, AliMCEvent * const mcEvent, AliESDfriend *const /*esdFriend*/)
{
  //
  // Select real events with large TPC dEdx signal
  //
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

  // get file name
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { 
    Printf("ERROR: Could not receive input chain");
    return;
  }
  TObjString fileName(chain->GetCurrentFile()->GetName());

  // trigger
  Bool_t isEventTriggered = kTRUE;
  AliPhysicsSelection *physicsSelection = NULL;
  AliTriggerAnalysis* triggerAnalysis = NULL;

  // 
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!inputHandler)
  {
    Printf("ERROR: Could not receive input handler");
    return;
  }

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
    vert[0] = vtxESD->GetXv();
    vert[1] = vtxESD->GetYv();
    vert[2] = vtxESD->GetZv();
    Int_t mult = vtxESD->GetNContributors();
    Double_t bz = esdEvent->GetMagneticField();
    Double_t runNumber = esdEvent->GetRunNumber();
    Double_t evtTimeStamp = esdEvent->GetTimeStamp();
    Int_t evtNumberInFile = esdEvent->GetEventNumberInFile();

    // large dEdx
    for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++)
    {
      AliESDtrack *track = esdEvent->GetTrack(iTrack);
      if(!track) continue;
      if(track->Charge()==0) continue;
      if(!esdTrackCuts->AcceptTrack(track)) continue;
      if(!accCuts->AcceptTrack(track)) continue;

      if(!IsHighDeDxParticle(track)) continue;
      TObjString triggerClass = esdEvent->GetFiredTriggerClasses().Data();

      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"dEdx"<<
      "fileName.="<<&fileName<<
      "runNumber="<<runNumber<<
      "evtTimeStamp="<<evtTimeStamp<<
      "evtNumberInFile="<<evtNumberInFile<<
	"triggerClass="<<&triggerClass<<      //  trigger
      "Bz="<<bz<<
      "vtxESD.="<<vtxESD<<
      "mult="<<mult<<
      "esdTrack.="<<track<<
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
  //if ((track0->GetITSclusters(0))<2) return 0;
  //if ((track1->GetITSclusters(0))<2) return 0; 
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
  mother = stack->Particle(motherLabel); 

return mother;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskFilteredTree::IsFromConversion(const Int_t label, AliStack *const stack) 
{
  Bool_t isFromConversion = kFALSE;

  if(stack) {
    TParticle* particle = stack->Particle(label);

    if(particle && particle->GetPDG() && particle->GetPDG()->Charge()!=0) 
    {
       Int_t mech = particle->GetUniqueID(); // production mechanism 
       Bool_t isPrim = stack->IsPhysicalPrimary(label);

       Int_t motherLabel = TMath::Abs(particle->GetMother(0));  
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
Bool_t AliAnalysisTaskFilteredTree::IsFromMaterial(const Int_t label, AliStack *const stack) 
{
  Bool_t isFromMaterial = kFALSE;

  if(stack) {
    TParticle* particle = stack->Particle(label);

    if(particle && particle->GetPDG() && particle->GetPDG()->Charge()!=0) 
    {
       Int_t mech = particle->GetUniqueID(); // production mechanism 
       Bool_t isPrim = stack->IsPhysicalPrimary(label);

       Int_t motherLabel = TMath::Abs(particle->GetMother(0));  
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
Bool_t AliAnalysisTaskFilteredTree::IsFromStrangeness(const Int_t label, AliStack *const stack) 
{
  Bool_t isFromStrangeness = kFALSE;

  if(stack) {
    TParticle* particle = stack->Particle(label);

    if(particle && particle->GetPDG() && particle->GetPDG()->Charge()!=0) 
    {
       Int_t mech = particle->GetUniqueID(); // production mechanism 
       Bool_t isPrim = stack->IsPhysicalPrimary(label);

       Int_t motherLabel = TMath::Abs(particle->GetMother(0));  
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

  //// must be deleted to store trees
  //if(fTreeSRedirector)  delete fTreeSRedirector; fTreeSRedirector=0;

  //// open temporary file and copy trees to the ouptut container

  //TChain* chain = 0;
  ////
  //chain = new TChain("highPt");
  //if(chain) { 
  //  chain->Add("jotwinow_Temp_Trees.root");
  //  fHighPtTree = chain->CopyTree("1");
  //  delete chain; chain=0; 
  //}
  //if(fHighPtTree) fHighPtTree->Print();

  ////
  //chain = new TChain("V0s");
  //if(chain) { 
  //  chain->Add("jotwinow_Temp_Trees.root");
  //  fV0Tree = chain->CopyTree("1");
  //  delete chain; chain=0; 
  //}
  //if(fV0Tree) fV0Tree->Print();

  ////
  //chain = new TChain("dEdx");
  //if(chain) { 
  //  chain->Add("jotwinow_Temp_Trees.root");
  //  fdEdxTree = chain->CopyTree("1");
  //  delete chain; chain=0; 
  //}
  //if(fdEdxTree) fdEdxTree->Print();

  ////
  //chain = new TChain("Laser");
  //if(chain) { 
  //  chain->Add("jotwinow_Temp_Trees.root");
  //  fLaserTree = chain->CopyTree("1");
  //  delete chain; chain=0; 
  //}
  //if(fLaserTree) fLaserTree->Print();

  ////
  //chain = new TChain("MCEffTree");
  //if(chain) { 
  //  chain->Add("jotwinow_Temp_Trees.root");
  //  fMCEffTree = chain->CopyTree("1");
  //  delete chain; chain=0; 
  //}
  //if(fMCEffTree) fMCEffTree->Print();

  ////
  //chain = new TChain("CosmicPairs");
  //if(chain) { 
  //  chain->Add("jotwinow_Temp_Trees.root");
  //  fCosmicPairsTree = chain->CopyTree("1");
  //  delete chain; chain=0; 
  //}
  //if(fCosmicPairsTree) fCosmicPairsTree->Print();  


  //OpenFile(1);

  // Post output data.
  //PostData(1, fHighPtTree);
  //PostData(2, fV0Tree);
  //PostData(3, fdEdxTree);
  //PostData(4, fLaserTree);
  //PostData(5, fMCEffTree);
  //PostData(6, fCosmicPairsTree);
}

//_____________________________________________________________________________
void AliAnalysisTaskFilteredTree::Terminate(Option_t *) 
{
  // Called one at the end 
  /*
  fOutputSummary = dynamic_cast<TTree*> (GetOutputData(1));
  if(fOutputSummary) delete fOutputSummary; fOutputSummary=0;
  TChain* chain = new TChain("highPt");
  if(!chain) return;
  chain->Add("jotwinow_HighPt_TrackAndV0_Trees.root");
  TTree *tree = chain->CopyTree("1");
  if (chain) { delete chain; chain=0; }
  if(!tree) return;
  tree->Print();
  fOutputSummary = tree;

  if (!fOutputSummary) {
    Printf("ERROR: AliAnalysisTaskFilteredTree::Terminate(): Output data not avaiable %p \n", GetOutputData(1));
    return;
  }
  */
  
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


