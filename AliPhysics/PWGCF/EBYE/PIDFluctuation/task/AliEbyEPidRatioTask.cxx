/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: ALICE Offline.                                                 *
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

//=========================================================================//
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                  Date: Wed Jul  9 18:38:30 CEST 2014                    // 
//          New approch to find particle ratio to reduce memory            //
//                             (Test Only)                                 //
//        Copied from NetParticle Classes
//        Origin: Authors: Jochen Thaeder <jochen@thaeder.de>
//                         Michael Weber <m.weber@cern.ch>
//=========================================================================//

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "THashList.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliTracker.h" 
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliESDtrackCuts.h"
#include "AliKineTrackCuts.h"
#include "AliMCParticle.h"
#include "AliESDVZERO.h"
#include "AliEbyEPidRatioTask.h"
#include "AliGenEventHeader.h"
#include "AliCentrality.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

using namespace std;
ClassImp(AliEbyEPidRatioTask)
//________________________________________________________________________
AliEbyEPidRatioTask::AliEbyEPidRatioTask(const char *name) :
  AliAnalysisTaskSE(name),
  fHelper(NULL),
  fEffCont(NULL),
  fEffContExtra(NULL),
  fDCA(NULL),
  fDist(NULL),
  fQA(NULL),

  fOutList(NULL),
  fOutListEff(NULL),
  fOutListCont(NULL),
  fOutListDCA(NULL),
  fOutListQA(NULL),

  fESD(NULL), 
  fESDHandler(NULL),

  fESDTrackCutsBase(NULL),
  fESDTrackCuts(NULL),
  fESDTrackCutsBkg(NULL),
  fESDTrackCutsEff(NULL),

  fAOD(NULL), 
  fAODHandler(NULL),

  fIsMC(kFALSE),
  fIsRatio(kFALSE),
  fIsPtBin(kFALSE),
  fIsDetectorWise(kFALSE),
  fIsAOD(kFALSE),
  fIsSub(kFALSE),
  fIsBS(kFALSE),
  fIsPer(kFALSE),

  fIsEffExtra(kFALSE),

  fESDTrackCutMode(0),
  fModeEffCreation(0),
  fModeDCACreation(0),
  fModeDistCreation(0),
  fModeQACreation(0),

  fMCEvent(NULL),
  fMCStack(NULL),

  fEtaMax(0.8),
  fEtaMaxEff(0.9),
  fPtRange(),
  fPtRangeEff(),

  fAODtrackCutBit(1024) {
  
  AliLog::SetClassDebugLevel("AliEbyEPidRatioTask",10);

  fPtRange[0]    = 0.4;
  fPtRange[1]    = 0.8;
  fPtRangeEff[0] = 0.2;
  fPtRangeEff[1] = 1.6;

  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
}

//________________________________________________________________________
AliEbyEPidRatioTask::~AliEbyEPidRatioTask() {
  // Destructor

  if (fESDTrackCutsBase) delete fESDTrackCutsBase;
  if (fESDTrackCuts)     delete fESDTrackCuts;
  if (fESDTrackCutsBkg)  delete fESDTrackCutsBkg;
  if (fESDTrackCutsEff)  delete fESDTrackCutsEff;

  if (fEffCont)          delete fEffCont;
  if (fEffContExtra)     delete fEffContExtra;
  if (fDCA)              delete fDCA;
  if (fDist)             delete fDist;
  if (fQA)               delete fQA;
  if (fHelper)           delete fHelper;
}

void AliEbyEPidRatioTask::SetIsRatio(Int_t i) {
  if      (i == 1) { fIsRatio = 1;  fIsPtBin = 0; fIsDetectorWise = 0; }
  else if (i == 2) { fIsRatio = 0;  fIsPtBin = 1; fIsDetectorWise = 0; }
  else if (i == 3) { fIsRatio = 1;  fIsPtBin = 1; fIsDetectorWise = 0; }         
  else if (i == 4) { fIsRatio = 0;  fIsPtBin = 0; fIsDetectorWise = 1; }         
  else if (i == 5) { fIsRatio = 0;  fIsPtBin = 1; fIsDetectorWise = 1; }         
  else if (i == 6) { fIsRatio = 1;  fIsPtBin = 1; fIsDetectorWise = 1; } 
  else if (i == 7) { fIsSub   = 1;  fIsBS    = 0; fIsPtBin        = 0; }
  else if (i == 8) { fIsSub   = 1;  fIsBS    = 1; fIsRatio        = 1; }
  else if (i == 9) { fIsSub   = 1;  fIsBS    = 1; fIsPtBin        = 0; }
  else if (i ==10) { fIsSub   = 1;  fIsBS    = 0; fIsPtBin        = 1; }
  else if (i ==11) { fIsSub   = 0;  fIsBS    = 1; fIsPtBin        = 1; }
  else if (i ==12) { fIsSub   = 1;  fIsBS    = 1; fIsPtBin        = 1; }
  else             { fIsRatio = 0;  fIsPtBin = 0; fIsDetectorWise = 0; }   

  if (fModeDistCreation == 0) 
    Printf(">>>> Task: No Physics Variable <<<<"); 
  if (fModeDistCreation > 0 && fIsRatio)       
    Printf(">>>> Task: Setting Ratios      : ON <<<<"); 
  else Printf(">>>> Task: Setting Ratios      : OFF <<<<");
  if (fModeDistCreation > 0 && fIsPtBin)        
    Printf(">>>> Task: Setting Binwise     : ON <<<<"); 
  else Printf(">>>> Task: Setting Binwise     : OFF <<<<");
  if (fModeDistCreation > 0 && fIsDetectorWise) 
    Printf(">>>> Task: Setting TOF-TPC wise: ON <<<<"); 
  else Printf(">>>> Task: Setting TOF-TPC wise: OFF <<<<");


}

//________________________________________________________________________
void AliEbyEPidRatioTask::UserCreateOutputObjects() {
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fOutList = new TList;
  fOutList->SetName(GetName()) ;
  fOutList->SetOwner(kTRUE);
 
  fOutListEff = new TList;
  fOutListEff->SetName(Form("%s_eff",GetName()));
  fOutListEff->SetOwner(kTRUE) ;

  fOutListCont = new TList;
  fOutListCont->SetName(Form("%s_cont",GetName()));
  fOutListCont->SetOwner(kTRUE) ;

  fOutListDCA = new TList;
  fOutListDCA->SetName(Form("%s_dca",GetName()));
  fOutListDCA->SetOwner(kTRUE) ;
 
  fOutListQA = new TList;
  fOutListQA->SetName(Form("%s_qa",GetName()));
  fOutListQA->SetOwner(kTRUE) ;

  Initialize();

  fOutList->Add(new TList);
  TList *list = static_cast<TList*>(fOutList->Last());
  list->SetName(Form("fStat"));
  list->SetOwner(kTRUE);

  list->Add(fHelper->GetHEventStat0());
  list->Add(fHelper->GetHEventStat1());
  list->Add(fHelper->GetHTriggerStat());
  list->Add(fHelper->GetHCentralityPercentile());
  list->Add(fHelper->GetHCentralityPercentileAll());

  if ((fIsAOD||fIsMC) && fModeEffCreation == 1) {
    if (fIsEffExtra) {
      for (Int_t i = 0; i < 4; i++) {
	for (Int_t j = 0; j < 2; j++) {
	  if (fEffContExtra->GetHnEff(i,j))
	    fOutListEff->Add(fEffContExtra->GetHnEff(i,j));
	}
	for (Int_t j = 2; j < 4; j++) {
	  if (fEffContExtra->GetHnEff(i,j))
	    fOutListCont->Add(fEffContExtra->GetHnEff(i,j));
	}

      }
    }
    else {
      fOutListEff->Add(fEffCont->GetHnEffMc());
      fOutListEff->Add(fEffCont->GetHnEffRec());
      
      fOutListCont->Add(fEffCont->GetHnContMc());
      fOutListCont->Add(fEffCont->GetHnContRec());

    }
  }
  
  if (fModeDCACreation == 1)
    fOutListDCA->Add(fDCA->GetHnDCA());

  if (fModeQACreation == 1) {
    fOutListQA->Add(fQA->GetHnQAPid());
    fOutListQA->Add(fQA->GetHnQADca());
  }

  TH1::AddDirectory(oldStatus);

  PostData(1,fOutList);
  PostData(2,fOutListEff);
  PostData(3,fOutListCont);
  PostData(4,fOutListDCA);
  PostData(5,fOutListQA);

  return;
}

//________________________________________________________________________
void AliEbyEPidRatioTask::UserExec(Option_t *) {
  
  if (SetupEvent() < 0) {
    PostData(1,fOutList);
    PostData(2,fOutListEff);
    PostData(3,fOutListCont);
    PostData(4,fOutListDCA);
    PostData(5,fOutListQA);
    return;
  }


  if ((fIsMC||fIsAOD) && fModeEffCreation == 1) {
    if (fIsEffExtra)  fEffContExtra->Process();
    else fEffCont->Process();
  }

  if (fModeDCACreation == 1)
    fDCA->Process();

  if (fModeDistCreation > 0)
    fDist->Process();

  if (fModeQACreation == 1)
    fQA->Process();

  PostData(1,fOutList);
  PostData(2,fOutListEff);
  PostData(3,fOutListCont);
  PostData(4,fOutListDCA);
  PostData(5,fOutListQA);

  return;
}      

//________________________________________________________________________
void AliEbyEPidRatioTask::Terminate(Option_t *){
  // Terminate
}

Int_t AliEbyEPidRatioTask::Initialize() {
  // Initialize event

  // ------------------------------------------------------------------
  // -- ESD TrackCuts
  // ------------------------------------------------------------------
  TString sModeName("");
  
  // -- Create ESD track cuts
  // --------------------------
  fESDTrackCutsBase = new AliESDtrackCuts;
  
  if (fESDTrackCutMode == 0) {
    fESDTrackCutsBase->SetMinNCrossedRowsTPC(70);                                             // TPC
    fESDTrackCutsBase->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);                    // TPC
  }
  else if (fESDTrackCutMode == 1) {
    fESDTrackCutsBase->SetMinNClustersTPC(70);                                                // TPC  2010
  }

  fESDTrackCutsBase->SetMaxChi2PerClusterTPC(4);                                              // TPC  2010
  fESDTrackCutsBase->SetAcceptKinkDaughters(kFALSE);                                          // TPC  2010
  fESDTrackCutsBase->SetRequireTPCRefit(kTRUE);                                               // TPC  2010

  if (fESDTrackCutMode == 0) {
    fESDTrackCutsBase->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff); // ITS
    fESDTrackCutsBase->SetClusterRequirementITS(AliESDtrackCuts::kSDD,AliESDtrackCuts::kOff); // ITS
    fESDTrackCutsBase->SetClusterRequirementITS(AliESDtrackCuts::kSSD,AliESDtrackCuts::kOff); // ITS
  } 
  else if (fESDTrackCutMode == 1) {
    fESDTrackCutsBase->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny); // ITS 2010
  //  fESDTrackCutsBase->SetMinNClustersITS(4);
  }

  fESDTrackCutsBase->SetRequireITSRefit(kTRUE);                                               // ITS 2010 
  fESDTrackCutsBase->SetMaxChi2PerClusterITS(36);                                             // ITS 2010

  fESDTrackCutsBase->SetDCAToVertex2D(kFALSE);                                                // VertexConstrained  2010
  fESDTrackCutsBase->SetRequireSigmaToVertex(kFALSE);                                         // VertexConstrained  2010
  fESDTrackCutsBase->SetMaxDCAToVertexZ(2);                                                   // VertexConstrained  2010
 
  fESDTrackCutsBase->SetEtaRange(-1.*fEtaMax, fEtaMax);                                       // Acceptance
  fESDTrackCutsBase->SetPtRange(fPtRange[0],fPtRange[1]);                                     // Acceptance

  // -- Mode : standard cuts
  if (fESDTrackCutMode == 0) 
    sModeName = "Std";
  // -- Mode : for comparison to LF
  else if (fESDTrackCutMode == 1)
    sModeName = "LF";
  // -- Mode : Default
  else
    sModeName = "Base";
  
  fESDTrackCutsBase->SetName(Form("NetParticleCuts2010_%s",sModeName.Data()));

  // -- Create ESD track cuts -> Base + DCA
  // ------------------------------
  fESDTrackCuts = static_cast<AliESDtrackCuts*>(fESDTrackCutsBase->Clone());
  fESDTrackCuts->SetName(Form("NetParticleCuts2010_%s",sModeName.Data()));
  if (fESDTrackCutMode == 0) 
    fESDTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");                       // 2010 VertexConstrained  ->  7*(0.0026+0.0050/pt^1.01)
    //    fESDTrackCuts->SetMaxDCAToVertexXY(0.3);
  else if (fESDTrackCutMode == 1)
    fESDTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");                       // 2010 VertexConstrained  ->  7*(0.0026+0.0050/pt^1.01)

  //  fESDTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);                                    // golden cut off

  // -- Create ESD BKG track cuts -> Base + Acceptance(Eff)
  // ------------------------------
  fESDTrackCutsBkg = static_cast<AliESDtrackCuts*>(fESDTrackCutsBase->Clone());
  fESDTrackCutsBkg->SetName(Form("NetParticleCuts2010_%s_Bkg",sModeName.Data()));
  fESDTrackCutsBkg->SetPtRange(fPtRangeEff[0],fPtRangeEff[1]);                              // Acceptance
  fESDTrackCutsBkg->SetEtaRange(-1.*fEtaMaxEff, fEtaMaxEff);                                // Acceptance
  
  // -- Create ESD Eff track cuts -> Base + DCA + Acceptance(Eff)
  // ------------------------------
  fESDTrackCutsEff = static_cast<AliESDtrackCuts*>(fESDTrackCuts->Clone());
  fESDTrackCutsEff->SetName(Form("NetParticleCuts2010_%s_Eff",sModeName.Data()));
  fESDTrackCutsEff->SetPtRange(fPtRangeEff[0],fPtRangeEff[1]);                              // Acceptance
  fESDTrackCutsEff->SetEtaRange(-1.*fEtaMaxEff, fEtaMaxEff);                                // Acceptance

  // ------------------------------------------------------------------
  // -- Initialize Helper
  // ------------------------------------------------------------------

  if (fHelper->Initialize(fESDTrackCutsEff, fIsMC,fIsRatio,fIsPtBin, fIsDetectorWise, fAODtrackCutBit, fModeDistCreation))
    return -1;

  // fHelper->SetIsRatio(fIsRatio);  
  // fHelper->SetIsPtBin(fIsPtBin);  

  // ------------------------------------------------------------------
  // -- Create / Initialize Efficiency/Contamination
  // ------------------------------------------------------------------
  if ((fIsMC||fIsAOD) && fModeEffCreation == 1) {
    if (fIsEffExtra) {
      fEffContExtra = new AliEbyEPidRatioEffContExtra;
      fEffContExtra->Initialize(fHelper, fESDTrackCutsEff);
      Printf(" >>>> AliEbyEPidRatioEffContExtra::Initialize()-ed  ");
    } else {
      fEffCont = new AliEbyEPidRatioEffCont;
      fEffCont->Initialize(fHelper, fESDTrackCutsEff);
      Printf(" >>>> AliEbyEPidRatioEffCont::Initialize()-ed  ");
    }
  }

  // ------------------------------------------------------------------
  // -- Create / Initialize DCA Determination
  // ------------------------------------------------------------------
  if (fModeDCACreation == 1) {
    fDCA = new AliEbyEPidRatioDCA;
    fDCA->SetESDTrackCutsBkg(fESDTrackCutsBkg);
    fDCA->Initialize(fHelper, fESDTrackCutsEff);
    Printf(" >>>> AliEbyEPidRatioDCA:Initialize()-ed  ");
  }

  // ------------------------------------------------------------------
  // -- Create / Initialize Phy Determination
  // ------------------------------------------------------------------
  if (fModeDistCreation > 0) {
    fDist = new AliEbyEPidRatioPhy;
    fDist->SetOutList(fOutList);
    if (fModeDistCreation == 2)  fDist->SetQA();
    if (fIsSub) fDist->SetSubRun();
    if (fIsBS) fDist->SetBSRun();
    if (fIsPer) fDist->SetIsPer();
    fDist->Initialize(fHelper, fESDTrackCuts);
    Printf(" >>>> AliEbyEPidRatioPhy:Initialize()-ed  ");
  }

  // ------------------------------------------------------------------
  // -- Create / Initialize QA Determination
  // ------------------------------------------------------------------
  if (fModeQACreation == 1) {
    fQA = new AliEbyEPidRatioQA();
    fQA->Initialize(fHelper, fESDTrackCutsEff);
    Printf(" >>>> AliEbyEPidRatioQA:Initialize()-ed  ");
  }

  // ------------------------------------------------------------------
  // -- Reset Event
  // ------------------------------------------------------------------
  ResetEvent();

  return 0;
}

//________________________________________________________________________
Int_t AliEbyEPidRatioTask::SetupEvent() {

  ResetEvent();

  // -- ESD Event
  // ------------------------------------------------------------------
  if (!fIsAOD && SetupESDEvent() < 0) {
    AliError("Setup ESD Event failed");
    return -1;
  }

  // -- AOD Event
  // ------------------------------------------------------------------
  if (fIsAOD && SetupAODEvent() < 0) {
    AliError("Setup AOD Event failed");
    return -1;
  }
  
  // -- Setup MC Event
  // ------------------------------------------------------------------
  if (fIsMC && SetupMCEvent() < 0) {
    AliError("Setup MC Event failed");
    return -1;
  }

  // -- Setup Event for Helper / EffCont  / DCA / Dist / QA classes
  // ------------------------------------------------------------------
  fHelper->SetupEvent(fESDHandler, fAODHandler, fMCEvent);

  
  if (fModeEffCreation && (fIsMC || fIsAOD) ) {
    if (fIsEffExtra) fEffContExtra->SetupEvent(); 
    else fEffCont->SetupEvent(); 
  }

  if (fModeDCACreation == 1)
    fDCA->SetupEvent();

  if (fModeDistCreation > 0)
    fDist->SetupEvent(); 

  if (fModeQACreation == 1)
    fQA->SetupEvent(); 

 
  // -- Evaluate Event cuts
  // ------------------------------------------------------------------
  return fHelper->IsEventRejected() ? -2 : 0;
}

//________________________________________________________________________
Int_t AliEbyEPidRatioTask::SetupESDEvent() {
  // -- Setup ESD Event
  // > return 0 for success 
  // > return -1 for failed setup


  

  fESDHandler= dynamic_cast<AliESDInputHandler*> 
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!fESDHandler) {
    AliError("Could not get ESD input handler");
    return -1;
  } 

  fESD = (AliESDEvent*)fESDHandler->GetEvent();
  if (!fESD) {
    AliError("Could not get ESD event");
    return -1;
  }

  // -- Check PID response
  // ------------------------------------------------------------------
  if (!fESDHandler->GetPIDResponse()) {
    AliError("Could not get PID response");
    return -1;
  } 

  // -- Check Vertex
  // ------------------------------------------------------------------
  if (!fESD->GetPrimaryVertexTracks()) {
    AliError("Could not get vertex from tracks");
    return -1;
  }

  // -- Check Centrality
  // ------------------------------------------------------------------
  if (!fESD->GetCentrality()) {
    AliError("Could not get centrality");
    return -1;
  }

  return 0;
}

//________________________________________________________________________
Int_t AliEbyEPidRatioTask::SetupAODEvent() {
  // -- Setup AOD Event
  // > return 0 for success 
  // > return -1 for failed setup

  fAODHandler= dynamic_cast<AliAODInputHandler*> 
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!fAODHandler) {
    AliError("Could not get AOD input handler");
    return -1;
  } 

  fAOD = fAODHandler->GetEvent();
  if (!fAOD) {
    AliError("Could not get AOD event");
    return -1;
  }

  // -- Check PID response
  // ------------------------------------------------------------------
  if (!fAODHandler->GetPIDResponse()) {
    AliError("Could not get PID response");
    return -1;
  } 

  // -- Check Vertex
  // ------------------------------------------------------------------
  if (!fAOD->GetPrimaryVertex()) {
    AliError("Could not get primary vertex");
    return -1;
  }

  // -- Check Centrality
  // ------------------------------------------------------------------
  if (!((AliVAODHeader*)fAOD->GetHeader())->GetCentralityP()) {
    AliError("Could not get centrality");
    return -1;
  }

  return 0;
}

//________________________________________________________________________
Int_t AliEbyEPidRatioTask::SetupMCEvent() {
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*> 
    (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  
  if (!mcH) {
    AliError("MC event handler not available");
    return -1;
  }

  fMCEvent = mcH->MCEvent();
  if (!fMCEvent) {
    AliError("MC event not available");
    return -1;
  }

  // -- Get MC header
  // ------------------------------------------------------------------
  AliHeader* header = fMCEvent->Header();
  if (!header) {
    AliError("MC header not available");
    return -1;
  }

  // -- Check Stack
  // ------------------------------------------------------------------
  fMCStack = fMCEvent->Stack(); 
  if (!fMCStack) {
    AliError("MC stack not available");
    return -1;
  }
    
  // -- Check GenHeader
  // ------------------------------------------------------------------
  if (!header->GenEventHeader()) {
    AliError("Could not retrieve genHeader from header");
    return -1;
  }

  // -- Check primary vertex
  // ------------------------------------------------------------------
  if (!fMCEvent->GetPrimaryVertex()){
    AliError("Could not get MC vertex");
    return -1;
  }

  return 0;
}

//________________________________________________________________________
void AliEbyEPidRatioTask::ResetEvent() {
  // -- Reset event
  
  // -- Reset ESD Event
  fESD       = NULL;

  // -- Reset AOD Event
  fAOD       = NULL;

  // -- Reset MC Event
  if (fIsMC)
    fMCEvent = NULL;

  // -- Reset Dist Creation 
  if (fModeDistCreation > 0)
    fDist->ResetEvent();

  return;
}


