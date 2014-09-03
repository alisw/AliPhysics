//-*- Mode: C++ -*-

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
#include "AliAnalysisTaskNetParticle.h"
#include "AliGenEventHeader.h"
#include "AliCentrality.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

using namespace std;

/**
 * Class for NetParticle Distributions
 * -- AnalysisTask
 * Authors: Jochen Thaeder <jochen@thaeder.de>
 *          Michael Weber <m.weber@cern.ch>
 */

ClassImp(AliAnalysisTaskNetParticle)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
AliAnalysisTaskNetParticle::AliAnalysisTaskNetParticle(const char *name) :
  AliAnalysisTaskSE(name),
  fHelper(NULL),
  fEffCont(NULL),
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
  fIsAOD(kFALSE),
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
  // Constructor   

  AliLog::SetClassDebugLevel("AliAnalysisTaskNetParticle",10);

  fPtRange[0] = 0.4;
  fPtRange[1] = 0.8;
  fPtRangeEff[0] = 0.2;
  fPtRangeEff[1] = 1.6;

  // -- Output slots
  // -------------------------------------------------
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskNetParticle::~AliAnalysisTaskNetParticle() {
  // Destructor

  if (fESDTrackCutsBase) delete fESDTrackCutsBase;
  if (fESDTrackCuts)     delete fESDTrackCuts;
  if (fESDTrackCutsBkg)  delete fESDTrackCutsBkg;
  if (fESDTrackCutsEff)  delete fESDTrackCutsEff;

  if (fEffCont)          delete fEffCont;
  if (fDCA)              delete fDCA;
  if (fDist)             delete fDist;
  if (fQA)               delete fQA;
  if (fHelper)           delete fHelper;
}

/*
 * ---------------------------------------------------------------------------------
 *                                 Public Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisTaskNetParticle::UserCreateOutputObjects() {
  // Create histograms

  // ------------------------------------------------------------------
  // -- Create Output Lists
  // ------------------------------------------------------------------
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

  // ------------------------------------------------------------------
  // -- Initialize all classes
  // ------------------------------------------------------------------
  Initialize();
 
  // ------------------------------------------------------------------
  // -- Get event / trigger statistics histograms
  // ------------------------------------------------------------------
  fOutList->Add(fHelper->GetHEventStat0());
  fOutList->Add(fHelper->GetHEventStat1());
  fOutList->Add(fHelper->GetHTriggerStat());
  fOutList->Add(fHelper->GetHCentralityStat());

  // ------------------------------------------------------------------
  // -- Add histograms from efficiency/contamination class
  // ------------------------------------------------------------------
  if ((fIsAOD||fIsMC) && fModeEffCreation == 1) {
    fOutListEff->Add(fEffCont->GetHnEff());
    fOutListCont->Add(fEffCont->GetHnCont());
  }

  // ------------------------------------------------------------------
  // -- Add histograms from DCA class
  // ------------------------------------------------------------------
  if (fModeDCACreation == 1)
    fOutListDCA->Add(fDCA->GetHnDCA());

  // ------------------------------------------------------------------
  // -- Add histograms from QA class
  // ------------------------------------------------------------------
  if (fModeQACreation == 1)
    fOutListQA->Add(fQA->GetHnQA());

  // ------------------------------------------------------------------

  TH1::AddDirectory(oldStatus);

  PostData(1,fOutList);
  PostData(2,fOutListEff);
  PostData(3,fOutListCont);
  PostData(4,fOutListDCA);
  PostData(5,fOutListQA);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskNetParticle::UserExec(Option_t *) {
  // Called for each event

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Setup Event
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  if (SetupEvent() < 0) {
    PostData(1,fOutList);
    PostData(2,fOutListEff);
    PostData(3,fOutListCont);
    PostData(4,fOutListDCA);
    PostData(5,fOutListQA);
    return;
  }

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Process Efficiency / Contamination Determination
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  if ((fIsMC||fIsAOD) && fModeEffCreation == 1)
    fEffCont->Process();

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Process DCA Determination
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  if (fModeDCACreation == 1)
    fDCA->Process();

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Process Distributions 
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  if (fModeDistCreation == 1)
    fDist->Process();

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Fill QA histograms
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  if (fModeQACreation == 1)
    fQA->Process();

  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  // -- Post output data
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

  PostData(1,fOutList);
  PostData(2,fOutListEff);
  PostData(3,fOutListCont);
  PostData(4,fOutListDCA);
  PostData(5,fOutListQA);

  return;
}      

//________________________________________________________________________
void AliAnalysisTaskNetParticle::Terminate(Option_t *){
  // Terminate
}

/*
 * ---------------------------------------------------------------------------------
 *                                 Public Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
Int_t AliAnalysisTaskNetParticle::Initialize() {
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
  if (fHelper->Initialize(fESDTrackCutsEff, fIsMC, fAODtrackCutBit, fModeDistCreation))
    return -1;

  // ------------------------------------------------------------------
  // -- Create / Initialize Efficiency/Contamination
  // ------------------------------------------------------------------
  if ((fIsMC||fIsAOD) && fModeEffCreation == 1) {
    fEffCont = new AliAnalysisNetParticleEffCont;
    fEffCont->Initialize(fHelper);
  }

  // ------------------------------------------------------------------
  // -- Create / Initialize DCA Determination
  // ------------------------------------------------------------------
  if (fModeDCACreation == 1) {
    fDCA = new AliAnalysisNetParticleDCA;
    fDCA->SetESDTrackCutsBkg(fESDTrackCutsBkg);
    fDCA->Initialize(fHelper);
  }

  // ------------------------------------------------------------------
  // -- Create / Initialize Distribution Determination
  // ------------------------------------------------------------------
  if (fModeDistCreation == 1) {
    fDist = new AliAnalysisNetParticleDistribution;
    fDist->SetOutList(fOutList);
    fDist->Initialize(fHelper, fESDTrackCuts);
  }

  // ------------------------------------------------------------------
  // -- Create / Initialize QA Determination
  // ------------------------------------------------------------------
  if (fModeQACreation == 1) {
    fQA = new AliAnalysisNetParticleQA();
    fQA->Initialize(fHelper);
  }

  // ------------------------------------------------------------------
  // -- Reset Event
  // ------------------------------------------------------------------
  ResetEvent();

  return 0;
}

/*
 * ---------------------------------------------------------------------------------
 *                            Setup/Reset Methods - private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
Int_t AliAnalysisTaskNetParticle::SetupEvent() {
  // Setup Reading of event
  // > return 0 for success / accepted event
  // > return -1 for failed setup
  // > return -2 for rejected event

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

  if (fModeEffCreation && (fIsMC || fIsAOD) )
    fEffCont->SetupEvent(); 

  if (fModeDCACreation == 1)
    fDCA->SetupEvent();

  if (fModeDistCreation == 1)
    fDist->SetupEvent(); 

  if (fModeQACreation == 1)
    fQA->SetupEvent(); 

  // -- Evaluate Event cuts
  // ------------------------------------------------------------------
  return fHelper->IsEventRejected() ? -2 : 0;
}

//________________________________________________________________________
Int_t AliAnalysisTaskNetParticle::SetupESDEvent() {
  // -- Setup ESD Event
  // > return 0 for success 
  // > return -1 for failed setup

  fESDHandler= dynamic_cast<AliESDInputHandler*> 
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!fESDHandler) {
    AliError("Could not get ESD input handler");
    return -1;
  } 

  fESD = fESDHandler->GetEvent();
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
Int_t AliAnalysisTaskNetParticle::SetupAODEvent() {
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
Int_t AliAnalysisTaskNetParticle::SetupMCEvent() {
  // -- Setup MC Event
  // > return 0 for success 
  // > return -1 for failed setup

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
void AliAnalysisTaskNetParticle::ResetEvent() {
  // -- Reset event
  
  // -- Reset ESD Event
  fESD       = NULL;

  // -- Reset AOD Event
  fAOD       = NULL;

  // -- Reset MC Event
  if (fIsMC)
    fMCEvent = NULL;

  // -- Reset Dist Creation 
  if (fModeDistCreation == 1)
    fDist->ResetEvent();

  return;
}


