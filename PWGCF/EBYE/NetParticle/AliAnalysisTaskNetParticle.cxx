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
 * Class for for NetParticle Distributions
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

  fHnQA(NULL),

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

  if (fHelper)           delete fHelper;
  if (fEffCont)          delete fEffCont;
  if (fDCA)              delete fDCA;
  if (fDist)             delete fDist;
}

/*
 * ---------------------------------------------------------------------------------
 *                                 Public Methods
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisTaskNetParticle::UserCreateOutputObjects() {
  // Create histograms

  // -- Initialize all classes
  Initialize();

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
  // -- Get event / trigger statistics histograms
  // ------------------------------------------------------------------
  fOutList->Add(fHelper->GetHEventStat0());
  fOutList->Add(fHelper->GetHEventStat1());
  fOutList->Add(fHelper->GetHTriggerStat());
  fOutList->Add(fHelper->GetHCentralityStat());

  // ------------------------------------------------------------------
  // -- Add histograms from distribution class
  // ------------------------------------------------------------------
  if (fModeDistCreation == 1)
    fDist->CreateHistograms(fOutList);

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
  // -- Create THnSparseF - QA
  // ------------------------------------------------------------------
  if (fModeQACreation == 1)
    CreateQAHists();

  // ------------------------------------------------------------------

  TH1::AddDirectory(oldStatus);

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
    FillQAHists();

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
  fESDTrackCutsBase->SetMinNCrossedRowsTPC(70);                                             // TPC
  fESDTrackCutsBase->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);                    // TPC

  fESDTrackCutsBase->SetMaxChi2PerClusterTPC(4);                                            // TPC
  fESDTrackCutsBase->SetAcceptKinkDaughters(kFALSE);                                        // TPC
  fESDTrackCutsBase->SetRequireTPCRefit(kTRUE);                                             // TPC

  fESDTrackCutsBase->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff); // ITS
  fESDTrackCutsBase->SetClusterRequirementITS(AliESDtrackCuts::kSDD,AliESDtrackCuts::kOff); // ITS
  fESDTrackCutsBase->SetClusterRequirementITS(AliESDtrackCuts::kSSD,AliESDtrackCuts::kOff); // ITS

  fESDTrackCutsBase->SetDCAToVertex2D(kFALSE);                                              // VertexConstrained 
  fESDTrackCutsBase->SetRequireSigmaToVertex(kFALSE);                                       // VertexConstrained 

  fESDTrackCutsBase->SetEtaRange(-1.*fEtaMax, fEtaMax);                                     // Acceptance
  fESDTrackCutsBase->SetPtRange(fPtRange[0],fPtRange[1]);                                   // Acceptance

  // -- Mode : clean cuts -> small contamination
  if (fESDTrackCutMode == 0) {
    sModeName = "Clean";
    fESDTrackCutsBase->SetRequireITSRefit(kTRUE);                                           // ITS
    fESDTrackCutsBase->SetMaxChi2PerClusterITS(36);                                         // ITS
  }  
  // -- Mode : dirty cuts -> high efficiency
  else if (fESDTrackCutMode == 1) {
    sModeName = "Dirty";
    fESDTrackCutsBase->SetRequireITSRefit(kFALSE);                                          // ITS
  }
  // -- Mode : Default
  else {
    sModeName = "Base";
  }

  fESDTrackCutsBase->SetName(Form("NetParticleCuts2010_%s",sModeName.Data()));

  // -- Create ESD BKG track cuts
  // ------------------------------
  fESDTrackCutsBkg = static_cast<AliESDtrackCuts*>(fESDTrackCutsBase->Clone());
  fESDTrackCutsBkg->SetName(Form("NetParticleCuts2010_%s_Bkg",sModeName.Data()));
  fESDTrackCutsBkg->SetMaxDCAToVertexZ(10.);                                                // VertexConstrained 
  
  // -- Create ESD track cuts
  // ------------------------------
  fESDTrackCuts = static_cast<AliESDtrackCuts*>(fESDTrackCutsBase->Clone());
  fESDTrackCuts->SetName(Form("NetParticleCuts2010_%s",sModeName.Data()));
  fESDTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");                         // VertexConstrained  ->  7*(0.0026+0.0050/pt^1.01)
  fESDTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);                                        // VertexConstrained
  fESDTrackCuts->SetMaxDCAToVertexZ(2);                                                     // VertexConstrained 

  // -- Create ESD Eff track cuts
  // ------------------------------
  fESDTrackCutsEff = static_cast<AliESDtrackCuts*>(fESDTrackCuts->Clone());
  fESDTrackCutsEff->SetName(Form("NetParticleCuts2010_%s_Eff",sModeName.Data()));
  fESDTrackCutsEff->SetPtRange(fPtRangeEff[0],fPtRangeEff[1]);                              // Acceptance
  fESDTrackCutsBase->SetEtaRange(-1.*fEtaMaxEff, fEtaMaxEff);                               // Acceptance

  // ------------------------------------------------------------------
  // -- Initialize Helper
  // ------------------------------------------------------------------
  if (fHelper->Initialize(fIsMC, fModeDistCreation))
    return -1;

  // ------------------------------------------------------------------
  // -- Create / Initialize Efficiency/Contamination
  // ------------------------------------------------------------------
  if ((fIsMC||fIsAOD) && fModeEffCreation == 1) {
    fEffCont = new AliAnalysisNetParticleEffCont;
    fEffCont->Initialize(fESDTrackCutsEff, fHelper, fAODtrackCutBit);
  }

  // ------------------------------------------------------------------
  // -- Create / Initialize DCA Determination
  // ------------------------------------------------------------------
  if (fModeDCACreation == 1) {
    fDCA = new AliAnalysisNetParticleDCA;
    fDCA->Initialize(fESDTrackCuts, fESDTrackCutsBkg, fHelper);
  }

  // ------------------------------------------------------------------
  // -- Create / Initialize DCA Determination
  // ------------------------------------------------------------------
  if (fModeDistCreation == 1) {
    fDist = new AliAnalysisNetParticleDistribution;
    fDist->Initialize(fHelper, fESDTrackCuts, fIsMC, fPtRange, fEtaMax, fAODtrackCutBit);
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

  // -- Setup Event for Helper / EffCont  / DCA / Dist classes
  // ------------------------------------------------------------------
  fHelper->SetupEvent(fESDHandler, fAODHandler, fMCEvent);

  if (fModeEffCreation && (fIsMC || fIsAOD) )
    fEffCont->SetupEvent(fESDHandler, fAODHandler, fMCEvent); 

  if (fModeDCACreation == 1)
    fDCA->SetupEvent(fESDHandler, fMCEvent);

  if (fModeDistCreation == 1)
    fDist->SetupEvent(fESDHandler, fAODHandler, fMCEvent);   // JMT need AOD?

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
  if (!fAOD->GetHeader()->GetCentralityP()) {
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

/*
 * ---------------------------------------------------------------------------------
 *                           Helper Methods - private
 * ---------------------------------------------------------------------------------
 */

//________________________________________________________________________
void AliAnalysisTaskNetParticle::CreateQAHists() {
  // -- Create QA histograms
  
  Double_t dCent[2] = {-0.5, 8.5};
  Int_t iCent       = 9;
  
  Double_t dEta[2]  = {-0.9, 0.9}; // -> 37 bins
  Int_t iEta        = Int_t((dEta[1]-dEta[0]) / 0.075) +1 ; 

  Double_t dRap[2]  = {-0.5, 0.5}; 
  Int_t iRap        = Int_t((dRap[1]-dRap[0]) / 0.075) +1 ; 

  Double_t dPhi[2]  = {0.0, TMath::TwoPi()};
  Int_t iPhi        = 42;

  Double_t dPt[2]   = {0.1, 3.0}; 
  Int_t iPt         = Int_t((dPt[1]-dPt[0]) / 0.075);

  Double_t dSign[2] = {-1.5, 1.5};
  Int_t iSign       = 3;

  //                      cent:isAccepted: pInner:     pt:     sign:tpcSignal:nSigmaTPC:nSigmaTOF:      eta:     phi:       y: dcar: dcaz: nSigmaCdd: nSigmaCzz  
  Int_t    bin[15] = {   iCent,       2  ,    iPt,    iPt,    iSign,      500,     50  ,     50  ,     iEta,    iPhi,    iRap,  50 ,  50 ,       50 ,       50 };
  Double_t min[15] = {dCent[0],      -0.5, dPt[0], dPt[0], dSign[0],       30,     -5.0,     -5.0,  dEta[0], dPhi[0], dRap[0], -10., -10.,      -10.,      -10.};
  Double_t max[15] = {dCent[1],       1.5, dPt[1], dPt[1], dSign[1],      500,      5.0,      5.0,  dEta[1], dPhi[1], dRap[1],  10.,  10.,       10.,       10.};
  
  fOutListQA->Add(new THnSparseF("fHnQA", "cent:isAccepted:pInner:pt:sign:tpcSignal:nSigmaTPC:nSigmaTOF:eta:phi:u:dcar:dcaz:nSigmaCdd:nSigmaCzz",
				 15, bin, min, max));
  
  fHnQA = static_cast<THnSparseF*>(fOutListQA->Last());
  fHnQA->Sumw2();
  
  fHnQA->GetAxis(0)->SetTitle("centrality");
  fHnQA->GetAxis(1)->SetTitle("isAccepted");
  fHnQA->GetAxis(2)->SetTitle("#it{p}_{Inner} (GeV/#it{c})");
  fHnQA->GetAxis(3)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fHnQA->GetAxis(4)->SetTitle("sign");
  
  fHnQA->GetAxis(5)->SetTitle("TPC signal");
  fHnQA->GetAxis(6)->SetTitle("n #sigma TPC");
  fHnQA->GetAxis(7)->SetTitle("n #sigma TOF");
  
  fHnQA->GetAxis(8)->SetTitle("#eta");
  fHnQA->GetAxis(9)->SetTitle("#varphi");
  fHnQA->GetAxis(10)->SetTitle("#it{y}");    
  
  fHnQA->GetAxis(11)->SetTitle("DCAr");
  fHnQA->GetAxis(12)->SetTitle("DCAz");
  
  fHnQA->GetAxis(13)->SetTitle("n #sigma #sqrt(Cdd)/DCAr");
  fHnQA->GetAxis(14)->SetTitle("n #sigma #sqrt(Czz)/DCAz");
  
  fHelper->BinLogAxis(fHnQA, 2);
  fHelper->BinLogAxis(fHnQA, 3);
}

//________________________________________________________________________
void AliAnalysisTaskNetParticle::FillQAHists() {
  // -- Process ESD tracks and fill QA histograms

  if (!fESD)
    return;

  for (Int_t idxTrack = 0; idxTrack < fESD->GetNumberOfTracks(); ++idxTrack) {
    AliESDtrack *track = fESD->GetTrack(idxTrack); 

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- Check track cuts
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    // -- Check if charged track is accepted for basic parameters
    if (!fHelper->IsTrackAcceptedBasicCharged(track))
      continue;
    
    // -- Check if accepted - BKG
    if (!fESDTrackCutsBkg->AcceptTrack(track)) 
      continue;

    // -- Check if accepted in rapidity window
    Double_t yP;
    if (!fHelper->IsTrackAcceptedRapidity(track, yP))
      continue;

    // -- Check if accepted bt PID from TPC or TPC+TOF
    Double_t pid[2];
    Bool_t isAcceptedPID = fHelper->IsTrackAcceptedPID(track, pid); 
    isAcceptedPID = isAcceptedPID; // JMT : NOT USED FOR NOW
    //    if (!isAcceptedPID)      // JMT : NOT USED FOR NOW
    //      continue;              // JMT : NOT USED FOR NOW

    // -- Check if accepted with thighter DCA cuts
    Bool_t isAcceptedDCA = fHelper->IsTrackAcceptedDCA(track);

    // -- Check track cuts
    Bool_t isAcceptedVertex = fESDTrackCuts->AcceptTrack(track);

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- Fill tracks in QA THnSparseF
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    if (fHnQA->GetEntries() > 5e5) 
      return;

    // -- Get dca r/z
    Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
    track->GetImpactParameters(dca,cov);
    
    Float_t dcaRoverCdd = ( TMath::Sqrt(cov[0]) != 0. )  ? dca[0]/TMath::Sqrt(cov[0]) : -9.99;
    Float_t dcaZoverCzz = ( TMath::Sqrt(cov[2]) != 0. )  ? dca[1]/TMath::Sqrt(cov[2]) : -9.99;
    
    Double_t aQA[15] = {
      Double_t(fHelper->GetCentralityBin()),        //  0 centrality 
      Double_t(isAcceptedVertex || isAcceptedDCA),  //  1 isAccepted -> Vertex || isAcceptedDCA
      track->GetTPCmomentum(),                      //  2 p at InnerParam
      track->Pt(),                                  //  3 pt
      track->GetSign(),                             //  4 sign
      track->GetTPCsignal(),                        //  5 TPC dE/dx
      pid[0],                                       //  6 n Sigma TPC
      pid[1],                                       //  7 n Sigma TOF
      track->Eta(),                                 //  8 eta
      track->Phi(),                                 //  9 phi
      yP,                                           // 10 rapidity   
      dca[0],                                       // 11 dca r
      dca[1],                                       // 12 dca z
      dcaRoverCdd,                                  // 13 sqrt(cov[dd])
      dcaZoverCzz,                                  // 14 sqrt(cov[zz])
    };
    
    fHnQA->Fill(aQA);
  
} // for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {

  return;
}

