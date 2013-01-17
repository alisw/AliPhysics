#include "TChain.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"                  
#include "TH3D.h"
#include "TArrayF.h"
#include "TF1.h"
#include "TRandom.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliEventplane.h"
#include "AliTHn.h"    

#include "AliEventPoolManager.h"           

#include "AliPID.h"                
#include "AliPIDResponse.h"        
#include "AliPIDCombined.h"        

#include "AliAnalysisTaskBFPsi.h"
#include "AliBalancePsi.h"
#include "AliAnalysisTaskTriggeredBF.h"


// Analysis task for the BF vs Psi code
// Authors: Panos.Christakoglou@nikhef.nl

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskBFPsi)

//________________________________________________________________________
AliAnalysisTaskBFPsi::AliAnalysisTaskBFPsi(const char *name) 
: AliAnalysisTaskSE(name), 
  fBalance(0),
  fRunShuffling(kFALSE),
  fShuffledBalance(0),
  fRunMixing(kFALSE),
  fRunMixingEventPlane(kFALSE),
  fMixingTracks(50000),
  fMixedBalance(0),
  fPoolMgr(0),
  fList(0),
  fListBF(0),
  fListBFS(0),
  fListBFM(0),
  fHistListPIDQA(0),
  fHistEventStats(0),
  fHistCentStats(0),
  fHistTriggerStats(0),
  fHistTrackStats(0),
  fHistVx(0),
  fHistVy(0),
  fHistVz(0),
  fHistEventPlane(0),
  fHistClus(0),
  fHistDCA(0),
  fHistChi2(0),
  fHistPt(0),
  fHistEta(0),
  fHistRapidity(0),
  fHistPhi(0),
  fHistPhiBefore(0),
  fHistPhiAfter(0),
  fHistPhiPos(0),
  fHistPhiNeg(0),
  fHistV0M(0),
  fHistRefTracks(0),
  fHistdEdxVsPTPCbeforePID(NULL),
  fHistBetavsPTOFbeforePID(NULL), 
  fHistProbTPCvsPtbeforePID(NULL), 
  fHistProbTOFvsPtbeforePID(NULL), 
  fHistProbTPCTOFvsPtbeforePID(NULL),
  fHistNSigmaTPCvsPtbeforePID(NULL), 
  fHistNSigmaTOFvsPtbeforePID(NULL), 
  fHistdEdxVsPTPCafterPID(NULL),
  fHistBetavsPTOFafterPID(NULL), 
  fHistProbTPCvsPtafterPID(NULL), 
  fHistProbTOFvsPtafterPID(NULL), 
  fHistProbTPCTOFvsPtafterPID(NULL),
  fHistNSigmaTPCvsPtafterPID(NULL), 
  fHistNSigmaTOFvsPtafterPID(NULL),  
  fPIDResponse(0x0),
  fPIDCombined(0x0),
  fParticleOfInterest(kPion),
  fPidDetectorConfig(kTPCTOF),
  fUsePID(kFALSE),
  fUsePIDnSigma(kTRUE),
  fUsePIDPropabilities(kFALSE), 
  fPIDNSigma(3.),
  fMinAcceptedPIDProbability(0.8),
  fESDtrackCuts(0),
  fCentralityEstimator("V0M"),
  fUseCentrality(kFALSE),
  fCentralityPercentileMin(0.), 
  fCentralityPercentileMax(5.),
  fImpactParameterMin(0.),
  fImpactParameterMax(20.),
  fUseMultiplicity(kFALSE),
  fNumberOfAcceptedTracksMin(0),
  fNumberOfAcceptedTracksMax(10000),
  fHistNumberOfAcceptedTracks(0),
  fUseOfflineTrigger(kFALSE),
  fVxMax(0.3),
  fVyMax(0.3),
  fVzMax(10.),
  nAODtrackCutBit(128),
  fPtMin(0.3),
  fPtMax(1.5),
  fEtaMin(-0.8),
  fEtaMax(-0.8),
  fDCAxyCut(-1),
  fDCAzCut(-1),
  fTPCchi2Cut(-1),
  fNClustersTPCCut(-1),
  fAcceptanceParameterization(0),
  fDifferentialV2(0),
  fUseFlowAfterBurner(kFALSE),
  fExcludeResonancesInMC(kFALSE),
  fUseMCPdgCode(kFALSE),
  fPDGCodeToBeAnalyzed(-1),
  fEventClass("EventPlane") {
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskBFPsi::~AliAnalysisTaskBFPsi() {

  // delete fBalance; 
  // delete fShuffledBalance; 
  // delete fList;
  // delete fListBF; 
  // delete fListBFS;

  // delete fHistEventStats; 
  // delete fHistTrackStats; 
  // delete fHistVx; 
  // delete fHistVy; 
  // delete fHistVz; 

  // delete fHistClus;
  // delete fHistDCA;
  // delete fHistChi2;
  // delete fHistPt;
  // delete fHistEta;
  // delete fHistPhi;
  // delete fHistV0M;
}

//________________________________________________________________________
void AliAnalysisTaskBFPsi::UserCreateOutputObjects() {
  // Create histograms
  // Called once

  // global switch disabling the reference 
  // (to avoid "Replacing existing TH1" if several wagons are created in train)
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  if(!fBalance) {
    fBalance = new AliBalancePsi();
    fBalance->SetAnalysisLevel("ESD");
    fBalance->SetEventClass(fEventClass);
    //fBalance->SetNumberOfBins(-1,16);
    //fBalance->SetInterval(-1,-0.8,0.8,16,0.,1.6,15.);
  }
  if(fRunShuffling) {
    if(!fShuffledBalance) {
      fShuffledBalance = new AliBalancePsi();
      fShuffledBalance->SetAnalysisLevel("ESD");
      fShuffledBalance->SetEventClass(fEventClass);
      //fShuffledBalance->SetNumberOfBins(-1,16);
      //fShuffledBalance->SetInterval(-1,-0.8,0.8,16,0.,1.6,15.);
    }
  }
  if(fRunMixing) {
    if(!fMixedBalance) {
      fMixedBalance = new AliBalancePsi();
      fMixedBalance->SetAnalysisLevel("ESD");
      fMixedBalance->SetEventClass(fEventClass);
    }
  }

  //QA list
  fList = new TList();
  fList->SetName("listQA");
  fList->SetOwner();

  //Balance Function list
  fListBF = new TList();
  fListBF->SetName("listBF");
  fListBF->SetOwner();

  if(fRunShuffling) {
    fListBFS = new TList();
    fListBFS->SetName("listBFShuffled");
    fListBFS->SetOwner();
  }

  if(fRunMixing) {
    fListBFM = new TList();
    fListBFM->SetName("listTriggeredBFMixed");
    fListBFM->SetOwner();
  }

  //PID QA list
  if(fUsePID) {
    fHistListPIDQA = new TList();
    fHistListPIDQA->SetName("listQAPID");
    fHistListPIDQA->SetOwner();
  }

  //Event stats.
  TString gCutName[5] = {"Total","Offline trigger",
                         "Vertex","Analyzed","sel. Centrality"};
  fHistEventStats = new TH2F("fHistEventStats",
                             "Event statistics;;Centrality percentile;N_{events}",
                             5,0.5,5.5,220,-5,105);
  for(Int_t i = 1; i <= 5; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fList->Add(fHistEventStats);

  TString gCentName[9] = {"V0M","FMD","TRK","TKL","CL0","CL1","V0MvsFMD","TKLvsV0M","ZEMvsZDC"};
  fHistCentStats = new TH2F("fHistCentStats",
                             "Centrality statistics;;Cent percentile",
			    9,-0.5,8.5,220,-5,105);
  for(Int_t i = 1; i <= 9; i++)
    fHistCentStats->GetXaxis()->SetBinLabel(i,gCentName[i-1].Data());
  fList->Add(fHistCentStats);

  fHistTriggerStats = new TH1F("fHistTriggerStats","Trigger statistics;TriggerBit;N_{events}",130,0,130);
  fList->Add(fHistTriggerStats);

  fHistTrackStats = new TH1F("fHistTrackStats","Event statistics;TrackFilterBit;N_{events}",130,0,130);
  fList->Add(fHistTrackStats);

  fHistNumberOfAcceptedTracks = new TH2F("fHistNumberOfAcceptedTracks",";N_{acc.};Centrality percentile;Entries",4001,-0.5,4000.5,220,-5,105);
  fList->Add(fHistNumberOfAcceptedTracks);

  // Vertex distributions
  fHistVx = new TH1F("fHistVx","Primary vertex distribution - x coordinate;V_{x} (cm);Entries",100,-0.5,0.5);
  fList->Add(fHistVx);
  fHistVy = new TH1F("fHistVy","Primary vertex distribution - y coordinate;V_{y} (cm);Entries",100,-0.5,0.5);
  fList->Add(fHistVy);
  fHistVz = new TH2F("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm);Centrality percentile;Entries",100,-20.,20.,220,-5,105);
  fList->Add(fHistVz);

  //Event plane
  fHistEventPlane = new TH2F("fHistEventPlane",";#Psi_{2} [deg.];Centrality percentile;Counts",100,0,360.,220,-5,105);
  fList->Add(fHistEventPlane);

  // QA histograms
  fHistClus = new TH2F("fHistClus","# Cluster (TPC vs. ITS)",10,0,10,200,0,200);
  fList->Add(fHistClus);
  fHistChi2 = new TH2F("fHistChi2","Chi2/NDF distribution;#chi^{2}/ndf;Centrality percentile",200,0,10,220,-5,105);
  fList->Add(fHistChi2);
  fHistDCA  = new TH2F("fHistDCA","DCA (xy vs. z)",400,-5,5,400,-5,5); 
  fList->Add(fHistDCA);
  fHistPt   = new TH2F("fHistPt","p_{T} distribution;p_{T} (GeV/c);Centrality percentile",200,0,10,220,-5,105);
  fList->Add(fHistPt);
  fHistEta  = new TH2F("fHistEta","#eta distribution;#eta;Centrality percentile",200,-2,2,220,-5,105);
  fList->Add(fHistEta);
  fHistRapidity  = new TH2F("fHistRapidity","y distribution;y;Centrality percentile",200,-2,2,220,-5,105);
  fList->Add(fHistRapidity);
  fHistPhi  = new TH2F("fHistPhi","#phi distribution;#phi (rad);Centrality percentile",200,0.0,2.*TMath::Pi(),220,-5,105);
  fList->Add(fHistPhi);
  fHistPhiBefore  = new TH2F("fHistPhiBefore","#phi distribution;#phi;Centrality percentile",200,0.,2*TMath::Pi(),220,-5,105);
  fList->Add(fHistPhiBefore);
  fHistPhiAfter  = new TH2F("fHistPhiAfter","#phi distribution;#phi;Centrality percentile",200,0.,2*TMath::Pi(),220,-5,105);
  fList->Add(fHistPhiAfter);
  fHistPhiPos  = new TH2F("fHistPhiPos","#phi distribution for positive particles;#phi;Centrality percentile",200,0.,2*TMath::Pi(),220,-5,105);
  fList->Add(fHistPhiPos);
  fHistPhiNeg  = new TH2F("fHistPhiNeg","#phi distribution for negative particles;#phi;Centrality percentile",200,0.,2.*TMath::Pi(),220,-5,105);
  fList->Add(fHistPhiNeg);
  fHistV0M  = new TH2F("fHistV0M","V0 Multiplicity C vs. A",500, 0, 20000, 500, 0, 20000);
  fList->Add(fHistV0M);
  TString gRefTrackName[6] = {"tracks","tracksPos","tracksNeg","tracksTPConly","clusITS0","clusITS1"};
  fHistRefTracks  = new TH2F("fHistRefTracks","Nr of Ref tracks/event vs. ref track estimator;;Nr of tracks",6, 0, 6, 400, 0, 20000);
  for(Int_t i = 1; i <= 6; i++)
    fHistRefTracks->GetXaxis()->SetBinLabel(i,gRefTrackName[i-1].Data());
  fList->Add(fHistRefTracks);

  // QA histograms for HBTinspired and Conversion cuts
  fList->Add(fBalance->GetQAHistHBTbefore());
  fList->Add(fBalance->GetQAHistHBTafter());
  fList->Add(fBalance->GetQAHistConversionbefore());
  fList->Add(fBalance->GetQAHistConversionafter());
  fList->Add(fBalance->GetQAHistPsiMinusPhi());

  // Balance function histograms
  // Initialize histograms if not done yet
  if(!fBalance->GetHistNp()){
    AliWarning("Histograms not yet initialized! --> Will be done now");
    AliWarning("--> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
    fBalance->InitHistograms();
  }

  if(fRunShuffling) {
    if(!fShuffledBalance->GetHistNp()) {
      AliWarning("Histograms (shuffling) not yet initialized! --> Will be done now");
      AliWarning("--> Add 'gBalance->InitHistograms()' in your configBalanceFunction");
      fShuffledBalance->InitHistograms();
    }
  }

  //for(Int_t a = 0; a < ANALYSIS_TYPES; a++){
  fListBF->Add(fBalance->GetHistNp());
  fListBF->Add(fBalance->GetHistNn());
  fListBF->Add(fBalance->GetHistNpn());
  fListBF->Add(fBalance->GetHistNnn());
  fListBF->Add(fBalance->GetHistNpp());
  fListBF->Add(fBalance->GetHistNnp());

  if(fRunShuffling) {
    fListBFS->Add(fShuffledBalance->GetHistNp());
    fListBFS->Add(fShuffledBalance->GetHistNn());
    fListBFS->Add(fShuffledBalance->GetHistNpn());
    fListBFS->Add(fShuffledBalance->GetHistNnn());
    fListBFS->Add(fShuffledBalance->GetHistNpp());
    fListBFS->Add(fShuffledBalance->GetHistNnp());
  }  

  if(fRunMixing) {
    fListBFM->Add(fMixedBalance->GetHistNp());
    fListBFM->Add(fMixedBalance->GetHistNn());
    fListBFM->Add(fMixedBalance->GetHistNpn());
    fListBFM->Add(fMixedBalance->GetHistNnn());
    fListBFM->Add(fMixedBalance->GetHistNpp());
    fListBFM->Add(fMixedBalance->GetHistNnp());
  }
  //}


  // Event Mixing
  if(fRunMixing){
    Int_t trackDepth = fMixingTracks; 
    Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
    
    // centrality bins
    Double_t centralityBins[] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,90.,100.}; // SHOULD BE DEDUCED FROM CREATED ALITHN!!!
    Double_t* centbins        = centralityBins;
    Int_t nCentralityBins     = sizeof(centralityBins) / sizeof(Double_t) - 1;
    
    // Zvtx bins
    Double_t vertexBins[] = {-10., -7., -5., -3., -1., 1., 3., 5., 7., 10.}; // SHOULD BE DEDUCED FROM CREATED ALITHN!!!
    Double_t* vtxbins     = vertexBins;
    Int_t nVertexBins     = sizeof(vertexBins) / sizeof(Double_t) - 1;
    
    // Event plane angle (Psi) bins
    Double_t psiBins[] = {0.,45.,135.,215.,305.,360.}; // SHOULD BE DEDUCED FROM CREATED ALITHN!!!
    Double_t* psibins     = psiBins;
    Int_t nPsiBins     = sizeof(psiBins) / sizeof(Double_t) - 1;
    
    // run the event mixing also in bins of event plane (statistics!)
    if(fRunMixingEventPlane){
      fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centbins, nVertexBins, vtxbins, nPsiBins, psibins);
    }
    else{
      fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centbins, nVertexBins, vtxbins);
    }
  }

  if(fESDtrackCuts) fList->Add(fESDtrackCuts);

  //====================PID========================//
  if(fUsePID) {
    fPIDCombined = new AliPIDCombined();
    fPIDCombined->SetDefaultTPCPriors();

    fHistdEdxVsPTPCbeforePID = new TH2D ("dEdxVsPTPCbefore","dEdxVsPTPCbefore", 1000, -10.0, 10.0, 1000, 0, 1000); 
    fHistListPIDQA->Add(fHistdEdxVsPTPCbeforePID); //addition 
    
    fHistBetavsPTOFbeforePID = new TH2D ("BetavsPTOFbefore","BetavsPTOFbefore", 1000, -10.0, 10., 1000, 0, 1.2); 
    fHistListPIDQA->Add(fHistBetavsPTOFbeforePID); //addition
    
    fHistProbTPCvsPtbeforePID = new TH2D ("ProbTPCvsPtbefore","ProbTPCvsPtbefore", 1000, -10.0,10.0, 1000, 0, 2.0); 
    fHistListPIDQA->Add(fHistProbTPCvsPtbeforePID); //addition 
    
    fHistProbTOFvsPtbeforePID = new TH2D ("ProbTOFvsPtbefore","ProbTOFvsPtbefore", 1000, -50, 50, 1000, 0, 2.0); 
    fHistListPIDQA->Add(fHistProbTOFvsPtbeforePID); //addition 

    fHistProbTPCTOFvsPtbeforePID =new TH2D ("ProbTPCTOFvsPtbefore","ProbTPCTOFvsPtbefore", 1000, -50, 50, 1000, 0, 2.0); 
    fHistListPIDQA->Add(fHistProbTPCTOFvsPtbeforePID); //addition 
    
    fHistNSigmaTPCvsPtbeforePID = new TH2D ("NSigmaTPCvsPtbefore","NSigmaTPCvsPtbefore", 1000, -10, 10, 1000, 0, 500); 
    fHistListPIDQA->Add(fHistNSigmaTPCvsPtbeforePID); //addition 
    
    fHistNSigmaTOFvsPtbeforePID = new TH2D ("NSigmaTOFvsPtbefore","NSigmaTOFvsPtbefore", 1000, -10, 10, 1000, 0, 500); 
    fHistListPIDQA->Add(fHistNSigmaTOFvsPtbeforePID); //addition 
    
    fHistdEdxVsPTPCafterPID = new TH2D ("dEdxVsPTPCafter","dEdxVsPTPCafter", 1000, -10, 10, 1000, 0, 1000); 
    fHistListPIDQA->Add(fHistdEdxVsPTPCafterPID); //addition 
    
    fHistBetavsPTOFafterPID = new TH2D ("BetavsPTOFafter","BetavsPTOFafter", 1000, -10, 10, 1000, 0, 1.2); 
    fHistListPIDQA->Add(fHistBetavsPTOFafterPID); //addition 
    
    fHistProbTPCvsPtafterPID = new TH2D ("ProbTPCvsPtafter","ProbTPCvsPtafter", 1000, -10, 10, 1000, 0, 2); 
    fHistListPIDQA->Add(fHistProbTPCvsPtafterPID); //addition 
  
    fHistProbTOFvsPtafterPID = new TH2D ("ProbTOFvsPtafter","ProbTOFvsPtafter", 1000,  -10, 10, 1000, 0, 2); 
    fHistListPIDQA->Add(fHistProbTOFvsPtafterPID); //addition  
    
    fHistProbTPCTOFvsPtafterPID =new TH2D ("ProbTPCTOFvsPtafter","ProbTPCTOFvsPtafter", 1000, -50, 50, 1000, 0, 2.0); 
    fHistListPIDQA->Add(fHistProbTPCTOFvsPtafterPID); //addition 

    fHistNSigmaTPCvsPtafterPID = new TH2D ("NSigmaTPCvsPtafter","NSigmaTPCvsPtafter", 1000, -10, 10, 1000, 0, 500); 
    fHistListPIDQA->Add(fHistNSigmaTPCvsPtafterPID); //addition  
    
    fHistNSigmaTOFvsPtafterPID = new TH2D ("NSigmaTOFvsPtafter","NSigmaTOFvsPtafter", 1000, -10, 10, 1000, 0, 500); 
    fHistListPIDQA->Add(fHistNSigmaTOFvsPtafterPID); //addition 
  }
  //====================PID========================//

  // Post output data.
  PostData(1, fList);
  PostData(2, fListBF);
  if(fRunShuffling) PostData(3, fListBFS);
  if(fRunMixing) PostData(4, fListBFM);
  if(fUsePID) PostData(5, fHistListPIDQA);       //PID

  TH1::AddDirectory(oldStatus);

}

//________________________________________________________________________
void AliAnalysisTaskBFPsi::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  TString gAnalysisLevel = fBalance->GetAnalysisLevel();
  Int_t gNumberOfAcceptedTracks = 0;
  Double_t fCentrality          = -1.;
  Double_t gReactionPlane       = -1.; 
  Float_t bSign = 0.;

  // get the event (for generator level: MCEvent())
  AliVEvent* eventMain = NULL;
  if(gAnalysisLevel == "MC") {
    eventMain = dynamic_cast<AliVEvent*>(MCEvent()); 
  }
  else{
    eventMain = dynamic_cast<AliVEvent*>(InputEvent()); 
    
    // for HBT like cuts need magnetic field sign
    bSign = (eventMain->GetMagneticField() > 0) ? 1 : -1;
  }
  if(!eventMain) {
    AliError("eventMain not available");
    return;
  }
  
  // PID Response task active?
  if(fUsePID) {
    fPIDResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
    if (!fPIDResponse) AliFatal("This Task needs the PID response attached to the inputHandler");
  }
  
  // check event cuts and fill event histograms
  if((fCentrality = IsEventAccepted(eventMain)) < 0){
    return;
  }
  
  //Compute Multiplicity or Centrality variable
  Double_t lMultiplicityVar = GetRefMultiOrCentrality( eventMain );

  // get the reaction plane
  gReactionPlane = GetEventPlane(eventMain);
  fHistEventPlane->Fill(gReactionPlane,fCentrality);
  if(gReactionPlane < 0){
    return;
  }
  
  // get the accepted tracks in main event
  TObjArray *tracksMain = GetAcceptedTracks(eventMain,fCentrality,gReactionPlane);
  gNumberOfAcceptedTracks = tracksMain->GetEntriesFast();

  //multiplicity cut (used in pp)
  fHistNumberOfAcceptedTracks->Fill(gNumberOfAcceptedTracks,fCentrality);
  if(fUseMultiplicity) {
    if((gNumberOfAcceptedTracks < fNumberOfAcceptedTracksMin)||(gNumberOfAcceptedTracks > fNumberOfAcceptedTracksMax))
      return;
  }

  // store charges of all accepted tracks, shuffle and reassign (two extra loops!)
  TObjArray* tracksShuffled = NULL;
  if(fRunShuffling){
    tracksShuffled = GetShuffledTracks(tracksMain);
  }
  
  // Event mixing 
  if (fRunMixing)
    {
      // 1. First get an event pool corresponding in mult (cent) and
      //    zvertex to the current event. Once initialized, the pool
      //    should contain nMix (reduced) events. This routine does not
      //    pre-scan the chain. The first several events of every chain
      //    will be skipped until the needed pools are filled to the
      //    specified depth. If the pool categories are not too rare, this
      //    should not be a problem. If they are rare, you could lose`
      //    statistics.
      
      // 2. Collect the whole pool's content of tracks into one TObjArray
      //    (bgTracks), which is effectively a single background super-event.
      
      // 3. The reduced and bgTracks arrays must both be passed into
      //    FillCorrelations(). Also nMix should be passed in, so a weight
      //    of 1./nMix can be applied.
      
      AliEventPool* pool = fPoolMgr->GetEventPool(fCentrality, eventMain->GetPrimaryVertex()->GetZ(),gReactionPlane);
      
      if (!pool){
	AliFatal(Form("No pool found for centrality = %f, zVtx = %f, psi = %f", fCentrality, eventMain->GetPrimaryVertex()->GetZ(),gReactionPlane));
      }
      else{
	
	//pool->SetDebug(1);

	if (pool->IsReady() || pool->NTracksInPool() > fMixingTracks / 10 || pool->GetCurrentNEvents() >= 5){ 
	  
	  
	  Int_t nMix = pool->GetCurrentNEvents();
	  //cout << "nMix = " << nMix << " tracks in pool = " << pool->NTracksInPool() << endl;
	  
	  //((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(2);
	  //((TH2F*) fListOfHistos->FindObject("mixedDist"))->Fill(centrality, pool->NTracksInPool());
	  //if (pool->IsReady())
	  //((TH1F*) fListOfHistos->FindObject("eventStat"))->Fill(3);
	  
	  // Fill mixed-event histos here  
	  for (Int_t jMix=0; jMix<nMix; jMix++) 
	    {
	      TObjArray* tracksMixed = pool->GetEvent(jMix);
	      fMixedBalance->CalculateBalance(gReactionPlane,tracksMain,tracksMixed,bSign,lMultiplicityVar);
	    }
	}
	
	// Update the Event pool
	pool->UpdatePool(tracksMain);
	//pool->PrintInfo();
	
      }//pool NULL check  
    }//run mixing
  
  // calculate balance function
  fBalance->CalculateBalance(gReactionPlane,tracksMain,NULL,bSign,lMultiplicityVar);
  
  // calculate shuffled balance function
  if(fRunShuffling && tracksShuffled != NULL) {
    fShuffledBalance->CalculateBalance(gReactionPlane,tracksShuffled,NULL,bSign,lMultiplicityVar);
  }
}      

//________________________________________________________________________
Double_t AliAnalysisTaskBFPsi::IsEventAccepted(AliVEvent *event){
  // Checks the Event cuts
  // Fills Event statistics histograms
  
  Bool_t isSelectedMain = kTRUE;
  Float_t fCentrality = -1.;
  TString gAnalysisLevel = fBalance->GetAnalysisLevel();

  fHistEventStats->Fill(1,fCentrality); //all events

  // Event trigger bits
  fHistTriggerStats->Fill(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());
  if(fUseOfflineTrigger)
    isSelectedMain = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
  if(isSelectedMain) {
    fHistEventStats->Fill(2,fCentrality); //triggered events
    
    //Centrality stuff 
    if(fUseCentrality) {
      if(gAnalysisLevel == "AOD") { //centrality in AOD header
	AliAODHeader *header = (AliAODHeader*) event->GetHeader();
	if(header){
	  fCentrality = header->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
	  
	  // QA for centrality estimators
	  fHistCentStats->Fill(0.,header->GetCentralityP()->GetCentralityPercentile("V0M"));
	  fHistCentStats->Fill(1.,header->GetCentralityP()->GetCentralityPercentile("FMD"));
	  fHistCentStats->Fill(2.,header->GetCentralityP()->GetCentralityPercentile("TRK"));
	  fHistCentStats->Fill(3.,header->GetCentralityP()->GetCentralityPercentile("TKL"));
	  fHistCentStats->Fill(4.,header->GetCentralityP()->GetCentralityPercentile("CL0"));
	  fHistCentStats->Fill(5.,header->GetCentralityP()->GetCentralityPercentile("CL1"));
	  fHistCentStats->Fill(6.,header->GetCentralityP()->GetCentralityPercentile("V0MvsFMD"));
	  fHistCentStats->Fill(7.,header->GetCentralityP()->GetCentralityPercentile("TKLvsV0M"));
	  fHistCentStats->Fill(8.,header->GetCentralityP()->GetCentralityPercentile("ZEMvsZDC"));
	  
	  // centrality QA (V0M)
	  fHistV0M->Fill(event->GetVZEROData()->GetMTotV0A(), event->GetVZEROData()->GetMTotV0C());
	  
	  // centrality QA (reference tracks)
	  fHistRefTracks->Fill(0.,header->GetRefMultiplicity());
	  fHistRefTracks->Fill(1.,header->GetRefMultiplicityPos());
	  fHistRefTracks->Fill(2.,header->GetRefMultiplicityNeg());
	  fHistRefTracks->Fill(3.,header->GetTPConlyRefMultiplicity());
	  fHistRefTracks->Fill(4.,header->GetNumberOfITSClusters(0));
	  fHistRefTracks->Fill(5.,header->GetNumberOfITSClusters(1));
	  fHistRefTracks->Fill(6.,header->GetNumberOfITSClusters(2));
	  fHistRefTracks->Fill(7.,header->GetNumberOfITSClusters(3));
	  fHistRefTracks->Fill(8.,header->GetNumberOfITSClusters(4));
	}//AOD header
      }//AOD

      else if(gAnalysisLevel == "ESD" || gAnalysisLevel == "MCESD"){ // centrality class for ESDs or MC-ESDs
 	AliCentrality *centrality = event->GetCentrality();
	fCentrality = centrality->GetCentralityPercentile(fCentralityEstimator.Data());

	// QA for centrality estimators
	fHistCentStats->Fill(0.,centrality->GetCentralityPercentile("V0M"));
	fHistCentStats->Fill(1.,centrality->GetCentralityPercentile("FMD"));
	fHistCentStats->Fill(2.,centrality->GetCentralityPercentile("TRK"));
	fHistCentStats->Fill(3.,centrality->GetCentralityPercentile("TKL"));
	fHistCentStats->Fill(4.,centrality->GetCentralityPercentile("CL0"));
	fHistCentStats->Fill(5.,centrality->GetCentralityPercentile("CL1"));
	fHistCentStats->Fill(6.,centrality->GetCentralityPercentile("V0MvsFMD"));
	fHistCentStats->Fill(7.,centrality->GetCentralityPercentile("TKLvsV0M"));
	fHistCentStats->Fill(8.,centrality->GetCentralityPercentile("ZEMvsZDC"));

	// centrality QA (V0M)
	fHistV0M->Fill(event->GetVZEROData()->GetMTotV0A(), event->GetVZEROData()->GetMTotV0C());
      }//ESD
      else if(gAnalysisLevel == "MC"){
	Double_t gImpactParameter = 0.;
	AliGenHijingEventHeader* headerH = dynamic_cast<AliGenHijingEventHeader*>(dynamic_cast<AliMCEvent*>(event)->GenEventHeader());
	if(headerH){
	  gImpactParameter = headerH->ImpactParameter();
	  fCentrality      = gImpactParameter;
	}//MC header
      }//MC
      else{
	fCentrality = -1.;
      }
    }
    
    // Event Vertex MC
    if(gAnalysisLevel == "MC"){
      AliGenEventHeader *header = dynamic_cast<AliMCEvent*>(event)->GenEventHeader();
      if(header){  
	TArrayF gVertexArray;
	header->PrimaryVertex(gVertexArray);
	//Printf("Vertex: %lf (x) - %lf (y) - %lf (z)",
	//gVertexArray.At(0),
	//gVertexArray.At(1),
	//gVertexArray.At(2));
	fHistEventStats->Fill(3,fCentrality); //events with a proper vertex
	if(TMath::Abs(gVertexArray.At(0)) < fVxMax) {
	  if(TMath::Abs(gVertexArray.At(1)) < fVyMax) {
	    if(TMath::Abs(gVertexArray.At(2)) < fVzMax) {
	      fHistEventStats->Fill(4,fCentrality); //analayzed events
	      fHistVx->Fill(gVertexArray.At(0));
	      fHistVy->Fill(gVertexArray.At(1));
	      fHistVz->Fill(gVertexArray.At(2),fCentrality);
	      
	      // take only events inside centrality class
	      if((fImpactParameterMin < fCentrality) && (fImpactParameterMax > fCentrality)){
		return fCentrality;	    
	      }//centrality class
	    }//Vz cut
	  }//Vy cut
	}//Vx cut
      }//header    
    }//MC
    
    // Event Vertex AOD, ESD, ESDMC
    else{
      const AliVVertex *vertex = event->GetPrimaryVertex();
      
      if(vertex) {
	Double32_t fCov[6];
	vertex->GetCovarianceMatrix(fCov);
	if(vertex->GetNContributors() > 0) {
	  if(fCov[5] != 0) {
	    fHistEventStats->Fill(3,fCentrality); //events with a proper vertex
	    if(TMath::Abs(vertex->GetX()) < fVxMax) {
	      if(TMath::Abs(vertex->GetY()) < fVyMax) {
		if(TMath::Abs(vertex->GetZ()) < fVzMax) {
		  fHistEventStats->Fill(4,fCentrality); //analyzed events
		  fHistVx->Fill(vertex->GetX());
		  fHistVy->Fill(vertex->GetY());
		  fHistVz->Fill(vertex->GetZ(),fCentrality);
		  
		  // take only events inside centrality class
		  if((fCentrality > fCentralityPercentileMin) && (fCentrality < fCentralityPercentileMax)){
		    fHistEventStats->Fill(5,fCentrality); //events with correct centrality
		    return fCentrality;		
		  }//centrality class
		}//Vz cut
	      }//Vy cut
	    }//Vx cut
	  }//proper vertex resolution
	}//proper number of contributors
      }//vertex object valid
    }//triggered event 
  }//AOD,ESD,ESDMC
  
  // in all other cases return -1 (event not accepted)
  return -1;
}


//________________________________________________________________________
Double_t AliAnalysisTaskBFPsi::GetRefMultiOrCentrality(AliVEvent *event){
    // Checks the Event cuts
    // Fills Event statistics histograms
  
  Float_t fCentrality = -1.;
  Double_t fMultiplicity = -100.;
  TString gAnalysisLevel = fBalance->GetAnalysisLevel();
  if(fEventClass == "Centrality"){
    Bool_t isSelectedMain = kTRUE;
    

    if(gAnalysisLevel == "AOD") { //centrality in AOD header
      AliAODHeader *header = (AliAODHeader*) event->GetHeader();
      if(header){
        fCentrality = header->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
      }//AOD header
    }//AOD
    
    else if(gAnalysisLevel == "ESD" || gAnalysisLevel == "MCESD"){ // centrality class for ESDs or MC-ESDs
      AliCentrality *centrality = event->GetCentrality();
      fCentrality = centrality->GetCentralityPercentile(fCentralityEstimator.Data());
    }//ESD
    else if(gAnalysisLevel == "MC"){
      Double_t gImpactParameter = 0.;
      AliGenHijingEventHeader* headerH = dynamic_cast<AliGenHijingEventHeader*>(dynamic_cast<AliMCEvent*>(event)->GenEventHeader());
      if(headerH){
        gImpactParameter = headerH->ImpactParameter();
        fCentrality      = gImpactParameter;
      }//MC header
    }//MC
    else{
      fCentrality = -1.;
    }
  }//End if "Centrality"
  if(fEventClass=="Multiplicity"&&gAnalysisLevel == "ESD"){
    fMultiplicity = fESDtrackCuts->GetReferenceMultiplicity(dynamic_cast<AliESDEvent*>(event), AliESDtrackCuts::kTrackletsITSTPC,0.5);
  }
  if(fEventClass=="Multiplicity"&&gAnalysisLevel != "ESD"){
    AliAODHeader *header = (AliAODHeader*) event->GetHeader();
    if(header){
      fMultiplicity = header->GetRefMultiplicity();
    }//AOD header
  }
  Double_t lReturnVal = -100;
  if(fEventClass=="Multiplicity"){
    lReturnVal = fMultiplicity;
  }else if(fEventClass=="Centrality"){
    lReturnVal = fCentrality;
  }
  return lReturnVal;
}

//________________________________________________________________________
Double_t AliAnalysisTaskBFPsi::GetEventPlane(AliVEvent *event){
  // Get the event plane

  TString gAnalysisLevel = fBalance->GetAnalysisLevel();

  Float_t gVZEROEventPlane    = -10.;
  Float_t gReactionPlane      = -10.;
  Double_t qxTot = 0.0, qyTot = 0.0;

  //MC: from reaction plane
  if(gAnalysisLevel == "MC"){
   
    AliGenHijingEventHeader* headerH = dynamic_cast<AliGenHijingEventHeader*>(dynamic_cast<AliMCEvent*>(event)->GenEventHeader());
    if (headerH) {
      gReactionPlane = headerH->ReactionPlaneAngle();
      //gReactionPlane *= TMath::RadToDeg();
    }
  }//MC
  
  // AOD,ESD,ESDMC: from VZERO Event Plane
  else{
   
    AliEventplane *ep = event->GetEventplane();
    if(ep){ 
      gVZEROEventPlane = ep->CalculateVZEROEventPlane(event,10,2,qxTot,qyTot);
      if(gVZEROEventPlane < 0.) gVZEROEventPlane += TMath::Pi();
      //gReactionPlane = gVZEROEventPlane*TMath::RadToDeg();
      gReactionPlane = gVZEROEventPlane;
    }
  }//AOD,ESD,ESDMC

  return gReactionPlane;
}

//________________________________________________________________________
TObjArray* AliAnalysisTaskBFPsi::GetAcceptedTracks(AliVEvent *event, Double_t fCentrality, Double_t gReactionPlane){
  // Returns TObjArray with tracks after all track cuts (only for AOD!)
  // Fills QA histograms

  TString gAnalysisLevel = fBalance->GetAnalysisLevel();

  //output TObjArray holding all good tracks
  TObjArray* tracksAccepted = new TObjArray;
  tracksAccepted->SetOwner(kTRUE);

  Double_t vCharge;
  Double_t vEta;
  Double_t vY;
  Double_t vPhi;
  Double_t vPt;


  if(gAnalysisLevel == "AOD") { // handling of TPC only tracks different in AOD and ESD
    // Loop over tracks in event
    for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++) {
      AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(event->GetTrack(iTracks));
      if (!aodTrack) {
	AliError(Form("Could not receive track %d", iTracks));
	continue;
      }
      
      // AOD track cuts
      
      // For ESD Filter Information: ANALYSIS/macros/AddTaskESDfilter.C
      // take only TPC only tracks 
      fHistTrackStats->Fill(aodTrack->GetFilterMap());
      if(!aodTrack->TestFilterBit(nAODtrackCutBit)) continue;
      
      vCharge = aodTrack->Charge();
      vEta    = aodTrack->Eta();
      vY      = aodTrack->Y();
      vPhi    = aodTrack->Phi();// * TMath::RadToDeg();
      vPt     = aodTrack->Pt();
      
      Float_t dcaXY = aodTrack->DCA();      // this is the DCA from global track (not exactly what is cut on)
      Float_t dcaZ  = aodTrack->ZAtDCA();   // this is the DCA from global track (not exactly what is cut on)
      
      
      // Kinematics cuts from ESD track cuts
      if( vPt < fPtMin || vPt > fPtMax)      continue;
      if( vEta < fEtaMin || vEta > fEtaMax)  continue;
      
      // Extra DCA cuts (for systematic studies [!= -1])
      if( fDCAxyCut != -1 && fDCAzCut != -1){
	if(TMath::Sqrt((dcaXY*dcaXY)/(fDCAxyCut*fDCAxyCut)+(dcaZ*dcaZ)/(fDCAzCut*fDCAzCut)) > 1 ){
	  continue;  // 2D cut
	}
      }
      
      // Extra TPC cuts (for systematic studies [!= -1])
      if( fTPCchi2Cut != -1 && aodTrack->Chi2perNDF() > fTPCchi2Cut){
	continue;
      }
      if( fNClustersTPCCut != -1 && aodTrack->GetTPCNcls() < fNClustersTPCCut){
	continue;
      }
      
      // fill QA histograms
      fHistClus->Fill(aodTrack->GetITSNcls(),aodTrack->GetTPCNcls());
      fHistDCA->Fill(dcaZ,dcaXY);
      fHistChi2->Fill(aodTrack->Chi2perNDF(),fCentrality);
      fHistPt->Fill(vPt,fCentrality);
      fHistEta->Fill(vEta,fCentrality);
      fHistRapidity->Fill(vY,fCentrality);
      if(vCharge > 0) fHistPhiPos->Fill(vPhi,fCentrality);
      else if(vCharge < 0) fHistPhiNeg->Fill(vPhi,fCentrality);
      fHistPhi->Fill(vPhi,fCentrality);
      
      // add the track to the TObjArray
      tracksAccepted->Add(new AliBFBasicParticle(vEta, vPhi, vPt, 1.*vCharge));
    }//track loop
  }// AOD analysis


  else if(gAnalysisLevel == "ESD" || gAnalysisLevel == "MCESD") { // handling of TPC only tracks different in AOD and ESD

    AliESDtrack *trackTPC   = NULL;
    AliMCParticle *mcTrack   = NULL;

    AliMCEvent*  mcEvent     = NULL;
    
    // for MC ESDs use also MC information
    if(gAnalysisLevel == "MCESD"){
      mcEvent = MCEvent(); 
      if (!mcEvent) {
	AliError("mcEvent not available");
	return tracksAccepted;
      }
    }
    
    // Loop over tracks in event
    for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++) {
      AliESDtrack* track = dynamic_cast<AliESDtrack *>(event->GetTrack(iTracks));
      if (!track) {
	AliError(Form("Could not receive track %d", iTracks));
	continue;
      }	

      // for MC ESDs use also MC information --> MC track not used further???
      if(gAnalysisLevel == "MCESD"){
	Int_t label = TMath::Abs(track->GetLabel());
	if(label > mcEvent->GetNumberOfTracks()) continue;
	if(label > mcEvent->GetNumberOfPrimaries()) continue;
	
	mcTrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(label));
	if(!mcTrack) continue;
      }

      // take only TPC only tracks
      trackTPC   = new AliESDtrack();
      if(!track->FillTPCOnlyTrack(*trackTPC)) continue;
      
      //ESD track cuts
      if(fESDtrackCuts) 
	if(!fESDtrackCuts->AcceptTrack(trackTPC)) continue;
      
      // fill QA histograms
      Float_t b[2];
      Float_t bCov[3];
      trackTPC->GetImpactParameters(b,bCov);
      if (bCov[0]<=0 || bCov[2]<=0) {
	AliDebug(1, "Estimated b resolution lower or equal zero!");
	bCov[0]=0; bCov[2]=0;
      }
      
      Int_t nClustersTPC = -1;
      nClustersTPC = trackTPC->GetTPCNclsIter1();   // TPC standalone
      //nClustersTPC = track->GetTPCclusters(0);   // global track
      Float_t chi2PerClusterTPC = -1;
      if (nClustersTPC!=0) {
	chi2PerClusterTPC = trackTPC->GetTPCchi2Iter1()/Float_t(nClustersTPC);      // TPC standalone
	//chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);     // global track
      }
      
      //===========================PID===============================//		    
      if(fUsePID) {
	Double_t prob[AliPID::kSPECIES]={0.};
	Double_t probTPC[AliPID::kSPECIES]={0.};
	Double_t probTOF[AliPID::kSPECIES]={0.};
	Double_t probTPCTOF[AliPID::kSPECIES]={0.};
	
	Double_t nSigma = 0.;
	UInt_t detUsedTPC = 0;
	UInt_t detUsedTOF = 0;
	UInt_t detUsedTPCTOF = 0;
	
	//Decide what detector configuration we want to use
	switch(fPidDetectorConfig) {
	case kTPCpid:
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
	  nSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)fParticleOfInterest));
	  detUsedTPC = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPC);
	  for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
	    prob[iSpecies] = probTPC[iSpecies];
	  break;
	case kTOFpid:
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
	  nSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)fParticleOfInterest));
	  detUsedTOF = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTOF);
	  for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
	    prob[iSpecies] = probTOF[iSpecies];
	  break;
	case kTPCTOF:
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
	  detUsedTPCTOF = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPCTOF);
	  for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
	    prob[iSpecies] = probTPCTOF[iSpecies];
	  break;
	default:
	  break;
	}//end switch: define detector mask
	
	//Filling the PID QA
	Double_t tofTime = -999., length = 999., tof = -999.;
	Double_t c = TMath::C()*1.E-9;// m/ns
	Double_t beta = -999.;
	Double_t  nSigmaTOFForParticleOfInterest = -999.;
	if ( (track->IsOn(AliESDtrack::kTOFin)) &&
	     (track->IsOn(AliESDtrack::kTIME))  ) { 
	  tofTime = track->GetTOFsignal();//in ps
	  length = track->GetIntegratedLength();
	  tof = tofTime*1E-3; // ns	
	  
	  if (tof <= 0) {
	    //Printf("WARNING: track with negative TOF time found! Skipping this track for PID checks\n");
	    continue;
	  }
	  if (length <= 0){
	    //printf("WARNING: track with negative length found!Skipping this track for PID checks\n");
	    continue;
	  }
	  
	  length = length*0.01; // in meters
	  tof = tof*c;
	  beta = length/tof;
	  
	  nSigmaTOFForParticleOfInterest = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)fParticleOfInterest);
	  fHistBetavsPTOFbeforePID ->Fill(track->P()*track->Charge(),beta);
	  fHistProbTOFvsPtbeforePID ->Fill(track->Pt(),probTOF[fParticleOfInterest]);
	  fHistNSigmaTOFvsPtbeforePID ->Fill(track->Pt(),nSigmaTOFForParticleOfInterest);
	}//TOF signal 
	
	
	Double_t  nSigmaTPCForParticleOfInterest = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)fParticleOfInterest);
	fHistdEdxVsPTPCbeforePID -> Fill(track->P()*track->Charge(),track->GetTPCsignal());
	fHistProbTPCvsPtbeforePID -> Fill(track->Pt(),probTPC[fParticleOfInterest]); 
	fHistNSigmaTPCvsPtbeforePID -> Fill(track->Pt(),nSigmaTPCForParticleOfInterest); 
	fHistProbTPCTOFvsPtbeforePID -> Fill(track->Pt(),probTPCTOF[fParticleOfInterest]);
	//end of QA-before pid
	
	if ((detUsedTPC != 0)||(detUsedTOF != 0)||(detUsedTPCTOF != 0)) {
	  //Make the decision based on the n-sigma
	  if(fUsePIDnSigma) {
	    if(nSigma > fPIDNSigma) continue;}
	  
	  //Make the decision based on the bayesian
	  else if(fUsePIDPropabilities) {
	    if(fParticleOfInterest != TMath::LocMax(AliPID::kSPECIES,prob)) continue;
	    if (prob[fParticleOfInterest] < fMinAcceptedPIDProbability) continue;      
	  }
	  
	  //Fill QA after the PID
	  fHistBetavsPTOFafterPID ->Fill(track->P()*track->Charge(),beta);
	  fHistProbTOFvsPtafterPID ->Fill(track->Pt(),probTOF[fParticleOfInterest]);
	  fHistNSigmaTOFvsPtafterPID ->Fill(track->Pt(),nSigmaTOFForParticleOfInterest);
	  
	  fHistdEdxVsPTPCafterPID -> Fill(track->P()*track->Charge(),track->GetTPCsignal());
	  fHistProbTPCvsPtafterPID -> Fill(track->Pt(),probTPC[fParticleOfInterest]); 
	  fHistProbTPCTOFvsPtafterPID -> Fill(track->Pt(),probTPCTOF[fParticleOfInterest]);
	  fHistNSigmaTPCvsPtafterPID -> Fill(track->Pt(),nSigmaTPCForParticleOfInterest); 
	}
      }
      //===========================PID===============================//
      vCharge = trackTPC->Charge();
      vY      = trackTPC->Y();
      vEta    = trackTPC->Eta();
      vPhi    = trackTPC->Phi();// * TMath::RadToDeg();
      vPt     = trackTPC->Pt();
      fHistClus->Fill(trackTPC->GetITSclusters(0),nClustersTPC);
      fHistDCA->Fill(b[1],b[0]);
      fHistChi2->Fill(chi2PerClusterTPC,fCentrality);
      fHistPt->Fill(vPt,fCentrality);
      fHistEta->Fill(vEta,fCentrality);
      fHistPhi->Fill(vPhi,fCentrality);
      fHistRapidity->Fill(vY,fCentrality);
      if(vCharge > 0) fHistPhiPos->Fill(vPhi,fCentrality);
      else if(vCharge < 0) fHistPhiNeg->Fill(vPhi,fCentrality);
      
      // add the track to the TObjArray
      tracksAccepted->Add(new AliBFBasicParticle(vEta, vPhi, vPt, vCharge));      

      delete trackTPC;
    }//track loop
  }// ESD analysis

  else if(gAnalysisLevel == "MC"){
    
    // Loop over tracks in event
    for (Int_t iTracks = 0; iTracks < dynamic_cast<AliMCEvent*>(event)->GetNumberOfPrimaries(); iTracks++) {
      AliMCParticle* track = dynamic_cast<AliMCParticle *>(event->GetTrack(iTracks));
      if (!track) {
	AliError(Form("Could not receive particle %d", iTracks));
	continue;
      }
	    
      //exclude non stable particles
      if(!(dynamic_cast<AliMCEvent*>(event)->IsPhysicalPrimary(iTracks))) continue;

      vEta    = track->Eta();
      vPt     = track->Pt();
      vY      = track->Y();
      
      if( vPt < fPtMin || vPt > fPtMax)      
	continue;
      if (!fUsePID) {
	if( vEta < fEtaMin || vEta > fEtaMax)  continue;
      }
      else if (fUsePID){
	if( vY < fEtaMin || vY > fEtaMax)  continue;
      }
      
      //analyze one set of particles
      if(fUseMCPdgCode) {
	TParticle *particle = track->Particle();
	if(!particle) continue;
	
	Int_t gPdgCode = particle->GetPdgCode();
	if(TMath::Abs(fPDGCodeToBeAnalyzed) != TMath::Abs(gPdgCode)) 
	  continue;
      }
      
      //Use the acceptance parameterization
      if(fAcceptanceParameterization) {
	Double_t gRandomNumber = gRandom->Rndm();
	if(gRandomNumber > fAcceptanceParameterization->Eval(track->Pt())) 
	  continue;
      }
      
      //Exclude resonances
      if(fExcludeResonancesInMC) {
	TParticle *particle = track->Particle();
	if(!particle) continue;
	
	Bool_t kExcludeParticle = kFALSE;
	Int_t gMotherIndex = particle->GetFirstMother();
	if(gMotherIndex != -1) {
	  AliMCParticle* motherTrack = dynamic_cast<AliMCParticle *>(event->GetTrack(gMotherIndex));
	  if(motherTrack) {
	    TParticle *motherParticle = motherTrack->Particle();
	    if(motherParticle) {
	      Int_t pdgCodeOfMother = motherParticle->GetPdgCode();
	      //if((pdgCodeOfMother == 113)||(pdgCodeOfMother == 213)||(pdgCodeOfMother == 221)||(pdgCodeOfMother == 223)||(pdgCodeOfMother == 331)||(pdgCodeOfMother == 333)) {
	      if(pdgCodeOfMother == 113) {
		kExcludeParticle = kTRUE;
	      }
	    }
	  }
	}
	
	//Exclude from the analysis decay products of rho0, rho+, eta, eta' and phi
	if(kExcludeParticle) continue;
      }
      
      vCharge = track->Charge();
      vPhi    = track->Phi();
      //Printf("phi (before): %lf",vPhi);
      
      fHistPt->Fill(vPt,fCentrality);
      fHistEta->Fill(vEta,fCentrality);
      fHistPhi->Fill(vPhi,fCentrality);
      //fHistPhi->Fill(vPhi*TMath::RadToDeg(),fCentrality);
      fHistRapidity->Fill(vY,fCentrality);
      //if(vCharge > 0) fHistPhiPos->Fill(vPhi*TMath::RadToDeg(),fCentrality);
      //else if(vCharge < 0) fHistPhiNeg->Fill(vPhi*TMath::RadToDeg(),fCentrality);
      if(vCharge > 0) fHistPhiPos->Fill(vPhi,fCentrality);
      else if(vCharge < 0) fHistPhiNeg->Fill(vPhi,fCentrality);
      
      //Flow after burner
      if(fUseFlowAfterBurner) {
	Double_t precisionPhi = 0.001;
	Int_t maxNumberOfIterations = 100;
	
	Double_t phi0 = vPhi;
	Double_t gV2 = fDifferentialV2->Eval(vPt);
	
	for (Int_t j = 0; j < maxNumberOfIterations; j++) {
	  Double_t phiprev = vPhi;
	  Double_t fl = vPhi - phi0 + gV2*TMath::Sin(2.*(vPhi - gReactionPlane*TMath::DegToRad()));
	  Double_t fp = 1.0 + 2.0*gV2*TMath::Cos(2.*(vPhi - gReactionPlane*TMath::DegToRad())); 
	  vPhi -= fl/fp;
	  if (TMath::AreEqualAbs(phiprev,vPhi,precisionPhi)) break;
	}
	//Printf("phi (after): %lf\n",vPhi);
	Double_t vDeltaphiBefore = phi0 - gReactionPlane*TMath::DegToRad();
	if(vDeltaphiBefore < 0) vDeltaphiBefore += 2*TMath::Pi();
	fHistPhiBefore->Fill(vDeltaphiBefore,fCentrality);
	
	Double_t vDeltaphiAfter = vPhi - gReactionPlane*TMath::DegToRad();
	if(vDeltaphiAfter < 0) vDeltaphiAfter += 2*TMath::Pi();
	fHistPhiAfter->Fill(vDeltaphiAfter,fCentrality);
      }
      
      //vPhi *= TMath::RadToDeg();
                
      tracksAccepted->Add(new AliBFBasicParticle(vEta, vPhi, vPt, vCharge));
      
    } //track loop
  }//MC
  
  return tracksAccepted;  
}
//________________________________________________________________________
TObjArray* AliAnalysisTaskBFPsi::GetShuffledTracks(TObjArray *tracks){
  // Clones TObjArray and returns it with tracks after shuffling the charges

  TObjArray* tracksShuffled = new TObjArray;
  tracksShuffled->SetOwner(kTRUE);

  vector<Short_t> *chargeVector = new vector<Short_t>;   //original charge of accepted tracks 

  for (Int_t i=0; i<tracks->GetEntriesFast(); i++)
  {
    AliVParticle* track = (AliVParticle*) tracks->At(i);
    chargeVector->push_back(track->Charge());
  }  
 
  random_shuffle(chargeVector->begin(), chargeVector->end());
  
  for(Int_t i = 0; i < tracks->GetEntriesFast(); i++){
    AliVParticle* track = (AliVParticle*) tracks->At(i);
    tracksShuffled->Add(new AliBFBasicParticle(track->Eta(), track->Phi(), track->Pt(),chargeVector->at(i)));
  }

  delete chargeVector;
   
  return tracksShuffled;
}


//________________________________________________________________________
void  AliAnalysisTaskBFPsi::FinishTaskOutput(){
  //Printf("END BF");

  if (!fBalance) {
    AliError("fBalance not available");
    return;
  }  
  if(fRunShuffling) {
    if (!fShuffledBalance) {
      AliError("fShuffledBalance not available");
      return;
    }
  }

}

//________________________________________________________________________
void AliAnalysisTaskBFPsi::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  // not implemented ...

}
