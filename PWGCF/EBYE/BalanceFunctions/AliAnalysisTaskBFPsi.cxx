#include "TChain.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h" 
#include "TH2D.h"                  
#include "TH3D.h"
#include "TArrayF.h"
#include "TF1.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TDatabasePDG.h"


#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h" 
#include "AliCollisionGeometry.h"
#include "AliGenEventHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliEventplane.h"
#include "AliTHn.h"    
#include "AliLog.h"
#include "AliAnalysisUtils.h"

#include "AliEventPoolManager.h"           

#include "AliPID.h"                
#include "AliPIDResponse.h"        
#include "AliPIDCombined.h"        

#include "AliAnalysisTaskBFPsi.h"
#include "AliBalanceEbyE.h"
#include "AliBalancePsi.h"
#include "AliAnalysisTaskTriggeredBF.h"
#include "TFile.h"
#include <iostream>


// Analysis task for the BF vs Psi code
// Authors: Panos.Christakoglou@nikhef.nl

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskBFPsi)

//________________________________________________________________________
AliAnalysisTaskBFPsi::AliAnalysisTaskBFPsi(const char *name) 
: AliAnalysisTaskSE(name),
  fDebugLevel(kFALSE),
  fArrayMC(0),
  fBalance(0),
  fBalanceEbyE(0),
  fRunShuffling(kFALSE),
  fShuffledBalance(0),
  fRunMixing(kFALSE),
  fRunMixingEventPlane(kFALSE),
  fRunEbyE(kFALSE),
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
  fHistCentStatsUsed(0),
  fHistTriggerStats(0),
  fHistTrackStats(0),
  fHistVx(0),
  fHistVy(0),
  fHistVz(0),
  fHistMixEvents(0),
  fHistMixTracks(0),
  fHistTPCvsVZEROMultiplicity(0),
  fHistCL1vsVZEROPercentile(0),
  fHistVZEROSignal(0),
  fHistEventPlane(0),
  fHistClus(0),
  fHistDCA(0),
  fHistChi2(0),
  fHistPt(0),
  fHistEta(0),
  fHistRapidity(0),
  fHistPhi(0),
  fHistEtaPhiPos(0), 	       	 
  fHistEtaPhiNeg(0), 
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
  fHistBetaVsdEdXbeforePID(NULL),
  fHistNSigmaTPCTOFvsPtbeforePID(NULL),
  fHistNSigmaTPCTOFPbefPID(NULL),
  fHistdEdxVsPTPCafterPID(NULL),
  fHistBetavsPTOFafterPID(NULL), 
  fHistProbTPCvsPtafterPID(NULL), 
  fHistProbTOFvsPtafterPID(NULL), 
  fHistProbTPCTOFvsPtafterPID(NULL),
  fHistNSigmaTPCvsPtafterPID(NULL), 
  fHistNSigmaTOFvsPtafterPID(NULL),  
  fHistBetaVsdEdXafterPID(NULL), 
  fHistNSigmaTPCTOFvsPtafterPID(NULL),
  fHistNSigmaTPCTOFPafterPID(NULL),
  fHistdEdxVsPTPCbeforePIDelectron(NULL),
  fHistNSigmaTPCvsPtbeforePIDelectron(NULL),
  fHistdEdxVsPTPCafterPIDelectron(NULL),
  fHistNSigmaTPCvsPtafterPIDelectron(NULL),
  fCentralityArrayBinsForCorrections(kCENTRALITY),
  fCentralityWeights(0x0),
  fPIDResponse(0x0),
  fPIDCombined(0x0),
  fParticleOfInterest(kPion),
  fPidDetectorConfig(kTPCTOF),
  fMassParticleOfInterest(0.13957),
  fUsePID(kFALSE),
  fUsePIDnSigma(kTRUE),
  fUsePIDPropabilities(kFALSE), 
  fUseRapidity(kFALSE),
  fPIDNSigma(3.),
  fMinAcceptedPIDProbability(0.8),
  fElectronRejection(kFALSE),
  fElectronOnlyRejection(kFALSE),
  fElectronRejectionNSigma(-1.),
  fElectronRejectionMinPt(0.),
  fElectronRejectionMaxPt(1000.),
  fESDtrackCuts(0),
  fCentralityEstimator("V0M"),
  fUseCentrality(kFALSE),
  fUseMultSelectionFramework(kFALSE),
  fUseUncheckedCentrality(kFALSE),
  fCentralityPercentileMin(0.), 
  fCentralityPercentileMax(5.),
  fImpactParameterMin(0.),
  fImpactParameterMax(20.),
  fMultiplicityEstimator("V0A"),
  fUseMultiplicity(kFALSE),
  fNumberOfAcceptedTracksMin(0),
  fNumberOfAcceptedTracksMax(10000),
  fHistNumberOfAcceptedTracks(0),
  fHistMultiplicity(0),
  fHistMultvsPercent(0),
  fUseOfflineTrigger(kFALSE),
  fCheckFirstEventInChunk(kFALSE),
  fCheckPileUp(kFALSE),
  fCheckPrimaryFlagAOD(kFALSE),
  fUseMCforKinematics(kFALSE),
  fUseAdditionalVtxCuts(kFALSE),
  fUseOutOfBunchPileUpCutsLHC15o(kFALSE),
  fVxMax(0.8),
  fVyMax(0.8),
  fVzMax(10.),
  fnAODtrackCutBit(128),
  fPtMin(0.3),
  fPtMax(1.5),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fPhiMin(0.),
  fPhiMax(360.),
  fDCAxyCut(-1),
  fDCAzCut(-1),
  fTPCchi2Cut(-1),
  fNClustersTPCCut(-1),
  fTPCsharedCut(-1),
  fAcceptanceParameterization(0),
  fDifferentialV2(0),
  fUseFlowAfterBurner(kFALSE),
  fIncludeSecondariesInMCgen(kFALSE),
  fExcludeSecondariesInMC(kFALSE),
  fExcludeWeakDecaysInMC(kFALSE),
  fExcludeResonancesInMC(kFALSE),
  fExcludeElectronsInMC(kFALSE),
  fExcludeParticlesExtra(kFALSE),
  fUseMCPdgCode(kFALSE),
  fPDGCodeToBeAnalyzed(-1),
  fExcludeResonancePDGInMC(-1),
  fIncludeResonancePDGInMC(-1),
  fEventClass("EventPlane"), 
  fCustomBinning(""),
  fHistVZEROAGainEqualizationMap(0),
  fHistVZEROCGainEqualizationMap(0),
  fHistVZEROChannelGainEqualizationMap(0),
  fUtils(0) {
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain

  //======================================================correction
  for (Int_t i=0; i<kCENTRALITY; i++){
    fHistCorrectionPlus[i] = NULL; 
    fHistCorrectionMinus[i] = NULL; 
    fCentralityArrayForCorrections[i] = -1.;
  }
  //=====================================================correction

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
  // delete fHistEtaPhiPos; 		 	 
  // delete fHistEtaPhiNeg;
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
  if(fRunEbyE){
    if(!fBalanceEbyE) {
      fBalanceEbyE = new AliBalanceEbyE();
      fBalanceEbyE->SetAnalysisLevel("AOD");
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
  if(fUsePID || fElectronRejection) {
    fHistListPIDQA = new TList();
    fHistListPIDQA->SetName("listQAPID");
    fHistListPIDQA->SetOwner();
  }

  //Event stats.
  TString gCutName[10] = {"Total","Offline trigger",
			  "Vertex","Analyzed","sel. Centrality","Not1stEvInChunk","No Pile-Up", "Add Vtx Cuts", "Rej OOB pile up LHC15o", "Rej TPC vs global LHC15o"};
  fHistEventStats = new TH2F("fHistEventStats",
                             "Event statistics;;Centrality percentile;N_{events}",
                             10,0.5,10.5,220,-5,105);
  for(Int_t i = 1; i <= 10; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fList->Add(fHistEventStats);

  TString gCentName[13] = {"V0M","V0A","V0C","FMD","TRK","TKL","CL0","CL1","ZNA","ZPA","V0MvsFMD","TKLvsV0M","ZEMvsZDC"};
  fHistCentStats = new TH2F("fHistCentStats",
                             "Centrality statistics;;Cent percentile",
			    13,-0.5,12.5,220,-5,105);
  for(Int_t i = 1; i <= 13; i++){
    fHistCentStats->GetXaxis()->SetBinLabel(i,gCentName[i-1].Data());
    //fHistCentStatsUsed->GetXaxis()->SetBinLabel(i,gCentName[i-1].Data());
  }
  fList->Add(fHistCentStats);

  fHistCentStatsUsed = new TH2F("fHistCentStatsUsed","Centrality statistics;;Cent percentile", 1,-0.5,0.5,220,-5,105);
  fHistCentStatsUsed->GetXaxis()->SetBinLabel(1,fCentralityEstimator.Data());
  fList->Add(fHistCentStatsUsed);

  fHistTriggerStats = new TH1F("fHistTriggerStats","Trigger statistics;TriggerBit;N_{events}",1025,0,1025);
  fList->Add(fHistTriggerStats);

  fHistTrackStats = new TH1F("fHistTrackStats","Event statistics;TrackFilterBit;N_{events}",16,0,16);
  fList->Add(fHistTrackStats);

  fHistNumberOfAcceptedTracks = new TH2F("fHistNumberOfAcceptedTracks",";N_{acc.};Centrality percentile;Entries",4001,-0.5,4000.5,220,-5,105);
  fList->Add(fHistNumberOfAcceptedTracks);

  fHistMultiplicity = new TH1F("fHistMultiplicity",";N_{ch.};Entries",30001,-0.5,30000.5);
  fList->Add(fHistMultiplicity);

  fHistMultvsPercent = new TH2F("fHistMultvsPercent",";N_{ch.};Centrality percentile;Entries",30001,-0.5,30000.5, 220,-5,105);
  fList->Add(fHistMultvsPercent);

  // Vertex distributions
  fHistVx = new TH1F("fHistVx","Primary vertex distribution - x coordinate;V_{x} (cm);Entries",100,-0.5,0.5);
  fList->Add(fHistVx);
  fHistVy = new TH1F("fHistVy","Primary vertex distribution - y coordinate;V_{y} (cm);Entries",100,-0.5,0.5);
  fList->Add(fHistVy);
  fHistVz = new TH2F("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm);Centrality percentile;Entries",100,-20.,20.,220,-5,105);
  fList->Add(fHistVz);

  // Event Mixing
  fHistMixEvents = new TH2F("fHistMixEvents","Number of mixed events;Centrality percentile;N_{mix,evts}",101, 0, 101, 200, 0, 200);
  fList->Add(fHistMixEvents);
  fHistMixTracks = new TH2F("fHistMixTracks","Number of mixed tracks;Centrality percentile;N_{mix,trks}",101, 0, 101, 200, 0, fMixingTracks * 1.5);
  fList->Add(fHistMixTracks);

  //TPC vs VZERO multiplicity
  fHistTPCvsVZEROMultiplicity = new TH2F("fHistTPCvsVZEROMultiplicity","VZERO vs TPC multiplicity",10001,-0.5,10000.5,4001,-0.5,4000.5);
  if(fMultiplicityEstimator == "V0A") 
    fHistTPCvsVZEROMultiplicity->GetXaxis()->SetTitle("VZERO-A multiplicity (a.u.)");
  else if(fMultiplicityEstimator == "V0C") 
    fHistTPCvsVZEROMultiplicity->GetXaxis()->SetTitle("VZERO-C multiplicity (a.u.)");
  else 
    fHistTPCvsVZEROMultiplicity->GetXaxis()->SetTitle("VZERO multiplicity (a.u.)");
  fList->Add(fHistTPCvsVZEROMultiplicity);

  fHistCL1vsVZEROPercentile = new TH2F("fHistCL1vsVZEROPercentile", "V0M vs CL1 centrality percentile", 101, 0, 101, 101, 0, 101);
  fList->Add(fHistCL1vsVZEROPercentile);

  fHistVZEROSignal = new TH2F("fHistVZEROSignal","VZERO signal vs VZERO channel;VZERO channel; Signal (a.u.)",64,0.5,64.5,3001,-0.5,30000.5);
  fList->Add(fHistVZEROSignal);

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
  fHistEtaPhiPos  = new TH3F("fHistEtaPhiPos","#eta-#phi distribution (+);#eta;#phi (rad);Centrality percentile",40,-1.6,1.6,72,0.,2.*TMath::Pi(),220,-5,105); 		 	 
  fList->Add(fHistEtaPhiPos); 			 
  fHistEtaPhiNeg  = new TH3F("fHistEtaPhiNeg","#eta-#phi distribution (-);#eta;#phi (rad);Centrality percentile",40,-1.6,1.6,72,0.,2.*TMath::Pi(),220,-5,105); 	       	 
  fList->Add(fHistEtaPhiNeg);
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

  // Balance function histograms
  // Initialize histograms if not done yet (including the custom binning)
  if(!fBalance->GetHistNp()){
    AliInfo("Histograms not yet initialized! --> Will be done now");
    fBalance->SetCustomBinning(fCustomBinning);
    fBalance->InitHistograms();
  }

  if(fRunShuffling) {
    if(!fShuffledBalance->GetHistNp()) {
      AliInfo("Histograms (shuffling) not yet initialized! --> Will be done now");
      fShuffledBalance->SetCustomBinning(fCustomBinning);
      fShuffledBalance->InitHistograms();
    }
  }

  if(fRunMixing) {
    if(!fMixedBalance->GetHistNp()) {
      AliInfo("Histograms (mixing) not yet initialized! --> Will be done now");
      fMixedBalance->SetCustomBinning(fCustomBinning);
      fMixedBalance->InitHistograms();
    }
  }

  if(fRunEbyE) {
    if(!fBalanceEbyE->GetHistBF()) {
      AliInfo("Histograms (EbyE) not yet initialized! --> Will be done now");
      fBalanceEbyE->InitHistograms();
    }
    fList->Add(fBalanceEbyE->GetHistBF());
  }

  // QA histograms for different cuts
  fList->Add(fBalance->GetQAHistHBTbefore());
  fList->Add(fBalance->GetQAHistHBTafter());
  fList->Add(fBalance->GetQAHistPhiStarHBTbefore());
  fList->Add(fBalance->GetQAHistPhiStarHBTafter());
  fList->Add(fBalance->GetQAHistConversionbefore());
  fList->Add(fBalance->GetQAHistConversionafter());
  fList->Add(fBalance->GetQAHistPsiMinusPhi());
  fList->Add(fBalance->GetQAHistResonancesBefore());
  fList->Add(fBalance->GetQAHistResonancesRho());
  fList->Add(fBalance->GetQAHistResonancesK0());
  fList->Add(fBalance->GetQAHistResonancesLambda());
  fList->Add(fBalance->GetQAHistQbefore());
  fList->Add(fBalance->GetQAHistQafter());

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
    Double_t* centbins = NULL;
    Int_t nCentralityBins;
    if(fBalance->IsUseVertexBinning()){
      centbins = fBalance->GetBinning(fBalance->GetBinningString(), "centralityVertex", nCentralityBins);
    }
    else{
      centbins = fBalance->GetBinning(fBalance->GetBinningString(), "centrality", nCentralityBins);
    }
    
    // multiplicity bins
    Double_t* multbins = NULL;
    Int_t nMultiplicityBins;
    multbins = fBalance->GetBinning(fBalance->GetBinningString(), "multiplicity", nMultiplicityBins);
    
    // Zvtx bins
    Double_t* vtxbins = NULL; 
    Int_t nVertexBins;
    if(fBalance->IsUseVertexBinning()){
      vtxbins = fBalance->GetBinning(fBalance->GetBinningString(), "vertexVertex", nVertexBins);
    }
    else{
      vtxbins = fBalance->GetBinning(fBalance->GetBinningString(), "vertex", nVertexBins);
    }

    // Event plane angle (Psi) bins
    Double_t* psibins = NULL;
    Int_t nPsiBins; 
    psibins = fBalance->GetBinning(fBalance->GetBinningString(), "eventPlane", nPsiBins);

  
    // run the event mixing also in bins of event plane (statistics!)
    if(fRunMixingEventPlane){
      if(fEventClass=="Multiplicity"){
	if(multbins && vtxbins && psibins){
	  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nMultiplicityBins, multbins, nVertexBins, vtxbins, nPsiBins, psibins);
	}
      }
      else{
	if(centbins && vtxbins && psibins){
	  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centbins, nVertexBins, vtxbins, nPsiBins, psibins);
	}
      }
    }
    else{
      if(fEventClass=="Multiplicity"){
	if(multbins && vtxbins){
	  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nMultiplicityBins, multbins, nVertexBins, vtxbins);
	}
      }
      else{
	if(centbins && vtxbins){
	  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centbins, nVertexBins, vtxbins);
	}
      }
    }
    
    if(centbins) delete [] centbins; 
    if(multbins) delete [] multbins; 
    if(vtxbins)  delete [] vtxbins; 
    if(psibins)  delete [] psibins; 

    // set minimum values for track depth, fraction, and number of events
    fPoolMgr->SetTargetValues(fMixingTracks, 0.1, 5);
    
    // check pool manager
    if(!fPoolMgr){
      AliError("Event Mixing required, but Pool Manager not initialized...");
      return;
    }
  }
  
  if(fESDtrackCuts) fList->Add(fESDtrackCuts);

  //====================PID========================//
  if(fUsePID) {
    fPIDCombined = new AliPIDCombined();
    fPIDCombined->SetDefaultTPCPriors();

    fHistdEdxVsPTPCbeforePID = new TH2D ("dEdxVsPTPCbefore","dEdxVsPTPCbefore", 1000, -10.0, 10.0, 1000, 0, 1000); 
    fHistListPIDQA->Add(fHistdEdxVsPTPCbeforePID);
    
    fHistBetavsPTOFbeforePID = new TH2D ("BetavsPTOFbefore","BetavsPTOFbefore", 1000, -10.0, 10., 1000, 0, 1.2); 
    fHistListPIDQA->Add(fHistBetavsPTOFbeforePID); 
    
    fHistProbTPCvsPtbeforePID = new TH2D ("ProbTPCvsPtbefore","ProbTPCvsPtbefore", 1000, -10.0,10.0, 1000, 0, 2.0); 
    fHistListPIDQA->Add(fHistProbTPCvsPtbeforePID); 
    
    fHistProbTOFvsPtbeforePID = new TH2D ("ProbTOFvsPtbefore","ProbTOFvsPtbefore", 1000, -50, 50, 1000, 0, 2.0); 
    fHistListPIDQA->Add(fHistProbTOFvsPtbeforePID);

    fHistProbTPCTOFvsPtbeforePID =new TH2D ("ProbTPCTOFvsPtbefore","ProbTPCTOFvsPtbefore", 1000, -50, 50, 1000, 0, 2.0); 
    fHistListPIDQA->Add(fHistProbTPCTOFvsPtbeforePID);
    
    fHistNSigmaTPCvsPtbeforePID = new TH2D ("NSigmaTPCvsPtbefore","NSigmaTPCvsPtbefore", 1000, -10, 10, 1000, -25, 25); 
    fHistListPIDQA->Add(fHistNSigmaTPCvsPtbeforePID);
    
    fHistNSigmaTOFvsPtbeforePID = new TH2D ("NSigmaTOFvsPtbefore","NSigmaTOFvsPtbefore", 1000, -10, 10, 1000, -25, 25); 
    fHistListPIDQA->Add(fHistNSigmaTOFvsPtbeforePID); 

    fHistBetaVsdEdXbeforePID = new TH2D ("BetaVsdEdXbefore","BetaVsdEdXbefore", 1000, 0., 1000, 1000, 0, 1.2); 
    fHistListPIDQA->Add(fHistBetaVsdEdXbeforePID);
    
    fHistNSigmaTPCTOFvsPtbeforePID = new TH2D ("NSigmaTPCTOFvsPtbefore","NSigmaTPCTOFvsPtbefore", 1000, -10., 10., 1000, -25, 25); 
    fHistListPIDQA->Add(fHistNSigmaTPCTOFvsPtbeforePID);
    
    //+++++++++++++++++//
    //p array
    const Int_t pBins = 36;
    Double_t nArrayP[pBins+1]={0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0};
    //nSigma Array
    const Int_t nSigmaBins = 250;
    Double_t nArrayS[nSigmaBins+1];
    for (Int_t i = 0; i <= nSigmaBins; i++){
      nArrayS[i]=i-125; //i+1
      //Printf("nS: %lf - i: %d", nSigmaArray[i], i);
    }
 
    fHistNSigmaTPCTOFPbefPID = new TH3D ("fHistNSigmaTPCTOFPbefPID","fHistNSigmaTPCTOFPbefPID;#sigma_{TPC};#sigma_{TOF};p_{T} (GeV/c)", nSigmaBins, nArrayS, nSigmaBins, nArrayS, pBins,nArrayP); 
    fHistListPIDQA->Add(fHistNSigmaTPCTOFPbefPID); 
    //++++++++++++++//

    fHistdEdxVsPTPCafterPID = new TH2D ("dEdxVsPTPCafter","dEdxVsPTPCafter", 1000, -10, 10, 1000, 0, 1000); 
    fHistListPIDQA->Add(fHistdEdxVsPTPCafterPID);
    
    fHistBetavsPTOFafterPID = new TH2D ("BetavsPTOFafter","BetavsPTOFafter", 1000, -10, 10, 1000, 0, 1.2); 
    fHistListPIDQA->Add(fHistBetavsPTOFafterPID); 
    
    fHistProbTPCvsPtafterPID = new TH2D ("ProbTPCvsPtafter","ProbTPCvsPtafter", 1000, -10, 10, 1000, 0, 2); 
    fHistListPIDQA->Add(fHistProbTPCvsPtafterPID);
  
    fHistProbTOFvsPtafterPID = new TH2D ("ProbTOFvsPtafter","ProbTOFvsPtafter", 1000,  -10, 10, 1000, 0, 2); 
    fHistListPIDQA->Add(fHistProbTOFvsPtafterPID); 
    
    fHistProbTPCTOFvsPtafterPID =new TH2D ("ProbTPCTOFvsPtafter","ProbTPCTOFvsPtafter", 1000, -50, 50, 1000, 0, 2.0); 
    fHistListPIDQA->Add(fHistProbTPCTOFvsPtafterPID);

    fHistNSigmaTPCvsPtafterPID = new TH2D ("NSigmaTPCvsPtafter","NSigmaTPCvsPtafter", 1000, -10, 10, 1000, -25, 25); 
    fHistListPIDQA->Add(fHistNSigmaTPCvsPtafterPID);
    
    fHistNSigmaTOFvsPtafterPID = new TH2D ("NSigmaTOFvsPtafter","NSigmaTOFvsPtafter", 1000, -10, 10, 1000, -25, 25); 
    fHistListPIDQA->Add(fHistNSigmaTOFvsPtafterPID);

    fHistBetaVsdEdXafterPID = new TH2D ("BetaVsdEdXafter","BetaVsdEdXafter", 1000, 0., 1000, 1000, 0, 1.2); 
    fHistListPIDQA->Add(fHistBetaVsdEdXafterPID);

    fHistNSigmaTPCTOFvsPtafterPID = new TH2D ("NSigmaTPCTOFvsPtafter","NSigmaTPCTOFvsPtafter", 1000, -10., 10., 1000, -25, 25); 
    fHistListPIDQA->Add(fHistNSigmaTPCTOFvsPtafterPID);

    fHistNSigmaTPCTOFPafterPID = new TH3D ("fHistNSigmaTPCTOFPafterPID","fHistNSigmaTPCTOFPafterPID;#sigma_{TPC};#sigma_{TOF};p_{T} (GeV/c)", nSigmaBins, nArrayS, nSigmaBins, nArrayS, pBins,nArrayP); 
    fHistListPIDQA->Add(fHistNSigmaTPCTOFPafterPID); //++++++++++++++
  }

  // for electron rejection only TPC nsigma histograms
  if(fElectronRejection) {
 
    fHistdEdxVsPTPCbeforePIDelectron = new TH2D ("dEdxVsPTPCbeforeelectron","dEdxVsPTPCbeforeelectron", 1000, -10.0, 10.0, 1000, 0, 1000); 
    fHistListPIDQA->Add(fHistdEdxVsPTPCbeforePIDelectron);
    
    fHistNSigmaTPCvsPtbeforePIDelectron = new TH2D ("NSigmaTPCvsPtbeforeelectron","NSigmaTPCvsPtbeforeelectron", 1000, -10, 10, 1000, 0, 500); 
    fHistListPIDQA->Add(fHistNSigmaTPCvsPtbeforePIDelectron);
    
    fHistdEdxVsPTPCafterPIDelectron = new TH2D ("dEdxVsPTPCafterelectron","dEdxVsPTPCafterelectron", 1000, -10, 10, 1000, 0, 1000); 
    fHistListPIDQA->Add(fHistdEdxVsPTPCafterPIDelectron);

    fHistNSigmaTPCvsPtafterPIDelectron = new TH2D ("NSigmaTPCvsPtafterelectron","NSigmaTPCvsPtafterelectron", 1000, -10, 10, 1000, 0, 500); 
    fHistListPIDQA->Add(fHistNSigmaTPCvsPtafterPIDelectron); 
  }
  //====================PID========================//

  // Post output data.
  PostData(1, fList);
  PostData(2, fListBF);
  if(fRunShuffling) PostData(3, fListBFS);
  if(fRunMixing) PostData(4, fListBFM);
  if(fUsePID || fElectronRejection) PostData(5, fHistListPIDQA);       //PID

  AliInfo("Finished setting up the Output");

  TH1::AddDirectory(oldStatus);
  
  fUtils = new AliAnalysisUtils();
}


//________________________________________________________________________
void AliAnalysisTaskBFPsi::SetInputCorrection(TString filename, 
					      Int_t nCentralityBins, 
					      Double_t *centralityArrayForCorrections) {
  //Open files that will be used for correction
  fCentralityArrayBinsForCorrections = nCentralityBins;
  for (Int_t i=0; i<nCentralityBins; i++)
    fCentralityArrayForCorrections[i] = centralityArrayForCorrections[i];

  // No file specified -> Abort
  if(!filename.Contains(".root")) {
    AliFatal(Form("No correction file specified (= %s) but correction requested ==> ABORT",filename.Data()));
    return;
  }

  //Open the input file
  TFile *f = TFile::Open(filename);
  if(!f->IsOpen()) {
    AliFatal(Form("File %s not found but correction requested ==> ABORT",filename.Data()));
    return;
  }
    
  //TString listEffName = "";
  for (Int_t iCent = 0; iCent < fCentralityArrayBinsForCorrections-1; iCent++) {    
    //Printf("iCent %d:",iCent);    
    TString histoName = "fHistCorrectionPlus";
    histoName += Form("%d-%d",(Int_t)(fCentralityArrayForCorrections[iCent]),(Int_t)(fCentralityArrayForCorrections[iCent+1]));
    fHistCorrectionPlus[iCent]= dynamic_cast<TH3F *>(f->Get(histoName.Data()));
    if(!fHistCorrectionPlus[iCent]) {
      AliFatal(Form("fHist %s not found but correction requested ==> ABORT",histoName.Data()));
      return;
    }
    
    histoName = "fHistCorrectionMinus";
    histoName += Form("%d-%d",(Int_t)(fCentralityArrayForCorrections[iCent]),(Int_t)(fCentralityArrayForCorrections[iCent+1]));
    fHistCorrectionMinus[iCent] = dynamic_cast<TH3F *>(f->Get(histoName.Data())); 
    if(!fHistCorrectionMinus[iCent]) {
      AliFatal(Form("fHist %s not found but correction requested ==> ABORT",histoName.Data()));
      return; 
    }
  }//loop over centralities: ONLY the PbPb case is covered
}

//________________________________________________________________________
void AliAnalysisTaskBFPsi::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  TString gAnalysisLevel = fBalance->GetAnalysisLevel();
  Int_t gNumberOfAcceptedTracks = 0;
  Double_t lMultiplicityVar     = -999.; //-1
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
  if(fUsePID || fElectronRejection) {
    fPIDResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
    if (!fPIDResponse) AliFatal("This Task needs the PID response attached to the inputHandler");
  }
 
  
  // check event cuts and fill event histograms
  if((lMultiplicityVar = IsEventAccepted(eventMain)) < 0){ 
    return;
  }
  // get the reaction plane
  if(fEventClass != "Multiplicity" && gAnalysisLevel!="AODnano") {
    gReactionPlane = GetEventPlane(eventMain);
    fHistEventPlane->Fill(gReactionPlane,lMultiplicityVar);
    if(gReactionPlane < 0){
      return;
    }
  }

  // get the accepted tracks in main event
  TObjArray *tracksMain = GetAcceptedTracks(eventMain,lMultiplicityVar,gReactionPlane);
  gNumberOfAcceptedTracks = tracksMain->GetEntriesFast();
  
  //multiplicity cut (used in pp)
  fHistNumberOfAcceptedTracks->Fill(gNumberOfAcceptedTracks,lMultiplicityVar);

  // store charges of all accepted tracks,shuffle and reassign(two extra loops)
  TObjArray* tracksShuffled = NULL;
  if(fRunShuffling){
    tracksShuffled = GetShuffledTracks(tracksMain,lMultiplicityVar);
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
      
      AliEventPool* pool = fPoolMgr->GetEventPool(lMultiplicityVar, eventMain->GetPrimaryVertex()->GetZ(),gReactionPlane);
      
      if (!pool){
	AliFatal(Form("No pool found for centrality = %f, zVtx = %f, psi = %f", lMultiplicityVar, eventMain->GetPrimaryVertex()->GetZ(),gReactionPlane));
      }
      else{
	
	//pool->SetDebug(1);

	if (pool->IsReady()){ 
	  
	  
	  Int_t nMix = pool->GetCurrentNEvents();
	  //cout << "nMix = " << nMix << " tracks in pool = " << pool->NTracksInPool() << endl;
	  
	  fHistMixEvents->Fill(lMultiplicityVar, nMix);
	  fHistMixTracks->Fill(lMultiplicityVar, pool->NTracksInPool());

	  // Fill mixed-event histos here  
	  for (Int_t jMix=0; jMix<nMix; jMix++) 
	    {
	      TObjArray* tracksMixed = pool->GetEvent(jMix);
	      fMixedBalance->CalculateBalance(gReactionPlane,tracksMain,tracksMixed,bSign,lMultiplicityVar,eventMain->GetPrimaryVertex()->GetZ());
	    }
	}
	
	// Update the Event pool
	pool->UpdatePool(tracksMain);
	//pool->PrintInfo();
	
      }//pool NULL check  
    }//run mixing
  
  // calculate balance function
  fBalance->CalculateBalance(gReactionPlane,tracksMain,NULL,bSign,lMultiplicityVar,eventMain->GetPrimaryVertex()->GetZ());
  
  // calculate shuffled balance function
  if(fRunShuffling && tracksShuffled != NULL) {
    fShuffledBalance->CalculateBalance(gReactionPlane,tracksShuffled,NULL,bSign,lMultiplicityVar,eventMain->GetPrimaryVertex()->GetZ());
  }

  // calculate balance function on an event-by-event basis
  if(fRunEbyE){
    fBalanceEbyE->CalculateBalance(gReactionPlane,tracksMain,NULL,bSign,lMultiplicityVar,eventMain->GetPrimaryVertex()->GetZ());
  }      
}

//________________________________________________________________________
Double_t AliAnalysisTaskBFPsi::IsEventAccepted(AliVEvent *event){
  // Checks the Event cuts
  // Fills Event statistics histograms
  
  Bool_t isSelectedMain = kTRUE;
  Float_t gRefMultiplicity = -1.;
  TString gAnalysisLevel = fBalance->GetAnalysisLevel();

  AliMCEvent *mcevent = dynamic_cast<AliMCEvent*>(event);

  fHistEventStats->Fill(1,gRefMultiplicity); //all events

  // check first event in chunk (is not needed for new reconstructions)
  if(fCheckFirstEventInChunk){
    if(fUtils->IsFirstEventInChunk(event)) 
      return -1.;
    fHistEventStats->Fill(6,gRefMultiplicity); 
  }
  // check for pile-up event
  AliAnalysisUtils ut;
  if(fCheckPileUp){
    fUtils->SetUseMVPlpSelection(kTRUE);
    fUtils->SetUseOutOfBunchPileUp(kTRUE);
    if(fUtils->IsPileUpEvent(event))
      return -1.;
    fHistEventStats->Fill(7,gRefMultiplicity); 
  }

  
  // Event trigger bits
  fHistTriggerStats->Fill(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());
  if(fUseOfflineTrigger)
    isSelectedMain = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
  if(isSelectedMain) {
    fHistEventStats->Fill(2,gRefMultiplicity); //triggered events
 
    // Event Vertex MC
    if(gAnalysisLevel == "MC") {
      if(!event) {
	AliError("mcEvent not available");
	return 0x0;
      }
      
      if(mcevent){
	AliGenEventHeader *header = dynamic_cast<AliGenEventHeader*>(mcevent->GenEventHeader());
	if(header){  
	  TArrayF gVertexArray;
	  header->PrimaryVertex(gVertexArray);
	  //Printf("Vertex: %lf (x) - %lf (y) - %lf (z)",
	  //gVertexArray.At(0),
	  //gVertexArray.At(1),
	  //gVertexArray.At(2));
	  fHistEventStats->Fill(3,gRefMultiplicity); //events with a proper vertex
	  if(TMath::Abs(gVertexArray.At(0)) < fVxMax) {
	    if(TMath::Abs(gVertexArray.At(1)) < fVyMax) {
	      if(TMath::Abs(gVertexArray.At(2)) < fVzMax) {
		fHistEventStats->Fill(4,gRefMultiplicity);//analyzed events

		// get the reference multiplicty or centrality
		gRefMultiplicity = GetRefMultiOrCentrality(event);

		fHistVx->Fill(gVertexArray.At(0));
		fHistVy->Fill(gVertexArray.At(1));
		fHistVz->Fill(gVertexArray.At(2),gRefMultiplicity);
		
		// take only events inside centrality class
		if(fUseCentrality) {
		  if((fImpactParameterMin < gRefMultiplicity) && (fImpactParameterMax > gRefMultiplicity)){
		    fHistEventStats->Fill(5,gRefMultiplicity); //events with correct centrality
		    return gRefMultiplicity;	    
		  }//centrality class
		}
		// take events only within the same multiplicity class
		else if(fUseMultiplicity) {
		  if((gRefMultiplicity > fNumberOfAcceptedTracksMin) && (gRefMultiplicity < fNumberOfAcceptedTracksMax)) {
		    fHistEventStats->Fill(5,gRefMultiplicity); //events with correct multiplicity
		    return gRefMultiplicity;
		  }
		}//multiplicity range
	      }//Vz cut
	    }//Vy cut
	  }//Vx cut
	}//header    
      }//MC event object
    }//MC
    
    // Event Vertex AOD, ESD, ESDMC
    else{
      const AliVVertex *vertex = event->GetPrimaryVertex();
      
      if(vertex) {
	Double32_t fCov[6];
	vertex->GetCovarianceMatrix(fCov);
	if(vertex->GetNContributors() > 0) {
	  if(fCov[5] != 0) {
	    fHistEventStats->Fill(3,gRefMultiplicity); //proper vertex
	    if(TMath::Abs(vertex->GetX()) < fVxMax) {
	      if(TMath::Abs(vertex->GetY()) < fVyMax) {
		if(TMath::Abs(vertex->GetZ()) < fVzMax) {
		  fHistEventStats->Fill(4,gRefMultiplicity);//analyzed events

		   if (fUseAdditionalVtxCuts){
		     const AliVVertex *vertex = event->GetPrimaryVertex();
		     if ((!vertex) || (vertex->GetNContributors()<1)){
		       Printf("No Primary Vertex");
		       return -1;
		     }
		     
		     const AliVVertex *vSPD = ((AliAODEvent*)event)->GetPrimaryVertexSPD();
		     if (!vSPD) {
		       Printf("No vertex SPD");
		       return -1;
		     }
		     Double_t dz = vSPD->GetZ()-vertex->GetZ();
		     double covTrc[6],covSPD[6];
		     vertex->GetCovarianceMatrix(covTrc);
		     vSPD->GetCovarianceMatrix(covSPD);
		     double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
		     double errTrc = TMath::Sqrt(covTrc[5]);
		     double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
		     if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20)
		       return -1;
		     
		     fHistEventStats->Fill(8,gRefMultiplicity); 
		   }
		   
		   // get the reference multiplicty or centrality for run1 data
		   if((event->GetRunNumber()<244824) && (fEventClass=="Multiplicity")&&(fMultiplicityEstimator.Contains("Utils"))) {
		     if ((fMultiplicityEstimator == "V0MUtils")) {
		       gRefMultiplicity = fUtils->GetMultiplicityPercentile(event,"V0MEq");
		    if ((fMultiplicityEstimator == "V0AUtils")) 
		      gRefMultiplicity = fUtils->GetMultiplicityPercentile(event,"V0AEq");
		    if ((fMultiplicityEstimator == "V0CUtils")) 
		      gRefMultiplicity = fUtils->GetMultiplicityPercentile(event,"V0CEq");
		    else 
		      AliError("The requested estimator from AliAnalysisUtils is not supported");
		    }//use the framework to define the multiplicity class
		  } 
		  else
		    gRefMultiplicity = GetRefMultiOrCentrality(event);
		  
		  fHistVx->Fill(vertex->GetX());
		  fHistVy->Fill(vertex->GetY());
		  fHistVz->Fill(vertex->GetZ(),gRefMultiplicity);
		  
		  // take only events inside centrality class
		  // if(fUseCentrality) {
		  if((gRefMultiplicity > fCentralityPercentileMin) && (gRefMultiplicity < fCentralityPercentileMax)){
		    
		    // centrality weighting (optional for 2011 if central and semicentral triggers are used)
		    if (fCentralityWeights && !AcceptEventCentralityWeight(gRefMultiplicity)){
		      AliInfo(Form("Rejecting event because of centrality weighting: %f", gRefMultiplicity));
		      return -1;
		    }

		    fHistEventStats->Fill(5,gRefMultiplicity); //events with correct centrality
		    return gRefMultiplicity;		
		  }//centrality class
		  
		  // take events only within the same multiplicity class RUN1! data! 
		  else if((fUseMultiplicity)&&(event->GetRunNumber()<244824)){
		    //if(fDebugLevel) 
		    //Printf("N(min): %.0f, N(max): %.0f - N(ref): %.0f",fNumberOfAcceptedTracksMin,
		    //fNumberOfAcceptedTracksMax,gRefMultiplicity);
		  
		    if((gRefMultiplicity > fNumberOfAcceptedTracksMin) && (gRefMultiplicity < fNumberOfAcceptedTracksMax)) {
		      fHistEventStats->Fill(5,gRefMultiplicity); //events with correct multiplicity
		      return gRefMultiplicity;
		    }
		  } //multiplicity range
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

  Float_t gCentrality = -1.;
  Double_t gMultiplicity = -1.;
  Double_t gMultiplicityFromAOD = -1.;
  TString gAnalysisLevel = fBalance->GetAnalysisLevel();
  
  // use AliMultSelection framework
  //if (fUseMultSelectionFramework) {
  
  AliMultSelection *multSelection = (AliMultSelection*) event->FindListObject("MultSelection");
  if (!multSelection)
    AliFatal("MultSelection not found in input event");
  
  if (fEventClass=="Multiplicity") {
    
    if (fUseUncheckedCentrality)
      gCentrality = multSelection->GetMultiplicityPercentile(fCentralityEstimator, kFALSE);
    else
      gCentrality = multSelection->GetMultiplicityPercentile(fCentralityEstimator, kTRUE);

    // error handling
    if (gCentrality > 100)
      gCentrality = -1;
    
    // QA for centrality estimators (only for checked centrality)
    fHistCentStats->Fill(0.,multSelection->GetMultiplicityPercentile("V0M", kTRUE));
    fHistCentStats->Fill(1.,multSelection->GetMultiplicityPercentile("V0A", kTRUE));
    fHistCentStats->Fill(2.,multSelection->GetMultiplicityPercentile("V0C", kTRUE));
    fHistCentStats->Fill(3.,multSelection->GetMultiplicityPercentile("FMD", kTRUE));
    fHistCentStats->Fill(4.,multSelection->GetMultiplicityPercentile("TRK", kTRUE));
    fHistCentStats->Fill(5.,multSelection->GetMultiplicityPercentile("TKL", kTRUE));
    fHistCentStats->Fill(6.,multSelection->GetMultiplicityPercentile("CL0", kTRUE));
    fHistCentStats->Fill(7.,multSelection->GetMultiplicityPercentile("CL1", kTRUE));
    fHistCentStats->Fill(8.,multSelection->GetMultiplicityPercentile("ZNA", kTRUE));
    fHistCentStats->Fill(9.,multSelection->GetMultiplicityPercentile("ZPA", kTRUE));
    fHistCentStats->Fill(10.,multSelection->GetMultiplicityPercentile("V0MvsFMD", kTRUE));
    fHistCentStats->Fill(11.,multSelection->GetMultiplicityPercentile("TKLvsV0M", kTRUE));
    fHistCentStats->Fill(12.,multSelection->GetMultiplicityPercentile("ZEMvsZDC", kTRUE));
    
    // Centrality estimator USED   ++++++++++++++++++++++++++++++
    fHistCentStatsUsed->Fill(0.,gCentrality);
    
    gMultiplicity = multSelection->GetEstimator(fCentralityEstimator)->GetValue();
    fHistMultiplicity->Fill(gMultiplicity);
    fHistMultvsPercent->Fill(gMultiplicity, gCentrality);

    if (fUseOutOfBunchPileUpCutsLHC15o) {
      if (TMath::Abs(multSelection->GetMultiplicityPercentile("V0M") - multSelection->GetMultiplicityPercentile("CL1")) > 7.5) {
      fHistEventStats->Fill(9, -1);
      return -1;
    }
      const Int_t nTracks = event->GetNumberOfTracks();
      Int_t multEsd = ((AliAODHeader*)event->GetHeader())->GetNumberOfESDTracks();
      Int_t multTPC = 0;
      for (Int_t it = 0; it < nTracks; it++) {
	AliAODTrack* AODTrk = (AliAODTrack*)event->GetTrack(it);
	if (!AODTrk){ delete AODTrk; continue; }
	if (AODTrk->TestFilterBit(128)) {multTPC++;}
      } // end of for (Int_t it = 0; it < nTracks; it++)

      if ((multEsd - 3.38*multTPC) > 15000) return -1;

    }
    
    fHistCL1vsVZEROPercentile->Fill(multSelection->GetMultiplicityPercentile("V0M"),multSelection->GetMultiplicityPercentile("CL1"));
    
    if(multSelection->GetEstimator("RefMult08"))
      fHistTPCvsVZEROMultiplicity->Fill( multSelection->GetEstimator("V0M")->GetValue(),multSelection->GetEstimator("RefMult08")->GetValue());
    else
      gMultiplicityFromAOD = GetReferenceMultiplicityFromAOD(event);
    
  }
  
  // use centrality framework
  else{ 
    
    // calculate centrality always (not only in centrality mode)
    if(gAnalysisLevel == "AOD"|| gAnalysisLevel == "MCAOD" || gAnalysisLevel == "MCAODrec" ) { //centrality in AOD header  //++++++++++++++
      AliAODHeader *header = (AliAODHeader*) event->GetHeader();
      
      if(header){
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
	
	//AOD header
	
	if(event->GetRunNumber()<244824) { //Run1 data. Old centrality framework. 
	  gCentrality = header->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
	  
	  // QA for centrality estimators
	  fHistCentStats->Fill(0.,header->GetCentralityP()->GetCentralityPercentile("V0M"));
	  fHistCentStats->Fill(1.,header->GetCentralityP()->GetCentralityPercentile("V0A"));
	  fHistCentStats->Fill(2.,header->GetCentralityP()->GetCentralityPercentile("V0C"));
	  fHistCentStats->Fill(3.,header->GetCentralityP()->GetCentralityPercentile("FMD"));
	  fHistCentStats->Fill(4.,header->GetCentralityP()->GetCentralityPercentile("TRK"));
	  fHistCentStats->Fill(5.,header->GetCentralityP()->GetCentralityPercentile("TKL")); 
	  fHistCentStats->Fill(6.,header->GetCentralityP()->GetCentralityPercentile("CL0"));
	  fHistCentStats->Fill(7.,header->GetCentralityP()->GetCentralityPercentile("CL1"));
	  fHistCentStats->Fill(8.,header->GetCentralityP()->GetCentralityPercentile("ZNA"));
	  fHistCentStats->Fill(9.,header->GetCentralityP()->GetCentralityPercentile("ZPA"));
	  fHistCentStats->Fill(10.,header->GetCentralityP()->GetCentralityPercentile("V0MvsFMD"));
	  fHistCentStats->Fill(11.,header->GetCentralityP()->GetCentralityPercentile("TKLvsV0M"));
	  fHistCentStats->Fill(12.,header->GetCentralityP()->GetCentralityPercentile("ZEMvsZDC"));
	  
	  // Centrality estimator USED   ++++++++++++++++++++++++++++++
	  fHistCentStatsUsed->Fill(0.,header->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data()));
	}//Run1 if
	
	else { //Run2 data. New multiplicity framework
	  
	  if (fUseUncheckedCentrality)
	    gCentrality = multSelection->GetMultiplicityPercentile(fCentralityEstimator, kFALSE);
	  else
	    gCentrality = multSelection->GetMultiplicityPercentile(fCentralityEstimator, kTRUE);
	  
	  // error handling
	  if (gCentrality > 100)
	    gCentrality = -1;
	  
	  // QA for centrality estimators (only for checked centrality)
	  fHistCentStats->Fill(0.,multSelection->GetMultiplicityPercentile("V0M", kTRUE));
	  fHistCentStats->Fill(1.,multSelection->GetMultiplicityPercentile("V0A", kTRUE));
	  fHistCentStats->Fill(2.,multSelection->GetMultiplicityPercentile("V0C", kTRUE));
	  fHistCentStats->Fill(3.,multSelection->GetMultiplicityPercentile("FMD", kTRUE));
	  fHistCentStats->Fill(4.,multSelection->GetMultiplicityPercentile("TRK", kTRUE));
	  fHistCentStats->Fill(5.,multSelection->GetMultiplicityPercentile("TKL", kTRUE));
	  fHistCentStats->Fill(6.,multSelection->GetMultiplicityPercentile("CL0", kTRUE));
	  fHistCentStats->Fill(7.,multSelection->GetMultiplicityPercentile("CL1", kTRUE));
	  fHistCentStats->Fill(8.,multSelection->GetMultiplicityPercentile("ZNA", kTRUE));
	  fHistCentStats->Fill(9.,multSelection->GetMultiplicityPercentile("ZPA", kTRUE));
	  fHistCentStats->Fill(10.,multSelection->GetMultiplicityPercentile("V0MvsFMD", kTRUE));
	  fHistCentStats->Fill(11.,multSelection->GetMultiplicityPercentile("TKLvsV0M", kTRUE));
	  fHistCentStats->Fill(12.,multSelection->GetMultiplicityPercentile("ZEMvsZDC", kTRUE));
	  
	  // Centrality estimator USED   ++++++++++++++++++++++++++++++
	  fHistCentStatsUsed->Fill(0.,gCentrality);
	  
	  if (fUseOutOfBunchPileUpCutsLHC15o) {
	    if (TMath::Abs(multSelection->GetMultiplicityPercentile("V0M") - multSelection->GetMultiplicityPercentile("CL1")) > 7.5) {
	      fHistEventStats->Fill(9, -1);
	      return -1;
	    }
	    const Int_t nTracks = event->GetNumberOfTracks();
	    Int_t multEsd = ((AliAODHeader*)event->GetHeader())->GetNumberOfESDTracks();
	    Int_t multTPC = 0;
	    for (Int_t it = 0; it < nTracks; it++) {
	      AliAODTrack* AODTrk = (AliAODTrack*)event->GetTrack(it);
	      if (!AODTrk){ delete AODTrk; continue; }
	      if (AODTrk->TestFilterBit(128)) {multTPC++;}
	    } // end of for (Int_t it = 0; it < nTracks; it++)
	    
	    if ((multEsd - 3.38*multTPC) > 15000) {
	      fHistEventStats->Fill(10, -1);
	      return -1;
	    }
	  }
	  
	  fHistCL1vsVZEROPercentile->Fill(multSelection->GetMultiplicityPercentile("V0M"),multSelection->GetMultiplicityPercentile("CL1"));
	  if(multSelection->GetEstimator("RefMult08"))
	    fHistTPCvsVZEROMultiplicity->Fill( multSelection->GetEstimator("V0M")->GetValue(),multSelection->GetEstimator("RefMult08")->GetValue());
	  else
	    gMultiplicityFromAOD = GetReferenceMultiplicityFromAOD(event);
	  
	  gMultiplicity = multSelection->GetEstimator(fCentralityEstimator)->GetValue(); 
	  fHistMultiplicity->Fill(gMultiplicity);
	  fHistMultvsPercent->Fill(gMultiplicity, gCentrality);
	}
	
	// centrality QA (V0M)
	fHistV0M->Fill(event->GetVZEROData()->GetMTotV0A(), event->GetVZEROData()->GetMTotV0C());
      } //AOD header      
    }//AOD
    
    // calculate centrality always (not only in centrality mode)
    else if(gAnalysisLevel == "AODnano" ) { //centrality via JF workaround
      
      AliAODHeader *header = (AliAODHeader*) event->GetHeader();
      if(header){
	gCentrality = (Float_t) gROOT->ProcessLine(Form("100.0 + 100.0 * ((AliNanoAODHeader*) %p)->GetCentrality(\"%s\")", header,fCentralityEstimator.Data())) / 100 - 1.0;
	
	// QA histogram
	fHistCentStatsUsed->Fill(0.,gCentrality);
	
      }//AOD header
    }//AODnano
    
    else if(gAnalysisLevel == "ESD" || gAnalysisLevel == "MCESD"){ // centrality class for ESDs or MC-ESDs //NOT SUPPORTED FOR NEW CENTRALITY FW (RUN2)
      AliCentrality *centrality = event->GetCentrality();
      gCentrality = centrality->GetCentralityPercentile(fCentralityEstimator.Data());
      
      // QA for centrality estimators
      fHistCentStats->Fill(0.,centrality->GetCentralityPercentile("V0M"));
      fHistCentStats->Fill(1.,centrality->GetCentralityPercentile("V0A"));
      fHistCentStats->Fill(2.,centrality->GetCentralityPercentile("V0C"));
      fHistCentStats->Fill(3.,centrality->GetCentralityPercentile("FMD"));
      fHistCentStats->Fill(4.,centrality->GetCentralityPercentile("TRK"));
      fHistCentStats->Fill(5.,centrality->GetCentralityPercentile("TKL"));
      fHistCentStats->Fill(6.,centrality->GetCentralityPercentile("CL0"));
      fHistCentStats->Fill(7.,centrality->GetCentralityPercentile("CL1"));
      fHistCentStats->Fill(8.,centrality->GetCentralityPercentile("ZNA"));
      fHistCentStats->Fill(9.,centrality->GetCentralityPercentile("ZPA"));
      fHistCentStats->Fill(10.,centrality->GetCentralityPercentile("V0MvsFMD"));
      fHistCentStats->Fill(11.,centrality->GetCentralityPercentile("TKLvsV0M"));
      fHistCentStats->Fill(12.,centrality->GetCentralityPercentile("ZEMvsZDC"));
      
      // Centrality estimator USED   ++++++++++++++++++++++++++++++
      fHistCentStatsUsed->Fill(0.,centrality->GetCentralityPercentile(fCentralityEstimator.Data()));
      
      // centrality QA (V0M)
      fHistV0M->Fill(event->GetVZEROData()->GetMTotV0A(), event->GetVZEROData()->GetMTotV0C());
    }//ESD
    
    else if(gAnalysisLevel == "MC"){
      Double_t gImpactParameter = 0.;
      AliMCEvent *gMCEvent = dynamic_cast<AliMCEvent*>(event);
      if(gMCEvent){
	AliCollisionGeometry* headerH = dynamic_cast<AliCollisionGeometry*>(gMCEvent->GenEventHeader());      
	if(headerH){
	  gImpactParameter = headerH->ImpactParameter();
	  gCentrality      = gImpactParameter;
	}//MC header
      }//MC event cast
    }//MC
    
    else{
      gCentrality = -1.;
    }

    if(event->GetRunNumber()<244824) { //Run1 data. Old centrality framework. 

      // calculate reference multiplicity always (not only in multiplicity mode) //only for run 1
      if(gAnalysisLevel == "ESD" || gAnalysisLevel == "MCESD"){
	AliESDEvent* gESDEvent = dynamic_cast<AliESDEvent*>(event);
	if(gESDEvent){
	  gMultiplicity = fESDtrackCuts->GetReferenceMultiplicity(gESDEvent, AliESDtrackCuts::kTrackletsITSTPC,0.5);
	  fHistMultiplicity->Fill(gMultiplicity);
	}//AliESDevent cast
      }//ESD mode
      
      else if(gAnalysisLevel == "AOD"|| gAnalysisLevel == "MCAOD" || gAnalysisLevel == "MCAODrec" ){
	AliAODHeader *header = (AliAODHeader*) event->GetHeader();
	if ((fMultiplicityEstimator == "V0M")||
	    (fMultiplicityEstimator == "V0A")||
	    (fMultiplicityEstimator == "V0C") ||
	    (fMultiplicityEstimator == "TPC")) {
	  gMultiplicity = GetReferenceMultiplicityFromAOD(event);
	  if(fDebugLevel) Printf("Reference multiplicity (calculated): %.0f",gMultiplicity);
	}
	else {
	  if(header)
	    gMultiplicity = header->GetRefMultiplicity();
	  if(fDebugLevel) Printf("Reference multiplicity (AOD header): %.0f",gMultiplicity);
	}
	
	fHistMultiplicity->Fill(gMultiplicity);
	fHistMultvsPercent->Fill(gMultiplicity, gCentrality);
      }//AOD mode
    }
    else if(gAnalysisLevel == "MC") {
      AliMCEvent* gMCEvent = dynamic_cast<AliMCEvent*>(event);
      //Calculating the multiplicity as the number of charged primaries
      //within \pm 0.8 in eta and pT > 0.1 GeV/c
      for(Int_t iParticle = 0; iParticle < gMCEvent->GetNumberOfPrimaries(); iParticle++) {
	AliMCParticle* track = dynamic_cast<AliMCParticle *>(gMCEvent->GetTrack(iParticle));
	if (!track) {
	  AliError(Form("Could not receive particle %d", iParticle));
	  continue;
	}
	
	//exclude non stable particles
	if(!(gMCEvent->IsPhysicalPrimary(iParticle))) continue;
	
	//++++++++++++++++
	if (fMultiplicityEstimator == "V0M") {
	  if((track->Eta() > 5.1 || track->Eta() < 2.8)&&(track->Eta() < -3.7 || track->Eta() > -1.7)) 
	    continue;}
	else if (fMultiplicityEstimator == "V0A") {
	  if(track->Eta() > 5.1 || track->Eta() < 2.8)  continue;}
	else if (fMultiplicityEstimator == "V0C") {
	  if(track->Eta() > -1.7 || track->Eta() < -3.7)  continue;}
	else if (fMultiplicityEstimator == "TPC") {
	  if(track->Eta() < fEtaMin || track->Eta() > fEtaMax)  continue;
	  if(track->Pt() < fPtMin || track->Pt() > fPtMax)  continue;
	}
	else{
	  if(track->Pt() < fPtMin || track->Pt() > fPtMax)  continue;
	  if(track->Eta() < fEtaMin || track->Eta() > fEtaMax)  continue;
	}
	//++++++++++++++++
	
	if(track->Charge() == 0) continue;
	
	gMultiplicity += 1;
      }//loop over primaries
      fHistMultiplicity->Fill(gMultiplicity);
    }//MC mode
    else{
      gMultiplicity = -1;
    }
  }//which centrality framework
  
  // decide what should be returned only here
  Double_t lReturnVal = -100;
  if(fEventClass=="Multiplicity"){
    lReturnVal = gCentrality;
    if (gAnalysisLevel == "MC")
      lReturnVal = gMultiplicity;
  }else if(fEventClass=="Centrality"){
    lReturnVal = gCentrality;
  }
  return lReturnVal;
}

//________________________________________________________________________
Double_t AliAnalysisTaskBFPsi::GetReferenceMultiplicityFromAOD(AliVEvent *event){
  //Function that returns the reference multiplicity from AODs (data or reco MC)
  //Different ref. mult. implemented: V0M, V0A, V0C, TPC
  Double_t gRefMultiplicity = 0., gRefMultiplicityTPC = 0.;
  Double_t gRefMultiplicityVZERO = 0., gRefMultiplicityVZEROA = 0., gRefMultiplicityVZEROC = 0.;

  AliAODHeader *header = dynamic_cast<AliAODHeader *>(event->GetHeader());
  if(!header) {
    Printf("ERROR: AOD header not available");
    return -999;
  }
  Int_t gRunNumber = header->GetRunNumber();

  // Loop over tracks in event
  for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++) {
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(event->GetTrack(iTracks));
    if (!aodTrack) {
      AliError(Form("Could not receive track %d", iTracks));
      continue;
    }
    
    // AOD track cuts
    if(!aodTrack->TestFilterBit(fnAODtrackCutBit)) continue;
            
    if(aodTrack->Charge() == 0) continue;
    // Kinematics cuts from ESD track cuts
    if( aodTrack->Pt() < fPtMin || aodTrack->Pt() > fPtMax)      continue;
    if( aodTrack->Eta() < fEtaMin || aodTrack->Eta() > fEtaMax)  continue;
    
    //=================PID (so far only for electron rejection)==========================//
    if(fElectronRejection) {
      // get the electron nsigma
      Double_t nSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kElectron));
      
      // check only for given momentum range
      if( aodTrack->Pt() > fElectronRejectionMinPt && aodTrack->Pt() < fElectronRejectionMaxPt ){
	//look only at electron nsigma
	if(!fElectronOnlyRejection) {
	  //Make the decision based on the n-sigma of electrons
	  if(nSigma < fElectronRejectionNSigma) continue;
	}
	else {
	  Double_t nSigmaPions   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kPion));
	  Double_t nSigmaKaons   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kKaon));
	  Double_t nSigmaProtons = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kProton));
	  
	  //Make the decision based on the n-sigma of electrons exclusively ( = track not in nsigma region for other species)
	  if(nSigma < fElectronRejectionNSigma
	     && nSigmaPions   > fElectronRejectionNSigma
	     && nSigmaKaons   > fElectronRejectionNSigma
	     && nSigmaProtons > fElectronRejectionNSigma ) continue;
	}
      }
    }//electron rejection
    
    gRefMultiplicityTPC += 1.0;
  }// track loop
  
  //VZERO segmentation in two detectors (0-31: VZERO-C, 32-63: VZERO-A)
  for(Int_t iChannel = 0; iChannel < 64; iChannel++) {
    fHistVZEROSignal->Fill(iChannel,event->GetVZEROEqMultiplicity(iChannel));
    
    if(iChannel < 32) 
      gRefMultiplicityVZEROC += event->GetVZEROEqMultiplicity(iChannel);
    else if(iChannel >= 32) 
      gRefMultiplicityVZEROA += event->GetVZEROEqMultiplicity(iChannel);
  }//loop over PMTs
  
  //Equalization of gain
  Double_t gFactorA = GetEqualizationFactor(gRunNumber,"A");
  if(gFactorA != 0)
    gRefMultiplicityVZEROA /= gFactorA;
  Double_t gFactorC = GetEqualizationFactor(gRunNumber,"C");
  if(gFactorC != 0)
    gRefMultiplicityVZEROC /= gFactorC;
  if((gFactorA != 0)&&(gFactorC != 0)) 
    gRefMultiplicityVZERO = (gRefMultiplicityVZEROA/gFactorA)+(gRefMultiplicityVZEROC/gFactorC);
  
  if(fDebugLevel) 
    Printf("VZERO multiplicity: %.0f - TPC multiplicity: %.0f",gRefMultiplicityVZERO,gRefMultiplicityTPC);

  fHistTPCvsVZEROMultiplicity->Fill(gRefMultiplicityVZERO,gRefMultiplicityTPC);

  if(fMultiplicityEstimator == "TPC") 
    gRefMultiplicity = gRefMultiplicityTPC;
  else if(fMultiplicityEstimator == "V0M")
    gRefMultiplicity = gRefMultiplicityVZERO;
  else if(fMultiplicityEstimator == "V0A")
    gRefMultiplicity = gRefMultiplicityVZEROA;
  else if(fMultiplicityEstimator == "V0C")
    gRefMultiplicity = gRefMultiplicityVZEROC;
  
  return gRefMultiplicity;
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
    if(!event) {
      AliError("mcEvent not available");
      return 0x0;
    }

    AliMCEvent *gMCEvent = dynamic_cast<AliMCEvent*>(event);
    if(gMCEvent){
      AliCollisionGeometry* headerH = dynamic_cast<AliCollisionGeometry*>(gMCEvent->GenEventHeader());    
      if (headerH) {
	gReactionPlane = headerH->ReactionPlaneAngle();
	//gReactionPlane *= TMath::RadToDeg();
      }//MC header
    }//MC event cast
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
Double_t AliAnalysisTaskBFPsi::GetTrackbyTrackCorrectionMatrix( Double_t vEta, 
								Double_t vPhi, 
								Double_t vPt, 
								Short_t vCharge, 
								Double_t gCentrality) {
  // -- Get efficiency correction of particle dependent on (eta, phi, pt, charge, centrality) 

  Double_t correction = 1.;
  Int_t gCentralityInt = -1;

  for (Int_t i=0; i<fCentralityArrayBinsForCorrections-1; i++){
    if((fCentralityArrayForCorrections[i] <= gCentrality)&&(gCentrality <= fCentralityArrayForCorrections[i+1])){
      gCentralityInt = i;
      break;
    }
  }  

  // centrality not in array --> no correction
  if(gCentralityInt < 0){
    correction = 1.;
  }
  else{
    
    //Printf("//=============CENTRALITY=============// %d:",gCentralityInt);
    
    if(fHistCorrectionPlus[gCentralityInt]){
      if (vCharge > 0) {
	correction = fHistCorrectionPlus[gCentralityInt]->GetBinContent(fHistCorrectionPlus[gCentralityInt]->FindBin(vEta,vPt,vPhi));
	//Printf("CORRECTIONplus: %.2f | Centrality %d",correction,gCentralityInt);  
      }
      if (vCharge < 0) {
	correction = fHistCorrectionMinus[gCentralityInt]->GetBinContent(fHistCorrectionMinus[gCentralityInt]->FindBin(vEta,vPt,vPhi));
	//Printf("CORRECTIONminus: %.2f | Centrality %d",correction,gCentralityInt); 
      }
    }
    else {
      correction = 1.;
    }
  }//centrality in array
  
  if (correction == 0.) { 
    AliError(Form("Should not happen : bin content = 0. >> eta: %.2f | phi : %.2f | pt : %.2f | cent %d",vEta, vPhi, vPt, gCentralityInt)); 
    correction = 1.; 
  } 
  
  return correction;
}

//________________________________________________________________________
TObjArray* AliAnalysisTaskBFPsi::GetAcceptedTracks(AliVEvent *event, Double_t gCentrality, Double_t gReactionPlane){
  // Returns TObjArray with tracks after all track cuts (only for AOD!)
  // Fills QA histograms

  TString gAnalysisLevel = fBalance->GetAnalysisLevel();

  //output TObjArray holding all good tracks
  TObjArray* tracksAccepted = new TObjArray;
  tracksAccepted->SetOwner(kTRUE);

  Short_t vCharge;
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
      //fHistTrackStats->Fill(aodTrack->GetFilterMap());
      for(Int_t iTrackBit = 0; iTrackBit < 16; iTrackBit++){
	fHistTrackStats->Fill(iTrackBit,aodTrack->TestFilterBit(1<<iTrackBit));
      }

      if(!aodTrack->TestFilterBit(fnAODtrackCutBit)) continue;


      // additional check on kPrimary flag
      if(fCheckPrimaryFlagAOD){
	if(aodTrack->GetType() != AliAODTrack::kPrimary)
	  continue;
      }

     
      vCharge = aodTrack->Charge();
      vEta    = aodTrack->Eta();
      vPhi    = aodTrack->Phi();// * TMath::RadToDeg();
      vPt     = aodTrack->Pt();
      vY = log( ( sqrt(fMassParticleOfInterest*fMassParticleOfInterest + vPt*vPt*cosh(vEta)*cosh(vEta)) + vPt*sinh(vEta) ) / sqrt(fMassParticleOfInterest*fMassParticleOfInterest + vPt*vPt) ); // convert eta to y; be aware that this works only for mass assumption of POI 
      
      

      //===========================PID (so far only for electron rejection)===============================//		    
      if(fElectronRejection) {

	// get the electron nsigma
	Double_t nSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kElectron));
	
	//Fill QA before the PID
	fHistdEdxVsPTPCbeforePIDelectron -> Fill(aodTrack->P()*aodTrack->Charge(),aodTrack->GetTPCsignal());
	fHistNSigmaTPCvsPtbeforePIDelectron -> Fill(aodTrack->P()*aodTrack->Charge(),nSigma); 
	//end of QA-before pid
	
	// check only for given momentum range
	if( vPt > fElectronRejectionMinPt && vPt < fElectronRejectionMaxPt ){
	  	  
	  //look only at electron nsigma
	  if(!fElectronOnlyRejection){
	    
	    //Make the decision based on the n-sigma of electrons
	    if(nSigma < fElectronRejectionNSigma) continue;
	  }
	  else{
	    
	    Double_t nSigmaPions   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kPion));
	    Double_t nSigmaKaons   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kKaon));
	    Double_t nSigmaProtons = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kProton));
	    
	    //Make the decision based on the n-sigma of electrons exclusively ( = track not in nsigma region for other species)
	    if(nSigma < fElectronRejectionNSigma
	       && nSigmaPions   > fElectronRejectionNSigma
	       && nSigmaKaons   > fElectronRejectionNSigma
	       && nSigmaProtons > fElectronRejectionNSigma ) continue;
	  }
	}
  
	//Fill QA after the PID
	fHistdEdxVsPTPCafterPIDelectron -> Fill(aodTrack->P()*aodTrack->Charge(),aodTrack->GetTPCsignal());
	fHistNSigmaTPCvsPtafterPIDelectron -> Fill(aodTrack->P()*aodTrack->Charge(),nSigma); 
	
      }
      //===========================end of PID (so far only for electron rejection)===============================//

      //+++++++++++++++++++++++++++++//
      //===========================PID===============================//		    
      if(fUsePID) {
	Double_t prob[AliPID::kSPECIES]={0.};
	Double_t probTPC[AliPID::kSPECIES]={0.};
	Double_t probTOF[AliPID::kSPECIES]={0.};
	Double_t probTPCTOF[AliPID::kSPECIES]={0.};
	
	Double_t nSigma = 0.;
	Double_t nSigmaTPC = 0.;
	Double_t nSigmaTOF = 0.; 
	Double_t nSigmaTPCTOF = 0.;
	Double_t nSigmaTPCTOFreq = 0.;
	UInt_t detUsedTPC = 0;
	UInt_t detUsedTOF = 0;
	UInt_t detUsedTPCTOF = 0;
	
	nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)fParticleOfInterest);
	nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(aodTrack,(AliPID::EParticleType)fParticleOfInterest);
	nSigmaTPCTOF = TMath::Sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF);
	nSigmaTPCTOFreq = nSigmaTPCTOF;
	if (nSigmaTOF == 999 ||  nSigmaTOF == -999){
	  nSigmaTPCTOF = nSigmaTPC;
	}


	//Decide what detector configuration we want to use
	switch(fPidDetectorConfig) {
	case kTPCpid:
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
	  nSigma = nSigmaTPC;
	  detUsedTPC = fPIDCombined->ComputeProbabilities(aodTrack, fPIDResponse, probTPC);
	  for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
	    prob[iSpecies] = probTPC[iSpecies];
	  break;
	case kTOFpid:
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
	  nSigma = nSigmaTOF;
	  detUsedTOF = fPIDCombined->ComputeProbabilities(aodTrack, fPIDResponse, probTOF);
	  for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
	    prob[iSpecies] = probTOF[iSpecies];
	  break;
	case kTPCTOF:
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
	  nSigma = nSigmaTPCTOF;
	  detUsedTPCTOF = fPIDCombined->ComputeProbabilities(aodTrack, fPIDResponse, probTPCTOF);
	  for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
	    prob[iSpecies] = probTPCTOF[iSpecies];
	  break;
	case kTPCTOFreq:
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
	  nSigma = nSigmaTPCTOFreq;
	  detUsedTPCTOF = fPIDCombined->ComputeProbabilities(aodTrack, fPIDResponse, probTPCTOF);
	  for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
	    prob[iSpecies] = probTPCTOF[iSpecies];
	  break;
	case kTPCTOFexcl:// as a starting point, more details below
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
	  nSigma = nSigmaTPC;
	  detUsedTPC = fPIDCombined->ComputeProbabilities(aodTrack, fPIDResponse, probTPC);
	  for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
	    prob[iSpecies] = probTPC[iSpecies];
	  break;
	default:
	  break;
	}//end switch: define detector mask
	
	//Filling the PID QA
	Double_t tofTime = -999., length = 999., tof = -999.;
	Double_t c = TMath::C()*1.E-9;// m/ns
	Double_t beta = -999.;


	if ((aodTrack->IsOn(AliAODTrack::kITSin)) &&  (aodTrack->IsOn(AliAODTrack::kTOFout)) ) { //check if track goes from ITS to matched TOF hit
	  
	  tofTime = aodTrack->GetTOFsignal();//in ps
	  length = aodTrack->GetIntegratedLength();
	  tof = tofTime*1E-3; // ns		      
	  if (tof <= 0) {
	    Printf("WARNING: track with negative TOF time found! Skipping this track for PID checks\n");
	    continue;
	  }
	  if (length <= 0){
	    // in old productions integrated track length is not stored in AODs -> need workaround
	    Double_t exptime[10];
	    aodTrack->GetIntegratedTimes(exptime);
	    length = exptime[0]*c*1E-3/0.01; //assume electrons are relativistic (and add all multiplication factors)
	    if (length <= 0){
	      Printf("WARNING: track with negative length found!Skipping this track for PID checks\n");
	      continue;
	    }
	  }	      
	  length = length*0.01; // in meters
	  tof = tof*c;
	  beta = length/tof;
	  
	  fHistBetavsPTOFbeforePID ->Fill(aodTrack->P()*aodTrack->Charge(),beta);
	  fHistProbTOFvsPtbeforePID ->Fill(aodTrack->Pt(),probTOF[fParticleOfInterest]);
	  fHistNSigmaTOFvsPtbeforePID ->Fill(aodTrack->Pt(),nSigmaTOF);
	}//TOF signal 
	
	fHistdEdxVsPTPCbeforePID -> Fill(aodTrack->GetTPCmomentum()*aodTrack->Charge(),aodTrack->GetTPCsignal()); //aodTrack->P()*aodTrack->Charge()
	fHistProbTPCvsPtbeforePID -> Fill(aodTrack->Pt(),probTPC[fParticleOfInterest]); 
	fHistNSigmaTPCvsPtbeforePID -> Fill(aodTrack->Pt(),nSigmaTPC);	
	fHistProbTPCTOFvsPtbeforePID -> Fill(aodTrack->Pt(),probTPCTOF[fParticleOfInterest]);
	
	//combined TPC&TOF
	fHistBetaVsdEdXbeforePID->Fill(aodTrack->GetTPCsignal(),beta); 	
	fHistNSigmaTPCTOFvsPtbeforePID -> Fill(aodTrack->Pt(),nSigmaTPCTOF);
	fHistNSigmaTPCTOFPbefPID ->Fill(nSigmaTPC,nSigmaTOF,aodTrack->P()); //++++++++++++++
	//Printf("NSIGMA: %lf - nSigmaTOF: %lf - nSigmaTPC: %lf - nSigmaTPCTOF: %lf",nSigma,nSigmaTOF,nSigmaTPC,nSigmaTPCTOF); 
	//end of QA-before pid

	Float_t probMis = fPIDResponse->GetTOFMismatchProbability(aodTrack);
	if (probMis < 0.01) { //if u want to reduce mismatch using also TPC
	  
	  if ((detUsedTPC != 0)||(detUsedTOF != 0)||(detUsedTPCTOF != 0)) {
	    
	    //Make the decision based on the n-sigma
	    if(fUsePIDnSigma) {
	      
	      // exclusive determination of PID (also taking into account other particle species, works only for pions so far)
	      if(fPidDetectorConfig == kTPCTOFexcl){
		
		Double_t nSigmaPionTPC   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kPion));
		Double_t nSigmaKaonTPC   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kKaon));
		Double_t nSigmaProtonTPC = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kProton));
		
		if( vPt > 0.2 && vPt < 0.6 ){
		  
		  if( nSigmaPionTPC > fPIDNSigma || nSigmaKaonTPC < fPIDNSigma || nSigmaProtonTPC < fPIDNSigma )
		    continue;
		}
		else if(vPt > 0.6){
		  
		  Double_t nSigmaPionTOF   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodTrack,AliPID::kPion));
		  Double_t nSigmaKaonTOF   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodTrack,AliPID::kKaon));
		  Double_t nSigmaProtonTOF = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodTrack,AliPID::kProton));
		  
		  if( nSigmaPionTPC > fPIDNSigma )
		    continue;	
		  
		  if( nSigmaPionTOF > fPIDNSigma || nSigmaKaonTOF < fPIDNSigma || nSigmaProtonTOF < fPIDNSigma )
		    continue;	
		  
		}
		else{
		  continue;
		}
	      }
	      // else just normal nSigma cut
	      else{
		if(TMath::Abs(nSigma) > fPIDNSigma) 
		  continue;
	      }	    
	      
	      fHistNSigmaTOFvsPtafterPID ->Fill(aodTrack->Pt(),nSigmaTOF);
	      fHistNSigmaTPCvsPtafterPID ->Fill(aodTrack->Pt(),nSigmaTPC);
	      fHistNSigmaTPCTOFvsPtafterPID ->Fill(aodTrack->Pt(),nSigmaTPCTOF);
	      fHistNSigmaTPCTOFPafterPID ->Fill(nSigmaTPC,nSigmaTOF,aodTrack->P());  //++++++++++++++
	    }
	    //Make the decision based on the bayesian
	    else if(fUsePIDPropabilities) {
	      if(fParticleOfInterest != TMath::LocMax(AliPID::kSPECIES,prob)) continue;
	      if (prob[fParticleOfInterest] < fMinAcceptedPIDProbability) continue;      
	      
	      fHistProbTOFvsPtafterPID ->Fill(aodTrack->Pt(),probTOF[fParticleOfInterest]);
	      fHistProbTPCvsPtafterPID ->Fill(aodTrack->Pt(),probTPC[fParticleOfInterest]); 
	      fHistProbTPCTOFvsPtafterPID ->Fill(aodTrack->Pt(),probTPCTOF[fParticleOfInterest]);
	    }	   

	    //Fill QA after the PID
	    fHistBetavsPTOFafterPID ->Fill(aodTrack->P()*aodTrack->Charge(),beta);
	    fHistdEdxVsPTPCafterPID ->Fill(aodTrack->P()*aodTrack->Charge(),aodTrack->GetTPCsignal());
	    fHistBetaVsdEdXafterPID ->Fill(aodTrack->GetTPCsignal(),beta);
	  
	  }
	  // if no detector flag remove track
	  else{
	    continue;
	  }
	}
	// if probably mismatch remove track
	else{
	  continue;
	}
      }
      //===========================PID===============================//
      //+++++++++++++++++++++++++++++//

       // Kinematics cuts from ESD track cuts
      if( vPt < fPtMin || vPt > fPtMax)      continue;
      if( vEta < fEtaMin || vEta > fEtaMax)  continue;

      // for extra DCA cuts
      Double_t pos[3];
      Double_t v[3];
      Float_t dcaXY = 0.;
      Float_t dcaZ  = 0.;
      
      // for constrained TPConly tracks
      if(fnAODtrackCutBit == 128){
	dcaXY = aodTrack->DCA();      // this is the DCA from global track (not exactly what is cut on)
	dcaZ  = aodTrack->ZAtDCA();   // this is the DCA from global track (not exactly what is cut on)
      }
      else{
	const AliVVertex *vertex = event->GetPrimaryVertex();
	vertex->GetXYZ(v);
	aodTrack->GetXYZ(pos);
	dcaXY  = TMath::Sqrt((pos[0] - v[0])*(pos[0] - v[0]) + (pos[1] - v[1])*(pos[1] - v[1]));
	dcaZ   = pos[2] - v[2];
      }

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

      // Extra cut on shared clusters
      if( fTPCsharedCut != -1 && aodTrack->GetTPCnclsS() > fTPCsharedCut){
	continue;
      }

      // fill QA histograms
      fHistClus->Fill(aodTrack->GetITSNcls(),aodTrack->GetTPCNcls());
      fHistDCA->Fill(dcaZ,dcaXY);
      fHistChi2->Fill(aodTrack->Chi2perNDF(),gCentrality);
      fHistPt->Fill(vPt,gCentrality);
      fHistEta->Fill(vEta,gCentrality);
      fHistRapidity->Fill(vY,gCentrality);
      if(vCharge > 0) fHistPhiPos->Fill(vPhi,gCentrality);
      else if(vCharge < 0) fHistPhiNeg->Fill(vPhi,gCentrality);
      fHistPhi->Fill(vPhi,gCentrality);
      if(vCharge > 0)      fHistEtaPhiPos->Fill(vEta,vPhi,gCentrality); 		 
      else if(vCharge < 0) fHistEtaPhiNeg->Fill(vEta,vPhi,gCentrality);
      
      //=======================================correction
      Double_t correction = GetTrackbyTrackCorrectionMatrix(vEta, vPhi, vPt, vCharge, gCentrality);  
      //Printf("CORRECTIONminus: %.2f | Centrality %lf",correction,gCentrality);
      
      // add the track to the TObjArray
      if(fUseRapidity){// use rapidity instead of pseudorapidity in correlation histograms
	tracksAccepted->Add(new AliBFBasicParticle(vY, vPhi, vPt, vCharge, correction)); 
      } 
      else{
	tracksAccepted->Add(new AliBFBasicParticle(vEta, vPhi, vPt, vCharge, correction)); 
      } 
    }//track loop
  }// AOD analysis


  // nano AODs
  else if(gAnalysisLevel == "AODnano") { // not fully supported yet (PID missing)
    // Loop over tracks in event
    
    for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++) {
      AliVTrack* aodTrack = dynamic_cast<AliVTrack *>(event->GetTrack(iTracks));
      if (!aodTrack) {
	AliError(Form("Could not receive track %d", iTracks));
	continue;
      }
      
      // AOD track cuts (not needed)
      //if(!aodTrack->TestFilterBit(fnAODtrackCutBit)) continue;
     
      vCharge = aodTrack->Charge();
      vEta    = aodTrack->Eta();
      vPhi    = aodTrack->Phi();// * TMath::RadToDeg();
      vPt     = aodTrack->Pt();
      vY = log( ( sqrt(fMassParticleOfInterest*fMassParticleOfInterest + vPt*vPt*cosh(vEta)*cosh(vEta)) + vPt*sinh(vEta) ) / sqrt(fMassParticleOfInterest*fMassParticleOfInterest + vPt*vPt) ); // convert eta to y; be aware that this works only for mass assumption of POI 
           
      
      // Kinematics cuts from ESD track cuts
      if( vPt < fPtMin || vPt > fPtMax)      continue;
      if( vEta < fEtaMin || vEta > fEtaMax)  continue;
      
       
      // fill QA histograms
      fHistPt->Fill(vPt,gCentrality);
      fHistEta->Fill(vEta,gCentrality);
      fHistRapidity->Fill(vY,gCentrality);
      if(vCharge > 0) fHistPhiPos->Fill(vPhi,gCentrality);
      else if(vCharge < 0) fHistPhiNeg->Fill(vPhi,gCentrality);
      fHistPhi->Fill(vPhi,gCentrality);
      if(vCharge > 0)      fHistEtaPhiPos->Fill(vEta,vPhi,gCentrality); 		 
      else if(vCharge < 0) fHistEtaPhiNeg->Fill(vEta,vPhi,gCentrality);
      
      //=======================================correction
      Double_t correction = GetTrackbyTrackCorrectionMatrix(vEta, vPhi, vPt, vCharge, gCentrality);  
      //Printf("CORRECTIONminus: %.2f | Centrality %lf",correction,gCentrality);
      
      // add the track to the TObjArray
      if(fUseRapidity){// use rapidity instead of pseudorapidity in correlation histograms
	tracksAccepted->Add(new AliBFBasicParticle(vY, vPhi, vPt, vCharge, correction)); 
      } 
      else{
	tracksAccepted->Add(new AliBFBasicParticle(vEta, vPhi, vPt, vCharge, correction)); 
      }
    }//track loop
  }// AOD nano analysis


  //==============================================================================================================
  else if(gAnalysisLevel == "MCAOD") {
    
    AliMCEvent* mcEvent = MCEvent();
    if (!mcEvent) {
      AliError("ERROR: Could not retrieve MC event");
    }
    else{
      
      for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {
	AliAODMCParticle *aodTrack = (AliAODMCParticle*) mcEvent->GetTrack(iTracks); 
	if (!aodTrack) {
	  AliError(Form("ERROR: Could not receive track %d (mc loop)", iTracks));
	  continue;
	}
	
	if(!aodTrack->IsPhysicalPrimary()) continue;
	
	vCharge = aodTrack->Charge();
	vEta    = aodTrack->Eta();
	vY      = aodTrack->Y();//true Y
	vPhi    = aodTrack->Phi();// * TMath::RadToDeg();
	vPt     = aodTrack->Pt();
	
	// Kinematics cuts from ESD track cuts
	if( vPt < fPtMin || vPt > fPtMax)      continue;
	if( vEta < fEtaMin || vEta > fEtaMax)  continue;
	
	// Remove neutral tracks
	if( vCharge == 0 ) continue;
	
	//Exclude resonances
	if(fExcludeResonancesInMC) {
	  
	  Bool_t kExcludeParticle = kFALSE;
	  Int_t gMotherIndex = aodTrack->GetMother();
	  if(gMotherIndex != -1) {
	    AliAODMCParticle* motherTrack = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(gMotherIndex));
	    if(motherTrack) {
	      Int_t pdgCodeOfMother = motherTrack->GetPdgCode();
	      //if((pdgCodeOfMother == 113)||(pdgCodeOfMother == 213)||(pdgCodeOfMother == 221)||(pdgCodeOfMother == 223)||(pdgCodeOfMother == 331)||(pdgCodeOfMother == 333)) {
	      //if(pdgCodeOfMother == 113) {
	      if(pdgCodeOfMother == 113  // rho0
		 || pdgCodeOfMother == 213 || pdgCodeOfMother == -213 // rho+
		 // || pdgCodeOfMother == 221  // eta
		 // || pdgCodeOfMother == 331  // eta'
		 // || pdgCodeOfMother == 223  // omega
		 // || pdgCodeOfMother == 333  // phi
		 || pdgCodeOfMother == 311  || pdgCodeOfMother == -311 // K0
		 // || pdgCodeOfMother == 313  || pdgCodeOfMother == -313 // K0*
		 // || pdgCodeOfMother == 323  || pdgCodeOfMother == -323 // K+*
		 || pdgCodeOfMother == 3122 || pdgCodeOfMother == -3122 // Lambda
		 || pdgCodeOfMother == 111  // pi0 Dalitz
		 || pdgCodeOfMother == 22  // photon
		 ) {
		kExcludeParticle = kTRUE;
	      }
	    }
	  }
	  
	  //Exclude from the analysis decay products of rho0, rho+, eta, eta' and phi
	  if(kExcludeParticle) continue;
	}

	//Exclude electrons with PDG
	if(fExcludeElectronsInMC) {
	  
	  if(TMath::Abs(aodTrack->GetPdgCode()) == 11) continue;
	  
	}
	
	// fill QA histograms
	fHistPt->Fill(vPt,gCentrality);
	fHistEta->Fill(vEta,gCentrality);
	fHistRapidity->Fill(vY,gCentrality);
	if(vCharge > 0) fHistPhiPos->Fill(vPhi,gCentrality);
	else if(vCharge < 0) fHistPhiNeg->Fill(vPhi,gCentrality);
	fHistPhi->Fill(vPhi,gCentrality);
	if(vCharge > 0)      fHistEtaPhiPos->Fill(vEta,vPhi,gCentrality); 		 
	else if(vCharge < 0) fHistEtaPhiNeg->Fill(vEta,vPhi,gCentrality);
	
	//=======================================correction
	Double_t correction = GetTrackbyTrackCorrectionMatrix(vEta, vPhi, vPt, vCharge, gCentrality);  
	//Printf("CORRECTIONminus: %.2f | Centrality %lf",correction,gCentrality);   
	
	// add the track to the TObjArray
	if(fUseRapidity){// use rapidity instead of pseudorapidity in correlation histograms
	  tracksAccepted->Add(new AliBFBasicParticle(vY, vPhi, vPt, vCharge, correction)); 
	} 
	else{
	  tracksAccepted->Add(new AliBFBasicParticle(vEta, vPhi, vPt, vCharge, correction)); 
	}
      }//aodTracks
    }//MC event
  }//MCAOD
  //==============================================================================================================

  //==============================================================================================================
  else if(gAnalysisLevel == "MCAODrec") {
    
    /* fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAOD) {
      printf("ERROR: fAOD not available\n");
      return;
      }*/
    
    fArrayMC = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName())); 
    if (!fArrayMC) { 
       AliError("No array of MC particles found !!!");
    }

    AliMCEvent* mcEvent = MCEvent();
    if (!mcEvent) {
       AliError("ERROR: Could not retrieve MC event");
       return tracksAccepted;
    }
     
    for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++) {
      AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(event->GetTrack(iTracks));
      if (!aodTrack) {
	AliError(Form("Could not receive track %d", iTracks));
	continue;
      }
      
      for(Int_t iTrackBit = 0; iTrackBit < 16; iTrackBit++){
	fHistTrackStats->Fill(iTrackBit,aodTrack->TestFilterBit(1<<iTrackBit));
      }
      if(!aodTrack->TestFilterBit(fnAODtrackCutBit)) continue;
      
      vCharge = aodTrack->Charge();
      vEta    = aodTrack->Eta();
      vPhi    = aodTrack->Phi();// * TMath::RadToDeg();
      vPt     = aodTrack->Pt();
      vY = log( ( sqrt(fMassParticleOfInterest*fMassParticleOfInterest + vPt*vPt*cosh(vEta)*cosh(vEta)) + vPt*sinh(vEta) ) / sqrt(fMassParticleOfInterest*fMassParticleOfInterest + vPt*vPt) ); // convert eta to y; be aware that this works only for mass assumption of POI 

      //analyze one set of particles
      if(fUseMCPdgCode) {

	Int_t label = TMath::Abs(aodTrack->GetLabel());
	AliAODMCParticle *AODmcTrackForPID = (AliAODMCParticle*) fArrayMC->At(label);

	if(!AODmcTrackForPID){
	  AliError(Form("No AliAODMCParticle for aodTrack with label %d ... skip",label));
	  continue;
	}

	Int_t gPdgCode = AODmcTrackForPID->PdgCode();
	if(TMath::Abs(fPDGCodeToBeAnalyzed) != TMath::Abs(gPdgCode)) 
	  continue;
      }
      
      //===========================use MC information for Kinematics===============================//		    
      if(fUseMCforKinematics){

	Int_t label = TMath::Abs(aodTrack->GetLabel());
	AliAODMCParticle *AODmcTrack = (AliAODMCParticle*) fArrayMC->At(label);

	if(AODmcTrack){
	  vCharge = AODmcTrack->Charge();
	  vEta    = AODmcTrack->Eta();
	  vY      = AODmcTrack->Y();//true Y
	  vPhi    = AODmcTrack->Phi();// * TMath::RadToDeg();
	  vPt     = AODmcTrack->Pt();
	}
	else{
	  AliDebug(1, "no MC particle for this track"); 
	  continue;
	}
      }
      //===========================end of use MC information for Kinematics========================//		    


      //===========================PID (so far only for electron rejection)===============================//		    
      if(fElectronRejection) {

	// get the electron nsigma
	Double_t nSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kElectron));

	//Fill QA before the PID
	fHistdEdxVsPTPCbeforePIDelectron -> Fill(aodTrack->P()*aodTrack->Charge(),aodTrack->GetTPCsignal());
	fHistNSigmaTPCvsPtbeforePIDelectron -> Fill(aodTrack->P()*aodTrack->Charge(),nSigma); 
	//end of QA-before pid
	
	// check only for given momentum range
	if( vPt > fElectronRejectionMinPt && vPt < fElectronRejectionMaxPt ){
	  	  
	  //look only at electron nsigma
	  if(!fElectronOnlyRejection){
	    
	    //Make the decision based on the n-sigma of electrons
	    if(nSigma < fElectronRejectionNSigma) continue;
	  }
	  else{
	    
	    Double_t nSigmaPions   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kPion));
	    Double_t nSigmaKaons   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kKaon));
	    Double_t nSigmaProtons = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kProton));
	    
	    //Make the decision based on the n-sigma of electrons exclusively ( = track not in nsigma region for other species)
	    if(nSigma < fElectronRejectionNSigma
	       && nSigmaPions   > fElectronRejectionNSigma
	       && nSigmaKaons   > fElectronRejectionNSigma
	       && nSigmaProtons > fElectronRejectionNSigma ) continue;
	  }
	}
  
	//Fill QA after the PID
	fHistdEdxVsPTPCafterPIDelectron -> Fill(aodTrack->P()*aodTrack->Charge(),aodTrack->GetTPCsignal());
	fHistNSigmaTPCvsPtafterPIDelectron -> Fill(aodTrack->P()*aodTrack->Charge(),nSigma); 
	
      }
      //===========================end of PID (so far only for electron rejection)===============================//
      
      // Kinematics cuts from ESD track cuts
      if( vPt < fPtMin || vPt > fPtMax)      continue;
      if( vEta < fEtaMin || vEta > fEtaMax)  continue;

      // for extra DCA cuts
      Double_t pos[3];
      Double_t v[3];
      Float_t dcaXY = 0.;
      Float_t dcaZ  = 0.;
      
      // for constrained TPConly tracks
      if(fnAODtrackCutBit == 128){
	dcaXY = aodTrack->DCA();      // this is the DCA from global track (not exactly what is cut on)
	dcaZ  = aodTrack->ZAtDCA();   // this is the DCA from global track (not exactly what is cut on)
      }
      else{
	const AliVVertex *vertex = event->GetPrimaryVertex();
	vertex->GetXYZ(v);
	aodTrack->GetXYZ(pos);
	dcaXY  = TMath::Sqrt((pos[0] - v[0])*(pos[0] - v[0]) + (pos[1] - v[1])*(pos[1] - v[1]));
	dcaZ   = pos[2] - v[2];
      }
      
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

     // Extra cut on shared clusters
      if( fTPCsharedCut != -1 && aodTrack->GetTPCnclsS() > fTPCsharedCut){
	continue;
      }

      //Exclude secondaries from material and weak decays
      if(fExcludeSecondariesInMC){
	
	Int_t label = TMath::Abs(aodTrack->GetLabel());
	AliAODMCParticle *AODmcTrack = (AliAODMCParticle*) fArrayMC->At(label);

	if (!AODmcTrack->IsPhysicalPrimary())
	  continue;   
      }

      //Exclude resonances
      if(fExcludeResonancesInMC) {
	
	Bool_t kExcludeParticle = kFALSE;

	Int_t label = TMath::Abs(aodTrack->GetLabel());
	AliAODMCParticle *AODmcTrack = (AliAODMCParticle*) fArrayMC->At(label);
      
        if (AODmcTrack){ 
	  
	  Int_t gMotherIndex = AODmcTrack->GetMother();
	  if(gMotherIndex != -1) {
	    AliAODMCParticle* motherTrack = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(gMotherIndex));
	    if(motherTrack) {
	      Int_t pdgCodeOfMother = motherTrack->GetPdgCode();
	      if(pdgCodeOfMother == 113  // rho0
		 // || pdgCodeOfMother == 213 || pdgCodeOfMother == -213 // rho+
		 // || pdgCodeOfMother == 221  // eta
		 // || pdgCodeOfMother == 331  // eta'
		 // || pdgCodeOfMother == 223  // omega
		 // || pdgCodeOfMother == 333  // phi
		 || pdgCodeOfMother == 310  || pdgCodeOfMother == -310 // K0
		 // || pdgCodeOfMother == 313  || pdgCodeOfMother == -313 // K0*
		 // || pdgCodeOfMother == 323  || pdgCodeOfMother == -323 // K+*
		 || pdgCodeOfMother == 3122 || pdgCodeOfMother == -3122 // Lambda
		 //|| pdgCodeOfMother == 111  // pi0 Dalitz
		 //|| pdgCodeOfMother == 22   // photon
		 ) {
		kExcludeParticle = kTRUE;
	      }
	    }
	  }
	}	
	//Exclude from the analysis decay products of rho0, rho+, eta, eta' and phi
	if(kExcludeParticle) continue;
      }

      //Include exclusively resonances with a specific PDG value
      if(fIncludeResonancePDGInMC > -1) {
	
	Bool_t kIncludeParticle = kFALSE;
	
	Int_t label = TMath::Abs(aodTrack->GetLabel());
	AliAODMCParticle *AODmcTrack = (AliAODMCParticle*) fArrayMC->At(label);
	
        if (AODmcTrack){ 
	  
	  Int_t gMotherIndex = AODmcTrack->GetMother();
	  if(gMotherIndex != -1) {
	    AliAODMCParticle* motherTrack = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(gMotherIndex));
	    if(motherTrack) {
	      Int_t pdgCodeOfMother = motherTrack->GetPdgCode();
	      if(TMath::Abs(pdgCodeOfMother) == fIncludeResonancePDGInMC) {
		kIncludeParticle = kTRUE;
	      }
	    }
	  }
	}	
	
	//Exclude from the analysis particle that are not decay products from this resonance
	if(!kIncludeParticle) continue;
      }
      
      //Exclude electrons with PDG
      if(fExcludeElectronsInMC) {
	
	Int_t label = TMath::Abs(aodTrack->GetLabel());
	AliAODMCParticle *AODmcTrack = (AliAODMCParticle*) fArrayMC->At(label);
	
        if (AODmcTrack){ 
	  if(TMath::Abs(AODmcTrack->GetPdgCode()) == 11) continue;
	}
      }
      
      // fill QA histograms
      fHistClus->Fill(aodTrack->GetITSNcls(),aodTrack->GetTPCNcls());
      fHistDCA->Fill(dcaZ,dcaXY);
      fHistChi2->Fill(aodTrack->Chi2perNDF(),gCentrality);
      fHistPt->Fill(vPt,gCentrality);
      fHistEta->Fill(vEta,gCentrality);
      fHistRapidity->Fill(vY,gCentrality);
      if(vCharge > 0) fHistPhiPos->Fill(vPhi,gCentrality);
      else if(vCharge < 0) fHistPhiNeg->Fill(vPhi,gCentrality);
      fHistPhi->Fill(vPhi,gCentrality);
      if(vCharge > 0)      fHistEtaPhiPos->Fill(vEta,vPhi,gCentrality); 		 
      else if(vCharge < 0) fHistEtaPhiNeg->Fill(vEta,vPhi,gCentrality);
      
      //=======================================correction
      Double_t correction = GetTrackbyTrackCorrectionMatrix(vEta, vPhi, vPt, vCharge, gCentrality);  
      //Printf("CORRECTIONminus: %.2f | Centrality %lf",correction,gCentrality);
      
      // add the track to the TObjArray
      if(fUseRapidity){// use rapidity instead of pseudorapidity in correlation histograms
	tracksAccepted->Add(new AliBFBasicParticle(vY, vPhi, vPt, vCharge, correction)); 
      } 
      else{
	tracksAccepted->Add(new AliBFBasicParticle(vEta, vPhi, vPt, vCharge, correction)); 
      } 
    }//track loop
  }//MCAODrec
  //==============================================================================================================

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
	Double_t nSigmaTPC = 0.;
	Double_t nSigmaTOF = 0.; 
	Double_t nSigmaTPCTOF = 0.;
	Double_t nSigmaTPCTOFreq = 0.;
	UInt_t detUsedTPC = 0;
	UInt_t detUsedTOF = 0;
	UInt_t detUsedTPCTOF = 0;
	
	nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)fParticleOfInterest);
	nSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)fParticleOfInterest);
	nSigmaTPCTOF = TMath::Sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF);
	nSigmaTPCTOFreq = nSigmaTPCTOF;
	if (nSigmaTOF == 999 ||  nSigmaTOF == -999){
	  nSigmaTPCTOF = nSigmaTPC;
	}
	
	//Decide what detector configuration we want to use
	switch(fPidDetectorConfig) {
	case kTPCpid:
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
	  nSigma = nSigmaTPC;
	  detUsedTPC = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPC);
	  for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
	    prob[iSpecies] = probTPC[iSpecies];
	  break;
	case kTOFpid:
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
	  nSigma = nSigmaTOF;
	  detUsedTOF = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTOF);
	  for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
	    prob[iSpecies] = probTOF[iSpecies];
	  break;
	case kTPCTOF:
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
	  nSigma = nSigmaTPCTOF;
	  detUsedTPCTOF = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPCTOF);
	  for(Int_t iSpecies = 0; iSpecies < AliPID::kSPECIES; iSpecies++)
	    prob[iSpecies] = probTPCTOF[iSpecies];
	  break;
	case kTPCTOFreq:
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
	  nSigma = nSigmaTPCTOFreq;
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
	    Printf("WARNING: track with negative TOF time found! Skipping this track for PID checks\n");
	    continue;
	  }
	  if (length <= 0){
	    Printf("WARNING: track with negative length found!Skipping this track for PID checks\n");
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
      vEta    = trackTPC->Eta();
      vPhi    = trackTPC->Phi();// * TMath::RadToDeg();
      vPt     = trackTPC->Pt();
      vY = log( ( sqrt(fMassParticleOfInterest*fMassParticleOfInterest + vPt*vPt*cosh(vEta)*cosh(vEta)) + vPt*sinh(vEta) ) / sqrt(fMassParticleOfInterest*fMassParticleOfInterest + vPt*vPt) ); // convert eta to y; be aware that this works only for mass assumption of POI 

      fHistClus->Fill(trackTPC->GetITSclusters(0),nClustersTPC);
      fHistDCA->Fill(b[1],b[0]);
      fHistChi2->Fill(chi2PerClusterTPC,gCentrality);
      fHistPt->Fill(vPt,gCentrality);
      fHistEta->Fill(vEta,gCentrality);
      fHistPhi->Fill(vPhi,gCentrality);
      if(vCharge > 0)      fHistEtaPhiPos->Fill(vEta,vPhi,gCentrality);
      else if(vCharge < 0) fHistEtaPhiNeg->Fill(vEta,vPhi,gCentrality);
      fHistRapidity->Fill(vY,gCentrality);
      if(vCharge > 0) fHistPhiPos->Fill(vPhi,gCentrality);
      else if(vCharge < 0) fHistPhiNeg->Fill(vPhi,gCentrality);
      
      //=======================================correction
      Double_t correction = GetTrackbyTrackCorrectionMatrix(vEta, vPhi, vPt, vCharge, gCentrality);  
      //Printf("CORRECTIONminus: %.2f | Centrality %lf",correction,gCentrality);

      // add the track to the TObjArray
      if(fUseRapidity){// use rapidity instead of pseudorapidity in correlation histograms
	tracksAccepted->Add(new AliBFBasicParticle(vY, vPhi, vPt, vCharge, correction)); 
      } 
      else{
	tracksAccepted->Add(new AliBFBasicParticle(vEta, vPhi, vPt, vCharge, correction)); 
      }
      
      delete trackTPC;
    }//track loop
  }// ESD analysis

  else if(gAnalysisLevel == "MC"){
    if(!event) {
      AliError("mcEvent not available");
      return 0x0;
    }

    AliMCEvent *gMCEvent = dynamic_cast<AliMCEvent*>(event);
    if(gMCEvent) {


      Int_t nTracks = gMCEvent->GetNumberOfPrimaries(); // default mode: running over primary particles only 
      if(fIncludeSecondariesInMCgen){
	nTracks = gMCEvent->GetNumberOfTracks(); // for weak decay studies, include also secondaries
      }
      
      // Loop over tracks in event
      for (Int_t iTracks = 0; iTracks < nTracks; iTracks++) {
	AliMCParticle* track = dynamic_cast<AliMCParticle *>(gMCEvent->GetTrack(iTracks));
	if (!track) {
	  AliError(Form("Could not receive particle %d", iTracks));
	  continue;
	}
	
	//exclude non stable particles (only if not secondaries included explicitly)
	if(!fIncludeSecondariesInMCgen && !(gMCEvent->IsPhysicalPrimary(iTracks))) continue;

	// exclude particles with strange behaviour in AMPT
	// - mothers that have physical primary daughters
	if(fExcludeParticlesExtra){

	  //exclude particles that are primary and have primary daughters
	  if(track->GetFirstDaughter()!=-1){
	    if(gMCEvent->IsPhysicalPrimary(track->GetFirstDaughter()))
	      continue;
	  }
	  if(track->GetLastDaughter()!=-1){
	    if(gMCEvent->IsPhysicalPrimary(track->GetLastDaughter()))
	      continue;
	  }
	}
	
	vCharge = track->Charge();
	vEta    = track->Eta();
	vPt     = track->Pt();
	vY      = track->Y();//true Y
	
	if( vPt < fPtMin || vPt > fPtMax)      
	  continue;
	if (!fUsePID) {
	  if( vEta < fEtaMin || vEta > fEtaMax)  continue;
	}
	else if (fUsePID){
	  if( vY < fEtaMin || vY > fEtaMax)  continue;
	}

	// Remove neutral tracks
	if( vCharge == 0 ) continue;
	
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

	//Exclude weak decay products (if not done by IsPhysicalPrimary)
	if(fExcludeWeakDecaysInMC) {
	  TParticle *particle = track->Particle();
	  if(!particle) continue;
	  
	  Bool_t kExcludeParticle = kFALSE;
	  Int_t gMotherIndex = particle->GetFirstMother();
	  if(gMotherIndex != -1) {
	    AliMCParticle* motherTrack = dynamic_cast<AliMCParticle *>(event->GetTrack(gMotherIndex));
	    if(motherTrack) {
	      TParticle *motherParticle = motherTrack->Particle();
	      if(motherParticle) {
		if(IsThisAWeakDecayingParticle(motherParticle)){
		  kExcludeParticle = kTRUE;
		}
	      }
	    }
	  }
	  
	  //Exclude from the analysis decay products of weakly decaying particles
	  if(kExcludeParticle) continue;
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
		if(pdgCodeOfMother == 113  // rho0
		   || pdgCodeOfMother == 213 || pdgCodeOfMother == -213 // rho+
		   // || pdgCodeOfMother == 221  // eta
		   // || pdgCodeOfMother == 331  // eta'
		   // || pdgCodeOfMother == 223  // omega
		   // || pdgCodeOfMother == 333  // phi
		   || pdgCodeOfMother == 311  || pdgCodeOfMother == -311 // K0
		   // || pdgCodeOfMother == 313  || pdgCodeOfMother == -313 // K0*
		   // || pdgCodeOfMother == 323  || pdgCodeOfMother == -323 // K+*
		   || pdgCodeOfMother == 3122 || pdgCodeOfMother == -3122 // Lambda
		   || pdgCodeOfMother == 111  // pi0 Dalitz
		   ) {
		  kExcludeParticle = kTRUE;
		}
	      }
	    }
	  }
	  
	  //Exclude from the analysis decay products of rho0, rho+, eta, eta' and phi
	  if(kExcludeParticle) continue;
	}



	//Exclude resonances with a specific PDG value
	if(fExcludeResonancePDGInMC > -1) {

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
		if( TMath::Abs(pdgCodeOfMother) == fExcludeResonancePDGInMC ){
		  kExcludeParticle = kTRUE;
		}
	      }
	    }
	  }
	  
	  //Exclude from the analysis decay products of rho0, rho+, eta, eta' and phi
	  if(kExcludeParticle) continue;
	}

	//Include exclusively resonances with a specific PDG value
	if(fIncludeResonancePDGInMC > -1) {

	  TParticle *particle = track->Particle();
	  if(!particle) continue;
	  
	  Bool_t kIncludeParticle = kFALSE;
	  Int_t gMotherIndex = particle->GetFirstMother();
	  if(gMotherIndex != -1) {
	    AliMCParticle* motherTrack = dynamic_cast<AliMCParticle *>(event->GetTrack(gMotherIndex));
	    if(motherTrack) {
	      TParticle *motherParticle = motherTrack->Particle();
	      if(motherParticle) {

		Int_t pdgCodeOfMother = motherParticle->GetPdgCode();
		if( TMath::Abs(pdgCodeOfMother) == fIncludeResonancePDGInMC ){
		  kIncludeParticle = kTRUE;
		}
	      }
	    }
	  }
	  
	  //Exclude from the analysis particle that are not decay products from this resonance
	  if(!kIncludeParticle) continue;
	}


	//Exclude electrons with PDG
	if(fExcludeElectronsInMC) {
	  
	  TParticle *particle = track->Particle();
	  
	  if (particle){ 
	    if(TMath::Abs(particle->GetPdgCode()) == 11) continue;
	  }
	}
      
	vPhi    = track->Phi();
	//Printf("phi (before): %lf",vPhi);
	
	fHistPt->Fill(vPt,gCentrality);
	fHistEta->Fill(vEta,gCentrality);
	fHistPhi->Fill(vPhi,gCentrality);
	if(vCharge > 0)      fHistEtaPhiPos->Fill(vEta,vPhi,gCentrality);
	else if(vCharge < 0) fHistEtaPhiNeg->Fill(vEta,vPhi,gCentrality);
	//fHistPhi->Fill(vPhi*TMath::RadToDeg(),gCentrality);
	fHistRapidity->Fill(vY,gCentrality);
	//if(vCharge > 0) fHistPhiPos->Fill(vPhi*TMath::RadToDeg(),gCentrality);
	//else if(vCharge < 0) fHistPhiNeg->Fill(vPhi*TMath::RadToDeg(),gCentrality);
	if(vCharge > 0) fHistPhiPos->Fill(vPhi,gCentrality);
	else if(vCharge < 0) fHistPhiNeg->Fill(vPhi,gCentrality);
	
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
	  fHistPhiBefore->Fill(vDeltaphiBefore,gCentrality);
	  
	  Double_t vDeltaphiAfter = vPhi - gReactionPlane*TMath::DegToRad();
	  if(vDeltaphiAfter < 0) vDeltaphiAfter += 2*TMath::Pi();
	  fHistPhiAfter->Fill(vDeltaphiAfter,gCentrality);
	  
	}
	
	//vPhi *= TMath::RadToDeg();
	
      //=======================================correction
      Double_t correction = GetTrackbyTrackCorrectionMatrix(vEta, vPhi, vPt, vCharge, gCentrality);  
      //Printf("CORRECTIONminus: %.2f | Centrality %lf",correction,gCentrality);

      if(fUseRapidity){// use rapidity instead of pseudorapidity in correlation histograms
	tracksAccepted->Add(new AliBFBasicParticle(vY, vPhi, vPt, vCharge, correction)); 
      } 
      else{
	tracksAccepted->Add(new AliBFBasicParticle(vEta, vPhi, vPt, vCharge, correction)); 
      } 
      } //track loop
    }//MC event object
  }//MC
  
  return tracksAccepted;  
}

//________________________________________________________________________
TObjArray* AliAnalysisTaskBFPsi::GetShuffledTracks(TObjArray *tracks, Double_t gCentrality){
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
    //==============================correction
    Double_t correction = GetTrackbyTrackCorrectionMatrix(track->Eta(), track->Phi(),track->Pt(), chargeVector->at(i), gCentrality);
    //Printf("CORRECTIONminus: %.2f | Centrality %lf",correction,gCentrality);
    tracksShuffled->Add(new AliBFBasicParticle(track->Eta(), track->Phi(), track->Pt(),chargeVector->at(i), correction));
  }

  delete chargeVector;
   
  return tracksShuffled;
}

//________________________________________________________________________
void  AliAnalysisTaskBFPsi::SetVZEROCalibrationFile(const char* filename,
						    const char* lhcPeriod) {
  //Function to setup the VZERO gain equalization
    //============Get the equilization map============//
  TFile *calibrationFile = TFile::Open(filename);
  if((!calibrationFile)||(!calibrationFile->IsOpen())) {
    Printf("No calibration file found!!!");
    return;
  }

  TList *list = dynamic_cast<TList *>(calibrationFile->Get(lhcPeriod));
  if(!list) {
    Printf("Calibration TList not found!!!");
    return;
  }

  fHistVZEROAGainEqualizationMap = dynamic_cast<TH1F *>(list->FindObject("gHistVZEROAGainEqualizationMap"));
  if(!fHistVZEROAGainEqualizationMap) {
    Printf("VZERO-A calibration object not found!!!");
    return;
  }
  fHistVZEROCGainEqualizationMap = dynamic_cast<TH1F *>(list->FindObject("gHistVZEROCGainEqualizationMap"));
  if(!fHistVZEROCGainEqualizationMap) {
    Printf("VZERO-C calibration object not found!!!");
    return;
  }

  fHistVZEROChannelGainEqualizationMap = dynamic_cast<TH2F *>(list->FindObject("gHistVZEROChannelGainEqualizationMap"));
  if(!fHistVZEROChannelGainEqualizationMap) {
    Printf("VZERO channel calibration object not found!!!");
    return;
  }
}

//________________________________________________________________________
void AliAnalysisTaskBFPsi::SetParticleOfInterest(kParticleOfInterest poi) {

  // Function to set the particle of interest (for PID analysis)
  // and the corresponding mass
  
  fParticleOfInterest = poi;

  if(fParticleOfInterest == kParticleOfInterest::kElectron){
    fMassParticleOfInterest = TDatabasePDG::Instance()->GetParticle(11)->Mass();
  }  
  else if(fParticleOfInterest == kParticleOfInterest::kMuon){
    fMassParticleOfInterest = TDatabasePDG::Instance()->GetParticle(13)->Mass();
  }
  else if(fParticleOfInterest == kParticleOfInterest::kPion){
    fMassParticleOfInterest = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  }
  else if(fParticleOfInterest == kParticleOfInterest::kKaon){
    fMassParticleOfInterest = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  }
  else if(fParticleOfInterest == kParticleOfInterest::kProton){
    fMassParticleOfInterest = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  }
  else{
    AliWarning("Particle type not known, set fMassParticleOfInterest to pion mass.");
    fMassParticleOfInterest = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  }
}


//________________________________________________________________________
Double_t AliAnalysisTaskBFPsi::GetChannelEqualizationFactor(Int_t run, 
							    Int_t channel) {
  //
  if(!fHistVZEROAGainEqualizationMap) return 1.0;

  for(Int_t iBinX = 1; iBinX <= fHistVZEROChannelGainEqualizationMap->GetNbinsX(); iBinX++) {
    Int_t gRunNumber = atoi(fHistVZEROChannelGainEqualizationMap->GetXaxis()->GetBinLabel(iBinX));
    if(gRunNumber == run)
      return fHistVZEROChannelGainEqualizationMap->GetBinContent(iBinX,channel+1);
  }

  return 1.0;
}

//________________________________________________________________________
Double_t AliAnalysisTaskBFPsi::GetEqualizationFactor(Int_t run, 
						     const char* side) {
  //
  if(!fHistVZEROAGainEqualizationMap) return 1.0;

  TString gVZEROSide = side;
  for(Int_t iBinX = 1; iBinX < fHistVZEROAGainEqualizationMap->GetNbinsX(); iBinX++) {
    Int_t gRunNumber = atoi(fHistVZEROAGainEqualizationMap->GetXaxis()->GetBinLabel(iBinX));
    //cout<<"Looking for run "<<run<<" - current run: "<<gRunNumber<<endl;
    if(gRunNumber == run) {
      if(gVZEROSide == "A") 
	return fHistVZEROAGainEqualizationMap->GetBinContent(iBinX);
      else if(gVZEROSide == "C") 
	return fHistVZEROCGainEqualizationMap->GetBinContent(iBinX);
    }
  }

  return 1.0;
}

//____________________________________________________________________
Bool_t AliAnalysisTaskBFPsi::AcceptEventCentralityWeight(Double_t centrality)
{
  // copied from AliAnalysisTaskPhiCorrelations
  //
  // rejects "randomly" events such that the centrality gets flat
  // uses fCentralityWeights histogram

  // TODO code taken and adapted from AliRDHFCuts; waiting for general class AliCentralityFlattening
  
  Double_t weight = fCentralityWeights->GetBinContent(fCentralityWeights->FindBin(centrality));
  Double_t centralityDigits = centrality*100. - (Int_t)(centrality*100.);
  
  Bool_t result = kFALSE;
  if (centralityDigits < weight) 
    result = kTRUE;
  
  AliInfo(Form("Centrality: %f; Digits: %f; Weight: %f; Result: %d", centrality, centralityDigits, weight, result));
  
  return result;
}

//____________________________________________________________________
Bool_t AliAnalysisTaskBFPsi::IsThisAWeakDecayingParticle(TParticle *thisGuy)
{
  // In order to prevent analyzing daughters from weak decays 
  // - AMPT does not only strong decays, so IsPhysicalPrimary does not catch it

 Int_t pdgcode = TMath::Abs( thisGuy->GetPdgCode() );

 Int_t myWeakParticles[7] = { 3322, 3312, 3222, // Xi0 Xi+- Sigma-+
			       3122, 3112, // Lambda0 Sigma+-
			       130, 310 // K_L0 K_S0
 };

 Bool_t found = kFALSE;
 for(Int_t i=0; i!=7; ++i)
   if( myWeakParticles[i] == pdgcode ) {
     found = kTRUE;
     break;
   }
 
 return found;
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

