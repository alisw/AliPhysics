#ifndef ALIANALYSISTASKSIGNEDBF_H
#define ALIANALYSISTASKSIGNEDBF_H

// Analysis task for the BF vs Psi code
// Authors: Panos Cristakoglou@cern.ch

class TList;
class TH1F;
class TH2F;
class TH3F; 
class TF1;
class TH3D;
class TParticle;

class AliBalanceEbyE;
class AliBalancePsi;
class AliESDtrackCuts;
class AliEventPoolManager;
class AliAnalysisUtils;
class AliPID;

#include "AliAnalysisTaskSE.h"
#include "AliBalancePsi.h"
#include "AliAODTrack.h"

#include "AliPID.h"  
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliTimeRangeCut.h"

 
//================================correction
#define kCENTRALITY 101
#define kNBRUN 200
//const Double_t centralityArrayForPbPb[kCENTRALITY+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};
//const TString centralityArrayForPbPb_string[kCENTRALITY] = {"0-5","5-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80"};
//================================correction

class AliAnalysisTaskSignedBF : public AliAnalysisTaskSE {
 public:
  enum etriggerSel{kMB, kCentral15, kCentral18, kINT7, kppHighMult};
  enum eCorrProcedure{kNoCorr, kDataDrivCorr, kMCCorr, kMC1DCorr};
  enum ePIDSel{kTrig, kAssoc, kBoth};
  
  AliAnalysisTaskSignedBF(const char *name = "AliAnalysisTaskSignedBF");
  virtual ~AliAnalysisTaskSignedBF();
   
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   FinishTaskOutput();
  virtual void   Terminate(Option_t *);

  //========================correction

  void SetCorrectionProcedure(AliAnalysisTaskSignedBF::eCorrProcedure corrProc) {fCorrProcedure = corrProc;}
  
  virtual void SetInputCorrection(TString filename, 
				    Int_t nCentralityBins, 
				    Double_t *centralityArrayForCorrections);
  //========================correction
  // void SetDebugLevel() {fDebugLevel = kTRUE;} //hides overloaded virtual function

  //======== New methods for data driven NUA(eta, phi, vz) run-by-run, NUE (pT, cont) from MC per centrality bins 

  void SetInputListForNUACorr(TString fileNUA);
  void SetInputListForNUECorr(TString fileNUE);

  void SetInputListForNUECorr3D(TString fileNUE);
 
  Double_t GetNUACorrection(Int_t gRun, Short_t vCharge, Double_t vVz, Float_t vEta, Float_t vPhi );
  Double_t GetNUECorrection(Int_t gCentrality, Short_t vCharge, Double_t vPt, Int_t poi);

  Bool_t SetSelectPID(AliAODTrack* track, Int_t poi);  
  Bool_t SetSelectPIDKaons(AliAODTrack* track, Int_t poi);
  Int_t GetIndexRun(Int_t runNb);
  Int_t GetIndexCentrality(Double_t gCentrality);
  
  void SetCentralityArrayBins(Int_t nCentralityBins, Double_t *centralityArrayForCorrections){
    fCentralityArrayBinsForCorrections = nCentralityBins;
    for (Int_t i=0; i<=nCentralityBins-1; i++)
      fCentralityArrayForCorrections[i] = centralityArrayForCorrections[i];
  }

  void SetArrayRuns(Int_t nRuns, Int_t *runsArrayForCorrections){
    fTotalNbRun = nRuns;
    for (Int_t i=0; i<=nRuns-1; i++)
      fRunNb[i] = runsArrayForCorrections[i];
  }
 
  void SetAnalysisObject(AliBalancePsi *const analysis) {
    fBalance         = analysis;
    }
  void SetShufflingObject(AliBalancePsi *const analysisShuffled) {
    fRunShuffling = kTRUE;
    fShuffledBalance = analysisShuffled;
  }
  void SetMixingObject(AliBalancePsi *const analysisMixed) {
    fRunMixing = kTRUE;
    fMixedBalance = analysisMixed;
  }
  void SetMixingWithEventPlane(Bool_t bMixingWithEventPlane = kTRUE) { fRunMixingEventPlane = bMixingWithEventPlane; }
  void SetMixingTracks(Int_t tracks) { fMixingTracks = tracks; }
  void SetMaxNbMixedEvents(Int_t ev) { fMaxNbMixedEvents = ev; }
  void SetEbyEObject(AliBalanceEbyE *const analysisEbyE){
    fRunEbyE = kTRUE;
    fBalanceEbyE = analysisEbyE;
  }
  void SetAnalysisCutObject(AliESDtrackCuts *const trackCuts) {
    fESDtrackCuts = trackCuts;}
  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {
    fVxMax = vx;
    fVyMax = vy;
    fVzMax = vz;
  }


  void SetRequireHighPtTrigger(Double_t gpTtrigger) {
    fRequireHighPtTrigger = kTRUE;
    fPtTriggerMin = gpTtrigger;
  }
  
  //==============AOD analysis==============//
  void SetAODtrackCutBit(Int_t bit){
    fnAODtrackCutBit = bit;
  }

  void SetKinematicsCutsAOD(Double_t ptmin, Double_t ptmax, Double_t etamin, Double_t etamax){
    fPtMin  = ptmin;  fPtMax  = ptmax;
    fEtaMin = etamin; fEtaMax = etamax;
  }
    
  void SetPtCutsCrossCorr(Double_t ptminTrig, Double_t ptmaxTrig, Double_t ptminAssoc, Double_t ptmaxAssoc){
    fPtCutsCrossCorr = kTRUE;
    fPtMinTrig  = ptminTrig;  fPtMaxTrig  = ptmaxTrig;
    fPtMinAssoc  = ptminAssoc;  fPtMaxAssoc  = ptmaxAssoc;
  }

  void SetExtraDCACutsAOD(Double_t DCAxy, Double_t DCAz){
    fDCAxyCut  = DCAxy;
    fDCAzCut = DCAz;
  }

  void SetExtraTPCCutsAOD(Double_t maxTPCchi2, Int_t minNClustersTPC, Int_t minNTPCCrossedRows, Float_t minNTPCFindableCls){
    fTPCchi2Cut      = maxTPCchi2;
    fNClustersTPCCut = minNClustersTPC;
    fMinTPCCrossedRows = minNTPCCrossedRows;
    fMinTPCRowsOverFindableCls =  minNTPCFindableCls;
  }

   void SetExtraTPCCutsSharedAOD(Int_t minTPCsharedCut){
    fTPCsharedCut = minTPCsharedCut;
  }

   void SetSphericityCut(Double_t sphericityMin, Double_t sphericityMax) {
     fSphericityMin = sphericityMin;
     fSphericityMax = sphericityMax;
     fUseSphericityCut = kTRUE;
   }
   
  //==============MC analysis==============//
  void SetKinematicsCutsMC(Double_t ptmin, Double_t ptmax,
                           Double_t etamin, Double_t etamax){
    fPtMin  = ptmin; fPtMax  = ptmax;
    fEtaMin = etamin; fEtaMax = etamax;
  }
  void UseFlowAfterBurner(TF1 *gDifferentialV2) {
    fDifferentialV2 = gDifferentialV2;
    fUseFlowAfterBurner = kTRUE;
  }

  void IncludeSecondariesInMCgen()  {fIncludeSecondariesInMCgen = kTRUE;}
  void ExcludeSecondariesInMC()  {fExcludeSecondariesInMC = kTRUE;}
  void ExcludeWeakDecaysInMC() {fExcludeWeakDecaysInMC = kTRUE;}
  void ExcludeResonancesInMC() {fExcludeResonancesInMC = kTRUE;}
  void ExcludeElectronsInMC()  {fExcludeElectronsInMC = kTRUE;}
  void ExcludeParticlesExtra() {fExcludeParticlesExtra = kTRUE;}
  void ExcludeResonancePDGInMC(Double_t pdgValue) {fExcludeResonancePDGInMC = pdgValue;}
  void IncludeResonancePDGInMC(Double_t pdgValue) {fIncludeResonancePDGInMC = pdgValue;}


  void ExcludeResonancesLabelCut(Int_t gPdgResonanceCode) {
     fExcludeResonancesLabel = kTRUE;
     fMotherPDGCodeToExclude = gPdgResonanceCode;
  }

  void SetPDGCode(Int_t gPdgCodeTrig, Int_t gPdgCodeAssoc, Bool_t setCrossCorr) {
    fUseMCPdgCode = kTRUE;
    fPDGCodeToBeAnalyzedTrig = gPdgCodeTrig;
    fPDGCodeToBeAnalyzedAssoc = gPdgCodeAssoc;
    fUseRapidity = kTRUE;
    fUsePIDMC = kTRUE;
    fCrossCorr = setCrossCorr;
    }
    
   void SetRejectInjectedSignals() {fExcludeInjectedSignals = kTRUE;}

   void SetRejectInjectedSignalsGenName(TString genToBeKept) {
    fGenToBeKept = genToBeKept; 
    fRejectCheckGenName=kTRUE;
    fExcludeInjectedSignals = kTRUE;
  }

  void SetUseNUADeep() {
    fUseNUADeep = kTRUE;
  }

  void SetUseRaaGeoCut(Float_t deadZoneWidth = 3, Float_t cutGeoNcrNclLength = 130, Float_t cutGeoNcrNclGeom1Pt = 1.5, Float_t cutGeoNcrNclFractionNcr = 0.85, Float_t cutGeoNcrNclFractionNcl = 0.7){
    fUseRaaGeoCut=kTRUE;
    fDeadZoneWidth = deadZoneWidth; 
    fCutGeoNcrNclLength = cutGeoNcrNclLength;
    fCutGeoNcrNclGeom1Pt = cutGeoNcrNclGeom1Pt;
    fCutGeoNcrNclFractionNcr = cutGeoNcrNclFractionNcr;
    fCutGeoNcrNclFractionNcl = cutGeoNcrNclFractionNcl;
  }

  void SetTimeRangeCutPbPb2018(){
    fUseTimeRangeCutForPbPb2018 = kTRUE;
  }

  void SetTOFBCPileUpCut(){
    fUseTOFBCPileUpCut = kTRUE;
  }

  void SetTPCInOutRowsCut(Int_t innermostRows = 2,  Int_t outermostRows = 20){
    fUseTPCInOutRowsCut = kTRUE;
    fInRows = innermostRows;
    fOutRows = outermostRows; 
  }
  
  //Centrality
  void SetCentralityEstimator(const char* centralityEstimator) {fCentralityEstimator = centralityEstimator;}
  const char* GetCentralityEstimator(void)  const              {return fCentralityEstimator;}
  void SetCentralityPercentileRange(Double_t min, Double_t max) { 
    fUseCentrality = kTRUE;
    fCentralityPercentileMin=min;
    fCentralityPercentileMax=max;
  }
  void SetImpactParameterRange(Double_t min, Double_t max) { 
    fUseCentrality = kTRUE;
    fImpactParameterMin=min;
    fImpactParameterMax=max;
  }
  void SetUseMultSelectionFramework(Bool_t multFramework) { fUseMultSelectionFramework = multFramework;} 
  void SetUseUncheckedCentrality(Bool_t unchecked) { fUseUncheckedCentrality = unchecked; }


  //multiplicity
  void SetMultiplicityEstimator(const char* multiplicityEstimator) {fMultiplicityEstimator = multiplicityEstimator;}
  const char* GetMultiplicityEstimator(void)  const              {return fMultiplicityEstimator;}
  void SetMultiplicityRange(Double_t min, Double_t max) {
    fUseMultiplicity = kTRUE;
    fNumberOfAcceptedTracksMin = min;
    fNumberOfAcceptedTracksMax = max;}

  //percentile
  void SetPercentileRange(Double_t min, Double_t max) { 
    fCentralityPercentileMin=min;
    fCentralityPercentileMax=max;
  }
  
  // additional event cuts (default = kFALSE)
  void UseOfflineTrigger() {fUseOfflineTrigger = kTRUE;}
  void CheckFirstEventInChunk() {fCheckFirstEventInChunk = kTRUE;}
  void CheckPileUp() {fCheckPileUp = kTRUE;}
  void UseSPDPileUpCuts(){fUsePileUpSPD = kTRUE;}

  void SetPileUpSPDParams(Int_t minVtxPileUpContrSPD, Float_t minPileUpZdistSPD){
    fModifySPDDefaultParams = kTRUE;
    fMinVtxPileUpContrSPD = minVtxPileUpContrSPD;
    fMinPileUpZdistSPD = minPileUpZdistSPD;
  }  
  void CheckPrimaryFlagAOD() {fCheckPrimaryFlagAOD = kTRUE;}
  void UseMCforKinematics() {fUseMCforKinematics = kTRUE;}
  void SetRebinnedCorrHistos() {fRebinCorrHistos = kTRUE;}
  void SetCentralityWeights(TH1* hist) { fCentralityWeights = hist; }
  Bool_t AcceptEventCentralityWeight(Double_t centrality);

  void SetUseAdditionalVtxCuts(Bool_t useAdditionalVtxCuts) {
    fUseAdditionalVtxCuts=useAdditionalVtxCuts;}

  void SetUseOutOfBunchPileUpCutsLHC15o(Bool_t useOutOfBunchPileUpCuts, Float_t slope=3.38, Float_t offset=15000) {
    fCheckOutOfBunchPileUp = kTRUE;
    fUseOOBPileUpCutsLHC15o = useOutOfBunchPileUpCuts;
    fPileupLHC15oSlope = slope;
    fPileupLHC15oOffset = offset;
  }
  
  void SetUseOutOfBunchPileUpCutsLHC15oJpsi(Bool_t useOutOfBunchPileUpCutsJpsi){
    fCheckOutOfBunchPileUp = kTRUE;
    fUseOOBPileUpCutsLHC15oJpsi = useOutOfBunchPileUpCutsJpsi;
  }

  void SetUseOutOfBunchPileUpCutsLHC18onTPCclus(Bool_t useOutOfBunchPileUpCutsnTPCclus, Float_t slope = 2000.0, Float_t param1 = 0.013, Float_t param2 = 1.25e-9){
    fCheckOutOfBunchPileUp = kTRUE;
    fUseOOBPileUpCutsLHC18nTPCclus = useOutOfBunchPileUpCutsnTPCclus;
    fOOBLHC18Slope = slope;
    fOOBLHC18Par1 = param1;
    fOOBLHC18Par2 = param2;
  }
  
  void SetUseDetailedTrackQA(Bool_t useDetailedTracksQA) {
    fDetailedTracksQA=useDetailedTracksQA;}
  
  // function to exclude the weak decay products
  Bool_t IsThisAWeakDecayingParticle(TParticle *thisGuy);
    
  //Acceptance filter
  void SetAcceptanceParameterization(TF1 *parameterization) {
    fAcceptanceParameterization = parameterization;}

  //pid
  enum kDetectorUsedForPID { kTPCpid, kTOFpid, kTPCTOF, kTPCTOFreq, kTPCTOFexcl }; // default TPC & TOF pid (via GetTPCpid & GetTOFpid)  
  enum kParticleOfInterest { kMuon, kElectron, kPion, kKaon, kProton };

  void SetUseBayesianPID(Double_t gMinProbabilityValue) {
    fUsePID = kTRUE; fUsePIDnSigma = kFALSE; fUsePIDPropabilities = kTRUE;
    fMinAcceptedPIDProbability = gMinProbabilityValue; }
  
  void SetUseNSigmaPID(Int_t gMaxNSigmaAcc, Int_t gMaxNSigmaExcl, Float_t pidMomCut) {
    fUsePID = kTRUE; fUsePIDPropabilities = kFALSE; fUsePIDnSigma = kTRUE;
    fPIDNSigmaAcc = gMaxNSigmaAcc;
    fPIDNSigmaExcl = gMaxNSigmaExcl;
    fPIDMomCut = pidMomCut;
    fUseRapidity = kTRUE;
  } // nsigma values for PID acceptance and rejection and pT threshold to move from TPC only and TPC+TOF for both methods: Bayes and nSigma Combined. usually 0.6 for pi and p and 0.4 for K.
    
  void SetUseNSigmaPIDKaons() {
    fUsePIDKaons = kTRUE; fUsePIDPropabilities = kFALSE;
    fUseRapidity = kTRUE;
  }
  
  void SetDetectorUsedForPID(kDetectorUsedForPID detConfig) {
    fPidDetectorConfig = detConfig;}
  void SetEventClass(TString receivedEventClass){
    fEventClass = receivedEventClass;
  }

  void SetCustomBinning(TString receivedCustomBinning) { fCustomBinning = receivedCustomBinning; }


    // electron rejection
    void SetElectronRejection(Double_t gMaxNSigma){
      fElectronRejection = kTRUE;
      fElectronRejectionNSigma = gMaxNSigma;
    }

    void SetElectronOnlyRejection(Double_t gMaxNSigma){
      fElectronRejection       = kTRUE;
      fElectronOnlyRejection   = kTRUE;
      fElectronRejectionNSigma = gMaxNSigma;
    }

    void SetElectronRejectionPt(Double_t minPt,Double_t maxPt){
      fElectronRejectionMinPt  = minPt;
      fElectronRejectionMaxPt  = maxPt;
    }

    void SetVZEROCalibrationFile(const char* filename, const char* lhcPeriod);
    void SetParticleOfInterest(AliPID::EParticleType trig, AliPID::EParticleType assoc, Bool_t setCrossCorr);
    

 private:
  Double_t    IsEventAccepted(AliVEvent* event);
  Double_t    GetRefMultiOrCentrality(AliVEvent* event);
  Double_t    GetReferenceMultiplicityFromAOD(AliVEvent* event);
  Double_t    GetEventPlane(AliVEvent* event);
  //===============================correction
  Double_t    GetTrackbyTrackCorrectionMatrix(Double_t vEta, 
					      Double_t vPhi, 
					      Double_t vPt, 
					      Short_t vCharge, 
					      Double_t gCentrality);
  //===============================correction
  TObjArray* GetAcceptedTracks(AliVEvent* event, Double_t gCentrality, Double_t gReactionPlane, Double_t &gSphericity, Int_t &nAcceptedTracksAboveHighPtThreshold);
  TObjArray* GetShuffledTracks(TObjArray* tracks, Double_t gCentrality, AliVEvent *event);

  Double_t GetChannelEqualizationFactor(Int_t run, Int_t channel);
  Double_t GetEqualizationFactor(Int_t run, const char *side);
 
  Bool_t fDebugLevel; // debug level

  TClonesArray* fArrayMC; //! AOD object  //+++++++++++++++++++++
  AliBalancePsi *fBalance; //BF object
  AliBalanceEbyE *fBalanceEbyE; //BF object (EbyE)
  Bool_t fRunShuffling;//run shuffling or not
  AliBalancePsi *fShuffledBalance; //BF object (shuffled)
  Bool_t fRunMixing;//run mixing or not
  Bool_t fRunMixingEventPlane;//run mixing with Event Plane
  Bool_t fRunEbyE;//run balance function on an event-by-event basis
  Int_t  fMixingTracks;
  Int_t  fMaxNbMixedEvents; 
  AliBalancePsi *fMixedBalance; //BF object (mixed)
  AliEventPoolManager*     fPoolMgr;         //! event pool manager

  TList *fList; //fList object
  TList *fListCrossCorr; //fList object
  TList *fListBF; //fList object
  TList *fListBFS; //fList object
  TList *fListBFM; //fList object
  TList *fHistListPIDQA;  //! list of histograms

  TList *fListNUA;  //fList of TH3F for NUA run-by-run corrections
  TList *fListNUE;   //fList of TH1F for NUE run-by-run corrections

  AliAnalysisTaskSignedBF::eCorrProcedure fCorrProcedure;

  //defualt kFALSE to be switch on for old correction method
  
  TH2F *fHistEventStats; //event stats
  TH2F *fHistCentStats; //centrality stats
  TH2F *fHistCentStatsUsed; //centrality stats USED
  TH1F *fHistTriggerStats; //trigger stats
  TH1F *fHistTrackStats; //Track filter bit stats
  TH1F *fHistVx; //x coordinate of the primary vertex
  TH1F *fHistVy; //y coordinate of the primary vertex
  TH2F *fHistVz; //z coordinate of the primary vertex
  TH1F *fHistCentrAfterEventSel; //event centrality distribution after all event selection

  TH2F *fHistMixEvents; //number of events that is mixed with in the current pool
  TH2F *fHistMixTracks; //number of tracks that is mixed with in the current pool

  TH2F *fHistTPCvsVZEROMultiplicity; //VZERO vs TPC reference multiplicity
  TH2F *fHistCL1vsVZEROPercentile; //VZERO vs TPC centrality to be used to monitor pileup 2015 data
  TH2F *fHistVZEROSignal; //VZERO channel vs signal

  TH2F *fHistEventPlane; //event plane distribution

  TH2F *fHistClus;//number of clusters (QA histogram)
  TH2F *fHistDCA;//DCA  (QA histogram)
  TH2F *fHistChi2;//track chi2 (QA histogram)
  TH2F *fHistPt;//transverse momentum (QA histogram)
  TH2F *fHistPtTrig;//transverse momentum for triggers (QA cross correlations histogram)
  TH2F *fHistPtAssoc;//transverse momentum for associated (QA cross correlations histogram)
  TH2F *fHistPtCorr;//transverse momentum after Corrrection (QA histogram)
  TH2F *fHistPtCorrTrig;//transverse momentum after Corrrection for triggers (QA cross correlations histogram)
  TH2F *fHistPtCorrAssoc;//transverse momentum after Corrrection for associated (QA cross correlations histogram)
  TH2F *fHistEta;//pseudorapidity (QA histogram)
  TH2F *fHistEtaCorr;//pseudorapidity after correction (QA histogram)
  TH2F *fHistRapidity;//rapidity (QA histogram)
  TH2F *fHistRapidityTrig;//rapidity for triggers (QA cross correlations histogram)
  TH2F *fHistRapidityAssoc;//rapidity for associated (QA cross correlations histogram)
  TH2F *fHistRapidityCorr;//rapidity after correction (QA histogram)
  TH2F *fHistRapidityCorrTrig;//rapidity after correction for triggers (QA cross correlations histogram)
  TH2F *fHistRapidityCorrAssoc;//rapidity after correction for associated (QA cross correlations histogram)
  TH2F *fHistPhi;//phi (QA histogram)
  TH2F *fHistPhiTrig;//phi for triggers (QA cross correlations histogram)
  TH2F *fHistPhiAssoc;//phi for associated (QA cross correlations histogram)
  TH2F *fHistPhiCorr;//phi after correction (QA histogram)
  TH2F *fHistPhiCorrTrig;//phi after correction for triggers (QA cross correlations histogram)
  TH2F *fHistPhiCorrAssoc;//phi after correction for associated (QA cross correlationshistogram)
  TH3F *fHistEtaVzPos;//eta vs Vz pos particles (QA histogram)
  TH3F *fHistEtaVzPosCorr;//eta vs Vz pos particles after correction(QA histogram) 
  TH3F *fHistEtaVzNeg;//eta vs Vz neg particles (QA histogram)
  TH3F *fHistEtaVzNegCorr;//eta vs Vz neg particles (QA histogram)
  TH3F *fHistEtaPhiPos;//eta-phi pos particles (QA histogram)
  TH3F *fHistEtaPhiPosCorr;//eta-phi pos particles after corrections  (QA histogram) 
  TH3F *fHistEtaPhiNeg;//eta-phi neg particles (QA histogram)
  TH3F *fHistEtaPhiNegCorr;//eta-phi neg particles after corrections (QA histogram)
  TH3F *fHistEtaPhiVzPlus;//eta-phi-Vz pos particles (QA histogram)
  TH3F *fHistEtaPhiVzMinus;//eta-phi-Vz neg particles (QA histogram)
  TH3F *fHistEtaPhiVzPlusCorr;//eta-phi-Vz pos particles after corrections  (QA histogram)
  TH3F *fHistEtaPhiVzMinusCorr;//eta-phi-Vz neg particles after corrections (QA histogram)
  TH2F *fHistPhiBefore;//phi before v2 afterburner (QA histogram)
  TH2F *fHistPhiAfter;//phi after v2 afterburner (QA histogram)
  TH2F *fHistPhiPos;//phi for positive particles (QA histogram)
  TH2F *fHistPhiNeg;//phi for negative particles (QA histogram)
  TH2F *fHistV0M;//V0 multiplicities (QA histogram)
  TH2F *fHistRefTracks;//reference track multiplicities (QA histogram)
  TH2F *fHistPhivZ;//phi vs Vz (QA histos) 
  TH2F *fHistEtavZ;//eta vs Vz (QA histos)
  TH2F *fHistPtPhi;//pt vs phi for GeOCut PbPb2018
  TH1F *fHistPdgMC;
  TH1F *fHistPdgMCAODrec;//pdg code of accepted tracks in MCAODrec
  TH1F *fHistSphericity; //sphericity of accepted tracks
  TH2F *fHistMultiplicityVsSphericity; //multiplicity vs sphericity of accepted tracks
  TH2F *fHistMeanPtVsSphericity; //mean pT vs sphericity of accepted tracks
  TH1F *fHistSphericityAfter; //sphericity of accepted tracks
  TH2F *fHistMultiplicityVsSphericityAfter; //multiplicity vs sphericity of accepted tracks
  TH2F *fHistMeanPtVsSphericityAfter; //mean pT vs sphericity of accepted tracks
  TH2F *fHistPhiNUADeep;
    
  //============PID============//
  TH2D *fHistdEdxVsPTPCbeforePID;//TPC dEdx vs momentum before PID cuts (QA histogram)
  TH2D *fHistBetavsPTOFbeforePID;//beta vs momentum before PID cuts (QA histogram)
  TH2D *fHistProbTPCvsPtbeforePID; //TPC probability vs pT before PID cuts (QA histogram)
  TH2D *fHistProbTOFvsPtbeforePID;//TOF probability vs pT before PID cuts (QA histogram)
  TH2D *fHistProbTPCTOFvsPtbeforePID;//TOF/TPC probability vs pT before PID cuts (QA histogram)
  TH2D *fHistNSigmaTPCvsPtbeforePID;//TPC nsigma vs pT before PID cuts (QA histogram)
  TH2D *fHistNSigmaTOFvsPtbeforePID;//TOF nsigma vs pT before PID cuts (QA histogram)
  TH2D *fHistBetaVsdEdXbeforePID;//TPCTOF  before PID cuts (QA histogram)
  TH2D *fHistNSigmaTPCTOFvsPtbeforePID;//TPCTOF  before PID cuts (QA histogram)
  TH3D *fHistNSigmaTPCTOFPbefPID;//+++++++++++++++

  TH2D *fHistdEdxVsPTPCafterPID;//TPC dEdx vs momentum after PID cuts (QA histogram)
  TH2D *fHistdEdxVsPTPCafterPIDTrig;//TPC dEdx vs momentum after PID cuts for triggers (QA cross correlations histogram)
  TH2D *fHistdEdxVsPTPCafterPIDAssoc;//TPC dEdx vs momentum after PID cuts for associated (QA cross correlations histogram)
  TH2D *fHistBetavsPTOFafterPID;//beta vs momentum after PID cuts (QA histogram)
  TH2D *fHistBetavsPTOFafterPIDTrig;//beta vs momentum after PID cuts for triggers (QA cross correlations histogram)
  TH2D *fHistBetavsPTOFafterPIDAssoc;//beta vs momentum after PID cuts for associated (QA cross correlations histogram)
  TH2D *fHistProbTPCvsPtafterPID; //TPC probability vs pT after PID cuts (QA histogram)
  TH2D *fHistProbTOFvsPtafterPID;//TOF probability vs pT after PID cuts (QA histogram)
  TH2D *fHistProbTPCTOFvsPtafterPID;//TOF/TPC probability vs pT after PID cuts (QA histogram)
  TH2D *fHistNSigmaTPCvsPtafterPID;//TPC nsigma vs pT after PID cuts (QA histogram)
  TH2D *fHistNSigmaTOFvsPtafterPID;//TOF nsigma vs pT after PID cuts (QA histogram)
  TH2D *fHistBetaVsdEdXafterPID;//TPCTOF  before PID cuts (QA histogram)
  TH2D *fHistNSigmaTPCTOFvsPtafterPID;//TPCTOF  before PID cuts (QA histogram)
  TH3D *fHistNSigmaTPCTOFPafterPID; //++++++++++++++++++

  TH2D *fHistdEdxVsPTPCbeforePIDelectron; //!
  TH2D *fHistNSigmaTPCvsPtbeforePIDelectron; //!
  TH2D *fHistdEdxVsPTPCafterPIDelectron; //!
  TH2D *fHistNSigmaTPCvsPtafterPIDelectron; //!
  
  TH3F *fHistCorrectionPlus[kCENTRALITY]; //====correction
  TH3F *fHistCorrectionMinus[kCENTRALITY]; //===correction

  TH3F *fHistNUACorrPlus[kNBRUN]; //====correction
  TH3F *fHistNUACorrMinus[kNBRUN]; //===correction

  Int_t fRunNb[kNBRUN]; //====correction
  Int_t fTotalNbRun; //total number of runs used in the analysis

  TH1F *fHistpTCorrPlus[kCENTRALITY]; //====correction
  TH1F *fHistpTCorrMinus[kCENTRALITY]; //===correction
  TH1F *fHistpTCorrPlusTrig[kCENTRALITY]; //====correction for triggers (cross correlations)
  TH1F *fHistpTCorrMinusTrig[kCENTRALITY]; //===correction for triggers (cross correlations)
  TH1F *fHistpTCorrPlusAssoc[kCENTRALITY]; //====correction for associated (cross correlations)
  TH1F *fHistpTCorrMinusAssoc[kCENTRALITY]; //===correction for associated (cross correlations)
  
  Double_t fCentralityArrayForCorrections[kCENTRALITY];
  Int_t fCentralityArrayBinsForCorrections;

  TH1* fCentralityWeights;		     // for centrality flattening

  AliPIDResponse *fPIDResponse;     //! PID response object
  AliPIDCombined       *fPIDCombined;     //! combined PID object
  
  AliPID::EParticleType fParticleOfInterest[2];//analyzed particle
  kDetectorUsedForPID   fPidDetectorConfig;//used detector for PID
  Double_t fMassParticleOfInterest[2];//particle mass (for rapidity calculation)

    
  Bool_t fUsePID; //flag to use PID
  Bool_t fUsePIDKaons; //flag to use PID trial for kaons
  Bool_t fUsePIDMC; //flag to use PID in generated MC ("MC" analysis option)
  Bool_t fUsePIDnSigma;//flag to use nsigma method for PID
  Bool_t fUsePIDPropabilities;//flag to use probability method for PID
  Bool_t fUseRapidity;//flag to use rapidity instead of pseudorapidity in correlation histograms
  Bool_t fCrossCorr;//cross correlations or not
  Bool_t fPtCutsCrossCorr;//flag to use different cuts on pt for trigger and associated particle in case of cross correlations
  Double_t fPIDNSigmaAcc;//nsigma cut for PID acceptance
  Double_t fPIDNSigmaExcl;//nsigma cut for exclusive PID
  Double_t fMinAcceptedPIDProbability;//probability cut for PID

  Bool_t   fElectronRejection;//flag to use electron rejection
  Bool_t   fElectronOnlyRejection;//flag to use electron rejection with exclusive electron PID (no other particle in nsigma range)
  Double_t fElectronRejectionNSigma;//nsigma cut for electron rejection
  Double_t fElectronRejectionMinPt;//minimum pt for electron rejection (default = 0.)
  Double_t fElectronRejectionMaxPt;//maximum pt for electron rejection (default = 1000.)

  Double_t fPIDMomCut; // pT value from which we switche from TPC only to TPC TOD PID (both nsigma and Bayes) 
  
  //============PID============//

  AliESDtrackCuts *fESDtrackCuts; //ESD track cuts

  TString fCentralityEstimator;      //"V0M","TRK","TKL","ZDC","FMD"
  Bool_t fUseCentrality;//use the centrality (PbPb) or not (pp)
  Bool_t fUseMultSelectionFramework;//use the AliMultSelection framework; default: kFALSE
  Bool_t fUseUncheckedCentrality; // use unchecked centrality; default: kFALSE
  Double_t fCentralityPercentileMin;//centrality percentile min
  Double_t fCentralityPercentileMax;//centrality percentile max
  Double_t fImpactParameterMin;//impact parameter min (used for MC)
  Double_t fImpactParameterMax;//impact parameter max (used for MC)  
  TString fMultiplicityEstimator;//"V0M","V0A","V0C","TPC"
  Bool_t fUseMultiplicity;//use the multiplicity cuts
  Double_t fNumberOfAcceptedTracksMin;//min. number of number of accepted tracks (used for the multiplicity dependence study - pp)
  Double_t fNumberOfAcceptedTracksMax;//max. number of number of accepted tracks (used for the multiplicity dependence study - pp)
  TH2F *fHistNumberOfAcceptedTracks;//histo to store the number of accepted tracks
  TH1F *fHistMultiplicity;//hisot to store the number of accepted tracks
  TH2F *fHistMultvsPercent;//hisot to store the multiplicity vs centrality percentile
    
  
  Bool_t fUseOfflineTrigger;//Usage of the offline trigger selection
  Bool_t fCheckFirstEventInChunk;//Usage of the "First Event in Chunk" check (not needed for new productions)
  Bool_t fCheckPileUp;//Usage of the "Pile-Up" event check
  Bool_t fUsePileUpSPD;//Usage of the pile-up rejection with SPD instead of MultiVertexer one
  Bool_t fCheckPrimaryFlagAOD;// Usage of check on AliAODtrack::kPrimary (default = OFF)
  Bool_t fUseMCforKinematics;//Usage of MC information for filling the kinematics information of particles (only in MCAODrec mode)
  Bool_t fRebinCorrHistos;//Rebinning of corrected plots
  Bool_t fUseAdditionalVtxCuts;//usage of additional clean up cuts for primary vertex.
  Bool_t fCheckOutOfBunchPileUp; //default kFALSE! 
  Bool_t fUseOOBPileUpCutsLHC15o;//usage of correlation cuts to exclude out of bunche pile up. To be used for 2015 PbPb data. multEsd - fPileupLHC15oSlope*multTPC) > fPileupLHC15oOffset
  Float_t fPileupLHC15oSlope; //parameters for LHC15o pile-up rejection  default: slope=3.35, offset 15000
  Float_t fPileupLHC15oOffset;
  
  Bool_t fUseOOBPileUpCutsLHC15oJpsi;//multVZERO < (-2200 + 2.5*ntrkTPCout + 1.2e-5*ntrkTPCout*ntrkTPCout
  
  Bool_t fUseOOBPileUpCutsLHC18nTPCclus; //multVZERO < (-fOOBLHC18Slope + fOOBLHC18Par1*nTPCclus + fOOBLHC18Par2*nTPCclus*nTPCclus 
  Float_t fOOBLHC18Slope;
  Float_t fOOBLHC18Par1;
  Float_t fOOBLHC18Par2;

  Bool_t  fModifySPDDefaultParams;
  Int_t   fMinVtxPileUpContrSPD;
  Float_t fMinPileUpZdistSPD;

  Bool_t fUseTOFBCPileUpCut;
  
  Bool_t fUseTPCInOutRowsCut;
  Int_t fInRows;
  Int_t fOutRows; 

  Bool_t fDetailedTracksQA; //fill Eta, Phi vs Vx histos to be used to check ME pools. 

  Double_t fVxMax;//vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax
  
  Bool_t fRequireHighPtTrigger;//use pT trigger
  Double_t fPtTriggerMin;//pT trigger min
  TH2F *fHistPtTriggerThreshold;//QA histo
  
  Int_t fnAODtrackCutBit;//track cut bit from track selection (only used for AODs)

  Bool_t fUseRaaGeoCut; //flag to switch on GeoCut for 2018PbPb data pass1
  Float_t fDeadZoneWidth; //parameters of the cut as implemented in AliESDtrackCuts.h, default values implemented as suggested by DPG and D mesons analysis
  Float_t fCutGeoNcrNclLength;
  Float_t fCutGeoNcrNclGeom1Pt;
  Float_t fCutGeoNcrNclFractionNcr;
  Float_t fCutGeoNcrNclFractionNcl;

  Float_t fUseTimeRangeCutForPbPb2018; 
  AliTimeRangeCut fTimeRangeCut;
 

  Double_t fPtMin;//only used for AODs
  Double_t fPtMax;//only used for AODs
  Double_t fPtMinTrig;//only used for AODs
  Double_t fPtMinAssoc;//only used for AODs
  Double_t fPtMaxTrig;//only used for AODs
  Double_t fPtMaxAssoc;//only used for AODs
  Double_t fEtaMin;//only used for AODs
  Double_t fEtaMax;//only used for AODs
  Double_t fPhiMin;//only used for AODs
  Double_t fPhiMax;//only used for AODs 

  Double_t fDCAxyCut;//only used for AODs
  Double_t fDCAzCut;//only used for AODs

  Double_t fTPCchi2Cut;//only used for AODs
  Int_t fNClustersTPCCut;//only used for AODs
  Int_t fMinTPCCrossedRows; // only used for AODs
  Float_t fMinTPCRowsOverFindableCls; //only used for AODs
  Int_t fTPCsharedCut;//only used for AODs

  Double_t fSphericityMin;//sphericity min cut (currently only for AODs)
  Double_t fSphericityMax;//sphericity max cut (currently only for AODs)
  Bool_t fUseSphericityCut;//sphericity cut (currently only for AODs)

  TF1 *fAcceptanceParameterization;//acceptance filter used for MC

  TF1 *fDifferentialV2;//pt-differential v2 (from real data)
  Bool_t fUseFlowAfterBurner;//Usage of a flow after burner
  Bool_t fUseNUADeep;//Usage of a deep in phi

  Bool_t fIncludeSecondariesInMCgen;//flag to include the secondaries from material and weak decays in the MC analysis (needed for fIncludeResonancePDGInMC)
  Bool_t fExcludeSecondariesInMC;//flag to exclude the secondaries from material and weak decays in the MCAODrec analysis
  Bool_t fExcludeWeakDecaysInMC;//flag to exclude the weak decay products (if not done by IsPhysicalPrimary) from the MC analysis
  Bool_t fExcludeResonancesInMC;//flag to exclude the resonances' decay products (and conversion) from the MC analysis
  Bool_t fExcludeResonancesLabel;//flag to exclude the resonances using mother's label;
  Bool_t fExcludeElectronsInMC;//flag to exclude the electrons from the MC analysis
  Bool_t fExcludeParticlesExtra;//flag to exclude particles from the MC analysis (extra)
  Bool_t fUseMCPdgCode; //Boolean to analyze a set of particles in MC and MCAODrec
  Int_t fPDGCodeToBeAnalyzedTrig;//Analyze a set of particles in MC and MCAODrec
  Int_t fPDGCodeToBeAnalyzedAssoc;//Analyze a set of particles in MC and MCAODrec
  Int_t fMotherPDGCodeToExclude; // exclude the resonance with this PDG with the label cut from the MC analysis
  Int_t fExcludeResonancePDGInMC;// exclude the resonance with this PDG from the MC analysis
  Int_t fIncludeResonancePDGInMC;// include excluvely this resonance with this PDG to the MC and MCAODrec analysis

  Bool_t fExcludeInjectedSignals; //Flag to reject MC injected signals from MC analysis
  Bool_t fRejectCheckGenName; // Flag to activate the injected signal rejection based on the name of the MC generator 
  TString fGenToBeKept; //String to select the generator name that has to be kept for analysis
  
  TString fEventClass; //Can be "EventPlane", "Centrality", "Multiplicity"
  TString fCustomBinning;//for setting customized binning (for output AliTHn of AliBalancePsi)
  
  //VZERO calibration
  TH1F *fHistVZEROAGainEqualizationMap;//VZERO calibration map
  TH1F *fHistVZEROCGainEqualizationMap;//VZERO calibration map
  TH2F *fHistVZEROChannelGainEqualizationMap; //VZERO calibration map

  TH2F *fHistGlobalvsESDBeforePileUpCuts; //histos to monitor Out of bunch pile up selection
  TH2F *fHistGlobalvsESDAfterPileUpCuts;

  TH2F *fHistV0MvsTPCoutBeforePileUpCuts; //histos to monitor pile up cuts J/psi
  TH2F *fHistV0MvsTPCoutAfterPileUpCuts;

  TH2F *fHistV0MvsnTPCclusBeforePileUpCuts; 
  TH2F *fHistV0MvsnTPCclusAfterPileUpCuts;

  TH1F *fHistCentrBeforePileUpCuts;
  TH1F *fHistCentrAfterPileUpCuts;
  
  //AliAnalysisUtils
  AliAnalysisUtils *fUtils;//AliAnalysisUtils

  AliAnalysisTaskSignedBF(const AliAnalysisTaskSignedBF&); // not implemented
  AliAnalysisTaskSignedBF& operator=(const AliAnalysisTaskSignedBF&); // not implemented
  
  ClassDef(AliAnalysisTaskSignedBF, 1); // example of analysis
};



#endif
