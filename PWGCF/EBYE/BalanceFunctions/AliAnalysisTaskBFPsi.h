#ifndef ALIANALYSISTASKBFPSI_H
#define ALIANALYSISTASKBFPSI_H

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


#include "AliAnalysisTaskSE.h"
#include "AliBalancePsi.h"

#include "AliPID.h"  
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
 
//================================correction
#define kCENTRALITY 101  
//const Double_t centralityArrayForPbPb[kCENTRALITY+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};
//const TString centralityArrayForPbPb_string[kCENTRALITY] = {"0-5","5-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80"};
//================================correction

class AliAnalysisTaskBFPsi : public AliAnalysisTaskSE {
 public:
  enum etriggerSel{kMB, kCentral, kINT7, kppHighMult};
  
  AliAnalysisTaskBFPsi(const char *name = "AliAnalysisTaskBFPsi");
  virtual ~AliAnalysisTaskBFPsi(); 
   
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   FinishTaskOutput();
  virtual void   Terminate(Option_t *);

  //========================correction
  virtual void   SetInputCorrection(TString filename, 
				    Int_t nCentralityBins, 
				    Double_t *centralityArrayForCorrections);
  //========================correction
  // void SetDebugLevel() {fDebugLevel = kTRUE;} //hides overloaded virtual function

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

  //==============AOD analysis==============//
  void SetAODtrackCutBit(Int_t bit){
    fnAODtrackCutBit = bit;
  }

  void SetKinematicsCutsAOD(Double_t ptmin, Double_t ptmax, Double_t etamin, Double_t etamax){
    fPtMin  = ptmin;  fPtMax  = ptmax;
    fEtaMin = etamin; fEtaMax = etamax;
  }

  void SetExtraDCACutsAOD(Double_t DCAxy, Double_t DCAz){
    fDCAxyCut  = DCAxy;
    fDCAzCut = DCAz;
  }

   void SetExtraTPCCutsAOD(Double_t maxTPCchi2, Int_t minNClustersTPC){
    fTPCchi2Cut      = maxTPCchi2;
    fNClustersTPCCut = minNClustersTPC;
  }

   void SetExtraTPCCutsSharedAOD(Int_t minTPCsharedCut){
    fTPCsharedCut = minTPCsharedCut;
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

  void SetPDGCode(Int_t gPdgCode) {
    fUseMCPdgCode = kTRUE;
    fPDGCodeToBeAnalyzed = gPdgCode;
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
  void CheckPrimaryFlagAOD() {fCheckPrimaryFlagAOD = kTRUE;}
  void UseMCforKinematics() {fUseMCforKinematics = kTRUE;}
  void SetCentralityWeights(TH1* hist) { fCentralityWeights = hist; }
  Bool_t AcceptEventCentralityWeight(Double_t centrality);

  void SetUseAdditionalVtxCuts(Bool_t useAdditionalVtxCuts) {
    fUseAdditionalVtxCuts=useAdditionalVtxCuts;}

  void SetUseOutOfBunchPileUpCutsLHC15o(Bool_t useOutOfBunchPileUpCuts) {
    fUseOutOfBunchPileUpCutsLHC15o=useOutOfBunchPileUpCuts;}
  
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

  void SetUseNSigmaPID(Double_t gMaxNSigma) {
    fUsePID = kTRUE; fUsePIDPropabilities = kFALSE; fUsePIDnSigma = kTRUE;
    fPIDNSigma = gMaxNSigma; }

  void SetDetectorUsedForPID(kDetectorUsedForPID detConfig) {
    fPidDetectorConfig = detConfig;}
  void SetEventClass(TString receivedEventClass){
    fEventClass = receivedEventClass;
  }

  void SetUseRapidity(Bool_t useRapidity = kTRUE){
    fUseRapidity = useRapidity;
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
    void SetParticleOfInterest(kParticleOfInterest poi);


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
  TObjArray* GetAcceptedTracks(AliVEvent* event, Double_t gCentrality, Double_t gReactionPlane);
  TObjArray* GetShuffledTracks(TObjArray* tracks, Double_t gCentrality);

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
  AliBalancePsi *fMixedBalance; //BF object (mixed)
  AliEventPoolManager*     fPoolMgr;         //! event pool manager

  TList *fList; //fList object
  TList *fListBF; //fList object
  TList *fListBFS; //fList object
  TList *fListBFM; //fList object
  TList *fHistListPIDQA;  //! list of histograms

  TH2F *fHistEventStats; //event stats
  TH2F *fHistCentStats; //centrality stats
  TH2F *fHistCentStatsUsed; //centrality stats USED
  TH1F *fHistTriggerStats; //trigger stats
  TH1F *fHistTrackStats; //Track filter bit stats
  TH1F *fHistVx; //x coordinate of the primary vertex
  TH1F *fHistVy; //y coordinate of the primary vertex
  TH2F *fHistVz; //z coordinate of the primary vertex

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
  TH2F *fHistEta;//pseudorapidity (QA histogram)
  TH2F *fHistRapidity;//rapidity (QA histogram)
  TH2F *fHistPhi;//phi (QA histogram)
  TH3F *fHistEtaPhiPos;//eta-phi pos particles (QA histogram) 		 	 
  TH3F *fHistEtaPhiNeg;//eta-phi neg particles (QA histogram)
  TH2F *fHistPhiBefore;//phi before v2 afterburner (QA histogram)
  TH2F *fHistPhiAfter;//phi after v2 afterburner (QA histogram)
  TH2F *fHistPhiPos;//phi for positive particles (QA histogram)
  TH2F *fHistPhiNeg;//phi for negative particles (QA histogram)
  TH2F *fHistV0M;//V0 multiplicities (QA histogram)
  TH2F *fHistRefTracks;//reference track multiplicities (QA histogram)

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
  TH2D *fHistBetavsPTOFafterPID;//beta vs momentum after PID cuts (QA histogram)
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
  Double_t fCentralityArrayForCorrections[kCENTRALITY];
  Int_t fCentralityArrayBinsForCorrections;

  TH1* fCentralityWeights;		     // for centrality flattening

  AliPIDResponse *fPIDResponse;     //! PID response object
  AliPIDCombined       *fPIDCombined;     //! combined PID object
  
  kParticleOfInterest  fParticleOfInterest;//analyzed particle
  kDetectorUsedForPID   fPidDetectorConfig;//used detector for PID
  Double_t fMassParticleOfInterest;//particle mass (for rapidity calculation) 

  Bool_t fUsePID; //flag to use PID 
  Bool_t fUsePIDnSigma;//flag to use nsigma method for PID
  Bool_t fUsePIDPropabilities;//flag to use probability method for PID
  Bool_t fUseRapidity;//flag to use rapidity instead of pseudorapidity in correlation histograms
  Double_t fPIDNSigma;//nsigma cut for PID
  Double_t fMinAcceptedPIDProbability;//probability cut for PID

  Bool_t   fElectronRejection;//flag to use electron rejection
  Bool_t   fElectronOnlyRejection;//flag to use electron rejection with exclusive electron PID (no other particle in nsigma range)
  Double_t fElectronRejectionNSigma;//nsigma cut for electron rejection
  Double_t fElectronRejectionMinPt;//minimum pt for electron rejection (default = 0.)
  Double_t fElectronRejectionMaxPt;//maximum pt for electron rejection (default = 1000.)
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
  TH2F *fHistNumberOfAcceptedTracks;//hisot to store the number of accepted tracks
  TH1F *fHistMultiplicity;//hisot to store the number of accepted tracks
  TH2F *fHistMultvsPercent;//hisot to store the multiplicity vs centrality percentile
    
  
  Bool_t fUseOfflineTrigger;//Usage of the offline trigger selection
  Bool_t fCheckFirstEventInChunk;//Usage of the "First Event in Chunk" check (not needed for new productions)
  Bool_t fCheckPileUp;//Usage of the "Pile-Up" event check
  Bool_t fCheckPrimaryFlagAOD;// Usage of check on AliAODtrack::kPrimary (default = OFF)
  Bool_t fUseMCforKinematics;//Usage of MC information for filling the kinematics information of particles (only in MCAODrec mode)

  Bool_t fUseAdditionalVtxCuts;//usage of additional clean up cuts for primary vertex.

  Bool_t fUseOutOfBunchPileUpCutsLHC15o;//usage of correlation cuts to exclude out of bunche pile up. To be used for 2015 PbPb data. 

  Double_t fVxMax;//vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax

  Int_t fnAODtrackCutBit;//track cut bit from track selection (only used for AODs)

  Double_t fPtMin;//only used for AODs
  Double_t fPtMax;//only used for AODs
  Double_t fEtaMin;//only used for AODs
  Double_t fEtaMax;//only used for AODs
  Double_t fPhiMin;//only used for AODs
  Double_t fPhiMax;//only used for AODs 

  Double_t fDCAxyCut;//only used for AODs
  Double_t fDCAzCut;//only used for AODs

  Double_t fTPCchi2Cut;//only used for AODs
  Int_t fNClustersTPCCut;//only used for AODs
  Int_t fTPCsharedCut;//only used for AODs

  TF1 *fAcceptanceParameterization;//acceptance filter used for MC

  TF1 *fDifferentialV2;//pt-differential v2 (from real data)
  Bool_t fUseFlowAfterBurner;//Usage of a flow after burner

  Bool_t fIncludeSecondariesInMCgen;//flag to include the secondaries from material and weak decays in the MC analysis (needed for fIncludeResonancePDGInMC)
  Bool_t fExcludeSecondariesInMC;//flag to exclude the secondaries from material and weak decays in the MCAODrec analysis
  Bool_t fExcludeWeakDecaysInMC;//flag to exclude the weak decay products (if not done by IsPhysicalPrimary) from the MC analysis
  Bool_t fExcludeResonancesInMC;//flag to exclude the resonances' decay products (and conversion) from the MC analysis
  Bool_t fExcludeElectronsInMC;//flag to exclude the electrons from the MC analysis
  Bool_t fExcludeParticlesExtra;//flag to exclude particles from the MC analysis (extra)
  Bool_t fUseMCPdgCode; //Boolean to analyze a set of particles in MC and MCAODrec
  Int_t fPDGCodeToBeAnalyzed; //Analyze a set of particles in MC and MCAODrec
  Int_t fExcludeResonancePDGInMC;// exclude the resonance with this PDG from the MC analysis
  Int_t fIncludeResonancePDGInMC;// include excluvely this resonance with this PDG to the MC and MCAODrec analysis
  TString fEventClass; //Can be "EventPlane", "Centrality", "Multiplicity"
  TString fCustomBinning;//for setting customized binning (for output AliTHn of AliBalancePsi)
  
  //VZERO calibration
  TH1F *fHistVZEROAGainEqualizationMap;//VZERO calibration map
  TH1F *fHistVZEROCGainEqualizationMap;//VZERO calibration map
  TH2F *fHistVZEROChannelGainEqualizationMap; //VZERO calibration map

  //AliAnalysisUtils
  AliAnalysisUtils *fUtils;//AliAnalysisUtils

  AliAnalysisTaskBFPsi(const AliAnalysisTaskBFPsi&); // not implemented
  AliAnalysisTaskBFPsi& operator=(const AliAnalysisTaskBFPsi&); // not implemented
  
  ClassDef(AliAnalysisTaskBFPsi, 9); // example of analysis
};



#endif
