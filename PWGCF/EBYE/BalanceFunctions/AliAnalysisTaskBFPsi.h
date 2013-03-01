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

class AliBalancePsi;
class AliESDtrackCuts;
class AliEventPoolManager;


#include "AliAnalysisTaskSE.h"
#include "AliBalancePsi.h"

#include "AliPID.h"  
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
 
//================================correction
#define kCENTRALITY 9  
const Double_t centralityArrayForPbPb[kCENTRALITY+1] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.};
//================================correction

class AliAnalysisTaskBFPsi : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskBFPsi(const char *name = "AliAnalysisTaskBFPsi");
  virtual ~AliAnalysisTaskBFPsi(); 
   
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   FinishTaskOutput();
  virtual void   Terminate(Option_t *);

  //========================correction
  virtual void   SetInputCorrection(const char* filename, const char* gSystem); 
  //========================correction

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
  void SetAnalysisCutObject(AliESDtrackCuts *const trackCuts) {
    fESDtrackCuts = trackCuts;}
  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {
    fVxMax = vx;
    fVyMax = vy;
    fVzMax = vz;
  }

  //==============AOD analysis==============//
  void SetAODtrackCutBit(Int_t bit){
    nAODtrackCutBit = bit;
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
  void ExcludeResonancesInMC() {fExcludeResonancesInMC = kTRUE;}

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

  //multiplicity
  void SetMultiplicityRange(Int_t min, Int_t max) {
    fUseMultiplicity = kTRUE;
    fNumberOfAcceptedTracksMin = min;
    fNumberOfAcceptedTracksMax = max;}
  
  void UseOfflineTrigger() {fUseOfflineTrigger = kTRUE;}
  
  //Acceptance filter
  void SetAcceptanceParameterization(TF1 *parameterization) {
    fAcceptanceParameterization = parameterization;}

  //pid
  enum kDetectorUsedForPID { kTPCpid, kTOFpid, kTPCTOF }; // default TPC & TOF pid (via GetTPCpid & GetTOFpid)  
  enum kParticleOfInterest { kMuon, kElectron, kPion, kKaon, kProton };  

  void SetUseBayesianPID(Double_t gMinProbabilityValue) {
    fUsePID = kTRUE; fUsePIDnSigma = kFALSE; fUsePIDPropabilities = kTRUE;
    fMinAcceptedPIDProbability = gMinProbabilityValue; }

  void SetUseNSigmaPID(Double_t gMaxNSigma) {
    fUsePID = kTRUE; fUsePIDPropabilities = kFALSE; fUsePIDnSigma = kTRUE;
    fPIDNSigma = gMaxNSigma; }

  void SetParticleOfInterest(kParticleOfInterest poi) {
    fParticleOfInterest = poi;}
  void SetDetectorUsedForPID(kDetectorUsedForPID detConfig) {
    fPidDetectorConfig = detConfig;}
    void SetEventClass(TString receivedEventClass){
        fEventClass = receivedEventClass;
    }

    //========================correction=============================//
    // void SetInputCorrection(const char* filename = "efficiencyCorrection.root"){}

 private:
  Double_t    IsEventAccepted(AliVEvent* event);
  Double_t    GetRefMultiOrCentrality(AliVEvent* event);
  Double_t    GetEventPlane(AliVEvent* event);
  //===============================correction
  Double_t    GetTrackbyTrackCorrectionMatrix(Double_t vEta, Double_t vPhi, Double_t vPt, Short_t vCharge, Double_t gCentrality);
  //===============================correction
  TObjArray* GetAcceptedTracks(AliVEvent* event, Double_t gCentrality, Double_t gReactionPlane);
  TObjArray* GetShuffledTracks(TObjArray* tracks, Double_t gCentrality);
  
  AliBalancePsi *fBalance; //BF object
  Bool_t fRunShuffling;//run shuffling or not
  AliBalancePsi *fShuffledBalance; //BF object (shuffled)
  Bool_t fRunMixing;//run mixing or not
  Bool_t fRunMixingEventPlane;//run mixing with Event Plane
  Int_t  fMixingTracks;
  AliBalancePsi *fMixedBalance; //TriggeredBF object (mixed)
  AliEventPoolManager*     fPoolMgr;         //! event pool manager

  TList *fList; //fList object
  TList *fListBF; //fList object
  TList *fListBFS; //fList object
  TList *fListBFM; //fList object
  TList *fHistListPIDQA;  //! list of histograms

  TH2F *fHistEventStats; //event stats
  TH2F *fHistCentStats; //centrality stats
  TH1F *fHistTriggerStats; //trigger stats
  TH1F *fHistTrackStats; //Track filter bit stats
  TH1F *fHistVx; //x coordinate of the primary vertex
  TH1F *fHistVy; //y coordinate of the primary vertex
  TH2F *fHistVz; //z coordinate of the primary vertex

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
  TH2D *fHistdEdxVsPTPCafterPID;//TPC dEdx vs momentum after PID cuts (QA histogram)
  TH2D *fHistBetavsPTOFafterPID;//beta vs momentum after PID cuts (QA histogram)
  TH2D *fHistProbTPCvsPtafterPID; //TPC probability vs pT after PID cuts (QA histogram)
  TH2D *fHistProbTOFvsPtafterPID;//TOF probability vs pT after PID cuts (QA histogram)
  TH2D *fHistProbTPCTOFvsPtafterPID;//TOF/TPC probability vs pT after PID cuts (QA histogram)
  TH2D *fHistNSigmaTPCvsPtafterPID;//TPC nsigma vs pT after PID cuts (QA histogram)
  TH2D *fHistNSigmaTOFvsPtafterPID;//TOF nsigma vs pT after PID cuts (QA histogram)
  
  TH3D *fHistMatrixCorrectionPlus[kCENTRALITY]; //======================================================correction
  TH3D *fHistMatrixCorrectionMinus[kCENTRALITY]; //=====================================================correction

  AliPIDResponse *fPIDResponse;     //! PID response object
  AliPIDCombined       *fPIDCombined;     //! combined PID object
  
  kParticleOfInterest  fParticleOfInterest;//analyzed particle
  kDetectorUsedForPID   fPidDetectorConfig;//used detector for PID

  Bool_t fUsePID; //flag to use PID 
  Bool_t fUsePIDnSigma;//flag to use nsigma method for PID
  Bool_t fUsePIDPropabilities;//flag to use probability method for PID
  Double_t fPIDNSigma;//nsigma cut for PID
  Double_t fMinAcceptedPIDProbability;//probability cut for PID
  //============PID============//

  AliESDtrackCuts *fESDtrackCuts; //ESD track cuts

  TString fCentralityEstimator;      //"V0M","TRK","TKL","ZDC","FMD"
  Bool_t fUseCentrality;//use the centrality (PbPb) or not (pp)
  Double_t fCentralityPercentileMin;//centrality percentile min
  Double_t fCentralityPercentileMax;//centrality percentile max
  Double_t fImpactParameterMin;//impact parameter min (used for MC)
  Double_t fImpactParameterMax;//impact parameter max (used for MC)

  Bool_t fUseMultiplicity;//use the multiplicity cuts
  Int_t fNumberOfAcceptedTracksMin;//min. number of number of accepted tracks (used for the multiplicity dependence study - pp)
  Int_t fNumberOfAcceptedTracksMax;//max. number of number of accepted tracks (used for the multiplicity dependence study - pp)
  TH2F *fHistNumberOfAcceptedTracks;//hisot to store the number of accepted tracks

  Bool_t fUseOfflineTrigger;//Usage of the offline trigger selection

  Double_t fVxMax;//vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax

  Int_t nAODtrackCutBit;//track cut bit from track selection (only used for AODs)

  Double_t fPtMin;//only used for AODs
  Double_t fPtMax;//only used for AODs
  Double_t fPtMinForCorrections;//only used for AODs
  Double_t fPtMaxForCorrections;//only used for AODs
  Double_t fPtBinForCorrections; //=================================correction
  Double_t fEtaMin;//only used for AODs
  Double_t fEtaMax;//only used for AODs
  Double_t fEtaMinForCorrections;//only used for AODs
  Double_t fEtaMaxForCorrections;//only used for AODs
  Double_t fEtaBinForCorrections; //=================================correction
  Double_t fPhiMin; //=================================correction 
  Double_t fPhiMax; //=================================correction
  Double_t fPhiMinForCorrections;//only used for AODs
  Double_t fPhiMaxForCorrections;//only used for AODs
  Double_t fPhiBinForCorrections; //=================================correction

  Double_t fDCAxyCut;//only used for AODs
  Double_t fDCAzCut;//only used for AODs

  Double_t fTPCchi2Cut;//only used for AODs
  Int_t fNClustersTPCCut;//only used for AODs

  TF1 *fAcceptanceParameterization;//acceptance filter used for MC

  TF1 *fDifferentialV2;//pt-differential v2 (from real data)
  Bool_t fUseFlowAfterBurner;//Usage of a flow after burner

  Bool_t fExcludeResonancesInMC;//flag to exclude the resonances' decay products from the MC analysis
  Bool_t fUseMCPdgCode; //Boolean to analyze a set of particles in MC
  Int_t fPDGCodeToBeAnalyzed; //Analyze a set of particles in MC
  TString fEventClass; //Can be "EventPlane", "Centrality", "Multiplicity"
  

  AliAnalysisTaskBFPsi(const AliAnalysisTaskBFPsi&); // not implemented
  AliAnalysisTaskBFPsi& operator=(const AliAnalysisTaskBFPsi&); // not implemented
  
  ClassDef(AliAnalysisTaskBFPsi, 5); // example of analysis
};



#endif
