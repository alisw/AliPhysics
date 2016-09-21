#ifndef ALI_ANALYSIS_TASK_PID_H
#define ALI_ANALYSIS_TASK_PID_H

/*
This task collects PID output from different detectors.
Only tracks fulfilling some standard quality cuts are taken into account.
At the moment, only data from TPC and TOF is collected. But in future,
data from e.g. HMPID is also foreseen.

Class written by Benjamin Hess.
Contact: bhess@cern.ch
*/

class TF1;
class TRandom3;
class AliAnalysisFilter;
class AliCFContainer;
class AliESDEvent;
class AliMCEvent;
class AliMCParticle;
class AliPID;
class AliPIDCombined;
class AliPIDResponse;
class AliPPVsMultUtils;
class AliTOFPIDResponse;
class AliVEvent;
class AliVTrack;

#include "TAxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TString.h"

#include "AliCentrality.h"
#include "AliCFContainer.h"

#include "AliPID.h"
#include "AliAnalysisTaskPIDV0base.h"

class AliAnalysisTaskPID : public AliAnalysisTaskPIDV0base {
 public:
  AliAnalysisTaskPID();
  AliAnalysisTaskPID(const char *name);
  virtual ~AliAnalysisTaskPID();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t*);
  
  enum ErrorCode { kNoErrors = 1, kWarning = 0, kError = -1};
  
  enum dataAxes { kDataMCID = 0, kDataSelectSpecies = 1, kDataPt = 2, kDataDeltaPrimeSpecies = 3, kDataCentrality = 4,
                  kDataJetPt = 5, kDataZ = 6, kDataXi = 7, kDataCharge = 8, kDataTOFpidInfo = 9, kDataDistance = 10, kDataJt = 11,
                  kDataNumAxes = 12 };

  enum genAxes { kGenMCID = 0, kGenSelectSpecies = 1, kGenPt = 2, kGenDeltaPrimeSpecies = 3, kGenCentrality = 4,
                 kGenJetPt = 5, kGenZ = 6, kGenXi = 7, kGenCharge = 8, kGenTOFpidInfo = 9, kGenDistance = 10, kGenJt = 11,
                 kGenNumAxes = 12 };
  
  enum genYieldAxes { kGenYieldMCID = 0, kGenYieldPt = 1, kGenYieldCentrality = 2, kGenYieldJetPt = 3, kGenYieldZ = 4, kGenYieldXi = 5,
                      kGenYieldCharge = 6, kGenYieldDistance = 7, kGenYieldJt = 8, kGenYieldNumAxes = 9 };
  
  enum ptResolutionAxes { kPtResJetPt = 0, kPtResGenPt = 1, kPtResRecPt = 2, kPtResCharge = 3, kPtResCentrality = 4, kPtResNumAxes = 5 };
  
  enum qaSharedClsAxes { kQASharedClsJetPt = 0, kQASharedClsPt = 1, kQASharedClsNumSharedCls = 2, kQASharedClsPadRow = 3,
                         kQASharedClsNumAxes = 4 };
  
  enum dEdxCheckAxes { kDeDxCheckPID = 0, kDeDxCheckP = 1, kDeDxCheckJetPt = 2, kDeDxCheckEtaAbs = 3 , kDeDxCheckDeDx = 4,
                       kDeDxCheckNumAxes = 5 };
  
  enum binZeroStudyAxes { kBinZeroStudyCentrality = 0, kBinZeroStudyGenPt = 1, kBinZeroStudyGenEta = 2, kBinZeroStudyNumAxes = 3 };
  
  enum efficiencyAxes { kEffMCID = 0, kEffTrackPt = 1, kEffTrackEta = 2, kEffTrackCharge = 3, kEffCentrality = 4, kEffJetPt = 5,
                        kEffZ = 6, kEffXi = 7, kEffDistance = 8, kEffJt = 9, kEffNumAxes = 10 };
  
  enum EffSteps { kStepGenWithGenCuts = 0, kStepRecWithGenCuts = 1, kStepRecWithGenCutsMeasuredObs = 2,
                  kStepRecWithRecCutsMeasuredObs = 3, kStepRecWithRecCutsMeasuredObsPrimaries = 4,
                  kStepRecWithRecCutsMeasuredObsStrangenessScaled = 5, kStepRecWithRecCutsPrimaries = 6, kNumSteps = 7};
  
  enum TOFpidInfo { kNoTOFinfo = -2, kNoTOFpid = -1, kTOFpion = 0, kTOFkaon = 1, kTOFproton = 2, kNumTOFspecies = 3,
                    kNumTOFpidInfoBins = 5 };
  
  enum EventCounterType { kTriggerSel = 0, kTriggerSelAndVtxCut = 1, kTriggerSelAndVtxCutAndZvtxCutNoPileUpRejection = 2,
                          kTriggerSelAndVtxCutAndZvtxCut = 3 };
  
  enum CutHistoType { kMCPtHardCut = 0 };
  
  enum EventGenerator { kPythia6Perugia0 = 0, kPythia6Perugia2011 = 0 };
  
  static Int_t PDGtoMCID(Int_t pdg);
  
  static void GetJetTrackObservables(Double_t trackPt, Double_t jetPt, Double_t& z, Double_t& xi, Bool_t storeXi = kTRUE);
  
  static Double_t GetMCStrangenessFactorCMS(Int_t motherPDG, Double_t motherGenPt);
  static Double_t GetMCStrangenessFactorCMS(AliMCEvent* mcEvent, AliMCParticle* daughter);
  
  static Bool_t IsSecondaryWithStrangeMotherMC(AliMCEvent* mcEvent, Int_t partLabel);
  
  static AliAnalysisTaskPID::EventGenerator GetEventGenerator() { return fgEventGenerator; };
  static void SetEventGenerator(AliAnalysisTaskPID::EventGenerator value) { fgEventGenerator = value; };
  
  virtual void ConfigureTaskForCurrentEvent(AliVEvent* event);
  
  Int_t GetIndexOfChargeAxisData() const
    { return fStoreAdditionalJetInformation ? kDataCharge : kDataCharge - 3; };
  Int_t GetIndexOfChargeAxisGen() const
    { return fStoreAdditionalJetInformation ? kGenCharge : kGenCharge - 3; };
  Int_t GetIndexOfChargeAxisGenYield() const
    { return fStoreAdditionalJetInformation ? kGenYieldCharge : kGenYieldCharge - 3; };
  
  Int_t GetIndexOfTOFpidInfoAxisData() const
    { return fStoreAdditionalJetInformation ? kDataTOFpidInfo : kDataTOFpidInfo - 3; };
  Int_t GetIndexOfTOFpidInfoAxisGen() const
    { return fStoreAdditionalJetInformation ? kGenTOFpidInfo : kGenTOFpidInfo - 3; };
  
  Bool_t FillXsec(Double_t xsection)
    { if (!fh1Xsec) return kFALSE; fh1Xsec->Fill("<#sigma>", xsection); return kTRUE; };
  Bool_t FillPythiaTrials(Double_t avgTrials)
    { if (!fh1Trials) return kFALSE; fh1Trials->Fill("#sum{ntrials}", avgTrials); return kTRUE; };
    
  Bool_t FillEfficiencyContainer(const Double_t* values, EffSteps step, Double_t weight = 1.0);
  
  Bool_t FillGeneratedYield(const Double_t* values, Double_t weight = 1.0);
  Bool_t FillPtResolution(Int_t mcID, const Double_t* values);
  Bool_t FillGenJets(Double_t centralityPercentile, Double_t jetPt, Double_t norm = -1.);
  Bool_t FillRecJets(Double_t centralityPercentile, Double_t jetPt, Double_t norm = -1.);
  
  Bool_t IncrementEventCounter(Double_t centralityPercentile, EventCounterType type);
  
  Bool_t FillCutHisto(Double_t value, CutHistoType type);
  
  void PostOutputData();
  
  void PrintSettings(Bool_t printSystematicsSettings = kFALSE) const;
  
  void PrintSystematicsSettings() const;
  
  Bool_t ProcessTrack(const AliVTrack* track, Int_t particlePDGcode, Double_t centralityPercentile, Double_t jetPt, Bool_t isMBSelected = kFALSE, 
                      Bool_t isMultSelected = kTRUE, Bool_t storeXi = kTRUE, Double_t radialDistanceToJet = -1, Double_t jT = -1);
  
  ErrorCode GenerateDetectorResponse(ErrorCode errCode, Double_t mean, Double_t sigma, Double_t* responses,
                                     Int_t nResponses,
                                     Bool_t usePureGaus = kFALSE);
  ErrorCode SetParamsForConvolutedGaus(Double_t gausMean, Double_t gausSigma);
  
  const TString GetCentralityEstimator() const { return fCentralityEstimator; };
  const TString GetPPCentralityEstimator() const {
    TString ppCentEstimator = fCentralityEstimator; ppCentEstimator = ppCentEstimator.ReplaceAll("ppMult", ""); return ppCentEstimator; }
  void SetCentralityEstimator(TString estimator) { fCentralityEstimator = estimator; };
  
  Double_t GetCentralityPercentile(AliVEvent* evt) const;
  
  inline Double_t GetConvolutedGaussTransitionPar(Int_t index) const;
  
  Bool_t SetConvolutedGaussLambdaParameter(Double_t lambda);

  Bool_t GetInputFromOtherTask() const { return fInputFromOtherTask; };
  void SetInputFromOtherTask(Bool_t flag) { fInputFromOtherTask = flag; };
  
  Bool_t GetDoPID() const { return fDoPID; };
  void SetDoPID(Bool_t flag) { fDoPID = flag; };
  
  Bool_t GetDoEfficiency() const { return fDoEfficiency; };
  void SetDoEfficiency(Bool_t flag) { fDoEfficiency = flag; };
  
  Bool_t GetDoPtResolution() const { return fDoPtResolution; };
  void SetDoPtResolution(Bool_t flag) { fDoPtResolution = flag; };
  
  Bool_t GetDoDeDxCheck() const { return fDoDeDxCheck; };
  void SetDoDeDxCheck(Bool_t flag) { fDoDeDxCheck = flag; };
  
  Bool_t GetDoBinZeroStudy() const { return fDoBinZeroStudy; };
  void SetDoBinZeroStudy(Bool_t flag) { fDoBinZeroStudy = flag; };
  
  Bool_t GetStoreCentralityPercentile() const { return fStoreCentralityPercentile; };
  void SetStoreCentralityPercentile(Bool_t flag) { fStoreCentralityPercentile = flag; };
  
  Bool_t GetStoreAdditionalJetInformation() const { return fStoreAdditionalJetInformation; };
  void SetStoreAdditionalJetInformation(Bool_t flag) { fStoreAdditionalJetInformation = flag; };
  
  Bool_t GetUseMCidForGeneration() const { return fUseMCidForGeneration; };
  void SetUseMCidForGeneration(Bool_t flag) { fUseMCidForGeneration = flag; };
  
  Bool_t GetUseConvolutedGaus() const { return fUseConvolutedGaus; };
  void SetUseConvolutedGaus(Bool_t flag) { fUseConvolutedGaus = flag; };
  
  Double_t GetAccuracyNonGaussianTail() const { return fAccuracyNonGaussianTail; };
  void SetAccuracyNonGaussianTail(Double_t value) { fAccuracyNonGaussianTail = value; };
  
  Bool_t GetTakeIntoAccountMuons() const { return fTakeIntoAccountMuons; };
  void SetTakeIntoAccountMuons(Bool_t flag) { fTakeIntoAccountMuons = flag; };
  
  Int_t GetTOFmode() const { return fTOFmode; };
  void SetTOFmode(Int_t tofMode) { fTOFmode = tofMode; };
  
  Bool_t GetUseTPCDefaultPriors() const { return fTPCDefaultPriors; };
  void SetUseTPCDefaultPriors(Bool_t flag) { fTPCDefaultPriors = flag; };
  
  Bool_t GetUsePriors() const { return fUsePriors; };
  void SetUsePriors(Bool_t flag) { fUsePriors = flag; };
  
  Bool_t GetUseITS() const { return fUseITS; };
  void SetUseITS(Bool_t flag) { fUseITS = flag; };
  
  Bool_t GetUseTOF() const { return fUseTOF; };
  void SetUseTOF(Bool_t flag) { fUseTOF = flag; };
  
  Bool_t GetStoreTOFInfo() const { return fStoreTOFInfo; }
  void SetStoreTOFInfo(Bool_t flag) { fStoreTOFInfo = flag; }
  
  Bool_t GetStoreCharge() const { return fStoreCharge; }
  void SetStoreCharge(Bool_t flag) { fStoreCharge = flag; }
  
  Double_t GetEtaAbsCutLow() const { return fEtaAbsCutLow; };
  Double_t GetEtaAbsCutUp() const { return fEtaAbsCutUp; };
  Bool_t SetEtaAbsCutRange(Double_t lowerLimit, Double_t upperLimit);
  
  Bool_t IsInAcceptedEtaRange(Double_t etaAbs) const { return (etaAbs >= fEtaAbsCutLow && etaAbs <= fEtaAbsCutUp); };
  
  Double_t GetSystematicScalingSplinesThreshold() const { return fSystematicScalingSplinesThreshold; };
  void SetSystematicScalingSplinesThreshold(Double_t threshold) { fSystematicScalingSplinesThreshold = threshold; };
  
  Double_t GetSystematicScalingSplinesBelowThreshold() const { return fSystematicScalingSplinesBelowThreshold; };
  void SetSystematicScalingSplinesBelowThreshold(Double_t scaleFactor) 
    { fSystematicScalingSplinesBelowThreshold = scaleFactor; CheckDoAnyStematicStudiesOnTheExpectedSignal(); };
    
  Double_t GetSystematicScalingSplinesAboveThreshold() const { return fSystematicScalingSplinesAboveThreshold; };
  void SetSystematicScalingSplinesAboveThreshold(Double_t scaleFactor) 
    { fSystematicScalingSplinesAboveThreshold = scaleFactor; CheckDoAnyStematicStudiesOnTheExpectedSignal(); };
  
  Double_t GetSystematicScalingEtaCorrectionMomentumThr() const { return fSystematicScalingEtaCorrectionMomentumThr; };
  void SetSystematicScalingEtaCorrectionMomentumThr(Double_t threshold) { fSystematicScalingEtaCorrectionMomentumThr = threshold; };
  
  Double_t GetSystematicScalingEtaCorrectionLowMomenta() const { return fSystematicScalingEtaCorrectionLowMomenta; };
  void SetSystematicScalingEtaCorrectionLowMomenta(Double_t scaleFactor) 
    { fSystematicScalingEtaCorrectionLowMomenta = scaleFactor; CheckDoAnyStematicStudiesOnTheExpectedSignal(); };
  
  Double_t GetSystematicScalingEtaCorrectionHighMomenta() const { return fSystematicScalingEtaCorrectionHighMomenta; };
  void SetSystematicScalingEtaCorrectionHighMomenta(Double_t scaleFactor) 
    { fSystematicScalingEtaCorrectionHighMomenta = scaleFactor; CheckDoAnyStematicStudiesOnTheExpectedSignal(); };
  
  Double_t GetSystematicScalingEtaSigmaParaThreshold() const { return fSystematicScalingEtaSigmaParaThreshold; };
  void SetSystematicScalingEtaSigmaParaThreshold(Double_t threshold) { fSystematicScalingEtaSigmaParaThreshold = threshold; };
  
  Double_t GetSystematicScalingEtaSigmaParaBelowThreshold() const { return fSystematicScalingEtaSigmaParaBelowThreshold; };
  void SetSystematicScalingEtaSigmaParaBelowThreshold(Double_t scaleFactor)
    { fSystematicScalingEtaSigmaParaBelowThreshold = scaleFactor; CheckDoAnyStematicStudiesOnTheExpectedSignal(); };
  
  Double_t GetSystematicScalingEtaSigmaParaAboveThreshold() const { return fSystematicScalingEtaSigmaParaAboveThreshold; };
  void SetSystematicScalingEtaSigmaParaAboveThreshold(Double_t scaleFactor)
    { fSystematicScalingEtaSigmaParaAboveThreshold = scaleFactor; CheckDoAnyStematicStudiesOnTheExpectedSignal(); };
  
  Double_t GetSystematicScalingMultCorrection() const { return fSystematicScalingMultCorrection; };
  void SetSystematicScalingMultCorrection(Double_t scaleFactor) 
    { fSystematicScalingMultCorrection = scaleFactor; CheckDoAnyStematicStudiesOnTheExpectedSignal(); };
  
  Double_t GetMaxEtaVariation(Double_t dEdxSplines);
  Bool_t CalculateMaxEtaVariationMapFromPIDResponse();
  
  void CleanupParticleFractionHistos();
  Bool_t GetParticleFraction(Double_t trackPt, Double_t jetPt, Double_t multiplicity,
                             AliPID::EParticleType species, Double_t& fraction, Double_t& fractionErrorStat,
                             Double_t& fractionErrorSys) const;
  Bool_t GetParticleFractions(Double_t trackPt, Double_t jetPt, Double_t centralityPercentile,
                              Double_t* prob, Int_t smearSpeciesByError, Int_t takeIntoAccountSpeciesSysError,
                              Bool_t uniformSystematicError = kFALSE) const;
  const TH3D* GetParticleFractionHisto(Int_t species, Bool_t sysError = kFALSE) const;
  Bool_t SetParticleFractionHisto(const TH3D* hist, Int_t species, Bool_t sysError = kFALSE);
  Int_t GetParticleFractionHistoNbinsTrackPt() const;
  Int_t GetParticleFractionHistoNbinsJetPt() const;
  Int_t GetParticleFractionHistoNbinsCentrality() const;
  Bool_t SetParticleFractionHistosFromFile(const TString filePathName, Bool_t sysError = kFALSE);
  Int_t GetRandomParticleTypeAccordingToParticleFractions(Double_t trackPt, Double_t jetPt, 
                                                          Double_t centralityPercentile,
                                                          Bool_t smearByError,
                                                          Bool_t takeIntoAccountSysError = kFALSE) const;
  
  TOFpidInfo GetTOFType(const AliVTrack* track, Int_t tofMode) const;
  
  //Underlying event
  Bool_t GetIsUEPID() const { return fIsUEPID; }
  void SetIsUEPID(Bool_t flag) { fIsUEPID = flag; }
  void FillUEDensity(Double_t cent, Double_t UEpt);
  void FillJetArea(Double_t cent, Double_t area);
  void NormalizeJetArea(Double_t jetParameter);
  
 protected:
  void CheckDoAnyStematicStudiesOnTheExpectedSignal();
  Double_t ConvolutedGaus(const Double_t* xx, const Double_t* par) const;
  inline Double_t FastGaus(Double_t x, Double_t mean, Double_t sigma) const;
  inline Double_t FastNormalisedGaus(Double_t x, Double_t mean, Double_t sigma) const;
  Int_t FindBinWithinRange(TAxis* axis, Double_t value) const;
  Int_t FindFirstBinAboveIn3dSubset(const TH3* hist, Double_t threshold, Int_t yValue, Int_t zValue) const;
  Int_t FindLastBinAboveIn3dSubset(const TH3* hist, Double_t threshold, Int_t yValue, Int_t zValue) const;
  virtual void SetUpGenHist(THnSparse* hist, Double_t* binsPt, Double_t* binsDeltaPrime, Double_t* binsCent, Double_t* binsJetPt,
                            Double_t* binsJt) const;
  virtual void SetUpGenYieldHist(THnSparse* hist, Double_t* binsPt, Double_t* binsCent, Double_t* binsJetPt, Double_t* binsJt) const;
  virtual void SetUpHist(THnSparse* hist, Double_t* binsPt, Double_t* binsDeltaPrime, Double_t* binsCent, Double_t* binsJetPt,
                         Double_t* binsJt) const;
  virtual void SetUpPtResHist(THnSparse* hist, Double_t* binsPt, Double_t* binsJetPt, Double_t* binsCent) const;
  virtual void SetUpSharedClsHist(THnSparse* hist, Double_t* binsPt, Double_t* binsJetPt) const;
  virtual void SetUpDeDxCheckHist(THnSparse* hist, const Double_t* binsPt, const Double_t* binsJetPt, const Double_t* binsEtaAbs) const;
  virtual void SetUpBinZeroStudyHist(THnSparse* hist, const Double_t* binsCent, const Double_t* binsPt) const;
  virtual void SetUpPIDcombined();
  
  static const Int_t fgkNumJetAxes; // Number of additional axes for jets
  static const Double_t fgkEpsilon; // Double_t threshold above zero
  static const Int_t fgkMaxNumGenEntries; // Maximum number of generated detector responses per track and delta(Prime) and associated species

  static const Double_t fgkOneOverSqrt2; // = 1. / TMath::Sqrt2();
  
 private:
  Int_t fRun; // Current run number
  AliPIDCombined* fPIDcombined; //! PID combined object
  AliPPVsMultUtils* fPPVsMultUtils; //! Utilities for pp vs. mult analysis
  
  Bool_t fInputFromOtherTask; // If set to kTRUE, no events are processed and the input must be fed in from another task. If set to kFALSE, normal event processing
  
  Bool_t fDoPID; // Do PID processing (and post the output), if flag is set to kTRUE
  Bool_t fDoEfficiency; // Do efficiency processing (and post the output), if flag is set to kTRUE
  Bool_t fDoPtResolution; // Do pT resolution processing (and post the output), if flag is set to kTRUE
  Bool_t fDoDeDxCheck; // Check dEdx, if flag set to kTRUE
  Bool_t fDoBinZeroStudy; // Do bin zero study, if flag is set to kTRUE
  
  static AliAnalysisTaskPID::EventGenerator fgEventGenerator;
  
  Bool_t fStoreCentralityPercentile; // If set to kTRUE, store centrality percentile for each event. In case of kFALSE (appropriate for pp), centrality percentile will be set to -1 for every event
  Bool_t fStoreAdditionalJetInformation; // If set to kTRUE, additional jet information like jetPt, z, xi will be stored in the THnSparses

  Bool_t fTakeIntoAccountMuons; // Also take into account muons for the generation of the expected response and the most probable PID
  Bool_t fUseITS; // Use ITS for PID combined probabilities
  Bool_t fUseTOF; // Use TOF for PID combined probabilities
  Bool_t fStoreTOFInfo;  //Store the TOF Info in the histogram
  Bool_t fUsePriors; // Use priors for PID combined probabilities
  Bool_t fTPCDefaultPriors; // Use TPC default priors for PID combined probabilities, if priors are enabled
  
  Bool_t fStoreCharge; // Store positive and negative charged tracks separately
    
  Bool_t fUseMCidForGeneration; // If MC, use MCid instead of PIDcombined to generate the signals
  
  Bool_t fUseConvolutedGaus; // Use convoluted gaus to generate detector response instead of pure gaus  
  const Int_t fkConvolutedGausNPar; // Number of parameters for convoluted gaus
  Double_t fAccuracyNonGaussianTail; // Accuracy of the non-gaussian tail (fraction of the maximum)
  const Double_t fkDeltaPrimeLowLimit; // Lower deltaPrime limit
  const Double_t fkDeltaPrimeUpLimit; // Upper deltaPrime limit
  TF1* fConvolutedGausDeltaPrime; // Gaus convoluted with exponential tail to generate detector response (deltaPrime)
  
  Int_t fTOFmode; // TOF mode used for TOF PID info (affects num sigma inclusion/exclusion)
  Double_t fConvolutedGaussTransitionPars[3]; // Parameter for transition from gaussian parameters to asymmetric shape
  static const Double_t fgkSigmaReferenceForTransitionPars; // Reference sigma chosen to calculate transition parameters
  
  Double_t fEtaAbsCutLow; // Lower cut value on |eta|
  Double_t fEtaAbsCutUp;  // Upper cut value on |eta|
  
  // For systematic studies
  Bool_t   fDoAnySystematicStudiesOnTheExpectedSignal; // Internal flag indicating whether any systematic studies are going to be performed
  Double_t fSystematicScalingSplinesThreshold;         // beta-gamma threshold for the systematic spline scale factor
  Double_t fSystematicScalingSplinesBelowThreshold;        // Systematic scale factor for the splines (1. = no systematics) below threshold
  Double_t fSystematicScalingSplinesAboveThreshold;        // Systematic scale factor for the splines (1. = no systematics) above threshold
  Double_t fSystematicScalingEtaCorrectionMomentumThr;  // Momentum threshold for the systematic scale factor for the eta correction (separates low-p from high-p
  Double_t fSystematicScalingEtaCorrectionLowMomenta;   // Systematic scale factor for the eta correction (1. = no systematics) at low momenta
  Double_t fSystematicScalingEtaCorrectionHighMomenta;  // Systematic scale factor for the eta correction (1. = no systematics) at high momenta
  
  Double_t fSystematicScalingEtaSigmaParaThreshold; // dEdx threshold for the systematic scale factor for the parametrisation of the relative signal width
  Double_t fSystematicScalingEtaSigmaParaBelowThreshold; // Systematic scale factor for the parametrisation of the relative signal width (1. = no systematics) below threshold
  Double_t fSystematicScalingEtaSigmaParaAboveThreshold; // Systematic scale factor for the parametrisation of the relative signal width (1. = no systematics) above threshold 
  Double_t fSystematicScalingMultCorrection; // Systematic scale factor for the multiplicity correction (1. = no systematics) 
  
  
  
  TH3D* fFractionHists[AliPID::kSPECIES]; // 3D histos of particle fraction as function  with trackPt, jetPt (-1 for inclusive spectra), centralityPercentile (-1 for pp)
  TH3D* fFractionSysErrorHists[AliPID::kSPECIES]; // 3D histos of sys. error of particle fraction as function  with trackPt, jetPt (-1 for inclusive spectra), centralityPercentile (-1 for pp)
  
  TString fCentralityEstimator; // Estimator for the centrality (e.g. V0A, V0M)
  
  THnSparseD* fhPIDdataAll; //! Data histo
  
  // Generated response information
  THnSparseD* fhGenEl; //! Generated response for el
  THnSparseD* fhGenKa; //! Generated response for ka
  THnSparseD* fhGenPi; //! Generated response for pi
  THnSparseD* fhGenMu; //! Generated response for mu
  THnSparseD* fhGenPr; //! Generated response for pr
  
  // Generated responses for a single track
  Double_t* fGenRespElDeltaPrimeEl; //! Generated responses for a single track
  Double_t* fGenRespElDeltaPrimeKa; //! Generated responses for a single track
  Double_t* fGenRespElDeltaPrimePi; //! Generated responses for a single track
  Double_t* fGenRespElDeltaPrimePr; //! Generated responses for a single track
  Double_t* fGenRespKaDeltaPrimeEl; //! Generated responses for a single track
  Double_t* fGenRespKaDeltaPrimeKa; //! Generated responses for a single track
  Double_t* fGenRespKaDeltaPrimePi; //! Generated responses for a single track
  Double_t* fGenRespKaDeltaPrimePr; //! Generated responses for a single track
  Double_t* fGenRespPiDeltaPrimeEl; //! Generated responses for a single track
  Double_t* fGenRespPiDeltaPrimeKa; //! Generated responses for a single track
  Double_t* fGenRespPiDeltaPrimePi; //! Generated responses for a single track
  Double_t* fGenRespPiDeltaPrimePr; //! Generated responses for a single track
  Double_t* fGenRespMuDeltaPrimeEl; //! Generated responses for a single track
  Double_t* fGenRespMuDeltaPrimeKa; //! Generated responses for a single track
  Double_t* fGenRespMuDeltaPrimePi; //! Generated responses for a single track
  Double_t* fGenRespMuDeltaPrimePr; //! Generated responses for a single track
  Double_t* fGenRespPrDeltaPrimeEl; //! Generated responses for a single track
  Double_t* fGenRespPrDeltaPrimeKa; //! Generated responses for a single track
  Double_t* fGenRespPrDeltaPrimePi; //! Generated responses for a single track
  Double_t* fGenRespPrDeltaPrimePr; //! Generated responses for a single track
  /*
  Double_t* fGenRespElDeltaEl; //! Generated responses for a single track
  Double_t* fGenRespElDeltaKa; //! Generated responses for a single track
  Double_t* fGenRespElDeltaPi; //! Generated responses for a single track
  Double_t* fGenRespElDeltaPr; //! Generated responses for a single track
  Double_t* fGenRespKaDeltaEl; //! Generated responses for a single track
  Double_t* fGenRespKaDeltaKa; //! Generated responses for a single track
  Double_t* fGenRespKaDeltaPi; //! Generated responses for a single track
  Double_t* fGenRespKaDeltaPr; //! Generated responses for a single track
  Double_t* fGenRespPiDeltaEl; //! Generated responses for a single track
  Double_t* fGenRespPiDeltaKa; //! Generated responses for a single track
  Double_t* fGenRespPiDeltaPi; //! Generated responses for a single track
  Double_t* fGenRespPiDeltaPr; //! Generated responses for a single track
  Double_t* fGenRespMuDeltaEl; //! Generated responses for a single track
  Double_t* fGenRespMuDeltaKa; //! Generated responses for a single track
  Double_t* fGenRespMuDeltaPi; //! Generated responses for a single track
  Double_t* fGenRespMuDeltaPr; //! Generated responses for a single track
  Double_t* fGenRespPrDeltaEl; //! Generated responses for a single track
  Double_t* fGenRespPrDeltaKa; //! Generated responses for a single track
  Double_t* fGenRespPrDeltaPi; //! Generated responses for a single track
  Double_t* fGenRespPrDeltaPr; //! Generated responses for a single track
  */
  
  TAxis* fDeltaPrimeAxis; //! Axis holding the deltaPrime binning
  TH1D* fhMaxEtaVariation; //! Histo holding the maximum deviation of the eta correction factor from unity vs. 1/dEdx(splines)
  
  TH1D* fhEventsProcessed; //! Histo holding the number of processed events (i.e. passing trigger selection, vtx and zvtx cuts and (if enabled) pile-up rejection)
  TH1D* fhEventsTriggerSel; //! Histo holding the number of events passing trigger selection
  TH1D* fhEventsTriggerSelVtxCut; //! Histo holding the number of events passing trigger selection and vtx cut
  TH1D* fhEventsProcessedNoPileUpRejection; //! Histo holding the number of processed events before pile-up rejection
  
  THnSparseD* fChargedGenPrimariesTriggerSel; //! Histo holding the generated charged primary yields for triggered events
  THnSparseD* fChargedGenPrimariesTriggerSelVtxCut; //! Histo holding the generated charged primary yields for triggered events passing vertex cuts
  THnSparseD* fChargedGenPrimariesTriggerSelVtxCutZ; //! Histo holding the generated charged primary yields for triggered events passing vertex cuts (including cut on z)
  THnSparseD* fChargedGenPrimariesTriggerSelVtxCutZPileUpRej; //! Histo holding the generated charged primary yields for triggered events passing vertex cuts (including cut on z) and pile-up rejection
  
  THnSparseD* fhMCgeneratedYieldsPrimaries; //! Histo holding the generated (no reco, no cuts) primary particle yields in considered eta range
  
  TH2D* fh2FFJetPtRec;            //! Number of reconstructed jets vs. jetPt and centrality
  TH2D* fh2FFJetPtGen;            //! Number of generated jets vs. jetPt and centrality
  
  TProfile* fh1Xsec;              //! pythia cross section and trials
  TH1D*     fh1Trials;            //! sum of trials
  TH1F*     fh1EvtsPtHardCut;     //! Number events before and after the cut on MC pT hard
  
  AliCFContainer* fContainerEff; //! Container for efficiency determination
  
  THnSparseD* fPtResolution[AliPID::kSPECIES]; //! Pt Resolution for the individual species
  THnSparseD* fQASharedCls; //! QA for shared clusters
  
  THnSparseD* fDeDxCheck; //! dEdx check
  
  TObjArray* fOutputContainer;  //! output data container
  
  TObjArray* fQAContainer; //! output data container for QA
  
  AliAnalysisTaskPID(const AliAnalysisTaskPID&); // not implemented
  AliAnalysisTaskPID& operator=(const AliAnalysisTaskPID&); // not implemented
  
  //For Underlying Event
  Bool_t fIsUEPID;
  TH2D* fh2UEDensity;
  TH1D* fh1JetArea;
  
  ClassDef(AliAnalysisTaskPID, 23);
};


//_____________________________________________________________________________
inline Bool_t AliAnalysisTaskPID::FillEfficiencyContainer(const Double_t* values, AliAnalysisTaskPID::EffSteps step,
                                                          Double_t weight) 
{
  // Fill efficiency container at step "step" with the values
  
  if (!fDoEfficiency)
    return kFALSE;
  
  if (!fContainerEff) {
    AliError("Efficiency container not initialised -> cannot be filled!");
    return kFALSE;
  }
  
  fContainerEff->Fill(values, step, weight);    
  
  return kTRUE;
}


//_____________________________________________________________________________
inline Bool_t AliAnalysisTaskPID::FillGeneratedYield(const Double_t* values, Double_t weight)
{
  // Fill histos with generated primary yields with provided values
  
  if (!fDoPID)
    return kFALSE;
  
  if (!fhMCgeneratedYieldsPrimaries) {
    AliError("Histo for generated primary yield not initialised -> cannot be filled!");
    return kFALSE;
  }
  
  fhMCgeneratedYieldsPrimaries->Fill(values, weight);
    
  return kTRUE;
}


//_____________________________________________________________________________
inline Bool_t AliAnalysisTaskPID::FillGenJets(Double_t centralityPercentile, Double_t jetPt, Double_t norm)
{
  if (!fDoPID && !fDoEfficiency)
    return kFALSE;
  
  if (!fh2FFJetPtGen)
    return kFALSE;
  
  if (norm > 0.)
    fh2FFJetPtGen->Fill(centralityPercentile, jetPt, 1. / norm);
  else
    fh2FFJetPtGen->Fill(centralityPercentile, jetPt);
  
  return kTRUE;
}


//_____________________________________________________________________________
inline Bool_t AliAnalysisTaskPID::FillRecJets(Double_t centralityPercentile, Double_t jetPt, Double_t norm)
{
  if (!fDoPID && !fDoEfficiency)
    return kFALSE;
  
  if (!fh2FFJetPtRec)
    return kFALSE;
  
  if (norm > 0.)
    fh2FFJetPtRec->Fill(centralityPercentile, jetPt, 1. / norm);
  else
    fh2FFJetPtRec->Fill(centralityPercentile, jetPt);
  
  return kTRUE;
}


//_____________________________________________________________________________
inline Bool_t AliAnalysisTaskPID::FillPtResolution(Int_t mcID, const Double_t* values)
{
  // Fill histos with pT resolution with provided values
  
  if (!fDoPtResolution || mcID < 0 || mcID >= AliPID::kSPECIES)
    return kFALSE;
  
  if (!fPtResolution[mcID]) {
    AliError(Form("Histo for pT resolution (species: %s) not initialised -> cannot be filled!", AliPID::ParticleName(mcID)));
    return kFALSE;
  }
  
  fPtResolution[mcID]->Fill(values);
    
  return kTRUE;
}
 

//_____________________________________________________________________________
inline Bool_t AliAnalysisTaskPID::IncrementEventCounter(Double_t centralityPercentile, AliAnalysisTaskPID::EventCounterType type)
{
  // Increment the number of events for the given centrality percentile and the corresponding counter type
  
  if (type == kTriggerSel) {
    if (!fhEventsTriggerSel) {
      AliError("Histogram for number of events (kTriggerSel) not initialised -> cannot be incremented!");
      return kFALSE;
    }
    
    fhEventsTriggerSel->Fill(centralityPercentile);
  }
  else if (type == kTriggerSelAndVtxCut) {
    if (!fhEventsTriggerSelVtxCut) {
      AliError("Histogram for number of events (kTriggerSelAndVtxCut) not initialised -> cannot be incremented!");
      return kFALSE;
    }
    
    fhEventsTriggerSelVtxCut->Fill(centralityPercentile);
  }
  else if (type == kTriggerSelAndVtxCutAndZvtxCut) {
    if (!fhEventsProcessed) {
      AliError("Histogram for number of events (kTriggerSelAndVtxCutAndZvtxCut) not initialised -> cannot be incremented!");
      return kFALSE;
    }
    
    fhEventsProcessed->Fill(centralityPercentile);
  }
  else if (type == kTriggerSelAndVtxCutAndZvtxCutNoPileUpRejection) {
    if (!fhEventsProcessedNoPileUpRejection) {
      AliError("Histogram for number of events (kTriggerSelAndVtxCutAndZvtxCutNoPileUpRejection) not initialised -> cannot be incremented!");
      return kFALSE;
    }
    
    fhEventsProcessedNoPileUpRejection->Fill(centralityPercentile);
  }
  
  
  return kTRUE;
};


//_____________________________________________________________________________
inline Bool_t AliAnalysisTaskPID::FillCutHisto(Double_t value, AliAnalysisTaskPID::CutHistoType type)
{
  // Fill cut histo of corresponding type at given value.
  if (type == kMCPtHardCut) {
    if (!fh1EvtsPtHardCut) {
      AliError("Histogram \"fh1EvtsPtHardCut\" not initialised -> cannot be filled!");
      return kFALSE;
    }
    
    fh1EvtsPtHardCut->Fill(value);
  }
  
  return kTRUE;
}


//_____________________________________________________________________________
inline Bool_t AliAnalysisTaskPID::SetEtaAbsCutRange(Double_t lowerLimit, Double_t upperLimit)
{
  if (lowerLimit >= upperLimit) {
    AliError(Form("Requested lower |eta| cut limit >= upper |eta| cut limit. Old eta cut range will be used (low %f, high %f).",
                  fEtaAbsCutLow, fEtaAbsCutUp));
    return kFALSE;
  }
  
  fEtaAbsCutLow = lowerLimit;
  fEtaAbsCutUp = upperLimit;
  
  return kTRUE;
};


//_____________________________________________________________________________
inline Double_t AliAnalysisTaskPID::GetConvolutedGaussTransitionPar(Int_t index) const
{
  if (index < 0 || index >= 3) {
    printf("Invalid index %d!\n", index);
    return -1;
  }
  return fConvolutedGaussTransitionPars[index];
}


//_____________________________________________________________________________
inline Int_t AliAnalysisTaskPID::GetParticleFractionHistoNbinsTrackPt() const
{
  if (!fFractionHists[AliPID::kPion])
    return -1;
  
  return fFractionHists[AliPID::kPion]->GetNbinsX();
}


//_____________________________________________________________________________
inline Int_t AliAnalysisTaskPID::GetParticleFractionHistoNbinsJetPt() const
{
  if (!fFractionHists[AliPID::kPion])
    return -1;
  
  return fFractionHists[AliPID::kPion]->GetNbinsY();
}


//_____________________________________________________________________________
inline Int_t AliAnalysisTaskPID::GetParticleFractionHistoNbinsCentrality() const
{
  if (!fFractionHists[AliPID::kPion])
    return -1;
  
  return fFractionHists[AliPID::kPion]->GetNbinsZ();
}


//_____________________________________________________________________________
inline Double_t AliAnalysisTaskPID::GetCentralityPercentile(AliVEvent* evt) const
{
  // WARNING: This function may not be used in case of special pp centrality estimators which require different handling
  // (and sometimes ESD events)
  if (!evt)
    return -1;
  
  Double_t centralityPercentile = -1.;
  if (fStoreCentralityPercentile)
    centralityPercentile = evt->GetCentrality()->GetCentralityPercentile(fCentralityEstimator.Data());
  
  return centralityPercentile;
}


//_____________________________________________________________________________
inline void AliAnalysisTaskPID::PostOutputData()
{
  PostData(1, fOutputContainer);
  
  if (fDoEfficiency)
    PostData(2, fContainerEff);
  
  if (fDoPtResolution || fDoDeDxCheck)
    PostData(3, fQAContainer);
}

#endif
