#ifndef ALITPCRECOPARAM_H
#define ALITPCRECOPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \class AliTPCRecoParam
/// \brief Class with TPC reconstruction parameters


#include "AliDetectorRecoParam.h"
#include "TVectorF.h"

class AliTPCRecoParam : public AliDetectorRecoParam
{
 public:
  enum {                       // methods used for correction maps time dependence
    kCorrMapInterpolation         // interpolate between 2 nearest timebins maps
    ,kCorrMapNoScaling            // no scaling, use just the single map matching to timestamp
    ,kCorrMapGlobalScalingLumi // scale current map by ratio of inst_lumi/<lumi_timebin>
  };
 public:
  AliTPCRecoParam();
  AliTPCRecoParam(const AliTPCRecoParam& src);
  AliTPCRecoParam& operator=(const AliTPCRecoParam& src);
  virtual ~AliTPCRecoParam();
  virtual void Print(const Option_t* option="") const;
  static   Bool_t  GetUseTimeCalibration();
  static   void    SetUseTimeCalibration(Bool_t useTimeCalibration);

  void     SetUseHLTClusters(Int_t useHLTClusters){fUseHLTClusters=useHLTClusters;}
  Int_t    GetUseHLTClusters() const {return fUseHLTClusters;}
  void     SetUseHLTPreSeeding(Int_t useHLTPreSeeding){fUseHLTPreSeeding=useHLTPreSeeding;}
  Int_t    GetUseHLTPreSeeding() const {return fUseHLTPreSeeding;}
  void     SetClusterSharing(Bool_t sharing){fBClusterSharing=sharing;}
  Bool_t   GetClusterSharing() const {return fBClusterSharing;}
  Double_t GetCtgRange() const     { return fCtgRange;}
  Double_t GetMaxSnpTracker() const{ return fMaxSnpTracker;}
  Double_t GetMaxSnpTrack() const  { return fMaxSnpTrack;}
  Bool_t   GetUseOuterDetectors() const { return fUseOuterDetectors;}
  void     SetUseOuterDetectors(Bool_t flag)  { fUseOuterDetectors=flag;}
  void     SetMaxChi2TPCTRD(Double_t maxChi2){fMaxChi2TPCTRD=maxChi2;}
  Double_t GetMaxChi2TPCTRD() const {return fMaxChi2TPCTRD;}
  void     SetMaxChi2TPCITS(Double_t maxChi2){fMaxChi2TPCITS=maxChi2;}
  Double_t GetMaxChi2TPCITS() const {return fMaxChi2TPCITS;}
  Double_t GetCutSharedClusters(Int_t index)const { return fCutSharedClusters[index];}
  void  SetCutSharedClusters(Int_t index, Float_t value){ fCutSharedClusters[index]=value;}
  Int_t GetClusterMaxRange(Int_t index)const { return fClusterMaxRange[index];}
  void     SetClusterMaxRange(Int_t index, Int_t value){ fClusterMaxRange[index]=value;}
  //
  Int_t    GetAccountDistortions()               const {return fAccountDistortions;}
  void     SetAccountDistortions(Int_t v)              {fAccountDistortions = v;}
  Bool_t   GetUseCorrectionMap()                 const {return fUseCorrectionMap;}
  void     SetUseCorrectionMap(Bool_t v=kTRUE)         {fUseCorrectionMap = v;}
  //
  // Outlier filtering configuration
  //
  Int_t   GetUseOulierClusterFilter() const { return fUseOulierClusterFilter;}  // swith to use outlier cluster filter
  void    SetUseOulierClusterFilter(Int_t value){ fUseOulierClusterFilter=value;}  // swith to use outlier cluster filter
  //
  Bool_t   DumpSignal()     const  { return fDumpSignal;}
  void     SetTimeInterval(Int_t first, Int_t last) { fFirstBin=first, fLastBin =last;}
  Int_t    GetFirstBin() const     { return fFirstBin;}
  Int_t    GetLastBin() const      { return fLastBin;}
  void     SetTimeBinRange(Int_t first, Int_t last){ fFirstBin = first; fLastBin = last;}
  Bool_t   GetCalcPedestal()       const  { return fBCalcPedestal;}
  Bool_t   GetDoUnfold()           const  { return fBDoUnfold;}
  void     SetDoUnfold(Bool_t unfold)     { fBDoUnfold = unfold;}
  Float_t  GetDumpAmplitudeMin()   const  { return fDumpAmplitudeMin;}
  Float_t  GetMaxNoise()           const  { return fMaxNoise;}
  //
  Int_t    GetUseOnePadCluster()   const  { return fUseOnePadCluster;}
  Bool_t   GetUseHLTOnePadCluster()const  { return fUseHLTOnePadCluster;}
  Float_t  GetMinMaxCutAbs()       const  { return fMinMaxCutAbs; }
  Float_t  GetMinLeftRightCutAbs() const  { return fMinLeftRightCutAbs;}  // minimal amplitude left right - PRF
  Float_t  GetMinUpDownCutAbs()    const  { return fMinUpDownCutAbs;}  // minimal amplitude up-down - TRF
  Float_t  GetMinMaxCutSigma()       const  { return fMinMaxCutSigma; }
  Float_t  GetMinLeftRightCutSigma() const  { return fMinLeftRightCutSigma;}  // minimal amplitude left right - PRF
  Float_t  GetMinUpDownCutSigma()    const  { return fMinUpDownCutSigma;}  // minimal amplitude up-down - TRF
  //
  void SetUseOnePadCluster(Int_t use)      {   fUseOnePadCluster = use;}
  void SetUseHLTOnePadCluster(Bool_t use)  {   fUseHLTOnePadCluster = use;}
  void SetMinMaxCutAbs(Float_t th)         {   fMinMaxCutAbs=th; }
  void SetMinLeftRightCutAbs(Float_t th)   {   fMinLeftRightCutAbs=th;}  // minimal amplitude left right - PRF
  void SetMinUpDownCutAbs(Float_t th)      {   fMinUpDownCutAbs=th;}  // minimal amplitude up-down - TRF
  void SetMinMaxCutSigma(Float_t th)       {   fMinMaxCutSigma=th; }
  void SetMinLeftRightCutSigma(Float_t th) {   fMinLeftRightCutSigma=th;}  // minimal amplitude left right - PRF
  void SetMinUpDownCutSigma(Float_t th)    {   fMinUpDownCutSigma=th;}  // minimal amplitude up-down - TRF
  void  SetUseTotCharge(Bool_t flag) {fUseTotCharge = flag;}
  void  SetCtgRange(Double_t ctgRange) {fCtgRange = ctgRange;}
  void  SetUseMultiplicityCorrectionDedx(Bool_t flag) {fUseMultiplicityCorrectionDedx = flag;}

  void  SetUseAlignmentTime(Bool_t flag) {fUseAlignmentTime = flag;}
  void  SetNeighborRowsDedx(Int_t nRows) {fNeighborRowsDedx = nRows;}
  void SetCorrectionHVandPTMode(Int_t value){ fGainCorrectionHVandPTMode =value;}
  void SetSkipTimeBins(Double_t value) {fSkipTimeBins=value;}
  //
  Int_t    GetLastSeedRowSec()       const  { return fLastSeedRowSec;}
  Int_t    GetSeedGapPrim()        const  { return fSeedGapPrim;}
  Int_t    GetSeedGapSec()         const  { return fSeedGapSec;}
  void     SetDoKinks(Bool_t on)   { fBKinkFinder=on; }
  Bool_t   GetDoKinks() const      { return fBKinkFinder;}
  Double_t GetKinkAngleCutChi2(Int_t index) const {return fKinkAngleCutChi2[index];}
  void     SetKinkAngleCutChi2(Int_t index,Double_t value) {fKinkAngleCutChi2[index]=value;}
  void     SetSeedGapPrim(Int_t seedGapPrim)         { fSeedGapPrim = seedGapPrim;}
  void     SetSeedGapSec(Int_t seedGapSec)          { fSeedGapSec  = seedGapSec;}
  Float_t  GetMaxC()    const      { return fMaxC;}
  Bool_t   GetSpecialSeeding() const { return fBSpecialSeeding;}
  //
  //

  //
  // Correction setup
  //
  void  SetUseFieldCorrection(Int_t flag){fUseFieldCorrection=flag;}
  void  SetUseComposedCorrection(Bool_t flag){fUseComposedCorrection=flag;}
  void  SetUseRPHICorrection(Int_t flag){fUseRPHICorrection=flag;}
  void  SetUseRadialCorrection(Int_t flag){fUseRadialCorrection=flag;}
  void  SetUseQuadrantAlignment(Int_t flag){fUseQuadrantAlignment=flag;}
  void  SetUseSectorAlignment(Int_t flag){fUseSectorAlignment=flag;}
  void  SetUseDriftCorrectionTime(Int_t flag){fUseDriftCorrectionTime=flag;}
  void  SetUseDriftCorrectionGY(Int_t flag){fUseDriftCorrectionGY=flag;}
  void  SetUseGainCorrectionTime(Int_t flag){fUseGainCorrectionTime=flag;}
  void  SetUseExBCorrection(Int_t flag){fUseExBCorrection=flag;}
  void  SetUseTOFCorrection(Bool_t flag) {fUseTOFCorrection = flag;}
  void  SetUseIonTailCorrection(Int_t flag) {fUseIonTailCorrection = flag;}
  void  SetCrosstalkCorrection(Float_t crosstalkCorrection) {fCrosstalkCorrection= crosstalkCorrection; }
  void  SetCrosstalkCorrectionMissingCharge(Float_t crosstalkCorrection) {fCrosstalkCorrectionMissingCharge= crosstalkCorrection; }
  //
  Int_t  GetCorrMapTimeDepMethod()      const {return fCorrMapTimeDepMethod;}
  void   SetCorrMapTimeDepMethod(int m)       {fCorrMapTimeDepMethod = m;}
  Int_t  GetUseLumiType()               const {return fUseLumiType;}
  void   SetUseLumiType(int tp)               {fUseLumiType  =tp;}
  //
  Int_t GetUseFieldCorrection() const {return fUseFieldCorrection;}
  Int_t GetUseComposedCorrection() const {return fUseComposedCorrection;}
  Int_t GetUseRPHICorrection() const {return fUseRPHICorrection;}
  Int_t GetUseRadialCorrection() const {return fUseRadialCorrection;}
  Int_t GetUseQuadrantAlignment() const {return fUseQuadrantAlignment;}
  Int_t GetUseSectorAlignment() const {return fUseSectorAlignment;}
  Int_t GetUseDriftCorrectionTime() const {return fUseDriftCorrectionTime;}
  Int_t GetUseDriftCorrectionGY() const {return fUseDriftCorrectionGY;}
  Int_t GetUseGainCorrectionTime() const {return fUseGainCorrectionTime;}
  Int_t GetUseExBCorrection() const {return fUseExBCorrection;}
  Bool_t GetUseTOFCorrection() {return fUseTOFCorrection;}
  Int_t GetUseIonTailCorrection() const {return fUseIonTailCorrection;}
  Double_t GetCrosstalkCorrection() const {return fCrosstalkCorrection;}
 Double_t GetCrosstalkCorrectionMissingCharge() const {return fCrosstalkCorrectionMissingCharge;}

  Bool_t GetUseMultiplicityCorrectionDedx() const {return fUseMultiplicityCorrectionDedx;}
  Int_t  GetGainCorrectionHVandPTMode() const  { return   fGainCorrectionHVandPTMode;}
  Double_t  GetSkipTimeBins() const {return fSkipTimeBins;}

  Bool_t GetUseAlignmentTime() const {return fUseAlignmentTime;}
  //
  Bool_t   GetUseTotCharge() const {return fUseTotCharge;}          // switch use total or max charge
  Float_t  GetMinFraction() const {return fMinFraction;}           // truncated mean - lower threshold
  Float_t  GetMaxFraction() const {return fMaxFaction;}            // truncated mean - upper threshold
  Int_t    GetNeighborRowsDedx() const {return fNeighborRowsDedx;}

  //
  void     SetSystematicError(Double_t *systematic){ for (Int_t i=0; i<5;i++) fSystematicErrors[i]=systematic[i];}
  void     SetSystematicErrorCluster(Double_t *systematic){ for (Int_t i=0; i<2;i++) fSystematicErrorCluster[i]=systematic[i];}
  Double_t GetUseDistortionFractionAsErrorY() const {return fDistortionFractionAsErrorYZ[0];}
  Double_t GetUseDistortionFractionAsErrorZ() const {return fDistortionFractionAsErrorYZ[1];}
  Double_t GetUseDistDispFractionAsErrorY() const {return fDistDispFractionAsErrorYZ[0];}
  Double_t GetUseDistDispFractionAsErrorZ() const {return fDistDispFractionAsErrorYZ[1];}
  void     SetUseDistortionFractionAsErrorY(double v) {fDistortionFractionAsErrorYZ[0] = v;}
  void     SetUseDistortionFractionAsErrorZ(double v) {fDistortionFractionAsErrorYZ[1] = v;}
  void     SetUseDistDispFractionAsErrorY(double v) {fDistDispFractionAsErrorYZ[0] = v;}
  void     SetUseDistDispFractionAsErrorZ(double v) {fDistDispFractionAsErrorYZ[1] = v;}
  const Double_t * GetSystematicError() const { return fSystematicErrors;}
  const Double_t * GetSystematicErrorClusterInner() const { return fSystematicErrorClusterInner;}
  const Double_t * GetSystematicErrorCluster() const { return fSystematicErrorCluster;}

  const TVectorF* GetSystErrClInnerRegZ()       const {return fSystErrClInnerRegZ;}
  const TVectorF* GetSystErrClInnerRegZSigInv() const {return fSystErrClInnerRegZSigInv;}
  void SetSystErrClInnerRegZ(TVectorF* zc)         {fSystErrClInnerRegZ = zc;}
  void SetSystErrClInnerRegZSigInv(TVectorF* zs)   {fSystErrClInnerRegZSigInv = zs;}
  
  void    SetUseSystematicCorrelation(Bool_t useCorrelation)  {fUseSystematicCorrelation=useCorrelation;}
  Bool_t  GetUseSystematicCorrelation() const { return fUseSystematicCorrelation;}

  static   AliTPCRecoParam *GetLowFluxParam();        // make reco parameters for low  flux env.
  static   AliTPCRecoParam *GetHighFluxParam();       // make reco parameters for high flux env.
  static   AliTPCRecoParam *GetHLTParam(); // special setting for HLT
  static   AliTPCRecoParam *GetLaserTestParam(Bool_t bPedestal);  // special setting for laser
  static   AliTPCRecoParam *GetCosmicTestParam(Bool_t bPedestal); // special setting for cosmic
  //
 protected:

  Int_t    fUseHLTClusters;  ///< allows usage of HLT clusters instead of RAW data
  Int_t    fUseHLTPreSeeding; ///< Usage of HLT pre-seeding
  Bool_t   fBClusterSharing; ///< allows or disable cluster sharing during tracking
  Double_t fCtgRange;        ///< +-fCtgRange is the ctg(Theta) window used for clusterization and tracking (MI)
  Double_t fMaxSnpTracker;   ///< max sin of local angle  - for TPC tracker
  Double_t fMaxSnpTrack;     ///< max sin of local angle  - for track
  Bool_t   fUseOuterDetectors; ///< switch - to use the outer detectors
  Double_t fMaxChi2TPCTRD;     ///< maximal allowed chi2 between the TRD in and TPC out to be accepted for refit
  Double_t fMaxChi2TPCITS;     ///< maximal allowed chi2 between the ITS in and TPC out to be accepted for backpropagation
  //
  // Outlier filtering configuration
  //
  Int_t   fUseOulierClusterFilter;  ///< swith to use outlier cluster filter

  Double_t fCutSharedClusters[2]; ///< cut value - maximal amount  of shared clusters
  Int_t fClusterMaxRange[2];   ///< neighborhood  - to define local maxima for cluster
  //
  //   clusterer parameters
  //
  Bool_t   fDumpSignal;      ///< Dump Signal information flag
  Int_t    fFirstBin;        ///< first time bin used by cluster finder
  Int_t    fLastBin;         ///< last time bin  used by cluster finder
  Bool_t   fBCalcPedestal;   ///< calculate Pedestal
  Bool_t   fBDoUnfold;       ///< do unfolding of clusters
  Float_t  fDumpAmplitudeMin; ///< minimal amplitude of signal to be dumped
  Float_t  fMaxNoise;        ///< maximal noise sigma on pad to be used in cluster finder
  Int_t    fUseOnePadCluster; ///< flag - use one pad cluster -0 not use >0 use
  Bool_t   fUseHLTOnePadCluster; ///< flag - use one HLT pad cluster for tracking
  Float_t  fMinMaxCutAbs;    ///< minimal amplitude at cluster maxima
  Float_t  fMinLeftRightCutAbs;  ///< minimal amplitude left right - PRF
  Float_t  fMinUpDownCutAbs;  ///< minimal amplitude up-down - TRF
  Float_t  fMinMaxCutSigma;    ///< minimal amplitude at cluster maxima
  Float_t  fMinLeftRightCutSigma;  ///< minimal amplitude left right - PRF
  Float_t  fMinUpDownCutSigma;  ///< minimal amplitude up-down - TRF
  //
  //
  Float_t  fMaxC;            ///< maximal curvature for tracking
  Bool_t   fBSpecialSeeding; ///< special seeding with big inclination angles allowed (for Cosmic and laser)
  Bool_t   fBKinkFinder;       ///< do kink finder reconstruction
  Double_t fKinkAngleCutChi2[2];   ///< angular cut for kinks
  Int_t    fLastSeedRowSec;     ///< Most Inner Row to make seeding for secondaries
  Int_t    fSeedGapPrim;   ///< seeding gap for primary tracks
  Int_t    fSeedGapSec;   ///< seeding gap for secondary tracks

  //
  // Correction switches
  //
  Int_t fUseFieldCorrection;     ///< use field correction
  Bool_t fUseComposedCorrection; ///< flag to use composed correction
  Int_t fUseRPHICorrection;      ///< use rphi correction
  Int_t fUseRadialCorrection;    ///< use radial correction
  Int_t fUseQuadrantAlignment;   ///< use quadrant alignment
  Int_t fUseSectorAlignment;     ///< use sector alignment
  Int_t fUseDriftCorrectionTime; ///< use drift correction time
  Int_t fUseDriftCorrectionGY;   ///< use drif correction global y
  Int_t fUseGainCorrectionTime;  ///< use gain correction time
  Int_t fUseExBCorrection;       ///< use ExB correction
  Bool_t fUseMultiplicityCorrectionDedx; ///< use Dedx multiplicity correction
  Bool_t fUseAlignmentTime;              ///< use time dependent alignment correction
  Int_t fUseIonTailCorrection;   ///< use ion tail correction
  Double_t fCrosstalkCorrection;   ///< crosstalk correction factor (fro each signal substracted by (mean signal in wite patch)xfCrosstalkCorrection) - Effect important only after removing oc capacitors in 2012
  Double_t fCrosstalkCorrectionMissingCharge;   ///< crosstalk correction factor - missing charge factor (from each signal substracted by (mean signal in wite patch)xfCrosstalkCorrection) - Effect important only after removing  capacitors in 2012
 //
  // dEdx switches
  //
  Bool_t   fUseTotCharge;          ///< switch use total or max charge
  Float_t fMinFraction;           ///< truncated mean - lower threshold
  Float_t fMaxFaction;            ///< truncated mean - upper threshold
  Int_t   fNeighborRowsDedx;      ///< number of neighboring rows to identify cluster below thres in dEdx calculation 0 -> switch off
  Int_t   fGainCorrectionHVandPTMode; ///< switch for the usage of GainCorrectionHVandPT (see AliTPCcalibDB::GetGainCorrectionHVandPT
  Int_t   fAccountDistortions;        ///< account for distortions in tracking
  Double_t fSkipTimeBins;        ///< number of time bins to be skiiped (corrupted signal druing gating opening)

  Bool_t fUseTOFCorrection;  ///< switch - kTRUE use TOF correction kFALSE - do not use
  //
  Bool_t fUseCorrectionMap;  ///< flag to use parameterized correction map (AliTPCChebCorr)
  Int_t  fCorrMapTimeDepMethod; ///< method used for correction time dependence
  Int_t  fUseLumiType;          ///< luminosity graph to be used for different lumi scalings
  //  misscalibration
  //
  TVectorF* fSystErrClInnerRegZ;        //< center of region in Z to apply extra systematic error
  TVectorF* fSystErrClInnerRegZSigInv;  //< inverse sigma forgaussian dumping aroung this center Z to apply extra syst error  
  Double_t fSystematicErrors[5];  ///< systematic errors in the track parameters - to be added to TPC covariance matrix
  Double_t fSystematicErrorClusterInner[2];  ///< systematic error of the cluster - used to downscale the information

  Double_t fSystematicErrorCluster[2];        ///< systematic error of the cluster - used e.g in OpenGG run to provide better cluster to track association efficiency
  Double_t fDistortionFractionAsErrorYZ[2];   ///< use fraction of distortion as additional error
  Double_t fDistDispFractionAsErrorYZ[2];   ///< use fraction of distortion dispersion as additional error
  Bool_t fUseSystematicCorrelation;         ///< switch to use the correlation for the sys


public:
  static Bool_t fgUseTimeCalibration; ///< flag usage the time dependent calibration
                                      // to be switched off for pass 0 reconstruction
                                      // Use static function, other option will be to use
                                      // additional specific storage ?
  /// \cond CLASSIMP
  ClassDef(AliTPCRecoParam, 28)
  /// \endcond
};


#endif
