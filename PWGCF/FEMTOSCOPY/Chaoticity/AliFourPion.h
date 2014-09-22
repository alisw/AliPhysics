#ifndef AliFourPion_cxx
#define AliFourPion_cxx
//
// Class AliFourPion
//
// AliFourPion
// author:
//        Dhevan Gangadharan (dhevan.raja.gangadharan@cern.ch)
//



class TH1F;
class TH3F;
class TH1D;
class TH2D;
class TH3D;

class TProfile;
class TProfile2D;
class TProfile3D;
class TRandom3;

class AliESDEvent;
class AliAODEvent;
class AliESDtrackCuts;
class AliESDpid;

#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliFourPionEventCollection.h"
#include "AliCentrality.h"

class AliFourPion : public AliAnalysisTaskSE {
 public:

  AliFourPion();
  AliFourPion(const Char_t *name);
  virtual ~AliFourPion();
  AliFourPion(const AliFourPion &obj); 
  AliFourPion &operator=(const AliFourPion &obj);
  

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  enum {
    kPairLimit = 10000,//10000
    kNormPairLimit = 45000,
    kMultLimitPbPb = 1800,//1800
    kMultLimitpp = 300,
    kMultBinspp = 10,
    kMCarrayLimit = 150000,// 110000
    kQbinsWeights = 40,
    kNDampValues = 16,
    kRmin = 5,// EW min radii 5 fm
    kDENtypes = 104,// was 44
  };

  static const Int_t fKbinsT   = 4;// Set fKstep as well !!!!
  static const Int_t fKbinsY   = 1;// Set fKstep as well !!!!
  static const Int_t fEDbins   = 2;
  static const Int_t fCentBins = 10;// 0-50%
  static const Int_t fMbinsMixing = 10;// 5% widths
  static const Int_t fRVALUES  = 7;// 7 EW radii (5-11) , was 8 Gaussian radii (3-10fm)


  Int_t GetNumKtBins() const {return AliFourPion::fKbinsT;}
  Int_t GetNumRValues() const {return AliFourPion::fRVALUES;}
  Int_t GetNumCentBins() const {return AliFourPion::fCentBins;}
  Int_t GetNumEDBins() const {return AliFourPion::fEDbins;}
  void SetWeightArrays(Bool_t legoCase=kTRUE, TH3F *histos[AliFourPion::fKbinsT][AliFourPion::fCentBins]=0x0);
  void SetMomResCorrections(Bool_t legoCase=kTRUE, TH2D *temp2DSC=0x0, TH2D *temp2DMC=0x0);
  void SetFSICorrelations(Bool_t legoCase=kTRUE, TH1D *tempss[12]=0x0, TH1D *tempos[12]=0x0);
  void SetMuonCorrections(Bool_t legoCase=kTRUE, TH2D *tempMuon=0x0);
  void Setc3FitEAs(Bool_t legoCase=kTRUE, TH3D *histoPbPb=0x0, TH3D *histopPb=0x0, TH3D *histopp=0x0);
  //
  void SetMCdecision(Bool_t mc) {fMCcase = mc;}
  void SetTabulatePairs(Bool_t tabulate) {fTabulatePairs = tabulate;}
  void SetInterpolationType(Bool_t linearInterp) {fLinearInterpolation = linearInterp;}
  void SetCollisionType(Bool_t ct) {fCollisionType = ct;}
  void SetGenerateSignal(Bool_t gen) {fGenerateSignal = gen;}
  void SetGeneratorOnly(Bool_t genOnly) {fGeneratorOnly = genOnly;}
  void SetCentBinRange(Int_t low, Int_t high) {fCentBinLowLimit = low; fCentBinHighLimit = high;}
  void SetLEGOCase(Bool_t lego) {fLEGO = lego;}
  void SetFilterBit(UInt_t filterbit) {fFilterBit = filterbit;}
  void SetMaxChi2NDF(Float_t MaxChi2NDF) {fMaxChi2NDF = MaxChi2NDF;}
  void SetMinTPCncls(Int_t MinTPCncls) {fMinTPCncls = MinTPCncls;}
  void SetPairSeparationCutEta(Float_t pairsep) {fMinSepPairEta = pairsep;}
  void SetPairSeparationCutPhi(Float_t pairsep) {fMinSepPairPhi = pairsep;}
  void SetNsigmaTPC(Float_t nsig) {fSigmaCutTPC = nsig;}
  void SetNsigmaTOF(Float_t nsig) {fSigmaCutTOF = nsig;}
  void SetRMax(Int_t rbin) {fRMax = rbin;}
  void SetfcSq(Float_t fcSq) {ffcSq = fcSq;}
  void SetMixedChargeCut(Bool_t mcCut) {fMixedChargeCut = mcCut;}
  void SetMinPt(Float_t minPt) {fMinPt = minPt;}
  void SetMaxPt(Float_t maxPt) {fMaxPt = maxPt;}
  void SetKT3transition(Float_t KT3trans) {fKT3transition = KT3trans;}
  void SetKT4transition(Float_t KT4trans) {fKT4transition = KT4trans;}
  void SetTriggerType(Int_t tt) {fTriggerType = tt;}
  //
 

 private:

  void ParInit();
  Bool_t AcceptPair(AliFourPionTrackStruct, AliFourPionTrackStruct);
  Bool_t AcceptPairPM(AliFourPionTrackStruct, AliFourPionTrackStruct);
  Float_t Gamov(Int_t, Int_t, Float_t);
  void Shuffle(Int_t*, Int_t, Int_t);
  Float_t GetQinv(Float_t[], Float_t[]);
  void GetQosl(Float_t[], Float_t[], Float_t&, Float_t&, Float_t&);
  void GetWeight(Float_t[], Float_t[], Float_t&, Float_t&);
  Float_t FSICorrelation(Int_t, Int_t, Float_t);
  Float_t MCWeight(Int_t[2], Float_t, Float_t, Float_t, Float_t);
  Float_t MCWeightOSL(Int_t, Int_t, Int_t, Int_t, Float_t, Float_t, Float_t, Float_t);
  Float_t MCWeight3(Int_t, Float_t, Float_t, Int_t[3], Float_t[3], Float_t[3]);
  Float_t MCWeightFSI3(Int_t, Float_t, Float_t, Int_t[3], Float_t[3]);
  Float_t MCWeight4(Int_t, Float_t, Float_t, Int_t[4], Float_t[6], Float_t[6]);
  Float_t MCWeightFSI4(Int_t, Float_t, Float_t, Int_t[4], Float_t[6]);
  //
  void SetFillBins2(Int_t, Int_t, Int_t&, Int_t&);
  void SetFillBins3(Int_t, Int_t, Int_t, Short_t, Int_t&, Int_t&, Int_t&, Bool_t&, Bool_t&, Bool_t&);
  void SetFillBins4(Int_t, Int_t, Int_t, Int_t, Int_t&, Int_t&, Int_t&, Int_t&, Int_t, Bool_t[13]);
  void SetFSIindex(Float_t);
  //
  Float_t cubicInterpolate(Float_t[4], Float_t);
  Float_t nCubicInterpolate(Int_t, Float_t*, Float_t[]);
  
  
  const char* fname;// name of class
  AliAODEvent            *fAOD; //!    // AOD object
  TList                  *fOutputList; //! Compact Output list
  AliPIDResponse         *fPIDResponse; //! PID response object; equivalent to AliAODpidUtil
  
  
  AliFourPionEventCollection ***fEC; //!
  AliFourPionEventStruct *fEvt; //!
  AliFourPionTrackStruct *fTempStruct; //!
  TRandom3* fRandomNumber; //!

  
 
  //////////////////////////////
  // histogram structures

  
  struct St6 {
    TH1D *fNorm3; //!
    TH1D *fTerms3; //!
    TH3D *fTerms33D; //!
    TProfile *fKfactor; //!
    TProfile3D *fKfactor3D; //!
    TProfile *fKfactorWeighted; //!
    TProfile *fMeanQinv; //!
    TH2D *fIdeal; //!
    TH2D *fSmeared; //!
    //
    TH3D *fMuonSmeared; //!
    TH3D *fMuonIdeal; //!
    TH3D *fMuonPionK3; //!
    TH3D *fPionPionK3; //!
    //
    TH2D *fTwoPartNorm; //!
    TH2D *fTwoPartNegNorm; //!
    TH2D *fTwoPartNormErr; //!
  };
  struct St7 {
    TH3D *fTerms2OSL; //!
    TH3D *fTerms2OSLQW; //!
  };
  struct St5 {
    TH2D *fTerms2; //!
    TH2D *fTerms2QW; //!
    TH3D *fTerms2ThreeD; //!
    TProfile2D *fAvgP; //!
    TH2D *fIdeal; //!
    TH2D *fSmeared; //!
    TH2D *fUnitMultBin; //!
    //
    TH2D *fMuonSmeared; //!
    TH2D *fMuonIdeal; //!
    TH2D *fMuonPionK2; //!
    TH2D *fPionPionK2; //!
    TH1D *fMCqinv; //!
    TH1D *fMCqinvQW; //!
    TH2D *fPIDpurityDen; //!
    TH3D *fPIDpurityNum; //!
    struct St7 OSL_ktbin[2];
  };
  struct StFourPT {
    TH1D *fNorm4; //!
    TH1D *fTerms4; //!
    TProfile *fKfactor; //!
    TProfile *fKfactorWeighted; //!
    TH2D *fIdeal; //!
    TH2D *fSmeared; //!
    //
    TH3D *fMuonSmeared; //!
    TH3D *fMuonIdeal; //!
    TH3D *fMuonPionK4; //!
    TH3D *fPionPionK4; //!
    //
    TH2D *fTwoPartNorm; //!
    TH2D *fTwoPartNegNorm; //!
    TH2D *fTwoPartNormErr; //!
    TH3D *fFullBuildFromFits; //!
    TH3D *fPartialBuildFromFits; //!
  };
  struct St_EDB {
    struct St5 TwoPT[2];
    struct St6 ThreePT[5];
    struct StFourPT FourPT[13];
  };
  struct St_M {
    struct St_EDB EDB[fEDbins];
  };
  struct St4 {
    struct St_M MB[fCentBins];
  };
  struct St3 {
    struct St4 Charge4[2];
    struct St_M MB[fCentBins];
  };
  struct St2 {
    struct St3 Charge3[2];
    struct St_M MB[fCentBins];
  };
  struct St1 {
    struct St2 Charge2[2];
  };
  struct St1 Charge1[2];//!


  /////////////////////
  // 4D r3 denominator
  struct St_Ky {
    struct St_M MB[fCentBins];
  };
  struct St_Kt {
    struct St_Ky KY[fKbinsY];
  };
  struct St_Kt KT[fKbinsT];//!
  
 
  Bool_t fLEGO;
  Bool_t fMCcase;
  Bool_t fAODcase;
  Short_t fCollisionType;
  Bool_t fGenerateSignal;
  Bool_t fGeneratorOnly;
  Bool_t fTabulatePairs;
  Bool_t fLinearInterpolation;
  Bool_t fMixedChargeCut;
  Int_t fRMax;
  Float_t fRstartMC;
  Float_t ffcSq;
  Float_t ffcSqMRC;
  UInt_t fFilterBit;
  Float_t fMaxChi2NDF;
  Int_t fMinTPCncls;
  Double_t fBfield;
  Int_t fMbin;
  Int_t fFSIindex;
  Int_t fEDbin;
  Int_t fMbins;
  Int_t fMultLimit;      
  Int_t fCentBinLowLimit;
  Int_t fCentBinHighLimit;
  Int_t fTriggerType;
  Int_t fEventCounter;
  Int_t fEventsToMix;
  Int_t fZvertexBins;
  Int_t fMultLimits[kMultBinspp+1];
  Float_t fMinPt;
  Float_t fMaxPt;
  Float_t fQcut;
  Float_t fQLowerCut;
  Float_t fNormQcutLow;
  Float_t fNormQcutHigh;
  Float_t fKupperBound;
  Float_t fQupperBoundQ2;
  Float_t fQupperBoundQ3;
  Float_t fQupperBoundQ4;
  Float_t fQbinsQ2;
  Float_t fQbinsQ3;
  Float_t fQbinsQ4;
  Float_t fQupperBoundWeights;
  Float_t fQbinsQinv3D;
  Float_t fQupperBoundQinv3D;
  Float_t fKstepT[fKbinsT];
  Float_t fKstepY[fKbinsY];
  Float_t fKmeanT[fKbinsT];
  Float_t fKmeanY[fKbinsY];
  Float_t fKmiddleT[fKbinsT];
  Float_t fKmiddleY[fKbinsY];
  Float_t fQstep;
  Float_t fQstepWeights;
  Float_t fQmean[kQbinsWeights];
  Float_t fDampStart;
  Float_t fDampStep;
  
  Float_t fTPCTOFboundry;
  Float_t fTOFboundry;
  Float_t fSigmaCutTPC;
  Float_t fSigmaCutTOF;
  
  Float_t fMinSepPairEta;
  Float_t fMinSepPairPhi;
  Float_t fShareQuality;
  Float_t fShareFraction;
  
  Float_t fTrueMassP, fTrueMassPi, fTrueMassK, fTrueMassKs, fTrueMassLam;
 
  Int_t fKtIndexL,fKtIndexH;
  //
  Int_t fQoIndexL,fQoIndexH;
  Int_t fQsIndexL,fQsIndexH;
  Int_t fQlIndexL,fQlIndexH;

  Bool_t fDummyB;

  Float_t fKT3transition;
  Float_t fKT4transition;
  
  Float_t farrP1[4][4][4];
  Float_t farrP2[4][4][4];
  
 
  //
  Char_t fDefaultsCharSwitch[kMultLimitPbPb];//!
  TArrayC *fLowQPairSwitch_E0E0[kMultLimitPbPb];//!
  TArrayC *fLowQPairSwitch_E0E1[kMultLimitPbPb];//!
  TArrayC *fLowQPairSwitch_E0E2[kMultLimitPbPb];//!
  TArrayC *fLowQPairSwitch_E0E3[kMultLimitPbPb];//!
  TArrayC *fLowQPairSwitch_E1E1[kMultLimitPbPb];//!
  TArrayC *fLowQPairSwitch_E1E2[kMultLimitPbPb];//!
  TArrayC *fLowQPairSwitch_E1E3[kMultLimitPbPb];//!
  TArrayC *fLowQPairSwitch_E2E3[kMultLimitPbPb];//!
  //
  TArrayC *fNormQPairSwitch_E0E0[kMultLimitPbPb];//!
  TArrayC *fNormQPairSwitch_E0E1[kMultLimitPbPb];//!
  TArrayC *fNormQPairSwitch_E0E2[kMultLimitPbPb];//!
  TArrayC *fNormQPairSwitch_E0E3[kMultLimitPbPb];//!
  TArrayC *fNormQPairSwitch_E1E1[kMultLimitPbPb];//!
  TArrayC *fNormQPairSwitch_E1E2[kMultLimitPbPb];//!
  TArrayC *fNormQPairSwitch_E1E3[kMultLimitPbPb];//!
  TArrayC *fNormQPairSwitch_E2E3[kMultLimitPbPb];//!

  

 public:
  TH2D *fMomResC2SC;
  TH2D *fMomResC2MC;
  TH2D *fWeightmuonCorrection;
  TH3D *fPbPbc3FitEA;
  TH3D *fpPbc3FitEA;
  TH3D *fppc3FitEA;
  TH1D *fFSIss[12];
  TH1D *fFSIos[12];
  TH3F *fNormWeight[fKbinsT][fCentBins];
  TF1 *ExchangeAmpPointSource[2][50];

  ClassDef(AliFourPion, 1); 
};

#endif
