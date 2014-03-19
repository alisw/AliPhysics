#ifndef AliThreePionRadii_cxx
#define AliThreePionRadii_cxx
//
// Class AliThreePionRadii
//
// AliThreePionRadii
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
class TRandom3;

class AliESDEvent;
class AliAODEvent;
class AliESDtrackCuts;
class AliESDpid;

#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliChaoticityEventCollection.h"

class AliThreePionRadii : public AliAnalysisTaskSE {
 public:

  AliThreePionRadii();
  AliThreePionRadii(const Char_t *name);
  virtual ~AliThreePionRadii();
  AliThreePionRadii(const AliThreePionRadii &obj); 
  AliThreePionRadii &operator=(const AliThreePionRadii &obj);
  

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  enum {
    kPairLimit = 15000,//15000
    kNormPairLimit = 45000,//45000
    kMultLimitPbPb = 2000,//2000
    kMultLimitPP = 300,//300
    kMCarrayLimit = 110000,//110000
    kQbins = 20,
    kQbinsWeights = 40,
    kQbinsPP = 50,// 40
    kQbinsWeightsPP = 40,
    kNDampValues = 16,
    kRmin = 2,// min radii for Momentum resolution calculations
    kDENtypes = 1,// was (kRVALUES)*kNDampValues
    kSCLimit2 = 1,// 1, 6
    kSCLimit3 = 1// 1, 10
  };

  static const Int_t fEDbins   = 3;
  static const Int_t fCentBins = 20;// 0-100% PbPb, pPb, pp
  static const Int_t fRVALUES  = 10;// 


  Int_t GetNumRValues() const {return AliThreePionRadii::fRVALUES;}
  Int_t GetNumCentBins() const {return AliThreePionRadii::fCentBins;}
  Int_t GetNumEDBins() const {return AliThreePionRadii::fEDbins;}
  void SetFSICorrelations(Bool_t legoCase=kTRUE, TH1D *temp1DSS[10]=0x0, TH1D *temp1DOS[10]=0x0);
  //
  void SetMCdecision(Bool_t mc) {fMCcase = mc;}
  void SetPbPbCase(Bool_t pbpb) {fPbPbcase = pbpb;}
  void SetGenerateSignal(Bool_t gen) {fGenerateSignal = gen;}
  void SetNumKt3Bins(Int_t kt3bins) {fKt3bins = kt3bins;}
  void SetV0Mbinning(Bool_t V0Mbinning) {fV0Mbinning = V0Mbinning;}
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
  void SetTriggerType(Int_t tt) {fTriggerType = tt;}
  //


 private:

  void ParInit();
  Bool_t AcceptPair(AliChaoticityTrackStruct*, AliChaoticityTrackStruct*);
  Float_t GamovFactor(Int_t, Int_t, Float_t);
  void Shuffle(Int_t*, Int_t, Int_t);
  Short_t FillIndex2part(Short_t);
  Short_t FillIndex3part(Short_t);
  Short_t SetQcutBin(Short_t);
  Short_t SetNormBin(Short_t);
  void SetFillBins2(Short_t, Short_t, Short_t, Int_t, Int_t, Int_t&, Int_t&);
  void SetFillBins3(Short_t, Short_t, Short_t, Short_t, Int_t, Int_t, Int_t, Short_t, Int_t&, Int_t&, Int_t&, Bool_t&, Bool_t&, Bool_t&);
  void ArrangeQs(Short_t, Short_t, Short_t, Short_t, Int_t, Int_t, Int_t, Float_t, Float_t, Float_t, Short_t, Short_t, Float_t&, Float_t&, Float_t&);
  Float_t GetQinv(Short_t, Float_t[], Float_t[]);
  void GetQosl(Float_t[], Float_t[], Float_t&, Float_t&, Float_t&);
  Float_t FSICorrelation2(Int_t, Int_t, Float_t);
  Float_t MCWeight(Int_t, Int_t, Float_t, Int_t, Float_t);
  Float_t MCWeight3D(Bool_t, Int_t, Int_t, Float_t, Float_t, Float_t);
  //
  
  
  const char* fname;// name of class
  AliAODEvent            *fAOD; //!    // AOD object
  TList                  *fOutputList; //! Compact Output list
  AliPIDResponse         *fPIDResponse; //! PID response object; equivalent to AliAODpidUtil
  
  
  AliChaoticityEventCollection ***fEC; //!
  AliChaoticityEventStruct *fEvt; //!
  AliChaoticityTrackStruct *fTempStruct; //!
  TRandom3* fRandomNumber; //!

  
 
  //////////////////////////////
  // histogram structures

  struct St6 {
    TH1D *fNorm3; //!
    TH3D *fTerms3; //!
    TH1D *fTermsQ3; //!
    TH1D *fIdeal; //!
    TH1D *fSmeared; //!
    TH1D *fMeanKt; //!
  };
  struct St5 {
    TH2D *fExplicit2; //!
    TH2D *fExplicit2QW; //!
    TProfile2D *fAvgP; //!
    TH2D *fIdeal; //!
    TH2D *fSmeared; //!
    TH1D *fMeanKt; //!
    //
    TH1D *fMCqinv; //!
    TH1D *fMCqinvQW; //!
    TH2D *fPIDpurityDen; //!
    TH3D *fPIDpurityNum; //!
  };
  struct St_EDB {// SC structure
    struct St5 TwoPT[2];
    struct St6 ThreePT[5];
  };
  struct St_M {
    struct St_EDB EDB[fEDbins];
  };
  struct St4 {
    struct St_M MB[fCentBins];
  };
  struct St3 {
    struct St4 SC[kSCLimit3];
  };
  struct St2 {
    struct St3 Charge3[2];
    struct St4 SC[kSCLimit2];
  };
  struct St1 {
    struct St2 Charge2[2];
  };
  struct St1 Charge1[2];//!


 
  Bool_t fLEGO;
  Bool_t fMCcase;
  Bool_t fAODcase;
  Bool_t fPbPbcase;
  Bool_t fGenerateSignal;
  Bool_t fGeneratorOnly;
  Bool_t fPdensityPairCut;
  Int_t fRMax;
  UInt_t fFilterBit;
  Float_t fMaxChi2NDF;
  Int_t fMinTPCncls;
  Double_t fBfield;
  Int_t fMbin;
  Int_t fFSIindex;
  Int_t fEDbin;
  Int_t fMbins;
  Int_t fMultLimit;  
  Int_t fKt3bins;
  Bool_t fV0Mbinning;
  Int_t fCentBinLowLimit;
  Int_t fCentBinHighLimit;
  Int_t fTriggerType;
  Int_t fEventCounter;
  Int_t fEventsToMix;
  Int_t fZvertexBins;
  Int_t fMultLimits[fCentBins+1];
  Float_t fQcut[3];
  Float_t fQLowerCut;
  Float_t fQlimitC2;
  Int_t fQbinsC2;
  Float_t fNormQcutLow[3];
  Float_t fNormQcutHigh[3];
  Float_t fKupperBound;
  Float_t fQupperBound;
  Int_t   fQbins;
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

  Bool_t fDummyB;

  
  Char_t fDefaultsCharMult[kMultLimitPbPb];//!
  Char_t fDefaultsCharSE[kPairLimit];//!
  Char_t fDefaultsCharME[2*kPairLimit];//!
  Int_t fDefaultsInt[kMultLimitPbPb];//!
  TArrayI *fPairLocationSE[kMultLimitPbPb];//!
  TArrayI *fPairLocationME[kMultLimitPbPb];//!
  TArrayC *fTripletSkip1[kPairLimit];//!
  TArrayC *fTripletSkip2[2*kPairLimit];//!
  TArrayI *fOtherPairLocation1[2][kMultLimitPbPb];//!
  TArrayI *fOtherPairLocation2[2][kMultLimitPbPb];//!
  TArrayC *fNormPairSwitch[3][kMultLimitPbPb];//!
  TArrayC *fPairSplitCut[4][kMultLimitPbPb];//!
  
  AliChaoticityNormPairStruct *fNormPairs[3];//!
  
 public:
  TH1D *fFSI2SS[10];
  TH1D *fFSI2OS[10];

  ClassDef(AliThreePionRadii, 1); 
};

#endif
