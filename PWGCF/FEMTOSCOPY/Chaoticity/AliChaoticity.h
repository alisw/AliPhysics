#ifndef AliChaoticity_cxx
#define AliChaoticity_cxx
//
// Class AliChaoticity
//
// AliChaoticity
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
#include "AliCentrality.h"

class AliChaoticity : public AliAnalysisTaskSE {
 public:

  AliChaoticity();
  AliChaoticity(const Char_t *name);
  virtual ~AliChaoticity();
  AliChaoticity(const AliChaoticity &obj); 
  AliChaoticity &operator=(const AliChaoticity &obj);
  

  virtual void   UserCreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  enum {
    kPairLimit = 15000,//15000
    kNormPairLimit = 45000,
    kMultLimitPbPb = 2000,//2000
    kMultLimitpp = 300,
    kMultBinspp = 11,//20 or 11
    kQbins = 20,
    kQbinsWeights = 40,
    kNDampValues = 16,
    kDENtypes = 1,// was (kRVALUES)*kNDampValues
    kSCLimit2 = 1,// 1, 6
    kSCLimit3 = 1// 1, 10
  };

  static const Int_t fKbinsT   = 3;// Set fKstep as well !!!!
  static const Int_t fKbinsY   = 1;// Set fKstep as well !!!!
  static const Int_t fEDbins   = 1;
  static const Int_t fCentBins = 10;// 0-50%
  static const Int_t fRVALUES  = 8;// 8 Gaussian radii (3-10fm)


  Int_t GetNumKtBins() const {return AliChaoticity::fKbinsT;}
  Int_t GetNumRValues() const {return AliChaoticity::fRVALUES;}
  Int_t GetNumCentBins() const {return AliChaoticity::fCentBins;}
  Int_t GetNumEDBins() const {return AliChaoticity::fEDbins;}
  void SetWeightArrays(Bool_t legoCase=kTRUE, TH3F *histos[AliChaoticity::fKbinsT][AliChaoticity::fCentBins]=0x0);
  void SetMomResCorrections(Bool_t legoCase=kTRUE, TH2D *temp2D=0x0);
  void SetFSICorrelations(Bool_t legoCase=kTRUE, TH2D *temp2DGaus[2]=0x0, TH2D *temp2DTherm[6]=0x0, TH3D *temp3Dos[6]=0x0, TH3D *temp3Dss[6]=0x0);
  //
  void SetMCdecision(Bool_t mc) {fMCcase = mc;}
  void SetTabulatePairs(Bool_t tabulate) {fTabulatePairs = tabulate;}
  void SetPbPbCase(Bool_t pbpb) {fPbPbcase = pbpb;}
  void SetGenerateSignal(Bool_t gen) {fGenerateSignal = gen;}
  void SetCentBinRange(Int_t low, Int_t high) {fCentBinLowLimit = low; fCentBinHighLimit = high;}
  void SetLEGOCase(Bool_t lego) {fLEGO = lego;}
  void SetFilterBit(UInt_t filterbit) {fFilterBit = filterbit;}
  void SetPairSeparationCut(Float_t pairsep) {fMinSepPair = pairsep;}
  void SetNsigmaTPC(Float_t nsig) {fSigmaCutTPC = nsig;}
  void SetNsigmaTOF(Float_t nsig) {fSigmaCutTOF = nsig;}
  void SetRBinMax(Int_t rbin) {fRBinMax = rbin;}
  void SetFixedLambdaBin(Int_t lbin) {fFixedLambdaBin = lbin;}
  
  //


 private:

  void ParInit();
  Bool_t AcceptPair(AliChaoticityTrackStruct, AliChaoticityTrackStruct);
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
  //void GetWeight(Float_t[], Float_t[], Float_t&, Float_t&);
  void GetWeight(Float_t[], Float_t[], Float_t[], Float_t[], Float_t&, Float_t&);
  void FourVectProdTerms(Float_t [], Float_t [], Float_t [], Float_t&, Float_t&, Float_t&, Float_t&, Float_t&);
  Float_t FSICorrelationGaus2(Int_t, Int_t, Int_t, Float_t);
  Float_t FSICorrelationTherm2(Int_t, Int_t, Float_t);
  Float_t MCWeight(Int_t, Int_t, Int_t, Int_t, Float_t);
  Float_t MCWeightOSL(Int_t, Int_t, Int_t, Int_t, Float_t, Float_t, Float_t, Float_t);
  Float_t MCWeight3D(Bool_t, Int_t, Int_t, Float_t, Float_t, Float_t);
  Double_t FSICorrelationOmega0(Bool_t, Double_t, Double_t, Double_t);
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

  struct St_DT {
    TH3D *fTwoPartNorm; //!
    //TH3D *fTwoPartNormErr; //!
    TH3D *f4VectProd1TwoPartNorm; //!
    TH3D *f4VectProd2TwoPartNorm; //!
    TH3D *f4VectProd1TwoPartNormIdeal; //!
    TH3D *f4VectProd2TwoPartNormIdeal; //!
    TH3D *f4VectProd1TwoPartNormSmeared; //!
    TH3D *f4VectProd2TwoPartNormSmeared; //!
  };  
  struct St6 {
    TH1D *fExplicit3; //!
    TH1D *fNormEx3; //!
    //
    TH1D *fNorm3; //!
    TH3D *fTerms3; //!
    TH3D *f4VectProd1Terms; //!
    TH3D *f4VectProd2Terms; //!
    TH3D *f4VectProd1TermsIdeal; //!
    TH3D *f4VectProd2TermsIdeal; //!
    TH3D *f4VectProd1TermsSmeared; //!
    TH3D *f4VectProd2TermsSmeared; //!
    TH3D *f4VectProd1TermsSumK3; //!
    TH3D *f4VectProd2TermsSumK3; //!
    TH3D *f4VectProd1TermsEnK3; //!
    TH3D *f4VectProd2TermsEnK3; //!
    TH3D *f4VectProd1TermsSumK2; //!
    TH3D *f4VectProd2TermsSumK2; //!
    TH3D *f4VectProd1TermsEnK2; //!
    TH3D *f4VectProd2TermsEnK2; //!
    TH3D *fIdeal; //!
    TH3D *fSmeared; //!
    TH3D *fQW12; //!
    TH3D *fQW13; //!
    TH3D *fSumK3; //!
    TH3D *fEnK3; //!
    TH3D *f4VectProd1Q3W; //!
    TH3D *f4VectProd2Q3W; //!
    //
    struct St_DT DT[kDENtypes];
  };
  struct St7 {
    TH3D *fExplicit2OSL; //!
    TH3D *fExplicit2OSLQW; //!
  };
  struct St5 {
    TH2D *fExplicit2; //!
    TH2D *fExplicit2QW; //!
    TH3D *fExplicit2ThreeD; //!
    TProfile2D *fAvgP; //!
    TH2D *fIdeal; //!
    TH2D *fSmeared; //!
    //
    TH1D *fMCqinv; //!
    TH1D *fMCqinvQW; //!
    TH2D *fPIDpurityDen; //!
    TH2D *fPIDpurityNum; //!
    struct St7 OSL_ktbin[2];
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
  Bool_t fPbPbcase;
  Bool_t fGenerateSignal;
  Bool_t fPdensityExplicitLoop;
  Bool_t fPdensityPairCut;
  Bool_t fTabulatePairs;
  Int_t fRBinMax;// 5 normally, (R=8fm)
  Int_t fFixedLambdaBin;// 11 normally, (lambda=0.52)
  UInt_t fFilterBit;
  Double_t fBfield;
  Int_t fMbin;
  Int_t fFSIbin;
  Int_t fEDbin;
  Int_t fMbins;
  Int_t fMultLimit;      
  Int_t fCentBinLowLimit;
  Int_t fCentBinHighLimit;
  Int_t fEventCounter;
  Int_t fEventsToMix;
  Int_t fZvertexBins;
  Int_t fMultLimits[kMultBinspp+1];
  Float_t fQcut[3];
  Float_t fQLowerCut;
  Float_t fNormQcutLow[3];
  Float_t fNormQcutHigh[3];
  Float_t fKupperBound;
  Float_t fQupperBound;
  Float_t fQupperBoundWeights;
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
  
  Float_t fMinSepPair;
  Float_t fShareQuality;
  Float_t fShareFraction;
  
  Float_t fTrueMassP, fTrueMassPi, fTrueMassK, fTrueMassKs, fTrueMassLam;
 
  Int_t fKtIndexL,fKtIndexH;
  //
  Int_t fQoIndexL,fQoIndexH;
  Int_t fQsIndexL,fQsIndexH;
  Int_t fQlIndexL,fQlIndexH;

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
  TH2D *fFSI2SS[2];
  TH2D *fFSI2OS[2];
  TH3D *fFSIOmega0SS[6];
  TH3D *fFSIOmega0OS[6];
  TH2D *fMomResC2;
  TH3F *fNormWeight[3][10];
    

  ClassDef(AliChaoticity, 1); 
};

#endif
