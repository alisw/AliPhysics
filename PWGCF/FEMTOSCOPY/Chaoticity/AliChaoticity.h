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
  AliChaoticity(const Char_t *name, Bool_t MCdecision=kFALSE, Bool_t Tabulatedecision=kFALSE, Bool_t PbPbdecision=kTRUE, Int_t lowCentBin=0, Int_t highCentBin=1.,  Bool_t lego=kTRUE);
  virtual ~AliChaoticity();
  AliChaoticity(const AliChaoticity &obj); 
  AliChaoticity &operator=(const AliChaoticity &obj);
  
 private:

  virtual void   UserCreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 enum {
    kPairLimit = 15000,//15000
    kNormPairLimit = 45000,
    kMultLimitPbPb = 2000,//2000
    kMultLimitpp = 300,
    kMultBinspp = 11,//20 or 11
    kKbinsT = 3,// Set fKstep as well !!!!
    kKbinsY = 1,// Set fKstep as well !!!!
    kQbins = 20,
    kQbinsWeights = 40,
    kEDbins = 1,
    kRVALUES = 8,
    kNDampValues = 16,// change to 11 soon
    kDENtypes = (kRVALUES)*kNDampValues,
    kCentBins=10,// 0-50%
    kSCLimit2 = 1,// 1, 6
    kSCLimit3 = 1,// 1, 10
    kNlinesCoulFile = 99
  };

  void ParInit();
  Bool_t AcceptPair(AliChaoticityTrackStruct, AliChaoticityTrackStruct);
  Float_t GamovFactor(Int_t, Int_t, Float_t);
  void Shuffle(Int_t*, Int_t, Int_t);
  short FillIndex2part(short);
  short FillIndex3part(short);
  short SetQcutBin(short);
  short SetNormBin(short);
  void SetFillBins2(short, short, short, Int_t, Int_t, Int_t&, Int_t&);
  void SetFillBins3(short, short, short, short, Int_t, Int_t, Int_t, short, Int_t&, Int_t&, Int_t&, Bool_t&, Bool_t&, Bool_t&);
  void ArrangeQs(short, short, short, short, Int_t, Int_t, Int_t, Float_t, Float_t, Float_t, short, short, Float_t&, Float_t&, Float_t&);
  Float_t GetQinv(short, Float_t[], Float_t[]);
  void GetQosl(Float_t[], Float_t[], Float_t&, Float_t&, Float_t&);
  void SetWeightArrays();
  void SetWeightArraysLEGO(TH3F *histos[kKbinsT][kCentBins]);
  void SetMomResCorrections();
  void SetMomResCorrectionsLEGO(TH2D *histo);
  Float_t GetMomRes(Int_t, Float_t);
  void GetWeight(Float_t[], Float_t[], Float_t&, Float_t&);
  Float_t CoulCorr(Int_t, Int_t, Int_t, Float_t);
  void SetCoulCorrections();
  void SetCoulCorrectionsLEGO(Float_t[kNlinesCoulFile], Float_t[kRVALUES][kNlinesCoulFile], Float_t[kRVALUES][kNlinesCoulFile]);
  Float_t MCWeight(Int_t, Int_t, Int_t, Int_t, Float_t);
  Float_t MCWeightOSL(Float_t, Float_t, Float_t, Float_t, Float_t, Float_t);
  Float_t MCWeightr3(Int_t, Int_t, Int_t, Float_t, Float_t, Float_t);
  //
  Int_t GetNumKtbins() const {return kKbinsT;}
  Int_t GetNumCoulLines() const {return kNlinesCoulFile;}
  Int_t GetNumRValues() const {return kRVALUES;}
  Int_t GetNumCentBins() const {return kCentBins;}

  
  const char* fname;// name of class
  AliAODEvent            *fAOD; //!    // AOD object
  AliESDEvent            *fESD; //!    // ESD object
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
  };  
  struct St6 {
    TH1D *fExplicit3; //!
    TH1D *fNormEx3; //!
    //
    TH1D *fNorm3; //!
    TH3D *fTerms3; //!
    //
    struct St_DT DT[kDENtypes];
  };
  struct St7 {
    TH3D *fExplicit2OSL; //!
    TH3D *fExplicit2OSLQW; //!
  };
  struct St5 {
    TH2D *fExplicit2; //!
    TH3I *fExplicit2ThreeD; //!
    TH2D *fIdeal; //!
    TH2D *fSmeared; //!
    struct St7 OSL_ktbin[2];
  };
  struct St_EDB {// SC structure
    struct St5 TwoPT[2];
    struct St6 ThreePT[5];
  };
  struct St_M {
    struct St_EDB EDB[kEDbins];
  };
  struct St4 {
    struct St_M MB[kCentBins];
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
  struct St1 Charge1[2];


  /////////////////////
  // 4D r3 denominator
  struct St_Ky {
    struct St_M MB[kCentBins];
  };
  struct St_Kt {
    struct St_Ky KY[kKbinsY];
  };
  struct St_Kt KT[kKbinsT];
  
 
  Bool_t fLEGO;
  Bool_t fMCcase;
  Bool_t fAODcase;
  Bool_t fPbPbcase;
  Bool_t fPdensityExplicitLoop;
  Bool_t fPdensityPairCut;
  Bool_t fTabulatePairs;
  Double_t fBfield;
  Int_t fMbin;
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
  Float_t fKstepT[kKbinsT];
  Float_t fKstepY[kKbinsY];
  Float_t fKmeanT[kKbinsT];
  Float_t fKmeanY[kKbinsY];
  Float_t fKmiddleT[kKbinsT];
  Float_t fKmiddleY[kKbinsY];
  Float_t fQstep;
  Float_t fQmean[kQbinsWeights];
  Float_t fDampStart;
  Float_t fDampStep;
  
  Float_t fQCoul[100];//! 2 MeV bins
  Float_t fCoulSS[kRVALUES][100];//! Radii, Q
  Float_t fCoulOS[kRVALUES][100];//! Radii, Q
  Float_t fMomResWeights[kDENtypes][kQbinsWeights];//!


  Float_t fTPCTOFboundry;
  Float_t fTOFboundry;
  Float_t fSigmaCutTPC;
  Float_t fSigmaCutTOF;
  
  Float_t fMinSepTPCEntrancePhi;
  Float_t fMinSepTPCEntranceEta;
  Float_t fShareQuality;
  Float_t fShareFraction;
  
  Float_t fTrueMassP, fTrueMassPi, fTrueMassK, fTrueMassKs, fTrueMassLam;
 
  Int_t fKtbinL,fKtbinH;
  Int_t fKybinL,fKybinH;
  //
  Int_t fQobinL,fQobinH;
  Int_t fQsbinL,fQsbinH;
  Int_t fQlbinL,fQlbinH;

  Bool_t fDummyB;

  Float_t *******fNormWeight;//! osl kt binning
  Float_t *******fNormWeightErr;//! osl kt binning

  
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
 
 
   
  ClassDef(AliChaoticity, 1); 
};

#endif
