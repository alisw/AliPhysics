/// \class AliAnalysisTaskKaonXiCorrelation

#ifndef __AliAnalysisTaskKaonXiCorrelation__
#define __AliAnalysisTaskKaonXiCorrelation__

#include "AliAnalysisTaskSE.h"
#include <Rtypes.h>
#include <TString.h>
#include "AliEventCuts.h"
#include "AliExternalBDT.h"
#include "AliVTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODcascade.h"

class AliPIDResponse;
class TH2F;
class TH3F;
class TList;
class TTree;

struct MiniXi {
  Double32_t fPt; //[-12.7,12.8,8]
  Double32_t fEta; //[-1.27,1.28,8]
  Double32_t fMass; //[1.29,1.35375,8]
  Double32_t fBdtOut; //[0.5,1.,16]
  unsigned char fRecFlag;
};

struct MiniXiMC : public MiniXi {
  Double32_t fPtMC; //[-12.7,12.8,8]
  Double32_t fEtaMC; //[-1.27,1.28,8]
  bool fIsReconstructed;
  unsigned char fFlag;
};

struct MiniKaon {
  Double32_t fPt; //[-3.175,3.2,8]
  Double32_t fEta; //[-1.27,1.28,8]
  Double32_t fNsigmaITS; //[-6.35,6.4,8]
  Double32_t fNsigmaTPC; //[-6.35,6.4,8]
  Double32_t fNsigmaTOF; //[-6.35,6.4,8]
  unsigned char fCutBitMap;
};

struct MiniKaonMC : public MiniKaon {
  Double32_t fPtMC; //[-3.175,3.2,8]
  Double32_t fEtaMC; //[-1.27,1.28,8]
  bool fIsReconstructed;
  unsigned char fFlag;
};

struct MiniCollision {
  Double32_t fZ; //[-12.7,12.8,8]
  unsigned char fCent;
  unsigned char fTrigger;
  unsigned short fNTrk;
  unsigned short fV0MAmp;
};

struct XiDaughter {
  XiDaughter(const double p_x, const double p_y, const double p_z, const int c) :
   px{p_x},
   py{p_y},
   pz{p_z},
   charge{c}
  {}
  bool IsSame(const XiDaughter other) { return std::abs(px - other.px) < 1.e-15 && std::abs(py - other.py) < 1.e-15 && std::abs(pz - other.pz) < 1.e-15 && charge == other.charge; }
  double px;
  double py;
  double pz;
  int charge;
};

class AliAnalysisTaskKaonXiCorrelation : public AliAnalysisTaskSE {
public:
  enum kStatusFlag {
    kPrimary = BIT(0),
    kSecondaryFromWD = BIT(1),
    kSecondaryFromMaterial = BIT(2)
  };

  enum kXiFlag {
    kHasTOFhitOrITSrefit = BIT(0),
    kCompetingMassCut = BIT(1)
  };

  enum kCutFlag {
    kDCAtightCut = BIT(0),   // 0.05 cm
    kDCAmidCut = BIT(1),     // 0.1 cm
    kTPCclsTightCut = BIT(2),// 80
    kTPCclsMidCut = BIT(3),  // 70
    kChi2TightCut = BIT(4),  // 2
    kChi2MidCut = BIT(5)     // 2.5
  };

  enum kReducedTrigger
  {
    kINT7 = BIT(0),
    kCentral = BIT(1),
    kSemiCentral = BIT(2),
    kPositiveB = BIT(3),
    kHasTwoXiFromSameDaughter = BIT(4)
  };

  AliAnalysisTaskKaonXiCorrelation(bool isMC = false, TString taskname = "KaonXiCorrelation");
  static AliAnalysisTaskKaonXiCorrelation* AddTask(bool isMC, TString tskname, TString suffix);
  virtual ~AliAnalysisTaskKaonXiCorrelation();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual Bool_t UserNotify();
  virtual void   Terminate(Option_t *) {}

  AliEventCuts  fEventCuts; ///<

  void SetFillCascades(bool toogle = true) { fFillCascades = toogle; }
  void SetApplyBdtToMC(bool toogle = true) { fApplyBdtToMC = toogle; }

  void SetEstimator(const int estimator = 0) { fEstimator = estimator; }

  // Setters for configurable cuts
  void SetRadiusCut(float cut = 1.2) { fCutRadius = cut; }
  void SetRadiusV0Cut(float cut = 3.0) { fCutRadiusV0 = cut; }
  void SetDCABachToPVCut(float cut = 0.1) { fCutDCABachToPV = cut; }
  void SetDCAV0toPVCut(float cut = 0.1) { fCutDCAV0toPV = cut; }
  void SetDCAV0piToPVCut(float cut = 0.2) { fCutDCAV0piToPV = cut; }
  void SetDCAV0prToPVCut(float cut = 0.2) { fCutDCAV0prToPV = cut; }
  void SetDCAV0tracksCut(float cut = 1.0) { fCutDCAV0tracks = cut; }
  void SetCutDCABachToV0XiCut(float cut = 1.0) { fCutDCABachToV0 = cut; }
  void SetCosPACut(float cut = 0.95) { fCutCosPA = cut; }
  void SetCosPAV0Cut(float cut = 0.95) { fCutCosPAV0 = cut; }
  void SetBachBarCosPACut(float cut = 0.99995) { fCutBachBarCosPA = cut; }
  void SetV0MassWindowCut(float cut = 0.005) { fCutV0MassWindow = cut; }
  void SetYCut(float cut = 0.5) { fCutY = cut; }
  void SetNsigmaTPCCut(float cut = 4.0) { fCutNsigmaTPC = cut; }
  void SetCtCut(float cut = 4) { fCutCt = cut; }
  void SetCtV0Cut(float cut = 30) { fCutCtV0 = cut; }
  void SetCompetingMassCut(float cut = 0.008) { fCutCompetingMass = cut; }
  void SetCascMassWindow(float cut = 0.006) { fCascMassWindow = cut; }
  void SetTPCcluCut(int cut = 70) { fCutTPCclu = cut; }
  void SetSaveOnlyTrueCandidates(bool cut = true) { fOnlyTrueCandidates = cut; }
  void SetTPCRowsCut(float cut = 80.) { fCutTPCrows = cut; }
  void SetTPCRowOvFCut(float cut = 0.8) { fCutRowsOvF = cut; }
  void UseOnTheFly(bool toggle = true) { fUseOnTheFly = toggle; }
  void SetMinCentrality(int minCentrality = 0) { fMinCentrality = minCentrality; }
  void SetMaxCentrality(int maxCentrality = 90) { fMaxCentrality = maxCentrality; }
  void SetMinPt(double minPt = 0.5) { fMinPt = minPt; }
  void SetMaxPt(double maxPt = 4.5) { fMaxPt = maxPt; }
  void SetRadiusOverflowCut(double cut = 10000.) { fRadiusOverflowCut = cut; }
  void SetRadiusV0OverflowCut(double cut = 10000.) { fRadiusV0OverflowCut = cut; }
  void SetDCAV0piToPVOverflowCut(double cut = 25.) { fDCAV0piToPVOverflowCut = cut; }
  void SetDCAV0prToPVOverflowCut(double cut = 12.7) { fDCAV0prToPVOverflowCut = cut; }
  void SetDCABachToPVOverflowCut(double cut = 2.5) { fDCABachToPVOverflowCut = cut; }
  void SetDCAV0toPVOverflowCut(double cut = 10.) { fDCAV0toPVOverflowCut = cut; }
  void SetBdtOutCut(double cut = 0.9) { fBdtOutCut = cut; }
  void SetNFeatures(int f = 9) { fNFeatures = f; }

  void SetFilterBit(double bit = BIT(4)) { fFilterBit = bit; }
  void SetTPCclsKaonCut(double cut = 70u) { fCutTPCclsKaon = cut; }
  void SetMaxChi2Cut(double cut = 4.) { fCutMaxChi2 = cut; }
  void SetMaxITSChi2Cut(double cut = 36.) { fCutMaxITSChi2 = cut; }
  void SetDCACut(double cutTight = .05, double cutMid = .1, double cutLoose = .5) { fCutDCA[0] = cutTight; fCutDCA[1] = cutMid; fCutDCA[2] = cutLoose; }
  void SetTPCclsCut(int cutTight = 80, int cutMid = 70, int cutLoose = 60) { fCutTPCcls[0] = cutTight; fCutTPCcls[1] = cutMid; fCutTPCcls[2] = cutLoose; }
  void SetChi2Cut(double cutTight = 2., double cutMid = 2.5, double cutLoose = 3.) { fCutChi2[0] = cutTight; fCutChi2[1] = cutMid; fCutChi2[2] = cutLoose; }
  void SetMaxPtKaon(double cut = 1.5) { fMaxPtKaon = cut; }
  void SetPtTofCut(double cut = 0.5) { fPtTofCut = cut; }
  void SetCutPtITSpid(double pt = 0.5) { fCutPtITSpid = pt; }
  void SetUseITSpid(bool toggle = true) { fUseITSpid = toggle; }

  void SetBDTPath(const char *path = "") { fBDTPath = path; }
  void SetPtBinsBDT(int nBins, double *ptBins) { fPtBinsBDT.Set(nBins+1,ptBins); }

  void SetCustomPidPath(const char* path = "") { fCustomPidPath = path; };
  void SetUseCustomPid(const bool toggle = true) { fUseCustomPid = toggle; };

  // Setters for Kaon track cuts

private:
  AliAnalysisTaskKaonXiCorrelation (const AliAnalysisTaskKaonXiCorrelation &source);
  AliAnalysisTaskKaonXiCorrelation &operator=(const AliAnalysisTaskKaonXiCorrelation &source);

  void PostAllData();

  TList*          fList = nullptr;             //!<! List of the output histograms
  TTree*          fTree = nullptr;             //!<! Tree for Xis

  MiniCollision fRecCollision;
  std::vector<MiniXi> fRecCascades;
  std::vector<MiniKaon> fRecKaons;
  std::vector<MiniXiMC> fGenCascades;
  std::vector<MiniKaonMC> fGenKaons;
  MiniXi* fXi = nullptr;
  MiniKaon* fKaon = nullptr;
  MiniXiMC fGenXi;
  MiniKaonMC fGenKaon;

  AliPIDResponse* fPID = nullptr;              //!<! ALICE PID framework
  std::vector<AliExternalBDT*> fBDT;           //!<! BDTs
  bool fMC;
  bool fOnlyTrueCandidates = true;  ///< Save only true Xi in MC
  bool fFillCascades = true;
  bool fFillKaons = true;
  bool fApplyBdtToMC = false;

  int fEstimator = 0;

  // configurable cuts
  float fCutRadius = 1.2;
  float fCutRadiusV0 = 3.0;
  float fCutDCABachToPV = 0.1;
  float fCutDCAV0toPV = 0.1;
  float fCutDCAV0piToPV = 0.2;
  float fCutDCAV0prToPV = 0.2;
  float fCutDCAV0tracks = 1.2;
  float fCutDCABachToV0 = 1.0;
  float fCutCosPA = 0.95;
  float fCutCosPAV0 = 0.95;
  float fCutBachBarCosPA = 0.99995;
  float fCutV0MassWindow = 0.005;
  float fCutY = 0.5;
  float fCutNsigmaTPC = 4.0;
  float fCutCt = 4;
  float fCutCtV0 = 30;
  float fCutCompetingMass = 0.008;
  float fCascMassWindow = 0.006;
  int fCutTPCclu = 70;
  float fCutTPCrows = 80.;
  float fCutRowsOvF = 0.8;
  double fCascLeastCRows;
  double fCascLeastCRowsOvF;
  double fV0LeastCRows;
  double fV0LeastCRowsOvF;
  double fMinCentrality = 0;
  double fMaxCentrality = 90;
  double fMinPt = 0.5;
  double fMaxPt = 4.5;
  double fRadiusOverflowCut = 51;
  double fRadiusV0OverflowCut = 100;
  double fDCAV0piToPVOverflowCut = 25;
  double fDCAV0prToPVOverflowCut = 12.7;
  double fDCABachToPVOverflowCut = 12.7;
  double fDCAV0toPVOverflowCut = 10.1;
  double fBdtOutCut = 0.9;
  int fNFeatures = 9;

  int fFilterBit = BIT(4);
  int fCutITSrecPoints = 2;
  int fCutSPDrecPoints = 1;
  int fCutSDDSSDrecPoints = 2;
  int fCutTPCclsKaon = 70;
  float fCutMaxChi2 = 2.5;
  float fCutMaxITSChi2 = 36.;
  float fCutDCA[3] = {0.05, 0.1, 0.5};
  int fCutTPCcls[3] = {80, 70, 60};
  float fCutChi2[3] = {1.5, 2., 2.5};
  double fCutKaonNsigmaTPC = 5.;
  double fCutKaonNsigmaTOF = 5.;
  double fCutKaonNsigmaITS = 5.;
  double fMaxPtKaon = 1.5;
  double fPtTofCut = 0.5;
  double fCutPtITSpid = 0.5;
  bool fUseITSpid = true;

  bool fUseOnTheFly = false;
  
  TArrayD fPtBinsBDT;
  std::string fBDTPath = "";
  std::string fCustomPidPath = "";
  TH3F* fCustomPidCalib[3] = {nullptr, nullptr, nullptr};
  bool fUseCustomPid = false;

  float Eta2y(float pt, float m, float eta) const;
  int WhichBDT(double pt);
  int GetITScls(AliAODTrack *track, int &nSPD, int &nSDD, int &nSSD);
  bool HasTOF(AliAODTrack *track);
  bool HasTwoXiFromSameDaughters(std::vector<XiDaughter> daughters);
  double GetCustomNsigma(AliAODTrack *t, double cent, int det = 0);

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskKaonXiCorrelation, 1);
  /// \endcond
};


#endif /* defined(__AliAnalysisTaskKaonXiCorrelation__) */
