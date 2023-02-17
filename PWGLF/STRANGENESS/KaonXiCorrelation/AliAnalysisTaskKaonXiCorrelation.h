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
class TList;
class TTree;

struct MiniXi {
  Double32_t fPt; //[-12.7,12.8,8]
  Double32_t fEta; //[-1.27,1.28,8]
  Double32_t fMass; //[1.29,1.35375,8]
  Double32_t fBdtOut; //[0.5,1.,16]
};

struct MiniXiMC : public MiniXi {
  Double32_t fPtMC; //[-12.7,12.8,8]
  Double32_t fEtaMC; //[-1.27,1.28,8]
  bool fIsReconstructed;
  unsigned char fFlag;
};

struct MiniKaon {
  Double32_t fPt; //[-12.7,12.8,8]
  Double32_t fEta; //[-1.27,1.28,8]
  Double32_t fNsigmaTPC; //[-6.35,6.4,8]
  Double32_t fNsigmaTOF; //[-6.35,6.4,8]
};

struct MiniKaonMC : public MiniKaon {
  Double32_t fPtMC; //[-12.7,12.8,8]
  Double32_t fEtaMC; //[-1.27,1.28,8]
  bool fIsReconstructed;
  unsigned char fFlag;
};

struct MiniCollision {
  Double32_t fZ; //[-12.7,12.8,8]
  unsigned char fCent;
  unsigned char fTrigger;
};

class AliAnalysisTaskKaonXiCorrelation : public AliAnalysisTaskSE {
public:
  enum kStatusFlag {
    kPrimary = BIT(0),
    kSecondaryFromWD = BIT(1),
    kSecondaryFromMaterial = BIT(2)
  };

  enum kReducedTrigger
  {
    kINT7 = BIT(0),
    kCentral = BIT(1),
    kSemiCentral = BIT(2),
    kPositiveB = BIT(3)
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
  void SetCompetingMassCut(float cut = 0.) { fCutCompetingMass = cut; }
  void SetTPCcluCut(int cut = 70) { fCutTPCclu = cut; }
  void SetSaveOnlyTrueCandidates(bool cut = true) { fOnlyTrueCandidates = cut; }
  void SetTPCRowsCut(float cut = 80.) { fCutTPCrows = cut; }
  void SetTPCRowOvFCut(float cut = 0.8) { fCutRowsOvF = cut; }
  void UseOnTheFly(bool toggle = true) { fUseOnTheFly = toggle; }
  void SetMinCentrality(int minCentrality = 0) { fMinCentrality = minCentrality; }
  void SetMaxCentrality(int maxCentrality = 90) { fMaxCentrality = maxCentrality; }
  void SetMinPt(double minPt = 0.5) { fMinPt = minPt; }
  void SetMaxPt(double maxPt = 4.5) { fMaxPt = maxPt; }
  void SetRadiusOverflowCut(double cut = 25.) { fRadiusOverflowCut = cut; }
  void SetRadiusV0OverflowCut(double cut = 25.) { fRadiusV0OverflowCut = cut; }
  void SetDCAV0piToPVOverflowCut(double cut = 2.5) { fDCAV0piToPVOverflowCut = cut; }
  void SetDCAV0prToPVOverflowCut(double cut = 2.5) { fDCAV0prToPVOverflowCut = cut; }
  void SetDCABachToPVOverflowCut(double cut = 2.5) { fDCABachToPVOverflowCut = cut; }
  void SetDCAV0toPVOverflowCut(double cut = 2.5) { fDCAV0toPVOverflowCut = cut; }
  void SetBdtOutCut(double cut = 0.9) { fBdtOutCut = cut; }

  void SetTPCsignalCut(double cut = 70u) { fCutTPCsignal = cut; }
  void SetMaxChi2Cut(double cut = 4.) { fCutMaxChi2 = cut; }
  void SetMaxITSChi2Cut(double cut = 36.) { fCutMaxITSChi2 = cut; }
  void SetDCAzCut(double cut = 1.) { fCutDCAz = cut; }
  void SetDCAxyCut(double cut = 0.1) { fCutDCAxy = cut; }
  void SetMaxPtKaon(double cut = 1.5) { fMaxPtKaon = cut; }

  void SetBDTPath(const char *path = "") { fBDTPath = path; }
  void SetPtBinsBDT(int nBins, double *ptBins) { fPtBinsBDT.Set(nBins+1,ptBins); }

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
  float fCutCompetingMass = 0.;
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
  double fRadiusOverflowCut = 25;
  double fRadiusV0OverflowCut = 25;
  double fDCAV0piToPVOverflowCut = 2.5;
  double fDCAV0prToPVOverflowCut = 2.5;
  double fDCABachToPVOverflowCut = 2.5;
  double fDCAV0toPVOverflowCut = 2.5;
  double fBdtOutCut = 0.9;

  int fCutTPCrecPoints = 70;
  int fCutITSrecPoints = 2;
  int fCutSPDrecPoints = 1;
  int fCutTPCsignal = 70;
  float fCutMaxChi2 = 4.;
  float fCutMaxITSChi2 = 36.;
  float fCutTPCfoundFraction = 0.8;
  float fCutMinEnergyLoss = 0.;
  float fCutDCAz = 1.;
  float fCutDCAxy = 0.1;
  double fCutKaonNsigmaTPC = 5.;
  double fCutKaonNsigmaTOF = 5.;
  double fMaxPtKaon = 1.5;

  bool fUseOnTheFly = false;
  
  TArrayD fPtBinsBDT;
  std::string fBDTPath = "";

  float Eta2y(float pt, float m, float eta) const;
  int WhichBDT(double pt);
  int GetITScls(AliAODTrack *track, int &nSPD);
  bool HasTOF(AliAODTrack *track);

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskKaonXiCorrelation, 1);
  /// \endcond
};


#endif /* defined(__AliAnalysisTaskKaonXiCorrelation__) */
