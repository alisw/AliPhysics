/// \class AliAnalysisTaskAlphaPiAODtrk

#ifndef __AliAnalysisTaskAlphaPiAODtrk__
#define __AliAnalysisTaskAlphaPiAODtrk__

#include <Rtypes.h>
#include <TString.h>

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliPID.h"

class THistManager;
class AliPIDResponse;
class TH2F;
class TList;
class TTree;

struct StructHypertrk {
  Double32_t pt;
  Double32_t Rapidity;
  Double32_t m;
  Double32_t ct;
  Double32_t V0CosPA;
  Double32_t V0radius;
  Double32_t Lrec;
  Double32_t fZ;
  Double32_t TPCnSigmaPi;
  Double32_t TPCnSigmaalpha;
  Double32_t TOFnSigmaalpha;
  Double32_t TPCmomalpha;
  Double32_t TPCsignalalpha;
  Double32_t alphaProngPvDCAXY;
  Double32_t PiProngPvDCAXY;
  Double32_t alphaProngPvDCA;
  Double32_t PiProngPvDCA;
  Double32_t ProngsDCA;
  unsigned char NpidClustersPion;
  unsigned char NpidClustersalpha;
  unsigned char NitsClustersalpha;
  unsigned char centrality;
  unsigned char trigger;
  bool Matter;
};

struct StructHypertrkMC : public StructHypertrk {
  float ptMC;
  float etaMC;
  float ctMC;
  float yMC;
  int pdg;
  bool isReconstructed;
  bool isDuplicated = false;
};

class AliAnalysisTaskAlphaPiAODtrk : public AliAnalysisTaskSE {
public:
  enum kReducedTrigger {
    kINT7 = BIT(0),
    kCentral = BIT(1),
    kSemiCentral = BIT(2),
    kPositiveB = BIT(3),
    kHighMultV0 = BIT(4)
  };
  AliAnalysisTaskAlphaPiAODtrk(bool isMC = false, TString taskname = "HyperAOD");
  static AliAnalysisTaskAlphaPiAODtrk *AddTask(bool isMC, TString tskname,
                                            TString suffix);
  virtual ~AliAnalysisTaskAlphaPiAODtrk();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *) {}

  AliEventCuts fEventCut; ///<

  void SetCustomBetheBloch(float resolution, const float bethe[5]);
  double customNsigma(double mom, double sig);
  void SaveOnlyTrueCandidates(bool toggle = true) {
    fOnlyTrueCandidates = toggle;
  }
  void UseOnTheFly(bool toggle = true) { fUseOnTheFly = toggle; }
  void UseCustomPID(bool toggle = true) { fUseCustomPID = toggle; }
  void SetMassRange(float min, float max) {
    fMassRange[0] = min;
    fMassRange[1] = max;
  }
  void SetFilterbitTrackCut(Double_t lParameter) { fFilterBit = lParameter; }
  void SetNucleusPID(AliPID::EParticleType pdg) { fNucleusPID = pdg; }
  void SetPIDrange(float min, float max) {
    fPIDrange[0] = min;
    fPIDrange[1] = max;
  }

private:
  AliAnalysisTaskAlphaPiAODtrk(const AliAnalysisTaskAlphaPiAODtrk &source);
  AliAnalysisTaskAlphaPiAODtrk &operator=(const AliAnalysisTaskAlphaPiAODtrk &source);

  void PostAllData();

  TTree *fTree = nullptr; //!<! Tree for Hyper

  StructHypertrk *fRecHyper = nullptr; //!<! Transient fRecHyper
  StructHypertrkMC fGenHyper;
  AliPIDResponse *fPID = nullptr; //!<! ALICE PID framework
  bool fMC;
  THistManager *fHistos;           //!
  bool fOnlyTrueCandidates = true; ///< Save only true Hyperhydrogens
  bool fUseOnTheFly = false;
  bool fUseCustomPID = false; //!

  float fCustomBethe[5] = {
      1.70184, 28.4426, 3.21871e-12, 2.06952,
      2.77971}; /// default values are from AliAnalysisTaskHelium3PiAOD.cxx
  float fCustomResolution = 0.06; /// default values are for LHC18qr
  double fMassRange[2] = {3.7, 4.1};
  UInt_t fFilterBit = 16; // Bit(4) 16: Loose StandardITSTPC2011 cut.
  AliPID::EParticleType fNucleusPID = AliPID::kAlpha;
  float fPIDrange[2] = {-5.f, 5.f};

  float Eta2y(float pt, float m, float eta) const;

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskAlphaPiAODtrk, 8);
  // 2: Use THistManager class, add QA histograms for TPC PID of alpha.
  // 3: Use default PID selection and add option for the customise.
  // 4: Add track loop option
  // 5: Add track cut with filterbit
  // 6: Add additional functionality for broad PID cut.
  // 7: Add TOF PID
  // 8: Update custom PID nsigma.
  /// \endcond
};

#endif /* defined(__AliAnalysisTaskAlphaPiAODtrk__) */
