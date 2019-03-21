#ifndef AliAnalysisTaskHyperTriton2He3piML_H
#define AliAnalysisTaskHyperTriton2He3piML_H

class TH1D;
class TGraph;
class TH2D;
class TList;
class TTree;

#include <TString.h>
#include <string>
#include <vector>

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "Math/Vector4D.h"


class AliPIDResponse;
class AliESDtrack;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> LVector_t;

struct RHyperTritonHe3pi {
  float fDecayX;
  float fDecayY;
  float fDecayZ;
  float fPxHe3;
  float fPyHe3;
  float fPzHe3;
  float fPxPi;
  float fPyPi;
  float fPzPi;
  Double32_t fChi2V0;               //[0.0,10.24,8] V0 fit chi2
  Double32_t fDcaHe32PrimaryVertex; //[0.0,0.256,8] DCA of the negative prong to the PV
  Double32_t fDcaPi2PrimaryVertex; //[0.0,0.256,8]  DCA of the positive prong to the PV
  Double32_t fDcaV0daughters;       //[0.0,2.56,8] DCA between the two prongs
  Double32_t fLeastXedOverFindable; //[0.36,1.0,8] Min xed roads/findable clusters
  Double32_t fMaxChi2PerCluster;    //[0,6.4,8] Max chi2 per cluster in TPC
  Double32_t fTPCnSigmaHe3;         //[-8.0,8.0,6] number of sigmas TPC pion for the positive prong
  Double32_t fTPCnSigmaPi;          //[-8.0,8.0,6] number of sigmas TPC 3He for the positive prong
  unsigned char fLeastNxedRows;     // Min number of xed roads
  unsigned char fLeastNpidClusters; // Min number of PID clusters in TPC
  bool fITSrefitHe3;
  bool fITSrefitPi;
  bool fSPDanyHe3;
  bool fSPDanyPi;
  bool fCowboy;
  bool fMatter;
};

struct RCollision {
  float fX;
  float fY;
  float fZ;
  float fCent;
};

struct SHyperTritonHe3pi {
  int   fRecoIndex;  /// To connect with the reconstructed information
  int   fPdgCode;
  float fDecayX;
  float fDecayY;
  float fDecayZ;
  float fPxHe3;
  float fPyHe3;
  float fPzHe3;
  float fPxPi;
  float fPyPi;
  float fPzPi;
  bool  fFake;
};

class AliAnalysisTaskHyperTriton2He3piML : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskHyperTriton2He3piML(bool mc = false, std::string name = "HyperTriton2He3piML");
  virtual ~AliAnalysisTaskHyperTriton2He3piML();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);

  static AliAnalysisTaskHyperTriton2He3piML* AddTask(bool isMC = false, TString suffix = "");

  void SetUseOnTheFlyV0s(bool toogle = true) { fUseOnTheFly = toogle; }

  void SetMinPt(float lMinPt) { fMinPtToSave = lMinPt; }
  void SetMaxPt(float lMaxPt) { fMaxPtToSave = lMaxPt; }

  void SetMinPtHe3(float min) { fMinHe3pt = min; }

  void SetCustomBetheBloch(float resolution, const float bethe[5]);

  void SetMaxTPCsigmas(float pi, float he3) {
    fMaxTPCpiSigma = pi;
    fMaxTPChe3Sigma = he3;
  }

  AliEventCuts fEventCuts;  /// Event cuts class

 private:
  TList* fListHist;  //! List of Cascade histograms
  TTree* fTreeV0;    //! Output Tree, V0s

  AliPIDResponse* fPIDResponse;  //! PID response object

  bool fMC;
  bool fUseOnTheFly;

  bool  fUseCustomBethe;
  float fCustomBethe[5];
  float fCustomResolution;

  /// Control histograms to monitor the filtering
  TH2D* fHistNsigmaHe3;          //! # sigma TPC proton for the positive prong
  TH2D* fHistNsigmaPi;           //! # sigma TPC pion for the negative prong
  TH2D* fHistInvMass;            //! # Invariant mass histogram

  float fMinPtToSave;  // minimum pt
  float fMaxPtToSave;  // maximum pt
  float fMaxTPCpiSigma;
  float fMaxTPChe3Sigma;
  float fMinHe3pt;

  std::vector<SHyperTritonHe3pi> fSHyperTriton;
  std::vector<RHyperTritonHe3pi> fRHyperTriton;
  RCollision fRCollision;

  AliAnalysisTaskHyperTriton2He3piML(
      const AliAnalysisTaskHyperTriton2He3piML&);  // not implemented
  AliAnalysisTaskHyperTriton2He3piML& operator=(
      const AliAnalysisTaskHyperTriton2He3piML&);  // not implemented

  ClassDef(AliAnalysisTaskHyperTriton2He3piML, 1);
};

#endif