#ifndef AliAnalysisTaskHyperTriton2He3piML_H
#define AliAnalysisTaskHyperTriton2He3piML_H

class TH1D;
class TGraph;
class TH2D;
class TList;
class TTree;

#include <TObjString.h>
#include <TString.h>
#include <string>
#include <vector>
#include "AliVertexerHyperTriton2Body.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "Math/Vector4D.h"

class TSpline3;
class AliPIDResponse;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> LVector_t;

struct RHyperTritonHe3pi
{
  float fDecayX;
  float fDecayY;
  float fDecayZ;
  float fPxHe3;
  float fPyHe3;
  float fPzHe3;
  float fPxPi;
  float fPyPi;
  float fPzPi;
  Double32_t fTPCmomHe3;              //[0.0,10.24,8]
  Double32_t fTPCmomPi;               //[0.0,10.24,8]
  Double32_t fChi2V0;                 //[0.0,10.24,8] V0 fit chi2
  Double32_t fDcaHe32PrimaryVertexXY; //[0.0,5.12,10] DCA of the negative prong to the PV
  Double32_t fDcaPi2PrimaryVertexXY;  //[0.0,20.48,12]  DCA of the positive prong to the PV
  Double32_t fDcaHe32PrimaryVertex;   //[0.0,10.24,10] DCA of the negative prong to the PV
  Double32_t fDcaPi2PrimaryVertex;    //[0.0,40.96,12]  DCA of the positive prong to the PV
  Double32_t fDcaV0daughters;         //[0.0,2.56,8] DCA between the two prongs
  Double32_t fLeastXedOverFindable;   //[0.36,1.0,8] Min xed roads/findable clusters
  Double32_t fMaxChi2PerCluster;      //[0,6.4,8] Max chi2 per cluster in TPC
  Double32_t fTPCnSigmaHe3;           //[-8.0,8.0,7] number of sigmas TPC pion
  Double32_t fTPCnSigmaPi;            //[-8.0,8.0,7] number of sigmas TPC 3He
  Double32_t fTOFnSigmaHe3;           //[-8.0,8.0,7] number of sigmas TOF pion
  Double32_t fTOFnSigmaPi;            //[-8.0,8.0,7] number of sigmas TOF 3He
  Double32_t fTPCsignalHe3;           //[0.,2048.,12] signal of the He3 track in the TPC
  Double32_t fTPCsignalPi;            //[0.,2048.,12] signal of the pion track in the TPC
  unsigned char fNpidClustersHe3;     // Number of PID clusters in TPC He3
  unsigned char fNpidClustersPi;      // Number of PID clusters in TPC Pion
  unsigned char fITSclusHe3;
  unsigned char fITSclusPi;
  bool fITSrefitHe3;
  bool fITSrefitPi;
  bool fTOFmatchHe3;
  bool fTOFmatchPi;
  bool fCowboy;
  bool fMatter;
};

struct RCollision
{
  float fX;
  float fY;
  float fZ;
  float fCent;
  unsigned char fTrigger;
};

struct RTracklet
{
  float fTheta;
  float fPhi;
  Double32_t fDeltaTheta; //[-0.12,0.12,8]
  Double32_t fDeltaPhi;   //[-0.12,0.12,8]
  bool fSharedCluster;
};

struct SHyperTritonHe3pi
{
  int fRecoIndex;    /// To connect with the reconstructed information
  int fRecoTracklet; /// To connect with the reconstructed information of the tracklets
  long fPdgCode;
  float fDecayX;
  float fDecayY;
  float fDecayZ;
  float fPxHe3;
  float fPyHe3;
  float fPzHe3;
  float fPxPi;
  float fPyPi;
  float fPzPi;
  bool fFake;
  bool fNegativeLabels;
};

struct SGenericV0
{ /// For the other V0s that are reconstructed
  int fRecoIndex;
  long fPdgCode;
  float fDecayX;
  float fDecayY;
  float fDecayZ;
  float fPx;
  float fPy;
  float fPz;
};

struct SGenericTracklet
{ /// For the other V0s that are reconstructed
  int fRecoIndex;
  long fPdgCode;
  float fPx;
  float fPy;
  float fPz;
};

class AliAnalysisTaskHyperTriton2He3piML : public AliAnalysisTaskSE
{
public:
  enum kReducedTrigger
  {
    kINT7 = BIT(0),
    kCentral = BIT(1),
    kSemiCentral = BIT(2),
    kPositiveB = BIT(3)
  };

  AliAnalysisTaskHyperTriton2He3piML(bool mc = false, std::string name = "HyperTriton2He3piML");
  virtual ~AliAnalysisTaskHyperTriton2He3piML();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  static AliAnalysisTaskHyperTriton2He3piML *AddTask(bool isMC = false, TString suffix = "");

  void SetUseOnTheFlyV0s(bool toogle = true) { fUseOnTheFly = toogle; }
  void SetUseNanoAODs(bool toogle = true) { fUseNanoAODs = toogle; }

  void SetMinPt(float lMinPt) { fMinPtToSave = lMinPt; }
  void SetMaxPt(float lMaxPt) { fMaxPtToSave = lMaxPt; }

  void SetMinPtHe3(float min) { fMinHe3pt = min; }

  void SetCustomBetheBloch(float resolution, const float bethe[5]);
  double customNsigma(double mom, double sig);

  template <class T, class M>
  void FillHyperCandidate(T *v0, AliVEvent *event, AliMCEvent *mcEvent, M mcMap, double *pP,
                          double *nP, int lKeyPos, int lKeyNeg, RHyperTritonHe3pi v0part, int he3index);

  void SetMaxTPCsigmas(float pi, float he3)
  {
    fMaxTPCpiSigma = pi;
    fMaxTPChe3Sigma = he3;
  }

  void SetMinTPCcluster(unsigned char minCl)
  {
    fMinTPCclusters = minCl;
  }

  void SetMinPIDcluster(unsigned char minClu)
  {
    fMinPIDclusters = minClu;
  }

  void SetMaxDeltaTheta(float maxDeltaTheta) { fMaxDeltaTheta = maxDeltaTheta; }
  void SetMaxDeltaPhi(float maxDeltaPhi) { fMaxDeltaPhi = maxDeltaPhi; }
  void SetMinTrackletCosP(float minTrackletCosP) { fMinTrackletCosP = minTrackletCosP; }
  void EnableLikeSign(bool enableIt = true)
  {
    fEnableLikeSign = enableIt;
    fV0Vertexer.fLikeSign = enableIt;
  }

  void SetCVMFSPath(std::string path) { fCVMFSPath = path; }

  AliEventCuts fEventCuts; /// Event cuts class
  bool fFillGenericV0s;
  bool fFillGenericTracklets; /// To check what is the background
  bool fFillTracklet;
  bool fStoreAllEvents;
  bool fSaveFileNames;
  bool fPropagetToPV;
  AliVertexerHyperTriton2Body fV0Vertexer; //
private:
  TList *fListHist; //! List of Cascade histograms
  TTree *fTreeV0;   //! Output Tree, V0s

  AliInputEventHandler *fInputHandler; //!
  AliPIDResponse *fPIDResponse;        //! PID response object
  std::string fCVMFSPath;
  bool fMC;
  bool fUseOnTheFly;
  bool fUseNanoAODs;

  bool fUseCustomBethe;
  float fCustomBethe[5];
  float fCustomResolution;

  /// Control histograms to monitor the filtering
  TH2D *fHistNsigmaHe3;          //! # sigma TPC proton for the positive prong
  TH2D *fHistNsigmaPi;           //! # sigma TPC pion for the negative prong
  TH2D *fHistInvMass;            //! # Invariant mass histogram
  TH2D *fHistTPCdEdx[2];         //! # TPC dE/dx for V0s
  TH2D *fHistTrackletThetaPhi;   //! # tracklet theta vs phi
  TH2D *fHistTrackletDThetaDPhi; //! # tracklet delta_theta vs delta_phi
  TH1D *fHistTrackletCosP;       //! # tracklet-V0 cosine of pointing angle

  float fMinPtToSave; // minimum pt
  float fMaxPtToSave; // maximum pt
  float fMaxTPCpiSigma;
  float fMaxTPChe3Sigma;
  float fMinHe3pt;
  unsigned char fMinTPCclusters;
  unsigned char fMinPIDclusters;

  float fMaxDeltaPhi;
  float fMaxDeltaTheta;
  float fMinTrackletCosP;

  bool fEnableLikeSign;

  TTree *fFileNameTree;        //!
  TObjString fCurrentFileName; //!

  std::vector<SHyperTritonHe3pi> fSHyperTriton;     //!
  std::vector<SGenericV0> fSGenericV0;              //!
  std::vector<RHyperTritonHe3pi> fRHyperTriton;     //!
  std::vector<RTracklet> fRTracklets;               //!
  std::vector<SGenericTracklet> fSGenericTracklets; //!
  RCollision fRCollision;                           //!

  AliAnalysisTaskHyperTriton2He3piML(
      const AliAnalysisTaskHyperTriton2He3piML &); // not implemented
  AliAnalysisTaskHyperTriton2He3piML &operator=(
      const AliAnalysisTaskHyperTriton2He3piML &); // not implemented

  ClassDef(AliAnalysisTaskHyperTriton2He3piML, 7);
};

#endif