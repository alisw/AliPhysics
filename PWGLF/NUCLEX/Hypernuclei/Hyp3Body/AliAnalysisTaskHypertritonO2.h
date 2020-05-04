#ifndef AliAnalysisTaskHypertritonO2_H
#define AliAnalysisTaskHypertritonO2_H

#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "Math/Vector4D.h"
#include "DCAFitterN.h"

#include <TString.h>

#include <list>
#include <map>
#include <string>
#include <vector>

class TH1D;
class TH2D;
class TList;
class TTree;
class TFile;
class TSpline3;

class AliPIDResponse;
class AliESDtrack;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> LVector_t;

struct RHyperTritonO2 {
  float fCent = -1.;
  float pt = -999.f;
  float phi = -999.f;
  float pz = -999.f;
  float ct = -1.f;
  float r = -1.f;
  float cosPA = -2.f;
  float m = -1;
  Double32_t dca_de = -1.0; //[0.0,8.0,8]
  Double32_t dca_pr = -1.0; //[0.0,8.0,8]
  Double32_t dca_pi = -1.0; //[0.0,8.0,8]
  Double32_t tpcNsig_de = -4.0; //[-4.0,4.0,8]
  Double32_t tpcNsig_pr = -4.0; //[-4.0,4.0,8]
  Double32_t tpcNsig_pi = -4.0; //[-4.0,4.0,8]
  Double32_t tofNsig_de = -4.0; //[-4.0,4.0,8]
  Double32_t tofNsig_pr = -4.0; //[-4.0,4.0,8]
  Double32_t tofNsig_pi = -4.0; //[-4.0,4.0,8]
  Double32_t dca_de_pr = -4.0; //[0.0,8.0,8]
  Double32_t dca_de_pi = -4.0; //[0.0,8.0,8]
  Double32_t dca_pr_pi = -4.0; //[0.0,8.0,8]
  Double32_t dca_de_sv = -4.0; //[0.0,8.0,8]
  Double32_t dca_pr_sv = -4.0; //[0.0,8.0,8]
  Double32_t dca_pi_sv = -4.0; //[0.0,8.0,8]
  Double32_t chi2 = -1.f;      //[0.0,16.,16]
  UChar_t tpcClus_de = 0u;
  UChar_t tpcClus_pr = 0u;
  UChar_t tpcClus_pi = 0u;
  UChar_t candidates = 0u;
  UChar_t fTrigger = 0u;
  bool hasTOF_de = false;
  bool hasTOF_pr = false;
  bool hasTOF_pi = false;
  bool positive = false;
};

struct SHyperTritonO2 : public RHyperTritonO2 {
  SHyperTritonO2() : RHyperTritonO2{} {}
  SHyperTritonO2(const RHyperTritonO2& other) : RHyperTritonO2{other} {}
  float gPt = -999.f;
  float gPhi = -999.f;
  float gPz = -999.f;
  float gCt = -1.f;
  float gT = -1.f;
  bool  gPositive = false;
  bool  gReconstructed = false;
};

class AliAnalysisTaskHypertritonO2 : public AliAnalysisTaskSE {

public:
  enum kReducedTrigger { kINT7 = BIT(0), kCentral = BIT(1), kSemiCentral = BIT(2), kPositiveB = BIT(3) };

  AliAnalysisTaskHypertritonO2(bool mc = false, std::string name = "HyperTriton3O2");
  virtual ~AliAnalysisTaskHypertritonO2();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  static AliAnalysisTaskHypertritonO2 *AddTask(bool isMC = false, TString suffix = "");

  void SetDownscaling(bool down) { fDownscaling = down; }

  void SetDownscalingFactorByEvent(float fraction) { fDownscalingFactorByEvent = fraction; }
  void SetDownscalingFactorByCandidate(float fraction) { fDownscalingFactorByCandidate = fraction; }

  void SetEnableEventMixing(bool enableEM) { fEnableEventMixing = enableEM; }
  void SetEventMixingPoolDepth(int maxDepth) { fEventMixingPoolDepth = maxDepth; }

  AliEventCuts fEventCuts;                  /// Event cuts class
  AliESDtrackCuts fTrackCuts = *AliESDtrackCuts::GetStandardV0DaughterCuts(); /// Track cuts Object

  o2::vertexing::DCAFitter3 fVertexer;
  enum kProng { kDeuteron = 0, kProton = 1, kPion = 2 };
  bool  fSwapSign = false;
  float fMassWindow[2] = {2.94f, 3.06f};
  float  fRequireTOFpid[3] = {10.,10.,10.}; /// momentum after which the TOF matching is required
  int   fMinTPCpidClusters[3] = {70, 70, 50};
  float fMinTrackDCA[3] = {0., 0., 0.};
  float fTPCsigmas[3] = {3.5f, 3.5f, 3.5f};
  float fTOFsigmas[3] = {4.f, 4.f, 4.f};
  float fCandidateCtRange[2] = {0.f, 35.f};
  float fCandidatePtRange[2] = {1.f, 9.f};
  float fTrackPtRange[3][2] = {{0.f, 7.f},{0.f, 4.f},{0.f, 1.f}};
  float fMinCosPA = 0.9;
  bool  fUseAbsCosPAcut = true;
  bool  fOnlyTrueCandidates = false;
  std::string fCosPAsplineName = "PWGLF/NUCLEX/HypertritonAnalysis/Cuts/spline3.root";


private:
  void FillGenHypertriton(int id, bool reco, AliMCEvent* mcEv);

  int FindEventMixingCentBin(const float centrality);
  int FindEventMixingZBin(const float zVtx);
  void FillEventMixingPool(const float centrality, const float xVtx, std::vector<AliESDtrack *> tracks);
  std::vector<AliESDtrack *> GetEventMixingTracks(const float centrality, const float zvtx);

  AliInputEventHandler* fInputHandler = nullptr; //!
  AliPIDResponse*       fPIDResponse = nullptr;  //!

  TList *fListHist = nullptr;    //! List of Cascade histograms
  TTree *fTreeHyp3 = nullptr;    //! Output Tree, V0s

  TFile *fCosPAsplineFile = nullptr;  //! Pointer to the spline file
  TSpline3 *fCosPAspline = nullptr;   //! Pointer to the cosPA cut spline

  bool fMC = false;
  bool fDownscaling = false;
  bool fEnableEventMixing = false;

  /// Control histograms to monitor the filtering
  TH2D *fHistNSigmaDeu = nullptr;    //! # sigma TPC for the deuteron
  TH2D *fHistNSigmaP = nullptr;      //! # sigma TPC proton for the positive prong
  TH2D *fHistNSigmaPi = nullptr;     //! # sigma TPC pion for the negative prong
  TH2D *fHistInvMass = nullptr;      //! # Invariant mass histogram

  float fDownscalingFactorByEvent = 1.;        // fraction of the events saved in the tree
  float fDownscalingFactorByCandidate = 1.;    // fraction of the candidates saved in the tree

  std::list<AliESDtrack> fEventMixingPool[10][10];    /// container for the ESD used fot event mixing
  int fEventMixingPoolDepth = 0;                      /// max depth of the event mixing pool

  SHyperTritonO2   fGenHyp;
  float            fGenRecDeutMom;
  float            fGenRecProtMom;
  float            fGenRecPiMom;
  RHyperTritonO2   fRecHyp;

  AliAnalysisTaskHypertritonO2(const AliAnalysisTaskHypertritonO2 &);               // not implemented
  AliAnalysisTaskHypertritonO2 &operator=(const AliAnalysisTaskHypertritonO2 &);    // not implemented

  ClassDef(AliAnalysisTaskHypertritonO2, 1);
};

#endif