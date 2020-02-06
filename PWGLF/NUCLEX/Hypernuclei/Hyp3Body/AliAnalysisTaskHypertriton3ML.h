#ifndef AliAnalysisTaskHypertriton3ML_H
#define AliAnalysisTaskHypertriton3ML_H

#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "AliMLResponse.h"
#include "AliVertexerHyperTriton3Body.h"
#include "Math/Vector4D.h"

#include <TObjString.h>
#include <TString.h>

#include <list>
#include <map>
#include <string>
#include <vector>

class TH1D;
class TH2D;
class TList;
class TTree;

class AliPIDResponse;
class AliESDtrack;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> LVector_t;

struct RHypertriton3 {
  float fDecayVtxX;
  float fDecayVtxY;
  float fDecayVtxZ;

  float fPxDeu;
  float fPyDeu;
  float fPzDeu;
  float fPxP;
  float fPyP;
  float fPzP;
  float fPxPi;
  float fPyPi;
  float fPzPi;

  float fPosXDeu;
  float fPosYDeu;
  float fPosZDeu;
  float fPosXP;
  float fPosYP;
  float fPosZP;
  float fPosXPi;
  float fPosYPi;
  float fPosZPi;

  Double32_t fDCAxyDeu;    // [0.0,5.12,9]
  Double32_t fDCAzDeu;     // [0.0,5.12,9]
  Double32_t fDCAxyP;      // [0.0,5.12,9]
  Double32_t fDCAzP;       // [0.0,5.12,9]
  Double32_t fDCAxyPi;     // [0.0,5.12,9]
  Double32_t fDCAzPi;      // [0.0,5.12,9]

  unsigned char fNClusterTPCDeu;
  unsigned char fNClusterTPCP;
  unsigned char fNClusterTPCPi;

  unsigned char fITSClusterMapDeu;
  unsigned char fITSClusterMapP;
  unsigned char fITSClusterMapPi;

  Double32_t fNSigmaTPCDeu;    // [-5.12,5.12,8]
  Double32_t fNSigmaTPCP;      // [-5.12,5.12,8]
  Double32_t fNSigmaTPCPi;     // [-5.12,5.12,8]

  Double32_t fNSigmaTOFDeu;    //  [-5.12,5.12,8]
  Double32_t fNSigmaTOFP;      //  [-5.12,5.12,8]
  Double32_t fNSigmaTOFPi;     //  [-5.12,5.12,8]

  bool fHasTOFDeu;
  bool fHasTOFP;
  bool fHasTOFPi;

  Double32_t fTrackChi2Deu;    // [0.0,10,24,10]
  Double32_t fTrackChi2P;      // [0.0,10,24,10]
  Double32_t fTrackChi2Pi;     // [0.0,10,24,10]

  Double32_t fDecayVertexChi2NDF;    // [0.0,102,4,10]

  bool fIsMatter;
};

struct REvent {
  float fX;
  float fY;
  float fZ;

  float fCent;

  unsigned char fTrigger;
};

struct SHypertriton3 {
  int fRecoIndex;    /// To connect with the reconstructed information

  long fPdgCode;

  float fDecayVtxX;
  float fDecayVtxY;
  float fDecayVtxZ;

  float fPxDeu;
  float fPyDeu;
  float fPzDeu;
  float fPxP;
  float fPyP;
  float fPzP;
  float fPxPi;
  float fPyPi;
  float fPzPi;
};

struct MLSelected {
  float score;
  float fInvMass;
  float fCt;
  float fCentrality;
  float fCandPt;
};

class AliAnalysisTaskHypertriton3ML : public AliAnalysisTaskSE {

public:
  enum kReducedTrigger { kINT7 = BIT(0), kCentral = BIT(1), kSemiCentral = BIT(2), kPositiveB = BIT(3) };

  AliAnalysisTaskHypertriton3ML(bool mc = false, std::string name = "HyperTriton3ML");
  virtual ~AliAnalysisTaskHypertriton3ML();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  static AliAnalysisTaskHypertriton3ML *AddTask(bool isMC = false, TString suffix = "");

  void SetDownscaling(bool down) { fDownscaling = down; }

  void SetDownscalingFactorByEvent(float fraction) { fDownscalingFactorByEvent = fraction; }
  void SetDownscalingFactorByCandidate(float fraction) { fDownscalingFactorByCandidate = fraction; }

  void SetApplyML(bool applyML) { fApplyML = applyML; }
  void SetMLResponseConfigfilePath(std::string configfilepath) { fMLResponseConfigfileRemotePath = configfilepath; }

  void SetEnableEventMixing(bool enableEM) { fEnableEventMixing = enableEM; }
  void SetEventMixingPoolDepth(int maxDepth) { fEventMixingPoolDepth = maxDepth; }

  void SetOnlyTrueCandidates(bool trueCand) { fOnlyTrueCandidates = trueCand; }

  void SetMinCandidatePt(float lPtMin) { fMinCanidatePtToSave = lPtMin; }
  void SetMaxCandidatePt(float lPtMax) { fMaxCanidatePtToSave = lPtMax; }

  void SetMinCandidateCt(float lCtMin) { fMinCanidateCtToSave = lCtMin; }
  void SetMaxCandidateCt(float lCtMax) { fMaxCanidateCtToSave = lCtMax; }

  void SetMaxEventCentrality(float lCentMax) { fMaxEventCentrality = lCentMax; }

  void SetMaxTPCsigmas(float nSigmaDeu, float nSigmaP, float nSigmaPi) {
    fMaxNSigmaTPCDeu = nSigmaDeu;
    fMaxNSigmaTPCP   = nSigmaP;
    fMaxNSigmaTPCPi  = nSigmaPi;
  }

  void SetMaxTOFsigmas(float nSigmaDeu, float nSigmaP, float nSigmaPi) {
    fMaxNSigmaTOFDeu = nSigmaDeu;
    fMaxNSigmaTOFP   = nSigmaP;
    fMaxNSigmaTOFPi  = nSigmaPi;
  }

  void SetMinDCA2PrimaryVtx(float dcaDeu, float dcaP, float dcaPi) {
    fMinDCA2PrimaryVtxDeu = dcaDeu;
    fMinDCA2PrimaryVtxP   = dcaP;
    fMinDCA2PrimaryVtxPi  = dcaPi;
  }

  void SetMinTPCcluster(unsigned char minCls) { fMinTPCNcluster = minCls; }

  void SetMaxPtDeu(float deuPtMax) { fMaxPtDeu = deuPtMax; }
  void SetMaxPtP(float pPtMax) { fMaxPtP = pPtMax; }
  void SetMaxPtPi(float piPtMax) { fMaxPtPi = piPtMax; }

  void SetVertexerToleranceGuessCompatibility(int nGuessCompatibility) {
    fVertexerToleranceGuessCompatibility = nGuessCompatibility;
  }
  void SetVertexerMaxDistanceInit(float vtxMaxDistanceInit) { fVertexerMaxDistanceInit = vtxMaxDistanceInit; }

  void SetMinCosPointingAngle(float minCosPA) { fMinCosPA = minCosPA; }

  AliEventCuts fEventCuts;                  /// Event cuts class
  AliVertexerHyperTriton3Body fVertexer;    /// custom 3-body decay vertexer
  AliMLResponse *fMLResponse;               /// object for the ML application
  AliESDtrackCuts fTrackCuts;               /// Track cuts Object

private:
  std::map<std::string, double> FeaturesMap(const RHypertriton3 &hypCand, const REvent &rEv);

  int FindEventMixingCentBin(const float centrality);
  int FindEventMixingZBin(const float zVtx);
  void FillEventMixingPool(const float centrality, const float xVtx, std::vector<AliESDtrack *> tracks);
  std::vector<AliESDtrack *> GetEventMixingTracks(const float centrality, const float zvtx);

  TList *fListHist;    //! List of Cascade histograms
  TTree *fTreeHyp3;    //! Output Tree, V0s

  AliInputEventHandler *fInputHandler;    //!
  AliPIDResponse *fPIDResponse;           //! PID response object

  bool fMC;
  bool fOnlyTrueCandidates;
  bool fDownscaling;
  bool fApplyML;
  bool fEnableEventMixing;

  /// Control histograms to monitor the filtering
  TH2D *fHistNSigmaDeu;    //! # sigma TPC for the deuteron
  TH2D *fHistNSigmaP;      //! # sigma TPC proton for the positive prong
  TH2D *fHistNSigmaPi;     //! # sigma TPC pion for the negative prong
  TH2D *fHistInvMass;      //! # Invariant mass histogram

  float fDownscalingFactorByEvent;        // fraction of the events saved in the tree
  float fDownscalingFactorByCandidate;    // fraction of the candidates saved in the tree

  float fMinCanidatePtToSave;    // min candidate pt to save
  float fMaxCanidatePtToSave;    // max candidate pt to save

  float fMinCanidateCtToSave;    // min candidate ct to save
  float fMaxCanidateCtToSave;    // max candidate ct to save

  float fMaxEventCentrality;    // max event centrality to save

  unsigned char fMinTPCNcluster;

  float fMaxNSigmaTPCDeu;    // nSigma TPC limit for deuteron
  float fMaxNSigmaTPCP;      // nSigma TPC limit for proton
  float fMaxNSigmaTPCPi;     // nSigma TPC limit for pion

  float fMaxNSigmaTOFDeu;    // nSigma TOF limit for deuteron
  float fMaxNSigmaTOFP;      // nSigma TOF limit for proton
  float fMaxNSigmaTOFPi;     // nSigma TOF limit for pion

  int fVertexerToleranceGuessCompatibility;    // minimum number of compatible tracks with the guess of the decay vertex
  float fVertexerMaxDistanceInit;              // max distance from the guess vertex for the compatible tracks

  float fMinCosPA;    // minimum cos(poninting angle) accepted

  float fMinDCA2PrimaryVtxDeu;
  float fMinDCA2PrimaryVtxP;
  float fMinDCA2PrimaryVtxPi;

  float fMaxPtDeu;
  float fMaxPtP;
  float fMaxPtPi;

  std::vector<SHypertriton3> fSHypertriton;    //!
  std::vector<RHypertriton3> fRHypertriton;    //!
  REvent fREvent;                              //!

  std::vector<MLSelected> fMLSelected;    //!

  std::vector<AliESDtrack *> fDeuVector;
  std::vector<AliESDtrack *> fPVector;
  std::vector<AliESDtrack *> fPiVector;

  std::string fMLResponseConfigfileRemotePath;    /// path for the remote ML config file
  std::string fMLResponseConfigfileLocalPath;     /// path for the local ML config file

  std::list<AliESDtrack> fEventMixingPool[10][10];    /// container for the ESD used fot event mixing
  int fEventMixingPoolDepth;                          /// max depth of the event mixing pool

  AliAnalysisTaskHypertriton3ML(const AliAnalysisTaskHypertriton3ML &);               // not implemented
  AliAnalysisTaskHypertriton3ML &operator=(const AliAnalysisTaskHypertriton3ML &);    // not implemented

  ClassDef(AliAnalysisTaskHypertriton3ML, 1);
};

#endif