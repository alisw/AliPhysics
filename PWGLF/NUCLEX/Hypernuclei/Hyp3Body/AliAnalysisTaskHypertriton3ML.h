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

  float fDCAxyDeu;
  float fDCAzDeu;
  float fDCAxyP;
  float fDCAzP;
  float fDCAxyPi;
  float fDCAzPi;

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
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);

  static AliAnalysisTaskHypertriton3ML* AddTask(bool isMC = false, TString suffix = "");

  void SetDownscaling(bool down) { fDownscaling = down; }

  void SetDownscalingFactorByEvent(float fraction) { fDownscalingFactorByEvent = fraction; }
  void SetDownscalingFactorByCandidate(float fraction) { fDownscalingFactorByCandidate = fraction; }

  void SetApplyML(bool applyML) { fApplyML = applyML; }
  void SetMLResponseConfigfilePath(std::string configfilepath) { fMLResponseConfigfilePath = configfilepath; }

  void SetEnableEventMixing(bool enableEM) { fEnableEventMixing = enableEM; }
  void SetEventMixingPoolDepth(int maxDepth) { fEventMixingPoolDepth = maxDepth; }

  void SetOnlyTrueCandidates(bool trueCand) { fOnlyTrueCandidates = trueCand; }

  void SetMinTPCclusterAll(unsigned char minCls) { fMinTPCNclusterAll = minCls; }
  void SetMinTPCcluster(unsigned char minClsDeu, unsigned char minClsP, unsigned char minClsPi) {
    fMinTPCNclusterDeuteron = minClsDeu;
    fMinTPCNclusterProton   = minClsP;
    fMinTPCNclusterPion     = minClsPi;
  }

  void SetMinDCA2PrimaryVtx(float dcaDeu, float dcaP, float dcaPi) {
    fMinDCA2PrimaryVtxDeuteron = dcaDeu;
    fMinDCA2PrimaryVtxProton   = dcaP;
    fMinDCA2PrimaryVtxPion     = dcaPi;
  }

  void SetPtRangeDeuteron(float deuPtMin, float deuPtMax) {
    fMinPtDeuteron = deuPtMin;
    fMaxPtDeuteron = deuPtMax;
  }
  void SetPtRangeProton(float pPtMin, float pPtMax) {
    fMinPtProton = pPtMin;
    fMaxPtProton = pPtMax;
  }
  void SetPtRangePion(float piPtMin, float piPtMax) {
    fMinPtPion = piPtMin;
    fMaxPtPion = piPtMax;
  }

  void SetTPCNSigmaRangeDeuteron(float nSigmaMinDeu, float nSigmaMaxDeu) {
    fMinNSigmaTPCDeuteron = nSigmaMinDeu;
    fMaxNSigmaTPCDeuteron = nSigmaMaxDeu;
  }

  void SetTPCNSigmaRangeProton(float nSigmaMinP, float nSigmaMaxP) {
    fMinNSigmaTPCProton = nSigmaMinP;
    fMaxNSigmaTPCProton = nSigmaMaxP;
  }

  void SetTPCNSigmaRangePion(float nSigmaMinPi, float nSigmaMaxPi) {
    fMinNSigmaTPCProton = nSigmaMinPi;
    fMaxNSigmaTPCProton = nSigmaMaxPi;
  }

  void SetTOFRequired(bool deuTOF, bool pTOF, bool piTOF) {
    fRequireDeuteronTOFpid = deuTOF;
    fRequireProtonTOFpid   = pTOF;
    fRequirePionTOFpid     = piTOF;
  }

  void SetMaxTOFsigmas(float nSigmaDeu, float nSigmaP, float nSigmaPi) {
    fMaxNSigmaTOFDeuteron = nSigmaDeu;
    fMaxNSigmaTOFProton   = nSigmaP;
    fMaxNSigmaTOFPion     = nSigmaPi;
  }

  void SetVertexerToleranceGuessCompatibility(int nGuessCompatibility) {
    fVertexerToleranceGuessCompatibility = nGuessCompatibility;
  }
  void SetVertexerMaxDistanceInit(float vtxMaxDistanceInit) { fVertexerMaxDistanceInit = vtxMaxDistanceInit; }

  void SetMinCandidatePt(float lPtMin) { fMinCanidatePtToSave = lPtMin; }
  void SetMaxCandidatePt(float lPtMax) { fMaxCanidatePtToSave = lPtMax; }

  void SetMinCandidateCt(float lCtMin) { fMinCanidateCtToSave = lCtMin; }
  void SetMaxCandidateCt(float lCtMax) { fMaxCanidateCtToSave = lCtMax; }

  void SetMinCosPointingAngle(float minCosPA) { fMinCosPA = minCosPA; }

  AliEventCuts fEventCuts;                  /// Event cuts class
  AliVertexerHyperTriton3Body fVertexer;    /// custom 3-body decay vertexer
  AliMLResponse* fMLResponse;               /// object for the ML application

private:
  bool AcceptDeuteronCandidate(const AliESDtrack* deuTrack, float dca, float nSigmaTPC);
  bool AcceptProtonCandidate(const AliESDtrack* pTrack, float dca, float nSigmaTPC);
  bool AcceptPionCandidate(const AliESDtrack* piTrack, float dca, float nSigmaTPC);

  std::map<std::string, double> FeaturesMap(const RHypertriton3& hypCand, const REvent& rEv);

  int FindEventMixingCentBin(const float centrality);
  int FindEventMixingZBin(const float zVtx);
  void FillEventMixingPool(const float centrality, const float xVtx, std::vector<AliESDtrack*> tracks);
  std::vector<AliESDtrack*> GetEventMixingTracks(const float centrality, const float zvtx);

  TList* fListHist;    //! List of Cascade histograms
  TTree* fTreeHyp3;    //! Output Tree, V0s

  AliInputEventHandler* fInputHandler;    //!
  AliPIDResponse* fPIDResponse;           //! PID response object

  bool fMC;
  bool fOnlyTrueCandidates;
  bool fDownscaling;
  bool fApplyML;
  bool fEnableEventMixing;

  bool fRequireDeuteronTOFpid;
  bool fRequireProtonTOFpid;
  bool fRequirePionTOFpid;

  float fDownscalingFactorByEvent;        // fraction of the events saved in the tree
  float fDownscalingFactorByCandidate;    // fraction of the candidates saved in the tree

  unsigned char fMinTPCNclusterAll;    // minimum number of TPC cluster for all daughters

  unsigned char fMinTPCNclusterDeuteron;    // minimum number of TPC cluster for Deuteron
  unsigned char fMinTPCNclusterProton;      // minimum number of TPC cluster for Proton
  unsigned char fMinTPCNclusterPion;        // minimum number of TPC cluster for Pion

  float fMinDCA2PrimaryVtxDeuteron;    // minimum DCA to the primary vertex for deuteron
  float fMinDCA2PrimaryVtxProton;      // minimum DCA to the primary vertex for proton
  float fMinDCA2PrimaryVtxPion;        // minimum DCA to the primary vertex for pion

  float fMinPtDeuteron;    // deuteron pT lower limit
  float fMaxPtDeuteron;    // deuteron pT upper limit
  float fMinPtProton;      // proton pT lower limit
  float fMaxPtProton;      // proton pT upper limit
  float fMinPtPion;        // pion pT lower limit
  float fMaxPtPion;        // pion pT upper limit

  float fMinNSigmaTPCDeuteron;    // nSigma TPC lower limit for deuteron
  float fMaxNSigmaTPCDeuteron;    // nSigma TPC upper limit for deuteron
  float fMinNSigmaTPCProton;      // nSigma TPC lower limit for proton
  float fMaxNSigmaTPCProton;      // nSigma TPC upper limit for proton
  float fMinNSigmaTPCPion;        // nSigma TPC lower limit for pion
  float fMaxNSigmaTPCPion;        // nSigma TPC upper limit for pion

  float fMaxNSigmaTOFDeuteron;    // nSigma TOF limit for deuteron
  float fMaxNSigmaTOFProton;      // nSigma TOF limit for proton
  float fMaxNSigmaTOFPion;        // nSigma TOF limit for pion

  int fVertexerToleranceGuessCompatibility;    // minimum number of compatible tracks with the guess of the decay vertex
  float fVertexerMaxDistanceInit;              // max distance from the guess vertex for the compatible tracks

  float fMinCanidatePtToSave;    // min candidate pt to save
  float fMaxCanidatePtToSave;    // max candidate pt to save

  float fMinCanidateCtToSave;    // min candidate ct to save
  float fMaxCanidateCtToSave;    // max candidate ct to save

  float fMinCosPA;    // minimum cos(poninting angle) accepted

  std::vector<SHypertriton3> fSHypertriton;    //!
  std::vector<RHypertriton3> fRHypertriton;    //!
  REvent fREvent;                              //!

  std::vector<MLSelected> fMLSelected;    //!

  std::vector<AliESDtrack*> fDeuteronVector;
  std::vector<AliESDtrack*> fProtonVector;
  std::vector<AliESDtrack*> fPionVector;

  std::string fMLResponseConfigfilePath;    /// path for the remote ML config file

  std::list<AliESDtrack> fEventMixingPool[10][10];    /// container for the ESD used fot event mixing
  int fEventMixingPoolDepth;                          /// max depth of the event mixing pool

  AliAnalysisTaskHypertriton3ML(const AliAnalysisTaskHypertriton3ML&);               // not implemented
  AliAnalysisTaskHypertriton3ML& operator=(const AliAnalysisTaskHypertriton3ML&);    // not implemented

  ClassDef(AliAnalysisTaskHypertriton3ML, 1);
};

#endif