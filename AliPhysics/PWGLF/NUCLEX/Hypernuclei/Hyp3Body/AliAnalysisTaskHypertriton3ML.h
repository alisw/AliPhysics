#ifndef AliAnalysisTaskHypertriton3ML_H
#define AliAnalysisTaskHypertriton3ML_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliVertexerHyperTriton3Body.h"
#include "Math/Vector4D.h"

#include <TObjString.h>
#include <TString.h>

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

  Double32_t fDCAxyDeu; // [0.0,5.12,9]
  Double32_t fDCAzDeu;  // [0.0,5.12,9]
  Double32_t fDCAxyP;   // [0.0,5.12,9]
  Double32_t fDCAzP;    // [0.0,5.12,9]
  Double32_t fDCAxyPi;  // [0.0,5.12,9]
  Double32_t fDCAzPi;   // [0.0,5.12,9]

  unsigned char fNClusterTPCDeu;
  unsigned char fNClusterTPCP;
  unsigned char fNClusterTPCPi;

  unsigned char fITSClusterMapDeu;
  unsigned char fITSClusterMapP;
  unsigned char fITSClusterMapPi;

  Double32_t fNSigmaTPCDeu; // [-5.12,5.12,8]
  Double32_t fNSigmaTPCP;   // [-5.12,5.12,8]
  Double32_t fNSigmaTPCPi;  // [-5.12,5.12,8]

  Double32_t fNSigmaTOFDeu; //  [-5.12,5.12,8]
  Double32_t fNSigmaTOFP;   //  [-5.12,5.12,8]
  Double32_t fNSigmaTOFPi;  //  [-5.12,5.12,8]

  bool fHasTOFDeu;
  bool fHasTOFP;
  bool fHasTOFPi;

  Double32_t fTrackChi2Deu; // [0.0,10,24,10]
  Double32_t fTrackChi2P;   // [0.0,10,24,10]
  Double32_t fTrackChi2Pi;  // [0.0,10,24,10]

  Double32_t fDecayVertexChi2NDF; // [0.0,102,4,10]

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
  int fRecoIndex; /// To connect with the reconstructed information

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

class AliAnalysisTaskHypertriton3ML : public AliAnalysisTaskSE {

public:
  AliAnalysisTaskHypertriton3ML(bool mc = false, std::string name = "HyperTriton2He3piML");
  virtual ~AliAnalysisTaskHypertriton3ML();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  static AliAnalysisTaskHypertriton3ML *AddTask(bool isMC = false, TString suffix = "");

  void SetDownscaling(bool down) { fDownscaling = down; }

  void SetDownscalingFactorByEvent(float fraction) { fDownscalingFactorByEvent = fraction; }
  void SetDownscalingFactorByCandidate(float fraction) { fDownscalingFactorByCandidate = fraction; }

  void SetOnlyTrueCandidates(bool trueCand) { fOnlyTrueCandidates = trueCand; }

  void SetMinCandidatePt(float lPtMin) { fMinCanidatePtToSave = lPtMin; }
  void SetMaxCandidatePt(float lPtMax) { fMaxCanidatePtToSave = lPtMax; }

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

  void SetMinITSnCluster(unsigned char nClsMin) { fMinITSNcluster = nClsMin; }
  void SetMinTPCcluster(unsigned char minCls) { fMinTPCNcluster = minCls; }

  void SetMaxPtPion(float piPtMax) { fMaxPtPion = piPtMax; }

  void SetVertexerToleranceGuessCompatibility(int nGuessCompatibility) {
    fVertexerToleranceGuessCompatibility = nGuessCompatibility;
  }
  void SetVertexerMaxDistanceInit(float vtxMaxDistanceInit) { fVertexerMaxDistanceInit = vtxMaxDistanceInit; }

  void SetMinCosPointingAngle(float minCosPA) { fMinCosPA = minCosPA; }

  AliEventCuts fEventCuts; /// Event cuts class

  AliVertexerHyperTriton3Body fVertexer; //

private:
  TList *fListHist; //! List of Cascade histograms
  TTree *fTreeHyp3; //! Output Tree, V0s

  AliInputEventHandler *fInputHandler; //!
  AliPIDResponse *fPIDResponse;        //! PID response object

  bool fMC;
  bool fOnlyTrueCandidates;
  bool fDownscaling;

  /// Control histograms to monitor the filtering
  TH2D *fHistNSigmaDeu; //! # sigma TPC for the deuteron
  TH2D *fHistNSigmaP;   //! # sigma TPC proton for the positive prong
  TH2D *fHistNSigmaPi;  //! # sigma TPC pion for the negative prong
  TH2D *fHistInvMass;   //! # Invariant mass histogram

  float fDownscalingFactorByEvent;     // fraction of the events saved in the tree
  float fDownscalingFactorByCandidate; // fraction of the candidates saved in the tree

  float fMinCanidatePtToSave; // min candidate pt to save
  float fMaxCanidatePtToSave; // max candidate pt to save

  unsigned char fMinITSNcluster;
  unsigned char fMinTPCNcluster;

  float fMaxNSigmaTPCDeu; // nSigma TPC limit for deuteron
  float fMaxNSigmaTPCP;   // nSigma TPC limit for proton
  float fMaxNSigmaTPCPi;  // nSigma TPC limit for pion

  float fMaxNSigmaTOFDeu; // nSigma TOF limit for deuteron
  float fMaxNSigmaTOFP;   // nSigma TOF limit for proton
  float fMaxNSigmaTOFPi;  // nSigma TOF limit for pion

  int fVertexerToleranceGuessCompatibility; // minimum number of compatible tracks with the guess of the decay vertex
  float fVertexerMaxDistanceInit;           // max distance from the guess vertex for the compatible tracks

  float fMinCosPA; // minimum cos(poninting angle) accepted

  float fMinDCA2PrimaryVtxDeu;
  float fMinDCA2PrimaryVtxP;
  float fMinDCA2PrimaryVtxPi;

  float fMaxPtPion;

  std::vector<SHypertriton3> fSHypertriton; //!
  std::vector<RHypertriton3> fRHypertriton; //!
  REvent fREvent;                           //!

  std::vector<AliESDtrack *> fDeuVector;
  std::vector<AliESDtrack *> fPVector;
  std::vector<AliESDtrack *> fPiVector;

  AliAnalysisTaskHypertriton3ML(const AliAnalysisTaskHypertriton3ML &);            // not implemented
  AliAnalysisTaskHypertriton3ML &operator=(const AliAnalysisTaskHypertriton3ML &); // not implemented

  ClassDef(AliAnalysisTaskHypertriton3ML, 1);
};

#endif