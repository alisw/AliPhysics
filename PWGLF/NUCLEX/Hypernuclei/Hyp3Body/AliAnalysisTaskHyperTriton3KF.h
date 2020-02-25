#ifndef AliAnalysisTaskHyperTriton3KF_H
#define AliAnalysisTaskHyperTriton3KF_H

#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "Math/Vector4D.h"

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

struct REvent3KF
{
  float fX;
  float fY;
  float fZ;
  float fCent;
  unsigned char fTrigger;
};

struct SHyperTriton3KF {
  float px = -999.f;
  float py = -999.f;
  float pz = -999.f;
  float l = -1.f;
  int reco = -1;
};


struct RHyperTriton3KF {
  float px = -999.f;
  float py = -999.f;
  float pz = -999.f;
  float l = -1.f;
  float r = -1.f;
  float chi2_deuprot = -1.f;
  float chi2_3prongs = -1.f;
  float chi2_topology = -1.f;
  float cosPA = -1.f;
  float m = -1;
};

class AliAnalysisTaskHyperTriton3KF : public AliAnalysisTaskSE {

public:
  enum kReducedTrigger { kINT7 = BIT(0), kCentral = BIT(1), kSemiCentral = BIT(2), kPositiveB = BIT(3) };

  AliAnalysisTaskHyperTriton3KF(bool mc = false, std::string name = "HyperTriton3KF");
  virtual ~AliAnalysisTaskHyperTriton3KF();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  static AliAnalysisTaskHyperTriton3KF *AddTask(bool isMC = false, TString suffix = "");

  void SetDownscaling(bool down) { fDownscaling = down; }

  void SetDownscalingFactorByEvent(float fraction) { fDownscalingFactorByEvent = fraction; }
  void SetDownscalingFactorByCandidate(float fraction) { fDownscalingFactorByCandidate = fraction; }

  void SetEnableEventMixing(bool enableEM) { fEnableEventMixing = enableEM; }
  void SetEventMixingPoolDepth(int maxDepth) { fEventMixingPoolDepth = maxDepth; }

  AliEventCuts fEventCuts;                  /// Event cuts class
  AliESDtrackCuts fTrackCuts = *AliESDtrackCuts::GetStandardV0DaughterCuts(); /// Track cuts Object

  enum kProng { kDeuteron = 0, kProton = 1, kPion = 2 };
  bool  fRequireTOFpid[3] = {false,false,false};
  int   fMinTPCpidClusters[3] = {100, 100, 70};
  float fMinTrackDCA[3] = {0., 0., 0.};
  float fTPCsigmas[3] = {3.5f, 3.5f, 3.5f};
  float fTOFsigmas[3] = {5.f, 5.f, 5.f};
  float fCandidateCtRange[2] = {0.f, 40.f};
  float fCandidatePtRange[2] = {0.f, 10.f};
  float fTrackPtRange[3][2] = {{0.f, 7.f},{0.f, 4.f},{0.f, 1.f}};
  float fMinCosPA = 0.99;
  float fMaxKFchi2[3] = {400.,400.,400.};
  bool  fOnlyTrueCandidates = false;


private:

  int FindEventMixingCentBin(const float centrality);
  int FindEventMixingZBin(const float zVtx);
  void FillEventMixingPool(const float centrality, const float xVtx, std::vector<AliESDtrack *> tracks);
  std::vector<AliESDtrack *> GetEventMixingTracks(const float centrality, const float zvtx);

  AliInputEventHandler* fInputHandler = nullptr; //!
  AliPIDResponse*       fPIDResponse = nullptr;  //!

  TList *fListHist = nullptr;    //! List of Cascade histograms
  TTree *fTreeHyp3 = nullptr;    //! Output Tree, V0s

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

  REvent3KF                    fREvent;
  std::vector<SHyperTriton3KF> fGenHyp;
  std::vector<RHyperTriton3KF> fRecHyp;

  AliAnalysisTaskHyperTriton3KF(const AliAnalysisTaskHyperTriton3KF &);               // not implemented
  AliAnalysisTaskHyperTriton3KF &operator=(const AliAnalysisTaskHyperTriton3KF &);    // not implemented

  ClassDef(AliAnalysisTaskHyperTriton3KF, 1);
};

#endif