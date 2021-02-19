#ifndef AliAnalysisTaskHypertriton3_H
#define AliAnalysisTaskHypertriton3_H

#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "Math/Vector4D.h"
#include "DCAFitterN.h"
#include "Hypertriton3structures.h"
#include "AliVertexerHyperTriton2Body.h"

#include <TString.h>
#include <AliMCEvent.h>

#include <list>
#include <map>
#include <string>
#include <vector>
#include <utility>

class TH1D;
class TH2D;
class TList;
class TTree;
class TFile;
class TSpline3;

class AliPIDResponse;
class AliESDtrack;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> LVector_t;

constexpr float kHyperTritonMass{2.99131};
/// helper functions
template <typename T>
double Sq(T a) { return a * a; }
template <typename F>
double Hypot(F a, F b, F c) { return std::sqrt(Sq(a) + Sq(b) + Sq(c)); }
template <typename F>
double Hypot(F a, F b, F c, F d) { return std::sqrt(Sq(a) + Sq(b) + Sq(c) + Sq(d)); }

struct EventMixingTrack {
  EventMixingTrack(AliESDtrack *t, float tpc, float tof, int us) : track{*t}, nSigmaTPC{tpc}, nSigmaTOF{tof}, used{us} {} 
  AliESDtrack track = AliESDtrack();
  float nSigmaTPC = -10.;
  float nSigmaTOF = -10.;
  int used = 0;
};

class AliAnalysisTaskHypertriton3 : public AliAnalysisTaskSE {

public:
  enum kReducedTrigger { kINT7 = BIT(0), kCentral = BIT(1), kSemiCentral = BIT(2), kPositiveB = BIT(3), kHighMultV0 = BIT(4) };

  AliAnalysisTaskHypertriton3(bool mc = false, std::string name = "HyperTriton3O2");
  virtual ~AliAnalysisTaskHypertriton3();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  static AliAnalysisTaskHypertriton3 *AddTask(bool isMC = false, TString suffix = "");

  void SetDownscaling(bool down) { fDownscaling = down; }

  void SetDownscalingFactorByEvent(float fraction) { fDownscalingFactorByEvent = fraction; }
  void SetDownscalingFactorByCandidate(float fraction) { fDownscalingFactorByCandidate = fraction; }

  void SetEnableEventMixing(bool enableEM) { fEnableEventMixing = enableEM; }
  void SetEventMixingPoolDepth(int maxDepth) { fEventMixingPoolDepth = maxDepth; }
  void SetEventMixingPoolMaxReuse(int maxDepth) { fEventMixingPoolMaxReuse = maxDepth; }

  int CheckPionCharge(std::vector<AliESDtrack *> tracks[2][2], AliESDv0 v0);

  AliEventCuts fEventCuts;                  /// Event cuts class
  AliESDtrackCuts fTrackCuts = *AliESDtrackCuts::GetStandardV0DaughterCuts(); /// Track cuts Object
  int fCounter;
  o2::vertexing::DCAFitter3 fVertexer;
  o2::vertexing::DCAFitter2 fVertexerLambda;
  enum kProng { kDeuteron = 0, kProton = 1, kPion = 2 };
  bool  fSwapSign = false;
  int   fMixingTrack = 0;
  float fMassWindow[2] = {2.94f, 3.06f};
  float  fRequireTOFpid[3] = {10.,10.,10.}; /// momentum after which the TOF matching is required
  int   fMinTPCpidClusters[3] = {70, 70, 50};
  float fMinTrackDCA[3] = {0., 0., 0.};
  float fMaxTrack2TrackDCA[3] = {8.,8.,8.};
  float fMaxTrack2SVDCA[3] = {8.,8.,8.};
  float fTPCsigmas[3] = {3.f, 3.f, 3.f};
  float fTOFsigmas[3] = {4.f, 4.f, 4.f};
  float fCandidateCtRange[2] = {0.f, 35.f};
  float fCandidatePtRange[2] = {1.f, 9.f};
  float fTrackPtRange[3][2] = {{0.f, 7.f},{0.f, 4.f},{0.f, 1.f}};
  float fMinCosPA = 0.99;
  bool  fUseAbsCosPAcut = true;
  bool  fOnlyTrueCandidates = false;
  bool  fLambdaCheck = true;
  bool  fKF = false;
  bool  fUseDoubleV0s = false;
  bool  fUseCovarianceCut = false;
  float fMaxKFchi2[3] = {40000.,40000.,40000.};
  std::string fCosPAsplineName = "PWGLF/NUCLEX/HypertritonAnalysis/Cuts/spline3.root";
  AliVertexerHyperTriton2Body fV0Vertexer;


private:
  template<class T>
  void FillGenHypertriton(T* ptr, int id, bool reco, AliMCEvent* mcEv);
  int FindEventMixingCentBin(const float centrality);
  int FindEventMixingZBin(const float zVtx);
  void FillEventMixingPool(const float centrality, const float xVtx, const std::vector<EventMixingTrack> &tracks);
  std::vector<EventMixingTrack*> GetEventMixingTracks(const float centrality, const float zvtx);

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
  TH1D *fHistDecVertexRes = nullptr; //! # decay vertex resolution

  float fDownscalingFactorByEvent = 1.;        // fraction of the events saved in the tree
  float fDownscalingFactorByCandidate = 1.;    // fraction of the candidates saved in the tree

  std::list<EventMixingTrack> fEventMixingPool[10][10];        //! container for the ESD used fot event mixing
  unsigned int fEventMixingPoolDepth = 10;                     /// max depth of the event mixing pool
  int fEventMixingPoolMaxReuse = 2;

  SHyperTriton3KF*   fGenHypKF = nullptr;
  SHyperTriton3O2*   fGenHypO2 = nullptr;
  RHyperTriton*   fRecHyp = nullptr;

  AliAnalysisTaskHypertriton3(const AliAnalysisTaskHypertriton3 &);               // not implemented
  AliAnalysisTaskHypertriton3 &operator=(const AliAnalysisTaskHypertriton3 &);    // not implemented

  ClassDef(AliAnalysisTaskHypertriton3, 1);
};

template<class T>
void AliAnalysisTaskHypertriton3::FillGenHypertriton(T* ptr, int id, bool reco, AliMCEvent *mcEvent)
{

  double mcVtx[3];
  mcEvent->GetPrimaryVertex()->GetXYZ(mcVtx);
  AliVParticle *part = mcEvent->GetTrack(id);

  if (!part)
  {
    ::Warning("AliAnalysisTaskHypertriton3::FillGenHypertriton",
              "Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skipping.", id);
    return;
  }

  double decayVtx[4]{0.0, 0.0, 0.0, 0.0};

  for (int iD = part->GetDaughterFirst(); iD <= part->GetDaughterLast(); ++iD)
  {
    AliVParticle *daughter = mcEvent->GetTrack(iD);

    if (mcEvent->IsSecondaryFromWeakDecay(iD) && daughter && std::abs(daughter->PdgCode()) != 11)
    {
      decayVtx[0] = daughter->Xv();
      decayVtx[1] = daughter->Yv();
      decayVtx[2] = daughter->Zv();
      decayVtx[3] = daughter->Tv() - part->Tv();
      break;
    }
  }
  ptr->gCt = Hypot(mcVtx[0] - decayVtx[0], mcVtx[1] - decayVtx[1], mcVtx[2] - decayVtx[2]) * kHyperTritonMass / part->P();
  ptr->gPt = part->Pt();
  ptr->gPhi = std::atan2(part->Py(), part->Px());
  ptr->gPz = part->Pz();
  ptr->gT = decayVtx[3];
  ptr->gPositive = part->PdgCode() > 0;
  ptr->gReconstructed = reco;
}

#endif
