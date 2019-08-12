#ifndef AliAnalysisTaskStrangenessLifetimes_H
#define AliAnalysisTaskStrangenessLifetimes_H

class TH1D;
class TGraph;
class TH2D;
class TList;
class TTree;

#include <string>
#include <vector>

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "Math/Vector4D.h"
#include "MCparticle.h"
#include "MiniV0.h"
#include "MiniEvent.h"
#include "HyperTriton2Body.h"


class AliPIDResponse;
class AliESDtrack;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> LVector_t;

class AliAnalysisTaskStrangenessLifetimes : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskStrangenessLifetimes(bool mc = false, std::string name = "TaskStrangenessLifetimes",float downscale=1,bool Hypertriton=true, bool V0s=true);
  virtual ~AliAnalysisTaskStrangenessLifetimes();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);

  // Task Configuration: trigger selection
  void SetUseOnTheFlyV0s(bool useThem = true) {
    fUseOnTheFly = useThem;
  }
  void SetDoV0Refit(bool lDoV0Refit = true) { fDoV0Refit = lDoV0Refit; }

  void SetMinPt(float lMinPt) { fMinPtToSave = lMinPt; }
  void SetMaxPt(float lMaxPt) { fMaxPtToSave = lMaxPt; }

  void SetCustomBetheBloch(float resolution, const float bethe[5]);

  void SetMaxTPCsigmas(float pi, float proton, float he3) {
    fMaxTPCpionSigma = pi;
    fMaxTPCprotonSigma = proton;
    fMaxTPChe3Sigma = he3;
  }
  // Functions for analysis Bookkeepinp
  // 1- Configure standard vertexing
  void SetupStandardVertexing();
  void SetupLooseVertexing();

  static LVector_t GetV0LorentzVector(int pdg, AliESDtrack* nTrack, AliESDtrack* pTrack, double alpha);

  AliEventCuts fEventCuts;  /// Event cuts class
  bool fHypertriton;
  bool fV0s;

 private:
  TList* fListHist;  //! List of Cascade histograms
  TTree* fTreeV0;    //! Output Tree, V0s

  AliPIDResponse* fPIDResponse;  //! PID response object

  bool fDoV0Refit;
  bool fMC;
  float fDownscale;
  bool fUseOnTheFly;

  bool  fUseCustomBethe;
  float fCustomBethe[5];
  float fCustomResolution;

  /// Control histograms to monitor the filtering
  TH1D* fHistMCct[2];               //! MC ct
  TH1D* fHistMCctPrimary[2];        //! MC ct only for primary particles
  TH1D* fHistMCctSecondaryFromMaterial[2]; //! MC ct for secondaries from material
  TH1D* fHistV0radius;              //! V0 decay vertex radius
  TH1D* fHistV0pt;                  //! V0 transverse momentum
  TH1D* fHistV0eta;                 //! V0 pseudorapidity
  TH2D* fHistInvMassK0s;            //! Invariant mass for K0s
  TH2D* fHistInvMassLambda;         //! Invariant mass for (anti-)Lambda
  TH1D* fHistDistOverTotMom;        //! L/p
  TH1D* fHistV0CosPA;               //! V0 cosine of pointing angle
  TH1D* fHistChi2V0;                //! V0 fit chi2
  TH1D* fHistDcaNeg2PrimaryVertex;  //! DCA of the negative prong to the PV
  TH1D* fHistDcaPos2PrimaryVertex;  //! DCA of the positive prong to the PV
  TH1D* fHistDcaV0daughters;        //! DCA between the two prongs
  TH1D* fHistV0armAlpha;            //! Armenteros alpha
  TH1D* fHistV0armPt;               //! Armenteros pt
  TH1D* fHistLeastNxedRows;         //! Min number of xed roads
  TH1D* fHistLeastXedOverFindable;  //! Min number of xed roads/findable clusters
  TH1D* fHistMaxChi2PerCluster;     //! Max chi2 per cluster in TPC
  TH1D* fHistNsigmaPosPion;         //! # sigma TPC pion for the positive prongedRo
  TH1D* fHistNsigmaPosProton;       //! # sigma TPC proton for the positive prong
  TH1D* fHistNsigmaNegPion;         //! # sigma TPC pion for the negative prong
  TH1D* fHistNsigmaNegProton;       //! # sigma TPC proton for the negative prong
  TH1D* fHistEtaPos;                //! Pseudorapidity of the positive prong
  TH1D* fHistEtaNeg;                //! Pseudorapidity of the negative prong
  TH2D* fHistArmenteros;            //! Pseudorapidity of the negative prong
  TH1D* fHistNsigmaPosHe;           //!
  TH1D* fHistNsigmaNegHe;           //!
  TH2D* fHistdEdxVsPt;              //!
  TH2D* fHistCtAnalysis;            //!  
  TH1D* fHistNhyp;                  //!
  float fMinPtToSave;  // minimum pt
  float fMaxPtToSave;  // maximum pt
  float fMaxTPCpionSigma;
  float fMaxTPCprotonSigma;
  float fMaxTPChe3Sigma;

  std::vector<Lifetimes::MiniV0 > fV0vector;
  std::vector<Lifetimes::MCparticle> fMCvector;
  std::vector<Lifetimes::HyperTriton2Body> fV0Hyvector;
  Lifetimes::MiniEvent fMiniEvent;

  AliAnalysisTaskStrangenessLifetimes(
      const AliAnalysisTaskStrangenessLifetimes&);  // not implemented
  AliAnalysisTaskStrangenessLifetimes& operator=(
      const AliAnalysisTaskStrangenessLifetimes&);  // not implemented

  ClassDef(AliAnalysisTaskStrangenessLifetimes, 7);
};

#endif