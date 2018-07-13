#ifndef AliAnalysisTaskStrangenessLifetimes_H
#define AliAnalysisTaskStrangenessLifetimes_H

class TH1D;
class TH2D;
class TList;
class TTree;

#include <string>
#include <vector>

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "MiniV0.h"

class AliPIDResponse;

class AliAnalysisTaskStrangenessLifetimes : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskStrangenessLifetimes(std::string name = "TaskStrangenessLifetimes");
  virtual ~AliAnalysisTaskStrangenessLifetimes();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);

  // Task Configuration: trigger selection
  void SetUseLightVertexers(bool lUseLightVertexers = true) {
    fUseLightVertexer = lUseLightVertexers;
  }
  void SetDoV0Refit(bool lDoV0Refit = true) { fDoV0Refit = lDoV0Refit; }

  // Setters for the V0 Vertexer Parameters
  void SetV0VertexerMaxChisquare(double lParameter) {
    fV0VertexerSels[0] = lParameter;
  }
  void SetV0VertexerDCAFirstToPV(double lParameter) {
    fV0VertexerSels[1] = lParameter;
  }
  void SetV0VertexerDCASecondtoPV(double lParameter) {
    fV0VertexerSels[2] = lParameter;
  }
  void SetV0VertexerDCAV0Daughters(double lParameter) {
    fV0VertexerSels[3] = lParameter;
  }
  void SetV0VertexerCosinePA(double lParameter) {
    fV0VertexerSels[4] = lParameter;
  }
  void SetV0VertexerMinRadius(double lParameter) {
    fV0VertexerSels[5] = lParameter;
  }
  void SetV0VertexerMaxRadius(double lParameter) {
    fV0VertexerSels[6] = lParameter;
  }

  void SetMinPt(float lMinPt) { fMinPtToSave = lMinPt; }
  void SetMaxPt(float lMaxPt) { fMaxPtToSave = lMaxPt; }
  void SetLambdaWindowParameters(double* fMeanPars, double* fSigmaPars) {
    for (Int_t ipar = 0; ipar < 5; ipar++)
      fLambdaMassMean[ipar] = fMeanPars[ipar];
    for (Int_t ipar = 0; ipar < 4; ipar++)
      fLambdaMassSigma[ipar] = fSigmaPars[ipar];
  }
  void SetLambdaWindowParametersStandard() {
    fLambdaMassMean[0] = 1.15768e+00;
    fLambdaMassMean[1] = -4.15945e-02;
    fLambdaMassMean[2] = -7.14294e-04;
    fLambdaMassMean[3] = -1.62793e-02;
    fLambdaMassMean[4] = -7.84067e+00;
    fLambdaMassSigma[0] = 1.30345e-03;
    fLambdaMassSigma[1] = 2.89679e-04;
    fLambdaMassSigma[2] = 1.52661e-03;
    fLambdaMassSigma[3] = -2.58251e+00;
  }

  void SetMaxTPCsigmas(float pi, float proton) {
    fMaxTPCpionSigma = pi;
    fMaxTPCprotonSigma = proton;
  }
  // Functions for analysis Bookkeepinp
  // 1- Configure standard vertexing
  void SetupStandardVertexing();
  void SetupLooseVertexing();

  AliEventCuts fEventCuts;  /// Event cuts class

 private:
  TList* fListHist;  //! List of Cascade histograms
  TTree* fTreeV0;    //! Output Tree, V0s

  AliPIDResponse* fPIDResponse;  //! PID response object

  bool fDoV0Refit;
  bool fUseLightVertexer;

  /// Control histograms to monitor the filtering
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
  TH1D* fHistNsigmaPosPion;         //! # sigma TPC pion for the positive prong
  TH1D* fHistNsigmaPosProton;       //! # sigma TPC proton for the positive prong
  TH1D* fHistNsigmaNegPion;         //! # sigma TPC pion for the negative prong
  TH1D* fHistNsigmaNegProton;       //! # sigma TPC proton for the negative prong
  TH1D* fHistEtaPos;                //! Pseudorapidity of the positive prong
  TH1D* fHistEtaNeg;                //! Pseudorapidity of the negative prong
  TH2D* fHistArmenteros;            //! Pseudorapidity of the negative prong

  double fV0VertexerSels[7];  // Array to store the 7 values for the different
                              // selections V0 related

  double fLambdaMassMean[5];  // Array to store the lambda mass mean
                              // parametrization
  //[0]+[1]*TMath::Exp([2]*x)+[3]*TMath::Exp([4]*x)

  double fLambdaMassSigma[4];  // Array to store the lambda mass sigma
                               // parametrization
  //[0]+[1]*x+[2]*TMath::Exp([3]*x)

  float fMinPtToSave;  // minimum pt
  float fMaxPtToSave;  // maximum pt
  float fMaxTPCpionSigma;
  float fMaxTPCprotonSigma;

  std::vector<Lifetimes::MiniV0> fV0vector;
  float fMultiplicity;

  AliAnalysisTaskStrangenessLifetimes(
      const AliAnalysisTaskStrangenessLifetimes&);  // not implemented
  AliAnalysisTaskStrangenessLifetimes& operator=(
      const AliAnalysisTaskStrangenessLifetimes&);  // not implemented

  ClassDef(AliAnalysisTaskStrangenessLifetimes, 1);
};

#endif
