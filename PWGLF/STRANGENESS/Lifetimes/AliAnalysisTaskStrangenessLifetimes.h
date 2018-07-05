#ifndef AliAnalysisTaskStrangenessLifetimes_H
#define AliAnalysisTaskStrangenessLifetimes_H

class TList;
class TTree;
class AliESDtrackCuts;

#include <vector>
#include <string>

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "MiniV0.h"

class AliPIDResponse;

class AliAnalysisTaskStrangenessLifetimes : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskStrangenessLifetimes(std::string name);
  virtual ~AliAnalysisTaskStrangenessLifetimes();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

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
  void SetLambdaWindowParameters(double *fMeanPars, double *fSigmaPars) {
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
  // Functions for analysis Bookkeepinp
  // 1- Configure standard vertexing
  void SetupStandardVertexing();
  void SetupLooseVertexing();

  AliEventCuts fEventCuts;  /// Event cuts class

 private:
  TList *fListHist;  //! List of Cascade histograms
  TTree *fTreeV0;    //! Output Tree, V0s

  AliPIDResponse *fPIDResponse;  // PID response object
  AliESDtrackCuts
      *fESDtrackCuts;  // ESD track cuts used for primary track definition

  bool fDoV0Refit;
  bool fUseLightVertexer;

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

  std::vector<Lifetimes::MiniV0> fV0vector;
  float fMultiplicity;

  AliAnalysisTaskStrangenessLifetimes(
      const AliAnalysisTaskStrangenessLifetimes &);  // not implemented
  AliAnalysisTaskStrangenessLifetimes &operator=(
      const AliAnalysisTaskStrangenessLifetimes &);  // not implemented

  ClassDef(AliAnalysisTaskStrangenessLifetimes, 1);
};

#endif
