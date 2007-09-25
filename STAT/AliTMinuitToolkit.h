#ifndef AliTMINUITTOOLKIT_H
#define AliTMINUITTOOLKIT_H


#include <TNamed.h>
#include <TVirtualFitter.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFormula.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMath.h>
#include <TString.h>


class AliTMinuitToolkit: public TNamed{
  //
  //
  //
public:
  AliTMinuitToolkit();
  virtual ~AliTMinuitToolkit();
  
  void FitHistogram(TH1F * his);
  void FitHistogram(TH2F * his);
  void SetInitialParam(TVectorD *param) { fParam=param;};
  void SetParamLimits(TMatrixD *paramLimits) { fParamLimits=paramLimits;};
  //
  void SetFitFunction(TFormula * formula) {fFormula=formula;};
  void SetWeightFunction(TFormula * weightFunction) {fWeightFunction=weightFunction;};
  void SetWeightFunction(Char_t * name, Float_t param1, Float_t param2 = 0.1);
  void SetFitAlgorithm(Char_t * name) {fFitAlgorithm=name;};
  void SetMaxCalls(Int_t calls) {fMaxCalls = calls;};
  void SetTolerance(Double_t tol) {fPrecision = tol;};
  void SetPoints(TMatrixD * points) {fPoints = points;};
  void SetWeights(TVectorD * weights) {fWeights = weights;};
  static void Test();
  static void FitterFCN(int &npar, double  *dummy, double &fchisq, double *gin, int iflag);
  void Fit();
  void EnableRobust(Bool_t b) {fUseRobust = b;};
  void EnableRobust(Bool_t b, Double_t sigma) {fUseRobust = b; fExpectedSigma = sigma;}; 
  //
  TMatrixD * GetPoints() {return fPoints;};
  TVectorD * GetWeights() {return fWeights;};
  TFormula * GetFormula() {return fFormula;};
  TVectorD * GetParameters() {return fParam;};
  TFormula * GetWeightFunction() {return fWeightFunction;};
  TMatrixD * GetCovarianceMatrix() {return fCovar;};
  Bool_t     GetStatus() {return fUseRobust;};
  Double_t   GetExpectedSigma() {return fExpectedSigma;};
 
private:
  //
  AliTMinuitToolkit(const AliTMinuitToolkit&); // fake copy constr. -- suppress warnings
  AliTMinuitToolkit& operator=(const AliTMinuitToolkit&); // fake -- suppress warnings
  //
  TFormula        * fFormula;            // formula of the fitted function
  TFormula        * fWeightFunction;     // weight function, must be defined between 0 and 1
  Char_t          * fFitAlgorithm;       // fit algorithm for TMinuit: migrad, simplex, ...
  //
  //
  TMatrixD        * fPoints;             // fitted points
  TVectorD        * fWeights;            // weights of the points
  TVectorD        * fParam;              // parameter values
  TMatrixD        * fParamLimits;        // limits of the parameters (up, down)
  TMatrixD        * fCovar;              // covariance matrix
  Double_t        * fChi2;               // chi-square
  Int_t             fMaxCalls;           // maximum number of calls, default value depends on fit algorithm
  Double_t          fPrecision;          // normalized tolerance, default value depends on fit algorithm
  Bool_t            fUseRobust;          // switch on/off robust option, default: kFalse
  Double_t          fExpectedSigma;      // expected sigma to normalize robust fitting
  //  
  ClassDef(AliTMinuitToolkit,1);
};



#endif
