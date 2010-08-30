#ifndef ALITMINUITTOOLKIT_H
#define ALITMINUITTOOLKIT_H


#include <TNamed.h>
class TH1F;
class TH2F;
class TFormula;
#include <TVectorDfwd.h>
#include <TMatrixDfwd.h>

class AliTMinuitToolkit: public TNamed{
public:
  //
  //
  //
  AliTMinuitToolkit();
  virtual ~AliTMinuitToolkit();
  
  void FitHistogram(TH1F *const his);
  void FitHistogram(TH2F *const his);
  void SetInitialParam(TVectorD *const param) { fParam=param;};
  void SetParamLimits(TMatrixD *const paramLimits) { fParamLimits=paramLimits;};
  //
  void SetFitFunction(TFormula *const formula) {fFormula=formula;};
  void SetWeightFunction(TFormula *const weightFunction) {fWeightFunction=weightFunction;};
  void SetWeightFunction(const Char_t *name, Float_t param1, Float_t param2 = 0.1);
  void SetFitAlgorithm(Char_t *const name) {fFitAlgorithm=name;};
  void SetMaxCalls(Int_t calls) {fMaxCalls = calls;};
  void SetTolerance(Double_t tol) {fPrecision = tol;};
  void SetPoints(TMatrixD *const points) {fPoints = points;};
  void SetWeights(TVectorD *const weights) {fWeights = weights;};
  static void Test();
  static void FitterFCN(int &npar, double *dummy, double &fchisq, double *gin, int iflag);
  void Fit();
  void EnableRobust(Bool_t b) {fUseRobust = b;};
  void EnableRobust(Bool_t b, Double_t sigma) {fUseRobust = b; fExpectedSigma = sigma;}; 
  //
  TMatrixD * GetPoints() const {return fPoints;};
  TVectorD * GetWeights() const {return fWeights;};
  TFormula * GetFormula() const {return fFormula;};
  TVectorD * GetParameters() const {return fParam;};
  TFormula * GetWeightFunction() const {return fWeightFunction;};
  TMatrixD * GetCovarianceMatrix() const {return fCovar;};
  Bool_t     GetStatus() const {return fUseRobust;};
  Double_t   GetExpectedSigma() const {return fExpectedSigma;};
 
private:
  //
  AliTMinuitToolkit(const AliTMinuitToolkit&); // fake copy constr. -- suppress warnings
  AliTMinuitToolkit& operator=(const AliTMinuitToolkit&); // fake -- suppress warnings
  //
  TFormula        * fFormula;            // formula of the fitted function
  TFormula        * fWeightFunction;     // weight function, must be defined between 0 and 1
  TString          fFitAlgorithm;       // fit algorithm for TMinuit: migrad, simplex, ...
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
