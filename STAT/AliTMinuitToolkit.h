#ifndef ALITMINUITTOOLKIT_H
#define ALITMINUITTOOLKIT_H


#include <TNamed.h>
class TH1F;
class TH2F;
class TFormula;
class TF1;
class TTree;
#include <TVectorDfwd.h>
#include <TMatrixDfwd.h>

class AliTMinuitToolkit: public TNamed{
public:
  //
  //
  //
  AliTMinuitToolkit();
  virtual ~AliTMinuitToolkit();  
  void  ClearData();
  void Fit();
  void  FitHistogram(TH1F *const his);
  //  Int_t UnbinnedFit(TTree * inputTree, TString values, TString variables, TString selection, Int_t firstEntry, Int_t nentries);
  Int_t FillFitter(TTree * inputTree, TString values, TString variables, TString selection, Int_t firstEntry, Int_t nentries);
  TString GetFitFunctionAsAlias();
  //
  void SetFitFunction(TFormula *const formula, Bool_t doReset);
  void SetInitialParam(TVectorD *const param);
  void SetParamLimits(TMatrixD *const paramLimits);  
  void SetLogLikelihoodFunction(TFormula *logLikelihood){fLogLikelihoodFunction=logLikelihood;}
  void SetFitAlgorithm(Char_t *const name) {fFitAlgorithm=name;};
  void SetMaxCalls(Int_t calls) {fMaxCalls = calls;};
  void SetTolerance(Double_t tol) {fPrecision = tol;};
  void SetPoints(TMatrixD *const points) {fPoints = points;};
  void SetValues(TMatrixD *const values) {fValues = values;};
  static void FitterFCN(int &npar, double *dummy, double &fchisq, double *gin, int iflag);
  void EnableRobust(Bool_t b) {fUseRobust = b;};
  //
  TMatrixD * GetPoints() const {return fPoints;};
  TMatrixD * GetValues() const {return fValues;};
  TFormula * GetFormula() const {return fFormula;};
  const TFormula *GetLogLikelihoodFunction() const { return fLogLikelihoodFunction;}
  const TVectorD * GetParameters() const {return fParam;};
  const TMatrixD * GetCovarianceMatrix() const {return fCovar;};
  Bool_t     GetStatus() const {return fUseRobust;};
  TObjArray  *GetListOfVariables(){return fVarNames;}
  void SetVerbose(Int_t verbose){fVerbose=verbose;}
  static Double_t RrndmGaus(Double_t mean=0, Double_t sigma=1);
  static Double_t RrndmLandau(Double_t mean=0, Double_t sigma=1);
private:
  //
  AliTMinuitToolkit(const AliTMinuitToolkit&);            // fake copy constr. -- suppress warnings
  AliTMinuitToolkit& operator=(const AliTMinuitToolkit&); // fake -- suppress warnings
  //
  Int_t           fVerbose;                 // verbosity flag
  TFormula        * fFormula;               // formula of the fitted function
  TFormula        * fLogLikelihoodFunction; // Log likelihood (const)  function
  TString           fFitAlgorithm;          // fit algorithm for TMinuit: migrad, simplex, ...
  //
  //
  TMatrixD        * fPoints;             // points -  IsOwner
  TMatrixD        * fValues;             // matrix of values and their errors -  IsOwner
  TObjArray       * fVarNames;           // variable names  
  TObjArray       * fValueNames;         // value  names  
                  
  //
  TVectorD        * fParam;              // parameter values
  TMatrixD        * fParamLimits;        // limits of the parameters (up, down)
  TMatrixD        * fCovar;              // covariance matrix
  Double_t          fChi2;               // chi-square
  Int_t             fMaxCalls;           // maximum number of calls, default value depends on fit algorithm
  Double_t          fPrecision;          // normalized tolerance, default value depends on fit algorithm
  Bool_t            fUseRobust;          // switch on/off robust option, default: kFalse
  //  
  ClassDef(AliTMinuitToolkit,1);
};



#endif
