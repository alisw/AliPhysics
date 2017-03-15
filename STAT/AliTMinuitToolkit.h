#ifndef ALITMINUITTOOLKIT_H
#define ALITMINUITTOOLKIT_H


#include <TNamed.h>
class TH1F;
class TH2F;
class TFormula;
class TF1;
class TTree;
class TTreeSRedirector;
#include <TVectorDfwd.h>
#include <TMatrixDfwd.h>
#include <TArrayI.h>


class AliTMinuitToolkit: public TNamed{
public: 
  enum EVerboseFlags{ // flags for verbosity
    kPrintMinuit=0x00001,   // flag:print Minuit messages 
    kSysInfo=0x00002        // flag:print SystemInfoamtion for profiling (Mempry, CPU...)
  };
  AliTMinuitToolkit(const char *streamerName=0);
  virtual ~AliTMinuitToolkit();  
  void  ClearData();
  void  Fit(Bool_t doRest=kFALSE);
  void  FitHistogram(TH1F *const his);
  //  Int_t UnbinnedFit(TTree * inputTree, TString values, TString variables, TString selection, Int_t firstEntry, Int_t nentries);
  Int_t FillFitter(TTree * inputTree, TString values, TString variables, TString selection, Int_t firstEntry, Int_t nentries);
  TString GetFitFunctionAsAlias();
  void Bootstrap(Int_t nIter, const char* reportName);
  void TwoFoldCrossValidation(Int_t nIter, const char*reportName);
  //
  void SetFitFunction(TFormula *const formula, Bool_t doReset);
  void SetInitialParam(TMatrixD *const initialParam);  
  void SetLogLikelihoodFunction(TFormula *logLikelihood){fLogLikelihoodFunction=logLikelihood;}
  void SetFitAlgorithm(Char_t *const name) {fFitAlgorithm=name;};
  void SetMaxCalls(Int_t calls) {fMaxCalls = calls;};
  void SetTolerance(Double_t tol) {fPrecision = tol;};
  void SetPoints(TMatrixD *const points) {fPoints = points;};
  void SetValues(TMatrixD *const values) {fValues = values;};
  static void FitterFCN(int &npar, double *info, double &fchisq, double *gin, int iflag);
  void EnableRobust(Bool_t b) {fIsHuberCost = b;};
  //
  const TMatrixD * GetPoints() const {return fPoints;};
  const TMatrixD * GetValues() const {return fValues;};
  TFormula * GetFormula() const {return fFormula;};
  const TFormula *GetLogLikelihoodFunction() const { return fLogLikelihoodFunction;}
  const TVectorD * GetParameters() const {return fParam;};
  const TMatrixD * GetCovarianceMatrix() const {return fCovar;};
  Int_t GetStatus(){return fStatus;}
  Bool_t     IsHuberCost() const {return fIsHuberCost;};
  TObjArray  *GetListOfVariables(){return fVarNames;}
  void SetVerbose(Int_t verbose){fVerbose=verbose;}
  static Double_t RrndmGaus(Double_t mean=0, Double_t sigma=1);
  static Double_t RrndmLandau(Double_t mean=0, Double_t sigma=1);  
  const TTreeSRedirector *GetStreamer(){return fStreamer;}
  static Double_t  GaussCachyLogLike(Double_t *x, Double_t *p);
private:
  //
  AliTMinuitToolkit(const AliTMinuitToolkit&);            // fake copy constr. -- suppress warnings
  AliTMinuitToolkit& operator=(const AliTMinuitToolkit&); // fake -- suppress warnings
  //
  TTreeSRedirector *fStreamer;              // !pointer to the streamer
  Int_t             fVerbose;               // verbosity flag
  TFormula        * fFormula;               // formula of the fitted function
  TFormula        * fLogLikelihoodFunction; // Log likelihood (const)  function
  TString           fFitAlgorithm;          // fit algorithm for TMinuit: migrad, simplex, ...  
  //
  TMatrixD        * fPoints;             // !points -  IsOwner
  TMatrixD        * fValues;             // !matrix of values and their errors -  IsOwner
  TObjArray       * fVarNames;           // !variable names  
  TObjArray       * fValueNames;         // value  names                    
  TArrayI           fPointIndex;         // point index - to specify points for likelihood calculation (used for Bottstrap and CrossValidation)             
  //
  TMatrixD        * fInitialParam;       // limits of the parameters (up, down)
  TVectorD        * fParam;              // parameter values
  TMatrixD        * fCovar;              // covariance matrix
  Int_t             fStatus;             // fstatus
  Double_t          fChi2;               // chi-square
  Int_t             fMaxCalls;           // maximum number of calls, default value depends on fit algorithm
  Double_t          fPrecision;          // normalized tolerance, default value depends on fit algorithm
  Bool_t            fIsHuberCost;        // switch on/off robust option, default: kFalse
  //  
  ClassDef(AliTMinuitToolkit,2);
};



#endif
