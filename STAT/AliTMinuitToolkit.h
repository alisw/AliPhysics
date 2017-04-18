#ifndef ALITMINUITTOOLKIT_H
#define ALITMINUITTOOLKIT_H


#include <TNamed.h>
class TH1;
class TH1F;
class TH2F;
class TFormula;
class TF1;
class TTree;
class TGraph;
class TTreeSRedirector;
#include <TVectorDfwd.h>
#include <TMatrixDfwd.h>
#include <TArrayI.h>
#include <map>


class AliTMinuitToolkit: public TNamed{
public: 
  enum EVerboseFlags{           // flags for verbosity
    kPrintMinuit=0x00001,       // flag:print Minuit messages 
    kSysInfo=0x00002,           // flag:print SystemInfoamtion for profiling (Mempry, CPU...)
    kPrintAll=0x00004,          // flag:print detailed information in minuit
    //
    kStreamFcn=0x00010,         // flag: stream fcn values   
    kStreamFcnPoint=0x00020     // flag: stream fcn points   
  };
  enum EFitFlags{
    kCovarFailed=0x0001     // flag - problem to extract covariance matrix 
  };

  AliTMinuitToolkit(const char *streamerName=0);
  virtual ~AliTMinuitToolkit();  
  // streaming intermediate results
  void  SetStreamer(const char *streamerName);
  TTreeSRedirector *GetStreamer(){return fStreamer;}
  
  //
  void  ClearData();
  void  Fit(Option_t* option = "");
  void  FitHistogram(TH1 *const his,  Option_t* option = "");
  void  FitGraph(TGraph *const gr, Option_t* option = "");
  //  Int_t UnbinnedFit(TTree * inputTree, TString values, TString variables, TString selection, Int_t firstEntry, Int_t nentries);
  Int_t FillFitter(TTree * inputTree, TString values, TString variables, TString selection, Int_t firstEntry, Int_t nentries, Bool_t doReset=kTRUE);
  TString GetFitFunctionAsAlias();
  void Bootstrap(Int_t nIter, const char* reportName, Option_t  *option=0);
  void TwoFoldCrossValidation(Int_t nIter, const char*reportName, Option_t *option=0);
  void MISAC(Int_t nFitPoints, Int_t nIter, const char*reportName, Option_t *option=0);
  // 
  void SetFitFunction(TF1 *const formula, Bool_t doReset);
  void SetInitialParam(TMatrixD *const initialParam);  
  TMatrixD *  GetInitialParam() {return fInitialParam;}
  void SetLogLikelihoodFunction(TF1 *logLikelihood){fLogLikelihoodFunction=logLikelihood;}
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
  TF1 * GetFormula() const {return fFormula;};
  const TF1 *GetLogLikelihoodFunction() const { return fLogLikelihoodFunction;}
  const TVectorD * GetParameters() const {return fParam;};
  const TVectorD * GetRMSEstimator() const {return fRMSEstimator;};
  const TMatrixD * GetCovarianceMatrix() const {return fCovar;};
  const TMatrixD * GetMISACEstimators() const {return fMISACEstimators;};
  
  Double_t GetChisquare(){return fChi2;}
  Int_t GetStatus(){return fStatus;}
  Bool_t     IsHuberCost() const {return fIsHuberCost;};
  TObjArray  *GetListOfVariables(){return fVarNames;}
  void SetVerbose(Int_t verbose){fVerbose=verbose;}
  Bool_t IsVerbose(Int_t flag){ return (fVerbose&flag)>0;}
  // static fittig
  static AliTMinuitToolkit *  Fit(TH1* his, const char *fitterName, Option_t* option = "", Option_t* goption = "", Option_t *foption="", Double_t xmin = 0, Double_t xmax = 0);
  static AliTMinuitToolkit *  Fit(TGraph* graph, const char *fitterName, Option_t* option = "", Option_t* goption = "", Option_t* foption="",Double_t xmin = 0, Double_t xmax = 0);
  static AliTMinuitToolkit *  GetPredefinedFitter(std::string fitterName){return fPredefinedFitters[fitterName];}
  static void   SetPredefinedFitter(std::string fitterName, AliTMinuitToolkit * fitter){fPredefinedFitters[fitterName]=fitter;}
  static  void  RegisterDefaultFitters();
  // helper functions for TFormula and TTreeFormuladrawing
  static Double_t  GaussCachyLogLike(Double_t *x, Double_t *p);
  static Double_t RrndmGaus(Double_t mean=0, Double_t sigma=1);
  static Double_t RrndmLandau(Double_t mean=0, Double_t sigma=1);  
  static void SetFunctionDrawOption(TF1 *fun, Option_t *option); // parse draw option and set function attributes
private:
  //
  AliTMinuitToolkit(const AliTMinuitToolkit&);            // fake copy constr. -- suppress warnings
  AliTMinuitToolkit& operator=(const AliTMinuitToolkit&); // fake -- suppress warnings
  //
  TTreeSRedirector *fStreamer;              // !pointer to the streamer
  Int_t             fVerbose;               // verbosity flag
  TF1        * fFormula;               // formula of the fitted function
  TF1        * fLogLikelihoodFunction; // Log likelihood (const)  function
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
  TVectorD        * fRMSEstimator;       // parameter spread as estimated in Bootstrap resp. TwoFold cross validation
  TMatrixD        * fMISACEstimators;    // MISAC estimators
  TMatrixD        * fCovar;              // covariance matrix
  Int_t             fStatus;             // fstatus
  Double_t          fChi2;               // chi-square
  Int_t             fMaxCalls;           // maximum number of calls, default value depends on fit algorithm
  Double_t          fPrecision;          // normalized tolerance, default value depends on fit algorithm
  Bool_t            fIsHuberCost;        // switch on/off robust option, default: kFalse
  static            std::map<std::string, AliTMinuitToolkit*> fPredefinedFitters;  
  ClassDef(AliTMinuitToolkit,2);
};



#endif
