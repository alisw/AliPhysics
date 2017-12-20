#ifndef ALITMINUITTOOLKIT_H
#define ALITMINUITTOOLKIT_H

/// \ingroup STAT
/// \class AliTMinuitToolkit
/// \brief The AliTMinuitToolkit serves as an easy to use interface for the TMinuit
///
///#### Usage:
///  * 1.)  Setting the fit formula:
///    * Set fit function:
///      \code void SetFitFunction(TFormula * formula, Bool_t doReset)".\endcode
///  * 2.)   Adding the data points
///    * Direct specification of the points is possible via
///      \code void SetPoints(TMatrixD * points) \endcode
///      * Each point corresponds to one row in the matrix.
///    * The fitter is then started with  the command
///      \code void Fit()\endcode
///    * In order to fit a histogram, graphs in once use:
///      \code AliTMinuitToolkit::FitHistogram(const TH1 * his, Option_t *option) \endcode
///      \code AliTMinuitToolkit::FitGraph(const TGraph * gr, Option_t *option) \endcode
///      \code  AliTMinuitToolkit::Fit(TH1* his, const char *fitterName, Option_t* option, Option_t* goption,Option_t* foption, Double_t xMin, Double_t xMax)\endcode
///      \code AliTMinuitToolkit::Fit(TGraph *graph, const char *fitterName, Option_t* option, Option_t* goption,  Option_t* foption, Double_t xMin, Double_t xMax)\endcode
///    * Filling input points using TTree queries (example)
///      \code tool1D->FillFitter(inputTree,"test1D+noise:1/sqrt(12.+0)","testx0", "", 0,fitEntries); \endcode
///  * 3.)  Accessing the fit results
///    * The N parameters of the formula are stored in a N-dimensional vector which is returned by
///     * \code TVectorD * GetParameters() \endcode
///    * The covariance (NxN)  matrix of the fit
///     * \code  TMatrixD * GetCovarianceMatrix() \endcode
///  * 4.)  Log likelihood minimization
///    * see example \code $ALICE_ROOT/../src/STAT/test/AliTMinuitToolkitTest.C:Test1D()\endcode
///>     Even a few outliers can lead to wrong results of a least-squares fitting
///>     procedure. In this case the use of robust(resistant)  methods can be
///>     helpful, but a stronger dependence on starting values or convergence to local minimum can occur.
///    * Cost (log likelihood)  functions
///      -# By default chi2 minimization is invoked.
///      -# Optionally more robust log likelihood- Huber cost function can be specified instead of chi^2.
///       For small deviations  the function behaves like x^2 and for larger deviations like |x| - the so  called Huber estimator:
///       \code
///           h(x) = x^2                              , for x < 2.5*sigma
///           h(x) = 2*(2.5*sigma)*x - (2.5*sigma)^2  , for x > 2.5*sigma
///       \endcode
///       * REFERENCE for non-linear robust fitting:
///         * Ekblom H. and Madsen K. (1988), Algorithms for non-linear Huber estimation,
///         * BIT Numerical Mathematics 29 (1989) 60-76.
///         * internet: http://www.springerlink.com/content/m277218542988344/
///      -# User defined log likelihood function can be specified
///       \code
///       TF1 * f1Cost = new TF1("1","abs(x)<10?-log(0.8*exp(-x**2)+0.2/(1+x**2)):-log(0.2/(1+x**2))",-20,20); // 80 % gaus + 20% cachy
///       tool1D->SetLogLikelihoodFunction(f1Cost);
///       tool1D->Fit();
///       \endcode
///  * 5.)   Additional fitting strategy (options)
///     * Bootstrap
///     * MISAC (similar to RANSAC algorithm - https://en.wikipedia.org/wiki/Random_sample_consensus)
///     * cross validation
///     * Example usage (.x $AliRoot_SRC/STAT/test/AliTMinuitToolkitTestLinear.C+(2000,3))
///     \code
///    AliTMinuitToolkit::RegisterDefaultFitters();
///    AliTMinuitToolkit*fitter4 = AliTMinuitToolkit::Fit(gr,"pol1H","", "", "funOption(2,3,2)");
///    fit4=*(fitter4->GetParameters());     rms4=*(fitter3->GetCovarianceMatrix());
///    AliTMinuitToolkit*fitter5 = AliTMinuitToolkit::Fit(gr,"pol1R","bootstrap20", "","funOption(6,3,4)");  // bootstrap20
///    fit5=*(fitter5->GetParameters());     rms5=*(fitter5->GetRMSEstimator());
///    AliTMinuitToolkit*fitter6 = AliTMinuitToolkit::Fit(gr,"pol1R","bootstrap50", "", "funOption(3,3,5)" ); // bootstrap50
///    fit6=*(fitter6->GetParameters());     rms6=*(fitter6->GetRMSEstimator());
///    AliTMinuitToolkit*fitter7 = AliTMinuitToolkit::Fit(gr,"pol1R","misac(4,20)","","funOption(3,3,6)");    // misac 20
///    fit7=*(fitter7->GetParameters()); rms7=*(fitter7->GetRMSEstimator());
///    AliTMinuitToolkit*fitter8 = AliTMinuitToolkit::Fit(gr,"pol1R","misac(4,50)","","funOption(3,3,6)");    // misac 50
///    \endcode
///  * 6.)  Work in progress
///        *  Multidimensional fits (more than one observable)
///          *    e.g fitting Er and Erphi in parallel
///        *  Regularization
///          * e.g positive definite, monotonous function
///\author Marian Ivanov marian.ivanov@cern.ch
///
/// * Example usage in the $AliRoot_SRC/STAT/test/
///   * TODO - save example images in the doxygen/html directory (STEERING part ?) - currently saved by hand - to check with OFFLINE
///   * $AliRoot_SRC/STAT/test/AliTMinuitToolkitTest.C+
///   * $AliRoot_SRC/STAT/test/AliTMinuitToolkitTestLinear.C+
///     *  compare fits example with different likelihood/resp fits strategies for data with outliers
///        \image html AliTMinuiToolkit.TestLinearFitExample.png
///     *  compare performance fits example with different likelihood/resp fits strategies for data with outliers
///        \image html AliTMinuiToolkit.TestLinearFitStatExample.png

//--------------------------------------------------------------------------------------


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
    kSysInfo=0x00002,           // flag:print system Information for profiling (Memory, CPU...)
    kPrintAll=0x00004,          // flag:print detailed information in minuit
    //
            kStreamFcn=0x00010,         // flag: stream fcn values
    kStreamFcnPoint=0x00020,    // flag: stream fcn points
    kFCNInfo=0x000400           // flag print FCN information
  };
  enum EFitFlags{
    kCovarFailed=0x0001     // flag - problem to extract covariance matrix 
  };
  AliTMinuitToolkit(const char *streamerName=nullptr, const char *mode="recreate");
  virtual ~AliTMinuitToolkit();
  // streaming intermediate results
  void  SetStreamer(const char *streamerName, const char *mode="recreate");
  TTreeSRedirector *GetStreamer(){return fStreamer;}
  //
  void  ClearData();
  void  Fit(Option_t* option = "");
  void  FitHistogram(const TH1 *his,  Option_t* option = "");
  void  FitGraph(const TGraph * gr, Option_t* option = "");
  Long64_t FillFitter(TTree * inputTree, TString values, TString variables, TString selection, Int_t firstEntry, Int_t nEntries, Bool_t doReset=kTRUE);
  TString GetFitFunctionAsAlias(Option_t *option="", TTree * tree= nullptr);
  void Bootstrap(ULong_t nIter, const char* reportName, Option_t  *option= nullptr);
  void TwoFoldCrossValidation(UInt_t nIter, const char*reportName, Option_t *option= nullptr);
  void MISAC(Int_t nFitPoints, UInt_t nIter, const char*reportName, Option_t *option= nullptr);
  // 
  void SetFitFunction(TF1 * formula, Bool_t doReset);
  void SetInitialParam(const TMatrixD * initialParam);
  TMatrixD *  GetInitialParam() {return fInitialParam;}
  void SetLogLikelihoodFunction(TF1 *logLikelihood){fLogLikelihoodFunction=logLikelihood;}
  void SetFitAlgorithm(Char_t *const name) {fFitAlgorithm=name;};
  void SetMaxCalls(Int_t calls) {fMaxCalls = calls;};
  void SetTolerance(Double_t tol) {fPrecision = tol;};
  void SetPoints(TMatrixD * points) {fPoints = points;};
  void SetValues(TMatrixD * values) {fValues = values;};
  static void FitterFCN(int &nParam, double *info, double &chi2, double *gin, int iflag);
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
  // static fitting
  static AliTMinuitToolkit *  Fit(TH1* his, const char *fitterName, Option_t* fitOption = "", Option_t* hisOption = "", Option_t *funOption="", Double_t xMin = 0, Double_t xMax = 0);
  static AliTMinuitToolkit *  Fit(TGraph* graph, const char *fitterName, Option_t* fitOption = "", Option_t* graphOption = "", Option_t* funOption="",Double_t xMin = 0, Double_t xMax = 0);
  static AliTMinuitToolkit *  GetPredefinedFitter(const std::string & fitterName){return fPredefinedFitters[fitterName];}
  static void   SetPredefinedFitter(const std::string & fitterName, AliTMinuitToolkit * fitter){fPredefinedFitters[fitterName]=fitter;}
  static  void  RegisterDefaultFitters();
  static AliTMinuitToolkit * RegisterPlaneFitter(Int_t nPlanes, Int_t logLikeType);
  // helper functions for TFormula and TTreeFormula drawing
  static Double_t  GaussCachyLogLike(const Double_t *x, const Double_t *p);
  static Double_t  HuberLogLike(const Double_t *x, const Double_t *p);
  static Double_t  PseudoHuberLogLike(const Double_t *x, const Double_t *p);
  static Double_t  GausKurtosisSkewness(const Double_t *x, const Double_t *p);
  static Double_t  RndmGaus(Double_t mean=0, Double_t sigma=1);
  static Double_t  RndmLandau(Double_t mean=0, Double_t sigma=1);
  static void SetFunctionDrawOption(TF1 *fun, Option_t *option); // parse draw option and set function attributes
private:
  //
  AliTMinuitToolkit(const AliTMinuitToolkit&);            //  private copy constructor -- suppress warnings
  AliTMinuitToolkit& operator=(const AliTMinuitToolkit&); //  private assignment operator
  //
  TTreeSRedirector *fStreamer;              // !pointer to the streamer
  Int_t             fVerbose;               // verbosity flag
  TF1             * fFormula;               // formula of the fitted function
  TF1             * fLogLikelihoodFunction; // Log likelihood (const)  function
  TString           fFitAlgorithm;          // fit algorithm for TMinuit: migrad, simplex, ...
  TMatrixD        * fPoints;             // !points -  IsOwner
  TMatrixD        * fValues;             // !matrix of values and their errors -  IsOwner
  TObjArray       * fVarNames;           // !variable names  
  TObjArray       * fValueNames;         // value  names                    
  TArrayI           fPointIndex;         // point index - to specify points for likelihood calculation (used for Bootstrap and CrossValidation)
  TMatrixD        * fInitialParam;       // limits of the parameters (up, down)
  TVectorD        * fParam;              // parameter values
  TVectorD        * fRMSEstimator;       // parameter spread as estimated in Bootstrap resp. TwoFold cross validation
  TMatrixD        * fMISACEstimators;    // MISAC estimators
  TMatrixD        * fCovar;              // covariance matrix
  Int_t             fStatus;             // fStatus
  Double_t          fChi2;               // chi-square
  Int_t             fMaxCalls;           // maximum number of calls, default value depends on fit algorithm
  Double_t          fPrecision;          // normalized tolerance, default value depends on fit algorithm
  Bool_t            fIsHuberCost;        // switch on/off robust option, default: kFalse
  Int_t            fFCNCounter;         // function call counter
  static            std::map<std::string, AliTMinuitToolkit*> fPredefinedFitters;
ClassDef(AliTMinuitToolkit,2);
};



#endif
