/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TMatrix.h>
#include <TRandom.h>
#include <TVector.h>
#include <TVectorD.h>
#include <TVirtualFitter.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTreeStream.h>
#include "AliTMinuitToolkit.h"
#include "TGraph.h"
#include "TPRegexp.h"
#include "TStatToolkit.h"
//#include "TFormulaPrimitive.h"
#include "TDataMember.h"

std::map<std::string, AliTMinuitToolkit*> AliTMinuitToolkit::fPredefinedFitters;


ClassImp(AliTMinuitToolkit)

AliTMinuitToolkit::AliTMinuitToolkit(const char *streamerName, const char *mode) :
        TNamed(),
        fStreamer(nullptr),
        fVerbose(0),
        fFormula(nullptr),
        fLogLikelihoodFunction(nullptr),
        fFitAlgorithm(""),
        fPoints(nullptr),
        fValues(nullptr),
        fVarNames(nullptr),           // variable names
        fValueNames(nullptr),         // value  names
        fPointIndex(0),              // point index - to specify points for likelihood calculation (used for Bootstrap and CrossValidation)
        fInitialParam(nullptr),
        fParam(nullptr),
        fRMSEstimator(nullptr),       // parameter spread as estimated in Bootstrap resp. TwoFold cross validation
        fMISACEstimators(nullptr),    // MISAC estimators - median, mean, rms
        fCovar(nullptr),
        fStatus(0),
        fChi2(0),
        fMaxCalls(0),
        fPrecision(0),
        fIsHuberCost(0),
        fFCNCounter(0)                // cost function call counter
{
  //
  // standard constructor
  //
  fMaxCalls = 500;
  fPrecision = 1;
  fIsHuberCost = false;
  if (streamerName!= nullptr) fStreamer= new TTreeSRedirector(streamerName,mode);
}



AliTMinuitToolkit::~AliTMinuitToolkit(){
  //
  // destructor
  //
  delete fPoints;
  delete fValues;
  delete fInitialParam;
  //   delete fFormula;
  delete fParam;
  delete fRMSEstimator;
  delete fMISACEstimators;
  delete fCovar;
  delete fStreamer;
}
/// \brief Set stream for  writing info into TTree.
///        see examples(STEER/STEERBase/TTreeStream.cxx)
/// \param streamerName - Name of .root file
/// \param mode - option for .root file open 
///        * recreate
///        * update
///        * ....
void  AliTMinuitToolkit::SetStreamer(const char *streamerName, const char *mode){
  //
  // set streamer
  delete fStreamer;
  fStreamer= new TTreeSRedirector(streamerName,mode);
}


void  AliTMinuitToolkit::ClearData(){
  delete fPoints;
  fPoints= nullptr;
  delete fValues;
  fValues= nullptr;
  //  delete fParam;
  //   fParam=0;
  //   delete fCovar;
  //   fCovar=0;
}
/// \brief Set formula of the fitted function and reset parameter. See doReset
/// \param formula - formula of the fitted function
/// \param doReset - if true reset parameters related to the fit function
void  AliTMinuitToolkit::SetFitFunction(TF1 * formula, Bool_t doReset) {
  //
  fFormula=formula;
  if (doReset){
    delete fParam;
    delete fRMSEstimator;
    delete fCovar;
    fParam=new TVectorD(formula->GetNpar());
    fRMSEstimator=new TVectorD(formula->GetNpar());
    fCovar=new TMatrixD(formula->GetNpar(),formula->GetNpar());
  }
}

/// \brief Set initial parameters matrix
/// \param paramLimits - Initialization matrix
///        * number of lines = number of  parameters
///        * 4 columns
///          * [0]      - initial value
///          * [1]      - initial error
///          * [2],[3]  - optional parameters to specify minimum and maximum of the parameter (see TMinuit)    
void   AliTMinuitToolkit::SetInitialParam(const TMatrixD * paramLimits) {
  fInitialParam=new TMatrixD(*paramLimits);
};


/// \brief Fit a one dimensional histogram
/// \param his    -   const TH1 object
/// \param option -   fit option e.g (default, bootstrap20,misac ... )
/// * In alternative usage more option and histogram range can be specified
///    \code Fit(TH1* his, const char *fitterName, Option_t* fitOption = "", Option_t* hisOption = "", Option_t *funOption="", Double_t xMin = 0, Double_t xMax = 0)\endcode
void AliTMinuitToolkit::FitHistogram(const TH1 * his, Option_t *option) {
  ClearData();
  fPoints  = new TMatrixD(his->GetNbinsX(), 1);
  fValues  = new TMatrixD(his->GetNbinsX(), 2);
  for(Int_t iBin=0; iBin < his->GetNbinsX(); iBin++) {
    Double_t x = his->GetXaxis()->GetBinCenter(iBin+1);
    Double_t y = his->GetBinContent(iBin+1);
    (*fPoints)(iBin, 0) = x;
    Double_t err=his->GetBinError(iBin+1);
    (*fValues)(iBin,0)=y;
    (*fValues)(iBin,1)=(err>0)? 1./err:0;
  }
  Fit(option);
}
/// \brief Fit a TGraph object
/// \param gr     - const TGraph object
/// \param option - fit option e.g (default, bootstrap20,misac .... )
/// * In alternative usage more option and graph range can be specified
///    \code Fit(TGraph* graph, const char *fitterName, Option_t* fitOption = "", Option_t* graphOption = "", Option_t* funOption="",Double_t xMin = 0, Double_t xMax = 0)\endcode
void AliTMinuitToolkit::FitGraph(const TGraph * gr, Option_t *option) {
  ClearData();
  Int_t nPoints=gr->GetN();
  fPoints  = new TMatrixD(nPoints, 1);
  fValues  = new TMatrixD(nPoints, 2);
  for(Int_t ip=0; ip < nPoints; ip++) {
    Double_t x = gr->GetX()[ip];
    Double_t y = gr->GetY()[ip];
    (*fPoints)(ip, 0) = x;
    Double_t err=gr->GetErrorY(ip);
    (*fValues)(ip,0)=y;
    (*fValues)(ip,1)=(err>0)? 1./err:0;
  }
  Fit(option);
}

/// \brief Internal function used in the TMinuit minimization
/// * see  documentation for FCN parameters g.g in ROOT- TMinuit::Eval
///
/// \param nPar          - not used in our interface
/// \param grad    - gradient of cost function
/// \param chi2(fVal))   - cost function to minimize (chi2/log-likelihood/)
/// \param gin(par)      - function parameters
/// \param iflag(flag)   - flag  used by root (1.) print/ 2.) derivative/ 3. final fit)
///                      - in our interface also possible to specify stream data (kStreamFcnPoint,kStreamFcn)
/// typedef void(* 	FCNFunc_t) (Int_t &nPar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag)
///
///*   Original documentation in TMinuit::Eval documentation in root
///  *  Input parameters:
///    *  nPar:    number of currently variable parameters
///    *  par:     array of (constant and variable) parameters
///    *  flag:    Indicates what is to be calculated (see example below) ()
///    *  grad:    array of gradients
///  * Output parameters:
///    *  fVal:    The calculated function value.
///    *  grad:    The (optional) vector of first derivatives).
void AliTMinuitToolkit::FitterFCN(int &/*nPar*/, double */*grad*/, double &chi2, double *gin, int iflag){
  AliTMinuitToolkit * fitter = (AliTMinuitToolkit*)TVirtualFitter::GetFitter()->GetObjectFit();
  if (!fitter) return;
  Int_t &fcnCounter=fitter->fFCNCounter;
  TF1 *logLike=fitter->fLogLikelihoodFunction;
  const Double_t* likeParam= (logLike!= nullptr) ? logLike->GetParameters(): nullptr;
  chi2 = 0;
  const TMatrixD & variables= (*fitter->GetPoints());
  const TMatrixD & values=    (*fitter->GetValues());
  Int_t nVars       = variables.GetNcols();
  Int_t nPoints     = variables.GetNrows();
  TVectorD param(fitter->fFormula->GetNpar(), gin);
  // calculate  log likelihood
  //
  if (fitter->fPointIndex.GetSize()>0) {
    nPoints=fitter->fPointIndex.GetSize();
  }
  for (Int_t jPoint=0; jPoint<nPoints; jPoint++){
    Int_t iPoint=(fitter->fPointIndex.GetSize()>0) ? fitter->fPointIndex[jPoint]:jPoint;
    Double_t x[100];
    for (Int_t ivar=0; ivar<nVars; ivar++){
      x[ivar] = variables(iPoint, ivar);
    }
    Double_t funX = (fitter->GetFormula())->EvalPar(x,gin);
    Double_t value=values(iPoint,0);
    Double_t weight=values(iPoint,1);
    Double_t delta = (value - funX);
    Double_t toAdd=0;
    if (logLike){
      Double_t normDelta = delta*weight;
      toAdd=logLike->EvalPar(&normDelta,likeParam);
    }else{
      if (fitter->IsHuberCost() != 0) {       //hubert norm
        delta = TMath::Abs(delta*weight);                   // normalization
        if (delta <= 2.5) toAdd= delta*delta; // new metric: Huber-k-estimator
        if (delta > 2.5)  toAdd= 2*(2.5)*delta - (2.5*2.5);
      } else {
        Double_t lChi2 = delta*weight;
        lChi2*=lChi2;
        toAdd=lChi2;   // chi2 (log likelihood)
      }
    }
    //
    if (fitter->IsVerbose(kStreamFcnPoint) || (iflag&kStreamFcnPoint)){
      TVectorD vecX(nVars, x);
      (*(fitter->fStreamer))<<"fcnDebugPoint"<<
                            "toAdd="<<toAdd<<  // log likelihood to add
                            "val="<<value<<    // "measurement"
                            "w="<<weight<<     // wight of point
                            "fun="<<funX<<     // function value
                            "x.="<<&vecX<<     // variable vector
                            "p.="<<&param<<    // parameter vector
                            "\n";
    }
    chi2+=toAdd;
  }
  if (fitter->IsVerbose(kStreamFcn) || (iflag&kStreamFcn )){
    (*(fitter->fStreamer))<<"fcnDebug"<<
                          "chi2="<<chi2<<
                          "p.="<<&param<<
                          "\n";
  }
  fcnCounter++;
}

/// \brief Internal function that calls the fitter
/// \param option:
///          * bootstrap[0-9]* -
///          * twofold[0-9]* -
///          * misac\([0-9]*,[0-9]*\) -
void AliTMinuitToolkit::Fit(Option_t *option) {
  TString sOption=(option== nullptr)?"":option;
  TPRegexp regBootstrap("bootstrap[0-9]*");
  TPRegexp regTwofold("twofold[0-9]*");
  TPRegexp regMISAC("misac\\([0-9]*,[0-9]*\\)");

  // run bootstrap if specified
  if (regBootstrap.Match(sOption)>0){      // run bootstrap
    TString newOption=sOption;
    regBootstrap.Substitute(newOption,"");
    UInt_t nPoints= static_cast<UInt_t>(TString(&(sOption.Data()[sOption.Index(regBootstrap) + 9])).Atoi());
    return Bootstrap(nPoints, nullptr,newOption.Data());
  }
  // run two fold minimization if specified
  if (regTwofold.Match(sOption)>0){      // run twofold minimization
    TString newOption=sOption;
    regTwofold.Substitute(newOption,"");
    UInt_t nPoints= static_cast<UInt_t>(TString(&(sOption.Data()[sOption.Index(regTwofold) + 7])).Atoi());
    return TwoFoldCrossValidation(nPoints, nullptr,newOption.Data());
  }
  // run MISAC
  if (regMISAC.Match(sOption)>0){      // run MISAC minimization  - mix of RANSAC and KFold
    TString newOption=sOption;
    regMISAC.Substitute(newOption,"");
    Int_t index1=sOption.Index(regMISAC)+6;
    Int_t index2=sOption.Index(",",index1+1)+1;
    Int_t nPoints=TString(&(sOption.Data()[index1])).Atoi();
    UInt_t nIter= static_cast<UInt_t>(TString(&(sOption.Data()[index2])).Atoi());
    return MISAC(nPoints,nIter, nullptr,newOption.Data());
  }
  fStatus=0;

  Int_t nParam  = fFormula->GetNpar();
  // create "initial param" in case not set before
  if (fParam== nullptr) fParam=new TVectorD(nParam);
  if (fInitialParam == nullptr) {
    ::Info("AliTMinuitToolkit::Fit","Default initial parameters set");
    fInitialParam=  new TMatrixD(nParam,4);
    for (Int_t iParam=0; iParam<nParam; iParam++) {
      (*fInitialParam)(iParam,0)=0;
      (*fInitialParam)(iParam,1)=1000;
      (*fInitialParam)(iParam,2)=0;
      (*fInitialParam)(iParam,3)=0;
    }
  }

  // migrad fit algorithm as default
  if (fFitAlgorithm == "") {
    fFitAlgorithm = "migrad";
  }

  // assign weights
  if (fValues == nullptr) {
    ::Error("AliTMinuitToolkit::Fit()","Invalid fit. Values not set");
    return ;
  }

  // set up the fitter
  TVirtualFitter *minuit = TVirtualFitter::Fitter(nullptr, nParam);
  minuit->SetObjectFit(this);
  minuit->SetFCN(AliTMinuitToolkit::FitterFCN);
  if ((fVerbose& kPrintMinuit)==0){  // MAKE minuit QUIET!!
    double p1 = -1;
    minuit->ExecuteCommand("SET PRINTOUT",&p1,1);
  }
  if ((fVerbose& kPrintAll)>0){  // MAKE minuit verbose!!  
    double pprint=2;
    minuit->ExecuteCommand("SET PRINTOUT",&pprint,1);
  }

  // initialize parameters (step size???)
  for (Int_t iParam=0; iParam<nParam; iParam++){
    if (sOption.Contains("rndmInit")){
      Double_t initValue=(*fInitialParam)(iParam, 0)+gRandom->Gaus()*(*fInitialParam)(iParam,1);
      minuit->SetParameter(iParam, Form("p[%d]",iParam), initValue, (*fInitialParam)(iParam,1), (*fInitialParam)(iParam, 2), (*fInitialParam)(iParam, 3));
    }else{
      minuit->SetParameter(iParam, Form("p[%d]",iParam), (*fInitialParam)(iParam,0), (*fInitialParam)(iParam,1), (*fInitialParam)(iParam, 2), (*fInitialParam)(iParam, 3));
    }
    /*
      if (doReset){
      minuit->SetParameter(iParam, Form("p[%d]",iParam), (*fInitialParam)(iParam,0), (*fInitialParam)(iParam,1), (*fInitialParam)(iParam, 2), (*fInitialParam)(iParam, 3));
      }else{
      minuit->SetParameter(iParam, Form("p[%d]",iParam), (*fParam)(iParam), TMath::Sqrt((*fCovar)(iParam,iParam)), (*fInitialParam)(iParam, 2), (*fInitialParam)(iParam, 3));
      }
    */
  }

  //
  Double_t argList[2];
  argList[0] = fMaxCalls; //maximal number of calls 
  argList[1] = fPrecision; //tolerance normalized to 0.001 
  //if (fMaxCalls == 500 && fPrecision == 1) minuit->ExecuteCommand(fFitAlgorithm, 0, 0); 
  //if (fMaxCalls != 500 || fPrecision != 1) 
  minuit->ExecuteCommand(fFitAlgorithm, argList, 2);
  // two additional arguments can be specified ExecuteCommand("migrad", argList, 2) - use 0,0 for default
  //  fStatus = (((TFitter*)minuit)->GetMinuit())->GetStatus(); // we should add separate status for Minuit and  AliTMinuitToolkit

  // fill parameter vector
  for (Int_t ivar=0; ivar<nParam; ivar++){
    (*fParam)(ivar) = minuit->GetParameter(ivar);
    fFormula->SetParameter(ivar, minuit->GetParameter(ivar));
    //fFormula->SetPar(ivar, minuit->GetParameter(ivar));
  }
  FitterFCN(nParam, nullptr,fChi2, fParam->GetMatrixArray(),0);

  // fill parameter vector
  for (Int_t ivar=0; ivar<nParam; ivar++){
    (*fParam)(ivar) = minuit->GetParameter(ivar);
    fFormula->SetParameter(ivar, minuit->GetParameter(ivar));
  }

  // fill covariance matrix
  if (fCovar== nullptr)   fCovar = new TMatrixD(nParam, nParam);
  if (minuit->GetCovarianceMatrix()){
    fCovar->SetMatrixArray(minuit->GetCovarianceMatrix());
  }else{
    fStatus|=kCovarFailed;
  }

}

/// \brief Fill input matrices of the fitter
/// * TODO - add better working example (using code/trees from the test directory)
/// \param inputTree      -  pointer to input tree
/// \param values         - values description (parameter:<error>)
/// \param variables      - variables array (var0:<var1>)  
/// \param selection      - tree selection 
/// \param firstEntry     - first entry to query
/// \param nEntries       - number of entries to query
/// \param doReset        - switch reset previous/or append
/// \return
/*!
   * Example: TPC occupancy 6 dimensional hyperplane fit
   \code
      AliTMinuitToolkit * hyp6R= AliTMinuitToolkit::GetPredefinedFitter("hyp6R");
      hyp6R->FillFitter(treeQA, "NoThreshold.fElements:0.1*NoThreshold.fElements", "1:dRNorm:Phi2Norm:QNorm:gxNorm:gyNorm", "region&&occFitCut&&rndm<0.2", 0, 10000000);
      hyp6R->Fit();
   \endcode

*/
Long64_t AliTMinuitToolkit::FillFitter(TTree * inputTree, TString values, TString variables, TString selection, Int_t firstEntry, Int_t nEntries, Bool_t doReset ){
  TString query=values;
  query+=":";
  query+=variables;
  if (inputTree== nullptr){
    ::Error("AliTMinuitToolkit::FillFitter","Input tree==NULL");
    return -1;
  }
  fValueNames=values.Tokenize(":");
  fVarNames=variables.Tokenize(":");
  Int_t nVal= fValueNames->GetEntries();
  Int_t nVars= fVarNames->GetEntries();
  if (doReset == 0 && fPoints != nullptr){
    if (fPoints->GetNrows()!=nVars){
      ::Error("AliTMinuitToolkit::FillFitter","Non compatible number of variables: %d instead of %d: variables:  %s",nVars, fPoints->GetNrows(), query.Data());
      return -1;
    }
  }

  Long64_t entries = inputTree->Draw(query.Data(),selection.Data(),"goff para",nEntries,firstEntry);
  if (entries<=0) {
    ::Error("AliTMinuitToolkit::FillFitter","badly formatted values or variables: %s",query.Data());
    ::Error("AliTMinuitToolkit::FillFitter","valueDescription: %s",values.Data());
    ::Error("AliTMinuitToolkit::FillFitter","variables: %s",variables.Data());
    return -1;
  }
  Int_t index0=0;
  if (doReset || fPoints== nullptr) {
    ClearData();
    fPoints=new TMatrixD(entries,nVars);
    fValues=new TMatrixD(entries,nVal);
  }else{
    index0= fPoints->GetNrows();
    fPoints->ResizeTo(index0+Int_t(entries),nVars);
    fValues->ResizeTo(index0+Int_t(entries),nVal);
  }
  for (Int_t iPoint=0; iPoint<entries; iPoint++){
    for (Int_t iVar=0; iVar<nVars; iVar++){
      (*fPoints)(index0+iPoint,iVar)=inputTree->GetVal(iVar+nVal)[iPoint];
    }
    for (Int_t iVal=0; iVal<nVal; iVal++){
      (*fValues)(index0+iPoint,iVal)=inputTree->GetVal(iVal)[iPoint];
    }
  }
  return entries;
}

/// \brief Format alias string for the for tree alias or for fit title
/// \param option  - use latex in case of title request
/// \return        - string description
TString AliTMinuitToolkit::GetFitFunctionAsAlias(Option_t *option, TTree * tree){
  //
  // construct string TTree alias for fit function
  TString inputString(fFormula->GetTitle());
  TString sOption(option);
  sOption.ToLower();

  if (fVarNames== nullptr){
    ::Error("AliTMinuitToolkit::GetFitFunctionAsAlias","Variable names not defined");
    return "";
  }

  for (Int_t iDim=0; iDim<fFormula->GetNdim(); iDim++){
    TString varName=GetListOfVariables()->At(iDim)->GetName();
    if (sOption.Contains("latex") &&tree){
      TNamed * meta = TStatToolkit::GetMetadata(tree,varName+".Title");
      if (meta) varName=meta->GetTitle();
    }
    inputString.ReplaceAll(TString::Format("x[%d]",iDim).Data(), varName.Data());
  }
  for (Int_t iPar=0; iPar<fFormula->GetNpar(); iPar++){
    if (sOption.Contains("latex")==0){
      inputString.ReplaceAll(TString::Format("[%d]",iPar).Data(), TString::Format("(%f)",(*fParam)[iPar]).Data());
    }else{
      if ((*GetRMSEstimator())[0]>0) {
        inputString.ReplaceAll(TString::Format("[%d]", iPar).Data(),
                               TString::Format("(%2.2f#pm%2.2f)", (*fParam)[iPar], (*GetRMSEstimator())[iPar]).Data());
      }else{
        inputString.ReplaceAll(TString::Format("[%d]", iPar).Data(),
                               TString::Format("(%2.2f#pm%2.2f)", (*fParam)[iPar], (*GetRMSEstimator())[iPar]).Data());
      }
    }
  }
  return inputString;
}

/// \brief Random gaus generator as static function (to use in root TFormula, TTreeFormula)
/// \param mean    - mean Gaus
/// \param sigma   - sigma gaus
/// \return
Double_t AliTMinuitToolkit::RndmGaus(Double_t mean, Double_t sigma){
  return gRandom->Gaus(mean,sigma);
}

/// \brief Random gaus generator as static function (to use in root TFormula, TTreeFormula)
/// \param mean
/// \param sigma
/// \return
Double_t  AliTMinuitToolkit::RndmLandau(Double_t mean, Double_t sigma){
  return gRandom->Landau(mean,sigma);
}


/// Log likelihood for the gaus+cauchy - used for Log likelihood minimization of distribution with tails
/// \param x   -
/// \param p   - p[0] - gaus fraction
///            - p[1] - cauchy 0.5*FWHM
/// \return    log likelihood for the mixture of gaus and cauchy distribution
///            representation as in https://en.wikipedia.org/wiki/Cauchy_distribution
Double_t  AliTMinuitToolkit::GaussCachyLogLike(const Double_t *x, const Double_t *p){
  Double_t vCauchy=p[1]/(TMath::Pi()*(p[1]*p[1]+(*x)*(*x)));
  Double_t vGaus=(TMath::Abs(*x)<20) ? TMath::Gaus(*x,0,1.,kTRUE):0.;
  Double_t p0= p[0]*TMath::Gaus(0,0,1,kTRUE)+(1-p[0])/(TMath::Pi()*p[1]);
  return -2.*TMath::Log((p[0]*vGaus+(1-p[0])*vCauchy)/p0);
}

/// Huber loss  is a loss function used in robust regression - representation as in https://en.wikipedia.org/w/index.php?title=Huber_loss&diff=809252384&oldid=798150659
/// \param x   -
/// \param p   - p[0] - gaus region boundary
/// \return    - loss  function (used as a log likelihood)
///
Double_t  AliTMinuitToolkit::HuberLogLike(const Double_t *x, const Double_t *p){
  Double_t absX=TMath::Abs(x[0]);
  if (absX<p[0]) return absX*absX;
  return p[0]*(2*absX-p[0]);
}

/// Pseudo Huber loss  is a loss function used in robust regression - representation as in https://en.wikipedia.org/w/index.php?title=Huber_loss&diff=809252384&oldid=798150659
/// \param x   -
/// \param p   - p[0] - gaus region boundary
/// \return    - loss  function (used as a log likelihood)
///
Double_t  AliTMinuitToolkit::PseudoHuberLogLike(const Double_t *x, const Double_t *p){
  Double_t p2=p[0]*p[0];
  Double_t x2=x[0]*x[0];
  Double_t result=2.*p2*(TMath::Sqrt(1.+x2/p2)-1.);
  return result;
}

/// \brief Gaus convoluted with Erf to emulate kurtosis and skewness
///
/// \param x
/// \param p
///          * p[0],p[1],p[2] - as for gauss
///          * p[3] - Kurtosis  - tails - proportional to 4 moment
///          * p[4] - Asymmetry  - proportional to 3 moment
/// \return
Double_t AliTMinuitToolkit::GausKurtosisSkewness(const Double_t *x, const Double_t *p){
  // TODO add protection against out of range
  static Float_t msqrt2=1./TMath::Sqrt(2.);
  Double_t gaussPart=p[0]*TMath::Exp(-(TMath::Power(TMath::Abs(x[0]-p[1])/p[2],p[3])));
  Double_t erfPart=(1.+TMath::Erf(p[4]*(x[0]-p[1])/p[2]*msqrt2));
  return gaussPart*erfPart;
}



/// \brief Bootstrap minimization
/// * fitting parameters done several times on modified data sample (random samples with replacement)
///  *   to emulate different data population
///  *   reduce problem with local minimum in the M-estimator estimated parameters - robust mean/median
///  *   obtained RMS used to estimate error (take into account also model error)
/// * TODO - parameters robust estimator setting
/// \param nIter        - number of iteration
/// \param reportName   - report name - in case results streamed
/// \param option       - fit option

void AliTMinuitToolkit::Bootstrap(ULong_t nIter, const char * reportName, Option_t *option){
  Int_t nPoints= fPoints->GetNrows();
  fPointIndex.Set(nPoints);
  // Double_t info[3]={0,0,0};
  Int_t nPar= fFormula->GetNpar();
  Int_t fcnP=-1;
  TObjString rName("bootstrap");
  if (reportName!= nullptr) rName=reportName;
  std::vector< std::vector<double> > vecPar(static_cast<unsigned long>(nPar), std::vector<double>(nIter));
  for (UInt_t iter=0; iter<nIter; iter++){
    for (Int_t iPoint=0; iPoint<nPoints; iPoint++){
      fPointIndex[iPoint]= static_cast<Int_t>(gRandom->Rndm() * nPoints);
    }
    Fit(option);
    FitterFCN(nPar, nullptr,fChi2, fParam->GetMatrixArray(),0);
    if (fcnP<0) fcnP=fFCNCounter;
    Int_t fcnCounter=fFCNCounter-fcnP;
    fcnP=fFCNCounter;
    for (Int_t iPar=0; iPar<nPar;iPar++)  vecPar[iPar][iter]=(*fParam)(iPar);
    if (fStreamer){
      (*fStreamer)<<"bootstrap"<<
                  "iter="<<iter<<
                  "fStatus="<<fStatus<<    // minuit status
                  "reportName.="<<&rName<<
                  "nPoints="<<nPoints<<
                  "fChi2="<<fChi2<<
                  "fcnCounterI="<<fFCNCounter<<
                  "fcnCounter="<<fcnCounter<<
                  "fParam.="<<fParam<<
                  "fCovar.="<<fCovar<<
                  "\n";
    }
  }
  if (fRMSEstimator== nullptr) fRMSEstimator=new TVectorD(*fParam);
  TVectorD &par=*fParam;
  TVectorD &rms=*fRMSEstimator;
  for (Int_t iPar=0; iPar<nPar;iPar++) {
    Double_t rmsF,meanF;
    TStatToolkit::EvaluateUni(static_cast<Int_t>(nIter), vecPar[iPar].data(), meanF, rmsF,
                              static_cast<Int_t>(TMath::Max(nIter * 0.75, 1.)));
    par[iPar]=TMath::Median(static_cast<Long64_t>(nIter), vecPar[iPar].data());
    rms[iPar]=rmsF;
  }
  fPointIndex.Set(0);
}


///\breief Two fold cross validation - https://en.wikipedia.org/wiki/Cross-validation_(statistics)
/// * Sample randomly divided to the independent training set and  complementary test set (Repeated random sub-sampling validation)
/// * log likelihood cost function calculated for training fit set
/// * TODO - Used for fit model validation. Report - Post processing using streamer not available for the moment
/// * TODO - Add:  The results are then averaged over the splits as described in the wiki section - (Repeated random sub-sampling validation)
/// \param nIter        - number of iterations
/// \param reportName   - report name (identifier in the streamer)
/// \param option       - fit option
void AliTMinuitToolkit::TwoFoldCrossValidation(UInt_t nIter, const char *reportName, Option_t *option) {
  Int_t nPoints = fPoints->GetNrows();
  TArrayI indexFold(nPoints);
  TArrayD rndmCV(nPoints);
  fPointIndex.Set(nPoints);
  Double_t chi2_00, chi2_01, /*chi2_10,*/ chi2_11;
  //Double_t info[3] = {0, 0, 0};
  Int_t nPar = fFormula->GetNpar();
  Int_t fcnP = -1;
  TObjString rName(reportName);
  std::vector<std::vector<double> > vecPar(static_cast<unsigned long>(2 * nPar + 4), std::vector<double>(nIter));  //store vectors and chi2
  for (UInt_t iter = 0; iter < nIter; iter++) {
    for (Int_t iPoint = 0; iPoint < nPoints; iPoint++) {
      rndmCV[iPoint] = gRandom->Rndm() * nPoints;
    }
    TMath::Sort(nPoints, rndmCV.GetArray(), indexFold.GetArray());
    //
    fPointIndex.Set(nPoints / 2, indexFold.GetArray());
    Fit(option);
    for (Int_t iPar = 0; iPar < nPar; iPar++) vecPar[iPar][iter] = (*fParam)(iPar);
    TVectorD param0(*fParam);
    TMatrixD covar0(*fCovar);
    FitterFCN(nPar, nullptr, chi2_00, fParam->GetMatrixArray(), 0);
    //
    fPointIndex.Set(nPoints - nPoints / 2, &(indexFold.GetArray()[nPoints / 2]));
    FitterFCN(nPar, nullptr, chi2_01, fParam->GetMatrixArray(), 0);
    Fit(option);
    for (Int_t iPar = 0; iPar < nPar; iPar++) vecPar[nPar + iPar][iter] = (*fParam)(iPar);
    FitterFCN(nPar, nullptr, chi2_11, fParam->GetMatrixArray(), 0);
    vecPar[2 * nPar][iter] = chi2_00 + chi2_01 + chi2_11;
    TVectorD param1(*fParam);
    TMatrixD covar1(*fCovar);
    // param distance
    TMatrixD covar(covar0);
    covar += covar1;
    covar.Invert();
    TMatrixD delta(nPar, 1, (param1 - param0).GetMatrixArray());
    TMatrixD delta2(covar, TMatrixD::kMult, delta);
    delta.T();
    TMatrixD chi2(delta, TMatrixD::kMult, delta2);
    chi2 *= vecPar[2 * nPar][iter];
    vecPar[2 * nPar + 1][iter] = chi2(0, 0);
    Double_t logLike = chi2_00 + chi2_01 + chi2_11;
    Double_t normDistance = TMath::Sqrt(chi2(0, 0));
    //
    if (fcnP < 0) fcnP = fFCNCounter;
    Int_t fcnCounter = fFCNCounter - fcnP;
    fcnP = fFCNCounter;
    if (fStreamer) {
      (*fStreamer) << "crossValidation" <<
                   "reportName.=" << &rName <<
                   "iter=" << iter <<          // iteration
                   "fStatus=" << fStatus <<    // minuit status
                   "nPoints=" << nPoints <<    // number of points
                   "chi2_00=" << chi2_00 <<    // log likelihood for training sample
                   "chi2_01=" << chi2_01 <<    // log likelihood for test sample using training fit
                   "chi2_11=" << chi2_11 <<    // log likelihood for test sample using test fit
                   "fcnCounterI=" << fFCNCounter <<   // fcn counter integral
                   "fcnCounter=" << fcnCounter << // fcn counter per fit
                   "param0.=" << &param0 <<    // parameters estimate training sample
                   "covar0.=" << &covar0 <<    // covariance for training sample
                   "param1.=" << &param1 <<    // parameters for complementary subsample
                   "covar1.=" << &covar1 <<    // covariance for complementary subsample
                   "logLike=" << logLike << "normDistance=" << normDistance << // parameters  normalized distance
                   "\n";
    }
  }
  fPointIndex.Set(0);
}

/// \brief  Fitting with outlier removal inspired by RANSAC  - see https://en.wikipedia.org/w/index.php?title=Random_sample_consensus&oldid=767657742
/// * TODO: Make robust n-dimensional estimator
/// \param nFitPoints  - number of points used in one fit
/// \param nIter       - number of iterations
/// \param reportName  - report name (!=NULL- specified stream to save intermediate data)
/// \param option      - fit option
/// * Algorithm:
///  * 0.) Generate permutations of data
///  * 1.) Calculate fit parameters for small subset of data
///  * 2.) Calculate normalized distance between each fit vector
///  * 3.) Reject outliers and calculate robust mean
///  * 4.) Extract combined statistic
///  * 5.) Dump results into tree in verbose mode
void AliTMinuitToolkit::MISAC(Int_t nFitPoints, UInt_t nIter, const char * reportName, Option_t *option){
  TObjString rName(reportName);
  const Double_t kRelMedDistanceCut=2;            // outlier rejection cut - should be parameter
  Int_t nPoints= fPoints->GetNrows();         //
  Int_t nPar=fFormula->GetNpar();             //
  Int_t nSamples=nFitPoints*nIter;            // number of fit samples
  Int_t nRepeat=(1+(nSamples/nPoints));       //
  Int_t *permutationIndex= new Int_t[nRepeat*nPoints];   // generate random permutations
  Double_t *xrndm=new Double_t[nPoints];
  std::vector<TMatrixD*>  paramArray(nIter); // fit parameters
  std::vector<TMatrixD*>  covarArray(nIter); // fit covariance
  std::vector<double> logLike(nIter);         // fit to points
  std::vector<double> medDistance(nIter);     // fit to other fits
  std::vector<int> fitStatus(nIter);         // fit status
  if (nullptr == fRMSEstimator){
    fRMSEstimator=new TVectorD(nPar);
  }
  // 0. Generate permutations of data
  for (Int_t iR=0; iR<nRepeat; iR++){         // make nRepeat random permutations of points
    for (Int_t iPoint=0; iPoint<nPoints; iPoint++){
      xrndm[iPoint]=gRandom->Rndm();
    }
    TMath::Sort(nPoints, xrndm, &(permutationIndex[iR*nPoints]));
  }
  // 1.) Calculate fit parameters for small subset of data 
  for (UInt_t iter=0; iter<nIter; iter++){
    fPointIndex.Set(nFitPoints, &(permutationIndex[iter*nFitPoints]));
    Fit(option);
    fitStatus[iter]=fStatus;
    FitterFCN(nPar, nullptr,fChi2, fParam->GetMatrixArray(),0);
    logLike[iter]=fChi2;
    paramArray[iter] = new TMatrixD(nPar,1,fParam->GetMatrixArray());
    covarArray[iter] = new TMatrixD(*fCovar);
  }
  // 2. Calculate normalized distance between each fits
  std::vector< std::vector<double> > logDistance(nIter, std::vector<double>(nIter));  //store vectors and chi2
  for (UInt_t iter0=0; iter0<nIter; iter0++){
    for (UInt_t iter1=0; iter1<nIter; iter1++){
      if ( ((fitStatus[iter0]&kCovarFailed)==0 && (fitStatus[iter1]&kCovarFailed))==0){
        TMatrixD covar(*(covarArray[iter0]));
        covar+=*(covarArray[iter1]);
        covar.Invert();
        TMatrixD delta(*(paramArray[iter0]));
        delta-=*(paramArray[iter1]);
        TMatrixD delta2(covar,TMatrixD::kMult,delta);
        delta.T();
        TMatrixD chi2(delta,TMatrixD::kMult,delta2);
        chi2*=logLike[iter0]+logLike[iter1];
        if (chi2(0,0)>0) {
          Double_t normDistance = TMath::Sqrt(chi2(0, 0));
          logDistance[iter0][iter1] = normDistance;
        }else{
          logDistance[iter0][iter1]=-1.;  ///TODO- check failure of the fit
        };
      }else{
        logDistance[iter0][iter1]=-1.;
      }
    }
  }
  // 3. Reject outliers and calculate robust mean
  for (UInt_t iter=0; iter<nIter; iter++){
    medDistance[iter]=TMath::Median(nIter, logDistance[iter].data());
  }
  Double_t medMedDistance=TMath::Median(nIter,medDistance.data());
  // 4.  Extract combined statistic
  // 4.a 1D median,mean,rms
  // 4.b nD ?
  if (fMISACEstimators== nullptr) fMISACEstimators = new TMatrixD(4,nPar);
  for (Int_t iPar=0;iPar<nPar; iPar++){
    std::vector<double>workingArray(nIter);
    Int_t nAccepted=0;
    for (UInt_t iter=0; iter<nIter; iter++){
      if (medDistance[iter]<0) continue;
      if (medDistance[iter]/medMedDistance>kRelMedDistanceCut) continue;
      workingArray[nAccepted]=(*(paramArray[iter]))(iPar,0);
      nAccepted++;
    }
    if (nAccepted>0){
      (*fMISACEstimators)(0,iPar)=TMath::Median(nAccepted,workingArray.data());
      (*fMISACEstimators)(1,iPar)=TMath::Mean(nAccepted,workingArray.data());
      (*fMISACEstimators)(2,iPar)=TMath::RMS(nAccepted,workingArray.data());
    }
    (*fParam)(iPar)=(*fMISACEstimators)(0,iPar);
    (*fRMSEstimator)(iPar)=(*fMISACEstimators)(2,iPar);

  }
  // 5. Dump results into tree in verbose mode
  if (fStreamer){
    for (UInt_t iter=0; iter<nIter; iter++){
      (*fStreamer)<<"misac"<<
                  "iter="<<iter<<          // iteration
                  "reportName.="<<&rName<<
                  "fStatus="<<fStatus<<    // minuit status
                  "nPoints="<<nPoints<<    // number of points
                  "param.="<<paramArray[iter]<<
                  "covar.="<<covarArray[iter]<<
                  "medDistance="<<medDistance[iter]<<
                  "medMedDistance="<<medMedDistance<<
                  "misacE.="<<fMISACEstimators<<   // estimator matrix (estimator(median,mean,rms),param)
                  "\n";
    }
  }
  fPointIndex.Set(0);
  delete [] xrndm;
  delete [] permutationIndex;
}


/// \brief Register default fitters
/// * Registered function can be used in standard minimization using sting identifier, e.g:
///   * \code AliTMinuitToolkit::Fit(TH1* his, const char *fitterName, Option_t* option, Option_t* goption,Option_t* foption, Double_t xMin, Double_t xMax) \endcode
///   * \code AliTMinuitToolkit::Fit(TH1* his, const char *fitterName, Option_t* option, Option_t* goption,Option_t* foption, Double_t xMin, Double_t xMax) \endcode
/// * Registered fitter can be accessed - eventually modified (e.g to set initial values)
///   * \code AliTMinuitToolkit * pol1R= AliTMinuitToolkit::GetPredefinedFitter("pol1R");\endocode
/// * Default fitters:
///   * formula: pol<N>, gauss
///   * log-likelihood: chi2, likeGausCachy, huber norm
void  AliTMinuitToolkit::RegisterDefaultFitters(){
  TF1 *likeGausCachy = new TF1("likeGausCachy", AliTMinuitToolkit::GaussCachyLogLike,-10,10,2);
  TF1 *likePseudoHuber = new TF1("likePseudoHuber", AliTMinuitToolkit::PseudoHuberLogLike,-10,10,1);
  likeGausCachy->SetParameters(0.8,1);
  likePseudoHuber->SetParameter(0,3);
  likeGausCachy->GetXaxis()->SetTitle("#Delta");
  likeGausCachy->GetYaxis()->SetTitle("-logLikelihood");
  likeGausCachy->SetLineColor(2);
  likePseudoHuber->GetXaxis()->SetTitle("#Delta");
  likePseudoHuber->GetYaxis()->SetTitle("-logLikelihood");
  likePseudoHuber->SetLineColor(4);
  //
  for (Int_t iPol=0; iPol<10; iPol++){ //register polynomial fitters
    TMatrixD initPar(iPol+1,4);
    for (Int_t iPar=0; iPar<iPol+1; iPar++) initPar(iPar,1)=1;
    //
    TF1 *fpol = new TF1(TString::Format("fpol%d",iPol).Data(),TString::Format("pol%d",iPol).Data());
    AliTMinuitToolkit * fitter1D = new AliTMinuitToolkit();
    fitter1D->SetVerbose(0x1); fitter1D->SetFitFunction(fpol,kTRUE);
    fitter1D->SetInitialParam(&initPar);
    fitter1D->SetName(TString::Format("pol%d",iPol).Data());
    AliTMinuitToolkit::SetPredefinedFitter(fitter1D->GetName(),fitter1D);
    // gaus log likelihood cost function
    AliTMinuitToolkit * fitter1DR = new AliTMinuitToolkit();
    fitter1DR->SetVerbose(0x1); fitter1DR->SetFitFunction(fpol,kTRUE);
    fitter1DR->SetLogLikelihoodFunction(likeGausCachy);
    fitter1DR->SetInitialParam(&initPar);
    fitter1DR->SetName(TString::Format("pol%dR",iPol).Data());
    AliTMinuitToolkit::SetPredefinedFitter(fitter1DR->GetName(),fitter1DR);
    // huber cost function
    AliTMinuitToolkit * fitter1DH = new AliTMinuitToolkit();
    fitter1DH->SetVerbose(0x1); fitter1DH->SetFitFunction(fpol,kTRUE);
    fitter1DH->SetInitialParam(&initPar);
    //fitter1DH->EnableRobust(true);
    fitter1DH->SetLogLikelihoodFunction(likePseudoHuber);
    fitter1DH->SetName(TString::Format("pol%dH",iPol).Data());
    AliTMinuitToolkit::SetPredefinedFitter(fitter1DH->GetName(),fitter1DH);
  }
  //
  TF1 *fGaus = new TF1("fgaus","gaus");
  TMatrixD initPar(3,4);
  initPar(0,0)=0; initPar(0,1)=1;
  initPar(1,0)=0; initPar(1,1)=1;
  initPar(2,0)=1; initPar(2,1)=1;
  AliTMinuitToolkit * fitterGR = new AliTMinuitToolkit();
  fitterGR->SetVerbose(0x1); fitterGR->SetFitFunction(fGaus,kTRUE);
  fitterGR->SetLogLikelihoodFunction(likeGausCachy);
  fitterGR->SetInitialParam(&initPar);
  AliTMinuitToolkit::SetPredefinedFitter("gausR",fitterGR);
  //
  //
  // TFormulaPrimitive::AddFormula(new TFormulaPrimitive("GausKS","GausKS",(TFormulaPrimitive::GenFuncG)AliTMinuitToolkit::GausKurtosisSkewness,5));
  // hack number of arguments in the primitive (protected member) - //TODO - fix ROOT
  // TFormulaPrimitive * p = TFormulaPrimitive::FindFormula("GausKS");
  // TDataMember *m = (TDataMember*)p->IsA()->GetListOfDataMembers()->FindObject("fNArguments");
  // *((Int_t*)(((char*)p)+m->GetOffset()))=1;
}

/// \breief Register plane fitters (As in the TLinearFitter)
/// Global registered fitters can be obtained e.g :
/// \code AliTMinuitToolkit * hyp2R= AliTMinuitToolkit::GetPredefinedFitter("hyp2R"); \endcode
/// \param nPlanes       - number of planes
/// \param logLikeType   - type of the log likelihood
///                        * 0) chi2, 1.) gaus+cauchy, 2.) pseudo-huber
/// \return return registered fitter
AliTMinuitToolkit * AliTMinuitToolkit::RegisterPlaneFitter(Int_t nPlanes, Int_t logLikeType){
  AliTMinuitToolkit *fitter = new AliTMinuitToolkit(TString::Format("AliTMinuitToolkitTest%d.root", nPlanes));
  TString sFormula = "[0]*x[0]";
  for (Int_t iPlane=1; iPlane<nPlanes; iPlane++) sFormula+=TString::Format("+[%d]*x[%d]", iPlane,iPlane);
  TFormula *formula = new TFormula(TString::Format("hyp%dR",nPlanes).Data(), sFormula.Data());
  fitter->SetFitFunction((TF1 *) formula, kTRUE);
  //
  if (logLikeType==0) { // standard minimization
    TF1 *likeGaus = new TF1("likeGausCachy", "x*x*0.5", -10, 10);
    fitter->SetLogLikelihoodFunction(likeGaus);
    fitter->SetName(TString::Format("hyp%d", nPlanes).Data());
  }
  if (logLikeType==1){
    TF1 *likeGausCachy = new TF1("likeGausCachy", AliTMinuitToolkit::GaussCachyLogLike, -10, 10, 2);
    likeGausCachy->SetParameters(0.8, 1);
    fitter->SetLogLikelihoodFunction(likeGausCachy);
    fitter->SetName(TString::Format("hyp%dR", nPlanes).Data());
  }
  if (logLikeType==2){
    TF1 *likePseudoHuber = new TF1("likePseudoHuber", AliTMinuitToolkit::PseudoHuberLogLike, -10, 10, 1);
    likePseudoHuber->SetParameter(0,3);
    fitter->SetLogLikelihoodFunction(likePseudoHuber);
    fitter->SetName(TString::Format("hyp%dH", nPlanes).Data());
  }
  TMatrixD *initParam=new TMatrixD(nPlanes,4);
  for (Int_t iPar=0; iPar<nPlanes; iPar++) (*initParam)(iPar,1)=1;
  fitter->SetInitialParam(initParam);
  AliTMinuitToolkit::SetPredefinedFitter(fitter->GetName(), fitter);
  return fitter;
}


/// \brief Make histogram fit
/// * TODO- add example code snippet - using the code from test directory
/// \param his              - pointer to histogram
/// \param fitterName       - fitter name ()
/// \param fitOption        - fit Option
/// \param hisOption        - graphic option for histogram
/// \param funOption        - graphic option for function drawing
/// \param xMin             - fit range min
/// \param xMax             - fit range max
/// \return                 - pointer to the fitter

AliTMinuitToolkit *  AliTMinuitToolkit::Fit(TH1* his, const char *fitterName, Option_t* fitOption, Option_t* hisOption,Option_t* funOption, Double_t xMin, Double_t xMax){
  AliTMinuitToolkit *fitter=GetPredefinedFitter(fitterName);
  if (fitter== nullptr){
    ::Error("AliTMinuitToolkit::Fit","Non supported fitter %s. \nfitter can be added to the list using. SetPredefinedFitter", fitterName);
    return nullptr;
  }
  if (his== nullptr){
    ::Error("AliTMinuitToolkit::Fit","NULL pointer");
    return nullptr;
  }
  fitter->FitHistogram(his,fitOption);
  TF1 * fitFun = (TF1*)(fitter->GetFormula()->Clone());
  fitFun->SetParameters(fitter->GetParameters()->GetMatrixArray());
  SetFunctionDrawOption(fitFun,funOption);
  his->GetListOfFunctions()->AddLast(fitFun);
  if (xMin<xMax) {
    fitFun->SetRange(xMin,xMax);
  }else{
    fitFun->SetRange(his->GetXaxis()->GetXmin(), his->GetXaxis()->GetXmax());
  }
  if (hisOption) his->Draw(hisOption);
  return fitter;
}


/// \brief Make graph fit
/// * TODO- add example code snippet - using the code from test directory
/// \param graph            - pointer to graph
/// \param fitterName       - fitter name ()
/// \param fitOption        - fitter setting option
/// \param graphOption      - graphic option
/// \param funOption        - graphic option for function
/// \param xMin             - fit range min
/// \param xMax             - fit range max
/// \return                 - pointer to the fitter

AliTMinuitToolkit *  AliTMinuitToolkit::Fit(TGraph *graph, const char *fitterName, Option_t* fitOption, Option_t* graphOption,  Option_t* funOption, Double_t xMin, Double_t xMax){
  //
  //
  AliTMinuitToolkit *fitter=GetPredefinedFitter(fitterName);
  if (fitter== nullptr){
    ::Error("AliTMinuitToolkit::Fit","Non supported fitter %s. \nfitter can be added to the list using. SetPredefinedFitter", fitterName);
    return nullptr;
  }
  if (graph== nullptr){
    ::Error("AliTMinuitToolkit::Fit","NULL pointer");
    return nullptr;
  }
  fitter->FitGraph(graph,fitOption);
  TF1 * fitFun = (TF1*)(fitter->GetFormula()->Clone());
  fitFun->SetParameters(fitter->GetParameters()->GetMatrixArray());
  fitFun->SetName(TString::Format("%s:%s",fitter->GetName(), fitOption).Data());
  fitFun->SetTitle(TString::Format("%s:%s",fitter->GetName(), fitOption).Data());
  SetFunctionDrawOption(fitFun,funOption);
  graph->GetListOfFunctions()->AddLast(fitFun);
  if (xMin<xMax) {
    fitFun->SetRange(xMin,xMax);
  }else{
    fitFun->SetRange(graph->GetXaxis()->GetXmin(), graph->GetXaxis()->GetXmax());
  }
  if (graphOption) graph->Draw(graphOption);
  return fitter;
}

/// Parse draw option for function
/// \param fun        - pointer to function to change
/// \param option     - draw option for function
/// * TODO- add example using the code from test directory
/// * Example usage:
/// \code AliTMinuitToolkit::SetFunctionDrawOption(fpol1,"funOption(1,1,1)")\endcode

void  AliTMinuitToolkit::SetFunctionDrawOption(TF1 *fun, Option_t *option){
  TString funOption=option; // e.g  "funOption(1,1,1)"
  TPRegexp regFunOption("funOption\\(");
  // Error message
  if (regFunOption.Match(funOption)){
    Int_t index0=funOption.Index(regFunOption)+9;
    Int_t index1=funOption.Index(")",index0+1);
    if (index1<index0) {
      ::Error("AliTMinuitToolkit::SetFunctionDrawOption","Invalid function draw option syntax %s", option);
      return;
    }
    Int_t index=index0+1;
    Int_t color=TString(&(funOption.Data()[index])).Atoi();
    if (color>0) fun->SetLineColor(static_cast<Color_t>(color));
    index=funOption.Index(",",index+1)+1;
    if (index>0) {
      Int_t width=TString(&(funOption.Data()[index])).Atoi();
      if (width>0) fun->SetLineWidth(static_cast<Width_t>(width));
      index=funOption.Index(",",index+1)+1;
      if (index>0) {
        Int_t style=TString(&(funOption.Data()[index])).Atoi();
        if (style>0) fun->SetLineStyle(static_cast<Style_t>(style));
      }
    }
  }
}
