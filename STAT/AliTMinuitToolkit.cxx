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
#include <TFormula.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TMatrix.h>
#include <TRandom.h>
#include <TString.h>
#include <TVector.h>
#include <TVectorD.h>
#include <TVirtualFitter.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTree.h>
#include <TRandom.h>
#include <TTreeStream.h>
#include "AliSysInfo.h"
#include "TFitter.h"
#include "TMinuit.h"
#include "AliTMinuitToolkit.h"
#include "TGraph.h"
#include <vector>
#include "TPRegexp.h"
#include "TStatToolkit.h"

std::map<std::string, AliTMinuitToolkit*> AliTMinuitToolkit::fPredefinedFitters;  

//--------------------------------------------------------------------------------------
//
// The AliTMinuitToolkit serves as an easy to use interface for the TMinuit 
// package: 

// A small example illustrating the usage of AliTMinuitToolkit is given in the function 
// "AliTMinuitToolkit::Test()" , respectively in  $ALICE_ROOT/../src/STAT/test/AliTMinuitToolkitTest.C+
// 
// 
// 1. Setting the formula:
//  The formula is simply set via "void SetFitFunction(TFormula * formula, Bool_t doReset)".
//
// 2. Adding the data points
//
//  - Direct specification of the points is possible via
//    "void SetPoints(TMatrixD * points)". 
//     Note, that the each point corresponds to one row in the matrix. The fitter is then started with  the command "void Fit()". 
//  - In order to fit a histogram, use "void FitHistogram(TH1F * his)"
//  - Filling input points useing TTree queries 
//     tool1D->FillFitter(inputTree,"test1D+noise:1/sqrt(12.+0)","testx0", "", 0,fitEntries);
//
// 3. Accessing the fit results
//
//   The N parameters of the formula are stored in a N-dimensional vector which
//   is returned by "TVectorD * GetParameters()". In a similar way the covariance 
//   matrix of the fit is returned via "TMatrixD * GetCovarianceMatrix()" which
//   is of the type N x N.
//
//
// 4. Non-linear robust fitting - see example $ALICE_ROOT/../src/STAT/test/AliTMinuitToolkitTest.C:Test1D()
//
//  Even a few outliers can lead to wrong results of a least-squares fitting 
//  procedure. In this case the use of robust(resistant) methods can be 
//  helpful, but a stronger dependence on starting values or convergence to
//  local minima can occur.
//
//  a.) By default chi2 minimization is invoked. 
//  b.) Optionally more robust log likelihood- Huber cost function can be specified instead of chi^2.
//   For small deviations  the function behaves like x^2 and for larger deviations like |x| - the so 
//   called Huber estimator:
//
//   h(x) = x^2                              , for x < 2.5*sigma
//   h(x) = 2*(2.5*sigma)*x - (2.5*sigma)^2  , for x > 2.5*sigma
//      REFERENCE for non-linear robust fitting:
//      Ekblom H. and Madsen K. (1988), Alogrithms for non-linear Huber estimation,
//      BIT Numerical Mathematics 29 (1989) 60-76.
//      internet: http://www.springerlink.com/content/m277218542988344/
//  c.) User defined log likelihood function can be specified 
//        TF1 * fcost = new TF1("1","abs(x)<10?-log(0.8*exp(-x**2)+0.2/(1+x**2)):-log(0.2/(1+x**2))",-20,20); // 80 % gaus + 20% cachy
//        tool1D->SetLogLikelihoodFunction(fcost);
//        tool1D->Fit();
//  
//  5.) Work in progress
//         *  Multidimensional fits (more than one observable) 
//         **    e.g fitting Er and Erphi in parallel
//         *  Regularization
//         ** e.g positive definite, monotonous function            
//         *  CrossValidation
//
//
// Comments and questions are always welcome: marian.ivanov@cern.ch, A.Kalweit@gsi.de
//
//--------------------------------------------------------------------------------------


ClassImp(AliTMinuitToolkit)

AliTMinuitToolkit::AliTMinuitToolkit(const char *streamerName, const char *mode) :
  TNamed(),
  fStreamer(0),
  fVerbose(0),
  fFormula(0),
  fLogLikelihoodFunction(0),
  fFitAlgorithm(""),
  fPoints(0),
  fValues(0),
  fVarNames(0),           // variable names  
  fValueNames(0),         // value  names                    
  fPointIndex(0),              // point index - to specify points for likelihood calculation (used for Bottstrap and CrossValidation)             
  fParam(0),
  fInitialParam(0),
  fRMSEstimator(0),       // parameter spread as estimated in Bootstrap resp. TwoFold cross validation
  fMISACEstimators(0),    // MISAC estimators - median, mean, rms
  fCovar(0),
  fChi2(0),
  fMaxCalls(0),
  fPrecision(0),
  fIsHuberCost(0)
{
 //
 // standard constructor
 //
 fMaxCalls = 500;
 fPrecision = 1;
 fIsHuberCost = false;
 if (streamerName!=NULL) fStreamer= new TTreeSRedirector(streamerName,mode);
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
  if (fStreamer) delete fStreamer;
}

void  AliTMinuitToolkit::SetStreamer(const char *streamerName, const char *mode){
  //
  // set streamer
  if (fStreamer) delete fStreamer;
  fStreamer= new TTreeSRedirector(streamerName,mode);
}


void  AliTMinuitToolkit::ClearData(){
  delete fPoints;
  fPoints=0;
  delete fValues;
  fValues=0;
  //  delete fParam;
  //   fParam=0;
  //   delete fCovar;
  //   fCovar=0;
}

void  AliTMinuitToolkit::SetFitFunction(TF1 *const formula, Bool_t doReset) {
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


void   AliTMinuitToolkit::SetInitialParam(TMatrixD *const paramLimits) { 
  fInitialParam=new TMatrixD(*paramLimits);
};  



void AliTMinuitToolkit::FitHistogram(TH1 *const his, Option_t *option) {
  //
  // Fit a one dimensional histogram
  //
  ClearData();
  fPoints  = new TMatrixD(his->GetNbinsX(), 1);
  fValues  = new TMatrixD(his->GetNbinsX(), 2); 
  for(Int_t ibin=0; ibin < his->GetNbinsX(); ibin++) {
    Double_t x = his->GetXaxis()->GetBinCenter(ibin+1);
    Double_t y = his->GetBinContent(ibin+1);  
    (*fPoints)(ibin, 0) = x;
    Double_t err=his->GetBinError(ibin+1);
    (*fValues)(ibin,0)=y;
    (*fValues)(ibin,1)=(err>0)? 1./err:0;
  }
  Fit(option);
}

void AliTMinuitToolkit::FitGraph(TGraph *const gr, Option_t *option) {
  //
  // Fit a one dimensional histogram
  //
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








void AliTMinuitToolkit::FitterFCN(int &/*npar*/, double */*info*/, double &fchisq, double *gin, int iflag){
  //
  // internal function which gives the specified function to the TMinuit function
  //  
  //
  static Int_t fcnCounter=0;
  //if (info) info[0]=fcnCounter;
  fcnCounter++;
  AliTMinuitToolkit * fitter = (AliTMinuitToolkit*)TVirtualFitter::GetFitter()->GetObjectFit();
  TF1 *logLike=fitter->fLogLikelihoodFunction;
  const Double_t* likeParam= (logLike!=NULL) ? logLike->GetParameters():NULL;
  fchisq = 0;
  const TMatrixD & variables= (*fitter->GetPoints());
  const TMatrixD & values=    (*fitter->GetValues());
  Int_t nvars       = variables.GetNcols();
  Int_t npoints     = variables.GetNrows();
  TVectorD param(fitter->fFormula->GetNpar(), gin);
  // calculate  log likelihood
  //
  if (fitter->fPointIndex.GetSize()>0) { 
    npoints=fitter->fPointIndex.GetSize();
  }
  for (Int_t jpoint=0; jpoint<npoints; jpoint++){
    Int_t ipoint=(fitter->fPointIndex.GetSize()>0) ? fitter->fPointIndex[jpoint]:jpoint;
    Double_t x[100];
    for (Int_t ivar=0; ivar<nvars; ivar++){
      x[ivar] = variables(ipoint, ivar);      
    }    
    Float_t funx = (fitter->GetFormula())->EvalPar(x,gin);   
    Double_t value=values(ipoint,0);
    Double_t weight=values(ipoint,1);
    Double_t delta = TMath::Abs(value - funx);
    Double_t toAdd=0;
    if (logLike){
      Double_t normDelta = delta*weight;
      toAdd=logLike->EvalPar(&delta,likeParam);
      //      continue;
    }else{
      if (fitter->IsHuberCost() == true) {       //hubert norm
	delta = delta*weight;                   // normalization
	if (delta <= 2.5) toAdd= delta*delta; // new metric: Huber-k-estimator
	if (delta > 2.5)  toAdd= 2*(2.5)*delta - (2.5*2.5);
      } else {
	Double_t chi2 = delta*weight;
	chi2*=chi2;
	toAdd=chi2;   // chi2 (log likelihood)
      }
    }
    //if (fitter-
    if (fitter->IsVerbose(kStreamFcnPoint) || (iflag&kStreamFcnPoint)){
      TVectorD vecX(nvars, x);
      (*(fitter->fStreamer))<<"fcnDebugPoint"<<
	"toAdd="<<toAdd<<  // log likelihood to add
	"val="<<value<<    // "measurement"
	"w="<<weight<<     // wight of point
	"fun="<<funx<<     // function value
	"x.="<<&vecX<<     // variable vector
	"p.="<<&param<<    // parameter vector
	"\n";
    }
    fchisq+=toAdd;    
  }
  if (fitter->IsVerbose(kStreamFcn) || (iflag&kStreamFcn )){
    (*(fitter->fStreamer))<<"fcnDebug"<<
      "fchisq="<<fchisq<<
      "p.="<<&param<<
      "\n";
  }
  fcnCounter++;
}
 

void AliTMinuitToolkit::Fit(Option_t *option) {
  //
  // internal function that calls the fitter
  //
  TString sOption=(option==NULL)?"":option;
  TPRegexp regBootstrap("bootstrap[0-9]*");
  TPRegexp regTwofold("twofold[0-9]*");
  TPRegexp regMISAC("misac\\([0-9]*,[0-9]*\\)");
  TPRegexp regFitName("rFitName:");

  // run bootstap if specified
  if (regBootstrap.Match(sOption)>0){      // run boootstrap
    TString newOption=sOption;
    regBootstrap.Substitute(newOption,"");
    Int_t nPoints=TString(&(sOption.Data()[sOption.Index(regBootstrap)+9])).Atoi();
    if (nPoints<0) nPoints=5;
    return Bootstrap(nPoints,0,newOption.Data());
  }
  // run two fold minimization if specified
  if (regTwofold.Match(sOption)>0){      // run twofold minimization
    TString newOption=sOption;
    regTwofold.Substitute(newOption,"");
    Int_t nPoints=TString(&(sOption.Data()[sOption.Index(regTwofold)+7])).Atoi();
    if (nPoints<0) nPoints=5;
    return TwoFoldCrossValidation(nPoints,0,newOption.Data());
  }
  // run MISAC
  if (regMISAC.Match(sOption)>0){      // run MISAC minimization  - mix of RAMSAC and KFold 
    TString newOption=sOption;
    regMISAC.Substitute(newOption,"");
    Int_t index1=sOption.Index(regMISAC)+6;
    Int_t index2=sOption.Index(",",index1+1)+1;
    Int_t nPoints=TString(&(sOption.Data()[index1])).Atoi();
    Int_t nIter=TString(&(sOption.Data()[index2])).Atoi();
    if (nPoints<0) nPoints=5;
    return MISAC(nPoints,nIter,0,newOption.Data());
  }
  fStatus=0;
  
  //  Int_t nParam  = fParam->GetNrows();
  Int_t nParam  = fFormula->GetNpar();
  Int_t npoints = fPoints->GetNrows();
  Int_t nvar    = fPoints->GetNcols()-1;
  
  // create "intial param" in case not set before
  if (fParam==NULL) fParam=new TVectorD(nParam);
  if (fInitialParam == 0) {
    ::Info("AliTMinuitToolkit::Fit","Default initial parameters set");
    fInitialParam=  new TMatrixD(nParam,4);
    for (Int_t iparam=0; iparam<nParam; iparam++) {
      (*fInitialParam)(iparam,0)=0;
      (*fInitialParam)(iparam,1)=1000;
      (*fInitialParam)(iparam,2)=0;
      (*fInitialParam)(iparam,3)=0;
    }
  }
    
  // migrad fit algorithm as default
  if (fFitAlgorithm == "") {
    fFitAlgorithm = "migrad";
  }
  
  // assign weights
  if (fValues == 0) {
    ::Error("AliTMinuitToolkit::Fit()","Invalid fit. Values not set");
    return ;
  }
  
  // set up the fitter
  TVirtualFitter *minuit = TVirtualFitter::Fitter(0, nParam);
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
  for (Int_t iparam=0; iparam<nParam; iparam++){
    if (sOption.Contains("rndmInit")){
      Double_t initValue=(*fInitialParam)(iparam, 0)+gRandom->Gaus()*(*fInitialParam)(iparam,1);
      minuit->SetParameter(iparam, Form("p[%d]",iparam), initValue, (*fInitialParam)(iparam,1), (*fInitialParam)(iparam, 2), (*fInitialParam)(iparam, 3));
    }else{
      minuit->SetParameter(iparam, Form("p[%d]",iparam), (*fInitialParam)(iparam,0), (*fInitialParam)(iparam,1), (*fInitialParam)(iparam, 2), (*fInitialParam)(iparam, 3));
    }
    /*
      if (doReset){
      minuit->SetParameter(iparam, Form("p[%d]",iparam), (*fInitialParam)(iparam,0), (*fInitialParam)(iparam,1), (*fInitialParam)(iparam, 2), (*fInitialParam)(iparam, 3));
      }else{
      minuit->SetParameter(iparam, Form("p[%d]",iparam), (*fParam)(iparam), TMath::Sqrt((*fCovar)(iparam,iparam)), (*fInitialParam)(iparam, 2), (*fInitialParam)(iparam, 3));
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
  //  fStatus = (((TFitter*)minuit)->GetMinuit())->GetStatus(); // we should add separate status for Minuit and  AliTMinuitTookit
  
  // fill parameter vector
  for (Int_t ivar=0; ivar<nParam; ivar++){
    (*fParam)(ivar) = minuit->GetParameter(ivar);
    fFormula->SetParameter(ivar, minuit->GetParameter(ivar));
    //fFormula->SetPar(ivar, minuit->GetParameter(ivar));
  }
  FitterFCN(nParam,0,fChi2, fParam->GetMatrixArray(),0); 
  
  // fill parameter vector
  for (Int_t ivar=0; ivar<nParam; ivar++){
   (*fParam)(ivar) = minuit->GetParameter(ivar);
   fFormula->SetParameter(ivar, minuit->GetParameter(ivar));
  }
  
  // fill covariance matrix
  if (fCovar==NULL)   fCovar = new TMatrixD(nParam, nParam);
  if (minuit->GetCovarianceMatrix()){
    fCovar->SetMatrixArray(minuit->GetCovarianceMatrix());
  }else{
    fStatus|=kCovarFailed;
  }
  
}


Int_t AliTMinuitToolkit::FillFitter(TTree * inputTree, TString values, TString variables, TString selection, Int_t firstEntry, Int_t nentries, Bool_t doReset ){
  //
  // Make unbinned fit
  //
  TString query=values;
  query+=":";
  query+=variables;
  if (inputTree==NULL){
    ::Error("AliTMinuitToolkit::UnbinnedFit","Zerro input tree");
    return -1;
  }
  fValueNames=values.Tokenize(":");
  fVarNames=variables.Tokenize(":");
  Int_t nVals= fValueNames->GetEntries();
  Int_t nVars= fVarNames->GetEntries();
  if (doReset==kFALSE  && fPoints!=NULL){
    if (fPoints->GetNrows()!=nVars){
      ::Error("AliTMinuitToolkit::UnbinnedFit","Not comatible number of variables: %d instead of %d: variables:  %s",nVars, fPoints->GetNrows(), query.Data());
      return -1;
    }
  }

  Int_t entries = inputTree->Draw(query.Data(),selection.Data(),"goffpara",nentries,firstEntry);
  if (entries<=0) {
    ::Error("AliTMinuitToolkit::UnbinnedFit","badly formatted values or variables: %s",query.Data());
    ::Error("AliTMinuitToolkit::UnbinnedFit","valueDescription: %s",values.Data());
    ::Error("AliTMinuitToolkit::UnbinnedFit","variables: %s",variables.Data());
    return -1;
  }  
  Int_t index0=0;
  if (doReset || fPoints==NULL) {
    ClearData();
    fPoints=new TMatrixD(entries,nVars);
    fValues=new TMatrixD(entries,nVals);
  }else{
    index0= fPoints->GetNrows();
    fPoints->ResizeTo(index0+entries,nVars);
    fValues->ResizeTo(index0+entries,nVals);
  }
  for (Int_t iPoint=0; iPoint<entries; iPoint++){
    for (Int_t iVar=0; iVar<nVars; iVar++){
      (*fPoints)(index0+iPoint,iVar)=inputTree->GetVal(iVar+nVals)[iPoint];
    }
    for (Int_t iVal=0; iVal<nVals; iVal++){
      (*fValues)(index0+iPoint,iVal)=inputTree->GetVal(iVal)[iPoint];
    }
  }
  return entries;
}

/// forrma alias string for the for tree alias or for fit title
/// \param option  - use latex in case of title request
/// \return        - string description
TString AliTMinuitToolkit::GetFitFunctionAsAlias(Option_t *option, TTree * tree){
  //
  // construct string TTree alias for fit function
  TString inputString(fFormula->GetTitle());
  TString sOption(option);
  sOption.ToLower();

  if (fVarNames==NULL){
    ::Error("AliTMinuitToolkit::GetFitFunctionAsAlias","Variable names not defined. Fucntion supported for TTree::UnbinnedFit");
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


Double_t AliTMinuitToolkit::RrndmGaus(Double_t mean, Double_t sigma){
  //
  return gRandom->Gaus(mean,sigma);
}

Double_t  AliTMinuitToolkit::RrndmLandau(Double_t mean, Double_t sigma){
  return gRandom->Landau(mean,sigma);
}

Double_t  AliTMinuitToolkit::GaussCachyLogLike(Double_t *x, Double_t *p){
  //
  // parameters
  Double_t vCauchy=p[1]/(TMath::Pi()*(p[1]*p[1]+(*x)*(*x)));
  Double_t vGaus=(TMath::Abs(*x)<20) ? TMath::Gaus(*x,0,1.,kTRUE):0.;
  return TMath::Abs(TMath::Log(p[0]*vGaus+(1-p[0])*vCauchy)-1);
  static Double_t norm= 1/TMath::Gaus(0,0,1.,kTRUE);
  return TMath::Abs(TMath::Log((p[0]*vGaus+(1-p[0])*vCauchy)*norm));
}

void AliTMinuitToolkit::Bootstrap(Int_t nIter, const char * reportName, Option_t *option){
  //
  // Bootstrap minimization 
  // fitting parameters done on several times on modified data (random samples with replacement)  (to be done togerher with M-stimator)
  //    to emulate different data population
  //    reduce problem with local minima in the M-estimator estimated parameters - robust mean/median
  //    obtained RMS can be used to estimate error of the estemator
  //   
  static Int_t counter=0;
  Int_t nPoints= fPoints->GetNrows();
  fPointIndex.Set(nPoints);
  Double_t info[3]={0,0,0};
  Int_t nPar=fFormula->GetNpar();
  Int_t fcnP=-1;
  TObjString rName("bootstrap");
  if (reportName!=NULL) rName=reportName;
  std::vector< std::vector<double> > vecPar(nPar, std::vector<double>(nIter));
  for (Int_t iter=0; iter<nIter; iter++){
    counter++;
    for (Int_t iPoint=0; iPoint<nPoints; iPoint++){
      fPointIndex[iPoint]=gRandom->Rndm()*nPoints;
    }
    Fit(option);
    FitterFCN(nPar,NULL,fChi2, fParam->GetMatrixArray(),0); 
    if (fcnP<0) fcnP=info[0];
    if (fVerbose&kSysInfo) AliSysInfo::AddStamp("Bootstrap",0,counter,iter,info[0]-fcnP);
    //Int_t fcnCounter=info[0]-fcnP;
    //fcnP=info[0];
    for (Int_t iPar=0; iPar<nPar;iPar++)  vecPar[iPar][iter]=(*fParam)(iPar);
    if (fStreamer){      
      (*fStreamer)<<"bootstrap"<<
	"iter="<<iter<<
	"fStatus="<<fStatus<<    // minuit status
	"reportName.="<<&rName<<
	"nPoints="<<nPoints<<   
	"fChi2="<<fChi2<<
	//	"fcnCounterI="<<info[0]<<
	//	"fcnCounter="<<fcnCounter<<
	"fParam.="<<fParam<<
	"fCovar.="<<fCovar<<
	"\n";
    }
  }
  if (fRMSEstimator==NULL) fRMSEstimator=new TVectorD(*fParam);
  TVectorD &par=*fParam;
  TVectorD &rms=*fRMSEstimator;
  for (Int_t iPar=0; iPar<nPar;iPar++) {
    Double_t rmsf,meanf;
    TStatToolkit::EvaluateUni(nIter,vecPar[iPar].data(),meanf,rmsf,TMath::Max(nIter*0.75,1.));
    par[iPar]=TMath::Median(nIter,vecPar[iPar].data());
    rms[iPar]=rmsf;
  }
  fPointIndex.Set(0);
  return ;
}


void AliTMinuitToolkit::TwoFoldCrossValidation(Int_t nIter, const char * reportName, Option_t *option){
  //
  // Two fold cross validation 
  // Sample randomly divided to the independent training set and  complementary test set
  //  - log likehood cost function calculated for training fit set 
  //     
  static Int_t counter=0;
  Int_t nPoints= fPoints->GetNrows();
  TArrayI indexFold(nPoints);
  TArrayF rndmCV(nPoints);  
  fPointIndex.Set(nPoints);
  Double_t chi2_00, chi2_01,  chi2_10, chi2_11;
  Double_t info[3]={0,0,0};
  Int_t nPar=fFormula->GetNpar();
  Int_t fcnP=-1;
  TObjString rName(reportName);
  std::vector< std::vector<double> > vecPar(2*nPar+4, std::vector<double>(nIter));  //store vectors and chi2
  for (Int_t iter=0; iter<nIter; iter++){
    for (Int_t iPoint=0; iPoint<nPoints; iPoint++){
      rndmCV[iPoint]=gRandom->Rndm()*nPoints;
    }
    TMath::Sort(nPoints,rndmCV.GetArray(), indexFold.GetArray());
    //
    fPointIndex.Set(nPoints/2,indexFold.GetArray());
    Fit(option);
    for (Int_t iPar=0; iPar<nPar;iPar++)  vecPar[iPar][iter]=(*fParam)(iPar);
    TVectorD param0(*fParam);
    TMatrixD covar0(*fCovar);
    FitterFCN(nPar,info,chi2_00, fParam->GetMatrixArray(),0); 
    //
    fPointIndex.Set(nPoints-nPoints/2,&(indexFold.GetArray()[nPoints/2]));
    FitterFCN(nPar,info,chi2_01, fParam->GetMatrixArray(),0); 
    Fit(option);
    for (Int_t iPar=0; iPar<nPar;iPar++)  vecPar[nPar+iPar][iter]=(*fParam)(iPar);    
    FitterFCN(nPar,info,chi2_11, fParam->GetMatrixArray(),0); 
    vecPar[2*nPar][iter]=chi2_00+chi2_01+chi2_11;
    TVectorD param1(*fParam);
    TMatrixD covar1(*fCovar);
    // param distance
    TMatrixD covar(covar0);
    covar+=covar1;
    covar.Invert();
    TMatrixD delta(nPar,1,(param1-param0).GetMatrixArray());
    TMatrixD delta2(covar,TMatrixD::kMult,delta);
    delta.T();
    TMatrixD chi2(delta,TMatrixD::kMult,delta2);
    chi2*=vecPar[2*nPar][iter];
    vecPar[2*nPar+1][iter]=chi2(0,0);
    Double_t logLike=chi2_00+chi2_01+chi2_11;
    Double_t normDistance=TMath::Sqrt(chi2(0,0));
    //
    if (fcnP<0) fcnP=info[0];
    if (fVerbose&kSysInfo) AliSysInfo::AddStamp("TwoFoldCrossValidation",0,counter,iter,info[0]-fcnP);
    counter++;
    Int_t fcnCounter=info[0]-fcnP;
    fcnP=info[0];
    if (fStreamer){      
      (*fStreamer)<<"crossValidation"<<
	"reportName.="<<&rName<<
	"iter="<<iter<<          // iteration 
	"fStatus="<<fStatus<<    // minuit status
	"nPoints="<<nPoints<<    // number of points
	"chi2_00="<<chi2_00<<    // log likelihodd for training sample
	"chi2_01="<<chi2_01<<    // log likelihood for test sample using traning fit
	"chi2_11="<<chi2_11<<    // log likelihood for test sample using test fit
	//
	"fcnCounterI="<<info[0]<<   // fcn counter integral
	"fcnCounter="<<fcnCounter<< // fcn counter per fit
	"param0.="<<&param0<<    // parameters estimate training sample
	"covar0.="<<&covar0<<    // covariance for training sample
	"param1.="<<&param1<<    // parameters for complementaray subsample
	"covar1.="<<&covar1<<    // covariance for complementaray subsample
	"logLike="<<logLike<<
	"normDistance="<<normDistance<< // parameters  normalizaed distance
	"\n";
    }
  }
  fPointIndex.Set(0);
  return ;
}

void AliTMinuitToolkit::MISAC(Int_t nFitPoints, Int_t nIter, const char * reportName, Option_t *option){
  //
  // Fitting with outlier removal inspired by RAMSAC  - see https://en.wikipedia.org/w/index.php?title=Random_sample_consensus&oldid=767657742
  //
  // nPoints - minimal number of points
  const Double_t kRelMedDistanceCut=2;            // otlier rejection cut - should be parameter
  //
  Int_t nPoints= fPoints->GetNrows();         //
  Int_t nPar=fFormula->GetNpar();             // 
  Int_t nSamples=nFitPoints*nIter;            // number of fit samples 
  Int_t nRepeat=(1+(nSamples/nPoints));       //
  Int_t *permutationIndex= new Int_t[nRepeat*nPoints];   // generate random permuations
  Double_t *xrndm=new Double_t[nPoints];
  std::vector<TMatrixD*>  paramArray(nIter); // fit parameters
  std::vector<TMatrixD*>  covarArray(nIter); // fit covariance
  std::vector<float> logLike(nIter);         // fit to points
  std::vector<float> medDistance(nIter);     // fit to other fits
  std::vector<int> fitStatus(nIter);         // fit status
  //
  // 0. Generate permuations of data
  for (Int_t iR=0; iR<nRepeat; iR++){         // make nRepeat random permutations of points
    for (Int_t iPoint=0; iPoint<nPoints; iPoint++){
      xrndm[iPoint]=gRandom->Rndm();
    }
    TMath::Sort(nPoints, xrndm, &(permutationIndex[iR*nPoints]));
  }
  // 1.) Calculate fit parameters for small subset of data 
  for (Int_t iter=0; iter<nIter; iter++){
    fPointIndex.Set(nFitPoints, &(permutationIndex[iter*nFitPoints]));
    Fit(option); 
    fitStatus[iter]=fStatus;
    FitterFCN(nPar,NULL,fChi2, fParam->GetMatrixArray(),0);
    logLike[iter]=fChi2;
    paramArray[iter] = new TMatrixD(nPar,1,fParam->GetMatrixArray());
    covarArray[iter] = new TMatrixD(*fCovar);    
  }
  // 2. Calculate normalized distance between each fits
  std::vector< std::vector<double> > logDistance(nIter, std::vector<double>(nIter));  //store vectors and chi2
  for (Int_t iter0=0; iter0<nIter; iter0++){    
    for (Int_t iter1=0; iter1<nIter; iter1++){
      if ( ((fitStatus[iter0]&kCovarFailed) || (fitStatus[iter1]&kCovarFailed))==0){
	TMatrixD covar(*(covarArray[iter0]));
	covar+=*(covarArray[iter1]);
	covar.Invert();
	TMatrixD delta(*(paramArray[iter0]));
	delta-=*(paramArray[iter1]);
	TMatrixD delta2(covar,TMatrixD::kMult,delta);
	delta.T();
	TMatrixD chi2(delta,TMatrixD::kMult,delta2);
	chi2*=logLike[iter0]+logLike[iter1];
	Double_t normDistance=TMath::Sqrt(chi2(0,0));      
	logDistance[iter0][iter1]=normDistance;
      }else{
	logDistance[iter0][iter1]=-1.;
      }
    }
  }
  // 3. Reject outliers and calculate robust mean
  for (Int_t iter=0; iter<nIter; iter++){    
    medDistance[iter]=TMath::Median(nIter, logDistance[iter].data());
  }
  Double_t medMedDistance=TMath::Median(nIter,medDistance.data());
  // 4.  Extract combined statistic
  // 4.a 1D median,mean,rms
  // 4.b nD ?
  if (fMISACEstimators==NULL) fMISACEstimators = new TMatrixD(4,nPar);
  for (Int_t iPar=0;iPar<nPar; iPar++){
    std::vector<double>workingArray(nIter);
    Int_t nAccepted=0;
    for (Int_t iter=0; iter<nIter; iter++){
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
    for (Int_t iter=0; iter<nIter; iter++){ 
      (*fStreamer)<<"misac"<<
	"iter="<<iter<<          // iteration 
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



void  AliTMinuitToolkit::RegisterDefaultFitters(){
  //
  //
  TF1 *likeGausCachy = new TF1("likeGausCachy", AliTMinuitToolkit::GaussCachyLogLike,-10,10,2);
  likeGausCachy->SetParameters(0.8,1);
  //
  for (Int_t ipol=0; ipol<10; ipol++){ //register polynomial fitters
    TMatrixD initPar(ipol+1,4);
    for (Int_t ipar=0; ipar<ipol+1; ipar++) initPar(ipar,1)=1; 
    //
    TF1 *fpol = new TF1(TString::Format("fpol%d",ipol).Data(),TString::Format("pol%d",ipol).Data());
    AliTMinuitToolkit * fitter1D = new AliTMinuitToolkit();
    fitter1D->SetVerbose(0x1); fitter1D->SetFitFunction(fpol,kTRUE);
    fitter1D->SetInitialParam(&initPar);
    fitter1D->SetName(TString::Format("pol%d",ipol).Data());
    AliTMinuitToolkit::SetPredefinedFitter(fitter1D->GetName(),fitter1D);
    // gaus log likelihood cost function
    AliTMinuitToolkit * fitter1DR = new AliTMinuitToolkit();
    fitter1DR->SetVerbose(0x1); fitter1DR->SetFitFunction(fpol,kTRUE);
    fitter1DR->SetLogLikelihoodFunction(likeGausCachy);
    fitter1DR->SetInitialParam(&initPar);
    fitter1DR->SetName(TString::Format("pol%dR",ipol).Data());
    AliTMinuitToolkit::SetPredefinedFitter(fitter1DR->GetName(),fitter1DR);
    // huber cost function
    AliTMinuitToolkit * fitter1DH = new AliTMinuitToolkit();
    fitter1DH->SetVerbose(0x1); fitter1DH->SetFitFunction(fpol,kTRUE);
    fitter1DH->SetInitialParam(&initPar);
    fitter1DH->EnableRobust(true);
    fitter1DH->SetName(TString::Format("pol%dH",ipol).Data());
    AliTMinuitToolkit::SetPredefinedFitter(fitter1DH->GetName(),fitter1DH);
  }
  //
  TF1 *fgaus = new TF1("fgaus","gaus");
  TMatrixD initPar(3,4);
  initPar(0,0)=0; initPar(0,1)=1; 
  initPar(1,0)=0; initPar(1,1)=1; 
  initPar(2,0)=1; initPar(2,1)=1; 
  AliTMinuitToolkit * fitterGR = new AliTMinuitToolkit();
  fitterGR->SetVerbose(0x1); fitterGR->SetFitFunction(fgaus,kTRUE);
  fitterGR->SetLogLikelihoodFunction(likeGausCachy);
  fitterGR->SetInitialParam(&initPar);
  AliTMinuitToolkit::SetPredefinedFitter("gausR",fitterGR);
 //
}


AliTMinuitToolkit *  AliTMinuitToolkit::Fit(TH1* his, const char *fitterName, Option_t* option, Option_t* goption,Option_t* foption, Double_t xmin, Double_t xmax){
  //
  //
  AliTMinuitToolkit *fitter=GetPredefinedFitter(fitterName);
  if (fitter==NULL){
    ::Error("AliTMinuitToolkit::Fit","Unsuported fitter %s. \nfitter can be added to the list using. SetPredefinedFitter", fitterName);
    return NULL;
  }
  if (his==NULL){
    ::Error("AliTMinuitToolkit::Fit","Zerro pointer");
    return NULL;
  }
  fitter->FitHistogram(his,foption);
  TF1 * fitFun = (TF1*)(fitter->GetFormula()->Clone());
  fitFun->SetParameters(fitter->GetParameters()->GetMatrixArray());
  SetFunctionDrawOption(fitFun,goption);
  his->GetListOfFunctions()->AddLast(fitFun);
  if (xmin<xmax) {
    fitFun->SetRange(xmin,xmax);
  }else{
    fitFun->SetRange(his->GetXaxis()->GetXmin(), his->GetXaxis()->GetXmax());
  }
  his->Draw(goption);
  return fitter;
}



AliTMinuitToolkit *  AliTMinuitToolkit::Fit(TGraph *graph, const char *fitterName, Option_t* option, Option_t* goption,  Option_t* foption, Double_t xmin, Double_t xmax){
  //
  //
  AliTMinuitToolkit *fitter=GetPredefinedFitter(fitterName);
  if (fitter==NULL){
    ::Error("AliTMinuitToolkit::Fit","Unsuported fitter %s. \nfitter can be added to the list using. SetPredefinedFitter", fitterName);
    return NULL;
  }
  if (graph==NULL){
    ::Error("AliTMinuitToolkit::Fit","Zerro pointer");
    return NULL;
  }
  fitter->FitGraph(graph,option);
  TF1 * fitFun = (TF1*)(fitter->GetFormula()->Clone());
  fitFun->SetParameters(fitter->GetParameters()->GetMatrixArray());
  fitFun->SetName(TString::Format("%s:%s",fitter->GetName(), option).Data());
  fitFun->SetTitle(TString::Format("%s:%s",fitter->GetName(), option).Data());
  SetFunctionDrawOption(fitFun,foption);
  graph->GetListOfFunctions()->AddLast(fitFun);
  if (xmin<xmax) {
    fitFun->SetRange(xmin,xmax);
  }else{
    fitFun->SetRange(graph->GetXaxis()->GetXmin(), graph->GetXaxis()->GetXmax());
  }
  graph->Draw();
  return fitter;
}

void  AliTMinuitToolkit::SetFunctionDrawOption(TF1 *fun, Option_t *option){
  TString funOption=option; // e.g  "funOption(1,1,1)"
  TPRegexp regFunOption("funOption\\(");
  // Error measges ???
  if (regFunOption.Match(funOption)){
    Int_t index0=funOption.Index(regFunOption)+9;
    Int_t index1=funOption.Index(")",index0+1);
    if (index1<index0) {
      ::Error("AliTMinuitToolkit::SetFunctionDrawOption","Invalide function draw option syntax %s", option);
      return;
    }
    Int_t index=index0+1;
    Int_t color=TString(&(funOption.Data()[index])).Atoi();
    if (color>0) fun->SetLineColor(color);
    index=funOption.Index(",",index+1)+1;
    if (index>0) {
      Int_t width=TString(&(funOption.Data()[index])).Atoi();
      if (width>0) fun->SetLineWidth(width);
      index=funOption.Index(",",index+1)+1;
      if (index>0) {
	Int_t style=TString(&(funOption.Data()[index])).Atoi();
	if (style>0) fun->SetLineStyle(style);
      }
    }
  }

}
