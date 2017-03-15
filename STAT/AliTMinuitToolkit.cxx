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

AliTMinuitToolkit::AliTMinuitToolkit(const char *streamerName) : 
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
 if (streamerName!=NULL) fStreamer= new TTreeSRedirector(streamerName,"update");
 
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
  delete fCovar;
  if (fStreamer) delete fStreamer;
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

void  AliTMinuitToolkit::SetFitFunction(TFormula *const formula, Bool_t doReset) {
  //
  fFormula=formula;
  if (doReset){
    delete fParam;
    delete fCovar;
    fParam=new TVectorD(formula->GetNpar());
    fCovar=new TMatrixD(formula->GetNpar(),formula->GetNpar());
  }
}


void   AliTMinuitToolkit::SetInitialParam(TMatrixD *const paramLimits) { 
  fInitialParam=new TMatrixD(*paramLimits);
};  



void AliTMinuitToolkit::FitHistogram(TH1F *const his) {
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
  Fit(kTRUE);
}






void AliTMinuitToolkit::FitterFCN(int &/*npar*/, double *info, double &fchisq, double *gin, int /*iflag*/){
  //
  // internal function which gives the specified function to the TMinuit function
  //  
  //
  static Int_t fcnCounter=0;
  if (info) info[0]=fcnCounter;
  fcnCounter++;
  AliTMinuitToolkit * fitter = (AliTMinuitToolkit*)TVirtualFitter::GetFitter()->GetObjectFit();
  TFormula *logLike=fitter->fLogLikelihoodFunction;
  const Double_t* likeParam= (logLike!=NULL) ? logLike->GetParameters():NULL;
  fchisq = 0;
  const TMatrixD & variables= (*fitter->GetPoints());
  const TMatrixD & values=    (*fitter->GetValues());
  Int_t nvars       = variables.GetNcols();
  Int_t npoints     = variables.GetNrows();
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
    if (logLike){
      Double_t normDelta = delta*weight;
      fchisq+=logLike->EvalPar(&delta,likeParam);
      continue;
    }
    if (fitter->IsHuberCost() == true) {       //hubert norm
     delta = delta*weight;                   // normalization
     if (delta <= 2.5) fchisq+= delta*delta; // new metric: Huber-k-estimator
     if (delta > 2.5) fchisq+= 2*(2.5)*delta - (2.5*2.5);
    } else {
     Double_t chi2 = delta*weight;
     chi2*=chi2;
     fchisq+= chi2;   // chi2 (log likelihood)
    }
  }
  fcnCounter++;
}
 

void AliTMinuitToolkit::Fit(Bool_t doReset) {
  //
  // internal function that calls the fitter
  //
  
  //  Int_t nParam  = fParam->GetNrows();
  Int_t nParam  = fFormula->GetNpar();
  Int_t npoints = fPoints->GetNrows();
  Int_t nvar    = fPoints->GetNcols()-1;
  
  // create "intial param" in case not set before
  if (fParam==NULL){
    fParam=new TVectorD(nParam);
  }
  if (fInitialParam == 0) {
    ::Info("AliTMinuitToolkit::Fit","Default initial parameters set");
    fInitialParam=  new TMatrixD(nParam,4);
    for (Int_t iparam=0; iparam<nParam; iparam++) (*fInitialParam)(iparam,1)=0;
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
  
  // initialize parameters (step size???)
  for (Int_t iparam=0; iparam<nParam; iparam++){
    if (doReset){
      minuit->SetParameter(iparam, Form("p[%d]",iparam), (*fInitialParam)(iparam,0), (*fInitialParam)(iparam,1), (*fInitialParam)(iparam, 2), (*fInitialParam)(iparam, 3));
    }else{
      minuit->SetParameter(iparam, Form("p[%d]",iparam), (*fParam)(iparam), TMath::Sqrt((*fCovar)(iparam,iparam)), (*fInitialParam)(iparam, 2), (*fInitialParam)(iparam, 3));
    }
  }
  
  //
  Double_t argList[2];
  argList[0] = fMaxCalls; //maximal number of calls 
  argList[1] = fPrecision; //tolerance normalized to 0.001 
  if (fMaxCalls == 500 && fPrecision == 1) minuit->ExecuteCommand(fFitAlgorithm, 0, 0); 
  if (fMaxCalls != 500 || fPrecision != 1) minuit->ExecuteCommand(fFitAlgorithm, argList, 2);
  // two additional arguments can be specified ExecuteCommand("migrad", argList, 2) - use 0,0 for default
  fStatus = (((TFitter*)minuit)->GetMinuit())->GetStatus(); // temporary hack XXXXXXXXXXXXXXXXXXx
  // fill parameter vector
  for (Int_t ivar=0; ivar<nParam; ivar++){
    (*fParam)(ivar) = minuit->GetParameter(ivar);
    fFormula->SetParameter(ivar, minuit->GetParameter(ivar));
  }
  
  
  // fill parameter vector
  for (Int_t ivar=0; ivar<nParam; ivar++){
   (*fParam)(ivar) = minuit->GetParameter(ivar);
   fFormula->SetParameter(ivar, minuit->GetParameter(ivar));
  }
  
  // fill covariance matrix
  fCovar = new TMatrixD(nParam, nParam);
  for(Int_t i=0; i < nParam; i++) {
   for(Int_t j=0; j < nParam; j++) {
    (*fCovar)(i,j) = minuit->GetCovarianceMatrixElement(i,j);
   }
  }
  
}


Int_t AliTMinuitToolkit::FillFitter(TTree * inputTree, TString values, TString variables, TString selection, Int_t firstEntry, Int_t nentries ){
  //
  // Make unbinned fit
  //
  ClearData();
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
  Int_t entries = inputTree->Draw(query.Data(),selection.Data(),"goffpara",nentries,firstEntry);
  if (entries<=0) {
    ::Error("AliTMinuitToolkit::UnbinnedFit","badly formatted values or variables: %s",query.Data());
    ::Error("AliTMinuitToolkit::UnbinnedFit","valueDescription: %s",values.Data());
    ::Error("AliTMinuitToolkit::UnbinnedFit","variables: %s",variables.Data());
    return -1;
  }  
  fPoints=new TMatrixD(entries,nVars);
  fValues=new TMatrixD(entries,nVals);
  for (Int_t iPoint=0; iPoint<entries; iPoint++){
    for (Int_t iVar=0; iVar<nVars; iVar++){
      (*fPoints)(iPoint,iVar)=inputTree->GetVal(iVar+nVals)[iPoint];
    }
    for (Int_t iVal=0; iVal<nVals; iVal++){
      (*fValues)(iPoint,iVal)=inputTree->GetVal(iVal)[iPoint];
    }
  }
}

TString AliTMinuitToolkit::GetFitFunctionAsAlias(){
  //
  // construct string TTree alias for fit function
  TString inputString(fFormula->GetTitle());

  if (fVarNames==NULL){
    ::Error("AliTMinuitToolkit::GetFitFunctionAsAlias","Variable names not defined. Fucntion supported for TTree::UnbinnedFit");
    return "";
  }
  for (Int_t iDim=0; iDim<fFormula->GetNdim(); iDim++){
    inputString.ReplaceAll(TString::Format("x[%d]",iDim).Data(), GetListOfVariables()->At(iDim)->GetName());
  }
  for (Int_t iPar=0; iPar<fFormula->GetNpar(); iPar++){
    inputString.ReplaceAll(TString::Format("[%d]",iPar).Data(), TString::Format("(%f)",(*fParam)[iPar]).Data());
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
}

void AliTMinuitToolkit::Bootstrap(Int_t nIter, const char * reportName){
  //
  // Bootstrap minimization 
  //
  static Int_t counter=0;
  Int_t nPoints= fPoints->GetNrows();
  fPointIndex.Set(nPoints);
  Double_t info[3]={0,0,0};
  Int_t npar=fFormula->GetNpar();
  Int_t fcnP=-1;
  TObjString rName(reportName);
  for (Int_t iter=0; iter<nIter; iter++){
    counter++;
    for (Int_t iPoint=0; iPoint<nPoints; iPoint++){
      fPointIndex[iPoint]=gRandom->Rndm()*nPoints;
    }
    Fit(kTRUE);
    FitterFCN(npar,info,fChi2, fParam->GetMatrixArray(),0); 
    if (fcnP<0) fcnP=info[0];
    if (fVerbose&kSysInfo) AliSysInfo::AddStamp("TwoFoldCrossValidation",0,counter,iter,info[0]-fcnP);
    Int_t fcnCounter=info[0]-fcnP;
    fcnP=info[0];
    if (fStreamer){      
      (*fStreamer)<<"bootstap"<<
	"reportName.="<<&rName<<
	"iter="<<iter<<
	"fStatus="<<fStatus<<    // minuit status
	"fChi2="<<fChi2<<
	"fcnCounterI="<<info[0]<<
	"fcnCounter="<<fcnCounter<<
	"fParam.="<<fParam<<
	"fCovar.="<<fCovar<<
	"\n";
    }
  }
  return ;
}


void AliTMinuitToolkit::TwoFoldCrossValidation(Int_t nIter, const char * reportName){
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
  Int_t npar=fFormula->GetNpar();
  Int_t fcnP=-1;
  TObjString rName(reportName);
  for (Int_t iter=0; iter<nIter; iter++){
    for (Int_t iPoint=0; iPoint<nPoints; iPoint++){
      rndmCV[iPoint]=gRandom->Rndm()*nPoints;
    }
    TMath::Sort(nPoints,rndmCV.GetArray(), indexFold.GetArray());
    //
    fPointIndex.Set(nPoints/2,indexFold.GetArray());
    Fit(kTRUE);
    TVectorD param0(*fParam);
    TMatrixD covar0(*fCovar);
    FitterFCN(npar,info,chi2_00, fParam->GetMatrixArray(),0); 
    //
    fPointIndex.Set(nPoints-nPoints/2,&(indexFold.GetArray()[nPoints/2]));
    FitterFCN(npar,info,chi2_01, fParam->GetMatrixArray(),0); 
    Fit(kTRUE);
    FitterFCN(npar,info,chi2_11, fParam->GetMatrixArray(),0); 
    TVectorD param1(*fParam);
    TMatrixD covar1(*fCovar);
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
	"chi2_00="<<chi2_00<<    // log likelihodd for training sample
	"chi2_01="<<chi2_01<<    // log likelihood for test sample using traning fit
	"chi2_11="<<chi2_11<<    // log likelihood for test sample using test fit
	"fcnCounterI="<<info[0]<<   // fcn counter integral
	"fcnCounter="<<fcnCounter<< // fcn counter per fit
	"param0.="<<&param0<<    // parameters estimate training sample
	"covar0.="<<&covar0<<    // covariance for training sample
	"param1.="<<&param1<<    // parameters for complementaray subsample
	"covar1.="<<&covar1<<    // covariance for complementaray subsample
	"\n";
    }
  }
  return ;
}


