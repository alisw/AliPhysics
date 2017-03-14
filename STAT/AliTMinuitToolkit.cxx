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

AliTMinuitToolkit::AliTMinuitToolkit() : 
   TNamed(),
   fVerbose(0),
   fFormula(0),
   fLogLikelihoodFunction(0),
   fFitAlgorithm(""),
   fPoints(0),
   fValues(0),
   fParam(0),
   fParamLimits(0),
   fCovar(0),
   fChi2(0),
   fMaxCalls(0),
   fPrecision(0),
   fUseRobust(0)
{
 //
 // standard constructor
 //
 fMaxCalls = 500;
 fPrecision = 1;
 fUseRobust = false;
}


AliTMinuitToolkit::~AliTMinuitToolkit(){
  //
  // destructor
  //
  delete fPoints;
  delete fValues;
  delete fParamLimits;
  delete fFormula;
  delete fParam;
  delete fCovar;
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

void   AliTMinuitToolkit::SetInitialParam(TVectorD *const param) { 
  fParam=new TVectorD(*param);
};

void   AliTMinuitToolkit::SetParamLimits(TMatrixD *const paramLimits) { 
  fParamLimits=paramLimits;
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
  Fit();
}






void AliTMinuitToolkit::FitterFCN(int &/*npar*/, double */*dummy*/, double &fchisq, double *gin, int /*iflag*/){
  //
  // internal function which gives the specified function to the TMinuit function
  //  

  //
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
  for (Int_t ipoint=0; ipoint<npoints; ipoint++){
    Double_t x[100];
    for (Int_t ivar=0; ivar<nvars; ivar++){
      x[ivar] = variables(ipoint, ivar);      
    }    
    Float_t funx = fitter->GetFormula()->EvalPar(x,gin);   
    Double_t value=values(ipoint,0);
    Double_t weight=values(ipoint,1);
    Double_t delta = TMath::Abs(value - funx);
    if (logLike){
      Double_t normDelta = delta*weight;
      fchisq+=logLike->EvalPar(&delta,likeParam);
      continue;
    }
    if (fitter->GetStatus() == true) {       //hubert norm
     delta = delta*weight;                   // normalization
     if (delta <= 2.5) fchisq+= delta*delta; // new metric: Huber-k-estimator
     if (delta > 2.5) fchisq+= 2*(2.5)*delta - (2.5*2.5);
    } else {
     Double_t chi2 = delta*weight;
     chi2*=chi2;
     fchisq+= chi2;   // chi2 (log likelihood)
    }
  }
 }
 

void AliTMinuitToolkit::Fit() {
  //
  // internal function that calls the fitter
  //
  Int_t nparam  = fParam->GetNrows();
  Int_t npoints = fPoints->GetNrows();
  Int_t nvar    = fPoints->GetNcols()-1;
  
  // set all paramter limits to infinity as default
  if (fParamLimits == 0) {
   fParamLimits = new TMatrixD(nparam ,2);
   for (Int_t iparam=0; iparam<nparam; iparam++){
    (*fParamLimits)(iparam, 0) = 0;
    (*fParamLimits)(iparam, 1) = 0;
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
  TVirtualFitter *minuit = TVirtualFitter::Fitter(0, nparam);
  minuit->SetObjectFit(this); 
  minuit->SetFCN(AliTMinuitToolkit::FitterFCN);
  
  // initialize paramters (step size???)
  for (Int_t iparam=0; iparam<nparam; iparam++){
    minuit->SetParameter(iparam, Form("p[%d]",iparam), (*fParam)(iparam), (*fParam)(iparam)/10+0.00000001, (*fParamLimits)(iparam, 0), (*fParamLimits)(iparam, 1));

  //   if (fParamLimits){
//       minuit->SetParameter(iparam, Form("p[%d]",iparam), (*fParam)(iparam), (*fParam)(iparam)/10, (*fParamLimits)(iparam, 0), (*fParamLimits)(iparam, 1));
//     }
//     //   else{
//     //       minuit->SetParameter(iparam, Form("p[%d]",iparam), (*fParam)(iparam), );
//     //     }

  }
  
  //
  Double_t argList[2];
  argList[0] = fMaxCalls; //maximal number of calls 
  argList[1] = fPrecision; //tolerance normalized to 0.001 
  if (fMaxCalls == 500 && fPrecision == 1) minuit->ExecuteCommand(fFitAlgorithm, 0, 0); 
  if (fMaxCalls != 500 || fPrecision != 1) minuit->ExecuteCommand(fFitAlgorithm, argList, 2);
  // two additional arguments can be specified ExecuteCommand("migrad", argList, 2) - use 0,0 for default
  
  // fill parameter vector
  for (Int_t ivar=0; ivar<nparam; ivar++){
    (*fParam)(ivar) = minuit->GetParameter(ivar);
    fFormula->SetParameter(ivar, minuit->GetParameter(ivar));
  }
  
  
  // fill parameter vector
  for (Int_t ivar=0; ivar<nparam; ivar++){
   (*fParam)(ivar) = minuit->GetParameter(ivar);
   fFormula->SetParameter(ivar, minuit->GetParameter(ivar));
  }
  
  // fill covariance matrix
  fCovar = new TMatrixD(nparam, nparam);
  for(Int_t i=0; i < nparam; i++) {
   for(Int_t j=0; j < nparam; j++) {
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


