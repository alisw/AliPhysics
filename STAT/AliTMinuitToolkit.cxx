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

#include "AliTMinuitToolkit.h"

//--------------------------------------------------------------------------------------
//
// The AliTMinuitToolkit serves as an easy to use interface for the TMinuit 
// package: 
// 
// - It allows to fit a curve to one and two dimensional histograms 
//   (TH2F::Fit() only allows to fit a hyperplane).
// - Or n points can be specified directly via a n x 2 matrix.
// - An option for robust fitting with non-linear functions is implemented.
//
// A small example illustrating the usage of AliTMinuitToolkit is given in the function 
// "AliTMinuitToolkit::Test()".
// 
// 
// 1. Setting the formula:
//
//  The formula is simply set via "void SetFitFunction(TFormula * formula)".
//
//
// 2. Adding the data points
//
//  - In order to fit a histogram, use "void FitHistogram(TH1F * his)" or 
//    "void FitHistogram(TH2F * his)". The fitter is started automatically
//  - Alternatively, the direct specification of the points is possible via
//    "void SetPoints(TMatrixD * points)". Note, that the each point 
//    corresponds to one row in the matrix. The fitter is then started with 
//    the command "void Fit()". The weight of each point can be specified
//    with an n-dimensional vector using "void SetWeights(TVectorD * weights)".
//
//
// 3. Accessing the fit results
//
//  The N parameters of the formula are stored in a N-dimensional vector which
//  is returned by "TVectorD * GetParameters()". In a similar way the covariance 
//  matrix of the fit is returned via "TMatrixD * GetCovarianceMatrix()" which
//  is of the type N x N.
//
//
// 4. Non-linear robust fitting:
//
//  Even a few outliers can lead to wrong results of a least-squares fitting 
//  procedure. In this case the use of robust(resistant) methods can be 
//  helpful, but a stronger dependence on starting values or convergence to
//  local minima can occur.
//
//  The robust option becomes active if EnableRobust(true, sigma) is called. It is
//  very much recommended that a normalization value (scale variable) corresponding 
//  to an expected deviation (sigma) is specified via 
//  "EnableRobust(Bool_t b, Double_t sigma)".
//
//  Performing the fit without knowledge of sigma is also possible if only
//  "EnableRobust(true)" is activated, but this is NOT RECOMMENDED.
//
//  The method is based on another estimator instead of chi^2. For small deviations 
//  the function behaves like x^2 and for larger deviations like |x| - the so 
//  called Huber estimator:
//
//   h(x) = x^2                              , for x < 2.5*sigma
//   h(x) = 2*(2.5*sigma)*x - (2.5*sigma)^2  , for x > 2.5*sigma
//
//  If a weighting function is specified in addition, a second run with the ordinary 
//  metric is started, but before entering the iteration every point is weighted 
//  according to its distance to the outcoming function of the first run. The weighting
//  function w(x) must be defined on the intervall x in [0,1]. w(0) then 
//  corresponds to the weight of the closest point and w(1) to the point with the
//  largest distance.
//
//  Some standard weighting functions are predefined in 
//  "SetWeightFunction(Char_t * name, Float_t param1, Float_t param2 = 0)":
//   - "BOX" equals to 1 if x < param1 and to 0 if x > param1.
//   - "EXPONENTIAL" corresponds to "Math::Exp(-TMath::Log(param1)*x)"
//   - "ERRORFUNCTION" corresponds to "TMath::Erfc((x-param1)/param2)"
//
//
//  REFERENCE for non-linear robust fitting:
//  Ekblom H. and Madsen K. (1988), Alogrithms for non-linear Huber estimation,
//  BIT Numerical Mathematics 29 (1989) 60-76.
//  internet: http://www.springerlink.com/content/m277218542988344/
//
//
// 5. examples:
//
//  A small example illustrating the working principles of AliTMinuitToolkit is given
//  in the function "AliTMinuitToolkit::Test()".
//
//
//
// Comments and questions are always welcome: A.Kalweit@gsi.de
//--------------------------------------------------------------------------------------


ClassImp(AliTMinuitToolkit)

AliTMinuitToolkit::AliTMinuitToolkit() : 
   TNamed(),
   fFormula(0),
   fWeightFunction(0),
   fFitAlgorithm(""),
   fPoints(0),
   fWeights(0),
   fParam(0),
   fParamLimits(0),
   fCovar(0),
   fChi2(0),
   fMaxCalls(0),
   fPrecision(0),
   fUseRobust(0),
   fExpectedSigma(0)
{
 //
 // standard constructor
 //
 fMaxCalls = 500;
 fPrecision = 1;
 fUseRobust = false;
 fExpectedSigma = 0;
}


AliTMinuitToolkit::AliTMinuitToolkit(const AliTMinuitToolkit&) :
   TNamed(),
   fFormula(0),
   fWeightFunction(0),
   fFitAlgorithm(""),
   fPoints(0),
   fWeights(0),
   fParam(0),
   fParamLimits(0),
   fCovar(0),
   fChi2(0),
   fMaxCalls(0),
   fPrecision(0),
   fUseRobust(0),
   fExpectedSigma(0)
{


}


AliTMinuitToolkit& AliTMinuitToolkit::operator=(const AliTMinuitToolkit&) {

 return *this;
}



AliTMinuitToolkit::~AliTMinuitToolkit(){
  //
  // destructor
  //
  delete fPoints;
  delete fWeights;
  delete fWeightFunction;
  delete fParamLimits;
  delete fFormula;
  delete fParam;
  delete fCovar;
  delete fChi2;
}

void AliTMinuitToolkit::FitHistogram(TH1F *const his) {
 //
 // Fit a one dimensional histogram
 //
 fPoints = new TMatrixD(his->GetNbinsX(), 2);
 
 for(Int_t ibin=0; ibin < his->GetNbinsX(); ibin++) {
  Double_t x = his->GetXaxis()->GetBinCenter(ibin+1);
  Double_t y = his->GetBinContent(ibin+1);
  
  (*fPoints)(ibin, 0) = x;
  (*fPoints)(ibin, 1) = y;
 }
 
 Fit();
}


void AliTMinuitToolkit::FitHistogram(TH2F *const his) {
 //
 // Fit a curve to a two dimensional histogram
 //
 fPoints = new TMatrixD((Long64_t)his->GetEntries(), 2);
 Long64_t entry = 0;
 
 for(Int_t ibin=0; ibin < his->GetNbinsX(); ibin++) {
  Double_t x = his->GetXaxis()->GetBinCenter(ibin);
  for(Int_t jbin=0; jbin < his->GetNbinsY(); jbin++) {   
   Long64_t n = his->GetBin(ibin, jbin);
   Double_t y = his->GetYaxis()->GetBinCenter(jbin);
   for(Int_t ientries=0; ientries < his->GetBinContent(n); ientries++) {
    (*fPoints)(entry,0) = x;
    (*fPoints)(entry,1) = y;
    entry++;
   }
   
  }
 }

 Fit();
}


void AliTMinuitToolkit::SetWeightFunction(const Char_t *name, Float_t param1, Float_t param2) {
 //
 // Set the weight function which must be defined on the interval [0,1].
 //
 TString FuncType(name);
 FuncType.ToUpper();
 
 if (FuncType == "EXPONENTIAL") fWeightFunction = new TFormula("exp", Form("TMath::Exp(-TMath::Log(%f)*x)", param1));
 if (FuncType == "BOX") fWeightFunction = new TFormula("box", Form("TMath::Erfc((x-%f)/0.0001)", param1));
 if (FuncType == "ERRORFUNCTION") fWeightFunction = new TFormula("err", Form("TMath::Erfc((x-%f)/%f)", param1, param2));
 
}


void AliTMinuitToolkit::FitterFCN(int &/*npar*/, double */*dummy*/, double &fchisq, double *gin, int /*iflag*/){
  //
  // internal function which gives the specified function to the TMinuit function
  //  

  //
  AliTMinuitToolkit * fitter = (AliTMinuitToolkit*)TVirtualFitter::GetFitter()->GetObjectFit();
  fchisq = 0;
  Int_t nvar       = fitter->GetPoints()->GetNcols()-1;
  Int_t npoints    = fitter->GetPoints()->GetNrows();
  
  // calculate mean deviation for normalization or use user-defined sigma
  Double_t dev = 0.;
  if (fitter->GetExpectedSigma() == 0 && fitter->GetStatus() == true) {
   for (Int_t ipoint=0; ipoint<npoints; ipoint++){
    Double_t x[100];
     for (Int_t ivar=0; ivar<nvar; ivar++){
      x[ivar] = (*fitter->GetPoints())(ipoint, ivar);      
     }
    Float_t funx = fitter->GetFormula()->EvalPar(x,gin);
    Double_t delta = (*fitter->GetPoints())(ipoint, nvar) - funx;
    dev += TMath::Sqrt(TMath::Abs(delta));
   }
   dev = dev/npoints; 
  } else {
   dev = fitter->GetExpectedSigma();
  }
  // calculate chisquare  
  for (Int_t ipoint=0; ipoint<npoints; ipoint++){
    Double_t x[100];
    for (Int_t ivar=0; ivar<nvar; ivar++){
      x[ivar] = (*fitter->GetPoints())(ipoint, ivar);      
    }
    Float_t funx = fitter->GetFormula()->EvalPar(x,gin);   
    Double_t delta = TMath::Abs((*fitter->GetPoints())(ipoint, nvar) - funx);
    if (fitter->GetStatus() == true) {
     delta = delta/dev; // normalization
     if (delta <= 2.5) fchisq+= delta*delta; // new metric: Huber-k-estimator
     if (delta > 2.5) fchisq+= 2*(2.5)*delta - (2.5*2.5);
    } else {
     Double_t weight = (*fitter->GetWeights())(ipoint);
     fchisq+= delta*delta*weight; //old metric
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
  
  // set all weights to 1 as default
  Bool_t weightFlag = false;
  if (fWeightFunction == 0) {
   fWeightFunction = new TFormula("constant", "1");
  } else {
   weightFlag = true;
  }
  
  // migrad fit algorithm as default
  if (fFitAlgorithm == "") {
   fFitAlgorithm = "migrad";
  }
  
  // assign weights
  if (fWeights == 0) {
   fWeights = new TVectorD(npoints);
   for (Int_t ipoint=0; ipoint<npoints; ipoint++) (*fWeights)(ipoint) = 1;
  }
  
  // set up the fitter
  TVirtualFitter *minuit = TVirtualFitter::Fitter(0, nparam);
  minuit->SetObjectFit(this); 
  minuit->SetFCN(AliTMinuitToolkit::FitterFCN);
  
  // initialize paramters (step size???)
  for (Int_t iparam=0; iparam<nparam; iparam++){
   minuit->SetParameter(iparam, Form("p[%d]",iparam), (*fParam)(iparam), (*fParam)(iparam)/10, (*fParamLimits)(iparam, 0), (*fParamLimits)(iparam, 1));
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

  // if a weight function is specified -> enter 2nd run with weights
  if (weightFlag == true && fUseRobust == true) {
   // sort points for weighting
   Double_t *sortList = new Double_t[npoints];
   Int_t  *indexList = new Int_t[npoints];   
   for (Int_t ipoint=0; ipoint<npoints; ipoint++){
    Double_t funx = fFormula->Eval((*fPoints)(ipoint, 0));
    Double_t delta = TMath::Abs((*fPoints)[ipoint][nvar] - funx);
    sortList[ipoint] = delta;
   } 
   TMath::Sort(npoints, sortList, indexList, false);
   for (Int_t ip=0; ip<npoints; ip++){
    Double_t t = ip/(Double_t)npoints;
    (*fWeights)(indexList[ip]) = fWeightFunction->Eval(t);
   }
   
   // set up the fitter
   fUseRobust = false;
   for (Int_t iparam=0; iparam<nparam; iparam++){
    minuit->SetParameter(iparam, Form("p[%d]",iparam), (*fParam)(iparam), (*fParam)(iparam)/10, (*fParamLimits)(iparam, 0), (*fParamLimits)(iparam, 1));
   }
   // start fitting
   if (fMaxCalls == 500 && fPrecision == 1) minuit->ExecuteCommand(fFitAlgorithm, 0, 0); 
   if (fMaxCalls != 500 || fPrecision != 1) minuit->ExecuteCommand(fFitAlgorithm, argList, 2);
   fUseRobust = true;
   
   delete [] sortList; 
   delete [] indexList;    
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
  
  if (weightFlag == false) fWeightFunction = 0;
}



void AliTMinuitToolkit::Test() {
 //
 // This test function shows the basic working principles of this class 
 // and illustrates how a robust fit can improve the results
 //
 
 // 1. provide some example histogram
 TH1F * hist = new TH1F("test", "with (red) and without (black) robust option", 20,0,4);
 TRandom * rand = new TRandom();
 for (Int_t i = 0; i < 10000; i++) {
  hist->Fill(rand->Exp(1));
  if (i < 1000) hist->Fill(3); //"outliers"
  if (i < 1070) hist->Fill(3.5);
  if (i < 670) hist->Fill(2);
  if (i < 770) hist->Fill(1.5);//"outliers"
  if (i < 740) hist->Fill(1);
 }
 TCanvas * canv = new TCanvas();
 canv->cd(1);
 hist->Draw();
 
 // 2. example fit without robust option
 AliTMinuitToolkit * tool = new AliTMinuitToolkit();
 TFormula *aFormExp = new TFormula("formExp", "[0]*TMath::Exp(-[1]*x)");
 tool->SetFitFunction(aFormExp);
 TVectorD *vec1 = new TVectorD(2); // Set initial values
 (*vec1)(0) = 1800;
 (*vec1)(1) = 1;
 tool->SetInitialParam(vec1);
 tool->FitHistogram(hist);
 
 // draw fit function
 TF1 *func = new TF1("test", "[0]*TMath::Exp(-[1]*x)", 0, 6);
 func->SetParameters((*tool->GetParameters())(0), (*tool->GetParameters())(1));
 func->Draw("same");
 
 // 3 . robust fit 
 TVectorD *vec2 = new TVectorD(2);
 (*vec2)(0) = 1800;
 (*vec2)(1) = 1;
 tool->SetInitialParam(vec2);
 tool->EnableRobust(true, 10);
 tool->SetWeightFunction("box", 0.75);
 tool->FitHistogram(hist);
 TF1 *func2 = new TF1("test2", "[0]*TMath::Exp(-[1]*x)", 0, 6);
 func2->SetParameter(0, (*tool->GetParameters())(0));
 func2->SetParameter(1, (*tool->GetParameters())(1));
 func2->SetLineColor(kRed);
 func2->Draw("same");
 
}

