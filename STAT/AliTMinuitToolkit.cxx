
#include "AliTMinuitToolkit.h"
#include <TNamed.h>
#include <TVirtualFitter.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFormula.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMath.h>
#include <TString.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TRandom.h>


//--------------------------------------------------------------------------
//
// The AliTMinuitToolkit serves as an easy to use interface for the TMinuit 
// package: 
// 
// - It allows to fit a curve to one and two dimensional histograms 
//   (TH2F::Fit() only allows to fit a hyperplane).
// - Or n points can be specified directly via a n x 2 matrix.
// - An option for robust fitting of non-linear functions is implemented.
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
//    the command "void Fit()".
//
//
// 3. Accessing the fit results
//
//  The N parameters of the formula are stored in a N-dimensional vector which
//  is returned by "TVectorD * GetParameters()". In a similar the covariance 
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
//  The robust option becomes active if a weighting function is specified. 
//  All points are sorted according to their distance to the curve and
//  weighted. The weighting function must be defined on the interval [0,1].
//
//  Some standard weighting functions are predefined in 
//  "SetWeightFunction(Char_t * name, Float_t param1, Float_t param2 = 0)":
//   - "BOX" equals to 1 if x < param1 and to 0 if x > param1.
//   - "EXPONENTIAL" corresponds to "Math::Exp(-TMath::Log(param1)*x)"
//   - "ERRORFUNCTION" corresponds to "TMath::Erfc((x-param1)/param2)"
//
//-------------------------------------------------------------------------


ClassImp(AliTMinuitToolkit)

AliTMinuitToolkit::AliTMinuitToolkit() : 
   TNamed(),
   fFormula(0),
   fWeightFunction(0),
   fFitAlgorithm(0),
   fPoints(0),
   fParam(0),
   fParamLimits(0),
   fCovar(0),
   fChi2(0)
{
 //
 // standard constructor
 //
 fMaxCalls = 500;
 fPrecision = 1;
}



AliTMinuitToolkit::~AliTMinuitToolkit(){
  //
  // destructor
  //
  delete fPoints;
  delete fWeightFunction;
  delete fParamLimits;
  delete fFormula;
  delete fParam;
  delete fCovar;
  delete fChi2;
}

void AliTMinuitToolkit::FitHistogram(TH1F * his) {
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


void AliTMinuitToolkit::FitHistogram(TH2F * his) {
 //
 // Fit a two dimensional histogram
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


void AliTMinuitToolkit::SetWeightFunction(Char_t * name, Float_t param1, Float_t param2) {
 //
 // Set the weight function which must be defined on the interval [0,1].
 //
 TString FuncType(name);
 FuncType.ToUpper();
 
 if (FuncType == "EXPONENTIAL") fWeightFunction = new TFormula("exp", Form("TMath::Exp(-TMath::Log(%f)*x)", param1));
 if (FuncType == "BOX") fWeightFunction = new TFormula("box", Form("TMath::Erfc((x-%f)/0.0001)", param1));
 if (FuncType == "ERRORFUNCTION") fWeightFunction = new TFormula("err", Form("TMath::Erfc((x-%f)/%f)", param1, param2));// !!!!!!!!!!!!!!!!!
 
}


void AliTMinuitToolkit::FitterFCN(int &npar, double *dummy, double &fchisq, double *gin, int iflag){
  //
  // internal function which gives the specified function to the TMinuit function
  //  

  // suppress warnings for unused variables:
  dummy = dummy;
  iflag = iflag;
  npar = npar;
  //
  AliTMinuitToolkit * fitter = (AliTMinuitToolkit*)TVirtualFitter::GetFitter()->GetObjectFit();
  fchisq = 0;
  Int_t nvar       = fitter->GetPoints()->GetNcols()-1;
  Int_t npoints    = fitter->GetPoints()->GetNrows();
  
  // sort points for weighting
  Double_t sortList[npoints];
  Int_t indexList[npoints];
  
  TVectorD *fWeight = new TVectorD(npoints);
  
  for (Int_t ipoint=0; ipoint<npoints; ipoint++){
    Double_t x[100];
    for (Int_t ivar=0; ivar<nvar; ivar++){
      x[ivar] = (*fitter->GetPoints())(ipoint, ivar);      
    }
    Float_t funx   = fitter->GetFormula()->EvalPar(x,gin);
    sortList[ipoint] = TMath::Abs((*fitter->GetPoints())(ipoint, nvar) - funx);
  }
  
  TMath::Sort(npoints, sortList, indexList, false);

  Double_t t;
  for (Int_t ipoint=0; ipoint<npoints; ipoint++){
   t = indexList[ipoint]/(Double_t)npoints;
   (*fWeight)(ipoint) = fitter->GetWeightFunction()->Eval(t);
  }
  //
  // calculate chisquare
  for (Int_t ipoint=0; ipoint<npoints; ipoint++){
    Double_t x[100];
    for (Int_t ivar=0; ivar<nvar; ivar++){
      x[ivar] = (*fitter->GetPoints())(ipoint, ivar);      
    }
    Float_t funx   = fitter->GetFormula()->EvalPar(x,gin);
    
    Double_t delta = (*fitter->GetPoints())(ipoint, nvar) - funx;
    fchisq+= delta*delta*(*fWeight)(ipoint);

  }
  delete fWeight;
}
 

void AliTMinuitToolkit::Fit() {
  //
  // internal function that calls the fitter
  //
  Int_t nparam = fParam->GetNrows();
  
  // set all paramter limits to infinity as default
  if (fParamLimits == 0) {
   fParamLimits = new TMatrixD(nparam ,2);
   for (Int_t iparam=0; iparam<nparam; iparam++){
    (*fParamLimits)(iparam, 0) = 0;
    (*fParamLimits)(iparam, 1) = 0;
   }
  }
  
  // set all weights to 1 as default
  if (fWeightFunction == 0) {
   fWeightFunction = new TFormula("constant", "1");
  }
  
  // migrad fit algorithm as default
  if (fFitAlgorithm == 0) {
   fFitAlgorithm = "migrad";
  }
  
  // set up the fitter
  TVirtualFitter *minuit = TVirtualFitter::Fitter(0, nparam);
  minuit->SetObjectFit(this); 
  minuit->SetFCN((void*)(AliTMinuitToolkit::FitterFCN));
  
  // initialize paramters (step size???)
  for (Int_t iparam=0; iparam<nparam; iparam++){
   minuit->SetParameter(iparam, Form("p[%d]",iparam), (*fParam)(iparam), (*fParam)(iparam)/10, (*fParamLimits)(iparam, 0), (*fParamLimits)(iparam, 1));
  }
  
  Double_t argList[2];
  argList[0] = fMaxCalls; //maximal number of calls 
  argList[1] = fPrecision; //tolerance normalized to 0.001 
  if (fMaxCalls == 500 && fPrecision == 1) minuit->ExecuteCommand(fFitAlgorithm, 0, 0); 
  if (fMaxCalls != 500 || fPrecision != 1) minuit->ExecuteCommand(fFitAlgorithm, argList, 2);
  // two additional arguments can be specified ExecuteCommand("migrad", argList, 2) - use 0,0 for default

  // fill parameter vector
  for (Int_t ivar=0; ivar<nparam; ivar++){
   (*fParam)(ivar) = minuit->GetParameter(ivar);
  }
  
  // fill covariance matrix
  fCovar = new TMatrixD(nparam, nparam);
  //TVirtualFitter *fitCov = TVirtualFitter::GetFitter();
  for(Int_t i=0; i < nparam; i++) {
   for(Int_t j=0; j < nparam; j++) {
    (*fCovar)(i,j) = minuit->GetCovarianceMatrixElement(i,j);
   }
  }
 
}



void AliTMinuitToolkit::Test() {
 //
 // This test function shows the basic working principles of this class 
 // and illustrates how a robust fit can improve the results
 //
 TFormula *FormExp = new TFormula("formExp", "[0]*TMath::Exp(-[1]*x)");
 SetFitFunction(FormExp);
 SetFitAlgorithm("simplex");
 // Set initial values
 TVectorD *vec1 = new TVectorD(2);
 (*vec1)(0) = 1800;
 (*vec1)(1) = 1;
 SetInitialParam(vec1);
 //provide some example histogram
 TH1F * hist = new TH1F("bla", "with (red) and without (black) robust option", 20,0,4);
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
 // fit it with the exponential decay
 FitHistogram(hist);
 // draw fit function
 TF1 *func = new TF1("test", "[0]*TMath::Exp(-[1]*x)", 0, 6);
 func->SetParameter(0, (*GetParameters())(0));
 func->SetParameter(1, (*GetParameters())(1));
 func->Draw("same");
 // robust fit
 TVectorD *vec2 = new TVectorD(2);
 (*vec2)(0) = 1800;
 (*vec2)(1) = 1;
 SetInitialParam(vec2);
 SetWeightFunction("Box", 0.7);
 FitHistogram(hist);
 TF1 *func2 = new TF1("test2", "[0]*TMath::Exp(-[1]*x)", 0, 6);
 func2->SetParameter(0, (*GetParameters())(0));
 func2->SetParameter(1, (*GetParameters())(1));
 func2->SetLineColor(kRed);
 func2->Draw("same");
 
}

