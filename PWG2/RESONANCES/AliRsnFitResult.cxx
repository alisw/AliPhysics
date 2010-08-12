//
//
//

#include <TF1.h>
#include <TH1.h>
#include <TMath.h>

#include "AliLog.h"
#include <AliRsnFitResult.h>

ClassImp(AliRsnFitResult)

//__________________________________________________________________________________________________
AliRsnFitResult::AliRsnFitResult(const char *name, ESignalType sType, EBackgroundType bType) :
  TNamed(name, ""),
  fSignalType(sType),
  fBackgroundType(bType)
{
//
// Default constructor
//

  InitFunctions();
}
  
//__________________________________________________________________________________________________
Double_t AliRsnFitResult::NormSquareRoot(Double_t *x, Double_t *par)
{
//
// Computes a square-root like function normalized 
// to the integral in the background range specified
// for this object.
// This is obtained by dividing the function value by
// its integral analytically computed in that range.
// Returns 0 in all cases where a floating point exception
// could be raised by computing square root of a negative number.
//

  if (x[0] < par[1]) return 0.0;
  
  Double_t fcnVal = TMath::Sqrt(x[0] - par[1]);

  Double_t fcnMax = fPeakRange[1] - par[1];
  Double_t fcnMin = fPeakRange[0] - par[1];
  Double_t fcnInt = (2.0 / 3.0) * ((fcnMax*TMath::Sqrt(fcnMax) - fcnMin*TMath::Sqrt(fcnMin)));
  
  return par[0] * fcnVal / fcnInt;
}

//__________________________________________________________________________________________________
Double_t AliRsnFitResult::NormLinear(Double_t *x, Double_t *par)
{
//
// Computes a 1st order polynomial function normalized 
// to the integral in the background range specified
// for this object.
// This is obtained by dividing the function value by
// its integral analytically computed in that range.
// Returns 0 in all cases where a floating point exception
// could be raised by computing square root of a negative number.
//
  
  Double_t fcnVal = 1.0 + x[0] * par[1];
  Double_t fcnInt = (fPeakRange[1] - fPeakRange[0]) + 0.5*par[1]*(fPeakRange[1]*fPeakRange[1] - fPeakRange[0]*fPeakRange[0]);
  
  return par[0] * fcnVal / fcnInt;
}
  
//__________________________________________________________________________________________________
Double_t AliRsnFitResult::NormPoly2(Double_t *x, Double_t *par)
{
//
// Computes a 2nd order polynomial function normalized 
// to the integral in the background range specified
// for this object.
// This is obtained by dividing the function value by
// its integral analytically computed in that range.
// Returns 0 in all cases where a floating point exception
// could be raised by computing square root of a negative number.
//
  
  Double_t fcnVal = 1.0 + x[0] * par[1] + x[0]*x[0] * par[2];
  Double_t fcnInt = (fPeakRange[1] - fPeakRange[0]) + par[1]*(fPeakRange[1]*fPeakRange[1] - fPeakRange[0]*fPeakRange[0])/2.0 + par[2]*(fPeakRange[1]*fPeakRange[1]*fPeakRange[1] - fPeakRange[0]*fPeakRange[0]*fPeakRange[0])/3.0;
  
  return par[0] * fcnVal / fcnInt;
}

//__________________________________________________________________________________________________
Double_t AliRsnFitResult::NormBreitWigner(Double_t *x, Double_t *par)
{
//
// Computes a Breit-Wigner function normalized 
// to the integral in the background range specified
// for this object.
// This is obtained by dividing the function value by
// its integral analytically computed in that range.
// Returns 0 in all cases where a floating point exception
// could be raised by computing square root of a negative number.
//
  
  Double_t fcnVal = 1.0 / ((x[0] - par[1])*(x[0] - par[1]) + 0.25*par[2]*par[2]);
  Double_t fcnInt = 2.0 / par[2] * (TMath::ATan((fPeakRange[1] - par[1]) / (0.5 * par[2])) - TMath::ATan((fPeakRange[0] - par[1]) / (0.5 * par[2])));
  
  return par[0] * fcnVal / fcnInt;
}

//__________________________________________________________________________________________________
Double_t AliRsnFitResult::NormGaus(Double_t *x, Double_t *par)
{
//
// Computes a Gaussian function normalized 
// to the integral in the background range specified
// for this object.
// This is obtained by dividing the function value by
// its integral analytically computed in that range.
// Returns 0 in all cases where a floating point exception
// could be raised by computing square root of a negative number.
//
  
  Double_t fcnVal      = TMath::Gaus(x[0], par[1], par[2], kTRUE);
  Double_t fcnIntLeft  = 0.5 * TMath::Erf(((par[1] - fPeakRange[0]) / par[2]) / TMath::Sqrt(2.0));
  Double_t fcnIntRight = 0.5 * TMath::Erf(((fPeakRange[1] - par[1]) / par[2]) / TMath::Sqrt(2.0));
  Double_t fcnInt      = fcnIntLeft + fcnIntRight;
  
  return par[0] * fcnVal / fcnInt;
}

//__________________________________________________________________________________________________
Double_t AliRsnFitResult::Signal(Double_t *x, Double_t *par)
{
//
// Computes the signal function according to chosen option
//

  switch (fSignalType)
  {
    case kBreitWigner:
      return NormBreitWigner(x, par);
    case kGaus:
      return NormGaus(x, par);
    default:
      AliError("Invalid signal function");
      return 0.0;
  }
}
  
//__________________________________________________________________________________________________
Double_t AliRsnFitResult::Background(Double_t *x, Double_t *par)
{
//
// Computes the background function according to chosen option
//

  switch (fBackgroundType)
  {
    case kSquareRoot:
      return NormSquareRoot(x, par);
    case kLinear:
      return NormLinear(x, par);
    case kPoly2:
      return NormPoly2(x, par);
    default:
      AliError("Invalid background function");
      return 0.0;
  }
}

//__________________________________________________________________________________________________
Double_t AliRsnFitResult::Sum(Double_t *x, Double_t *par)
{
//
// Computes the sum of the signal and the background chosen.
// First 3 parameters are for the signal (they are always 3),
// and other parameters, up to the necessary amount, are for
// the background.
//

  Double_t signal     = Signal(x, par);
  Double_t background = Background(x, &par[3]);

  return signal + background;
}

//__________________________________________________________________________________________________
Int_t AliRsnFitResult::GetNParBackground()
{
//
// Tells how many parameters has the background function
//

  switch (fBackgroundType)
  {
    case kSquareRoot: return 2;
    case kLinear    : return 2;
    case kPoly2     : return 3;
    default         : return 0;
  }
}

//__________________________________________________________________________________________________
Bool_t AliRsnFitResult::InitFunctions()
{
//
// Initialize functions
//

  //if (!fSignalFcn) fSignalFcn = new TF1(Form("sg%s", GetName()), Signal, fFullRange[0], fFullRange[1], 3);
  //if (!fBackgroundFcn) fBackgroundFcn = new TF1(Form("bg%s", GetName()), Background, fFullRange[0], fFullRange[1], GetNParBackground());
  //if (!fSumFcn) fSumFcn = new TF1(Form("sum%s", GetName()), Sum, fFullRange[0], fFullRange[1], 3+GetNParBackground());
}

//__________________________________________________________________________________________________
void AliRsnFitResult::SetHistogram(TH1F *histogram)
{
//
// Clones the passed histogram to proceed with fits
//

  if (fHistogram) delete fHistogram;
  
  fHistogram = (TH1D*)histogram->Clone();
  fHistogram->SetName(Form("h_%s_%s", GetName(), histogram->GetName()));
  fHistogram->SetTitle(histogram->GetTitle());
}

//__________________________________________________________________________________________________
Bool_t AliRsnFitResult::SingleFit(const char *opt, Double_t mass, Double_t width)
{
//
// Executes a single fit of the function
//
/*
  // require histogram
  if (!fHistogram)
  {
    AliError("Require an initialized histogram!");
    return kFALSE;
  }
  
  // if necessary, initialize functions
  if (!fSignalFcn || !fBackgroundFcn || !fSumFcn) InitFunctions();
  
  // step #0: know how many parameter we have
  Int_t npar = GetNParBackground();
  
  // step #1:
  // fit outside peak to roughly initialize the background
  for (Int_t i = 0; i < npar; i++) fBackgroundFcn->SetParameter(i, 0.5);
  status = (Int_t)h->Fit(fcnOut, opt);
  if (status) return kFALSE;

  // step #2:
  // estimate signal/background in the peak
  Int_t    imax = h->GetMaximumBin();
  Double_t xmax = h->GetBinCenter(imax);
  Double_t ymax = h->GetMaximum() - fBackgroundFcn->Eval(xmax);

  // step #3:
  // fit whole range with signal + background
  // which is initialized to the reference values for mass and width
  fSignalFcn->SetParameter(0, ymax);
  fSignalFcn->SetParameter(1, mass);
  fSignalFcn->SetParameter(2, width);
  for (Int_t i = 0; i < npar; i++) fcnSum->SetParameter(i + 3, fBackgroundFcn->GetParameter(i));
  status = (Int_t)h->Fit(fcnSum, opt);
  if (status) return kFALSE;
  
  // get integral and its error here
  fValue[kSumIntegralError] = fSumFcn->IntegralError(fitPeakRangeMin, fitPeakRangeMax);
  fValue[kSumIntegral]      = fSumFcn->Integral     (fitPeakRangeMin, fitPeakRangeMax);
  fValue[kChi2]             = fSumFcn->GetChisquare();
  fValue[kNDF]              = fSumFcn->GetNDF();
*/  
  return kTRUE;
}
