#include "AliForwardUtil.h"
#include <AliAnalysisManager.h>
#include "AliAODForwardMult.h"
#include <AliLog.h>
#include <AliInputEventHandler.h>
#include <AliESDEvent.h>
#include <AliPhysicsSelection.h>
#include <AliTriggerAnalysis.h>
#include <AliMultiplicity.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TMath.h>
#include <TError.h>

//====================================================================
Int_t    AliForwardUtil::fgConvolutionSteps  = 100;
Double_t AliForwardUtil::fgConvolutionNSigma = 5;
namespace {
  /** 
   * The shift of the most probable value for the ROOT function TMath::Landau 
   */
  const Double_t  mpshift  = -0.22278298;
  /** 
   * Integration normalisation 
   */
  const Double_t  invSq2pi = 1. / TMath::Sqrt(2*TMath::Pi());

  /** 
   * Utility function to use in TF1 defintition 
   */
  Double_t landauGaus1(Double_t* xp, Double_t* pp) 
  {
    Double_t x        = xp[0];
    Double_t constant = pp[0];
    Double_t delta    = pp[1];
    Double_t xi       = pp[2];
    Double_t sigma    = pp[3];
    Double_t sigma_n  = pp[4];

    return constant * AliForwardUtil::LandauGaus(x, delta, xi, sigma, sigma_n);
  }

  /** 
   * Utility function to use in TF1 defintition 
   */
  Double_t landauGausN(Double_t* xp, Double_t* pp) 
  {
    Double_t  x        = xp[0];
    Double_t  constant = pp[0];
    Double_t  delta    = pp[1];
    Double_t  xi       = pp[2];
    Double_t  sigma    = pp[3];
    Double_t  sigma_n  = pp[4];
    Int_t     n        = Int_t(pp[5]);
    Double_t* a        = &(pp[6]);

    return constant * AliForwardUtil::NLandauGaus(x, delta, xi, sigma, sigma_n,
						  n, a);
  }


}
//____________________________________________________________________
Double_t 
AliForwardUtil::Landau(Double_t x, Double_t delta, Double_t xi)
{
  return TMath::Landau(x, delta - xi * mpshift, xi);
}
//____________________________________________________________________
Double_t 
AliForwardUtil::LandauGaus(Double_t x, Double_t delta, Double_t xi,
			   Double_t sigma, Double_t sigma_n)
{
  Double_t deltap = delta - xi * mpshift;
  Double_t sigma2 = sigma_n*sigma_n + sigma*sigma;
  Double_t sigma1 = TMath::Sqrt(sigma2);
  Double_t xlow   = x - fgConvolutionNSigma * sigma1;
  Double_t xhigh  = x - fgConvolutionNSigma * sigma1;
  Double_t step   = (xhigh - xlow) / fgConvolutionSteps;
  Double_t sum    = 0;
  
  for (Int_t i = 0; i <= fgConvolutionSteps/2; i++) { 
    Double_t x1 = xlow + (i - .5) * step;
    Double_t x2 = xlow - (i - .5) * step;
    
    sum += TMath::Landau(x1, deltap, xi, kTRUE) * TMath::Gaus(x, x1, sigma1);
    sum += TMath::Landau(x2, deltap, xi, kTRUE) * TMath::Gaus(x, x2, sigma1);
  }
  return step * sum * invSq2pi / sigma1;
}

//____________________________________________________________________
Double_t 
AliForwardUtil::NLandauGaus(Double_t x, Double_t delta, Double_t xi, 
			    Double_t sigma, Double_t sigma_n, Int_t n, 
			    Double_t* a)
{
  Double_t result = LandauGaus(x, delta, xi, sigma, sigma_n);
  for (Int_t i = 2; i <= n; i++) { 
    Double_t delta_i =  i * (delta + xi * TMath::Log(i));
    Double_t xi_i    =  i * xi;
    Double_t sigma_i =  TMath::Sqrt(Double_t(n)*sigma);
    Double_t a_i     =  a[i];
    result           += a_i * AliForwardUtil::LandauGaus(x, delta_i, xi_i, 
							 sigma_i, sigma_n);
  }
  return result;
}

//====================================================================
AliForwardUtil::ELossFitter::ELossFitter(Double_t lowCut, 
					 Double_t maxRange, 
					 UShort_t minusBins) 
  : fLowCut(lowCut), fMaxRange(maxRange), fMinusBins(minusBins), 
    fFitResults(0), fFunctions(0)
{
  fFitResults.SetOwner();
  fFunctions.SetOwner();
}
//____________________________________________________________________
AliForwardUtil::ELossFitter::~ELossFitter()
{
  fFitResults.Delete();
  fFunctions.Delete();
}
//____________________________________________________________________
void
AliForwardUtil::ELossFitter::Clear()
{
  fFitResults.Clear();
  fFunctions.Clear();
}
//____________________________________________________________________
TF1*
AliForwardUtil::ELossFitter::Fit1Particle(TH1* dist, Double_t sigman)
{
  // Clear the cache 
  Clear();
  
  // Find the fit range 
  dist->GetXaxis()->SetRangeUser(fLowCut, fMaxRange);
  
  // Normalize peak to 1 
  Double_t max = dist->GetMaximum(); 
  dist->Scale(1/max);
  
  // Get the bin with maximum 
  Int_t    maxBin = dist->GetMaximumBin();
  Double_t maxE   = dist->GetBinLowEdge(maxBin);
  
  // Get the low edge 
  dist->GetXaxis()->SetRangeUser(fLowCut, maxE);
  Int_t    minBin = maxBin - fMinusBins; // dist->GetMinimumBin();
  Double_t minE   = TMath::Max(dist->GetBinCenter(minBin),fLowCut);
  Double_t maxEE  = dist->GetBinCenter(maxBin+2*fMinusBins);

  // Restore the range 
  dist->GetXaxis()->SetRangeUser(0, fMaxRange);
  
  // Define the function to fit 
  TF1*          landau1 = new TF1("landau1", landauGaus1, minE, maxEE, 5);

  // Set initial guesses, parameter names, and limits  
  landau1->SetParameters(5,.5,.07,1,sigman);
  landau1->SetParNames("C","#Delta_{p}","#xi", "#sigma", "#sigma_{n}");
  landau1->SetParLimits(1, minE, fMaxRange);
  landau1->SetParLimits(2, 0.00, fMaxRange);
  landau1->SetParLimits(3, 0.01, fMaxRange);
  if (sigman <= 0)  landau1->FixParameter(4, 0);
  else              landau1->SetParLimits(4, 0, fMaxRange);

  // Do the fit, getting the result object 
  TFitResultPtr r = dist->Fit(landau1, "RNQS", "", minE, maxEE);

  fFitResults.AddAtAndExpand(new TFitResult(*r), 0);
  fFunctions.AddAtAndExpand(landau1, 0);

  return landau1;
}
//____________________________________________________________________
TF1*
AliForwardUtil::ELossFitter::FitNParticle(TH1* dist, UShort_t n, 
					  Double_t sigman)
{
  // Get the seed fit result 
  TFitResult* r = static_cast<TFitResult*>(fFitResults.At(0));
  TF1*        f = static_cast<TF1*>(fFunctions.At(0));
  if (!r || !f) { 
    f = Fit1Particle(dist, sigman);
    r = static_cast<TFitResult*>(fFitResults.At(0));
    if (!r || !f) { 
      ::Warning("FitNLandau", "No first shot at landau fit");
      return 0;
    }
  }

  // Get some parameters from seed fit 
  Double_t delta1  = r->Parameter(1);
  Double_t xi1     = r->Parameter(2);
  Double_t maxEi   = n * (delta1 + xi1 * TMath::Log(n)) + 2 * n * xi1;
  Double_t minE    = f->GetXmin();

  // Make the fit function 
  TF1* landaun     = new TF1(Form("landau%d", n), &landauGausN,minE,maxEi,5+n);
  landaun->SetLineStyle((n % 10)+1);
  landaun->SetLineWidth(2);
  landaun->SetParNames("C","#Delta_{p}","#xi", "#sigma", "#sigma_{n}", "N");

  // Set the initial parameters from the seed fit 
  landaun->SetParameter(0, r->Parameter(0)); // Constant
  landaun->SetParameter(1, r->Parameter(1)); // Delta 
  landaun->SetParameter(2, r->Parameter(2)); // xi
  landaun->SetParameter(3, r->Parameter(3)); // sigma
  landaun->SetParameter(4, r->Parameter(4)); // sigma_n
  landaun->SetParLimits(1, minE, fMaxRange); // Delta
  landaun->SetParLimits(2, 0.00, fMaxRange); // xi
  landaun->SetParLimits(3, 0.01, fMaxRange); // sigma
  landaun->SetParLimits(4, 0.00, fMaxRange); // sigma_n
  // Fix the number parameter 
  landaun->FixParameter(5, n);               // N

  // Set the range and name of the scale parameters 
  for (UShort_t i = 2; i <= n; i++) {// Take parameters from last fit 
    landaun->SetParameter(5+i-2, n == 2 ? 0.05 : 0);
    landaun->SetParLimits(5+i-2, 0,1);
    landaun->SetParName(5+i-2, Form("a_{%d}", i));
  }

  // Do the fit 
  TFitResultPtr tr = dist->Fit(landaun, "RSQN", "", minE, maxEi);
  
  fFitResults.AddAtAndExpand(new TFitResult(*tr), n-1);
  fFunctions.AddAtAndExpand(landaun, n-1);
  
  return landaun;
}  

//====================================================================
AliForwardUtil::Histos::~Histos()
{
  if (fFMD1i) delete fFMD1i;
  if (fFMD2i) delete fFMD2i;
  if (fFMD2o) delete fFMD2o;
  if (fFMD3i) delete fFMD3i;
  if (fFMD3o) delete fFMD3o;
}

//____________________________________________________________________
TH2D*
AliForwardUtil::Histos::Make(UShort_t d, Char_t r, 
			     const TAxis& etaAxis) const
{
  Int_t ns = (r == 'I' || r == 'i') ? 20 : 40;
  TH2D* hist = new TH2D(Form("FMD%d%c_cache", d, r), 
			Form("FMD%d%c cache", d, r),
			etaAxis.GetNbins(), etaAxis.GetXmin(), 
			etaAxis.GetXmax(), ns, 0, 2*TMath::Pi());
  hist->SetXTitle("#eta");
  hist->SetYTitle("#phi [radians]");
  hist->SetZTitle("d^{2}N_{ch}/d#etad#phi");
  hist->Sumw2();
  hist->SetDirectory(0);

  return hist;
}
//____________________________________________________________________
void
AliForwardUtil::Histos::Init(const TAxis& etaAxis)
{
  fFMD1i = Make(1, 'I', etaAxis);
  fFMD2i = Make(2, 'I', etaAxis);
  fFMD2o = Make(2, 'O', etaAxis);
  fFMD3i = Make(3, 'I', etaAxis);
  fFMD3o = Make(3, 'O', etaAxis);
}
//____________________________________________________________________
void
AliForwardUtil::Histos::Clear(Option_t* option)
{
  fFMD1i->Reset(option);
  fFMD2i->Reset(option);
  fFMD2o->Reset(option);
  fFMD3i->Reset(option);
  fFMD3o->Reset(option);
}

//____________________________________________________________________
TH2D*
AliForwardUtil::Histos::Get(UShort_t d, Char_t r) const
{
  switch (d) { 
  case 1: return fFMD1i;
  case 2: return (r == 'I' || r == 'i' ? fFMD2i : fFMD2o);
  case 3: return (r == 'I' || r == 'i' ? fFMD3i : fFMD3o);
  }
  return 0;
}
//====================================================================
TList*
AliForwardUtil::RingHistos::DefineOutputList(TList* d) const
{
  if (!d) return 0;
  TList* list = new TList;
  list->SetName(fName.Data());
  d->Add(list);
  return list;
}
//____________________________________________________________________
TList*
AliForwardUtil::RingHistos::GetOutputList(TList* d) const
{
  if (!d) return 0;
  TList* list = static_cast<TList*>(d->FindObject(fName.Data()));
  return list;
}

//____________________________________________________________________
TH1*
AliForwardUtil::RingHistos::GetOutputHist(TList* d, const char* name) const
{
  return static_cast<TH1*>(d->FindObject(name));
}

//
// EOF
//
