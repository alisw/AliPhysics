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
    Double_t constant = pp[AliForwardUtil::ELossFitter::kC];
    Double_t delta    = pp[AliForwardUtil::ELossFitter::kDelta];
    Double_t xi       = pp[AliForwardUtil::ELossFitter::kXi];
    Double_t sigma    = pp[AliForwardUtil::ELossFitter::kSigma];
    Double_t sigma_n  = pp[AliForwardUtil::ELossFitter::kSigmaN];

    return constant * AliForwardUtil::LandauGaus(x, delta, xi, sigma, sigma_n);
  }

  /** 
   * Utility function to use in TF1 defintition 
   */
  Double_t landauGausN(Double_t* xp, Double_t* pp) 
  {
    Double_t  x        = xp[0];
    Double_t constant  = pp[AliForwardUtil::ELossFitter::kC];
    Double_t delta     = pp[AliForwardUtil::ELossFitter::kDelta];
    Double_t xi        = pp[AliForwardUtil::ELossFitter::kXi];
    Double_t sigma     = pp[AliForwardUtil::ELossFitter::kSigma];
    Double_t sigma_n   = pp[AliForwardUtil::ELossFitter::kSigmaN];
    Int_t     n        = Int_t(pp[AliForwardUtil::ELossFitter::kN]);
    Double_t* a        = &(pp[AliForwardUtil::ELossFitter::kA]);

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
  Double_t sigma1 = sigma_n == 0 ? sigma : TMath::Sqrt(sigma2);
  Double_t xlow   = x - fgConvolutionNSigma * sigma1;
  Double_t xhigh  = x + fgConvolutionNSigma * sigma1;
  Double_t step   = (xhigh - xlow) / fgConvolutionSteps;
  Double_t sum    = 0;
  
  for (Int_t i = 0; i <= fgConvolutionSteps/2; i++) { 
    Double_t x1 = xlow  + (i - .5) * step;
    Double_t x2 = xhigh - (i - .5) * step;
    
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
    Double_t sigma_i =  TMath::Sqrt(Double_t(n))*sigma;
    Double_t a_i     =  a[i-2];
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
  TF1*          landau1 = new TF1("landau1", landauGaus1, minE,maxEE,kSigmaN+1);

  // Set initial guesses, parameter names, and limits  
  landau1->SetParameters(1,0.5,0.07,0.1,sigman);
  landau1->SetParNames("C","#Delta_{p}","#xi", "#sigma", "#sigma_{n}");
  landau1->SetNpx(500);
  landau1->SetParLimits(kDelta, minE, fMaxRange);
  landau1->SetParLimits(kXi,    0.00, fMaxRange);
  landau1->SetParLimits(kSigma, 0.01, 0.1);
  if (sigman <= 0)  landau1->FixParameter(kSigmaN, 0);
  else              landau1->SetParLimits(kSigmaN, 0, fMaxRange);

  // Do the fit, getting the result object 
  TFitResultPtr r = dist->Fit(landau1, "RNQS", "", minE, maxEE);
  landau1->SetRange(minE, fMaxRange);
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
  Double_t delta1  = r->Parameter(kDelta);
  Double_t xi1     = r->Parameter(kXi);
  Double_t maxEi   = n * (delta1 + xi1 * TMath::Log(n)) + 2 * n * xi1;
  Double_t minE    = f->GetXmin();

  // Make the fit function 
  TF1* landaun     = new TF1(Form("landau%d", n), &landauGausN,minE,maxEi,kN+n);
  landaun->SetLineStyle(((n-2) % 10)+2); // start at dashed
  landaun->SetLineColor(((n-2) % 10)+2); // start at red
  landaun->SetLineWidth(1);
  landaun->SetNpx(500);
  landaun->SetParNames("C","#Delta_{p}","#xi", "#sigma", "#sigma_{n}", "N");

  // Set the initial parameters from the seed fit 
  landaun->SetParameter(kC,      r->Parameter(kC));      // Constant
  landaun->SetParameter(kDelta,  r->Parameter(kDelta));  // Delta 
  landaun->SetParameter(kXi,     r->Parameter(kXi));     // xi
  landaun->SetParameter(kSigma,  r->Parameter(kSigma));  // sigma
  landaun->SetParameter(kSigmaN, r->Parameter(kSigmaN)); // sigma_n
  landaun->SetParLimits(kDelta,  minE, fMaxRange);       // Delta
  landaun->SetParLimits(kXi,     0.00, fMaxRange);       // xi
  landaun->SetParLimits(kSigma,  0.01, 1);            // sigma
  // Check if we're using the noise sigma 
  if (sigman <= 0)  landaun->FixParameter(kSigmaN, 0);
  else              landaun->SetParLimits(kSigmaN, 0, fMaxRange);
  // Fix the number parameter 
  landaun->FixParameter(kN, n);               // N

  // Set the range and name of the scale parameters 
  for (UShort_t i = 2; i <= n; i++) {// Take parameters from last fit 
    landaun->SetParameter(kA+i-2, n == 2 ? 0.05 : 0.000001);
    landaun->SetParLimits(kA+i-2, 0,1);
    landaun->SetParName(kA+i-2, Form("a_{%d}", i));
  }

  // Do the fit 
  TFitResultPtr tr = dist->Fit(landaun, "RSQN", "", minE, maxEi);
  
  landaun->SetRange(minE, fMaxRange);
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
