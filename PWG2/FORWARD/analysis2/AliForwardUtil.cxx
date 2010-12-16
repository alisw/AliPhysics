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
UShort_t
AliForwardUtil::ParseCollisionSystem(const char* sys)
{
  TString s(sys);
  s.ToLower();
  if (s.Contains("p-p")   || s.Contains("pp"))   return AliForwardUtil::kPP; 
  if (s.Contains("pb-pb") || s.Contains("pbpb")) return AliForwardUtil::kPbPb;
  if (s.Contains("a-a")   || s.Contains("aa"))   return AliForwardUtil::kPbPb;
  return AliForwardUtil::kUnknown;
}
//____________________________________________________________________
const char*
AliForwardUtil::CollisionSystemString(UShort_t sys)
{
  switch (sys) { 
  case AliForwardUtil::kPP:   return "pp";
  case AliForwardUtil::kPbPb: return "PbPb";
  }
  return "unknown";
}
//____________________________________________________________________
UShort_t
AliForwardUtil::ParseCenterOfMassEnergy(UShort_t sys, Float_t v)
{
  Float_t energy = v; 
  if (sys != AliForwardUtil::kPP) energy = energy / 208 * 82;
  if (TMath::Abs(energy - 900.)   < 10)  return 900;
  if (TMath::Abs(energy - 2400.)  < 10)  return 2400;
  if (TMath::Abs(energy - 2750.)  < 10)  return 2750;
  if (TMath::Abs(energy - 5500.)  < 40)  return 5500;
  if (TMath::Abs(energy - 7000.)  < 10)  return 7000;
  if (TMath::Abs(energy - 10000.) < 10)  return 10000;
  if (TMath::Abs(energy - 14000.) < 10)  return 14000;
  return 0;
}
//____________________________________________________________________
const char* 
AliForwardUtil::CenterOfMassEnergyString(UShort_t cms)
{
  return Form("%04dGeV", cms);
}
//____________________________________________________________________
Short_t
AliForwardUtil::ParseMagneticField(Float_t v)
{
  if (TMath::Abs(v - 5.) < 1 ) return +5;
  if (TMath::Abs(v + 5.) < 1 ) return -5;
  if (TMath::Abs(v) < 1)       return 0;
  return 999;
}
//____________________________________________________________________
const char* 
AliForwardUtil::MagneticFieldString(Short_t f)
{
  return Form("%01dkG", f);
}


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
  /** 
   * Utility function to use in TF1 defintition 
   */
  Double_t landauGausI(Double_t* xp, Double_t* pp) 
  {
    Double_t x         = xp[0];
    Double_t constant  = pp[AliForwardUtil::ELossFitter::kC];
    Double_t delta     = pp[AliForwardUtil::ELossFitter::kDelta];
    Double_t xi        = pp[AliForwardUtil::ELossFitter::kXi];
    Double_t sigma     = pp[AliForwardUtil::ELossFitter::kSigma];
    Double_t sigma_n   = pp[AliForwardUtil::ELossFitter::kSigmaN];
    Int_t    i         = Int_t(pp[AliForwardUtil::ELossFitter::kN]);

    return constant * AliForwardUtil::ILandauGaus(x,delta,xi,sigma,sigma_n,i);
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
AliForwardUtil::ILandauGaus(Double_t x, Double_t delta, Double_t xi, 
			    Double_t sigma, Double_t sigma_n, Int_t i)
{
  Double_t delta_i =  (i == 1 ? delta : i * (delta + xi * TMath::Log(i)));
  Double_t xi_i    =  i * xi;
  Double_t sigma_i =  (i == 1 ? sigma : TMath::Sqrt(Double_t(i))*sigma);
  if (sigma_i < 1e-10) { 
    // Fall back to landau 
    return AliForwardUtil::Landau(x, delta_i, xi_i);
  }
  return AliForwardUtil::LandauGaus(x, delta_i, xi_i, sigma_i, sigma_n);
}

//____________________________________________________________________
Double_t 
AliForwardUtil::IdLandauGausdPar(Double_t x, 
				 UShort_t par,   Double_t dPar, 
				 Double_t delta, Double_t xi, 
				 Double_t sigma, Double_t sigma_n, 
				 Int_t    i)
{
  if (dPar == 0) return 0;
  Double_t dp      = dPar;
  Double_t d2      = dPar / 2;
  Double_t delta_i =  i * (delta + xi * TMath::Log(i));
  Double_t xi_i    =  i * xi;
  Double_t si      =  TMath::Sqrt(Double_t(i));
  Double_t sigma_i =  si*sigma;
  Double_t y1      = 0;
  Double_t y2      = 0;
  Double_t y3      = 0;
  Double_t y4      = 0;
  switch (par) {
  case 0: 
    y1 = ILandauGaus(x, delta_i+i*dp, xi_i, sigma_i, sigma_n, i);
    y2 = ILandauGaus(x, delta_i+i*d2, xi_i, sigma_i, sigma_n, i);
    y3 = ILandauGaus(x, delta_i-i*d2, xi_i, sigma_i, sigma_n, i);
    y4 = ILandauGaus(x, delta_i-i*dp, xi_i, sigma_i, sigma_n, i);
    break;
  case 1: 
    y1 = ILandauGaus(x, delta_i, xi_i+i*dp, sigma_i, sigma_n, i);
    y2 = ILandauGaus(x, delta_i, xi_i+i*d2, sigma_i, sigma_n, i);
    y3 = ILandauGaus(x, delta_i, xi_i-i*d2, sigma_i, sigma_n, i);
    y4 = ILandauGaus(x, delta_i, xi_i-i*dp, sigma_i, sigma_n, i);
    break;
  case 2: 
    y1 = ILandauGaus(x, delta_i, xi_i, sigma_i+si*dp, sigma_n, i);
    y2 = ILandauGaus(x, delta_i, xi_i, sigma_i+si*d2, sigma_n, i);
    y3 = ILandauGaus(x, delta_i, xi_i, sigma_i-si*d2, sigma_n, i);
    y4 = ILandauGaus(x, delta_i, xi_i, sigma_i-si*dp, sigma_n, i);
    break;
  case 3: 
    y1 = ILandauGaus(x, delta_i, xi_i, sigma_i, sigma_n+dp, i);
    y2 = ILandauGaus(x, delta_i, xi_i, sigma_i, sigma_n+d2, i);
    y3 = ILandauGaus(x, delta_i, xi_i, sigma_i, sigma_n-d2, i);
    y4 = ILandauGaus(x, delta_i, xi_i, sigma_i, sigma_n-dp, i);
    break;
  default:
    return 0;
  } 
  
  Double_t d0  = y1 - y4;
  Double_t d1  = 2 * (y2 - y3);
  
  Double_t g   = 1/(2*dp) * (4*d1 - d0) / 3;
   
  return g;
}

//____________________________________________________________________
Double_t 
AliForwardUtil::NLandauGaus(Double_t x, Double_t delta, Double_t xi, 
			    Double_t sigma, Double_t sigma_n, Int_t n, 
			    Double_t* a)
{
  Double_t result = ILandauGaus(x, delta, xi, sigma, sigma_n, 1);
  for (Int_t i = 2; i <= n; i++) 
    result += a[i-2] * AliForwardUtil::ILandauGaus(x,delta,xi,sigma,sigma_n,i);
  return result;
}
namespace { 
  const Int_t kColors[] = { kRed+1, 
			    kPink+3, 
			    kMagenta+2, 
			    kViolet+2, 
			    kBlue+1, 
			    kAzure+3, 
			    kCyan+1, 
			    kTeal+2, 
			    kGreen+2, 
			    kSpring+3, 
			    kYellow+2, 
			    kOrange+2 };
}

//____________________________________________________________________
TF1*
AliForwardUtil::MakeNLandauGaus(Double_t  c, 
				Double_t  delta, Double_t xi, 
				Double_t  sigma, Double_t sigma_n, Int_t n, 
				Double_t* a, 
				Double_t  xmin, Double_t xmax)
{
  Int_t npar       = AliForwardUtil::ELossFitter::kN+n;
  TF1* landaun     = new TF1(Form("nlandau%d", n), &landauGausN,xmin,xmax,npar);
  // landaun->SetLineStyle(((n-2) % 10)+2); // start at dashed
  landaun->SetLineColor(kColors[((n-1) % 12)]); // start at red
  landaun->SetLineWidth(2);
  landaun->SetNpx(500);
  landaun->SetParNames("C","#Delta_{p}","#xi", "#sigma", "#sigma_{n}", "N");

  // Set the initial parameters from the seed fit 
  landaun->SetParameter(AliForwardUtil::ELossFitter::kC,      c);       
  landaun->SetParameter(AliForwardUtil::ELossFitter::kDelta,  delta);   
  landaun->SetParameter(AliForwardUtil::ELossFitter::kXi,     xi);      
  landaun->SetParameter(AliForwardUtil::ELossFitter::kSigma,  sigma);   
  landaun->SetParameter(AliForwardUtil::ELossFitter::kSigmaN, sigma_n); 
  landaun->FixParameter(AliForwardUtil::ELossFitter::kN,      n);       

  // Set the range and name of the scale parameters 
  for (UShort_t i = 2; i <= n; i++) {// Take parameters from last fit 
    landaun->SetParameter(AliForwardUtil::ELossFitter::kA+i-2, a[i-2]);
    landaun->SetParName(AliForwardUtil::ELossFitter::kA+i-2, Form("a_{%d}", i));
  }
  return landaun;
}
//____________________________________________________________________
TF1*
AliForwardUtil::MakeILandauGaus(Double_t  c, 
				Double_t  delta, Double_t xi, 
				Double_t  sigma, Double_t sigma_n, Int_t i, 
				Double_t  xmin, Double_t xmax)
{
  Int_t npar       = AliForwardUtil::ELossFitter::kN+1;
  TF1* landaui     = new TF1(Form("ilandau%d", i), &landauGausI,xmin,xmax,npar);
  // landaui->SetLineStyle(((i-2) % 10)+2); // start at dashed
  landaui->SetLineColor(kColors[((i-1) % 12)]); // start at red
  landaui->SetLineWidth(1);
  landaui->SetNpx(500);
  landaui->SetParNames("C","#Delta_{p}","#xi", "#sigma", "#sigma_{n}", "i");

  // Set the initial parameters from the seed fit 
  landaui->SetParameter(AliForwardUtil::ELossFitter::kC,      c);       
  landaui->SetParameter(AliForwardUtil::ELossFitter::kDelta,  delta);   
  landaui->SetParameter(AliForwardUtil::ELossFitter::kXi,     xi);      
  landaui->SetParameter(AliForwardUtil::ELossFitter::kSigma,  sigma);   
  landaui->SetParameter(AliForwardUtil::ELossFitter::kSigmaN, sigma_n); 
  landaui->FixParameter(AliForwardUtil::ELossFitter::kN,      i);       

  return landaui;
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
  TF1* landau1 = new TF1("landau1", landauGaus1, minE,maxEE,kSigmaN+1);

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

  // Array of weights 
  TArrayD a(n-1);
  for (UShort_t i = 2; i <= n; i++) 
    a.fArray[i-2] = (n == 2 ? 0.05 : 0.000001);
  // Make the fit function 
  TF1* landaun     = MakeNLandauGaus(r->Parameter(kC),
				     r->Parameter(kDelta),
				     r->Parameter(kXi),
				     r->Parameter(kSigma),
				     r->Parameter(kSigmaN),
				     n,a.fArray,minE,maxEi);
  landaun->SetParLimits(kDelta,  minE, fMaxRange);       // Delta
  landaun->SetParLimits(kXi,     0.00, fMaxRange);       // xi
  landaun->SetParLimits(kSigma,  0.01, 1);            // sigma
  // Check if we're using the noise sigma 
  if (sigman <= 0)  landaun->FixParameter(kSigmaN, 0);
  else              landaun->SetParLimits(kSigmaN, 0, fMaxRange);

  // Set the range and name of the scale parameters 
  for (UShort_t i = 2; i <= n; i++) {// Take parameters from last fit 
    landaun->SetParLimits(kA+i-2, 0,1);
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
