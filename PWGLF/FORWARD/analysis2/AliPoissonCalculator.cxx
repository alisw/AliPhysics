#include "AliPoissonCalculator.h"
#include "AliForwardCorrectionManager.h"
#include <TH2D.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <TMath.h>
#include <TList.h>
#include <iostream>
#include <TAxis.h>

// 
// A class to calculate the multiplicity in @f$(x,y)@f$ bins
// using Poisson statistics. 
//
// The input is assumed to be binned in @f$(x,y)@f$ as described by
// the 2D histogram passwd to the Reset member function.
//
// The data is grouped in to regions as defined by the parameters
// fXLumping and fYLumping.  The total number of cells and number of
// empty cells is then calculate in each region.  The mean
// multiplicity over the region is then determined as
//
// @f[
// \lange m\rangle = -\log\left(\frac{e}{t}\right)
// @f]
// where @f$ e@f$ is the number of empty cells and @f$t@f$ is the
// total number of cells in the region.  A correction for counting
// statistics, is then applied 
// @f{eqnarray*}
//    c &=& \frac{1}{1 - \exp{-\lange m\rangle}}
//      &=& \frac{1}{1 - \frac{e}{t}}
// @f{eqnarray*}
// and the final number in each cell is then 
// @f[
//   h_i c \lange m\rangle 
// @f] 
// where @f$h_i@f$ is the number of hits in the cell @f$i@f$ 
// 
//

namespace {
  const char* kBasicN = "basic";
  const char* kEmptyN = "empty";
  const char* kTotalN = "total";
  const char* kBasicT = "Basic number of hits";
  const char* kEmptyT = "Empty number of bins/region";
  const char* kTotalT = "Total number of bins/region";
}




//____________________________________________________________________
AliPoissonCalculator::AliPoissonCalculator()
  : TNamed(),
    fXLumping(32), 
    fYLumping(4), 
    fTotal(0), 
    fEmpty(0), 
    fBasic(0),
    fEmptyVsTotal(0),
    fMean(0), 
    fOcc(0),
    fCorr(0)
{
  //
  // CTOR
  // 
}

//____________________________________________________________________
AliPoissonCalculator::AliPoissonCalculator(const char*)
  : TNamed("poissonCalculator", "Calculate N_ch using Poisson stat"),
    fXLumping(32), 
    fYLumping(4), 
    fTotal(0), 
    fEmpty(0), 
    fBasic(0),
    fEmptyVsTotal(0),
    fMean(0), 
    fOcc(0),
    fCorr(0)
{
  //
  // CTOR
  //
}
//____________________________________________________________________
AliPoissonCalculator::AliPoissonCalculator(const AliPoissonCalculator& o)
  : TNamed(o),
    fXLumping(o.fXLumping),
    fYLumping(o.fYLumping),
    fTotal(0), 
    fEmpty(0),
    fBasic(0), 
    fEmptyVsTotal(0),
    fMean(0), 
    fOcc(0),
    fCorr(0)
{
  Init();
  Reset(o.fBasic);
}

//____________________________________________________________________
AliPoissonCalculator::~AliPoissonCalculator()
{
  // CleanUp();
}

//____________________________________________________________________
void
AliPoissonCalculator::CleanUp()
{
  if (fTotal)        { delete fTotal;        fTotal        = 0; }
  if (fEmpty)        { delete fEmpty;        fEmpty        = 0; }
  if (fBasic)        { delete fBasic;        fBasic        = 0; }
  if (fEmptyVsTotal) { delete fEmptyVsTotal; fEmptyVsTotal = 0; }
  if (fMean)         { delete fMean;         fMean         = 0; }
  if (fOcc)          { delete fOcc;          fOcc          = 0; } 
  if (fCorr)         { delete fCorr;         fCorr         = 0; }
}
//____________________________________________________________________
AliPoissonCalculator&
AliPoissonCalculator::operator=(const AliPoissonCalculator& o)
{
  if (&o == this) return *this;
  TNamed::operator=(o);
  fXLumping = o.fXLumping;
  fYLumping = o.fYLumping;
  CleanUp();
  Init();
  Reset(o.fBasic);
  return *this;
}

//____________________________________________________________________
void
AliPoissonCalculator::Init(Int_t xLumping, Int_t yLumping)
{
  // 
  // Initialize 
  // 
  if (xLumping > 0) fXLumping = xLumping;
  if (yLumping > 0) fYLumping = yLumping;

  //Create diagnostics if void
  if (fEmptyVsTotal) return;
  
  MakeOutput();  
}
//____________________________________________________________________
void
AliPoissonCalculator::Define(const TAxis& xaxis, const TAxis& yaxis)
{
  // 
  // Initialize 
  // 
  const Double_t* xBins = xaxis.GetXbins()->GetArray();
  const Double_t* yBins = yaxis.GetXbins()->GetArray();
  Int_t           nX    = xaxis.GetNbins();
  Int_t           nY    = yaxis.GetNbins();
  Double_t        lX    = xaxis.GetXmin();
  Double_t        hX    = xaxis.GetXmax();
  Double_t        lY    = yaxis.GetXmin();
  Double_t        hY    = yaxis.GetXmax();
  
  if (xBins) { 
    if (yBins) fBasic = new TH2D(kBasicN, kBasicT, nX, xBins, nY, yBins);
    else       fBasic = new TH2D(kBasicN, kBasicT, nX, xBins, nY, lY, hY);
  }
  else { 
    if (yBins) fBasic = new TH2D(kBasicN, kBasicT, nX, lX, hX, nY, yBins);
    else       fBasic = new TH2D(kBasicN, kBasicT, nX, lX, hX, nY, lY, hY);
  }
  fBasic->SetXTitle(xaxis.GetTitle());
  fBasic->SetYTitle(yaxis.GetTitle());

  Reset(fBasic);
}
//____________________________________________________________________
void AliPoissonCalculator::MakeOutput() {

  Int_t n = fXLumping * fYLumping + 1;
  fEmptyVsTotal = new TH2D("emptyVsTotal", 
			   "# of empty # bins vs total # bins", 
			   n, -.5, n-.5, n, -.5, n-.5);
  fEmptyVsTotal->SetXTitle("Total # bins");
  fEmptyVsTotal->SetYTitle("# empty bins");
  fEmptyVsTotal->SetZTitle("Correlation");
  fEmptyVsTotal->SetOption("colz");
  fEmptyVsTotal->SetDirectory(0);

  n = (fXLumping + fYLumping);
  fMean = new TH1D("poissonMean", "Mean N_{ch} as calculated by Poisson",
		   10*n+1, -.1, n+.1);
  fMean->SetXTitle("#bar{N_{ch}}=log(empty/total)");
  fMean->SetYTitle("Events");
  fMean->SetFillColor(kRed+1);
  fMean->SetFillStyle(3001);
  fMean->SetLineColor(kBlack);
  fMean->SetDirectory(0);

  fOcc = new TH1D("occupancy", "Occupancy = #int_{1}^{#infty}dN P(N)",
		  1000, 0, 100);
  fOcc->SetXTitle("#int_{1}^{#infty}dN P(N) [%]");
  fOcc->SetYTitle("Events");
  fOcc->SetFillColor(kBlue+1);
  fOcc->SetFillStyle(3001);
  fOcc->SetLineColor(kBlack);
  fOcc->SetDirectory(0);

  fCorr = new TH2D("correction", "Correction as function of mean N_{ch}", 
		   10*n+1, -.1, n+.1, 100, 0, 10);
  fCorr->SetXTitle("#bar{N_{ch}}");
  fCorr->SetYTitle("Correction 1/(1-e^{#bar{N_{c}}})");
  fCorr->SetZTitle("Events");
  fCorr->SetOption("colz");
  fCorr->SetDirectory(0);
  
 
}
//____________________________________________________________________
void
AliPoissonCalculator::Output(TList* d)
{
  if (!d) return;
  if (!fEmptyVsTotal) MakeOutput();
  d->Add(fEmptyVsTotal);
  d->Add(fMean);
  d->Add(fOcc);
  d->Add(fCorr);
}

//____________________________________________________________________
void 
AliPoissonCalculator::SetLumping(UShort_t nx, UShort_t ny) 
{ 
  if (nx == fXLumping && ny == fYLumping && 
      fEmptyVsTotal && fTotal) 
    // Check if we have something to do. 
    return;
  CleanUp(); 
  Init(nx, ny);
}

//____________________________________________________________________
Int_t
AliPoissonCalculator::CheckLumping(char which, Int_t nBins, Int_t lumping) const
{
  if ((nBins % lumping) == 0) return lumping;
  Int_t l = lumping;
  do { 
    l--;
  } while (l > 0 && ((nBins % l) != 0));
  Warning("CheckLumping", "%c lumping %d is not a divisor of %d, set to %d", 
	  which, lumping, nBins, l);
  return l;
}

//____________________________________________________________________
void
AliPoissonCalculator::Reset(const TH2D* base)
{
  // 
  // Reset histogram 
  // 
  if (fBasic && fTotal && fEmpty) {
    fBasic->Reset();
    fTotal->Reset();
    fEmpty->Reset();
    return;
  }

  if (!base) return;

  Int_t    nXF   = base->GetNbinsX();
  Double_t xMin  = base->GetXaxis()->GetXmin();
  Double_t xMax  = base->GetXaxis()->GetXmax();
  Int_t    nYF   = base->GetNbinsY();
  Double_t yMin  = base->GetYaxis()->GetXmin();
  Double_t yMax  = base->GetYaxis()->GetXmax();

  fXLumping = CheckLumping('X', nXF, fXLumping);
  fYLumping = CheckLumping('Y', nYF, fYLumping);
  
  Int_t    nY    = nYF / fYLumping;
  Int_t    nX    = nXF / fXLumping;

  if (fBasic != base) { 
    fBasic = static_cast<TH2D*>(base->Clone(kBasicN));
    fBasic->SetTitle(kBasicT);
    fBasic->SetDirectory(0);
    fBasic->Sumw2();
  }

  fTotal = new TH2D(kTotalN, kTotalT, nX, xMin, xMax, nY, yMin, yMax);
  fTotal->SetDirectory(0);
  fTotal->SetXTitle(fBasic->GetXaxis()->GetTitle());
  fTotal->SetYTitle(fBasic->GetYaxis()->GetTitle());
  fTotal->Sumw2();
    
  fEmpty = static_cast<TH2D*>(fTotal->Clone(kEmptyN));
  fEmpty->SetTitle(kEmptyT);
  fEmpty->SetDirectory(0);
  // fEmpty->Sumw2();
}

//____________________________________________________________________
void
AliPoissonCalculator::Fill(UShort_t x, UShort_t y, Bool_t hit, Double_t weight)
{
  // 
  // Fill in an observation 
  // 
  // Parameters:
  //    x       X value 
  //    Y       Y value
  //    hit     True if hit 
  //    weight  Weight if this 
  //
  fTotal->Fill(x, y);
  if (hit) fBasic->Fill(x, y, weight);
  else     fEmpty->Fill(x, y);
}

//____________________________________________________________________
Double_t 
AliPoissonCalculator::CalculateMean(Double_t empty, Double_t total) const
{
  if (total <= 0) return 0;
  if (empty < .001) empty = .001;
  return -TMath::Log(empty/total);
}
//____________________________________________________________________
Double_t 
AliPoissonCalculator::CalculateCorrection(Double_t empty, Double_t total) const
{
  if (total <= 0) return 0;
  if (TMath::Abs(empty-total) < .001) empty = total - .001;
  return 1 / (1 - empty / total);
}

//____________________________________________________________________
Int_t
AliPoissonCalculator::GetReducedXBin(Int_t ix) const
{
  if (!fBasic) return 0;
  Double_t mx = fBasic->GetXaxis()->GetBinCenter(ix);
  return GetReducedXBin(mx);
}
//____________________________________________________________________
Int_t
AliPoissonCalculator::GetReducedXBin(Double_t x) const
{
  if (!fEmpty) return 0;
  return fEmpty->GetXaxis()->FindBin(x);
}
//____________________________________________________________________
Int_t
AliPoissonCalculator::GetReducedYBin(Int_t iy) const
{
  if (!fBasic) return 0;
  Double_t my = fBasic->GetYaxis()->GetBinCenter(iy);
  return GetReducedYBin(my);
}
//____________________________________________________________________
Int_t
AliPoissonCalculator::GetReducedYBin(Double_t y) const
{
  if (!fEmpty) return 0;
  return fEmpty->GetYaxis()->FindBin(y);
}



//____________________________________________________________________
TH2D*
AliPoissonCalculator::Result(Bool_t correct)
{
  // 
  // Calculate result and store in @a output
  // 
  // Return:
  //    The result histogram (fBase overwritten)
  //
  
  // Double_t total = fXLumping * fYLumping;
  
  for (Int_t ix = 1; ix <= fBasic->GetNbinsX(); ix++) { 
    // Double_t x        = fBasic->GetXaxis()->GetBinCenter(ix);
    Int_t    jx       = GetReducedXBin(ix); // fEmpty->GetXaxis()->FindBin(x);
    for (Int_t iy = 1; iy <= fBasic->GetNbinsY(); iy++) { 
      // Double_t y        = fBasic->GetYaxis()->GetBinCenter(iy);
      Int_t    jy       = GetReducedYBin(iy); // fEmpty->GetYaxis()->FindBin(y);
      Double_t empty    = fEmpty->GetBinContent(jx, jy);
      Double_t total    = fTotal->GetBinContent(jx, jy);
      Double_t hits     = fBasic->GetBinContent(ix,iy);
      // Mean in region of interest 
      Double_t poissonM = CalculateMean(empty, total);
      Double_t poissonC = (correct ? CalculateCorrection(empty, total) : 1);
      
      Double_t poissonV = hits * poissonM * poissonC;
      Double_t poissonE = TMath::Sqrt(poissonV);
      if(poissonV > 0) poissonE = TMath::Sqrt(poissonV);
	  
      fBasic->SetBinContent(ix,iy,poissonV);
      fBasic->SetBinError(ix,iy,poissonE);
    }
  }
  for (Int_t ix = 1; ix <= fEmpty->GetNbinsX(); ix++) { 
    for (Int_t iy = 1; iy <= fEmpty->GetNbinsY(); iy++) { 
      Double_t empty    = fEmpty->GetBinContent(ix, iy);
      Double_t total    = fTotal->GetBinContent(ix, iy);
      Double_t mean     = CalculateMean(empty, total);
      Double_t corr     = CalculateCorrection(empty, total);
      fEmptyVsTotal->Fill(total, empty);
      fMean->Fill(mean);
      fOcc->Fill(100 * (1 - empty/total));
      //Old fOcc->Fill(100 * (1 - TMath::PoissonI(0,mean)));
      fCorr->Fill(mean, corr);
    }
  }
  return fBasic;
}
  
//____________________________________________________________________
void
AliPoissonCalculator::Print(const Option_t*) const
{
  // 
  // Print information 
  // 
  // Parameters:
  //    option Not used
  //
  char ind[gROOT->GetDirLevel()+3];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  std::cout << ind << ClassName() << ": " << GetName() << '\n'
	    << std::boolalpha 
	    << ind << " X lumping:              " << fXLumping << '\n'
	    << ind << " Y lumping:              " << fYLumping << '\n'
	    << std::noboolalpha << std::endl;
}
//____________________________________________________________________
void
AliPoissonCalculator::Browse(TBrowser* b)
{
  // 
  // Browse this object
  // 
  // Parameters:
  //    b Object to browse 
  //
  if (fTotal)        b->Add(fTotal);
  if (fEmpty)        b->Add(fEmpty);
  if (fBasic)        b->Add(fBasic);
  if (fEmptyVsTotal) b->Add(fEmptyVsTotal);
  if (fMean)         b->Add(fMean);
  if (fOcc)          b->Add(fOcc);
  if (fCorr)         b->Add(fCorr);
}
// 
// EOF
//
