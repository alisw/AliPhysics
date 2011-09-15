#include "AliPoissonCalculator.h"
#include <TH2D.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <TMath.h>
#include <TList.h>
#include <iostream>

// 
// A class to calculate the multiplicity in @f$(\eta,\varphi)@f$ bins
// using Poisson statistics. 
//
// The input is assumed to be binned in @f$(\eta,\varphi)@f$ as
// described by the 2D histogram passwd to the Reset member function.  
//
// The data is grouped in to regions as defined by the parameters
// fEtaLumping and fPhiLumping.  The total number of cells and number
// of empty cells is then calculate in each region.  The mean
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

//____________________________________________________________________
AliPoissonCalculator::AliPoissonCalculator()
  : TNamed(),
    fEtaLumping(32), 
    fPhiLumping(4), 
    fTotal(0), 
    fEmpty(0), 
    fBasic(0),
    fEmptyVsTotal(0),
    fMean(0), 
    fOcc(0),
    fCorr(0),
    fTotalList(),
    fEmptyList(),
    fRunningAverage(false)
{
  //
  // CTOR
  // 
}

//____________________________________________________________________
AliPoissonCalculator::AliPoissonCalculator(const char*/*, UShort_t d, Char_t r*/)
  : TNamed("poissonCalculator", "Calculate N_ch using Poisson stat"),
    fEtaLumping(32), 
    fPhiLumping(4), 
    fTotal(0), 
    fEmpty(0), 
    fBasic(0),
    fEmptyVsTotal(0),
    fMean(0), 
    fOcc(0),
    fCorr(0),
    fTotalList(),
    fEmptyList(),
    fRunningAverage(false)
{
  //
  // CTOR
  //
  fEmptyList.SetOwner();
  fTotalList.SetOwner();
  

}
//____________________________________________________________________
AliPoissonCalculator::AliPoissonCalculator(const AliPoissonCalculator& o)
  : TNamed(o),
    fEtaLumping(o.fEtaLumping),
    fPhiLumping(o.fPhiLumping),
    fTotal(0), 
    fEmpty(0),
    fBasic(0), 
    fEmptyVsTotal(0),
    fMean(0), 
    fOcc(0),
    fCorr(0),
    fTotalList(),
    fEmptyList(),
    fRunningAverage(o.fRunningAverage)
{
  Init();
  Reset(o.fBasic);
}

//____________________________________________________________________
AliPoissonCalculator::~AliPoissonCalculator()
{
  CleanUp();
}

//____________________________________________________________________
void
AliPoissonCalculator::CleanUp()
{
  if (fTotal)        delete fTotal;        fTotal        = 0;
  if (fEmpty)        delete fEmpty;        fEmpty        = 0;
  if (fBasic)        delete fBasic;        fBasic        = 0;
  if (fEmptyVsTotal) delete fEmptyVsTotal; fEmptyVsTotal = 0;
  if (fMean)         delete fMean;         fMean         = 0;
  if (fOcc)          delete fOcc;          fOcc          = 0;
  if (fCorr)         delete fCorr;         fCorr         = 0;
}
//____________________________________________________________________
AliPoissonCalculator&
AliPoissonCalculator::operator=(const AliPoissonCalculator& o)
{
  TNamed::operator=(o);
  fEtaLumping = o.fEtaLumping;
  fPhiLumping = o.fPhiLumping;
  fRunningAverage = o.fRunningAverage;
  CleanUp();
  Init();
  Reset(o.fBasic);
  return *this;
}

//____________________________________________________________________
void
AliPoissonCalculator::Init(UShort_t d, Char_t r, Int_t etaLumping, Int_t phiLumping)
{
  // 
  // Initialize 
  // 
  if (etaLumping > 0) SetEtaLumping(etaLumping);
  if (phiLumping > 0) SetPhiLumping(phiLumping);
  if(d > 0) {

    Int_t    nEtaF   = (r == 'I' ? 512 : 256);
    Int_t    nEta    = nEtaF / fEtaLumping;
    Double_t etaMin  = -0.5;
    Double_t etaMax  = nEtaF-0.5 ;
    Int_t    nPhiF   = (r == 'I' ? 20 : 40);
    Int_t    nPhi    = nPhiF / fPhiLumping;
    Double_t phiMin  = -0.5;
    Double_t phiMax  = nPhiF - 0.5;
    
    fBasic = new TH2D("basic", "Basic number of hits",
		      nEtaF, etaMin, etaMax, nPhiF, phiMin, phiMax);
    fBasic->SetDirectory(0);
    fBasic->SetXTitle("#eta");
    fBasic->SetYTitle("#varphi [radians]");
    fBasic->Sumw2();
    
    for(Int_t v = 1 ; v < 11; v++)  { //CHC bins
      for(Int_t centbin = 0 ; centbin < 13; centbin++) {
	
	TH2D* hTotal = new TH2D(Form("totalFMD%d%c_vertex%d_cent%d",d,r,v,centbin),"Total number of bins/region",
			      nEta, etaMin, etaMax, nPhi, phiMin, phiMax);
	TH2D* hEmpty = new TH2D(Form("emptyFMD%d%c_vertex%d_cent%d",d,r,v,centbin), "Empty number of bins/region",
				nEta, etaMin, etaMax, nPhi, phiMin, phiMax);
	hEmpty->Sumw2();
	hTotal->Sumw2();
	fEmptyList.Add(hEmpty);
	fTotalList.Add(hTotal);
	
      }
    }
  }
  //Create diagnostics if void
  if (fEmptyVsTotal) return;
  
  Int_t n = fEtaLumping * fPhiLumping + 1;
  fEmptyVsTotal = new TH2D("emptyVsTotal", 
			   "# of empty # bins vs total # bins", 
			   n, -.5, n-.5, n, -.5, n-.5);
  fEmptyVsTotal->SetXTitle("Total # bins");
  fEmptyVsTotal->SetYTitle("# empty bins");
  fEmptyVsTotal->SetZTitle("Correlation");
  fEmptyVsTotal->SetOption("colz");
  fEmptyVsTotal->SetDirectory(0);

  n = (fEtaLumping + fPhiLumping);
  fMean = new TH1D("poissonMean", "Mean N_{ch} as calculated by Poisson",
		   10*n+1, -.1, n+.1);
  fMean->SetXTitle("#bar{N_{ch}}=#log(empty/total)");
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
void AliPoissonCalculator::SetObject(UShort_t d, Char_t r, UShort_t v, Double_t cent) {
  
  Int_t centbin = 0;
  if(cent > 0) {
    if(cent > 0 && cent <5) centbin = 1;
    if(cent > 5 && cent <10) centbin = 2;
    else if (cent>10) centbin = (Int_t)(cent/10.) + 2;
  }
  
  fTotal = static_cast<TH2D*>(fTotalList.FindObject(Form("totalFMD%d%c_vertex%d_cent%d",d,r,v,centbin)));
  fEmpty = static_cast<TH2D*>(fEmptyList.FindObject(Form("emptyFMD%d%c_vertex%d_cent%d",d,r,v,centbin)));
  
  return;
  
}
//____________________________________________________________________
void
AliPoissonCalculator::Output(TList* d)
{
  if (!d) return;
  if (!fEmptyVsTotal) Init();
  d->Add(fEmptyVsTotal);
  d->Add(fMean);
  d->Add(fOcc);
  d->Add(fCorr);
}

//____________________________________________________________________
void
AliPoissonCalculator::Reset(const TH2D* base)
{
  // 
  // Reset histogram 
  // 
  if (!base) return;
  if (fBasic /* && fTotal && fEmpty*/) {
    fBasic->Reset();
    if(!fRunningAverage) {
      fTotal->Reset();
      fEmpty->Reset();
    }
    return;
  }
  /*  
  Int_t    nEtaF   = base->GetNbinsX();
  Int_t    nEta    = nEtaF / fEtaLumping;
  Double_t etaMin  = base->GetXaxis()->GetXmin();
  Double_t etaMax  = base->GetXaxis()->GetXmax();
  Int_t    nPhiF   = base->GetNbinsY();
  Int_t    nPhi    = nPhiF / fPhiLumping;
  Double_t phiMin  = base->GetYaxis()->GetXmin();
  Double_t phiMax  = base->GetYaxis()->GetXmax();
  

  
  
  
  //fTotal = new TH2D("total", "Total number of bins/region",
  //		    nEta, etaMin, etaMax, nPhi, phiMin, phiMax);
  //fEmpty = new TH2D("empty", "Empty number of bins/region",
  //nEta, etaMin, etaMax, nPhi, phiMin, phiMax);
  fBasic = new TH2D("basic", "Basic number of hits",
		    nEtaF, etaMin, etaMax, nPhiF, phiMin, phiMax);
    
  //fTotal->SetDirectory(0);
  //fEmpty->SetDirectory(0);
  fBasic->SetDirectory(0);
  //fTotal->SetXTitle("#eta");
  //fEmpty->SetXTitle("#eta");
  fBasic->SetXTitle("#eta");
  //fTotal->SetYTitle("#varphi [radians]");
  //fEmpty->SetYTitle("#varphi [radians]");
  fBasic->SetYTitle("#varphi [radians]");
  //fTotal->Sumw2();
  //fEmpty->Sumw2();
  fBasic->Sumw2();
  */

}

//____________________________________________________________________
void
AliPoissonCalculator::Fill(UShort_t strip, UShort_t sec, Bool_t hit, 
			   Double_t weight)
{
  // 
  // Fill in an observation 
  // 
  // Parameters:
  //    eta     Eta value 
  //    phi     Phi value
  //    hit     True if hit 
  //    weight  Weight if this 
  //
  
  fTotal->Fill(strip, sec);
  if (hit) fBasic->Fill(strip, sec, weight);
  else     fEmpty->Fill(strip, sec);
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
TH2D*
AliPoissonCalculator::Result()
{
  // 
  // Calculate result and store in @a output
  // 
  // Return:
  //    The result histogram (fBase overwritten)
  //
  
  // Double_t total = fEtaLumping * fPhiLumping;
  
  for (Int_t ieta = 1; ieta <= fBasic->GetNbinsX(); ieta++) { 
    for (Int_t iphi = 1; iphi <= fBasic->GetNbinsY(); iphi++) { 
      Double_t eta      = fBasic->GetXaxis()->GetBinCenter(ieta);
      Double_t phi      = fBasic->GetYaxis()->GetBinCenter(iphi);
      Int_t    jeta     = fEmpty->GetXaxis()->FindBin(eta);
      Int_t    jphi     = fEmpty->GetYaxis()->FindBin(phi);
      Double_t empty    = fEmpty->GetBinContent(jeta, jphi);
      Double_t total    = fTotal->GetBinContent(jeta, jphi);
      Double_t hits     = fBasic->GetBinContent(ieta,iphi);
      // Mean in region of interest 
      Double_t poissonM = CalculateMean(empty, total);
      Double_t poissonC = CalculateCorrection(empty, total);
      
      Double_t poissonV = hits * poissonM * poissonC;
      Double_t poissonE = TMath::Sqrt(poissonV);
      if(poissonV > 0) poissonE = TMath::Sqrt(poissonV);
	  
      fBasic->SetBinContent(ieta,iphi,poissonV);
      fBasic->SetBinError(ieta,iphi,poissonE);
    }
  }
  for (Int_t ieta = 1; ieta <= fEmpty->GetNbinsX(); ieta++) { 
    for (Int_t iphi = 1; iphi <= fEmpty->GetNbinsY(); iphi++) { 
      Double_t empty    = fEmpty->GetBinContent(ieta, iphi);
      Double_t total    = fTotal->GetBinContent(ieta, iphi);
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
	    << ind << " Eta lumping:            " << fEtaLumping << '\n'
	    << ind << " Phi lumping:            " << fPhiLumping << '\n'
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
