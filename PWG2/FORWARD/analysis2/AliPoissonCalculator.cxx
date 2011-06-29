#include "AliPoissonCalculator.h"
#include <TH2D.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <iostream>

//____________________________________________________________________
AliPoissonCalculator::AliPoissonCalculator()
  : TNamed(),
    fEtaLumping(5), 
    fPhiLumping(5), 
    fTotal(0), 
    fEmpty(0), 
    fBasic(0)
{}

//____________________________________________________________________
AliPoissonCalculator::AliPoissonCalculator(const char*)
  : TNamed("poissonCalculator", "Calculate N_ch using Poisson stat"),
    fEtaLumping(5), 
    fPhiLumping(5), 
    fTotal(0), 
    fEmpty(0), 
    fBasic(0)
{}

//____________________________________________________________________
void
AliPoissonCalculator::Reset(const TH2D* base)
{
  if (fBasic && fTotal && fEmpty) {
    fBasic->Reset();
    fTotal->Reset();
    fEmpty->Reset();
    return;
  }
  
  Int_t    nEtaF   = base->GetNbinsX();
  Int_t    nEta    = nEtaF / fEtaLumping;
  Double_t etaMin  = base->GetXaxis()->GetXmin();
  Double_t etaMax  = base->GetXaxis()->GetXmax();
  Int_t    nPhiF   = base->GetNbinsY();
  Int_t    nPhi    = nPhiF / phiLumping;
  Double_t phiMin  = base->GetYaxis()->GetXmin();
  Double_t phiMax  = base->GetYaxis()->GetXmax();

  fTotal = new TH2D("total", "Total number of bins/region",
		    nEta, etaMin, etaMax, nPhi, phiMin, phiMax);
  fEmpty = new TH2D("empty", "Empty number of bins/region",
		    nEta, etaMin, etaMax, nPhi, phiMin, phiMax);
  fBasic = new TH2D("basic", "Basic number of hits",
		    nEtaF, etaMin, etaMax, nPhiF, phiMin, phiMax);
  
  fTotal->SetDirectory(0);
  fEmpty->SetDirectory(0);
  fBasic->SetDirectory(0);
  fTotal->SetXTitle("#eta");
  fEmpty->SetXTitle("#eta");
  fBasic->SetXTitle("#eta");
  fTotal->SetYTitle("#varphi [radians]");
  fEmpty->SetYTitle("#varphi [radians]");
  fBasic->SetYTitle("#varphi [radians]");
  fTotal->Sumw2();
  fEmpty->Sumw2();
  fBasic->Sumw2();
}

//____________________________________________________________________
void
AliPoissonCalculator::Fill(Double_t eta, Double_t phi, Bool_t hit, 
			   Double_t weight)
{
  fTotal->Fill(eta, phi);
  if (hit) fBasic->Fill(eta, phi, weight);
  else     fEmpty->Fill(eta, phi);
}

//____________________________________________________________________
void 
AliPoissonCalculator::Result(TH2D* output)
{
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
      Double_t poissonM = (total <= 0 || empty <= 0 ? 0 : 
			   -TMath::Log(empty / total));

      // Note, that given filled=total-empty, and 
      //
      //     m = -log(empty/total)
      //       = -log(1 - filled/total)
      // 
      //     v = m / (1 - exp(-m))
      //       = -total/filled * (log(total-filled)-log(total))
      //       = -total / (total-empty) * log(empty/total)
      //       = total (log(total)-log(empty)) / (total-empty)
      //  
      Double_t poissonV = hits;
      if(poissonM > 0)
	// Correct for counting statistics and weight by counts 
	poissonV *= poissonM / (1 - TMath::Exp(-1*poissonM));
      Double_t poissonE = TMath::Sqrt(hits);
      if(poissonV > 0) poissonE = TMath::Sqrt(poissonV);
	  
      output->SetBinContent(ieta,iphi,poissonV);
      output->SetBinError(ieta,iphi,poissonE);
    }
  }
}
  
//____________________________________________________________________
void
AliPoissonCalculator::Print(const Option_t*) const
{
  char ind[gROOT->GetDirLevel()+3];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  std::cout << ind << ClassName() << ": " << GetName() << '\n'
	    << std::boolalpha 
	    << ind << " Eta lumping:            " << fEtaLumping << '\n'
	    << ind << " Phi lumping:            " << fPhiLumping << '\n'
	    << std::endl;
}
//____________________________________________________________________
void
AliPoissonCalculator::Browse(TBrowser* b)
{
  if (fBasic) b->Add(fBasic);
  if (fTotal) b->Add(fTotal);
  if (fEmpty) b->Add(fEmpty);

}
// 
// EOF
//
