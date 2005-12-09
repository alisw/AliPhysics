//
// $Id$
//
// Script that contains a class to draw hits, using the
// AliFMDInputHits class in the util library. 
//
// It draws the energy loss versus the p/(mq^2).  It can be overlayed
// with the Bethe-Bloc curve to show how the simulation behaves
// relative to the expected. 
//
// Use the script `Compile.C' to compile this class using ACLic. 
//
#include <TH2D.h>
#include <AliFMDHit.h>
#include <AliFMDInput.h>
#include <iostream>
#include <TStyle.h>
#include <TArrayF.h>

class DrawHits : public AliFMDInputHits
{
private:
  TH2D* fElossVsPMQ; // Histogram 
  Int_t fPdg;
public:
  TArrayF MakeLogScale(Int_t n, Double_t min, Double_t max) 
  {
    TArrayF bins(n+1);
    Float_t dp   = n / TMath::Log10(max / min);
    Float_t pmin = TMath::Log10(min);
    bins[0]      = min;
    for (Int_t i = 1; i < n+1; i++) {
      Float_t p = pmin + i / dp;
      bins[i]   = TMath::Power(10, p);
    }
    return bins;
  }
  DrawHits(Int_t m=1000, Double_t emin=1, Double_t emax=1000, 
	   Int_t n=900, Double_t tmin=1e-2, Double_t tmax=1e3, 
	   Int_t pdg=211) 
    : fPdg(pdg)
  { 
    TArrayF tkine(MakeLogScale(n, tmin, tmax));
    TArrayF eloss(MakeLogScale(m, emin, emax));
    fElossVsPMQ = new TH2D("bad", "#Delta E/#Delta x vs. p/(mq^{2})", 
			   tkine.fN-1, tkine.fArray, eloss.fN-1, eloss.fArray);
    fElossVsPMQ->SetXTitle(Form("p%s", (fPdg == 0 ? "/(mq^{2})" : " [GeV]")));
    fElossVsPMQ->SetYTitle("#Delta E/#Delta x [MeV/cm]");
  }
  Bool_t ProcessHit(AliFMDHit* hit, TParticle*) 
  {
    if (!hit) {
      std::cout << "No hit" << std::endl;
      return kFALSE;
    }

    if (hit->IsStop()) return kTRUE;
    Float_t x = hit->P();
    Float_t y = hit->Edep()/hit->Length();
    if (fPdg != 0) {
      if (TMath::Abs(hit->Pdg()) != fPdg) return kTRUE;
    }
    else {
      Float_t q = hit->Q() / 3.;
      if (hit->M() == 0 || q == 0) return kTRUE;
      x /= hit->M();
      y /= q * q;
    }
    fElossVsPMQ->Fill(x, y);
    return kTRUE;
  }
  Bool_t Finish()
  {
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderSize(0);
    fElossVsPMQ->SetStats(kFALSE);
    fElossVsPMQ->Draw("COLZ");
    return kTRUE;
  }
};

//____________________________________________________________________
//
// EOF
//
