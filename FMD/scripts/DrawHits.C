//____________________________________________________________________
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
#include <TParticle.h>

/** @class DrawHits
    @brief Draw hit energy loss
    @code 
    Root> .L Compile.C
    Root> Compile("DrawHits.C")
    Root> DrawHits c
    Root> c.Run();
    @endcode
    @ingroup FMD_script
 */
class DrawHits : public AliFMDInput
{
private:
  TH2D* fElossVsPMQ; // Histogram 
  Int_t fPdg;
public:
  //__________________________________________________________________
  TArrayF MakeLogScale(Int_t n, Double_t min, Double_t max) 
  {
    TArrayF bins(n+1);
    bins[0]      = min;
    if (n <= 20) {
      for (Int_t i = 1; i < n+1; i++) bins[i] = bins[i-1] + (max-min)/n;
      return bins;
    }
    Float_t dp   = n / TMath::Log10(max / min);
    Float_t pmin = TMath::Log10(min);
    for (Int_t i = 1; i < n+1; i++) {
      Float_t p = pmin + i / dp;
      bins[i]   = TMath::Power(10, p);
    }
    return bins;
  }
  //__________________________________________________________________
  DrawHits(Int_t m=1000, Double_t emin=1, Double_t emax=1000, 
	   Int_t n=900, Double_t tmin=1e-2, Double_t tmax=1e3, 
	   Int_t pdg=211) 
    : fPdg(pdg)
  { 
    AddLoad(kKinematics);
    AddLoad(kHits);
    TArrayF tkine(MakeLogScale(n, tmin, tmax));
    TArrayF eloss(MakeLogScale(m, emin, emax));
    TString name(Form("elossVsP%s", (fPdg == 0 ? "MQ" : "")));
    TString title(Form("#Delta E/#Delta x vs. p%s", 
		       (fPdg == 0 ? "/(mq^{2})" : "")));
    fElossVsPMQ = new TH2D(name.Data(), title.Data(), 
			   tkine.fN-1, tkine.fArray, 
			   eloss.fN-1, eloss.fArray);
    fElossVsPMQ->SetXTitle(Form("p%s", (fPdg == 0 ? "/(mq^{2})" : " [GeV]")));
    fElossVsPMQ->SetYTitle("#Delta E/#Delta x [MeV/cm]");
  }
  //__________________________________________________________________
  Bool_t ProcessHit(AliFMDHit* hit, TParticle* p) 
  {
    if (!hit) {
      std::cout << "No hit" << std::endl;
      return kFALSE;
    }

    if (!p) {
      std::cout << "No track" << std::endl;
      return kFALSE;
    }
    if (!p->IsPrimary()) return kTRUE;
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
  //__________________________________________________________________
  Bool_t Finish()
  {
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderSize(0);
    fElossVsPMQ->SetStats(kFALSE);
    fElossVsPMQ->Draw("COLZ box");
    return kTRUE;
  }
  ClassDef(DrawHits,0);
};

//____________________________________________________________________
//
// EOF
//
