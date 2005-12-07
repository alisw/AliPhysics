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
public:
  DrawHits() 
  { 
    Double_t emax = 1e4;
    Double_t emin = 1e-5;
    Int_t    n    = 901;
    TArrayF tkine(n);
    Float_t dp   = 1/TMath::Log10(emax/emin)/10;
    Float_t pmin = TMath::Log10(emin);
    tkine[0]     = emin;
    for (Int_t i=1; i < tkine.fN; i++) {
      Float_t el = pmin + i * dp;
      tkine[i]   = TMath::Power(10, el);
    }
    Double_t dmin = .1;
    Double_t dmax = 10;
    Int_t    nd = 1000;
    TArrayF eloss(nd+1);
    eloss[0] = dmin;
    for (Int_t i = 1; i < eloss.fN; i++){
      eloss[i] = dmin + i * (dmax-dmin)/(eloss.fN-1);
    }
    fElossVsPMQ = new TH2D("bad", "#Delta E vs. p/(mq^{2})", 
			   tkine.fN-1, tkine.fArray, eloss.fN-1, eloss.fArray);
    fElossVsPMQ->SetXTitle("p/(mq^{2})");
    fElossVsPMQ->SetYTitle("#Delta E [MeV]");
  }
  Bool_t ProcessHit(AliFMDHit* hit, TParticle*) 
  {
    if (!hit) {
      std::cout << "No hit" << std::endl;
      return kFALSE;
    }
    if (TMath::Abs(hit->Pdg()) != 211) return kTRUE;
    // Float_t pmq = 0;
    // Float_t q = hit->Q() / 3.;
    // if (hit->M() != 0 && hit->Q() != 0) 
    //  pmq = hit->P() / (hit->M()*q*q);
    Float_t pmq = hit->P();
    fElossVsPMQ->Fill(pmq, hit->Edep());
    return kTRUE;
  }
  Bool_t Finish()
  {
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);
    fElossVsPMQ->SetStats(kFALSE);
    fElossVsPMQ->Draw("COLZ");
    return kTRUE;
  }
};

//____________________________________________________________________
//
// EOF
//
