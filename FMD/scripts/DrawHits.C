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

class DrawHits : public AliFMDInputHits
{
private:
  TH2D* fElossVsPMQ; // Histogram 
public:
  DrawHits() 
  { 
    fElossVsPMQ = new TH2D("bad", "#Delta E vs. p/(mq^{2})>1GeV", 
		     1000, 1, 100, 50, 0.00001, 10);
    fElossVsPMQ->SetXTitle("p/(mq^{2}) [GeV/GeV]");
    fElossVsPMQ->SetYTitle("#Delta E [MeV]");
  }
  Bool_t ProcessHit(AliFMDHit* hit, TParticle*) 
  {
    if (!hit) {
      std::cout << "No hit" << std::endl;
      return kFALSE;
    }
    Float_t pmq = 0;
    if (hit->M() != 0 && hit->Q() != 0) 
      pmq = hit->P() / hit->M() / TMath::Power(hit->Q()/3, 2);
    fElossVsPMQ->Fill(pmq, hit->Edep());
    return kTRUE;
  }
  Bool_t Finish()
  {
    gStyle->SetPalette(1);
    fElossVsPMQ->SetStats(kFALSE);
    fElossVsPMQ->Draw("COLZ");
    return kTRUE;
  }
};

//____________________________________________________________________
//
// EOF
//
