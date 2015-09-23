// Histogram check of DIME output of p+p -> p + X + p, X(rhorho) -> 4pi from the simulation after GEANT
// Reads data using galice.root and goes through the kinematics tree, and selects generator level pions
//
// Mikael.Mieskolainen@cern.ch, 22.9.2015

#include "TMath.h"

void plots() {

  gSystem->Load("libEVGEN"); // Needs to be!

  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  rl->LoadKinematics();
  rl->LoadHeader();

  // 4pi histograms
  TH1* hM    = new TH1D("hM",  "DIME #rho#rho;M_{4#pi} #(){GeV/#it{c}^{2}}", 100,  1.0, 3.0);
  TH1* hPt   = new TH1D("hPt", "DIME #rho#rho;p_{T}#(){4#pi} #(){GeV/#it{c}}", 100,  0.0, 3.0);

  // pi+- histograms
  TH1* hPt1   = new TH1D("hPt1", "DIME #rho#rho;p_{T}#(){#pi^{#pm}} #(){Gev/#it{c}}",  100, 0.0, 3.0);

  AliStack* stack = NULL;
  TParticle* part = NULL;
  TLorentzVector v[4];
  TLorentzVector vSum;

  // Loop over events
  for (Int_t i = 0; i < rl->GetNumberOfEvents(); ++i) {

    rl->GetEvent(i);
    stack = rl->Stack();
    Int_t nPrimary = 0;

    // Loop over all particles
    for (Int_t j = 0; j < stack->GetNtrack(); ++j) {
      part = stack->Particle(j);           // Get particle
      part->Print();                       // Print contents

      if (abs(part->GetPdgCode()) == 211   // Is pi+ or pi-
          & part->GetStatusCode() == 1     // Is stable final state
          & stack->IsPhysicalPrimary(j)) { // Is from generator level

        part->Momentum(v[nPrimary]);       // Set content of v
        ++nPrimary;
      }
    }
    if (nPrimary != 4) {
      printf("Error: nPrimary=%d != 4 \n", nPrimary);
      continue;
    }

    // 4-vector sum
    vSum = v[0] + v[1] + v[2] + v[3];

    // Fill 4pi histograms
    hM->Fill(vSum.M());
    hPt->Fill(vSum.Perp());

    // Fill pi+- histograms
    for (Int_t k = 0; k < 4; ++k) {
      hPt1->Fill(v[k].Perp());
    }
    printf("\n");
  }

  // Save plots as pdf
  hM->Draw();    c1->SaveAs("plotM.pdf");
  hPt->Draw();   c1->SaveAs("plotPt.pdf");
  hPt1->Draw();  c1->SaveAs("plotPt1.pdf");

}
