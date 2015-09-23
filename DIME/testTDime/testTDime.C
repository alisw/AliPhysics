// Test TDime class
//
// Mikael.Mieskolainen@cern.ch, 22.9.2015

//#include "TMath.h"

void testTDime(Int_t nev = 100) {

  gSystem->Load("libEVGEN");
  gSystem->Load("libTDime");
  gSystem->Load("libdime");

  TDime* dime = new TDime();
  dime->SetEnergyCMS(7000.0);
  dime->SetYRange(-2.0, 2.0);   // Set rapidity range of mesons
  dime->SetMinPt(0.1);          // Minimum pT of mesons
  dime->Initialize();

  // (pi+pi-) histograms
  TH1* hM = new TH1D("hM", "DIME #pi^{+}#pi^{-};M_{#pi^{+}#pi^{-}} #[]{GeV/#it{c}^{2}}", 100,  0.0, 5.0);

  TClonesArray* particles = new TClonesArray("TParticle", 6);
  TParticle* part = NULL;
  TLorentzVector v[2];
  TLorentzVector vSum;

  // Event loop
  for (Int_t i = 0; i < nev; ++i) {

    dime->GenerateEvent();
    Int_t np = dime->ImportParticles(particles, "All");
    printf("\n DIME Event %d: Imported %3d particles \n", i, np);

    Int_t nPrimary = 0;

    // Loop over pion (j = 4,5) tracks
    for (Int_t j = 4; j < 6; ++j) {
      part = (TParticle*) particles->At(j); // Choose the particle
      part->Print();
      part->Momentum(v[nPrimary]);          // Copy content to v
      nPrimary++;
    }
    //particles.Clear();

    // 4-vector sum
    vSum = v[0] + v[1];

    // Fill pi+pi- histograms
    hM->Fill(vSum.M());
  }

  // Save plots as pdf
  hM->Draw();    c1->SaveAs("massTDime.pdf");

}
