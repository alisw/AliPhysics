// Test TDime class
//
// Mikael.Mieskolainen@cern.ch

//#include "TMath.h"

void testTDime(Int_t nev = 500) {

  gSystem->Load("libEG.so");
  gSystem->Load("libdime.so");
  gSystem->Load("libTDime.so");

  TDime* dime = new TDime();
  dime->SetEnergyCMS(7000.0);
  dime->Initialize();

  // (pi+pi-) histograms
  TH1* hM = new TH1D("hM", "DIME;M_{#pi^{+}#pi^{-}} #[]{GeV/#it{c}^{2}}", 100,  0.0, 5.0);

  TClonesArray* particles = new TClonesArray("TParticle", 6);
  TParticle* part = NULL;
  TLorentzVector v[2];
  TLorentzVector vSum;

  // Event loop
  for (Int_t i = 0; i < nev; i++) {

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
    if (nPrimary != 2) {
      Printf("Error: nPrimary=%d != 2", nPrimary);
      continue;
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
