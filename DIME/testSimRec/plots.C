// Histograms and text file output check of DIME output from Sim run
//
// Mikael.Mieskolainen@cern.ch

#include "TMath.h"

void plots() {

  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  rl->LoadKinematics();
  rl->LoadHeader();

  // (pi+pi-) histograms
  TH1* hM    = new TH1D("hM",    "DIME;M_{#pi^{+}#pi^{-}} #(){GeV/#it{c}^{2}}", 100,  0.0, 20.0);
  TH1* hPt   = new TH1D("hPt",   "DIME;p_{T}#(){#pi^{+}#pi^{-}} #(){GeV/#it{c}}", 100,  0.0, 3.0);

  // pi+- histograms
  TH1* hPt1   = new TH1D("hPt1",  "DIME;p_{T}#(){#pi^{#pm}} #(){Gev/#it{c}}",  100, 0.0, 5.0);

  std::ofstream ofs("pipm.txt");

  AliStack* stack = NULL;
  TParticle* part = NULL;
  TLorentzVector v[2];
  TLorentzVector vSum;

  // Loop over events
  for (Int_t i = 0, n(rl->GetNumberOfEvents()); i < n; ++i) {

    rl->GetEvent(i);
    stack = rl->Stack();
    Int_t nPrimary(0);

    // Loop over tracks
    for (Int_t j = 0, m(stack->GetNtrack()); j < m; ++j) {
      part = stack->Particle(j);     // Get particle
      if (part->IsPrimary()) {
        part->Momentum(v[nPrimary]); // Set content of v
        nPrimary++;
      }
    }
    if (nPrimary != 2) {
      Printf("Error: nPrimary=%d != 2", nPrimary);
      continue;
    }

    // 4-vector sum
    vSum = v[0] + v[1];

    // Fill pi+pi- histograms
    hM->Fill(vSum.M());
    hPt->Fill(vSum.Perp());

    // Fill pi(+-) histograms
    for (Int_t k = 0; k < 2; ++k) {
      hPt1->Fill(v[k].Perp());
    }

    // Text file output
    ofs << std::fixed << std::setprecision(4)
	<< vSum.M() << " " << vSum.Perp() << " " << vSum.Rapidity() << " "
	<< v[0].Eta() << " " << v[0].Px() << " " << v[0].Py() << " " << v[0].Pz() << " "
	<< v[1].Eta() << " " << v[1].Px() << " " << v[1].Py() << " " << v[1].Pz()
	<< std::endl;
  }

  // Save plots as pdf
  hM->Draw();    c1->SaveAs("plotM.pdf");
  hPt->Draw();   c1->SaveAs("plotPt.pdf");
  hPt1->Draw();  c1->SaveAs("plotPt1.pdf");

}
