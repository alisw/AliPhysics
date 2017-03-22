// -*- C++ -*-
// $Id$

// creates histogram of generated rho0 kinematical data
// writes a text file "rho0.txt" which should be identical to ../startlight/TSTarlight/rho0.txt

void Kin2Txt() {
  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  rl->LoadKinematics();
  rl->LoadHeader();

  TH1* hM  = new TH1D("hM",  "TreeK;M#(){#pi^{+}#pi^{-}} #(){GeV/#it{c}^{2}}", 100, 0., 2.);
  TH1* hPt = new TH1D("hPt", "TreeK;P_{T}#(){#pi^{+}#pi^{-}} #(){GeV/#it{c}}", 100, 0., 1.);
  TH1* hY  = new TH1D("hY",  "TreeK;Y#(){#pi^{+}#pi^{-}}",                     160,-8., 8.);
  std::ofstream ofs("rho0.txt");
  
  AliStack *stack = NULL;
  TParticle *part = NULL;
  TLorentzVector v[2], vSum;
  for (Int_t i=0, n(rl->GetNumberOfEvents()); i<n; ++i) {
    rl->GetEvent(i);
    stack = rl->Stack();

    Int_t nPrimary(0);
    for (Int_t j(0), m(stack->GetNtrack()); j<m; ++j) {
      part = stack->Particle(j);
      if (part->IsPrimary())
	part->Momentum(v[nPrimary++]);
    }
    if (nPrimary != 2) {
      Printf("Error: nPrimary=%d != 2", nPrimary);
      continue;
    }
    vSum = v[0] + v[1];
    hM->Fill(vSum.M());
    hPt->Fill(vSum.Perp());
    hY->Fill(vSum.Rapidity());
    ofs << std::fixed << std::setprecision(4)
	<< vSum.M() << " " << vSum.Perp() << " " << vSum.Rapidity() << " "
	<< v[0].Eta() << " " << v[0].Px() << " " << v[0].Py() << " " << v[0].Pz() << " "
	<< v[1].Eta() << " " << v[1].Px() << " " << v[1].Py() << " " << v[1].Pz()
	<< std::endl;
  }
  hM->Draw();
  c1->SaveAs("TreeK.pdf(");
  hPt->Draw();
  c1->SaveAs("TreeK.pdf");
  hY->Draw();
  c1->SaveAs("TreeK.pdf)");
}
