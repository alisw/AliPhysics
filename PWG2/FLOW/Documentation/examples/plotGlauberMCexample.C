{
  // open file
  TFile* f = new TFile("GlauberMC_PbPb_ntuple.root");
  TTree* t = (TTree*)gDirectory->Get("nt_Pb_Pb"); // get the tree

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetStatBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetPalette(1);

  TCanvas* c = new TCanvas("GlauberMC","GlauberMC");
  c->Divide(2,3);
  c->cd(1);
  // plot number of collisions
  t->Draw("Ncoll:B","","prof");
  c->cd(2);
  // plot number of wounded nucleons
  t->Draw("Npart:B","","prof");
  c->cd(3);
  t->Draw("VarE:Npart","","prof");
  c->cd(4);
  t->Draw("VarEpart:Npart","","prof");
  c->cd(5);
  t->Draw("Mult:B");
  c->cd(6);
  t->Draw("MultGBW:B");

  c->cd(0);

  c->Update();
}
