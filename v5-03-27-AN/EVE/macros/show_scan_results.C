// Display some histograms from scanning.
//
// BIT(1) stores the original selection.
// BIT(0) stores the user selection (set to same value as b1 at init).
//
// This allows to check all possible combinations.


void show_scan_results()
{
  TFile *f = TFile::Open("scan_results.root");

  TTree* t = (TTree*) gDirectory->Get("SR");

  if (t == 0)
    Error("show_scan_results", "Tree 'SR' with scan results not found.");

  TCanvas *c = 0;


  //----------------------------------------------------------------------------
  // Tracks
  //----------------------------------------------------------------------------

  c = new TCanvas("Tracks", "Track Scanning Results", 800, 600);
  c->Divide(2, 3);

  c->cd(1);
  t->Draw("Sum$(T.fLabel & 1)");

  c->cd(2);
  t->Draw("T.GetSign()", "T.fLabel & 1");

  c->cd(3);
  t->Draw("T.Pt()", "T.fLabel & 1");

  c->cd(4);
  t->Draw("T.Eta()", "T.fLabel & 1");

  c->cd(5);
  t->Draw("T.Phi()", "T.fLabel & 1");

  c->Modified();
  c->Update();


  //----------------------------------------------------------------------------
  // Trackelts
  //----------------------------------------------------------------------------

  c = new TCanvas("Tracklets", "Tracklet Scanning Results", 800, 600);
  c->Divide(2, 3);

  c->cd(1);
  t->Draw("Sum$(M.fLabels & 1)");

  c->cd(2);
  t->Draw("M.fNsingle");

  c->cd(3);
  t->Draw("M.fFiredChips[1]:Sum$(M.fLabels & 1)");

  c->cd(4);
  t->Draw("M.fDeltTh", "M.fLabels & 1");

  c->cd(5);
  t->Draw("M.fDeltPhi", "M.fLabels & 1");

  c->cd(6);
  t->Draw("M.fPhi", "M.fLabels & 1");

  c->Modified();
  c->Update();


  //----------------------------------------------------------------------------
  // Vertices
  //----------------------------------------------------------------------------

  c = new TCanvas("Vertices", "Vertex Scanning Results", 800, 600);
  c->Divide(3, 3);

  c->cd(1);
  t->Draw("VT.GetX()", "VT.GetNContributors()>0");

  c->cd(2);
  t->Draw("VT.GetY()", "VT.GetNContributors()>0");

  c->cd(3);
  t->Draw("VT.GetZ()", "VT.GetNContributors()>0");

  c->cd(4);
  t->Draw("VSPD.GetX()", "VSPD.GetNContributors()>0 && VSPD.fTitle.Contains(\"3D\")");

  c->cd(5);
  t->Draw("VSPD.GetY()", "VSPD.GetNContributors()>0 && VSPD.fTitle.Contains(\"3D\")");

  c->cd(6);
  t->Draw("VSPD.GetZ()", "VSPD.GetNContributors()>0 && VSPD.fTitle.Contains(\"3D\")");

  c->cd(7);
  t->Draw("VTPC.GetX()", "VTPC.GetNContributors()>0");

  c->cd(8);
  t->Draw("VTPC.GetY()", "VTPC.GetNContributors()>0");

  c->cd(9);
  t->Draw("VTPC.GetZ()", "VTPC.GetNContributors()>0");

  //----------------------------------------------------------------------------
  // End
  //----------------------------------------------------------------------------

  f->Close();
  delete f;
}
