// -*- C++ -*-
// $Id$


TGraph* ConfHist(TGraph *g, const char* title) {
  g->GetHistogram()->SetTitle(title);
  gPad->Update();
  TPaveStats *ps = (TPaveStats*)g->FindObject("stats");
  ps->SetFillColor(kYellow);
  ps->Draw();
  return g;
}

void MakeLumiRegion_4690_LSC() {
  gROOT->LoadMacro("AliLuminousRegionFit.cxx+");
  gROOT->LoadMacro("Util.C");

  const Int_t t0 = 1449188165;
  const ScanData sd[] = {
    { 4690,
      11,
      "root/4690/AnalysisResults_101116.root",
      "txt/4690/SepVStime_steps.out",
      "Scan1",
      "Scan1X", 0, t0, t0+1020+90, 0,
    },
  };

  const Int_t n = sizeof(sd)/sizeof(ScanData);

  // fill 4269
  const Int_t bcs[1] = {
    -1,   //all
  };

#if 1
  for (Int_t i=0; i<n; ++i) {
    CheckCopyFile(sd[i].vtxFileName);
    for (Int_t j=0; j<1; ++j) {
      Int_t bcid = bcs[j];
      AliLuminousRegionFit f(sd[i].fillNumber,
			     sd[i].minNumberOfTracks,
			     sd[i].vtxFileName,
			     sd[i].sepFileName);
      f.DoFit(sd[i].scanName1, sd[i].t1, sd[i].t2, sd[i].scanType1, sd[i].offset1, bcid);
    }
  }
#else

  gStyle->SetOptFit(111);
  TFile::Open("root/4690/lumiRegion_Scan1X.root");

  TF1 *f = new TF1("f", "[0]+[1]*x+[2]*TMath::Floor((x+0.5)/60*5)");
  f->SetLineColor(kBlue);
  f->SetNpx(1000);

  TF1 *g = new TF1("f", "[0]+[1]*x+[2]*TMath::Floor((x+0.5)/60)");
  g->SetLineColor(kBlue);
  g->SetNpx(1000);

  TCanvas *c1 = new TCanvas("c1", "", 900, 500);
  c1->Divide(5,2);

  c1->cd(1);
  geX->Draw("APE");
  geX->Fit(f, "", "", 0, 59);
  ConfHist(geX, "<X>;step index;<X> (cm)");

  c1->cd(2);
  geY->Draw("APE");
  geY->Fit(f, "", "", 60, 119);
  ConfHist(geY, "<Y>;step index;<Y> (cm)");

  c1->cd(3);
  geZ->Draw("APE");
  geZ->Fit(g, "", "", 0, 119);
  ConfHist(geZ, "<Z>;step index;<Z> (cm)");

  c1->cd(6);
  geSX->Draw("APE");
  geSX->Fit(g, "", "", 0, 119);
  ConfHist(geSX, "#sigma_{X};step index;#sigma_{X} (cm)");

  c1->cd(7);
  geSY->Draw("APE");
  geSY->Fit(g, "", "", 0, 119);
  ConfHist(geSY, "#sigma_{Y};step index;#sigma_{Y} (cm)");

  c1->cd(8);
  geSZ->Draw("APE");
  geSZ->Fit(g, "", "", 0, 119);
  ConfHist(geSZ, "#sigma_{Z};step index;#sigma_{Z} (cm)");

  c1->cd(9);
  geCXY->Draw("APE");
  geCXY->Fit(g, "", "", 0, 119);
  ConfHist(geCXY, "#rho_{XY};step index;#rho_{XY}");

  c1->cd(10);
  gek->Draw("APE");
  gek->Fit(g, "", "", 0, 119);
  ConfHist(gek, "k;step index;k");

  c1->cd(4);
  gesX->Draw("APE");
  gesX->Fit(g, "", "", 0, 119);
  ConfHist(gesX, "#mu'_{X};step index;#mu'_{X}");

  c1->cd(5);
  gesY->Draw("APE");
  gesY->Fit(g, "", "", 0, 119);
  ConfHist(gesY, "#mu'_{Y};step index;#mu'_{Y}");

  c1->SaveAs("pdf/4690/LSC.pdf");
#endif

}
