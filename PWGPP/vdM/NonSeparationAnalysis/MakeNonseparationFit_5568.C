// -*- C++ -*-


const char* momentFileNames[] = {
  "root/5568/lumiRegion_Scan1X.root",
  "root/5568/lumiRegion_Scan1Y.root",

  "root/5568/lumiRegion_Scan2X.root",
  "root/5568/lumiRegion_Scan2Y.root",
};

void MakeFit(TString options) {
  TVirtualFitter::SetMaxIterations(100*1000);

  AliNonseparationModelFit f;

  f.SetVar( 0, 8.0e-3, 0.0001, 0.001,  0.0142);
  f.SetVar( 1, 7.0e-3, 0.0001, 0.001,  0.0142);
  f.SetVar( 2, 7.0   , 0.01, 6.5, 8.5);
  f.SetVar( 3, 0.0   , 0.01, -1, 1);

  f.SetVar( 4, 1.3   , 0.01, 0.5, 1.75);
  f.SetVar( 5, 1.3   , 0.01, 0.5, 1.75);
  f.SetVar( 6, 1.1   , 0.01, 0.5, 1.75);
  f.SetVar( 7, 0.0   , 0.01, -1, 1);

  f.SetVar( 8, 0.7   , 0.01, 0.1, 0.9);

  f.SetVar( 9, 8.0e-3, 0.0001, 0.001,  0.0142);
  f.SetVar(10, 7.0e-3, 0.0001, 0.001,  0.0142);
  f.SetVar(11, 7.0   , 0.01, 6.5, 8.5);
  f.SetVar(12, 0.0   , 0.01, -1, 1);

  f.SetVar(13, 1.3   , 0.01, 0.5, 1.75);
  f.SetVar(14, 1.3   , 0.01, 0.5, 1.75);
  f.SetVar(15, 1.1   , 0.01, 0.5, 1.75);
  f.SetVar(16, 0.0   , 0.01, -1, 1);

  f.SetVar(17, 0.7   , 0.01, 0.1, 0.9);

  f.SetVar(18, 0.0   , 1e-7, -0.1, +0.1);
  f.SetVar(19, -62.6e-6, 1e-6, -0.1, +0.1);

  // f.SetVar(20, 0.080 , 0.001, 0.078, 0.084);
  // f.SetVar(21, 0.360 , 0.001, 0.356, 0.364);
  // f.SetVar(22, 2.0 , 0.1, -2, 4);

  // f.SetVar(23, 0.080 , 0.001, 0.078, 0.084);
  // f.SetVar(24, 0.360 , 0.001, 0.356, 0.364);
  // f.SetVar(25, 2.0 , 0.1, -2, 4);

  // f.SetVar(26, 2.2e-6 , 1e-4, 0.0, 0.001);

  // f.SetVar(27, 1 , 0.001, 0.8, 1.2);
  // f.SetVar(28, 1 , 0.001, 0.8, 1.2);

  f.SetVar(20, 0 , 0.001, -0.078, 0.084);
  f.SetVar(21, 0 , 0.001, -0.356, 0.364);
  f.SetVar(22, 0 , 0.1, -2, 4);

  f.SetVar(23, 0 , 0.001, -0.078, 0.084);
  f.SetVar(24, 0 , 0.001, -0.356, 0.364);
  f.SetVar(25, 0 , 0.1, -2, 4);

  f.SetVar(26, 3.7e-09, 1e-9, 0.0, 1e-3);

  f.SetVar(27, 1 , 0.001, 0.8, 1.3);
  f.SetVar(28, 1 , 0.001, 0.8, 1.3);

  TCut cut = "modelPar.k>1.0 && modelPar.k<1.13 && abs(beamSep.X)<0.019 && abs(beamSep.Y)<0.019";
  TFile *fr = TFile::Open("root/5568/rates_5568.root");
  for (Int_t i=0,j=0; i<4; ++i) {
    if (options.Contains("Scan1") && i<2) {
      TFile::Open(momentFileNames[i]);
      TTree *t = (TTree*)gFile->Get("TBeamSpot");
      TGraph *g = (TGraph*)fr->Get(Form("gRateT0_scan%d", i));
      f.Add(j++, t, cut, g);
    }
    if (options.Contains("Scan2") && i>=2) {
      TFile::Open(momentFileNames[i]);
      TTree *t = (TTree*)gFile->Get("TBeamSpot");
      TGraph *g = (TGraph*)fr->Get(Form("gRateT0_scan%d", i));
      f.Add(j++, t, cut, g);
    }
  }
  f.SetFitToRates(!kFALSE);
  // f.GetMinimizer().FixVariable(3);
  // f.GetMinimizer().FixVariable(7);
  // f.GetMinimizer().FixVariable(12);
  // f.GetMinimizer().FixVariable(16);

  f.GetMinimizer().FixVariable(20);
  f.GetMinimizer().FixVariable(21);
  f.GetMinimizer().FixVariable(22);
  f.GetMinimizer().FixVariable(23);
  f.GetMinimizer().FixVariable(24);
  f.GetMinimizer().FixVariable(25);

  f.GetMinimizer().FixVariable(27);
  f.GetMinimizer().FixVariable(28);

  f.DoFit(TString::Format("root/%d/fit_%s.root", 5568, options.Data()),
          options);

}

void MakeNonseparationFit_5568() {
  gROOT->LoadMacro("AliDoubleGaussianBeamProfile.cxx+");
  gROOT->LoadMacro("AliNonseparationModelFit.cxx+");
  gROOT->LoadMacro("MakePlots.C");
  gROOT->LoadMacro("ExtractFromCanvas.C");

  const char* scans[3] = {
    "Scan1",
    "Scan2",
    "Scan1_Scan2",
  };
  const char* var[3]= {
    "",
    "_fix_alpha_XZ",
    "_fix_alpha_XZ_fix_alpha_YZ"
  };

#if 0
  for (Int_t i=0; i<3; ++i) {
    MakeFit(Form("Scan1%s", var[i]));
    MakeFit(Form("Scan2%s", var[i]));
    MakeFit(Form("Scan1_Scan2%s", var[i]));
  }
#endif
  TCanvas *c1 = new TCanvas;
  c1->Divide(3,3);
  for (Int_t i=0; i<3; ++i) {
    for (Int_t j=0; j<3; ++j) {
      MakePlots(Form("root/5568/fit_%s%s.root", scans[i], var[j]),
                "ALICE Pb-p #sqrt{#it{s}_{NN}}=8.2 TeV",
                0.2,
                kTRUE);
      TString line = ExtractFromCanvas(Form("pdf/5568/fit_%s%s.pdf_canvas.root", scans[i], var[j]));
      std::ofstream ofs(Form("pdf/5568/fit_%s%s_XZ.pdf_canvas.txt", scans[i], var[j]));
      ofs << line.Data() << std::endl;
       c1->cd(1+3*i+j);
      TFile::Open(Form("pdf/5568/fit_%s%s.pdf_canvas.root", scans[i], var[j]));
      Double_t *r = R->GetMatrixArray();
      Double_t r0 = TMath::Mean(20, r);
      Double_t sr0 = TMath::RMS(20, r);
      TH1* h = new TH1D(TString::Format("h_%s%s", scans[i], var[j]), TString::Format("%s%s;R", scans[i],var[j]), 20, r0-10*sr0, r0+10*sr0);
      for (Int_t k=0; k<R->GetNoElements(); ++k)
        h->Fill(r[k]);
      h->SetLabelSize(0.06, "XY");
      h->SetTitleSize(0.06, "XY");
      h->SetTitleOffset(0.7, "XY");
      h->Draw();
      gPad->Update();
      TPaveStats *ps = (TPaveStats*)h->FindObject("stats");
      ps->SetX1NDC(0.6);
      ps->SetX2NDC(0.9);
      ps->SetY1NDC(0.6);
      ps->SetY2NDC(0.9);
      ps->Draw();
      TLine *lline = new TLine;
      lline->SetLineColor(kRed);
      lline->SetLineWidth(2);
      //      line->SetLineStyle(2);
      lline->DrawLine(r[0], 0, r[0], h->GetMaximum());
    }
  }
  c1->SaveAs(Form("pdf/5568/R.pdf"));
}
