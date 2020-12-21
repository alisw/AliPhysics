// -*- C++ -*-

const char* momentFileNames[] = {
  "root/6380/lumiRegion_Scan1X.root",
  "root/6380/lumiRegion_Scan1Y.root",

  "root/6380/lumiRegion_Scan2X.root",
  "root/6380/lumiRegion_Scan2Y.root",
};

void MakeFit(TString options, Int_t bc) {
  TString bcName = (bc == -1
                    ? ""
                    : TString::Format("_bcid%d", bc));

  TVirtualFitter::SetMaxIterations(100*1000);

  AliNonseparationModelFit f;
  f.SetScaleRateError(bc == -1 ? 1.0 : 1.0/(10*10));

#if 0
  f.SetVar( 0, 8e-3, 0.0001, 0.001,  0.02);
  f.SetVar( 1, 8e-3, 0.0001, 0.001,  0.02);
  f.SetVar( 2, 6.5   , 0.01, 3.5, 10.5);
  f.SetVar( 3, 0.0   , 0.01, -1, 1);

  f.SetVar( 4, 1.0   , 0.01, 0.5, 1.75);
  f.SetVar( 5, 1.0   , 0.01, 0.5, 1.75);
  f.SetVar( 6, 0.8   , 0.01, 0.5, 1.75);
  f.SetVar( 7, 0.0   , 0.01, -1, 1);

  f.SetVar( 8, 0.7   , 0.01, 0.1, 0.9);

  f.SetVar( 9, 8.0e-3, 0.0001, 0.001,  0.02);
  f.SetVar(10, 8.0e-3, 0.0001, 0.001,  0.02);
  f.SetVar(11, 6.5   , 0.01, 3.5, 10.5);
  f.SetVar(12, 0.0   , 0.01, -1, 1);

  f.SetVar(13, 1.2   , 0.01, 0.5, 1.75);
  f.SetVar(14, 1.2   , 0.01, 0.5, 1.75);
  f.SetVar(15, 0.8   , 0.01, 0.5, 1.75);
  f.SetVar(16, 0.0   , 0.01, -1, 1);

  f.SetVar(17, 0.7   , 0.01, 0.1, 0.9);
#else
  f.SetVar( 0, 12.4e-3, 0.0001, 0.001,  0.02);
  f.SetVar( 1, 10.6e-3, 0.0001, 0.001,  0.02);
  f.SetVar( 2, 8.0   , 0.01, 3.5, 10.5);
  f.SetVar( 3, 0.0   , 0.01, -1, 1);

  f.SetVar( 4, 0.86  , 0.01, 0.5, 1.75);
  f.SetVar( 5, 1.15  , 0.01, 0.5, 1.75);
  f.SetVar( 6, 0.90   , 0.01, 0.5, 1.75);
  f.SetVar( 7, 0.0   , 0.01, -1, 1);

  f.SetVar( 8, 0.71  , 0.01, 0.1, 0.9);

  f.SetVar( 9, 11.6e-3, 0.0001, 0.001,  0.02);
  f.SetVar(10, 11.3e-3, 0.0001, 0.001,  0.02);
  f.SetVar(11, 7.4   , 0.01, 3.5, 10.5);
  f.SetVar(12, 0.0   , 0.01, -1, 1);

  f.SetVar(13, 0.92  , 0.01, 0.5, 1.75);
  f.SetVar(14, 1.14  , 0.01, 0.5, 1.75);
  f.SetVar(15, 0.82  , 0.01, 0.5, 1.75);
  f.SetVar(16, 0.0   , 0.01, -1, 1);

  f.SetVar(17, 0.88  , 0.01, 0.1, 0.9);
#endif
  f.SetVar(18, 0.0   , 1e-7, -0.1, +0.1);
  f.SetVar(19, -340e-6, 1e-6, -0.1, +0.1);

  f.SetVar(20, 0 , 0.001, -0.078, 0.084);
  f.SetVar(21, 0 , 0.001, -0.356, 0.364);
  f.SetVar(22, 0 , 0.1, -2, 4);

  f.SetVar(23, 0 , 0.001, -0.078, 0.084);
  f.SetVar(24, 0 , 0.001, -0.356, 0.364);
  f.SetVar(25, 0 , 0.1, -2, 4);

  f.SetVar(26, 3.79500e-06*(bc==-1 ? 1 : 22), 1e-6, 0.0, 1e-3);

  f.SetVar(27, 1 , 0.001, 0.8, 1.3);
  f.SetVar(28, 1 , 0.001, 0.8, 1.3);

  TCut cut = (bc == -1
              ? "modelPar.k>0.9 && modelPar.k<1.2 && abs(beamSep.X)<0.06 && abs(beamSep.Y)<0.06"
              : "modelPar.k>0.7 && modelPar.k<1.4 && abs(beamSep.X)<0.06 && abs(beamSep.Y)<0.06 && abs(modelPar.rhoXY)<0.9");
  TFile *fr = TFile::Open("~/cernbox/data/vdM/6380/RatesVsSep.root");
  for (Int_t i=0,j=0; i<4; ++i) {
    if (options.Contains("Scan1") && i<2) {
      TString s = momentFileNames[i];
      if (bc != -1)
        s.ReplaceAll(".root", bcName+".root");
      TFile::Open(s);
      TTree *t = (TTree*)gFile->Get("TBeamSpot");
      TGraph *g = fr ? (TGraph*)fr->Get(Form("T0Rate_Scan%d%s", i, bcName.Data())) : NULL;
      f.Add(j++, t, cut, g);
    }
    if (options.Contains("Scan2") && i>=2) {
      TString s = momentFileNames[i];
      if (bc != -1)
        s.ReplaceAll(".root", bcName+".root");
      TFile::Open(s);
      TTree *t = (TTree*)gFile->Get("TBeamSpot");
      TGraph *g = fr ? (TGraph*)fr->Get(Form("T0Rate_Scan%d%s", i, bcName.Data())) : NULL;
      f.Add(j++, t, cut, g);
    }
  }
  f.SetFitToRates(kTRUE);

  f.GetMinimizer().FixVariable(20);
  f.GetMinimizer().FixVariable(21);
  f.GetMinimizer().FixVariable(22);
  f.GetMinimizer().FixVariable(23);
  f.GetMinimizer().FixVariable(24);
  f.GetMinimizer().FixVariable(25);

  f.GetMinimizer().FixVariable(27);
  f.GetMinimizer().FixVariable(28);
  Printf("calling dofit");
  f.DoFit(TString::Format("root/%d/fit_%s%s.root", 6380, options.Data(), bcName.Data()),
          options);
}

void MakeNonseparationFit_6380(Int_t bc=-1) {
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
    // MakeFit(Form("Scan1%s", var[i]), bc);
    // MakeFit(Form("Scan2%s", var[i]), bc);
    MakeFit(Form("Scan1_Scan2%s", var[i]), bc);
  }
#else
  TString bcName = (bc == -1
                    ? ""
                    : TString::Format("_bcid%d", bc));
  TCanvas *c1 = new TCanvas;
  c1->Divide(3,3);
  for (Int_t i=2; i<3; ++i) {
    for (Int_t j=0; j<1; ++j) {
      MakePlots(Form("root/6380/fit_%s%s%s.root", scans[i], var[j], bcName.Data()),
                "ALICE p-p #sqrt{#it{s}_{NN}}=5 TeV",
                0.75,
                !kTRUE);
      TString line = ExtractFromCanvas(Form("pdf/6380/fit_%s%s%s.pdf_canvas.root", scans[i], var[j], bcName.Data()));
      std::ofstream ofs(Form("pdf/6380/fit_%s%s%s_XZ.pdf_canvas.txt", scans[i], var[j], bcName.Data()));
      ofs << line.Data() << std::endl;
      c1->cd(1+3*i+j);
      TFile::Open(Form("pdf/6380/fit_%s%s%s.pdf_canvas.root", scans[i], var[j], bcName.Data()));
      TVectorD *R = NULL;
      gFile->GetObject("R", R);
      Double_t const* r = R->GetMatrixArray();
      Double_t r0 = TMath::Mean(R->GetNoElements(), R->GetMatrixArray());
      Double_t sr0 = TMath::RMS(R->GetNoElements(), R->GetMatrixArray());
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
  c1->SaveAs(Form("pdf/6380/R%s.pdf", bcName.Data()));
#endif
  // TString fn = TString::Format("pdf/%d/par_with0TVX.pdf_canvas.root", 6380, s.Data(), wRange, fixCrossingAngles);
  // TString line = ExtractFromCanvas(fn);
  // ofs << line.Data() << std::endl;
}
