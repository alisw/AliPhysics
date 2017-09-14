// -*- C++ -*-


const char* momentFileNames[] = {
  "root/5533/lumiRegion_Scan1X.root",
  "root/5533/lumiRegion_Scan1Y.root",

  "root/5533/lumiRegion_Scan2X.root",
  "root/5533/lumiRegion_Scan2Y.root",
};

void MakeNonseparationFit_5533() {
  gROOT->LoadMacro("AliDoubleGaussianBeamProfile.cxx+");
  gROOT->LoadMacro("AliNonseparationModelFit.cxx+");
  gROOT->LoadMacro("MakePlots.C");
  gROOT->LoadMacro("ExtractFromCanvas.C");
#if 1
  AliNonseparationModelFit f;

  f.SetVar( 0, 3.5e-3, 0.0001, 0.001,  0.0142);
  f.SetVar( 1, 2.5e-3, 0.0001, 0.001,  0.0142);
  f.SetVar( 2, 8.5   , 0.01, 4.5, 10.5);
  f.SetVar( 3, 0.0   , 0.01, -1, 1);

  f.SetVar( 4, 1.1   , 0.01, 0.5, 1.75);
  f.SetVar( 5, 1.1   , 0.01, 0.5, 1.75);
  f.SetVar( 6, 1.1   , 0.01, 0.5, 1.75);
  f.SetVar( 7, 0.0   , 0.01, -1, 1);

  f.SetVar( 8, 0.7   , 0.01, 0.1, 0.9);

  f.SetVar( 9, 4.0e-3, 0.0001, 0.001,  0.0142);
  f.SetVar(10, 2.5e-3, 0.0001, 0.001,  0.0142);
  f.SetVar(11, 8.5   , 0.01, 4.5, 10.5);
  f.SetVar(12, 0.0   , 0.01, -1, 1);

  f.SetVar(13, 1.1   , 0.01, 0.5, 1.75);
  f.SetVar(14, 1.1   , 0.01, 0.5, 1.75);
  f.SetVar(15, 1.1   , 0.01, 0.5, 1.75);
  f.SetVar(16, 0.0   , 0.01, -1, 1);

  f.SetVar(17, 0.7   , 0.01, 0.1, 0.9);

  f.SetVar(18, 0.0   , 1e-7, -0.1, +0.1);
  f.SetVar(19, -40e-6, 1e-6, -0.1, +0.1);

  f.SetVar(20, 0.080 , 0.001, 0.078, 0.084);
  f.SetVar(21, 0.360 , 0.001, 0.356, 0.364);
  f.SetVar(22, 2.0 , 0.1, -2, 4);

  f.SetVar(23, 0.080 , 0.001, 0.078, 0.084);
  f.SetVar(24, 0.360 , 0.001, 0.356, 0.364);
  f.SetVar(25, 2.0 , 0.1, -2, 4);

  f.SetVar(26, 2.2e-6 , 1e-4, 0.0, 0.001);

  f.SetVar(27, 1.0 , 0.001, 0.8, 1.2);
  f.SetVar(28, 1.0 , 0.001, 0.8, 1.2);

  TCut cut = "modelPar.k>1.0 && modelPar.k<1.13 && abs(beamSep.X)<0.015 && abs(beamSep.Y)<0.015";
  for (Int_t i=0; i<4; ++i) {
    TFile::Open(momentFileNames[i]);
    TTree *t = (TTree*)gFile->Get("TBeamSpot");
    f.Add(i, t, cut, NULL);
  }

  f.SetFitToRates(kFALSE);
  f.GetMinimizer().FixVariable(3);
  f.GetMinimizer().FixVariable(7);
  f.GetMinimizer().FixVariable(12);
  f.GetMinimizer().FixVariable(16);

  f.DoFit(TString::Format("root/%d/par_with0TVX_0_1.root", 5533));
#endif
  TString s = "";

  Int_t bcid=-1, wRange=0;
  Bool_t fixCrossingAngles = kTRUE;
  MakePlots("root/5533/par_with0TVX_0_1.root",
            "ALICE p-Pb #sqrt{#it{s}_{NN}}=8.2 TeV",
            0.2,
            kTRUE);
  TString fn = TString::Format("pdf/%d/par_with0TVX.pdf_canvas.root", 5533, s.Data(), wRange, fixCrossingAngles);
  TString line = ExtractFromCanvas(fn);
  ofs << line.Data() << std::endl;
}
