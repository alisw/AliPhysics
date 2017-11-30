// -*- C++ -*-


const char* momentFileNames[] = {
  "root/4269/lumiRegion_Scan1X.root",
  "root/4269/lumiRegion_Scan1Y.root",

  "root/4269/lumiRegion_Scan2X.root",
  "root/4269/lumiRegion_Scan2Y.root",

  "root/4269/lumiRegion_ScanOffsetX.root",
  "root/4269/lumiRegion_ScanOffsetY.root"
};

void MakeNonseparationFit_4269() {
  gROOT->LoadMacro("AliDoubleGaussianBeamProfile.cxx+");
  gROOT->LoadMacro("AliNonseparationModelFit.cxx+");
  gROOT->LoadMacro("MakePlots.C");
  gROOT->LoadMacro("ExtractFromCanvas.C");
#if 1
  AliNonseparationModelFit f;

  f.SetVar( 0, 8.0e-3, 0.0001, 0.001,  0.0142);
  f.SetVar( 1, 7.0e-3, 0.0001, 0.001,  0.0142);
  f.SetVar( 2, 7.0   , 0.01, 6.5, 8.5);
  f.SetVar( 3, 0.0   , 0.01, -1, 1);

  f.SetVar( 4, 1.3   , 0.01, 0.75, 1.5);
  f.SetVar( 5, 1.3   , 0.01, 0.75, 1.5);
  f.SetVar( 6, 1.1   , 0.01, 0.75, 1.5);
  f.SetVar( 7, 0.0   , 0.01, -1, 1);

  f.SetVar( 8, 0.7   , 0.01, 0.1, 0.9);

  f.SetVar( 9, 8.0e-3, 0.0001, 0.001,  0.0142);
  f.SetVar(10, 7.0e-3, 0.0001, 0.001,  0.0142);
  f.SetVar(11, 7.0   , 0.01, 6.5, 8.5);
  f.SetVar(12, 0.0   , 0.01, -1, 1);

  f.SetVar(13, 1.3   , 0.01, 0.75, 1.5);
  f.SetVar(14, 1.3   , 0.01, 0.75, 1.5);
  f.SetVar(15, 1.1   , 0.01, 0.75, 1.5);
  f.SetVar(16, 0.0   , 0.01, -1, 1);

  f.SetVar(17, 0.7   , 0.01, 0.1, 0.9);

  f.SetVar(18, 0.0   , 1e-7, -0.1, +0.1);
  f.SetVar(19, 195e-6, 1e-6, -0.1, +0.1);

  f.SetVar(20, 0.076 , 0.001, 0.072, 0.080);
  f.SetVar(21, 0.537 , 0.001, 0.532, 0.541);
  f.SetVar(22, -0.9 , 0.1, -2, 2);

  f.SetVar(23, 0.076 , 0.001, 0.072, 0.080);
  f.SetVar(24, 0.537 , 0.001, 0.532, 0.541);
  f.SetVar(25, -0.9 , 0.1, -2, 2);

  f.SetVar(26, 2.2e-6 , 1e-4, 0.0, 0.001);

  f.SetVar(27, 1 , 0.001, 0.8, 1.2);
  f.SetVar(28, 1 , 0.001, 0.8, 1.2);

  TCut cut = "modelPar.k>0.9 && modelPar.k<1.15";
  for (Int_t i=0; i<6; ++i) {
    TFile::Open(momentFileNames[i]);
    TTree *t = (TTree*)gFile->Get("TBeamSpot");
    if (i<4)
      f.Add(i, t, cut, NULL);
    else
      f.Add(i, t, cut*TCut(), NULL);
  }

  f.SetFitToRates(kFALSE);
  f.GetMinimizer().FixVariable(3);
  f.GetMinimizer().FixVariable(7);
  f.GetMinimizer().FixVariable(12);
  f.GetMinimizer().FixVariable(16);

  f.DoFit(TString::Format("root/%d/par_with0TVX_0_1.root", 4269));
#endif
  TString s = "";

  Int_t bcid=-1, wRange=0;
  Bool_t fixCrossingAngles = kTRUE;
  MakePlots("root/4269/par_with0TVX_0_1.root",
            "ALICE pp #sqrt{#it{s}}=13 TeV",
            0.6,
            kTRUE);
  TString fn = TString::Format("pdf/%d/par_with0TVX.pdf_canvas.root", 4269, s.Data(), wRange, fixCrossingAngles);
  TString line = ExtractFromCanvas(fn);
  ofs << line.Data() << std::endl;
}
