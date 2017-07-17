// -*- C++ -*-


const char* momentFileNames[] = {
  "root/4690/lumiRegion_Scan1X.root",
  "root/4690/lumiRegion_Scan1Y.root",

  "root/4690/lumiRegion_Scan2X.root",
  "root/4690/lumiRegion_Scan2Y.root",

  "root/4690/lumiRegion_ScanOffsetX.root",
  "root/4690/lumiRegion_ScanOffsetY.root",
};

void MakeNonseparationFit_4690() {
  gROOT->LoadMacro("AliDoubleGaussianBeamProfile.cxx+");
  gROOT->LoadMacro("AliNonseparationModelFit.cxx+");
  gROOT->LoadMacro("MakePlots.C");
  gROOT->LoadMacro("ExtractFromCanvas.C");
#if 1
  AliNonseparationModelFit f;

  f.SetVar( 0, 15e-4, 1e-4, 1e-4, 80e-4);
  f.SetVar( 1, 15e-4, 1e-4, 1e-4, 80e-4);
  f.SetVar( 2, 8.0   , 0.01, 2.5, 10.5);
  f.SetVar( 3, 0.0   , 0.01, -1, 1);

  f.SetVar( 4, 1.1   , 0.01, 0.75, 1.75);
  f.SetVar( 5, 1.1   , 0.01, 0.75, 1.75);
  f.SetVar( 6, 1.1   , 0.01, 0.75, 1.75);
  f.SetVar( 7, 0.0   , 0.01, -1, 1);

  f.SetVar( 8, 0.7   , 0.01, 0.5, 0.9);

  f.SetVar( 9, 15e-4, 1e-4, 1e-4, 80e-4);
  f.SetVar(10, 15e-4, 1e-4, 1e-4, 80e-4);
  f.SetVar(11, 8.0   , 0.01, 2.5, 10.5);
  f.SetVar(12, 0.0   , 0.01, -1, 1);

  f.SetVar(13, 1.1   , 0.01, 0.75, 1.75);
  f.SetVar(14, 1.1   , 0.01, 0.75, 1.75);
  f.SetVar(15, 1.1   , 0.01, 0.75, 1.75);
  f.SetVar(16, 0.0   , 0.01, -1, 1);

  f.SetVar(17, 0.7   , 0.01, 0.5, 0.9);

  f.SetVar(18, 0.0   , 1e-7, -0.1, +0.1);
  f.SetVar(19, 60e-6, 1e-6, -0.1, +0.1);

  f.SetVar(20, 0.0691 , 0.001, 0.063, 0.074);
  f.SetVar(21, 0.3300 , 0.001, 0.328, 0.335);
  f.SetVar(22, 0.2 , 0.1, -2, 2);

  f.SetVar(23, 0.0691 , 0.001, 0.063, 0.074);
  f.SetVar(24, 0.3300 , 0.001, 0.328, 0.335);
  f.SetVar(25, 0.13 , 0.01, -2, 2);

  f.SetVar(26, 2.2e-6 , 1e-4, 0.0, 0.001);

  f.SetVar(27, 1 , 0.001, 0.8, 1.2);
  f.SetVar(28, 1 , 0.001, 0.8, 1.2);

  TCut cut("modelPar.k>0.8 && modelPar.k<1.5");
  TCut cutSep[6] = {
    TCut("abs(beamSep.X)<0.011"),
    TCut("abs(beamSep.Y)<0.010"),
    TCut("abs(beamSep.X)<0.011"),
    TCut("abs(beamSep.Y)<0.010"),
    TCut("abs(beamSep.X)<0.009"),
    TCut("abs(beamSep.Y)<0.008")
  };
  for (Int_t i=0; i<6; ++i) {
    TFile::Open(momentFileNames[i]);
    TTree *t = (TTree*)gFile->Get("TBeamSpot");
    f.Add(i, t, cut*cutSep[i], NULL);
  }

  f.SetFitToRates(kFALSE);
  f.GetMinimizer().FixVariable(3);
  f.GetMinimizer().FixVariable(7);
  f.GetMinimizer().FixVariable(12);
  f.GetMinimizer().FixVariable(16);

  f.DoFit(TString::Format("root/%d/par_with0TVX_0_1.root", 4690));
#endif
  TString s = "";

  Int_t bcid=-1, wRange=0;
  Bool_t fixCrossingAngles = kTRUE;
  MakePlots("root/4690/par_with0TVX_0_1.root",
            "ALICE Pb-Pb #sqrt{#it{s}_{NN}}=5.02 TeV",
            0.15,
            kTRUE);
  TString fn = TString::Format("pdf/%d/par_with0TVX.pdf_canvas.root", 4690, s.Data(), wRange, fixCrossingAngles);
  TString line = ExtractFromCanvas(fn);
  ofs << line.Data() << std::endl;
}
