// -*- C++ -*-


const char* momentFileNames[] = {
  "root/4634/lumiRegion_Scan1X.root",
  "root/4634/lumiRegion_Scan1Y.root",

  "root/4634/lumiRegion_Scan2X.root",
  "root/4634/lumiRegion_Scan2Y.root",

  "root/4634/lumiRegion_ScanOffsetX.root",
  "root/4634/lumiRegion_ScanOffsetY.root"
};

TVectorD MakeStartParameters() {
  TVectorD par(26);

  par[ 0] = 0.008;  par[ 1] = 0.008;  par[ 2] = 4.0; par[ 3] = -0.1;
  par[ 4] = 1.2;    par[ 5] = 1.2;    par[ 6] = 1.1; par[ 7] =  0.1;
  par[ 8] = 0.9;

  par[ 9] = 0.008;  par[10] = 0.008;  par[11] = 4.0; par[12] =  0.1;
  par[13] = 1.2;    par[14] = 1.2;    par[15] = 1.1; par[16] = -0.1;
  par[17] = 0.9;

  par[18] = 0.0; par[19] = 360e-6;

  par[20] = 0.075; par[21] = 0.335; par[22] = 0.467;

  par[23] = 2.3e-6;

  par[24] = 1; par[25] = 1;
  return par;
}

void MakeNonseparationFit_4634(Bool_t useRates=kTRUE) {
  gROOT->LoadMacro("AliDoubleGaussianBeamProfile.cxx+");
  gROOT->LoadMacro("AliNonseparationModelFit.cxx+");

  // Scans 1 and 2
  TGraph *gRate[6];
  TFile::Open("root/4634/244369_C0TVX-B.root");
  for (Int_t i=0; i<4; ++i)
    gRate[i] = (TGraph*)gFile->Get(Form("gRateVsSep_Scan%d", i));

  // Offset scan
  TFile::Open("root/4634/244375_C0TVX-B.root");
  for (Int_t i=4; i<6; ++i)
    gRate[i] = (TGraph*)gFile->Get(Form("gRateVsSep_Scan%d", i-4));

  AliNonseparationModelFit f;
  for (Int_t i=0; i<6; ++i) {
    TFile::Open(momentFileNames[i]);
    TTree *t = (TTree*)gFile->Get("TBeamSpot");
    f.Add(i, t, useRates ? gRate[i] : NULL);
  }

  TVectorD par = MakeStartParameters();
  f.DoFit(par, TString::Format("root/%d/par_%s.root", 4634, useRates ? "with0TVX" : "withoutRates"));
}
