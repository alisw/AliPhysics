// -*- C++ -*-


TVectorD MakeStartParameters() {
  TVectorD par(26);

  par[ 0] = 0.0118; par[ 1] = 0.0100; par[ 2] = 8.1; par[ 3] = 0.0;
  par[ 4] = 1.37;   par[ 5] = 1.58;   par[ 6] = 0.8; par[ 7] = 0.0;
  par[ 8] = 0.82;
  
  par[ 9] = 0.0112; par[10] = 0.0107; par[11] = 6.8; par[12] = 0.0;
  par[13] = 1.33;   par[14] = 1.53;   par[15] = 1.1; par[16] = 0.0;
  par[17] = 0.71;

  par[18] = 0.0; par[19] = 360e-6;

  par[20] = 0.0745; par[21] = 0.335; par[22] = 0.451;

  par[23] = 9.1e-8;

  par[24] = 1; par[25] = 1;
  return par;
}

const Int_t bcs[] = {
  -1,   // -B
  344,  // -I1 BCM5   344H1L3219H
  464,  // -I2 BCM6   464H1L3099H
  827,  // -I3 BCM7   827H1L2736H
  1187, // -I4 BCM8  1187H1L2376H
  1558, // -I5 BCM9  1558H1L2005H
  1678, // -I6 BCM10 1678H1L1885H
  3177, // -I7 BCM11 3177H1L386H
  3297  // -I8 BCM12 3297H1L266H
};


const char* bcNames[] = {
  "-B",
  "-I1", "-I2", "-I3", "-I4", "-I5", "-I6", "-I7", "-I8"
};
const char* momentFileNames[] = {
  "root/4634/lumiRegion_Scan1X%s.root",
  "root/4634/lumiRegion_Scan1Y%s.root",

  "root/4634/lumiRegion_Scan2X%s.root",
  "root/4634/lumiRegion_Scan2Y%s.root",

  "root/4634/lumiRegion_ScanOffsetX%s.root",
  "root/4634/lumiRegion_ScanOffsetY%s.root"
};


void MakeNonseparationFit_4634(Bool_t useRates, Int_t bcIndex) {
  gROOT->LoadMacro("AliDoubleGaussianBeamProfile.cxx+");
  gROOT->LoadMacro("AliNonseparationModelFit.cxx+");

  // Scans 1 and 2
  TGraph *gRate[6];
  TFile::Open(Form("root/4634/244369_C0TVX%s.root", bcNames[bcIndex]));
  for (Int_t i=0; i<4; ++i)
    gRate[i] = (TGraph*)gFile->Get(Form("gRateVsSep_Scan%d", i));

  // Offset scan
  TFile::Open(Form("root/4634/244375_C0TVX%s.root", bcNames[bcIndex]));
  for (Int_t i=4; i<6; ++i)
    gRate[i] = (TGraph*)gFile->Get(Form("gRateVsSep_Scan%d", i-4));

  const TString bcString = (bcIndex == 0 ? "" : Form("_bcid%d", bcs[bcIndex]));

  AliNonseparationModelFit f;
  TTree *t(NULL);
  for (Int_t i=0; i<6; ++i) {
    TFile::Open(Form(momentFileNames[i], bcString.Data()));
    t = (TTree*)gFile->Get("TBeamSpot");
    f.Add(i, t, useRates ? gRate[i] : NULL);
  }

  TVectorD par = MakeStartParameters();
  f.DoFit(par, TString::Format("root/%d/par_%s%s.root", 4634, useRates ? "with0TVX" : "withoutRates", bcString.Data()));
}
