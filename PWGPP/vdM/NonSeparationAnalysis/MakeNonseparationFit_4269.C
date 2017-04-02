// -*- C++ -*-


const char* momentFileNames[] = {
  "root/4269/lumiRegion_Scan1X.root",
  "root/4269/lumiRegion_Scan1Y.root",

  "root/4269/lumiRegion_Scan2X.root",
  "root/4269/lumiRegion_Scan2Y.root",

  "root/4269/lumiRegion_ScanOffsetX.root",
  "root/4269/lumiRegion_ScanOffsetY.root"
};

TVectorD MakeStartParameters() {
  TVectorD par(26);

  par[ 0] = 0.008; par[ 1] = 0.007; par[ 2] = 7.0; par[ 3] =  0.0;
  par[ 4] = 1.3;   par[ 5] = 1.3;   par[ 6] = 1.1; par[ 7] = -0.0;

  par[ 8] = 0.7;

  par[ 9] = 0.008; par[10] = 0.007; par[11] = 7.0; par[12] = -0.0;
  par[13] = 1.3;   par[14] = 1.3;   par[15] = 1.1; par[16] =  0.0;

  par[17] = 0.7;

  par[18] = 0.0; par[19] = 195e-6;

  par[20] = 0.076; par[21] = 0.537; par[22] = -0.918;

  par[23] = 2.2e-06;

  par[24] = 1; par[25] = 1;

  return par;
}

void MakeNonseparationFit_4269() {
  gROOT->LoadMacro("AliDoubleGaussianBeamProfile.cxx+");
  gROOT->LoadMacro("AliNonseparationModelFit.cxx+");

//   TVectorD muOffsetsX(3);
//   muOffsetsX[0] =  +3.76; // #mum
//   muOffsetsX[1] = +12.70; // #mum
//   muOffsetsX[2] =  -9.75; // #mum

//   TVectorD muOffsetsY(3);
//   muOffsetsY[0] =  -5.80; // #mum
//   muOffsetsY[1] =  -7.56; // #mum
//   muOffsetsY[2] =  -7.87; // #mum

  AliNonseparationModelFit f;
  // without T0 rates there is only a good fit only without the offset scans
  for (Int_t i=0; i<4; ++i) {
    TFile::Open(momentFileNames[i]);
    TTree *t = (TTree*)gFile->Get("TBeamSpot");
    f.Add(i, t, NULL);
  }

  TVectorD par = MakeStartParameters();
  f.DoFit(par, TString::Format("root/%d/par.root", 4269));
}
