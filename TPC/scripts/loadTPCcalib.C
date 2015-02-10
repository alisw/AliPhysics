/// \file loadTPCcalib.C

void loadTPCcalib(){
  ///

  gROOT->Macro("~/NimStyle.C");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  printf("LOAD TPC calibration libraries\n\n");
}
