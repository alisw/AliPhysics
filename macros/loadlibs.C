void loadlibs () 
{
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");

  // Uncomment the following line for macosx
  // Waiting for a better solution
  // gSystem->Load("libg2c_sh");
  gSystem->Load("libmicrocern");
  gSystem->Load("libpdf");
  gSystem->Load("libpythia6");

  gSystem->Load("libEG");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");

  gSystem->Load("libEGPythia6");
  gSystem->Load("libSTEER");
  gSystem->Load("libEVGEN");
  gSystem->Load("libFASTSIM");
  gSystem->Load("libAliPythia6");

  gSystem->Load("libRAW");

  gSystem->Load("libSTRUCT");
  gSystem->Load("libPHOS");
  gSystem->Load("libMUON");
  gSystem->Load("libFMD");
  gSystem->Load("libPMD");
  gSystem->Load("libRICH");
  gSystem->Load("libSTART");
  gSystem->Load("libZDC");
  gSystem->Load("libCRT");
  gSystem->Load("libVZERO");
  gSystem->Load("libEMCAL");
  gSystem->Load("libCONTAINERS");

  // The following lines have to be commented on Darwin
  // for the moment due to cross dependencies
  gSystem->Load("libTPCbase");
  gSystem->Load("libTPCrec");
  gSystem->Load("libTPCsim");
  gSystem->Load("libTPCfast");
  gSystem->Load("libITS");
  gSystem->Load("libTPCBarrel");
  gSystem->Load("libTRDbase");
  gSystem->Load("libTRDsim");
  gSystem->Load("libTRDrec");
  gSystem->Load("libTRDfast");
  gSystem->Load("libTOF");
}
