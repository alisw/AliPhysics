void loadlibs () 
{
  gSystem->Load("libPhysics");

  // Uncomment the following line for Darwin
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
  //  gSystem->Load("libEMCAL");
  gSystem->Load("libCONTAINERS");

  // The following lines have to be commented on Darwin
  // for the moment due to cross dependencies
  gSystem->Load("libTPC");
  gSystem->Load("libITS");
  gSystem->Load("libTPCBarrel");
  gSystem->Load("libTRD");
  gSystem->Load("libTOF");
}
