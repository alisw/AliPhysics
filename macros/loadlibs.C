void loadlibs () 
{
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");

  // libraries required by EVGEN
  // (commented libraries are already loaded in prevoius sequence)

  gSystem->Load("libmicrocern");
  gSystem->Load("libEG"); 
  gSystem->Load("libSTEER");
  gSystem->Load("libEVGEN");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libpdf");
  gSystem->Load("libpythia6");
  gSystem->Load("libAliPythia6");

  gSystem->Load("libPhysics");

  gSystem->Load("libCONTAINERS");
  gSystem->Load("libFMD");
  gSystem->Load("libMUON");
  gSystem->Load("libPHOS");
  gSystem->Load("libPMD");
  gSystem->Load("libRICH");
  gSystem->Load("libSTRUCT");
  gSystem->Load("libTPC");
  gSystem->Load("libTRD");
  gSystem->Load("libTOF");
  gSystem->Load("libZDC");
  gSystem->Load("libITS");
  gSystem->Load("libCRT");
  gSystem->Load("libSTART");
//  gSystem->Load("libEMCAL");
  gSystem->Load("libVZERO");
}
