void loadlibs () 
{
  gSystem->Load("libMC");

  // libraries required by EVGEN
  // (commented libraries are already loaded in prevoius sequence)

  gSystem->Load("libmicrocern");
  gSystem->Load("libSTEER");
  gSystem->Load("libEG"); 
  gSystem->Load("libEGPythia6");
  gSystem->Load("libdummypythia6");
  gSystem->Load("libdummyhijing");
  gSystem->Load("libTHijing");
  gSystem->Load("libdummyHBTP");
  gSystem->Load("libTHbtp");
  gSystem->Load("libdummymevsim");
  gSystem->Load("libTMevSim");
  gSystem->Load("libEVGEN");

  gSystem->Load("libPhysics");

  gSystem->Load("libCONTAINERS");
  gSystem->Load("libFMD");
  gSystem->Load("libMUON");
  gSystem->Load("libPHOS");
  gSystem->Load("libPMD");
  gSystem->Load("libRICH");
  gSystem->Load("libSTRUCT");
  gSystem->Load("libTOF");
  gSystem->Load("libTPC");
  gSystem->Load("libTRD");
  gSystem->Load("libZDC");
  gSystem->Load("libITS");
  gSystem->Load("libCRT");
  gSystem->Load("libSTART");
  gSystem->Load("libEMCAL");
  gSystem->Load("libVZERO");
  gSystem->Load("libdummyherwig");
  gSystem->Load("libTHerwig");
}
