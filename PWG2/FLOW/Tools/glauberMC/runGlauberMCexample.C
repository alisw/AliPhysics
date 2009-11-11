{
  //load libraries
  gROOT->Macro("initGlauberMC.C");

  //run the example code:
  AliGlauberMC::runAndSaveNucleons(100);
  AliGlauberMC::runAndSaveNtuple(100);
}
