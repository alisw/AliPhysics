{
  //load libraries
  gROOT->Macro("initGlauberMC.C");

  //run the example code:
  //AliGlauberMC::runAndSaveNucleons(10000,"Pb","Pb",72);
  AliGlauberMC::runAndSaveNtuple(10000,"Pb","Pb",72);
}
