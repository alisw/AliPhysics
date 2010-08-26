{
  //load libraries
  gROOT->LoadMacro("AliGlauberNucleon.cxx+");
  gROOT->LoadMacro("AliGlauberNucleus.cxx+");
  gROOT->LoadMacro("AliGlauberMC.cxx+");


  Int_t nevents = 10000; // number of events to simulate 
  // supported systems are e.g. "p", "d", "Si", "Au", "Pb", "U" 
  Option_t *sysA="Pb"; 
  Option_t *sysB="Pb";
  Double_t signn=72; // inelastic nucleon nucleon cross section
  const char *fname="GlauberMC_PbPb_ntuple.root"; // name output file

  // run the code to produce an ntuple:
  //  AliGlauberMC::runAndSaveNucleons(10000,"Pb","Pb",72);
  Double_t mind=0.4;
  AliGlauberMC::runAndSaveNtuple(nevents,sysA,sysB,signn,mind,fname);

}
