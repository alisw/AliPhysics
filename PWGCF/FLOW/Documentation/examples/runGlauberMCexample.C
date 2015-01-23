{
  //load libraries
 gSystem->Load("libPWGGlauber");
 //  gSystem->SetBuildDir("/tmp");
 //  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/Glauber/AliGlauberNucleon.cxx+");
 //  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/Glauber/AliGlauberNucleus.cxx+");
 //  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/Glauber/AliGlauberMC.cxx+");

  //set the random seed from current time
  TTimeStamp time;
  Int_t seed = time->GetSec();
  gRandom->SetSeed(seed);

  Int_t nevents = 1000; // number of events to simulate 
  // supported systems are e.g. "p", "d", "Si", "Au", "Pb", "U" 
  Option_t *sysA="Pb"; 
  Option_t *sysB="Pb";
  Double_t signn=64; // inelastic nucleon nucleon cross section
  //const char *fname="GlauberMC_PbPb_ntuple.root"; // name output file

  // run the code to produce an ntuple:
  //  AliGlauberMC::runAndSaveNucleons(10000,"Pb","Pb",72);
  Double_t mind=0.4;
  AliGlauberMC::RunAndSaveNtuple(nevents,sysA,sysB,signn,mind);

}
