{
  //load libraries
 gSystem->Load("libPWG2flowTools");
 //  gSystem->SetBuildDir("/tmp");
 //  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FLOW/AliFlowTools/glauberMC/AliGlauberNucleon.cxx+");
 //  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FLOW/AliFlowTools/glauberMC/AliGlauberNucleus.cxx+");
 //  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FLOW/AliFlowTools/glauberMC/AliGlauberMC.cxx+");

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
  //  AliGlauberMC::RunAndSaveNtuple(nevents,sysA,sysB,signn,mind);
  Double_t r=6.62;
  Double_t a=0.546;
  const char *fname="glau_pbpb_ntuple.root";

  AliGlauberMC mcg(sysA,sysB,signn);
  mcg.SetMinDistance(mind);
  mcg.Setr(r);
  mcg.Seta(a);
  mcg.SetDoPartProduction(kFALSE);
  //mcg.SetDoPartProduction(kTRUE);
  mcg.Run(nevents);

  TNtuple  *nt=mcg.GetNtuple();
  TFile out(fname,"recreate",fname,9);
  if(nt) nt->Write();
  printf("total cross section with a nucleon-nucleon cross section \t%f is \t%f",signn,mcg.GetTotXSect());
  out.Close();

}
