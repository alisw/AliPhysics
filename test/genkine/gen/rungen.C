void rungen(Int_t nev=1){
  // Simulation and reconstruction
  TStopwatch timer;
  timer.Start();
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT");
  gROOT->LoadMacro("fastGen.C+");
  fastGen(nev);
  timer.Stop();
  timer.Print();
}
