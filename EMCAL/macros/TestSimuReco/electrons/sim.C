///
/// Example of simulation macro to generate electrons in EMCAL 
/// and TCP/ITS. Adapted from $ALICE_ROOT/test/gun
///
void sim( Int_t nev = 10 ) 
{
  // libraries required by geant321
  gSystem->Load("liblhapdf");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libpythia6");
  gSystem->Load("libAliPythia6");
  gSystem->Load("libgeant321");

  AliSimulation simulator;
  simulator.SetMakeSDigits("EMCAL");
  simulator.SetMakeDigitsFromHits("ITS TPC");

  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
