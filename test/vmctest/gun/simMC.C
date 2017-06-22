// $Id$
//
// Macro for running MC simulation only in test/vmctest/gun.
// Adapted from from sim.C. 

void simMC(Int_t nev=3, const TString& config) {
  if (gSystem->Getenv("EVENT"))
   nev = atoi(gSystem->Getenv("EVENT")) ;   
  
  AliSimulation simulator(config);
  simulator.SetMakeSDigits("");
  simulator.SetMakeDigits("");
  simulator.SetMakeDigitsFromHits("");
  simulator.SetRunQA(":") ; 
  simulator.SetRunHLT("") ; 

  simulator.SetDefaultStorage("local://$ALIROOT_OCDB_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));
  
  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
