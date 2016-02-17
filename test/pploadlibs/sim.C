void sim(Int_t nev=20) {
  new AliRun("gAlice","The ALICE Off-line Simulation Framework");

  AliSimulation simulator;
  simulator.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");
    simulator.SetMakeDigitsFromHits("ITS TPC");
 
  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));
  simulator.SetRunQA(":");

  simulator.SetRunHLT("default"); // In case we do not have ancored production

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
