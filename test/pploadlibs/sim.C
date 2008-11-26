void sim(Int_t nev=20) {
  gROOT->Macro("loadlibssim.C");
  new AliRun("gAlice","The ALICE Off-line Simulation Framework");

  AliSimulation simulator;
  simulator.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");
  simulator.SetMakeDigitsFromHits("ITS TPC");
 
  simulator.SetDefaultStorage("local://$ALICE_ROOT");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));
  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
