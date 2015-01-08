void simTGeoUpgrade(Int_t nev=10) {

  gSystem->Exec(" rm *.root ");
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");
  AliSimulation simulator("ConfigTgeoUpgrade.C");
  simulator.SetMakeSDigits("");

  simulator.SetMakeDigits("");

  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));
  simulator.SetRunHLT("");
  simulator.SetRunQA(":");

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
