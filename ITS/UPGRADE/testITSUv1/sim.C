void sim(Int_t nev=3) {

  gSystem->Exec(" rm itsSegmentations.root ");
  AliSimulation simulator;
  //  simulator.SetMakeSDigits("");

  //  simulator.SetMakeDigits("");

  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));
  simulator.SetSpecificStorage("ITS/Align/Data",
			       Form("local://%s",gSystem->pwd()));
  simulator.SetSpecificStorage("ITS/Calib/SimuParam",
			       Form("local://%s",gSystem->pwd()));
  simulator.SetRunHLT("");
  simulator.SetRunQA(":");

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
