void rec() {
  new AliRun("gAlice","The ALICE Off-line Simulation Framework");

  AliReconstruction reco;


  reco.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  reco.SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s",gSystem->pwd()));
  reco.SetRunPlaneEff(kTRUE);

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
