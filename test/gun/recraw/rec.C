void rec() {

  AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1"); // Do nkt forget to keep consistency with Config.C
  
  AliReconstruction reco;

  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();

  reco.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  reco.SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s/..",gSystem->pwd()));
  reco.SetRunPlaneEff(kTRUE);

  reco.SetInput("raw.root");

  reco.SetNumberOfEventsPerFile(-1); // all events in one single file

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
