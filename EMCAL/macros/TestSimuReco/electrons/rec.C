///
/// Example of reconstruction macro for electrons in EMCAL
/// and TCP/ITS.
/// Adapted from $ALICE_ROOT/test/gun.
///
void rec() 
{
  AliReconstruction reco;

  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();

  reco.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  reco.SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s",gSystem->pwd()));
  reco.SetRunPlaneEff(kTRUE);

  reco.SetFractionFriends(1.);

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
