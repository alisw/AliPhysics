void rec() {
  AliReconstruction reco;
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();

  // Use the GRP from the backgr
  reco.SetDefaultStorage("local://$ALICE_ROOT");
  reco.SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s/../backgr",gSystem->pwd()));
  reco.SetRecoParam("ITS",AliITSRecoParam::GetHighFluxParam()); // to change the default vertexer

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
