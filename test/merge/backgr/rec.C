void rec() {

  AliReconstruction reco;
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();

  reco.SetDefaultStorage("local://$ALICE_ROOT");
  reco.SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s",gSystem->pwd()));
  reco.SetRecoParam("ITS",AliITSRecoParam::GetHighFluxParam()); // to change the default vertexer
 
  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
