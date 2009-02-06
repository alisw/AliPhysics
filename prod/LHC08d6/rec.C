void rec() {

  AliReconstruction reco;
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();

  reco.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/?cacheFold=/tmp/CDBCache?operateDisconnected=kFALSE");

  // No write access to the OCDB => specific storage
  reco.SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s",gSystem->pwd()));

  reco.SetRunQA("ALL:ALL");

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
