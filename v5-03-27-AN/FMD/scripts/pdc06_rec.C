void rec() {
  AliReconstruction reco;
  reco.SetUniformFieldTracking(kFALSE);
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  AliTPCReconstructor::SetCtgRange(2.); // for pp events
  AliTPCReconstructor::SetStreamLevel(1);
  reco.SetDefaultStorage("alien://Folder=/alice/data/2006/LHC06a/CDB");
  //  reco.SetInput("raw.root");

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
