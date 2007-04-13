void rec() {
  AliReconstruction reco;
  reco.SetUniformFieldTracking(kFALSE);
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  AliTPCReconstructor::SetStreamLevel(1);
  AliTPCReconstructor::SetRecoParam(AliTPCRecoParam::GetLowFluxParam());
  //  reco.SetInput("raw.root");

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
