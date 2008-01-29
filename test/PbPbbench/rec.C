void rec() {
  AliReconstruction reco;

  reco.SetUniformFieldTracking(kFALSE);
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetHighFluxParam();
  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  //  AliTPCReconstructor::SetStreamLevel(1);
  reco.SetRunVertexFinderTracks(kFALSE);

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
