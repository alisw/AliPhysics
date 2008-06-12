void rec(Int_t embrun=0) {
  new AliRun("gAlice","The ALICE Off-line Simulation Framework");
  AliReconstruction reco;
  reco.SetUniformFieldTracking(kFALSE);
  reco.SetWriteESDfriend(kTRUE);
  reco.SetWriteAlignmentData(kFALSE);
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetLowFluxParam();
  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  AliTPCReconstructor::SetStreamLevel(0);
  reco.SetRunReconstruction("ITS TPC TRD TOF");
  reco.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-12-Release/Ideal/");
  reco.SetRunQA(kFALSE);
  reco.SetRunGlobalQA(kFALSE);

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
