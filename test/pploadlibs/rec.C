void rec() {
  gSystem->Load("libTree");
  gROOT->Macro("loadlibsrec.C");
  new AliRun("gAlice","The ALICE Off-line Simulation Framework");
  AliReconstruction reco;
  reco.SetUniformFieldTracking(kFALSE);
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetLowFluxParam();
  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  AliTPCReconstructor::SetStreamLevel(1);
  //   reco.SetInput("raw.root");

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
