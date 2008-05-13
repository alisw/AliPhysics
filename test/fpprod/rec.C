void rec() {
  AliReconstruction reco;
  reco.SetUniformFieldTracking(kFALSE);
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetLowFluxParam();
  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  AliTPCReconstructor::SetStreamLevel(1);
  //   reco.SetInput("raw.root");
  reco.SetRunReconstruction("ITS TPC TRD TOF HMPID PHOS EMCAL MUON VZERO T0 FMD PMD ZDC");
  reco.SetRunQA(kFALSE);
  reco.SetRunGlobalQA(kFALSE);

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
