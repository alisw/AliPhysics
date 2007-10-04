void rec() {
  AliReconstruction reco;
  reco.SetUniformFieldTracking(kFALSE);
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  AliTPCReconstructor::SetStreamLevel(1);
  //  AliTPCReconstructor::SetRecoParam(AliTPCRecoParam::GetLowFluxParam());
  //  reco.SetInput("raw.root");
  //  reco.SetRunReconstruction("ITS TPC TRD TOF HMPID PHOS EMCAL MUON T0 VZERO FMD PMD ZDC");
 
  AliCDBManager::Instance()->SetCacheFlag(kFALSE);
 

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
