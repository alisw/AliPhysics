void rec() {

  AliCDBManager::Instance()->SetRun(0);

  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., AliMagFMaps::k5kG);
  AliTracker::SetFieldMap(field,kTRUE);

  AliReconstruction reco;

  reco.SetUniformFieldTracking(kFALSE);
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetHighFluxParam();
  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  //  AliTPCReconstructor::SetStreamLevel(1);

  reco.SetRunReconstruction("ITS TPC TRD TOF HMPID PHOS EMCAL MUON FMD PMD ZDC T0 VZERO");
  reco.SetInput("raw.root");

  reco.SetNumberOfEventsPerFile(-1); // all events in one single file
  reco.SetRunVertexFinderTracks(kFALSE);

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
