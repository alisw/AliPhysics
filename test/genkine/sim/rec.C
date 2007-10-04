void rec() {
  AliReconstruction reco;

  //  reco.SetRunReconstruction("ITS TPC TRD TOF HMPID FMD PMD VZERO START MUON ZDC");

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
