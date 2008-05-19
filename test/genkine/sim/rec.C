void rec() {
  AliReconstruction reco;

  //  reco.SetRunReconstruction("ITS TPC TRD TOF HMPID FMD PMD VZERO START MUON ZDC");

// **** The field map settings must be the same as in Config.C !
  AliMagFMaps *field=new AliMagFMaps("Maps","Maps",2,1.,10.,AliMagFMaps::k5kG);
  Bool_t uniform=kFALSE;
  AliTracker::SetFieldMap(field,uniform);

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
