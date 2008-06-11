void rec() {
  AliReconstruction reco;

  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetLowFluxParam();
  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  AliTPCReconstructor::SetStreamLevel(1);
  //   reco.SetInput("raw.root");
  reco.SetRunReconstruction("ITS TPC TRD TOF HMPID PHOS EMCAL MUON VZERO T0 FMD PMD ZDC");
  reco.SetRunQA("");
  reco.SetRunGlobalQA(kFALSE);
  reco.SetMeanVertexConstraint(kFALSE);

// **** The field map settings must be the same as in Config.C !
  AliMagWrapCheb* field = 0x0;
  //field = new AliMagWrapCheb("Maps","Maps", 2, 0., 10., AliMagWrapCheb::k2kG);
  //Bool_t uniform=kTRUE;
  //AliTracker::SetFieldMap(field,uniform); // tracking with the uniform field
  field = new AliMagWrapCheb("Maps","Maps", 2, 1., 10., AliMagWrapCheb::k5kG);
  Bool_t uniform=kFALSE;
  AliTracker::SetFieldMap(field,uniform);  // tracking with the real map

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
