void recMC(){

  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetLowFluxParam();  
  AliTPCReconstructor::SetStreamLevel(10);
  //

  AliReconstruction rec;
  rec.SetRecoParam("TPC",tpcRecoParam);
  //
  rec.SetLoadAlignData("");
  rec.SetWriteESDfriend(kTRUE);
  rec.SetFillTriggerESD(kFALSE);
  rec.SetRunVertexFinder(kTRUE);
  rec.SetWriteAlignmentData(kTRUE);
  rec.SetRunReconstruction("ITS TPC TRD TOF");
  //rec.SetRunReconstruction("ITS TPC");
  rec.SetCleanESD(kFALSE);
  //
  rec.SetWriteESDfriend(kTRUE);
  rec.SetWriteAlignmentData(kTRUE);
  rec.SetRunQA(kFALSE);
  rec.SetUniformFieldTracking(kFALSE);
  rec.SetWriteAlignmentData(kTRUE);
  //
  rec.SetEventRange(0,20);

  AliMagWrapCheb* field = 0x0;
  AliMagF* field = new AliMagWrapCheb("Maps","Maps", 2, 1, 10., AliMagWrapCheb::k5kG,kTRUE,"$(ALICE_ROOT)/data/maps/mfchebKGI_sym.root");
  Bool_t uniform=kFALSE;
  AliTracker::SetFieldMap(field,uniform);  // tracking with the real map
  //

  //
  rec.Run();
}
