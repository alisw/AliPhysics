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

  //

  //
  rec.Run();
}
