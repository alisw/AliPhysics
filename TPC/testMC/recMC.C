/// \file recMC.C

void recMC(const char* tpcDBpath="local://$ALICE_ROOT", Int_t nevents){
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetLowFluxParam();  
  AliTPCReconstructor::SetStreamLevel(10);
  //

  AliReconstruction rec;
  rec.SetSpecificStorage("TPC/*/*",tpcDBpath);
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
  //disable HLT and QA
  //
  rec.SetRunHLTTracking(kFALSE);
  rec.SetUseHLTData(" ");
  rec.SetRunQA(kFALSE);
  //
  rec.SetWriteESDfriend(kTRUE);
  rec.SetWriteAlignmentData(kTRUE);
  rec.SetUniformFieldTracking(kFALSE);
  rec.SetWriteAlignmentData(kTRUE);
  //
  rec.SetEventRange(0,nevents);
  //
  rec.Run();
}
