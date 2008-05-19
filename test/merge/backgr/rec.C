void rec() {

  AliReconstruction reco;
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  AliTPCReconstructor::SetStreamLevel(1);
  AliTPCReconstructor::SetRecoParam(AliTPCRecoParam::GetLowFluxParam());
  //  reco.SetInput("raw.root");

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
