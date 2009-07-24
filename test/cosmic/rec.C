void rec(const char *filename="raw.root")
{
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  // First version of the reconstruction
  // script for the 2009 data (LHC09b) 
  //
  /////////////////////////////////////////////////////////////////////////////////////////
  gSystem->Load("libRAliEn.so");
  gSystem->Load("libNet.so");
  gSystem->Load("libMonaLisa.so");
  new TMonaLisaWriter(0, "GridAliRoot-rec.C", 0, 0, "global");
  gSystem->Setenv("APMON_INTERVAL", "120");

  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  
  // Reconstruction settings
  AliReconstruction rec;

  // QA options
  rec.SetRunQA("ALL:ALL") ;
  rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  AliGRPRecoParam *grpRecoParam = AliGRPRecoParam::GetLowFluxParam();
  grpRecoParam->SetVertexerTracksConstraintITS(kFALSE);
  grpRecoParam->SetVertexerTracksConstraintTPC(kFALSE);
  grpRecoParam->SetMostProbablePt(3.0);
  rec.SetRecoParam("GRP",grpRecoParam);

  // AliReconstruction settings
  rec.SetWriteESDfriend(kTRUE);
  rec.SetWriteAlignmentData();
  rec.SetInput(filename);
  rec.SetRunReconstruction("ALL");
  rec.SetUseTrackingErrorsForAlignment("ITS");

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  AliLog::Flush();
  rec.Run();

}
