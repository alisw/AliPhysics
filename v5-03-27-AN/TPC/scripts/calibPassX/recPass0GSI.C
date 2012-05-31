//
//   rec.C to be used for pass0
//   

void rec(const char *filename="raw.root",Int_t nevents=-1)
{
  // Load some system libs for Grid and monitoring
  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  //man->SetDefaultStorage("raw://");
  man->SetDefaultStorage("local:///lustre/alice/alien/alice/data/2010/OCDB/");
  // Reconstruction settings
  AliReconstruction rec;

  // Set protection against too many events in a chunk (should not happen)
  if (nevents>0) rec.SetEventRange(0,nevents);

  // Switch off HLT until the problem with schema evolution resolved
  //rec.SetRunReconstruction("ALL-HLT");
  //
  // QA options
  //
  AliQAManager *qam = AliQAManager::QAManager(AliQAv1::kRECMODE) ;
  rec.SetRunQA(":");
  rec.SetRunGlobalQA(kFALSE);

  // AliReconstruction settings
  rec.SetWriteESDfriend(kTRUE);
  rec.SetWriteAlignmentData();
  rec.SetInput(filename);
  rec.SetUseTrackingErrorsForAlignment("ITS");
  rec.SetRunReconstruction("ITS TPC TRD TOF");
  rec.SetFillESD("ITS TPC TRD TOF");

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  AliLog::Flush();
  rec.Run();

}


