/*
   rec.C to be used for pass0
   - reconstruction of raw data  
   - QA information switched off
   - store all friends
   - default OCDB storage set to "raw://"

   Example:
   aliroot -b -q 'recPass0.C("raw.root",100)'
*/

void recPass0(const char *filename="raw.root",Int_t nevents=-1, const char *ocdb="raw://")
{
  // Load some system libs for Grid and monitoring
  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdb);
 
  // Reconstruction settings
  AliReconstruction rec;

  // All friends
  rec.SetFractionFriends(1.0);

  // Set protection against too many events in a chunk (should not happen)
  if (nevents>0) rec.SetEventRange(0,nevents);

  //
  // QA options - all QA is off
  //
  rec.SetRunQA(":");
  rec.SetRunGlobalQA(kFALSE);

  // AliReconstruction settings
  rec.SetWriteESDfriend(kTRUE);
  rec.SetWriteAlignmentData();
  rec.SetInput(filename);
  rec.SetUseTrackingErrorsForAlignment("ITS");
  rec.SetRunReconstruction("ALL");
  rec.SetFillESD("ALL");
  rec.SetCleanESD(kFALSE);

  //Ignore SetStopOnError
  rec.SetStopOnError(kFALSE);

  AliLog::Flush();
  rec.Run();
}


