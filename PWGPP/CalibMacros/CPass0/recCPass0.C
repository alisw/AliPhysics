/*
   rec.C to be used for pass0
   - reconstruction of raw data
   - QA information switched off
   - store all friends
   - default OCDB storage set to "raw://"

   Example:
   aliroot -b -q 'recCPass0.C("raw.root",100)'
*/

void recCPass0(const char *filename="raw.root",Int_t nevents=-1, const char *ocdb="raw://")
{
  // Load some system libs for Grid and monitoring
  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdb);
  // Reconstruction settings
  AliReconstruction rec;
  // Upload CDB entries from the snapshot (local root file) if snapshot exist
  if (gSystem->AccessPathName("OCDB.root", kFileExists)==0) {        
    rec.SetFromCDBSnapshot("OCDB.root");
  }
  // All friends
  rec.SetFractionFriends(1.0);

 // AliReconstruction settings - hardwired MB trigger for calibration

  TString newfilename = filename;
  newfilename += "?Trigger=kCalibBarrel";
  rec.SetInput(newfilename.Data());

  // Set protection against too many events in a chunk (should not happen)
  if (nevents>0) rec.SetEventRange(0,nevents);

  // Remove recpoints after each event
  rec.SetDeleteRecPoints("TPC TRD ITS");

  // Switch off the V0 finder - saves time!
  //  rec.SetRunMultFinder(kFALSE);
  rec.SetRunV0Finder(kFALSE); 

  //
  // QA options - all QA is off
  //
  rec.SetRunQA(":");
  rec.SetRunGlobalQA(kFALSE);

  // AliReconstruction settings
  rec.SetWriteESDfriend(kTRUE);
  rec.SetWriteAlignmentData();
  rec.SetUseTrackingErrorsForAlignment("ITS");
  rec.SetRunReconstruction("ALL");
  rec.SetFillESD("ALL");
  rec.SetCleanESD(kFALSE);

  // Specific reco params for ZDC (why isn't this automatic?)
//  rec.SetRecoParam("ZDC",AliZDCRecoParamPbPb::GetHighFluxParam(2760));

  //Ignore SetStopOnError
  rec.SetStopOnError(kFALSE);

  AliLog::Flush();
  rec.Run();
}

