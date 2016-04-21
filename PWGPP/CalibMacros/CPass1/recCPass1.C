/*
   rec.C to be used for pass0
   - reconstruction of raw data
   - QA information switched off
   - store all friends
   - default OCDB storage set to "raw://"

   Example:
   aliroot -b -q 'recCPass1.C("raw.root",100)'
*/

void recCPass1(const char *filename="raw.root",Int_t nevents=-1, const char *ocdb="raw://", const char* options="?Trigger=kCalibBarrel")
{
  if (gSystem->Getenv("ALIROOT_FORCE_COREDUMP"))
  {
    printf("ALIROOT_FORCE_COREDUMP set\n");
    gSystem->ResetSignal(kSigFloatingException);
    gSystem->ResetSignal(kSigSegmentationViolation);
  }

  // removing apparently pile-up clusters to speadup reconstruction
  const double kZOutSectorCut = 3.; // cut on clusters on wrong side of CE (added to extendedRoadZ)
  AliTPCReconstructor::SetZOutSectorCut(kZOutSectorCut);

  // Load some system libs for Grid and monitoring
  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdb);
  
  // Reconstruction settings
  AliReconstruction rec;
  if (gSystem->Getenv("disableOuter")){
    rec.SetRunLocalReconstruction("ITS TPC TRD TOF T0");
    rec.SetRunReconstruction("ITS TPC TRD TOF T0");
    rec.SetRunTracking("ITS TPC TRD TOF T0");
  } else {
    rec.SetRunReconstruction("ALL");
  }
  // Upload CDB entries from the snapshot (local root file) if snapshot exist
  if (gSystem->AccessPathName("OCDB.root", kFileExists)==0) {        
    rec.SetCDBSnapshotMode("OCDB.root");
  }

  if (gSystem->AccessPathName("localOCDBaccessConfig.C", kFileExists)==0) {        
    gROOT->LoadMacro("localOCDBaccessConfig.C");
    localOCDBaccessConfig();
  }

  // All friends
  rec.SetFractionFriends(2.0);

 // AliReconstruction settings - hardwired MB trigger for calibration

  TString newfilename = filename;
  newfilename += options;
  rec.SetInput(newfilename.Data());

  // Set protection against too many events in a chunk (should not happen)
  if (nevents>0) rec.SetEventRange(0,nevents);

  // Remove recpoints after each event
  rec.SetDeleteRecPoints("TPC TRD ITS");


  // Switch off the V0 finder - saves time!
  rec.SetRunMultFinder(kTRUE);
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
  rec.SetCleanESD(kFALSE);

  //Ignore SetStopOnError
  rec.SetStopOnError(kFALSE);

  AliLog::Flush();
  rec.Run();
}

