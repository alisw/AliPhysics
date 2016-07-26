/*
   rec.C to be used for pass1
   - reconstruction of raw data
   - QA information switched off
   - store all friends
   - default OCDB storage set to "raw://"

   Example:
   aliroot -b -q 'recCPass1.C("raw.root",100)'
*/

void recCPass1(const char *filename="raw.root", Int_t nevents=-1, const char *ocdb="raw://", const char* options="?Trigger=kCalibBarrelMB", TString additionalRecOptions="")
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
  //
  // do we extract the TPC recpoints in advance
  Bool_t noTPCLocalRec = gSystem->Getenv("preclusterizeTPC")!=NULL;
  if (noTPCLocalRec) printf("TPC local reconstruction assumed to be already done\n");
  //
  if (gSystem->Getenv("disableOuter")!=NULL){
    TString disOuter = gSystem->Getenv("disableOuter");
    TString disOuterLoc = disOuter;
    if (noTPCLocalRec) {
      disOuterLoc.ReplaceAll("TPC","");
      disOuterLoc.ReplaceAll("HLT","");
    }
    rec.SetRunReconstruction(disOuter.Data());
    rec.SetRunLocalReconstruction(disOuterLoc.Data());
  } 
  else if (noTPCLocalRec) {
    rec.SetRunReconstruction("ALL -HLT");
    rec.SetRunLocalReconstruction("ALL -TPC -HLT");
  }
  else {
    rec.SetRunLocalReconstruction("ALL");
  }

  // Set additional reconstruction options.
  // They are in the form "Detector:value;Detector2:value" in a single string.
  // For instance: additionalRecOptions="TPC:useHLTorRAW"
  {
    TIter nexttok( additionalRecOptions.Tokenize(";") );
    while (( os = (TObjString *)nexttok() )) {
      TString detOpt = os->String();
      Ssiz_t idx = detOpt.Index(":");
      if (idx < 0) continue;
      TString detOptKey = detOpt(0,idx);
      TString detOptVal = detOpt(idx+1,999);
      if (detOptKey.IsNull() || detOptVal.IsNull()) continue;
      printf("Setting additional reconstruction option: %s --> %s\n", detOptKey.Data(),
                                                                      detOptVal.Data());
      rec.SetOption(detOptKey.Data(), detOptVal.Data());
    }
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
  TString delRecPoints="TPC TRD ITS";
  if (noTPCLocalRec) delRecPoints.ReplaceAll("TPC","");
  rec.SetDeleteRecPoints(delRecPoints.Data()); 
  //

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

