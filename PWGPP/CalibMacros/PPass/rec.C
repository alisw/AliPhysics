void rec(const char *filename="raw.root", int nevents=-1, const char* ocdb="/cvmfs/alice-ocdb.cern.ch/calibration/data", TString additionalRecOptions="")
{
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  // Reconstruction script for 2010 RAW data
  //
  /////////////////////////////////////////////////////////////////////////////////////////

  
  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdb);

  AliReconstruction rec;

  // Set reconstruction flags (skip detectors here if neded with -<detector name>

  rec.SetRunReconstruction("ALL");

  // QA options
  //  rec.SetRunQA("Global:ESDs") ;
  //  rec.SetRunQA(":") ;
  //  rec.SetRunQA("ALL:ALL") ;
  rec.SetRunQA("Global MUON:ALL") ;

  rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  // AliReconstruction settings
  rec.SetWriteESDfriend(kTRUE);
  rec.SetWriteAlignmentData();
  rec.SetInput(filename);
  rec.SetUseTrackingErrorsForAlignment("ITS");

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

  // Set protection against too many events in a chunk (should not happen)
  if (nevents>0) rec.SetEventRange(0,nevents);

  // Upload CDB entries from the snapshot (local root file) if snapshot exist
  if (gSystem->AccessPathName("OCDB.root", kFileExists)==0) {
    rec.SetCDBSnapshotMode("OCDB.root");
  }

  // Specific AD storage, see https://alice.its.cern.ch/jira/browse/ALIROOT-6056
  //  rec.SetSpecificStorage("AD/Calib/TimeSlewing", "alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal");

  //load any specific local storages
  if (gSystem->AccessPathName("localOCDBaccessConfig.C", kFileExists)==0) {        
    gROOT->LoadMacro("localOCDBaccessConfig.C");
    localOCDBaccessConfig();
  }

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  //Ignore SetStopOnError
  rec.SetStopOnError(kFALSE);

  // Delete recpoints
  rec.SetDeleteRecPoints("TPC TRD");

  // Set 100% of friends
  // rec.SetFractionFriends(2.0);

  AliLog::Flush();
  rec.Run();

}
