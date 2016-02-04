void rec(const char *filename="raw.root")
{
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  // Reconstruction script for 2010 RAW data
  //
  /////////////////////////////////////////////////////////////////////////////////////////

  
  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("alien://folder=/alice/data/2015/OCDB");

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
  rec->SetDeleteRecPoints("TPC TRD");

  // Set 100% of friends
  // rec.SetFractionFriends(2.0);

  AliLog::Flush();
  rec.Run();

}
