void raw2clust(const char *filename="raw.root", Int_t nevents=-1,const char *ocdb="raw://", const char* options="")
{
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  // Reconstruction script for 2010 RAW data
  //
  /////////////////////////////////////////////////////////////////////////////////////////

  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();

  AliReconstruction rec;

  man->SetDefaultStorage(ocdb);
  // Upload CDB entries from the snapshot (local root file) if snapshot exist
  if (gSystem->AccessPathName("OCDB.root", kFileExists)==0) {        
    rec.SetCDBSnapshotMode("OCDB.root");
  }

  if (gSystem->AccessPathName("localOCDBaccessConfig.C", kFileExists)==0) {        
    gROOT->LoadMacro("localOCDBaccessConfig.C");
    localOCDBaccessConfig();
  }

  rec.SetRunReconstruction("");
  rec.SetRunLocalReconstruction("TPC HLT");

  // use compact clusters
  AliTPCReconstructor::SetCompactClusters(kTRUE); 

  // QA options
  rec.SetRunQA(":") ;
  rec.SetRunGlobalQA(kFALSE);

  // AliReconstruction settings
  TString newfilename = filename;
  newfilename += options;
  rec.SetInput(newfilename.Data());

  //Ignore SetStopOnError
  rec.SetStopOnError(kFALSE);

  // Set protection against too many events in a chunk (should not happen)
  if (nevents>0) rec.SetEventRange(0,nevents);

  AliLog::Flush();
  rec.Run();

}
