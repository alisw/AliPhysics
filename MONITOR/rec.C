void rec(const char *filename="raw.root")
{
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  // Script for the online reconstruction/visualization
  //
  /////////////////////////////////////////////////////////////////////////////////////////

  // Setting CDB
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("local:///local/cdb");
  man->SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s/..",gSystem->pwd()));
  man->SetSpecificStorage("GRP/CTP/Config",
			  Form("local://%s/..",gSystem->pwd()));
  
  // Reconstruction settings
  AliReconstruction rec;

  // QA options
  rec.SetRunQA(":") ;
  rec.SetRunGlobalQA(kFALSE);
  rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

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
