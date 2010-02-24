void rec(Int_t runNumber)
{
  gSystem->Load("libRAliEn.so");
  gSystem->Load("libNet.so");

  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");
  
  // Reconstruction settings
  AliReconstruction rec;

  // QA options
  rec.SetRunQA(":") ;
  rec.SetRunGlobalQA(kFALSE);
  rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;
  rec.SetRunPlaneEff(kFALSE);

  // AliReconstruction settings
  rec.SetWriteESDfriend(kTRUE);
  rec.SetWriteAlignmentData();
  rec.SetInput(Form("raw://run%d",runNumber));
  rec.SetRunReconstruction("ALL -HLT");
  rec.SetUseTrackingErrorsForAlignment("ITS");

  rec.SetEventRange(0,10000);

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  rec.SetOutput(Form("root_archive.zip#AliESDs.root:AliESDs.root,AliESDfriends.root@dataset://run%d",runNumber));

  rec.Run();
}
