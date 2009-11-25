void rec(Int_t runNumber, Int_t nev=10000, Int_t firstev=0)
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
  rec.SetRunPlaneEff(kTRUE);

  // AliReconstruction settings
  rec.SetWriteESDfriend(kTRUE);
  rec.SetWriteAlignmentData();
  rec.SetInput(Form("raw://run%d",runNumber));
//  rec.SetRunReconstruction("ALL -HLT -PMD -MUON -PHOS");
  rec.SetRunReconstruction("ALL");
  rec.SetUseTrackingErrorsForAlignment("ITS");

  Int_t lastev=nev+firstev-1;
  rec.SetEventRange(0,lastev);

  // high flux settings (for injection test runs)
  // rec.SetRecoParam("ITS", AliITSRecoParam::GetHighFluxParam());

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  rec.SetOutput(Form("root_archive.zip#AliESDs.root:AliESDs.root,AliESDfriends.root@dataset://run%d",runNumber));

  rec.Run();
}
