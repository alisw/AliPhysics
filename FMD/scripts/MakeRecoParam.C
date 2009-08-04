void
MakeRecoParam()
{
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);

  AliFMDRecoParam param;
  AliCDBId        id("FMD/Calib/RecoParam",0,999999);
  AliCDBMetaData  meta;
  meta = new AliCDBMetaData;					    
  meta.SetResponsible(gSystem->GetUserInfo()->fRealName.Data());	
  meta.SetAliRootVersion(gROOT->GetVersion()); 
  meta.SetBeamPeriod(1);						
  meta.SetComment("Dummy data for testing"); 
  meta.SetProperty("key1", &param);
  
  cdb->Put(&param, id, &meta);
}
  
 
