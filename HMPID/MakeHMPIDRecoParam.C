void MakeHMPIDRecoParam(AliRecoParam::EventSpecie_t default=AliRecoParam::kLowMult)
{
  //
  // Create HMPID Reco Parameters in OCDB
  //  
  AliCDBManager* man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet()) man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  TObject *recoParamArray = new TObject();
  
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Domenico DiBari");
  md->SetComment("Reconstruction parameters of HMPID, v1");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("HMPID/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  man->GetDefaultStorage()->Put(recoParamArray,id, md);
  
  
  
  return;
}
