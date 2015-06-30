void MakeMFTRecoParam(TString Storage = "alien://folder=/alice/cern.ch/user/a/auras/OCDB/") {
  
  const char* macroname = "MakeMFTRecoParam.C";

  TObjArray *obj = new TObjArray();
  AliMFTRecoParam *param = new AliMFTRecoParam();
  obj -> AddAt(param, 0);
  
  // save in CDB storage
  if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
    Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
    return;
  }
  Info(macroname,"Saving Reconstruction Parameters objects in CDB storage %s", Storage.Data());
  AliCDBManager* cdb = AliCDBManager::Instance();
  AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
  if(!storage){
    Error(macroname,"Unable to open storage %s\n",Storage.Data());
    return;
  }
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Antonio Uras");
  md->SetComment("MFT Calibration Data");
  md->SetAliRootVersion(gROOT->GetVersion());
  AliCDBId id("MFT/Calib/Data",0,AliCDBRunRange::Infinity());
  //  AliCDBId id("MFT/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  storage->Put(obj,id,md);

}
