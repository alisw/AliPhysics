void MakePHOSCpvBadChannels()
{

  const char* macroname = "MakePHOSCpvBadChannels.C";

  AliPHOSCpvBadChannelsMap* badmap = new AliPHOSCpvBadChannelsMap();

  TString Storage = gSystem->Getenv("STORAGE");
  if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
    Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
    return;
  }
  Info(macroname,"Saving CPV bad channel map in CDB storage %s",
       Storage.Data());
  AliCDBManager* cdb = AliCDBManager::Instance();
  AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
  if(!storage){
    Error(macroname,"Unable to open storage %s\n",Storage.Data());
    return;
  }
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Yuri Kharlov");
  md->SetComment("CPV default bad channel map");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  AliCDBId id("PHOS/Calib/CpvBadChannels",0,AliCDBRunRange::Infinity());
  storage->Put(badmap,id,md);
}
  
