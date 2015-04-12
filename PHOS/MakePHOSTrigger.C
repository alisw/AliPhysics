void MakePHOSTrigger()
{

  const char* macroname = "MakePHOSTrigger.C";

  gSystem->SetIncludePath("-I$ALICE_ROOT/PHOS/PHOSbase");
  gROOT->SetMacroPath("$ALICE_ROOT/PHOS/macros/Trigger/OCDB");
  gROOT->LoadMacro("readTRUPedestals.C+g");

  AliPHOSTriggerParameters* parameters = new AliPHOSTriggerParameters("phosTParametersDat");
  readTRUPedestals(parameters);

  TString Storage = gSystem->Getenv("STORAGE");
  if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
    Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
    return;
  }
  Info(macroname,"Saving altro mapping objects in CDB storage %s",
       Storage.Data());
  AliCDBManager* cdb = AliCDBManager::Instance();
  AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
  if(!storage){
    Error(macroname,"Unable to open storage %s\n",Storage.Data());
    return;
  }
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Yuri Kharlov");
  md->SetComment("PHOS trigger parameters for Run1");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  AliCDBId id("PHOS/Trigger/Parameters",0,AliCDBRunRange::Infinity());
  storage->Put(parameters,id,md);

  drawTRUPedestals(parameters);
}
  
