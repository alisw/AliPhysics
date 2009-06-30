Bool_t MakeCosmicTriggersEntry(const char *fileName)
{
  const char* macroname = "MakeCosmicTriggersEntry.C";

  if (gSystem->AccessPathName(fileName)) {
    Error(macroname,Form("file (%s) not found", fileName));
    return kFALSE;
  }

  ifstream *file = new ifstream(fileName);
  if (!*file) {
    Error(macroname,Form("Error opening file (%s) !",fileName));
    file->close();
    delete file;
    return kFALSE;
  }

  THashTable *table = new THashTable();
  table->SetName("List of defined cosmic triggers");

  TString strLine;
  while (strLine.ReadLine(*file)) {

    if (strLine.BeginsWith("#")) continue;

    strLine.ReplaceAll(" ","");
    strLine.ReplaceAll("\t","");
    if (strLine.IsNull()) continue;

    TObjString *obj = new TObjString(strLine.Data());
    table->Add(obj);
  }

  file->close();
  delete file;


  // save in CDB storage
  TString Storage = gSystem->Getenv("STORAGE");
  if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
    Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
    return;
  }
  Info(macroname,"Saving alignment objects in CDB storage %s",Storage.Data());
  AliCDBManager* cdb = AliCDBManager::Instance();
  AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
  if(!storage){
    Error(macroname,"Unable to open storage %s\n",Storage.Data());
    return;
  }
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Federico Antinori");
  md->SetComment("List of the defined cosmic triggers. It is used in order to steer the reconstruction, namely in the selection of the proper event specie. It is maintained and updated by the trigger coordinator.");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  AliCDBId id("GRP/Calib/CosmicTriggers",0,AliCDBRunRange::Infinity());
  storage->Put(table,id,md);

  table->Delete();
  delete table;

  return kTRUE;
}
