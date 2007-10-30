void MakeTPCAltroMapping(){
  // Create TObjArray of TPC altro mapping objects and
  // store it in the CDB
  //
  const char* macroname = "MakeTPCAltroMapping.C";

  TObjArray mappingsArray(6);
  
  TString path = gSystem->Getenv("ALICE_ROOT");
  path += "/TPC/mapping/Patch";
  TString path2;
  for(Int_t i = 0; i < 6; i++) {
    path2 = path;
    path2 += i;
    path2 += ".data";
    AliAltroMapping *mapping = new AliTPCAltroMapping(path2.Data());
    mappingsArray.Add(mapping);
  }

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "TPCAltroMapping.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving altro mapping objects to the file %s", filename);
    f.cd();
    f.WriteObject(&mappingsArray,"TPCAtroMappings","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
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
    md->SetResponsible("Marian Ivanov");
    md->SetComment("Default ALTRO mapping for TPC");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("TPC/Calib/Mapping",0,AliCDBRunRange::Infinity());
    storage->Put(&mappingsArray,id,md);

    delete md;
  }

  mappingsArray.Delete();

}

