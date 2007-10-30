void MakePHOSAltroMapping(){
  // Create TObjArray of PHOS altro mapping objects and
  // store it in the CDB
  //
  const char* macroname = "MakePHOSAltroMapping.C";

  TObjArray mappingsArray(6);
  
  TString path = gSystem->Getenv("ALICE_ROOT");
  path += "/PHOS/mapping/RCU";
  TString path2;
  for(Int_t i = 0; i < 4; i++) {
    path2 = path;
    path2 += i;
    path2 += ".data";
    AliAltroMapping *mapping = new AliCaloAltroMapping(path2.Data());
    mappingsArray.Add(mapping);
  }

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "PHOSAltroMapping.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving altro mapping objects to the file %s", filename);
    f.cd();
    f.WriteObject(&mappingsArray,"PHOSAtroMappings","kSingleKey");
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
    md->SetResponsible("Yuri Kharlov");
    md->SetComment("Default ALTRO mapping for PHOS");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("PHOS/Calib/Mapping",0,AliCDBRunRange::Infinity());
    storage->Put(&mappingsArray,id,md);

    delete md;
  }

  mappingsArray.Delete();

}

