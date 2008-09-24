void MakeEMCALAltroMapping(){
  // Create TObjArray of EMCAL altro mapping objects and
  // store it in the CDB
  //
  const char* macroname = "MakeEMCALAltroMapping.C";

  TObjArray mappingsArray(6);
  
  TString path = gSystem->Getenv("ALICE_ROOT");
  path += "/EMCAL/mapping/RCU";
  TString path2;
  TString side[] = {"A","C"};//+ and - pseudarapidity supermodules
  for(Int_t j = 0; j < 2; j++){
    for(Int_t i = 0; i < 2; i++) {
      path2 = path;
      path2 += i;
      path2 +=side[j]; 
      path2 += ".data";
      printf("File: %s\n",path2.Data());
      AliAltroMapping *mapping = new AliCaloAltroMapping(path2.Data());
      mappingsArray.Add(mapping);
    }
  }
  
  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "EMCALAltroMapping.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving altro mapping objects to the file %s", filename);
    f.cd();
    f.WriteObject(&mappingsArray,"EMCALAtroMappings","kSingleKey");
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
    md->SetResponsible("Jennifer Klay");
    md->SetComment("Default ALTRO mapping for EMCAL");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("EMCAL/Calib/Mapping",0,AliCDBRunRange::Infinity());
    storage->Put(&mappingsArray,id,md);

    delete md;
  }

  mappingsArray.Delete();

}

