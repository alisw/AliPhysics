void MakePHOSAltroMapping(){
  // Create TObjArray of PHOS altro mapping objects.
  // Set the environment variables in order to store it in the OCDB:
  // export TOCDB=kTRUE
  // export STORAGE=local://$ALICE_ROOT/OCDB
  // Then the newly created root file $ALICE_ROOT/OCDB/PHOS/Calib/Mapping
  // should be committed to SVN and submitted to AliEn production manager
  //
  // Yuri Kharlov. 10 September 2009
  // $Id$

  const char* macroname = "MakePHOSAltroMapping.C";

  TObjArray mappingsArray(20);
  
  nModules = 5;
  nRCU     = 4;
  TString path = gSystem->Getenv("ALICE_ROOT");
  path += "/PHOS/mapping/";
  TString path1, path2;
  for(Int_t m = 0; m < nModules; m++) {
    path1 = path;
    path1 += "Mod";
    path1 += m;
    path1 += "RCU";
    for(Int_t i = 0; i < nRCU; i++) {
      path2 = path1;
      path2 += i;
      path2 += ".data";
      Info(macroname,"Mapping file: %s",path2.Data());
      AliAltroMapping *mapping = new AliCaloAltroMapping(path2.Data());
      mappingsArray.Add(mapping);
    }
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
    md->SetComment("Default ALTRO mapping for PHOS: 20 mapping objects, one per modules per RCU");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("PHOS/Calib/Mapping",0,AliCDBRunRange::Infinity());
    storage->Put(&mappingsArray,id,md);

    delete md;
  }

  mappingsArray.Delete();

}

