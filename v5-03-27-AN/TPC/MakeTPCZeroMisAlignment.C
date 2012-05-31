void MakeTPCZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for TPC
  //
  const char* macroname = "MakeTPCZeroMisAlignment.C";

  TClonesArray *array = new TClonesArray("AliAlignObjParams",100);
  TClonesArray &alobj = *array;
  
  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  
  AliCDBStorage* storage = NULL;

  if(TString(gSystem->Getenv("TOCDB")) == TString("kTRUE")){
    TString Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    storage = cdb->GetStorage(Storage.Data());
    if(!storage){
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }
    AliCDBPath path("GRP","Geometry","Data");
    AliCDBEntry *entry = storage->Get(path.GetPath(),cdb->GetRun());
    if(!entry) Fatal(macroname,"Could not get the specified CDB entry!");
    entry->SetOwner(0);
    TGeoManager* geom = (TGeoManager*) entry->GetObject();
    AliGeomManager::SetGeometry(geom);
  }else{
    AliGeomManager::LoadGeometry("geometry.root"); //load geom from default CDB storage
  }    
  
  Double_t dx=-0.159, dy=-0.05, dz=0.034, dpsi=-0.00183, dtheta=0.01835, dphi=0.02865;
  Int_t j = 0;

  new(alobj[j++]) AliAlignObjParams("ALIC_1/TPC_M_1", 0, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  
  // RS = local
  dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;
  for (Int_t iLayer = AliGeomManager::kTPC1; iLayer <= AliGeomManager::kTPC2; iLayer++) {
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer); iModule++) {

      UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iModule);
      const char *symname = AliGeomManager::SymName(volid);
      new(alobj[j++]) AliAlignObjParams(symname, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
    }
  }

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "TPCzeroMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"TPCAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    Info(macroname,"Saving alignment objects in CDB storage %s",
	 Storage.Data());
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Marian Ivanov");
    md->SetComment("Zero misalignment for TPC");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("TPC/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}

