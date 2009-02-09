void MakeACORDEZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for ACORDE
  //
  const char* macroname = "MakeACORDEZeroMisAlignment.C";
  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB/ACORDE/Align/Data");
  cdb->SetRun(0);
  
  AliCDBStorage* storage;
  
  //load geom from local file till ACORDE is not switched on by default in standard config-files
  if( TString(gSystem->Getenv("TOCDB")) == TString("kTRUE") ){
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
    AliGeomManager::LoadGeometry(); //load geom from default CDB storage
    
  }
  //  AliGeomManager::LoadGeometry("geometry.root");  

  TClonesArray *array = new TClonesArray("AliAlignObjParams",64);
  TClonesArray &alobj = *array;
  
  TRandom *rnd = new TRandom(4321);
  Int_t j = 0;
  Double_t dx=0, dy=0, dz=0, dpsi=0, dtheta=0, dphi=0;

  // RS = local
  // sigma translation 
  // sigma rotation 
  
  TString symname;
  TString basename = "ACORDE/Array";
  Int_t iIndex=0;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

  for (Int_t imod=1; imod<61; imod++){
    symname = basename;
    symname += imod;    
    new(alobj[j++]) AliAlignObjParams(symname, volid, dx, dy, dz,dpsi, dtheta, dphi, kFALSE);
  }

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "ACORDEZeroMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"ACORDEAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("E. Cuautle & M. Rodriguez");
    md->SetComment("Zero misalignment for ACORDE");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("ACORDE/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}

