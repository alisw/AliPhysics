void MakeITSZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for ITS
  //
  TClonesArray *array = new TClonesArray("AliAlignObjParams",4000);
  TClonesArray &alobj = *array;
  const char* macroname = "MakeITSZeroMisAlignment.C";

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

  Int_t j = 0;

  //=****************************************
  // overall ITS misalignment according to survey as reported by Werner Riegler (18/07/2008) 
  //=****************************************
  Float_t its_dx     = -0.12;
  Float_t its_dy     = -0.07;
  Float_t its_dz     = 0.29;
  Float_t its_dpsi   = 0.;  
  Float_t its_dtheta = 0.03;
  Float_t its_dphi   = 0.04;

  new(alobj[j++]) AliAlignObjParams("ITS", 0, its_dx, its_dy, its_dz, its_dpsi, its_dtheta, its_dphi, kTRUE);

  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0., globalZ=0.;

  for ( Int_t iLayer = AliGeomManager::kSPD1; iLayer <= AliGeomManager::kSSD2; iLayer++) {
    
    printf("%i modules in layer %i\n", AliGeomManager::LayerSize(iLayer), iLayer);
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer); iModule++) {

      UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iModule);
      const char *symname = AliGeomManager::SymName(volid);

      new(alobj[j++]) AliAlignObjParams(symname, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);

    }
  }

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "ITSzeroMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"ITSAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    Info(macroname,"Saving alignment objects in CDB storage %s", Storage.Data());
    AliCDBMetaData *md= new AliCDBMetaData();
    md->SetResponsible("Ludovic Gaudichet");
    md->SetComment("Alignment objects with zero ITS misalignment");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("ITS/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id, md);
  }

  array->Delete();

}


