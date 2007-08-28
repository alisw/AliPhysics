void MakeTPCZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for TPC
  //
  if(!AliGeomManager::GetGeometry()){
    if(!(AliCDBManager::Instance())->IsDefaultStorageSet())
      AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
      AliCDBManager::Instance()->SetRun(0);
    AliGeomManager::LoadGeometry();
  }
  // needed for the constructors with local coordinates not to fail

  TClonesArray *array = new TClonesArray("AliAlignObjParams",100);
  TClonesArray &alobj = *array;
  
  AliAlignObjParams o;
  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;
  Int_t j = 0;

  // RS = local
  for (Int_t iLayer = AliGeomManager::kTPC1; iLayer <= AliGeomManager::kTPC2; iLayer++) {
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer); iModule++) {

      UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iModule);
      const char *symname = AliGeomManager::SymName(volid);
      new(alobj[j]) AliAlignObjParams(symname, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
      j++;
    }
  }


  const char* macroname = "MakeTPCZeroMisAlignment.C";
  if( gSystem->Getenv("TOCDB") != TString("kTRUE") ){
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
    TString Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    Info(macroname,"Saving alignment objects in CDB storage %s",
	 Storage.Data());
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
    if(!storage){
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Marian Ivanov");
    md->SetComment("Zero misalignment for TPC");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("TPC/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}

