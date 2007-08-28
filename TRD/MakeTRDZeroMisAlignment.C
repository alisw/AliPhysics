void MakeTRDZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for TRD
  //
  TClonesArray *array = new TClonesArray("AliAlignObjParams",1000);
  TClonesArray &alobj = *array;
   
  if(!AliGeomManager::GetGeometry()){
    if(!(AliCDBManager::Instance())->IsDefaultStorageSet())
      AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
      AliCDBManager::Instance()->SetRun(0);
    AliGeomManager::LoadGeometry();
  }
  // needed for the constructors with local coordinates not to fail

  AliAlignObjParams a;

  Double_t dx=0.,dy=0.,dz=0.,rx=0.,ry=0.,rz=0.;

  Int_t j=0;
  UShort_t volid;
  const char *symname;

  // create the chambers' alignment objects
  for (Int_t iLayer = AliGeomManager::kTRD1; iLayer <= AliGeomManager::kTRD6; iLayer++) {
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer); iModule++) {
      volid = AliGeomManager::LayerToVolUID(iLayer,iModule);
      symname = AliGeomManager::SymName(volid);
      new(alobj[j++]) AliAlignObjParams(symname,volid,dx,dy,dz,rx,ry,rz,kTRUE);
    }
  }

  const char* macroname = "MakeTRDZeroMisAlignment.C";
  if( gSystem->Getenv("TOCDB") != TString("kTRUE") ){
    // save on file
    const char* filename = "TRDzeroMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"TRDAlignObjs","kSingleKey");
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
    md->SetResponsible("Dariusz Miskowiec");
    md->SetComment("Zero misalignment for TRD");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("TRD/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}


