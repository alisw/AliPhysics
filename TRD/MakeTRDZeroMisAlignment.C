void MakeTRDZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for TRD
  //
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",1000);
  TClonesArray &alobj = *array;
   
  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

  AliAlignObjAngles a;

  Double_t dx=0.,dy=0.,dz=0.,rx=0.,ry=0.,rz=0.;

  Int_t j=0;
  UShort_t volid;
  const char *symname;

  // create the chambers' alignment objects
  for (Int_t iLayer = AliAlignObj::kTRD1; iLayer <= AliAlignObj::kTRD6; iLayer++) {
    for (Int_t iModule = 0; iModule < AliAlignObj::LayerSize(iLayer); iModule++) {
      volid = AliAlignObj::LayerToVolUID(iLayer,iModule);
      symname = AliAlignObj::SymName(volid);
      new(alobj[j++]) AliAlignObjAngles(symname,volid,dx,dy,dz,rx,ry,rz,kTRUE);
    }
  }

  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("TRDzeroMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"TRDAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Dariusz Miskowiec");
    md->SetComment("Zero misalignment for TRD");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("TRD/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}


