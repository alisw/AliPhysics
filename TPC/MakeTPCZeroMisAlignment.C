void MakeTPCZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for TPC
  //
  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

  TClonesArray *array = new TClonesArray("AliAlignObjAngles",100);
  TClonesArray &alobj = *array;
  
  AliAlignObjAngles o;
  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;
  Int_t j = 0;

  // RS = local
  for (Int_t iLayer = AliAlignObj::kTPC1; iLayer <= AliAlignObj::kTPC2; iLayer++) {
    for (Int_t iModule = 0; iModule < AliAlignObj::LayerSize(iLayer); iModule++) {

      UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iModule);
      const char *symname = AliAlignObj::SymName(volid);
      new(alobj[j]) AliAlignObjAngles(symname, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
      j++;
    }
  }


  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("TPCzeroMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"TPCAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Marian Ivanov");
    md->SetComment("Zero misalignment for TPC");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("TPC/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}

