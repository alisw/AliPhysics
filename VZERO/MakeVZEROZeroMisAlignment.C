void MakeVZEROZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for VZERO
  //
  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

  TClonesArray *array = new TClonesArray("AliAlignObjAngles",10);
  TClonesArray &alobj = *array;

  AliAlignObjAngles a;

  Double_t dx, dy, dz, dpsi, dtheta, dphi;

  Int_t iIndex=0;
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iIndex);

  TString V0right("VZERO/V0C");
  new(alobj[0]) AliAlignObjAngles(V0right.Data(), volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  TString V0left("VZERO/V0A");
  new(alobj[1]) AliAlignObjAngles(V0left.Data(), volid, dx, dy, dz, dpsi, dtheta, dphi,kTRUE);

  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("V0zeroMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"V0ZeroAlObjs ","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Brigitte Cheynis");
    md->SetComment("Alignment objects for V0 zero-misalignment");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("VZERO/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}

