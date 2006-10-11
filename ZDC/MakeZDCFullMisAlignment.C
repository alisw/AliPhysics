void MakeZDCFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for ZDC
  //
  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

  TClonesArray *array = new TClonesArray("AliAlignObjAngles",10);
  TClonesArray &alobj = *array;

  AliAlignObjAngles a;

  Double_t dx=0., dy=2., dz=0.;
  Double_t dpsi=0., dtheta=0., dphi=0.;

  const char *ZDC_N="ALIC_1/ZDC_1/ZNEU_1";
  const char *ZDC_P="ALIC_1/ZDC_1/ZPRO_1";

  Int_t iIndex=0; //let all modules have index=0 in a layer with no LUT
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iIndex);

  new(alobj[0]) AliAlignObjAngles(ZDC_N, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[1]) AliAlignObjAngles(ZDC_P, volid, dx, dy, dz, dpsi, dtheta, dphi,kTRUE);

  if(!gSystem->Getenv("$TOCDB")){
    // save in file
    TFile f("ZDCfullMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"ZDCFullObjs ","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Chiara Oppedisano");
    md->SetComment("Alignment objects for ZDC full misalignment");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("ZDC/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}

