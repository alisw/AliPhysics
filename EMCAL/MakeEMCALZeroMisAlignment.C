void MakeEMCALZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for EMCAL
  //
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",10);
  TClonesArray &alobj = *array;

  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

  AliAlignObjAngles a;

  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;

  const TString fbasepath = "EMCAL/FullSupermodule";
  const TString hbasepath = "EMCAL/HalfSupermodule";
  TString pathstr;

  Int_t iIndex=0; // let all modules have index=0 in a layer with no LUT
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iIndex);

  Int_t i;
  Int_t j=0;

  for(i=0; i<10; i++){
    pathstr=fbasepath;
    pathstr+=(i+1);
    new(alobj[j++]) AliAlignObjAngles(pathstr, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }

  for(i=0; i<2; i++){
    pathstr=hbasepath;
    pathstr+=(i+1);
    new(alobj[j++]) AliAlignObjAngles(pathstr, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }

  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("EMCALzeroMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"EMCALAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Jennifer Clay");
    md->SetComment("Zero misalignment for EMCAL");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("EMCAL/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}

