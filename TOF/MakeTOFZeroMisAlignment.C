void MakeTOFZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for TOF
  //
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",2000);
  TClonesArray &alobj = *array;
   
  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

  AliAlignObjAngles a;

  AliAlignObj::ELayerID idTOF = AliAlignObj::kTOF;
  Int_t i;
  Int_t j=0;
  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;

  for(i=0; i<AliAlignObj::LayerSize(idTOF); i++) {
    new(alobj[j++]) AliAlignObjAngles(AliAlignObj::SymName(idTOF,i), AliAlignObj::LayerToVolUID(idTOF,i), dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }

  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("TOFzeroMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"TOFAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Silvia Arcelli");
    md->SetComment("Zero misalignment for TOF");
    md->SetAliRootVersion("HEAD");
    AliCDBId id("TOF/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}


