void MakeVZEROResMisAlignment(){
  // Create TClonesArray of residual misalignment objects for VZERO
  //
  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

  TClonesArray *array = new TClonesArray("AliAlignObjAngles",10);
  TClonesArray &alobj = *array;

  AliAlignObjAngles a;

  Double_t dx, dy, dz, dpsi, dtheta, dphi;
  TRandom *rnd   = new TRandom(4321);
  Double_t sigmatr = 0.1; // max shift in cm
  Double_t sigmarot = 0.5; // max rot in degrees

  const char *V0right="VZERO/V0C";
  const char *V0left="VZERO/V0A";

  Int_t iIndex=0;
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iIndex);

  dx = rnd->Gaus(0.,sigmatr);
  dy = rnd->Gaus(0.,sigmatr);
  dz = rnd->Gaus(0.,sigmatr);
  dpsi = rnd->Gaus(0.,sigmarot);
  dtheta = rnd->Gaus(0.,sigmarot);
  dphi = rnd->Gaus(0.,sigmarot);
  new(alobj[0]) AliAlignObjAngles(V0right, volid, dx, dy, dz, dpsi, dtheta, dphi, kFALSE);
  dx = rnd->Gaus(0.,sigmatr);
  dy = rnd->Gaus(0.,sigmatr);
  dz = rnd->Gaus(0.,sigmatr);
  dpsi = rnd->Gaus(0.,sigmarot);
  dtheta = rnd->Gaus(0.,sigmarot);
  dphi = rnd->Gaus(0.,sigmarot);
  new(alobj[1]) AliAlignObjAngles(V0left, volid, dx, dy, dz, dpsi, dtheta, dphi,kFALSE);

  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("V0residualMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"V0ResidualObjs ","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Brigitte Cheynis");
    md->SetComment("Alignment objects for V0 residual misalignment");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("VZERO/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}

