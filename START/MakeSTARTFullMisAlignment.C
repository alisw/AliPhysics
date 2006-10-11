void MakeSTARTFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for T0
  //
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",10);
  TClonesArray &alobj = *array;

  if(!gGeoManager) TGeoManager::Import("geometry.root");

  AliAlignObjAngles a;

  Double_t dx, dy, dz, dpsi, dtheta, dphi;
  TRandom *rnd   = new TRandom(4321);
  Double_t sigmatr = 0.05; // max shift in cm w.r.t. RS
  Double_t sigmarot = 0.3; // max rot in degrees w.r.t. RS

  // null shifts and rotations

  const char *T0right="/ALIC_1/0STR_1";
  const char *T0left="/ALIC_1/0STL_1";

  Int_t iIndex=0; //let all modules have index=0 in a layer with no LUT
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iIndex);

  dx = rnd->Gaus(0.,sigmatr);
  dy = rnd->Gaus(0.,sigmatr);
  dz = rnd->Gaus(0.,sigmatr);
  dpsi = rnd->Gaus(0.,sigmarot);
  dtheta = rnd->Gaus(0.,sigmarot);
  dphi = rnd->Gaus(0.,sigmarot);
  cout<<dx<<"  "<<dy<<"  "<<dz<<"  "<<dpsi<<"  "<<dtheta<<"  "<<dphi<<"\n";
  new(alobj[0]) AliAlignObjAngles(T0right, volid, dx, dy, dz, dpsi, dtheta, dphi, kFALSE);
  dx = rnd->Gaus(0.,sigmatr);
  dy = rnd->Gaus(0.,sigmatr);
  dz = rnd->Gaus(0.,sigmatr);
  dpsi = rnd->Gaus(0.,sigmarot);
  dtheta = rnd->Gaus(0.,sigmarot);
  dphi = rnd->Gaus(0.,sigmarot);
  cout<<dx<<"  "<<dy<<"  "<<dz<<"  "<<dpsi<<"  "<<dtheta<<"  "<<dphi<<"\n";
  new(alobj[1]) AliAlignObjAngles(T0left, volid, dx, dy, dz, dpsi, dtheta, dphi, kFALSE);

  // just save in a file you will send me
  // or you will put in your afs public dir ...

  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("T0fullMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"T0FullObjs ","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    // save in CDB storage
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Tomasz Malkiewicz");
    md->SetComment("Full misalignment for T0, produced with sigmatr=0.05 and sigmarot=0.3 in the local RS");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("START/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}

