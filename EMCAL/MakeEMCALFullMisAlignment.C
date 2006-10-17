void MakeEMCALFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for EMCAL
  //
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",10);
  TClonesArray &alobj = *array;

  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

  AliAlignObjAngles a;

  Double_t dx, dy, dz, dpsi, dtheta, dphi;
  TRandom *rnd   = new TRandom(4321);
  Double_t sigmatr = 0.1; // max shift in cm w.r.t. local RS
  Double_t sigmarot = 0.1; // max rot in degrees w.r.t. local RS

  // null shifts and rotations

  const TString basepath = "EMCAL/FullSupermodule";
  TString pathstr;

  Int_t iIndex=0; //let all modules have index=0 in a layer with no LUT
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iIndex); //dummy volume identity
  Int_t i;

  for(i=0; i<10; i++){

    dx = rnd->Gaus(0.,sigmatr);
    dy = rnd->Gaus(0.,sigmatr);
    dz = rnd->Gaus(0.,sigmatr);
    dpsi = rnd->Gaus(0.,sigmarot);
    dtheta = rnd->Gaus(0.,sigmarot);
    dphi = rnd->Gaus(0.,sigmarot);
    pathstr=basepath;
    pathstr+=(i+1);
    cout<<pathstr<<"  "<<dx<<"  "<<dy<<"  "<<dz<<"  "<<dpsi<<"  "<<dtheta<<"  "<<dphi<<"\n";
    new(alobj[i]) AliAlignObjAngles(pathstr, volid, dx, dy, dz, dpsi, dtheta, dphi, kFALSE);
  }

  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("EMCALfullMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"T0FullObjs ","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Jennifer Clay");
    md->SetComment("Full misalignment for EMCAL");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("EMCAL/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}

