void MakeT0ResMisAlignment(){
  // Create TClonesArray of residual misalignment objects for T0
  //
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",30);
  TClonesArray &alobj = *array;

  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

  AliAlignObjAngles a;

  Double_t dx, dy, dz, dpsi, dtheta, dphi;
  TRandom *rnd   = new TRandom(4321);
  Double_t sigmatr = 0.05; // max shift in cm
  Double_t sigmarot = 0.3; // max rot in degrees

  TString symName, sn;

  Int_t iIndex=0;
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iIndex);

  Int_t j=0;
  for (Int_t imod=0; imod<24; imod++){
    if (imod < 12){
      sn="T0/C/PMT";
    }else{
      sn="T0/A/PMT";
    }
    symName = sn;
    symName += imod+1;
    
    dx = rnd->Gaus(0.,sigmatr);
    dy = rnd->Gaus(0.,sigmatr);
    dz = rnd->Gaus(0.,sigmatr);
    dpsi = rnd->Gaus(0.,sigmarot);
    dtheta = rnd->Gaus(0.,sigmarot);
    dphi = rnd->Gaus(0.,sigmarot);
    
    new(alobj[j++]) AliAlignObjAngles(symName.Data(), volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }
  
  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("T0residualMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"T0ResidualObjs ","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Tomasz Malkiewicz");
    md->SetComment("Residual misalignment for T0, produced with sigmatr=0.05 and sigmarot=0.3 in the local RS");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("T0/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}

