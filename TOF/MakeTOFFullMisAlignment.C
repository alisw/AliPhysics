void MakeTOFFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for TOF
  // 
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",100);
  TClonesArray &alobj = *array;
   
  if(!gGeoManager) TGeoManager::Import("geometry.root"); //needed for
  // the constructors with local coordinates not to fail

  AliAlignObjAngles a;

  Int_t nSMTOF = 18;
  Int_t iIndex=0; //let all modules have index=0 in a layer with no LUT
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iIndex);
  
  Int_t i;
  Int_t j=0;
  Double_t dx, dy, dz, dpsi, dtheta, dphi;
  TRandom *rnd   = new TRandom(2345);
  Double_t sigmatr = 0.4; // max shift in cm w.r.t. local ideal RS
  Double_t sigmarot = 0.06; // max rot in deg w.r.t. local ideal RS (~ 1 mrad)

  for(i=0; i<18; i++) {
    Char_t  path[100];
    sprintf(path,"TOF/sm%02d",i);
    dx = rnd->Gaus(0.,sigmatr);
    dy = 0;
    dz = rnd->Gaus(0.,sigmatr);
    dpsi = 0;
    dtheta = rnd->Gaus(0.,sigmarot);
    dphi = 0.;
    new(alobj[j]) AliAlignObjAngles(path, volid, dx, dy, dz, dpsi, dtheta, dphi,kFALSE);
    alobj[j]->Print();
    j++;
  }

  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("TOFfullMisalignment.root","RECREATE");
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
    md->SetComment("Full misalignment for TOF");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("TOF/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}


