void MakeTOFResMisAlignment(){
  // Create TClonesArray of residual misalignment objects for TOF
  //
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",2000);
  TClonesArray &alobj = *array;
   
  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

  AliAlignObjAngles a;

  AliAlignObj::ELayerID idTOF = AliAlignObj::kTOF;
  Int_t i;
  Int_t j=0;
  Double_t dx=0.; 
  Double_t dy=0.; 
  Double_t dz=0.;
  Double_t  dpsi, dtheta, dphi;
  TRandom *rnd   = new TRandom(4357);
  Double_t sigmatr = 0.1; // max shift in cm w.r.t. local ideal RS

  for(i=0; i<AliAlignObj::LayerSize(idTOF); i++) {
    dx = 0;
    dy = rnd->Gaus(0.,sigmatr);
    dz = rnd->Gaus(0.,sigmatr);
    dpsi = 0.;
    dtheta = 0.;
    dphi = 0.;
    new(alobj[j]) AliAlignObjAngles(AliAlignObj::SymName(idTOF,i), AliAlignObj::LayerToVolUID(idTOF,i), dx, dy, dz, dpsi, dtheta, dphi, kFALSE);
    j++;
  }

  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("TOFresidualMisalignment.root","RECREATE");
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
    md->SetComment("Residual misalignment for TOF, sigmatr=1mm in the local RS");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("TOF/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}


