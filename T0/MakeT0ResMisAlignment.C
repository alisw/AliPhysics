void MakeT0ResMisAlignment(){
  // Create TClonesArray of residual misalignment objects for T0
  //
  TClonesArray *array = new TClonesArray("AliAlignObjParams",30);
  TClonesArray &alobj = *array;

  if(!AliGeomManager::GetGeometry()){
    if(!(AliCDBManager::Instance())->IsDefaultStorageSet())
      AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
      AliCDBManager::Instance()->SetRun(0);
    AliGeomManager::LoadGeometry();
  }
  // needed for the constructors with local coordinates not to fail

  AliAlignObjParams a;

  Double_t dx, dy, dz, dpsi, dtheta, dphi;
  TRandom *rnd   = new TRandom(4321);
  Double_t sigmatr = 0.05; // max shift in cm
  Double_t sigmarot = 0.3; // max rot in degrees

  TString symName, sn;

  Int_t iIndex=0;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

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
    
    new(alobj[j++]) AliAlignObjParams(symName.Data(), volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }
  
  const char* macroname = "MakeT0ResMisAlignment.C";
  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "T0residualMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"T0AlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    TString Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    Info(macroname,"Saving alignment objects in CDB storage %s",
	 Storage.Data());
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
    if(!storage){
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Tomasz Malkiewicz");
    md->SetComment("Residual misalignment for T0, produced with sigmatr=0.05 and sigmarot=0.3 in the local RS");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("T0/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}

