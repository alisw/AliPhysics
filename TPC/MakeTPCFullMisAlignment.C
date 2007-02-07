void MakeTPCFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for TPC
  //
  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

  TClonesArray *array = new TClonesArray("AliAlignObjAngles",100);
  TClonesArray &alobj = *array;
  
  TRandom *rnd   = new TRandom(4357);
  AliAlignObjAngles o;
  Int_t j = 0;
  Double_t dx, dy, dz, dpsi, dtheta, dphi;

  // RS = local
  // sigma translation = 0.1 mm
  // sigma rotation = 0.1 mrad
  Float_t sigmatr=0.01;
  Float_t sigmarot = 0.006;

  for (Int_t iLayer = AliAlignObj::kTPC1; iLayer <= AliAlignObj::kTPC2; iLayer++) {
    for (Int_t iModule = 0; iModule < AliAlignObj::LayerSize(iLayer); iModule++) {

      dx = (rnd->Uniform()-0.5)*sigmatr;
      dy = (rnd->Uniform()-0.5)*sigmatr;
      dz = (rnd->Uniform()-0.5)*sigmatr;
      dpsi = (rnd->Uniform()-0.5)*sigmarot;
      dtheta = (rnd->Uniform()-0.5)*sigmarot;
      dphi = (rnd->Uniform()-0.5)*sigmarot;

      UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iModule);
      const char *symname = AliAlignObj::SymName(volid);
      new(alobj[j++]) AliAlignObjAngles(symname, volid, dx, dy, dz, dpsi, dtheta, dphi, kFALSE);
    }
  }


  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("TPCfullMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"TPCAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Marian Ivanov");
    md->SetComment("Full misalignment for TPC, sigmatr=0.01 and sigmarot=0.6 in the local RS");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("TPC/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}

