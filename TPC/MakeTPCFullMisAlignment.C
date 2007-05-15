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

  for (Int_t iLayer = AliGeomManager::kTPC1; iLayer <= AliGeomManager::kTPC2; iLayer++) {
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer); iModule++) {

      dx = rnd->Gaus(0,sigmatr);
      dy = rnd->Gaus(0,sigmatr);
      dz = rnd->Gaus(0,sigmatr);
      dpsi = rnd->Gaus(0,sigmarot);
      dtheta = rnd->Gaus(0,sigmarot);
      dphi = rnd->Gaus(0,sigmarot);

      UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iModule);
      const char *symname = AliGeomManager::SymName(volid);
      new(alobj[j++]) AliAlignObjAngles(symname, volid, dx, dy, dz, dpsi, dtheta, dphi, kFALSE);
    }
  }


  if( gSystem->Getenv("TOCDB") != TString("kTRUE") ){
    // save on file
    TFile f("TPCfullMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"TPCAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Marian Ivanov");
    md->SetComment("Full misalignment for TPC, sigmatr=0.01 and sigmarot=0.6 in the local RS");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("TPC/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}

