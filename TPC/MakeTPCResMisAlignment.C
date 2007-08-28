void MakeTPCResMisAlignment(){
  // Create TClonesArray of residual misalignment objects for TPC
  //
  if(!AliGeomManager::GetGeometry()){
    if(!(AliCDBManager::Instance())->IsDefaultStorageSet())
      AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
      AliCDBManager::Instance()->SetRun(0);
    AliGeomManager::LoadGeometry();
  }
  // needed for the constructors with local coordinates not to fail

  TClonesArray *array = new TClonesArray("AliAlignObjParams",100);
  TClonesArray &alobj = *array;
  
  TRandom *rnd   = new TRandom(4357);
  AliAlignObjParams o;
  Double_t dx, dy, dz, dpsi, dtheta, dphi;
  Int_t j = 0;

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
      new(alobj[j]) AliAlignObjParams(symname, volid, dx, dy, dz, dpsi, dtheta, dphi, kFALSE);
      j++;
    }
  }


  const char* macroname = "MakeTPCResMisAlignment.C";
  if( gSystem->Getenv("TOCDB") != TString("kTRUE") ){
    // save on file
    const char* filename = "TPCresidualMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"TPCAlignObjs","kSingleKey");
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
    md->SetResponsible("Marian Ivanov");
    md->SetComment("Residual misalignment for TPC, sigmatr=0.01 and sigmarot=0.6 in the local RS");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("TPC/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}

