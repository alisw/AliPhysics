void MakeTRDResMisAlignment(){
  // Create TClonesArray of residual misalignment objects for TRD
  //
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",1000);
  TClonesArray &alobj = *array;
   
  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

  AliAlignObjAngles a;

  // sigmas for the chambers
  Double_t chdx=0.002; // 20 microns
  Double_t chdy=0.003; // 30 microns
  Double_t chdz=0.007; // 70 microns
  Double_t chrx=0.3/1000/TMath::Pi()*180; // 0.3 mrad
  Double_t chry=0.3/1000/TMath::Pi()*180; // 0.3 mrad
  Double_t chrz=0.1/1000/TMath::Pi()*180; // 0.1 mrad

  Double_t dx,dy,dz,rx,ry,rz;

  Int_t j=0;
  TRandom *ran = new TRandom(4357);
  UShort_t volid;
  const char *path;

  // create the chambers' alignment objects
  for (Int_t iLayer = AliGeomManager::kTRD1; iLayer <= AliGeomManager::kTRD6; iLayer++) {
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(iLayer); iModule++) {
      ran.Rannor(dx,rx);
      ran.Rannor(dy,ry);
      ran.Rannor(dz,rz);
      dx*=chdx;
      dy*=chdy;
      dz*=chdz;
      rx*=chrx;
      ry*=chry;
      rz*=chrz;
      volid = AliGeomManager::LayerToVolUID(iLayer,iModule);
      symname = AliGeomManager::SymName(volid);
      new(alobj[j++]) AliAlignObjAngles(symname,volid,dx,dy,dz,rx,ry,rz,kFALSE);
    }
  }

  if( gSystem->Getenv("TOCDB") != TString("kTRUE") ){
    // save on file
    TFile f("TRDresidualMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"TRDAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Dariusz Miskowiec");
    md->SetComment("Residual misalignment for TRD");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("TRD/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}


