void MakeTRDFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for TRD
  // 
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",1000);
  TClonesArray &alobj = *array;
   
  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

  AliAlignObjAngles a;

  // sigmas for the chambers
  Double_t chdx=1; // 1 cm
  Double_t chdy=1; // 1 cm
  Double_t chdz=1; // 1 cm
  Double_t chrx=1.0/1000/TMath::Pi()*180; // 1 mrad
  Double_t chry=1.0/1000/TMath::Pi()*180; // 1 mrad
  Double_t chrz=0.7/1000/TMath::Pi()*180; // 0.7 mrad

  Double_t dx,dy,dz,rx,ry,rz;

  Int_t j=0;
  TRandom *ran = new TRandom(4357);
  UShort_t volid;
  const char *path;

  // create the chambers' alignment objects
  for (Int_t iLayer = AliAlignObj::kTRD1; iLayer <= AliAlignObj::kTRD6; iLayer++) {
    for (Int_t iModule = 0; iModule < AliAlignObj::LayerSize(iLayer); iModule++) {
      ran.Rannor(dx,rx);
      ran.Rannor(dy,ry);
      ran.Rannor(dz,rz);
      dx*=chdx;
      dy*=chdy;
      dz*=chdz;
      rx*=chrx;
      ry*=chry;
      rz*=chrz;
      volid = AliAlignObj::LayerToVolUID(iLayer,iModule);
      symname = AliAlignObj::SymName(volid);
      new(alobj[j++]) AliAlignObjAngles(symname,volid,dx,dy,dz,rx,ry,rz,kFALSE);
    }
  }

  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("TRDfullMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"TRDAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Dariusz Miskowiec");
    md->SetComment("Full misalignment for TRD");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("TRD/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}


