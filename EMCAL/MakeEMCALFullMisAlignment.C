void MakeEMCALFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for EMCAL
  //
  TClonesArray *array = new TClonesArray("AliAlignObjParams",10);
  TClonesArray &alobj = *array;
  const char* macroname = "MakeEMCALFullMisAlignment.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT");
  cdb->SetRun(0);
  
  AliCDBStorage* storage;
  
  if( TString(gSystem->Getenv("TOCDB")) == TString("kTRUE") ){
    TString Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    storage = cdb->GetStorage(Storage.Data());
    if(!storage){
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }
    AliCDBPath path("GRP","Geometry","Data");
    AliCDBEntry *entry = storage->Get(path.GetPath(),cdb->GetRun());
    if(!entry) Fatal(macroname,"Could not get the specified CDB entry!");
    entry->SetOwner(0);
    TGeoManager* geom = (TGeoManager*) entry->GetObject();
    AliGeomManager::SetGeometry(geom);
  }else{
    AliGeomManager::LoadGeometry(); //load geom from default CDB storage
  }    

  Double_t dx, dy, dz, dpsi, dtheta, dphi;
  TRandom *rnd   = new TRandom(4321);
  Double_t sigmatr = 0.1; // max shift in cm w.r.t. local RS
  Double_t sigmarot = 0.1; // max rot in degrees w.r.t. local RS

  // null shifts and rotations

  const TString basepath = "EMCAL/FullSupermodule";
  const TString hbasepath = "EMCAL/HalfSupermodule";
  TString pathstr;

  Int_t iIndex=0; //let all modules have index=0 in a layer with no LUT
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex); //dummy volume identity
  Int_t i;
  Int_t j=0;

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
    new(alobj[j++]) AliAlignObjParams(pathstr, volid, dx, dy, dz, dpsi, dtheta, dphi, kFALSE);
  }

  for(i=0; i<2; i++){
    dx = rnd->Gaus(0.,sigmatr);
    dy = rnd->Gaus(0.,sigmatr);
    dz = rnd->Gaus(0.,sigmatr);
    dpsi = rnd->Gaus(0.,sigmarot);
    dtheta = rnd->Gaus(0.,sigmarot);
    dphi = rnd->Gaus(0.,sigmarot);
    pathstr=hbasepath;
    pathstr+=(i+1);
    new(alobj[j++]) AliAlignObjParams(pathstr, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "EMCALfullMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"EMCALAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Jennifer Clay");
    md->SetComment("Full misalignment for EMCAL");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("EMCAL/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}

