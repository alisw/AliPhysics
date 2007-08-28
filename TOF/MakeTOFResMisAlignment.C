void MakeTOFResMisAlignment(){
  // Create TClonesArray of residual misalignment objects for TOF
  //
  TClonesArray *array = new TClonesArray("AliAlignObjParams",2000);
  TClonesArray &alobj = *array;
   
  const char* macroname = "MakeTOFResMisAlignment.C";

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT");
  cdb->SetRun(0);
  
  AliCDBStorage* storage;
  
  if( gSystem->Getenv("TOCDB") == TString("kTRUE") ){
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

  AliGeomManager::ELayerID idTOF = AliGeomManager::kTOF;
  Int_t i;
  Int_t j=0;
  Double_t dx=0.; 
  Double_t dy=0.; 
  Double_t dz=0.;
  Double_t  dpsi, dtheta, dphi;
  TRandom *rnd   = new TRandom(4357);
  Double_t sigmatr = 0.1; // sigma (in cm) for shift w.r.t. local ideal RS

  for(i=0; i<AliGeomManager::LayerSize(idTOF); i++) {
    dx = 0;
    dy = rnd->Gaus(0.,sigmatr);
    dz = rnd->Gaus(0.,sigmatr);
    dpsi = 0.;
    dtheta = 0.;
    dphi = 0.;
    new(alobj[j]) AliAlignObjParams(AliGeomManager::SymName(idTOF,i), AliGeomManager::LayerToVolUID(idTOF,i), dx, dy, dz, dpsi, dtheta, dphi, kFALSE);
    j++;
  }

  if( gSystem->Getenv("TOCDB") != TString("kTRUE") ){
    // save on file
    const char* filename = "TOFresidualMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"TOFAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Silvia Arcelli");
    md->SetComment("Residual misalignment for TOF, sigmatr=1mm in the local RS");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("TOF/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}


