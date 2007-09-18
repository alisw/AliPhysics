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

  AliGeomManager::ELayerID idTOF = AliGeomManager::kTOF;
  Int_t j=0;
  Int_t strId=-1;
  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;
  TRandom *rnd   = new TRandom(4357);
  Double_t sigmatr = 0.1; // sigma (in cm) for shift w.r.t. local ideal RS

  Int_t sActive[18]={0,1,1,0,0,0,1,1,0,1,1,1,1,0,0,1,1,1};

  Int_t nstrA=15;
  Int_t nstrB=19;
  Int_t nstrC=19;
  Int_t nSectors=18;
  Int_t nStrips=nstrA+2*nstrB+2*nstrC;

  for (Int_t isect = 0; isect < nSectors; isect++) {
    for (Int_t istr = 1; istr <= nStrips; istr++) {
    dy = rnd->Gaus(0.,sigmatr);
    dz = rnd->Gaus(0.,sigmatr);
      strId++;
      if( (TString(gSystem->Getenv("PARTGEOM")) == TString("kTRUE")) && !sActive[iSect] ) continue;
      new(alobj[j++]) AliAlignObjParams(AliGeomManager::SymName(idTOF,strId),AliGeomManager::LayerToVolUID(idTOF,strId), dx, dy, dz, dpsi, dtheta, dphi, kFALSE);
    }
  }

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
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


