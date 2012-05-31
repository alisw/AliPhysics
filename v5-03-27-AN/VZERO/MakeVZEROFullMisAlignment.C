void MakeVZEROFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for VZERO
  // 
  const char* macroname = "MakeVZEROFullMisAlignment.C";

  TClonesArray *array = new TClonesArray("AliAlignObjParams",10);
  TClonesArray &alobj = *array;

  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
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
  Double_t sigmatr = 0.1; // max shift in cm
  Double_t sigmarot = 0.5; // max rot in degrees

  const char *V0right="VZERO/V0C";
  const char *V0left="VZERO/V0A";

  Int_t iIndex=0; //let all modules have index=0 in a layer with no LUT
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

  dx = rnd->Gaus(0.,sigmatr);
  dy = rnd->Gaus(0.,sigmatr);
  dz = rnd->Gaus(0.,sigmatr);
  dpsi = rnd->Gaus(0.,sigmarot);
  dtheta = rnd->Gaus(0.,sigmarot);
  dphi = rnd->Gaus(0.,sigmarot);
  new(alobj[0]) AliAlignObjParams(V0right, volid, dx, dy, dz, dpsi, dtheta,
                                  dphi, kFALSE);				  
  AliAlignObjParams* itsalobj = (AliAlignObjParams*) alobj.UncheckedAt(0);
  itsalobj->ApplyToGeometry();	
  			  
  dx = rnd->Gaus(0.,sigmatr);
  dy = rnd->Gaus(0.,sigmatr);
  dz = rnd->Gaus(0.,sigmatr);
  dpsi = rnd->Gaus(0.,sigmarot);
  dtheta = rnd->Gaus(0.,sigmarot);
  dphi = rnd->Gaus(0.,sigmarot);
  new(alobj[1]) AliAlignObjParams(V0left, volid, dx, dy, dz, dpsi, dtheta,
                                  dphi,kFALSE);
  AliAlignObjParams* itsalobj = (AliAlignObjParams*) alobj.UncheckedAt(1);
  itsalobj->ApplyToGeometry();

  if(TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
     // save on file
     const char* filename = "VZEROfullMisalignment.root";
     TFile f(filename,"RECREATE");
     if(!f){
       Error(macroname,"cannot open file for output\n");
       return;
     }
     Info(macroname,"Saving alignment objects to the file %s", filename);
     f.cd();
     f.WriteObject(array,"VZEROAlignObjs","kSingleKey");
     f.Close();
  }else{
     // save in CDB storage
     AliCDBMetaData* md = new AliCDBMetaData();
     md->SetResponsible("Brigitte Cheynis");
     md->SetComment("Alignment objects for V0 full misalignment");
     md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
     AliCDBId id("VZERO/Align/Data",0,AliCDBRunRange::Infinity());
     storage->Put(array,id,md);
  }

  array->Print();
  array->Delete();

}

