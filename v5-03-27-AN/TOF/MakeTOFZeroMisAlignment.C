void MakeTOFZeroMisAlignment() {
  //
  // Create TClonesArray of zero misalignment objects for TOF
  //

  const char* macroname = "MakeTOFZeroMisAlignment.C";

  TClonesArray *array = new TClonesArray("AliAlignObjParams",2000);
  TClonesArray &alobj = *array;

  // Activate CDB storage to load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  
  AliCDBStorage* storage;

  if ( TString(gSystem->Getenv("TOCDB")) == TString("kTRUE") ) {
    TString Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
      return;
    }
    storage = cdb->GetStorage(Storage.Data());
    if (!storage) {
      Error(macroname,"Unable to open storage %s\n",Storage.Data());
      return;
    }
    AliCDBPath path("GRP","Geometry","Data");
    AliCDBEntry *entry = (AliCDBEntry*)storage->Get(path.GetPath(),cdb->GetRun());
    if (!entry)
      Fatal(macroname,"Could not get the specified CDB entry!");

    entry->SetOwner(0);
    TGeoManager* geom = (TGeoManager*) entry->GetObject();
    AliGeomManager::SetGeometry(geom);
  } else
    AliGeomManager::LoadGeometry(); //load geom from default CDB storage


  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;
  AliGeomManager::ELayerID idTOF = AliGeomManager::kTOF;
  Int_t j=0;
  Int_t nSectors=18;

  Int_t iIndex=0; //let all modules have index=0 in a layer with no LUT
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t dvoluid = AliGeomManager::LayerToVolUID(iLayer,iIndex); //dummy vol id

  for (Int_t isect = 0; isect < nSectors; isect++) {
    TString symname(Form("TOF/sm%02d",isect));
    new(alobj[j++]) AliAlignObjParams(symname.Data(),
				      dvoluid,
				      dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }

  Int_t strId=-1;
  Int_t nstrA=15;
  Int_t nstrB=19;
  Int_t nstrC=19;
  Int_t nStrips=nstrA+2*nstrB+2*nstrC;

  for (Int_t isect = 0; isect < nSectors; isect++) {
    for (Int_t istr = 1; istr <= nStrips; istr++) {
      strId++;
      if ((isect==13 || isect==14 || isect==15) && (istr >= 39 && istr <= 53)) continue;
      new(alobj[j++]) AliAlignObjParams(AliGeomManager::SymName(idTOF,strId),
					AliGeomManager::LayerToVolUID(idTOF,strId),
					dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
    }
  }

  if ( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ) {
    // save on file
    const char* filename = "TOFzeroMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"TOFAlignObjs","kSingleKey");
    f.Close();
  } else {
    // save in CDB storage
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Silvia Arcelli");
    md->SetComment("Zero misalignment for TOF");
    md->SetAliRootVersion("HEAD");
    AliCDBId id("TOF/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}
