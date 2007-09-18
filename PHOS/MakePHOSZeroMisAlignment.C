void MakePHOSZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for PHOS
  //
  const char* macroname = "MakePHOSZeroMisAlignment.C";
  const AliPHOSGeometry *phosGeom = AliPHOSGeometry::GetInstance("IHEP", "IHEP");
  if (!phosGeom) {
    Error("MakePHOSFullMisAlignment", "Cannot obtain AliPHOSGeometry singleton\n");
    return;
  }

  AliPHOSEMCAGeometry *emca = phosGeom->GetEMCAGeometry();
  TClonesArray *array = new TClonesArray("AliAlignObjParams", 16 + phosGeom->GetNModules() * 
                                         emca->GetNStripX() * emca->GetNStripZ());
  TClonesArray &alobj = *array;
   
  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;
  Int_t iIndex=0; // let all modules have index=0 in a layer with no LUT
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);
 
  Int_t i=0 ;
  // Alignment for 5 PHOS modules
  new(alobj[i++]) AliAlignObjParams("PHOS/Module1",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module2",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module3",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module4",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module5",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);

  // Alignment of CPV modules
  new(alobj[i++]) AliAlignObjParams("PHOS/Module1/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module2/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module3/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module4/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Module5/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
 

  // Alignment for PHOS cradle
  new(alobj[i++]) AliAlignObjParams("PHOS/Cradle0",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Cradle1",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);

  // Alignment for cradle wheels
  new(alobj[i++])  AliAlignObjParams("PHOS/Wheel0",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++])  AliAlignObjParams("PHOS/Wheel1",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++])  AliAlignObjParams("PHOS/Wheel2",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjParams("PHOS/Wheel3",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);

  AliPHOSSurvey geodesicData;
  geodesicData.CreateNullObjects(alobj, phosGeom);

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "PHOSzeroMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"PHOSAlignObjs","kSingleKey");
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
    AliCDBMetaData *md= new AliCDBMetaData();
    md->SetResponsible("Yuri Kharlov");
    md->SetComment("Zero misalignment objects");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("PHOS/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id, md);
  }

  array->Delete();

}
