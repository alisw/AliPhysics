void MakeEMCALZeroMisAlignment(TString geoname = "EMCAL_COMPLETE12SMV1_DCAL_8SM"){
  // Create TClonesArray of zero misalignment objects for EMCAL
  //
  const char* macroname = "MakeEMCALZeroMisAlignment.C";
  if(geoname=="")geoname=AliEMCALGeometry::GetDefaultGeometryName();
  const AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance(geoname,"");
  if(!geom) {
    Error("MakeEMCALZeroMisAlignment","Cannot obtain AliEMCALGeometry singleton\n");
    return;
  }

  TClonesArray *array = new TClonesArray("AliAlignObjParams",geom->GetNumberOfSuperModules());
  TClonesArray &alobj = *array;

  /*
  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;

  const TString fbasepath = "EMCAL/FullSupermodule";
  const TString hbasepath = "EMCAL/HalfSupermodule";
  TString pathstr;

  Int_t iIndex=0; // let all modules have index=0 in a layer with no LUT
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

  Int_t i;
  Int_t j=0;

  for(i=0; i<10; i++){
    pathstr=fbasepath;
    pathstr+=(i+1);
    new(alobj[j++]) AliAlignObjParams(pathstr, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }

  for(i=0; i<2; i++){
    pathstr=hbasepath;
    pathstr+=(i+1);
    new(alobj[j++]) AliAlignObjParams(pathstr, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }
  */

  AliEMCALSurvey emcalSurvey;
  emcalSurvey.CreateNullObjects(alobj,geom);

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "EMCALzeroMisalignment.root";
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
    md->SetResponsible("Jennifer Klay");
    md->SetComment("Zero misalignment for EMCAL");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("EMCAL/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}

