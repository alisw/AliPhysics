void MakeZDCZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for ZDC
  // 
  const char* macroname = "MakeZDCZeroMisAlignment.C";

  TClonesArray *array = new TClonesArray("AliAlignObjParams",10);
  TClonesArray &alobj = *array;

  Double_t dx=0., dy=0., dz=0.;
  Double_t dpsi=0., dtheta=0., dphi=0.;

  const char *ZDCn="ZDC/NeutronZDC";
  const char *ZDCp="ZDC/ProtonZDC";

  UShort_t iIndex=0;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

  new(alobj[0]) AliAlignObjParams(ZDCn, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[1]) AliAlignObjParams(ZDCp, volid, dx, dy, dz, dpsi, dtheta, dphi,kTRUE);

  if( gSystem->Getenv("TOCDB") != TString("kTRUE") ){
    // save in file
    const char* filename = "ZDCzeroMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"ZDCAlignObjs","kSingleKey");
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
    md->SetResponsible("Chiara Oppedisano");
    md->SetComment("Alignment objects for ZDC zero misalignment");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("ZDC/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}

