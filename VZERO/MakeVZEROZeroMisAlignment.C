void MakeVZEROZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for VZERO
  //
  const char* macroname = "MakeVZEROZeroMisAlignment.C";

  TClonesArray *array = new TClonesArray("AliAlignObjParams",10);
  TClonesArray &alobj = *array;

  Double_t dx, dy, dz, dpsi, dtheta, dphi;

  Int_t iIndex=0;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

  TString V0right("VZERO/V0C");
  new(alobj[0]) AliAlignObjParams(V0right.Data(), volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  TString V0left("VZERO/V0A");
  new(alobj[1]) AliAlignObjParams(V0left.Data(), volid, dx, dy, dz, dpsi, dtheta, dphi,kTRUE);

  if( gSystem->Getenv("TOCDB") != TString("kTRUE") ){
    // save on file
    const char* filename = "VZEROzeroMisalignment.root";
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
    md->SetResponsible("Brigitte Cheynis");
    md->SetComment("Alignment objects for V0 zero-misalignment");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("VZERO/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}

