void MakeMFTZeroMisAlignment(TString Storage = "alien://folder=/alice/cern.ch/user/a/auras/OCDB/") {

  // Create TClonesArray of zero misalignment objects for MFT

  const char* macroname = "MakeMFTZeroMisAlignment.C";

  TClonesArray *array = new TClonesArray("AliAlignObjParams",10);
  TClonesArray &alobj = *array;

  Double_t dx=0, dy=0, dz=0, dpsi=0, dtheta=0, dphi=0;

  Int_t iIndex=0;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

  TString MFT("MFT");
  new (alobj[0]) AliAlignObjParams(MFT.Data(), volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);

  // save in CDB storage
  if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
    Error(macroname,"STORAGE variable set to %s is not valid. Exiting\n",Storage.Data());
    return;
  }
  Info(macroname,"Saving alignment objects in CDB storage %s", Storage.Data());
  AliCDBManager* cdb = AliCDBManager::Instance();
  AliCDBStorage* storage = cdb->GetStorage(Storage.Data());
  if(!storage){
    Error(macroname,"Unable to open storage %s\n",Storage.Data());
    return;
  }
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Antonio Uras");
  md->SetComment("Alignment objects for MFT zero-misalignment");
  md->SetAliRootVersion(gROOT->GetVersion());
  AliCDBId id("MFT/Align/Data",0,AliCDBRunRange::Infinity());
  storage->Put(array,id,md);

  array->Delete();

}

