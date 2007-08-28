void MakeTOFZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for TOF
  //
  TClonesArray *array = new TClonesArray("AliAlignObjParams",2000);
  TClonesArray &alobj = *array;
  const char* macroname = "MakeTOFZeroMisAlignment.C";

   
  AliGeomManager::ELayerID idTOF = AliGeomManager::kTOF;
  Int_t i;
  Int_t j=0;
  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;

  for(i=0; i<AliGeomManager::LayerSize(idTOF); i++) {
    new(alobj[j++]) AliAlignObjParams(AliGeomManager::SymName(idTOF,i), AliGeomManager::LayerToVolUID(idTOF,i), dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }

  if( gSystem->Getenv("TOCDB") != TString("kTRUE") ){
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
    md->SetResponsible("Silvia Arcelli");
    md->SetComment("Zero misalignment for TOF");
    md->SetAliRootVersion("HEAD");
    AliCDBId id("TOF/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}


