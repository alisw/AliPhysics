void MakeITSZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for ITS
  //
  TClonesArray *array = new TClonesArray("AliAlignObjParams",4000);
  TClonesArray &alobj = *array;
  const char* macroname = "MakeITSZeroMisAlignment.C";


  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0., globalZ=0.;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer; 

  Int_t j = 0;

  new(alobj[j]) AliAlignObjParams("ITS", 0, dx, dy, globalZ, dpsi, dtheta, dphi, kTRUE);
  j++;

  for ( Int_t l = AliGeomManager::kSPD1; l <= AliGeomManager::kSSD2; l++) {
    
    printf("%i modules in layer %i\n", AliGeomManager::LayerSize(l), l);
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(l); iModule++) {

      iLayer = AliGeomManager::kInvalidLayer; 

      switch (l) {
      case 1: {
	iLayer = AliGeomManager::kSPD1;
      }; break;
      case 2: {
	iLayer = AliGeomManager::kSPD2;
      }; break;
      case 3: {
	iLayer = AliGeomManager::kSDD1;
      }; break;
      case 4: {
	iLayer = AliGeomManager::kSDD2;
      }; break;
      case 5: {
	iLayer = AliGeomManager::kSSD1;
      }; break;
      case 6: {
	iLayer = AliGeomManager::kSSD2;
      }; break;
      default: Printf("Wrong layer index in ITS (%d) !",l);
      };
      UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iModule);
      const char *symname = AliGeomManager::SymName(volid);

      new(alobj[j]) AliAlignObjParams(symname, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
      j++;

    }
  }

  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "ITSzeroMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"ITSAlignObjs","kSingleKey");
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
    md->SetResponsible("Ludovic Gaudichet");
    md->SetComment("Alignment objects with zero ITS misalignment");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("ITS/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id, md);
  }

  array->Delete();

}


