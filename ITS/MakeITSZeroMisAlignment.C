void MakeITSZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for ITS
  //
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",4000);
  TClonesArray &alobj = *array;
   
  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

  AliAlignObjAngles a;

  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0., globalZ=0.;
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer; 

  Int_t j = 0;

  new(alobj[j]) AliAlignObjAngles("ITS", 0, dx, dy, globalZ, dpsi, dtheta, dphi, kTRUE);
  j++;

  for ( Int_t l = AliAlignObj::kSPD1; l <= AliAlignObj::kSSD2; l++) {
    
    printf("%i modules in layer %i\n", AliAlignObj::LayerSize(l), l);
    for (Int_t iModule = 0; iModule < AliAlignObj::LayerSize(l); iModule++) {

      iLayer = AliAlignObj::kInvalidLayer; 

      switch (l) {
      case 1: {
	iLayer = AliAlignObj::kSPD1;
      }; break;
      case 2: {
	iLayer = AliAlignObj::kSPD2;
      }; break;
      case 3: {
	iLayer = AliAlignObj::kSDD1;
      }; break;
      case 4: {
	iLayer = AliAlignObj::kSDD2;
      }; break;
      case 5: {
	iLayer = AliAlignObj::kSSD1;
      }; break;
      case 6: {
	iLayer = AliAlignObj::kSSD2;
      }; break;
      default: Printf("Wrong layer index in ITS (%d) !",l);
      };
      UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iModule);
      const char *symname = AliAlignObj::SymName(volid);

      new(alobj[j]) AliAlignObjAngles(symname, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
      j++;

    }
  }

  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("ITSzeroMisalignment.root","RECREATE");
    if(!f) {cerr<<"cannot open file for output\n";}
    f.WriteObject(array,"ITSAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager *CDB = AliCDBManager::Instance();
    AliCDBStorage* storage = CDB->GetStorage(Storage);
    AliCDBMetaData *md= new AliCDBMetaData();
    md->SetResponsible("Ludovic Gaudichet");
    md->SetComment("Alignment objects with zero ITS misalignment");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("ITS/Align/Data",0,9999999);
    storage->Put(array,id, md);
  }

  array->Delete();

}


