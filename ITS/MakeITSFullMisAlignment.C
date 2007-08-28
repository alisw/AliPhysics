void MakeITSFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for ITS
  //
  TClonesArray *array = new TClonesArray("AliAlignObjParams",4000);
  TClonesArray &alobj = *array;
  const char* macroname = "MakeITSFullMisAlignment.C";
   
  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT");
  cdb->SetRun(0);
  
  AliCDBStorage* storage;
  
  if( gSystem->Getenv("TOCDB") == TString("kTRUE") ){
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

  Double_t globalZ = 0.015; // in cm, = 150 microns
  Double_t mecanicalPrec = 0.0020;

  Double_t resFact = 0.;
  Double_t spdXY   = 0.0015*resFact;
  Double_t sddXYZ  = 0.0030*resFact;
  Double_t ssdXY   = 0.0020*resFact;
  Double_t rot     = 0.018;
 
  Double_t spdZ    = 0.002;
  Double_t ssdZ    = 0.010;


  TRandom *rnd   = new TRandom(65416087);

  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;

  Int_t j = 0;
  new(alobj[j++]) AliAlignObjParams("ITS", 0, dx, dy, globalZ, dpsi, dtheta, dphi, kTRUE);
  AliAlignObjParams* its_alobj = (AliAlignObjParams*) array->UncheckedAt(0);
  its_alobj->ApplyToGeometry();

  for ( Int_t l = AliGeomManager::kSPD1; l <= AliGeomManager::kSSD2; l++) {
    
    printf("%i modules in layer %i\n", AliGeomManager::LayerSize(l), l);
    for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(l); iModule++) {

      dpsi   = rnd->Gaus(0., rot);
      dtheta = rnd->Gaus(0., rot);
      dphi   = rnd->Gaus(0., rot);

      AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer; 
      switch (l) {
      case 1: {
	iLayer = AliGeomManager::kSPD1;
	dx = rnd->Gaus(0., spdXY + mecanicalPrec);
	dy = rnd->Gaus(0., spdXY + mecanicalPrec);
	dz = rnd->Gaus(0., spdZ + mecanicalPrec);
      }; break;
      case 2: {
	iLayer = AliGeomManager::kSPD2;
	dx = rnd->Gaus(0., spdXY + mecanicalPrec);
	dy = rnd->Gaus(0., spdXY + mecanicalPrec);
	dz = rnd->Gaus(0., spdZ + mecanicalPrec);
      }; break;
      case 3: {
	iLayer = AliGeomManager::kSDD1;
	dx = rnd->Gaus(0., sddXYZ + mecanicalPrec);
	dy = rnd->Gaus(0., sddXYZ + mecanicalPrec);
	dz = rnd->Gaus(0., sddXYZ + mecanicalPrec);
      }; break;
      case 4: {
	iLayer = AliGeomManager::kSDD2;
	dx = rnd->Gaus(0., sddXYZ + mecanicalPrec);
	dy = rnd->Gaus(0., sddXYZ + mecanicalPrec);
	dz = rnd->Gaus(0., sddXYZ + mecanicalPrec);
      }; break;
      case 5: {
	iLayer = AliGeomManager::kSSD1;
	dx = rnd->Gaus(0., ssdXY + mecanicalPrec);
	dy = rnd->Gaus(0., ssdXY + mecanicalPrec);
	dz = rnd->Gaus(0., ssdZ + mecanicalPrec);
      }; break;
      case 6: {
	iLayer = AliGeomManager::kSSD2;
	dx = rnd->Gaus(0., ssdXY + mecanicalPrec);
	dy = rnd->Gaus(0., ssdXY + mecanicalPrec);
	dz = rnd->Gaus(0., ssdZ + mecanicalPrec);
      }; break;
      default: Printf("Wrong layer index in ITS (%d) !",l);
      };
      UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iModule);
      const char *symname = AliGeomManager::SymName(volid);

      new(alobj[j++]) AliAlignObjParams(symname, volid, dx, dy, dz, dpsi, dtheta, dphi, kFALSE);

    }
  }

  if( gSystem->Getenv("TOCDB") != TString("kTRUE") ){
    // save on file
    const char* filename = "ITSfullMisalignment.root";
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
    AliCDBMetaData *md= new AliCDBMetaData();
    md->SetResponsible("Ludovic Gaudichet");
    md->SetComment("Alignment objects with actual ITS misalignment");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("ITS/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id, md);
  }

  array->Delete();

}



