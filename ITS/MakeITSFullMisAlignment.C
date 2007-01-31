void MakeITSFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for ITS
  //
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",4000);
  TClonesArray &alobj = *array;
   
  if(!gGeoManager) TGeoManager::Import("geometry.root");
  // needed for the constructors with local coordinates not to fail

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
  AliAlignObjAngles a;

  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;

  Int_t j = 0;
  new(alobj[j++]) AliAlignObjAngles("ITS", 0, dx, dy, globalZ, dpsi, dtheta, dphi, kTRUE);

  for ( Int_t l = AliAlignObj::kSPD1; l <= AliAlignObj::kSSD2; l++) {
    
    printf("%i modules in layer %i\n", AliAlignObj::LayerSize(l), l);
    for (Int_t iModule = 0; iModule < AliAlignObj::LayerSize(l); iModule++) {

      dpsi   = rnd->Gaus(0., rot);
      dtheta = rnd->Gaus(0., rot);
      dphi   = rnd->Gaus(0., rot);

      AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer; 
      switch (l) {
      case 1: {
	iLayer = AliAlignObj::kSPD1;
	dx = rnd->Gaus(0., spdXY + mecanicalPrec);
	dy = rnd->Gaus(0., spdXY + mecanicalPrec);
	dz = rnd->Gaus(0., spdZ + mecanicalPrec);
      }; break;
      case 2: {
	iLayer = AliAlignObj::kSPD2;
	dx = rnd->Gaus(0., spdXY + mecanicalPrec);
	dy = rnd->Gaus(0., spdXY + mecanicalPrec);
	dz = rnd->Gaus(0., spdZ + mecanicalPrec);
      }; break;
      case 3: {
	iLayer = AliAlignObj::kSDD1;
	dx = rnd->Gaus(0., sddXYZ + mecanicalPrec);
	dy = rnd->Gaus(0., sddXYZ + mecanicalPrec);
	dz = rnd->Gaus(0., sddXYZ + mecanicalPrec);
      }; break;
      case 4: {
	iLayer = AliAlignObj::kSDD2;
	dx = rnd->Gaus(0., sddXYZ + mecanicalPrec);
	dy = rnd->Gaus(0., sddXYZ + mecanicalPrec);
	dz = rnd->Gaus(0., sddXYZ + mecanicalPrec);
      }; break;
      case 5: {
	iLayer = AliAlignObj::kSSD1;
	dx = rnd->Gaus(0., ssdXY + mecanicalPrec);
	dy = rnd->Gaus(0., ssdXY + mecanicalPrec);
	dz = rnd->Gaus(0., ssdZ + mecanicalPrec);
      }; break;
      case 6: {
	iLayer = AliAlignObj::kSSD2;
	dx = rnd->Gaus(0., ssdXY + mecanicalPrec);
	dy = rnd->Gaus(0., ssdXY + mecanicalPrec);
	dz = rnd->Gaus(0., ssdZ + mecanicalPrec);
      }; break;
      default: Printf("Wrong layer index in ITS (%d) !",l);
      };
      UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iModule);
      const char *symname = AliAlignObj::SymName(volid);

      new(alobj[j++]) AliAlignObjAngles(symname, volid, dx, dy, dz, dpsi, dtheta, dphi, kFALSE);

    }
  }

  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("ITSfullMisalignment.root","RECREATE");
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
    md->SetComment("Alignment objects with actual ITS misalignment");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("ITS/Align/Data",0,9999999);
    storage->Put(array,id, md);
  }

  array->Delete();

}



