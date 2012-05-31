void MakeT0FullMisAlignment(){
  // Create TClonesArray of full misalignment objects for T0
  //
  TClonesArray *array = new TClonesArray("AliAlignObjParams",4);
  TClonesArray &alobj = *array;

  Double_t dx, dy, dz, dpsi, dtheta, dphi;
  TRandom3 *rnd   = new TRandom3(4321);
  Double_t sigmatr = 0.006; // sigma for shifts in cm
  Double_t sigmarot = 0.001; // sigma for tilts in degrees

  TString symName, sn;

  Int_t iIndex=0;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

  for (Int_t imod=0; imod<2; imod++)
  {
    symName="/ALIC_1/0STR_1";
    if(imod==1) symName="/ALIC_1/0STL_1";	
    dx = rnd->Gaus(0.,sigmatr);
    dy = rnd->Gaus(0.,sigmatr);
    dz = rnd->Gaus(0.,sigmatr);
    dpsi = rnd->Gaus(0.,sigmarot);
    dtheta = rnd->Gaus(0.,sigmarot);
    dphi = rnd->Gaus(0.,sigmarot);
    new(alobj[imod]) AliAlignObjParams(symName.Data(), volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  }
    
  const char* macroname = "MakeT0FullMisAlignment.C";
  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "T0fullMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"T0AlignObjs","kSingleKey");
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
    md->SetResponsible("Tomasz Malkiewicz");
    md->SetComment("Full misalignment for T0, produced with sigmatr=0.05 and sigmarot=0.3 in the local RS");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("T0/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

  array->Delete();

}

