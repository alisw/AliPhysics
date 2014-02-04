void MakeADZeroMisAlignment()
{
  // Create TClonesArray of zero misalignment objects for AD
  // 
  const char* macroname = "MakeADZeroMisAlignment.C";

  TClonesArray *array = new TClonesArray("AliAlignObjParams",10);
  TClonesArray &alobj = *array;

  Double_t dx=0., dy=0., dz=0.;
  Double_t dpsi=0., dtheta=0., dphi=0.;

  const char * ADC1 = "AD/ADC1";
  const char * ADC2 = "AD/ADC2";
  const char * ADA1 = "AD/ADA1";
  const char * ADA2 = "AD/ADA2";

  UShort_t iIndex=0;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

  new( alobj[0] ) AliAlignObjParams( ADC1, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE );
  new( alobj[1] ) AliAlignObjParams( ADC2, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE );
  new( alobj[2] ) AliAlignObjParams( ADA1, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE );
  new( alobj[3] ) AliAlignObjParams( ADA2, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE );

  if ( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ) {
     
    // save in file
    const char* filename = "ADzeroMisalignment.root";
    TFile f( filename, "RECREATE");
    if (!f) {
      Error( macroname, "cannot open file for output\n");
      return;
    }
    Info( macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject( array, "ADAlignObjs", "kSingleKey");
    f.Close();
    
  } else {
     
    // save in CDB storage
    TString storage = gSystem->Getenv("STORAGE");
    if ( !storage.BeginsWith("local://") && !storage.BeginsWith("alien://") ) {
      Error( macroname, "STORAGE variable set to %s is not valid. Exiting\n", storage.Data());
      return;
    }
    Info( macroname, "Saving alignment objects in CDB storage %s", storage.Data() );
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storageCDB = cdb->GetStorage(storage.Data());
    if ( !storageCDB ) {
      Error( macroname,"Unable to open storage %s\n", storageCDB.Data());
      return;
    }
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible( "Ernesto Lopez" );
    md->SetComment( "Alignment objects for AD zero misalignment" );
    md->SetAliRootVersion( gSystem->Getenv("ARVERSION") );
    AliCDBId id( "AD/Align/Data", 0, AliCDBRunRange::Infinity() );
    storageCDB->Put( array, id, md );
  }

  array->Delete();

}

