void MakeFITZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for FIT
  //
  TClonesArray *array = new TClonesArray("AliAlignObjParams",4);
  TClonesArray &alobj = *array;

  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;

  TString symName;

  Int_t iIndex=0;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);
  cout<<" volid "<<volid<<" iLayer "<<iLayer<<" iIndex "<< iIndex<<endl;
  for (Int_t imod=0; imod<2; imod++)
  {
    symName="/ALIC_1/0STR_1";
    if(imod==1) 
    {
      symName="/ALIC_1/0STL_1";
    }
    new(alobj[imod]) AliAlignObjParams(symName.Data(), volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
    cout<<" align "<<symName.Data()<<" "<<volid<<" "<< dx <<" "<<dy<<" "<< dz<<" "<< dpsi<<" "<< dtheta<<" "<< dphi<<endl;
  }

  const char* macroname = "MakeFITZeroMisAlignment.C";
  /*  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "FITzeroMisalignment.root";
    TFile f(filename,"RECREATE");
    if(!f){
      Error(macroname,"cannot open file for output\n");
      return;
    }
    Info(macroname,"Saving alignment objects to the file %s", filename);
    f.cd();
    f.WriteObject(array,"FITAlignObjs","kSingleKey");
    f.Close();
  }else{
  */
    // save in CDB storage
  //   TString Storage = gSystem->Getenv("STORAGE");
  /* AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Alla");
  TString fPath="FIT/Align/Data";
  cout<<fPath.Data()<<endl;
 AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
  if(storage) {
     AliCDBId id(fPath.Data(),firstRun, lastRun );
    storage->Put(clb, id, &md);
  }
  */
   TString Storage = "local://$ALICE_ROOT/OCDB";
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
    md->SetComment("Zero misalignment for FIT, produced with sigmatr=0.05 and sigmarot=0.3 in the local RS");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("FIT/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
    // }

  array->Delete();

}

