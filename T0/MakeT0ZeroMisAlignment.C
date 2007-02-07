void MakeT0ZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for T0
  //
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",10);
  TClonesArray &alobj = *array;

  AliAlignObjAngles a;

  Double_t dx=0, dy=0, dz=0, dpsi=0, dtheta=0, dphi=0;

  TString symName, sn;

  Int_t iIndex=0;
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iIndex);

  Int_t j=0;
  for (Int_t imod=0; imod<24; imod++)
    {
      if (imod < 12){
	sn="T0/C/PMT";
      }else{
	sn="T0/A/PMT";
      }
      symName = sn;
      symName += imod+1;

      new(alobj[j++]) AliAlignObjAngles(symName.Data(), volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
    }

  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("T0zeroMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"T0ZeroObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBStorage* storage = cdb->GetStorage(Storage);
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Tomasz Malkiewicz");
    md->SetComment("Zero misalignment for T0, produced with sigmatr=0.05 and sigmarot=0.3 in the local RS");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("T0/Align/Data",0,9999999);
    storage->Put(array,id,md);
  }

  array->Delete();

}

