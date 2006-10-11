void MakePHOSResMisAlignment(){
  // Create TClonesArray of residual misalignment objects for PHOS
  //
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",11);
  TClonesArray &alobj = *array;
   
  AliAlignObjAngles a;

  Double_t dpsi=0., dtheta=0., dphi=0.;
  Double_t displacement = 0.2;

  Int_t iIndex=0; // let all modules have index=0 in a layer with no LUT
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iIndex);

  // Alignment for 5 PHOS modules
  new(alobj[0]) AliAlignObjAngles("/ALIC_1/PHOS_1",
        volid, -0.20, -0.1, +0.0, dpsi, dtheta, 0.2, kTRUE);
  new(alobj[1]) AliAlignObjAngles("/ALIC_1/PHOS_2",
	volid, -0.10, +0.0, -0.2, dpsi, dtheta, 0.2, kTRUE);
  new(alobj[2]) AliAlignObjAngles("/ALIC_1/PHOS_3",
	volid,  0.05, -0.1,  0.2, dpsi, dtheta, 0.0, kTRUE);
  new(alobj[3]) AliAlignObjAngles("/ALIC_1/PHOS_4",
	volid, +0.10, -0.0, -0.1, dpsi, dtheta, 0.1, kTRUE);
  new(alobj[4]) AliAlignObjAngles("/ALIC_1/PHOS_5",
	volid, +0.20, -0.1,  0.1, dpsi, dtheta, 0.2, kTRUE);

  // Alignment for PHOS cradle
  new(alobj[5]) AliAlignObjAngles("/ALIC_1/PCRA_0",
	volid, 0., 0., -displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[6]) AliAlignObjAngles("/ALIC_1/PCRA_1",
	volid, 0., 0., +displacement, dpsi, dtheta, dphi, kTRUE);

  // Alignment for cradle wheels
  new(alobj[7])  AliAlignObjAngles("/ALIC_1/PWHE_0",
	volid, 0., 0., -displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[8])  AliAlignObjAngles("/ALIC_1/PWHE_1",
	volid, 0., 0., -displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[9])  AliAlignObjAngles("/ALIC_1/PWHE_2",
	volid, 0., 0., +displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[10]) AliAlignObjAngles("/ALIC_1/PWHE_3",
	volid, 0., 0., +displacement, dpsi, dtheta, dphi, kTRUE);


  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("PHOSresidualMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"PHOSResidualObjs ","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager *CDB = AliCDBManager::Instance();
    AliCDBStorage* storage = CDB->GetStorage(Storage);
    AliCDBMetaData *md= new AliCDBMetaData();
    md->SetResponsible("Yuri Kharlov");
    md->SetComment("Alignment objects for slightly misaligned geometry (residual misalignment");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("PHOS/Align/Data",0,9999999);
    storage->Put(array,id, md);
  }

  array->Delete();

}
