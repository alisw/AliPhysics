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
  new(alobj[0]) AliAlignObjAngles("PHOS/Module1",
        volid, -0.20, -0.1, +0.0, dpsi, dtheta, 0.2, kTRUE);
  new(alobj[1]) AliAlignObjAngles("PHOS/Module2",
	volid, -0.10, +0.0, -0.2, dpsi, dtheta, 0.2, kTRUE);
  new(alobj[2]) AliAlignObjAngles("PHOS/Module3",
	volid,  0.05, -0.1,  0.2, dpsi, dtheta, 0.0, kTRUE);
  new(alobj[3]) AliAlignObjAngles("PHOS/Module4",
	volid, +0.10, -0.0, -0.1, dpsi, dtheta, 0.1, kTRUE);
  new(alobj[4]) AliAlignObjAngles("PHOS/Module5",
	volid, +0.20, -0.1,  0.1, dpsi, dtheta, 0.2, kTRUE);

  // Alignment for PHOS cradle
  new(alobj[5]) AliAlignObjAngles("PHOS/Cradle0",
	volid, 0., 0., -displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[6]) AliAlignObjAngles("PHOS/Cradle1",
	volid, 0., 0., +displacement, dpsi, dtheta, dphi, kTRUE);

  // Alignment for cradle wheels
  new(alobj[7])  AliAlignObjAngles("PHOS/Wheel0",
	volid, 0., 0., -displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[8])  AliAlignObjAngles("PHOS/Wheel1",
	volid, 0., 0., -displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[9])  AliAlignObjAngles("PHOS/Wheel2",
	volid, 0., 0., +displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[10]) AliAlignObjAngles("PHOS/Wheel3",
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
