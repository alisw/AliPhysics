void MakePHOSZeroMisAlignment(){
  // Create TClonesArray of zero misalignment objects for PHOS
  //
  TClonesArray *array = new TClonesArray("AliAlignObjAngles",11);
  TClonesArray &alobj = *array;
   
  AliAlignObjAngles a;

  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;

  Int_t iIndex=0; // let all modules have index=0 in a layer with no LUT
  AliAlignObj::ELayerID iLayer = AliAlignObj::kInvalidLayer;
  UShort_t volid = AliAlignObj::LayerToVolUID(iLayer,iIndex);
 
  // Alignment for 5 PHOS modules
  new(alobj[0]) AliAlignObjAngles("PHOS/Module1",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[1]) AliAlignObjAngles("PHOS/Module2",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[2]) AliAlignObjAngles("PHOS/Module3",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[3]) AliAlignObjAngles("PHOS/Module4",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[4]) AliAlignObjAngles("PHOS/Module5",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);

  // Alignment for PHOS cradle
  new(alobj[5]) AliAlignObjAngles("PHOS/Cradle0",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[6]) AliAlignObjAngles("PHOS/Cradle1",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);

  // Alignment for cradle wheels
  new(alobj[7])  AliAlignObjAngles("PHOS/Wheel0",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[8])  AliAlignObjAngles("PHOS/Wheel1",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[9])  AliAlignObjAngles("PHOS/Wheel2",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[10]) AliAlignObjAngles("PHOS/Wheel3",
	volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);


  if(!gSystem->Getenv("$TOCDB")){
    // save on file
    TFile f("PHOSzeroMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"PHOSAlignObjs","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("$STORAGE");
    AliCDBManager *CDB = AliCDBManager::Instance();
    AliCDBStorage* storage = CDB->GetStorage(Storage);
    AliCDBMetaData *md= new AliCDBMetaData();
    md->SetResponsible("Yuri Kharlov");
    md->SetComment("Zero misalignment objects");
    md->SetAliRootVersion(gSystem->Getenv("$ARVERSION"));
    AliCDBId id("PHOS/Align/Data",0,9999999);
    storage->Put(array,id, md);
  }

  array->Delete();

}
