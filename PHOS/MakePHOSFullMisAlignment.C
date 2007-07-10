void MakePHOSFullMisAlignment(){
  // Create TClonesArray of full misalignment objects for PHOS
  //
  const AliPHOSGeometry *phosGeom = AliPHOSGeometry::GetInstance("IHEP", "IHEP");
  if (!phosGeom) {
    Error("MakePHOSFullMisAlignment", "Cannot obtain AliPHOSGeometry singleton\n");
    return;
  }

  AliPHOSEMCAGeometry *emca = phosGeom->GetEMCAGeometry();
  TClonesArray *array = new TClonesArray("AliAlignObjAngles", 16 + phosGeom->GetNModules() * 
                                         emca->GetNStripX() * emca->GetNStripZ());
  TClonesArray &alobj = *array;
   
  Double_t dpsi=0., dtheta=0., dphi=0.;
  Double_t displacement = 10;
  Int_t iIndex=0; //let all modules have index=0 in a layer with no LUT
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);
  Int_t i=0 ;

  // Alignment for 5 PHOS modules
  new(alobj[i++]) AliAlignObjAngles("PHOS/Module1",
	  volid, -20., -10.,   0., dpsi, dtheta, 5, kTRUE);
  new(alobj[i++]) AliAlignObjAngles("PHOS/Module2",
	  volid, -10.,   0., -10., dpsi, dtheta, 2, kTRUE);
  new(alobj[i++]) AliAlignObjAngles("PHOS/Module3",
	  volid,   5., -10.,  10., dpsi, dtheta, 0, kTRUE);
  new(alobj[i++]) AliAlignObjAngles("PHOS/Module4",
	  volid, +10.,  -0., -10., dpsi, dtheta, 2, kTRUE);
  new(alobj[i++]) AliAlignObjAngles("PHOS/Module5",
	  volid, +20., -10.,   0., dpsi, dtheta, 5, kTRUE);

  Double_t dx=0., dy=0., dz=0. ;
  // Alignment of CPV modules
  new(alobj[i++]) AliAlignObjAngles("PHOS/Module1/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjAngles("PHOS/Module2/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjAngles("PHOS/Module3/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjAngles("PHOS/Module4/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjAngles("PHOS/Module5/CPV",
        volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
 
  // Alignment for PHOS cradle
  new(alobj[i++]) AliAlignObjAngles("PHOS/Cradle0",
	  volid, 0., 0., -displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjAngles("PHOS/Cradle1",
	  volid, 0., 0., +displacement, dpsi, dtheta, dphi, kTRUE);

  // Alignment for cradle wheels
  new(alobj[i++]) AliAlignObjAngles("PHOS/Wheel0",
	  volid, 0., 0., -displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjAngles("PHOS/Wheel1",
	  volid, 0., 0., -displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjAngles("PHOS/Wheel2",
	  volid, 0., 0., +displacement, dpsi, dtheta, dphi, kTRUE);
  new(alobj[i++]) AliAlignObjAngles("PHOS/Wheel3",
	  volid, 0., 0., +displacement, dpsi, dtheta, dphi, kTRUE);

//  AliPHOSSurvey geodesicData("phos_mod3_survey.txt");
//  geodesicData.CreateAliAlignObjAngles(alobj);

  AliPHOSSurvey1 geodesicData("phos_mod3_survey_EDMS.txt", "T1_");
  geodesicData.CreateAliAlignObjAngles(alobj);

  // *************************    2nd step    ***************

  if(!gSystem->Getenv("TOCDB")){
    // save on file
    TFile f("PHOSfullMisalignment.root","RECREATE");
    if(!f) cerr<<"cannot open file for output\n";
    f.cd();
    f.WriteObject(array,"PHOSFullObjs ","kSingleKey");
    f.Close();
  }else{
    // save in CDB storage
    const char* Storage = gSystem->Getenv("STORAGE");
    AliCDBManager *CDB = AliCDBManager::Instance();
    AliCDBStorage* storage = CDB->GetStorage(Storage);
    AliCDBMetaData *md= new AliCDBMetaData();
    md->SetResponsible("Yuri Kharlov");
    md->SetComment("Alignment objects for fully misaligned geometry");
    md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
    AliCDBId id("PHOS/Align/Data",0,9999999);
    storage->Put(array,id, md);
  }

  array->Delete();

}
