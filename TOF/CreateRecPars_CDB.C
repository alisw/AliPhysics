void CreateRecPars_CDB(){
  // Create TOF Calibration Object for Ideal calibration and 
  // write it on CDB
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE");
  AliTOFGeometry *geom = new AliTOFGeometryV5(); 
  AliTOFcalib *tofcalib = new AliTOFcalib(geom);
  AliTOFRecoParam *param = new AliTOFRecoParam();
  tofcalib->WriteRecParOnCDB("TOF/Calib",0,0,param);
}


