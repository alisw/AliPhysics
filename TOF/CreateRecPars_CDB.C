void CreateRecPars_CDB(){
  // Create TOF Calibration Object for Ideal calibration and 
  // write it on CDB
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliTOFcalib *tofcalib = new AliTOFcalib();
  AliTOFRecoParam *param = new AliTOFRecoParam();
  tofcalib->WriteRecParOnCDB("TOF/Calib",0,999999999,param);
}


