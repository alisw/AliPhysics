void CreateIdealOfflineCalibPars(){
  // Create TOF Dummy Offline Calibration Object for reconstruction
  // and write it on CDB
  AliTOFcalib *tofcalib = new AliTOFcalib();
  tofcalib->CreateCalArrays();
  TObjArray *tofCalOffline = (TObjArray*) tofcalib->GetTOFCalArrayOffline(); 
  // Write the dummy offline calibration object on CDB

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE");
  tofcalib->WriteParOfflineOnCDB("TOF/Calib","ideal");
  return;
}
