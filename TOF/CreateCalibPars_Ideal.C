void CreateCalibPars_Ideal(){
  // Create TOF Calibration Object for Ideal calibration and 
  // write it on CDB
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliTOFcalib *tofcalib = new AliTOFcalib();
  tofcalib->CreateCalArrays();
  TObjArray *tofCalOffline = (TObjArray*) tofcalib->GetTOFCalArrayOffline(); 
  TH1F *hToT= new TH1F(); //"empty" ToT histo as a default for ideal 
  tofcalib->WriteParOfflineOnCDB("TOF/Calib","valid",0,AliCDBRunRange::Infinity());
  tofcalib->WriteSimHistoOnCDB("TOF/Calib",0,AliCDBRunRange::Infinity(),hToT);
}


