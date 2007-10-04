void CreateCalibPars_Ideal(){
  // Create TOF Calibration Object for Ideal calibration and 
  // write it on CDB
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE");
  AliTOFcalib *tofcalib = new AliTOFcalib();
  tofcalib->CreateSimCalArrays();
  TObjArray *tofCalOnline = (TObjArray*) tofcalib->GetTOFSimCalArrayOnline(); 
  TObjArray *tofCalOffline = (TObjArray*) tofcalib->GetTOFSimCalArrayOffline(); 
  TH1F *hToT= new TH1F(); //"empty" ToT histo as a default for ideal 
  tofcalib->WriteSimParOnlineOnCDB("TOF/Calib",0,0,tofCalOnline);
  tofcalib->WriteSimParOfflineOnCDB("TOF/Calib","valid",0,0,tofCalOffline,hToT);
}


