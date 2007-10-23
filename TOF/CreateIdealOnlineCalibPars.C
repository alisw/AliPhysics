void CreateIdealOnlineCalibPars(){
  // Create TOF Dummy (delay=0) Offline Calibration Object for reconstruction
  // and write it on CDB
  AliTOFcalib *tofcalib = new AliTOFcalib();
  tofcalib->CreateCalArrays();
  TObjArray *tofCalOnline = (TObjArray*) tofcalib->GetTOFCalArrayOnline(); 
  // Write the dummy offline calibration object on CDB

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE");
  Int_t nChannels = AliTOFGeometry::NSectors()*(2*(AliTOFGeometry::NStripC()+AliTOFGeometry::NStripB())+AliTOFGeometry::NStripA())*AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX();
  for (Int_t ipad = 0 ; ipad<nChannels; ipad++){
    AliTOFChannelOnline *calChannelOnline = (AliTOFChannelOnline*)tofCalOnline->At(ipad);
    Float_t delay = 0.;
    calChannelOnline->SetDelay(delay);
  }
  tofcalib->WriteParOnlineOnCDB("TOF/Calib");
  return;
}
