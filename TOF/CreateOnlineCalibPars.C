void CreateOnlineCalibPars(){
  // Create TOF Offline Calibration Object for reconstruction
  // and write it on CDB
  AliTOFcalib *tofcalib = new AliTOFcalib();
  tofcalib->CreateCalArrays();
  TObjArray *tofCalOnline = (TObjArray*) tofcalib->GetTOFCalArrayOnline(); 
  // Write the offline calibration object on CDB

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE");
  Int_t nChannels = AliTOFGeometry::NSectors()*(2*(AliTOFGeometry::NStripC()+AliTOFGeometry::NStripB())+AliTOFGeometry::NStripA())*AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX();
  Float_t delay=0.;
  Float_t meanDelay=0.3;
  Float_t sigmaDelay=0.08;
  TRandom *rnd   = new TRandom(4357);
  for (Int_t ipad = 0 ; ipad<nChannels; ipad++){
    AliTOFChannelOnline *calChannelOnline = (AliTOFChannelOnline*)tofCalOnline->At(ipad);
    delay = rnd->Gaus(meanDelay,sigmaDelay);
    calChannelOnline->SetDelay(delay);
  }
  tofcalib->WriteParOnlineOnCDB("TOF/Calib");
  return;
}
