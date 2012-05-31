void CreateOnlineCalibPars(){

  // Create TOF Online Calibration Object for reconstruction
  // and write it on CDB;
  // NB: only delay set, status still ok
  // Both old type objects (using TObjArrays) and new type objects (using AliTOFChannelOnlineArray/
  // AliTOFChannelOnlineStatusArray are written

  AliTOFcalib *tofcalib = new AliTOFcalib();
  tofcalib->CreateCalArrays();
  TObjArray *tofCalOnline = (TObjArray*) tofcalib->GetTOFCalArrayOnline(); 
  TObjArray *tofCalOnlinePulser = (TObjArray*) tofcalib->GetTOFCalArrayOnlinePulser(); 
  TObjArray *tofCalOnlineNoise = (TObjArray*) tofcalib->GetTOFCalArrayOnlineNoise(); 
  TObjArray *tofCalOnlineHW = (TObjArray*) tofcalib->GetTOFCalArrayOnlineHW(); 
  AliTOFChannelOnlineArray *delayObj = (AliTOFChannelOnlineArray*) tofcalib->GetTOFOnlineDelay();
  AliTOFChannelOnlineStatusArray *status = (AliTOFChannelOnlineStatusArray*) tofcalib->GetTOFOnlineStatus();
  // Write the online calibration object on CDB

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  Int_t nChannels = AliTOFGeometry::NSectors()*(2*(AliTOFGeometry::NStripC()+AliTOFGeometry::NStripB())+AliTOFGeometry::NStripA())*AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX();
  Float_t delay=0.;
  Float_t meanDelay=0.3;
  Float_t sigmaDelay=0.08;
  TRandom *rnd   = new TRandom(4357);
  for (Int_t ipad = 0 ; ipad<nChannels; ipad++){
    AliTOFChannelOnline *calChannelOnline = (AliTOFChannelOnline*)tofCalOnline->At(ipad);
    AliTOFChannelOnlineStatus *calChannelOnlinePulser = (AliTOFChannelOnlineStatus*)tofCalOnlinePulser->At(ipad);
    AliTOFChannelOnlineStatus *calChannelOnlineNoise = (AliTOFChannelOnlineStatus*)tofCalOnlineNoise->At(ipad);
    AliTOFChannelOnlineStatus *calChannelOnlineHW = (AliTOFChannelOnlineStatus*)tofCalOnlineHW->At(ipad);
    delay = rnd->Gaus(meanDelay,sigmaDelay);
    delayObj->SetDelay(ipad,delay);
    calChannelOnline->SetDelay(delay);
    calChannelOnlinePulser->SetStatus(AliTOFChannelOnlineStatus::kTOFPulserOk);
    calChannelOnlineNoise->SetStatus(AliTOFChannelOnlineStatus::kTOFNoiseOk);
    calChannelOnlineHW->SetStatus(AliTOFChannelOnlineStatus::kTOFHWOk);
    status->SetHWStatus(ipad,AliTOFChannelOnlineStatusArray::kTOFHWOk);
    status->SetPulserStatus(ipad,AliTOFChannelOnlineStatusArray::kTOFPulserOk);
    status->SetNoiseStatus(ipad,AliTOFChannelOnlineStatusArray::kTOFNoiseOk);
  }
  tofcalib->WriteParOnlineDelayOnCDB("TOF/Calib");   // new obj
  tofcalib->WriteParOnlineStatusOnCDB("TOF/Calib");  // new obj
  tofcalib->WriteParOnlineOnCDB("TOF/Calib");        // old object
  tofcalib->WriteParOnlinePulserOnCDB("TOF/Calib");  // old obj
  tofcalib->WriteParOnlineNoiseOnCDB("TOF/Calib");   // old obj
  tofcalib->WriteParOnlineHWOnCDB("TOF/Calib");      // old obj
  return;
}
