void CreateOnlineCalibPars_CalibHisto(){
  // Create TOF Calibration Object from AliTOFcalibHisto class
  // and write it on CDB

  AliTOFcalib *tofcalib = new AliTOFcalib();
  tofcalib->CreateCalArrays();
  AliTOFChannelOnlineArray *delayArray = (AliTOFChannelOnlineArray*) tofcalib->GetTOFOnlineDelay();

  /* get calib histo andl and load params */
  AliTOFcalibHisto calibHisto;
  calibHisto.LoadCalibPar();

  /* turn time-slewing correction off to only retrieve constants */
  calibHisto.SetFullCorrectionFlag(AliTOFcalibHisto::kTimeSlewingCorr, kFALSE);

  /* OCDB init */
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  Int_t nChannels = AliTOFGeometry::NSectors()*(2*(AliTOFGeometry::NStripC()+AliTOFGeometry::NStripB())+AliTOFGeometry::NStripA())*AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX();
  
  /* channel-related params */
  Double_t delay;
  for (Int_t ipad = 0 ; ipad<nChannels; ipad++){
    AliTOFChannelOnline *calChannelOnline = (AliTOFChannelOnline *)tofCalOnline->At(ipad);
    delay = calibHisto.GetFullCorrection(ipad);
    delayArray->SetDelay(ipad, delay);
  }

  /* write */
  tofcalib->WriteParOnlineDelayOnCDB("TOF/Calib",0,AliCDBRunRange::Infinity());
}


