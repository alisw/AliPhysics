void CreateCalibPars_CalibHisto(){
  // Create TOF Calibration Object from AliTOFcalibHisto class
  // and write it on CDB

  AliTOFcalib *tofcalib = new AliTOFcalib();
  tofcalib->CreateCalArrays();
  TObjArray *tofCalOffline = (TObjArray*) tofcalib->GetTOFCalArrayOffline(); 

  /* get calib histo andl and load params */
  AliTOFcalibHisto calibHisto;
  calibHisto.LoadCalibPar();

  /* turn time-slewing correction off to only retrieve constants */
  calibHisto.SetFullCorrectionFlag(AliTOFcalibHisto::kTimeSlewingCorr, kFALSE);

  /* OCDB init */
  Float_t par[6] = {0.,0.,0.,0.,0.,0.};
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  Int_t nChannels = AliTOFGeometry::NSectors()*(2*(AliTOFGeometry::NStripC()+AliTOFGeometry::NStripB())+AliTOFGeometry::NStripA())*AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX();
  
  /* common time-slewing params */
  for (Int_t iSlew = 1; iSlew < 6; iSlew++) 
    par[iSlew] = calibHisto.GetCalibPar(AliTOFcalibHisto::kTimeSlewingPar, iSlew);
  
  /* channel-related params */
  for (Int_t ipad = 0 ; ipad<nChannels; ipad++){
    AliTOFChannelOffline *calChannelOffline = (AliTOFChannelOffline*)tofCalOffline->At(ipad);
    par[0] = calibHisto.GetFullCorrection(ipad);
    calChannelOffline->SetSlewPar(par);
  }

  /* write */
  tofcalib->WriteParOfflineOnCDB("TOF/Calib","valid",0,AliCDBRunRange::Infinity());
}


