void MakeCalibrationSSD(Int_t firstRun=0,Int_t lastRun=999999999 ){
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Macro to generate and store the calibration files for SSD         
  // Modified: Enrico Fragiacomo - 04/05/2015                          
  // Generates:                                                        
  //  1 file with calibration objects (AliITSNoiseSSDv2, AliITSPedestalSSDv2, AliITSGainSSDv2, AliITSBadChannelsSSDv2
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
    AliCDBManager::Instance()->SetDefaultStorage("local://OCDB/");
  }
  
  AliCDBMetaData *md1= new AliCDBMetaData(); // metaData describing the object
  md1->SetObjectClassName("AliITSCalibration");
  md1->SetResponsible("Enrico Fragiacomo");
  md1->SetBeamPeriod(0);
 
  AliCDBId idNoiseSSD("ITS/Calib/NoiseSSD",firstRun, lastRun);
  AliCDBId idPedestalSSD("ITS/Calib/PedestalSSD",firstRun, lastRun);
  AliCDBId idGainSSD("ITS/Calib/GainSSD",firstRun, lastRun);
  AliCDBId idBadChannelsSSD("ITS/Calib/BadChannelsSSD",firstRun, lastRun);
  
  AliITSNoiseSSDv2 *noiseSSD = new AliITSNoiseSSDv2();
  AliITSPedestalSSDv2 *pedestalSSD = new AliITSPedestalSSDv2();
  AliITSGainSSDv2 *gainSSD = new AliITSGainSSDv2();
  AliITSBadChannelsSSDv2 *badchannelsSSD = new AliITSBadChannelsSSDv2();

  for(Int_t i=0; i<1698; i++) {
    for(Int_t j=0; j<768; j++) {
      noiseSSD->AddNoiseP(i,j,3.);
      noiseSSD->AddNoiseN(i,j,5.);
      gainSSD->AddGainP(i,j,0.8);
      gainSSD->AddGainN(i,j,1.2);
    }
  }

  AliCDBManager::Instance()->GetDefaultStorage()->Put( (TObject*) noiseSSD, idNoiseSSD, md1);
  AliCDBManager::Instance()->GetDefaultStorage()->Put( (TObject*) gainSSD, idGainSSD, md1);
  AliCDBManager::Instance()->GetDefaultStorage()->Put( (TObject*) badchannelsSSD, idBadChannelsSSD, md1);
  AliCDBManager::Instance()->GetDefaultStorage()->Put( (TObject*) pedestalSSD, idPedestalSSD, md1);
  
}
