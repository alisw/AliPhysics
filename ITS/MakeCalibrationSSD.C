void MakeCalibrationSSD(Int_t firstRun=0,Int_t lastRun=999999999 ){
  
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB/");
  }
  
  AliCDBMetaData *md1= new AliCDBMetaData(); // metaData describing the object
  md1->SetObjectClassName("AliITSCalibration");
  md1->SetResponsible("Enrico Fragiacomo");
  md1->SetBeamPeriod(0);
  md1->SetAliRootVersion("trunk090708"); //root version

  AliCDBId idNoiseSSD("ITS/Calib/NoiseSSD",firstRun, lastRun);
  AliCDBId idPedestalSSD("ITS/Calib/PedestalSSD",firstRun, lastRun);
  AliCDBId idGainSSD("ITS/Calib/GainSSD",firstRun, lastRun);
  AliCDBId idBadChannelsSSD("ITS/Calib/BadChannelsSSD",firstRun, lastRun);
  
  AliITSNoiseSSD *noiseSSD = new AliITSNoiseSSD();
  AliITSPedestalSSD *pedestalSSD = new AliITSPedestalSSD();
  AliITSGainSSD *gainSSD = new AliITSGainSSD();
  AliITSBadChannelsSSD *badchannelsSSD = new AliITSBadChannelsSSD();

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
