void MakeCalibrationSSD(Int_t firstRun=0,Int_t lastRun=9999999 ){
  
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/");
  }
  
  AliCDBMetaData *md1= new AliCDBMetaData(); // metaData describing the object
  md1->SetObjectClassName("AliITSCalibration");
  md1->SetResponsible("Enrico Fragiacomo");
  md1->SetBeamPeriod(0);
  md1->SetAliRootVersion("head 23/08/07"); //root version

  AliCDBId idNoiseSSD("ITS/Calib/NoiseSSD",firstRun, lastRun);
  AliCDBId idGainSSD("ITS/Calib/GainSSD",firstRun, lastRun);
  AliCDBId idBadChannelsSSD("ITS/Calib/BadChannelsSSD",firstRun, lastRun);
  
  TObjArray noiseSSD(1698);
  TObjArray gainSSD(1698);
  TObjArray badchannelsSSD(1698);

  noiseSSD.SetOwner(kFALSE);
  gainSSD.SetOwner(kFALSE);
  badchannelsSSD.SetOwner(kFALSE);
  
  Double_t noiseP[768];
  Double_t noiseN[768];

  Double_t gainP[768];
  Double_t gainN[768];
  
  // loop over SSD modules
  for(Int_t i=0;i<1698;i++){
    
    AliITSNoiseSSD* noise = new AliITSNoiseSSD();
    AliITSGainSSD* gain = new AliITSGainSSD();
    AliITSBadChannelsSSD* badchannels = new AliITSBadChannelsSSD();

    // 768 strips on P- and N-side
    noise->SetNNoiseP(768);
    noise->SetNNoiseN(768);
    gain->SetNGainP(768);
    gain->SetNGainN(768);
    badchannels->SetNBadPChannelsList(10);
    badchannels->SetNBadNChannelsList(10);

    // take a reasonable averaged value for the noise on P- and N-side strips 
    for(Int_t j=0; j<768; j++) {
      noise->AddNoiseP(j,2.);
      gain->AddGainP(j,1.);
      noise->AddNoiseN(j,4.);
      gain->AddGainN(j,1.);
    }

    // 10 random strips per module tagged as "bad"
    for(Int_t j=0; j<10; j++) {
      badchannels->AddBadPChannel(j,gRandom->Uniform(0,767));
      badchannels->AddBadNChannel(j,gRandom->Uniform(0,767));
    }

    noiseSSD.Add(noise);
    gainSSD.Add(gain);
    badchannelsSSD.Add(badchannels);

  }

  AliCDBManager::Instance()->GetDefaultStorage()->Put(&noiseSSD, idNoiseSSD, md1);
  AliCDBManager::Instance()->GetDefaultStorage()->Put(&gainSSD, idGainSSD, md1);
  AliCDBManager::Instance()->GetDefaultStorage()->Put(&badchannelsSSD, idBadChannelsSSD, md1);

}
