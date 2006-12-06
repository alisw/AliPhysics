void AliTRDCheckPreprocessor()
{
  // load library
  gSystem->Load("libTestShuttle.so");


  AliTestShuttle::SetOCDBStorage("local://TestCDB_New");
  AliTestShuttle::SetReferenceStorage("local://TestReference_New");
  printf("Test OCDB storage Uri: %s\n", AliTestShuttle::GetOCDBStorage().Data());
  printf("Test Reference storage Uri: %s\n", AliTestShuttle::GetReferenceStorage().Data());

  //Check the reference data 
  //***************************

  

  //Test reference data gain
  //***************************
  Int_t ErrorRefDataGain = 0;
  AliCDBEntry* entry = AliCDBManager::Instance()->GetStorage(AliTestShuttle::GetReferenceStorage())->Get("TRD/HLTData/Gain", 7);  
  if(!entry) ErrorRefDataGain = 1;
  else{
    TH2I *histogain = (TH2I *) entry->GetObject();
    if(!histogain) ErrorRefDataGain = 2;
    else{
      Int_t NbinsX = histogain->GetNbinsX();
      if(NbinsX != 540) ErrorRefDataGain = 3;
    }
  }


  //Test reference data vdriftt0
  //***************************
  Int_t ErrorRefDataVdriftT0 = 0;
  AliCDBEntry* entry = AliCDBManager::Instance()->GetStorage(AliTestShuttle::GetReferenceStorage())->Get("TRD/HLTData/VdriftT0", 7);
  if(!entry) ErrorRefDataVdriftT0 = 1;
  else{
    TProfile2D *histovdriftt0 = (TProfile2D *) entry->GetObject();
    if(!histovdriftt0) ErrorRefDataVdriftT0 = 2;
    else{
      Int_t NbinsX = histovdriftt0->GetNbinsX();
      if(NbinsX != 540) ErrorRefDataVdriftT0 = 3;
    }
  }
  
  
  //Test reference data PRF
  //***************************
  Int_t ErrorRefDataPRF = 0;
  AliCDBEntry* entry = AliCDBManager::Instance()->GetStorage(AliTestShuttle::GetReferenceStorage())->Get("TRD/HLTData/PRF", 7);
  if(!entry) ErrorRefDataPRF = 1;
  else{
    TProfile2D *histoprf = (TProfile2D *) entry->GetObject();
    if(!histoprf) ErrorRefDataPRF = 2;
    else{
      Int_t NbinsX = histoprf->GetNbinsX();
      if(NbinsX != 540) ErrorRefDataPRF = 3;
    }
  }
  

  //Check the detector OCDB values
  //********************************
 
  //Test for pads
  //******************
  //Gain
  //*****
  Int_t ErrorGainPad = 0;
  AliCDBEntry* entry = AliCDBManager::Instance()->GetStorage(AliTestShuttle::GetOCDBStorage())->Get("TRD/Calib/LocalGainFactor", 7);
  AliTRDCalPad *calPad = (AliTRDCalPad *) entry->GetObject();

  for(Int_t det = 0; det < 540; det++){
    
    AliTRDCalROC *calROC = calPad->GetCalROC(det);
    for(Int_t channel =0; channel < calROC->GetNchannels(); channel++){
      if(calROC->GetValue(channel) != 1.0) ErrorGainPad++;
    }//channel loop
  }//det loop

  //Vdrift
  //*****
  Int_t ErrorVdriftPad = 0;
  entry = AliCDBManager::Instance()->GetStorage(AliTestShuttle::GetOCDBStorage())->Get("TRD/Calib/LocalVdrift", 7);
  calPad = (AliTRDCalPad *) entry->GetObject();

  for(Int_t det = 0; det < 540; det++){
    
    AliTRDCalROC *calROC = calPad->GetCalROC(det);
    for(Int_t channel =0; channel < calROC->GetNchannels(); channel++){
      if(calROC->GetValue(channel) != 1.0) ErrorVdriftPad++;
    }//channel loop
  }//det loop

 
  //PRFWidth
  //********
  Int_t ErrorPRFWidthPad = 0;
  entry = AliCDBManager::Instance()->GetStorage(AliTestShuttle::GetOCDBStorage())->Get("TRD/Calib/PRFWidth", 7);
  calPad = (AliTRDCalPad *) entry->GetObject();
  Float_t value = 0.0;

  for(Int_t plane = 0; plane < 6; plane++){

    if(plane == 0) value = 0.515;
    if(plane == 1) value = 0.502;
    if(plane == 2) value = 0.491;
    if(plane == 3) value = 0.481;
    if(plane == 4) value = 0.471;
    if(plane == 5) value = 0.463;

    for(Int_t chamber = 0; chamber < 5; chamber++){
      for(Int_t sector = 0; sector < 18; sector++){
    
	AliTRDCalROC *calROC = calPad->GetCalROC(plane,chamber,sector);
	for(Int_t channel =0; channel < calROC->GetNchannels(); channel++){
	  if((calROC->GetValue(channel) > 1.1*value) || (calROC->GetValue(channel) < 0.9*value)) ErrorPRFWidth++;
	}//channel loop
      }//sector loop
    }//chamber loop
  }//plane loop


  //Test for detector values
  //*************************

  //Gain
  //******
  Int_t ErrorGainDetector = 0;
  entry = AliCDBManager::Instance()->GetStorage(AliTestShuttle::GetOCDBStorage())->Get("TRD/Calib/ChamberGainFactor", 7);  
  AliTRDCalDet *object = (AliTRDCalDet *) entry->GetObject();
  for(Int_t det = 0; det < 540; det++){
    if((object->GetValue(det)> 1.2) || (object->GetValue(det) < 0.8)) ErrorGainDetector++;
  }
  

  //Vdrift
  //******
  Int_t ErrorVdriftDetector = 0;
  entry = AliCDBManager::Instance()->GetStorage(AliTestShuttle::GetOCDBStorage())->Get("TRD/Calib/ChamberVdrift", 7);  
  object = (AliTRDCalDet *) entry->GetObject();
  for(Int_t det = 0; det < 540; det++){
    if((object->GetValue(det)> 1.6) || (object->GetValue(det) < 1.4)) ErrorVdriftDetector++;
  }

  

  //Extract case of T0
  //********************


  //T0
  //*****
  Int_t ErrorT0Pad = 0;
  Int_t ErrorT0Detector = 0;
  Int_t ErrorT0 = 0;
  entry = AliCDBManager::Instance()->GetStorage(AliTestShuttle::GetOCDBStorage())->Get("TRD/Calib/LocalT0", 7);
  calPad = (AliTRDCalPad *) entry->GetObject();
  AliCDBEntry *entry1 = AliCDBManager::Instance()->GetStorage(AliTestShuttle::GetOCDBStorage())->Get("TRD/Calib/ChamberT0", 7);  
  object = (AliTRDCalDet *) entry->GetObject();

  for(Int_t det = 0; det < 540; det++){

    if((object->GetValue(det)> 0.05) || (object->GetValue(det) < 0.0)) ErrorT0Detector++;
    
    AliTRDCalROC *calROC = calPad->GetCalROC(det);
    Float_t valuedetector = calROC->GetValue(0);
    if(((object->GetValue(det)*valuedetector)> 0.05) || ((object->GetValue(det)*valuedetector) < 0.0)) ErrorT0++;
    for(Int_t channel =0; channel < calROC->GetNchannels(); channel++){
      if(calROC->GetValue(channel) != valuedetector) {
	ErrorT0Pad++;
      }
    }//channel loop
  }//det loop

  
  //Bilan
  //**************

  //OCDB values
  //************

  printf("For the local gain factor there are %d strange values\n",ErrorGainPad);
  printf("For the local vdrift there are %d strange values\n",ErrorVdriftPad);
  printf("For the local t0 there are %d strange values\n",ErrorT0Pad);
  printf("For the chamber gain factor there are %d strange values\n",ErrorGainDetector);
  printf("For the chamber vdrift there are %d strange values\n",ErrorVdriftDetector);
  printf("For the chamber t0 there are %d strange values and %d strange product values\n",ErrorT0Detector,ErrorT0);
  printf("For the prf width there are %d strange values\n",ErrorPRFWidthPad);




 //Reference data
 //****************
  
  if(ErrorRefDataGain == 1) printf("There is no reference data entry for the gain!\n");
  
  if(ErrorRefDataGain == 2) printf("There is no reference data histogram for the gain!\n");
  
  if(ErrorRefDataGain == 3) printf("The reference data histogram has not the good number of Xbins for the gain!\n");
  
  if(ErrorRefDataVdriftT0 == 1) printf("There is no reference data entry for the vdriftt0!\n");
  
  if(ErrorRefDataVdriftT0 == 2) printf("There is no reference data histogram for the vdriftt0!\n");
  
  if(ErrorRefDataVdriftT0 == 3) printf("The reference data profile has not the good number of Xbins for the gain!\n");
  
  if(ErrorRefDataPRF == 1) printf("There is no reference data entry for the prf!\n");
  
  if(ErrorRefDataPRF == 2) printf("There is no reference data profile for the prf!\n");
  
  if(ErrorRefDataPRF == 3) printf("The reference data profile has not the good number of Xbins for the prf!\n");
  
 
}
