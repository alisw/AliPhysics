void MakeADThresholdCalibEntry(const char *outputCDB = "local://$ALICE_ROOT/OCDB")
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(outputCDB);

  // Creation of the threshold calibration OCDB object
  Double_t a[16] = {3.43814, 3.07214, 3.98605, 3.84432, 3.88496, 3.71621, 3.86557, 3.21873, 3.05334, 1.17287, 2.83189, 3.90145, 3.07584, 2.59038, 3.21505, 4.14664}; 
     
  Double_t b[16] = {7.31771, 7.05281, 9.64445, 8.4491,  8.67408, 8.33752, 8.97269, 7.40062, 7.10964, 2.77174, 6.79971, 9.32446, 7.37046, 6.14143, 7.71174, 9.83657}; 
		    
  TH2F *thrCalib = new TH2F("ADThresholdCalib", "AD threshold calibration parameters", 16, -0.5, 15.5, 2, -0.5, 1.5);
  for(Int_t channel = 0; channel < 16; ++channel) {
    thrCalib->SetBinContent(channel+1,1,a[channel]);
    thrCalib->SetBinContent(channel+1,2,b[channel]);
  }
	
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Michal Broz");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Threshold calibration parameters channel by channel");
  md->AddDateToComment();
  md->PrintMetaData();

  AliCDBId id("AD/Calib/Thresholds", 0, AliCDBRunRange::Infinity());
  man->Put(thrCalib, id, md);

  delete md;

}
