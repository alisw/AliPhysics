void MakeADPMGainsEntry(const char *outputCDB = "local://$ALICE_ROOT/../AliRoot/OCDB")
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(outputCDB);

  // Creation of the PM gains OCDB object
  Double_t a[16] = {875.876, 861.803, 824.828, 949.700, 947.877, 934.311, 821.803, 781.871,
  		    957.520, 946.115, 866.103, 902.057, 1028.79, 976.087, 961.054, 927.773};
  Double_t b[16] = {5.991, 5.986, 5.408, 5.452, 5.189, 5.761, 4.727, 4.737,
		    6.102, 6.508, 5.537, 5.582, 7.934, 6.481, 5.990, 5.286};

  Double_t da[16] = {9.6220, 8.5490, 13.163, 16.494, 17.426, 10.404, 17.290, 15.457,
  		     15.055, 13.412, 11.509, 10.308, 12.092, 10.445, 20.887, 32.637};
  Double_t db[16] = {0.123, 0.111, 0.145, 0.168, 0.156, 0.116, 0.151, 0.141,
		     0.183, 0.187, 0.134, 0.115, 0.212, 0.139, 0.249, 0.304};		  

  TH2F *gains = new TH2F("ADPMGains", "AD PM gain factors", 16, -0.5, 15.5, 2, -0.5, 1.5);
  for(Int_t channel = 0; channel < 16; ++channel) {
    gains->SetBinContent(channel+1,1,a[channel]);
    gains->SetBinContent(channel+1,2,b[channel]);
    gains->SetBinError(channel+1,1,da[channel]);
    gains->SetBinError(channel+1,2,db[channel]);
  }
	
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Michal Broz");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("PM gain factors channel by channel");
  md->AddDateToComment();
  md->PrintMetaData();

  AliCDBId id("AD/Calib/PMGains", 0, AliCDBRunRange::Infinity());
  man->Put(gains, id, md);

  delete md;

}
