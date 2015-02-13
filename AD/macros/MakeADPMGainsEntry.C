void MakeADPMGainsEntry()
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  // Creation of the PM gains OCDB object
  Double_t a[16] = {1012, 1203, 1021, 1083, 835, 1067, 1163, 1050,
  		    1033, 1026, 1121, 1135, 1050, 1312, 957, 1013};
  Double_t b[16] = {7.32, 8.91, 6.05, 6.57, 5.34, 6.64, 8.25, 6.44,
		    7.24, 6.77, 6.85, 7.67, 4.86, 7.38, 5.48 ,7.13};

  TH2F *gains = new TH2F("ADPMGains","AD PM gain factors",16,-0.5,15.5,2,-0.5,1.5);
  for(Int_t channel = 0; channel < 16; ++channel) {
    gains->SetBinContent(channel+1,1,a[channel]);
    gains->SetBinContent(channel+1,2,b[channel]);
  }
	
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Michal Broz");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("PM gain factors channel by channel");
  md->PrintMetaData();

  AliCDBStorage *storLoc = man->GetDefaultStorage();
  AliCDBId id("AD/Calib/PMGains",0,AliCDBRunRange::Infinity());

  storLoc->Put(gains, id, md);

  storLoc->Delete();
  delete md;

}
