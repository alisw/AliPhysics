void MakeVZEROPMGainsEntry()
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  // Creation of the PM gains OCDB object
  Double_t a[64] = {-39.68,-35.83,-36.92,-36.42,-37.02,-37.50,-43.05,-39.39,
		    -36.62,-36.93,-37.30,-36.46,-39.51,-40.32,-39.92,-39.20,
		    -35.39,-37.95,-38.85,-42.76,-40.68,-40.32,-39.00,-37.36,
		    -39.64,-38.86,-37.59,-39.90,-37.97,-36.32,-38.88,-41.35,
		    -36.01,-36.82,-39.48,-36.86,-38.22,-32.55,-39.44,-35.08,
		    -29.91,-37.88,-33.25,-36.49,-37.25,-35.89,-40.31,-39.15,
		    -41.71,-37.07,-38.94,-36.04,-36.62,-32.96,-36.99,-30.71,
		    -36.66,-37.23,-35.98,-36.56,-35.64,-36.97,-35.88,-38.78};
  Double_t b[66] = {7.40,  6.83,  7.02,  6.94,  7.03,  7.04,  7.79,  7.27,
		    6.92,  6.96,  7.01,  6.90,  7.28,  7.38,  7.33,  7.23,
		    6.71,  7.05,  7.17,  7.69,  7.41,  7.38,  7.21,  7.11,
		    7.26,  7.12,  6.98,  7.35,  6.99,  6.79,  7.13,  7.58,
		    6.95,  7.01,  7.33,  7.01,  7.21,  6.01,  7.34,  6.44,
		    5.68,  7.12,  6.07,  6.92,  7.04,  6.82,  7.04,  7.24,
		    7.53,  6.99,  7.10,  6.89,  7.07,  6.35,  6.88,  5.77,
		    6.81,  7.01,  6.89,  6.84,  6.68,  6.95,  6.73,  7.14};

  TH2F *gains = new TH2F("VZEROPMGains","VZERO PM gain factors",64,-0.5,63.5,2,-0.5,1.5);
  for(Int_t channel = 0; channel < 64; ++channel) {
    gains->SetBinContent(channel+1,1,a[channel]);
    gains->SetBinContent(channel+1,2,b[channel]);
  }
	
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("PM gain factors channel by channel");
  md->PrintMetaData();

  AliCDBStorage *storLoc = man->GetDefaultStorage();
  AliCDBId id("VZERO/Calib/PMGains",0,AliCDBRunRange::Infinity());

  storLoc->Put(gains, id, md);

  storLoc->Delete();
  delete md;

}
