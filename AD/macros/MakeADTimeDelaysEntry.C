void MakeADTimeDelaysEntry(const char *outputCDB = "local://$ALICE_ROOT/../AliRoot/OCDB")
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(outputCDB);

  // Creation of the time delays OCDB object
  //const Double_t timeShift[18] = {0.0, 203.2, 203.4, 203.5, 203.0, 203.4, 203.5, 203.1, 203.2, 194.2, 194.4, 194.5, 194.2, 194.7, 194.5, 194.3, 192.8, 0.0};
  const Double_t timeShift[18] = {0.0, 61.6091, 61.1891, 60.5191, 61.3591, 60.7691, 62.0291, 61.1091, 61.4591, 62.3491, 62.7891, 59.7791, 60.0991, 63.3091, 62.7691, 59.6491, 61.5091, 0.0};
  TH1F *delays = new TH1F("ADTimeDelays", "AD Time delays", 16, -0.5, 15.5);
  delays->SetContent(timeShift);
	
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Michal Broz");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Time delays channel by channel");
  md->PrintMetaData();

  AliCDBStorage *storLoc = man->GetDefaultStorage();
  AliCDBId id("AD/Calib/TimeDelays", 0, AliCDBRunRange::Infinity());

  man->Put(delays, id, md);

  delete md;

}
