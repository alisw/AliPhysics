void MakeADTimeDelaysEntry()
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  // Creation of the time delays OCDB object
  const Double_t timeShift[18] = {0.0 ,203.2, 203.4, 203.5, 203.0, 203.4, 203.5, 203.1, 203.2, 194.2, 194.4, 194.5, 194.2, 194.7, 194.5, 194.3, 192.8, 0.0};
  TH1F *delays = new TH1F("ADTimeDelays","AD Time delays",16,-0.5,15.5);
  delays->SetContent(timeShift);
	
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Michal Broz");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Time delays channel by channel");
  md->PrintMetaData();

  AliCDBStorage *storLoc = man->GetDefaultStorage();
  AliCDBId id("AD/Calib/TimeDelays",0,AliCDBRunRange::Infinity());

  storLoc->Put(delays, id, md);

  storLoc->Delete();
  delete md;

}
