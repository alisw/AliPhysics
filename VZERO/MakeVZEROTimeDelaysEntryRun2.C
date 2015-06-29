void MakeVZEROTimeDelaysEntryRun2()
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://./OCDB");

  // Creation of the time delays OCDB object

  const Double_t timeShift[66] = {0.0 , 263.453366 , 263.554201 , 263.455896 , 263.705908 , 263.178068 , 263.260511 , 263.352692 , 263.371615 , 264.275609 , 263.022590 , 263.923256 , 263.477724 , 263.328535 , 263.338143 , 263.684347 , 263.636735 , 264.240256 , 264.704787 , 263.341181 , 265.875077 , 264.362158 , 264.009584 , 263.985062 , 264.507631 , 264.711630 , 264.702209 , 264.983539 , 265.156149 , 265.340929 , 265.185957 , 265.402229 , 267.060006 , 260.193899 , 260.831277 , 260.870380 , 260.231921 , 259.848971 , 263.069287 , 262.829099 , 261.297212 , 260.468547 , 260.962657 , 260.754787 , 260.782074 , 260.244392 , 263.248285 , 262.224661 , 259.936356 , 262.307604 , 262.698634 , 262.259535 , 262.425491 , 262.041006 , 264.711811 , 264.537483 , 262.019158 , 263.940932 , 263.309260 , 263.819921 , 264.324985 , 263.419804 , 265.954522 , 266.150658 , 263.942502 , 0.0};
  TH1F *delays = new TH1F("VZEROTimeDelays","VZERO Time delays",64,-0.5,63.5);
  delays->SetContent(timeShift);
	
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Time delays channel by channel for Run2");
  md->PrintMetaData();

  AliCDBId id("VZERO/Calib/TimeDelays",215011,AliCDBRunRange::Infinity());

  man->Put(delays, id, md);

  delete md;

}
