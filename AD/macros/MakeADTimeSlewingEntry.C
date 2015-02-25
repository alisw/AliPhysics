void MakeADTimeSlewingEntry(const char *outputCDB = "local://$ALICE_ROOT/../AliRoot/OCDB")
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(outputCDB);

  // Creation of the time slewing OCDB object
  TF1 *slew = new TF1("TimeSlewing","[0]*TMath::Power(x,[1])",1,1024);
  slew->SetParameter(0,286.8);
  slew->SetParameter(1,-1.437);
	
  TObjString str("AD Time-slewing correction");

  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Michal Broz");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Time-slewing correction used in reconstruction and MC simulation");
  md->PrintMetaData();

  AliCDBStorage *storLoc = man->GetDefaultStorage();
  AliCDBId id("AD/Calib/TimeSlewing",0,AliCDBRunRange::Infinity());

  storLoc->Put(slew, id, md);

  storLoc->Delete();
  delete md;

}
