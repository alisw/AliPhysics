

void TestPreprocessor()
{
  gSystem->Load("libTestShuttle.so");

  AliCDBManager::Instance()->SetDefaultStorage("local://./TestCDB");

  AliTestShuttle* shuttle = new AliTestShuttle(0, 0, 1);

  TMap* dcsAliasMap = CreateDCSAliasMap();

  shuttle->SetDCSInput(dcsAliasMap);

  shuttle->AddInputFile(AliTestShuttle::kDAQ, "T00", "TIME", "LDC0", "DAQfile.root");

  AliPreprocessor* start = new AliSTARTPreprocessor("T00", shuttle);

  shuttle->Process();

  AliCDBEntry* entry = AliCDBManager::Instance()->Get("T00/Calib/Data", 0);
  if (!entry)
  {
    printf("The file is not there. Something went wrong.\n");
    return;
  }

  AliSTARTCalc* output = dynamic_cast<AliSTARTCalc*> (entry->GetObject());

   // output->Print();
}

TMap* CreateDCSAliasMap()
{
  TMap* aliasMap = new TMap;
  aliasMap->SetOwner(1);

  for(int nAlias=0;nAlias<24;nAlias++)
  {
    TObjArray* valueSet = new TObjArray;
    valueSet->SetOwner(1);

    TString aliasName="T0HV";
    aliasName += nAlias;

    for (int timeStamp=0;timeStamp<1;timeStamp++)
    {
      AliDCSValue* dcsVal = new AliDCSValue((Float_t) nAlias, timeStamp);
      valueSet->Add(dcsVal);
//    printf("hello! dcsVal= %d %d\n" ,dcsVal->GetFloat(), dcsVal->GetTimeStamp());
    }
    aliasMap->Add(new TObjString(aliasName), valueSet);
  }

  return aliasMap;
}

TMap* ReadDCSAliasMap()
{
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("DET/DCS/Data", 0);
  return dynamic_cast<TMap*> (entry->GetObject());
}

void WriteDCSAliasMap()
{
  TMap* dcsAliasMap = CreateDCSAliasMap();

  AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("Responsible person");
	metaData.SetComment("Test object for TestPreprocessor.C");

  AliCDBId id("DET/DCS/Data", 0, 0);

  // initialize location of CDB
  AliCDBManager::Instance()->SetDefaultStorage("local://./TestCDB");

  AliCDBManager::Instance()->Put(dcsAliasMap, id, &metaData);
}
