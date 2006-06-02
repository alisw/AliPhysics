// This class runs the test preprocessor
// It uses AliTestShuttle to simulate a full Shuttle process

// The input data is created in the functions
//   CreateDCSAliasMap() creates input that would in the same way come from DCS
//   ReadDCSAliasMap() reads from a file
//   CreateInputFilesMap() creates a list of local files, that can be accessed by the shuttle

void TestPreprocessor()
{
  gSystem->Load("libTestShuttle");

  // initialize location of CDB
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://TestCDB");

  // TODO decide where to read the data
  TMap* dcsAliasMap = 0;
  //dcsAliasMap = CreateDCSAliasMap();
  dcsAliasMap = ReadDCSAliasMap();

  // create TMap of available files
  TMap* inputFiles = CreateInputFilesMap();

  // create AliTestShuttle instance
  AliTestShuttle* shuttle = new AliTestShuttle(inputFiles);
  TestShuttle(shuttle);

  // create preprocessor
  AliPreprocessor* pp = new AliTestPreprocessor("TPC", shuttle);

  // call preprocessor
  pp->Initialize(0, 0, 1);
  pp->Process(dcsAliasMap);
}

TMap* CreateDCSAliasMap()
{
  // fill fake DCS object

  TObjArray* valueSet1 = new TObjArray;


  valueSet1->Add();

  TObjArray* valueSet2 = new TObjArray;
  valueSet2->Add();

  TMap* aliasMap = new TMap;
  aliasMap->Add(new TObjString("DCSAlias1"), valueSet1);
  aliasMap->Add(new TObjString("DCSAlias2"), valueSet2);

  return aliasMap;
}

TMap* ReadDCSAliasMap()
{
  // open the file
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("TPC/DCS/Data",21200);
  return dynamic_cast<TMap*> (entry->GetObject());
}

TMap* CreateInputFilesMap()
{
  // create a list of files which will be available from the AliTestShuttle

  inputFile1 = new TMap;
  inputFile1->Add(new TObjString("GDC"), new TObjString("file1.root"));

  inputFile2 = new TMap;
  inputFile2->Add(new TObjString("LDC0"), new TObjString("file2a.root"));
  inputFile2->Add(new TObjString("LDC1"), new TObjString("file2b.root"));
  inputFile2->Add(new TObjString("LDC2"), new TObjString("file2c.root"));

  TMap* inputFiles = new TMap;
  inputFiles->Add(new TObjString("DAQ-TPC-PEDESTALS"), inputFile1);
  inputFiles->Add(new TObjString("DAQ-TPC-DRIFTVELOCITY"), inputFile2);

  return inputFiles;
}

void TestShuttle(AliShuttleInterface* shuttle)
{
  const char* file = shuttle->GetFile(AliShuttleInterface::kDAQ, "TPC", "PEDESTALS", "GDC");
  cout << "GetFile: " << file << endl;
  cout << "GetFileSources: " << endl;
  shuttle->GetFileSources(AliShuttleInterface::kDAQ, "TPC", "DRIFTVELOCITY")->Print();

  shuttle->Log("TPC", "Log test");

  //shuttle->Store();
}
