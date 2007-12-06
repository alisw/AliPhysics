void TestPreprocessor()
{
  gSystem->Load("libT0shuttle.so");
  gSystem->Load("$ALICE_ROOT/SHUTTLE/TestShuttle/libTestShuttle.so");
  gSystem->Load("libSpectrum");


  AliTestShuttle::SetMainCDB("local://./TestCDB");
  AliTestShuttle::SetLocalCDB("local://./TestCDB");

  AliTestShuttle::SetMainRefStorage("local://./TestRef");
  AliTestShuttle::SetLocalRefStorage("local://./TestRef");

  AliTestShuttle* shuttle = new AliTestShuttle(0, 0, 1);

  TMap* dcsAliasMap = CreateDCSAliasMap();

  shuttle->SetDCSInput(dcsAliasMap);

  shuttle->SetInputRunType("T0_STANDALONE_LASER");

  shuttle->AddInputFile(AliTestShuttle::kDAQ, "T00", "LASER", "LDC0","daLaser.root");
 
  //shuttle->SetInputRunType("PHYSICS");

  shuttle->AddInputFile(AliTestShuttle::kDAQ, "T00", "PHYSICS", "LDC0", "daPhys.root");

  AliPreprocessor* start = new AliT0Preprocessor(shuttle);

  shuttle->Process();
  
  /* AliCDBManager::Instance()->SetDefaultStorage("local://TestCDB");

  AliCDBEntry* entry = AliCDBManager::Instance()->Get("T00/Calib/Data", 0);
  if (!entry)
  {
    printf("The file is not there. Something went wrong.\n");
    return;
  }

  AliT0Calc* output = dynamic_cast<AliT0Calc*> (entry->GetObject());

   // output->Print();
  */
}

TMap* CreateDCSAliasMap()
{
  TMap* aliasMap = new TMap;
  aliasMap->SetOwner(1);

  Int_t n_ac_scalers=32;
  
  for(int nAlias=0;nAlias<n_ac_scalers;nAlias++)
  {
    TObjArray* valueSet = new TObjArray;
    valueSet->SetOwner(1);

    TString aliasName=Form("t00_ac_scaler_%02d",nAlias);
    
    Int_t nValues=10;

    for (int timeStamp=0;timeStamp<nValues;timeStamp++)
    {
      AliDCSValue* dcsVal = new AliDCSValue((Float_t) gRandom->Gaus(3.0e8,50), timeStamp);
      valueSet->Add(dcsVal);
      printf("Alias: %s - value n. %d: (val=%d timestamp=%d)\n" ,
    	    aliasName.Data(), timeStamp, dcsVal->GetFloat(), dcsVal->GetTimeStamp());
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
