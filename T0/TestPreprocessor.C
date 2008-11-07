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

  //shuttle->SetInputRunType("STANDALONE");

  shuttle->AddInputFile(AliTestShuttle::kDAQ, "T00", "LASER", "LDC0","daLaser.root");
 
  shuttle->SetInputRunType("PHYSICS");

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

  TString aliasName;
  Int_t n_T0aliases=184;
  Int_t nValues=10;	
  
  for(int nAlias=0;nAlias<n_T0aliases;nAlias++)
  {
    TObjArray* valueSet = new TObjArray;
    valueSet->SetOwner(1);
    if(nAlias < 32)
    {			
      aliasName=Form("t00_ac_scaler_%02d",nAlias);
    }
    else if(nAlias < 64)
    {
      aliasName=Form("t00_ac_scaler_sec_%02d",nAlias);
    }
    else if(nAlias < 76)
    {
      aliasName=Form("t00_a_hv_imon_%02d",nAlias);
    }
    else if(nAlias < 88)
    {
      aliasName=Form("t00_a_hv_vmon_%02d",nAlias);
    }
    else if(nAlias < 90)
    {
      aliasName=Form("t00_a_lv_imon_%02d",nAlias);
    }
    else if(nAlias < 92)
    {
      aliasName=Form("t00_a_lv_vmon_%02d",nAlias);
    }
    else if(nAlias < 104)
    {
      aliasName=Form("t00_c_hv_imon_%02d",nAlias);
    }
    else if(nAlias < 116)
    {
      aliasName=Form("t00_c_hv_vmon_%02d",nAlias);
    }
    else if(nAlias < 118)
    {
      aliasName=Form("t00_c_lv_imon_%02d",nAlias);
    }
    else if(nAlias < 120)
    {
      aliasName=Form("t00_c_lv_vmon_%02d",nAlias);
    }
    else if(nAlias < 132)
    {
      aliasName=Form("t00_a_cfd_thre_%02d",nAlias);
    }
    else if(nAlias < 144)
    {
      aliasName=Form("t00_a_cfd_walk_%02d",nAlias);
    }
    else if(nAlias < 156)
    {
      aliasName=Form("t00_c_cfd_thre_%02d",nAlias);
    }
    else if(nAlias < 168)
    {
      aliasName=Form("t00_c_cfd_walk_%02d",nAlias);
    }
    else if(nAlias < 178)
    {
      aliasName=Form("t00_ac_trm_%02d",nAlias);
    }
    else if(nAlias < 183)
    {
      aliasName=Form("t00_ac_drm_%02d",nAlias);
    }
    else
    {
      aliasName=Form("t00_ac_atten",nAlias);
    }

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
