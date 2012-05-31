void TestPreprocessor()
{
  gSystem->Load("libT0shuttle.so");
  gSystem->Load("$ALICE_ROOT/SHUTTLE/TestShuttle/libTestShuttle.so");
  gSystem->Load("libSpectrum");


  AliTestShuttle::SetMainCDB("local://./TestCDB");
  AliTestShuttle::SetLocalCDB("local://./TestCDB");

  AliTestShuttle::SetMainRefStorage("local://./TestRef");
  AliTestShuttle::SetLocalRefStorage("local://./TestRef");

  AliTestShuttle* shuttle = new AliTestShuttle(104890, 0, 1);

 TMap* dcsAliasMap = CreateDCSAliasMap();



  shuttle->SetDCSInput(dcsAliasMap);

 shuttle->SetInputRunType("AMPLITUDE_CALIBRATION");
 
 //shuttle->AddInputFile(AliTestShuttle::kDAQ, "T00", "AMPLITUDE_CALIBRATION", "LDC0","daLaser.root");
   shuttle->SetInputRunType("PHYSICS");

  shuttle->AddInputFile(AliTestShuttle::kDAQ, "T00", "PHYSICS", "LDC0", "daPhys.root");

  AliPreprocessor* start = new AliT0Preprocessor(shuttle);

  shuttle->Process();
}

TMap* CreateDCSAliasMap()
{
  TMap* aliasMap = new TMap;
  aliasMap->SetOwner(1);

  TString aliasName;
  Int_t n_T0aliases=201;
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
      aliasName=Form("t00_ac_scaler_sec_%02d",nAlias-32);
    }
    else if(nAlias < 76)
    {
      aliasName=Form("t00_a_hv_imon_%02d",nAlias-64);
    }
    else if(nAlias < 88)
    {
      aliasName=Form("t00_a_hv_vmon_%02d",nAlias-76);
    }
    else if(nAlias < 90)
    {
      aliasName=Form("t00_a_lv_imon_%01d",nAlias-88);
    }
    else if(nAlias < 92)
    {
      aliasName=Form("t00_a_lv_vmon_%01d",nAlias-90);
    }
    else if(nAlias < 104)
    {
      aliasName=Form("t00_c_hv_imon_%02d",nAlias-92);
    }
    else if(nAlias < 116)
    {
      aliasName=Form("t00_c_hv_vmon_%02d",nAlias-104);
    }
    else if(nAlias < 118)
    {
      aliasName=Form("t00_c_lv_imon_%01d",nAlias-116);
    }
    else if(nAlias < 120)
    {
      aliasName=Form("t00_c_lv_vmon_%01d",nAlias-118);
    }
    else if(nAlias < 132)
    {
      aliasName=Form("t00_a_cfd_thre_%02d",nAlias-120);
    }
    else if(nAlias < 144)
    {
      aliasName=Form("t00_a_cfd_walk_%02d",nAlias-132);
    }
    else if(nAlias < 156)
    {
      aliasName=Form("t00_c_cfd_thre_%02d",nAlias-144);
    }
    else if(nAlias < 168)
    {
      aliasName=Form("t00_c_cfd_walk_%02d",nAlias-156);
    }
    else if(nAlias < 188)
    {
      aliasName=Form("t00_ac_trm_%02d",nAlias-168);
    }
    else if(nAlias < 193)
    {
      aliasName=Form("t00_ac_drm_%02d",nAlias-188);
    }
    else if (nAlias < 194)
    {
      aliasName=Form("t00_ac_atten");
    }
    else if(nAlias < 195)
    {
      aliasName=Form("t00_a_mpd_cent");
    }
    else if(nAlias < 196)
    {
      aliasName=Form("t00_c_mpd_cent");
    }
    else if(nAlias < 197)
    {
      aliasName=Form("t00_a_mpd_scent");
    }
    else if(nAlias < 198)
    {
      aliasName=Form("t00_c_mpd_scent");
    }
    else if(nAlias < 199)
    {
      aliasName=Form("t00_ac_tvdc_top");
    }
    else if(nAlias < 200)
    {
      aliasName=Form("t00_ac_tvdc_bottom");
    }
    else
    {
      aliasName=Form("t00_ac_mpd_mode");		
    }

    for (int timeStamp=0;timeStamp<nValues;timeStamp++)
    {
      //CHIARA's original // AliDCSValue* dcsVal = new AliDCSValue((Float_t) gRandom->Gaus(3.0e8,50), timeStamp);
      AliDCSValue* dcsVal = new AliDCSValue((Float_t) gRandom->Gaus(3.0e3,50), timeStamp);
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
	metaData.SetResponsible("Tomek");
	metaData.SetComment("Test object for TestPreprocessor.C");

  AliCDBId id("DET/DCS/Data", 0, 0);

  // initialize location of CDB
  AliCDBManager::Instance()->SetDefaultStorage("local://./TestCDB");

  AliCDBManager::Instance()->Put(dcsAliasMap, id, &metaData);
}
