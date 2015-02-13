void MakeADLightYieldsEntry(const char *outputCDB = "local://$ALICE_ROOT/../AliRoot/OCDB")
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(outputCDB);

  // Creation of the light yields OCDB object
  const Double_t lightYieldCorr[66] = {0.0,
				       0.01051, 0.00955, 0.00861, 0.00948, 0.01082, 0.00870, 0.01023, 0.01012,
				       0.01270, 0.01184, 0.01110, 0.01266, 0.00956, 0.00826, 0.00966, 0.00891,
				       0.0};

  TH1F *yields = new TH1F("ADLightYields", "AD Light Yields", 16, -0.5, 15.5);
  yields->SetContent(lightYieldCorr);
	
  AliCDBMetaData *md = new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Michal Broz");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Light Yields channel by channel");
  md->PrintMetaData();

  AliCDBId id("AD/Calib/LightYields", 0, AliCDBRunRange::Infinity());
  man->Put(yields, id, md);

  delete md;

}
