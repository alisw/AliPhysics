void MakeADLightYieldsEntry(const char *outputCDB = "local://$ALICE_ROOT/../AliRoot/OCDB")
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(outputCDB);

  // Creation of the light yields OCDB object
  const Double_t lightYieldCorr[18] = {0.0,
				       2.2e-4,2.2e-4,2.2e-4,2.2e-4, 2.2e-4,2.2e-4,2.2e-4,2.2e-4,
                                       2.4e-4,2.4e-4,2.6e-4,2.6e-4, 2.4e-4,2.4e-4,2.6e-4,2.6e-4,
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
