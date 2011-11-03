void MakeVZEROLightYieldsEntry()
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  // Creation of the light yields OCDB object
  const Double_t lightYieldCorr[66] = {0.0,
				       0.01051, 0.00955, 0.00861, 0.00948, 0.01082, 0.00870, 0.01023, 0.01012,
				       0.01270, 0.01184, 0.01110, 0.01266, 0.00956, 0.00826, 0.00966, 0.00891,
				       0.01358, 0.01543, 0.01516, 0.01337, 0.01908, 0.01641, 0.01767, 0.01512,
				       0.01664, 0.01326, 0.01536, 0.01500, 0.01439, 0.01445, 0.01504, 0.01079,
				       0.00105, 0.00110, 0.00143, 0.00093, 0.00072, 0.01919, 0.00073, 0.02580,
				       0.02911, 0.00148, 0.03176, 0.00126, 0.00158, 0.00111, 0.02804, 0.00109,
				       0.00157, 0.00158, 0.00104, 0.00120, 0.00123, 0.00188, 0.00193, 0.02133,
				       0.00200, 0.00185, 0.00143, 0.00257, 0.00201, 0.00119, 0.00197, 0.00282,
				       0.0};

  TH1F *yields = new TH1F("VZEROLightYields","VZERO Light Yields",64,-0.5,63.5);
  yields->SetContent(lightYieldCorr);
	
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Light Yields channel by channel");
  md->PrintMetaData();

  AliCDBStorage *storLoc = man->GetDefaultStorage();
  AliCDBId id("VZERO/Calib/LightYields",0,AliCDBRunRange::Infinity());

  storLoc->Put(yields, id, md);

  storLoc->Delete();
  delete md;

}
