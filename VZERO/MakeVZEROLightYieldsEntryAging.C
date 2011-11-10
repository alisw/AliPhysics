void MakeVZEROLightYieldsEntryAging()
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  // Creation of the light yields OCDB object
  Double_t lightYieldCorr[66] = {0.0,
				       0.01051, 0.00955, 0.00861, 0.00948, 0.01082, 0.00870, 0.01023, 0.01012,
				       0.01270, 0.01184, 0.01110, 0.01266, 0.00956, 0.00826, 0.00966, 0.00891,
				       0.01358, 0.01543, 0.01516, 0.01337, 0.01908, 0.01641, 0.01767, 0.01512,
				       0.01664, 0.01326, 0.01536, 0.01500, 0.01439, 0.01445, 0.01504, 0.01079,
				       0.00105, 0.00110, 0.00143, 0.00093, 0.00072, 0.01919, 0.00073, 0.02580,
				       0.02911, 0.00148, 0.03176, 0.00126, 0.00158, 0.00111, 0.02804, 0.00109,
				       0.00157, 0.00158, 0.00104, 0.00120, 0.00123, 0.00188, 0.00193, 0.02133,
				       0.00200, 0.00185, 0.00143, 0.00257, 0.00201, 0.00119, 0.00197, 0.00282,
				       0.0};

  Double_t aging[66]={0.0, 0.74232, 0.86317, 0.93388, 0.77380, 0.71052, 0.63083, 0.50076, 0.46546, 0.72292, 0.82868, 0.83174, 0.86246, 0.60871, 0.54256, 0.60856, 0.55087, 0.73893, 0.88906, 0.71472, 0.59136, 0.53206, 0.61151, 0.62547, 0.64436, 0.74800, 0.83493, 0.77828, 0.68809, 0.75009, 1.00658, 0.71355, 0.69793, 0.59247, 0.56017, 0.52624, 0.73453, 0.71224, 0.48992, 0.37814, 0.61555, 0.37054, 0.82003, 0.69987, 0.84906, 0.68629, 0.63420, 0.61967, 0.49794, 0.44471, 0.82987, 0.87021, 0.70299, 0.53928, 0.70221, 0.65760, 0.81541, 1.05422, 0.80957, 0.62292, 0.80441, 0.88877, 0.80181, 0.74406, 0.73023, 0.0};
  for(Int_t i = 0; i < 66; ++i) lightYieldCorr[i] *= aging[i];

  TH1F *yields = new TH1F("VZEROLightYields","VZERO Light Yields",64,-0.5,63.5);
  yields->SetContent(lightYieldCorr);
	
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Light Yields channel by channel (corrected for aging in 2011)");
  md->PrintMetaData();

  AliCDBStorage *storLoc = man->GetDefaultStorage();
  AliCDBId id("VZERO/Calib/LightYields",0,AliCDBRunRange::Infinity());

  storLoc->Put(yields, id, md);

  storLoc->Delete();
  delete md;

}
