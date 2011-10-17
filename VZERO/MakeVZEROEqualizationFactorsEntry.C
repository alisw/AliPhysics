void MakeVZEROEqualizationFactorsEntry()
{

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  // Creation of the light yields OCDB object
  const Double_t alpha[66] = {0.0,
			      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
			      0.0};

  TH1F *eqFactors = new TH1F("VZEROEqualizationFactors","VZERO Equalization Factors for Pb-Pb",64,-0.5,63.5);
  eqFactors->SetContent(alpha);
	
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Default entry for VZERO Equalization Factors object");
  md->PrintMetaData();

  AliCDBStorage *storLoc = man->GetDefaultStorage();
  AliCDBId id("VZERO/Calib/EqualizationFactors",0,AliCDBRunRange::Infinity());

  storLoc->Put(eqFactors, id, md);

  storLoc->Delete();
  delete md;

}
