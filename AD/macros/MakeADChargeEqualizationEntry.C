// -*- C++ -*-

void MakeADChargeEqualizationEntry(const char *outputCDB = "local://$ALICE_ROOT/../AliRoot/OCDB")
{
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage(outputCDB);

  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("C. Mayer");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("AD charge equalization factors");
  md->PrintMetaData();

  AliCDBId id("AD/Calib/ChargeEqualization", 0, AliCDBRunRange::Infinity());

  // make a default charge equalization object
  AliADChargeEqualization *eq = new AliADChargeEqualization;

  man->Put(eq, id, md);

  delete md;
}
