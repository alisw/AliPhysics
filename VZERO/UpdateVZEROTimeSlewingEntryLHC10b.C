void UpdateVZEROTimeSlewingEntryLHC10b()
{
  // Macro repairs the OCDB entry for the
  // run range 106031 - 116353 (from LHC10b)
  // where the discri thr was not stored correctly
  // in VZERO/Calib/Data.
  // It is just ascale factor of TMath::Power(2.5/4.0,-0.5);


  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  // Creation of the time slewing OCDB object
  // x = ADC charge / TDC threshold
  TF1 *slew = new TF1("TimeSlewing","[0]*TMath::Power(x,[1])",1,1024);
  slew->SetParameter(0,1.05e1*TMath::Power(2.5/4.0,-0.5));
  slew->SetParameter(1,-5.0e-1);

  TObjString str("VZERO Time-slewing correction (corrected for wrong VZERO/Calib/Data entry in runs 106031-116353)");

  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Time-slewing correction used in reconstruction and MC simulation");
  md->PrintMetaData();

  AliCDBStorage *storLoc = man->GetDefaultStorage();
  AliCDBId id("VZERO/Calib/TimeSlewing",106031,116353);

  storLoc->Put(slew, id, md);

  storLoc->Delete();
  delete md;

}
