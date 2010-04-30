MakeCDBEntryRunParams(Int_t startRun = 0, Int_t endRun = AliCDBRunRange::Infinity())
{

  const Int_t nPoints = 1;
  UInt_t timestamp[nPoints] = {0.};
  Float_t t0[nPoints] = {0.};
  Float_t tofReso[nPoints] = {80.};
  Float_t t0Spread[nPoints] = {0.};

  /* create object */
  AliTOFRunParams *obj = new AliTOFRunParams(1);
  obj->SetTimestamp(timestamp);
  obj->SetT0(t0);
  obj->SetTOFResolution(tofReso);
  obj->SetT0Spread(t0Spread);

  /* create cdb info */
  AliCDBId id("TOF/Calib/RunParams", startRun, endRun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Roberto Preghenella");
  md->SetComment("RunParams");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);

  /* put object in cdb */
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->GetDefaultStorage()->Put(obj, id, md);

}
