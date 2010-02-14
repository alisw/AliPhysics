MakeCDBEntryDeltaBCOffset(Int_t value = 0, Int_t startRun = 0, Int_t endRun = AliCDBRunRange::Infinity())
{

  /* create object */
  AliTOFDeltaBCOffset *obj = new AliTOFDeltaBCOffset();
  obj->SetDeltaBCOffset(value);

  /* create cdb info */
  AliCDBId id("TOF/Calib/DeltaBCOffset", startRun, endRun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Roberto Preghenella");
  md->SetComment("DeltaBCOffset (BC bins)");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);

  /* put object in cdb */
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->GetDefaultStorage()->Put(obj, id, md);

}
