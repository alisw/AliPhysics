MakeCDBEntryT0Fill(Float_t value = 0., Int_t startRun = 0, Int_t endRun = AliCDBRunRange::Infinity())
{

  /* create object */
  AliTOFT0Fill *obj = new AliTOFT0Fill();
  obj->SetT0Fill(value);

  /* create cdb info */
  AliCDBId id("TOF/Calib/T0Fill", startRun, endRun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Roberto Preghenella");
  md->SetComment("T0Fill (ps)");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);

  /* put object in cdb */
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->GetDefaultStorage()->Put(obj, id, md);

}
