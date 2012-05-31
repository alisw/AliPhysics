MakeCDBEntryCTPLatency(Float_t value = 0., Int_t startRun = 0, Int_t endRun = AliCDBRunRange::Infinity())
{

  /* create object */
  AliTOFCTPLatency *obj = new AliTOFCTPLatency();
  obj->SetCTPLatency(value);

  /* create cdb info */
  AliCDBId id("TOF/Calib/CTPLatency", startRun, endRun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Roberto Preghenella");
  md->SetComment("CTPLatency (ps)");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);

  /* put object in cdb */
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->GetDefaultStorage()->Put(obj, id, md);

}
