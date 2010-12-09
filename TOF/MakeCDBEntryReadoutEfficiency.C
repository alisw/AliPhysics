MakeCDBEntryReadoutEfficiency(Int_t startRun = 0, Int_t endRun = AliCDBRunRange::Infinity())
{

  /* create object */
  TH1F *obj = new TH1F("hReadoutEfficiency", "", 157248, 0., 157248.);
  for (Int_t i = 0; i < 157248; i++)
    obj->SetBinContent(i + 1, 1.);

  /* create cdb info */
  AliCDBId id("TOF/Calib/ReadoutEfficiency", startRun, endRun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Roberto Preghenella");
  md->SetComment("ReadoutEfficiency");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);

  /* put object in cdb */
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->GetDefaultStorage()->Put(obj, id, md);

}
