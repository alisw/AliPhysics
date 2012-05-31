MakeCDBEntryT0FillOnlineCalib(Float_t offset = 0., Float_t coefficient = 0., Int_t startRun = 0, Int_t endRun = AliCDBRunRange::Infinity())
{

  /* create object */
  AliTOFT0FillOnlineCalib *obj = new AliTOFT0FillOnlineCalib();
  obj->SetOffset(offset);
  obj->SetCoefficient(coefficient);

  /* create cdb info */
  AliCDBId id("TOF/Calib/T0FillOnlineCalib", startRun, endRun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Roberto Preghenella");
  md->SetComment("Online T0-fill calibration parameters");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);

  /* put object in cdb */
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->GetDefaultStorage()->Put(obj, id, md);

}
