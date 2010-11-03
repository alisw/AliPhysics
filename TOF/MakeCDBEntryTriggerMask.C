MakeCDBEntryTriggerMask(Int_t startRun = 0, Int_t endRun = AliCDBRunRange::Infinity())
{

  UInt_t triggerMask[72];
  for (Int_t i = 0; i < 72; i++)
    triggerMask[i] = 0xffffff;

  /* create object */
  AliTOFTriggerMask *obj = new AliTOFTriggerMask();
  obj->SetTriggerMaskArray(triggerMask);

  /* create cdb info */
  AliCDBId id("TRIGGER/TOF/TriggerMask", startRun, endRun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Roberto Preghenella");
  md->SetComment("TOF Trigger Mask");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);

  /* put object in cdb */
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->GetDefaultStorage()->Put(obj, id, md);

}
