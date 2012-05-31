AliTOFTriggerMask *
ReadCDBEntryTriggerMask(Int_t run, const Char_t *defaultStorage = "raw://", const Char_t *specificStorage = NULL)
{
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(defaultStorage);
  if (specificStorage)
    cdb->SetSpecificStorage("TRIGGER/TOF/TriggerMask", specificStorage);
  cdb->SetRun(run);
  AliCDBEntry *cdbe = cdb->Get("TRIGGER/TOF/TriggerMask");
  AliTOFTriggerMask *triggerMaskObj = (AliTOFTriggerMask *)cdbe->GetObject();

  return triggerMaskObj;
}
