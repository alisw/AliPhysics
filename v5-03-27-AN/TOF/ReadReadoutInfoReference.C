AliTOFReadoutInfo *
ReadReadoutInfoReference(Int_t run, Int_t year = 2010)
{

  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(Form("alien://folder=/alice/data/%d/Reference", year));
  cdb->SetRun(run);
  AliCDBEntry *cdbe = cdb->Get("TOF/Calib/ReadoutInfo");
  AliTOFReadoutInfo *ri = (AliTOFReadoutInfo *)cdbe->GetObject();
  
  return ri;
}
