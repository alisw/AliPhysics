TH1F *
ReadT0FillReference(Int_t run, Int_t year = 2010)
{

  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(Form("alien://folder=/alice/data/%d/Reference", year));
  cdb->SetRun(run);
  AliCDBEntry *cdbe = cdb->Get("TOF/Calib/T0Fill");
  TH1F *h = (TH1F *)cdbe->GetObject();
  h->DrawCopy();
  return h;

}
