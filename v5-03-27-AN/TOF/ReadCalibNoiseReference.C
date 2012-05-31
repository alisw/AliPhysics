TH1F *
ReadCalibNoiseReference(Int_t run, Int_t year = 2010)
{

  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(Form("alien://folder=/alice/data/%d/Reference", year));
  cdb->SetRun(run);
  AliCDBEntry *cdbe = cdb->Get("TOF/Calib/CalibNoise");
  TH1F *h = (TH1F *)cdbe->GetObject();
  h->DrawCopy();
  return h;

}
