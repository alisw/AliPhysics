CheckCalibStatus(const Char_t *fileName)
{

  TFile *file = TFile::Open(fileName);
  AliCDBEntry *cdbe = (AliCDBEntry *)file->Get("AliCDBEntry");
  AliTOFChannelOnlineStatusArray *array = (AliTOFChannelOnlineStatusArray *)cdbe->GetObject();

  TH1F *h = new TH1F("h", "Channel status;index;status", array->GetSize(), 0., array->GetSize(););
  for (Int_t i = 0; i <  array->GetSize(); i++)
    h->SetBinContent(i + 1, array->GetStatus(i));
  
  h->Draw();

}
