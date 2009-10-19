CheckCalibStatus(const Char_t *fileName)
{

  TFile *file = TFile::Open(fileName);
  AliCDBEntry *cdbe = (AliCDBEntry *)file->Get("AliCDBEntry");
  AliTOFChannelOnlineStatusArray *array = (AliTOFChannelOnlineStatusArray *)cdbe->GetObject();

  TH1F *hStatus = new TH1F("hStatus", "Channel status;index;status", array->GetSize(), 0., array->GetSize(););
  TH2F *hNoiseMap = new TH2F("hNoiseMap", "Noise map;sector;strip", 72, 0., 18., 91, 0., 91.);
  TH2F *hEnableMap = new TH2F("hEnableMap", "Enable map;sector;strip", 72, 0., 18., 91, 0., 91.);

  AliTOFcalibHisto calib;
  calib.LoadCalibHisto();

  Int_t sector, sectorStrip, padx, fea;
  Float_t hitmapx, hitmapy;
  for (Int_t i = 0; i <  array->GetSize(); i++) {
    hStatus->SetBinContent(i + 1, array->GetStatus(i));
    sector = calib.GetCalibMap(AliTOFcalibHisto::kSector, i);
    sectorStrip = calib.GetCalibMap(AliTOFcalibHisto::kSectorStrip, i);
    padx = calib.GetCalibMap(AliTOFcalibHisto::kPadX, i);
    fea = padx / 12;
    hitmapx = sector + ((Double_t)(3 - fea) + 0.5) / 4.;
    hitmapy = sectorStrip;
    if (array->GetHWStatus(i) == AliTOFChannelOnlineStatusArray::kTOFHWOk) hEnableMap->Fill(hitmapx, hitmapy);
    if (array->GetNoiseStatus(i) == AliTOFChannelOnlineStatusArray::kTOFNoiseBad) hNoiseMap->Fill(hitmapx, hitmapy);
  }

  TFile *fout = TFile::Open("CheckCalibStatus.root", "RECREATE");
  hStatus->Write();
  hNoiseMap->Write();
  hEnableMap->Write();
  fout->Close();

}
