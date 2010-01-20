CheckCalibStatus(Int_t run)
{

  TGrid::Connect("alien");
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  cdb->SetRun(run);
  AliCDBEntry *cdbe = cdb->Get("TOF/Calib/Status");
  CheckCalibStatus(cdbe);

}

CheckCalibStatus(const Char_t *fileName)
{

  TFile *file = TFile::Open(fileName);
  AliCDBEntry *cdbe = (AliCDBEntry *)file->Get("AliCDBEntry");
  CheckCalibStatus(cdbe);
}

CheckCalibStatus(AliCDBEntry *cdbe)
{

  if (!cdbe) {
    printf("invalid CDB entry\n");
    return;
  }

  AliTOFChannelOnlineStatusArray *array = (AliTOFChannelOnlineStatusArray *)cdbe->GetObject();

  TH1F *hStatus = new TH1F("hStatus", "Channel status;index;status", array->GetSize(), 0., array->GetSize(););
  TH1F *hChEnabled = new TH1F("hChEnabled", "Channel enabled;index;enabled", array->GetSize(), 0., array->GetSize(););
  TH1F *hChNoisy = new TH1F("hChNoisy", "Channel noise flag;index;noise flag", array->GetSize(), 0., array->GetSize(););
  TH2F *hNoiseMap = new TH2F("hNoiseMap", "Noise map;sector;strip", 72, 0., 18., 91, 0., 91.);
  TH2F *hEnableMap = new TH2F("hEnableMap", "Enable map;sector;strip", 72, 0., 18., 91, 0., 91.);
  TH2F *hStatusMap = new TH2F("hStatusMap", "Status map;sector;strip", 72, 0., 18., 91, 0., 91.);

  AliTOFcalibHisto calib;
  calib.LoadCalibHisto();
  calib.LoadCalibStat(); /* temp */

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

    if (array->GetHWStatus(i) == AliTOFChannelOnlineStatusArray::kTOFHWOk) 
      hEnableMap->Fill(hitmapx, hitmapy);
    if (calib.GetCalibStat(AliTOFcalibHisto::kStripStat, i) == 1) 
      hStatusMap->Fill(hitmapx, hitmapy);
    if (array->GetNoiseStatus(i) == AliTOFChannelOnlineStatusArray::kTOFNoiseBad) {
      hChNoisy->SetBinContent(i + 1, 1);
      hNoiseMap->Fill(hitmapx, hitmapy);
    }

    if (array->GetHWStatus(i) == AliTOFChannelOnlineStatusArray::kTOFHWOk &&
	array->GetNoiseStatus(i) != AliTOFChannelOnlineStatusArray::kTOFNoiseBad &&
	calib.GetCalibStat(AliTOFcalibHisto::kStripStat, i) == 1) {
      hChEnabled->SetBinContent(i + 1, 1);
    }
   
  }

  TFile *fout = TFile::Open("CheckCalibStatus.root", "RECREATE");
  hStatus->Write();
  hChNoisy->Write();
  hChEnabled->Write();
  hNoiseMap->Write();
  hEnableMap->Write();
  hStatusMap->Write();
  fout->Close();

}
