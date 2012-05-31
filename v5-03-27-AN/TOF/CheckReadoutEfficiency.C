TH1F *
CheckReadoutEfficiency(Int_t run)
{

  TGrid::Connect("alien");
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  cdb->SetRun(run);
  AliCDBEntry *cdbe = cdb->Get("TOF/Calib/ReadoutEfficiency");
  CheckReadoutEfficiency(cdbe);

}

TH1F *
CheckReadoutEfficiency(const Char_t *fileName)
{

  TFile *file = TFile::Open(fileName);
  AliCDBEntry *cdbe = (AliCDBEntry *)file->Get("AliCDBEntry");
  CheckReadoutEfficiency(cdbe);
}

TH1F *
CheckReadoutEfficiency(AliCDBEntry *cdbe)
{

  if (!cdbe) {
    printf("invalid CDB entry\n");
    return;
  }

  TH1F *data = (TH1F *)cdbe->GetObject();

  TH2F *hEfficiencyMap = new TH2F("hEfficiencyMap", "Readout efficiency map;sector;strip", 72, 0., 18., 91, 0., 91.);
  TH1F *hEfficiencyFlag = new TH1F("hEfficiencyFlag", "Readout efficiency flag;index;flag", 157248, 0., 157248.);

  AliTOFcalibHisto calib;
  calib.LoadCalibHisto();
  calib.LoadCalibStat(); /* temp */

  Int_t sector, sectorStrip, padx, fea;
  Float_t efficiency, hitmapx, hitmapy;
  for (Int_t i = 0; i <  data->GetNbinsX(); i++) {
    efficiency = data->GetBinContent(i + 1);
    sector = calib.GetCalibMap(AliTOFcalibHisto::kSector, i);
    sectorStrip = calib.GetCalibMap(AliTOFcalibHisto::kSectorStrip, i);
    padx = calib.GetCalibMap(AliTOFcalibHisto::kPadX, i);
    fea = padx / 12;
    hitmapx = sector + ((Double_t)(3 - fea) + 0.5) / 4.;
    hitmapy = sectorStrip;
    hEfficiencyMap->Fill(hitmapx, hitmapy, efficiency / 24.);
    if (efficiency >= 0.95)
      hEfficiencyFlag->SetBinContent(i + 1, 1);
  }

  hEfficiencyMap->DrawCopy("colz");
  
  return hEfficiencyFlag;
}
