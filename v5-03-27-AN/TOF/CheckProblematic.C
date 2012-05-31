CheckProblematic(Int_t run, const Char_t *dbString = "raw://")
{

  TGrid::Connect("alien");
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(dbString);
  cdb->SetRun(run);
  AliCDBEntry *cdbe = cdb->Get("TOF/Calib/Problematic");
  CheckProblematic(cdbe);

}

CheckProblematic(const Char_t *fileName)
{

  TFile *file = TFile::Open(fileName);
  AliCDBEntry *cdbe = (AliCDBEntry *)file->Get("AliCDBEntry");
  CheckProblematic(cdbe);
}

CheckProblematic(AliCDBEntry *cdbe)
{

  if (!cdbe) {
    printf("invalid CDB entry\n");
    return;
  }

  TH1C *data = (TH1C *)cdbe->GetObject();

  TH2F *hProblematicMap = new TH2F("hProblematicMap", "Problematic map;sector;strip", 72, 0., 18., 91, 0., 91.);

  AliTOFcalibHisto calib;
  calib.LoadCalibHisto();
  calib.LoadCalibStat(); /* temp */

  Int_t sector, sectorStrip, padx, fea;
  Float_t efficiency, hitmapx, hitmapy;
  for (Int_t i = 0; i <  data->GetNbinsX(); i++) {
    if (data->GetBinContent(i + 1) == 0) continue;
    sector = calib.GetCalibMap(AliTOFcalibHisto::kSector, i);
    sectorStrip = calib.GetCalibMap(AliTOFcalibHisto::kSectorStrip, i);
    padx = calib.GetCalibMap(AliTOFcalibHisto::kPadX, i);
    fea = padx / 12;
    hitmapx = sector + ((Double_t)(3 - fea) + 0.5) / 4.;
    hitmapy = sectorStrip;

    hProblematicMap->Fill(hitmapx, hitmapy);
  }

  hProblematicMap->DrawCopy("colz");
  TFile *fout = TFile::Open("CheckProblematic.root", "RECREATE");
  hProblematicMap->Write();
  fout->Close();

}
