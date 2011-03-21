CheckEnabledChannels(Int_t run, const Char_t *dbString)
{

  /* init */
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(dbString);
  cdb->SetRun(run);
  AliTOFcalib calib;
  calib.Init();

  TH2F *hEnabledMap = new TH2F("hEnabledMap", "Enabled channel map;sector;strip", 72, 0., 18., 91, 0., 91.);

  AliTOFcalibHisto calibhisto;
  calibhisto.LoadCalibHisto();
  calibhisto.LoadCalibStat(); /* temp */

  Int_t sector, sectorStrip, padx, fea;
  Float_t hitmapx, hitmapy;
  /* loop over channels */
  for (Int_t ich = 0; ich < 157248; ich++) {
    if (!calib.IsChannelEnabled(ich, kTRUE, kTRUE)) continue;
    sector = calibhisto.GetCalibMap(AliTOFcalibHisto::kSector, ich);
    sectorStrip = calibhisto.GetCalibMap(AliTOFcalibHisto::kSectorStrip, ich);
    padx = calibhisto.GetCalibMap(AliTOFcalibHisto::kPadX, ich);
    fea = padx / 12;
    hitmapx = sector + ((Double_t)(3 - fea) + 0.5) / 4.;
    hitmapy = sectorStrip;
    hEnabledMap->Fill(hitmapx, hitmapy);
  }
  
  hEnabledMap->DrawCopy("colz");
  TFile *fileout = TFile::Open("CheckEnabledChannels.root", "RECREATE");
  hEnabledMap->Write();
  fileout->Close();

}
