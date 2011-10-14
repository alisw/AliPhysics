CheckEnabledChannels(const Char_t *runlist)
{

  ifstream is(runlist);
  Char_t buf[4096];
  Int_t run[1024];
  Int_t nrun = 0;
  while(!is.eof()) {
    is.getline(buf, 4096);
    if (is.eof()) break;
    run[nrun] = atoi(buf);
    printf("added run number %d\n", run[nrun]);
    nrun++;
  }
  printf("%d runs added\n", nrun);
  is.close();

  TH1F *hActive = new TH1F("hActive", "active channels;run;fraction", nrun, 0, nrun);
  TH1F *hReadout = new TH1F("hReadout", "good readout;run;fraction", nrun, 0, nrun);
  for (Int_t irun = 0; irun < nrun; irun++) {
    hr = CheckEnabledChannels(run[irun], kTRUE);
    ha = CheckEnabledChannels(run[irun], kFALSE);
    hReadout->SetBinContent(irun + 1, hr->Integral());
    hActive->SetBinContent(irun + 1, ha->Integral());
    hReadout->GetXaxis()->SetBinLabel(irun + 1, Form("%d", run[irun]));
    delete hr; delete ha;
  }
  
  hReadout->SetMarkerStyle(20);
  hReadout->SetMarkerColor(4);
  hActive->SetMarkerStyle(25);
  hActive->SetMarkerColor(2);
  hReadout->Sumw2();
  hActive->Sumw2();
  hReadout->Divide(hReadout, hActive, 1., 1., "B");
  hActive->Scale(1. / 152928.);
  hReadout->SetMinimum(0.);
  hReadout->SetMaximum(1.);
  hReadout->Draw("E");
  hActive->Draw("E, same");
  TLegend *l = gPad->BuildLegend();
  l->SetFillStyle(0);

}

TH1F *
CheckEnabledChannels(Int_t run, Bool_t checkROEff = kTRUE, const Char_t *dbString = "raw://")
{

  /* init */
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(dbString);
  cdb->SetRun(run);
  AliTOFcalib calib;
  calib.Init();

  TH2F *hEnabledMap = new TH2F("hEnabledMap", "Enabled channel map;sector;strip", 72, 0., 18., 91, 0., 91.);
  TH1F *hEnabledFlag = new TH1F("hEnabledFlag", "Enabled channel flag;index;flag", 157248, 0., 157248.);

  AliTOFcalibHisto calibhisto;
  calibhisto.LoadCalibHisto();
  calibhisto.LoadCalibStat(); /* temp */

  Int_t sector, sectorStrip, padx, fea;
  Float_t hitmapx, hitmapy;
  /* loop over channels */
  for (Int_t ich = 0; ich < 157248; ich++) {
    if (!calib.IsChannelEnabled(ich, checkROEff)) continue;
    sector = calibhisto.GetCalibMap(AliTOFcalibHisto::kSector, ich);
    sectorStrip = calibhisto.GetCalibMap(AliTOFcalibHisto::kSectorStrip, ich);
    padx = calibhisto.GetCalibMap(AliTOFcalibHisto::kPadX, ich);
    fea = padx / 12;
    hitmapx = sector + ((Double_t)(3 - fea) + 0.5) / 4.;
    hitmapy = sectorStrip;
    hEnabledMap->Fill(hitmapx, hitmapy);
    hEnabledFlag->SetBinContent(ich + 1, 1);
  }
  
  hEnabledMap->DrawCopy("colz");
  return hEnabledFlag;

}
