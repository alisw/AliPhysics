CheckTOFFEElight(Char_t *fileName)
{

  AliTOFFEEReader r;
  r.LoadFEElightConfig(fileName);
  Int_t nch = r.ParseFEElightConfig();
  printf("found %d channels enabled\n", nch);

  TH1F *hDO = new TH1F("hDO", "detector-oriented (DO);index;enabled", 157248, 0, 157248);
  TH1F *hEO = new TH1F("hEO", "electronics-oriented (EO);index;enabled", 172800, 0, 172800);

  AliTOFcalibHisto ch;
  ch.LoadCalibHisto();

  Int_t drm, trm, chain, tdc, channel, indexEO;

  for (Int_t indexDO = 0; indexDO < 157248; indexDO++)
    if (r.IsChannelEnabled(indexDO)) {
      drm = (Int_t)ch.GetCalibMap(AliTOFcalibHisto::kDDL, indexDO);
      trm = (Int_t)ch.GetCalibMap(AliTOFcalibHisto::kTRM, indexDO);
      chain = (Int_t)ch.GetCalibMap(AliTOFcalibHisto::kChain, indexDO);
      tdc = (Int_t)ch.GetCalibMap(AliTOFcalibHisto::kTDC, indexDO);
      channel = (Int_t)ch.GetCalibMap(AliTOFcalibHisto::kChannel, indexDO);
      indexEO = (Int_t)ch.GetIndexEO(drm, trm, chain, tdc, channel);
      hDO->Fill(indexDO);
      hEO->Fill(indexEO);
    }

  TCanvas *c = new TCanvas("c");
  c->Divide(1,2);
  c->cd(1);
  hDO->Draw();


  c->cd(2);
  hEO->SetLineColor(2);
  hEO->SetLineWidth(2);
  hEO->Draw();
  for (Int_t i = 0; i < 720; i++) {
    TLine *l = new TLine(i * 240, 0., i * 240, 0.25);
    l->Draw("same");
  }
  for (Int_t i = 0; i < 72; i++) {
    TLine *l = new TLine(i * 2400, 0., i * 2400, 0.5);
    l->SetLineColor(4);
    l->SetLineWidth(2);
    l->Draw("same");
  }
  
    
}
