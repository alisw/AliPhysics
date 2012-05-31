TH1F *
ReadT0FillReference(Int_t run, Int_t year = 2010)
{

  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(Form("alien://folder=/alice/data/%d/Reference", year));
  cdb->SetRun(run);
  AliCDBEntry *cdbe = cdb->Get("TOF/Calib/T0Fill");
  TH1F *hT0Fill = (TH1F *)cdbe->GetObject();

  /* rebin until maximum bin has required minimum entries */
  Int_t maxBin = hT0Fill->GetMaximumBin();
  Float_t maxBinContent = hT0Fill->GetBinContent(maxBin);
  Float_t binWidth = hT0Fill->GetBinWidth(maxBin);
  while (maxBinContent < 400 && binWidth < 90.) {
    hT0Fill->Rebin(2);
    maxBin = hT0Fill->GetMaximumBin();
    maxBinContent = hT0Fill->GetBinContent(maxBin);
    binWidth = hT0Fill->GetBinWidth(maxBin);
  }
  Float_t maxBinCenter = hT0Fill->GetBinCenter(maxBin);

  /* rough fit of the edge */
  TF1 *gaus = (TF1 *)gROOT->GetFunction("gaus");
  gaus->SetParameter(1, maxBinCenter);
  Float_t fitMin = maxBinCenter - 100.; /* fit from 0.1 ns before max */
  Float_t fitMax = maxBinCenter + 100.; /* fit until 0.1 ns above max */
  hT0Fill->Fit("gaus", "q0", "", fitMin, fitMax);
  /* better fit of the edge */
  Float_t mean, sigma;
  for (Int_t istep = 0; istep < 10; istep++) {
    mean = gaus->GetParameter(1);
    sigma = gaus->GetParameter(2);
    fitMin = mean - 1. * sigma;
    fitMax = mean + 0.1 * sigma;
    hT0Fill->Fit("gaus", "q", "", fitMin, fitMax);
  }
  /* print params */
  mean = gaus->GetParameter(1);
  sigma = gaus->GetParameter(2);
  Float_t meane = gaus->GetParError(1);
  Float_t sigmae = gaus->GetParError(2);
  printf("edge fit: mean  = %f +- %f ps\n", mean, meane);
  printf("edge fit: sigma = %f +- %f ps\n", sigma, sigmae);

  hT0Fill->DrawCopy();
  gaus->Draw("same");

  return hT0Fill;

}
