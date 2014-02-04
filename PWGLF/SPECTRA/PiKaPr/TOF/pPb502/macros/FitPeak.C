TObjArray *
FitPeak(TH2 *h, Float_t startSigma, Float_t nSigmaMin, Float_t nSigmaMax, const Char_t *name = "hsigma", TF1 *fitFunc = NULL)
{
  /*
   * fit peak
   */

  TH1 *hpy = NULL;
  TH1 *hout[2];
  TH1 *hmean = h->ProjectionX(Form("%s_mean", name));
  hmean->Reset();
  TH1 *hsigma = h->ProjectionX(Form("%s_sigma", name));
  hsigma->Reset();
  for (Int_t i = 0; i < h->GetNbinsX(); i++) {
    hpy = h->ProjectionY("hpy", i + 1, i + 1);
    if (hpy->Integral() <= 0.) {
      delete hpy;
      continue;
    }
    fitFunc = FitPeak(hpy, startSigma, nSigmaMin, nSigmaMax, fitFunc);
    hmean->SetBinContent(i + 1, fitFunc->GetParameter(1));
    hmean->SetBinError(i + 1, fitFunc->GetParError(1));
    hsigma->SetBinContent(i + 1, fitFunc->GetParameter(2));
    hsigma->SetBinError(i + 1, fitFunc->GetParError(2));
    delete hpy;
  }
  //  new TCanvas("cmean");
  //  hmean->DrawClone();
  //  new TCanvas("csigma");
  //  hsigma->DrawClone();
  TObjArray *oa = new TObjArray();
  oa->Add(hmean);
  oa->Add(hsigma);
  return oa;
}

TF1 *
FitPeak(TH1 *h, Float_t startSigma, Float_t nSigmaMin, Float_t nSigmaMax, TF1 *fitFunc = NULL)
{
  /*
   * fit peak
   */

  if (!fitFunc)
    fitFunc = (TF1 *)gROOT->GetFunction("gaus");

  Double_t fitCent = h->GetBinCenter(h->GetMaximumBin());
  Double_t fitMin = fitCent - nSigmaMin * startSigma;
  Double_t fitMax = fitCent + nSigmaMax * startSigma;
  if (fitMin < h->GetXaxis()->GetXmin()) fitMin = h->GetXaxis()->GetXmin();
  if (fitMax > h->GetXaxis()->GetXmax()) fitMax = h->GetXaxis()->GetXmax();
  fitFunc->SetParameter(1, fitCent);
  fitFunc->SetParameter(2, startSigma);
  Int_t fitres = h->Fit(fitFunc, "WWq0", "", fitMin, fitMax);
  if (fitres != 0) return NULL;
  /* refit with better range */
  for (Int_t i = 0; i < 3; i++) {
    fitCent = fitFunc->GetParameter(1);
    fitMin = fitCent - nSigmaMin * fitFunc->GetParameter(2);
    fitMax = fitCent + nSigmaMax * fitFunc->GetParameter(2);
    if (fitMin < h->GetXaxis()->GetXmin()) fitMin = h->GetXaxis()->GetXmin();
    if (fitMax > h->GetXaxis()->GetXmax()) fitMax = h->GetXaxis()->GetXmax();
    fitres = h->Fit(fitFunc, "q0", "", fitMin, fitMax);
    if (fitres != 0) return NULL;
  }
  return fitFunc;
}
