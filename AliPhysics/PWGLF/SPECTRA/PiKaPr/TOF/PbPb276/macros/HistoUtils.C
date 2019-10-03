TH1 *
HistoUtils_weightedmean(TH1 *h1, TH1 *h2)
{

  TH1 *ho = (TH1 *)h1->Clone("ho");
  ho->Reset();
  Double_t val1, val2, w1, w2, mean, meane;
  for (Int_t i = 0; i < ho->GetNbinsX(); i++) {
    val1 = h1->GetBinContent(i + 1);
    w1 = 1. / (h1->GetBinError(i + 1) * h1->GetBinError(i + 1));
    val2 = h2->GetBinContent(i + 1);
    w2 = 1. / (h2->GetBinError(i + 1) * h2->GetBinError(i + 1));

    if (val1 == 0 && val2 == 0) continue;

    mean = (w1 * val1 + w2 * val2) / (w1 + w2);
    meane = TMath::Sqrt(1. / (w1 + w2));

    ho->SetBinContent(i + 1, mean);
    ho->SetBinError(i + 1, meane);
  }

  return ho;
}

//__________________________________________________________________

TH1 *
HistoUtils_smartdifference(TH1 *hnum, TH1 *hden)
{

  TH1 *hr = (TH1 *)hnum->Clone("hr");
  hr->Reset();
  Double_t ref;
  for (Int_t i = 0; i < hr->GetNbinsX(); i++) {
    if (hnum->GetBinError(i + 1) <= 0.) continue;
    ref = hden->Interpolate(hr->GetBinCenter(i + 1));
    if (ref <= 0.) continue;
    hr->SetBinContent(i + 1, (hnum->GetBinContent(i + 1) - ref) / ref);
    hr->SetBinError(i + 1, hnum->GetBinError(i + 1) / ref);
  }
  return hr;
}

//__________________________________________________________________

TH1 *
HistoUtils_smartratio(TH1 *hnum, TH1 *hden)
{

  TH1 *hr = (TH1 *)hnum->Clone("hr");
  hr->Reset();
  Double_t ref;
  for (Int_t i = 0; i < hr->GetNbinsX(); i++) {
    if (hnum->GetBinError(i + 1) <= 0.) continue;
    ref = hden->Interpolate(hr->GetBinCenter(i + 1));
    if (ref <= 0.) continue;
    hr->SetBinContent(i + 1, hnum->GetBinContent(i + 1) / ref);
    hr->SetBinError(i + 1, hnum->GetBinError(i + 1) / ref);
  }
  return hr;
}

//__________________________________________________________________

TH1 *
HistoUtils_smartratio(TH1 *hnum, TGraph *hden)
{

  TH1 *hr = (TH1 *)hnum->Clone("hr");
  hr->Reset();
  Double_t ref;
  for (Int_t i = 0; i < hr->GetNbinsX(); i++) {
    if (hnum->GetBinError(i + 1) <= 0.) continue;
    ref = hden->Eval(hr->GetBinCenter(i + 1));
    if (ref <= 0.) continue;
    hr->SetBinContent(i + 1, hnum->GetBinContent(i + 1) / ref);
    hr->SetBinError(i + 1, hnum->GetBinError(i + 1) / ref);
  }
  return hr;
}

//__________________________________________________________________

HistoUtils_drawthemall(const Char_t *filename, Int_t sleepms = 100)
{
  TFile *f1 = TFile::Open(filename);
  TList *l1 = f1->GetListOfKeys();
  TObject *o;
  Char_t *name;
  for (Int_t i = 0; i < l1->GetEntries(); i++) {
    name = l1->At(i)->GetName();
    o = f1->Get(name);
    o->Draw();
    gPad->Update();
    gSystem->Sleep(sleepms);
  }
  f1->Close();
}

//__________________________________________________________________

HistoUtils_autoratio(const Char_t *f1name, const Char_t *f2name, const Char_t *outname, const Char_t *title = NULL)
{

  TFile *f1 = TFile::Open(f1name);
  TFile *f2 = TFile::Open(f2name);
  TFile *fo = TFile::Open(outname, "RECREATE");
  TList *l1 = f1->GetListOfKeys();


  TH1D *h1, *h2, *hr;
  Char_t *name;
  for (Int_t i = 0; i < l1->GetEntries(); i++) {
    name = l1->At(i)->GetName();
    h1 = (TH1D *)f1->Get(name);
    h2 = (TH1D *)f2->Get(name);
    if (!h1 || !h2) continue;
    hr = new TH1D(*h1);
    hr->Divide(h2);
    if (title)
      hr->SetTitle(title);
    fo->cd();
    hr->Write();
  }

  delete hr;
  f1->Close();
  f2->Close();
  fo->Close();

}

//__________________________________________________________________

HistoUtils_autosystematics(const Char_t *f1name, const Char_t *f2name, const Char_t *outname, const Char_t *title = NULL)
{

  TFile *f1 = TFile::Open(f1name);
  TFile *f2 = TFile::Open(f2name);
  if (!f1 || !f1->IsOpen() || !f2 || !f2->IsOpen()) return;
  TFile *fo = TFile::Open(outname, "RECREATE");
  TList *l1 = f1->GetListOfKeys();


  TH1D *h1, *h2, *hd;
  Char_t *name;
  for (Int_t i = 0; i < l1->GetEntries(); i++) {
    name = l1->At(i)->GetName();
    h1 = (TH1D *)f1->Get(name);
    h2 = (TH1D *)f2->Get(name);
    if (!h1 || !h2) continue;

    hd = new TH1D(*h1);
    Double_t val1, val2, vald, vale1, vale2, valde;
    for (Int_t ii = 0; ii < hd->GetNbinsX(); ii++) {
      val1 = h1->GetBinContent(ii + 1);
      vale1 = h1->GetBinError(ii + 1);
      val2 = h2->GetBinContent(ii + 1);
      vale2 = h2->GetBinError(ii + 1);
      if (val2 == 0.) continue;
      vald = (val1 - val2) / val2;
      valde = TMath::Sqrt(TMath::Abs((vale1 * vale1 - vale2 * vale2))) / val2;
      hd->SetBinContent(ii + 1, vald);
      hd->SetBinError(ii + 1, valde);
    }

    if (title)
      hd->SetTitle(title);
    fo->cd();
    hd->Write();
    delete hd;
  }

  f1->Close();
  f2->Close();
  fo->Close();

}

//__________________________________________________________________

TH1 *
HistoUtils_ratio(const Char_t *f1name, const Char_t *f2name, const Char_t *h1name, const Char_t *h2name)
{

  if (!f2name) f2name = f1name;
  TFile *f1 = TFile::Open(f1name);
  TFile *f2 = TFile::Open(f2name);

  if (!h2name) h2name = h1name;
  TH1 *h1 = (TH1 *)f1->Get(h1name);
  TH1 *h2 = (TH1 *)f2->Get(h2name);

  TH1 *hr = h1->Clone("hr");
  hr->Sumw2();
  hr->Divide(h2);

  return hr;
  
}

//__________________________________________________________________

TH1 *
HistoUtils_systematics(const Char_t *f1name, const Char_t *f2name, const Char_t *h1name, const Char_t *h2name)
{

  if (!f2name) f2name = f1name;
  TFile *f1 = TFile::Open(f1name);
  TFile *f2 = TFile::Open(f2name);

  if (!h2name) h2name = h1name;
  TH1 *h1 = (TH1 *)f1->Get(h1name);
  TH1 *h2 = (TH1 *)f2->Get(h2name);

  TH1 *hd = h1->Clone("hd");
  Double_t val1, val2, vald, vale1, vale2, valde;
  for (Int_t i = 0; i < h1->GetNbinsX(); i++) {
    val1 = h1->GetBinContent(i + 1);
    vale1 = h1->GetBinError(i + 1);
    val2 = h2->GetBinContent(i + 1);
    vale2 = h2->GetBinError(i + 1);
    if (val2 == 0.) continue;
    vald = (val1 - val2) / val2;
    valde = TMath::Sqrt(TMath::Abs((vale1 * vale1 - vale2 * vale2))) / val2;
    hd->SetBinContent(i + 1, vald);
    hd->SetBinError(i + 1, valde);
  }

  return hd;
  
}

//__________________________________________________________________

HistoUtils_BinLogX(TH1 *h)
{
  TAxis *axis = h->GetXaxis();
  Int_t bins = axis->GetNbins();
  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];
  for (int i = 0; i <= bins; i++) new_bins[i] = TMath::Power(10, from + i * width);
  axis->Set(bins, new_bins);
  delete new_bins;
}

//__________________________________________________________________

HistoUtils_BinNormX(TH1 *h)
{
  TAxis *axis = h->GetXaxis();
  Int_t bins = axis->GetNbins();
  Double_t c, ec;
  Double_t w;
  for (Int_t i = 0; i < bins; i++) {
    c = h->GetBinContent(i + 1);
    ec = h->GetBinError(i + 1);
    w = axis->GetBinWidth(i + 1);
    c /= w;
    ec /= w;
    h->SetBinContent(i + 1, c);
    h->SetBinError(i + 1, ec);
  }
}

//__________________________________________________________________

TObjArray *
HistoUtils_FitPeak(TF1 *fitFunc, TH2 *h, Float_t startSigma, Float_t nSigmaMin, Float_t nSigmaMax, Int_t minIntegral = 100., const Char_t *basename = "hParam", Bool_t monitor = kFALSE)
{

  /* gaus function if not specified */
  if (!fitFunc)
    fitFunc = (TF1*)gROOT->GetFunction("gaus");

  /* prepare output histos */
  TObjArray *outArray = new TObjArray();
  Int_t npars = fitFunc->GetNpar();
  TH1D *hpx = h->ProjectionX("hpx");
  TH1D *hParam;
  for (Int_t ipar = 0; ipar < npars; ipar++) {
    hParam = new TH1D(*hpx);
    hParam->SetName(Form("%s_%d", basename, ipar));
    hParam->Reset();
    outArray->Add(hParam);
  }

  /* loop over x-bins */
  for (Int_t ibin = 0; ibin < hpx->GetNbinsX(); ibin++) {

    /* check integral */
    if (hpx->GetBinContent(ibin + 1) < minIntegral) continue;
    /* projection y */
    TH1D *hpy = h->ProjectionY("hpy", ibin + 1, ibin + 1);
    /* fit peak */
    if (HistoUtils_FitPeak(fitFunc, hpy, startSigma, nSigmaMin, nSigmaMax) != 0) {
      delete hpy;
      continue;
    }
    /* setup output histos */
    for (Int_t ipar = 0; ipar < npars; ipar++) {
      hParam = (TH1D *)outArray->At(ipar);
      hParam->SetBinContent(ibin + 1, fitFunc->GetParameter(ipar));
      hParam->SetBinError(ibin + 1, fitFunc->GetParError(ipar));
    }
    /* monitor */
    if (monitor) {
      hpy->SetMarkerStyle(20);
      hpy->SetMarkerColor(4);
      hpy->Draw("E1");
      fitFunc->Draw("same");
      gPad->Update();
      getchar();
    }
    /* delete */
    delete hpy;

  }
  /* delete */
  delete hpx;
  /* return output array */
  return outArray;
}

//__________________________________________________________________

Int_t
HistoUtils_FitPeak(TF1 *fitFunc, TH1 *h, Float_t startSigma, Float_t nSigmaMin, Float_t nSigmaMax)
{

  Double_t fitCent = h->GetBinCenter(h->GetMaximumBin());
  Double_t fitMin = fitCent - nSigmaMin * startSigma;
  Double_t fitMax = fitCent + nSigmaMax * startSigma;
  if (fitMin < h->GetXaxis()->GetXmin()) fitMin = h->GetXaxis()->GetXmin();
  if (fitMax > h->GetXaxis()->GetXmax()) fitMax = h->GetXaxis()->GetXmax();
  fitFunc->SetParameter(1, fitCent);
  fitFunc->SetParameter(2, startSigma);
  Int_t fitres = h->Fit(fitFunc, "WWq0", "", fitMin, fitMax);
  if (fitres != 0) return fitres;
  /* refit with better range */
  for (Int_t i = 0; i < 3; i++) {
    fitCent = fitFunc->GetParameter(1);
    fitMin = fitCent - nSigmaMin * fitFunc->GetParameter(2);
    fitMax = fitCent + nSigmaMax * fitFunc->GetParameter(2);
    if (fitMin < h->GetXaxis()->GetXmin()) fitMin = h->GetXaxis()->GetXmin();
    if (fitMax > h->GetXaxis()->GetXmax()) fitMax = h->GetXaxis()->GetXmax();
    fitres = h->Fit(fitFunc, "q0", "", fitMin, fitMax);
    if (fitres != 0) return fitres;
  }
  return fitres;

}

//__________________________________________________________________

void
HistoUtils_Function2Profile(TF1 *fin, TProfile *p, Int_t ntry = 10000)
{

  TF1 *f = new TF1(*fin);
  Int_t npars = f->GetNpar();
  Double_t par[1000];
  Double_t pare[1000];
  for (Int_t ipar = 0; ipar < npars; ipar++) {
    par[ipar] = f->GetParameter(ipar);
    pare[ipar] = f->GetParError(ipar);
  }

  Double_t x;
  for (Int_t ibin = 0; ibin < p->GetNbinsX(); ibin++) {
    for (Int_t itry = 0; itry < ntry; itry++) {
      for (Int_t ipar = 0; ipar < npars; ipar++)
	f->SetParameter(ipar, gRandom->Gaus(par[ipar], pare[ipar]));
      x = gRandom->Uniform(p->GetXaxis()->GetBinLowEdge(ibin + 1), p->GetXaxis()->GetBinUpEdge(ibin + 1));
      p->Fill(x, f->Eval(x));
    }
  }
}
