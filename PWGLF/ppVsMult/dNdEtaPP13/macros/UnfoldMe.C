#include "InputData.C"

Float_t integr_eff = 0.;
Float_t integr_eff_err = 0.;
Float_t dndeta_eff = 0.;
Float_t dndeta_eff_err = 0.;

playV0M()
{
  
  TH1 *heff_ev_pho_bay = new TH1F("heff_ev_pho_bay", "", 12, 0, 12);
  TH1 *heff_ch_pho_bay = new TH1F("heff_ch_pho_bay", "", 12, 0, 12);
  TH1 *heff_ev_phy_bay = new TH1F("heff_ev_phy_bay", "", 12, 0, 12);
  TH1 *heff_ch_phy_bay = new TH1F("heff_ch_phy_bay", "", 12, 0, 12);
  TH1 *heff_ev_pho_chi = new TH1F("heff_ev_pho_chi", "", 12, 0, 12);
  TH1 *heff_ch_pho_chi = new TH1F("heff_ch_pho_chi", "", 12, 0, 12);
  TH1 *heff_ev_phy_chi = new TH1F("heff_ev_phy_chi", "", 12, 0, 12);
  TH1 *heff_ch_phy_chi = new TH1F("heff_ch_phy_chi", "", 12, 0, 12);

  for (Int_t i = 0; i < 12; i++) {
  
    TH1 *hpho_bay = UnfoldMe_TAG("LHC10d_ae.pass2", "LHC10d_ae.LHC10f6", "V0M", i, kTRUE, kFALSE, kFALSE, 0.001, 500, AliUnfolding::kPowerLaw, 100., kTRUE);
    heff_ev_pho_bay->SetBinContent(i + 1, integr_eff);
    heff_ev_pho_bay->SetBinError(i + 1, integr_eff_err);
    heff_ch_pho_bay->SetBinContent(i + 1, dndeta_eff);
    heff_ch_pho_bay->SetBinError(i + 1, dndeta_eff_err);
    
    TH1 *hphy_bay = UnfoldMe_TAG("LHC10d_ae.pass2", "LHC10d_ae.LHC10f6a", "V0M", i, kTRUE, kFALSE, kFALSE, 0.001, 500, AliUnfolding::kPowerLaw, 100., kTRUE);
    heff_ev_phy_bay->SetBinContent(i + 1, integr_eff);
    heff_ev_phy_bay->SetBinError(i + 1, integr_eff_err);
    heff_ch_phy_bay->SetBinContent(i + 1, dndeta_eff);
    heff_ch_phy_bay->SetBinError(i + 1, dndeta_eff_err);
    
    TH1 *hpho_chi = UnfoldMe_TAG("LHC10d_ae.pass2", "LHC10d_ae.LHC10f6", "V0M", i, kTRUE, kFALSE, kFALSE, 0.001, 500, AliUnfolding::kPowerLaw, 100., kFALSE);
    heff_ev_pho_chi->SetBinContent(i + 1, integr_eff);
    heff_ev_pho_chi->SetBinError(i + 1, integr_eff_err);
    heff_ch_pho_chi->SetBinContent(i + 1, dndeta_eff);
    heff_ch_pho_chi->SetBinError(i + 1, dndeta_eff_err);
    
    TH1 *hphy_chi = UnfoldMe_TAG("LHC10d_ae.pass2", "LHC10d_ae.LHC10f6a", "V0M", i, kTRUE, kFALSE, kFALSE, 0.001, 500, AliUnfolding::kPowerLaw, 100., kFALSE);
    heff_ev_phy_chi->SetBinContent(i + 1, integr_eff);
    heff_ev_phy_chi->SetBinError(i + 1, integr_eff_err);
    heff_ch_phy_chi->SetBinContent(i + 1, dndeta_eff);
    heff_ch_phy_chi->SetBinError(i + 1, dndeta_eff_err);

  }

  TCanvas *c1 = new TCanvas("c1");
  c1->DrawFrame(0., 0., 12., 1.);
  heff_ev_pho_bay->Draw("same");
  heff_ev_phy_bay->Draw("same");
  heff_ev_pho_chi->Draw("same");
  heff_ev_phy_chi->Draw("same");
}

testMBcorr(Char_t *datatag, Char_t *mctag, Bool_t bayes = kFALSE)
{

  TCanvas *c = new TCanvas("c");
  c->Divide(4, 3);
  for (Int_t i = 0; i < 12; i++) {
    TH1 *hnomb = UnfoldMe_TAG(datatag, mctag, "V0M", i, kFALSE, kFALSE, kFALSE, 1., 4, AliUnfolding::kPowerLaw, 5000., bayes);
    TH1 *hmb = UnfoldMe_TAG(datatag, mctag, "V0M", i, kTRUE, kFALSE, kFALSE, 1., 4, AliUnfolding::kPowerLaw, 5000., bayes);
    hnomb->Divide(hmb);
    c->cd(i+1)->DrawFrame(1., 0.9, 100., 1.1);
    c->cd(i+1)->SetLogx();
    hnomb->Draw("same");
    c->Update();
  }
  
}

playsmoothiter(Char_t *datatag, Char_t *mctag, Char_t *anatag, Int_t bin, Bool_t ismc = kFALSE)
{

  TCanvas *c = new TCanvas("c");
  c->DrawFrame(0., 0.5, 100., 1.5);
  
  TH1 *href = UnfoldMe_TAG(datatag, mctag, anatag, bin, kTRUE, kFALSE, ismc, 1., 4);
  TH1 *hout = NULL;
  TH1 *heff_ev = new TH1F("heff_ev", "", 12, 0, 12);
  TH1 *heff_ch = new TH1F("heff_ch", "", 12, 0, 12);
  for (Int_t j = 2; j < 10; j+=2) {
  for (Int_t i = 5; i < 16; i+=2) {
    printf("%d %d\n", j, i);
    hout = UnfoldMe_TAG(datatag, mctag, anatag, bin, kTRUE, kFALSE, ismc, 0.1 * i, j);
    heff_ev->Fill(integr_eff);
    heff_ch->Fill(dndeta_eff);

    hout->Divide(href);
    c->cd();
    hout->Draw("same");
    c->Update();
  }}

}

playV0M(Char_t *datatag, Char_t *mctag, Bool_t ismc = kFALSE)
{

  TH1 *heff_ev = new TH1F("heff_ev", "", 12, 0, 12);
  TH1 *heff_ch = new TH1F("heff_ch", "", 12, 0, 12);
  for (Int_t i = 0; i < 12; i++) {
    UnfoldMe_TAG(datatag, mctag, "V0M", i, 0, 0, ismc);
    heff_ev->SetBinContent(i + 1, integr_eff);
    heff_ch->SetBinContent(i + 1, dndeta_eff);
  }
  new TCanvas("cEff");
  heff_ev->SetMarkerStyle(21);
  heff_ev->SetMarkerColor(2);
  heff_ev->SetLineColor(2);
  heff_ch->SetMarkerStyle(20);
  heff_ch->SetMarkerColor(4);
  heff_ch->SetLineColor(4);
  heff_ev->Draw("histoc");
  heff_ch->Draw("samehistoc");
  
}

TH1 *
UnfoldMe_TAG(Char_t *datatag, Char_t *mctag, Char_t *anatag, Int_t bin, Bool_t useMBcorr = kTRUE, Bool_t usecorrfit = kFALSE, Bool_t ismc = kFALSE, Float_t smooth = 0.001, Int_t iter = 50, Int_t regul = AliUnfolding::kPowerLaw, Float_t weight = 100., Bool_t bayesian = kTRUE, Int_t nloop = 1)
{

  if (ismc)
    return UnfoldMe(InputMC(datatag), InputMC(mctag), anatag, bin, useMBcorr, usecorrfit, ismc, smooth, iter, regul, weight, bayesian, nloop);
  else
    return UnfoldMe(InputData(datatag), InputMC(mctag), anatag, bin, useMBcorr, usecorrfit, ismc, smooth, iter, regul, weight, bayesian, nloop);
      
}

TH1 *
UnfoldMe(Char_t *data, Char_t *mc, Char_t *anatag, Int_t bin, Bool_t useMBcorr = kTRUE, Bool_t usecorrfit = kFALSE, Bool_t ismc = kFALSE, Float_t smooth = 0.001, Int_t iter = 50, Int_t regul = AliUnfolding::kPowerLaw, Float_t weight = 100., Bool_t bayesian = kTRUE, Int_t nloop = 1)
{

  if (ismc)
    TFile *fdt = TFile::Open(data);
  else
    TFile *fdt = TFile::Open(data);
  TFile *fmc = TFile::Open(mc);
  
  TList *ldt = (TList *)fdt->Get(Form("clist_%s", anatag));
  TList *lmc = (TList *)fmc->Get(Form("clist_%s", anatag));
  
  TH2 *hmatdt = (TH2 *)ldt->FindObject(Form("b%d_corrMatrix", bin));
  if (useMBcorr)
    TH2 *hmatmc = (TH2 *)lmc->FindObject("effMatrix");
  else
    TH2 *hmatmc = (TH2 *)lmc->FindObject(Form("b%d_corrMatrix", bin));
 
  TH1 *hdata = hmatdt->ProjectionY("hdata");
  hdata->Sumw2();
  hdata->SetBinContent(1, 0.);
  hdata->SetBinError(1, 0.);
  //  hdata->Scale(1. / hdata->Integral());
  hdata->SetMarkerStyle(25);
  TH1 *htrue = hmatdt->ProjectionX("htrue");
  htrue->Sumw2();
  //  htrue->Scale(1. / htrue->Integral());
  htrue->SetMarkerStyle(7);
  htrue->SetMarkerColor(2);
  htrue->SetBinContent(1, 0.);
  htrue->SetBinError(1, 0.);
  TH2 *hcorr = (TH2 *)hmatmc->Clone("hcorr");
  TH1 *hinit = (TH1 *)hdata->Clone("hinit");
  TH1 *hresu = (TH1 *)hdata->Clone("hresu");
  TH1 *hbias = (TH1 *)hdata->Clone("hbias");
  hresu->SetMarkerStyle(20);
  hresu->SetMarkerColor(4);
  hresu->Reset();

  TH1 *hnum = hcorr->ProjectionY("hnum");
  TH1 *hden = hcorr->ProjectionY("hden");
  TH1 *heff = hcorr->ProjectionY("heff");
  hnum->Reset();
  hnum->Sumw2();
  hden->Reset();
  hden->Sumw2();
  heff->Reset();
  for (Int_t i = 0; i < heff->GetNbinsX(); i++) {
    Float_t int1 = hcorr->Integral(i + 1, i + 1, 0, -1);
    if (int1 <= 0.) continue;
    Float_t int2 = hcorr->Integral(i + 1, i + 1, 2, -1);
    hnum->SetBinContent(i + 1, int2);
    hnum->SetBinError(i + 1, TMath::Sqrt(int2));
    hden->SetBinContent(i + 1, int1);
    hden->SetBinError(i + 1, TMath::Sqrt(int1));
  }
  new TCanvas("cEfficiency");
  heff->Divide(hnum, hden, 1., 1., "B");
  heff->Draw();
#if 0
  for (Int_t ii = 0; ii < heff->GetNbinsX(); ii++) {
    heff->SetBinContent(ii + 1, 1.);
    heff->SetBinError(ii + 1, 0.);
  }
#endif
  
  for (Int_t i = 0; i < hcorr->GetNbinsX(); i++) {
    hcorr->SetBinContent(i + 1, 1, 0.);
    hcorr->SetBinError(i + 1, 1, 0.);
  }
  for (Int_t i = 0; i < hcorr->GetNbinsY(); i++) {
    hcorr->SetBinContent(1, i + 1, 0.);
    hcorr->SetBinError(1, i + 1, 0.);
  }
  TH2 *hcorrfit = ReturnCorrFromFit(hcorr);

  for (Int_t iloop = 0; iloop < nloop; iloop++) {
    if (bayesian) {
      AliUnfolding::SetUnfoldingMethod(AliUnfolding::kBayesian);
      AliUnfolding::SetBayesianParameters(smooth, iter);
    } else {
      AliUnfolding::SetUnfoldingMethod(AliUnfolding::kChi2Minimization);
      AliUnfolding::SetChi2Regularization(regul, weight);
    }
    AliUnfolding::SetSkip0BinInChi2(kTRUE);
    AliUnfolding::SetSkipBinsBegin(1);
    AliUnfolding::SetNbins(150, 150);
    AliUnfolding::Unfold(usecorrfit ? hcorrfit : hcorr, heff, hdata, hinit, hresu);
    hinit = (TH1 *)hresu->Clone(Form("hinit_%d", iloop));
  }

  printf("hdata->Integral(2, -1) = %f\n", hdata->Integral(2, -1));
  printf("hresu->Integral(2, -1) = %f\n", hresu->Integral(2, -1));
  
  
  TCanvas *cUnfolded = new TCanvas ("cUnfolded", "", 400, 800);
  cUnfolded->Divide(1, 2);
  cUnfolded->cd(1)->SetLogx();
  cUnfolded->cd(1)->SetLogy();
  hdata->Draw();
  hresu->Draw("same");
  htrue->Draw("same");
  cUnfolded->cd(2)->SetLogx();
  cUnfolded->cd(2)->DrawFrame(1., 0.75, 300., 1.25);
  TH1 *hrat = (TH1 *)hresu->Clone("hrat");
  hrat->Divide(htrue);
  hrat->Draw("same");

  TH1 *htrig = (TH1 *)hresu->Clone("htrig");
  htrig->Multiply(heff);

  Float_t dndeta_resu = 0.;
  Float_t integr_resu = 0.;
  Float_t dndeta_trig = 0.;
  Float_t integr_trig = 0.;
  for (Int_t i = 1; i < hresu->GetNbinsX(); i++) {
    dndeta_resu += hresu->GetBinContent(i + 1) * hresu->GetBinLowEdge(i + 1);
    integr_resu += hresu->GetBinContent(i + 1);
    dndeta_trig += htrig->GetBinContent(i + 1) * htrig->GetBinLowEdge(i + 1);
    integr_trig += htrig->GetBinContent(i + 1);
  }
  //  dndeta_resu /= integr_resu;
  //  dndeta_trig /= integr_trig;

  integr_eff = integr_trig / integr_resu;
  integr_eff_err = TMath::Sqrt(integr_eff * (1. - integr_eff) / integr_resu);
  dndeta_eff = dndeta_trig / dndeta_resu;
  dndeta_eff_err = TMath::Sqrt(dndeta_eff * (1. - dndeta_eff) / dndeta_resu);
  
  printf("INEL > 0 efficiency: %.3f +- %.3f\n", integr_eff, integr_eff_err);
  printf("dN/dEta correction:  %.3f +- %.3f\n", dndeta_eff, dndeta_eff_err);

  return hresu;
}

TH2 *
ReturnCorrFromFit(TH2 *hcorr)
{

  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptFit(kTRUE);
  
  TObjArray *oa = new TObjArray();
  hcorr->FitSlicesX(0, 0, -1, 0, "QNR", oa);

  TCanvas *cFit = new TCanvas("cFit");
  cFit->Divide(2, 2);
  
  TF1 *fMean = new TF1("fMean", "[0] + [1] * TMath::Power(x, [2])", 1., 100.);
  fMean->SetParameter(0, 0.);
  fMean->SetParameter(1, 1.);
  fMean->SetParameter(2, 1.);
  TH1 *hMean = (TH1 *)oa->At(1);
  hMean->Fit(fMean, "0q", "I", 1., 100.);
  cFit->cd(1)->SetLogx();
  cFit->cd(1)->SetLogy();
  hMean->Draw();
  fMean->Draw("same");
  
  TF1 *fSigma = new TF1("fSigma", "[0] + [1] * TMath::Power(x, [2])", 1., 100.);
  fSigma->SetParameter(0, 0.);
  fSigma->SetParameter(1, 1.);
  fSigma->SetParameter(2, 0.5);
  TH1 *hSigma = (TH1 *)oa->At(2);
  hSigma->Fit(fSigma, "0q", "", 1., 100.);
  cFit->cd(3)->SetLogx();
  cFit->cd(3)->SetLogy();
  hSigma->Draw();
  fSigma->Draw("same");

  cFit->cd(2)->SetLogx();
  cFit->cd(2)->SetLogy();
  cFit->cd(2)->SetLogz();
  hcorr->Draw("colz");

  TH2 *hcorrfit = (TH2 *)hcorr->Clone("hcorrfit");
  //  hcorrfit->Reset();
  for (Int_t i = 0; i < hcorr->GetNbinsX(); i++) {
    Float_t cent = hcorr->GetXaxis()->GetBinCenter(i + 1);
    Float_t mean = fMean->Eval(cent);
    Float_t sigma = fSigma->Eval(cent);
    if (cent < 25 || cent > 100) continue;
    for (Int_t j = 0; j < 10000; j++) {
      Float_t val = gRandom->Gaus(mean, sigma);
      if (val <= 0.) continue;
      hcorrfit->Fill(val, cent);
    }
  }
  
  cFit->cd(4)->SetLogx();
  cFit->cd(4)->SetLogy();
  cFit->cd(4)->SetLogz();
  hcorrfit->Draw("colz");

  return hcorrfit;
}
