/////////////////////////////////////////////////
//
// Macro for plotting rates of identified Non-photonic electrons
// for the EMCAL PPR
//
// J.L. Klay (Cal Poly)
//
/////////////////////////////////////////////////

TLegend* leg;

void plotNPERates(const char* which = "EMC") {

  gROOT->LoadMacro("makeCombinedData.C");
  makeData("data/scaled25Oct09/TOTALhistosscaled-LHC09b2-0.root",
           "data/scaled25Oct09/histosscaledLHC09b4AODc.root",
           "data/scaled25Oct09/histosWboson.root");

  //define common legend
  leg = new TLegend(0.5,0.6,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetTextSize(leg->GetTextSize()*1.2);
  //leg->AddEntry(alltte,"All N-P e candidates","l");
  leg->AddEntry(sumtte,"All N-P electrons","l");
  leg->AddEntry(btte,"Bottom e","l");
  leg->AddEntry(ctte,"Charm e","l");
  leg->AddEntry(cbtte,"B-->C e","l");
  leg->AddEntry(daltte,"Dalitz e","l");
  leg->AddEntry(convtte,"Conversion e","l");
  leg->AddEntry(wztte,"W Boson e","l");
  //  leg->AddEntry(htte,"Misidentified hadrons","l");

  gStyle->SetOptStat(0);
  //  drawAnnualYields(which);
  // drawPtCutRates(which);
  //drawCompareTruth();

}

void drawAnnualYields(char* which = "EMC") {

  TCanvas* crates = new TCanvas();
  crates->SetFillColor(0);
  crates->SetBorderMode(0);
  crates->SetBorderSize(2);
  crates->SetFrameBorderMode(0);
  crates->SetFrameBorderMode(0);

  crates->cd();
  gPad->SetLogy();

  if(strcmp(which,"EMC")==0) {
    allemc->SetXTitle("p_{T} (GeV/c)");
    allemc->SetYTitle("Annual yield");
    allemc->SetTitle("Annual yield of non-phot. electron candidates (EMCAL pid)");
    allemc->Rebin(5); allemc->Scale(0.2);

    sumemc->SetXTitle("p_{T} (GeV/c)");
    sumemc->SetYTitle("Annual yield");
    sumemc->SetTitle("Annual yield of non-phot. electrons (EMCAL pid)");
    sumemc->Rebin(5); sumemc->Scale(0.2);

    bemc->Rebin(5); bemc->Scale(0.2);
    cemc->Rebin(5); cemc->Scale(0.2);
    cbemc->Rebin(5); cbemc->Scale(0.2);
    convemc->Rebin(5); convemc->Scale(0.2);
    dalemc->Rebin(5); dalemc->Scale(0.2);
    wzemc->Rebin(5); wzemc->Scale(0.2);
    hemc->Rebin(5); hemc->Scale(0.2);

    allemc->GetYaxis()->SetRangeUser(1.,2.e6);
    allemc->GetXaxis()->SetRangeUser(10.,49.);
    //allemc->Draw();
    sumemc->GetYaxis()->SetRangeUser(1.,2.e6);
    sumemc->GetXaxis()->SetRangeUser(10.,49.);
    sumemc->Draw();
    bemc->Draw("same");
    cemc->Draw("same");
    cbemc->Draw("same");
    convemc->Draw("same");
    dalemc->Draw("same");
    wzemc->Draw("same");
    //    hemc->Draw("same");
    leg->Draw();
    crates->Print("NPERates_EMC_all.pdf");
  }
  if(strcmp(which,"TRK")==0) {
    alltrk->SetXTitle("p_{T} (GeV/c)");
    alltrk->SetYTitle("Annual yield");
    alltrk->SetTitle("Annual yield of non-phot. electron candidates (TPC+TRD pid)");
    alltrk->Rebin(5); alltrk->Scale(0.2);
    sumtrk->SetXTitle("p_{T} (GeV/c)");
    sumtrk->SetYTitle("Annual yield");
    sumtrk->SetTitle("Annual yield of non-phot. electrons (TPC+TRD pid)");
    sumtrk->Rebin(5); sumtrk->Scale(0.2);

    btrk->Rebin(5); btrk->Scale(0.2);
    ctrk->Rebin(5); ctrk->Scale(0.2);
    cbtrk->Rebin(5); cbtrk->Scale(0.2);
    convtrk->Rebin(5); convtrk->Scale(0.2);
    daltrk->Rebin(5); daltrk->Scale(0.2);
    wztrk->Rebin(5); wztrk->Scale(0.2);
    htrk->Rebin(5); htrk->Scale(0.2);

    alltrk->GetYaxis()->SetRangeUser(1.,6.e6);
    alltrk->GetXaxis()->SetRangeUser(10.,49.);
    //    alltrk->Draw();
    sumtrk->GetYaxis()->SetRangeUser(1.,6.e6);
    sumtrk->GetXaxis()->SetRangeUser(10.,49.);
    sumtrk->Draw();
    btrk->Draw("same");
    ctrk->Draw("same");
    cbtrk->Draw("same");
    convtrk->Draw("same");
    daltrk->Draw("same");
    wztrk->Draw("same");
    //    htrk->Draw("same");
    leg->Draw();
    crates->Print("NPERates_TRK_all.pdf");
  }
  if(strcmp(which,"TTE")==0) {
    alltte->SetXTitle("p_{T} (GeV/c)");
    alltte->SetYTitle("Annual yield");
    alltte->SetTitle("Annual yield of non-phot. electron candidates (Tracking+EMCAL pid)");
    alltte->Rebin(5); alltte->Scale(0.2);
    sumtte->SetXTitle("p_{T} (GeV/c)");
    sumtte->SetYTitle("Annual yield");
    sumtte->SetTitle("Annual yield of non-phot. electrons (Tracking+EMCAL pid)");
    sumtte->Rebin(5); sumtte->Scale(0.2);

    btte->Rebin(5); btte->Scale(0.2);
    ctte->Rebin(5); ctte->Scale(0.2);
    cbtte->Rebin(5); cbtte->Scale(0.2);
    convtte->Rebin(5); convtte->Scale(0.2);
    daltte->Rebin(5); daltte->Scale(0.2);
    wztte->Rebin(5); wztte->Scale(0.2);
    htte->Rebin(5); htte->Scale(0.2);

    alltte->GetYaxis()->SetRangeUser(1.,2.e6);
    alltte->GetXaxis()->SetRangeUser(10.,49.);
    //    alltte->Draw();
    sumtte->GetYaxis()->SetRangeUser(1.,2.e6);
    sumtte->GetXaxis()->SetRangeUser(10.,49.);
    sumtte->Draw();
    btte->Draw("same");
    ctte->Draw("same");
    cbtte->Draw("same");
    convtte->Draw("same");
    daltte->Draw("same");
    wztte->Draw("same");
    htte->Draw("same");
    leg->Draw();
    crates->Print("NPERates_TTE_all.pdf");
  }

}

void drawPtCutRates(char* which = "EMC") {

  TCanvas* cptcut = new TCanvas();
  cptcut->SetFillColor(0);
  cptcut->SetBorderMode(0);
  cptcut->SetBorderSize(2);
  cptcut->SetFrameBorderMode(0);
  cptcut->SetFrameBorderMode(0);
  cptcut->cd();
  gPad->SetLogy();
  if(strcmp(which,"EMC")==0) {
    //    TH1F* alleptcut = GetPtCutHisto(allemc);
    TH1F* alleptcut = GetPtCutHisto(sumemc);
    TH1F* beleptcut = GetPtCutHisto(bemc);
    TH1F* celeptcut = GetPtCutHisto(cemc);
    TH1F* cbeleptcut = GetPtCutHisto(cbemc);
    TH1F* dalitzptcut = GetPtCutHisto(dalemc);
    TH1F* convptcut = GetPtCutHisto(convemc);
    TH1F* wzptcut = GetPtCutHisto(wzemc);
    TH1F* misidptcut = GetPtCutHisto(hemc);
    alleptcut->SetXTitle("p_{T}^{cut} (GeV/c)");
    alleptcut->SetYTitle("Annual Yield in EMCAL for p_{T}>p_{T}^{cut}");
    alleptcut->SetTitle("Pb+Pb, 5.5 TeV reconstructed N-P electrons (EMCAL pid)");
    alleptcut->GetXaxis()->SetRangeUser(10.,49.);
    alleptcut->GetYaxis()->SetRangeUser(1,4.e5);
    alleptcut->Draw();
    beleptcut->Draw("same");
    celeptcut->Draw("same");
    cbeleptcut->Draw("same");
    dalitzptcut->Draw("same");
    convptcut->Draw("same");
    wzptcut->Draw("same");
    //misidptcut->Draw("same");
    leg->Draw();
    cptcut->Print("NPERates_EMC_ptcut_all.pdf");
  }
  if(strcmp(which,"TRK")==0) {
    //    TH1F* alleptcut = GetPtCutHisto(alltrk);
    TH1F* alleptcut = GetPtCutHisto(sumtrk);
    TH1F* beleptcut = GetPtCutHisto(btrk);
    TH1F* celeptcut = GetPtCutHisto(ctrk);
    TH1F* cbeleptcut = GetPtCutHisto(cbtrk);
    TH1F* dalitzptcut = GetPtCutHisto(daltrk);
    TH1F* convptcut = GetPtCutHisto(convtrk);
    TH1F* wzptcut = GetPtCutHisto(wztrk);
    TH1F* misidptcut = GetPtCutHisto(htrk);
    alleptcut->SetXTitle("p_{T}^{cut} (GeV/c)");
    alleptcut->SetYTitle("Annual Yield in EMCAL for p_{T}>p_{T}^{cut}");
    alleptcut->SetTitle("Pb+Pb, 5.5 TeV reconstructed N-P electrons (TPC+TRD pid)");
    alleptcut->GetXaxis()->SetRangeUser(10.,49.);
    alleptcut->GetYaxis()->SetRangeUser(1,6.e6);
    alleptcut->Draw();
    beleptcut->Draw("same");
    celeptcut->Draw("same");
    cbeleptcut->Draw("same");
    dalitzptcut->Draw("same");
    convptcut->Draw("same");
    wzptcut->Draw("same");
    //misidptcut->Draw("same");
    leg->Draw();
    cptcut->Print("NPERates_TRK_ptcut_all.pdf");
  }
  if(strcmp(which,"TTE")==0) {
    //    TH1F* alleptcut = GetPtCutHisto(alltte);
    TH1F* alleptcut = GetPtCutHisto(sumtte);
    TH1F* beleptcut = GetPtCutHisto(btte);
    TH1F* celeptcut = GetPtCutHisto(ctte);
    TH1F* cbeleptcut = GetPtCutHisto(cbtte);
    TH1F* dalitzptcut = GetPtCutHisto(daltte);
    TH1F* convptcut = GetPtCutHisto(convtte);
    TH1F* wzptcut = GetPtCutHisto(wztte);
    TH1F* misidptcut = GetPtCutHisto(htte);
    alleptcut->SetXTitle("p_{T}^{cut} (GeV/c)");
    alleptcut->SetYTitle("Annual Yield in EMCAL for p_{T}>p_{T}^{cut}");
    alleptcut->SetTitle("Pb+Pb, 5.5 TeV reconstructed N-P electrons (Tracking+EMCAL pid)");
    alleptcut->GetXaxis()->SetRangeUser(10.,49.);
    alleptcut->GetYaxis()->SetRangeUser(1,4.e5);
    alleptcut->Draw();
    beleptcut->Draw("same");
    celeptcut->Draw("same");
    cbeleptcut->Draw("same");
    dalitzptcut->Draw("same");
    convptcut->Draw("same");
    wzptcut->Draw("same");
    misidptcut->Draw("same");

    leg->Draw();
    cptcut->Print("NPERates_TTE_ptcut_all.pdf");
  }

}

TH1F* drawCompareTruth() {

  TH1F* mctruth = (TH1F*)belemc->Clone();
  mctruth->SetName("mctruth");
  mctruth->Add(celemc);
  mctruth->Add(candbmc);
  mctruth->Add(wzmc);
  mctruth->Rebin(2); mctruth->Scale(0.5);

  TFile* effic = new TFile("elec_eff.root");
  TH1F* heff = (TH1F*)effic->Get("h111");

  TH1F* hcorr = (TH1F*)sumHFemc->Clone();
  hcorr->SetName("hcorr");
  for(Int_t i = 1; i < heff->GetNbinsX(); i++) {
    Double_t pt = heff->GetBinCenter(i);
    Double_t eff = heff->GetBinContent(i);    
    Double_t corr = 0.;
    if(eff > 0.) corr = hcorr->GetBinContent(i)/eff;
    hcorr->SetBinContent(i,corr);
  }
  hcorr->Rebin(2); hcorr->Scale(0.5);
  sumHFemc->Rebin(2); sumHFemc->Scale(0.5);

  Double_t efinal = 0.258;
  TGraphErrors* eerr = new TGraphErrors();
  eerr->SetName("emcErr");
  for(Int_t i = 1; i <= hcorr->GetNbinsX(); i++) {
    eerr->SetPoint(i-1,hcorr->GetBinCenter(i),hcorr->GetBinContent(i));
    eerr->SetPointError(i-1,0.,efinal*hcorr->GetBinContent(i));
  }
  eerr->SetFillColor(kRed-8);

  TCanvas* ctruth = new TCanvas();
  ctruth->SetFillColor(0);
  ctruth->SetBorderMode(0);
  ctruth->SetBorderSize(2);
  ctruth->SetFrameBorderMode(0);
  ctruth->SetFrameBorderMode(0);

  ctruth->cd();
  gPad->SetLogy();
  mctruth->SetTitle("Comparison of MC and reco HF+W electrons");
  mctruth->SetMarkerColor(kBlack); mctruth->SetLineColor(kBlack);
  mctruth->SetYTitle("Annual yield in EMCAL dN/dp_{T} (GeV/c)^{-1}");
  mctruth->SetXTitle("p_{T} (GeV/c)");
  mctruth->GetXaxis()->SetRangeUser(10.,49.);
  mctruth->Draw();
  eerr->Draw("3same");
  hcorr->SetMarkerColor(kRed);
  hcorr->SetLineColor(kRed);
  hcorr->Draw("same");
  sumHFemc->SetMarkerColor(kBlue);
  sumHFemc->SetLineColor(kBlue);
  sumHFemc->Draw("same");
  TLegend *legy = new TLegend(0.3,0.7,0.9,0.9);
  legy->SetFillColor(0);
  legy->AddEntry(mctruth,"MC HF+W electrons","l");
  legy->AddEntry(sumHFemc,"Rec (EMCAL) HF+W electrons","l");
  legy->AddEntry(hcorr,"Eff. corrected reco HF+W electrons","l");
  legy->AddEntry(eerr,"Systematic uncertainty","f");
  legy->Draw();

  ctruth->Print("NPERates_TruthComparison.pdf");

}
