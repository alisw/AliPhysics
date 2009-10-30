/////////////////////////////////////////////////
//
// Macro for plotting MC rates of electrons
// for the EMCAL PPR
//
// J.L. Klay (Cal Poly)
//
/////////////////////////////////////////////////

TLegend* leg;

void plotMCRates() {

  gROOT->LoadMacro("makeCombinedData.C");
  makeData("data/scaled25Oct09/histosLHC08d6.root",
           "data/scaled25Oct09/TOTALhistosscaled-LHC09b2-0.root",
           "data/scaled25Oct09/histosscaledLHC09b4AODc.root",
           "data/scaled25Oct09/histosWboson.root");

  gStyle->SetOptStat(0);
  //drawXSRates();
  drawAnnualYields();
  drawPtCutRates();
  drawHadEleRatios();
  drawSigBkg();

}

void drawAnnualYields() {

  TCanvas* crates = new TCanvas();
  crates->SetFillColor(0);
  crates->SetBorderMode(0);
  crates->SetBorderSize(2);
  crates->SetFrameBorderMode(0);
  crates->SetFrameBorderMode(0);

  crates->cd();
  gPad->SetLogy();
  allmc->SetXTitle("p_{T} (GeV/c)");
  allmc->SetTitle("MC electrons in Pb+Pb, 5.5 TeV");
  allmc->SetYTitle("Annual yield in EMCAL dN/dp_{T} (GeV/c)^{-1}");
  allmc->GetYaxis()->SetRangeUser(1,2.E6);
  allmc->GetXaxis()->SetRangeUser(10.,50.);
  allmc->Draw();
  belemc->Draw("same");  
  celemc->Draw("same");  
  candbmc->Draw("same");  
  convmc->Draw("same");  
  dalmc->Draw("same");  
  wzmc->Draw("same");

  leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetTextSize(leg->GetTextSize()*1.2);
  leg->AddEntry(allmc,"All MC electrons","l");
  leg->AddEntry(belemc,"Bottom e","l");
  leg->AddEntry(celemc,"Charm e","l");
  leg->AddEntry(candbmc,"B-->C e","l");
  leg->AddEntry(dalmc,"Dalitz e","l");
  leg->AddEntry(convmc,"Conversion e","l");
  leg->AddEntry(wzmc,"W Boson e","l");
  leg->Draw();
  crates->Print("MCRates_all.pdf");

}

void drawSigBkg() {

  TCanvas* csigbkg = new TCanvas();
  csigbkg->SetFillColor(0);
  csigbkg->SetBorderMode(0);
  csigbkg->SetBorderSize(2);
  csigbkg->SetFrameBorderMode(0);
  csigbkg->SetFrameBorderMode(0);

  csigbkg->cd();
  gPad->SetLogy();
  allmc->SetXTitle("p_{T} (GeV/c)");
  allmc->SetTitle("MC electrons in Pb+Pb, 5.5 TeV");
  allmc->SetYTitle("Annual Yield in EMCAL dN/dp_{T} (GeV/c)^{-1}");
  allmc->GetYaxis()->SetRangeUser(1.,6.E8);
  allmc->GetXaxis()->SetRangeUser(0.,50.);
  allmc->Draw();
  sigemc->Draw("same");  
  bkgemc->Draw("same");  
  hijemc->Draw("same");  
  wallemc->Draw("same");

  TLegend* leg1 = new TLegend(0.6,0.6,0.9,0.9);
  leg1->SetFillColor(0);
  leg1->SetTextSize(leg->GetTextSize()*1.2);
  leg1->AddEntry(allmc,"All MC electrons","l");
  leg1->AddEntry(sigemc,"B-Jet Events","l");
  leg1->AddEntry(hijemc,"Pb+Pb Underlying Event","l");
  leg1->AddEntry(bkgemc,"Jet-Jet Events","l");
  leg1->AddEntry(wallemc,"W-decay Events","l");
  leg1->Draw();
  csigbkg->Print("MCRates_byEventSource.pdf");

}

void drawPtCutRates() {

  TCanvas* cptcut = new TCanvas();
  cptcut->SetFillColor(0);
  cptcut->SetBorderMode(0);
  cptcut->SetBorderSize(2);
  cptcut->SetFrameBorderMode(0);
  cptcut->SetFrameBorderMode(0);

  cptcut->cd();
  gPad->SetLogy();
  TH1F* alleptcut = GetPtCutHisto(allmc);
  TH1F* beleptcut = GetPtCutHisto(belemc);
  TH1F* celeptcut = GetPtCutHisto(celemc);
  TH1F* cbeleptcut = GetPtCutHisto(candbmc);
  TH1F* dalitzptcut = GetPtCutHisto(dalmc);
  TH1F* convptcut = GetPtCutHisto(convmc);
  TH1F* wzptcut = GetPtCutHisto(wzmc);
  alleptcut->GetXaxis()->SetRangeUser(10,50);
  alleptcut->GetYaxis()->SetRangeUser(10,2.e6);
  alleptcut->SetXTitle("p_{T}^{cut} (GeV/c)");
  alleptcut->SetYTitle("Annual Yield in EMCAL for p_{T}>p_{T}^{cut}");
  alleptcut->SetTitle("MC electrons in Pb+Pb, 5.5 TeV");
  alleptcut->Draw();
  beleptcut->Draw("same");
  celeptcut->Draw("same");
  cbeleptcut->Draw("same");
  dalitzptcut->Draw("same");
  convptcut->Draw("same");
  wzptcut->Draw("same");
  leg->Draw();
  cptcut->Print("MCRates_ptcut_all.pdf");

}

void drawHadEleRatios() {

  TCanvas* ceh = new TCanvas();
  ceh->SetFillColor(0);
  ceh->SetBorderMode(0);
  ceh->SetBorderSize(2);
  ceh->SetFrameBorderMode(0);
  ceh->SetFrameBorderMode(0);

  ceh->cd();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  TH1F* allratio = (TH1F*)allmc->Clone();
  TH1F* behratio = (TH1F*)belemc->Clone();
  allratio->SetTitle("MC hadrons and electrons in Pb+Pb, 5.5 TeV");
  allratio->SetXTitle("p_{T} (GeV/c)");
  allratio->SetYTitle("Hadrons/Electrons");
  for(Int_t i = 1; i < allmc->GetNbinsX(); i++) {
    Double_t vale = allmc->GetBinContent(i);
    Double_t valb = belemc->GetBinContent(i);
    Double_t valh = mchad->GetBinContent(i);
    //printf("pT %.2f, Hadron %.1f, Electron %.1f, B-electron %.1f\n",all->GetBinCenter(i),valh,vale,valb);
    if(vale>0) allratio->SetBinContent(i,valh/vale);
    else allratio->SetBinContent(i,0.);

    if(valb>0) behratio->SetBinContent(i,valh/valb);
    else behratio->SetBinContent(i,0.);

    allratio->SetBinError(i,0.);
    behratio->SetBinError(i,0.);
  }
  allratio->Rebin(5); allratio->Scale(1./5.);
  behratio->Rebin(5); behratio->Scale(1./5.);
  allratio->GetYaxis()->SetRangeUser(50,1e4);
  allratio->GetXaxis()->SetRangeUser(10.,49.);
  behratio->GetXaxis()->SetRangeUser(10.,49.);
  allratio->SetMarkerStyle(20);
  behratio->SetMarkerStyle(24);
  behratio->SetMarkerColor(1);
  allratio->Draw();
  behratio->Draw("psame");

  TLegend *heleg = new TLegend(0.4,0.75,0.75,0.9);
  heleg->SetFillColor(0);
  heleg->SetTextSize(heleg->GetTextSize()*1.5);
  heleg->AddEntry(allratio,"All electrons","l");
  heleg->AddEntry(behratio,"Bottom electrons","p");
  heleg->Draw();
  ceh->Print("MCRates_heratio.pdf");
}

