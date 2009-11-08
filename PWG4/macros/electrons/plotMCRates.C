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
  makeData("data/scaled25Oct09/TOTALhistosscaled-LHC09b2-0.root",
           "data/scaled25Oct09/histosscaledLHC09b4AODc.root",
	   "data/scaled25Oct09/histosWboson.root");

  gStyle->SetOptStat(0);
  drawAnnualYields();
  drawPtCutRates();
  drawHadEleRatios();

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
  allMC->SetXTitle("p_{T} (GeV/c)");
  allMC->SetTitle("MC electrons in Pb+Pb, 5.5 TeV");
  allMC->SetYTitle("Annual yield in EMCAL dN/dp_{T} (GeV/c)^{-1}");
  allMC->GetYaxis()->SetRangeUser(1,2.E6);
  allMC->GetXaxis()->SetRangeUser(10.,50.);
  allMC->Draw();
  bMC->Draw("same");
  cMC->Draw("same");  
  cbMC->Draw("same");  
  convMC->Draw("same");  
  dalMC->Draw("same");  
  wzMC->Draw("same");

  leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetTextSize(leg->GetTextSize()*1.2);
  leg->AddEntry(allMC,"All MC electrons","l");
  leg->AddEntry(bMC,"Bottom e","l");
  leg->AddEntry(cMC,"Charm e","l");
  leg->AddEntry(cbMC,"B-->C e","l");
  leg->AddEntry(dalMC,"Dalitz e","l");
  leg->AddEntry(convMC,"Conversion e","l");
  leg->AddEntry(wzMC,"W Boson e","l");
  leg->Draw();
  crates->Print("MCRates_all.pdf");

  TCanvas* crates2 = new TCanvas();
  crates2->Divide(2,4);
  crates2->cd(1); gPad->SetLogy(); allMC->Draw();
  crates2->cd(2); gPad->SetLogy(); bMC->Draw();
  crates2->cd(3); gPad->SetLogy(); cMC->Draw();
  crates2->cd(4); gPad->SetLogy(); cbMC->Draw();
  crates2->cd(5); gPad->SetLogy(); convMC->Draw();
  crates2->cd(6); gPad->SetLogy(); dalMC->Draw();
  crates2->cd(7); gPad->SetLogy(); wzMC->Draw();
  crates2->cd(8); gPad->SetLogy(); mchad->Draw();

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
  TH1F* alleptcut = GetPtCutHisto(allMC);
  TH1F* beleptcut = GetPtCutHisto(bMC);
  TH1F* celeptcut = GetPtCutHisto(cMC);
  TH1F* cbeleptcut = GetPtCutHisto(cbMC);
  TH1F* dalitzptcut = GetPtCutHisto(dalMC);
  TH1F* convptcut = GetPtCutHisto(convMC);
  TH1F* wzptcut = GetPtCutHisto(wzMC);
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
  allheratio->Rebin(2); allheratio->Scale(1./2.);
  behratio->Rebin(2); behratio->Scale(1./2.);
  allheratio->SetLineWidth(2);
  allheratio->GetYaxis()->SetRangeUser(10,2e3);
  allheratio->GetXaxis()->SetRangeUser(10.,49.);
  behratio->GetXaxis()->SetRangeUser(10.,49.);
  allheratio->SetMarkerStyle(20);
  behratio->SetMarkerStyle(24);
  behratio->SetMarkerColor(1);
  allheratio->Draw();
  behratio->Draw("psame");

  TLegend *heleg = new TLegend(0.15,0.15,0.5,0.35);
  heleg->SetFillColor(0);
  heleg->SetTextSize(heleg->GetTextSize()*1.5);
  heleg->AddEntry(allheratio,"All electrons","l");
  heleg->AddEntry(behratio,"Bottom electrons","p");
  heleg->Draw();
  ceh->Print("MCRates_heratio.pdf");
}

