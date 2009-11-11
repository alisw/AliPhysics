/////////////////////////////////////////////////
//
// Macro for plotting rates of identified Non-photonic electrons
// for the EMCAL PPR
//
// J.L. Klay (Cal Poly)
//
/////////////////////////////////////////////////

void plotPIDCompare(char* which = "EMC") {

  gROOT->LoadMacro("makeCombinedData.C");
  makeData("data/scaled25Oct09/TOTALhistosscaled-LHC09b2-0.root",
	   "data/scaled25Oct09/histosscaledLHC09b4AODc.root",
	   "data/scaled25Oct09/histosWboson.root");

  //define common legend
  TLegend* leg = new TLegend(0.5,0.6,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetTextSize(leg->GetTextSize()*1.2);
  leg->AddEntry(alltte,"All N-P e candidates","l");
  //leg->AddEntry(sumtte,"All N-P electrons","l");
  leg->AddEntry(btte,"Bottom e","l");
  leg->AddEntry(ctte,"Charm e","l");
  leg->AddEntry(cbtte,"B-->C e","l");
  leg->AddEntry(daltte,"Dalitz e","l");
  leg->AddEntry(convtte,"Conversion e","l");
  leg->AddEntry(wztte,"W Boson e","l");
  //leg->AddEntry(htte,"Mis-identified hadrons","l");

  gStyle->SetOptStat(0);
  
  TCanvas* crates = new TCanvas("crates","",0,0,1600,600);
  crates->SetFillColor(0);
  crates->SetBorderMode(0);
  crates->SetBorderSize(2);
  crates->SetFrameBorderMode(0);
  crates->SetFrameBorderMode(0);
  crates->Divide(2,1);
  crates->cd(1);
  gPad->SetLogy();
  alltrk->SetXTitle("p_{T} (GeV/c)");
  alltrk->SetYTitle("Annual yield in EMCAL dN/dp_{T} (GeV/c)^{-1}");
  alltrk->SetTitle("PID comparison: Tracking only vs. EMCAL only");
  alltrk->Rebin(5); alltrk->Scale(0.2);
  alltrk->GetYaxis()->SetRangeUser(10.,alltrk->GetMaximum()*2.);
  alltrk->GetXaxis()->SetRangeUser(10.,49.);
  alltrk->Draw();
  htrk->Rebin(5); htrk->Scale(0.2);
  htrk->Draw("same");
  TH1F* tempallemc = (TH1F*)allemc->Clone();
  tempallemc->SetNameTitle("tempallemc","tempallemc");
  tempallemc->SetLineColor(kBlue);
  tempallemc->Rebin(5); tempallemc->Scale(0.2);
  tempallemc->Draw("same");
  
  TH1F* temphemc = (TH1F*)hemc->Clone();
  temphemc->SetNameTitle("temphemc","temphemc");
  temphemc->SetLineColor(kOrange-3);
  temphemc->Rebin(5); temphemc->Scale(0.2);
  temphemc->Draw("same");

  TLegend* leg2 = new TLegend(0.35,0.6,0.9,0.9);
  leg2->SetFillColor(0);
  leg2->SetTextSize(leg->GetTextSize()*1.2);
  leg2->AddEntry(alltrk,"Electron Candidates (Tracking PID only)","l");
  leg2->AddEntry(htrk,"Hadron Contamination (Tracking PID only)","l");
  leg2->AddEntry(tempallemc,"Electron Candidates (EMCAL PID)","l");
  leg2->AddEntry(temphemc,"Hadron Contamination (EMCAL PID)","l");
  leg2->Draw();

  crates->cd(2);
  gPad->SetLogy();
  TH1F* subtrk = (TH1F*)alltrk->Clone();
  for(Int_t i = 1; i <= alltrk->GetNbinsX(); i++) {
    Double_t diff = alltrk->GetBinContent(i) - htrk->GetBinContent(i);
    Double_t unc = 0.;
    if(diff < 0) diff = 0.;
    if(diff > 0) 
      unc = diff/TMath::Sqrt(diff+2.*htrk->GetBinContent(i));
    printf("<%d> Cand %d, Contam %d, diff %d, unc %d\n",i,alltrk->GetBinContent(i),htrk->GetBinContent(i), diff, unc);
    subtrk->SetBinContent(i,diff);
    subtrk->SetBinError(i,unc);
  }
  subtrk->SetYTitle("Annual yield in EMCAL dN/dp_{T} (GeV/c)^{-1}");
  subtrk->SetLineColor(kRed);
  subtrk->SetMarkerStyle(20); subtrk->SetMarkerColor(kRed);
  subtrk->Draw();

  TH1F* subemc = (TH1F*)tempallemc->Clone();
  for(Int_t i = 1; i <= tempallemc->GetNbinsX(); i++) {
    Double_t diff = tempallemc->GetBinContent(i) - temphemc->GetBinContent(i);
    Double_t unc = 0.;
    if(diff < 0.) diff = subemc->GetBinContent(i-1);
    if(diff > 0)
      unc = diff/TMath::Sqrt(diff+2.*temphemc->GetBinContent(i));
    //    printf("Cand %d, Contam %d, diff %d, unc %d\n",tempallemc->GetBinContent(i),temphemc->GetBinContent(i), diff, unc);
    subemc->SetBinContent(i,diff);
    subemc->SetBinError(i,unc);
  }

  TFile* effic = new TFile("elec_eff.root");
  TH1F* heff = (TH1F*)effic->Get("h111");
  heff->Rebin(5); heff->Scale(0.2);
  TH1F* hcorr = (TH1F*)subemc->Clone();
  hcorr->SetName("hcorr");
  for(Int_t i = 1; i < heff->GetNbinsX(); i++) {
    Double_t pt = heff->GetBinCenter(i);
    Double_t eff = heff->GetBinContent(i);
    Double_t corr = 0.;
    if(eff > 0.) corr = hcorr->GetBinContent(i)/eff;
    hcorr->SetBinContent(i,corr);
  }

  Double_t efinal = 0.258;
  TGraphErrors* eerr = new TGraphErrors();
  eerr->SetName("emcErr");
  int count=0;
  for(Int_t i = 1; i <= hcorr->GetNbinsX(); i++) {
    if (hcorr->GetBinCenter(i) <10.0) continue;
    if (hcorr->GetBinCenter(i) >50.0) break;
    cout <<"bin:"<< i << ", bin-center:"<< hcorr->GetBinCenter(i)<< endl;
    eerr->SetPoint(count,hcorr->GetBinCenter(i),hcorr->GetBinContent(i));
    eerr->SetPointError(count,0.,efinal*hcorr->GetBinContent(i));
    if(hcorr->GetBinCenter(i) <20.0) eerr->SetPointError(count,0.,1.5*efinal*hcorr->GetBinContent(i));
    count++;
  }
  eerr->SetFillColor(kRed-8);
  eerr->Draw("3same");
  subtrk->Draw("same");
  hcorr->SetMarkerStyle(20); hcorr->SetMarkerColor(kBlue);
  hcorr->SetLineColor(kBlue);
  hcorr->Draw("same");
  
  allMC->Draw("same");
  mchad->Draw("same");

  TLegend *legx = new TLegend(0.3,0.7,0.9,0.9);
  legx->SetFillColor(0);
  legx->AddEntry(subtrk,"Signal (candidates - contam) (Tracking PID only)","pl");
  legx->AddEntry(hcorr,"Eff. corrected signal (EMCAL PID)","pl");
  legx->AddEntry(eerr,"Systematic uncertainty","f");
  legx->AddEntry(allMC,"MC Electrons","l");
  legx->AddEntry(mchad,"MC Hadrons","l");
  legx->Draw();

  
  TLatex* latex = new TLatex(0.5,0.6,"Unc = #frac{S}{#sqrt{S+2B}}");
  latex->SetNDC();
  //  latex->Draw();

  crates->Print("NPERates_PIDCompare_all.pdf");

  TCanvas* ccomp = new TCanvas("ccomp","",0,0,600,800);
  ccomp->SetFillColor(0);
  ccomp->SetBorderMode(0);
  ccomp->SetBorderSize(2);
  ccomp->SetFrameBorderMode(0);
  ccomp->SetFrameBorderMode(0);
  TPad*    upperPad = new TPad("upperPad", "upperPad", 
			       .005, .2525, .995, .995);
  upperPad->SetFillColor(0);
  upperPad->SetBorderMode(0);
  upperPad->SetBorderSize(2);
  upperPad->SetFrameBorderMode(0);
  upperPad->SetFrameBorderMode(0);
  TPad*    lowerPad = new TPad("lowerPad", "lowerPad", 
			       .005, .005, .995, .2525);
  lowerPad->SetFillColor(0);
  lowerPad->SetFillColor(0);
  lowerPad->SetBorderMode(0);
  lowerPad->SetBorderSize(2);
  lowerPad->SetFrameBorderMode(0);
  lowerPad->SetFrameBorderMode(0);
  upperPad->Draw();        
  lowerPad->Draw(); 
  upperPad->cd();
  gPad->SetLogy();
  hcorr->SetTitle("Efficiency corrected non-photonic electrons");
  hcorr->SetYTitle("Annual yield in EMCAL dN/dp_{T} (GeV/c)^{-1}");
  hcorr->SetXTitle("p_{T} (GeV/c)");
  hcorr->GetXaxis()->SetRangeUser(10.,49.);
  hcorr->GetYaxis()->SetTitleOffset(1.2);
  hcorr->Draw();
  eerr->Draw("3same");
  hcorr->Draw("same");
  allMC->Draw("same");
  TLegend *myleg = new TLegend(0.3,0.7,0.9,0.9);
  myleg->SetFillColor(0);
  myleg->AddEntry(allMC,"MC Electrons","l");
  myleg->AddEntry(hcorr,"Eff. corrected EMCAL N-P electrons","pl");
  myleg->AddEntry(eerr,"Systematic uncertainty","f");
  myleg->Draw();
  lowerPad->cd();
  TH1F* hmc = (TH1F*)allMC->Clone(); hmc->SetName("hmc");
  hmc->Rebin(5); hmc->Scale(0.2);
  TH1F* hmcratio = (TH1F*)hcorr->Clone();  hmcratio->SetName("hmcratio");
  hmcratio->Divide(hmc);
  hmcratio->SetTitle();
  hmcratio->SetYTitle("MC/Data");
  hmcratio->GetYaxis()->SetRangeUser(0.,2.);
  hmcratio->GetYaxis()->SetNdivisions(505);
  hmcratio->GetYaxis()->SetLabelSize(2.*hmcratio->GetYaxis()->GetLabelSize());
  hmcratio->GetYaxis()->SetTitleSize(3.*hmcratio->GetYaxis()->GetTitleSize());
  hmcratio->GetYaxis()->SetTitleOffset(0.3);
  hmcratio->SetXTitle("");
  hmcratio->GetXaxis()->SetLabelSize(2.*hmcratio->GetXaxis()->GetLabelSize());
  hmcratio->Draw();
  TGraphErrors* rerr = new TGraphErrors();
  rerr->SetName("ratioErr");
  count=0;
  for(Int_t i = 1; i <= hmcratio->GetNbinsX(); i++) {
    if (hmcratio->GetBinCenter(i) <10.0) continue;
    if (hmcratio->GetBinCenter(i) >50.0) break;
    rerr->SetPoint(count,hmcratio->GetBinCenter(i),hmcratio->GetBinContent(i));
    rerr->SetPointError(count,0.,efinal*hmcratio->GetBinContent(i));
    if(hmcratio->GetBinCenter(i) <20.0) rerr->SetPointError(count,0.,1.5*efinal*hmcratio->GetBinContent(i));
    count++;
  }
  rerr->SetFillColor(kRed-8);
  rerr->Draw("3same");
  hmcratio->Draw("same");
  TLine* line = new TLine(10.,1.,50.,1.);
  line->SetLineStyle(2);
  line->Draw();
  ccomp->Print("EMCAL_NPE.pdf");
  crates->cd();
}

