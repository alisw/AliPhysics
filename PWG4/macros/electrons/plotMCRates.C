/////////////////////////////////////////////////
//
// Macro for plotting MC rates of electrons
// for the EMCAL PPR
//
// J.L. Klay (Cal Poly)
//
/////////////////////////////////////////////////

TH2F* mcele;
TH1F* all;
TH1F* bele;
TH1F* cele;
TH1F* candb;
TH1F* conv;
TH1F* dalitz;
TH1F* wz;
TH1F* other;
TH1F* mchad;
TLegend* leg;

void plotMCRates() {

  bool IsScaled = kTRUE; //PYTHIA
  //bool IsScaled = kFALSE; //HIJING
  double myscale;

  if( IsScaled ){
    myscale = (1.E6)*0.5*208*208;
    printf("Processing histosscaled.root\n");
    TFile* f = new TFile("histosscaledLHC09b4v2.root");
    TList* list = (TList*)f->Get("histosscaled");
    mcele = (TH2F*)histosscaled->FindObject("AnaElectron_hPtMCElectronScaled");
    mchad = (TH1F*)histosscaled->FindObject("AnaElectron_hPtMCHadronScaled");
  }
  else{
    myscale = 0.050*(1.E6)*0.5*208*208/(1.E4);// events per year divided by number events
    printf("Processing histos.root\n");
    TFile* f = new TFile("histos.root");
    TList* list = (TList*)f->Get("histos");
    mcele = (TH2F*)histos->FindObject("AnaElectron_hPtMCElectron");
    mchad = (TH1F*)histosscaled->FindObject("AnaElectron_hPtMCHadronScaled");
  }

  printf("myscale factor = %g\n",myscale);

  all = (TH1F*)mcele->ProjectionX("all",1,1);
  bele = (TH1F*)mcele->ProjectionX("b",2,2);
  cele = (TH1F*)mcele->ProjectionX("c",3,3);
  candb = (TH1F*)mcele->ProjectionX("candb",4,4);
  conv = (TH1F*)mcele->ProjectionX("conv",5,5);
  dalitz = (TH1F*)mcele->ProjectionX("dalitz",6,6);
  wz = (TH1F*)mcele->ProjectionX("wz",7,7);
  other = (TH1F*)mcele->ProjectionX("other",8,8);

  ScaleAndConfigure(all,myscale,kBlack,kFALSE);
  ScaleAndConfigure(bele,myscale,kRed,kFALSE);
  ScaleAndConfigure(cele,myscale,kBlue,kFALSE);
  ScaleAndConfigure(candb,myscale,kViolet,kFALSE);
  ScaleAndConfigure(conv,myscale,kOrange-3,kFALSE);
  ScaleAndConfigure(dalitz,myscale,kGreen-3,kFALSE);
  ScaleAndConfigure(wz,myscale,kOrange-7,kFALSE);
  ScaleAndConfigure(mchad,myscale,kGreen+2,kFALSE);

  gStyle->SetOptStat(0);
  //drawXSRates();
  drawAnnualYields();
  drawPtCutRates();
  drawHadEleRatios();

}

void ScaleAndConfigure(TH1F* hist,Double_t scale, Int_t color,Bool_t keepErr)
{
  hist->Scale(scale);
  hist->SetLineColor(color);
  hist->SetLineWidth(2);
  if(keepErr == kFALSE) {
    //remove the error bars - useful for MC rates
    for(Int_t i = 1; i < hist->GetNbinsX(); i++) {
      hist->SetBinError(i,0.);
    }
  }
}

void drawAnnualYields() {

  TCanvas* crates = new TCanvas("crates","crates",20,20,600,400);
  crates->cd();
  gPad->SetLogy();
  all->SetXTitle("p_{T} (GeV/c)");
  all->SetYTitle("Annual yield");
  all->GetYaxis()->SetRangeUser(1.E1,1.E9);
  all->Draw();
  bele->Draw("same");  
  cele->Draw("same");  
  candb->Draw("same");  
  conv->Draw("same");  
  dalitz->Draw("same");  
  wz->Draw("same");

  leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->SetTextSize(leg->GetTextSize()*1.2);
  leg->AddEntry(all,"All MC electrons","l");
  leg->AddEntry(bele,"Bottom e","l");
  leg->AddEntry(cele,"Charm e","l");
  leg->AddEntry(candb,"B-->C e","l");
  leg->AddEntry(dalitz,"Dalitz e","l");
  leg->AddEntry(conv,"Conversion e","l");
  leg->AddEntry(wz,"W Boson e","l");
  leg->Draw();

}

void drawPtCutRates() {

  TCanvas* cptcut = new TCanvas();
  cptcut->cd();
  gPad->SetLogy();
  TH1F* alleptcut = GetPtCutHisto(all);
  TH1F* beleptcut = GetPtCutHisto(bele);
  TH1F* celeptcut = GetPtCutHisto(cele);
  TH1F* cbeleptcut = GetPtCutHisto(candb);
  TH1F* dalitzptcut = GetPtCutHisto(dalitz);
  TH1F* convptcut = GetPtCutHisto(conv);
  alleptcut->GetXaxis()->SetRangeUser(0,50);
  alleptcut->GetYaxis()->SetRangeUser(1,alleptcut->GetMaximum()*5);
  alleptcut->SetXTitle("p_{T}^{cut} (GeV/c)");
  alleptcut->SetYTitle("Annual Yield in EMCAL for p_{T}>p_{T}^{cut}");
  alleptcut->SetTitle("Annual electron yield in Pb+Pb for p_{T}>p_{T}^{cut}");
  alleptcut->Draw();
  beleptcut->Draw("same");
  celeptcut->Draw("same");
  cbeleptcut->Draw("same");
  dalitzptcut->Draw("same");
  convptcut->Draw("same");
  leg->Draw();

}

void drawHadEleRatios() {

  TCanvas* ceh = new TCanvas();
  ceh->cd();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  TH1F* allratio = (TH1F*)all->Clone();
  TH1F* behratio = (TH1F*)bele->Clone();
  allratio->SetTitle("PYTHIA p+p, 5.5 TeV");
  allratio->SetXTitle("p_{T} (GeV/c)");
  allratio->SetYTitle("N_{hadron}/N_{electron}");
  for(Int_t i = 1; i < all->GetNbinsX(); i++) {
    Double_t vale = all->GetBinContent(i);
    Double_t valb = bele->GetBinContent(i);
    Double_t valh = mchad->GetBinContent(i);
    printf("pT %.2f, Hadron %.1f, Electron %.1f, B-electron %.1f\n",all->GetBinCenter(i),valh,vale,valb);
    if(vale>0) {
      allratio->SetBinContent(i,valh/vale);
      behratio->SetBinContent(i,valh/valb);
    } else {
      allratio->SetBinContent(i,0.);
      behratio->SetBinContent(i,0.);
    }
    allratio->SetBinError(i,0.);
    behratio->SetBinError(i,0.);
  }
  allratio->GetYaxis()->SetRangeUser(1,10000);
  allratio->GetXaxis()->SetRangeUser(0,50);
  allratio->SetMarkerStyle(20);
  behratio->SetMarkerStyle(24);
  //  behratio->Draw();
  allratio->Draw();

  TLegend *heleg = new TLegend(0.5,0.65,0.85,0.85);
  heleg->SetTextSize(heleg->GetTextSize()*1.5);
  heleg->AddEntry(allratio,"All electrons","p");
  heleg->AddEntry(behratio,"Bottom electrons","p");
  heleg->Draw();

}

TH1F* GetPtCutHisto(TH1F* input) 
{
  //Given a rate histogram vs pt, return the histogram with yield
  //above a given pTcut

  TH1F* result = (TH1F*)input->Clone();
  char name[100];
  sprintf(name,"%s_ptCut",result->GetName());
  result->SetNameTitle(name,name);
  for(Int_t i = 1; i <= result->GetNbinsX(); i++) {
    Double_t val = input->Integral(i,result->GetNbinsX());
    result->SetBinContent(i,val);
    result->SetBinError(i,0.);
  }

  return result;

}

