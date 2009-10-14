/////////////////////////////////////////////////
//
// Macro for plotting MC rates of electrons
// for the EMCAL PPR
//
// J.L. Klay (Cal Poly)
//
/////////////////////////////////////////////////

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

void plotMCRates(char* hijfname = "data/histos-merged-LHC0bd6.root",
		 char* jjfname = "data/histosscaledLHC09b2ESDb.root",
		 char* bfname = "data/histosscaledLHC09b4v2.root",
		 char* wfname = "data/histos_wboson_pt29b.root") {

  //For HIJING need to divide by the number of events, which we
  //can get from the file and do when we perform scaling
  double hijscale = 0.05*(1.E6)*0.5*7700; //0-5% * seconds*lumi*PbPb x-section
  //For bjet and jet-jet events
  double pyscale = (1.E6)*0.5*208*208; //seconds*lumi*Pb*Pb
  double wscale = pyscale*6.29e-05; //JLK: This is temporary X-sec
				     //info from 2-29 GeV bin until we get pyxsec files

  TFile* hijfile = new TFile(hijfname);
  if(!hijfile) { printf("NO HIJING FILE\n"); return; }
  TList* hijlist = (TList*)hijfile->Get("histos");
  TH2F* hijmcele = (TH2F*)hijlist->FindObject("AnaElectron_hPtMCElectron");
  TH1F* hijmchad = (TH1F*)hijlist->FindObject("AnaElectron_hPtMCHadron");
  TH1F* refmult = (TH1F*)hijlist->FindObject("AnaElectron_hRefMult");
  Int_t nEvt = refmult->GetEntries();
  if(nEvt == 0) { printf("NO HIJING EVENTS\n"); return; }
  hijmcele->Scale(hijscale/nEvt);
  hijmchad->Scale(hijscale/nEvt);

  TFile* jjfile = new TFile(jjfname);
  if(!jjfile) { printf("NO JET-JET FILE\n"); return; }
  TList* jjlist = (TList*)jjfile->Get("histos");
  TH2F* jjmcele = (TH2F*)histosscaled->FindObject("AnaElectron_hPtMCElectronScaled");
  TH1F* jjmchad = (TH1F*)histosscaled->FindObject("AnaElectron_hPtMCHadronScaled");
  jjmcele->Scale(pyscale);
  jjmchad->Scale(pyscale);

  TFile* bfile = new TFile(bfname);
  if(!bfile) { printf("NO B-JET FILE\n"); return; }
  TList* blist = (TList*)bfile->Get("histos");
  TH2F* bmcele = (TH2F*)histosscaled->FindObject("AnaElectron_hPtMCElectronScaled");
  TH1F* bmchad = (TH1F*)histosscaled->FindObject("AnaElectron_hPtMCHadronScaled");
  bmcele->Scale(pyscale);
  bmchad->Scale(pyscale);

  TFile* wfile = new TFile(wfname);
  if(!wfile) { printf("NO W-BOSON FILE\n"); return; }
  TList* wlist = (TList*)wfile->Get("histos");
  TH2F* wmcele = (TH2F*)histos->FindObject("AnaElectron_hPtMCElectron");
  TH1F* wmchad = (TH1F*)histos->FindObject("AnaElectron_hPtMCHadron");
  wmcele->Scale(wscale);
  wmchad->Scale(wscale);

  TH2F* combined = (TH2F*)hijmcele->Clone();
  combined->Add(jjmcele);
  combined->Add(bmcele);
  combined->SetTitle("MC electrons in Pb+Pb in EMCAL acceptance");
  combined->SetName("CombinedMCEle");
  combined->SetXTitle("p_T (GeV/c)");

  mchad = (TH1F*)hijmchad->Clone();
  mchad->Add(jjmchad);
  mchad->Add(bmchad);
  mchad->Add(wmchad);
  mchad->SetTitle("MC hadrons in Pb+Pb in EMCAL acceptance");
  mchad->SetName("CombinedMCHad");
  mchad->SetXTitle("p_T (GeV/c)");

  all = (TH1F*)combined->ProjectionX("all",1,1);
  bele = (TH1F*)combined->ProjectionX("b",2,2);
  cele = (TH1F*)combined->ProjectionX("c",3,3);
  candb = (TH1F*)combined->ProjectionX("candb",4,4);
  conv = (TH1F*)combined->ProjectionX("conv",5,5);
  dalitz = (TH1F*)combined->ProjectionX("dalitz",6,6);
  other = (TH1F*)combined->ProjectionX("other",8,8);

  wz = (TH1F*)wmcele->ProjectionX("wz",7,7);

  double myscale = 1.; //we already scaled them
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
  all->GetYaxis()->SetRangeUser(1.E1,3.E9);
  all->GetXaxis()->SetRangeUser(0.,50.);
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
  alleptcut->GetYaxis()->SetRangeUser(100,alleptcut->GetMaximum()*2);
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
  allratio->SetTitle("Pb+Pb, 5.5 TeV");
  allratio->SetXTitle("p_{T} (GeV/c)");
  allratio->SetYTitle("Hadrons/Electrons");
  for(Int_t i = 1; i < all->GetNbinsX(); i++) {
    Double_t vale = all->GetBinContent(i);
    Double_t valb = bele->GetBinContent(i);
    Double_t valh = mchad->GetBinContent(i);
    //printf("pT %.2f, Hadron %.1f, Electron %.1f, B-electron %.1f\n",all->GetBinCenter(i),valh,vale,valb);
    if(vale>0) allratio->SetBinContent(i,valh/vale);
    else allratio->SetBinContent(i,0.);

    if(valb>0) behratio->SetBinContent(i,valh/valb);
    else behratio->SetBinContent(i,0.);

    allratio->SetBinError(i,0.);
    behratio->SetBinError(i,0.);
  }
  allratio->Rebin();
  behratio->Rebin();
  allratio->GetYaxis()->SetRangeUser(1,10000);
  allratio->GetXaxis()->SetRangeUser(0,50);
  behratio->GetXaxis()->SetRangeUser(0,50);
  allratio->SetMarkerStyle(20);
  behratio->SetMarkerStyle(24);
  allratio->Draw();
  behratio->Draw("psame");


  TLegend *heleg = new TLegend(0.4,0.75,0.75,0.9);
  heleg->SetTextSize(heleg->GetTextSize()*1.5);
  heleg->AddEntry(allratio,"All electrons","l");
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

