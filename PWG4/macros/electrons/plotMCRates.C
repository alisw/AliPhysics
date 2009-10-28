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
TH1F* sige;
TH1F* bkge;
TH1F* walle;
TH1F* hije;

TLegend* leg;

void plotMCRates(char* hijfname = "data/scaled25Oct09/histosLHC08d6.root",
		 char* jjfname = "data/scaled25Oct09/TOTALhistosscaled-LHC09b2-0.root",
		 char* bfname = "data/scaled25Oct09/histosscaledLHC09b4AODc.root",
		 char* wfname = "data/scaled25Oct09/histosWBoson.root") {

  //For HIJING need to divide by the number of events, which we
  //can get from the file and do when we perform scaling
  double hijscale = 0.05*(1.E6)*0.5*7700; //0-5% * seconds*lumi*PbPb x-section
  //For bjet and jet-jet events
  double pyscale = (1.E6)*0.5*208*208*100/360; //seconds*lumi*Pb*Pb*acceptance
  double bscale = pyscale; //Do we need the Branching ratio for forced
				//semi-leptonic decays?
  double wscale = pyscale; 
  
  TFile* hijfile = new TFile(hijfname);
  if(!hijfile) { printf("NO HIJING FILE\n"); return; }
  TList* hijlist = (TList*)hijfile->Get("histos");
  TH2F* hijmcele = (TH2F*)histos->FindObject("AnaElectron_hPtMCElectron");
  TH1F* hijmchad = (TH1F*)histos->FindObject("AnaElectron_hPtMCHadron");
  TH1F* hijmult = (TH1F*)histos->FindObject("AnaElectron_hRefMult");
  Int_t nEvt = hijmult->GetEntries();
  if(nEvt == 0) { printf("NO HIJING EVENTS\n"); return; }
  hijmcele->Scale(hijscale/nEvt);
  hijmchad->Scale(hijscale/nEvt);

  TFile* jjfile = new TFile(jjfname);
  if(!jjfile) { printf("NO JET-JET FILE\n"); return; }
  TH2F* jjmcele = (TH2F*)jjfile->Get("AnaElectron_hPtMCElectronScaled");
  TH1F* jjmchad = (TH1F*)jjfile->Get("AnaElectron_hPtMCHadronScaled");
  TH1F* jjmult = (TH1F*)jjfile->Get("AnaElectron_hRefMultScaled");
  Int_t nEvtJJ = jjmult->GetEntries();
  jjmcele->Scale(pyscale);
  jjmchad->Scale(pyscale);

  TFile* bfile = new TFile(bfname);
  if(!bfile) { printf("NO B-JET FILE\n"); return; }
  TH2F* bmcele = (TH2F*)bfile->Get("AnaElectron_hPtMCElectronScaled");
  TH1F* bmchad = (TH1F*)bfile->Get("AnaElectron_hPtMCHadronScaled");
  TH1F* bmult = (TH1F*)bfile->Get("AnaElectron_hRefMultScaled");
  Int_t nEvtB = bmult->GetEntries();
  bmcele->Scale(bscale);
  bmchad->Scale(bscale);

  TFile* wfile = new TFile(wfname);
  if(!wfile) { printf("NO W-BOSON FILE\n"); return; }
  TH2F* wmcele = (TH2F*)wfile->Get("AnaElectron_hPtMCElectron");
  TH1F* wmchad = (TH1F*)wfile->Get("AnaElectron_hPtMCHadron");
  TH1F* wmult = (TH1F*)wfile->Get("AnaElectron_hRefMult");
  Int_t nEvtW = wmult->GetEntries();
  wmcele->Scale(wscale);
  wmchad->Scale(wscale);

  printf("Event statistics: %d (HIJING)  %d (JET-JET)  %d (B-JET)  %d (W-Boson)\n",nEvt,nEvtJJ,nEvtB,nEvtW);

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
  
  all->Add(wz); //because it had to be done separately

  //For comparing contributions
  walle = (TH1F*)wmcele->ProjectionX("walle",7,7);
  sige = (TH1F*)bmcele->ProjectionX("sige",1,1);
  bkge = (TH1F*)jjmcele->ProjectionX("bkge",1,1);
  hije = (TH1F*)hijmcele->ProjectionX("hije",1,1);

  double myscale = 1.; //we already scaled them
  ScaleAndConfigure(all,myscale,kBlack,kFALSE);
  ScaleAndConfigure(bele,myscale,kRed,kFALSE);
  ScaleAndConfigure(sige,myscale,kRed,kFALSE);
  ScaleAndConfigure(cele,myscale,kBlue,kFALSE);
  ScaleAndConfigure(candb,myscale,kViolet,kFALSE);
  ScaleAndConfigure(conv,myscale,kOrange-3,kFALSE);
  ScaleAndConfigure(bkge,myscale,kOrange-3,kFALSE);
  ScaleAndConfigure(dalitz,myscale,kGreen-3,kFALSE);
  ScaleAndConfigure(wz,myscale,kOrange-7,kFALSE);
  ScaleAndConfigure(walle,myscale,kOrange-7,kFALSE);
  ScaleAndConfigure(mchad,myscale,kGreen+2,kFALSE);
  ScaleAndConfigure(hije,myscale,kGreen+2,kFALSE);

  gStyle->SetOptStat(0);
  //drawXSRates();
  drawAnnualYields();
  drawPtCutRates();
  drawHadEleRatios();
  drawSigBkg();

  TFile* mcout = new TFile("MCElectrons.root","RECREATE");
  all->Write();
  bele->Write();
  cele->Write();
  candb->Write();
  conv->Write();
  dalitz->Write();
  wz->Write();
  other->Write();
  mchad->Write();
  sige->Write();
  bkge->Write();
  walle->Write();
  hije->Write();
  mcout->Close();

}

void ScaleAndConfigure(TH1F* hist,Double_t scale, Int_t color,Bool_t keepErr=kFALSE)
{
  hist->Scale(scale);
  hist->SetLineColor(color);
  hist->SetLineWidth(2);
  if(keepErr == kFALSE) {
    //remove the error bars - useful for MC rates
    for(Int_t i = 1; i <= hist->GetNbinsX(); i++) {
      if(hist->GetBinContent(i) > 0.) {
        if(hist->GetBinError(i)/hist->GetBinContent(i) > 0.5) {
          Double_t avg = 0.;
          if(i > 1 && i < hist->GetNbinsX())
            avg = (hist->GetBinContent(i-1) + hist->GetBinContent(i+1))/2.;
          hist->SetBinContent(i,avg);
        }
      }
      hist->SetBinError(i,0.);
    }
  }
}

void drawAnnualYields() {

  TCanvas* crates = new TCanvas();
  crates->cd();
  gPad->SetLogy();
  all->SetXTitle("p_{T} (GeV/c)");
  all->SetTitle("MC electrons in Pb+Pb, 5.5 TeV");
  all->SetYTitle("Annual yield in EMCAL dN/dp_{T} (GeV/c)^{-1}");
  all->GetYaxis()->SetRangeUser(1,2.E6);
  all->GetXaxis()->SetRangeUser(10.,50.);
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
  crates->Print("MCRates_all.pdf");

}

void drawSigBkg() {

  TCanvas* csigbkg = new TCanvas();
  csigbkg->cd();
  gPad->SetLogy();
  all->SetXTitle("p_{T} (GeV/c)");
  all->SetTitle("MC electrons in Pb+Pb, 5.5 TeV");
  all->SetYTitle("Annual Yield in EMCAL dN/dp_{T} (GeV/c)^{-1}");
  all->GetYaxis()->SetRangeUser(1.,6.E8);
  all->GetXaxis()->SetRangeUser(0.,50.);
  all->Draw();
  sige->Draw("same");  
  bkge->Draw("same");  
  hije->Draw("same");  
  walle->Draw("same");

  TLegend* leg1 = new TLegend(0.6,0.6,0.9,0.9);
  leg1->SetTextSize(leg->GetTextSize()*1.2);
  leg1->AddEntry(all,"All MC electrons","l");
  leg1->AddEntry(sige,"B-Jet Events","l");
  leg1->AddEntry(hije,"Pb+Pb Underlying Event","l");
  leg1->AddEntry(bkge,"Jet-Jet Events","l");
  leg1->AddEntry(walle,"W-decay Events","l");
  leg1->Draw();
  csigbkg->Print("MCRates_byEventSource.pdf");

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
  TH1F* wzptcut = GetPtCutHisto(wz);
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
  ceh->cd();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
  TH1F* allratio = (TH1F*)all->Clone();
  TH1F* behratio = (TH1F*)bele->Clone();
  allratio->SetTitle("MC hadrons and electrons in Pb+Pb, 5.5 TeV");
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
  allratio->Rebin(5); allratio->Scale(1./5.);
  behratio->Rebin(5); behratio->Scale(1./5.);
  allratio->GetYaxis()->SetRangeUser(50,1e4);
  allratio->GetXaxis()->SetRangeUser(10.,49.);
  behratio->GetXaxis()->SetRangeUser(10.,49.);
  allratio->SetMarkerStyle(20);
  behratio->SetMarkerStyle(24);
  allratio->Fit("pol0");
  allratio->Draw();
  behratio->Draw("psame");

  TLegend *heleg = new TLegend(0.4,0.75,0.75,0.9);
  heleg->SetTextSize(heleg->GetTextSize()*1.5);
  heleg->AddEntry(allratio,"All electrons","l");
  heleg->AddEntry(behratio,"Bottom electrons","p");
  heleg->Draw();
  ceh->Print("MCRates_heratio.pdf");
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

