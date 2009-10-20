/////////////////////////////////////////////////
//
// Macro for plotting rates of identified Non-photonic electrons
// for the EMCAL PPR
//
// J.L. Klay (Cal Poly)
//
/////////////////////////////////////////////////

TH1F *alltte,      *alltrk,      *allemc;
TH1F *bpttte,      *bpttrk,      *bptemc;
TH1F *cpttte,      *cpttrk,      *cptemc;
TH1F *candbpttte,  *candbpttrk,  *candbptemc;
TH1F *convpttte,   *convpttrk,   *convptemc;
TH1F *dalitzpttte, *dalitzpttrk, *dalitzptemc;
TH1F *wzpttte,     *wzpttrk,     *wzptemc;
TH1F *otherpttte,  *otherpttrk,  *otherptemc;
TH1F *misidpttte,  *misidpttrk,  *misidptemc;

TLegend* leg;

void plotNPERates(char* which = "TTE", 
		  char* hijfname = "scale/histosLHC08d6.root",
		  char* jjfname = "scale/histosscaledLHC09b2ESD.root",
		  char* bfname = "scale/histosscaledLHC09b4.root",
		  char* wfname = "scale/histosWBoson.root") {

  //For HIJING need to divide by the number of events, which we
  //can get from the file and do when we perform scaling
  double hijscale = 0.05*(1.E6)*0.5*7700; //0-5% * seconds*lumi*PbPb x-section
  //For bjet and jet-jet events
  double pyscale = (1.E6)*0.5*208*208*100/360; //seconds*lumi*Pb*Pb*acceptance
  double bscale = pyscale*0.10; //Branching ratio for forced
				//semi-leptonic decays
  double wscale = pyscale*6.29e-05; //JLK: This is temporary X-sec
				     //info from 2-29 GeV bin until we
				     //get pyxsec files; also need to
				     //divide by nTrials.  For now,
				     //use nEvt
  
  TFile* hijfile = new TFile(hijfname);
  if(!hijfile) { printf("NO HIJING FILE\n"); return; }
  TList* hijlist = (TList*)hijfile->Get("histos");
  TH2F* hijtte = (TH2F*)histos->FindObject("AnaElectron_hPtNPEleTTE");
  TH2F* hijemc = (TH2F*)histos->FindObject("AnaElectron_hPtNPEleEMCAL");
  TH2F* hijtrk = (TH2F*)histos->FindObject("AnaElectron_hPtNPEleTPCTRD");
  TH1F* hijmult = (TH1F*)histos->FindObject("AnaElectron_hRefMult");
  Int_t nEvt = hijmult->GetEntries();
  if(nEvt == 0) { printf("NO HIJING EVENTS\n"); return; }
  hijtte->Scale(hijscale/nEvt);
  hijemc->Scale(hijscale/nEvt);
  hijtrk->Scale(hijscale/nEvt);

  TFile* jjfile = new TFile(jjfname);
  if(!jjfile) { printf("NO JET-JET FILE\n"); return; }
  TH2F* jjtte = (TH2F*)jjfile->Get("AnaElectron_hPtNPEleTTEScaled");
  TH2F* jjemc = (TH2F*)jjfile->Get("AnaElectron_hPtNPEleEMCALScaled");
  TH2F* jjtrk = (TH2F*)jjfile->Get("AnaElectron_hPtNPEleTPCTRDScaled");
  TH1F* jjmult = (TH1F*)jjfile->Get("AnaElectron_hRefMultScaled");
  Int_t nEvtJJ = jjmult->GetEntries();
  jjtte->Scale(pyscale);
  jjemc->Scale(pyscale);
  jjtrk->Scale(pyscale);

  TFile* bfile = new TFile(bfname);
  if(!bfile) { printf("NO B-JET FILE\n"); return; }
  TH2F* btte = (TH2F*)bfile->Get("AnaElectron_hPtNPEleTTEScaled");
  TH2F* bemc = (TH2F*)bfile->Get("AnaElectron_hPtNPEleEMCALScaled");
  TH2F* btrk = (TH2F*)bfile->Get("AnaElectron_hPtNPEleTPCTRDScaled");
  TH1F* bmult = (TH1F*)bfile->Get("AnaElectron_hRefMultScaled");
  Int_t nEvtB = bmult->GetEntries();
  btte->Scale(bscale);
  bemc->Scale(bscale);
  btrk->Scale(bscale);

  TFile* wfile = new TFile(wfname);
  if(!wfile) { printf("NO W-BOSON FILE\n"); return; }
  TList* wlist = (TList*)wfile->Get("histos");
  TH2F* wtte = (TH2F*)histos->FindObject("AnaElectron_hPtNPEleTTE");
  TH2F* wemc = (TH2F*)histos->FindObject("AnaElectron_hPtNPEleEMCAL");
  TH2F* wtrk = (TH2F*)histos->FindObject("AnaElectron_hPtNPEleTPCTRD");
  TH1F* wmult = (TH1F*)histos->FindObject("AnaElectron_hRefMult");
  Int_t nEvtW = wmult->GetEntries();
  wtte->Scale(wscale/nEvtW);
  wemc->Scale(wscale/nEvtW);
  wtrk->Scale(wscale/nEvtW);

  printf("Event statistics: %d (HIJING)  %d (JET-JET)  %d (B-JET)  %d (W-Boson)\n",nEvt,nEvtJJ,nEvtB,nEvtW);

  TH2F* combTTE = (TH2F*)hijtte->Clone();
  combTTE->Add(jjtte);
  combTTE->Add(btte);
  combTTE->SetTitle("Identified non-phot. electrons (TPC+TRD+EMCAL)");
  combTTE->SetName("CombinedEleTTE");
  combTTE->SetXTitle("p_T (GeV/c)");

  alltte = (TH1F*)combTTE->ProjectionX("alltte",1,1);
  bpttte = (TH1F*)combTTE->ProjectionX("btte",2,2);
  cpttte = (TH1F*)combTTE->ProjectionX("ctte",3,3);
  candbpttte = (TH1F*)combTTE->ProjectionX("candbtte",4,4);
  convpttte = (TH1F*)combTTE->ProjectionX("convtte",5,5);
  dalitzpttte = (TH1F*)combTTE->ProjectionX("dalitztte",6,6);
  otherpttte = (TH1F*)combTTE->ProjectionX("othertte",8,8);
  misidpttte = (TH1F*)combTTE->ProjectionX("misidtte",9,9);

  wzpttte = (TH1F*)wtte->ProjectionX("wz",7,7);
  alltte->Add(wzpttte);

  double myscale = 1.; //we already scaled them
  ScaleAndConfigure(alltte,myscale,kBlack,kFALSE);
  ScaleAndConfigure(bpttte,myscale,kRed,kFALSE);
  ScaleAndConfigure(cpttte,myscale,kBlue,kFALSE);
  ScaleAndConfigure(candbpttte,myscale,kViolet,kFALSE);
  ScaleAndConfigure(convpttte,myscale,kOrange-3,kFALSE);
  ScaleAndConfigure(dalitzpttte,myscale,kGreen-3,kFALSE);
  ScaleAndConfigure(misidpttte,myscale,kGreen+2,kFALSE);
  ScaleAndConfigure(wzpttte,myscale,kOrange-7,kFALSE);

  TH2F* combEMC = (TH2F*)hijemc->Clone();
  combEMC->Add(jjemc);
  combEMC->Add(bemc);
  combEMC->SetTitle("Identified non-phot. electrons (EMCAL)");
  combEMC->SetName("CombinedEleEMC");
  combEMC->SetXTitle("p_T (GeV/c)");

  allemc = (TH1F*)combEMC->ProjectionX("allemc",1,1);
  bptemc = (TH1F*)combEMC->ProjectionX("bemc",2,2);
  cptemc = (TH1F*)combEMC->ProjectionX("cemc",3,3);
  candbptemc = (TH1F*)combEMC->ProjectionX("candbemc",4,4);
  convptemc = (TH1F*)combEMC->ProjectionX("convemc",5,5);
  dalitzptemc = (TH1F*)combEMC->ProjectionX("dalitzemc",6,6);
  otherptemc = (TH1F*)combEMC->ProjectionX("otheremc",8,8);
  misidptemc = (TH1F*)combEMC->ProjectionX("misidemc",9,9);

  wzptemc = (TH1F*)wemc->ProjectionX("wz",7,7);
  allemc->Add(wzptemc);

  double myscale = 1.; //we already scaled them
  ScaleAndConfigure(allemc,myscale,kBlack,kFALSE);
  ScaleAndConfigure(bptemc,myscale,kRed,kFALSE);
  ScaleAndConfigure(cptemc,myscale,kBlue,kFALSE);
  ScaleAndConfigure(candbptemc,myscale,kViolet,kFALSE);
  ScaleAndConfigure(convptemc,myscale,kOrange-3,kFALSE);
  ScaleAndConfigure(dalitzptemc,myscale,kGreen-3,kFALSE);
  ScaleAndConfigure(misidptemc,myscale,kGreen+2,kFALSE);
  ScaleAndConfigure(wzptemc,myscale,kOrange-7,kFALSE);

  TH2F* combTRK = (TH2F*)hijtrk->Clone();
  combTRK->Add(jjtrk);
  combTRK->Add(btrk);
  combTRK->SetTitle("Identified non-phot. electrons (TPC+TRD)");
  combTRK->SetName("CombinedEleTRK");
  combTRK->SetXTitle("p_T (GeV/c)");

  alltrk = (TH1F*)combTRK->ProjectionX("alltrk",1,1);
  bpttrk = (TH1F*)combTRK->ProjectionX("btrk",2,2);
  cpttrk = (TH1F*)combTRK->ProjectionX("ctrk",3,3);
  candbpttrk = (TH1F*)combTRK->ProjectionX("candbtrk",4,4);
  convpttrk = (TH1F*)combTRK->ProjectionX("convtrk",5,5);
  dalitzpttrk = (TH1F*)combTRK->ProjectionX("dalitztrk",6,6);
  otherpttrk = (TH1F*)combTRK->ProjectionX("othertrk",8,8);
  misidpttrk = (TH1F*)combTRK->ProjectionX("misidtrk",9,9);

  wzpttrk = (TH1F*)wtrk->ProjectionX("wztrk",7,7);
  alltrk->Add(wzpttrk);

  double myscale = 1.; //we already scaled them
  ScaleAndConfigure(alltrk,myscale,kBlack,kFALSE);
  ScaleAndConfigure(bpttrk,myscale,kRed,kFALSE);
  ScaleAndConfigure(cpttrk,myscale,kBlue,kFALSE);
  ScaleAndConfigure(candbpttrk,myscale,kViolet,kFALSE);
  ScaleAndConfigure(convpttrk,myscale,kOrange-3,kFALSE);
  ScaleAndConfigure(dalitzpttrk,myscale,kGreen-3,kFALSE);
  ScaleAndConfigure(misidpttrk,myscale,kGreen+2,kFALSE);
  ScaleAndConfigure(wzpttrk,myscale,kOrange-7,kFALSE);

  //define common legend
  leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->SetTextSize(leg->GetTextSize()*1.2);
  leg->AddEntry(alltte,"All Non-phot. electrons","l");
  leg->AddEntry(bpttte,"Bottom e","l");
  leg->AddEntry(cpttte,"Charm e","l");
  leg->AddEntry(candbpttte,"B-->C e","l");
  leg->AddEntry(dalitzpttte,"Dalitz e","l");
  leg->AddEntry(convpttte,"Conversion e","l");
  leg->AddEntry(wzpttte,"W Boson e","l");
  leg->AddEntry(misidpttte,"Mis-identified hadrons","l");

  gStyle->SetOptStat(0);
  drawAnnualYields(which);
  drawComparePID();
  drawPtCutRates(which);

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
  } else {
    //Set the error bars to statistics of the bin content
    for(Int_t i = 1; i < hist->GetNbinsX(); i++) {
      hist->SetBinError(i,TMath::Sqrt(hist->GetBinContent(i)));
    }
  }

}

void drawAnnualYields(char* which = "TTE") {

  TCanvas* crates = new TCanvas();
  crates->cd();
  gPad->SetLogy();

  if(strcmp(which,"TTE")==0) {
    alltte->SetXTitle("p_{T} (GeV/c)");
    alltte->SetYTitle("Annual yield");
    alltte->SetTitle("Annual yield of non-phot. electrons identified with TPC+TRD+EMCAL");
    alltte->GetYaxis()->SetRangeUser(1.,alltte->GetMaximum()*2.);
    alltte->GetXaxis()->SetRangeUser(0.,50.);
    alltte->Draw();
    bpttte->Draw("same");
    cpttte->Draw("same");
    candbpttte->Draw("same");
    convpttte->Draw("same");
    dalitzpttte->Draw("same");
    wzpttte->Draw("same");
    misidpttte->Draw("same");
    leg->Draw();
  }

}

void drawComparePID() {

  TCanvas* crates = new TCanvas();
  crates->cd();
  gPad->SetLogy();
  alltrk->SetXTitle("p_{T} (GeV/c)");
  alltrk->SetYTitle("Annual yield");
  alltrk->SetTitle("PID comparison: Tracking only vs. with EMCAL");
  alltrk->GetYaxis()->SetRangeUser(1.,alltrk->GetMaximum()*2.);
  alltrk->GetXaxis()->SetRangeUser(0.,50.);
  alltrk->Draw();
  misidpttrk->Draw("same");
  TH1F* tempalltte = (TH1F*)alltte->Clone();
  tempalltte->SetNameTitle("tempalltte","tempalltte");
  tempalltte->SetLineColor(kBlue);
  tempalltte->Draw("same");
  
  TH1F* tempmisidtte = (TH1F*)misidpttte->Clone();
  tempmisidtte->SetNameTitle("tempmisidtte","tempmisidtte");
  tempmisidtte->SetLineColor(kOrange-3);
  tempmisidtte->Draw("same");

  TLegend* leg2 = new TLegend(0.6,0.6,0.9,0.9);
  leg2->SetTextSize(leg->GetTextSize()*1.2);
  leg2->AddEntry(alltrk,"All electrons (Tracking PID only)","l");
  leg2->AddEntry(misidpttrk,"Hadron contamination (Tracking PID only)","l");
  leg2->AddEntry(tempalltte,"All electrons (with EMCAL PID)","l");
  leg2->AddEntry(tempmisidtte,"Hadron contamination (with EMCAL PID)","l");
  leg2->Draw();

}

void drawPtCutRates(char* which = "TTE") {

  TCanvas* cptcut = new TCanvas();
  cptcut->cd();
  gPad->SetLogy();
  if(strcmp(which,"TTE")==0) {
    TH1F* alleptcut = GetPtCutHisto(alltte);
    TH1F* beleptcut = GetPtCutHisto(bpttte);
    TH1F* celeptcut = GetPtCutHisto(cpttte);
    TH1F* cbeleptcut = GetPtCutHisto(candbpttte);
    TH1F* dalitzptcut = GetPtCutHisto(dalitzpttte);
    TH1F* convptcut = GetPtCutHisto(convpttte);
    TH1F* wzptcut = GetPtCutHisto(wzpttte);
    TH1F* misidptcut = GetPtCutHisto(misidpttte);
    alleptcut->GetXaxis()->SetRangeUser(0,50);
    alleptcut->GetYaxis()->SetRangeUser(1,alleptcut->GetMaximum()*2);
    alleptcut->SetXTitle("p_{T}^{cut} (GeV/c)");
    alleptcut->SetYTitle("Annual Yield in EMCAL for p_{T}>p_{T}^{cut}");
    alleptcut->SetTitle("Annual non-phot. electron yield for p_{T}>p_{T}^{cut} (TPC+TRD+EMCAL)");
    alleptcut->Draw();
    beleptcut->Draw("same");
    celeptcut->Draw("same");
    cbeleptcut->Draw("same");
    dalitzptcut->Draw("same");
    convptcut->Draw("same");
    wzptcut->Draw("same");
    misidptcut->Draw("same");

    leg->Draw();
  }

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
    result->SetBinError(i,TMath::Sqrt(val));
  }

  return result;

}

