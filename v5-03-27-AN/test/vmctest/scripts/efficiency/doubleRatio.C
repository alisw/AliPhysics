// plot macro for double ratio of reconstruction efficiencies of different particles
// Author: Eva Sicking

void doubleRatio(Int_t trackType=0, Int_t particle=3)
{

  //particle names
  TString partName[4]={"Deu", "Tri", "He3", "He4"};

  // track type names
  TString trackName[4]={"Global", "TPC", "ITS_SA", "ITS"};

  //open file and lists
  //output of AliAnalysisTaskEfficiency
  TFile *f = new TFile("AnalysisResults.root");
  TList *list = (TList *)f->Get(Form("QAHists/QAHists_%s",trackName[trackType].Data()));

  //style
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  const Int_t NRGBs = 5;
  const Int_t NCont = 500;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
    
  // Pt     
  TH1F* fHistRECpt  =  (TH1F*) list->
    FindObject(Form("fHistRECptChargeAnti%s",partName[particle].Data()));
  TH1F* fHistMCpt   =  (TH1F*) list->
    FindObject(Form("fHistMCptChargeAnti%s",partName[particle].Data()));
  TH1F* fHistFAKEpt =  (TH1F*) list->
    FindObject(Form("fHistFAKEptChargeAnti%s",partName[particle].Data()));
  TH1F* fHistMCNRpt =  (TH1F*) list->
    FindObject(Form("fHistMCNRptChargeAnti%s",partName[particle].Data()));
  
    
  c1 = new TCanvas(Form("pT_antinuclei_%s",trackName[trackType].Data()),
		   Form("pT_antinuclei_%s",trackName[trackType].Data()),
		   100, 100, 1000,600 );
  c1->Divide(2,2);
  c1->cd(1);
  c1->GetPad(1)->SetLogy();
    
  fHistMCpt->SetXTitle("p_{T} [GeV]");
  fHistMCpt->SetMinimum(1);
  fHistMCpt->Draw();
  fHistMCpt->SetTitle("");
  fHistRECpt->SetLineColor(2);
  fHistRECpt->Draw("same");
  fHistFAKEpt->SetLineColor(4);
  fHistFAKEpt->Draw("same");
  fHistMCNRpt->SetLineColor(3);
  fHistMCNRpt->Draw("same");

  c1->cd(4);
  fHistPtInEff = (TH1F*) fHistMCNRpt->Clone();
  fHistPtInEff->Divide(fHistMCpt);
  fHistPtInEff->SetXTitle("p_{T} [GeV]");
  fHistPtInEff->SetTitle("Inefficiency from non-rec. MC tracks");
  fHistPtInEff->Draw();
  fHistPtInEff->SetMinimum(0);
  fHistPtInEff->SetMaximum(1);

  c1->cd(3);
  fHistPtEff = (TH1F*) fHistRECpt->Clone();
  fHistPtEff->Add(fHistFAKEpt, -1.);
  fHistPtEff->Divide(fHistMCpt);
  //fHistPtEff->Add(fHistPtInEff);
  fHistPtEff->SetXTitle("p_{T} [GeV]");
  fHistPtEff->SetTitle("Efficiency");
  fHistPtEff->Draw();
  fHistPtEff->SetMinimum(0);
  fHistPtEff->SetMaximum(1);

  c1->cd(2);
  TLine *line = new TLine(0.07, 0.75, 0.2, 0.75);
  line->SetLineWidth(2);
  line->Draw();
  tex = new TLatex(0.25, 0.73, "MC");
  tex->SetLineWidth(2);
  tex->Draw();

  line = new TLine(0.07, 0.55, 0.2, 0.55);
  line->SetLineColor(2);
  line->SetLineWidth(2);
  line->Draw();
  tex = new TLatex(0.25, 0.53, "Reconstructed");
  tex->SetLineWidth(2);
  tex->Draw();

  line = new TLine(0.07, 0.35, 0.2, 0.35);
  line->SetLineColor(4);
  line->SetLineWidth(2);
  line->Draw();
  tex = new TLatex(0.25, 0.33, "Fake Tracks");
  tex->SetLineWidth(2);
  tex->Draw();


  line = new TLine(0.07, 0.15, 0.2, 0.15);
  line->SetLineColor(3);
  line->SetLineWidth(2);
  line->Draw();
  tex = new TLatex(0.25, 0.13, "MC not reconstructed");
  tex->SetLineWidth(2);
  tex->Draw();



  // Pt     
  TH1F* fHistRECpt2  =  (TH1F*) list->
    FindObject(Form("fHistRECptCharge%s",partName[particle].Data()));
  TH1F* fHistMCpt2   =  (TH1F*) list->
    FindObject(Form("fHistMCptCharge%s",partName[particle].Data()));
  TH1F* fHistFAKEpt2 =  (TH1F*) list->
    FindObject(Form("fHistFAKEptCharge%s",partName[particle].Data()));
  TH1F* fHistMCNRpt2 =  (TH1F*) list->
    FindObject(Form("fHistMCNRptCharge%s",partName[particle].Data()));
    
  c2 = new TCanvas(Form("pT_nuclei_%s",trackName[trackType].Data()),
		   Form("pT_nuclei_%s",trackName[trackType].Data()),
		   100, 100, 1000,600 );
  c2->Divide(2,2);
  c2->cd(1);
  c2->GetPad(1)->SetLogy();
    
  fHistMCpt2->SetXTitle("p_{T} [GeV]");
  fHistMCpt2->SetMinimum(1);
  fHistMCpt2->Draw();
  fHistMCpt2->SetTitle("");
  fHistRECpt2->SetLineColor(2);
  fHistRECpt2->Draw("same");
  fHistFAKEpt2->SetLineColor(4);
  fHistFAKEpt2->Draw("same");
  fHistMCNRpt2->SetLineColor(3);
  fHistMCNRpt2->Draw("same");

  c2->cd(4);
  fHistPtInEff2 = (TH1F*) fHistMCNRpt2->Clone();
  fHistPtInEff2->Divide(fHistMCpt2);
  fHistPtInEff2->SetXTitle("p_{T} [GeV]");
  fHistPtInEff2->SetTitle("Inefficiency from non-rec. MC tracks");
  fHistPtInEff2->Draw();
  fHistPtInEff2->SetMinimum(0);
  fHistPtInEff2->SetMaximum(1);

  c2->cd(3);
  fHistPtEff2 = (TH1F*) fHistRECpt2->Clone();
  fHistPtEff2->Add(fHistFAKEpt2, -1.);
  fHistPtEff2->Divide(fHistMCpt2);
  //    fHistPtEff->Add(fHistPtInEff);
  fHistPtEff2->SetXTitle("p_{T} [GeV]");
  fHistPtEff2->SetTitle("Efficiency");
  fHistPtEff2->Draw();
  fHistPtEff2->SetMinimum(0);
  fHistPtEff2->SetMaximum(1);

  c2->cd(2);
  TLine *line = new TLine(0.07, 0.75, 0.2, 0.75);
  line->SetLineWidth(2);
  line->Draw();
  tex = new TLatex(0.25, 0.73, "MC");
  tex->SetLineWidth(2);
  tex->Draw();

  line = new TLine(0.07, 0.55, 0.2, 0.55);
  line->SetLineColor(2);
  line->SetLineWidth(2);
  line->Draw();
  tex = new TLatex(0.25, 0.53, "Reconstructed");
  tex->SetLineWidth(2);
  tex->Draw();

  line = new TLine(0.07, 0.35, 0.2, 0.35);
  line->SetLineColor(4);
  line->SetLineWidth(2);
  line->Draw();
  tex = new TLatex(0.25, 0.33, "Fake Tracks");
  tex->SetLineWidth(2);
  tex->Draw();


  line = new TLine(0.07, 0.15, 0.2, 0.15);
  line->SetLineColor(3);
  line->SetLineWidth(2);
  line->Draw();
  tex = new TLatex(0.25, 0.13, "MC not reconstructed");
  tex->SetLineWidth(2);
  tex->Draw();




  new TCanvas;
  
  fHistPtEffPosNeg = (TH1F*) fHistPtEff->Clone();
  fHistPtEffPosNeg->Divide(fHistPtEff2);
  fHistPtEffPosNeg->SetLineWidth(2);
  fHistPtEffPosNeg->Draw();
  fHistPtEffPosNeg->SetMinimum(0);
  fHistPtEffPosNeg->SetMaximum(1.2);
  fHistPtEff->Draw("same");//anti-particle
  fHistPtEff->SetLineColor(kRed);
  fHistPtEff->SetLineWidth(2);
  fHistPtEff2->Draw("same");
  fHistPtEff2->SetLineColor(kBlack);
  fHistPtEff2->SetLineWidth(2);
  TF1 *func = new TF1("func", "1",-5,20);
  func->SetLineWidth(1);
  func->Draw("same");






  TLegend *legp= new TLegend(0.65,0.78,0.9,0.98);
  legp->SetFillColor(kWhite);
  legp->SetBorderSize(0);
  Int_t number =4;
  TF1 *fun[2];
  for(Int_t i=0;i<2;i++){
    fun[i]= new TF1(Form("fun%d",i),"gaus",-5.0,5.0);
    fun[i]->SetLineColor(2**i);
    fun[i]->SetLineStyle(1);
  }
  legp->AddEntry(fun[0],Form("%s",partName[particle].Data() ),"l");   
  legp->AddEntry(fun[1],Form("Anti-%s",partName[particle].Data() ),"l");   



  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);


  TCanvas *c = new TCanvas("c","PtAll" , 100, 100, 880, 680);
  c->Draw();
  c ->cd();
  TPad *pad_spectrum = new TPad( "pad_spectrum","pad_spectrum", 0.0,0.39,1.0,0.95);
  c ->cd();
  TPad *pad_ratio = new TPad( "pad_ratio","pad_ratio", 0.0,0.03,1.0,0.39);  
  
  
  c ->cd();
  pad_spectrum->SetTopMargin(0.);
  pad_spectrum->SetBottomMargin(0.);
  pad_spectrum->SetLeftMargin(0.12);
  pad_spectrum->SetRightMargin(0.055);
  pad_spectrum->Draw();
  pad_spectrum->cd(); 
  TH1F * histo = new TH1F("histo", "",150, -0.5,149.5 );
  histo->GetXaxis()->SetRangeUser(0,11);
  //  histo->GetXaxis()->SetRangeUser(25,35);
  histo->GetYaxis()->SetRangeUser(0.01,1.3);
  histo->SetXTitle("N_{charged}");
  histo->SetTitle("");
  histo->SetYTitle("efficiency  ");
  histo->Draw();
  histo->GetYaxis()->SetNdivisions(10);
  histo->GetYaxis()->SetLabelSize(0.07);
  histo->GetYaxis()->SetTitleSize(0.08);
  histo->GetYaxis()->SetTitleOffset(0.6);
 
  fHistPtEff->Draw("same");
  fHistPtEff2->Draw("same");

  legp->Draw();
  
  

  c ->cd();
  pad_ratio->SetTopMargin(0.);
  pad_ratio->SetBottomMargin(0.36);
  pad_ratio->SetLeftMargin(0.12);
  pad_ratio->SetRightMargin(0.055);
  pad_ratio->Draw();
  pad_ratio->cd(); 
  TH1F * histo2 = new TH1F("histo2", "",150, -0.5,149.5 );
  histo2->GetXaxis()->SetRangeUser(0,90);
  //  histo->GetXaxis()->SetRangeUser(25,35);
  histo2->SetXTitle("p_{T} (GeV/c)");
  histo2->SetTitle("");
  histo2->SetYTitle("Ratio Neg / Pos");
  histo2->Draw();
  histo2->GetXaxis()->SetRangeUser(0.0,11);
  histo2->SetLabelSize(0.1,"X");
  histo2->SetLabelSize(0.1,"Y");
  histo2->SetTitleSize(0.09,"X");
  histo2->SetTitleSize(0.09,"Y");
  histo2->SetTitleOffset(0.3,"Y");
  histo2->GetYaxis()->SetRangeUser(0.59,1.16);
  histo2->GetYaxis()->SetNdivisions(4);
  histo2->GetYaxis()->SetLabelSize(0.12);
  histo2->GetXaxis()->SetLabelSize(0.12);
  histo2->GetXaxis()->SetTitleSize(0.15);
  histo2->GetYaxis()->SetTitleSize(0.11);

  fHistPtEffPosNeg->Draw("same");
  func->Draw("same");
  func->SetLineStyle(2);

  c->Print(Form("%s.png",partName[particle].Data()));

}
