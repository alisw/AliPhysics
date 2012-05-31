// monitor reconstruction efficiencies for different track types
// use output of AliAnalysisTaskEfficiency
// Authors: Veronica Canoa Roman, Eva Sicking

void macrosalida(Int_t prod=0,Int_t trackType=0 )
{

  //postfix of analysis output files ("X")
  TString prodName[6]={"A","B","C","D","E","F"};
  //names of physics lists using in analysed productions (LHC11d6x)
  TString prodNameDetail[6]={"Geant4_QGSP_BERT_EMV_p02",
			     "Geant4_QGSP_BERT_CHIPS_p02",
			     "Geant4_QGSP_BERT_EMV_b01",
			     "Geant4_QGSP_BERT_CHIPS_b01",
			     "Geant4_QGSP_FTFP_BERT_b01",
			     "Geant3"};
  // name of track types
  // ITS = left over tracks from global tracking, efficiency does not make sense
  TString trackName[4]={"Global", "TPC", "ITS_SA", "ITS"};
  
  //open file
  // TFile *f = new TFile(Form("AnalysisResults%s.root",prodName[prod].Data()));
  TFile *f = new TFile(Form("AnalysisResults.root",prodName[prod].Data()));
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
    

  
  // Pt efficiency histos     
  TH1F* fHistRECpt  =  (TH1F*) list->FindObject("fHistRECpt");
  TH1F* fHistMCpt   =  (TH1F*) list->FindObject("fHistMCpt");
  TH1F* fHistFAKEpt =  (TH1F*) list->FindObject("fHistFAKEpt");
  TH1F* fHistMCNRpt =  (TH1F*) list->FindObject("fHistMCNRpt");
  
    
  c1 = new TCanvas(Form("pT_%s_%s",prodNameDetail[prod].Data(),
			trackName[trackType].Data()),
		   Form("pT_%s_%s",prodNameDetail[prod].Data(),
			trackName[trackType].Data()),
		   100, 100, 1000,800 );
  c1->Divide(2,2);
  c1->cd(1);
  c1->GetPad(1)->SetLogy();
    
  fHistMCpt->SetXTitle("p_{T} [GeV]");
  fHistMCpt->SetMinimum(1);
  fHistMCpt->Draw();
  fHistMCpt->SetTitle("");
  fHistMCpt->SetLineWidth(2);
  fHistRECpt->SetLineColor(2);
  fHistRECpt->SetLineWidth(2);
  fHistRECpt->Draw("same");
  fHistFAKEpt->SetLineColor(4);
  fHistFAKEpt->SetLineWidth(2);
  fHistFAKEpt->Draw("same");
  fHistMCNRpt->SetLineColor(3);
  fHistMCNRpt->SetLineWidth(2);
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
  //    fHistPtEff->Add(fHistPtInEff);
  fHistPtEff->SetXTitle("p_{T} [GeV]");
  fHistPtEff->SetTitle("Efficiency");
  fHistPtEff->Draw();
  fHistPtEff->SetMinimum(0);
  fHistPtEff->SetMaximum(1);

  c1->cd(2);
  TLine *line = new TLine(0.07, 0.75, 0.2, 0.75);
  line->SetLineWidth(3);
  line->Draw();
  tex = new TLatex(0.25, 0.73, "MC");
  tex->SetLineWidth(2);
  tex->SetTextSize(0.08);
  tex->Draw();

  line = new TLine(0.07, 0.55, 0.2, 0.55);
  line->SetLineColor(2);
  line->SetLineWidth(3);
  line->Draw();
  tex = new TLatex(0.25, 0.53, "Reconstructed");
  tex->SetLineWidth(2);
  tex->SetTextSize(0.08);
  tex->Draw();

  line = new TLine(0.07, 0.35, 0.2, 0.35);
  line->SetLineColor(4);
  line->SetLineWidth(3);
  line->Draw();
  tex = new TLatex(0.25, 0.33, "Fake Tracks");
  tex->SetLineWidth(3);
  tex->SetTextSize(0.08);
  tex->Draw();


  line = new TLine(0.07, 0.15, 0.2, 0.15);
  line->SetLineColor(3);
  line->SetLineWidth(3);
  line->Draw();
  tex = new TLatex(0.25, 0.13, "MC not reconstructed");
  tex->SetLineWidth(2);
  tex->SetTextSize(0.08);
  tex->Draw();

  // c1->Print(Form("pT_%s_%s.png",prodNameDetail[prod].Data(),trackName[trackType].Data()));

  // Eta efficiency histos     
  c2 = new TCanvas(Form("eta_%s_%s",prodNameDetail[prod].Data(),
			trackName[trackType].Data()),
		   Form("eta_%s_%s",prodNameDetail[prod].Data(),
			trackName[trackType].Data()),
		   100, 100, 1000,600 );  
  c2->Divide(2,2);
  c2->cd(1);

  TH1F* fHistRECeta  =  (TH1F*) list->FindObject("fHistRECeta");
  TH1F* fHistMCeta   =  (TH1F*) list->FindObject("fHistMCeta");
  TH1F* fHistFAKEeta =  (TH1F*) list->FindObject("fHistFAKEeta");
  TH1F* fHistMCNReta =  (TH1F*) list->FindObject("fHistMCNReta");
  fHistMCeta->SetMinimum(0.);
  fHistMCeta->SetXTitle("#eta");
  fHistMCeta->Draw();
  fHistMCeta->SetTitle("");
  fHistRECeta->SetLineColor(2);
  fHistRECeta->Draw("same");
  fHistFAKEeta->SetLineColor(4);
  fHistFAKEeta->Draw("same");
  fHistMCNReta->SetLineColor(3);
  fHistMCNReta->Draw("same");
  c2->cd(4);
  fHistEtaInEff = (TH1F*) fHistMCNReta->Clone();
  fHistEtaInEff->Divide(fHistMCeta);
  fHistEtaInEff->SetXTitle("#eta");
  fHistEtaInEff->SetTitle("Inefficiency from non-rec. MC tracks");
  fHistEtaInEff->Draw();
  fHistEtaInEff->SetMinimum(0);
  fHistEtaInEff->SetMaximum(1);

  c2->cd(3);
  fHistEtaEff = (TH1F*) fHistRECeta->Clone();
  fHistEtaEff->Add(fHistFAKEeta, -1.);
  fHistEtaEff->Divide(fHistMCeta);
  fHistEtaEff->SetXTitle("#eta");
  fHistEtaEff->SetTitle("Efficiency");
  fHistEtaEff->Draw();
  fHistEtaEff->SetMinimum(0);
  fHistEtaEff->SetMaximum(1);

  c2->cd(2);
  TLine *line = new TLine(0.07, 0.75, 0.2, 0.75);
  line->Draw();
  line->SetLineWidth(2);
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
  
  //c2->Print(Form("eta_%s_%s.png",prodNameDetail[prod].Data(),trackName[trackType].Data()));


  // Phi efficiency histos     
  c3 = new TCanvas(Form("phi_%s_%s",prodNameDetail[prod].Data(),
			trackName[trackType].Data()),
		   Form("phi_%s_%s",prodNameDetail[prod].Data(),
			trackName[trackType].Data()),
		   100, 100, 1000,600 );  
  c3->Divide(2,2);
  c3->cd(1);

  TH1F* fHistRECphi  =  (TH1F*) list->FindObject("fHistRECphi");
  TH1F* fHistMCphi   =  (TH1F*) list->FindObject("fHistMCphi");
  TH1F* fHistFAKEphi =  (TH1F*) list->FindObject("fHistFAKEphi");
  TH1F* fHistMCNRphi =  (TH1F*) list->FindObject("fHistMCNRphi");
    
  fHistMCphi->SetMinimum(0.);
  fHistMCphi->SetXTitle("#phi");
  fHistMCphi->Draw();
  fHistMCphi->SetTitle("");
  fHistRECphi->SetLineColor(2);
  fHistRECphi->Draw("same");
  fHistFAKEphi->SetLineColor(4);
  fHistFAKEphi->Draw("same");
  fHistMCNRphi->SetLineColor(3);
  fHistMCNRphi->Draw("same");

  c3->cd(4);
    
  fHistPhiInEff = (TH1F*) fHistMCNRphi->Clone();
  fHistPhiInEff->Divide(fHistMCphi);
  fHistPhiInEff->SetXTitle("#phi");
  fHistPhiInEff->SetTitle("Inefficiency from non-rec. MC tracks");
  fHistPhiInEff->Draw();
  fHistPhiInEff->SetMinimum(0);
  fHistPhiInEff->SetMaximum(1);

  c3->cd(3);
  fHistPhiEff = (TH1F*) fHistRECphi->Clone();
  fHistPhiEff->Add(fHistFAKEphi, -1.);
  fHistPhiEff->Divide(fHistMCphi);
  fHistPhiEff->SetXTitle("#phi");
  fHistPhiEff->SetTitle("Efficiency");
  fHistPhiEff->Draw();
  fHistPhiEff->SetMinimum(0);
  fHistPhiEff->SetMaximum(1);


  c3->cd(2);
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

  // c3->Print(Form("phi_%s_%s.png",prodNameDetail[prod].Data(),trackName[trackType].Data()));


}
