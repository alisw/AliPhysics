const Double_t ptmin = 0.2;
const Double_t ptmax = 5.;

void DefineTPChisto(TH2 *h2, const char *yname);
void ALICEWorkInProgress(TCanvas *c, TString date = "today");
void PerformanceSpectrumUncorr(const Char_t *fname = "HFEtask.root"){

  gROOT->SetStyle("Plain");
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.1);
  gStyle->SetTitleY(0.96);

  TFile *in = TFile::Open(fname);
  TList *res = (TList *)in->Get("HFE_Results");
  TList *qa = (TList *)in->Get("HFE_QA");
  gROOT->cd();
  AliHFEcontainer *tcont = dynamic_cast<AliHFEcontainer *>(res->FindObject("trackContainer"));
  AliCFContainer *c = tcont->GetCFContainer("recTrackContReco");
  AliCFContainer *cb = tcont->GetCFContainer("hadronicBackground");

  TH1 *spec = c->Project(c->GetNStep() - 1, 0);
  spec->GetXaxis()->SetTitle("p_{T} / GeV/c");
  spec->GetYaxis()->SetTitle("#frac{dN}{dp_{T}} / (GeV/c)^{-1}");
  spec->GetXaxis()->SetRangeUser(ptmin, ptmax);
  spec->GetYaxis()->SetTitleOffset(1.1);
  spec->SetTitle();
  spec->SetStats(kFALSE);

  spec->SetLineColor(kBlue);
  spec->SetLineWidth(1);
  spec->SetMarkerColor(kBlue);
  spec->SetMarkerStyle(22);

  // Produce background subtracted spectrum
  AliCFDataGrid tracks("tracks", "track grid", *c, c->GetNStep() - 1);
  AliCFDataGrid background("background", "background grid", *cb, 1);
  tracks.ApplyBGCorrection(background);
  TH1 *spec_subtracted = tracks.Project(0);
  spec_subtracted->GetXaxis()->SetTitle("p_{T} / GeV/c");
  spec_subtracted->GetYaxis()->SetTitle("#frac{dN}{dp_{T}} / (GeV/c)^{-1}");
  spec_subtracted->GetXaxis()->SetRangeUser(ptmin, ptmax);
  spec_subtracted->GetYaxis()->SetTitleOffset(1.1);
  spec_subtracted->SetTitle();
  spec_subtracted->SetStats(kFALSE);
  spec_subtracted->SetLineColor(kRed);
  spec_subtracted->SetLineWidth(1);
  spec_subtracted->SetMarkerColor(kRed);
  spec_subtracted->SetMarkerStyle(22);

  TLegend *leg = new TLegend(0.2, 0.25, 0.4, 0.35);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(spec, "Raw Spectrum", "p");
  leg->AddEntry(spec_subtracted, "Spectrum after background subtraction", "p");
 
  TCanvas *c1 = new TCanvas("cspec", "Single-inclusive electron spectrum", 1200, 750);
  c1->cd();
  c1->SetLogy();
  c1->SetGridx(kFALSE);
  c1->SetGridy(kFALSE);
  spec->Draw("ep");
  spec_subtracted->Draw("epsame");
  leg->Draw();
  ALICEWorkInProgress(c1, "today");

  // PID
  TList *pidqa = (TList *)qa->FindObject("HFEpidQA");
  AliHFEtpcPIDqa *tpcqa = (AliHFEtpcPIDqa *)pidqa->FindObject("TPCQA");
  AliHFEtofPIDqa *tofqa = (AliHFEtofPIDqa *)pidqa->FindObject("TOFQA");

  // Make Plots for TPC
  // Create histograms by projecting the THnSparse
  TH2 *hTPCall = tpcqa->MakeSpectrumdEdx(AliHFEdetPIDqa::kBeforePID);
  TH2 *hTPCselected = tpcqa->MakeSpectrumdEdx(AliHFEdetPIDqa::kAfterPID);
  TH2 *hTPCsigmaAll = tpcqa->MakeSpectrumNSigma(AliHFEdetPIDqa::kBeforePID);
  TH2* hTPCsigmaSelected = tpcqa->MakeSpectrumNSigma(AliHFEdetPIDqa::kAfterPID);
  // Make Plots for TOF
  TH2 *hTOFsigmaAll = tofqa->MakeSpectrumNSigma(AliHFEdetPIDqa::kBeforePID);
  TH2 *hTOFsigmaSelected = tofqa->MakeSpectrumNSigma(AliHFEdetPIDqa::kAfterPID);

  hTPCsigmaAll->SetTitle("TPC n#sigma around the electron line");
  hTPCsigmaSelected->SetTitle("TPC n#sigma around the electron line for selected tracks");
  hTOFsigmaAll->SetTitle("TOF n#sigma around the electron line");
  hTOFsigmaSelected->SetTitle("TOF n#sigma around the electron line for selected tracks");
  DefineTPChisto(hTPCall, "TPC Signal / a.u");
  DefineTPChisto(hTPCselected, "TPC Signal / a.u.");
  DefineTPChisto(hTPCsigmaAll, "TPC Sigma");
  DefineTPChisto(hTPCsigmaSelected, "TPC Sigma");

  // Also make nice histograms for TOF
  DefineTPChisto(hTOFsigmaAll, "TOF Sigma");
  DefineTPChisto(hTOFsigmaSelected, "TOF Sigma");

  // Plot them
  TCanvas *c2 = new TCanvas("cTPCall", "TPC Signal for all tracks", 640, 480);
  c2->cd();
  c2->SetGridx(kFALSE);
  c2->SetGridy(kFALSE);
  c2->SetLogx();
  c2->SetLogz();
  hTPCall->GetYaxis()->SetRangeUser(40., 100.);
  hTPCall->Draw("colz");
  ALICEWorkInProgress(c2, "today");

  TCanvas *c3 = new TCanvas("cTPCsel", "TPC Signal for selected tracks", 640, 480);
  c3->cd();
  c3->SetGridx(kFALSE);
  c3->SetGridy(kFALSE);
  c3->SetLogx();
  c3->SetLogz();
  hTPCselected->GetYaxis()->SetRangeUser(40., 100.);
  hTPCselected->Draw("colz");
  ALICEWorkInProgress(c3, "today");

  TCanvas *c4 = new TCanvas("cTPCsigAll", "TPC Sigma for all tracks", 640, 480);
  c4->cd();
  c4->SetGridx(kFALSE);
  c4->SetGridy(kFALSE);
  c4->SetLogx();
  c4->SetLogz();
  //hTPCsigmaAll->GetYaxis()->SetRangeUser(-3.5, 5.);
  hTPCsigmaAll->Draw("colz");
  ALICEWorkInProgress(c4, "today");

  TCanvas *c5 = new TCanvas("cTPCsigSel", "TPC Sigma for selected tracks", 640, 480);
  c5->cd();
  c5->SetGridx(kFALSE);
  c5->SetGridy(kFALSE);
  c5->SetLogx();
  c5->SetLogz();
  hTPCsigmaSelected->GetYaxis()->SetRangeUser(-3.5, 5.);
  hTPCsigmaSelected->Draw("colz");
  ALICEWorkInProgress(c5, "today");

  TCanvas *c6 = new TCanvas("cTOFsigAll", "TOF Sigma for all tracks", 640, 480);
  c6->cd();
  c6->SetGridx(kFALSE);
  c6->SetGridy(kFALSE);
  c6->SetLogx();
  c6->SetLogz();
  hTOFsigmaAll->Draw("colz");
  ALICEWorkInProgress(c6, "today");

  TCanvas *c7 = new TCanvas("cTOFsigSel", "TOF Sigma for selected tracks", 640, 480);
  c7->cd();
  c7->SetGridx(kFALSE);
  c7->SetGridy(kFALSE);
  c7->SetLogx();
  c7->SetLogz();
  //hTOFsigmaSelected->GetYaxis()->SetRangeUser(-3, 3);
  hTOFsigmaSelected->Draw("colz");
  ALICEWorkInProgress(c7, "today");

  TFile *output = new TFile("Performance.root", "RECREATE");
  output->cd();
  spec->Write();
  hTPCall->Write();
  hTPCselected->Write();
  hTPCsigmaAll->Write();
  hTPCsigmaSelected->Write();
  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();
  c5->Write();
  c6->Write();
  c7->Write();
  output->Close();
  delete output;
}

void DefineTPChisto(TH2 *h, const char *yaxis){
  h->SetStats(kFALSE);
  h->GetXaxis()->SetRangeUser(ptmin, ptmax);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent();
  h->GetXaxis()->SetTitle("p / GeV/c");
  h->GetYaxis()->SetTitle(yaxis);
}

void ALICEWorkInProgress(TCanvas *c,TString today){
 //date must be in the form: 04/05/2010
 if(today=="today"){
   TDatime startt;                                                                                                                                                        
   int date=startt.GetDate();

   int y=date/10000;
   int m=(date%10000)/100;
   int d=date%100;


   today="";
   today+=d;
   if(m<10)
     today.Append("/0");
   else today.Append("/");
   today+=m;
   today.Append("/");
   today+=y;  

 }
 TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",0.67,0.65,0.82,0.89);
 //  myPadLogo->SetFillColor(2); 
 myPadLogo->SetBorderMode(0);
 myPadLogo->SetBorderSize(2);
 myPadLogo->SetFrameBorderMode(0);
 myPadLogo->SetLeftMargin(0.0);
 myPadLogo->SetTopMargin(0.0);
 myPadLogo->SetBottomMargin(0.0);
 myPadLogo->SetRightMargin(0.0);
 myPadLogo->Draw();
 myPadLogo->cd();
 TASImage *myAliceLogo = new TASImage("/u/mfasel/work/electron/Spectrum/alice_logo.png");
 myAliceLogo->Draw();
 c->cd();
 TPaveText* t1=new TPaveText(0.59,0.59,0.89,0.66,"NDC");
 t1->SetFillStyle(0);
 t1->SetBorderSize(0);
 t1->AddText(0.,0.,"ALICE Performance");
 t1->SetTextColor(kRed);
 t1->SetTextFont(42);
 t1->Draw();
 TPaveText* t2=new TPaveText(0.59,0.54,0.89,0.60,"NDC");
 t2->SetFillStyle(0);
 t2->SetBorderSize(0);
 t2->SetTextColor(kRed);
 t2->SetTextFont(52);
 t2->AddText(0.,0.,today.Data());
 t2->Draw();
}

