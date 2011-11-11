void checkThr(const char *filename, const char *what="")
{
  TFile *f = TFile::Open(filename);
  TList *list = (TList*)f->Get("coutput");

  TH1F *h3 = (TH1F*)list->FindObject("fV0Percent");

  TH1F *hCent;
  TH1F *hCentAll;
  TH1F *hSemiCent;
  TH1F *hSemiCentAll;
  Double_t nentCent;
  Double_t nentCentall;
  Double_t nentSemiCent;
  Double_t nentSemiCentall;
  TF1 *ffCent;
  TF1 *ffSemiCent;

  hCent = (TH1F*)list->FindObject(Form("fV0Cent%s",what));
  hCentAll = (TH1F*)list->FindObject(Form("fV0Cent%sAll",what));
  hSemiCent = (TH1F*)list->FindObject(Form("fV0SemiCent%s",what));
  hSemiCentAll = (TH1F*)list->FindObject(Form("fV0SemiCent%sAll",what));

  nentCent = hCent->GetEntries();
  Double_t nCentEvts = (Double_t)hCent->Integral(hCent->FindBin(0.0),hCent->FindBin(91.0));
  nentCentall = hCentAll->GetEntries();
  nentSemiCent = hSemiCent->GetEntries();
  Double_t nSemiCentEvts = (Double_t)hSemiCent->Integral(hSemiCent->FindBin(0.0),hSemiCent->FindBin(91.0));
  nentSemiCentall = hSemiCentAll->GetEntries();

  hCent->Sumw2();
  hCent->Divide(hCent,h3,1,1,"B");
  hSemiCent->Sumw2();
  hSemiCent->Divide(hSemiCent,h3,1,1,"B");

  ffCent = new TF1("ffCent","1.-1./(1.+TMath::Exp(-(x-[0])/[1]))",0,100);
  ffCent->SetLineColor(kBlue);
  Double_t par0 = hCent->GetBinCenter(hCent->FindLastBinAbove(0.95));
  ffCent->SetParameter(0,par0);
  ffCent->SetParameter(1,0.3);
  hCent->Fit(ffCent,"R+");
  hCent->Fit(ffCent,"R+");
  hCent->Fit(ffCent,"R+");

  ffSemiCent = new TF1("ffSemiCent","1.-1./(1.+TMath::Exp(-(x-[0])/[1]))",0,100);
  ffSemiCent->SetLineColor(kRed);
  par0 = hSemiCent->GetBinCenter(hSemiCent->FindLastBinAbove(0.95));
  ffSemiCent->SetParameter(0,par0);
  ffSemiCent->SetParameter(1,0.3);
  hSemiCent->Fit(ffSemiCent,"R+");
  hSemiCent->Fit(ffSemiCent,"R+");
  hSemiCent->Fit(ffSemiCent,"R+");

  new TCanvas;
  hSemiCent->SetLineColor(kBlue);
  hSemiCent->SetMarkerColor(kBlue);
  hCent->Draw("e");
  hSemiCent->SetLineColor(kRed);
  hSemiCent->SetMarkerColor(kRed);
  hSemiCent->Draw("e same");

  printf("\n\n\n");

  printf("Purity: Centr      %.3f %.3f   Cuts ( 0.0-%.1f)%%   rate=%.1f Hz\n",
	 ffCent->Integral(0,10)/ffCent->Integral(0,100),
	 ffCent->Integral(0,10)/ffCent->Integral(0,100)*nentCent/nentCentall,
	 ffCent->GetX(0.99),
	 30.*10./(ffCent->Integral(0,10)/ffCent->Integral(0,100)*nentCent/nentCentall));
  printf("Purity: Semi-Centr %.3f %.3f   Cuts ( 0.0-%.1f)%%   rate=%.1f Hz\n",
	 ffSemiCent->Integral(0,50)/ffSemiCent->Integral(0,100),
	 ffSemiCent->Integral(0,50)/ffSemiCent->Integral(0,100)*nentSemiCent/nentSemiCentall,
	 ffSemiCent->GetX(0.99),
	 30.*50./(ffSemiCent->Integral(0,50)/ffSemiCent->Integral(0,100)*nentSemiCent/nentSemiCentall));

  TH1F *hand = (TH1F*)list->FindObject("fV0Percent");
  TH1F *handall = (TH1F*)list->FindObject("fV0PercentAll");
  if (handall) {
    TH1F *hmb63 = (TH1F*)list->FindObject("fV0Percent63");
    TH1F *hmb63all = (TH1F*)list->FindObject("fV0Percent63All");
    printf("V0AND: purity=%f %%\n",hand->Integral(0,hand->FindBin(90.))/(Float_t)handall->GetEntries()*100.);
    printf(" MB63: purity=%f %% eff=%f %%\n",hmb63->Integral(0,hmb63->FindBin(90.))/hmb63all->GetEntries()*100.,hmb63->Integral(0,hmb63->FindBin(90.))/hand->Integral(0,hand->FindBin(90.))*100.);
    hmb63->Sumw2();
    hmb63->Divide(hmb63,hand,1,1,"B");
    new TCanvas;
    hmb63->SetStats(0);
    hmb63->SetTitle("MultA+MultC>=63 efficiency as a function of centrality percentile");
    hmb63->Draw();
  }
}
