void extractThr(const char *filename, const char *fparname, Float_t centrValue = 10.0, Float_t semiCentrValue = 50.0)
{
  gROOT->LoadMacro(Form("%s/AliAnaVZEROTrigger.cxx+",
  			gSystem->pwd()));

  AliAnaVZEROTrigger *task = new AliAnaVZEROTrigger;
  task->Setup(fparname);

  TFile *f = TFile::Open(filename);

  TList *list = (TList*)f->Get("coutput");

  TH2F *h1 = (TH2F*)list->FindObject("fV0Charge2d");
  TH2F *h2 = (TH2F*)list->FindObject("fV0Charge2dPercent");
  h2->Divide(h1);
  TProfile *h1_pfx = h1->ProfileX("h1_pfx");
  TProfile *h1_pfy = h1->ProfileY("h1_pfy");
  new TCanvas;
  gStyle->SetPalette(1,0);
  h2->Draw("colz");
  new TCanvas;
  TF1 *f1 = new TF1("f1","[0]*x",0,task->GetCentCuts(0));
  f1->SetParameter(0,task->GetRatio());
  h1_pfx->Fit(f1,"WR+");
  h1_pfx->Fit(f1,"WR+");
  h1_pfx->Fit(f1,"WR+");
  h1_pfx->Draw();

  new TCanvas;
  TF1 *f2 = new TF1("f2","x/[0]",0,task->GetCentCuts(0)*f1->GetParameter(0));
  f2->SetParameter(0,1./task->GetRatio());
  h1_pfy->Fit(f2,"WR+");
  h1_pfy->Fit(f2,"WR+");
  h1_pfy->Fit(f2,"WR+");
  h1_pfy->Draw();

  TH1F *h3 = (TH1F*)list->FindObject("fV0Percent");
  new TCanvas;
  h3->Draw();
  Double_t nTotalEvts = (Double_t)h3->Integral(h3->FindBin(0.0),h3->FindBin(91.0))/0.9;

  Int_t nb = task->GetNThr();
  new TCanvas;
  h3->Sumw2();
  TH1F **hh = new TH1F*[nb];
  TH1F **hhall = new TH1F*[nb];
  TF1 **ff = new TF1*[nb];
  Double_t *nent = new Double_t[nb];
  Double_t *nBinEvts = new Double_t[nb];
  Double_t *nentall = new Double_t[nb];
  for(Int_t i = 0; i < nb; ++i) {
    hh[i] = (TH1F*)list->FindObject(Form("fV0PercentBins_%d",i));
    hhall[i] = (TH1F*)list->FindObject(Form("fV0PercentBinsAll_%d",i));

    nent[i] = hh[i]->GetEntries();
    nBinEvts[i] = (Double_t)hh[i]->Integral(hh[i]->FindBin(0.0),hh[i]->FindBin(91.0));
    nentall[i] = hhall[i]->GetEntries();

    hh[i]->Sumw2();
    hh[i]->Divide(hh[i],h3,1,1,"B");
    //    TF1 *ff = new TF1(Form("ff_%d",i),"x > [0] ? TMath::Exp(-(x-[0])*(x-[0])/(2.*[1]*[1])) : 1.0",0,100);
    ff[i] = new TF1(Form("ff_%d",i),"1.-1./(1.+TMath::Exp(-(x-[0])/[1]))",0,100);
    ff[i]->SetLineColor(i%9+1);
    printf("%d\n",i);
    if (nent[i] > 100) {
      Double_t par0 = hh[i]->GetBinCenter(hh[i]->FindLastBinAbove(0.9));
      printf("%f\n",par0);
      ff[i]->SetParameter(0,par0);
      ff[i]->SetParameter(1,0.3);
      hh[i]->Fit(ff[i],"R+");
      hh[i]->Fit(ff[i],"R+");
      hh[i]->Fit(ff[i],"R+");
    }
  }

  new TCanvas;
  Double_t *eff99 = new Double_t[nb];
  Double_t *eff05 = new Double_t[nb];
  Double_t *purity = new Double_t[nb];
  Double_t *purity2 = new Double_t[nb];
  Bool_t first = kTRUE;
  for(Int_t i = 0; i < nb; ++i) {
    hh[i]->SetLineColor(i%9+1);
    hh[i]->SetMarkerColor(i%9+1);
    hh[i]->SetLineWidth(2.0);
    hh[i]->SetStats(0);
    if (first) {
      first = kFALSE;
      hh[i]->Draw("e");
      hh[i]->GetXaxis()->SetTitle("Centrality percentile");
      hh[i]->GetYaxis()->SetTitle("Efficiency");
    }
    else
      hh[i]->Draw("e same");
    eff99[i] = ff[i]->GetX(0.99);
    eff05[i] = ff[i]->GetX(0.05)-eff99[i];
    purity[i] = ff[i]->Integral(0,eff99[i])/ff[i]->Integral(0,100);
    purity2[i] = purity[i]*nent[i]/nentall[i];
    //    printf("%d  %d %d\n",i,nent[i],nentall[i]);
  }

  new TCanvas;
  TGraph *gr99 = new TGraph(nb,eff99,eff05);
  gr99->SetMarkerStyle(21);
  gr99->GetHistogram()->GetXaxis()->SetTitle("Centrality percentile @ Eff=99%");
  gr99->GetHistogram()->GetYaxis()->SetTitle("#Delta(Centrality percentile @ Eff=5% - @ Eff=99%)");
  gr99->GetHistogram()->GetYaxis()->SetRangeUser(-7,10);
  gr99->SetTitle("");
  gr99->Draw("AP");

  TLegend *leg0 = new TLegend(0.75,0.55,0.89,0.89);
  leg0->AddEntry(gr99,"A and C sides","P");
  leg0->Draw();

  new TCanvas;
  TGraph *pur = new TGraph(nb,eff99,purity);
  pur->SetMarkerStyle(21);
  pur->GetHistogram()->GetXaxis()->SetTitle("Centrality percentile @ Eff=99%");
  pur->GetHistogram()->GetYaxis()->SetTitle("Purity");
  pur->SetTitle("");
  TF1 *fitPur = new TF1("fitPur","[0]-[1]*TMath::Exp(-x/[2])",6,50);
  fitPur->SetParameter(0,0.9);
  fitPur->SetParameter(1,0.4);
  fitPur->SetParameter(2,9.0);
  pur->Fit(fitPur,"R+");
  pur->Fit(fitPur,"R+");
  pur->Fit(fitPur,"R+");
  pur->Draw("AP");

  TGraph *pur2 = new TGraph(nb,eff99,purity2);
  pur2->SetMarkerStyle(22);
  pur2->GetHistogram()->GetXaxis()->SetTitle("Centrality percentile @ Eff=99%");
  pur2->GetHistogram()->GetYaxis()->SetTitle("Purity (no phys & centr selection)");
  pur2->SetTitle("");
  TF1 *fitPur2 = new TF1("fitPur2","[0]-[1]*TMath::Exp(-x/[2])",6,50);
  fitPur2->SetParameter(0,0.9);
  fitPur2->SetParameter(1,0.4);
  fitPur2->SetParameter(2,9.0);
  pur2->Fit(fitPur2,"R+");
  pur2->Fit(fitPur2,"R+");
  pur2->Fit(fitPur2,"R+");
  pur2->Draw("Psame");

  TLegend *leg = new TLegend(0.75,0.55,0.89,0.89);
  leg->AddEntry(pur,"A and C, Phys & centr selection","P");
  leg->AddEntry(pur2,"A and C, No event selection","P");
  leg->Draw();

  new TCanvas;
  Double_t *thrA = new Double_t[nb];
  Double_t *thrC = new Double_t[nb];
  Double_t ratio = task->GetRatio();
  for(Int_t i = 0; i < nb; ++i) {
    thrA[i] = task->GetThrA(i);
    thrC[i] = task->GetThrC(i);
  }
  TGraph *thr = new TGraph(nb,eff99,thrA);
  thr->SetMarkerStyle(21);
  TF1 *fThr = new TF1("fThr","[0]*x*x*x+[1]*x*x+[2]*x+[3]",1.5,50);
  thr->Fit(fThr,"R");
  thr->Draw("AP");

  new TCanvas;
  Double_t *intThr = new Double_t[nb];
  for(Int_t i = 0; i < nb; ++i) {
    intThr[i] = 100.*nBinEvts[i]/nTotalEvts;
  }
  TGraph *grInt = new TGraph(nb,intThr,thrA);
  grInt->SetMarkerStyle(21);
  TF1 *fIntThr = new TF1("fIntThr","[0]*x*x*x+[1]*x*x+[2]*x+[3]",1.5,50);
  grInt->Fit(fIntThr,"R");
  grInt->Draw("AP");

  new TCanvas;
  TH2F *hall = (TH2F*)list->FindObject("fV0Charge2dAll");
  h2->Draw("colz");
  TLine *line0 = new TLine(fThr->Eval(centrValue),fThr->Eval(centrValue)*task->GetRatio(),25000,fThr->Eval(centrValue)*task->GetRatio());
  line0->Draw("same");
  TLine *line1 = new TLine(fThr->Eval(centrValue),fThr->Eval(centrValue)*task->GetRatio(),fThr->Eval(centrValue),25000);
  line1->Draw("same");

  TLine *line01 = new TLine(fIntThr->Eval(centrValue),fIntThr->Eval(centrValue)*task->GetRatio(),25000,fIntThr->Eval(centrValue)*task->GetRatio());
  line01->Draw("same");
  TLine *line11 = new TLine(fIntThr->Eval(centrValue),fIntThr->Eval(centrValue)*task->GetRatio(),fIntThr->Eval(centrValue),25000);
  line11->Draw("same");

  TLine *line2 = new TLine(fThr->Eval(semiCentrValue),fThr->Eval(semiCentrValue)*task->GetRatio(),25000,fThr->Eval(semiCentrValue)*task->GetRatio());
  line2->Draw("same");
  TLine *line3 = new TLine(fThr->Eval(semiCentrValue),fThr->Eval(semiCentrValue)*task->GetRatio(),fThr->Eval(semiCentrValue),25000);
  line3->Draw("same");

  printf("Slope param: %f %f\n",f1->GetParameter(0),f2->GetParameter(0));
  printf("Mean slope param: %f\n",0.5*(f1->GetParameter(0)+f2->GetParameter(0)));
  printf("Delta slope param: %f %%\n",100.*(f1->GetParameter(0)-f2->GetParameter(0))/(0.5*(f1->GetParameter(0)+f2->GetParameter(0))));

  printf("A&C @ %.1f %%   ThrA = %f   ThrC = %f\n",centrValue,fThr->Eval(centrValue),fThr->Eval(centrValue)*task->GetRatio());
  printf("A&C @ %.1f %%   ThrA = %f   ThrC = %f\n",semiCentrValue,fThr->Eval(semiCentrValue),fThr->Eval(semiCentrValue)*task->GetRatio());

  printf("A&C (Integral) @ %.1f %%   ThrA = %f   ThrC = %f\n",centrValue,fIntThr->Eval(centrValue),fIntThr->Eval(centrValue)*task->GetRatio());
  printf("A&C (Intergal) @ %.1f %%   ThrA = %f   ThrC = %f\n",semiCentrValue,fIntThr->Eval(semiCentrValue),fIntThr->Eval(semiCentrValue)*task->GetRatio());

  printf("\n\n");
  FILE *fout=fopen("./trigger.txt","w");
  if (!fout) {
    printf("Failed to open local result file\n");
    return;
  }
  printf("%.0f %.0f %f %d %.0f %.0f %.0f %.0f\n\n",
	 task->GetMinThr(),
	 task->GetMaxThr(),
	 0.5*(f1->GetParameter(0)+f2->GetParameter(0)),
	 task->GetNThr(),
	 fThr->Eval(semiCentrValue),fThr->Eval(semiCentrValue)*task->GetRatio(),
	 fIntThr->Eval(centrValue),fIntThr->Eval(centrValue)*task->GetRatio());
  fprintf(fout,"%.0f %.0f %f %d %.0f %.0f %.0f %.0f\n",
	  task->GetMinThr(),
	  task->GetMaxThr(),
	  0.5*(f1->GetParameter(0)+f2->GetParameter(0)),
	  task->GetNThr(),
	  fThr->Eval(semiCentrValue),fThr->Eval(semiCentrValue)*task->GetRatio(),
	  fIntThr->Eval(centrValue),fIntThr->Eval(centrValue)*task->GetRatio());
  fclose(fout);
}
