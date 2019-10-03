void CompareCorrectedSpectra(const char *datafilea, const char *datafileb){
  
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);

  // Take files

  TFile *filea = TFile::Open(datafilea);
  TGraphErrors *correctedSpectrum_a = (TGraphErrors *) filea->Get("AlltogetherSpectrum");
 
  TH1D *histo = (TH1D*) filea->Get("RatioUnfoldingAlltogetherSpectrum");
  histo->SetName("historatio");
  TH1D *histoa = (TH1D*) histo->Clone();
  histoa->Sumw2();
  histoa->SetName("a");
  TH1D *histob = (TH1D*) histo->Clone();
  histob->Sumw2();
  histob->SetName("b");
  
 
  TFile *fileb = TFile::Open(datafileb);
  TGraphErrors *correctedSpectrum_b = (TGraphErrors *) fileb->Get("AlltogetherSpectrum");
    
  // Style
  
 correctedSpectrum_a->SetMarkerStyle(24);
 correctedSpectrum_a->SetMarkerColor(1);
 correctedSpectrum_a->SetLineColor(1);

 correctedSpectrum_b->SetMarkerStyle(27);
 correctedSpectrum_b->SetMarkerColor(4);
 correctedSpectrum_b->SetLineColor(4);

  //

  TCanvas *c1 = new TCanvas("CorrectedSpectrum","CorrectedSpectrum",800,800);
  c1->cd(1);
  TH1D *total = new TH1D("total","",1,0.38,4.3);
  total->SetMaximum(1.0);
  total->SetMinimum(1.0e-09);
  total->SetXTitle("p_{T} [GeV/c]");
  total->SetYTitle("1/2#pip_{T} d^{2}N/dp_{T}dy [GeV/c]^{-2}, |#eta| < 0.8");
  total->SetTitleOffset(1.5,"Y");
  total->Draw();
  gPad->SetLeftMargin(0.13);
  gPad->SetLogy();
  gPad->SetTicks();
  gPad->SetFillColor(10);
  gPad->SetFrameFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);
  correctedSpectrum_a->Draw("Psame");
  correctedSpectrum_b->Draw("Psame");
  TLegend *leg = new TLegend(0.25,0.825,0.48,0.875);
  leg->AddEntry(correctedSpectrum_a,"Corrected spectrum a","lep");
  leg->AddEntry(correctedSpectrum_b,"Corrected spectrum b","lep");
  leg->SetFillStyle(0);
  leg->Draw("same");
    
  // Ratio

  Double_t x[300];
  Double_t ry[300];
  Double_t y[300];
  Double_t rex[300];
  Double_t rey[300];

  double xa,ya,xb,yb,eya,exa,eyb,exb;
  Int_t npointsa = correctedSpectrum_a->GetN();
  Int_t npointsb = correctedSpectrum_b->GetN();
  if(npointsa != npointsb) {
    printf("Problem the two spectra have not the same number of points");
    return;
  }
  for(Int_t k = 0; k < npointsa; k++){
    correctedSpectrum_a->GetPoint(k,xa,ya);
    correctedSpectrum_b->GetPoint(k,xb,yb);
    //
    Double_t centerhisto = histoa->GetBinCenter(k+1);
    //printf("bin center %f and center %f\n",centerhisto,xa);
    histoa->SetBinContent(k+1,ya);
    histob->SetBinContent(k+1,yb);
    //
    if(TMath::Abs(xa-xb) > 0.0001) {
      printf("Problem the two spectra have not the same number of points");
      return;
    }
    eya = correctedSpectrum_a->GetErrorY(k);
    exa = correctedSpectrum_a->GetErrorX(k);
    eyb = correctedSpectrum_b->GetErrorY(k);
    exb = correctedSpectrum_b->GetErrorX(k);
    x[k] = xa;
    rex[k] = exa;
    if(yb > 0.0) y[k] = ya/yb;
    ry[k] = ya-yb;
    Double_t error = 0.0;
    if((yb > 0.0) && (ya > 0.0)) {
      error = TMath::Sqrt(TMath::Abs(eya*eya-eyb*eyb));
    }
    rey[k] = error;
    histoa->SetBinError(k+1,eya);
    histob->SetBinError(k+1,eyb);
  }
  
  histo->Sumw2();
  histo->Divide(histoa,histob,1.0,1.0,"B");

  TGraph *gratio = new TGraph(npointsa,&x[0],&y[0]);
  gratio->SetName("RatioOfSpectra");

  TGraphErrors *gdiff = new TGraphErrors(npointsa,&x[0],&ry[0],&rex[0],&rey[0]);
  gdiff->SetName("Difference");

  TCanvas *c2 = new TCanvas("ratio","ratio",800,800);
  c2->cd(1);
  TH1D *ratio = new TH1D("ratio","",1,0.38,10.0);
  ratio->SetMaximum(1.5);
  ratio->SetMinimum(0.5);
  ratio->SetXTitle("p_{T} [GeV/c]");
  ratio->SetYTitle("Ratio of corrected spectra a/b");
  ratio->SetTitleOffset(1.5,"Y");
  //ratio->Draw();
  histo->SetLineColor(2);
  histo->Draw();
  gPad->SetLeftMargin(0.13);
  gPad->SetTicks();
  //gPad->SetGridx();
  //gPad->SetGridy();
  gPad->SetFillColor(10);
  gPad->SetFrameFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);
  gratio->SetMarkerStyle(24);
  gratio->SetMarkerColor(4);
  gratio->SetLineColor(4);
  gratio->Draw("P");

  TCanvas *c3 = new TCanvas("Difference","Difference",800,800);
  c3->cd(1);
  TH1D *diff = new TH1D("difference","",1,0.38,4.3);
  diff->SetMaximum(1.0e-04);
  diff->SetMinimum(-1.0e-04);
  diff->SetXTitle("p_{T} [GeV/c]");
  diff->SetYTitle("Difference of corrected spectra a-b");
  diff->SetTitleOffset(1.5,"Y");
  diff->Draw();
  gPad->SetLeftMargin(0.13);
  gPad->SetTicks();
  //gPad->SetGridx();
  //gPad->SetGridy();
  gPad->SetFillColor(10);
  gPad->SetFrameFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);
  gdiff->SetMarkerStyle(24);
  gdiff->SetMarkerColor(4);
  gdiff->SetLineColor(4);
  gdiff->Draw("P");
  

 return;

}

