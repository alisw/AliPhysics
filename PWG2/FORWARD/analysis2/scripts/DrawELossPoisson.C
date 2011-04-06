Double_t
DrawRingELossPoisson(TList* p, UShort_t d, Char_t r)
{
  if (!p) return;

  TList* ring = static_cast<TList*>(p->FindObject(Form("FMD%d%c",d,r)));
  if (!ring) { 
    Error("DrawELossPoisson", "List FMD%d%c not found in %s",d,r,p->GetName());
    return;
  }
  
  TH2* corr = static_cast<TH2D*>(ring->FindObject("elossVsPoisson"));
  if (!corr) { 
    Error("DrawRingELossPoisson", "Histogram esdEloss not found in FMD%d%c",
	  d, r);
    return;
  }
  gPad->SetGridy();
  gPad->SetGridx();
  gPad->SetLogz();
  gPad->SetFillColor(0);
  corr->SetTitle(Form("FMD%d%c",d,r));
  corr->Draw("colz");

  TLine* l = new TLine(0,0,corr->GetXaxis()->GetXmax(),
		       corr->GetYaxis()->GetXmax());
  l->SetLineWidth(1);
  l->SetLineColor(kBlack);
  l->Draw();
  // corr->GetXaxis()->SetRangeUser(0, 2);

  Info("", "FMD%d%c correlation coefficient: %9.5f%%", 
       d, r, 100*corr->GetCorrelationFactor());
  gPad->cd();
  return corr->GetCorrelationFactor();
}


void
DrawELossPoisson(const char* filename="forward.root")
{
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleW(.4);
  gStyle->SetTitleH(.1);
  gStyle->SetTitleColor(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(.6);
  
  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("DrawELossPoisson", "failed to open %s", filename);
    return;
  }

  TList* forward = static_cast<TList*>(file->Get("Forward"));
  if (!forward) { 
    Error("DrawELossPoisson", "List Forward not found in %s", filename);
    return;
  }

  TList* dc = static_cast<TList*>(forward->FindObject("fmdDensityCalculator"));
  if (!dc) { 
    Error("DrawELossPoisson", "List fmdDensityCalculator not found in Forward");
    return;
  }
  
  TCanvas* c = new TCanvas("elossVsPoisson", 
			   "N_ch from ELoss versus from Poisson", 900, 700);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->SetHighLightColor(0);
  c->Divide(3, 2, 0, 0);

  Double_t corrs[5];
  c->cd(1); corrs[0] = DrawRingELossPoisson(dc, 1, 'I');
  c->cd(2); corrs[1] = DrawRingELossPoisson(dc, 2, 'I');
  c->cd(5); corrs[2] = DrawRingELossPoisson(dc, 2, 'O');
  c->cd(3); corrs[3] = DrawRingELossPoisson(dc, 3, 'I');
  c->cd(6); corrs[4] = DrawRingELossPoisson(dc, 3, 'O');

  TVirtualPad* p = c->cd(4);
  p->SetTopMargin(0.05);
  p->SetRightMargin(0.10);
  p->SetLeftMargin(0.15);
  p->SetBottomMargin(0.15);
  p->SetFillColor(0);

  TH1D* hc = new TH1D("corrs", "Correlation factors", 5, .5, 5.5);
  hc->SetFillColor(kRed+1);
  hc->SetFillStyle(3001);
  hc->SetMinimum(0.0);
  hc->SetMaximum(1.1);
  hc->GetXaxis()->SetBinLabel(1,"FMD1i"); hc->SetBinContent(1,corrs[0]);
  hc->GetXaxis()->SetBinLabel(2,"FMD2i"); hc->SetBinContent(2,corrs[1]);
  hc->GetXaxis()->SetBinLabel(3,"FMD2o"); hc->SetBinContent(3,corrs[2]);
  hc->GetXaxis()->SetBinLabel(4,"FMD3i"); hc->SetBinContent(4,corrs[3]);
  hc->GetXaxis()->SetBinLabel(5,"FMD3o"); hc->SetBinContent(5,corrs[4]);
  hc->GetXaxis()->SetLabelSize(0.08);
  hc->GetYaxis()->SetTitle("r");
  hc->SetMarkerSize(1.5);
  hc->Draw("text hist");

  // TH2D* highCuts = static_cast<TH2D*>(dc->FindObject("highCuts"));
  // if (highCuts) highCuts->Draw("colz");
  c->cd();
  
}

  
  
 
