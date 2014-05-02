void
TestSPD(const TString& which, Double_t nVar=2)
{
  TFile* file = TFile::Open("forward.root", "READ");
  if (!file) return;

  Bool_t spd = which.EqualTo("spd", TString::kIgnoreCase);
  
  TList* l = 0;
  if (spd) l = static_cast<TList*>(file->Get("CentralSums"));
  else     l = static_cast<TList*>(file->Get("ForwardSums"));
  if (!l) { 
    Warning("", "%sSums not found", spd ? "Central" : "Forward");
    return;
  }

  TList* ei = static_cast<TList*>(l->FindObject("fmdEventInspector"));
  if (!l) { 
    Warning("", "fmdEventInspector not found");
    return;
  }
  
  TObject* run = ei->FindObject("runNo");
  if (!run) 
    Warning("", "No run number found");
  ULong_t runNo = run ? run->GetUniqueID() : 0;

  TH2* h = 0;
  if (spd) h = static_cast<TH2*>(l->FindObject("nClusterVsnTracklet"));
  else { 
    TList* den = static_cast<TList*>(l->FindObject("fmdDensityCalculator"));
    if (!den) { 
      Error("", "fmdDensityCalculator not found");
      return;
    }
    TList* rng = static_cast<TList*>(den->FindObject(which));
    if (!rng) { 
      Error("", "%s not found", which.Data());
      return;
    }
    h = static_cast<TH2*>(rng->FindObject("elossVsPoisson"));
  }
  if (!h) { 
    Warning("", "%s not found", spd ? nClusterVsnTracklet : "elossVsPoisson");
    return;
  }

  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c", Form("Run %u", runNo));
  c->Divide(2,2);
  
  TVirtualPad* p = c->cd(1);
  if (spd) {
    p->SetLogx();
    p->SetLogy();
  }
  p->SetLogz();
  h->Draw("colz");

  TObjArray* fits = new TObjArray;
  h->FitSlicesY(0, 1, -1, 0, "QN", fits);

  TF1* mean = new TF1("mean", "pol1");
  TF1* var  = new TF1("var", "pol1");
  // mean->FixParameter(0, 0);
  // var->FixParameter(0, 0);
  for (Int_t i = 0; i < 3; i++) { 
    p = c->cd(2+i);
    if (spd) { 
      p->SetLogx();
      p->SetLogy();
    }
    TH1* hh = static_cast<TH1*>(fits->At(i));
    hh->Draw();

    if (i == 0) continue;
    
    hh->Fit((i == 1? mean : var), "+Q");
    
  }

  TGraphErrors* g1 = new TGraphErrors(h->GetNbinsX());
  g1->SetFillColor(kBlue-10);
  g1->SetFillStyle(3001);
  g1->SetLineStyle(1);
  TGraph* u1 = new TGraph(h->GetNbinsX());
  TGraph* l1 = new TGraph(h->GetNbinsX());
  u1->SetLineColor(kBlue+1);
  l1->SetLineColor(kBlue+1);
  u1->SetName("u1");
  l1->SetName("l1");
  TGraphErrors* g2 = new TGraphErrors(h->GetNbinsX());
  g2->SetFillColor(kRed-10);
  g2->SetFillStyle(3001);
  g2->SetLineStyle(2);
  TGraph* u2 = new TGraph(h->GetNbinsX());
  TGraph* l2 = new TGraph(h->GetNbinsX());
  u2->SetLineColor(kRed+1);
  l2->SetLineColor(kRed+1);
  u2->SetName("u2");
  l2->SetName("l2");
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    Double_t x  = hh->GetXaxis()->GetBinCenter(i);
    Double_t y  = mean->Eval(x);
    Double_t e  = var->Eval(y);
    Double_t e1 = nVar * e;
    if (spd) e1 *= TMath::Log10(e);
    // Printf("%10e -> %10e +/- %10e", x, y, ee);
    g1->SetPoint(i-1, x, y);
    g1->SetPointError(i-1, 0, e1);
    u1->SetPoint(i-1, x, y+e1);
    l1->SetPoint(i-1, x, y-e1);
    // Printf("%3d: %f -> %f +/- %f", i, x, y, ee);

    Double_t e2 = nVar*0.05*x;
    g2->SetPoint(i-1, x, x);
    g2->SetPointError(i-1, 0, e2);
    u2->SetPoint(i-1, x, x+e2);
    l2->SetPoint(i-1, x, x-e2);
  }

  p = c->cd(1);
  c->Clear();
  c->cd();
  c->SetLogz();
  h->Draw("colz");
  g1->Draw("3 same");
  u1->Draw("l same");
  l1->Draw("l same");
  g2->Draw("3 same");
  u2->Draw("l same");
  l2->Draw("l same");

  Double_t ly = 0.9;
  Double_t dy = 0.06;
  TLatex* ltx = new TLatex(0.15, ly, Form("#LTy#GT = %f + %f x",
					   mean->GetParameter(0),
					   mean->GetParameter(1)));
  ltx->SetNDC();
  ltx->SetTextSize(dy);
  ltx->SetTextAlign(13);
  ltx->Draw();

  ly -= dy + 0.01;
  ltx->DrawLatex(0.15, ly, Form("#sigma_{y} = %f + %f x", 
				var->GetParameter(0),
				var->GetParameter(1)));
  
  ly -= dy + 0.01;
  ltx->DrawLatex(0.15, ly, Form("#delta = %f #sigma %s", 
				nVar, (spd ? "log_{10}(#sigma" : "")));
	    
					   
}

   
  
