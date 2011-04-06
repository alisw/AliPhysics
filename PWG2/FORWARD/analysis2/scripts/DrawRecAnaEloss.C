void
DrawRingRecAnaEloss(TList* p, UShort_t d, Char_t r)
{
  if (!p) return;

  TList* ring = static_cast<TList*>(p->FindObject(Form("FMD%d%c",d,r)));
  if (!ring) { 
    Error("DrawRecAnaEloss", "List FMD%d%c not found in %s",d,r,p->GetName());
    return;
  }
  
  TH1* before = static_cast<TH2D*>(ring->FindObject("esdEloss"));
  if (!before) { 
    Error("DrawRingRecAnaEloss", "Histogram esdEloss not found in FMD%d%c",
	  d, r);
    return;
  }
  TH1* after = static_cast<TH2D*>(ring->FindObject("anaEloss"));
  if (!after) { 
    Error("DrawRingRecAnaEloss", "Histogram anaEloss not found in FMD%d%c",
	  d, r);
    return;
  }
  gPad->SetLogy();
  gPad->SetFillColor(0);
  before->SetTitle(Form("FMD%d%c",d,r));
  before->Draw("");
  after->Draw("same");

  before->GetXaxis()->SetRangeUser(0, 2);
  gPad->cd();
}


void
DrawRecAnaEloss(const char* filename="forward.root")
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
    Error("DrawRecAnaEloss", "failed to open %s", filename);
    return;
  }

  TList* forward = static_cast<TList*>(file->Get("Forward"));
  if (!forward) { 
    Error("DrawRecAnaEloss", "List Forward not found in %s", filename);
    return;
  }

  TList* sf = static_cast<TList*>(forward->FindObject("fmdSharingFilter"));
  if (!sf) { 
    Error("DrawRecAnaEloss", "List fmdSharingFilter not found in Forward");
    return;
  }
  
  TCanvas* c = new TCanvas("recAnaELoss", 
			   "Reconstructed and Analysed energy loss", 900, 700);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->Divide(3, 2, 0, 0);
  
  c->cd(1); DrawRingRecAnaEloss(sf, 1, 'I');
  c->cd(2); DrawRingRecAnaEloss(sf, 2, 'I');
  c->cd(5); DrawRingRecAnaEloss(sf, 2, 'O');
  c->cd(3); DrawRingRecAnaEloss(sf, 3, 'I');
  c->cd(6); DrawRingRecAnaEloss(sf, 3, 'O');
  TVirtualPad* p = c->cd(4);
  // p->SetTopMargin(0.05);
  p->SetRightMargin(0.15);
  p->SetFillColor(0);
  TH2D* highCuts = static_cast<TH2D*>(sf->FindObject("highCuts"));
  if (highCuts) highCuts->Draw("colz");
  c->cd();
  
}

  
  
 
