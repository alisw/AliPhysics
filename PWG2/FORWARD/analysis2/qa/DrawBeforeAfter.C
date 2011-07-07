/** 
 * Draw the before/after merging image for a single ring
 * 
 * @param p 
 * @param d 
 * @param r 
 *
 * @ingroup pwg2_forward_scripts_qa
 */
void
DrawRingBeforeAfter(TList* p, UShort_t d, Char_t r)
{
  if (!p) return;

  TList* ring = static_cast<TList*>(p->FindObject(Form("FMD%d%c",d,r)));
  if (!ring) { 
    Error("DrawBeforeAfter", "List FMD%d%c not found in %s",d,r,p->GetName());
    return;
  }
  
  TH2* corr = static_cast<TH2D*>(ring->FindObject("beforeAfter"));
  if (!corr) { 
    Error("DrawRingBeforeAfter", "Histogram esdEloss not found in FMD%d%c",
	  d, r);
    return;
  }
  // gPad->SetLogz();
  gPad->SetFillColor(0);
  corr->SetTitle(Form("FMD%d%c",d,r));
  corr->Draw("colz");

  corr->GetXaxis()->SetRangeUser(-.5, 4);
  corr->GetYaxis()->SetRangeUser(-.5, 4);
  gPad->cd();
}


/** 
 * Draw the before/after sharing image for all rings 
 * 
 * @param filename 
 *
 * @ingroup pwg2_forward_scripts_qa
 */
void
DrawBeforeAfter(const char* filename="forward.root")
{
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleW(.4);
  gStyle->SetTitleH(.1);
  gStyle->SetTitleColor(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(.6);
  
  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("DrawBeforeAfter", "failed to open %s", filename);
    return;
  }

  TList* forward = static_cast<TList*>(file->Get("Forward"));
  if (!forward) { 
    Error("DrawBeforeAfter", "List Forward not found in %s", filename);
    return;
  }

  TList* sf = static_cast<TList*>(forward->FindObject("fmdSharingFilter"));
  if (!sf) { 
    Error("DrawBeforeAfter", "List fmdSharingFilter not found in Forward");
    return;
  }
  
  TCanvas* c = new TCanvas("beforeAfter", 
			   "Signals before and after merging", 900, 700);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.02);
  c->Divide(3, 2, 0, 0);
  
  c->cd(1); DrawRingBeforeAfter(sf, 1, 'I');
  c->cd(2); DrawRingBeforeAfter(sf, 2, 'I');
  c->cd(5); DrawRingBeforeAfter(sf, 2, 'O');
  c->cd(3); DrawRingBeforeAfter(sf, 3, 'I');
  c->cd(6); DrawRingBeforeAfter(sf, 3, 'O');
  TVirtualPad* p = c->cd(4);
  // p->SetTopMargin(0.05);
  p->SetRightMargin(0.15);
  p->SetFillColor(0);
  TH2D* highCuts = static_cast<TH2D*>(sf->FindObject("highCuts"));
  if (highCuts) highCuts->Draw("colz");
  c->cd();
  
}

  
  
 
//
// EOF
//
