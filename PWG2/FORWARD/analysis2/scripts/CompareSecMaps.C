/** 
 * Make ratio of two specific maps 
 * 
 * @param d        Detector
 * @param r        Ring
 * @param v        Vertex bin (1 based)
 * @param first    First correction
 * @param second   Second correction
 * 
 * @return Ratio of the two, or null
 */
TH2*
Compare2Maps(UShort_t d, Char_t r, UShort_t v, 
	     const AliFMDCorrSecondaryMap& first, 
	     const AliFMDCorrSecondaryMap& second)
{
  TH2* h1 = first.GetCorrection(d, r, v);
  TH2* h2 = second.GetCorrection(d, r, v);
  
  if (!h1) { 
    Error("Compare2Maps", "Map for FMD%d%c, vtxbin %3d not found in first", 
	  d, r, v);
    return 0;
  }
  if (!h1) { 
    Error("Compare2Maps", "Map for FMD%d%c, vtxbin %3d not found in second", 
	  d, r, v);
    return 0;
  }
  
  Double_t vl    = first.GetVertexAxis().GetBinLowEdge(v);
  Double_t vh    = first.GetVertexAxis().GetBinUpEdge(v);
  TH2*     ratio = static_cast<TH2*>(h1->Clone(Form("tmpFMD%d%c_%3d",d,r,v)));
  ratio->SetName(Form("FMD%d%c_vtx%03d_ratio", d, r, v));
  ratio->SetTitle(Form("%+5.1f<v_{z}<%-+5.1f", vl, vh));
  ratio->Divide(h2);
  ratio->SetStats(0);
  ratio->SetDirectory(0);
  ratio->SetZTitle("ratio");
  // ratio->SetMinimum(0.9);
  // ratio->SetMaximum(1.5);

  return ratio;
}

//____________________________________________________________________
TPad* 
ClearCanvas(TCanvas* c, UShort_t nVtx, UShort_t d, Char_t r, 
	    const char* n1, const char* n2)
{
  c->Clear();
  TPad* p1 = new TPad("top", "Top", 0, .95, 1, 1, 0, 0);
  p1->Draw();
  p1->cd();

  TLatex* l = new TLatex(.5, .5, Form("Ratio of secondary maps for "
				       "FMD%d%c (%s / %s)", d, r, n1, n2));
  l->SetNDC();
  l->SetTextAlign(22);
  l->SetTextSize(0.3);
  l->Draw();
  
  c->cd();
  TPad* body = new TPad("body", "Body", 0, 0, 1, .95, 0, 0);
  body->SetTopMargin(0.05);
  body->SetRightMargin(0.05);
  body->Divide(2, (nVtx+1)/2, 0.001, 0.001);
  body->Draw();

  return body;
}

//____________________________________________________________________
void
CompareSecMaps(const char* fn1,   const char* fn2, 
	       const char* n1=0,  const char* n2=0)
{

  // --- Load libraries ----------------------------------------------
  // gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");
  
  // --- Open files --------------------------------------------------
  const char* nam1 = n1;
  const char* nam2 = n2;
  if (!n1) nam1 = fn1;				
  if (!n2) nam2 = fn2;

  TFile* file1 = TFile::Open(fn1, "READ");
  TFile* file2 = TFile::Open(fn2, "READ");

  if (!file1) { 
    Error("CompareSecMaps", "File %s cannot be opened for %s", fn1, n1);
    return;
  }

  if (!file2) { 
    Error("CompareSecMaps", "File %s cannot be opened for %s", fn2, n2);
    return;
  }
  
  // --- Find Objects ------------------------------------------------
  const char* objName = AliForwardCorrectionManager::Instance()
    .GetObjectName(AliForwardCorrectionManager::kSecondaryMap);
  
  AliFMDCorrSecondaryMap* obj1 = 
    static_cast<AliFMDCorrSecondaryMap*>(file1->Get(objName));
  AliFMDCorrSecondaryMap* obj2 = 
    static_cast<AliFMDCorrSecondaryMap*>(file2->Get(objName));
  
  if (!obj1) {
    Error("CompareSecMaps", "File %s does not contain an object named %s", 
	  fn1, objName);
    return;
  }
  if (!obj2) {
    Error("CompareSecMaps", "File %s does not contain an object named %s", 
	  fn2, objName);
    return;
  }
  UShort_t nVtx = obj1->GetVertexAxis().GetNbins();

  // --- Make canvas -------------------------------------------------
  const char* pdfName = "secMapComparison.pdf";
  gStyle->SetPalette(1);
  gStyle->SetTitleX(.10);
  gStyle->SetTitleY(.99);
  gStyle->SetTitleW(.85);
  gStyle->SetTitleH(.085);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);

  TCanvas* c = new TCanvas("c", "c", 800, TMath::Sqrt(2)*800);
  c->SetFillColor(0);
  
  c->Print(Form("%s[", pdfName), "pdf");

  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nR = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nR; q++) { 
      Char_t   r  = (q == 0 ? 'I' : 'O');
      UShort_t nS = (q == 0 ?  20 :  40);

      TPad* body = ClearCanvas(c, nVtx, d, r, nam1, nam2);
      TList hists;
      for (UShort_t v=1; v <= nVtx; v++) { 
	TVirtualPad* p = body->cd(v);
	// p->SetTopMargin(0.1);
	// p->SetBottomMargin(0.05);
	// p->SetRightMargin(0.05);
	
	TH2* ratio = Compare2Maps(d, r, v, *obj1, *obj2);
	if (ratio->GetMaximum()-ratio->GetMinimum() > 10) 
	  p->SetLogz();

	ratio->Draw("colz");
	hists.AddAt(ratio, v-1);
      }
      c->Print(pdfName, Form("Title:FMD%d%c", d, r));
      
      body = ClearCanvas(c, nVtx, d, r, nam1, nam2);

      for (UShort_t v=1; v <= nVtx; v++) { 
	TVirtualPad* p    = body->cd(v);
	TH2*         hist = static_cast<TH2*>(hists.At(v-1));
	TH1*         prof = hist->ProjectionX();
	prof->Scale(1. / nS);
	prof->SetStats(0);
	prof->SetMinimum(0.8);
	prof->SetMaximum(1.2);

	// prof->Draw("hist");
	// prof->DrawCopy("e same");
	prof->Draw();
	prof->Fit("pol0","Q");

	TF1* f = prof->GetFunction("pol0");

	TLatex* l = new TLatex(0.5, 0.4, Form("A = %f #pm %f", 
					      f->GetParameter(0), 
					      f->GetParError(0)));
	l->SetTextAlign(22);
	l->SetNDC();
	l->Draw();
	l->DrawLatex(0.5, 0.3, Form("#chi^2/NDF = %f / %d = %f", 
				    f->GetChisquare(), 
				    f->GetNDF(), 
				    f->GetChisquare() / f->GetNDF()));
	Double_t dist = TMath::Abs(1 - f->GetParameter(0));
	l->DrawLatex(0.5, 0.35, Form("|1 - A| = %f %s #deltaA", 
				     dist, dist <= f->GetParError(0) ? 
				     "#leq" : ">")); 

	TLine* l1 = new TLine(-4, 1, 6, 1);
	l1->SetLineColor(kRed);
	l1->SetLineStyle(2);
	l1->Draw();
      }

      c->Print(pdfName, Form("Title:FMD%d%c profiles", d, r));
    }
  }

  c->Print(Form("%s]", pdfName), "pdf");
  file1->Close();
  file2->Close();
}

  
