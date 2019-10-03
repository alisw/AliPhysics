//____________________________________________________________________
/** 
 * 
 * 
 * @param fn1 
 * @param fn2 
 * @param n1 
 * @param n2 
 * @param load 
 *
 * @ingroup pwglf_forward_scripts_corr
 */
void
CompareSecMaps(const char* fn1,   const char* fn2, 
	       const char* n1=0,  const char* n2=0,
	       bool load=true)
{

  // --- Load Utilities ----------------------------------------------
  if (load) {
    gROOT->Macro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/scripts/LoadLibs.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/corrs/CompareCorrs.C");
  }

  // --- Get Objects -------------------------------------------------
  TObject* o1 = 0;
  TObject* o2 = 0;
  GetObjects("AliFMDCorrSecondaryMap", fn1, fn2, o1, o2);
  if (!o1 || !o2) return; 
  AliFMDCorrSecondaryMap* obj1 = static_cast<AliFMDCorrSecondaryMap*>(o1);
  AliFMDCorrSecondaryMap* obj2 = static_cast<AliFMDCorrSecondaryMap*>(o2);
  UShort_t nVtx = obj1->GetVertexAxis().GetNbins();

  // --- Make canvas -------------------------------------------------
  Canvas* c = new Canvas("secMapComparison", "Ratio of secondary maps", n1, n2);
  c->Open();

  TObjArray* allVtx = new TObjArray(nVtx);
  // --- Loop over the data ------------------------------------------
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nR = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nR; q++) { 
      Char_t   r   = (q == 0 ? 'I' : 'O');
      UShort_t nS  = (q == 0 ?  20 :  40);
      Color_t  col = ((d == 1 ? kRed : (d == 2 ? kGreen : kBlue))
		      + ((r == 'I' || r == 'i') ? -3 : -9));

      // --- Make 2D ratios ------------------------------------------
      c->Clear(nVtx, d, r, true);
      TList hists;
      for (UShort_t v=1; v <= nVtx; v++) { 
	TVirtualPad* p = c->cd(v);
	p->SetGridx();
	p->SetGridy();
	// p->SetRightMargin(0.01);
	// p->SetTopMargin(0.01);
	
	TH2* h1 = obj1->GetCorrection(d, r, v);
	TH2* h2 = obj2->GetCorrection(d, r, v);
  
	if (!h1) { 
	  Error("CompareSecMaps", 
		"Map for FMD%d%c, vtxbin %3d not found in first", 
		d, r, v);
	  continue;
	}
	if (!h2) { 
	  Error("CompareSecMaps", 
		"Map for FMD%d%c, vtxbin %3d not found in second", 
		d, r, v);
	  continue;
	}
  
	Double_t vl    = obj1->GetVertexAxis().GetBinLowEdge(v);
	Double_t vh    = obj1->GetVertexAxis().GetBinUpEdge(v);
	TH2*     ratio = 
	  static_cast<TH2*>(h1->Clone(Form("tmpFMD%d%c_%3d",d,r,v)));
	ratio->SetName(Form("FMD%d%c_vtx%03d_ratio", d, r, v));
	ratio->SetTitle(Form("%+5.1f<v_{z}<%-+5.1f", vl, vh));
	ratio->Divide(h2);
	ratio->SetStats(0);
	ratio->SetDirectory(0);
	ratio->SetZTitle("ratio");
	ratio->GetXaxis()->SetTitleSize(0.01/p->GetHNDC());
	ratio->GetXaxis()->SetLabelSize(0.01/p->GetHNDC());
	ratio->GetYaxis()->SetTitleSize(0.01/p->GetHNDC());
	ratio->GetYaxis()->SetLabelSize(0.01/p->GetHNDC());
	
	if (ratio->GetMaximum()-ratio->GetMinimum() > 10) 
	  p->SetLogz();

	ratio->Draw("colz");
	hists.AddAt(ratio, v-1);
      }
      c->Print(d, r);
      
      // --- Make 1D profiles ----------------------------------------
      c->Clear(nVtx, d, r, true);
      for (UShort_t v=1; v <= nVtx; v++) { 
	TVirtualPad* p = c->cd(v);
	p->SetGridx();
	p->SetGridy();


	TH2* hist = static_cast<TH2*>(hists.At(v-1));
	TH1* prof = hist->ProjectionX();
	prof->SetDirectory(0);
	prof->Scale(1. / nS);
	// TProfile* prof = hist->ProfileX();
	prof->SetStats(0);
	prof->SetMinimum(0.8);
	prof->SetMaximum(1.2);
	prof->SetLineColor(col);
	prof->SetMarkerColor(col);
	prof->SetMarkerStyle(20);
	prof->GetXaxis()->SetTitleSize(0.01/p->GetHNDC());
	prof->GetXaxis()->SetLabelSize(0.01/p->GetHNDC());
	prof->GetYaxis()->SetTitleSize(0.01/p->GetHNDC());
	prof->GetYaxis()->SetLabelSize(0.01/p->GetHNDC());

	THStack* stack = static_cast<THStack*>(allVtx->At(v-1));
	if (!stack) {
	  stack = new THStack(prof->GetName(), prof->GetTitle());
	  stack->SetMinimum(prof->GetMinimum());
	  stack->SetMaximum(prof->GetMaximum());
	  allVtx->AddAt(stack, v-1);
	}
	stack->Add(static_cast<TH1*>(prof->Clone()));

	prof->Draw();
	prof->Fit("pol0","Q");

	TF1* f = prof->GetFunction("pol0");
	if (!f) continue;

	f->SetLineColor(kBlack);
	Double_t chi2nu = (f->GetNDF() != 0 ?
			   f->GetChisquare() / f->GetNDF() : 0);
	Double_t dist = TMath::Abs(1 - f->GetParameter(0));

	Double_t ly = .8;
	TLatex* l = new TLatex(0.5, ly,
			       Form("A = %f #pm %f",
				    f->GetParameter(0), f->GetParError(0)));
	l->SetTextAlign(22);
	l->SetNDC();
	l->Draw();
	ly -= 1.1*l->GetTextSize();
	l->DrawLatex(0.5, ly, Form("|1 - A| = %f %s #deltaA", 
				   dist, dist <= f->GetParError(0) ? 
				   "#leq" : ">")); 
	ly -= 1.1*l->GetTextSize();
	l->DrawLatex(0.5, ly, Form("#chi^2/NDF = %f / %d = %f", 
				   f->GetChisquare(), f->GetNDF(), chi2nu));

	TLine* l1 = new TLine(prof->GetXaxis()->GetXmin(), 1,
			      prof->GetXaxis()->GetXmax(), 1);
	l1->SetLineColor(kGray+2);
	l1->SetLineStyle(2);
	l1->Draw();

      }

      c->Print(d, r, "profiles");
    }
  }
  c->Clear(nVtx, 0, '*', true);
  for (UShort_t v=1; v <= nVtx; v++) { 
    TVirtualPad* p = c->cd(v);
    p->SetGridx();
    p->SetGridy();
    THStack* stack = static_cast<THStack*>(allVtx->At(v-1));
    if (!stack) continue;
    stack->Draw("nostack");
    TH1* h = stack->GetHistogram();
    if (!h) continue;
    TLine* l1 = new TLine(h->GetXaxis()->GetXmin(), 1,
			  h->GetXaxis()->GetXmax(), 1);
    l1->SetLineColor(kGray+2);
    l1->SetLineStyle(2);
    l1->Draw();
  }
  c->Print(0,'*',"summary");
    
  
  // --- Close stuff -------------------------------------------------
  c->Close();
  // file1->Close();
  // file2->Close();
}

  
//____________________________________________________________________
//
// EOF
// 
