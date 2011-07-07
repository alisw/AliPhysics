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

//____________________________________________________________________
void
CompareSecMaps(const char* fn1,   const char* fn2, 
	       const char* n1=0,  const char* n2=0,
	       bool load=true)
{

  // --- Load Utilities ----------------------------------------------
  if (load) {
    gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/CompareCorrs.C");
  }

  // --- Get Objects -------------------------------------------------
  TObject* o1 = 0;
  TObject* o2 = 0;
  GetObjects(AliForwardCorrectionManager::kSecondaryMap, fn1, fn2, o1, o2);
  if (!o1 || !o2) return; 
  AliFMDCorrSecondaryMap* obj1 = static_cast<AliFMDCorrSecondaryMap*>(o1);
  AliFMDCorrSecondaryMap* obj2 = static_cast<AliFMDCorrSecondaryMap*>(o2);
  UShort_t nVtx = obj1->GetVertexAxis().GetNbins();

  // --- Make canvas -------------------------------------------------
  Canvas* c = new Canvas("secMapComparison", "Ratio of secondary maps", n1, n2);
  c->Open();

  // --- Loop over the data ------------------------------------------
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nR = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nR; q++) { 
      Char_t   r  = (q == 0 ? 'I' : 'O');
      UShort_t nS = (q == 0 ?  20 :  40);

      // --- Make 2D ratios ------------------------------------------
      c->Clear(nVtx, d, r);
      TList hists;
      for (UShort_t v=1; v <= nVtx; v++) { 
	TVirtualPad* p = c->cd(v);
	
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

	if (ratio->GetMaximum()-ratio->GetMinimum() > 10) 
	  p->SetLogz();

	ratio->Draw("colz");
	hists.AddAt(ratio, v-1);
      }
      c->Print(d, r);
      
      // --- Make 1D profiles ----------------------------------------
      c->Clear(nVtx, d, r);
      for (UShort_t v=1; v <= nVtx; v++) { 
	c->cd(v);
	TH2* hist = static_cast<TH2*>(hists.At(v-1));
	TH1* prof = hist->ProjectionX();
	prof->Scale(1. / nS);
	prof->SetStats(0);
	prof->SetMinimum(0.8);
	prof->SetMaximum(1.2);

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

      c->Print(d, r, "profiles");
    }
  }

  // --- Close stuff -------------------------------------------------
  c->Close();
  // file1->Close();
  // file2->Close();
}

  
//____________________________________________________________________
//
// EOF
// 
