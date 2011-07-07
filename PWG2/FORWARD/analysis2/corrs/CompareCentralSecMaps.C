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
CompareCentralSecMaps(const char* fn1,   const char* fn2, 
		      const char* n1=0,  const char* n2=0,
		      bool load=true)
{

  // --- Load Utilities ----------------------------------------------
  if (load) {
    gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");
    //   gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/CompareCorrs.C");
  }

  // --- Get Objects -------------------------------------------------
  TObject* o1 = 0;
  TObject* o2 = 0;
  TFile file1(fn1);
  TFile file2(fn2);
  
  AliCentralMultiplicityTask task("tmp");
  
  o1 = file1.Get(task.GetManager().GetSecMapName());
  o2 = file2.Get(task.GetManager().GetSecMapName());
  if (!o1 || !o2) return; 
  AliCentralCorrSecondaryMap* obj1 = static_cast<AliCentralCorrSecondaryMap*>(o1);
  AliCentralCorrSecondaryMap* obj2 = static_cast<AliCentralCorrSecondaryMap*>(o2);
  UShort_t nVtx = obj1->GetVertexAxis().GetNbins();

  // --- Make canvas -------------------------------------------------
  // Canvas* c = new Canvas("CentralSecMapComparison", "Ratio of central secondary maps", n1, n2);
  // c->Open();

  // --- Loop over the data ------------------------------------------
  //for (UShort_t d = 1; d <= 3; d++) { 
  //  UShort_t nR = (d == 1 ? 1 : 2);
  //  for (UShort_t q = 0; q < nR; q++) { 
  //   Char_t   r  = (q == 0 ? 'I' : 'O');
  //    UShort_t nS = (q == 0 ?  20 :  40);

      // --- Make 2D ratios ------------------------------------------
      //c->Clear(nVtx, d, r);
  TCanvas* c2D = new TCanvas("2Dcomparison","2Dcomparison",800,1000);
  TCanvas* c1D = new TCanvas("1Dcomparison","1Dcomparison",800,1000);
  c2D->Divide(2,5);
  c1D->Divide(2,5);
  c2D->cd();
  
  TList hists;
  for (UShort_t v=1; v <= nVtx; v++) { 
    //	TVirtualPad* p = c->cd(v);
    c2D->cd(v);
    TH2* h1 = obj1->GetCorrection( v);
    TH2* h2 = obj2->GetCorrection( v);
    //std::cout<<h1<<"  "<<h2<<std::endl;
    /*	if (!h1) { 
	Error("CompareSecMaps", 
	"Map for SPD, vtxbin %3d not found in first", 
	v);
	continue;
	}
	if (!h2) { 
	Error("CompareSecMaps", 
	"Map for SPD, vtxbin %3d not found in second", 
	v);
	continue;
	}
    */
    Double_t vl    = obj1->GetVertexAxis().GetBinLowEdge(v);
    Double_t vh    = obj1->GetVertexAxis().GetBinUpEdge(v);
    TH2*     ratio = 
      static_cast<TH2*>(h1->Clone(Form("tmpSPD_%3d",v)));
    ratio->SetName(Form("SPD_vtx%03d_ratio",  v));
    ratio->SetTitle(Form("%+5.1f<v_{z}<%-+5.1f", vl, vh));
    ratio->Divide(h2);
    ratio->SetStats(0);
    ratio->SetDirectory(0);
    ratio->SetZTitle("ratio");
    
	//if (ratio->GetMaximum()-ratio->GetMinimum() > 10) 
	//  p->SetLogz();

    ratio->DrawCopy("colz");
    hists.AddAt(ratio, v-1);
    
  }
      // c->Print(d, r);
      
      // --- Make 1D profiles ----------------------------------------
      //c->Clear(nVtx, d, r);
    
  c1D->cd();
  for (UShort_t v=1; v <= nVtx; v++) { 
    //c->cd(v);
    c1D->cd(v);
    TH2* hist = static_cast<TH2*>(hists.At(v-1));
    TH1D* prof = hist->ProjectionX();
    prof->Clear();
    for(Int_t i=1; i<=hist->GetNbinsX();i++) {
      
      Float_t sum = 0;
      Float_t error = 0;
      Int_t   n   = 0;
      for(Int_t j=1; j<=hist->GetNbinsY();j++) {
	if(hist->GetBinContent(i,j) > 0) {
	  sum = sum+hist->GetBinContent(i,j);
	  error = error + TMath::Power(hist->GetBinError(i,j),2);
	  n++;
	}
	
      }
      if(n>0) {
	sum = sum/(Float_t)n;
	error = TMath::Sqrt(error) / (Float_t)n;
	prof->SetBinContent(i,sum);
	prof->SetBinError(i,error);
      }
    }
    
    //    prof->Scale(0.05);
    prof->SetStats(0);
    prof->SetMinimum(0.8);
    prof->SetMaximum(1.2);
    
    
    prof->Fit("pol0","Q");
    prof->DrawCopy();
    
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
  
  //c->Print(d, r, "profiles");
  c2D->Print("comparisonSPD.pdf(");
  c1D->Print("comparisonSPD.pdf)");
  

  // --- Close stuff -------------------------------------------------
  //c->Close();
  // file1->Close();
  // file2->Close();
}

  
//____________________________________________________________________
//
// EOF
// 
