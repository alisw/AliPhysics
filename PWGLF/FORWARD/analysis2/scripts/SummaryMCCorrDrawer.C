#include "SummaryDrawer.C"

/**
 * Class to draw a summary of the AOD production
 *
 * @par Input: 
 * - The merged <tt>forward.root</tt> file.
 *   If the file isn't merged, it should still work. 
 *
 * @par Output:
 * - A PDF file named after the input, but with <tt>.root</tt>
 *   replaced with <tt>pdf</tt>
 * 
 */
class SummaryMCCorrDrawer : public SummaryDrawer
{
public:
  enum EFlags { 
    kEventInspector    = 0x001, 
    kTrackDensity      = 0x002,
    kVertexBins        = 0x004,
    kResults           = 0x008, 
    kCentral           = 0x010,
    kNormal            = 0x01F
  };
  SummaryMCCorrDrawer() 
    : SummaryDrawer(),
      fSums(0),
      fResults(0)
  {}
  virtual ~SummaryMCCorrDrawer() {}
  
  //__________________________________________________________________
  /** 
   * 
   * 
   * @param fname 
   * @param what 
   */
  void Run(const char* fname, UShort_t what=kNormal)
  {
    // --- Open the file ---------------------------------------------
    TString filename(fname);
    TFile*  file = TFile::Open(filename, "READ");
    if (!file) { 
      Error("Run", "Failed to open \"%s\"", filename.Data());
      return;
    }


    // --- Get top-level collection ----------------------------------
    fSums = GetCollection(file, "ForwardCorrSums");
    if (!fSums) return;

    // --- Make our canvas -------------------------------------------
    TString pdfName(filename);
    pdfName.ReplaceAll(".root", ".pdf");
    CreateCanvas(pdfName, what & kLandscape);

    // --- Make a Title page -------------------------------------------
    DrawTitlePage(file);
    
    // --- Possibly make a chapter here ------------------------------
    if (what & kCentral && GetCollection(file, "CentralCorrSums")) 
      MakeChapter("Forward");
    
    // --- Set pause flag --------------------------------------------
    fPause = what & kPause;

    // --- Do each sub-algorithm -------------------------------------
    if (what & kEventInspector) DrawEventInspector(fSums);
    if (what & kTrackDensity)   DrawTrackDensity(fSums);
    if (what & kVertexBins)     DrawVertexBins(true);
  
    // --- Do the results ----------------------------------------------
    fResults = GetCollection(file, "ForwardCorrResults");
    if (!fResults) fResults = fSums; // Old-style
    
    if (what & kResults) DrawResults(true);

    // --- SPD clusters ----------------------------------------------
    if (what & kCentral) { 
      // --- Get top-level collection --------------------------------
      fSums = GetCollection(file, "CentralCorrSums");
      if (fSums) {
	MakeChapter("Central");
	if (what & kEventInspector) DrawEventInspector(fSums);
	if (what & kTrackDensity)   DrawTrackDensity(fSums);
	if (what & kVertexBins)     DrawVertexBins(false);
      }
      else 
	Warning("", "No CentralCorrSums found");

      fResults = GetCollection(file, "CentralCorrResults");
      if (!fResults) fResults = fSums; // Old-style
      if (!fResults) 
	Warning("", "No CentralCorrResults found");

      if (what & kResults) DrawResults(false);
    }

    CloseCanvas();
  }
protected:
  //____________________________________________________________________
  void DrawTitlePage(TFile* file)
  {
    TCollection* c   = GetCollection(file, "ForwardCorrSums");

    fBody->cd();
    
    Double_t y = .9;
    TLatex* ltx = new TLatex(.5, y, "ESD+MC #rightarrow Corrections");
    ltx->SetTextSize(0.07);
    ltx->SetTextFont(62);
    ltx->SetTextAlign(22);
    ltx->SetNDC();
    ltx->Draw();
    y -= .075;

    TCollection* ei = GetCollection(c, "fmdEventInspector");
    if (ei) { 
      UShort_t sys=0;
      UShort_t sNN=0;
      Int_t field=0;
      ULong_t runNo=0;
      GetParameter(ei, "sys", sys);
      GetParameter(ei, "sNN", sNN);
      GetParameter(ei, "field", field);
      GetParameter(ei, "runNo", runNo);

      TString tS; SysString(sys, tS); DrawParameter(y, "System", tS);
      TString tE; SNNString(sNN, tE); DrawParameter(y, "#sqrt{s_{NN}}", tE);
      DrawParameter(y, "L3 B field", Form("%+2dkG", field));
      DrawParameter(y, "Run #", Form("%6lu", runNo));
    }

    PrintCanvas("MC Corrections");
  }
  //____________________________________________________________________
  TCollection* GetVertexList(TCollection* parent, const TAxis& axis, Int_t bin)
  {
    TString folder = TString::Format("vtx%+05.1f_%+05.1f", 
				     axis.GetBinLowEdge(bin), 
				     axis.GetBinUpEdge(bin));
    folder.ReplaceAll(".", "d");
    folder.ReplaceAll("+", "p");
    folder.ReplaceAll("-", "m");
    
    TCollection* c = GetCollection(parent, folder);
    if (!c) 
      Warning("GetVertexList", "List %s not found", folder.Data());
	
    return c;
  }
  //____________________________________________________________________
  void DrawVertexBins(Bool_t forward)
  {
    TH1* vtxHist = GetH1(fSums, "vtxAxis");
    if (!vtxHist) return;
    
    Info("DrawVertexBins", "Drawing %s vertex bins - %d bins",
	 forward ? "Forward" : "Central", vtxHist->GetNbinsX());

    TAxis* vtxAxis = vtxHist->GetXaxis();
    for (Int_t i = 1; i <= vtxAxis->GetNbins(); i++) {
      Info("DrawVertexBins", " - Bin %d (%+5.1f - %+5.1f)", i,
	   vtxAxis->GetBinLowEdge(i), vtxAxis->GetBinUpEdge(i));
      TCollection* c = GetVertexList(fSums, *vtxAxis, i);
      if (!c) continue;
      
      if (forward) {
	DivideForRings(true, true);
	for (UShort_t d = 1; d <= 3; d++) { 
	  for (UShort_t q = 0; q < (d == 1 ? 1 : 2); q++) { 
	    Char_t r = q == 0 ? 'I' : 'O';
	    DrawInRingPad(d,r, GetH2(c, Form("FMD%d%c_cache",d,r)), "colz");
	  }
	}
	DrawInPad(fBody, (fLandscape ? 4: 2), GetH2(c, "primary"), "colz");
      }
      else {
	fBody->Divide(1,3, 0, 0);
	
	DrawInPad(fBody, 1, GetH2(c, "hits"),     "colz");
	DrawInPad(fBody, 2, GetH2(c, "clusters"), "colz");
	DrawInPad(fBody, 3, GetH2(c, "primary"), "colz");
      }
      PrintCanvas(Form("%s sums - IP_{z} bin %+5.1f - %+5.1f",
		       (forward ? "Forward" : "Central"),
		       vtxAxis->GetBinLowEdge(i), 
		       vtxAxis->GetBinUpEdge(i)));
    }
  }
  //____________________________________________________________________
  void DrawResults(Bool_t forward)
  {
    Info("DrawResults", "Drawing resulting %s vertex bins",
	 forward ? "Forward" : "Central");
    TH1* vtxHist = GetH1(fSums, "vtxAxis");
    if (!vtxHist) return;
    
    TAxis* vtxAxis = vtxHist->GetXaxis();
    for (Int_t i = 1; i <= vtxAxis->GetNbins(); i++) {
      Info("DrawResults", " - Bin %d (%+5.1f - %+5.1f)", i,
	   vtxAxis->GetBinLowEdge(i), vtxAxis->GetBinUpEdge(i));
      TCollection* c = GetVertexList(fResults, *vtxAxis, i);
      if (!c) continue;

      if (forward) {
	THStack* all = new THStack("all", 
				   "2^{nd} correction averaged over #phi");
	DivideForRings(true, true);
	for (UShort_t d = 1; d <= 3; d++) { 
	  for (UShort_t q = 0; q < (d == 1 ? 1 : 2); q++) { 
	    Char_t r = q == 0 ? 'I' : 'O';
	    // TVirtualPad* p = RingPad(d, r);
	    // p->cd();
	    // p->Divide(1,2,0,0);
	    TH2* h = GetH2(c, Form("FMD%d%c_vtxbin%03d",d,r,i));
	    DrawInRingPad(d,r, h,"colz");
	    // TVirtualPad* pp = p->cd(1);
	    // TVirtualPad* ppp = p->cd(2);
	    // ppp->SetRightMargin(pp->GetRightMargin());
	    TH1* hh = h->ProjectionX();
	    hh->Scale(1./ h->GetNbinsY());
	    hh->SetFillColor(RingColor(d,r));
	    hh->SetLineColor(RingColor(d,r));
	    hh->SetMarkerColor(RingColor(d,r));
	    hh->SetFillStyle(3001);
	    hh->SetTitle(Form("#LT%s#GT", hh->GetTitle()));
	    // DrawInPad(p,2, hh, "hist e");
	    all->Add(hh, "hist e");
	  }
	  TVirtualPad* p = RingPad(0, '0');
	  p->SetBottomMargin(0.10);
	  p->SetLeftMargin(0.10);
	  p->SetRightMargin(0.05);
	  DrawInRingPad(0, 'O', all, "nostack");
	}
      }
      else {
	fBody->Divide(2,4,0,0);
	TH2* secMap               = GetH2(c, "secMap");
	TH2* secMapAlt            = GetH2(c, "secMapEff");
	if (!secMapAlt) secMapAlt = GetH2(c, "secMapHit");
	TH1* acc                  = GetH1(c, "acc");
	TH1* accAlt               = GetH1(c, "accEff");
	if (!accAlt) accAlt       = GetH1(c, "accHit");

	TH1* secMapProj = secMap->ProjectionX();
	secMapProj->Scale(1./ secMap->GetNbinsY());
	secMapProj->Divide(acc);
	secMapProj->SetFillColor(kRed+1);
	secMapProj->SetFillStyle(3001);
	secMapProj->SetTitle(Form("#LT%s#GT/Acceptance", secMap->GetTitle()));

	TH1* secMapAltProj = secMapAlt->ProjectionX();
	secMapAltProj->Scale(1./ secMapAlt->GetNbinsY());
	secMapAltProj->Divide(accAlt);
	secMapAltProj->SetFillColor(kBlue+1);
	secMapAltProj->SetFillStyle(3001);
	secMapAltProj->SetTitle(Form("#LT%s#GT/Acceptance", 
				     secMapAlt->GetTitle()));

	Double_t secMapMax = TMath::Max(secMap->GetMaximum(),
					secMapAlt->GetMaximum());
	secMap->SetMaximum(secMapMax);
	secMapAlt->SetMaximum(secMapMax);

	Double_t secMapProjMax = TMath::Max(secMapProj->GetMaximum(),
					    secMapAltProj->GetMaximum());
	secMapProj->SetMaximum(secMapProjMax);
	secMapAltProj->SetMaximum(secMapProjMax);

	acc->SetFillColor(kRed+1);
	acc->SetFillStyle(3001);
	accAlt->SetFillColor(kBlue+1);
	accAlt->SetFillStyle(3001);

	Double_t accMax = TMath::Max(acc->GetMaximum(),accAlt->GetMaximum());
	acc->SetMaximum(accMax);
	accAlt->SetMaximum(accMax);

	if (secMap->GetMean(1) > 0 && secMap->GetMean(2) > 0)
	  DrawInPad(fBody, 1, secMap,    "colz", kGridx);
	if (secMapAlt->GetMean(1) > 0 && secMapAlt->GetMean(2) > 0)
	  DrawInPad(fBody, 2, secMapAlt, "colz", kGridx);

	TVirtualPad* p = fBody;
	TVirtualPad* pp = p->cd(1);
	TVirtualPad* ppp = p->cd(3);
	ppp->SetRightMargin(pp->GetRightMargin());
	DrawInPad(p,3, secMapProj, "hist", kGridx|kGridy);

	ppp = p->cd(5);
	ppp->SetRightMargin(pp->GetRightMargin());
	DrawInPad(fBody, 5, acc, "", kGridx|kGridy);

	pp = p->cd(2);
	pp->SetLeftMargin(0.10);
	ppp = p->cd(4);
	ppp->SetRightMargin(pp->GetRightMargin());
	ppp->SetLeftMargin(0.10);
	DrawInPad(p,4, secMapAltProj, "hist", kGridx|kGridy);

	pp = p->cd(4);
	pp->SetLeftMargin(0.10);
	ppp = p->cd(6);
	ppp->SetRightMargin(pp->GetRightMargin());
	ppp->SetLeftMargin(0.10);
	DrawInPad(fBody, 6, accAlt, "", kGridx|kGridy);

	TH2* diag = GetH2(c, "diagnostics");
	if (diag->GetMean(1) > 0 && diag->GetMean(2) > 0) 
	  DrawInPad(fBody, 7, diag, "colz");

	if (secMap->GetMean(1) > 0 && secMap->GetMean(2) > 0 &&
	    secMapAlt->GetMean(1) > 0 && secMapAlt->GetMean(2) > 0) {
	  TH2* ratio = static_cast<TH2*>(secMap->Clone("ratio"));
	  ratio->Add(secMapAlt, -1);
	  ratio->Divide(secMapAlt);
	  ratio->SetTitle("Relative difference between maps");
	  ratio->SetZTitle("#frac{S - S_{alt}}{S_{alt}}");
	
	  pp = p->cd(8);
	  pp->SetLeftMargin(0.10);
	  DrawInPad(fBody, 8, ratio, "colz");
	}
      }
      PrintCanvas(Form("%s results - IP_{z} bin %+5.1f - %+5.1f",
		       (forward ? "Forward" : "Central"),
		       vtxAxis->GetBinLowEdge(i), 
		       vtxAxis->GetBinUpEdge(i)));
    }
  }
  TCollection* fSums;
  TCollection* fResults;
};

// #endif
