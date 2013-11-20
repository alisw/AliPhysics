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
class SummaryAODDrawer : public SummaryDrawer
{
public:
  enum EFlags { 
    kEventInspector    = 0x001, 
    kSharingFilter     = 0x002, 
    kDensityCalculator = 0x004,
    kCorrector         = 0x008,
    kHistCollector     = 0x010,
    kSteps             = 0x020, 
    kResults           = 0x040, 
    kCentral           = 0x080,
    kNormal            = 0x0FF
  };
  SummaryAODDrawer() 
    : SummaryDrawer(),
      fSums(0),
      fResults(0)
  {}
  
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
    fSums = GetCollection(file, "ForwardSums");
    if (!fSums) {
      Info("Run", "Trying old name Forward");
      fSums = GetCollection(file, "Forward");
      if (!fSums) return;
    }

    // --- Do the results ----------------------------------------------
    fResults = GetCollection(file, "ForwardResults");
    if (!fResults) fResults = fSums; // Old-style

    // --- Make our canvas -------------------------------------------
    TString pdfName(filename);
    pdfName.ReplaceAll(".root", ".pdf");
    CreateCanvas(pdfName, what & kLandscape);
    DrawTitlePage(file);

    // --- Possibly make a chapter here ------------------------------
    TCollection* centralSums = GetCollection(file, "CentralSums", false);
    if (!centralSums) 
      centralSums = GetCollection(file, "Central", false);
    if (what & kCentral && centralSums) 
      MakeChapter("Forward");
    
    // --- Set pause flag --------------------------------------------
    fPause = what & kPause;

    // --- Do each sub-algorithm -------------------------------------
    if (what & kEventInspector)    DrawEventInspector(fSums);
    if (what & kSharingFilter)     DrawSharingFilter();
    if (what & kDensityCalculator) DrawDensityCalculator();
    if (what & kCorrector)         DrawCorrector();
    if (what & kHistCollector)     DrawHistCollector();
  
    
    if (what & kSteps) DrawSteps();
    if (what & kResults) DrawResults();

    // --- SPD clusters ----------------------------------------------
    if (what & kCentral) { 
      // --- Get top-level collection --------------------------------
      fSums = GetCollection(file, "CentralSums");
      if (fSums) {
	MakeChapter("Central");
	DrawCentral();
	if (what & kEventInspector) DrawEventInspector(fSums);
      }
      fResults = GetCollection(file, "CentralResults");
      if (fResults && (what & kResults)) {
	DrawCentralResults();
      }

      if (what & kResults) DrawBoth(file);
    }

    
    CloseCanvas();
  }
protected:
  //____________________________________________________________________
  void DrawTitlePage(TFile* f)
  {
    fBody->cd();

    TLatex* ltx = new TLatex(.5, .7, "ESD #rightarrow AOD filtering");
    ltx->SetNDC();
    ltx->SetTextSize(0.07);
    ltx->SetTextAlign(22);
    ltx->Draw();

    TCollection* fwd = fSums; // GetCollection(f, "ForwardSums");
    TCollection* cen = GetCollection(f, "CentralSums");
    Double_t y = .6;
    
    Double_t save = fParName->GetTextSize();
    fParName->SetTextSize(0.03);
    fParVal->SetTextSize(0.03);

    DrawParameter(y, "Tasks", (fwd ? "Forward" : ""));
    DrawParameter(y, "",      (cen ? "Central" : ""));

    if (fwd) { 
      TCollection* ei = GetCollection(fwd, "fmdEventInspector");
      if (ei) { 

	UShort_t sys=0, sNN=0;
	Int_t field=0;
	ULong_t runNo=0;
	Bool_t mc=false;
	GetParameter(ei, "sys", sys);
	GetParameter(ei, "sNN", sNN);
	GetParameter(ei, "field", field);
	GetParameter(ei, "runNo", runNo);
	if (!GetParameter(ei, "mc", mc, false)) mc = false;
	
	TString sysString;    SysString(sys, sysString);
	TString sNNString;    SNNString(sNN, sNNString);
	
	DrawParameter(y, "System", sysString);
	DrawParameter(y, "#sqrt{s_{NN}}", sNNString);
	DrawParameter(y, "L3 B field", Form("%+2dkG", field));
	DrawParameter(y, "Run #", Form("%6lu", runNo));
	DrawParameter(y, "Simulation", (mc ? "yes" : "no"));	
      }
    }
    PrintCanvas("Title page");
    fParName->SetTextSize(save);
    fParVal->SetTextSize(save);
  }
  //____________________________________________________________________
  void DrawSharingFilter()
  {
    Info("DrawSharingFilter", "Drawing sharing filter");
    TCollection* c = GetCollection(fSums, "fmdSharingFilter");
    if (!c) return;
    TCollection* rc = GetCollection(fResults, "fmdSharingFilter");
    if (!rc) rc = c;

    fBody->Divide(1, 3);
    fBody->cd(1);
  
    Double_t y = .8;
    Bool_t   angle=false, lowSignal=false, simple=false;

    if (GetParameter(c, "angle", angle))
      DrawParameter(y, "Angle correct", (angle ? "yes" : "no")); 
    if (GetParameter(c, "lowSignal", lowSignal))
      DrawParameter(y, "Lower signal",  (lowSignal ? "yes" : "no"));
    
    if (GetParameter(c, "simple", simple))
      DrawParameter(y, "Simple method", (simple ? "yes" : "no"));
    TParameter<int>* nFiles = 
      static_cast<TParameter<int>*>(GetObject(c, "nFiles"));
    if (nFiles)
      DrawParameter(y, "# files merged", Form("%d", nFiles->GetVal()));

    
    TCollection* lc = GetCollection(c, "lCuts");
    TCollection* hc = GetCollection(c, "hCuts");
    Int_t lm, hm;
    if (GetParameter(lc, "method", lm)) 
      DrawParameter(y, "Low cut method", (lm == 0 ? "fixed" : 
					  lm == 1 ? "fraction of MPV" : 
					  lm == 2 ? "fit range" : 
					  lm == 3 ? "Landau width" : 
					  "unknown"));
    if (GetParameter(hc, "method", hm)) 
      DrawParameter(y, "High cut method", (hm == 0 ? "fixed" : 
					  hm == 1 ? "fraction of MPV" : 
					  hm == 2 ? "fit range" : 
					  hm == 3 ? "Landau width" : 
					  "unknown"));

					  
    TH2* hLow  = GetH2(c, "lowCuts");
    TH2* hHigh = GetH2(c, "highCuts");
    if (hLow  && nFiles) hLow->Scale(1. / nFiles->GetVal());
    if (hHigh && nFiles) hHigh->Scale(1. / nFiles->GetVal());
    DrawInPad(fBody, 2, hLow,  "colz");
    DrawInPad(fBody, 3, hHigh, "colz");
  
    PrintCanvas("Sharing filter");

    const char* subs[] = { "FMD1I", "FMD2I", "FMD2O", "FMD3O", "FMD3I", 0 };
    const char** ptr   = subs;
    while (*ptr) { 
      TCollection* sc = GetCollection(c, *ptr);
      if (!sc) { ptr++; continue; }
    
      fBody->Divide(2,3);
      DrawInPad(fBody, 1, GetH1(sc, "esdEloss"),       "",     kLogy,
		"#Delta/#Delta_{mip} reconstructed and merged");
      DrawInPad(fBody, 1, GetH1(sc, "anaEloss"),       "same", kLogy|kLegend);
      DrawInPad(fBody, 2, GetH1(sc, "singleEloss"),    "",     kLogy,
		"#Delta/#Delta_{mip} for single, double, and tripple hits");
      DrawInPad(fBody, 2, GetH1(sc, "doubleEloss"),    "same", kLogy);
      DrawInPad(fBody, 2, GetH1(sc, "tripleEloss"),    "same", kLogy|kLegend);  
      DrawInPad(fBody, 3, GetH2(sc, "singlePerStrip"), "colz", kLogz);
      // DrawInPad(fBody, 4, GetH1(sc, "distanceBefore"), "",     0x2);
      // DrawInPad(fBody, 4, GetH1(sc, "distanceAfter"),  "same", 0x12);
      DrawInPad(fBody, 4, GetH2(sc, "summed"),         "colz", 0x0);

      TH2* nB = GetH2(sc, "neighborsBefore");
      if (nB) { 
	nB->GetXaxis()->SetRangeUser(0,8); 
	nB->GetYaxis()->SetRangeUser(0,8); 
      }
      DrawInPad(fBody, 5, nB, "colz", kLogz);
      DrawInPad(fBody, 5, GetH2(sc, "neighborsAfter"), "p same", kLogz,
		"Correlation of neighbors before and after merging");
      DrawInPad(fBody, 6, GetH2(sc, "beforeAfter"),    "colz",   kLogz);

      PrintCanvas(Form("Sharing filter - %s", *ptr));
      ptr++;
    }

    // --- MC --------------------------------------------------------
    TCollection* cc = GetCollection(c, "esd_mc_comparion", false); // Spelling!
    if (!cc) return; // Not MC 

    DivideForRings(false, false);
    DrawInRingPad(1, 'I', GetH2(cc, "FMD1i_corr"), "colz", kLogz);
    DrawInRingPad(2, 'I', GetH2(cc, "FMD2i_corr"), "colz", kLogz);
    DrawInRingPad(2, 'O', GetH2(cc, "FMD2o_corr"), "colz", kLogz);
    DrawInRingPad(3, 'I', GetH2(cc, "FMD3i_corr"), "colz", kLogz);
    DrawInRingPad(3, 'O', GetH2(cc, "FMD3o_corr"), "colz", kLogz);

    PrintCanvas("Sharing filter - MC vs Reco");

    // --- MC --------------------------------------------------------
    DrawTrackDensity(c);
  }

  //____________________________________________________________________
  void DrawDensityCalculator()
  {
    Info("DrawDensityCalculator", "Drawing density calculator");
    TCollection* c = GetCollection(fSums, "fmdDensityCalculator");
    if (!c) return;

    fBody->Divide(2, 2);
    fBody->cd(1);
  
    Double_t y = .8;
    Int_t maxParticles=0, phiAcceptance=0, etaLumping=0, phiLumping=0;
    Bool_t method=false, recalcEta=false, recalcPhi=false;
    Double_t size = fLandscape ? 0.06 : 0.04;
  
    GetParameter(c, "maxParticle", maxParticles);

    if (GetParameter(c, "phiAcceptance", phiAcceptance))
      DrawParameter(y, "#phi acceptance method", 
		    (phiAcceptance == 1 ? "N_{ch}" : 
		     phiAcceptance == 2 ? "#DeltaE" : "none"),       size);
    
    if (GetParameter(c, "etaLumping", etaLumping) &&
	GetParameter(c, "phiLumping", phiLumping))
      DrawParameter(y, "Region size (sector#timesstrip)", 
		    Form("%2d #times %2d", phiLumping, etaLumping),  size);
    
    if (GetParameter(c, "method", method))
      DrawParameter(y, "Method", (method ? "Poisson" : "#DeltaE"),   size); 

    if (GetParameter(c, "recalcEta", recalcEta))
      DrawParameter(y, "Recalculate #eta",(recalcEta ? "yes" : "no"),size); 

    if (GetParameter(c, "recalcPhi", recalcPhi))
      DrawParameter(y, "Recalculate #phi",(recalcPhi ? "yes" : "no"),size); 
    TParameter<int>* nFiles = 
      static_cast<TParameter<int>*>(GetObject(c, "nFiles"));
    if (nFiles)
      DrawParameter(y, "# files merged", Form("%d", nFiles->GetVal()));

    TCollection* lc = GetCollection(c, "lCuts");
    Int_t lm;
    if (GetParameter(lc, "method", lm)) 
      DrawParameter(y, "Low cut method", (lm == 0 ? "fixed" : 
					  lm == 1 ? "fraction of MPV" : 
					  lm == 2 ? "fit range" : 
					  lm == 3 ? "Landau width" : 
					  "unknown"));


    TVirtualPad* p = fBody; // fBody->cd(2);
    // p->Divide(3,1);

    TH1* accI = GetH1(c, "accI");
    TH1* accO = GetH1(c, "accO");
    if (accI) { 
      Double_t scale = 1./accI->GetMaximum();
      accI->Scale(scale); 
      accO->Scale(scale);
      accI->SetMinimum(0); 
    }
    TH2* lCuts = GetH2(c, "lowCuts");
    TH2* maxW  = GetH2(c, "maxWeights");
    if (nFiles && lCuts) lCuts->Scale(1. / nFiles->GetVal());
    if (nFiles && maxW)  maxW->Scale(1. / nFiles->GetVal());
    DrawInPad(p, 2, accI); 
    DrawInPad(p, 2, accO,  "same", kLegend); 
    DrawInPad(p, 3, lCuts, "colz");
    DrawInPad(p, 4, maxW,  "colz");
  
    PrintCanvas("Density calculator");

    const char* subs[] = { "FMD1I", "FMD2I", "FMD2O", "FMD3O", "FMD3I", 0 };
    const char** ptr   = subs;
    while (*ptr) { 
      TCollection* sc = GetCollection(c, *ptr);
      if (!sc) { ptr++; continue; }
    
      fBody->Divide(2,3);
    
      DrawInPad(fBody, 1, GetH2(sc, "elossVsPoisson"),   "colz",   kLogz);
      DrawInPad(fBody, 2, GetH1(sc, "diffElossPoisson"), "HIST E", kLogy);
      DrawInPad(fBody, 3, GetH1(sc, "occupancy"),        "",       kLogy);
      DrawInPad(fBody, 4, GetH1(sc, "eloss"),            "",       kLogy,
		"#Delta/#Delta_{mip} before and after cuts");
      DrawInPad(fBody, 4, GetH1(sc, "elossUsed"),        "same",kLogy|kLegend);
      TH1* phiB = GetH1(sc, "phiBefore");
      TH1* phiA = GetH1(sc, "phiAfter");
      if (phiB && phiA) { 
	phiA->Add(phiB, -1);
	phiA->Divide(phiB);
	phiA->SetTitle("#Delta#phi from Ip (x,y) correction");
	phiA->SetYTitle("(#phi_{after}-#phi_{before})/#phi_{before}");
      }
      DrawInPad(fBody, 5, phiA);
      DrawInPad(fBody, 6, GetH2(sc, "phiAcc"), "colz",   kLogz);
    
      PrintCanvas(Form("Density calculator - %s", *ptr));
      ptr++;    
    }

    TCollection* cc = GetCollection(c, "esd_mc_comparison", false); 
    if (!cc) return; // Not MC 

    fBody->Divide(2,5);
    DrawInPad(fBody, 1, GetH2(cc, "FMD1I_corr_mc_esd"), "colz", kLogz);
    DrawInPad(fBody, 3, GetH2(cc, "FMD2I_corr_mc_esd"), "colz", kLogz);
    DrawInPad(fBody, 5, GetH2(cc, "FMD2O_corr_mc_esd"), "colz", kLogz);
    DrawInPad(fBody, 7, GetH2(cc, "FMD3O_corr_mc_esd"), "colz", kLogz);
    DrawInPad(fBody, 9, GetH2(cc, "FMD3I_corr_mc_esd"), "colz", kLogz);
    DrawInPad(fBody, 2,  GetH1(cc, "FMD1I_diff_mc_esd"), "", kLogy);
    DrawInPad(fBody, 4,  GetH1(cc, "FMD2I_diff_mc_esd"), "", kLogy);
    DrawInPad(fBody, 6,  GetH1(cc, "FMD2O_diff_mc_esd"), "", kLogy);
    DrawInPad(fBody, 8,  GetH1(cc, "FMD3O_diff_mc_esd"), "", kLogy);
    DrawInPad(fBody, 10, GetH1(cc, "FMD3I_diff_mc_esd"), "", kLogy);

    PrintCanvas("Density calculator - MC vs Reco");
  }

  //____________________________________________________________________
  void DrawCorrector()
  {
    Info("DrawCorrector", "Drawing corrector"); 
    TCollection* c = GetCollection(fSums, "fmdCorrector"); 
    if (!c) return;
  
    fBody->cd();
  
    Double_t y = .8;  
    Bool_t secondary=false, vertexBias=false, acceptance=false, merging=false;  
    if (GetParameter(c, "secondary", secondary))
      DrawParameter(y, "Secondary corr.", secondary ? "yes" : "no");
    if (GetParameter(c, "acceptance", acceptance))
      DrawParameter(y, "Acceptance corr.", acceptance ? "yes" : "no");
    if (GetParameter(c, "vertexBias", vertexBias))
      DrawParameter(y, "Vertex bias corr.", vertexBias ? "yes" : "no");
    if (GetParameter(c, "merging", merging))  
      DrawParameter(y, "Merging eff.", merging ? "yes" : "no");
    
    PrintCanvas("Corrector");

    TCollection* cc = GetCollection(c, "esd_mc_comparison", false); 
    if (!cc) return; // Not MC 
    
    DivideForRings(false, false);

    DrawInRingPad(1, 'I', GetH2(cc, "FMD1I_esd_vs_mc"), "colz", 0x0);
    DrawInRingPad(2, 'I', GetH2(cc, "FMD2I_esd_vs_mc"), "colz", 0x0);
    DrawInRingPad(2, 'O', GetH2(cc, "FMD2O_esd_vs_mc"), "colz", 0x0);
    DrawInRingPad(3, 'O', GetH2(cc, "FMD3O_esd_vs_mc"), "colz", 0x0);
    DrawInRingPad(3, 'I', GetH2(cc, "FMD3I_esd_vs_mc"), "colz", 0x0);

    PrintCanvas("Corrector - MC vs Reco");
  }

  //____________________________________________________________________
  void DrawHistCollector()
  {
    Info("DrawHistCollector", "Drawing histogram collector");  
    TCollection* c = GetCollection(fSums, "fmdHistCollector");
    if (!c) return;

    fBody->Divide(1, 3);
    fBody->cd(1);

    Double_t y = .8;
    Int_t nCutBins=0, fiducial=0, merge=0, skipRings=0;
    Double_t fiducialCut=0.;
    Bool_t  bgAndHits=false;

    if (GetParameter(c, "nCutBins", nCutBins))
      DrawParameter(y, "# of bins to cut",      Form("%d", nCutBins));

    if (GetParameter(c, "skipRings", skipRings)) {
      TString skipped;
      if (skipRings & 0x05) skipped.Append("FMD1i ");
      if (skipRings & 0x09) skipped.Append("FMD2i ");
      if (skipRings & 0x0a) skipped.Append("FMD2o ");
      if (skipRings & 0x11) skipped.Append("FMD3i ");
      if (skipRings & 0x12) skipped.Append("FMD3o ");
      DrawParameter(y, "Skipped rings", skipped);
    }

    if (GetParameter(c, "bgAndHits", bgAndHits))
      DrawParameter(y, "Bg & hit maps stored.", bgAndHits?"yes":"no");

    if (GetParameter(c, "merge", merge))
      DrawParameter(y, "Merge method", 
		    (merge == 0 ? "straight mean" :
		     merge == 1 ? "straight mean, no zeroes" : 
		     merge == 2 ? "weighted mean" : 
		     merge == 3 ? "least error" : 
		     merge == 4 ? "sum" : "unknown"));

    if (GetParameter(c, "fiducial", fiducial))
      DrawParameter(y, "Fiducial method.", 
		    fiducial == 0 ? "cut" : "distance");

    if (GetParameter(c, "fiducialCut", fiducialCut))
      DrawParameter(y, "Fiducial cut.", Form("%f", fiducialCut));

		 
    DrawInPad(fBody, 2, GetH2(c, "sumRings"), "colz"); 
    DrawInPad(fBody, 3, GetH2(c, "coverage"), "colz");

    fBody->cd(1)->Modified();
    fBody->cd(2)->Modified();
    fBody->cd(3)->Modified();
    fBody->cd(1)->Update();
    fBody->cd(2)->Update();
    fBody->cd(3)->Update();
    PrintCanvas("Histogram collector");
		
    
    TIter next(c);
    TObject* o = 0;
    TRegexp regexp("[pm][0-9]+_[pm][0-9]+");
    while ((o = next())) { 
      TString name(o->GetName());
      if (name.Index(regexp) == kNPOS) continue;
      
      TList* vl = static_cast<TList*>(o);

      DivideForRings(false, false);
      
      DrawInRingPad(1, 'I', GetH2(vl, "secMapFMD1I"), "colz", 0x0);
      DrawInRingPad(2, 'I', GetH2(vl, "secMapFMD2I"), "colz", 0x0);
      DrawInRingPad(2, 'O', GetH2(vl, "secMapFMD2O"), "colz", 0x0);
      DrawInRingPad(3, 'O', GetH2(vl, "secMapFMD3O"), "colz", 0x0);
      DrawInRingPad(3, 'I', GetH2(vl, "secMapFMD3I"), "colz", 0x0);
      DrawInRingPad(1, 'I', GetH2(vl, "hitMapFMD1I"), "box same", 0x0);
      DrawInRingPad(2, 'I', GetH2(vl, "hitMapFMD2I"), "box same", 0x0);
      DrawInRingPad(2, 'O', GetH2(vl, "hitMapFMD2O"), "box same", 0x0);
      DrawInRingPad(3, 'O', GetH2(vl, "hitMapFMD3O"), "box same", 0x0);
      DrawInRingPad(3, 'I', GetH2(vl, "hitMapFMD3I"), "box same", 0x0);

      PrintCanvas(Form("Histogram Collector - Vertex bin %s", vl->GetName()));
    }

    o = c->FindObject("byCentrality");
    if (!o) return;
    TList* bc = static_cast<TList*>(o);

    DrawInPad(fBody, GetH3(bc, "FMD1I"), "box", 0);
    DrawInPad(fBody, GetH3(bc, "FMD2I"), "box same", 0);
    DrawInPad(fBody, GetH3(bc, "FMD2O"), "box same", 0);
    DrawInPad(fBody, GetH3(bc, "FMD3O"), "box same", 0);
    DrawInPad(fBody, GetH3(bc, "FMD3I"), "box same", kLegend);
  }

  //____________________________________________________________________
  void DrawCentral()
  {
    Info("DrawCentral", "Drawing central (SPD)");  
    TCollection* c = fSums; 
    if (!c) return;

    fBody->Divide(1, 3);
    fBody->cd(1);

		 
    DrawInPad(fBody, 1, GetH2(c, "coverage"), "col", 0,
	      "#eta coverage per v_{z}");
    DrawInPad(fBody, 2, GetH2(c, "nClusterVsnTracklet"), "colz", kLogx|kLogy,
	      "Correlation of # of tracklets and clusters"); 
    DrawInPad(fBody, 3, GetH2(c, "clusterPerTracklet"), "colz", 0x0,
	      "# clusters per tracklet vs #eta"); 

    fBody->cd(1)->Modified();
    fBody->cd(2)->Modified();
    fBody->cd(3)->Modified();
    fBody->cd(1)->Update();
    fBody->cd(2)->Update();
    fBody->cd(3)->Update();
    PrintCanvas("Central - overview");
		
    
    TIter next(c);
    TObject* o = 0;
    TRegexp regexp("[pm][0-9]+_[pm][0-9]+");
    while ((o = next())) { 
      TString name(o->GetName());
      if (name.Index(regexp) == kNPOS) continue;
      
      TList* vl = static_cast<TList*>(o);

      fBody->Divide(1, 3);
    
      DrawInPad(fBody, 1, GetH1(vl, "acceptance"), "", 0);

      TH1* sec = GetH1(vl, "secondary");
      sec->SetMarkerStyle(21);
      sec->SetMarkerSize(1.2);
      DrawInPad(fBody, 2, sec, "", 0);
      DrawInPad(fBody, 2, GetH1(vl, "secondaryFiducial"),    "same", 0x0);
      DrawInPad(fBody, 3, GetH2(vl, "secondaryMapFiducial"), "colz", 0);
      DrawInPad(fBody, 3, GetH2(vl, "hitMap"),               "box same", 0x0);

      fBody->cd(1)->Modified();
      fBody->cd(2)->Modified();
      fBody->cd(3)->Modified();
      fBody->cd(1)->Update();
      fBody->cd(2)->Update();
      fBody->cd(3)->Update();
      PrintCanvas(Form("Central - Vertex bin %s", vl->GetName()));
    }
  }

  
  //____________________________________________________________________
  void AddToAll(THStack* all, const THStack* stack, Int_t curr, Int_t step)
  {
    if (!stack) return;

    TIter   next(stack->GetHists());
    TH1*    h = 0;
    while ((h = static_cast<TH1*>(next()))) {
      TH1* copy = static_cast<TH1*>(h->Clone(Form("%s_copy", h->GetName())));
      copy->SetDirectory(0);
      if (curr != step) {
	copy->SetMarkerColor(kGray);
	copy->SetLineColor(kGray);
      }
      all->Add(copy);
    }
  }
  //____________________________________________________________________
  void AddToAll(THStack* all, const THStack* stack)
  {
    if (!stack) return;

    TIter   next(stack->GetHists());
    TH1*    h = 0;
    while ((h = static_cast<TH1*>(next()))) {
      TH1* copy = static_cast<TH1*>(h->Clone(Form("%s_copy", h->GetName())));
      copy->SetDirectory(0);
      copy->SetMarkerColor(kGray);
      copy->SetLineColor(kGray);
      all->Add(copy);
    }
  }

  //____________________________________________________________________
  void DrawStep(Int_t        step,
		THStack*     all,
		TObject*     cur,
		TLegend*     leg,
		const char*  title,
		TVirtualPad* can)
  {
    if (all->GetHists()->GetEntries() <= 0 || !cur) return;

    // Info("", "Drawing step # %d", step);
    Bool_t left = (step % 2) == 1; 
    TVirtualPad* p = can->cd(step);
    gStyle->SetOptTitle(0);
    p->SetTitle(Form("Step # %d", step));
    p->SetFillColor(kWhite);
    p->SetRightMargin(left ? 0 : 0.02);
    p->SetTopMargin(0); // 0.02);

    p->cd();
    all->Draw("nostack");
    all->GetHistogram()->SetXTitle("#eta");
    all->GetHistogram()->SetYTitle("signal");

    // p->cd();
    gROOT->SetSelectedPad(p);
    cur->DrawClone("same nostack");
    leg->DrawClone("");

    TLatex* ltx = new TLatex(.97, .97, title);
    ltx->SetNDC();
    ltx->SetTextSize(.06);
    ltx->SetTextAlign(33);
    ltx->Draw();

    ltx = new TLatex((left ? .12 : .02), .97, p->GetTitle());
    ltx->SetNDC();
    ltx->SetTextSize(.06);
    ltx->SetTextAlign(13);
    ltx->Draw();

    p->Modified();
    p->Update();
    p->cd();

    gStyle->SetOptTitle(1);
  }

  //____________________________________________________________________
  void FixStack(THStack* stack, const TString& title, 
		const TString& extra, Int_t marker)
  {
    if (!stack) return;
    stack->SetTitle(title);
    TIter next(stack->GetHists());
    TH1*  h = 0;
    while ((h = static_cast<TH1*>(next())))  {
      h->SetMarkerStyle(marker);
      TString tit(h->GetTitle());
      tit.ReplaceAll("cache", "");
      tit.Append(extra);
      h->SetTitle(tit);
    }
  }
  void AddLegendEntry(TLegend* l, 
		      const TH1* h, 
		      const TString& title)
  {
    if (!h) return;

    TLegendEntry* e = l->AddEntry("dummy", title.Data(), "pl");
    e->SetMarkerStyle(h->GetMarkerStyle());
    e->SetMarkerColor(kGray);
    e->SetLineColor(kGray);
    e->SetTextColor(kGray);
  }
		      

  //____________________________________________________________________
  void DrawSteps()
  {
    // MakeChapter(can, "Steps");

    THStack* deltas  = GetStack(GetCollection(fResults, "fmdSharingFilter"), 
				"sums", "summed");
    THStack* nchs    = GetStack(GetCollection(fResults, 
					      "fmdDensityCalculator"), 
				"sums", "inclDensity");
    THStack* prims   = GetStack(GetCollection(fResults, "fmdCorrector"), 
				"sums", "primaryDensity");
    THStack* rings   = GetStack(GetCollection(fResults, "ringResults"), "all");
    THStack* mcRings = GetStack(GetCollection(fResults, "mcRingResults", false),
				"all","dndeta_eta", false);
    TH1*     dndeta  = GetH1(fResults, "dNdeta");
    if (dndeta) dndeta->SetMarkerColor(kBlack);

    FixStack(deltas, "#sum_{} #Delta/#Delta_{mip}",  "",     20);
    FixStack(nchs,   "#sum_{} N_{ch,incl}", 	     "",     21);
    FixStack(prims,  "#sum_{} N_{ch,primary}",       "",     22);
    FixStack(rings,  "dN/d#eta per ring",            "",     23);
    FixStack(mcRings,"dN/d#eta per ring (MC)",       "(MC)", 34);

    THStack* all = new THStack;
    AddToAll(all, mcRings);
    AddToAll(all, deltas);
    AddToAll(all, nchs);
    AddToAll(all, prims);
    AddToAll(all, rings);

    TH1* res = 0;
    if (dndeta) {
      res = static_cast<TH1*>(dndeta->Clone("dNdeta"));
      res->SetTitle("dN/d#eta");
      res->SetMarkerColor(kGray);
      res->SetLineColor(kGray);
      res->SetDirectory(0);
      all->Add(res);
    }

    TLegend* l = new TLegend(.35, .2, .55, .9);
    l->SetFillColor(kWhite);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    TLegendEntry* e = 0;

    TH1* h = 0;
    if (mcRings) {
      h = static_cast<TH1*>(mcRings->GetHists()->At(0));
      AddLegendEntry(l, h, mcRings->GetTitle());
    }

    if (deltas) {
      h = static_cast<TH1*>(deltas->GetHists()->At(0));
      AddLegendEntry(l, h, deltas->GetTitle());    
    }

    if (nchs) {
      h = static_cast<TH1*>(nchs->GetHists()->At(0));
      AddLegendEntry(l, h, nchs->GetTitle());    
    }

    if (prims) {
      h = static_cast<TH1*>(prims->GetHists()->At(0));
      AddLegendEntry(l, h, prims->GetTitle());    
    }

    if (rings) {
      h = static_cast<TH1*>(rings->GetHists()->At(0));
      AddLegendEntry(l, h, rings->GetTitle());    
    }

    if (res) {
      h = res;
      AddLegendEntry(l, h, h->GetTitle());
    }
    
    TObject* objs[] = { mcRings, 
			deltas, 
			nchs, 
			prims, 
			rings, 
			dndeta };
    const char* titles[] = { "MC", 
			     "After merging", 
			     "After particle counting", 
			     "After corrections", 
			     "After normalization", 
			     "After combining" };
    
    fBody->Divide(2, 3, 0, 0); // mcRings ? 4 : 3, 0, 0);
    Int_t step = 0;
    for (Int_t i = 0; i < 6; i++) { 
      TObject* obj = objs[i];
      if (!obj) continue;
      e = static_cast<TLegendEntry*>(l->GetListOfPrimitives()->At(step));
      e->SetMarkerColor(kBlack);
      e->SetLineColor(kBlack);
      e->SetTextColor(kBlack);
      step++;
      
      DrawStep(step, all, obj, l, titles[i], fBody);
      e->SetMarkerColor(kGray);
      e->SetLineColor(kGray);
      e->SetTextColor(kGray);
    }

    if (!mcRings && deltas) { 
      fBody->cd(6);
      TLegend* ll = new TLegend(0.01, 0.11, 0.99, 0.99);
      // ll->SetNDC();
      ll->SetFillColor(kWhite);
      ll->SetFillStyle(0);
      ll->SetBorderSize(0);

      TIter next(deltas->GetHists());
      TH1*  hh = 0;
      while ((hh = static_cast<TH1*>(next()))) {
	e = ll->AddEntry("dummy", hh->GetTitle(), "pl");
	e->SetMarkerColor(hh->GetMarkerColor());
	e->SetMarkerStyle(hh->GetMarkerStyle());
	e->SetLineColor(kBlack);
      }
      ll->Draw();
    }
    PrintCanvas("Steps");
  }


  //____________________________________________________________________
  void DrawResults()
  {
    // MakeChapter(can, "Results");

    fBody->Divide(2,2);

    TCollection* c = GetCollection(fResults, "ringResults");
    if (!c) return;
  
    THStack* mcRings = GetStack(GetCollection(fResults, "mcRingResults", false),
				"all", "dndeta_eta", false);

    TH1* dndeta_phi = GetH1(fResults, "dNdeta");
    TH1* dndeta_eta = GetH1(fResults, "dNdeta_");
    dndeta_phi->SetTitle("1/N_{ev}dN_{ch}/d#eta (#varphi norm)");
    dndeta_eta->SetTitle("1/N_{ev}dN_{ch}/d#eta (#eta norm)");
    dndeta_eta->SetMarkerSize(0.7);

    THStack* allPhi = new THStack("phiAcc", "#varphi Acceptance");
    THStack* allEta = new THStack("etaCov", "#eta Coverage");
    const char*  rings[] = { "FMD1I", "FMD2I", "FMD2O", "FMD3I", "FMD3O", 0 };
    const char** pring   = rings;
    
    while ((*pring)) { 
      TCollection* cc     = GetCollection(c, *pring);
      TH1*         etaCov = GetH1(cc, "etaCov");
      TH1*         phiAcc = GetH1(cc, "phiAcc");
      TH1*         dndeta = GetH1(cc, "dndeta_phi");
      Int_t        color  = kBlack;
      if (dndeta)  color  = dndeta->GetMarkerColor();
      if (etaCov) { 
	etaCov->SetTitle(*pring);
	etaCov->SetFillColor(color);
	etaCov->SetLineColor(color);
	allEta->Add(etaCov);
      }
      if (phiAcc) { 
	phiAcc->SetFillColor(color);
	phiAcc->SetLineColor(color);
	allPhi->Add(phiAcc);
      }
      pring++;
    }
    DrawInPad(fBody, 1, GetStack(c, "all"), "nostack", mcRings ? 0 : kLegend,
	      "Individual ring results");
    DrawInPad(fBody, 1, mcRings, "nostack same", kLegend);
    DrawInPad(fBody, 2, dndeta_phi, 
	      "1/#it{N}_{ev} d#it{N}_{ch}/d#it{#eta}");
    DrawInPad(fBody, 2, dndeta_eta, "Same", kLegend);
    DrawInPad(fBody, 3, allEta, "nostack hist", kLegend,
	      "#phi acceptance and #eta coverage per ring");
    DrawInPad(fBody, 3, allPhi, "nostack hist same", 0x0);
    DrawInPad(fBody, 4, GetH1(fResults, "norm"), "", 0x0, 
	      "Total #phi acceptance and #eta coverage");
    DrawInPad(fBody, 4, GetH1(fResults, "phi"), "same", kLegend);
    // DrawInPad(fBody, 4, GetH1(fSums,    "d2Ndetadphi"), "colz");

    // fBody->cd(1);
    // TLatex* l = new TLatex(.5, .2, "Ring results");
    // l->SetNDC();
    // l->SetTextAlign(21);
    // l->Draw();

    // fBody->cd(2);
    // l->DrawLatex(.5, .2, "1/N_{ev}dN_{ch}/d#eta");

    // fBody->cd(3);
    // l->DrawLatex(.5, .2, "1/N_{ev}dN_{ch}/d#eta (#vta norm.)");
    
    PrintCanvas("Results");
  }

  //____________________________________________________________________
  void DrawCentralResults()
  {
    // MakeChapter(can, "Results");
    Info("DrawCentralResults", "Drawing central results");

    fBody->Divide(1,2,0,0);

    TH1* dndeta_ = GetH1(fResults, "dNdeta_");
    TH1* dndeta  = GetH1(fResults, "dNdeta");
    THStack* stack = new THStack("dndetas", 
				 "d#it{N}_{ch}/d#it{#eta} - central");
    stack->Add(dndeta_);
    stack->Add(dndeta);
    
    DrawInPad(fBody, 1, stack, "nostack");
    stack->GetHistogram()->SetXTitle("#it{#eta}");
    stack->GetHistogram()->SetYTitle("#frac{d#it{N}_{ch}}{d#it{#eta}}");
    fBody->cd(1);
    TLegend* l = new TLegend(.3, .05, .7, .4);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->AddEntry(dndeta_, "Normalized to coverage",        "lp");
    l->AddEntry(dndeta,  "Normalized to #phi acceptance", "lp");
    l->Draw();

    DrawInPad(fBody, 2, GetH1(fResults, "norm"));
    DrawInPad(fBody, 2, GetH1(fResults, "phi"), "same", kLegend);

    PrintCanvas("Central Results");

  }
  void DrawBoth(TFile* file)
  {
    Info("DrawBoth", "Drawing central & forward results");
    TCollection* central = GetCollection(file, "CentralResults");
    TCollection* forward = GetCollection(file, "ForwardResults");
    
    if (!central || !forward) {
      Warning("DrawBoth", "central %p or forward %p results not found", 
	      central, forward);
      return;
    }

    TH1* f1 = GetH1(forward, "dNdeta_");
    TH1* c1 = GetH1(central, "dNdeta_");
    TH1* f2 = GetH1(forward, "dNdeta");
    TH1* c2 = GetH1(central, "dNdeta");
    f1->SetLineColor(kBlack);
    f2->SetLineColor(kBlack);
    c1->SetLineColor(kBlack);
    c2->SetLineColor(kBlack);
    f1->SetMarkerColor(f2->GetMarkerColor());
    f1->SetMarkerStyle(24);
    c1->SetMarkerStyle(24);
    c2->SetMarkerStyle(20);
    c2->SetMarkerColor(c1->GetMarkerColor());
    THStack* s = new THStack("dndetas", "d#it{N}_{ch}/d#it{#eta}");
    s->Add(f1);
    s->Add(c1);
    s->Add(f2);
    s->Add(c2);

    fBody->Divide(1, 2, 0, 0);
    DrawInPad(fBody, 1, s, "nostack");
    s->GetHistogram()->SetXTitle("#it{#eta}");
    s->GetHistogram()->SetYTitle("#frac{d#it{N}_{ch}}{d#it{#eta}}");
    
    fBody->cd(1);
    TLegend* l = new TLegend(.4, .05, .8, .4);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    TLegendEntry* entry = l->AddEntry("dummy", "Forward", "f");
    entry->SetFillColor(f1->GetMarkerColor());
    entry->SetLineColor(f1->GetMarkerColor());
    entry->SetFillStyle(1001);
    entry->SetLineWidth(0);
    entry = l->AddEntry("dummy", "Central", "f");
    entry->SetFillColor(c1->GetMarkerColor());
    entry->SetLineColor(c1->GetMarkerColor());
    entry->SetLineWidth(0);
    entry->SetFillStyle(1001);
    entry = l->AddEntry("dummy", "Normalized to coverage", "lp");
    entry->SetMarkerStyle(f1->GetMarkerStyle());
    entry = l->AddEntry("dummy", "Normalized to #phi acceptance", "lp");
    entry->SetMarkerStyle(f2->GetMarkerStyle());
    l->Draw();

    TH1* f3 = GetH1(forward, "norm");
    TH1* c3 = GetH1(central, "norm");
    TH1* f4 = GetH1(forward, "phi");
    TH1* c4 = GetH1(central, "phi");
    f3->SetFillColor(f1->GetMarkerColor());
    f4->SetFillColor(f1->GetMarkerColor());
    c3->SetFillColor(c1->GetMarkerColor());
    c4->SetFillColor(c1->GetMarkerColor());
    f3->SetLineColor(f1->GetMarkerColor());
    f4->SetLineColor(f1->GetMarkerColor());
    c3->SetLineColor(c1->GetMarkerColor());
    c4->SetLineColor(c1->GetMarkerColor());
    
    THStack* a = new THStack("norms", "Normalizations");
    a->Add(f3);
    a->Add(c3);
    a->Add(f4);
    a->Add(c4);
    
    a->SetMaximum(a->GetMaximum("nostack")*1.2);
    DrawInPad(fBody, 2, a, "nostack");
    a->GetHistogram()->SetXTitle("#it{#eta}");
    a->GetHistogram()->SetYTitle("Normalization (coverage or acceptance)");
    
    fBody->cd(2);
    l = new TLegend(.2, .94, .9, .99);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetNColumns(2);
    // entry = l->AddEntry("dummy", "Forward", "f");
    // entry->SetFillColor(f1->GetMarkerColor());
    // entry->SetLineColor(f1->GetMarkerColor());
    // entry->SetFillStyle(1001);
    // entry->SetLineWidth(0);
    // entry = l->AddEntry("dummy", "Central", "f");
    // entry->SetFillColor(c1->GetMarkerColor());
    // entry->SetLineColor(c1->GetMarkerColor());
    // entry->SetLineWidth(0);
    // entry->SetFillStyle(1001);
    entry = l->AddEntry("dummy", "#eta Coverage", "f");
    entry->SetFillStyle(f3->GetFillStyle());
    entry->SetFillColor(kBlack);
    entry = l->AddEntry("dummy", "#phi Acceptance", "f");
    entry->SetFillStyle(f4->GetFillStyle());
    entry->SetFillColor(kBlack);
    l->Draw();

    PrintCanvas("Both results");
  }
  TCollection* fSums;
  TCollection* fResults;
};

// #endif
