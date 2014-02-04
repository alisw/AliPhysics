#include "SummaryDrawer.C"
#ifndef __CINT__
# include <TGraph.h>
# include <TGraphErrors.h>
# include <TF1.h>
# include <TArrow.h>
#else
class TGraph;
#endif

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
    if (!centralSums) {
      Info("Run", "Trying old name \"Central\"");
      centralSums = GetCollection(file, "Central", false);
    }
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
      if (!fSums) 
	fSums = GetCollection(file, "Central");
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
  TGraph* CreateCutGraph(Int_t method, Int_t iy, TH2* cuts, TH1* eloss,
			 Int_t color)
  {
    TGraph*  ret = new TGraph(4);
    Double_t y0  = TMath::Max(eloss->GetMinimum(),1.);
    Double_t y1  = eloss->GetMaximum();
    Double_t min = 1000;
    Double_t max = 0;
    if (method == 0) { // Fixed value
      max = cuts->GetBinContent(1, iy);
      min = eloss->GetXaxis()->GetXmin();
    }
    else {
      for (Int_t ix=1; ix <= cuts->GetNbinsX(); ix++) {
	Double_t c = cuts->GetBinContent(ix, iy);
	if (c <= 0.0001) continue;
	min = TMath::Min(c, min);
	max = TMath::Max(c, max);
      }
    }
    // Printf("Cuts between %f,%f @%f, %f", min, max, y0,y1);
    ret->SetPoint(0, min, y0); 
    ret->SetPoint(1, min, y1); 
    ret->SetPoint(2, max, y1);
    ret->SetPoint(3, max, y0);
    ret->SetFillColor(color);
    ret->SetFillStyle(3002);
    ret->SetLineColor(kBlack);
    ret->SetLineStyle(2);
    ret->SetName(Form("g%s", cuts->GetName()));
    ret->SetTitle(cuts->GetTitle());
    
    return ret;
  }
  //____________________________________________________________________
  const Char_t* CutMethodName(Int_t lm) const
  {
    switch (lm) {
    case 0: return "fixed";
    case 1: return "fraction of MPV";
    case 2: return "fit range";
    case 3: return "Landau width";
    case 4: return "Probability";
    }
    return "unknown";
  }
  //____________________________________________________________________
  Int_t PrintCut(const TCollection* c, Double_t& y, const Char_t* name,
		 Double_t size=0)
  {
    if (!c) return -1;

    Int_t method = 0;
    if (!GetParameter(c, "method", method)) return -1;

    Double_t val = 0;
    Bool_t   sig = false;
    switch (method) {
    case 0: // Fixed 
      DrawParameter(y, name, "fixed", size); 
      break;
    case 1: // MPV
      GetParameter(c, "frac", val);
      DrawParameter(y, name, Form("Fraction of #Delta_{p} (%3d%%)",
				  Int_t(val*100)), size);
      break;
    case 2: // Fit range
      DrawParameter(y, name, "Fit range"); 
      break;
    case 3: // Landau width
      GetParameter(c, "nXi", val);
      GetParameter(c, "sigma", sig);
      DrawParameter(y, name, Form("N#times%s#xi%s (N=%4.1f)",
				  (sig ? "(" : ""),
				  (sig ? "+#sigma)" : ""), 
				  val), size);
      break;
    case 4: // Probability;
      GetParameter(c, "probability", val);
      DrawParameter(y, name, Form("P(#Delta)<%f", val), size);
      break;
    default:
      DrawParameter(y, name, "Unknown", size);
      break;
    }
    return method;
  }
  //____________________________________________________________________
  void DrawSharingFilter()
  {
    Info("DrawSharingFilter", "Drawing sharing filter");
    TCollection* c = GetCollection(fSums, "fmdSharingFilter");
    if (!c) return;
    TCollection* rc = GetCollection(fResults, "fmdSharingFilter");
    if (!rc) rc = c;

    // --- Draw summary information ----------------------------------
    fBody->Divide(1, 3, 0, 0);
    fBody->cd(1);
  
    Double_t y = .8;
    Bool_t   angle=false, lowSignal=false, disabled=false;

    if (GetParameter(c, "angle", angle))
      DrawParameter(y, "Angle correct", (angle ? "yes" : "no")); 
    if (GetParameter(c, "lowSignal", lowSignal))
      DrawParameter(y, "Lower signal",  (lowSignal ? "yes" : "no"));
    TParameter<int>* nFiles = 
      static_cast<TParameter<int>*>(GetObject(c, "nFiles"));
    if (nFiles)
      DrawParameter(y, "# files merged", Form("%d", nFiles->GetVal()));    
    if (GetParameter(c, "disabled", disabled)) 
      DrawParameter(y, "Merging disabled", (disabled ? "yes" : "no"));

    Int_t lm    = 0;
    Int_t hm    = 0;
    TH2*  hLow  = 0;
    TH2*  hHigh = 0;
    if (!disabled) {
      Bool_t simple=false, three=false;
      if (GetParameter(c, "simple", simple))
	DrawParameter(y, "Simple method", (simple ? "yes" : "no"));
      if (GetParameter(c, "sumThree", three)) 
	DrawParameter(y, "3-strip merging", (three ? "yes" : "no"));
    
      TCollection* lc = GetCollection(c, "lCuts");
      TCollection* hc = GetCollection(c, "hCuts");
      lm              = PrintCut(lc, y, "Low cut");
      hm              = PrintCut(hc, y, "High cut");
      hLow            = GetH2(c, "lowCuts");
      hHigh           = GetH2(c, "highCuts");
      // if (hLow  && nFiles) hLow->Scale(1. / nFiles->GetVal());
      // if (hHigh && nFiles) hHigh->Scale(1. / nFiles->GetVal());
      DrawInPad(fBody, 2, hLow,  "colz");
      DrawInPad(fBody, 3, hHigh, "colz");
    }
    PrintCanvas("Sharing filter");

    if (!disabled) {
      // --- Draw rings individually -----------------------------------
      Double_t savX = fParVal->GetX();
      Double_t savY = fParVal->GetY();
      fParVal->SetX(0.6);
      fParVal->SetY(0.6);
      const char** ptr   = GetRingNames(false);
      UShort_t     iq    = 1;
      while (*ptr) { 
	TCollection* sc = GetCollection(c, *ptr);
	if (!sc) { ptr++; iq++; continue; }
	
	if (fLandscape) fBody->Divide(3, 2);
	else            fBody->Divide(2,3);
	
	TH1*    esdELoss = GetH1(sc, "esdEloss");
	esdELoss->GetXaxis()->SetRangeUser(-.1, 2);
	TGraph* lowCut   = CreateCutGraph(lm, iq,  hLow,  esdELoss, kYellow+1);
	TGraph* highCut  = CreateCutGraph(hm, iq,  hHigh, esdELoss, kCyan+1);
	
	DrawInPad(fBody, 1, esdELoss, "", kLogy,
		  "#Delta/#Delta_{mip} reconstructed and merged");
	DrawInPad(fBody, 1, GetH1(sc, "anaEloss"), "same");
	DrawInPad(fBody, 1, lowCut,  "lf same"); 
	DrawInPad(fBody, 1, highCut, "lf same", kLogy|kLegend); 
	
	DrawInPad(fBody, 2, GetH1(sc, "singleEloss"),    "",    kLogy,
		  "#Delta/#Delta_{mip} for single, double, and tripple hits");
	DrawInPad(fBody, 2, GetH1(sc, "doubleEloss"),    "same",kLogy);
	DrawInPad(fBody, 2, GetH1(sc, "tripleEloss"),    "same",kLogy|kLegend);
	
	DrawInPad(fBody, 3, GetH2(sc, "singlePerStrip"), "colz",kLogz);
	// DrawInPad(fBody, 4, GetH1(sc, "distanceBefore"), "",     0x2);
	// DrawInPad(fBody, 4, GetH1(sc, "distanceAfter"),  "same", 0x12);
	DrawInPad(fBody, 4, GetH2(sc, "summed"),         "colz",0x0);
	
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
	iq++;
      }
      fParVal->SetX(savX);
      fParVal->SetY(savY);
    }

    // --- MC --------------------------------------------------------
    TCollection* cc = GetCollection(c, "esd_mc_comparion", false); // Spelling!
    if (!cc) return; // Not MC 

    DivideForRings(false, false);
    const char** ptr = GetRingNames(false);
    while (*ptr) { 
      DrawInRingPad(GetH2(cc, Form("%s_corr", *ptr)), "colz", kLogz);
      ptr++;
    }

    PrintCanvas("Sharing filter - MC vs Reco");

    // --- MC --------------------------------------------------------
    DrawTrackDensity(c);
  }

  void ShowSliceFit(Bool_t inY, TH2* h, Double_t nVar, 
		    TVirtualPad* p, Int_t sub, UShort_t flags=0,
		    Double_t cut=-1)
  {
    if (!h) return;

    TObjArray* fits = new TObjArray;
    fits->SetOwner();
    if (inY) h->FitSlicesY(0, 1, -1, 10, "QN", fits);
    else     h->FitSlicesX(0, 1, -1, 10, "QN", fits);
    if (!fits) { 
      Warning("ShowSliceFit", "No fits returned");
      return;
    }
    TH1* mean = static_cast<TH1*>(fits->At(1));
    TH1* var  = static_cast<TH1*>(fits->At(2));
    if (!mean || !var) {
      Warning("ShowSliceFit", "Didn't get histograms");
      fits->Delete();
      return;
    }
    TF1* fmean = new TF1("mean", "pol1");
    TF1* fvar  = new TF1("var",  "pol1");
    mean->Fit(fmean, "Q0+");
    var->Fit(fvar, "Q0+");
    if (!fmean || !fvar) {
      Warning("ShowSliceFit", "No functions returned");
      fits->Delete();
      return;
    }

    TGraphErrors* g = new TGraphErrors(h->GetNbinsX());
    g->SetName(Form("g%s", h->GetName()));
    TString xTit = h->GetXaxis()->GetTitle();
    TString yTit = h->GetYaxis()->GetTitle();
    g->SetTitle(Form("Correlation of %s and %s",
		     inY  ? xTit.Data() : yTit.Data(), 
		     !inY ? xTit.Data() : yTit.Data()));
    g->SetFillColor(kBlue-10);
    g->SetFillStyle(3001);
    TGraph* up  = (cut > 0 ? new TGraph(h->GetNbinsX()) : 0);
    TGraph* low = (cut > 0 ? new TGraph(h->GetNbinsX()) : 0);
    if (up)  { 
      up ->SetLineColor(kBlack); 
      up ->SetLineWidth(2); 
      up ->SetLineStyle(2); 
    }
    if (low) { 
      low->SetLineColor(kBlack); 
      low->SetLineWidth(2); 
      low->SetLineStyle(2); 
    }
    for (Int_t i = 1; i <= h->GetNbinsX(); i++) { 
      Double_t x  = h->GetXaxis()->GetBinCenter(i);
      Double_t y  = fmean->Eval(x);
      Double_t e  = fvar->Eval(x);
      Double_t ee = nVar * e;
      if (flags & 0x8000) ee *= e > 0 ? TMath::Log10(e) : 1;
      g->SetPoint(i-1, x, y);
      g->SetPointError(i-1, 0, ee);

      if (up)  up ->SetPoint(i-1,x,x+cut*x);
      if (low) low->SetPoint(i-1,x,x-cut*x);
    }
    DrawInPad(p, sub, g, "3", flags);
    if (up)  DrawInPad(p, sub, up, "l", flags);
    if (low) DrawInPad(p, sub, low, "l", flags);
    fmean->SetRange(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    fmean->SetLineWidth(2);
    DrawInPad(p, sub, fmean, "same", flags);

    TVirtualPad* pp = p->GetPad(sub);
    Double_t y = 1-pp->GetTopMargin()-.01;
    TLatex* l = new TLatex(.15, y, 
			   Form("#LT%s#GT(%s) = "
				"%f + %f %s", 
				yTit.Data(), xTit.Data(), 
				fmean->GetParameter(0), 
				fmean->GetParameter(1), xTit.Data()));
    l->SetNDC();
    l->SetTextAlign(13);
    l->SetTextSize(0.04);
    l->SetTextFont(42);
    l->Draw();
    l->DrawLatex(0.15, y-0.07, 
		 Form("#sigma_{%s}(%s) = "
		      "%f + %f %s", 
		      yTit.Data(), xTit.Data(),
		      fvar->GetParameter(0), 
		      fvar->GetParameter(1),  xTit.Data()));
    l->DrawLatex(0.15, y-0.14, Form("#delta = %3.1f %s #sigma",
				    nVar, 
				    flags & 0x8000 ? "log_{10}(#sigma)" : ""));
    fits->Delete();
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
    Double_t maxOutliers=0, outlierCut=0;
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
    if (GetParameter(c, "maxOutliers", maxOutliers))
      DrawParameter(y, "Max relative N_{outlier}",
		    Form("%5.3f",maxOutliers),size);
    if (GetParameter(c, "outlierCut", outlierCut))
      DrawParameter(y, "Max relative deviation",Form("%5.3f",outlierCut),size);


    TParameter<int>* nFiles = 
      static_cast<TParameter<int>*>(GetObject(c, "nFiles"));
    if (nFiles)
      DrawParameter(y, "# files merged", Form("%d", nFiles->GetVal()), size);

    TCollection* lc = GetCollection(c, "lCuts");
    PrintCut(lc, y, "Threshold", size);

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
    if (lCuts)           lCuts->SetTitle("Thresholds");
    if (nFiles && lCuts) lCuts->Scale(1. / nFiles->GetVal());
    if (nFiles && maxW)  maxW->Scale(1. / nFiles->GetVal());
    DrawInPad(p, 2, accI); 
    DrawInPad(p, 2, accO,  "same", kLegend); 
    DrawInPad(p, 3, lCuts, "colz");
    DrawInPad(p, 4, maxW,  "colz");
  
    PrintCanvas("Density calculator");

    const char** ptr   = GetRingNames(false);
    while (*ptr) { 
      TCollection* sc = GetCollection(c, *ptr);
      if (!sc) { ptr++; continue; }
    
      if (fLandscape) fBody->Divide(3,2);
      else            fBody->Divide(2,3);
    
      TH2* corr      = GetH2(sc, "elossVsPoisson");
      TH2* corrOut   = GetH2(sc, "elossVsPoissonOutlier");
      TH1* diff      = GetH1(sc, "diffElossPoisson");
      TH1* diffOut   = GetH1(sc, "diffElossPoissonOutlier");
      TH1* eloss     = GetH1(sc, "eloss");
      TH1* elossUsed = GetH1(sc, "elossUsed");
      TH1* occ       = GetH1(sc, "occupancy");
      if (eloss)     eloss    ->SetLineWidth(1);
      if (elossUsed) elossUsed->SetLineWidth(1);
      
      DrawInPad(fBody, 1, corr,    "colz",        kLogz);
      DrawInPad(fBody, 1, corrOut, "same",        kLogz);
      DrawInPad(fBody, 2, diff,    "HIST E",      kLogy);
      DrawInPad(fBody, 2, diffOut, "HIST E SAME", kLogy|kLegend);
      DrawInPad(fBody, 3, occ,      "",           kLogy);
      DrawInPad(fBody, 4, eloss,    "",           kLogy, 
		"#Delta/#Delta_{mip} before and after cuts");
      DrawInPad(fBody, 4, elossUsed, "same",      kLogy|kLegend);
      TH1* phiB = GetH1(sc, "phiBefore");
      TH1* phiA = GetH1(sc, "phiAfter");
      TH1* outliers = GetH1(sc, "outliers");
      if (outliers) 
	DrawInPad(fBody, 5, outliers, "hist", kLogy);
      else if (phiB && phiA) { 
	phiA->Add(phiB, -1);
	phiA->Divide(phiB);
	phiA->SetTitle("#Delta#phi from Ip (x,y) correction");
	phiA->SetYTitle("(#phi_{after}-#phi_{before})/#phi_{before}");
	DrawInPad(fBody, 5, phiA);
      }
      else {
	fBody->cd(5);
	TLatex* ltx = new TLatex(0.5, 0.5, "No outliers or #phi corrections");
	ltx->SetTextAlign(22);
	ltx->SetTextSize(0.07);
	ltx->SetNDC();
	ltx->Draw();
      }
      DrawInPad(fBody, 6, GetH2(sc, "phiAcc"), "colz",   kLogz);
    
      ShowSliceFit(true, corr, 10, fBody, 1, kLogz, outlierCut);

      if (diff && diffOut) { 
	fBody->cd(2);
	Double_t in  = diff->GetEntries();
	Double_t out = diffOut->GetEntries();
	if ((in+out) > 0) {
	  TLatex*  ltx = new TLatex(0.11, 0.89, 
				    Form("Fraction: %7.3f%%", 
					 100*out/(in+out)));
	  ltx->SetNDC();
	  ltx->SetTextAlign(13);
	  ltx->SetTextSize(0.06);
	  ltx->Draw();
	}
      }
      PrintCanvas(Form("Density calculator - %s", *ptr));
      ptr++;    
    }

    TCollection* cc = GetCollection(c, "esd_mc_comparison", false); 
    if (!cc) return; // Not MC 

    fBody->Divide(2,5);
    ptr   = GetRingNames(false);
    Int_t cnt = 0;
    while (*ptr) { 
      DrawInPad(fBody, 2*cnt+1, GetH2(cc, Form("%s_corr_mc_esd", *ptr)),
		"colz", kLogz);
      DrawInPad(fBody, 2*(cnt+1), GetH2(cc, Form("%s_diff_mc_esd", *ptr)),
		"", kLogz);
      ptr++;
      cnt++;
    }

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
    const char** ptr = GetRingNames(false);
    while (*ptr) { 
      DrawInRingPad(GetH2(cc, Form("%s_esd_vs_mc", *ptr)), "colz", 0x0);
      ptr++;
    }

    PrintCanvas("Corrector - MC vs Reco");
  }

  //____________________________________________________________________
  void DrawHistCollector()
  {
    Info("DrawHistCollector", "Drawing histogram collector");  
    TCollection* c = GetCollection(fSums, "fmdHistCollector");
    if (!c) return;

    fBody->Divide(2, 1);
    TVirtualPad* p = fBody->cd(1);
    p->Divide(1,2);
    p->cd(1);

    Double_t y = .8;
    Int_t nCutBins=0, fiducial=0, merge=0, skipRings=0;
    Double_t fiducialCut=0.;
    Bool_t  bgAndHits=false;
    Double_t size = fLandscape ? 0.06 : 0.04;
    if (GetParameter(c, "nCutBins", nCutBins))
      DrawParameter(y, "# of bins to cut", Form("%d", nCutBins),size);

    if (GetParameter(c, "skipRings", skipRings)) {
      TString skipped;
      if (skipRings & 0x05) skipped.Append("FMD1i ");
      if (skipRings & 0x09) skipped.Append("FMD2i ");
      if (skipRings & 0x0a) skipped.Append("FMD2o ");
      if (skipRings & 0x11) skipped.Append("FMD3i ");
      if (skipRings & 0x12) skipped.Append("FMD3o ");
      if (skipped.IsNull()) skipped = "none";
      DrawParameter(y, "Skipped rings", skipped, size);
    }
    if (GetParameter(c, "bgAndHits", bgAndHits))
      DrawParameter(y, "Bg & hit maps stored.", bgAndHits?"yes":"no",size);
    if (GetParameter(c, "merge", merge))
      DrawParameter(y, "Merge method", 
		    (merge == 0 ? "straight mean" :
		     merge == 1 ? "straight mean, no zeroes" : 
		     merge == 2 ? "weighted mean" : 
		     merge == 3 ? "least error" : 
		     merge == 4 ? "sum" : "unknown"),size);
    if (GetParameter(c, "fiducial", fiducial))
      DrawParameter(y, "Fiducial method.", 
		    fiducial == 0 ? "cut" : "distance", size);
    if (GetParameter(c, "fiducialCut", fiducialCut))
      DrawParameter(y, "Fiducial cut.", Form("%f", fiducialCut), size);

    // p->cd(2);
    // Printf("Drawing skipped");
    TH1* skipped = GetH1(c, "skipped");
    if (skipped) { 
      skipped->SetFillColor(kRed+1);
      skipped->SetFillStyle(3001);
    }
    DrawInPad(p, 2, skipped, "hist");
		 
    p = fBody->cd(2);
    p->Divide(1,2,0,0);

    // Printf("Drawing sumRings");
    DrawInPad(p, 1, GetH2(c, "sumRings"), "colz"); 
    // Printf("Drawing coverage");
    DrawInPad(p, 2, GetH2(c, "coverage"), "colz");
    // Printf("Done drawing for now");
    PrintCanvas("Histogram collector");
		
    
    TIter next(c);
    TObject* o = 0;
    TRegexp regexp("[pm][0-9]+_[pm][0-9]+");
    while ((o = next())) { 
      TString name(o->GetName());
      if (name.Index(regexp) == kNPOS) continue;
      
      TList* vl = static_cast<TList*>(o);
      if (!vl) continue;

      DivideForRings(false, false);
      const char** ptr = GetRingNames(false);
      while (*ptr) { 
	DrawInRingPad(GetH2(vl, Form("secMap%s", *ptr)), "colz", 0x0);
	DrawInRingPad(GetH2(vl, Form("hitMap%s", *ptr)), "box same", 0x0);
	ptr++;
      }
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

    fBody->Divide(2, 2);
    fBody->cd(1);
    Double_t y = .7;  
    Bool_t secondary=false, acceptance=false;  
    if (GetParameter(c, "secondary", secondary))
      DrawParameter(y, "Secondary corr.", secondary ? "yes" : "no");
    if (GetParameter(c, "acceptance", acceptance))
      DrawParameter(y, "Acceptance corr.", acceptance ? "yes" : "no");

		 
    DrawInPad(fBody, 2, GetH2(c, "coverage"), "col", 0,
	      "#eta coverage per v_{z}");
    TH2* cvst = GetH2(c, "nClusterVsnTracklet");
    if (cvst) {
      // cvst->Scale(1, "width");
      cvst->GetXaxis()->SetTitle("N_{free cluster}");
      cvst->GetYaxis()->SetTitle("N_{tracklet}");
      cvst->GetXaxis()->SetRangeUser(1,10000);
      cvst->GetYaxis()->SetRangeUser(1,10000);
    }
    DrawInPad(fBody, 3, cvst, "colz", kLogx|kLogy|kLogz,
	      "Correlation of # of tracklets and clusters"); 
    DrawInPad(fBody, 4, GetH2(c, "clusterPerTracklet"), "colz", 0x0,
	      "# clusters per tracklet vs #eta"); 
    ShowSliceFit(true, cvst, 3, fBody, 3, 0x8000|kLogz);

    fBody->cd(1)->Modified();
    fBody->cd(2)->Modified();
    fBody->cd(3)->Modified();
    fBody->cd(4)->Modified();
    fBody->cd(1)->Update();
    fBody->cd(2)->Update();
    fBody->cd(3)->Update();
    fBody->cd(4)->Update();
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
		TVirtualPad* can,
		Int_t        sub,
		Int_t        nCol)
  {
    if (all->GetHists()->GetEntries() <= 0 || !cur) return;

    // Info("", "Drawing step # %d", step);
    Bool_t       left = sub % nCol == 1; 
    Bool_t       right= sub % nCol == 0;
    Bool_t       top  = (sub-1) / nCol == 0;
    TVirtualPad* p    = can->cd(sub);
    gStyle->SetOptTitle(0);
    p->SetTitle(Form("Step # %d", step));
    p->SetFillColor(kWhite);
    p->SetRightMargin(right ? 0.02 : 0);
    p->SetTopMargin(top ? 0.02 : 0); // 0.02);
    // Info("", "Drawing step %d in sub-pad %d (%s)", 
    //      step, sub, (left?"left":"right"));

    p->cd();
    all->Draw("nostack");
    all->GetHistogram()->SetXTitle("#eta");
    all->GetHistogram()->SetYTitle("signal");

    TLegendEntry* e = 
      static_cast<TLegendEntry*>(leg->GetListOfPrimitives()->At(step-1));
    if (e) {
      e->SetMarkerColor(kBlack);
      e->SetLineColor(kBlack);
      e->SetTextColor(kBlack);
    }

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

    if (step > 1) { 
      Double_t x1 = 0.5*(p->GetUxmax()+p->GetUxmin());
      Double_t x2 = x1;
      Double_t y1 = p->GetUymax();
      Double_t y2 = 0.92*y1;
      Double_t sz = 0.05;
      if (fLandscape) { 
	x1 = 0.99*p->GetUxmin();
	x2 = 0.80*x1;
	y1 = .5*(p->GetUymax()+p->GetUymin());
	y2 = y1;
	sz = 0.034;
      }
      // Info("", "Arrow at (x1,y1)=%f,%f (x2,y2)=%f,%f", x1, y1, x2, y2);
      TArrow* a = new TArrow(x1, y1, x2, y2, sz, "|>");
      // (fLandscape ? "<|" : "|>"));
      a->SetFillColor(kGray+1);
      a->SetLineColor(kGray+1);
      a->Draw();
    }
    p->Modified();
    p->Update();
    p->cd();

    if (e) {
      e->SetMarkerColor(kGray);
      e->SetLineColor(kGray);
      e->SetTextColor(kGray);
    }
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

    THStack* esds    = GetStack(GetCollection(fResults, "fmdSharingFilter"), 
				"sumsESD", "summedESD");
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

    FixStack(esds,   "#sum_{s} #Delta/#Delta_{mip}", "",     20);
    FixStack(deltas, "#sum_{c} #Delta/#Delta_{mip}", "",     21);
    FixStack(nchs,   "#sum_{b} N_{ch,incl}", 	     "",     22);
    FixStack(prims,  "#sum_{b} N_{ch,primary}",      "",     23);
    FixStack(rings,  "dN/d#eta per ring",            "",     33);
    FixStack(mcRings,"dN/d#eta per ring (MC)",       "(MC)", 34);

    THStack* all = new THStack;
    AddToAll(all, mcRings);
    AddToAll(all, esds);
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

    if (esds) {
      h = static_cast<TH1*>(esds->GetHists()->At(0));
      AddLegendEntry(l, h, esds->GetTitle());
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
			esds, 
			deltas, 
			nchs, 
			prims, 
			rings, 
			dndeta };
    const char* titles[] = { /* 1 */ "MC",  
			     /* 2 */ "ESD input",
			     /* 3 */ "After merging", 
			     /* 4 */ "After particle counting", 
			     /* 5 */ "After corrections", 
			     /* 6 */ "After normalization", 
			     /* 7 */ "After combining" };
    Int_t nY = mcRings ? 4 : 3;
    Int_t nX = 2;
    if (fLandscape) {
      Int_t tmp = nX;
      nX        = nY;
      nY        = tmp;
    }
    fBody->Divide(nX, nY, 0, 0);
    
    Int_t step = 0;
    for (Int_t i = 0; i < 7; i++) { 
      TObject* obj = objs[i];
      if (!obj) continue;

      step++;
      Int_t padNo = step;
      if (!fLandscape) {
	switch (step) { 
	case 1: padNo = 1; break; 
	case 2: padNo = 3; break; 
	case 3: padNo = 5; break; 
	case 4: padNo = (mcRings ? 7 : 2); break; 
	case 5: padNo = (mcRings ? 2 : 4); break; 
	case 6: padNo = (mcRings ? 4 : 6); break; 
	case 7: padNo = (mcRings ? 6 : 8); break; 
    }
      }
      //Printf("Drawing step %d in sub-pad %d (%s)",step,padNo,obj->GetTitle());
      DrawStep(step, all, obj, l, titles[i], fBody, padNo, nX);
    }

    if (!esds && !mcRings && deltas) { 
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
    // Printf("Done drawing steps");
    PrintCanvas("Steps");
  }


  //____________________________________________________________________
  void DrawResults()
  {
    // MakeChapter(can, "Results");

    fBody->Divide(2,1);

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
    const char** pring   = GetRingNames(false);
    
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
    Double_t savX = fParVal->GetX();
    Double_t savY = fParVal->GetY();
    fParVal->SetX(.3);
    fParVal->SetY(.2);
    TVirtualPad* p = fBody->cd(1);
    p->Divide(1,2,0,0);
    DrawInPad(p, 1, GetStack(c, "all"), "nostack", mcRings ? 0 : kLegend,
	      "Individual ring results");
    DrawInPad(p, 1, mcRings, "nostack same", kLegend|kSilent);
    DrawInPad(p, 2, allEta, "nostack hist", kLegend,
	      "#phi acceptance and #eta coverage per ring");
    DrawInPad(p, 2, allPhi, "nostack hist same", 0x0);

    p = fBody->cd(2);
    p->Divide(1,2,0,0);
    DrawInPad(p, 1, dndeta_phi, "", 0x0, 
	      "1/#it{N}_{ev} d#it{N}_{ch}/d#it{#eta}");
    DrawInPad(p, 1, dndeta_eta, "Same", kLegend);
    DrawInPad(p, 2, GetH1(fResults, "norm"), "", 0x0, 
	      "Total #phi acceptance and #eta coverage");
    DrawInPad(p, 2, GetH1(fResults, "phi"), "same", kLegend);
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
    
    fParVal->SetX(savX);
    fParVal->SetY(savY);
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
    TH1* h = stack->GetHistogram();
    if (h) {
      h->SetXTitle("#it{#eta}");
      h->SetYTitle("#frac{d#it{N}_{ch}}{d#it{#eta}}");
    }
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
