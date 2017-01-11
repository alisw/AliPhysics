#include "SummaryDrawer.C"
#ifndef __CINT__
# include <TGraph.h>
# include <TGraphErrors.h>
# include <TF1.h>
# include <TArrow.h>
#else
class TGraph;
class TFile;
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
    kESDFixer          = 0x100,
    kNormal            = 0x1FF
  };
  SummaryAODDrawer() 
    : SummaryDrawer(),
      fSums(0),
      fResults(0)
  {}
  //__________________________________________________________________
  TFile* Init(const char* fname)
  {
    // --- Open the file ---------------------------------------------
    TString filename(fname);
    TFile*  file = TFile::Open(filename, "READ");
    if (!file) { 
      Error("Run", "Failed to open \"%s\"", filename.Data());
      return 0;
    }

    // --- Get top-level collection ----------------------------------
    fSums = GetCollection(file, "ForwardSums");
    if (!fSums) {
      Info("Run", "Trying old name Forward");
      fSums = GetCollection(file, "Forward");
      if (!fSums) return 0;
    }

    // --- Do the results ----------------------------------------------
    fResults = GetCollection(file, "ForwardResults");
    if (!fResults) fResults = fSums; // Old-style

    return file;
  }
  void SummarizeSharing(const char* fname, UShort_t what=0)
  {
    // --- Initialize ------------------------------------------------
    TFile* file = 0;
    if (!(file = Init(fname))) return;

    // --- Make our canvas -------------------------------------------
    TString pdfName("sharing.pdf");
    pdfName.ReplaceAll(".root", ".pdf");
    CreateCanvas(pdfName, what & kLandscape, true, false);
    
    TCollection* c = GetCollection(fSums, "fmdSharingFilter");
    if (!c) return;
    TCollection* rc = GetCollection(fResults, "fmdSharingFilter");
    if (!rc) rc = c;

    Int_t    nFiles = 0;
    TParameter<int>* pnFiles = 
      static_cast<TParameter<int>*>(GetObject(c, "nFiles"));
    if (pnFiles) {
      nFiles = pnFiles->GetVal();
    }

    TCollection* lc    = GetCollection(c, "lCuts");
    TCollection* hc    = GetCollection(c, "hCuts");
    Int_t        lm    = 0;
    Int_t        hm    = 0;
    TH2*         hLow  = GetH2(c, "lowCuts");
    TH2*         hHigh = GetH2(c, "highCuts");
    GetParameter(lc, "method", lm);
    GetParameter(hc, "method", hm);
    if (hLow  && nFiles > 0 && !hLow->TestBit(BIT(20)))
      hLow->Scale(1. / nFiles);
    if (hHigh && nFiles > 0 && !hHigh->TestBit(BIT(20))) 
      hHigh->Scale(1. / nFiles);
    
    DivideForRings(true,true);
    const char** ptr   = GetRingNames(false);
    UShort_t     iq    = 1;
    while (*ptr) { 
      TCollection* sc = GetCollection(c, *ptr);
      if (!sc) { ptr++; iq++; continue; }
      UShort_t d = Int_t((*ptr)[3])-48;
      Char_t   r = (*ptr)[4];
      
      TH1*     esdELoss = GetH1(sc, "esdEloss");
      TH1*     anaELoss = GetH1(sc, "anaEloss");
      TGraph*  lowCut   = CreateCutGraph(lm, iq,  hLow,  esdELoss, kYellow+1);
      TGraph*  highCut  = CreateCutGraph(hm, iq,  hHigh, esdELoss, kCyan+1);
      // Double_t ignCut   = TMath::Max(lowCut->GetX()[3],0.);
      // Int_t    esdLow   = esdELoss->FindBin(ignCut);
      // Int_t    anaLow   = anaELoss->FindBin(ignCut);
      // Double_t esdInt   = esdELoss->Integral(esdLow,esdELoss->GetNbinsX()+1);
      // Double_t anaInt   = anaELoss->Integral(anaLow,anaELoss->GetNbinsX()+1);
      // Double_t frac     = esdInt > 0 ? (esdInt-anaInt)/esdInt : 1;
      esdELoss->GetXaxis()->SetRangeUser(-.1, 2);
      
      DrawInRingPad(d,r, esdELoss, "", kLogy, *ptr);
      // 		    "#Delta/#Delta_{mip} reconstructed and merged");
      DrawInRingPad(d, r, anaELoss, "same");
      DrawInRingPad(d, r, lowCut,  "lf same"); 
      DrawInRingPad(d, r, highCut, "lf same"); 
      ptr++;
      iq++;
    }
    TVirtualPad* p = RingPad(0,0);
    p->cd();
    TLegend* l = new TLegend(0.1, 0.1, 0.98, 0.98, "");
    l->SetFillStyle(0);
    l->SetFillColor(0);
    l->SetBorderSize(0);
    TLegendEntry* e = 0;
    e = l->AddEntry("dummy", "ESD signal", "f");
    e->SetFillStyle(3002);
    e->SetFillColor(kBlack);
    e = l->AddEntry("dummy", "Merged signal", "f");
    e->SetFillStyle(3001);
    e->SetFillColor(kBlack);
    e = l->AddEntry("dummy", "Low cut", "f");
    e->SetFillStyle(3002);
    e->SetFillColor(kYellow+1);
    e->SetLineWidth(0);
    e->SetLineColor(kWhite);
    e = l->AddEntry("dummy", "High cut", "f");
    e->SetFillStyle(3002);
    e->SetFillColor(kCyan+1);
    e->SetLineWidth(0);
    e->SetLineColor(kWhite);
    l->Draw();
    
    PrintCanvas("Summary of sharing filter");
    CloseCanvas();
    
  }
  void SummarizeSteps(const char* fname, UShort_t what=0)
  {
    // --- Initialize ------------------------------------------------
    TFile* file = 0;
    if (!(file = Init(fname))) return;

    // --- Make our canvas -------------------------------------------
    TString pdfName("steps.pdf");
    pdfName.ReplaceAll(".root", ".pdf");
    CreateCanvas(pdfName, what & kLandscape, true, false);
    DrawSteps();
    CloseCanvas();
  }
  
  //__________________________________________________________________
  /** 
   * 
   * 
   * @param fname 
   * @param what 
   */
  void Run(const char* fname, UShort_t what=kNormal)
  {
    // --- Initialize ------------------------------------------------
    TFile* file = 0;
    if (!(file = Init(fname))) return;

    // --- Make our canvas -------------------------------------------
    TString pdfName(fname);
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

    // Plot status if found 
    TH1* hStatus = GetH1(fSums, "status", false);
    if (hStatus) { 
      hStatus->SetMaximum(hStatus->GetMaximum()*1.2);
      fBody->SetRightMargin(0.10);
      DrawInPad(fBody,0,hStatus, "hist text30");
      PrintCanvas("Status");
    }

    // --- Do each sub-algorithm -------------------------------------
    if (what & kEventInspector)    DrawEventInspector(fSums);
    if (what & kESDFixer)          DrawESDFixer(fSums);
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
    Int_t    iiy = (iy == 4 ? 5 : iy == 5 ? 4 : iy);
    if (method == 0) { // Fixed value
      max = cuts->GetBinContent(1, iiy);
      min = eloss->GetXaxis()->GetXmin();
    }
    else {
      for (Int_t ix=1; ix <= cuts->GetNbinsX(); ix++) {
	Double_t c = cuts->GetBinContent(ix, iiy);
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
    case 0: return "c=X";
    case 1: return "c=X#times#Delta_{p}";
    case 2: return "c:Lower bound of fit range";
    case 3: return "c=#Delta_{p}-X#times#xi";
    case 4: return "c=#Delta_{p}-X#times(#xi+#sigma)";
    case 5: return "c:P(#Delta<c)<X";
    case 6: return "c:#Delta_{p}-X#times#bar#xi+#sigma}";
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
    DrawParameter(y, name, CutMethodName(method), size);

    TString params;
    const char*  cuts[] = { "fmd1i", "fmd2i", "fmd2o", "fmd3i", "fmd3o", 0 };
    const char** pcut   = cuts;
    while (*pcut) { 
      Double_t cut;
      GetParameter(c, *pcut, cut, false);
      if (pcut != cuts) params.Append(", ");
      params.Append(Form("%5.2f", cut));
      pcut++;
    }
    if (params.IsNull()) {
      Double_t frac = 0;
      GetParameter(c, "frac", frac);
      params = Form("%f", frac);
    }
    DrawParameter(y, "Parameters", params, size);
    return method;
  }
  //____________________________________________________________________
  void DrawCut(TVirtualPad* parent, Int_t sub, TH2* cuts)
  {
    if (!cuts) return;
    THStack* stack = new THStack(cuts,"x");
    stack->SetTitle(cuts->GetTitle());
    for (Int_t i = 1; i <= cuts->GetNbinsY(); i++) {
      TH1*     hist = static_cast<TH1*>(stack->GetHists()->At(i-1));
      TString  name(cuts->GetYaxis()->GetBinLabel(i));
      UShort_t det = UShort_t(name[3]-48);
      Char_t   rng = name[4];
      Color_t  col = RingColor(det, rng);
      hist->SetDirectory(0);
      hist->SetTitle(name);
      hist->SetMarkerStyle(20);
      hist->SetMarkerColor(col);
      hist->SetLineColor(col);
      hist->SetFillColor(col);
      hist->SetLineWidth(0);
      hist->SetFillStyle(0);
      hist->SetXTitle("#eta");
      hist->SetYTitle(cuts->GetZaxis()->GetTitle());
    }
    DrawInPad(parent, sub, stack, "nostack p", kLegend|kCenter|kSouth);
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
  
    Double_t y = .95;
    Bool_t   angle=false, lowSignal=false, disabled=false;
    Int_t    nFiles = 0;
    if (GetParameter(c, "angle", angle))
      DrawParameter(y, "Angle correct", (angle ? "yes" : "no")); 
    if (GetParameter(c, "lowSignal", lowSignal))
      DrawParameter(y, "Lower signal",  (lowSignal ? "yes" : "no"));
    TParameter<int>* pnFiles = 
      static_cast<TParameter<int>*>(GetObject(c, "nFiles"));
    if (pnFiles) {
      nFiles = pnFiles->GetVal();
      DrawParameter(y, "# files merged", Form("%d", nFiles));    
    }
    if (GetParameter(c, "disabled", disabled, false)) 
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
      if (hLow  && nFiles > 0 && !hLow->TestBit(BIT(20)))
	hLow->Scale(1. / nFiles);
      if (hHigh && nFiles > 0 && !hHigh->TestBit(BIT(20))) 
	hHigh->Scale(1. / nFiles);
      DrawCut(fBody, 2, hLow);
      DrawCut(fBody, 3, hHigh);
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
	TH1*    anaELoss = GetH1(sc, "anaEloss");
	TGraph* lowCut   = CreateCutGraph(lm, iq,  hLow,  esdELoss, kYellow+1);
	TGraph* highCut  = CreateCutGraph(hm, iq,  hHigh, esdELoss, kCyan+1);
	Double_t ignCut  = TMath::Max(lowCut->GetX()[3],0.);
	Int_t    esdLow  = esdELoss->FindBin(ignCut);
	Int_t    anaLow  = anaELoss->FindBin(ignCut);
	Double_t esdInt  = esdELoss->Integral(esdLow,esdELoss->GetNbinsX()+1);
	Double_t anaInt  = anaELoss->Integral(anaLow,anaELoss->GetNbinsX()+1);
	Double_t frac    = esdInt > 0 ? (esdInt-anaInt)/esdInt : 1;
	esdELoss->GetXaxis()->SetRangeUser(-.1, 2);
	
	DrawInPad(fBody, 1, esdELoss, "", kLogy,
		  "#Delta/#Delta_{mip} reconstructed and merged");
	DrawInPad(fBody, 1, anaELoss, "same");
	DrawInPad(fBody, 1, lowCut,  "lf same"); 
	DrawInPad(fBody, 1, highCut, "lf same", kLogy|kLegend|kNorth|kWest); 
	TVirtualPad* p = fBody->GetPad(1);
	p->cd();
	TLatex* l = new TLatex(1-p->GetRightMargin(), 
			       0.5, Form("Loss: %5.1f%%", frac*100));
	l->SetNDC();
	l->SetTextAlign(32);
	l->Draw();
	l->DrawLatex(1-p->GetRightMargin(), 0.45,
		     Form("%f #rightarrow #infty", ignCut));

	TH1*     singles  = GetH1(sc, "singleEloss");
	TH1*     doubles  = GetH1(sc, "doubleEloss");
	TH1*     tripples = GetH1(sc, "tripleEloss");
	Double_t int1     = singles->Integral(0,singles->GetNbinsX()+1);
	Double_t int2     = doubles->Integral(0,doubles->GetNbinsX()+1);
	Double_t int3     = tripples->Integral(0,tripples->GetNbinsX()+1);
	Double_t intT     = int1 + int2 + int3;
	Double_t f1       = intT > 0 ? int1 / intT : 0;
	Double_t f2       = intT > 0 ? int2 / intT : 0;
	Double_t f3       = intT > 0 ? int3 / intT : 0;

	singles->GetXaxis()->SetRangeUser(-.1, 2);
	DrawInPad(fBody, 2, singles,    "",    kLogy,
		  "#Delta/#Delta_{mip} for single, double, and tripple hits");
	DrawInPad(fBody, 2, doubles,    "same",    kLogy);
	DrawInPad(fBody, 2, tripples,   "same",    kLogy);
	DrawInPad(fBody, 2, lowCut,     "lf same", kLogy); 
	DrawInPad(fBody, 2, highCut,    "lf same", kLogy|kLegend|kNorth|kWest); 
	
	fBody->cd(2);
	Double_t nameX = fParName->GetX();
	Double_t valX  = fParVal->GetX();
	Double_t intY  = 0.4;
	fParName->SetX(0.5);
	fParVal->SetX(0.7);
	DrawParameter(intY, "Singles",  Form("%5.1f%%", 100*f1), 0.05);
	DrawParameter(intY, "Doubles",  Form("%5.1f%%", 100*f2), 0.05);
	DrawParameter(intY, "Tripples", Form("%5.1f%%", 100*f3), 0.05);
	fParName->SetX(nameX);
	fParVal->SetX(valX);

	DrawInPad(fBody, 3, GetH2(sc, "singlePerStrip"), "colz",kLogz);
	// DrawInPad(fBody, 4, GetH1(sc, "distanceBefore"), "",     0x2);
	// DrawInPad(fBody, 4, GetH1(sc, "distanceAfter"),  "same", 0x12);
	DrawInPad(fBody, 4, GetH2(sc, "summed"),         "colz",0x0);
	
	TH2* nB = GetH2(sc, "neighborsBefore");
	if (nB) { 
	  nB->GetXaxis()->SetRangeUser(0,2); 
	  nB->GetYaxis()->SetRangeUser(0,2); 
	}
	DrawInPad(fBody, 5, nB, "colz cont3", kLogz,
		  "Correlation of neighbors before merging");

	TH2* nA = GetH2(sc, "neighborsAfter");
	if (nA) { 
	  nA->GetXaxis()->SetRangeUser(0,2); 
	  nA->GetYaxis()->SetRangeUser(0,2); 
	}
	DrawInPad(fBody, 6, nA, "colz cont3", kLogz,
		  "Correlation of neighbors after merging");

	// DrawInPad(fBody, 6, GetH2(sc, "beforeAfter"),    "colz",   kLogz);
	
	PrintCanvas(Form("Sharing filter - %s", *ptr));
	ptr++;
	iq++;
      }
      fParVal->SetX(savX);
      fParVal->SetY(savY);
    }

    // --- MC --------------------------------------------------------
    TCollection* cc = GetCollection(c, "esd_mc_comparion", false); // Spelling!
    if (!cc) {
      cc = GetCollection(c, "esd_mc_comparison", false); // Spelling!
      if (!cc) return; // Not MC 
    }

    DivideForRings(false, false);
    const char** ptr = GetRingNames(false);
    while (*ptr) { 
      TString nam(Form("%s_corr", *ptr));
      TH2* hc = GetH2(cc, nam, false);
      if (!hc) {
	nam[4] = (nam[4] == 'I' ? 'i' : 'o');
	hc = GetH2(cc, nam, false);
	if (!hc) { 
	  ptr++;
	  continue;
	}
      }
      DrawInRingPad(hc, "colz", kLogz);
      ptr++;
    }

    PrintCanvas("Sharing filter - MC vs Reco");

    // --- MC --------------------------------------------------------
    DrawTrackDensity(c);
  }
  //__________________________________________________________________
  /** 
   * Draw a slice fit on a 2D histogram
   * 
   * @param inY   Whether to slice in Y
   * @param h     2D histogram
   * @param nVar  Number of variances
   * @param p     Master pad to draw in 
   * @param sub   Sub pad number 
   * @param flags Flags 
   * @param cut   Cut value 
   */
  void ShowSliceFit(Bool_t inY, TH2* h, Double_t nVar, 
		    TVirtualPad* p, Int_t sub, UShort_t flags=0,
		    Double_t cut=-1)
  {
    if (!h) return;
    TAxis* indep = (inY ? h->GetXaxis() : h->GetYaxis());
    TAxis* dep   = (inY ? h->GetYaxis() : h->GetXaxis());
    

    TObjArray* fits = new TObjArray;
    fits->SetOwner();
    // FitSlicesX sets onX to true.  If onX is true, then it calls
    // ProjectionX for each Y bin.  That is, we fit each X slice. 
    Int_t minB = dep->FindBin(0.)+1;
    Int_t maxB = dep->GetNbins();
    Info("", "Bin range: %d - %d", minB, maxB);
    if (inY) h->FitSlicesY(0, minB, maxB, 10, "QN", fits);
    else     h->FitSlicesX(0, minB, maxB, 10, "QN", fits);
    if (!fits) { 
      Warning("ShowSliceFit", "No fits returned");
      return;
    }

    // If inY is true then this is the mean,variance as a function of X
    TH1* mean = static_cast<TH1*>(fits->At(1));
    TH1* var  = static_cast<TH1*>(fits->At(2));
    if (!mean || !var) {
      Warning("ShowSliceFit", "Didn't get histograms");
      fits->Delete();
      return;
    }
#if 0
    TH1* hh[] = { mean, var, 0 };
    for (Int_t ib=1; ib<=mean->GetNbinsX();ib++) { 
      TH1** hp = hh;
      while (*hp) { 
	if ((*hp)->GetBinContent(ib) <= 0) { 
	  (*hp)->SetBinContent(ib,0);
	  (*hp)->SetBinError(ib,0);
	}
	hp++;
      }
    }
#endif

    // If inY is true then this is the mean,variance as a function of X
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
    TString xTit = indep->GetTitle();
    TString yTit = dep->GetTitle();
    g->SetTitle(Form("Correlation of %s and %s",xTit.Data(),yTit.Data()));
    g->SetFillColor(kBlue-10);
    g->SetFillStyle(3001);
    TGraph* up  = (cut > 0 ? new TGraph(indep->GetNbins()) : 0);
    TGraph* low = (cut > 0 ? new TGraph(indep->GetNbins()) : 0);
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
      Double_t x  = indep->GetBinCenter(i);
      Double_t y  = fmean->Eval(x);
      Double_t e  = fvar->Eval(x);
      // Double_t y  = mean->GetBinContent(i);
      // Double_t e  = var->GetBinContent(i); 
      Double_t ee = nVar * e;
      if (flags & 0x8000) ee *= e > 0 ? TMath::Log10(e) : 1;
      if (inY) {
	g->SetPoint(i-1, x, y);
	g->SetPointError(i-1, 0, ee);
      }
      else { 
	g->SetPoint(i-1, y, x);
	g->SetPointError(i-1, ee, 0);
      }
      if (up)  {
	if (inY) up->SetPoint(i-1,x,x+cut*x);
	else     up->SetPoint(i-1,y+cut*y,y);
      }
      if (low) {
	if (inY) low->SetPoint(i-1,x,x-cut*x);
	else     low->SetPoint(i-1,y-cut*y,y);
      }
    }
    DrawInPad(p, sub, g, "3", flags);
    if (up)  DrawInPad(p, sub, up, "l", flags);
    if (low) DrawInPad(p, sub, low, "l", flags);
    fmean->SetRange(indep->GetXmin(), indep->GetXmax());
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
    l->SetTextSize(fParVal->GetTextSize());
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
  
    Double_t y = .9;
    Int_t maxParticles=0, phiAcceptance=0, etaLumping=0, phiLumping=0;
    Bool_t method=false, recalcPhi=false;
    Double_t maxOutliers=0, outlierCut=-1;
    Double_t size = fLandscape ? 0.05 : 0.03;
  
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
    if (GetParameter(c, "recalcPhi", recalcPhi))
      DrawParameter(y, "Recalculate #phi",(recalcPhi ? "yes" : "no"),size); 
    if (GetParameter(c, "maxOutliers", maxOutliers, false))
      DrawParameter(y, "Max relative N_{outlier}",
		    Form("%5.3f",maxOutliers),size);
    Bool_t hasOutCut = false;
    if ((hasOutCut = GetParameter(c, "outlierCut", outlierCut, false)))
      DrawParameter(y, "Max relative deviation",Form("%5.3f",outlierCut),size);


    TParameter<int>* nFiles = 
      static_cast<TParameter<int>*>(GetObject(c, "nFiles"));
    if (nFiles)
      DrawParameter(y, "# files merged", Form("%d", nFiles->GetVal()), size);

    TCollection* lc = GetCollection(c, "lCuts");
    Int_t tm = PrintCut(lc, y, "Threshold", size);

    TVirtualPad* p = fBody; // fBody->cd(2);
    // p->Divide(3,1);

    TH1* accI = GetH1(c, "accI");
    TH1* accO = GetH1(c, "accO");
    if (accI) { 
      Double_t scale = 1./accI->GetMaximum();
      accI->Scale(scale); 
      accO->Scale(scale);
      accI->SetMinimum(0); 
      accI->SetMaximum(1.3);
    }
    TH2* lCuts = GetH2(c, "lowCuts");
    TH2* maxW  = GetH2(c, "maxWeights");
    if (lCuts)           lCuts->SetTitle("Thresholds");
    if (nFiles && lCuts) lCuts->Scale(1. / nFiles->GetVal());
    if (nFiles && maxW)  maxW->Scale(1. / nFiles->GetVal());
    DrawInPad(p, 2, accI); 
    DrawInPad(p, 2, accO,  "same", kLegend|kNorth|kCenter); 
    DrawCut(p, 3, lCuts);
    DrawCut(p, 4, maxW);
  
    PrintCanvas("Density calculator");
    Float_t save = fParVal->GetTextSize();
    fParVal->SetTextSize(0.03);

    UShort_t iq = 1;
    const char** ptr   = GetRingNames(false);
    while (*ptr) { 
      TCollection* sc = GetCollection(c, *ptr);
      if (!sc) { ptr++; continue; }
    
      if (fLandscape) fBody->Divide(3,2);
      else            fBody->Divide(2,3);
    
      TH2* corr      = GetH2(sc, "elossVsPoisson");
      TH2* corrOut   = GetH2(sc, "elossVsPoissonOutlier", false);
      TH1* diff      = GetH1(sc, "diffElossPoisson");
      TH1* diffOut   = GetH1(sc, "diffElossPoissonOutlier", false);
      TH1* eloss     = GetH1(sc, "eloss");
      TH1* elossUsed = GetH1(sc, "elossUsed");
      TH1* occ       = GetH1(sc, "occupancy");
      if (eloss)     eloss    ->SetLineWidth(1);
      if (elossUsed) elossUsed->SetLineWidth(1);
      if (eloss)     eloss->GetXaxis()->SetRangeUser(0.05, 2);
      
      corr->SetXTitle("N_{ch,#Delta}");
      corr->SetYTitle("N_{ch,Poisson}");
      DrawInPad(fBody,  1, corr,    "colz",        kLogz);
      if (corrOut) 
	DrawInPad(fBody,1, corrOut, "same",        kLogz);
      ShowSliceFit(true, corr, 10, fBody, 1, kLogz, outlierCut);

      DrawInPad(fBody,  2, diff,    "HIST E",      kLogy);
      if (diffOut) 
	DrawInPad(fBody,2, diffOut, "HIST E SAME", kLogy|kLegend|kNorth|kWest);
      DrawInPad(fBody,  3, occ,      "",           kLogy);
      DrawInPad(fBody,  4, eloss,    "",           kLogy, 
		"#Delta/#Delta_{mip} before and after cuts");
      DrawInPad(fBody, 4, elossUsed, "same",      kLogy);
      TGraph* thres = CreateCutGraph(tm, iq,  lCuts,  eloss, kYellow+1);
      DrawInPad(fBody, 4, thres, "lf same", kLogy|kLegend|kNorth|kWest); 

      if (eloss && elossUsed) {
	Int_t    lowBin    = eloss->GetXaxis()->FindBin(0.)+1;
	Int_t    upBin     = eloss->GetNbinsX()+1;
	Double_t beforeInt = eloss->Integral(lowBin,upBin);
	Double_t afterInt  = elossUsed->Integral(lowBin,upBin);
	Double_t frac      = beforeInt > 0 ? (beforeInt-afterInt)/beforeInt : 1;
	TVirtualPad* pp = fBody->GetPad(4);
	pp->cd();
	TLatex* l = new TLatex(1-pp->GetRightMargin(), 
			       0.5, Form("Loss: %5.1f%%", frac*100));
	l->SetNDC();
	l->SetTextAlign(32);
	l->Draw();
      }


      TH1* phiB = GetH1(sc, "phiBefore");
      TH1* phiA = GetH1(sc, "phiAfter");
      TH1* outliers = GetH1(sc, "outliers", false);
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
      iq++;
    }

    TCollection* cc = GetCollection(c, "esd_mc_comparison", false); 
    if (!cc) {
      fParVal->SetTextSize(save);
      return; // Not MC 
    }

    ptr   = GetRingNames(false);
    THStack* profiles = new THStack("profiles", "MC-ESD");
    // Int_t cnt = 0;
    while (*ptr) { 
      fBody->Divide(2,1);

      TH2* corr = GetH2(cc, Form("%s_corr_mc_esd", *ptr));
      if (corr) { 
	// corr->GetXaxis()->SetRangeUser(-1,51);
	// corr->GetYaxis()->SetRangeUser(-1,51);
	corr->SetXTitle("MC");
	corr->SetYTitle("ANA");
	TH1* h = new TH1F("",Form("Correlation of N_{ch,incl} for %s", *ptr),
			  corr->GetNbinsX(), 
			  corr->GetXaxis()->GetXmin(), 
			  corr->GetXaxis()->GetXmax());
	h->SetMinimum(corr->GetYaxis()->GetXmin());
	h->SetMaximum(1.2*corr->GetYaxis()->GetXmax());
	h->SetXTitle("MC"); // corr->GetXaxis()->GetTitle());
	h->SetYTitle(corr->GetYaxis()->GetTitle());
	DrawInPad(fBody, 1, h, "", kLogz);
	DrawInPad(fBody, 1, corr,  "colz same", kLogz);
	ShowSliceFit(true, corr, 1, fBody, 1, 0, -1);
      }
      TH1* diff = GetH1(cc, Form("%s_diff_mc_esd", *ptr));
      if (diff && diff->GetEntries() > 0) {
	// diff->Rebin(2);
	diff->Smooth(3);
	diff->Scale(1./diff->GetEntries(), "width");
	TF1* f = new TF1(Form("%s_tri", *ptr), 
			 "[0]*TMath::Exp(-TMath::Abs(x-[1])/[2])",-20,20);
	f->SetParNames("A", "x_{0}", "w");
	f->SetParameters(1, diff->GetMean(), diff->GetRMS());
	diff->Fit(f,"RQ0+");
	DrawInPad(fBody, 2, diff, "", kLogy);
	DrawInPad(fBody, 2, f, "same", kLogy);
	Double_t py = .88;
	TLatex* l = new TLatex(.89, py, "Ae^{-|x-x_{0}|/w}");
	l->SetTextAlign(33);
	l->SetTextSize(fParVal->GetTextSize());
	l->SetTextFont(42);
	l->SetNDC();
	l->Draw();
	py -= fParVal->GetTextSize();
	l->DrawLatex(.89, py, Form("#chi^{2}/#nu=%f",  f->GetNDF()>0 ? 
				   f->GetChisquare()/f->GetNDF() : 0));
	for (Int_t i = 0; i < 3; i++) { 
	  py -= fParVal->GetTextSize();
	  l->DrawLatex(.89, py, Form("%s: %f #pm %f", f->GetParName(i), 
				     f->GetParameter(i), f->GetParError(i)));
	}	       
      }
      TH2* vs = GetH2(cc, Form("%s_esd_vs_mc", *ptr));
      if (vs) { 
	const Double_t lVs = -0.3;
	const Double_t hVs = +0.4;
	Int_t nX = vs->GetNbinsX();
	Int_t nY = vs->GetNbinsY();
	TProfile* vsp = new TProfile(Form("%s_esd_vs_mc_px", *ptr), *ptr,
				   nX, vs->GetXaxis()->GetXmin(),
				   vs->GetXaxis()->GetXmax());
	vsp->SetXTitle(vs->GetXaxis()->GetTitle());
	vsp->SetYTitle("#LT#deltaN_{ch}#GT");
	vsp->SetLineColor(diff->GetLineColor());
	vsp->SetMarkerColor(diff->GetLineColor());
	vsp->SetFillColor(diff->GetLineColor());
	vsp->SetDirectory(0);
	Bool_t hasSome = false;
	for (Int_t ix = 1; ix <= nX; ix++) { 	
	  Double_t vsx = vs->GetXaxis()->GetBinCenter(ix);
	  for (Int_t iy = 1; iy <= nY; iy++) { 
	    Double_t vsc = vs->GetBinContent(ix, iy);
	    Double_t vse = vs->GetBinError(ix, iy);
	    if (vsc < lVs || vsc > hVs) continue;
	    hasSome = true;
	    vsp->Fill(vsx, vsc, vse);
	  }
	}
	if (hasSome) profiles->Add(vsp);
	else delete vsp;
      }
      PrintCanvas(Form("Density calculator - MC vs Reco - %s", *ptr));
      ptr++;
      // cnt++;
    }
    if (profiles->GetHists() && profiles->GetHists()->GetEntries() > 0) {
      Double_t pmax = profiles->GetMaximum("nostack")*1.3;
      profiles->SetMaximum(+pmax);
      profiles->SetMinimum(-pmax);
      DrawInPad(fBody, 0, profiles, "nostack");
      PrintCanvas("Density calculator - MC vs Reco - All");
    }
      
    fParVal->SetTextSize(save);
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
    TH1* skipped = GetH1(c, "skipped", false);
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
				"sumsESD", "summedESD",false);
    THStack* deltas  = GetStack(GetCollection(fResults, "fmdSharingFilter"), 
				"sums", "summed",false);
    THStack* nchs    = GetStack(GetCollection(fResults, 
					      "fmdDensityCalculator"), 
				"sums", "inclDensity",false);
    THStack* prims   = GetStack(GetCollection(fResults, "fmdCorrector"), 
				"sums", "primaryDensity",false);
    THStack* rings   = GetStack(GetCollection(fResults, "ringResults"), 
				"all",0, false);
    THStack* mcRings = GetStack(GetCollection(fResults, "mcRingResults", false),
				"all","dndeta_eta", false);
    TH1*     dndeta  = GetH1(fResults, "dNdeta", false);
    if (dndeta) dndeta->SetMarkerColor(kBlack);
    if (!(esds || deltas || nchs || prims || rings)) return;

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

    THStack* stacks[] = { mcRings, 
			  esds, 
			  deltas, 
			  nchs, 
			  prims, 
			  rings };
    
    Int_t nHist = 0;
    for (Int_t i = 0; i < 6; i++) {
      if (!stacks[i]) continue;
      TH1* h = static_cast<TH1*>(stacks[i]->GetHists()->At(0));
      AddLegendEntry(l, h, stacks[i]->GetTitle());
      nHist++;
    }
    if (res) {
      AddLegendEntry(l, res, res->GetTitle());
      nHist++;
    }
    
    TObject* objs[] = { stacks[0], 
			stacks[1], 
			stacks[2], 
			stacks[3], 
			stacks[4], 
			stacks[5], 
			dndeta };
    const char* titles[] = { /* 1 */ "MC",  
			     /* 2 */ "ESD input",
			     /* 3 */ "After merging", 
			     /* 4 */ "After particle counting", 
			     /* 5 */ "After corrections", 
			     /* 6 */ "After normalization", 
			     /* 7 */ "After combining" };
    Int_t nY = nHist > 6 ? 4 : 3;
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
    if (!f1 || !c1 || !f2 || !c2) return; 
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
