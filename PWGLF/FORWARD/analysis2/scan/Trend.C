#ifndef __TREND_C_
#define __TREND_C_
#include "SummaryDrawer.C"
#include <TH1.h>
#include <THStack.h>
#include <TMultiGraph.h>
#include <TGraphAsymmErrors.h>
#include <TTree.h>
#include <TList.h>
#include <TLegend.h>

/**
 * Structure to make trendind plots from cut scans 
 * 
 */
struct Trend : public SummaryDrawer
{
  //__________________________________________________________________
  /** 
   * Constructor
   */
  Trend() 
  {
    fSLCuts.SetName("sl");
    fSHCuts.SetName("sh");
    fDCCuts.SetName("dc");
    fOrder[0] = &fSLCuts;
    fOrder[1] = &fSHCuts;
    fOrder[2] = &fDCCuts;
  }
  //=== Customization member function ================================
  /** 
   * Add a sharing filter low cut 
   * 
   * @param name    Name of the cut type (fix,mpv,sig,xi,prob)
   * @param values  Space-separated cut parameters 
   */
  void AddSLCut(const TString& name, const TString& values) 
  {
    fSLCuts.Add(new TNamed(name, values));
  }
  //__________________________________________________________________
  /** 
   * Add a sharing filter high cut 
   * 
   * @param name    Name of the cut type (fix,mpv,sig,xi,prob)
   * @param values  Space-separated cut parameters 
   */
  void AddSHCut(const TString& name, const TString& values) 
  {
    fSHCuts.Add(new TNamed(name, values));
  }
  //__________________________________________________________________
  /** 
   * Add a density calculator threshold 
   * 
   * @param name    Name of the cut type (fix,mpv,sig,xi,prob)
   * @param values  Space-separated cut parameters 
   */
  void AddDCCut(const TString& name, const TString& values) 
  {
    fDCCuts.Add(new TNamed(name, values));
  }
  //__________________________________________________________________
  /** 
   * Add a run to analyse 
   * 
   * @param name Run number
   */
  void AddRun(const TString& name)
  {
    fRuns.Add(new TNamed(name, name));
  }
  //__________________________________________________________________
  /** 
   * Add a centrality bin to analyse 
   * 
   * @param l  Lower bound
   * @param h  Upper bound
   */
  void AddCentrality(UShort_t l, UShort_t h)
  {
    TObjString* n = 0;
    if (l <= 0 && h >= 100) n = new TObjString("all");
    else                    n = new TObjString(Form("cent%03d_%03d", l, h));
    n->SetUniqueID((h & 0xFFFF) << 16 | (l & 0xFFFF));
    fCents.Add(n);
  }
  //__________________________________________________________________
  /** 
   * Set the order in which to do the comparison.  Most important
   * should be last.
   * 
   * @param order Permutation of sl, sh, and dc, space separated
   */
  void SetOrder(const TString& order) 
  {
    TObjArray* a = order.Tokenize(" ");
    if (a->GetEntries() < 3) {
      Error("SetOrder", "Order must contain a pertubation of "
	    "\"sl\", \"sh\", and \"dc\", separated by spaces");
      fOrder[0] = fOrder[1] = fOrder[2] = 0;
      return;
    }
    for (Int_t i = 0; i < 3; i++) { 
      const TString& n = static_cast<TObjString*>(a->At(i))->String();
      TList* l = 0;
      if      (n.EqualTo("sl", TString::kIgnoreCase)) l = &fSLCuts;
      else if (n.EqualTo("sh", TString::kIgnoreCase)) l = &fSHCuts;
      else if (n.EqualTo("dc", TString::kIgnoreCase)) l = &fDCCuts;
      if (!l) { 
	Error("SetOrder", "Unknown cut \"%s\"", n.Data());
	fOrder[0] = fOrder[1] = fOrder[2] = 0;
	return;
      }
      fOrder[i] = l;
    }
    a->Delete();
  }
  //=== Processing member functions ==================================
  /** 
   * Run it
   * 
   * @param output
   * @param sys 
   * @param sNN 
   * @param trg 
   */
  void Run(const char* output="trending.root", 
	   UShort_t sys=2, 
	   UShort_t sNN=2760, 
	   UShort_t trg=1)
  {
    TString outName(output);
    outName.ReplaceAll(".pdf", ".root");
    if (!outName.EndsWith(".root")) outName.Append(".root");
    TString pdfName(outName);
    pdfName.ReplaceAll(".root", ".pdf");
    
    CreateCanvas(pdfName, true);
    
    TFile*      out    = TFile::Open(outName, "RECREATE");
    TDirectory* refDir = out->mkdir("reference");

    TIter iCent(&fCents);
    TObject* pCent = 0;
    while ((pCent = iCent())) {
      UShort_t low   = (pCent->GetUniqueID() & 0xFFFF);
      UShort_t high  = (pCent->GetUniqueID() >> 16) & 0xFFFF;
      TGraph*  graph = GetOther(sys, sNN, trg, low, high);
      if (!graph) break;
      
      TString name(graph->GetName());
      name.Append(Form("_%s", pCent->GetName()));
      graph->SetName(name);

      fOthers.Add(graph);
      refDir->Add(graph);
    }

    TIter    iRun(&fRuns);
    TObject* pRun = 0;
    while ((pRun = iRun())) {
      TString     sRun = pRun->GetName();
      TDirectory* dRun = out->mkdir(sRun);
      
      MakeChapter(sRun);
      
      for (Int_t i = 0; i < 1/*2*/; i++) { 
	TString     three(Form("%dstrip", i+2));
	TDirectory* dStrip = dRun->mkdir(three);

	MakeChapter(Form("_%s", three.Data()));

	three.Append("_slX_shX_dcX");
	
	LoopCuts(sRun, three, fOrder[0], fOrder[1], fOrder[2], dStrip);
      }
    }
    out->Write();
    out->Close();
    CloseCanvas();
  }
  //__________________________________________________________________
  /** 
   * Loop over cuts.  Recursively calls self. 
   * 
   * @param run      Run number
   * @param cur      Current path
   * @param cuts1    Cuts
   * @param cuts2    Cuts (or null)
   * @param cuts3    Cuts (or null)
   * @param out      Output diretory
   * @param upBin    Parents bin number 
   * @param upMean   Parents mean
   * @param upRatios Parents ratios
   */    
  void LoopCuts(const TString& run, 
		const TString& cur, 
		const TList*   cuts1,
		const TList*   cuts2,
		const TList*   cuts3,
		TDirectory*    out,
		Int_t          upBin=0,
		THStack*       upMean=0, 
		THStack*       upRatios=0)
  {
    TString pre(Form("%s%s",(cuts3 ? "" : "_"), (cuts2 ? "" : "_")));
    Printf("%s%s", pre.Data(), cur.Data());

    TDirectory* dCut = out->mkdir(cuts1->GetName()); 
    TIter       iCut(cuts1);
    TObject*    pCut = 0;
    dCut->cd();
    while ((pCut = iCut())) { 
      TString     sMethod = pCut->GetName();
      TDirectory* dMethod = dCut->mkdir(sMethod);
      TString     sValues = pCut->GetTitle();
      TObjArray*  aValues = sValues.Tokenize(" ");
      TIter       iValue(aValues);
      TObjString* pValue = 0;
      TString     templ(Form("%s%s", cuts1->GetName(), sMethod.Data()));
      TString     base(cur);
      base.ReplaceAll(Form("%sX", cuts1->GetName()), templ);
      aValues->SetName(cuts1->GetName());

      THStack*     mean   = 0;
      THStack*     ratios = 0;
      dMethod->cd();
      MakeStacks(base, aValues, mean, ratios);
      dMethod->Add(mean);
      dMethod->Add(ratios);

      Int_t xbin = 1;
      while ((pValue = static_cast<TObjString*>(iValue()))) {
	TString     sValue(Form("%05.2f", pValue->String().Atof()));
	sValue.ReplaceAll(".", "d");
	TDirectory* dValue = dMethod->mkdir(sValue);
	TString now(base);
	TString sub(Form("%s%s", templ.Data(), sValue.Data()));
	now.ReplaceAll(templ.Data(), sub.Data());
	// now.Append(sValue);
	dValue->cd();

	if (cuts2) {
	  // Loop over sub-cut.
	  LoopCuts(run, now, cuts2, cuts3, 0, 
		   dValue, xbin, mean, ratios);
	}
	else {
	  // Process for a given cut with the current value of that cut. 
	  if (!NextFile(run, now, dValue, xbin, 
			mean, ratios)) {
	    // xbin++;
	    dMethod->cd();
	    continue;
	  }
	  // After this, we have points for the current cut value
	  // stored in the stacks mean and wSpread, and added ratios
	  // for the current cut to the stack ratios.  Since we're not
	  // done with the cut values just yet, we should wait to
	  // update the parent stacks.
	}

	if (upMean /*&& upWSpread*/) {
	  // For each centrality extract 
	  // - mean of mean of ... 
	  // - var of var of ... 
	  // - spread of spread of ... 
	  TIter    iCent(&fCents);
	  TObject* pCent = 0;
	  Int_t    jCent = 0;
	  while ((pCent = iCent())) {
	    TH1* h = static_cast<TH1*>(mean->GetHists()->At(jCent));
	    // Info("", "Updating parent %s with %s @ %d", 
	    //      upMean->GetTitle(), h->GetTitle(),upBin);
	    UpdateStacks(h, jCent, upBin, upMean);
	    jCent++;
	  }
	} // if ups
	
	xbin++;
	dMethod->cd();
      } // for values 
      
      if (upRatios) {
	TIter nextHist(ratios->GetHists());
	TH1*  hist = 0;
	while ((hist = static_cast<TH1*>(nextHist()))) 
	  upRatios->Add(hist);
      }      
      dCut->cd();

      FixMinMax(ratios);
      FixMinMax(mean);

      if (mean && mean->GetHists()->GetEntries() > 0 && 
	  ratios && ratios->GetHists() && 
	  ratios->GetHists()->GetEntries() > 0) {
	fBody->Divide(2,1);
	DrawStacks(fBody->GetPad(1), mean, kNorth|kWest);
	DrawStacks(fBody->GetPad(2), ratios, kNorth|kCenter, kSouth|kCenter);
	PrintCanvas(Form("%s_%s", pre.Data(), base.Data()), 0.5);
      }
    } // for methods 
    out->cd();
  }
  //__________________________________________________________________
  /** 
   * Process the next file 
   * 
   * @param run        Run 
   * @param now        Time
   * @param out        Output directory 
   * @param binx       Bin number 
   * @param mean       Graphs of mean +/- variance 
   * @param ratios     Ratios stack
   * 
   * @return true on success 
   */
  Bool_t NextFile(const TString& run, const TString& now,
		  TDirectory* out, Int_t binx, 
		  THStack* mean, THStack* ratios)
  {
    TString dir(Form("%s_dndeta_%s", run.Data(), now.Data()));
    TString path(Form("%s/forward_dndeta.root", dir.Data()));
    if (gSystem->AccessPathName(path.Data())) {
      Warning("NextFile", "%s not found", path.Data());
      return false;
    }
    TFile* file = TFile::Open(path, "READ");
    if (!file) { 
      Warning("NextFile", "Failed to open %s", path.Data());
      return false;
    }
    // Info("NextFile", "Opened %s", path.Data());
    
    TCollection* results = GetCollection(file, "ForwarddNdetaResults");
    if (!results) return false;

    TCollection* mcResults = GetCollection(file, "MCTruthdNdetaResults");

    THStack* all = new THStack("all",    "All");
    THStack* rat = new THStack("ratios", "Ratios");
    
    TIter    iCent(&fCents);
    TObject* pCent = 0;
    Int_t    jCent = 0;
    while ((pCent = iCent())) {
      TString folderName(pCent->GetName());
      TCollection* centFolder = GetCollection(results, folderName);
      if (!centFolder) {
	Warning("", "Didn't get the centrality %s folder from %s",
		folderName.Data(), results->GetName());
	// results->ls();
	break;
      }

      TCollection* mcCentFolder = 0;
      if (mcResults) mcCentFolder = GetCollection(mcResults, folderName);
      
      TH1* dNdeta = GetH1(centFolder, Form("dndetaForward%s", 
					   fRebinned ? "_rebin05" : ""));
      if (!dNdeta) {
	Warning("", "Didn't get histogram for jCent=%d in %s", 
		jCent, path.Data());
	// results->ls();
	break;
      }
      dNdeta->SetDirectory(out);
      dNdeta->SetName(folderName);
      dNdeta->SetMarkerColor(jCent+1);
      dNdeta->SetLineColor(jCent+1);

      TH1* other = 0; 
      if (mcCentFolder) 
	other = GetH1(mcCentFolder, Form("dndetaMCTruth%s", 
					 fRebinned ? "_rebin05" : ""));
      if (!other) {
	TGraph*  graph = static_cast<TGraph*>(fOthers.At(jCent));
	if (!graph) break;

	other = G2H(graph, *(dNdeta->GetXaxis()));
      }

      if (!other) {
	Warning("", "No other data found for %s", path.Data());
	break;
      }

      other->SetMarkerColor(dNdeta->GetMarkerColor());
      other->SetMarkerSize(dNdeta->GetMarkerSize());
      other->SetLineColor(dNdeta->GetLineColor());
      other->SetDirectory(out);
      folderName.ReplaceAll("cent", "other");
      folderName.ReplaceAll("all",  "other");
      other->SetName(folderName);

      all->Add(other);
      all->Add(dNdeta);
      
      folderName.ReplaceAll("other", "ratio");
      TH1* ratio = static_cast<TH1*>(dNdeta->Clone(folderName));
      ratio->SetTitle(Form("%s %s", dir.Data(), pCent->GetName()));
      ratio->SetDirectory(out);
      ratio->Divide(other);
      ratio->SetMarkerStyle(20+binx-1);
      ratio->SetMarkerColor(jCent+1);
      ratio->SetLineColor(kGray);
      ratio->SetLineWidth(0);
      ratios->Add(ratio);

      rat->Add(ratio);
      
      UpdateStacks(ratio, jCent, binx, mean);
      jCent++;
    }
    out->Add(all);
    out->Add(rat);

    FixMinMax(all);
    FixMinMax(rat);
    
    fBody->Divide(1,2,0,0);
    DrawStacks(fBody->GetPad(1), all/*, 21*/);
    DrawStacks(fBody->GetPad(2), rat, kNorth|kCenter);
    PrintCanvas(Form("____%s", now.Data()), .4);

    file->Close();
    if (jCent <= 0) return false;

    return true;
  }
  //=== Stack functions ==============================================
  /** 
   * Make stacks 
   * 
   * @param run      Run number
   * @param values   Cut parameters
   * @param mean     Stack of means
   * @param ratios   Stack of ratios 
   */
  void MakeStacks(const TString&   run,
		  const TObjArray* values, 
		  THStack*&        mean, 
		  THStack*&        ratios)
  {
    mean    = new THStack("mean", run);
    ratios  = new THStack("ratios", run);

    // --- Create histograms and graphs ------------------------
    Int_t    nValues = values->GetEntriesFast();
    TIter    nextCent(&fCents);
    TObject* pcent = 0;
    Int_t    col   = 1;
    while ((pcent = nextCent())) {
      UShort_t low   = (pcent->GetUniqueID() & 0xFFFF);
      UShort_t high  = (pcent->GetUniqueID() >> 16) & 0xFFFF;
      
      TH1* hMean = new TH1D(pcent->GetName(), 
			    Form("%s %d%%-%d%% central", run.Data(), 
				 low, high),
			    nValues, .5, nValues+.5);
      hMean->SetMarkerColor(col);
      hMean->SetMarkerStyle(20);
      hMean->SetLineColor(col);
      hMean->SetXTitle(Form("Cut parameter X_{%s}", values->GetName()));
      hMean->SetYTitle("Average, spread, and min/max of ratio");
      hMean->SetDirectory(0);
      TIter nextV(values);
      TObjString* pvalue = 0;
      Int_t xbin = 1;
      while ((pvalue = static_cast<TObjString*>(nextV()))) {
	TString& value = pvalue->String();
	hMean->GetXaxis()->SetBinLabel(xbin, value);
	xbin++;
      }

      TGraphAsymmErrors* gWSpread = new TGraphAsymmErrors(nValues);
      gWSpread->SetName("minmax");
      gWSpread->SetTitle(hMean->GetTitle());
      gWSpread->SetMarkerColor(col);
      gWSpread->SetLineColor(col);
      gWSpread->SetFillColor(0);
      gWSpread->SetFillStyle(0);

      hMean->GetListOfFunctions()->Add(gWSpread, "[]pl same");
      mean->Add(hMean, "x0 e1");
      col++;
    }
  }
  //__________________________________________________________________
  /** 
   * Update stacks 
   * 
   * @param h       Histogram
   * @param i       Index 
   * @param binx    Bin number 
   * @param stack   Stack to update
   */
  void UpdateStacks(TH1*         h, 
		    Int_t        i, 
		    Int_t        binx,
		    THStack*     stack)
  {
    Double_t avg = 0;
    Double_t var = 0;
    Double_t min = +100000;
    Double_t max = -100000;
    HistStatistics(h, avg, var, min, max);
    h->SetMinimum(0.95*min);
    h->SetMaximum(1.1*max);
    
    TH1* hMean = static_cast<TH1*>(stack->GetHists()->At(i));
    hMean->SetBinContent(binx, avg);
    hMean->SetBinError(binx,var);

    TObject* pG = hMean->GetListOfFunctions()->FindObject("minmax");
    if (!pG) {
      hMean->GetListOfFunctions()->ls();
      return;
    }
    TGraphAsymmErrors* gWSpread = static_cast<TGraphAsymmErrors*>(pG);
    gWSpread->SetPoint(binx-1, binx, avg);
    gWSpread->SetPointError(binx-1,0, 0, /*0.5,0.5,*/
			    TMath::Abs(avg-min), TMath::Abs(max-avg));
  }
  //=== Graphics member functions ====================================
  //__________________________________________________________________
  /** 
   * Build a centrality legend
   * 
   * @param p     Pad
   * @param where See above 
   * @param stack Stack to make legend for 
   */
  void BuildCentLegend(TVirtualPad* p, UInt_t where, THStack* stack=0)
  {
    TLegend* l     = MakeLegend(p, where, false);
    TIter    iCent(&fCents);
    TObject* pCent = 0;
    Int_t    col   = 1;
    while ((pCent = iCent())) {
      UShort_t low   = (pCent->GetUniqueID() & 0xFFFF);
      UShort_t high  = (pCent->GetUniqueID() >> 16) & 0xFFFF;

      TLegendEntry* e = l->AddEntry("dummy", 
				    Form("%2d%% - %2d%% central", low, high),
				    "pl");
      e->SetFillColor(col);
      e->SetMarkerColor(col);
      e->SetLineColor(col);
      e->SetMarkerStyle(20);
      col++;
    }
    if (stack && stack->GetHistogram()) { 
      stack->GetHistogram()->GetListOfFunctions()->Add(l);
    }
    else 
      l->Draw();
    
  }
  //__________________________________________________________________
  /** 
   * Make a cut legend. 
   * 
   * @param p      Pad
   * @param where  See above 
   * @param stack  Stack to make legend for
   */
  void BuildCutLegend(TVirtualPad* p, UInt_t where, THStack* stack)
  {
    TLegend* l     = MakeLegend(p, where, false);
    l->SetX1(p->GetLeftMargin());
    l->SetX2(1-p->GetRightMargin());
    // l->SetBorderSize(1);
    l->SetNColumns(1);

    TList    seen;
    TIter    iHist(stack->GetHists());
    TH1*     pHist = 0;
    while ((pHist = static_cast<TH1*>(iHist()))) {
      TString title(pHist->GetTitle());

      Int_t   idx = title.Index(" cent");
      if (idx != kNPOS) title.Remove(idx, 5+3+3+1);

      idx = title.Index("_dndeta");
      if (idx != kNPOS) title.Remove(0, idx+6+1+1);
      
      TObject* before = seen.FindObject(title);
      if (before) continue;

      seen.Add(new TObjString(title));
      
      TLegendEntry* e = l->AddEntry("dummy", title, "p");
      e->SetMarkerColor(kBlack);
      e->SetMarkerStyle(pHist->GetMarkerStyle());
    }
    seen.IsOwner();
    if (stack->GetHistogram()) { 
      stack->GetHistogram()->GetListOfFunctions()->Add(l);
    }
    else 
      l->Draw();
    
  }
  //__________________________________________________________________
  /** 
   * Draw stacks
   *  
   * @param p      Pad
   * @param stack  Stacks
   * @param cent   Centrality legend location
   * @param cuts   Cut legend location
   */
  void DrawStacks(TVirtualPad* p, 
		  THStack* stack, 
		  UInt_t cent=0, 
		  UInt_t cuts=0)
  {
    if (!stack) {
      Warning("DrawStacks", "Stack is missing!");
      return;
    }
    if (!stack->GetHists() || stack->GetHists()->GetEntries() <= 0) { 
      Warning("DrawStacks", "Stack is empty");
      return;
    }
    TH1*    first = static_cast<TH1*>(stack->GetHists()->At(0));
    TString xT    = first->GetXaxis()->GetTitle();
    
    p->cd();
    stack->Draw("nostack");
    // DrawInPad(p, 0, stack, "nostack", flags);
    stack->GetXaxis()->SetTitle(xT);
    FixMinMax(stack);

    if (cent > 0) BuildCentLegend(p, cent, stack);
    if (cuts > 0) BuildCutLegend(p, cuts, stack);

    p->Modified();
    p->Update();
    p->cd();
  }
  //=== Utility static member functions ==============================
  /** 
   * Update statistics
   * 
   * @param y    Current observation
   * @param cnt  On entry, last count, on return current count
   * @param mean On entry, last mean, on return current mean
   * @param var  On entry, last variance, on return current variance 
   */
  static void Statistics(Double_t  y,
			 Int_t&    cnt,
			 Double_t& mean, 
			 Double_t& var) 
  {
    cnt        += 1;
    mean       += (y - mean) / cnt;
    var        += (cnt > 1 ? (TMath::Power(y-mean,2)/(cnt-1)-var/cnt) : 0);
  }
  //__________________________________________________________________
  /** 
   * Calculate histogram statistics
   *  
   * @param h     Histogram
   * @param mean  On return, the Y-mean
   * @param var   On return, the Y variance 
   * @param min   On return, the least Y 
   * @param max   On return, the largest Y 
   */
  static void HistStatistics(const TH1* h, 
			     Double_t& mean, 
			     Double_t& var, 
			     Double_t& min, 
			     Double_t& max)
  {
    mean      = 0;
    var       = 0;
    min       = +100000;
    max       = -100000;
    Int_t cnt = 0;
    for (Int_t i = 1; i <= h->GetNbinsX(); i++) { 
      Double_t y = h->GetBinContent(i);
      if (TMath::Abs(y) <= 1e-9) continue;
      min        =  TMath::Min(min, y);
      max        =  TMath::Max(max, y);
      Statistics(y, cnt, mean, var);
    }
    // Info("", "Stats for %s:  mean=%f +/- %f [%f,%f]",
    //      h->GetTitle(), mean, var, min, max);
  }
  //__________________________________________________________________
  /** 
   * Fix the min/max of a stack
   * 
   * @param stack 
   */
  static void FixMinMax(THStack* stack)
  {
    TIter iHist(stack->GetHists());
    TH1*  pHist = 0;
    Double_t m1 = 10000000;
    Double_t m2 = -10000000;
    while((pHist = static_cast<TH1*>(iHist()))) { 
      m1 = TMath::Min(m1, pHist->GetMinimum(1e-6));
      m2 = TMath::Max(m2, pHist->GetMaximum());
    }
    // Double_t m1 = stack->GetMinimum("nostack e");
    // Double_t m2 = stack->GetMaximum("nostack e");
    // Printf("Stack %s minimum: %f", stack->GetTitle(), m1);
    stack->SetMinimum((m1 < 0 ? 1.05 : 0.95)  * m1);
    stack->SetMaximum((m2 > 0 ? 1.05 : 0.95)  * m2);
  }
  //__________________________________________________________________
  /** 
   * Turn a graph into a histogram
   * 
   * @param g    Graph
   * @param axis Axis to use 
   * 
   * @return Newly allocated histogram 
   */
  static TH1* G2H(const TGraph* g, const TAxis& axis)
  {
    TH1* h = 0;
    if (axis.GetXbins()->GetArray()) 
      h = new TH1D(g->GetName(), g->GetTitle(), 
		   axis.GetNbins(), axis.GetXbins()->GetArray());
    else 
      h = new TH1D(g->GetName(), g->GetTitle(), 
		   axis.GetNbins(), axis.GetXmin(), axis.GetXmax());
    h->SetMarkerColor(g->GetMarkerColor());
    h->SetMarkerStyle(g->GetMarkerStyle());
    h->SetMarkerSize(g->GetMarkerSize());
    h->SetLineColor(g->GetLineColor());
    h->SetLineStyle(g->GetLineStyle());
    h->SetLineWidth(g->GetLineWidth());
    h->SetFillColor(g->GetFillColor());
    h->SetFillStyle(g->GetFillStyle());
    h->SetDirectory(0);
    
    for (Int_t i = 1; i <= h->GetNbinsX(); i++) { 
      Double_t x = h->GetXaxis()->GetBinCenter(i);
      Double_t y = g->Eval(x); // , 0, "S");
      h->SetBinContent(i, y);
    }
    return h;
  }
  //__________________________________________________________________
  /** 
   * Get other data
   * 
   * @param sys    Collision system
   * @param sNN    Collision energy 
   * @param trg    Trigger 
   * @param lowC   Low centrality 
   * @param highC  High centrality 
   * 
   * @return Newly allocated histogram
   */
  static TGraph* GetOther(UShort_t      sys, 
			  UShort_t      sNN, 
			  UShort_t      trg,
			  UShort_t      lowC, 
			  UShort_t      highC)
  {
    return 0;
#if 0
    // --- Set the macro pathand load other data script --------------
    // Always recompile 
    if (!gROOT->GetClass("RefData")) {
      TString savPath(gROOT->GetMacroPath());
      gROOT->SetMacroPath(Form("%s:$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2",
			       gROOT->GetMacroPath()));
      gROOT->LoadMacro("OtherData.C++");
      gROOT->SetMacroPath(savPath);
    }
    Long_t ret = gROOT->ProcessLine(Form("RefData::GetData(%d,%d,%d,%d,%d,0x4)",
					 sys, sNN, trg, lowC, highC));
    if (!ret) { 
      Warning("GetOther", "No other data for %d %d %d %d%%-%d%% central", 
	      sys, sNN, trg, lowC, highC);
      return 0;
    }
    TMultiGraph* others = reinterpret_cast<TMultiGraph*>(ret);    
    TGraph*      other  =static_cast<TGraph*>(others->GetListOfGraphs()->At(0));
    if (!other) {
      Warning("GetOther", "No ALICE data for %d %d %d %d%%-%d%% central", 
	      sys, sNN, trg, lowC, highC);
      return 0;
    }
    
    // Info("", "Got graph %s/%s", other->GetName(), other->GetTitle());

    return other;
#endif 
  }

  //__________________________________________________________________
  TList     fSLCuts;
  TList     fSHCuts;
  TList     fDCCuts;
  TList     fRuns;
  TList     fCents;
  TList     fOthers;
  TList*    fOrder[3];	
  Bool_t    fRebinned;
};

#endif
