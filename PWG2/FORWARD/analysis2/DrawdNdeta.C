/**
 * @file   DrawdNdeta.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:07:10 2011
 * 
 * @brief  Script to visualise the dN/deta 
 *
 * This script is independent of any AliROOT code - and should stay
 * that way.
 * 
 * 
 * @ingroup pwg2_forward_dndeta
 */
#include <TH1.h>
#include <THStack.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TList.h>
#include <TString.h>
#include <TError.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TImage.h>

/**
 * Class to draw dN/deta results 
 * 
 * @ingroup pwg2_forward_tasks_dndeta
 * @ingroup pwg2_forward_dndeta
 */
struct dNdetaDrawer 
{
  /**
   * POD of data for range zooming 
   */
  struct RangeParam 
  {
    TAxis*       fMasterAxis; // Master axis 
    TAxis*       fSlave1Axis; // First slave axis 
    TVirtualPad* fSlave1Pad;  // First slave pad 
    TAxis*       fSlave2Axis; // Second slave axis 
    TVirtualPad* fSlave2Pad;  // Second slave pad 
  };
  //__________________________________________________________________
  /** 
   * Constructor 
   * 
   */
  dNdetaDrawer()
    : fShowOthers(false),	// Bool_t 
      fShowAlice(false),        // Bool_t
      fShowRatios(false),	// Bool_t 
      fShowLeftRight(false),	// Bool_t
      fRebin(5),		// UShort_t
      fCutEdges(false),		// Bool_t
      fTitle(""),		// TString
      fTrigString(0),		// TNamed*
      fNormString(0),           // TNamed*
      fSNNString(0),		// TNamed*
      fSysString(0),		// TNamed*
      fVtxAxis(0),		// TAxis*
      fCentAxis(0),             // TAxis*
      fTriggers(0),		// TH1*
      fRangeParam(0)

  {
    fRangeParam = new RangeParam;
    fRangeParam->fMasterAxis = 0;
    fRangeParam->fSlave1Axis = 0;
    fRangeParam->fSlave1Pad  = 0;
    fRangeParam->fSlave2Axis = 0;
    fRangeParam->fSlave2Pad  = 0;
  }
  //__________________________________________________________________
  virtual ~dNdetaDrawer()
  {
    if (fRatios  && fRatios->GetHists())  fRatios->GetHists()->Delete();
    if (fResults && fResults->GetHists()) fResults->GetHists()->Delete();

    if (fTrigString) { delete fTrigString; fTrigString = 0; }
    if (fSNNString)  { delete fSNNString;  fSNNString  = 0; }
    if (fSysString)  { delete fSysString;  fSysString  = 0; }
    if (fVtxAxis)    { delete fVtxAxis;    fVtxAxis    = 0; }
    if (fCentAxis)   { delete fCentAxis;   fCentAxis   = 0; }
    if (fResults)    { delete fResults;    fResults    = 0; }
    if (fRatios)     { delete fRatios;     fRatios     = 0; }
    if (fOthers)     { delete fOthers;     fOthers     = 0; }
    if (fTriggers)   { delete fTriggers;   fTriggers   = 0; } 
    fRangeParam = 0;
  }

  //==================================================================
  /** 
   * @{ 
   * @name Set parameters 
   */
  /** 
   * Show other (UA5, CMS, ...) data 
   * 
   * @param x Whether to show or not 
   */
  void SetShowOthers(Bool_t x)    { fShowOthers = x; }
  //__________________________________________________________________
  /** 
   * Show ALICE published data 
   * 
   * @param x Wheter to show or not 
   */
  void SetShowAlice(Bool_t x)     { fShowAlice = x; }
  //__________________________________________________________________
  /** 
   * Whether to show ratios or not.  If there's nothing to compare to,
   * the ratio panel will be implicitly disabled
   * 
   * @param x Whether to show or not 
   */
  void SetShowRatios(Bool_t x)    { fShowRatios = x; }
  //__________________________________________________________________
  /** 
   * 
   * Whether to show the left/right asymmetry 
   *
   * @param x To show or not 
   */
  void SetShowLeftRight(Bool_t x) { fShowLeftRight = x; }
  //__________________________________________________________________
  /** 
   * Set the rebinning factor 
   * 
   * @param x Rebinning factor (must be a divisor in the number of bins) 
   */
  void SetRebin(UShort_t x)       { fRebin = x; }
  //__________________________________________________________________
  /** 
   * Wheter to cut away the edges 
   * 
   * @param x Whether or not to cut away edges 
   */
  void SetCutEdges(Bool_t x)      { fCutEdges = x; }
  //__________________________________________________________________
  /** 
   * Set the title of the plot
   * 
   * @param x Title
   */
  void SetTitle(TString x)        { fTitle = x; }
  /* @} */
  //==================================================================  
  /** 
   * @{ 
   * @name Override settings from input 
   */
  /** 
   * Override setting from file 
   * 
   * @param sNN Center of mass energy per nucleon pair (GeV)
   */
  void SetSNN(UShort_t sNN) 
  {
    fSNNString = new TNamed("sNN", Form("%04dGeV", sNN));
    fSNNString->SetUniqueID(sNN);
  }
  //__________________________________________________________________
  /** 
   * Set the collision system 
   * - 1: pp 
   * - 2: PbPb
   * 
   * @param sys collision system
   */
  void SetSys(UShort_t sys)
  {
    fSysString = new TNamed("sys", (sys == 1 ? "pp" : 
				    sys == 2 ? "PbPb" : "unknown"));
    fSysString->SetUniqueID(sys);
  }
  //__________________________________________________________________
  /** 
   * Set the vertex range in centimeters 
   * 
   * @param vzMin Min @f$ v_z@f$
   * @param vzMax Max @f$ v_z@f$
   */
  void SetVertexRange(Double_t vzMin, Double_t vzMax) 
  {
    fVtxAxis = new TAxis(10, vzMin, vzMax);
    fVtxAxis->SetName("vtxAxis");
    fVtxAxis->SetTitle(Form("v_{z}#in[%+5.1f,%+5.1f]cm", vzMin, vzMax));
  }
  //__________________________________________________________________
  void SetTrigger(UShort_t trig)
  {
    fTrigString = new TNamed("trigString", (trig & 0x1 ? "INEL" : 
					    trig & 0x2 ? "INEL>0" : 
					    trig & 0x4 ? "NSD" : 
					    "unknown"));
    fTrigString->SetUniqueID(trig);
  }


  //==================================================================
  /** 
   * @{ 
   * @name Main steering functions 
   */
  /** 
   * Run the code to produce the final result. 
   * 
   * @param filename  File containing the data 
   */
  void Run(const char* filename="forward_dndeta.root") 
  {
    Double_t max = 0, rmax=0, amax=0;

    gStyle->SetPalette(1);

    // --- Open input file -------------------------------------------
    TFile* file = TFile::Open(filename, "READ");
    if (!file) { 
      Error("Open", "Cannot open %s", filename);
      return;
    }
    // --- Get forward list ------------------------------------------
    TList* forward = static_cast<TList*>(file->Get("ForwardResults"));
    if (!forward) { 
      Error("Open", "Couldn't find list ForwardResults");
      return;
    }
    // --- Get information on the run --------------------------------
    FetchInformation(forward);
    // --- Set the macro pathand load other data script --------------
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWG2/FORWARD/analysis2",
			     gROOT->GetMacroPath()));
    gROOT->LoadMacro("OtherData.C");

    // --- Get the central results -----------------------------------
    TList* clusters = static_cast<TList*>(file->Get("CentralResults"));
    if (!clusters) Warning("Open", "Couldn't find list CentralResults");

    // --- Make our containtes ---------------------------------------
    fResults   = new THStack("results", "Results");
    fRatios    = new THStack("ratios",  "Ratios");
    fLeftRight = new THStack("asymmetry", "Left-right asymmetry");
    fOthers    = new TMultiGraph();
    
    // --- Loop over input data --------------------------------------
    FetchResults(forward,  "Forward", max, rmax, amax);
    FetchResults(clusters, "Central", max, rmax, amax);

    // --- Get trigger information -----------------------------------
    TList* sums = static_cast<TList*>(file->Get("ForwardSums"));
    if (sums) {
      TList* all = static_cast<TList*>(sums->FindObject("all"));
      if (all) {
	fTriggers = FetchResult(all, "triggers");
	if (!fTriggers) all->ls();
      }
      else  {
	Warning("Run", "List all not found in ForwardSums");
	sums->ls();
      }
    }
    else { 
      Warning("Run", "No ForwardSums directory found in %s", file->GetName());
      file->ls();
    }
    
    // --- Check our stacks ------------------------------------------
    if (!fResults->GetHists() || 
	fResults->GetHists()->GetEntries() <= 0) { 
      Error("Run", "No histograms in result stack!");
      return;
    }
    if (!fOthers->GetListOfGraphs() || 
	fOthers->GetListOfGraphs()->GetEntries() <= 0) { 
      Warning("Run", "No other data found - disabling that");
      fShowOthers = false;
    }
    if (!fRatios->GetHists() || 
	fRatios->GetHists()->GetEntries() <= 0) { 
      Warning("Run", "No ratio data found - disabling that");
      // fRatios->ls();
      fShowRatios = false;
    }
    if (!fLeftRight->GetHists() || 
	fLeftRight->GetHists()->GetEntries() <= 0) { 
      Warning("Run", "No left/right data found - disabling that");
      // fLeftRight->ls();
      fShowLeftRight = false;
    }
    
    // --- Close the input file --------------------------------------
    file->Close();

    

    // --- Plot results ----------------------------------------------
    Plot(max, rmax, amax);
  }

  //__________________________________________________________________
  /** 
   * Fetch the information on the run from the results list
   * 
   * @param results  Results list
   */
  void FetchInformation(const TList* results)
  {
    if (!fTrigString) 
      fTrigString = static_cast<TNamed*>(results->FindObject("trigger"));
    if (!fNormString) 
      fNormString = static_cast<TNamed*>(results->FindObject("scheme"));
    if (!fSNNString) 
      fSNNString  = static_cast<TNamed*>(results->FindObject("sNN"));
    if (!fSysString) 
      fSysString  = static_cast<TNamed*>(results->FindObject("sys"));
    if (!fVtxAxis)
      fVtxAxis    = static_cast<TAxis*>(results->FindObject("vtxAxis"));
    if (!fCentAxis) 
      fCentAxis   = static_cast<TAxis*>(results->FindObject("centAxis"));

    TNamed* options = static_cast<TAxis*>(results->FindObject("options"));
    if (!fTrigString) fTrigString = new TNamed("trigger", "unknown");
    if (!fNormString) fNormString = new TNamed("scheme", "unknown");
    if (!fSNNString)  fSNNString  = new TNamed("sNN", "unknown");
    if (!fSysString)  fSysString  = new TNamed("sys", "unknown");
    if (!fVtxAxis) { 
      fVtxAxis    = new TAxis(1,0,0);
      fVtxAxis->SetName("vtxAxis");
      fVtxAxis->SetTitle("v_{z} range unspecified");
    }

    TString centTxt("none");
    if (fCentAxis) { 
      Int_t nCent = fCentAxis->GetNbins();
      centTxt = Form("%d bins", nCent);
      for (Int_t i = 0; i <= nCent; i++) 
	centTxt.Append(Form("%c%d", i == 0 ? ' ' : '-', 
			    int(fCentAxis->GetXbins()->At(i))));
    }
    Info("FetchInformation", 
	 "Initialized for\n"
	 "   Trigger:       %-30s  (%d)\n"
	 "   sqrt(sNN):     %-30s  (%dGeV)\n"
	 "   System:        %-30s  (%d)\n"
	 "   Vz range:      %-30s  (%f,%f)\n"
	 "   Normalization: %-30s  (%d)\n"
	 "   Centrality:    %s\n"
	 "   Options:       %s",
	 fTrigString->GetTitle(), fTrigString->GetUniqueID(), 
	 fSNNString->GetTitle(),  fSNNString->GetUniqueID(), 
	 fSysString->GetTitle(),  fSysString->GetUniqueID(), 
	 fVtxAxis->GetTitle(), fVtxAxis->GetXmin(), fVtxAxis->GetXmax(),
	 fNormString->GetTitle(), fNormString->GetUniqueID(),
	 centTxt.Data(), (options ? options->GetTitle() : "none"));
  }
  //__________________________________________________________________
  TMultiGraph* FetchOthers(UShort_t centLow, UShort_t centHigh)
  {
    TMultiGraph* thisOther = 0;
    if (!fShowOthers && !fShowAlice) return 0;

    Bool_t   onlya = (fShowOthers ? false : true);
    UShort_t sys   = (fSysString  ? fSysString->GetUniqueID() : 0);
    UShort_t trg   = (fTrigString ? fTrigString->GetUniqueID() : 0);
    UShort_t snn   = (fSNNString  ? fSNNString->GetUniqueID() : 0);
    Long_t   ret   = gROOT->ProcessLine(Form("GetData(%d,%d,%d,%d,%d,%d);",
					     sys,snn,trg,
					     centLow,centHigh,onlya));
    if (!ret) { 
      Warning("FetchResults", "No other data found for sys=%d, sNN=%d, "
	      "trigger=%d %d%%-%d%% central %s", 
	      sys, snn, trg, centLow, centHigh, 
	      onlya ? " (ALICE results)" : "all");
      return 0;
    }
    thisOther = reinterpret_cast<TMultiGraph*>(ret);
    
    return thisOther;
  }
  //__________________________________________________________________
  /** 
   * Get the results from the top-level list 
   * 
   * @param list  List 
   * @param name  name 
   * @param max   On return, maximum of data 
   * @param rmax  On return, maximum of ratios
   * @param amax  On return, maximum of left-right comparisons
   */
  void FetchResults(const TList* list, 
		    const char*  name, 
		    Double_t&    max,
		    Double_t&    rmax,
		    Double_t&    amax)
  {
    UShort_t n = fCentAxis ? fCentAxis->GetNbins() : 0;
    if (n == 0) {
      TList* all = static_cast<TList*>(list->FindObject("all"));
      if (!all)
	Error("FetchResults", "Couldn't find list 'all' in %s", 
	      list->GetName());
      else 
	FetchResults(all, name, FetchOthers(0,0), -1, 0, max, rmax, amax);
      return;
    }
    
    Int_t   nCol = gStyle->GetNumberOfColors();
    for (UShort_t i = 0; i < n; i++) { 
      UShort_t centLow  = fCentAxis->GetXbins()->At(i);
      UShort_t centHigh = fCentAxis->GetXbins()->At(i+1);
      TString  lname    = Form("cent%03d_%03d", centLow, centHigh);
      TList*   thisCent = static_cast<TList*>(list->FindObject(lname));

      Float_t fc   = (centLow+double(centHigh-centLow)/2) / 100;
      Int_t   icol = TMath::Min(nCol-1,int(fc * nCol + .5));
      Int_t   col  = gStyle->GetColorPalette(icol);
      Info("FetchResults","Centrality %d-%d color index %d (=%d*%f) -> %d", 
	   centLow, centHigh, icol, nCol, fc, col);

      TString centTxt = Form("%3d%%-%3d%% central", centLow, centHigh);
      if (!thisCent)
	Error("FetchResults", "Couldn't find list '%s' in %s", 
	      lname.Data(), list->GetName());
      else 
	FetchResults(thisCent, name, FetchOthers(centLow, centHigh), 
		     col, centTxt.Data(), max, rmax, amax);
    }
  } 
  //__________________________________________________________________
  void SetAttributes(TH1* h, Int_t color)
  {
    if (!h) return;
    if (color < 0) return;
    // h->SetLineColor(color);
    h->SetMarkerColor(color);
    // h->SetFillColor(color);
  }
  //__________________________________________________________________
  void SetAttributes(TGraph* g, Int_t color)
  {
    if (!g) return;
    if (color < 0) return;
    // g->SetLineColor(color);
    g->SetMarkerColor(color);
    // g->SetFillColor(color);
  }
  //__________________________________________________________________
  void ModifyTitle(TNamed* h, const char* centTxt)
  {
    if (!centTxt || !h) return;
    h->SetTitle(Form("%s, %s", h->GetTitle(), centTxt));
  }

  //__________________________________________________________________
  /** 
   * Fetch results for a particular centrality bin
   * 
   * @param list       List 
   * @param name       Name 
   * @param thisOther  Other graphs 
   * @param color      Color 
   * @param centTxt    Centrality text
   * @param max        On return, data maximum
   * @param rmax       On return, ratio maximum 
   * @param amax       On return, left-right maximum 
   */
  void FetchResults(const TList* list, 
		    const char*  name, 
		    TMultiGraph* thisOther,
		    Int_t        color,
		    const char*  centTxt,
		    Double_t&    max,
		    Double_t&    rmax,
		    Double_t&    amax)
  {
    TH1* dndeta      = FetchResult(list, Form("dndeta%s", name));
    TH1* dndetaMC    = FetchResult(list, Form("dndeta%sMC", name));
    TH1* dndetaTruth = FetchResult(list, "dndetaTruth");
    TH1* dndetaSym   = 0;
    TH1* dndetaMCSym = 0;
    TH1* tracks      = FetchResult(list, "tracks");
    if (tracks) tracks->SetTitle("ALICE Tracks");
    SetAttributes(dndeta,     color);
    SetAttributes(dndetaMC,   color+2);
    SetAttributes(dndetaTruth,color);
    SetAttributes(dndetaSym,  color);
    SetAttributes(dndetaMCSym,color+2);
    SetAttributes(tracks,     color+3);
    ModifyTitle(dndeta,     centTxt);
    ModifyTitle(dndetaMC,   centTxt);
    ModifyTitle(dndetaTruth,centTxt);
    ModifyTitle(dndetaSym,  centTxt);
    ModifyTitle(dndetaMCSym,centTxt);
    ModifyTitle(tracks,     centTxt);
      

    max = TMath::Max(max, AddHistogram(fResults, dndetaTruth, "e5 p"));
    max = TMath::Max(max, AddHistogram(fResults, dndetaMC,    dndetaMCSym));
    max = TMath::Max(max, AddHistogram(fResults, dndeta,      dndetaSym));
    max = TMath::Max(max, AddHistogram(fResults, tracks));
    
    // Info("FetchResults", "Got %p, %p, %p from %s with name %s, max=%f", 
    //      dndeta, dndetaMC, dndetaTruth, list->GetName(), name, max);

    if (fShowLeftRight) {
      fLeftRight->Add(Asymmetry(dndeta,    amax));
      fLeftRight->Add(Asymmetry(dndetaMC,  amax));
      fLeftRight->Add(Asymmetry(tracks,    amax));
    }

    if (thisOther) {
      TIter next(thisOther->GetListOfGraphs());
      TGraph* g = 0;
      while ((g = static_cast<TGraph*>(next()))) {
	fRatios->Add(Ratio(dndeta,    g, rmax));
	fRatios->Add(Ratio(dndetaSym, g, rmax));
	SetAttributes(g, color);
	ModifyTitle(g, centTxt);
	if (!fOthers->GetListOfGraphs() || 
	    !fOthers->GetListOfGraphs()->FindObject(g->GetName())) {
	  max = TMath::Max(max,TMath::MaxElement(g->GetN(), g->GetY()));
	  fOthers->Add(g);
	}
      }
      // fOthers->Add(thisOther);
    }
    if (tracks) { 
      if (!fRatios->GetHists() || 
	  !fRatios->GetHists()->FindObject(tracks->GetName()))
	fRatios->Add(Ratio(dndeta, tracks, rmax));
    }

    if (dndetaTruth) { 
      fRatios->Add(Ratio(dndeta,      dndetaTruth, rmax));
      fRatios->Add(Ratio(dndetaSym,   dndetaTruth, rmax));
      fRatios->Add(Ratio(dndetaMC,    dndetaTruth, rmax));
      fRatios->Add(Ratio(dndetaMCSym, dndetaTruth, rmax));
    }
  }
  //__________________________________________________________________
  /** 
   * Plot the results
   * @param max        Max value 
   * @param rmax       Maximum diviation from 1 of ratios 
   * @param amax       Maximum diviation from 1 of asymmetries 
   */
  void Plot(Double_t     max, 
	    Double_t     rmax,
	    Double_t     amax)
  {
    gStyle->SetOptTitle(0);
    gStyle->SetTitleFont(132, "xyz");
    gStyle->SetLabelFont(132, "xyz");
    
    Int_t    h  = 800;
    Int_t    w  = 800; // h / TMath::Sqrt(2);
    Double_t y1 = 0;
    Double_t y2 = 0;
    Double_t y3 = 0;
    if (!fShowRatios)    w  *= 1.4;
    else                 y1 =  0.3;
    if (!fShowLeftRight) w  *= 1.4;
    else { 
      Double_t y11 = y1;
      y1 = (y11 > 0.0001 ? 0.4 : 0.2);
      y2 = (y11 > 0.0001 ? 0.2 : 0.3);
    }
    TCanvas* c = new TCanvas("Results", "Results", w, h);
    c->SetFillColor(0);
    c->SetBorderSize(0);
    c->SetBorderMode(0);

#if 1
    Info("Plot", "y1=%f, y2=%f, y3=%f extra: %s %s", 
	 y1, y2, y2, (fShowRatios ? "ratios" : ""), 
	 (fShowLeftRight ? "right/left" : ""));
    Info("Plot", "Maximum is %f", max);
#endif
    PlotResults(max, y1);
    c->cd();

    PlotRatios(rmax, y2, y1);
    c->cd( );

    PlotLeftRight(amax, y3, y2);
    c->cd();

    
    Int_t   vMin = fVtxAxis->GetXmin();
    Int_t   vMax = fVtxAxis->GetXmax();    
    TString trg(fTrigString->GetTitle());
    Int_t   nev  = 0;
    if (fTriggers) nev = fTriggers->GetBinContent(1);
    trg          = trg.Strip(TString::kBoth);
    TString base(Form("dndeta_%s_%s_%s_%c%02d%c%02dcm_%09dev",
		      fSysString->GetTitle(), 
		      fSNNString->GetTitle(), 
		      trg.Data(),
		      vMin < 0 ? 'm' : 'p',  TMath::Abs(vMin),
		      vMax < 0 ? 'm' : 'p',  TMath::Abs(vMax),
		      nev));
    c->SaveAs(Form("%s.png",  base.Data()));
    c->SaveAs(Form("%s.root", base.Data()));
    c->SaveAs(Form("%s.C",    base.Data()));
  }
  //__________________________________________________________________
  /** 
   * Build main legend 
   * 
   * @param stack   Stack to include 
   * @param mg      (optional) Multi graph to include 
   * @param x1      Lower X coordinate in the range [0,1]
   * @param y1      Lower Y coordinate in the range [0,1]
   * @param x2      Upper X coordinate in the range [0,1]
   * @param y2 	    Upper Y coordinate in the range [0,1]
   */
  void BuildLegend(THStack* stack, TMultiGraph* mg, 
		   Double_t x1, Double_t y1, Double_t x2, Double_t y2)
  {
    TLegend* l = new TLegend(x1,y1,x2,y2);
    l->SetNColumns(2);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextFont(132);

    TIter    next(stack->GetHists());
    TObject* hist = 0;
    while ((hist = next())) { 
      TString n(hist->GetTitle());
      if (n.Contains("mirrored")) continue;
      l->AddEntry(hist, hist->GetTitle(), "pl");
    }
    if (mg) {
      TIter nexto(mg->GetListOfGraphs());
      while ((hist = nexto())) { 
	TString n(hist->GetTitle());
	if (n.Contains("mirrored")) continue;
	l->AddEntry(hist, hist->GetTitle(), "pl");
      }
    }
    TLegendEntry* d1 = l->AddEntry("d1", "Data", "lp");
    d1->SetLineColor(kBlack);
    d1->SetMarkerColor(kBlack);
    d1->SetMarkerStyle(20);
    TLegendEntry* d2 = l->AddEntry("d2", "Mirrored data", "lp");
    d2->SetLineColor(kBlack);
    d2->SetMarkerColor(kBlack);
    d2->SetMarkerStyle(24);
    
    l->Draw();
  }
  //__________________________________________________________________
  /** 
   * Plot the results
   *    
   * @param max       Maximum 
   * @param yd        Bottom position of pad 
   */
  void PlotResults(Double_t max, Double_t yd) 
  {
    // Make a sub-pad for the result itself
    TPad* p1 = new TPad("p1", "p1", 0, yd, 1.0, 1.0, 0, 0, 0);
    p1->SetTopMargin(0.05);
    p1->SetBorderSize(0);
    p1->SetBorderMode(0);
    p1->SetBottomMargin(yd > 0.001 ? 0.001 : 0.1);
    p1->SetRightMargin(0.05);
    p1->SetGridx();
    p1->SetTicks(1,1);
    p1->SetNumber(1);
    p1->Draw();
    p1->cd();
    
    // Info("PlotResults", "Plotting results with max=%f", max);
    fResults->SetMaximum(1.15*max);
    fResults->SetMinimum(yd > 0.00001 ? -0.1 : 0);

    FixAxis(fResults, 1/(1-yd)/1.7, "#frac{1}{N} #frac{dN_{ch}}{d#eta}");

    p1->Clear();
    fResults->DrawClone("nostack e1");

    fRangeParam->fSlave1Axis = fResults->GetXaxis();
    fRangeParam->fSlave1Pad  = p1;

    // Draw other data
    if (fShowOthers || fShowAlice) {
      TGraphAsymmErrors* o      = 0;
      TIter              nextG(fOthers->GetListOfGraphs());
      while ((o = static_cast<TGraphAsymmErrors*>(nextG())))
        o->DrawClone("same p");
    }

    // Make a legend in the main result pad
    BuildLegend(fResults, fOthers, .15,p1->GetBottomMargin()+.01,.90,.35);
#if 0
    TLegend* l = p1->BuildLegend(.15,p1->GetBottomMargin()+.01,.90,.35);
    l->SetNColumns(2);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextFont(132);
#endif

    // Put a title on top
    fTitle.ReplaceAll("@", " ");
    TLatex* tit = new TLatex(0.10, 0.95, fTitle.Data());
    tit->SetNDC();
    tit->SetTextFont(132);
    tit->SetTextSize(0.05);
    tit->Draw();

    // Put a nice label in the plot
    TString     eS;
    UShort_t    snn = fSNNString->GetUniqueID();
    const char* sys = fSysString->GetTitle();
    if (snn > 1000) eS = Form("%4.2fTeV", float(snn)/1000);
    else            eS = Form("%3dGeV", snn);
    TLatex* tt = new TLatex(.93, .93, Form("%s #sqrt{s}=%s, %s", 
					   sys, 
					   eS.Data(), 
					   fTrigString->GetTitle()));
    tt->SetNDC();
    tt->SetTextFont(132);
    tt->SetTextAlign(33);
    tt->Draw();

    // Put number of accepted events on the plot
    Int_t nev = 0;
    if (fTriggers) nev = fTriggers->GetBinContent(1);
    TLatex* et = new TLatex(.93, .83, Form("%d events", nev));
    et->SetNDC();
    et->SetTextFont(132);
    et->SetTextAlign(33);
    et->Draw();

    // Put number of accepted events on the plot
    if (fVtxAxis) { 
      TLatex* vt = new TLatex(.93, .88, fVtxAxis->GetTitle());
      vt->SetNDC();
      vt->SetTextFont(132);
      vt->SetTextAlign(33);
      vt->Draw();
    }
    // results->Draw("nostack e1 same");

    fRangeParam->fSlave1Axis = FindXAxis(p1, fResults->GetName());
    fRangeParam->fSlave1Pad  = p1;


    // Mark the plot as preliminary
    TLatex* pt = new TLatex(.12, .93, "Preliminary");
    pt->SetNDC();
    pt->SetTextFont(22);
    pt->SetTextSize(0.07);
    pt->SetTextColor(kRed+1);
    pt->SetTextAlign(13);
    pt->Draw();

    if (!gSystem->AccessPathName("ALICE.png")) { 
      TPad* logo = new TPad("logo", "logo", .12, .65, .25, .85, 0, 0, 0);
      logo->SetFillStyle(0);
      logo->Draw();
      logo->cd();
      TImage* i = TImage::Create();
      i->ReadImage("ALICE.png");
      i->Draw();
    }
    p1->cd();
  }
  //__________________________________________________________________
  /** 
   * Plot the ratios 
   * 
   * @param max     Maximum diviation from 1 
   * @param y1      Lower y coordinate of pad
   * @param y2      Upper y coordinate of pad
   */
  void PlotRatios(Double_t max, Double_t y1, Double_t y2) 
  {
    if (!fShowRatios) return;

    bool isBottom = (y1 < 0.0001);
    Double_t yd = y2 - y1;
    // Make a sub-pad for the result itself
    TPad* p2 = new TPad("p2", "p2", 0, y1, 1.0, y2, 0, 0, 0);
    p2->SetTopMargin(0.001);
    p2->SetRightMargin(0.05);
    p2->SetBottomMargin(isBottom ? 1/yd * 0.07 : 0.0001);
    p2->SetGridx();
    p2->SetTicks(1,1);
    p2->SetNumber(2);
    p2->Draw();
    p2->cd();

    // Fix up axis
    FixAxis(fRatios, 1/yd/1.7, "Ratios", 7);

    fRatios->SetMaximum(1+TMath::Max(.22,1.05*max));
    fRatios->SetMinimum(1-TMath::Max(.32,1.05*max));
    p2->Clear();
    fRatios->DrawClone("nostack e1");

    
    // Make a legend
    BuildLegend(fRatios, 0, .15,p2->GetBottomMargin()+.01,.9,
		isBottom ? .6 : .4);
#if 0
    TLegend* l2 = p2->BuildLegend(.15,p2->GetBottomMargin()+.01,.9,
				  isBottom ? .6 : .4);
    l2->SetNColumns(2);
    l2->SetFillColor(0);
    l2->SetFillStyle(0);
    l2->SetBorderSize(0);
    l2->SetTextFont(132);
#endif
    // Make a nice band from 0.9 to 1.1
    TGraphErrors* band = new TGraphErrors(2);
    band->SetPoint(0, fResults->GetXaxis()->GetXmin(), 1);
    band->SetPoint(1, fResults->GetXaxis()->GetXmax(), 1);
    band->SetPointError(0, 0, .1);
    band->SetPointError(1, 0, .1);
    band->SetFillColor(kYellow+2);
    band->SetFillStyle(3002);
    band->SetLineStyle(2);
    band->SetLineWidth(1);
    band->Draw("3 same");
    band->DrawClone("X L same");
    
    // Replot the ratios on top
    fRatios->DrawClone("nostack e1 same");

    if (isBottom) {
      fRangeParam->fMasterAxis = FindXAxis(p2, fRatios->GetName());
      p2->AddExec("range", Form("RangeExec((dNdetaDrawer::RangeParam*)%p)", 
				fRangeParam));
    }
    else { 
      fRangeParam->fSlave2Axis = FindXAxis(p2, fRatios->GetName());
      fRangeParam->fSlave2Pad  = p2;
    }
  }
  //__________________________________________________________________
  /** 
   * Plot the asymmetries
   * 
   * @param max     Maximum diviation from 1 
   * @param y1      Lower y coordinate of pad
   * @param y2      Upper y coordinate of pad
   */
  void PlotLeftRight(Double_t max, Double_t y1, Double_t y2) 
  {
    if (!fShowLeftRight) return;

    bool isBottom = (y1 < 0.0001);
    Double_t yd = y2 - y1;
    // Make a sub-pad for the result itself
    TPad* p3 = new TPad("p3", "p3", 0, y1, 1.0, y2, 0, 0, 0);
    p3->SetTopMargin(0.001);
    p3->SetRightMargin(0.05);
    p3->SetBottomMargin(isBottom ? 1/yd * 0.07 : 0.0001);
    p3->SetGridx();
    p3->SetTicks(1,1);
    p3->SetNumber(2);
    p3->Draw();
    p3->cd();

    // Fix up axis
    FixAxis(fLeftRight, 1/yd/1.7, "Right/Left", 4);

    fLeftRight->SetMaximum(1+TMath::Max(.12,1.05*max));
    fLeftRight->SetMinimum(1-TMath::Max(.15,1.05*max));
    p3->Clear();
    fLeftRight->DrawClone("nostack e1");

    
    // Make a legend
    TLegend* l2 = p3->BuildLegend(.15,p3->GetBottomMargin()+.01,.9,.5);
    l2->SetNColumns(2);
    l2->SetFillColor(0);
    l2->SetFillStyle(0);
    l2->SetBorderSize(0);
    l2->SetTextFont(132);

    // Make a nice band from 0.9 to 1.1
    TGraphErrors* band = new TGraphErrors(2);
    band->SetPoint(0, fResults->GetXaxis()->GetXmin(), 1);
    band->SetPoint(1, fResults->GetXaxis()->GetXmax(), 1);
    band->SetPointError(0, 0, .05);
    band->SetPointError(1, 0, .05);
    band->SetFillColor(kYellow+2);
    band->SetFillStyle(3002);
    band->SetLineStyle(2);
    band->SetLineWidth(1);
    band->Draw("3 same");
    band->DrawClone("X L same");

    fLeftRight->DrawClone("nostack e1 same");
    if (isBottom) {
      fRangeParam->fMasterAxis = FindXAxis(p3, fLeftRight->GetName());
      p3->AddExec("range", Form("RangeExec((dNdetaDrawer::RangeParam*)%p)", 
				fRangeParam));
    }
    else { 
      fRangeParam->fSlave2Axis = FindXAxis(p3, fLeftRight->GetName());
      fRangeParam->fSlave2Pad  = p3;
    }
  }
  /** @} */
  //==================================================================
  /** 
   * @{ 
   * @name Data utility functions 
   */
  /** 
   * Get a result from the passed list
   * 
   * @param list List to search 
   * @param name Object name to search for 
   * 
   * @return 
   */
  TH1* FetchResult(const TList* list, const char* name) const 
  {
    if (!list) return 0;
    
    TH1* ret = static_cast<TH1*>(list->FindObject(name));
    if (!ret) {
      // all->ls();
      Warning("GetResult", "Histogram %s not found", name);
    }

    return ret;
  }
  //__________________________________________________________________
  /** 
   * Add a histogram to the stack after possibly rebinning it  
   * 
   * @param stack   Stack to add to 
   * @param hist    histogram
   * @param option  Draw options 
   * 
   * @return Maximum of histogram 
   */
  Double_t AddHistogram(THStack* stack, TH1* hist, Option_t* option="") const 
  {
    // Check if we have input 
    if (!hist) return 0;

    // Rebin if needed 
    Rebin(hist);

    // Info("AddHistogram", "Adding %s to %s", 
    //      hist->GetName(), stack->GetName());
    stack->Add(hist, option);
    // stack->ls();
    return hist->GetMaximum();
  }
  //__________________________________________________________________
  /** 
   * Add a histogram to the stack after possibly rebinning it  
   * 
   * @param stack   Stack to add to 
   * @param hist    histogram
   * @param option  Draw options 
   * @param sym     On return, the data symmetriced (added to stack)
   * 
   * @return Maximum of histogram 
   */
  Double_t AddHistogram(THStack* stack, TH1* hist, TH1*& sym,
			Option_t* option="") const 
  {
    // Check if we have input 
    if (!hist) return 0;

    // Rebin if needed 
    Rebin(hist);
    stack->Add(hist, option);

    // Now symmetrice the histogram 
    sym = Symmetrice(hist);
    stack->Add(sym, option);

    // Info("AddHistogram", "Adding %s and %s to %s", 
    //      hist->GetName(), sym->GetName(), stack->GetName());
    return hist->GetMaximum();
  }

  //__________________________________________________________________
  /** 
   * Rebin a histogram 
   * 
   * @param h     Histogram to rebin
   * 
   * @return 
   */
  virtual void Rebin(TH1* h) const
  { 
    if (fRebin <= 1) return;

    Int_t nBins = h->GetNbinsX();
    if(nBins % fRebin != 0) {
      Warning("Rebin", "Rebin factor %d is not a devisor of current number "
	      "of bins %d in the histogram %s", fRebin, nBins, h->GetName());
      return;
    }
    
    // Make a copy 
    TH1* tmp = static_cast<TH1*>(h->Clone("tmp"));
    tmp->Rebin(fRebin);
    tmp->SetDirectory(0);
    tmp->Reset();
    // The new number of bins 
    Int_t nBinsNew = nBins / fRebin;
    for(Int_t i = 1;i<= nBinsNew; i++) {
      Double_t content = 0;
      Double_t sumw    = 0;
      Double_t wsum    = 0;
      Int_t    nbins   = 0;
      for(Int_t j = 1; j<=fRebin;j++) {
	Int_t    bin = (i-1)*fRebin + j;
	Double_t c   =  h->GetBinContent(bin);

	if (c <= 0) continue;

	if (fCutEdges) {
	  if (h->GetBinContent(bin+1)<=0 || 
	      h->GetBinContent(bin-1)) {
	    Warning("Rebin", "removing bin %d=%f of %s (%d=%f,%d=%f)", 
		    bin, c, h->GetName(), 
		    bin+1, h->GetBinContent(bin+1), 
		    bin-1, h->GetBinContent(bin-1));
	    continue;
	  }	
	}
	Double_t e =  h->GetBinError(bin);
	Double_t w =  1 / (e*e); // 1/c/c
	content    += c;
	sumw       += w;
	wsum       += w * c;
	nbins++;
      }
      
      if(content > 0 && nbins > 1 ) {
	tmp->SetBinContent(i, wsum / sumw);
	tmp->SetBinError(i,1./TMath::Sqrt(sumw));
      }
    }

    // Finally, rebin the histogram, and set new content
    h->Rebin(fRebin);
    h->Reset();
    for(Int_t i = 1; i<= nBinsNew; i++) {
      h->SetBinContent(i,tmp->GetBinContent(i));
      h->SetBinError(i,  tmp->GetBinError(i));
    }
    
    delete tmp;
  }
  //__________________________________________________________________
  /** 
   * Make an extension of @a h to make it symmetric about 0 
   * 
   * @param h Histogram to symmertrice 
   * 
   * @return Symmetric extension of @a h 
   */
  TH1* Symmetrice(const TH1* h) const
  {
    Int_t nBins = h->GetNbinsX();
    TH1* s = static_cast<TH1*>(h->Clone(Form("%s_mirror", h->GetName())));
    s->SetTitle(Form("%s (mirrored)", h->GetTitle()));
    s->Reset();
    s->SetBins(nBins, -h->GetXaxis()->GetXmax(), -h->GetXaxis()->GetXmin());
    s->SetMarkerColor(h->GetMarkerColor());
    s->SetMarkerSize(h->GetMarkerSize());
    s->SetMarkerStyle(h->GetMarkerStyle()+4);
    s->SetFillColor(h->GetFillColor());
    s->SetFillStyle(h->GetFillStyle());
    s->SetDirectory(0);

    // Find the first and last bin with data 
    Int_t first = nBins+1;
    Int_t last  = 0;
    for (Int_t i = 1; i <= nBins; i++) { 
      if (h->GetBinContent(i) <= 0) continue; 
      first = TMath::Min(first, i);
      last  = TMath::Max(last,  i);
    }
    
    Double_t xfirst = h->GetBinCenter(first-1);
    Int_t    f1     = h->GetXaxis()->FindBin(-xfirst);
    Int_t    l2     = s->GetXaxis()->FindBin(xfirst);
    for (Int_t i = f1, j=l2; i <= last; i++,j--) { 
      s->SetBinContent(j, h->GetBinContent(i));
      s->SetBinError(j, h->GetBinError(i));
    }
    // Fill in overlap bin 
    s->SetBinContent(l2+1, h->GetBinContent(first));
    s->SetBinError(l2+1, h->GetBinError(first));
    return s;
  }
  //__________________________________________________________________
  /** 
   * Calculate the left-right asymmetry of input histogram 
   * 
   * @param h   Input histogram
   * @param max On return, the maximum distance from 1 of the histogram
   * 
   * @return Asymmetry 
   */
  TH1* Asymmetry(TH1* h, Double_t& max)
  {
    if (!h) return 0;

    TH1* ret = static_cast<TH1*>(h->Clone(Form("%s_leftright", h->GetName())));
    // Int_t    oBins = h->GetNbinsX();
    // Double_t high  = h->GetXaxis()->GetXmax();
    // Double_t low   = h->GetXaxis()->GetXmin();
    // Double_t dBin  = (high - low) / oBins;
    // Int_t    tBins = Int_t(2*high/dBin+.5);
    // ret->SetBins(tBins, -high, high);
    ret->SetDirectory(0);
    ret->Reset();
    ret->SetTitle(Form("%s (+/-)", h->GetTitle()));
    ret->SetYTitle("Right/Left");
    Int_t nBins = h->GetNbinsX();
    for (Int_t i = 1; i <= nBins; i++) { 
      Double_t x = h->GetBinCenter(i);
      if (x > 0) break;
      
      Double_t c1 = h->GetBinContent(i);
      Double_t e1 = h->GetBinError(i);
      if (c1 <= 0) continue; 
      
      Int_t    j  = h->FindBin(-x);
      if (j <= 0 || j > nBins) continue;

      Double_t c2 = h->GetBinContent(j);
      Double_t e2 = h->GetBinError(j);

      Double_t c12 = c1*c1;
      Double_t e   = TMath::Sqrt((e2*e2*c1*c1+e1*e1*c2*c2)/(c12*c12));
      
      Int_t    k   = ret->FindBin(x);
      ret->SetBinContent(k, c2/c1);
      ret->SetBinError(k, e);
    }
    max = TMath::Max(max, RatioMax(ret));

    return ret;
  }
  //__________________________________________________________________
  /** 
   * Transform a graph into a histogram 
   * 
   * @param g 
   * 
   * @return 
   */
  TH1* Graph2Hist(const TGraphAsymmErrors* g) const
  {
    Int_t    nBins = g->GetN();
    TArrayF  bins(nBins+1);
    Double_t dx = 0;
    for (Int_t i = 0; i < nBins; i++) { 
      Double_t x   = g->GetX()[i];
      Double_t exl = g->GetEXlow()[i];
      Double_t exh = g->GetEXhigh()[i];
      bins.fArray[i]   = x-exl;
      bins.fArray[i+1] = x+exh;
      Double_t dxi = exh+exl;
      if (i == 0) dx  = dxi;
      else if (dxi != dx) dx = 0;
    }
    TString name(g->GetName());
    TString title(g->GetTitle());
    TH1D* h = 0;
    if (dx != 0) {
      h = new TH1D(name.Data(), title.Data(), nBins, bins[0], bins[nBins]);
    }
    else {
      h = new TH1D(name.Data(), title.Data(), nBins, bins.fArray);
    }
    h->SetMarkerStyle(g->GetMarkerStyle());
    h->SetMarkerColor(g->GetMarkerColor());
    h->SetMarkerSize(g->GetMarkerSize());
    
    return h;
  }
  /* @} */
  //==================================================================
  /** 
   * @{ 
   * @name Ratio utility functions 
   */
  /** 
   * Get the maximum diviation from 1 in the passed ratio
   * 
   * @param h Ratio histogram
   * 
   * @return Max diviation from 1 
   */
  Double_t RatioMax(TH1* h) const
  {
    Int_t    nBins = h->GetNbinsX();
    Double_t ret   = 0;
    for (Int_t i = 1; i <= nBins; i++) { 
      Double_t c = h->GetBinContent(i);
      if (c == 0) continue;
      Double_t e = h->GetBinError(i);
      Double_t d = TMath::Abs(1-c-e);
      ret        = TMath::Max(d, ret);
    }
    return ret;
  }
  //__________________________________________________________________
  /** 
   * Make the ratio of h1 to h2 
   * 
   * @param o1  First object (numerator) 
   * @param o2  Second object (denominator)
   * @param max Maximum diviation from 1 
   * 
   * @return o1 / o2
   */
  TH1* Ratio(const TObject* o1, const TObject* o2, Double_t& max) const
  {
    if (!o1 || !o2) return 0;

    TH1* r = 0;
    const TH1* h1 = dynamic_cast<const TH1*>(o1); 
    if (h1) { 
      // Info("Ratio", "First is a TH1");
      const TH1* h2 = dynamic_cast<const TH1*>(o2); 
      if (h2) { 
	// Info("Ratio", "Second is a TH1");
	r = RatioHH(h1,h2,max);
      }
      else {
	const TGraph* g2 = dynamic_cast<const TGraph*>(o2);
	if (g2) { 
	  // Info("Ratio", "Second os a TGraph");
	  r = RatioHG(h1,g2,max);      
	}
      }
    }
    else {
      const TGraphAsymmErrors* g1 = dynamic_cast<const TGraphAsymmErrors*>(o1);
      if (g1) { 
	// Info("Ratio", "First is a TGraphAsymmErrors");
	const TGraphAsymmErrors* g2 = 
	  dynamic_cast<const TGraphAsymmErrors*>(o2);
	if (g2) {
	  // Info("Ratio", "Second is a TGraphAsymmErrors");
	  r = RatioGG(g1, g2, max);
	}
      }
    }
    if (!r) {
      Warning("Ratio", "Don't know how to divide a %s (%s) with a %s (%s)", 
	      o1->ClassName(), o1->GetName(), o2->ClassName(), o2->GetName());
      return 0;
    }
    // Check that the histogram isn't empty
    if (r->GetEntries() <= 0) { 
      delete r; 
      r = 0; 
    }
    // for (Int_t bin = 1; bin <= r->GetNbinsX(); bin++) 
    //   if (r->GetBinContent(bin) != 0) return r;
      
    return r;
  }
  //__________________________________________________________________
  /** 
   * Compute the ratio of @a h to @a g.  @a g is evaluated at the bin
   * centers of @a h 
   * 
   * @param h  Numerator 
   * @param g  Divisor 
   * @param max Maximum diviation from 1 
   * 
   * @return h/g 
   */
  TH1* RatioHG(const TH1* h, const TGraph* g, Double_t& max) const 
  {
    if (!h || !g) return 0;

    TH1* ret = static_cast<TH1*>(h->Clone("tmp"));
    ret->SetName(Form("%s_over_%s", h->GetName(), g->GetName()));
    ret->SetTitle(Form("%s / %s", h->GetTitle(), g->GetTitle()));
    ret->Reset();
    ret->SetMarkerStyle(h->GetMarkerStyle());
    ret->SetMarkerColor(g->GetMarkerColor());
    ret->SetMarkerSize(0.9*h->GetMarkerSize());
    // ret->SetMarkerStyle(g->GetMarkerStyle());
    // ret->SetMarkerColor(h->GetMarkerColor());
    // ret->SetMarkerSize(0.9*g->GetMarkerSize());
    ret->SetDirectory(0);
    Double_t xlow  = g->GetX()[0];
    Double_t xhigh = g->GetX()[g->GetN()-1];
    if (xlow > xhigh) { Double_t t = xhigh; xhigh = xlow; xlow = t; }

    for (Int_t i = 1; i <= h->GetNbinsX(); i++) { 
      Double_t c = h->GetBinContent(i);
      if (c <= 0) continue;

      Double_t x = h->GetBinCenter(i);
      if (x < xlow || x > xhigh) continue; 

      Double_t f = g->Eval(x);
      if (f <= 0) continue; 

      ret->SetBinContent(i, c / f);
      ret->SetBinError(i, h->GetBinError(i) / f);
    }
    if (ret->GetEntries() > 0) 
      max = TMath::Max(RatioMax(ret), max);
    return ret;
  }
  //__________________________________________________________________
  /** 
   * Make the ratio of h1 to h2 
   * 
   * @param h1 First histogram (numerator) 
   * @param h2 Second histogram (denominator)
   * @param max Maximum diviation from 1 
   * 
   * @return h1 / h2
   */
  TH1* RatioHH(const TH1* h1, const TH1* h2, Double_t& max) const
  {
    if (!h1 || !h2) return 0;
    TH1* t1 = static_cast<TH1*>(h1->Clone(Form("%s_%s", 
					       h1->GetName(), 
					       h2->GetName())));
    t1->SetTitle(Form("%s / %s", h1->GetTitle(), h2->GetTitle()));
    t1->Divide(h2);
    // t1->SetMarkerColor(h1->GetMarkerColor());
    // t1->SetMarkerStyle(h2->GetMarkerStyle());
    t1->SetMarkerColor(h2->GetMarkerColor());
    t1->SetMarkerStyle(h1->GetMarkerStyle());
    t1->SetMarkerSize(0.9*h1->GetMarkerSize());
    t1->SetDirectory(0);
    max = TMath::Max(RatioMax(t1), max);
    return t1;
  }
  //__________________________________________________________________
  /** 
   * Calculate the ratio of two graphs - g1 / g2
   * 
   * @param g1 Numerator 
   * @param g2 Denominator
   * @param max Maximum diviation from 1 
   * 
   * @return g1 / g2 in a histogram 
   */
  TH1* RatioGG(const TGraphAsymmErrors* g1, 
	       const TGraphAsymmErrors* g2, Double_t& max) const
  {
    Int_t    nBins = g1->GetN();
    TArrayF  bins(nBins+1);
    Double_t dx   = 0;
    for (Int_t i = 0; i < nBins; i++) {
      Double_t x   = g1->GetX()[i];
      Double_t exl = g1->GetEXlow()[i];
      Double_t exh = g1->GetEXhigh()[i];
      bins.fArray[i]   = x-exl;
      bins.fArray[i+1] = x+exh;
      Double_t dxi = exh+exl;
      if (i == 0) dx  = dxi;
      else if (dxi != dx) dx = 0;
    }
    TString name(Form("%s_%s", g1->GetName(), g2->GetName()));
    TString title(Form("%s / %s", g1->GetTitle(), g2->GetTitle()));
    TH1* h = 0;
    if (dx != 0) {
      h = new TH1F(name.Data(), title.Data(), nBins, bins[0], bins[nBins]);
    }
    else {
      h = new TH1F(name.Data(), title.Data(), nBins, bins.fArray);
    }
    h->SetMarkerStyle(g1->GetMarkerStyle());
    h->SetMarkerColor(g2->GetMarkerColor());
    h->SetMarkerSize(0.9*g1->GetMarkerSize());
    // h->SetMarkerStyle(g2->GetMarkerStyle());
    // h->SetMarkerColor(g1->GetMarkerColor());
    // h->SetMarkerSize(0.9*g2->GetMarkerSize());
    h->SetDirectory(0);

    Double_t low  = g2->GetX()[0];
    Double_t high = g2->GetX()[g2->GetN()-1];
    if (low > high) { Double_t t = low; low = high; high = t; }
    for (Int_t i = 0; i < nBins; i++) { 
      Double_t x  = g1->GetX()[i];
      if (x < low-dx || x > high+dx) continue;
      Double_t c1 = g1->GetY()[i];
      Double_t e1 = g1->GetErrorY(i);
      Double_t c2 = g2->Eval(x);
      
      h->SetBinContent(i+1, c1 / c2);
      h->SetBinError(i+1, e1 / c2);
    }
    max = TMath::Max(RatioMax(h), max);
    return h;
  }  
  /* @} */
  //==================================================================
  /** 
   * @{ 
   * @name Graphics utility functions 
   */
  /** 
   * Find an X axis in a pad 
   * 
   * @param p     Pad
   * @param name  Histogram to find axis for 
   * 
   * @return Found axis or null
   */
  TAxis* FindXAxis(TVirtualPad* p, const char* name)
  {
    TObject* o = p->GetListOfPrimitives()->FindObject(name);
    if (!o) { 
      Warning("FindXAxis", "%s not found in pad", name);
      return 0;
    }
    THStack* stack = dynamic_cast<THStack*>(o);
    if (!stack) { 
      Warning("FindXAxis", "%s is not a THStack", name);
      return 0;
    }
    if (!stack->GetHistogram()) { 
      Warning("FindXAxis", "%s has no histogram", name);
      return 0;
    }
    TAxis* ret = stack->GetHistogram()->GetXaxis();
    return ret;
  }

  //__________________________________________________________________
  /**
   * Fix the apperance of the axis in a stack
   *
   * @param stack  stack of histogram
   * @param s      Scaling factor
   * @param ytitle Y axis title
   * @param force  Whether to draw the stack first or not
   * @param ynDiv  Divisions on Y axis
   */
  void FixAxis(THStack* stack, Double_t s, const char* ytitle,
               Int_t ynDiv=210, Bool_t force=true)
  {
    if (!stack) { 
      Warning("FixAxis", "No stack passed for %s!", ytitle);
      return;
    }
    if (force) stack->Draw("nostack e1");

    TH1* h = stack->GetHistogram();
    if (!h) { 
      Warning("FixAxis", "Stack %s has no histogram", stack->GetName());
      return;
    }
    
    h->SetXTitle("#eta");
    h->SetYTitle(ytitle);
    TAxis* xa = h->GetXaxis();
    TAxis* ya = h->GetYaxis();
    if (xa) {
      xa->SetTitle("#eta");
      // xa->SetTicks("+-");
      xa->SetTitleSize(s*xa->GetTitleSize());
      xa->SetLabelSize(s*xa->GetLabelSize());
      xa->SetTickLength(s*xa->GetTickLength());

      if (stack != fResults) {
	TAxis* rxa = fResults->GetXaxis();
	xa->Set(rxa->GetNbins(), rxa->GetXmin(), rxa->GetXmax());
      }
    }
    if (ya) {
      ya->SetTitle(ytitle);
      ya->SetDecimals();
      // ya->SetTicks("+-");
      ya->SetNdivisions(ynDiv);
      ya->SetTitleSize(s*ya->GetTitleSize());
      ya->SetTitleOffset(ya->GetTitleOffset()/s);
      ya->SetLabelSize(s*ya->GetLabelSize());
    }
  }
  /* @} */



  //__________________________________________________________________
  Bool_t       fShowOthers;   // Show other data
  Bool_t       fShowAlice;    // Show ALICE published data
  Bool_t       fShowRatios;   // Show ratios 
  Bool_t       fShowLeftRight;// Show asymmetry 
  UShort_t     fRebin;        // Rebinning factor 
  Bool_t       fCutEdges;     // Whether to cut edges
  TString      fTitle;        // Title on plot
  TNamed*      fTrigString;   // Trigger string (read, or set)
  TNamed*      fNormString;   // Normalisation string (read, or set)
  TNamed*      fSNNString;    // Energy string (read, or set)
  TNamed*      fSysString;    // Collision system string (read or set)
  TAxis*       fVtxAxis;      // Vertex cuts (read or set)
  TAxis*       fCentAxis;     // Centrality axis
  THStack*     fResults;      // Stack of results 
  THStack*     fRatios;       // Stack of ratios 
  THStack*     fLeftRight;    // Left-right asymmetry
  TMultiGraph* fOthers;       // Older data 
  TH1*         fTriggers;     // Number of triggers
  RangeParam*  fRangeParam;   // Parameter object for range zoom 
  
};

//=== Stuff for auto zooming =========================================
void UpdateRange(dNdetaDrawer::RangeParam* p)
{
  if (!p) { 
    Warning("UpdateRange", "No parameters %p", p);
    return;
  }
  if (!p->fMasterAxis) { 
    Warning("UpdateRange", "No master axis %p", p->fMasterAxis);
    return;
  }
  Int_t    first = p->fMasterAxis->GetFirst();
  Int_t    last  = p->fMasterAxis->GetLast();
  Double_t x1    = p->fMasterAxis->GetBinCenter(first);
  Double_t x2    = p->fMasterAxis->GetBinCenter(last);
  //Info("UpdateRange", "Range set to [%3d,%3d]->[%f,%f]", first, last, x1,x2);

  if (p->fSlave1Axis) { 
    Int_t i1 = p->fSlave1Axis->FindBin(x1);
    Int_t i2 = p->fSlave1Axis->FindBin(x2);
    p->fSlave1Axis->SetRange(i1, i2);
    p->fSlave1Pad->Modified();
    p->fSlave1Pad->Update();
  }
  if (p->fSlave2Axis) { 
    Int_t i1 = p->fSlave2Axis->FindBin(x1);
    Int_t i2 = p->fSlave2Axis->FindBin(x2);
    p->fSlave2Axis->SetRange(i1, i2);
    p->fSlave2Pad->Modified();
    p->fSlave2Pad->Update();
  }
  TCanvas*  c = gPad->GetCanvas();
  c->cd();
}
  
//____________________________________________________________________
void RangeExec(dNdetaDrawer::RangeParam* p)
{
  // Event types: 
  //  51:   Mouse move 
  //  53:   
  //  1:    Button down 
  //  21:   Mouse drag
  //  11:   Mouse release 
  // dNdetaDrawer::RangeParam* p = 
  //   reinterpret_cast<dNdetaDrawer::RangeParam*>(addr);
  Int_t event     = gPad->GetEvent();
  TObject *select = gPad->GetSelected();
  if (event == 53) { 
    UpdateRange(p);
    return;
  }
  if (event != 11 || !select || select != p->fMasterAxis) return;
  UpdateRange(p);
}

//=== Steering function ==============================================  
/** 
 * Draw @f$ dN/d\eta@f$ 
 * 
 * @param filename  File name 
 * @param flags     Flags 
 * @param title     Title 
 * @param rebin     Rebinning factor 
 * @param cutEdges  Whether to cut edges when rebinning
 * @param sNN       (optional) Collision energy [GeV]
 * @param sys       (optional) Collision system (1: pp, 2: PbPb)
 * @param trg       (optional) Trigger (1: INEL, 2: INEL>0, 4: NSD)   
 * @param vzMin     Least @f$ v_z@f$
 * @param vzMax     Largest @f$ v_z@f$
 *
 * @ingroup pwg2_forward_dndeta
 */
void
DrawdNdeta(const char* filename="forward_dndeta.root", 
	   Int_t       flags=0xf,
	   const char* title="",
	   UShort_t    rebin=5, 
	   Bool_t      cutEdges=false,
	   UShort_t    sNN=0, 
	   UShort_t    sys=0,
	   UShort_t    trg=0,
	   Float_t     vzMin=999, 
	   Float_t     vzMax=-999)
{
  dNdetaDrawer* pd = new dNdetaDrawer;
  dNdetaDrawer& d = *pd;
  d.SetRebin(rebin);
  d.SetCutEdges(cutEdges);
  d.SetTitle(title);
  d.SetShowOthers(flags & 0x1);
  d.SetShowAlice(flags & 0x2);
  d.SetShowRatios(flags & 0x4);
  d.SetShowLeftRight(flags & 0x8);
  // Do the below if your input data does not contain these settings 
  if (sNN > 0) d.SetSNN(sNN);     // Collision energy per nucleon pair (GeV)
  if (sys > 0) d.SetSys(sys);     // Collision system (1:pp, 2:PbPB)
  if (trg > 0) d.SetTrigger(trg); // Collision trigger (1:INEL, 2:INEL>0, 4:NSD)
  if (vzMin < 999 && vzMax > -999) 
    d.SetVertexRange(vzMin,vzMax); // Collision vertex range (cm)
  d.Run(filename);
}
//____________________________________________________________________
//
// EOF
//

