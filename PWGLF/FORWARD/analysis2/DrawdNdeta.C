/**
 * @file   DrawdNdeta.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 14:07:10 2011
 * 
 * @brief  Script to visualise the dN/deta for pp and PbPb
 *
 * This script is independent of any AliROOT code - and should stay
 * that way.
 * 
 * The script is <i>very</i> long - sigh - the joy of drawing
 * things nicely in ROOT
 * 
 * @ingroup pwglf_forward_dndeta
 */
#include <TH1.h>
#include <TColor.h>
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
#include <TRandom.h>
#include <TParameter.h>
#include <TGClient.h>
#include <fstream>
#include <iostream>
/** Systematic error color */
// #define SYSERR_COLOR kGray;
// #define SYSERR_COLOR kBlue-10
// #define SYSERR_COLOR kCyan-10
#define SYSERR_COLOR TColor::GetColor(220, 220, 255)
/** Systematic error style */
#define SYSERR_STYLE 1001

Double_t myFunc(Double_t* xp, Double_t* pp);

/**
 * Class to draw dN/deta results 
 * 
 * @ingroup pwglf_forward_tasks_dndeta
 * @ingroup pwglf_forward_dndeta
 */
struct dNdetaDrawer 
{
  enum EFlags { 
    kShowRatios    = 0x0001, 
    kShowLeftRight = 0x0002, 
    kShowSysError  = 0x0004, 
    kShowRings     = 0x0008,
    kCutEdges      = 0x0010,
    kRemoveOuters  = 0x0020, 
    kUseFinalMC    = 0x0040,
    kUseEmpirical  = 0x0080,
    kForceMB       = 0x0100,
    kMirror        = 0x0200,
    kExport        = 0x0400, 
    kAddExec       = 0x0800,
    kOldFormat     = 0x1000,
    kVerbose       = 0x2000,
    kHiRes         = 0x4000,
    kExtraWhite    = 0x8000,
    kDefaultOptions = 0xCE07
  };
  enum EOutFormat { 
    kPNG     = 0x1, 
    kPDF     = 0x2, 
    kROOT    = 0x4, 
    kScript  = 0x8,
    kAllFormats = 0xF
  };
  
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
    : fOptions(kDefaultOptions),
      fFormats(kAllFormats),
      fShowOthers(0),        // Show other data
      // Settings 
      fRebin(0),             // Rebinning factor 
      fFwdSysErr(0.076),     // Systematic error in forward range
      fCenSysErr(0),         // Systematic error in central range 
      fTitle(""),            // Title on plot
      fBase(""),             // Optional base name of output files
      fClusterScale(""),     // Scaling of clusters to tracklets      
      // Read (or set) information 
      fTrigString(0),        // Trigger string (read, or set)
      fNormString(0),        // Normalisation string (read, or set)
      fSNNString(0),         // Energy string (read, or set)
      fSysString(0),         // Collision system string (read or set)
      fVtxAxis(0),           // Vertex cuts (read or set)
      fCentAxis(0),          // Centrality axis
      fTriggerEff(1),        // Trigger efficency 
      fExtTriggerEff(false), // True if fTriggerEff was read 
      fCentMin(0),           // Least centrality to plot
      fCentMax(100),         // Largest centrality to plot
      // Resulting plots 
      fResults(0),           // Stack of results 
      fRatios(0),            // Stack of ratios 
      fLeftRight(0),         // Left-right asymmetry
      fOthers(0),            // Older data 
      fTriggers(0),          // Number of triggers
      fTruth(0),             // Pointer to truth 
      fRangeParam(0)         // Parameter object for range zoom 
  {
    fRangeParam = new RangeParam;
    fRangeParam->fMasterAxis = 0;
    fRangeParam->fSlave1Axis = 0;
    fRangeParam->fSlave1Pad  = 0;
    fRangeParam->fSlave2Axis = 0;
    fRangeParam->fSlave2Pad  = 0;

    TColor* sysErr = gROOT->GetColor(kSysErrColor);
    sysErr->SetAlpha(0.7);
  }
  /** 
   * Cpoy constructor 
   */
  dNdetaDrawer(const dNdetaDrawer&) {}
  /** 
   * Assignment operator
   * 
   * 
   * @return Reference to this object
   */
  dNdetaDrawer& operator=(const dNdetaDrawer&) { return *this; }

  //__________________________________________________________________
  /** 
   * Destructor 
   */
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
  void SetShowOthers(UInt_t others) { fShowOthers = others; }
  /** 
   * Set the rebinning factor 
   * 
   * @param x Rebinning factor (must be a divisor in the number of bins) 
   */
  void SetRebin(UShort_t x)       { fRebin = x; }
  //__________________________________________________________________
  /** 
   * Set the title of the plot
   * 
   * @param x Title
   */
  void SetTitle(TString x)        { fTitle = x; }
  //__________________________________________________________________
  /** 
   * Set the base name of the output files 
   * 
   * @param x Base name 
   */
  void SetBase(TString x) { fBase = x; }
  //__________________________________________________________________
  /** 
   * Set the systematic error in the forward region
   * 
   * @param e Systematic error in the forward region 
   */
  void SetForwardSysError(Double_t e=0) { fFwdSysErr = e; }
  //__________________________________________________________________
  /** 
   * Set the systematic error in the forward region
   * 
   * @param e Systematic error in the forward region 
   */
  void SetCentralSysError(Double_t e=0) { fCenSysErr = e; }
  /** 
   * Set the 'Final MC' correction file.  This is needed if the
   * secondary maps where produced using the old code
   * 
   * @param file Filename 
   */
  void SetFinalMC(const TString& file) { fFinalMC = file; }
  /** 
   * Set the file that contains the empirical correction.  This is
   * needed when the secondary maps was generated with an in-accurate
   * geometry, and when we're analysing data at nominal interaction
   * points
   * 
   * @param file Filename 
   */
  void SetEmpirical(const TString& file) { fEmpirical = file; }
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
				    sys == 2 ? "PbPb" : 
				    sys == 3 ? "pPb" : 
				    "unknown"));
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
  /** 
   * Set the centrality range in centimeters 
   * 
   * @param centMin Min @f$ v_z@f$
   * @param centMax Max @f$ v_z@f$
   */
  void SetCentralityRange(UShort_t centMin, UShort_t centMax) 
  {
    fCentMin = centMin;
    fCentMax = centMax;
  }
  //__________________________________________________________________
  /** 
   * Set the trigger mask (overrides what's in the file)
   * 
   * @param trig Trigger mask (0x1: INEL, 0x2: INEL>0, 0x4: NSD)
   */
  void SetTrigger(UShort_t trig)
  {
    fTrigString = new TNamed("trigString", (trig & 0x1 ? "INEL" : 
					    trig & 0x2 ? "INEL>0" : 
					    trig & 0x4 ? "NSD" : 
					    "unknown"));
    fTrigString->SetUniqueID(trig);
  }
  //__________________________________________________________________
  /** 
   * Set the trigger efficiency - if set, then scale result histograms 
   * by this factor 
   * 
   * @param eff @f$\varepsilon_{T}@f$ 
   */
  void SetTriggerEfficiency(Float_t eff)  
  { 
    fTriggerEff = eff; 
    fExtTriggerEff = false;
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
  void Run(const char* filename="forward_dndeta.root",
	   UInt_t flags=kDefaultOptions, UInt_t formats=kAllFormats) 
  {
    fOptions          = flags;
    fFormats          = formats;
    SetForwardSysError(flags & kShowSysError ? 0.076 : 0);
    SetFinalMC        (flags & kUseFinalMC ? "forward_dndetamc.root" : "");
    SetEmpirical      (flags & kUseEmpirical ? "EmpiricalCorrection.root" : "");
    // SetBase(base);

    Double_t max = 0, rmax=0, amax=0;

    gStyle->SetPalette(1);

    // --- Open input file -------------------------------------------
    TFile* file = TFile::Open(filename, "READ");
    if (!file) { 
      Error("Run", "Cannot open %s", filename);
      return;
    }
    Info("Run", "Drawing results from %s", file->GetName());

    // --- Get forward list ------------------------------------------
    TList* forward = static_cast<TList*>(file->Get("ForwarddNdetaResults"));
    if (!forward) { 
      Error("Run", "Couldn't find list ForwarddNdetaResults");
      return;
    }
    TList* sums = static_cast<TList*>(file->Get("ForwarddNdetaSums"));
    if (!sums) { 
      Error("Run", "Couldn't find list ForwarddNdetaSums");
      return;
    }
    TParameter<bool>* p = 
      static_cast<TParameter<bool>*>(sums->FindObject("empirical"));
    if (p && p->GetVal() && !fEmpirical.IsNull()) {
      Warning("Run", "Empirical correction already applied");
      fEmpirical = "__task__";
    }
    // --- Get information on the run --------------------------------
    FetchInformation(forward);

    // --- Print settings --------------------------------------------
    Info("Run", "Settings for the drawer:\n"
	 "   Show ratios:                      %5s\n"
	 "   Show Left/right:                  %5s\n"
	 "   Show rings:                       %5s\n"
	 "   Export to file:                   %5s\n"
	 "   Cut edges when rebinning:         %5s\n"
	 "   Remove outer rings:               %5s\n"
	 "   Mirror to un-covered regions:     %5s\n"
	 "   Force minimum bias:               %5s\n"
	 "   Show other results:               0x%03x\n"
	 "   Rebinning factor:                 %5d\n"
	 "   Forward systematic error:         %5.1f%%\n"
	 "   Central systematic error:         %5.1f%%\n"
	 "   Title on plot:                    %s\n"
	 "   Scaling of clusters to tracklets: %s\n"
	 "   Final MC correction file:         %s\n"
	 "   Empirical correction file:        %s",
	 ((fOptions & kShowRatios)    ? "yes" : "no"), 
	 ((fOptions & kShowLeftRight) ? "yes" : "no"),
	 ((fOptions & kShowRings)     ? "yes" : "no"),
	 ((fOptions & kExport)        ? "yes" : "no"), 
	 ((fOptions & kCutEdges)      ? "yes" : "no"),
	 ((fOptions & kRemoveOuters)  ? "yes" : "no"),
	 ((fOptions & kMirror)        ? "yes" : "no"),
	 ((fOptions & kForceMB)       ? "yes" : "no"),
	 fShowOthers, fRebin, (100*fFwdSysErr), (100*fCenSysErr), 
	 fTitle.Data(), fClusterScale.Data(), fFinalMC.Data(), 
	 fEmpirical.Data());

    // --- Set the macro pathand load other data script --------------
    TString savPath(gROOT->GetMacroPath());
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));
    // Always recompile 
    if (!gROOT->GetClass("RefData"))
      gROOT->LoadMacro("OtherData.C++");
    gROOT->SetMacroPath(savPath);

    // --- Get the central results -----------------------------------
    TList* clusters = static_cast<TList*>(file->Get("CentraldNdetaResults"));
    if (!clusters) Warning("Run", "Couldn't find list CentraldNdetaResults");

    // --- Get the central results -----------------------------------
    TList* mcTruth = static_cast<TList*>(file->Get("MCTruthdNdetaResults"));
    if (!mcTruth) Warning("Run", "Couldn't find list MCTruthdNdetaResults");

    // --- Make our containtes ---------------------------------------
    fResults   = new THStack("results", "Results");
    fRatios    = new THStack("ratios",  "Ratios");
    fLeftRight = new THStack("asymmetry", "Left-right asymmetry");
    fOthers    = new TMultiGraph();

    // --- Try to open the final MC file, and find relevant lists ----
    TList* forwardMC = 0;
    // TList* centralMC = 0;
    if (!fFinalMC.IsNull()) { 
      TFile* finalMC = TFile::Open(fFinalMC, "READ");
      if (!finalMC) { 
	Warning("Run", "Failed to open file %s for final MC corrections", 
		fFinalMC.Data());
      }
      else { 
	forwardMC = static_cast<TList*>(finalMC->Get("ForwarddNdetaResults"));
	if (!forwardMC) 
	  Warning("Run","Couldn't find list ForwarddNdetaResults for final MC");
#if 0
	centralMC = static_cast<TList*>(finalMC->Get("CentradNdetalResults"));
	if (!centralMC) 
	  Warning("Run","Couldn't find list CentraldNdetaResults for final MC");
#endif
      }
    }
    if (!forwardMC) fFinalMC = "";

    // --- Try to get the emperical correction -----------------------
    TGraphErrors* empCorr = 0;
    if (!fEmpirical.IsNull() && !fEmpirical.EqualTo("__task__")) {
      if (gSystem->AccessPathName(fEmpirical.Data())) { // Not found here
	fEmpirical = 
	  gSystem->ExpandPathName(Form("$ALICE_ROOT/PWGLF/FORWARD/"
				       "corrections/Empirical/%s", 
				       fEmpirical.Data()));
	if (gSystem->AccessPathName(fEmpirical.Data())) { // Not found here
	  Warning("Run", "Couldn't get empirical correction file");
	  fEmpirical = "";
	}
      }
      if (!fEmpirical.IsNull()) {
	TFile* empirical = TFile::Open(fEmpirical, "READ");
	if (!empirical) { 
	  Warning("Run", "couldn't open empirical correction file: %s",
		  fEmpirical.Data());
	  fEmpirical = "";
	}
	const char* empPath = "fmdfull/average";
	empCorr = static_cast<TGraphErrors*>(empirical->Get(empPath));
	if (!empCorr) {
	  Warning("Run", "Didn't find the graph %s in %s", 
		  empPath, fEmpirical.Data());
	  fEmpirical = "";
	}
      }
    }
    if (!empCorr && !fEmpirical.EqualTo("__task__")) fEmpirical = "";

    // --- Loop over input data --------------------------------------
    TObjArray truths;
    FetchTopResults(mcTruth,  0, 0, "MCTruth", max, rmax, amax,truths);
    TObjArray* fwdA = FetchTopResults(forward,  forwardMC, empCorr, "Forward", 
				   max, rmax, amax,truths);
    TObjArray* cenA = FetchTopResults(clusters, 0, 0, "Central", 
				   max, rmax, amax,truths);

    // --- Get trigger information -----------------------------------
    // TList* sums = static_cast<TList*>(file->Get("ForwardSums"));
    if (sums) {
      TList* all = (fOptions & kOldFormat ? sums : 
		    static_cast<TList*>(sums->FindObject("all")));
      if (all) {
	fTriggers = FetchHistogram(all, "triggers");
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
      fShowOthers = 0;
    }
    if (!fRatios->GetHists() || 
	fRatios->GetHists()->GetEntries() <= 0) { 
      Warning("Run", "No ratio data found - disabling that");
      // fRatios->ls();
      fOptions &= ~kShowRatios;
    }
    if (!fLeftRight->GetHists() || 
	fLeftRight->GetHists()->GetEntries() <= 0) { 
      Warning("Run", "No left/right data found - disabling that");
      // fLeftRight->ls();
      fOptions &= ~kShowLeftRight;
    }
    if (fFwdSysErr > 0) { 
      if (fCenSysErr <= 0) fCenSysErr = fFwdSysErr;
      for (Int_t i = 0; i < fwdA->GetEntriesFast(); i++) {
	TH1* fwd = static_cast<TH1*>(fwdA->At(i));
	TH1* cen = static_cast<TH1*>(cenA->At(i));
	CorrectForward(fwd);
	CorrectCentral(cen);
	Double_t low, high;
	TH1* tmp = Merge(cen, fwd, low, high);
	TF1* f   = FitMerged(tmp, low, high);
	MakeSysError(tmp, cen, fwd, f);
	delete f;
	if (fOptions & kVerbose) 
	  Info("", "Adding systematic error histogram %s", tmp->GetName());
	fResults->GetHists()->AddFirst(tmp, "e5");

	if (!(fOptions & kMirror)) continue;

	TH1* tmp2 = Symmetrice(tmp);
	tmp2->SetFillColor(tmp->GetFillColor());
	tmp2->SetFillStyle(tmp->GetFillStyle());
	tmp2->SetMarkerStyle(tmp->GetMarkerStyle());
	tmp2->SetLineWidth(tmp->GetLineWidth());
	fResults->GetHists()->AddFirst(tmp2, "e5");
	fResults->Modified();
      }
    }
    delete fwdA;
    delete cenA;
    
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
    if (fTriggerEff <= 0 || (1-fTriggerEff)<1e-6) {
      TParameter<double>* eff = 
	static_cast<TParameter<double>*>(results->FindObject("triggerEff"));
      if (eff) {
	fTriggerEff = eff->GetVal();
	fExtTriggerEff = true;
      }
      if (fTriggerEff <= 0) SetTriggerEfficiency(1);
    }

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
    if (fCentAxis) {
      TArrayD  bins(fCentAxis->GetNbins()+1);
      Int_t    nBins = 0;
      Double_t high  = -1;
      for (Int_t i = 1; i <= fCentAxis->GetNbins(); i++) {
	Double_t binLow  = fCentAxis->GetBinLowEdge(i);
	Double_t binHigh = fCentAxis->GetBinUpEdge(i);
	if (binLow  < fCentMin-.5) continue;
	if (binHigh > fCentMax+.5) continue;
	high = binHigh;
	bins[nBins] = binLow;
	nBins++;
      }
      bins[nBins] = high;
      fCentAxis->Set(nBins, bins.GetArray());
    }
	

    if (fOptions & kVerbose) {
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
	   "   Trigger:       %-30s  (0x%x)\n"
	   "   Efficiency:    %-6.4f\n"
	   "   sqrt(sNN):     %-30s  (%dGeV)\n"
	   "   System:        %-30s  (%d)\n"
	   "   Vz range:      %-30s  (%f,%f)\n"
	   "   Normalization: %-30s  (%d)\n"
	   "   Centrality:    %s\n"
	   "   Options:       %s",
	   fTrigString->GetTitle(), fTrigString->GetUniqueID(), 
	   fTriggerEff,
	   fSNNString->GetTitle(),  fSNNString->GetUniqueID(), 
	   fSysString->GetTitle(),  fSysString->GetUniqueID(), 
	   fVtxAxis->GetTitle(), fVtxAxis->GetXmin(), fVtxAxis->GetXmax(),
	   fNormString->GetTitle(), fNormString->GetUniqueID(),
	   centTxt.Data(), (options ? options->GetTitle() : "none"));
    }
    if (fSysString->GetUniqueID() == 3) {
      Info("FetchTopResults", "Left/Right assymmetry, mirror, and systematic "
	   "errors explicitly disabled for pPb");
      fOptions   &= ~kShowLeftRight;
      fOptions   &= ~kMirror;
      fFwdSysErr =  0;
      fCenSysErr =  0;
    }
  }
  //__________________________________________________________________
  TMultiGraph* FetchOthers(UShort_t centLow, UShort_t centHigh)
  {
    TMultiGraph* thisOther = 0;
    if (fShowOthers == 0) return 0;

    UShort_t sys   = (fSysString  ? fSysString->GetUniqueID() : 0);
    UShort_t trg   = (fTrigString ? fTrigString->GetUniqueID() : 0);
    UShort_t snn   = (fSNNString  ? fSNNString->GetUniqueID() : 0);
    Long_t   ret   = 
      gROOT->ProcessLine(Form("RefData::GetData(%d,%d,%d,%d,%d,%d);",
			      sys,snn,trg,centLow,centHigh,fShowOthers));
    if (!ret) {
      Warning("FetchOthers", 
	      "No other data for %s %s %s %3d%%-%3d%% central (0x%x)", 
	      fSysString  ? fSysString->GetTitle()  : "unknown", 
	      fTrigString ? fTrigString->GetTitle() : "unknown", 
	      fSNNString  ? fSNNString->GetTitle()  : "unknown", 
	      centLow, centHigh, fShowOthers);
      return 0;
    }

    thisOther = reinterpret_cast<TMultiGraph*>(ret);    
    return thisOther;
  }
  //__________________________________________________________________
  /** 
   * Get the results from the top-level list (MC, SPD, FMD)
   * 
   * @param list    List 
   * @param mcList  List of histograms from MC
   * @param empCorr Emperical correction if any
   * @param name    name 
   * @param max     On return, maximum of data 
   * @param rmax    On return, maximum of ratios
   * @param amax    On return, maximum of left-right comparisons
   * @param truths  List of MC truths to compare to. 
   *
   * @return Array of results
   */
  TObjArray* 
  FetchTopResults(const TList*  list, 
		  const TList*  mcList,
		  TGraphErrors* empCorr,
		  const char*   name, 
		  Double_t&     max,
		  Double_t&     rmax,
		  Double_t&     amax,
		  TObjArray&    truths)
  {
    if (!list) return 0;
    UShort_t   n = HasCent() ? fCentAxis->GetNbins() : 0;
    // Info("FetchTopResults","got %d centrality bins", n);
    if (n == 0) {
      TH1*  h = FetchOne(list, mcList, empCorr, name, "all",
			 FetchOthers(0,0), -1000, 0, 
			 max, rmax, amax, fTruth);
      if (!h) return 0;
      TObjArray* a = new TObjArray;
      // Info("FetchTopResults", "Adding %s to result stack", h->GetName());
      a->AddAt(h, 0);
      return a;
    }
    
    TObjArray* a = new TObjArray;
    truths.Expand(n);
    for (UShort_t i = 0; i < n; i++) { 
      UShort_t centLow  = fCentAxis->GetBinLowEdge(i+1);
      UShort_t centHigh = fCentAxis->GetBinUpEdge(i+1);
      TString  lname    = Form("cent%03d_%03d", centLow, centHigh);
      Int_t    col      = GetCentralityColor(i+1);
      TString  centTxt  = Form("%3d%%-%3d%% central", centLow, centHigh);

      TH1* tt = static_cast<TH1*>(truths.At(i));
      TH1* ot = tt;
      TH1* h  = FetchOne(list, mcList, empCorr, name, lname,
			 FetchOthers(centLow,centHigh), col, 
			 centTxt.Data(), max, rmax, amax, fTruth);
      if (!h) continue;
      if (ot != tt) { 
	//Info("FetchTopResults", "old truth=%p new truth=%p (%s)", 
	//     ot, tt, name);
	truths.AddAt(tt, i);
      }
      // Info("FetchTopResults", "Adding %p to result stack", h);
      a->AddAt(h, i);
    }
    return a;
  } 
  //__________________________________________________________________
  /** 
   * Steer retrieval one centrality bin results
   * 
   * @param list          Input list
   * @param mcList        Possible MC list
   * @param empCorr       Possible empirical correction
   * @param name          Name of bing 
   * @param folderName    What sub-folder to get 
   * @param others        What else to plot 
   * @param col           Color 
   * @param txt           Centrality text 
   * @param max           Current maximum, on return new maximum 
   * @param rmax          Current range maximum, on return new maximum 
   * @param amax          Current A maximum, on return new maximum 
   * @param truth         Possible MC truth histogram 
   * 
   * @return 
   */
  TH1* FetchOne(const TList*  list, 
		const TList*  mcList,
		TGraphErrors* empCorr,
		const char*   name, 
		const char*   folderName,
		TMultiGraph*  others, 
		Int_t         col,
		const char*   /* txt */,
		Double_t&     max,
		Double_t&     rmax,
		Double_t&     amax,
		TH1*&         truth)
  {
    TList* folder = (fOptions & kOldFormat ? const_cast<TList*>(list) :
		     static_cast<TList*>(list->FindObject(folderName)));
    if (!folder) {
      Error("FetchResults", "Couldn't find list '%s' in %s", 
	    folderName, list->GetName());
      return 0;
    }
    TList* mcFolder = 0;
    if (mcList) {
      mcFolder = static_cast<TList*>(mcList->FindObject(folderName));
      if (!mcFolder) 
	Warning("FetchResults", 
		"Didn't find the list '%s' in %s for final MC correction", 
		folderName, mcList->GetName());
    }
    if (fOptions & kVerbose) {
      TObject* normCalc = folder->FindObject("normCalc");
      if (normCalc) Info("FetchOne", "%s:\n%s", 
			 folderName, normCalc->GetTitle());
    }
    TH1* h = FetchCentResults(folder, mcFolder, empCorr, name, 
			      others, col, folderName, max, rmax, amax, truth);
    return h;
  }
  //__________________________________________________________________
  /** 
   * Fetch results for a particular centrality bin
   * 
   * @param list       List 
   * @param mcList     List of MC results
   * @param empCorr    Emperical correction if any 
   * @param name       Name 
   * @param thisOther  Other graphs 
   * @param color      Color 
   * @param centTxt    Centrality text
   * @param max        On return, data maximum
   * @param rmax       On return, ratio maximum 
   * @param amax       On return, left-right maximum 
   * @param truth      MC truth to compare to or possibly update
   *
   * @return Histogram of results 
   */
  TH1* FetchCentResults(const TList*  list, 
			const TList*  mcList, 
			TGraphErrors* empCorr,
			const char*   name, 
			TMultiGraph*  thisOther,
			Int_t         color,
			const char*   centTxt,
			Double_t&     max,
			Double_t&     rmax,
			Double_t&     amax, 
			TH1*&         truth)
  {
    
    TH1* dndeta      = FetchHistogram(list, Form("dndeta%s", name));
    TH1* dndetaMC    = FetchHistogram(list, Form("dndeta%sMC", name));
    TH1* dndetaTruth = FetchHistogram(list, "dndetaTruth");

    if (mcList && FetchHistogram(mcList, "finalMCCorr")) 
      Warning("FetchCentResults", "dNdeta already corrected for final MC");
    else 
      CorrectFinalMC(dndeta, mcList);
      
    CorrectEmpirical(dndeta, empCorr);
    CorrectTriggerEff(dndeta);
    CorrectTriggerEff(dndetaMC);

    TH1* dndetaSym   = 0;
    TH1* dndetaMCSym = 0;
    SetAttributes(dndeta,     color);
    SetAttributes(dndetaMC,   HasCent() ? color : color+2);
    SetAttributes(dndetaTruth,color);
    SetAttributes(dndetaSym,  color);
    SetAttributes(dndetaMCSym,HasCent() ? color : color+2);
    if (dndetaMC && HasCent()) 
      dndetaMC->SetMarkerStyle(dndetaMC->GetMarkerStyle()+2);
    if (dndetaMCSym && HasCent()) 
      dndetaMCSym->SetMarkerStyle(dndetaMCSym->GetMarkerStyle()+2);
    if (dndetaTruth && HasCent()) {
      dndetaTruth->SetMarkerStyle(34);
      dndetaTruth->SetMarkerColor(kYellow-1);
    }
    if (dndetaTruth) { 
      dndetaTruth->SetLineColor(kBlack); 
      dndetaTruth->SetFillColor(kBlack); 
      dndetaTruth->SetFillStyle(3002); 
      // dndetaTruth->SetLineColor(kBlack); 
    }
    ModifyTitle(dndeta,     centTxt);
    ModifyTitle(dndetaMC,   centTxt);
    ModifyTitle(dndetaTruth,centTxt);
    ModifyTitle(dndetaSym,  centTxt);
    ModifyTitle(dndetaMCSym,centTxt);


    max = TMath::Max(max, AddHistogram(fResults, dndetaTruth, "e5"));
    max = TMath::Max(max, AddHistogram(fResults, dndetaMC,    dndetaMCSym));
    max = TMath::Max(max, AddHistogram(fResults, dndeta,      dndetaSym));

    if (dndetaTruth) {
      truth = dndetaTruth;
    }
    else {
      if ((fOptions & kShowRings)) {
	THStack* rings = static_cast<THStack*>(list->FindObject("dndetaRings"));
	if (rings) { 
	  TIter next(rings->GetHists());
	  TH1*  hist = 0;
	  while ((hist = static_cast<TH1*>(next()))) 
	    max = TMath::Max(max, AddHistogram(fResults, hist));
	}
      }
      // Info("FetchCentResults", "Got %p, %p, %p from %s with name %s, max=%f",
      //      dndeta, dndetaMC, dndetaTruth, list->GetName(), name, max);
      
      if ((fOptions & kShowLeftRight)) {
	fLeftRight->Add(Asymmetry(dndeta,    amax));
	fLeftRight->Add(Asymmetry(dndetaMC,  amax));
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
    }
    if (dndetaMC) { 
      fRatios->Add(Ratio(dndeta,    dndetaMC,    rmax));
      fRatios->Add(Ratio(dndetaSym, dndetaMCSym, rmax));
    }
    if (truth) {
      fRatios->Add(Ratio(dndeta,      truth, rmax));
      fRatios->Add(Ratio(dndetaSym,   truth, rmax));
    }
    return dndeta;
  }
  //__________________________________________________________________
  void CorrectFinalMC(TH1* dndeta, const TList* mcList)
  {
    if (!dndeta) return;
    if (!mcList) return;

    TH1* dndetaMC    = FetchHistogram(mcList, dndeta->GetName());
    TH1* dndetaTruth = FetchHistogram(mcList, "dndetaTruth");
    if (!dndetaMC || !dndetaTruth) return;
    
    TH1* corr = static_cast<TH1*>(dndetaMC->Clone("finalMCCorr"));
    corr->Divide(dndetaTruth);
    
    Info("CorrectFinalMC", "Correcting dN/deta with final MC correction");
    dndeta->Divide(corr);
  }
  //__________________________________________________________________
  void CorrectEmpirical(TH1* dndeta, const TGraphErrors* empCorr) 
  {
    if (!dndeta) return;
    if (!empCorr) return;
   
    Info("CorrectEmpirical", "Doing empirical correction of dN/deta");
    TAxis* xAxis = dndeta->GetXaxis();
    for (Int_t i = 1; i <= xAxis->GetNbins(); i++) {
      Double_t x = xAxis->GetBinCenter(i);
      Double_t y = dndeta->GetBinContent(i);
      Double_t c = empCorr->Eval(x);
      dndeta->SetBinContent(i, y / c);
    }
  }
  //__________________________________________________________________
  void CorrectTriggerEff(TH1* dndeta)
  {
    if (fExtTriggerEff) return;
    if (!dndeta) return;
    if (fTriggerEff <= 0 || fTriggerEff >= 1) return;
    dndeta->Scale(fTriggerEff);
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
    gStyle->SetTitleFont(kFont, "xyz");
    gStyle->SetLabelFont(kFont, "xyz");
    
    Int_t    h  = (gROOT->IsBatch() ? 
		   ((fOptions & kHiRes) ? 1600 : 900) :
		   gClient->GetDisplayHeight());
    Int_t    w  = h; // h / TMath::Sqrt(2);
    Double_t y1 = 0;
    Double_t y2 = 0;
    Double_t y3 = 0;
    if (!(fOptions & kShowRatios))    w  *= 1.3;
    else                 y1 =  0.3;
    if (!(fOptions & kShowLeftRight)) w  *= 1.3;
    else { 
      Double_t y11 = y1;
      y1 = (y11 > 0.0001 ? 0.4 : 0.2);
      y2 = (y11 > 0.0001 ? 0.2 : 0.3);
    }
    TCanvas* c = new TCanvas("Results", "Results", w, h);
    c->SetFillColor(0);
    c->SetBorderSize(0);
    c->SetBorderMode(0);

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
    if (HasCent()) trg = "CENT";
    trg          = trg.Strip(TString::kBoth);
    trg.ReplaceAll(" ", "_");
    trg.ReplaceAll(">", "Gt");
    trg.ReplaceAll("&", "AND");
    trg.ReplaceAll("|", "OR");
    if (fBase.IsNull()) 
      fBase = "dndeta_<sys>_<snn>_<trig>_<ipmin><ipmax>cm_<nev>ev";
    fBase.ReplaceAll("<sys>",   fSysString->GetTitle());
    fBase.ReplaceAll("<snn>",   fSNNString->GetTitle());
    fBase.ReplaceAll("<trig>",  trg.Data());
    fBase.ReplaceAll("<ipmin>", Form("%c%02d",vMin<0?'m':'p',TMath::Abs(vMin)));
    fBase.ReplaceAll("<ipmax>", Form("%c%02d",vMax<0?'m':'p',TMath::Abs(vMax)));
    fBase.ReplaceAll("<nev>",   Form("%09d",  nev));
    if ((fFormats & kPNG))   c->SaveAs(Form("%s.png",  fBase.Data()));
    if ((fFormats & kROOT))  c->SaveAs(Form("%s.root", fBase.Data()));
    if ((fFormats & kScript))c->SaveAs(Form("%s.C",    fBase.Data()));
    if ((fFormats & kPDF))   c->SaveAs(Form("%s.pdf",  fBase.Data()));
    if (fOptions & kExport) {
      TString exp(fBase);
      exp.ReplaceAll("dndeta", "export");
      exp.ReplaceAll("dNdeta", "export");
      Export(exp);
    }
  }
  //__________________________________________________________________
  /** 
   * Build main legend 
   * 
   * @param stack    Stack to include 
   * @param mg       (optional) Multi graph to include 
   * @param x1       Lower X coordinate in the range [0,1]
   * @param y1       Lower Y coordinate in the range [0,1]
   * @param x2       Upper X coordinate in the range [0,1]
   * @param y2 	     Upper Y coordinate in the range [0,1]
   * @param forceCol If non-zero, force this many columns
   */
  void BuildLegend(THStack* stack, TMultiGraph* mg, 
		   Double_t x1, Double_t y1, Double_t x2, Double_t y2,
		   Int_t forceCol=0)
  {
    TLegend* l = new TLegend(x1,y1,x2,y2);
    Int_t nCol = forceCol;
    if (nCol <= 0) nCol = HasCent() ? 1 : 2;
    l->SetNColumns(nCol);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextFont(kFont);
    l->SetTextColor(kAliceBlue);

    // Loop over items in stack and get unique items, while ignoring
    // mirrored data and systematic error bands 
    TIter    next(stack->GetHists());
    TH1*     hist = 0;
    TObjArray unique;
    unique.SetOwner();
    Bool_t   sysErrSeen = false;
    while ((hist = static_cast<TH1*>(next()))) { 
      TString t(hist->GetTitle());
      TString n(hist->GetName());
      n.ToLower();
      if (t.Contains("mirrored")) continue;
      if (n.Contains("syserror")) { sysErrSeen = true; continue; }
      if (unique.FindObject(t.Data())) continue;
      TObjString* s1 = new TObjString(hist->GetTitle());
      s1->SetUniqueID(((hist->GetMarkerStyle() & 0xFFFF) << 16) |
		      ((hist->GetMarkerColor() & 0xFFFF) <<  0));
      unique.Add(s1);
      // l->AddEntry(hist, hist->GetTitle(), "pl");
    }
    if (mg) {
      // If we have other stuff, scan for unique names 
      TIter nexto(mg->GetListOfGraphs());
      TGraph* g = 0;
      while ((g = static_cast<TGraph*>(nexto()))) { 
	TString n(g->GetTitle());
	if (n.Contains("mirrored")) continue;
	if (unique.FindObject(n.Data())) continue;
	TObjString* s2 = new TObjString(n);
	s2->SetUniqueID(((g->GetMarkerStyle() & 0xFFFF) << 16) |
			 ((g->GetMarkerColor() & 0xFFFF) <<  0));
	unique.Add(s2);
	// l->AddEntry(hist, hist->GetTitle(), "pl");
      }
    }

    // Add legend entries for unique items only
    TIter nextu(&unique);
    TObject* s = 0;
    Int_t i = 0;
    while ((s = nextu())) { 
      TLegendEntry* dd = l->AddEntry(Form("data%2d", i++), 
				     s->GetName(), "lp");
      Int_t style = (s->GetUniqueID() >> 16) & 0xFFFF;
      Int_t color = (s->GetUniqueID() >>  0) & 0xFFFF;
      dd->SetLineColor(kBlack);
      if (HasCent()) dd->SetMarkerColor(kBlack);
      else           dd->SetMarkerColor(color);
      dd->SetMarkerStyle(style);
    }
    if (sysErrSeen) {
      // Add entry for systematic errors 
      TLegendEntry* d0 = l->AddEntry("d0", Form("%4.1f%% Systematic error", 
						100*fFwdSysErr), "f");
      d0->SetLineColor(kSysErrColor);
      d0->SetMarkerColor(kSysErrColor);
      d0->SetFillColor(kSysErrColor);
      d0->SetFillStyle(SYSERR_STYLE);
      d0->SetMarkerStyle(0);
      d0->SetLineWidth(0);
      i++;
    }
    if (nCol == 2 && i % 2 == 1)  {
      // To make sure the 'data' and 'mirrored' entries are on a line
      // by themselves 
      TLegendEntry* dd = l->AddEntry("dd", "   ", "");
      dd->SetTextSize(0);
      dd->SetFillStyle(0);
      dd->SetFillColor(0);
      dd->SetLineWidth(0);
      dd->SetLineColor(0);
      dd->SetMarkerSize(0);
    }
    if ((fOptions & kMirror)) {
      // Add entry for 'data'
      TLegendEntry* d1 = l->AddEntry("d1", "Data", "lp");
      d1->SetLineColor(kBlack);
      d1->SetMarkerColor(kBlack);
      d1->SetMarkerStyle(20);

      // Add entry for 'mirrored data'
      TLegendEntry* d2 = l->AddEntry("d2", "Mirrored data", "lp");
      d2->SetLineColor(kBlack);
      d2->SetMarkerColor(kBlack);
      d2->SetMarkerStyle(24);
    }
    
    l->Draw();
  }
  //__________________________________________________________________
  /** 
   * Build centrality legend 
   * 
   * @param x1      Lower X coordinate in the range [0,1]
   * @param y1      Lower Y coordinate in the range [0,1]
   * @param x2      Upper X coordinate in the range [0,1]
   * @param y2 	    Upper Y coordinate in the range [0,1]
   */
  void BuildCentLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
  {
    if (!HasCent()) return;

    TLegend* l = new TLegend(x1,y1,x2,y2);
    l->SetNColumns(1);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextFont(kFont);
    l->SetTextColor(kAliceBlue);

    Int_t n = fCentAxis->GetNbins();
    for (Int_t i = 1; i <= n; i++) { 
      Double_t low = fCentAxis->GetBinLowEdge(i);
      Double_t upp = fCentAxis->GetBinUpEdge(i);
      TLegendEntry* e = l->AddEntry(Form("dummy%02d", i),
				    Form("%3d%% - %3d%%", 
					 int(low), int(upp)), "pl");
      e->SetMarkerColor(GetCentralityColor(i));
    }
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
    p1->SetRightMargin(kRightMargin);
    if ((fOptions & kShowLeftRight) || (fOptions & kShowRatios)) p1->SetGridx();
    p1->SetTicks(1,1);
    p1->SetNumber(1);
    p1->Draw();
    p1->cd();

    // Info("PlotResults", "Plotting results with max=%f", max);
    fResults->SetMaximum((fOptions & kExtraWhite ? 1.4 : 1.15)*max);
    fResults->SetMinimum(yd > 0.00001 ? -0.02*max : 0);
    // fResults->SetMinimum(yd > 0.00001 ? -0.02*max : 0);

    FixAxis(fResults, (1-yd)*(yd > .001 ? 1 : .9 / 1.2), 
	    "#frac{1}{#it{N}}#kern[.1]{#frac{d#it{N}_{ch}}{d#it{#eta}}}");

    p1->Clear();
    fResults->DrawClone("nostack e1");

    fRangeParam->fSlave1Axis = fResults->GetXaxis();
    fRangeParam->fSlave1Pad  = p1;

    // Draw other data
    if (fShowOthers != 0) {
      TGraphAsymmErrors* o      = 0;
      TIter              nextG(fOthers->GetListOfGraphs());
      while ((o = static_cast<TGraphAsymmErrors*>(nextG())))
        o->DrawClone("same p");
    }

    // Make a legend in the main result pad
    BuildCentLegend(.12, 1-p1->GetTopMargin()-.01-.5,  
		    .35, 1-p1->GetTopMargin()-.01-.1);
    Double_t x1 = .15;
    Double_t x2 = .90;
    Double_t y1 = p1->GetBottomMargin()+.01;
    Double_t y2 = .35;
    if (HasCent()) { 
      if (fOptions & kExtraWhite) { 
	x1 = .40;
	x2 = .60;
	y1 = .70;
	y2 = 1-p1->GetTopMargin()-0.06;
      }
      else {
	x1 = .70;
	x2 = 1-p1->GetRightMargin()-.01;
	y1 = .50;
	y2 = 1-p1->GetTopMargin()-0.16;
      }
    }
		   
    BuildLegend(fResults, fOthers, x1, y1, x2, y2);

    // Put a title on top
    fTitle.ReplaceAll("@", " ");
    TLatex* tit = new TLatex(0.10, .99, fTitle.Data());
    tit->SetNDC();
    tit->SetTextFont(kFont);
    tit->SetTextAlign(13);
    tit->SetTextSize(0.045);
    tit->SetTextColor(kAlicePurple);
    tit->Draw();

    // TColor::GetColor(41,73,156);
    Double_t x         = .95;
    Double_t y         = .93;
    // Put a nice label in the plot
    TString     eS;
    UShort_t    snn = fSNNString->GetUniqueID();
    const char* sys = fSysString->GetTitle();
    if (snn == 2750) snn = 2760;
    if (snn > 1000) eS = Form("%4.2fTeV", float(snn)/1000);
    else            eS = Form("%3dGeV", snn);
    TLatex* tt = new TLatex(x, y, Form("%s #sqrt{s%s}=%s, %s", 
				       sys, 
				       (HasCent() ? "_{NN}" : ""),
				       eS.Data(), 
				       HasCent() ? "by centrality" : 
				       fTrigString->GetTitle()));
    tt->SetTextColor(kAliceBlue);
    tt->SetNDC();
    tt->SetTextFont(kFont);
    tt->SetTextAlign(33);
    tt->Draw();
    y -= tt->GetTextSize() + .01;
    
    // Put number of accepted events on the plot
    Int_t nev = 0;
    if (fTriggers) nev = fTriggers->GetBinContent(1);
    TLatex* et = new TLatex(x, y, Form("%d events", nev));
    et->SetTextColor(kAliceBlue);
    et->SetNDC();
    et->SetTextFont(kFont);
    et->SetTextAlign(33);
    et->Draw();
    y -= et->GetTextSize() + .01;

    // Put number of accepted events on the plot
    if (fVtxAxis) { 
      TLatex* vt = new TLatex(x, y, fVtxAxis->GetTitle());
      vt->SetNDC();
      vt->SetTextFont(kFont);
      vt->SetTextAlign(33);
      vt->SetTextColor(kAliceBlue);
      vt->Draw();
      y -= vt->GetTextSize() + .01;
    }
    // results->Draw("nostack e1 same");

    TString corrs;
    if (!fEmpirical.IsNull()) corrs.Append("Emperical");
    if (!fFinalMC.IsNull())   {
      if (!corrs.IsNull()) corrs.Append("+");
      corrs.Append("Final MC");
    }

    if (!corrs.IsNull()) {
      corrs.Append(" correction");
      if (corrs.Index("+") != kNPOS) corrs.Append("s");
      TLatex* em = new TLatex(x, y, corrs);
      em->SetNDC();
      em->SetTextFont(kFont);
      em->SetTextAlign(33);
      em->SetTextColor(kAliceBlue);
      em->Draw();
      y -= em->GetTextSize() + .01;
    }
      
    if (fTriggerEff > 0 && fTriggerEff <= 1 && !HasCent()) { 
      TLatex* ef = new TLatex(x, y, Form("#varepsilon_{%s} = %5.3f", 
					 fTrigString->GetTitle(), 
					 fTriggerEff));
      ef->SetNDC();
      ef->SetTextFont(kFont);
      ef->SetTextAlign(33);
      ef->SetTextColor(kAliceBlue);
      ef->Draw();
      y -= ef->GetTextSize() + .01;
    }
    
    fRangeParam->fSlave1Axis = FindXAxis(p1, fResults->GetName());
    fRangeParam->fSlave1Pad  = p1;


    // Mark the plot as preliminary
    TLatex* pt = new TLatex(.12, .93, "Work in progress");
    pt->SetNDC();
    pt->SetTextFont(62);
    // pt->SetTextSize();
    pt->SetTextColor(kAliceRed);
    pt->SetTextAlign(13);
    pt->Draw();
    TDatime now;
    TLatex* dt = new TLatex(.12, .88, now.AsSQLString());
    dt->SetNDC();
    dt->SetTextFont(42);
    dt->SetTextSize(0.04);
    dt->SetTextColor(kAliceBlue); // kAliceRed);
    dt->SetTextAlign(13);
    dt->Draw();

    TString savPath(gROOT->GetMacroPath());
    gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWGLF/FORWARD/analysis2/scripts",
			     gROOT->GetMacroPath()));
    // Always recompile 
    if (!gROOT->GetClass("AliceLogo"))
      gROOT->LoadMacro("AliceLogo.C++");
    gROOT->SetMacroPath(savPath);

    if (gROOT->GetClass("AliceLogo")) {
      p1->cd();
      p1->Range(0,0,1,1);
      gROOT->ProcessLine("AliceLogo* al = new AliceLogo();");
      gROOT->ProcessLine("al->Draw(0,.88,.5,.2, 0, 0);");
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
    if ((fOptions & kShowRatios) == 0) return;

    bool isBottom = (y1 < 0.0001);
    Double_t yd = y2 - y1;
    // Make a sub-pad for the result itself
    TPad* p2 = new TPad("p2", "p2", 0, y1, 1.0, y2, 0, 0, 0);
    p2->SetTopMargin(0.001);
    p2->SetRightMargin(kRightMargin);
    p2->SetBottomMargin(isBottom ? 1/yd * 0.07 : 0.0001);
    p2->SetGridx();
    p2->SetTicks(1,1);
    p2->SetNumber(2);
    p2->Draw();
    p2->cd();

    // Fix up axis
    FixAxis(fRatios, yd, "Ratios", 7);

    fRatios->SetMaximum(1+TMath::Max(.22,1.05*max));
    fRatios->SetMinimum(1-TMath::Max(.32,1.05*max));
    p2->Clear();
    fRatios->DrawClone("nostack e1");

    
    // Make a legend
    BuildLegend(fRatios, 0, .15,p2->GetBottomMargin()+.01,.9,
		isBottom ? .6 : .4, 2);
#if 0
    TLegend* l2 = p2->BuildLegend(.15,p2->GetBottomMargin()+.01,.9,
				  isBottom ? .6 : .4);
    l2->SetNColumns(2);
    l2->SetFillColor(0);
    l2->SetFillStyle(0);
    l2->SetBorderSize(0);
    l2->SetTextFont(kFont);
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

    if ((fOptions & kAddExec)) {
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
    if (!(fOptions & kShowLeftRight)) return;

    bool isBottom = (y1 < 0.0001);
    Double_t yd = y2 - y1;
    // Make a sub-pad for the result itself
    TPad* p3 = new TPad("p3", "p3", 0, y1, 1.0, y2, 0, 0, 0);
    p3->SetTopMargin(0.001);
    p3->SetRightMargin(kRightMargin);
    p3->SetBottomMargin(isBottom ? 1/yd * 0.07 : 0.0001);
    p3->SetGridx();
    p3->SetTicks(1,1);
    p3->SetNumber(2);
    p3->Draw();
    p3->cd();

    // Fix up axis
    FixAxis(fLeftRight, yd, "Right/Left", 4);

    fLeftRight->SetMaximum(1+TMath::Max(.12,1.05*max));
    fLeftRight->SetMinimum(1-TMath::Max(.15,1.05*max));
    p3->Clear();
    fLeftRight->DrawClone("nostack e1");

    
    // Make a legend
    Double_t xx1 = (HasCent() ? .7                           : .15); 
    Double_t xx2 = (HasCent() ? 1-p3->GetRightMargin()-.01   : .90);
    Double_t yy1 = p3->GetBottomMargin()+.01;
    Double_t yy2 = (HasCent() ? 1-p3->GetTopMargin()-.01-.15 : .5);
    BuildLegend(fLeftRight, 0, xx1, yy1, xx2, yy2);
    // TLegend* l2 = p3->BuildLegend(.15,p3->GetBottomMargin()+.01,.9,.5);
    // l2->SetNColumns(2);
    // l2->SetFillColor(0);
    // l2->SetFillStyle(0);
    // l2->SetBorderSize(0);
    // l2->SetTextFont(132);

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
  //__________________________________________________________________
  /** 
   * Get the color for a centrality bin
   * 
   * @param bin Centrality bin 
   * 
   * @return Color 
   */
  Int_t GetCentralityColor(Int_t bin) const
  {
    if (fCentAxis->GetNbins() < 6) { 
      switch (bin) { 
      case 1: return kRed+2;
      case 2: return kGreen+2;
      case 3: return kBlue+1;
      case 4: return kCyan+1;
      case 5: return kMagenta+1;
      case 6: return kYellow+2;
      }
    }
    UShort_t centLow  = fCentAxis->GetBinLowEdge(bin);
    UShort_t centHigh = fCentAxis->GetBinUpEdge(bin);
    Float_t  fc       = (centLow+double(centHigh-centLow)/2) / 100;
    Int_t    nCol     = gStyle->GetNumberOfColors();
    Int_t    icol     = TMath::Min(nCol-1,int(fc * nCol + .5));
    Int_t    col      = gStyle->GetColorPalette(icol);
    //Info("GetCentralityColor","%3d: %3d-%3d -> %3d",bin,centLow,centHigh,col);
    return col;
  }
  //__________________________________________________________________
  /** 
   * Set attributed on a histogram 
   * 
   * @param h     Histogram
   * @param color Color 
   */
  void SetAttributes(TH1* h, Int_t color)
  {
    if (!h) return;
    if (color < 0) return;
    // h->SetLineColor(color);
    h->SetMarkerColor(color);
    // h->SetFillColor(color);
  }
  //__________________________________________________________________
  /** 
   * Set attributed on a graph 
   * 
   * @param g     Graph
   * @param color Color 
   */
  void SetAttributes(TGraph* g, Int_t color)
  {
    if (!g) return;
    if (color < 0) return;
    // g->SetLineColor(color);
    g->SetMarkerColor(color);
    // g->SetFillColor(color);
  }
  //__________________________________________________________________
  /** 
   * Modify the title 
   * 
   */
  void ModifyTitle(TNamed* h, const char* centTxt)
  {
    if (!h) return;

    TString title(h->GetTitle());
    title.ReplaceAll("ALICE ","");
    if (title.Contains("Central")) 
      title.ReplaceAll("CentraldNdeta", "SPD clusters");
    if (title.Contains("Forward"))
      title.ReplaceAll("ForwarddNdeta", "FMD");
    h->SetTitle(title);

    if (centTxt && centTxt[0] != '\0') {
      TString name(h->GetName());
      name.Append(Form("_%s", centTxt));
      h->SetName(name);
    }
    
    return;
    // if (!centTxt || !h) return;
    // h->SetTitle(Form("%s, %s", h->GetTitle(), centTxt));
  }
  /** 
   * Get a result from the passed list
   * 
   * @param list List to search 
   * @param name Object name to search for 
   * 
   * @return Histogram
   */
  TH1* FetchHistogram(const TList* list, const char* name) const 
  {
    if (!list) return 0;
    
    TH1* ret = static_cast<TH1*>(list->FindObject(name));
#if 0
    if (!ret) {
      // all->ls();
      Warning("GetResult", "Histogram %s not found", name);
    }
#endif

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

    stack->Add(hist, option);
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
    if ((fOptions & kMirror)) {
      sym = Symmetrice(hist);
      stack->Add(sym, option);
    }

    return hist->GetMaximum();
  }

  //__________________________________________________________________
  /** 
   * Rebin a histogram 
   * 
   * @param h     Histogram to rebin
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

	if ((fOptions & kCutEdges)) {
	  if (h->GetBinContent(bin+1)<=0 || 
	      h->GetBinContent(bin-1)<=0) {
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
   * @param g Graph
   * 
   * @return Histogram
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

    TH1*        r  = 0;
    const TAttMarker* m1 = 0;
    const TAttMarker* m2 = 0;
    const TH1* h1 = dynamic_cast<const TH1*>(o1); 
    if (h1) { 
      m1 = h1;
      const TH1* h2 = dynamic_cast<const TH1*>(o2); 
      if (h2) { 
	m2 = h2;
	r = RatioHH(h1,h2);
      }
      else {
	const TGraph* g2 = dynamic_cast<const TGraph*>(o2);
	if (g2) { 
	  m2 = g2;
	  r = RatioHG(h1,g2);      
	}
      }
    }
    else {
      const TGraphAsymmErrors* g1 = dynamic_cast<const TGraphAsymmErrors*>(o1);
      if (g1) { 
	m1 = g1;
	const TGraphAsymmErrors* g2 = 
	  dynamic_cast<const TGraphAsymmErrors*>(o2);
	if (g2) {
	  m2 = g2;
	  r = RatioGG(g1, g2);
	}
      }
    }
    if (!r) {
      // Warning("Ratio", "Don't know how to divide a %s (%s) with a %s (%s)", 
      //         o1->ClassName(),o1->GetName(),o2->ClassName(),o2->GetName());
      return 0;
    }
    // Check that the histogram isn't empty
    if (r->GetEntries() <= 0) { 
      delete r; 
      r = 0; 
    }
    if (r) {
      r->SetMarkerStyle(m2->GetMarkerStyle());
      r->SetMarkerColor(m1->GetMarkerColor());
      if (TString(o2->GetName()).Contains("truth", TString::kIgnoreCase)) 
	r->SetMarkerStyle(m1->GetMarkerStyle());
      r->SetMarkerSize(0.9*m1->GetMarkerSize());
      r->SetName(Form("%s_over_%s", o1->GetName(), o2->GetName()));
      r->SetTitle(Form("%s / %s", o1->GetTitle(), o2->GetTitle()));
      r->SetDirectory(0);
      max = TMath::Max(RatioMax(r), max);
    }

    return r;
  }
  //__________________________________________________________________
  /** 
   * Compute the ratio of @a h to @a g.  @a g is evaluated at the bin
   * centers of @a h 
   * 
   * @param h  Numerator 
   * @param g  Divisor 
   * 
   * @return h/g 
   */
  TH1* RatioHG(const TH1* h, const TGraph* g) const 
  {
    if (!h || !g) return 0;

    TH1* ret = static_cast<TH1*>(h->Clone("tmp"));
    ret->Reset();
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
    return ret;
  }
  //__________________________________________________________________
  /** 
   * Make the ratio of h1 to h2 
   * 
   * @param h1 First histogram (numerator) 
   * @param h2 Second histogram (denominator)
   * 
   * @return h1 / h2
   */
  TH1* RatioHH(const TH1* h1, const TH1* h2) const
  {
    if (!h1 || !h2) return 0;
    Bool_t bad = false;
    if (h1->GetNbinsX() != h2->GetNbinsX()) {
      Error("RatioHH", "They have differnet number of bins");
      bad = true;
    }
    for (Int_t i = 1; i <= h1->GetNbinsX(); i++) {
      if (h1->GetXaxis()->GetBinLowEdge(i) != 
	  h2->GetXaxis()->GetBinLowEdge(i)) {
	// Error("RatioHH", "They have incompatible variable bins");
	bad = true;
	break;
      }
    }
    if (bad) return 0;
    
    TH1* t1 = static_cast<TH1*>(h1->Clone("tmp"));
    // Printf("Dividing %s with %s", h1->GetName(), h2->GetName());
    t1->Divide(h2);
    return t1;
  }
  //__________________________________________________________________
  /** 
   * Calculate the ratio of two graphs - g1 / g2
   * 
   * @param g1 Numerator 
   * @param g2 Denominator
   * 
   * @return g1 / g2 in a histogram 
   */
  TH1* RatioGG(const TGraphAsymmErrors* g1, 
	       const TGraphAsymmErrors* g2) const
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
    TH1* h = 0;
    if (dx != 0) {
      h = new TH1F("tmp", "tmp", nBins, bins[0], bins[nBins]);
    }
    else {
      h = new TH1F("tmp", "tmp", nBins, bins.fArray);
    }

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
   * @param yd     How the canvas is cut
   * @param ytitle Y axis title
   * @param ynDiv  Divisions on Y axis
   * @param force  Whether to draw the stack first or not
   */
  void FixAxis(THStack* stack, Double_t yd, const char* ytitle,
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
    Double_t s = 1/yd/1.2;
    // Info("FixAxis", "for %s, s=1/%f=%f", stack->GetName(), yd, s);

    h->SetXTitle("#it{#eta}");
    h->SetYTitle(ytitle);
    TAxis* xa = h->GetXaxis();
    TAxis* ya = h->GetYaxis();

    // Int_t   npixels = 20
    // Float_t dy = gPad->AbsPixeltoY(0) - gPad->AbsPixeltoY(npixels);
    // Float_t ts = dy/(gPad->GetY2() - gPad->GetY1());

    if (xa) {
      // xa->SetTitle(h->GetXTitle());
      // xa->SetTicks("+-");
      xa->SetTitleSize(s*xa->GetTitleSize());
      xa->SetLabelSize(s*xa->GetLabelSize());
      xa->SetTickLength(s*xa->GetTickLength());
      // xa->SetTitleOffset(xa->GetTitleOffset()/s);

      if (stack != fResults) {
	TAxis* rxa = fResults->GetXaxis();
	xa->Set(rxa->GetNbins(), rxa->GetXmin(), rxa->GetXmax());
      }
    }
    if (ya) {
      // ya->SetTitle(h->GetYTitle());
      ya->SetDecimals();
      // ya->SetTicks("+-");
      ya->SetNdivisions(ynDiv);
      ya->SetTitleSize(s*ya->GetTitleSize());
      ya->SetTitleOffset(1.15*ya->GetTitleOffset()/s);
      ya->SetLabelSize(s*ya->GetLabelSize());
    }
  }
  //__________________________________________________________________
  /** 
   * Merge two histograms into one 
   * 
   * @param cen    Central part
   * @param fwd    Forward part
   * @param xlow   On return, lower eta bound
   * @param xhigh  On return, upper eta bound
   * 
   * @return Newly allocated histogram or null
   */
  TH1* 
  Merge(const TH1* cen, const TH1* fwd, Double_t& xlow, Double_t& xhigh)
  {
    TH1* tmp = static_cast<TH1*>(fwd->Clone("tmp"));
    TString name(fwd->GetName());
    name.ReplaceAll("Forward", "Merged");
    tmp->SetName(name);

    // tmp->SetMarkerStyle(28);
    // tmp->SetMarkerColor(kBlack);
    tmp->SetDirectory(0);
    xlow  = 100;
    xhigh = -100;
    for (Int_t i = 1; i <= tmp->GetNbinsX(); i++) {
      Double_t cc = cen->GetBinContent(i);
      Double_t cf = fwd->GetBinContent(i);
      Double_t ec = cen->GetBinError(i);
      Double_t ef = fwd->GetBinError(i);
      Double_t nc = cf;
      Double_t ne = ef;
      if (cc < 0.001 && cf < 0.01) continue;
      xlow  = TMath::Min(tmp->GetXaxis()->GetBinLowEdge(i),xlow);
      xhigh = TMath::Max(tmp->GetXaxis()->GetBinUpEdge(i),xhigh);
      if (cc > 0.001) {
	nc = cc;
	ne = ec;
	if (cf > 0.001) {
	  nc  = (cf + cc) / 2;
	  ne  = TMath::Sqrt(ec*ec + ef*ef);
	}
      }
      tmp->SetBinContent(i, nc);
      tmp->SetBinError(i, ne);
    }
    return tmp;
  }
  //____________________________________________________________________
  /** 
   * Fit  @f$g(x;A_1,A_2,\sigma_1,\sigma_2)@f$ to histogram data 
   * 
   * @param tmp    Histogram
   * @param xlow   Lower x bound
   * @param xhigh  Upper x bound 
   *
   * @return Fitted function 
   */
  TF1* 
  FitMerged(TH1* tmp, Double_t xlow, Double_t xhigh)
  {
    TF1* tmpf  = new TF1("tmpf", "gaus", xlow, xhigh);
    tmp->Fit(tmpf, "NQ", "");
    tmp->SetDirectory(0);

    TF1* fit = new TF1("f", myFunc, xlow, xhigh, 4);
    fit->SetParNames("a_{1}", "a_{2}", "#sigma_{1}", "#sigma_{2}");
    fit->SetParameters(tmpf->GetParameter(0), 
		       .2, 
		       tmpf->GetParameter(2), 
		       tmpf->GetParameter(2)/4);
    fit->SetParLimits(3, 0, 100);
    fit->SetParLimits(4, 0, 100);
    const char* fitOpts = (fOptions & kVerbose ? "0W" : "Q0W");
    tmp->Fit(fit,fitOpts,"");

    delete tmpf;
    return fit;
  }
  //____________________________________________________________________
  /** 
   * Make band of systematic errors 
   * 
   * @param tmp Histogram
   * @param cen Central 
   * @param fwd Forward 
   * @param fit Fit 
   */
  void
  MakeSysError(TH1* tmp, TH1* cen, TH1* fwd, TF1* fit)
  {
    for (Int_t i = 1; i <= tmp->GetNbinsX(); i++) {
      Double_t tc = tmp->GetBinContent(i);
      if (tc < 0.01) continue;
      Double_t fc = fwd->GetBinContent(i);
      Double_t cc = cen->GetBinContent(i);
      Double_t sysErr = fFwdSysErr;
      if (cc > .01 && fc > 0.01) 
	sysErr = (fFwdSysErr+fCenSysErr) / 2;
      else if (cc > .01) 
	sysErr = fCenSysErr;
      Double_t x = tmp->GetXaxis()->GetBinCenter(i);
      Double_t y = fit->Eval(x);
      tmp->SetBinContent(i, y);
      tmp->SetBinError(i,sysErr*y);
    }
    TString name(tmp->GetName());
    name.ReplaceAll("Merged", "SysError");
    tmp->SetName(name);
    tmp->SetMarkerColor(kSysErrColor);
    tmp->SetLineColor(kSysErrColor);
    tmp->SetFillColor(kSysErrColor);
    tmp->SetFillStyle(SYSERR_STYLE);
    tmp->SetMarkerStyle(0);
    tmp->SetLineWidth(0);
  }
  void CorrectForward(TH1* h) const
  {
    if (!(fOptions & kRemoveOuters)) return;
    
    for (Int_t i = 1; i <= h->GetNbinsX(); i++) { 
      Double_t eta = h->GetBinCenter(i);
      if (TMath::Abs(eta) < 2.3) { 
	h->SetBinContent(i, 0);
	h->SetBinError(i, 0);
      }
    }
  }
  void CorrectCentral(TH1* h) const 
  {
    if (fClusterScale.IsNull()) return;
    TString t(h->GetTitle());
    Info("CorrectCentral", "Replacing Central with Tracklets in %s", t.Data());
    t.ReplaceAll("Central", "Tracklets");
    h->SetTitle(t);

    TF1* cf = new TF1("clusterScale", fClusterScale);
    for (Int_t i = 1; i <= h->GetNbinsX(); i++) { 
      Double_t eta = h->GetBinCenter(i);
      Double_t f   = cf->Eval(eta);
      Double_t c   = h->GetBinContent(i);
      if (f < .1) f = 1;
      h->SetBinContent(i, c / f);
    }
    delete cf;
  }
  //____________________________________________________________________
  void Export(const char* basename)
  {
    TString bname(basename);
    bname.ReplaceAll(" ", "_");
    bname.ReplaceAll("-", "_");
    TString fname(Form("%s.C", bname.Data()));

    std::ofstream outf(fname.Data());
    if (!outf) { 
      Error("Export", "Failed to open output file %s", fname.Data());
      return;
    }
    Info("Export", "Exporting data to %s", fname.Data());
    outf << "// Create by dNdetaDrawer\n"
	 << "void " << bname << "(THStack* stack, TLegend* l, Int_t m)\n"
	 << "{"
	 << "   Int_t ma[] = { 24, 25, 26, 32,\n"
	 << "                  20, 21, 22, 33,\n"
	 << "                  34, 30, 29, 0, \n"
	 << "                  23, 27 };\n"
	 << "   Int_t mm = ((m < 20 || m > 34) ? 0 : ma[m-20]);\n\n";
    TList* hists = fResults->GetHists();
    TIter  next(hists);
    TH1*   hist = 0;
    while ((hist = static_cast<TH1*>(next()))) { 
      TString hname = hist->GetName();
      hname.Append(Form("_%04x", (gRandom->Integer(0xffff) & 0xffff)));
      hist->SetName(hname);
      hist->GetListOfFunctions()->Clear();
      hist->SavePrimitive(outf, "nodraw");
      bool mirror = hname.Contains("mirror");
      bool syserr = hname.Contains("SysError");
      if (!syserr) 
	outf << "   " << hname << "->SetMarkerStyle(" 
	     << (mirror ? "mm" : "m") << ");\n";
      else 
	outf << "   " << hname << "->SetMarkerStyle(1);\n";
      outf << "   stack->Add(" << hname 
	   << (syserr ? ",\"e5\"" : "") << ");\n\n";
    }
    UShort_t    snn = fSNNString->GetUniqueID();
    // const char* sys = fSysString->GetTitle();
    TString eS;
    if      (snn == 2750)     snn = 2760;
    if      (snn < 1000)      eS = Form("%3dGeV", snn);
    else if (snn % 1000 == 0) eS = Form("%dTeV", snn/1000);
    else                      eS = Form("%4.2fTeV", float(snn)/1000);
    outf << "  if (l) {\n"
	 << "    TLegendEntry* e = l->AddEntry(\"\",\"" << eS << "\",\"pl\");\n"
	 << "    e->SetMarkerStyle(m);\n"
	 << "    e->SetMarkerColor(kBlack);\n"
	 << "  }\n"
	 << "}\n" << std::endl;
  }
  /* @} */ 
  /** 
   * Check if we have centrality dependent information, and we're not
   * forcing to use minimum bias
   * 
   * 
   * @return True if we should do centrality dependent ploting 
   */
  Bool_t HasCent() const 
  { 
    return fCentAxis && fCentAxis->GetNbins() > 0 && !(fOptions & kForceMB); 
  }



  //__________________________________________________________________
  /** 
   * @{ 
   * @name Options 
   */
  UInt_t       fOptions;      // Options 
  UInt_t       fFormats;      // Output formats
  UShort_t     fShowOthers;   // Show other data
  /* @} */
  /** 
   * @{ 
   * @name Settings 
   */
  UShort_t     fRebin;        // Rebinning factor 
  Double_t     fFwdSysErr;    // Systematic error in forward range
  Double_t     fCenSysErr;    // Systematic error in central range 
  TString      fTitle;        // Title on plot
  TString      fBase;         // Base name of output 
  TString      fClusterScale; // Scaling of clusters to tracklets      
  TString      fFinalMC;      // Final MC correction file name
  TString      fEmpirical;    // Empirical correction file name
  /* @} */
  /** 
   * @{ 
   * @name Read (or set) information 
   */
  TNamed*      fTrigString;   // Trigger string (read, or set)
  TNamed*      fNormString;   // Normalisation string (read, or set)
  TNamed*      fSNNString;    // Energy string (read, or set)
  TNamed*      fSysString;    // Collision system string (read or set)
  TAxis*       fVtxAxis;      // Vertex cuts (read or set)
  TAxis*       fCentAxis;     // Centrality axis
  Float_t      fTriggerEff;   // Trigger efficiency 
  Bool_t       fExtTriggerEff;// True if read externally 
  UShort_t     fCentMin;      // Least centrality to plot
  UShort_t     fCentMax;      // Largest centrality to plot
  /* @} */
  /** 
   * @{ 
   * @name Resulting plots 
   */
  THStack*     fResults;      // Stack of results 
  THStack*     fRatios;       // Stack of ratios 
  THStack*     fLeftRight;    // Left-right asymmetry
  TMultiGraph* fOthers;       // Older data 
  TH1*         fTriggers;     // Number of triggers
  TH1*         fTruth;        // Pointer to truth 
  /* @} */
  RangeParam*  fRangeParam;   // Parameter object for range zoom 

  static const Float_t kRightMargin;
  static const Int_t   kFont;
  static const Int_t   kAliceBlue;
  static const Int_t   kAliceRed;
  static const Int_t   kAlicePurple;
  static const Int_t   kAliceYellow;
  static const Int_t   kSysErrColor;
};

const Float_t dNdetaDrawer::kRightMargin = 0.02;
const Int_t   dNdetaDrawer::kFont        = 42; // 132 for serif
const Int_t   dNdetaDrawer::kAliceBlue   = TColor::GetColor(40,   58, 68);
const Int_t   dNdetaDrawer::kAliceRed    = TColor::GetColor(226,   0, 26);
const Int_t   dNdetaDrawer::kAlicePurple = TColor::GetColor(202,  71, 67);
const Int_t   dNdetaDrawer::kAliceYellow = TColor::GetColor(238, 125, 17);
const Int_t   dNdetaDrawer::kSysErrColor = SYSERR_COLOR;


//____________________________________________________________________
/** 
 * Function to calculate 
 * @f[
 *  g(x;A_1,A_2,\sigma_1,\sigma_2) = 
 *       A_1\left(\frac{1}{2\pi\sigma_1}e^{(x/\sigma_1)^2} - 
 *           A_2\frac{1}{2\pi\sigma_2}e^{(x/\sigma_2)^2}\right)
 * @f]
 * 
 * @param xp Pointer to x array
 * @param pp Pointer to parameter array 
 * 
 * @return @f$g(x;A_1,A_2,\sigma_1,\sigma_2)@f$
 */
Double_t myFunc(Double_t* xp, Double_t* pp)
{
  Double_t x  = xp[0];
  Double_t a1 = pp[0];
  Double_t a2 = pp[1];
  Double_t s1 = pp[2];
  Double_t s2 = pp[3];
  return a1*(TMath::Gaus(x, 0, s1) - a2 * TMath::Gaus(x, 0, s2));
}

//=== Stuff for auto zooming =========================================
/** 
 * Update canvas range 
 * 
 * @param p Parameter 
 */
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
/** 
 * Called when user changes X range 
 * 
 * @param p Parameter
 */
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

//=== Steering functions
//==============================================  
/** 
 * Display usage information
 * 
 */
void
Usage()
{
  std::ostream& o = std::cout;
  o << "Usage: DrawdNdeta(FILE,TITLE,REBIN,OTHERS,FLAGS,"
    << "SNN,SYS,TRIG,IPZMIN,IPZMAX,BASE,FMT)\n"
    << "  const char* FILE   File name to open (\"forward_dndeta.root\")\n"
    << "  const char* TITLE  Title to put on plot (\"\")\n"
    << "  UShort_t    REBIN  Rebinning factor (1)\n"
    << "  UShort_t    OTHERS Other data to draw - more below (0x7)\n"
    << "  UShort_t    FLAGS  Visualisation flags - more below (0x7)\n"
    << "  UShort_t    SYS    (optional) 1:pp, 2:PbPb, 3:pPb\n"
    << "  UShort_t    SNN    (optional) sqrt(s_NN) in GeV\n"
    << "  UShort_t    TRIG   (optional) 1: INEL, 2: INEL>0, 4: NSD, ...\n"
    << "  Float_t     EFF    (optional) Trigger efficiency\n"
    << "  Float_t     IPZMIN (optional) Least z coordinate of IP\n"
    << "  Float_t     IPZMAX (optional) Largest z coordinate of IP\n"
    << "  const char* BASE   (optional) base name of output files\n"
    << "  UShort_t    FMT    (optional) Output formats\n"
    << "\n";
  o << " OTHERS is a bit mask of\n"
    << "  0x1     Show UA5 data (INEL,NSD, ppbar, 900GeV)\n"
    << "  0x2     Show CMS data (NSD, pp)\n"
    << "  0x4     Show published ALICE data (INEL,INEL>0,NSD, pp)\n"
    << "  0x8     Show event genertor data\n"
    << "\n";
  o << " FMT is a bit mask of\n"
    << "  0x1     Make PNG output\n"
    << "  0x2     Make PDF output\n"
    << "  0x4     Make ROOT file output\n"
    << "  0x8     Make ROOT script output\n"
    << "\n";
  o << " FLAGS is a bit mask of\n"
    << "  0x1     Show ratios of data to other data and possibly MC\n"
    << "  0x2     Show left-right asymmetry\n"
    << "  0x4     Show systematic error band\n"
    << "  0x8     Show individual ring results (INEL only)\n"
    << "  0x10    Cut edges when rebinning\n"
    << "  0x20    Remove FMDxO points\n"
    << "  0x40    Apply `final MC' correction\n"
    << "  0x80    Apply `Emperical' correction\n"
    << "  0x100   Force use of MB\n"
    << "  0x200   Mirror data\n"
    << "  0x400   Export results to script\n"
    << "  0x800   Add code to do combined zooms on eta axis\n"
    << "  0x1000  Assume old-style input\n"
    << "  0x2000  Be verbose\n"
    << "  0x4000  Hi-res batch output\n"
    << "  0x8000  Add aditional white-space above results\n"
    << "\n";
  o << "0x200 requires the file forward_dndetamc.root\n"
    << "0x400 requires the file EmpiricalCorrection.root\n"
    << "To specify that you want ratios, force MB, apply empirical "
    << "correction, and export to script, set flags to\n\n"
    << "   0x1|0x100|0x80|0x400=0x581\n"
    << std::endl;

}

//____________________________________________________________________
/** 
 * Draw @f$ dN/d\eta@f$ 
 * 
 * @param filename  File name 
 * @param title     Title 
 * @param rebin     Rebinning factor 
 * @param others    What other data to show 
 * @param flags     Flags 
 * @param sNN       (optional) Collision energy [GeV]
 * @param sys       (optional) Collision system (1: pp, 2: PbPb)
 * @param trg       (optional) Trigger (1: INEL, 2: INEL>0, 4: NSD)   
 * @param eff       (optional) Trigger efficiency 
 * @param vzMin     Least @f$ v_z@f$
 * @param vzMax     Largest @f$ v_z@f$
 * @param base      Base name 
 * @param outflg    Output flags 
 *
 * @ingroup pwglf_forward_dndeta
 */
void
DrawdNdeta(const char* filename="forward_dndeta.root", 
	   const char* title="",
	   UShort_t    rebin=5, 
	   UShort_t    others=0x7,
	   UShort_t    flags=dNdetaDrawer::kDefaultOptions,
	   UShort_t    sNN=0, 
	   UShort_t    sys=0,
	   UShort_t    trg=0,
	   Float_t     eff=0,
	   UShort_t    centMin=0,
	   UShort_t    centMax=100,
	   Float_t     vzMin=999, 
	   Float_t     vzMax=-999,
	   const char* base="", 
	   UShort_t    outflg=dNdetaDrawer::kAllFormats)
{
  TString fname(filename);
  fname.ToLower();
  if (fname.CompareTo("help") == 0 || 
      fname.CompareTo("--help") == 0) { 
    Usage();
    return;
  }
  dNdetaDrawer* pd = new dNdetaDrawer;
  dNdetaDrawer& d = *pd;
  d.SetRebin(rebin);
  d.SetTitle(title);
  d.SetShowOthers(others);
  d.SetBase(base);
  // d.fClusterScale = "1.06 -0.003*x +0.0119*x*x";
  // Do the below if your input data does not contain these settings 
  if (sNN > 0) d.SetSNN(sNN);     // Collision energy per nucleon pair (GeV)
  if (sys > 0) d.SetSys(sys);     // Collision system (1:pp, 2:PbPB)
  if (trg > 0) d.SetTrigger(trg); // Collision trigger (1:INEL, 2:INEL>0, 4:NSD)
  if (eff > 0) d.SetTriggerEfficiency(eff); // Trigger efficiency
  if (vzMin < 999 && vzMax > -999) 
    d.SetVertexRange(vzMin,vzMax); // Collision vertex range (cm)
  d.SetCentralityRange(centMin,centMax); // Collision vertex range (cm)
  d.Run(filename, flags, outflg);
}
//____________________________________________________________________
//
// EOF
//

