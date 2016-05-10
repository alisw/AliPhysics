/**
 * @file   AliTrackletdNdeta2.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr 27 16:50:04 2016
 * 
 * @brief  To post processing 2nd version
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */
#ifndef ALITRACKLETDNDETA_H
#define ALITRACKLETDNDETA_H
#include <AliTrackletAODUtils.C>
#ifndef __CINT__
#include <THStack.h>
#include <TFile.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TError.h>
#include <TParameter.h>
#include <TMath.h>
#include <TLine.h>
#include <TLegendEntry.h>
#include <TStyle.h>
#include <TBrowser.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TFitResult.h>
#include <TF1.h>
#include <TSystem.h>
#else
class TPad;
class TLatex;
class TObject;
class TSeqCollection;
class TH1;
class TH2;
class TF1;
class TFitResultPtr;
class THStack;
class TCanvas;
class TVirtualPad;
class TFile;
class TAxis;
class TLegend;
class TDirectory;
#endif

//====================================================================
/**
 * Post processing 
 * 
 * @ingroup pwglf_forward_tracklets
 */
struct AliTrackletdNdeta2 : public AliTrackletAODUtils
{
  typedef AliTrackletAODUtils::Container Container;
  /**
   * Display options 
   * 
   */
  enum {
    /** Draw general information */
    kGeneral      = 0x0001,
    /** Draw parameters */
    kParameters   = 0x0002,
    /** Draw weights */
    kWeights      = 0x0004,    
    /** Draw dNch/deta */
    kdNdetas      = 0x0010,
    /** Draw delta information */   
    kDeltas       = 0x0020,
    /** Draw backgrounds */
    kBackgrounds  = 0x0040,
    /** Draw alphas */
    kAlphas       = 0x0080,
    /** Whether to make a PDF */
    kPDF          = 0x0100,
    kPNG          = 0x0100,
    /** Whether to pause after each plot */
    kPause        = 0x0200,
    /** Draw in landscape */
    kLandscape    = 0x0400,
    /** Alternative markers */
    kAltMarker    = 0x0800,
    /** Default options */
    kDefaultViz   = 0x07ff
  };
  /**
   * Calculation options 
   * 
   */
  enum {
    /** Do scaling by unity */
    kUnitScale = 0x0001,
    /** Do scaling by full average */
    kAverageScale = 0x0002,
    /** Do scaling by eta differential */
    kEtaScale = 0x0004,
    /** Do scaling by fully differential */
    kFullScale = 0x0008,
    /** MC closure test */
    kClosure      = 0x01000,
    /** Default processing options */
    kDefaultProc  = 0x00083 
  };
  enum {
    kTopBackground = kAzure-8
  };
  //==================================================================
  /** Vizualisation options */
  UInt_t   fViz;
  /** Processing options */
  UInt_t   fProc;
  /** Cut on @f$ \Delta@f$ read from real data input file */
  Double_t fDeltaCut;
  /** Lower cut on @f$\Delta@f$ tail read from real data input file */
  Double_t fTailDelta;
  /** Highest value of @f$\Delta@f$ read from real data input file */
  Double_t fMaxDelta;
  Double_t fMinK;
  Double_t fMaxK;
  Double_t fMinAlpha;
  Double_t fMaxAlpha;
  //==================================================================
  /** The canvas to draw in */
  TCanvas* fCanvas;
  /** Top bar to put title in */
  TPad*    fTop;
  /** Main body of plots */
  TPad*    fBody;
  /** Cache of last title */
  TString  fLastTitle;
  /** Cache of centrality bin title */
  TString fLastBin;
  /** Cache of formated centrality bin */
  TString fLastShort;
  /** The page header text */
  TLatex*  fHeader;

  /** 
   * Constructor
   */
  AliTrackletdNdeta2();


  /** 
   * Run it 
   * 
   * @param proc     Process mask 
   * @param viz      Visualisation mask
   * @param maxBins  Maximum number of bins to look at 
   * @param dataName Name of file from real data 
   * @param simName  Name of file from simulated data 
   * @param output   Output file name 
   */
  void Run(UInt_t      proc     = kDefaultProc,
	   UInt_t      viz      = kDefaultViz,
	   UShort_t    maxBins  = 9,
	   const char* dataName = "data.root",
	   const char* simName  = "sim.root",
	   const char* output   = 0);

  Bool_t Process(Container*  realTop,
		 Container*  simTop,
		 TDirectory* outTop,
		 Int_t       maxBins);
  /** 
   * Process a single centrality bin 
   * 
   * @param c1         Lower centrality bound
   * @param c2         Upper centrality bound 
   * @param realTop    Real container 
   * @param simTop     Simulation container 
   * @param outTop     Top-level output directory 
   * 
   * @return true on success 
   */
  Bool_t ProcessBin(Double_t    c1,
		    Double_t    c2,
		    Container*  realTop,
		    Container*  simTop,
		    TDirectory* outTop);
  /** 
   * Process a single centrality bin 
   * 
   * @param c1         Lower centrality bound
   * @param c2         Upper centrality bound 
   * @param realCont   Real container 
   * @param simCont    Simulation container 
   * @param outTop     Top-level output directory 
   * @param outDir     Current output directory 
   * @param dimen      Dimensions of k depdendence 
   * 
   * @return true on success 
   */
  Bool_t ProcessBin(Double_t    c1,
		    Double_t    c2,
		    Container*  realCont,
		    Container*  simCont,
		    TDirectory* outTop, 
		    TDirectory* outDir,
		    Int_t       dimen);
//__________________________________________________________________
  /** 
   * @{ 
   * @name @f$\Delta@f$ distributions 
   */
  /** 
   * Get @f$ \Delta@f$ distributions and scale appropriately
   * 
   * @param realCont  Container of real distributions 
   * @param simCont   Container of simulated distributions 
   * @param outParent Parent directory for output 
   * @param dim       Dimensions of dependence 
   *
   * @return true on success
   */
  Bool_t Deltas(Container*  realCont,
		Container*  simCont,
		TDirectory* outParent,
		Int_t       dim);
  /** 
   * Get @f$ \Delta@f$ distributions but only scale those from
   * injection. In this member function, we use a fixed scaling factor
   * of 1 except for the injection distribtions, where we use the
   * integrated scaling factors. That is,
   * we extract the scaling factors as
   *
   * @f[
   k_{XY} = \frac{\int_{\mathrm{tail}}d\Delta_{X}\frac{dN}{d\Delta_{X}}}{
   \int_{\mathrm{tail}}d\Delta_{Y}\frac{dN}{d\Delta_{Y}}}
   @f] 
   * where @f$ X,Y@f$ are either @f$ M, I@f$ or @f$ M', I'@f$
   * 
   * @param realCont  Container of real distributions 
   * @param simCont   Container of simulated distributions 
   * @param outParent Parent directory for output 
   *
   * @return true on success
   */
  Bool_t Deltas0D(Container*  realCont,
		  Container*  simCont,
		  TDirectory* outParent);
  /** 
   * Get @f$ \Delta@f$ distributions and scale them appropriately.  In
   * this member function, we use the integrated scaling factors
   * (e.g., as done in the past). That is, we extract the scaling
   * factors as
   *
   * @f[
   k_{XY} = \frac{\int_{\mathrm{tail}}d\Delta_{X}\frac{dN}{d\Delta_{X}}}{
   \int_{\mathrm{tail}}d\Delta_{Y}\frac{dN}{d\Delta_{Y}}}
   @f] 
   * where @f$ X,Y@f$ are one of @f$ M, M', I, I'@f$ 
   * 
   * @param realCont  Container of real distributions 
   * @param simCont   Container of simulated distributions 
   * @param outParent Parent directory for output 
   *
   * @return true on success
   */
  Bool_t Deltas1D(Container*  realCont,
		   Container*  simCont,
		   TDirectory* outParent);
  /** 
   * Get @f$ \Delta@f$ distributions and scale them appropriately.  In
   * this member function, we use the @f$\eta@f$ differential scaling
   * factors. That is, we extract the scaling factors as
   *
   * @f[
   k_{XY}(\eta) = 
   \frac{\int_{\mathrm{tail}}d\Delta_{X}\frac{d^2N}{d\Delta_{X}d\eta}}{
   \int_{\mathrm{tail}}d\Delta_{Y}\frac{d^2N}{d\Delta_{Y}d\eta}}
   @f] 
   * where @f$ X,Y@f$ are one of @f$ M, M', I, I'@f$ 
   * 
   * @param realCont  Container of real distributions 
   * @param simCont   Container of simulated distributions 
   * @param outParent Parent directory for output 
   *
   * @return true on success
   */
  Bool_t Deltas2D(Container*  realCont,
		  Container*  simCont,
		  TDirectory* outParent);
  /** 
   * Get @f$ \Delta@f$ distributions and scale them appropriately.  In
   * this member function, we use the @f$\eta,\mathrm{IP}_z@f$
   * differential scaling factors. That is, we extract the scaling
   * factors as
   *
   * @f[
   k_{XY}(\eta) = 
   \frac{
   \int_{\mathrm{tail}}d\Delta_{X}\frac{d^3N}{d\Delta_{X}d\eta d\mathrm{IP}_Z}}{
   \int_{\mathrm{tail}}d\Delta_{Y}\frac{d^3N}{d\Delta_{Y}d\eta d\mathrm{IP}_Z}}
   @f] 
   * where @f$ X,Y@f$ are one of @f$ M, M', I, I'@f$ 
   * 
   * @param realCont  Container of real distributions 
   * @param simCont   Container of simulated distributions 
   * @param outParent Parent directory for output 
   *
   * @return true on success
   */
  Bool_t Deltas3D(Container*  realCont,
		  Container*  simCont,
		  TDirectory* outParent);
  /** 
   * Write 1-dimensional @f$\Delta@f$ distributions to disk
   * 
   * @param outDir      Output directory 
   * @param realDeltaM  Real measured @f$ dN/d\Delta@f$ 
   * @param realDeltaI  Real injection @f$ dN/d\Delta@f$ 
   * @param simDeltaM   Simulated measured @f$ dN/d\Delta@f$ 
   * @param simDeltaI   Simulated injection @f$ dN/d\Delta@f$  
   * @param simDeltaC   Simulated combinatorial @f$ dN/d\Delta@f$  
   * @param simDeltaP   Simulated primary @f$ dN/d\Delta@f$  
   * @param simDeltaS   Simulated secondary @f$ dN/d\Delta@f$  
   */
  void WriteDeltas(TDirectory* outDir,
		   TH1* realDeltaM, TH1* realDeltaI,
		   TH1* simDeltaM,  TH1* simDeltaI,
		   TH1* simDeltaC,  TH1* simDeltaP,
		   TH1* simDeltaS);
  /* @} */

  /** 
   * @{ 
   * @name Result extraction 
   */
  TH1* Results(Container*  realCont,
	       Container*  simCont,
	       TDirectory* outParent,
	       Int_t       deltaDimen);
  
  /* @} */
  //____________________________________________________________________
  /**
   * @{ 
   * @name Canvas functions 
   */
  /** 
   * Clear our canvas 
   * 
   */  
  void ClearCanvas();
  /** 
   * Create our canvas 
   * 
   * @param outputName Output file name 
   */
  void CreateCanvas(const TString& outputName);
  /** 
   * Close the canvas for good 
   * 
   */
  void CloseCanvas();
  /** 
   * Print the canvas 
   * 
   * @param title       Title of this page 
   * @param shortTitle  Short title of page 
   * @param size        Size of title 
   */
  void PrintCanvas(const char* title,
		   const char* shortTitle="page",
		   Float_t     size=.7);
  /** 
   * Draw an object in a sub-pad 
   * 
   * @param c   Canvas (or pad)
   * @param pad Sub-pad number 
   * @param o   Object to draw 
   * @param opt Drawing options 
   *
   * @return if a legend is drawn, return that
   */
  TLegend* DrawInPad(TVirtualPad* c, Int_t pad, TObject* o, Option_t* opt);
  /** 
   * Modify placement of legend 
   * 
   * @param p   Pad it's drawn in 
   * @param l   The legend 
   * @param x1  New X1 coordiante in pad NDC
   * @param y1  New Y1 coordiante in pad NDC
   * @param x2  New X2 coordiante in pad NDC
   * @param y2  New Y2 coordiante in pad NDC
   */
  void ModLegend(TVirtualPad* p, TLegend* l,
		 Double_t x1, Double_t y1,
		 Double_t x2, Double_t y2);
  /** 
   * Make a stack of histograms for real and simulated data. 
   * 
   * @param name     Name of stack 
   * @param title    Title of stack
   * @param realList Container of real data 
   * @param simList  Container of simulated data 
   * @param dataOpt  Options for real data 
   * @param simOpt   Options for simulated data 
   * 
   * @return The created stack 
   */
  THStack* Make2Stack(const char*      name,
		      const char*      title,
		      Container*       realList,
		      Container*       simList,
		      Option_t*        dataOpt="",
		      Option_t*        simOpt="");
  /* @} */
  //____________________________________________________________________
  /** 
   * @{ 
   * @name Visualization functions 
   */
  Bool_t Visualize(Container*  realSums,
		   Container*  simSums,
		   Container*  realRess,
		   Container*  simRess,
		   TDirectory* outTop,
		   Int_t       maxBins);
  void   VisualizeGeneral(Container* realList, Container* simList);
  /** 
   * Draw the used simulation weights
   * 
   * @param simList Simulation list 
   */
  void VisualizeWeights(Container* simList);
  void VisualizeFinal(TDirectory* outDir, Int_t i);
  Bool_t VisualizeBin(Double_t    c1,
		      Double_t    c2,
		      Container*  simList, 
		      TDirectory* outTop);
  Bool_t VisualizeSpecies(Container* simCont);
  Bool_t VisualizePrimary(Container* simCont);
  Bool_t VisualizeDelta(TDirectory* outTop, Int_t dimen);
  Bool_t VisualizeResult(TDirectory* outTop,  Int_t       dimen);
  /* @} */
  //____________________________________________________________________
  /** 
   * @{ 
   * @name Drawing parameters 
   */
  /** 
   * Draw parameter name and value 
   * 
   * @param name  Name of parameter 
   * @param y     Current y coordinate 
   * @param val   Value 
   */
  void VisualizeParam(const char* name, Double_t& y,  const char* val);
  /** 
   * Draw a real valued paraemeter 
   * 
   * @param name  Name of parameter 
   * @param y     Current y coordinate 
   * @param val   Value 
   */  
  void VisualizeParam(const char* name, Double_t& y,  Double_t val);
  /** 
   * Draw an integer valued paraemeter 
   * 
   * @param name  Name of parameter 
   * @param y     Current y coordinate 
   * @param val   Value 
   */  
  void VisualizeParam(const char* name, Double_t& y,  Int_t val);
  /** 
   * Draw a boolean valued paraemeter 
   * 
   * @param name  Name of parameter 
   * @param y     Current y coordinate 
   * @param val   Value 
   */  
  void VisualizeParam(const char* name, Double_t& y,  Bool_t val);
  /** 
   * Draw all parameters in a container 
   * 
   * @param pars 
   * @param title 
   * @param comb if true, write for combinatorial background 
   */
  void VisualizeParams(Container* pars, const char* title);
  /** 
   * Draw parameters from both real and simulated analysis 
   * 
   * @param realSums Real data 
   * @param simSums  simulated data 
   */
  void VisualizeParams(Container* realSums, Container* simSums);
  /* @} */
  
  //____________________________________________________________________
  /** 
   * @{
   * @name Some utilities 
   */
  /** 
   * Open a file 
   * 
   * @param filename File name 
   * 
   * @return Opened file handle or null 
   */
  TFile* OpenFile(const char* filename);
  /** 
   * Set attributes on a histogram 
   * 
   * @param h      Histogram 
   * @param color  Color of histogram 
   * @param marker Marker style 
   * @param size   Marker size 
   * @param fill   Fill style 
   * @param line   Line style 
   * @param width  Line width 
   * 
   * @return The modified histogram
   */
  TH1* SetAttr(TH1* h,
	       Color_t  color,
	       Style_t  marker=20,
	       Double_t size=1.,
	       Style_t  fill=0,
	       Style_t  line=1,
	       Width_t  width=1);
  /** 
   * The observable title 
   * 
   * @return Observable title 
   */
  const char* ObsTitle() const { return "d#it{N}_{ch}/d#eta"; }
  /** 
   * Get name of centrality folder 
   * 
   * @param c1 Lower bound 
   * @param c2 Upper bound 
   * 
   * @return String 
   */
  const char* CentName(Double_t c1, Double_t c2);
  /** 
   * Get mean of a histogram content by fitting with a constant. 
   * 
   * @param h The histogram to get the mean off
   * @param e On return, the fit error 
   * 
   * @return The fit parameter
   */
  Double_t MeanY(TH1* h, Double_t& e);
  /** 
   * Get mean of a histogram content by fitting with a constant. 
   * 
   * @param h The histogram to get the mean off
   * @param e On return, the fit error 
   * 
   * @return The fit parameter
   */
  Double_t MeanZ(TH2* h, Double_t& e);
  /* @} */
};

//====================================================================  
struct SuppressGuard
{
  Int_t save = 0;
  SuppressGuard(Int_t lvl=2000)
  {
    save = gErrorIgnoreLevel;
    gErrorIgnoreLevel = lvl;
  }
  ~SuppressGuard()
  {
    gErrorIgnoreLevel = save;
  }
};

//====================================================================
AliTrackletdNdeta2::AliTrackletdNdeta2()
  : AliTrackletAODUtils(),
    fProc(0),
    fViz(0),
    fDeltaCut(0),
    fTailDelta(0),
    fMaxDelta(0),
    fMinK(.7),
    fMaxK(1.5),
    fMinAlpha(0),
    fMaxAlpha(2.5),
    fCanvas(0),
    fTop(0),
    fBody(0)    
{}

//====================================================================
void AliTrackletdNdeta2::Run(UInt_t      proc,
			     UInt_t      viz, 
			     UShort_t    maxBins,
			     const char* dataName,
			     const char* simName,
			     const char* outName)
{
  // Store options 
  fProc = proc;
  fViz  = viz;

  // Open the input files
  TFile* dataFile = 0;
  TFile* simFile  = 0;
  if (!(dataFile = OpenFile(dataName))) return;
  if (!(simFile  = OpenFile(simName)))  return;

  // Get some top-level contianers 
  const char* base     = "MidRapidity";
  Container*  realSums = GetC(dataFile, Form("%sSums",      base));
  Container*  realRess = GetC(dataFile, Form("%sResults",   base));
  Container*  simSums  = GetC(simFile,  Form("%sMCSums",    base));
  Container*  simRess  = GetC(simFile,  Form("%sMCResults", base));
  if (!realSums || !realRess || !simSums || !simRess) return;

  // Get parameters from the real data file 
  Container* params   = GetC(realSums, "parameters");
  fDeltaCut  = GetD(params, "DeltaCut");
  fTailDelta = GetD(params, "TailDelta");
  fMaxDelta  = GetD(params, "MaxDelta");

  // Create output file name 
  TString outBase(outName);
  if (outBase.IsNull())          outBase.Form("MiddNdeta_0x%04x", fProc);
  if (outBase.EndsWith(".root")) outBase.ReplaceAll(".root", "");  
  // Open the output file 
  TFile* out = TFile::Open(Form("%s.root", outBase.Data()), "RECREATE");
  
  
  Process(realRess, simRess, out, maxBins);
  out->Write();

  Visualize(realSums, simSums, realRess, simRess, out, maxBins);
}

//====================================================================
Bool_t AliTrackletdNdeta2::Process(Container*  realTop,
				   Container*  simTop,
				   TDirectory* outDir,
				   Int_t       maxBins)
{
  TH1* realCent = CopyH1(realTop, "cent", "realCent");
  TH1* simCent  = CopyH1(simTop,  "cent", "simCent");
  TH1* realIPz  = GetH1(realTop,  "ipz");
  TH1* simIPz   = GetH1(simTop,   "ipz");

  // Check consistency of found histograms 
  if (!CheckConsistency(realCent, simCent)) {
    Warning("Post", "Centrality bins are incompatible, giving up");
    return false;
  }
  if (!CheckConsistency(realIPz, simIPz)) {
    Warning("Post", "Centrality bins are incompatible, giving up");
    return false;
  }
  
  // Check if we're doing a closure test 
  if (fProc & kClosure) realTop = simTop;


  THStack* mids = new THStack("mids", ""); 
  {
    Int_t    nbin   = 9;
    Double_t bins[] = { 0., 5., 10., 20., 30., 40., 50., 60., 70., 80.,  };
    Double_t vals[] = { 1948, 1590, 1180, 786, 512, 318, 183, 96.3, 44.9, };
    Double_t errs[] = {   38,   32,   31,  20,  15,  12,   8,  5.8,  3.4, };
    TH1*     publ   = new TH1D("published", "dNch/deta in |eta|<0.5",
			       nbin, bins);
    publ->SetMarkerStyle(24);
    publ->SetMarkerColor(kBlack);
    publ->SetMarkerSize(1.3);
    for (Int_t i = 0; i < nbin; i++) {
      publ->SetBinContent(i+1, vals[i]);
      publ->SetBinError  (i+1, errs[i]);
    }
    mids->Add(publ);
  }
  
  // Make histogram for mid-rapidiy results
  for (Int_t d = 0; d < 4; d++) {
    if ((fProc & (1 << d)) == 0) continue;
    TDirectory* dd  = outDir->mkdir(Form("final%dd", d));
    TH1*        mid = Make1D(0,"mid",Form("%s|_{|#eta|<0.5}", ObsTitle()),
			     2+d, 20, *(realCent->GetXaxis()));    
    mid->SetDirectory(dd);
    mid->SetXTitle("Centrality [%]");
    mid->SetYTitle(mid->GetTitle());
    mids->Add(mid);
    THStack* full = new THStack("full","");
    dd->cd();
    full->Write();
    outDir->cd();
  }

  // Write centralities to file 
  realCent->Write();
  simCent ->Write();

  // Loop over defined centrality bins 
  for (Int_t i = 1; i <= realCent->GetNbinsX() && i <= maxBins ; i++) {
    Double_t c1 = realCent->GetXaxis()->GetBinLowEdge(i);
    Double_t c2 = realCent->GetXaxis()->GetBinUpEdge (i);
      
    ProcessBin(c1, c2, realTop, simTop, outDir);
  }
  // Close output file
  outDir->cd();
  mids->Write();

  return true;
}


//____________________________________________________________________
Bool_t AliTrackletdNdeta2::ProcessBin(Double_t    c1,
				      Double_t    c2,
				      Container*  realTop,
				      Container*  simTop,
				      TDirectory* outTop)
{
  // Form the folder name
  TString centName(CentName(c1,c2));

  // Get centrality containers 
  Container* realCont = GetC(realTop, centName);
  Container* simCont  = GetC(simTop,  centName);
  if (!realCont || !simCont) return false;

  TDirectory* outDir = outTop->mkdir(centName);

  Printf("%5.1f - %5.1f%%", c1, c2);
  for (Int_t i = 0; i < 4; i++) {
    if ((fProc & (1 << i)) == 0) continue;
    if (!ProcessBin(c1, c2, realCont, simCont, outTop, outDir, i))
      return false;
  }
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletdNdeta2::ProcessBin(Double_t    c1,
				      Double_t    c2,
				      Container*  realCont,
				      Container*  simCont,
				      TDirectory* outTop, 
				      TDirectory* outDir,
				      Int_t       dimen)
{
  if (!Deltas(realCont, simCont, outDir, dimen)) return false;
  TH1* dndeta = Results(realCont, simCont, outDir, dimen);
  if (!dndeta) {
    Warning("ProcessBin", "Failed on Deltas for %f - %f", c1, c2);
    return false;
  }

  TDirectory* final = outTop->GetDirectory(Form("final%dd", dimen));
  if (!final) {
    Warning("ProcessBin", "Failed on results for %f - %f", c1, c2);
    return false;
  }

  TH1*     mid  = static_cast<TH1*>    (GetO(final, "mid"));
  THStack* full = static_cast<THStack*>(GetO(final, "full"));
  if (!mid || !final) {
    Warning("ProcessBin", "Missing one of mid (%p) or full (%p)", mid, full);
    return false;
  }

  TF1* f = static_cast<TF1*>(dndeta->GetListOfFunctions()->At(0));
  if (!f) {
    Warning("ProcessBin", "No fit found on %s", dndeta->GetTitle());
    return false;
  }

  Double_t c = (c1+c2)/2;
  Int_t    b = mid->GetXaxis()->FindBin(c);
  if (b < 1 || b > mid->GetNbinsX()) {
    Warning("ProcessBin", "Centrality %f - %f out of range", c1, c2);
    return false;
  }

  mid->SetBinContent(b, f->GetParameter(0));
  mid->SetBinError  (b, f->GetParError (0));
  
  const Color_t cc[] = { kMagenta+2, // 0
			 kBlue+2,    // 1
			 kAzure-1,   // 2 // 10,
			 kCyan+2,    // 3
			 kGreen+1,   // 4 
			 kSpring+5,  // 5 //+10,
			 kYellow+1,  // 6
			 kOrange+5,  // 7 //+10,
			 kRed+1,     // 8
			 kPink+5,    // 9 //+10,
			 kBlack };   // 10
  Color_t tc = cc[b % 10];
  TH1*    copy = static_cast<TH1*>(dndeta->Clone(outDir->GetName()));
  copy->SetDirectory(final);
  copy->GetListOfFunctions()->Clear();
  copy->SetTitle(Form("%5.1f#minus%5.1f%%", c1, c2));
  SetAttr(copy, tc);
  full->Add(copy);
  final->cd();
  full->Write(full->GetName(), TObject::kOverwrite);
  
  return true;
}
  

//====================================================================
Bool_t AliTrackletdNdeta2::Deltas(Container*  realCont,
				  Container*  simCont,
				  TDirectory* outParent,
				  Int_t       dim)
{
  switch (dim) {
  case 0: return Deltas0D(realCont, simCont, outParent);
  case 1: return Deltas1D(realCont, simCont, outParent);
  case 2: return Deltas2D(realCont, simCont, outParent);
  case 3: return Deltas3D(realCont, simCont, outParent);
  }
  return false;
}
//____________________________________________________________________
Bool_t AliTrackletdNdeta2::Deltas0D(Container*  realCont,
				  Container*  simCont,
				    TDirectory* outParent)
{
  // Make an output directory 
  TDirectory* outDir = outParent->mkdir("delta0d");

  // Get the real and simulated measured folders 
  Container* realMeas  = GetC(realCont, "measured");
  Container* simMeas   = GetC(simCont,  "measured");
  if (!realMeas || !simMeas) return false;

  // Create a flat 2D scaling histogram for later user 
  TH2*       h         = CopyH2(realMeas, "etaIPz", "scale");
  TH1*       hp        = h->ProjectionX("scaleProj");
  hp->Reset();
  hp->SetMinimum(fMinK);
  hp->SetMaximum(fMaxK);
  hp->SetYTitle("k");
  h->SetZTitle("k");
  h->SetTitle("k=1");
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    hp->SetBinContent(i, 1);
    hp->SetBinError  (i, 0);
    for (Int_t j = 1; j <= h->GetNbinsY(); j++) {
      if (h->GetBinContent(i,j) < 1e-6) {
	h->SetBinContent(i,j,0);
	h->SetBinError  (i,j,0);
      }
      else {
	h->SetBinContent(i,j,1);
	h->SetBinError  (i,j,0);
      }
    }
  }
  h->SetDirectory(outDir);
  h->SetMinimum(fMinK);
  h->SetMaximum(fMaxK);
  hp->SetDirectory(outDir);

  // Get the raw projected Delta distributions for each component 
  TH1* realDeltaM = CopyH1(realMeas,                      "delta","realDeltaM");
  TH1* realDeltaI = CopyH1(GetC(realCont,"injected"),     "delta","realDeltaI");
  TH1* simDeltaM  = CopyH1(simMeas,                       "delta","simDeltaM");
  TH1* simDeltaI  = CopyH1(GetC(simCont, "injected"),     "delta","simDeltaI");
  TH1* simDeltaC  = CopyH1(GetC(simCont, "combinatorics"),"delta","simDeltaC");
  TH1* simDeltaP  = CopyH1(GetC(simCont, "primaries"),    "delta","simDeltaP");
  TH1* simDeltaS  = CopyH1(GetC(simCont, "secondaries"),  "delta","simDeltaS");

  // Get integrated scaling factor for injections, and scale the
  // injection distributions by that 
  Double_t realScaleI  = GetD(GetC(realCont,"injected"), "scale");
  Double_t realScaleIE = GetD(GetC(realCont,"injected"), "scaleError");
  Double_t simScaleI   = GetD(GetC(simCont, "injected"), "scale");
  Double_t simScaleIE  = GetD(GetC(simCont, "injected"), "scaleError");
  Scale(realDeltaI, realScaleI, realScaleIE);
  Scale(simDeltaI,  simScaleI,  simScaleIE);
  realDeltaI->SetTitle(Form("k_{I}#times%s",realScaleI,realDeltaI->GetTitle()));
  simDeltaI ->SetTitle(Form("k_{I'}#times%s",simScaleI,simDeltaI ->GetTitle()));
  
  WriteDeltas(outDir, realDeltaM, realDeltaI, simDeltaM, simDeltaI,
	      simDeltaC, simDeltaP, simDeltaS);

  outParent->cd();
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletdNdeta2::Deltas1D(Container*  realCont,
				    Container*  simCont,
				    TDirectory* outParent)
{
  // Make an output directory 
  TDirectory* outDir = outParent->mkdir("delta1d");

  // Get the real and simulated measured folders 
  Container* realMeas  = GetC(realCont, "measured");
  Container* simMeas   = GetC(simCont,  "measured");
  if (!realMeas || !simMeas) return false;

  // Get the integrated tails of the real and simulated observed
  // distribtutions, and calcualte scaling factor.
  Double_t   realTail  = GetD(realMeas, "deltaTail");
  Double_t   realTailE = GetD(realMeas, "deltaTailError");
  Double_t   simTail   = GetD(simMeas,  "deltaTail");
  Double_t   simTailE  = GetD(simMeas,  "deltaTailError");
  Double_t   scaleE    = 0;
  Double_t   scale     = RatioE(realTail, realTailE, simTail, simTailE, scaleE);

  // Create a flat 2D scaling histogram for later user 
  TH2*       h         = CopyH2(realMeas, "etaIPz", "scale");
  TH1*       hp        = h->ProjectionX("scaleProj");
  hp->Reset();
  hp->SetMinimum(fMinK);
  hp->SetMaximum(fMaxK);
  hp->SetYTitle("k");
  h->SetZTitle("k");
  h->SetTitle(Form("k=%5.3f#pm%5.3f", scale, scaleE));
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    hp->SetBinContent(i, scale);
    hp->SetBinError  (i, scaleE);
    for (Int_t j = 1; j <= h->GetNbinsY(); j++) {
      if (h->GetBinContent(i,j) < 1e-6) {
	h->SetBinContent(i,j,0);
	h->SetBinError  (i,j,0);
      }
      else {
	h->SetBinContent(i,j,scale);
	h->SetBinError  (i,j,scaleE);
      }
    }
  }
  h->SetDirectory(outDir);
  h->SetMinimum(fMinK);
  h->SetMaximum(fMaxK);
  hp->SetDirectory(outDir);

  // Get the raw projected Delta distributions for each component 
  TH1* realDeltaM = CopyH1(realMeas,                      "delta","realDeltaM");
  TH1* realDeltaI = CopyH1(GetC(realCont,"injected"),     "delta","realDeltaI");
  TH1* simDeltaM  = CopyH1(simMeas,                       "delta","simDeltaM");
  TH1* simDeltaI  = CopyH1(GetC(simCont, "injected"),     "delta","simDeltaI");
  TH1* simDeltaC  = CopyH1(GetC(simCont, "combinatorics"),"delta","simDeltaC");
  TH1* simDeltaP  = CopyH1(GetC(simCont, "primaries"),    "delta","simDeltaP");
  TH1* simDeltaS  = CopyH1(GetC(simCont, "secondaries"),  "delta","simDeltaS");

  // Get integrated scaling factor for injections, and scale the
  // injection distributions by that 
  Double_t realScaleI  = GetD(GetC(realCont,"injected"), "scale");
  Double_t realScaleIE = GetD(GetC(realCont,"injected"), "scaleError");
  Double_t simScaleI   = GetD(GetC(simCont, "injected"), "scale");
  Double_t simScaleIE  = GetD(GetC(simCont, "injected"), "scaleError");
  Scale(realDeltaI, realScaleI, realScaleIE);
  Scale(simDeltaI,  simScaleI,  simScaleIE);
  realDeltaI->SetTitle(Form("k_{I}#times%s", realDeltaI->GetTitle()));
  simDeltaI ->SetTitle(Form("k_{I'}#times%s",simDeltaI ->GetTitle()));
  
  TH1* toScale[]  = { simDeltaM,simDeltaI,simDeltaC,simDeltaP,simDeltaS,0};
  TH1**    pScale     = toScale;
  while ((*pScale)) { 
    Scale(*pScale, scale, scaleE);
    (*pScale)->SetTitle(Form("k_{M}#times%s", (*pScale)->GetTitle()));
    pScale++;
  }

  WriteDeltas(outDir, realDeltaM, realDeltaI, simDeltaM, simDeltaI,
	      simDeltaC, simDeltaP, simDeltaS);

  outParent->cd();
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletdNdeta2::Deltas2D(Container*  realCont,
				    Container*  simCont,
				    TDirectory* outParent)
{
  // Make an output directory 
  TDirectory* outDir = outParent->mkdir("delta2d");

  // Get the real and simulated measured folders 
  Container* realMeas  = GetC(realCont, "measured");
  Container* simMeas   = GetC(simCont,  "measured");
  if (!realMeas || !simMeas) return false;

  // Get the eta-differential tails of the real and simulated observed
  // distribtutions, and calcualte scaling factor.
  TH1*     realTail  = GetH1(realMeas, "etaDeltaTail");
  TH1*     simTail   = GetH1(simMeas,  "etaDeltaTail");
  TH1*     scale     = static_cast<TH1*>(realTail->Clone("scaleProj"));
  scale->Divide(simTail);
  scale->SetYTitle("k");
  Double_t sE, s     = MeanY(scale, sE);
  TGraphErrors* g = new TGraphErrors(2);
  g->SetLineStyle(2);
  g->SetLineColor(kBlack);
  g->SetFillColor(kYellow);
  g->SetFillStyle(3002);
  g->SetPoint(0, scale->GetXaxis()->GetXmin(), s); g->SetPointError(0, 0, sE);
  g->SetPoint(1, scale->GetXaxis()->GetXmax(), s); g->SetPointError(1, 0, sE);
  scale->GetListOfFunctions()->Add(g, "le3");
  
  // Create a flat 2D scaling histogram for later user 
  TH2*       h         = CopyH2(realMeas, "etaIPz", "scale");  
  h->SetZTitle("k");
  h->SetTitle(Form("#LTk(#eta)#GT_{#eta}=%5.3f#pm%5.3f", s, sE));
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) { // eta 
    for (Int_t j = 1; j <= h->GetNbinsY(); j++) { // ipz
      if (h->GetBinContent(i,j) < 1e-6) {
	h->SetBinContent(i,j,0);
	h->SetBinError  (i,j,0);
      }
      else {
	h->SetBinContent(i,j,scale->GetBinContent(i));
	h->SetBinError  (i,j,scale->GetBinError  (i));
      }
    }
  }
  h->SetMinimum(fMinK);
  h->SetMaximum(fMaxK);
  h->SetDirectory(outDir);
  scale->SetMinimum(fMinK);
  scale->SetMaximum(fMaxK);
  scale->SetDirectory(outDir);

  // Get the raw projected Delta distributions for each component 
  TH2* r2DeltaM = CopyH2(realMeas,                      "etaDelta","r2DeltaM");
  TH2* r2DeltaI = CopyH2(GetC(realCont,"injected"),     "etaDelta","r2DeltaI");
  TH2* s2DeltaM = CopyH2(simMeas,                       "etaDelta","s2DeltaM");
  TH2* s2DeltaI = CopyH2(GetC(simCont, "injected"),     "etaDelta","s2DeltaI");
  TH2* s2DeltaC = CopyH2(GetC(simCont, "combinatorics"),"etaDelta","s2DeltaC");
  TH2* s2DeltaP = CopyH2(GetC(simCont, "primaries"),    "etaDelta","s2DeltaP");
  TH2* s2DeltaS = CopyH2(GetC(simCont, "secondaries"),  "etaDelta","s2DeltaS");

  // Get eta-differential scaling factor for injections, and scale the
  // injection distributions by that
  TH1* rScaleI  = GetH1(GetC(realCont,"injected"), "etaScale");
  TH1* sScaleI  = GetH1(GetC(simCont, "injected"), "etaScale");
  Double_t rIE, rI = MeanY(rScaleI, rIE);  
  Double_t sIE, sI = MeanY(sScaleI, sIE);  
  Scale(r2DeltaI, rScaleI);
  Scale(s2DeltaI, sScaleI);
  r2DeltaI ->SetTitle(Form("#LTk_{I}#GT_{#eta}#times%s",
			   r2DeltaI->GetTitle()));
  s2DeltaI ->SetTitle(Form("#LTk_{I'}#GT_{#eta}#times%s",
			   s2DeltaI->GetTitle()));
  
  TH2*  toScale[] = { s2DeltaM,s2DeltaI,s2DeltaC,s2DeltaP,s2DeltaS,0};
  TH2** pScale    = toScale;
  while ((*pScale)) { 
    Scale(*pScale, scale);
    (*pScale)->SetTitle(Form("#LTk_{M}#GT_{#eta}#times%s",
			     (*pScale)->GetTitle()));
    pScale++;
  }

  TH1* rDeltaM = ProjectDelta(r2DeltaM); 
  TH1* rDeltaI = ProjectDelta(r2DeltaI); 
  TH1* sDeltaM = ProjectDelta(s2DeltaM); 
  TH1* sDeltaI = ProjectDelta(s2DeltaI); 
  TH1* sDeltaC = ProjectDelta(s2DeltaC); 
  TH1* sDeltaP = ProjectDelta(s2DeltaP); 
  TH1* sDeltaS = ProjectDelta(s2DeltaS); 
  rDeltaM->SetTitle(r2DeltaM->GetTitle()); rDeltaM->SetName("realDeltaM");
  rDeltaI->SetTitle(r2DeltaI->GetTitle()); rDeltaI->SetName("realDeltaI");
  sDeltaM->SetTitle(s2DeltaM->GetTitle()); sDeltaM->SetName("simDeltaM");
  sDeltaI->SetTitle(s2DeltaI->GetTitle()); sDeltaI->SetName("simDeltaI");
  sDeltaC->SetTitle(s2DeltaC->GetTitle()); sDeltaC->SetName("simDeltaC");
  sDeltaP->SetTitle(s2DeltaP->GetTitle()); sDeltaP->SetName("simDeltaP");  
  sDeltaS->SetTitle(s2DeltaS->GetTitle()); sDeltaS->SetName("simDeltaS");
  
  WriteDeltas(outDir,rDeltaM,rDeltaI,sDeltaM,sDeltaI,sDeltaC,sDeltaP,sDeltaS);

  TDirectory* full = outDir->mkdir("full");
  r2DeltaM->SetDirectory(full); r2DeltaM->SetName("realDeltaM");
  r2DeltaI->SetDirectory(full); r2DeltaI->SetName("realDeltaI");
  s2DeltaM->SetDirectory(full); s2DeltaM->SetName("simDeltaM");
  s2DeltaI->SetDirectory(full); s2DeltaI->SetName("simDeltaI");
  s2DeltaC->SetDirectory(full); s2DeltaC->SetName("simDeltaC");
  s2DeltaP->SetDirectory(full); s2DeltaP->SetName("simDeltaP");
  s2DeltaS->SetDirectory(full); s2DeltaS->SetName("simDeltaS");

  outParent->cd();
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletdNdeta2::Deltas3D(Container*  realCont,
				    Container*  simCont,
				    TDirectory* outParent)
{
  // Make an output directory 
  TDirectory* outDir = outParent->mkdir("delta3d");

  // Get the real and simulated measured folders 
  Container* realMeas  = GetC(realCont, "measured");
  Container* simMeas   = GetC(simCont,  "measured");
  if (!realMeas || !simMeas) return false;

  // Get the eta-differential tails of the real and simulated observed
  // distribtutions, and calcualte scaling factor.
  TH2*     realTail  = GetH2(realMeas, "etaIPzDeltaTail");
  TH2*     simTail   = GetH2(simMeas,  "etaIPzDeltaTail");
  TH2*     scale     = static_cast<TH2*>(realTail->Clone("scale"));
  scale->SetDirectory(0);
  scale->Divide(simTail);
  Double_t sE, s     = MeanZ(scale, sE);  
  scale->SetZTitle("k");
  scale->SetTitle(Form("#LTk(#eta)#GT_{#eta,IP_{#it{z}}}=%5.3f#pm%5.3f",
		       s, sE));
  scale->SetDirectory(outDir);
  scale->SetMinimum(fMinK);
  scale->SetMaximum(fMaxK);
  TH1* etaScale = AverageOverIPz(scale, "scaleProj", 1, 0, 0);
  etaScale->SetYTitle("#LTk#GT_{IP_{z}}");
  etaScale->SetDirectory(outDir);
  etaScale->SetMinimum(fMinK);
  etaScale->SetMaximum(fMaxK);
  TGraphErrors* g = new TGraphErrors(2);
  g->SetLineStyle(2);
  g->SetLineColor(kBlack);
  g->SetFillColor(kYellow);
  g->SetFillStyle(3002);
  g->SetPoint(0, etaScale->GetXaxis()->GetXmin(), s); 
  g->SetPoint(1, etaScale->GetXaxis()->GetXmax(), s);
  g->SetPointError(0, 0, sE);
  g->SetPointError(1, 0, sE);
  etaScale->GetListOfFunctions()->Add(g, "le3");
  
  // Get the raw projected Delta distributions for each component
  const char* nm = "etaDeltaIPz";
  TH3* r3DeltaM = CopyH3(realMeas,                      nm,"r3DeltaM");
  TH3* r3DeltaI = CopyH3(GetC(realCont,"injected"),     nm,"r3DeltaI");
  TH3* s3DeltaM = CopyH3(simMeas,                       nm,"s3DeltaM");
  TH3* s3DeltaI = CopyH3(GetC(simCont, "injected"),     nm,"s3DeltaI");
  TH3* s3DeltaC = CopyH3(GetC(simCont, "combinatorics"),nm,"s3DeltaC");
  TH3* s3DeltaP = CopyH3(GetC(simCont, "primaries"),    nm,"s3DeltaP");
  TH3* s3DeltaS = CopyH3(GetC(simCont, "secondaries"),  nm,"s3DeltaS");

  // Get eta-differential scaling factor for injections, and scale the
  // injection distributions by that
  TH2* rScaleI  = GetH2(GetC(realCont,"injected"), "etaIPzScale");
  TH2* sScaleI  = GetH2(GetC(simCont, "injected"), "etaIPzScale");
  Double_t rIE, rI = MeanZ(rScaleI, rIE);  
  Double_t sIE, sI = MeanZ(sScaleI, sIE);  
  ScaleDelta(r3DeltaI, rScaleI);
  ScaleDelta(s3DeltaI, sScaleI);
  r3DeltaI ->SetTitle(Form("#LTk_{I}#GT_{#eta,IP_{#it{z}}}#times%s",
			   r3DeltaI->GetTitle()));
  s3DeltaI ->SetTitle(Form("#LTk_{I'}#GT_{#eta,IP_{#it{z}}}#times%s",
			   s3DeltaI->GetTitle()));
  
  TH3*  toScale[] = { s3DeltaM,s3DeltaI,s3DeltaC,s3DeltaP,s3DeltaS,0};
  TH3** pScale    = toScale;
  while ((*pScale)) { 
    ScaleDelta(*pScale, scale);
    (*pScale)->SetTitle(Form("#LTk_{M}#GT_{#eta,IP_{#it{z}}}#times%s",
			     (*pScale)->GetTitle()));
    pScale++;
  }

  TH1* rDeltaM = ProjectDeltaFull(r3DeltaM); 
  TH1* rDeltaI = ProjectDeltaFull(r3DeltaI); 
  TH1* sDeltaM = ProjectDeltaFull(s3DeltaM); 
  TH1* sDeltaI = ProjectDeltaFull(s3DeltaI); 
  TH1* sDeltaC = ProjectDeltaFull(s3DeltaC); 
  TH1* sDeltaP = ProjectDeltaFull(s3DeltaP); 
  TH1* sDeltaS = ProjectDeltaFull(s3DeltaS); 
  rDeltaM->SetTitle(r3DeltaM->GetTitle()); rDeltaM->SetName("realDeltaM");
  rDeltaI->SetTitle(r3DeltaI->GetTitle()); rDeltaI->SetName("realDeltaI");
  sDeltaM->SetTitle(s3DeltaM->GetTitle()); sDeltaM->SetName("simDeltaM");
  sDeltaI->SetTitle(s3DeltaI->GetTitle()); sDeltaI->SetName("simDeltaI");
  sDeltaC->SetTitle(s3DeltaC->GetTitle()); sDeltaC->SetName("simDeltaC");
  sDeltaP->SetTitle(s3DeltaP->GetTitle()); sDeltaP->SetName("simDeltaP");  
  sDeltaS->SetTitle(s3DeltaS->GetTitle()); sDeltaS->SetName("simDeltaS");
    
  WriteDeltas(outDir,rDeltaM,rDeltaI,sDeltaM,sDeltaI,sDeltaC,sDeltaP,sDeltaS);

  TDirectory* full = outDir->mkdir("full");
  r3DeltaM->SetDirectory(full); r3DeltaM->SetName("realDeltaM");
  r3DeltaI->SetDirectory(full); r3DeltaI->SetName("realDeltaI");
  s3DeltaM->SetDirectory(full); s3DeltaM->SetName("simDeltaM");
  s3DeltaI->SetDirectory(full); s3DeltaI->SetName("simDeltaI");
  s3DeltaC->SetDirectory(full); s3DeltaC->SetName("simDeltaC");
  s3DeltaP->SetDirectory(full); s3DeltaP->SetName("simDeltaP");
  s3DeltaS->SetDirectory(full); s3DeltaS->SetName("simDeltaS");

  outParent->cd();
  return true;
}

//____________________________________________________________________
void AliTrackletdNdeta2::WriteDeltas(TDirectory* outDir,
				     TH1* realDeltaM, TH1* realDeltaI,
				     TH1* simDeltaM,  TH1* simDeltaI,
				     TH1* simDeltaC,  TH1* simDeltaP,
				     TH1* simDeltaS)
{
  THStack* all    = new THStack("all", "");
  SetAttr(realDeltaM, kRed+2,    20, 1.0);
  SetAttr(realDeltaI, kOrange+2, 21, 1.0);
  realDeltaM->SetDirectory(outDir);
  realDeltaI->SetDirectory(outDir);
  all->Add(realDeltaM);
  all->Add(realDeltaI);
  
  TH1*     toScale[]  = { simDeltaM,simDeltaI,simDeltaC,simDeltaP,simDeltaS,0};
  Color_t  toColor[]  = { kRed,     kOrange,  kMagenta, kGreen,   kBlue };
  Style_t  toStyle[]  = { 24,       25,       30,       26,       32    };
  TH1**    pScale     = toScale;
  Color_t* pColor     = toColor;
  Style_t* pStyle     = toStyle;
  while ((*pScale)) { 
    (*pScale)->SetDirectory(outDir);
    all->Add((*pScale));
    SetAttr(*pScale, (*pColor)+2, *pStyle, 1.2);
    pScale++;
    pColor++;
    pStyle++;
  }
  outDir->cd();
  all->Write();

  THStack* ratios = new THStack("ratios", "");
  TH1*     ratioM = static_cast<TH1*>(simDeltaM->Clone("ratioM"));
  ratioM->SetTitle("#Delta_{M'}/#Delta_{M}");
  ratioM->Divide(realDeltaM);
  ratioM->SetDirectory(outDir);
  ratios->Add(ratioM);

  TH1*     ratioI = static_cast<TH1*>(simDeltaI->Clone("ratioI"));
  ratioI->SetTitle("#Delta_{I'}/#Delta_{I}");
  ratioI->Divide(realDeltaI);
  ratioI->SetDirectory(outDir);
  ratios->Add(ratioI);

  TH1*     ratioIC = static_cast<TH1*>(simDeltaC->Clone("ratioIC"));
  ratioIC->SetTitle("#Delta_{C'}/#Delta_{I}");
  ratioIC->Divide(realDeltaI);
  ratioIC->SetDirectory(outDir);
  ratios->Add(ratioIC);

  ratios->Write();
}


//====================================================================
TH1* AliTrackletdNdeta2::Results(Container*  realCont,
				 Container*  simCont,
				 TDirectory* outParent,
				 Int_t       deltaDimen)
{
  TDirectory* outDir = outParent->mkdir(Form("results%dd", deltaDimen));
  TDirectory* delDir = outParent->GetDirectory(Form("delta%dd", deltaDimen));

  TH2* scale = static_cast<TH2*>(delDir->Get("scale"));
  TH2* realM = CopyH2(GetC(realCont, "measured"),      "etaIPz", "realM");
  TH2* realS = CopyH2(GetC(realCont, "measured"),      "etaIPz", "realS");
  TH2* realC = CopyH2(GetC(realCont, "measured"),      "etaIPz", "realC");
  TH2* simM  = CopyH2(GetC(simCont,  "measured"),      "etaIPz", "simM");
  TH2* simC  = CopyH2(GetC(simCont,  "combinatorics"), "etaIPz", "simC");
  TH2* simG  = CopyH2(GetC(simCont,  "generated"),     "etaIPz", "simG");
  TH2* simS  = CopyH2(GetC(simCont,  "measured"),      "etaIPz", "simS");
  TH2* simA  = CopyH2(GetC(simCont,  "generated"),     "etaIPz", "simA");
  TH2* simB  = CopyH2(GetC(simCont,  "combinatorics"), "etaIPz", "simB");
  TH1* realZ = CopyH1(realCont, "ipz", "realZ");
  TH1* simZ  = CopyH1(simCont,  "ipz", "simZ");
  
  // Scale combinatorial background to measured to get beta
  simB->Divide(simM);
  simB->SetTitle("#beta'");
  simB->SetZTitle("#beta'");

  // Copy simulated beta to real beta, and scale by scalar
  TH2* realB = static_cast<TH2*>(simB->Clone("realB"));
  realB->SetDirectory(0);
  realB->SetTitle("#beta");
  realB->SetZTitle("#beta");
  realB->Multiply(scale);

  // Multiply real beta onto real measured to get background
  realC->Multiply(realB);
  realC->SetTitle("C");
  realC->SetZTitle("C");

  // Substract the real background off the real measured
  realS->Add(realC, -1);
  realS->SetTitle("S");
  realS->SetZTitle("S");
  
  // Substract combinatorial background from measured to get signal 
  simS->Add(simC, -1);
  simS->SetTitle("S'");
  simS->SetZTitle("S'");

  // Scale MC truth primaries by signal to get correction 
  simA->Divide(simS);
  simA->SetTitle("A'");
  simA->SetZTitle("#alpha'");

  // Make a fiducial distribution, and coerce the others to fit this
  TH2* fiducial = static_cast<TH2*>(simA->Clone("fiducial"));
  fiducial->SetTitle("F");
  fiducial->SetZTitle("F");
  for (Int_t i = 1; i <= fiducial->GetNbinsX(); i++) {
    for (Int_t j = 1; j <= fiducial->GetNbinsY(); j++) {
      Double_t c = fiducial->GetBinContent(i,j);
      fiducial->SetBinContent(i,j,c > fMinAlpha && c <= fMaxAlpha);
      fiducial->SetBinError(i,j,0);
    }
  }
  realM->Multiply(fiducial);
  realC->Multiply(fiducial);
  realS->Multiply(fiducial);
  realB->Multiply(fiducial);
  simM ->Multiply(fiducial);
  simC ->Multiply(fiducial);
  simS ->Multiply(fiducial);
  simB ->Multiply(fiducial);
  simA ->Multiply(fiducial);
	
  // We can now make our result
  TH2* result = static_cast<TH2*>(realS->Clone("result"));
  result->Multiply(simA);
  result->SetTitle("R");
  result->SetZTitle("R");

  // Output directory for full stuff
  TDirectory* full = outDir->mkdir("full");
  
  // Make a stack of 1D projections
  struct Rec {
    TH2* h;
    TH1* s;
    TH2* c;
    Style_t sty;
    Color_t col;
    Float_t siz;
    const char* tit;
  };
  Rec      sC = { simC,  simZ,  result, 32, kMagenta+2, 1.4, "background" };
  Rec      sS = { simS,  simZ,  result, 27, kGreen+2,   1.8, "signal" };
  Rec      sM = { simM,  simZ,  result, 26, kBlue+2,    1.4, "measured" };
  Rec      sG = { simG,  simZ,  0,      24, kRed+2,     1.4, "generated" };
  Rec      rC = { realC, realZ, result, 23, kMagenta+2, 1.2, "background" };
  Rec      rS = { realS, realZ, result, 33, kGreen+2,   1.6, "signal" };
  Rec      rM = { realM, realZ, result, 22, kBlue+2,    1.2, "measured" };
  Rec      rR = { result,realZ, 0,      20, kRed+2,     1.3, ObsTitle() };
  Rec*     recs[]  = { &rR, &sG, &rS, &sS, &rM, &sM, &rC, &sC, 0 };
  Rec**    ptr     = recs;
  TH1*     dndeta  = 0;
  THStack* all     = new THStack("all", "");
  if (fViz & kAltMarker) {
    rR.sty = 21;
    rR.siz = 1.2;
    sG.sty = 25;
    sG.siz = 1.3;
  }
  while ((*ptr)) {
    Rec*  src = *ptr;
    src->h->SetDirectory(full);
    TH1*  proj = AverageOverIPz(src->h, src->h->GetName(), 1,
				src->s, src->c);
    proj->SetYTitle(src->h->GetZaxis()->GetTitle());
    proj->SetTitle(Form("%s - %s", src->h->GetTitle(), src->tit));
    proj->SetDirectory(outDir);
    all->Add(proj);
    SetAttr(proj, src->col, src->sty, src->siz);
    if (src->h == result) {
      dndeta = proj;
      dndeta->SetYTitle(ObsTitle());
    }
    ptr++;
  }
  TF1* tmp = new TF1("mid", "pol0", -.5, +.5);
  dndeta->Fit(tmp, "Q0R+");
  TLatex* ltx = new TLatex(0,tmp->GetParameter(0)/2,
			   Form("%s|_{|#eta|<0.5}=%.1f#pm%.1f",
				ObsTitle(), 
				tmp->GetParameter(0),
				tmp->GetParError(0)));
  Printf("  %dD: %6.1f +/- %6.1f  (%5.2f)",
	 deltaDimen, tmp->GetParameter(0), tmp->GetParError(0),
	 tmp->GetChisquare()/tmp->GetNDF());
  
  ltx->SetTextAlign(22);
  ltx->SetTextFont(42);
  dndeta->GetListOfFunctions()->Add(ltx);
    
  outDir->cd();
  all->Write();
  realB   ->SetDirectory(full);
  simB    ->SetDirectory(full);
  simA    ->SetDirectory(full);
  fiducial->SetDirectory(full);
  
  outParent->cd();

  return dndeta;
}


//====================================================================
void AliTrackletdNdeta2::ClearCanvas()
{
  fTop->Clear();
  fTop->SetNumber(1);
  fTop->SetFillColor(kGray); // kTopBackground);
  fTop->SetFillStyle(1001);
  fTop->SetBorderSize(0);
  fTop->SetBorderMode(0);

  fBody->Clear();
  fBody->SetNumber(2);
  fBody->SetFillColor(0);
  fBody->SetFillStyle(0);
  fBody->SetBorderSize(0);
  fBody->SetBorderMode(0);
  fBody->SetTopMargin(0.01);
  fBody->SetRightMargin(0.01);
  fBody->SetBottomMargin(0.10);
  fBody->SetLeftMargin(0.10);
  fBody->SetTicks();
  
  fCanvas->cd();
}
//____________________________________________________________________
void AliTrackletdNdeta2::CreateCanvas(const TString& outputName)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  Int_t h    = 1000;
  Int_t w    = h / TMath::Sqrt(2);
  if (fViz & kLandscape) {
    Int_t t = h;
    h       = w;
    w       = t;
  }

  fCanvas = new TCanvas("c",outputName,w,h);
  fCanvas->SetFillColor(0);
  fCanvas->SetBorderSize(0);
  fCanvas->SetBorderMode(0);

  if (fViz & kPDF) {
    SuppressGuard g;    
    fCanvas->Print(Form("%s.pdf[", outputName.Data()),
		   Form("pdf %s", (fViz & kLandscape) ? "Landscape" : ""));
  }
  if (fViz & kPNG) {
    gSystem->mkdir(outputName);
  }
  fCanvas->SetLeftMargin  (0.10);
  fCanvas->SetRightMargin (0.05);
  fCanvas->SetTopMargin   (0.05);
  fCanvas->SetBottomMargin(0.10);
    
  Float_t dy = 0.05;
  fTop = new TPad("top","Top",0,1-dy,1,1,0,0);
  // fTop->SetNumber(1);
  // fTop->SetFillColor(kTopBackground);
  // fTop->SetBorderSize(0);
  // fTop->SetBorderMode(0);
  fCanvas->cd();
  fTop->Draw();

  fBody = new TPad("body","Body",0,0,1,1-dy,0,0);
  fBody->SetNumber(2);
  // fBody->SetFillColor(0);
  // fBody->SetFillStyle(0);
  // fBody->SetBorderSize(0);
  // fBody->SetBorderMode(0);
  fCanvas->cd();
  fBody->Draw();

  ClearCanvas();
  
  fHeader = new TLatex(.5, .5, "Title");
  fHeader->SetNDC();
  fHeader->SetTextAlign(22);
  // fHeader->SetTextColor(kWhite);
  fHeader->SetTextFont(62);
  fHeader->SetTextSize(0.7);
    
  fCanvas->cd();
}
//____________________________________________________________________
void AliTrackletdNdeta2::CloseCanvas()
{
  if ((fViz & kPDF) && fCanvas) {
    SuppressGuard g;
    fCanvas->Print(Form("%s.pdf]", fCanvas->GetTitle()),
		   Form("pdf %s Title:%s",
			(fViz & kLandscape) ? "Landscape" : "",
			fLastTitle.Data()));
    Printf("PDF %s written", fCanvas->GetTitle());
  }
  if (fCanvas) fCanvas->Close();
  fCanvas = 0;
}

//____________________________________________________________________
void AliTrackletdNdeta2::PrintCanvas(const char* title,
				     const char* shortTitle,
				     Float_t     size)
{
  if (fTop) {
    fTop->cd();
    fHeader->SetTextSize(size);
    fHeader->DrawLatex(.5,.5,title);
  }
  fCanvas->Modified();
  fCanvas->Update();
  fCanvas->cd();

  if (fViz & kPDF) {
    TString tit;
    tit.Form("pdf %s Title:%s", (fViz & kLandscape) ? "Landscape" : "",
	     title);
    // Suppress prints
    SuppressGuard g;
    fCanvas->Print(Form("%s.pdf",fCanvas->GetTitle()), tit);
  }
  static Int_t cnt = 1;
  if (fViz & kPNG) {
    SuppressGuard g;
    fCanvas->Print(Form("%s/%03d_%s.png",fCanvas->GetTitle(),cnt,shortTitle));
    cnt++;
  }
  fLastTitle = title;
  if (fViz & kPause) fCanvas->WaitPrimitive();
}
//____________________________________________________________________
TLegend* AliTrackletdNdeta2::DrawInPad(TVirtualPad* c,
				       Int_t        pad,
				       TObject*     o,
				       Option_t*    opt)
{
  if (!o) {
    // Warning("", "Nothing to draw in pad %d", pad);
    return 0;
  }
  TLegend*     l = 0;
  TVirtualPad* p = c->cd(pad);
  TString option(opt);
  option.ToLower();
  if (option.Contains("logx")) { p->SetLogx(); option.ReplaceAll("logx",""); }
  if (option.Contains("logy")) { p->SetLogy(); option.ReplaceAll("logy",""); }
  if (option.Contains("logz")) { p->SetLogz(); option.ReplaceAll("logz",""); }
  if (option.Contains("grid")) { p->SetGridx(); p->SetGridy();
    option.ReplaceAll("grid",""); }
  Int_t leg = 0;
  if (option.Contains("leg3")) { leg = 3; option.ReplaceAll("leg3",""); }
  if (option.Contains("leg2")) { leg = 2; option.ReplaceAll("leg2",""); }
  if (option.Contains("leg")) {  leg = 1; option.ReplaceAll("leg",""); }
  // Printf("Drawing %p %s with %s", o, o->GetName(), option.Data());
  o->Draw(option);
  if (leg) {
    l = p->BuildLegend(0.5, 0.73, .98, .98);
    l->SetNColumns(leg);
    TObject* frame = 0;
    TIter next(l->GetListOfPrimitives());
    TLegendEntry* ent = 0;
    while ((ent = static_cast<TLegendEntry*>(next()))) {
      if (TString(ent->GetLabel()).EqualTo("frame")) frame = ent;
    }
    if (frame) l->GetListOfPrimitives()->Remove(frame);
    // l->GetListOfPrimitives()->Print();
  }
  p->Modified();
  // p->Update();
  // p->cd();
  return l;
}
//____________________________________________________________________
void AliTrackletdNdeta2::ModLegend(TVirtualPad* p, TLegend* l,
				   Double_t x1, Double_t y1,
				   Double_t x2, Double_t y2)
{
  Double_t px1 = p->GetX1();
  Double_t px2 = p->GetX2();
  Double_t py1 = p->GetY1();
  Double_t py2 = p->GetY2();
  l->SetX1(px1+(px2-px1)*x1);
  l->SetX2(px1+(px2-px1)*x2);
  l->SetY1(py1+(py2-py1)*y1);
  l->SetY2(py1+(py2-py1)*y2);
  p->Modified();
}
//____________________________________________________________________
THStack* AliTrackletdNdeta2::Make2Stack(const char*      name,
					const char*      title,
					Container*       realList,
					Container*       simList,
					Option_t*        realOpt,
					Option_t*        simOpt)
{
  TString  nme(name);
  THStack* stack = new THStack(name, title);
  TH1*     real  = CopyH1(realList, name, Form("real%s",name));
  TH1*     sim   = CopyH1(simList,  name, Form("sim%s",name));
  real->SetMarkerStyle(20);
  sim ->SetMarkerStyle(24);
  real->SetFillStyle(3004);
  sim ->SetFillStyle(3005);
  real->SetBarWidth(0.4);
  sim ->SetBarWidth(0.4);
  real->SetBarOffset(0.1);
  sim ->SetBarOffset(0.5);
  TString dtit(real->GetTitle());
  if (dtit.Contains("\\")) dtit.Form("%s\\hbox{ - real}", real->GetTitle());
  else                     dtit.Form("%s - real", real->GetTitle());
  real->SetTitle(dtit);
  TString stit(sim->GetTitle());
  if (stit.Contains("\\")) stit.Form("%s\\hbox{ - sim.}", sim->GetTitle());
  else                     stit.Form("%s - sim.", sim->GetTitle());
  sim->SetTitle(stit);
  stack->Add(real, realOpt);
  stack->Add(sim,  simOpt);
  return stack;
}


//====================================================================
Bool_t AliTrackletdNdeta2::Visualize(Container*  realSums,
				     Container*  simSums,
				     Container*  realRess,
				     Container*  simRess, 
				     TDirectory* outDir,
				     Int_t       maxBins)
{
  // --- Visualization -----------------------------------------------
  TH1* realCent = static_cast<TH1*>(GetO(outDir, "realCent"));
  TString outName(outDir->GetName());
  outName.ReplaceAll(".root", "");
  CreateCanvas(outName);
  VisualizeParams(realSums, simSums);
  VisualizeGeneral(realRess, simRess);
  VisualizeWeights(simRess);
  for (Int_t i = 0; i < 4; i++) {
    if ((fProc & (1 << i)) == 0) continue;
    VisualizeFinal(outDir, i);
  }
  for (Int_t i = 1; i <= realCent->GetNbinsX() && i <= maxBins ; i++) {
    Double_t c1 = realCent->GetXaxis()->GetBinLowEdge(i);
    Double_t c2 = realCent->GetXaxis()->GetBinUpEdge (i);
      
    VisualizeBin(c1, c2, simRess, outDir);
  }
  CloseCanvas();
}

//====================================================================
void AliTrackletdNdeta2::VisualizeGeneral(Container* realList,
					  Container* simList)
{
  THStack*   ipz    = Make2Stack("ipz",   "IP_{#it{z}}",      realList,simList);
  THStack*   cent   = Make2Stack("cent",  "Centrality [%]",   realList,simList);
  THStack*   status = Make2Stack("status","Task status",      realList,simList,
				 "B text90", "B text90");
  THStack*   centTr = Make2Stack("centTracklets", "#LTtracklets#GT",
				 realList, simList, "E", "E");
  ClearCanvas();
  fBody->Divide(1,3);
  for (Int_t i = 1; i <= 3; i++) {
    if (i < 3) fBody->GetPad(i)->SetRightMargin(0.01);
    fBody->GetPad(i)->SetTopMargin(0.01);
  }
  TLegend* l = 0;
  TVirtualPad* q = fBody->GetPad(1);
  q->Divide(2,1);
  l = DrawInPad(q,1,ipz,     "nostack leg");
  ModLegend(q->GetPad(1),l,.4,.1,.75,.4);
  l = DrawInPad(q,2,cent,    "nostack leg");
  ModLegend(q->GetPad(2),l,.6,.1,.99,.4);
  q = fBody->GetPad(2);
  q->Divide(2,1);
  l = DrawInPad(q,1,status,  "nostack hist text90 leg");
  ModLegend(q->GetPad(1),l,.5,.7,.99,.99);
  l = DrawInPad(q,2,centTr,  "nostack leg");
  ModLegend(q->GetPad(2),l,.5,.7,.99,.99);

  TH2* real = GetH2(realList, "etaPhi");
  TH2* sim  = GetH2(simList,  "etaPhi");
  if (sim) {
    sim->SetMarkerColor(kBlack);
    sim->SetMarkerStyle(0);
    sim->SetMarkerSize(1);
    sim->SetLineColor(kBlack);
    sim->SetFillColor(kBlack);
    sim->SetFillStyle(0);
    sim->SetName("simEtaPhi");
  }
  DrawInPad(fBody, 3, real, "colz"); 
  DrawInPad(fBody, 3, sim,  "box same");
  
  PrintCanvas("General information","general");
}

namespace {
  void SetCentColors(THStack* s, TH1* dist=0)
  {
    if (!s->GetHists()) return;

    const Color_t cc[] = { kMagenta+2, // 0
			   kBlue+2,    // 1
			   kAzure-1,   // 2 // 10,
			   kCyan+2,    // 3
			   kGreen+1,   // 4 
			   kSpring+5,  // 5 //+10,
			   kYellow+1,  // 6
			   kOrange+5,  // 7 //+10,
			   kRed+1,     // 8
			   kPink+5,    // 9 //+10,
			   kBlack };   // 10    
    TIter next(s->GetHists());
    TH1*  h = 0;
    Int_t i = 0;
    Double_t min = +10000;
    Double_t max = -10000;
    while ((h = static_cast<TH1*>(next()))) {
      Color_t c = cc[i % 10];
      h->SetMarkerColor(c);
      h->SetFillColor(c);
      h->SetLineColor(c);
      h->SetMarkerStyle(20+(i%4));
      h->SetDirectory(0);
      min = TMath::Min(h->GetMinimum(), min);
      max = TMath::Max(h->GetMaximum(), max);
      i++;
      if (!dist) continue;
      h->SetTitle(Form("%5.1f-%5.1f%%",
		       dist->GetXaxis()->GetBinLowEdge(i),
		       dist->GetXaxis()->GetBinUpEdge(i)));
		       
    }
    s->SetMinimum(min*.9);
    s->SetMaximum(max*1.1);
  }
  THStack* GetPdgStack(AliTrackletAODUtils::Container* w, const char* name)
  {
    AliTrackletAODUtils::Container* c = AliTrackletAODUtils::GetC(w, name);
    if (!c) return 0;

    THStack* s = new THStack(name, "");
    TH1*     h = 0;
    TIter    n(c);
    while ((h = static_cast<TH1*>(n()))) s->Add(h);

    return s;
  }
}
//____________________________________________________________________
void AliTrackletdNdeta2::VisualizeWeights(Container* simList)
{
  TH1*       c = GetH1(simList, "cent");
  Container* w = GetC(simList,"weights", false);
  if (!w) return;

  Float_t right = .75;
  ClearCanvas();
  fBody->SetTopMargin(0.01);
  fBody->Divide(1,3);
  TVirtualPad* pp[] = { fBody->GetPad(1), fBody->GetPad(2) };
  for (Int_t i = 0; i < 2; i++) {
    pp[i]->SetRightMargin(0.01);
    pp[i]->SetTopMargin(0.01);
    pp[i]->SetPad(pp[i]->GetXlowNDC(), pp[i]->GetYlowNDC(),
		  right, pp[i]->GetYlowNDC()+pp[i]->GetHNDC());
    pp[i]->Modified();
  }
  
  THStack* ef = new THStack(GetP2(simList,"etaWeight"),"x","effWeights","");
  THStack* pt = new THStack(GetH2(w,      "centPt"),   "y","pt","");
  THStack* ab = GetPdgStack(w, "abundance");
  THStack* st = GetPdgStack(w, "strangeness");
  SetCentColors(ef, c);
  SetCentColors(pt);
  SetCentColors(ab);
  SetCentColors(st);

  ef->SetMinimum(0.98);
  ef->SetMaximum(1.02);
  TLegend* l = DrawInPad(fBody, 1, ef, "nostack leg");
  DrawInPad(fBody, 2, pt, "nostack");
  ef->GetHistogram()->SetYTitle("Average weight");
  ef->GetHistogram()->SetXTitle("#eta");
  pt->GetHistogram()->SetYTitle("Weight");
  pt->GetHistogram()->SetXTitle("#it{p}_{T}");

  fBody->GetPad(1)->GetListOfPrimitives()->Remove(l);
  fBody->GetPad(1)->Modified();

  fBody->cd();
  l->Draw();
  ModLegend(fBody, l, right, pp[1]->GetYlowNDC(),
	    .99, 1-fBody->GetTopMargin());
  fBody->Modified();
  
  TVirtualPad* p3 = fBody->GetPad(3);
  p3->SetTopMargin(0.01);
  p3->SetRightMargin(0.01);
  p3->Divide(2,1);
  p3->GetPad(1)->SetRightMargin(0.01);
  p3->GetPad(2)->SetRightMargin(0.01);
  p3->GetPad(1)->SetTopMargin(0.01);
  p3->GetPad(2)->SetTopMargin(0.01);

  DrawInPad(p3, 1, ab, "nostack leg");
  DrawInPad(p3, 2, st, "nostack leg");

  ab->GetHistogram()->SetYTitle("Weight");
  ab->GetHistogram()->SetXTitle("Centrality [%]");
  st->GetHistogram()->SetYTitle("Weight");
  st->GetHistogram()->SetXTitle("Centrality [%]");

  p3->GetPad(1)->Modified();
  p3->GetPad(2)->Modified(); 
  
  PrintCanvas("Simulation weights","weights");
}

//____________________________________________________________________
void AliTrackletdNdeta2::VisualizeFinal(TDirectory* outDir, Int_t i)
{
  TDirectory* dd  = outDir->GetDirectory(Form("final%dd", i));
  if (!dd) return;

  THStack* all = static_cast<THStack*>(GetO(dd, "full"));
  TH1*     mid = static_cast<TH1*>    (GetO(dd, "mid"));
  TH1*     pub = static_cast<TH1*>    (GetO(outDir, "published"));
  Double_t max = TMath::Max(all->GetMaximum("nostack"),mid->GetMaximum());
  Double_t min = TMath::Min(all->GetMinimum("nostack"),mid->GetMinimum());

  all->SetMinimum(.9*min);
  mid->SetMinimum(.9*min);
  all->SetMaximum(1.2*max);
  mid->SetMaximum(1.2*max);
  mid->SetLineColor(kBlack);
  mid->SetFillColor(kBlack);
  mid->SetMarkerColor(kBlack);
  
  ClearCanvas();

  TPad* p1 = new TPad("p1","p1",0,0,.4,1);
  p1->SetTopMargin(0.01);
  p1->SetRightMargin(0.0);
  p1->SetLeftMargin(0.12);
  p1->SetBottomMargin(0.15);
  p1->SetTicks();
  fBody->cd();
  p1->Draw();
  p1->SetNumber(2);

  Double_t right = .7;
  TPad* p2 = new TPad("p2","p2",.4,0,1,1);
  p2->SetTopMargin(0.01);
  p2->SetRightMargin(1-right);
  p2->SetLeftMargin(0.0);
  p2->SetBottomMargin(0.15);
  p2->SetTicks();
  fBody->cd();
  p2->Draw();
  p2->SetNumber(2);

  
  DrawInPad(p1,0, mid, "logy grid");
  DrawInPad(p1,0, pub, "logy grid same");
  TLegend* l = DrawInPad(p2,0, all, "nostack logy grid leg");
  all->GetHistogram()->SetXTitle("#eta");

  ModLegend(p2, l, right, .15, .99, .99);
  l->SetMargin(0.2);
  l->SetEntrySeparation(0.1);
  l->SetTextSize(0.04);
  p1->Modified();
  p2->Modified();
  fBody->Modified();
  
  const char* what = (i == 3 ? "d^{3}N/(d#Deltad#etadIP_{z})" :
		      i == 2 ? "d^{2}N/(d#Deltad#eta)" :
		      i == 1 ? "dN/d#Delta" : "dN/d#Delta (k#equiv1)");
  PrintCanvas(Form("Results #topbar %s", what), "results");
}

//====================================================================
void AliTrackletdNdeta2::VisualizeParam(const char* name,
				  Double_t&   y,
				  const char* val)
{
  TLatex* ln = new TLatex(.49, y, name);
  ln->SetTextAlign(31);
  ln->SetTextSize(0.02/gPad->GetHNDC());
  ln->SetNDC();
  ln->Draw();
  TLatex* lv = new TLatex(.51, y, val);
  lv->SetTextAlign(11);
  lv->SetTextSize(0.02/gPad->GetHNDC());
  lv->SetNDC();
  lv->Draw();
  y -= 0.025/gPad->GetHNDC();
}
//____________________________________________________________________  
void AliTrackletdNdeta2::VisualizeParam(const char* name,
					Double_t&   y,
					Double_t    val)
{
  VisualizeParam(name, y, Form("%f", val));
}
//____________________________________________________________________  
void AliTrackletdNdeta2::VisualizeParam(const char* name,
					Double_t&   y,
					Int_t       val)
{
  VisualizeParam(name, y, Form("%d", val));
}
//____________________________________________________________________  
void AliTrackletdNdeta2::VisualizeParam(const char* name,
					Double_t&   y,
					Bool_t      val)
{
  VisualizeParam(name, y, val ? "yes" : "no");
}
//____________________________________________________________________
void AliTrackletdNdeta2::VisualizeParams(Container*  pars,
					 const char* title)
{
  // c->Clear();
  Double_t y = .9;
  TLatex* latex = new TLatex(.5, y, title);
  latex->SetTextAlign(21);
  latex->SetTextSize(0.023/gPad->GetHNDC());
  latex->SetNDC();
  latex->Draw();
  y -= 0.028/gPad->GetHNDC();
  if (!pars) return;
  VisualizeParam("#delta#phi shift",      y, GetD(pars, "DPhiShift"));
  VisualizeParam("Shifted #delta#phi cut",y, GetD(pars, "ShiftedDPhiCut"));
  VisualizeParam("#Delta cut",            y, GetD(pars, "DeltaCut"));
  VisualizeParam("max#Delta",             y, GetD(pars, "MaxDelta"));
  VisualizeParam("min#Delta_{tail}",      y, GetD(pars,"TailDelta"));
  VisualizeParam("abs.min#it{c}",         y, GetD(pars,"AbsMinCent"));
}
//____________________________________________________________________
void AliTrackletdNdeta2::VisualizeParams(Container* realSums,
					 Container* simSums)
{
  ClearCanvas();
  fBody->Divide(1,3,0,0);

  Double_t yr = .1;
  TVirtualPad* p1 = fBody->GetPad(1);  p1->SetPad(0,1-yr,1,1);
  TVirtualPad* p2 = fBody->GetPad(2);  p2->SetPad(0,(1-yr)/2,1,1-yr);
  TVirtualPad* p3 = fBody->GetPad(3);  p3->SetPad(0,0,1,(1-yr)/2);
  
  // Post-processing stuff 
  fBody->cd(1);
  Double_t y = .80;
  TLatex* latex = new TLatex(.5, y, "Post-processing");
  latex->SetTextAlign(21);
  latex->SetTextSize(0.023/p1->GetHNDC());
  latex->SetNDC();
  latex->Draw();
  y -= 0.028/p1->GetHNDC();
  // TString scaleM("none");
  // if      (fProc & kScaleNone)    scaleM = "none";
  // else if (fProc & kScaleFix)     scaleM.Form("%4.2f", fCombinatoricsScale);
  // else if (fProc & kScaleAverage) scaleM = "average";
  // else if (fProc & kScaleFull)    scaleM = "full";
  // if      (fProc & kScaleDouble)  scaleM.Append(" (double)");  
  // VisualizeParam("Scaling of comb. bg.", y, scaleM);	     
  VisualizeParam("min#alpha", y, fMinAlpha);
  VisualizeParam("max#alpha", y, fMaxAlpha);
  
  // From tasks 
  fBody->cd(2);
  VisualizeParams(GetC(realSums, "parameters"), "Real data");
  fBody->cd(3);
  VisualizeParams(GetC(simSums, "parameters"), "Simulated data");
  PrintCanvas("Parameters", "parameters");
}

//====================================================================
Bool_t AliTrackletdNdeta2::VisualizeBin(Double_t    c1,
					Double_t    c2,
					Container*  simList, 
					TDirectory* outTop)
{
  // Form the folder name
  TString centName(CentName(c1,c2));
  fLastBin.Form("%.1f#minus%.1f%%", c1, c2);
  fLastShort = centName;
  
  TDirectory* outDir = outTop->GetDirectory(centName);
  if (!outDir) {
    Warning("VisualizeBin", "Directory %s not found in %s",
	    centName.Data(), outTop->GetName());
    return false;
  }

  Printf("%5.1f - %5.1f%%", c1, c2);
  VisualizeSpecies(GetC(simList, centName));
  VisualizePrimary(GetC(GetC(simList, centName), "generated"));
  for (Int_t i = 0; i < 4; i++) {
    if ((fProc & (1 << i)) == 0) continue;
    VisualizeDelta(outDir, i);
    VisualizeResult(outDir, i);
  }
  
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletdNdeta2::VisualizeSpecies(Container* simCont)
{
  if (!simCont) return true;
  ClearCanvas();
  Container* meas = GetC(simCont, "measured");
  Container* comb = GetC(simCont, "combinatorics");
  Container* prim = GetC(simCont, "primaries");
  Container* seco = GetC(simCont, "secondaries");
  Container* cont[] = { meas, comb, prim, seco };
  const char* tit[] = { "M' by primary mother's specie",
			"C' by primary mother's specie",
			"P' by species",
			"S' by primary mother's specie" };
  fBody->SetTopMargin(0.01);
  fBody->SetRightMargin(0.01);
  fBody->Divide(2,2,0.01,0.01);

  for (Int_t i = 0; i < 4; i++) {
    if (!cont[i]) continue;

    Container* species = GetC(cont[i], "types");
    if (!species) continue;

    THStack* all    = static_cast<THStack*>(GetO(species, "all"));
    THStack* toPion = static_cast<THStack*>(GetO(species, "toPion"));
    
    
    TVirtualPad* p = fBody->GetPad(i+1);
    p->SetTopMargin(0.10);
    p->SetRightMargin(0.15);
    p->SetBottomMargin(0.15);
    p->cd();
    TLatex* ltx = new TLatex(.5, .99, tit[i]);
    ltx->SetTextAlign(23);
    ltx->SetTextSize(0.04);
    ltx->SetTextFont(42);
    ltx->Draw();
    
    p->Divide(1,2,0,0);
    TLegend* l = DrawInPad(p, 1, all, "nostack grid leg");
    DrawInPad(p, 2, toPion, "nostack grid");
    all->GetHistogram()->SetYTitle("dN_{X}/d#eta");
    all->GetHistogram()->SetXTitle("#eta");
    all->GetHistogram()->GetYaxis()->SetTitleSize(0.08);
    all->GetHistogram()->GetYaxis()->SetTitleOffset(0.6);
    all->GetHistogram()->GetYaxis()->SetLabelSize(0.08);
    toPion->GetHistogram()->SetYTitle("Relative to #pi^{#pm} mothers");
    toPion->GetHistogram()->SetXTitle("#eta");
    toPion->GetHistogram()->GetYaxis()->SetTitleSize(0.08);
    toPion->GetHistogram()->GetYaxis()->SetTitleOffset(0.6);
    toPion->GetHistogram()->GetYaxis()->SetLabelSize(0.08);
    toPion->GetHistogram()->GetXaxis()->SetTitleSize(0.08);
    toPion->GetHistogram()->GetXaxis()->SetTitleOffset(0.6);
    toPion->GetHistogram()->GetXaxis()->SetLabelSize(0.08);
    p->GetPad(1)->GetListOfPrimitives()->Remove(l);
    p->GetPad(1)->Modified();
    p->GetPad(2)->Modified();
    p->cd();
    l->Draw();
    
    ModLegend(p, l, .85, p->GetBottomMargin(), .99, 1-p->GetTopMargin());
    p->Modified();
  }
  
  PrintCanvas(Form("Species #topbar %s", fLastBin.Data()),
	      Form("%s_species", fLastShort.Data()));
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletdNdeta2::VisualizePrimary(Container* simCont)
{
  if (!simCont) return true;
  ClearCanvas();

  TH2* etaPdg = GetH2(simCont, "etaPdg");
  TH2* etaPt  = GetH2(simCont, "etaPt");


  TH1*      pion = 0;
  THStack*  all     = new THStack("all",    "All primaries");
  THStack*  ratios  = new THStack("ratios", "Ratios to pions, all primaries");
  TAxis*    pdgAxis = etaPdg->GetYaxis();
  
  for (Int_t i = 1; i <= pdgAxis->GetNbins(); i++) {
    TString lbl  = pdgAxis->GetBinLabel(i);
    Int_t   apdg = lbl.Atoi();
    switch (apdg) {
    case 22: continue; // ignore photons
    }
    TString prty;
    Color_t colo;
    Style_t styl;
    PdgAttr(apdg, prty, colo, styl);
    TH1*    proj = static_cast<TH1*>(etaPdg->ProjectionX(Form("h%d",apdg),
							 i, i));
    proj->SetTitle(prty);
    proj->SetMarkerColor(colo);
    proj->SetFillColor(colo);
    proj->SetMarkerColor(colo);
    proj->SetLineColor(colo);
    proj->SetMarkerStyle(styl);
    proj->SetFillStyle(0);
    proj->SetFillColor(0);
    switch (apdg) {
    case 321:  proj->SetBinContent(0, 0.15);   break; // PRC88,044910
    case 2212: proj->SetBinContent(0, 0.05);   break; // PRC88,044910
    case 310:  proj->SetBinContent(0, 0.075);  break; // PRL111,222301
    case 3122: proj->SetBinContent(0, 0.018);  break; // PRL111,222301
    case 3212: proj->SetBinContent(0, 0.0055); break; // NPA904,539
    case 3322: proj->SetBinContent(0, 0.005);  break; // PLB734,409
    case 211:  proj->SetBinContent(0, 1);      break; // it self 
    default:   proj->SetBinContent(0, -1);     break; // Unknown
    }
    if (apdg == 211) pion = proj;
    all->Add(proj);
  }
    
  TIter next(all->GetHists());
  TH1*  proj = 0;
  Double_t rmin = +1000000000;
  Double_t rmax = -1000000000;
  while ((proj = static_cast<TH1*>(next()))) {
    if (proj == pion) continue;
    Double_t r276 = proj->GetBinContent(0);
    if (r276 < 0 || r276 >= 1) continue;
    TH1*     copy = static_cast<TH1*>(proj->Clone(Form("r%s",proj->GetName())));
    copy->Divide(pion);
    copy->SetFillStyle(0);
    copy->SetFillColor(0);
    copy->SetTitle(Form("%s / %s", copy->GetTitle(), pion->GetTitle()));
    ratios->Add(copy);

    if (r276 > 0 && r276 < 1) {
      TGraphErrors* g = new TGraphErrors(1);
      g->SetName(Form("%s_2760", copy->GetName()));
      g->SetTitle(Form("%s in #sqrt{s_{NN}}=2.76TeV", copy->GetTitle()));
      g->SetPoint(0,0,r276);
      g->SetPointError(0,.5,0);
      g->SetLineColor(proj->GetLineColor());
      g->SetLineStyle(proj->GetLineStyle());
      g->SetMarkerColor(proj->GetMarkerColor());
      g->SetMarkerStyle(proj->GetMarkerStyle());
      g->SetMarkerSize(1.2*proj->GetMarkerSize());
      copy->GetListOfFunctions()->Add(g,"p");
      copy->SetMaximum(TMath::Max(copy->GetMaximum(),r276));
      copy->SetMinimum(TMath::Min(copy->GetMinimum(),r276));
      rmin = TMath::Min(copy->GetMinimum(),rmin);
      rmax = TMath::Max(copy->GetMaximum(),rmax);
    }
  }
  // ratios->SetMinimum(rmin);
  ratios->SetMaximum(rmax*1.1);

  Int_t etaM = etaPt->GetXaxis()->FindBin(-.5);
  Int_t etaP = etaPt->GetXaxis()->FindBin(+.5);
  TH1*  pt   = etaPt->ProjectionY("pt", etaM, etaP);
  pt->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T}");
  pt->GetYaxis()->SetTitleSize(0.08);
  pt->GetYaxis()->SetLabelSize(0.08);
  pt->GetYaxis()->SetTitleOffset(0.6);  
  pt->GetXaxis()->SetTitleSize(0.08);
  pt->GetXaxis()->SetLabelSize(0.08);
  pt->GetXaxis()->SetTitleOffset(0.6);  

  TPad* p1 = new TPad("p1","p1",0,.3,1,1);
  p1->SetTopMargin(.01);
  p1->SetRightMargin(.01);
  fBody->cd();
  p1->Draw(); 
  p1->SetNumber(1);
  p1->Divide(1,2,0,0);
  p1->GetPad(1)->SetRightMargin(0.2);
  p1->GetPad(2)->SetRightMargin(0.2);
  TLegend* l = DrawInPad(p1,1,all,    "leg2 nostack logy");
  all->GetHistogram()->GetYaxis()->SetTitle("d#it{N}_{X}/d#eta");
  all->GetHistogram()->GetYaxis()->SetTitleSize(0.08);
  all->GetHistogram()->GetYaxis()->SetLabelSize(0.08);
  all->GetHistogram()->GetYaxis()->SetTitleOffset(0.6);  
  ModLegend(p1->GetPad(1), l,
	    1-p1->GetPad(1)->GetRightMargin(),
	    p1->GetPad(1)->GetBottomMargin(),
	    1-p1->GetPad(1)->GetTopMargin(),
	    .99);
  l->SetBorderSize(0);
  l = DrawInPad(p1,2,ratios, "nostack leg");
  ratios->GetHistogram()->GetYaxis()->SetTitle("Ratios to #pi");
  ratios->GetHistogram()->GetYaxis()->SetTitleSize(0.08);
  ratios->GetHistogram()->GetYaxis()->SetLabelSize(0.08);
  ratios->GetHistogram()->GetYaxis()->SetTitleOffset(0.6);  
  ratios->GetHistogram()->GetXaxis()->SetTitle("#eta");
  ratios->GetHistogram()->GetXaxis()->SetTitleSize(0.08);
  ratios->GetHistogram()->GetXaxis()->SetLabelSize(0.08);
  ratios->GetHistogram()->GetXaxis()->SetTitleOffset(0.6);  
  ModLegend(p1->GetPad(2), l,
	    1-p1->GetPad(2)->GetRightMargin(),
	    p1->GetPad(2)->GetBottomMargin(),
	    1-p1->GetPad(2)->GetTopMargin(),
	    .99);
  l->SetBorderSize(0);
  p1->Modified();
  p1->Update();
  p1->cd();
  
  TPad* p2 = new TPad("p2","p2",0,0,1,.3);
  p2->SetTopMargin(.01);
  p2->SetRightMargin(.01);
  p2->SetBottomMargin(0.15);
  fBody->cd();
  p2->Draw(); 
  p2->SetNumber(2);
  DrawInPad(p2,0,pt,     "logx logy");
  p2->Modified();
  p2->Update();
  p2->cd();

  fBody->Modified();
  
  PrintCanvas(Form("Primary species #topbar %s", fLastBin.Data()),
	      Form("%s_primary_species", fLastShort.Data()));
		
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletdNdeta2::VisualizeDelta(TDirectory* outTop,
					  Int_t       dimen)
{
  TDirectory* outDir = outTop->GetDirectory(Form("delta%dd", dimen));
  if (!outDir) {
    Warning("VisualizeDelta", "Directory detla%dd not found in %s",
	    dimen, outTop->GetName());
    return false;
  }

  ClearCanvas();
  fBody->cd();
  TPad* pq = new TPad("p1","p1",0, .3, 1, 1);
  pq->SetNumber(1);
  pq->Draw();

  TVirtualPad* q = fBody->cd(1);
  q->SetTopMargin(0.01);
  q->SetRightMargin(0.01);
  q->Divide(1,2,0,0);
  q->GetPad(1)->SetRightMargin(0.15);
  q->GetPad(2)->SetRightMargin(0.15);
  TVirtualPad* qq = q->GetPad(1);
  THStack* all = static_cast<THStack*>(GetO(outDir,"all"));
  TLegend* l = DrawInPad(q,1,all,"nostack logx logy grid leg");
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetMargin(0.2);  
  l->SetEntrySeparation(0.1);
  all->GetHistogram()->GetYaxis()->SetTitle("d#it{N}/d#Delta");
  all->GetHistogram()->GetYaxis()->SetLabelSize(0.06);
  all->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  all->GetHistogram()->GetYaxis()->SetTitleOffset(0.6);
  ModLegend(qq, l, 1-qq->GetRightMargin(), qq->GetBottomMargin(), .99,
	    1-qq->GetTopMargin()-.01);

  
  THStack* ratios = static_cast<THStack*>(GetO(outDir,"ratios"));
  ratios->SetMinimum(.6);
  ratios->SetMaximum(1.4);
  qq = q->GetPad(2);
  l = DrawInPad(q,2,ratios,"nostack logx grid leg");
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetMargin(0.2);
  l->SetEntrySeparation(0.1);
  ratios->GetHistogram()->GetXaxis()->SetTitle("#Delta");
  ratios->GetHistogram()->GetYaxis()->SetTitle("Ratio");
  ratios->GetHistogram()->GetYaxis()->SetLabelSize(0.06);
  ratios->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
  ratios->GetHistogram()->GetYaxis()->SetTitleOffset(0.6);
  ratios->GetHistogram()->GetXaxis()->SetLabelSize(0.06);
  ratios->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
  ratios->GetHistogram()->GetXaxis()->SetTitleOffset(0.6);
  
  ModLegend(qq, l, 1-qq->GetRightMargin(), qq->GetBottomMargin(), .99,
	    1-qq->GetTopMargin()-.01);

  fBody->cd();
  pq = new TPad("p2","p2",0, 0, 1, .3);
  pq->SetBottomMargin(0.25);
  pq->SetNumber(2);
  pq->Draw();

  q = fBody->cd(2);
  q->SetTopMargin(0.01);
  q->SetRightMargin(0.01);
  q->Divide(1,2,0,0);
  q->GetPad(1)->SetRightMargin(0.10);
  q->GetPad(2)->SetRightMargin(0.10);
  TH2* scale     = static_cast<TH2*>(GetO(outDir,"scale"));
  TH1* scaleProj = static_cast<TH2*>(GetO(outDir,"scaleProj"));
  DrawInPad(q,1,scale,    "colz");
  DrawInPad(q,2,scaleProj,"");
  scale->SetYTitle("IP_{#it{z}}");
  scaleProj->SetYTitle("#it{k}");
  scaleProj->SetXTitle("#eta");
  scale->GetYaxis()->SetLabelSize(0.12);
  scale->GetYaxis()->SetTitleSize(0.12);
  scale->GetYaxis()->SetTitleOffset(0.4);
  scaleProj->GetYaxis()->SetLabelSize(0.12);
  scaleProj->GetYaxis()->SetTitleSize(0.12);
  scaleProj->GetYaxis()->SetTitleOffset(0.4);
  scaleProj->GetYaxis()->SetNdivisions(207);
  scaleProj->GetXaxis()->SetLabelSize(0.12);
  scaleProj->GetXaxis()->SetTitleSize(0.12);
  scaleProj->GetXaxis()->SetTitleOffset(0.6);
  
  const char* what = (dimen == 3 ? "d^{3}N/(d#Deltad#etadIP_{z})" :
		      dimen == 2 ? "d^{2}N/(d#Deltad#eta)" :
		      dimen == 1 ? "dN/d#Delta" : "dN/d#Delta (k#equiv1)");
  PrintCanvas(Form("%s #topbar %s", what, fLastBin.Data()),
	      Form("%s_deltas", fLastShort.Data()));

  return true;
}

//____________________________________________________________________
Bool_t AliTrackletdNdeta2::VisualizeResult(TDirectory* outTop,
					   Int_t       dimen)
{
  TDirectory* outDir = outTop->GetDirectory(Form("results%dd", dimen));
  if (!outDir) {
    Warning("VisualizeDelta", "Directory results%dd not found in %s",
	    dimen, outTop->GetName());
    return false;
  }

  ClearCanvas();
  Double_t yr = .2;
  fBody->SetTopMargin(yr);
  fBody->SetRightMargin(0.01);
  fBody->SetLeftMargin(0.15);

  THStack* all = static_cast<THStack*>(GetO(outDir, "all"));
  TLegend* l = DrawInPad(fBody,0,all, "nostack leg2");
  all->GetHistogram()->SetXTitle("#eta");
  all->GetHistogram()->SetYTitle(ObsTitle());
  all->GetHistogram()->GetYaxis()->SetTitleOffset(1.7);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetMargin(0.12);
  l->SetEntrySeparation(0.1);
  l->SetHeader("R=G'/(M'-C')#times(M-C) #kern[2]{ } [C=kC'/M'#timesM]");
  // l->SetTextSize(0.04);
  // l->SetTextAlign(12);
  ModLegend(fBody, l, fBody->GetLeftMargin()-.01, 1-yr,
	    1-fBody->GetRightMargin(), .99);
  fBody->Modified();
  // fBody->GetListOfPrimitives()->Print();
  
  const char* what = (dimen == 3 ? "k(#eta,IP_{z})" :
		      dimen == 2 ? "k(#eta)" :
		      dimen == 1 ? "k=const." : "k#equiv1");
  PrintCanvas(Form("%s [%s] #topbar %s", ObsTitle(), what, fLastBin.Data()),
	      Form("%s_summary", fLastShort.Data()));

  return true;
}

//====================================================================
TFile* AliTrackletdNdeta2::OpenFile(const char* filename)
{
  TFile* file = TFile::Open(filename, "READ");
  if (!file) {
    Warning("OpenFile", "Failed to open \"%s\"", filename);
    return 0;
  }
  return file;
}
//____________________________________________________________________
const char* AliTrackletdNdeta2::CentName(Double_t c1, Double_t c2)
{
  static TString tmp;
  tmp.Form("cent%06.2f_%06.2f", c1, c2);
  tmp.ReplaceAll(".", "d");
  return tmp.Data();
}
  
//____________________________________________________________________
TH1* AliTrackletdNdeta2::SetAttr(TH1* h,
				Color_t  color,
				Style_t  marker,
				Double_t size,
				Style_t  fill,
				Style_t  line,
				Width_t  width)
{
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(size);
  h->SetFillColor(color);
  h->SetFillStyle(fill);
  h->SetLineColor(color);
  h->SetLineStyle(line);
  h->SetLineWidth(width);
  h->GetXaxis()->SetNdivisions(210);
  h->GetYaxis()->SetNdivisions(210);
}

//____________________________________________________________________
Double_t AliTrackletdNdeta2::MeanY(TH1* h, Double_t& e)
{
  Double_t sum = 0;
  Double_t sumw = 0;
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    Double_t c =  h->GetBinContent(i);
    Double_t e =  h->GetBinError  (i);
    if (c < 1e-6 || e/c > 0.1) continue;
    Double_t w =  1/e/e;
    sum        += w * c;
    sumw       += w;
  }
  e = TMath::Sqrt(1/sumw);
  return sum / sumw;
}
//____________________________________________________________________
Double_t AliTrackletdNdeta2::MeanZ(TH2* h, Double_t& e)
{
  Double_t sum = 0;
  Double_t sumw = 0;
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    for (Int_t j = 1; j <= h->GetNbinsY(); j++) { 
      Double_t c =  h->GetBinContent(i,j);
      Double_t e =  h->GetBinError  (i,j);
      if (c < 1e-6 || e/c > 0.1) continue;
      Double_t w =  1/e/e;
      sum        += w * c;
      sumw       += w;
    }
  }
  e = TMath::Sqrt(1/sumw);
  return sum / sumw;
}
  



#endif
//____________________________________________________________________
//
// EOF
//
