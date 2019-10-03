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
#include <TFractionFitter.h>
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
    kdNdetas      = 0x0008,
    /** Draw alphas */
    kSpecies      = 0x0010,
    /** Draw delta information */   
    kDeltas       = 0x0020,
    /** Draw backgrounds */
    kDetails      = 0x0040,
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
    kDefaultViz   = 0x03ff
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
    /** Correct for decay of strange to secondary  */
    kDTFS = 0x0010,
    /** Correct for fakes even if doing DTFS */
    kFake = 0x0020,
    /** Correct for centrality assumetry */
    kCentAverage = 0x0040,
    /** MC closure test */
    kClosure      = 0x01000,
    /** Debug flat */
    kDebug = 0x100,
    /** Default processing options */
    kDefaultProc  = 0x0002 
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
  /** Upper cut on @f$\Delta@f$ tail read from real data input file */
  Double_t fTailMax;
  /** Highest value of @f$\Delta@f$ read from real data input file */
  Double_t fMaxDelta;
  /** Least background scalar */
  Double_t fMinK;
  /** Largest background scalar */
  Double_t fMaxK;
  /** Least primary correction */
  Double_t fMinAlpha;
  /** Largest primary correction */
  Double_t fMaxAlpha;
  /** Fudge factor */
  Double_t fFudge;
  /** Strangeness enhancement factor - scale sim to real */
  Double_t fSEF;
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
  Bool_t fRealIsSim;
  
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
   * @param fudge    Fudge parameter 
   */
  void Run(UInt_t      proc     = kDefaultProc,
	   UInt_t      viz      = kDefaultViz,
	   UShort_t    maxBins  = 9,
	   const char* dataName = "data.root",
	   const char* simName  = "sim.root",
	   const char* output   = 0,
	   Double_t    fudge    = 1);
  //__________________________________________________________________
  /** 
   * @{ 
   * @name Steer processing 
   */
  /** 
   * Process the data 
   * 
   * @param realTop Top-level container of real data 
   * @param simTop  Top-level container of simulated data 
   * @param outTop  Top-level output directory 
   * @param maxBins Maximum number of bins to process 
   * 
   * @return true on success
   */
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
   * @param realEff    Trigger efficiency for real data
   * @param simEff     Trigger efficiency for simulated data
   * @param outTop     Top-level output directory 
   * 
   * @return true on success 
   */
  Bool_t ProcessBin(Double_t    c1,
		    Double_t    c2,
		    Container*  realTop,
		    Container*  simTop,
		    Double_t    realEff,
		    Double_t    simEff,
		    TDirectory* outTop);
  /** 
   * Process a single centrality bin 
   * 
   * @param c1         Lower centrality bound
   * @param c2         Upper centrality bound 
   * @param realCont   Real container 
   * @param simCont    Simulation container 
   * @param realEff    Trigger efficiency for real data
   * @param simEff     Trigger efficiency for simulated data
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
		    Double_t    realEff,
		    Double_t    simEff,
		    TDirectory* outTop, 
		    TDirectory* outDir,
		    Int_t       dimen);
  /* @} */
  //__________________________________________________________________
  /** 
   * @{ 
   * @name Strangeness 
   */
  /** 
   * Calculate the strangeness enhancement factor by looking at the
   * species distributions relative to pions.
   * 
   * @param simCont  Simulation container 
   * @param centAxis Centrality axis 
   * @param out      Output directory 
   * 
   * @return true on success 
   */
  Bool_t CalculateSEF(Container*   simCont,
		      const TAxis* centAxis,
		      TDirectory*  out);
  /* @} */
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
   * @param simDeltaD   Simulated distinct combinatorial @f$ dN/d\Delta@f$  
   * @param fit         Template fit 
   */
  void WriteDeltas(TDirectory* outDir,
		   TH1* realDeltaM, TH1* realDeltaI,
		   TH1* simDeltaM,  TH1* simDeltaI,
		   TH1* simDeltaC,  TH1* simDeltaP,
		   TH1* simDeltaS,  TH1* simDeltaD,
		   TH1* fit);
  /** 
   * Do a template fit of 
   *
   *@f[ 
   \frac{dN}{d\Delta_M} = a\frac{dN}{d\Delta_C'} +
   b\frac{dN}{d\Delta_P'} + c\frac{dN}{d\Delta_S'} @f]
   *
   * That is, find parameters @f$ a, b, c@f$ that scales the
   * simulation combinatorics @f$ C'@f$, primaries @f$ P'@f$, and
   * secondaries @f$ S'@f$ to match the data @f$ M@f$
   * 
   * @param outDir 
   * @param realDeltaM 
   * @param simDeltaC 
   * @param simDeltaP 
   * @param simDeltaS 
   *
   * @return fitted sum 
   */
  TH1* FractionFit(TDirectory* outDir,
		   TH1* realDeltaM,
		   TH1* simDeltaC,
		   TH1* simDeltaP,
		   TH1* simDeltaS);
  /* @} */

  /** 
   * @{ 
   * @name Result extraction 
   */
  /** 
   * Calculate the result as 
   * @f[
   R_{\eta,\mathrm{IP}_z} = 
   \frac{G\prime_{\eta,\mathrm{IP}_z}}{
     (1-\beta\prime)M\prime_{\eta,\mathrm{IP}_z}}
     (1-\beta)M_{\eta,\mathrm{IP}_z}
   @f] 
   * where 
   *
   * - @f$ G'@f$ is the generated, primary, charged particles 
   * - @f$ M'@f$ is the measured tracklets in simulated data 
   * - @f$ \beta'=C'/M'@f$ is the fraction of background tracklets in
   *   simulated data
   * - @f$ C'@f$ is the tracklets from combinatorics in simulated data
   * - @f$ M@f$  is the measured tracklets in real data 
   * - @f$ \beta=k\beta'@f$ is the fraction of background tracklets in
   *   real data
   * - @f$ k@f$ is the scalar to real, background, tracklets from
   *   simulated, background, tracklets
   * 
   * @param realCont   Container of real data 
   * @param simCont    Container of simulated data 
   * @param realEff    Trigger efficiency for real data
   * @param simEff     Trigger efficiency for simulated data
   * @param outParent  Output directory 
   * @param deltaDimen Dimensionality of @f$\Delta@f$ distribution. 
   * 
   * @return The result histogram 
   */
  TH1* Results(Container*  realCont,
	       Container*  simCont,
	       Double_t    realEff,
	       Double_t    simEff,
	       TDirectory* outParent,
	       Int_t       deltaDimen);
  /** 
   * Flatten the result with respect to the mean centrality of the
   * bin.  This scales each IPz row by the mean centrality in that row
   * divided by the mean centrality of the full bin.  This only
   * happens if the appropriate flag is set.
   * 
   * @param h Result histogram 
   * @param c Real data container.
   */
  void FlattenCentrality(TH2* h, Container* c);
  
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
  /** 
   * Visualize the results 
   * 
   * @param realSums Top-level container of real data sums 
   * @param simSums  Top-level container of simulated data sums 
   * @param realRess Top-level container of real data results
   * @param simRess  Top-level container of simulated data results 
   * @param outTop   Top-level directory of output
   * @param maxBins  Maximum number of bins to visualize 
   * 
   * @return true on success 
   */
  Bool_t Visualize(Container*  realSums,
		   Container*  simSums,
		   Container*  realRess,
		   Container*  simRess,
		   TDirectory* outTop,
		   Int_t       maxBins);
  /** 
   * Visualize the general stuff 
   * 
   * @param realList Real data output list 
   * @param simList  Simulated data output list
   */
  void   VisualizeGeneral(Container* realList, Container* simList);
  /** 
   * Draw the used simulation weights
   * 
   * @param simList Simulation list 
   */
  void VisualizeWeights(Container* simList);
  /** 
   * Visualize the final output 
   * 
   * @param outDir Output directory 
   * @param i      Bin number 
   */
  void VisualizeFinal(TDirectory* outDir, Int_t i);
  /** 
   * Visualize the closure test
   * 
   * @param outDir Output directory 
   * @param i      Bin number 
   */
  void VisualizeClosure(TDirectory* outDir, Int_t i);
  /** 
   * Visualize a centrality bin 
   * 
   * @param c1       Lowest centrality limit 
   * @param c2       Highest centrality limit 
   * @param simList  Simulated data output 
   * @param outTop   Output top directory 
   * 
   * @return true on success 
   */
  Bool_t VisualizeBin(Double_t    c1,
		      Double_t    c2,
		      Container*  simList, 
		      TDirectory* outTop);
  /** 
   * Visualize species from simulation
   * 
   * @param simCont Centrality container of simulated data
   * 
   * @return true on success 
   */
  Bool_t VisualizeSpecies(Container* simCont);
  /** 
   * Visualize specie deltas from simulation
   * 
   * @param simCont Centrality container of simulated data
   * @param outDir  Output directory
   *
   * @return true on success 
   */
  Bool_t VisualizeSpeciesDelta(Container* simCont, TDirectory* outDir);
  /** 
   * Visualize primary particles 
   * 
   * @param simCont Centrality container of simulated data
   * 
   * @return true on success 
   */
  Bool_t VisualizePrimary(Container* simCont);
  /** 
   * Visualize the @f$\Delta@f$ distributions 
   * 
   * @param outTop Output top directory 
   * @param dimen  The dimension to show 
   * 
   * @return true on success 
   */
  Bool_t VisualizeDelta(TDirectory* outTop, Int_t dimen);
  /** 
   * Visualize the details of the calculations 
   * 
   * @param outTop Output directory 
   * @param dimen  The dimension to show 
   * 
   * @return true on success
   */
  Bool_t VisualizeDetails(TDirectory* outTop, Int_t dimen);
  /** 
   * Visualize the @f$ dN_{\mathrm{ch}}/d\eta@f$ and components
   * 
   * @param outTop Output top directory 
   * @param dimen  The dimension to show 
   * 
   * @return true on success 
   */
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
  static void StackMinMax(THStack* stack, Double_t& min, Double_t& max);
  /** 
   * The observable title 
   * 
   * @return Observable title 
   */
  const char* ObsTitle() const { return "d#it{N}_{ch}/d#eta"; }
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
namespace {
  /** 
   * A guard to suppress messages 
   */
  struct SuppressGuard2
  {
    /** The previous message level */
    Int_t save = 0;
    /** 
     * Constructor 
     * 
     * @param lvl Level to suppress to 
     */
    SuppressGuard2(Int_t lvl=2000)
    {
      save = gErrorIgnoreLevel;
      gErrorIgnoreLevel = lvl;
    }
    /** 
     * Destructor 
     */
    ~SuppressGuard2()
    {
      gErrorIgnoreLevel = save;
    }
  };
}

//====================================================================
AliTrackletdNdeta2::AliTrackletdNdeta2()
  : AliTrackletAODUtils(),
    fProc(0),
    fViz(0),
    fDeltaCut(0),
    fTailDelta(0),
    fTailMax(0),
    fMaxDelta(0),
    fMinK(.7),
    fMaxK(2),
    fMinAlpha(0),
    fMaxAlpha(2.5),
    fFudge(0),
    fSEF(1),
    fCanvas(0),
    fTop(0),
    fBody(0),
    fRealIsSim(false)
{}

//====================================================================
void AliTrackletdNdeta2::Run(UInt_t      proc,
			     UInt_t      viz, 
			     UShort_t    maxBins,
			     const char* dataName,
			     const char* simName,
			     const char* outName,
			     Double_t    fudge)
{
  // Store options 
  fProc  = proc;
  fViz   = viz;
  fFudge = fudge;

  Printf("***************************************************\n"
	 "  Data file:       %s\n"
	 "  Simulation file: %s\n"
	 "  Output file:     %s\n"
	 "  Processing:      0x%x\n"
	 "  Visualize:       0x%x\n"
	 "  Max(Nbins):      %d\n"
	 "  Fudge:           %f\n",	
	 dataName, simName, outName, proc, viz, maxBins, fudge);
  

  // Open the input files
  TFile* dataFile = 0;
  TFile* simFile  = 0;
  if (!(dataFile = OpenFile(dataName))) return;
  if (!(simFile  = OpenFile(simName)))  return;

  // Get some top-level contianers 
  fRealIsSim = false;
  const char* base     = "MidRapidity";
  Container*  realSums = GetC(dataFile, Form("%sSums",      base));
  Container*  realRess = GetC(dataFile, Form("%sResults",   base));
  Container*  simSums  = GetC(simFile,  Form("%sMCSums",    base));
  Container*  simRess  = GetC(simFile,  Form("%sMCResults", base));
  if (!realSums || !realRess) {
    realSums = GetC(dataFile, Form("%sMCSums",      base));
    realRess = GetC(dataFile, Form("%sMCResults",   base));
    if (realSums && realRess)
      Warning("Run","\n"
	      "*********************************************\n"
	      "* Testing MC vs. MC:                        *\n"
	      "*  'Data' file:      %23s *\n"
	      "*  Simulation file:  %23s *\n"
	      "*********************************************\n",
	      dataName, simName);
    fRealIsSim = true;
  }
  if (!realSums || !realRess || !simSums || !simRess) return;

  // Get parameters from the real data file 
  Container* params   = GetC(realSums, "parameters");
  fDeltaCut  = GetD(params, "DeltaCut");
  fTailDelta = GetD(params, "TailDelta");
  fTailMax   = GetD(params, "TailMax");
  fMaxDelta  = GetD(params, "MaxDelta");

  // Create output file name 
  TString outBase(outName);
  if (outBase.IsNull())          outBase.Form("MiddNdeta_0x%04x", fProc);
  if (outBase.EndsWith(".root")) outBase.ReplaceAll(".root", "");

  // Make output directory
  gSystem->mkdir(outBase);
  
  // Open the output file 
  TFile* out = TFile::Open(Form("%s/result.root", outBase.Data()), "RECREATE");
  
  
  if (!Process(realRess, simRess, out, maxBins)) return;
  
  out->Write();

  if (fViz == 0) return; // Do not visualize anything
  
  Visualize(realSums, simSums, realRess, simRess, out, maxBins);
}

//====================================================================
Bool_t AliTrackletdNdeta2::Process(Container*  realTop,
				   Container*  simTop,
				   TDirectory* outDir,
				   Int_t       maxBins)
{
  DebugGuard g(fProc&kDebug,1,"Processing");
  TH1* realCent = CopyH1(realTop, "cent", "realCent");
  TH1* simCent  = CopyH1(simTop,  "cent", "simCent");
  TH1* realIPz  = GetH1(realTop,  "ipz");
  TH1* simIPz   = GetH1(simTop,   "ipz");

  // Check consistency of found histograms 
  if (!CheckConsistency(realCent, simCent)) {
    Warning("Process", "Centrality bins are incompatible, giving up");
    return false;
  }
  if (!CheckConsistency(realIPz, simIPz)) {
    Warning("Process", "IPz bins are incompatible, giving up");
    return false;
  }
  
  // Check if we're doing a closure test 
  if (fProc & kClosure) realTop = simTop;

  PrintAxis(*realCent->GetXaxis(), 2, "Real centrality");
  PrintAxis(*simCent ->GetXaxis(), 2, "Simulated centrality");

  THStack* mids = 0;
  TH1*     publ = 0;
  if (realCent->GetNbinsX() > 1) {
    Printf("Creating stack for mid rapidity results");
    mids            = new THStack("mids", "");
    Int_t    nbin   = 9;
    Double_t bins[] = { 0., 5., 10., 20., 30., 40., 50., 60., 70., 80.,  };
    Double_t vals[] = { 1948, 1590, 1180, 786, 512, 318, 183, 96.3, 44.9, };
    Double_t errs[] = {   38,   32,   31,  20,  15,  12,   8,  5.8,  3.4, };
    publ            = new TH1D("published", "PRL116,222302 (2016)",
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
    Printf("Making directory for final result");
    TDirectory* dd  = outDir->mkdir(Form("final%dd", d));
    if (realCent->GetNbinsX() > 1) {
      TH1*        mid = Make1D(0,"mid",Form("%s|_{|#eta|<0.5}", ObsTitle()),
			       2+d, 20, *(realCent->GetXaxis()));    
      mid->SetDirectory(dd);
      mid->SetXTitle("Centrality [%]");
      mid->SetYTitle(mid->GetTitle());
      mids->Add(mid);
    }
    THStack* full = new THStack("full","");
    dd->cd();
    full->Write();

    if (fRealIsSim) {
      THStack* clos = new THStack("closure","");
      clos->Write();
    }
    outDir->cd();
  }

  // Write centralities to file 
  realCent->Write();
  simCent ->Write();

  // Calculate strangeness enhancement factor
  CalculateSEF(simTop,simCent->GetXaxis(),outDir);

  // Get trigger efficiencies
  Double_t realEff  = GetD(realTop, "triggerEfficiency", 1);
  Double_t simEff   = GetD(simTop,  "triggerEfficiency", 1);
  Printf("Trigger efficiencies:  real=%6.4f  sim=%6.4f", realEff, simEff);

  // Process MB bin 
  if (realCent->GetNbinsX() <= 1)
    ProcessBin(0, 0, realTop, simTop, realEff, simEff, outDir);
  
  // Loop over defined centrality bins 
  for (Int_t i = 1; i <= realCent->GetNbinsX() && i <= maxBins ; i++) {
    Double_t c1 = realCent->GetXaxis()->GetBinLowEdge(i);
    Double_t c2 = realCent->GetXaxis()->GetBinUpEdge (i);
      
    ProcessBin(c1, c2, realTop, simTop, 1, 1, outDir);
  }
  // Close output file
  outDir->cd();
  if (mids) {
    mids->Write();
    
    TH1* mid = static_cast<TH1*>(mids->GetHists()->At(1));
    for (Int_t i = 1; i <= mid->GetNbinsX(); i++) {
      Double_t c1 = mid->GetXaxis()->GetBinLowEdge(i);
      Double_t c2 = mid->GetXaxis()->GetBinUpEdge(i);
      Int_t    j  = publ->GetXaxis()->FindBin((c1+c2)/2);
      if (j < 1 || j > publ->GetNbinsX()) continue;
      Double_t vh = mid ->GetBinContent(i);
      Double_t eh = mid ->GetBinError  (i);
      Double_t vp = publ->GetBinContent(j);
      Double_t ep = publ->GetBinError  (j);
      Double_t er, r = RatioE(vh,eh, vp, ep, er);
      Printf("%5.1f - %5.1f%%: "
	     "Here %6.1f +/- %4.1f "
	     "Published %6.1f +/- %4.1f "
	     "Ratio %5.3f +/- %5.3f",
	     c1, c2, vh, eh, vp, ep, r, er);
    }
  }
  
  return true;
}


//____________________________________________________________________
Bool_t AliTrackletdNdeta2::ProcessBin(Double_t    c1,
				      Double_t    c2,
				      Container*  realTop,
				      Container*  simTop,
				      Double_t    realEff,
				      Double_t    simEff,
				      TDirectory* outTop)
{
  // Form the folder name
  TString centName(CentName(c1,c2));
  if (TMath::Abs(c1 - c2) < 1e-6) centName = "all";
  DebugGuard g(fProc&kDebug,1,"Processing bin %s", centName.Data());
  
  // Get centrality containers 
  Container* realCont = GetC(realTop, centName);
  Container* simCont  = GetC(simTop,  centName);
  if (!realCont || !simCont) return false;

  // Test if we have any data at all
  TH1*       realDist  = GetH1(GetC(realCont, "measured"), "delta");
  TH1*       simDist   = GetH1(GetC(realCont, "measured"), "delta");
  if (realDist->GetEntries() <= 0 ||
      simDist ->GetEntries() <= 0) {
    Warning("ProcessBin", "No entries for bin %s", centName.Data());
    return false;
  }
  TDirectory* outDir = outTop->mkdir(centName);

  printf("%5.1f - %5.1f%%:", c1, c2);
  for (Int_t i = 0; i < 4; i++) {
    if ((fProc & (1 << i)) == 0) continue;
    if (!ProcessBin(c1, c2, realCont, simCont, realEff, simEff,
		    outTop, outDir, i))
      return false;
  }
  printf("\n");
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletdNdeta2::ProcessBin(Double_t    c1,
				      Double_t    c2,
				      Container*  realCont,
				      Container*  simCont,
				      Double_t    realEff,
				      Double_t    simEff,
				      TDirectory* outTop, 
				      TDirectory* outDir,
				      Int_t       dimen)
{
  if (!outDir) {
    Warning("ProcessBin", "No directory passed for %s and %s",
	    realCont->GetName(), simCont->GetName());
    return false;
  }  
  DebugGuard g(fProc&kDebug,1,"Processing bin %s", realCont->GetName());
  if (!Deltas(realCont, simCont, outDir, dimen)) return false;
  TH1* dndeta = Results(realCont, simCont, realEff, simEff, outDir, dimen);
  if (!dndeta) {
    Warning("ProcessBin", "Failed on Deltas for %f - %f", c1, c2);
    return false;
  }

  TDirectory* final = GetT(outTop,Form("final%dd", dimen));
  if (!final) {
    Warning("ProcessBin", "Failed on results for %f - %f", c1, c2);
    return false;
  }

  TH1*     mid  = GetH1(final, "mid");
  THStack* full = GetHS(final, "full");
  if (!final) {
    Warning("ProcessBin", "Missing full (%p)", full);
    return false;
  }

  Int_t b = 10;
  if (mid) {
    TF1* f = static_cast<TF1*>(dndeta->GetListOfFunctions()->At(0));
    if (!f) {
      Warning("ProcessBin", "No fit found on %s", dndeta->GetTitle());
      return false;
    }
    
    Double_t c = (c1+c2)/2;
    b          = mid->GetXaxis()->FindBin(c);
    if (b < 1 || b > mid->GetNbinsX()) {
      Warning("ProcessBin", "Centrality %f - %f out of range", c1, c2);
      return false;
    }   
    mid->SetBinContent(b, f->GetParameter(0));
    mid->SetBinError  (b, f->GetParError (0));
  }
  
  Color_t tc   = CentColor(b);
  TH1*    copy = static_cast<TH1*>(dndeta->Clone(outDir->GetName()));
  copy->SetDirectory(final);
  copy->GetListOfFunctions()->Clear();
  copy->SetTitle(Form("%5.1f#minus%5.1f%%", c1, c2));
  SetAttr(copy, tc);
  full->Add(copy);
  final->cd();
  full->Write(full->GetName(), TObject::kOverwrite);

  THStack* clos = GetHS(final, "closure", false);
  TH1*     clss = GetH1(GetT(outDir, Form("results%dd",dimen)),
			"closure", false);
  if (clos && clss) {
    copy = static_cast<TH1*>(clss->Clone(outDir->GetName()));
    copy->SetDirectory(0);
    copy->GetListOfFunctions()->Clear();
    copy->SetTitle(Form("%5.1f#minus%5.1f%%", c1, c2));
    SetAttr(copy, tc);
    clos->Add(copy);
    clos->Write(clos->GetName(), TObject::kOverwrite);
  }
    
  
  return true;
}

//====================================================================
Bool_t AliTrackletdNdeta2::CalculateSEF(Container*   simCont,
					const TAxis* centAxis,
					TDirectory*  out)
{
  DebugGuard g(fProc&kDebug,1,"Calculating strangeness enhancement factor");
  Double_t c1 = centAxis->GetBinLowEdge(1);
  Double_t c2 = centAxis->GetBinUpEdge(1);
  TString  centName(CentName(c1,c2));
  if (TMath::Abs(c1 - c2) < 1e-6) centName = "all";

  // Find the most central bin 
  Container* centBin = GetC(simCont, centName);
  if (!centBin) return false;

  Container* generated = GetC(centBin, "generated");
  if (!generated) return false;

  Container* mix   = GetC(generated, "mix");
  if (!mix) return false;

  THStack* toPion = GetHS(mix, "toPion");
  if (!toPion) return false;

  TIter    next(toPion->GetHists());
  TH1*     hist = 0;
  Double_t sum  = 0;
  Double_t sumw = 0;
  while ((hist = static_cast<TH1*>(next()))) {
    Double_t r2760 = hist->GetBinContent(0);
    Double_t e2760 = 0.07*r2760; // Fixed 7% error 

    TF1* f = new TF1("f", "pol0", -.5, +.5);
    hist->Fit(f, "QN", "", -.5, +.5);
    Double_t rHere = f->GetParameter(0);
    Double_t eHere = f->GetParError (0);

    Double_t eCh, rCh = RatioE(r2760, e2760, rHere, eHere, eCh);
#if 0    
    Printf("%20s: @ 2.76TeV=%6.4f+/-%6.4f  here=%6.4f+/-%6.4f -> %6.4f+/-%6.4f",
	   hist->GetTitle(), r2760, e2760, rHere, eHere, rCh, eCh);
#endif 
    delete f;

    sum  += rHere*rCh;
    sumw += rHere;
  }  
  Double_t avg = 1;
  if (sumw>0) avg = sum / sumw;
  Printf(" Weighted average of factor: %6.4f", avg);
  Printf(" Preset:                     %6.4f", fSEF);
  if (fSEF == 1) fSEF = avg;
  Printf( "Strangeness enhancement factor set to %6.4f", fSEF);
  
  return true;
}

//====================================================================
Bool_t AliTrackletdNdeta2::Deltas(Container*  realCont,
				  Container*  simCont,
				  TDirectory* outParent,
				  Int_t       dim)
{
  DebugGuard g(fProc&kDebug,1,"Doing Delta calculations");
  if (!outParent) {
    Warning("Deltas", "No directory passed!");
    return false;
  }
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
  DebugGuard g(fProc&kDebug,1,"Doing 0D Delta calculations");
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
  TH1* simDeltaD  = CopyH1(GetC(simCont, "distinct"),     "delta","simDeltaD");

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
  
  TH1* fit = FractionFit(outDir, realDeltaM, simDeltaC, simDeltaP, simDeltaS);
  WriteDeltas(outDir, realDeltaM, realDeltaI, simDeltaM, simDeltaI,
	      simDeltaC, simDeltaP, simDeltaS, simDeltaD, fit);

  outParent->cd();
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletdNdeta2::Deltas1D(Container*  realCont,
				    Container*  simCont,
				    TDirectory* outParent)
{
  DebugGuard g(fProc&kDebug,1,"Doing 1D Delta calculations");
  // Make an output directory
  if (!outParent) {
    Warning("Deltas1D", "No directory passed!");
    return false;
  }
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
  TH1* simDeltaD  = CopyH1(GetC(simCont, "distinct"),     "delta","simDeltaD");

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
  
  TH1* toScale[]  = { simDeltaM,simDeltaI,simDeltaC,
		      simDeltaP,simDeltaS,simDeltaD, 0};
  TH1**    pScale     = toScale;
  while ((*pScale)) { 
    Scale(*pScale, scale, scaleE);
    (*pScale)->SetTitle(Form("k_{M}#times%s", (*pScale)->GetTitle()));
    pScale++;
  }

  TH1* fit = FractionFit(outDir, realDeltaM, simDeltaC, simDeltaP, simDeltaS);
  WriteDeltas(outDir, realDeltaM, realDeltaI, simDeltaM, simDeltaI,
	      simDeltaC, simDeltaP, simDeltaS, simDeltaD, fit);

  outParent->cd();
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletdNdeta2::Deltas2D(Container*  realCont,
				    Container*  simCont,
				    TDirectory* outParent)
{
  DebugGuard g(fProc&kDebug,1,"Doing 2D Delta calculations");
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
  TGraphErrors* avgScale = new TGraphErrors(2);
  avgScale->SetLineStyle(2);
  avgScale->SetLineColor(kBlack);
  avgScale->SetFillColor(kYellow);
  avgScale->SetFillStyle(3002);
  avgScale->SetPoint     (0, scale->GetXaxis()->GetXmin(), s); 
  avgScale->SetPoint     (1, scale->GetXaxis()->GetXmax(), s);
  avgScale->SetPointError(0, 0, sE);
  avgScale->SetPointError(1, 0, sE);
  scale->GetListOfFunctions()->Add(avgScale, "le3");
  
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
  TH2* s2DeltaD = CopyH2(GetC(simCont, "distinct"),     "etaDelta","s2DeltaD");

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
  
  TH2*  toScale[] = { s2DeltaM,s2DeltaI,s2DeltaC,s2DeltaP,s2DeltaS,s2DeltaD,0 };
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
  TH1* sDeltaD = ProjectDelta(s2DeltaD); 
  rDeltaM->SetTitle(r2DeltaM->GetTitle()); rDeltaM->SetName("realDeltaM");
  rDeltaI->SetTitle(r2DeltaI->GetTitle()); rDeltaI->SetName("realDeltaI");
  sDeltaM->SetTitle(s2DeltaM->GetTitle()); sDeltaM->SetName("simDeltaM");
  sDeltaI->SetTitle(s2DeltaI->GetTitle()); sDeltaI->SetName("simDeltaI");
  sDeltaC->SetTitle(s2DeltaC->GetTitle()); sDeltaC->SetName("simDeltaC");
  sDeltaP->SetTitle(s2DeltaP->GetTitle()); sDeltaP->SetName("simDeltaP");  
  sDeltaS->SetTitle(s2DeltaS->GetTitle()); sDeltaS->SetName("simDeltaS");
  sDeltaD->SetTitle(s2DeltaD->GetTitle()); sDeltaD->SetName("simDeltaD");
  
  TH1* f2  = FractionFit(outDir, r2DeltaM, s2DeltaC, s2DeltaP, s2DeltaS);
  TH1* fit = ProjectDelta(static_cast<TH2*>(f2));
  WriteDeltas(outDir,rDeltaM,rDeltaI,
	      sDeltaM,sDeltaI,sDeltaC,
	      sDeltaP,sDeltaS,sDeltaD,
	      fit);

  TDirectory* full = outDir->mkdir("full");
  r2DeltaM->SetDirectory(full); r2DeltaM->SetName("realDeltaM");
  r2DeltaI->SetDirectory(full); r2DeltaI->SetName("realDeltaI");
  s2DeltaM->SetDirectory(full); s2DeltaM->SetName("simDeltaM");
  s2DeltaI->SetDirectory(full); s2DeltaI->SetName("simDeltaI");
  s2DeltaC->SetDirectory(full); s2DeltaC->SetName("simDeltaC");
  s2DeltaP->SetDirectory(full); s2DeltaP->SetName("simDeltaP");
  s2DeltaS->SetDirectory(full); s2DeltaS->SetName("simDeltaS");
  s2DeltaD->SetDirectory(full); s2DeltaD->SetName("simDeltaD");

  outParent->cd();
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletdNdeta2::Deltas3D(Container*  realCont,
				    Container*  simCont,
				    TDirectory* outParent)
{
  DebugGuard g(fProc&kDebug,1,"Doing 3D Delta calculations");
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
  TGraphErrors* avgScale = new TGraphErrors(2);
  avgScale->SetLineStyle(2);
  avgScale->SetLineColor(kBlack);
  avgScale->SetFillColor(kYellow);
  avgScale->SetFillStyle(3002);
  avgScale->SetPoint(0, etaScale->GetXaxis()->GetXmin(), s); 
  avgScale->SetPoint(1, etaScale->GetXaxis()->GetXmax(), s);
  avgScale->SetPointError(0, 0, sE);
  avgScale->SetPointError(1, 0, sE);
  etaScale->GetListOfFunctions()->Add(avgScale, "le3");
  
  // Get the raw projected Delta distributions for each component
  const char* nm = "etaDeltaIPz";
  TH3* r3DeltaM = CopyH3(realMeas,                      nm,"r3DeltaM");
  TH3* r3DeltaI = CopyH3(GetC(realCont,"injected"),     nm,"r3DeltaI");
  TH3* s3DeltaM = CopyH3(simMeas,                       nm,"s3DeltaM");
  TH3* s3DeltaI = CopyH3(GetC(simCont, "injected"),     nm,"s3DeltaI");
  TH3* s3DeltaC = CopyH3(GetC(simCont, "combinatorics"),nm,"s3DeltaC");
  TH3* s3DeltaP = CopyH3(GetC(simCont, "primaries"),    nm,"s3DeltaP");
  TH3* s3DeltaS = CopyH3(GetC(simCont, "secondaries"),  nm,"s3DeltaS");
  TH3* s3DeltaD = CopyH3(GetC(simCont, "distinct"),     nm,"s3DeltaD");

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
#if 0
  TH2* scale2 = static_cast<TH2*>(scale->Clone("scaleMain"));
  scale2->SetDirectory(0);
  scale2->Reset();
  Int_t sigBin = r3DeltaM->GetYaxis()->FindBin(1.5);
  for (Int_t i = 1; i <= r3DeltaM->GetNbinsX(); i++) {
    for (Int_t j = 1; j <= r3DeltaM->GetNbinsZ(); j++) {
      // Integrate over Delta 
      Double_t rintg = 0, reintg = 0;
      Double_t sintg = 0, seintg = 0;
      rintg = r3DeltaM->IntegralAndError(i,i,1,sigBin,j,j,reintg);
      sintg = s3DeltaM->IntegralAndError(i,i,1,sigBin,j,j,seintg);
      Double_t re, r = RatioE(rintg, reintg, sintg, seintg, re);
      
      scale2->SetBinContent(i, j, r);
      scale2->SetBinError  (i, j, re);
    }
  }
  Double_t rS2, rS = MeanZ(scale2, rS2);  
  Printf("Scalar of Inj %6.3f +/- %6.3f", rS, rS2);
#endif 
  
  TH3*  toScale[] = { s3DeltaM,s3DeltaI,s3DeltaC,s3DeltaP,s3DeltaS,s3DeltaD,0};
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
  TH1* sDeltaD = ProjectDeltaFull(s3DeltaD); 
  rDeltaM->SetTitle(r3DeltaM->GetTitle()); rDeltaM->SetName("realDeltaM");
  rDeltaI->SetTitle(r3DeltaI->GetTitle()); rDeltaI->SetName("realDeltaI");
  sDeltaM->SetTitle(s3DeltaM->GetTitle()); sDeltaM->SetName("simDeltaM");
  sDeltaI->SetTitle(s3DeltaI->GetTitle()); sDeltaI->SetName("simDeltaI");
  sDeltaC->SetTitle(s3DeltaC->GetTitle()); sDeltaC->SetName("simDeltaC");
  sDeltaP->SetTitle(s3DeltaP->GetTitle()); sDeltaP->SetName("simDeltaP");  
  sDeltaS->SetTitle(s3DeltaS->GetTitle()); sDeltaS->SetName("simDeltaS");
  sDeltaD->SetTitle(s3DeltaD->GetTitle()); sDeltaD->SetName("simDeltaD");
    
  TH1* f3  = FractionFit(outDir, r3DeltaM, s3DeltaC, s3DeltaP, s3DeltaS);
  TH3* ff3 = static_cast<TH3*>(f3);

  TDirectory* full = outDir->mkdir("full");
  r3DeltaM->SetDirectory(full); r3DeltaM->SetName("realDeltaM");
  r3DeltaI->SetDirectory(full); r3DeltaI->SetName("realDeltaI");
  s3DeltaM->SetDirectory(full); s3DeltaM->SetName("simDeltaM");
  s3DeltaI->SetDirectory(full); s3DeltaI->SetName("simDeltaI");
  s3DeltaC->SetDirectory(full); s3DeltaC->SetName("simDeltaC");
  s3DeltaP->SetDirectory(full); s3DeltaP->SetName("simDeltaP");
  s3DeltaS->SetDirectory(full); s3DeltaS->SetName("simDeltaS");
  s3DeltaD->SetDirectory(full); s3DeltaD->SetName("simDeltaD");
  TH1* fit = 0;
  if (ff3){    
    ff3->SetDirectory(full);
    ff3->SetName("simDeltaF");
    TH2* fetaDelta = static_cast<TH2*>(ff3->Project3D("yx e"));
    fetaDelta->SetName("simDeltaFF");
    fetaDelta->SetDirectory(full);
    fetaDelta->Scale(1./ff3->GetNbinsZ());
    outDir->cd();
    fit = fetaDelta->ProjectionY("simDeltaF");
    fit->SetTitle("#Delta_{F}");
    fit->SetDirectory(outDir);
    fit->Scale(1. / fetaDelta->GetNbinsX());
    // delete fetaDelta;
  }
  WriteDeltas(outDir,rDeltaM,rDeltaI,
	      sDeltaM,sDeltaI,sDeltaC,
	      sDeltaP,sDeltaS,sDeltaD,
	      fit);
  
  outParent->cd();
  return true;
}

//____________________________________________________________________
TH1* AliTrackletdNdeta2::FractionFit(TDirectory* outDir,
				     TH1*        rDeltaM,
				     TH1*        sDeltaC,
				     TH1*        sDeltaP,
				     TH1*        sDeltaS)
{
  DebugGuard g(fProc&kDebug,1,"Fraction fitting in %s", outDir->GetName());
  // We don't do this, as it doesn't seem to do much. 
  return 0;
  if (!rDeltaM || !sDeltaC || !sDeltaP || !sDeltaS) {
    Warning("FractionFit", "Missing M=%p, C'=%s, P'=%s, or S'=%p",
	    rDeltaM, sDeltaC, sDeltaP, sDeltaS);
    return 0;
  }
  TDirectory* savDir = gDirectory;
  gROOT->cd();
  Double_t intg  = rDeltaM->Integral();
  Double_t mintg = sDeltaP->Integral()+sDeltaS->Integral()+sDeltaC->Integral();
  TH1*     dat   = static_cast<TH1*>(rDeltaM->Clone("tmpM"));
  dat->SetDirectory(0);
  dat->Scale(1./intg);
  TH1*  sig   = static_cast<TH1*>(sDeltaP->Clone("tmpPS"));
  sig->SetDirectory(0);
  sig->Add(sDeltaS);
  sig->Scale(1./sig->Integral()); // mintg);
  TH1*  bg    = static_cast<TH1*>(sDeltaC->Clone("tmpC"));
  bg->SetDirectory(0);
  bg->Scale(1./bg->Integral()); // mintg);

  TObjArray mc;
  mc.SetOwner();
  mc.Add(sig);
  mc.Add(bg);
  // mc.Add(arr[3]);
  TFractionFitter f(dat, &mc, "Q");
  Int_t status = f.Fit();
  savDir->cd();
  if (status != 0) {
    Warning("FractionFit", "Fit failed w/status=%d", status);
    return 0;
  }
  Printf("\nTemplate fits");
  for (Int_t i = 0; i < 2; i++) {
    Double_t v, e;
    f.GetResult(i, v, e);
    Printf("%30s=%6.4f +/- %6.4f",
	   mc.At(i)->GetName(), e, v);
  }
  TH1* ret = f.GetPlot();
  ret->Scale(intg);
  delete dat;
  return ret;
}

//____________________________________________________________________
void AliTrackletdNdeta2::WriteDeltas(TDirectory* outDir,
				     TH1* rDeltaM, TH1* rDeltaI,
				     TH1* sDeltaM, TH1* sDeltaI,
				     TH1* sDeltaC, TH1* sDeltaP,
				     TH1* sDeltaS, TH1* sDeltaD,
				     TH1* fit)
{
  DebugGuard g(fProc&kDebug,1,"Writing Deltas to %s", outDir->GetName());
  
  THStack* all    = new THStack("all", "");
  SetAttr(rDeltaM, kRed+2,    20, 1.0);
  SetAttr(rDeltaI, kOrange+2, 21, 1.0);
  rDeltaM->SetDirectory(outDir);
  rDeltaI->SetDirectory(outDir);
  all->Add(rDeltaM);
  all->Add(rDeltaI);
  
  TH1*     toScale[] = {sDeltaM,sDeltaI,sDeltaC, sDeltaP,sDeltaS,sDeltaD,fit,0};
  Color_t  toColor[] = {kRed,   kOrange,kMagenta,kGreen, kBlue,  kPink, kBlack};
  Style_t  toStyle[] = {24,     25,     30,      26,     32,     30, 24    };
  TH1**    pScale    = toScale;
  Color_t* pColor    = toColor;
  Style_t* pStyle    = toStyle;
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
  TH1*     ratioM = static_cast<TH1*>(sDeltaM->Clone("ratioM"));
  ratioM->SetTitle("#Delta_{M'}/#Delta_{M}");
  ratioM->Divide(rDeltaM);
  ratioM->SetDirectory(outDir);
  ratios->Add(ratioM);

  TH1*     ratioI = static_cast<TH1*>(sDeltaI->Clone("ratioI"));
  ratioI->SetTitle("#Delta_{I'}/#Delta_{I}");
  ratioI->Divide(rDeltaI);
  ratioI->SetDirectory(outDir);
  ratios->Add(ratioI);

  TH1*     ratioIC = static_cast<TH1*>(sDeltaC->Clone("ratioIC"));
  ratioIC->SetTitle("#Delta_{C'}/#Delta_{I}");
  ratioIC->Divide(rDeltaI);
  ratioIC->SetDirectory(outDir);
  ratios->Add(ratioIC);

  if (!fit) { ratios->Write(); return; }
  
  TH1*     ratioF = static_cast<TH1*>(fit->Clone("ratioF"));
  ratioF->SetTitle("#Delta_{fit}/#Delta_{M}");
  ratioF->Divide(rDeltaM);
  ratioF->SetDirectory(outDir);
  ratios->Add(ratioF);
  
  ratios->Write();
}


//====================================================================
TH1* AliTrackletdNdeta2::Results(Container*  realCont,
				 Container*  simCont,
				 Double_t    realEff,
				 Double_t    simEff,
				 TDirectory* outParent,
				 Int_t       deltaDimen)
{
  DebugGuard g(fProc&kDebug,1,"Calculating results for %s",realCont->GetName());
  
  TDirectory* outDir = outParent->mkdir(Form("results%dd", deltaDimen));
  TDirectory* delDir = outParent->GetDirectory(Form("delta%dd", deltaDimen));

  Bool_t fk = (fProc & kFake);
  Bool_t sc = (fProc & kDTFS);
  Bool_t rs = fRealIsSim;
  
  TH2* scale = static_cast<TH2*>(delDir->Get("scale"));
  TH2* realM =      CopyH2(GetC(realCont,"measured"),     "etaIPz", "realM");
  TH2* realS =      CopyH2(GetC(realCont,"measured"),     "etaIPz", "realS");
  TH2* realC = fk ? CopyH2(GetC(realCont,"measured"),     "etaIPz", "realC"):0;
  TH2* simM  =      CopyH2(GetC(simCont, "measured"),     "etaIPz", "simM");
  TH2* simC  = fk ? CopyH2(GetC(simCont, "combinatorics"),"etaIPz", "simC") :0;
  TH2* simG  =      CopyH2(GetC(simCont, "generated"),    "etaIPz", "simG");
  TH2* simS  =      CopyH2(GetC(simCont, "measured"),     "etaIPz", "simS");
  TH2* simA  =      CopyH2(GetC(simCont, "generated"),    "etaIPz", "simA");
  TH2* simB  = fk ? CopyH2(GetC(simCont, "combinatorics"),"etaIPz", "simB") :0;
  TH2* simT  = sc ? CopyH2(GetC(simCont, "secondaries"),  "dtfs",   "simT") :0;
  TH2* realG = rs ? CopyH2(GetC(realCont,"generated"),    "etaIPz", "realG"):0;
  TH1* realZ = CopyH1(realCont, "ipz", "realZ");
  TH1* simZ  = CopyH1(simCont,  "ipz", "simZ");

  Double_t realEt = GetD(realCont, "triggerEfficiency"); 
  Double_t simEt  = GetD(simCont,  "triggerEfficiency"); 
  Double_t realEz = GetD(realCont, "ipEfficiency"); 
  Double_t simEz  = GetD(simCont,  "ipEfficiency"); 
  
  if (scale->GetEntries() <= 0) {
    Warning("Results", "Scale has no entries");
  }
  
  // Scale combinatorial background to measured to get beta
  if (simB) {
    simB->Divide(simM);
    simB->SetTitle("#beta'");
    simB->SetZTitle("#beta'");
  }
  
  // Substract combinatorial background from measured to get signal
  if (simC) {
    simS->Add(simC, -1);
    simS->SetTitle("S'");
    simS->SetZTitle("S'");
  }
  
  // Possibly subtract secondaries from strange 
  if (simT) {
    simT->SetTitle("T'");
    simT->SetZTitle("T'");
    simS->Add(simT, -1);
  }
  
  // Scale MC truth primaries by signal to get correction 
  simA->Divide(simS);
  simA->SetTitle("A'");
  simA->SetZTitle("#alpha'");

  // Copy simulated beta to real beta, and scale by scalar
  TH2* realB = 0;
  if (simB) {
    realB = static_cast<TH2*>(simB->Clone("realB"));
    realB->SetDirectory(0);
    realB->SetTitle("#beta");
    realB->SetZTitle("#beta");
    realB->Multiply(scale);
    Scale(realB, fFudge, 0);
  }
  
  // Multiply real beta onto real measured to get background
  if (realC) { 
    realC->Multiply(realB);
    realC->SetTitle("C");
    realC->SetZTitle("C");
  }
  
  // Substract the real background off the real measured
  if (realC) {
    realS->Add(realC, -1);
    realS->SetTitle("S");
    realS->SetZTitle("S");
  }
  
  // Possibly subtract secondaries from strange
  TH2* realT = 0;
  if (simT) {
    realT = static_cast<TH2*>(simT->Clone("realT"));
    realT->SetDirectory(0);
    realT->SetTitle("T");
    realT->SetZTitle("T");
    realT->Scale(fSEF);
    realS->Add(realT, -1);
  }
  
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
  /*      */ realM->Multiply(fiducial);
  if (realC) realC->Multiply(fiducial);
  /*      */ realS->Multiply(fiducial);
  if (realB) realB->Multiply(fiducial);
  /*      */ simM ->Multiply(fiducial);
  if (simC)  simC ->Multiply(fiducial);
  /*      */ simS ->Multiply(fiducial);
  if (simB)  simB ->Multiply(fiducial);
  /*      */ simA ->Multiply(fiducial);
  if (simT)  simT ->Multiply(fiducial);
  if (realT) realT->Multiply(fiducial);
  
  // We can now make our result
  TH2* result = static_cast<TH2*>(realS->Clone("result"));
  result->Multiply(simA);
  result->SetTitle("R");
  result->SetZTitle("R");

  FlattenCentrality(result,realCont);
		    
  // Output directory for full stuff
  TDirectory* full = outDir->mkdir("full");
  
  // Make a stack of 1D projections
  struct Rec {
    TH2*     h;  // histogram to average 
    TH1*     s;  // IPz distribution
    Double_t e;  // IP efficiency 
    TH2*     c;  // Mask for bins to use 
    TH1*     p;  // The result 
    Style_t sty; // Marker style 
    Color_t col; // Color 
    Float_t siz; // Marker size 
    const char* tit; // Title 
  };
  Rec   sT = { simT,    simZ, simEz, 0,     0,30,kSpring+2, 1.4,"strangeness"};
  Rec   sC = { simC,    simZ, simEz, result,0,32,kMagenta+2,1.4,"background"};
  Rec   sS = { simS,    simZ, simEz, result,0,27,kGreen+2,  1.8,"signal" };
  Rec   sM = { simM,    simZ, simEz, result,0,26,kBlue+2,   1.4,"measured" };
  Rec   sG = { simG,    simZ, simEz, 0,     0,24,kRed+2,    1.4,"generated" };
  Rec   sB = { simB,    simZ, simEz, 0,     0,28,kPink+2,   1.4,"#beta" };
  Rec   rC = { realC,   realZ,realEz,result,0,23,kMagenta+2,1.2,"background"};
  Rec   rT = { realT,   realZ,realEz,0,     0,29,kSpring+2, 1.4,"strangeness"};
  Rec   rS = { realS,   realZ,realEz,result,0,33,kGreen+2,  1.6,"signal" };
  Rec   rM = { realM,   realZ,realEz,result,0,22,kBlue+2,   1.2,"measured" };
  Rec   rR = { result,  realZ,realEz,0,     0,20,kRed+2,    1.3,ObsTitle() };
  Rec   rB = { realB,   realZ,realEz,0,     0,34,kPink+2,   1.4,"#beta" };
  Rec   sA = { simA,    simZ, simEz, 0,     0,30,kSpring+2, 1.4,"#alpha" };
  Rec   sF = { fiducial,simZ, simEz, 0,     0,31,kSpring+2, 1.4,"fiducial" };
  Rec   rG = { realG,   realZ,realEz,0,     0,24,kBlack,    1.4,"truth" };
  Rec*  recs[]={ &rR, &sG, &rS, &sS, &rM, &sM, &rC, &sC,
		 &rB, &sB, &sA, &sF, &rG, &rT, &sT, 0 };
  Rec** ptr   = recs;

  TH1*     dndeta  = 0;
  THStack* all     = new THStack("all", "");
  if (fViz & kAltMarker) {
    rR.sty = 21;
    rR.siz = 1.2;
    sG.sty = 25;
    sG.siz = 1.3;
    rG.sty = 25;
    rG.siz = 1.3;
  }
  while ((*ptr)) {
    Rec*  src = *ptr;
    if (!src->h) { ptr++; continue; }
    src->h->SetDirectory(full);
    src->p = AverageOverIPz(src->h, src->h->GetName(), 1,
			    src->s, src->e, src->c);
    src->p->SetYTitle(src->h->GetZaxis()->GetTitle());
    src->p->SetTitle(Form("%s - %s", src->h->GetTitle(), src->tit));
    src->p->SetDirectory(outDir);
    if (src->e > 1e-3) src->p->Scale(src->e);
    if (src->h != simB && src->h != realB &&
	src->h != simA && src->h != fiducial) all->Add(src->p);
    SetAttr(src->p, src->col, src->sty, src->siz);
    if (src->h == result) {
      dndeta = src->p;
      dndeta->SetYTitle(ObsTitle());
    }
    ptr++;
  }
  // Show example calculation
  Int_t mi = result->GetXaxis()->FindBin(0.);
  Int_t mj = result->GetYaxis()->FindBin(0.);
  TString simST;
  TString realST;
  if (simT) {
    simST.Form("-%4.1f", simT->GetBinContent(mi,mj));
    realST.Form("-%4.1f", realT->GetBinContent(mi,mj));
  }
  printf("R=G'/[(1-beta')M'](1-beta)M="
	 "%6.1f /((1-%6.3f)*%6.1f%s)*(1-%6.3f)*(%6.1f%s)="
	 "%6.3f * %6.3f="
	 "%6.1f [%6.1f]",
	 /*  */ sG.h->GetBinContent(mi,mj),   
	 sB.h ? sB.h->GetBinContent(mi,mj) : 0,   
	 /*  */ sM.h->GetBinContent(mi,mj),   
	 /*  */ simST.Data(),
	 rB.h ? rB.h->GetBinContent(mi,mj) : 0,
	 /*  */ rM.h->GetBinContent(mi,mj),   
	 /*  */ realST.Data(),
	 /*  */ sA.h->GetBinContent(mi,mj),   
	 /*  */ rS.h->GetBinContent(mi,mj),   
	 /*  */ rR.h->GetBinContent(mi,mj),   
	 /*  */ rG.h ? rG.h->GetBinContent(mi,mj) : -1);

  if (rG.p) {
    TH1* ratio = RatioH(dndeta, rG.p, "closure");
    ratio->SetYTitle("Closure test");
    ratio->SetDirectory(outDir);
  }

  
  THStack* ratios  = new THStack("ratios", "");
  TH1*     scaleC  = 0;
  ratios->Add(RatioH(rM.p, sM.p, "rMeaured"));
  ratios->Add(RatioH(rC.p, sC.p, "rBackground"));
  ratios->Add(RatioH(rS.p, sS.p, "rSignal"));
  ratios->Add(RatioH(rR.p, sG.p, "rResult"));
  ratios->Add((scaleC = AverageOverIPz(scale,scale->GetName(),1,realZ,0,0)));
  scaleC->SetMarkerColor(kBlack);
  scaleC->SetMarkerColor(kGray);
  scaleC->SetMarkerStyle(31);
  scaleC->SetLineStyle(2);
  scaleC->SetLineColor(kGray);
  
  TF1* tmp = new TF1("mid", "pol0", -.5, +.5);
  dndeta->Fit(tmp, "Q0R+");
  TLatex* ltx = new TLatex(0,tmp->GetParameter(0)/2,
			   Form("%s|_{|#eta|<0.5}=%.1f#pm%.1f",
				ObsTitle(), 
				tmp->GetParameter(0),
				tmp->GetParError(0)));
#if 0
  printf(" %dD: %6.1f +/- %4.1f (%4.2f)",
	 deltaDimen, tmp->GetParameter(0), tmp->GetParError(0),
	 tmp->GetChisquare()/tmp->GetNDF());
#endif
  
  ltx->SetTextAlign(22);
  ltx->SetTextFont(42);
  dndeta->GetListOfFunctions()->Add(ltx);
    
  outDir->cd();
  all->Write();
  ratios->Write();
  if (realB) realB   ->SetDirectory(full);
  if (simB)  simB    ->SetDirectory(full);
  /*      */ simA    ->SetDirectory(full);
  /*      */ fiducial->SetDirectory(full);

  outParent->cd();
  dndeta->SetBinContent(0,                   realEff);
  dndeta->SetBinContent(dndeta->GetNbinsX()+1,simEff);
  
  return dndeta;
}

//____________________________________________________________________
void AliTrackletdNdeta2::FlattenCentrality(TH2* h, Container* c)
{
  if (!(fProc & kCentAverage)) return;
  DebugGuard g(fProc&kDebug,1,"Correcting for non-flat centrality in bin");

  TH1* centIpz = GetH1(c, "centIpz");
  if (!centIpz) {
    Warning("FlattenCentrality", "No average cent vs IPz data");
    return;
  }
  Double_t meanC = centIpz->GetMean(2);
  if (meanC < 0 || meanC > 100) {
    Warning("FlattenCentrality", "Average centrality %f funny", meanC);
    return;
  }
  for (Int_t iz = 1; iz <= centIpz->GetXaxis()->GetNbins(); iz++) {
    Double_t cc = centIpz->GetBinContent(iz);
    Double_t rc = meanC / cc;
#if 0
    DebugGuard g2(fProc&kDebug,1,"<c>=%7.4f%% <ci>=%7.4f%% -> %6.4f",
		  meanC, cc, rc);
#endif 
    for (Int_t ieta = 1; ieta <= h->GetXaxis()->GetNbins(); ieta++) {
      h->SetBinContent(ieta, iz, rc * h->GetBinContent(ieta, iz));
      h->SetBinError  (ieta, iz, rc * h->GetBinError  (ieta, iz)); 
    }
  }
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
    SuppressGuard2 g;    
    fCanvas->Print(Form("%s/summary.pdf[", outputName.Data()),
		   Form("pdf %s", (fViz & kLandscape) ? "Landscape" : ""));
  }
  // if (fViz & kPNG) {
  // }
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
    SuppressGuard2 g;
    fCanvas->Print(Form("%s/summary.pdf]", fCanvas->GetTitle()),
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
    SuppressGuard2 g;
    fCanvas->Print(Form("%s/summary.pdf",fCanvas->GetTitle()), tit);
  }
  static Int_t cnt = 1;
  if (fViz & kPNG) {
    SuppressGuard2 g;
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
  if (!p) {
    Warning("DrawInPad", "Sub-pad %d does not exist", pad);
    return 0;
  }
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
  if (!p || !l) return;
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
  TH1*     real  = CopyH1(realList, name, Form("real%s",name));
  TH1*     sim   = CopyH1(simList,  name, Form("sim%s",name));
  if (real->GetNbinsX() <= 1 && sim->GetNbinsX() <= 1) return 0;
  THStack* stack = new THStack(name, title);
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
  DebugGuard g(fProc&kDebug,1,"Visualizing");
  // --- Visualization -----------------------------------------------
  TH1* realCent = GetH1(outDir, "realCent");
  if (!realCent) {
    Warning("Visualize", "realCent histogram not found");
    return false;
  }
  TString outName(gSystem->DirName(outDir->GetName()));
  outName.ReplaceAll(".root", "");
  CreateCanvas(outName);
  if (fViz & kParameters) VisualizeParams(realSums, simSums);
  if (fViz & kGeneral)    VisualizeGeneral(realRess, simRess);
  if (fViz & kWeights)    VisualizeWeights(simRess);

  if (fViz &  kdNdetas) { 
    for (Int_t i = 0; i < 4; i++) {
      if ((fProc & (1 << i)) == 0) continue;
      VisualizeFinal(outDir, i);
      VisualizeClosure(outDir, i);
    }
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
  DebugGuard g(fProc&kDebug,1,"Visualizing general things");
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
    if (!s || !s->GetHists()) return;

    TIter next(s->GetHists());
    TH1*  h = 0;
    Int_t i = 0;
    Double_t min = +10000;
    Double_t max = -10000;
    while ((h = static_cast<TH1*>(next()))) {
      Color_t c = AliTrackletAODUtils::CentColor(i);
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

    if (!s->GetHists() || s->GetHists()->GetEntries() <= 0) {
      delete s;
      s = 0;
    }
    return s;
  }
}
//____________________________________________________________________
void AliTrackletdNdeta2::VisualizeWeights(Container* simList)
{
  DebugGuard g(fProc&kDebug,1,"Visualizing weights");
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
  TH2*     hp = GetH2(w,      "centPt");
  THStack* ef = new THStack(GetP2(simList,"etaWeight"),"x","effWeights","");
  THStack* ab = GetPdgStack(w, "abundance");
  THStack* st = GetPdgStack(w, "strangeness");
  SetCentColors(ef, c);
  SetCentColors(ab);
  SetCentColors(st);

  Double_t eMin = +1e9;
  Double_t eMax = -1e9;
  TIter next(ef->GetHists());
  TH1* h = 0;
  while ((h = static_cast<TH1*>(next()))) {
    eMin = TMath::Min(h->GetMinimum(), eMin);
    eMax = TMath::Max(h->GetMaximum(), eMax);
  }
  ef->SetMinimum(eMax);
  ef->SetMinimum(eMin);
  // ef->SetMaximum(1.02);
  TLegend* l = DrawInPad(fBody, 1, ef, "nostack p leg");
  ef->GetHistogram()->SetYTitle("Average weight");
  ef->GetHistogram()->SetXTitle("#eta");
  fBody->GetPad(1)->GetListOfPrimitives()->Remove(l);
  fBody->GetPad(1)->Modified();

  if (hp) {
    if (hp->GetNbinsY() > 1) {
      THStack* pt = new THStack(hp,   "y","pt","");      
      SetCentColors(pt);
      DrawInPad(fBody, 2, pt, "nostack");
      pt->GetHistogram()->SetYTitle("Weight");
      pt->GetHistogram()->SetXTitle("#it{p}_{T}");
    }
    else {
      TArrayD bins(hp->GetNbinsX()+1);
      bins[0] = hp->GetXaxis()->GetBinLowEdge(1);
      for (Int_t i = 1; i <= hp->GetNbinsX(); i++)
	bins[i] = hp->GetXaxis()->GetBinUpEdge(1);
      TH1* pt = new TH1D("pt","", bins.GetSize()-1,bins.GetArray());
      for (Int_t i = 1; i <= hp->GetNbinsX(); i++) {
	pt->SetBinContent(i, hp->GetBinContent(i,1));
	pt->SetBinError  (i, hp->GetBinError  (i,1));
      }
      pt->SetYTitle("Weight");
      pt->SetXTitle("Centrality [%]");
      pt->SetMarkerStyle(2);
      pt->SetMarkerColor(kRed+2);
      DrawInPad(fBody, 2, pt, "");      
    }
  }
  
  fBody->cd();
  if (l) {
    l->Draw();
    ModLegend(fBody, l, right, pp[1]->GetYlowNDC(),
	      .99, 1-fBody->GetTopMargin());
    fBody->Modified();
  }
  
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

  if (ab) {
    ab->GetHistogram()->SetYTitle("Weight");
    ab->GetHistogram()->SetXTitle("Centrality [%]");
  }
  if (st) {
    st->GetHistogram()->SetYTitle("Weight");
    st->GetHistogram()->SetXTitle("Centrality [%]");
  }
  p3->GetPad(1)->Modified();
  p3->GetPad(2)->Modified(); 
  
  PrintCanvas("Simulation weights","weights");
}

//____________________________________________________________________
void AliTrackletdNdeta2::StackMinMax(THStack* stack,
				     Double_t& min,
				     Double_t& max)
{
  if (!stack->GetHists()) return;
  min = +1e6;
  max = -1e6;
  TIter next(stack->GetHists());
  TH1*  hist = 0;
  while ((hist = static_cast<TH1*>(next()))) {
    if (hist->GetEntries() <= 0) continue;
    min = TMath::Min(min, hist->GetMinimum());
    max = TMath::Max(max, hist->GetMaximum());
  }
  stack->SetMinimum(min);
  stack->SetMaximum(max);
}

    
 
//____________________________________________________________________
void AliTrackletdNdeta2::VisualizeFinal(TDirectory* outDir, Int_t i)
{
  DebugGuard g(fProc&kDebug,1,"Visualizing final @ %d", i);
  TDirectory* dd  = outDir->GetDirectory(Form("final%dd", i));
  if (!dd) return;

  THStack* all = GetHS(dd, "full");
  TH1*     mid = GetH1(dd, "mid");
  TH1*     pub = GetH1(outDir, "published");
  if (!mid) return;

  Double_t min, max;
  StackMinMax(all, min, max);
  if (all->GetHists() && all->GetHists()->GetEntries() > 1 && min <= 0) min=2;
  if (mid->GetMinimum() <= 0) mid->SetMinimum(min);
  if (fProc&kDebug) all->GetHists()->Print();
  
  max = TMath::Max(max,mid->GetMaximum());
  min = TMath::Min(min,mid->GetMinimum());
  all->SetMinimum(.9*min);
  all->SetMaximum(1.2*max);
  mid->SetMinimum(.9*min);
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
  TLegend* l = DrawInPad(p1,0, pub, "logy grid same leg");
  ModLegend(p1, l, .4, .90, .99, .99);
  l->SetMargin(0.2);
  l->SetEntrySeparation(0.1);
  l->SetTextSize(0.04);
  
  l = DrawInPad(p2,0, all, "nostack logy grid leg");
  if (all && all->GetHistogram()) all->GetHistogram()->SetXTitle("#eta");

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

//____________________________________________________________________
void AliTrackletdNdeta2::VisualizeClosure(TDirectory* outDir, Int_t i)
{
  DebugGuard g(fProc&kDebug,1,"Visualizing closure @ %d", i);
  TDirectory* dd  = outDir->GetDirectory(Form("final%dd", i));
  if (!dd) return;

  THStack* all = GetHS(dd, "closure", false);
  if (!all) return;
  
  Double_t max = all->GetMaximum("nostack");
  Double_t min = all->GetMinimum("nostack");
  all->SetMinimum(.95*min);
  all->SetMaximum(1.05*max);
  
  ClearCanvas();
  Double_t right = .3;
  fBody->SetRightMargin(right);
  TLegend* l = DrawInPad(fBody,0, all, "grid nostack leg");
  ModLegend(fBody, l, 1-right+.01, fBody->GetBottomMargin(), .99,
	    1-fBody->GetTopMargin());
  l->SetMargin(0.2);
  l->SetEntrySeparation(0.1);
  l->SetTextSize(0.04);
  l->SetBorderSize(0);
  fBody->Modified();
  
  const char* what = (i == 3 ? "d^{3}N/(d#Deltad#etadIP_{z})" :
		      i == 2 ? "d^{2}N/(d#Deltad#eta)" :
		      i == 1 ? "dN/d#Delta" : "dN/d#Delta (k#equiv1)");
  PrintCanvas(Form("Closure #topbar %s", what), "closure");
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
  VisualizeParam("min#Delta_{tail}",      y, GetD(pars, "TailDelta"));
  VisualizeParam("max#Delta_{tail}",      y, GetD(pars, "TailMax"));
  VisualizeParam("abs.min#it{c}",         y, GetD(pars, "AbsMinCent"));
}
//____________________________________________________________________
void AliTrackletdNdeta2::VisualizeParams(Container* realSums,
					 Container* simSums)
{
  DebugGuard g(fProc&kDebug,1,"Visualizing parameters");
  ClearCanvas();
  fBody->Divide((fViz & kLandscape ? 3 : 1),
		(fViz & kLandscape ? 1 : 3), 0, 0);

  TVirtualPad* p1 = fBody->GetPad(1);  
  TVirtualPad* p2 = fBody->GetPad(2);  
  TVirtualPad* p3 = fBody->GetPad(3);
#if 0
  if (!(fViz & kLandscape)) {
    Double_t yr = .2;
    p1->SetPad(0,1-yr,1,1);
    p2->SetPad(0,(1-yr)/2,1,1-yr);
    p3->SetPad(0,0,1,(1-yr)/2);
  }
#endif
  
  // Post-processing stuff 
  fBody->cd(1);
  Double_t y = .80;
  TLatex* latex = new TLatex(.5, y, "Post-processing");
  latex->SetTextAlign(21);
  latex->SetTextSize(0.023/p1->GetHNDC());
  latex->SetNDC();
  latex->Draw();
  y -= 0.028/p1->GetHNDC();
  TString scaleM("");
  if (fProc & kUnitScale)    scaleM.Append(" unit");
  if (fProc & kAverageScale) scaleM.Append(" const");
  if (fProc & kEtaScale)     scaleM.Append(" dN/d#eta");
  if (fProc & kFullScale)    scaleM.Append(" d^{2}N/(d#etadIP_{z})");
  TString tmp = scaleM.Strip(TString::kBoth, ' ');
  VisualizeParam("Scaling of comb. bg.", y, tmp);	     
  VisualizeParam("min#alpha", y, fMinAlpha);
  VisualizeParam("max#alpha", y, fMaxAlpha);
  VisualizeParam("min#it{k}", y, fMinK);
  VisualizeParam("max#it{k}", y, fMaxK);
  VisualizeParam("Flatten c", y, (fProc & kCentAverage) != 0);
  if (fFudge != 1) VisualizeParam("Fudge",     y, fFudge);
  
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
  if (TMath::Abs(c1 - c2) < 1e-6) centName = "all";
  DebugGuard g(fProc&kDebug,1,"Visualizing bin %s", centName.Data());

  fLastBin.Form("%.1f#minus%.1f%%", c1, c2);
  fLastShort = centName;
  
  TDirectory* outDir = outTop->GetDirectory(centName);
  if (!outDir) {
    Warning("VisualizeBin", "Directory %s not found in %s",
	    centName.Data(), outTop->GetName());
    return false;
  }

  Printf("%5.1f - %5.1f%%", c1, c2);
  if (fViz & kSpecies) {
    VisualizeSpecies(GetC(simList, centName));
    if (fViz & kDeltas) 
      VisualizeSpeciesDelta(GetC(simList, centName), outDir);
    VisualizePrimary(GetC(GetC(simList, centName), "generated"));
  }
  if (fViz & kdNdetas) {
    for (Int_t i = 0; i < 4; i++) {
      if ((fProc & (1 << i)) == 0) continue;
      if (fViz & kDeltas)      VisualizeDelta(outDir, i);
      if (fViz & kDetails)     VisualizeDetails(outDir, i);
      if (fViz & kdNdetas)     VisualizeResult(outDir, i);
    }
  }
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletdNdeta2::VisualizeSpecies(Container* simCont)
{
  DebugGuard g(fProc&kDebug,1,"Visualizing species");
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

    THStack* all    = GetHS(species, "all");
    THStack* toPion = GetHS(species, "toPion");
    
    
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
Bool_t AliTrackletdNdeta2::VisualizeSpeciesDelta(Container* simCont,
						 TDirectory* outDir)
{
  DebugGuard g(fProc&kDebug,1,"Visualizing species deltas");
  if (!simCont) return true;
  ClearCanvas();
  Container*  comb = GetC(GetC(simCont, "combinatorics"), "specie", false);
  Container*  prim = GetC(GetC(simCont, "primaries"),     "specie", false);
  Container*  seco = GetC(GetC(simCont, "secondaries"),   "specie", false);
  Container*  cont[] = { comb, prim, seco, 0 };
  if (!comb || !prim || !seco) return true;
  
  const char* tit[]  = { "C' by primary mother's specie",
			 "P' by species",
			 "S' by primary mother's specie" };

  // --- Draw the delta distributions --------------------------------
  THStack* hp = 0;
  THStack* hs = 0;
  THStack* hc = 0;
  Double_t rr = .8;
  for (Int_t i = 0; i < 2; i++) {
    TString  sub = (i == 0 ? "mid" : "fwd");
    TPad*    pad = new TPad("pad","pad",i*rr/2,0,(i+1)*rr/2,1);
    pad->SetTopMargin(0.10);
    pad->SetRightMargin(0.01);
    pad->SetNumber(i+1);
    fBody->cd();
    pad->Draw();
    pad->Divide(1,3,0,0);
    pad->GetPad(1)->SetRightMargin(0.01);
    pad->GetPad(2)->SetRightMargin(0.01);
    pad->GetPad(3)->SetRightMargin(0.01);
    DrawInPad(pad, 1, hp = GetHS(GetC(prim, sub), "all"),
	      "nostack grid logy logx");
    DrawInPad(pad, 2, hs = GetHS(GetC(seco, sub), "all"),
	      "nostack grid logy logx");
    DrawInPad(pad, 3, hc = GetHS(GetC(comb, sub), "all"),
	      "nostack grid logy logx");
    hp->GetHistogram()->SetYTitle("Primaries");
    hs->GetHistogram()->SetYTitle("Secondaries");
    hc->GetHistogram()->SetYTitle("Combinatorial");
    hp->GetHistogram()->SetXTitle("#Delta");
    hs->GetHistogram()->SetXTitle("#Delta");
    hc->GetHistogram()->SetXTitle("#Delta");
    TLatex* txt = new TLatex(pad->GetLeftMargin()+
			     (1-pad->GetLeftMargin()-pad->GetRightMargin())/2,
			     1-pad->GetTopMargin()+.01, hp->GetTitle());
    txt->SetNDC();
    txt->SetTextAlign(21);
    txt->SetTextSize(0.07);
    DrawInPad(pad,0,txt,"");
  }
  // --- Make a legend -----------------------------------------------
  const TAxis& pdgs = PdgAxis();
  TLegend* l = new TLegend(rr, fBody->GetBottomMargin(),
			   .99, 1-fBody->GetTopMargin());
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  TLegendEntry* e = 0;
  for (Int_t i = 0; i < 3; i++) {
    THStack* tmp = i == 0 ? hs : i == 1 ? hc : hp;
    TIter next(tmp->GetHists());
    TH1*  h = 0;
    while ((h = static_cast<TH1*>(next()))) {
      Int_t   bin = h->GetUniqueID();
      TString lbl = pdgs.GetBinLabel(bin);
      Int_t   pdg = lbl.Atoi();
      TString nme;
      Color_t col;
      Style_t sty;
      PdgAttr(pdg, nme, col, sty);
      if (nme.EqualTo("0")) {
	nme = "All";
	sty = 24;
	col = kBlack;
	h->SetMarkerStyle(sty);
      }
      nme.Append("#rightarrowX_{ch}");
      TIter         nextE(l->GetListOfPrimitives());
      while ((e = static_cast<TLegendEntry*>(nextE()))) {
	if (nme.EqualTo(e->GetLabel())) break;
      }
      if (e) continue;
      e = l->AddEntry(lbl,nme,"p");
      e->SetMarkerStyle(sty);
      e->SetMarkerColor(col);
      e->SetLineColor  (col);
    }
  }
  fBody->cd();
  l->Draw();
  
  PrintCanvas(Form("Species #Delta #topbar %s", fLastBin.Data()),
	      Form("%s_species_delta", fLastShort.Data()));

  // --- Draw relative contributions to signal/background ------------
  ClearCanvas();
  fBody->SetTopMargin(0.2);
  fBody->Divide(2,3,0,0);
  gStyle->SetPaintTextFormat("9.7f");
  Bool_t drawLog = true;
  for (Int_t i = 0; i < 2; i++) {
    TString  sub = (i == 0 ? "mid" : "fwd");
    for (Int_t j = 0; j < 3; j++) {
      Container*   par   = (j == 0 ? prim :
			    j == 1 ? seco : comb);
      TString      parT  = (j == 0 ? "Primaries" :
			    j == 1 ? "Secondaries" : "Combinatorics");
      Int_t        padNo = j * 2 + i + 1;
      TVirtualPad* pad   = fBody->cd(padNo);
      if ((padNo % 2) == 0) pad->SetRightMargin(0.01);
      pad->SetTicks();
      pad->SetGridx();
      pad->SetGridy();
      THStack*     rhs   = GetHS(GetC(par, sub), "ratios");
      TObjLink*    ptr   = rhs->GetHists()->FirstLink();
      while (ptr) {
	ptr->SetOption("hist bar0");
	TH1* h = static_cast<TH1*>(ptr->GetObject());
	h->SetMarkerSize(2);
	ptr = ptr->Next();
      }
      rhs->SetMaximum(drawLog ? 2    : 1.1);
      rhs->SetMinimum(drawLog ? 1e-8 : 0);
      // Printf("Draw %s/%s in %d", sub.Data(), parT.Data(), padNo);
      DrawInPad(fBody, padNo, hp = rhs, Form("nostack %s",
					     drawLog ? "logy" : ""));
      hp->GetHistogram()->SetYTitle(parT);
    }
  }
  // --- Make a legend -----------------------------------------------
  l = new TLegend(fBody->GetLeftMargin(),
		  1-fBody->GetTopMargin(),
		  .99, .99);
  l->SetNColumns(2);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  e = l->AddEntry("dummy", "Signal #minus #Delta<1.5", "f");
  e->SetFillColor(kGreen+1);
  e->SetFillStyle(1001);
  e = l->AddEntry("dummy", "Background #minus 5<#Delta<25", "f");
  e->SetFillColor(kRed+1);
  e->SetFillStyle(1001);
  fBody->cd();
  l->Draw();
  
  PrintCanvas(Form("Species contribution #topbar %s", fLastBin.Data()),
	      Form("%s_species_contrib", fLastShort.Data()));


  // --- Draw selected particles -------------------------------------
  TH1*     lt1p       = GetH1(GetC(prim, "mid"), "totalIntegrals");
  TH1*     lt1s       = GetH1(GetC(seco, "mid"), "totalIntegrals");
  TH1*     lt1c       = GetH1(GetC(comb, "mid"), "totalIntegrals");
  TH1*     gt1p       = GetH1(GetC(prim, "fwd"), "totalIntegrals");
  TH1*     gt1s       = GetH1(GetC(seco, "fwd"), "totalIntegrals");
  TH1*     gt1c       = GetH1(GetC(comb, "fwd"), "totalIntegrals");
  Double_t totalA     = (lt1p->GetBinContent(1)+
			 lt1s->GetBinContent(1)+
			 lt1c->GetBinContent(1)+
			 gt1p->GetBinContent(1)+
			 gt1s->GetBinContent(1)+
			 gt1c->GetBinContent(1));
  
  // --- Draw relative contributions to signal/background ------------
  ClearCanvas();
  fBody->SetTopMargin(0.2);
  fBody->SetRightMargin(0.01);
  fBody->Divide(2,3,0,0);
  gStyle->SetPaintTextFormat("6.2f");
  drawLog = false;
  for (Int_t i = 0; i < 2; i++) {
    TString  sub = (i == 0 ? "mid" : "fwd");
    TString  bin = (i == 0 ? "|#eta|<1" : "|#eta|>1");
    for (Int_t j = 0; j < 3; j++) {
      Container*   par   = (j == 0 ? prim :
			    j == 1 ? seco : comb);
      TString      parT  = (j == 0 ? "Primaries" :
			    j == 1 ? "Secondaries" : "Combinatorics");
      Int_t        padNo = j * 2 + i + 1;
      TVirtualPad* pad   = fBody->cd(padNo);
      if ((padNo % 2) == 0) pad->SetRightMargin(0.10);
      else                  pad->SetLeftMargin(0.20);
      pad->SetTicks();
      pad->SetGridx();
      pad->SetGridy();
      TH1*         itg   = GetH1(GetC(par, sub), "totalIntegrals");
      THStack*     rhs   = GetHS(GetC(par, sub), "rows");
      TObjLink*    ptr   = rhs->GetHists()->FirstLink();
      while (ptr) {
	ptr->SetOption("hist bar0 text30");
	TH1* h = static_cast<TH1*>(ptr->GetObject());
	h->SetMarkerSize(3);
	ptr = ptr->Next();
      }
      rhs->SetMaximum(drawLog ? 200  : 110);
      rhs->SetMinimum(drawLog ? 1e-3 : 0);
      // Printf("Draw %s/%s in %d", sub.Data(), parT.Data(), padNo);
      DrawInPad(fBody, padNo, rhs, Form("nostack %s",
					drawLog ? "logy" : ""));
      rhs->GetHistogram()->SetYTitle(parT);
      rhs->GetHistogram()->GetXaxis()->SetLabelSize(0.03/pad->GetHNDC());
      rhs->GetHistogram()->GetYaxis()->SetLabelSize(0.03/pad->GetHNDC());
      rhs->GetHistogram()->GetYaxis()->SetTitleSize(0.03/pad->GetHNDC());
      rhs->GetHistogram()->GetYaxis()->SetNdivisions(207);
      rhs->GetHistogram()->GetYaxis()->SetTitleOffset(1);
      
      pad->cd();
      Double_t total = itg->GetBinContent(1);
      TLatex* txt = new TLatex(pad->GetLeftMargin()+.05,
			       .99, Form("%5.2f%% of all",
					 100*total/totalA));
      txt->SetTextAlign(13);
      txt->SetNDC();
      txt->SetTextSize(0.06);
      txt->Draw();

      txt = new TLatex(pad->GetLeftMargin()+.05,
		       .99-txt->GetTextSize(),
		       Form("Signal %5.2f%% of all",
			    100*itg->GetBinContent(2)/totalA));
      txt->SetTextAlign(13);
      txt->SetNDC();
      txt->SetTextSize(0.06);
      txt->SetTextColor(kGreen+1);
      txt->Draw();

      txt = new TLatex(pad->GetLeftMargin()+.05,
		       .99-2*txt->GetTextSize(),
		       Form("Background %5.2f%% of all",
			    100*itg->GetBinContent(3)/totalA));
      txt->SetTextAlign(13);
      txt->SetNDC();
      txt->SetTextSize(0.06);
      txt->SetTextColor(kRed+1);
      txt->Draw();
    }
  }
  // --- Make a legend -----------------------------------------------
  l = new TLegend(fBody->GetLeftMargin(),
		  1-fBody->GetTopMargin(),
		  .99, .99);
  l->SetNColumns(2);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  e = l->AddEntry("dummy", "Signal #minus #Delta<1.5", "f");
  e->SetFillColor(kGreen+1);
  e->SetFillStyle(1001);
  e = l->AddEntry("dummy", "Background #minus 5<#Delta<25", "f");
  e->SetFillColor(kRed+1);
  e->SetFillStyle(1001);
  fBody->cd();
  l->Draw();
  
  PrintCanvas(Form("Species strangeness #topbar %s", fLastBin.Data()),
	      Form("%s_species_strange", fLastShort.Data()));


    
}  
  

//____________________________________________________________________
Bool_t AliTrackletdNdeta2::VisualizePrimary(Container* simCont)
{
  DebugGuard g(fProc&kDebug,1,"Visualizing primaries");
  if (!simCont) return true;
  ClearCanvas();

  THStack*  all     = GetHS(GetC(simCont,"mix"),"all");
  THStack*  toPion  = GetHS(GetC(simCont,"mix"),"toPion");
  THStack*  toAll   = GetHS(GetC(simCont,"mix"),"toAll");
  TH2*      etaPt   = GetH2(simCont, "etaPt");
  Int_t     etaM    = etaPt->GetXaxis()->FindBin(-.5);
  Int_t     etaP    = etaPt->GetXaxis()->FindBin(+.5);
  TH1*      pt      = etaPt->ProjectionY("pt", etaM, etaP);
  pt->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T}");
  pt->GetYaxis()->SetTitleSize(0.08);
  pt->GetYaxis()->SetLabelSize(0.08);
  pt->GetYaxis()->SetTitleOffset(0.6);  
  pt->GetXaxis()->SetTitleSize(0.08);
  pt->GetXaxis()->SetLabelSize(0.08);
  pt->GetXaxis()->SetTitleOffset(0.6);  

  if (false) {
    TIter next(toPion->GetHists());
    TH1*  rat = 0;
    while ((rat = static_cast<TH1*>(next()))) {
      TGraphErrors* g =
	static_cast<TGraphErrors*>(rat->GetListOfFunctions()
				   ->FindObject(Form("%s_2760",
						     rat->GetName())));
      TF1* f = new TF1("fit", "pol0", -.5, +.5);
      rat->Fit(f, "Q0R+", "", -.5, +.5);
      Double_t re, r = RatioE(f->GetParameter(0), f->GetParError(0),
			      g->GetY()[0], g->GetEY()[0], re);
      Printf("%10s:  2760: %6.4f +/- %6.4f  Here: %6.4f +/- %6.4f  "
	     "Ratio: %6.4f +/- %6.4f",
	     rat->GetName(), g->GetY()[0], g->GetEY()[0],
	     f->GetParameter(0), f->GetParError(0), r, re);
    
      f->SetLineColor(rat->GetLineColor());
      f->SetLineStyle(7);
    }
  }
  
  TPad* p1 = new TPad("p1","p1",0,.3,1,1);
  p1->SetTopMargin(.01);
  p1->SetRightMargin(.01);
  fBody->cd();
  p1->Draw(); 
  p1->SetNumber(1);
  p1->Divide(1,3,0,0);
  p1->GetPad(1)->SetRightMargin(0.2);
  p1->GetPad(2)->SetRightMargin(0.2);
  p1->GetPad(3)->SetRightMargin(0.2);
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

  l = DrawInPad(p1,2,toAll, "nostack leg2 logy");
  toAll->GetHistogram()->GetYaxis()->SetTitle("Ratio to all");
  toAll->GetHistogram()->GetYaxis()->SetTitleSize(0.08);
  toAll->GetHistogram()->GetYaxis()->SetLabelSize(0.08);
  toAll->GetHistogram()->GetYaxis()->SetTitleOffset(0.6);  
  toAll->GetHistogram()->GetXaxis()->SetTitle("#eta");
  toAll->GetHistogram()->GetXaxis()->SetTitleSize(0.08);
  toAll->GetHistogram()->GetXaxis()->SetLabelSize(0.08);
  toAll->GetHistogram()->GetXaxis()->SetTitleOffset(0.6);  
  ModLegend(p1->GetPad(2), l,
	    1-p1->GetPad(2)->GetRightMargin(),
	    p1->GetPad(2)->GetBottomMargin(),
	    1-p1->GetPad(2)->GetTopMargin(),
	    .99);
  l->SetBorderSize(0);

  l = DrawInPad(p1,3,toPion, "nostack leg");
  toPion->GetHistogram()->GetYaxis()->SetTitle("Ratio to #pi");
  toPion->GetHistogram()->GetYaxis()->SetTitleSize(0.08);
  toPion->GetHistogram()->GetYaxis()->SetLabelSize(0.08);
  toPion->GetHistogram()->GetYaxis()->SetTitleOffset(0.6);  
  toPion->GetHistogram()->GetXaxis()->SetTitle("#eta");
  toPion->GetHistogram()->GetXaxis()->SetTitleSize(0.08);
  toPion->GetHistogram()->GetXaxis()->SetLabelSize(0.08);
  toPion->GetHistogram()->GetXaxis()->SetTitleOffset(0.6);  
  ModLegend(p1->GetPad(3), l,
	    1-p1->GetPad(3)->GetRightMargin(),
	    p1->GetPad(3)->GetBottomMargin(),
	    1-p1->GetPad(3)->GetTopMargin(),
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
  DebugGuard g(fProc&kDebug,1,"Visualizing delta in %s/%dD",
	       outTop->GetName(), dimen);
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
  THStack* all = GetHS(outDir,"all");
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

  
  THStack* ratios = GetHS(outDir,"ratios");
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
  TH2* scale     = GetH2(outDir,"scale");
  TH1* scaleProj = GetH1(outDir,"scaleProj");
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
Bool_t AliTrackletdNdeta2::VisualizeDetails(TDirectory* outTop,
					    Int_t       dimen)
{
  DebugGuard g(fProc&kDebug,1,"Visualize details for %s %dD",
	       outTop->GetName(), dimen);
  TDirectory* outDir = outTop->GetDirectory(Form("results%dd", dimen));
  if (!outDir) {
    Warning("VisualizeDelta", "Directory results%dd not found in %s",
	    dimen, outTop->GetName());
    return false;
  }
  Double_t tbase = 0.015;
  ClearCanvas();
  fBody->cd();
  fBody->SetTopMargin(0.01);
  fBody->SetRightMargin(0.01);
  fBody->SetBottomMargin(0.15);
  fBody->SetLeftMargin(0.15);
  fBody->Divide(2,6,0,0);
  
  const char* names[] = { "realM",  "simM",
			  "realC",  "simC",
			  "realB",  "simB",
			  "realS",  "simS",
			  "simA",   "fiducial",
			  "result", "simG" };
  const char* titles[] = { "M",     "M'",
			   "C",     "C'",
			   "#beta", "#beta'",
			   "S",     "S'",
			   "#alpha", "F'",
			   "R",      "G'" };

  for (Int_t i = 0; i < 12; i++) {
    const char*  name  = names[i];
    const char*  title = titles[i];
    TVirtualPad* pad   = fBody->GetPad(i+1);
    pad->SetLeftMargin(0.15);
    pad->SetRightMargin(0);
    pad->Divide(1,2,0,0);
    TVirtualPad* q2 = pad->GetPad(1); q2->SetRightMargin(0.15);
    TVirtualPad* q1 = pad->GetPad(2); q1->SetRightMargin(0.15);
    
    TH2* h2 = GetH2(outDir, Form("full/%s", name));
    TH1* h1 = GetH1(outDir, name);
    if (!h2 || !h1) {
      Warning("VisualizeDetails", "Didn't find full/%s (%p) or %s (%p)",
	      name, name);
      continue;
    }
    h2 = static_cast<TH2*>(h2->Clone());
    h1 = static_cast<TH1*>(h1->Clone());
    h2->SetDirectory(0);
    h1->SetDirectory(0);
    h2->SetXTitle("#eta");
    h1->SetXTitle("#eta");
    h2->SetYTitle(title);
    h1->SetYTitle("");
    TAxis*        axis[] = { h2->GetXaxis(), h2->GetYaxis(), h2->GetZaxis(),
			     h1->GetXaxis(), h1->GetYaxis(), 0 };
    TVirtualPad*  pads[] = { q2, q2, q2, q1, q1, 0 };
    TAxis**       pa     = axis;
    TVirtualPad** pq     = pads;
    while ((*pa)) {
      (*pa)->SetTitleSize(tbase/pad->GetHNDC()/(*pq)->GetHNDC());
      (*pa)->SetLabelSize(tbase/pad->GetHNDC()/(*pq)->GetHNDC());
      (*pa)->SetTitleOffset(0.4);
      (*pa)->SetNdivisions(207);
      pa++;
      pq++;
    }
    DrawInPad(pad, 1, h2, "colz");
    DrawInPad(pad, 2, h1, "");
    pad->Modified();
  }
  fBody->Modified();
  
  PrintCanvas(Form("Details #topbar %s", fLastBin.Data()),
	      Form("%s_details", fLastShort.Data()));
  return true;
}

//____________________________________________________________________
Bool_t AliTrackletdNdeta2::VisualizeResult(TDirectory* outTop,
					   Int_t       dimen)
{
  DebugGuard g(fProc&kDebug,1,"Visualize results for  %s %dD",
	       outTop->GetName(), dimen);
  TDirectory* outDir = outTop->GetDirectory(Form("results%dd", dimen));
  if (!outDir) {
    Warning("VisualizeDelta", "Directory results%dd not found in %s",
	    dimen, outTop->GetName());
    return false;
  }

  ClearCanvas();
  Double_t tbase = 0.03;
  Double_t yr = .2;
  Double_t yf = .2;
  TPad* p1 = new TPad("p1","p1",0,yf,1,1);
  p1->SetTopMargin(yr);
  p1->SetRightMargin(0.01);
  p1->SetLeftMargin(0.15);
  p1->SetBottomMargin(0);
  p1->SetTicks();
  fBody->cd();
  p1->Draw();
  p1->SetNumber(1);

  TPad* p2 = new TPad("p2","p2",0,0,1,yf);
  p2->SetTopMargin(0.01);
  p2->SetRightMargin(.01);
  p2->SetLeftMargin(0.15);
  p2->SetBottomMargin(0.20);
  p2->SetTicks();
  fBody->cd();
  p2->Draw();
  p2->SetNumber(2);
  
  THStack* all = GetHS(outDir, "all");
  TLegend* l = DrawInPad(fBody,1,all, "nostack leg2");
  all->GetHistogram()->SetXTitle("#eta");
  all->GetHistogram()->SetYTitle(ObsTitle());
  all->GetHistogram()->GetYaxis()->SetTitleOffset(1.7);
  all->GetHistogram()->GetYaxis()->SetTitleSize(tbase/(1-yf));
  all->GetHistogram()->GetYaxis()->SetLabelSize(tbase/(1-yf));
  all->GetHistogram()->GetYaxis()->SetNdivisions(205);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetMargin(0.12);
  l->SetEntrySeparation(0.1);
  l->SetHeader("R=G'/(M'-C')#times(M-C) #kern[2]{ } [C=kC'/M'#timesM]");
  // l->SetTextSize(0.04);
  // l->SetTextAlign(12);
  ModLegend(p1, l, p1->GetLeftMargin()-.01, 1-yr,
	    1-p1->GetRightMargin(), .99);

  THStack* ratios = GetHS(outDir, "ratios");
  FixMinMax(ratios, true);
  DrawInPad(fBody,2,ratios, "nostack");
  ratios->GetHistogram()->SetXTitle("#eta");
  ratios->GetHistogram()->SetYTitle("Ratios");
  ratios->GetHistogram()->GetYaxis()->SetTitleOffset(.45);
  ratios->GetHistogram()->GetXaxis()->SetTitleOffset(.7);
  ratios->GetHistogram()->GetYaxis()->SetTitleSize(tbase/yr);
  ratios->GetHistogram()->GetYaxis()->SetLabelSize(tbase/yr);
  ratios->GetHistogram()->GetXaxis()->SetTitleSize(tbase/yr);
  ratios->GetHistogram()->GetXaxis()->SetLabelSize(tbase/yr);
  ratios->GetHistogram()->GetYaxis()->SetNdivisions(205);
  
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
