#ifndef ALITRACKLETDNDETAPOST_H
#define ALITRACKLETDNDETAPOST_H
#include <AliTrackletdNdetaUtils.C>
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
#include <THashList.h>
#include <TStyle.h>
#include <TBrowser.h>
#include <TGraphAsymmErrors.h>
#include <TProfile.h>
#else
class TPad;
class TLatex;
class TObject;
class TSeqCollection;
class TH1;
class TH2;
class THStack;
class TCanvas;
class TVirtualPad;
class TFile;
class TAxis;
#endif

//====================================================================
struct AliTrackletdNdetaPost : public AliTrackletdNdetaUtils
{
  /**
   * Processing options 
   * 
   */
  enum {
    /** Use combinatorial background for real data */
    kRealComb    = 0x001,
    /** Use combinatorial background for simulated data */
    kSimComb     = 0x002,
    /** Draw dNch/deta */
    kdNdetas     = 0x004,
    /** Draw general information */
    kGeneral     = 0x008,
    /** Draw parameters */
    kParameters  = 0x010,
    /** Draw simulation information */
    kSimInfo     = 0x020,
    /** Draw species information */
    kSpecies     = 0x040,
    /** Draw delta information */   
    kDeltas      = 0x080,
    /** Draw backgrounds */
    kBackgrounds = 0x100,
    /** Draw alphas */
    kAlphas      = 0x200,
    /** MC closure test */
    kClosure     = 0x400,
    /** Print to PDF */
    kPDF         = 0x1000,
    /** Whether to pause after each plot */
    kPause       = 0x2000,
    /** Draw in landscape */
    kLandscape   = 0x4000,
    /** Default options */
    kDefault     = 0x33ff
  };
  enum {
    kTopBackground = kAzure-8
  };
  //_________________________________________________________________________
  /** 
   * @{
   * @name Service functions to check histogram consistencies 
   */
  /** 
   * Check if both axis have the same number of bins 
   * 
   * @param which Which axis is being checked 
   * @param a1    First axis
   * @param a2    Second axis 
   * 
   * @return true of both axis have equal number of bins 
   */
  static Bool_t CheckAxisNBins(const char*  which,
			       const TAxis* a1,
			       const TAxis* a2);
  /** 
   * Check axis limits (min,max) are the same 
   * 
   * @param which Which axis we're checking 
   * @param a1    First axis 
   * @param a2    Second axis 
   * 
   * @return true if axis min/max are indetical 
   */
  static Bool_t CheckAxisLimits(const char*  which,
				const TAxis* a1,
				const TAxis* a2);
  /** 
   * In case of non-uniform bins, check each bin edge to see if they
   * are the same
   *  
   * @param which Which axis we're checking  
   * @param a1    First axis 
   * @param a2    Second axis 
   * 
   * @return true if all bin edges are the same 
   */
  static Bool_t CheckAxisBins(const char*  which,
			      const TAxis* a1,
			      const TAxis* a2);
  /** 
   * Check if all bin labels (if specified) are the same 
   * 
   * @param which Which axis we're checking 
   * @param a1    First axis 
   * @param a2    Second axis 
   * 
   * @return True if all bin labels are the same 
   */
  static Bool_t CheckAxisLabels(const char*  which,
				const TAxis* a1,
				const TAxis* a2);
  /** 
   * Check that axis definitions are compatible 
   * 
   * @param which    Which axis we're checking 
   * @param a1       First axis 
   * @param a2       Second axis 
   * @param alsoLbls If true, also check labels 
   *  
   * @return true if the axes are compatible 
   */
  static Bool_t CheckAxis(const char*  which, 
			  const TAxis* a1,
			  const TAxis* a2,
			  Bool_t       alsoLbls);
  /** 
   * Check if two histograms are compatible by checking each defined
   * axis.
   * 
   * @param h1 First histogram 
   * @param h2 Second histogram 
   * 
   * @return true of they are compatible 
   */
  static Bool_t CheckConsistency(const TH1* h1, const TH1* h2);
  /* @} */
  //==================================================================
  /** The canvas to draw in */
  TCanvas* fCanvas;
  /** Top bar to put title in */
  TPad*    fTop;
  /** Main body of plots */
  TPad*    fBody;
  /** Whether we should be in landscape mode */
  Bool_t   fLandscape;
  /** If we should print to a PDF */
  Bool_t   fPDF;
  /** Whether we should pause */
  Bool_t   fPause;
  /** Cache of last title */
  TString  fLastTitle;
  /** Cache of centrality bin title */
  TString fLastBin;
  /** The page header text */
  TLatex*  fHeader;
  /** Cut on @f$ \Delta@f$ read from real data input file */
  Double_t fDeltaCut;
  /** Lower cut on @f$\Delta@f$ tail read from real data input file */
  Double_t fTailDelta;
  /** Highest value of @f$\Delta@f$ read from real data input file */
  Double_t fMaxDelta;
  /** The scalar on combinatorial background for real data */
  Double_t fCombinatoricsScale;
  /** Least value of @f$\alpha@f$ */
  Double_t fAlphaMin;
  /** Largest value of @f$\alpha@f$ */
  Double_t fAlphaMax;
  /** 
   * Constructor 
   */
  AliTrackletdNdetaPost();
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
   * @param landscape  Should it be landscape? 
   * @param pdf        Should we store in PDF?
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
   * @param title Title of this page 
   * @param size  Size of title 
   */
  void PrintCanvas(const TString& title, Float_t size=.7);
  //____________________________________________________________________
  /** 
   * Draw an object in a sub-pad 
   * 
   * @param c   Canvas (or pad)
   * @param pad Sub-pad number 
   * @param o   Object to draw 
   * @param opt Drawing options 
   */
  void DrawInPad(TVirtualPad* c, Int_t pad, TObject* o, Option_t* opt);
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
  void DrawParam(const char* name, Double_t& y,  const char* val);
  /** 
   * Draw a real valued paraemeter 
   * 
   * @param name  Name of parameter 
   * @param y     Current y coordinate 
   * @param val   Value 
   */  
  void DrawParam(const char* name, Double_t& y,  Double_t val);
  /** 
   * Draw an integer valued paraemeter 
   * 
   * @param name  Name of parameter 
   * @param y     Current y coordinate 
   * @param val   Value 
   */  
  void DrawParam(const char* name, Double_t& y,  Int_t val);
  /** 
   * Draw a boolean valued paraemeter 
   * 
   * @param name  Name of parameter 
   * @param y     Current y coordinate 
   * @param val   Value 
   */  
  void DrawParam(const char* name, Double_t& y,  Bool_t val);
  /** 
   * Draw all parameters in a container 
   * 
   * @param pars 
   * @param title 
   */
  void DrawParams(Container* pars, const char* title);
  /** 
   * Draw parameters from both real and simulated analysis 
   * 
   * @param realSums Real data 
   * @param simSums  simulated data 
   */
  void DrawParams(Container* realSums, Container* simSums);
  /* @} */
  //____________________________________________________________________
  /** 
   * @{ 
   * @name General information 
   */
  /** 
   * Make a stack of histograms for real and simulated data. 
   * 
   * @param name     Name of stack 
   * @param title    Title of stack
   * @param dataList Container of real data 
   * @param simList  Container of simulated data 
   * @param dataOpt  Options for real data 
   * @param simOpt   Options for simulated data 
   * 
   * @return The created stack 
   */
  THStack* Make2Stack(const char*      name,
		      const char*      title,
		      Container*       dataList,
		      Container*       simList,
		      Option_t*        dataOpt="",
		      Option_t*        simOpt="");
  /** 
   * Draw general information 
   * 
   * @param realSums Sums from real data 
   * @param simSums  Sums from simulated data 
   * @param realRess Results from real data 
   * @param simRess  Results from simulated data 
   */
  void
  DrawGeneral(Container* realSums, Container* simSums,
	      Container* realRess, Container* simRess);
  /* @} */
  //____________________________________________________________________
  /** 
   * Functions related to single centrality bins 
   */
  /** 
   * Service function to find @f$\Delta@f$ distribution 
   * 
   * @param ress    Results 
   * @param sub     Sub set 
   * @param scaled  Whether to find the scaled value 
   * 
   * @return Histogram or null 
   */
  TH1* FindDelta(Container* ress, const char* sub, Bool_t scaled=false);
  /** 
   * utility function to find a sub-element 
   * 
   * @param ress Container 
   * @param sub  Sub-path 
   * @param name Name 
   * 
   * @return Histogram or null
   */
  TH1* FindSub(Container* ress, const char* sub, const char* name);
  /** 
   * Draw species information from simulated data 
   * 
   * @param mcRess Results from simulated data 
   */
  void DrawSpecies(Container* mcRess);
  /** 
   * Draw @f$\Delta@f$ distributions for a single bin 
   * 
   * @param realRess 
   * @param simRess 
   */
  void DrawDeltas(Container* realRess, Container* simRess);
  /** 
   * Draw background estimates in a single bin for a single estimate 
   * 
   * @param c     Top pad 
   * @param pad   Sub pad 
   * @param ress  Results 
   * @param name  Name of the results 
   */
  void DrawBackground(TVirtualPad* c,
		      Int_t        pad,
		      Container*   ress,
		      const char*  name,
		      const char*  pretty);
  /** 
   * Draw all background estimates for a centrality bin 
   * 
   * @param realRess Results from real data 
   * @param simRess  Results from simulated data 
   */
  void DrawBackground(Container* realRess, Container* simRess);
  /** 
   * Cut @f$\alpha@f$ edges between fMinAlpha and fMaxAlpha
   * 
   * @param alpha Alpha histogram to modify 
   * 
   * @return @a alpha
   */
  TH2* CutAlpha(TH2* alpha);
  /** 
   * Get the alpha and cut off edges 
   * 
   * @param c    Container 
   * @param name Sub container name 
   * @param mc   MC or not 
   * 
   * @return Copy of the retrieved alpha
   */
  TH2* GetAlpha(Container* c, const char* name, Bool_t sel);
  /** 
   * Draw background estimates in a single bin for a single estimate 
   * 
   * @param c     Top pad 
   * @param pad   Sub pad 
   * @param ress  Results 
   * @param name  Name of the results 
   */
  void DrawAlpha(TVirtualPad* c,
		 Int_t        pad,
		 Container*   ress,
		 const char*  name,
		 const char*  pretty);
  /** 
   * Draw all background estimates for a centrality bin 
   * 
   * @param realRess Results from real data 
   * @param simRess  Results from simulated data 
   */
  void DrawAlpha(Container* realRess, Container* simRess);
  /** 
   * Draw specific stuff from the simulation
   * 
   * @param simRess Container of results from simulated data
   */
  void DrawSim(Container* simSums, Container* simRess);
  /** 
   */
  //__________________________________________________________________
  /** 
   * @{ 
   * @name dN/deta stuff 
   */
  /** 
   * The observable title 
   * 
   * @return Observable title 
   */
  const char* ObsTitle() const
  {
    return "\\mathrm{d}N_{\\mathrm{ch}}/\\mathrm{d}\\eta";
  }
  /** 
   * Get the background from a container.  
   * 
   * @param comb If true, get combinatorial background 
   * @param list The container to get background from 
   * 
   * @return Background histogram
   */
  TH2* GetBg(Bool_t comb, Container* list);
  /** 
   * Correct measured distribution for background.  If we do not use
   * MC-labels then we simply substract the background.  
   *
   * If we use MC labels, then multiply the measures signal by
   * @f$ 1-\beta@f$, where @f$ \beta@f$ is 
   *
   * @f[
   *   C/M\prime
   * @f]
   *
   * where @f$ C@f$ is the combinatorial background, and @f$ M\prime@f$
   * is the measured distribution in simulations.
   * 
   * @param comb     true if we use combinatorial background 
   * @param meas     The measured distribution to correct (in-place)
   * @param bg       The background (or combinatorial background)
   * @param simMeas  The measured distribution (in MC)
   * 
   * @return Pointer to @f$ meas@f$ after correcting 
   */
  TH2* CorrectSignal(Bool_t comb, TH2* meas, TH2* bg, TH2* simMeas);
  /** 
   * Calculates the result as 
   *
   * @f[
   *  R &=& \frac{P}{(1-\beta\prime)M\prime} (1-\beta)M
   * @f] 
   *
   * where 
   * 
   * - @f$ P@f$ is the input primary particle distribution 
   * - @f$ \beta\prime@f$ is background correction for simulated data 
   * - @f$ M\prime@f$ is the observed distribution in simulations 
   * - @f$ \beta@f$ is background correction for real data 
   * - @f$ M@f$ is the observed distribution in real data 
   * 
   * In case we use MC-labels for the background correction in
   * simulations, we have that
   *
   * @f[
   *   (1-\beta\prime)M\prime = (1-C\prime/M\prime)M\prime = M\prime-C\prime
   * @f] 
   * 
   * where 
   *
   * - @f$ C\prime@f$ is the observed number of tracklets for which the 
   *   constituent clusters have different parent simulated traks.
   *
   * In case we use MC-labels for the background correction for real
   * data, we have that
   *
   * @f[ 
   *  (1-\beta)M = (1-k\beta\prime)M = (1-kC/M\prime)M
   * @f] 
   *
   * where @f$ k@f$ is a constant scaling factor (1.3). 
   *
   * If we use injected clusters as an estimator for the background
   * correction in simulations, we have 
   *
   * @f[
   *  (1-\beta\prime)M\prime = (1-B\prime/M\prime)M\prime = M\prime-B\prime
   * @f] 
   *
   * where 
   *
   * - @f$ B\prime@f$ is the observed number of tracklets when injecting
   *   fake clusters, scaled so that the corresponding @f$\Delta@f$
   *   distribution match the measured (in simulations) @f$\Delta@f$
   *   distribution in the tails (@f$\Delta\in[5,20]@f$)
   *
   * The same expression holds for the case when we use injected
   * clusters to estimate the background in real data.
   *
   *  (1-\beta)M = (1-B/M)M = M-B
   *
   * where @f$ B@f$ has the same meaning as @f$ B\prime@f$ except we use
   * real data. 
   *
   * Note, if we use MC-labels for both the simulated and real data
   * backgrounds, and @f$ k=1@f$ we find that 
   *
   * @f[
   *  R = \frac{P}{M\prime}M
   * @f] 
   *
   * @param realList  List of real data histograms 
   * @param realComb  Whether to use MC-labels for real data 
   * @param simList   List of simulated data histograms 
   * @param simComb   Whether to use MC-labels for simulated data 
   * @param color     Color to use 
   * @param stack     Stack to add results to 
   */
  void DrawdNdeta(Container* realRess,
		  Bool_t     dataComb,
		  Container* simRess,
		  Bool_t     simComb,
		  Color_t    color=kBlack,
		  THStack*   stack=0);
  /** 
   * Process a single bin 
   * 
   * @param what     Bit mask of processing options 
   * @param bin      Centrality bin number 
   * @param c1       Least centrality 
   * @param c2       Largest centrality 
   * @param realSums Sums from real data 
   * @param simSums  Sums from simulated data 
   * @param realRess Results from real data 
   * @param simRess  Results from simulated data 
   * @param stack    Possible stack to add dN/deta to 
   */
  void ProcessBin(UInt_t     what,     Int_t      bin,
		  Double_t   c1,       Double_t   c2,
		  Container* realSums, Container* simSums,
		  Container* realRess, Container* simRess,
		  THStack*   stack=0);
  //====================================================================
  /** 
   * Run it 
   * 
   * @param dataName File from real data 
   * @param simName  File from simulated data 
   * @param what     Bit mask of processing options 
   */
  void Run(UInt_t      what     = kDefault, 
	   UShort_t    maxBins  = 9,
	   const char* dataName = "data.root",
	   const char* simName  = "sim.root");
};

//____________________________________________________________________
Bool_t AliTrackletdNdetaPost::CheckAxisNBins(const char*  which,
					     const TAxis* a1,
					     const TAxis* a2)
{
  if (a1->GetNbins() != a2->GetNbins()) {
    ::Warning("CheckAxisNBins", "Incompatible number %s bins: %d vs %d",
	      which, a1->GetNbins(), a2->GetNbins());
    return false;
  }
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaPost::CheckAxisLimits(const char*  which,
					      const TAxis* a1,
					      const TAxis* a2)
{
  if (!TMath::AreEqualRel(a1->GetXmin(), a2->GetXmin(),1.E-12) ||
      !TMath::AreEqualRel(a1->GetXmax(), a2->GetXmax(),1.E-12)) {
    Warning("CheckAxisLimits",
	    "Limits of %s axis incompatible [%f,%f] vs [%f,%f]", which,
	    a1->GetXmin(), a1->GetXmax(), a2->GetXmin(), a2->GetXmax());
    return false;
  }
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaPost::CheckAxisBins(const char*  which,
					    const TAxis* a1,
					    const TAxis* a2)
{
  const TArrayD * h1Array = a1->GetXbins();
  const TArrayD * h2Array = a2->GetXbins();
  Int_t fN = h1Array->fN;
  if ( fN == 0 ) return true;
  if (h2Array->fN != fN) {
    // Redundant?
    Warning("CheckAxisBins", "Not equal number of %s bin limits: %d vs %d",
	    which, fN, h2Array->fN);
    return false;
  }
  else {
    for (int i = 0; i < fN; ++i) {
      if (!TMath::AreEqualRel(h1Array->GetAt(i),h2Array->GetAt(i),1E-10)) {
	Warning("CheckAxisBins",
		"%s limit # %3d incompatible: %f vs %f",
		which, i, h1Array->GetAt(i),h2Array->GetAt(i));
	return false;
      }
    }
  }
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaPost::CheckAxisLabels(const char*  which,
					      const TAxis* a1,
					      const TAxis* a2)
{
  // check that axis have same labels
  THashList *l1 = (const_cast<TAxis*>(a1))->GetLabels();
  THashList *l2 = (const_cast<TAxis*>(a2))->GetLabels();
  
  if (!l1 && !l2) return true;
  if (!l1 ||  !l2) {
    Warning("CheckAxisLabels", "Difference in %s labels: %p vs %p",
	    which, l1, l2);
    return false;
  }
  // check now labels sizes  are the same
  if (l1->GetSize() != l2->GetSize()) {
    Warning("CheckAxisLabels", "Different number of %s labels: %d vs %d",
	    which, l1->GetSize(), l2->GetSize());
    return false;
  }
  for (int i = 1; i <= a1->GetNbins(); ++i) {
    TString label1 = a1->GetBinLabel(i);
    TString label2 = a2->GetBinLabel(i);
    if (label1 != label2) {
      Warning("CheckAxisLabels", "%s label # %d not the same: '%s' vs '%s'",
	      which, i, label1.Data(), label2.Data());
      return false;
    }
  }
    
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaPost::CheckAxis(const char*  which, 
					const TAxis* a1,
					const TAxis* a2,
					Bool_t       alsoLbls)
{
  if (!CheckAxisNBins (which, a1, a2)) return false;
  if (!CheckAxisLimits(which, a1, a2)) return false;
  if (!CheckAxisBins  (which, a1, a2)) return false;
  if (alsoLbls && !CheckAxisLabels(which, a1, a2)) return false;
  return true;
}
//____________________________________________________________________
Bool_t AliTrackletdNdetaPost::CheckConsistency(const TH1* h1, const TH1* h2)
{
  // Check histogram compatibility
  if (h1 == h2) return true;
    
  if (h1->GetDimension() != h2->GetDimension() ) {
    Warning("CheckConsistency",
	    "%s and %s have different dimensions %d vs %d",
	    h1->GetName(), h2->GetName(),
	    h1->GetDimension(), h2->GetDimension());
    return false;
  }
  Int_t dim = h1->GetDimension(); 
    
  Bool_t alsoLbls = (h1->GetEntries() != 0 && h2->GetEntries() != 0);
  if (!CheckAxis("X", h1->GetXaxis(), h2->GetXaxis(), alsoLbls)) return false;
  if (dim > 1 &&
      !CheckAxis("Y", h1->GetYaxis(), h2->GetYaxis(), alsoLbls)) return false;
  if (dim > 2 &&
      !CheckAxis("Z", h1->GetZaxis(), h2->GetZaxis(), alsoLbls)) return false;
    
  return true;
}

//____________________________________________________________________
AliTrackletdNdetaPost::AliTrackletdNdetaPost()
  : AliTrackletdNdetaUtils(),
    fCanvas(0),
    fTop(0),
    fBody(0),
    fLandscape(false),
    fPDF(false),
    fPause(false),
    fLastTitle(""),
    fLastBin(""),
    fHeader(0),
    fDeltaCut(0),
    fTailDelta(0),
    fMaxDelta(0),
    fCombinatoricsScale(1.3),
    fAlphaMin(0),
    fAlphaMax(2.5)    
{}
//====================================================================  
void AliTrackletdNdetaPost::ClearCanvas()
{
  fTop->Clear();
  fTop->SetNumber(1);
  fTop->SetFillColor(kTopBackground);
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

  fCanvas->cd();
}
//____________________________________________________________________
void AliTrackletdNdetaPost::CreateCanvas(const TString& outputName)
{
  Int_t h    = 1000;
  Int_t w    = h / TMath::Sqrt(2);
  if (fLandscape) {
    Int_t t = h;
    h       = w;
    w       = t;
  }

  fCanvas = new TCanvas("c",outputName,w,h);
  fCanvas->SetFillColor(0);
  fCanvas->SetBorderSize(0);
  fCanvas->SetBorderMode(0);

  if (fPDF) 
    fCanvas->Print(Form("%s[", outputName.Data()),
		   Form("pdf %s", fLandscape ? "Landscape" : ""));
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
void AliTrackletdNdetaPost::CloseCanvas()
{
  if (fPDF && fCanvas) {
    fCanvas->Print(Form("%s]", fCanvas->GetTitle()),
		   Form("pdf %s Title:%s",
			fLandscape ? "Landscape" : "",
			fLastTitle.Data()));
  }
  if (fCanvas) fCanvas->Close();
  fCanvas = 0;
}
//____________________________________________________________________
void AliTrackletdNdetaPost::PrintCanvas(const TString& title, Float_t size)
{
  if (fTop) {
    fTop->cd();
    fHeader->SetTextSize(size);
    fHeader->DrawLatex(.5,.5,title);
  }
  fCanvas->Modified();
  fCanvas->Update();
  fCanvas->cd();

  if (fPDF) {
    TString tit;
    tit.Form("pdf %s Title:%s",
	     fLandscape ? "Landscape" : "",
	     title.Data());
    fCanvas->Print(fCanvas->GetTitle(), tit);
  }
  fLastTitle = title;
  if (fPause) fCanvas->WaitPrimitive();
}
//====================================================================
TFile* AliTrackletdNdetaPost::OpenFile(const char* filename)
{
  TFile* file = TFile::Open(filename, "READ");
  if (!file) {
    Warning("OpenFile", "Failed to open \"%s\"", filename);
    return 0;
  }
  return file;
}
//____________________________________________________________________
TH1* AliTrackletdNdetaPost::SetAttr(TH1* h,
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
void AliTrackletdNdetaPost::DrawInPad(TVirtualPad* c,
				      Int_t        pad,
				      TObject*     o,
				      Option_t*    opt)
{
  if (!o) {
    Warning("", "Nothing to draw in pad %d", pad);
    return;
  }
  TVirtualPad* p = c->cd(pad);
  TString option(opt);
  option.ToLower();
  if (option.Contains("logx")) { p->SetLogx(); option.ReplaceAll("logx",""); }
  if (option.Contains("logy")) { p->SetLogy(); option.ReplaceAll("logy",""); }
  if (option.Contains("logz")) { p->SetLogz(); option.ReplaceAll("logz",""); }
  if (option.Contains("grid")) { p->SetGridx(); p->SetGridy();
    option.ReplaceAll("grid",""); }
  Int_t leg = 0;
  if (option.Contains("leg2")) { leg = 2; option.ReplaceAll("leg2",""); }
  if (option.Contains("leg")) {  leg = 1; option.ReplaceAll("leg",""); }
  o->Draw(opt);
  if (leg) {
    TLegend* l = p->BuildLegend(0.5, 0.73, .98, .98);
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
}
//====================================================================
void AliTrackletdNdetaPost::DrawParam(const char* name,
				      Double_t&   y,
				      const char* val)
{
  TLatex* ln = new TLatex(.49, y, name);
  ln->SetTextAlign(31);
  ln->SetTextSize(0.04);
  ln->SetNDC();
  ln->Draw();
  TLatex* lv = new TLatex(.51, y, val);
  lv->SetTextAlign(11);
  lv->SetTextSize(0.04);
  lv->SetNDC();
  lv->Draw();
  y -= 0.045;
}
//____________________________________________________________________  
void AliTrackletdNdetaPost::DrawParam(const char* name,
				      Double_t&   y,
				      Double_t    val)
{
  DrawParam(name, y, Form("%f", val));
}
//____________________________________________________________________  
void AliTrackletdNdetaPost::DrawParam(const char* name,
				      Double_t&   y,
				      Int_t       val)
{
  DrawParam(name, y, Form("%d", val));
}
//____________________________________________________________________  
void AliTrackletdNdetaPost::DrawParam(const char* name,
				      Double_t&   y,
				      Bool_t      val)
{
  DrawParam(name, y, val ? "yes" : "no");
}
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawParams(Container* pars, const char* title)
{
  // c->Clear();
  Double_t y = .9;
  TLatex* latex = new TLatex(.5, y, title);
  latex->SetTextAlign(21);
  latex->SetTextSize(0.04);
  latex->SetNDC();
  latex->Draw();
  y -= 0.045;
  if (!pars) return;
  Int_t mode = GetI(pars, "RecMode");
  DrawParam("\\hbox{Reconstruction mode}",y,
	    Form("%s%s%s",
		 (mode & 0x2 ? "recon" : ""),
		 (mode & 0x4 ? " inj" : ""),
		 (mode & 0x8 ? " rot" : "")));
  DrawParam("\\hbox{Scale by }\\sin^2(\\theta)",y, GetB(pars, "ScaleDTheta"));
  DrawParam("\\delta\\phi\\hbox{ shift}", y, GetD(pars, "DPhiShift"));
  DrawParam("\\hbox{Shifted }\\delta\\phi\\hbox{ cut}",
	    y, GetD(pars, "ShiftedDPhiCut"));
  DrawParam("\\hbox{Scaled }\\delta\\theta\\hbox{ cut}",
	    y, GetD(pars, "ScaledDThetaCut"));
  DrawParam("\\Delta\\hbox{ cut}",        y, GetD(pars, "DeltaCut"));
  DrawParam("\\mathrm{max}\\Delta",       y, GetD(pars, "MaxDelta"));
  DrawParam("\\mathrm{min}\\Delta_{\\mathrm{tail}}", y,
	    GetD(pars,"TailDelta"));
  DrawParam("d\\theta\\hbox{ window}",    y, GetD(pars, "DThetaWindow"));
  DrawParam("d\\phi\\hbox{ window}",      y, GetD(pars, "DPhiWindow"));
  DrawParam("\\phi\\hbox{ overlap cut}",  y, GetD(pars, "PhiOverlapCut"));
  DrawParam("z-\\eta\\hbox{ overlap cut}",y, GetD(pars, "ZEtaOverlapCut"));
  DrawParam("\\phi\\hbox{ rotation}",     y, TMath::RadToDeg()*
	    GetD(pars, "PhiRotation"));
}
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawParams(Container* realSums, Container* simSums)
{
  ClearCanvas();
  fBody->Divide(1,2,0,0);
  fBody->cd(1);
  DrawParams(GetC(realSums, "parameters"), "Real data");
  fBody->cd(2);
  DrawParams(GetC(simSums, "parameters"), "Simulated data");
  PrintCanvas("Parameters");
}
//====================================================================
THStack* AliTrackletdNdetaPost::Make2Stack(const char*      name,
					   const char*      title,
					   Container*       dataList,
					   Container*       simList,
					   Option_t*        dataOpt,
					   Option_t*        simOpt)
{
  TString  nme(name);
  THStack* stack = new THStack(name, title);
  TH1*     data  = GetH1(dataList, name);
  TH1*     sim   = GetH1(simList,  name);
  TString dtit(data->GetTitle());
  if (dtit.Contains("\\")) dtit.Form("%s\\hbox{ - real}", data->GetTitle());
  else                     dtit.Form("%s - real", data->GetTitle());
  data->SetTitle(dtit);
  TString stit(sim->GetTitle());
  if (stit.Contains("\\")) stit.Form("%s\\hbox{ - sim.}", sim->GetTitle());
  else                     stit.Form("%s - sim.", sim->GetTitle());
  sim->SetTitle(stit);
  data->SetBarOffset(0.1); // @todo OK
  data->SetBarWidth (0.4); // @todo OK
  sim ->SetBarOffset(0.5); // @todo OK
  sim ->SetBarWidth (0.4); // @todo OK
  stack->Add(data, dataOpt);
  stack->Add(sim,  simOpt);
  return stack;
}
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawGeneral(Container* realSums,
					Container* simSums,
					Container* realRess,
					Container* simRess)
{
  gStyle->SetOptStat(0);
  Container* dataC  = realRess;
  Container* simC   = simRess;
  THStack* ipz    = Make2Stack("ipz",   "\\mathrm{IP}_{z}", dataC, simC);
  THStack* cent   = Make2Stack("cent",  "Centrality [%]",   dataC, simC);
  THStack* status = Make2Stack("status","Task status",      dataC, simC,
			       "B text90", "B text90");
  THStack* all0   = Make2Stack("allClusters0", "All clusters layer 0",
			       dataC, simC, "colz", "box");
  THStack* all1   = Make2Stack("allClusters1", "All clusters layer 1",
			       dataC, simC, "colz", "box");
  THStack* used0  = Make2Stack("usedClusters0", "Used clusters layer 0",
			       dataC, simC, "colz", "box");
  THStack* used1  = Make2Stack("usedClusters1", "Used clusters layer 1",
			       dataC, simC, "colz", "box");

  TH1* mcStatus = GetH1(simRess, "statusMC");
  mcStatus->SetMaximum(1.2*mcStatus->GetMaximum());
  mcStatus->SetMarkerSize(2); // @todo OK
  ClearCanvas();
  fBody->Divide(2,4);
  DrawInPad(fBody,1,ipz,     "nostack leg");
  DrawInPad(fBody,2,cent,    "nostack leg");
  DrawInPad(fBody,3,status,  "nostack hist text90 leg");
  DrawInPad(fBody,4,mcStatus,"hist text90");
  DrawInPad(fBody,5,all0,    "nostack");
  DrawInPad(fBody,6,all1,    "nostack");
  DrawInPad(fBody,7,used0,   "nostack");
  DrawInPad(fBody,8,used1,   "nostack");
  PrintCanvas("General information");
}
//====================================================================
TH1* AliTrackletdNdetaPost::FindDelta(Container*  ress,
				      const char* sub,
				      Bool_t      scaled)
{
  if (!ress) return 0;
  Container* subCont = GetC(ress, sub);
  if (!subCont) return 0;
  TH1* h = GetH1(subCont, Form("delta%s", scaled ? "Bg" : ""));
  if (!h) return 0;
  if (!scaled) return h;
  TString tit(h->GetTitle());
  if (tit.Contains("\\")) tit.Append(Form("\\hbox{ - %s}",sub));
  else                    tit.Append(Form(" - %s", sub));
  return h;
}
//____________________________________________________________________
TH1* AliTrackletdNdetaPost::FindSub(Container*  ress,
				    const char* sub,
				    const char* name)
{
  if (!ress)  return 0;
  TObjArray*  tokens = TString(sub).Tokenize("/");
  Container*  subCont = ress;
  TIter       next(tokens);
  TObjString* dir = 0;
  while ((dir = static_cast<TObjString*>(next()))) {
    subCont = GetC(ress, dir->GetName());
    if (!subCont) break;
  }
  TH1* ret = 0;
  if (subCont) ret = GetH1(subCont, name);
  tokens->Delete();
  return ret;
}
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawSpecies(Container* mcRess)
{
  THStack* particlePdgs = new THStack("particlePdg","Particle species");
  particlePdgs->Add(GetH1(mcRess, "allPrimaryPdg"));
  particlePdgs->Add(GetH1(mcRess, "cutPrimaryPdg"));  
  particlePdgs->Add(GetH1(mcRess, "allSecondaryPdg"));
  particlePdgs->Add(GetH1(mcRess, "cutSecondaryPdg"));
  THStack* parentPdgs = new THStack("parentPdg","Parent species");
  parentPdgs->Add(GetH1(mcRess, "allPrimaryParentPdg"));
  parentPdgs->Add(GetH1(mcRess, "cutPrimaryParentPdg"));  
  parentPdgs->Add(GetH1(mcRess, "allSecondaryParentPdg"));
  parentPdgs->Add(GetH1(mcRess, "cutSecondaryParentPdg"));
  
  ClearCanvas();
  fBody->Divide(1,2,0,0);
  fBody->GetPad(1)->SetRightMargin(0.01);
  fBody->GetPad(2)->SetRightMargin(0.01);

  DrawInPad(fBody, 1, particlePdgs, "nostack logy leg2 grid");
  DrawInPad(fBody, 2, parentPdgs,   "nostack logy leg2 grid");

  particlePdgs->GetHistogram()->GetXaxis()->LabelsOption("v");
  parentPdgs  ->GetHistogram()->GetXaxis()->LabelsOption("v");
  particlePdgs->SetMaximum(1.4*particlePdgs->GetMaximum("nostack"));
  parentPdgs  ->SetMaximum(1.4*parentPdgs  ->GetMaximum("nostack"));
  
  fBody->GetPad(1)->Modified();
  fBody->GetPad(2)->Modified();

  PrintCanvas(Form("%s\\hbox{ - species}", fLastBin.Data()));
}
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawDeltas(Container* realRess, Container* simRess)
{

  ClearCanvas();
  fBody->Divide(1,2,0,0);
  fBody->GetPad(1)->SetRightMargin(0.01);
  fBody->GetPad(2)->SetRightMargin(0.01);
      
  
  THStack* orig = new THStack("orig", "\\Delta");
  orig->Add(FindDelta(realRess, "data",                     false));
  orig->Add(FindDelta(realRess, "injection",                false));
  orig->Add(FindDelta(simRess,  "data",                     false));
  orig->Add(FindDelta(simRess,  "injection",                false));
  orig->Add(FindDelta(simRess,  "primaries",                false));
  orig->Add(FindDelta(simRess,  "combinatorics",            false));
  orig->Add(FindDelta(simRess,  "secondaries",              false));
  orig->Add(FindDelta(simRess,  "uncorrelatedCombinatorics",false));  
  
  THStack* scaled = new THStack("scaled", "\\Delta (scaled)");
  scaled->Add(FindDelta(realRess, "data",                     false));
  scaled->Add(FindDelta(realRess, "injection",                true));
  scaled->Add(FindDelta(simRess,  "data",                     false));
  scaled->Add(FindDelta(simRess,  "injection",                true));
  scaled->Add(FindDelta(simRess,  "primaries",                false));
  scaled->Add(FindDelta(simRess,  "combinatorics",            true));
  scaled->Add(FindDelta(simRess,  "secondaries",              false));
  scaled->Add(FindDelta(simRess,  "uncorrelatedCombinatorics",true));

  Double_t max = 1.5*orig->GetMaximum("nostack");
  TH1* templ = FindDelta(realRess, "data", false);
  TH2* frame = new TH2D("frame","",
			templ->GetNbinsX(),
			templ->GetXaxis()->GetXmin(),
			templ->GetXaxis()->GetXmax(),
			500, 0.001, max);
  frame->SetXTitle(templ->GetXaxis()->GetTitle());
  frame->SetYTitle(templ->GetYaxis()->GetTitle());
  frame->SetStats(0);
  frame->SetDirectory(0);
  TLine* cuts = new TLine(fDeltaCut, 0.001, fDeltaCut, .99*max);
  cuts->SetLineColor(kGreen+1);
  frame->GetListOfFunctions()->Add(cuts);
  TLine* cutt = new TLine(fTailDelta, 0.001, fTailDelta, .99*max);
  cutt->SetLineColor(kRed+1);
  frame->GetListOfFunctions()->Add(cutt);


  TH1* f1 = static_cast<TH1*>(frame->Clone());
  f1->SetDirectory(0);
  DrawInPad(fBody, 1, f1, "axis logx logy grid");
  DrawInPad(fBody, 1, orig, "nostack same leg2");
  
  
  TH1* f2 = static_cast<TH1*>(frame->Clone());
  f2->SetDirectory(0);
  DrawInPad(fBody, 2, f2, "axis logx logy grid");    
  DrawInPad(fBody, 2, scaled, "nostack same leg2");

  PrintCanvas(Form("%s - \\Delta", fLastBin.Data()));
}
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawBackground(TVirtualPad* c,
					   Int_t        pad,
					   Container*   ress,
					   const char*  name,
					   const char*  pretty)
{
  Container*   sc    = GetC(ress, name);
  TH2*         mask  = GetH2(sc,  "alphaMask");// fiducial
  TH2*         data  = GetH2(sc,  "etaVsIPz");
  TH2*         bgEst = GetH2(sc,  "backgroundEst");
  if (!bgEst)  bgEst = GetH2(sc,  "backgroundEstimate");
  TH2*         mBeta = GetH2(sc,  "oneMinusBeta");
  TH2*         sigEst= GetH2(sc,  "signalEst");
  if (data    && mask) data  ->Multiply(mask);
  if (bgEst   && mask) bgEst ->Multiply(mask);
  if (mBeta   && mask) mBeta ->Multiply(mask);
  if (sigEst  && mask) sigEst->Multiply(mask);

  TVirtualPad* q     = c->cd(pad);
  q->SetRightMargin(0.10);
  q->Divide(4,1,0,0);
  DrawInPad(q,1,data,      "colz");
  DrawInPad(q,2,bgEst,     "colz");
  DrawInPad(q,3,mBeta,     "colz");
  DrawInPad(q,4,sigEst,    "colz");
  TVirtualPad* r = q->GetPad(4);
  r->SetRightMargin(0.12);
  r->SetPad(r->GetXlowNDC(), r->GetYlowNDC(),
	    .99, r->GetYlowNDC()+r->GetHNDC());

  TLatex* ltx = new TLatex(.01,.5, pretty);
  ltx->SetNDC();
  ltx->SetTextFont(62);
  ltx->SetTextAngle(90);
  ltx->SetTextAlign(23);
  ltx->SetTextSize(0.07);
  DrawInPad(q,1,ltx,"");
}
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawBackground(Container* realRess,
					   Container* simRess)
{
  ClearCanvas();

  fBody->SetLeftMargin(0.2);
  fBody->SetTopMargin(0.10);
  fBody->Divide(1,4,0,0);
  DrawBackground(fBody, 1, realRess, "injection", "Injection - Real");
  DrawBackground(fBody, 2, simRess,  "injection", "Injection - Sim.");
  DrawBackground(fBody, 3, simRess,  "combinatorics", "MC Labels");
  DrawBackground(fBody, 4, simRess,  "uncorrelatedCombinatorics",
		 "MC Labels - Uncorr.");

  const char* headings[] = { "\\eta\\hbox{ vs }\\mathrm{IP}_{z}",
			     "\\beta", "1-\\beta",
			     "signal", 0 };
  const char** ptr = headings;
  Int_t i = 0;
  while (*ptr) {
    TLatex* head = new TLatex(fBody->GetRightMargin()+
			      (1-fBody->GetRightMargin())/4*(i+.5),
			      .99, *ptr);
    head->SetNDC();
    head->SetTextAlign(23);
    head->SetTextFont(62);
    head->SetTextSize(0.025);
    DrawInPad(fBody, 0, head, "");
    ptr++;
    i++;
  }
  PrintCanvas(Form("%s\\hbox{ - backgrounds}", fLastBin.Data()));
}
//____________________________________________________________________
TH2* AliTrackletdNdetaPost::CutAlpha(TH2* r)
{
  for (Int_t i = 1; i <= r->GetNbinsX(); i++) {
    for (Int_t j = 1; j <= r->GetNbinsY(); j++) {
      Double_t c = r->GetBinContent(i,j);
      if (c >= fAlphaMin && c <= fAlphaMax) continue;
      r->SetBinContent(i,j,0);
      r->SetBinError  (i,j,0);
    }
  }
  return r;
}
//____________________________________________________________________
TH2* AliTrackletdNdetaPost::GetAlpha(Container* c, const char* name, Bool_t sel)
{
  TH2* a = static_cast<TH2*>(FindSub(c,name,Form("alpha%s",sel ? "Sel" :"")));
  TH2* r = static_cast<TH2*>(a->Clone());
  r->SetDirectory(0);
  CutAlpha(r);
  return r;
}
	
			       
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawAlpha(TVirtualPad* c,
				      Int_t        pad,
				      Container*   ress,
				      const char*  name,
				      const char*  pretty)
{
  Container*   sc    = GetC(ress, name);
  TH2*         mask  = GetH2(sc,  "alphaMask");
  TH2*         al    = GetH2(sc,  "alpha");
  TH2*         als   = GetH2(sc,  "alphaSel");
  if (al  && mask) al ->Multiply(mask);
  if (als && mask) als->Multiply(mask);
    
  TVirtualPad* q = c->cd(pad);
  q->SetRightMargin(0.10);
  q->Divide(2,1,0,0);

  DrawInPad(q,1,al,  "colz");
  DrawInPad(q,2,als, "colz");
  TVirtualPad* r = q->GetPad(2);
  r->SetRightMargin(0.12);
  r->SetPad(r->GetXlowNDC(), r->GetYlowNDC(),
	    .99, r->GetYlowNDC()+r->GetHNDC());

  TLatex* ltx = new TLatex(.01,.5, pretty);
  ltx->SetNDC();
  ltx->SetTextFont(62);
  ltx->SetTextAngle(90);
  ltx->SetTextAlign(23);
  ltx->SetTextSize(0.07);
  DrawInPad(q,1,ltx,"");
}
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawAlpha(Container* realRess, Container* simRess)
{
  ClearCanvas();

  fBody->SetLeftMargin(0.2);
  fBody->SetTopMargin(0.10);
  fBody->Divide(1,3,0,0);
  DrawAlpha(fBody, 1, simRess,  "injection", "Injection - Sim.");
  DrawAlpha(fBody, 2, simRess,  "combinatorics", "MC Labels");
  DrawAlpha(fBody, 3, simRess,  "uncorrelatedCombinatorics",
	    "MC Labels - Uncorr.");

  const char* headings[] = { "\\alpha\\hbox{ all MC}",
			     "\\alpha\\hbox{ selected MC}",
			     0 };
  const char** ptr = headings;
  Int_t i = 0;
  while (*ptr) {
    TLatex* head = new TLatex(fBody->GetRightMargin()+
			      (1-fBody->GetRightMargin())/2*(i+.5),
			      .99, *ptr);
    head->SetNDC();
    head->SetTextAlign(23);
    head->SetTextFont(62);
    head->SetTextSize(0.025);
    DrawInPad(fBody, 0, head, "");
    ptr++;
    i++;
  }
  PrintCanvas(Form("%s - \\alpha", fLastBin.Data()));
}
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawSim(Container* simSums, Container* simRess)
{
  ClearCanvas();

  fBody->Divide(1,2);
  THStack* dNdeta =
    new THStack("dndeta","\\mathrm{d}N_{\\mathrm{ch}}/\\mathrm{d}\\eta");
  TH1* all = GetH1(simRess, "dNdeta");
  TH1* sel = GetH1(simRess, "dNdetaSel");
  all->SetMarkerSize(1.4*all->GetMarkerSize()); // @todo OK
  dNdeta->Add(all);
  dNdeta->Add(sel);
  DrawInPad(fBody, 1, dNdeta, "nostack leg");
  dNdeta->SetMaximum(1.5*dNdeta->GetMaximum("nostack"));

  TVirtualPad* p = fBody->cd(2);
  p->SetRightMargin(0.01);
  p->Divide(1,3,0,0);
  p->GetPad(1)->SetRightMargin(0.01);
  p->GetPad(2)->SetRightMargin(0.01);
  p->GetPad(3)->SetRightMargin(0.01);
  THStack* ipz = new THStack("ipz", "\\mathrm{IP}_{z}");
  ipz->Add(GetH1(simRess, "ipz"));
  ipz->Add(GetH1(simRess, "ipzGen"));
  ipz->Add(GetH1(simRess, "ipzSel"));
  DrawInPad(p, 1, ipz, "nostack grid leg");
  TGraphAsymmErrors* ipzEff =
    static_cast<TGraphAsymmErrors*>(GetO(simRess, "ipzEff",
					 TGraphAsymmErrors::Class()));
  TH2* frame = new TH2D("frame","",100,ipz->GetXaxis()->GetXmin(),
			ipz->GetXaxis()->GetXmax(),
			100, 0.75, 1.1);
  frame->SetDirectory(0);
  frame->SetStats(0);
  frame->SetYTitle("\\epsilon_{\\mathrm{IP}_{z}}");
  DrawInPad(p, 2, frame, "axis");
  DrawInPad(p, 2, ipzEff, "p5 grid");
  DrawInPad(p, 3, GetH2(simSums, "ipzVsGenIPz")->ProfileX(), "grid");

  PrintCanvas(Form("%s\\hbox{ - MC information}",fLastBin.Data()));
}
//____________________________________________________________________
TH2* AliTrackletdNdetaPost::GetBg(Bool_t     comb,
				  Container* cont)
{
  if (comb)    return  GetH2(cont, "etaVsIPz");
  TH2*         bgEst = GetH2(cont, "backgroundEst");
  if (!bgEst)  bgEst = GetH2(cont, "backgroundEstimate");
  return bgEst;
}
//____________________________________________________________________
TH2* AliTrackletdNdetaPost::CorrectSignal(Bool_t comb,
					  TH2*   meas,
					  TH2*   bg,
					  TH2*   simMeas)
{
  if (!comb) {
    // Just return measured minus background 
    meas->Add(bg, -1);
    return meas;
  }
  // Otherwise, we need to calculate 1-beta 
  TH2* beta = static_cast<TH2*>(bg->Clone("beta"));
  beta->Divide(simMeas);
  TH2* one = static_cast<TH2*>(bg->Clone("one"));
  one->Reset();
  for (Int_t i = 1; i <= one->GetNbinsX(); i++) 
    for (Int_t j = 1; j <= one->GetNbinsY(); j++)
      one->SetBinContent(i,j,1);
  one->Add(beta,-1);
  meas->Multiply(one);
  delete one;
  delete beta;    
  return meas;
}				  
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawdNdeta(Container* realList,
				       Bool_t     realComb,
				       Container* simList,
				       Bool_t     simComb,
				       Color_t    color,
				       THStack*   stack)
{
  Double_t   k        = (realList == simList ? 1 : fCombinatoricsScale);
  Container* realData = GetC(realList, "data");
  Container* simData  = GetC(simList,  "data");
  Container* simCombC = GetC(simList, "combinatorics");
  Container* simBgC   = simComb ? simCombC : GetC(simList, "injection");

  
  TH2* realMeas = GetH2(realData, "etaVsIPz");
  TH1* realIPz  = GetH1(realList, "ipz");
  TH2* realBg   = 0;
  if (realComb) {
    // @todo Get directly.
    //
    // If we use MC combinatorics as real background, we simply copy
    // that and scale by constant - note, we need to clone so we do
    // not modify the original
    realBg = static_cast<TH2*>(GetBg(true, simCombC)->Clone("realBg"));
    Printf("Scale real background by %6.4f", k);
    realBg->Scale(k);
  }
  else
    realBg = GetBg(false, realList);
  if (!realMeas || !realBg || !realIPz) {
    Warning("DrawdNdeta", "One or more real data histograms missing: "
	    "meas=%p bg=%p ipz=%p", realMeas, realBg, realIPz);
    return;
  }

  TH2* simMeas  = GetH2(simData,  "etaVsIPz");
  TH1* simIPz   = GetH1(simList,  "ipz");
  TH2* simBg    = GetBg(simComb,  simBgC); // @todo Get directly
  TH2* alpha    = CutAlpha(GetH2(simBgC,   "alpha"));
  if (!simMeas || !simBg || !simIPz || !alpha) {
    Warning("DrawdNdeta", "One or more simuluated data histograms missing: "
	    "meas=%p bg=%p ipz=%p alpha=%p", simMeas, simBg, simIPz, alpha);
    return;
  }

  TH2* trueGen  = GetH2(simList,  "etaVsIPzMC");
  TH1* trueIPz  = GetH1(simList,  "ipzGen");
  if (!trueGen || !trueIPz) {
    Warning("DrawdNdeta", "One or more generator data histograms missing: "
	    "gen=%p ipz=%p", trueGen, trueIPz);
    return;
  }

  realMeas->SetName("realMeas");realMeas->SetTitle("Measured (real)");  
  realBg  ->SetName("realBg");  realBg  ->SetTitle("Background (real)");
  realIPz ->SetName("realIPz"); realIPz ->SetTitle("IP_{z} (real)");
  simMeas ->SetName("simMeas"); simMeas ->SetTitle("Measured (simulated)");
  simBg   ->SetName("simBg");   simBg   ->SetTitle("Background (simulated)");
  simIPz  ->SetName("simIPz");  simIPz  ->SetTitle("IP_{z} (simulated)");
  trueGen ->SetName("trueGen"); trueGen ->SetTitle("Generated");
  trueIPz ->SetName("trueIPz"); trueIPz ->SetTitle("IP_{z} (generated)");
  
  // Calculate real signal as measured minus background
  TH2* realSig = static_cast<TH2*>(realMeas->Clone("realSignal"));
  realSig->SetTitle("Signal (real)");
  realSig->SetDirectory(0);
  // if realComb is true, then realBg points to scaled copy of simBg 
  CorrectSignal(realComb, realSig, realBg, simMeas);

  // Calculate simulated signal as measured minus background - for drawing
  // @todo We should get this from the input once I correct the task 
  TH2* simSig = static_cast<TH2*>(simMeas->Clone("simSig"));
  simSig->SetTitle("Signal (simulated)");
  simSig->SetDirectory(0);
  // CorrectSignal(simComb, simSig, simBg, simMeas);
  simSig->Add(simBg, -1);
  
  // Create result as alpha times real signal 
  TH2* result = static_cast<TH2*>(realSig->Clone("result"));
  result->SetTitle("Result");
  result->SetDirectory(0);
  result->Multiply(alpha);  

  UShort_t mode = 2;// @todo Should be 1 when task is fixed
  // Calculate the primary dN/deta
  TH1* truth  = AverageOverIPz(trueGen, "truth", mode, trueIPz);
  truth->SetYTitle(ObsTitle());
  SetAttr(truth, color, 24, 1.5, 0, 7, 1);
  
  // Calculate the real dN/deta 
  TH1* dndeta = AverageOverIPz(result, "dndeta", mode, realIPz);  
  dndeta->SetYTitle(ObsTitle());
  SetAttr(dndeta, color, 20, 1.2, 0, 1, 1);

  // Get some backgrounds
  TH1* realAvgMeas = AverageOverIPz(realMeas, "realAvgMeas", mode, realIPz);
  TH1* realAvgSig  = AverageOverIPz(realSig,  "realAvgSig",  mode, realIPz);
  TH1* realAvgBg   = AverageOverIPz(realBg,   "realAvgBg",   mode, realIPz);
  TH1* simAvgMeas  = AverageOverIPz(simMeas,  "simAvgMeas",  mode, simIPz);
  TH1* simAvgSig   = AverageOverIPz(simSig,   "simAvgSig",   mode, simIPz);
  TH1* simAvgBg    = AverageOverIPz(simBg,    "simAvgBg",    mode, simIPz);
  SetAttr(realAvgMeas, kGreen+2, 21);
  SetAttr(realAvgSig,  kGreen+2, 22, 1.6);
  SetAttr(realAvgBg,   kGreen+2, 23, 1.6);
  SetAttr(simAvgMeas,  kBlue+2,  21);
  SetAttr(simAvgSig,   kBlue+2,  22, 1.4);
  SetAttr(simAvgBg,    kBlue+2,  23, 1.4);
	   
  
  // Put everything together 
  THStack* summary = new THStack("summary", dndeta->GetYaxis()->GetTitle());
  summary->Add(dndeta, "e2");
  summary->Add(truth,  "e");
  summary->Add(realAvgMeas);
  summary->Add(simAvgMeas);
  summary->Add(realAvgSig);
  summary->Add(simAvgSig);
  summary->Add(realAvgBg);
  summary->Add(simAvgBg);
  
  ClearCanvas();
  DrawInPad(fBody, 0, summary, "nostack leg2 grid");
  if (summary->GetHistogram()->GetMinimum() < 1)
    summary->SetMinimum(1);
  summary->SetMaximum(1.5*summary->GetMaximum("nostack"));

  if (stack) { 
    dndeta->SetTitle(fLastBin);   
    stack->Add(dndeta, "e2");
    stack->Add(truth,  "e");
  }
  
  PrintCanvas(Form("%s - %s", fLastBin.Data(),ObsTitle()));  
}
//____________________________________________________________________
void AliTrackletdNdetaPost::ProcessBin(UInt_t     what,     Int_t      bin,
				       Double_t   c1,       Double_t   c2,
				       Container* realSums, Container* simSums,
				       Container* realRess, Container* simRess,
				       THStack*   stack)
{
  fLastBin.Form("%5.1f\\hbox{%%} - %5.1f\\hbox{%%}", c1, c2);
  Printf("Centrality bin %5.1f%% - %5.1f%%", c1, c2);
  TString name;
  name.Form("cent%03dd%02d_%03dd%02d",
	    Int_t(c1), Int_t(c1*100)%100, 
	    Int_t(c2), Int_t(c2*100)%100);
  Container* realRessBin = GetC(realRess, name);
  Container* simRessBin  = GetC(simRess,  name);
  Container* realSumsBin = GetC(realSums, name);
  Container* simSumsBin  = GetC(simSums,  name);
  if (!realRessBin || !simRessBin) {
    Warning("PostBin", "Missing bin for %5.1f%% - %5.1f%% (%s): %p %p",
	    c1, c2, name.Data(), realRessBin, simRessBin);
    return;
  }
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
  Color_t color = (bin <= 0 ? cc[10] : cc[(bin-1)%11]);

  if (what & kSimInfo)     DrawSim       (simSumsBin,  simRessBin);
  if (what & kSpecies)     DrawSpecies   (simRessBin);
  if (what & kDeltas)      DrawDeltas    (realRessBin, simRessBin);
  if (what & kBackgrounds) DrawBackground(realRessBin, simRessBin);
  if (what & kAlphas)      DrawAlpha     (realRessBin, simRessBin);
  if (what & kdNdetas)     DrawdNdeta    (realRessBin, what & kRealComb,
					  simRessBin,  what & kSimComb,
					  color, stack);
}
//====================================================================
void AliTrackletdNdetaPost::Run(UInt_t      what, 
				UShort_t    maxBins,
				const char* dataName,
				const char* simName)
{
  fLandscape = what & kLandscape;
  fPause     = what & kPause;
  fPDF       = what & kPDF;
    
  TFile* dataFile = 0;
  TFile* simFile  = 0;
  if (!(dataFile = OpenFile(dataName))) return;
  if (!(simFile  = OpenFile(simName)))  return;
    
  Container* realSums = GetC(dataFile, "MiddNdetaSums");
  Container* realRess = GetC(dataFile, "MiddNdetaResults");
  Container* simSums  = GetC(simFile,  "MidMCdNdetaSums");
  Container* simRess  = GetC(simFile,  "MidMCdNdetaResults");
  if (!realSums || !realRess || !simSums || !simRess) return;

  Container* params   = GetC(realSums, "parameters");
  fDeltaCut  = GetD(params, "DeltaCut");
  fTailDelta = GetD(params, "TailDelta");
  fMaxDelta  = GetD(params, "MaxDelta");
  
  CreateCanvas("MiddNdeta.pdf");
  if (what & kGeneral)    DrawGeneral(realSums, simSums, realRess, simRess);
  if (what & kParameters) DrawParams(realSums, simSums);

  TH1* dataCent = GetH1(realSums, "cent");
  TH1* simCent  = GetH1(simSums,  "cent");
  TH1* dataIpz  = GetH1(realSums, "ipz");
  TH1* simIpz   = GetH1(simSums,  "ipz");
    
  if (!CheckConsistency(dataCent, simCent)) {
    Warning("Post", "Centrality bins are incompatible, giving up");
    return;
  }
  if (!CheckConsistency(dataIpz, simIpz)) {
    Warning("Post", "Centrality bins are incompatible, giving up");
    return;
  }
  if (what & kClosure) {
    realSums = simSums;
    realRess = simRess;
  }
  THStack* stack = new THStack("all", ObsTitle());
  // Draw "min-bias" bin
  ProcessBin(what,0,0,100,realSums,simSums,realRess,simRess, stack);

  for (Int_t i = 1; i <= dataCent->GetNbinsX() && i <= maxBins ; i++) {
    Double_t c1 = dataCent->GetXaxis()->GetBinLowEdge(i);
    Double_t c2 = dataCent->GetXaxis()->GetBinUpEdge (i);
      
    ProcessBin(what, i, c1, c2, realSums, simSums, realRess, simRess, stack);
  }

  if (stack->GetHists() && stack->GetHists()->GetEntries() > 0) {
    ClearCanvas();
    TFile* out = TFile::Open("dndeta.root", "RECREATE");
    DrawInPad(fBody, 0, stack, "nostack grid e2");
    if (stack->GetHistogram()->GetMinimum() < 1)
      stack->SetMinimum(1);
    stack->Write();
    out->Write();
    TFile* res = TFile::Open("result.root","READ");
    if (res) {
      THStack* other = static_cast<THStack*>(res->Get("result"));
      if (other)
	DrawInPad(fBody, 0, other, "nostack same");
    }
    PrintCanvas(ObsTitle());
  }
  CloseCanvas();
}

#endif
// Local Variables:
//  mode: C++
// End:

