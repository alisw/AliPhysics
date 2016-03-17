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
#include <TStyle.h>
#include <TBrowser.h>
#include <TGraphAsymmErrors.h>
#include <TProfile.h>
#include <TFitResult.h>
#include <TF1.h>
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
struct AliTrackletdNdetaPost : public AliTrackletdNdetaUtils
{
  /**
   * Processing options 
   * 
   */
  enum {
    /** Use combinatorial background for real data */
    kRealComb     = 0x00001,
    /** Use combinatorial background for simulated data */
    kSimComb      = 0x00002,
    /** Draw dNch/deta */
    kdNdetas      = 0x00004,
    /** Draw general information */
    kGeneral      = 0x00008,
    /** Draw parameters */
    kParameters   = 0x00010,
    /** Draw simulation information */
    kSimInfo      = 0x00020,
    /** Draw species information */
    kSpecies      = 0x00040,
    /** Draw delta information */   
    kDeltas       = 0x00080,
    /** Draw backgrounds */
    kBackgrounds  = 0x00100,
    /** Draw alphas */
    kAlphas       = 0x00200,
    /** MC closure test */
    kClosure      = 0x00400,
    /** Print to PDF */
    kPDF          = 0x01000,
    /** Whether to pause after each plot */
    kPause        = 0x02000,
    /** Draw in landscape */
    kLandscape    = 0x04000,
    /** Fixed scaling */
    kFixedScale   = 0x08000,
    /** Dynamic cent scaling */
    kDynCentScale = 0x10000,
    /** Alternative markers */
    kAltMarker    = 0x20000,
    /** Default options */
    kDefault      = 0x133ff
  };
  enum ECombinatoricsScaleMode {
    kFixed,
    kPerCent,
    kPerEtaCent
  };
  enum {
    kTopBackground = kAzure-8
  };
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
  /** Mode for scaling the combinatorics background */
  ECombinatoricsScaleMode fCombinatoricsScaleMode;
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
   *
   * @return if a legend is drawn, return that
   */
  TLegend* DrawInPad(TVirtualPad* c, Int_t pad, TObject* o, Option_t* opt);
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
  void DrawParams(Container* pars, const char* title, Bool_t comb);
  /** 
   * Draw parameters from both real and simulated analysis 
   * 
   * @param realSums Real data 
   * @param simSums  simulated data 
   */
  void DrawParams(Container* realSums, Bool_t realComb,
		  Container* simSums,  Bool_t simComb);
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
  /** 
   * Draw clusters in sub-pad 
   * 
   * @param pad       Which sub-pad 
   * @param name      Which clusters
   * @param realList  Real data list
   * @param simList   Simulated data list
   */
  void DrawClusters(Int_t pad, const char* name,
		    Container* realList, Container* simList);
  /** 
   * Draw general information 
   * 
   * @param realList Results from real data 
   * @param simList  Results from simulated data 
   */
  void DrawGeneral(Container* realList, Container* simList);
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
   * @param mcList Results from simulated data 
   */
  void DrawSpecies(Container* mcList);
  /** 
   * Draw @f$\Delta@f$ distributions for a single bin 
   * 
   * @param realList 
   * @param simList 
   */
  void DrawDeltas(Container* realList, Container* simList);
  /** 
   * Draw background estimates in a single bin for a single estimate 
   * 
   * @param c     Top pad 
   * @param pad   Sub pad 
   * @param ress  Results 
   * @param name  Name of the results 
   *
   * @return Maximum of histograms
   */
  Double_t DrawBackground(TVirtualPad* c,
			  Int_t        pad,
			  Container*   ress,
			  const char*  name,
			  const char*  pretty);
  /** 
   * Draw all background estimates for a centrality bin 
   * 
   * @param realList Results from real data 
   * @param simList  Results from simulated data 
   */
  void DrawBackground(Container* realList, Container* simList);
  /** 
   * Cut @f$\alpha@f$ edges between fMinAlpha and fMaxAlpha
   * 
   * @param alpha Alpha histogram to modify 
   * 
   * @return @a alpha
   */
  TH2* CutAlpha(TH2* alpha);
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
   * @param realList Results from real data 
   * @param simList  Results from simulated data 
   */
  void DrawAlpha(Container* realList, Container* simList);
  /** 
   * Draw specific stuff from the simulation
   * 
   * @param simList Container of results from simulated data
   */
  void DrawSim(Container* simList);
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
  const char* ObsTitle() const { return "d#it{N}_{ch}/d#eta"; }
  /** 
   * Get the scalar between data and injection 
   * 
   * @param c Top container 
   * @param e Error on scalar 
   * 
   * @return The scalar
   */
  Double_t GetDeltaScale(Container* c, Double_t&  e) const;
  /** 
   * Get marker style to use 
   * 
   * @param what For what 
   * @param sim  For sim?
   * @param alt  Alternate set?
   * 
   * @return Marker style
   */
  Style_t MS(Int_t what, Bool_t sim, Bool_t alt) const;
  /** 
   * Get the scalar distribution between data and injection 
   * 
   * @param c Top container 
   * 
   * @return The scalar histogram
   */
  TH1* GetDeltaScaleDist(Container* c) const;
  /** 
   * Get the integrated delta distribution as a function of eta
   * 
   * @param c Container 
   * 
   * @return Histogram (either created or from the file)
   */
  TH1* GetDeltaIntDist(Container* c) const;
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
   * @param out       (optional) output directory 
   * @param alt       Set alternative marker set 
   */
  void DrawdNdeta(Container*  realList,
		  Bool_t      dataComb,
		  Container*  simList,
		  Bool_t      simComb,
		  Color_t     color=kBlack,
		  THStack*    stack=0,
		  TDirectory* out=0,
		  Bool_t      alt=false);
  //====================================================================
  /** 
   * Process a single bin 
   * 
   * @param what     Bit mask of processing options 
   * @param bin      Centrality bin number 
   * @param c1       Least centrality 
   * @param c2       Largest centrality 
   * @param realList Results from real data 
   * @param simList  Results from simulated data 
   * @param stack    Possible stack to add dN/deta to 
   */
  void ProcessBin(UInt_t     what,     Int_t       bin,
		  Double_t   c1,       Double_t    c2,
		  Container* realList, Container*  simList,
		  THStack*   stack=0,  TDirectory* out=0);
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
    fCombinatoricsScaleMode(kPerEtaCent),
    fAlphaMin(0),
    fAlphaMax(2.5)    
{}
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

//____________________________________________________________________
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
  fBody->SetTicks();
  
  fCanvas->cd();
}
//____________________________________________________________________
void AliTrackletdNdetaPost::CreateCanvas(const TString& outputName)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
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

  if (fPDF) {
    SuppressGuard g;    
    fCanvas->Print(Form("%s[", outputName.Data()),
		   Form("pdf %s", fLandscape ? "Landscape" : ""));
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
void AliTrackletdNdetaPost::CloseCanvas()
{
  if (fPDF && fCanvas) {
    SuppressGuard g;
    fCanvas->Print(Form("%s]", fCanvas->GetTitle()),
		   Form("pdf %s Title:%s",
			fLandscape ? "Landscape" : "",
			fLastTitle.Data()));
    Printf("PDF %s written", fCanvas->GetTitle());
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
    // Suppress prints
    SuppressGuard g;
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
TLegend* AliTrackletdNdetaPost::DrawInPad(TVirtualPad* c,
					  Int_t        pad,
					  TObject*     o,
					  Option_t*    opt)
{
  if (!o) {
    Warning("", "Nothing to draw in pad %d", pad);
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
  if (option.Contains("leg2")) { leg = 2; option.ReplaceAll("leg2",""); }
  if (option.Contains("leg")) {  leg = 1; option.ReplaceAll("leg",""); }
  o->Draw(opt);
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
//====================================================================
void AliTrackletdNdetaPost::DrawParam(const char* name,
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
void AliTrackletdNdetaPost::DrawParams(Container*  pars,
				       const char* title,
				       Bool_t      comb)
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
  Int_t mode = GetI(pars, "RecMode");
  DrawParam("Reconstruction mode",y,
	    Form("%s%s%s",
		 (mode & 0x2 ? "recon" : ""),
		 (mode & 0x4 ? " inj" : ""),
		 (mode & 0x8 ? " rot" : "")));
  DrawParam("Scale by sin^{2}(#it{#theta})",y,
	    bool(GetB(pars, "ScaleDTheta")));
  DrawParam("#delta#phi shift", y, GetD(pars, "DPhiShift"));
  DrawParam("Shifted #delta#phi cut",
	    y, GetD(pars, "ShiftedDPhiCut"));
  Double_t tmp = GetD(pars, "ScaledDThetaCut");
  DrawParam("Scaled #delta#theta cut",
	    y, (tmp < 0 ? "n/a" : Form("%f", tmp)));
  DrawParam("#Delta cut",        y, GetD(pars, "DeltaCut"));
  DrawParam("max#Delta",       y, GetD(pars, "MaxDelta"));
  DrawParam("min#Delta_{tail}", y,
	    GetD(pars,"TailDelta"));
  DrawParam("d#theta window",    y, GetD(pars, "DThetaWindow"));
  DrawParam("d#phi window",      y, GetD(pars, "DPhiWindow"));
  DrawParam("#phi overlap cut",  y, GetD(pars, "PhiOverlapCut"));
  DrawParam("#it{z}-#eta overlap cut",y, GetD(pars, "ZEtaOverlapCut"));
  DrawParam("#phi rotation",     y, TMath::RadToDeg()*
	    GetD(pars, "PhiRotation"));
  y -= 0.045;
  DrawParam("Background method", y, comb ? "Combinatorics" : "Injection");
}
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawParams(Container* realSums, Bool_t realComb,
				       Container* simSums,  Bool_t simComb)
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
  DrawParam("Scaling of comb. bg.", y,
	    (fCombinatoricsScaleMode == kFixed ?
	     Form("%f", fCombinatoricsScale) :
	     fCombinatoricsScaleMode == kPerCent ?
	     "per centrality" :
	     "per eta,centrality"));
  DrawParam("min#alpha", y, fAlphaMin);
  DrawParam("max#alpha", y, fAlphaMax);
  
  // From tasks 
  fBody->cd(2);
  DrawParams(GetC(realSums, "parameters"), "Real data", realComb);
  fBody->cd(3);
  DrawParams(GetC(simSums, "parameters"), "Simulated data", simComb);
  PrintCanvas("Parameters");
}
//====================================================================
THStack* AliTrackletdNdetaPost::Make2Stack(const char*      name,
					   const char*      title,
					   Container*       realList,
					   Container*       simList,
					   Option_t*        dataOpt,
					   Option_t*        simOpt)
{
  TString  nme(name);
  THStack* stack = new THStack(name, title);
  TH1*     data  = CopyH1(realList, name, Form("real%s",name));
  TH1*     sim   = CopyH1(simList,  name, Form("sim%s",name));
  TString dtit(data->GetTitle());
  if (dtit.Contains("\\")) dtit.Form("%s\\hbox{ - real}", data->GetTitle());
  else                     dtit.Form("%s - real", data->GetTitle());
  data->SetTitle(dtit);
  TString stit(sim->GetTitle());
  if (stit.Contains("\\")) stit.Form("%s\\hbox{ - sim.}", sim->GetTitle());
  else                     stit.Form("%s - sim.", sim->GetTitle());
  sim->SetTitle(stit);
  stack->Add(data, dataOpt);
  stack->Add(sim,  simOpt);
  return stack;
}
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawClusters(Int_t pad,
					 const char* name,
					 Container* realList,
					 Container* simList)
{
  fBody->cd(pad);
  TH2* real = GetH2(realList, name);
  TH2* sim  = GetH2(simList,  name);
  real->Draw("colz");
  sim->Draw("box same");
  fBody->GetPad(pad)->Modified();
}

//____________________________________________________________________
void ModLegend(TVirtualPad* p, TLegend* l,
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
void AliTrackletdNdetaPost::DrawGeneral(Container* realList,
					Container* simList)
{
  THStack*   ipz    = Make2Stack("ipz",   "IP_{#it{z}}", realList,simList);
  THStack*   cent   = Make2Stack("cent",  "Centrality [%]",   realList,simList);
  THStack*   status = Make2Stack("status","Task status",      realList,simList,
				 "B text90", "B text90");
  TH1* mcStatus = GetH1(simList, "statusMC");
  mcStatus->SetMaximum(1.2*mcStatus->GetMaximum());
  ClearCanvas();
  fBody->Divide(2,4);
  for (Int_t i = 1; i <= 8; i++) {
    fBody->GetPad(i)->SetRightMargin(0.01);
    fBody->GetPad(i)->SetTopMargin(0.01);
  }
  TLegend* l = 0;
  l = DrawInPad(fBody,1,ipz,     "nostack leg");
  ModLegend(fBody->GetPad(1),l,.4,.1,.75,.4);
  l = DrawInPad(fBody,2,cent,    "nostack leg");
  ModLegend(fBody->GetPad(2),l,.6,.1,.99,.4);
  l = DrawInPad(fBody,3,status,  "nostack hist text90 leg");
  ModLegend(fBody->GetPad(3),l,.5,.7,.99,.99);
  DrawInPad(fBody,4,mcStatus,"hist text90");
  DrawClusters(5, "allClusters0",  realList, simList);
  DrawClusters(6, "allClusters1",  realList, simList);
  DrawClusters(7, "usedClusters0", realList, simList);
  DrawClusters(8, "usedClusters1", realList, simList);
  PrintCanvas("General information");
}
//====================================================================
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
TH1* AliTrackletdNdetaPost::FindDelta(Container*  ress,
				      const char* sub,
				      Bool_t      scaled)
{
  if (!ress) return 0;
  Container* subCont = GetC(ress, sub);
  if (!subCont) return 0;
  TH1* h = GetH1(subCont, Form("delta%s", scaled ? "Scaled" : ""));
  if (!h) return 0;
  if (!scaled) return h;
  TString tit(h->GetTitle());
  if (tit.Contains("\\")) tit.Append(Form("\\hbox{ - %s}",sub));
  else                    tit.Append(Form(" - %s", sub));
  return h;
}
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawSpecies(Container* mcList)
{
  THStack* particlePdgs = new THStack("particlePdg","Particle species");
  particlePdgs->Add(GetH1(mcList, "allPrimaryPdg"));
  particlePdgs->Add(GetH1(mcList, "cutPrimaryPdg"));  
  particlePdgs->Add(GetH1(mcList, "allSecondaryPdg"));
  particlePdgs->Add(GetH1(mcList, "cutSecondaryPdg"));
  particlePdgs->SetMaximum(particlePdgs->GetMaximum("nostack")*50);
  particlePdgs->SetMinimum(8*1e-7);
  THStack* parentPdgs = new THStack("parentPdg","Parent species");
  parentPdgs->Add(GetH1(mcList, "allPrimaryParentPdg"));
  parentPdgs->Add(GetH1(mcList, "cutPrimaryParentPdg"));  
  parentPdgs->Add(GetH1(mcList, "allSecondaryParentPdg"));
  parentPdgs->Add(GetH1(mcList, "cutSecondaryParentPdg"));
  parentPdgs->SetMaximum(parentPdgs->GetMaximum("nostack")*50);
  parentPdgs->SetMinimum(8*1e-5);
  
  ClearCanvas();
  fBody->Divide(1,2,0,0);
  fBody->GetPad(1)->SetRightMargin(0.01);
  fBody->GetPad(2)->SetRightMargin(0.01);

  TLegend* p1Leg = DrawInPad(fBody, 1, particlePdgs, "nostack logy leg2 grid");
  TLegend* p2Leg = DrawInPad(fBody, 2, parentPdgs,   "nostack logy leg2 grid");
  p1Leg->SetHeader(particlePdgs->GetTitle());
  p2Leg->SetHeader(parentPdgs  ->GetTitle());
  ModLegend(fBody->GetPad(1),p1Leg,.2,.8,.99,1);
  ModLegend(fBody->GetPad(2),p2Leg,.2,.8,.99,1);
  
  particlePdgs->GetHistogram()->GetXaxis()->LabelsOption("v");
  parentPdgs  ->GetHistogram()->GetXaxis()->LabelsOption("v");
  // particlePdgs->SetMaximum(1.4*particlePdgs->GetMaximum("nostack"));
  // parentPdgs  ->SetMaximum(1.4*parentPdgs  ->GetMaximum("nostack"));
  
  fBody->GetPad(1)->Modified();
  fBody->GetPad(2)->Modified();

  PrintCanvas(Form("%s - species", fLastBin.Data()));
}
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawDeltas(Container* realList, Container* simList)
{

  ClearCanvas();
  fBody->Divide(1,2,0,0);
  fBody->GetPad(1)->SetRightMargin(0.01);
  fBody->GetPad(2)->SetRightMargin(0.01);
      
  
  THStack* orig = new THStack("orig", "#Delta");
  orig->Add(FindDelta(realList, "data",                     false));
  orig->Add(FindDelta(realList, "injection",                false));
  orig->Add(FindDelta(simList,  "data",                     false));
  orig->Add(FindDelta(simList,  "injection",                false));
  orig->Add(FindDelta(simList,  "primaries",                false));
  orig->Add(FindDelta(simList,  "combinatorics",            false));
  orig->Add(FindDelta(simList,  "secondaries",              false));
  orig->Add(FindDelta(simList,  "uncorrelatedCombinatorics",false));  
  
  THStack* scaled = new THStack("scaled", "#Delta (scaled)");
  scaled->Add(FindDelta(realList, "data",                     false));
  scaled->Add(FindDelta(realList, "injection",                true));
  scaled->Add(FindDelta(simList,  "data",                     false));
  scaled->Add(FindDelta(simList,  "injection",                true));
  scaled->Add(FindDelta(simList,  "primaries",                false));
  scaled->Add(FindDelta(simList,  "combinatorics",            true));
  scaled->Add(FindDelta(simList,  "secondaries",              false));
  scaled->Add(FindDelta(simList,  "uncorrelatedCombinatorics",true));

  Double_t max   = /*5*/1.2*orig->GetMaximum("nostack");
  TH1*     templ = FindDelta(realList, "data", false);
  TH2*     frame = new TH2D("frame","",
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

  TLegend* l = 0;
  TH1* f1 = static_cast<TH1*>(frame->Clone());
  f1->SetDirectory(0);
  DrawInPad(fBody, 1, f1,   "axis logx logy grid");
  l = DrawInPad(fBody, 1, orig, "nostack same leg2");
  // ModLegend(fBody->GetPad(1), l, .2, .8, .99, .99);
  ModLegend(fBody->GetPad(1), l,
	    fBody->GetPad(1)->GetLeftMargin(),
	    fBody->GetPad(1)->GetBottomMargin(),
	    .6, .35);
  
  
  TH1* f2 = static_cast<TH1*>(frame->Clone());
  f2->SetDirectory(0);
  DrawInPad(fBody, 2, f2,     "axis logx logy grid");    
  l = DrawInPad(fBody, 2, scaled, "nostack same leg2");
  // ModLegend(fBody->GetPad(1), l, .2, .8, .99, 1);
  ModLegend(fBody->GetPad(2), l,
	    fBody->GetPad(2)->GetLeftMargin(),
	    fBody->GetPad(2)->GetBottomMargin(),
	    .6, .45);

  PrintCanvas(Form("%s - #Delta", fLastBin.Data()));
}
//____________________________________________________________________
Double_t AliTrackletdNdetaPost::DrawBackground(TVirtualPad* c,
					       Int_t        pad,
					       Container*   ress,
					       const char*  name,
					       const char*  pretty)
{
  SuppressGuard g(kError);
  Container*   sc    = GetC(ress, name);
  TH2*         mask  = GetH2(sc,   "fiducial");
  TH2*         data  = CopyH2(sc,  "etaVsIPz",      "raw");
  TH2*         bgEst = CopyH2(sc,  "backgroundEst", "rawBg");
  TH2*         mBeta = CopyH2(sc,  "beta"); // oneMinusBeta
  TH2*         sigEst= CopyH2(sc,  "signalEst",     "raw");
  if (data    && mask) data  ->Multiply(mask);
  if (bgEst   && mask) bgEst ->Multiply(mask);
  if (mBeta   && mask) mBeta ->Multiply(mask);
  if (sigEst  && mask) sigEst->Multiply(mask);
  mBeta->SetMinimum(0);
  mBeta->SetMaximum(.7);
  Double_t  max = TMath::Max(data->GetMaximum(), sigEst->GetMaximum());
  data  ->SetMinimum(0);
  sigEst->SetMinimum(0);
  bgEst ->SetMinimum(0);
  data  ->SetMaximum(max);
  sigEst->SetMaximum(max);
  
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

  return max;
}
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawBackground(Container* realList,
					   Container* simList)
{
  ClearCanvas();

  fBody->SetLeftMargin(0.2);
  fBody->SetTopMargin(0.10);
  fBody->Divide(1,4,0,0);
  Double_t max = 0;
  max = TMath::Max(max, DrawBackground(fBody, 1, realList,
				       "injection", "Injection - Real"));
  max = TMath::Max(max, DrawBackground(fBody, 2, simList,
				       "injection", "Injection - Sim."));
  max = TMath::Max(max, DrawBackground(fBody, 3, simList,
				       "combinatorics", "MC Labels"));
  max = TMath::Max(max, DrawBackground(fBody, 4, simList,
				       "uncorrelatedCombinatorics",
				       "MC Labels - Uncorr."));
  Int_t i = 0;
  for (i = 1; i <= 4; i++) {
    TVirtualPad* p = fBody->GetPad(i);
    if (!p) continue;
    for (Int_t j = 1; j <= 4; j++) {
      TVirtualPad* q = p->GetPad(j);
      TH1*         h = static_cast<TH1*>(q->FindObject("raw"));
      if (h) {
	h->SetMaximum(max);
	continue;
      }
      h  = static_cast<TH1*>(q->FindObject("rawBg"));
      if (h) h->SetMaximum(max/5);      
    }
  }
  const char* headings[] = { "Measured",
			     "Background est.",
			     "#beta",
			     "Signal", 0 };
  const char** ptr = headings;
  i = 0;
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
  PrintCanvas(Form("%s - backgrounds", fLastBin.Data()));
}
//____________________________________________________________________
void PrintH(TH2* h, Int_t prec=2)
{
  Printf("Content of %s - %s", h->GetName(), h->GetTitle());
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    printf("%3d: ", i);
    for (Int_t j = 1; j <= h->GetNbinsY(); j++) {
      Double_t c = h->GetBinContent(i,j);
      Double_t e = h->GetBinError  (i,j);
      if (TMath::IsNaN(c) || TMath::IsNaN(e))
	printf("*** NAN ***");
      else 
	printf("%.*f+/-%.*f ", prec, c, prec, e);
    }
    printf("\n");
  }
}
//____________________________________________________________________
TH2* AliTrackletdNdetaPost::CutAlpha(TH2* h)
{
  TH2* r = static_cast<TH2*>(h->Clone("fiducial"));
  r->SetDirectory(0);
  r->Reset();
  r->SetTitle("Fiducial cut");
  for (Int_t i = 1; i <= r->GetNbinsX(); i++) {
    for (Int_t j = 1; j <= r->GetNbinsY(); j++) {
      Double_t c = h->GetBinContent(i,j);
      if (c <= fAlphaMin || c > fAlphaMax) continue;
      r->SetBinContent(i,j,1);
      r->SetBinError  (i,j,0);
    }
  }
  // PrintH(h, 1);
  // PrintH(r, 0);
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
  TH2*         mask  = GetH2(sc,  "fiducial");
  TH2*         al    = CopyH2(sc,  "alpha");
  TH2*         als   = CopyH2(sc,  "alphaSel");
  if (al  && mask) al ->Multiply(mask);
  if (als && mask) als->Multiply(mask);
  al ->SetMinimum(fAlphaMin); al ->SetMaximum(fAlphaMax);
  als->SetMinimum(fAlphaMin); als->SetMaximum(fAlphaMax);
    
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
void AliTrackletdNdetaPost::DrawAlpha(Container* realList, Container* simList)
{
  ClearCanvas();

  fBody->SetLeftMargin(0.2);
  fBody->SetTopMargin(0.10);
  fBody->Divide(1,3,0,0);
  DrawAlpha(fBody, 1, simList,  "injection", "Injection - Sim.");
  DrawAlpha(fBody, 2, simList,  "combinatorics", "MC Labels");
  DrawAlpha(fBody, 3, simList,  "uncorrelatedCombinatorics",
	    "MC Labels - Uncorr.");

  const char* headings[] = { "#alpha all MC",
			     "#alpha selected MC",
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
  PrintCanvas(Form("%s - #alpha", fLastBin.Data()));
}
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawSim(Container* simList)
{
  ClearCanvas();

  fBody->Divide(1,2);
  THStack* dNdeta =  new THStack("dndeta",ObsTitle());
  TH1* all = GetH1(simList, "truedNdeta");
  TH1* sel = GetH1(simList, "truedNdetaSel");
  dNdeta->Add(all);
  dNdeta->Add(sel);
  fBody->GetPad(1)->SetTopMargin(0.01);
  fBody->GetPad(1)->SetRightMargin(0.01);
  TLegend* l = DrawInPad(fBody, 1, dNdeta, "nostack grid leg");
  ModLegend(fBody->GetPad(1), l, .4, .1, .99, .4);
  // dNdeta->SetMaximum(1.5*dNdeta->GetMaximum("nostack"));

  TVirtualPad* p = fBody->cd(2);
  p->SetRightMargin(0.0);
  p->SetBottomMargin(0.15);
  p->Divide(1,3,0,0);
  p->GetPad(1)->SetRightMargin(0.01);
  p->GetPad(2)->SetRightMargin(0.01);
  p->GetPad(3)->SetRightMargin(0.01);
  THStack* ipz = new THStack("ipz", "IP_{#it{z}}");
  ipz->Add(GetH1(simList, "ipz"));
  ipz->Add(GetH1(simList, "ipzGen"));
  ipz->Add(GetH1(simList, "ipzSel"));
  l = DrawInPad(p, 1, ipz, "nostack grid leg");
  ipz->GetHistogram()->SetYTitle("Fraction");
  ipz->GetHistogram()->GetYaxis()->SetTitleSize(0.12);
  ipz->GetHistogram()->GetYaxis()->SetTitleOffset(0.4);
  ipz->GetHistogram()->GetYaxis()->SetLabelSize(0.09);
  ModLegend(p->GetPad(1), l, .4, 0, .75, .4);  
  TGraphAsymmErrors* ipzEff =
    static_cast<TGraphAsymmErrors*>(GetO(simList, "ipzEff",
					 TGraphAsymmErrors::Class()));
  DrawInPad(p, 2, ipzEff, "ap5 grid");
  ipzEff->GetHistogram()->SetYTitle("#epsilon_{IP_{#it{z}}}");
  ipzEff->GetHistogram()->GetYaxis()->SetTitleSize(0.12);
  ipzEff->GetHistogram()->GetYaxis()->SetLabelSize(0.09);
  ipzEff->GetHistogram()->GetYaxis()->SetTitleOffset(0.4);
  ipzEff->GetHistogram()->GetXaxis()->SetRangeUser(ipz->GetXaxis()->GetXmin(),
						   ipz->GetXaxis()->GetXmax());
  p->GetPad(2)->Modified();
  TH1* ipzRes = 0;
  DrawInPad(p, 3, ipzRes  = GetH1(simList, "ipzVsGenIPz"), "grid");
  ipzRes->GetYaxis()->SetTitleSize(0.12);
  ipzRes->GetYaxis()->SetTitleOffset(0.4);
  ipzRes->GetYaxis()->SetLabelSize(0.09);
  ipzRes->GetXaxis()->SetTitleSize(0.09);
  ipzRes->GetXaxis()->SetLabelSize(0.09);
  ipzRes->GetXaxis()->SetTitleOffset(0.4);

  PrintCanvas(Form("%s - MC information",fLastBin.Data()));
}
//____________________________________________________________________
Double_t AliTrackletdNdetaPost::GetDeltaScale(Container* c, Double_t&  e) const
{
  Container* data = GetC(c, "data");
  Container* inj  = GetC(c, "injection");
  // Get integral and error from data 
  Double_t dataIntg = GetD(data, "deltaTailIntg");
  Double_t dataInte = GetD(data, "deltaTailIntE");
  // Get integral and error from injection
  Double_t injIntg  = GetD(inj,  "deltaTailIntg");
  Double_t injInte  = GetD(inj,  "deltaTailIntE");
  // Return the ratio of data to injection
  Double_t r = RatioE(dataIntg, dataInte, injIntg, injInte, e);
  // Printf("%20s: data=%9g+/-%9g  inj=%9g+/-%9g  ->  %9g+/-%9g",
  //        c->GetName(), dataIntg, dataInte, injIntg, injInte, r, e);
  return r;
}

//____________________________________________________________________
TH1* AliTrackletdNdetaPost::GetDeltaIntDist(Container* c) const
{
  TH1* deltaIntg = 0;
  {
    // SuppressGuard g(kError);
    deltaIntg = GetH1(c, "deltaTailInt");
  }
  if (deltaIntg) return deltaIntg;
  
  TH2* etaVsDelta = GetH2(c, "etaVsDelta");
  Double_t maxDelta = etaVsDelta->GetYaxis()->GetXmax();
  deltaIntg         = etaVsDelta->ProjectionX("deltaTailInt");
  deltaIntg->SetDirectory(0);
  deltaIntg->Reset();
  deltaIntg->SetYTitle(Form("#int_{%3.1f}^{%4.1f}d#Delta",
			    fTailDelta, maxDelta));
  for (Int_t i = 1; i <= deltaIntg->GetNbinsX(); i++) {
    TH1* tmp = etaVsDelta->ProjectionY("tmp",i,i);
    Double_t eintg;
    Double_t intg = Integrate(tmp, fTailDelta, maxDelta, eintg);
    deltaIntg->SetBinContent(i, intg);
    deltaIntg->SetBinError  (i, eintg);
    delete tmp;
  }
  return deltaIntg;
}
//____________________________________________________________________
TH1* AliTrackletdNdetaPost::GetDeltaScaleDist(Container* c) const
{
  Container* data    = GetC(c, "data");
  Container* inj     = GetC(c, "injection");
  TH1*       dataInt = GetDeltaIntDist(data);
  TH1*       injInt  = GetDeltaIntDist(inj);
  TH1*       ret     = static_cast<TH1*>(dataInt->Clone("deltaTailRatio"));
  ret->SetDirectory(0);
  ret->SetYTitle("Ratio of tails (data/injection)");
  ret->Divide(injInt);
  return ret;
}

//____________________________________________________________________
Style_t AliTrackletdNdetaPost::MS(Int_t what, Bool_t sim, Bool_t alt) const
{
  if (!alt) {
    switch (what) {
    case 0: /* dndeta */ return (sim ? 24 : 20); // cirlce
    case 1: /* M      */ return (sim ? 24 : 20); // circle
    case 2: /* B      */ return (sim ? 25 : 21); // square
    case 3: /* S      */ return (sim ? 26 : 22); // up-triangle
    default:             return (sim ?  5 :  2); //
    }
  }
  switch (what) {
  case 0: /* dndeta */ return (sim ? 25 : 21); // square
  case 1: /* M      */ return (sim ? 28 : 34); // cross
  case 2: /* B      */ return (sim ? 27 : 33); // diamond
  case 3: /* S      */ return (sim ? 30 : 29); // star
  }
  return (sim ?  5 :  2); //
}
//____________________________________________________________________
void AliTrackletdNdetaPost::DrawdNdeta(Container*  realList,
				       Bool_t      realComb,
				       Container*  simList,
				       Bool_t      simComb,
				       Color_t     color,
				       THStack*    stack,
				       TDirectory* out,
				       Bool_t      alt)
{
  Double_t   k        = (realList == simList ? 1 : fCombinatoricsScale);
  Container* realData = GetC(realList, "data");
  Container* simData  = GetC(simList,  "data");
  Container* simCombC = GetC(simList, "combinatorics");
  Container* simBgC   = simComb  ? simCombC : GetC(simList,  "injection");
  Container* realBgC  = realComb ? simCombC : GetC(realList, "injection");

  
  TH2* realMeas = CopyH2(realData, "etaVsIPz",      "realMeas");
  TH1* realIPz  = CopyH1(realList, "ipz",           "realIPz");
  TH2* realBg   = CopyH2(realBgC,  "backgroundEst", "realBg");
  TH2* realSig  = 0;
  Double_t scale  = 0;
  Double_t scaleE = 0;
  TH1*     scaleH = 0;
  if (realComb) {
    // Scale combinatorial background
    if (fCombinatoricsScaleMode == kFixed) {
      // fixed scalar
      Scale(realBg, fCombinatoricsScale, 0);
      scale = k;
    }
    else if (fCombinatoricsScaleMode == kPerCent) {
      // Printf("");
      // Get scalar and error from data
      Double_t eReal = 0;
      Double_t rReal = GetDeltaScale(realList, eReal);
      // Get integral and error from simulations
      Double_t eSim  = 0;
      Double_t rSim  = GetDeltaScale(simList,  eSim);
      // Calculate the scalar of the combinatorial background 
      scale             = RatioE(rReal, eReal, rSim, eSim, scaleE);
      Scale(realBg, scale, scaleE);
    }
    else {
      TH1* realTail = GetDeltaScaleDist(realList);
      TH1* simTail  = GetDeltaScaleDist(simList);
      scaleH        = static_cast<TH1*>(realTail->Clone("scaleBg"));
      scaleH->SetDirectory(0);
      scaleH->SetYTitle("k_{#eta}");
      scaleH->SetTitle("Scalar");
      SetAttr(scaleH, scaleH->GetMarkerColor(), 34);
      scaleH->Divide(simTail);
      TF1* scaleF = new TF1("ff","pol0",
			    scaleH->GetXaxis()->GetXmin(),
			    scaleH->GetXaxis()->GetXmax());
      TFitResultPtr scaleR = scaleH->Fit(scaleF, "QN0SR");
      scale  = scaleR->Parameter(0);
      scaleE = scaleR->ParError(0);
      Scale(realBg, scaleH);
    }
    // Get the background beta, and scale it
    TH2* beta = CopyH2(realBgC, "beta");
    Scale(beta, scale, scaleE);
    // Create unit histogram
    TH2* one  = CopyH2(realBgC, "beta", "one");
    one->Reset();
    for (Int_t i = 1; i <= one->GetNbinsX(); i++) 
      for (Int_t j = 1; j <= one->GetNbinsY(); j++)
	one->SetBinContent(i,j,1);
    // Subtract k times beta    
    one->Add(beta,-1);
    // Multiply on to the measured distribution
    realSig = CopyH2(realData, "etaVsIPz", "realSig");
    realSig->SetMarkerStyle(beta->GetMarkerStyle());
    realSig->Multiply(one);
    // Clean-up temporary histograms
    delete beta;
    delete one;
  }
  else {
    // We can get the signal estimate directory from the background
    // container.
    realSig = CopyH2(realBgC, "signalEst");
  }
  if (!realMeas || !realSig || !realBg || !realIPz) {
    Warning("DrawdNdeta", "One or more real data histograms missing: "
	    "meas=%p sig=%p bg=%p ipz=%p", realMeas, realSig, realBg, realIPz);
    return;
  }
  if (realBg ->GetMarkerStyle() == 30) realBg ->SetMarkerStyle(29);
  if (realBg ->GetMarkerStyle() == 27) realBg ->SetMarkerStyle(33);
  if (realSig->GetMarkerStyle() == 30) realSig->SetMarkerStyle(29);
  if (realSig->GetMarkerStyle() == 27) realSig->SetMarkerStyle(33);

  TH2* simMeas  = CopyH2(simData,  "etaVsIPz",      "simMeas");
  TH1* simIPz   = CopyH1(simList,  "ipz",           "simIPz");
  TH2* simBg    = CopyH2(simBgC,   "backgroundEst", "simBg");
  TH2* simSig   = CopyH2(simBgC,   "signalEst",     "simSig");
  TH2* alpha    = CopyH2(simBgC,   "alpha");
  TH2* fiducial = CutAlpha(alpha);
  if (!simMeas || !simSig || !simBg || !simIPz || !alpha) {
    Warning("DrawdNdeta", "One or more simuluated data histograms missing: "
	    "meas=%p sig=%p bg=%p ipz=%p alpha=%p",
	    simMeas, simSig, simBg, simIPz, alpha);
    return;
  }
  alpha->Multiply(fiducial);
  TH2* trueGen  = CopyH2(simList,  "etaVsIPzMC", "trueGen");
  TH1* trueIPz  = CopyH1(simList,  "ipzGen",     "trueIPz");
  if (!trueGen || !trueIPz) {
    Warning("DrawdNdeta", "One or more generator data histograms missing: "
	    "gen=%p ipz=%p", trueGen, trueIPz);
    return;
  }

  realMeas->SetTitle("Measured (real)");  
  realBg  ->SetTitle("Background (real)");
  realSig ->SetTitle("Signal (real)");
  realIPz ->SetTitle("IP_{z} (real)");
  simMeas ->SetTitle("Measured (simulated)");
  simBg   ->SetTitle("Background (simulated)");
  simSig  ->SetTitle("Signal (simulated)");
  simIPz  ->SetTitle("IP_{z} (simulated)");
  trueGen ->SetTitle("Generated");
  trueIPz ->SetTitle("IP_{z} (generated)");
  
  // Create result as alpha times real signal 
  TH2* result = static_cast<TH2*>(realSig->Clone("result"));
  result->SetTitle("Result");
  result->SetDirectory(0);
  result->Multiply(alpha);  

  UShort_t mode = 1;
  // Calculate the primary dN/deta
  TH1* truth  = AverageOverIPz(trueGen, "truth", mode, trueIPz, 0);
  truth->SetYTitle(ObsTitle());
  SetAttr(truth, color, MS(0,true,alt), 1.5, 0, 7, 1);
  
  // Calculate the real dN/deta 
  TH1* dndeta = AverageOverIPz(result, "dndeta", mode, realIPz, 0);  
  dndeta->SetYTitle(ObsTitle());
  SetAttr(dndeta, color, MS(0,false,alt), 1.2, 0, 1, 1);

  // Multiply histograms by fiducial cut
  realMeas->Multiply(fiducial);
  simMeas ->Multiply(fiducial);
  realSig ->Multiply(fiducial);
  simSig  ->Multiply(fiducial); 
  realBg  ->Multiply(fiducial);
  simBg   ->Multiply(fiducial);
  realSig ->Multiply(fiducial);
  simSig  ->Multiply(fiducial);
  // Get some step histograms.  Note, we use the result histogram to
  // select which bins we project so that we get the used sub-sample.
  TH1* realAvgMeas = AverageOverIPz(realMeas,"realAvgMeas",mode,realIPz,result);
  TH1* realAvgSig  = AverageOverIPz(realSig, "realAvgSig", mode,realIPz,result);
  TH1* realAvgBg   = AverageOverIPz(realBg,  "realAvgBg",  mode,realIPz,result);
  TH1* simAvgMeas  = AverageOverIPz(simMeas, "simAvgMeas", mode,simIPz, result);
  TH1* simAvgSig   = AverageOverIPz(simSig,  "simAvgSig",  mode,simIPz, result);
  TH1* simAvgBg    = AverageOverIPz(simBg,   "simAvgBg",   mode,simIPz, result);
  SetAttr(realAvgMeas,realAvgMeas->GetMarkerColor(), MS(1,false,alt));
  SetAttr(realAvgSig, realAvgSig ->GetMarkerColor(), MS(2,false,alt));
  SetAttr(realAvgBg,  realAvgSig ->GetMarkerColor(), MS(3,false,alt));
  SetAttr(simAvgMeas, simAvgSig  ->GetMarkerColor(), MS(1,true, alt));
  SetAttr(simAvgSig,  simAvgSig  ->GetMarkerColor(), MS(2,true, alt));
  SetAttr(simAvgBg,   simAvgBg   ->GetMarkerColor(), MS(3,true, alt));
  realAvgMeas->SetYTitle("#LT M#GT");
  realAvgSig ->SetYTitle("#LT S#GT");
  realAvgBg  ->SetYTitle("#LT B#GT");
  simAvgMeas ->SetYTitle("#LT M'#GT");
  simAvgSig  ->SetYTitle("#LT S'#GT");
  simAvgBg   ->SetYTitle("#LT B'#GT");
  
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
  if (scaleH) {
    TH1* tmp = static_cast<TH1*>(scaleH->Clone("scaleBg"));
    Double_t sf = TMath::Power(10,Int_t(TMath::Log10(realBg->GetMaximum())+.3));
    tmp->SetTitle(Form("%s#times%3.0f",scaleH->GetTitle(),sf));
    tmp->Scale(sf);
    tmp->SetDirectory(0);
    summary->Add(tmp);
  }
  
  ClearCanvas();
  fBody->SetLeftMargin(0.15);
  TLegend* l = DrawInPad(fBody, 0, summary, "nostack leg2 grid");
  summary->GetHistogram()->SetXTitle("#eta");
  summary->GetHistogram()->SetYTitle(ObsTitle());
  summary->GetHistogram()->GetYaxis()->SetTitleOffset(1.7);
  fBody->Modified();
  fBody->Update();
  fBody->cd();
  l->SetName("legend");
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  ModLegend(fBody, l, fBody->GetLeftMargin(), .7, .99, .99);
  Double_t lim = .5;
  TF1* f = new TF1("ff", "pol0", -lim, lim);
  TFitResultPtr r = dndeta->Fit(f, "N0QSR", "", -lim, lim);
  // r->Print();
  TLatex* ltx = new TLatex(fBody->GetLeftMargin()+
			   (1-fBody->GetLeftMargin()-fBody->GetRightMargin())/2,
			   .3,
			   Form("%s|_{|#eta|<%3.1f}"
				"=%6.1f#pm%6.1f (#chi^{2}/#nu=%5.2f)",
				ObsTitle(), lim, r->Parameter(0),
				r->ParError(0), r->Chi2()/r->Ndf()));
  ltx->SetTextAlign(22);
  ltx->SetNDC();
  ltx->SetTextFont(42);
  ltx->SetTextSize(0.035);
  ltx->Draw();
  printf("%6.1f +/- %6.1f (%5.2f - %5.3f+/-%5.3f)",	 
	 r->Parameter(0), r->ParError(0),
	 r->Chi2()/r->Ndf(), scale, scaleE);
  if (scale > 0)
    ltx->DrawLatex(fBody->GetLeftMargin()+
		   (1-fBody->GetLeftMargin()-fBody->GetRightMargin())/2,
		   .25,
		   Form("%sk%s=%4.2f #pm %4.2f",
			scaleH ? "#LT" : "",
			scaleH ? "#GT"  : "",
			scale, scaleE));
  summary->GetHistogram()->GetListOfFunctions()->Add(l);
  summary->GetHistogram()->GetListOfFunctions()->Add(ltx);
  
  fBody->Modified();
  if (summary->GetHistogram()->GetMinimum() < 1)
    summary->SetMinimum(1);
  summary->SetMaximum(1.5*summary->GetMaximum("nostack"));

  if (stack) { 
    dndeta->SetTitle(fLastBin);   
    stack->Add(dndeta, "e2");
    stack->Add(truth,  "e");
  }

  if (out) {
    out->cd();
    result  ->Write();
    dndeta  ->Write();
    summary ->Write();
    if (scaleH) scaleH->Write();
    (new TParameter<double>("scale",  scale)) ->Write();
    (new TParameter<double>("scaleE", scaleE))->Write();
    
    TDirectory* details = out->mkdir("details");
    details->cd();
    realMeas->Write();
    realBg  ->Write();
    realSig ->Write();
    simMeas ->Write();
    simBg   ->Write();
    simSig  ->Write();
    trueGen ->Write();
    alpha   ->Write();
    fiducial->Write();
    GetH1(realData, "delta")                   ->Write("realDataDelta");
    GetH1(GetC(realList, "injection"), "delta")->Write("realInjDelta");
    GetH1(simData, "delta")                    ->Write("simDataDelta");
    GetH1(GetC(simList, "injection"), "delta") ->Write("simInjDelta");
    
    TDirectory* averages = out->mkdir("averages");
    averages->cd();
    realAvgMeas ->Write(); 
    realAvgSig  ->Write(); 
    realAvgBg   ->Write(); 
    simAvgMeas  ->Write(); 
    simAvgSig   ->Write(); 
    simAvgBg    ->Write(); 
    truth       ->Write();
    realIPz     ->Write();
    simIPz      ->Write();
    trueIPz     ->Write();
  }
  PrintCanvas(Form("%s - %s", fLastBin.Data(),ObsTitle()));  
}
//____________________________________________________________________
void AliTrackletdNdetaPost::ProcessBin(UInt_t     what,     Int_t       bin,
				       Double_t   c1,       Double_t    c2,
				       Container* realList, Container*  simList,
				       THStack*   stack,    TDirectory* out)
{
  fLastBin.Form("%5.1f%% - %5.1f%%", c1, c2);
  printf("Centrality bin %5.1f%% - %5.1f%%: ", c1, c2);
  TString name;
  name.Form("cent%03dd%02d_%03dd%02d",
	    Int_t(c1), Int_t(c1*100)%100, 
	    Int_t(c2), Int_t(c2*100)%100);
  Container* realListBin = GetC(realList, name);
  Container* simListBin  = GetC(simList,  name);
  if (!realListBin || !simListBin) {
    Warning("PostBin", "Missing bin for %5.1f%% - %5.1f%% (%s): %p %p",
	    c1, c2, name.Data(), realListBin, simListBin);
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
  Color_t     color  = (bin <= 0 ? cc[10] : cc[(bin-1)%11]);
  TDirectory* binOut = 0;
  Bool_t      alt    = (what & kAltMarker);
  if (out) binOut = out->mkdir(name);
    
  if (what & kSimInfo)     DrawSim       (simListBin);
  if (what & kSpecies)     DrawSpecies   (simListBin);
  if (what & kDeltas)      DrawDeltas    (realListBin, simListBin);
  if (what & kBackgrounds) DrawBackground(realListBin, simListBin);
  if (what & kAlphas)      DrawAlpha     (realListBin, simListBin);
  if (what & kdNdetas)     DrawdNdeta    (realListBin, what & kRealComb,
					  simListBin,  what & kSimComb,
					  color, stack, binOut, alt);
  printf(" done \n");
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
  if      (what & kFixedScale)   fCombinatoricsScaleMode = kFixed;
  else if (what & kDynCentScale) fCombinatoricsScaleMode = kPerCent;
  else                           fCombinatoricsScaleMode = kPerEtaCent;
    
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
  
  CreateCanvas(Form("MiddNdeta_0x%04x.pdf", what));
  if (what & kGeneral)    DrawGeneral(realRess, simRess);
  if (what & kParameters) DrawParams(realSums, what & kRealComb,
				     simSums,  what & kSimComb);

  TH1* realCent = CopyH1(realSums, "cent", "realCent");
  TH1* simCent  = CopyH1(simSums,  "cent", "simCent");
  TH1* realIPz  = GetH1(realSums, "ipz");
  TH1* simIPz   = GetH1(simSums,  "ipz");
    
  if (!CheckConsistency(realCent, simCent)) {
    Warning("Post", "Centrality bins are incompatible, giving up");
    return;
  }
  if (!CheckConsistency(realIPz, simIPz)) {
    Warning("Post", "Centrality bins are incompatible, giving up");
    return;
  }
  if (what & kClosure) realRess = simRess;

  THStack* stack = new THStack("all", ObsTitle());
  // Draw "min-bias" bin
  ProcessBin(what, 0, 0, 100, realRess, simRess, stack);

  TFile* out = TFile::Open(Form("MiddNdeta_0x%1x.root", 0x3 & what),
			   "RECREATE");
  realCent->Write();
  simCent ->Write();
  for (Int_t i = 1; i <= realCent->GetNbinsX() && i <= maxBins ; i++) {
    Double_t c1 = realCent->GetXaxis()->GetBinLowEdge(i);
    Double_t c2 = realCent->GetXaxis()->GetBinUpEdge (i);
      
    ProcessBin(what, i, c1, c2, realRess, simRess, stack, out);
  }

  if (stack->GetHists() && stack->GetHists()->GetEntries() > 0) {
    ClearCanvas();
    DrawInPad(fBody, 0, stack, "nostack grid e2");
    if (stack->GetHistogram()->GetMinimum() < 1)
      stack->SetMinimum(1);
    out->cd();
    stack->Write();
    TFile* res = 0;
    {
      SuppressGuard g(kFatal);
      res = TFile::Open("result.root","READ");
    }
    if (res) {
      THStack* other = static_cast<THStack*>(res->Get("result"));
      if (other)
	DrawInPad(fBody, 0, other, "nostack same");
      out->cd();
      other->Write();
    }
    PrintCanvas(ObsTitle());
  }
  Printf("Results stored in %s", out->GetName());
  out->Write();
  CloseCanvas();
}

#endif
// Local Variables:
//  mode: C++
// End:

