/**
 * @file   DrawdNdeta.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
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
#include <TH2.h>
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
#include <TArrow.h>
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
 * @deprecated Use new GSE based drawing 
 * @ingroup pwglf_forward_dndeta
 */
struct dNdetaDrawer 
{
  enum EFlags { 
    kShowRatios     = 0x00001, 
    kShowLeftRight  = 0x00002, 
    kShowSysError   = 0x00004, 
    kShowRings      = 0x00008,
    kCutEdges       = 0x00010,
    kRemoveOuters   = 0x00020, 
    kUseFinalMC     = 0x00040,
    kUseEmpirical   = 0x00080,
    kForceMB        = 0x00100,
    kMirror         = 0x00200,
    kExport         = 0x00400, 
    kAddExec        = 0x00800,
    kOldFormat      = 0x01000,
    kVerbose        = 0x02000,
    kHiRes          = 0x04000,
    kExtraWhite     = 0x08000,
    kLogo           = 0x10000,
    kNoCentral      = 0x20000,
    kNoLabels       = 0x40000,
    kDefaultOptions = 0x1CE07
  };
  enum EOutFormat { 
    kPNG        = 0x1, 
    kPDF        = 0x2, 
    kROOT       = 0x4, 
    kScript     = 0x8,
    kAllFormats = 0xF
  };
  struct MarkerUtil 
  {
    /**
     * Marker styles 
     */
    enum { 
      kSolid        = 0x000, 
      kHollow       = 0x001, 
      kCircle       = 0x002,
      kSquare       = 0x004, 
      kUpTriangle   = 0x006, 
      kDownTriangle = 0x008, 
      kDiamond      = 0x00a,
      kCross        = 0x00c,
      kStar         = 0x00e
    };
    /** 
     * Get the marker style from option bits
     * 
     * @param bits Option bits 
     * 
     * @return Marker style 
     */
    static Int_t GetMarkerStyle(UShort_t bits)
    {
      Int_t  base   = bits & (0xFE);
      Bool_t hollow = bits & kHollow;
      switch (base) { 
      case kCircle:       return (hollow ? 24 : 20);
      case kSquare:       return (hollow ? 25 : 21);
      case kUpTriangle:   return (hollow ? 26 : 22);
      case kDownTriangle: return (hollow ? 32 : 23);
      case kDiamond:      return (hollow ? 27 : 33); 
      case kCross:        return (hollow ? 28 : 34); 
      case kStar:         return (hollow ? 30 : 29); 
      }
      return 1;
    }
    /** 
     * Get the marker option bits from a style 
     * 
     * @param style Style
     * 
     * @return option bits
     */
    static UShort_t GetMarkerBits(Int_t style)
    { 
      UShort_t bits = 0;
      switch (style) { 
      case 24: case 25: case 26: case 27: case 28: case 30: case 32: 
	bits |= kHollow; break;
      }
      switch (style) { 
      case 20: case 24: bits |= kCircle;       break;
      case 21: case 25: bits |= kSquare;       break;
      case 22: case 26: bits |= kUpTriangle;   break;
      case 23: case 32: bits |= kDownTriangle; break;
      case 27: case 33: bits |= kDiamond;      break;
      case 28: case 34: bits |= kCross;        break;
      case 29: case 30: bits |= kStar;         break;
      }
      return bits;
    }
    /** 
     * Flip an option bit 
     * 
     * @param style Style parameter
     * 
     * @return New style 
     */
    static Int_t FlipHollowStyle(Int_t style)
    {
      UShort_t bits = GetMarkerBits(style);
      Int_t    ret  = GetMarkerStyle(bits ^ kHollow);
      return ret;
    }    
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
      fFwdSysErr(0.07),      // Systematic error in forward range
      fCenSysErr(0.02),      // Systematic error in central range 
      fTitle(""),            // Title on plot
      fBase(""),             // Optional base name of output files
      fClusterScale(""),     // Scaling of clusters to tracklets      
      fFinalMC(""),          // Final MC correction file name
      fEmpirical(""),        // Empirical correction file name
      fDelta(0.),            // IP_delta wrt to (0.005,0.184)
      // Read (or set) information 
      fTrigString(0),        // Trigger string (read, or set)
      fNormString(0),        // Normalisation string (read, or set)
      fSNNString(0),         // Energy string (read, or set)
      fSysString(0),         // Collision system string (read or set)
      fVtxAxis(0),           // Vertex cuts (read or set)
      fCentAxis(0),          // Centrality axis
      fCentMeth(0),          // Centrality axis
      fTriggerEff(1),        // Trigger efficency 
      fExtTriggerEff(false), // True if fTriggerEff was read 
      fCentMin(0),           // Least centrality to plot
      fCentMax(100),         // Largest centrality to plot
      fCentSeen(0x0),        // 32bits of centtraliy seen flags
      // Resulting plots 
      fResults(0),           // Stack of results 
      fRatios(0),            // Stack of ratios 
      fLeftRight(0),         // Left-right asymmetry
      fOthers(0),            // Older data 
      fTriggers(0),          // Number of triggers
      fTruth(0),             // Pointer to truth
      // Other stuff 
      fRangeParam(0),        // Parameter object for range zoom 
      fEmpCorr(0)
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
    if (fCentMeth)   { delete fCentMeth;   fCentMeth   = 0; }
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
  /** 
   * Set shift of IP in (x,y) wrt to reference (x,y)=(0.005,0.184)
   * position. Thie shift must be given in milimeters.
   * 
   * @param delta Shift wrt to refernce in milimeter
   */
  void SetDelta(Double_t delta) { fDelta = delta; }
  void SetDelta(Double_t meanIpx, Double_t meanIpy)
  {
    Info("","Setting delta from %f,%f",meanIpx,meanIpy);
    if (meanIpx < 0 && meanIpy < 0) return;
    const Double_t refX = -0.004;
    const Double_t refY = 0.184;
    Double_t       dx   = (meanIpx - refX);
    Double_t       dy   = (meanIpy - refY);
    Info("","Shifts (%f-%f)=%f, (%f-%f)=%f",
	 meanIpx, refX, dx,
	 meanIpy, refY, dy);
    if (TMath::Abs(dx) < 1e-3 && TMath::Abs(dy) < 1e-3) return;
    Double_t       delta = TMath::Sqrt(dx*dx+dy*dy);
    fDelta = delta;
  }
		  
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
				    sys == 4 ? "Pbp" : 
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
    fTrigString = new TNamed("trigString", (trig & 0x1 ? "MBOR" : 
					    trig & 0x2 ? "INEL>0" : 
					    trig & 0x4 ? "MBAND5" :
					    trig & 0x2000 ? "V0-AND" :
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
  /** 
   * Try to get empirical correction from a path or prefix path 
   * 
   * @param prx     Path or prefix 
   * @param fwdEmp  On return, pounter to the objects or null 
   * @param cenEmp  On return, pounter to the objects or null 
   * @param empName On return, the full path
   * 
   * @return true on success 
   */
  Bool_t GetEmpirical(const TString& prx,
		      Bool_t         useCen,
		      TObject*&      fwdEmp,
		      TObject*&      cenEmp,
		      TString&       empName)   
  {
    empName = "";
    TString path(prx);
    if (!path.Contains("empirical.root")) {
      path = gSystem->ConcatFileName(gSystem->ExpandPathName(prx.Data()),
				     "empirical.root");
      path.Append("#default");
    }
    TUrl   empUrl(path);
    TFile* empirical = TFile::Open(empUrl.GetUrl(), "READ");
    if (!empirical) return false;

    const char* empPath = empUrl.GetAnchor();
    TObject*    fwdObj  = empirical->Get(Form("Forward/%s", empPath));
    TObject*    cenObj  = empirical->Get(Form("Central/%s", empPath));
    if (!(fwdObj &&
	  (fwdObj->IsA()->InheritsFrom(TH1::Class()) || 
	   fwdObj->IsA()->InheritsFrom(TGraphAsymmErrors::Class())))) { 
      Warning("Run", "Didn't get the object Forward/%s from %s", 
	      empPath, empUrl.GetUrl());
    }
    if (useCen &&
	!(cenObj &&
	  (cenObj->IsA()->InheritsFrom(TH1::Class()) || 
	   cenObj->IsA()->InheritsFrom(TGraphAsymmErrors::Class())))) { 
      Warning("Run", "Didn't get the object Central/%s from %s", 
	      empPath, empUrl.GetUrl());
    }
    else {
      fwdEmp  = fwdObj;
      cenEmp  = fwdObj;
      empName = empUrl.GetUrl();
      if (fwdEmp->IsA()->InheritsFrom(TH1::Class())) {
	TH1* h = static_cast<TH1*>(fwdEmp);
	h->SetDirectory(0);
	if (fDelta) {
#if 0
	  TF1* f = new TF1("corr",
	  		   "1+[0]"
	  		   "+(x<0)*sqrt(2)*[0]"
	  		   "+(x<[1])*[2]*[0]*TMath::Power(x-[1],2)",
	  		   -6, 6);
#else 
	  TF1* f = new TF1("corr", "1+[2]*([0]+(x<[1])*pow([0]*(x-[1]),2))");
#endif
	  f->SetParNames("delta","eta0","a");
	  f->SetParameter(0,fDelta);
	  f->SetParameter(1,-2.0);
	  f->SetParameter(2,.10); //TMath::Sqrt(2)); // 0.5);
	  // f->Print();
	  Info("", "Applying correction for IP_delta=%f", fDelta);
	  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
	    Double_t c   = h->GetBinContent(i);
	    if (c < 1e-6) continue;
	    
	    Double_t e   = h->GetBinError(i);
	    Double_t eta = h->GetXaxis()->GetBinCenter(i);
	    Double_t cor = f->Eval(eta);
	    // Info("", "%5.2f -> %7.4f", eta, cor);
	    h->SetBinContent(i, c*cor);
	    h->SetBinError(i, e*cor);
	  }
	}
      }
      if (cenEmp->IsA()->InheritsFrom(TH1::Class()))
	static_cast<TH1*>(cenEmp)->SetDirectory(0);
    }
    empirical->Close();
    return !empName.IsNull();    
  }
  
  //==================================================================
  /** 
   * @{ 
   * @name Main steering functions 
   */
  void Run(const char* filename="forward_dndeta.root", 
	   const char* title="", 
	   const char* others="all",
	   const char* options="default",
	   const char* formats="all",
	   UShort_t    rebin=5,
	   Float_t     eff=0,
	   UShort_t    centMin=0, 
	   UShort_t    centMax=0, 
	   Float_t     vzMin=+999,
	   Float_t     vzMax=-999, 
	   const char* base="")
  {
    Warning("Run","\n"
	    "============================================================\n"
	    "\n"
	    "This script is deprecated.  Please use new GSE based drawing\n"
	    "\n"
	    "============================================================\n");
    TString  ostr(others); ostr.ToUpper();
    UShort_t obits = 0x0;
    if (ostr.EqualTo("ALL")) obits = 0xf;
    else { 
      if (ostr.Contains("UA5"))   obits |= 0x1;
      if (ostr.Contains("CMS"))   obits |= 0x2;
      if (ostr.Contains("ALICE")) obits |= 0x4;
      if (ostr.Contains("WIP"))   obits |= 0x8;
    }
   
    TString fstr(options);
    UInt_t fbits  = 0;
    if (fstr.EqualTo("default", TString::kIgnoreCase)) fbits = kDefaultOptions;
    else {
      TObjArray* farr = fstr.Tokenize(" ,");
      TIter next(farr);
      TObjString* ftoken = 0;
      while ((ftoken = static_cast<TObjString*>(next()))) {
	TString& token = ftoken->String();
	token.ToLower();
	if      (token.BeginsWith("ratio"))         fbits |= kShowRatios;
        else if (token.BeginsWith("asym"))          fbits |= kShowLeftRight;
        else if (token.BeginsWith("left"))          fbits |= kShowLeftRight;
        else if (token.BeginsWith("syse"))          fbits |= kShowSysError;
        else if (token.BeginsWith("rings"))         fbits |= kShowRings;
        else if (token.BeginsWith("noedge"))        fbits |= kCutEdges;
        else if (token.BeginsWith("noout"))         fbits |= kRemoveOuters;
        else if (token.BeginsWith("finalmc"))       fbits |= kUseFinalMC;
        else if (token.BeginsWith("mb"))            fbits |= kForceMB;
        else if (token.BeginsWith("mirror"))        fbits |= kMirror;
        else if (token.BeginsWith("export"))        fbits |= kExport;
        else if (token.BeginsWith("exec"))          fbits |= kAddExec;
        else if (token.BeginsWith("old"))           fbits |= kOldFormat;
        else if (token.BeginsWith("verbose"))       fbits |= kVerbose;
        else if (token.BeginsWith("hires"))         fbits |= kHiRes;
        else if (token.BeginsWith("extraw"))        fbits |= kExtraWhite;
        else if (token.BeginsWith("logo"))          fbits |= kLogo;
        else if (token.BeginsWith("nocentr"))       fbits |= kNoCentral;
        else if (token.BeginsWith("nolabels"))      fbits |= kNoLabels;
        else if (token.BeginsWith("empirical"))  {
	  fbits |= kUseEmpirical;
	  TObjArray* parts=token.Tokenize("=");
	  if (parts->GetEntriesFast() > 1) 
	    SetEmpirical(parts->At(1)->GetName());
	  delete parts;
	}
      }  
      delete farr;
    }
    TString estr(formats); estr.ToUpper();
    UShort_t ebits = 0x0;
    if (ostr.EqualTo("ALL")) ebits = kAllFormats;
    else { 
      if (ostr.Contains("PNG"))  ebits |= kPNG;
      if (ostr.Contains("PDF"))  ebits |= kPDF;
      if (ostr.Contains("ROOT")) ebits |= kROOT;
      if (ostr.Contains("C"))    ebits |= kScript;
    }

    Run(filename, title, rebin, obits, fbits, 0, 0, 0, eff, 
	centMin, centMax, vzMin, vzMax, base, ebits);
  }
  /** 
   * Run the job
   * 
   * @param filename  Input file name  
   * @param title     Title to put on plot
   * @param rebin     Rebinning factor 
   * @param others    Which other to draw 
   * @param flags     Flags for the drawing 
   * @param sys       (optional) Collision system
   * @param sNN       (optional) Collision energy 
   * @param trg       (optional) Trigger 
   * @param eff       (optional) Efficiency 
   * @param centMin   (optional) Least centrality 
   * @param centMax   (optional) Largest centrality 
   * @param vzMin     (optional) Least @f$ IP_z@f$ 
   * @param vzMax     (optional) Largest @f$ IP_z@f$
   * @param base      Basename for output files
   * @param formats   Formats to export to 
   */
  void Run(const char* filename, 
	   const char* title, 
	   UShort_t    rebin, 
	   UShort_t    others=0x7, 
	   UInt_t      flags=kDefaultOptions,
	   UShort_t    sys=0,
	   UShort_t    sNN=0, 
	   UShort_t    trg=0, 
	   Float_t     eff=0, 
	   UShort_t    centMin=0, 
	   UShort_t    centMax=0, 
	   Float_t     vzMin=+999, 
	   Float_t     vzMax=-999,
	   const char* base="", 
	   UShort_t    formats=kAllFormats)
  {
    SetRebin(rebin);
    SetTitle(title);
    SetShowOthers(others);
    SetBase(base);
    // d.fClusterScale = "1.06 -0.003*x +0.0119*x*x";
    // Do the below if your input data does not contain these settings 
    if (sNN > 0) SetSNN(sNN);     // Collision energy per nucleon pair (GeV)
    if (sys > 0) SetSys(sys);     // Collision system (1:pp, 2:PbPB)
    if (trg > 0) SetTrigger(trg); // Collision trigger (1:INEL, 2:INEL>0, 4:NSD)
    if (eff > 0) SetTriggerEfficiency(eff); // Trigger efficiency
    if (vzMin < 999 && vzMax > -999) 
      SetVertexRange(vzMin,vzMax); // Collision vertex range (cm)
    SetCentralityRange(centMin,centMax); // Collision vertex range (cm)

    fCentSeen         = 0;
    fOptions          = flags;
    fFormats          = formats;
    SetForwardSysError(flags & kShowSysError ? 0.07 : 0);
    SetFinalMC        (flags & kUseFinalMC ? "forward_dndetamc.root" : "");
    // "EmpiricalCorrection.root"
    SetEmpirical      (flags & kUseEmpirical ? fEmpirical.Data() : "");
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
    
    if (!fEmpirical.IsNull()) {
      TParameter<bool>* p = 
	static_cast<TParameter<bool>*>(sums->FindObject("empirical"));
      if (p && p->GetVal() && !fEmpirical.IsNull()) {
	Warning("Run", "Empirical correction already applied");
	fEmpirical = "__task__";
      }
      else if (forward->FindObject("dndetaEmp")) {
	Warning("Run", "Empirical correction already applied");
	fEmpirical = "__task__";
      }
    } 

    if (!forward->FindObject("deltaIP")) {
      TH2* vertexXY = static_cast<TH2*>(sums->FindObject("vertexAccXY"));
      if (vertexXY && fDelta <= 0)
	SetDelta(vertexXY->GetMean(1), vertexXY->GetMean(2));
    }
    
    // --- Get information on the run --------------------------------
    FetchInformation(forward);

    // --- Print settings --------------------------------------------
    Info("Run", "Settings for the drawer:\n"
	 "   Show ratios:                      %5s\n"
	 "   Show Left/right:                  %5s\n"
	 "   Show rings:                       %5s\n"
	 "   Cut edges when rebinning:         %5s\n"
	 "   Remove outer rings:               %5s\n"
	 "   Force minimum bias:               %5s\n"
	 "   Mirror to un-covered regions:     %5s\n"
	 "   Export to file:                   %5s\n"
	 "   Add Zoom code:                    %5s\n"
	 "   Assume old format:                %5s\n"
	 "   Be verbose:                       %5s\n"
	 "   Hi-resolution plot:               %5s\n"
	 "   Extra whitespace:                 %5s\n"
	 "   Show logo:                        %5s\n"
	 "   Show clusters:                    %5s\n"
	 "   Show y-axis labels:               %5s\n"
	 "   Show other results:               0x%03x\n"
	 "   Rebinning factor:                 %5d\n"
	 "   Forward systematic error:         %5.1f%%\n"
	 "   Central systematic error:         %5.1f%%\n"
	 "   Trigger efficiency:               %5.1f%%\n"
	 "   Title on plot:                    %s\n"
	 "   Scaling of clusters to tracklets: %s\n"
	 "   Final MC correction file:         %s\n"
	 "   Empirical correction file:        %s",
	 ((fOptions & kShowRatios)    ? "yes" : "no"), 
	 ((fOptions & kShowLeftRight) ? "yes" : "no"),
	 ((fOptions & kShowRings)     ? "yes" : "no"),
	 ((fOptions & kCutEdges)      ? "yes" : "no"),
	 ((fOptions & kRemoveOuters)  ? "yes" : "no"),
	 ((fOptions & kForceMB)       ? "yes" : "no"),
	 ((fOptions & kMirror)        ? "yes" : "no"),
	 ((fOptions & kExport)        ? "yes" : "no"), 
	 ((fOptions & kAddExec)       ? "yes" : "no"), 
	 ((fOptions & kOldFormat)     ? "yes" : "no"), 
	 ((fOptions & kVerbose)       ? "yes" : "no"), 
	 ((fOptions & kHiRes)         ? "yes" : "no"), 
	 ((fOptions & kExtraWhite)    ? "yes" : "no"), 
	 ((fOptions & kLogo)          ? "yes" : "no"), 
	 ((fOptions & kNoCentral)     ? "no"  : "yes"),
	 ((fOptions & kNoLabels)      ? "no"  : "yes"),
	 fShowOthers, 
	 fRebin, 
	 (100*fFwdSysErr), 
	 (100*fCenSysErr), 
	 (100*fTriggerEff),
	 fTitle.Data(), 
	 fClusterScale.Data(), 
	 fFinalMC.Data(), 
	 fEmpirical.Data());

    // --- Set the macro pathand load other data script --------------
    TString savPath(gROOT->GetMacroPath());
    gROOT->SetMacroPath(Form("%s:$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2",
			     gROOT->GetMacroPath()));
    // Always recompile 
    if (!gROOT->GetClass("RefData"))
      gROOT->LoadMacro("OtherData.C+");
    gROOT->SetMacroPath(savPath);

    Bool_t useCen = !(fOptions & kNoCentral);
    // --- Get the central results -----------------------------------
    TList* clusters = 0;
    if (useCen) {
      clusters = static_cast<TList*>(file->Get("CentraldNdetaResults"));
      if (!clusters) Warning("Run", "Couldn't find list CentraldNdetaResults");
    }

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
    TObject* fwdEmp = 0;
    TObject* cenEmp = 0;
    TUrl     empUrl(fEmpirical);
    TFile*   empirical = 0;
    if (!fEmpirical.IsNull() && !fEmpirical.EqualTo("__task__")) {
      // Try some other locations
      const char* test[] = {
	fEmpirical.Data(),
	"file://${ANA_SRC}",
	"file://${FWD}",
	"file://${OADB_PATH}/PWGLF/FORWARD/EMPIRICAL",
	"file://${ALICE_PHYSICS}/OADB/PWGLF/FORWARD/EMPIRICAL",
	0 };
      const char** ptr = test;
      Bool_t ok = false;
      while (*ptr) {
	const char*  testf[] = { "", "empirical_000138190.root", 0 };
	const char** ptr2 = testf;
	while (*ptr2) {
	  TString path(gSystem->ConcatFileName(*ptr, *ptr2));
	  if (GetEmpirical(*ptr, useCen, fwdEmp, cenEmp, fEmpirical)) {
	    ok = true;
	    break;
	  }
	  ptr2++;
	}
	if (ok) break;
	ptr++;
      }
    }

    // --- Loop over input data --------------------------------------
    TObjArray truths;
    FetchTopResults(mcTruth,  0, 0, "MCTruth", max, rmax, amax,truths);
    TObjArray* fwdA = FetchTopResults(forward, forwardMC, fwdEmp, "Forward", 
				   max, rmax, amax,truths);
    TObjArray* cenA = FetchTopResults(clusters, 0, cenEmp, "Central", 
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
	TH1* cen = (cenA ? static_cast<TH1*>(cenA->At(i)) : 0);
	CorrectForward(fwd);
	CorrectCentral(cen);
	Double_t low, high;
	TH1* tmp = Merge(cen, fwd, low, high);
	TF1* f   = 0; // FitMerged(tmp, low, high);
	MakeSysError(tmp, cen, fwd, f);
	delete f;

	if (fOptions & kNoCentral && tmp) { 
	  // Split Sys error into two histograms 
	  const char* nme  = tmp->GetName();
	  TH1* tmpp = static_cast<TH1*>(tmp->Clone(Form("%s_a", nme)));
	  tmp->SetName(Form("%s_c", nme));
	  for (Int_t k = 1; k <= tmp->GetNbinsX(); k++) { 
	    Double_t x = tmp->GetXaxis()->GetBinCenter(k);
	    TH1* tmppp = (x < 0 ? tmpp : tmp);
	    tmppp->SetBinError(k, 0);
	    tmppp->SetBinContent(k, 0);
	  }
	  fResults->GetHists()->AddFirst(tmpp, (f ? "e5" : "e2"));
	}

	if (fOptions & kVerbose) 
	  Info("", "Adding systematic error histogram %s", tmp->GetName());
	fResults->GetHists()->AddFirst(tmp, (f ? "e5" : "e2"));

	if (!(fOptions & kMirror)) continue;

	TH1* tmp2 = Symmetrice(tmp);
	if (tmp2) {
	  tmp2->SetFillColor(tmp->GetFillColor());
	  tmp2->SetFillStyle(tmp->GetFillStyle());
	  tmp2->SetMarkerStyle(tmp->GetMarkerStyle());
	  tmp2->SetLineWidth(tmp->GetLineWidth());
	  fResults->GetHists()->AddFirst(tmp2, "e5");
	  fResults->Modified();
	}
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
    if (!fCentAxis) {
      TObject* cO = results->FindObject("centAxis");
      if (cO) {
	if (cO->IsA()->InheritsFrom(TH1::Class())) {
	  TH1* cH     = static_cast<TH1*>(cO);
	  fCentAxis   = cH->GetXaxis();
	}
	else if (cO->IsA()->InheritsFrom(TAxis::Class())) 
	  fCentAxis   = static_cast<TAxis*>(cO);
      }
    }
    if (!fCentMeth) 
      fCentMeth = results->FindObject("centEstimator");
    if (!fCentMeth && HasCent()) {
      TString tmp(fTitle);
      tmp.Remove(0,tmp.Index("CENT")+4);
      fCentMeth = new TNamed("centEstimator", tmp.Data()); // default");
      fCentMeth->SetUniqueID(1);
    }
    if (fSysString && fSysString->GetUniqueID() != 1 &&
	fCentMin == 0 && fCentMax == 0) {
      fCentMin = 0;
      fCentMax = 100;
    }

    if (fTriggerEff < 0) { 
      // Allow complete overwrite by passing negative number 
      SetTriggerEfficiency(TMath::Abs(fTriggerEff));
    }
    else if (fTriggerEff <= 0 || TMath::Abs(1-fTriggerEff)<1e-6) {
      TParameter<double>* eff = 
	static_cast<TParameter<double>*>(results->FindObject("triggerEff"));
      if (eff) {
	fTriggerEff = eff->GetVal();
	fExtTriggerEff = true;
	Info("FetchInformation", "External trigger efficicency: %5.3f",
	     fTriggerEff);
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
      Printf("Centrality axis: %d [%f,%f] vs [%f,%f]",
	     fCentAxis->GetNbins(),
	     fCentAxis->GetXmin(),
	     fCentAxis->GetXmax(),
	     fCentMin,
	     fCentMax);
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
	
    if (fTrigString) { 
      UInt_t mask = 0x200f & fTrigString->GetUniqueID();
      if      (mask == 0x2)      fTrigString->SetTitle("INEL>0");
      else if (fTriggerEff != 1) {
	if      (mask == 0x1)    fTrigString->SetTitle("INEL");
	else if (mask == 0x2000) fTrigString->SetTitle("NSD");
	else if (mask == 0x4)    fTrigString->SetTitle("NSD");
      }
    }

    if (true /*fOptions & kVerbose*/) {
      TString centTxt("none");
      TString centMeth(fCentMeth ? fCentMeth->GetTitle() : "V0M");
      if (fCentAxis) {
	Int_t nCent = fCentAxis->GetNbins();
	centTxt = Form("%s %d bins", centMeth.Data(), nCent);
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
    if (fSysString->GetUniqueID() == 3 ||
	fSysString->GetUniqueID() == 4) {
      Info("FetchInformation", "Left/Right asymmetry, and mirroring "
	   "explicitly disabled for pPb/Pbp");
      fOptions   &= ~kShowLeftRight;
      fOptions   &= ~kMirror;
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
    UShort_t cen   = (fCentMeth   ? fCentMeth->GetUniqueID() : 0);
    if (centLow < centHigh) {
      // Possibly modify trg according to method
      // Printf("We have centrality cen=%d meth=%s", cen,
      //        (fCentMeth ? fCentMeth->GetTitle() : "?"));
      UShort_t msk = 0x01; // Default to V0M comparison
      if (cen == 0 && fCentMeth) {
	TString cm(fCentMeth->GetTitle());
	if      (cm.EqualTo("V0M", TString::kIgnoreCase))    msk = 0x01;
	else if (cm.EqualTo("V0A", TString::kIgnoreCase))    msk = 0x02;
	else if (cm.EqualTo("ZNA", TString::kIgnoreCase))    msk = 0x04;
	else if (cm.EqualTo("ZNC", TString::kIgnoreCase))    msk = 0x08;
	else if (cm.EqualTo("V0C", TString::kIgnoreCase))    msk = 0x10;
	else if (cm.EqualTo("Cl1", TString::kIgnoreCase))    msk = 0x20;
      }
      else {
	switch (cen) {
	case 1:   msk = 0x01; break; // default V0M
	case 2:   msk = 0x01; break; // V0M
	case 3:   msk = 0x02; break; // V0A
	case 4:               break; // V0A123
	case 5:   msk = 0x10; break; // V0C
	case 6:               break; // FMD
	case 7:               break; // Tracks
	case 8:               break; // Tracklets
	case 9:               break; // CL0
	case 10:  msk = 0x20; break; // CL1
	case 11:              break; // CND
	case 12:  msk = 0x04; break; // ZNA
	case 13:  msk = 0x08; break; // ZNC
	case 14:              break; // ZPA
	case 15:              break; // ZPC
	case 16:              break; // NPA
	case 17:              break; // V0MvsFMD
	case 18:              break; // V0MvsTracklets
	case 19:              break; // ZEMvsZDC
	case 20:              break; // RefMult
	case 21:              break; // HMTF V0A
	case 22:              break; // HMTF V0M
	case 23:              break; // HMTF V0C
	default:              break;
	}
      }
      trg = (trg & 0x200f) | (msk << 4);
    }
    Long_t   ret   = 
      gROOT->ProcessLine(Form("RefData::GetData(%d,%d,%d,%d,%d,%d);",
			      sys,snn,trg,centLow,centHigh,fShowOthers));
    if (!ret) {
#if 0
      Warning("", "RefData::GetData(%d,%d,0x%x,%d,%d,0x%x);",
	      sys,snn,trg,centLow,centHigh,fShowOthers);
      Warning("FetchOthers", 
	      "No other data for %s %s %s %3d%%-%3d%% central (0x%x)", 
	      fSysString  ? fSysString->GetTitle()  : "unknown", 
	      fTrigString ? fTrigString->GetTitle() : "unknown", 
	      fSNNString  ? fSNNString->GetTitle()  : "unknown", 
	      centLow, centHigh, fShowOthers);
#endif
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
		  TObject*      empCorr,
		  const char*   name, 
		  Double_t&     max,
		  Double_t&     rmax,
		  Double_t&     amax,
		  TObjArray&    truths)
  {
    if (!list) return 0;
    UShort_t   n = HasCent() ? fCentAxis->GetNbins() : 1;
    // // Info("FetchTopResults","got %d centrality bins", n);
    // if (n == 0) {
    //   TH1* h  = FetchOne(list, mcList, empCorr, name, "all",
    // 			 FetchOthers(0,0), -1000, 0, 
    // 			 max, rmax, amax, truths);
    //   if (!h) return 0;
    //   TObjArray* a = new TObjArray;
    //   // Info("FetchTopResults", "Adding %s to result stack", h->GetName());
    //   a->AddAt(h, 0);
    //   return a;
    // }
    
    TObjArray* a = new TObjArray;
    truths.Expand(n);
    for (UShort_t i = 0; i < n; i++) { 
      Float_t  centLow  = 0;
      Float_t  centHigh = 0;
      TString  lname    = "all";
      Int_t    col      = -1000;
      TString  centTxt  = "";
      if (HasCent()) {
	centLow  = fCentAxis->GetBinLowEdge(i+1);
	centHigh = fCentAxis->GetBinUpEdge(i+1);
	lname    = Form("cent%03dd%02d_%03dd%02d",
			Int_t(centLow),  Int_t(centLow *100)%100,
			Int_t(centHigh), Int_t(centHigh*100)%100);
	col      = GetCentralityColor(i+1);
        centTxt  = Form("%6.2f%%-%6.2f%% central", centLow, centHigh);
      }
      TH1* tt = static_cast<TH1*>(truths.At(i));
      TH1* ot = tt;
      TH1* h  = FetchOne(list, mcList, empCorr, name, lname,
			 FetchOthers(centLow,centHigh), col, 
			 centTxt.Data(), max, rmax, amax, tt);
      if (!h) continue;
      fCentSeen |= (1 << i);

      if (tt != ot) { 
	truths.AddAt(tt, i);
      }
      // Info("FetchTopResults", "Adding %pto result stack", h);
      a->AddAt(h, i);
    }
    if (a->GetEntries() <= 0) {
      delete a;
      a = 0;
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
		TObject*      empCorr,
		const char*   name, 
		const char*   folderName,
		TMultiGraph*  others, 
		Int_t         col,
		const char*   txt,
		Double_t&     max,
		Double_t&     rmax,
		Double_t&     amax,
		TH1*&         truth)
  {
    TString foldName(folderName);
    TList* folder = (fOptions & kOldFormat ? const_cast<TList*>(list) :
		     static_cast<TList*>(list->FindObject(foldName)));
    if (!folder) {
      Warning("FetchOne",
	      "Couldn't find list '%s' in %s, trying w/o decimals (%s)", 
	      foldName.Data(), list->GetName(), txt);
      foldName.ReplaceAll("d00", "");
      foldName.ReplaceAll("d10", "");
      folder = (fOptions & kOldFormat ? const_cast<TList*>(list) :
		static_cast<TList*>(list->FindObject(foldName)));
      if (!folder) {
	Error("FetchOne", "Couldn't find list '%s' in %s", 
	      foldName.Data(), list->GetName());	
	return 0;
      }
    }
    TList* mcFolder = 0;
    if (mcList) {
      mcFolder = static_cast<TList*>(mcList->FindObject(folderName));
      if (!mcFolder) 
	Warning("FetchOne", 
		"Didn't find the list '%s' in %s for final MC correction", 
		folderName, mcList->GetName());
    }
    if (fOptions & kVerbose) {
      TObject* normCalc = folder->FindObject("normCalc");
      if (normCalc) 
	Info("FetchOne", "%s:\n%s", 
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
			TObject*      empCorr,
			const char*   name, 
			TMultiGraph*  thisOther,
			Int_t         color,
			const char*   centTxt,
			Double_t&     max,
			Double_t&     rmax,
			Double_t&     amax, 
			TH1*&         truth)
  {
    TH1* norm        = FetchHistogram(list, Form("norm%s", name));
    if (!norm) return 0;
    if (norm->GetMaximum() < 1000) {
      Warning("FetchCentResults", "Too few events in %s: %ld",
	      list->GetName(), Long_t(norm->GetMaximum()));
      // return 0;
    }
    
    TH1* dndeta      = FetchHistogram(list, Form("dndeta%s", name));
    TH1* dndetaMC    = FetchHistogram(list, Form("dndeta%sMC", name));
    TH1* dndetaTruth = FetchHistogram(list, "dndetaTruth");
    TH1* dndetaEmp   = FetchHistogram(list, Form("dndeta%sEmp", name));
    // Info("", "dN/deta truth from %s: %p", list->GetName(), dndetaTruth);
    // Info("", "dN/deta truth from external: %p", truth);
    if (dndetaEmp) dndeta = dndetaEmp;

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
      } // If show rings
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
      } // if others for this 
    } // if not truth 
    if (dndetaMC) { 
      fRatios->Add(Ratio(dndeta,    dndetaMC,    rmax));
      fRatios->Add(Ratio(dndetaSym, dndetaMCSym, rmax));
    }
    if (truth) {
      Info("", "Forming ratio to truth:\n\t%s\n\t%s",
	   dndeta->GetName(), 
	   truth->GetName());
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
  void CorrectEmpirical(TH1* dndeta, TObject* empObj) 
  {
    if (!dndeta) return;
    if (!empObj) return;
   
    
    if (empObj->IsA()->InheritsFrom(TGraphAsymmErrors::Class())) {
      Info("CorrectEmpirical", "Doing empirical correction of dN/deta");
      TGraphAsymmErrors* empCorr = static_cast<TGraphAsymmErrors*>(empObj);
      TAxis* xAxis = dndeta->GetXaxis();
      for (Int_t i = 1; i <= xAxis->GetNbins(); i++) {
	Double_t x = xAxis->GetBinCenter(i);
	Double_t y = dndeta->GetBinContent(i);
	Double_t c = empCorr->Eval(x);
	dndeta->SetBinContent(i, y / c);
      }
    }
    else if (empObj->IsA()->InheritsFrom(TH1::Class())) {
      Info("CorrectEmpirical", "Doing empirical correction of dN/deta");
      TH1* empCorr = static_cast<TH1*>(empObj);
      dndeta->Divide(empCorr);
    }
    else { 
      Warning("CorrectEmpirical", 
	      "Don't know how to apply a %s as an empirical correction",
	      empObj->IsA()->GetName());
    }
  }
  //__________________________________________________________________
  void CorrectTriggerEff(TH1* dndeta)
  {
    if (fExtTriggerEff) return;
    if (!dndeta) return;
    if (fTriggerEff <= 0) return; //  || fTriggerEff >= 1) return;
    Info("CorrectTriggerEff", "Correcting with trigger efficiency %5.3f",
	 fTriggerEff);
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
      y2 = (y11 > 0.0001 ? 0.2 : 0.2);
    }
    TCanvas* c = new TCanvas("Results", "Results", w, h);
    c->SetFillColor(0);
    c->SetBorderSize(0);
    c->SetBorderMode(0);

    Double_t s = Float_t(h) / 900;

    PlotResults(max, y1, s);
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
   * Plot the title on a pad 
   * 
   * @param p       Pad to draw in
   * @param yd      Division height
   * @param bottom  Bottom or top of pad 
   */
  void PlotTitle(TVirtualPad* p, Double_t yd, Bool_t bottom=true) 
  {
    // Put a title on top
    p->cd();
    fTitle.ReplaceAll("@", " ");
    Double_t s = 1/yd/1.2;
    TLatex* tit = new TLatex((bottom ? kRightMargin : p->GetLeftMargin()),
			     (bottom ? 0.01 : .99), fTitle.Data());
    tit->SetNDC();
    tit->SetTextFont(kFont);
    tit->SetTextAlign(bottom ? 11 : 13);
    tit->SetTextSize(s*0.045);
    tit->SetTextColor(kAlicePurple);
    tit->Draw();
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
    Bool_t   mirrorSeen = false;
    while ((hist = static_cast<TH1*>(next()))) { 
      TString t(hist->GetTitle());
      TString n(hist->GetName());
      n.ToLower();
      if (t.Contains("mirrored")) { mirrorSeen = true; continue; }
      if (n.Contains("syserror")) { sysErrSeen = true; continue; }
      if (unique.FindObject(t.Data())) continue;
      // TObjString* s1 = new TObjString(hist->GetTitle());
      TParameter<float>* s1 = new TParameter<float>(t,hist->GetMarkerSize());
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
	// TObjString* s2 = new TObjString(n);
	TParameter<float>* s2 = new TParameter<float>(n,g->GetMarkerSize());
	s2->SetUniqueID(((g->GetMarkerStyle() & 0xFFFF) << 16) |
			 ((g->GetMarkerColor() & 0xFFFF) <<  0));
	unique.Add(s2);
	// l->AddEntry(hist, hist->GetTitle(), "pl");
      }
    }

    // Add legend entries for unique items only
    TIter nextu(&unique);
    // TObject* s = 0;
    TParameter<float>* s = 0;
    Int_t i = 0;
    while ((s = static_cast<TParameter<float>*>(nextu()))) { 
      TLegendEntry* dd = l->AddEntry(Form("data%2d", i++), 
				     s->GetName(), "lp");
      Int_t style = (s->GetUniqueID() >> 16) & 0xFFFF;
      Int_t color = (s->GetUniqueID() >>  0) & 0xFFFF;
      dd->SetLineColor(kBlack);
      if (HasCent()) dd->SetMarkerColor(kBlack);
      else           dd->SetMarkerColor(color);
      dd->SetMarkerStyle(style);
      dd->SetMarkerSize(s->GetVal());
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
    if (mirrorSeen /* (fOptions & kMirror) */) {
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
  const char* CentLimitName(Bool_t isMult, Float_t v)
  {
    if (isMult) return Form("%3d", Int_t(v));
    if ((Int_t(v*100) % 100) == 0) return Form("%3d%%", Int_t(v));
    if ((Int_t(v*100)  % 10) == 0) return Form("%5.1f%%", v);
    return Form("%6.2f%%", v);
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

    Printf("Centralities seen 0x%x", fCentSeen);
    
    if (fCentAxis->GetNbins() <= 4) y1 += .15;
    TLegend* l = new TLegend(x1,y1,x2,y2);
    l->SetNColumns(1);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextFont(kFont);
    l->SetTextColor(kAliceBlue);

    TString     centMeth(fCentMeth->GetTitle());
    Bool_t      isMult   = centMeth.EqualTo("MULT");
    Int_t       lowOff   = (isMult ? 1 : 0);
    Int_t       nextOff  = 0;
    Int_t n = fCentAxis->GetNbins();
    for (Int_t i = 1; i <= n; i++) {
      if (!(fCentSeen & (1 << (i-1)))) { nextOff = 0; continue; }
      Double_t low  = fCentAxis->GetBinLowEdge(i) + lowOff + nextOff;
      Double_t upp  = fCentAxis->GetBinUpEdge(i);
      TString  txt  = Form("%s - %s",
			   CentLimitName(isMult,low),
			   CentLimitName(isMult,upp));
      if (isMult && upp == low) {
	nextOff = -1;
	txt     = CentLimitName(isMult,low-1);
      }
      else
	nextOff = 0;
      TLegendEntry* e = l->AddEntry(Form("dummy%02d", i), txt, "pl");
      e->SetMarkerColor(GetCentralityColor(i));
    }
    l->Draw();
  }
  //__________________________________________________________________
  void AttachExec(TVirtualPad* p, const char* plot, UShort_t id, 
		  Bool_t isBottom) 
  {
    if (!(fOptions & kAddExec)) return;

    if (isBottom) {
      fRangeParam->fMasterAxis = FindXAxis(p, plot);
      p->AddExec("range", Form("RangeExec((dNdetaDrawer::RangeParam*)%p)", 
				fRangeParam));
    }
    else { 
      if (id == 1) {
	fRangeParam->fSlave1Axis = FindXAxis(p, plot);
	fRangeParam->fSlave1Pad  = p;
      }
      else if (id == 2) {
	fRangeParam->fSlave2Axis = FindXAxis(p, plot);
	fRangeParam->fSlave2Pad  = p;
      }
    }
  }
  //__________________________________________________________________
  /** 
   * Plot the results
   *    
   * @param max       Maximum 
   * @param yd        Bottom position of pad 
   * @param s         Scaling 
   */
  void PlotResults(Double_t max, Double_t yd, Double_t s) 
  {
    // --- Make a sub-pad for the result itself ----------------------
    TPad* p1 = new TPad("p1", "p1", 0, yd, 1.0, 1.0, 0, 0, 0);
    p1->SetTopMargin(kRightMargin);
    p1->SetBorderSize(0);
    p1->SetBorderMode(0);
    p1->SetBottomMargin(yd > 0.001 ? 0.001 : 0.1);
    p1->SetRightMargin(kRightMargin);
    if ((fOptions & kShowLeftRight) || (fOptions & kShowRatios)) p1->SetGridx();
    p1->SetTicks(1,1);
    p1->SetNumber(1);
    p1->Draw();
    p1->cd();

    // --- Figure out min/max ----------------------------------------
    // Info("PlotResults", "Plotting results with max=%f", max);
    fResults->SetMaximum((fOptions & kExtraWhite ? 1.4 : 1.15)*max);
    fResults->SetMinimum(yd > 0.00001 ? -0.02*max : 0);
    // fResults->SetMinimum(yd > 0.00001 ? -0.02*max : 0);

    // --- Fix up axis -----------------------------------------------
    Double_t yyd  = (1-yd)*(yd > .001 ? 1 : .9 / 1.2);
    FixAxis(fResults, yyd, 
	    "1/#it{N}#kern[.1]{d#it{N}_{ch}/d#it{#eta}}"
	    //"#frac{1}{#it{N}}#kern[.1]{#frac{d#it{N}_{ch}}{d#it{#eta}}}"
	    );

    // --- Fix up marker size ---------------------------------------
    TIter next(fResults->GetHists());
    TH1*  h = 0;
    while ((h = static_cast<TH1*>(next()))) 
      h->SetMarkerSize(h->GetMarkerSize()*s);
    // --- Clear pad and re-draw ------------------------------------
    p1->Clear();
    fResults->DrawClone("nostack e1");

    // --- Draw other data -------------------------------------------
    if (fShowOthers != 0) {
      TGraphAsymmErrors* o      = 0;
      TIter              nextG(fOthers->GetListOfGraphs());
      while ((o = static_cast<TGraphAsymmErrors*>(nextG()))) {
	Double_t gS = s;
	switch (o->GetMarkerStyle()) { 
	case 29: 
	case 30: gS *= 1.2; break; // Star		
	case 27: 
	case 33: gS *= 1.2; break; // Diamond
	}
	o->SetMarkerSize(o->GetMarkerSize()*gS);
	  
        o->DrawClone("same p");
      }
    }

    // --- Make a legend in the result pad ---------------------------
    Double_t x1 = p1->GetLeftMargin()+.08;
    Double_t x2 = 1-p1->GetRightMargin()-.08;
    Double_t y1 = p1->GetBottomMargin()+.01;
    Double_t y2 = .35;
    Int_t    nC = 2;
    if (HasCent()) { 
      if (fCentAxis->GetNbins() <= 4) {
	x1 = p1->GetLeftMargin()+.15;
	x2 = 1-p1->GetRightMargin()-.15;
	y2 = .2;
      }
      else {
	nC = 1;
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
    }	
    BuildLegend(fResults, fOthers, x1, y1, x2, y2, nC);

    // --- Parameters for stuff on the right -------------------------
    Double_t yTop      = 1-p1->GetTopMargin()-.02;
    Double_t xR        = .95;
    Double_t yR        = yTop;

    // --- Put a nice label in the plot ------------------------------
    TString     eS;
    UShort_t    snn = fSNNString->GetUniqueID();
    const char* sys = fSysString->GetTitle();
    if (snn == 2750) snn = 2760;
    if (snn < 1000)           eS = Form("%3dGeV", snn);
    else if (snn % 1000 == 0) eS = Form("%dTeV", snn/1000);
    else                      eS = Form("%.2fTeV", float(snn)/1000);
    Bool_t nn = (fSysString->GetUniqueID() != 1);
    TString tS(fTrigString->GetTitle());
    if (HasCent()) {
      TString centMeth(fCentMeth ? fCentMeth->GetTitle() : "CENTV0M");
      if (centMeth.EqualTo("MULT"))
	tS = "by N_{#lower[-.2]{ch}} |#it{#eta}|<0.8";
      else 
	tS  = "by centrality";
      if (fCentMeth && fCentMeth->GetTitle()[0] != '\0')
	tS.Append(Form(" (%s)", fCentMeth->GetTitle()));
#if 0
      UShort_t trg = fTrigString->GetUniqueID();
      switch (trg) { 
      case 0x10: tS.Append(" (V0M)"); break;
      case 0x20: tS.Append(" (V0A)"); break;
      case 0x40: tS.Append(" (ZNA)"); break;
      case 0x80: tS.Append(" (ZNC)"); break;
      }
#endif
    }
    
    TLatex* tt = new TLatex(xR, yR, Form("%s #sqrt{s%s}=%s, %s", 
					 sys, (nn ? "_{NN}" : ""),
					 eS.Data(), tS.Data()));
    tt->SetTextColor(kAliceBlue);
    tt->SetNDC();
    tt->SetTextFont(kFont);
    tt->SetTextAlign(33);
    tt->Draw();
    yR -= tt->GetTextSize() + .01;
    
    if (fSysString->GetUniqueID() == 3) { 
      // pPb - put arrow at y_CM
      UShort_t a1  = 1;
      UShort_t a2  = 208;
      UShort_t z1  = 1;
      UShort_t z2  = 82;
      Double_t yCM = .5 * TMath::Log(Float_t(z1*a2)/z2/a1);
      TArrow*  a   = new TArrow(yCM, 0, yCM, 0.05*max, 0.01, "<|");
      a->SetAngle(30);
      a->SetFillColor(kAliceBlue);
      a->SetLineColor(kAliceBlue);
      a->Draw();
    }

    // --- Put number of accepted events on the plot -----------------
    Int_t nev = 0;
    if (fTriggers) nev = fTriggers->GetBinContent(1);
    TLatex* et = new TLatex(xR, yR, Form("%d events", nev));
    et->SetTextColor(kAliceBlue);
    et->SetNDC();
    et->SetTextFont(kFont);
    et->SetTextAlign(33);
    et->Draw();
    yR -= et->GetTextSize() + .01;

    // --- Put vertex axis on the plot -------------------------------
    if (fVtxAxis) { 
      TLatex* vt = new TLatex(xR, yR, fVtxAxis->GetTitle());
      vt->SetNDC();
      vt->SetTextFont(kFont);
      vt->SetTextAlign(33);
      vt->SetTextColor(kAliceBlue);
      vt->Draw();
      yR -= vt->GetTextSize() + .01;
    }
    // results->Draw("nostack e1 same");

    // --- Put statement on corrections used on the plot -------------
    TString corrs;
    if (!fEmpirical.IsNull()) corrs.Append("Emperical");
    if (fDelta > 0) {
      if (!corrs.IsNull()) corrs.Append("+");
      corrs.Append("IP_{xy}");
    }
    if (!fFinalMC.IsNull())   {
      if (!corrs.IsNull()) corrs.Append("+");
      corrs.Append("Final MC");
    }

    if (!corrs.IsNull()) {
      corrs.Append(" correction");
      if (corrs.Index("+") != kNPOS) corrs.Append("s");
      TLatex* em = new TLatex(xR, yR, corrs);
      em->SetNDC();
      em->SetTextFont(kFont);
      em->SetTextAlign(33);
      em->SetTextColor(kAliceBlue);
      em->Draw();
      yR -= em->GetTextSize() + .01;
    }

    // --- Put trigger efficiency on the plot (only pp and if != 1) --
    if (fTriggerEff > 0 && fTriggerEff < 1 && !HasCent()) { 
      TLatex* ef = new TLatex(xR, yR, Form("#varepsilon_{%s} = %5.3f", 
					 fTrigString->GetTitle(), 
					 fTriggerEff));
      ef->SetNDC();
      ef->SetTextFont(kFont);
      ef->SetTextAlign(33);
      ef->SetTextColor(kAliceBlue);
      ef->Draw();
      yR -= ef->GetTextSize() + .01;
    }

    // --- Put logo on the plot --------------------------------------
    if (fOptions & kLogo) {
      TString savPath(gROOT->GetMacroPath());
      TString fwd("$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2");
      gROOT->SetMacroPath(Form("%s:%s/scripts",
			       gROOT->GetMacroPath(), fwd.Data()));
      // Always recompile 
      if (!gROOT->GetClass("AliceLogo"))
	gROOT->LoadMacro("AliceLogo.C+");
      gROOT->SetMacroPath(savPath);
      
      if (gROOT->GetClass("AliceLogo")) {
	yR -= .22;
	p1->cd();
	p1->Range(0,0,1,1);
	gROOT->ProcessLine("AliceLogo* al = new AliceLogo();");
	gROOT->ProcessLine(Form("al->Draw(0,.88,%f,.2, 0, 0);", yR));
	yR -= .01;
      }
    }

    
    // --- Parameter for stuff on the left ---------------------------
    Double_t yL = yTop;
    Double_t xL = .12;

    // --- Mark as work in progress ----------------------------------
    TLatex* pt = new TLatex(xL, yL, "Work in progress");
    pt->SetNDC();
    pt->SetTextFont(62);
    // pt->SetTextSize();
    pt->SetTextColor(kAliceRed);
    pt->SetTextAlign(13);
    pt->Draw();
    yL -= pt->GetTextSize()+.01;

    TDatime now;    
    TLatex* dt = new TLatex(xL, yL, now.AsSQLString());
    dt->SetNDC();
    dt->SetTextFont(42);
    dt->SetTextSize(0.04);
    dt->SetTextColor(kAliceBlue); // kAliceRed);
    dt->SetTextAlign(13);
    dt->Draw();
    yL -= dt->GetTextSize()+.01;

    // --- Possible centrality legend --------------------------------
    if (fSysString->GetUniqueID() == 1) { xL += .2; yL = y1; }
    BuildCentLegend(xL, yL-.4, xL+.23, yL);

    // --- Attach Zoom executor --------------------------------------
    AttachExec(p1, fResults->GetName(), 1, false);

    // --- Possibly add title ----------------------------------------
    if (yd < 0.0001) PlotTitle(p1, yyd, true);

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

    // --- Make a sub-pad for the ratios -----------------------------
    TPad* p2 = new TPad("p2", "p2", 0, y1, 1.0, y2, 0, 0, 0);
    p2->SetTopMargin(0.001);
    p2->SetRightMargin(kRightMargin);
    p2->SetBottomMargin(isBottom ? 1/yd * 0.07 : 0.0001);
    p2->SetGridx();
    p2->SetTicks(1,1);
    p2->SetNumber(2);
    p2->Draw();
    p2->cd();

    // --- Fix up axis -----------------------------------------------
    FixAxis(fRatios, yd, "Ratios", 7);

    // --- Set up min/max --------------------------------------------
    fRatios->SetMaximum(1+TMath::Max(.22,1.05*max));
    fRatios->SetMinimum(1-TMath::Max(.32,1.05*max));

    // --- Clear pad and draw ----------------------------------------
    p2->Clear();
    fRatios->DrawClone("nostack e1");

    
    // --- Make a legend ---------------------------------------------
    BuildLegend(fRatios, 0, .15,p2->GetBottomMargin()+.01,.9,
		isBottom ? .6 : .4, 2);

    // --- Make a nice band from 0.9 to 1.1 --------------------------
    TGraphErrors* band = new TGraphErrors(2);
    band->SetPoint(0, fResults->GetXaxis()->GetXmin(), 1);
    band->SetPoint(1, fResults->GetXaxis()->GetXmax(), 1);
    band->SetPointError(0, 0, (fFwdSysErr > 0 ? fFwdSysErr : .1));
    band->SetPointError(1, 0, (fFwdSysErr > 0 ? fFwdSysErr : .1));
    band->SetFillColor(kYellow+2);
    band->SetFillStyle(3002);
    band->SetLineStyle(2);
    band->SetLineWidth(1);
    band->Draw("3 same");
    band->DrawClone("X L same");
    
    // --- Replot the ratios on top ----------------------------------
    fRatios->DrawClone("nostack e1 same");

    // --- Some more stuff -------------------------------------------
    AttachExec(p2, fRatios->GetName(), 2, isBottom);
    if (isBottom) PlotTitle(p2, yd, true);
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
    // --- Make a sub-pad for the asymmetry --------------------------
    TPad* p3 = new TPad("p3", "p3", 0, y1, 1.0, y2, 0, 0, 0);
    p3->SetTopMargin(0.001);
    p3->SetRightMargin(kRightMargin);
    p3->SetBottomMargin(isBottom ? 1/yd * 0.07 : 0.0001);
    p3->SetGridx();
    p3->SetTicks(1,1);
    p3->SetNumber(2);
    p3->Draw();
    p3->cd();

    // --- Fix up axis -----------------------------------------------
    FixAxis(fLeftRight, yd, "Right/Left", 4);

    // ---- Setup min/max --------------------------------------------
    fLeftRight->SetMaximum(1+TMath::Max(.12,1.05*max));
    fLeftRight->SetMinimum(1-TMath::Max(.15,1.05*max));

    // --- Clear pad and redraw --------------------------------------
    p3->Clear();
    fLeftRight->DrawClone("nostack e1");

    
    // --- Make a legend ---------------------------------------------
    Double_t xx1 = (HasCent() ? .7                           : .15); 
    Double_t xx2 = (HasCent() ? 1-p3->GetRightMargin()-.01   : .90);
    Double_t yy1 = p3->GetBottomMargin()+.01;
    Double_t yy2 = (HasCent() ? 1-p3->GetTopMargin()-.01-.15 : .5);
    BuildLegend(fLeftRight, 0, xx1, yy1, xx2, yy2);

    // --- Make a nice band from 0.9 to 1.1 --------------------------
    TGraphErrors* band = new TGraphErrors(2);
    band->SetPoint(0, fResults->GetXaxis()->GetXmin(), 1);
    band->SetPoint(1, fResults->GetXaxis()->GetXmax(), 1);
    band->SetPointError(0, 0, (fFwdSysErr > 0 ? fFwdSysErr : .05));
    band->SetPointError(1, 0, (fFwdSysErr > 0 ? fFwdSysErr : .05));
    band->SetFillColor(kYellow+2);
    band->SetFillStyle(3002);
    band->SetLineStyle(2);
    band->SetLineWidth(1);
    band->Draw("3 same");
    band->DrawClone("X L same");

    // --- Re-draw over ----------------------------------------------
    fLeftRight->DrawClone("nostack e1 same");

    // --- Misc stuff ------------------------------------------------
    AttachExec(p3, fLeftRight->GetName(), 2, isBottom);
    if (isBottom) PlotTitle(p3, yd, true);
  }
  /** @} */
  //==================================================================
  /** 
   * @{ 
   * @name Data utility functions 
   */
  Color_t Brighten(Color_t origNum, Int_t nTimes=2) const
  {
    TColor* col   = gROOT->GetColor(origNum);
    if (!col) return origNum;
    Int_t   origR = Int_t(0xFF * col->GetRed());
    Int_t   origG = Int_t(0xFF * col->GetGreen());
    Int_t   origB = Int_t(0xFF * col->GetBlue());
    // TColor::Pixel2RGB(TColor::Number2Pixel(origNum), origR, origG, origB);
    Int_t off    = nTimes*0x33;
    Int_t newR   = TMath::Min((origR+off),0xff);
    Int_t newG   = TMath::Min((origG+off),0xff);
    Int_t newB   = TMath::Min((origB+off),0xff);
    Int_t newNum = TColor::GetColor(newR, newG, newB);
    return newNum;
  }
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
    // Double_t centLow  = fCentAxis->GetBinLowEdge(bin);
    Double_t max      = fCentAxis->GetXmax();
    Double_t min      = fCentAxis->GetXmin();
    Double_t centHigh = fCentAxis->GetBinUpEdge(bin);
    Double_t centLow  = fCentAxis->GetBinUpEdge(bin);
    Float_t  fc       = ((centHigh+centLow)/2-min) / (max-min);
    Int_t    nCol     = gStyle->GetNumberOfColors();
    Int_t    icol     = TMath::Min(nCol-1,int(fc * nCol + .5));
    Int_t    col      = gStyle->GetColorPalette(icol);
    //Info("GetCentralityColor","%3d: %3d-%3d -> %3d",bin,centLow,centHigh,col);
    return col;
    // return Brighten(col);
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
    h->SetLineColor(color);
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
    if (!h) return 0;
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
  /** 
   * Get error bands from vanilla graph
   *
   * @param g 
   * @param low
   * @param up
   */
  static void ErrorGraphs(const TGraph* g, TGraph*& low, TGraph*& up)
  {
    Warning("ErrorGraphs", "Called with vanila TGraph (%s)",
	    g->IsA()->GetName());
    Int_t n = g->GetN();
    low = new TGraph(n);
    up  = new TGraph(n);
    for (Int_t i = 0; i < n; i++) {
      low->SetPoint(i, 0, 0);
      up->SetPoint(i, 0, 0);
    }
  }
  /** 
   * Get error bands 
   */
  static void ErrorGraphs(const TGraphErrors* g, TGraph*& low, TGraph*& up)
  {
    // Info("ErrorGraphs", "Called with TGraphErrors");
    Int_t n = g->GetN();
    low = new TGraph(n);
    up  = new TGraph(n);
    for (Int_t i = 0; i < n; i++) {
      low->SetPoint(i, g->GetX()[i], g->GetY()[i] - g->GetEY()[i]);
      up->SetPoint(i, g->GetX()[i], g->GetY()[i] + g->GetEY()[i]);
    }
  }
  /** 
   * Get error bands 
   */
  static void ErrorGraphs(const TGraphAsymmErrors* g, TGraph*& low, TGraph*& up)
  {
    // Info("ErrorGraphs", "Called with TGraphAsymmErrors");
    Int_t n = g->GetN();
    low = new TGraph(n);
    up  = new TGraph(n);
    for (Int_t i = 0; i < n; i++) {
      low->SetPoint(i, g->GetX()[i], g->GetY()[i] - g->GetEYlow()[i]);
      up->SetPoint(i, g->GetX()[i], g->GetY()[i] + g->GetEYhigh()[i]);
    }
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

    bool mirror = false;
    TString n1(o1->GetName());
    if (n1.Contains("mirror")) mirror = true;

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
      UShort_t m2bits = MarkerUtil::GetMarkerBits(m2->GetMarkerStyle());
      if (mirror) m2bits |= MarkerUtil::kHollow;
      else        m2bits &= ~MarkerUtil::kHollow;
      r->SetMarkerStyle(MarkerUtil::GetMarkerStyle(m2bits));
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
    Int_t    n     = g->GetN();
    Double_t xlow  = g->GetX()[0];
    Double_t xhigh = g->GetX()[n-1];
    TGraph*  glow  = 0;
    TGraph*  ghigh = 0;
    if (g->IsA()->InheritsFrom(TGraphErrors::Class())) { 
      const TGraphErrors* ge = static_cast<const TGraphErrors*>(g);
      xlow  -= ge->GetErrorXlow(0);
      xhigh += ge->GetErrorXhigh(n-1);
      ErrorGraphs(ge, glow, ghigh);
    }
    if (g->IsA()->InheritsFrom(TGraphAsymmErrors::Class())) { 
      const TGraphAsymmErrors* ge = static_cast<const TGraphAsymmErrors*>(g);
      xlow  -= ge->GetErrorXlow(0);
      xhigh += ge->GetErrorXhigh(n-1);
      ErrorGraphs(ge, glow, ghigh);
    }
    if (xlow > xhigh) { Double_t t = xhigh; xhigh = xlow; xlow = t; }

    for (Int_t i = 1; i <= h->GetNbinsX(); i++) { 
      Double_t c = h->GetBinContent(i);
      if (c <= 0) continue;

      Double_t x = h->GetBinCenter(i);
      if (x < xlow || x > xhigh) continue; 

      Double_t f = g->Eval(x);
      if (f <= 0) continue;

      Double_t efl = f - glow->Eval(x);
      Double_t efh = ghigh->Eval(x) - f;
      Double_t ec  = h->GetBinError(i);
      Double_t e   = TMath::Max(ec, TMath::Max(efl, efh));

      ret->SetBinContent(i, c / f);
      ret->SetBinError(i, e / f);
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
    Bool_t supLabel = (fResults == stack && ((fOptions & kNoLabels)));
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
      ya->SetTitleOffset(/*1.15*/1.4*ya->GetTitleOffset()/s);
      ya->SetLabelSize(supLabel ? 0 : s*ya->GetLabelSize());
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
    if (!fwd) return 0;
    TH1* tmp = static_cast<TH1*>(fwd->Clone("tmp"));
    TString name(fwd->GetName());
    name.ReplaceAll("Forward", "Merged");
    tmp->SetName(name);

    // tmp->SetMarkerStyle(28);
    // tmp->SetMarkerColor(kBlack);
    tmp->SetDirectory(0);
    if (!cen) return tmp;

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
    if (!tmp) return 0;
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
    if (!tmp || !fwd) return;
    for (Int_t i = 1; i <= tmp->GetNbinsX(); i++) {
      Double_t tc = tmp->GetBinContent(i);
      if (tc < 0.01) continue;
      Double_t fc     = fwd->GetBinContent(i);
      Double_t cc     = cen ? cen->GetBinContent(i) : 0;
      Double_t sysErr = fFwdSysErr;
      Double_t mc     = fc;
      if (cc > .01 && fc > 0.01) {
	sysErr = (fFwdSysErr+fCenSysErr) / 2;
	mc     = (fc+cc) / 2;
      }
      else if (cc > .01) {
	sysErr = fCenSysErr;
	mc     = cc;
      }
      Double_t x = tmp->GetXaxis()->GetBinCenter(i);
      Double_t y = (fit ? fit->Eval(x) : mc);
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
    if (!h) return;
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
	 << "// Setting for this draw\n"
	 << "//\n" 
	 << "//  Draw(\"\",\"" << fTitle << "\"," << fRebin << ",0x"
	 << std::hex << fShowOthers << ",0x" << fOptions << ","
	 << std::dec << fSNNString->GetUniqueID() << ","
	 << fSysString->GetUniqueID() << ","
	 << std::hex << fTrigString->GetUniqueID() << ","
	 << std::dec << fTriggerEff << "," << fCentMin << "," << fCentMax << ","
	 << fVtxAxis->GetXmin() << "," << fVtxAxis->GetXmax() << ",\""
	 << fBase << "\",0x" << std::hex << fFormats
	 << std::dec << ");\n"
	 << "//\n"
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

    TString tgt(fRebin > 1 ? "rebin" : "full");
    TString clean(Form("sed -e 's/_[0-9a-f]\\{4\\}"
		       "\\(\\|_\\{1,4\\}[0-9a-f]\\{1,4\\}\\)\\b//g'"
		       " -e 's/%s/%s/g' < %s > %s.C",
		       bname.Data(), tgt.Data(), fname.Data(), tgt.Data()));
    Printf("Execute \"%s\"", clean.Data());
    gSystem->Exec(clean);
    TString trg(fTrigString->GetTitle());
    TString mth;
    if (HasCent()) {
      trg = "CENT"; // Form("CENT%s", fCentMeth->GetTitle());
      mth = fCentMeth->GetTitle();
    }
    Printf("Copy %s.C to %s/%s/%05d/%s%s/%s.C",tgt.Data(),
	   (fEmpirical.IsNull() ? "normal" : "nosec"),
	   fSysString->GetTitle(), snn, trg.Data(), mth.Data(),
	   tgt.Data());

    if      (trg.EqualTo("INEL") || trg.EqualTo("MBOR"))      trg = "INEL";
    else if (trg.EqualTo("INEL>0") || trg.EqualTo("INELGt0")) trg = "INELGt0";
    else if (trg.EqualTo("NSD") || trg.EqualTo("V0AND"))      trg = "NSD";
    

    std::ofstream outs("gse.C");
    if (!outs) { 
      Error("Export", "Failed to open output file %s", fname.Data());
      return;
    }
    outs << "// \n"
	 << "TList* gse() {\n"
	 << "  TString mkLib = gSystem->GetMakeSharedLib();\n"
	 << "  mkLib.ReplaceAll(\"-std=c++14\", \"-std=c++98\");\n"
	 << "  gSystem->SetMakeSharedLib(mkLib);\n"      
	 << "  TString gseDir(gSystem->ExpandPathName(\"~/GraphSysErr\"));\n"
	 << "  TString phyDir(gSystem->ExpandPathName(\"${ALICE_PHYSICS}\"));\n"
	 << "  gROOT->SetMacroPath(Form(\"%s:%s/PWGLF/FORWARD/analysis2/dndeta:%s\",gseDir.Data(),phyDir.Data(),gROOT->GetMacroPath()));\n"
	 << "  gSystem->AddIncludePath(Form(\"-I%s\",gseDir.Data()));\n"
	 << "  gROOT->LoadMacro(\"GraphSysErr.C+g\");\n"
	 << "  gROOT->LoadMacro(\"SysErrorAdder.C+g\");\n"
	 << "  TString  t = \"" << trg << "\";\n"
	 << "  TString  c = \"" << mth <<"\";\n"
	 << "  TString  s = \"" << fSysString->GetTitle() << "\";\n"
	 << "  UShort_t e = " << snn << ";\n"
	 << "  TString  m = \"" << tgt << ".C\";\n"
	 << "  TList*   r = new TList;\n"
	 << "  THStack* l = new THStack(\"l\",\"l\");\n"
	 << "  gROOT->Macro(Form(\"%s((THStack*)%p,0,20)\",m.Data(),l));\n"
	 << "  SysErrorAdder* a = SysErrorAdder::Create(t,s,e,c);\n"
	 << "  TIter  n(l->GetHists());\n"
	 << "  TH1*   h = 0;\n"
	 << "  Bool_t f = true;\n"
	 << "  while ((h = static_cast<TH1*>(n()))) {\n"
	 << "    TString nme(h->GetName());\n"
	 << "    if (nme.Contains(\"mirror\")||nme.Contains(\"SysError\"))\n"
	 << "      continue;\n"
	 << "    h->SetMarkerColor(kBlack);\n"
	 << "    h->SetFillColor(kBlack);\n"
	 << "    h->SetLineColor(kBlack);\n"
	 << "    GraphSysErr* g = a->Make(h,0);\n"
	 << "    g->SetSumOption(GraphSysErr::kBox);\n"
	 << "    g->SetSumFillColor(g->GetMarkerColor());\n"
	 << "    g->SetSumFillStyle(3002);\n"
	 << "    g->SetCommonSumOption(GraphSysErr::kBox);\n"
	 << "    g->SetCommonSumFillColor(g->GetMarkerColor());\n"
	 << "    g->SetCommonSumFillStyle(3001);\n"
	 << "    if (f) g->Draw(\"SUM QUAD AXIS\");\n"
	 << "    else   g->Draw(\"SUM QUAD\");\n"
	 << "    f = false;\n"
	 << "    r->Add(g);\n"
	 << "  }\n"
	 << "  TFile* file = TFile::Open(Form(\"%s_%05d_%s%s.root\",s.Data(),e,t.Data(),c.Data()),\"RECREATE\");\n"
	 << "  r->Write(\"container\",TObject::kSingleKey);\n"
	 << "  file->Write();\n"
	 << "  return r;\n"
	 << "}\n"
	 << "// EOF" << std::endl;
    outs.close();
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
    return (!(fOptions & kForceMB) && 
	    fCentAxis && 
	    (fCentAxis->GetNbins() > 1 ||
	     (fCentAxis->GetNbins() == 1 && 
	      fCentAxis->GetXmin() <= fCentAxis->GetXmax()))); 
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
  Double_t     fDelta;        // IP_delta wrt to (0.005,0.184)
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
  TObject*     fCentMeth;     // Centrality method 
  Float_t      fTriggerEff;   // Trigger efficiency 
  Bool_t       fExtTriggerEff;// True if read externally 
  UShort_t     fCentMin;      // Least centrality to plot
  UShort_t     fCentMax;      // Largest centrality to plot
  UInt_t       fCentSeen;     // List of seen centralities
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
  TObject*     fEmpCorr;      // Empirical correction 

  static const Float_t kRightMargin;
  static const Int_t   kFont;
  static const Int_t   kAliceBlue;
  static const Int_t   kAliceRed;
  static const Int_t   kAlicePurple;
  static const Int_t   kAliceYellow;
  static const Int_t   kSysErrColor;
};

const Float_t dNdetaDrawer::kRightMargin = 0.01;
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
 * @deprecated Use new GSE based drawing 
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
 * @param centMin   Least centrality 
 * @param centMax   Largest centrality 
 *
 * @ingroup pwglf_forward_dndeta
 * @deprecated Use new GSE based drawing 
 */
void
DrawdNdeta(const char* filename="forward_dndeta.root", 
	   const char* title="",
	   UShort_t    rebin=5, 
	   UShort_t    others=0x7,
	   UInt_t      flags=dNdetaDrawer::kDefaultOptions,
	   Double_t    meanIpX=-1,
	   Double_t    meanIpY=-1,
	   Float_t     eff=0,
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
  pd->SetEmpirical("file://../empirical.root#default");
  pd->SetDelta(meanIpX, meanIpY);
  // d.fClusterScale = "1.06 -0.003*x +0.0119*x*x";
  pd->Run(filename, title, rebin, others, flags, 0, 0, 0, eff,
	  0, 0, +999, -999, base, outflg);
}

/** 
 * Display usage information
 * 
 */
void
UsageS()
{
  std::ostream& o = std::cout;
  o << "Usage: DrawdNdeta(FILE,TITLE,OTHERS,OPTIONS,FORMATS,REBIN,EFF,"
    << "CMIN,CMAX,IPZMIN,IPZMAX,BASE)\n"
    << "  const char* FILE   File name to open (\"forward_dndeta.root\")\n"
    << "  const char* TITLE  Title to put on plot (\"\")\n"
    << "  const char* OTHERS Other data to draw - more below (\"all\")\n"
    << "  const char* FLAGS  Visualisation flags - more below (\"default\")\n"
    << "  const char* FMT    (optional) Output formats (\"all\")\n"
    << "  UShort_t    REBIN  (optional) Rebinning factor (5)\n"
    << "  Float_t     EFF    (optional) Trigger efficiency\n"
    << "  Float_t     IPZMIN (optional) Least z coordinate of IP\n"
    << "  Float_t     IPZMAX (optional) Largest z coordinate of IP\n"
    << "  UShort_t    CMIN   (optional) Least centrality\n"
    << "  UShort_t    CMAX   (optional) Largest centrality\n"
    << "  const char* BASE   (optional) base name of output files\n"
    << "\n";
  o << " OTHERS space separated list of\n"
    << "  UA5     Show UA5 data (INEL,NSD, ppbar, 900GeV)\n"
    << "  CMS     Show CMS data (NSD, pp)\n"
    << "  ALICE   Show published ALICE data (INEL,INEL>0,NSD, pp)\n"
    << "  WIP     Show event genertor data/work-in-progress\n"
    << "\n";
  o << " FMT space separated list of \n"
    << "  PNG     Make PNG output\n"
    << "  PDF     Make PDF output\n"
    << "  ROOT    Make ROOT file output\n"
    << "  C       Make ROOT script output\n"
    << "\n";
  o << " FLAGS is a bit mask of\n"
    << "  ratio           Show ratios of data to other data and possibly MC\n"
    << "  asymmetry       Show left-right asymmetry\n"
    << "  syserror        Show systematic error band\n"
    << "  rings           Show individual ring results (INEL only)\n"
    << "  noedges         Cut edges when rebinning\n"
    << "  noouters        Remove FMDxO points\n"
    << "  finalmc         Apply `final MC' correction\n"
    << "  empirical[=URL] Apply `Emperical' correction\n"
    << "  mb              Force use of MB\n"
    << "  mirror          Mirror data\n"
    << "  export          Export results to script\n"
    << "  exec            Add code to do combined zooms on eta axis\n"
    << "  old             Assume old-style input\n"
    << "  verbose         Be verbose\n"
    << "  hires           Hi-res batch output\n"
    << "  extrawhite      Add aditional white-space above results\n"
    << "  logo            Add ALICE logo\n"
    << "  nocentral       Do not plot cluster data\n"
    << "  nolabels        No labels on y-axis\n"
    << "\n";
  o << "finalmc requires the file forward_dndetamc.root\n"
    << "empirical requires the histogram at URL \n"
    << std::endl;

}

void 
Draw(const char* filename,
     const char* title="",
     const char* others="ALL",
     const char* options="DEFAULT",
     const char* outFlg="ALL",
     UShort_t    rebin=5,
     Float_t     eff=0, 
     const char* base="")
{
  TString fname(filename);
  fname.ToLower();
  if (fname.CompareTo("help") == 0 || 
      fname.CompareTo("--help") == 0) { 
    UsageS();
    return;
  }

  dNdetaDrawer* pd = new dNdetaDrawer;
  pd->SetEmpirical("file://../empirical.root#default");

  pd->Run(filename, title, others, options, outFlg, rebin, eff,
	  0, 0, +999, -999, base);
}
//____________________________________________________________________
//
// EOF
//

