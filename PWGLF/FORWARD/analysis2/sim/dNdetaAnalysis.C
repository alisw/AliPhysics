#ifndef __CINT__
# include <TH1.h>
# include <TH2.h>
# include <TMath.h>
# include <TParticle.h>
# include <TObjArray.h>
# include <TString.h>
#else
class TH1;
class TParticle;
#endif
#include "FastAnalysis.C"

//====================================================================
/** 
 * Base class for making @f$ 1/N dN_{ch}/d\eta@f$ 
 */
struct dNdetaAnalysis : public FastAnalysis
{
  TH1* fdNdeta; //!
  
  dNdetaAnalysis(Bool_t verbose=false, Int_t monitor=0)
    : FastAnalysis(verbose,monitor), fdNdeta(0)
  {}
  /** 
   * Static member function to create a histogram 
   * 
   * @return Newly created histogram
   */  
  static TH1D* CreatedNdeta()
  {
    Double_t maxEta = 5; // 10;
    Double_t dEta   = 10./200 * 5;
    TH1D* eta = new TH1D("dNdeta", "Charged particle pseudo-rapidity density",
	 		 Int_t(2*maxEta/dEta+.5), -maxEta, +maxEta);
    eta->Sumw2();
    eta->SetXTitle("#it{#eta}");
    eta->SetYTitle("1/N d#it{N}_{ch}/d#it{#eta}");
    eta->SetMarkerColor(kRed+2);
    eta->SetMarkerStyle(20);
    eta->SetDirectory(0);
    
    return eta;
  }
  /** 
   * Called on each slave at start of processing
   */
  virtual void SlaveBegin(TTree*)
  {
    Info("SlaveBegin", "Making dN/deta histogram");
    fdNdeta = CreatedNdeta();
    fdNdeta->SetMarkerColor(kBlack);
    fdNdeta->SetMarkerStyle(21);
    fdNdeta->SetTitle(GetName());
    fOutput->Add(fdNdeta);
  }
  /** 
   * Process a particle.  
   * 
   * @param p Particle to process
   */    
  virtual Bool_t ProcessParticle(const TParticle* p)
  {
    Double_t   pT    = p->Pt();
    Double_t   pZ    = p->Pz();
    Double_t   theta = TMath::ATan2(pT, pZ);
    Double_t   eta   = (pT < 1e-10 ? 1024 : -TMath::Log(TMath::Tan(theta/2)));
    if (TMath::Abs(eta) > 1000) return false;

    Fill(eta);
    return true;
  }
  virtual void Fill(Double_t eta) { fdNdeta->Fill(eta); }
  /** 
   * Final processing.  Scales the histogram to the nubmer of events
   * and the bin width.
   */
  virtual void Terminate()
  {
    fOK = GetEventCount();
    if (fOK <= 0) {
      SetStatus(-1);
      Warning("Terminate", "No events selected");
      return;
    }
    fdNdeta = static_cast<TH1*>(GetOutputObject("dNdeta", TH1::Class()));
    if (!fdNdeta) {
      SetStatus(-1);
      Warning("Terminate", "No dN/deta histogram found");
      return;
    }
    Printf("A total of %ld events", fOK);
    TH1* copy = static_cast<TH1*>(fdNdeta->Clone("before"));
    fOutput->Add(copy);
						
    fdNdeta->Scale(1. / fOK, "width");
    fdNdeta->Draw();
    SetStatus(0);
  }
  /** 
   * Get the list of monitored objects 
   * 
   * @return The list of monitored objects 
   */
  virtual TList* GetMonitorObjects()
  {
    TObject* m1 = new TNamed("dNdeta", "");
    m1->SetUniqueID(0x8); // Scale 

    TList* ret = new TList;
    ret->Add(m1);
    
    return ret;
  }
  /** 
   * Create a new analysis object
   * 
   * @param type The type 
   * @param verbose Whether to be verbose 
   * @param monitor Monitor frequency 
   * 
   * @return newly created object, or null
   */
  static FastAnalysis* Make(const char* type,
			    Bool_t verbose,
			    Int_t monitor);
  ClassDef(dNdetaAnalysis,1);
};

//====================================================================
/** 
 * Select NSD events and builds the @f$ 1/N dN_{ch}/d\eta@f$ 
 * 
 */
struct NSDAnalysis : public dNdetaAnalysis
{
  /** 
   * Constructor 
   * 
   * @param filename File to open 
   * @param verbose  Whether to be verbose 
   */
  NSDAnalysis(Bool_t verbose=false, Int_t monitor=0)
    : dNdetaAnalysis(verbose, monitor)
  {}
  /** 
   * Process the header.  Return 
   * 
   * @return True in case the event is flagged as NSD, false
   * otherwise.
   */
  virtual Bool_t ProcessHeader()
  {
    if (fHeader->fType & 0x1) return false;
    return true;
  }
  ClassDef(NSDAnalysis,1);
};

//====================================================================
/**
 * Processes INEL events and build the @f$ 1/N dN_{ch}/d\eta@f$ 
 * 
 */
struct INELAnalysis : public dNdetaAnalysis
{
  /** 
   * Constructor 
   * 
   * @param verbose  Whether to be verbose 
   */  
  INELAnalysis(Bool_t verbose=false, Int_t monitor=0)
    : dNdetaAnalysis(verbose,monitor)
  {}
  /** 
   * Process the header.  
   * 
   * @return Always true - by definition INEL is all events. 
   */
  virtual Bool_t ProcessHeader() { return true; }
  ClassDef(INELAnalysis,1);
};

//====================================================================
/**
 * Processes INEL events and build the @f$ 1/N dN_{ch}/d\eta@f$ 
 * 
 */
struct INELGt0Analysis : public dNdetaAnalysis
{
  /** 
   * Constructor 
   * 
   * @param verbose  Whether to be verbose 
   */  
  INELGt0Analysis(Bool_t verbose=false, Int_t monitor=0)
    : dNdetaAnalysis(verbose,monitor)
  {}
  /** 
   * Process the header.  
   * 
   * @return Always true - by definition INEL is all events. 
   */
  virtual Bool_t ProcessHeader()
  {
    if (fHeader->fType & 0x4) return true;
    return false;
  }
  ClassDef(INELGt0Analysis,1);
};

//====================================================================
struct CentAnalysis : public dNdetaAnalysis
{
  TAxis*       fCentAxis;
  /// THStack*     fCentStack; //!
  TList*       fCentList;
  Int_t        fCentBin;
  const Long_t fMinEvents;
  TH1*         fCentAll; //!
  TH1*         fCentAcc; //!
  TH2*         fCentMult; //!
  TH2*         fCentNPart; //!
  TH2*         fCentNBin; //!
  TH2*         fCentB; //!
  TH1*         fBtoC; //! 
  /** 
   * Constructor 
   * 
   * @param method   Centrality method 
   * @param verbose  Verbosity flag 
   */
  CentAnalysis(const char* method="V0M", Bool_t verbose=true, Int_t monitor=0)
    : dNdetaAnalysis(verbose,monitor),
      fCentAxis(0),
      fCentList(0),
      fCentBin(0),
      fMinEvents(100),
      fCentAll(0),
      fCentAcc(0),
      fCentMult(0),
      fCentNPart(0),
      fCentNBin(0),
      fCentB(0),
      fBtoC(0)
  {
    TString    axis("default");
    TString    meth(method);
    TObjArray* tokens = meth.Tokenize(":");
    
    fCentMethod = tokens->At(0)->GetName();
    if (tokens->GetEntriesFast() > 1)
      axis = tokens->At(1)->GetName();
    
    if (fCentMethod.Contains("RefMult")) SetCentAxis("mult");
    else                                 SetCentAxis(axis);
  }
  /** 
   * Set the centrality axis (equi-distant)
   * 
   * @param n     Number of bin
   * @param low   Lowest bound 
   * @param high  Highest bound 
   */
  void SetCentAxis(Int_t n, Double_t low, Double_t high)
  {
    if (fCentAxis) {
      delete fCentAxis;
      fCentAxis = 0;
    }
    fCentAxis = new TAxis(n, low, high);
  }
  /** 
   * Set the centrality axis.
   * 
   * @param n       Number of bins 
   * @param edges   Bin edges 
   */
  void SetCentAxis(Int_t n, Double_t* edges)
  {
    if (fCentAxis) {
      delete fCentAxis;
      fCentAxis = 0;
    }
    fCentAxis = new TAxis(n, edges);
  }
  /** 
   * Set the centrality axis.  Spec is either a pre-defined string, or
   * a list of colon separated bin edges. Pre-defined settings are 
   *
   * - pbpb, aa, default: For PbPb
   * - ppb, pbp, pa, ap: For pPb/Pbp 
   * - pp: for pp centrality 
   * 
   * Example of specs could be 
   *
   * - 0:5:10:20:30:40:50:60:80:100
   * - 0:0.1:1:5:10:20:40:70:100 
   * 
   * @param spec 
   */
  void SetCentAxis(const char* spec)
  {
    Info("SetCentAxis", "Setting centrality axis from %s", spec);
    TString    s(spec);
    s.ToLower();
    if (s.IsNull()) return;
    if (s.EqualTo("pbpb") || s.EqualTo("aa") || s.EqualTo("default")) {
      Printf("Setting centrality axis Pb-Pb");
      Double_t aa[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
      SetCentAxis(11, aa);
      fCentAxis->SetUniqueID(2);
      return;
    }
    if (s.EqualTo("ppb") || s.EqualTo("pbp") ||
	s.EqualTo("pa") || s.EqualTo("ap")) {
      Printf("Setting centrality axis p-Pb/Pb-p");
      Double_t pa[] = { 0, 5, 10, 20, 40, 60, 80, 100 };
      SetCentAxis(7, pa);
      fCentAxis->SetUniqueID(3);
      return;
    }
    if (s.EqualTo("pp")) {
      Printf("Setting centrality axis pp");
      Double_t pp[] = { 0, 0.01, 0.1, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100 };
      SetCentAxis(12, pp);
      fCentAxis->SetUniqueID(1);
      return;
    }
    
    TObjArray* tokens  = s.Tokenize(":");
    Int_t      nTokens = tokens->GetEntriesFast();
    TArrayD    edges(nTokens);
    for (Int_t i = 0; i < nTokens; i++) {
      TObjString* token = static_cast<TObjString*>(tokens->At(i));
      TString&    edge  = token->String();
      edges[i]          = edge.Atof();
    }
    SetCentAxis(edges.GetSize()-1, edges.GetArray());
    delete tokens;
  }
  /** 
   * Get the color associated with a centrality bin. 
   * 
   * @param low  Low edge 
   * @param high High edge 
   * 
   * @return Color identifier. 
   */
  Int_t GetCentralityColor(Double_t, Double_t high) const
  {
    Float_t  fc       = high / 100;
    Int_t    nCol     = gStyle->GetNumberOfColors();
    Int_t    icol     = TMath::Min(nCol-1,int(fc * nCol + .5));
    Int_t    col      = gStyle->GetColorPalette(icol);
    //Info("GetCentralityColor","%3d: %3d-%3d -> %3d",bin,centLow,centHigh,col);
    return col;
  }
  /** 
   * Get the histogram name 
   * 
   * @param low  Low edge 
   * @param high High edge
   * 
   * @return String
   */
  virtual const char* HistName(Double_t low, Double_t high) const
  {
    return Form("h%03dd%02d_%03dd%02d",
		Int_t(low),  Int_t(100*low) %100,
		Int_t(high), Int_t(100*high)%100);
  }
  /** 
   * Get the histogram title
   * 
   * @param low  Low edge 
   * @param high High edge
   * 
   * @return String 
   */
  virtual const char* HistTitle(Double_t low, Double_t high) const
  {
    return Form("%6.2f-%6.2f%%", low, high);
  }
  /** 
   * Modify a histogram.  Sets the line, marker, and fill styles and
   * colors, as well as the name and title. 
   * 
   * @param h    Histogram to modify 
   * @param low  Low edge of our bin
   * @param up   High edge of our bin 
   */
  void ModHist(TH1* h, Double_t low, Double_t high)
  {
    Color_t col = GetCentralityColor(low, high);
    h->SetLineColor(col);
    h->SetLineStyle(1);
    h->SetMarkerColor(col);
    h->SetMarkerStyle(24);
    h->SetFillColor(kWhite);
    h->SetFillStyle(0);
    h->SetName(HistName(low, high));
    h->SetTitle(HistTitle(low, high));
  }
  /** 
   * Called on each slave at start of processing
   * 
   * @param t Ignored
   */
  virtual void SlaveBegin(TTree* t)
  {
    dNdetaAnalysis::SlaveBegin(t);
    // We use the parent fdNdeta histogram as a cache 
    fdNdeta->SetName("cache");
    fdNdeta->SetTitle("Cache");

    Info("SlaveBegin", "Making stack of dN/deta histograms");
    
    // Create our stack 
    // fCentStack = new THStack("stack", "Stack of all dN_{ch}/d#eta");
    fCentList  = new TList;
    fCentList->SetName("byCent");
    for (Int_t i = 1; i <= fCentAxis->GetNbins(); i++) {
      Double_t low  = fCentAxis->GetBinLowEdge(i);
      Double_t high = fCentAxis->GetBinUpEdge(i);
      TH1*     hist = CreatedNdeta();
      ModHist(hist, low, high);
      hist->SetMinimum(0);
      hist->Sumw2();
      fCentList->Add(hist);
    }
    fOutput->Add(fCentList);

    if (fCentAxis->GetXbins()->GetArray()) 
      fCentAll = new TH1D("cent", "All Centralities",
			   fCentAxis->GetNbins(),
			   fCentAxis->GetXbins()->GetArray());
    else
      fCentAll = new TH1D("cent", "All Centralities",
			   fCentAxis->GetNbins(),
			   fCentAxis->GetXmin(),
			   fCentAxis->GetXmax());
    fCentAll->SetXTitle("Centrality [%]");
    fCentAll->SetYTitle("Events");
    fCentAll->SetFillColor(kRed+2);
    fCentAll->SetFillStyle(3002);
    fCentAll->SetMarkerStyle(20);
    fCentAll->SetMarkerColor(kRed+2);
    fCentAll->SetLineColor(kRed+2);
    fCentAll->SetDirectory(0);
    fCentAll->SetMinimum(0);
    fOutput->Add(fCentAll);

    fCentAcc = static_cast<TH1*>(fCentAll->Clone("centAcc"));
    fCentAcc->SetTitle("Accepted centralities");
    fCentAcc->SetFillColor(kGreen+2);
    fCentAcc->SetMarkerColor(kGreen+2);
    fCentAcc->SetLineColor(kGreen+2);
    fCentAcc->SetDirectory(0);
    fOutput->Add(fCentAcc);


    fOutput->ls();
  }
  virtual Bool_t SetupEstimator()
  {
    if (!dNdetaAnalysis::SetupEstimator()) return false;
    
    if (fCentMult) return true;

    fCentMult  = 0;
    fCentNPart = 0;
    fCentNBin  = 0;
    Int_t maxNPart = 2*210;
    Int_t maxNBin  = 3*210;
    if (fCentAxis->GetXbins()->GetArray()) {
      fCentNPart = new TH2D("centNPart", "Centrality vs. N_{part}",
			    fCentAxis->GetNbins(), 
			    fCentAxis->GetXbins()->GetArray(),
			    maxNPart, 0, maxNPart);
      fCentNBin  = new TH2D("centNBin", "Centrality vs. N_{bin}", 
			    fCentAxis->GetNbins(), 
			    fCentAxis->GetXbins()->GetArray(),
			    maxNBin, 0, maxNBin);
      fCentB     = new TH2D("centB", "Centrality vs. b",
			    fCentAxis->GetNbins(), 
			    fCentAxis->GetXbins()->GetArray(),
			    200, 0, 20);
    } else {
      fCentNPart = new TH2D("centNPart", "Centrality vs. N_{part}",
			    fCentAxis->GetNbins(), fCentAxis->GetXmin(), 
			    fCentAxis->GetXmax(),
			    maxNPart, 0, maxNPart);
      fCentNBin  = new TH2D("centNBin", "Centrality vs. N_{bin}", 
			    fCentAxis->GetNbins(), fCentAxis->GetXmin(), 
			    fCentAxis->GetXmax(),
			    maxNBin, 0, maxNBin);
      fCentB     = new TH2D("centB", "Centrality vs. b",
			    fCentAxis->GetNbins(), fCentAxis->GetXmin(), 
			    fCentAxis->GetXmax(),
			    200, 0, 20);
    }

    fCentNPart->SetXTitle(Form("Centrality (%s) %%", fCentMethod.Data()));
    fCentNPart->SetYTitle("N_{part}");
    fCentNPart->SetDirectory(0);

    fCentNBin->SetXTitle(Form("Centrality (%s) %%", fCentMethod.Data()));
    fCentNBin->SetYTitle("N_{bin}");
    fCentNBin->SetDirectory(0);

    fCentB->SetXTitle(Form("Centrality (%s) %%", fCentMethod.Data()));
    fCentB->SetYTitle("b [fm]");
    fCentB->SetDirectory(0);
    
    fOutput->Add(fCentNPart);
    fOutput->Add(fCentNBin);
    fOutput->Add(fCentB);

    if (fCentAxis->GetUniqueID() == 2) {
      Double_t bin[25];
      Double_t cent[24];
      Int_t    i = 0;
      bin[i] =  0.00; cent[i] = Double_t( 0+  1)/2; i++; //  0 -   1 ( 0.00 -  1.57)
      bin[i] =  1.57; cent[i] = Double_t( 1+  2)/2; i++; //  1 -   2 ( 1.57 -  2.22)
      bin[i] =  2.22; cent[i] = Double_t( 2+  3)/2; i++; //  2 -   3 ( 2.22 -  2.71)
      bin[i] =  2.71; cent[i] = Double_t( 3+  4)/2; i++; //  3 -   4 ( 2.71 -  3.13)
      bin[i] =  3.13; cent[i] = Double_t( 4+  5)/2; i++; //  4 -   5 ( 3.13 -  3.50)
      bin[i] =  3.50; cent[i] = Double_t( 5+ 10)/2; i++; //  5 -  10 ( 3.50 -  4.94)
      bin[i] =  4.94; cent[i] = Double_t(10+ 15)/2; i++; // 10 -  15 ( 4.94 -  6.05)
      bin[i] =  6.05; cent[i] = Double_t(15+ 20)/2; i++; // 15 -  20 ( 6.05 -  6.98)
      bin[i] =  6.98; cent[i] = Double_t(20+ 25)/2; i++; // 20 -  25 ( 6.98 -  7.81)
      bin[i] =  7.81; cent[i] = Double_t(25+ 30)/2; i++; // 25 -  30 ( 7.81 -  8.55)
      bin[i] =  8.55; cent[i] = Double_t(30+ 35)/2; i++; // 30 -  35 ( 8.55 -  9.23)
      bin[i] =  9.23; cent[i] = Double_t(35+ 40)/2; i++; // 35 -  40 ( 9.23 -  9.88)
      bin[i] =  9.88; cent[i] = Double_t(40+ 45)/2; i++; // 40 -  45 ( 9.88 - 10.47)
      bin[i] = 10.47; cent[i] = Double_t(45+ 50)/2; i++; // 45 -  50 (10.47 - 11.04)
      bin[i] = 11.04; cent[i] = Double_t(50+ 55)/2; i++; // 50 -  55 (11.04 - 11.58)
      bin[i] = 11.58; cent[i] = Double_t(55+ 60)/2; i++; // 55 -  60 (11.58 - 12.09)
      bin[i] = 12.09; cent[i] = Double_t(60+ 65)/2; i++; // 60 -  65 (12.09 - 12.58)
      bin[i] = 12.58; cent[i] = Double_t(65+ 70)/2; i++; // 65 -  70 (12.58 - 13.05)
      bin[i] = 13.05; cent[i] = Double_t(70+ 75)/2; i++; // 70 -  75 (13.05 - 13.52)
      bin[i] = 13.52; cent[i] = Double_t(75+ 80)/2; i++; // 75 -  80 (13.52 - 13.97)
      bin[i] = 13.97; cent[i] = Double_t(80+ 85)/2; i++; // 80 -  85 (13.97 - 14.43)
      bin[i] = 14.43; cent[i] = Double_t(85+ 90)/2; i++; // 85 -  90 (14.43 - 14.96)
      bin[i] = 14.96; cent[i] = Double_t(90+ 95)/2; i++; // 90 -  95 (14.96 - 15.67)
      bin[i] = 15.67; cent[i] = Double_t(95+100)/2; i++; // 95 - 100 (15.67 - 20.00)
      bin[i] = 20;
      fBtoC  = new TH1D("bToC", "Centrality from b", 24, bin);
      fBtoC->SetXTitle("b [fm]");
      fBtoC->SetYTitle("Centrality [%]");
      for (Int_t i = 0; i < 24; i++) fBtoC->SetBinContent(i+1, cent[i]);
      fOutput->Add(fBtoC);
    }
    else if (fCentAxis->GetUniqueID()) {
#if 0
      const Double_t cl[] = {   0,    5,   10,   20,   40,   60,   80, -1 };
      const Double_t ch[] = {   5,   10,   20,   40,   60,   80,  100, -1 };
      const Double_t bm[] = {3.12, 3.50, 3.85, 4.54, 5.57, 6.63, 7.51, 30 };
      const Double_t bs[] = {1.39, 1.48, 1.57, 1.69, 1.69, 1.45, 1.11, -1 };
      Double_t       bl[] = {   0,   -1,   -1,   -1,   -1,   -1,   -1, -1, -1 };
      for (Int_t i = 1; i <= 7; i++) bl[i] = bm[i-1] + (bm[i]-bm[i-1])/2;
      fBtoC  = new TH1D("bToC", "Centrality from b", 7, bl);
      fBtoC->SetXTitle("b [fm]");
      fBtoC->SetYTitle("Centrality [%]");
      for (Int_t i = 0; i < 7; i++) fBtoC->SetBinContent(i+1, (ch[i]+cl[i])/2);
#else 
      const Double_t cl[] = {   0,    5,   10,   20,   40,   60,   80, -1 };
      const Double_t ch[] = {   5,   10,   20,   40,   60,   80,  100, -1 };
      const Double_t bm[] = { 0, 1.83675, 2.59375,  3.66875, 5.18625, 6.35475,
			      7.40225, 13.8577, -1};
      fBtoC  = new TH1D("bToC", "Centrality from b", 7, bm);
      fBtoC->SetXTitle("b [fm]");
      fBtoC->SetYTitle("Centrality [%]");
      for (Int_t i = 1; i <= 7; i++) {
	fBtoC->SetBinContent(i, (ch[i-1]+cl[i-1]) / 2);
      }      
#endif
      fOutput->Add(fBtoC);
    }
    
    if (!fCentHist) return true;
    if (fCentAxis->GetXbins()->GetArray()) {
      fCentMult  = new TH2D("centMult","Event multiplicity vs. centrality",
			    fCentHist->GetXaxis()->GetNbins(),
			    fCentHist->GetXaxis()->GetXmin(),
			    fCentHist->GetXaxis()->GetXmax(),
			    fCentAxis->GetNbins(),
			    fCentAxis->GetXbins()->GetArray());
    } else {
      fCentMult  = new TH2D("centMult","Event multiplicity vs. centrality",
			    fCentHist->GetXaxis()->GetNbins(),
			    fCentHist->GetXaxis()->GetXmin(),
			    fCentHist->GetXaxis()->GetXmax(),
			    fCentAxis->GetNbins(),
			    fCentAxis->GetXmin(), fCentAxis->GetXmax());
    }
    fCentMult->SetXTitle(Form("Event multiplicity (%s)", fCentMethod.Data()));
    fCentMult->SetYTitle(Form("Centrality (%s) %%", fCentMethod.Data()));
    fCentMult->SetDirectory(0);
    
    fOutput->Add(fCentMult);
    return true;
  }    
		     
  /** 
   * Clear our internal caches
   */
  virtual void Clear(Option_t* option="") 
  {
    dNdetaAnalysis::Clear(option);
    fdNdeta->Reset(); // Clear the cache
  }
  /** 
   * Process the header.  Accepts events in range 
   * 
   * @return true 
   */
  virtual Bool_t ProcessHeader()
  {
    Double_t cent = 0;
    if (fCentMethod.EqualTo("B"))
      cent = fBtoC->GetBinContent(fBtoC->FindBin(fHeader->fB));
    else
      cent = GetCentrality();
    fCentAll->Fill(cent);
    if (fCentMult) fCentMult->Fill(fEventMult, cent);
    fCentB->Fill(cent, fHeader->fB);
    if (cent < 0 || cent > 999) {
      Warning("ProcessHeader",
	      "Centrality is unreasonable: %f -> %f",fEventMult, cent);
      return false;
    }
    Int_t nBin = fCentAxis->GetNbins();
    fCentBin = fCentAxis->FindBin(cent);
    if (fCentBin-1 == nBin && cent == fCentAxis->GetXmax())
      fCentBin = nBin;
    if (fCentBin < 1 || fCentBin > nBin) {
      Warning("ProcessHeader", "Centrality %f -> %f -> bin # %d",
	      fEventMult, cent, fCentBin);
      fCentBin = -1;
      return false;
    }
    fCentNPart->Fill(cent, fHeader->fNtgt+fHeader->fNproj);
    fCentNBin->Fill(cent, fHeader->fNbin);
    if (fCentBin == nBin) cent -= 0.001;
    fCentAcc->Fill(cent);
    return true;
  }
  /** 
   * Process a single event.  
   *
   * First we fill the internal cache using the base class methods.
   * Then we find the histogram for this particular reference
   * multiplicity and add our event cache to that bin.  The number of
   * events in each bin is counted in the unique ID of each bin.
   */
  virtual void ProcessParticles()
  {
    // Check we got a bin 
    if (fCentBin < 0) return;

    // Find the histogram to update 
    TH1*     out  = static_cast<TH1*>(fCentList->/*GetHists()->*/
				      At(fCentBin-1));
    // If we still have no histogram, return immediately 
    if (!out) return;

    // Use parent function to fill cache 
    dNdetaAnalysis::ProcessParticles();
    if (fVerbose) {
      Double_t n0   = fdNdeta->GetBinContent(fdNdeta->GetNbinsX()/2);
      Double_t d0   = fdNdeta->GetXaxis()->GetBinWidth(fdNdeta->GetNbinsX()/2);
      Double_t eta0 = (n0 / d0);
      Printf("Centrality %6.2f-%6.2f (bin %4d) "
	     "Nch=%8.1f/%4.2f=%8.1f out=%s (%d)",
	     fCentAxis->GetBinLowEdge(fCentBin),
	     fCentAxis->GetBinUpEdge(fCentBin),
	     fCentBin,
	     n0, d0, eta0,
	     out->GetName(),
	     Int_t(out->GetBinContent(0)));
    }

    // Make sure we have nothing in the underflow bin. 
    fdNdeta->SetBinContent(0,0);

    // Add our cache to the appropriate bin 
    out->Add(fdNdeta);
    // Increment underflow bin for the event count 
    out->SetBinContent(0, out->GetBinContent(0)+1);
  }
  /** 
   * Final processing. 
   * 
   * Normalize each bin to the number of events in each bin (stored in
   * the unique ID of that bin).
   */
  virtual void Terminate()
  {
    fOK = GetEventCount();
    if (fOK <= 0) {
      SetStatus(-1);
      Warning("Terminate", "No events selected");
      return;
    }
    
    fCentAll   = static_cast<TH1*>(GetOutputObject("cent",      TH1::Class()));
    fCentAcc   = static_cast<TH1*>(GetOutputObject("centAcc",   TH1::Class()));
    fCentMult  = static_cast<TH2*>(GetOutputObject("centMult",  TH2::Class()));
    fCentNPart = static_cast<TH2*>(GetOutputObject("centNPart", TH2::Class()));
    fCentNBin  = static_cast<TH2*>(GetOutputObject("centNBin",  TH2::Class()));
	 
    // fCentStack = static_cast<THStack*>(GetOutputObject("stack",
    // THStack::Class()));
    fCentList = static_cast<TList*>(GetOutputObject("byCent",
						    TList::Class()));
    if (!fCentList || !fCentAll || !fCentAcc) {
      Warning("Terminate", "Missing stack and histograms");
      SetStatus(-1);
      return;
    }
    // fCentAll->Scale(1./fOK, "width");
    // fCentAcc->Scale(1./fOK, "width");
    Info("Terminate", "Accepted %d/%d=%6.2f%% events",
	 int(fCentAcc->GetEntries()), int(fCentAll->GetEntries()),
	 100*fCentAcc->GetEntries()/fCentAll->GetEntries());

    THStack*  stack = new THStack("all", "All");
    TList*    hists = fCentList; // ->GetHists();
    TObjLink* link  = hists->FirstLink();
    Int_t     bin   = 1;
    Long64_t  sum   = 0;
    Long64_t  total = 0;
    Long64_t  all   = fCentAll->GetEntries();
    while (link) {
      TObject* o = link->GetObject();
      if (!o) {
	link = link->Next();
	bin++;
	continue;
      }
      TH1*  h = static_cast<TH1*>(o);
      Int_t n = h->GetBinContent(0);
      Int_t m = fCentAcc->GetBinContent(bin);
      total += m;
      h->SetBinContent(0,0);
      printf("%9d (%9d) events in bin %s ...", n, m, o->GetTitle());
      if (n < fMinEvents) {
	// Too few event, remove this
	TObjLink* tmp = link->Next();
	hists->Remove(link);
	link = tmp;
	delete o;
	Printf(" removed");
	bin++;
	continue;
      }
      sum += m;
      // Scale
      h->Scale(1. / n, "width");
      stack->Add(h);
      Printf(" scaled");
      link = link->Next();
      bin++;
    }
    Printf("ana/acc/all: %9lld/%9lld/%9lld [%6.2f%%/%6.2f%%]",
	   sum, total, all, float(100*sum)/total, float(100*total)/all);
    fOutput->Add(stack);
  }
  /** 
   * Get the list of monitored objects 
   * 
   * @return The list of monitored objects 
   */
  virtual TList* GetMonitorObjects()
  {
    TObject* m1 = new TNamed("cent",     "hist text30");
    TObject* m2 = new TNamed("centAcc",  "hist text30");
    TObject* m3 = new TNamed("byCent",   "e");
    
    m3->SetUniqueID(0x8); // Scale 
    TList* ret = new TList;
    ret->Add(m1);
    ret->Add(m2);
    ret->Add(m3);
    
    return ret;
  }
  ClassDef(CentAnalysis,1);
};
  
//====================================================================
/**
 * Processes events and build the @f$ 1/N dN_{ch}/d\eta@f$ for each
 * bin in reference multiplicity.  The reference multiplicity is the
 * number of charged particles with @f$|\eta|\le0.8@f$
 */
struct MultAnalysis : public CentAnalysis
{
  /** 
   * Constructor. 
   * 
   * @param filename File to open
   * @param verbose  Whether to verbose
   */
  MultAnalysis(const char* method="RefMult00d80",
	       Bool_t verbose=false,
	       Int_t monitor=0)
    : CentAnalysis(method, verbose,monitor)
  {
    //              +1 +2 +3 +3  +5, +5, +5, +5,+10,+10,+10,+10,+10,+10,+10
    Double_t bins[]={ 0, 3, 6, 9, 14, 19, 24, 29, 39, 49, 59, 69, 79, 89, 99 };
    CentAnalysis::SetCentAxis(14, bins);
  }
  /** 
   * Get the histogram title
   * 
   * @param low  Low edge 
   * @param high High edge
   * 
   * @return String 
   */
  virtual const char* HistTitle(Double_t low, Double_t high) const
  {
    Int_t iLow  = low;
    Int_t iHigh = high;
    if (iLow == iHigh) {
      if (iLow == 0) return " 0";
      else           return Form("%2d+", iLow);
    }
    return Form("%2d-%2d", iLow+1, iHigh);
  }
  /** 
   * Called on each slave at start of processing
   * 
   * @param t Ignored
   */
  virtual void SlaveBegin(TTree* t)
  {
    CentAnalysis::SlaveBegin(t);

    // Create null-bin 
    TH1* first = CreatedNdeta();
    ModHist(first, 0, 0);
    fCentList->/*GetHists()->*/AddFirst(first);

    // Create overflow-bin
    TH1* last = CreatedNdeta();
    ModHist(last, 100, 100);
    fCentList->/*GetHists()->*/AddLast(last);
  }
  /** 
   * Process the header.  Accepts events in range 
   * 
   * @return true 
   */
  virtual Bool_t ProcessHeader()
  {
    if (fEventMult < 0) return false;
    fCentBin = fCentAxis->FindBin(Int_t(fEventMult)-.1)+1;
    Printf("Event multiplicity: %d -> bin %d", Int_t(fEventMult), fCentBin);
    return true;
  }
  ClassDef(MultAnalysis,1);
};

//====================================================================
/*
 * The function to make our analyser 
 */
FastAnalysis*
dNdetaAnalysis::Make(const char* type,
		     Bool_t verbose,
		     Int_t monitor)
{
  TString t(type);
  if      (t.EqualTo("INEL"))    return new INELAnalysis(verbose,monitor);
  else if (t.EqualTo("NSD"))     return new NSDAnalysis(verbose,monitor);
  else if (t.EqualTo("INELGt0")) return new INELGt0Analysis(verbose,monitor);
  else if (t.BeginsWith("MULT") || t.BeginsWith("CENT")) {
    TString w(t(4, t.Length()-4));
    if (!(w.BeginsWith("RefMult") ||
	  w.BeginsWith("ZNA") || 
	  w.BeginsWith("ZNC") || 
	  w.BeginsWith("ZPA") || 
	  w.BeginsWith("ZPC") || 
	  w.BeginsWith("V0M") ||
	  w.BeginsWith("V0A") ||
	  w.BeginsWith("V0C") ||
	  w.BeginsWith("B"))) {
      Printf("Warning: dNdetaAnalysis::Make: Unknown estimator: %s", w.Data());
      return 0;
    }
    if (t.BeginsWith("MULT"))
      return new MultAnalysis(w, verbose, monitor);
    else
      return new CentAnalysis(w, verbose, monitor);
  }
  Printf("Error: dNdetaAnalysis::Run: Invalid spec: %s", t.Data());
  return 0;
}
//
// EOF
//
