/**
 * @file   SummaryDrawer.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Sun Nov 25 11:36:41 2012
 * 
 * @brief  Base class for classes to draw summaries 
 * 
 * 
 */
#ifndef SUMMARYDRAWER_C
# define SUMMARYDRAWER_C
# ifndef __CINT__
#  include <THStack.h>
#  include <TH1.h>
#  include <TH2.h>
#  include <TH3.h>
#  include <TParameter.h>
#  include <TCanvas.h>
#  include <TList.h>
#  include <TFile.h>
#  include <TError.h>
#  include <TLatex.h>
#  include <TLegend.h>
#  include <TLegendEntry.h>
#  include <TMath.h>
#  include <TString.h>
#  include <TStyle.h>
#  include <TSystem.h>
#  include <TProfile.h>
#  include <TGaxis.h>
#  include <TPad.h>
#  include <TRegexp.h>
#  include <TGraph.h>
#  include <sstream>
#  include <iomanip>
# else 
#  ifdef __ACLIC__
class THStack;
class TH1;
class TH2;
class TH3;
class TCollection;
class TCanvas;
class TVirtualPad;
class TPad;
class TLatex;
class TAxis;
#  endif
# endif

/**
 * Base class for summary drawers
 * 
 */
class SummaryDrawer 
{
public:
  enum { 
    kLogx   = 0x1, 
    kLogy   = 0x2, 
    kLogz   = 0x4, 
    kLegend = 0x10, 
    kGridx  = 0x100, 
    kGridy  = 0x200, 
    kGridz  = 0x400,
    kSilent = 0x800,
    kNorth  = 0x1000, 
    kMiddle = 0x2000,
    kSouth  = 0x3000, 
    kEast   = 0x10000, 
    kCenter = 0x20000,
    kWest   = 0x30000
  };
  enum { 
    kLandscape         = 0x100, 
    kPause             = 0x200
  };    
  SummaryDrawer() 
    : fCanvas(0), 
      fTop(0), 
      fBody(0),
      fHeader(0),
      fParName(0),
      fParVal(0),
      fPause(false),
      fLandscape(false), 
      fRingMap(0), 
      fPDF(true),
      fLastTitle("")
  {
    fRingMap = new TVirtualPad*[6];
    fRingMap[0] = 0;
    fRingMap[1] = 0;
    fRingMap[2] = 0;
    fRingMap[3] = 0;
    fRingMap[4] = 0;
    fRingMap[5] = 0;
  }
  virtual ~SummaryDrawer() {}

protected:
  //____________________________________________________________________
  /** 
   * Get null terminated array of ring names 
   * 
   * @param lower If true, return in the form FMD[1-3][io], otherwise 
   *              in the form FMD[1-3][IO]
   * 
   * @return Null terminated array of ring names 
   */
  static const Char_t** GetRingNames(Bool_t lower=false) 
  {
    static const Char_t* lN[]={ "FMD1i", "FMD2i", "FMD2o", "FMD3o", "FMD3i", 0};
    static const Char_t* uN[]={ "FMD1I", "FMD2I", "FMD2O", "FMD3O", "FMD3I", 0};
    return (lower ? lN : uN);
  }
  //____________________________________________________________________
  /** 
   * Get the standard color for a ring  
   *
   * @param d Detector
   * @param r Ring 
   * 
   * @return 
   */
  static Color_t RingColor(UShort_t d, Char_t r)
  { 
    return ((d == 1 ? kRed : (d == 2 ? kGreen : kBlue))
	    + ((r == 'I' || r == 'i') ? 2 : -3));
  }
  //____________________________________________________________________
  TLegend* DrawRingLegend(TVirtualPad* p, UInt_t flags)
  {
    TLegend* l = MakeLegend(p, flags, false);

    for (UShort_t i = 0; i < 5; i++) {
      UShort_t      d = (i+1)/2+1;
      Char_t        r = (i/2 == 1) ? 'o' : 'i';
      TLegendEntry* e = l->AddEntry("dummy", Form("FMD%d%c", d, r), "f");
      e->SetFillColor(RingColor(d, r));
      e->SetFillStyle(1001);
      e->SetLineColor(kBlack);
    }
   
    l->Draw();
    return l;
  }
  //____________________________________________________________________
  static void SysString(UShort_t sys, TString& str)
  {
    str = "?";
    switch (sys) { 
    case 1: str = "pp";   break;
    case 2: str = "PbPb"; break;
    case 3: str = "pPb";  break;
    case 4: str = "Pbp";  break;
    }
  }
  //____________________________________________________________________
  static void SNNString(UShort_t sNN, TString& str)
  {
    str = "?";
    if      (sNN < 1000) str = Form("%dGeV", sNN);
    else if (sNN < 3000) str = Form("%4.2fTeV", 0.001*sNN);
    else                 str = Form("%dTeV", sNN/1000);
  }
  //____________________________________________________________________
  /** 
   * Append an & to a string and the next term.
   * 
   * @param trg  Output string
   * @param what Term
   */
  static void AppendAnd(TString& trg, const TString& what)
  {
    if (!trg.IsNull()) trg.Append(" & ");
    trg.Append(what);
  }
  //____________________________________________________________________
  static void TriggerString(ULong_t trigger, TString& str)
  {
    str = "";
        /** 
     * Bits of the trigger pattern
     */
    enum { 
      /** In-elastic collision */
      kInel        = 0x0001, 
      /** In-elastic collision with at least one SPD tracklet */
      kInelGt0     = 0x0002, 
      /** Non-single diffractive collision */
      kNSD         = 0x0004, 
      /** Empty bunch crossing */
      kEmpty       = 0x0008, 
      /** A-side trigger */
      kA           = 0x0010, 
      /** B(arrel) trigger */
      kB           = 0x0020, 
      /** C-side trigger */
      kC           = 0x0080,  
      /** Empty trigger */
      kE           = 0x0100,
      /** pileup from SPD */
      kPileUp      = 0x0200,    
      /** true NSD from MC */
      kMCNSD       = 0x0400,    
      /** Offline MB triggered */
      kOffline     = 0x0800,
      /** At least one SPD cluster */ 
      kNClusterGt0 = 0x1000,
      /** V0-AND trigger */
      kV0AND       = 0x2000, 
      /** Satellite event */
      kSatellite   = 0x4000
    };
    if ((trigger & kInel)        != 0x0) AppendAnd(str, "INEL");
    if ((trigger & kInelGt0)     != 0x0) AppendAnd(str, "INEL>0");
    if ((trigger & kNSD)         != 0x0) AppendAnd(str, "NSD");
    if ((trigger & kV0AND)       != 0x0) AppendAnd(str, "V0AND");
    if ((trigger & kA)           != 0x0) AppendAnd(str, "A");
    if ((trigger & kB)           != 0x0) AppendAnd(str, "B");
    if ((trigger & kC)           != 0x0) AppendAnd(str, "C");
    if ((trigger & kE)           != 0x0) AppendAnd(str, "E");
    if ((trigger & kMCNSD)       != 0x0) AppendAnd(str, "MCNSD");
    if ((trigger & kNClusterGt0) != 0x0) AppendAnd(str, "NCluster>0");
    if ((trigger & kSatellite)   != 0x0) AppendAnd(str, "Satellite");
  }
    
  //__________________________________________________________________
  /** 
   * Find an object in a collection
   * 
   * @param parent Parent directory
   * @param name   Name of object
   * @param verb   Be verbose 
   * 
   * @return Pointer to object or null 
   */
  static TObject* GetObject(const TObject* parent, 
			    const TString& name, 
			    Bool_t verb=true)
  {
    if (!parent) {
      if (verb) Warning("GetObject", "No parent given");
      return 0;
    }
    if (name.IsNull()) { 
      if (verb) Warning("GetObject", "No name specified");
      return 0;
    }
    TObject* o = 0;
    if (parent->IsA()->InheritsFrom(TCollection::Class())) {
      const TCollection* p = static_cast<const TCollection*>(parent);
      o = p->FindObject(name);
    }
    else if (parent->IsA()->InheritsFrom(TDirectory::Class())) {
      const TDirectory* d = static_cast<const TDirectory*>(parent);
      o = const_cast<TDirectory*>(d)->Get(name);
    }
    else 
      Warning("GetObject", "Do not know how to find an object (%s) in "
	      "%s (of class %s)", name.Data(), 
	      parent ? parent->GetName() : "?", 
	      parent ? parent->ClassName() : "?");
    if (!o) {
      if (verb) Warning("GetObject", "Object \"%s\" not found in parent \"%s\"",
			name.Data(), parent->GetName());
      return 0;
    }
    return o;    
  }
  //____________________________________________________________________
  /** 
   * Check the type of a found object 
   * 
   * @param o   Object 
   * @param cl  Class 
   * @param src Source of object
   * 
   * @return true on success, false otherwise 
   */
  static Bool_t CheckType(const TObject* o, 
			  const TClass*  cl, 
			  const TString& src)
  {
    if (!o->IsA()->InheritsFrom(cl)) { 
      Warning("CheckType", "Object \"%s\" retrieved from \"%s\" is not a "
	      "%s but a %s", o->GetName(), src.Data(), cl->GetName(), 
	      o->ClassName());
      return false;
    }
    return true;
  }
  //__________________________________________________________________
  /** 
   * Check a possibly returned object. 
   * 
   * @param o  Object found - if any
   * @param p  Parent of object
   * 
   * @return A pointer to the object cast to the right type
   */
  template <typename T>
  static T* DoGetObject(TObject* o, const TObject* p) 
  {
    if (!o) return 0;
    if (!CheckType(o, T::Class(), p->GetName())) return 0;
    return static_cast<T*>(o);
  }
  //__________________________________________________________________
  /** 
   * Check a returned parameter from a parent
   * 
   * @param o     Possibly found object
   * @param p     Parent object
   * @param value Value 
   * 
   * @return true on success, false otherwise 
   */
  template <typename T>
  static Bool_t DoGetParameter(TObject* o, const TObject* p, T& value) 
  {
    TParameter<T>* r = DoGetObject<TParameter<T> >(o, p);
    if (!r) return false;
    // if (r->TestBit(TParameter<T>::kFirst)) value = r->GetVal();
    // else                                   value = r->GetUniqueID();
    value = r->GetVal();
    if (!r->TestBit(BIT(19))) {
      TObject* oc = GetObject(p, "count", false);
      if (oc) {
	TParameter<int>* pc = static_cast<TParameter<int>*>(oc);
	int cnt = pc->GetVal();
	value /= cnt;
      }
      else 
	value = r->GetUniqueID();
    }
    // value = r->GetUniqueID();
    return true;
  }
    
  //___________________________________________________________________
  /** 
   * Get a Short_t parameter value 
   * 
   * @param c      Parent collection
   * @param name   Name of parameter
   * @param value  On return the value
   * @param verb   If true, complain if not found 
   */
  static Bool_t GetParameter(const TObject*  c, 
			     const TString&  name, 
			     Short_t&        value,
			     Bool_t          verb=true)
  {
    int v;
    Bool_t r = DoGetParameter(GetObject(c, name, verb), c, v); 
    value = v;
    return r;
  }
  //___________________________________________________________________
  /** 
   * Get a UShort_t parameter value 
   * 
   * @param c      Parent collection
   * @param name   Name of parameter
   * @param value  On return the value
   * @param verb   If true, complain if not found 
   */
  static Bool_t GetParameter(const TObject*  c, 
			     const TString&  name, 
			     UShort_t&       value,
			     Bool_t          verb=true)
  {
    int v = 0;
    Bool_t r = DoGetParameter(GetObject(c, name, verb), c, v); 
    value = v;
    return r;
  }
  //___________________________________________________________________
  /** 
   * Get a ULong_t parameter value 
   * 
   * @param c      Parent collection
   * @param name   Name of parameter
   * @param value  On return the value
   * @param verb   If true, complain if not found 
   */
  static Bool_t GetParameter(const TObject*  c, 
			     const TString&  name, 
			     ULong_t&        value,
			     Bool_t          verb=true)
  {
    Long_t v = 0;
    Bool_t r = DoGetParameter(GetObject(c, name, verb), c, v); 
    value = v;
    return r;
  }
  //_____________________________________________________________________
  /** 
   * Get a Int_t parameter value 
   * 
   * @param c      Parent collection
   * @param name   Name of parameter
   * @param value  On return the value
   * @param verb   If true, complain if not found 
   */
  static Bool_t GetParameter(const TObject*  c, 
			     const TString&  name, 
			     Int_t&          value,
			     Bool_t          verb=true)
  {
    return DoGetParameter(GetObject(c, name, verb), c, value);
  }
  //_____________________________________________________________________
  /** 
   * Get a Double_t parameter value 
   * 
   * @param c      Parent collection
   * @param name   Name of parameter
   * @param value  On return the value
   * @param verb   If true, complain if not found 
   */
  static Bool_t GetParameter(const TObject*      c, 
			     const TString&      name, 
			     Double_t&           value,
			     Bool_t              verb=true);
  //_____________________________________________________________________
  /** 
   * Get a Bool_t parameter value 
   * 
   * @param c      Parent collection
   * @param name   Name of parameter
   * @param value  On return the value
   * @param verb   If true, complain if not found 
   */
  static Bool_t GetParameter(const TObject*  c, 
			     const TString&      name, 
			     Bool_t&             value,
			     Bool_t              verb=true)
  {
    return DoGetParameter(GetObject(c, name, verb), c, value);
  }
  //____________________________________________________________________
  /** 
   * Find a collection in another collection 
   * 
   * @param parent Parent collection 
   * @param name   Name of the collection 
   * @param verb   If true and not found, complain
   *
   * @return pointer to collection on success, otherwise null 
   */
  static TCollection* GetCollection(const TObject*     parent, 
				    const TString&     name,
				    Bool_t             verb=true)
  {
    return DoGetObject<TCollection>(GetObject(parent, name, verb), parent);
  }
  //____________________________________________________________________
  /** 
   * Check a 1D histogram object from a parent
   * 
   * @param parent Parent collection 
   * @param name   Name of histogram
   * @param verb   Possibly be verbose
   * 
   * @return pointer or null
   */
  static TH1* GetH1(const TObject*     parent, 
		    const TString&     name,
		    Bool_t             verb=true)
  {
    return DoGetObject<TH1>(GetObject(parent, name, verb), parent);
  }
  //____________________________________________________________________
  /** 
   * Get a 2D histogram from a collection
   * 
   * @param parent Parent collection 
   * @param name   Name of histogram 
   * @param verb   If true and not found, complain
   * 
   * @return pointer or null
   */
  static TH2* GetH2(const TObject*     parent, 
		    const TString&     name, 
		    Bool_t             verb=true)
  {
    return DoGetObject<TH2>(GetObject(parent, name, verb), parent);
  }
  //____________________________________________________________________
  /** 
   * Get a 2D histogram from a collection
   * 
   * @param parent Parent collection 
   * @param name   Name of histogram 
   * @param verb   If true and not found, complain
   * 
   * @return pointer or null
   */
  static TH3* GetH3(const TCollection* parent, 
		    const TString&     name, 
		    Bool_t             verb=true)
  {
    // Info("GetH2", "Getting 2D histogram of %s from %p", name.Data(), c);
    // --- Find the object -------------------------------------------
    return DoGetObject<TH3>(GetObject(parent, name, verb), parent);
  }
  //__________________________________________________________________
  /** 
   * Get a histogram stack from a collection
   * 
   * @param parent Parent collection 
   * @param name   Name of histogram 
   * @param sub    If set, fill from sub-component 
   * @param verb   If true and not found, complain
   * 
   * @return pointer or null
   */
  static THStack* GetStack(const TObject*  parent, 
			   const TString&  name,
			   const char*     sub=0,
			   Bool_t          verb=true)
  {
    THStack* stack = DoGetObject<THStack>(GetObject(parent,name,verb),parent);
    if (!stack) return 0;
    if (sub == 0) return stack;
  
    if (stack->GetHists()->GetEntries() <= 0 ||stack->GetMaximum() < 1) { 
      // Info("GetStack", "No entries in %s", name.Data());
      stack->GetHists()->Delete();
      const char** ptr   = GetRingNames(false);
      while (*ptr) { 
	TCollection* sc = GetCollection(parent, *ptr, true);
	if (!sc) { ptr++; continue; }

	TObject* obj = GetObject(sc, sub);
	if (!obj) {
	  continue;
	  ptr++;
	}

	if (obj->IsA()->InheritsFrom(TH2::Class())) {
	  TH2* h = static_cast<TH2*>(obj);
	  TH1* p = h->ProjectionX(*ptr, 1, h->GetNbinsY(), "e");
	  p->Scale(1., "width");
	  p->SetTitle(*ptr);
	  p->SetDirectory(0);
	  stack->Add(p);
	}
	else if (obj->IsA()->InheritsFrom(TH1::Class())) {
	  TH1* hh = static_cast<TH1*>(obj);
	  hh->SetTitle(*ptr);
	  stack->Add(hh);
	}
	ptr++;
      }
    }
    // --- Return the collection -------------------------------------
    return stack;
  }
  //____________________________________________________________________
  /** 
   * Clear canvas 
   * 
   */
  void ClearCanvas()
  {
    if (fTop) {
      fTop->Clear();
      fTop->SetNumber(1);
      fTop->SetFillColor(kBlue-5);
      fTop->SetBorderSize(0);
      fTop->SetBorderMode(0);
    }

    fBody->Clear();
    fBody->SetNumber(2);
    fBody->SetFillColor(0);
    fBody->SetFillStyle(0);
    fBody->SetBorderSize(0);
    fBody->SetBorderMode(0);
    fBody->SetTopMargin(0.01);
    fBody->SetLeftMargin(0.10);
    fBody->SetRightMargin(0.01);
    fBody->SetBottomMargin(0.10);

    fRingMap[0] = 0;
    fRingMap[1] = 0;
    fRingMap[2] = 0;
    fRingMap[3] = 0;
    fRingMap[4] = 0;
    fRingMap[5] = 0;
    
    fCanvas->cd();    
  }
  //____________________________________________________________________
  /** 
   * Create a canvas 
   * 
   * @param pname     Name of PDF file to make 
   * @param landscape If true, print in landscape 
   * @param pdf       Make PDF
   * @param useTop    Make top area
   *
   * @return Created canvas 
   */
  void CreateCanvas(const TString& pname, 
		    Bool_t landscape=false, 
		    Bool_t pdf=true,
		    Bool_t useTop=true)
  {
    // Info("CreateCanvas", "Creating canvas");
    fLandscape = landscape;
    fPDF       = pdf;
    Int_t height = 1000;
    Int_t width  = height / TMath::Sqrt(2);
    if (fLandscape) {
      Int_t tmp = height; 
      height    = width;
      width     = tmp;
    }
    fCanvas = new TCanvas("c", pname.Data(), width, height);
    fCanvas->SetFillColor(0);
    fCanvas->SetBorderSize(0);
    fCanvas->SetBorderMode(0);
    if (fPDF) 
      fCanvas->Print(Form("%s[", pname.Data()), 
		     Form("pdf %s", fLandscape ? "Landscape" : ""));
    fCanvas->SetLeftMargin(.1);
    fCanvas->SetRightMargin(.05);
    fCanvas->SetBottomMargin(.1);
    fCanvas->SetTopMargin(.05);
  
    gStyle->SetOptStat(0);
    gStyle->SetTitleColor(0);
    gStyle->SetTitleStyle(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleX(.5);
    gStyle->SetTitleY(1);
    gStyle->SetTitleW(.8);
    gStyle->SetTitleH(.09);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetFrameBorderSize(1);
    gStyle->SetFrameBorderMode(1);
    gStyle->SetPalette(1);

    Float_t dy = useTop ? .05 : 0;
    if (useTop) {
      fTop = new TPad("top", "Top", 0, 1-dy, 1, 1, 0, 0);
      fTop->SetNumber(1);
      fTop->SetFillColor(kBlue-5);
      fTop->SetBorderSize(0);
      fTop->SetBorderMode(0);
      fCanvas->cd();
      fTop->Draw();
    }

    fBody = new TPad("body", "Body", 0, 0, 1, 1-dy, 0, 0);
    fBody->SetNumber(2);
    fBody->SetFillColor(0);
    fBody->SetFillStyle(0);
    fBody->SetBorderSize(0);
    fBody->SetBorderMode(0);
    fCanvas->cd();
    fBody->Draw();

    fHeader = new TLatex(.5, .5, "Title");
    fHeader->SetNDC();
    fHeader->SetTextAlign(22);
    fHeader->SetTextSize(.7);
    fHeader->SetTextColor(kWhite);
    fHeader->SetTextFont(62);

    Double_t x1 = .1;
    Double_t x2 = .6;
    Double_t y  = .8;
    Double_t s  = fLandscape ? 0.08 : 0.05;
    fParName = new TLatex(x1, y, "");
    fParName->SetTextAlign(13);
    fParName->SetNDC();
    fParName->SetTextSize(s);
    fParName->SetTextFont(62);
    
    fParVal = new TLatex(x2, y, "");
    fParVal->SetTextAlign(13);
    fParVal->SetNDC();
    fParVal->SetTextSize(s);
    fParVal->SetTextFont(42);

    fCanvas->cd();
  }

  //____________________________________________________________________
  /** 
   * Close the PDF
   * 
   */
  void CloseCanvas()
  {
    // Info("CloseCanvas", "Closing canvas");
    // ClearCanvas();
    if (fPDF && fCanvas) {
      // Printf("Closing canvas with last title %s", fLastTitle.Data());
      fCanvas->Print(Form("%s]", fCanvas->GetTitle()),
		     Form("pdf %s Title:%s", 
			  fLandscape ? "Landscape" : "",
			  fLastTitle.Data()));
    }
    if (fCanvas)
      fCanvas->Close();
    fCanvas = 0;
  }

  //__________________________________________________________________
  /** 
   * Print the canvas 
   * 
   * @param title  Title 
   * @param size   Size of text 
   */
  void PrintCanvas(const TString& title, Float_t size=.7)
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
      tit.Form("pdf %s Title:%s", fLandscape ? "Landscape" : "",
	       title.Data());

#ifdef DEBUG
      Info("PrintCanvas", "Printing to %s (%s)", 
	   fCanvas->GetTitle(), tit.Data());
#else
      gSystem->RedirectOutput("/dev/null");
#endif
      fCanvas->Print(fCanvas->GetTitle(), tit);
#ifndef DEBUG
      gSystem->RedirectOutput(0);
#endif
      fLastTitle = title;
      Pause();
      
      ClearCanvas();
    }
  }
  //__________________________________________________________________
  /** 
   * Make a chapter page 
   * 
   * @param title Title 
   */
  void MakeChapter(const TString& title)
  {
    fBody->cd();

    TLatex* ltx = new TLatex(.5, .5, title);
    ltx->SetNDC();
    ltx->SetTextAlign(22);
    ltx->Draw();
    
    PrintCanvas(title);
  }
  //__________________________________________________________________
  /** 
   * Draw an object in pad 
   * 
   * @param c       PArent pad 
   * @param padNo   Sub-pad number (0 is self)
   * @param h       Object to draw 
   * @param opts    Options
   * @param flags   Flags
   * @param title   Title on plot
   *
   * @return Drawn object - if any
   */
  TObject* DrawInPad(TVirtualPad* c, 
		     Int_t        padNo, 
		     TObject*     h, 
		     Option_t*    opts="",
		     UInt_t       flags=0x0,
		     const char*  title="")
  {
    TVirtualPad* p = c->GetPad(padNo);
    if (!p) { 
      Warning("DrawInPad", "Pad # %d not found in %s", padNo, c->GetName());
      return 0;
    }
    return DrawInPad(p, h, opts, flags, title);
  }
  /** 
   * Draw a clone of an object
   * 
   * @param o       Object
   * @param options Draw options
   * @param title   Title of object
   *
   * @return Drawn object - if any
   */
  virtual TObject* DrawObjClone(TObject* o, Option_t* options, 
				const char* title)
  {
    if (o->IsA()->InheritsFrom(TH1::Class())) 
      return DrawObjClone(static_cast<TH1*>(o), options, title);
    else if (o->IsA()->InheritsFrom(THStack::Class())) 
      return DrawObjClone(static_cast<THStack*>(o), options, title);
    else if (o->IsA()->InheritsFrom(TGraph::Class()))
      return o->DrawClone(options);
    else 
      o->Draw(options);
    return o;
  }
  /** 
   * Draw an object clone 
   * 
   * @param o        Stack object
   * @param options  Draw options 
   * @param title    Title on plot
   *
   * @return Drawn object - if any
   */
  virtual TObject* DrawObjClone(THStack* o, Option_t* options, 
				const char* title)
  {
    // THStack* tmp = static_cast<THStack*>(o->Clone());
    o->Draw(options);
    if (title && title[0] != '\0') o->GetHistogram()->SetTitle(title);
    TAxis*   xAxis = o->GetXaxis();
    if (!xAxis) {
      Warning("DrawObjClone", "No X-axis for drawn stack %s", o->GetName());
      return o;
    }
    TH1*     h     = 0;
    Int_t    nBins = xAxis->GetNbins();
    Double_t xMin  = xAxis->GetXmin();
    Double_t xMax  = xAxis->GetXmax();
    TIter    next(o->GetHists());
    while ((h = static_cast<TH1*>(next()))) {
      TAxis* a = h->GetXaxis();
      nBins    = TMath::Max(nBins, a->GetNbins()); 
      xMin     = TMath::Min(xMin, a->GetXmin());
      xMax     = TMath::Max(xMax, a->GetXmax());
    }
    if (nBins != xAxis->GetNbins() || 
	xMin  != xAxis->GetXmin() || 
	xMax  != xAxis->GetXmax()) {
      xAxis->Set(nBins, xMin, xMax);
      o->GetHistogram()->Rebuild();
    }
    o->GetHistogram()->SetYTitle(title);
    return o;
  }
  /** 
   * Draw an object clone 
   * 
   * @param o        Histogram
   * @param options  Draw options 
   * @param title    Title on plot 
   *
   * @return Drawn object - if any
   */
  virtual TObject* DrawObjClone(TH1* o, Option_t* options, const char* title)
  {
    TH1* tmp = o->DrawCopy(options);
    if (title && title[0] != '\0') tmp->SetTitle(title);
    return tmp;
  }    
  //__________________________________________________________________
  static void GetLegendPosition(UInt_t    flags, TVirtualPad* p, 
				Double_t& x1,    Double_t&    y1, 
				Double_t& x2,    Double_t&    y2)
  {
    UInt_t   horiz = (flags & 0xF0000);
    UInt_t   verti = (flags & 0xF000);
    Double_t eps   = .01;
    Double_t dY    = .4;
    Double_t dX    = .4;
    Double_t yB    = p->GetBottomMargin()+eps;
    Double_t yT    = 1-p->GetTopMargin()-eps;
    Double_t xL    = p->GetLeftMargin()+eps;
    Double_t xR    = 1-p->GetRightMargin()-eps;
    switch (verti) { 
    case kNorth:  y1 = yT-dY;           break;
    case kSouth:  y1 = yB;              break;
    case kMiddle: y1 = (yB+yT-dY)/2;    break;
    }
    y2 = TMath::Min(y1 + dY, yT);
    
    switch (horiz) { 
    case kEast:   x1 = xL;              break;
    case kWest:   x1 = xR-dX;           break;
    case kCenter: x1 = (xL+xR-dX)/2;    break;
    }
    x2 = TMath::Min(x1 + dX, xR);
  } 
  //__________________________________________________________________
  /** 
   * Make a legend 
   * 
   * @param p 
   * @param flags 
   * @param autoFill 
   * 
   * @return 
   */
  TLegend* MakeLegend(TVirtualPad* p, UInt_t flags, Bool_t autoFill)
  {
    Double_t x1 = fParVal->GetX();
    Double_t y1 = fParVal->GetY();
    Double_t x2 = 0;
    Double_t y2 = 0;
    GetLegendPosition(flags, p, x1, y1, x2, y2);

    //Printf("Legend at (%f,%f)x(%f,%f)", x1, y1, x2, y2);
    TLegend* l = 0;
    p->cd();
    if (autoFill) l = p->BuildLegend(x1, y1, x2, y2);
    else          l = new TLegend(x1, y1, x2, y2);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);

    return l;
  }
  //__________________________________________________________________
  /** 
   * Draw an object in pad 
   * 
   * @param p       Pad
   * @param h       Object to draw 
   * @param opts    Options
   * @param flags   Flags
   * @param title   Title on plot
   *
   * @return Drawn object - if any
   */
  TObject* DrawInPad(TVirtualPad* p, 
		     TObject*     h, 
		     Option_t*    opts="",
		     UInt_t       flags=0x0,
		     const char*  title="")
  {
    if (!p) { 
      Warning("DrawInPad", "No pad specified");
      return 0;
    }
    p->cd();
    // Info("DrawInPad", "Drawing in pad %p", p);
    // fBody->ls();
    if (flags & kLogx) p->SetLogx();
    if (flags & kLogy) p->SetLogy();
    if (flags & kLogz) p->SetLogz();
    if (flags & kGridx) p->SetGridx();
    if (flags & kGridy) p->SetGridy();
    // if (flags & kGridz) p->SetGridz();
    p->SetFillColor(0);
    TString o(opts);
    if (o.Contains("colz", TString::kIgnoreCase)) 
      p->SetRightMargin(0.15);
    if (!h) {
      if (!(flags & kSilent))
	Warning("DrawInPad", "Nothing to draw in pad # %s", p->GetName());
      return 0;
    }
    if (o.Contains("text", TString::kIgnoreCase)) {
      TH1* hh = static_cast<TH1*>(h);
      hh->SetMaximum(1.1*hh->GetMaximum());
      hh->SetMarkerSize(2);
      o.Append("30");
    }
    TObject* ret = DrawObjClone(h, o, title);
    
    if (flags & kLegend) {
      MakeLegend(p, flags, true);
    }
    p->Modified();
    p->Update();
    p->cd();

    return ret;
  }
  //__________________________________________________________________
  /** 
   * Draw two graphs in the same frame, but with separate y-axis 
   * 
   * @param c      Mother pad 
   * @param padNo  Sub-pad number (0 is self)
   * @param h1     First histogram
   * @param h2     Second histogram
   * @param opts   Options
   * @param flags  Flags
   */
  void DrawTwoInPad(TVirtualPad* c, Int_t padNo, TH1* h1, TH1* h2,
		    Option_t* opts="", UShort_t flags=0x0)
  {
    TVirtualPad* p = c->cd(padNo);
    if (!p) { 
      Warning("DrawInPad", "Pad # %d not found in %s", padNo, c->GetName());
      return;
    }
    if (flags & kLogx) p->SetLogx();
    if (flags & kLogy) p->SetLogy();
    if (flags & kLogz) p->SetLogz();
    if (flags & kGridx) p->SetGridx();
    if (flags & kGridy) p->SetGridy();
    // if (flags & kGridz) p->SetGridz();
    p->SetFillColor(0);

    TString o(opts);
    o.ToLower();
    TString fopts(o);
    Bool_t e3 = o.Contains("e3");
    if (e3) {
      fopts.ReplaceAll("e3", " same");
    }
    
    h1->GetYaxis()->SetLabelSize(0);
    h1->GetYaxis()->SetTicks("");
    h1->GetYaxis()->SetNdivisions(0);
    h1->DrawCopy(o); // First draw with opts 
    if (e3) h1->DrawCopy(fopts);
    p->Update();

    Double_t m1 = 1.05 * h1->GetMaximum();
    if (m1 > 0) {
      TGaxis*  a1 = new TGaxis(p->GetUxmin(), p->GetUymin(), 
			       p->GetUxmin(), p->GetUymax(), 
			       0, m1, 510);
      a1->SetLineColor(h1->GetLineColor());
      a1->Draw();
    }

    o.Append(" same");
    Double_t m2    = 1.1 * h2->GetMaximum();
    Double_t scale = m2 > 0 ? m1 / m2 : 1;
    h2->Scale(scale);
    h2->DrawCopy(o);
    if (e3) h2->DrawCopy(fopts);

    if (m2 > 0) {
      TGaxis*  a2 = new TGaxis(p->GetUxmax(), p->GetUymin(), 
			       p->GetUxmax(), p->GetUymax(), 
			       0, m2, 510, "+L");
      a2->SetLineColor(h2->GetLineColor());
      a2->Draw();
    }
    if (flags& kLegend) { 
      MakeLegend(p, flags, true);
    }
    p->Modified();
    p->Update();
    p->cd();
  }  
  
  //____________________________________________________________________
  /** 
   * Draw a parameter.  
   * 
   * @param y       Current y position. On return new y position
   * @param name    Parameter name
   * @param value   Parameter value 
   * @param size    Optional text size
   */
  void DrawParameter(Double_t&      y, 
		     const TString& name, 
		     const TString& value,
		     Double_t       size=0)
  {
    Double_t s = fParName->GetTextSize();
    Double_t t = fParVal->GetTextSize();
    if (name.IsNull() && value.IsNull()) return;
    if (size > 0) { 
      fParName->SetTextSize(size);
      fParVal->SetTextSize(size);
    }
    if (!name.IsNull())
      fParName->DrawLatex(fParName->GetX(), y, Form("%s:", name.Data()));
    if (!value.IsNull())
      fParVal->DrawLatex(fParVal->GetX(), y, value.Data());
    if (!name.IsNull())
      y -= 1.2 * fParName->GetTextSize();
    else if (!value.IsNull())
      y -= 1.2 * fParVal->GetTextSize();

    fParName->SetTextSize(s);
    fParVal->SetTextSize(t);
  }  
  template <typename T>
  void DrawTParameter(Double_t&      y,
		      TList*         list, 
		      const TString& name) {
    T value;
    if (!GetParameter(list, name, value)) 
      return;
    std::stringstream s;
    s << std::boolalpha << value;
    DrawParameter(y, name, s.str().c_str(), 0);
  }

#ifndef __CINT__
  //__________________________________________________________________
  /**
   * Structure to hold a dived pad 
   */
  struct DividedPad { 
    TVirtualPad*  fParent;
    TVirtualPad** fSubs;
    Bool_t        fLandscape;
    Int_t         fNCol;
    Int_t         fNRow;

    DividedPad(TVirtualPad* p, Bool_t landscape, Int_t nCol, Int_t nRow) 
      : fParent(p), 
	fSubs(0),
	fLandscape(landscape),
	fNCol(landscape ? nRow : nCol),
	fNRow(landscape ? nCol : nRow)
    {
      Int_t nPad = fNCol * fNRow;
      fSubs      = new TVirtualPad*[nPad];
    }
    void Divide(Bool_t commonX, Bool_t commonY) {
      if ((!commonX && !commonY) || (commonX && commonY)) {
	// In case we have no common axis or do have both to be common,
	// we directly use the TVirtualPad::Divide member function 
	fParent->Divide(fNCol, fNRow, commonX ? 0 : 0.01, commonY ? 0 : 0.01);
	for (Int_t iPad = 1; iPad <= fNRow*fNCol; iPad++) 
	  fSubs[iPad-1] = fParent->GetPad(iPad);
      }
      else if (commonX && !commonY) {
	// We need to have common X axis, but not common Y axis. We first
	// divide the pad in fNCol columns, and then each in to fNRow rows
	fParent->Divide(fNCol, 1);
	for (Int_t iCol = 1; iCol <= fNCol; iCol++) { 
	  TVirtualPad* q = fParent->GetPad(iCol);

	  if (fNRow == 1) {
	    fSubs[GetIdx(iCol,0)] = q;
	    continue;
	  }

	  q->Divide(1,fNRow,0,0);
	  for (Int_t iRow = 1; iRow <= fNRow; iRow++) 
	    fSubs[GetIdx(iCol, iRow)] = q->GetPad(iRow);
	}
      }
      else if (!commonX && commonY) { 
	// We need to have common Y axis, but not common X axis. We first
	// divide the pad in fNRow rows, and then each in to fNCol columns
	fParent->Divide(1, fNRow);
	for (Int_t iRow = 1; iRow <= fNRow; iRow++) { 
	  TVirtualPad* q = fParent->GetPad(iRow);

	  if (fNCol == 1) {
	    fSubs[GetIdx(0,iRow)] = q;
	    continue;
	  }
	  
	  q->Divide(fNCol,1,0,0);
	  for (Int_t iCol = 1; iCol <= fNCol; iCol++) 
	    fSubs[GetIdx(iCol, iRow)] = q->GetPad(iCol);
	}
      }
    }
    virtual ~DividedPad() { if (fSubs) delete [] fSubs; }
    /** 
     * Get a sub-pad 
     * 
     * @param idx Index (0 based)
     * 
     * @return Pad or null
     */
    TVirtualPad* GetPad(Int_t idx) {
      if (!fSubs) {
	::Warning("GetPad","No sub-pads");
	return 0;
      }
      if (idx < 0 || idx >= (fNRow*fNCol)) {
	::Warning("GetPad", "Inded %d out of bounds [%d,%d]", 
		  idx, 0, fNRow*fNCol);
	return 0;
      }
      return fSubs[idx];
    }
    Int_t GetIdx(Int_t iCol, Int_t iRow) const 
    {
      return (iRow-1) * fNCol + iCol;
    }
    /** 
     * Get a sub-pad 
     * 
     * @param iRow  Row number (1-based)
     * @param iCol  Column number (1-based)
     * 
     * @return Pad or null
     */
    TVirtualPad* GetPad(Int_t iCol, Int_t iRow) { 
      if (iRow < 0 || iRow > fNRow) return 0;
      if (iCol < 0 || iRow > fNCol) return 0;
      return GetPad(GetIdx(iCol, iRow));
    }
  };
#endif
  
  //__________________________________________________________________
  void DivideForRings(Bool_t commonX, Bool_t commonY)
  {
    // 
    // Divide canvas for rings 
    // 
    if ((!commonX && !commonY) || 
	(commonX  && commonY)) {
    // Portrait:
    //    +----------+----------+
    //    | 1: FMD1i | 2: Free  |
    //    +----------+----------+
    //    | 3: FMD2i | 4: FMD2o |
    //    +----------+----------+
    //    | 5: FMD3i | 6: FMD3o |
    //    +----------+----------+
    // 
    // Landscape:
    //    +----------+----------+----------+
    //    | 1: FMD1i | 2: FMD2i | 3: FMD3i |
    //    +----------+----------+----------+
    //    | 4: Free  | 5: FMD2o | 6: FMD3o |
    //    +----------+----------+----------+
    // 
      fBody->Divide(fLandscape ? 3 : 2, fLandscape ? 2 : 3,
		    commonX ? 0 : 0.01, commonY ? 0 : 0.01);
      fRingMap[0] = fBody->GetPad(1);                  // FMD1i;
      fRingMap[1] = fBody->GetPad(fLandscape ? 2 : 3); // FMD2i;
      fRingMap[2] = fBody->GetPad(fLandscape ? 5 : 4); // FMD2o;
      fRingMap[3] = fBody->GetPad(fLandscape ? 3 : 5); // FMD3i;
      fRingMap[4] = fBody->GetPad(6);                  // FMD3o;
      fRingMap[5] = fBody->GetPad(fLandscape ? 4 : 2); // Free
    }
    else if (commonX && !commonY) {
      // Divide into two  - left/right
      // Portrait:
      //    +----------++----------+
      //    | 1: FMD1i || 1: Free  |
      //    +----------++----------+
      //    | 2: FMD2i || 2: FMD2o |
      //    +----------++----------+
      //    | 3: FMD3i || 3: FMD3o |
      //    +----------++----------+
      // 
      // Landscape:
      //    +----------++----------++----------+
      //    | 1: FMD1i || 1: FMD2i || 1: FMD3i |
      //    +----------++----------++----------+
      //    | 2: Free  || 2: FMD2o || 2: FMD3o |
      //    +----------++----------++----------+
      // 
      fBody->Divide(fLandscape ? 3 : 2, 1);
      TVirtualPad* left = fBody->cd(1);
      left->Divide(fLandscape ? 2 : 3);
      TVirtualPad* middle = fBody->cd(2);
      middle->Divide(fLandscape ? 2 : 3);

      Info("left","%p",left); left->ls();
      Info("middle","%p",middle); middle->ls();

      fRingMap[0] = left->GetPad(1); // FMD1i;
      if (!fLandscape) {
	fRingMap[1] = left->GetPad(2);   // FMD2i
	fRingMap[2] = middle->GetPad(2); // FMD2o
	fRingMap[3] = left->GetPad(3);   // FMD3i
	fRingMap[4] = middle->GetPad(3); // FMD3o
	fRingMap[5] = middle->GetPad(1); // Free
      }
      else {
	TVirtualPad* right = fBody->cd(3);
	right->Divide(fLandscape ? 2 : 3);
	fRingMap[1] = middle->GetPad(1); // FMD2i
	fRingMap[2] = middle->GetPad(2); // FMD2o
	fRingMap[3] = right->GetPad(1);  // FMD3i
	fRingMap[4] = right->GetPad(2);  // FMD3o
	fRingMap[5] = left->GetPad(2);   // Free
      }
    }
    else { 
      // Divide into two  - left/right
      // Portrait:
      //    +----------+----------+
      //    | 1: FMD1i | 2: Free  |
      //    +----------+----------+
      //    +----------+----------+
      //    | 1: FMD2i | 2: FMD2o |
      //    +----------+----------+
      //    +----------+----------+
      //    | 1: FMD3i | 2: FMD3o |
      //    +----------+----------+
      // 
      // Landscape:
      //    +----------+----------+----------+
      //    | 1: FMD1i | 2: FMD2i | 3: FMD3i |
      //    +----------+----------+----------+
      //    +----------+----------+----------+
      //    | 1: Free  | 2: FMD2o | 3: FMD3o |
      //    +----------+----------+----------+
      // 
      fBody->Divide(1, fLandscape ? 2 : 3);
      TVirtualPad* top = fBody->cd(1);
      top->Divide(fLandscape ? 3 : 2);
      TVirtualPad* middle = fBody->cd(2);
      middle->Divide(fLandscape ? 3 : 2);

      fRingMap[0] = top->GetPad(1); // FMD1i;
      if (!fLandscape) {
	TVirtualPad* bottom = fBody->cd(2);
	bottom->Divide(2);

	fRingMap[1] = middle->GetPad(1); // FMD2i
	fRingMap[2] = middle->GetPad(2); // FMD2o
	fRingMap[3] = bottom->GetPad(1); // FMD3i
	fRingMap[4] = bottom->GetPad(2); // FMD3o
	fRingMap[5] = top->GetPad(2);    // Free
      }
      else {
	fRingMap[1] = top->GetPad(2);    // FMD2i
	fRingMap[2] = middle->GetPad(2); // FMD2o
	fRingMap[3] = top->GetPad(3);    // FMD3i
	fRingMap[4] = middle->GetPad(3); // FMD3o
	fRingMap[5] = middle->GetPad(1); // Free
      }
    }
    if (fRingMap[0]) fRingMap[0]->SetTitle("FMD1i");
    if (fRingMap[1]) fRingMap[1]->SetTitle("FMD2i");
    if (fRingMap[2]) fRingMap[2]->SetTitle("FMD2o");
    if (fRingMap[3]) fRingMap[3]->SetTitle("FMD3i");
    if (fRingMap[4]) fRingMap[4]->SetTitle("FMD3o");
    if (fRingMap[5]) fRingMap[5]->SetTitle("Other");
  }
  //__________________________________________________________________
  TVirtualPad* RingPad(UShort_t d, Char_t r) const
  {
    Int_t idx = 0;
    switch (d) { 
    case 0:     idx = 5; break;
    case 1:     idx = 0; break;
    case 2:     idx = 1 + ((r == 'I' || r == 'i') ? 0 : 1); break;
    case 3:     idx = 3 + ((r == 'I' || r == 'i') ? 0 : 1); break;
    default: return 0;
    }
    return fRingMap[idx];
    // return fBody->GetPad(no);
  }
  //__________________________________________________________________
  TVirtualPad* RingPad(const char* name) const
  {
    TString n(name);
    Int_t idx = n.Index("FMD");
    if (n == kNPOS) return 0;
    n.Remove(0, idx+3);
    Int_t det = n.Atoi();
    n.Remove(0,1);
    Char_t rng = n[0];
    return RingPad(det, rng);
  }
  //__________________________________________________________________
  /** 
   * Draw an object in pad 
   * 
   * @param d       Detector 
   * @param r       Ring 
   * @param h       Object to draw 
   * @param opts    Options
   * @param flags   Flags
   * @param title   Title on plot
   */
  void DrawInRingPad(UShort_t    d, 
		     Char_t      r, 
		     TObject*    h, 
		     Option_t*   opts="",
		     UShort_t    flags=0x0,
		     const char* title="")
  {
    TVirtualPad* p = RingPad(d, r);
    if (!p) {
      Warning("DrawInRingPad", "No pad found for FMD%d%c", d, r);
      return;
    }
    DrawInPad(p, h, opts, flags, title);
  }
  /** 
   * Draw object in a ring pad
   * 
   * @param name    Name of ring
   * @param h       Object to draw
   * @param opts    Options
   * @param flags   Flags
   * @param title   Possible new title
   */
  void DrawInRingPad(const char* name, 
		     TObject*    h, 
		     Option_t*   opts="", 
		     UShort_t    flags=0x0, 
		     const char* title="")
  {
    TVirtualPad* p = RingPad(name);
    if (!p) {
      Warning("DrawInRingPad", "No pad found for \"%s\"", name);
      return;
    }
    DrawInPad(p, h, opts, flags, title);
  }
  /** 
   * Draw object in a ring pad. Which pad to draw in depends on the
   * name or title of the drawn object (must contain the ring name as
   * a sub-string).
   * 
   * @param h       Object to draw
   * @param opts    Options
   * @param flags   Flags
   * @param title   Possible new title
   */
  void DrawInRingPad(TObject*    h, 
		     Option_t*   opts="", 
		     UShort_t    flags=0x0, 
		     const char* title="")
  {
    if (!h) return;
    TVirtualPad* p = RingPad(h->GetName());
    if (!p) {
      p = RingPad(h->GetTitle());
      if (!p) {
	Warning("DrawInRingPad", "No pad found for %s/%s", 
		h->GetName(), h->GetTitle());
	return;
      }
    }
    DrawInPad(p, h, opts, flags, title);
  }
    

  //__________________________________________________________________
  /** 
   * Pause after each plot
   * 
   */
  void Pause()
  {
    if (!fPause) return;    
    // printf("Press space to continue");
    fCanvas->WaitPrimitive();
    // std::cin.get();
  }
  static void CompileScript(const TString& name, 
			    const TString& sub, 
			    const TString& check,
			    Bool_t         force)
  {
    if (!check.IsNull() && gROOT->GetClass(check)) return;

    TString fwd =gSystem->ExpandPathName("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2");
    TString macPath(gROOT->GetMacroPath());
    TString incPath(gSystem->GetIncludePath());
    if (!macPath.Contains(fwd)) macPath.Append(Form(":%s", fwd.Data()));
    if (!incPath.Contains(fwd)) gSystem->AddIncludePath(Form("-I%s",
							     fwd.Data()));
    if (!sub.IsNull()) { 
      TObjArray* subs = sub.Tokenize(": ");
      TObject*   pSub = 0;
      TIter      iSub(subs);
      while ((pSub = iSub())) {
	TString subDir = gSystem->ConcatFileName(fwd, pSub->GetName());
	if (!macPath.Contains(subDir))
	  macPath.Append(Form(":%s", subDir.Data()));
	if (!incPath.Contains(subDir)) 
	  gSystem->AddIncludePath(Form("-I%s", subDir.Data()));
      }
    }
    gROOT->SetMacroPath(macPath);
    gROOT->LoadMacro(Form("%s%s", name.Data(), (force ? "++g" : "+")));
  }
  //____________________________________________________________________
  virtual void DrawEventInspector(TCollection* parent)
  {
    Info("DrawEventInspector", "Drawing event inspector");
    TCollection* c = GetCollection(parent, "fmdEventInspector");
    if (!c) return;

    UShort_t sys=0, sNN=0;
    Int_t field=0;
    ULong_t runNo=0;
    Int_t lowFlux=0, nPileUp=0, ipMethod=0;
    ULong_t aliRev=0, aliBra=0;    
    Bool_t v0and=false;
    Double_t dPileUp=0.;
    Double_t y = .8;

    fBody->cd();

    Double_t save = fParName->GetTextSize();
    fParName->SetTextSize(0.03);
    fParVal->SetTextSize(0.03);
    
    GetParameter(c, "sys", sys);
    GetParameter(c, "sNN", sNN);
    GetParameter(c, "field", field);
    GetParameter(c, "runNo", runNo);
    GetParameter(c, "lowFlux", lowFlux);
    GetParameter(c, "ipMethod", ipMethod, false);
    GetParameter(c, "v0and", v0and);
    GetParameter(c, "nPileUp", nPileUp);
    GetParameter(c, "dPileup", dPileUp);
    GetParameter(c, "alirootRev", aliRev);
    GetParameter(c, "alirootBranch", aliBra);
    
    TString tS; SysString(sys, tS); DrawParameter(y, "System", tS);
    TString tE; SNNString(sNN, tE); DrawParameter(y, "#sqrt{s_{NN}}", tE);
    DrawParameter(y, "L3 B field", Form("%+2dkG", field));
    DrawParameter(y, "Run #", Form("%lu", runNo));
    DrawParameter(y, "Low flux cut", Form("%d", lowFlux));
    TString sIpMeth("unknown");
    switch(ipMethod) { 
    case 0: sIpMeth = "Normal"; break;
    case 1: sIpMeth = "pA in 2012"; break;
    case 2: sIpMeth = "pA in 2013"; break;
    case 3: sIpMeth = "PWG-UD"; break;
    case 4: sIpMeth = "Satellite"; break;
    }
    DrawParameter(y, "Use PWG-UD vertex", sIpMeth);
    DrawParameter(y, "Use V0AND for NSD", (v0and ? "yes" : "no"));
    DrawParameter(y, "Least # of pile-up vertex", Form("%d", nPileUp));
    DrawParameter(y, "Least distance of pile-up vertex",
		  Form("%fcm", dPileUp));
    DrawParameter(y, "AliROOT", Form("%lu/0x%08lx", ULong_t(aliRev), 
				     ULong_t(aliBra)));

    TH1*    triggers     = GetH1(c, "triggers");
    TH1*    vertex       = GetH1(c, "vertex", false);
    Bool_t  mc           = (vertex != 0);
    if (mc) { 
      Int_t nInelMC = vertex->GetEntries();
      Int_t nInel   = triggers->GetBinContent(1);
      Int_t nNSDMC  = triggers->GetBinContent(11);
      Int_t nNSD    = triggers->GetBinContent(4);
      DrawParameter(y, 
		    Form("#varepsilon_{INEL} = #bf{%d/%d}", nInel, nInelMC),
		    Form("%5.3f", float(nInel)/nInelMC));
      DrawParameter(y, 
		    Form("#varepsilon_{NSD} = #bf{%d/%d}", nNSD, nNSDMC),
		    Form("%5.3f", float(nNSD)/nNSDMC));
    }

    PrintCanvas("Event Inspector");
    fParName->SetTextSize(save);
    fParVal->SetTextSize(save);

    if (fLandscape) fBody->Divide(4,2);
    else            fBody->Divide(2,4);
    
    TH1*    nEventsTr    = GetH1(c, "nEventsTr");
    TH1*    nEventsTrVtx = GetH1(c, "nEventsTrVtx");
    TH1*    nEventsAcc   = GetH1(c, "nEventsAccepted");
    if (nEventsTr)    nEventsTr->Rebin(2);
    if (nEventsTrVtx) nEventsTrVtx->Rebin(2);
    if (vertex) {
      // vertex->Rebin(2);
      vertex->SetFillColor(kMagenta+2);
    }
    DrawInPad(fBody, 1, nEventsTr, "", kLogy, 
	      "Events w/trigger, trigger+vertex, accepted");
    if (vertex) DrawInPad(fBody, 1, vertex, "same");
    DrawInPad(fBody, 1, nEventsTrVtx, "same"); 
    DrawInPad(fBody, 1, nEventsAcc, "same", kLegend);


    DrawInPad(fBody, 2, GetH2(c, "nEventsAcceptedXY"), "colz", kLogz);
    DrawInPad(fBody, 3, triggers,          "hist text");
    if (GetH1(c, "trgStatus"))
      DrawInPad(fBody, 4, GetH1(c, "trgStatus"),       "hist text");
    else  // Old one 
      DrawInPad(fBody, 4, GetH2(c, "triggerCorr"),     "colz", kLogz);
    DrawInPad(fBody, 5, GetH1(c, "status"),            "hist text");
    if (GetH1(c, "vtxStatus"))
      DrawInPad(fBody, 6, GetH1(c, "vtxStatus"),       "hist text");
    else // old 
      DrawInPad(fBody, 6, GetH1(c, "type"),            "hist text");

    TH1* cent     = GetH1(c, "cent");
    if (cent) { 
      cent->Scale(1, "width");
      DrawInPad(fBody, 7, cent, "", kLogy);
    }

    TH1* pileupStatus = GetH1(c, "pileupStatus", false);
    if (pileupStatus) DrawInPad(fBody, 8, pileupStatus, "hist text30");
    else {
      TH2* centQual = GetH2(c, "centVsQuality");
      if (centQual) { 
	centQual->Scale(1, "width");
	DrawInPad(fBody, 8, centQual, "colz", kLogz);
      }
    }
    
    PrintCanvas("EventInspector - Histograms");  

    if (!mc) return; // not MC 
  
    TH1* phiR         = GetH1(c, "phiR");
    TH1* b            = GetH1(c, "b");
    TH2* bVsNpart     = GetH2(c, "bVsParticipants");
    TH2* bVsNbin      = GetH2(c, "bVsBinary");
    TH2* bVsCent      = GetH2(c, "bVsCentrality");
    TH2* vzComparison = GetH2(c, "vzComparison");
    TH2* centVsNpart  = GetH2(c, "centralityVsParticipans");// Spelling!
    TH2* centVsNbin   = GetH2(c, "centralityVsBinary");
  
    fBody->Divide(2,3);

    DrawInPad(fBody, 1, phiR);
    DrawInPad(fBody, 2, vzComparison, "colz", kLogz);
    DrawInPad(fBody, 3, b);

    TProfile* nPartB = bVsNpart->ProfileX("nPartB",1,-1,"s");
    TProfile* nBinB  = bVsNbin->ProfileX("nBinB",1,-1,"s");
    nPartB->SetMarkerColor(kBlue+2);
    nPartB->SetMarkerStyle(20);
    nPartB->SetLineColor(kBlue+2);
    nPartB->SetFillColor(kBlue-10);
    nPartB->SetFillStyle(1001);
    nPartB->SetMarkerSize(0.7);
    nBinB->SetMarkerColor(kRed+2);
    nBinB->SetMarkerStyle(21);
    nBinB->SetLineColor(kRed+2);
    nBinB->SetFillColor(kRed-10);
    nBinB->SetMarkerSize(0.7);
    nBinB->SetFillStyle(1001);

    DrawTwoInPad(fBody, 4, nPartB, nBinB, "e3 p", kLegend);

    DrawInPad(fBody, 5, bVsCent, "colz", kLogz);

    TProfile* nPartC = centVsNpart->ProfileY("nPartC",1,-1,"s");
    TProfile* nBinC  = centVsNbin->ProfileY("nBinC",1,-1,"s");
    nPartC->SetMarkerColor(kBlue+2);
    nPartC->SetMarkerStyle(20);
    nPartC->SetLineColor(kBlue+2);
    nPartC->SetFillColor(kBlue-10);
    nPartC->SetFillStyle(1001);
    nPartC->SetMarkerSize(0.7);
    nBinC->SetMarkerColor(kRed+2);
    nBinC->SetMarkerStyle(21);
    nBinC->SetLineColor(kRed+2);
    nBinC->SetFillColor(kRed-10);
    nBinC->SetMarkerSize(0.7);
    nBinC->SetFillStyle(1001);

    DrawTwoInPad(fBody, 6, nPartC, nBinC, "e3 p", kLegend);

    PrintCanvas("EventInspector - Monte-Carlo");  
  }
  //____________________________________________________________________
  virtual void DrawESDFixer(TCollection* parent)
  {
    Info("DrawESDFixer", "Drawing ESD fixer");
    TCollection* c = GetCollection(parent, "fmdESDFixer");
    if (!c) return;

    Int_t  recoFactor = 0;
    Bool_t recalcEta = false;
    Bool_t invalidIsEmpty = false;

    fBody->cd();

    Double_t save = fParName->GetTextSize();
    fParName->SetTextSize(0.05);
    fParVal->SetTextSize(0.05);

    fBody->Divide(2,2);
    fBody->cd(1);
    
    Double_t y = .8;
    if (GetParameter(c, "recoFactor", recoFactor))
	DrawParameter(y, "Noise factor used in reco",
		      Form("%d (assumed)", recoFactor));
    if (GetParameter(c, "recalcEta", recalcEta))
	DrawParameter(y, "Recalculate #eta",
		      Form("%s", (recalcEta ? "yes" : "no")));
    if (GetParameter(c, "invalidIsEmpty", invalidIsEmpty))
	DrawParameter(y, "Assume invalid strips are empty",
		      Form("%s", (invalidIsEmpty ? "yes" : "no")));

    TCollection* xd = GetCollection(c, "extraDead");
    if (xd) 
      DrawParameter(y, "# extra dead strips", 
		    Form("%d", xd->GetEntries()));
    
    DrawInPad(fBody, 2, GetH1(c, "noiseChange"), "", kLogy);
    DrawInPad(fBody, 3, GetH1(c, "etaChange"), "", kLogy);
    DrawInPad(fBody, 4, GetH1(c, "deadChange"), "", kLogy);
	
    PrintCanvas("ESD Fixer");
    fParName->SetTextSize(save);
    fParVal->SetTextSize(save);
  }
  //____________________________________________________________________
  void DrawTrackDensity(TCollection* parent, 
			const char* folderName="mcTrackDensity")
  {
    Info("DrawTrackDensity", "Drawing track density");

    // --- MC --------------------------------------------------------
    TCollection* mc = GetCollection(parent, folderName, false);
    if (!mc) return; // Not MC 

    fBody->Divide(2,3);
    DrawInPad(fBody, 1, GetH2(mc, "binFlow"),    "colz", kLogz);
    DrawInPad(fBody, 2, GetH2(mc, "binFlowEta"), "colz", kLogz);
    DrawInPad(fBody, 3, GetH2(mc, "binFlowPhi"), "colz", kLogz);
    DrawInPad(fBody, 4, GetH1(mc, "nRefs"),       "",    kLogy,
	      "# of references");
    DrawInPad(fBody, 4, GetH1(mc, "clusterRefs",   false), "same");
    DrawInPad(fBody, 4, GetH1(mc, "clusterSize",   false), "same");
    DrawInPad(fBody, 4, GetH1(mc, "nClusters",     false), "same", kLegend);
    DrawInPad(fBody, 5, GetH2(mc, "clusterVsRefs", false),"colz", kLogz);

    PrintCanvas("Track density");  
  }
  //__________________________________________________________________
  TCanvas* fCanvas;  // Our canvas 
  TPad*    fTop;     // Top part 
  TPad*    fBody;    // Body part 
  TLatex*  fHeader;  // Header text 
  TLatex*  fParName; // Parameter name 
  TLatex*  fParVal;  // Parameter value 
  Bool_t   fPause;   // Whether to pause after drawing a canvas
  Bool_t   fLandscape; // Landscape or Portrait orientation
  TVirtualPad** fRingMap;
  Bool_t   fPDF;
  TString  fLastTitle;
};
#if 0
template <> 
inline Bool_t 
SummaryDrawer::DoGetParameter<Double_t>(TObject* o, const TObject* p, 
					Double_t& value)
{
  TParameter<Double_t>* r = DoGetObject<TParameter<Double_t>(o, p);
  UInt_t  i = o->GetUniqueID();
  Float_t v = *reinterpret_cast<Float_t*>(&i);
  value = v;
  return true;
}
#endif
inline Bool_t 
SummaryDrawer::GetParameter(const TObject*   c, 
			    const TString&   name, 
			    Double_t&        value,
			    Bool_t           verb)
  
{
  return DoGetParameter(GetObject(c, name, verb), c, value);
}

#endif
//
// EOF
//

