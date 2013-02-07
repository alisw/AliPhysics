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
# else
class THStack;
class TH1;
class TH2;
class TCollection;
class TCanvas;
class TVirtualPad;
class TPad;
class TLatex;
class TAxis;
# endif


class SummaryDrawer 
{
public:
  SummaryDrawer() 
    : fCanvas(0), 
      fTop(0), 
      fBody(0),
      fHeader(0),
      fParName(0),
      fParVal(0),
      fPause(false),
      fLandscape(false), 
      fRingMap(0)
  {
    fRingMap = new TVirtualPad*[6];
    fRingMap[0] = 0;
    fRingMap[1] = 0;
    fRingMap[2] = 0;
    fRingMap[3] = 0;
    fRingMap[4] = 0;
    fRingMap[5] = 0;
  }
  

protected:
  //__________________________________________________________________
  /** 
   * Find an object in a collection
   * 
   * @param parent Parent list
   * @param name   Name of object
   * @param verb   Be verbose 
   * 
   * @return Pointer to object or null 
   */
  static TObject* GetObject(const TCollection* parent, 
			    const TString&     name,
			    Bool_t             verb=true)
  {
    if (!parent) {
      if (verb) Warning("GetObject", "No parent list");
      return 0;
    }
    if (name.IsNull()) { 
      if (verb) Warning("GetObject", "No name specified");
      return 0;
    }
    TObject* o = parent->FindObject(name);
    if (!o) {
      if (verb) Warning("GetObject", "Object \"%s\" not found in parent \"%s\"",
			name.Data(), parent->GetName());
      return 0;
    }
    return o;
  }
  //__________________________________________________________________
  /** 
   * Find an object in a directory
   * 
   * @param parent Parent directory
   * @param name   Name of object
   * @param verb   Be verbose 
   * 
   * @return Pointer to object or null 
   */
  static TObject* GetObject(const TDirectory* parent, 
			    const TString& name, 
			    Bool_t verb=true)
  {
    if (!parent) {
      if (verb) Warning("GetObject", "No parent directory");
      return 0;
    }
    if (name.IsNull()) { 
      if (verb) Warning("GetObject", "No name specified");
      return 0;
    }
    TObject* o = const_cast<TDirectory*>(parent)->Get(name);
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
  //___________________________________________________________________
  /** 
   * Get a UShort_t parameter value 
   * 
   * @param c      Parent collection
   * @param name   Name of parameter
   * @param value  On return the value
   * @param verb   If true, complain if not found 
   */
  static void GetParameter(const TCollection*  c, 
			   const TString&      name, 
			   UShort_t&           value,
			   Bool_t              verb=true)
  {
    TObject* o = GetObject(c, name, verb);
    if (!o) return;
    value = o->GetUniqueID();
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
  static void GetParameter(const TCollection*  c, 
			   const TString&      name, 
			   Int_t&              value,
			   Bool_t              verb=true)

  {
    TObject* o = GetObject(c, name, verb);
    if (!o) return;
    value = o->GetUniqueID();
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
  static void GetParameter(const TCollection*  c, 
			   const TString&      name, 
			   Double_t&           value,
			   Bool_t              verb=true)

  {
    TObject* o = GetObject(c, name, verb);
    if (!o) return;
    UInt_t  i = o->GetUniqueID();
    Float_t v = *reinterpret_cast<Float_t*>(&i);
    value = v;
  }
  //_____________________________________________________________________
  /** 
   * Get a Bool_t parameter value 
   * 
   * @param c      Parent collection
   * @param name   Name of parameter
   * @param value  On return the value
   * @param verb   If true, complain if not found 
   */
  static void GetParameter(const TCollection*  c, 
			   const TString&      name, 
			   Bool_t&             value,
			   Bool_t              verb=true)
  {
    TObject* o = GetObject(c, name, verb);
    if (!o) return;
    value = o->GetUniqueID();
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
  static TCollection* GetCollection(const TCollection* parent, 
				    const TString&     name,
				    Bool_t             verb=true)
  {
    TObject* o = GetObject(parent, name, verb);
    if (!o) return 0;
    if (!CheckType(o, TCollection::Class(), parent->GetName())) return 0;
    return static_cast<TCollection*>(o);
  }
  //____________________________________________________________________
  /** 
   * Find a collection in a directory
   * 
   * @param parent Parent directory
   * @param name   Name of the collection 
   * @param verb   If true and not found, complain
   * 
   * @return pointer to collection on success, otherwise null 
   */
  static TCollection* GetCollection(const TDirectory* parent, 
				    const TString&    name,
				    Bool_t            verb=true)
  {
    TObject* o = GetObject(parent, name, verb);
    if (!o) return 0;
    if (!CheckType(o, TCollection::Class(), parent->GetName())) return 0;
    return static_cast<TCollection*>(o);
  }

  //____________________________________________________________________
  /** 
   * Get a 1D histogram from a collection
   * 
   * @param parent Parent collection 
   * @param name   Name of histogram 
   * @param verb   If true and not found, complain
   * 
   * @return pointer or null
   */
  static TH1* GetH1(const TCollection* parent, 
		    const TString&     name,
		    Bool_t             verb=true)
  {
    TObject* o = GetObject(parent, name, verb);
    if (!o) return 0;
    if (!CheckType(o, TH1::Class(), parent->GetName())) return 0;
    return static_cast<TH1*>(o);
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
  static TH2* GetH2(const TCollection* parent, 
		    const TString&     name, 
		    Bool_t             verb=true)
  {
    // Info("GetH2", "Getting 2D histogram of %s from %p", name.Data(), c);
    // --- Find the object -------------------------------------------
    TObject* o = GetObject(parent, name, verb);
    if (!o) return 0;

    // --- Check the type of object ----------------------------------
    if (!CheckType(o, TH2::Class(), parent->GetName())) return 0;
  
    // --- Return the collection -------------------------------------
    return static_cast<TH2*>(o);
  }
  //__________________________________________________________________
  /** 
   * Get a histogram stack from a collection
   * 
   * @param parent Parent collection 
   * @param name   Name of histogram 
   * @param sub    Sub-component 
   * @param verb   If true and not found, complain
   * 
   * @return pointer or null
   */
  static THStack* GetStack(const TCollection* parent, 
			   const TString&     name,
			   const char*        sub=0,
			   Bool_t             verb=true)
  {
    // --- Find the object -------------------------------------------
    TObject* o = GetObject(parent, name, verb);
    if (!o) return 0;

    // --- Check the type of object -----------------------------------
    if (!CheckType(o, THStack::Class(), parent->GetName())) return 0;
  
    THStack* stack = static_cast<THStack*>(o);
    if (sub == 0) return stack;
  
    if (stack->GetHists()->GetEntries() <= 0 ||stack->GetMaximum() < 1) { 
      // Info("GetStack", "No entries in %s", name.Data());
      stack->GetHists()->Delete();
      const char* subs[] = { "FMD1I", "FMD2I", "FMD2O", "FMD3O", "FMD3I", 0 };
      const char** ptr   = subs;
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
    fTop->Clear();
    fTop->SetNumber(1);
    fTop->SetFillColor(kBlue-5);
    fTop->SetBorderSize(0);
    fTop->SetBorderMode(0);

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
   * @param pname Name of PDF file to make 
   * @param landscape if true, print in landscape 
   *
   * @return Created canvas 
   */
  void CreateCanvas(const TString& pname, Bool_t landscape=false)
  {
    // Info("CreateCanvas", "Creating canvas");
    fLandscape = landscape;
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

    Float_t dy = .05;
    fTop = new TPad("top", "Top", 0, 1-dy, 1, 1, 0, 0);
    fTop->SetNumber(1);
    fTop->SetFillColor(kBlue-5);
    fTop->SetBorderSize(0);
    fTop->SetBorderMode(0);
    fCanvas->cd();
    fTop->Draw();
    
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
    
    fParVal = new TLatex(x2, y, "");
    fParVal->SetTextAlign(13);
    fParVal->SetNDC();
    fParVal->SetTextSize(s);

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
    fCanvas->Print(Form("%s]", fCanvas->GetTitle()),
		   Form("pdf %s", fLandscape ? "Landscape" : ""));
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
    TString tit;
    tit.Form("pdf %s Title:%s", fLandscape ? "Landscape" : "",
	     title.Data());

    fTop->cd();
    fHeader->SetTextSize(size);
    fHeader->DrawLatex(.5,.5,title);
  
    fCanvas->Modified();
    fCanvas->Update();
    fCanvas->cd();

    gSystem->RedirectOutput("/dev/null");
    fCanvas->Print(fCanvas->GetTitle(), tit);
    gSystem->RedirectOutput(0);
  
    Pause();

    ClearCanvas();
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
   */
  void DrawInPad(TVirtualPad* c, 
		 Int_t        padNo, 
		 TObject*     h, 
		 Option_t*    opts="",
		 UShort_t     flags=0x0)
  {
    TVirtualPad* p = c->GetPad(padNo);
    if (!p) { 
      Warning("DrawInPad", "Pad # %d not found in %s", padNo, c->GetName());
      return;
    }
    DrawInPad(p, h, opts, flags);
  }
  //__________________________________________________________________
  /** 
   * Draw an object in pad 
   * 
   * @param p       Pad
   * @param h       Object to draw 
   * @param opts    Options
   * @param flags   Flags
   */
  void DrawInPad(TVirtualPad* p, 
		 TObject*     h, 
		 Option_t*    opts="",
		 UShort_t     flags=0x0)
  {
    if (!p) { 
      Warning("DrawInPad", "No pad specified");
      return;
    }
    p->cd();
    // Info("DrawInPad", "Drawing in pad %p", p);
    // fBody->ls();
    if (flags & 0x1) p->SetLogx();
    if (flags & 0x2) p->SetLogy();
    if (flags & 0x4) p->SetLogz();
    p->SetFillColor(0);
    TString o(opts);
    if (o.Contains("colz", TString::kIgnoreCase)) 
      p->SetRightMargin(0.15);
    if (o.Contains("text", TString::kIgnoreCase)) {
      TH1* hh = static_cast<TH1*>(h);
      hh->SetMaximum(1.1*hh->GetMaximum());
      hh->SetMarkerSize(2);
      o.Append("30");
    }
    if (!h) {
      // Warning("DrawInPad", "Nothing to draw in pad # %s", p->GetName());
      return;
    }
    h->Draw(o);
    
    if (flags& 0x10) { 
      TLegend* l = p->BuildLegend(0.33, .67, .66, .99-p->GetTopMargin());
      l->SetFillColor(0);
      l->SetFillStyle(0);
      l->SetBorderSize(0);
    }
    p->Modified();
    p->Update();
    p->cd();
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
    if (flags & 0x1) p->SetLogx();
    if (flags & 0x2) p->SetLogy();
    if (flags & 0x4) p->SetLogz();
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
    TGaxis*  a1 = new TGaxis(p->GetUxmin(), p->GetUymin(), 
			     p->GetUxmin(), p->GetUymax(), 
			     0, m1, 510);
    a1->SetLineColor(h1->GetLineColor());
    a1->Draw();
    
    o.Append(" same");
    Double_t m2    = 1.1 * h2->GetMaximum();
    Double_t scale = m1 / m2;
    h2->Scale(scale);
    h2->DrawCopy(o);
    if (e3) h2->DrawCopy(fopts);
    
    TGaxis*  a2 = new TGaxis(p->GetUxmax(), p->GetUymin(), 
			     p->GetUxmax(), p->GetUymax(), 
			     0, m2, 510, "+L");
    a2->SetLineColor(h2->GetLineColor());
    a2->Draw();
    
    if (flags& 0x10) { 
      TLegend* l = p->BuildLegend();
      l->SetFillColor(0);
      l->SetFillStyle(0);
      l->SetBorderSize(0);
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
  /** 
   * Draw an object in pad 
   * 
   * @param d       Detector 
   * @param r       Ring 
   * @param h       Object to draw 
   * @param opts    Options
   * @param flags   Flags
   */
  void DrawInRingPad(UShort_t  d, 
		     Char_t    r, 
		     TObject*  h, 
		     Option_t* opts="",
		     UShort_t  flags=0x0)
  {
    TVirtualPad* p = RingPad(d, r);
    if (!p) {
      Warning("DrawInRingPad", "No pad found for FMD%d%c", d, r);
      return;
    }
    DrawInPad(p, h, opts, flags);
  }
    

  //__________________________________________________________________
  /** 
   * Pause after each plot
   * 
   */
  void Pause()
  {
    if (!fPause) return;
    printf("Press enter to continue");
    std::cin.get();
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
};
#endif

