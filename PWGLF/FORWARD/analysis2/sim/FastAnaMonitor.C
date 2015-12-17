#ifndef __CINT__
# include <TQObject.h>
# include <TObject.h>
# include <TSelector.h>
# include <TCanvas.h>
# include <TROOT.h>
# include <TProof.h>
# include <TTimer.h>
# include <TH1.h>
# include <TObjArray.h>
# include <TClass.h>
#else
class TCanvas;
class TSelector;
class TVirtualPad;
class TTimer;
class TH1;
class TObjArray;
class TClass;
class TQObject;
class TProof;
#endif

//====================================================================
struct FastMonitor : public TObject, public TQObject 
{
  /** 
   * Execute a PROOF command. Short hand convinience 
   * 
   * @param cmd Command, or empty string. 
   * 
   * @return If cmd is empty, test if gProof is defined, other wise
   * result of command.
   */
  static Long_t ProofExec(const char* cmd=0)
  {
    Bool_t hasCmd = (cmd && cmd[0] != '\0');
    TString lne;
    lne.Form("gProof%s%s", (hasCmd ? "->" : ""), (hasCmd ? cmd : ""));
    Printf("FastMonitor::ProofExec: %s", lne.Data());
    return gROOT->ProcessLine(lne);
  }
  /** 
   * Prepare a draw 
   * 
   * @param c      List to add to 
   * @param name   Name (path) of object 
   * @param title  Drawing options  
   * @param flags  Flags 
   */
  static void PrepDraw(TCollection* c,
		       const char* name,
		       const char* title,
		       UInt_t      flags)
  {
    TNamed* n = new TNamed(name, title);
    n->SetUniqueID(flags);
    c->Add(n);
  }
  /** 
   * Constructor 
   * 
   * 
   * @return 
   */
  FastMonitor(TSelector* s=0,
		 const TString& name="FastMonitor",
		 TCollection* names=0)
    : fName(name),
      fCanvas(0),
      fSelector(s),
      fNPads(0)
  {
    Construct(names);
  }
  /** 
   * Copy constructor 
   */
  FastMonitor(const FastMonitor& m)
    : TObject(m),
      TQObject(),
      fCanvas(0),
      fSelector(m.fSelector)
  {
    Construct(0);
  }
  /** 
   * Desctructor 
   */
  virtual ~FastMonitor() 
  {
    if (ProofExec() == 0) return;
    ProofExec(Form("Disconnect(\"Feedback(TList *objs)\","
		   "%p,\"Feedback(TList* objs)\"", this));
  }
  /** 
   * Assignment operator 
   * 
   * @param m Object to assign from 
   * 
   * @return Reference to this 
   */
  FastMonitor& operator=(const FastMonitor& m)
  {
    if (&m == this) return *this;
    fSelector = m.fSelector;
    if (fCanvas) delete fCanvas;
    Construct();
    return *this;
  }
  /** 
   * Setup connection, canvas, and registered object
   * 
   * @param names Names of objects to monitor 
   */
  void Construct(TCollection* names)
  {
    if (!names) return;
    
    if (gROOT->IsBatch()) {
      Warning("FastMonitor", "Batch processing, no monitoring");
      return;
    }
    if (ProofExec()) {
      // We're on Proof
      gROOT->ProcessLine(Form("((FastMontor*)%p)->SetName("
			      "gProot->GetSessionTag())", this));
      Long_t ret = ProofExec(Form("Connect(\"Feedback(TList *objs)\","
				  "        \"FastMonitor\",%p,"
				  "        \Feedback(TList *objs)\")", this));
      if (!ret) {
	Warning("FastMonitor", "Failed to connect to Proof");
	return;
      }	
    }
    if (fCanvas) return;
    
    fCanvas = new TCanvas(fName, Form("Monitor %s", fName.Data()), 1000, 800);
    fCanvas->SetFillColor(0);
    fCanvas->SetFillStyle(0);
    fCanvas->SetTopMargin(0.01);
    fCanvas->SetRightMargin(0.01);
    
    Int_t nTotal = names->GetEntries();
    Int_t nRow   = Int_t(TMath::Sqrt(nTotal)+.5);
    Int_t nCol   = nRow;
    if (nCol * nRow < nTotal) nCol++;
    fNPads = nTotal;
    fCanvas->Divide(nCol,nRow);
    Info("FastMonitor","Create canvas with (%dx%d) [%d] pads",nCol,nRow,nTotal);
    TIter next(names);
    TObject* o = 0;
    Int_t    i = 1;
    while ((o = next()))
      RegisterDraw(i++, o->GetName(), o->GetTitle(), o->GetUniqueID());
  }
  /** 
   * Register a draw of a an object 
   * 
   * @param i      Pad number 
   * @param name   Name of object 
   * @param option Drawing option
   * @param flags  Flags 
   *
   *  - 0x1   Log(x)
   *  - 0x2   Log(y)
   *  - 0x4   Log(z)
   *  - 0x8   Scale to events and bin width 
   */
  void RegisterDraw(Int_t i,
		    const char* name,
		    const char* option,
		    UShort_t    flags=0)
  {
    TVirtualPad* p = fCanvas->GetPad(i);
    if (!p) {
      Warning("RegisterDraw", "Not enough sub-pads (%d)", i);
      return;
    }
    Info("RegisterDraw", "Adding draw # %d %s [%s] (0x%x)",
	 i, name, option, flags);
    p->SetFillColor(0);
    p->SetFillStyle(0);
    p->SetTopMargin(0.01);
    p->SetRightMargin(0.01);
    p->SetName(Form("p_%s", name));
    p->SetTitle(option);
    if (flags & 0x1) p->SetLogx();
    if (flags & 0x2) p->SetLogy();
    if (flags & 0x4) p->SetLogz();
    if (flags & 0x8) p->SetBit(BIT(15));

    fCanvas->Modified();
  }
  /** 
   * Set name of this object 
   * 
   * @param name Name 
   */
  void SetName(const char* name) { fName = name; }
  /** 
   * Get the name of this object 
   * 
   * @return Name 
   */
  const char* GetName() const { return fName.Data(); }
  /** 
   * Find pad corresponding to an object
   * 
   * @param name Name of object 
   * 
   * @return Pointer to pad or null
   */
  TVirtualPad* FindPad(const TString& name)
  {
    TVirtualPad* p = 0;
    Int_t        i = 1;
    TString      t = Form("p_%s", name.Data());
    while ((p = fCanvas->GetPad(i))) {
      if (t.EqualTo(p->GetName())) return p;
      i++;
    }
    return 0;
  }
  /** 
   * Find an object in the list @a l which corresponds to a registered
   * pad.
   * 
   * @param padName Pad's name 
   * @param l       Input collection 
   * 
   * @return The found object, or null
   */
  TObject* FindPadObject(const Char_t* padName, TCollection* l)
  {
    TString path(padName);
    path.Remove(0,2);
    if (path.Index("/") == kNPOS)
      return l->FindObject(path);
    TObjArray*   tokens  = path.Tokenize("/");
    Int_t        nTokens = tokens->GetEntriesFast();
    TCollection* current = l;
    TObject*     ret     = 0;
    for (Int_t i = 0; i < nTokens; i++) {
      TObject* o = current->FindObject(tokens->At(i)->GetName());
      if (!o) break;
      if (i == nTokens - 1) {
	ret = o;
	break;
      }
      if (!o->IsA()->InheritsFrom(TCollection::Class())) {
	Warning("FindPadObject", "Path object %s of %s is not a collection "
		"but a %s", o->GetName(), path.Data(), o->ClassName());
	break;
      }
      current = static_cast<TCollection*>(o);
    }
    delete tokens;
    if (!ret) l->ls();
    return ret;
  }
  /** 
   * Called when we get notified of 
   * 
   * @param objs List of monitored objects
   */
  void Feedback(TList* objs)
  {
    // Info("FeedBack", "List is %p", objs);
    // if (objs) objs->ls();
    if (!fCanvas) return;

    if (!objs) {
      Warning("Feedback", "No list");
      return;
    }
    // Info("FeedBack", "Looping over %d pads", fNPads);	 
    for (Int_t iPad = 1; iPad <= fNPads; iPad++) {
      TVirtualPad* p = fCanvas->cd(iPad);
      // Info("Feedback", "Drawing in sub-pad # %d: %s", iPad, p->GetName());
      TObject* o = FindPadObject(p->GetName(), objs);
      if (!o) {
	// Warning("Feedback", "Object correspondig to pad %s (%d) not found",
	//         p->GetName(), iPad);
	iPad++;
	continue; 
      }
      DrawObject(o, p->GetTitle(), p->TestBit(BIT(15)), false);
      p->cd();
      p->Modified();
      // iPad++;
    }
    fCanvas->Modified();
    fCanvas->Update();
    fCanvas->cd();
  }
  /** 
   * Draw an object.
   * 
   * @param o 
   * @param same 
   */
  void DrawObject(TObject* o, Option_t* opt, Bool_t scale, Bool_t same=false)
  {
    if (o->IsA()->InheritsFrom(TH1::Class())) {
      TH1* h = static_cast<TH1*>(o);
      TH1* c = static_cast<TH1*>(h->Clone(Form("cpy_%s", h->GetName())));
      c->SetDirectory(0);
      if (scale) {
	// Info("Feedback", "Scaling %s by 1./%d and width",
	//      c->GetName(), nEvents);	
	Int_t nEvents = c->GetBinContent(0);
	if (nEvents <= 0) return;
	c->Scale(1./nEvents, "width");
	c->SetMinimum(0);
      }
      if (p->TestBit(BIT(16)))
	// suppress stats
	c->SetStats(0);
      c->Draw(Form("%s %s",opt, (same ? "same" : "")));
      c->SetBit(TObject::kCanDelete);
      // Info("DrawObject","Drawing histogram '%s' with \"%s %s\"",
      //      c->GetName(), opt, (same ? "same" : ""));
    }
    else if (o->IsA()->InheritsFrom(TGraph::Class())) {
      TGraph* g = static_cast<TGraph*>(o);
      TGraph* c = static_cast<TGraph*>(g->Clone(Form("cpy_%s",g->GetName())));
      c->Draw(Form("%s %s",opt, (same ? "" : "a")));
      c->SetBit(TObject::kCanDelete);
      // Info("DrawObject","Drawing Graph '%s' with \"%s %s\"",
      //      c->GetName(), opt, (same ? "" : "a"));
    }
    else if (o->IsA()->InheritsFrom(TCollection::Class())) {
      TCollection* c = static_cast<TCollection*>(o);
      TIter        n(c);
      TObject*     co = 0;
      Bool_t       first = true;
      // Info("DrawObject","Drawing collection '%s' with \"%s\"",
      //      c->GetName(), opt);
      while ((co = n())) {
	DrawObject(co, opt, scale, !first);
	first = false;
      }
    }
    else {
      TObject* c = o->DrawClone(opt);
      c->SetBit(TObject::kCanDelete);
    }
  }
  /** 
   * Function to handle connect signals 
   * 
   */
  void Handle()
  {
    HandleTimer(0);
  }
  /**
   * Function to handle timer events 
   */
  Bool_t HandleTimer(TTimer*)
  {
    // Info("HandleTimer", "Selector=%p", fSelector);
    if (!fSelector) return false;
    Feedback(fSelector->GetOutputList());
    return true;
  }
  /** Our name */
  TString fName;
  /** Our canvas */
  TCanvas* fCanvas;
  /** Possibly link to selector */
  TSelector* fSelector;
  Int_t fNPads;
  ClassDef(FastMonitor,1);
};
