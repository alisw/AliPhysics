#ifndef FASTMONITOR_C
#define FASTMONITOR_C
#ifndef __CINT__
# include <TQObject.h>
# include <TObject.h>
# include <TSelector.h>
# include <TCanvas.h>
# include <TROOT.h>
// # include <TProof.h>
# include <TGraph.h>
# include <TMath.h>
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
// class TProof;
#endif

//====================================================================
struct FastMonitor : public TObject, public TQObject 
{
  /** 
   * Flags to pass 
   */
  enum {
    kLogx    = 0x01,
    kLogy    = 0x02,
    kLogz    = 0x04,
    kScale   = 0x08,
    kNoStats = 0x10
  };
  /** 
   * Internal bits set on pads 
   */
  enum {
    kBitScale   = (1<<17), // BIT(15),
    kBitNoStats = (1<<16)  // BIT(16)
  };
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
    // Printf("FastMonitor::ProofExec: %s", lne.Data());
    return gROOT->ProcessLine(lne);
  }
  FastMonitor()
    : TObject(),
      TQObject(),
      fName("FastMonitor"),
      fPaths(),
      fCanvas(0),
      fSelector(0),
      fNPads(0),
      fTimer(0)
  {
    fPaths.SetOwner();
  }
  /** 
   * Constructor 
   * 
   * 
   * @return 
   */
  FastMonitor(TSelector* s, const TString& name)
    : TObject(),
      TQObject(),
      fName(name),
      fPaths(),
      fCanvas(0),
      fSelector(s),
      fNPads(0),
      fTimer(0)
  {
    fPaths.SetOwner();
  }
  /** 
   * Copy constructor 
   */
  FastMonitor(const FastMonitor& m)
    : TObject(m),
      TQObject(),
      fPaths(),
      fCanvas(0),
      fSelector(m.fSelector),
      fTimer(0)
  {
    fPaths.SetOwner();
    TIter next(&m.fPaths);
    TObject* obj = 0;
    while ((obj = next())) Register(obj);
  }
  /** 
   * Desctructor 
   */
  virtual ~FastMonitor() { Disconnect(); }
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
    fPaths.Clear();
    TIter next(&m.fPaths);
    TObject* obj = 0;
    while ((obj = next())) Register(obj);
    return *this;
  }
  /** 
   * Connect to the source 
   * 
   */
  void Connect(Int_t freq=-1)
  {
    if (gROOT->IsBatch()) {
      Warning("FastMonitor", "Batch processing, no monitoring");
      return;
    }
    if (freq <= 0) return;
    if (ProofExec()) {
      // We're on Proof
      // Info("Construct", "Attaching to PROOF");
      gROOT->ProcessLine(Form("((FastMonitor*)%p)->SetName("
			      "gProof->GetSessionTag())", this));
      Long_t ret = ProofExec(Form("Connect(\"Feedback(TList *objs)\","
				  "        \"FastMonitor\",(void*)%p,"
				  "        \"Feedback(TList *objs)\")", this));
      if (!ret) {
	Warning("FastMonitor", "Failed to connect to Proof");
	return;
      }	
      ProofExec(Form("SetParameter(\"PROOF_FeedbackPeriod\",%d)",
		     freq*1000));
    }
    else {
      fTimer = new TTimer(freq*1000);
      fTimer->Connect("Timeout()","FastMonitor",this, "Handle()");
      fTimer->Start(-1,false);
    }
    // Info("Connect", "Monitor %s connected", fName.Data());
  }
  /** 
   * Disconnect from the source 
   * 
   */
  void Disconnect()
  {
    if (ProofExec())
      ProofExec(Form("Disconnect(\"Feedback(TList *objs)\","
		     "(void*)%p,\"Feedback(TList* objs)\"", this));
    else if (fTimer)
      fTimer->Stop();
  }
  /** 
   * Register an object.  Note the object passed here is a descripter,
   * that gives the name (descr->GetName()) and options
   * (descr->GetTitle()) and the pad options (descr->GetUniqueID()) -
   * not the actual object to draw.
   * 
   * @param descr 
   * @param proof  IF true, register for Proof(Lite)
   */
  void Register(TObject* descr, Bool_t proof=true)
  {
    Register(descr->GetName(), descr->GetTitle(), descr->GetUniqueID(), proof);
  }
  /** 
   * Prepare a draw 
   * 
   * @param name   Name (path) of object 
   * @param title  Drawing options  
   * @param flags  Flags 
   * @param proof  Register for proof 
   */
  void Register(const char* name,
		const char* title="",
		UInt_t      flags=0,
		Bool_t      proof=true)
  {
    TNamed* n = new TNamed(name, title);
    n->SetUniqueID(flags);
    fPaths.Add(n);
    if (ProofExec() && proof) {
      ProofExec(Form("AddFeedback(\"%s\")", name));
    }    
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
   * Called when we get notified of 
   * 
   * @param objs List of monitored objects
   */
  void Feedback(TList* objs)
  {
    // Info("FeedBack", "List is %p", objs);
    // if (objs) objs->ls();
    if (!fCanvas && !SetupCanvas()) return;
    
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
	Warning("Feedback", "Object correspondig to pad %s (%d) not found",
	         p->GetName(), iPad);
	iPad++;
	continue; 
      }
      DrawObject(o,
		 p->GetTitle(),
		 p->TestBit(kBitScale),
		 p->TestBit(kBitNoStats),
		 false);
      p->cd();
      p->Modified();
      // iPad++;
    }
    fCanvas->Modified();
    fCanvas->Update();
    fCanvas->cd();
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
protected:
  /** 
   * Setup canvas, and registered object
   */
  Bool_t SetupCanvas()
  {
    if (gROOT->IsBatch()) {
      Warning("FastMonitor", "Batch processing, no monitoring");
      return false;
    }
    if (fCanvas) return true;
    
    // Info("SetupCanvas", "Creating canvas");
    fCanvas = new TCanvas(fName, Form("Monitor %s", fName.Data()), 1000, 800);
    fCanvas->SetFillColor(0);
    fCanvas->SetFillStyle(0);
    fCanvas->SetTopMargin(0.01);
    fCanvas->SetRightMargin(0.01);
    
    Int_t nTotal = fPaths.GetEntries();
    Int_t nRow   = Int_t(TMath::Sqrt(nTotal)+.5);
    Int_t nCol   = nRow;
    if (nCol * nRow < nTotal) nCol++;
    fNPads = nTotal;
    fCanvas->Divide(nCol,nRow);
    Info("FastMonitor","Create canvas with (%dx%d) [%d] pads",nCol,nRow,nTotal);
    TIter next(&fPaths);
    TObject* o = 0;
    Int_t    i = 1;
    while ((o = next()))
      SetupDraw(i++, o->GetName(), o->GetTitle(), o->GetUniqueID());

    return true;
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
  void SetupDraw(Int_t i,
		 const char* name,
		 const char* option,
		 UInt_t    flags=0)
  {
    TVirtualPad* p = fCanvas->GetPad(i);
    if (!p) {
      Warning("RegisterDraw", "Not enough sub-pads (%d)", i);
      return;
    }
    //Info("RegisterDraw",
    //     "Adding draw # %d %s [%s] (0x%x)",i,name,option,flags);
    p->SetFillColor(0);
    p->SetFillStyle(0);
    p->SetTopMargin(0.01);
    p->SetRightMargin(0.01);
    p->SetName(Form("p_%s", name));
    p->SetTitle(option);
    if (flags & kLogx)    p->SetLogx();
    if (flags & kLogy)    p->SetLogy();
    if (flags & kLogz)    p->SetLogz();
    if (flags & kScale)   p->SetBit(kBitScale);
    if (flags & kNoStats) p->SetBit(kBitNoStats);

    fCanvas->Modified();
  }
  /** 
   * Draw an object.
   * 
   * @param o       Object to draw 
   * @param opt     Options 
   * @param scale   Scalar 
   * @param nostats Do not show stats 
   * @param same    Draw on current canvas 
   */
  void DrawObject(TObject*  o,
		  Option_t* opt,
		  Bool_t    scale,
		  Bool_t    nostats,
		  Bool_t    same=false)
  {
    // Info("DrawObject","Drawing %s '%s' with \"%s %s\" (%s with%s stats)",
    //      o->ClassName(), o->GetName(), opt, (same ? "same" : ""),
    //      scale ? "scaled" : "raw", nostats ? "out" : "");
    // Printf("M: Draw object: %s (%s)", o->GetName(), o->ClassName());
    if (o->IsA()->InheritsFrom(TH1::Class())) {
      TH1* h = static_cast<TH1*>(o);
      TH1* c = static_cast<TH1*>(h->Clone(Form("cpy_%s", h->GetName())));
      // Printf("M:  Is a histogram");
      c->SetDirectory(0);
      if (scale) {
	Int_t nEvents = c->GetBinContent(0);
	// Info("Feedback", "Scaling %s by 1./%d and width",
	//      c->GetName(), nEvents);	
	if (nEvents <= 0) return;
	c->Scale(1./nEvents, "width");
	c->SetMinimum(0);
      }
      if (nostats)
	// suppress stats
	c->SetStats(0);
      c->Draw(Form("%s %s",opt, (same ? "same" : "")));
      c->SetBit(TObject::kCanDelete);
    }
    else if (o->IsA()->InheritsFrom(TGraph::Class())) {
      TGraph* g = static_cast<TGraph*>(o);
      TGraph* c = static_cast<TGraph*>(g->Clone(Form("cpy_%s",g->GetName())));
      // Printf("M:  Is a graph");
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
      // Printf("M:  Is a collection");
      // Info("DrawObject","Drawing collection '%s' with \"%s\"",
      //      c->GetName(), opt);
      while ((co = n())) {
	DrawObject(co, opt, scale, nostats, !first);
	first = false;
      }
    }
    else {
      Printf("M:  Is a generic object");
      TObject* c = o->DrawClone(opt);
      c->SetBit(TObject::kCanDelete);
    }
  }
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
    // Info("FindPadObject", "pad=%s -> %p", padName, l);
    delete tokens;
    if (!ret) l->ls();
    return ret;
  }
  /** Our name */
  TString fName;
  /** List of things to monitor */
  TList fPaths;
  /** Our canvas */
  TCanvas* fCanvas;
  /** Possibly link to selector */
  TSelector* fSelector;
  /** Number of pads registered */
  Int_t fNPads;
  /** Possible timer */
  TTimer* fTimer;
  ClassDef(FastMonitor,1);
};
#endif
//
// EOF
//
