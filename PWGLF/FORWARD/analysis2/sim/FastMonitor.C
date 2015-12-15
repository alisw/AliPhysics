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
/** 
 * Monitor output objects
 */
struct FastMonitor : public TObject, public TQObject 
{
  /** 
   * Constructor 
   * 
   * 
   * @return 
   */
  FastMonitor(TSelector* s=0)
    : TObject(),
      TQObject(),
      fName("FastMonitor"),
      fCanvas(0),
      fSelector(s)
  {
    Construct();
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
    Construct();
  }
  FastMonitor& operator=(const FastMonitor& m)
  {
    if (&m == this) return *this;
    fSelector = m.fSelector;
    if (fCanvas) delete fCanvas;
    Construct();
    return *this;
  }
protected:
  void Construct()
  {
    if (gROOT->IsBatch()) {
      Warning("FastMonitor", "Batch processing, no monitoring");
      return;
    }
    if (gProof) {
      fName = gProof->GetSessionTag();
      gDirectory->Add(this);
      Bool_t ret = gProof->Connect("Feedback(TList *objs)",
				   "FastMonitor", this, 
				   "Feedback(TList *objs)");
      if (!ret) {
	Warning("FastMonitor", "Failed to connect to Proof");
	return;
      }
    }
    else if (!fSelector) return;

    if (!fCanvas) {
      fCanvas = new TCanvas(fName, Form("Monitor %s", fName.Data()), 1000, 800);
      fCanvas->SetFillColor(0);
      fCanvas->SetFillStyle(0);
      fCanvas->SetTopMargin(0.01);
      fCanvas->SetRightMargin(0.01);

      fCanvas->Divide(3,2);
      RegisterDraw(1, "histograms/type",            "", 0);
      RegisterDraw(2, "histograms/b",               "", 0);
      RegisterDraw(3, "histograms/cent",            "", 0);
      RegisterDraw(4, "histograms/dNdeta",          "", 0x8);
      RegisterDraw(5, "estimators/rawV0M",          "", 0x2);
      RegisterDraw(6, "estimators/rawRefMult00d80", "", 0x2);
    }
  }
public:
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
  }
  /** 
   * Desctructor 
   */
  virtual ~FastMonitor() 
  {
    if (!gProof) return;
    gProof->Disconnect("Feedback(TList *objs)",this, 
		       "Feedback(TList* objs)");
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
    // if (!ret) l->ls();
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

    // objs->ls();
    TList* l = static_cast<TList*>(objs->FindObject("list"));
    if (!l) {
      Warning("Feedback", "No list");
      return;
    }
    TList* hs = static_cast<TList*>(l->FindObject("histograms"));
    Int_t nEvents = 1;
    TObject* oIpz = hs->FindObject("ipZ");
    if (oIpz && oIpz->IsA()->InheritsFrom(TH1::Class())) 
      nEvents = static_cast<TH1*>(oIpz)->GetEntries();
    else 
      Warning("Feedback", "Histogram ipZ not found");

    Int_t        iPad = 1;
    TVirtualPad* p    = 0;
    while ((p = fCanvas->GetPad(iPad))) {
      TObject* o = FindPadObject(p->GetName(), l);
      if (!o) {
	Warning("Feedback", "Object correspondig to pad %s (%d) not found",
		p->GetName(), iPad);
	iPad++;
      }
      p->cd();
      if (o->IsA()->InheritsFrom(TH1::Class())) {
	TH1* h = static_cast<TH1*>(o);
	TH1* c = h->DrawCopy(p->GetTitle());
	c->SetDirectory(0);
	c->SetBit(TObject::kCanDelete);
	if (p->TestBit(BIT(15))) {
	  // Info("Feedback", "Scaling %s by 1./%d and width",
	  //      c->GetName(), nEvents);
	  c->Scale(1./nEvents, "width");
	}
      }
      else {
	TObject* c = o->DrawClone(p->GetTitle());
	c->SetBit(TObject::kCanDelete);
      }
      p->Modified();
      iPad++;
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
  /** Our name */
  TString fName;
  /** Our canvas */
  TCanvas* fCanvas;
  /** Possibly link to selector */
  TSelector* fSelector;
  ClassDef(FastMonitor,1);
};
