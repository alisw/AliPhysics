/**
 * @file   FastAnalysis.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Fri Mar 20 12:13:28 2015
 * 
 * @brief This script defines classes for looping over the data
 * produced by FastSim.C
 * 
 */

#include <TSelector.h>
#include <TQObject.h>
#ifndef __CINT__
# include <TParticle.h>
# include <TMath.h>
# include <TTree.h>
# include <TFile.h>
# include <THStack.h>
# include <TH1.h>
# include <TH2.h>
# include <TClonesArray.h>
# include <TAxis.h>
# include <TCanvas.h>
# include <TStopwatch.h>
# include <TStyle.h>
# include <TLegend.h>
# include <TLegendEntry.h>
# include <TMultiGraph.h>
# include <TGraph.h>
# include <TArrayI.h>
# include <TChain.h>
# include <TParameter.h>
# include <TSystem.h>
# include <TUrl.h>
# include <TGraph.h>
// # include <TProof.h>
#else
class TTree;
class TH1;
class TH2;
class TH1D;
class THStack;
class TCanvas;
class TAxis;
class TStopwatch;
class TLegend;
class TMultiGraph;
class TParticle;
class TArrayI;
class TProof;
class TUrl;
class TVirtualPad;
#endif


//====================================================================
/** 
 * Header in the EG tree 
 * 
 */
struct ShortHeader {
  UInt_t   fRunNo;
  UInt_t   fEventId;
  UInt_t   fNtgt;
  UInt_t   fNproj;
  UInt_t   fNbin;
  UInt_t   fType;
  Double_t fIpX;
  Double_t fIpY;
  Double_t fIpZ;
  Double_t fB;
  Double_t fC;
  Double_t fPhiR;
  UInt_t   fNSpecNproj;  // # of spectator neutrons in projectile
  UInt_t   fNSpecNtgt;   // # of spectator neutrons in target 
  UInt_t   fNSpecPproj;  // # of spectator protons in projectile
  UInt_t   fNSpecPtgt;   // # of spectator protons in target
  
  void Print()
  {
    Printf(" Run #/Event:          %9d/%9d", fRunNo, fEventId);
    Printf(" Participants/binary:  %4d/%4d/%3d", fNtgt, fNproj, fNbin);
    Printf(" Event type:           %7s%12s",(fType==1?"Non":
					     fType==2?"Single":
					     "Double"), "-diffractive");
    Printf(" IP:                   (%-5.1f,%-5.1f,%-5.1f)",fIpX,fIpY,fIpZ);
    Printf(" Impact par./cent.:    (%13f/%-3d)", fB, Int_t(fC));
    Printf(" Reaction plane:       %19f", fPhiR);
    Printf(" Specs (Nt,Np,Pt,Pp):  %4d/%4d/%4d/%4d",
	   fNSpecNtgt, fNSpecNproj, fNSpecPtgt, fNSpecPproj);
  }
  void Clear(Option_t* option="")
  {
    Reset(0,0);
  }
  void Reset(UInt_t runNo, UInt_t eventNo)
  {
    fRunNo      = runNo;
    fEventId    = eventNo;
    fIpX        = 1024;
    fIpY        = 1024;
    fIpZ        = 1024;
    fNtgt       = -1;
    fNproj      = -1;
    fNbin       = -1;
    fPhiR       = -1;
    fB          = -1;
    fC          = -1;
    fNSpecNtgt  = -1;
    fNSpecNproj = -1;
    fNSpecPtgt  = -1;
    fNSpecPproj = -1;
  }
};

//====================================================================
struct FastAnaMonitor : public TObject, public TQObject 
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
    Printf("FastAnaMonitor::ProofExec: %s", lne.Data());
    return gROOT->ProcessLine(lne);
  }
  /** 
   * Constructor 
   * 
   * 
   * @return 
   */
  FastAnaMonitor(TSelector* s=0,
		 const TString& name="FastAnaMonitor",
		 TCollection* names=0)
    : fName(name),
      fCanvas(0),
      fSelector(s),
      fNPads(0)
  {
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
    Info("FastAnaMonitor", "Create canvas with (%dx%d) [%d] pads",
	 nCol, nRow, nTotal);
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
   * Desctructor 
   */
  virtual ~FastAnaMonitor() 
  {
    if (ProofExec() == 0) return;
    ProofExec(Form("Disconnect(\"Feedback(TList *objs)\","
		   "%p,\"Feedback(TList* objs)\"", this));
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
    // objs->ls();
    // fCanvas->ls();

    // Int_t        iPad = 1;
    // TVirtualPad* p    = 0;
    // while ((p = fCanvas->GetPad(iPad))) {
    // Info("FeedBack", "Looping over %d pads", fNPads);	 
    for (Int_t iPad = 1; iPad <= fNPads; iPad++) {
      TVirtualPad* p = fCanvas->cd(iPad);
      // Info("Feedback", "Drawing in sub-pad # %d: %s", iPad, p->GetName());
      TObject* o = FindPadObject(p->GetName(), objs);
      if (!o) {
	Warning("Feedback", "Object correspondig to pad %s (%d) not found",
		p->GetName(), iPad);
	// iPad++;
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
  ClassDef(FastAnaMonitor,1);
};

//====================================================================
/** 
 * Base class for processors 
 */
struct FastAnalysis : public TSelector
{
  /** Pointer to tree being analysed */
  TTree*        fTree; //! 
  /** Cache of our header */
  ShortHeader*  fHeader;     //!
  /** List of particles */
  TClonesArray* fParticles;  //!
  /** Whether to be verbose */
  Bool_t        fVerbose;    //
  /** A simple timer */
  TStopwatch    fTimer;      //!
  /** Name */
  TString       fName;       //!
  /** Cache of event discriminator */
  ULong64_t     fEventMult;
  /** The centrality method to use */
  TString       fCentMethod;
  /** The centrality histogram to use - if any */
  TH1*          fCentHist; //!
  /** Number of good events */
  ULong_t fOK;
  /** Monitor frequency in seconds */
  Int_t fMonitor;
  /** 
   * Constructor.  Opens the file passed and sets internal pointers to
   * tree, header, and particle list.
   * 
   * @param verb     Whether to be verbose 
   */  
  FastAnalysis(Bool_t verb=false, Int_t monitor=0)
    : fTree(0),
      fHeader(0),
      fParticles(0),
      fVerbose(verb),
      fEventMult(0),
      fCentMethod(""),
      fCentHist(0),
      fOK(0),
      fMonitor(monitor)
  {
  }
  /**
   * Destructor 
   */
  virtual ~FastAnalysis()
  {
    if (fHeader)    delete fHeader;
    if (fParticles) delete fParticles;
  }
  /** 
   * Set the verbosity flag 
   * 
   * @param verb If true, be verbose
   */
  void SetVerbose(Bool_t verb) { fVerbose = verb; }
  /** 
   * Set up a monitor 
   */
  void SetupMonitor()
  {
    // We disable monitoring in batch mode
    if (gROOT->IsBatch()) {
      Warning("SetupMonitor", "Batch processing, no monitoring");
      return;
    }

    // If the monitor frequency is 0 or negative, do nothing 
    if (fMonitor <= 0) {
      Info("SetupMonitor", "Monitoring frequency too low");
      return;
    }

    // Get list of names of monitor objects. 
    TList* objs = GetMonitorObjects();
    // If there's no monitored objects defined, do nothing 
    if (!objs || objs->GetEntries() < 1) {
      Info("SetupMonitor", "No monitored objects defined");
      return;
    }

    TString name("FastAnaMonitor");
    if (ProofExec() != 0) {
      Long_t ret = ProofExec("GetSessionTag()");
      name = reinterpret_cast<const char*>(ret);
    }
    
    // Create our monitor 
    FastAnaMonitor* monitor = new FastAnaMonitor(this, name, objs);
    Info("SetupMonitor", "Created monitor %s with objects",
	 monitor->GetName());
    objs->Print();
    if (ProofExec() != 0) {
      gDirectory->Add(monitor);
      Long_t ret = ProofExec(Form("Connect(\"Feedback(TList *objs)\","
				  "\"FastAnaMonitor\", (void*)%p, "
				  "\"Feedback(TList *objs)\")",monitor));
      if (!ret) {
	Warning("FastMonitor", "Failed to connect to Proof");
	delete monitor;
	return;
      }
      ProofExec(Form("SetParameter(\"PROOF_FeedbackPeriod\",%d)",
		     fMonitor*1000));
      TIter    next(objs);
      TObject* obj = 0;
      while ((obj = next()))
	ProofExec(Form("AddFeedback(\"%s\")", obj->GetName()));
    }
    else {
      TTimer* timer = new TTimer(fMonitor*1000);
      timer->Connect("Timeout()","FastAnaMonitor",monitor, "Handle()");
      // ::Info("Run", "Turning on monitoring");
      timer->Start(-1,false);
    }
  }
  /** 
   * @{ 
   * @name Selector interface 
   */
  /** 
   * Set-up our branches 
   * 
   * @return true on success, false otherwise 
   */
  Bool_t SetupBranches()
  {
    fHeader = new ShortHeader;
    fParticles = new TClonesArray("TParticle");
    fTree->SetBranchAddress("header",    &(fHeader->fRunNo));
    fTree->SetBranchAddress("particles", &fParticles);
    if (!fCentMethod.IsNull() && !fCentMethod.EqualTo("B"))
      fTree->SetBranchAddress(fCentMethod, &fEventMult);
    return true;
   
  }
  /** 
   * Set-up the estimator is one is chosen.  We get a pointer to the
   * current file, and retrieve the relevant histogram from that.
   * 
   * @return true on success, false otherwise 
   */
  virtual Bool_t SetupEstimator()
  {
    if (fCentMethod.IsNull()) return true;
    if (fCentMethod.EqualTo("B",TString::kIgnoreCase)) return true;

    if (!fTree) {
      Warning("SetupEstimator", "No tree defined!");
      return false;
    }
    TFile* f = fTree->GetCurrentFile();
    if (!f) {
      Warning("SetupEstimator", "Failed to get current file from tree");
      return false;
    }
    // TList* l = static_cast<TList*>(f->Get("estimators"));
    // if (!l) {
    //   Warning("SetupEstimator",
    //          "Failed to get list of estimators from file %s",
    // 		f->GetPath());
    //   f->ls();
    //   return false;
    // }
    // TObject* o = l->FindObject(fCentMethod);
    TObject* o = f->Get(Form("estimators/%s",fCentMethod.Data()));
    if (!o) {
      Warning("SetupEstimator", "Failed to get %s from estimator list:",
	      fCentMethod.Data());
      // l->ls();
      f->ls();
      return false;
    }
    if (!o->IsA()->InheritsFrom(TH1::Class())) {
      Warning("SetupEstimator", "Estimator object %s is not a histogram "
	      "but a %s", o->GetName(), o->ClassName());
      return false;
    }    
    fCentHist = static_cast<TH1*>(o);
    fCentHist->SetDirectory(0);
    fTree->SetBranchAddress(fCentMethod, &fEventMult);
    return true;
  }
  /** 
   * Called on initialization 
   * 
   * @param tree Tree to analyse 
   */
  void Init(TTree* tree)
  {
    Info("Init", "Initializing with tree %p (%s)",
	 tree, (tree ? tree->ClassName() : ""));
    if (!tree) return;

    TFile* file = tree->GetCurrentFile();
    Info("Init", "Current file: (%p) %s", file,
	 (file ? file->GetName() : ""));
    
    fTree = tree;
    if (!SetupBranches())
      Fatal("Init", "Failed to set-up branches");
    // if (!SetupEstimator())
    //   Fatal("Init", "Failed to set-up estimator");
  }
  virtual void Begin(TTree*)
  {
    SetupMonitor();
  }
  /** 
   * Called when the file changes 
   * 
   * @return true on success, false otherwise 
   */
  Bool_t Notify()
  {
    TFile* file = fTree->GetCurrentFile();
    Info("Notify", "processing file: (%p) %s", file,
	 (file ? file->GetName() : ""));
    if (!file) return true;
    if (!SetupBranches()) {
      Warning("Notify", "Failed to set-up branches");
      return false;
    }
    if (!SetupEstimator()) {
      Warning("Notify", "Failed to set-up estimator");
      return false;
    }
    return true;
  }
  /** 
   * Called on each event 
   * 
   * @param entry Entry in chain 
   * 
   * @return true on success, false otherwise 
   */
  virtual Bool_t Process(Long64_t entry)
  {
    Clear();
    Int_t read = fTree->GetTree()->GetEntry(entry);
    if (read <= 0) return false;
    
    if (fVerbose)
      printf("Process event %7lld (0x%x) ",
	     entry, fHeader->fType);
    if (!ProcessHeader()) {
      if (fVerbose) Printf("No accepted");      
      return true;
    }
    fOK++;

    ProcessParticles();

    return true;
  }
  /** 
   * Called at the end of the slave processing.  Puts the event count
   * into a parameter that is stored in the output list.  Users can
   * overload this to do more stuff should it be needed.
   */
  virtual void SlaveTerminate()
  {
    Printf("A total of %ld events accepted", fOK);
    TParameter<Long_t>* cnt = new TParameter<Long_t>("count",fOK,'+');
    fOutput->Add(cnt);
    fOutput->ls();
  }
  Int_t Version() const { return 2; }
  /* @} */

  /** 
   * @{ 
   * @name Some service functions 
   */
  /** 
   * Get the event centrality.  
   * 
   * If the centrality estimator is not set, we use the centrality
   * stored in the header.
   * 
   * Otherwise, we find the bin corresponding to the event
   * multiplicity (as stored in the selected branch), of the estimator
   * histogram, and return the corresponding centrality (bin content).
   * 
   * If the event multiplcity is below the range given by the
   * estimator, we return -1. If the event multiplicity is above the
   * range defined by the estimator, we return 1000.
   * 
   * @return The event centrality. 
   */
  Double_t GetCentrality() const
  {
    if (!fCentHist) return fHeader->fC;
    Int_t nBin = fCentHist->GetNbinsX();
    Int_t bin  = fCentHist->GetXaxis()->FindBin(fEventMult);
    if (bin <= 0)    return -1;
    if (bin-1 == nBin && fEventMult == fCentHist->GetXaxis()->GetXmax()) {
      bin        =  nBin;
    }
    if (bin >  nBin) return 1000;
    return fCentHist->GetBinContent(bin);
  }
  /** 
   * Get the event count from internal counter, or if that is 0 or
   * smaller, from the output
   * 
   * @return Event count 
   */
  Long_t GetEventCount()
  {
    if (fOK > 0) return fOK;

    typedef TParameter<Long_t> Param_t;
    Param_t* p = static_cast<Param_t*>(GetOutputObject("count",
						       Param_t::Class()));
    if (!p) return 0;
    return fOK = p->GetVal();
  }
  /** 
   * Get an object from the output list, possibly checking the type 
   * 
   * @param name Name of object 
   * @param cls  Possible class pointer 
   * 
   * @return Found object, or null
   */
  TObject* GetOutputObject(const char* name, TClass* cls)
  {
    TObject* o = fOutput->FindObject(name);
    if (!o) {
      Warning("GetOutputObject", "Object %s not found in output", name);
      fOutput->ls();
      return 0;
    }
    if (cls && !o->IsA()->InheritsFrom(cls)) {
      Warning("GetOutputObject", "Object %s from output is not a %s, but a %s",
	      o->GetName(), cls->GetName(), o->ClassName());
      return o;
    }
    return o;
  }
  /* @} */

  /** 
   * @{
   * @name overloadable behaviour 
   */
  /** 
   * Get the list of monitor objects.  Sub-classes should overload
   * this to return a TList of TNamed objects.  The name of each entry
   * must corresponding to the path of an object in the output list.
   * The title encodes the draw option used for the object, and the
   * unique id sets pad options, such as log scale, and scaling to
   * number of events.  If scaling to number of events is requested,
   * then the number of events is assumed to be encode in the
   * underflow bin.
   * 
   * @return A TList of monitor object identifiers, or null
   */
  virtual TList* GetMonitorObjects() { return 0; }
  /** 
   * Clear internal caches.  Called at start of each event.  Can be
   * overloaded to do some more stuff if needed.
   */
  virtual void Clear(Option_t* option="")
  {
    fHeader->Clear(option);
    fParticles->Clear();
    fEventMult = 0;
  }
  /** 
   * Process the particles.  Called once for each event.
   *
   * This in turn calls ProcessParticle for each particle. 
   *
   * By default secondaries and neutral particles are not processed.
   * If a derived class need to look at these, that class should
   * overwrite AcceptSecondaries and/or AcceptNeutrals to return true.
   */
  virtual void ProcessParticles()
  {
    Int_t nParticles = fParticles->GetEntriesFast();
    Int_t nOK        = 0;
    Int_t nSec       = 0;
    Int_t nNeut      = 0;
    Int_t nOut       = 0;
    for (Int_t iPart = 0; iPart < nParticles; iPart++) {
      TObject*  oPart = fParticles->At(iPart);
      if (!oPart) continue;      
      if (!oPart->TestBit((1<<14)) || oPart->TestBit(1<<15)) {
	// If this is a secondary, increment counter 
	nSec++;
	if (!AcceptSecondaries()) continue;
      }
      if (!oPart->TestBit((1<<16))) {
	// If this is neutral, increment counter 
	nNeut++;
	if (!AcceptNeutrals()) continue;
      }

      if (!ProcessParticle(static_cast<TParticle*>(oPart))) {
	nOut++;
	continue;
      }
      nOK++;
    }
    if (fVerbose) 
      Printf("ok:%7d sec:%7d neu=%7d out:%7d total:%7d",
	     nOK, nSec, nNeut, nOut, nParticles);

  }
  /** 
   * Whether to accept secondary particles.  By default this returns
   * false, meaning we do not process secondary particles.  A derived
   * class can overload this to return true in case one want to look
   * at secondary particles.  
   *
   * If one need finer control, this can be overwritten to return
   * true, and one can inspect bit 14 of the object bits to see if a
   * particle is a secondary.
   * 
   * @return false 
   */
  virtual Bool_t AcceptSecondaries() const { return false; }
  /** 
   * Whether to accept neutral particles.  By default this returns
   * false, meaning we do not process neutral particles.  A derived
   * class can overload this to return true in case one want to look
   * at neutral particles.
   * 
   * If one need finer control, this can be overwritten to return
   * true, and one can inspect bit 15 of the object bits to see if a
   * particle is a secondary.
   * 
   * @return false 
   */
  virtual Bool_t AcceptNeutrals() const { return false; }
  /** 
   * Process a single particle.  
   *
   * In a derived class one can inspect bit 14 (15) to test of the
   * particle is a secondary (neutral) particle.  By default secondary
   * (netrual) particles are not processes, unless the derived class
   * overwrites AcceptSecondaries (AcceptNeutrals) to return true.
   * 
   * @param p Pointer to TParticle object 
   * 
   * @return true if the particle was accepted.  
   */
  virtual Bool_t ProcessParticle(const TParticle* p) = 0;
  /** 
   * Process the header.  Shall return true if the event is accepted,
   * false otherwise.  Must be overloaded by derived class. 
   * 
   * @return True if the event is to be taken. 
   */
  virtual Bool_t ProcessHeader() = 0;
  /* @} */

  /** 
   * @{ 
   * @name Static interface for running 
   */
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
    Printf("FastAnalysis::ProofExec: %s", lne.Data());
    return gROOT->ProcessLine(lne);
  }
  /** 
   * Extract key value pair from string 
   * 
   * @param in  Input string 
   * @param key On return, the key 
   * @param val On return, the value
   * @param sep Separator between key an value 
   * 
   * @return false of separator isn't found in input 
   */
  static Bool_t Str2KeyVal(const TString& in,
			   TString&       key,
			   TString&       val,
			   const char     sep='=')
  {
    Int_t idx = in.Index(sep);
    if (idx == kNPOS) return false;

    key = in(0,idx);
    val = in(idx+1, in.Length()-idx-1);
    return true;
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
  /** 
   * Set-up PROOF 
   * 
   * @return true on success, false otherwise 
   */
  static Bool_t SetupProof(TUrl& url, const char* opt)
  {
    Long_t ret = 0;
    gROOT->LoadClass("TProof", "libProof");
    ret = gROOT->ProcessLine(Form("TProof::Reset(\"%s\")",url.GetUrl()));
    ret = gROOT->ProcessLine(Form("TProof::Open(\"%s\")",url.GetUrl()));
    if (!ret) {
      Printf("Error: FastAnalysis::SetupProof: Failed to connect");
      return false;
    }
    ret = ProofExec("ClearCache()");

    TString phy = gSystem->ExpandPathName("$(ALICE_PHYSICS)");
    TString ali = gSystem->ExpandPathName("$(ALICE_ROOT)");
    TString fwd = phy + "/PWGLF/FORWARD/analysis2";
    TString mkLib = gSystem->GetMakeSharedLib();
    mkLib.ReplaceAll("-std=c++14", "-std=c++98");

    ProofExec(Form("Exec(\"gSystem->SetMakeSharedLib(\\\"%s\\\")\")",
		   mkLib.Data()));

    ret = ProofExec(Form("Load(\"%s/sim/FastAnalysis.C+%s\",true)",
			 fwd.Data(), opt));
    if (ret != 0) {
      Printf("Error: FastAnalysis::SetupProof: Failed to load");
      return false;
    }
    return true;
  }
  
  /** 
   * Run a selector 
   * 
   * @param url Url to process 
   * 
   * @return true on success, false on error
   */
  static Bool_t Run(const char* url, const char* output, const char* opt="g")
  {
    Long64_t     nev     = -1;
    TString      type    = "INEL";
    Int_t        monitor = -1;
    Bool_t       verbose = false;
    TUrl         u(url);
    TString      uout    = "";
    TObjArray*   opts    = TString(u.GetOptions()).Tokenize("&");
    TObjString*  token   = 0;
    TIter        nextToken(opts);
    while ((token = static_cast<TObjString*>(nextToken()))) {
      TString& str = token->String();
      if (str.IsNull()) continue;

      if (str.EqualTo("verbose")) {
	verbose = true;
	continue;
      }
      TString  key, val;
      if (!Str2KeyVal(str,key,val)) {
	if (!uout.IsNull()) uout.Append("&");
	uout.Append(str);
	continue;
      }
      
      if      (key.EqualTo("events")) nev     = val.Atoll();
      else if (key.EqualTo("type"))   type    = val;
      else if (key.EqualTo("monitor"))monitor = val.Atoi();
      else {
	if (!uout.IsNull()) uout.Append("&");
	uout.Append(str);
      }
    }
    opts->Delete();
    u.SetOptions(uout);
    if (!u.IsValid()) {
      Printf("Error: FastAnalysis::Run: URL %s is invalid", u.GetUrl());
      return false;
    }
    FastAnalysis* analysis = Make(type,verbose,monitor);
    if (!analysis) return false;

    TString treeName = u.GetAnchor();
    if (treeName.IsNull()) treeName = "T";
    TFile*  file = TFile::Open(u.GetFile(), "READ");
    if (!file) {
      Printf("Error: FastAnalysis::Run: Failed to open %s",
	     u.GetFile());
      return false;
    }

    TChain*   chain = new TChain(treeName, treeName);
    TObject*  o = file->Get(treeName);
    if (!o) {
      Printf("Error: FastAnalysis::Run: Couldn't get %s from %s",
	     treeName.Data(), u.GetFile());
      file->Close();
      return false;
    }
    Int_t cret = 0;
    if (o->IsA()->InheritsFrom(TChain::Class()))
      cret = chain->Add(static_cast<TChain*>(o));
    else if (o->IsA()->InheritsFrom(TTree::Class()))
      cret = chain->AddFile(u.GetFile());
    else if (o->IsA()->InheritsFrom(TCollection::Class()))
      cret = chain->AddFileInfoList(static_cast<TCollection*>(o));
    file->Close();
    if (cret <= 0 || chain->GetListOfFiles()->GetEntries() <= 0) {
      Printf("Error: FastAnalysis::Run: Failed to create chain");
      return false;
    }
    
    TString       proto    = u.GetProtocol();
    Bool_t        isProof  = (proto.EqualTo("proof") || proto.EqualTo("lite"));
    if (isProof) {
      if (!SetupProof(u,opt)) return false;
      chain->SetProof();
    }

    Printf("===================================================\n"
	   "\n"
	   " Processing chain %s with selector %p\n"
	   " Type: %s, maxEvents: %lld\n"
	   " URL: %s\n"
	   "\n"
	   "===================================================",
	   chain->GetName(), analysis, type.Data(), nev, u.GetUrl());
    if (nev < 0) nev = TChain::kBigNumber;
    Long64_t ret = chain->Process(analysis, "", nev, 0);
    
    TFile* out = TFile::Open(output, "RECREATE");
    analysis->GetOutputList()->Write("out",TObject::kSingleKey);
    out->Write();
    Printf("Saved in %s", out->GetName());

    return ret > 0;
  }    
  
  ClassDef(FastAnalysis,1);
};

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
FastAnalysis::Make(const char* type,
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
      Printf("Warning: FastAnalysis::Make: Unknown estimator: %s", w.Data());
      return 0;
    }
    if (t.BeginsWith("MULT"))
      return new MultAnalysis(w, verbose, monitor);
    else
      return new CentAnalysis(w, verbose, monitor);
  }
  Printf("Error: FastAnalysis::Run: Invalid spec: %s", t.Data());
  return 0;
}



//
//  EOF
//  
	  
