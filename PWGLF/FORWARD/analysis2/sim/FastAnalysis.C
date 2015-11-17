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




//
//  EOF
//  
	  
