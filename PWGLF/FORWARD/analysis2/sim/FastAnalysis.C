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
# include "FastShortHeader.C"
# include "FastMonitor.C"
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
class FastShortHeader;
#endif

//====================================================================
/** 
 * Base class for processors 
 */
struct FastAnalysis : public TSelector
{
  /** Pointer to tree being analysed */
  TTree*        fTree; //! 
  /** Cache of our header */
  FastShortHeader*  fHeader;     //!
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
  Bool_t fCompatB;
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
      fMonitor(monitor),
      fCompatB(false)
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

    FastMonitor* monitor = new FastMonitor(this, "FastMonitor");
    TObject* obj = 0;
    TIter    next(objs);
    while ((obj = next())) monitor->Register(obj);
    monitor->Connect(fMonitor);
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
    fHeader = new FastShortHeader;
    fParticles = new TClonesArray("TParticle");
    fTree->SetBranchAddress("header",    &(fHeader->fRunNo));
    fTree->SetBranchAddress("particles", &fParticles);
    if (!fCentMethod.IsNull()) {
      if (fCentMethod.EqualTo("B") && !fTree->GetBranch("B"))
	fCompatB = true;
      else 
	fTree->SetBranchAddress(fCentMethod, &fEventMult);
    }
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
    if (fCompatB) {
      Warning("SetupEstimator", "Using fixed estimator");
      Double_t bs[] = { 0,      1.57,  2.22,  2.71,  3.13,
			3.50,   4.94,  6.05,  6.98,  7.81,
			8.55,   9.23,  9.88, 10.47, 11.04,
			11.58, 12.09, 12.58, 13.05, 13.52,
			13.97, 14.43, 14.96, 15.67, 20.00 };
      Double_t cs[] = { 0.5,   1.5,   2.5,   3.5,   4.5,
			7.5,   12.5,  17.5,  22.5,  27.5,
			32.5,  37.5,  42.5,  47.5,  52.5,
			57.5,  62.5,  67.5,  72.5,  77.5,
			82.5,  87.5,  92.5,  97.5 };
      TH1* h = new TH1D("rawB", "B to Centrality", 28, bs);
      h->SetDirectory(0);
      h->SetXTitle("b\\hbox{ [fm]}");
      h->SetYTitle("c\\hbox{ [\\%]}");
      h->SetBinContent(0,1);
      for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
	h->SetBinContent(i, cs[i-1]);
      }
      // fCentHist = h;
      return true;
    }
    
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
  /** 
   * At beginning of job
   * 
   */
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
    CopyEgHistogram();
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
#if 0
    Info("GetCentrality", "fCentHist=%p, fCentMethod=%s cent=%f",
	 fCentHist, fCentMethod.Data(), fHeader->fC);
#endif 
    if (!fCentHist) return fHeader->fC;
    Int_t nBin = fCentHist->GetNbinsX();
    Int_t bin  = fCentHist->GetXaxis()->FindBin(fEventMult);
    if (bin <= 0)    {
      Warning("", "Look-up of %f failed (%d)", fEventMult, bin);
      return -1;
    }
    if (bin-1 == nBin && fEventMult == fCentHist->GetXaxis()->GetXmax()) {
      bin        =  nBin;
    }
    if (bin >  nBin) {
      Warning("", "Look-up of %f failed (%d > %d)", fEventMult, bin, nBin);
      return 200;
    }
    Double_t cent = fCentHist->GetBinContent(bin);
    // Info("", "Look-up of %ld -> %d -> %5.1f%%", fEventMult, bin, cent);
    return cent;
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
  void CopyEgHistogram()
  {
    if (fOutput->FindObject("eg")) return;
    TFile* file = (fTree ? fTree->GetCurrentFile() : 0);
    Info("CopyEgHistogram", "Current file: (%p) %s", file,
	 (file ? file->GetName() : ""));
    if (!file) return;
    TH1* egTitle = static_cast<TH1*>(file->Get("eg"));
    if (!egTitle) {
      Warning("CopyEgHistogram", "No EG histogram");
      return;
    }
    TH1* copy = static_cast<TH1*>(egTitle->Clone());
    copy->SetDirectory(0);
    copy->Reset();
    copy->Fill(0.5);
    fOutput->Add(copy);
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
  static Bool_t ProofLoad(const char* file, const char* opt)
  {
    const char* real = gSystem->Which(gROOT->GetMacroPath(),
				      file, kReadPermission);
    if (!real) {
      ::Warning("ProofLoad", "%s not found in load path: %s",
		file, gROOT->GetMacroPath());
      delete real;
      return false;
    }
    TString cc(opt); if (cc == "-") cc = "" ; else cc.Prepend("+");
    Long_t ret = ProofExec(Form("Load(\"%s%s\",true)",real, cc.Data()));
    if (ret != 0) {
      ::Warning("ProofLoad", "Failed to load %s%s", real, cc.Data());
      delete real;
      return false;
    }
    delete real;
    return true;
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
  static Bool_t SetupProof(TUrl& url, const char* opt,
			   const char* extra)
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
    gROOT->SetMacroPath(Form("%s:%s/sim", 
			     gROOT->GetMacroPath(), fwd.Data()));
    ProofExec(Form("Exec(\"gSystem->SetMakeSharedLib(\\\"%s\\\")\")",
		   mkLib.Data()));

    if (!ProofLoad("FastMonitor.C", opt)) return false;
    if (!ProofLoad("FastShortHeader.C", opt)) return false;
    if (!ProofLoad("FastAnalysis.C", opt)) return false;
    if (!ProofLoad("FastCentHelper.C", opt)) return false;
    if (!extra || extra[0] == '\0') return true;
    if (!ProofLoad(extra, opt)) return false;

    return true;
  }

  //==================================================================
  /**
   * A class that implements construction of specific analysers.
   * Sub-classes should implement the Make function to either return a
   * fast analyser of the requested type, or null.
   * 
   */
  struct Maker : public TNamed
  {
    /** DTOR */
    virtual ~Maker() {}
    /** 
     * Pure virtual function to implement 
     * 
     * @param subtype Type of analyser 
     * @param monitor Monitor period 
     * @param verbose Verbosity 
     * @param uout    Possible options 
     * 
     * @return Analyer or null
     */
    virtual FastAnalysis* Make(const TString& subtype,
			       Int_t          monitor,
			       Bool_t         verbose,
			       TString&       uout) = 0;
    /** 
     * List available sub-types
     */
    virtual void List() const = 0;
    /** 
     * Script to load 
     */
    virtual const char* Script() const = 0;
  protected:
    /** CTOR */
    Maker(const char* type="");    
    ClassDef(Maker,0);
  };
  //==================================================================
#ifndef __CINT__ 
  /** 
   * A class that can create analysers from different makers. 
   */
  struct Factory : public TObject
  {
    /** 
     * Singleton function 
     */
    static Factory& Instance()
    {
      static Factory* instance = 0;
      if (!instance) instance = new Factory;
      return *instance;
    }
    /** 
     * Make an analysis. 
     * 
     * @param type     Type of analyser
     * @param subtype  Sub-type of analyser
     * @param monitor  Monitor period in seconds 
     * @param verbose  Whether to be verbose 
     *
     * @param uout     Possible options.  A maker can modify this to
     *                 extract optoins for the analyser .
     * 
     * @return analyser or null
     */
    FastAnalysis* Make(const TString& type,
		       const TString& subtype,
		       Int_t          monitor,
		       Bool_t         verbose,
		       TString&       uout)
    {
      Bool_t help = (type.EqualTo("help",TString::kIgnoreCase) ||
		     type.EqualTo("list",TString::kIgnoreCase));
      TIter next(&fList);
      Maker* maker = 0;
      while ((maker = static_cast<Maker*>(next()))) {
	if (!help && !type.EqualTo(maker->GetName(), TString::kIgnoreCase))
	  continue;
	if (help) Printf("%s", maker->GetName());
	if (help ||
	    subtype.EqualTo("help", TString::kIgnoreCase) ||
	    subtype.EqualTo("list", TString::kIgnoreCase)) {
	  maker->List();
	  continue;
	}
	    
	FastAnalysis* a = maker->Make(subtype,monitor,verbose,uout);
	if (a) return a;
      }
      return 0;      
    }
    const char* Script(const TString& type) const
    {
      TIter next(&fList);
      Maker* maker = 0;
      while ((maker = static_cast<Maker*>(next()))) {
	if (!type.EqualTo(maker->GetName(), TString::kIgnoreCase))
	  continue;
	return maker->Script();
      }
      return 0;
    }
    /** 
     * Register maker 
     * 
     * @param m Maker 
     */
    void Register(Maker* m)
    {
      fList.Add(m);
    }
  private: 
    Factory() : fList() {}
    TList fList;
  };
#endif
  /** 
   * Run an analysis.  Argument URL has format of 
   *
   *  protocol://[[user@[password:]]host]/file?[options]#treeName
   * 
   * where 
   *
   * - protocol can be @c local, @c lite, or @c proof 
   * - user@password:host gives Proof credentials and master 
   * - file is either a data file, or contains a list of file names 
   * - options are options for the execution environment 
   * - treeName is the input tree name 
   * 
   * @param url     Processing and input URL
   * @param output  Output file name 
   * @param a       Analyser 
   * @param nev     Max number of events
   * @param offset  Offset in events
   * @param monitor Monitor period in seconds (<0 disables)
   * @param verbose Whether to be verbose 
   * @param opt     Optimization used 
   * 
   * @return        true on success
   */
  static Bool_t Run(const char*   url,
		    const char*   output,
		    FastAnalysis* a,
		    const char*   script,
		    Long64_t      nev=-1,
		    Long64_t      offset=0,
		    Int_t         monitor=-1,
		    Bool_t        verbose=false,
		    const char*   opt="")
    
  {
    if (!a) {
      Printf("Error: FastAnalysis::Run: No analyser given");
      return false;
    }
    TUrl u(url);
    if (!u.IsValid()) {
      Printf("Error: FastAnalysis::Run: URL %s is invalid", u.GetUrl());
      return false;
    }

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
      if (!SetupProof(u,opt,script)) return false;
      chain->SetProof();
    }

    Printf("===================================================\n"
	   "\n"
	   " Processing chain %s with selector %s\n"
	   " Max events: %lld\n"
	   " Event offset: %lld\n"
	   " URL: %s\n"
	   "\n"
	   "===================================================",
	   chain->GetName(), a->GetName(), nev, offset, u.GetUrl());
    if (nev < 0) nev = TChain::kBigNumber;
    Long64_t ret = chain->Process(a, "", nev, offset);

    if (output && output[0] != '\0') {
      TFile* out = TFile::Open(output, "RECREATE");
      a->GetOutputList()->Write("out",TObject::kSingleKey);
      out->Write();
      Printf("Saved in %s", out->GetName());
    }

    return ret > 0;
  }
  
  /** 
   * Run this.  
   * 
   * @param url     Url to process 
   * @param output  Output file 
   * @param opt     Compilation options 
   * 
   * @return true on success
   */
  static Bool_t Run(const char* url,
		    const char* output,
		    const char* opt="g")
  {
    Long64_t     nev     = -1;
    Long64_t     off     = 0;
    TString      type    = "";
    TString      sub     = "";
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
      
      if      (key.EqualTo("events"))  nev     = val.Atoll();
      else if (key.EqualTo("offset"))  off     = val.Atoll();
      else if (key.EqualTo("type"))    type    = val;
      else if (key.EqualTo("subtype")) sub     = val;
      else if (key.EqualTo("monitor")) monitor = val.Atoi();
      else {
	if (!uout.IsNull()) uout.Append("&");
	uout.Append(str);
      }
    }
    opts->Delete();
    FastAnalysis* a = Factory::Instance().Make(type,sub, monitor,verbose,uout);
    const char*   s = Factory::Instance().Script(type);
    if (type.EqualTo("help",TString::kIgnoreCase) ||
	type.EqualTo("list",TString::kIgnoreCase) ||
	sub .EqualTo("help",TString::kIgnoreCase) ||
	sub .EqualTo("list",TString::kIgnoreCase)) return false;
    u.SetOptions(uout);

    return Run(u.GetUrl(), output, a, s, nev, off, monitor, verbose, opt);
  }    
  
  ClassDef(FastAnalysis,1);
};

FastAnalysis::Maker::Maker(const char* type)
  : TNamed(type, "")
{
  // Automatically register 
  FastAnalysis::Factory::Instance().Register(this);
}


//
//  EOF
//  
	  
