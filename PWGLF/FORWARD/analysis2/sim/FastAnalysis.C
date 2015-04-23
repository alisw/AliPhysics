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
#endif

//====================================================================
/** 
 * Header in the EG tree 
 * 
 */
struct Header {
  /** The run number */
  UInt_t   fRunNo;
  /** Event number */
  UInt_t   fEventId;
  /** Number of target participants */
  UInt_t   fNtgt;
  /** Number of projectile participants */
  UInt_t   fNproj;
  /** Number of binary collisions */
  UInt_t   fNbin;
  /** Event type */
  UInt_t   fType;
  /** Collision point X coordinate */
  Double_t fIpX;
  /** Collision point Y coordinate */
  Double_t fIpY;
  /** Collision point Z coordinate */
  Double_t fIpZ;
  /** Impact parameter */
  Double_t fB;
  /** Centratlity */
  Double_t fC;
  /** Reaction plane angle */
  Double_t fPhiR;
  /** 
   * Clear the header 
   */
  void Clear(Option_t* ) {
    fRunNo   = 0;
    fEventId = 0;
    fNtgt    = 0;
    fNproj   = 0;
    fNbin    = 0;
    fType    = 0;
    fIpX     = 0;
    fIpY     = 0;
    fIpZ     = 9999;
    fB       = -1;
    fC       = -1;
    fPhiR    = -1;
  }
  /** 
   * Print the header 
   */
  void Print() {
    printf("Run No:          %9u\n"
	   "Event Id:        %9u\n"
	   "N_target:        %9u\n"
	   "N_projectile:    %9u\n"
	   "Type:            0x%07x\n"
	   "IP:          (%3.1f,%3.1f,%3.1f)\n"
	   "b:               %6f fm\n"
	   "c:               %7f %%\n"
	   "phi_R:           %9f\n",
	   fRunNo, fEventId, fNtgt, fNproj, fType,
	   fIpX, fIpY, fIpZ, fB, fC, fPhiR);
  }
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
  Header*       fHeader;     //!
  /** List of particles */
  TClonesArray* fParticles;  //!
  /** Whether to be verbose */
  Bool_t        fVerbose;    //
  /** A simple timer */
  TStopwatch    fTimer;      //!
  /** Name */
  TString       fName;       //!
  /** Cache of event discriminator */
  Double_t      fEventMult;
  /** The centrality method to use */
  TString       fCentMethod;
  /** The centrality histogram to use - if any */
  TH1*          fCentHist; //!
  /** Number of good events */
  ULong_t fOK;
  
  /** 
   * Constructor.  Opens the file passed and sets internal pointers to
   * tree, header, and particle list.
   * 
   * @param verb     Whether to be verbose 
   */  
  FastAnalysis(Bool_t verb=false)
    : fTree(0),
      fHeader(0),
      fParticles(0),
      fVerbose(verb),
      fEventMult(0),
      fCentMethod(""),
      fCentHist(0),
      fOK(0)
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
  void SetVerbose(Bool_t verb) { fVerbose = verb; }
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
    fHeader = new Header;
    fParticles = new TClonesArray("TParticle");
    fTree->SetBranchAddress("header",    &(fHeader->fRunNo));
    fTree->SetBranchAddress("particles", &fParticles);
    if (!fCentMethod.IsNull())
      fTree->SetBranchAddress(fCentMethod, &fEventMult);
    return true;
   
  }
  /** 
   * Set-up the estimator is one is chosen.  We get a pointer to the
   * current file, and retrieve the relevant histogram from that.
   * 
   * @return true on success, false otherwise 
   */
  Bool_t SetupEstimator()
  {
    if (fCentMethod.IsNull()) return true;

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
    TObject* o = f->Get(fCentMethod);
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
    Int_t bin = fCentHist->GetXaxis()->FindBin(fEventMult);
    if (bin <= 0)                      return -1;
    if (bin >= fCentHist->GetNbinsX()) return 1000;
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
    return p->GetVal();
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
      if (!oPart->TestBit((1<<14))) {
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
   * 
   * @return newly created object, or null
   */
  static FastAnalysis* Make(const char* type);
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
    ret = gROOT->ProcessLine("gProof->ClearCache()");

    TString phy = gSystem->ExpandPathName("$(ALICE_PHYSICS)");
    TString ali = gSystem->ExpandPathName("$(ALICE_ROOT)");
    TString fwd = phy + "/PWGLF/FORWARD/analysis2";

    TString load(Form("gProof->Load(\"%s/sim/FastAnalysis.C+%s\",true)",
		      fwd.Data(), opt));
    ret = gROOT->ProcessLine(load.Data());
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
    TUrl         u(url);
    TString      uout    = "";
    TObjArray*   opts    = TString(u.GetOptions()).Tokenize("&");
    TObjString*  token   = 0;
    TIter        nextToken(opts);
    while ((token = static_cast<TObjString*>(nextToken()))) {
      TString& str = token->String();
      if (str.IsNull()) continue;
      
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

    FastAnalysis* analysis = Make(type);
    if (!analysis) return false;
    
    TChain*       chain    = new TChain("T", "");
    if (!chain->AddFile(u.GetFile())) {
      Printf("Error: FastAnalysis::Run: Failed to add file %s",
	     u.GetFile());
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
  
  dNdetaAnalysis(Bool_t verbose=false)
    : FastAnalysis(verbose), fdNdeta(0)
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
    fdNdeta->Scale(1. / fOK, "width");
    SetStatus(0);
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
  NSDAnalysis(Bool_t verbose=false)
    : dNdetaAnalysis(verbose)
  {}
  /** 
   * Process the header.  Return 
   * 
   * @return True in case the event is flagged as NSD, false
   * otherwise.
   */
  virtual Bool_t ProcessHeader()
  {
    if (fHeader->fType == 0x1) return false;
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
  INELAnalysis(Bool_t verbose=false)
    : dNdetaAnalysis(verbose)
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
struct CentAnalysis : public dNdetaAnalysis
{
  TAxis*       fCentAxis;
  /// THStack*     fCentStack; //!
  TList*       fCentList;
  Int_t        fCentBin;
  const Long_t fMinEvents;
  TH1*         fCentAll; //!
  TH1*         fCentAcc; //!
  /** 
   * Constructor 
   * 
   * @param method   Centrality method 
   * @param verbose  Verbosity flag 
   */
  CentAnalysis(const char* method="V0M", Bool_t verbose=true)
    : dNdetaAnalysis(verbose),
      fCentAxis(0),
      fCentList(0),
      fCentBin(0),
      fMinEvents(100),
      fCentAll(0),
      fCentAcc(0)
  {
    fCentMethod = method;
    if (fCentMethod.Contains("RefMult")) SetCentAxis("mult");
    else                                 SetCentAxis("default");
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
    TString    s(spec);
    s.ToLower();
    if (s.IsNull()) return;
    if (s.EqualTo("pbpb") || s.EqualTo("aa") || s.EqualTo("default")) {
      Double_t aa[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
      SetCentAxis(11, aa);
      return;
    }
    if (s.EqualTo("ppb") || s.EqualTo("pbp") ||
	s.EqualTo("pa") || s.EqualTo("ap")) {
      Double_t pa[] = { 0, 5, 10, 20, 40, 60, 80, 100 };
      SetCentAxis(7, pa);
      return;
    }
    if (s.EqualTo("pp")) {
      Double_t pp[] = { 0, 0.01, 0.1, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100 };
      SetCentAxis(12, pp);
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
    fCentAll->SetDirectory(0);
    fOutput->Add(fCentAll);

    fCentAcc = static_cast<TH1*>(fCentAll->Clone("centAcc"));
    fCentAcc->SetTitle("Accepted centralities");
    fCentAcc->SetFillColor(kGreen+2);
    fCentAcc->SetDirectory(0);
    fOutput->Add(fCentAcc);

    fOutput->ls();
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
    Double_t cent = GetCentrality();
    fCentAll->Fill(cent);
    if (cent < 0 || cent > 999) return false;
    fCentBin = fCentAxis->FindBin(cent);
    if (fCentBin < 1 || fCentBin > fCentAxis->GetNbins()) {
      fCentBin = -1;
      return false;
    }
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
    TH1*     out  = static_cast<TH1*>(fCentList->/*GetHists()->*/At(fCentBin-1));
    // If we still have no histogram, return immediately 
    if (!out) return;

    // Use parent function to fill cache 
    dNdetaAnalysis::ProcessParticles();

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
    
    fCentAll   = static_cast<TH1*>(GetOutputObject("cent", TH1::Class()));
    fCentAcc   = static_cast<TH1*>(GetOutputObject("centAcc", TH1::Class()));
    // fCentStack = static_cast<THStack*>(GetOutputObject("stack",
    // THStack::Class()));
    fCentList = static_cast<TList*>(GetOutputObject("byCent",
						    TList::Class()));
    if (!fCentList || !fCentAll || !fCentAcc) {
      Warning("Terminate", "Missing stack and histograms");
      SetStatus(-1);
      return;
    }
    fCentAll->Scale(1./fOK, "width");
    fCentAcc->Scale(1./fOK, "width");

    THStack*  stack = new THStack("all", "All");
    TList*    hists = fCentList; // ->GetHists();
    TObjLink* link  = hists->FirstLink();
    while (link) {
      TObject* o = link->GetObject();
      if (!o) {
	link = link->Next();
	continue;
      }
      TH1*  h = static_cast<TH1*>(o);
      Int_t n = h->GetBinContent(0);
      printf("%9d events in bin %s ...", n, o->GetTitle());
      if (n < fMinEvents) {
	// Too few event, remove this
	TObjLink* tmp = link->Next();
	hists->Remove(link);
	link = tmp;
	delete o;
	Printf(" removed");
	continue;
      }
      // Scale
      h->Scale(1. / n, "width");
      stack->Add(h);
      Printf(" scaled");
      link = link->Next();
    }
    fOutput->Add(stack);
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
  MultAnalysis(const char* method="RefMult00d80", Bool_t verbose=false)
    : CentAnalysis(method, verbose)
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
FastAnalysis::Make(const char* type)
{
  TString t(type);
  if      (t.EqualTo("INEL"))    return new INELAnalysis;
  else if (t.EqualTo("NSD"))     return new NSDAnalysis;
  else if (t.BeginsWith("MULT") || t.BeginsWith("CENT")) {
    TString w(t(4, t.Length()-4));
    if (!(w.EqualTo("RefMult00d80") ||
	  w.EqualTo("RefMult00d50") ||
	  w.EqualTo("V0M") ||
	  w.EqualTo("V0A") ||
	  w.EqualTo("V0C"))) {
      Printf("Warning: FastAnalysis::Make: Unknown estimator: %s", w.Data());
      return 0;
    }
      
    if (t.BeginsWith("MULT"))
      return new MultAnalysis(w);
    else
      return new CentAnalysis(w, true);
  }
  Printf("Error: FastAnalysis::Run: Invalid spec: %s", t.Data());
  return 0;
}



//
//  EOF
//  
	  
