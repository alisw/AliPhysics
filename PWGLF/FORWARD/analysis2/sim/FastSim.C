#ifndef FASTSIM_H
#define FASTSIM_H
#include <TSelector.h>
#include <TQObject.h>
#ifndef __CINT__
# include "AliGenerator.h"
# include "AliRunLoader.h"
# include "AliStack.h"
# include "AliHeader.h"
# include "AliGenEventHeader.h"
# include "AliRun.h"
# include "AliCollisionGeometry.h"
# include "AliGenPythiaEventHeader.h"
# include "AliGenDPMjetEventHeader.h"
# include "AliGenGeVSimEventHeader.h"
# include "AliGenHerwigEventHeader.h"
# include <TROOT.h>
# include <TString.h>
# include <TMath.h>
# include <TParticle.h>
# include <TH1.h>
# include <TTree.h>
# include <TClonesArray.h>
# include <TList.h>
# include <TProof.h>
# include <TParticlePDG.h>
# include <TStopwatch.h>
# include <TFile.h>
# include <TProofOutputFile.h>
# include <TCanvas.h>
# include <TTimer.h>
# include <fstream>
#else
class AliGenerator;
class AliRunLoader;
class AliStack;
class AliHeader;
class AliGenEventHeader;
class TH1;
class TTree;
class TClonesArray;
class TBrowser;
class TList;
class TFile;
class TProofOutputFile;
class TCanvas;
class TVirtualPad;
class TTimer;
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
    : fName("FastMonitor"),
      fCanvas(0),
      fSelector(s)
  {
    if (gROOT->IsBatch()) {
      Warning("FastMonitor", "Batch processing, no monitoring");
      return;
    }

    if (gProof) {
      fName = gProof->GetSessionTag();
      gDirectory->Add(this);
      Bool_t ret = gProof->Connect("Feedback(TList *objs)", "FastMonitor", this, 
				   "Feedback(TList *objs)");
      if (!ret) {
	Warning("FastMonitor", "Failed to connect to Proof");
	return;
      }
    }
    else if (!s) return;
    
    fCanvas = new TCanvas(fName, Form("Monitor %s", fName.Data()), 1000, 800);
    fCanvas->SetFillColor(0);
    fCanvas->SetFillStyle(0);
    fCanvas->SetTopMargin(0.01);
    fCanvas->SetRightMargin(0.01);

    fCanvas->Divide(2,2);
    RegisterDraw(1, "type",   "", 0);
    RegisterDraw(2, "b",      "", 0);
    RegisterDraw(3, "cent",   "", 0);
    RegisterDraw(4, "dNdeta", "", 0x8);
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
   * Called when we get notified of 
   * 
   * @param objs List of monitored objects
   */
  void Feedback(TList* objs)
  {
    // Info("FeedBack", "List is %p", objs);
    // if (objs) objs->ls();
    if (!fCanvas) return;

    TList* l = static_cast<TList*>(objs->FindObject("histograms"));
    if (!l) {
      Warning("Feedback", "No histograms");
      return;
    }
    Int_t nEvents = 1;
    TObject* oIpz = l->FindObject("ipZ");
    if (oIpz && oIpz->IsA()->InheritsFrom(TH1::Class())) 
      nEvents = static_cast<TH1*>(oIpz)->GetEntries();
    else 
      Warning("Feedback", "Histogram ipZ not found");
    
    TIter next(l);
    TObject* o = 0;
    while ((o = next())) {
      TVirtualPad* p = FindPad(o->GetName());
      if (!p) 
	// Info("FeedBack", "no pad for %s", o->GetName());
	continue;

      p->cd();
      if (o->IsA()->InheritsFrom(TH1::Class())) {
	TH1* h = static_cast<TH1*>(o);
	TH1* c = h->DrawCopy(p->GetTitle());
	c->SetDirectory(0);
	c->SetBit(TObject::kCanDelete);
	if (p->TestBit(BIT(15))) {
	  Info("Feedback", "Scaling %s by 1./%d and width",
	       c->GetName(), nEvents);
	  c->Scale(1./nEvents, "width");
	}
      }
      else {
	TObject* c = o->DrawClone(p->GetTitle());
	c->SetBit(TObject::kCanDelete);
      }
      p->Modified();
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
    Info("HandleTimer", "Selector=%p", fSelector);
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


//====================================================================
/** 
 * Run a event generator simulation 
 */
struct FastSim : public TSelector
{
  /** 
   * Constructor 
   * 
   * @param eg     Event generator 
   * @param runNo  Run number to simulate 
   * @param bMin   Lease impact parameter 
   * @param bMax   Largest impact parameter 
   */
  FastSim(const char* eg="",
	  ULong_t runNo=0,
	  Double_t bMin=0,
	  Double_t bMax=20,
	  Long64_t nEvents=0)
    : TSelector(),
      fEGName(eg),
      fRunNo(runNo),
      fBMin(bMin),
      fBMax(bMax),
      fGRP(0),
      fNEvents(nEvents),
      fGenerator(0),
      fRunLoader(0),
      fStack(0),
      fHeader(0),
      fTree(0),
      fParticles(0),
      fList(0),
      fHEta(0),
      fHIpz(0),
      fHType(0),
      fHCent(0),
      fHB(0),
      fHPhiR(0),
      fHTime(0),
      fProofFile(0),
      fFile(0),
      fFileName("")
  {}
  const char* FileName() const
  {
    static TString fn;
    if (fn.IsNull()) {
      if (!fFileName.IsNull())  fn = fFileName;
      else {
	const char* egName = (fGenerator ?
			      fGenerator->GetName() :
			      fEGName.Data());
	fn = Form("%s_%09d", egName, fRunNo);
	if (fNEvents > 0) {
	  if (fNEvents >= 1000000)
	    fn.Append(Form("_%lldM", fNEvents/1000000));
	  else if (fNEvents >= 1000)
	    fn.Append(Form("_%lldk", fNEvents/1000));
	  else
	    fn.Append(Form("_%lld", fNEvents));
	}
	fn.Append(".root");
	fFileName = fn;
      }
    }
    return fn.Data();
    /*
      if (fFileName.IsNull())
      fFileName = Form("%s_%09d.root", fEGName.Data(), fRunNo);
      return fFileName.Data();*/
  }
  const char* GetName() const { return "FastSim"; }
  const char* GetTitle() const { return "ALICE Event Generator simulation"; }
  /** 
   * Create our outputs 
   * 
   * 
   * @return true on success 
   */
  Bool_t SetupOutput()
  {
    Info("SetupOutput", "First the file");
    Bool_t isProof = false;
    if (fInput && fInput->FindObject("PROOF_Ordinal"))
      isProof = true;
    if (isProof) {
      Info("SetupOutput", "Making Proof File");
      fProofFile = new TProofOutputFile(FileName(), "M");
      // TProofOutputFile::kMerge,
      // TProofOutputFile::kRemote);
      fFile = fProofFile->OpenFile("RECREATE");
    }
    else
      fFile = TFile::Open(FileName(), "RECREATE");

    Info("SetupOutput", "Making our tree");
    fTree      = new TTree("T", "T");
    fParticles = new TClonesArray("TParticle");
    fTree->Branch("header", &fShortHead,
		  "run/i:event:npart:nbin:type:ipx/D:ipy:ipz:b:c:phir");
    fTree->Branch("particles", &fParticles);
    fTree->AutoSave();
    fTree->SetDirectory(fFile);
    fTree->SetAlias("primary", "(particles.fBits&(1<<14))");
    fTree->SetAlias("weak",    "(particles.fBits&(1<<15))");
    fTree->SetAlias("charged", "(particles.fBits&(1<<16))");
    fTree->SetAlias("pt",      "(sqrt(pow(particles.fPx,2)+"
		    /*       */"pow(particles.fPy,2)))");
    fTree->SetAlias("eta",     "(pt<1e-10?1024:"
		    "-log(tan(atan2(particles.Pt(),particles.fPz)/2)))");
    fTree->SetAlias("good",    "(primary&&charged&&abs(eta)<1000)");
    fTree->SetAlias("sd",      "(header.fType & 0x1)");
    fTree->SetAlias("dd",      "(header.fType & 0x2)");
    fTree->SetAlias("pion",    "(abs(particles.fPdgCode)==211)");
    fTree->SetAlias("kaon",    "(abs(particles.fPdgCode)==321)");
    fTree->SetAlias("proton",  "(abs(particles.fPdgCode)==2212)");
    fTree->SetAlias("electron","(abs(particles.fPdgCode)==11)");
    fTree->SetAlias("other",   "(!pion&&!kaon&&!proton&&!electron)");
    fTree->SetAlias("beta",    "(particles.P()/particle.Energy())");
    fTree->SetAlias("gamma",   "(1./sqrt(1-beta*beta))");

    Info("SetupOutput", "Making histograms");
    Double_t maxEta = 10;
    Double_t dEta   = 10./200;
    fHEta = new TH1D("dNdeta", "Charged particle pseudo-rapidity density",
		     Int_t(2*maxEta/dEta+.5), -maxEta, +maxEta);
    fHEta->Sumw2();
    fHEta->SetXTitle("#it{#eta}");
    fHEta->SetYTitle("1/N d#it{N}_{ch}/d#it{#eta}");
    fHEta->SetMarkerColor(kRed+2);
    fHEta->SetMarkerStyle(20);
    fHEta->SetDirectory(0);
    
    fHIpz = new TH1D("ipZ", "Z-coordinate of interaction point",
		     10, -10, 10);
    fHIpz->SetMarkerColor(kGreen+2);
    fHIpz->SetFillColor(kGreen+2);
    fHIpz->SetFillStyle(3001);
    fHIpz->SetXTitle("IP_{#it{z}} [cm]");
    fHIpz->SetDirectory(0);

    fHType = new TH1D("type", "Diffractive", 3, .5, 3.5);
    fHType->SetMarkerColor(kOrange+2);
    fHType->SetFillColor(kOrange+2);
    fHType->SetFillStyle(3001);
    fHType->SetDirectory(0);
    fHType->GetXaxis()->SetBinLabel(1, "Non");
    fHType->GetXaxis()->SetBinLabel(2, "Single");
    fHType->GetXaxis()->SetBinLabel(3, "Double");

    fHCent = new TH1D("cent", "Centrality", 20, 0, 100);
    fHCent->SetMarkerColor(kPink+2);
    fHCent->SetFillColor(kPink+2);
    fHCent->SetFillStyle(3001);
    fHCent->SetDirectory(0);
    fHCent->SetXTitle("Centrality [%]");
    
    fHB = new TH1D("b", "Impact parameter", 20, 0, 20);
    fHB->SetMarkerColor(kYellow+2);
    fHB->SetFillColor(kYellow+2);
    fHB->SetFillStyle(3001);
    fHB->SetXTitle("#it{b} [fm]");
    fHB->SetDirectory(0);

    fHPhiR = new TH1D("phiR", "Event plane angle", 360, 0, 360);
    fHPhiR->SetMarkerColor(kMagenta+2);
    fHPhiR->SetFillColor(kMagenta+2);
    fHPhiR->SetFillStyle(3001);
    fHPhiR->SetXTitle("#it{#Phi}_{R} [degrees]");
    fHPhiR->SetDirectory(0);

    fHTime = new TH1D("timing", "Timing of processing", 5,0.5,5.5);
    fHTime->SetMarkerColor(kBlue+2);
    fHTime->SetFillColor(kBlue+2);
    fHTime->SetFillStyle(3001);
    fHTime->GetXaxis()->SetBinLabel(1,"Reset");
    fHTime->GetXaxis()->SetBinLabel(2,"Generate");
    fHTime->GetXaxis()->SetBinLabel(3,"Header");
    fHTime->GetXaxis()->SetBinLabel(4,"Particles");
    fHTime->GetXaxis()->SetBinLabel(5,"Filling");
    fHTime->SetDirectory(0);
				    
    fList = new TList;
    fList->SetName("histograms");
    fList->SetOwner(true);
    fList->Add(fHEta);
    fList->Add(fHIpz);
    fList->Add(fHType);
    fList->Add(fHCent);
    fList->Add(fHB);
    fList->Add(fHPhiR);
    fList->Add(fHTime);

    Info("SetupOutput", "Adding list ot outputs");
    fOutput->Add(fList);

    return true;
  }
  Bool_t SetupGen()
  {
    Printf(" === Setup ==============================");
    Printf("  Run #:                          %6d", fRunNo);
    Printf("  EG:     %30s", fEGName.Data());
    Printf("  B range:             %5.1ffm - %5.1ffm", fBMin, fBMax);
    Printf(" ========================================");
    Printf("Macro path: %s", gROOT->GetMacroPath());

    // --- Check if we shoud get the GRP line ------------------------
    if (!fGRP && fInput) {
      fGRP = fInput->FindObject("GRP");
      std::ofstream* pout = new std::ofstream("grp.dat");
      if (pout) {
	Info("SetupGen", "Writing GRP line '%s' to \"grp.dat\"",
	     fGRP->GetTitle());
	std::ostream& out = *pout;
	out << fGRP->GetTitle() << std::endl;
	pout->close();
      }
    }

    // --- Load our settings -----------------------------------------
    Info("SetupGen", "Loading scripts");
    gROOT->Macro(Form("GRP.C(%d)", fRunNo));
    gROOT->Macro("BaseConfig.C");
    gROOT->Macro("EGConfig.C");

    gROOT->ProcessLine(Form("VirtualEGCfg::LoadGen(\"%s\")",fEGName.Data()));

    // --- Make our generator ----------------------------------------
    Info("SetupGen", "Creating generator");
    TString egMk = Form("egCfg->MakeGenerator(\"%s\",%f,%f)",
			fEGName.Data(), fBMin, fBMax);
    Long64_t egPtr = gROOT->ProcessLine(egMk);
    if (egPtr == 0) {
      Error("Setup", "Failed to make generator");
      return false;
    }
    fGenerator = reinterpret_cast<AliGenerator*>(egPtr);


    if (fFileName.IsNull()) FileName();
    Info("SetupRun", "File name is '%s'", fFileName.Data());

    return true;
  }    
  /** 
   * Setup the generator etc. of the job 
   * 
   * @param nev Maximum number of events per file 
   * 
   * @return true on success 
   */
  Bool_t SetupRun()
  {
    // --- gAlice (bare ROOT) ----------------------------------------
    if (!gAlice)
      new AliRun("gAlice", "The ALICE Off-line framework");

    Long64_t nev = (fNEvents <= 0 ? 0xFFFFFFFF : fNEvents);
    // --- Run-loader, stack, etc  -----------------------------------
    Info("SetupRun", "Set-up run Loader");    
    fRunLoader = AliRunLoader::Open("galice.root", "FASTRUN", "RECREATE");
    fRunLoader->SetCompressionLevel(2);
    fRunLoader->SetNumberOfEventsPerFile(nev);
    fRunLoader->LoadKinematics("RECREATE");
    fRunLoader->MakeTree("E");
    gAlice->SetRunLoader(fRunLoader);
    fRunLoader->MakeStack();
    fStack  = fRunLoader->Stack();
    fHeader = fRunLoader->GetHeader();

    // --- Initialize generator --------------------------------------
    Info("SetupRun", "Initializing generator");
    fGenerator->Init();
    fGenerator->SetStack(fStack);

    return true;
  }
  /** 
   * Read the previously created grp.dat file 
   * 
   */
  Bool_t ReadGRPLine()
  {
    std::ifstream* pin = new std::ifstream("grp.dat");
    if (!pin) {
      Warning("ReadGRPLine", "Failed to open \"grp.dat\"");
      return false;
    }
    std::istream&  in  = *pin;
    TString line;
    TString env;
    do {
      line.ReadLine(in);
      if (line.IsNull()) continue;
      if (line.BeginsWith("#")) continue;
      env = line;
      break;
    } while (!in.eof());
    pin->close();

    if (env.IsNull()) {
      Warning("ReadGRPLine", "Got no line from \"grp.dat\"");
      return false;
    }
    
    fGRP = new TNamed("GRP",env.Data());
    return true;
  }
    
  /** 
   * Set up job 
   * 
   */
  void Init(TTree*)
  {
  }
  /** 
   * Set up job 
   * 
   */
  void Begin(TTree*)
  {
    // Make a monitor
    Info("Begin", "gProof=%p Nomonitor=%p",
	 gProof, (gProof ? gProof->GetParameter("NOMONITOR") : 0));

    if (gProof && !gProof->GetParameter("NOMONITOR")) { 
      new FastMonitor;
      gProof->AddFeedback("histograms");
      Info("Begin", "Adding monitoring");
    }
    gROOT->Macro(Form("GRP.C(%d)", fRunNo));
    if (ReadGRPLine()) {
      if(gProof) {
	gProof->AddInput(fGRP);
      }
    }
  }
  /** 
   * Set-up this sub-job 
   * 
   */
  void SlaveBegin(TTree*)
  {
    SetupGen();
    SetupOutput();
    SetupRun();
  }
  /* Reset internal caches etc. 
   * 
   * @param iEv Event number 
   * 
   * @return true on success
   */
  Bool_t PreEvent(Long64_t iEv)
  {
    // --- Reset header ----------------------------------------------
    fShortHead.fRunNo   = fRunNo;
    fShortHead.fEventId = iEv;
    fShortHead.fIpX     = 1024;
    fShortHead.fIpY     = 1024;
    fShortHead.fIpZ     = 1024;
    fShortHead.fNpart   = -1;
    fShortHead.fNbin    = -1;
    fShortHead.fPhiR    = -1;
    fShortHead.fB       = -1;
    fShortHead.fC       = -1; 
    fParticles->Clear();
    // --- Reset header, etc.  ---------------------------------------
    fHeader->Reset(fRunNo, iEv);
    fRunLoader->SetEventNumber(iEv);
    fStack->Reset();
    fRunLoader->MakeTree("K");

    return true;
  }
  /** 
   * Process the event header 
   * 
   * @return true if the event should be diagnosed
   */
  Bool_t ProcessHeader()
  {
    // --- Copy to short header --------------------------------------
    fShortHead.fRunNo   = fHeader->GetRun();
    fShortHead.fEventId = fHeader->GetEvent();
    TArrayF ip;
    fHeader->GenEventHeader()->PrimaryVertex(ip);
    fShortHead.fIpX     = ip[0];
    fShortHead.fIpY     = ip[1];
    fShortHead.fIpZ     = ip[2];

    // --- Check header type -----------------------------------------
    AliGenEventHeader* genHeader = fHeader->GenEventHeader();
    AliCollisionGeometry* geometry = 
      dynamic_cast<AliCollisionGeometry*>(genHeader);
    AliGenPythiaEventHeader* pythia    = 
      dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
    AliGenDPMjetEventHeader* dpm       = 
      dynamic_cast<AliGenDPMjetEventHeader*>(genHeader);
    AliGenGeVSimEventHeader* gev       = 
      dynamic_cast<AliGenGeVSimEventHeader*>(genHeader);
    AliGenHerwigEventHeader* herwig    = 
      dynamic_cast<AliGenHerwigEventHeader*>(genHeader);
    if (geometry) {
      fShortHead.fB     = geometry->ImpactParameter();
      fShortHead.fNpart = (geometry->ProjectileParticipants() + 
			   geometry->TargetParticipants());
      fShortHead.fNbin  = geometry->NN();
      fShortHead.fPhiR  = geometry->ReactionPlaneAngle();
    }
    // --- Determine diffraction flags -------------------------------
    Bool_t sd = false;
    Bool_t dd = false;
    if (pythia) {
      Int_t type = pythia->ProcessType();
      if (type < 100) { // pythia6
	switch (type) {
	case 92: case 93: sd = true; break;
	case 94:          dd = true; break;
	}
      }
      else {
	switch (type) { // Pythia8
	case 103: case 104: sd = true; break;
	case 105:           dd = true; break;
	}
      }
      fShortHead.fB     = pythia->GetImpactParameter();
      fShortHead.fNpart = 2;
      fShortHead.fNbin  = 1;
    }
    if (dpm) {
      Int_t type = dpm->ProcessType();
      switch (type) {
      case 5: case 6: sd = true;

      case 7:         dd = true;
      }
    }
    if (gev) fShortHead.fPhiR = gev->GetEventPlane();
    if (herwig) {
      Int_t type = herwig->ProcessType();
      switch (type) {
      case 5: case 6: sd = true; break;
      }
      fShortHead.fNpart = 2;
      fShortHead.fNbin  = 1;
    }
    fShortHead.fType = (sd ? 0x1 : 0) | (dd ? 0x2 : 0);

    // --- Check centrality -----------------------------------------
    // PbPb only 
    // Updated 4th of November 2014 from 
    // cern.ch/twiki/bin/view/ALICE/CentStudies#Tables_with_centrality_bins_AN1
    Float_t  np = 0;
    UInt_t   nc = 0;
    Double_t c  = -1;
    Double_t b  = fShortHead.fB;
    if (b >= 0) {
      if      (0.00 >= b  && b < 1.57)  { c=0.5;  np=403.8; nc=1861; } 
      else if (1.57 >= b  && b < 2.22)  { c=1.5;  np=393.6; nc=1766; } 
      else if (2.22 >= b  && b < 2.71)  { c=2.5;  np=382.9; nc=1678; } 
      else if (2.71 >= b  && b < 3.13)  { c=3.5;  np=372;   nc=1597; }  
      else if (3.13 >= b  && b < 3.50)  { c=4.5;  np=361.1; nc=1520; } 
      else if (3.50 >= b  && b < 4.94)  { c=7.5;  np=329.4; nc=1316; } 
      else if (4.94 >= b  && b < 6.05)  { c=12.5; np=281.2; nc=1032; } 
      else if (6.05 >= b  && b < 6.98)  { c=17.5; np=239;   nc=809.8; }
      else if (6.98 >= b  && b < 7.81)  { c=22.5; np=202.1; nc=629.6; }
      else if (7.81 >= b  && b < 8.55)  { c=27.5; np=169.5; nc=483.7; }
      else if (8.55 >= b  && b < 9.23)  { c=32.5; np=141;   nc=366.7; }
      else if (9.23 >= b  && b < 9.88)  { c=37.5; np=116;   nc=273.4; }
      else if (9.88 >= b  && b < 10.47) { c=42.5; np=94.11; nc=199.4; } 
      else if (10.47 >= b && b < 11.04) { c=47.5; np=75.3;  nc=143.1; } 
      else if (11.04 >= b && b < 11.58) { c=52.5; np=59.24; nc=100.1; }
      else if (11.58 >= b && b < 12.09) { c=57.5; np=45.58; nc=68.46; }
      else if (12.09 >= b && b < 12.58) { c=62.5; np=34.33; nc=45.79; }
      else if (12.58 >= b && b < 13.05) { c=67.5; np=25.21; nc=29.92; }
      else if (13.05 >= b && b < 13.52) { c=72.5; np=17.96; nc=19.08; }
      else if (13.52 >= b && b < 13.97) { c=77.5; np=12.58; nc=12.07; }
      else if (13.97 >= b && b < 14.43) { c=82.5; np=8.812; nc=7.682; }
      else if (14.43 >= b && b < 14.96) { c=87.5; np=6.158; nc=4.904; }
      else if (14.96 >= b && b < 15.67) { c=92.5; np=4.376; nc=3.181; }
      else if (15.67 >= b && b < 20.00) { c=97.5; np=3.064; nc=1.994; }
      fShortHead.fC = c;
      // Be careful to round off
      if (fShortHead.fNpart <= 0) fShortHead.fNpart = Int_t(np+.5);
      if (fShortHead.fNbin  <= 0) fShortHead.fNbin  = Int_t(nc+.5)/2;
    }
    
    // --- Check if within vertex cut -------------------------------
    Bool_t selected = (fShortHead.fIpZ <= fHIpz->GetXaxis()->GetXmax() &&
		       fShortHead.fIpZ >= fHIpz->GetXaxis()->GetXmin());

    // --- Only update histograms if within IPz cut ------------------
    if (selected) {
      fHPhiR->Fill(fShortHead.fPhiR*TMath::RadToDeg());
      fHB->Fill(fShortHead.fB);
      fHIpz->Fill(fShortHead.fIpZ);
      fHType->Fill(dd ? 3 : sd ? 2 : 1);
      fHCent->Fill(c);
    }
    return selected;
  }
  /** 
   * Process all particles 
   * 
   * @param selected True if particle information should be diagnosed
   * 
   * @return true on success
   */
  Bool_t ProcessParticles(Bool_t selected)
  {
    Int_t nPart = fStack->GetNprimary();
    for (Int_t iPart = 0; iPart < nPart; iPart++) {
      TParticle*    particle  = fStack->Particle(iPart);
      TParticlePDG* pdg       = particle->GetPDG();
      Bool_t        primary   = fStack->IsPhysicalPrimary(iPart);
      Bool_t        weakDecay = fStack->IsSecondaryFromWeakDecay(iPart);
      Bool_t        charged   = (pdg &&  TMath::Abs(pdg->Charge()) > 0);
      if (primary)   particle->SetBit(BIT(14));
      if (weakDecay) particle->SetBit(BIT(15));
      if (charged)   particle->SetBit(BIT(16));

      new ((*fParticles)[iPart]) TParticle(*particle);

      if (!selected || !charged || !primary) continue;
      Double_t pT    = particle->Pt();
      if (pT < 1e-10) continue; /// Along beam axis 
      Double_t pZ    = particle->Pz();
      Double_t theta = TMath::ATan2(pT, pZ);
      Double_t eta   = -TMath::Log(TMath::Tan(theta/2));
      fHEta->Fill(eta);
    }
    return true;
  }
  /** 
   * Do final event processing (fill output)
   * 
   */
  void PostEvent()
  {
    fHeader->SetNprimary(fStack->GetNprimary());
    fHeader->SetNtrack(fStack->GetNtrack());

    fTree->Fill();
    
    fStack->FinishEvent();
    fHeader->SetStack(fStack);
    
    fRunLoader->TreeE()->Fill();
    fRunLoader->WriteKinematics("OVERWRITE");
  }    
  /** 
   * Process one event 
   * 
   * @param iEv Event number 
   * 
   * @return true on success, false otherwize 
   */
  Bool_t Process(Long64_t iEv)
  {
    // --- The stopwatch ---------------------------------------------
    TStopwatch timer;
    timer.Start();
    PreEvent(iEv);
    fHTime->Fill(1, timer.RealTime());
    
    // --- Generate event --------------------------------------------
    timer.Start();
    fGenerator->Generate();
    fHTime->Fill(2, timer.RealTime());

    // --- Process the header ----------------------------------------
    timer.Start();
    Bool_t selected = ProcessHeader();
    fHTime->Fill(3, timer.RealTime());
    
    // --- Loop over particles ---------------------------------------
    timer.Start();
    ProcessParticles(selected);
    fHTime->Fill(4, timer.RealTime());

    // --- Do final stuff --------------------------------------------    
    timer.Start();
    PostEvent();
    fHTime->Fill(5, timer.RealTime());

    return true;
  }
  /** 
   * Finalize this sub-job  
   * 
   */
  void SlaveTerminate()
  {
    fGenerator->FinishRun();
    fRunLoader->WriteHeader("OVERWRITE");
    fGenerator->Write();
    fRunLoader->Write();

    if (fFile) {
      if (fProofFile) {
	fOutput->Add(fProofFile);
	fOutput->Add(new TH1F("filename", fFileName.Data(),1,0,1));
      }
      // Flush out tree 
      fFile->cd();
      fTree->Write(0, TObject::kOverwrite);
      fFile->Close();
      fFile->Delete();
      fFile = 0;
    }
  }
  /** 
   * Final processing of the data 
   * 
   */
  void Terminate()
  {
    if (gProof) gProof->ClearFeedback();
    
    if (!fList)
      fList = static_cast<TList*>(fOutput->FindObject("histograms"));
    if (!fList) {
      Error("Terminate", "No output list");
      return;
    }
    
    if (!fProofFile) {
      TObject* fn = fOutput->FindObject("filename");
      if (fn) fFileName  = fn->GetTitle();
      fProofFile =
	static_cast<TProofOutputFile*>(fOutput->FindObject(FileName()));
    }
    if (fProofFile) 
      fFile = fProofFile->OpenFile("UPDATE");
    if (!fFile)
      fFile = TFile::Open(FileName(),"UPDATE");
    
	
    fHEta  = static_cast<TH1*>(fList->FindObject("dNdeta"));
    fHIpz  = static_cast<TH1*>(fList->FindObject("ipZ"));
    fHType = static_cast<TH1*>(fList->FindObject("type"));
    fHCent = static_cast<TH1*>(fList->FindObject("cent"));
    fHB    = static_cast<TH1*>(fList->FindObject("b"));
    fHPhiR = static_cast<TH1*>(fList->FindObject("phiR"));
    fHTime = static_cast<TH1*>(fList->FindObject("timing"));

    if (!(fHEta && fHIpz && fHType && fHB && fHPhiR && fHTime)) {
      Warning("Terminate", "Missing histograms (%p,%p,%p,%p,%p,%p)",
	      fHEta, fHIpz, fHType, fHB, fHPhiR, fHTime);
      return;
    }

    Int_t nTotal = fHIpz->GetEntries();
    fHEta ->Scale(1./nTotal, "width");
    fHB   ->Scale(1./nTotal, "width");
    fHPhiR->Scale(1./nTotal, "width");
    fHTime->Scale(1./nTotal, "width");

    if (!fFile){
      Warning("Terminate", "No file to write to");
      return;
    }

    fHEta ->Write();
    fHIpz ->Write();
    fHType->Write();
    fHCent->Write();
    fHB   ->Write();
    fHPhiR->Write();
    fHTime->Write();

    fTree = static_cast<TTree*>(fFile->Get("T"));
    if (!fTree)  Warning("Terminate", "No tree");
    
    fFile->Close();
  }
  /** 
   * Interface version used 
   * 
   * @return 1
   */
  Int_t Version() const { return 1; }

  /**
   * @{ 
   * @name Parameters 
   */
  TString  fEGName;               // Name of event generator
  Int_t    fRunNo;                // Run to simulate 
  Double_t fBMin;                 // Least impact parameter 
  Double_t fBMax;                 // Largest impact parameter
  TObject* fGRP;                  //! GRP in one line
  Long64_t fNEvents;              //  Number of requested events
  /* @} */
  /** 
   * @{ 
   * @name ALICE EG interface 
   */
  AliGenerator* fGenerator;       //! Event generator
  AliRunLoader* fRunLoader;       //! Loader of trees
  AliStack*     fStack;           //! Stack of particles
  AliHeader*    fHeader;          //! Header handler
  /* @} */
  /** 
   * @{ 
   * @name Custom output 
   */
  TTree*        fTree;            //! Custom tree 
  TClonesArray* fParticles;       //! List of particles
  /**
   * @{ 
   * @name Diagnostics 
   */
  TList* fList;                   //! List of outputs
  TH1*   fHEta;                   //! dN/deta
  TH1*   fHIpz;                   //! IPz histogram
  TH1*   fHType;                  //! Event type histogram
  TH1*   fHCent;                  //! Event type histogram
  TH1*   fHB;                     //! B histogram
  TH1*   fHPhiR;                  //! Reaction plane
  TH1*   fHTime;                  //! Timing 
  /* @} */
  /**
   * @{ 
   * @name Output files 
   */
  TProofOutputFile* fProofFile;   //! Proof output file 
  TFile*            fFile;        //! Output file
  mutable TString   fFileName;    //! Output file name 
  /* @} */

  // Hide from CINT 
#ifndef __CINT__
  struct ShortHeader {
    UInt_t   fRunNo;
    UInt_t   fEventId;
    UInt_t   fNpart;
    UInt_t   fNbin;
    UInt_t   fType;
    Double_t fIpX;
    Double_t fIpY;
    Double_t fIpZ;
    Double_t fB;
    Double_t fC;
    Double_t fPhiR;
  } fShortHead;
#endif
  /** 
   * Run this selector as a normal process
   * 
   * @param nev        Number of events
   * @param run        Run number to anchor in
   * @param gen        Generator 
   * @param bMin       Least impact parameter [fm]
   * @param bMax       Largest impact parameter [fm]
   * @param monitor    Monitor frequency [s]
   * 
   * @return true on succes
   */
  static Bool_t  Run(Long64_t    nev,
		     UInt_t      run,
		     const char* gen,
		     Double_t    bMin,
		     Double_t    bMax,
		     Int_t       monitor)
  {
    TStopwatch stopwatch;
    stopwatch.Start();
    
    FastSim* sim = new FastSim(gen,run,bMin,bMax,nev);
    sim->Begin(0);
    sim->SlaveBegin(0);

    TTimer* timer = 0;
    if (monitor > 0) {
      // timer = new TTimer(new FastMonitor(sim), monitor*1000,true);
      timer = new TTimer(1000);
      timer->Connect("Timeout()","FastMonitor",
		     new FastMonitor(sim), "Handle()");
      ::Info("Run", "Turning on monitoring");
      timer->Start(-1,false);
    }
      
    for (Long64_t i=0; i <nev; i++) {
      Printf("=== Event # %6lld/%6lld ==========================",
	     i+1, nev);
      sim->Process(i);
      if (timer && (i > 0) && (i % 500 == 0)) {
	if (timer->CheckTimer(gSystem->Now()))
	  Printf("Fired timer");
      }
    }
    if (timer) timer->TurnOff();
    sim->SlaveTerminate();
    sim->Terminate();

    stopwatch.Print();
    return true;
  }
  static void ProofLoadLibs()
  {
    if (!gProof) return;

    // Remember to copy changes to RunFast.C
    TList clsLib;
    clsLib.Add(new TNamed("TVirtualMC",              "libVMC"));
    clsLib.Add(new TNamed("TLorentzVector",          "libPhysics"));
    clsLib.Add(new TNamed("TLinearFitter",           "libMinuit"));
    clsLib.Add(new TNamed("TTree",                   "libTree"));
    clsLib.Add(new TNamed("TProof",                  "libProof"));
    clsLib.Add(new TNamed("TGFrame",                 "libGui"));
    clsLib.Add(new TNamed("TSAXParser",              "libXMLParser"));
    clsLib.Add(new TNamed("AliVEvent",               "libSTEERBase"));
    clsLib.Add(new TNamed("AliESDEvent",             "libESD"));
    clsLib.Add(new TNamed("AliAODEvent",             "libAOD"));
    clsLib.Add(new TNamed("AliAnalysisManager",      "libANALYSIS"));
    clsLib.Add(new TNamed("AliCDBManager",           "libCDB"));
    clsLib.Add(new TNamed("AliRawVEvent",            "libRAWDatabase"));
    clsLib.Add(new TNamed("AliHit",                  "libSTEER"));
    clsLib.Add(new TNamed("AliGenMC",                "libEVGEN"));
    clsLib.Add(new TNamed("AliFastEvent",            "libFASTSIM"));

    TIter next(&clsLib);
    TObject* obj = 0;
    while ((obj = next())) {
      gProof->Exec(Form("gROOT->LoadClass(\"%s\",\"%s\");",
			obj->GetName(), obj->GetTitle()));
    }
  }
  /** 
   * Run this selector in PROOF(Lite)
   * 
   * @param url        Proof URL
   * @param nev        Number of events
   * @param run        Run number to anchor in
   * @param gen        Generator 
   * @param bMin       Least impact parameter [fm]
   * @param bMax       Largest impact parameter [fm]
   * @param monitor    Monitor frequency [s]
   * @param opt        Compilation options
   * 
   * @return true on succes
   */
  static Bool_t Proof(const char*  url,
		      Long64_t     nev,
		      UInt_t       run,
		      const char*  gen,
		      Double_t     bMin,
		      Double_t     bMax,
		      Int_t        monitor=-1,
		      const char*  opt="")
  {
    Printf("# events:  %lld", nev);
    Printf("Run #:     %u",   run);
    Printf("Generator: %s",   gen);
    Printf("b range:   %5.1f-%5.1f", bMin, bMax);
    Printf("monitor:   %ds", monitor);
	   
    TStopwatch timer;
    timer.Start();
    
    TProof::Reset(url);
    TProof::Open(url);
    gProof->ClearCache();

    TString ali = gSystem->ExpandPathName("$(ALICE_ROOT)");
    // TString fwd = gSystem->ExpandPathName("$ANA_SRC");
    TString fwd = ali + "/PWGLF/FORWARD/analysis2";

    gProof->AddIncludePath(Form("%s/include", ali.Data()));
    ProofLoadLibs();
    gProof->Load(Form("%s/sim/GRP.C",fwd.Data()), true);
    gProof->Load(Form("%s/sim/BaseConfig.C",fwd.Data()), true);
    gProof->Load(Form("%s/sim/EGConfig.C",fwd.Data()), true);

    // gROOT->ProcessLine("gProof->SetLogLevel(5);");
    gProof->Load(Form("%s/sim/FastSim.C+%s", fwd.Data(), opt),true);

    if (monitor <= 0) gProof->SetParameter("NOMONITOR", true/*ignored*/);
    else              gProof->SetParameter("PROOF_FeedbackPeriod",
					   monitor*1000/*ms*/);

    FastSim* sim = new FastSim(gen,run,bMin,bMax,nev);
    gProof->Process(sim, nev, "");

    timer.Print();
    return true; // status >= 0;
  }
  ClassDef(FastSim,1); 
};

#endif
//
// EOF
// 
