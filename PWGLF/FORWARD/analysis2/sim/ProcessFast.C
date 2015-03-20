/**
 * @file   ProcessFast.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Fri Mar 20 12:13:28 2015
 * 
 * @brief This script defines classes for looping over the data
 * produced by FastSim.C
 * 
 */

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
#else
class TTree;
class TH1;
class TH2;
class THStack;
class TCanvas;
class TAxis;
class TStopwatch;
class TLegend;
class TMultiGraph;
class TParticle;
class TArrayI;
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
  void Clear() {
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
struct Processor
{
  /** The tree */
  TTree*        fTree;       //!
  /** Cache of our header */
  Header*       fHeader;     //!
  /** List of particles */
  TClonesArray* fParticles;  //!
  /** Whether to be verbose */
  Bool_t        fVerbose;    //
  /** Number of OK events */
  Long64_t      fOK;         //
  /** A simple timer */
  TStopwatch    fTimer;      //!
  /** Name */
  TString       fName;       //! 
  /** 
   * Constructor.  Opens the file passed and sets internal pointers to
   * tree, header, and particle list.
   * 
   * @param fileName File to read from 
   * @param verb     Whether to be verbose 
   */  
  Processor(const char* fileName, Bool_t verb=false)
    : fTree(0), fHeader(0), fParticles(0), fVerbose(verb), fOK(0),
      fName(fileName)
  {
    Int_t idx = fName.Index("_");
    if (idx != kNPOS)
      fName.Remove(idx, fName.Length()-idx);
			   
    
    TFile* file = TFile::Open(fileName, "READ");
    if (!file) {
      Printf("E: Failed to open \"%s\"", fileName);
      return;
    }
    
    fTree = static_cast<TTree*>(file->Get("T"));
    if (!fTree) {
      Printf("E: Failed to get tree from %s", fileName);
      file->Close();
      return;
    }
    fHeader = new Header;
    fParticles = new TClonesArray("TParticle");
    fTree->SetBranchAddress("header",    &(fHeader->fRunNo));
    fTree->SetBranchAddress("particles", &fParticles);
  }
  /**
   * Destructor 
   */
  virtual ~Processor()
  {
    if (fHeader)    delete fHeader;
    if (fParticles) delete fParticles;
    if (fTree) {
      if (fTree->GetCurrentFile()) 
	fTree->GetCurrentFile()->Close();
      fTree = 0;
    }
  }
  /** 
   * Get the name of this processor.  Must be overloaded to return a
   * string.
   * 
   * @return The name of the processor 
   */
  virtual const char* GetName() const = 0;
  /** 
   * Clear internal caches.  Called at start of each event.
   */
  virtual void Clear()
  {
    fHeader->Clear();
    fParticles->Clear();
  }
  /** 
   * Set-up.  Called before event processing.  Here we should declare
   * histograms and the like.  Thist must be overloaded by derived
   * classes.
   */
  virtual void Setup() = 0;
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
  /** 
   * Process a single event.  
   * 
   * @param event Entry number to get from the tree 
   * 
   * @return True on success, false when there's no more to process 
   */
  virtual Bool_t ProcessEvent(Long64_t event)
  {
    Clear();
    Int_t read = fTree->GetEntry(event);
    if (read <= 0) return false;
    
    if (fVerbose)
      printf("Process event %7lld (0x%x) ",
	     event, fHeader->fType);

    if (!ProcessHeader()) {
      if (fVerbose) Printf("No accepted");
      return true;
    }
    fOK++;

    ProcessParticles();

    return true;
  }
  /** 
   * Final processing.  Must be overloaded to do what we need to do
   * here.
   */
  virtual void Terminate() = 0;
  /** 
   * Report progress
   * 
   * @param event The event number we're looking at. 
   */
  virtual void Progress(Int_t event, Long64_t max)
  {
    if (event == 0) {
      fTimer.Start(true);
      return;
    }
    Bool_t final = event < 0;
    if (!final && (event % 1000) != 0) return;
    Long64_t ntot = (max < 0 ? fTree->GetEntries() : max);
    Double_t elap = fTimer.RealTime();
    Int_t    proj = elap;
    Double_t rati = 0;
    if (!final) {
      Double_t rate = elap / event;
      Double_t ttot = rate*ntot;
      rati          = Double_t(event) / ntot;
      proj          = int(ttot - elap);
    }
    Int_t    hour = proj/3600;
    Int_t    minu = (proj/60)%60;
    Int_t    seco = proj%60;

    if (!final) {
      printf("\r%9d/%9lld (%3d%%) %d:%02d:%02d remaining", event, ntot,
	     int(100*rati), hour, minu, seco);
      fTimer.Continue();
    }
    else {
      printf("\nProcessed %d events in %d:%02d:%02d\n",
	     -event, hour, minu, seco);
      fTimer.Stop();
    }
    fflush(stdout);
  }
  /** 
   * Run the processing of the tree.  
   * 
   * @param max Maximum number of entries to read. 
   */
  virtual Bool_t Run(Long64_t max)
  {
    if (!fTree) {
      Printf("E: No tree opened");
      return false;
    }
    Setup();
    Long64_t nEvents = fTree->GetEntries();
    if (max > 0) nEvents = TMath::Min(nEvents, max);
    Long64_t iEvent = 0;
    while ((max < 0) || (iEvent < max)) {
      if (!ProcessEvent(iEvent)) break;
      Progress(iEvent, nEvents);
      iEvent++;
    }
    Progress(-iEvent, nEvents);
    Terminate();

    return fOK > 0;
  }
  /** 
   * Write stuff to the output.  Must be called after terminate. x
   * 
   * @param dir Output directory. 
   */
  virtual void Write(TDirectory* dir) = 0;
};

//====================================================================
/** 
 * Base class for making @f$ 1/N dN_{ch}/d\eta@f$ 
 */
struct dNdetaProcessor : public Processor
{
  TH1* fdNdeta; //!
  
  dNdetaProcessor(const char* filename, Bool_t verbose=false)
    : Processor(filename, verbose), fdNdeta(0)
  {}
  virtual const char* GetName() const { return "dummy"; }
  /** 
   * Static member function to create a histogram 
   * 
   * @return Newly created histogram
   */  
  static TH1* CreatedNdeta()
  {
    Double_t maxEta = 5; // 10;
    Double_t dEta   = 10./200 * 5;
    TH1* eta = new TH1D("dNdeta", "Charged particle pseudo-rapidity density",
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
   * Create our histogram 
   * 
   */
  virtual void Setup()
  {
    fdNdeta = CreatedNdeta();
    fdNdeta->SetMarkerColor(kBlack);
    fdNdeta->SetMarkerStyle(21);
    fdNdeta->SetTitle(GetName());
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
    if (fOK <= 0) {
      delete fdNdeta;
      fdNdeta = 0;
      Printf("W: No events selected");
      return;
    }
    Printf("A total of %lld %s events", fOK, GetName());
    fdNdeta->Scale(1. / fOK, "width");
  }
  /** 
   * Write histogram to output 
   * 
   * @param dir Output directory 
   */
  virtual void Write(TDirectory* dir)
  {
    TDirectory* out = dir->mkdir(GetName());
    out->cd();
    fdNdeta->Write();
    dir->cd();
  }
  
};

//====================================================================
/** 
 * Select NSD events and builds the @f$ 1/N dN_{ch}/d\eta@f$ 
 * 
 */
struct NSDProcessor : public dNdetaProcessor
{
  /** 
   * Constructor 
   * 
   * @param filename File to open 
   * @param verbose  Whether to be verbose 
   */
  NSDProcessor(const char* filename, Bool_t verbose=false)
    : dNdetaProcessor(filename, verbose)
  {}
  /** 
   * Get the name 
   * 
   * @return The name 
   */
  virtual const char* GetName() const { return "nsd"; }
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
};

//====================================================================
/**
 * Processes INEL events and build the @f$ 1/N dN_{ch}/d\eta@f$ 
 * 
 */
struct INELProcessor : public dNdetaProcessor
{
  /** 
   * Constructor 
   * 
   * @param filename File to open 
   * @param verbose  Whether to be verbose 
   */  
  INELProcessor(const char* filename, Bool_t verbose=false)
    : dNdetaProcessor(filename, verbose)
  {}
  /** 
   * Get the name 
   * 
   * @return The name 
   */
  virtual const char* GetName() const { return "inel"; }
  /** 
   * Process the header.  
   * 
   * @return Always true - by definition INEL is all events. 
   */
  virtual Bool_t ProcessHeader() { return true; }
};

//====================================================================
/**
 * Processes events and build the @f$ 1/N dN_{ch}/d\eta@f$ for each
 * bin in reference multiplicity.  The reference multiplicity is the
 * number of charged particles with @f$|\eta|\le0.8@f$
 */
struct MultProcessor : public dNdetaProcessor
{
  /** Stack of histograms */
  THStack* fStack;
  /** Upper limits on our bins */
  TArrayI  fBins;
  /** Least number of events */
  Int_t fMinEvents;
  /** Reference multiplicity */
  Int_t fRefMult;
  /** 
   * Constructor. 
   * 
   * @param filename File to open
   * @param verbose  Whether to verbose
   */
  MultProcessor(const char* filename, Bool_t verbose=false, Int_t min=1000)
    : dNdetaProcessor(filename, verbose), fStack(0),
      fMinEvents(min), fRefMult(0)
  {
    //              +1 +2 +3 +3  +5, +5, +5, +5,+10,+10,+10,+10,+10,+10,+10
    Int_t bins[] = { 0, 3, 6, 9, 14, 19, 24, 29, 39, 49, 59, 69, 79, 89, 99 };
    fBins.Set(15, bins);
  }
  /** 
   * Get the name 
   * 
   * @return The name 
   */
  virtual const char* GetName() const { return "mult"; }
  /** 
   * Get the color associated with a centrality bin. 
   * 
   * @param low  Low edge 
   * @param high High edge 
   * 
   * @return Color identifier. 
   */
  Int_t GetCentralityColor(Int_t low, Int_t high) const
  {
    UShort_t centLow  = low;
    UShort_t centHigh = high;
    Float_t  fc       = (centLow+double(centHigh-centLow)/2) / 100;
    Int_t    nCol     = gStyle->GetNumberOfColors();
    Int_t    icol     = TMath::Min(nCol-1,int(fc * nCol + .5));
    Int_t    col      = gStyle->GetColorPalette(icol);
    //Info("GetCentralityColor","%3d: %3d-%3d -> %3d",bin,centLow,centHigh,col);
    return col;
  }
  /** 
   * Modify a histogram.  Sets the line, marker, and fill styles and
   * colors, as well as the name and title. 
   * 
   * @param h    Histogram to modify 
   * @param low  Low edge of our bin
   * @param up   High edge of our bin 
   */
  void ModHist(TH1* h, Int_t low, Int_t up)
  {
    Color_t col = GetCentralityColor(low, up);
    h->SetLineColor(col);
    h->SetLineStyle(1);
    h->SetMarkerColor(col);
    h->SetMarkerStyle(24);
    h->SetFillColor(kWhite);
    h->SetFillStyle(0);
    h->SetName(Form("h%03d_%03d", low, up));
    if (low == up) {
      if (low > 0) h->SetTitle(Form("%3d+", up));
      else         h->SetTitle(Form("%2d", up));
    }
    else           h->SetTitle(Form("%2d-%2d", low, up));
  }
  
  /** 
   * Create our histogram and histogram stack 
   */
  void Setup()
  {
    // We use the base class object as a cache 
    dNdetaProcessor::Setup();
    fdNdeta->SetName("cache");
    fdNdeta->SetTitle("Cache");

    // Create our stack 
    fStack = new THStack("pythia", "pythia");
    Int_t last = 0;
    for (Int_t i = 0; i < fBins.GetSize(); i++) {
      Int_t up = fBins[i];
      TH1*  h  = CreatedNdeta();
      ModHist(h, last, up);
      fStack->Add(h);
      last = up+1;
    }
    // Add the final bin that goes to infinity 
    TH1* t = CreatedNdeta();
    ModHist(t, last, last);
    fStack->Add(t);
  }
  /** 
   * Process the header.  Accepts all events. 
   * 
   * @return true 
   */
  virtual Bool_t ProcessHeader() { return true; }
  /** 
   * Clear our internal caches
   */
  virtual void Clear() 
  {
    dNdetaProcessor::Clear();
    fdNdeta->Reset(); // Clear the cache
    fRefMult = 0;
  }
  virtual void Fill(Double_t eta)
  {
    dNdetaProcessor::Fill(eta);
    if (TMath::Abs(eta) <= 0.8) fRefMult++;
  }
  /** 
   * Process a single event.  
   *
   * First we fill the internal cache using the base class methods.
   * Then we find the histogram for this particular reference
   * multiplicity and add our event cache to that bin.  The number of
   * events in each bin is counted in the unique ID of each bin.
   * 
   * @param event Event number 
   * 
   * @return True on success, false when there's no more to process 
   */
  virtual Bool_t ProcessEvent(Long64_t event)
  {
    if (!dNdetaProcessor::ProcessEvent(event)) return false;

    // Make sure we have nothing in the underflow bin. 
    fdNdeta->SetBinContent(0,0);

    // Now integrate between eta=-0.8 and +0.8
    // Int_t    b1   = fdNdeta->GetXaxis()->FindBin(-0.8);
    // Int_t    b2   = fdNdeta->GetXaxis()->FindBin(+0.8);
    // Double_t mult = fdNdeta->Integral(b1, b2);
    Int_t mult = fRefMult;
    
    // Find the histogram to update 
    TH1*     out  = 0; 
    for (Int_t i = 0; i < fBins.GetSize(); i++) {
      if (mult > fBins[i]) continue;
      out = static_cast<TH1*>(fStack->GetHists()->At(i));
      break;
    }
    // If we didn't find any histogram, take the last one 
    if (!out) out = static_cast<TH1*>(fStack->GetHists()->Last());

    // If we still have no histogram, return immediately 
    if (!out) return true;

    // Printf("Integral (%d,%d): %f -> %s", b1, b2, mult, out->GetTitle());
    // Add our cache to the appropriate bin 
    out->Add(fdNdeta);
    // Increment the event count for this bin in the unique id
    out->SetUniqueID(out->GetUniqueID()+1);

    return true;
  }
    
  /** 
   * Final processing. 
   * 
   * Normalize each bin to the number of events in each bin (stored in
   * the unique ID of that bin).
   */
  virtual void Terminate()
  {
    if (fOK <= 0) {
      delete fdNdeta;
      fdNdeta = 0;
      delete fStack;
      fStack = 0;
      Printf("W: No events selected");
      return;
    }
    Printf("A total of %lld Mult events", fOK);

    TList*    hists = fStack->GetHists();
    TObjLink* link  = hists->FirstLink();
    while (link) {
      TObject* o = link->GetObject();
      if (!o) {
	link = link->Next();
	continue;
      }
      Int_t n = o->GetUniqueID();
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
      TH1* h = static_cast<TH1*>(o);
      h->Scale(1. / n, "width");
      Printf(" scaled");
      link = link->Next();
    }
#if 0
    // Loop over each bin 
    TIter next(fStack->GetHists());
    TH1* h = 0;
    while ((h = static_cast<TH1*>(next()))) {
      Int_t n = h->GetUniqueID();
      Printf("%9d events in bin %s", n, h->GetTitle());
      h->Scale(1. / n, "width");
    }
#endif
  }
  /** 
   * Write or data to the output 
   * 
   * @param dir 
   */
  virtual void Write(TDirectory* dir)
  {
    TDirectory* out = dir->mkdir(GetName());
    out->cd();
    fStack->Write();
    dir->cd();
  }
};

//====================================================================
/**
 * A steering class for running processors. 
 */
struct Master 
{
  /** The processor to use */
  dNdetaProcessor* fProcessor;
  /**
   * Constructor 
   */
  Master() : fProcessor(0) {}
  /** 
   * Run this job
   * 
   * @param maxEvents     Maximum number of events
   * @param type          What type of processor 
   * @param ppFileName    Input file name 
   * @param dataFileName  Optional data file name 
   * @param title         Title of job
   */   
  void Run(Long64_t maxEvents,
	   UShort_t    type, 
	   const char* ppFileName,
	   const char* dataFileName,
	   const char* title)
  {
    Printf("pp file: %s", ppFileName);
    // Create our processor 
    switch (type) {
    case 0: fProcessor = new INELProcessor(ppFileName); break;
    case 1: fProcessor = new NSDProcessor(ppFileName); break;
    case 2: fProcessor = new MultProcessor(ppFileName); break;
    default:
      Error("Run", "Unknown type: %d", type);
      return;
    }

    // Run the processor 
    if (!fProcessor->Run(maxEvents)) return;

    // Prepare for the output 
    TString outName(title);
    outName.ReplaceAll("@", "");
    outName.ReplaceAll(" ", "_");
    outName.ReplaceAll("__", "_");
    outName.ReplaceAll("#sqrt{#it{s}}=", "");

    // Write the output 
    TFile* out = TFile::Open(Form("%s.root", outName.Data()), "RECREATE");
    fProcessor->Write(out);
    out->Write();
    Printf("Saved in %s", out->GetName());
    
    TMultiGraph* data = 0;
    if (dataFileName && dataFileName[0] != '\0') {
      TFile* dataFile = TFile::Open(dataFileName, "READ");
      if (dataFile) {
	data = static_cast<TMultiGraph*>(dataFile->Get("all"));
	if (data) {
	  TList* l = data->GetListOfGraphs();
	  TObject* old = l->FindObject(Form("alice_%s", fProcessor->GetName()));
	  if (old) l->Remove(old);
	}
	dataFile->Close();
      }
    }
    DrawProcessor(data, title, outName);
  }
  /** 
   * Draw the result.  Delegates to appropriate member function. 
   * 
   * @param data      Optional data results 
   * @param title     Title 
   * @param outName   Output name 
   */
  void DrawProcessor(TMultiGraph* data,
		     const char*  title,
		     const char*  outName)
  {
    MultProcessor* mult = dynamic_cast<MultProcessor*>(fProcessor);
    if (mult) DrawProcessor(mult,       data, title, outName);
    else      DrawProcessor(fProcessor, data, title, outName);
  }
  /** 
   * Set some styles on a histogram 
   * 
   * @param h 
   */
  void SetStyle(TH1* h)
  {
    h->GetXaxis()->SetTitleFont(42);
    h->GetYaxis()->SetTitleFont(42);
    h->GetZaxis()->SetTitleFont(42);
    h->GetXaxis()->SetLabelFont(42);
    h->GetYaxis()->SetLabelFont(42);
    h->GetZaxis()->SetLabelFont(42);
    h->GetYaxis()->SetTitleOffset(1.2);
  }
  /** 
   * Draw a 'normal processor'
   *  
   * @param p      Processor 
   * @param data   Optional data
   * @param title  Title 
   * @param out    Output name 
   */
  void DrawProcessor(dNdetaProcessor* p,
		     TMultiGraph*     data,
		     const char*      title,
		     const char*      out)
  {
    Printf("Drawing a simple dNdetaProcessor");
    TH1*     pp = p->fdNdeta;
    if (!pp) return;
    
    TCanvas* c  = MakeCanvas("dndeta");

    SetStyle(pp);
    pp->SetStats(0);
    pp->SetTitle(title);
    pp->Smooth(3);
    pp->SetLineWidth(2);
    pp->SetLineColor(kBlack);
    pp->Draw("hist l");
    pp->SetMaximum(1.2*pp->GetMaximum());
    pp->SetMinimum(0);
    if (data) DrawData(data, pp);

    c->SetName(out);    
    PrintCanvas(c);
  }
  /** 
   * Draw a 'normal processor'
   *  
   * @param p      Processor 
   * @param data   Optional data
   * @param title  Title 
   * @param out    Output name 
   */
  void DrawProcessor(MultProcessor*   p,
		     TMultiGraph*     data,
		     const char*      title,
		     const char*      out)
  {
    Printf("Drawing a mult dNdetaProcessor");
    THStack* pp = p->fStack;
    if (!pp) return;
    
    TCanvas* c  = MakeCanvas("dndeta");
    pp->SetTitle("");
    pp->Draw("nostack" /*" hist l"*/);
    pp->GetHistogram()->SetXTitle("#it{#eta}");
    pp->GetHistogram()->SetYTitle("1/#it{N} d#it{N}_{ch}/d#it{#eta}");
    pp->SetMaximum(1.2*pp->GetMaximum("nostack"));
    pp->SetMinimum(-.5);
    if (data) DrawData(data, p->fdNdeta);

    TLegend* l = c->BuildLegend(.7, .6, .95, .95, title);
    l->SetFillColor(0);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextFont(42);
    
    c->SetName(out);    
    PrintCanvas(c);
  }
  /** 
   * Draw (optional) data results 
   * 
   * @param data List of graphs 
   */
  void DrawData(TMultiGraph* data, TH1* pp)
  {
    if (!data) return;

    data->Draw("p");
    TLegend* l = new TLegend(0.2, 0.12, .95, .4);
    l->SetNColumns(1);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
    TLegendEntry* e = l->AddEntry("dummy", fProcessor->fName.Data(), "l");
    e->SetLineColor(pp->GetLineColor());
    e->SetLineWidth(pp->GetLineWidth());
    e->SetMarkerColor(pp->GetMarkerColor());
    e->SetMarkerSize(pp->GetMarkerSize());
    e->SetMarkerStyle(pp->GetMarkerStyle());
    TIter nextD(data->GetListOfGraphs());
    TGraph* g = 0;
    while ((g = static_cast<TGraph*>(nextD()))) {
      TString name(g->GetName());
      if (name.Contains("Forward"))
	e = l->AddEntry("dummy", "PWGLF/GEO work-in-progess", "pl");
      else 
	e = l->AddEntry("dummy", "PWGUD/MULT work-in-progess", "pl");
      e->SetLineColor(g->GetLineColor());
      e->SetLineWidth(g->GetLineWidth());
      e->SetMarkerColor(g->GetMarkerColor());
      e->SetMarkerSize(g->GetMarkerSize());
      e->SetMarkerStyle(g->GetMarkerStyle());
    }
    l->Draw();    
  }
  /** 
   * utility to make a canvas. 
   * 
   * @param name Sur-name of canvas 
   * 
   * @return The created canvas 
   */
  TCanvas* MakeCanvas(const char* name)
  {
    const char* egPP  = fProcessor->fName.Data();
    
    TCanvas* c = new TCanvas(Form("%s_%s", egPP, name),
			     "Canvas", 800, 900);
    c->SetTopMargin(0.01);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.02);
    return c;
  }
  /** 
   * Print the canvas to a (set of) file(s)
   * 
   * @param c     Canvas to print
   * @param types Space separated list of types 
   */
  void PrintCanvas(TCanvas* c, const char* types="png")
  {
    c->Modified();
    c->Update();
    c->cd();
    TString    t(types);
    TObjArray* tokens = t.Tokenize(" ,+:");
    TIter      next(tokens);
    TObject*   token = 0;
    while ((token = next())) {
      c->Print(Form("%s.%s", c->GetName(), token->GetName()));
    }
    if (tokens) delete tokens;
  }
};

//====================================================================
void
ProcessFast(Long64_t    maxEvents,
	    Int_t       type, 
	    const char* ppFileName,
	    const char* dataFileName=0,
	    const char* title=0)
{

  Master* m = new Master();
  m->Run(maxEvents, type, ppFileName, dataFileName, title);
}

//
//  EOF
//  
	  
