// EventTimeTask.C
#ifndef __CINT__
# include <TTree.h>
# include <TError.h>
# include <TChain.h>
#else
class TTree;
// Force load of libGui
class TGWindow;
#endif
class AliESDEvent;

/**
 * Data structure to be filled by task - one for each event
 * 
 * @ingroup pwglf_forward_eventtime
 */
struct EventTimeData 
{
  /**
   * Widths of field in the full timestamp 
   */
  enum { 
    kBCWidth      = 12, 
    kOrbitWidth   = 24, 
    kPeriodWidth  = 28
  };
  /** Full time stamp */
  ULong64_t fFull;
  /** LDC creation time - standard time_t */
  UInt_t    fTime;
  /** Mask of detectors present in the event - FMD is 0x1000 */
  UInt_t    fDetectors;
  /** Type of event - 7 is physics */
  UShort_t  fType;
  /** GUID */
  UChar_t   fGUID[42];
  /** 
   * Create a branch in a tree 
   * 
   * @param tree Tree to create the branch in 
   */
  void CreateBranch(TTree* tree)
  {
    tree->Branch("event", &(this->fFull), 
		 "full/l:time/i:detector:type/s:guid/C");
  }
  /** 
   * Set the address of a branch for reading back objects from the tree
   * 
   * @param tree Tree to read back from 
   */
  void ReadBranch(TTree* tree) 
  {
    tree->SetBranchAddress("event", &(this->fFull));
  }
  /** 
   * Utility function to encode the full time stamp from components. 
   * 
   * @param period Period counter (overflow of orbit counter)
   * @param orbit  Orbit counter (period of 88us)
   * @param bc     Bunch crossing number (period of 25ns)
   * 
   * @return Encoded full time stamp 
   */
  static ULong64_t EncodeFull(ULong64_t period, 
			      ULong64_t orbit, 
			      ULong64_t bc) 
  {
    const ULong64_t bcMask     = (1 << kBCWidth) - 1;
    const ULong64_t orbitMask  = (1 << kOrbitWidth) - 1;
    const ULong64_t periodMask = (1 << kPeriodWidth) - 1;
    ULong64_t ret = ((bc     & bcMask)                  | 
		     (orbit  & orbitMask)  <<  kBCWidth | 
		     (period & periodMask) << (kBCWidth+kOrbitWidth));
    return ret;
  }
  /** 
   * Fill information from ESD into this data structure. 
   * 
   * @param esd  Event 
   * @param dets List of active detectors in this event. 
   * @param guid Current file GUID
   */
  void Fill(AliESDEvent* esd, UInt_t dets, const TString& guid);
};

#ifndef NO_TASK
# ifndef __CINT__
#  include <AliESDEvent.h>
# endif
inline void EventTimeData::Fill(AliESDEvent* esd, UInt_t dets, 
				const TString& guid)
{
  ULong64_t period = esd->GetPeriodNumber();
  ULong64_t orbit  = esd->GetOrbitNumber();
  ULong64_t bc     = esd->GetBunchCrossNumber();
  fType            = esd->GetEventType();
  fTime            = esd->GetTimeStamp();//LDC time
  fDetectors       = dets; // esd->GetDAQDetectorPattern();
  fFull            = EncodeFull(period, orbit, bc);
  Int_t i = 0;
  for (i = 0; i < guid.Length() && i < 42; i++) fGUID[i] = guid[i];
  for (; i < 42; i++) fGUID[i] = '\0';
  fGUID[41] = '\0';
}
# else
inline void EventTimeData::Fill(AliESDEvent*, UInt_t, const TString&) 
{
  Warning("Fill", "Calling empty method - shouldn't happen");
}
#endif

#ifndef NO_TASK
# ifndef __CINT__
#  include <AliAnalysisManager.h>
#  include <AliVEventHandler.h>
#  include <AliESDEvent.h>
#  include <TTree.h>
#  include <TH2.h>
#  include <TList.h>
#  include <AliCDBManager.h>
#  include <AliCDBEntry.h>
#  include <AliTriggerConfiguration.h>
#  include <AliTriggerClass.h>
#  include <AliTriggerCluster.h>
#  include <AliDAQ.h>
#  include <TObjArray.h>
#  include <TDirectory.h>
#  include <TUrl.h>
# else
class AliAnalysisManager;
class TTree;
class AliTimeStamp;
class TList;
class TH2;
# endif
# include <AliAnalysisTaskSE.h>
# include <vector>
 
/**
 * A task to record the unique timestamp of each event. 
 *
 * @par Input: ESD 
 * @par Output:  A tree with a single branch 
 * @ingroup pwglf_forward_eventtime
 */
class EventTimeTask : public AliAnalysisTaskSE
{
public:
  enum { 
    kListSlot = 1,
    kTreeSlot = 2
  };
  /** 
   * Default CTOR - for I/O only.
   */
  EventTimeTask() 
    : AliAnalysisTaskSE(), 
      fTree(0),
      fHistograms(0),
      fDetVsType(0)
  {
  }
  /** 
   * User constructor 
   * 
   * @param name Name of task 
   */
  EventTimeTask(const char* name) 
    : AliAnalysisTaskSE(name), 
      fTree(0),
      fHistograms(0),
      fDetVsType(0),
      fGUID("")
  {
    DefineOutput(kListSlot, TList::Class());
    DefineOutput(kTreeSlot, TTree::Class());
    fBranchNames = "ESD:AliESDRun.,AliESDHeader.";
  }
  /** 
   * Create user output objects.
   *
   * Called on each slave at start of job
   */
  void UserCreateOutputObjects()
  {
    Printf("Creating tree and histogram");
    fGUID = "";
    fHistograms = new TList();
    fHistograms->SetOwner();
    fHistograms->SetName("L");
    
    fDetVsType = new TH2D("detVsType", "Detector vs type", 
		     16, -.5, 15.5, 31, -.5, 30.5);
    fDetVsType->SetXTitle("Type");
    fDetVsType->SetYTitle("Detector");
    fDetVsType->SetDirectory(0);
    fHistograms->Add(fDetVsType);
    Printf("Histogram (%d,%f,%f)x(%d,%f,%f)",
	   fDetVsType->GetXaxis()->GetNbins(), 
	   fDetVsType->GetXaxis()->GetXmin(), 
	   fDetVsType->GetXaxis()->GetXmax(),
	   fDetVsType->GetYaxis()->GetNbins(), 
	   fDetVsType->GetYaxis()->GetXmin(), 
	   fDetVsType->GetYaxis()->GetXmax());

    // TDirectory* savdir = gDirectory;
    // Printf("Opening file at slot %d", kTreeSlot);
    // OpenFile(kTreeSlot);
    Printf("Make tree and disassociate from file");
    fTree = new TTree("T", "T");
    fTree->SetDirectory(0);
    Printf("Create branch");
    fData.CreateBranch(fTree);
    // savdir->cd();
    

    PostData(kListSlot, fHistograms);
    PostData(kTreeSlot, fTree);
  }
  /** 
   * Analyse a single event 
   * 
   */
  void UserExec(Option_t*)
  {
    static Bool_t first = true;
    LoadBranches();

    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!esd) return;
    if (!esd->GetHeader()) return;     

    if (first) {
      LoadTriggerConfig(esd->GetRunNumber());
      first = false;
    }
    ULong64_t mask = esd->GetTriggerMask();
    UInt_t    dets = 0;
    for (UShort_t i = 0; i < fDets.size(); i++) {
      if ((1 << i) & mask) dets |= fDets[i];
    }
    // Printf("Event mask 0x%016llx -> 0x%08x", mask, dets);
    fData.Fill(esd, dets, fGUID);
    fTree->Fill();

    UInt_t type      = esd->GetEventType();
    UInt_t detectors = dets;
    for (UInt_t i = 0; i < 31; i++) { 
      if ((1 << i) & detectors) fDetVsType->Fill(type, i);
    }

    PostData(kListSlot, fHistograms);
    PostData(kTreeSlot, fTree);
  }
  Bool_t UserNotify()
  {
    fGUID = "";
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    AliVEventHandler*   inp = mgr->GetInputEventHandler();
    if (!inp) { 
      Warning("UserNotify", "No input handler");
      return true;
    }
    TTree* tree = inp->GetTree();
    if (!tree) { 
      Warning("UserNotify", "No input tree");
      return true;
    }
    TFile* file = tree->GetCurrentFile();
    if (!file) { 
      Warning("UserNotify", "No current file for tree");
      return true;
    }
    const TUrl* url = file->GetEndpointUrl();
    if (!url) { 
      Warning("UserNotify", "No end point for file");
      return false;
    }
    fGUID = gSystem->BaseName(url->GetFile());
    Printf("Got GUID=%s from %s", fGUID.Data(), url->GetUrl());
    
    return true;
  }
  void LoadTriggerConfig(Int_t runNo)
  {
    Printf("Loading trigger configuration for run %d",  runNo);
    // --- Connect to CDB --------------------------------------------
    AliCDBManager* cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorageFromRun(runNo);
    cdb->SetRun(runNo);
    
    // --- Get entry -------------------------------------------------
    AliCDBEntry* entry = cdb->Get("GRP/CTP/Config");
    if (!entry || !entry->GetObject()) { 
      Warning("LoadTriggerConfig", "Couldn't get trigger configuration");
      return;
    }
    AliTriggerConfiguration* config = 
      static_cast<AliTriggerConfiguration*>(entry->GetObject());

    // --- Get the classes, and resize cache -------------------------
    const TObjArray& clss = config->GetClasses();
    fDets.resize(clss.GetEntries());

    // --- Loop over configurations ----------------------------------
    TIter next(&clss);
    AliTriggerClass* cls = 0;
    while ((cls = static_cast<AliTriggerClass*>(next()))) {
      Int_t              mask    = cls->GetMask();
      AliTriggerCluster* cluster = cls->GetCluster();
      if (!cluster) { 
	Warning("LoadTriggerConfig", 
		"Failed to get trigger cluster for %s", cls->GetName());
	continue;
      }
      TString names = cluster->GetDetectorsInCluster();
      if (names.IsNull()) { 
	Warning("LoadTriggerConfig", "No detectors for cluster %s",
		cls->GetName());
	continue;
      }
      UInt_t dets = AliDAQ::DetectorPattern(names);
      UInt_t id   = UInt_t(TMath::Log2(mask));
      fDets[id]   = dets;
      Printf("Trigger mask 0x%08x (%3d): 0x%08x (%s)", 
	     mask, id, dets, names.Data());
    }
    // cdb->Destroy();
  }
  /** 
   * Register with manager and connect output containers
   * 
   * @param mgr Manager 
   */
  void Connect(AliAnalysisManager* mgr)
  {
    mgr->AddTask(this);
    mgr->ConnectInput(this, 0, mgr->GetCommonInputContainer());
    AliAnalysisDataContainer* contTree = 
      mgr->CreateContainer("T", TTree::Class(), 
			   AliAnalysisManager::kOutputContainer,
			   "time.root");
    AliAnalysisDataContainer* contList = 
      mgr->CreateContainer("L", TList::Class(), 
			   AliAnalysisManager::kOutputContainer,
			   "hist.root");
    mgr->ConnectOutput(this, kTreeSlot, contTree);
    mgr->ConnectOutput(this, kListSlot, contList);
  }
  /** 
   * Create an instance of this task, and register it and it's
   * outputs.
   */
  static void Create()
  {
    EventTimeTask* task = new EventTimeTask("eventTime");
    task->Connect(AliAnalysisManager::GetAnalysisManager());
  }
  TTree*              fTree;       // Our tree
  EventTimeData       fData;       // Our data 
  TList*              fHistograms; // List
  TH2D*               fDetVsType;  // Histogram
  std::vector<UInt_t> fDets;       // Per-trigger-bit detector mask 
  TString             fGUID;
  
  ClassDef(EventTimeTask,3);
};
#endif // NO_TASK 

#ifndef NO_MAP
#include <utility>
#include <map>

typedef std::pair<ULong64_t,ULong64_t> EventTimeMapPair;

/**
 * A map of event time-stamp to distance to previous event
 * 
 * @ingroup pwglf_forward_eventtime
 */
struct EventTimeMap : public TObject
{
  /** Map type */
  typedef std::map<ULong64_t,ULong64_t> Map;
  /** Key type */
  typedef Map::key_type Key;
  /** Mapped value type */
  typedef Map::mapped_type Value;
  /** Element type */
  typedef Map::value_type Pair;
  /** Iterator type */
  typedef Map::iterator   Iterator;
  /** Constant iterator type */
  typedef Map::const_iterator ConstIterator;
  /** 
   * Constructor
   */
  EventTimeMap() : TObject(), fMap() {}
  /** 
   * Destructor 
   */
  virtual ~EventTimeMap() {}
  /** 
   * Get name of this object - always the same 
   * 
   * @return The string "eventTimeMap"
   */
  const char* GetName() const { return "eventTimeMap"; }
  /** 
   * Element access.  If the key @a k doesn't already exist, it is
   * created
   * 
   * @param k Key 
   * 
   * @return A pair of key and value 
   */
  Value& operator[](const Key& k)
  {
    return fMap[k];
  } 
  /** 
   * Find the element whos key is @a k 
   * 
   * @param k Key to look for
   * 
   * @return Iterator pointing to element, or end of container if not found 
   */
  Iterator Find(const Key& k)
  {
    return fMap.find(k);
  }
  /** 
   * Find the element whos key is @a k 
   * 
   * @param k Key to look for
   * 
   * @return Iterator pointing to element, or end of container if not found 
   */
  ConstIterator Find(const Key& k) const
  {
    return fMap.find(k);
  }
  /** 
   * Get forward iterator pointing beginning of the container
   *
   * @return Iterator to start of container
   */
  Iterator Begin()
  {
    return fMap.begin();
  }
  /** 
   * Get forward iterator pointing beginning of the container
   *
   * @return Iterator to start of container
   */
  ConstIterator Begin() const
  {
    return fMap.begin();
  }
  /** 
   * Get forward iterator pointing just beyond the end of the container
   *
   * @return Iterator just beyond container
   */
  Iterator End()
  {
    return fMap.end();
  }
  /** 
   * Get forward iterator pointing just beyond the end of the container
   *
   * @return Iterator just beyond container
   */
  ConstIterator End() const
  {
    return fMap.end();
  }
  enum { 
    kInvalidTime = 0xFFFFFFFFFFFFFFFF
  };
  /**
   * Get the time difference to previous event from a event with a
   * given time stamp.
   *
   * @param timestamp Time stamp of the event 
   *
   * @return time difference or kInvalidTime if not found 
   */
  Value Get(const Key& timestamp) const 
  {
    ConstIterator i = Find(timestamp);
    if (i == End()) return kInvalidTime;
    return i->second;
  }
  /** 
   * Set the time to previous event for a given event time stamp
   * 
   * @param timestamp Event time stamp 
   * @param diff      Time to previous event 
   */
  void Set(const Key& timestamp, const Value& diff)
  {
    this->operator[](timestamp) = diff;
  }
  ULong64_t Size() const { return fMap.size(); }
  /** Our embeded map */
  Map fMap;
  ClassDef(EventTimeMap,1)
};
#endif // NO_MAP

#ifndef NO_SORTER
# ifndef __CINT__
#  include <TFile.h>
#  include <TArrayL64.h>
#  include <TMath.h>
#  include <TParameter.h>
#  include <TCanvas.h>
#  include <TH1.h>
# else 
class TFile;
class TTree;
class TCanvas;
# endif

/** 
 * A class to sort the tree and generate our timestamp->dT map.
 *
 * @ingroup pwglf_forward_eventtime
 */
struct EventTimeSorter 
{
  TTree*        fTree; 
  EventTimeData fData;
  ULong64_t     fMin;
  ULong64_t     fMax;

  /** 
   * Constructor 
   */
  EventTimeSorter() : fTree(0), fData(), fMin(0xFFFFFFFFFFFFFFFF), fMax(0) {}
  /** 
   * Progress meter 
   * 
   * @param cur    Current entry
   * @param total  Total number of entries
   */
  void Progress(Long64_t cur, Long64_t total) const
  {
    static UShort_t old = 0;
    if (cur == 0) old = 0;
    UShort_t now = UShort_t(100 * (cur + 1 == total ? 1 : 
				   Double_t(cur)/total));
    if (now > old) 
      std::cout << "\rLooping over " << total << " entries: ... "
		<< now << '%' << std::flush;
    if (now >= 100) std::cout << std::endl;
    old = now;
  }
  /** 
   * Connect a tree 
   * 
   * @param inputName 
   * @param treeName 
   * 
   * @return 
   */
  Bool_t OpenInput(const char* inputName, const char* treeName)
  {
    CloseInput();

    // --- Get input -------------------------------------------------
    TChain* chain = new TChain(treeName);
    chain->Add(inputName);
    if (chain->GetListOfFiles()->GetEntries() < 1) { 
      Error("Run", "Failed to add \"%s\" to chain", inputName);
      return false;
    }
    fTree = chain;

    // --- Set branch address ---------------------------------------
    fData.ReadBranch(fTree);

    return true;
  }
  /** 
   * Disconnect tree
   */		
  void CloseInput()
  {
    if (!fTree) return;
    TFile* file = fTree->GetCurrentFile();
    if (file) file->Close();
    fTree = 0;
  }
  /** 
   * Run the sorter 
   * 
   * @param inputName   Input file name 
   * @param outputName  Output file name 
   * @param treeName    Tree name 
   * 
   * @return true on success
   */
  Bool_t Run(const char* inputName, const char* outputName, 
	     const char* treeName="T") 
  {
    if (!OpenInput(inputName, treeName)) return false;

    // --- Loop over the data ----------------------------------------
    Long64_t   nEnt    = fTree->GetEntries();
    Long64_t   n       = 0;
    ULong64_t* values  = new ULong64_t[nEnt];
    UInt_t*    dets    = new UInt_t[nEnt];
    UShort_t*  types   = new UShort_t[nEnt];
    UInt_t*    times   = new UInt_t[nEnt];
    for (Long64_t iEnt = 0; iEnt < nEnt; iEnt++) { 
      Progress(iEnt, nEnt);

      fTree->GetEntry(iEnt);

      if (!(fData.fDetectors & 0x1000)) 
	// No FMD read-out 
	continue;
      if (fData.fType != 7) { 
	// Ignore all but physics events 
	Warning("Run", "Non-PHYSICS event (%d) with FMD in it", fData.fType);
	// continue;
      }

      // Store values 
      Int_t j = n;
      values[j] = fData.fFull;
      dets[j]   = fData.fDetectors;
      types[j]  = fData.fType;
      times[j]  = fData.fTime;
      n++;
    }
    
    // --- Now sort all values ---------------------------------------
    TArrayL64 index(n);
    TMath::Sort(n, values, index.fArray, false);
    
    // Maps the unique time to the distance to previous event 
    EventTimeMap  eventTimeMap;
    ULong64_t     key   = values[index[0]];
    ULong64_t     prev  = 0;
    ULong64_t     dt    = key-prev;
    ULong64_t     min   = 0xFFFFFFFFFFFFFFFF;
    ULong64_t     max   = 0;
    eventTimeMap[key]   = dt; 
    for (Long64_t i = 1; i < n; i++) { 
      Long64_t j        = index[i];
      Long64_t k        = index[i-1];
      key               = values[j];
      prev              = values[k];
      dt                = key - prev;
      if (dt <= 0) {
	Printf("0x%016llx==0x%016llx -> dt=%10llu [%10lld %10lld] "
	       "(det: 0x%08x 0x%08x  type: 0x%04x 0x%04x  time: %08d %08d)",
	       key, prev, dt, j, k, dets[j], dets[k], 
	       types[j], types[k], times[j], times[k]);
	// continue;
      }
      eventTimeMap[key] = dt;
      min               = TMath::Min(dt, min);
      max               = TMath::Max(dt, max);
    }
    std::cout << "Range is [" << min << "," << max << ']' << std::endl;

    TFile* output = TFile::Open(outputName, "RECREATE");
    if (!output) { 
      Error("Run", "Failed to create output file \"%s\"", outputName);
      return false;
    }
    fMin = min;
    fMax = max;
    eventTimeMap.Write();
    output->Write();
    output->Close();
    
    delete [] values; 

    CloseInput();

    return true;
  }
  Bool_t Test(const char* inputName, const char* outputName, 
	      const char* treeName="T") 
  {
    if (!OpenInput(inputName, treeName)) return false;

    // --- Get our map -----------------------------------------------
    TFile* output = TFile::Open(outputName, "UPDATE");
    if (!output) { 
      Error("Test", "Failed to open \"%s\"", outputName);
      return false;
    }

    EventTimeMap* eventTimeMap = 
      static_cast<EventTimeMap*>(output->Get("eventTimeMap"));
    if (!eventTimeMap) { 
      Error("Test", "Failed to get \"%s\" from \"%s\"",
	    "eventTimeMap", outputName);
      return false;
    }

    // --- Histograms --------------------------------------------------
    ULong64_t mmin = TMath::Min(25*fMin, 900000ULL);
    ULong64_t mmax = TMath::Min(25*fMax, 110000000ULL);
    TH1*      diff = new TH1D("diff", "Time-to-last-event (10#mus bins w/FMD)", 
			      (mmax-mmin)/10000, mmin, mmax);
    diff->SetStats(0);
    diff->SetXTitle("#Deltat [ns]");
    diff->SetFillColor(kRed+1);
    diff->SetFillStyle(3001);
    diff->SetDirectory(0);
    

    TH1* ldiff = new TH1D("ldiff", "log(#Deltat) - Events w/FMD",
			  300, 0, 15);
    ldiff->SetStats(0);
    ldiff->SetXTitle("log_{10}(#Deltat) [ns]");
    ldiff->SetFillColor(kGreen+1);
    ldiff->SetFillStyle(3001);
    ldiff->SetDirectory(0);
			 
    // --- Loop over the data ----------------------------------------
    Long64_t   nEnt    = fTree->GetEntries();
    Long64_t   nZero   = 0;
    Long64_t   nMiss   = 0;
    Long64_t   nGood   = 0;
    Long64_t   nNoFMD  = 0;
    for (Long64_t iEnt = 0; iEnt < /*10*/ nEnt; iEnt++) { 
      Progress(iEnt, nEnt);
      fTree->GetEntry(iEnt);
      
      if (!(fData.fDetectors & 0x1000)) {
	// No FMD read-out 
	nNoFMD++;
	continue;
      }
      if (fData.fType != 7) { 
	// Ignore all but physics events 
	Warning("Run", "Non-PHYSICS event (%d) with FMD in it", fData.fType);
	// continue;
      }
      
      // Look-up the timestamp 
      ULong64_t value = fData.fFull;
      ULong64_t dT    = eventTimeMap->Get(value);
      if (dT == EventTimeMap::kInvalidTime) {
	Warning("Test", "Value %llu not found", value);
	ldiff->Fill(1);
	nMiss++;
	continue;
      }
      if (dT <= 0) { 
#if 0
	Warning("Test", "Impossible dt=%llu found for %llu (0x%0x %2d)", 
		dT, value, fData.fDetectors, fData.fType);
#endif
	ldiff->Fill(1);
	nZero++;
	continue;
      }
      nGood++;
      Double_t  dt    = 25.*dT;
      diff->Fill(dt);
      Double_t  logDt = TMath::Log10(dt);
      ldiff->Fill(logDt);
    }
    CloseInput();
    Printf("missing: %llu, bad: %llu, good: %llu, no FMD: %llu, all: %llu", 
	   nMiss,nZero,nGood,nNoFMD,nEnt);

    // --- Draw results ----------------------------------------------
    TCanvas* c = new TCanvas("c", "C");
    c->SetTopMargin(0.01);
    c->SetRightMargin(0.01);
    c->Divide(2,1); // ,0,0);
    TVirtualPad* p = c->cd(2);
    // p->SetRightMargin(0.10);
    p->SetLogy();
    ldiff->DrawCopy();
    
    p = c->cd(1);
    p->SetLogy();
    p->SetLogx();
    diff->DrawCopy();

    c->Print("dt.png");
    c->Print("dt.C");

    // --- Write histogram to file -----------------------------------
    output->cd();
    diff->Write();
    output->Write();
    output->Close();
    
    return true;
  }
};
#endif // NO_SORTER 
#ifdef __MAKECINT__
# ifndef NO_MAP
#  pragma link C++ class std::pair<ULong64_t,ULong64_t>+;
#  pragma link C++ class std::map<ULong64_t,ULong64_t>+;
# endif
# ifndef NO_TASK
#  pragma link C++ class std::vector<UInt_t>+;
# endif
#endif

// EOF
        
        
