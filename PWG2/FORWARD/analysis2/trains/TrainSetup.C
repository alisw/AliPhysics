/**
 * @defgroup pwg2_forward_trains Trains
 * 
 * @ingroup pwg2_forward
 */
/**
 * @file   TrainSetup.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 12:12:00 2011
 * 
 * @brief  
 * 
 * @ingroup pwg2_forward_trains
 * 
 */

#ifndef __CINT__
#include <fstream>
#include <iostream>

#include <TAlienCollection.h>
#include <TArrayI.h>
#include <TChain.h>
#include <TDatime.h>
#include <TEnv.h>
#include <TFile.h>
#include <TGrid.h>
#include <TList.h>
#include <TObjString.h>
#include <TProof.h>
#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TROOT.h>

#include <AliAODHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisDataContainer.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisAlien.h>
#include <AliESDInputHandler.h>
#include <AliMCEventHandler.h>
#include <AliVEventHandler.h>
#include <AliPhysicsSelection.h>
#include <AliCentralitySelectionTask.h>
#else
class TArrayI;
class TChain;
class AliAnalysisManager;
#endif

//====================================================================
/** 
 * Generic set-up of an analysis train using the grid-handler (AliEn plugin). 
 * 
 * Users should define a class that derives from this.  The class
 * should implement the member function CreateTasks to add needed
 * tasks to the train
 * 
 * @code 
 * // MyTrain.C 
 * class MyTrain : public TrainSetup
 * {
 * public:
 *   MyTrain(Bool_t   dateTime = false, 
 *           UShort_t year     = 0, 
 *           UShort_t month    = 0, 
 *           UShort_t day      = 0, 
 *           UShort_t hour     = 0, 
 *           UShort_t min      = 0) 
 *     : TrainSetup("My train", dateTime, year, month, day, hour, min)
 *   {}
 *   void Run(const char* type, const char* mode, const char* oper, 
 *            Int_t nEvents=-1, Bool_t mc=false,
 *            Bool_t usePar=false)
 *   {
 *     Exec(type, mode, oper, nEvents, mc, usePar);
 *   }
 * protected:
 *   void CreateTasks(EMode mode, Bool_t par, AliAnalysisManager* mgr)
 *   {
 *     AliAnalysisManager::SetCommonFileName("my_analysis.root");
 *     LoadLibrary("MyAnalysis", mode, par, true);
 *     Bool_t mc = mgr->GetMCtruthEventHandler() != 0;
 *     gROOT->Macro("MyAnalysis.C");
 *   }
 * };
 * @endcode 
 * 
 * This can then be run like 
 * 
 * @verbatim 
 * > aliroot 
 * Root> .L TrainSetup.C 
 * Root> .L MyTrain.C 
 * Root> MyTrain t;
 * Root> t.Run();
 * @endverbatim 
 * 
 * or as a script 
 * 
 * @code 
 * {
 *   gROOT->LoadMacro("TrainSetup.C");
 *   gROOT->LoadMacro("MyTrain.C");
 *   MyTrain t;
 *   t.Run();
 * }
 * @endcode 
 * 
 * To byte compile this, you need to 
 * - load the ROOT AliEn library
 * - load the analysis libraries 
 * - add $ALICE_ROOT/include to header search 
 * first 
 *
 * @verbatim 
 * > aliroot 
 * Root> gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWG2/FORWARD/analysis2:"
 * Root>                          "$ALICE_ROOT/ANALYSIS/macros",
 * Root> 			  gROOT->GetMacroPath()));
 * Root> gSystem->AddIncludePath("-I${ALICE_ROOT}/include");
 * Root> gSystem->Load("libRAliEn");
 * Root> gSystem->Load("libANALYSIS");
 * Root> gSystem->Load("libANALYSISalice");
 * Root> gROOT->LoadMacro("TrainSetup.C+");
 * @endverbatim 
 * 
 * 
 * @ingroup pwg2_forward_trains
 * 
 */
struct TrainSetup
{
  /** 
   * Data type to process 
   */
  enum EType { 
    /** Event Summary Data */
    kESD, 
    /** Analysis Object Data */
    kAOD
  };
  /**
   * How to run the analysis
   * 
   */
  enum EMode {
    /** Locally */ 
    kLocal = 1, 
    /** In PROOF(-Lite) cluster */
    kProof, 
    /** On the Grid */
    kGrid 
  };
  /**
   * What stage of the analysis to run 
   * 
   */
  enum EOper { 
    /** Testing.  Local processing, a single copied from Grid */
    kTest, 
    /** Off-line */
    kOffline, 
    /** Submit to queue */
    kSubmit, 
    /** Merge and terminate */
    kTerminate, 
    /** Full run */
    kFull
  };

  //__________________________________________________________________
  /** 
   * Constructor 
   * 
   * @param name         Name of analysis (free-form)
   * @param useDateTime  Whether to append date and time to the name 
   * @param year         Year - if not specified, taken from current date
   * @param month        Month - if not specified, taken from current date 
   * @param day          Day - if not specified, taken from current date 
   * @param hour         Hour - if not specified, taken from current time  
   * @param min          Minute - if not specified, taken from current time  
   */
  TrainSetup(const char* name, Bool_t useDateTime=true, 
	     UShort_t year=0, UShort_t month=0, 
	     UShort_t day=0, UShort_t hour=0, UShort_t min=0) 
    : fName(name),
      fRootVersion("v5-28-00a"),
      fAliRootVersion("v4-21-18-AN"),
      fAliEnAPIVersion("V1.1x"),
      fProofServer("alicecaf.cern.ch"),
      fDataDir("/alice/data/2010/LHC10c"),
      fDataSet("/COMMON/COMMON/LHC09a4_run8100X#/esdTree"),
      fXML(""), 
      fRunNumbers(0),
      fListOfPARs(),
      fListOfSources(),
      fListOfLibraries(),
      fListOfExtras(),
      fNReplica(4),
      fESDPass(3),
      fPassPostfix(""),
      fEscapedName(name),
      fAllowOverwrite(kFALSE),
      fUseGDB(kFALSE)
  {
    char  c[] = { ' ', '/', '@', 0 };
    char* p   = c;
    while (*p) { 
      fEscapedName.ReplaceAll(Form("%c", *p), "_");
      p++;
    }

    if (useDateTime) { 
      if (year == 0 || month == 0 || day == 0) {
	TDatime now;
	year  = now.GetYear();
	month = now.GetMonth();
	day   = now.GetDay();
	hour  = now.GetHour();
	min   = now.GetMinute();
      }
      fEscapedName.Append(Form("_%04d%02d%02d_%02d%02d", 
			       year, month, day, hour, min));
    }

  }

  //__________________________________________________________________
  /** 
   * Parse a string into a type enum
   * 
   * @param type String to pass
   * 
   * @return Enumaration value 
   */
  static EType ParseType(const char* type, Bool_t& /*mc*/)
  {
    // mc = false;
    TString sType(type);
    sType.ToUpper();
    EType eType = kESD;
    // if      (sType.Contains("MC"))    mc    = true;
    if      (sType.Contains("ESD"))   eType = kESD; 
    else if (sType.Contains("AOD"))   eType = kAOD;
    else 
      Fatal("Run", "Unknown type '%s'", type);
    
    return eType;
  }
  //__________________________________________________________________
  /** 
   * Return a string that reflects the passed mode
   * 
   * @param eMode Mode 
   * 
   * @return String representation of mode 
   */
  static const char* ModeString(EMode eMode) 
  {
    switch (eMode) {
    case kLocal:	return "LOCAL";
    case kProof:	return "PROOF";
    case kGrid:		return "GRID";
    }
    return 0;
  }
  //__________________________________________________________________
  /** 
   * Parse a string for mode specifier 
   * 
   * @param mode Mode string
   * 
   * @return EMode value
   */
  static EMode ParseMode(const char* mode)
  {
    TString sMode(mode);
    sMode.ToUpper();
    EMode eMode = kLocal;
    if      (sMode == "LOCAL") eMode = kLocal;
    else if (sMode == "PROOF") eMode = kProof;
    else if (sMode == "GRID")  eMode = kGrid;
    else 
      Fatal("Run", "Unknown mode '%s'", mode);
    return eMode;
  }

  //__________________________________________________________________
  /** 
   * Return a string that reflects the passed operation
   * 
   * @param eOper Operation
   * 
   * @return String representation of operation 
   */
  static const char* OperString(EOper eOper) 
  {
    switch (eOper) {
    case kTest:		return "TEST";
    case kOffline:	return "OFFLINE";
    case kSubmit:	return "SUBMIT";
    case kTerminate:	return "TERMINATE";
    case kFull:		return "FULL";
    }
    return 0;
  }
  //__________________________________________________________________
  /** 
   * Parse an operation string 
   * 
   * @param oper Operation 
   * 
   * @return An EOper value
   */
  static EOper ParseOperation(const char* oper)
  {
    TString sOper(oper);
    sOper.ToUpper();
    EOper eOper = kFull;
    if      (sOper == "TEST")      eOper = kTest;
    else if (sOper == "OFFLINE")   eOper = kOffline;
    else if (sOper == "SUBMIT")    eOper = kSubmit;
    else if (sOper == "TERMINATE") eOper = kTerminate;
    else if (sOper == "FULL")      eOper = kFull;
    else 
      Fatal("Run", "unknown operation '%s'", oper);
    return eOper;
  }

  //__________________________________________________________________
  /** 
   * Set ROOT version to use 
   * 
   * @param v Version string of ROOT 
   */
  void SetROOTVersion(const char* v)    { fRootVersion = v; }
  //__________________________________________________________________
  /** 
   * Set AliROOT version to use 
   * 
   * @param v Version string of AliROOT 
   */
  void SetAliROOTVersion(const char* v) { fAliRootVersion = v; }
  //__________________________________________________________________
  /** 
   * Set the PROOF server URL
   * 
   * @param s PROOF server URL 
   */
  void SetProofServer(const char* s)    { fProofServer = s; }
  //__________________________________________________________________
  /** 
   * Set the GRID/Local data dir 
   * 
   * @param d Directory with data 
   */
  void SetDataDir(const char* d) { fDataDir = d; }
  //__________________________________________________________________
  /** 
   * Set the PROOF data set 
   * 
   * @param d PROOF registered data set 
   */
  void SetDataSet(const char* d) { fDataSet = d; }
  //__________________________________________________________________
  /** 
   * Set the XML file to use 
   * 
   * @param x XML file 
   */
  void SetXML(const char* x) { fXML = x; }
  //__________________________________________________________________
  /** 
   * Set how many replicas of the output we want 
   * 
   * @param n Number of replicas requested 
   */
  void SetNReplica(Int_t n) { fNReplica = n; }
  /** 
   * Set the ESD pass to use 
   * 
   * @param pass Pass number 
   */
  void SetESDPass(Int_t pass) { fESDPass = pass; }
  /** 
   * Set the ESD pass to use 
   * 
   * @param postfix Post fix to pass number 
   */
  void SetPassPostfix(const char* postfix) { fPassPostfix = postfix; }
  //__________________________________________________________________
  /** 
   * Use GDB to wrap PROOF slaves 
   * 
   * @param use Whether to use GDB or not 
   */
  void SetUseGDB(Bool_t use=kTRUE) { fUseGDB = use; }
  //__________________________________________________________________
  /** 
   * Add a source file to be copied and byte compiled on slaves 
   * 
   * @param src          Sources 
   * @param addToExtra   If false, do not copy 
   */
  void AddSource(const char* src, bool addToExtra=true) 
  { 
    fListOfSources.Add(new TObjString(src)); 
    if (addToExtra) AddExtraFile(src); // Source code isn't copied!
  }
  //__________________________________________________________________
  /** 
   * Add binary data to be uploaded to slaves 
   * 
   * @param lib Name of binary file 
   */
  void AddLibrary(const char* lib) { fListOfLibraries.Add(new TObjString(lib));}
  //__________________________________________________________________
  /** 
   * Add a run to be analysed
   *  
   * @param run Run number
   */
  void AddRun(Int_t run) 
  {
    Int_t i = fRunNumbers.fN; fRunNumbers.Set(i+1); fRunNumbers[i] = run;
  }
  //__________________________________________________________________
  /** 
   * Read run numbers from a file 
   * 
   * @param filename File name 
   */
  void ReadRunNumbers(const char* filename)
  {
    std::ifstream file(filename);
    if (!file) 
      Fatal("ReadRunNumbers", "Cannot read from %s", filename);
    
    while (!file.eof()) { 
      Int_t run;
      file >> run;
      AddRun(run);
      Char_t c;
      file >> c;
      if (file.bad()) break;
    }
    file.close();
  }
  //__________________________________________________________________
  /** 
   * Add an extra file to be uploaded to slave 
   * 
   * @param file Extra file to be uploaded 
   */
  void AddExtraFile(const char* file)
  {
    if (!file || file[0] == '\0') return;
    fListOfExtras.Add(new TObjString(file));
  }
  //__________________________________________________________________
  /** 
   * Set whether to allow overwritting existing files/directories 
   * 
   * @param allow If true, allow overwritting files/directories
   */
  void SetAllowOverwrite(Bool_t allow) { fAllowOverwrite = allow; }
  //__________________________________________________________________
  /** 
   * Print the setup 
   * 
   */
  void Print() const 
  {
    bool mc = AliAnalysisManager::GetAnalysisManager()
      ->GetMCtruthEventHandler();
    std::cout << fName << " train setup\n"
	      << std::boolalpha
	      << "  ROOT version:         " << fRootVersion    << "\n"
	      << "  AliROOT version:      " << fAliRootVersion << "\n"
	      << "  Name of proof server: " << fProofServer    << "\n"
	      << "  Grid Input directory: " << fDataDir        << "\n"
	      << "  Proof data set name:  " << fDataSet        << "\n"
	      << "  XML collection:       " << fXML            << "\n"
	      << "  Monte-Carlo input:    " << mc              << "\n"
	      << "  Storage replication:  " << fNReplica       << "\n"
	      << "  Run numbers:          " << std::flush;
    for (Int_t i = 0; i < fRunNumbers.GetSize(); i++) 
      std::cout << (i == 0 ? "" : ", ") << fRunNumbers.At(i);

    std::cout << "\n"
	      << "  PAR files:            " << std::flush;
    Bool_t first = true;
    TObject* obj = 0;
    TIter nextPar(&fListOfPARs);
    while ((obj = nextPar())) {
      std::cout << (first ? "" : ", ") << obj->GetName();
      first = false;
    }

    std::cout << "\n"
	      << "  Script sources:       " << std::flush;
    first = true;
    TIter nextSrc(&fListOfSources);
    while ((obj = nextSrc())) {
      std::cout << (first ? "" : ", ") << obj->GetName();
      first = false;
    }

    std::cout << "\n"
	      << "  Libraries to load:    " << std::flush;
    first = true;
    TIter nextLib(&fListOfLibraries);
    while ((obj = nextLib())) {
      std::cout << (first ? "" : ", ") << obj->GetName();
      first = false;
    }
    std::cout << std::noboolalpha << std::endl;

    AliAnalysisGrid* plugin = 
      AliAnalysisManager::GetAnalysisManager()->GetGridHandler();
    if (!plugin) return;
    
  }

protected:
  //__________________________________________________________________
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from
   */
  TrainSetup(const TrainSetup& o)
    : fName(o.fName),
      fRootVersion(o.fRootVersion),
      fAliRootVersion(o.fAliRootVersion),
      fProofServer(o.fProofServer),
      fDataDir(o.fDataDir),	
      fDataSet(o.fDataSet),	
      fXML(o.fXML),	
      fRunNumbers(o.fRunNumbers),
      fListOfPARs(),
      fListOfSources(),
      fListOfLibraries(),
      fListOfExtras(),
      fNReplica(o.fNReplica),
      fESDPass(o.fESDPass)
  {
    if (isdigit(fName[0])) { 
      Warning("TrainSetup", "Name starts with a digit, prepending 'a' to name");
      fName = Form("a%s", fName.Data());
    }
    TObject* obj = 0;
    TIter nextPar(&o.fListOfPARs);
    while ((obj = nextPar())) fListOfPARs.Add(obj->Clone());
    TIter nextSrc(&o.fListOfSources);
    while ((obj = nextSrc())) fListOfSources.Add(obj->Clone());
    TIter nextLib(&o.fListOfLibraries);
    while ((obj = nextLib())) fListOfLibraries.Add(obj->Clone());
    TIter nextExa(&o.fListOfExtras);
    while ((obj = nextExa())) fListOfExtras.Add(obj->Clone());
  }
  //__________________________________________________________________
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object. 
   */
  TrainSetup& operator=(const TrainSetup& o)
  {
    fName		= o.fName;
    fRootVersion	= o.fRootVersion;
    fAliRootVersion	= o.fAliRootVersion;
    fProofServer	= o.fProofServer;
    fDataDir		= o.fDataDir;	
    fDataSet		= o.fDataSet;	
    fXML		= o.fXML;	
    fNReplica		= o.fNReplica;	
    fESDPass            = o.fESDPass;
    fRunNumbers         = o.fRunNumbers;
    TObject* obj = 0;
    TIter nextPar(&o.fListOfPARs);
    while ((obj = nextPar())) fListOfPARs.Add(obj->Clone());
    TIter nextSrc(&o.fListOfSources);
    while ((obj = nextSrc())) fListOfSources.Add(obj->Clone());
    TIter nextLib(&o.fListOfLibraries);
    while ((obj = nextLib())) fListOfLibraries.Add(obj->Clone());
    TIter nextExa(&o.fListOfExtras);
    while ((obj = nextExa())) fListOfExtras.Add(obj->Clone());

    return *this;
  }

  //__________________________________________________________________
  /** 
   * Run this analysis 
   * 
   * @param type    Type of input for analysis  (kESD, kAOD)
   * @param mode    Mode of job (kLocal, kProof, kGrid)
   * @param oper    Operation 
   * @param nEvents Number of events to analyse (<0 means all)
   * @param mc      Whether to connect MC data 
   * @param usePar  Whether to use PARs  
   * @param dbg     Debug level
   */
  void Exec(const char*  type, 
	    const char*  mode="GRID", 
	    const char*  oper="FULL", 
	    Int_t        nEvents=-1, 
	    Bool_t       mc=false, 
	    Bool_t       usePar=false, 
	    Int_t        dbg=0)
  {
    Info("Exec", "Doing exec with type=%s, mode=%s, oper=%s, events=%d "
	 "mc=%d, usePar=%d", type, mode, oper, nEvents, mc, usePar);
    EType eType = ParseType(type, mc);
    EMode eMode = ParseMode(mode);
    EOper eOper = ParseOperation(oper);

    Exec(eType, eMode, eOper, nEvents, mc, usePar, dbg);
  }

  //__________________________________________________________________
  /** 
   * Run this analysis 
   * 
   * @param type    Type of input for analysis  (kESD, kAOD)
   * @param mode    Mode of job (kLocal, kProof, kGrid)
   * @param oper    Operation 
   * @param nEvents Number of events to analyse (<0 means all)
   * @param mc      Whether to connect MC data 
   * @param usePar  Whether to use PARs  
   * @param dbg     Debug level
   */
  void Exec(EType  type, 
	    EMode  mode, 
	    EOper  oper, 
	    Int_t  nEvents, 
	    Bool_t mc, 
	    Bool_t usePar, 
	    Int_t  dbg=0)
  {
    Info("Exec", "Doing exec with type=%d, mode=%d, oper=%d, events=%d "
	 "mc=%d, usePar=%d", type, mode, oper, nEvents, mc, usePar);

    if (mode == kProof) usePar    = true;

    if (!Connect(mode)) return;

    TString cwd = gSystem->WorkingDirectory();
    TString nam = EscapedName();
    if (oper != kTerminate) { 
      if (!fAllowOverwrite && !gSystem->AccessPathName(nam.Data())) {
	Error("Exec", "File/directory %s already exists", nam.Data());
	return;
      }
      if (gSystem->AccessPathName(nam.Data())) {
	if (gSystem->MakeDirectory(nam.Data())) {
	  Error("Exec", "Failed to make directory %s", nam.Data());
	  return;
	}
      }
    }
    else {
      if (gSystem->AccessPathName(nam.Data())) {
	Error("Exec", "File/directory %s does not exists", nam.Data());
	return;
      }
    }
      
    if (!gSystem->ChangeDirectory(nam.Data())) { 
      Error("Exec", "Failed to change directory to %s", nam.Data());
      return;
    }
    Info("Exec", "Made subdirectory %s, and cd'ed there", nam.Data());
      
    if (!LoadCommonLibraries(mode, usePar)) return;
    
    // --- Create analysis manager -----------------------------------
    AliAnalysisManager *mgr  = new AliAnalysisManager(fName,"Analysis Train");

    // In test mode, collect system information on every event 
    // if (oper == kTest)  mgr->SetNSysInfo(1); 
    if (dbg  >  0)      mgr->SetDebugLevel(dbg);
    if (mode == kLocal) mgr->SetUseProgressBar(kTRUE, 100);
   
    // --- ESD input handler ------------------------------------------
    AliVEventHandler*  inputHandler = CreateInputHandler(type);
    if (inputHandler) mgr->SetInputEventHandler(inputHandler);
    
    // --- Monte-Carlo ------------------------------------------------
    AliVEventHandler*  mcHandler = CreateMCHandler(type,mc);
    if (mcHandler) mgr->SetMCtruthEventHandler(mcHandler);
    
    // --- AOD output handler -----------------------------------------
    AliVEventHandler*  outputHandler = CreateOutputHandler(type);
    if (outputHandler) mgr->SetOutputEventHandler(outputHandler);
    
    // --- Include analysis macro path in search path ----------------
    gROOT->SetMacroPath(Form("%s:%s:$ALICE_ROOT/ANALYSIS/macros",
			     cwd.Data(), gROOT->GetMacroPath()));

    // --- Physics selction - only for ESD ---------------------------
    if (type == kESD) CreatePhysicsSelection(mc, mgr);
    
    // --- Create centrality task ------------------------------------
    CreateCentralitySelection(mc, mgr);

    // --- Create tasks ----------------------------------------------
    CreateTasks(mode, usePar, mgr);

    // --- Create Grid handler ----------------------------------------
    // _must_ be done after all tasks has been added
    AliAnalysisAlien* gridHandler = CreateGridHandler(type, mode, oper);
    if (gridHandler) mgr->SetGridHandler(gridHandler);
    
    // --- Create the chain ------------------------------------------
    TChain* chain = CreateChain(type, mode, oper, mc);
    if (mode == kLocal && !chain) {
      Error("Exec", "No chain defined in local mode!");
      return;
    }

    // --- Print setup -----------------------------------------------
    Print();

    // --- Initialise the train --------------------------------------
    if (!mgr->InitAnalysis())  {
      gSystem->ChangeDirectory(cwd.Data());
      Error("Run","Failed to initialise train");
      return;
    }

    // --- Show status -----------------------------------------------
    mgr->PrintStatus();

    Long64_t ret = StartAnalysis(mgr, mode, chain, nEvents);

    // Make sure we go back 
    gSystem->ChangeDirectory(cwd.Data());

    if (ret < 0) Error("Exec", "Analysis failed");
  }
  //__________________________________________________________________
  /** 
   * Start the analysis 
   * 
   * @param mgr       Analysis manager
   * @param mode      Run mode
   * @param chain     Input data (local and proof only)
   * @param nEvents   Number of events to analyse 
   */
  Long64_t StartAnalysis(AliAnalysisManager* mgr, 
			 EMode               mode, 
			 TChain*             chain,
			 Int_t               nEvents)
  {
    // --- Run the analysis ------------------------------------------
    switch (mode) { 
    case kLocal: 
      if (!chain) {
	Error("StartAnalysis", "No chain defined");
	return -1;
      }
      if (nEvents < 0) nEvents = chain->GetEntries();
      return mgr->StartAnalysis(ModeString(mode), chain, nEvents);
    case kProof: 
      if (fDataSet.IsNull()) {
	if (!chain) { 
	  Error("StartAnalysis", "No chain defined");
	  return -1;
	}
	if (nEvents < 0) nEvents = chain->GetEntries();
	return mgr->StartAnalysis(ModeString(mode), chain, nEvents);
      }
      return mgr->StartAnalysis(ModeString(mode), fDataSet);
    case kGrid: 
      if (nEvents < 0)
	return mgr->StartAnalysis(ModeString(mode));
      return mgr->StartAnalysis(ModeString(mode), nEvents);
    }
    // We should never get  here 
    return -1;
  }
  //__________________________________________________________________
  /** 
   * Return the escaped name 
   * 
   * 
   * @return Escaped name 
   */
  const TString& EscapedName() const 
  {
    return fEscapedName;
  }
  //__________________________________________________________________
  /** 
   * Create a grid handler 
   * 
   * @param type Data type
   * @param mode Run mode 
   * @param oper Operation 
   * 
   * @return Grid handler 
   */
  virtual AliAnalysisAlien* 
  CreateGridHandler(EType type, EMode mode, EOper oper)
  {
    if (mode != kGrid) return 0;

    TString name = EscapedName();

    // Create the plug-in object, and set run mode 
    AliAnalysisAlien* plugin = new AliAnalysisAlien();
    plugin->SetRunMode(OperString(oper));
    
    // Production mode - not used here 
    // plugin->SetProductionMode();
    
    // Set output to be per run 
    plugin->SetOutputToRunNo();

    // Set the job tag 
    plugin->SetJobTag(fName);

    // Set number of test files - used in test mode only 
    plugin->SetNtestFiles(1);
    
    // Set required version of software 
    plugin->SetAPIVersion(fAliEnAPIVersion);
    plugin->SetROOTVersion(fRootVersion);
    plugin->SetAliROOTVersion(fAliRootVersion);

    // Keep log files 
    plugin->SetKeepLogs();

    // Declare root of input data directory 
    plugin->SetGridDataDir(fDataDir);

    // Data search patterns 
    TString pat;
    if (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()) {
      pat = "*/";
      plugin->SetRunPrefix("");
    }
    else {
      pat = Form("*ESDs/pass%d%s/*/", fESDPass, fPassPostfix.Data());
      plugin->SetRunPrefix("000");
    }
    pat.Append(Form("*%s.root", type == kESD ? "ESDs" : "AOD"));
    plugin->SetDataPattern(pat);

    // Add the run numbers 
    for (Int_t i = 0; i < fRunNumbers.fN; i++) {
      if (fRunNumbers[i] < 0) continue; 
      plugin->AddRunNumber(fRunNumbers[i]);
    }
    
    // Set the working directory to be the trains name (with special
    // characters replaced by '_' and the date appended), and also set
    // the output directory (relative to working directory)
    plugin->SetGridWorkingDir(name.Data());
    plugin->SetGridOutputDir("output");

    // Enable configured PARs 
    TIter nextPar(&fListOfPARs);
    TObject* parName;
    while ((parName = nextPar()))
      plugin->EnablePackage(parName->GetName());
    
    // Add sources that need to be compiled on the workers using
    // AcLIC. 
    TString addSources = SetupSources();
    if (!addSources.IsNull()) plugin->SetAnalysisSource(addSources.Data());

    // Add binary libraries that should be uploaded to the workers 
    TString addLibs = SetupLibraries();
    if (!addLibs.IsNull()) plugin->SetAdditionalLibs(addLibs.Data());
    
    // Disable default outputs 
    plugin->SetDefaultOutputs(true);

    // Merge parameters 
    plugin->SetMaxMergeFiles(20);
    plugin->SetMergeExcludes("AliAOD.root "
			    "*EventStat*.root "
			    "*event_stat*.root");

    // Set number of runs per master - set to one to per run
    plugin->SetNrunsPerMaster(1);

    // Loop over defined containers in the analysis manager, 
    // and declare these as outputs 
    TString listOfAODs  = "";
    TString listOfHists = "";
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    AliAnalysisDataContainer* cont = 0;
    TIter nextCont(mgr->GetOutputs());
    while ((cont = static_cast<AliAnalysisDataContainer*>(nextCont()))) {
      TString outName(cont->GetFileName());
      TString& list = (outName == "default" ? listOfAODs : listOfHists);
      if (outName == "default") { 
	if (!mgr->GetOutputEventHandler()) continue; 

	outName = mgr->GetOutputEventHandler()->GetOutputFileName();
      }
      if (list.Contains(outName)) continue;
      if (!list.IsNull()) list.Append(",");
      list.Append(outName);
    }
    if (!mgr->GetExtraFiles().IsNull()) { 
      if (!listOfAODs.IsNull()) listOfAODs.Append("+");
      TString extra = mgr->GetExtraFiles();
      extra.ReplaceAll(" ", ",");
      listOfAODs.Append(extra);
   }
    TString outArchive = Form("stderr, stdout@disk=%d", fNReplica);
    if (!listOfHists.IsNull()) 
      outArchive.Append(Form(" hist_archive.zip:%s@disk=%d", 
			     listOfHists.Data(), fNReplica));
    if (!listOfAODs.IsNull()) 
      outArchive.Append(Form(" aod_archive.zip:%s@disk=%d", 
			     listOfAODs.Data(), fNReplica));
    if (listOfAODs.IsNull() && listOfHists.IsNull()) 
      Fatal("CreateGridHandler", "No outputs defined");
    // Disabled for now 
    // plugin->SetOutputArchive(outArchive);
    
    // Set name of generated analysis macro 
    plugin->SetAnalysisMacro(Form("%s.C", name.Data()));
    
    // Maximum number of sub-jobs 
    // plugin->SetSplitMaxInputFileNumber(25);
    
    // Set the Time-To-Live 
    plugin->SetTTL(70000);
    
    // Re-submit failed jobs as long as the ratio of failed jobs is
    // below this percentage. 
    plugin->SetMasterResubmitThreshold(95);

    // Set the input format
    plugin->SetInputFormat("xml-single");

    // Set the name of the generated jdl 
    plugin->SetJDLName(Form("%s.jdl", name.Data()));

    // Set the name of the generated executable 
    plugin->SetExecutable(Form("%s.sh", name.Data()));
    
    // Set the job price !?
    plugin->SetPrice(1);

    // Set whether to merge via JDL 
    plugin->SetMergeViaJDL(true);
    
    // Fast read otion 
    plugin->SetFastReadOption(false);

    // Whether to overwrite existing output 
    plugin->SetOverwriteMode(true);

    // Set the executable binary name and options 
    plugin->SetExecutableCommand("aliroot -b -q -x");
    
    // Split by storage element - must be lower case!
    plugin->SetSplitMode("se");

    return plugin;
  }
  //__________________________________________________________________
  /** 
   * Create input handler 
   * 
   * @param type 
   * 
   * @return 
   */
  virtual AliVEventHandler* CreateInputHandler(EType type)
  {
    switch (type) {
    case kESD: return new AliESDInputHandler(); 
    case kAOD: return new AliAODInputHandler(); 
    }
    return 0;
  }
  //__________________________________________________________________
  /** 
   * Create input handler 
   * 
   * @param type  Run type (ESD or AOD)
   * @param mc    Assume monte-carlo input 
   * 
   * @return 
   */
  virtual AliVEventHandler* CreateMCHandler(EType type, bool mc)
  {
    if (!mc)          return 0;
    if (type != kESD) return 0;
    Info("CreateMCHandler", "Making MC handler");
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mcHandler->SetReadTR(true); 
    return mcHandler;
  }
  //__________________________________________________________________
  /** 
   * Create output event handler 
   * 
   * @param type 
   * 
   * @return 
   */
  virtual AliVEventHandler* CreateOutputHandler(EType type)
  {
    AliAODHandler* ret = new AliAODHandler();
    switch (type) { 
    case kESD: 
      ret->SetOutputFileName("AliAOD.root");
      break;
    case kAOD: 
      ret->SetOutputFileName("AliAOD.pass2.root");
      break;
    }
    return ret;
  }
  //__________________________________________________________________
  /** 
   * Create physics selection , and add to manager
   * 
   * @param mc Whether this is for MC 
   * @param mgr Manager
   */
  virtual void CreatePhysicsSelection(Bool_t mc,
				      AliAnalysisManager* mgr)
  {
    gROOT->Macro(Form("AddTaskPhysicsSelection.C(%d)", mc));
    mgr->RegisterExtraFile("event_stat.root");
  }
  //__________________________________________________________________
  /** 
   * Create physics selection , and add to manager
   * 
   * @param mc Whether this is for MC 
   * @param mgr Manager
   */
  virtual void CreateCentralitySelection(Bool_t mc, AliAnalysisManager* mgr)
  {
    gROOT->Macro("AddTaskCentrality.C");
    AliCentralitySelectionTask* ctask = 
      dynamic_cast<AliCentralitySelectionTask*>(mgr->GetTask("CentralitySelection"));
    if (!ctask) return;
    ctask->SetPass(fESDPass);
    if (mc) ctask->SetMCInput();
  }
  //__________________________________________________________________
  /** 
   * Create analysis tasks 
   * 
   * @param mode Run mode
   * @param mgr  Manager
   * @param par  Whether to use pars 
   */
  virtual void CreateTasks(EMode mode, Bool_t par, AliAnalysisManager* mgr)=0;
  //__________________________________________________________________
  /** 
   * Connect to external services (Proof and/or grid)
   * 
   * @param mode Running mode 
   * 
   * @return true on success 
   */
  virtual Bool_t Connect(EMode mode)
  {
    if (mode == kLocal) return true;
			  
    // --- Set-up connections to Proof cluster and alien -------------
    if (mode == kProof) { 
      // --- Find user name ------------------------------------------
      TString userName(gSystem->Getenv("alien_API_USER"));
      if (userName.IsNull()) {
	userName = gSystem->GetUserInfo()->fUser;
	Warning("Connect", 
		"environment variable 'alien_API_USER' not set, using %s", 
		userName.Data());
      }

      // --- Set prefered GSI method ---------------------------------
      gEnv->SetValue("XSec.GSI.DelegProxy", "2");
      
      // --- Now open connection to PROOF cluster --------------------
      TString serv = "";
      Bool_t  lite = false;
      if (fProofServer.BeginsWith("workers=") || fProofServer.IsNull()) {
	lite = true;
	serv = fProofServer;
      }
      else 
	serv = Form("%s@%s", userName.Data(), fProofServer.Data());
      TProof::Open(serv);
      if (!gProof) { 
	Error("Connect", "Failed to connect to Proof cluster %s as %s",
	      fProofServer.Data(), userName.Data());
	return false;
      }
      if (fUseGDB) { 
	TProof::AddEnvVar("PROOF_WRAPPERCMD", 
			  "gdb --batch -ex run -ex bt --args ");
      }
      if (lite) return true;
    }

    // --- Open a connection to the grid -----------------------------
    TGrid::Connect("alien://");
    if (!gGrid || !gGrid->IsConnected()) { 
      // This is only fatal in grid mode 
      Error("Connect", "Failed to connect to AliEN");
      if (mode == kGrid) return false; 
      return true;
    }
    if (mode == kGrid) return true;

    
    // --- Set and make output directory -----------------------------
    TString name = EscapedName();
    TString homeDir(gGrid->GetHomeDirectory());
    TString workDir(homeDir);
    workDir.Append("/");
    workDir.Append(name);
    
    // Make working directory 
    if (!gGrid->Cd(workDir)) { 
      gGrid->Cd(homeDir);
      if (gGrid->Mkdir(workDir)) {
	gGrid->Cd(name);
	Info("Connect", "Directory %s created", workDir.Data());
      }
    }
    // Make output directory 
    gGrid->Mkdir("proof_output");
    gGrid->Cd("proof_output");

    return true;
  }	  
  //__________________________________________________________________
  /** 
   * Load common libraries 
   * 
   * @param mode Running mode			
   * @param par  If true, load as PARs 
   * 
   * @return true on success 
   */
  Bool_t LoadCommonLibraries(EMode mode, Bool_t par) 
  {
    if (!gSystem->Getenv("ALICE_ROOT")) { 
      Error("LoadCommonLibraries", "Local AliROOT not available");
      return false;
    }
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libMinuit.so");

    Bool_t ret   = true;
    Bool_t basic = mode == kGrid ? false : par;
    
    ret &= LoadLibrary("STEERBase",     mode, basic, false);
    ret &= LoadLibrary("ESD",           mode, basic, false);
    ret &= LoadLibrary("AOD",           mode, basic, false);
    ret &= LoadLibrary("ANALYSIS",      mode, basic, true);
    ret &= LoadLibrary("ANALYSISalice", mode, basic, true);

    return ret;
  }
  //__________________________________________________________________
  /** 
   * Load a library 
   * 
   * @param what What library to load
   * @param mode Mode (local, proof, grid)
   * @param par  If true, load as PAR
   * @param rec  If true, also load on slaves
   * 
   * @return true on success 
   */
  Bool_t LoadLibrary(const char* what, EMode mode, Bool_t par, Bool_t rec=false)
  {
    if (!what || what[0] == '\0') return true;
    
    TString module(what);
    TString libName(what);
    if (!libName.BeginsWith("lib")) libName = Form("lib%s", libName.Data());
    if (!libName.EndsWith(".so"))   libName.Append(".so");

    Int_t ret = 0;

    switch (mode) { 
    case kLocal: // Just load and exit 
      gSystem->Load(libName.Data());
      break;
    case kGrid: 
      if (par) { 
	ret = SetupPAR(what) ? 0 : -1;
	if (rec) fListOfPARs.Add(new TObjString(what));
      } else  {
	ret = gSystem->Load(libName.Data());
	if (rec) fListOfLibraries.Add(new TObjString(libName));
      }
      break;
    case kProof: 
      ret = gProof->UploadPackage(what);
      if (ret < 0)  {	
	  ret = gProof->UploadPackage(gSystem->ExpandPathName(Form("../%s.par",
								   what)));
	if (ret < 0) {	
	  ret = 
	    gProof->UploadPackage(gSystem
				  ->ExpandPathName(Form("$ALICE_ROOT/%s.par", 
							what)));
	  if (ret < 0) {
	    Error("LoadLibrary", 
		  "Could not find module %s.par in current directory nor "
		  "in $ALICE_ROOT", module.Data());
	    return false;
	  }
	}
      }
      ret = gProof->EnablePackage(what);
      break;
    }
    if (ret < 0) { 
      Error("LoadLibrary", "Couldn't load %s", what);
      return false;
    }
    return true;
  }
          
  //__________________________________________________________________
  Bool_t SetupPAR(const char* what)
  {
    if (!what || what[0] == '\0') return -1;
    
    TString parFile(Form("%s.par", what));
    if (gSystem->AccessPathName(parFile.Data())) { 
      if (gSystem->AccessPathName(Form("../%s.par", what))) { 
	// If not found 
	TString aliParFile = 
	  gSystem->ExpandPathName(Form("$(ALICE_ROOT)/%s.par", what));
	if (gSystem->AccessPathName(aliParFile.Data())) { 
	  Error("SetupPAR", "PAR file %s not found in current directory or "
		"$(ALICE_ROOT)", what);
	  return false;
	}
	// Copy to current directory 
	TFile::Cp(aliParFile, parFile);
      }
      else 
	gSystem->Exec(Form("ln -s ../%s.par .", what));
    }
    
    // Extract archive 
    gSystem->Exec(Form("tar xvzf %s", parFile.Data()));
    
    // Change directory into par archive
    TString cwd = gSystem->WorkingDirectory();
    
    if (!gSystem->ChangeDirectory(what)) { 
      Error("SetupPAR", "Failed to change directory to %s", what);
      return false;
    }
    
    // Test the build 
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      Info("SetupPar", "Building in PAR archive %s", what);
      if (gSystem->Exec("PROOF-INF/BUILD.sh")) { 
	Error("SetupPar", "Failed to build in PAR directory %s", what);
	gSystem->ChangeDirectory(cwd.Data());
	return false;
      }
    }
    
    // Check for setup script
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      Info("SetupPAR", "Setting up for PAR %s", what);
      gROOT->Macro("PROOF-INF/SETUP.C");
    }
    if (!gSystem->ChangeDirectory(cwd.Data())) return false;

    return true;
  }
  //__________________________________________________________________
  TString SetupExtras()
  {
    TString ret;
    TIter next(&fListOfExtras);
    TObjString* obj = 0;
    while ((obj = static_cast<TObjString*>(next()))) {
      TString path = gSystem->ExpandPathName(obj->GetName());
      if (!path.BeginsWith("/")) 
	// If not an absolute path, prepend to up-one
	path = Form("../%s", path.Data());
      if (gSystem->AccessPathName(path.Data())) { 
	// File not accessible 
	Warning("SetupExtras", "File %s not accessible", path.Data());
	continue;
      }
      ret.Append(Form("%s ", gSystem->BaseName(obj->GetName())));
      gSystem->Exec(Form("ln -s %s .", path.Data()));
    }
    ret = ret.Strip();
    return ret;
  }
  //__________________________________________________________________
  TString SetupSources()
  {
    TString nam = EscapedName();
    TString ret;
    TIter next(&fListOfSources); 
    TObject* src;
    while ((src = next())) {
      TString path = gSystem->ExpandPathName(src->GetName());
      if (!path.BeginsWith("/")) 
	// If not an absolute path, prepend to up-one
	path = Form("../%s", path.Data());
      if (gSystem->AccessPathName(path.Data())) { 
	// File not accessible 
	Warning("SetupSources", "File %s not accessible", path.Data());
	continue;
      }
      ret.Append(Form("%s ", gSystem->BaseName(src->GetName())));
      gSystem->Exec(Form("ln -s %s .", path.Data()));
    }
    ret = ret.Strip();
    return ret;
  }
  //__________________________________________________________________
  TString SetupLibraries()
  {
    TString ret;
    TIter next(&fListOfLibraries); 
    TObject* lib;
    while ((lib = next())) {
      ret.Append(lib->GetName());
      ret.Append(" ");
    }
    // Also add extra files to this variable 
    ret.Append(SetupExtras());
    ret = ret.Strip();
    return ret;
  }
  //__________________________________________________________________
  /** 
   * Scan directory @a dir (possibly recursive) for tree files to add
   * to the chain.    This does not follow sym-links
   * 
   * @param dir        Directory to scan
   * @param chain      Chain to add to
   * @param type       Type of tree (ESD or AOD)
   * @param recursive  Whether to scan recursively 
   * @param mc         Look also for MC files if true 
   *
   * @return true if any files where added 
   */
  Bool_t ScanDirectory(TSystemDirectory* dir, TChain* chain, 
		       EType type, bool recursive, bool mc)
  {
    TString fnPattern;
    switch (type) { 
    case kESD:  fnPattern = "AliESD"; break;
    case kAOD:  fnPattern = "AliAOD"; break;
    }

    // Assume failure 
    Bool_t ret = false;

    // Get list of files, and go back to old working directory
    TString oldDir(gSystem->WorkingDirectory());
    TList*  files = dir->GetListOfFiles();
    if (!gSystem->ChangeDirectory(oldDir)) { 
      Error("ScanDirectory", "Failed to go back to %s", oldDir.Data());
      return false;
    }
    if (!files) return false;

    TList toAdd;
    toAdd.SetOwner();
    Bool_t hasGAlice = (!mc ? true : false);
    Bool_t hasKine   = (!mc ? true : false);
    Bool_t hasTrRef  = (!mc ? true : false);
    
    // Sort list of files and check if we should add it 
    files->Sort();
    TIter next(files);
    TSystemFile* file = 0;
    while ((file = static_cast<TSystemFile*>(next()))) {
      TString name(file->GetName());
      TString title(file->GetTitle());
      TString full(gSystem->ConcatFileName(file->GetTitle(), name.Data()));
      if (dynamic_cast<TSystemDirectory*>(file)) full = title;
      // Ignore special links 
      if (name == "." || name == "..") { 
	// Info("ScanDirectory", "Ignoring %s", name.Data());
	continue;
      }

      FileStat_t fs;
      if (gSystem->GetPathInfo(full.Data(), fs)) {
	Warning("ScanDirectory", "Cannot stat %s (%s)", full.Data(),
                gSystem->WorkingDirectory());
	continue;
      }
      // Check if this is a directory 
      if (file->IsDirectory(full)) { 
	if (recursive) {
	  // if (title[0] == '/') 
	  TSystemDirectory* d = new TSystemDirectory(file->GetName(),
						     full.Data());
          if (ScanDirectory(d,chain,type,recursive,mc))
	    ret = true;
	  delete d;
	}
        continue;
      }
    
      // If this is not a root file, ignore 
      if (!name.EndsWith(".root")) continue;

      // If this file does not contain AliESDs, ignore 
      if (!name.Contains(fnPattern)) { 
	// Info("ScanDirectory", "%s does not match pattern %s", 
	//      name.Data(), fnPattern.Data());
	if (mc) { 
	  if (name.CompareTo("galice.root") == 0)     hasGAlice = true;
	  if (name.CompareTo("Kinematics.root") == 0) hasKine   = true;
	  if (name.CompareTo("TrackRefs.root")  == 0) hasTrRef = true;
	}
	continue;
      }
    
      // Add 
      // Info("ScanDirectory", "Adding %s", full.Data());
      toAdd.Add(new TObjString(full));
    }

    if (mc && toAdd.GetEntries() > 0 && 
	(!hasGAlice || !hasKine || !hasTrRef)) { 
      Warning("ScanDirectory", 
	      "one or more of {galice,Kinematics,TrackRefs}.root missing from "
	      "%s, not adding anything from this directory", 
	      dir->GetTitle());
      toAdd.Delete();
    }

    TIter nextAdd(&toAdd);
    TObjString* s = 0;
    while ((s = static_cast<TObjString*>(nextAdd()))) {
      // Info("ScanDirectory", "Adding %s", s->GetString().Data());
      chain->Add(s->GetString());
    }
    if (toAdd.GetEntries() > 0) ret = true;

    gSystem->ChangeDirectory(oldDir);
    return ret;
  }
  //__________________________________________________________________
  /** 
   * Create a chain from an XML containing an collection
   * 
   * @param treeName Name of tree's 
   * @param xmlFile  XML collection
   * 
   * @return Newly allocated chain or null
   */
  TChain* CreateChainFromXML(const char* treeName, 
			     const char* xmlFile) 
  {
    TGridCollection* collection = TAlienCollection::Open(xmlFile);
    if (!collection) { 
      Error("CreateChainFromXML", "Cannot create AliEn collection from "
	    "XML file %s", xmlFile);
      return 0;
    }

    TChain* chain = new TChain(treeName);
    collection->Reset();
    while (collection->Next()) chain->Add(collection->GetTURL(""));
    
    return chain;
  }
  //__________________________________________________________________
  /** 
   * Create a chain of data 
   * 
   * @param type Type of data
   * @param mode Operation mode 
   * @param mc   Assume MC input if true
   *
   * @return TChain of data 
   */    
  TChain* CreateChain(EType type, EMode mode, EOper /* oper */, Bool_t mc)
  {
    TString treeName;
    switch (type) { 
    case kESD:  treeName = "esdTree"; break;
    case kAOD:  treeName = "aodTree"; break;
    }

    TChain* chain = 0;
    switch (mode) { 
    case kProof: 
      if (!fDataSet.IsNull()) break; 
      // Otherwise fall through
    case kLocal:
      if (fXML.IsNull()) {
	chain = new TChain(treeName.Data());
	TString dir(fDataDir);
	if (dir == ".") dir = "";
	if (!dir.BeginsWith("/")) dir = Form("../%s", dir.Data());
	TString savdir(gSystem->WorkingDirectory());
	TSystemDirectory d(gSystem->BaseName(dir.Data()), dir.Data());
	if (!ScanDirectory(&d, chain, type, true, mc)) { 
	  delete chain;
	  chain = 0;
	}
	gSystem->ChangeDirectory(savdir);
      }
      else 
	chain = CreateChainFromXML(treeName.Data(), fXML.Data());
      break;
    case kGrid:  break; // Do nothing - we use plugin
    }
    
    if (chain && chain->GetNtrees() <= 0) { 
      delete chain;
      return 0;
    }
    return chain;
  }
  //__________________________________________________________________
  TString fName;             // Name of analysis
  TString fRootVersion;      // ROOT version to use 
  TString fAliRootVersion;   // AliROOT version to use 
  TString fAliEnAPIVersion;  // AliEn API version to use 
  TString fProofServer;      // Name of proof server
  TString fDataDir;          // Grid Input directory 
  TString fDataSet;          // Proof data set name 
  TString fXML;              // XML collection for local/proof mode
  TArrayI fRunNumbers;       // List of run number 
  TList   fListOfPARs;       // List of PAR files to use 
  TList   fListOfSources;    // List of sources to upload and AcLIC
  TList   fListOfLibraries;  // List of libraries to load
  TList   fListOfExtras;     // List of extra files to upload
  Int_t   fNReplica;         // Storage replication
  Int_t   fESDPass;          // ESD pass number 
  TString fPassPostfix;      // Possible pass postfix
  TString fEscapedName;      // Name escaped for special chars
  Bool_t  fAllowOverwrite;   // Allow overwriting output dir
  Bool_t  fUseGDB;           // Wrap PROOF slaves in GDB 
};


void
BuildTrainSetup()
{
  gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWG2/FORWARD/analysis2:"
			   "$ALICE_ROOT/ANALYSIS/macros",
			   gROOT->GetMacroPath()));
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include");
  gSystem->Load("libRAliEn");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  TString path = gSystem->Which(gROOT->GetMacroPath(), "TrainSetup.C");
  Info("BuildTrainSetup", "Path=%s", path.Data());
  TString tmp("TrainSetup");
  FILE* fp = gSystem->TempFileName(tmp, ".");
  fclose(fp);
  gSystem->Unlink(tmp);
  tmp.Append(".C");
  Info("BuildTrainSetup", "Copy %s -> %s", path.Data(), tmp.Data());
  gSystem->CopyFile(path, tmp);
  gROOT->LoadMacro(Form("%s+g", tmp.Data()));
  gSystem->Unlink(tmp);
  tmp.ReplaceAll(".C", "_C.so");
  gSystem->Unlink(tmp);
  tmp.ReplaceAll("_C.so", "_C.d");
  gSystem->Unlink(tmp);
}

  
//____________________________________________________________________
//
// EOF
//
