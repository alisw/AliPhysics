/**
 * @defgroup pwglf_forward_trains Trains
 * 
 * Train specifications 
 *
 * @ingroup pwglf_forward
 */
/**
 * @file   TrainSetup.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Wed Mar 23 12:12:00 2011
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_trains
 * 
 */

#ifndef __CINT__
#include <fstream>
#include <iostream>
#include <cstdlib>

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
class TDatime;
class TString;
class TSystemDirectory;
#endif



//====================================================================
/** 
 * Generic set-up of an analysis train using the grid-handler (AliEn plugin). 
 *
 * See also @ref train_setup_doc
 *
 * @ingroup pwglf_forward_trains
 * 
 */
struct TrainSetup
{
  // Forward declaration 
  class Runner;
  /** 
   * Data type to process 
   */
  enum EType { 
    /** Event Summary Data */
    kESD, 
    /** Analysis Object Data */
    kAOD,
    /** User defined */
    kUser
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
    kFull,
    /** Only intialize */
    kInitialize
  };


  //__________________________________________________________________
  virtual ~TrainSetup() {}

  //__________________________________________________________________
  /** 
   * Constructor 
   * 
   * @param name         Name of analysis (free-form)
   */
  TrainSetup(const char* name)
    : fName(name),
      fEscapedName(name),
      fRootVersion("v5-28-00a"),
      fAliRootVersion("v4-21-18-AN"),
      fAliEnAPIVersion("V1.1x"),
      fProofServer("alicecaf.cern.ch"),
      fDataDir("/alice/data/2010/LHC10c"),
      fDataPattern("*"),
      fDataSet("/COMMON/COMMON/LHC09a4_run8100X#/esdTree"),
      fXML(""), 
      fNReplica(4),
      fAllowOverwrite(kFALSE),
      fUseGDB(kFALSE), 
      fMaxSplit(50),
      fRunNumbers(0),
      fListOfPARs(),
      fListOfSources(),
      fListOfLibraries(),
      fListOfExtras(),
      fDatime(1995, 0, 0, 0, 0, 0), 
      fExecType(kUser), 
      fExecMode(kLocal), 
      fExecOper(kFull),
      fUsePar(false), 
      fMC(false), 
      fPerRunMerge(false),
      fVerbose(0)
  {
    fEscapedName = EscapeName(fName, fDatime);
  }
  //__________________________________________________________________
  /** 
   * @{ 
   * @name Software environment 
   */
  /** 
   * Set ROOT version to use 
   * 
   * @param v Version string of ROOT 
   */
  void SetROOTVersion(const char* v)    { fRootVersion = v; }
  /** 
   * Set AliROOT version to use 
   * 
   * @param v Version string of AliROOT 
   */
  void SetAliROOTVersion(const char* v) { fAliRootVersion = v; }
  /** 
   * Set the AliEn API version to use 
   * 
   * @param v AliEn API version 
   */
  void SetAliEnAPIVersion(const char* v) { fAliEnAPIVersion = v; }
  /** 
   * Wether to use par files through-out. Mandetory and enforced in
   * case of a PROOF job,
   * 
   * @param usePar If true, use PAR files - even for base libraries 
   */
  void SetUsePar(Bool_t usePar) { fUsePar = usePar; }
  /* @} */

  //__________________________________________________________________
  /** 
   * @{ 
   * @name Input data  
   */
  /** 
   * Set the GRID/Local data dir 
   * 
   * @param d Directory with data 
   */
  void SetDataDir(const char* d) { fDataDir = d; }
  /** 
   * Set the glob pattern to search input files in - Grid only
   * 
   * @param pattern Glob pattern
   */
  void SetDataPattern(const char* pattern) { fDataPattern = pattern; }
  /** 
   * Set the PROOF data set 
   * 
   * @param d PROOF registered data set 
   */
  void SetDataSet(const char* d) { fDataSet = d; }
  /** 
   * Set the XML file to use 
   * 
   * @param x XML file 
   */
  void SetXML(const char* x) { fXML = x; }
  /** 
   * Wether to assume the input comes from MC.  If this is set to
   * true, and the CreateMCHandler member function isn't overloaded to
   * return null, then the files @c galice.root, @c Kinematics.root,
   * and @c TrackRefs.root must be present for each input file (@c
   * AliESDs.root or @c AliAOD.root)
   * 
   * @param isMC If true, assume MC input 
   */
  void SetMC(Bool_t isMC) { fMC = isMC; } 
  /* @} */

  //__________________________________________________________________
  /** 
   * @{ 
   * @name Grid storage and splitting 
   */
  /** 
   * Set how many replicas of the output we want 
   * 
   * @param n Number of replicas requested 
   */
  void SetNReplica(Int_t n) { fNReplica = n; }
  /** 
   * Set the maximum number of files per sub-job.  
   * 
   * @param max Maximum number of files per sub-job
   */  
  void SetMaxSplit(UShort_t max=50) { fMaxSplit = max; }
  /* @} */

  //__________________________________________________________________
  /** 
   * @{ 
   * @name Name and local working directory 
   */
  /** 
   * Set whether to allow overwritting existing files/directories 
   * 
   * @param allow If true, allow overwritting files/directories
   */
  void SetAllowOverwrite(Bool_t allow) { fAllowOverwrite = allow; }
  /** 
   * Set the date and time 
   * 
   * @param year     Year (>1994) 
   * @param month    Month
   * @param day      Day
   * @param hour     Hour
   * @param minutes  Minute
   */
  void SetDateTime(UShort_t year, UShort_t month, UShort_t day, 
		   UShort_t hour, UShort_t minutes)
  {
    fDatime.Set((year<1995?1995:year), month, day, hour, minutes, 0);
    fEscapedName = EscapeName(fName, fDatime);
  }
  /** 
   * Set the date and time from a string.
   * 
   * @param date Formatted like YYYY/MM/DD HH:MM:SS
   */
  void SetDateTime(const TString& date)
  {
    if (date.IsNull())
      fDatime.Set(1985,0,0,0,0,0);
    else if (date.EqualTo("now", TString::kIgnoreCase)) 
      fDatime.Set();
    else 
      fDatime.Set(date);
    fEscapedName = EscapeName(fName, fDatime);
  }
  /** 
   * Return the escaped name 
   * 
   * @return Escaped name 
   */
  const TString& EscapedName() const 
  {
    return fEscapedName;
  }
  /* @} */
  //__________________________________________________________________
  /** 
   * @{ 
   * @name Execution parameters 
   */
  /** 
   * Set the PROOF server URL
   * 
   * @param s PROOF server URL 
   */
  void SetProofServer(const char* s)    { fProofServer = s; }
  /** 
   * Set the type of analysis
   * 
   * @param type AOD or ESD
   */
  void SetType(EType type) { fExecType = type; }
  /** 
   * Set the type of analysis
   * 
   * @param type AOD or ESD
   */
  void SetType(const char* type) { SetType(ParseType(type)); }
  /** 
   * Set the execution mode of analysis
   * 
   * @param mode LOCAL, PROOF, GRID 
   */
  void SetMode(EMode mode) { fExecMode = mode; }
  /** 
   * Set the execution mode of analysis
   * 
   * @param mode LOCAL, PROOF, GRID 
   */
  void SetMode(const char* mode) { SetMode(ParseMode(mode)); }
  /** 
   * Set the execution operation of analysis
   * 
   * @param oper FULL, TERMINATE, INIT 
   */
  void SetOperation(EOper oper) { fExecOper = oper; }
  /** 
   * Set the execution operation of analysis
   * 
   * @param oper FULL, TERMINATE, INIT 
   */
  void SetOperation(const char* oper) { SetOperation(ParseOperation(oper)); }
  /** 
   * Use GDB to wrap PROOF slaves 
   * 
   * @param use Whether to use GDB or not 
   */
  void SetUseGDB(Bool_t use=kTRUE) { fUseGDB = use; }
  /* @} */

  //__________________________________________________________________
  /** 
   * @{ 
   * @name Stuff to upload 
   */
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
  /** 
   * Add binary data to be uploaded to slaves 
   * 
   * @param lib Name of binary file 
   */
  void AddLibrary(const char* lib) { fListOfLibraries.Add(new TObjString(lib));}
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
  /* @} */

  //__________________________________________________________________
  /** 
   * @{ 
   * @name Run numbers 
   */
  /** 
   * Add a run to be analysed
   *  
   * @param run Run number
   */
  void AddRun(Int_t run) 
  {
    Int_t i = fRunNumbers.fN; fRunNumbers.Set(i+1); fRunNumbers[i] = run;
  }
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
  /** 
   * Set the runs to read from a string.  The parts should be
   * delimited by a character in the string @a delim.  If a non-number
   * part is seen, it is assumed to be the name of a file containing
   * run numbers.
   * 
   * @param runs   String of runs
   * @param delim  Delimiters 
   */
  void SetRuns(const TString& runs, const char* delim=":, \t") 
  {
    TIter next(runs.Tokenize(delim));
    TObjString* os = 0;
    while ((os = static_cast<TObjString*>(next()))) {
      TString s(os->String());
      if (s.IsNull()) continue;
      if (!s.IsDigit()) ReadRunNumbers(s);
      else              AddRun(s.Atoi());
    }
  }
  /** 
   * Whether final merge should be done over all runs (argument true),
   * or for each run individually. 
   * 
   * @param perRun If true, do final merge over all runs 
   */
  void SetPerRunMerge(Bool_t perRun) { fPerRunMerge = perRun; }
  /* @} */
  //__________________________________________________________________
  /**
   * @{ 
   * @name Execution 
   */
  /** 
   * Initialize the job 
   * 
   * @return true on success, false otherwise
   */
  Bool_t Init()
  {
    if (fExecMode == kProof) fUsePar    = true;

    // Info("Init", "Connecting in mode=%d", mode);
    if (!Connect()) return false;

    // --- Get current directory and set-up sub-directory ------------
    TString cwd = gSystem->WorkingDirectory();
    if (!SetupWorkingDirectory()) return false;

    // --- Load the common libraries ---------------------------------
    if (!LoadCommonLibraries()) return false;
    
    // --- Create analysis manager -----------------------------------
    AliAnalysisManager *mgr  = new AliAnalysisManager(fName,"Analysis Train");

    // In test mode, collect system information on every event 
    // if (oper == kTest)  mgr->SetNSysInfo(1); 
    if (fVerbose  >  0)      mgr->SetDebugLevel(fVerbose);
    if (fExecMode == kLocal) mgr->SetUseProgressBar(kTRUE, 100);
   
    // --- ESD input handler ------------------------------------------
    AliVEventHandler*  inputHandler = CreateInputHandler(fExecType);
    if (inputHandler) mgr->SetInputEventHandler(inputHandler);
    
    // --- Monte-Carlo ------------------------------------------------
    AliVEventHandler*  mcHandler = CreateMCHandler(fExecType,fMC);
    if (mcHandler) mgr->SetMCtruthEventHandler(mcHandler);
    
    // --- AOD output handler -----------------------------------------
    AliVEventHandler*  outputHandler = CreateOutputHandler(fExecType);
    if (outputHandler) mgr->SetOutputEventHandler(outputHandler);
    
    // --- Include analysis macro path in search path ----------------
    gROOT->SetMacroPath(Form("%s:%s:$ALICE_ROOT/ANALYSIS/macros",
			     cwd.Data(), gROOT->GetMacroPath()));

    // --- Physics selction - only for ESD ---------------------------
    if (fExecType == kESD) CreatePhysicsSelection(fMC, mgr);
    
    // --- Create centrality task ------------------------------------
    CreateCentralitySelection(fMC, mgr);

    // --- Create tasks ----------------------------------------------
    CreateTasks(fExecMode, fUsePar, mgr);

    // --- Create Grid handler ----------------------------------------
    // _must_ be done after all tasks has been added
    AliAnalysisAlien* gridHandler = CreateGridHandler();
    if (gridHandler) mgr->SetGridHandler(gridHandler);
    
    // --- Print setup -----------------------------------------------
    Print();
    // if (mode == kProof) {
    // Info("Run", "Exported environment variables to PROOF slaves:");
    // TProof::GetEnvVars()->ls();
    // Info("Run", "Environment variables for this session:");
    // gSystem->Exec("printenv");
    // }

    // --- Initialise the train --------------------------------------
    if (!mgr->InitAnalysis())  {
      gSystem->ChangeDirectory(cwd.Data());
      Error("Run","Failed to initialise train");
      return false;
    }

    // --- Show status -----------------------------------------------
    mgr->PrintStatus();

    return true;
  }
  //------------------------------------------------------------------
  /** 
   * Run the analysis. 
   * 
   * @param nEvents Number of events to analyse 
   * @param r       Possible runner object 
   * @param asShell Passed to SaveSetup
   */
  virtual void Run(Int_t nEvents, Runner* r=0, Bool_t asShell=false)
  {
    // Info("Exec", "Doing exec with type=%d, mode=%d, oper=%d, events=%d "
    //      "mc=%d, usePar=%d", type, mode, oper, nEvents, mc, usePar);

    TString cwd = gSystem->WorkingDirectory();
    
    if (!Init()) { 
      Error("Run", "Failed to intialize the train");
      return;
    }
    if (r) SaveSetup(*r, nEvents, asShell);
    if (fExecOper == kInitialize) return;
    
    // --- Create the chain ------------------------------------------
    TChain* chain = CreateChain();
    if (fExecMode == kLocal) {
      if (!chain) {
	Error("Run", "No chain defined in local mode!");
	return;
      }
      if (chain->GetListOfFiles()->GetEntries() < 1) { 
	Error("Run", "Empty chain in local mode!");
	return;
      }
    }

    // --- Get manager and execute -----------------------------------
    AliAnalysisManager *mgr  =AliAnalysisManager::GetAnalysisManager();
    Long64_t ret = StartAnalysis(mgr, chain, nEvents);

    // Make sure we go back 
    gSystem->ChangeDirectory(cwd.Data());

    // Return. 
    if (ret < 0) Error("Exec", "Analysis failed");
  }
  //------------------------------------------------------------------
  /** 
   * Print the setup 
   * 
   */
  virtual void Print() const 
  {
    bool mc=AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
    std::cout << fName << " train setup\n"
	      << std::boolalpha;
    PrintField(std::cout, "Escaped name",               fEscapedName);
    PrintField(std::cout, "ROOT version",		fRootVersion);
    PrintField(std::cout, "AliROOT version",		fAliRootVersion);
    PrintField(std::cout, "AliEn API version",		fAliEnAPIVersion);
    PrintField(std::cout, "Name of proof server",	fProofServer);
    PrintField(std::cout, "Input directory",		fDataDir);
    PrintField(std::cout, "Data pattern",		fDataPattern);
    PrintField(std::cout, "Proof data set name",	fDataSet);
    PrintField(std::cout, "XML collection",		fXML);
    PrintField(std::cout, "Storage replication",	fNReplica);
    PrintField(std::cout, "Allow overwrite",            fAllowOverwrite);
    PrintField(std::cout, "Do GDB debugging",           fUseGDB);
    PrintField(std::cout, "Max # files per split",      fMaxSplit);
    PrintField(std::cout, "Monte-Carlo input",		fMC);
    PrintField(std::cout, "Monte-Carlo handler",        mc);
    PrintField(std::cout, "Per run merge",              fPerRunMerge);
    PrintFieldName(std::cout, "Run numbers");
    for (Int_t i = 0; i < fRunNumbers.GetSize(); i++) 
      std::cout << (i == 0 ? "" : ", ") << fRunNumbers.At(i);
    std::cout << std::endl;

    PrintFieldList(std::cout, "PAR files", 		fListOfPARs);
    PrintFieldList(std::cout, "Script sources", 	fListOfSources);
    PrintFieldList(std::cout, "Libraries", 		fListOfLibraries);
    PrintFieldList(std::cout, "Extras", 		fListOfExtras, "\n  ");

    std::cout << std::noboolalpha << std::endl;

    AliAnalysisGrid* plugin = 
      AliAnalysisManager::GetAnalysisManager()->GetGridHandler();
    if (!plugin) return;
    
  }
  /** 
   * Whether to be verbosity level.  0 means no messages, while higher
   * numbers increase the verbosity
   * 
   * @param verb Verbosity level 
   */
  void SetVerbose(Int_t verb) { fVerbose = verb; }
  /* @} */
protected:
  //__________________________________________________________________
  /** 
   * @{ 
   * @name Copying 
   */
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from
   */
  TrainSetup(const TrainSetup& o)
  : fName(o.fName),
    fEscapedName(o.fEscapedName),
    fRootVersion(o.fRootVersion),
    fAliRootVersion(o.fAliRootVersion),
    fAliEnAPIVersion(o.fAliEnAPIVersion),
    fProofServer(o.fProofServer),
    fDataDir(o.fDataDir),	
    fDataPattern(o.fDataPattern),
    fDataSet(o.fDataSet),	
    fXML(o.fXML),	
    fNReplica(o.fNReplica),
    fAllowOverwrite(o.fAllowOverwrite),
    fUseGDB(o.fUseGDB),
    fMaxSplit(o.fMaxSplit),
    fRunNumbers(o.fRunNumbers),
    fListOfPARs(),
    fListOfSources(),
    fListOfLibraries(),
    fListOfExtras(),
    fDatime(o.fDatime),
    fExecType(o.fExecType), 
    fExecMode(o.fExecMode), 
    fExecOper(o.fExecOper),
    fUsePar(o.fUsePar), 
    fMC(o.fMC), 
    fPerRunMerge(o.fPerRunMerge),
    fVerbose(o.fVerbose)
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
  //------------------------------------------------------------------
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
    fDataPattern        = o.fDataPattern;
    fDataSet		= o.fDataSet;	
    fXML		= o.fXML;	
    fNReplica		= o.fNReplica;	
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
  /* @} */

  //__________________________________________________________________
  /** 
   * @{ 
   * @name Utility functions 
   */
  /** 
   * Escape bad elements of the name 
   * 
   * @param name   Name to escape 
   * @param datime Date and Time 
   * 
   * @return escaped name 
   */  
  static TString EscapeName(const char* name, const TDatime& datime)
  {
    TString escaped = name;
    char  c[] = { ' ', '/', '@', 0 };
    char* p   = c;
    while (*p) { 
      escaped.ReplaceAll(Form("%c", *p), "_");
      p++;
    }
    if (datime.GetYear() <= 1995 ||
	datime.GetMonth() == 0 || 
	datime.GetDay() == 0) return escaped;
    escaped.Append(Form("_%04d%02d%02d_%02d%02d", 
			datime.GetYear(), 
			datime.GetMonth(), 
			datime.GetDay(), 
			datime.GetHour(), 
			datime.GetMinute()));
    return escaped;
  }    
  //------------------------------------------------------------------
  static void PrintFieldName(std::ostream& o, const char* name)
  {
    o << "  " << std::left << std::setw(25) << name << ": " << std::flush;
  }
  //------------------------------------------------------------------
  static void PrintFieldList(std::ostream& o, const char* name, 
			     const TCollection& c, const char* sep=", ")
  {
    PrintFieldName(o, name);
    Bool_t   first = true;
    TObject* obj = 0;
    TIter    next(&c);
    while ((obj = next())) {
      o << (first ? "" : sep) << obj->GetName();
      first = false;
    }
    std::cout << std::endl;
  }
  //------------------------------------------------------------------
  template <typename T>
  static void PrintField(std::ostream& o, const char* name, T& value) 
  {
    PrintFieldName(o, name);
    o << value << std::endl;
  }
  //------------------------------------------------------------------
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
    case kInitialize:   return "INIT";
    }
    return 0;
  }
  //------------------------------------------------------------------
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
    if      (sOper.Contains("TEST"))      eOper = kTest;
    else if (sOper.Contains("OFFLINE"))   eOper = kOffline;
    else if (sOper.Contains("SUBMIT"))    eOper = kSubmit;
    else if (sOper.Contains("TERMINATE")) eOper = kTerminate;
    else if (sOper.Contains("FULL"))      eOper = kFull;
    else if (sOper.Contains("INIT"))      eOper = kInitialize;
    else 
      Fatal("Run", "unknown operation '%s'", oper);
    return eOper;
  }
  //------------------------------------------------------------------
  /** 
   * Return a string that reflects the passed mode
   * 
   * @param eType Type of analysis 
   * 
   * @return String representation of execution type
   */
  static const char* TypeString(EType eType) 
  {
    switch (eType) {
    case kESD:	return "ESD";
    case kAOD:	return "AOD";
    case kUser:	return "USER";
    }
    return 0;
  }
  //------------------------------------------------------------------
  /** 
   * Parse a string into a type enum
   * 
   * @param type String to pass
   * 
   * @return Enumaration value 
   */
  static EType ParseType(const char* type)
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
  //------------------------------------------------------------------
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
  //------------------------------------------------------------------
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
    if      (sMode.Contains("LOCAL")) eMode = kLocal;
    else if (sMode.Contains("PROOF")) eMode = kProof;
    else if (sMode.Contains("GRID"))  eMode = kGrid;
    else 
      Fatal("Run", "Unknown mode '%s'", mode);
    return eMode;
  }
  /* @} */


  //__________________________________________________________________
  /** 
   * @{ 
   * @name Overloadable creators 
   */
  /** 
   * Create a grid handler 
   * 
   * @return Grid handler 
   */
  virtual AliAnalysisAlien* 
  CreateGridHandler()
  {
    if (fExecMode != kGrid) return 0;

    TString name = EscapedName();

    // Create the plug-in object, and set run mode 
    AliAnalysisAlien* plugin = new AliAnalysisAlien();
    plugin->SetRunMode(OperString(fExecOper == kInitialize ? 
				  kFull : fExecOper));
    
    // Production mode - not used here 
    // plugin->SetProductionMode();
    
    // Set output to be per run 
    plugin->SetOutputToRunNo(true); 

    // Set the job tag 
    plugin->SetJobTag(fName);

    // Set number of test files - used in test mode only 
    plugin->SetNtestFiles(1);

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
    plugin->SetSplitMaxInputFileNumber(fMaxSplit);

    // Disable default outputs 
    plugin->SetDefaultOutputs(true);

    // Merge parameters 
    plugin->SetMaxMergeFiles(20);
    plugin->SetMergeExcludes("AliAOD.root "
			    "*EventStat*.root "
			    "*event_stat*.root");
    
    // Keep log files 
    plugin->SetKeepLogs();

    // Set the working directory to be the trains name (with special
    // characters replaced by '_' and the date appended), and also set
    // the output directory (relative to working directory)
    plugin->SetGridWorkingDir(name.Data());
    plugin->SetGridOutputDir("output");

    // Set required version of software 
    if (!fAliEnAPIVersion.IsNull()) plugin->SetAPIVersion(fAliEnAPIVersion);
    if (!fRootVersion.IsNull())     plugin->SetROOTVersion(fRootVersion);
    if (!fAliRootVersion.IsNull())  plugin->SetAliROOTVersion(fAliRootVersion);

    // Declare root of input data directory 
    TString dataDir(fDataDir);
    if (dataDir.BeginsWith("alien://")) 
      dataDir.ReplaceAll("alien://", "");
    plugin->SetGridDataDir(dataDir);

    // Data search patterns 
    TString pat;
    if (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())
      plugin->SetRunPrefix("");
    else {
      plugin->SetRunPrefix("000");
    }
    pat = fDataPattern;
    if (!pat.EndsWith("/")) pat.Append("/");
    pat.Append(Form("*%s.root", fExecType == kESD ? "ESDs" : "AOD"));
    plugin->SetDataPattern(pat);

    // Add the run numbers 
    Int_t nRun = 0;
    for (Int_t i = 0; i < fRunNumbers.fN; i++) {
      if (fRunNumbers[i] < 0) continue; 
      plugin->AddRunNumber(fRunNumbers[i]);
      nRun++;
    }
    // Set number of runs per master - set to one to per run
    if (fPerRunMerge) plugin->SetNrunsPerMaster(1);
    else              plugin->SetNrunsPerMaster(nRun+1);
    
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
    

    return plugin;
  }
  //------------------------------------------------------------------
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
    case kUser: return 0;
    }
    return 0;
  }
  //------------------------------------------------------------------
  /** 
   * Create input handler 
   * 
   * @param type  Run type (ESD or AOD)
   * @param mc    Assume monte-carlo input 
   * 
   * @return 
   */
  virtual AliVEventHandler* CreateMCHandler(EType /*type*/, bool mc)
  {
    if (!mc)          return 0;
    // if (type != kESD) return 0;
    Info("CreateMCHandler", "Making MC handler");
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mcHandler->SetReadTR(true); 
    return mcHandler;
  }
  //------------------------------------------------------------------
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
    case kUser: 
      break;
    }
    
    return ret;
  }
  //------------------------------------------------------------------
  /** 
   * Create physics selection, and add to manager
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
  //------------------------------------------------------------------
  /** 
   * Create centrality selection, and add to manager
   * 
   * @param mc Whether this is for MC 
   * @param mgr Manager
   */
  virtual void CreateCentralitySelection(Bool_t mc, AliAnalysisManager* mgr)
  {
    gROOT->Macro("AddTaskCentrality.C");
    const char* name = "CentralitySelection";
    AliCentralitySelectionTask* ctask = 
      dynamic_cast<AliCentralitySelectionTask*>(mgr->GetTask(name));
    if (!ctask) return;
    if (mc) ctask->SetMCInput();
  }
  //------------------------------------------------------------------
  /** 
   * Create analysis tasks.  Must be overloaded by sub-class
   * 
   * @param mode Run mode
   * @param mgr  Manager
   * @param par  Whether to use pars 
   */
  virtual void CreateTasks(EMode mode, Bool_t par, AliAnalysisManager* mgr)=0;
  /* @} */

  //__________________________________________________________________
  /** 
   * @{ 
   * @name Library loading 
   */
  //------------------------------------------------------------------
  /** 
   * Load common libraries 
   * 
   * @return true on success 
   */
  Bool_t LoadCommonLibraries() 
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
    if (fExecMode == kProof) { 
      gProof->Exec("gSystem->Load(\"libTree.so\");");
      gProof->Exec("gSystem->Load(\"libGeom.so\");");
      gProof->Exec("gSystem->Load(\"libMinuit.so\");");
      gProof->Exec("gSystem->Load(\"libVMC.so\");");

      
    }

    Bool_t ret   = true;
    Bool_t basic = fExecMode == kGrid ? false : fUsePar;
    
    ret &= LoadLibrary("STEERBase",     basic, false);
    ret &= LoadLibrary("ESD",           basic, false);
    ret &= LoadLibrary("AOD",           basic, false);
    ret &= LoadLibrary("ANALYSIS",      basic, true);
    ret &= LoadLibrary("OADB",          basic, true);
    ret &= LoadLibrary("ANALYSISalice", basic, true);

    return ret;
  }
  //------------------------------------------------------------------
  /** 
   * Load a library 
   * 
   * @param what What library to load
   * @param par  If true, load as PAR
   * @param rec  If true, also load on slaves
   * 
   * @return true on success 
   */
  Bool_t LoadLibrary(const char* what, Bool_t par, Bool_t rec=false)
  {
    if (!what || what[0] == '\0') return true;
    
    TString module(what);
    TString libName(what);
    if (!libName.BeginsWith("lib")) { 
      // Check if the library corresponds to a compiled macro 
      if (!gSystem->AccessPathName(Form("%s_C.so", libName.Data()))) {
	libName.Append("_C");
      }
      else if (!gSystem->AccessPathName(Form("../%s_C.so", libName.Data()))) {
	libName = Form("../%s_C", what);
      }
      else 
	libName = Form("lib%s", libName.Data());
    }
    if (!libName.EndsWith(".so"))   libName.Append(".so");

    Int_t ret = 0;

    switch (fExecMode) { 
    case kLocal: // Just load and exit 
      if (gSystem->Load(libName.Data()) < 0) {
	Error("LoadLibrary", "Failed to load library %s", libName.Data());
	return false;
      }
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
      Info("LoadLibrary", "Uploading %s", what);
      ret = gProof->UploadPackage(what, TProof::kRemoveOld);
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
      Info("LoadLibrary", "Enabling package %s", what);
      ret = gProof->EnablePackage(what);
      break;
    }
    if (ret < 0) { 
      Error("LoadLibrary", "Couldn't load %s", what);
      return false;
    }
    return true;
  }
  /* @} */

  //__________________________________________________________________
  /** 
   * @{ 
   * @name PAR generation from script 
   */
  /** 
   * Service function to make a PAR out of a script.  
   * 
   * The script should contain can contain a sub-class of AliAnalysisTask. 
   * The script will be compiled on the slaves before loading the 
   * AliAnalysisManager.  Parts to (not) be compiled can be protected like 
   * 
   * @code 
   * #ifdef BUILD_PAR
   * // This will _only_ be compiled in the servers 
   * #endif
   * #ifndef BUILD_PAR
   * // This will not be compiled in the servers 
   * #endif
   * @endcode
   * 
   * @param mode   Execution mode (Grid, PROOF, Local)
   * @param script Script to upload and compile in the PAR
   * @param deps   Dependency pars 
   * 
   * @return true on success. 
   */
  static Bool_t MakeScriptPAR(EMode mode, const char* script, const char* deps)
  {
    // Get the base name 
    Info("MakeScriptPAR", "Making par file for %s", script);
    TString base(gSystem->BaseName(script));
    Int_t   idx = base.Last('.');
    if (idx != kNPOS) base.Remove(idx);
    Bool_t retval = true;
    // Info("MakeScriptPAR", "script=%s, base=%s", script, base.Data());

    if (mode == kLocal) { 
      if (gROOT->LoadMacro(Form("%s.C++g", base.Data())) < 0)
	return false;
      return true;
    }

    TString tmpdir(gSystem->TempDirectory());
    int   ltempl = tmpdir.Length() + 1 + 5 + 6 + 1;
    char* templ  = new char[ltempl];
    snprintf(templ, ltempl, "%s/trainXXXXXX", tmpdir.Data());
    if (!mkdtemp(templ)) {
      Error("MakeScriptPAR", 
	    "Failed to generate temporary directory from template %s", 
	    templ);
      return false;
    }

    try {
      // Check name of script file 
      TString scr(script);
      TString ext;
      if      (scr.EndsWith(".C"))   ext = "C"; 
      else if (scr.EndsWith(".cxx")) ext = "cxx";
      else                           { ext = "C"; scr.Append(".C"); }
      
      // Check if we can access the file 
      TString path = TString::Format(".:%s", TROOT::GetMacroPath());
      char* loc = gSystem->Which(path, scr);
      if (!loc) throw TString::Format("Script %s not found in %s", 
				      scr.Data(), path.Data());
      TString full(loc);
      
      TString dir = TString::Format("%s/%s", templ, base.Data());
      // Set-up directories 
      if (gSystem->MakeDirectory(dir) < 0) 
	throw TString::Format("Could not make directory '%s'", base.Data());
      
      if (gSystem->MakeDirectory(Form("%s/PROOF-INF", dir.Data()))) 
	throw TString::Format("Could not make directory %s/PROOF-INF", 
			      base.Data());
      
      // Copy the script to the setup directory 
      TString dest = TString::Format("%s/%s.%s", dir.Data(),
				     base.Data(), ext.Data());
      Int_t ret = gSystem->CopyFile(full, dest, true);
      switch (ret) { 
      case -1: throw TString::Format("Couldn't open %s for copy", scr.Data());
      case -2: throw TString::Format("File %s exists", dest.Data());
      case -3: throw TString::Format("Error while copying %s", scr.Data());
      }
      
      // Make our build file 
      // Info("MakeScriptPAR", "Making build script %s/PROOF-INF/BUILD.sh", dir.Data());
      std::ofstream b(Form("%s/PROOF-INF/BUILD.sh", dir.Data()));
      if (!b) 
	throw TString::Format("Failed to open b shell script");
      b << "#!/bin/sh\n"
	<< "echo BUILD.sh@`hostname`: Building " << base << "\n"
	<< "root.exe -l -b -q PROOF-INF/BUILD.C 2>&1 | tee " << base << ".log\n"
	<< "echo BUILD.sh@`hostname`: done: $?\n"
	<< std::endl;
      b.close();
      if (gSystem->Chmod(Form("%s/PROOF-INF/BUILD.sh", dir.Data()), 0755) != 0)
	throw TString::Format("Failed to set exectuable flags on "
			      "%s/PROOF-INF/BUILD.sh", dir.Data());
      
      // Info("MakeScriptPAR", "Making utility script %s/PROOF-INF/UTIL.C", dir.Data());
      std::ofstream u(Form("%s/PROOF-INF/UTIL.C", dir.Data()));
      if (!u) 
	throw TString::Format("Failed to open utility script");
      u << "void LoadROOTLibs() {\n"
	<< "  gSystem->Load(\"libVMC\");\n"
	<< "  gSystem->Load(\"libNet\");\n"
	<< "  gSystem->Load(\"libTree\");\n"
	<< "  gSystem->Load(\"libPhysics\");\n"
	<< "  gSystem->Load(\"libMinuit\");\n"
	<< "}\n\n"
	<< "void AddAliROOT() {\n"
	<< "  TString val(gSystem->Getenv(\"ALICE_ROOT\"));\n"
	<< "  if (val.IsNull())\n"
	<< "    Warning(\"Add\",\"ALICE_ROOT not defined\");\n"
	<< "  else\n"
	<< "    gSystem->AddIncludePath(Form(\"-I%s/include\",val.Data()));\n"
	<< "}\n\n"
	<< "void AddDep(const char* env) {\n"
	<< "  TString val(gSystem->Getenv(Form(\"%s_INCLUDE\",env)));\n"
	<< "  if (val.IsNull())\n"
	<< "    Warning(\"Add\",\"%s_INCLUDE not defined\", env);\n"
	<< "  else {\n"
	<< "    gSystem->AddIncludePath(Form(\"-I../%s\",val.Data()));\n"
	<< "  }\n"
	<< "}\n\n"
	<< "void LoadDep(const char* name) {\n"
	<< "  gSystem->AddDynamicPath(Form(\"../%s\",name));\n"
	<< "  char* full = gSystem->DynamicPathName(name,true);\n"
	<< "  if (!full) \n"
	<< "   full = gSystem->DynamicPathName(Form(\"lib%s\",name),true);\n"
	<< "  if (!full) \n"
	<< "   full = gSystem->DynamicPathName(Form(\"lib%s.so\",name),true);\n"
	<< "  if (!full) {\n"
	<< "    Warning(\"LoadDep\",\"Module %s not found\", name);\n"
	<< "    return;\n"
	<< "  }\n"
	<< "  gSystem->Load(full);\n"
	<< "}\n"
	<< std::endl;
      u.close();

      // Info("MakeScriptPAR", "Making utility script %s/PROOF-INF/BUILD.C", dir.Data());
      std::ofstream cbuild(Form("%s/PROOF-INF/BUILD.C", dir.Data()));
      if (!cbuild) 
	throw TString::Format("Failed to open build script");
      cbuild << "void BUILD() {\n"
	     << "  gSystem->AddIncludePath(\"-DBUILD_PAR=1\");\n"
	     << "  gROOT->LoadMacro(\"PROOF-INF/UTIL.C\");\n"
	     << "  LoadROOTLibs();\n"
	     << "  AddAliROOT();\n";
      TObjArray*  depList = TString(deps).Tokenize(",");
      TIter       next(depList);
      TObject*    dep = 0;
      while ((dep = next())) {
	cbuild << "  AddDep(\"" << dep->GetName() << "\");\t"
	       << "  LoadDep(\"" << dep->GetName() << "\");\n";
      }
      cbuild << "  // gDebug = 5;\n"
	     << "  int ret = gROOT->LoadMacro(\"" 
	     << base << "." << ext << "++g\");\n"
	     << "  if (ret != 0) Fatal(\"BUILD\",\"Failed to build\");\n"
	     << "  // else Info(\"BUILD\", \"Made " << base << "\");\n"
	     << "}\n"
	     << std::endl;
      cbuild.close();
      
      // Make our set-up script 
      // Info("MakeScriptPAR", "Making setup script %s/PROOF-INF/SETUP.C", dir.Data());
      std::ofstream setup(Form("%s/PROOF-INF/SETUP.C", dir.Data()));
      if (!setup) 
	throw TString::Format("Failed to open setup script");
      setup << "void SETUP() {\n"
	    << "  gROOT->LoadMacro(\"PROOF-INF/UTIL.C\");\n"
	    << "  LoadROOTLibs();\n"
	    << "  // Info(\"SETUP\",\"Loading libraries\");\n";
      next.Reset();
      dep = 0;
      while ((dep = next())) 
	setup << "  LoadDep(\"" << dep->GetName() << "\");\n";
      setup << "  // gDebug = 5;\n"
	    << "  // Info(\"SETUP\",\"Loading " << base << "_" << ext << ".so\");\n"
	    << "  gSystem->Load(\"" << base << "_" << ext << ".so\");\n"
	    << "  // gDebug = 0;\n"
	    << "  gROOT->ProcessLine(\".include " << base << "\");\n"
	    << "  gSystem->Setenv(\"" << base << "_INCLUDE\",\"" 
	    << base << "\");\n"
	    << "  // Info(\"SETUP\", \"Done\");\n"
	    << "}\n"
	    << std::endl;
      setup.close();

      // Info("MakeScriptPAR", "Packing up tar-archive");
      ret = gSystem->Exec(Form("(cd %s && tar -czf %s.par %s)", 
			       templ, base.Data(),base.Data()));
      if (ret != 0) 
	throw TString::Format("Failed to create PAR file %s.PAR from %s", 
			      base.Data(), dir.Data());

      // Info("MakeScriptPAR", "Moving PAR archive");
      ret = gSystem->Exec(Form("mv -f %s/%s.par %s.par", templ, base.Data(), 
			       base.Data()));
      if (ret != 0) 
	throw TString::Format("Failed to rename %s/%s.par to %s.par: %s", 
			      templ, base.Data(), base.Data(), 
			      gSystem->GetError());
    }
    catch (TString& e) { 
      Error("MakeScriptPAR", "%s", e.Data()); 
      retval = false;
    }
    // Info("MakeScriptPAR", "Removing temperary directory %s", templ);
    gSystem->Exec(Form("rm -rf %s", templ));
    return retval;
  }
  /* @} */

  //__________________________________________________________________
  /** 
   * @{ 
   * @name Execution implementation
   */
  /** 
   * Start the analysis 
   * 
   * @param mgr       Analysis manager
   * @param chain     Input data (local and proof only)
   * @param nEvents   Number of events to analyse 
   */
  Long64_t StartAnalysis(AliAnalysisManager* mgr, 
			 TChain*             chain,
			 Int_t               nEvents)
  {
    // --- Run the analysis ------------------------------------------
    TString mode = ModeString(fExecMode);
    switch (fExecMode) { 
    case kLocal: 
      if (!chain) {
	Error("StartAnalysis", "No chain defined");
	return -1;
      }
      if (nEvents < 0) nEvents = chain->GetEntries();
      return mgr->StartAnalysis(mode, chain, nEvents);
    case kProof: 
      if (fDataSet.IsNull()) {
	if (!chain) { 
	  Error("StartAnalysis", "No chain defined");
	  return -1;
	}
	if (nEvents < 0) nEvents = chain->GetEntries();
	return mgr->StartAnalysis(mode, chain, nEvents);
      }
      return mgr->StartAnalysis(mode, fDataSet);
    case kGrid: 
      if (nEvents < 0)
	return mgr->StartAnalysis(mode);
      return mgr->StartAnalysis(mode, nEvents);
    }
    // We should never get  here 
    return -1;
  }
  //------------------------------------------------------------------
  /** 
   * Connect to external services (Proof and/or grid)
   * 
   * @return true on success 
   */
  virtual Bool_t Connect()
  {
    if (fExecMode == kLocal) return true;
			  
    // --- Set-up connections to Proof cluster and alien -------------
    if (fExecMode == kProof) { 
      Info("Connect", "Opening connection to proof server");
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

      // --- Figure out some server settings -------------------------
      TString serv = "";
      Bool_t  lite = false;
      if (fProofServer.BeginsWith("workers=") || fProofServer.IsNull()) {
	lite = true;
	serv = fProofServer;
      }
      else 
	serv = Form("%s@%s", userName.Data(), fProofServer.Data());

      // --- Possibly debug slave sessions with GDB ------------------
      if (fUseGDB) { 
	TString gdbCmd("/usr/bin/gdb --batch -ex run -ex bt --args");
	// TString gdbCmd("\"gdb --batch -ex run -ex bt --args\"");
	Info("Connect", "Using GDB to wrap slaves: %s", gdbCmd.Data());
	TProof::AddEnvVar("PROOF_WRAPPERCMD", gdbCmd);
      }
      
      // --- Add ALICE_ROOT directory to search path for packages ----
      Info("Connect", "Set location of packages");
      gEnv->SetValue("Proof.GlobalPackageDirs", 
		     Form("%s:%s", 
			  gEnv->GetValue("Proof.GlobalPackageDirs", "."), 
			  gSystem->Getenv("ALICE_ROOT")));
						     
      // --- Set OADB path on workers --------------------------------
      const char* oadbPath = AliAnalysisManager::GetOADBPath();
      TProof::AddEnvVar("OADB_PATH", oadbPath);
      // if (lite) gSystem->Setenv("OADB_PATH", oadbPath);
      // Info("Connect", "OADB_PATH=%s", gSystem->Getenv("OADB_PATH"));

      // --- Now open connection to PROOF cluster --------------------
      TProof::Open(serv);
      if (!gProof) { 
	Error("Connect", "Failed to connect to Proof cluster %s as %s",
	      fProofServer.Data(), userName.Data());
	return false;
      }
      Info("Connect", "Now connected to Proof");
      // gProof->SetParameter("PROOF_LookupOpt", "all");
      if (lite) return true;
    }

    // --- Open a connection to the grid -----------------------------

    TGrid::Connect("alien://");
    if (!gGrid || !gGrid->IsConnected()) { 
      // This is only fatal in grid mode 
      Error("Connect", "Failed to connect to AliEN");
      if (fExecMode == kGrid) return false; 
      return true;
    }
    if (fExecMode == kGrid) return true;

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
  //------------------------------------------------------------------
  /** 
   * Get the output directory (local or Grid)
   * 
   * @param mode Mode of execution 
   * 
   * @return Path to output directory 
   */
  TString GetOutputDirectory(EMode mode) const 
  {
    TString ret(fEscapedName);
    if (mode != kGrid) return ret;
    
    AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
    if (!am) { 
      Warning("GetOutputDirectory", "No analysis manager defined yet");
      return ret;
    }
    AliAnalysisGrid* ag = am->GetGridHandler();
    if (!ag) { 
      Warning("GetOutputDirectory", "No grid handler defined yet");
      return ret;
    }
    AliAnalysisAlien* aa = dynamic_cast<AliAnalysisAlien*>(ag);
    if (!aa) { 
      Warning("GetOutputDirectory", "Grid handler isn't for AliEn");
      return ret;
    }
    ret = aa->GetGridOutputDir();
    if (!ret.BeginsWith("/")) {
      if (gGrid)
	ret = Form("%s/%s/%s", gGrid->GetHomeDirectory(), 
		   fEscapedName.Data(), aa->GetGridOutputDir());
      else 
	ret = Form("%s/%s", fEscapedName.Data(), aa->GetGridOutputDir());
    }
    return ret;
  }
  /* @} */

  //__________________________________________________________________
  /** 
   * @{ 
   * @name Setup 
   */
  /** 
   * Make our working directory if so requested 
   * 
   * @return true on success
   */
  Bool_t SetupWorkingDirectory()
  {
    TString nam = EscapedName();
    //Info("Init","Current dir=%s, escaped name=%s",cwd.Data(),nam.Data());
    Bool_t exists = gSystem->AccessPathName(nam.Data()) == 0;
    if (fExecOper == kTerminate && !exists) {
      Error("SetupWorkingDirectory", "File/directory %s does not exists", 
	    nam.Data());
      return false;
    }

	
    if (!fAllowOverwrite && exists) {
      Error("SetupWorkingDirectory", "File/directory %s already exists", 
	    nam.Data());
      return false;
    }

    if (!exists) {
      if (gSystem->MakeDirectory(nam.Data())) {
	Error("SetupWorkingDirectory", "Failed to make directory '%s'", 
	      nam.Data());
	return false;
      }
    }
      
    if (!gSystem->ChangeDirectory(nam.Data())) { 
      Error("SetupWorkingDirectory", "Failed to change directory to %s", 
	    nam.Data());
      return false;
    }
    Info("SetupWorkingDirectory", "Made subdirectory %s, and cd'ed there", 
	 nam.Data());
    return true;
  }
  //------------------------------------------------------------------
  /** 
   * Set-up a PAR file 
   * 
   * @param what PAR file 
   * 
   * @return true on success
   */
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
    gSystem->Exec(Form("tar xzf %s", parFile.Data()));
    
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
      // Info("SetupPAR", "Setting up for PAR %s", what);
      gROOT->Macro("PROOF-INF/SETUP.C");
    }
    if (!gSystem->ChangeDirectory(cwd.Data())) return false;

    return true;
  }
  //------------------------------------------------------------------
  /** 
   * Set-up extra sources. 
   * 
   * @return true on success
   */
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
  //------------------------------------------------------------------
  /** 
   * Set-up sources for upload 
   * 
   * 
   * @return String of sources 
   */
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
  //------------------------------------------------------------------
  /** 
   * Set-up extra libraries to upload 
   * 
   * @return String of libraries 
   */
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
  /* @} */

  //__________________________________________________________________
  /** 
   * @{ 
   * @name Chain building 
   */
  /** 
   * Check if we can add a file to the chain 
   * 
   * @param path   Full path to file 
   * @param chain  Chain 
   * 
   * @return true on success, false otherwise
   */
  Bool_t CheckFile(const TString& path, TChain* chain)
  {
    TFile* test = TFile::Open(path, "READ");
    if (!test) { 
      Warning("CheckFile", "Failed to open %s", path.Data());
      return false;
    }

    Bool_t ok = false; // Assume failure
    TObject* o = test->Get(chain->GetName());
    if (!o) 
      Warning("CheckFile", "The file %s does not contain the object %s", 
	      path.Data(), chain->GetName());
    else if (!dynamic_cast<TTree*>(o)) 
      Warning("CheckFile", "Object %s found in %s is not a TTree", 
	      o->GetName(), path.Data());
    else 
      ok = true;
    test->Close();
    if (ok) chain->AddFile(path);

    return ok;
  }
  /** 
   * Scan directory @a dir (possibly recursive) for tree files to add
   * to the chain.    This does not follow sym-links
   * 
   * @param dir        Directory to scan
   * @param chain      Chain to add to
   * @param recursive  Whether to scan recursively 
   *
   * @return true if any files where added 
   */
  Bool_t ScanDirectory(TSystemDirectory* dir, TChain* chain, 
		       bool recursive)
  {
    TString fnPattern;
    switch (fExecType) { 
    case kESD:  fnPattern = "AliESD"; break;
    case kAOD:  fnPattern = "AliAOD"; break;
    case kUser: fnPattern = "";       break;
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
    Bool_t hasGAlice = (!fMC ? true : false);
    Bool_t hasKine   = (!fMC ? true : false);
    Bool_t hasTrRef  = (!fMC ? true : false);
    
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
          if (ScanDirectory(d,chain,recursive))
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
	if (fMC) { 
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

    if (fMC && toAdd.GetEntries() > 0 && 
	(!hasGAlice || !hasKine || !hasTrRef)) { 
      Warning("ScanDirectory", 
	      "one or more of {galice,Kinematics,TrackRefs}.root missing from "
	      "%s, not adding anything from this directory", 
	      dir->GetTitle());
      toAdd.Delete();
    }

    TIter nextAdd(&toAdd);
    TObjString* s = 0;
    Int_t added = 0;
    while ((s = static_cast<TObjString*>(nextAdd()))) {
      // Info("ScanDirectory", "Adding %s", s->GetString().Data());
      TString fn = s->GetString();
      if (!CheckFile(fn, chain)) continue;

      added++;
    }
    if (added > 0) ret = true;

    gSystem->ChangeDirectory(oldDir);
    return ret;
  }
  //------------------------------------------------------------------
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
  //------------------------------------------------------------------
  /** 
   * Create a chain of data 
   *
   * @return TChain of data 
   */    
  TChain* CreateChain()
  {
    TString treeName;
    switch (fExecType) { 
    case kESD:  treeName = "esdTree"; break;
    case kAOD:  treeName = "aodTree"; break;
    case kUser: treeName = "";        break;
    }

    TChain* chain = 0;
    switch (fExecMode) { 
    case kProof: 
      if (!fDataSet.IsNull()) break; 
      // Otherwise fall through
    case kLocal:
      if (fXML.IsNull()) {
	chain = new TChain(treeName.Data());
	TString dir(fDataDir);
	if (dir == ".") dir = "";
	if (!dir.BeginsWith("/")) dir = Form("../%s", dir.Data());
	FileStat_t stat; 
	gSystem->GetPathInfo(dir, stat);
	if (!R_ISDIR(stat.fMode)) { // A file, check it 
	  if (!CheckFile(dir, chain)) { 
	    delete chain;
	    chain = 0;
	  }
	  break;
	}
	TString savdir(gSystem->WorkingDirectory());
	TSystemDirectory d(gSystem->BaseName(dir.Data()), dir.Data());
	if (!ScanDirectory(&d, chain, true)) { 
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
  /* @} */

public:
  //====================================================================
  /**
   * Option class 
   * 
   */
  struct Option
  {
    /** 
     * Constructor 
     * 
     * @param name Option name
     * @param desc Option description 
     * @param arg  Option argument, if any
     */
    Option(const char* name, const char* desc, const char* arg="")
      : fName(name),
	fDesc(desc),
	fArg(arg),
	fIsSet(false),
	fValue("")
    {}
    /** 
     * Process an option string. 
     * 
     * @param opt String to process
     * 
     * @return true, if this handled the option, false otherwise
     */
    Bool_t Process(const TString& opt)
    {
      // Info("Option::Process", "Option %s processing %s",
      //      fName.Data(), opt.Data());
      if (!opt.BeginsWith(fName, TString::kIgnoreCase)) return false;

      // We've got a value 
      fIsSet = true;
    
      // No argument options shouldn't set a value 
      if (fArg.IsNull()) return true;

      // Parse out value 
      Int_t eq = opt.Index("=");

      // Empty string is an OK value 
      if (eq == kNPOS) return true;

      TString tmp = opt(eq+1,opt.Length()-eq-1);
      fValue = tmp.Strip();
    
      return true;
    }
    /** 
     * Get option value as a string 
     * 
     * @return Option value
     */
    const TString& AsString() const { return fValue; } 
    /** 
     * Get option value as a double
     * 
     * @return Option value
     */
    Double_t AsDouble() const { return fValue.Atof(); }
    /** 
     * Get option value as an integer 
     * 
     * @return Option value
     */
    Int_t AsInt() const { return fValue.Atoi(); } 
    /** 
     * Get option value as a boolean
     * 
     * @return Option value
     */
    Bool_t   AsBool() const { return fIsSet; }
    /** 
     * Test if the option has been set.
     * 
     * @return True if the option was given
     */
    Bool_t   IsSet() const { return fIsSet; }
    /** 
     * Print help. 
     * 
     * @param o Stream to write on
     * @param prefix Prefix
     */
    void PrintHelp(std::ostream& o, const char* prefix) const 
    {
      TString arg(fName);
      if (!fArg.IsNull()) arg.Append(Form("=%s", fArg.Data()));
      o << "  " << (prefix ? prefix : "") 
	<< std::left << std::setw(30) << arg 
	<< " " << fDesc << std::endl;
    }
    /** 
     * Print the setting 
     * 
     * @param o Stream to write on
     */
    void PrintSettings(std::ostream& o) const 
    {
      o << "  " << std::left << std::setw(30) << fName << ": ";
      if (fArg.IsNull()) o << (IsSet() ? "true" : "false");
      else               o << fValue;
      o << std::endl;
    }
    /** 
     * Save the setting 
     * 
     * @param str  object nmae 
     * @param val  Value
     * @param o Stream to write on
     */
    void Save(std::ostream& o, const char* str, bool val)
    {
      if (!val) return;
      if (str[0] == '-') {
	o << "  " << str << fName << " \\" << std::endl;
	return;
      }
      o << "  " << str << ".Append(\"" << fName <<  ",\");" << std::endl;
    }
    /** 
     * Save the setting 
     * 
     * @param str  object nmae 
     * @param val  Value
     * @param o Stream to write on
     */
    void Save(std::ostream& o, const char* str, Int_t val)
    {
      if (str[0] == '-') { 
	o << "  " << str << fName << "=" << val << " \\" << std::endl;
	return;
      }
      o << "  " << str << ".Append(\"" << fName <<  "=" << val 
	<< ",\");" << std::endl;
    }
    /** 
     * Save the setting 
     * 
     * @param str  object nmae 
     * @param val  Value
     * @param o Stream to write on
     */
    void Save(std::ostream& o, const char* str, Double_t val)
    {
      if (str[0] == '-') { 
	o << "  " << str << fName << "=" << val << "  \\" << std::endl;
	return;
      }
      o << "  " << str << ".Append(\"" << fName <<  "=" << val 
	<< ",\");" << std::endl;
    }
    /** 
     * Save the setting 
     * 
     * @param str  object nmae 
     * @param val  Value
     * @param o Stream to write on
     */
    void Save(std::ostream& o, const char* str, const char* val)
    {
      if (str[0] == '-') { 
	TString sval(val);
	sval.ReplaceAll(" ", "\\ ");
	o << "  " << str << fName << "=" << sval << " \\" << std::endl;
	return;
      }
      o << "  " << str << ".Append(\"" << fName <<  "=" << val 
	<< ",\");" << std::endl;
    }

    TString fName;  // Name of the option 
    TString fDesc;  // Decription 
    TString fArg;   // Argument, if any
    Bool_t  fIsSet; // Whether the option has been set.
    TString fValue; // Value of the option. 
  };
  //====================================================================
  /**
   * Run a train setup
   * 
   */
  struct Runner 
  {
    /** 
     * Constructor 
     * 
     * @param train Train to run 
     * @param max   Maximum number of options
     */
    Runner(TrainSetup& train, UShort_t max=30)
      : fTrain(&train), fOptions(0), fN(0), fMax(max)
    {
      fOptions = new Option*[fMax];
      for (Int_t i = 0; i < fMax; i++) fOptions[i] = 0;
    }
    /** 
     * Add an option 
     * 
     * @param opt Option to add 
     */
    void Add(Option* opt)
    {
      if (fN >= fMax) {
	Warning("AddOption", "No room for option %s", opt->fName.Data());
	return;
      }
      fOptions[fN++] = opt;
    }
    /** 
     * Remove an option
     * 
     * @param name 
     */
    void Remove(const TString& name)
    {
      Option** ptr = fOptions;
      Option** tmp = 0;
      while (*ptr) { 
	if (name.EqualTo((*ptr)->fName)) {
	  tmp = ptr;
	  break;
	}
	ptr++;
      }
      if (!tmp) // nothing found, returning 
	return;
      
      ptr = tmp;
      delete *tmp;
      tmp++;
      fN--;
      while (*tmp) { 
	*ptr = *tmp;
	ptr++;
	tmp++;
      }
      *ptr = 0;
    }

    /** 
     * Parse option string 
     * 
     * @param options Option string. 
     * @param delim   Delimiters 
     * 
     * @return true on success. 
     */
    Bool_t Parse(const TString& options, const char* delim=",;")
    {
      TObjArray* a = options.Tokenize(delim);
      return Parse(*a);
    }
    /** 
     * Parse options 
     * 
     * @param options 
     * 
     * @return true on success
     */
    Bool_t Parse(TObjArray& options) 
    {
      TIter next(&options);
      TObjString* os = 0;
      while ((os = static_cast<TObjString*>(next()))) {
	TString s(os->String());
	// Info("Runner::Parse", "Processing option %s", s.Data());
	if (s.IsNull()) continue;

	Bool_t   ok  = false;
	Option** ptr = fOptions;
	while (*ptr && !ok) { 
	  Option* o = *ptr;
	  if (o->Process(s)) ok = true;
	  ptr++;
	}
	
	if (!ok) 
	  Warning("Parse", "Unknown option %s", s.Data());
      }
      return true;
    }
    /** 
     * Check if we asked for help 
     * 
     * 
     * @return 
     */
    Bool_t IsHelpAsked() const 
    {
      Option* help = FindOption("help");
      return (help && help->IsSet());
    }
    /** 
     * Print help
     * 
     * @param out Stream to write on  
     * @param prefix Prefix
     */
    void PrintHelp(std::ostream& out, const char* prefix="") const
    {
      Option** ptr = fOptions;
      while (*ptr) { 
	(*ptr)->PrintHelp(out, prefix);
	ptr++;
      }
    }
    /** 
     * Print the settings 
     * 
     * @param out Stream to write on. 
     */
    void PrintSettings(std::ostream& out) const
    {
      Option** ptr = fOptions;
      while (*ptr) { 
	(*ptr)->PrintSettings(out);
	ptr++;
      }
    }
    /** 
     * Find an option by name 
     * 
     * @param name Name of option to find
     * 
     * @return Pointer to option, or null
     */
    Option* FindOption(const TString& name) const
    {
      Option** ptr = fOptions;
      while (*ptr) { 
	if (name.EqualTo((*ptr)->fName)) return *ptr;
	ptr++;
      }
      return 0;
    }
    /** 
     * Initialize the train
     * 
     * @param options  Execution options 
     * 
     * @return true on success
     */
    Bool_t Init(const TString& options)
    {
      fTrain->MakeOptions(*this);
      if (!Parse(options)) return false;
      return true;
    }
    /** 
     * Run the train
     * 
     * @param runs     Run numbers 
     * @param nEvents  Number of events
     * @param asShell  Save set-up as shell script 
     * 
     * @return 
     */
    Bool_t Run(const TString& runs, Int_t nEvents, Bool_t asShell=false)
    {
      PrintSettings(std::cout);

      fTrain->SetOptions(*this);
      fTrain->SetRuns(runs);
      // fTrain->SaveSetup(*this, nEvents, asShell);
      
      fTrain->Run(nEvents, this, asShell);
      return true;
    }
      
    TrainSetup* fTrain;
    Option** fOptions;  // Our options 
    UShort_t fN;        // Current number of options 
    UShort_t fMax;      // Maximum number of options
  };
protected:
  //__________________________________________________________________
  /** 
   * @{ 
   * @name Options 
   */
  /** 
   * Class name of this train setup.  Sub-classes must define this. 
   * 
   * @return Class name of this setup 
   */
  virtual const char* ClassName() const = 0;
  //------------------------------------------------------------------
  /** 
   * Make the options for this train.  Sub-classes can overload this
   * to define new options, or append to the set of default option.
   */    
  virtual void MakeOptions(Runner& r)
  {
    r.Add(new Option("help",   "Show this help"));
    r.Add(new Option("par",    "Use PAR files (PROOF and Grid)"));
    r.Add(new Option("mc",     "Assume simulation input"));
    r.Add(new Option("debug",  "Execute in debugger"));
    r.Add(new Option("type",   "Type of train",    "AOD|ESD"));
    r.Add(new Option("mode",   "Execution mode",   "LOCAL|PROOF|GRID"));
    r.Add(new Option("oper",   "Operation mode",   "TEST|TERMINATE|FULL|INIT"));
    r.Add(new Option("date",   "Set date string",  "YYYY-MM-DD HH:MM:SS"));
    r.Add(new Option("cluster","PROOF cluster",         "HOST"));
    r.Add(new Option("dataSet","Data set (PROOF only)", "NAME"));
    r.Add(new Option("dataDir","Data directory",        "DIRECTORY"));
    r.Add(new Option("pattern","Data pattern (grid only)", "GLOB")); 
    r.Add(new Option("verb",   "Verbosity",             "NUMBER")); 
    r.Add(new Option("root",   "ROOT version (Grid)",   "TAG")); 
    r.Add(new Option("aliroot","AliROOT version (Grid)","TAG")); 
    r.Add(new Option("alien",  "AliEn API version (Grid)","TAG")); 
    r.Add(new Option("overwrite", "Allow overwrite"));
    r.Add(new Option("per-run", "Per run merge"));
  } 
  //------------------------------------------------------------------
  /** 
   * Set the option values on the train.  Sub-classes can overload 
   * this to set custom options on the train. 
   */
  virtual void SetOptions(Runner& r)
  {
    Option* debug	= r.FindOption("debug");
    Option* date	= r.FindOption("date");
    Option* cluster	= r.FindOption("cluster");
    Option* dataSet	= r.FindOption("dataSet");
    Option* dataDir	= r.FindOption("dataDir");
    Option* pattern	= r.FindOption("pattern");
    Option* par         = r.FindOption("par");
    Option* type        = r.FindOption("type");
    Option* mode        = r.FindOption("mode");
    Option* oper        = r.FindOption("oper");
    Option* mc          = r.FindOption("mc");
    Option* verb        = r.FindOption("verb");
    Option* root        = r.FindOption("root");
    Option* aliroot     = r.FindOption("aliroot");
    Option* alien       = r.FindOption("alien");
    Option* overwrite   = r.FindOption("overwrite");
    Option* run_merge   = r.FindOption("per-run");
    
    if (date && date->IsSet()) SetDateTime(date->AsString());
    if (cluster)               SetProofServer(cluster->AsString());
    if (dataSet)               SetDataSet(dataSet->AsString());
    if (dataDir)               SetDataDir(dataDir->AsString());
    if (pattern)               SetDataPattern(pattern->AsString());
    if (debug)                 SetUseGDB(debug->AsBool());
    if (type && type->IsSet()) SetType(type->AsString());
    if (mode && mode->IsSet()) SetMode(mode->AsString());
    if (oper && oper->IsSet()) SetOperation(oper->AsString());
    if (par)                   SetUsePar(par->AsBool());
    if (mc)                    SetMC(mc->AsBool());
    if (verb)                  SetVerbose(verb->AsInt());
    if (root)                  SetROOTVersion(root->AsString());
    if (aliroot)               SetAliROOTVersion(aliroot->AsString());
    if (alien)                 SetAliEnAPIVersion(alien->AsString());
    if (overwrite)             SetAllowOverwrite(overwrite->AsBool());
    if (run_merge)             SetPerRunMerge(run_merge->AsBool());
  }
  //------------------------------------------------------------------
  /** 
   * Set the option values on the train.  Sub-classes can overload 
   * this to set custom options on the train. 
   */
  virtual void SaveOptions(std::ostream& o, const char* str, Runner& r)
  {
    Option* debug	= r.FindOption("debug");
    Option* date	= r.FindOption("date");
    Option* cluster	= r.FindOption("cluster");
    Option* dataSet	= r.FindOption("dataSet");
    Option* dataDir	= r.FindOption("dataDir");
    Option* pattern	= r.FindOption("pattern");
    Option* par         = r.FindOption("par");
    Option* type        = r.FindOption("type");
    Option* mode        = r.FindOption("mode");
    Option* oper        = r.FindOption("oper");
    Option* mc          = r.FindOption("mc");
    Option* verb        = r.FindOption("verb");
    Option* root        = r.FindOption("root");
    Option* aliroot     = r.FindOption("aliroot");
    Option* alien       = r.FindOption("alien");
    Option* overwrite   = r.FindOption("overwrite");
    Option* run_merge   = r.FindOption("per-run");
    
    if (date) date->Save(o, str, 
			 Form("%04d-%02d-%02d %02d:%02d:00",
			      fDatime.GetYear(), 
			      fDatime.GetMonth(), 
			      fDatime.GetDay(), 
			      fDatime.GetHour(),
			      fDatime.GetMinute()));
    if (cluster)  cluster->Save(o, str, fProofServer);
    if (dataSet)  dataSet->Save(o, str, fDataSet);
    if (dataDir)  dataDir->Save(o, str, fDataDir);
    if (pattern)  pattern->Save(o, str, fDataPattern);
    if (debug)    debug->Save(o, str, fUseGDB);
    if (type)     type->Save(o, str, TypeString(fExecType));
    if (mode)     mode->Save(o, str, ModeString(fExecMode));
    if (oper)     oper->Save(o, str, OperString(fExecOper));
    if (par)      par->Save(o, str, fUsePar);
    if (mc)       mc->Save(o, str, fMC);
    if (verb)     verb->Save(o, str, fVerbose);
    if (root)     root->Save(o, str, fRootVersion);
    if (aliroot)  aliroot->Save(o, str, fAliRootVersion);
    if (alien)    alien->Save(o, str, fAliEnAPIVersion);
    if (overwrite)overwrite->Save(o, str, fAllowOverwrite);
    if (run_merge)run_merge->Save(o, str, fPerRunMerge);
  }
  /** 
   * Save the setup to file for later re-execution 
   * 
   * @param r         Runner object 
   * @param nEvents   Number of events 
   * @param asShell   If true, save as shell script - otherwise ROOT script
   */
  virtual void SaveSetup(Runner& r, Int_t nEvents, Bool_t asShell=false)
  {
    if (asShell) SaveSetupShell(r, nEvents); 
    /* else */   SaveSetupROOT(r, nEvents);
  }
  /** 
   * Save the setup to shell script for later re-execution 
   * 
   * @param r         Runner object 
   * @param nEvents   Number of events 
   */
  virtual void SaveSetupShell(Runner& r, Int_t nEvents)
  {
    std::ofstream o("rerun.sh");
    o << "#!/bin/bash\n\n"
      << "oper=$1\n"
      << "if test x$oper = x ; then oper=full ; fi \n\n"
      << "class=\"" << ClassName() << "\"\n"
      << "name=\"" << fName << "\"\n"
      << "nev=" << nEvents << "\n\n"
      << "opts=(--class=$class \\\n"
      << "  --name=$name \\\n"
      << "  --events=$nev \\" << std::endl;
    for (Int_t i = 0; i < fRunNumbers.GetSize(); i++) 
      o << "  --run=" << fRunNumbers.At(i) << " \\\n";
    SaveOptions(o, "--", r);
    o << "  --oper=$oper)\n\n"
      << "echo \"Running runTrain ${opts[@]}\"\n"
      << "runTrain \"${opts[@]}\"\n\n"
      << "#EOF\n" 
      << std::endl;
    o.close();
    gSystem->Exec("chmod a+x rerun.sh");
  }
  /** 
   * Save the setup to shell script for later re-execution 
   * 
   * @param r         Runner object 
   * @param nEvents   Number of events 
   */
  virtual void SaveSetupROOT(Runner& r, Int_t nEvents) 
  {
    std::ofstream o("rerun.C");
    o << "void rerun(bool terminate=false)\n"
      << "{\n" 
      << "  TString opts;" << std::endl;
    SaveOptions(o, "opts", r);
      
    o << "  if (terminate) opts.Append(\"mode=terminate;\");\n\n"
      << "  TString runs(\"";
    for (Int_t i = 0; i < fRunNumbers.GetSize(); i++) 
      o << (i == 0 ? "" : ", ") << fRunNumbers.At(i);
    o << "\");\n\n"
      << "  Int_t   nEvents = " << nEvents << ";\n\n"
      << "  gROOT->LoadMacro(\"$ALICE_ROOT/PWGLF/FORWARD/analysis2/trains/RunTrain.C\");\n"
      << "  RunTrain(\"" << ClassName() << "\",\"" 
      << fName << "\",opts,runs,nEvents);\n"
      << "}\n"
      << "// EOF" << std::endl;
    o.close();
  }
  /* @} */

  //__________________________________________________________________
  TString fName;             // Name of analysis
  TString fEscapedName;      // Name escaped for special chars
  TString fRootVersion;      // ROOT version to use 
  TString fAliRootVersion;   // AliROOT version to use 
  TString fAliEnAPIVersion;  // AliEn API version to use 
  TString fProofServer;      // Name of proof server
  TString fDataDir;          // Grid Input directory 
  TString fDataPattern;      // Data directory pattern
  TString fDataSet;          // Proof data set name 
  TString fXML;              // XML collection for local/proof mode
  Int_t   fNReplica;         // Storage replication
  Bool_t  fAllowOverwrite;   // Allow overwriting output dir
  Bool_t  fUseGDB;           // Wrap PROOF slaves in GDB 
  Int_t   fMaxSplit;         // Maximum number of files per split
  TArrayI fRunNumbers;       // List of run number 
  TList   fListOfPARs;       // List of PAR files to use 
  TList   fListOfSources;    // List of sources to upload and AcLIC
  TList   fListOfLibraries;  // List of libraries to load
  TList   fListOfExtras;     // List of extra files to upload
  TDatime fDatime;           // Date and time 
  EType   fExecType;         // Execution type (ESD, AOD)
  EMode   fExecMode;         // Execution mode (PROOF, local, Grid)
  EOper   fExecOper;         // Execution operation (full, terminate, ...)
  Bool_t  fUsePar;           // Wether to use PAR files 
  Bool_t  fMC;               // Whether to assume MC input 
  Bool_t  fPerRunMerge;      // Whether to merge per run or over-all
  Int_t   fVerbose;          // Verbosity level 
};
  
//____________________________________________________________________
//
// EOF
//
