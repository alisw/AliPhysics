/**
 * @defgroup pwglf_forward_trains Trains
 * 
 * Train specifications 
 *
 * @ingroup pwglf_forward
 */
/**
 * @file   TrainSetup.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 17:56:57 2012
 * 
 * @brief  Base classs for train specifications 
 * 
 * @ingroup pwglf_forward_trains
 */
#ifndef TRAINSETUP_C
#define TRAINSETUP_C
#ifndef __CINT__
# include "Helper.C"
# include "Option.C"
# include <TDatime.h>
# include <TUrl.h>
# include <TString.h>
# include <TApplication.h>
# include <AliAnalysisManager.h>
# include <AliVEventHandler.h>
# include <AliPhysicsSelection.h>
# include <AliPhysicsSelectionTask.h>
# include <AliCentralitySelectionTask.h>
# include <AliESDInputHandler.h>
# include <AliAODInputHandler.h>
# include <AliAODHandler.h>
# include <AliMCEventHandler.h>
#else 
struct Helper;
struct OptionList;
class TDatime;
class TUrl;
class TString;
class AliVEventHandler;
class AliAnalysisManager;
class AliInputEventHandler;
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
  /** 
   * Constructor 
   * 
   * @param name Name of the train
   */
  TrainSetup(const TString& name)
    : fName(name), 
      fEscapedName(name),
      fHelper(0)
  {
    fOptions.Add("help", "Show help");
    fOptions.Add("date", "YYYY-MM-DD HH:MM", "Set date", "now");
    fOptions.Add("mc", "Assume MC input");
    fOptions.Add("bare-ps", "Use bare physics selection w/o task");
    fOptions.Add("verbose", "LEVEL", "Set verbosity level", "0");
    fOptions.Add("url", "URL", "Job location & input URL", "");
    fOptions.Add("overwrite", "Allow overwrite");
    fOptions.Add("events", "N", "Number of events to analyse", "-1");
    fOptions.Add("type", "ESD|AOD|USER", "Input data stype", "");
    fEscapedName = EscapeName(fName, "");
  }
  /** 
   * Destructor
   */
  virtual ~TrainSetup() {}
  /* @} */
  //__________________________________________________________________
  /** 
   * @{ 
   * @name Execution 
   */
  /** 
   * Initialize 
   * 
   * @return true on success 
   */
  Bool_t Init()
  {
    // --- Create the helper -----------------------------------------
    TString  url     = fOptions.Get("url");
    Int_t    verbose = fOptions.AsInt("verbose");
    Bool_t   mc      = fOptions.AsBool("mc");

    fHelper = Helper::Create(url.Data(), verbose);
    if (!fHelper) { 
      Error("Init", "Failed to make the worker for URL %s", url.Data());
      return false;
    }

    UShort_t type    = fHelper->InputType();
    if (fOptions.Has("type")) { 
      const TString& it = fOptions.Get("type");
      if      (it.EqualTo("ESD",TString::kIgnoreCase)) type = Helper::kESD;
      else if (it.EqualTo("AOD",TString::kIgnoreCase)) type = Helper::kAOD;
      else if (it.EqualTo("user",TString::kIgnoreCase)) 
	type = Helper::kUser;
    }

    // --- Get current directory and set-up sub-directory ------------
    TString cwd = gSystem->WorkingDirectory();
    if (!SetupWorkingDirectory()) return false;

    // --- Do initial helper setup -----------------------------------
    if (!fHelper->PreSetup()) return false;

    // --- Load ROOT libraries ---------------------------------------
    if (!fHelper->LoadROOT()) return false;
    
    // --- Load AliROOT libraries ------------------------------------
    if (!fHelper->LoadAliROOT()) return false;

    // --- Create analysis manager -----------------------------------
    AliAnalysisManager *mgr  = CreateAnalysisManager(fName);

    // In test mode, collect system information on every event 
    // if (oper == kTest)  mgr->SetNSysInfo(1); 
    if (verbose  >  0)      mgr->SetDebugLevel(verbose);
    if (fHelper->Mode() == Helper::kLocal) 
      mgr->SetUseProgressBar(kTRUE, 100);
   
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
    if (type == Helper::kESD) CreatePhysicsSelection(mc, mgr);
    
    // --- Create centrality task ------------------------------------
    CreateCentralitySelection(mc, mgr);

    // --- Create tasks ----------------------------------------------
    CreateTasks(mgr);

    // --- Post set-up initialization of helper ----------------------
    if (!fHelper->PostSetup()) return false;

    // --- Set debug level on defined tasks --------------------------
    if (verbose > 0) {
      TIter next(mgr->GetTasks());
      AliAnalysisTask* sub = 0;
      while ((sub = static_cast<AliAnalysisTask*>(next()))) { 
	AliAnalysisTaskSE* se = dynamic_cast<AliAnalysisTaskSE*>(sub);
	if (!se) continue;
	se->SetDebugLevel(verbose);
      }
    }

    // --- Print this setup ------------------------------------------
    Print();

    // --- Initialise the train --------------------------------------
    if (!mgr->InitAnalysis())  {
      gSystem->ChangeDirectory(cwd.Data());
      Error("Init","Failed to initialise train");
      return false;
    }

    // --- Enable progress bar ---------------------------------------
    if (fHelper->Mode() != Helper::kGrid) 
      mgr->SetUseProgressBar(true, 100);

    // --- Save setup to disk ----------------------------------------
    SaveSetup(true);

    // Some information
    mgr->PrintStatus();

    return true;
  }
  Bool_t Run(Bool_t doExit=false)
  {
    TString cwd = gSystem->WorkingDirectory();
    Bool_t status = false;
    try {
      if (!Init()) throw TString("Failed to intialize the train");

      // if (r) SaveSetup(*r, nEvents, asShell);
      
      Long64_t nEvents = fOptions.AsLong("events", -1);
      Long64_t ret     = fHelper->Run(nEvents);
      
      // Make sure we go back 
      gSystem->ChangeDirectory(cwd.Data());
      
      // Return. 
      if (ret < 0) throw TString("Analysis failed");

      status = true;
    }
    catch (TString& e) {
      Error("Main", e);
      status = false;
    }
    if (gApplication && doExit) {
      gSystem->Sleep(3);
      gApplication->Terminate(status ? 0 : 1);
    }
    return true;
  }
  /** 
   * Get the options 
   * 
   * 
   * @return Reference ot the options 
   */
  OptionList& Options() { return fOptions; }
  /** 
   * Print information to standard output 
   * 
   */
  void Print(Option_t* ="") const
  {
    std::cout << "Train: " << fName << " (" << fEscapedName << ")" 
	      << std::endl;
    fOptions.Show(std::cout);
    if (fHelper) fHelper->Print();
  }
  /** 
   * Show the help 
   * 
   * @param o 
   * 
   * @return 
   */
  Bool_t Help(std::ostream& o=std::cout, bool asProg=false)
  {
    if (!fOptions.Has("help")) return true;

    if (asProg) 
      o << "Usage: runTrain --name=NAME --class=CLASS [OPTIONS]";
    else 
      o << "Usage: RunTrain(NAME, CLASS, OPTIONS)";
    
    o << "\n\nOptions:\n\n";
    if (asProg) {
      OptionList tmp(fOptions);
      tmp.Add("name", "STRING", "Name of train", fName);
      tmp.Add("class", "NAME", "Name of setup class", "");
      tmp.Help(o,"  --");
    }
    else 
      fOptions.Help(o, "  ");
    
    o << "\n";

    if (!fHelper && fOptions.Has("url")) {
      TString url = fOptions.Get("url");
      fHelper = Helper::Create(url.Data());
    }
    if (fHelper) { 
      o << fHelper->Desc() << " URL form:\n\n"
	<< "    " << fHelper->UrlHelp() << "\n\n"
	<< "Options:\n";
      fHelper->Options().Help(o, "    ");
      o << "\n";
    }
    else { 
      o << "Possible URL forms:\n\n";
      Helper::ShowUrlHelp("LocalHelper");
      Helper::ShowUrlHelp("ProofHelper");
      Helper::ShowUrlHelp("LiteHelper");
      Helper::ShowUrlHelp("AAFHelper");
      Helper::ShowUrlHelp("AAFPluginHelper");
      Helper::ShowUrlHelp("GridHelper");
      o << "\n";
    }
    return false;
  }
  /** 
   * Run train.  This will AcLic compile the setup script, create
   * an object of that type with the given name, and then pass the 
   * options to it.  Then, it will run the setup.
   * 
   * @param name  Train name 
   * @param cls   Class name 
   * @param opts  Comma seperated list of options
   * 
   * @return true on success
   */
  static Bool_t Main(const TString& name, const TString& cls, 
		     const TCollection* opts, 
		     Bool_t asProg=true)
  {
    if (cls.IsNull()) { 
      Error("Main", "No class name specified");
      return false;
    }
    if (name.IsNull()) { 
      Error("Main", "No train name specified");
      return false;
    }
    Int_t error = 0;
    Int_t r1 = gROOT->LoadMacro(Form("%s.C++g", cls.Data()), &error);
    if (r1 < 0 || error) { 
      Error("Main", "Failed to load setup %s: %d", cls.Data(), error);
      return false;
    }

    // Make our object using the interpreter 
    TString create = TString::Format("new %s(\"%s\")", 
				   cls.Data(), name.Data());
    gROOT->ProcessLine("gSystem->RedirectOutput(\"/dev/null\",\"w\");");
    Long_t ret = gROOT->ProcessLine(create, &error);
    gROOT->ProcessLine("gSystem->RedirectOutput(0);");
    if (!ret || error) { 
      Error("Main", "Failed to make object of class %s: 0x%08lx/%d\n\t%s", 
	    cls.Data(), ret, error, create.Data());
      return false;
    }
    TrainSetup* train = reinterpret_cast<TrainSetup*>(ret);
    
    // Now parse the options 
    if (!train->Options().Parse(opts)) { 
      Error("Main", "Failed to parse options");
      return false;
    }

    // Check if we got a help request
    if (train->Options().Has("help")) { 
      train->Help(std::cout, asProg);
      return true;
    }
    // return train->Init();
    return train->Run(asProg);
  }
protected:
  //__________________________________________________________________
  /** 
   * @{ 
   * @Name Overloadable behaviour 
   */
  //------------------------------------------------------------------
  /** 
   * Create the analysis manager 
   * 
   * @param name Name of the analysis 
   * 
   * @return Created analysis manager 
   */
  virtual AliAnalysisManager* CreateAnalysisManager(const char* name)
  {
    return new AliAnalysisManager(name,"Analysis Train");
  }
  //------------------------------------------------------------------
  /** 
   * Create input handler 
   * 
   * @param type 
   * 
   * @return 
   */
  virtual AliVEventHandler* CreateInputHandler(UShort_t type)
  {
    switch (type) {
    case Helper::kESD:  return new AliESDInputHandler(); 
    case Helper::kAOD:  return new AliAODInputHandler(); 
    case Helper::kUser: return 0;
    }
    return 0;
  }
  //------------------------------------------------------------------
  /** 
   * Create MC input handler 
   * 
   * @param type  Run type (ESD or AOD)
   * @param mc    Assume monte-carlo input 
   * 
   * @return 
   */
  virtual AliVEventHandler* CreateMCHandler(UShort_t /*type*/, bool mc)
  {
    if (!mc) return 0;
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
  virtual AliVEventHandler* CreateOutputHandler(UShort_t type)
  {
    AliAODHandler* ret = new AliAODHandler();
    switch (type) { 
    case Helper::kESD: 
      ret->SetOutputFileName("AliAOD.root");
      break;
    case Helper::kAOD: 
      ret->SetOutputFileName("AliAOD.pass2.root");
      break;
    case Helper::kUser: 
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
  virtual void CreatePhysicsSelection(Bool_t mc, AliAnalysisManager* mgr)
  {
    if (fOptions.Has("bare-ps")) {
      AliInputEventHandler* input = 
	dynamic_cast<AliInputEventHandler*> (mgr->GetInputEventHandler());
      if (!input) return;

      AliPhysicsSelection* ps = new AliPhysicsSelection();
      if (mc) ps->SetAnalyzeMC();

      input->SetEventSelection(ps);

      return;
    }
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
   * @param mgr  Manager
   */
  virtual void CreateTasks(AliAnalysisManager* mgr)=0;
  virtual const Char_t* ClassName() const = 0;
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
  static TString EscapeName(const char* name, const TString& datimeStr)
  {
    TString escaped = name;
    char  c[] = { ' ', '/', '@', 0 };
    char* p   = c;
    while (*p) { 
      escaped.ReplaceAll(Form("%c", *p), "_");
      p++;
    }
    if (!datimeStr.IsNull()) {
      TDatime datime;
      if (datimeStr.EqualTo("now", TString::kIgnoreCase)) 
	datime.Set();
      else 
	datime.Set(datimeStr.Data());
      if (datime.GetYear() <= 1995 ||
	  datime.GetMonth() == 0 || 
	  datime.GetDay() == 0) return escaped;
      escaped.Append(Form("_%04d%02d%02d_%02d%02d", 
			  datime.GetYear(), 
			  datime.GetMonth(), 
			  datime.GetDay(), 
			  datime.GetHour(), 
			  datime.GetMinute()));
    }
    return escaped;
  }    
  /** 
   * Make our working directory if so requested 
   * 
   * @return true on success
   */
  Bool_t SetupWorkingDirectory()
  {
    // Get the name of the target directory 
    TString& nam = fEscapedName;

    // Check if the directory exists already 
    Bool_t exists = gSystem->AccessPathName(nam.Data()) == 0;
    if (fHelper->Operation() == Helper::kTerminate && !exists) {
      Error("SetupWorkingDirectory", "File/directory %s does not exists", 
	    nam.Data());
      return false;
    }

    Bool_t overwrite = fOptions.Has("overwrite");
    // If we're not allowed to overwrite, then complain
    if (!overwrite && exists) {
      Error("SetupWorkingDirectory", "File/directory %s already exists", 
	    nam.Data());
      return false;
    }

    // Make the target directory if it doesn't exists 
    if (!exists) {
      if (gSystem->MakeDirectory(nam.Data())) {
	Error("SetupWorkingDirectory", "Failed to make directory '%s'", 
	      nam.Data());
	return false;
      }
    }
      
    // Change directory to target directory 
    if (!gSystem->ChangeDirectory(nam.Data())) { 
      Error("SetupWorkingDirectory", "Failed to change directory to %s", 
	    nam.Data());
      return false;
    }
    Info("SetupWorkingDirectory", "Made subdirectory %s, and cd'ed there", 
	 nam.Data());
    return true;
  }
  /** 
   * Save the setup as a ROOT script and possibly also a shell script
   * 
   * @param asShellScript If true, also save as shell script
   */
  virtual void SaveSetup(Bool_t asShellScript)
  {
    if (asShellScript) 
      SaveSetupShell("rerun", ClassName(), fName, fOptions);
    SaveSetupROOT("ReRun", ClassName(), fName, fOptions);
  }
  /** 
   * Save a setup as a shell script 
   * 
   * @param out   Output name of shell script 
   * @param cls   Class of the train 
   * @param name  Name of the train
   * @param opts  Option list
   */
  static void SaveSetupShell(const TString& out, const TString& cls,
			      const TString& name, const OptionList& opts)
  {
    std::ofstream o(Form("%s.sh", out.Data()));
    o << "#!/bin/bash\n\n"
      << "class=\"" << cls << "\"\n"
      << "name=\"" << name << "\"\n\n"
      << "# Available options\n"
      << "# \n";
    opts.Help(o, "#    --");
    o << "#\n"
      << "opts=(--class=$class \\\n"
      << "  --name=$name";
    opts.Store(o, " \\\n  --", "", true);
    o << ")\n\n"
      << "echo \"Running runTrain ${opts[@]} $@\"\n"
      << "runTrain \"${opts[@]}\" $@\n\n"
      << "# EOF" << std::endl;
    o.close();
    gSystem->Exec(Form("chmod a+x %s.sh", out.Data()));
  }
  /** 
   * Save a setup as a ROOT script 
   * 
   * @param out   Output name of shell script 
   * @param cls   Class of the train 
   * @param name  Name of the train
   * @param opts  Option list
   */
  static void SaveSetupROOT(const TString& out, const TString& cls,
			    const TString& name, const OptionList& opts)
  {
    OptionList tmp(opts);
    tmp.Remove("url");
    std::ofstream o(Form("%s.C", out.Data()));
    o << "/* Available options:\n"
      << " *\n";
    tmp.Help(o, " *    ");
    o << " */\n"
      << "Bool_t " << out << "()\n"
      << "{\n"
      << "  TString name(\"" << name << "\");\n"
      << "  TString cls(\"" << cls << "\");\n"
      << "  TString uri(\"" << opts.Get("url") << "\");\n"
      << "  TString opts(";
    tmp.Store(o, "\"", ",\"\n               ", false);
    o << ");\n\n"
      << "  gROOT->LoadMacro(\"RunTrain.C\");\n\n"
      << "  return RunTrain(name, cls, uri, opts);\n"
      << "}\n" 
      << "//\n"
      << "// EOF\n" 
      << "//" << std::endl;
    o.close();
  }    
  /* @} */
  TString      fName;
  TString      fEscapedName;
  OptionList   fOptions;
  Helper*      fHelper;
};
#endif
