/**
 * @defgroup pwglf_forward_trains Trains
 * 
 * Train specifications 
 *
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
# include <TStopwatch.h>
# include <AliAnalysisManager.h>
# include <AliVEventHandler.h>
# include <AliPhysicsSelection.h>
# include <AliPhysicsSelectionTask.h>
# include <AliCentralitySelectionTask.h>
# include <AliESDInputHandler.h>
# include <AliAODInputHandler.h>
# include <AliAODHandler.h>
# include <AliMCEventHandler.h>
# include <ctime>
#else 
struct Helper;
struct OptionList;
class TDatime;
class TUrl;
class TString;
class TStopwatch;
class AliVEventHandler;
class AliAnalysisManager;
class AliInputEventHandler;
#endif

//====================================================================
/** 
 * Generic set-up of an analysis train 
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
      fOptions(),
      fHelper(0)
  {
    fOptions.Add("help", "Show help", false);
    fOptions.Add("date", "YYYY-MM-DD HH:MM", "Set date", "now");
    fOptions.Add("bare-ps", "Use bare physics selection w/o task", false);
    fOptions.Add("verbose", "LEVEL", "Set verbosity level", 0);
    fOptions.Add("url", "URL", "Job location & input URL", "");
    fOptions.Add("overwrite", "Allow overwrite", false);
    fOptions.Add("events", "N", "Number of events to analyse", -1);
    fOptions.Add("type", "ESD|AOD|USER", "Input data stype", "");
    fOptions.Add("setup", "Only do the setup", false);
    fOptions.Add("branches", "Load only requested branches", false);
    fEscapedName = EscapeName(fName, "");
  }
  TrainSetup(const TrainSetup& o) 
    : fName(o.fName), 
      fEscapedName(o.fEscapedName), 
      fOptions(o.fOptions), 
      fHelper(o.fHelper)
  {}
  TrainSetup& operator=(const TrainSetup& o) 
  {
    if (&o == this) return *this;
    fName        = o.fName;
    fEscapedName = o.fEscapedName;
    fOptions     = o.fOptions;
    fHelper      = o.fHelper;
    return *this;
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

    fHelper = Helper::Create(url.Data(), verbose);
    if (!fHelper) { 
      Error("Init", "Failed to make the worker for URL %s", url.Data());
      return false;
    }

    // --- Check the type, if possible -------------------------------
    UShort_t type    = fHelper->InputType();
    Bool_t   mc      = fHelper->IsMC();
    if (fOptions.Has("type")) { 
      const TString& it = fOptions.Get("type");
      if      (it.EqualTo("ESD",TString::kIgnoreCase)) type = Helper::kESD;
      else if (it.EqualTo("AOD",TString::kIgnoreCase)) type = Helper::kAOD;
      else if (it.EqualTo("user",TString::kIgnoreCase)) 
	type = Helper::kUser;
    }

    // --- Rewrite the escpaed name ----------------------------------
    if (fOptions.Has("date")) { 
      TString date = fOptions.Get("date");
      fEscapedName = EscapeName(fName, date);
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
    AliAnalysisManager *mgr  = CreateAnalysisManager(fEscapedName);

    // In test mode, collect system information on every event 
    // if (oper == kTest)  mgr->SetNSysInfo(1); 
    if (verbose  >  0)      mgr->SetDebugLevel(verbose);
    mgr->SetAutoBranchLoading(!fOptions.Has("branches"));
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

    // --- Some information ------------------------------------------
    mgr->PrintStatus();
    if (fHelper->Mode() != Helper::kLocal) {
      TIter next(mgr->GetTasks());
      AliAnalysisTask* sub = 0;
      while ((sub = static_cast<AliAnalysisTask*>(next()))) { 
	sub->Print();
      }
    }
    return true;
  }
  void PrintTimer(TStopwatch& timer, const char* where)
  {
    timer.Stop();
    Double_t t = timer.RealTime();
    Int_t    h = Int_t(t / 3600); t -= h * 3600;
    Int_t    m = Int_t(t /   60); t -= m *   60;
    if (t < 0) t = 0;
    Info(where, "took %4d:%02d:%06.3f", h, m, t);
  }
    
  Bool_t Run()
  {
    TString cwd = gSystem->WorkingDirectory();
    Bool_t status = false;
    TStopwatch timer;
    timer.Start();
    try {
      if (!Init()) throw TString("Failed to intialize the train");
      PrintTimer(timer, "Initialization");
      timer.Continue();
      
      // Check if we're asked to only do the setup 
      if (fOptions.Has("setup")) {
	status = true;
	throw TString("Only did setup, no running");
      }

      // if (r) SaveSetup(*r, nEvents, asShell);
      
      Long64_t nEvents = fOptions.AsLong("events", -1);
      Long64_t ret     = fHelper->Run(nEvents);
      PrintTimer(timer, "Processing");
      timer.Continue();
      
      // Make sure we go back 
      gSystem->ChangeDirectory(cwd.Data());
      
      // Return. 
      if (ret < 0) throw TString("Analysis failed");

      status = true;
    }
    catch (TString& e) {
      if (status) Warning("Run", "%s", e.Data());
      else    	  Error("Run", "%s", e.Data());
    }
    if (fOptions.Has("date")) {
      TString escaped = EscapeName(fName, "");
      gSystem->Exec(Form("rm -f last_%s", escaped.Data()));
      gSystem->Exec(Form("ln -sf %s last_%s", 
			 fEscapedName.Data(), escaped.Data()));
    }
    PrintTimer(timer, "Finish");
    timer.Continue();
    return status;
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
   * @param o      Output stream
   * @param asProg If true, output as program options 
   * 
   * @return 
   */
  Bool_t Help(std::ostream& o=std::cout, bool asProg=false)
  {
    if (!fOptions.Has("help")) return true;

    if (!asProg) 
      o << "Usage: RunTrain(NAME, CLASS, OPTIONS)";
    
    o << "\n\nTrain Options:\n";
    fOptions.Help(o, asProg ? "  --" : "  ");
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
   * @param name   Train name 
   * @param cls    Class name 
   * @param opts   Comma seperated list of options
   * @param asProg Run as program 
   * @param spawn  Spawn ROOT shell after execution 
   * 
   * @return true on success
   */
  static Bool_t Main(const TString& name, const TString& cls, 
		     const TCollection* opts, 
		     Bool_t asProg=true,
		     Bool_t spawn=false)
  {
    Bool_t ret = false;
    try {
      if (cls.IsNull()) 
	throw TString("No class name specified");
      if (name.IsNull()) 
	throw TString("No train name specified");

      gROOT->ProcessLine("gSystem->RedirectOutput(\"build.log\",\"w\");");
      Int_t error = 0;
      Int_t r1 = gROOT->LoadMacro(Form("%s.C++g", cls.Data()), &error);
      gROOT->ProcessLine("gSystem->RedirectOutput(0);");
      if (r1 < 0 || error) 
	throw TString::Format("Failed to load setup %s: %d - see build.log", 
			      cls.Data(), error);

      // Make our object using the interpreter 
      TString create = TString::Format("new %s(\"%s\")", 
				       cls.Data(), name.Data());
      gROOT->ProcessLine("gSystem->RedirectOutput(\"build.log\",\"a\");");
      Long_t retP = gROOT->ProcessLine(create, &error);
      gROOT->ProcessLine("gSystem->RedirectOutput(0);");
      if (!retP || error) 
	throw TString::Format("Failed to make object of class %s "
			      "(see build.log): 0x%08lx/%d\n\t%s", 
			      cls.Data(), retP, error, create.Data());

      TrainSetup* train = reinterpret_cast<TrainSetup*>(retP);
    
      // Now parse the options 
      if (!train->Options().Parse(opts)) 
	throw TString("Failed to parse options");

      // Info("", "URL=%s", train->Options().Get("url").Data());

      // Check if we got a help request
      if (train->Options().Has("help")) { 
	train->Help(std::cout, asProg);
	ret = true;
	throw TString("");
      }

      // return train->Init();
      ret = train->Run();
    }
    catch (TString& e) { 
      if (!e.IsNull()) Error("Main", "%s", e.Data());
    }
    // Info("Main", "End of main loop (app=%p, asProg=%s, spawn=%s)",
    //	 gApplication, asProg ? "true" : "false", spawn ? "true" : "false");
    if (gApplication && asProg) {
      if (!spawn) {
	gSystem->Sleep(3);
	gApplication->Terminate(ret ? 0 : 1);
      }
    }
    return ret;
  }
protected:
  //__________________________________________________________________
  /** 
   * @{ 
   * @name Overloadable behaviour 
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
    mgr->AddStatisticsTask(AliVEvent::kAny);
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
    gROOT->Macro("AddTaskCentrality.C(true)");
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
  /** 
   * Add a task using a script and possibly some arguments 
   * 
   * @param macro Script to execute 
   * @param args  Optional arguments to the script 
   * 
   * @return Created task or null
   */
  virtual AliAnalysisTask* AddTask(const TString& macro, 
				   const TString& args)
  {
    TString p = gSystem->Which(gROOT->GetMacroPath(), macro.Data());
    if (p.IsNull()) { 
      Error("AddTask", "Macro %s not found", macro.Data());
      return 0;
    }
    TString cmd(p);
    if (!args.IsNull()) 
      cmd.Append(TString::Format("(%s)", args.Data()));
    
    Int_t err;
    Long_t ret = gROOT->Macro(cmd.Data(), &err, false);
    if (!ret) { 
      Error("AddTask", "Failed to execute %s", cmd.Data());
      return 0;
    }
    return reinterpret_cast<AliAnalysisTask*>(ret);
  }
  /** 
   * Add a task to the train with no arguments passed to the script 
   * 
   * @param macro The <b>AddTask</b> macro. 
   * 
   * @return The added task, if any
   */
  virtual AliAnalysisTask* AddTask(const TString& macro)
  {
    TString args;
    return AddTask(macro, args);
  }
  /** 
   * Add a single event analysis task to the train, passing the
   * specified arguments to the macro.
   * 
   * @param macro The <b>AddTask</b> macro 
   * @param args  Arguments to pass the macro 
   * 
   * @return The added task, if any 
   */
  virtual AliAnalysisTaskSE* AddSETask(const TString& macro, 
				       const TString& args)
  {
    return dynamic_cast<AliAnalysisTaskSE*>(AddTask(macro, args));
  }
  /** 
   * Add a single event task to the train with no arguments passed to
   * the script
   * 
   * @param macro The <b>AddTask</b> macro. 
   * 
   * @return The added task, if any
   */
  virtual AliAnalysisTaskSE* AddSETask(const TString& macro)
  {
    TString args;
    return AddSETask(macro, args);
  }

  /** 
   * Set the name of the train - should be name of the class.  Must be
   * overloaded.
   * 
   * @return Class name as a constant C string 
   */
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
   * @param name      Name to escape 
   * @param datimeStr Date and Time string 
   * 
   * @return escaped name 
   */  
  static TString EscapeName(const char* name, const TString& datimeStr)
  {
    TString escaped = name;
    char  c[] = { ' ', '/', '@', 0 };
    char* p   = c;
    while (*p) { 
      char tmp[] = { *p, '\0' };
      escaped.ReplaceAll(tmp, "_");
      p++;
    }
    if (!datimeStr.IsNull()) {
      TDatime datime;
      if (datimeStr.EqualTo("now", TString::kIgnoreCase)) 
	datime.Set();
      else {
	// Try various formats 
	struct tm t;
	const char* formats[] = { "%Ec", // Locale 
				  "%c", // Locale 
				  "%Ex EX", // Locale 
				  "%x %X", // Locale 
				  "%Y%m%d_%H%M", // YYYYMMDD_HHMM
				  "%F %R", // ISO standard, no seconds 
				  0 };
	const char** f = formats;
	Bool_t found = false;
	while (*f && !found) { 
	  // Reset needed fields 
	  t.tm_year  = 0;
	  t.tm_mon   = 0;
	  t.tm_mday  = 0;
	  t.tm_hour  = 0;
	  t.tm_min   = 0;
	  // Stop processing on first match 
	  if (strptime(datimeStr.Data(), *f, &t) != 0) found = true;
	  f++;
	}
	if (found) {
	  t.tm_mon += 1; // Return 0-based month
	  datime.Set(t.tm_year, t.tm_mon, t.tm_mday, t.tm_hour, t.tm_min, 0); 
	}
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
    // Info("SetupWorkingDirectory", "Made subdirectory %s, and cd'ed there", 
    //      nam.Data());
    return true;
  }
  /** 
   * Save the setup as a ROOT script and possibly also a shell script
   * 
   * @param asShellScript If true, also save as shell script
   */
  virtual void SaveSetup(Bool_t asShellScript)
  {
    OptionList tmp(fOptions);
    const OptionList* uopts = (fHelper ? &fHelper->Options() : 0);
    if (tmp.Find("overwrite")) tmp.Set("overwrite");
    if (tmp.Find("date") && fEscapedName.Length() > fName.Length()+1) {
      Int_t n = fName.Length()+1;
      tmp.Set("date", fEscapedName(n, fEscapedName.Length()-n));
    }
    if (asShellScript) 
      SaveSetupShell("rerun", ClassName(), fName, tmp, uopts);
    SaveSetupROOT("ReRun", ClassName(), fName, tmp, uopts);
    if (fHelper) fHelper->AuxSave(fEscapedName, asShellScript);
  }
  /** 
   * Save a setup as a shell script 
   * 
   * @param out   Output name of shell script 
   * @param cls   Class of the train 
   * @param name  Name of the train
   * @param opts  Option list
   * @param uopts Url options 
   */
  static void SaveSetupShell(const TString& out, const TString& cls,
			     const TString& name, const OptionList& opts,
			     const OptionList* uopts)
  {
    std::ofstream o(Form("%s.sh", out.Data()));
    o << "#!/bin/bash\n\n"
      << "class=\"" << cls << "\"\n"
      << "name=\"" << name << "\"\n\n"
      << "# Available options\n"
      << "# \n";
    opts.Help(o, "#    --");
    if (uopts) {
      o << "#\n"
	<< "# Available URI options\n"
	<< "# \n";
      uopts->Help(o, "#      ");
    }
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
   * @param uopts Url options 
   */
  static void SaveSetupROOT(const TString& out, const TString& cls,
			    const TString& name, const OptionList& opts,
			    const OptionList* uopts)
  {
    OptionList tmp(opts);
    tmp.Remove("url");

    std::ofstream o(Form("%s.C", out.Data()));
    o << "// Available options:\n"
      << "// \n";
    tmp.Help(o, "//     ");
    if (uopts) {
      o << "// \n"
	<< "// Available URI options\n";
      uopts->Help(o, "//      ");
    }
    o << "//\n"
      << "Bool_t " << out << "()\n" 
      << "{\n"
      << "  TString name(\"" << name << "\");\n"
      << "  TString cls(\"" << cls << "\");\n"
      << "  TUrl    uri(\"" << opts.Get("url") << "\");\n"
      << "  \n"
      << "  TString opts(";
    tmp.Store(o, "\"", ",\"\n               ", false);
    o << ");\n\n"
      << "  TString path(";
    TString     path(gROOT->GetMacroPath());
    TObjArray*  elements = path.Tokenize(":");
    TObjString* element = 0;
    TIter       next(elements);
    while ((element = static_cast<TObjString*>(next()))) {
      if (element->String().IsNull()) continue;
      o << "\n               \"" << element->GetName() << ":\"";
    }
    elements->Delete();
    o << ");\n"
      << "  path.Append(\"$ALICE_ROOT/PWGLF/FORWARD/trains\");\n"
      << "  gROOT->SetMacroPath(path);\n\n"
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
