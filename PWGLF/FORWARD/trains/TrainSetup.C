/**
 * @defgroup pwglf_forward_trains Trains.
 * 
 * Train specifications. 
 * See also @ref train_setup_doc
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
# include "Railway.C"
# include "Option.C"
# include <TDatime.h>
# include <TUrl.h>
# include <TString.h>
# include <TApplication.h>
# include <TStopwatch.h>
# include <AliAnalysisManager.h>
# include <AliVEvent.h>
# include <AliVEventHandler.h>
# include <AliPhysicsSelection.h>
# include <AliPhysicsSelectionTask.h>
# include <AliCentralitySelectionTask.h>
# include <AliESDInputHandler.h>
# include <AliESDInputHandlerRP.h>
# include <AliAODInputHandler.h>
# include <AliAODHandler.h>
# include <AliMCEventHandler.h>
# include <ctime>
#else 
struct Railway;
struct OptionList;
class TDatime;
class TUrl;
class TString;
class TStopwatch;
class AliVEventHandler;
class AliAnalysisManager;
class AliInputEventHandler;
#endif
class AliAnalysisTask;

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
   * Version of this framework 
   */
  enum {
    kVersion = 2
  }; 
  /** 
   * Constructor 
   * 
   * @param name Name of the train
   */
  TrainSetup(const TString& name)
    : fName(name), 
      fEscapedName(name),
      fDatimeString(""),
      fOptions(),
      fRailway(0)
  {
    fOptions.Add("help",   "Show help",                                  false);
    fOptions.Add("date",   "YYYY-MM-DD HH:MM", "Set date",               "now");
    fOptions.Add("ps",     "MODE",             "Physics selection mode", "");
    fOptions.Add("verbose","LEVEL",            "Set verbosity level",    0);
    fOptions.Add("url",    "URL",              "Job location & input URL","");
    fOptions.Add("overwrite","Allow overwrite",                          false);
    fOptions.Add("events", "N",             "Number of events to analyse",-1);
    fOptions.Add("type",   "ESD|AOD|USER",     "Input data stype",       "");
    fOptions.Add("setup",  "Only do the setup",                          false);
    fOptions.Add("branches","Load only requested branches",              false);
    fOptions.Add("version","Print version and exit",                     false);
    fOptions.Add("tender", "WHICH",            "Specify tender supplies","");
    fOptions.Add("ocdb",   "(TENDER_SNAPSHOT)","Enable OCDB",            "");
    fOptions.Add("friends","(AOD_FRIENDS)","Enable friends (list of files)","");
    fOptions.Add("cent-oadb","PERIOD","Alternative OADB for centrality","");
    fOptions.Add("no-link","Do not make symlink to output",              false);
    fDatimeString = "";
    fEscapedName  = EscapeName(fName, fDatimeString);
  }
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  TrainSetup(const TrainSetup& o) 
    : fName(o.fName), 
      fEscapedName(o.fEscapedName), 
      fDatimeString(o.fDatimeString),
      fOptions(o.fOptions), 
      fRailway(o.fRailway)
  {}
  /** 
   * Assignment operator
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this object
   */
  TrainSetup& operator=(const TrainSetup& o) 
  {
    if (&o == this) return *this;
    fName         = o.fName;
    fEscapedName  = o.fEscapedName;
    fDatimeString = o.fDatimeString;
    fOptions      = o.fOptions;
    fRailway      = o.fRailway;
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
    // --- Print the version number ----------------------------------
    Info("Init", "Running with TrainSetup version %d", kVersion);
    
    // --- Create the helper -----------------------------------------
    TString  url     = fOptions.Get("url");
    Int_t    verbose = fOptions.AsInt("verbose");

    fRailway = Railway::Create(url.Data(), verbose);
    if (!fRailway) { 
      Error("Init", "Failed to make the worker for URL %s", url.Data());
      return false;
    }

    // --- Check the type, if possible -------------------------------
    UShort_t type    = fRailway->InputType();
    Bool_t   mc      = fRailway->IsMC();
    if (fOptions.Has("type")) { 
      const TString& it = fOptions.Get("type");
      if      (it.EqualTo("ESD",TString::kIgnoreCase)) type = Railway::kESD;
      else if (it.EqualTo("AOD",TString::kIgnoreCase)) type = Railway::kAOD;
      else if (it.EqualTo("user",TString::kIgnoreCase)) 
	type = Railway::kUser;
    }

    // --- Rewrite the escpaed name ----------------------------------
    if (fOptions.Has("date")) { 
      fDatimeString = fOptions.Get("date");
      fEscapedName  = EscapeName(fName, fDatimeString);
    }
    
    // --- Get current directory and set-up sub-directory ------------
    TString cwd = gSystem->WorkingDirectory();
    if (!SetupWorkingDirectory()) return false;

    // --- Do initial helper setup -----------------------------------
    if (!fRailway->PreSetup()) return false;

    // --- Load ROOT libraries ---------------------------------------
    if (!fRailway->LoadROOT()) return false;
    
    // --- Load AliROOT libraries ------------------------------------
    if (!fRailway->LoadAliROOT()) return false;

    // --- Load AliROOT libraries ------------------------------------
    if (!fRailway->LoadAliPhysics()) return false;

    // --- Create analysis manager -----------------------------------
    AliAnalysisManager *mgr  = CreateAnalysisManager(fEscapedName);

    // In test mode, collect system information on every event 
    // if (oper == kTest)  mgr->SetNSysInfo(1); 
    if (verbose  >  0)      mgr->SetDebugLevel(verbose);
    mgr->SetAutoBranchLoading(!fOptions.Has("branches"));
    if (fRailway->Mode() == Railway::kLocal) 
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
    gROOT->SetMacroPath(Form("%s:%s:$ALICE_ROOT/ANALYSIS/macros:"
			     "$ALICE_PHYSICS/OADB/macros",
			     cwd.Data(), gROOT->GetMacroPath()));

    // --- Tender/OCDB -----------------------------------------------
    if (type == Railway::kESD) {
      TString supplies = fOptions.Get("tender");
      if (!supplies.IsNull()) {
	AddTender(supplies);
      }
      else if (fOptions.AsBool("ocdb")) {
	AddOCDBConnect();
      }
    }
    
    // --- Physics selction - only for ESD ---------------------------
    if (type == Railway::kESD) CreatePhysicsSelection(mc, mgr);
    
    // --- Create centrality task ------------------------------------
    CreateCentralitySelection(mc);

    // --- Create tasks ----------------------------------------------
    CreateTasks(mgr);

    // --- Create monitor objects ------------------------------------
    CreateMonitors();

    // --- Post set-up initialization of helper ----------------------
    if (!fRailway->PostSetup()) return false;

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
    if (fRailway->Mode() != Railway::kGrid) 
      mgr->SetUseProgressBar(true, 100);

    // --- Save setup to disk ----------------------------------------
    SaveSetup(true);

    // --- Some information ------------------------------------------
    mgr->PrintStatus();
    if (fRailway->Mode() != Railway::kLocal) {
      TIter next(mgr->GetTasks());
      AliAnalysisTask* sub = 0;
      while ((sub = static_cast<AliAnalysisTask*>(next()))) { 
	sub->Print();
      }
    }
    return true;
  }
  /** 
   * Print timer information
   * 
   * @param timer The timer
   * @param where Where this was called from 
   */
  void PrintTimer(TStopwatch& timer, const char* where)
  {
    timer.Stop();
    Double_t t = timer.RealTime();
    Int_t    h = Int_t(t / 3600); t -= h * 3600;
    Int_t    m = Int_t(t /   60); t -= m *   60;
    if (t < 0) t = 0;
    Info(where, "took %4d:%02d:%06.3f", h, m, t);
  }
  /** 
   * Run this train 
   * 
   * @return true on success 
   */    
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
      Long64_t ret     = fRailway->Run(nEvents);
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
    if (fOptions.Has("date") &&
	!fOptions.Get("date").EqualTo("none") &&
	!fOptions.Has("no-link")) {
      TString tmp     = "";
      TString escaped = EscapeName(fName, tmp);
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
    if (fRailway) fRailway->Print();
  }
  /** 
   * Show the version number 
   *
   * @param o      Output stream 
   * 
   * @return true help wasn't requested
   */
  Bool_t Version(std::ostream& o=std::cout)
  {
    if (!fOptions.Has("version")) return true;

    o << "TrainSetup version " << kVersion << std::endl;
    return false;
  }
  /** 
   * Show the help 
   * 
   * @param o      Output stream
   * @param asProg If true, output as program options 
   * 
   * @return true help wasn't requested
   */
  Bool_t Help(std::ostream& o=std::cout, bool asProg=false)
  {
    if (!fOptions.Has("help")) return true;
    
    if (!asProg) 
      o << "Usage: RunTrain(NAME, CLASS, OPTIONS)";
    
    o << "\n\nTrain Options:\n";
    fOptions.Help(o, asProg ? "  --" : "  ");
    o << "\n";

    if (!fRailway && fOptions.Has("url")) {
      TString url = fOptions.Get("url");
      fRailway = Railway::Create(url.Data());
    }
    if (fRailway) { 
      o << fRailway->Desc() << " URL form:\n\n"
	<< "    " << fRailway->UrlHelp() << "\n\n"
	<< "Options:\n";
      fRailway->Options().Help(o, "    ");
      o << "\n";
    }
    else { 
      o << "Possible URL forms:\n\n";
      Railway::ShowUrlHelp("LocalRailway");
      Railway::ShowUrlHelp("ProofRailway");
      Railway::ShowUrlHelp("LiteRailway");
      Railway::ShowUrlHelp("VAFRailway");
      Railway::ShowUrlHelp("AAFRailway");
      Railway::ShowUrlHelp("AAFPluginRailway");
      Railway::ShowUrlHelp("GridRailway");
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
      TString mkLib = gSystem->GetMakeSharedLib();
      mkLib.ReplaceAll("-std=c++14", "-std=c++98");
      // mkLib.Append("-ffpe-trap=invalid,zero,overflow,underflow,inexact");
      gSystem->SetMakeSharedLib(mkLib);
  
      gROOT->ProcessLine("gSystem->RedirectOutput(\"build.log\",\"w\");");
      Int_t error = 0;
      Int_t r1 = gROOT->LoadMacro(Form("%s.C+g", cls.Data()), &error);
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
      if (!train->Help(std::cout, asProg)) {
	ret = true;
	throw TString("");
      }
      // Check if we got a version request
      if (!train->Version(std::cout)) {
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
  /** 
   * Create input handler 
   * 
   * @param type         Type of input (ESD, AOD, or user)
   * @param esdRecPoints if type is ESD and this is true, create input
   * handler for rec-points (clusters).
   * 
   * @return 
   */
  virtual AliVEventHandler* CreateInputHandler(UShort_t type,
					       Bool_t   esdRecPoints=false)
  {
    Info("CreateInputHandler", "Making handler for %d (%d)",
	 type, esdRecPoints);
    AliVEventHandler* ret = 0;
    switch (type) {
    case Railway::kESD:  {
      AliESDInputHandler* input = 0;
      if (!esdRecPoints) input = new AliESDInputHandler();
      else {
	Info("CreateInputHandler", "Special handler for rec-points");
	AliESDInputHandlerRP* esd = new AliESDInputHandlerRP();
	// esd->ReadFromDirectory();
	input = esd;	
      }
      input->SetReadFriends(fOptions.AsBool("friends"));
      ret = input;
    }
      break;
    case Railway::kAOD:  {
      AliAODInputHandler* input = new AliAODInputHandler();
      TString    fr  = fOptions.Get("friends");
      TObjArray* afr = fr.Tokenize(",+:");
      TObject*   ofr = 0;
      TIter      nfr(afr);
      while (ofr = nfr()) input->AddFriend(const_cast<char*>(ofr->GetName()));
      afr->Delete();
      ret = input;
    }
      break;
    case Railway::kUser: return 0;
    }
    // Info("CreateInput", "Returning input handler %p", ret);
    return ret;
  }
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
    case Railway::kESD: 
      ret->SetOutputFileName("AliAOD.root");
      break;
    case Railway::kAOD: 
      ret->SetOutputFileName("AliAOD.pass2.root");
      break;
    case Railway::kUser: 
      break;
    }
    
    return ret;
  }
  /** 
   * Create physics selection, and add to manager
   *
   * The physis selectuion used can be specified in the option @c ps. 
   * 
   * - NONE            Do not set a physics selection 
   * - BARE            Set physics selection on input handler directly 
   * - CUSTOM[=macro]  Read custom physics selection, @c macro.  
   *                   If not macro is given, then assume it to be 
   *                   @c CustomPS.C 
   * - ALL             Do not ignore background triggers. 
   *
   * The macro passed to CUSTOM must accept a single argument of a
   * pointer to an AliPhysicsSelection object.  The macro should set
   * the collision and background candidates, as well as the filling
   * scheme, and trigger analysis.
   *
   * @code 
   void CustomPS(AliPhysicsSelection* ps)
   {
     // --- Defaults filling scheme ----------------------------------
     AliOADBFillingScheme * cFS = new AliOADBFillingScheme("Default");
     cFS->SetFillingSchemeName("Default");
     cFS->SetBXIDs("B",  "");
     cFS->SetBXIDs("A",  "");
     cFS->SetBXIDs("AC", "");
     cFS->SetBXIDs("ACE","");
     cFS->SetBXIDs("C",  "");
     cFS->SetBXIDs("E",  "");

     // --- Define hardware triggers, and how to replay them offline -
     // This redefines Min.Bias. by an or of VZERO-OR and SPD-FASTOR 
     AliOADBPhysicsSelection * cPS = new AliOADBPhysicsSelection("cPS");
     Int_t id = 0;
     // B (bunch-crossing) trigger - or of CINT5 (V0-OR) 
     // and SPD COSMB (FASTOR)
     cPS->AddCollisionTriggerClass(AliVEvent::kMB,
                                   "+CINT5-B-NOPF-ALLNOTRD,"
				   "+C0SMB-B-NOPF-ALLNOTRD","B",    id);
     // Control triggers 
     cPS->AddBGTriggerClass       (AliVEvent::kMB,
                                   "+CINT5-A-NOPF-ALLNOTRD,"
                                   "+C0SMB-A-NOPF-ALLNOTRD","A",    id);
     cPS->AddBGTriggerClass       (AliVEvent::kMB,
                                   "+CINT5-C-NOPF-ALLNOTRD,"
                                   "+C0SMB-C-NOPF-ALLNOTRD","C",    id);
     cPS->AddBGTriggerClass       (AliVEvent::kMB,
                                   "+CINT5-E-NOPF-ALLNOTRD,"
                                   "+C0SMB-E-NOPF-ALLNOTRD","E",    id);
     cPS->AddBGTriggerClass       (AliVEvent::kMB,
                                   "+CINT5-ACE-NOPF-ALLNOTRD,"
                                   "+C0SMB-ACE-NOPF-ALLNOTRD","ACE",id);
     // How to replay hardware trigger 
     cPS->SetHardwareTrigger      (id,"SPDGFO >= 1 || V0A || V0C");
     // Additional off-line trigger 
     cPS->SetOfflineTrigger       (id,"(SPDGFO >= 1 || V0A || V0C) "
                                  "&& !V0ABG && !V0CBG && !TPCLaserWarmUp");
  
     // --- Trigger analysis defaults --------------------------------
     AliOADBTriggerAnalysis * cTA = new AliOADBTriggerAnalysis("Default");
     cTA->SetZDCCorrParameters(0.5, 0, 4*0.7, 4*0.7);
   }
   @endcode 
   *
   * @c AliPhysicsSelection::AddCollisionTriggerClass sets a hardware
   * trigger (combination) to be a collision of a given kind (e.g.,
   * kMB). First argument is the collision kind (see AliVEvent), second is the 
   *
   * @c AliPhysicsSelection::AddBGTrigerClass sets a hardware trigger
   * (combination) to be a background trigger of the quivilent
   * collisions trigger.
   *
   * @c AliPhysicsSelection::SetHardwareTrigger sets how to replay the
   * hardware trigger off-line.
   *
   * @c AliPhysicsSelection::SetOfflineTrigger sets the off-line
   * conditions to meet.
   * 
   * @param mc Whether this is for MC 
   * @param mgr Manager
   */
  virtual void CreatePhysicsSelection(Bool_t mc, AliAnalysisManager* mgr)
  {
    TString opt = fOptions.Get("ps");
    // opt.ToUpper();

    if (opt.EqualTo("NONE",TString::kIgnoreCase)) return;

    // Possibly load as PAR
    LoadOADB();

    AliPhysicsSelection* ps = 0;
    AliInputEventHandler* input = 
      dynamic_cast<AliInputEventHandler*> (mgr->GetInputEventHandler());
    if (!input) return;

    // --- Create object, either directory or via task ---------------
    if (opt.Contains("BARE",TString::kIgnoreCase)) {
      // --- Just create object and set on input handler -------------
      ps = new AliPhysicsSelection();
      if (mc) ps->SetAnalyzeMC();

      input->SetEventSelection(ps);
    }
    else {
      // --- Create and add the task ---------------------------------
      CoupleSECar("AddTaskPhysicsSelection.C", Form("%d", mc), AliVEvent::kAny);
      mgr->RegisterExtraFile("event_stat.root");
      // mgr->AddStatisticsTask(AliVEvent::kAny);
      
      // --- Retrive object from input handler -----------------------
      ps = dynamic_cast<AliPhysicsSelection*>(input->GetEventSelection());
    }
    if (opt.Contains("CUSTOM",TString::kIgnoreCase)) {
      // --- Load custom trigger definitions -------------------------
      TString macro("CustomPS.C");
      Int_t eq = opt.Index("custom=",7,0,TString::kIgnoreCase);
      if (eq != kNPOS) {
	Int_t end = opt.Index(".C",2,eq+7,TString::kIgnoreCase);
	if (end != kNPOS) {
	  macro = opt(eq+7,end+2-eq-7);
	}
      }
      fRailway->LoadAux(macro);
      TString rmacro = gSystem->Which(gROOT->GetMacroPath(), macro);
      if (rmacro.IsNull()) {
	Error("CreatePhysicsSelection", "Custom PS script %s not found",
	      macro.Data());
	return;
      }
      Info("CreatePhysicsSelection", "Loading custom PS from %s",rmacro.Data());
      TString base(gSystem->BaseName(rmacro.Data()));
      gROOT->Macro(Form("%s((AliPhysicsSelection*)%p)", base.Data(), ps));
    }

    if (opt.Contains("ALL",TString::kIgnoreCase)) {
      // --- Ignore trigger class when selecting events.  This means -
      // --- that we get offline+(A,C,E) events too ------------------
      Info("CreatePhysicsSelection", "Skipping trigger selection");
      ps->SetSkipTriggerClassSelection(true);
    }      
  }
  /** 
   * Create centrality selection, and add to manager
   * 
   * @param mc Whether this is for MC 
   */
  virtual void CreateCentralitySelection(Bool_t mc)
  {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    AliVEventHandler*   inp = mgr->GetInputEventHandler();
    if (!inp) return;

    // Possibly load as PAR
    LoadOADB();
    
    gROOT->SetMacroPath(Form("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros"
			     ":%s", gROOT->GetMacroPath()));
    AliAnalysisTaskSE* task = CoupleSECar("AddTaskMultSelection.C","false");
    FromOption(task, "AlternateOADBforEstimators", "cent-oadb", "");
    if (!task->HasBranches()) {
      // Everything except tracks since that slows does things and is
      // really only needed for reference multiplicities
      task->SetBranches("ESD:AliESDRun.,AliESDHeader.,AliESDZDC.,"
			"AliESDVZERO.,AliESDTZERO.,TPCVertex.,"
			"SPDVertex.,PrimaryVertex.,AliMultiplicity."
			"SPDPileupVertices,TrkPileupVertices,"
			"AliESDAD. "
			"AOD:header,vertices,AliAODTZERO,AliAODVZERO,"
			"AliAODZDC,AliAODAD");
    }
    return;

    // Ignore the rest - just kept for historical reasons 
    Bool_t isAOD = inp->IsA()->InheritsFrom(AliAODInputHandler::Class());
    if (isAOD) return;
    
    task = CoupleSECar("AddTaskCentrality.C",
		       Form("true,%s", isAOD ? "true" : "false"));
    AliCentralitySelectionTask* ctask = 
      dynamic_cast<AliCentralitySelectionTask*>(task);
    if (!ctask) return;
    if (mc) ctask->SetMCInput();
  }
  /** 
   * Create analysis tasks.  Must be overloaded by sub-class
   * 
   * @param mgr  Manager
   */
  virtual void CreateTasks(AliAnalysisManager* mgr)=0;
  /* @} */
  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Functions for adding Cars (tasks) 
   */
  /** 
   * Add a task using a script and possibly some arguments 
   * 
   * @param macro Script to execute 
   * @param args  Optional arguments to the script 
   * 
   * @return Created task or null
   */
  virtual AliAnalysisTask* CoupleCar(const TString& macro, 
				     const TString& args)
  {
    TString p = gSystem->Which(gROOT->GetMacroPath(), macro.Data());
    if (p.IsNull()) { 
      Error("CoupleCar", "Macro %s not found", macro.Data());
      return 0;
    }
    TString cmd(p);
    if (!args.IsNull()) 
      cmd.Append(TString::Format("(%s)", args.Data()));
    Info("CoupleCar", "Execute %s", cmd.Data());
    Int_t err;
    Long_t ret = gROOT->Macro(cmd.Data(), &err, false);
    if (!ret) { 
      Error("CoupleCar", "Failed to execute %s (%ld)", cmd.Data(), ret);
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
  virtual AliAnalysisTask* CoupleCar(const TString& macro)
  {
    TString args;
    return CoupleCar(macro, args);
  }
  /** 
   * Add a single event analysis task to the train, passing the
   * specified arguments to the macro.
   * 
   * @param macro The <b>AddTask</b> macro 
   * @param args  Arguments to pass the macro 
   * @param mask  Possible trigger mask  (if 0, no mask is set)
   * 
   * @return The added task, if any 
   */
  virtual AliAnalysisTaskSE* CoupleSECar(const TString& macro, 
					 const TString& args,
					 UInt_t         mask=0)
  {
    AliAnalysisTaskSE* task =
      dynamic_cast<AliAnalysisTaskSE*>(CoupleCar(macro, args));
    // Explicitly set mask 
    if (mask > 0) task->SelectCollisionCandidates(mask);
    return task;
  }
  /** 
   * Add a single event task to the train with no arguments passed to
   * the script
   * 
   * @param macro The <b>AddTask</b> macro. 
   * @param mask  Possible trigger mask  (if 0, no mask is set)
   *
   * @return The added task, if any
   */
  virtual AliAnalysisTaskSE* CoupleSECar(const TString& macro,
					 UInt_t         mask=0)
  {
    TString args;
    return CoupleSECar(macro, args, mask);
  }
  /** 
   * Find an already added task 
   * 
   * @param name    Name of the task 
   * @param verbose If true, 
   * 
   * @return 
   */
  virtual AliAnalysisTask* FindCar(const TString& name, 
				    Bool_t verbose=true) const
  {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
      ::Warning("FindCar", "No manager defined");
      return 0;
    }
    AliAnalysisTask* task = mgr->GetTask(name);
    if (!task && verbose)
      ::Warning("FindCar", "Task \"%s\" not found in train", 
		name.Data());
    return task;
  }
  /** 
   * Load OADB
   * 
   * @param asPar If true, possibly as par
   */
  void LoadOADB(Bool_t asPar=false)
  {
    Bool_t usePar = asPar;
    fRailway->UsePar(usePar); // Disabled pars temporarily
    fRailway->LoadLibrary("OADB");
    fRailway->UsePar(usePar); // Set to old value 
  }
  /* @} */
  /** 
   * @{ 
   * @name Tender 
   */
  virtual void AddOCDBConnect()
  {
    fRailway->LoadLibrary("PWGPP");
    Long_t ret = gROOT->ProcessLine("new AliTaskCDBconnect(\"cdb\")");
    if (!ret) {
      Fatal("AddOCDBConnect", "Failed to add OCDB connection task");
      return;
    }
    AliAnalysisTask* task = reinterpret_cast<AliAnalysisTask*>(ret);
    if (!task->HasBranches()) task->SetBranches("ESD:AliESDRun. AOD:header");
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  }  
  /** 
   * Enumeration of tenders 
   */
  enum {
    kTenderV0    = 0x0001,
    kTenderTPC   = 0x0002,
    kTenderPtFix = 0x0004,
    kTenderT0    = 0x0008,
    kTenderTOF   = 0x0010,
    kTenderTRD   = 0x0020,
    kTenderVTX   = 0x0040,
    kTenderEMCAL = 0x0080,
    kTenderPID   = 0x0100,
    kTenderHMPID = 0x0200,
    kTenderPHOS  = 0x0400
  };
  /** 
   * Add individual tender supplies 
   * 
   * @param tender Tender task 
   * @param flag   Which supply to add
   * @param debug  Debug flag for supply
   * 
   * @return Pointer to tender supply 
   */
  virtual void* AddTenderSupply(void*     tender,
				UShort_t  flag,
				Int_t     debug)
  {
    if (flag == 0) return 0;

    TString c;
    TString n;
    switch (flag) {
    case kTenderV0:    n = "VZERO";    c = "VZERO";    break;
    case kTenderTPC:   n = "TPC";      c = "TPC";      break;
    case kTenderPtFix: n = "TrackFix"; c = "PtInvFix"; break;
    case kTenderT0:    n = "T0";       c = "TZERO";    break;
    case kTenderTOF:   n = "TOF";      c = "TOF";      break;
    case kTenderTRD:   n = "TRD";      c = "TRD";      break;
    case kTenderVTX:   n = "Vtx";      c = "IP";       break;
    case kTenderEMCAL: n = "EMCAL";    c = "EMCAL";    break;
    case kTenderPID:   n = "PID";      c = "PID";      break;
    case kTenderHMPID: n = "HMPID";    c = "HMPID";    break;
    case kTenderPHOS:  n = "PHOS";     c = "PHOS";     break;
    default:     
      Warning("AddTenderSupply", "Unknown tender flag: 0x%08x", flag);
      return 0;
    }
    TString m;
    switch (flag) {
    case kTenderV0:    m = Form("s->SetDebug(%d);", debug);     break;
    case kTenderEMCAL: // Fall through 
    case kTenderTOF:   // Fall through 
    case kTenderTPC:   // Fall through 
    case kTenderPtFix: // Fall through 
    case kTenderTRD:   m = Form("s->SetDebugLevel(%d);", debug); break;
    }

    TString cls(Form("Ali%sTenderSupply", c.Data()));
    Long_t ret =
      gROOT->ProcessLine(Form("{ %s* s = new %s(\"%s\");"
			      "((AliTender*)%p)->AddSupply(s);%s"
			      "} s;",
			      cls.Data(),
			      cls.Data(),
			      n.Data(),
			      tender,
			      m.Data()));
    void* ptr = reinterpret_cast<void*>(ret);
    Info("AddTenderSupply", "Adding supply %s (an %s object): %p",
	 n.Data(), cls.Data(), ret);
    return ptr;
  }
  /** 
   * Add tender to train.  
   * 
   */
  virtual void AddTender(const TString& sup)
  {
    TString supplies = sup;

    UShort_t which = 0;
    supplies.ToUpper();
    if (supplies.Contains("V0") ||
	supplies.Contains("VZERO")) which |= kTenderV0    ;
    if (supplies.Contains("TPC"))   which |= kTenderTPC   ;
    if (supplies.Contains("PTFIX")) which |= kTenderPtFix ;
    if (supplies.Contains("T0"))    which |= kTenderT0    ;
    if (supplies.Contains("TOF"))   which |= kTenderTOF   ;
    if (supplies.Contains("TRD"))   which |= kTenderTRD   ;
    if (supplies.Contains("VTX"))   which |= kTenderVTX   ;
    if (supplies.Contains("EMCAL")) which |= kTenderEMCAL ;
    if (supplies.Contains("PID"))   which |= kTenderPID   ;

    AddTender(which);
  }
  /** 
   * Add a tender with supplies specified in argument 
   *
   *
   */
  virtual void AddTender(UShort_t which)
  {
    if (which == 0) return;
    
    Int_t  vrb      = fOptions.AsInt("verbose", 3);

    fRailway->LoadLibrary("Tender");
    fRailway->LoadLibrary("TenderSupplies");
    
    Long_t ret = gROOT->ProcessLine("new AliTender(\"Tender\")");
    if (!ret) {
      Warning("AddTender", "Failed to make tender");
      return;
    }
    void* tender = reinterpret_cast<void*>(ret);
    gROOT->ProcessLine(Form("((AliTender*)%p)->SetCheckEventSelection(%d)",
			    tender, (which & kTenderV0)));
    gROOT->ProcessLine(Form("((AliTender*)%p)->SetDebugLevel(%d)",
			    tender, vrb));
    
    // OCDB settings for tender 
    TString ocdb = fOptions.Get("ocdb");
    if (ocdb.IsNull())
      ocdb = "raw://";
    else if (ocdb.EndsWith(".root")) {
      fRailway->LoadLibrary("CDB");
      fRailway->LoadAux(ocdb);
      gROOT->ProcessLine(Form("AliCDBManager::Instance()->"
			      "SetSnapshotMode(\"%s\");",
			      ocdb.Data()));
      ocdb = "raw://";
    }
    gROOT->ProcessLine(Form("((AliTender*)%p)->SetDefaultCDBStorage(\"%s\")",
			    tender, ocdb.Data()));

    AddTenderSupply(tender, which & kTenderV0,     vrb);
    AddTenderSupply(tender, which & kTenderTPC,    vrb);
    AddTenderSupply(tender, which & kTenderPtFix,  vrb);
    AddTenderSupply(tender, which & kTenderT0,     vrb);
    AddTenderSupply(tender, which & kTenderTOF,    vrb);
    AddTenderSupply(tender, which & kTenderTRD,    vrb);
    AddTenderSupply(tender, which & kTenderVTX,    vrb);
    AddTenderSupply(tender, which & kTenderEMCAL,  vrb);
    AddTenderSupply(tender, which & kTenderPID,    vrb);

    gROOT->ProcessLine(Form("((AliTender*)%p)->GetSupplies()->Print()",tender));
    
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    AliAnalysisTask*    tsk = reinterpret_cast<AliAnalysisTask*>(tender);
    mgr->AddTask(tsk);
    AliAnalysisDataContainer* cnt =
      mgr->CreateContainer("tender_event",
			   AliESDEvent::Class(),
			   AliAnalysisManager::kExchangeContainer,
			   "default_tender");
    mgr->ConnectInput (tsk, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(tsk, 1, cnt);
    
  }
  /* @} */
  //------------------------------------------------------------------
  /** 
   * @{ 
   * @name Set parameters on a task 
   */
  /** 
   * Set a unsigned integer parameter on the task 
   * 
   * @param task Task 
   * @param what What to set 
   * @param val  Value 
   */
  void SetOnTask(AliAnalysisTaskSE* task,
		 const char* what,
		 UInt_t val)
  {
    const char* cls = task->ClassName();
    gROOT->ProcessLine(Form("((%s*)%p)->Set%s(%u)",cls,task,what,val));
  }    
  /** 
   * Set a integer parameter on the task 
   * 
   * @param task Task 
   * @param what What to set 
   * @param val  Value 
   */
  void SetOnTask(AliAnalysisTaskSE* task,
		 const char* what,
		 Int_t val)
  {
    const char* cls = task->ClassName();
    gROOT->ProcessLine(Form("((%s*)%p)->Set%s(%d)",cls,task,what,val));
  }    
  /** 
   * Set a integer parameter on the task 
   * 
   * @param task Task 
   * @param what What to set 
   * @param val  Value 
   */
  void SetOnTask(AliAnalysisTaskSE* task,
		 const char* what,
		 Long64_t val)
  {
    const char* cls = task->ClassName();
    gROOT->ProcessLine(Form("((%s*)%p)->Set%s(%lld)",cls,task,what,val));
  }    
  /** 
   * Set a real parameter on the task 
   * 
   * @param task Task 
   * @param what What to set 
   * @param val  Value 
   */
  void SetOnTask(AliAnalysisTaskSE* task,
		 const char* what,
		 Double_t val)
  {
    const char* cls = task->ClassName();
    gROOT->ProcessLine(Form("((%s*)%p)->Set%s(%lg)",cls,task,what,val));
  }    
  /** 
   * Set a boolean parameter on the task 
   * 
   * @param task Task 
   * @param what What to set 
   * @param val  Value 
   */
  void SetOnTask(AliAnalysisTaskSE* task,
		 const char* what,
		 Bool_t val)
  {
    const char* cls = task->ClassName();
    gROOT->ProcessLine(Form("((%s*)%p)->Set%s(%d)",cls,task,what,val));
  }    
  /** 
   * Set a string parameter on the task 
   * 
   * @param task Task 
   * @param what What to set 
   * @param val  Value 
   */
  void SetOnTask(AliAnalysisTaskSE* task,
		 const char* what,
		 const char* val)
  {
    const char* cls = task->ClassName();
    gROOT->ProcessLine(Form("((%s*)%p)->Set%s(\"%s\")",cls,task,what,val));
  }
  /* @} */
  //------------------------------------------------------------------
  /** 
   * @{
   * @name Set parameters on tasks based on options 
   */
  /** 
   * Set a real parameter on the task based on option value 
   * 
   * @param task    Task 
   * @param what    What to set 
   * @param opt     Option name 
   * @param defval  Default value if option not given
   */
  void FromOption(AliAnalysisTaskSE* task,
		  const char* what,
		  const char* opt,
		  Double_t defval)
  {
    Double_t val = fOptions.AsDouble(opt,defval);;
    SetOnTask(task, what, val);
  }
  /** 
   * Set a boolean parameter on the task based on option value 
   * 
   * @param task    Task 
   * @param what    What to set 
   * @param opt     Option name 
   */
  void FromOption(AliAnalysisTaskSE* task,
		  const char* what,
		  const char* opt,
		  Bool_t      /*defval*/)
  {
    Bool_t val = fOptions.AsBool(opt);;
    SetOnTask(task, what, val);
  }
  /** 
   * Set a string parameter on the task based on option value 
   * 
   * @param task    Task 
   * @param what    What to set 
   * @param opt     Option name 
   * @param defval  Default value if option not given
   */
  void FromOption(AliAnalysisTaskSE* task,
		 const char* what,
		 const char* opt,
		 const char* defval)
  {
    TString val = fOptions.AsString(opt,defval);
    SetOnTask(task, what, val.Data());
  }
  /** 
   * Set a integer parameter on the task based on option value 
   * 
   * @param task    Task 
   * @param what    What to set 
   * @param opt     Option name 
   * @param defval  Default value if option not given
   */
  void FromOption(AliAnalysisTaskSE* task,
		  const char* what,
		  const char* opt,
		  Int_t      defval)
  {
    Int_t   val = fOptions.AsInt(opt, defval);
    SetOnTask(task, what, val);
  }
  /** 
   * Set a integer parameter on the task based on option value 
   * 
   * @param task    Task 
   * @param what    What to set 
   * @param opt     Option name 
   * @param defval  Default value if option not given
   */
  void FromOption(AliAnalysisTaskSE* task,
		  const char* what,
		  const char* opt,
		  Long64_t    defval)
  {
    Long64_t   val = fOptions.AsLong(opt, defval);
    SetOnTask(task, what, val);
  }
  /* @} */

  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name Train query utilities 
   */
  /** 
   * Check if we have an MC handler attached 
   * 
   * @return True if MC handler is found in a valid manager.  False if
   * manager is not defined, or has no MC handler.
   */
  virtual Bool_t HasMCHandler() const 
  {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return false;
    return mgr->GetMCtruthEventHandler() != 0;
  }
  /* @} */
  /** 
   * Set the name of the train - should be name of the class.  Must be
   * overloaded.
   * 
   * @return Class name as a constant C string 
   */
  virtual const Char_t* ClassName() const = 0;
  /* @} */
  //__________________________________________________________________
  virtual void AddMonitor(const TString& name)
  {
    if (fRailway->Mode() != Railway::kProof) return;
    Warning("CreateMonitors", "Monitoring not supported yet");
  }
  //__________________________________________________________________
  /** 
   * Create the monitors
   * 
   */
  virtual void CreateMonitors() 
  {
  }
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
  static TString EscapeName(const char* name, TString& datimeStr)
  {
    TString escaped = name;
    char  c[] = { ' ', '/', '@', 0 };
    char* p   = c;
    while (*p) { 
      char tmp[] = { *p, '\0' };
      escaped.ReplaceAll(tmp, "_");
      p++;
    }
    if (!datimeStr.IsNull() && 
	!datimeStr.EqualTo("none", TString::kIgnoreCase)) {
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
      datimeStr = Form("%04d%02d%02d_%02d%02d", 
		       datime.GetYear(), 
		       datime.GetMonth(), 
		       datime.GetDay(), 
		       datime.GetHour(), 
		       datime.GetMinute());
      escaped.Append(Form("_%s", datimeStr.Data()));
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
    if (fRailway->Operation() == Railway::kTerminate && !exists) {
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
  /* @} */
  // -----------------------------------------------------------------
  /** 
   * @{ 
   * @name For post-processing 
   */
  /** 
   * Save the setup as a ROOT script and possibly also a shell script
   * 
   * @param asShellScript If true, also save as shell script
   */
  virtual void SaveSetup(Bool_t asShellScript)
  {
    OptionList tmp(fOptions);
    const OptionList* uopts = (fRailway ? &fRailway->Options() : 0);
    if (tmp.Find("overwrite")) tmp.Set("overwrite");
    if (tmp.Find("date") && fEscapedName.Length() > fName.Length()+1) {
      Int_t n = fName.Length()+1;
      tmp.Set("date", fEscapedName(n, fEscapedName.Length()-n));
    }
    if (asShellScript) 
      SaveSetupShell("rerun", ClassName(), fName, tmp, uopts);
    SaveSetupROOT("ReRun", ClassName(), fName, tmp, uopts);
    if (fRailway) fRailway->AuxSave(fEscapedName, asShellScript);
    SavePostShellScript();
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
    TString url = opts.Get("url");
    OptionList tmp(opts);
    tmp.Set("url", "${url}");
    o << "#!/bin/bash\n\n"
      << "class=\"" << cls << "\"\n"
      << "name=\"" << name << "\"\n"
      << "url=\"" << url << "\"\n\n" 
      << "# Available options\n"
      << "# \n";
    tmp.Help(o, "#    --");
    if (uopts) {
      o << "#\n"
	<< "# Available URI options\n"
	<< "# \n";
      uopts->Help(o, "#      ");
    }
    o << "#\n"
      << "opts=(--class=$class \\\n"
      << "  --name=$name";
    tmp.Store(o, " \\\n  --", "", true);
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
      << "  path.Append(\"$ALICE_PHYSICS/PWGLF/FORWARD/trains\");\n"
      << "  gROOT->SetMacroPath(path);\n\n"
      << "  gROOT->LoadMacro(\"RunTrain.C\");\n\n"
      << "  return RunTrain(name, cls, uri, opts);\n"
      << "}\n" 
      << "//\n"
      << "// EOF\n" 
      << "//" << std::endl;
    o.close();
  }    
  /** 
   * Write shell code to do post processing after terminate.  This
   * code should deal with a single run (or run range).  The following
   * shell variables are available to the code:
   *
   * - @c $prefix  Relative path to job directory or empty 
   * - @c $dest    Destination for output to be stored
   * 
   * Note, the code is injected into a shell function, and should
   * therefor not define new functions or the like.
   * 
   * @param o The output stream.  
   */
  virtual void PostShellCode(std::ostream& o)
  {
    o << "  echo \"Nothing to do for " << ClassName() 
      << " train\"" << std::endl;
  }
  /** 
   * Save a script to do post processing after terminate on each run
   * or run-range. 
   * 
   * The shell script will execute the train defined code (in
   * PostShellCode) for each run or run-range.  The train defined code
   * and call drawing and summarizing macros or the like.
   *
   * In case of Grid analysis, this script will download and extract
   * the appropriate ZIP files to separate directories, and then
   * change directory to these directories and execute the train
   * defined shell code there.  In this case, the script defines the
   * shell variable @c $prefix as the relative path to the job
   * directory.
   * 
   */
  void SavePostShellScript()
  {
    std::ofstream f("post.sh");
    if (!f) { 
      Error("SavePostAll", "Failed to open post.sh script");
      return;
    }
    f << "#!/bin/bash\n"
      << "# Generated by " << ClassName() << "\n"
      << "set -e\n"
      << "\n"
      << "dest=$1\n"
      << "prefix=\n"
      << "\n"
      << "doall() {"
      << std::endl;
    PostShellCode(f);
    f << "}\n"
      << "\n"
      << "if test ! -f Download.C ;then\n"
      << "  doall\n"
      << "  exit\n"
      << "fi\n"
      << "\n"
      << "if test ! -f .download ; then\n"
      << "  aliroot -l -b -q Download.C\\(1\\)\n"
      << "  touch .download\n"
      << "fi\n"
      << "prefix=../\n"
      << "\n"
      << "for i in root_archive_*.zip ; do\n"
      << "  d=`basename $i .zip` \n"
      << "  if test ! -d $d ; then\n"
      << "    echo \"Directory $d missing\"\n"
      << "    continue\n"
      << "  fi\n"
      << "  \n"
      << "  (cd $d && doall)\n"
      << "done\n"
      << "# EOF"
      << std::endl;
    f.close();
    gSystem->Exec("chmod a+x post.sh");
  }    
  /* @} */
  TString      fName;
  TString      fEscapedName;
  TString      fDatimeString;
  OptionList   fOptions;
  Railway*     fRailway;
};
#endif
