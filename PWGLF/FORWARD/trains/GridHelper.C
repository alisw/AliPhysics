/**
 * @file   GridHelper.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 19:01:27 2012
 * 
 * @brief  Grid Analysis Helper
 * 
 * @ingroup pwglf_forward_trains_helper
 * 
 */
#ifndef GRIDHELPER_C
#define GRIDHELPER_C
#include "PluginHelper.C"
#ifndef __CINT__
# include <TUrl.h>
# include <TString.h>
# include <TGrid.h>
# include <AliAnalysisManager.h>
# include <AliAnalysisAlien.h>
#else
class TUrl;
class AliAnalysisAlien;
#endif

// ===================================================================
/**
 * Handle analysis on an the Grid
 * 
 * This helper is triggered by a URL of the form 
 *
 * @code
 * alien:///<directory>[?<options>][#<pattern>]
 * @endcode 
 * where 
 * <dl>
 *   <dt>&lt;directory&gt;</dt>
 *   <dd>Grid directory that holds the data</dd>
 *   <dt>&lt;treeName&gt;</dt>
 *   <dd>Tree to loop over</dd>
 *   <dt>&lt;options&gt;</dt>
 *   <dd>List of options separated by an &amp;
 *     <dl>
 *       <dt><tt>storage=&lt;url&gt;</tt></dt>
 *       <dd>Specify a non-default storage location for special output
 *         (e.g., AOD trees).  &lt;url&gt; should be a valid XRootd 
 *         server URI accessible to the slaves - e.g., 
 *         <tt>root://lxplus.cern.ch:10930//tmp</tt>.</dd>
 *       <dt><tt>mode=[default,rec,sim,train,custom]</tt></dt>
 *       <dd>Set the AliROOT mode.  If not specified <tt>default</tt> 
 *         is assumed.  See also CreateAliROOTPar</dd>
 *       <dt><tt>par</tt></dt>
 *       <dd> Use PAR files</dd>
 *       <dt><tt>runs=[list or file]</tt></dt>
 *       <dd>Comma separated list of run numbers, or file(s) containing 
 *         run numbers</dd> 
 *       <dt><tt>oper=[FULL,TERMINATE,SUBMIT,OFFLINE,TEST]</tt></dt>
 *       <dd>How to run the analysis</dd>
 *       <dt><tt>split=&lt;N&gt;</tt></dt>
 *       <dd>Maximum number of files per split</dd>
 *       <dt><tt>merge=&lt;N&gt;</tt></dt>
 *       <dd>Maximum number of files per merger</dd>
 *       <dt><tt>mc</tt></dt>
 *       <dd>Scan also for MC files (<tt>galice.root</tt>, 
 *          <tt>Kinematics.root</tt>, and <tt>TrackRefs.root</tt>) when 
 *          scanning &lt;datadir&gt;</dd>
 *       <dt><tt>pattern=&lt;GLOB&gt;</tt></dt>
 *       <dd>Shell glob pattern that files must check when scanning 
 *         &lt;datadir&gt;</dd>
 *     </dl>
 *   </dd>
 * </dl>  
 *
 * @ingroup pwglf_forward_trains_helper
 */
struct GridHelper : public PluginHelper
{
  /** 
   * Constructor 
   * 
   * @param url  Url 
   * @param verbose Verbosity level
   */
  GridHelper(const TUrl& url, Int_t verbose)
    : PluginHelper(url, verbose), fRuns()
  {
    fOptions.Add("oper", "FULL|TERMINATE|SUBMIT", "Analysis operation", "FULL");
    fOptions.Add("split",  "N|max",  "Max number of files before split","max");
    fOptions.Add("merge",  "N|max",  "Max number of files for merge",   "max");
    fOptions.Add("run",    "RUNS",   "Range, list, and/or file of runs");
    fOptions.Add("pattern","GLOB",   "File/directory name pattern");
    fOptions.Add("alien",  "VERSION","Alien API version",              "V1.1x");
    fOptions.Add("concat", "Concatenate all runs");
    fOptions.Add("ttl",    "N|max",  "Time to live",                    "max");
  }
  GridHelper(const GridHelper& o)
    : PluginHelper(o), fRuns()
  {}
  GridHelper& operator=(const GridHelper& o)
  {
    if (&o == this) return *this;
    PluginHelper::operator=(o);
    return *this;
  }
  virtual ~GridHelper() {}
  /** 
   * Get the mode identifier 
   * 
   * @return Always kProof
   */
  virtual UShort_t Mode() const { return kGrid; }
  /**
   * Get the mode string used for AliAnalysisManager::StartAnalysis
   */
  virtual const char* ModeString() const { return "grid"; }
  /** 
   * Set-up done before task set-ups 
   * 
   * @return true on success 
   */
  virtual UShort_t Operation() const 
  {
    if (!fOptions.Has("oper")) return kFull;
    const TString& oper = fOptions.Get("oper");
    if      (oper.EqualTo("FULL",      TString::kIgnoreCase)) return kFull;
    else if (oper.EqualTo("OFFLINE",   TString::kIgnoreCase)) return kOffline;
    else if (oper.EqualTo("SUBMIT",    TString::kIgnoreCase)) return kSubmit;
    else if (oper.EqualTo("TERMINATE", TString::kIgnoreCase)) return kTerminate;
    else if (oper.EqualTo("TEST",      TString::kIgnoreCase)) return kTest;
    return kFull;
  }
  void StoreRun(Int_t r)
  {
    TObject* o = new TObject;
    o->SetUniqueID(r);
    fRuns.Add(o);
  }
  /**
   * Read run numbers 
   *
   * @return Number of registered runs 
   */
  virtual Int_t RegisterRuns()
  {
    if (!fOptions.Find("run")) {
      Error("GridHelper::RegisterRuns", "No runs specified");
      return -1;
    }
    Int_t       nRuns  = 0;
    TString     runs   = fOptions.Get("run");
    TObjArray*  tokens = runs.Tokenize(",+.");
    TObjString* part   = 0;
    TIter       next(tokens);
    Bool_t      range  = false;
    Bool_t      individual = false;
    // Info("GridHelper::RegisterRuns", "Runs specified are %s", runs.Data());
    while ((part = static_cast<TObjString*>(next()))) {
      TString& s = part->String();
      if (s.Contains("-")) { // Run range 
	if (range) { 
	  Warning("GridHelper::RegisterRuns", "Run range already specified, "
		  "ignoring %s", s.Data());
	  continue;
	}
	if (individual) { 
	  Warning("GridHelper::RegisterRuns", 
		  "Run ranges and individual run specs do not mix, "
		  "ignoring %s", s.Data());
	  continue;
	}
	TObjArray* ranges = s.Tokenize("-");
	if (ranges->GetEntriesFast() > 2) { 
	  Warning("GridHelper::RegisterRuns", "Invalid run range: %s", 
		  s.Data());
	  ranges->Delete();
	  continue;
	}
	Int_t first = static_cast<TObjString*>(ranges->At(0))->String().Atoi();
	Int_t last  = static_cast<TObjString*>(ranges->At(1))->String().Atoi();
	nRuns       = last-first+1;
	// Info("GridHelper::RegisterRuns", "Run range %d -> %d", first, last);
	fHandler->SetRunRange(first, last);
	ranges->Delete();
	range = true;
	for (Int_t r = first; r <= last; r++) StoreRun(r);
	continue;
      }
      if (s.IsDigit()) { // single run
	if (range) { 
	  Warning("GridHelper::RegisterRuns", 
		  "Run ranges and individual run specs do not mix, "
		  "ignoring %s", s.Data());
	  continue;
	}
	// Info("GridHandler::RegisterRuns", "Adding run %s", s.Data());
	fHandler->AddRunNumber(s.Atoi());
	StoreRun(s.Atoi());
	nRuns++;
	individual = true;
	continue;
      }
      if (range) { 
	Warning("GridHelper::RegisterRuns", "Run ranges and list file "
		"do not mix, ignoring %s", s.Data());
	continue;
      }

      // We assume this part is a file 
      // Info("GridHelper::RegisterRuns", "Reading runs from %s", s.Data());
      std::ifstream in(s.Data());
      if (!in) { 
	s.Prepend("../");
	in.open(s.Data());
	if (!in) {
	  Warning("GridHelper::RegisterRuns", "Failed to open %s", s.Data());
	  continue;
	}
      }
      while (!in.eof()) { 
	Int_t r;
	in >> r;
	// Info("GridHelper::RegisterRuns", "Read %d, adding", r);
	fHandler->AddRunNumber(r);
	StoreRun(r);
	nRuns++;
	Char_t c;
	in >> c;
	if (in.bad()) break;
      }
      individual = true;
      in.close();
    }
    return nRuns;
  }
  /** 
   * Executed before setting up tasks 
   * 
   * @return true on success 
   */
  virtual Bool_t PreSetup() 
  {
    if (!PluginHelper::PreSetup()) return false;
    
    // --- Add system library dir to load path -----------------------
    gSystem->AddDynamicPath("/usr/lib");

    // --- Open a connection to the grid -----------------------------
    if (!TGrid::Connect(Form("%s://", fUrl.GetProtocol()))) { 
      Error("GridHelper::PreSetup", "Failed to connect to AliEN");
      return false;
    }
    if (!gGrid || !gGrid->IsConnected()) { 
      Error("GridHelper::PreSetup", "Failed to connect to AliEN");
      return false;
    }

    return true;
  }
  /** 
   * Set-up done after the task set-ups 
   *
   * @return true on success 
   */
  virtual Bool_t PostSetup() 
  {
    // Info("GridHelper::PostSetup", "Calling super.PostSetup");
    if (!PluginHelper::PostSetup()) return false;

    // --- API version -----------------------------------------------
    fHandler->SetAPIVersion(fOptions.Get("alien"));
    
    // --- Get the name ----------------------------------------------
    // Info("GridHelper", "Proceeding with plugin setup");
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    TString name(mgr->GetName());

    // --- Set the operation to do (TEST, SUBMIT, TERMINATE, FULL) ---
    TString operation("FULL");
    if (fOptions.Has("oper")) operation = fOptions.Get("oper");
    fHandler->SetRunMode(operation);

    // --- Add the run numbers ---------------------------------------
    fHandler->SetRunPrefix(fOptions.Has("mc") ? "%d" : "%09d");
    Int_t nRun = RegisterRuns();

    // --- Do not test copying ---------------------------------------
    fHandler->SetCheckCopy(false);
    
    // --- Set output to be per run ----------------------------------
    fHandler->SetOutputToRunNo(true); 

    // --- Set the job tag -------------------------------------------
    fHandler->SetJobTag(name);

    // --- Set number of test files - used in test mode only ---------
    fHandler->SetNtestFiles(1);

    // --- Set the Time-To-Live --------------------------------------
    if (fOptions.Has("ttl")) { 
      if (!fOptions.Get("ttl").EqualTo("max")) {
	fHandler->SetTTL(fOptions.AsInt("ttl"));
      }
    }
    
    // --- Re-submit failed jobs as long as the ratio of failed jobs -
    // --- is this percentage.
    fHandler->SetMasterResubmitThreshold(95);

    // --- Set the input format --------------------------------------
    fHandler->SetInputFormat("xml-single");

    // --- Set names of generated files ------------------------------
    fHandler->SetAnalysisMacro(Form("%s.C", name.Data()));
    fHandler->SetJDLName(Form("%s.jdl", name.Data()));
    fHandler->SetExecutable(Form("%s.sh", name.Data()));
    
    // ---- Set the job price !? -------------------------------------
    fHandler->SetPrice(1);

    // --- Set whether to merge via JDL ------------------------------
    fHandler->SetMergeViaJDL(true);
    
    // --- Fast read otion -------------------------------------------
    fHandler->SetFastReadOption(false);

    // --- Whether to overwrite existing output ----------------------
    fHandler->SetOverwriteMode(true);

    // --- Set the executable binary name and options ----------------
    fHandler->SetExecutableCommand("aliroot -b -q -x");

    // --- Split by storage element - must be lower case! ------------
    fHandler->SetSplitMode("se");

    // --- How much to split -----------------------------------------
    if (fOptions.Has("split")) { 
      if (!fOptions.Get("split").EqualTo("max")) {
	fHandler->SetSplitMaxInputFileNumber(fOptions.AsInt("split"));
      }
    }
    // --- Merge parameters ------------------------------------------
    if (fOptions.Has("merge")) { 
      if (!fOptions.Get("merge").EqualTo("max")) { 
	fHandler->SetMaxMergeFiles(fOptions.AsInt("merge"));
      }
    }
    fHandler->SetMergeExcludes("AliAOD.root *EventStat*.root "
			       "*event_stat*.root");
    
    // --- Set number of runs per master - 1 or all ------------------
    fHandler->SetNrunsPerMaster(fOptions.Has("concat") ? nRun+1 : 1);


    // --- Enable default outputs ------------------------------------
    fHandler->SetDefaultOutputs(true);

    // --- Keep log files ------------------------------------------
    fHandler->SetKeepLogs();

    // --- Set the working directory to be the trains name (with -----
    // --- special characters replaced by '_' and the date appended),
    // --- and also set the output directory (relative to working
    // --- directory)
    fHandler->SetGridWorkingDir(name.Data());
    fHandler->SetGridOutputDir("output");
    fHandler->SetGridDataDir(fUrl.GetFile());

    // --- Get the tree name and set the file pattern ----------------
    TString pattern;
    if (fOptions.Has("pattern")) pattern = fOptions.Get("pattern");
    else {
      TString treeName(fUrl.GetAnchor());
      if (treeName.IsNull()) { 
	Warning("GridHelper::PreSetup", "No tree name specified, assuming T");
	treeName = "T";
      }
      if      (treeName.EqualTo("esdTree")) pattern = "AliESD";
      else if (treeName.EqualTo("aodTree")) pattern = "AliAOD";
    }
    fHandler->SetDataPattern(pattern);

    // --- Loop over defined containers in the analysis manager, and -
    // --- declare these as outputs
    TString listOfAODs  = "";
    TString listOfHists = "";
    TString listOfTerms = "";

    TObjArray*  outs[] = { mgr->GetOutputs(), mgr->GetParamOutputs(), 0 };
    TObjArray** out    = outs;
    while (*out) {
      AliAnalysisDataContainer* cont = 0;
      TIter nextCont(*out);
      while ((cont = static_cast<AliAnalysisDataContainer*>(nextCont()))) {
	TString outName(cont->GetFileName());
	Bool_t   term = (*out == outs[1]);
	TString& list = (outName == "default" ? listOfAODs : 
			 !term ? listOfHists : listOfTerms);
	if (outName == "default") { 
	  if (!mgr->GetOutputEventHandler()) continue; 
	  
	  outName = mgr->GetOutputEventHandler()->GetOutputFileName();
	}
	if (list.Contains(outName)) continue;
	if (!list.IsNull()) list.Append(",");
	list.Append(outName);
      }
      out++;
    }
    TString extra = mgr->GetExtraFiles();
    if (!extra.IsNull()) { 
      if (!listOfAODs.IsNull()) listOfAODs.Append("+");
      extra.ReplaceAll(" ", ",");
      listOfAODs.Append(extra);
   }

#if 0
    Int_t nReplica = 2;
    TString outArchive = Form("stderr, stdout@disk=%d", nReplica);
    if (!listOfHists.IsNull()) 
      outArchive.Append(Form(" hist_archive.zip:%s@disk=%d", 
			     listOfHists.Data(), nReplica));
    if (!listOfAODs.IsNull()) 
      outArchive.Append(Form(" aod_archive.zip:%s@disk=%d", 
			     listOfAODs.Data(), nReplica));
    // Disabled for now 
    // plugin->SetOutputArchive(outArchive);
#endif 

    if (listOfAODs.IsNull() && listOfHists.IsNull()) 
      Fatal("PostSetup", "No outputs defined");
    if (!listOfTerms.IsNull()) 
      fHandler->SetTerminateFiles(listOfTerms);
    
    return true;
  };
  /** 
   * Start the analysis 
   * 
   * @param nEvents Number of events to analyse 
   * 
   * @return The return value of AliAnalysisManager::StartAnalysis
   */
  virtual Long64_t Run(Long64_t nEvents=-1) 
  {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (nEvents == 0) return 0;
    Long64_t ret = mgr->StartAnalysis("grid", nEvents);

#if 1
    std::ofstream outJobs(Form("%s.jobid", mgr->GetName()));
    outJobs << fHandler->GetGridJobIDs() << std::endl;
    outJobs.close();

    std::ofstream outStages(Form("%s.stage", mgr->GetName()));
    outStages << fHandler->GetGridStages() << std::endl;
    outStages.close();
#endif
    return ret;
  }
  /** 
   * Link an auxilary file to working directory 
   * 
   * @param name Name of the file
   * @param copy  Whether to copy or not 
   * 
   * @return true on success
   */
  virtual Bool_t AuxFile(const TString& name, bool copy=false)
  {
    if (!Helper::AuxFile(name, copy)) return false;
    // We need to add this file as an additional 'library', so that the 
    // file is uploaded to the users Grid working directory. 
    fHandler->AddAdditionalLibrary(gSystem->BaseName(name.Data()));
    return true;
  }
  /** 
   * Get the output (directory)
   *
   */
  virtual TString OutputPath() const
  {
    TString ret;
    if (!fHandler) {
      Warning("GridHelper::OutputLocation", "No AliEn handler");
      return ret;
    }
    ret = fHandler->GetGridOutputDir();
    if (ret.BeginsWith("/")) return ret;

    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) { 
      Warning("GridHelper::OutputLocation", "No analysis manager");
      return ret;
    }
    ret.Prepend(Form("%s/",  mgr->GetName()));
    if (gGrid) 
      ret.Prepend(Form("%s/", gGrid->GetHomeDirectory())); 
    
    return ret;
  }
  /** 
   * @return URL help string
   */
  virtual const Char_t* UrlHelp() const 
  {
    return "alien:///<datadir>[?<options>][#<treeName>]";
  }
  /** 
   * @return Short description
   */
  virtual const char* Desc() const { return "AliEn"; }
  /** 
   * Write auxillary ROOT (and possible shell) script for more 
   * (post-)processing e.g., terminate
   * 
   * @param escaped        Escaped name  
   */
  void AuxSave(const TString& escaped, 
	       Bool_t /*asShellScript*/) 
  {
    // Write plug-in to file 
    TFile* plug = TFile::Open(Form("%s_plugin.root", escaped.Data()), 
			      "RECREATE");
    fHandler->Write("plugin");
    plug->Close();
    
    TIter       nextLib(&fExtraLibs);
    TObjString* lib = 0;
    TString     libs;
    while ((lib = static_cast<TObjString*>(nextLib()))) {
      if (!libs.IsNull()) libs.Append(" ");
      libs.Append(lib->String());
    }
    TIter       nextPar(&fExtraPars);
    TObjString* par = 0;
    TString     pars;
    while ((par = static_cast<TObjString*>(nextPar()))) {
      if (!pars.IsNull()) pars.Append(" ");
      pars.Append(par->String());
    }
    TIter       nextSrc(&fExtraSrcs);
    TObjString* src = 0;
    TString     srcs;
    while ((src = static_cast<TObjString*>(nextSrc()))) {
      if (!srcs.IsNull()) srcs.Append(" ");
      srcs.Append(src->String());
    }
    TString macDir("$ALICE_ROOT/PWGLF/FORWARD/trains");
    std::ofstream t("Terminate.C");
    if (!t) { 
      Error("GridHelper::AuxSave", "Failed to make terminate ROOT script");
      return;
    }

    t << "// Generated by GridHelper\n"
      << "Bool_t Terminate()\n"
      << "{\n"
      << "  TString name = \"" << escaped << "\";\n"
      << "  TString libs = \"" << libs << "\";\n"
      << "  TString pars = \"" << pars << "\";\n"
      << "  TString srcs = \"" << srcs << "\";\n\n"
      << "  gROOT->LoadMacro(\"" << macDir << "/GridTerminate.C\");\n\n"
      << "  return GridTerminate(name,libs,pars,srcs);\n"
      << "}\n"
      << "// EOF\n"
      << std::endl;
    t.close();

    TString runs;
    TString format(fOptions.Has("mc") ? "%d" : "%09d");
    if (fOptions.Has("concat")) {
      Int_t first = fRuns.First()->GetUniqueID();
      Int_t last  = fRuns.Last()->GetUniqueID();
      TString fmt(format); 
      fmt.Append("_"); 
      fmt.Append(format);
      if (!runs.IsNull()) runs.Append(" ");
      runs.Append(TString::Format(fmt, first, last));
    }
    else {
      TIter next(&fRuns);
      TObject* o = 0;
      while ((o = next())) { 
	if (!runs.IsNull()) runs.Append(" ");
	runs.Append(Form(format, o->GetUniqueID()));
      }
    }

    std::ofstream d("Download.C");
    if (!d) { 
      Error("GridHelper::AuxSave", "Failed to make ROOT script Download.C");
      return;
    }
    d << "// Generated by GridHelper\n"
      << "void Download()\n"
      << "{\n"
      << "  TString base = \"" << fUrl.GetProtocol() << "://" 
      << OutputPath() << "\";\n"
      << "  TString runs = \"" << runs << "\";\n\n"
      << "  gROOT->LoadMacro(\"" << macDir << "/GridDownload.C\");\n\n"
      << "  GridDownload(base, runs);\n"
      << "}\n"
      << "// EOF\n"
      << std::endl;
    d.close();

    std::ofstream w("Watch.C");
    if (!w) {
      Error("GridHelper::AuxSave", "Failed to make ROOT script Watch.C");
      return;
    }
    w << "// Generated by GridHelper\n"
      << "void Watch(Int_t delay=5*60)\n"
      << "{\n"
      << "  TString name = \"" << escaped << "\";\n"
      << "  gROOT->LoadMacro(\"" << macDir << "/GridWatch.C+g\");\n\n"
      << "  GridWatch(name,delay);\n"
      << "}\n"
      << "// EOF\n"
      << std::endl;
    w.close();

  }
  TList fRuns;
};
#endif
//
// EOF
//
