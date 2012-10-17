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
 *   <dt>&lt;directory@gt;</dt>
 *   <dd>Grid directory that holds the data</dd>
 *   <dt>&lt;treeName@gt;</dt>
 *   <dd>Tree to loop over</dd>
 *   <dt>&lt;options@gt;</dt>
 *   <dd>List of options separated by an &amp;
 *     <dl>
 *       <dt><tt>storage=&lt;url&gt;</tt></dt>
 *       <dd>Specify a non-default storage location for special output
 *         (e.g., AOD trees).  &lt;url&gt; should be a valid XRootd 
 *         server URI accessible to the slaves - e.g., 
 *         <tt>root://lxplus.cern.ch:10930//tmp</tt>.</dd>
 *       <dt><tt>mode=[default,rec,sim,train,custom]</tt></dt>
 *       <dd>Set the AliROOT mode.  If not specified <tt>default</tt> 
 *         is assumed</tt>.  See also CreateAliROOTPar</dd>
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
   * @param opts Options 
   */
  GridHelper(const TUrl& url, Int_t verbose)
    : PluginHelper(url, verbose)
  {
    fOptions.Add("oper", "FULL|TERMINATE|SUBMIT", "Analysis operation", "FULL");
    fOptions.Add("split",  "N", "Maximum number of files before split", "max");
    fOptions.Add("merge",  "N", "Maximum number of files for merge", "max");
    fOptions.Add("run",    "RUNS", "Range, list, and/or file of runs", "");
    fOptions.Add("pattern","GLOB", "File/directory name pattern", "");
    fOptions.Add("mc", "Assume MC input");
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
    TString     runs   = fOptions.Get("runs");
    TObjArray*  tokens = runs.Tokenize(",");
    TObjString* part   = 0;
    TIter       next(tokens);
    Bool_t      range  = false;
    Bool_t      individual = false;
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
	fHandler->SetRunRange(first, last);
	ranges->Delete();
	range = true;
	continue;
      }
      if (s.IsDigit()) { // single run
	if (range) { 
	  Warning("GridHelper::RegisterRuns", 
		  "Run ranges and individual run specs do not mix, "
		  "ignoring %s", s.Data());
	  continue;
	}
	fHandler->AddRunNumber(s.Data());
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
      std::ifstream in(s.Data());
      if (!in) { 
	Warning("GridHelper::RegisterRuns", "Failed to open %s", s.Data());
	continue;
      }
      while (!in.eof()) { 
	Int_t r;
	in >> r;
	fHandler->AddRunNumber(r);
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

    // --- Open a connection to the grid -----------------------------
    TGrid::Connect(Form("%s://", fUrl.GetProtocol()));
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
    if (!PluginHelper::PostSetup()) return false;

    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    TString name(mgr->GetName());

    // --- Set the operation to do (TEST, SUBMIT, TERMINATE, FULL) ---
    TString operation("FULL");
    if (fOptions.Has("oper")) operation = fOptions.Get("oper");
    fHandler->SetRunMode(operation);

    // --- Do not test copying ---------------------------------------
    fHandler->SetCheckCopy(false);
    
    // --- Set output to be per run ----------------------------------
    fHandler->SetOutputToRunNo(true); 

    // --- Set the job tag -------------------------------------------
    fHandler->SetJobTag(name);

    // --- Set number of test files - used in test mode only ---------
    fHandler->SetNtestFiles(1);

    // --- Set the Time-To-Live --------------------------------------
    fHandler->SetTTL(70000);
    
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
	fHandler->SetSplitMaxInputFileNumber(fOptions.Get("split").Atoi());
      }
    }
    
    // --- Enable default outputs ------------------------------------
    fHandler->SetDefaultOutputs(true);

    // --- Merge parameters ------------------------------------------
    if (fOptions.Has("merge")) { 
      if (!fOptions.Get("merge").EqualTo("max")) { 
	fHandler->SetMaxMergeFiles(fOptions.Get("merge").Atoi());
      }
    }
    fHandler->SetMergeExcludes("AliAOD.root *EventStat*.root "
			       "*event_stat*.root");

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
    fHandler->SetRunPrefix(mgr->GetMCtruthEventHandler() ? "" : "000");

    // --- Add the run numbers ---------------------------------------
    Int_t nRun = RegisterRuns();

    // --- Set number of runs per master - set to one to per run -----
    fHandler->SetNrunsPerMaster(fOptions.Has("run-merge") ? 1 : nRun+1);

    // --- Loop over defined containers in the analysis manager, and -
    // --- declare these as outputs
    TString listOfAODs  = "";
    TString listOfHists = "";

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

    Int_t nReplica = 2;
    TString outArchive = Form("stderr, stdout@disk=%d", nReplica);
    if (!listOfHists.IsNull()) 
      outArchive.Append(Form(" hist_archive.zip:%s@disk=%d", 
			     listOfHists.Data(), nReplica));
    if (!listOfAODs.IsNull()) 
      outArchive.Append(Form(" aod_archive.zip:%s@disk=%d", 
			     listOfAODs.Data(), nReplica));
    if (listOfAODs.IsNull() && listOfHists.IsNull()) 
      Fatal("PostSetup", "No outputs defined");
    
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
    
    return mgr->StartAnalysis("grid", nEvents);
  }
  /** 
   * Link an auxilary file to working directory 
   * 
   * @param name Name of the file
   * 
   * @return true on success
   */
  virtual Bool_t AuxFile(const TString& name)
  {
    if (!Helper::AuxFile(name)) return false;
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
};
#endif
//
// EOF
//
