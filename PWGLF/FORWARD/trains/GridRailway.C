/**
 * @file   GridRailway.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 19:01:27 2012
 * 
 * @brief  Grid Analysis Railway
 * 
 * @ingroup pwglf_forward_trains_helper
 * 
 */
#ifndef GRIDHELPER_C
#define GRIDHELPER_C
#include "PluginRailway.C"
#ifndef __CINT__
# include <TUrl.h>
# include <TString.h>
# include <TGrid.h>
# include <AliAnalysisManager.h>
# include <AliAnalysisAlien.h>
# include <sstream>
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
struct GridRailway : public PluginRailway
{
  /** 
   * Constructor 
   * 
   * @param url  Url 
   * @param verbose Verbosity level
   */
  GridRailway(const TUrl& url, Int_t verbose)
    : PluginRailway(url, verbose), fRuns()
  {
    // Note, split, merge, and ttl are by default set to values
    // optimized for AOD production on real PbPb data.
    //
    // TTL shouldn't be much smaller than 4h10m.  Split and merge
    // shouldn't be much larger than 75, but probably not smaller than
    // 50.
    fOptions.Add("oper", "FULL|TERMINATE|SUBMIT", "Analysis operation", "FULL");
    fOptions.Add("split",  "N|max",  "Max number of files before split","50");
    fOptions.Add("merge",  "N|max",  "Max number of files for merge",   "50");
    fOptions.Add("run",    "RUNS",   "Range, list, and/or file of runs", "");
    fOptions.Add("alien",  "VERSION","Alien API version",              "V1.1x");
    fOptions.Add("ttl",    "N|max",  "Time to live",                   "6h");
    fOptions.Add("pattern","GLOB",   "File/directory name pattern", "");
    fOptions.Add("concat", "Concatenate all runs");
    fOptions.Add("exclude", "GLOB","Comma separated list of merge excludes","");
    fOptions.Add("files",   "FILE", "file containing list of files", "");
  }
  GridRailway(const GridRailway& o)
    : PluginRailway(o), fRuns()
  {}
  GridRailway& operator=(const GridRailway& o)
  {
    if (&o == this) return *this;
    PluginRailway::operator=(o);
    return *this;
  }
  virtual ~GridRailway() {}
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
      Error("GridRailway::RegisterRuns", "No runs specified");
      return -1;
    }
    Int_t       nRuns  = 0;
    TString     runs   = fOptions.Get("run");
    TObjArray*  tokens = runs.Tokenize(",+:");
    TObjString* part   = 0;
    TIter       next(tokens);
    Bool_t      range  = false;
    Bool_t      individual = false;
    // Info("GridRailway::RegisterRuns", "Runs specified are %s", runs.Data());
    while ((part = static_cast<TObjString*>(next()))) {
      TString& s = part->String();
      if (s.Contains("-")) { // Run range 
	if (range) { 
	  Warning("GridRailway::RegisterRuns", "Run range already specified, "
		  "ignoring %s", s.Data());
	  continue;
	}
	if (individual) { 
	  Warning("GridRailway::RegisterRuns", 
		  "Run ranges and individual run specs do not mix, "
		  "ignoring %s", s.Data());
	  continue;
	}
	TObjArray* ranges = s.Tokenize("-");
	if (ranges->GetEntriesFast() > 2) { 
	  Warning("GridRailway::RegisterRuns", "Invalid run range: %s", 
		  s.Data());
	  ranges->Delete();
	  continue;
	}
	Int_t first = static_cast<TObjString*>(ranges->At(0))->String().Atoi();
	Int_t last  = static_cast<TObjString*>(ranges->At(1))->String().Atoi();
	nRuns       = last-first+1;
	// Info("GridRailway::RegisterRuns", "Run range %d -> %d", first, last);
	fHandler->SetRunRange(first, last);
	ranges->Delete();
	range = true;
	for (Int_t r = first; r <= last; r++) StoreRun(r);
	continue;
      }
      if (s.IsDigit()) { // single run
	if (range) { 
	  Warning("GridRailway::RegisterRuns", 
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
	Warning("GridRailway::RegisterRuns", "Run ranges and list file "
		"do not mix, ignoring %s", s.Data());
	continue;
      }

      // We assume this part is a file 
      // Info("GridRailway::RegisterRuns", "Reading runs from %s", s.Data());
      std::ifstream in(s.Data());
      if (!in) { 
	s.Prepend("../");
	in.open(s.Data());
	if (!in) {
	  Warning("GridRailway::RegisterRuns", "Failed to open %s", s.Data());
	  continue;
	}
      }
      while (!in.eof()) {
	TString lne;
	lne.ReadLine(in);

	TString bare = lne.Strip(TString::kBoth);
	if (bare[0] == '#') continue;

	TObjArray* ltokens = bare.Tokenize(" \t,");
	TIter lnext(ltokens);
	TObjString* str = 0;
	while ((str = static_cast<TObjString*>(lnext()))) {
	  const TString& token = str->String();
	  if (!token.IsDigit()) continue;
	  
	  int r = token.Atoi();
	  fHandler->AddRunNumber(r);
	  StoreRun(r);
	  nRuns++;
	}
	ltokens->Delete();
      }
#if 0
      while (!in.eof()) { 
	Int_t r;
	in >> r;
	// Info("GridRailway::RegisterRuns", "Read %d, adding", r);
	fHandler->AddRunNumber(r);
	StoreRun(r);
	nRuns++;
	Char_t c;
	in >> c;
	if (in.bad()) break;
      }
#endif
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
    if (!PluginRailway::PreSetup()) return false;
    
    // --- Add system library dir to load path -----------------------
    gSystem->AddDynamicPath("/usr/lib");

    // --- Open a connection to the grid -----------------------------
    if (!TGrid::Connect(Form("%s://", fUrl.GetProtocol()))) { 
      Error("GridRailway::PreSetup", "Failed to connect to AliEN");
      return false;
    }
    if (!gGrid || !gGrid->IsConnected()) { 
      Error("GridRailway::PreSetup", "Failed to connect to AliEN");
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
    // Info("GridRailway::PostSetup", "Calling super.PostSetup");
    if (!PluginRailway::PostSetup()) return false;

    // --- API version -----------------------------------------------
    fHandler->SetAPIVersion(fOptions.Get("alien"));
    
    // --- Get the name ----------------------------------------------
    // Info("GridRailway", "Proceeding with plugin setup");
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    TString name(mgr->GetName());

    // --- Set the operation to do (TEST, SUBMIT, TERMINATE, FULL) ---
    TString operation("FULL");
    if (fOptions.Has("oper")) operation = fOptions.Get("oper");
    fHandler->SetRunMode(operation);

    // --- Add the run numbers ---------------------------------------
    fHandler->SetRunPrefix(fOptions.Has("mc") ? "%d" : "%09d");
    Int_t nRun = RegisterRuns();

    // --- Possibly add files from input file ------------------------
    AddFiles();
    
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
      TString sttl = fOptions.Get("ttl");
      if (!sttl.EqualTo("max")) {
	Int_t   ttl  = 0;
	if (sttl.IsDigit()) ttl = sttl.Atoi();
	else { 
	  // Parse string of the form <DAYS>d<HOURS>h<MINUTES>m<SECONDS>s
	  Int_t id = sttl.Index("d", 0);
	  if (id == kNPOS) id = -1;
	  else { 
	    TString sdays(sttl(0,id));
	    ttl += 24 * 60 * 60 * sdays.Atoi();
	  }
	  Int_t ih = sttl.Index("h", id+1);
	  if (ih == kNPOS) ih = id;
	  else { 
	    TString shour(sttl(id+1,ih-id-1));
	    ttl += 60 * 60 * shour.Atoi();
	  }
	  Int_t im = sttl.Index("m", ih+1);
	  if (im == kNPOS) im = ih;
	  else { 
	    TString smin(sttl(ih+1, im-ih-1));
	    ttl += 60 * smin.Atoi();
	  }
	  Int_t is = sttl.Index("s", im+1);
	  if (is != kNPOS) { 
	    TString ssec(sttl(im+1, is-im-1));
	    ttl += ssec.Atoi();
	  }
	}
	if (ttl != 0) fHandler->SetTTL(ttl);
	else 
	  Warning("", "Option ttl given but no value found");
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
    TString exclude="AliAOD.root *EventStat*.root *event_stat*.root";
    if (fOptions.Has("exclude")) { 
      TString exOpt = fOptions.Get("exclude");
      exOpt.ReplaceAll(",", " ");
      exclude.Append(" ");
      exclude.Append(exOpt);
    }
    fHandler->SetMergeExcludes(exclude);
    
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
	Warning("GridRailway::PreSetup", "No tree name specified, assuming T");
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
  void ScanFiles()
  {
    // Check if runs where registered, and if so that the first run
    // registered isn't 0
    if (fRuns.GetEntries() > 0 && 
	fRuns.At(0)->GetUniqueID() != 0) return;

    TString path    = fUrl.GetFile();
    TString pattern = fOptions.Get("pattern");

    if (path.IsNull() || pattern.IsNull()) {
      Warning("ScanFiles", "No search path (%s) or pattern (%s) specified",
	      path.Data(), pattern.Data());
      return;
    }

    TString cmd = Form("alien_find %s %s", path.Data(), pattern.Data());
    TString ret = gSystem->GetFromPipe(cmd);
    if (ret.IsNull()) {
      Warning("ScanFiles", "Command %s failed", cmd.Data());
      return;
    }

    std::stringstream str(ret.Data());
    AddFiles(str);
    
  }
  void AddFiles()
  {
    TString files = fOptions.Get("files");
    // Info("AddFiles", "Getting list of files from '%s'", files.Data());
    if (files.IsNull()) {
      ScanFiles();
      return;
    }
    
    std::ifstream in(files.Data());
    if (!in) {
      Warning("", "Failed to open the file %s", files.Data());
      in.open(Form("../%s", files.Data()));
      if (!in) {
	Warning("", "Failed to open the file ../%s - giving up", files.Data());
	return;
      }
    }
    AddFiles(in);
    in.close();
  }
  void AddFiles(std::istream& in)
  {
    do {
      TString l;
      l.ReadLine(in);
      
      TString tmp = l.Strip(TString::kBoth, ' ');
      if (!tmp.EndsWith(".root")) {
	// Info("AddFiles", "'%s' is not a ROOT file!", l.Data());
	continue;
      }
      // Info("AddFiles", "Adding %s to list of inputs", tmp.Data());
      fHandler->AddDataFile(tmp);
    } while (!in.eof());
  }
  
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
    if (nEvents < 0)  nEvents = 123456789;
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
  virtual Bool_t AuxFile(TString& name, bool copy=false)
  {
    TString local = name;
    if (!Railway::AuxFile(local, copy)) return false;
    // We need to add this file as an additional 'library', so that the 
    // file is uploaded to the users Grid working directory. 
    fHandler->AddAdditionalLibrary(local);
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
      Warning("GridRailway::OutputLocation", "No AliEn handler");
      return ret;
    }
    ret = fHandler->GetGridOutputDir();
    if (ret.BeginsWith("/")) return ret;

    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) { 
      Warning("GridRailway::OutputLocation", "No analysis manager");
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
    TString macDir("$ALICE_PHYSICS/PWGLF/FORWARD/trains");
    std::ofstream t("Terminate.C");
    if (!t) { 
      Error("GridRailway::AuxSave", "Failed to make terminate ROOT script");
      return;
    }

    t << "// Generated by GridRailway\n"
      << "Bool_t Terminate(Bool_t localMerge=false)\n"
      << "{\n"
      << "  TString name = \"" << escaped << "\";\n"
      << "  TString libs = \"" << libs << "\";\n"
      << "  TString pars = \"" << pars << "\";\n"
      << "  TString srcs = \"" << srcs << "\";\n\n"
      << "  gSystem->Load(\"libANALYSIS\");\n"
      << "  gSystem->Load(\"libANALYSISalice\");\n"
      << "  gSystem->AddIncludePath(\"-I$ALICE_ROOT/include\");\n\n"
      << "  gSystem->AddIncludePath(\"-I$ALICE_PHYSICS/include\");\n\n"
      << "  gROOT->LoadMacro(\"" << macDir << "/GridTerminate.C+g\");\n\n"
      << "  return GridTerminate(name,libs,pars,srcs,localMerge);\n"
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
      Error("GridRailway::AuxSave", "Failed to make ROOT script Download.C");
      return;
    }
    d << "// Generated by GridRailway\n"
      << "void Download(Bool_t unpack=true)\n"
      << "{\n"
      << "  TString base = \"" << fUrl.GetProtocol() << "://" 
      << OutputPath() << "\";\n"
      << "  TString runs = \"" << runs << "\";\n\n"
      << "  gROOT->LoadMacro(\"" << macDir << "/GridDownload.C\");\n\n"
      << "  GridDownload(base, runs, unpack);\n"
      << "}\n"
      << "// EOF\n"
      << std::endl;
    d.close();

    std::ofstream w("Watch.C");
    if (!w) {
      Error("GridRailway::AuxSave", "Failed to make ROOT script Watch.C");
      return;
    }
    w << "// Generated by GridRailway\n"
      << "void Watch(Bool_t batch=false, Int_t delay=5*60)\n"
      << "{\n"
      << "  TString name = \"" << escaped << "\";\n"
      << "  gROOT->LoadMacro(\"" << macDir << "/GridWatch.C+g\");\n\n"
      << "  GridWatch(name,batch,delay);\n"
      << "}\n"
      << "// EOF\n"
      << std::endl;
    w.close();

    if (!AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler())
      // No AODs generated - no need for DownloadAODs.C script 
      return; 

    std::ofstream a("DownloadAODs.C");
    if (!a) { 
      Error("SaveDownloadAODs", "Failed to open DownloadAODs.C");
      return;
    }
    a << "// Generated by GridRailway\n"
      << "void DownloadAODs(Bool_t force=false)\n"
      << "{\n"
      << "  if (!TGrid::Connect(\"alien://\")) {\n"
      << "    Error(\"DownloadAODs\",\"Failed to connect to AliEn\");\n"
      << "    return;\n"
      << "  }\n\n"
      << "  TString dir(\"" << OutputPath() << "\");\n"
      << "  TString pat(\"*/AliAOD.root\");\n"
      << "  TGridResult* r = gGrid->Query(dir,pat);\n"
      << "  if (!r) {\n"
      << "    Error(\"DownloadAODs\",\"No result from query\");\n"
      << "    return;\n"
      << "  }\n\n"
      << "  Int_t n = r->GetEntries();\n"
      << "  Printf(\"=== Got a total of %d AOD files\",n);\n"
      << "  for (Int_t i = 0; i < n; i++) {\n"
      << "     TString path(r->GetKey(i, \"turl\"));\n"
      << "     TString dir(gSystem->DirName(path));\n"
      << "     TString sub(gSystem->BaseName(dir));\n"
      << "     TString subsub(gSystem->BaseName(gSystem->DirName(dir)));\n"
      << "     TString out = TString::Format(\"AliAOD_%s_%s.root\",\n"
      << "                                   subsub.Data(),sub.Data());\n"
      << "     if (!gSystem->AccessPathName(out.Data()) && !force) {\n"
      << "       Printf(\"=== Already have %s\",out.Data());\n"
      << "       continue;\n"
      << "     }\n"
      << "     Printf(\"=== Getting %s %s (%3d/%3d)\",\n"
      << "            subsub.Data(),sub.Data(),i,n);\n"
      << "     if (!TFile::Cp(path, out)) {\n"
      << "       Warning(\"DownloadAODs\",\"Failed to copy %s -> %s\",\n"
      << "               path.Data(), out.Data());\n"
      << "       continue;\n"
      << "     }\n"
      << "   }\n"
      << "}\n"
      << "// EOF\n"
      << std::endl;
    a.close();
    
  }
  TList fRuns;
};
#endif
//
// EOF
//
