/**
 * @defgroup pwglf_forward_trains_helper Analysis helpers
 *
 * @ingroup pwglf_forward_trains
 * 
 * 
 */
/**
 * @file   Railway.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 19:00:17 2012
 * 
 * @brief  Base class for analysis helpers
 * 
 * @ingroup pwglf_forward_trains_helper
 * 
 */
#ifndef TRAIN_HELPER_C
#define TRAIN_HELPER_C
#ifndef __CINT__
# include "Option.C"
# include "ChainBuilder.C"
# include <TUrl.h>
# include <TString.h>
# include <TMap.h>
# include <TObjString.h>
# include <TSystem.h>
# include <TROOT.h>
# include <TError.h>
# include <TObjArray.h>
# include <TFile.h>
# include <TChain.h>
# include <AliAnalysisManager.h>
# include <iostream>
#else 
class TString;
class TUrl;
class TMap;
class Option;
class OptionList;
class TChain;
#endif

/**
 * Helper class to set-up an analysis using a train 
 *
 * @ingroup pwglf_forward_trains_helper
 */
struct Railway 
{
  enum EMode { 
    kLocal, 
    kProof, 
    kGrid
  };
  enum EOperation { 
    kTest,
    kOffline,
    kSubmit,
    kTerminate, 
    kFull
  };
  enum EInput { 
    kESD, 
    kAOD, 
    kUser
  };
  Railway(const Railway& o) 
    : fUrl(o.fUrl), fOptions(o.fOptions), fVerbose(o.fVerbose)
  {}
  Railway& operator=(const Railway&) { return *this; }
  /** 
   * @{ 
   * @name Create a railway 
   */
  /** 
   * Create a helper object. 
   *
   * @param verbose Verbosity 
   * @param url Url describing the job. 
   *
   * - Local mode: 
   *
   * @code
   * local://<path>[#<treeName>][?[recursive[&]]]
   * @endcode 
   *
   * &lt;path&gt; can be a single ROOT file, a text file with one
   * entry per file to add to the chain, or a directory containing the
   * files. If &lt;path&gt; does not start with a '/' then it is
   * interpreted as a relative path.
   *
   * - Grid mode: 
   *
   * @code 
   * alien:///<path>#<pattern>
   * @endcode
   *
   * - PROOF mode: 
   * 
   * Several options 
   *
   * @code 
   * lite:///<path>[?[recursive[&]][workers=<n>]][#treeName]
   * proof:///<path>[?[recursive[&]][workers=<n>]][#treeName]
   * @endcode 
   *
   * @code 
   * proof://<host>/<dsname>[?[workers=<n>[&]][dsname[=<outname>]]][#treeName]
   * @endcode 
   *
   * Note, if &lt;host&gt; is recognised as an Alice Analysis
   * Facility, then the Grid handler (AliAnalysisAlien) is used unless
   * the option <tt>plain</tt> was given. 
   * 
   * @return Newly allocated helper or null 
   */
  static Railway* Create(const TUrl& url, Int_t verbose=0);
  /** 
   * Create an instance of a helper class 
   */
  static Railway* CreateObject(const TString& cl, 
			      const TUrl&    url, 
			      Int_t verbose=0)
  {
    if (verbose < 3) gSystem->RedirectOutput("/dev/null","w");
    if (cl.Contains("proof", TString::kIgnoreCase) || 
	cl.Contains("vaf",   TString::kIgnoreCase) || 
	cl.Contains("lite",  TString::kIgnoreCase) || 
	cl.Contains("aaf",   TString::kIgnoreCase)) {
      gSystem->Load("libProof");
      gSystem->Load("libProofPlayer");
    }
    // (Always) recompile and with debug symbols 
    gROOT->LoadMacro(Form("%s.C++g",cl.Data()));
    Long_t ptr = gROOT->ProcessLine(Form("new %s(\"%s\", %d);", 
					 cl.Data(), url.GetUrl(), verbose));
    if (verbose < 3) gSystem->RedirectOutput(0);
    if (!ptr) { 
      Warning("Railway::CreateObject", "Failed to instantize a %s", cl.Data());
      return 0;
    }
    Railway* h = reinterpret_cast<Railway*>(ptr);
    return h;
  }
  /** 
   * Show help on URL using the interpreter 
   * 
   * @param cl Railway class 
   */
  static void ShowUrlHelp(const TString& cl)
  {
    Railway* h = CreateObject(cl, "", true);
    if (!h) return;

    std::cout << "   " << h->UrlHelp() << std::endl;
  }
  /** 
   * Show help on URL and options using the interpreter 
   * 
   * @param cl Railway class 
   */
  static void ShowFullHelp(const TString& cl) 
  {
    Railway* h = CreateObject(cl, "", true);
    if (!h) return;
    
    std::cout << h->Desc() << ":\n" 
	      << "==============\n"
	      << "  " << h->UrlHelp() << "\n\n"
	      << "Options: " << std::endl;
    h->Options().Help(std::cout);
    std::cout << std::endl;
  }
  /* @} */
  /** 
   * @{ 
   * @name Loading stuff 
   */
  /** 
   * Set whether to use pars.  On return, the argument is set to the
   * old value.  So to temporaruly turn off pars, do
   * 
   * @code 
   Bool_t usePar = false;
   fRailway->UsePar(usePar); // usePar is now old value 
   // Load real libraries 
   fRailway->UsePar(usePar); // Restore old value 
   * @endcode 
   * 
   * @param use Whether to use pars or not.  On return contains old value 
   */
  virtual void UsePar(Bool_t& use) {}
  /** 
   * Add an include path
   * 
   * @param path The path to add for header search 
   * 
   * @return true on success
   */
  virtual Bool_t AddIncludePath(const TString& path)
  {
    TString p(path);
    if (!p.BeginsWith("-")) {
      // IF the path does not begin with a -, we assume its a path
      p.Prepend("-I");
    }
    gSystem->AddIncludePath(p);
    return true;
  }
  /** 
   * Load a library 
   * 
   * @param name   Name of library 
   * @param slave  If true also load on slaves
   * @param forcePar if true, force load as PAR 
   * 
   * @return true on success, false otherwise 
   */
  virtual Bool_t LoadLibrary(const TString& name, 
			     Bool_t slave=true,
			     Bool_t forcePar=false) = 0;
  /** 
   * Load a source file, and compile it 
   * 
   * @param name Name of the source file 
   * @param copy Copy rather than link 
   * 
   * @return true on success
   */
  virtual Bool_t LoadSource(const TString& name, bool copy=false)
  {
    TString local = name;
    if (!AuxFile(local, copy)) {
      Warning("LoadSource", "Failed to add aux source file %s", name.Data());
      return false;
    }
    // TString base(gSystem->BaseName(name));
    Info("LoadSource", "Building %s", local.Data());
    // gROOT->ProcessLine("gSystem->RedirectOutput(\"build.log\",\"a\");");
    gROOT->LoadMacro(Form("%s+g", local.Data()));
    // gROOT->ProcessLine("gSystem->RedirectOutput(0);");
    Info("LoadSource", "End of loading source");
    return true;
  }
  /** 
   * Load auxillary file - not compiled or sourced.  Just copied to
   * working directory
   * 
   * @param name Extra file name 
   * @param copy Copy rather than link 
   *
   * @return true on success 
   */
  virtual Bool_t LoadAux(const TString& name, Bool_t copy=false)
  {
    TString local = name;
    if (!AuxFile(local, copy)) return false;
    return true;
  }

  /** 
   * Load needed ROOT libaries 
   */
  virtual Bool_t LoadROOT()
  {
    if (gSystem->Load("libTree")    < 0) return false;
    if (gSystem->Load("libGeom")    < 0) return false;
    if (gSystem->Load("libVMC")     < 0) return false;
    if (gSystem->Load("libPhysics") < 0) return false;
    if (gSystem->Load("libMinuit")  < 0) return false;
    return true;
  }
  /** 
   * Set-up to load the AliROOT libraries 
   * 
   * @return true on success
   */
  virtual Bool_t LoadAliROOT()
  {
    if (!gSystem->Getenv("ALICE_ROOT")) { 
      Error("Railway::LoadAliROOT", "Local AliROOT not available");
      return false;
    }
    if (!LoadLibrary("STEERBase"))     return false;
    if (!LoadLibrary("ESD"))           return false;
    if (!LoadLibrary("AOD"))           return false;
    if (!LoadLibrary("ANALYSIS"))      return false;
    if (!LoadLibrary("OADB"))          return false;
    if (!LoadLibrary("ANALYSISalice")) return false;
    return true;
  }
  /** 
   * Set-up to load the AliROOT libraries 
   * 
   * @return true on success
   */
  virtual Bool_t LoadAliPhysics()
  {
    if (!gSystem->Getenv("ALICE_PHYSICS")) { 
      Error("Railway::LoadAliPhysics", "Local AliPhysics not available");
      return false;
    }
    if (!LoadLibrary("OADB"))          return false;
    return true;
  }
  /* @} */
  /** 
   * @{ 
   * @name Mode of exection 
   */
  /** 
   * Get the execution mode 
   * 
   * @return Execution mode set in set-up URL
   */
  virtual UShort_t Mode() const = 0;
  /**
   * Get the mode string used for AliAnalysisManager::StartAnalysis
   */
  virtual const char* ModeString() const { return "unknown"; }
  /** 
   * Get the operation - this only makes sense for Grid jobs
   * 
   * @return Operation type
   */
  virtual UShort_t Operation() const { return kFull; }
  /* @} */
  /** 
   * @{ 
   * @name Input/Output
   */
  /** 
   * Get the input data type 
   *
   * @return Input data type 
   */
  virtual Short_t InputType() const
  {
    UShort_t ret = DeduceType(fUrl.GetAnchor());
    if (ret != kUser) return ret;

    if (fOptions.Has("pattern")) ret = DeduceType(fOptions.Get("pattern"));
    if (ret != kUser) return ret;
    
    ret = DeduceType(fUrl.GetFile());
    return ret;
  }
  /** 
   * Check if the MC option was set
   * 
   * @return true if the MC option was given 
   */
  virtual Bool_t IsMC() const { return fOptions.Has("mc"); }
  /** 
   * The file part of tehe output URL - overwritten by derived classes. 
   * 
   * 
   * @return File part of output URL
   */
  virtual TString OutputPath() const { return ""; }
  /** 
   * Get the location of the output data.  Use ful to define second pass 
   * scripts, etc. 
   * 
   * 
   * @return Url string 
   */
  virtual TString OutputLocation() const 
  {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr || !mgr->GetOutputEventHandler()) return "";

    TString path(OutputPath());
    if (path.IsNull()) {
      path = gSystem->WorkingDirectory(); 
      // mgr->GetName();
    }

    TUrl u(fUrl);
    u.SetFile(path);
    u.SetAnchor("aodTree");
    TString opt(u.GetOptions());
    // if (opt.Contains("AliESDs")) 
    opt.ReplaceAll("AliESDs", "AliAOD");
    // else                         opt.Append("&pattern=AliAOD*.root");
    u.SetOptions(opt);

    return u.GetUrl();
  }
  /* @} */
  /** 
   * @{ 
   * @name Processing 
   */
  /** 
   * Set-up done before task setup 
   * 
   * @return true on success
   */
  virtual Bool_t PreSetup() = 0;
  /** 
   * Set-up done after task setup 
   * 
   * @return true on success
   */
  virtual Bool_t PostSetup() = 0;
  /** 
   * Run the analysis
   * 
   * @param nEvents Number of events to analyse 
   * 
   * @return The return value of AliAnalysisManager::StartAnalysis
   */
  virtual Long64_t Run(Long64_t nEvents=-1) = 0;
  /**
   * Add a monitor object - only for PROOF 
   */
  virtual Bool_t AddMonitor(const TString&) { return true; }
  /* @} */
  /** 
   * @{ 
   * @name Various help functions
   */
  /** 
   * Print information to standard output 
   * 
   * @return 
   */
  virtual void Print(Option_t* ="") const
  {
    std::cout << "Url: " << fUrl.GetUrl() << std::endl;
    fOptions.Show(std::cout);
  }
  /** 
   * @return URL help string 
   */
  virtual const Char_t* UrlHelp() const  = 0;
  /** 
   * @return Short description string 
   */
  virtual const char* Desc() const = 0;
  /** 
   * Get the input URL 
   * 
   * @return Input URL 
   */
  const TUrl& Url() const { return fUrl; }
  /** 
   * Get the list of options 
   * 
   * @return Reference to list of options 
   */
  const OptionList& Options() const { return fOptions; }
  /** 
   * Write auxillary ROOT (and possible shell) script for more 
   * (post-)processing e.g., terminate
   * 
   */
  virtual void AuxSave(const TString& /*escaped*/, 
		       Bool_t /*asShellScript*/) {}
protected:
  /** 
   * Constructor 
   * 
   * @param url Set-up URL
   * @param verbose Verbosity level 
   */
  Railway(const TUrl& url, Int_t verbose) 
    : fUrl(url), fOptions(), fVerbose(verbose)
  {
    fOptions.Add("mc", "Assume simulation input");
  }

  virtual Bool_t ParseOptions()
  {
    //std::cout << "Url options: \"" << fUrl.GetOptions() << "\"" << std::endl;
    return fOptions.Parse(fUrl.GetOptions(), "&");
  }
  /** 
   * Destructor 
   */
  virtual ~Railway()
  {
  }
  /** 
   * @{ 
   * @name Loading things 
   */
  /** 
   * Normalize a library name
   * 
   * @param name 
   * 
   * @return 
   */
  const TString& MakeLibraryName(const TString& name)
  {
    static TString libName;

    libName = name;

    if (!libName.BeginsWith("lib")) { 
      // Check if the library corresponds to a compiled macro 
      if (!gSystem->AccessPathName(Form("%s_C.so", libName.Data()))) 
	libName.Append("_C");
      else if (!gSystem->AccessPathName(Form("../%s_C.so", libName.Data()))) 
	libName = Form("../%s_C", libName.Data());
      else 
	libName = Form("lib%s", libName.Data());
    }
    if (!libName.EndsWith(".so"))   libName.Append(".so");

    return libName;
  }
  /** 
   * Link an auxilary file to working directory 
   * 
   * @param name Name of the file
   * @param copy Copy rather than link
   *
   * @return true on success
   */
  virtual Bool_t AuxFile(TString& name, bool copy=false)
  {
    TString path(gSystem->ExpandPathName(name.Data()));
    TString search(gSystem->GetIncludePath());
    search.ReplaceAll("-I", ":");
    search.ReplaceAll(" ", "");
    search.Append(":");
    search.Append(gROOT->GetMacroPath());    
    // If not absolute, prepend up-one
    if (!path.BeginsWith("/")) search.Prepend("../:");
    
    // Remove any trailing AcLic instructions 
    if (path.EndsWith("+") || path.EndsWith("+g")) {
      do {
	char c = path[path.Length()-1];
	if (c != '+' && c != 'g') break;
	path.Remove(path.Length()-1);
      } while (path.Length());
    }
    TString query = gSystem->Which(search, path);
    if (gSystem->AccessPathName(query.Data())) { 
      // File not accessible
      Warning("Railway::AuxFile", "File %s not accessible", path.Data());
      return false;
    }
    TString base(gSystem->BaseName(query.Data()));
    name = base;
    if (gSystem->AccessPathName(base.Data()) == 0) { 
      // File or link exists - remove it 
      if (gSystem->Unlink(base) != 0) { 
	Error("Railway::AuxFile", "Failed to remove old %s", base.Data());
	return false;
      }
    }
    if (copy) 
      TFile::Cp(path, base);
    else 
      gSystem->Exec(Form("ln -s %s .", query.Data()));
    return true;
  }
  /* @} */
  /** 
   * @{ 
   * @name Input/output 
   */
  /** 
   * Create a local chain based on URL and options 
   * 
   * 
   * @return Chain or null
   */
  TChain* LocalChain()
  {
    // -- Check the source -------------------------------------------
    TString  src       = fUrl.GetFile();
    if (src.IsNull()) {
      Error("LocalChain", "No input source specified");
      return 0;
    }

    // --- Check possible pattern ------------------------------------
    TString  pattern   = (fOptions.Has("pattern") ?fOptions.Get("pattern") :"");
    pattern.ReplaceAll("@", "#");

    // --- Get the tree name -----------------------------------------
    TString  treeName  = fUrl.GetAnchor();

    // --- Create flags for the chain builder ------------------------
    UShort_t flags     = 0;
    if (fOptions.Has("mc") &&
	AliAnalysisManager::GetAnalysisManager()
	->GetMCtruthEventHandler() != 0) flags |= ChainBuilder::kMC;
    if (fOptions.Has("recursive"))       flags |= ChainBuilder::kRecursive;
    if (fOptions.Has("trackref"))        flags |= ChainBuilder::kTrRef;
    if (fOptions.Has("clean"))           flags |= ChainBuilder::kClean;
    if (fOptions.Has("scan"))            flags |= ChainBuilder::kScan;
    if (fVerbose > 5)                    flags |= ChainBuilder::kVerbose;

    // --- Check input -----------------------------------------------
    UShort_t type      = ChainBuilder::CheckSource(src, flags);
    if (type == ChainBuilder::kInvalid) {
      Error("LocalChain", "Cannot generate TChain from %s", src.Data());
      return 0;
    }

    // --- Create the chain ------------------------------------------
    TChain* chain = ChainBuilder::Create(type, src, treeName, pattern, flags);
    if (!chain) { 
      Error("LocalChain", "No chain defined "
	    "(src=%s, treeName=%s, pattern=%s, flags=0x%x)", 
	    src.Data(), treeName.Data(), pattern.Data(), flags);
      return 0;
    }
    
    return chain;
  }
  /* @} */
  /** 
   * @{ 
   * @name Helper to figure out what we're doing
   */
  /** 
   * Deduce the top of job from a string 
   * 
   * @param str String 
   * 
   * @return Job type
   */
  static UShort_t DeduceType(const TString& str)
  {
    if (str.IsNull()) return kUser;
    if (str.Contains("aod", TString::kIgnoreCase)) return kAOD;
    if (str.Contains("esd", TString::kIgnoreCase)) return kESD;
    return kUser;
  }
  /* @} */
  // --- Data members ------------------------------------------------
  TUrl        fUrl;     // The URI
  OptionList  fOptions; 
  Int_t       fVerbose;
};




// ===================================================================
Railway* 
Railway::Create(const TUrl& url, Int_t verbose)
{
  if (!url.IsValid()) { 
    Warning("Railway::Create", "URL is invalid");
    return 0;
  }

  TString prot(url.GetProtocol());
  prot.ToLower();

  TUrl tmp(url);
  TString opts(tmp.GetOptions());
  TString host(url.GetHost());
  TString cl = "";
  if (prot.EqualTo("alien")) { 
    // Create an AliEn helper 
    cl = "GridRailway";
  }
  else if (prot.EqualTo("local")) { 
    // Create Lite helper 
    cl = "LocalRailway";
  }
  else if (prot.EqualTo("proof")) { 
    // Create a Proof helper 
    if (host.IsNull()) 
      cl = "LiteRailway";
    else if (host.BeginsWith("alice-caf")) { 
      // AAF
      ::Warning("Railway::Create", "CAF has been decommissioned");
      cl = opts.Contains("plugin") ? "AAFPluginRailway" : "AAFRailway";
    }
    else if (host.BeginsWith("alivaf")) {
      // VAF
      cl = "VAFRailway";
    }
    else 
      cl = "ProofRailway";
  }
  else if (prot.EqualTo("lite")) { 
    // Create a Proof helper 
    cl = "LiteRailway";
  }
  else if (prot.EqualTo("help")) {
    // Special HELP protocol
    if (host.Contains("options")) {
      std::cout << "Possible URL types and options are:" << std::endl;
      ShowFullHelp("LocalRailway");
      ShowFullHelp("ProofRailway");
      ShowFullHelp("LiteRailway");
      ShowFullHelp("VAFRailway");
      ShowFullHelp("AAFRailway");
      ShowFullHelp("AAFPluginRailway");
      ShowFullHelp("GridRailway");
      return 0;
    }
    std::cout << "Possible URL types are:" << std::endl;
    ShowUrlHelp("LocalRailway");
    ShowUrlHelp("ProofRailway");
    ShowUrlHelp("LiteRailway");
    ShowUrlHelp("VAFRailway");
    ShowUrlHelp("AAFRailway");
    ShowUrlHelp("AAFPluginRailway");
    ShowUrlHelp("GridRailway");
    return 0;
  }
  // --- Check if we got a scheme ------------------------------------
  if (cl.IsNull()) {
    Error("Railway::Create", "Unknown scheme: %s", prot.Data());
    return 0;
  }

  // --- Use interpreter to make our object --------------------------
  Railway* helper = CreateObject(cl, url, verbose);
  if (!helper) {
    Error("Railway::Create", "Failed to make object of class %s", cl.Data());
    return 0;
  }

  // --- Parse options -----------------------------------------------
  if (!helper->ParseOptions()) {
    delete helper;
    helper = 0;
  }

  return helper;
}

#endif
