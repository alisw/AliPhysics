/**
 * @defgroup pwglf_forward_trains_helper Analysis helpers
 *
 * @ingroup pwglf_forward_trains
 * 
 * 
 */
/**
 * @file   Helper.C
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
# include <TUrl.h>
# include <TString.h>
# include <TMap.h>
# include <TObjString.h>
# include <TSystem.h>
# include <TROOT.h>
# include <TError.h>
# include <TObjArray.h>
# include <TFile.h>
# include <AliAnalysisManager.h>
# include <iostream>
#else 
class TString;
class TUrl;
class TMap;
class Option;
class OptionList;
#endif

/**
 * Helper class to set-up an analysis using a train 
 *
 * @ingroup pwglf_forward_trains_helper
 */
struct Helper 
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
  Helper(const Helper& o) 
    : fUrl(o.fUrl), fOptions(o.fOptions), fVerbose(o.fVerbose)
  {}
  Helper& operator=(const Helper&) { return *this; }
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
  static Helper* Create(const TUrl& url, Int_t verbose=0);
  /** 
   * Load a library 
   * 
   * @param name   Name of library 
   * @param slave  If true also load on slaves
   * 
   * @return true on success, false otherwise 
   */
  virtual Bool_t LoadLibrary(const TString& name, 
			     Bool_t slave=true) = 0;
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
    if (!AuxFile(name, copy)) return false;
    TString base(gSystem->BaseName(name));
    gROOT->ProcessLine("gSystem->RedirectOutput(\"build.log\",\"a\");");
    gROOT->LoadMacro(Form("%s++g", base.Data()));
    gROOT->ProcessLine("gSystem->RedirectOutput(0);");
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
    if (!AuxFile(name, copy)) return false;
    return true;
  }

  /** 
   * Load needed ROOT libaries 
   */
  virtual Bool_t LoadROOT()
  {
    if (gSystem->Load("libTree.so")    < 0) return false;
    if (gSystem->Load("libGeom.so")    < 0) return false;
    if (gSystem->Load("libVMC.so")     < 0) return false;
    if (gSystem->Load("libPhysics.so") < 0) return false;
    if (gSystem->Load("libMinuit.so")  < 0) return false;
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
      Error("Helper::LoadAliROOT", "Local AliROOT not available");
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
   * Get the execution mode 
   * 
   * @return Execution mode set in set-up URL
   */
  virtual UShort_t Mode() const = 0;
  /** 
   * Get the operation - this only makes sense for Grid jobs
   * 
   * @return Operation type
   */
  virtual UShort_t Operation() const { return kFull; }
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
    opt.ReplaceAll("AliESDs", "AliAOD");
    u.SetOptions(opt);

    return u.GetUrl();
  }
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
  /** 
   * Create an instance of a helper class 
   */
  static Helper* CreateObject(const TString& cl, 
			      const TUrl&    url, 
			      Int_t verbose=0)
  {
    if (verbose < 3) gSystem->RedirectOutput("/dev/null","w");
    if (cl.Contains("proof", TString::kIgnoreCase) || 
	cl.Contains("lite",  TString::kIgnoreCase) || 
	cl.Contains("aaf",   TString::kIgnoreCase)) {
      gSystem->Load("libProof");
      gSystem->Load("libProofPlayer");
    }
    // Always recompile and with debug symbols 
    gROOT->LoadMacro(Form("%s.C++g",cl.Data()));
    Long_t ptr = gROOT->ProcessLine(Form("new %s(\"%s\", %d);", 
					 cl.Data(), url.GetUrl(), verbose));
    if (verbose < 3) gSystem->RedirectOutput(0);
    if (!ptr) { 
      Warning("Helper::CreateObject", "Failed to instantize a %s", cl.Data());
      return 0;
    }
    Helper* h = reinterpret_cast<Helper*>(ptr);
    return h;
  }
  /** 
   * Show help on URL using the interpreter 
   * 
   * @param cl Helper class 
   */
  static void ShowUrlHelp(const TString& cl)
  {
    Helper* h = CreateObject(cl, "", true);
    if (!h) return;

    std::cout << "   " << h->UrlHelp() << std::endl;
  }
  /** 
   * Show help on URL and options using the interpreter 
   * 
   * @param cl Helper class 
   */
  static void ShowFullHelp(const TString& cl) 
  {
    Helper* h = CreateObject(cl, "", true);
    if (!h) return;
    
    std::cout << h->Desc() << ":\n" 
	      << "==============\n"
	      << "  " << h->UrlHelp() << "\n\n"
	      << "Options: " << std::endl;
    h->Options().Help(std::cout);
    std::cout << std::endl;
  }
protected:
  /** 
   * Constructor 
   * 
   * @param url Set-up URL
   * @param verbose Verbosity level 
   */
  Helper(const TUrl& url, Int_t verbose) 
    : fUrl(url), fOptions(), fVerbose(verbose)
  {
    fOptions.Add("mc", "Assume simulation input");
  }

  virtual Bool_t ParseOptions()
  {
    return fOptions.Parse(fUrl.GetOptions(), "&");
  }
  /** 
   * Destructor 
   */
  virtual ~Helper()
  {
  }
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
  virtual Bool_t AuxFile(const TString& name, bool copy=false)
  {
    TString path(gSystem->ExpandPathName(name.Data()));
    // If not absolute, prepend up-one
    if (!path.BeginsWith("/")) path.Prepend("../");
    if (gSystem->AccessPathName(path.Data())) { 
      // File not accessible
      Warning("Helper::AuxFile", "File %s not accessible", path.Data());
      return false;
    }
    TString base(gSystem->BaseName(path.Data()));
    if (gSystem->AccessPathName(base.Data()) == 0) { 
      // File or link exists - remove it 
      if (gSystem->Unlink(base) != 0) { 
	Error("Helper::AuxFile", "Failed to remove old %s", base.Data());
	return false;
      }
    }
    if (copy) 
      TFile::Cp(path, base);
    else 
      gSystem->Exec(Form("ln -s %s .", path.Data()));
    return true;
  }
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
  // --- Data members ------------------------------------------------
  TUrl        fUrl;     // The URI
  OptionList  fOptions; 
  Int_t       fVerbose;
};




// ===================================================================
Helper* 
Helper::Create(const TUrl& url, Int_t verbose)
{
  if (!url.IsValid()) { 
    Warning("Helper::Create", "URL is invalid");
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
    cl = "GridHelper";
  }
  else if (prot.EqualTo("local")) { 
    // Create Lite helper 
    cl = "LocalHelper";
  }
  else if (prot.EqualTo("proof")) { 
    // Create a Proof helper 
    if (host.IsNull()) 
      cl = "LiteHelper";
    else if (host.BeginsWith("alice-caf")) { 
      // AAF 
      cl = opts.Contains("plugin") ? "AAFPluginHelper" : "AAFHelper";
    }
    else 
      cl = "ProofHelper";
  }
  else if (prot.EqualTo("lite")) { 
    // Create a Proof helper 
    cl = "LiteHelper";
  }
  else if (prot.EqualTo("help")) {
    // Special HELP protocol
    if (host.Contains("options")) {
      std::cout << "Possible URL types and options are:" << std::endl;
      ShowFullHelp("LocalHelper");
      ShowFullHelp("ProofHelper");
      ShowFullHelp("LiteHelper");
      ShowFullHelp("AAFHelper");
      ShowFullHelp("AAFPluginHelper");
      ShowFullHelp("GridHelper");
      return 0;
    }
    std::cout << "Possible URL types are:" << std::endl;
    ShowUrlHelp("LocalHelper");
    ShowUrlHelp("ProofHelper");
    ShowUrlHelp("LiteHelper");
    ShowUrlHelp("AAFHelper");
    ShowUrlHelp("AAFPluginHelper");
    ShowUrlHelp("GridHelper");
    return 0;
  }
  // --- Check if we got a scheme ------------------------------------
  if (cl.IsNull()) {
    Error("Helper::Create", "Unknown scheme: %s", prot.Data());
    return 0;
  }

  // --- Use interpreter to make our object --------------------------
  Helper* helper = CreateObject(cl, url, verbose);
  if (!helper) {
    Error("Helper::Create", "Failed to make object of class %s", cl.Data());
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
