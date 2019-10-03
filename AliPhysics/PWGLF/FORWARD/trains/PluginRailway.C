/**
 * @file   PluginRailway.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 18:57:18 2012
 * 
 * @brief  Base class for helpers using the AliAnalysisAlien plugin
 *
 * @ingroup pwglf_forward_trains_helper
 * 
 */
#ifndef PLUGINHELPER_C
#define PLUGINHELPER_C
#include "Railway.C"
#ifndef __CINT__
# include "AvailableSoftware.C"
# include "ParUtilities.C"
# include "OutputUtilities.C"
# include <TUrl.h>
# include <TString.h>
# include <TEnv.h>
# include <AliAnalysisManager.h>
# include <AliAnalysisAlien.h>
#else
class TUrl;
class AliAnalysisAlien;
#endif

// ===================================================================
/**
 * Handle analysis on using the AliAnalysisAlien plugin - i.e., AAF or Grid 
 * 
 * This helper is triggered by a URL of the form 
 *
 * @code
 * <protocol>://[<user>@][<host>][:<port>]/[<file>][?<options>][#<anchor>]
 * @endcode 
 * where &lt;options&gt; contains <tt>plugin</tt>
 * <dl>
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
 *       <dd>Use par files</dd>
 *       <dt><tt>aliroot=&lt;version&gt;</tt></dt>
 *       <dd>Set AliROOT version to use </dd>
 *       <dt><tt>root=&lt;version&gt;</tt></dt>
 *       <dd>Set ROOT version to use </dd>
 *     </dl>
 *   </dd>
 * </dl>  
 *
 * @ingroup pwglf_forward_trains_helper
 */
struct PluginRailway : public Railway
{
  /** 
   * Constructor 
   * 
   * @param url  Url 
   * @param verbose Verbosity level 
   */
  PluginRailway(const TUrl& url, Int_t verbose)
    : Railway(url, verbose), fHandler(0), fUsePars(false), 
      fTestBuild(true), fExtraLibs(), fExtraPars(), fExtraSrcs()
  {
    fHandler = new AliAnalysisAlien();

    fOptions.Add("aliphysics", "VERSION", "AliPhysics version", "last");
    fOptions.Add("aliroot", "VERSION", "AliROOT version", "last");
    fOptions.Add("root",    "VERSION", "ROOT version", "last");
    fOptions.Add("par", "Use par files");
    fOptions.Add("mode", "default|rec|sim", "AliROOT mode", "default");
    fOptions.Add("storage", "URL", "Location for external storage", "");    
    fOptions.Add("plugin", "Use AliEn handler");
    fOptions.Add("testpar", "Test build PARs");

    fExtraLibs.SetOwner();
    fExtraPars.SetOwner();
    fExtraSrcs.SetOwner();
  }
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  PluginRailway(const PluginRailway& o) 
    : Railway(o), fHandler(o.fHandler), fUsePars(o.fUsePars), 
      fTestBuild(o.fTestBuild), fExtraLibs(), fExtraPars(), fExtraSrcs()
  {}
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this 
   */
  PluginRailway& operator=(const PluginRailway& o) 
  {
    if (&o == this) return *this;
    Railway::operator=(o);
    fHandler   = o.fHandler;
    fUsePars   = o.fUsePars;
    fTestBuild = o.fTestBuild;
    return *this;
  }
  /** 
   * Destructor 
   */
  virtual ~PluginRailway() {}
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
  void UsePar(Bool_t& use)
  {
    Bool_t tmp = fUsePars;
    fUsePars   = use;
    use        = tmp;
  }
  /** 
   * Add an include path
   * 
   * @param path The path to add for header search 
   * 
   * @return true on success
   */
  virtual Bool_t AddIncludePath(const TString& path)
  {
    Railway::AddIncludePath(path); // Also do on client 
    TString p(path);
    if (!p.BeginsWith("-")) {
      // IF the path does not begin with a -, we assume its a path
      p.Prepend("-I");
    }
    fHandler->AddIncludePath(p);
    return true;
  }
  /** 
   * Load a library/PAR/script 
   * 
   * @param name   Name 
   * @param slaves If true, also load on slaves
   * @param forcePar if true, force load as PAR 
   * 
   * @return true on success 
   */
  virtual Bool_t LoadLibrary(const TString& name, 
			     Bool_t slaves=true,
			     Bool_t forcePar=false)
  {
    Info("LoadLibrary", "Loading %s (slaves=%d, forcePar=%d)",
	 name.Data(), slaves, forcePar);
    if (!fUsePars && !forcePar) {
      TString fullName(MakeLibraryName(name));
      Int_t ret = gSystem->Load(fullName);
      if (ret < 0) return false;
      if (slaves) {	
	fHandler->AddAdditionalLibrary(fullName);
	fExtraLibs.Add(new TObjString(fullName));
      }
    }
    else {
      if (fExtraPars.FindObject(name)) return true;
      Info("LoadLibrary", "Looking for PAR %s", name.Data());
      if (!ParUtilities::Find(name)) { 
	Error("PluginRailway::LoadLibrary", "Failed to find PAR file %s", 
	      name.Data());
	return false;
      }
      Info("LoadLibrary", "Building PAR %s", name.Data());
      if (fTestBuild && !ParUtilities::Build(name)) { 
	Error("PluginRailway::LoadLibrary", "Failed to build PAR file %s", 
	      name.Data());
	return false;
      }
      fHandler->EnablePackage(name);
      fExtraPars.Add(new TObjString(name));
    }
    return true;
  }
  /** 
   * Load a source file, and compile it 
   * 
   * @param name Name of the source file 
   * @param copy If true, copy here instead of link
   * 
   * @return true on success
   */
  virtual Bool_t LoadSource(const TString& name, bool copy=false)
  {
    static TString s;
    if (!Railway::LoadSource(name, copy)) return false;
    s.Append(Form(" %s", gSystem->BaseName(name.Data())));
    fHandler->SetAnalysisSource(s);
    fExtraSrcs.Add(new TObjString(name));
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
      Error("PluginRailway::LoadAliROOT", "Local AliROOT not available");
      return false;
    }
    Bool_t tmp = fUsePars;
    fUsePars   = false;
    
    if (!LoadLibrary("STEERBase"))     return false;
    if (!LoadLibrary("ESD"))           return false;
    if (!LoadLibrary("AOD"))           return false;
    if (!LoadLibrary("ANALYSIS"))      return false;
    if (!LoadLibrary("ANALYSISalice")) return false;
    fUsePars = tmp;

    return true;
  }
  /** 
   * Set-up to load the AliPHYSICS libraries 
   * 
   * @return true on success
   */
  virtual Bool_t LoadAliPhysics()
  {
    if (!gSystem->Getenv("ALICE_PHYSICS")) { 
      Error("PluginRailway::LoadAliPhysics", "Local AliPhysics not available");
      return false;
    }
    Bool_t tmp = fUsePars;
    fUsePars   = false;
    
    // if (!LoadLibrary("OADB"))          return false;
    fUsePars = tmp;

    return true;
  }
  /** 
   * Set-up done before task set-ups 
   * 
   * @return true on success 
   */
  virtual Bool_t PreSetup() 
  {
    // --- Set prefered GSI method ---------------------------------
    gEnv->SetValue("XSec.GSI.DelegProxy", "2");

    TString aliroot("last");
    TString root("last");
    TString aliphysics("last");
    if (fOptions.Has("aliphysics")) aliphysics = fOptions.Get("aliphysics");
    if (fOptions.Has("aliroot"))    aliroot    = fOptions.Get("aliroot");
    if (fOptions.Has("root"))       root       = fOptions.Get("root");

    AvailableSoftware::Check(aliphysics, aliroot, root);
    fOptions.Set("aliphysics", aliphysics);
    fOptions.Set("aliroot", aliroot);
    fOptions.Set("root", root);

    fUsePars   = fOptions.Has("par");
    fTestBuild = fOptions.Has("testpar");

    // fHandler->SetAPIVersion("V1.1x");
    fHandler->SetROOTVersion(root);
    fHandler->SetAliROOTVersion(aliroot);
    fHandler->SetAliPhysicsVersion(aliphysics);
    // Add AliPhysics include path
    fHandler->AddIncludePath("-I.");
    fHandler->AddIncludePath("-I$ALICE_PHYSICS/include");
    fHandler->AddIncludePath("-I$ALICE_ROOT/include");
    // Execute through interpreter until patch is applied
    fHandler->SetDropToShell(false);
    if (fOptions.Has("mode"))
      fHandler->SetAliRootMode(fOptions.Get("mode"));
    else 
      fHandler->SetAliRootMode("default");
    
    return true;
  }
  /** 
   * Set-up done after the task set-ups 
   *
   * @return true on success 
   */
  virtual Bool_t PostSetup() 
  {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) { 
      Error("PluginRailway::PostSetup", "No analysis manager defined");
      return false;
    }
    mgr->SetGridHandler(fHandler);

    fHandler->SetJobTag(mgr->GetName());
    fHandler->SetAnalysisMacro(Form("%s.C", mgr->GetName()));

    if (fOptions.Has("storage"))
      OutputUtilities::RegisterStorage(fOptions.Get("storage"));

    return true;
  };
  /** 
   * Pure virtual overload
   * 
   * @return 
   */
  virtual Long64_t Run(Long64_t) = 0;
  /** 
   * Pure virtual overload 
   * 
   * @return Short description
   */
  virtual const char* Desc() const = 0;
  /** 
   * Pure virtual overload 
   * 
   * @return URI help string 
   */
  virtual const Char_t* UrlHelp() const = 0;
  /** 
   * Overload 
   * 
   * @param option Options 
   */
  virtual void Print(Option_t* option="") const 
  {
    Railway::Print(option);
    fHandler->Print(option);
  }
  AliAnalysisAlien* fHandler;
  Bool_t fUsePars;
  Bool_t fTestBuild;
  TList  fExtraLibs;
  TList  fExtraPars;
  TList  fExtraSrcs;
};
#endif
//
// EOF
//
