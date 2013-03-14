/**
 * @file   PluginHelper.C
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
#include "Helper.C"
#ifndef __CINT__
# include "AvailableSoftware.C"
# include "ParUtilities.C"
# include "OutputUtilities.C"
# include <TUrl.h>
# include <TString.h>
# include <TEnv.h>
# include <TProof.h>
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
struct PluginHelper : public Helper
{
  /** 
   * Constructor 
   * 
   * @param url  Url 
   * @param verbose Verbosity level 
   */
  PluginHelper(const TUrl& url, Int_t verbose)
    : Helper(url, verbose), fHandler(0), fUsePars(false), 
      fExtraLibs(), fExtraPars(), fExtraSrcs()
  {
    fHandler = new AliAnalysisAlien();

    fOptions.Add("aliroot", "VERSION", "AliROOT version", "last");
    fOptions.Add("root",    "VERSION", "ROOT version", "last");
    fOptions.Add("par", "Use par files");
    fOptions.Add("mode", "default|rec|sim", "AliROOT mode", "default");
    fOptions.Add("storage", "URL", "Location for external storage", "");    
    fOptions.Add("plugin", "Use AliEn handler");

    fExtraLibs.SetOwner();
    fExtraPars.SetOwner();
    fExtraSrcs.SetOwner();
  }
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  PluginHelper(const PluginHelper& o) 
    : Helper(o), fHandler(o.fHandler), fUsePars(o.fUsePars), 
      fExtraLibs(), fExtraPars(), fExtraSrcs()
  {}
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this 
   */
  PluginHelper& operator=(const PluginHelper& o) 
  {
    if (&o == this) return *this;
    Helper::operator=(o);
    fHandler = o.fHandler;
    fUsePars = o.fUsePars;
    return *this;
  }
  /** 
   * Destructor 
   */
  virtual ~PluginHelper() {}
  /** 
   * Load a library/PAR/script 
   * 
   * @param name   Name 
   * @param slaves If true, also load on slaves
   * 
   * @return true on success 
   */
  virtual Bool_t LoadLibrary(const TString& name, 
			     Bool_t slaves=true)
  {
    if (!fUsePars) {
      TString fullName(MakeLibraryName(name));
      Int_t ret = gSystem->Load(fullName);
      if (ret < 0) return false;
      if (slaves) {
	fHandler->AddAdditionalLibrary(fullName);
	fExtraLibs.Add(new TObjString(fullName));
      }
    }
    else { 
      if (!ParUtilities::Find(name)) { 
	Error("PluginHelper::LoadLibrary", "Failed to find PAR file %s", 
	      name.Data());
	return false;
      }
      if (!ParUtilities::Build(name)) { 
	Error("PluginHelper::LoadLibrary", "Failed to build PAR file %s", 
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
    if (!Helper::LoadSource(name, copy)) return false;
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
      Error("PluginHelper::LoadAliROOT", "Local AliROOT not available");
      return false;
    }
    Bool_t tmp = fUsePars;
    fUsePars   = false;
    
    if (!LoadLibrary("STEERBase"))     return false;
    if (!LoadLibrary("ESD"))           return false;
    if (!LoadLibrary("AOD"))           return false;
    if (!LoadLibrary("ANALYSIS"))      return false;
    if (!LoadLibrary("OADB"))          return false;
    if (!LoadLibrary("ANALYSISalice")) return false;
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
    if (fOptions.Has("aliroot")) aliroot = fOptions.Get("aliroot");
    if (fOptions.Has("root"))    root    = fOptions.Get("root");

    AvailableSoftware::Check(aliroot, root);
    fOptions.Set("aliroot", aliroot);
    fOptions.Set("root", root);

    fUsePars = fOptions.Has("par");

    fHandler->SetROOTVersion(root);
    fHandler->SetAliROOTVersion(aliroot);
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
      Error("PluginHelper::PostSetup", "No analysis manager defined");
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
    Helper::Print(option);
    fHandler->Print(option);
  }
  AliAnalysisAlien* fHandler;
  Bool_t fUsePars;
  TList  fExtraLibs;
  TList  fExtraPars;
  TList  fExtraSrcs;
};
#endif
//
// EOF
//
