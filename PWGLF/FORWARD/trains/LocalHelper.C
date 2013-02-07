/**
 * @file   LocalHelper.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 18:59:42 2012
 * 
 * @brief  Local analysis helper
 * 
 * @ingroup pwglf_forward_trains_helper
 * 
 */
#ifndef LOCALHELPER_C
#define LOCALHELPER_C
#ifndef __CINT__
# include "Helper.C"
# include "ChainBuilder.C"
# include <TUrl.h>
# include <TString.h>
# include <TSystem.h>
# include <AliAnalysisManager.h>
#else
class TChain;
class Helper;
class TUrl;
#endif

// ===================================================================
/** 
 * Handle local analysis jobs 
 *
 * This is triggered by URIs of the form 
 *
 * @code 
 * local:///<datadir>[?<options>][#treeName]
 * local:///<collection>[?<options>][#treeName]
 * local:///<file>[?<options>][#treeName]
 * @endcode 
 *
 * where 
 *
 * <dl>
 *   <dt><tt>&lt;datadir&gt;</tt></dt>
 *   <dd>is the base directory holding data files </dd>
 *   <dt><tt>&lt;collection&gt;</tt></dt>
 *   <dd>is an ASCII or XML list of input sources</dd>
 *   <dt><tt>&lt;file&gt;</tt></dt>
 *   <dd>is a single ROOT file</dd>
 *   <dt><tt>&lt;options&gt;</tt></dt>
 *   <dd>A &amp; separated list of options
 *     <dl>
 *       <dt><tt>recursive</tt></dt>
 *       <dd>Scan &lt;datadir&gt; recursively</dd>
 *       <dt><tt>mc</tt></dt>
 *       <dd>Scan also for MC files (<tt>galice.root</tt>, 
 *          <tt>Kinematics.root</tt>, and <tt>TrackRefs.root</tt>) when 
 *          scanning &lt;datadir&gt;</dd>
 *       <dt><tt>pattern=&lt;GLOB&gt;</tt></dt>
 *       <dd>Shell glob pattern that files must check when scanning 
 *         &lt;datadir&gt;</dd>
 *      </dl>
 *   </dd>
 * </dl>
 *       
 * @ingroup pwglf_forward_trains_helper
 */
struct LocalHelper : public Helper
{
  /** 
   * Constructor 
   * 
   * @param url   Url 
   * @param verbose Verbosity level 
   */
  LocalHelper(const TUrl& url, Int_t verbose)
    : Helper(url, verbose), fChain(0)
  {
    fOptions.Add("recursive","Scan recursive");
    fOptions.Add("pattern",  "GLOB", "File name pattern", "*.root");
  }
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  LocalHelper(const LocalHelper& o) 
    : Helper(o), fChain(o.fChain)
  {}
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this 
   */
  LocalHelper& operator=(const LocalHelper& o) 
  {
    if (&o == this) return *this;
    Helper::operator=(o);
    fChain = o.fChain;
    return *this;
  }
  /** 
   * Destructor 
   */
  virtual ~LocalHelper() {}
  /** 
   * Load a library 
   * 
   * @param name    Name of library 
   * 
   * @return true on success
   */
  virtual Bool_t LoadLibrary(const TString& name, Bool_t)
  {
    Int_t ret = gSystem->Load(MakeLibraryName(name));
    return ret >= 0;
  }
  /** 
   * Get the execution mode 
   *
   * @return Always kLocal
   */
  virtual UShort_t Mode() const { return kLocal; }
  /**
   * Get the mode string used for AliAnalysisManager::StartAnalysis
   */
  virtual const char* ModeString() const { return "local"; }
  /** 
   * Set-up done before task set-ups 
   * 
   * @return true on success 
   */
  virtual Bool_t PreSetup() 
  {
    return true;
  }
  /** 
   * Set-up done after the task set-ups 
   *
   * @return true on success 
   */
  virtual Bool_t PostSetup() 
  {
    TString treeName(fUrl.GetAnchor());
    TString pattern(fOptions.Has("pattern") ? fOptions.Get("pattern") : "");
    Bool_t  recursive = fOptions.Has("recursive");
    Bool_t  mc        = fOptions.Has("mc");

    fChain = ChainBuilder::Create(fUrl.GetFile(), treeName, 
				  pattern, mc, recursive);

    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) { 
      Error("PostSetup", "No analysis manager defined");
      return false;
    }
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
    
    if (nEvents < 0) nEvents = fChain->GetEntries();
    return mgr->StartAnalysis(fUrl.GetProtocol(), fChain, nEvents);
  }
  /** 
   * @return URL help string
   */
  virtual const Char_t* UrlHelp() const 
  {
    return "local://<datadir or list>[?<options>][#<treeName>]";
  }
  /** 
   * @return Short description 
   */
  virtual const char* Desc() const { return "Local"; }
  TChain* fChain; // Our chain
};
#endif
//
// EOF
//
