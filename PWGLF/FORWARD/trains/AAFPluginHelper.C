/**
 * @file   AAFPluginHelper.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 19:01:45 2012
 * 
 * @brief  AAF (using AliAnalysisAlien) analysis helper
 * 
 * @ingroup pwglf_forward_trains_helper
 * 
 */
#ifndef AAFPLUGINHELPER_C
#define AAFPLUGINHELPER_C
#include "PluginHelper.C"
#ifndef __CINT__
# include <TUrl.h>
# include <TString.h>
# include <AliAnalysisManager.h>
# include <AliAnalysisAlien.h>
#else
class TUrl;
class AliAnalysisAlien;
#endif

// ===================================================================
/**
 * Handle analysis on an Alice Analysis Facility (AAF)
 * 
 * This helper is triggered by a URL of the form 
 *
 * @code
 * proof://[<user>@]<host>[:<port>]/<dsname>[?<options>][#<treename>]
 * @endcode 
 * where &lt;host&gt; is a known AAF (e.g., <tt>alice-caf.cern.ch</tt>),
 * and the &lt;options&gt; contains <tt>plugin</tt>
 * <dl>
 *   <dt>&lt;user&gt;</dt>
 *   <dd>Optional user name</dd>
 *   <dt>&lt;host&gt;</dt>
 *   <dd>PROOF cluster master host</dd>
 *   <dt>&lt;port&gt;</dt>
 *   <dd>Optional PROOF cluster port on master host</dd>
 *   <dt>&lt;dsname&gt;</dt>
 *   <dd>Data set name</dd>
 *   <dt>&lt;treename&gt;</dt>
 *   <dd>Optional tree name in data set, often <tt>esdTree</tt> or
 *   <tt>aodTree</tt></dd>
 *   <dt>&lt;options&gt;</dt>
 *   <dd>List of options separated by an &amp;
 *     <dl>
 *       <dt><tt>dsname</tt>[=&lt;output dataset&gt;]</dt>
 *       <dd>Register tree output (e.g., AOD) as a new data set on the
 *         PROOF cluster. If &lt;output dataset&gt; is not specified, take
 *         the name of the train.</dd>
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
 *       <dt><tt>workers=</tt><i>N</i><tt>[x]</tt></dt>
 *       <dd>Set the number of workers to use.  If <tt>x</tt> is appended, 
 *         then it's maximum number of workers per slave</dd>
 *     </dl>
 *   </dd>
 * </dl>  
 *
 * @ingroup pwglf_forward_trains_helper
 */
struct AAFPluginHelper : public PluginHelper
{
  /** 
   * Constructor 
   * 
   * @param url  Url 
   * @param verbose Verbosity level
   */
  AAFPluginHelper(const TUrl& url, Int_t verbose)
    : PluginHelper(url, verbose)
  {
    fOptions.Add("workers", "N[x]", "Number of workers to use", 0);
    fOptions.Add("dsname",  "NAME", "Make output dataset", "");
    fOptions.Add("wrapper", "CMD", "Wrapper command", "");
    fOptions.Add("clear",   "Clear all packages");
    fOptions.Add("reset",   "soft|hard", "Reset cluster", "hard");

  }
  /** 
   * Destructor
   */
  virtual ~AAFPluginHelper() {}
  /** 
   * Called before setting up 
   * 
   * @return true on success 
   */
  virtual Bool_t PreSetup()
  {
    // --- Handle software options -----------------------------------
    TString root = fOptions.Get("root");
    fHandler->SetRootVersionForProof(Form("VO_ALICE@ROOT::%s", root.Data()));
    fHandler->SetProofCluster(fUrl.GetHost());
    fHandler->SetProofDataSet(fUrl.GetFile());

    // --- Handle worker options -------------------------------------
    if (fOptions.Has("workers")) {
      TString nwork = fOptions.Get("workers");
      if (nwork.EndsWith("x")) 
	fHandler->SetNproofWorkersPerSlave(nwork.Atoi());
      else 
	fHandler->SetNproofWorkers(nwork.Atoi());
    }
    
    // --- Check if we're using a wrapper ----------------------------
    if (fOptions.Has("wrapper")) { 
      TString wrapper = fOptions.Get("wrapper");
      if (wrapper.IsNull()) 
	// In case of no argument, use GDB 
	// Just run and backtrace 
	wrapper = "/usr/bin/gdb --batch -ex run -ex bt --args";
      Info("ProofHelper::PreSetup", "Using wrapper command: %s", 
	   wrapper.Data());
      TProof::AddEnvVar("PROOF_WRAPPERCMD", wrapper);
    }
    
    // --- Check if we need to clear packages ------------------------
    fHandler->SetClearPackages(fOptions.Has("clear"));

    // --- Check if we need to reset first ---------------------------
    if (fOptions.Has("reset")) { 
      TString reset = fOptions.Get("reset");
      Bool_t  hard  = (reset.IsNull() || 
		       reset.EqualTo("hard", TString::kIgnoreCase));
      Info("AAFPluginHelper::PreSetup", "Will do a %s reset of %s", 
	   hard ? "hard" : "soft", fUrl.GetHost());
      fHandler->SetProofReset(hard ? 2 : 1);
    }
    
    return PluginHelper::PreSetup();
  }
  /** 
   * Set-up done after the task set-ups 
   *
   * @return true on success 
   */
  virtual Bool_t PostSetup() 
  {
    if (!PluginHelper::PostSetup()) return false;
    if (fOptions.Has("dsname")) 
      OutputUtilities::RegisterDataset(fOptions.Get("dsname"));

    return true;
  };
  /** 
   * Get the mode identifier 
   * 
   * @return Always kProof
   */
  virtual UShort_t Mode() const { return kProof; }
  /**
   * Get the mode string used for AliAnalysisManager::StartAnalysis
   */
  virtual const char* ModeString() const { return "proof"; }
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
    
    TString dsName(fUrl.GetFile());
    if (fUrl.GetAnchor() && fUrl.GetAnchor()[0] != '\0') 
      dsName.Append(Form("#%s", fUrl.GetAnchor()));
    return mgr->StartAnalysis(fUrl.GetProtocol(), dsName, nEvents);
  }
  /** 
   * @return URI help string 
   */
  virtual const Char_t* UrlHelp() const 
  {
    return "proof://<host>/<dataset>?plugin[&<options>][#<treename>]";
  }
  /**
   * @return Short description
   */
  virtual const char* Desc() const { return "CAF w/plugin"; }
};
#endif
//
// EOF
//
