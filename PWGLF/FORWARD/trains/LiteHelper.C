/**
 * @file   LiteHelper.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 18:59:59 2012
 * 
 * @brief  Proof-Lite analysis helper
 * 
 * @ingroup pwglf_forward_trains_helper
 * 
 */
#ifndef LITEHELPER_C
#define LITEHELPER_C
#include "ProofHelper.C"
#ifndef __CINT__
# include "ChainBuilder.C"
# include <TUrl.h>
# include <TString.h>
# include <TChain.h>
# include <AliAnalysisManager.h>
#else
class TChain;
class TUrl;
#endif

// ===================================================================
/**
 * Handler of analysis in Proof-Lite.  This is triggered by URIs of the 
 * form 
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
 *       <dt><tt>workers=N[x]</tt></dt>
 *       <dd>Set the number of workers to use.  If <tt>x</tt> is appended, 
 *         then it's maximum number of workers per slave</dd>
 *       <dt><tt>par[=all]</tt></dt>
 *       <dd>Use PAR files.  If the value <tt>all</tt> is given, then also 
 *         PAR files of STEERBase, ESD, AOD, ANALYSIS, OADB, ANALYSISalice 
 *         are used. </dd>
 *       <dt><tt>mode=[default,rec,sim,train,custom]</tt></dt>
 *       <dd>Set the AliROOT mode.  If not specified <tt>default</tt> 
 *         is assumed.  See also CreateAliROOTPar</dd>
 *      </dl>
 *   </dd>
 * </dl>
 *       
 * @ingroup pwglf_forward_trains_helper
 */
struct LiteHelper : public ProofHelper
{
  /** 
   * Constructor 
   * 
   * @param url     Url 
   * @param verbose Verbosity
   */
  LiteHelper(const TUrl& url, Int_t verbose)
    : ProofHelper(url, verbose), fChain(0)
  {
    fOptions.Add("recursive","Recursive scan");
    fOptions.Add("pattern",  "GLOB", "File name pattern", "*.root");
    fOptions.Remove("dsname");
    fOptions.Remove("storage");
  }
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  LiteHelper(const LiteHelper& o) 
    : ProofHelper(o), fChain(o.fChain)
  {}
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this 
   */
  LiteHelper& operator=(const LiteHelper& o) 
  {
    if (&o == this) return *this;
    ProofHelper::operator=(o);
    fChain = o.fChain;
    return *this;
  }
  /** 
   * Destructor 
   */
  virtual ~LiteHelper() {}
  /** 
   * Set-up done before task set-ups 
   * 
   * @return true on success 
   */
  virtual Bool_t PreSetup() 
  {
    fUrl.SetProtocol("lite");
    Bool_t ret = ProofHelper::PreSetup();
    return ret;
  }
  /** 
   * Set-up done after task set-ups
   * 
   * @return true on success
   */
  virtual Bool_t PostSetup()
  {
    // -- Check for local chain --------------------------------------
    TString  pattern   = (fOptions.Has("pattern") ?fOptions.Get("pattern") :"");
    TString  treeName  = fUrl.GetAnchor();
    Bool_t   recursive = fOptions.Has("recursive");
    Bool_t   mc        = fOptions.Has("mc");
    TString  src       = fUrl.GetFile();
    UShort_t type      = ChainBuilder::CheckSource(src);
    if (type == ChainBuilder::kInvalid) {
      Error("LiteHelper", "Cannot generate TChain from %s", src.Data());
      return false;
    }

    // --- Create the chain ------------------------------------------
    fChain = ChainBuilder::Create(type, src, treeName, pattern, mc, recursive);
    if (!fChain) { 
      Error("PreSetup", "No chain defined");
      return false;
    }

    return ProofHelper::PostSetup();
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
    gProof->SetLogLevel(TMath::Max(fVerbose-2,0), 
			/* TProofDebug::kPacketizer| */
			TProofDebug::kLoop|
			/* TProofDebug::kSelector|
			TProofDebug::kOutput|
			TProofDebug::kInput|
			TProofDebug::kGlobal|*/
			TProofDebug::kPackage);
    if (nEvents < 0) nEvents = fChain->GetEntries();
    Long64_t ret = mgr->StartAnalysis("proof", fChain, nEvents);
    
    if (fVerbose > 2) 
      TProof::Mgr(fUrl.GetUrl())->GetSessionLogs()->Save("*","lite.log");
    return ret;
  }

  /** 
   * @return URL help string
   */
  virtual const Char_t* UrlHelp() const 
  {
    return "lite://<datadir_or_list>[?<options>][#<treeName]";
  }
  /** 
   * @return The short description
   */
  virtual const char* Desc() const { return "PROOF-lite"; }
  /** Our chain */
  TChain* fChain;
};
#endif
//
// EOF
//

