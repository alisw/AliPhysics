/**
 * @file   LiteRailway.C
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
#include "ProofRailway.C"
#ifndef __CINT__
# include <TUrl.h>
# include <TString.h>
# include <TChain.h>
# include <TDSet.h>
# include <AliAnalysisManager.h>
# include <AliVEventHandler.h>
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
struct LiteRailway : public ProofRailway
{
  /** 
   * Constructor 
   * 
   * @param url     Url 
   * @param verbose Verbosity
   */
  LiteRailway(const TUrl& url, Int_t verbose)
    : ProofRailway(url, verbose), fChain(0)
  {
    fOptions.Add("recursive","Recursive scan");
    fOptions.Add("trackref", "For MC input, check TrackRef.root presence");
    fOptions.Add("scan",     "Scan for number of events in chain");
    fOptions.Add("clean",    "Clean chain elements");
    fOptions.Add("pattern",  "GLOB", "File name pattern", "*.root");
    fOptions.Remove("dsname");
    fOptions.Remove("storage");
  }
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  LiteRailway(const LiteRailway& o) 
    : ProofRailway(o), fChain(o.fChain)
  {}
  /** 
   * Assignment operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this 
   */
  LiteRailway& operator=(const LiteRailway& o) 
  {
    if (&o == this) return *this;
    ProofRailway::operator=(o);
    fChain = o.fChain;
    return *this;
  }
  /** 
   * Destructor 
   */
  virtual ~LiteRailway() {}
  /** 
   * Set-up done before task set-ups 
   * 
   * @return true on success 
   */
  virtual Bool_t PreSetup() 
  {
    fUrl.SetProtocol("lite");
    Bool_t ret = ProofRailway::PreSetup();
    return ret;
  }
  /** 
   * Set-up done after task set-ups
   * 
   * @return true on success
   */
  virtual Bool_t PostSetup()
  {
    // --- Create the chain ------------------------------------------
    fChain = LocalChain();
    if (!fChain) return false;

    return ProofRailway::PostSetup();
  }
  /** 
   * Load extra sources.  Since we're on a single host, we might as
   * well load it directly from the working directory rather than by
   * uploading to the slaves cache.
   * 
   */
  virtual Bool_t LoadExtraSrcs()
  {
    TString    tmp2 = fExtraSrcs.Strip(TString::kBoth, ':');
    TObjArray* srcs = tmp2.Tokenize(":");
    TIter      next2(srcs);
    TObject*   obj = 0;
    while ((obj = next2())) {
      TString full = gSystem->ConcatFileName(gSystem->WorkingDirectory(),
					     obj->GetName());
      Info("LoadExtraSrcs", "Will load %s", full.Data());
      Int_t ret = gProof->Load(Form("%s+g", full.Data()), true,true);
      if (ret < 0) { 
	Error("ProofRailway::PostSetup", "Failed to compile %s",obj->GetName());
	return false;
      }
    }
    return true;
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
    Long64_t off = fOptions.AsLong("offset", 0);
    if (nEvents > 0 && nEvents < off) {
      Warning("Run", "Number of events %lld < offset (%lld), stopping", 
	      nEvents, off);
      return 0;
    }
    Long64_t ret = mgr->StartAnalysis("proof", fChain, nEvents, off);
    
    if (fVerbose > 2) 
      TProof::Mgr(fUrl.GetUrl())->GetSessionLogs()->Save("*","lite.log");
    return ret;
  }
  /** 
   * Path of output 
   * 
   * @return Path to output - possibly a data set
   */
  virtual TString OutputPath() const 
  {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return "";

    AliVEventHandler* outH = mgr->GetOutputEventHandler();
    if (!outH) return "";
    
    TString ret = gSystem->ConcatFileName(gSystem->WorkingDirectory(),
					  outH->GetOutputFileName());
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

