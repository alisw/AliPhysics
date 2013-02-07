/**
 * @file   AAFHelper.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 19:02:14 2012
 * 
 * @brief  AAF analysis helper
 * 
 * @ingroup pwglf_forward_trains_helper
 * 
 */
#ifndef AAFHELPER_C
#define AAFHELPER_C
#include "ProofHelper.C"
#ifndef __CINT__
# include "AvailableSoftware.C"
# include <TUrl.h>
# include <TString.h>
# include <TProof.h>
# include <AliAnalysisManager.h>
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
 * where &lt;host@&t; is a known AAF (e.g., <tt>alice-caf.cern.ch</tt>)
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
 *       <dd> Use par files </dd>
 *     </dl>
 *   </dd>
 * </dl>  
 *
 * Note, this helper does not use the AliAnalysisAlien plugin
 *
 * @ingroup pwglf_forward_trains_helper
 */
struct AAFHelper : public ProofHelper
{
  /** 
   * Constructor 
   * 
   * @param url  Url 
   * @param verbose Verbosity 
   */
  AAFHelper(const TUrl& url, Int_t verbose)
    : ProofHelper(url, verbose)
  {
    fOptions.Add("aliroot", "VERSION", "AliROOT version", "last");
    fOptions.Add("root",    "VERISON", "ROOT version", "last");
    fOptions.Add("nocache", "Disable tree cache");
  }
  virtual ~AAFHelper() {}
  /** 
   * Get the name of the AliROOT par file to use 
   * 
   * @return String 
   */
  virtual const char* AliROOTParName() const
  {
    return Form("VO_ALICE@AliRoot::%s", fOptions.Get("aliroot").Data());
  }
  virtual Bool_t CreateAliROOTPar()
  {
    return true;
  }
  /** 
   * Set-up done before task set-ups. Overload ProofHelper::PreSetup
   * to specify the ROOT version using TProofMgr::SetROOTVersion
   * 
   * @return true on success 
   */
  virtual Bool_t PreSetup() 
  {
    TString aliroot("last");
    TString root("last");
    if (fOptions.Has("aliroot")) aliroot = fOptions.Get("aliroot");
    if (fOptions.Has("root"))    root    = fOptions.Get("root");

    AvailableSoftware::Check(aliroot, root);
    fOptions.Set("aliroot", aliroot);
    fOptions.Set("root", root);

    fBasePars = false;

    // Set this before we try to access the cluster
    gEnv->SetValue("XSec.GSI.DelegProxy", "2");

    TProof::Mgr(fUrl.GetHost())
      ->SetROOTVersion(Form("VO_ALICE@ROOT::%s", root.Data()));

    if (!ProofHelper::PreSetup()) return false;

    if (fOptions.Has("nocache")) 
      gProof->SetParameter("PROOF_UseTreeCache", 0);
    return true;
  }
  /** 
   * @return URI help string 
   */
  virtual const Char_t* UrlHelp() const 
  {
    return "proof://<host>/<dataset>?[&<options>][#<treename>]";
  }
  /** 
   * @return short description
   */
  virtual const char* Desc() const { return "AAF"; }
};
#endif
//
// EOF
//
