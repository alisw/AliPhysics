/**
 * @file   VAFRailway.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 19:02:14 2012
 * 
 * @brief  AAF analysis helper
 * 
 * @ingroup pwglf_forward_trains_helper
 * 
 */
#ifndef VAFHELPER_C
#define VAFHELPER_C
#include "ProofRailway.C"
#ifndef __CINT__
# include <TUrl.h>
# include <TString.h>
# include <TProof.h>
# include <TSystem.h>
#else
class TUrl;
#endif

// ===================================================================
/**
 * Handle analysis on Proof-on-Demand
 * 
 * This helper is triggered by a URL of the form 
 *
 * @code
 * proof://host/<datadir>[?<options>][#<treename>]
 * @endcode 
 * where &lt;host@&t; is a known AAF (e.g., <tt>alivafh</tt>)
 * <dl>
 *   <dt>&lt;datadir&gt;</dt>
 *   <dd>Base directory to search for data</dd>
 *   <dt>&lt;pattern&gt;</dt>
 *   <dd>Filename pattern to search for. '%' is a wild-card</dd>
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
 *       <dt><tt>pattern=</tt>&lt;file pattern&gt;</dt>
 *       <dd>The file name pattern to search for</dd>
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
struct VAFRailway : public ProofRailway
{
  /** 
   * Constructor 
   * 
   * @param url  Url 
   * @param verbose Verbosity 
   */
  VAFRailway(const TUrl& url, Int_t verbose)
    : ProofRailway(url, verbose)
  {
    fOptions.Add("nocache", "Disable tree cache");
    fOptions.Add("pattern","SEARCH","Search pattern", "");
    fOptions.Add("alien","Enable ALIEN file access",true); //MUST be true
  }
  virtual ~VAFRailway() {}
  /** 
   * Do some sanity checks, then load 
   * 
   * @return See ProofRailway::LoadAliROOT 
   */
  virtual Bool_t LoadAliROOT()
  {
    Bool_t enabPhys = true;
    TString aliVer  = gSystem->Getenv("VafAliRootVersion");
    TString aliPhys = gSystem->Getenv("VafAliPhysicsVersion");
    if (!aliVer.IsNull()) {
      Info("VAFRailway::LoadAliROOT", "Using AliROOT=%s", aliVer.Data());
      if (!aliPhys.IsNull()) {
	Warning("VAFRailway::LoadAliROOT",
		"AliPhysics not loaded, even though version %s was requested",
		aliPhys.Data());
      }
      enabPhys = false;
    }
    else if (aliPhys.IsNull()) {
      Error("VAFRailway::LoadAliROOT", "Neither AliROOT nor AliPhysics "
	    "versions specified, giving up");
      return false;
    }
    else
      Info("VAFRailway::LoadAliROOT", "Using AliPhysics=%s", aliPhys.Data());

    fOptions.Set("alien");
    return ProofRailway::LoadAliROOT();
  }        
  /** 
   * Get the name of the Aliphysics par file to use 
   * 
   * @return String 
   */
  virtual const char* AliPhysicsParName() const { return 0; }
  virtual Bool_t CreateAliPhysicsPar() { return true; }
  virtual Bool_t EnableAliPhysics() { return true; }
  /** 
   * Get the name of the AliROOT par file to use 
   * 
   * @return String 
   */
  virtual const char* AliROOTParName() const
  {
    return "/afs/cern.ch/alice/offline/vaf/AliceVaf.par";
  }
  virtual Bool_t CreateAliROOTPar()
  {
    return true;
  }
  /** 
   * Set-up done before task set-ups. Overload ProofRailway::PreSetup
   * to specify the ROOT version using TProofMgr::SetROOTVersion
   * 
   * @return true on success 
   */
  virtual Bool_t PreSetup() 
  {
    if (!gSystem->Getenv("VafPodRemoteEnv")) {
      Printf("\n"
	     "=== ERROR:\n"
	     "Are you sure you're on alivaf-XXX and did 'vaf-enter'? "
	     "'cause I'm not so I have to give up!\n");
      return false;
    }

    fBasePars = false;

    if (!ProofRailway::PreSetup()) return false;

    if (fOptions.Has("nocache")) 
      gProof->SetParameter("PROOF_UseTreeCache", 0);
    return true;
  }
    /** 
   * Get the data-set name. On VAF, this is a search description
   * 
   * @param dsName On return, must contain the data set name 
   */
  virtual void GetDataSet(TString& dsName)
  {
    dsName = "Find;";
    dsName.Append(Form("BasePath=%s;", fUrl.GetFile()));
    TString pat = fOptions.Get("pattern");
    Int_t   idx = pat.Index("@");
    if (idx != kNPOS) {
      TString fn = pat(0,idx);
      TString an = pat(idx+1,pat.Length()-idx-1);
      dsName.Append(Form("FileName=%s;",fn.Data()));
      dsName.Append(Form("Anchor=%s;", an.Data()));
    }
    else
      dsName.Append(Form("FileName=%s;",pat.Data()));
    dsName.Append(Form("Tree=/%s;", fUrl.GetAnchor()));
    dsName.Append("Mode=remote;");
      
  }
  /** 
   * Override base class to replace the proof protocol with pod. 
   * 
   * @param url Connection URL
   * @param opt Possible options 
   * 
   * @return true on success 
   */
  virtual Bool_t Connect(const TUrl& url, const TString& opt)
  {
    TUrl tmp(url);
    tmp.SetProtocol("pod");
    return ProofRailway::Connect(tmp, opt);
  }
  /** 
   * @return URI help string 
   */
  virtual const Char_t* UrlHelp() const 
  {
    return "proof://alivaf/<datadir>?[&<options>][#<treename>]";
  }
  /** 
   * @return short description
   */
  virtual const char* Desc() const { return "VAF"; }
};
#endif
//
// EOF
//
