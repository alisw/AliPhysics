/**
 * @file   ProofHelper.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 18:58:37 2012
 * 
 * @brief  
 * 
 *
 * @ingroup pwglf_forward_trains_helper
 * 
 */
#ifndef PROOFHELPER_C
#define PROOFHELPER_C
#include "Helper.C"
#ifndef __CINT__
# include "OutputUtilities.C"
# include "ParUtilities.C"
# include "ChainBuilder.C"
# include <TUrl.h>
# include <TString.h>
# include <TProof.h>
# include <TProofLog.h>
# include <TProofDebug.h>
# include <AliAnalysisManager.h>
# include <TEnv.h>
# include <TChain.h>
// For SendFile
# include <TSystem.h>
# include <TSlave.h>
# include <TSocket.h>
# include <cerrno>
#else
class TUrl;
class TChain;
#endif

// ===================================================================
/**
 * Handle analysis on a Proof farm. 
 * 
 * This helper is triggered by URIs of the form 
 *
 * @code
 * proof://[<user>@]<host>[:<port>]/<dsname>[?<options>][#<treename>]
 * @endcode 
 * where 
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
 *       <dt><tt>workers=N[x]</tt></dt>
 *       <dd>Set the number of workers to use.  If <tt>x</tt> is appended, 
 *         then it's maximum number of workers per slave</dd>
 *       <dt><tt>dsname</tt>[=&lt;output dataset&gt;]</dt>
 *       <dd>Register tree output (e.g., AOD) as a new data set on the
 *         PROOF cluster. If &lt;output dataset&gt; is not specified, take
 *         the name of the train.</dd>
 *       <dt><tt>par[=all]</tt></dt>
 *       <dd>Use PAR files.  If the value <tt>all</tt> is given, then also 
 *         PAR files of STEERBase, ESD, AOD, ANALYSIS, OADB, ANALYSISalice 
 *         are used. </dd>
 *       <dt><tt>mode=[default,rec,sim,train,custom]</tt></dt>
 *       <dd>Set the AliROOT mode.  If not specified <tt>default</tt> 
 *         is assumed.  See also CreateAliROOTPar</dd>
 *       <dt><tt>storage=&lt;url&gt;</tt></dt>
 *       <dd>Specify a non-default storage location for special output
 *         (e.g., AOD trees).  &lt;url&gt; should be a valid XRootd 
 *         server URI accessible to the slaves - e.g., 
 *         <tt>root://lxplus.cern.ch:10930//tmp</tt>.</dd>
 *     </dl>
 *   </dd>
 * </dl>  
 *
 * @ingroup pwglf_forward_trains_helper
 */
struct ProofHelper : public Helper
{
  /** 
   * Constructor 
   * 
   * @param url  Url 
   * @param verbose Verbosity level
   */
  ProofHelper(const TUrl& url, Int_t verbose)
    : Helper(url, verbose), 
      fExtraLibs(""),
      fExtraPars(""),
      fExtraSrcs(""),
      fUsePars(false), 
      fBasePars(false),
      fAuxFiles()
  {
    fOptions.Add("workers",  "N[x]", "Number of workers to use", 0);
    fOptions.Add("dsname",   "NAME", "Make output dataset", "");
    fOptions.Add("par",      "tasks|all", "Use par files",           "tasks");
    fOptions.Add("mode",     "default|rec|sim", "AliROOT mode",      "default");
    fOptions.Add("storage",  "URL", "Location for external storage", "");    
    fOptions.Add("wrapper",  "CMD", "Wrapper command", "");
    fOptions.Add("clear",    "PKGS", "Clear packages ','-separated", "");
    fOptions.Add("reset",    "soft|hard", "Reset cluster", "hard");
    fOptions.Add("feedback", "Enable feedback mechanism");
    fOptions.Add("env",      "SCRIPT", "Script to set-up environment","-none-");
    fOptions.Add("offset",   "EVENTS", "Skip this number of events", 0);
    if (!fUrl.GetUser() || fUrl.GetUser()[0] == '\0') 
      fUrl.SetUser(gSystem->GetUserInfo()->fUser);
    fAuxFiles.SetOwner();
  }
  ProofHelper(const ProofHelper& o) 
    : Helper(o),
      fExtraLibs(""),
      fExtraPars(""),
      fExtraSrcs(""),
      fUsePars(false), 
      fBasePars(false),
      fAuxFiles()
  {}
  ProofHelper& operator=(const ProofHelper& o) 
  {
    if (&o == this) return *this;
    Helper::operator=(o);
    fExtraLibs = o.fExtraLibs;
    fExtraPars = o.fExtraPars;
    fExtraSrcs = o.fExtraSrcs;
    fUsePars   = o.fUsePars;
    fBasePars  = o.fBasePars;
    // fAuxFiles;
    return *this;
  }
  /** 
   * Destructor 
   */
  virtual ~ProofHelper() {}
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
    Bool_t isBase = false;
    if (!fBasePars) { 
      if (name.EqualTo("STEERBase")      || 
	  name.EqualTo("ESD")            || 
	  name.EqualTo("AOD")            || 
	  name.EqualTo("ANALYSIS")       || 
	  name.EqualTo("OADB")           || 
	  name.EqualTo("ANALYSISalice")) 
	isBase = true;
    }
    if (!fUsePars || isBase) {
      Int_t ret = gSystem->Load(MakeLibraryName(name));
      if (ret < 0) return false;
      if (slaves) fExtraLibs.Append(Form(":%s", name.Data()));
    }
    else { 
      if (!ParUtilities::Find(name)) { 
	Error("ProofHelper::LoadLibrary", "Failed to find PAR file %s", 
	      name.Data());
	return false;
      }
      if (!ParUtilities::Build(name)) { 
	Error("ProofHelper::LoadLibrary", "Failed to build PAR file %s", 
	      name.Data());
	return false;
      }
      if (gProof->UploadPackage(name.Data(), TProof::kRemoveOld) < 0) {
	Error("ProofHelper::LoadLibrary", "Failed to upload PAR file %s", 
	      name.Data());
	return false;
      }
      fExtraPars.Append(Form(":%s", name.Data()));
    }
    return true;
  }
  /** 
   * Load a source file, and compile it 
   * 
   * @param name Name of the source file 
   * @param copy If true, copy not link 
   * 
   * @return true on success
   */
  virtual Bool_t LoadSource(const TString& name, bool copy=false)
  {
    if (!Helper::LoadSource(name, copy)) return false;
    fExtraSrcs.Append(Form(":%s", gSystem->BaseName(name.Data())));
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
      Error("ProofHelper::LoadAliROOT", "Local AliROOT not available");
      return false;
    }

    Bool_t tmp = fUsePars;
    fUsePars   = fBasePars;
    if (!LoadLibrary("STEERBase"))     return false;
    if (!LoadLibrary("ESD"))           return false;
    if (!LoadLibrary("AOD"))           return false;
    if (!LoadLibrary("ANALYSIS"))      return false;
    if (!LoadLibrary("OADB"))          return false;
    if (!LoadLibrary("ANALYSISalice")) return false;
    fUsePars = tmp;

    return CreateAliROOTPar();
  }
  /** 
   * Get the name of the AliROOT par file to use 
   * 
   * @return String 
   */
  virtual const char* AliROOTParName() const
  {
    return "ALIROOT";
  }
  /** 
   * Create an AliROOT par file from the executing AliROOT.  This PAR
   * file basically uses the environment of the client - that is, we
   * assume that the used AliROOT is accessible on the slaves - e.g.,
   * via an NFS export.
   * 
   * Note, the SETUP.C script take one argument - a TList of TNamed
   * parameters.  Parameters processed are	
   *
   * - ALIROOT_MODE=[default,aliroot,rec,sim,train]
   *   - default: Load base analysis libraries 
   *   - aliroot: Load $ALICE_ROOT/macros/loadlibs.C
   *   - rec:     Load $ALICE_ROOT/macros/loadlibsrec.C
   *   - sim:     Load $ALICE_ROOT/macros/loadlibssim.C
   * - ALIROOT_EXTRA_LIBS Colon separated list of additional (Ali)ROOT
   *   libraries to load on the slaves.
   * 
   * The generated PAR file is uploaded but not enabled until we have 
   * populated fExtraLibs.  The enabling takes place at the end of the 
   * set-up. 
   * 
   * @return true on success, false otherwise.     */
  virtual Bool_t CreateAliROOTPar()
  {
    if (fBasePars) return true;

    TString parName(AliROOTParName());
    TString parFile(Form("%s.par", parName.Data()));

    // --- Check if we have the drirectory already -------------------
    if (gSystem->AccessPathName(parName.Data()) == 0) { 
      // Let's remove it to get a clean slate 
      if (gSystem->Exec(Form("rm -rf %s", parName.Data())) != 0) {
	Error("ProofHelper", "Failed to remove %s", parName.Data());
	return false;
      }
    }
    // --- Check if the PAR file is there, and remove it if so -------
    if (gSystem->AccessPathName(parFile.Data()) == 0) { 
      if (gSystem->Unlink(parFile.Data()) != 0) { 
	Error("ProofHelper::CreateAliROOTPar", "Failed to remove %s", 
	      parFile.Data());
	return false;
      }
    }
      

    // Set-up directories 
    if (gSystem->MakeDirectory(parName) < 0) {
      Error("ProofHelper::CreateAliROOTPar", "Could not make directory '%s'", 
	    parName.Data());
      return false;
    }
    
    if (gSystem->MakeDirectory(Form("%s/PROOF-INF", parName.Data()))) {
      Error("ProofHelper::CreateAliROOTPar", 
	    "Could not make directory %s/PROOF-INF", 
	    parName.Data());
      return false;
    }

    std::ofstream b(Form("%s/PROOF-INF/BUILD.sh",parName.Data()));
    if (!b) { 
      Error("ProofHelper::CreateAliROOTPar", 
	    "Failed to make BUILD.sh shell script");
      return false;
    }
    b << "#!/bin/sh\n\n"
      << "# echo Nothing to do\n"
      << "exit 0\n"
      << std::endl;
    b.close();
    gSystem->Exec(Form("chmod a+x %s/PROOF-INF/BUILD.sh",parName.Data()));

    TString envScript = fOptions.Get("env");
    if (envScript.EqualTo("-none-", TString::kIgnoreCase)) 
      envScript = "";
    if (!envScript.IsNull()) { 
      // If an environment script was specified, copy that to the par
      if (gSystem->AccessPathName(envScript.Data()) == 0) { 
	// Copy script 
	if (gSystem->Exec(Form("cp %s %s/PROOF-INF/", envScript.Data(), 
			       parName.Data())) != 0) {
	  Error("ProofHelper", "Failed to copy %s", envScript.Data());
	  return false;
	}
      }
      else {
	Warning("CreateALIROOTPar", "Couldn't read %s", envScript.Data());
	envScript = "";
      }
    }
    std::ofstream s(Form("%s/PROOF-INF/SETUP.C", parName.Data()));
    if (!s) { 
      Error("ProofHelper::CreateAliROOTPar", 
	    "Failed to make SETUP.C ROOT script");
      return false;
    }
    s << "void SETUP(TList* opts) {\n";
    if (envScript.IsNull()) {
      s << "  gSystem->Setenv(\"ALICE\",\"" 
      << gSystem->Getenv("ALICE") << "\");\n"
      << "  gSystem->Setenv(\"ALICE_ROOT\",\"" 
      << gSystem->Getenv("ALICE_ROOT") << "\");\n"
      << "  gSystem->Setenv(\"ALICE_TARGET\",\"" 
	<< gSystem->Getenv("ALICE_TARGET") << "\");\n";
    if (gSystem->Getenv("OADB_PATH")) 
      s << "  gSystem->Setenv(\"OADB_PATH\",\"" 
	<< gSystem->Getenv("OADB_PATH") << "\");\n";
    }
    else { 
      s << "  gROOT->Macro(\"PROOF-INF/" << gSystem->BaseName(envScript.Data())
	<< "\");\n";
    }
    s  	<< "  gSystem->AddDynamicPath("
	<< "\"$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)\");\n";
    s << "  \n"
      << "  // Info(\"SETUP\",\"Loading ROOT libraries\");\n"
      << "  gSystem->Load(\"libTree\");\n"
      << "  gSystem->Load(\"libGeom\");\n"
      << "  gSystem->Load(\"libVMC\");\n"
      << "  gSystem->Load(\"libPhysics\");\n"
      << "  gSystem->Load(\"libMinuit\");\n"
      << "  \n";
    s << "  // Info(\"SETUP\",\"Parameter list:\");\n"
      << "  if (!opts) return;\n"
      << "  //opts->ls();\n"
      << "  \n";
    s << "  TObject* par = opts->FindObject(\"ALIROOT_MODE\");\n"
      << "  if (par) {\n"
      << "    // Info(\"SETUP\",\"ALIROOT mode: %s\", par->GetTitle());\n"
      << "    TString mode(par->GetTitle());\n"
      << "    if (mode.EqualTo(\"default\",TString::kIgnoreCase)) {\n"
      << "      gSystem->Load(\"libSTEERBase\");\n"
      << "      gSystem->Load(\"libESD\");\n"
      << "      gSystem->Load(\"libAOD\");\n"
      << "      gSystem->Load(\"libANALYSIS\");\n"
      << "      gSystem->Load(\"libOADB\");\n"
      << "      gSystem->Load(\"libANALYSISalice\");\n"
      << "    }\n"
      << "    else if (mode.EqualTo(\"aliroot\",TString::kIgnoreCase)) \n"
      << "      gROOT->Macro(\"$ALICE_ROOT/macros/loadlibs.C\");\n"
      << "    else if (mode.EqualTo(\"rec\",TString::kIgnoreCase)) \n"
      << "      gROOT->Macro(\"$ALICE_ROOT/macros/loadlibsrec.C\");\n"
      << "    else if (mode.EqualTo(\"sim\",TString::kIgnoreCase)) \n"
      << "      gROOT->Macro(\"$ALICE_ROOT/macros/loadlibssim.C\");\n"
      << "    else if (mode.EqualTo(\"train\",TString::kIgnoreCase)) \n"
      << "      gROOT->Macro(\"$ALICE_ROOT/macros/loadlibstrain.C\");\n"
      << "    else if (mode.EqualTo(\"custom\",TString::kIgnoreCase)) \n"
      << "      gROOT->Macro(\"$ALICE_ROOT/macros/loadlibstrain.C\");\n"
      << "  }\n"
      << "  \n";
    s << "  par = opts->FindObject(\"ALIROOT_EXTRA_LIBS\");\n"
      << "  if (par) {\n"
      << "    Info(\"SETUP\",\"Libaries to load: %s\\n\",par->GetTitle());\n"
      << "    TString tit(par->GetTitle());\n"
      << "    TObjArray* tokens = tit.Tokenize(\":\");\n"
      << "    TObject*   lib    = 0;\n"
      << "    TIter      next(tokens);\n"
      << "    while ((lib = next())) {\n"
      << "      TString libName(lib->GetName());\n"
      << "      if (!libName.BeginsWith(\"lib\")) libName.Prepend(\"lib\");\n"
      << "      // Info(\"SETUP\",\"Loading %s ...\",libName.Data());\n"
      << "      gSystem->Load(libName.Data());\n"
      << "    }\n"
      << "  }\n"
      << "}\n"
      << std::endl;
    s.close();

    Int_t ret = gSystem->Exec(Form("tar -czf %s %s",
				   parFile.Data(), parName.Data()));
    if (ret != 0) { 
      Error("ProofHelper::CreateAliROOTPar", "Failed to pack up PAR file %s",
	    parFile.Data());
      return false;
    }

    ret = gProof->UploadPackage(parFile.Data(),TProof::kRemoveOld);
    if (ret != 0) { 
      Error("ProofHelper::CreateAliROOTPar", 
	    "Failed to upload the AliROOT PAR file");
      return false;
    }
    // Note, the PAR isn't enabled until much later when we've
    // collected all the needed libraries in fExtraLibs
    return true;
  }
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
   * Set-up done before task set-ups 
   * 
   * @return true on success 
   */
  virtual Bool_t PreSetup()
  {
    // --- Set prefered GSI method ---------------------------------
    gEnv->SetValue("XSec.GSI.DelegProxy", "2");

    // --- Add ALICE_ROOT directory to search path for packages ----
    // Info("ProofHelper::PreSetup", "Set location of packages");
    gEnv->SetValue("Proof.GlobalPackageDirs", 
		   Form("%s:%s", 
			gEnv->GetValue("Proof.GlobalPackageDirs", "."), 
			gSystem->Getenv("ALICE_ROOT")));

    // --- Forming the URI we use to connect with --------------------
    TUrl connect(fUrl);
    connect.SetAnchor("");
    connect.SetFile("");
    connect.SetOptions("");

    // --- Check if we need to reset first ---------------------------
    if (fOptions.Has("reset")) { 
      TString reset = fOptions.Get("reset");
      Bool_t  hard  = (reset.IsNull() || 
		       reset.EqualTo("hard", TString::kIgnoreCase));
      Info("ProofHelper::PreSetup", "Doing a %s reset of %s", 
	   hard ? "hard" : "soft", connect.GetUrl());
      TProof::Reset(connect.GetUrl(), hard);
      Int_t secs = 3;
      Info("ProofHelper::PreSetup", 
	   "Waiting for %d second%s for things to settle", secs,
	   secs > 1 ? "s" : "");
      gSystem->Sleep(1000*secs);
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

    // --- PAR parameters --------------------------------------------
    fUsePars  = fOptions.Has("par");
    fBasePars = (fUsePars && 
		 fOptions.Get("par").EqualTo("all",TString::kIgnoreCase));

    // --- Connect to the cluster ------------------------------------
    TString opts;
    if (fOptions.Has("workers")) 
      opts.Append(Form("workers=%s", fOptions.Get("workers").Data()));
      
    Info("ProofHelper::PreSetup", "Connecting to %s with %soptions %s", 
	 connect.GetUrl(), 
	 opts.IsNull() ? "no " : "", 
	 opts.Data());
    TString proto(connect.GetProtocol());
    if (proto.BeginsWith("lite") && fOptions.Has("workers")) 
      TProof::Open(opts);
    else 
      TProof::Open(connect.GetUrl(), opts);
    // TProof::Open(connect.GetHost(), opts);
    if (!gProof) { 
      Error("ProofHelper::PreSetup", "Failed to open Proof connection %s", 
	    connect.GetUrl());
      return false;
    }
    
    // --- Check if we need to clear packages ------------------------
    if (fOptions.Has("clear")) {
      TString pkgs = fOptions.Get("clear");
      if (pkgs.IsNull() || pkgs.EqualTo("all", TString::kIgnoreCase)) { 
	// No value given, clear all 
	if (gProof->ClearPackages() != 0) 
	  Warning("ProofHelper::PreSetup", "Failed to lear all packages");
      }
      else { 
	// Tokenize on ',' and clear each package 
	TObjArray* pars = pkgs.Tokenize(",");
	TObject*   pkg  = 0;
	TIter      next(pars); 
	while ((pkg = next())) { 
	  if (gProof->ClearPackage(pkg->GetName()) != 0)
	    Warning("ProofHelper::PreSetup", "Failed to clear package %s", 
		    pkg->GetName());
	}
	pars->Delete();
      }
    }
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
      Error("ProofHelper::PostSetup", "No analysis manager defined");
      return false;
    }

    // --- Check for output ------------------------------------------
    if (fOptions.Has("dsname")) 
      OutputUtilities::RegisterDataset(fOptions.Get("dsname"));
    if (fOptions.Has("storage"))
      OutputUtilities::RegisterStorage(fOptions.Get("storage"));

    // --- Check for feedback mechanism ------------------------------
    if (!fOptions.Has("feedback"))
      gProof->ClearFeedback();

    // --- If we are not using PARs for Base, enable special PAR -----
    if (!fBasePars) {
      TString tmp(fExtraLibs.Strip(TString::kBoth,':'));
      TList* params = new TList;
      params->SetOwner(true);
      params->Add(new TNamed("ALIROOT_EXTRA_LIBS", tmp.Data()));
      if (fOptions.Has("mode"))
	params->Add(new TNamed("ALIROOT_MODE", fOptions.Get("mode").Data()));
      else
	params->Add(new TNamed("ALIROOT_MODE", "default"));
      Int_t ret = gProof->EnablePackage(AliROOTParName(), params, true);
      if (ret < 0) {
	Error("ProofHelper::EnableAliROOT", "Failed to enable AliROOT PAR %s", 
	      AliROOTParName());
	return false;
      }
    }

    // --- Make PAR file of Aux Files --------------------------------
    if (fAuxFiles.GetEntries() > 0) { 
      TString name = TString::Format("%s_auxfiles", mgr->GetName());
      ParUtilities::MakeAuxFilePAR(fAuxFiles, name);

      if (gProof->UploadPackage(name.Data(), TProof::kRemoveOld) < 0) 
	Error("ProofHelper::PostSetup", "Failed to upload PAR file %s", 
	      name.Data());
      else 
	fExtraPars.Append(Form(":%s", name.Data()));
    }
    
    // --- Load par files --------------------------------------------
    TString    tmp  = fExtraPars.Strip(TString::kBoth,':');
    TObjArray* pars = tmp.Tokenize(":");
    TObject*   obj  = 0;
    TIter      next(pars);
    while ((obj = next())) { 
      // Enable the package, but do not build on client - already done
      Int_t ret = gProof->EnablePackage(obj->GetName(), true);
      if (ret < 0) { 
	Error("ProofHelper::PostSetup", "Failed to enable PAR %s",
	      obj->GetName());
	return false;
      }
    }
    
    // --- Load extra sources ----------------------------------------
    TString    tmp2 = fExtraSrcs.Strip(TString::kBoth, ':');
    TObjArray* srcs = tmp2.Tokenize(":");
    TIter      next2(srcs);
    while ((obj = next())) { 
      Int_t ret = gProof->Load(Form("%s++g", obj->GetName()), true);
      if (ret < 0) { 
	Error("ProofHelper::PostSetup", "Failed to compile %s", obj->GetName());
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
    TString dsName(fUrl.GetFile());
    // if (fUrl.GetAnchor() && fUrl.GetAnchor()[0] != '\0') 
    //   dsName.Append(Form("#%s", fUrl.GetAnchor()));
    // Info("Run", "Output objects registered with PROOF:");
    // gProof->GetOutputList()->ls();
    Long64_t off = fOptions.AsLong("offset", 0);
    if (nEvents > 0 && nEvents < off) {
      Warning("Run", "Number of events %lld < offset (%lld), stopping", 
	      nEvents, off);
      return 0;
    }
    Long64_t ret = mgr->StartAnalysis(fUrl.GetProtocol(), dsName, nEvents, off);
    
    if (fVerbose > 10) 
      TProof::Mgr(fUrl.GetUrl())->GetSessionLogs()->Save("*","proof.log");
    return ret;
  }
#if 0
  Bool_t AddMonitor(const TString& path)
  {
    if (path.IsNull()) return true;

    TObjArray* tokens  = path.Tokenize("/");
    Int_t      nTokens = tokens->GetEntries();
    if (nTokens < 2) { 
      Error("AddMonitor", "Monitors must be of the form:\n"
	    "  <task>[:<slot>]/<name>\n"
	    "  <task>[:<slot>]/<path>/<name>");
      return false;
    }
    // --- Get the manager 
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();

    // --- Extract task and possibly slot number 
    TString& sTask  = static_cast<TObjString*>(tokens->At(0))->String();
    Int_t    slotNo = 0;
    Ssiz_t   colon  = sTask.Index(":");
    if (colon != kNPOS) { 
      TString sSlot = sTask(colon+1, sTask.Length()-colon-1);
      if (!sSlot.IsNull()) slotNo = sSlot.Atoi();
      sTask.Remove(colon, sTask.Length()-colon);
    }
    
    AliAnalysisTask* task = mgr->GetTask(sTask);
    if (!task) { 
      Error("AddMonitor", "Task \"%s\" not registered with manager", 
	    sTask.Data());
      return false;
    }
    AliAnalysisDataSlot* slot = task->GetOutputSlot(slotNo);
    if (!slot) { 
      Error("AddMonitor", "Task \"%s\" does not have an output slot at %d",
	    task->GetName(), slotNo);
      return false;
    }
    AliAnalysisDataContainer* cont = slot->GetContainer();
    if (!cont) {
      Error("AddMonitor", "Output slot %d of task \"%s\" has no container",
	    slotNo, task->GetName());
      return false;
    }
    Int_t    idx   = 1;
    TString& first = static_cast<TObjString*>(tokens->At(idx))->String(); 
    if (first.EqualTo(cont->GetName())) {
      idx++;
    }
    TObject* data = cont->GetData();
    TObject* obj  = data; 
    for (; idx < nTokens; idx++) {
    }
    return true;
  }
#endif
  /** 
   * Print information to standard output
   * 
   * @param option 
   */
  virtual void Print(Option_t* option="") const 
  {
    Helper::Print(option);
    std::cout << std::boolalpha 
	      << "  --- Other settings -------\n"
	      << "  Extra libraries  : " << fExtraLibs << "\n"
	      << "  Extra PARs       : " << fExtraPars << "\n"
	      << "  Extra sources    : " << fExtraSrcs << "\n"
	      << "  Use PARs of tasks: " << fUsePars   << "\n"
	      << "  Use PARs of base : " << fBasePars  
	      << std::noboolalpha << std::endl;
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
    Bool_t ret = Helper::AuxFile(name, copy);
    if (!name.BeginsWith("/")) {
      fAuxFiles.Add(new TObjString(name));
    }
#if 0
    if (ret && name.EndsWith(".root")) { 
      TFile* file = TFile::Open(name, "READ");
      if (file) {
	Info("AuxFile", "Adding input file %s", name.Data());
	gProof->AddInputData(file, true);
      }
    }
#endif
    return ret;
  }
  Int_t SendFile(const TString& fileName) 
  {
    Int_t    bufSize = 32768;
    Char_t   buf[bufSize];
    Long64_t size = 0;
    Long_t   id = 0, flags = 0, modtime = 0;
    if (gSystem->GetPathInfo(fileName.Data(), &id, &size, &flags, &modtime)==1 
	|| size <= 0) {
      Error("SendFile", "Cannot stat %s", fileName.Data());
      return -1;
    }
    TString fn(gSystem->BaseName(fileName.Data()));
    TList*  slaves = 0; // gProof->GetListOfActiveSlaves(); - protected
    TIter   next(slaves);
    TSlave* sl   = 0;
    Int_t   ret  = 0;
    Int_t   fd = open(fileName.Data(), O_RDONLY);
    while ((sl = static_cast<TSlave*>(next()))) {
      if (!sl->IsValid()) continue;
      if (sl->GetSlaveType() != TSlave::kSlave) continue;
      
      // Always binary (first 1), never forward (last 0).
      snprintf(buf,bufSize,"%s %d %lld %d", fn.Data(), 1, size, 0);
      if (sl->GetSocket()->Send(buf, kPROOF_SENDFILE) == -1) {
	Warning("SendFile", "Could not send kPROOF_SENDFILE request");
	continue;
      }

      // Go to the beginning of the file 
      lseek(fd, 0, SEEK_SET);
      Int_t len = 0;
      do { 
	while ((len = read(fd, buf, bufSize)) < 0 && 
	       TSystem::GetErrno() == EINTR)
	  TSystem::ResetErrno();
	if (len < 0) { 
	  Error("SendFile", "error reading input");
	  close(fd);
	  return -1;
	}
	if (len > 0 && sl->GetSocket()->SendRaw(buf, len) == -1) {
	  Error("SendFile", "error writing to slave");
	  sl = 0;
	  break;
	}
      } while (len > 0);
      ret ++;

      // Wait for slave - private
      // if (sl) gProof->Collect(sl,gEnv->GetValue("Proof.CollectTimeout",-1));
    }

    // Close the file 
    close(fd);

    return ret;
  }
  /** 
   * Path of output 
   * 
   * @return Path to output - possibly a data set
   */
  virtual TString OutputPath() const 
  {
    TString ret;
    if (fOptions.Has("dsname")) {
      ret = Form("/%s/%s/", gProof->GetGroup(), gProof->GetUser());
      ret.Append(OutputUtilities::RegisteredDataset());
    }
    return ret;
  }
  /** 
   * @return URL help string
   */
  virtual const Char_t* UrlHelp() const 
  {
    return "proof://<host>[:<port>]/[<dataset>|<path>][?<options>][#<treeName>]";
  }
  /** 
   * @return Short description 
   */
  virtual const char* Desc() const { return "PROOF"; }
  TString fExtraLibs;
  TString fExtraPars;
  TString fExtraSrcs;
  Bool_t  fUsePars;
  Bool_t  fBasePars;
  TList   fAuxFiles;
};
#endif
//
// EOF
//
