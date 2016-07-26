/**
 * @file   ProofRailway.C
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
#include "Railway.C"
#ifndef __CINT__
# include "OutputUtilities.C"
# include "ParUtilities.C"
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
struct ProofRailway : public Railway
{
  /** 
   * Constructor 
   * 
   * @param url  Url 
   * @param verbose Verbosity level
   */
  ProofRailway(const TUrl& url, Int_t verbose)
    : Railway(url, verbose), 
      fExtraLibs(""),
      fExtraPars(""),
      fExtraSrcs(""),
      fUsePars(false), 
      fBasePars(false),
      fTestBuild(true),
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
    fOptions.Add("testpar",  "Test build PARs");
    if (!fUrl.GetUser() || fUrl.GetUser()[0] == '\0') 
      fUrl.SetUser(gSystem->GetUserInfo()->fUser);
    fAuxFiles.SetOwner();
  }
  ProofRailway(const ProofRailway& o) 
    : Railway(o),
      fExtraLibs(""),
      fExtraPars(""),
      fExtraSrcs(""),
      fUsePars(false), 
      fBasePars(false),
      fTestBuild(true),
      fAuxFiles()
  {}
  ProofRailway& operator=(const ProofRailway& o) 
  {
    if (&o == this) return *this;
    Railway::operator=(o);
    fExtraLibs = o.fExtraLibs;
    fExtraPars = o.fExtraPars;
    fExtraSrcs = o.fExtraSrcs;
    fUsePars   = o.fUsePars;    
    fBasePars  = o.fBasePars;
    fTestBuild = o.fTestBuild;
    // fAuxFiles;
    return *this;
  }
  /** 
   * Destructor 
   */
  virtual ~ProofRailway() {}
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
    gProof->AddIncludePath(p);
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
    if ((!fUsePars || isBase) && !forcePar) {
      Int_t ret = gSystem->Load(MakeLibraryName(name));
      if (ret < 0) return false;
      if (slaves) fExtraLibs.Append(Form(":%s", name.Data()));
    }
    else { 
      if (!ParUtilities::Find(name)) { 
	Error("ProofRailway::LoadLibrary", "Failed to find PAR file %s", 
	      name.Data());
	return false;
      }
      if (fTestBuild && !ParUtilities::Build(name)) { 
	Error("ProofRailway::LoadLibrary", "Failed to build PAR file %s", 
	      name.Data());
	return false;
      }
      if (gProof->UploadPackage(name.Data(), TProof::kRemoveOld) < 0) {
	Error("ProofRailway::LoadLibrary", "Failed to upload PAR file %s", 
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
  virtual Bool_t LoadSource(const TString& name, bool copy=true)
  {
    if (!Railway::LoadSource(name, copy)) return false;
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
      Error("ProofRailway::LoadAliROOT", "Local AliROOT not available");
      return false;
    }

    Bool_t tmp = fUsePars;
    fUsePars   = fBasePars;
    if (!LoadLibrary("STEERBase"))     return false;
    if (!LoadLibrary("ESD"))           return false;
    if (!LoadLibrary("AOD"))           return false;
    if (!LoadLibrary("ANALYSIS"))      return false;
    if (!LoadLibrary("ANALYSISalice")) return false;
    fUsePars = tmp;

    return CreateAliROOTPar();
  }
  /** 
   * Set-up to load the AliPHYSICS libraries 
   * 
   * @return true on success
   */
  virtual Bool_t LoadAliPhysics()
  {
    if (!gSystem->Getenv("ALICE_PHYSICS")) { 
      Error("ProofRailway::LoadAliPhysics",
	    "Local AliPhysics not available");
      return false;
    }

    Bool_t tmp = fUsePars;
    fUsePars   = fBasePars;
    if (!LoadLibrary("OADB"))          return false;
    fUsePars = tmp;

    return CreateAliPhysicsPar();
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
   * Get the name of the AliPHYSICS par file to use 
   * 
   * @return String 
   */
  virtual const char* AliPhysicsParName() const
  {
    return "ALIPHYSICS";
  }
  /** 
   * Create an AliROOT/AliPhysics par file from the executing AliROOT.
   * This PAR file basically uses the environment of the client - that
   * is, we assume that the used AliROOT is accessible on the slaves -
   * e.g., via an NFS export.
   * 
   * Note, the SETUP.C script take one argument - a TList of TNamed
   * parameters.  Parameters processed are	
   *
   * - ALIROOT_MODE=[default,aliroot,rec,sim,train]
   *   - default: Load base analysis libraries 
   *   - aliroot: Load $ALICE_ROOT/macros/loadralibs.C
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
  virtual Bool_t CreatePseudoPar(const TString& parName,
				 const TString& env,
				 const TString& setup)

  {
    if (fBasePars) return true;

    TString parFile(Form("%s.par", parName.Data()));

    // --- Check if we have the drirectory already -------------------
    if (gSystem->AccessPathName(parName.Data()) == 0) { 
      // Let's remove it to get a clean slate 
      if (gSystem->Exec(Form("rm -rf %s", parName.Data())) != 0) {
	Error("ProofRailway", "Failed to remove %s", parName.Data());
	return false;
      }
    }
    // --- Check if the PAR file is there, and remove it if so -------
    if (gSystem->AccessPathName(parFile.Data()) == 0) { 
      if (gSystem->Unlink(parFile.Data()) != 0) { 
	Error("ProofRailway::CreatePseudoPar", "Failed to remove %s", 
	      parFile.Data());
	return false;
      }
    }
      

    // Set-up directories 
    if (gSystem->MakeDirectory(parName) < 0) {
      Error("ProofRailway::CreatePseudoPar", "Could not make directory '%s'", 
	    parName.Data());
      return false;
    }
    
    if (gSystem->MakeDirectory(Form("%s/PROOF-INF", parName.Data()))) {
      Error("ProofRailway::CreatePseudoPar", 
	    "Could not make directory %s/PROOF-INF", 
	    parName.Data());
      return false;
    }

    // --- Build script does nothing ---------------------------------
    std::ofstream b(Form("%s/PROOF-INF/BUILD.sh",parName.Data()));
    if (!b) { 
      Error("ProofRailway::CreatePseudoPar", 
	    "Failed to make BUILD.sh shell script");
      return false;
    }
    b << "#!/bin/sh\n\n"
      << "# echo Nothing to do\n"
      << "exit 0\n"
      << std::endl;
    b.close();
    gSystem->Exec(Form("chmod a+x %s/PROOF-INF/BUILD.sh",parName.Data()));

    // --- Possible environment script -------------------------------
    TString envScript = env;
    if (!envScript.EndsWith(".C"))
      envScript = "";
    if (!envScript.IsNull()) {
      // If an environment script was specified, copy that to the par
      if (gSystem->AccessPathName(envScript.Data()) == 0) { 
	// Copy script 
	if (gSystem->Exec(Form("cp %s %s/PROOF-INF/", envScript.Data(), 
			       parName.Data())) != 0) {
	  Error("ProofRailway", "Failed to copy %s", envScript.Data());
	  return false;
	}
      }
      else {
	Warning("CreateALIROOTPar", "Couldn't read %s", envScript.Data());
	envScript = "";
      }
    }
    // --- Write SETUP script ----------------------------------------
    std::ofstream s(Form("%s/PROOF-INF/SETUP.C", parName.Data()));
    if (!s) { 
      Error("ProofRailway::CreatePseudoPar", 
	    "Failed to make SETUP.C ROOT script");
      return false;
    }
    s << "void SETUP(TList* opts) {\n";
    if (envScript.IsNull()) {
      s << env << std::endl;
    }
    else { 
      s << "  gROOT->Macro(\"PROOF-INF/" << gSystem->BaseName(envScript.Data())
	<< "\");\n";
    }
    s  << setup << "}\n"
      << std::endl;
    s.close();

    // --- TAR up everything -----------------------------------------
    Int_t ret = gSystem->Exec(Form("tar -czf %s %s",
				   parFile.Data(), parName.Data()));
    if (ret != 0) { 
      Error("ProofRailway::CreatePseudoPar", "Failed to pack up PAR file %s",
	    parFile.Data());
      return false;
    }

    // --- And upload the package ------------------------------------
    ret = gProof->UploadPackage(parFile.Data(),TProof::kRemoveOld);
    if (ret != 0) { 
      Error("ProofRailway::CreatePseudoPar", 
	    "Failed to upload the AliROOT PAR file");
      return false;
    }
    // Note, the PAR isn't enabled until much later when we've
    // collected all the needed libraries in fExtraLibs
    return true;
  }
  static void ExportEnvVar(TString& out, const TString& name)
  {
    TString env(gSystem->Getenv(name.Data()));
    if (env.IsNull()) return;

    out.Append(Form("  gSystem->Setenv(\"%s\",\"%s\");\n",
		    name.Data(), env.Data()));
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
    TString parName(AliROOTParName());
    TString envScript = fOptions.Get("env");
    if (envScript.EqualTo("-none-", TString::kIgnoreCase) ||
	!envScript.EndsWith(".C"))
      envScript = "";

    TString env       = "";
    ExportEnvVar(env, "ALICE");
    ExportEnvVar(env, "ALICE_ROOT");
    
    TString setup("gSystem->AddDynamicPath(\"$(ALICE_ROOT)/lib\");\n"
                  "gSystem->AddIncludePath(\"-I${ALICE_ROOT}/include\");\n"
		  "// Info(\"SETUP\",\"Loading ROOT libraries\");\n"
		  "gSystem->Load(\"libTree\");\n"
		  "gSystem->Load(\"libGeom\");\n"
		  "gSystem->Load(\"libVMC\");\n"
		  "gSystem->Load(\"libPhysics\");\n"
		  "gSystem->Load(\"libMinuit\");\n"
		  "\n"
		  "// Info(\"SETUP\",\"Parameter list:\");\n"
		  "if (!opts) return;\n"
		  "//opts->ls();\n"
		  "\n"
		  "TObject* par = opts->FindObject(\"ALIROOT_MODE\");\n"
		  "if (par) {\n"
		  " //Info(\"SETUP\",\"ALIROOT mode: %s\",par->GetTitle());\n"
		  "  TString mode(par->GetTitle());mode.ToLower();\n"
		  "  if (mode.EqualTo(\"default\")) {\n"
		  "    gSystem->Load(\"libSTEERBase\");\n"
		  "    gSystem->Load(\"libESD\");\n"
		  "    gSystem->Load(\"libAOD\");\n"
		  "    gSystem->Load(\"libANALYSIS\");\n"
		  "    gSystem->Load(\"libANALYSISalice\");\n"
		  "  }\n"
		  "  else if (mode.EqualTo(\"aliroot\")) \n"
		  "    gROOT->Macro(\"$ALICE_ROOT/macros/loadlibs.C\");\n"
		  "  else if (mode.EqualTo(\"rec\")) \n"
		  "    gROOT->Macro(\"$ALICE_ROOT/macros/loadlibsrec.C\");\n"
		  "  else if (mode.EqualTo(\"sim\")) \n"
		  "    gROOT->Macro(\"$ALICE_ROOT/macros/loadlibssim.C\");\n"
		  "  else if (mode.EqualTo(\"train\")) \n"
		  "   gROOT->Macro(\"$ALICE_ROOT/macros/loadlibstrain.C\");\n"
		  "  else if (mode.EqualTo(\"custom\")) \n"
		  "   gROOT->Macro(\"$ALICE_ROOT/macros/loadlibstrain.C\");\n"
		  "}\n"
		  "\n"
		  "par = opts->FindObject(\"ALIROOT_EXTRA_LIBS\");\n"
		  "if (par) {\n"
		  "  Info(\"SETUP\",\"Libaries to load: %s\\n\",\n"
		  "       par->GetTitle());\n"
		  "  TString tit(par->GetTitle());\n"
		  "  TObjArray* tokens = tit.Tokenize(\":\");\n"
		  "  TObject*   lib    = 0;\n"
		  "  TIter      next(tokens);\n"
		  "  while ((lib = next())) {\n"
		  "    TString libName(lib->GetName());\n"
		  "    if (!libName.BeginsWith(\"lib\"))\n"
		  "      libName.Prepend(\"lib\");\n"
		  "    // Info(\"SETUP\",\"Loading %s ...\",libName.Data());\n"
		  "    gSystem->Load(libName.Data());\n"
		  "  }\n"
		  "}");
    // Note, the PAR isn't enabled until much later when we've
    // collected all the needed libraries in fExtraLibs
    return CreatePseudoPar(parName,envScript.IsNull() ? env : envScript,setup);
  }
  virtual Bool_t CreateAliPhysicsPar()
  {
    TString parName(AliPhysicsParName());
    TString env       = "";
    ExportEnvVar(env, "ALICE_PHYSICS");
    ExportEnvVar(env, "OADB_PATH");
    
    TString setup("gSystem->AddDynamicPath(\"$(ALICE_PHYSICS)/lib\");\n"
                  "gSystem->AddIncludePath(\"-I${ALICE_PHYSICS}/include\");\n"
		  "// Info(\"SETUP\",\"Parameter list:\");\n"
		  "if (!opts) return;\n"
		  "//opts->ls();\n"
		  "\n"
		  "TObject* par = opts->FindObject(\"ALIPHYSICS_MODE\");\n"
		  "if (par) {\n"
		  " //Info(\"SETUP\",\"ALIPHYSICS mode:%s\",par->GetTitle());\n"
		  "  TString mode(par->GetTitle());mode.ToLower();\n"
		  "  if (mode.EqualTo(\"default\")) {\n"
		  "    gSystem->Load(\"libOADB\");\n"
		  "  }\n"
		  "}\n"
		  "\n"
		  "par = opts->FindObject(\"ALIPHYSICS_EXTRA_LIBS\");\n"
		  "if (par) {\n"
		  "  Info(\"SETUP\",\"Libaries to load: %s\\n\",\n"
		  "       par->GetTitle());\n"
		  "  TString tit(par->GetTitle());\n"
		  "  TObjArray* tokens = tit.Tokenize(\":\");\n"
		  "  TObject*   lib    = 0;\n"
		  "  TIter      next(tokens);\n"
		  "  while ((lib = next())) {\n"
		  "    TString libName(lib->GetName());\n"
		  "    if (!libName.BeginsWith(\"lib\"))\n"
		  "      libName.Prepend(\"lib\");\n"
		  "    // Info(\"SETUP\",\"Loading %s ...\",libName.Data());\n"
		  "    gSystem->Load(libName.Data());\n"
		  "  }\n"
		  "}");
    // Note, the PAR isn't enabled until much later when we've
    // collected all the needed libraries in fExtraLibs
    return CreatePseudoPar(parName,env,setup);
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
   * Connect to PROOF(,-Lite,-on-Demand) cluster 
   * 
   * @param url Connection URL
   * @param opt Possible options 
   * 
   * @return true on success 
   */
  virtual Bool_t Connect(const TUrl& url, const TString& opts)
  {
    Info("ProofRailway::Connect", "Connecting to %s with %soptions %s", 
	 url.GetUrl(), 
	 opts.IsNull() ? "no " : "", 
	 opts.Data());
    TProof::Open(url.GetUrl(), opts);
    // TProof::Open(connect.GetHost(), opts);
    if (!gProof) { 
      Error("ProofRailway::Connect", "Failed to open Proof connection %s", 
	    url.GetUrl());
      return false;
    }
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

    // --- Add ALICE_ROOT directory to search path for packages ----
    // Info("ProofRailway::PreSetup", "Set location of packages");
    gEnv->SetValue("Proof.GlobalPackageDirs", 
		   Form("%s:%s:%s", 
			gEnv->GetValue("Proof.GlobalPackageDirs", "."),
			gSystem->Getenv("ALICE_PHYSICS"),
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
      Info("ProofRailway::PreSetup", "Doing a %s reset of %s", 
	   hard ? "hard" : "soft", connect.GetUrl());
      TProof::Reset(connect.GetUrl(), hard);
      Int_t secs = 3;
      Info("ProofRailway::PreSetup", 
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
      Info("ProofRailway::PreSetup", "Using wrapper command: %s", 
	   wrapper.Data());
      TProof::AddEnvVar("PROOF_WRAPPERCMD", wrapper);
    }

    // --- PAR parameters --------------------------------------------
    fUsePars   = fOptions.Has("par");
    fBasePars  = (fUsePars && 
		  fOptions.Get("par").EqualTo("all",TString::kIgnoreCase));
    fTestBuild = fOptions.Has("testpar");

    // --- Connect to the cluster ------------------------------------
    TString opts;
    if (fOptions.Has("workers")) 
      opts.Append(Form("workers=%s", fOptions.Get("workers").Data()));
      
    if (!Connect(connect, opts)) return false;
    
    // --- Check if we need to clear packages ------------------------
    if (fOptions.Has("clear")) {
      TString pkgs = fOptions.Get("clear");
      if (pkgs.IsNull() || pkgs.EqualTo("all", TString::kIgnoreCase)) { 
	// No value given, clear all 
	if (gProof->ClearPackages() != 0) 
	  Warning("ProofRailway::PreSetup", "Failed to lear all packages");
      }
      else { 
	// Tokenize on ',' and clear each package 
	TObjArray* pars = pkgs.Tokenize(",");
	TObject*   pkg  = 0;
	TIter      next(pars); 
	while ((pkg = next())) { 
	  if (gProof->ClearPackage(pkg->GetName()) != 0)
	    Warning("ProofRailway::PreSetup", "Failed to clear package %s", 
		    pkg->GetName());
	}
	pars->Delete();
      }
    }
    return true;
  }
  virtual Bool_t EnableSpecial(const TString& parName,
			       const TString& prefix)
  {
    if (prefix.IsNull()) {
      Warning("EnableSpecial", "No prefix specified");
      return false;
    }
    const char* prf = prefix.Data();
    if (parName.IsNull()) {
      Warning("EnableSpecial", "No par name specified for %s", prf);
      return false;
    }
    // List of parameters for PAR x
    TList* params = new TList;
    params->SetOwner(true);
    
    // Extra libraries 
    TString tmp(fExtraLibs.Strip(TString::kBoth,':'));
    params->Add(new TNamed(Form("%s_EXTRA_LIBS", prf), tmp.Data()));

    // Check for mode 
    if (fOptions.Has("mode"))
      params->Add(new TNamed(Form("%s_MODE", prf),
			     fOptions.Get("mode").Data()));
    else
      params->Add(new TNamed(Form("%s_MODE", prf), "default"));

    // Check for AliEn 
    if (fOptions.Has("alien"))
      params->Add(new TNamed(Form("%s_ENABLE_ALIEN", prf), "1"));

    // Try to enable the package - do we need to load first? 
    Int_t ret = gProof->EnablePackage(parName.Data(), params, true);
    if (ret < 0) {
      Error("ProofRailway::EnableSpecial", "Failed to enable %s PAR %s", 
	    prf, parName.Data());
      return false;
    }
    return true;
  }
      
  /** 
   * Enable the special AliROOT package on the cluster 
   * 
   * @return true on success
   */
  virtual Bool_t EnableAliROOT()
  {
    return EnableSpecial(AliROOTParName(), "ALIROOT");
  }
  /** 
   * Enable the special AliROOT package on the cluster 
   * 
   * @return true on success
   */
  virtual Bool_t EnableAliPhysics()
  {
    return EnableSpecial(AliPhysicsParName(), "ALIPHYSICS");
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
      Error("ProofRailway::PostSetup", "No analysis manager defined");
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
      if (!EnableAliROOT()) return false;
      if (!EnableAliPhysics()) return false;
    }

    // --- Make PAR file of Aux Files --------------------------------
    if (fAuxFiles.GetEntries() > 0) { 
      TString name = TString::Format("%s_auxfiles", mgr->GetName());
      ParUtilities::MakeAuxFilePAR(fAuxFiles, name);

      if (gProof->UploadPackage(name.Data(), TProof::kRemoveOld) < 0) 
	Error("ProofRailway::PostSetup", "Failed to upload PAR file %s", 
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
	Error("ProofRailway::PostSetup", "Failed to enable PAR %s",
	      obj->GetName());
	return false;
      }
    }
    
    // --- Load extra sources ----------------------------------------
    return LoadExtraSrcs();
  }
  virtual Bool_t LoadExtraSrcs()
  {
    TString    tmp2 = fExtraSrcs.Strip(TString::kBoth, ':');
    TObjArray* srcs = tmp2.Tokenize(":");
    TIter      next2(srcs);
    TObject*   obj = 0;
    while ((obj = next2())) { 
      Int_t ret = gProof->Load(Form("%s+g", obj->GetName()), true);
      if (ret < 0) { 
	Error("ProofRailway::PostSetup", "Failed to compile %s",obj->GetName());
	return false;
      }
    }
    return true;
  }    
  /** 
   * Get the data-set name 
   * 
   * @param dsname On return, must contain the data set name 
   */
  virtual void GetDataSet(TString& dsname)
  {
    dsname = fUrl.GetFile();
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
    TString dsName;
    GetDataSet(dsName);
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

    TString store = fOptions.Get("storage");
    if (!store.IsNull() && store.EqualTo("auto",TString::kIgnoreCase))
      OutputUtilities::StopXrootd();
    
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
    Railway::Print(option);
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
  virtual Bool_t AuxFile(TString& name, bool copy=false)
  {
    TString local = name;
    Bool_t ret = Railway::AuxFile(local, copy);
    if (!name.BeginsWith("/")) {
      fAuxFiles.Add(new TObjString(local));
    }
    name = local;
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
    return
      "proof://<host>[:<port>]/[<dataset>|<path>][?<options>][#<treeName>]";
  }
  /** 
   * @return Short description 
   */
  virtual const char* Desc() const { return "PROOF"; }
  /** 
   * Write auxillary ROOT (and possible shell) script for more 
   * (post-)processing e.g., terminate
   * 
   * @param escaped        Escaped name  
   */
  void AuxSave(const TString& escaped, 
	       Bool_t /*asShellScript*/) 
  {
    // Write train to file - for terminate actions 
    TFile* mgr = TFile::Open(Form("%s.root", escaped.Data()), "RECREATE");
    AliAnalysisManager::GetAnalysisManager()->Write();
    mgr->Write();

    TString     libs = fExtraLibs; libs.ReplaceAll(":", " ");
    TString     pars = fExtraPars; pars.ReplaceAll(":", " ");
    TString     srcs = fExtraSrcs; srcs.ReplaceAll(":", " ");
    TString     macDir("$ALICE_PHYSICS/PWGLF/FORWARD/trains");
    std::ofstream t("Terminate.C");
    if (!t) { 
      Error("ProofRailway::AuxSave", "Failed to make terminate ROOT script");
      return;
    }

    t << "// Generated by ProofRailway\n"
      << "Bool_t Terminate()\n"
      << "{\n"
      << "  TString name = \"" << escaped << "\";\n"
      << "  TString libs = \"" << libs << "\";\n"
      << "  TString pars = \"" << pars << "\";\n"
      << "  TString srcs = \"" << srcs << "\";\n\n"
      << "  gSystem->Load(\"libANALYSIS\");\n"
      << "  gSystem->Load(\"libANALYSISalice\");\n"
      << "  gSystem->AddIncludePath(\"-I$ALICE_ROOT/include\");\n\n"
      << "  gSystem->AddIncludePath(\"-I$ALICE_PHYSICS/include\");\n\n"
      << "  gROOT->LoadMacro(\"" << macDir << "/ParUtilities.C+g\");\n\n"
      << "  gROOT->LoadMacro(\"" << macDir << "/ProofTerminate.C+g\");\n\n"
      << "  return ProofTerminate(name,libs,pars,srcs);\n"
      << "}\n"
      << "// EOF\n"
      << std::endl;
    t.close();    
  }
  TString fExtraLibs;
  TString fExtraPars;
  TString fExtraSrcs;
  Bool_t  fUsePars;
  Bool_t  fBasePars;
  Bool_t  fTestBuild;
  TList   fAuxFiles;
};
#endif
//
// EOF
//
