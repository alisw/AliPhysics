/**
 * @file   ParUtilities.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Tue Oct 16 17:51:10 2012
 * 
 * @brief  PAR file utilities 
 * 
 * @ingroup pwglf_forward_trains_util
 * 
 */
#ifndef PARHELPER_C
#define PARHELPER_C
#ifndef __CINT__
# include <TString.h>
# include <TProof.h>
# include <TSystem.h>
# include <TError.h>
# include <TFile.h>
# include <TSystem.h>
# include <TROOT.h>
# include <fstream>
# include <cstdlib>
#else
class TString;
#endif

// ===================================================================
/**
 * Helper to set-up and load PARs
 *
 * @ingroup pwglf_forward_trains_util
 */
struct ParUtilities
{
  /** 
   * Find PAR file (either in current or parent directory or directly 
   * in $ALICE_ROOT), and link it here
   * 
   * @param what PAR file name (sans .par)
   * 
   * @return true on success
   */
  static Bool_t Find(const TString& what)
  {
    if (what.IsNull()) return false;
    
    TString parFile(what);
    if (!parFile.EndsWith(".par")) parFile.Append(".par");
    if (gSystem->AccessPathName(parFile.Data())) { 
      // If not found
      TString src;
      if (gSystem->AccessPathName(Form("../%s.par", parFile.Data())) == 0) 
	src.Form("../%s", parFile.Data());
      else {
	// If not found 
	TString aliParFile = 
	  gSystem->ExpandPathName(Form("$(ALICE_ROOT)/%s", parFile.Data()));
	if (gSystem->AccessPathName(aliParFile.Data()) == 0) 
	  src = aliParFile;
      }
      if (src.IsNull()) {
	  Error("ParUtilities::Find", 
		"PAR file %s not found in current or parent "
		"directory nor in $(ALICE_ROOT)", parFile.Data());
	  return false;
      }
      // Copy to current directory 
      // TFile::Copy(aliParFile, parFile);
      if (gSystem->Exec(Form("ln -s %s %s", src.Data(), parFile.Data())) != 0){
	Error("ParUtilities::Find", "Failed to symlink %s to %s", 
	      src.Data(), parFile.Data());
	return false;
      }
    }
    return true;
  }
  /** 
   * Unpack and load a PAR file previously found with Find.
   * 
   * @param name PAR file name 
   * @deprecated Use Find  and Build instead
   * @return true on success
   */
  static Bool_t Load(const TString& name) 
  {
    if (name.IsNull()) return true;
    if (!gProof) { 
      Error("ParUtilities::Load", "No connection to a Proof cluster");
      return false;
    }

    // Load par library 
    TString fn(name);
    Info("ParUtilities::LoadLibrary", "Uploading %s", name.Data());
    
    // First check in current directory
    Int_t ret = gProof->UploadPackage(fn, TProof::kRemoveOld);
    
    if (ret < 0)  {
      // IF not found there, then check parent directory 
      fn.Prepend("../");
      gSystem->ExpandPathName(fn);
      ret = gProof->UploadPackage(fn);
    }

    if (ret < 0) {	
      // If not found in current or parent directory, try the 
      // the ALICE_ROOT directory 
      fn  = Form("$ALICE_ROOT/%s.par", name.Data());
      gSystem->ExpandPathName(fn);
      ret = gProof->UploadPackage(fn);
    }
    
    if (ret < 0) {
      // IF not found, bark 
      Error("ParUtilities::Load", 
	    "Could not find module %s.par in current or parent directory "
	    "nor in $ALICE_ROOT", name.Data());
      return false;
    }
    
    ret = gProof->EnablePackage(name);
    Info("ParUtilities::Load", "Enabled package %s (from %s)", 
	 name.Data(), fn.Data());
    
    return true;
  }
  /** 
   * Unpack, build, and load a PAR file. 
   * 
   * @param what Which PAR file 
   * 
   * @return 
   */
  static Bool_t Build(const TString& what)
  {
    if (what.IsNull()) return false;
    
    TString parFile(what);
    if (!parFile.EndsWith(".par")) parFile.Append(".par");

    // Extract archive 
    gSystem->Exec(Form("tar xzf %s", parFile.Data()));
    
    // Change directory into par archive
    TString cwd = gSystem->WorkingDirectory();
    
    TString dir(what);
    if (dir.EndsWith(".par")) dir.ReplaceAll(".par", "");
    if (!gSystem->ChangeDirectory(dir)) { 
      Error("ParUtilities::Setup", "Failed to change directory to %s", 
	    dir.Data());
      return false;
    }
    
    // Test the build 
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      Info("ParUtilities::Setup", "Building in PAR archive %s", parFile.Data());
      if (gSystem->Exec("PROOF-INF/BUILD.sh")) { 
	Error("ParUtilities::Setup", "Failed to build in PAR directory %s", 
	      dir.Data());
	gSystem->ChangeDirectory(cwd.Data());
	return false;
      }
    }
    
    // Check for setup script
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      // Info("ParUtilities::SetupPAR", "Setting up for PAR %s", what);
      gROOT->Macro("PROOF-INF/SETUP.C");
    }
    if (!gSystem->ChangeDirectory(cwd.Data())) return false;

    return true;
  }
  //__________________________________________________________________
  /** 
   * @{ 
   * @name PAR generation from script 
   */
  /** 
   * Service function to make a PAR out of a script.  
   * 
   * The script should contain can contain a sub-class of AliAnalysisTask. 
   * The script will be compiled on the slaves before loading the 
   * AliAnalysisManager.  Parts to (not) be compiled can be protected like 
   * 
   * @code 
   * #ifdef BUILD_PAR
   * // This will _only_ be compiled in the servers 
   * #endif
   * #ifndef BUILD_PAR
   * // This will not be compiled in the servers 
   * #endif
   * @endcode
   * 
   * @param script Script to upload and compile in the PAR
   * @param deps   Dependency pars 
   * @param isLocal Local build 
   * 
   * @return true on success. 
   */
  static Bool_t MakeScriptPAR(Bool_t isLocal, 
			      const TString& script, 
			      const TString& deps)
  {
    // --- In local mode, just AcLic and load ------------------------
    if (isLocal) { 
      if (gROOT->LoadMacro(Form("%s++g", script.Data())) < 0)
	return false;
      return true;
    }

    // --- Get the base name -----------------------------------------
    Info("ParUtilities::MakeScriptPAR", "Making par file for %s", 
	 script.Data());
    TString base(gSystem->BaseName(script));
    Int_t   idx = base.Last('.');
    if (idx != kNPOS) base.Remove(idx);

    // --- Check name of script file ---------------------------------
    TString scr(script);
    TString ext;
    if      (script.EndsWith(".C"))   ext = "C"; 
    else if (script.EndsWith(".cxx")) ext = "cxx";
    else                              { ext = "C"; scr.Append(".C"); }
    
    // --- Check if we can access the file ---------------------------
    TString path = TString::Format(".:%s", TROOT::GetMacroPath());
    char*   loc  = gSystem->Which(path, scr);
    if (!loc) {
      Error("ParUtilities::MakeScriptPAR", 
	    "Script %s not found in %s", scr.Data(), path.Data());
      return false;
    }
    TString full(loc);

    // --- Create our temporary directory ----------------------------
    TString tmpdir(gSystem->TempDirectory());
    int     ltempl = tmpdir.Length() + 1 + 5 + 6 + 1;
    char*   templ  = new char[ltempl];
    snprintf(templ, ltempl, "%s/trainXXXXXX", tmpdir.Data());
    if (!mkdtemp(templ)) {
      Error("ParUtilities::MakeScriptPAR", 
	    "Failed to generate temporary directory from template %s", 
	    templ);
      return false;
    }

    Bool_t retVal = false;
    try {
      // --- Make directories for package ------------------------------
      TString dir = TString::Format("%s/%s", templ, base.Data());
      // Set-up directories 
      if (gSystem->MakeDirectory(dir) < 0) 
	throw TString::Format("Could not make directory '%s'", base.Data());
      if (gSystem->MakeDirectory(Form("%s/PROOF-INF", dir.Data()))) 
	throw TString::Format("Could not make directory %s/PROOF-INF", 
			      base.Data());
      
      // --- Copy the script to the setup directory --------------------
      TString dest = TString::Format("%s/%s.%s", dir.Data(),
				     base.Data(), ext.Data());
      Int_t ret = gSystem->CopyFile(full, dest, true);
      switch (ret) { 
      case -1: throw TString::Format("Couldn't open %s for copy", scr.Data());
      case -2: throw TString::Format("File %s exists", dest.Data());
      case -3: throw TString::Format("Error while copying %s", scr.Data());
      }
      
      // --- Make scripts, etc. ----------------------------------------
      TObjArray* depList = deps.Tokenize(", ");
      if (!MakeBuildScript(dir, base)) 
	throw TString::Format("Failed to make build script");
      if (!MakeUtilityScript(dir)) 
	throw TString::Format("Failed to make utility script");
      if (!MakeBuildMacro(dir, base, ext, depList)) 
	throw TString::Format("Failed to make build macro");
      if (!MakeSetupMacro(dir, base, ext, depList)) 
	throw TString::Format("Failed to setup macro");

      // --- Pack up the archive ---------------------------------------
      ret = gSystem->Exec(Form("(cd %s && tar -czf %s.par %s)", 
			       templ, base.Data(),base.Data()));
      if (ret != 0) 
	throw TString::Format("Failed to create PAR file %s.PAR from %s", 
			      base.Data(), dir.Data());

      // --- Move PAR file to here -------------------------------------
      ret = gSystem->Exec(Form("mv -f %s/%s.par %s.par", templ, base.Data(), 
			       base.Data()));
      if (ret != 0) 
	throw TString::Format("Failed to rename %s/%s.par to %s.par: %s", 
			      templ, base.Data(), base.Data(), 
			      gSystem->GetError());
      retVal = true;
    }
    catch (TString& e) {
      Error("ParUtilities::MakeScriptPAR", "%s", e.Data());
      retVal = false;
    }
    
    // --- Remove temporary directory --------------------------------
    gSystem->Exec(Form("rm -rf %s", templ));

    return retVal;
  }
  /** 
   * Write a build script
   * 
   * @param dir Directory to put it in
   * @param base Base name 
   * 
   * @return true on success
   */
  static Bool_t MakeBuildScript(const TString& dir, 
				const TString& base)
  {
    // Make our build file 
    std::ofstream out(Form("%s/PROOF-INF/BUILD.sh", dir.Data()));
    if (!out) {
      Error("ParUtilities::MakeBuildScript", "Failed to open out shell script");
      return false;
    }
    out << "#!/bin/sh\n"
	<< "if test x$ALICE_ROOT != x ; then\n"
	<< "  if test x$ALICE_TARGET = x ; then\n"
	<< "    export ALICE_TARGET=`$ROOTSYS/bin/root-config --arch`\n"
	<< "  fi\n"
	<< "  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"
	<< "${ALICE_ROOT}/lib/tgt_${ALICE_TARGET}\n"
	<< "fi\n"
	<< "echo BUILD.sh@`hostname`: Building " << base << "\n"
	<< "root.exe -l -out -q PROOF-INF/BUILD.C 2>&1 | tee " 
	<< base << ".log\n"
	<< "echo BUILD.sh@`hostname`: done: $?\n"
	<< std::endl;
    out.close();
    if (gSystem->Chmod(Form("%s/PROOF-INF/BUILD.sh", dir.Data()), 0755) != 0) {
      Error("ParUtilities::MakeBuildScript", 
	    "Failed to set exectuable flags on %s/PROOF-INF/BUILD.sh", 
	    dir.Data());
      return false;
    }
    return true;
  }
  /** 
   * Write a build macro 
   * 
   * @param dir   Directory to put macro in
   * @param deps  Dependencies
   * @param base  Base name of script to compile
   * @param ext   `extension' - last part of file name 
   * 
   * @return true on success
   */
  static Bool_t MakeBuildMacro(const TString& dir, 
			       const TString& base, 
			       const TString& ext,
			       const TCollection* deps)  {
    std::ofstream out(Form("%s/PROOF-INF/BUILD.C", dir.Data()));
    if (!out) {
      Error("ParUtilities::MakeBuildMacro", "Failed to open build script");
      return false;
    }
    out << "void BUILD() {\n"
	<< "  gSystem->AddIncludePath(\"-DBUILD_PAR=1\");\n"
	<< "  gROOT->LoadMacro(\"PROOF-INF/UTIL.C\");\n"
	<< "  LoadROOTLibs();\n"
	<< "  AddAliROOT();\n";
    if (deps) {
      TIter       next(deps);
      TObject*    dep = 0;
      while ((dep = next())) {
	out << "  AddDep(\"" << dep->GetName() << "\");\t"
	    << "  LoadDep(\"" << dep->GetName() << "\");\n";
      }
    }
    out << "  // gDebug = 5;\n"
	<< "  int ret = gROOT->LoadMacro(\"" 
	<< base << "." << ext << "++g\");\n"
	<< "  if (ret != 0) Fatal(\"BUILD\",\"Failed to build\");\n"
	<< "  // else Info(\"BUILD\", \"Made " << base << "\");\n"
	<< "}\n"
	<< std::endl;
    out.close();

    return true;
  }
  /** 
   * Make a utility macro 
   * 
   * @param dir Directory to put the macro in
   * 
   * @return true on success
   */
  static Bool_t MakeUtilityScript(const TString& dir)
  {
    std::ofstream out(Form("%s/PROOF-INF/UTIL.C", dir.Data()));
    if (!out) {
      Error("ParUtilities::MakeUtilityScript", "Failed to open utility script");
      return false;
    }
    out << "void LoadROOTLibs() {\n"
	<< "  gSystem->Load(\"libVMC\");\n"
	<< "  gSystem->Load(\"libNet\");\n"
	<< "  gSystem->Load(\"libTree\");\n"
	<< "  gSystem->Load(\"libPhysics\");\n"
	<< "  gSystem->Load(\"libMinuit\");\n"
	<< "}\n\n"
	<< "void AddAliROOT() {\n"
	<< "  TString val(gSystem->Getenv(\"ALICE_ROOT\"));\n"
	<< "  if (val.IsNull())\n"
	<< "    Warning(\"Add\",\"ALICE_ROOT not defined\");\n"
	<< "  else\n"
	<< "    gSystem->AddIncludePath(Form(\"-I%s/include\",val.Data()));\n"
	<< "}\n\n"
	<< "void AddDep(const char* env) {\n"
	<< "  TString val(gSystem->Getenv(Form(\"%s_INCLUDE\",env)));\n"
	<< "  if (val.IsNull())\n"
	<< "    Warning(\"Add\",\"%s_INCLUDE not defined\", env);\n"
	<< "  else {\n"
	<< "    gSystem->AddIncludePath(Form(\"-I../%s\",val.Data()));\n"
	<< "  }\n"
	<< "}\n\n"
	<< "void LoadDep(const char* name) {\n"
	<< "  gSystem->AddDynamicPath(Form(\"../%s\",name));\n"
	<< "  char* full = gSystem->DynamicPathName(name,true);\n"
	<< "  if (!full) \n"
	<< "   full = gSystem->DynamicPathName(Form(\"lib%s\",name),true);\n"
	<< "  if (!full) \n"
	<< "   full = gSystem->DynamicPathName(Form(\"lib%s.so\",name),true);\n"
	<< "  if (!full) {\n"
	<< "    Warning(\"LoadDep\",\"Module %s not found\", name);\n"
	<< "    return;\n"
	<< "  }\n"
	<< "  gSystem->Load(full);\n"
	<< "}\n"
	<< std::endl;
    out.close();
    return true;
  }
  /** 
   * Make a setup script 
   * 
   * @param dir   Directory to put it in 
   * @param base  Base name of target script 
   * @param ext   Extension of target script 
   * @param deps  Dependencies 
   * 
   * @return true on success
   */
  static Bool_t MakeSetupMacro(const TString& dir, 
			       const TString& base, 
			       const TString& ext,
			       const TCollection* deps)
  {
    // Make our set-up script 
    std::ofstream out(Form("%s/PROOF-INF/SETUP.C", dir.Data()));
    if (!out) {
      Error("ParUtilities::MakeSetupMacro", "Failed to open setup script");
      return false;
    }
    out << "void SETUP() {\n"
	<< "  gROOT->LoadMacro(\"PROOF-INF/UTIL.C\");\n"
	<< "  LoadROOTLibs();\n"
	<< "  // Info(\"SETUP\",\"Loading libraries\");\n";
    if (deps) {
      TIter next(deps);
      TObject* dep = 0;
      while ((dep = next())) 
	out << "  LoadDep(\"" << dep->GetName() << "\");\n";
    }
    out << "  // gDebug = 5;\n"
	<< "  // Info(\"SETUP\",\"Loading " << base <<"_"<< ext << ".so\");\n"
	<< "  gSystem->Load(\"" << base << "_" << ext << ".so\");\n"
	<< "  // gDebug = 0;\n"
	<< "  gROOT->ProcessLine(\".include " << base << "\");\n"
	<< "  gSystem->Setenv(\"" << base << "_INCLUDE\",\"" 
	<< base << "\");\n"
	<< "  // Info(\"SETUP\", \"Done\");\n"
	<< "}\n"
	<< std::endl;
    out.close();
    return true;
  }
  /* @} */
};
#endif
// 
// EOF
//
