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
# include "Railway.C"
#else
class TString;
class Railway;
class TProof; // Autoload 
#endif

// ===================================================================
/**
 * Railway to set-up and load PARs
 *
 * @ingroup pwglf_forward_trains_util
 */
struct ParUtilities
{
  static Bool_t DoFind(const TString& what, TString& src)
  {
    src = "";
    if (what.IsNull()) return false;

    TString aliRoot = gSystem->ExpandPathName("$ALICE_ROOT");
    TString aliPhys = gSystem->ExpandPathName("$ALICE_PHYSICS");
    TString search = Form(".:..:%s:%s/PARfiles:%s:%s/PARfiles",
			  aliPhys.Data(), aliPhys.Data(),
 	         	  aliRoot.Data(), aliRoot.Data());
    char* found = gSystem->Which(search.Data(), what.Data());
    if (!found || found[0] == '\0') {
      Error("ParUtilities::Find", "PAR file %s not found in %s",
	    what.Data(), search.Data());
      return false;
    }
    src = found;
    return true;
  }
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
    TString parFile(what);
    if (!parFile.EndsWith(".par")) parFile.Append(".par");
    TString src;
    if (!DoFind(parFile, src)) return false;
    
    // Copy to current directory 
    // TFile::Copy(aliParFile, parFile);
    Info("", "Found PAR %s at %s", what.Data(), src.Data());
    if (gSystem->Exec(Form("ln -s %s %s", src.Data(), parFile.Data())) != 0){
      if (!gSystem->AccessPathName(parFile.Data())) return true;
      Error("ParUtilities::Find", "Failed to symlink %s to %s", 
	      src.Data(), parFile.Data());
      return false;
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

    TString parFile(name);
    if (!parFile.EndsWith(".par")) parFile.Append(".par");
    TString src;
    if (!DoFind(parFile, src)) return false;
    
    // First check in current directory
    Int_t ret = gProof->UploadPackage(src, TProof::kRemoveOld);
    
    if (ret < 0) {
      // IF not found, bark 
      Error("ParUtilities::Load", 
	    "Could not upload module %s.par", name.Data());
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

    // We need to make sure the current directory is in the load path 
    gSystem->SetDynamicPath(Form("%s:%s", gSystem->WorkingDirectory(), 
				 gSystem->GetDynamicPath()));
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
   * @param script  Script to upload and compile in the PAR
   * @param deps    Dependency pars 
   * @param isLocal Local build 
   * @param helper  Railway 
   * 
   * @return true on success. 
   */
  static Bool_t MakeScriptPAR(Bool_t         isLocal, 
			      const TString& script, 
			      const TString& deps, 
			      Railway*        helper)
  {
    TObjArray* depList = deps.Tokenize(", ");
    
    // --- In local mode, just AcLic and load ------------------------
    if (isLocal) { 
      // Load dependencies 
      TIter       next(depList);
      TObject*    dep = 0;
      while ((dep = next())) 
	helper->LoadLibrary(dep->GetName());
      
      // AcLic and load 
      Info("ParUtilities::MakeScriptPAR", "Loading macro %s", script.Data());
      if (gROOT->LoadMacro(Form("%s++g", script.Data())) < 0) {
	Error("ParUtilities::MakeScriptPAR", 
	      "Failed to build local library %s", script.Data());
	return false;
      }
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
    Info("", "Building PAR in %s", templ);

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
      Info("", "Made directory %s", dir.Data());

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
      if (!MakeScriptBuildScript(dir, base)) 
	throw TString::Format("Failed to make build script");
      if (!MakeScriptUtilityScript(dir)) 
	throw TString::Format("Failed to make utility script");
      if (!MakeScriptBuildMacro(dir, base, ext, depList)) 
	throw TString::Format("Failed to make build macro");
      if (!MakeScriptSetupMacro(dir, base, ext, depList)) 
	throw TString::Format("Failed to setup macro");

      // --- Pack up the archive ---------------------------------------
      ret = gSystem->Exec(Form("(cd %s && tar -czf %s.par %s)", 
			       templ, base.Data(),base.Data()));
      if (ret != 0) 
	throw TString::Format("Failed to create PAR file %s.PAR from %s", 
			      base.Data(), dir.Data());
      
      Info("", "Made par archive %s/%s.par - moving here", 
	   templ, base.Data());
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
  static Bool_t MakeScriptBuildScript(const TString& dir, 
				      const TString& base)
  {
    // Make our build file 
    std::ofstream out(Form("%s/PROOF-INF/BUILD.sh", dir.Data()));
    if (!out) {
      Error("ParUtilities::MakeScriptBuildScript", 
	    "Failed to open out shell script");
      return false;
    }
    out << "#!/bin/sh\n"
	<< "if test x$ALICE_ROOT != x ; then\n"
	<< "  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ALICE_ROOT}/lib\n"
	<< "fi\n"
	<< "if test x$ALICE_PHYSICS != x ; then\n"
	<< "  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ALICE_PHYSICS}/lib\n"
	<< "fi\n"
	<< "# printenv | sort -u\n"
	<< "echo BUILD.sh@`hostname`: Building " << base << "\n"
	<< "root.exe -l -out -q PROOF-INF/BUILD.C 2>&1 | tee " 
	<< base << ".log\n"
	<< "echo BUILD.sh@`hostname`: done: $?\n"
	<< std::endl;
    out.close();
    if (gSystem->Chmod(Form("%s/PROOF-INF/BUILD.sh", dir.Data()), 0755) != 0) {
      Error("ParUtilities::MakeScriptBuildScript", 
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
  static Bool_t MakeScriptBuildMacro(const TString& dir, 
				     const TString& base, 
				     const TString& ext,
				     const TCollection* deps)  {
    std::ofstream out(Form("%s/PROOF-INF/BUILD.C", dir.Data()));
    if (!out) {
      Error("ParUtilities::MakeScriptBuildMacro","Failed to open build script");
      return false;
    }
    out << "void BUILD() {\n"
	<< "  gSystem->AddIncludePath(\"-DBUILD_PAR=1\");\n"
	<< "  gROOT->LoadMacro(\"PROOF-INF/UTIL.C\");\n"
	<< "  LoadROOTLibs();\n"
	<< "  AddAliROOT();\n"
	<< "  AddAliPhysics();\n";
    if (deps) {
      TIter       next(deps);
      TObject*    dep = 0;
      while ((dep = next())) {
	out << "  AddDep(\"" << dep->GetName() << "\");\t"
	    << "  LoadDep(\"" << dep->GetName() << "\");\n";
      }
    }
    out << "  // gDebug = 5;\n"
	<< "  // gSystem->ListLibraries();\n"
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
  static Bool_t MakeScriptUtilityScript(const TString& dir)
  {
    std::ofstream out(Form("%s/PROOF-INF/UTIL.C", dir.Data()));
    if (!out) {
      Error("ParUtilities::MakeScriptUtilityScript", 
	    "Failed to open utility script");
      return false;
    }
    out << "void LoadROOTLibs() {\n"
	<< "  gSystem->Load(\"libVMC\");\n"
	<< "  gSystem->Load(\"libNet\");\n"
	<< "  gSystem->Load(\"libTree\");\n"
	<< "  gSystem->Load(\"libPhysics\");\n"
	<< "  gSystem->Load(\"libMinuit\");\n"
	<< "}\n\n";
    out << "void AddAliROOT() {\n"
	<< "  TString val(gSystem->Getenv(\"ALICE_ROOT\"));\n"
	<< "  if (val.IsNull())\n"
	<< "    Warning(\"Add\",\"ALICE_ROOT not defined\");\n"
	<< "  else\n"
	<< "    gSystem->AddIncludePath(Form(\"-I%s/include\",val.Data()));\n"
	<< "}\n\n";
    out << "void AddAliPhysics() {\n"
	<< "  TString val(gSystem->Getenv(\"ALICE_PHYSICS\"));\n"
	<< "  if (val.IsNull())\n"
	<< "    Warning(\"Add\",\"ALICE_PHYSICS not defined\");\n"
	<< "  else\n"
	<< "    gSystem->AddIncludePath(Form(\"-I%s/include\",val.Data()));\n"
	<< "}\n\n";
    out << "void AddDep(const char* env) {\n"
	<< "  TString val(gSystem->Getenv(Form(\"%s_INCLUDE\",env)));\n"
	<< "  if (val.IsNull())\n"
	<< "    Warning(\"Add\",\"%s_INCLUDE not defined\", env);\n"
	<< "  else {\n"
	<< "    // gSystem->AddIncludePath(Form(\"-I../%s\",val.Data()));\n"
	<< "    // Prepend to include path\n"
	<< "    TString incPath(gSystem->GetIncludePath());\n"
	<< "    incPath.Prepend(Form(\"-I../%s \",val.Data()));\n"
	<< "    gSystem->SetIncludePath(incPath);\n"
	<< "    // Printf(\"Include path: %s\",incPath.Data());\n"
	<< "  }\n"
	<< "}\n\n";
    out << "void LoadDep(const char* name) {\n"
	<< "  // Printf(\"Loading dependency \\\"%s\\\" ...\",name);\n"
	<< "  // gSystem->AddDynamicPath(Form(\"../%s\",name));\n"
	<< "  // Prepend to dynamic path\n"
	<< "  TString dynPath(gSystem->GetDynamicPath());\n"
	<< "  dynPath.Prepend(Form(\"../%s:\",name));\n"
	<< "  gSystem->SetDynamicPath(dynPath);\n"
	<< "  // Printf(\"Dynamic path: %s\",dynPath.Data());\n"
	<< "  char* full = gSystem->DynamicPathName(name,true);\n"
	<< "  if (!full) \n"
	<< "   full = gSystem->DynamicPathName(Form(\"lib%s\",name),true);\n"
	<< "  if (!full) \n"
	<< "   full = gSystem->DynamicPathName(Form(\"lib%s.so\",name),true);\n"
	<< "  if (!full) \n"
	<< "   full = gSystem->DynamicPathName(Form(\"%s_C.so\",name),true);\n"
	<< "  if (!full) \n"
	<< "   full = gSystem->DynamicPathName(Form(\"%s_h.so\",name),true);\n"
	<< "  if (!full) {\n"
	<< "    Warning(\"LoadDep\",\"Module %s not found\", name);\n"
	<< "    return;\n"
	<< "  }\n"
	<< "  Printf(\"Loading \\\"%s\\\" for \\\"%s\\\"\",full,name);\n"
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
  static Bool_t MakeScriptSetupMacro(const TString& dir, 
				     const TString& base, 
				     const TString& ext,
				     const TCollection* deps)
  {
    // Make our set-up script 
    std::ofstream out(Form("%s/PROOF-INF/SETUP.C", dir.Data()));
    if (!out) {
      Error("ParUtilities::MakeScriptSetupMacro", 
	    "Failed to open setup script");
      return false;
    }
    out << "void SETUP() {\n"
	<< "  gROOT->LoadMacro(\"PROOF-INF/UTIL.C\");\n"
	<< "  Printf(\"Loading \\\"" << base << "\\\" ...\");\n"
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
  //__________________________________________________________________
  /** 
   * @{ 
   * @name PAR generation from aux file list
   */
  static Bool_t MakeAuxFilePAR(const TList& files, 
			       const TString& name,
			       Bool_t verbose=false)
  {
    // --- Check input -----------------------------------------------
    if (files.GetEntries() <= 0) return true;

    // --- Create our temporary directory ----------------------------
    Bool_t  retval = true;
    TString tmpdir(gSystem->TempDirectory());
    int     ltempl = tmpdir.Length() + 1 + 5 + 6 + 1;
    char*   templ  = new char[ltempl];
    snprintf(templ, ltempl, "%s/trainXXXXXX", tmpdir.Data());
    if (!mkdtemp(templ)) {
      Error("ParUtilities::MakeAuxFilePAR", 
	    "Failed to generate temporary directory from template %s", 
	    templ);
      return false;
    }
    if (verbose) Printf("Preparing PAR file in %s", templ);
    
    try {
      // --- Make directories for package ------------------------------
      TString dir = TString::Format("%s/%s", templ, name.Data());
      // Set-up directories 
      if (gSystem->MakeDirectory(dir) < 0) 
	throw TString::Format("Could not make directory '%s'", name.Data());
      if (gSystem->MakeDirectory(Form("%s/PROOF-INF", dir.Data()))) 
	throw TString::Format("Could not make directory %s/PROOF-INF", 
			      name.Data());

      TIter next(&files);
      TObject* o = 0;
      while ((o = next())) { 
	TString fn(o->GetName());
	if (verbose) Printf("Got %s", fn.Data());
	if (fn.BeginsWith("/")) {
	  Warning("MakeAuxFilePAR", "Will not include absolute path %s",
		  fn.Data());
	  continue; // absolute path 
	}
	
	if (gSystem->AccessPathName(fn.Data())) {
	  Warning("MakeAuxFilePAR", "Cannot access %s", fn.Data());
	  continue; // non-exist
	}
	// Loop over path components and make directories as needed 
	TObjArray*  comps = fn.Tokenize("/");
	TString     cur   = dir;
	Int_t       n     = comps->GetEntriesFast();
	if (verbose) Printf("Got %d path components in %s", n-1, fn.Data());
	Int_t       lvl   = 0;
	for (Int_t i = 0; i < n-1; i++) {
	  TObjString* comp = static_cast<TObjString*>(comps->At(i));
	  TString&    c    = comp->String();
	  if (c.IsNull()) continue;
	  if (c.EqualTo(".")) continue;
	  
	  Bool_t doMake = true;
	  if (c.EqualTo("..")) { doMake = false; lvl--; }
	  
	  cur = gSystem->ConcatFileName(cur, c);
	  if (lvl < 0) {
	    Warning("MakeAuxFilePAR", "Path %s points outside archive, ignored",
		    cur.Data());
	    break;
	  }

	  if (doMake) { 
	    lvl++;
	    if (!gSystem->AccessPathName(cur)) continue;
	    if (verbose) Printf("Making directory %s", cur.Data());
	    gSystem->MakeDirectory(cur);
	  }
	} // for(i)
	if (verbose) Printf("cur=%s for %s lvl=%d", cur.Data(), fn.Data(), lvl);
	comps->Delete();
	if (lvl < 0) continue;

	TString dest = TString::Format("%s/%s", cur.Data(), 
				       gSystem->BaseName(fn.Data()));
	if (verbose) Printf("%s -> %s", fn.Data(), dest.Data());
	// Int_t ret = gSystem->CopyFile(fn, dest, true);
	Int_t ret = gSystem->Exec(Form("cp -f %s %s", fn.Data(), dest.Data()));
	switch (ret) { 
	case -1: throw TString::Format("Couldn't open %s for copy", fn.Data());
	case -2: throw TString::Format("File %s exists", dest.Data());
	case -3: throw TString::Format("Error while copying %s", fn.Data());
	}
      }

      {
	// Make our build file 
	if (verbose) Printf("Making build script");
	std::ofstream out(Form("%s/PROOF-INF/BUILD.sh", dir.Data()));
	if (!out) {
	  Error("ParUtilities::MakeAuxFilePAR", 
		"Failed to open out shell script");
	  return false;
	}
	out << "#!/bin/sh\n\n"
	    << "echo \"Nothing to be done\"\n\n"
	    << "# EOF" << std::endl;
	out.close();
	if (gSystem->Chmod(Form("%s/PROOF-INF/BUILD.sh", dir.Data()), 0755)) {
	  Error("ParUtilities::MakeAuxFilePAR", 
		"Failed to set exectuable flags on %s/PROOF-INF/BUILD.sh", 
		dir.Data());
	  return false;
	}
      }
      {
	if (verbose) Printf("Making setup script");
	// Make our setup file 
	std::ofstream out(Form("%s/PROOF-INF/SETUP.C", dir.Data()));
	if (!out) {
	  Error("ParUtilities::MakeAuxFilePAR", 
		"Failed to open out ROOT script");
	  return false;
	}
	// The SETUP script is executed in the package's directory in
	// the package cache - not in the session directory.  Hence,
	// we take special care to get a link to the session directory
	// from the package cache directory
	out << "void SETUP()\n"
	    << "{\n"
	    << "  if (!gProofServ) {\n"
	    << "    Info(\"SETUP\",\"Not in Proof Session, nothing to do\");\n"
	    << "    return;\n"
	    << "  }\n"
	    << "  TString oldDir(gSystem->WorkingDirectory());\n"
	    << "  TSystemDirectory* dir = new TSystemDirectory(\".\",\".\");\n"
	    << "  TList*  files = dir->GetListOfFiles();\n"
	    << "  if (!gSystem->ChangeDirectory(oldDir)) {\n"
	    << "    Error(\"SETUP\", \"Failed to go back to %s\",\n"
	    << "          oldDir.Data());\n"
	    << "    return;\n"
	    << "  }\n"
	    << "  if (!files) {\n"
	    << "    Warning(\"SETUP\", \"No files\");\n"
	    << "    gSystem->Exec(\"pwd; ls -al\");\n"
	    << "    return;\n"
	    << "  }\n"
	    << "  files->Sort();\n"
	    << "  TString pkgDir = gSystem->WorkingDirectory();\n"
	    << "  TString sesDir = gProofServ->GetSessionDir();\n"
	    << "  Info(\"\",\"Session dir: %s\",sesDir);\n"
	    << "  TIter next(files);\n"
	    << "  TSystemFile* file = 0;\n"
	    << "  while ((file = static_cast<TSystemFile*>(next()))) {\n"
	    << "    TString name(file->GetName());\n"
	    << "    if (name == \".\" || name == \"..\") continue;\n"
	    << "    TString title(file->GetTitle());\n"
	    << "    TString full(gSystem->ConcatFileName(pkgDir.Data(),\n"
	    << "                                         name.Data()));\n"
	    << "    TString tgt(gSystem->ConcatFileName(sesDir.Data(),\n"
	    << "                                        name.Data()));\n"
	    << "    if(file->IsA()->InheritsFrom(TSystemDirectory::Class())){\n"
	    << "      gSystem->mkdir(tgt.Data(), true);\n"
	    << "      continue;\n"
	    << "    }\n"
	    << "    Info(\"\",\"Linking %s to %s\",full.Data(),tgt.Data());\n"
	    << "    gSystem->Symlink(full, tgt);\n"
	    << "  }\n"
	    << "}\n"
	    << "// EOF " << std::endl;
	out.close();
      }
      if (verbose) Printf("Packing up");
      Int_t ret = 0;
      ret = gSystem->Exec(Form("(cd %s && tar -c%szf %s.par %s)", 
			       templ, (verbose ? "v" : ""), 
			       name.Data(),name.Data()));
      if (ret != 0) 
	throw TString::Format("Failed to create PAR file %s.PAR from %s", 
			      name.Data(), name.Data());

      // --- Move PAR file to here -------------------------------------
      if (verbose) Printf("Move here");
      ret = gSystem->Exec(Form("mv -f %s/%s.par %s.par", templ, name.Data(), 
			       name.Data()));
      if (ret != 0) 
	throw TString::Format("Failed to rename %s/%s.par to %s.par: %s", 
			      templ, name.Data(), name.Data(), 
			      gSystem->GetError());


      if (verbose) {
	Printf("List content");
	gSystem->Exec(Form("tar tzf %s.par", name.Data()));
      }
      retval = true;
    }
    catch (TString& e) {
      Error("ParUtilities::MakeAuxFilePAR", "%s", e.Data());
      retval = false;
    }
    
    // --- Remove temporary directory --------------------------------
    gSystem->Exec(Form("rm -rf %s", templ));
    
    return retval;
  }
};
#endif
// 
// EOF
//
