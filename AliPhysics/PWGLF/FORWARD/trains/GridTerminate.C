/**
 * @file   GridTerminate.C
 * @author Christian Holm Christensen <cholm@master.hehi.nbi.dk>
 * @date   Wed Jan 23 21:22:28 2013
 * 
 * @brief  Worker script to run terminate part for Grid  
 * 
 * 
 * @ingroup pwglf_forward_trains_helper
 */
#ifndef __CINT__
# include <TString.h>
# include <TSystem.h>
# include <TGrid.h>
# include <TFile.h>
# include <TObjArray.h>
# include <TObjString.h>
# include <TError.h>
# include <TEnv.h>
# include <TROOT.h>
# include <AliAnalysisManager.h>
# include <AliAnalysisAlien.h>
# include <fstream>
#else
class TString;
#endif
/** 
 * Load a library 
 * 
 * @param libName Library name 
 * 
 * @return true on success 
 */
Bool_t LoadLib(const char* libName)
{
  if (gSystem->Load(libName) < 0) {
    Error("GridTerminate", "Failed to load library %s",libName);
    return false;
  }
  Info("GridTerminate","Loaded library %s",libName);
  return true;
}
/** 
 * Load a PAR 
 * 
 * @param parName PAR file name 
 * 
 * @return true on success
 */
Bool_t LoadPar(const char* parName)
{
  if (!AliAnalysisAlien::SetupPar(parName)) {
    Error("GridTerminate","Failed to load PAR %s",parName);
    return false;
  }
  Info("GridTerminate","Loaded package %s",parName);
  return true;
}

/** 
 * Load the analysis handler if needed 
 * 
 * @param name 
 * 
 * @return 
 */
Bool_t LoadHandler(const TString& name)
{
  // Load plugin
  TFile* file = TFile::Open(Form("%s_plugin.root",name.Data()),"READ");
  // TFile* plug = TFile::Open("plugin.root","READ");
  if (!file) {
    // Error("GridTerminate","Failed to open %s_plugin.root",name.Data());
    Error("GridTerminate","Failed to open %s_plugin.root",
	  name.Data());
    return false;
  }
  AliAnalysisAlien* handler = 
    static_cast<AliAnalysisAlien*>(file->Get("plugin"));
  if (!handler) {
    Error("GridTerminate","Failed to load plugin");
    return false;
  }
  Info("GridTerminate","Setting grid handler");
  handler->SetRunMode("terminate");
  handler->SetMergeViaJDL(true);
  AliAnalysisManager::GetAnalysisManager()->SetGridHandler(handler);

  return true;
}

/** 
 * Setup our manager et al 
 * 
 * @param name            Name of the job
 * @param libs            Libraries to load (space separated string)
 * @param pars            PARs to load (space separated string)
 * @param srcs            Sources to load (space separated string)
 * @param local           If true, run local terminate job on already 
 *                        merged and downloaded files.
 * @param localLibsNotPar if @a local is true and this is true, then
 *                        we load the corresponding locally compiled 
 *                        library for each specified PAR file.
 * 
 * @return true on success 
 */
Bool_t Setup(const TString& name, 
	     const TString& libs, 
	     const TString& pars, 
	     const TString& srcs,
	     Bool_t         local=false,
	     Bool_t         localLibsNotPar=true)
{
  // Load basic ROOT libraries
  gSystem->AddDynamicPath("/usr/lib");
  if (gSystem->Load("libTree")       < 0) return false;
  if (gSystem->Load("libGeom")       < 0) return false;
  if (gSystem->Load("libVMC")        < 0) return false;
  if (gSystem->Load("libPhysics")    < 0) return false;
  if (gSystem->Load("libMinuit")     < 0) return false;

  // Load basic AliROOT libraries
  if (gSystem->Load("libSTEERBase")     < 0) return false;
  if (gSystem->Load("libESD")           < 0) return false;
  if (gSystem->Load("libAOD")           < 0) return false;
  if (gSystem->Load("libANALYSIS")      < 0) return false;
  if (gSystem->Load("libOADB")          < 0) return false;
  if (gSystem->Load("libANALYSISalice") < 0) return false;

  // Load libraries
  if (!libs.IsNull()) {
    TObjArray*  libsArray = libs.Tokenize(" ");
    TObjString* lib       = 0;
    TIter       nextLib(libsArray);
    while ((lib = static_cast<TObjString*>(nextLib()))) {
      const TString& libName = lib->String();
      if (libName.Contains("libSTEERBase") ||
	  libName.Contains("libESD")       ||
	  libName.Contains("libAOD")       ||
	  libName.Contains("libANALYSIS")  ||
	  libName.Contains("libOADB")      ||
	  libName.Contains("libANALYSISalice")) continue;
      if (!libName.Contains(".so")) continue;
      if (!LoadLib(libName.Data())) return false;
    }
    libsArray->Delete();
  }
  
  // Load packages
  if (!pars.IsNull()) {
    TObjArray*  parArray = pars.Tokenize(" ");
    TObjString* par      = 0;
    TIter       nextPar(parArray);
    while ((par = static_cast<TObjString*>(nextPar()))) { 
      TString parName(par->String());
      if (parName.EndsWith(".par")) parName.ReplaceAll(".par", "");
      if (parName.EqualTo("STEERBase") ||
	  parName.EqualTo("ESD")       ||
	  parName.EqualTo("AOD")       ||
	  parName.EqualTo("ANALYSIS")  ||
	  parName.EqualTo("OADB")      ||
	  parName.EqualTo("ANALYSISalice")) continue;
      Bool_t ret = true;
      if (local && localLibsNotPar) 
	ret = LoadLib(Form("lib%s.so", parName.Data()));
      else 
	ret = LoadPar(parName.Data());
      if (!ret) return false;
    }
  }

  // Load sources
  if (!srcs.IsNull()) {
    TObjArray*  srcArray = srcs.Tokenize(" ");
    TObjString* src      = 0;
    TIter       nextSrc(srcArray);
    while ((src = static_cast<TObjString*>(nextSrc()))) { 
      const TString& srcName = src->String();
      gROOT->ProcessLine(Form(".L %s+g", srcName.Data()));
    }
  }

  // Connect to the grid
  gEnv->SetValue("XSec.GSI.DelegProxy", "2");
  // TGrid::Connect("alien://");
  // if (!gGrid) {
  //   Error("GridTerminate", "Failed to connect to AliEn");
  //   return false;
  // }

  // Load the analysis manager from file
  TString base(name);
  base.Append(".root");
  if (gSystem->AccessPathName(base.Data())) {
    // Couldn't read from current directory, try sub-dir
    TString sub(gSystem->ConcatFileName(name, base));
    if (gSystem->AccessPathName(sub)) {
      Error("GridTerminate","Couldn't find manager file %s",base.Data());
      return false;
    }
    base = sub;
  }
  AliAnalysisManager* mgr= AliAnalysisAlien::LoadAnalysisManager(base);
  if (!mgr) {
    Error("GridTerminate", "Failed to load manager from %s",base.Data());
    return false;
  }
  if (!name.EqualTo(mgr->GetName())) {
    Error("GridTerminate","Read manager %s is not %s",
	  mgr->GetName(),name.Data());
    return false;
  }
  Info("GridTerminate","Loaded analysis manager");

  // If we do a local merge, do not do any else 
  if (local) return true;

  return LoadHandler(name);
}

/** 
 * Submit the terminate job 
 * 
 * @param name            Name of the job
 * @param libs            Libraries to load (space separated string)
 * @param pars            PARs to load (space separated string)
 * @param srcs            Sources to load (space separated string)
 * @param local           If true, run local terminate job on already 
 *                        merged and downloaded files.
 * @param localLibsNotPar if @a local is true and this is true, then
 *                        we load the corresponding locally compiled 
 *                        library for each specified PAR file.
 * 
 * @return true on success 
 */
Bool_t GridTerminate(const TString& name, 
		     const TString& libs, 
		     const TString& pars, 
		     const TString& srcs,
		     Bool_t         local=false,
		     Bool_t         localLibsNotPar=true)
{
  if (!Setup(name, libs, pars, srcs, local, localLibsNotPar)) return false; 

  // Run the terminate job
  Info("GridTerminate","Starting terminate job - %s", local ? "locally" : "on the grid");

  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!local) {

    AliAnalysisAlien* handler = static_cast<AliAnalysisAlien*>(mgr->GetGridHandler());
    if (!handler) { 
      Error("GridTerminate", "Manager does now have an AliEn handler");
      return false;
    }
    if (mgr->StartAnalysis("grid") < 0) return false;

    std::ofstream outJobs(Form("%s_merge.jobid", mgr->GetName()));
    outJobs << handler->GetGridJobIDs() << std::endl;
    outJobs.close();
    
    std::ofstream outStages(Form("%s_merge.stage", mgr->GetName()));
    outStages << handler->GetGridStages() << std::endl;
    outStages.close();

    return true;
  }
  
  // mgr->SetDebugLevel(2);
  mgr->SetSkipTerminate(false);
  TTree* dummy = 0;
  if (mgr->StartAnalysis("gridterminate", dummy, -1, 0) < 0) return false;
  
  return true;
}

// 
// EOF
//

