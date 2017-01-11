/**
 * @file   CreateIndex.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Mon Nov 14 20:38:50 2016
 * 
 * @brief  Create an index of downloaded files 
 *
 * @ingroup pwglf_forward_trains_helper
 */

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


void
CreateIndex(const TString& dir,
	    const TString& tree="esdTree",
	    const char*    remote=0,
	    const char*    plibs=0,
	    const char*    ppars=0,
	    const char*    psrcs=0)
{
  gROOT->SetMacroPath(Form("$ALICE_PHYSICS/PWGLF/FORWARD/trains:%s", 
			   gROOT->GetMacroPath()));
  gROOT->LoadMacro("ChainBuilder.C+");
  gROOT->Macro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/scripts/LoadLibs.C");

  gSystem->AddDynamicPath("/usr/lib");

  // Load basic ROOT libraries
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
  TString libs(plibs);
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
  Bool_t local           = (!remote || remote[0] == '\0');
  Bool_t localLibsNotPar = false;
  
  // Load packages
  TString pars(ppars);
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
  TString srcs(psrcs);
  if (!srcs.IsNull()) {
    TObjArray*  srcArray = srcs.Tokenize(" ");
    TObjString* src      = 0;
    TIter       nextSrc(srcArray);
    while ((src = static_cast<TObjString*>(nextSrc()))) { 
      const TString& srcName = src->String();
      gROOT->ProcessLine(Form(".L %s+g", srcName.Data()));
    }
  }
  
  Bool_t mc = false;
  Bool_t zip = false;
  if (tree.BeginsWith("mc")) { 
    mc = true;
    zip = true;
    tree.Remove(0,2);
  }
  if (tree.BeginsWith("zip")) { 
    zip = true;
    tree.Remove(0,3);
  }
  
  TString pat("*.root");
  if      (tree.EqualTo("esdTree",  TString::kIgnoreCase)) pat="AliESDs*";
  else if (tree.EqualTo("aodTree",  TString::kIgnoreCase)) pat="AliAOD*";
  else    Warning("", "Unknown tree: %s, pattern set to *.root", tree.Data());
  if (zip) {
    pat.Prepend("root_archive.zip@");
    pat.ReplaceAll("*", ".root");
    pat.ReplaceAll(".zip", "*.zip");
  }


  TString opts;
  opts.Append(Form("pattern=%s", pat.Data()));
  opts.Append("&check");
  opts.Append("&clean");
  opts.Append("&recursive");
  opts.Append("&verbose");
  if (mc) opts.Append("&mc");

  TString realDir(dir);
  if (!remote)              realDir = gSystem->ExpandPathName(dir.Data());
  if (realDir.EqualTo(".")) realDir = gSystem->WorkingDirectory();

  TUrl url;
  url.SetProtocol("local");
  url.SetPort(0);
  url.SetFile(realDir);
  url.SetAnchor(tree);
  url.SetOptions(opts);
  
  Printf("Running ChainBuilder::CreateCollection(\"%s/index.root\",\"%s\")",
	 realDir.Data(), url.GetUrl());
  TString out(Form("%s/%s.root", realDir.Data(),
		   !remote ? "index" : "remote"));
  ChainBuilder::CreateCollection(out, url, remote);
}

				 
  
  
