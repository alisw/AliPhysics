/**
 * @file   ProcessFast.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Fri Mar 20 12:13:28 2015
 * 
 * @brief This script defines classes for looping over the data
 * produced by FastSim.C
 * 
 */
Bool_t
ProcessFast(const char* url,
	    const char* out,
	    const char* opt="g",
	    const char* extra=0)
{
  TString mkLib = gSystem->GetMakeSharedLib();
  mkLib.ReplaceAll("-std=c++14", "-std=c++98");
  gSystem->SetMakeSharedLib(mkLib);

  TString fwd = ""; // gSystem->Getenv("ANA_SRC");
  if (fwd.IsNull()) 
    fwd = gSystem->ExpandPathName("${ALICE_PHYSICS}/PWGLF/FORWARD/analysis2");
  gSystem->AddIncludePath(Form("-I${ALICE_ROOT}/include "
			       "-I${ALICE_PHYSICS}/include "
			       "-I%s/sim",
			       fwd.Data()));
  gROOT->SetMacroPath(Form("%s:%s/sim", gROOT->GetMacroPath(), fwd.Data()));
  gROOT->LoadMacro(Form("FastMonitor.C+%s",opt));
  gROOT->LoadMacro("FastShortHeader.C");
  gROOT->LoadMacro(Form("FastAnalysis.C+%s",opt));
  // gDebug = 3; // Show compile steps
  // gDebug = 7; // Keep generated files 
  gROOT->LoadMacro(Form("FastCentHelper.C+%s",opt));
  gROOT->LoadMacro(Form("dNdetaAnalysis.C+%s",opt));
  gROOT->LoadMacro(Form("dNdyAnalysis.C+%s",opt));
  gROOT->LoadMacro(Form("MidNchAnalysis.C+%s",opt));
  gDebug = 0;
  gROOT->LoadMacro(Form("spectraAnalysis.C+%s",opt));
  if (extra && extra[0] != '\0') gROOT->LoadMacro(Form("%s+%s",extra, opt));
    
  // new dNdetaMaker;
  // new spectraMaker;

  return FastAnalysis::Run(url, out, opt);
}

//
//  EOF
//  
	  
