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
	    const char* opt="g")
{
  TString fwd = ""; // gSystem->Getenv("ANA_SRC");
  if (fwd.IsNull()) 
    fwd = gSystem->ExpandPathName("${ALICE_PHYSICS}/PWGLF/FORWARD/analysis2");
  gSystem->AddIncludePath(Form("-I${ALICE_ROOT}/include "
			       "-I${ALICE_PHYSICS}/include "
			       "-I%s/include",
			       fwd.Data()));
  gROOT->SetMacroPath(Form("%s:%s/sim", gROOT->GetMacroPath(), fwd.Data()));
  gROOT->LoadMacro(Form("%s/sim/FastAnalysis.C+%s",fwd.Data(),opt));


  return FastAnalysis::Run(url, out, opt);
}

//
//  EOF
//  
	  
