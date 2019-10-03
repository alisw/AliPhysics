/**
 * @file   TryBuild.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Apr 27 16:52:29 2016
 * 
 * @brief  
 * 
 * @ingroup pwglf_forward_tracklets
 * 
 */
/** 
 * Try to build the code 
 * 
 * @ingroup pwglf_forward_tracklets
 */
void
TryBuild()
{
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include "
			  "-I${ALICE_PHYSICS}/include "
			  "-I${ANA_SRC}/dndeta/tracklets3");
  gROOT->SetMacroPath(Form("${ANA_SRC}/dndeta/tracklets3:%s",
			   gROOT->GetMacroPath()));
  gROOT->LoadMacro("AliAODSimpleHeader.C+g");
  gROOT->LoadMacro("AliSimpleHeaderTask.C+g");
  gROOT->LoadMacro("AliAODTracklet.C+g");
  gROOT->LoadMacro("AliTrackletWeights.C+g");
  gROOT->LoadMacro("AliTrackletAODUtils.C+g");
  gROOT->LoadMacro("AliTrackletAODTask.C+g");
  gROOT->LoadMacro("AliTrackletAODdNdeta.C+g");
  gROOT->LoadMacro("AliTrackletdNdeta.C+g");
  gROOT->LoadMacro("AliTrackletdNdeta2.C+g");
  gROOT->LoadMacro("MakeDeltaWeights.C+g");
}
//
// EOF
//
