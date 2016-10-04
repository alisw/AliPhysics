/**
 * @file   MakeDeltas.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Wed Aug 31 10:57:32 2016
 * 
 * @brief  
 * 
 * 
 * @ingroup pwglf_forward_tracklets
 */

/** 
 * 
 * 
 * @param realFileName 
 * @param simFileName 
 *
 * @ingroup pwglf_forward_tracklets
 *
 * @relates MakeDeltaWeights
 */
void
MakeDeltas(const char* realFileName,
	   const char* simFileName,
	   Int_t       dimen=2,
	   Bool_t      scaleToTail=true)
{
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include "
			  "-I${ALICE_PHYSICS}/include "
			  "-I${ANA_SRC}/dndeta/tracklets3");
  gROOT->SetMacroPath(Form("${ANA_SRC}/dndeta/tracklets3:%s",
			   gROOT->GetMacroPath()));
  gROOT->LoadMacro("AliAODTracklet.C+g");
  gROOT->LoadMacro("AliTrackletWeights.C+g");
  gROOT->LoadMacro("AliTrackletAODUtils.C+g");
  gROOT->LoadMacro("MakeDeltaWeights.C+g");

  TCanvas* c1 = new TCanvas("c1","c1");
  MakeDeltaWeights* mdw = new MakeDeltaWeights;
  AliTrackletBaseWeights* w = mdw->Run(realFileName,simFileName,
				       dimen, scaleToTail);
  TCanvas* c2 = new TCanvas("c2","c2");
  w->Draw();
}
//
// EOF
// 
