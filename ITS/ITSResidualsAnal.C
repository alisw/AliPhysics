#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TArrayI.h>
#include <TSelector.h>
#include <TString.h>
#include <TTree.h>
#include "AliGeomManager.h"
#include "AliTrackPoints.h"
#include "AliITSResidualsAnalysis.h"
#endif
void ITSResidualsAnal(Int_t layer = 1,
		      TString AliTrackPoints = "AliTrackPoints.root",
		      TString geometry = "geometry.root")
{
  //
  // Sample Macro for AliITSResidualsAnalysis Class
  //
  // Parameters:
  //   layer          = layer to analyze
  //   AliTrackPoints = file with the AliTrackPoints (including path if needed)
  //   geometry       = file with the geometry (including path if needed)
  //
  // Results are stored in ResidualsAnalysisTree.root in a tree called "analysisTree"
  //

  AliITSResidualsAnalysis *res = new AliITSResidualsAnalysis(AliTrackPoints,geometry);

  TArrayI *volids = (TArrayI*)res->GetSingleLayerVolids(layer);

  res->InitHistograms(volids);

  AliGeomManager::ELayerID lay = (AliGeomManager::ELayerID)layer;
  res->ProcessPoints("fast",0,lay,lay,"");

  return;
}
