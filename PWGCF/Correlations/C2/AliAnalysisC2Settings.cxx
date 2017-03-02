#include "TString.h"
#include "TMath.h"
#include "AliAnalysisC2Settings.h"

ClassImp(AliAnalysisC2Settings)
//________________________________________________________________________
AliAnalysisC2Settings::AliAnalysisC2Settings() :
  fDataType(-1),
  fEtaEdgesBwd(0),
  fEtaEdgesIts(0),
  fEtaEdgesFwd(0),
  fPhiAcceptanceLowEdge(0),
  fPhiAcceptanceUpEdge(2*TMath::Pi()),
  fNPhiBins(26),
  fPtBinEdges(0),
  fMultBinEdges(0),
  fZVtxAcceptanceLowEdge(-10),
  fZVtxAcceptanceUpEdge(10),
  fNZvtxBins(20),
  fMaxDcaLong(3.0),
  fMaxDcaTang(2.4),
  fMultEstimatorRefMult08("RefMult08"),
  fMultEstimatorV0M("V0M"),
  fMultEstimatorValidTracks("ValidTracks"),
  fMultEstimator(""),
  fTriggerCint7("CINT7-B-"),
  fTriggerVhmV0M("CVHMV0M-B-"),
  fTriggerMask(0),
  fTrigger_str(""),
  fUseFMD(false)
{
  Double_t _ptbins[] = {3.0, 4.0, 6.0, 8.0, 15.0};
  fPtBinEdges = edgeContainer(_ptbins, _ptbins + sizeof(_ptbins) / sizeof(_ptbins[0]));
  Double_t _multbins[] = {0, 20, 40, 90};
  fMultBinEdges = edgeContainer(_multbins, _multbins + sizeof(_multbins) / sizeof(_multbins[0]));

  const Float_t eta_bin_width = 0.10;
  for (float etaEdge= -3.1; etaEdge <= -1.7 + 0.0001; etaEdge += eta_bin_width) {
    this->fEtaEdgesBwd.push_back(etaEdge);
  }
  for (float etaEdge= -1.7; etaEdge <= 1.7 + 0.0001; etaEdge += eta_bin_width) {
    this->fEtaEdgesIts.push_back(etaEdge);
  }
  for (float etaEdge= 1.7; etaEdge <= 4.9 + 0.0001; etaEdge += eta_bin_width) {
    this->fEtaEdgesFwd.push_back(etaEdge);
  }

}
// this comment serves no purpose
