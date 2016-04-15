#include "AliAnalysisC2Settings.h"

ClassImp(AliAnalysisC2Settings)
//________________________________________________________________________
AliAnalysisC2Settings::AliAnalysisC2Settings() :
  fDataType(AliAnalysisC2Settings::kRECON),
  fEtaAcceptanceLowEdge(-0.8),
  fEtaAcceptanceUpEdge(0.8),
  fNEtaBins(15),
  fPhiAcceptanceLowEdge(0),
  fPhiAcceptanceUpEdge(2*TMath::Pi()),
  fNPhiBins(26),
  fPtBinEdges(0),
  fMultBinEdges(0),
  fZVtxAcceptanceLowEdge(-10),
  fZVtxAcceptanceUpEdge(10),
  fNZvtxBins(20),
  fMaxDcaLong(3.0),
  fMaxDcaTang(2.4)
{
  Double_t _ptbins[] = {3.0, 4.0, 6.0, 8.0, 15.0};
  fPtBinEdges = edgeContainer(_ptbins, _ptbins + sizeof(_ptbins) / sizeof(_ptbins[0]));
  Double_t _multbins[] = {0, 20, 40, 90};
  fMultBinEdges = edgeContainer(_multbins, _multbins + sizeof(_multbins) / sizeof(_multbins[0]));
}
