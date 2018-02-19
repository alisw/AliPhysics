#include "TString.h"
#include "TMath.h"
#include "AliForwardFlowRun2Settings.h"
#include "TFile.h"

//________________________________________________________________________
AliForwardFlowRun2Settings::AliForwardFlowRun2Settings() :
  fDataType(-1),
  fPhiAcceptanceLowEdge(0),
  fPhiAcceptanceUpEdge(2*TMath::Pi()),
  fNPhiBins(20),
  fZVtxAcceptanceLowEdge(-10),
  fZVtxAcceptanceUpEdge(10),
  fNZvtxBins(20),
  qctype("std"),
  fnoSamples(10),
  fNRefEtaBins(1),
  fNDiffEtaBins(24),
  fCentBins(10),
  nuacentral(),
  nuaforward(),
  doNUA(false),
  gap(0.0),
  mc(false)
{


}
