#include "TString.h"
#include "TMath.h"
#include "AliForwardSettings.h"
#include "TFile.h"

//________________________________________________________________________
AliForwardSettings::AliForwardSettings() :
  fPhiAcceptanceLowEdge(0),
  fPhiAcceptanceUpEdge(2*TMath::Pi()),
  fNPhiBins(20),
  fEtaLowEdge(-4.0),
  fEtaUpEdge(6.0),
  fZVtxAcceptanceLowEdge(-10),
  fZVtxAcceptanceUpEdge(10),
  fNZvtxBins(10),
  fnoSamples(10),
  fNRefEtaBins(1),
  fNDiffEtaBins(50),
  fCentBins(60),
  nuacentral(),
  nuaforward(),
  seccorr_fwd(),
  seccorr_cen(),
  doNUA(false),
  gap(0.0),
  minpt(0.2),
  maxpt(5),
  mc(kFALSE),
  esd(kFALSE),
  tracktype(kHybrid),
  nua_mode(kNormal),
  ref_mode(kTPCref),
  useTPC{kTRUE},
  useSPD(kFALSE),
  useITS(kFALSE),
  use_primaries_cen(kFALSE),
  use_primaries_fwd(kFALSE),
  useEventcuts(kTRUE),
  centrality_estimator("V0M"),//CL0, V0M
  etagap(kTRUE),
  makeFakeHoles(kFALSE),
  fnoClusters(70),
  fCutChargedDCAxyMax(0.),
  fCutChargedDCAzMax(0.),
  doPt(kFALSE),
  stdQC(kFALSE),
  sec_corr(kFALSE),
  a5(kFALSE)
{
}
