#include <TString.h>

#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliVEventHandler.h"

#include "AliPicoHeaderCJ.h"

ClassImp(AliPicoHeaderCJ)

//_____________________________________________________________________________
AliPicoHeaderCJ::AliPicoHeaderCJ() :
TNamed(),
fPhysSelMask(0),
fFiredTriggerClass(""),
fCentralityV0M(-1.),
fCentralityV0A(-1.),
fCentralityCL1(-1.),
fCentralityZNA(-1.),
fEventPlane(999.),
fBackgroundRhoRD02(0.),
fBackgroundRhoRD03(0.),
fBackgroundRhoRD04(0.),
fBackgroundRhoMC02(0.),
fBackgroundRhoMC03(0.),
fBackgroundRhoMC04(0.)
{
//
// AliPicoHeaderCJ::AliPicoHeaderCJ
//

  for (Int_t i=3; i--;) fVtx[i] = 0.;
}

//_____________________________________________________________________________
AliPicoHeaderCJ::AliPicoHeaderCJ(const AliPicoHeaderCJ &src) :
TNamed(src),
fPhysSelMask(src.fPhysSelMask),
fFiredTriggerClass(src.fFiredTriggerClass),
fCentralityV0M(src.fCentralityV0M),
fCentralityV0A(src.fCentralityV0A),
fCentralityCL1(src.fCentralityCL1),
fCentralityZNA(src.fCentralityZNA),
fEventPlane(src.fEventPlane),
fBackgroundRhoRD02(src.fBackgroundRhoRD02),
fBackgroundRhoRD03(src.fBackgroundRhoRD03),
fBackgroundRhoRD04(src.fBackgroundRhoRD04),
fBackgroundRhoMC02(src.fBackgroundRhoMC02),
fBackgroundRhoMC03(src.fBackgroundRhoMC03),
fBackgroundRhoMC04(src.fBackgroundRhoMC04)
{
//
// AliPicoHeaderCJ::AliPicoHeaderCJ
//

  for (Int_t i=3; i--;) fVtx[i] = src.fVtx[i];
}

//_____________________________________________________________________________
AliPicoHeaderCJ& AliPicoHeaderCJ::operator=(const AliPicoHeaderCJ &src)
{
//
// AliPicoHeaderCJ::operator=
//

  if (&src==this) return *this;

  TNamed::operator=(src);

  fPhysSelMask       = src.fPhysSelMask;
  fFiredTriggerClass = src.fFiredTriggerClass;
  fCentralityV0M     = src.fCentralityV0M;
  fCentralityV0A     = src.fCentralityV0A;
  fCentralityCL1     = src.fCentralityCL1;
  fCentralityZNA     = src.fCentralityZNA;

  fEventPlane        = src.fEventPlane;

  fBackgroundRhoRD02 = src.fBackgroundRhoRD02;
  fBackgroundRhoRD03 = src.fBackgroundRhoRD03;
  fBackgroundRhoRD04 = src.fBackgroundRhoRD04;

  fBackgroundRhoMC02 = src.fBackgroundRhoMC02;
  fBackgroundRhoMC03 = src.fBackgroundRhoMC03;
  fBackgroundRhoMC04 = src.fBackgroundRhoMC04;

  for (Int_t i=3; i--;) fVtx[i] = src.fVtx[i];

  return *this;
}

//_____________________________________________________________________________
AliPicoHeaderCJ::~AliPicoHeaderCJ()
{
//
// AliPicoHeaderCJ::~AliPicoHeaderCJ
//
}

//_____________________________________________________________________________
void AliPicoHeaderCJ::SetEventInfo(AliVEventHandler* const pEH)
{
//
// AliPicoHeaderCJ::SetEventInfo
//

  AliVEvent   *pEV = pEH->GetEvent();
  AliAODEvent *pA = dynamic_cast<AliAODEvent*>(pEV);
  AliESDEvent *pE = dynamic_cast<AliESDEvent*>(pEV);

  fPhysSelMask = pEH->IsEventSelected(); 
  if (pA) fFiredTriggerClass = pA->GetFiredTriggerClasses();
  if (pE) fFiredTriggerClass = pE->GetFiredTriggerClasses();

  const AliVVertex *pVtx = pEV->GetPrimaryVertex();
  this->SetTitle(pVtx->GetTitle());
  pVtx->GetXYZ(fVtx);

  AliCentrality *pCent = pEV->GetCentrality();
  if (pCent) {
    fCentralityV0M = pCent->GetCentralityPercentile("V0M");
    fCentralityV0A = pCent->GetCentralityPercentile("V0A");
    fCentralityCL1 = pCent->GetCentralityPercentile("CL1");
    fCentralityZNA = pCent->GetCentralityPercentile("ZNA");
  }
 
  AliEventplane   *pEventPlane = pEV->GetEventplane();
  if (pEventPlane) fEventPlane = pEventPlane->GetEventplane("Q");

  return;
}

//_____________________________________________________________________________
Double_t AliPicoHeaderCJ::BackgroundRho(const TString sJet) const
{
//
//  AliPicoHeaderCJ::BackgroundRho
//

  if (sJet.IsNull()) return 0.;
//=============================================================================

  if (sJet=="RD02") return fBackgroundRhoRD02;
  if (sJet=="RD03") return fBackgroundRhoRD03;
  if (sJet=="RD04") return fBackgroundRhoRD04;

  if (sJet=="MC02") return fBackgroundRhoMC02;
  if (sJet=="MC03") return fBackgroundRhoMC03;
  if (sJet=="MC04") return fBackgroundRhoMC04;
//=============================================================================

  return 0.;
}

//_____________________________________________________________________________
void AliPicoHeaderCJ::Reset()
{
//
//  AliPicoHeaderCJ::Reset
//

  fPhysSelMask       = 0,
  fFiredTriggerClass = "";

  fCentralityV0M = -1.;
  fCentralityV0A = -1.;
  fCentralityCL1 = -1.;
  fCentralityZNA = -1.;

  fEventPlane       = 999.;

  fBackgroundRhoRD02 = 0.;
  fBackgroundRhoRD03 = 0.;
  fBackgroundRhoRD04 = 0.;

  fBackgroundRhoMC02 = 0.;
  fBackgroundRhoMC03 = 0.;
  fBackgroundRhoMC04 = 0.;

  for (Int_t i=3; i--;) fVtx[i] = 0.;

  return;
}
