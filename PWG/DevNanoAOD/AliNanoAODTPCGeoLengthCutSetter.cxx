#include "AliNanoAODTPCGeoLengthCutSetter.h"

#include "AliAODEvent.h"
#include "AliESDtrack.h"
#include "AliNanoAODHeader.h"
#include "AliVTrack.h"

ClassImp(AliNanoAODTPCGeoLengthCutSetter)

    AliNanoAODTPCGeoLengthCutSetter::AliNanoAODTPCGeoLengthCutSetter(
        const char* name)
    : AliNanoAODCustomSetter(name),
      fMode(0),
      fDeltaY(3.0),
      fDeltaZ(220.0),
      fMagField(-1),
      fRequireCutGeoNcrNclLength(130),
      fRequireCutGeoNcrNclGeom1Pt(1.5),
      fCutGeoNcrNclFractionNcr(0.85),
      fCutGeoNcrNclFractionNcl(0.7),
      fIndex(-1) {}

AliNanoAODTPCGeoLengthCutSetter::~AliNanoAODTPCGeoLengthCutSetter() {}

void AliNanoAODTPCGeoLengthCutSetter::SetNanoAODHeader(const AliAODEvent* event,
                                                       AliNanoAODHeader* head,
                                                       TString varListHeader) {
  fMagField = event->GetMagneticField();
}

void AliNanoAODTPCGeoLengthCutSetter::SetNanoAODTrack(
    const AliAODTrack* aodTrack,
    AliNanoAODTrack* spTrack) {
  if (fIndex == -1)
    fIndex =
        AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstTPCGeoLength");
  if (fIndex == -1)
    fIndex = -2;  // prevent useless run.

  if (fIndex >= 0) {
    auto checkResult = true;
    AliESDtrack fESDTrack(aodTrack);
    fESDTrack.SetTPCClusterMap(aodTrack->GetTPCClusterMap());
    fESDTrack.SetTPCSharedMap(aodTrack->GetTPCSharedMap());
    fESDTrack.SetTPCPointsF(aodTrack->GetTPCNclsF());

    auto nCrossedRowsTPC = fESDTrack.GetTPCCrossedRows();
    auto lengthInActiveZoneTPC = fESDTrack.GetLengthInActiveZone(
        fMode, fDeltaY, fDeltaZ, fMagField);
    auto cutGeoNcrNclLength = fRequireCutGeoNcrNclLength -
                              TMath::Power(TMath::Abs(fESDTrack.GetSigned1Pt()),
                                           fRequireCutGeoNcrNclGeom1Pt);

    if (lengthInActiveZoneTPC < cutGeoNcrNclLength)
      checkResult = false;
    if (nCrossedRowsTPC < fCutGeoNcrNclFractionNcr * cutGeoNcrNclLength)
      checkResult = false;
    if (fESDTrack.GetTPCncls() < fCutGeoNcrNclFractionNcl * cutGeoNcrNclLength)
      checkResult = false;

    spTrack->SetVar(fIndex, (checkResult) ? 1. : 0.);
  }
};