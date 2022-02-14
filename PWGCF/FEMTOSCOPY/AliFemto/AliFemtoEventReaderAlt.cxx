///
/// \file AliFemtoEventReaderAlt.cxx
///


#include "AliFemtoEventReaderAlt.h"

#include "AliFemtoEvent.h"
#include "AliFemtoTrack.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliFemtoModelHiddenInfo.h"

#include <TRandom3.h>


AliFemtoEventReaderAlt::AliFemtoEventReaderAlt()
  : AliFemtoEventReaderAODMultSelection()
  , fRng(nullptr)
  , fEnhanceSmearing(0.0)
  , fDistributeMCParticles(false)
{
}

AliFemtoEventReaderAlt::AliFemtoEventReaderAlt(const AliFemtoEventReaderAlt &orig)
  : AliFemtoEventReaderAODMultSelection(orig)
  , fRng(nullptr)
  , fEnhanceSmearing(0.0)
  , fDistributeMCParticles(false)
{
  SetEnhanceSmearing(orig.GetEnhanceSmearing());
}

AliFemtoEventReaderAlt AliFemtoEventReaderAlt::operator=(const AliFemtoEventReaderAlt& rhs)
{
  if (&rhs != this) {
    AliFemtoEventReaderAODMultSelection::operator=(rhs);
    SetEnhanceSmearing(rhs.GetEnhanceSmearing());
  }

  return *this;
}

AliFemtoEventReaderAlt::~AliFemtoEventReaderAlt()
{
  delete fRng;
}

void
AliFemtoEventReaderAlt::SetEnhanceSmearing(double n)
{
  if (fRng == nullptr) {
    fRng = new TRandom3();
  }
  fEnhanceSmearing = n;
}

void
AliFemtoEventReaderAlt::SetShouldDistribute(bool val)
{
  if (fRng == nullptr) {
    fRng = new TRandom3();
  }
  fDistributeMCParticles = val;
}

AliFemtoEvent*
AliFemtoEventReaderAlt::CopyAODtoFemtoEvent()
{
  auto *femto_event = AliFemtoEventReaderAODMultSelection::CopyAODtoFemtoEvent();
  if (!femto_event) {
    return nullptr;
  }

  if (fDistributeMCParticles) {
    RandomlyDistributeParticles(*femto_event);
  }

  return femto_event;
}


AliFemtoTrack*
AliFemtoEventReaderAlt::CopyAODtoFemtoTrack(AliAODTrack *aod_trk)
{
  auto *femto_trk = AliFemtoEventReaderAODMultSelection::CopyAODtoFemtoTrack(aod_trk);

  if (!femto_trk) {
    return nullptr;
  }

  femto_trk->SetITSchi2(aod_trk->GetITSchi2());
  femto_trk->SetITSncls(aod_trk->GetITSNcls());
  femto_trk->SetTPCchi2(aod_trk->GetTPCchi2());
  femto_trk->SetTPCncls(aod_trk->GetTPCNcls());
  femto_trk->SetTPCnclsF(aod_trk->GetTPCNclsF());

  if (fEnhanceSmearing != 0.0) {
    auto p = femto_trk->P();
    p.SetX(p.x() * fRng->Gaus(1, fEnhanceSmearing));
    p.SetY(p.y() * fRng->Gaus(1, fEnhanceSmearing));
    p.SetZ(p.z() * fRng->Gaus(1, fEnhanceSmearing));
    femto_trk->SetP(p);
  }

  return femto_trk;
}


void
AliFemtoEventReaderAlt::CopyPIDtoFemtoTrack(AliAODTrack *aod_trk, AliFemtoTrack *femto_trk)
{
  AliFemtoEventReaderAODMultSelection::CopyPIDtoFemtoTrack(aod_trk, femto_trk);

  femto_trk->SetITSchi2(aod_trk->GetITSchi2());
  femto_trk->SetITSncls(aod_trk->GetITSNcls());
  femto_trk->SetTPCchi2(aod_trk->GetTPCchi2());
  femto_trk->SetTPCncls(aod_trk->GetTPCNcls());
  femto_trk->SetTPCnclsF(aod_trk->GetTPCNclsF());
}

void
AliFemtoEventReaderAlt::RandomlyDistributeParticles(AliFemtoEvent &femto_event)
{
  const double R = 7.406240171018426;

  if (auto *mc_header = dynamic_cast<AliAODMCHeader*>(fEvent->FindListObject(AliAODMCHeader::StdBranchName()))) {
    double v[3];
    mc_header->GetVertex(v);

    double pv[3];
    fEvent->GetPrimaryVertex()->GetXYZ(pv);

    const double
      psi = mc_header->GetReactionPlaneAngle(),
      impact_parameter = mc_header->GetImpactParameter();

    const double
      wx = 2.0 * R - impact_parameter,
      wy = std::sqrt(std::max(R*R -  impact_parameter * impact_parameter / 4.0, 1e-5)),
      wz = 1.0;

    for (auto *track : *femto_event.TrackCollection()) {
      auto *hi = track->GetHiddenInfo();

      if (auto *info = dynamic_cast<AliFemtoModelHiddenInfo*>(hi)) {

        if (info->GetOrigin() != 0) {
          continue;
        }

        AliFemtoThreeVector x(fRng->Gaus(0, wx),
                              fRng->Gaus(0, wy),
                              fRng->Gaus(0, wz));
        x.RotateZ(psi);

        AliFemtoLorentzVector emission_point(x, 0);
        emission_point += *info->GetEmissionPoint();

        info->SetEmissionPoint(emission_point);
      }
    }
  }
}
