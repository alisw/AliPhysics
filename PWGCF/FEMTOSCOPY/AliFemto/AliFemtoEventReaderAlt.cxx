
#include "AliFemtoEventReaderAlt.h"

#include "AliFemtoTrack.h"
#include "AliAODTrack.h"

#include <TRandom3.h>


AliFemtoEventReaderAlt::AliFemtoEventReaderAlt()
  : AliFemtoEventReaderAODMultSelection()
  , fRng(nullptr)
  , fEnhanceSmearing(0.0)
{
}

AliFemtoEventReaderAlt::AliFemtoEventReaderAlt(const AliFemtoEventReaderAlt &orig)
  : AliFemtoEventReaderAODMultSelection(orig)
  , fRng(nullptr)
  , fEnhanceSmearing(0.0)
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
}

void
AliFemtoEventReaderAlt::SetEnhanceSmearing(double n)
{
  if (fRng == nullptr) {
    fRng = new TRandom3();
  }
  fEnhanceSmearing = n;
}

AliFemtoTrack*
AliFemtoEventReaderAlt::CopyAODtoFemtoTrack(AliAODTrack *aod_trk)
{
  auto *femto_trk = AliFemtoEventReaderAODMultSelection::CopyAODtoFemtoTrack(aod_trk);

  if (femto_trk) {
    femto_trk->SetITSchi2(aod_trk->GetITSchi2());
    femto_trk->SetITSncls(aod_trk->GetITSNcls());
    femto_trk->SetTPCchi2(aod_trk->GetTPCchi2());
    femto_trk->SetTPCncls(aod_trk->GetTPCNcls());
    femto_trk->SetTPCnclsF(aod_trk->GetTPCNclsF());
  }

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

