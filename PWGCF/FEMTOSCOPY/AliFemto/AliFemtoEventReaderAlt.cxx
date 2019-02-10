
#include "AliFemtoEventReaderAlt.h"

#include "AliFemtoTrack.h"
#include "AliAODTrack.h"

AliFemtoEventReaderAlt::AliFemtoEventReaderAlt()
  : AliFemtoEventReaderAODMultSelection()
{
}

AliFemtoEventReaderAlt::~AliFemtoEventReaderAlt()
{
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

  return femto_trk;
}
