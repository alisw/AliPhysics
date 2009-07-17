//
// Class AliRsnCutMgr
//
// The cut manager: contains a complete set of cut definitions
// to be applied to all possible targets (one for each target),
// in order to ease the set-up procedure of cuts and allow to
// pass them at once to each object which must use them
//
// author: Martin Vala (martin.vala@cern.ch)
//

#include "AliLog.h"

#include "AliRsnCutSet.h"
#include "AliRsnCutMgr.h"

ClassImp(AliRsnCutMgr)

//_____________________________________________________________________________
AliRsnCutMgr::AliRsnCutMgr() :
    TNamed("defaultName", "defaultTitle")
{
//
// Constructor without arguments.
//

  Int_t i;
  for (i = 0; i < AliRsnCut::kLastCutTarget; i++) {
    fCutSets[i] = 0;
  }
}

//_____________________________________________________________________________
AliRsnCutMgr::AliRsnCutMgr(const char *name, const char *title) :
    TNamed(name, title)
{
//
// Constructor with name and title.
//

  Int_t i;
  for (i = 0; i < AliRsnCut::kLastCutTarget; i++) {
    fCutSets[i] = 0;
  }
}

//_____________________________________________________________________________
AliRsnCutMgr::~AliRsnCutMgr()
{
//
// Destructor.
// Deletes all cut definitions.
//

  Int_t i;
  for (i = 0; i < AliRsnCut::kLastCutTarget; i++) {
    delete fCutSets[i];
  }
}

//_____________________________________________________________________________
void AliRsnCutMgr::SetCutSet(AliRsnCut::ETarget type, AliRsnCutSet* const cutset)
{
//
// Assign a cut set to a given target
//

  if (!fCutSets[type]) fCutSets[type] = (AliRsnCutSet*) cutset->Clone();
  AliDebug(AliLog::kDebug, Form("DatasetName %s", fCutSets[type]->GetName()));
}

//_____________________________________________________________________________
Bool_t AliRsnCutMgr::IsSelected(AliRsnCut::ETarget type, TObject*const obj)
{
//
// Check if a given object passes the cuts defined for it.
// The target of the check is here a TObject, in order to allo generality
// but then the kind of cut to be used is defined as first argument, and
// in the single cut it will be checked if it is appropriate for passed target
//

  AliDebug(AliLog::kDebug, "<-");
  if (!fCutSets[type]) return kTRUE;

  switch (type) {
  case AliRsnCut::kParticle:
    return fCutSets[type]->IsSelected(type, (AliRsnDaughter*)obj);
    break;
  case AliRsnCut::kPair:
    return fCutSets[type]->IsSelected(type, (AliRsnPairParticle*)obj);
    break;
  case AliRsnCut::kEvent:
    return fCutSets[type]->IsSelected(type, (AliRsnEvent*)obj);
    break;
  default:
    AliWarning("Wrong target selected.");
    return kTRUE;
    break;
  }

  return kTRUE;
}
