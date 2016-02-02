///
/// \file AliFemtoEventReaderAODMultSelection.cxx
///


#include "AliFemtoEventReaderAODMultSelection.h"
#include "AliFemtoEvent.h"
#include "../../../OADB/COMMON/MULTIPLICITY/AliMultSelection.h"  // "AliMultSelection.h"

AliFemtoEvent* AliFemtoEventReaderAODMultSelection::CopyAODtoFemtoEvent()
{
  AliFemtoEvent *femto_event = AliFemtoEventReaderAODChain::CopyAODtoFemtoEvent();

  AliMultSelection *mult_selection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
  if (mult_selection == NULL) {
    cout << "W-AliFemtoEventReaderAODMultSelection: No AliMultSelection found."
            "Was the AddTaskMultSelection run? ignoring.\n";
    return femto_event;
  }

  femto_event->SetCentralityV0(mult_selection->GetMultiplicityPercentile("V0M"));
  femto_event->SetCentralityCL1(mult_selection->GetMultiplicityPercentile("CL1"));

  switch (fEstEventMult) {
  case kCentrality:
    femto_event->SetNormalizedMult(lrint(10 * mult_selection->GetMultiplicityPercentile("V0M")));
    break;
  case kCentralityCL1:
    femto_event->SetNormalizedMult(lrint(10 * mult_selection->GetMultiplicityPercentile("CL1")));
    break;
  default:
    break;
  }

  return femto_event;
}
