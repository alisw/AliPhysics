///
/// \file AliFemtoEventReaderAODMultSelection.cxx
///


#include "AliFemtoEventReaderAODMultSelection.h"
#include "AliFemtoEvent.h"
#include "AliMultSelection.h"

AliFemtoEvent* AliFemtoEventReaderAODMultSelection::CopyAODtoFemtoEvent()
{
  AliFemtoEvent *femto_event = AliFemtoEventReaderAODChain::CopyAODtoFemtoEvent();

  if (femto_event == NULL) {
    return femto_event;
  }

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
  case kCentralityV0A:
    femto_event->SetNormalizedMult(lrint(10 * mult_selection->GetMultiplicityPercentile("V0A")));
    break;
  case kCentralityV0C:
    femto_event->SetNormalizedMult(lrint(10 * mult_selection->GetMultiplicityPercentile("V0C")));
    break;
  case kCentralityZNA:
    femto_event->SetNormalizedMult(lrint(10 * mult_selection->GetMultiplicityPercentile("ZNA")));
    break;
  case kCentralityZNC:
    femto_event->SetNormalizedMult(lrint(10 * mult_selection->GetMultiplicityPercentile("ZNC")));
    break;
  case kCentralityCL0:
    femto_event->SetNormalizedMult(lrint(10 * mult_selection->GetMultiplicityPercentile("CL0")));
    break;
  case kCentralityCL1:
    femto_event->SetNormalizedMult(lrint(10 * mult_selection->GetMultiplicityPercentile("CL1")));
    break;
  default:
    break;
  }

  return femto_event;
}
