//
// Class AliRsnValueMult
//
// Inherits from AliRsnValue, and computer multiplicity of the event
// in several ways
//
// Author: A. Pulvirenti
// Email : alberto.pulvirenti@ct.infn.it
//

#include "AliESDEvent.h"
#include "AliRsnEvent.h"
#include "AliRsnValueMult.h"

ClassImp(AliRsnValueMult)

//_____________________________________________________________________________
AliRsnValueMult::AliRsnValueMult() :
  AliRsnValue(),
  fMode(kESDcuts),
  fESDcuts()
{
//
// Main constructor (version 1)
// This can also be created without any argument.
//
}

//_____________________________________________________________________________
AliRsnValueMult::AliRsnValueMult
(const char *name, EValueType type, Int_t nbins, Double_t min, Double_t max) :
  AliRsnValue(name, type, nbins, min, max),
  fMode(kESDcuts),
  fESDcuts()
{
//
// Main constructor (version 1)
// This can also be created without any argument.
//

  SetBins(nbins, min, max);
}

//_____________________________________________________________________________
AliRsnValueMult::AliRsnValueMult
(const char *name, EValueType type, Double_t min, Double_t max, Double_t step) :
  AliRsnValue(name, type, min, max, step),
  fMode(kESDcuts),
  fESDcuts()
{
//
// Main constructor (version 2)
//

  SetBins(min, max, step);
}

//_____________________________________________________________________________
AliRsnValueMult::AliRsnValueMult
(const char *name, EValueType type, Int_t nbins, Double_t *array) :
  AliRsnValue(name, type, nbins, array),
  fMode(kESDcuts),
  fESDcuts()
{
//
// Main constructor (version 2)
//

  SetBins(nbins, array);
}

//_____________________________________________________________________________
Bool_t AliRsnValueMult::Eval(AliRsnMother * const, AliRsnPairDef * const, AliRsnEvent * const event)
{
//
// Evaluation of the required value.
// It is implemented like this for compatibility with mother class
// but it uses just the passed event
//

  // avoid segfaults
  if (!event) return kFALSE;
  
  // if using method 1, we need to convert it into ESD event
  // and count how many tracks do pass the ESD selections
  AliESDEvent *esd = event->GetRefESD();
  if (!esd && fMode == kESDcuts)
  {
    AliError("Cannot use method based on ESD cuts when input is not ESD.");
    return kFALSE;
  }
  
  // loopf on modes
  switch (fMode)
  {
    case kESDcuts:
      fValue = fESDcuts.CountAcceptedTracks(esd);
      break;
    case kNTracks:
      fValue = event->GetRef()->GetNumberOfTracks();
      break;
    case kNTracklets:
      fValue = -1.0;
      AliWarning("Not yet implemented counting of tracklets");
      break;
    default:
      fValue = -1.0;
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnValueMult::Eval(AliRsnDaughter * const, AliRsnEvent * const event)
{
//
// Works as above, using only the event
//

  return Eval(0x0, 0x0, event);
}
