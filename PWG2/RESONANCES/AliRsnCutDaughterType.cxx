//
// Class AliRsnCutDaughterType
//
// General implementation of a single cut strategy, which can be:
// - a value contained in a given interval  [--> IsBetween()   ]
// - a value equal to a given reference     [--> MatchesValue()]
//
// In all cases, the reference value(s) is (are) given as data members
// and each kind of cut requires a given value type (Int, UInt, Double),
// but the cut check procedure is then automatized and chosen thanks to
// an enumeration of the implemented cut types.
// At the end, the user (or any other point which uses this object) has
// to use the method IsSelected() to check if this cut has been passed.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliRsnDaughter.h"
#include "AliRsnCutDaughterType.h"

ClassImp(AliRsnCutDaughterType)

//_________________________________________________________________________________________________
AliRsnCutDaughterType::AliRsnCutDaughterType() :
  AliRsnCut(),
  fRefType(kTypes)
{
//
// Default constructor.
//

  SetTargetType(AliRsnTarget::kDaughter);
}

//_________________________________________________________________________________________________
AliRsnCutDaughterType::AliRsnCutDaughterType
(const char *name, EType type) :
  AliRsnCut(name, AliRsnCut::kDaughter, 0.0, 0.0),
  fRefType(type)
{
//
// Main constructor.
//
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutDaughterType::IsSelected(TObject *object)
{
//
// Cut checker.
//

  // coherence check
  if (!TargetOK(object)) return kFALSE;
  
  // check the daughter according to the selected type
  // in some cases this means to retrieve the track status
  AliRsnDaughter *daughter = dynamic_cast<AliRsnDaughter*>(object);
  AliVTrack   *track  = dynamic_cast<AliVTrack*>(daughter->GetRef());
  AliESDtrack *esdT   = dynamic_cast<AliESDtrack*>(daughter->GetRef());
  ULong_t      status = 0x0;
  if (track) status = (ULong_t)track->GetStatus();
  
  switch (fRefType)
  {
    case kTrackTPC:
      return ((status & AliESDtrack::kTPCin)  != 0);
    case kTrackITSSA:
      if (esdT && track)
      {
        UChar_t itsCluMap = track->GetITSClusterMap();
        Int_t   k, nITS   = 0;
        for(k = 2; k < 6; k++) if(itsCluMap & (1 << k)) ++nITS;
        if (nITS < 3) return kFALSE;
      }
      return ((status & AliESDtrack::kTPCin)  == 0 && (status & AliESDtrack::kITSrefit) != 0 && (status & AliESDtrack::kITSpureSA) == 0 && (status & AliESDtrack::kITSpid) != 0);
    case kV0:
      return daughter->IsV0();
    default:
      AliError("No good reference type is chosen. Cut skipped");
      return kTRUE;
  }
}
