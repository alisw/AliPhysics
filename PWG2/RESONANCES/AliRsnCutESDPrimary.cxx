//
// Class AliRsnCutESDPrimary
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
#include "AliRsnCutESDPrimary.h"

ClassImp(AliRsnCutESDPrimary)

//_________________________________________________________________________________________________
AliRsnCutESDPrimary::AliRsnCutESDPrimary() :
    AliRsnCut(),
    fCuts()
{
//
// Default constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutESDPrimary::AliRsnCutESDPrimary
(const char *name) :
    AliRsnCut(name, 0.0, 0.0),
    fCuts()
{
//
// Main constructor.
//
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutESDPrimary::IsSelected(ETarget tgt, AliRsnDaughter * const track)
{
//
// Cut checker.
//

  // coherence check
  if (tgt != AliRsnCut::kParticle) {
    AliError(Form("[%s] Wrong target. Skipping cut", GetName()));
    return kTRUE;
  }

  // retrieve the TPC signal
  AliVParticle *vpart = track->GetRef();
  AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*>(vpart);
  if (!esdTrack) {
    AliError("ESD information unavailable");
    return kTRUE;
  }

  // check cut
  return fCuts.IsSelected(esdTrack);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutESDPrimary::IsSelected(ETarget /*tgt*/, AliRsnPairParticle* /*pair*/)
{
//
// Cut checker
//

  AliWarning("Cannot apply this cut to pairs");
  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutESDPrimary::IsSelected(ETarget /*tgt*/, AliRsnEvent* /*event*/)
{
//
// Cut checker
//

  AliWarning("Cannot apply this cut to events");
  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutESDPrimary::IsSelected(ETarget /*tgt*/, AliRsnEvent* /*ev1*/, AliRsnEvent* /*ev2*/)
{
//
// Cut checker
//

  AliWarning("Cannot apply this cut to event mixing");
  return kTRUE;
}
