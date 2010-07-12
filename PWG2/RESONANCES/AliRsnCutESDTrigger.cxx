//
// Class AliRsnCutESDTrigger
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

#include "AliESDEvent.h"
#include "AliRsnEvent.h"
#include "AliRsnCutESDTrigger.h"

ClassImp(AliRsnCutESDTrigger)

//_________________________________________________________________________________________________
AliRsnCutESDTrigger::AliRsnCutESDTrigger() :
    AliRsnCut(),
    fTrigger()
{
//
// Default constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutESDTrigger::AliRsnCutESDTrigger
(const char *name, const char *mask) :
    AliRsnCut(name, 0.0, 0.0),
    fTrigger(mask)
{
//
// Main constructor.
//
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutESDTrigger::IsSelected(ETarget /*tgt*/, AliRsnDaughter* /*const track*/)
{
//
// Cut checker.
//

  AliWarning("Cannot apply this cut to tracks");
  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutESDTrigger::IsSelected(ETarget /*tgt*/, AliRsnPairParticle* /*pair*/)
{
//
// Cut checker
//

  AliWarning("Cannot apply this cut to pairs");
  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutESDTrigger::IsSelected(ETarget tgt, AliRsnEvent* const event)
{
//
// Cut checker
//

  // coherence check
  if (tgt != AliRsnCut::kEvent) {
    AliError(Form("[%s] Wrong target. Skipping cut", GetName()));
    return kTRUE;
  }

  // retrieve the event trigger mask
  AliVEvent *vevent = event->GetRef();
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(vevent);
  if (!esd) {
    AliError("ESD information unavailable");
    return kTRUE;
  }
  
  // check trigger mask
  return esd->IsTriggerClassFired(fTrigger.Data());
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutESDTrigger::IsSelected(ETarget /*tgt*/, AliRsnEvent* /*ev1*/, AliRsnEvent* /*ev2*/)
{
//
// Cut checker
//

  AliWarning("Cannot apply this cut to event mixing");
  return kTRUE;
}
