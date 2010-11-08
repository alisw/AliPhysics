//
// Class AliRsnCutESDCutMultiplicity
//
// Cuts on event multiplicity computed from number o tracks passing
// the ESDtrackCuts defined as data member.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliRsnEvent.h"
#include "AliRsnCutESDCutMultiplicity.h"

ClassImp(AliRsnCutESDCutMultiplicity)

//_________________________________________________________________________________________________
AliRsnCutESDCutMultiplicity::AliRsnCutESDCutMultiplicity() :
  AliRsnCut(AliRsnCut::kEvent),
  fCuts()
{
//
// Default constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutESDCutMultiplicity::AliRsnCutESDCutMultiplicity
(const char *name, Int_t min, Int_t max) :
  AliRsnCut(name, AliRsnCut::kEvent, min, max),
  fCuts()
{
//
// Main constructor.
//
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutESDCutMultiplicity::IsSelected(TObject *obj1, TObject* /*obj2*/)
{
//
// Cut checker.
//

  // coherence check
  AliRsnEvent *event = dynamic_cast<AliRsnEvent*>(obj1);
  if (!event) return kFALSE;
  AliESDEvent *esd   = event->GetRefESD();
  if (!esd) return kFALSE;
  
  // count the tracks passing the cut
  fCutValueI = fCuts.CountAcceptedTracks(esd);
  return OkRangeI();
}
