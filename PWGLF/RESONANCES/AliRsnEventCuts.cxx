//
// Class AliRsnEventCuts
//
// This class constitutes an interface to the AliEventCuts class 
// supported and maintained by the DPG starting from LHC Run 2
// for period-by-period standard cuts on events, including 
// vertex cuts, background and pileup rejection, etc...
//
// authors: P. Ganoti, F. Bellini
//

#include "AliRsnEventCuts.h"
#include "AliEventCuts.h"

ClassImp(AliRsnEventCuts)

//_________________________________________________________________________________________________
AliRsnEventCuts::AliRsnEventCuts(const char *name) :
  AliRsnCut(name, AliRsnCut::kEvent),
  fEvCuts()
{
  //
  // Main constructor.
  //
  
}

//_________________________________________________________________________________________________
Bool_t AliRsnEventCuts::Init(TObject *object)
{ 
  // fills data member objects without applying cuts
  if (!TargetOK(object)) return kFALSE;
  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnEventCuts::IsSelected(TObject *object)
{
//
// Cut checker
//
   // coherence check
   // which also fills data member objects
   if (!Init(object)) return kFALSE;
   // retrieve event
   AliVEvent *vevt = dynamic_cast<AliVEvent *>(fEvent->GetRef());
   if (!vevt) return kFALSE;
   return fEvCuts.AcceptEvent(vevt);
}
