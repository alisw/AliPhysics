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
//#include "AliESDtrackCuts.h"
#include "AliMultSelection.h"
//#include "AliMCEvent.h"
//#include "AliMCParticle.h"
//#include "AliAODMCParticle.h"
//#include <AliHeader.h>
//#include <AliAODMCHeader.h>
//#include <AliGenDPMjetEventHeader.h>

ClassImp(AliRsnEventCuts)

//_________________________________________________________________________________________________
AliRsnEventCuts::AliRsnEventCuts(const char *name) :
  AliRsnCut(name, AliRsnCut::kEvent),
  fEvCuts(0x0)
{
  //
  // Main constructor.
  //
  
}

//_________________________________________________________________________________________________
AliRsnEventCuts::AliRsnEventCuts(const AliRsnEventCuts &copy) :
  AliRsnCut(copy),
  fEvCuts(copy.fEvCuts)
{
  //
  // Copy constructor.
  //
}

//-------------------------------------------------------------------------------------------------
AliRsnEventCuts &AliRsnEventCuts::operator=(const AliRsnEventCuts &copy)
{
  //
  // Assignment operator.
  // Works like copy constructor.
  //
  AliRsnCut::operator=(copy);
  if (this == &copy)
    return *this;
 
  fEvCuts = copy.fEvCuts;
	
  return (*this);
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
   fEvCuts = new AliEventCuts();
   if (fUsePbPb2018) fEvCuts->SetupPbPb2018();

   Bool_t accept = kTRUE;
   if(!IsAcceptedMultSelection()) return kFALSE;
   if (!fEvCuts->AcceptEvent(vevt)) return kFALSE;
  
   return accept;
}

//_________________________________________________________________________________________________
Bool_t AliRsnEventCuts::IsAcceptedMultSelection() {

  AliMultSelection *MultSelection=0;

  if(!fEvent) return kFALSE;

  AliESDEvent* esdEvt=0;
  Bool_t isESD=fEvent->IsESD();
  if(isESD) esdEvt=dynamic_cast<AliESDEvent *>(fEvent->GetRef());
  if(isESD && esdEvt){
    MultSelection=(AliMultSelection*) esdEvt->FindListObject("MultSelection");
    if(!MultSelection) return kTRUE;
    if(MultSelection->IsEventSelected()) return kTRUE;
    return kFALSE;
  }

  return kFALSE;
}
