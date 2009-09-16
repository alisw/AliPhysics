//
// Class AliRsnCutPrimaryVertex
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

// #include "AliLog.h"
// #include "AliESDEvent.h"
// #include "AliESDVertex.h"
// 
// #include "AliRsnEvent.h"
// #include "AliRsnCutPrimaryVertex.h"

// #include "Riostream.h"
// #include "TMath.h"
// 
// #include "AliLog.h"
// #include "AliESDEvent.h"
// #include "AliESDVertex.h"
// 
// #include "AliRsnEvent.h"
#include "AliRsnCutPrimaryVertex.h"

ClassImp(AliRsnCutPrimaryVertex)

//_________________________________________________________________________________________________
AliRsnCutPrimaryVertex::AliRsnCutPrimaryVertex() :
    AliRsnCut()
{
//
// Default constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutPrimaryVertex::AliRsnCutPrimaryVertex
(const char *name, Int_t nContributors) :
    AliRsnCut(name, 0, nContributors)
{
//
// Main constructor.
// the cut range is outside the interval 0 - min number of contributors.
//
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPrimaryVertex::IsSelected(AliRsnCut::ETarget /*tgt*/, AliRsnDaughter*/*const track*/)
{
//
// Cut checker.
//
  AliWarning("Cannot apply this cut to particles");
  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPrimaryVertex::IsSelected(AliRsnCut::ETarget, AliRsnPairParticle*/*const pair*/)
{
//
// Cut checker
//

  AliWarning("Cannot apply this cut to pairs");
  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPrimaryVertex::IsSelected(AliRsnCut::ETarget, AliRsnEvent*event)
{
//
// Cut checker
//

  // retrieve ESD event
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(event->GetRef());
  if (!esd) {
    AliDebug(AliLog::kDebug+2, "NO ESD");
    return kTRUE;
  }

  // check primary vertex and eventually fill step 1
  // if a vertex with tracks was successfully reconstructed,
  // it is used for computing DCA;
  // otherwise, the one computed with SPD is used.
  // This is known from the "Status" parameter of the vertex itself.
  const AliESDVertex *v = esd->GetPrimaryVertex();
  if (!v->GetStatus()) v = esd->GetPrimaryVertexSPD();
  if (!v->GetStatus()) {
    AliDebug(AliLog::kDebug+2, "Bad vertex status");
    return kFALSE;
  }

  fCutValueI = v->GetNContributors();

  return !OkRange();
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPrimaryVertex::IsSelected(AliRsnCut::ETarget, AliRsnEvent*/*ev1*/, AliRsnEvent*/*ev2*/)
{
//
// Cut checker
//

  AliWarning("Cannot apply this cut to event mixing");
  return kTRUE;
}
