//
// Class AliRsnCutPrimaryVertex
//
// This cut implementation checks the quality of event primary vertex.
// It currently works only with ESD events (not AOD).
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliRsnCutPrimaryVertex.h"

ClassImp(AliRsnCutPrimaryVertex)

//_________________________________________________________________________________________________
AliRsnCutPrimaryVertex::AliRsnCutPrimaryVertex() :
    AliRsnCut(),
    fAcceptTPC(kFALSE)
{
//
// Default constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutPrimaryVertex::AliRsnCutPrimaryVertex
(const char *name, Int_t nContributors, Bool_t acceptTPC) :
    AliRsnCut(name, 0, nContributors - 1),
    fAcceptTPC(acceptTPC)
{
//
// Main constructor.
// Defines the cut range between 0 and
// the minimum required number of contributors.
// The cut will be passed when if the event has a
// primary vertex with number of contributors outside this interval.
// ---
// If the 'acceptTPC' argument is true, events with TPC
// primary vertex will be checked, otherwise they will be
// rejected by default.
// ---
// Since the range check uses the '>=' and '<=', the high edge
// must be decreased by 1 to get the right behaviour, since this is integer.
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
    return kFALSE;
  }

  // get the three possible primary vertexes of the event:
  // 0 = default one
  // 1 = SPD
  // 2 = TPC
  // then, if the default vertex is TPC, the event is rejected,
  // otherwise, the event is rejected only if its vertex status is 'false'
  // get primary vertexes
  const AliESDVertex *vert[3];
  vert[0] = esd->GetPrimaryVertex();
  vert[1] = esd->GetPrimaryVertexSPD();
  vert[2] = esd->GetPrimaryVertexTPC();

  // if TPC vertexes are rejected by default do this now
  if (!fAcceptTPC && (vert[0] == vert[2])) {
    AliDebug(AliLog::kDebug+2, "Rejecting TPC vertex");
    return kFALSE;
  }

  // if we are here, vert[0] contains the default primary vertex
  // in case it is with tracks or SPD, its status must be OK
  // because otherwise the ESD event returns the lower level vertex
  // then, we can check just the first element in the array
  if (!vert[0]->GetStatus()) {
    AliDebug(AliLog::kDebug+2, "Bad vertex status");
    return kFALSE;
  }

  // if we are here, the status of default vertex is OK
  // and then we check the number of contributors:
  // it must be *outside* the 'allowed' range
  fCutValueI = vert[0]->GetNContributors();
  return /*there is a NOT operator */!/*here*/OkRange();
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
