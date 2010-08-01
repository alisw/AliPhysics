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
  AliRsnCut(AliRsnCut::kEvent),
  fAcceptTPC(kFALSE)
{
//
// Default constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutPrimaryVertex::AliRsnCutPrimaryVertex
(const char *name, Double_t maxVz, Int_t nContributors, Bool_t acceptTPC) :
  AliRsnCut(name, AliRsnCut::kEvent, 0, nContributors - 1),
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

  fMinD = 0.0;
  fMaxD = maxVz;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPrimaryVertex::IsSelected(TObject *obj1, TObject* /*obj2*/)
{
//
// Cut checker
//

  // retrieve ESD event
  AliRsnEvent *rsn = dynamic_cast<AliRsnEvent*>(obj1);
  if (!rsn) return kFALSE;
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(rsn->GetRef());
  if (!esd) 
  {
    AliDebug(AliLog::kDebug+2, "NO ESD");
    return kFALSE;
  }

  // get the best primary vertex:
  // first try the one with tracks
  const AliESDVertex *vTrk  = esd->GetPrimaryVertexTracks();
  const AliESDVertex *vSPD  = esd->GetPrimaryVertexSPD();
  const AliESDVertex *vTPC  = esd->GetPrimaryVertexTPC();
  Double_t            vzTrk = 2.0 * fMaxD;
  Double_t            vzSPD = 2.0 * fMaxD;
  Double_t            vzTPC = 2.0 * fMaxD;
  if (vTrk) vzTrk = TMath::Abs(vTrk->GetZv());
  if (vSPD) vzSPD = TMath::Abs(vSPD->GetZv());
  if (vTPC) vzTPC = TMath::Abs(vTPC->GetZv());
  AliDebug(AliLog::kDebug + 1, Form("Vertex with tracks: contributors = %d, abs(vz) = %f", vTrk->GetNContributors(), vzTrk));
  AliDebug(AliLog::kDebug + 1, Form("Vertex with SPD,    contributors = %d, abs(vz) = %f", vSPD->GetNContributors(), vzSPD));
  AliDebug(AliLog::kDebug + 1, Form("Vertex with TPC,    contributors = %d, abs(vz) = %f", vTPC->GetNContributors(), vzTPC));
  if(vTrk->GetStatus())
  {
    fCutValueI = vTrk->GetNContributors();
    fCutValueD = vzTrk;
  }
  else if (vSPD->GetStatus())
  {
    fCutValueI = vSPD->GetNContributors();
    fCutValueD = vzSPD;
  }
  else if (vTPC->GetStatus() && fAcceptTPC)
  {
    fCutValueI = vTPC->GetNContributors();
    fCutValueD = vzTPC;
  }
  else
    return kFALSE;
  
  return OkRangeI() && OkRangeD();
}

