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
  fMaxD = maxVz + 1E-6;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPrimaryVertex::IsSelected(TObject *obj1, TObject* /*obj2*/)
{
//
// Cut checker
//

  static Int_t evNum = 0;
  evNum++;
  
  // retrieve ESD event
  AliRsnEvent *rsn = dynamic_cast<AliRsnEvent*>(obj1);
  if (!rsn) return kFALSE;
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(rsn->GetRef());
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(rsn->GetRef());
  
  if (esd)
  {
    // get the best primary vertex:
    // first try the one with tracks
    const AliESDVertex *vTrk  = esd->GetPrimaryVertexTracks();
    const AliESDVertex *vSPD  = esd->GetPrimaryVertexSPD();
    const AliESDVertex *vTPC  = esd->GetPrimaryVertexTPC();
    Int_t               ncTrk = -1;
    Int_t               ncSPD = -1;
    Int_t               ncTPC = -1;
    Double_t            vzTrk = 2.0 * fMaxD;
    Double_t            vzSPD = 2.0 * fMaxD;
    Double_t            vzTPC = 2.0 * fMaxD;
    if (vTrk) vzTrk = TMath::Abs(vTrk->GetZv());
    if (vSPD) vzSPD = TMath::Abs(vSPD->GetZv());
    if (vTPC) vzTPC = TMath::Abs(vTPC->GetZv());
    if (vTrk) ncTrk = (Int_t)vTrk->GetNContributors();
    if (vSPD) ncSPD = (Int_t)vSPD->GetNContributors();
    if (vTPC) ncTPC = (Int_t)vTPC->GetNContributors();
    if(vTrk && ncTrk > 0)
    {
      fCutValueI = ncTrk;
      fCutValueD = vzTrk;
    }
    else if (vSPD && ncSPD > 0)
    {
      fCutValueI = ncSPD;
      fCutValueD = vzSPD;
    }
    else if (vTPC && fAcceptTPC && ncTPC > 0)
    {
      fCutValueI = ncTPC;
      fCutValueD = vzTPC;
    }
    else
    {
      fCutValueI = -1;
      fCutValueD = 2.0 * fMaxD;
    }
  }
  else if (aod)
  {
    // lines suggested by Andrea to reject TPC-only events
    if(!aod->GetPrimaryVertexSPD()) return kFALSE;
    if(aod->GetPrimaryVertexSPD()->GetNContributors() < 1) return kFALSE;
    
    AliAODVertex *prim = (AliAODVertex*)aod->GetPrimaryVertex();
    fCutValueI = prim->GetNContributors();
    fCutValueD = prim->GetZ();
  }
  else
    return kFALSE;
    
  // output
  Bool_t result = ((!OkRangeI()) && OkRangeD());
  return result;
}

