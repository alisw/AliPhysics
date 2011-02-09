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
AliRsnCutPrimaryVertex::AliRsnCutPrimaryVertex
(const char *name, Double_t maxVz, Int_t nContributors, Bool_t acceptTPC) :
   AliRsnCut(name, AliRsnCut::kEvent, 0, nContributors - 1, 0.0, maxVz),
   fAcceptTPC(acceptTPC),
   fCheckPileUp(kFALSE)
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
Bool_t AliRsnCutPrimaryVertex::IsSelected(TObject *object)
{
//
// Cut checker
//

   // retrieve ESD event
   AliRsnEvent *rsn = dynamic_cast<AliRsnEvent*>(object);
   if (!rsn) return kFALSE;
   AliESDEvent *esd = rsn->GetRefESD();
   AliAODEvent *aod = rsn->GetRefAOD();

   if (esd) {
      // pile-up check
      if (fCheckPileUp) {
         if (esd->IsPileupFromSPD()) return kFALSE;
      }

      // get the best primary vertex:
      // first try the one with tracks
      const AliESDVertex *vTrk  = esd->GetPrimaryVertexTracks();
      const AliESDVertex *vSPD  = esd->GetPrimaryVertexSPD();
      const AliESDVertex *vTPC  = esd->GetPrimaryVertexTPC();
      Int_t               ncTrk = -1;
      Int_t               ncSPD = -1;
      Int_t               ncTPC = -1;
      Double_t            vzTrk = 1000000.0;
      Double_t            vzSPD = 1000000.0;
      Double_t            vzTPC = 1000000.0;
      if (vTrk) vzTrk = TMath::Abs(vTrk->GetZv());
      if (vSPD) vzSPD = TMath::Abs(vSPD->GetZv());
      if (vTPC) vzTPC = TMath::Abs(vTPC->GetZv());
      if (vTrk) ncTrk = (Int_t)vTrk->GetNContributors();
      if (vSPD) ncSPD = (Int_t)vSPD->GetNContributors();
      if (vTPC) ncTPC = (Int_t)vTPC->GetNContributors();
      if (vTrk && ncTrk > 0) {
         fCutValueI = ncTrk;
         fCutValueD = vzTrk;
      } else if (vSPD && ncSPD > 0) {
         fCutValueI = ncSPD;
         fCutValueD = vzSPD;
      } else if (vTPC && ncTPC > 0) {
         if (!fAcceptTPC)
            return kFALSE;
         else {
            fCutValueI = ncTPC;
            fCutValueD = vzTPC;
         }
      } else
         return kFALSE;
   } else if (aod) {
      // pile-up check is not yet available for AODs

      // lines suggested by Andrea to reject TPC-only events
      if (!fAcceptTPC) {
         if (!aod->GetPrimaryVertexSPD()) return kFALSE;
         else if (aod->GetPrimaryVertexSPD()->GetNContributors() < 1) return kFALSE;
      }

      AliAODVertex *prim = (AliAODVertex*)aod->GetPrimaryVertex();
      if (!prim) return kFALSE;

      fCutValueI = prim->GetNContributors();
      fCutValueD = prim->GetZ();
   } else
      return kFALSE;

   // output
   Bool_t result = ((!OkRangeI()) && OkRangeD());
   return result;
}

//_________________________________________________________________________________________________
void AliRsnCutPrimaryVertex::Print(const Option_t *) const
{
//
// Print information on this cut
//

   AliInfo(Form("Cut name                     : %s", GetName()));
   AliInfo(Form("Accepting TPC primary vertex : %s", (fAcceptTPC ? "YES" : "NO")));
   AliInfo(Form("Contributors range (outside) : %d - %d", fMinI, fMaxI));
   AliInfo(Form("Z-vertex     range (inside)  : %f - %f", fMinD, fMaxD));
}
