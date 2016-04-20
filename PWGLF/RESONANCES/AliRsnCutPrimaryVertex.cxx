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
#include "AliAnalysisUtils.h"

ClassImp(AliRsnCutPrimaryVertex)

//_________________________________________________________________________________________________
AliRsnCutPrimaryVertex::AliRsnCutPrimaryVertex
(const char *name, Double_t maxVz, Int_t nContributors, Bool_t acceptTPC, Bool_t acceptSPD) :
   AliRsnCut(name, AliRsnCut::kEvent, 0, nContributors - 1, -maxVz, maxVz + 1E-6),
   fAcceptTPC(acceptTPC),
   fAcceptSPD(acceptSPD),
   fCheckPileUp(kFALSE),
   fCheckZResolutionSPD(kFALSE),
   fMaxZResolutionSPD(0.25),
   fCheckDispersionSPD(kFALSE),
   fMaxDispersionSPD(0.04),
   fCheckZDifferenceSPDTrack(kFALSE),
   fMaxZDifferenceSPDTrack(0.5)
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
Bool_t AliRsnCutPrimaryVertex::IsSelected(TObject *object)
{
//
// Cut checker
//
   // coherence check
   // which also fills data member objects
   if (!TargetOK(object)) return kFALSE;

   // retrieve ESD event
   AliESDEvent *esd = dynamic_cast<AliESDEvent *>(fEvent->GetRef());
   AliAODEvent *aod = dynamic_cast<AliAODEvent *>(fEvent->GetRef());
   AliVEvent *vevt = dynamic_cast<AliVEvent *>(fEvent->GetRef());
   // pile-up check
   if (fCheckPileUp) {
     AliAnalysisUtils * utils = new AliAnalysisUtils();
     if (utils->IsPileUpSPD(vevt)) return kFALSE;
   }
   
   if (esd) {
      // get the best primary vertex:
      // first try the one with tracks
      const AliESDVertex *vTrk  = esd->GetPrimaryVertexTracks();
      const AliESDVertex *vSPD  = esd->GetPrimaryVertexSPD();
      const AliESDVertex *vTPC  = esd->GetPrimaryVertexTPC();

      if(fCheckZResolutionSPD && vSPD && !GoodZResolutionSPD(vSPD)) return kFALSE;
      if(fCheckDispersionSPD && vSPD && !GoodDispersionSPD(vSPD)) return kFALSE;
      if(fCheckZDifferenceSPDTrack && vTrk && vSPD && !GoodZDifferenceSPDTrack(vTrk, vSPD)) return kFALSE;

      Int_t               ncTrk = -1;
      Int_t               ncSPD = -1;
      Int_t               ncTPC = -1;
      Double_t            vzTrk = 1000000.0;
      Double_t            vzSPD = 1000000.0;
      Double_t            vzTPC = 1000000.0;
      if (vTrk) vzTrk = TMath::Abs(vTrk->GetZ());
      if (vSPD) vzSPD = TMath::Abs(vSPD->GetZ());
      if (vTPC) vzTPC = TMath::Abs(vTPC->GetZ());
      if (vTrk) ncTrk = (Int_t)vTrk->GetNContributors();
      if (vSPD) ncSPD = (Int_t)vSPD->GetNContributors();
      if (vTPC) ncTPC = (Int_t)vTPC->GetNContributors();
      if (vTrk && ncTrk > 0) {
         fCutValueI = ncTrk;
         fCutValueD = vzTrk;
      } else if (vSPD && ncSPD > 0) {
         if (!fAcceptSPD)
            return kFALSE;
         else {
         fCutValueI = ncSPD;
         fCutValueD = vzSPD;
	 }
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
      // in this case, as suggested by Andrea Dainese,
      // we first check if the SPD primary vertex is there
      // if it is not, then the only available is the TPC
      // stand-alone primary vertex, which is rejected
      
      if(fAcceptSPD){
      	AliAODVertex *aodv_SPD = aod->GetPrimaryVertexSPD();
      	if (!aodv_SPD) {
         	AliDebugClass(1, "Not found SPD vertex --> TPC only available, skipped");
         	return kFALSE;
      	}

	if(fCheckZResolutionSPD && !GoodZResolutionSPD(aodv_SPD)) return kFALSE;

      	// now check primary vertex
      	AliAODVertex *aodv_Trk = (AliAODVertex *)aod->GetPrimaryVertex();
      	if (CheckVertex(aodv_Trk)) {
         	AliDebugClass(1, "Vertex TRK is OK");
		if(fCheckZDifferenceSPDTrack && !GoodZDifferenceSPDTrack(aodv_Trk, aodv_SPD)) return kFALSE;
         	fCutValueD = aodv_Trk->GetZ();
         	fCutValueI = aodv_Trk->GetNDaughters(); //aodv->GetNContributors();
      	}
      	else {
         	if (CheckVertex(aodv_SPD)) {
            	AliDebugClass(1, "Vertex TRK is BAD, but vertex SPD is OK");
            	fCutValueD = aodv_SPD->GetZ();
            	fCutValueI = aodv_SPD->GetNDaughters(); //aodv->GetNContributors();
         	} else {
            	AliDebugClass(1, "Vertex TRK is BAD, and vertex SPD is BAD");
            	return kFALSE;
         	}
	 
	 	} 
     }
     else{     
 	 const AliVVertex *vertex = aod->GetPrimaryVertex();
  	 if(!vertex) return kFALSE;
  	 else{
    	 TString title=vertex->GetTitle();
    	 if(title.Contains("Z") ) return kFALSE;
    	 else if(title.Contains("3D") ) return kFALSE;
	 fCutValueI = vertex->GetNContributors();
	 fCutValueD = TMath::Abs(vertex->GetZ());
 	 }
     }
     
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
   AliInfo(Form("Accepting SPD primary vertex : %s", (fAcceptSPD ? "YES" : "NO")));
   AliInfo(Form("Contributors range (outside) : %d - %d", fMinI, fMaxI));
   AliInfo(Form("Z-vertex     range (inside)  : %f - %f", fMinD, fMaxD));
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutPrimaryVertex::CheckVertex(AliVVertex *vertex)
{
//
// Checks if a candidate primary vertex is good,
// which is true if it is not null and has at
// least one contributor
//

   if (!vertex) return kFALSE;
   if (vertex->GetNContributors() < 1) return kFALSE;
   return kTRUE;
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutPrimaryVertex::GoodZResolutionSPD(const AliESDVertex* v)
{
   if(!v) return kTRUE;
   if(v->IsFromVertexerZ() && v->GetZRes()>=fMaxZResolutionSPD) return kFALSE;
   return kTRUE;
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutPrimaryVertex::GoodZResolutionSPD(AliAODVertex* v)
{
   if(!v) return kTRUE;
   Double_t sigma[3];
   v->GetSigmaXYZ(sigma);
   if(sigma[2]<0.) return kFALSE;
   if(TMath::Sqrt(sigma[2])>=fMaxZResolutionSPD) return kFALSE;
   return kTRUE;
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutPrimaryVertex::GoodZResolutionSPD()
{
   if(!fEvent) return kTRUE;

   AliESDEvent* esd=dynamic_cast<AliESDEvent *>(fEvent->GetRef());
   if(esd){
     const AliESDVertex* vesd=esd->GetPrimaryVertexSPD();
     return GoodZResolutionSPD(vesd);
   }

   AliAODEvent* aod=dynamic_cast<AliAODEvent *>(fEvent->GetRef());
   if(aod){
     AliAODVertex* vaod=aod->GetPrimaryVertexSPD();
     return GoodZResolutionSPD(vaod);
   }

   return kTRUE;
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutPrimaryVertex::GoodDispersionSPD(const AliESDVertex *v)
{
   //works for ESDs only
   if(!v) return kTRUE;
   if(v->IsFromVertexerZ() && v->GetDispersion()>=fMaxDispersionSPD) return kFALSE;
   return kTRUE;
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutPrimaryVertex::GoodDispersionSPD()
{
   //works for ESDs only
   if(!fEvent) return kTRUE;
   AliESDEvent* esd=dynamic_cast<AliESDEvent *>(fEvent->GetRef());
   if(!esd) return kTRUE;

   const AliESDVertex* vSPD=esd->GetPrimaryVertexSPD();
   return GoodDispersionSPD(vSPD);
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutPrimaryVertex::GoodZDifferenceSPDTrack(const AliESDVertex *vTrk, const AliESDVertex *vSPD)
{
   if(!vTrk || !vSPD) return kTRUE;
   Double_t vzTrk = vTrk->GetZ();
   Double_t vzSPD = vSPD->GetZ();
   if(TMath::Abs(vzTrk - vzSPD)>fMaxZDifferenceSPDTrack) return kFALSE;
   return kTRUE;
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutPrimaryVertex::GoodZDifferenceSPDTrack(AliAODVertex *vTrk, AliAODVertex *vSPD)
{
   if(!vTrk || !vSPD) return kTRUE;
   Double_t vzTrk = vTrk->GetZ();
   Double_t vzSPD = vSPD->GetZ();
   if(TMath::Abs(vzTrk - vzSPD)>fMaxZDifferenceSPDTrack) return kFALSE;
   return kTRUE;
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutPrimaryVertex::GoodZDifferenceSPDTrack()
{
   if(!fEvent) return kTRUE;
   AliESDEvent* esd=dynamic_cast<AliESDEvent *>(fEvent->GetRef());
   if(esd){
     const AliESDVertex* vTrk_esd=esd->GetPrimaryVertexTracks();
     const AliESDVertex* vSPD_esd=esd->GetPrimaryVertexSPD();
     return GoodZDifferenceSPDTrack(vTrk_esd, vSPD_esd);
   }

   AliAODEvent* aod=dynamic_cast<AliAODEvent *>(fEvent->GetRef());
   if(aod){
     AliAODVertex* vTrk_aod=aod->GetPrimaryVertex();
     AliAODVertex* vSPD_aod=aod->GetPrimaryVertexSPD();
     return GoodZDifferenceSPDTrack(vTrk_aod, vSPD_aod);
   }

   return kTRUE;
}
