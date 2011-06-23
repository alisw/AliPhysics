//
// This cut implements all the checks done to accept a track as a Kaon
// for the PbPb analysis using 2010 runs. 
// It is based on standard cuts on track quality and nsigma cuts
// with respect to the TPC and TOF signals for the PID.
//

#include <Riostream.h>

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliRsnCutKaonForPhi2010.h"

ClassImp(AliRsnCutKaonForPhi2010)

//__________________________________________________________________________________________________
AliRsnCutKaonForPhi2010::AliRsnCutKaonForPhi2010
(const char *name, Double_t nSigmaTPC, Double_t nSigmaTOF, Double_t tofLimit) :
   AliRsnCut(name, AliRsnTarget::kDaughter),
   fOnlyQuality(kFALSE),
   fOnlyTPC(kFALSE),
   fOnlyTOF(kFALSE),
   fCutTPC(nSigmaTPC),
   fCutTOF(nSigmaTOF),
   fTOFthreshold(tofLimit),
   fCutQuality(Form("%s_quality", name))
{
//
// Constructor
// Initialize the contained cuts and sets defaults
//

   // track quality
   //fCutQuality.AddStatusFlag(AliESDtrack::kTPCin   , kTRUE);
   //fCutQuality.AddStatusFlag(AliESDtrack::kTPCrefit, kTRUE);
   //fCutQuality.AddStatusFlag(AliESDtrack::kITSrefit, kTRUE);
   fCutQuality.SetPtRange(0.15, 1E+20);
   fCutQuality.SetEtaRange(-0.8, 0.8);
   fCutQuality.SetDCARPtFormula("0.0182+0.0350/pt^1.01");
   fCutQuality.SetDCAZmax(2.0);
   fCutQuality.SetSPDminNClusters(1);
   fCutQuality.SetITSminNClusters(0);
   fCutQuality.SetITSmaxChi2(1E+20);
   fCutQuality.SetTPCminNClusters(70);
   fCutQuality.SetTPCmaxChi2(4.0);
   fCutQuality.SetRejectKinkDaughters();
   fCutQuality.SetAODTestFilterBit(5);
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutKaonForPhi2010::IsSelected(TObject *obj)
{
//
// Global check
//

   // coherence check
   if (!TargetOK(obj)) return kFALSE;
   
   // check track
   AliVTrack *track = fDaughter->Ref2Vtrack();
   if (!track) {
      if (!fDaughter->GetRef()) AliWarning("NULL ref");
      return kFALSE;
   }
   
   // check flags
   if ((track->GetStatus() & AliESDtrack::kTPCin   ) == 0) return kFALSE;
   if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) return kFALSE;
   if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) return kFALSE;
   
   // quality
   if (!fCutQuality.IsSelected(obj)) return kFALSE;
   
   // if not PID is done, exit here
   if (fOnlyQuality) return kTRUE;
   
   // check initialization of PID object
   AliPIDResponse *pid = fEvent->GetPIDResponse();
   if (!pid) {
      AliFatal("NULL PID response");
      return kFALSE;
   }
   
   // check if TOF is matched
   // and computes all values used in the PID cut
   Bool_t   accept = kFALSE;
   Bool_t   isTOF  = MatchTOF(track);
   Double_t nsTPC  = TMath::Abs(pid->NumberOfSigmasTPC(track, AliPID::kKaon));
   Double_t nsTOF  = TMath::Abs(pid->NumberOfSigmasTOF(track, AliPID::kKaon));
   
   // if only one detector is chosen, do this here
   if (fOnlyTPC) {
      AliDebugClass(1, Form("Checking only TPC: nsigma = %f - cut = %f", nsTPC, fCutTPC));
      accept = (nsTPC <= fCutTPC);
   } else if (fOnlyTOF) {
      if (!isTOF) {
         AliDebugClass(1, "Checking only TOF: rejecting non-TOF track");
         accept = kFALSE;
      } else {
         AliDebugClass(1, Form("Checking only TOF: nsigma = %f - cut = %f", nsTOF, fCutTOF));
         accept = (nsTOF <= fCutTOF);
      }
   } else {
      // combined PID:
      // below momentum threshold, start with TPC
      if (track->P() < fTOFthreshold) {
         AliDebugClass(1, Form("Checking both PID: p = %f below threshold (TOF not required)", track->P())); 
         if (isTOF) {
            AliDebugClass(1, Form("TOF present: nsigmaTPC: = %f - cut = %f", nsTPC, fCutTPC));
            AliDebugClass(1, Form("TOF present: nsigmaTOF: = %f - cut = %f", nsTOF, fCutTOF));
            accept = ((nsTPC <= fCutTPC) && (nsTOF <= fCutTOF));
         } else {
            AliDebugClass(1, Form("TOF absent : nsigmaTPC: = %f - cut = %f", nsTPC, fCutTPC));
            accept = (nsTPC <= fCutTPC);
         }
      } else {
         AliDebugClass(1, Form("Checking both PID: p = %f above threshold (TOF required)", track->P())); 
         if (isTOF) {
            AliDebugClass(1, Form("TOF present: nsigmaTPC: = %f - cut = %f", nsTPC, fCutTPC));
            AliDebugClass(1, Form("TOF present: nsigmaTOF: = %f - cut = %f", nsTOF, fCutTOF));
            accept = ((nsTPC <= fCutTPC) && (nsTOF <= fCutTOF));
         } else {
            AliDebugClass(1, "TOF absent : track rejected");
            accept = kFALSE;
         }
      }
   }
   
   AliDebugClass(1, Form("Track %s", (accept ? "accepted" : "rejected")));
   return accept;
}
