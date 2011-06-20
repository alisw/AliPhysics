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
   Bool_t   isTOF  = MatchTOF(track);
   Double_t nsTPC  = TMath::Abs(pid->NumberOfSigmasTPC(track, AliPID::kKaon));
   Double_t nsTOF  = TMath::Abs(pid->NumberOfSigmasTOF(track, AliPID::kKaon));
   
   // if only one detector is chosen, do this here
   if (fOnlyTPC) {
      return (nsTPC <= fCutTPC);
   } else if (fOnlyTOF) {
      return (isTOF && (nsTOF <= fCutTOF));
   } else {
      // combined PID:
      // below momentum threshold, start with TPC
      if (track->P() < fTOFthreshold) {
         if (isTOF)
            return ((nsTPC <= fCutTPC) && (nsTOF <= fCutTOF));
         else
            return (nsTPC <= fCutTPC);
      } else {
         if (isTOF) 
            return ((nsTPC <= fCutTPC) && (nsTOF <= fCutTOF));
         else
            return kFALSE;
      }
   }
}
