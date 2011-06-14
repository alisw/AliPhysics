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
AliRsnCutKaonForPhi2010::AliRsnCutKaonForPhi2010(const char *name) :
   AliRsnCut(name, AliRsnTarget::kDaughter, 0.0, 3.0),
   fOnlyQuality(kFALSE),
   fOnlyTPC(kFALSE),
   fOnlyTOF(kFALSE),
   fCutTPC(2.0),
   fCutTOF(2.0),
   fTOFthreshold(0.8),
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
   
   // PID 
   Double_t nsigmaTPC = TMath::Abs(pid->NumberOfSigmasTPC(track, AliPID::kKaon));
   Double_t nsigmaTOF = TMath::Abs(pid->NumberOfSigmasTOF(track, AliPID::kKaon));
   Bool_t   matchTOF  = MatchTOF(track);
   
   // if only one detector is chosen, do this here
   if (fOnlyTPC) {
      return (nsigmaTPC <= fCutTPC);
   } else if (fOnlyTOF) {
      return (matchTOF && (nsigmaTOF <= fCutTOF));
   } else {
      // combined PID:
      // below momentum threshold, start with TPC
      if (track->P() < fTOFthreshold) {
         if (nsigmaTPC > fCutTPC) return kFALSE;
         if (matchTOF && (nsigmaTOF > fCutTOF)) return kFALSE;
         return kTRUE;
      } else {
         if (!matchTOF) return kFALSE;
         if (nsigmaTOF > fCutTOF) return kFALSE;
         if (nsigmaTPC > 4.0) return kFALSE;
      }
      return kTRUE;
   }
}
