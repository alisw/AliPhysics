//
// All cuts for single kaons in phi analysis 2010
//

#include <Riostream.h>

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliRsnCutKaonForPhi2010.h"

ClassImp(AliRsnCutKaonForPhi2010)

//__________________________________________________________________________________________________
AliRsnCutKaonForPhi2010::AliRsnCutKaonForPhi2010(const char *name) :
   AliRsnCut(name, AliRsnTarget::kDaughter, -3.0, 3.0),
   fCutQuality(Form("%sQuality", name))
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
   
   // check initialization of PID object
   AliPIDResponse *pid = fEvent->GetPIDResponse();
   if (!pid) {
      AliFatal("NULL PID response");
      return kFALSE;
   }
   
   // PID ITS :
   // depends on momentum
   //SetRangeD(0.0, 4.0);
   //fCutValueD = TMath::Abs(pid->NumberOfSigmasITS(track, AliPID::kKaon));
   //if (!OkRangeD()) return kFALSE;
   
   // PID TPC :
   // depends on momentum
   //SetRangeD(0.0, 3.0);
   //if (track->GetTPCmomentum() < 0.350) SetRangeD(0.0, 5.0);
   fCutValueD = TMath::Abs(pid->NumberOfSigmasTPC(track, AliPID::kKaon));
   if (!OkRangeD()) return kFALSE;
   
   // if TOF is not matched, end here
   // otherwise check TOF
   if (!MatchTOF(track)) 
      return kTRUE;
   else {
      //SetRangeD(0.0, 3.0);
      fCutValueD = TMath::Abs(pid->NumberOfSigmasTOF(track, AliPID::kKaon));
      return OkRangeD();
   }
}
