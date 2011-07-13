//
// This cut implements all the checks done to accept a track as a Kaon
// for the PbPb analysis using 2010 runs. 
// It is based on standard cuts on track quality and nsigma cuts
// with respect to the TPC and TOF signals for the PID.
//

#include <Riostream.h>

#include "AliLog.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliESDpid.h"
#include "AliAODpidUtil.h"
#include "AliRsnCutKaonForPhi2010.h"

ClassImp(AliRsnCutKaonForPhi2010)

//__________________________________________________________________________________________________
AliRsnCutKaonForPhi2010::AliRsnCutKaonForPhi2010
(const char *name, Double_t nSigmaTPC, Double_t nSigmaTOF, Double_t tofLimit) :
   AliRsnCut(name, AliRsnTarget::kDaughter),
   fMode(kDefaultPID),
   fCutTPC(nSigmaTPC),
   fCutTOF(nSigmaTOF),
   fTOFthreshold(tofLimit),
   fMyPID(0x0),
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
AliRsnCutKaonForPhi2010::AliRsnCutKaonForPhi2010(const AliRsnCutKaonForPhi2010 &copy) :
   AliRsnCut(copy),
   fMode(copy.fMode),
   fCutTPC(copy.fCutTPC),
   fCutTOF(copy.fCutTOF),
   fTOFthreshold(copy.fTOFthreshold),
   fMyPID(copy.fMyPID),
   fCutQuality(copy.fCutQuality)
{
//
// Copy constructor
//
}

//__________________________________________________________________________________________________
AliRsnCutKaonForPhi2010& AliRsnCutKaonForPhi2010::operator=(const AliRsnCutKaonForPhi2010 &copy)
{
//
// Assignment operator
//

   AliRsnCut::operator=(copy);
   fMode = copy.fMode;
   fCutTPC = copy.fCutTPC;
   fCutTOF = copy.fCutTOF;
   fTOFthreshold = copy.fTOFthreshold;
   fMyPID = copy.fMyPID;
   fCutQuality = copy.fCutQuality;
   
   return *this;
}

//__________________________________________________________________________________________________
void AliRsnCutKaonForPhi2010::InitMyPID(Bool_t isMC, Bool_t isESD)
{
//
// Initialize manual PID object
//

   if (isESD) 
      fMyPID = new AliESDpid(isMC);
   else
      fMyPID = new AliAODpidUtil(isMC);
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
   
   // check quality and reject always bad quality tracks
   if (!fCutQuality.IsSelected(obj)) {
      AliDebugClass(1, Form("[%s] Track quality is bad", GetName()));
      return kFALSE;
   }
   
   // initialize check variables
   Bool_t   accept = kFALSE;
   Bool_t   isTOF  = MatchTOF(track);
   Double_t nsTPC  = 1E20;
   Double_t nsTOF  = 1E20;
   
   // if PID is required, compute it check initialization of PID object
   if (fMode >= kOnlyTPC && fMode <= kDefaultPID) {
      AliPIDResponse *pid = fEvent->GetPIDResponse();
      if (!pid) {
         AliFatal("NULL PID response");
         return kFALSE;
      }
      // TPC PID
      if (fMyPID) 
         nsTPC = TMath::Abs(fMyPID->NumberOfSigmasTPC(track, AliPID::kKaon));
      else
         nsTPC = TMath::Abs(pid->NumberOfSigmasTPC(track, AliPID::kKaon));
      // TOF PID
      nsTOF = TMath::Abs(pid->NumberOfSigmasTOF(track, AliPID::kKaon));
   }
   
   // decide cut result depending on mode
   switch (fMode) {
      case kQuality:
         // in this case, since bad quality tracks are rejected above,
         // all tracks arrived here can be accepted
         AliDebugClass(1, Form("[%s] Required only track quality", GetName()));
         accept = kTRUE;
         break;
      case kOnlyTPC:
         // in this case, only TPC PID is checked
         // all tracks have one, so nothing additional is checked
         AliDebugClass(1, Form("[%s] Checking only TPC: nsigma = %f - cut = %f", GetName(), nsTPC, fCutTPC));
         accept = (nsTPC <= fCutTPC);
         break;
      case kOnlyTOF:
         // in this case, only TOF PID is checked
         // additional check: we want that TOF is matched
         AliDebugClass(1, Form("[%s] Checking only TOF: nsigma = %f - cut = %f", GetName(), nsTOF, fCutTOF));
         if (isTOF) {
            accept = (nsTOF <= fCutTOF);
         } else {
            AliDebugClass(1, Form("[%s] TOF not matched", GetName()));
            accept = kFALSE;
         }
         break;
      case kDefaultPID:
         // in this case, default PID check is done
         // TPC PID is checked and tracks are rejected if it is not passed
         // if their momentum is below the TOF threshold, they are required
         // to be matched in TOF, otherwise TPC only is OK
         AliDebugClass(1, Form("[%s] Default PID TPC: nsigma = %f - cut = %f", GetName(), nsTPC, fCutTPC));
         AliDebugClass(1, Form("[%s] Default PID TOF: nsigma = %f - cut = %f", GetName(), nsTOF, fCutTOF));
         // step 0: check TPC
         if (nsTPC > fCutTPC) {
            AliDebugClass(1, Form("[%s] TPC PID cut not passed", GetName()));
            accept = kFALSE;
         } else {
            if (isTOF) {
               accept = (nsTOF <= fCutTOF);
            } else {
               if (track->P() >= fTOFthreshold) {
                  AliDebugClass(1, Form("[%s] p = %f above threshold: TOF is required but missing", GetName(), track->P()));
                  accept = kFALSE;
               } else {
                  AliDebugClass(1, Form("[%s] p = %f below threshold: TOF is not required", GetName(), track->P()));
                  accept = kTRUE;
               }
            }
         }
         break;
      default:
         AliDebugClass(1, Form("[%s] Wrong mode", GetName()));
         accept = kFALSE;
   }
   
   AliDebugClass(1, Form("[%s] Track %s", GetName(), (accept ? "accepted" : "rejected")));
   return accept;
}
