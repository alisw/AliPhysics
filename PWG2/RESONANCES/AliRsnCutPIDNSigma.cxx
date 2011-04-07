//
// Class AliRsnCutPIDNSigma
//
// General implementation of a single cut strategy, which can be:
// - a value contained in a given interval  [--> IsBetween()   ]
// - a value equal to a given reference     [--> MatchesValue()]
//
// In all cases, the reference value(s) is (are) given as data members
// and each kind of cut requires a given value type (Int, UInt, Double),
// but the cut check procedure is then automatized and chosen thanks to
// an enumeration of the implemented cut types.
// At the end, the user (or any other point which uses this object) has
// to use the method IsSelected() to check if this cut has been passed.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliPIDResponse.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMultiInputEventHandler.h"

#include "AliRsnCutPIDNSigma.h"

ClassImp(AliRsnCutPIDNSigma)

//_________________________________________________________________________________________________
AliRsnCutPIDNSigma::AliRsnCutPIDNSigma
(const char *name, AliPID::EParticleType species, EDetector det, Double_t nsigma) :
   AliRsnCut(name, AliRsnCut::kDaughter, -nsigma, nsigma),
   fSpecies(species),
   fDetector(det),
   fMomMin(0.0),
   fMomMax(0.0),
   fRejectOutside(kFALSE),
   fRejectUnmatched(kFALSE)
{
//
// Main constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutPIDNSigma::AliRsnCutPIDNSigma
(const AliRsnCutPIDNSigma& copy) :
   AliRsnCut(copy),
   fSpecies(copy.fSpecies),
   fDetector(copy.fDetector),
   fMomMin(copy.fMomMin),
   fMomMax(copy.fMomMax),
   fRejectOutside(copy.fRejectOutside),
   fRejectUnmatched(copy.fRejectUnmatched)
{
//
// Copy constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutPIDNSigma& AliRsnCutPIDNSigma::operator=(const AliRsnCutPIDNSigma& copy)
{
//
// Assignment operator
//

   AliRsnCut::operator=(copy);

   fSpecies = copy.fSpecies;
   fDetector = copy.fDetector;
   fMomMin = copy.fMomMin;
   fMomMax = copy.fMomMax;
   fRejectOutside = copy.fRejectOutside;
   fRejectUnmatched = copy.fRejectUnmatched;

   return (*this);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPIDNSigma::IsSelected(TObject *object)
{
//
// Cut checker.
//

   // coherence check
   if (!TargetOK(object)) return kFALSE;
   
   // check initialization of PID object
   AliPIDResponse *pid = fEvent->GetPIDResponse();
   if (!pid) {
      AliFatal("NULL PID response");
      return kFALSE;
   }

   // get reference momentum, for range cut
   Double_t momentum = -1.0;
   if (!fDaughter->GetRefVtrack()) {
      AliDebugClass(2, "Referenced daughter is not a track");
      return kFALSE;
   }
   if (fDetector == kTPC)
      momentum = fDaughter->GetRefVtrack()->GetTPCmomentum();
   else
      momentum = fDaughter->GetRef()->P();
      
   // check momentum range, if required
   if (fRejectOutside) {
      if (momentum < fMomMin || momentum > fMomMax) {
         AliDebugClass(2, Form("Track momentum = %.5f, outside allowed range [%.2f - %.2f]", momentum, fMomMin, fMomMax));
         return kFALSE;
      }
   }
   
   // matching check, if required
   if (fRejectUnmatched) {
      switch (fDetector) {
         case kITS:
            if (!IsITS()) {
               AliDebug(3, "Rejecting track not matched in ITS");
               return kFALSE;
            }
            break;
         case kTPC:
            if (!IsTPC()) {
               AliDebug(3, "Rejecting track not matched in TPC");
               return kFALSE;
            }
            break;
         case kTOF:
            if (!IsTOF()) {
               AliDebug(3, "Rejecting track not matched in TOF");
               return kFALSE;
            }
            break;
         default:
            AliWarning("Required to reject unmatched tracks, but no detector defined. Track rejected");
            return kFALSE;
      }
   }
   
   // check PID
   // the number of sigmas is set as cut value, which is then checked
   // using the basic functions available in AliRsnCut
   switch (fDetector) {
      case kITS:
         fCutValueD = pid->NumberOfSigmasITS(fDaughter->GetRef(), fSpecies);
         break;
      case kTPC:
         fCutValueD = pid->NumberOfSigmasTPC(fDaughter->GetRef(), fSpecies);
         break;
      case kTOF:
         fCutValueD = pid->NumberOfSigmasTOF(fDaughter->GetRef(), fSpecies);
         break;
      default:
         fCutValueD = 1E20;
         return kFALSE;
   }
   
   return OkRangeD();
}

//_________________________________________________________________________________________________
void AliRsnCutPIDNSigma::Print(const Option_t *) const
{
//
// Print information on this cut
//

   Char_t mom[200], det[10], match[200];
   
   if (fRejectOutside)
      snprintf(mom, 200, "Tracks are accepted only in the momentum range %.2f --> %.2f", fMomMin, fMomMax);
   else
      snprintf(mom, 200, "No check in momentum range");
      
   if (fRejectUnmatched)
      snprintf(match, 200, "Unmatched tracks are rejected");
   else
      snprintf(match, 200, "No check on track matching");
      
   switch (fDetector) {
      case kITS: sprintf(det, "ITS"); break;
      case kTPC: sprintf(det, "TPC"); break;
      case kTOF: sprintf(det, "TOF"); break;
      default  : sprintf(det, "undefined");
   }

   AliInfo(Form("Cut name          : %s", GetName()));
   AliInfo(Form("--> PID detector  : %s", det));
   AliInfo(Form("--> match criteria: %s", match));
   AliInfo(Form("--> momentum range: %s", mom));   
}
