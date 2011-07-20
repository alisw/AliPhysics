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
AliRsnCutPIDNSigma::AliRsnCutPIDNSigma() : 
   AliRsnCut("cut", AliRsnTarget::kDaughter),
   fSpecies(AliPID::kUnknown),
   fDetector(kDetectors),
   fRejectUnmatched(kFALSE),
   fMomMin(0.0),
   fMomMax(1E20),
   fNSigma(1E20)
{
//
// Main constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutPIDNSigma::AliRsnCutPIDNSigma
(const char *name, AliPID::EParticleType species, EDetector det, Double_t nsigma) :
   AliRsnCut(name, AliRsnTarget::kDaughter),
   fSpecies(species),
   fDetector(det),
   fRejectUnmatched(kFALSE),
   fMomMin(0.0),
   fMomMax(1E20),
   fNSigma(nsigma)
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
   fRejectUnmatched(copy.fRejectUnmatched),
   fMomMin(copy.fMomMin),
   fMomMax(copy.fMomMax),
   fNSigma(copy.fNSigma)
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
   fRejectUnmatched = copy.fRejectUnmatched;
   fMomMin = copy.fMomMin;
   fMomMax = copy.fMomMax;
   fNSigma = copy.fNSigma;

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
   AliVTrack *vtrack = fDaughter->Ref2Vtrack();
   if (!vtrack) {
      AliDebugClass(2, "Referenced daughter is not a track");
      return kFALSE;
   }
   if (fDetector == kTPC)
      momentum = vtrack->GetTPCmomentum();
   else
      momentum = vtrack->P();
      
   // check momentum range
   if (momentum < fMomMin || momentum > fMomMax) {
      AliDebugClass(2, Form("Track momentum = %.5f, outside allowed range [%.2f - %.2f]", momentum, fMomMin, fMomMax));
      return kFALSE;
   }
   
   // check PID
   Bool_t matched;
   Double_t nsigma;
   switch (fDetector) {
      case kITS:
         matched = IsITS(vtrack);
         nsigma  = pid->NumberOfSigmasITS(vtrack, fSpecies);
         break;
      case kTPC:
         matched = IsTPC(vtrack);
         nsigma  = pid->NumberOfSigmasTPC(vtrack, fSpecies);
         break;
      case kTOF:
         matched = IsTOF(vtrack);
         nsigma  = pid->NumberOfSigmasTOF(vtrack, fSpecies);
         break;
      default:
         AliError("Bad detector chosen. Rejecting track");
         return kFALSE;
   }
   
   // determine cut result
   if (fRejectUnmatched && (matched == kFALSE)) {
      AliDebugClass(2, "Required to reject unmatched traks, and this track is not matched in the detector");
      return kFALSE;
   } else {
      AliDebugClass(2, Form("Nsigma = %.5f, maximum allowed = %.2f", nsigma, fNSigma));
      if (TMath::Abs(nsigma) <= fNSigma) {
         AliDebugClass(2, "Track accepted");
         return kTRUE;
      } else {
         AliDebugClass(2, "Track rejected");
         return kFALSE;
      }
   }
}

//_________________________________________________________________________________________________
void AliRsnCutPIDNSigma::Print(const Option_t *) const
{
//
// Print information on this cut
//

   Char_t mom[200], det[100], match[200];
      
   if (fRejectUnmatched)
      snprintf(match, 200, "Unmatched tracks are rejected");
   else
      snprintf(match, 200, "No check on track matching");
      
   switch (fDetector) {
      case kITS: snprintf(det, 3, "ITS"); break;
      case kTPC: snprintf(det, 3, "TPC"); break;
      case kTOF: snprintf(det, 3, "TOF"); break;
      default  : snprintf(det, 3, "undefined");
   }

   AliInfo(Form("Cut name          : %s", GetName()));
   AliInfo(Form("--> PID detector  : %s", det));
   AliInfo(Form("--> match criteria: %s", match));
   AliInfo(Form("--> momentum range: %s", mom));   
}
