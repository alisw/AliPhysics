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
#include "AliESDpid.h"
#include "AliAODpidUtil.h"

#include "AliRsnCutPIDNSigma.h"

ClassImp(AliRsnCutPIDNSigma)
ClassImp(AliRsnCutPIDNSigma::AliRsnPIDRange)

//_________________________________________________________________________________________________
AliRsnCutPIDNSigma::AliRsnCutPIDNSigma() : 
   AliRsnCut("cut", AliRsnTarget::kDaughter),
   fSpecies(AliPID::kUnknown),
   fDetector(kDetectors),
   fRejectUnmatched(kFALSE),
   fTrackNSigma(0.0),
   fTrackMom(0.0),
   fMyPID(0x0),
   fRanges("AliRsnCutPIDNSigma::AliRsnPIDRange", 0)
{
//
// Main constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutPIDNSigma::AliRsnCutPIDNSigma
(const char *name, AliPID::EParticleType species, EDetector det) :
   AliRsnCut(name, AliRsnTarget::kDaughter),
   fSpecies(species),
   fDetector(det),
   fRejectUnmatched(kFALSE),
   fTrackNSigma(0.0),
   fTrackMom(0.0),
   fMyPID(0x0),
   fRanges("AliRsnCutPIDNSigma::AliRsnPIDRange", 0)
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
   fTrackNSigma(0.0),
   fTrackMom(0.0),
   fMyPID(copy.fMyPID),
   fRanges(copy.fRanges)
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
   fMyPID = copy.fMyPID;
   fRanges = copy.fRanges;

   return (*this);
}

//__________________________________________________________________________________________________
void AliRsnCutPIDNSigma::InitMyPID(Bool_t isMC, Bool_t isESD)
{
//
// Initialize manual PID object
//

   if (isESD) 
      fMyPID = new AliESDpid(isMC);
   else
      fMyPID = new AliAODpidUtil(isMC);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPIDNSigma::IsSelected(TObject *object)
{
//
// Cut checker.
// As usual, there are 'kFALSE' exit points whenever one of the conditions is not passed,
// and at the end, it returns kTRUE since it bypassed all possible exit points.
//

   // coherence check
   if (!TargetOK(object)) return kFALSE;
   
   // check initialization of PID object
   // if manual PID is used, use that, otherwise get from source event
   AliPIDResponse *pid = 0x0;
   if (fMyPID) 
      pid = fMyPID; 
   else 
      pid = fEvent->GetPIDResponse();
   if (!pid) {
      AliFatal("NULL PID response");
      return kFALSE;
   }

   // convert input object into AliVTrack
   // if this fails, the cut cannot be checked
   AliVTrack *vtrack = fDaughter->Ref2Vtrack();
   if (!vtrack) {
      AliDebugClass(2, "Referenced daughter is not a track");
      return kFALSE;
   }
   
   // check matching, if required
   // a) if detector is not matched and matching is required, reject the track
   // b) if detector is not matched and matching is not required, accept blindly the track
   //    since missing the matching causes one not to be able to rely that detector
   if (!MatchDetector(vtrack)) {
      AliDebugClass(2, Form("Detector not matched. fRejectUnmatched = %s --> track is %s", (fRejectUnmatched ? "true" : "false"), (fRejectUnmatched ? "rejected" : "accepted")));
      return (!fRejectUnmatched);
   }
   
   // get reference momentum
   fTrackMom = (fDetector == kTPC) ? vtrack->GetTPCmomentum() : vtrack->P();
   
   // get number of sigmas
   switch (fDetector) {
      case kITS:
         fTrackNSigma = TMath::Abs(pid->NumberOfSigmasITS(vtrack, fSpecies));
         break;
      case kTPC:
         fTrackNSigma = TMath::Abs(pid->NumberOfSigmasTPC(vtrack, fSpecies));
         break;
      case kTOF:
         fTrackNSigma = TMath::Abs(pid->NumberOfSigmasTOF(vtrack, fSpecies));
         break;
      default:
         AliError("Bad detector chosen. Rejecting track");
         return kFALSE;
   }
   
   // loop on all ranges, and use the one which contains this momentum
   // if none is found, the cut is not passed
   Bool_t accept = kFALSE;
   Int_t  i, goodRange = -1, nRanges = fRanges.GetEntriesFast();
   for (i = 0; i < nRanges; i++) {
      AliRsnPIDRange *range = (AliRsnPIDRange*)fRanges[i];
      if (!range) continue;
      if (!range->IsInRange(fTrackMom)) continue;
      else {
         goodRange = i;
         accept = range->CutPass(fTrackNSigma);
         AliDebugClass(2, Form("[%s] NSigma = %.3f, max = %.3f, track %s", GetName(), fTrackNSigma, range->NSigmaCut(), (accept ? "accepted" : "rejected")));
         break;
      }
   }
   if (goodRange < 0) {
      AliDebugClass(2, Form("[%s] No good range found. Rejecting track", GetName()));
      return kFALSE;
   } else {
      AliDebugClass(2, Form("[%s] Mom = %.3f, good range found (#%d), track was %s", GetName(), fTrackMom, goodRange, (accept ? "accepted" : "rejected")));
      return accept;
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
