//
// Class AliRsnCutPIDTOF
//
// Implements the PID check with TOF detector,
// computed as a compatibility within a given range
// expressed in number of sigmas w.r. to expected time.
// Uses the default cut checking facilities of AliRsnCut
// to check this.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliAODTrack.h"
#include "AliAODpidUtil.h"

#include "AliRsnCutPIDTOF.h"

ClassImp(AliRsnCutPIDTOF)

//_________________________________________________________________________________________________
AliRsnCutPIDTOF::AliRsnCutPIDTOF
(const char *name, AliPID::EParticleType ref, Double_t min, Double_t max, Bool_t rejectUnmatched) :
   AliRsnCut(name, AliRsnCut::kDaughter, min, max),
   fRejectUnmatched(rejectUnmatched),
   fRefType(AliPID::kUnknown),
   fRefMass(0.0),
   fESDpid(0x0),
   fAODpid(0x0)
{
//
// Default constructor.
// To set the reference PID type, calls the SetRefType method,
// which sets the mass accordingly and coherently.
//

   SetRefType(ref);
}

//_________________________________________________________________________________________________
AliRsnCutPIDTOF::AliRsnCutPIDTOF(const AliRsnCutPIDTOF& copy) :
   AliRsnCut(copy),
   fRejectUnmatched(copy.fRejectUnmatched),
   fRefType(AliPID::kUnknown),
   fRefMass(0.0),
   fESDpid(copy.fESDpid),
   fAODpid(copy.fAODpid)
{
//
// Copy constructor.
// To set the reference PID type, calls the SetRefType method,
// which sets the mass accordingly and coherently.
//

   SetRefType(copy.fRefType);
}

//_________________________________________________________________________________________________
AliRsnCutPIDTOF& AliRsnCutPIDTOF::operator=(const AliRsnCutPIDTOF& copy)
{
//
// Assignment operator.
// To set the reference PID type, calls the SetRefType method,
// which sets the mass accordingly and coherently.
//

   fRejectUnmatched = copy.fRejectUnmatched;
   fESDpid          = copy.fESDpid;
   fAODpid          = copy.fAODpid;

   SetRefType(copy.fRefType);

   return (*this);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPIDTOF::IsSelected(TObject *object)
{
//
// Cut checker.
//

   // coherence check
   if (!TargetOK(object)) return kFALSE;

   // reject always non-track objects
   AliVTrack *vtrack = fDaughter->Ref2Vtrack();
   if (!vtrack) {
      AliDebug(AliLog::kDebug + 2, Form("Impossible to process an object of type '%s'. Cut applicable only to ESD/AOD tracks", fDaughter->GetRef()->ClassName()));
      return kFALSE;
   }

   // checks that track is matched in TOF:
   // if not, the track is accepted or rejected
   // depending on the 'fRejectUnmatched' data member:
   // -- kTRUE  --> all unmatched tracks are rejected
   // -- kFALSE --> all unmatched tracks are accepted (it is assumed that other PIDs are done)
   if (!IsMatched(vtrack)) {
      AliDebug(AliLog::kDebug + 2, "Track is not matched with TOF");
      if (fRejectUnmatched) return kFALSE;
   }

   // retrieve real object type and
   // prepare some useful variables
   Double_t     tof, sigma, times[5];
   Double_t    &ref = times[(Int_t)fRefType];
   AliESDtrack *esdTrack = fDaughter->Ref2ESDtrack();
   AliAODTrack *aodTrack = fDaughter->Ref2AODtrack();

   // cut check depends on the object type
   if (esdTrack) {
      // setup the ESD PID object
      AliESDEvent *esd = 0x0;
      if (fEvent) esd = fEvent->GetRefESD();
      if (!esd) {
         AliError("Processing an ESD track, but target is not an ESD event");
         return kFALSE;
      }
      if (!fESDpid) fESDpid = new AliESDpid;
      fESDpid->SetTOFResponse(esd, AliESDpid::kTOF_T0);

      // get time of flight, reference times and sigma
      esdTrack->GetIntegratedTimes(times);
      tof   = (Double_t)(esdTrack->GetTOFsignal() - fESDpid->GetTOFResponse().GetStartTime(esdTrack->P()));
      sigma = (Double_t)fESDpid->GetTOFResponse().GetExpectedSigma(esdTrack->P(), ref, fRefMass);

      // port values to standard AliRsnCut checker
      fCutValueD = (tof - ref) / sigma;
      return OkRangeD();
   } else if (aodTrack) {
      // for AOD tracks, all operations are done by the AOD PID utility
      if (!fAODpid) fAODpid = new AliAODpidUtil;
      fCutValueD = (Double_t)fAODpid->NumberOfSigmasTOF(aodTrack, fRefType);
      return OkRangeD();
   } else {
      AliDebug(AliLog::kDebug + 2, Form("Impossible to process an object of type '%s'. Cut applicable only to ESD/AOD tracks", fDaughter->GetRef()->ClassName()));
      return kFALSE;
   }
}

//_________________________________________________________________________________________________
void AliRsnCutPIDTOF::Print(const Option_t *) const
{
//
// Print information on this cut
//

   AliInfo(Form("Cut name, type            : %s %s", GetName(), ClassName()));
   AliInfo(Form("TOF PID cut range (sigmas): %.3f %.3f", fMinD, fMaxD));
   AliInfo(Form("Unmatched tracks are      : %s", (fRejectUnmatched ? "rejected" : "accepted")));
}

