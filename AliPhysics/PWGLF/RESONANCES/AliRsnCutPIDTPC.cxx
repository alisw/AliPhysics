//
// Class AliRsnCutPIDTPC
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

#include "AliPID.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliAODpidUtil.h"

#include "AliRsnCutPIDTPC.h"

ClassImp(AliRsnCutPIDTPC)

//_________________________________________________________________________________________________
AliRsnCutPIDTPC::AliRsnCutPIDTPC
(const char *name, AliPID::EParticleType type, Double_t min, Double_t max, Bool_t rejectOutside) :
   AliRsnCut(name, AliRsnCut::kDaughter, min, max),
   fRejectOutside(rejectOutside),
   fMomMin(0.0),
   fMomMax(1E+20),
   fRefType(type),
   fESDpid(0x0),
   fAODpid(0x0)
{
//
// Main constructor.
//

   fBB[0] = fBB[1] = fBB[2] = fBB[3] = fBB[4] = 0.0;
}

//_________________________________________________________________________________________________
AliRsnCutPIDTPC::AliRsnCutPIDTPC
(const AliRsnCutPIDTPC &copy) :
   AliRsnCut(copy),
   fRejectOutside(copy.fRejectOutside),
   fMomMin(copy.fMomMin),
   fMomMax(copy.fMomMax),
   fRefType(copy.fRefType),
   fESDpid(copy.fESDpid),
   fAODpid(copy.fAODpid)
{
//
// Copy constructor.
//

   Int_t i;
   for (i = 0; i < 5; i++) fBB[i] = copy.fBB[i];
}

//_________________________________________________________________________________________________
AliRsnCutPIDTPC &AliRsnCutPIDTPC::operator=(const AliRsnCutPIDTPC &copy)
{
//
// Assignment operator
//

   AliRsnCut::operator=(copy);
   if (this == &copy)
      return *this;

   fRejectOutside = copy.fRejectOutside;
   fMomMin        = copy.fMomMin;
   fMomMax        = copy.fMomMax;
   fRefType       = copy.fRefType;
   fESDpid        = copy.fESDpid;
   fAODpid        = copy.fAODpid;

   Int_t i;
   for (i = 0; i < 5; i++) fBB[i] = copy.fBB[i];

   return (*this);
}

//_________________________________________________________________________________________________
void AliRsnCutPIDTPC::SetBBParam(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4)
{
//
// Properly set the Bethe-Bloch parameters in all places where it is needed.
//

   fBB[0] = p0;
   fBB[1] = p1;
   fBB[2] = p2;
   fBB[3] = p3;
   fBB[4] = p4;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPIDTPC::IsSelected(TObject *object)
{
//
// Cut checker.
//

   // coherence check
   if (!TargetOK(object)) return kFALSE;

   // common evaluation variables
   Double_t     mom;
   AliESDtrack *esdTrack = fDaughter->Ref2ESDtrack();
   AliAODTrack *aodTrack = fDaughter->Ref2AODtrack();

   // get inner momentum, needed for BB computation
   if (esdTrack) {
      if (!esdTrack->GetInnerParam()) {
         AliDebug(AliLog::kDebug + 2, "No inner param");
         return kFALSE;
      }
      mom = esdTrack->GetInnerParam()->P();
   } else if (aodTrack) {
      if (!aodTrack->GetDetPid()) {
         AliDebug(AliLog::kDebug + 2, "No def-pid object");
         return kFALSE;
      }
      mom = aodTrack->GetDetPid()->GetTPCmomentum();
      if (mom < 1E-6) return kFALSE;
   } else {
      AliDebug(AliLog::kDebug + 2, Form("Impossible to process an object of type '%s'. Cut applicable only to ESD/AOD tracks", fDaughter->GetRef()->ClassName()));
      return kFALSE;
   }

   // assign PID nsigmas to default cut check value
   // since bad object types are rejected before, here we have an ESD track or AOD track
   if (esdTrack) {
      if (!fESDpid) {
         fESDpid = new AliESDpid;
         fESDpid->GetTPCResponse().SetBetheBlochParameters(fBB[0], fBB[1], fBB[2], fBB[3], fBB[4]);
      }
      fCutValueD = fESDpid->GetTPCResponse().GetNumberOfSigmas(mom, esdTrack->GetTPCsignal(), esdTrack->GetTPCsignalN(), fRefType);
   } else {
      if (!fAODpid) {
         fAODpid = new AliAODpidUtil;
         fAODpid->GetTPCResponse().SetBetheBlochParameters(fBB[0], fBB[1], fBB[2], fBB[3], fBB[4]);
      }
      if (aodTrack->GetTPCsignalN() == 0) aodTrack->GetDetPid()->SetTPCsignalN(aodTrack->GetTPCNcls());
      fCutValueD = fAODpid->NumberOfSigmasTPC(aodTrack, fRefType);
   }

   // use AliRsnCut default method to check cut
   Bool_t cutCheck = OkRangeD();

   // now check the momentum:
   // -- if it stays inside the accepted range, track just checked
   //    with respect to the nsigma band
   // -- if it stays outside the accepted range and 'fRejectOutside' is kTRUE,
   //    track is always rejected, while if 'fRejectOutside' is kFALSE,
   //    track is accepted if it stays inside the nsigma band
   if ((mom >= fMomMin && mom <= fMomMax))
      return cutCheck;
   else {
      AliDebug(AliLog::kDebug + 2, Form("Track momentum = %.5f, outside allowed range", mom));
      return ((!fRejectOutside) && cutCheck);
   }
}

//_________________________________________________________________________________________________
void AliRsnCutPIDTPC::Print(const Option_t *) const
{
//
// Print information on this cut
//

   AliInfo(Form("Cut name                    : %s", GetName()));
   AliInfo(Form("--> cut range (nsigma)      : %.3f %.3f", fMinD, fMaxD));
   AliInfo(Form("--> momentum range          : %.3f %.3f", fMomMin, fMomMax));
   AliInfo(Form("--> tracks outside range are: %s", (fRejectOutside ? "rejected" : "accepted")));
}
