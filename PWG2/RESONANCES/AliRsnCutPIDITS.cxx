//
// Class AliRsnCutPIDITS
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

#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"

#include "AliRsnCutPIDITS.h"

ClassImp(AliRsnCutPIDITS)

//_________________________________________________________________________________________________
AliRsnCutPIDITS::AliRsnCutPIDITS
(const char *name, AliPID::EParticleType ref, Double_t min, Double_t max, Bool_t rejectOutside) :
   AliRsnCut(name, AliRsnCut::kDaughter, min, max),
   fIsMC(kFALSE),
   fRejectOutside(rejectOutside),
   fMomMin(0.0),
   fMomMax(1E+20),
   fRefType(ref),
   fESDpid(),
   fAODpid()
{
//
// Main constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutPIDITS::AliRsnCutPIDITS
(const AliRsnCutPIDITS& copy) :
   AliRsnCut(copy),
   fIsMC(copy.fIsMC),
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
}

//_________________________________________________________________________________________________
AliRsnCutPIDITS& AliRsnCutPIDITS::operator=(const AliRsnCutPIDITS& copy)
{
//
// Assignment operator
//

   AliRsnCut::operator=(copy);

   fIsMC          = copy.fIsMC;
   fRejectOutside = copy.fRejectOutside;
   fMomMin        = copy.fMomMin;
   fMomMax        = copy.fMomMax;
   fRefType       = copy.fRefType;
   fESDpid        = copy.fESDpid;
   fAODpid        = copy.fAODpid;

   return (*this);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPIDITS::IsSelected(TObject *object)
{
//
// Cut checker.
//

   // coherence check
   if (!TargetOK(object)) return kFALSE;

   // reject not ITS tracks
   // status is checked in the same way for all tracks
   AliVTrack *vtrack = fDaughter->Ref2Vtrack();
   if (!vtrack) {
      AliDebug(AliLog::kDebug + 2, Form("Impossible to process an object of type '%s'. Cut applicable only to ESD/AOD tracks", fDaughter->GetRef()->ClassName()));
      return kFALSE;
   }

   // check status, to know it track is an ITS+TPC or ITS standalone
   // and reject it if it is of none of them
   Bool_t isSA = kFALSE;
   if (IsTPC(vtrack)) isSA = kFALSE;
   else if (IsSA(vtrack)) isSA = kTRUE;
   else {
      AliDebug(AliLog::kDebug + 2, "Status flags unmatched");
      return kFALSE;
   }

   // common evaluation variables
   Int_t        k, nITSpidLayers = 0;
   Double_t     mom      = vtrack->P();
   AliESDtrack *esdTrack = fDaughter->Ref2ESDtrack();
   AliAODTrack *aodTrack = fDaughter->Ref2AODtrack();

   // count number of PID layers...
   if (esdTrack) {
      UChar_t itsCluMap = esdTrack->GetITSClusterMap();
      for (k = 2; k < 6; k++) if (itsCluMap & (1 << k)) ++nITSpidLayers;
   } else if (aodTrack) {
      for (k = 2; k < 6; k++) if (TESTBIT(aodTrack->GetITSClusterMap(), k)) ++nITSpidLayers;
   } else {
      AliDebug(AliLog::kDebug + 2, Form("Impossible to process an object of type '%s'. Cut applicable only to ESD/AOD tracks", fDaughter->GetRef()->ClassName()));
      return kFALSE;
   }
   // ...and reject tracks where it is smaller than 3
   if (nITSpidLayers < 3) {
      AliDebug(AliLog::kDebug + 2, "Rejecting track with too few ITS pid layers");
      return kFALSE;
   }

   // assign PID nsigmas to default cut check value
   // since bad object types are rejected before, here we have an ESD track or AOD track
   if (esdTrack) {
      if (!fESDpid) fESDpid = new AliESDpid(fIsMC);
      fCutValueD = fESDpid->GetITSResponse().GetNumberOfSigmas(mom, esdTrack->GetITSsignal(), fRefType, nITSpidLayers, isSA);
   } else {
      if (!fAODpid) fAODpid = new AliAODpidUtil(fIsMC);
      fCutValueD = fAODpid->NumberOfSigmasITS(aodTrack, fRefType);
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
void AliRsnCutPIDITS::Print(const Option_t *) const
{
//
// Print information on this cut
//

   AliInfo(Form("Cut name                    : %s", GetName()));
   AliInfo(Form("--> cut range (nsigma)      : %.3f %.3f", fMinD, fMaxD));
   AliInfo(Form("--> momentum range          : %.3f %.3f", fMomMin, fMomMax));
   AliInfo(Form("--> tracks outside range are: %s", (fRejectOutside ? "rejected" : "accepted")));
}
