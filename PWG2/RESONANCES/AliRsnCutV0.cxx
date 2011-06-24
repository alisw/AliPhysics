//
// Class AliRsnCutV0
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

#include <Riostream.h>
#include <TFormula.h>
#include <TBits.h>

#include "AliLog.h"
#include "AliESDtrackCuts.h"

#include "AliRsnEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnCutV0.h"

ClassImp(AliRsnCutV0)

//_________________________________________________________________________________________________
AliRsnCutV0::AliRsnCutV0(const char *name, Int_t hypothesis) :
   AliRsnCut(name, AliRsnTarget::kDaughter),
   fHypothesis(0),
   fMass(0.0),
   fTolerance(0.2),
   fMaxDCAVertex(0.3),
   fMinCosPointAngle(0.95),
   fMaxDaughtersDCA(0.1),
   fESDtrackCuts(0x0)
{
//
// Default constructor.
// Initializes all cuts in such a way that all of them are disabled.
//

   SetHypothesis(hypothesis);
}

//_________________________________________________________________________________________________
AliRsnCutV0::AliRsnCutV0(const AliRsnCutV0 &copy) :
   AliRsnCut(copy),
   fHypothesis(copy.fHypothesis),
   fMass(copy.fMass),
   fTolerance(copy.fTolerance),
   fMaxDCAVertex(copy.fMaxDCAVertex),
   fMinCosPointAngle(copy.fMinCosPointAngle),
   fMaxDaughtersDCA(copy.fMaxDaughtersDCA),
   fESDtrackCuts(copy.fESDtrackCuts)
{
//
// Copy constructor.
// Just copy all data member values.
//
}

//_________________________________________________________________________________________________
AliRsnCutV0& AliRsnCutV0::operator=(const AliRsnCutV0 &copy)
{
//
// Assignment operator.
// Just copy all data member values.
//


   fHypothesis = copy.fHypothesis;
   fMass = copy.fMass;
   fTolerance = copy.fTolerance;
   fMaxDCAVertex = copy.fMaxDCAVertex;
   fMinCosPointAngle = copy.fMinCosPointAngle;
   fMaxDaughtersDCA = copy.fMaxDaughtersDCA;
   fESDtrackCuts = copy.fESDtrackCuts;

   return (*this);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutV0::IsSelected(TObject *object)
{
//
// Cut checker.
// Checks the type of object being evaluated
// and then calls the appropriate sub-function (for ESD or AOD)
//

   // coherence check
   if (!TargetOK(object)) return kFALSE;
   
   // check cast
   AliESDv0 *v0esd = fDaughter->Ref2ESDv0();
   AliAODv0 *v0aod = fDaughter->Ref2AODv0();
   //cout << fDaughter->GetRef()->ClassName() << ' ' << v0esd << ' ' << v0aod << endl;
   
   // operate depending on cast
   if (v0esd) {
      return CheckESD(v0esd);
   } else if (v0aod) {
      return CheckAOD(v0aod);
   } else {
      AliDebugClass(1, "Object is not a V0");
      return kFALSE;
   }
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutV0::CheckESD(AliESDv0 *v0)
{
//
// Check an ESD V0.
// This is done using the default track checker for ESD.
// It is declared static, not to recreate it every time.
//

   AliDebugClass(1, "Check ESD");
   if (v0->GetOnFlyStatus()) {
      AliDebugClass(1, "Rejecting V0 in 'on fly' status");
      return kFALSE; // if kTRUE, then this V0 is recontructed
   }
   
   // retrieve pointer to owner event
   AliESDEvent *lESDEvent = fEvent->GetRefESD();
   Double_t xPrimaryVertex = lESDEvent->GetPrimaryVertex()->GetX();
   Double_t yPrimaryVertex = lESDEvent->GetPrimaryVertex()->GetY();
   Double_t zPrimaryVertex = lESDEvent->GetPrimaryVertex()->GetZ();
   AliDebugClass(2, Form("Primary vertex: %f %f %f", xPrimaryVertex, yPrimaryVertex, zPrimaryVertex));
   
   // retrieve the V0 daughters
   UInt_t lIdxPos      = (UInt_t) TMath::Abs(v0->GetPindex());
   UInt_t lIdxNeg      = (UInt_t) TMath::Abs(v0->GetNindex());
   AliESDtrack *pTrack = lESDEvent->GetTrack(lIdxPos);
   AliESDtrack *nTrack = lESDEvent->GetTrack(lIdxNeg);
   
   // check quality cuts
   if (fESDtrackCuts) {
      AliDebugClass(2, "Checking quality cuts");
      if (!fESDtrackCuts->IsSelected(pTrack)) {
         AliDebugClass(2, "Positive daughter failed quality cuts");
         return kFALSE;
      }
      if (!fESDtrackCuts->IsSelected(nTrack)) {
         AliDebugClass(2, "Negative daughter failed quality cuts");
         return kFALSE;
      }
   }
   
   // filter like-sign V0
   //if ((TMath::Abs(pTrack->GetSign()) - TMath::Abs(nTrack->GetSign()) ) < 0.1) {
   //   AliDebugClass(2, "Failed like-sign V0 check");
   //   return kFALSE;
   //}
   
   // check compatibility with expected species hypothesis
   v0->ChangeMassHypothesis(fHypothesis);
   if ((TMath::Abs(v0->GetEffMass() - fMass)) > fTolerance) {
      AliDebugClass(2, "V0 is not in the expected inv mass range");
      return kFALSE;
   }
   
   // topological checks
   if (TMath::Abs(v0->GetD(xPrimaryVertex, yPrimaryVertex, zPrimaryVertex)) > fMaxDCAVertex) {
      AliDebugClass(2, "Failed check on DCA to primary vertes");
      return kFALSE;
   }
   if (TMath::Abs(v0->GetV0CosineOfPointingAngle()) < fMinCosPointAngle) {
      AliDebugClass(2, "Failed check on cosine of pointing angle");
      return kFALSE;
   }
   if (TMath::Abs(v0->GetDcaV0Daughters()) > fMaxDaughtersDCA) {
      AliDebugClass(2, "Failed check on DCA between daughters");
      return kFALSE;
   }
   
   // if we reach this point, all checks were successful
   AliDebugClass(2, "Good V0 (hallelujah)");
   return kTRUE;   
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutV0::CheckAOD(AliAODv0 *v0)
{
//
// Check an AOD V0.
// This is done doing directly all checks, since there is not
// an equivalend checker for AOD tracks
//

   AliWarning("Cuts is not yet implemented for AOD");
   
   return kTRUE;
}

//_________________________________________________________________________________________________
void AliRsnCutV0::Print(const Option_t *) const
{
//
// Print information on this cut
//
}
