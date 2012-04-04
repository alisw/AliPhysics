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
//AliRsnCutV0::AliRsnCutV0(const char *name, Int_t hypothesis) :
AliRsnCutV0::AliRsnCutV0(const char *name, Int_t hypothesis, AliPID::EParticleType pid, AliPID::EParticleType pid2) :
   AliRsnCut(name, AliRsnTarget::kDaughter),
   fHypothesis(0),
   fMass(0.0),
   fTolerance(0.01),
   fMaxDCAVertex(0.3),
   fMinCosPointAngle(0.95),
   fMaxDaughtersDCA(0.5),
   fMaxRapidity(0.8),
   fPID(pid),
   fPID2(pid2),
   fPIDCut1(4.0),
   fPIDCut2(3.0),
   fPIDCut3(3.0),
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
   fMaxRapidity(copy.fMaxRapidity),
   fPID(copy.fPID),
   fPID2(copy.fPID2),
   fPIDCut1(copy.fPIDCut1),
   fPIDCut2(copy.fPIDCut2),
   fPIDCut3(copy.fPIDCut3),
   fESDtrackCuts(copy.fESDtrackCuts)
{
//
// Copy constructor.
// Just copy all data member values.
//
}

//_________________________________________________________________________________________________
AliRsnCutV0 &AliRsnCutV0::operator=(const AliRsnCutV0 &copy)
{
//
// Assignment operator.
// Just copy all data member values.
//

   if (this == &copy)
      return *this;
   fHypothesis = copy.fHypothesis;
   fMass = copy.fMass;
   fTolerance = copy.fTolerance;
   fMaxDCAVertex = copy.fMaxDCAVertex;
   fMinCosPointAngle = copy.fMinCosPointAngle;
   fMaxDaughtersDCA = copy.fMaxDaughtersDCA;
   fMaxRapidity = copy.fMaxRapidity;
   fPID = copy.fPID;
   fPID2 = copy.fPID2;
   fPIDCut1 = copy.fPIDCut1;
   fPIDCut2 = copy.fPIDCut2;
   fPIDCut3 = copy.fPIDCut3;
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
   if ( TMath::Abs( ((pTrack->GetSign()) - (nTrack->GetSign())) ) < 0.1) {
      AliDebugClass(2, "Failed like-sign V0 check");
      return kFALSE;
   }


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
   if (TMath::Abs(v0->Y(fHypothesis)) > fMaxRapidity) {
      AliDebugClass(2, "Failed check on V0 rapidity");
      return kFALSE;
   }


   // check PID on proton or antiproton from V0

   // check initialization of PID object
   AliPIDResponse *pid = fEvent->GetPIDResponse();
   if (!pid) {
      AliFatal("NULL PID response");
      return kFALSE;
   }

   // check if TOF is matched
   // and computes all values used in the PID cut
   //Bool_t   isTOFpos  = MatchTOF(ptrack);
   //Bool_t   isTOFneg  = MatchTOF(ntrack);
   Double_t pospTPC   = pTrack->GetTPCmomentum();
   Double_t negpTPC   = nTrack->GetTPCmomentum();
   //Double_t posp      = pTrack->P();
   //Double_t negp      = nTrack->P();
   Double_t posnsTPC   = TMath::Abs(pid->NumberOfSigmasTPC(pTrack, fPID));
   Double_t posnsTPC2  = TMath::Abs(pid->NumberOfSigmasTPC(pTrack, fPID2));
   //Double_t posnsTOF  = TMath::Abs(pid->NumberOfSigmasTOF(ptrack, fPID));
   Double_t negnsTPC   = TMath::Abs(pid->NumberOfSigmasTPC(nTrack, fPID));
   Double_t negnsTPC2  = TMath::Abs(pid->NumberOfSigmasTPC(nTrack, fPID2));
   //Double_t negnsTOF  = TMath::Abs(pid->NumberOfSigmasTOF(ntrack, fPID));
   Double_t maxTPC = 1E20;
   Double_t maxTPC2 = 1E20;
   //Double_t maxTOF = 1E20;

   // applies the cut differently depending on the PID and the momentum

   if(fHypothesis==kLambda0) {
      //if (isTOFpos) {
      // TPC: 5sigma cut for all
      //if (posnsTPC > 5.0) return kFALSE;
      // TOF: 3sigma
      // maxTOF = 3.0;
      //return (posnsTOF <= maxTOF);
      //} else {
      // TPC:
      // below 600 MeV: 4sigma
      // above 600 MeV: 3sigma

      if (pospTPC <= 0.6 && fPID==AliPID::kProton)
         //maxTPC = 4.0;
         maxTPC = fPIDCut1;
      else if (pospTPC > 0.6 && fPID==AliPID::kProton)
         //maxTPC = 3.0;
         maxTPC = fPIDCut2;
      //else
      //return kFALSE;

      //maxTPC2 = 3.0;
      maxTPC2 = fPIDCut3;

      if (! ((posnsTPC <= maxTPC) && (negnsTPC2 <= maxTPC2)) ) {
         AliDebugClass(2, "Failed check on V0 PID");
         return kFALSE;
      }
   }

   //}

   if(fHypothesis==kLambda0Bar) {
      //if (isTOFneg) {
      // TPC: 5sigma cut for all
      //if (negnsTPC > 5.0) return kFALSE;
      // TOF: 3sigma
      // maxTOF = 3.0;
      //return (negnsTOF <= maxTOF);
      //} else {
      // TPC:
      // below 600 MeV: 4sigma
      // above 600 MeV: 3sigma

      if (negpTPC <= 0.6 && fPID==AliPID::kProton)
         //maxTPC = 4.0;
         maxTPC = fPIDCut1;
      else if (negpTPC > 0.6 && fPID==AliPID::kProton)
         //maxTPC = 3.0;
         maxTPC = fPIDCut2;
      else
         return kFALSE;

      //maxTPC2 = 3.0;
      maxTPC2 = fPIDCut3;

      if(! ((negnsTPC <= maxTPC) && (posnsTPC2 <= maxTPC2)) ) {
         AliDebugClass(2, "Failed check on V0 PID");
         return kFALSE;
      }
   }
   //}


   // if we reach this point, all checks were successful
   AliDebugClass(2, "Good V0 (hallelujah)");
   return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutV0::CheckAOD(AliAODv0 *)
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
