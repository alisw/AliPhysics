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
// authors: Massimo Venaruzzo (massimo.venaruzzo@ts.infn.it)
// modified: Enrico Fragiacomo (enrico.fragiacomo@ts.infn.it)
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
AliRsnCutV0::AliRsnCutV0(const char *name, Int_t hypothesis, AliPID::EParticleType pid, AliPID::EParticleType pid2) :
   AliRsnCut(name, AliRsnTarget::kDaughter),
   fHypothesis(0),
   fMass(0.0),
   fTolerance(0.01),
   fMaxDCAVertex(0.3),
   fMinCosPointAngle(0.95),
   fMaxDaughtersDCA(0.5),
   fMinTPCcluster(70),
   fMaxRapidity(0.8),
   fDCARptFormula(""),
   fPID(pid),
   fPID2(pid2),
   fPIDCutProton(0),
   fPIDCutPion(0),
   fESDtrackCuts(0x0),
   fCutQuality(Form("%sDaughtersQuality", name)),
   fAODTestFilterBit(5)
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
   fMinTPCcluster(copy.fMinTPCcluster),
   fMaxRapidity(copy.fMaxRapidity),
   fDCARptFormula(copy.fDCARptFormula),   
   fPID(copy.fPID),
   fPID2(copy.fPID2),
   fPIDCutProton(copy.fPIDCutProton),
   fPIDCutPion(copy.fPIDCutPion),
   fESDtrackCuts(copy.fESDtrackCuts),
   fCutQuality(copy.fCutQuality),
   fAODTestFilterBit(copy.fAODTestFilterBit)
{
//
// Copy constructor.
// Just copy all data member values.:IsSelected: Object is not a V0 (RESONANCES/AliRsnCutV0.cxx:149)

//
   fCutQuality.SetPtRange(0.15, 1E+20);
   fCutQuality.SetEtaRange(-0.8, 0.8);
   fCutQuality.SetDCARPtFormula(fDCARptFormula);
   fCutQuality.SetDCAZmax(2.0);
   fCutQuality.SetSPDminNClusters(1);
   fCutQuality.SetITSminNClusters(0);
   fCutQuality.SetITSmaxChi2(1E+20);
   fCutQuality.SetTPCminNClusters(fMinTPCcluster);
   fCutQuality.SetTPCmaxChi2(4.0);
   fCutQuality.SetRejectKinkDaughters();
   fCutQuality.SetAODTestFilterBit(5);

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
   fMinTPCcluster = copy.fMinTPCcluster;
   fMaxRapidity = copy.fMaxRapidity;
   fDCARptFormula = copy.fDCARptFormula;
   fCutQuality = copy.fCutQuality;
   fPID = copy.fPID;
   fPID2 = copy.fPID2;
   fPIDCutProton = copy.fPIDCutProton;
   fPIDCutPion = copy.fPIDCutPion;
   fESDtrackCuts = copy.fESDtrackCuts;
   fCutQuality = copy.fCutQuality;
   fAODTestFilterBit = copy.fAODTestFilterBit;

   return (*this);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutV0::IsSelected(TObject *object)
{
//:IsSelected: Object is not a V0 (RESONANCES/AliRsnCutV0.cxx:149)

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

   // operate depending on cast:IsSelected: Object is not a V0 (RESONANCES/AliRsnCutV0.cxx:149)

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
   //Double_t pospTPC   = pTrack->GetTPCmomentum();
   //Double_t negpTPC   = nTrack->GetTPCmomentum();
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

      maxTPC = fPIDCutProton;
      maxTPC2 = fPIDCutPion;

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

     
         maxTPC = fPIDCutProton;
         maxTPC2 = fPIDCutPion;

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
Bool_t AliRsnCutV0::CheckAOD(AliAODv0 *v0)
{
//
// Check an AOD V0.
// This is done doing directly all checks, since there is not
// an equivalend checker for AOD tracks
//

   AliDebugClass(2, "Check AOD");
   if (v0->GetOnFlyStatus()) {
      AliDebugClass(2, "Rejecting V0 in 'on fly' status");
      return kFALSE; // if kTRUE, then this V0 is recontructed
   }



   // retrieve pointer to owner event
   AliAODEvent *lAODEvent = fEvent->GetRefAOD();
   Double_t xPrimaryVertex = lAODEvent->GetPrimaryVertex()->GetX();
   Double_t yPrimaryVertex = lAODEvent->GetPrimaryVertex()->GetY();
   Double_t zPrimaryVertex = lAODEvent->GetPrimaryVertex()->GetZ();
   AliDebugClass(2, Form("Primary vertex: %f %f %f", xPrimaryVertex, yPrimaryVertex, zPrimaryVertex));

   // retrieve the V0 daughters
   AliAODTrack *pTrack = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(0));
   AliAODTrack *nTrack = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(1));


   // check quality cuts
   UInt_t  filtermapP = 9999;
   UInt_t  filtermapN = 9999;
   filtermapP = pTrack->GetFilterMap();
   filtermapN = nTrack->GetFilterMap();


   if ( !pTrack->TestFilterBit(fAODTestFilterBit)   ) {
      AliDebugClass(2, Form("Positive daughter failed quality cuts filtermapP=%d",filtermapP));
      return kFALSE;
   }
   if ( !nTrack->TestFilterBit(fAODTestFilterBit)   ) {
      AliDebugClass(2, Form("Negative daughter failed quality cuts filtermapN=%d",filtermapN));
      return kFALSE;
   }



   /*AliDebugClass(1, Form("fESDtrackCuts=%p",fESDtrackCuts));
   if (fESDtrackCuts) { // use fESDtrackCuts to retrieve cuts values


      AliDebugClass(2, "Checking quality cuts");



      const Bool_t  CutAcceptKinkDaughters  = fESDtrackCuts->GetAcceptKinkDaughters(); // 0 = kFalse
      const Float_t CutMaxChi2PerClusterTPC = fESDtrackCuts->GetMaxChi2PerClusterTPC();
      const Int_t   CutMinNClusterTPC       = fESDtrackCuts->GetMinNClusterTPC();
      const Bool_t  CutRequireTPCRefit      = fESDtrackCuts->GetRequireTPCRefit();
      //AliDebugClass(2, Form("accept kink=%d maxchi2=%f minnclSTPC=%d requireTPCrefit=%d", CutAcceptKinkDaughters,
      //        CutMaxChi2PerClusterTPC, CutMinNClusterTPC, CutRequireTPCRefit));
      //AliDebugClass(2, Form("pTrack->TPCNcls=%d", pTrack->GetTPCNcls() ));
      //AliDebugClass(2, Form("nTrack->TPCNcls=%d", nTrack->GetTPCNcls() ));


      if(pTrack->GetTPCNcls() < CutMinNClusterTPC) {AliDebugClass(2, "Positive daughter not MinNclsTPC"); return kFALSE;}
      if(nTrack->GetTPCNcls() < CutMinNClusterTPC) {AliDebugClass(2, "Negative daughter not MinNclsTPC"); return kFALSE;}

   }*/

   // filter like-sign V0
   if ( TMath::Abs( ((pTrack->Charge()) - (nTrack->Charge())) ) < 0.1) {
      AliDebugClass(2, "Failed like-sign V0 check");
      return kFALSE;
   }

   // check compatibility with expected species hypothesis
   Double_t mass = 0.0;
   if(fHypothesis==kLambda0) {
      mass = v0->MassLambda();
   }
   else if (fHypothesis==kLambda0Bar) {
      mass = v0->MassAntiLambda();
   }
   if ((TMath::Abs(mass - fMass)) > fTolerance) {
      AliDebugClass(2, Form("V0 is not in the expected inv mass range  Mass: %d %f %f", fHypothesis, fMass, mass));
      return kFALSE;
   }
   AliDebugClass(2, Form("Mass: %d %f %f", fHypothesis, fMass, mass));


   // topological checks
   if (TMath::Abs(v0->DcaV0ToPrimVertex()) > fMaxDCAVertex) {
      AliDebugClass(2, Form("Failed check on DCA to primary vertes dca=%f maxdca=%f",TMath::Abs(v0->DcaV0ToPrimVertex()),fMaxDCAVertex));
      return kFALSE;
   }
   if (TMath::Abs( TMath::Cos(v0->OpenAngleV0()) ) < fMinCosPointAngle) {
      AliDebugClass(2, "Failed check on cosine of pointing angle");
      return kFALSE;
   }
   if (TMath::Abs(v0->DcaV0Daughters()) > fMaxDaughtersDCA) {
      AliDebugClass(2, "Failed check on DCA between daughters");
      return kFALSE;
   }
   if (TMath::Abs(v0->RapLambda()) > fMaxRapidity) {
      AliDebugClass(2, "Failed check on V0 rapidity");
      return kFALSE;
   }





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
   //Double_t pospTPC   = pTrack->GetTPCmomentum();
   //Double_t negpTPC   = nTrack->GetTPCmomentum();
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
     
     
         maxTPC = fPIDCutProton; 
         maxTPC2 = fPIDCutPion; 

      if (! ((posnsTPC <= maxTPC) && (negnsTPC2 <= maxTPC2)) ) {
         AliDebugClass(2, "Failed check on V0 PID");
         return kFALSE;
      }
   }

   if(fHypothesis==kLambda0Bar) {

      //if (isTOFneg) {
      // TPC: 5sigma cut for all
      //if (negnsTPC > 5.0) return kFALSE;
      // TOF: 3sigma
      // maxTOF = 3.0;
      //return (negnsTOF <= maxTOF);
      //} else {
      // TPC:

      
         maxTPC = fPIDCutProton;
         maxTPC2 = fPIDCutPion;

      if(! ((negnsTPC <= maxTPC) && (posnsTPC2 <= maxTPC2)) ) {
         AliDebugClass(2, "Failed check on V0 PID");
         return kFALSE;
      }
   }



   // if we reach this point, all checks were successful
   AliDebugClass(1, "Good AOD V0 (hallelujah)");
   AliDebugClass(1, Form("Mass: %d %f %f %d %d", fHypothesis, fMass, mass, filtermapP, filtermapN));
   return kTRUE;

}

//_________________________________________________________________________________________________
void AliRsnCutV0::Print(const Option_t *) const
{
//
// Print information on this cut
//
}
