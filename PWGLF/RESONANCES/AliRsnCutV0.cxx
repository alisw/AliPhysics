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

// modified: Kunal Garg (kgarg@cern.ch)
//  Modifications: Added Competing V0 Rejection, Lifetime cut and a switch for using or not using Competing V0 Rejection
//Modification (13 March 2018): Added pT dependent Mass Tolerance cut for K* charged

/* Note: Competing V0 Rejection
For selection of V0 particle, we typically set a wide range around the mass of the particle. With Competing V0 rejection, we reject the other V0 particle in the same region with a small range.
Lifetime cut should be quite intuitive. Set the value to a high value if not needed in the analysis.
*/

/*Switch to a pT dependent Mass Tolerance Cut has been added. fpT_Tolerance is a flag to switch to this pT dependent version of the cut. fMassTolSigma is as it says, how strict you want the cut to be.
 fpT_Tolerance takes values from 0-3.
 The values go as follows:
 fpT_Tolerance = 0 ; Default. Fixed cut for mass tolerance
 fpT_Tolerance = 1, 2, 3; for 13 TeV, 8 TeV and 5.02 TeV respectively for K* charged analysis
 P.S: these values have been calculated for K* charged at different energies. For other analysis, it's recommended to stick to the fix Mass Tolerance Cut (Date Modified: 4 May 2018 by Kunal Garg)

*/

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
   fpT_Tolerance(0),
   fMassTolSigma(0),
   fMass(0.0),
   fTolerance(0.01),
   fToleranceVeto(0.01),
   fSwitch(0),
   fLife(0),
   fLowRadius(0),
   fHighRadius(0),
   fMaxDCAVertex(0.3),
   fMinCosPointAngle(0.95),
   fMaxDaughtersDCA(0.5),
   fMinTPCcluster(70),
   fMaxRapidity(0.8),
   fMaxPseudorapidity(1e20),
   fPID(pid),
   fPID2(pid2),
   fPIDCutProton(0),
   fPIDCutPion(0),
   fESDtrackCuts(0x0),
   fCutQuality(Form("%sDaughtersQuality", name)),
   fAODTestFilterBit(5),
   fCustomTrackDCACuts(kFALSE),
   fMinDCAPositiveTrack(0.001),
   fMinDCANegativeTrack(0.001),
   fCheckOOBPileup(kFALSE)
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
   fpT_Tolerance(copy.fpT_Tolerance),
   fMassTolSigma(copy.fMassTolSigma),
   fMass(copy.fMass),
   fTolerance(copy.fTolerance),
   fToleranceVeto(copy.fToleranceVeto),
   fSwitch(copy.fSwitch),
   fLife(copy.fLife),
   fLowRadius(copy.fLowRadius),
   fHighRadius(copy.fHighRadius),
   fMaxDCAVertex(copy.fMaxDCAVertex),
   fMinCosPointAngle(copy.fMinCosPointAngle),
   fMaxDaughtersDCA(copy.fMaxDaughtersDCA),
   fMinTPCcluster(copy.fMinTPCcluster),
   fMaxRapidity(copy.fMaxRapidity),
   fMaxPseudorapidity(copy.fMaxPseudorapidity),
   fPID(copy.fPID),
   fPID2(copy.fPID2),
   fPIDCutProton(copy.fPIDCutProton),
   fPIDCutPion(copy.fPIDCutPion),
   fESDtrackCuts(copy.fESDtrackCuts),
   fCutQuality(copy.fCutQuality),
   fAODTestFilterBit(copy.fAODTestFilterBit),
   fCustomTrackDCACuts(copy.fCustomTrackDCACuts),
   fMinDCAPositiveTrack(copy.fMinDCAPositiveTrack),
   fMinDCANegativeTrack(copy.fMinDCANegativeTrack),
   fCheckOOBPileup(copy.fCheckOOBPileup)
{
//
// Copy constructor.
// Just copy all data member values.:IsSelected: Object is not a V0 (RESONANCES/AliRsnCutV0.cxx:149)

//
   fCutQuality.SetPtRange(0.15, 1E+20);
   fCutQuality.SetEtaRange(-0.8, 0.8);
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
   fpT_Tolerance = copy.fpT_Tolerance;
   fMassTolSigma = copy.fMassTolSigma;
   fMass = copy.fMass;
   fTolerance = copy.fTolerance;
   fToleranceVeto = copy.fToleranceVeto;
   fSwitch = copy.fSwitch;
   fLife  = copy.fLife;
   fLowRadius = copy.fLowRadius;
   fHighRadius  = copy.fHighRadius;
   fMaxDCAVertex = copy.fMaxDCAVertex;
   fMinCosPointAngle = copy.fMinCosPointAngle;
   fMaxDaughtersDCA = copy.fMaxDaughtersDCA;
   fMinTPCcluster = copy.fMinTPCcluster;
   fMaxRapidity = copy.fMaxRapidity;
   fMaxPseudorapidity = copy.fMaxPseudorapidity;
   fCutQuality = copy.fCutQuality;
   fPID = copy.fPID;
   fPID2 = copy.fPID2;
   fPIDCutProton = copy.fPIDCutProton;
   fPIDCutPion = copy.fPIDCutPion;
   fESDtrackCuts = copy.fESDtrackCuts;
   fCutQuality = copy.fCutQuality;
   fAODTestFilterBit = copy.fAODTestFilterBit;
   fCustomTrackDCACuts = copy.fCustomTrackDCACuts;
   fMinDCAPositiveTrack = copy.fMinDCAPositiveTrack;
   fMinDCANegativeTrack = copy.fMinDCANegativeTrack;
   fCheckOOBPileup = copy.fCheckOOBPileup;

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

   // The Cascade classes inherit from the V0 classes. Cascades will therefore be counted as V0s unless we remove them.
   AliESDcascade *Xesd = fDaughter->Ref2ESDcascade();
   AliAODcascade *Xaod = fDaughter->Ref2AODcascade();
   if (Xesd || Xaod) return kFALSE;

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

    // filter like-sign V0
   if ( TMath::Abs( ((pTrack->GetSign()) - (nTrack->GetSign())) ) < 0.1) {
      AliDebugClass(2, "Failed like-sign V0 check");
      return kFALSE;
   }

   if (fCheckOOBPileup) {
      Double_t bfield = lESDEvent->GetMagneticField();
      if(!TrackPassesOOBPileupCut(pTrack, bfield) &&
         !TrackPassesOOBPileupCut(nTrack, bfield)) return kFALSE;
   }

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

      // Apply different DCAxy cut for positive and negative V0 daughters

   if (fCustomTrackDCACuts) {
     Float_t impParPos[2], impParNeg[2];
     Float_t covMatPos[3], covMatNeg[3];
     pTrack->GetImpactParameters(impParPos,covMatPos);
     nTrack->GetImpactParameters(impParNeg,covMatNeg);

     if (covMatPos[0]<=0 || covMatPos[2]<=0) {
       Printf("Estimated b resolution lower or equal zero!");
       covMatPos[0]=0; covMatPos[2]=0;
     }

     if (covMatNeg[0]<=0 || covMatNeg[2]<=0) {
       Printf("Estimated b resolution lower or equal zero!");
       covMatNeg[0]=0; covMatNeg[2]=0;
     }
     
     if(TMath::Abs(impParPos[0])<fMinDCAPositiveTrack) return kFALSE;
     if(TMath::Abs(impParNeg[0])<fMinDCANegativeTrack) return kFALSE;
     
   }

   // topological checks
   if (TMath::Abs(v0->GetDcaV0Daughters()) > fMaxDaughtersDCA) {
     AliDebugClass(2, "Failed check on DCA between daughters");
     return kFALSE;
   }
   if (TMath::Abs(v0->GetD(xPrimaryVertex, yPrimaryVertex, zPrimaryVertex)) > fMaxDCAVertex) {
     AliDebugClass(2, "Failed check on DCA to primary vertes");
     return kFALSE;
   }
   //if (TMath::Abs(v0->GetV0CosineOfPointingAngle()) < fMinCosPointAngle) {
   if ( (TMath::Abs(v0->GetV0CosineOfPointingAngle()) < fMinCosPointAngle) || (TMath::Abs(v0->GetV0CosineOfPointingAngle()) >= 1 ) ) {
     AliDebugClass(2, "Failed check on cosine of pointing angle");
     return kFALSE;
   }
   if (TMath::Abs(v0->Y(fHypothesis)) > fMaxRapidity) {
      AliDebugClass(2, "Failed check on V0 rapidity");
      return kFALSE;
   }
   if (TMath::Abs(v0->Eta()) > fMaxPseudorapidity) {
      AliDebugClass(2, "Failed check on V0 pseudorapidity");
      return kFALSE;
   }
   //
   Double_t v0Position[3]; // from $ALICE_ROOT/ANALYSIS/AliESDV0Cuts.cxx
   v0->GetXYZ(v0Position[0],v0Position[1],v0Position[2]);
   Double_t radius = TMath::Sqrt(TMath::Power(v0Position[0],2) + TMath::Power(v0Position[1],2));
   if ( ( radius < fLowRadius ) || ( radius > fHighRadius ) ) {
     AliDebugClass(2, "Failed fiducial volume");
     return kFALSE;
   }

    // Lifetime cut for negative and positive track
    
    // Total Momentum
  Double_t tV0mom[3];
  v0->GetPxPyPz( tV0mom[0],tV0mom[1],tV0mom[2] );
  Double_t lV0TotalMomentum =
      TMath::Sqrt(tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );
    
    
   Double_t fLength = TMath::Sqrt(TMath::Power(v0Position[0]- xPrimaryVertex,2) + TMath::Power(v0Position[1] - yPrimaryVertex,2)+ TMath::Power(v0Position[2]- zPrimaryVertex,2));
   
   if( TMath::Abs(fMass*fLength/lV0TotalMomentum) > fLife)
   {
       AliDebugClass(2, "Failed Lifetime Cut on positive track V0");
       return kFALSE;
    }
    
   Double_t v0pT = TMath::Abs(TMath::Sqrt(tV0mom[0]*tV0mom[0] + tV0mom[1]*tV0mom[1]));
	
	if(fpT_Tolerance==0)
	{
		
		v0->ChangeMassHypothesis(fHypothesis);
		if ((TMath::Abs(v0->GetEffMass() - fMass)) > fTolerance) {
			AliDebugClass(2, "V0 is not in the expected inv mass range");
			return kFALSE;
		}
	}

	Double_t upper_limit =0, lower_limit =0;
	
    if(fpT_Tolerance==1 || fpT_Tolerance==2)		//Read the Note at the top of the code for more information .
	 {
		 v0->ChangeMassHypothesis(fHypothesis);

		 if(v0pT<=1.5)
		 {
			 if(v0pT<0.15){v0pT = 0.15;}
			 upper_limit = 0.49722 + 9.83285e-04* log(v0pT) + fMassTolSigma*(3.47153e-03 + 2.70453e-04* v0pT);
			 lower_limit = 0.49722 + 9.83285e-04* log(v0pT) - fMassTolSigma*(3.47153e-03 + 2.70453e-04* v0pT);
		}
		 
		 else if(v0pT > 1.5)
		 {
			 upper_limit = 0.49761 + fMassTolSigma*(3.47153e-03 + 2.70453e-04* v0pT);
			 lower_limit = 0.49761 - fMassTolSigma*(3.47153e-03 + 2.70453e-04* v0pT);
		 }
		 
		 if ((v0->GetEffMass() < lower_limit ) || (v0->GetEffMass() > upper_limit))
		 {
			 AliDebugClass(2, "V0 is not in the expected inv mass range");
			 return kFALSE;
		 }
	 }
	
	if(fpT_Tolerance==3)
	{
		v0->ChangeMassHypothesis(fHypothesis);
		Double_t upper_limit =0, lower_limit =0;
		
		if(v0pT<=1.5)
		{
			if(v0pT<0.15) {v0pT = 0.15;}
			upper_limit = 0.49797 + 9.99398e-04* log(v0pT) + fMassTolSigma*(3.47153e-03 + 2.70453e-04* v0pT);
			lower_limit = 0.49797 + 9.99398e-04* log(v0pT) - fMassTolSigma*(3.47153e-03 + 2.70453e-04* v0pT);
		}
		else if(v0pT > 1.5)
		{
			upper_limit = 0.498375 + fMassTolSigma*(3.47153e-03 + 2.70453e-04* v0pT);
			lower_limit = 0.498375 - fMassTolSigma*(3.47153e-03 + 2.70453e-04* v0pT);
		}
		
		if ((v0->GetEffMass() < lower_limit ) || (v0->GetEffMass() > upper_limit))
		{
			AliDebugClass(2, "V0 is not in the expected inv mass range");
			return kFALSE;
		}
	}

    //Set Switch to kTRUE to use Competing V0 Rejection
    if(fSwitch){
        
    if(fHypothesis == kK0Short)
    {v0->ChangeMassHypothesis(kLambda0);
     if ((TMath::Abs(v0->GetEffMass() - 1.115683)) < fToleranceVeto) {
         v0->ChangeMassHypothesis(kK0Short);
      return kFALSE;
   
     }
      v0->ChangeMassHypothesis(kK0Short);
    
    }

        if(fHypothesis == kLambda0)
    {v0->ChangeMassHypothesis(kK0Short);
     if ((TMath::Abs(v0->GetEffMass() - 0.497614)) < fToleranceVeto) {
         v0->ChangeMassHypothesis(kLambda0);
      return kFALSE;
   
     }
      v0->ChangeMassHypothesis(kLambda0);
    }
        
        if(fHypothesis == kLambda0Bar)
    {v0->ChangeMassHypothesis(kK0Short);
     if ((TMath::Abs(v0->GetEffMass() - 0.497614)) < fToleranceVeto) {
         v0->ChangeMassHypothesis(kLambda0Bar);
      return kFALSE;
   
     }
      v0->ChangeMassHypothesis(kLambda0Bar);
    
    }
        
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
   else if(fHypothesis==kLambda0Bar) {
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
   else if(fHypothesis==kK0Short) {
      //if (isTOFneg) {
      // TPC: 5sigma cut for all
      //if (negnsTPC > 5.0) return kFALSE;
      // TOF: 3sigma
      // maxTOF = 3.0;
      //return (negnsTOF <= maxTOF);
      //} else {
      // TPC:

         maxTPC = fPIDCutPion;
         maxTPC2 = fPIDCutPion;

      if(! ((negnsTPC <= maxTPC) && (posnsTPC2 <= maxTPC2)) ) {
         AliDebugClass(2, "Failed check on V0 PID");
         return kFALSE;
      }
   }

   // if we reach this point, all checks were successful
   AliDebugClass(2, "Good V0");
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

   //-----
   // next lines commented out by EF - 17/01/2014
   // NOTE that the filter bit test on V0 daughters removes a huge amount of V0 candidates, including good ones.
   //      Likely wrong -> requires a DCA max!
   //      Removing the test, there's a little gain in efficiency in the
   //      final search for Sigma* candidates
   // NOTE that further constrains (e.g. DCA of daughters greater than xxx),
   //      necessary to remove background, are already in V0s. (see also below)
   /*
   if ( !pTrack->TestFilterBit(fAODTestFilterBit)   ) {
      AliDebugClass(2, Form("Positive daughter failed quality cuts filtermapP=%d",filtermapP));
      return kFALSE;
   }
   if ( !nTrack->TestFilterBit(fAODTestFilterBit)   ) {
      AliDebugClass(2, Form("Negative daughter failed quality cuts filtermapN=%d",filtermapN));
      return kFALSE;
   }
   */

   //----
   // next lines are not necessary. Just left there (commented-out) to remind that the requirement on the DCA of V0 daughters
   //      is already in the V0, so requiring dca>0.050 (with 0.050 cm the default value from the Lambda analysis)
   //      does not remove V0s candidates
   /*
   Double_t dca = v0->DcaPosToPrimVertex() ;
     AliDebugClass(2, Form("DCA of Lambda positive daughter %f",dca));
   if(dca<0.050) {
     AliDebugClass(2, Form("DCA of Lambda positive daughter (%f) less than 0.05",dca));
      return kFALSE;
   }
   dca = v0->DcaNegToPrimVertex();
   if(dca<0.050) {
     AliDebugClass(2, Form("DCA of Lambda negative daughter (%f) less than 0.05",dca));
      return kFALSE;
   }
   */

   // EF - 17/01/2014 - next check apparently not effective!? Already in V0s?
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
   else if (fHypothesis==kK0Short) {
      mass = v0->MassK0Short();
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

   // next cut is effective (should it be in AODV0?)
   AliAODVertex *vertex = lAODEvent->GetPrimaryVertex();
   Double_t cospointangle = v0->CosPointingAngle(vertex);
   if (TMath::Abs( cospointangle )  < fMinCosPointAngle) {
     AliDebugClass(2, "Failed check on cosine of pointing angle");
     return kFALSE;
   }
 
  // next cut is effective (should it be in AODV0?)
   if (TMath::Abs(v0->DcaV0Daughters()) > fMaxDaughtersDCA) {
      AliDebugClass(2, "Failed check on DCA between daughters");
      return kFALSE;
   }

   if (TMath::Abs(v0->RapLambda()) > fMaxRapidity) {
      AliDebugClass(2, "Failed check on V0 rapidity");
      return kFALSE;
   }

   if (TMath::Abs(v0->Eta()) > fMaxPseudorapidity) {
      AliDebugClass(2, "Failed check on V0 pseusorapidity");
      return kFALSE;
   }

   Double_t radius = v0->RadiusV0();
   if ( ( radius < fLowRadius ) || ( radius > fHighRadius ) ){
     AliDebugClass(2, "Failed fiducial volume");
     return kFALSE;
   }

   //-----------------------------------------------------------
   // check initialization of PID object
   AliPIDResponse *pid = fEvent->GetPIDResponse();
   if (!pid) {
      AliFatal("NULL PID response");
      return kFALSE;
   }

   Double_t posnsTPC   = TMath::Abs(pid->NumberOfSigmasTPC(pTrack, fPID));
   Double_t posnsTPC2  = TMath::Abs(pid->NumberOfSigmasTPC(pTrack, fPID2));
   Double_t negnsTPC   = TMath::Abs(pid->NumberOfSigmasTPC(nTrack, fPID));
   Double_t negnsTPC2  = TMath::Abs(pid->NumberOfSigmasTPC(nTrack, fPID2));
   Double_t maxTPC = 1E20;
   Double_t maxTPC2 = 1E20;

   // applies the cut differently depending on the PID and the momentum
   if(fHypothesis==kLambda0) {
     maxTPC = fPIDCutProton;
     maxTPC2 = fPIDCutPion;
     if (! ((posnsTPC <= maxTPC) && (negnsTPC2 <= maxTPC2)) ) {
       AliDebugClass(2, "Failed check on V0 PID");
       return kFALSE;
     }
   }
   else if(fHypothesis==kLambda0Bar) {
     maxTPC = fPIDCutProton;
     maxTPC2 = fPIDCutPion;
     if(! ((negnsTPC <= maxTPC) && (posnsTPC2 <= maxTPC2)) ) {
       AliDebugClass(2, "Failed check on V0 PID");
       return kFALSE;
     }
   }
   else if(fHypothesis==kK0Short) {
     maxTPC = fPIDCutPion;
     maxTPC2 = fPIDCutPion;
     if(! ((negnsTPC <= maxTPC) && (posnsTPC2 <= maxTPC2)) ) {
       AliDebugClass(2, "Failed check on V0 PID");
       return kFALSE;
     }
   }
   
   //---------------------------------------------------------------
   // if we reach this point, all checks were successful
   AliDebugClass(1, "Good AOD V0");
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

//_________________________________________________________________________________________________
Bool_t AliRsnCutV0::TrackPassesOOBPileupCut(AliESDtrack* t, Double_t b){
   if (!t) return true;
   if ((t->GetStatus() & AliESDtrack::kITSrefit) == AliESDtrack::kITSrefit) return true;
   if (t->GetTOFExpTDiff(b, true) + 2500 > 1e-6) return true;
   return false;
}
