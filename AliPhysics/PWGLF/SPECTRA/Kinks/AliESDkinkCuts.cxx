#include "AliESDkinkCuts.h"
#include <limits>
#include <TCanvas.h>
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TH1.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TMCProcess.h>
#include <TPDGCode.h>
#include <TTree.h>
#include <TDirectory.h>
#include <AliAnalysisManager.h>
#include <AliESDEvent.h>
#include <AliESDInputHandler.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliESDtrackCuts.h>
#include <AliCFParticleGenCuts.h>
#include <AliPIDResponse.h>
#include <AliCFParticleGenCuts.h>
#include <AliESDpidCuts.h>
#include <AliESDkink.h>
#include <AliStack.h>
#include <AliMCParticle.h>
#include <AliMCVertex.h>
#include <AliCFTrackKineCuts.h>
#include <AliKineTrackCuts.h>
#include <AliAODEvent.h>
#include <AliAODInputHandler.h>
#include <AliAODTrack.h>
#include <AliAODMCParticle.h>
#include <AliMultiplicity.h>

// Kink and resonance analysis task
// Author: Filimon Roukoutakis, University of Athens

ClassImp(AliESDkinkCuts)

//________________________________________________________________________
Bool_t AliESDkinkCuts::IsGoodKaonKink(const AliESDkink* kink, const AliESDEvent* esdEvent, /* Future general function AliESDkinkCuts* esdKinkCuts,*/ AliVTrack*& kinkMother, AliVTrack*& kinkDaughter) const
{
  
		if (!kink) return kFALSE;
		//if ( (kink->GetQt() < fMinQt) || (kink->GetQt() > fMaxQt) ) return kFALSE;
		if ( !IsInQtLimits(kink->GetQt()) ) return kFALSE;
		if ( !GetKinkTracks(kink, esdEvent, kinkMother, kinkDaughter) ) return kFALSE;
		if ( fEsdTrackCutsKinkMother && !fEsdTrackCutsKinkMother->IsSelected(kinkMother) ) return kFALSE;
		if ( fEsdTrackCutsKinkDaughter && !fEsdTrackCutsKinkDaughter->IsSelected(kinkDaughter) ) return kFALSE;
    Float_t kinkAngle = TMath::RadToDeg() * kink->GetAngle(2);
		if (kinkAngle < fMinAngleKaon) return kFALSE;
		/*AliExternalTrackParam* kinkMotherTPCParam = (AliExternalTrackParam *)kinkMother->GetTPCInnerParam();
    if (!kinkMotherTPCParam) return kFALSE;*/
	  TVector3 mother3MomentumAtKinkVtx(kink->GetMotherP());
	  TVector3 daughter3MomentumAtKinkVtx(kink->GetDaughterP());
		TVector3 transferedMom = mother3MomentumAtKinkVtx-daughter3MomentumAtKinkVtx;
		Float_t energyDaughterMu = TMath::Sqrt(daughter3MomentumAtKinkVtx.Mag()*daughter3MomentumAtKinkVtx.Mag()+fMuonMass*fMuonMass);
    Float_t invariantMassKmu = (energyDaughterMu+transferedMom.Mag())*(energyDaughterMu+transferedMom.Mag())-mother3MomentumAtKinkVtx.Mag()*mother3MomentumAtKinkVtx.Mag();
    if ( /*(mother3MomentumAtKinkVtx.Mag() < fMinKaonMomentum) &&*/ ( (invariantMassKmu < 0) || (TMath::Sqrt(invariantMassKmu) > fInvMassKmuMaxAngle) ) ) return kFALSE; // Kmu cut	
		//if(!IsKink(esd, kinkMother->GetKinkIndex(0), kinkMother3Momentum)) return kFALSE;
		Double_t maxDecAngKmu = fMaxDecayAngleCurveKmu->Eval(mother3MomentumAtKinkVtx.Mag(), 0., 0., 0.);
		Double_t maxDecAngpimu = fMaxDecayAngleCurvePimu->Eval(mother3MomentumAtKinkVtx.Mag(), 0., 0., 0.);
		if ( kinkAngle < maxDecAngpimu * 1.2 ) return kFALSE; // Possibly pion, reject above theoretical pion kink curve
		if( (mother3MomentumAtKinkVtx.Mag() > fMinKaonMomentum) && ( kinkAngle > maxDecAngKmu * 0.98) ) return kFALSE; // Reject above theoretical kaon kink curve
		// Here try various cuts for mother3MomentumAtKinkVtx.Mag() < fMinKaonMomentum (mother momentum at kink vertex below good kink kinematical region)
		//if( (mother3MomentumAtKinkVtx.Mag() < fMinKaonMomentum) && ( kink->GetQt() < 0.03) ) return kFALSE; // Hard qt cut for p<1.1 only, bad
		if ( /*(mother3MomentumAtKinkVtx.Mag() < fMinKaonMomentum) &&*/ !IsSelected(kinkMother, AliPID::kKaon) ) return kFALSE; // Optional dEdx cut, MSS compatible yields at 2.76
		//if ( /*(mother3MomentumAtKinkVtx.Mag() < fMinKaonMomentum) &&*/ (TMath::Abs(pidResponse->NumberOfSigmasTPC(kinkDaughter, AliPID::kElectron)) < 3.5) ) return kTRUE; // NEW! Daughter electron cut
		//if ( (mother3MomentumAtKinkVtx.Mag() > 0.5) && (mother3MomentumAtKinkVtx.Mag() < fMinKaonMomentum) && (TMath::RadToDeg()*kink->GetAngle(2) > 3) ) return kFALSE; // Bad
		//if ( (mother3MomentumAtKinkVtx.Mag() < fMinKaonMomentum) && ((kinkMother->GetTPCNcls() + kinkDaughter->GetTPCNcls()) < 120) ) return kFALSE; // To test
		//if (TMath::Abs(pidResponse->NumberOfSigmasTPC(kinkDaughter, AliPID::kMuon)) > 3.0) return kFALSE;

		return kTRUE;
}

//________________________________________________________________________
Bool_t AliESDkinkCuts::IsGoodPionKink(const AliESDkink* kink, const AliESDEvent* esdEvent, /* Future general function AliESDkinkCuts* esdKinkCuts,*/ AliVTrack*& kinkMother, AliVTrack*& kinkDaughter) const
{
  
		if (!kink) return kFALSE;
		//if ( (kink->GetQt() < fMinQt) || (kink->GetQt() > fMaxQt) ) return kFALSE;
		if ( !IsInQtLimitsPion(kink->GetQt()) ) return kFALSE;
		if ( !GetKinkTracks(kink, esdEvent, kinkMother, kinkDaughter) ) return kFALSE;
		if ( fEsdTrackCutsKinkMother && !fEsdTrackCutsKinkMother->IsSelected(kinkMother) ) return kFALSE;
		if ( fEsdTrackCutsKinkDaughter && !fEsdTrackCutsKinkDaughter->IsSelected(kinkDaughter) ) return kFALSE;
    Float_t kinkAngle = TMath::RadToDeg() * kink->GetAngle(2);
		if (kinkAngle < fMinAnglePion) return kFALSE;
		/*AliExternalTrackParam* kinkMotherTPCParam = (AliExternalTrackParam *)kinkMother->GetTPCInnerParam();
    if (!kinkMotherTPCParam) return kFALSE;*/
	  TVector3 mother3MomentumAtKinkVtx(kink->GetMotherP());
	  TVector3 daughter3MomentumAtKinkVtx(kink->GetDaughterP());
		TVector3 transferedMom = mother3MomentumAtKinkVtx-daughter3MomentumAtKinkVtx;
		Float_t energyDaughterMu = TMath::Sqrt(daughter3MomentumAtKinkVtx.Mag()*daughter3MomentumAtKinkVtx.Mag()+fMuonMass*fMuonMass);
    //Float_t invariantMassKmu = (energyDaughterMu+transferedMom.Mag())*(energyDaughterMu+transferedMom.Mag())-mother3MomentumAtKinkVtx.Mag()*mother3MomentumAtKinkVtx.Mag();
    //if ( /*(mother3MomentumAtKinkVtx.Mag() < fMinKaonMomentum) &&*/ ( (invariantMassKmu < 0) || (TMath::Sqrt(invariantMassKmu) > fInvMassKmuMaxAngle) ) ) return kFALSE; // Kmu cut	
		//if(!IsKink(esd, kinkMother->GetKinkIndex(0), kinkMother3Momentum)) return kFALSE;
		//Double_t maxDecAngKmu = fMaxDecayAngleCurveKmu->Eval(mother3MomentumAtKinkVtx.Mag(), 0., 0., 0.);
		Double_t maxDecAngpimu = fMaxDecayAngleCurvePimu->Eval(mother3MomentumAtKinkVtx.Mag(), 0., 0., 0.);
		if ( kinkAngle > maxDecAngpimu ) return kFALSE; // Possibly kaon, reject above theoretical pion kink curve
		//if( (mother3MomentumAtKinkVtx.Mag() > fMinKaonMomentum) && ( kinkAngle > maxDecAngKmu) ) return kFALSE; // Reject above theoretical kaon kink curve
		// Here try various cuts for mother3MomentumAtKinkVtx.Mag() < fMinKaonMomentum (mother momentum at kink vertex below good kink kinematical region)
		//if( (mother3MomentumAtKinkVtx.Mag() < fMinKaonMomentum) && ( kink->GetQt() < 0.03) ) return kFALSE; // Hard qt cut for p<1.1 only, bad
		if ( /*(mother3MomentumAtKinkVtx.Mag() < fMinKaonMomentum) &&*/ !IsSelected(kinkMother, AliPID::kPion) ) return kFALSE; // Optional dEdx cut, MSS compatible yields at 2.76
		//if ( /*(mother3MomentumAtKinkVtx.Mag() < fMinKaonMomentum) &&*/ (TMath::Abs(pidResponse->NumberOfSigmasTPC(kinkDaughter, AliPID::kElectron)) < 3.5) ) return kTRUE; // NEW! Daughter electron cut
		//if ( (mother3MomentumAtKinkVtx.Mag() > 0.5) && (mother3MomentumAtKinkVtx.Mag() < fMinKaonMomentum) && (TMath::RadToDeg()*kink->GetAngle(2) > 3) ) return kFALSE; // Bad
		//if ( (mother3MomentumAtKinkVtx.Mag() < fMinKaonMomentum) && ((kinkMother->GetTPCNcls() + kinkDaughter->GetTPCNcls()) < 120) ) return kFALSE; // To test
		//if (TMath::Abs(pidResponse->NumberOfSigmasTPC(kinkDaughter, AliPID::kMuon)) > 3.0) return kFALSE;

		return kTRUE;
}
