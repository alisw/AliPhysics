#ifndef AliESDkinkCuts_cxx
#define AliESDkinkCuts_cxx

// Generic kink kinematics cut class
// Author: Filimon Roukoutakis, University of Athens
// Helper class to identify K and pi from their kink decays to mu. 

#include <AliAnalysisTaskSE.h>
#include <AliVEvent.h>
#include <AliInputEventHandler.h>
#include <AliPID.h>
#include <THnSparse.h>
#include <TF1.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TVector3.h>
#include <AliAnalysisCuts.h>
#include <AliTrackReference.h>
#include <AliMCParticle.h>
#include <AliESDVertex.h>
#include <AliESDEvent.h>
#include <AliESDpidCuts.h>
#include <AliESDkink.h>
#include <AliPIDResponse.h>
#include "AliESDkinkCuts.h"

//class TH1F;
class TH3;
class AliESDtrackCuts;
class AliCFParticleGenCuts;
class AliKineTrackCuts;
class AliVVertex;
class AliESDtrack;
//class TF1;
//class AliESDkink;
//class THnSparseF;
//class AliESDEvent;
//class AliPIDResponse;
//class AliESDpidCuts;

	class AliESDkinkCuts : public AliAnalysisCuts {
	public:
    enum EFiducialCuts { kProcessedKinksFiducial = 1, kMinAbsZ, kMaxAbsZ, kMinR, kMaxR, kFiducialCutsMax };
    enum EKineCuts { kProcessedKinksKine = 1, kMinPt, kAbsY, kAbsEta, kKineCutsMax };
    enum EKinkCuts { kProcessedKinksKink = 1, kMinQt, kMaxQt, kMinAngle, kMaxAngle, kKinkCutsMax };
		AliESDkinkCuts() : AliAnalysisCuts(), fEsdTrackCutsKinkMother(0), fEsdTrackCutsKinkDaughter(0), fPIDResponse(0), fEsdPIDCutsArray(), fMuonMass(0), fMinKaonMomentum(0), fMinPionMomentum(0), fMaxKaonMomentum(0), fMaxPionMomentum(0), fMinPt(0), fMaxAbsY(0), fMaxAbsEta(0), fMinQt(0), fMaxQt(0), fMinQtPion(0), fMaxQtPion(0), fMinKinkR(0), fMaxKinkR(0), fMinAbsKinkZ(0), fMaxAbsKinkZ(0), fMinAngleKaon(0), fMinAnglePion(0), fInvMassKmuMaxAngle(0), fMaxDecayAngleCurveKmu(0), fMaxDecayAngleCurvePimu(0), fFillCutHist(kFALSE), fhCutFiducial(0), fhCutFiducialKine(0), fhCutKinkKaon(0), fhCutKinkPion(0) {};
		AliESDkinkCuts(const char* name, const char* title) : AliAnalysisCuts(name, title), fEsdTrackCutsKinkMother(0), fEsdTrackCutsKinkDaughter(0), fPIDResponse(0), fEsdPIDCutsArray(), fMuonMass(0), fMinKaonMomentum(0), fMinPionMomentum(0), fMaxKaonMomentum(0), fMaxPionMomentum(0), fMinPt(0), fMaxAbsY(0), fMaxAbsEta(0), fMinQt(0), fMaxQt(0), fMinQtPion(0), fMaxQtPion(0), fMinKinkR(0), fMaxKinkR(0), fMinAbsKinkZ(0), fMaxAbsKinkZ(0), fMinAngleKaon(0), fMinAnglePion(0), fInvMassKmuMaxAngle(0), fMaxDecayAngleCurveKmu(0), fMaxDecayAngleCurvePimu(0), fFillCutHist(kFALSE), fhCutFiducial(0), fhCutFiducialKine(0), fhCutKinkKaon(0), fhCutKinkPion(0) {};
		AliESDkinkCuts(const char* name, const char* title, AliESDtrackCuts* esdTrackCutsKinkMother, AliESDtrackCuts* esdTrackCutsKinkDaughter/*, AliKineCuts* fiducialDecayCut*/, const AliPIDResponse* pidResponse, /*AliESDpidCuts**/ const TArrayF esdPIDcutsKinks, const Bool_t fillCutHist, const Double_t minPt, const Double_t maxAbsY, const Double_t minQt, const Double_t maxQt, const Double_t minQtPion, const Double_t maxQtPion, const Double_t minKinkR=120.0, const Double_t maxKinkR=210.0, const Double_t minAbsKinkZ=0.5, const Double_t maxAbsKinkZ=225.0, const Double_t minAngleKaon=2.0, const Double_t minAnglePion=1.0) : AliAnalysisCuts(name, title), fEsdTrackCutsKinkMother(esdTrackCutsKinkMother), fEsdTrackCutsKinkDaughter(esdTrackCutsKinkDaughter), fPIDResponse(pidResponse), fEsdPIDCutsArray(esdPIDcutsKinks), fMuonMass(TDatabasePDG::Instance()->GetParticle("mu+")->Mass()), fMinKaonMomentum(1.20), fMinPionMomentum(0.01), fMaxKaonMomentum(100.0), fMaxPionMomentum(100.0), fMinPt(minPt), fMaxAbsY(maxAbsY), fMaxAbsEta(1.5), fMinQt(minQt), fMaxQt(maxQt), fMinQtPion(minQtPion), fMaxQtPion(maxQtPion), fMinKinkR(minKinkR), fMaxKinkR(maxKinkR), fMinAbsKinkZ(minAbsKinkZ), fMaxAbsKinkZ(maxAbsKinkZ), fMinAngleKaon(minAngleKaon), fMinAnglePion(minAnglePion), fInvMassKmuMaxAngle(0.6), fMaxDecayAngleCurveKmu(0), fMaxDecayAngleCurvePimu(0), fFillCutHist(fillCutHist), fhCutFiducial(0), fhCutFiducialKine(0), fhCutKinkKaon(0), fhCutKinkPion(0) {
			   fMaxDecayAngleCurveKmu=new TF1("fMaxDecayAngleCurveKmu","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]", fMinKaonMomentum, 1.1*fMaxKaonMomentum);
         fMaxDecayAngleCurveKmu->SetParameter(0, TDatabasePDG::Instance()->GetParticle("K+")->Mass());
         fMaxDecayAngleCurveKmu->SetParameter(1, 0.9127037);
         fMaxDecayAngleCurveKmu->SetParameter(2, TMath::Pi());
         fMaxDecayAngleCurvePimu=new TF1("fMaxDecayAngleCurvePimu","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]", fMinPionMomentum, 1.1*fMaxPionMomentum);
         fMaxDecayAngleCurvePimu->SetParameter(0, TDatabasePDG::Instance()->GetParticle("pi+")->Mass());
         fMaxDecayAngleCurvePimu->SetParameter(1, 0.2731374);
         fMaxDecayAngleCurvePimu->SetParameter(2, TMath::Pi());
/*
//kaon to muon maximum kink angle
TF1 *f1=new TF1("f1","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",1.1,15.0);
f1->SetParameter(0,0.493677);
f1->SetParameter(1,0.9127037);
f1->SetParameter(2,TMath::Pi());

//kaon to pion maximum kink angle
TF1 *f2=new TF1("f2","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",0.8,15.0);
f2->SetParameter(0,0.493677);
f2->SetParameter(1,0.826607291);
f2->SetParameter(2,TMath::Pi());

//kaon to electron maximum kink angle
TF1 *f3=new TF1("f3","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",0.3,15.0);
f3->SetParameter(0,0.493677);
f3->SetParameter(1,0.435193217);
f3->SetParameter(2,TMath::Pi());

//pion to muon maximum kink angle
TF1 *f4=new TF1("f4","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",0.1,15.0);
f4->SetParameter(0,0.13957018);
f4->SetParameter(1,0.2731374);
f4->SetParameter(2,TMath::Pi());
*/
	// Cuts
	fhCutFiducial = new TH1F(Form("%s::fhCutFiducial", GetName()), Form("%s; Cut type; Number of events", GetTitle()), kFiducialCutsMax, kProcessedKinksFiducial, kFiducialCutsMax);
	fhCutFiducial->SetOption("b text0");
	TAxis* axis = fhCutFiducial->GetXaxis();
	axis->SetBinLabel(kProcessedKinksFiducial, "all");
	axis->SetBinLabel(kMinAbsZ, "minabs(z)");
	axis->SetBinLabel(kMaxAbsZ, "maxabs(z)");
	axis->SetBinLabel(kMinR, "minR");
	axis->SetBinLabel(kMaxR, "maxR");
	fhCutFiducialKine = new TH1F(Form("%s::fhCutFiducialKine", GetName()), Form("%s; Cut type; Number of events", GetTitle()), kKineCutsMax, kProcessedKinksKine, kKineCutsMax);
	fhCutFiducialKine->SetOption("b text0");
	axis = fhCutFiducialKine->GetXaxis();
	axis->SetBinLabel(kProcessedKinksKine, "all");
	axis->SetBinLabel(kMinPt, "minPt");
	axis->SetBinLabel(kAbsY, "abs(Y)");
	axis->SetBinLabel(kAbsEta, "abs(Eta)");
	fhCutKinkKaon = new TH1F(Form("%s::fhCutKinkKaon", GetName()), Form("%s; Cut type; Number of events", GetTitle()), kKinkCutsMax, kProcessedKinksKink, kKinkCutsMax);
	fhCutKinkKaon->SetOption("b text0");
	axis = fhCutKinkKaon->GetXaxis();
	axis->SetBinLabel(kProcessedKinksKink, "all");
	axis->SetBinLabel(kMinQt, "minQt");
	axis->SetBinLabel(kMaxQt, "maxQt");
	axis->SetBinLabel(kMinAngle, "minAngle");
	axis->SetBinLabel(kMaxAngle, "maxAngle");
	fhCutKinkPion = new TH1F(Form("%s::fhCutKinkPion", GetName()), Form("%s; Cut type; Number of events", GetTitle()), kKinkCutsMax, kProcessedKinksKink, kKinkCutsMax);
	fhCutKinkPion->SetOption("b text0");
	axis = fhCutKinkPion->GetXaxis();
	axis->SetBinLabel(kProcessedKinksKink, "all");
	axis->SetBinLabel(kMinQt, "minQt");
	axis->SetBinLabel(kMaxQt, "maxQt");
	axis->SetBinLabel(kMinAngle, "minAngle");
	axis->SetBinLabel(kMaxAngle, "maxAngle");
		};
    virtual Bool_t IsSelected(const AliVParticle *track, const AliPID::EParticleType type) const { return( (fPIDResponse && (fEsdPIDCutsArray.GetSize() == AliPID::kSPECIES) ) ? (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, type)) < fEsdPIDCutsArray[type]) : kTRUE); };
		virtual Bool_t IsSelected(TObject*) { return kFALSE; };
		virtual Bool_t IsSelected(TList*) { return kFALSE; };
		virtual Bool_t IsSelected(const AliESDtrack* /*esdTrack*/) const { return kFALSE; };
		virtual Bool_t IsSelected(const AliESDkink* /*esdKink*/) const { return kFALSE; };		
	inline Bool_t IsInQtLimits(const Double_t qt) const {
		if (fFillCutHist) fhCutKinkKaon->Fill(kMinQt);
		if ( qt < fMinQt ) return kFALSE;
		if (fFillCutHist) fhCutKinkKaon->Fill(kMaxQt);
		if ( qt > fMaxQt ) return kFALSE;
		return kTRUE;
	};
	inline Bool_t IsInQtLimitsPion(const Double_t qt) const {
		if (fFillCutHist) fhCutKinkPion->Fill(kMinQt);
		if ( qt < fMinQtPion ) return kFALSE;
		if (fFillCutHist) fhCutKinkPion->Fill(kMaxQt);
		if ( qt > fMaxQtPion ) return kFALSE;
		return kTRUE;
	};
	inline Bool_t IsInFiducial(const TVector3& position) const {
		if (fFillCutHist) fhCutFiducial->Fill(kProcessedKinksFiducial);
		if (TMath::Abs(position.Z()) < fMinAbsKinkZ) return kFALSE;
		if (fFillCutHist) fhCutFiducial->Fill(kMinAbsZ);
		if (TMath::Abs(position.Z()) > fMaxAbsKinkZ) return kFALSE;
		if (fFillCutHist) fhCutFiducial->Fill(kMaxAbsZ);
		if (position.Perp() < fMinKinkR) return kFALSE;
		if (fFillCutHist) fhCutFiducial->Fill(kMinR);
		if (position.Perp() > fMaxKinkR) return kFALSE; 
		if (fFillCutHist) fhCutFiducial->Fill(kMaxR);
		return kTRUE;
	};
	inline Bool_t IsInFiducialKine(const TLorentzVector& momentum) const { 
		if (fFillCutHist) fhCutFiducialKine->Fill(kProcessedKinksKine);
		if (momentum.Pt() < fMinPt) return kFALSE;
		if (fFillCutHist) fhCutFiducialKine->Fill(kMinPt);
		if (TMath::Abs(momentum.Rapidity()) > fMaxAbsY) return kFALSE;
		if (fFillCutHist) fhCutFiducialKine->Fill(kAbsY);
		if (TMath::Abs(momentum.Eta()) > fMaxAbsEta) return kFALSE; 
		if (fFillCutHist) fhCutFiducialKine->Fill(kAbsEta);
		return kTRUE;
	};
	inline Bool_t GetKinkTracks(const AliESDkink* kink, const AliESDEvent* esdEvent, AliVTrack*& kinkMother, AliVTrack*& kinkDaughter) const {
		if(!kink || !esdEvent) return kFALSE;
		Int_t motherIndex = kink->GetIndex(0);
		Int_t daughterIndex = kink->GetIndex(1);
		kinkMother = esdEvent->GetTrack(TMath::Abs(motherIndex));
		kinkDaughter = esdEvent->GetTrack(TMath::Abs(daughterIndex));
		if ( !kinkMother || !kinkDaughter || (kinkMother == kinkDaughter) ) return kFALSE;
		return kTRUE;
	}
	Bool_t IsGoodKaonKink(const AliESDkink* kink, const AliESDEvent* esdEvent, /* Future general function AliESDkinkCuts* esdKinkCuts,*/ AliVTrack*& kinkMother, AliVTrack*& kinkDaughter) const;
	Bool_t IsGoodPionKink(const AliESDkink* kink, const AliESDEvent* esdEvent, /* Future general function AliESDkinkCuts* esdKinkCuts,*/ AliVTrack*& kinkMother, AliVTrack*& kinkDaughter) const;
	inline void SetPIDResponse(AliPIDResponse* pidResponse) { fPIDResponse=pidResponse; };
	private:
		Bool_t IsKinematicallyKaonKink(const AliESDkink* /*esdKink*/) const { return kFALSE; };
		Bool_t IsKinkMotherPrimary(const AliESDtrack* /*esdTrack*/) const { return kFALSE; };
		AliESDtrackCuts* fEsdTrackCutsKinkMother;
		AliESDtrackCuts* fEsdTrackCutsKinkDaughter;
	  const AliPIDResponse* fPIDResponse;
		const TArrayF fEsdPIDCutsArray;
    AliESDkinkCuts(const AliESDkinkCuts& rhs);
    AliESDkinkCuts& operator=(const AliESDkinkCuts& rhs);
	public:
    // Move to private
	  const Float_t fMuonMass;
		const Double_t fMinKaonMomentum;
		const Double_t fMinPionMomentum;
		const Double_t fMaxKaonMomentum;
		const Double_t fMaxPionMomentum;
    Double_t fMinPt;
		Double_t fMaxAbsY;
		Double_t fMaxAbsEta;
		Double_t fMinQt;
		Double_t fMaxQt;
		Double_t fMinQtPion;
		Double_t fMaxQtPion;
		Double_t fMinKinkR;
		Double_t fMaxKinkR;
		Double_t fMinAbsKinkZ;
		Double_t fMaxAbsKinkZ;
		Double_t fMinAngleKaon;
		Double_t fMinAnglePion;
		Double_t fInvMassKmuMaxAngle;
		TF1* fMaxDecayAngleCurveKmu;
		TF1* fMaxDecayAngleCurvePimu;
	// Cut hist
	Bool_t fFillCutHist;
	TH1* fhCutFiducial;
	TH1* fhCutFiducialKine;
	TH1* fhCutKinkKaon;
	TH1* fhCutKinkPion; 
  ClassDef(AliESDkinkCuts, 1); // example of analysis
	};

#endif
