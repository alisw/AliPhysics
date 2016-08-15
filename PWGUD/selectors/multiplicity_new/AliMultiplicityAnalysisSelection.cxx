#include "AliMultiplicityAnalysisSelection.h"
#include "AliMultiplicityHelper.h"

ClassImp ( AliMultiplicityAnalysisSelection )

AliMultiplicityAnalysisSelection::AliMultiplicityAnalysisSelection ( const char* name, const char* title ) : AliMultiplicityEventSelection ( name, title ),
fEtaCut(1.0),
fZvCut(5.5),
useZcutTS(kFALSE),
fCorrelateTracks(kFALSE),
checkIfPhysSelected(kFALSE),
fPhiCut(0),
fSPDInCuts(kFALSE),
fGlobalInCuts(kFALSE),
fMCVInCuts(kFALSE),
fPreSelected(kFALSE) {
	fZvCutTS[0] = - fZvCut;
	fZvCutTS[1] = fZvCut;
}

Bool_t AliMultiplicityAnalysisSelection::IsVertexInCuts ( const AliESDEvent* esd, UInt_t mask ) {
	
	fZcache[0] = esd->GetPrimaryVertexSPD()->GetZ();
	fZcache[1] = esd->GetPrimaryVertexTracks()->GetZ();
	fZcache[2] = esd->GetPrimaryVertex()->GetZ();
	
	if ( useZcutTS ) {
		fGlobalInCuts = ( (Bool_t) ( ( mask & AliMultiplicityHelper::kVGlobal ) != 0 ) ) && ( (AliMultiplicityHelper::GetVertex ( AliMultiplicityHelper::kVertexTracks, esd )->GetZ() > GetZvCutTS(0)) && ( AliMultiplicityHelper::GetVertex ( AliMultiplicityHelper::kVertexTracks, esd )->GetZ() < GetZvCutTS(1)) );
		fGlobalInCuts = ( (Bool_t) ( ( mask & AliMultiplicityHelper::kVSPD ) != 0 ) ) && ( (AliMultiplicityHelper::GetVertex ( AliMultiplicityHelper::kVertexSPD, esd )->GetZ() > GetZvCutTS(0)) && ( AliMultiplicityHelper::GetVertex ( AliMultiplicityHelper::kVertexSPD, esd )->GetZ() < GetZvCutTS(1)) );
	} else {
		fGlobalInCuts = ( (Bool_t) ( ( mask & AliMultiplicityHelper::kVGlobal ) != 0 ) ) && ( TMath::Abs ( AliMultiplicityHelper::GetVertex ( AliMultiplicityHelper::kVertexTracks, esd )->GetZ() ) < GetZvCut() );
		fSPDInCuts = ( (Bool_t) ( ( mask & AliMultiplicityHelper::kVSPD ) != 0) ) && ( TMath::Abs ( AliMultiplicityHelper::GetVertex ( AliMultiplicityHelper::kVertexSPD, esd )->GetZ() ) < GetZvCut() );
	}
	
	return fGlobalInCuts || fSPDInCuts;
}

Bool_t AliMultiplicityAnalysisSelection::IsVertexInCuts ( const AliMCEvent* mc ) {
	fZcacheGen = mc->GetPrimaryVertex()->GetZ();
	if ( useZcutTS ) {
		fMCVInCuts = ( ( ( mc->GetPrimaryVertex()->GetZ() > GetZvCutTS(0) ) ) && ( mc->GetPrimaryVertex()->GetZ() < GetZvCutTS(1) ) );
	} else {
		fMCVInCuts = ( TMath::Abs ( mc->GetPrimaryVertex()->GetZ() ) < GetZvCut() ) ;
	}
	return fMCVInCuts;
}

Double_t AliMultiplicityAnalysisSelection::GetEtaCut() {
	return fEtaCut;
}

Double_t AliMultiplicityAnalysisSelection::GetZvCut() {
	return fZvCut;
}

Double_t AliMultiplicityAnalysisSelection::GetZvCutTS ( Int_t side ) {
	if ( ( side > -1 ) && ( side < 2 ) ) return fZvCutTS[side];
	AliWarning ( "Wrong side" );
	return -100.0;
}


void AliMultiplicityAnalysisSelection::SetEtaCut ( Double_t etaCut ) {
	if ( etaCut > 0 ) fEtaCut = etaCut;
}

void AliMultiplicityAnalysisSelection::SetZvCut ( Double_t ZvCut ) {
	if ( ZvCut > 0 ) fZvCut = ZvCut;
	useZcutTS = kFALSE;
}

void AliMultiplicityAnalysisSelection::SetZvCutTS ( Double_t zlow, Double_t zhi ) {
	useZcutTS = kTRUE;
	fZvCutTS[0] = zlow;
	fZvCutTS[1] = zhi;
}


Bool_t AliMultiplicityAnalysisSelection::GetQueryPhysSelection() {
	return checkIfPhysSelected;
}

void AliMultiplicityAnalysisSelection::SetQueryPhysSelection ( Bool_t query ) {
	checkIfPhysSelected = query;
}

Bool_t AliMultiplicityAnalysisSelection::GetCorrelateTracks() {
	return fCorrelateTracks;
}

void AliMultiplicityAnalysisSelection::SetCorrelateTracks ( Bool_t corr ) {
	fCorrelateTracks = corr;
}

void AliMultiplicityAnalysisSelection::ResetAccepted() {
	fSPDInCuts = kFALSE;
	fGlobalInCuts = kFALSE;
	fMCVInCuts = kFALSE;
	AliMultiplicityEventSelection::ResetAccepted();
	fPreSelected = kFALSE;
	ResetCounters();
}
