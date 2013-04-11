// -------------------------------------------------------------------------
// Copies all info that is needed for the DiHadronPID analysis.
// Possible Extension: At this moment the object is protected for returning
// pointers to the original tracks. It could at some point be beneficial to
// be able to access this information.
// -------------------------------------------------------------------------
// Author: Misha Veldhoen (misha.veldhoen@cern.ch)

#include <iostream>
using namespace std;

// AOD includes.
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODMCParticle.h"

// PID includes.
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"

// Objects own include.
#include "AliTrackDiHadronPID.h"

ClassImp(AliTrackDiHadronPID);

// -------------------------------------------------------------------------
AliTrackDiHadronPID::AliTrackDiHadronPID():
	TObject(),
	fAODTrack(0x0),
	fAODGlobalTrack(0x0),
	fAODEvent(0x0),
	fAODMCParticle(0x0),
	fPIDResponse(0x0),
	fBasicInfoAvailable(kFALSE),
	fFlagsAvailable(kFALSE),
	fDCAInfoAvailable(kFALSE),
	fITSInfoAvailable(kFALSE),
	fTPCInfoAvailable(kFALSE),
	fTOFInfoAvailable(kFALSE),
	fMCInfoAvailable(kFALSE),
	fPt(-999.),
	fEta(-999.),
	fPhi(-999.),
	fFlags(0),
	fFilterMap(0),
	fID(0),
	fLabel(0),
	fCharge(0),
	fDCAz(-999.),
	fDCAxy(-999.),
	fTOFsignal(-999.),
	fIsTOFmismatch(kFALSE),
	fTPCsignal(-999.),
	fTPCmomentum(-999.),
	fMCPt(-999.),
	fMCEta(-999.),
	fMCPhi(-999.),
	fMCY(-999.),
	fPdgCode(0),
	fIsPhysicalPrimary(kFALSE),
	fIsSecondaryFromWeakDecay(kFALSE),
	fIsSecondaryFromMaterial(kFALSE),
	fDebug(0)

{

	//
	// Default Constructor.
	//

	if (fDebug > 2) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {
		fTOFsignalMinusExpected[iSpecies] = -999.;
		fTOFNsigma[iSpecies] = -999.;
		fTPCsignalMinusExpected[iSpecies] = -999.;
		fTPCNsigma[iSpecies] = -999.;
		fY[iSpecies] = -999.;
	}

	for (Int_t iITSlayer = 0; iITSlayer < 6; iITSlayer++) {
		fITSHits[iITSlayer] = kFALSE;
	}

}

// -------------------------------------------------------------------------
AliTrackDiHadronPID::AliTrackDiHadronPID(AliAODTrack* track, AliAODTrack* globaltrack, AliAODMCParticle* mcparticle, AliPIDResponse* pidresponse):
	TObject(),
	fAODTrack(0x0),
	fAODGlobalTrack(0x0),
	fAODEvent(0x0),
	fAODMCParticle(0x0),
	fPIDResponse(0x0),
	fBasicInfoAvailable(kFALSE),
	fFlagsAvailable(kFALSE),
	fDCAInfoAvailable(kFALSE),
	fITSInfoAvailable(kFALSE),
	fTPCInfoAvailable(kFALSE),
	fTOFInfoAvailable(kFALSE),
	fMCInfoAvailable(kFALSE),
	fPt(-999.),
	fEta(-999.),
	fPhi(-999.),
	fFlags(0),
	fFilterMap(0),
	fID(0),
	fLabel(0),
	fCharge(0),
	fDCAz(-999.),
	fDCAxy(-999.),
	fTOFsignal(-999.),
	fIsTOFmismatch(kFALSE),
	fTPCsignal(-999.),
	fTPCmomentum(-999.),
	fMCPt(-999.),
	fMCEta(-999.),
	fMCPhi(-999.),
	fMCY(-999.),	
	fPdgCode(0),
	fIsPhysicalPrimary(kFALSE),
	fIsSecondaryFromWeakDecay(kFALSE),
	fIsSecondaryFromMaterial(kFALSE),
	fDebug(0)
{

	//
	// Constructor.
	//

	if (fDebug > 2) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {
		fTOFsignalMinusExpected[iSpecies] = -999.;
		fTOFNsigma[iSpecies] = -999.;
		fTPCsignalMinusExpected[iSpecies] = -999.;
		fTPCNsigma[iSpecies] = -999.;
		fY[iSpecies] = -999.;	
	}

	for (Int_t iITSlayer = 0; iITSlayer < 6; iITSlayer++) {
		fITSHits[iITSlayer] = kFALSE;
	}

	if (track) {
		fAODTrack = track;
		fAODEvent = track->GetAODEvent();
	}
	if (globaltrack) fAODGlobalTrack = globaltrack;
	if (mcparticle) fAODMCParticle = mcparticle;
	if (pidresponse) fPIDResponse = pidresponse;

	// Copy AOD Track info.
	if (fAODTrack) {
		CopyBasicTrackInfo();
	} else {
		AliError("No Track Supplied.");
	}

	// Find the Global Track.
	if (fID >= 0) fAODGlobalTrack = fAODTrack;

	// Copy DCA and PID info.
	if (fAODGlobalTrack) {
		CopyFlags();
		if (fAODEvent) CopyDCAInfo();
		else AliError("Couln't find AOD Event.");
		CopyITSInfo();
		if (fPIDResponse) CopyTPCInfo();
		CopyTOFInfo();
	} else {
		AliError("Couldn't find Global Track.");
	} 

	// Copy MC info.
	if (fAODMCParticle) {
		CopyMCInfo();
	} 

	// Test 
	/*	Double_t sigmaTOFProton = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(fAODTrack, AliPID::kProton));
		if ( sigmaTOFProton < 1.0) {cout<<"tofsigmabelowone: "<<sigmaTOFProton<<endl;}
	
		Double_t sigmaTPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fAODTrack, AliPID::kProton));
		if ( sigmaTPCProton < 1.0) {cout<<"tpcsigmabelowone: "<<sigmaTPCProton<<endl;}*/
}

// -------------------------------------------------------------------------
Bool_t AliTrackDiHadronPID::CopyBasicTrackInfo() {

	//
	// Copies everything available in every AOD track.
	//

	if (fDebug > 2) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

 	fPt = fAODTrack->Pt();
	fEta = fAODTrack->Eta();
	fPhi = fAODTrack->Phi();

	fY[0] = fAODTrack->Y(AliAODTrack::kPion);
	fY[1] = fAODTrack->Y(AliAODTrack::kKaon);
	fY[2] = fAODTrack->Y(AliAODTrack::kProton);

	//fFlags = fAODTrack->GetFlags(); // FLAGS MUST BE TAKEN FROM GLOBAL TRACKS.
	fFilterMap = fAODTrack->GetFilterMap();

	fID = fAODTrack->GetID();
	fLabel = fAODTrack->GetLabel();

	fCharge = fAODTrack->Charge();

	fBasicInfoAvailable = kTRUE;
	return fBasicInfoAvailable;

}

// -------------------------------------------------------------------------
Bool_t AliTrackDiHadronPID::CopyFlags() {

	//
	// Copies Flags (properly stored in global track)
	//

	if (fDebug > 2) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Copy Flags
	fFlags = fAODGlobalTrack->GetFlags();

	// Is TOF mismatch?
	if (AliAODTrack::kTOFmismatch&fFlags) {
		fIsTOFmismatch = kTRUE;
		//cout<<"Found TOF mismatch!"<<endl;
	}
	else fIsTOFmismatch = kFALSE; 

	fFlagsAvailable = kTRUE;
	return fFlagsAvailable;

}

// -------------------------------------------------------------------------
Bool_t AliTrackDiHadronPID::CopyDCAInfo() {

	//
	// Copies DCA info. (only stored in a global track)
	//

	if (fDebug > 2) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Propagate track to DCA.
	Double_t PosAtDCA[2] = {-999,-999};
    Double_t covar[3] = {-999,-999,-999};
    //cout<<fAODTrack<<" "<<fAODGlobalTrack<<endl;
    Bool_t propagate = fAODGlobalTrack->PropagateToDCA(fAODEvent->GetPrimaryVertex(),fAODEvent->GetMagneticField(),100.,PosAtDCA,covar);
    	
    if (propagate) {
    	fDCAxy = PosAtDCA[0];
    	fDCAz = PosAtDCA[1];
    } else {
    	//AliError("Could not propagate track to DCA.");
    }

    if (propagate) fDCAInfoAvailable = kTRUE;
    return fDCAInfoAvailable;

}

// -------------------------------------------------------------------------
Bool_t AliTrackDiHadronPID::CopyITSInfo() {

	//
	// Copies ITS info.
	//

	if (fDebug > 2) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

    // Get the ITS clustermap
    fITSClusterMap = fAODGlobalTrack->GetITSClusterMap();

    // Copy the ITS hits.
    for (Int_t iITSlayer = 0; iITSlayer < 6; iITSlayer++) {
		fITSHits[iITSlayer] = fAODGlobalTrack->HasPointOnITSLayer(iITSlayer);
	}

    fITSInfoAvailable = kTRUE;
    return fITSInfoAvailable;

}

// -------------------------------------------------------------------------
Bool_t AliTrackDiHadronPID::CopyTPCInfo() {

	//
	// Copies TPC info. (needs global track and pid response)
	//

	if (fDebug > 2) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}	

    // Get TPC signal.
    fTPCsignal = fAODGlobalTrack->GetTPCsignal();

    // Compute expected TPC signal under pi/K/p mass assumption.
    AliTPCPIDResponse& TPCPIDResponse = fPIDResponse->GetTPCResponse();
    fTPCmomentum = fAODGlobalTrack->GetTPCmomentum();

	fTPCsignalMinusExpected[0] = fTPCsignal - TPCPIDResponse.GetExpectedSignal(fTPCmomentum,AliPID::kPion);
	fTPCsignalMinusExpected[1] = fTPCsignal - TPCPIDResponse.GetExpectedSignal(fTPCmomentum,AliPID::kKaon);
	fTPCsignalMinusExpected[2] = fTPCsignal - TPCPIDResponse.GetExpectedSignal(fTPCmomentum,AliPID::kProton);

	fTPCNsigma[0] = fPIDResponse->NumberOfSigmasTPC(fAODGlobalTrack, AliPID::kPion);
	fTPCNsigma[1] = fPIDResponse->NumberOfSigmasTPC(fAODGlobalTrack, AliPID::kKaon);
	fTPCNsigma[2] = fPIDResponse->NumberOfSigmasTPC(fAODGlobalTrack, AliPID::kProton);

    fTPCInfoAvailable = kTRUE;
    return fTPCInfoAvailable;

}

// -------------------------------------------------------------------------
Bool_t AliTrackDiHadronPID::CopyTOFInfo() {

	//
	// Copies TOF info. (needs global track)
	//

	if (fDebug > 2) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

    // Get TOF signal.
    fTOFsignal = fAODGlobalTrack->GetTOFsignal();

    // Compute expected TOF signal under pi/K/p mass assumption.
    Double_t times[AliPID::kSPECIES];
    fAODGlobalTrack->GetIntegratedTimes(times);
    fTOFsignalMinusExpected[0] = fTOFsignal - times[AliPID::kPion];
	fTOFsignalMinusExpected[1] = fTOFsignal - times[AliPID::kKaon];
	fTOFsignalMinusExpected[2] = fTOFsignal - times[AliPID::kProton];

	fTOFNsigma[0] = fPIDResponse->NumberOfSigmasTOF(fAODGlobalTrack, AliPID::kPion);
	fTOFNsigma[1] = fPIDResponse->NumberOfSigmasTOF(fAODGlobalTrack, AliPID::kKaon);
	fTOFNsigma[2] = fPIDResponse->NumberOfSigmasTOF(fAODGlobalTrack, AliPID::kProton);	

	fTOFInfoAvailable = kTRUE;
	return fTOFInfoAvailable;

}

// -------------------------------------------------------------------------
Bool_t AliTrackDiHadronPID::CopyMCInfo() {

	// Copies MC info (needs an MC track with the same label)

	// Check if the label of the current track matches the label of the
	// generated particle. Note that the label of the AOD track can be
	// negative. This means that the quality of this track is not awesome,
	// but that it does correspond to the MC particle.

	if (fDebug > 2) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}
	
	/*
	if (fAODMCParticle->Label() != TMath::Abs(fAODTrack->GetLabel())) {
		cout<<"Label of supplied MC particle and reconstructed track do not match."<<endl;	
		return kFALSE;
	}
	*/
	// Note: It seems like the Label of the AOD track points to the INDEX of the
	// MCPArticle, not to the label (See AliAnalysisTaskCompareAODTrackCuts.cxx)

	fMCPt = fAODMCParticle->Pt();
	fMCEta = fAODMCParticle->Eta();
	fMCPhi = fAODMCParticle->Phi();
	fMCY = fAODMCParticle->Y();
	fPdgCode = fAODMCParticle->PdgCode();

	TClonesArray* mcArray = 0x0;
	mcArray = dynamic_cast<TClonesArray*>(fAODEvent->FindListObject(AliAODMCParticle::StdBranchName()));
   	if (!mcArray) {
        AliFatal("No MC array found in the AOD.");
    }

	// Primary particle
	if ( fAODMCParticle->IsPhysicalPrimary() ){
		fIsPhysicalPrimary = kTRUE;
	} else {
		// Safety check for mother existence.
		if (fAODMCParticle->GetMother() >= 0){

			Int_t mcMotherPDG = -999;
			Int_t firstInt = -999;

			AliAODMCParticle* mcMother = (AliAODMCParticle*) mcArray->At(TMath::Abs(fAODMCParticle->GetMother()));
			mcMotherPDG = TMath::Abs(mcMother->GetPdgCode());

			// Need a way to get the first intiger, for now Marek's method:
			firstInt = Int_t (mcMotherPDG/ TMath::Power(10, Int_t(TMath::Log10(mcMotherPDG))));
			// cout<<"Mother PDG: "<<mcMotherPDG<<"; Firt integer: "<<firstInt<<endl;

			// Weak decay
			if( firstInt == 3){
				fIsSecondaryFromWeakDecay = kTRUE;
			// Material decay
			} else {
				fIsSecondaryFromMaterial = kTRUE;
			}
		}
	}

	fMCInfoAvailable = kTRUE;
	return fMCInfoAvailable;

}
