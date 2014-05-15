/************************************************************************* 
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. * 
*                                                                        * 
* Author: The ALICE Off-line Project.                                    * 
* Contributors are mentioned in the code where appropriate.              * 
*                                                                        * 
* Permission to use, copy, modify and distribute this software and its   * 
* documentation strictly for non-commercial purposes is hereby granted   * 
* without fee, provided that the above copyright notice appears in all   * 
* copies and that both the copyright notice and this permission notice   * 
* appear in the supporting documentation. The authors make no claims     * 
* about the suitability of this software for any purpose. It is          * 
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/

// -----------------------------------------------------------------------
//  Track class for the DiHadronPID analysis.
// -----------------------------------------------------------------------
//  Author: Misha Veldhoen (misha.veldhoen@cern.ch)

#include "AliTrackDiHadronPID.h"

#include "AliAODVertex.h"
#include "AliPID.h"
#include "AliTPCPIDResponse.h"

ClassImp(AliTrackDiHadronPID);

Double_t AliTrackDiHadronPID::fSigmaTOFStd = 80.;	// Should perhaps be replaced with a 
Double_t AliTrackDiHadronPID::fSigmaTPCStd = 3.5;   // function later.

// -----------------------------------------------------------------------
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
	fNclsTPC(-999),
	fDCAz(-999.),
	fDCAxy(-999.),
	fTOFsignal(-999.),
	fTOFMatchingStatus(-1),
	fTPCsignal(-999.),
	fTPCmomentum(-999.),
	fITSClusterMap(0),
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

	for (Int_t iN = 0; iN < 3; ++iN) {
		fTOFLabel[iN] = -1;	// Same convention as in ESDs
	}

}

// -----------------------------------------------------------------------
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
	fNclsTPC(-999),
	fDCAz(-999.),
	fDCAxy(-999.),
	fTOFsignal(-999.),
	fTOFMatchingStatus(-1),
	fTPCsignal(-999.),
	fTPCmomentum(-999.),
	fITSClusterMap(0),	
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

	for (Int_t iN = 0; iN < 3; ++iN) {
		fTOFLabel[iN] = -1;	// Same convention as in ESDs
	}

	if (track) {
		fAODTrack = track;
		fAODEvent = const_cast<AliAODEvent*>(track->GetAODEvent());
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

	// Copy the rest of the track parameters if the filtermap is nonzero.
	// If fFiltermap == 0, then propagation to the DCA will result in a floating point error.
	if (fFilterMap) {

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

	}	

	// Test 
	/*	Double_t sigmaTOFProton = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(fAODTrack, AliPID::kProton));
		if ( sigmaTOFProton < 1.0) {cout<<"tofsigmabelowone: "<<sigmaTOFProton<<endl;}
	
		Double_t sigmaTPCProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fAODTrack, AliPID::kProton));
		if ( sigmaTPCProton < 1.0) {cout<<"tpcsigmabelowone: "<<sigmaTPCProton<<endl;}*/
}

// -----------------------------------------------------------------------
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
	fNclsTPC = fAODTrack->GetTPCNcls();

	fBasicInfoAvailable = kTRUE;
	return fBasicInfoAvailable;

}

// -----------------------------------------------------------------------
Bool_t AliTrackDiHadronPID::CopyFlags() {

	//
	// Copies Flags (properly stored in global track)
	//

	if (fDebug > 2) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Copy Flags
	fFlags = fAODGlobalTrack->GetFlags();
/*
	// Is TOF mismatch?
	if (AliAODTrack::kTOFmismatch&fFlags) {
		fTOFMatchingStatus = kTRUE;
		//cout<<"Found TOF mismatch!"<<endl;
	}
	else fTOFMatchingStatus = kFALSE; 
*/
	fFlagsAvailable = kTRUE;
	return fFlagsAvailable;

}

// -----------------------------------------------------------------------
Bool_t AliTrackDiHadronPID::CopyDCAInfo() {

	//
	// Copies DCA info. (only stored in a global track)
	//

	if (fDebug > 2) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	if (fAODGlobalTrack->IsMuonTrack()) return kFALSE;

	// Propagate track to DCA.
	Double_t PosAtDCA[2] = {-999,-999};
    Double_t covar[3] = {-999,-999,-999};
    //cout<<fAODTrack<<" "<<fAODGlobalTrack<<endl;
    AliAODTrack* clone = (AliAODTrack*) fAODGlobalTrack->Clone("trk_clone"); //need clone, in order not to change track parameters
    Bool_t propagate = clone->PropagateToDCA(fAODEvent->GetPrimaryVertex(),fAODEvent->GetMagneticField(),100.,PosAtDCA,covar);
    delete clone;    	

    if (propagate) {
    	fDCAxy = PosAtDCA[0];
    	fDCAz = PosAtDCA[1];
    } else {
    	//AliError("Could not propagate track to DCA.");
    }

    if (propagate) fDCAInfoAvailable = kTRUE;
    return fDCAInfoAvailable;

}

// -----------------------------------------------------------------------
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

// -----------------------------------------------------------------------
Bool_t AliTrackDiHadronPID::CopyTPCInfo() {

	//
	// Copies TPC info. (needs global track and pid response).
	// See https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PIDInAnalysis#Signal_Deltas
	// for more info!
	//

	if (fDebug > 2) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}	

	// Get TPC signal and momentum.
	fTPCsignal = fAODGlobalTrack->GetTPCsignal();
	fTPCmomentum = fAODGlobalTrack->GetTPCmomentum();

	// Obtaining (signal - expected).
	fPIDResponse->GetSignalDelta(AliPIDResponse::kTPC, fAODGlobalTrack, AliPID::kPion, fTPCsignalMinusExpected[0], kFALSE);
	fPIDResponse->GetSignalDelta(AliPIDResponse::kTPC, fAODGlobalTrack, AliPID::kKaon, fTPCsignalMinusExpected[1], kFALSE);
	fPIDResponse->GetSignalDelta(AliPIDResponse::kTPC, fAODGlobalTrack, AliPID::kProton, fTPCsignalMinusExpected[2], kFALSE);

	// Obtaining nSigma.
	fTPCNsigma[0] = fPIDResponse->NumberOfSigmasTPC(fAODGlobalTrack, AliPID::kPion);
	fTPCNsigma[1] = fPIDResponse->NumberOfSigmasTPC(fAODGlobalTrack, AliPID::kKaon);
	fTPCNsigma[2] = fPIDResponse->NumberOfSigmasTPC(fAODGlobalTrack, AliPID::kProton);

    fTPCInfoAvailable = kTRUE;
    return fTPCInfoAvailable;

}

// -----------------------------------------------------------------------
Bool_t AliTrackDiHadronPID::CopyTOFInfo() {

	//
	// Copies TOF info. (needs global track)
	// See https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PIDInAnalysis#Signal_Deltas
	// for more info!	
	//

	if (fDebug > 2) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

    // Get TOF signal minus the start time.
	fTOFsignal = fAODGlobalTrack->GetTOFsignal();

	// Get the expected times.
	Double_t expectedTimes[AliPID::kSPECIES];
	fAODGlobalTrack->GetIntegratedTimes(expectedTimes);
/*
	// Get the exptected TOF resolution.
	AliTOFHeader* tofH = (AliTOFHeader*)ev->GetTOFHeader();
	Double_t TOFpidRes[AliPID::kSPECIES];
	tr->GetDetPid()->GetTOFpidResolution(TOFpidRes);
*/
	// Obtaining (signal - expected).
	fPIDResponse->GetSignalDelta(AliPIDResponse::kTOF, fAODGlobalTrack, AliPID::kPion, fTOFsignalMinusExpected[0], kFALSE);
	fPIDResponse->GetSignalDelta(AliPIDResponse::kTOF, fAODGlobalTrack, AliPID::kKaon, fTOFsignalMinusExpected[1], kFALSE);
	fPIDResponse->GetSignalDelta(AliPIDResponse::kTOF, fAODGlobalTrack, AliPID::kProton, fTOFsignalMinusExpected[2], kFALSE);

	// Obtaining nSigma.
	fTOFNsigma[0] = fPIDResponse->NumberOfSigmasTOF(fAODGlobalTrack, AliPID::kPion);
	fTOFNsigma[1] = fPIDResponse->NumberOfSigmasTOF(fAODGlobalTrack, AliPID::kKaon);
	fTOFNsigma[2] = fPIDResponse->NumberOfSigmasTOF(fAODGlobalTrack, AliPID::kProton);	

	// Q: what do the different TOF labels mean?
	// It seems that in AOD090 the TOF labels aren't copied properly.
	//Int_t TOFlabeltmp[3] = {0};
	fAODGlobalTrack->GetTOFLabel(fTOFLabel);
	//for (Int_t iN = 0; iN < 3; ++iN) {fTOFLabel[iN] = TOFlabeltmp[iN];}
	/*
	if (fTOFLabel[1] == fLabel || fTOFLabel[2] == fLabel) {
		cout<<"fLabel = " << fLabel << " fTOFLabel =  {" << fTOFLabel[0] << "," << fTOFLabel[1] << "," << fTOFLabel[2] <<"}"<<endl; 
	}
	*/
	// The following will only work in an AOD production with the fTOFlabels set.
	// If it wasn't set, then every track will be labeled as no match.
	if (fTOFLabel[0] == -1) {fTOFMatchingStatus = 2;} 			// TPC Track was not matched to any TOF hit.
	else if (fLabel == fTOFLabel[0]) {fTOFMatchingStatus = 0;}	// TPC Track was correctly matched to a TOF hit.
	else {fTOFMatchingStatus = 1;}								// TPC Track was mismatched.

	fTOFInfoAvailable = kTRUE;
	return fTOFInfoAvailable;

}

// -----------------------------------------------------------------------
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

// -----------------------------------------------------------------------
Bool_t AliTrackDiHadronPID::UnknownSpecies(Int_t species) const {

	if (fDebug > 2) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}
	if (species < 0 || species > 2) {
		cout<<"ERROR: Unknown species"<<endl;
		return kTRUE;
	} else {
		return kFALSE;
	}

}
