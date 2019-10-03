/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *																								*
 * Authors: Friederike Bock															*
 * Version 1.0																				*
 *																								*
 * Permission to use, copy, modify and distribute this software and its	 *
 * documentation strictly for non-commercial purposes is hereby granted	 *
 * without fee, provided that the above copyright notice appears in all	 *
 * copies and that both the copyright notice and this permission notice	 *
 * appear in the supporting documentation. The authors make no claims	 *
 * about the suitability of this software for any purpose. It is			*
 * provided "as is" without express or implied warranty.						*
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// QA Task for V0 Reader V1
//---------------------------------------------
////////////////////////////////////////////////

#include "AliAnalysisTaskResolution.h"
#include "TChain.h"
#include "AliAnalysisManager.h"
#include "TParticle.h"
#include "TVectorF.h"
#include "AliPIDResponse.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "TFile.h"

class iostream;

using namespace std;

ClassImp(AliAnalysisTaskResolution)

//________________________________________________________________________
AliAnalysisTaskResolution::AliAnalysisTaskResolution() : AliAnalysisTaskSE(),
	fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
	fConversionGammas(NULL),
	fEventCuts(NULL),
	fConversionCuts(NULL),
	fTreeEvent(NULL),
	fTreeResolution(NULL),
	fPrimVtxZ(0.),
	fNContrVtx(0),
	fNESDtracksEta09(0),
	fNESDtracksEta0914(0),
	fNESDtracksEta14(0),
	fGammaRecCoords(5),
	fGammaMCCoords(5),
	fChi2ndf(0),
	fIsHeavyIon(0),
	fIsMC(kFALSE),
	fOutputList(NULL),
	fEventList(NULL),
	fResolutionList(NULL),
	fESDEvent(NULL),
	fMCEvent(NULL)
{

   
}

//________________________________________________________________________
AliAnalysisTaskResolution::AliAnalysisTaskResolution(const char *name) : AliAnalysisTaskSE(name),
	fV0Reader(NULL),
    fV0ReaderName("V0ReaderV1"),
	fConversionGammas(NULL),
	fEventCuts(NULL),
	fConversionCuts(NULL),
	fTreeEvent(NULL),
	fTreeResolution(NULL),
	fPrimVtxZ(0.),
	fNContrVtx(0),
	fNESDtracksEta09(0),
	fNESDtracksEta0914(0),
	fNESDtracksEta14(0),
	fGammaRecCoords(5),
	fGammaMCCoords(5),
	fChi2ndf(0),
	fIsHeavyIon(0),
	fIsMC(kFALSE),
	fOutputList(NULL),
	fEventList(NULL),
	fResolutionList(NULL),
	fESDEvent(NULL),
	fMCEvent(NULL)
{
	// Default constructor

	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskResolution::~AliAnalysisTaskResolution()
{
	// default deconstructor
   
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskResolution::Notify()
{
    if (fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){        
        fEventCuts->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      fEventCuts->SetPeriodEnum(fV0Reader->GetPeriodName());
    }  
    return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskResolution::UserCreateOutputObjects()
{
	// Create User Output Objects

	if(fOutputList != NULL){
		delete fOutputList;
		fOutputList = NULL;
	}
	if(fOutputList == NULL){
		fOutputList = new TList();
		fOutputList->SetOwner(kTRUE);
	}
	
	fEventList = new TList();
	fEventList->SetName("EventList");
	fEventList->SetOwner(kTRUE);
	fOutputList->Add(fEventList);
	
	fTreeEvent = new TTree("Event","Event");   
	fTreeEvent->Branch("primVtxZ",&fPrimVtxZ,"fPrimVtxZ/F");
	fTreeEvent->Branch("nContrVtx",&fNContrVtx,"fNContrVtx/I");
	fTreeEvent->Branch("nGoodTracksEta09",&fNESDtracksEta09,"fNESDtracksEta09/I");
	fTreeEvent->Branch("nGoodTracksEta14",&fNESDtracksEta14,"fNESDtracksEta14/I");
	fEventList->Add(fTreeEvent);
	
	if (fIsMC){		
		fResolutionList = new TList();
		fResolutionList->SetName("ResolutionList");
		fResolutionList->SetOwner(kTRUE);
		fOutputList->Add(fResolutionList);
								
		fTreeResolution = new TTree("Resolution","Resolution");
		fTreeResolution->Branch("nGoodTracksEta09",&fNESDtracksEta09,"fNESDtracksEta09/I");
		fTreeResolution->Branch("RecCoords",&fGammaRecCoords);
		fTreeResolution->Branch("MCCoords",&fGammaMCCoords);
		fTreeResolution->Branch("chi2ndf",&fChi2ndf,"fChi2ndf/F");
		fResolutionList->Add(fTreeResolution);
	}
	PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskResolution::UserExec(Option_t *){

    fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());

	Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
	if(eventQuality != 0){// Event Not Accepted
		return;
	}
	fESDEvent = (AliESDEvent*) InputEvent();
	if (fESDEvent==NULL) return;
	if(MCEvent()){
		fMCEvent = MCEvent();
	}
 
	if(MCEvent()){
		// Process MC Particle
		if(fEventCuts->GetSignalRejection() != 0){
// 		if(fESDEvent->IsA()==AliESDEvent::Class()){
			fEventCuts->GetNotRejectedParticles(	fEventCuts->GetSignalRejection(),
													fEventCuts->GetAcceptedHeader(),
													fMCEvent);
// 		}
// 		else if(fInputEvent->IsA()==AliAODEvent::Class()){
// 			fEventCuts->GetNotRejectedParticles(	fEventCuts->GetSignalRejection(),
// 													fEventCuts->GetAcceptedHeader(),
// 													fInputEvent);
// 		}
		}
	}

   
    if(fIsHeavyIon > 0 && !fEventCuts->IsCentralitySelected(fESDEvent,fMCEvent)) return;
	fNESDtracksEta09 = CountTracks09(); // Estimate Event Multiplicity
	fNESDtracksEta0914 = CountTracks0914(); // Estimate Event Multiplicity
	fNESDtracksEta14 = fNESDtracksEta09 + fNESDtracksEta0914;
 	if(fESDEvent){
		if(fESDEvent->GetPrimaryVertexTracks()->GetNContributors()>0) {
			fNContrVtx = fESDEvent->GetPrimaryVertexTracks()->GetNContributors();
		} else {
			fNContrVtx = 0;
		}	
// 		else if(fESDEvent->GetPrimaryVertexTracks()->GetNContributors()<1) {
// 			if(fESDEvent->GetPrimaryVertexSPD()->GetNContributors()>0) {
// 				fNContrVtx = fESDEvent->GetPrimaryVertexSPD()->GetNContributors();
// 			} else if(fESDEvent->GetPrimaryVertexSPD()->GetNContributors()<1) {
// 				fNContrVtx = 0;
// 			}
// 		}
	}
	fPrimVtxZ = fESDEvent->GetPrimaryVertex()->GetZ();
	
// 	if (fIsHeavyIon == 2){
// 		if (!(fNESDtracksEta09 > 20 && fNESDtracksEta09 < 80)) return;
// 	}	

	
	if (fTreeEvent){
		fTreeEvent->Fill();
	}

	fConversionGammas=fV0Reader->GetReconstructedGammas();
	ProcessPhotons();
	PostData(1, fOutputList);
}


///________________________________________________________________________
void AliAnalysisTaskResolution::ProcessPhotons(){

	// Fill Histograms for QA and MC
	for(Int_t firstGammaIndex=0;firstGammaIndex<fConversionGammas->GetEntriesFast();firstGammaIndex++){
		AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(firstGammaIndex));
		if (gamma ==NULL) continue;
		if(!fConversionCuts->PhotonIsSelected(gamma,fESDEvent)) continue;
		fNESDtracksEta09 = CountTracks09(); // Estimate Event Multiplicity
		fGammaRecCoords(0) = gamma->GetPhotonPt();
		fGammaRecCoords(1) = gamma->GetPhotonPhi();
		fGammaRecCoords(2) = gamma->GetPhotonEta();
		fGammaRecCoords(3) = gamma->GetConversionRadius();
		fGammaRecCoords(4) = gamma->GetConversionZ();
		fChi2ndf = gamma->GetChi2perNDF();
		if(MCEvent()){
            TParticle *posDaughter = gamma->GetPositiveMCDaughter(fMCEvent);
            TParticle *negDaughter = gamma->GetNegativeMCDaughter(fMCEvent);
// 			cout << "generate Daughters: "<<posDaughter << "\t" << negDaughter << endl;
            if(fMCEvent && fEventCuts->GetSignalRejection() != 0){
                Int_t isPosFromMBHeader = fEventCuts->IsParticleFromBGEvent(gamma->GetMCLabelPositive(), fMCEvent, fESDEvent);
                Int_t isNegFromMBHeader = fEventCuts->IsParticleFromBGEvent(gamma->GetMCLabelNegative(), fMCEvent, fESDEvent);
				if( (isNegFromMBHeader < 1) || (isPosFromMBHeader < 1)) continue;
			}
			
			if(posDaughter == NULL || negDaughter == NULL){ 
				continue;
			} else if(posDaughter->GetMother(0) != negDaughter->GetMother(0) || (posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1)){ 
				continue;
			} else {
// 				cout << "same mother" << endl;
				Int_t pdgCodePos; 
				if (posDaughter->GetPdgCode()) pdgCodePos = posDaughter->GetPdgCode(); else continue;
				Int_t pdgCodeNeg; 
				if (negDaughter->GetPdgCode()) pdgCodeNeg = negDaughter->GetPdgCode(); else continue;
// 				cout << "PDG codes daughters: " << pdgCodePos << "\t" << pdgCodeNeg << endl;
				Int_t pdgCode; 
                if (gamma->GetMCParticle(fMCEvent)->GetPdgCode()) pdgCode = gamma->GetMCParticle(fMCEvent)->GetPdgCode(); else continue;
// 				cout << "PDG code: " << pdgCode << endl;
				if(TMath::Abs(pdgCodePos)!=11 || TMath::Abs(pdgCodeNeg)!=11)
					continue;
				else if ( !(pdgCodeNeg==pdgCodePos)){
                    TParticle *truePhotonCanditate = gamma->GetMCParticle(fMCEvent);
					if(pdgCode == 111) 
						continue;
					else if (pdgCode == 221) 
						continue;
					else if (!(negDaughter->GetUniqueID() != 5 || posDaughter->GetUniqueID() !=5)){
						if(pdgCode == 22){
							fGammaMCCoords(0) = truePhotonCanditate->Pt();
                            fGammaMCCoords(1) = gamma->GetNegativeMCDaughter(fMCEvent)->Phi();
                            fGammaMCCoords(2) = gamma->GetNegativeMCDaughter(fMCEvent)->Eta();
                            fGammaMCCoords(3) = gamma->GetNegativeMCDaughter(fMCEvent)->R();
                            fGammaMCCoords(4) = gamma->GetNegativeMCDaughter(fMCEvent)->Vz();
							
							if (fTreeResolution){
								fTreeResolution->Fill();
							}
						}		
					} else continue; //garbage
				} else continue; //garbage
			}
		}
	}
}

//________________________________________________________________________
Int_t AliAnalysisTaskResolution::CountTracks09(){
	Int_t fNumberOfESDTracks = 0;
	if(fInputEvent->IsA()==AliESDEvent::Class()){
	// Using standard function for setting Cuts
			
		Bool_t selectPrimaries=kTRUE;
		AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
		EsdTrackCuts->SetMaxDCAToVertexZ(2);
		EsdTrackCuts->SetEtaRange(-0.9, 0.9);
		EsdTrackCuts->SetPtRange(0.15);
		
		for(Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++){
			AliESDtrack* curTrack = (AliESDtrack*) fInputEvent->GetTrack(iTracks);
			if(!curTrack) continue;
			if(EsdTrackCuts->AcceptTrack(curTrack) ){
				if (fMCEvent){
                    if(fEventCuts->GetSignalRejection() != 0){
                        Int_t isFromMBHeader = fEventCuts->IsParticleFromBGEvent(TMath::Abs(curTrack->GetLabel()), fMCEvent, fESDEvent);
						if( (isFromMBHeader < 1) ) continue;
					}					
				}	
				fNumberOfESDTracks++;
			}	
		}
		delete EsdTrackCuts;
		EsdTrackCuts=0x0;
	} else if(fInputEvent->IsA()==AliAODEvent::Class()){
		for(Int_t iTracks = 0; iTracks<fInputEvent->GetNumberOfTracks(); iTracks++){
			AliAODTrack* curTrack = (AliAODTrack*) fInputEvent->GetTrack(iTracks);
			if(!curTrack->IsPrimaryCandidate()) continue;
            if(TMath::Abs(curTrack->Eta())>0.9) continue;
			if(curTrack->Pt()<0.15) continue;
            if(TMath::Abs(curTrack->ZAtDCA())>2) continue;
			fNumberOfESDTracks++;
		}
	}

	return fNumberOfESDTracks;
}

//________________________________________________________________________
Int_t AliAnalysisTaskResolution::CountTracks0914(){
	Int_t fNumberOfESDTracks = 0;
	if(fInputEvent->IsA()==AliESDEvent::Class()){
		// Using standard function for setting Cuts
		
		AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
		EsdTrackCuts->SetMaxDCAToVertexZ(5);
		EsdTrackCuts->SetEtaRange(0.9, 1.4);
		EsdTrackCuts->SetPtRange(0.15);
		
		for(Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++){
			AliESDtrack* curTrack = (AliESDtrack*) fInputEvent->GetTrack(iTracks);
			if(!curTrack) continue;
			if(EsdTrackCuts->AcceptTrack(curTrack) ){
				if (fMCEvent){
                    if(fEventCuts->GetSignalRejection() != 0){
                        Int_t isFromMBHeader = fEventCuts->IsParticleFromBGEvent(TMath::Abs(curTrack->GetLabel()), fMCEvent, fESDEvent);
						if( (isFromMBHeader < 1) ) continue;
					}					
				}	
				fNumberOfESDTracks++;
			}
		}
		EsdTrackCuts->SetEtaRange(-1.4, -0.9);
		for(Int_t iTracks = 0; iTracks < fESDEvent->GetNumberOfTracks(); iTracks++){
			AliESDtrack* curTrack = fESDEvent->GetTrack(iTracks);
			if(!curTrack) continue;
			if(EsdTrackCuts->AcceptTrack(curTrack) ){
				if (fMCEvent){
                    if(fEventCuts->GetSignalRejection() != 0){
                        Int_t isFromMBHeader = fEventCuts->IsParticleFromBGEvent(TMath::Abs(curTrack->GetLabel()), fMCEvent, fESDEvent);
						if( (isFromMBHeader < 1) ) continue;
					}					
				}	
				fNumberOfESDTracks++;
			}	
		}
		delete EsdTrackCuts;
		EsdTrackCuts=0x0;
	} else if(fInputEvent->IsA()==AliAODEvent::Class()){
		for(Int_t iTracks = 0; iTracks<fInputEvent->GetNumberOfTracks(); iTracks++){
			AliAODTrack* curTrack = (AliAODTrack*) fInputEvent->GetTrack(iTracks);
            if(TMath::Abs(curTrack->Eta())<0.9 || TMath::Abs(curTrack->Eta())>1.4 ) continue;
			if(curTrack->Pt()<0.15) continue;
            if(TMath::Abs(curTrack->ZAtDCA())>5) continue;
			fNumberOfESDTracks++;
		}
	}
	
	return fNumberOfESDTracks;
}

//________________________________________________________________________
void AliAnalysisTaskResolution::Terminate(Option_t *)
{
	
}
