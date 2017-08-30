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

#include "AliAnalysisTaskMaterial.h"
#include "TChain.h"
#include "AliAnalysisManager.h"
#include "TParticle.h"
#include "TVectorF.h"
#include "AliPIDResponse.h"
#include "AliESDtrackCuts.h"
#include "TFile.h"

class iostream;

using namespace std;

ClassImp(AliAnalysisTaskMaterial)

AliAnalysisTaskMaterial::AliAnalysisTaskMaterial() : AliAnalysisTaskSE(),
	fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
	fConversionGammas(NULL),
	fConversionCuts(NULL),
	fEventCuts(NULL),
	fOutputList(NULL),
	fEventList(NULL),
	fRecGammaList(NULL),
	fAllMCGammaList(NULL),
	fAllMCConvGammaList(NULL),
	fTreeEvent(NULL),
	fTreeMaterialRec(NULL),
	fTreeMaterialAllGamma(NULL),
	fTreeMaterialConvGamma(NULL),
	fPrimVtxZ(0.),
	fNContrVtx(0),
	fNESDtracksEta09(0),
	fNESDtracksEta0914(0),
	fNESDtracksEta14(0),
	fGammaMCPt(0.),
	fGammaMCTheta(0.),
	fGammaMCConvPt(0.),
	fGammaMCConvTheta(0.),
	fMCConvCords(5),
	fMCConvDaughterProp(4),
	fGammaPt(0.),
	fGammaTheta(0.),
	fGammaChi2NDF(0.),
	fRecCords(5),
	fDaughterProp(4),
	fKind(0),      
	fIsHeavyIon(0),
	fIsMC(kFALSE),
	fESDEvent(NULL),
	fMCEvent(NULL)
{
   
}


//________________________________________________________________________
AliAnalysisTaskMaterial::AliAnalysisTaskMaterial(const char *name) : AliAnalysisTaskSE(name),
	fV0Reader(NULL),
    fV0ReaderName("V0ReaderV1"),
	fConversionGammas(NULL),
	fConversionCuts(NULL),
	fEventCuts(NULL),
	fOutputList(NULL),
	fEventList(NULL),
	fRecGammaList(NULL),
	fAllMCGammaList(NULL),
	fAllMCConvGammaList(NULL),
	fTreeEvent(NULL),
	fTreeMaterialRec(NULL),
	fTreeMaterialAllGamma(NULL),
	fTreeMaterialConvGamma(NULL),
	fPrimVtxZ(0.),
	fNContrVtx(0),
	fNESDtracksEta09(0),
	fNESDtracksEta0914(0),
	fNESDtracksEta14(0),
	fGammaMCPt(0.),
	fGammaMCTheta(0.),
	fGammaMCConvPt(0.),
	fGammaMCConvTheta(0.),
	fMCConvCords(5),
	fMCConvDaughterProp(4),
	fGammaPt(0.),
	fGammaTheta(0.),
	fGammaChi2NDF(0.),
	fRecCords(5),
	fDaughterProp(4),
	fKind(0),      
	fIsHeavyIon(0),
	fIsMC(kFALSE),
	fESDEvent(NULL),
	fMCEvent(NULL)
{
	// Default constructor
	
	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskMaterial::~AliAnalysisTaskMaterial()
{
	// default deconstructor
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskMaterial::Notify()
{
  if (fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){        
    fEventCuts->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
  } else if (fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
    fEventCuts->SetPeriodEnum(fV0Reader->GetPeriodName());
  } 
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskMaterial::UserCreateOutputObjects()
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
	
	fRecGammaList= new TList();
	fRecGammaList->SetName("RecGammaList");
	fRecGammaList->SetOwner(kTRUE);
	fOutputList->Add(fRecGammaList);
		
	fTreeMaterialRec = new TTree("ConvPointRec","ConvPointRec");   
	fTreeMaterialRec->Branch("nGoodTracksEta09",&fNESDtracksEta09,"fNESDtracksEta09/I");
	fTreeMaterialRec->Branch("recCords",&fRecCords);
	fTreeMaterialRec->Branch("daughterProp",&fDaughterProp);
	fTreeMaterialRec->Branch("pt",&fGammaPt,"fGammaPt/F");
	fTreeMaterialRec->Branch("theta",&fGammaTheta,"fGammaTheta/F");
	fTreeMaterialRec->Branch("chi2ndf",&fGammaChi2NDF,"fGammaChi2NDF/F");
	if (fIsMC) {
		fTreeMaterialRec->Branch("kind",&fKind,"fKind/b");
	}   
	fRecGammaList->Add(fTreeMaterialRec);
	
// 	if (fIsMC) {
// 		fAllMCGammaList = new TList();
// 		fAllMCGammaList->SetName("AllMCGammaList");
// 		fAllMCGammaList->SetOwner(kTRUE);
// 		fOutputList->Add(fAllMCGammaList);
// 		
// 		fTreeMaterialAllGamma = new TTree("AllGamma","AllGamma");   
// 		fTreeMaterialAllGamma->Branch("pt",&fGammaMCPt,"fGammaMCPt/F");
// 		fTreeMaterialAllGamma->Branch("theta",&fGammaMCTheta,"fGammaMCTheta/F");
// 		fAllMCGammaList->Add(fTreeMaterialAllGamma);
// 		
// 		fAllMCConvGammaList = new TList();
// 		fAllMCConvGammaList->SetName("AllMCGammaConvList");
// 		fAllMCConvGammaList->SetOwner(kTRUE);
// 		fOutputList->Add(fAllMCConvGammaList);
// 
// 	//       fMCConvCords = new Float_t[5];
// 	//       fMCConvDaughterProp = new Float_t[4];
// 
// 		
// 		fTreeMaterialConvGamma = new TTree("ConvGammaMC","ConvGammaMC");   
// 		fTreeMaterialConvGamma->Branch("Cords",&fMCConvCords);
// 		fTreeMaterialConvGamma->Branch("daughterProp",&fMCConvDaughterProp);
// 		fTreeMaterialConvGamma->Branch("Pt",&fGammaMCConvPt,"fGammaMCConvPt/F");
// 		fTreeMaterialConvGamma->Branch("Theta",&fGammaMCConvTheta,"fGammaMCConvTheta/F");   
// 		fAllMCConvGammaList->Add(fTreeMaterialConvGamma);
// 	}
    fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
	if(fV0Reader && fV0Reader->GetProduceV0FindingEfficiency())
		if (fV0Reader->GetV0FindingEfficiencyHistograms())
			fOutputList->Add(fV0Reader->GetV0FindingEfficiencyHistograms());


	
	PostData(1, fOutputList);
   
}

//________________________________________________________________________
void AliAnalysisTaskMaterial::UserExec(Option_t *){


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
// 	if(MCEvent()){
// 		ProcessMCPhotons();
// 	}
	PostData(1, fOutputList);
}

// ///________________________________________________________________________
// void AliAnalysisTaskMaterial::FillMCTree(Int_t eventPos){
// 	TParticle* candidate = (TParticle *)fMCEvent->Particle(eventPos);
// 	
// 	if(fConversionCuts->PhotonIsSelectedMC(candidate,fMCEvent,kFALSE)){
// 		fGammaMCPt = candidate->Pt();
// 		fGammaMCTheta = candidate->Theta();
// 		if (fTreeMaterialAllGamma){
// 			fTreeMaterialAllGamma->Fill();
// 		}	
// 	}
// 	if(fConversionCuts->PhotonIsSelectedMC(candidate,fMCEvent,kTRUE)){
// 		fGammaMCConvPt = candidate->Pt();
// 		fGammaMCConvTheta = candidate->Theta();
// 		TParticle* daughter1 = (TParticle *)fMCEvent->Particle(candidate->GetFirstDaughter());
// 		TParticle* daughter2 = (TParticle *)fMCEvent->Particle(candidate->GetLastDaughter());
// 		fMCConvCords(0) = (Float_t)daughter1->Vx();
// 		fMCConvCords(1) = (Float_t)daughter1->Vy();
// 		fMCConvCords(2) = (Float_t)daughter1->Vz();
// 		fMCConvCords(3) = (Float_t)daughter1->R();
// 		fMCConvCords(4) = (Float_t)daughter1->Phi();
// 		
// 		fMCConvDaughterProp(0) = (Float_t)daughter1->Pt();
// 		fMCConvDaughterProp(1) = (Float_t)daughter1->Theta();
// 		fMCConvDaughterProp(2) = (Float_t)daughter2->Pt();
// 		fMCConvDaughterProp(3) = (Float_t)daughter2->Theta();      
// 			
// 		if (fTreeMaterialConvGamma){
// 			fTreeMaterialConvGamma->Fill();
// 		}
// 	} // Converted MC Gamma		
// }
// 
// ///________________________________________________________________________
// void AliAnalysisTaskMaterial::ProcessMCPhotons(){
// 	// Loop over all primary MC particle
// 	for(Int_t i = 0; i < fMCEvent->GetNumberOfPrimaries(); i++) {
// 		TParticle* particle = (TParticle *)fMCEvent->Particle(i);
// 		if (!particle) continue;
// 		Int_t isMCFromMBHeader = -1;
// 		if(fConversionCuts->GetSignalRejection() != 0){
// 			isMCFromMBHeader
// 				= fConversionCuts->IsParticleFromBGEvent(i, fMCEvent, fESDEvent);
// 			if(isMCFromMBHeader == 0) continue;
// 		}		
// 		if (particle->GetPdgCode() == 111 && particle->GetFirstDaughter() >= fMCEvent->GetNumberOfPrimaries()){
// // 			cout << "Undecayed pi0 found with mother: " << particle->GetMother(0) << endl;
// 			for (Int_t j = 0; j < 2 ; j++){
//  				FillMCTree(particle->GetDaughter(j));
// 			}
// 		} else {
// 			FillMCTree(i);
// 		}	
// 	}	
// }

///________________________________________________________________________
void AliAnalysisTaskMaterial::ProcessPhotons(){

		// Fill Histograms for QA and MC
	for(Int_t firstGammaIndex=0;firstGammaIndex<fConversionGammas->GetEntriesFast();firstGammaIndex++){
		AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(firstGammaIndex));
		if (gamma ==NULL) continue;
		if(!fConversionCuts->PhotonIsSelected(gamma,fESDEvent)) continue;

		fNESDtracksEta09 = CountTracks09();
		fGammaPt = gamma->GetPhotonPt();
		fGammaTheta = gamma->GetPhotonTheta();
		fGammaChi2NDF = gamma->GetChi2perNDF();
		fRecCords(0) = (Float_t)gamma->GetConversionX();
		fRecCords(1) = (Float_t)gamma->GetConversionY();
		fRecCords(2) = (Float_t)gamma->GetConversionZ();
		fRecCords(3) = (Float_t)gamma->GetConversionRadius();
		fRecCords(4) = (Float_t)gamma->GetPhotonPhi();
		
		AliESDtrack * negTrack = fConversionCuts->GetESDTrack(fESDEvent, gamma->GetTrackLabelNegative());
		AliESDtrack * posTrack = fConversionCuts->GetESDTrack(fESDEvent, gamma->GetTrackLabelPositive());
		fDaughterProp(0) = (Float_t)posTrack->Pt();
		fDaughterProp(1) = (Float_t)posTrack->Theta();
		fDaughterProp(2) = (Float_t)negTrack->Pt();
		fDaughterProp(3) = (Float_t)negTrack->Theta();

		fKind = 9;	
		
		if(fMCEvent){
			
// 			cout << "generating MC stack"<< endl;
			
			const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
			Double_t mcProdVtxX 	= primVtxMC->GetX();
			Double_t mcProdVtxY 	= primVtxMC->GetY();
			Double_t mcProdVtxZ 	= primVtxMC->GetZ();

            TParticle *posDaughter = gamma->GetPositiveMCDaughter(fMCEvent);
            TParticle *negDaughter = gamma->GetNegativeMCDaughter(fMCEvent);
// 			cout << "generate Daughters: "<<posDaughter << "\t" << negDaughter << endl;
		
            if(fMCEvent && fEventCuts->GetSignalRejection() != 0){
				Int_t isPosFromMBHeader
                = fEventCuts->IsParticleFromBGEvent(gamma->GetMCLabelPositive(), fMCEvent, fESDEvent);
				Int_t isNegFromMBHeader
                = fEventCuts->IsParticleFromBGEvent(gamma->GetMCLabelNegative(), fMCEvent, fESDEvent);
				if( (isNegFromMBHeader < 1) || (isPosFromMBHeader < 1)) continue;
			}
				
			if(posDaughter == NULL || negDaughter == NULL){ 
				fKind = 9; // garbage
// 				cout << "one of the daughters not available" << endl;
			} else if(posDaughter->GetMother(0) != negDaughter->GetMother(0) || (posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1)){ 
				// Not Same Mother == Combinatorial Bck
				fKind = 1;
// 				cout << "not the same mother" << endl;
				Int_t pdgCodePos; 
				if (posDaughter->GetPdgCode()) pdgCodePos = posDaughter->GetPdgCode(); else continue;
				Int_t pdgCodeNeg; 
				if (negDaughter->GetPdgCode()) pdgCodeNeg = negDaughter->GetPdgCode(); else continue;
// 				cout << "PDG codes daughters: " << pdgCodePos << "\t" << pdgCodeNeg << endl;
				if(TMath::Abs(pdgCodePos)==11 && TMath::Abs(pdgCodeNeg)==11)
					fKind = 10; //Electron Combinatorial
				if(TMath::Abs(pdgCodePos)==11 && TMath::Abs(pdgCodeNeg)==11 && (posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1))
					fKind = 15; //direct Electron Combinatorial
				
				if(TMath::Abs(pdgCodePos)==211 && TMath::Abs(pdgCodeNeg)==211)
					fKind = 11; //Pion Combinatorial
				if((TMath::Abs(pdgCodePos)==211 && TMath::Abs(pdgCodeNeg)==2212) ||
					(TMath::Abs(pdgCodePos)==2212 && TMath::Abs(pdgCodeNeg)==211))
					fKind = 12; //Pion, Proton Combinatorics
				if((TMath::Abs(pdgCodePos)==211 && TMath::Abs(pdgCodeNeg)==11) ||
					(TMath::Abs(pdgCodePos)==11 && TMath::Abs(pdgCodeNeg)==211))
					fKind = 13; //Pion, Electron Combinatorics
				if (TMath::Abs(pdgCodePos)==321 || TMath::Abs(pdgCodeNeg)==321)	
					fKind = 14; //Kaon combinatorics

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
					fKind = 2; // combinatorics from hadronic decays
				else if ( !(pdgCodeNeg==pdgCodePos)){
                    TParticle *truePhotonCanditate = gamma->GetMCParticle(fMCEvent);
					Int_t motherLabelPhoton = truePhotonCanditate->GetMother(0);
                    Bool_t gammaIsPrimary = fEventCuts->IsConversionPrimaryESD( fMCEvent, posDaughter->GetMother(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
					if(pdgCode == 111) 
						fKind = 3; // pi0 Dalitz
					else if (pdgCode == 221) 
						fKind = 4; // eta Dalitz
					else if (!(negDaughter->GetUniqueID() != 5 || posDaughter->GetUniqueID() !=5)){
						if(pdgCode == 22 && gammaIsPrimary){
							fKind = 0; // primary photons
						} else if (pdgCode == 22){
							fKind = 5; //secondary photons
						}		
					} else 	fKind = 9; //garbage
				} else fKind = 9; //garbage
			}					
		}
				
		if (fTreeMaterialRec){
			fTreeMaterialRec->Fill();
		}
	}
}

//________________________________________________________________________
Int_t AliAnalysisTaskMaterial::CountTracks09(){
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
	} 
	return fNumberOfESDTracks;
}

//________________________________________________________________________
Int_t AliAnalysisTaskMaterial::CountTracks0914(){
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
	} 
	
	return fNumberOfESDTracks;
}

//________________________________________________________________________
void AliAnalysisTaskMaterial::Terminate(Option_t *)
{
//    if (fStreamMaterial){
//       fStreamMaterial->GetFile()->Write();
//    }
}
