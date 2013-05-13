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

//________________________________________________________________________
AliAnalysisTaskMaterial::AliAnalysisTaskMaterial(const char *name) : AliAnalysisTaskSE(name),
   fV0Reader(NULL),
   fConversionGammas(NULL),
   fConversionCuts(NULL),
   fStreamMaterial(NULL),
   fStreamResolution(NULL),
   fIsHeavyIon(kFALSE),
   fOutputList(NULL),
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
   if(fStreamMaterial){
      delete fStreamMaterial;
      fStreamMaterial = 0x0;
   }
   if(fStreamResolution){
      delete fStreamResolution;
      fStreamResolution = 0x0;
   }
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

   // V0 Reader Cuts
   TString cutnumber = fConversionCuts->GetCutNumber();

   fStreamMaterial = new TTreeSRedirector(Form("GammaConvV1_Material_%s.root",cutnumber.Data()),"recreate");
	 fStreamResolution = new TTreeSRedirector(Form("GammaConvV1_Resolution_%s.root",cutnumber.Data()),"recreate");
   PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskMaterial::UserExec(Option_t *){

   fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");

   Int_t eventQuality = ((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetEventQuality();
	if(eventQuality != 0){// Event Not Accepted
      return;
   }
   fESDEvent = (AliESDEvent*) InputEvent();
   if (fESDEvent==NULL) return;
   if(fIsHeavyIon && !fConversionCuts->IsCentralitySelected(fESDEvent)) return;
	Int_t nESDtracksEta09 = CountESDTracks09(); // Estimate Event Multiplicity
	Int_t nESDtracksEta0914 = CountESDTracks0914(); // Estimate Event Multiplicity
	//	Int_t nESDtracksEta14 = CountESDTracks14(); // Estimate Event Multiplicity
	Int_t nESDtracksEta14; // Estimate Event Multiplicity
  nESDtracksEta14= nESDtracksEta09 + nESDtracksEta0914;
	Int_t nContrVtx;
	if(fESDEvent){
		if(fESDEvent->GetPrimaryVertexTracks()->GetNContributors()>0) {
			nContrVtx = fESDEvent->GetPrimaryVertexTracks()->GetNContributors();
		} else if(fESDEvent->GetPrimaryVertexTracks()->GetNContributors()<1) {
			if(fESDEvent->GetPrimaryVertexSPD()->GetNContributors()>0) {
				nContrVtx = fESDEvent->GetPrimaryVertexSPD()->GetNContributors();

			} else if(fESDEvent->GetPrimaryVertexSPD()->GetNContributors()<1) {
				nContrVtx = 0;
			}
		}
   }
	Float_t primVtxZ = fESDEvent->GetPrimaryVertex()->GetZ();
	
	if (fStreamMaterial){
		(*fStreamMaterial)<<"Event"
						<< "primVtxZ=" << primVtxZ
						<< "nContrVtx=" << nContrVtx
						<< "nGoodTracksEta09=" << nESDtracksEta09
						<< "nGoodTracksEta0914=" << nESDtracksEta0914
						<< "nGoodTracksEta14=" << nESDtracksEta14
						<< "\n";
	}
	
   fConversionGammas=fV0Reader->GetReconstructedGammas();
	if(MCEvent()){
		fMCEvent = MCEvent();
	}
   ProcessPhotons();
	if(MCEvent()){
		ProcessMCPhotons();
	}
   PostData(1, fOutputList);
}

///________________________________________________________________________
void AliAnalysisTaskMaterial::FillMCTree(Int_t stackPos){
	AliStack *MCStack = fMCEvent->Stack();
	TParticle* candidate = (TParticle *)MCStack->Particle(stackPos);
	if(fConversionCuts->PhotonIsSelectedMC(candidate,MCStack,kFALSE)){
		Float_t gammaPt = candidate->Pt();
		Float_t gammaTheta = candidate->Theta();
		if (fStreamMaterial){
			(*fStreamMaterial)<<"AllGamma"
							<< "pt=" << gammaPt
// 							<< "p=" << gammaP
							<< "theta=" << gammaTheta
							<< "\n";
		}	
	}
	if(fConversionCuts->PhotonIsSelectedMC(candidate,MCStack,kTRUE)){
		Float_t gammaPt = candidate->Pt();
		Float_t gammaTheta = candidate->Theta();
		TParticle* daughter1 = (TParticle *)MCStack->Particle(candidate->GetFirstDaughter()); 
		TParticle* daughter2 = (TParticle *)MCStack->Particle(candidate->GetLastDaughter()); 
		TVectorF coord(5);	
		coord(0) = daughter1->Vx();
		coord(1) = daughter1->Vy();
		coord(2) = daughter1->Vz();
		coord(3) = daughter1->R();
		coord(4) = candidate->Phi();
		TVectorF daughterProp(4);	
		if (daughter1->	GetPdgCode() < 0){
			daughterProp(0) = daughter2->Pt();
			daughterProp(1) = daughter2->Theta();
			daughterProp(2) = daughter1->Pt();
			daughterProp(3) = daughter1->Theta();
		} else {
			daughterProp(0) = daughter1->Pt();
			daughterProp(1) = daughter1->Theta();
			daughterProp(2) = daughter2->Pt();
			daughterProp(3) = daughter2->Theta();
		}
		if (fStreamMaterial){
			(*fStreamMaterial)<<"ConvGammaMC"
							<< "pt=" << gammaPt
// 							<< "p=" << gammaP
 							<< "theta=" << gammaTheta
							<< "coord.=" << &coord
							<< "daughterProp.=" << &daughterProp
							<< "\n";
		}
	} // Converted MC Gamma		
}

///________________________________________________________________________
void AliAnalysisTaskMaterial::ProcessMCPhotons(){
	// Loop over all primary MC particle
	AliStack *ffMCStack = fMCEvent->Stack();
	for(Int_t i = 0; i < ffMCStack->GetNprimary(); i++) {
		TParticle* particle = (TParticle *)ffMCStack->Particle(i);
		if (!particle) continue;		
		if (particle->GetPdgCode() == 111 && particle->GetFirstDaughter() >= ffMCStack->GetNprimary()){
// 			cout << "Undecayed pi0 found with mother: " << particle->GetMother(0) << endl;
			for (Int_t j = 0; j < 2 ; j++){
 				FillMCTree(particle->GetDaughter(j));
			}
		} else {
			FillMCTree(i);
		}	
	}	
}

///________________________________________________________________________
void AliAnalysisTaskMaterial::ProcessPhotons(){

	// Fill Histograms for QA and MC
   for(Int_t firstGammaIndex=0;firstGammaIndex<fConversionGammas->GetEntriesFast();firstGammaIndex++){
      AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(firstGammaIndex));
      if (gamma ==NULL) continue;
      if(!fConversionCuts->PhotonIsSelected(gamma,fESDEvent)) continue;
// 		cout << "i=  " <<firstGammaIndex << " of "<< fConversionGammas->GetEntriesFast() << endl;
      Float_t gammaPt = gamma->GetPhotonPt();
// 		Float_t gammaP = gamma->GetPhotonP();
		Float_t gammaTheta = gamma->GetPhotonTheta();
		Float_t gammaEta = gamma->GetPhotonEta();
		Float_t gammaChi2NDF = gamma->GetChi2perNDF();
// 		Float_t gammaX = gamma->GetConversionX();
// 		Float_t gammaY = gamma->GetConversionY();
 		Float_t gammaZ = gamma->GetConversionZ();
      Float_t gammaR = gamma->GetConversionRadius();
 		Float_t gammaPhi = gamma->GetPhotonPhi();
		
		TVectorF coord(5);	
		coord(0) = gamma->GetConversionX();
		coord(1) = gamma->GetConversionY();
		coord(2) = gamma->GetConversionZ();
		coord(3) = gamma->GetConversionRadius();
		coord(4) = gamma->GetPhotonPhi();
		TVectorF daughterProp(4);	
		AliESDtrack * negTrack = fConversionCuts->GetESDTrack(fESDEvent, gamma->GetTrackLabelNegative());
      AliESDtrack * posTrack = fConversionCuts->GetESDTrack(fESDEvent, gamma->GetTrackLabelPositive());
		daughterProp(0) = posTrack->Pt();
      daughterProp(1) = posTrack->Theta();
      daughterProp(2) = negTrack->Pt();
      daughterProp(3) = negTrack->Theta();
				
		UInt_t kind = 9;
		if(fMCEvent){
// 			cout << "generating MC stack"<< endl;
			AliStack *fMCStack = fMCEvent->Stack();
			if (!fMCStack) continue;
			TParticle *posDaughter = gamma->GetPositiveMCDaughter(fMCStack);
			TParticle *negDaughter = gamma->GetNegativeMCDaughter(fMCStack);
// 			cout << "generate Daughters: "<<posDaughter << "\t" << negDaughter << endl;
			
			if(posDaughter == NULL || negDaughter == NULL){ 
				kind = 9; // garbage
// 				cout << "one of the daughters not available" << endl;
			} else if(posDaughter->GetMother(0) != negDaughter->GetMother(0) || (posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1)){ 
				// Not Same Mother == Combinatorial Bck
				kind = 1;
// 				cout << "not the same mother" << endl;
				Int_t pdgCodePos; 
				if (posDaughter->GetPdgCode()) pdgCodePos = posDaughter->GetPdgCode(); else continue;
				Int_t pdgCodeNeg; 
				if (negDaughter->GetPdgCode()) pdgCodeNeg = negDaughter->GetPdgCode(); else continue;
// 				cout << "PDG codes daughters: " << pdgCodePos << "\t" << pdgCodeNeg << endl;
				if(TMath::Abs(pdgCodePos)==11 && TMath::Abs(pdgCodeNeg)==11)
					kind = 10; //Electron Combinatorial
				if(TMath::Abs(pdgCodePos)==11 && TMath::Abs(pdgCodeNeg)==11 && (posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1))
					kind = 15; //direct Electron Combinatorial
				
				if(TMath::Abs(pdgCodePos)==211 && TMath::Abs(pdgCodeNeg)==211)
					kind = 11; //Pion Combinatorial
				if((TMath::Abs(pdgCodePos)==211 && TMath::Abs(pdgCodeNeg)==2212) ||
					(TMath::Abs(pdgCodePos)==2212 && TMath::Abs(pdgCodeNeg)==211))
					kind = 12; //Pion, Proton Combinatorics
				if((TMath::Abs(pdgCodePos)==211 && TMath::Abs(pdgCodeNeg)==11) ||
					(TMath::Abs(pdgCodePos)==11 && TMath::Abs(pdgCodeNeg)==211))
					kind = 13; //Pion, Electron Combinatorics
				if (TMath::Abs(pdgCodePos)==321 || TMath::Abs(pdgCodeNeg)==321)	
					kind = 14; //Kaon combinatorics

			} else {
// 				cout << "same mother" << endl;
				Int_t pdgCodePos; 
				if (posDaughter->GetPdgCode()) pdgCodePos = posDaughter->GetPdgCode(); else continue;
				Int_t pdgCodeNeg; 
				if (negDaughter->GetPdgCode()) pdgCodeNeg = negDaughter->GetPdgCode(); else continue;
// 				cout << "PDG codes daughters: " << pdgCodePos << "\t" << pdgCodeNeg << endl;
				Int_t pdgCode; 
				if (gamma->GetMCParticle(fMCStack)->GetPdgCode()) pdgCode = gamma->GetMCParticle(fMCStack)->GetPdgCode(); else continue;
// 				cout << "PDG code: " << pdgCode << endl;
				if(TMath::Abs(pdgCodePos)!=11 || TMath::Abs(pdgCodeNeg)!=11)
					kind = 2; // combinatorics from hadronic decays
				else if ( !(pdgCodeNeg==pdgCodePos)){
					TParticle *truePhotonCanditate = gamma->GetMCParticle(fMCStack);
					Int_t motherLabelPhoton = truePhotonCanditate->GetMother(0);
					if(pdgCode == 111) 
						kind = 3; // pi0 Dalitz
					else if (pdgCode == 221) 
						kind = 4; // eta Dalitz
					else if (!(negDaughter->GetUniqueID() != 5 || posDaughter->GetUniqueID() !=5)){
						if(pdgCode == 22 && motherLabelPhoton < fMCStack->GetNprimary()){
							kind = 0; // primary photons
						} else if (pdgCode == 22){
							kind = 5; //secondary photons
						}
						if(pdgCode == 22){
							Float_t mcPt   = truePhotonCanditate->Pt();
							Float_t mcR = gamma->GetNegativeMCDaughter(fMCStack)->R();
							Float_t mcZ = gamma->GetNegativeMCDaughter(fMCStack)->Vz();
							Float_t mcPhi = gamma->GetNegativeMCDaughter(fMCStack)->Phi();
							Float_t mcEta = gamma->GetNegativeMCDaughter(fMCStack)->Eta();
							if (fStreamResolution){
								(*fStreamResolution)<<"Resolution"
								<< "ESDpt=" << gammaPt
								<< "ESDphi=" << gammaPhi
								<< "ESDeta=" << gammaEta
								<< "ESDr="<< gammaR
								<< "ESDz="<< gammaZ
								<< "MCpt=" << mcPt
								<< "MCphi=" << mcPhi
								<< "MCeta=" << mcEta
								<< "MCr="<< mcR
								<< "MCz="<< mcZ
								<< "chi2ndf=" << gammaChi2NDF
								<< "\n";
							}
						}		
					} else 	kind = 9; //garbage
				} else kind = 9; //garbage
			}
// 			cout << gammaPt << "\t" << gammaP<< "\t" << gammaEta<<  "\t" <<gammaChi2NDF << "\t" << gammaX<<  "\t" <<gammaY << "\t" << gammaZ<< "\t" << gammaR<<  "\t" <<  gammaPhi<< "\t" <<kind << endl;
		
			if (fStreamMaterial){
				(*fStreamMaterial)<<"ConvPointRec"
								<< "pt=" << gammaPt
								<< "theta=" << gammaTheta
								<< "chi2ndf=" << gammaChi2NDF
								<< "kind=" << kind
								<< "coord.=" << &coord
								<< "daughterProp.=" << &daughterProp
								<< "\n";
			}		
		} else {
				if (fStreamMaterial){
				(*fStreamMaterial)<<"ConvPointRec"
								<< "pt=" << gammaPt
								<< "theta=" << gammaTheta
								<< "chi2ndf=" << gammaChi2NDF
								<< "coord.=" << &coord
								<< "daughterProp.=" << &daughterProp
								<< "\n";
			}
		}
	}
}

//________________________________________________________________________
Int_t AliAnalysisTaskMaterial::CountESDTracks09(){
   
   // Using standard function for setting Cuts
   Bool_t selectPrimaries=kTRUE;
   AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
   EsdTrackCuts->SetMaxDCAToVertexZ(2);
   EsdTrackCuts->SetEtaRange(-0.9, 0.9);
   EsdTrackCuts->SetPtRange(0.15);

   Int_t fNumberOfESDTracks = 0;
   for(Int_t iTracks = 0; iTracks < fESDEvent->GetNumberOfTracks(); iTracks++){
      AliESDtrack* curTrack = fESDEvent->GetTrack(iTracks);
      if(!curTrack) continue;
      if(EsdTrackCuts->AcceptTrack(curTrack) ) fNumberOfESDTracks++;
   }
   delete EsdTrackCuts;
   EsdTrackCuts=0x0;

   return fNumberOfESDTracks;
}

Int_t AliAnalysisTaskMaterial::CountESDTracks0914(){

   // Using standard function for setting Cuts ; We use TPCOnlyTracks for outer eta region
   //Bool_t selectPrimaries=kTRUE;
	 //   EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
   AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
	 EsdTrackCuts->SetMaxDCAToVertexXY(5);
	 //	 EsdTrackCuts->SetMaxDCAToVertexXYPtDep("sqrt(0.15^2+(0.4/pt)^2");
   EsdTrackCuts->SetEtaRange(0.9, 1.4);
   EsdTrackCuts->SetPtRange(0.15);

   Int_t fNumberOfESDTracks = 0;
   for(Int_t iTracks = 0; iTracks < fESDEvent->GetNumberOfTracks(); iTracks++){
      AliESDtrack* curTrack = fESDEvent->GetTrack(iTracks);
      if(!curTrack) continue;
      if(EsdTrackCuts->AcceptTrack(curTrack) ) fNumberOfESDTracks++;
   }
   EsdTrackCuts->SetEtaRange(-1.4, -0.9);
   for(Int_t iTracks = 0; iTracks < fESDEvent->GetNumberOfTracks(); iTracks++){
      AliESDtrack* curTrack = fESDEvent->GetTrack(iTracks);
      if(!curTrack) continue;
      if(EsdTrackCuts->AcceptTrack(curTrack) ) fNumberOfESDTracks++;
   }
   delete EsdTrackCuts;
   EsdTrackCuts=0x0;

   return fNumberOfESDTracks;
}

Int_t AliAnalysisTaskMaterial::CountESDTracks14(){

   // Using standard function for setting Cuts
   Bool_t selectPrimaries=kTRUE;
   AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
   EsdTrackCuts->SetMaxDCAToVertexZ(2);
   EsdTrackCuts->SetEtaRange(-1.4, 1.4);
   EsdTrackCuts->SetPtRange(0.15);

   Int_t fNumberOfESDTracks = 0;
   for(Int_t iTracks = 0; iTracks < fESDEvent->GetNumberOfTracks(); iTracks++){
      AliESDtrack* curTrack = fESDEvent->GetTrack(iTracks);
      if(!curTrack) continue;
      if(EsdTrackCuts->AcceptTrack(curTrack) ) fNumberOfESDTracks++;
   }
   delete EsdTrackCuts;
   EsdTrackCuts=0x0;

   return fNumberOfESDTracks;
}


//________________________________________________________________________
void AliAnalysisTaskMaterial::Terminate(Option_t *)
{
//    if (fStreamMaterial){
//       fStreamMaterial->GetFile()->Write();
//    }
//    if (fStreamResolution){
//       fStreamResolution->GetFile()->Write();
//    }
}
