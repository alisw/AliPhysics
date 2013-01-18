/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
//
//
//             Base class for DStar - Hadron Correlations Analysis
//
//-----------------------------------------------------------------------
//          
//
//						   Author S.Bjelogrlic
//                         Utrecht University 
//                      sandro.bjelogrlic@cern.ch
//
//-----------------------------------------------------------------------

/* $Id$ */

//#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TVector3.h>
#include <TChain.h>
#include "TROOT.h"

#include "AliAnalysisTaskDStarCorrelations.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODPidHF.h"
#include "AliVParticle.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliNormalizationCounter.h"
#include "AliReducedParticle.h"
#include "AliHFCorrelator.h"
#include "AliAODMCHeader.h"
#include "AliEventPoolManager.h"

using std::cout;
using std::endl;


ClassImp(AliAnalysisTaskDStarCorrelations)


//__________________________________________________________________________
AliAnalysisTaskDStarCorrelations::AliAnalysisTaskDStarCorrelations() :
AliAnalysisTaskSE(),
fhandler(0x0),
fmcArray(0x0),
fCounter(0x0),
fCorrelator(0x0),
fselect(0),
fmontecarlo(kFALSE),
fmixing(kFALSE),
fSystem(kFALSE),
fReco(kTRUE),
fEvents(0),
fDebugLevel(0),
fDisplacement(0),

fOutput(0x0), 
fOutputMC(0x0),
fCuts(0),
fCuts2(0)
{
// default constructor	
}

//__________________________________________________________________________
AliAnalysisTaskDStarCorrelations::AliAnalysisTaskDStarCorrelations(const Char_t* name,AliRDHFCutsDStartoKpipi* cuts, AliHFAssociatedTrackCuts *AsscCuts) :
AliAnalysisTaskSE(name),

fhandler(0x0),
fmcArray(0x0),
fCounter(0x0),
fCorrelator(0x0),
fselect(0),
fmontecarlo(kFALSE),
fmixing(kFALSE),
fSystem(kFALSE),
fReco(kTRUE),
fEvents(0),
fDebugLevel(0),
fDisplacement(0),

fOutput(0x0),   
fOutputMC(0x0), 
fCuts(0),
fCuts2(AsscCuts)
{
	fCuts=cuts;
	Info("AliAnalysisTaskDStarCorrelations","Calling Constructor");
	DefineInput(0, TChain::Class());
	
	DefineOutput(1,TList::Class()); // histos from data and MC
	DefineOutput(2,TList::Class()); // histos from MC
	DefineOutput(3,AliRDHFCutsDStartoKpipi::Class()); // my D meson cuts
	DefineOutput(4,AliHFAssociatedTrackCuts::Class()); // my associated tracks cuts
	DefineOutput(5,AliNormalizationCounter::Class());   // normalization
}

//__________________________________________________________________________

AliAnalysisTaskDStarCorrelations::~AliAnalysisTaskDStarCorrelations() {
	//
	// destructor
	//
	
	Info("AliAnalysisTaskDStarCorrelations","Calling Destructor");  
	
	if(fhandler) {delete fhandler; fhandler = 0;}
	//if(fPoolMgr) {delete fPoolMgr; fPoolMgr = 0;}    
	if(fmcArray) {delete fmcArray; fmcArray = 0;}
	if(fCounter) {delete fCounter; fCounter = 0;}
	if(fCorrelator) {delete fCorrelator; fCorrelator = 0;}
	if(fOutput) {delete fOutput; fOutput = 0;}
	if(fOutputMC) {delete fOutputMC; fOutputMC = 0;}
	if(fCuts) {delete fCuts; fCuts = 0;}
	if(fCuts2) {delete fCuts2; fCuts2=0;}

}

//___________________________________________________________
void AliAnalysisTaskDStarCorrelations::Init(){
	//
	// Initialization
	//
	if(fDebugLevel > 1) printf("AliAnalysisTaskDStarCorrelations::Init() \n");
	
	AliRDHFCutsDStartoKpipi* copyfCuts=new AliRDHFCutsDStartoKpipi(*fCuts);
	
	
	
	
	// Post the D* cuts
	PostData(3,copyfCuts);
	
	// Post the hadron cuts
	PostData(4,fCuts2);
	

	
	return;
}


//_________________________________________________
void AliAnalysisTaskDStarCorrelations::UserCreateOutputObjects(){
	Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
	
	//slot #1  
	//OpenFile(0);
	fOutput = new TList();
	fOutput->SetOwner();
	
	fOutputMC = new TList();
	fOutputMC->SetOwner();
	
	// define histograms
	DefineHistoForAnalysis();
	fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(5)->GetContainer()->GetName()));
	fCounter->Init();
	
    Double_t Pi = TMath::Pi();
	fCorrelator = new AliHFCorrelator("Correlator",fCuts2,fSystem); // fCuts2 is the hadron cut object, fSystem to switch between pp or PbPb
	fCorrelator->SetDeltaPhiInterval((-0.5-1./32)*Pi,(1.5-1./32)*Pi); // set correct phi interval
	fCorrelator->SetEventMixing(fmixing); //set kFALSE/kTRUE for mixing Off/On
	fCorrelator->SetAssociatedParticleType(fselect); // set 1/2/3 for hadron/kaons/kzeros
	fCorrelator->SetApplyDisplacementCut(fDisplacement); //set kFALSE/kTRUE for using the displacement cut
	fCorrelator->SetUseMC(fmontecarlo);
	fCorrelator->SetUseReco(fReco);
	Bool_t pooldef = fCorrelator->DefineEventPool();
	
	if(!pooldef) AliInfo("Warning:: Event pool not defined properly");


	
	PostData(1,fOutput); // set the outputs
	PostData(2,fOutputMC); // set the outputs
	PostData(5,fCounter); // set the outputs
}
//_________________________________________________
void AliAnalysisTaskDStarCorrelations::UserExec(Option_t *){

	
	if(fDebugLevel){
		
		if(fReco) std::cout << "USING RECONSTRUCTION" << std::endl;
		if(!fReco) std::cout << "USING MC TRUTH" << std::endl;
        std::cout << " " << std::endl;
        std::cout << "=================================================================================" << std::endl;
        if(!fmixing){
            if(fselect==1) std::cout << "TASK::Correlation with hadrons on SE "<< std::endl;
            if(fselect==2) std::cout << "TASK::Correlation with kaons on SE "<< std::endl;
            if(fselect==3) std::cout << "TASK::Correlation with kzeros on SE "<< std::endl;
        }
        if(fmixing){
            if(fselect==1) std::cout << "TASK::Correlation with hadrons on ME "<< std::endl;
            if(fselect==2) std::cout << "TASK::Correlation with kaons on ME "<< std::endl;
            if(fselect==3) std::cout << "TASK::Correlation with kzeros on ME "<< std::endl;
        }
        
    }// end if debug
    
	if (!fInputEvent) {
		Error("UserExec","NO EVENT FOUND!");
		return;
	}
	
	AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
	if(!aodEvent){
	  AliError("AOD event not found!");
	  return;
	}
	
	
	
	fEvents++; // event counter
	((TH1D*)fOutput->FindObject("NofEvents"))->Fill(0);
    fCounter->StoreEvent(aodEvent,fCuts,fmontecarlo);
  
	// load MC array
	fmcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	if(fmontecarlo && !fmcArray){
	  AliError("Array of MC particles not found");
	  return;
	}
	

	
	
	Bool_t isEvSel=fCuts->IsEventSelected(aodEvent);
	if(!isEvSel) return;
	
    fCorrelator->SetAODEvent(aodEvent); // set the event to be processed
    
    ((TH1D*)fOutput->FindObject("NofEvents"))->Fill(1);
	//
	Bool_t correlatorON = fCorrelator->Initialize(); //define the pool for mixing
	if(!correlatorON) {
		AliInfo("AliHFCorrelator didn't initialize the pool correctly or processed a bad event");
		return;
	}
	((TH1D*)fOutput->FindObject("NofEvents"))->Fill(2);
	
	if(fmontecarlo) fCorrelator->SetMCArray(fmcArray);
	
	
	// check the event type
	// load MC header
	
	if(fmontecarlo){
		AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
		if (fmontecarlo && !mcHeader) {
			AliError("Could not find MC Header in AOD");
			return;
		}
	
		Bool_t isMCeventgood = kFALSE;
       
		
		Int_t eventType = mcHeader->GetEventType();
		Int_t NMCevents = fCuts2->GetNofMCEventType();
		
		for(Int_t k=0; k<NMCevents; k++){
			Int_t * MCEventType = fCuts2->GetMCEventType();
			
			if(eventType == MCEventType[k]) isMCeventgood= kTRUE;
			((TH1D*)fOutputMC->FindObject("EventTypeMC"))->Fill(eventType);
		}
		
		if(NMCevents && !isMCeventgood){
		if(fDebugLevel)	std::cout << "The MC event " << eventType << " not interesting for this analysis: skipping" << std::endl;
			return;	
		}
		
	} // end if montecarlo
	
	// D* reconstruction
	TClonesArray *arrayDStartoD0pi=0;
	if(!aodEvent && AODEvent() && IsStandardAOD()) {
		// In case there is an AOD handler writing a standard AOD, use the AOD 
		// event in memory rather than the input (ESD) event.    
		aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
		// in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
		// have to taken from the AOD event hold by the AliAODExtension
		AliAODHandler* aodHandler = (AliAODHandler*) 
		((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
		if(aodHandler->GetExtensions()) {
			AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
			AliAODEvent *aodFromExt = ext->GetAOD();
			arrayDStartoD0pi=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
		}
	} else {
		arrayDStartoD0pi=(TClonesArray*)aodEvent->GetList()->FindObject("Dstar");
	}
	
	if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;
	
	// initialize variables you will need for the D*
	
	Double_t ptDStar;//
	Double_t phiDStar;//
	Double_t etaDStar;//
	Bool_t isInPeak, isInSideBand, isDStarMCtag;
	Double_t invMassDZero;
	Double_t deltainvMDStar;

	
	Double_t mPDGD0=1.8648;//TDatabasePDG::Instance()->GetParticle(421)->Mass();
	Double_t mPDGDstar=2.01022;//TDatabasePDG::Instance()->GetParticle(413)->Mass();
		
	
	//MC tagging for DStar
	//D* and D0 prongs needed to MatchToMC method
	Int_t pdgDgDStartoD0pi[2]={421,211};
	Int_t pdgDgD0toKpi[2]={321,211};

	Bool_t isDStarCand = kFALSE;
	//loop on D* candidates
	for (Int_t iDStartoD0pi = 0; iDStartoD0pi<arrayDStartoD0pi->GetEntriesFast(); iDStartoD0pi++) {
		isInPeak = kFALSE;
		isInSideBand = kFALSE;
		isDStarMCtag = kFALSE;
		ptDStar = -123.4;
		phiDStar = -999;
		etaDStar = -56.;
		invMassDZero = - 999;
		deltainvMDStar = -998;
		
		AliAODRecoCascadeHF* dstarD0pi = (AliAODRecoCascadeHF*)arrayDStartoD0pi->At(iDStartoD0pi);
		if(!dstarD0pi->GetSecondaryVtx()) continue;
		AliAODRecoDecayHF2Prong* theD0particle = (AliAODRecoDecayHF2Prong*)dstarD0pi->Get2Prong();
		if (!theD0particle) continue;

		
		// track quality cuts
		Int_t isTkSelected = fCuts->IsSelected(dstarD0pi,AliRDHFCuts::kTracks); // quality cuts on tracks
		// region of interest + topological cuts + PID
		Int_t isSelected=fCuts->IsSelected(dstarD0pi,AliRDHFCuts::kCandidate); //selected
		//apply selections
		if(!isTkSelected) continue;
		if(!isSelected) continue;
		if(!fCuts->IsInFiducialAcceptance(dstarD0pi->Pt(),dstarD0pi->YDstar())) continue;    
		Int_t mcLabelDStar = -999;
		if(fmontecarlo){
			// find associated MC particle for D* ->D0toKpi
			mcLabelDStar = dstarD0pi->MatchToMC(413,421,pdgDgDStartoD0pi,pdgDgD0toKpi,fmcArray,kFALSE);
			if(mcLabelDStar>=0) isDStarMCtag = kTRUE;
		}
		
		ptDStar = dstarD0pi->Pt();
		phiDStar = dstarD0pi->Phi(); 
		etaDStar = dstarD0pi->Eta();
		
		phiDStar = fCorrelator->SetCorrectPhiRange(phiDStar); // set the phi of the D meson in the correct range
		
		Int_t ptbin=fCuts->PtBin(dstarD0pi->Pt());
		
		Double_t dmDStarWindow =0.0019;// 0.0019 = 3 sigma
		Double_t mD0Window=0.074;
		
		if (!fSystem){ // pp
			if (ptbin==1)  mD0Window = 0.026; //0.5-1
			if (ptbin==2)  mD0Window = 0.022; //1-2
			if (ptbin==3)  mD0Window = 0.024; //2-3
			if (ptbin==4)  mD0Window = 0.032;
			if (ptbin==5)  mD0Window = 0.032;
			if (ptbin==6)  mD0Window = 0.036;
			if (ptbin==7)  mD0Window = 0.036;
			if (ptbin==8)  mD0Window = 0.036;
			if (ptbin==9)  mD0Window = 0.058;
			if (ptbin==10) mD0Window = 0.058;
			if (ptbin>10)  mD0Window = 0.074;
		}
		if(fSystem){// PbPb
			if (ptbin==0)  mD0Window = 0.032; //1-1
			if (ptbin==1)  mD0Window = 0.032; //2-3
			if (ptbin==2)  mD0Window = 0.032; //3-4
			if (ptbin==3)  mD0Window = 0.032; //4-5
			if (ptbin==4)  mD0Window = 0.036; //5-6
			if (ptbin==5)  mD0Window = 0.036; //6-8
			if (ptbin==6)  mD0Window = 0.055; //8-12
			if (ptbin==7)  mD0Window = 0.074; //12-16
			if (ptbin==8)  mD0Window = 0.074; //16-24
			if (ptbin==9)  mD0Window = 0.074; //24-35 
		 }
		
		invMassDZero = dstarD0pi->InvMassD0();
		((TH2F*)fOutput->FindObject("D0InvMass"))->Fill(ptDStar,invMassDZero);
		
		deltainvMDStar = dstarD0pi->DeltaInvMass();
		
	
		//good candidates
		if (TMath::Abs(invMassDZero-mPDGD0)<mD0Window){
			
			((TH2F*)fOutput->FindObject("DeltaInvMass"))->Fill(ptDStar,deltainvMDStar);
			if(TMath::Abs(deltainvMDStar-(mPDGDstar-mPDGD0))<dmDStarWindow){ // is in DStar peak region?
				
				((TH1F*)fOutput->FindObject("RecoPtDStar"))->Fill(ptDStar);
				isInPeak = kTRUE;
				((TH2F*)fOutput->FindObject("PhiEtaTrigger"))->Fill(phiDStar,etaDStar);
			}
		}// end if good candidates
		
		//sidebands
		if (TMath::Abs(invMassDZero-mPDGD0)>1.3*mD0Window && TMath::Abs(invMassDZero-mPDGD0)<4.*mD0Window ){
			((TH2F*)fOutput->FindObject("bkgDeltaInvMass"))->Fill(ptDStar,deltainvMDStar);
			((TH2F*)fOutput->FindObject("D0InvMassinSB"))->Fill(ptDStar,invMassDZero);
			
			if(TMath::Abs(deltainvMDStar-(mPDGDstar-mPDGD0))<dmDStarWindow){ // is in DStar peak region?
				((TH1F*)fOutput->FindObject("RecoPtBkg"))->Fill(ptDStar);
				isInSideBand = kTRUE;
				((TH2F*)fOutput->FindObject("PhiEtaSideBand"))->Fill(phiDStar,etaDStar);
			}
			
		}//end if sidebands
		// getting the number of triggers in the MCtag D* case 
		

        if(fmontecarlo && isDStarMCtag) ((TH1F*)fOutput->FindObject("MCtagPtDStar"))->Fill(ptDStar);
		
		
		if(!isInPeak && !isInSideBand) continue; // skip if it is not side band or peak event - SAVE CPU TIME
		
		isDStarCand = kTRUE;
		
	    fCorrelator->SetTriggerParticleProperties(ptDStar,phiDStar,etaDStar); // pass to the object the necessary trigger part parameters
		
		Short_t daughtercharge = ((AliAODTrack*)theD0particle->GetDaughter(0))->Charge();
		fCorrelator->SetTriggerParticleDaughterCharge(daughtercharge);
		
		
		Int_t trackiddaugh0 = ((AliAODTrack*)theD0particle->GetDaughter(0))->GetID();
		Int_t trackiddaugh1 = ((AliAODTrack*)theD0particle->GetDaughter(1))->GetID();
		Int_t trackidsoftPi = ((AliAODTrack*)dstarD0pi->GetBachelor())->GetID();
		
		Bool_t execPool = fCorrelator->ProcessEventPool();
		if(fmixing && !execPool) {
			AliInfo("Mixed event analysis: pool is not ready");
			continue;
		}
		
		Int_t NofEventsinPool = 1;
		if(fmixing) NofEventsinPool = fCorrelator->GetNofEventsInPool();
				
		for (Int_t jMix =0; jMix < NofEventsinPool; jMix++){// loop on events in the pool; if it is SE analysis, stops at one
		
			Bool_t analyzetracks = fCorrelator->ProcessAssociatedTracks(jMix);
			
			if(!analyzetracks) {
				AliInfo("AliHFCorrelator::Cannot process the track array");
				continue;
			}
		
            //initialization of variables for correlations with leading particles
            Double_t DeltaPhiLeading = -999.;
			Double_t DeltaEtaLeading = -999.;
			//Double_t ptleading = -999.;
            Int_t labelleading = -999;
		
			Int_t NofTracks = fCorrelator->GetNofTracks();
			
			for(Int_t iTrack = 0; iTrack<NofTracks; iTrack++){ // looping on track candidates
			
				Bool_t runcorrelation = fCorrelator->Correlate(iTrack);
				if(!runcorrelation) continue;
			
				Double_t DeltaPhi = fCorrelator->GetDeltaPhi();
				Double_t DeltaEta = fCorrelator->GetDeltaEta();
			
				AliReducedParticle * hadron = fCorrelator->GetAssociatedParticle();
				
				Double_t ptHad = hadron->Pt();
				Double_t phiHad = hadron->Phi();
				Double_t etaHad = hadron->Eta(); 
				Int_t label = hadron->GetLabel(); 
				Int_t trackid = hadron->GetID();
				
				phiHad = fCorrelator->SetCorrectPhiRange(phiHad);
				
				if(!fmixing){ // skip D* Daughetrs
					if(trackid == trackiddaugh0) continue;
					if(trackid == trackiddaugh1) continue;
					if(trackid == trackidsoftPi) continue;
				}
			
				// from here on it is up to the user to decide what object to fill
				
				if(fmontecarlo && isDStarMCtag){ // check correlations of MC tagged DStars in MonteCarlo
				
					Bool_t* PartSource = fCuts2->IsMCpartFromHF(label,fmcArray); // check source of associated particle (hadron/kaon/K0)
					FillMCTagCorrelations(ptDStar,DeltaPhi,DeltaEta,ptHad,PartSource);
					
					
					((TH3F*)fOutputMC->FindObject("MCPhiEtaPart"))->Fill(phiHad,etaHad,0);
					if(PartSource[0]) ((TH3F*)fOutputMC->FindObject("MCPhiEtaPart"))->Fill(phiHad,etaHad,1);
					if(PartSource[1]) ((TH3F*)fOutputMC->FindObject("MCPhiEtaPart"))->Fill(phiHad,etaHad,2);
					if(PartSource[2]&&PartSource[0]) ((TH3F*)fOutputMC->FindObject("MCPhiEtaPart"))->Fill(phiHad,etaHad,3);
					if(PartSource[2]&&PartSource[1]) ((TH3F*)fOutputMC->FindObject("MCPhiEtaPart"))->Fill(phiHad,etaHad,4);
					if(PartSource[3]) ((TH3F*)fOutputMC->FindObject("MCPhiEtaPart"))->Fill(phiHad,etaHad,5);
				
				
				}
			
				if(isInPeak)  {
				
					if(fselect==1) ((TH3D*)fOutput->FindObject("DPhiDStarHadron"))->Fill(DeltaPhi,ptDStar,DeltaEta);
					if(fselect==2) ((TH3D*)fOutput->FindObject("DPhiDStarKaon"))->Fill(DeltaPhi,ptDStar,DeltaEta);
					if(fselect==3) ((TH3D*)fOutput->FindObject("DPhiDStarKZero"))->Fill(DeltaPhi,ptDStar,DeltaEta);
				    ((TH2F*)fOutput->FindObject("PhiEtaPart"))->Fill(phiHad,etaHad);
					//counterPeak++; // count tracks per peak per event
				
				}
			
				if(isInSideBand) {
				
					if(fselect==1) ((TH3D*)fOutput->FindObject("bkgDPhiDStarHadron"))->Fill(DeltaPhi,ptDStar,DeltaEta);
					if(fselect==2) ((TH3D*)fOutput->FindObject("bkgDPhiDStarKaon"))->Fill(DeltaPhi,ptDStar,DeltaEta);
					if(fselect==3) ((TH3D*)fOutput->FindObject("bkgDPhiDStarKZero"))->Fill(DeltaPhi,ptDStar,DeltaEta);
				
					
				
					//counterSB++;
				}
			
			
			} // end loop on track candidates
            
            // fill the leading particle histograms
            
            if(isInPeak) ((TH3D*)fOutput->FindObject("LeadingCand"))->Fill(DeltaPhiLeading,ptDStar,DeltaEtaLeading);
			if(isInSideBand) ((TH3D*)fOutput->FindObject("LeadingSB"))->Fill(DeltaPhiLeading,ptDStar,DeltaEtaLeading);
            
			if(fmontecarlo && isDStarMCtag){
                Bool_t* LeadPartSource = fCuts2->IsMCpartFromHF(labelleading,fmcArray);
                FillMCTagLeadingCorrelations(ptDStar,DeltaPhiLeading,DeltaEtaLeading,LeadPartSource);
                
            }
            
		} // end loop on events in the pool
				
	}// end loop on D* candidates		
	
	
	Bool_t updated = fCorrelator->PoolUpdate();
	
	if(updated) EventMixingChecks(aodEvent);
	if(!updated) AliInfo("Pool was not updated");
	
	
	
		
} //end the exec

//________________________________________ terminate ___________________________
void AliAnalysisTaskDStarCorrelations::Terminate(Option_t*)
{    
	// The Terminate() function is the last function to be called during
	// a query. It always runs on the client, it can be used to present
	// the results graphically or save the results to file.
	
	AliAnalysisTaskSE::Terminate();
	
	fOutput = dynamic_cast<TList*> (GetOutputData(1));
	if (!fOutput) {     
		printf("ERROR: fOutput not available\n");
		return;
	}

	return;
}


//_____________________________________________________
void AliAnalysisTaskDStarCorrelations::DefineHistoForAnalysis(){
	
	Double_t Pi = TMath::Pi();
	Int_t nbinscorr = 32;
	Double_t lowcorrbin = -0.5*Pi - Pi/32; // shift the bin by half the width so that at 0 is it the bin center
	Double_t upcorrbin = 1.5*Pi - Pi/32;	
	
	// ========================= histograms for both Data and MonteCarlo
	
	
	TH1D * NofEvents = new TH1D("NofEvents","NofEvents",11,0,11);
	fOutput->Add(NofEvents);
	
	
	
	
	
	TH2F *D0InvMass = new TH2F("D0InvMass","K#pi invariant mass distribution",300,0,30,1500,0.5,3.5);
	fOutput->Add(D0InvMass);
	
	TH2F *D0InvMassinSB = new TH2F("D0InvMassinSB","K#pi invariant mass distribution in sb",300,0,30,1500,0.5,3.5);
	fOutput->Add(D0InvMassinSB);
	
	TH2F *DeltaInvMass = new TH2F("DeltaInvMass","K#pi#pi - K#pi invariant mass distribution",300,0,30,750,0.1,0.2);
	fOutput->Add(DeltaInvMass);
	
	TH2F *bkgDeltaInvMass = new TH2F("bkgDeltaInvMass","K#pi#pi - K#pi invariant mass distribution",300,0,30,750,0.1,0.2);
	fOutput->Add(bkgDeltaInvMass);
	
	TH1F *RecoPtDStar = new TH1F("RecoPtDStar","RECO DStar pt distribution",50,0,50);
	fOutput->Add(RecoPtDStar);
	
	TH1F *RecoPtBkg = new TH1F("RecoPtBkg","RECO pt distribution side bands",50,0,50);
	fOutput->Add(RecoPtBkg);
	
	TH1F *MCtagPtDStar = new TH1F("MCtagPtDStar","RECO pt of MCtagged DStars side bands",50,0,50);
	fOutput->Add(MCtagPtDStar);
	
	TH2F *KZeroSpectra = new TH2F("KZeroSpectra","Spectra of K0s",500,0.3,0.8,250,0,25);
	if(fselect==3) fOutput->Add(KZeroSpectra);
	
	TH2F *KZeroSpectraifHF = new TH2F("KZeroSpectraifHF","Spectra of K0s in association with a D*",500,0.3,0.8,250,0,25);
	if(fselect==3) fOutput->Add(KZeroSpectraifHF);
	
	TH1D * NofTracksInPeak = new TH1D("NofTracksInPeak","NofTracksInPeak",500,0.5,500.5);
	fOutput->Add(NofTracksInPeak);
	
	TH1D * NofTracksInSB = new TH1D("NofTracksInSB","NofTracksInSB",500,0.5,500.5);
	fOutput->Add(NofTracksInSB);
	
	TH2I * EventMixingCheck = new TH2I("EventMixingCheck","EventMixingCheck",5,-0.5,4.5,7,-0.5,6.5);
	if(fmixing) fOutput->Add(EventMixingCheck);
	
	

	
	
	TH2F * PhiEtaTrigger = new TH2F("PhiEtaTrigger","#phi distribution of the trigger particle",36,-0.5*Pi,1.5*Pi,18,-0.9,0.9);
	fOutput->Add(PhiEtaTrigger);
	
	TH2F * PhiEtaSideBand = new TH2F("PhiEtaSideBand","#phi distribution of the sideband particle",36,-0.5*Pi,1.5*Pi,18,-0.9,0.9);
	fOutput->Add(PhiEtaSideBand);
	
	TH2F * PhiEtaPart = new TH2F("PhiEtaPart","#phi distribution of the associated particle",36,-0.5*Pi,1.5*Pi,18,-0.9,0.9);
	fOutput->Add(PhiEtaPart);
	
	
	//correlations histograms
	TString histoname1 = "DPhiDStar";
	if(fselect==1) histoname1 += "Hadron";
	if(fselect==2) histoname1 += "Kaon";
	if(fselect==3) histoname1 += "KZero";
	
	
	TH3D * DPhiDStar = new TH3D(histoname1.Data(),histoname1.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
	
	TH3D * DPhiDStarKZero1 = new TH3D("DPhiDStarKZero1","DPhiDStarKZero1",nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
	
	//side band background histograms
	TString histoname2 = "bkg";
	histoname2 += histoname1;
	TH3D * bkgDPhiDStar = new TH3D(histoname2.Data(),histoname2.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
	TH3D * bkgDPhiDStarKZero1 = new TH3D("bkgDPhiDStarKZero1","bkgDPhiDStarKZero1",nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
	
	
	fOutput->Add(DPhiDStar);
	
	if(fselect==3){fOutput->Add(DPhiDStarKZero1);}
	
	fOutput->Add(bkgDPhiDStar);
	
	if(fselect==3){fOutput->Add(bkgDPhiDStarKZero1);}
	
	
	// leading particle
	TH3D * leadingcand = new TH3D("LeadingCand","LeadingCand",nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
	TH3D * leadingsidebands = new TH3D("LeadingSB","LeadingSB",nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
	
	fOutput->Add(leadingcand);
	fOutput->Add(leadingsidebands);
	
	// ========================= histos for analysis on MC only
	
	TH1D * EventTypeMC = new TH1D("EventTypeMC","EventTypeMC",100,-0.5,99.5);
	if(fmontecarlo) fOutputMC->Add(EventTypeMC);
	
	TH1F * MCSources = new TH1F("MCSources","Origin of associated particles in MC", 10, -0.5, 9.5);
	MCSources->GetXaxis()->SetBinLabel(1,"All ");
	MCSources->GetXaxis()->SetBinLabel(2," from hadron Heavy flavour");
	MCSources->GetXaxis()->SetBinLabel(3," from c->D");
	MCSources->GetXaxis()->SetBinLabel(4," from b->D");
	MCSources->GetXaxis()->SetBinLabel(5," from b->B");
	MCSources->GetXaxis()->SetBinLabel(6," from quark Heavy flavour");
	MCSources->GetXaxis()->SetBinLabel(7," from c");
	MCSources->GetXaxis()->SetBinLabel(8," from b");
	
	if(fmontecarlo) fOutputMC->Add(MCSources);
    
    // leading particle from mc source
    TH1F * LeadingMCSources = new TH1F("LeadingMCSources","Origin of associated leading particles in MC", 10, -0.5, 9.5);
	LeadingMCSources->GetXaxis()->SetBinLabel(1,"All ");
	LeadingMCSources->GetXaxis()->SetBinLabel(2," from hadron Heavy flavour");
	LeadingMCSources->GetXaxis()->SetBinLabel(3," from c->D");
	LeadingMCSources->GetXaxis()->SetBinLabel(4," from b->D");
	LeadingMCSources->GetXaxis()->SetBinLabel(5," from b->B");
	LeadingMCSources->GetXaxis()->SetBinLabel(6," from quark Heavy flavour");
	LeadingMCSources->GetXaxis()->SetBinLabel(7," from c");
	LeadingMCSources->GetXaxis()->SetBinLabel(8," from b");
	
	if(fmontecarlo) fOutputMC->Add(LeadingMCSources);
	
    // all hadrons
	TString histoname3 = "MCTag";
	histoname3 += histoname1;
	TH3D * MCTagDPhiDStar = new TH3D(histoname3.Data(),histoname3.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
	
	TString histoname44 = "CharmDOrigin";
	histoname44 += histoname1;
	histoname44 += "MC";
	
	TH3D * CharmDOriginDPhiDStar = new TH3D(histoname44.Data(),histoname44.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
	
	
	TString histoname54 = "BeautyDOrigin";
	histoname54 += histoname1;
	histoname54 += "MC";
	TH3D * BeautyDOriginDPhiDStar = new TH3D(histoname54.Data(),histoname54.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
	
	TString histoname55 = "BeautyBOrigin";
	histoname55 += histoname1;
	histoname55 += "MC";
	TH3D * BeautyBOriginDPhiDStar = new TH3D(histoname55.Data(),histoname55.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
	
	TString histoname4 = "CharmQuarkOrigin";
	histoname4 += histoname1;
	histoname4 += "MC";
	TH3D * CharmQuarkOriginDPhiDStar = new TH3D(histoname4.Data(),histoname4.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
	
	TString histoname5 = "BeautyQuarkOrigin";
	histoname5 += histoname1;
	histoname5 += "MC";
	TH3D * BeautyQuarkOriginDPhiDStar = new TH3D(histoname5.Data(),histoname5.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
    
    if(fmontecarlo){
        
        fOutputMC->Add(MCTagDPhiDStar);
        fOutputMC->Add(CharmDOriginDPhiDStar);
        fOutputMC->Add(BeautyDOriginDPhiDStar);
        fOutputMC->Add(BeautyBOriginDPhiDStar);
        fOutputMC->Add(CharmQuarkOriginDPhiDStar);
        fOutputMC->Add(BeautyQuarkOriginDPhiDStar);
        
	}
    
    // ========================= histos for analysis on MC
    // all leading hadron
	TString Leadinghistoname3 = "LeadingMCTag";
	Leadinghistoname3 += histoname1;
	TH3D * LeadingMCTagDPhiDStar = new TH3D(Leadinghistoname3.Data(),Leadinghistoname3.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
    
	TString Leadinghistoname44 = "LeadingCharmDOrigin";
	Leadinghistoname44 += histoname1;
	Leadinghistoname44 += "MC";
	
	TH3D * LeadingCharmDOriginDPhiDStar = new TH3D(Leadinghistoname44.Data(),Leadinghistoname44.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
	
	
	TString Leadinghistoname54 = "LeadingBeautyDOrigin";
	Leadinghistoname54 += histoname1;
	Leadinghistoname54 += "MC";
	TH3D * LeadingBeautyDOriginDPhiDStar = new TH3D(Leadinghistoname54.Data(),Leadinghistoname54.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
	
	TString Leadinghistoname55 = "LeadingBeautyBOrigin";
	Leadinghistoname55 += histoname1;
	Leadinghistoname55 += "MC";
	TH3D * LeadingBeautyBOriginDPhiDStar = new TH3D(Leadinghistoname55.Data(),Leadinghistoname55.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
	
	TString Leadinghistoname4 = "LeadingCharmQuarkOrigin";
	Leadinghistoname4 += histoname1;
	Leadinghistoname4 += "MC";
	TH3D * LeadingCharmQuarkOriginDPhiDStar = new TH3D(Leadinghistoname4.Data(),Leadinghistoname4.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
	
	TString Leadinghistoname5 = "LeadingBeautyQuarkOrigin";
	Leadinghistoname5 += histoname1;
	Leadinghistoname5 += "MC";
	TH3D * LeadingBeautyQuarkOriginDPhiDStar = new TH3D(Leadinghistoname5.Data(),Leadinghistoname5.Data(),nbinscorr,lowcorrbin,upcorrbin,50,0,50,39,-1.95,1.95);
    
    
	
	
	if(fmontecarlo){
		
		fOutputMC->Add(LeadingMCTagDPhiDStar);
		fOutputMC->Add(LeadingCharmDOriginDPhiDStar);
		fOutputMC->Add(LeadingBeautyDOriginDPhiDStar);
		fOutputMC->Add(LeadingBeautyBOriginDPhiDStar);
		fOutputMC->Add(LeadingCharmQuarkOriginDPhiDStar);
		fOutputMC->Add(LeadingBeautyQuarkOriginDPhiDStar);
		
	}
	
	TH3F * MCPhiEtaPart = new TH3F("MCPhiEtaPart","#phi distribution of the associated particle",36,-0.5*Pi,1.5*Pi,50,-2.5,2.5,6,-0.5,6.5);
	MCPhiEtaPart->GetZaxis()->SetBinLabel(1,"All particles");
	MCPhiEtaPart->GetZaxis()->SetBinLabel(2,"from c quark");
	MCPhiEtaPart->GetZaxis()->SetBinLabel(3,"from b quark");
	MCPhiEtaPart->GetZaxis()->SetBinLabel(4,"from D from c");
	MCPhiEtaPart->GetZaxis()->SetBinLabel(5,"from D from b");
	MCPhiEtaPart->GetZaxis()->SetBinLabel(6,"from B from b");
	if(fmontecarlo) fOutputMC->Add(MCPhiEtaPart);
	
	// ============================= EVENT MIXING CHECKS ======================================
	
	Int_t MaxNofEvents = fCuts2->GetMaxNEventsInPool();
	Int_t MinNofTracks = fCuts2->GetMinNTracksInPool();
	Int_t NofCentBins = fCuts2->GetNCentPoolBins();
	Double_t * CentBins = fCuts2->GetCentPoolBins();
	Int_t NofZVrtxBins = fCuts2->GetNZvtxPoolBins();
	Double_t *ZVrtxBins = fCuts2->GetZvtxPoolBins();
	
	Int_t k =0;
	
	if(fSystem) k = 100; // PbPb centrality
	if(!fSystem) k = NofCentBins; // pp multiplicity
	
	
	Double_t minvalue = CentBins[0];
	Double_t maxvalue = CentBins[NofCentBins+1];
	Double_t Zminvalue = ZVrtxBins[0];
	Double_t Zmaxvalue = ZVrtxBins[NofCentBins+1];
	
	

	Double_t Nevents[]={0,2*MaxNofEvents/10,4*MaxNofEvents/10,6*MaxNofEvents/10,8*MaxNofEvents/10,MaxNofEvents};
	Double_t * events = Nevents;
	
	TH3D * EventsPerPoolBin = new TH3D("EventsPerPoolBin","Number of events in bin pool",NofCentBins,CentBins,NofZVrtxBins,ZVrtxBins,5,events);
	EventsPerPoolBin->GetXaxis()->SetTitle("Centrality/multiplicity ");
	EventsPerPoolBin->GetYaxis()->SetTitle("Z vertex [cm]");
	EventsPerPoolBin->GetZaxis()->SetTitle("Number of events in pool bin");
	if(fmixing) fOutput->Add(EventsPerPoolBin);
	
	Int_t MaxNofTracks = (MaxNofEvents+1)*MinNofTracks;
	Int_t Diff = MaxNofTracks-MinNofTracks;
	
	Double_t Ntracks[]={MinNofTracks,MinNofTracks+Diff/5,MinNofTracks+2*Diff/5,MinNofTracks+3*Diff/5,MinNofTracks+4*Diff/5,MaxNofTracks};
	Double_t  * trackN = Ntracks;
	
	TH3D * NofTracksPerPoolBin = new TH3D("NofTracksPerPoolBin","Number of tracks in bin pool",NofCentBins,CentBins,NofZVrtxBins,ZVrtxBins,5,trackN);
	NofTracksPerPoolBin->GetXaxis()->SetTitle("Centrality/multiplicity ");
	NofTracksPerPoolBin->GetYaxis()->SetTitle("Z vertex [cm]");
	NofTracksPerPoolBin->GetZaxis()->SetTitle("Number of tracks per bin");
	
	if(fmixing) fOutput->Add(NofTracksPerPoolBin);
	
	TH2D * NofPoolBinCalls = new TH2D("NofPoolBinCalls","Number of tracks in bin pool",NofCentBins,CentBins,NofZVrtxBins,ZVrtxBins);
	NofPoolBinCalls->GetXaxis()->SetTitle("Centrality/multiplicity ");
	NofPoolBinCalls->GetYaxis()->SetTitle("Z vertex [cm]");
	if(fmixing) fOutput->Add(NofPoolBinCalls);
	

	
	TH2D * EventProps = new TH2D("EventProps","Number of tracks in bin pool",k,minvalue,maxvalue,100,Zminvalue,Zmaxvalue);
	EventProps->GetXaxis()->SetTitle("Centrality/multiplicity ");
	EventProps->GetYaxis()->SetTitle("Z vertex [cm]");
	if(fmixing) fOutput->Add(EventProps);
	
}



//____________________________  Function for MC correlations ___________________________________________________
void AliAnalysisTaskDStarCorrelations::FillMCTagCorrelations(Double_t ptTrig, Double_t DelPhi,  Double_t DelEta, Double_t ptTrack, Bool_t *mcSource){

	
	
	
			
		if(fselect==1) ((TH3D*)fOutputMC->FindObject("MCTagDPhiDStarHadron"))->Fill(DelPhi,ptTrig,DelEta);
		if(fselect==2 && ptTrack <1.5) ((TH3D*)fOutputMC->FindObject("MCTagDPhiDStarKaon"))->Fill(DelPhi,ptTrig,DelEta);
		if(fselect==3) ((TH3D*)fOutputMC->FindObject("MCTagDPhiDStarKZero"))->Fill(DelPhi,ptTrig,DelEta);
		
	
		
		((TH1F*)fOutputMC->FindObject("MCSources"))->Fill(0);
		
		if(fDebugLevel){
			std::cout << "MC source " << mcSource[0] << " "  << mcSource[1] << " " << mcSource[2] << " " << mcSource[3] << std::endl;
		
			if(mcSource[0]) std::cout << "mcSource 0 " << std::endl;
			if(mcSource[1]) std::cout << "mcSource 1 " << std::endl;
			if(mcSource[2]) std::cout << "mcSource 2 " << std::endl;
			if(mcSource[3]) std::cout << "mcSource 3 " << std::endl;
		
		}
		if(mcSource[0]){ // is from charm quark
			((TH1F*)fOutputMC->FindObject("MCSources"))->Fill(5); // all HF quarks
			((TH1F*)fOutputMC->FindObject("MCSources"))->Fill(6); //  charm quarks
			if(fselect==1) ((TH3D*)fOutputMC->FindObject("CharmQuarkOriginDPhiDStarHadronMC"))->Fill(DelPhi,ptTrig,DelEta);
			if(fselect==2 && ptTrack <1.5) ((TH3D*)fOutputMC->FindObject("CharmQuarkOriginDPhiDStarKaonMC"))->Fill(DelPhi,ptTrig,DelEta);
			if(fselect==3) ((TH3D*)fOutputMC->FindObject("CharmQuarkOriginDPhiDStarKZeroMC"))->Fill(DelPhi,ptTrig,DelEta);
		}
		
		if(mcSource[1]){ // is from b quark
			((TH1F*)fOutputMC->FindObject("MCSources"))->Fill(5); // all HF quarks
			((TH1F*)fOutputMC->FindObject("MCSources"))->Fill(7); // beauty quarks
			if(fselect==1) ((TH3D*)fOutputMC->FindObject("BeautyQuarkOriginDPhiDStarHadronMC"))->Fill(DelPhi,ptTrig,DelEta);
			if(fselect==2 && ptTrack <1.5) ((TH3D*)fOutputMC->FindObject("BeautyQuarkOriginDPhiDStarKaonMC"))->Fill(DelPhi,ptTrig,DelEta);
			if(fselect==3) ((TH3D*)fOutputMC->FindObject("BeautyQuarkOriginDPhiDStarKZeroMC"))->Fill(DelPhi,ptTrig,DelEta);
			
		}
		
		if(mcSource[2]&&mcSource[0]){ // is from D meson and charm quark
			((TH1F*)fOutputMC->FindObject("MCSources"))->Fill(1); // all HF mesons
			((TH1F*)fOutputMC->FindObject("MCSources"))->Fill(2); //  charm + D
			if(fselect==1) ((TH3D*)fOutputMC->FindObject("CharmDOriginDPhiDStarHadronMC"))->Fill(DelPhi,ptTrig,DelEta);
			if(fselect==2 && ptTrack <1.5) ((TH3D*)fOutputMC->FindObject("CharmDOriginDPhiDStarKaonMC"))->Fill(DelPhi,ptTrig,DelEta);
			if(fselect==3) ((TH3D*)fOutputMC->FindObject("CharmDOriginDPhiDStarKZeroMC"))->Fill(DelPhi,ptTrig,DelEta);
		}
		
		if(mcSource[2]&&mcSource[1]){ // is from D meson and b quark
			((TH1F*)fOutputMC->FindObject("MCSources"))->Fill(1); // all HF mesons
			((TH1F*)fOutputMC->FindObject("MCSources"))->Fill(3); //  beauty + D
			if(fselect==1) ((TH3D*)fOutputMC->FindObject("BeautyDOriginDPhiDStarHadronMC"))->Fill(DelPhi,ptTrig,DelEta);
			if(fselect==2 && ptTrack <1.5) ((TH3D*)fOutputMC->FindObject("BeautyDOriginDPhiDStarKaonMC"))->Fill(DelPhi,ptTrig,DelEta);
			if(fselect==3) ((TH3D*)fOutputMC->FindObject("BeautyDOriginDPhiDStarKZeroMC"))->Fill(DelPhi,ptTrig,DelEta);
		}
		
	return;
}

//____________________________  Function for MC leading part correlations ___________________________________________________
void AliAnalysisTaskDStarCorrelations::FillMCTagLeadingCorrelations(Double_t ptTrig, Double_t DelPhi,  Double_t DelEta, Bool_t *mcSource){
    // correlations with leading hadron on MC
    
	if(fselect==1) ((TH3D*)fOutputMC->FindObject("LeadingMCTagDPhiDStarHadron"))->Fill(DelPhi,ptTrig,DelEta);
	
    
    
	((TH1F*)fOutputMC->FindObject("LeadingMCSources"))->Fill(0);
	
	if(fDebugLevel){ std::cout << "MC source " << mcSource[0] << " "  << mcSource[1] << " " << mcSource[2] << " " << mcSource[3] << std::endl;
    
		if(mcSource[0]) std::cout << "mcSource 0 " << std::endl;
		if(mcSource[1]) std::cout << "mcSource 1 " << std::endl;
		if(mcSource[2]) std::cout << "mcSource 2 " << std::endl;
		if(mcSource[3]) std::cout << "mcSource 3 " << std::endl;
    }
    
	if(mcSource[0]){ // is from charm quark
		((TH1F*)fOutputMC->FindObject("LeadingMCSources"))->Fill(5); // all HF quarks
		((TH1F*)fOutputMC->FindObject("LeadingMCSources"))->Fill(6); //  charm quarks
		if(fselect==1) ((TH3D*)fOutputMC->FindObject("LeadingCharmQuarkOriginDPhiDStarHadronMC"))->Fill(DelPhi,ptTrig,DelEta);
		
	}
    
	if(mcSource[1]){ // is from b quaLeadingrk
		((TH1F*)fOutputMC->FindObject("LeadingMCSources"))->Fill(5); // all HF quarks
		((TH1F*)fOutputMC->FindObject("LeadingMCSources"))->Fill(7); // beauty quarks
		if(fselect==1) ((TH3D*)fOutputMC->FindObject("LeadingBeautyQuarkOriginDPhiDStarHadronMC"))->Fill(DelPhi,ptTrig,DelEta);
		
		
	}
    
	if(mcSource[2]&&mcSource[0]){ // is from D meson and charm quark
		((TH1F*)fOutputMC->FindObject("LeadingMCSources"))->Fill(1); // all HF mesons
		((TH1F*)fOutputMC->FindObject("LeadingMCSources"))->Fill(2); //  charm + D
		if(fselect==1) ((TH3D*)fOutputMC->FindObject("LeadingCharmDOriginDPhiDStarHadronMC"))->Fill(DelPhi,ptTrig,DelEta);
		
	}
    
	if(mcSource[2]&&mcSource[1]){ // is from D meson and b quark
		((TH1F*)fOutputMC->FindObject("LeadingMCSources"))->Fill(1); // all HF mesons
		((TH1F*)fOutputMC->FindObject("LeadingMCSources"))->Fill(3); //  beauty + D
		if(fselect==1) ((TH3D*)fOutputMC->FindObject("LeadingBeautyDOriginDPhiDStarHadronMC"))->Fill(DelPhi,ptTrig,DelEta);
		
	}
	
    
	return;
}


//____________________________  Run checks on event mixing ___________________________________________________
void AliAnalysisTaskDStarCorrelations::EventMixingChecks(AliAODEvent* AOD){
	
	AliCentrality *centralityObj = 0;
	Int_t multiplicity = -1;
	Double_t MultipOrCent = -1;
	
	// get the pool for event mixing
	if(!fSystem){ // pp
		multiplicity = AOD->GetNTracks();
		MultipOrCent = multiplicity; // convert from Int_t to Double_t
	}
	if(fSystem){ // PbPb
		
		centralityObj = AOD->GetHeader()->GetCentralityP();
		MultipOrCent = centralityObj->GetCentralityPercentileUnchecked("V0M");
		AliInfo(Form("Centrality is %f", MultipOrCent));
	}
	
	AliAODVertex *vtx = AOD->GetPrimaryVertex();
	Double_t zvertex = vtx->GetZ(); // zvertex
	
	
	
	
	AliEventPool * pool = fCorrelator->GetPool();
	

	
	
	((TH2D*)fOutput->FindObject("NofPoolBinCalls"))->Fill(MultipOrCent,zvertex); // number of calls of pool
	((TH2D*)fOutput->FindObject("EventProps"))->Fill(MultipOrCent,zvertex); // event properties
	
	((TH3D*)fOutput->FindObject("EventsPerPoolBin"))->Fill(MultipOrCent,zvertex,pool->NTracksInPool()); // number of events in the pool
	((TH3D*)fOutput->FindObject("NofTracksPerPoolBin"))->Fill(MultipOrCent,zvertex,pool->GetCurrentNEvents()); // number of calls of pool
}
	




