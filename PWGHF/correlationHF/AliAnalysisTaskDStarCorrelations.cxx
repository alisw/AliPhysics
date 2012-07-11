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




ClassImp(AliAnalysisTaskDStarCorrelations)


//__________________________________________________________________________
AliAnalysisTaskDStarCorrelations::AliAnalysisTaskDStarCorrelations() :
AliAnalysisTaskSE(),
fhandler(0x0),
//fPoolMgr(0x0),      
fmcArray(0x0),
fCounter(0x0),
fCorrelator(0x0),
fselect(0),
fmontecarlo(kFALSE),
fmixing(kFALSE),
fSystem(kFALSE),
fEvents(0),
fDebug(0),
fDisplacement(0),

fOutput(0x0),            
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
fEvents(0),
fDebug(0),
fDisplacement(0),

fOutput(0x0),                
fCuts(0),
fCuts2(AsscCuts)
{
	fCuts=cuts;
	Info("AliAnalysisTaskDStarCorrelations","Calling Constructor");
	DefineInput(0, TChain::Class());
	DefineOutput(1,TList::Class()); // histos from data
	DefineOutput(2,AliRDHFCutsDStartoKpipi::Class()); // my cuts
	DefineOutput(3,AliNormalizationCounter::Class());   // normalization
	DefineOutput(4,AliHFAssociatedTrackCuts::Class()); // my cuts
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
	if(fCuts) {delete fCuts; fCuts = 0;}
	if(fCuts2) {delete fCuts2; fCuts2=0;}

}

//___________________________________________________________
void AliAnalysisTaskDStarCorrelations::Init(){
	//
	// Initialization
	//
	if(fDebug > 1) printf("AliAnalysisTaskDStarCorrelations::Init() \n");
	
	AliRDHFCutsDStartoKpipi* copyfCuts=new AliRDHFCutsDStartoKpipi(*fCuts);
	
	
	
	
	// Post the D* cuts
	PostData(2,copyfCuts);
	
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
	
	// define histograms
	DefineHistoForAnalysis();
	fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(3)->GetContainer()->GetName()));
	fCounter->Init();
	
    Double_t Pi = TMath::Pi();
	fCorrelator = new AliHFCorrelator("Correlator",fCuts2,fSystem); // fCuts2 is the hadron cut object, fSystem to switch between pp or PbPb
	fCorrelator->SetDeltaPhiInterval((-0.5-1./32)*Pi,(1.5-1./32)*Pi); // set correct phi interval
	fCorrelator->SetEventMixing(fmixing); //set kFALSE/kTRUE for mixing Off/On
	fCorrelator->SetAssociatedParticleType(fselect); // set 1/2/3 for hadron/kaons/kzeros
	fCorrelator->SetApplyDisplacementCut(fDisplacement); //set kFALSE/kTRUE for using the displacement cut
	fCorrelator->SetUseMC(fmontecarlo);
	Bool_t pooldef = fCorrelator->DefineEventPool();
	
	if(!pooldef) AliInfo("Warning:: Event pool not defined properly");

	
	PostData(1,fOutput); // set the outputs
	PostData(3,fCounter); // set the outputs
}
//_________________________________________________
void AliAnalysisTaskDStarCorrelations::UserExec(Option_t *){

	
	
	cout << " " << endl;
	cout << "=================================================================================" << endl;
	if(!fmixing){
	if(fselect==1) cout << "TASK::Correlation with hadrons on SE "<< endl;
	if(fselect==2) cout << "TASK::Correlation with kaons on SE "<< endl;
	if(fselect==3) cout << "TASK::Correlation with kzeros on SE "<< endl;
	}
	if(fmixing){
		if(fselect==1) cout << "TASK::Correlation with hadrons on ME "<< endl;
		if(fselect==2) cout << "TASK::Correlation with kaons on ME "<< endl;
		if(fselect==3) cout << "TASK::Correlation with kzeros on ME "<< endl;
	}
	if (!fInputEvent) {
		Error("UserExec","NO EVENT FOUND!");
		return;
	}
	
	AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
	if(!aodEvent){
	  AliError("AOD event not found!");
	  return;
	}
	
	fCorrelator->SetAODEvent(aodEvent); // set the event to be processed
	
	fEvents++; // event counter
	((TH1D*)fOutput->FindObject("NofEvents"))->Fill(0);
	fmcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	
	if(fmontecarlo && !fmcArray){
	  AliError("Array of MC particles not found");
	  return;
	}
	Bool_t isEvSel=fCuts->IsEventSelected(aodEvent);
	if(!isEvSel) return;
	((TH1D*)fOutput->FindObject("NofEvents"))->Fill(1);
	//
	Bool_t correlatorON = fCorrelator->Initialize(); //define the pool for mixing
	if(!correlatorON) {
		AliInfo("AliHFCorrelator didn't initialize the pool correctly or processed a bad event");
		return;
	}
	((TH1D*)fOutput->FindObject("NofEvents"))->Fill(2);
	
	if(fmontecarlo) fCorrelator->SetMCArray(fmcArray);
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
				Double_t label = hadron->GetLabel(); 
				Double_t trackid = hadron->GetID();
				
				phiHad = fCorrelator->SetCorrectPhiRange(phiHad);
				
				if(!fmixing){ // skip D* Daughetrs
					if(trackid == trackiddaugh0) continue;
					if(trackid == trackiddaugh1) continue;
					if(trackid == trackidsoftPi) continue;
				}
			
				// from here on it is up to the user to decide what object to fill
				
				if(fmontecarlo && isDStarMCtag){ // check correlations of MC tagged DStars in MonteCarlo
				
					Int_t PartSource = fCuts2->IsMCpartFromHF(label,fmcArray); // check source of associated particle (hadron/kaon/K0)
				
					FillMCTagCorrelations(ptDStar,DeltaPhi,DeltaEta,ptHad,PartSource);
				
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
		} // end loop on events in the pool
				
	}// end loop on D* candidates		
	Bool_t updated = fCorrelator->PoolUpdate();
	if(!updated) AliInfo("Pool was not updated");
	
		//cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> END OF THE EVENT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
	
	
		
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
	
	
	TH1D * NofEvents = new TH1D("NofEvents","NofEvents",3,0,3);
	fOutput->Add(NofEvents);
		
	TH2F *D0InvMass = new TH2F("D0InvMass","K#pi invariant mass distribution",300,0,30,1500,0.5,3.5);
	fOutput->Add(D0InvMass);
	
	TH2F *D0InvMassinSB = new TH2F("D0InvMassinSB","K#pi invariant mass distribution in sb",300,0,30,1500,0.5,3.5);
	fOutput->Add(D0InvMassinSB);
	
	TH2F *DeltaInvMass = new TH2F("DeltaInvMass","K#pi#pi - K#pi invariant mass distribution",300,0,30,750,0.1,0.2);
	fOutput->Add(DeltaInvMass);
	
	TH2F *bkgDeltaInvMass = new TH2F("bkgDeltaInvMass","K#pi#pi - K#pi invariant mass distribution",300,0,30,750,0.1,0.2);
	fOutput->Add(bkgDeltaInvMass);
	
	TH1F *RecoPtDStar = new TH1F("RecoPtDStar","RECO DStar pt distribution",30,0,30);
	fOutput->Add(RecoPtDStar);
	
	TH1F *RecoPtBkg = new TH1F("RecoPtBkg","RECO pt distribution side bands",30,0,30);
	fOutput->Add(RecoPtBkg);
	
	TH1F *MCtagPtDStar = new TH1F("MCtagPtDStar","RECO pt of MCtagged DStars side bands",30,0,30);
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
	

	
	TH1F * MCSources = new TH1F("MCSources","Origin of associated particles in MC", 10, -0.5, 9.5);
	MCSources->GetXaxis()->SetBinLabel(1,"All ");
	MCSources->GetXaxis()->SetBinLabel(2," from Heavy flavour");
	MCSources->GetXaxis()->SetBinLabel(3," from c->D");
	MCSources->GetXaxis()->SetBinLabel(4," from b->D");
	MCSources->GetXaxis()->SetBinLabel(5," from b->B");
	if(fmontecarlo) fOutput->Add(MCSources);
	
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

	
	TH3D * DPhiDStar = new TH3D(histoname1.Data(),histoname1.Data(),nbinscorr,lowcorrbin,upcorrbin,30,0,30,19,-0.95,0.95);
	
	TH3D * DPhiDStarKZero1 = new TH3D("DPhiDStarKZero1","DPhiDStarKZero1",nbinscorr,lowcorrbin,upcorrbin,30,0,30,19,-0.95,0.95);
	
	//side band background histograms
	TString histoname2 = "bkg";
	histoname2 += histoname1;
	TH3D * bkgDPhiDStar = new TH3D(histoname2.Data(),histoname2.Data(),nbinscorr,lowcorrbin,upcorrbin,30,0,30,19,-0.95,0.95);
	TH3D * bkgDPhiDStarKZero1 = new TH3D("bkgDPhiDStarKZero1","bkgDPhiDStarKZero1",nbinscorr,lowcorrbin,upcorrbin,30,0,30,19,-0.95,0.95);
	
	
	fOutput->Add(DPhiDStar);

	if(fselect==3){fOutput->Add(DPhiDStarKZero1);}
	
	fOutput->Add(bkgDPhiDStar);

	if(fselect==3){fOutput->Add(bkgDPhiDStarKZero1);}


	// ========================= histos for analysis on MC
	TString histoname3 = "MCTag";
	histoname3 += histoname1;
	TH3D * MCTagDPhiDStar = new TH3D(histoname3.Data(),histoname3.Data(),nbinscorr,lowcorrbin,upcorrbin,30,0,30,19,-0.95,0.95);
	 
	TString histoname44 = "CharmDOrigin";
	histoname44 += histoname1;
	histoname44 += "MC";
	
	TH3D * CharmDOriginDPhiDStar = new TH3D(histoname44.Data(),histoname44.Data(),nbinscorr,lowcorrbin,upcorrbin,30,0,30,19,-0.95,0.95);
	
	
	TString histoname54 = "BeautyDOrigin";
	histoname54 += histoname1;
	histoname54 += "MC";
	TH3D * BeautyDOriginDPhiDStar = new TH3D(histoname54.Data(),histoname54.Data(),nbinscorr,lowcorrbin,upcorrbin,30,0,30,19,-0.95,0.95);
	
	TString histoname55 = "BeautyBOrigin";
	histoname55 += histoname1;
	histoname55 += "MC";
	TH3D * BeautyBOriginDPhiDStar = new TH3D(histoname55.Data(),histoname55.Data(),nbinscorr,lowcorrbin,upcorrbin,30,0,30,19,-0.95,0.95);
	
	if(fmontecarlo){
	
	fOutput->Add(MCTagDPhiDStar);
	fOutput->Add(CharmDOriginDPhiDStar);
	fOutput->Add(BeautyDOriginDPhiDStar);
	fOutput->Add(BeautyBOriginDPhiDStar);
	
	}
	

	
}



//____________________________  Function for MC correlations ___________________________________________________
void AliAnalysisTaskDStarCorrelations::FillMCTagCorrelations(Double_t ptTrig, Double_t DelPhi,  Double_t DelEta, Double_t ptTrack, Int_t mcSource){


	if(fselect==1) ((TH3D*)fOutput->FindObject("MCTagDPhiDStarHadron"))->Fill(DelPhi,ptTrig,DelEta);
	if(fselect==2 && ptTrack <1.5) ((TH3D*)fOutput->FindObject("MCTagDPhiDStarKaon"))->Fill(DelPhi,ptTrig,DelEta);
	if(fselect==3) ((TH3D*)fOutput->FindObject("MCTagDPhiDStarKZero"))->Fill(DelPhi,ptTrig,DelEta);



((TH1F*)fOutput->FindObject("MCSources"))->Fill(0);

if (mcSource==44){ // is from charm ->D
	if(fselect==1) ((TH3D*)fOutput->FindObject("CharmDOriginDPhiDStarHadronMC"))->Fill(DelPhi,ptTrig,DelEta);
	if(fselect==2 && ptTrack <1.5) ((TH3D*)fOutput->FindObject("CharmDOriginDPhiDStarKaonMC"))->Fill(DelPhi,ptTrig,DelEta);
	if(fselect==3) ((TH3D*)fOutput->FindObject("CharmDOriginDPhiDStarKZeroMC"))->Fill(DelPhi,ptTrig,DelEta);
	
	
	((TH1F*)fOutput->FindObject("MCSources"))->Fill(1);
	((TH1F*)fOutput->FindObject("MCSources"))->Fill(2);
	}

if (mcSource==54){ // is from beauty -> D
	if(fselect==1) ((TH3D*)fOutput->FindObject("BeautyDOriginDPhiDStarHadronMC"))->Fill(DelPhi,ptTrig,DelEta);
	if(fselect==2 && ptTrack <1.5) ((TH3D*)fOutput->FindObject("BeautyDOriginDPhiDStarKaonMC"))->Fill(DelPhi,ptTrig,DelEta);
	if(fselect==3) ((TH3D*)fOutput->FindObject("BeautyDOriginDPhiDStarKZeroMC"))->Fill(DelPhi,ptTrig,DelEta);
	if(fselect==3) ((TH1F*)fOutput->FindObject("MCSources"))->Fill(1);
	if(fselect==3) ((TH1F*)fOutput->FindObject("MCSources"))->Fill(3);
	}

if (mcSource==55){ // is from beauty ->B
	if(fselect==1) ((TH3D*)fOutput->FindObject("BeautyBOriginDPhiDStarHadronMC"))->Fill(DelPhi,ptTrig,DelEta);
	if(fselect==2 && ptTrack <1.5) ((TH3D*)fOutput->FindObject("BeautyBOriginDPhiDStarKaonMC"))->Fill(DelPhi,ptTrig,DelEta);
	if(fselect==3) ((TH3D*)fOutput->FindObject("BeautyOriginBDPhiDStarKZeroMC"))->Fill(DelPhi,ptTrig,DelEta);
	if(fselect==3) ((TH1F*)fOutput->FindObject("MCSources"))->Fill(1);
	if(fselect==3) ((TH1F*)fOutput->FindObject("MCSources"))->Fill(4);
	}
	return;
}





