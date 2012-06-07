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

#include <TDatabasePDG.h>
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
#include "AliEventPoolManager.h"
#include "AliVParticle.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliNormalizationCounter.h"
#include "AliReducedParticle.h"



ClassImp(AliAnalysisTaskDStarCorrelations)


//__________________________________________________________________________
AliAnalysisTaskDStarCorrelations::AliAnalysisTaskDStarCorrelations() :
AliAnalysisTaskSE(),
fhandler(0x0),
fPoolMgr(0x0),      
fmcArray(0x0),
fCounter(0x0),
fselect(0),
fmontecarlo(kFALSE),
fmixing(kFALSE),
fEvents(0),
fDebug(0),

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
fPoolMgr(0x0),   
fmcArray(0x0),
fCounter(0x0),
fselect(0),
fmontecarlo(kFALSE),
fmixing(kFALSE),
fEvents(0),
fDebug(0),

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
	if(fPoolMgr) {delete fPoolMgr; fPoolMgr = 0;}    
	if(fmcArray) {delete fmcArray; fmcArray = 0;}
	if(fCounter) {delete fCounter; fCounter = 0;}
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
	
   
	// definition of the Pool Manager for Event Mixing
	
	Int_t MaxNofEvents = 200;
	Int_t MinNofTracks = 1000;
	
	Int_t NofMultiplicityBins = 5;
	Double_t MBins[]={0,20,40,60,80,500};
	Double_t * MultiplicityBins = MBins;
	
	Int_t NofZVrtxBins = 5;
	Double_t ZBins[]={-10,-5,-2.5,2.5,5,10};
	Double_t *ZVrtxBins = ZBins;
	
	
	fPoolMgr = new AliEventPoolManager(MaxNofEvents, MinNofTracks, NofMultiplicityBins, MultiplicityBins, NofZVrtxBins, ZVrtxBins);
	
	
}
//_________________________________________________
void AliAnalysisTaskDStarCorrelations::UserExec(Option_t *){

	cout << " " << endl;
	cout << "=================================================================================" << endl;
	
	if(fselect==1) cout << "TASK::Correlation with hadrons "<< endl;
	if(fselect==2) cout << "TASK::Correlation with kaons "<< endl;
	if(fselect==3) cout << "TASK::Correlation with kzeros "<< endl;
	
	if (!fInputEvent) {
		Error("UserExec","NO EVENT FOUND!");
		return;
	}
	
	AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
	if(!aodEvent){
	  AliError("AOD event not found!");
	  return;
	}
	Double_t pi = TMath::Pi();
	
	fEvents++; // event counter
	((TH1D*)fOutput->FindObject("NofEvents"))->Fill(0);
	fmcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	if(fmontecarlo && !fmcArray){
	  AliError("Array of MC particles not found");
	  return;
	}
	
	// initialize the pool for event mixing
	Int_t multiplicity = aodEvent->GetNTracks();
	AliAODVertex *vtx = aodEvent->GetPrimaryVertex();
	Double_t zvertex = vtx->GetZ();
	Double_t multip = multiplicity;
	
	if(TMath::Abs(zvertex)>=10 || multip>500 || multip == 0) {
		AliInfo(Form("Event with Zvertex = %.2f cm and multiplicity = %.0f out of pool bounds, SKIPPING",zvertex,multip));
		return;
	}
	
	
	AliEventPool* pool = fPoolMgr->GetEventPool(multip, zvertex);
	if (!pool) AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f cm", multip, zvertex));
	
	

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

	
	Double_t mPDGD0=TDatabasePDG::Instance()->GetParticle(421)->Mass();
	Double_t mPDGDstar=TDatabasePDG::Instance()->GetParticle(413)->Mass();
	
	
	
	if(fselect ==3){// check the K0 invariant mass
	TObjArray * fillkzerospectra = AcceptAndReduceKZero(aodEvent, 0,0);
	delete fillkzerospectra;
	}
	
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

		if(fmontecarlo){
			// find associated MC particle for D* ->D0toKpi
			Int_t mcLabelDStar = dstarD0pi->MatchToMC(413,421,pdgDgDStartoD0pi,pdgDgD0toKpi,fmcArray);
			if(mcLabelDStar>=0) isDStarMCtag = kTRUE;
		}
		
		ptDStar = dstarD0pi->Pt();
		phiDStar = dstarD0pi->Phi(); 
		etaDStar = dstarD0pi->Eta();
		
		Int_t ptbin=fCuts->PtBin(dstarD0pi->Pt());
		
		Double_t dmDStarWindow =0.0019;// 0.0019 = 3 sigma
		Double_t mD0Window=0.074;
		
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

		
		invMassDZero = dstarD0pi->InvMassD0();
		((TH2F*)fOutput->FindObject("D0InvMass"))->Fill(ptDStar,invMassDZero);
		
		deltainvMDStar = dstarD0pi->DeltaInvMass();
		
		
		
		//good candidates
		if (TMath::Abs(invMassDZero-mPDGD0)<mD0Window){
			
			((TH2F*)fOutput->FindObject("DeltaInvMass"))->Fill(ptDStar,deltainvMDStar);
			if(TMath::Abs(deltainvMDStar-(mPDGDstar-mPDGD0))<dmDStarWindow){ // is in DStar peak region?
				
				((TH1F*)fOutput->FindObject("RecoPtDStar"))->Fill(ptDStar);
				isInPeak = kTRUE;
			}
		}// end if good candidates
		
		//sidebands
		if (TMath::Abs(invMassDZero-mPDGD0)>1.3*mD0Window && TMath::Abs(invMassDZero-mPDGD0)<4.*mD0Window ){
			((TH2F*)fOutput->FindObject("bkgDeltaInvMass"))->Fill(ptDStar,deltainvMDStar);
			((TH2F*)fOutput->FindObject("D0InvMassinSB"))->Fill(ptDStar,invMassDZero);
			
			if(TMath::Abs(deltainvMDStar-(mPDGDstar-mPDGD0))<dmDStarWindow){ // is in DStar peak region?
				((TH1F*)fOutput->FindObject("RecoPtBkg"))->Fill(ptDStar);
				isInSideBand = kTRUE;
			}
			
		}//end if sidebands
        

		
		if(!isInPeak && !isInSideBand) continue; // skip if it is not side band or peak event - SAVE CPU TIME
		
	
		Int_t trackiddaugh0 = ((AliAODTrack*)theD0particle->GetDaughter(0))->GetID();
		Int_t trackiddaugh1 = ((AliAODTrack*)theD0particle->GetDaughter(1))->GetID();
		Int_t trackidsoftPi = ((AliAODTrack*)dstarD0pi->GetBachelor())->GetID();
		
		ptDStar = dstarD0pi->Pt();
		phiDStar = dstarD0pi->Phi(); 
		etaDStar = dstarD0pi->Eta();
				
		if(!fmixing){ // analysis on Single Event
						
			TObjArray* selectedtracks = new TObjArray();
			if(fselect==1 || fselect ==2)	selectedtracks = AcceptAndReduceTracks(aodEvent);
			if(fselect==3) {cout << " 2 "<< endl; selectedtracks = AcceptAndReduceKZero(aodEvent,iDStartoD0pi,1);}	
			Int_t noftracks = selectedtracks->GetEntriesFast(); 
			Int_t counterPeak =0;
			Int_t counterSB = 0;

			for(Int_t i =0; i<noftracks; i++){ // loop on tracks/k0 candidates in aod event

				AliReducedParticle *redpart = (AliReducedParticle*)selectedtracks->At(i);
				Double_t phiHad=redpart->Phi();

				Double_t ptHad=redpart->Pt();
				Double_t etaHad=redpart->Eta();

				Int_t label = redpart->GetLabel();

				Int_t trackid = redpart->GetID();


				// check that the track is not a D* daughter
				if(trackid == trackiddaugh0) continue;
				if(trackid == trackiddaugh1) continue;
				if(trackid == trackidsoftPi) continue;
				
				if(fmontecarlo && isDStarMCtag){ // check correlations of MC tagged DStars in MonteCarlo
					
					Int_t PartSource = fCuts2->IsMCpartFromHF(label,fmcArray); // check source of associated particle (hadron/kaon/K0)

					FillMCTagCorrelations(ptDStar,phiDStar,etaDStar,ptHad,phiHad,etaHad,PartSource);
				}
				if(isInPeak)  {
					FillCorrelations(ptDStar,phiDStar,etaDStar,phiHad,etaHad);// function for correlations
					counterPeak++;
					if (phiDStar > 1.5*pi) phiDStar = phiDStar - 2*pi;
					if (phiDStar < -0.5*pi) phiDStar = phiDStar + 2*pi;

					((TH1F*)fOutput->FindObject("PhiTrigger"))->Fill(phiDStar);

					
						
						if (phiHad > 1.5*pi) phiHad = phiHad - 2*pi;
						if (phiHad < -0.5*pi) phiHad = phiHad + 2*pi;
						((TH1F*)fOutput->FindObject("PhiPart"))->Fill(phiHad);
					
				}
				cout << "fill6" << endl;
				if(isInSideBand) {

					FillSideBandCorrelations(ptDStar,phiDStar,etaDStar,phiHad,etaHad); // function for sidebands
					if (phiDStar > 1.5*pi) phiDStar = phiDStar - 2*pi;
					if (phiDStar < -0.5*pi) phiDStar = phiDStar + 2*pi;
					((TH1F*)fOutput->FindObject("PhiSideBand"))->Fill(phiDStar);

					counterSB++;
				}
				cout << "fill7" << endl;
			} // end loop on tracks	
			
			if(counterPeak) ((TH1D*)fOutput->FindObject("NofTracksInPeak"))->Fill(counterPeak);
			if(counterSB) ((TH1D*)fOutput->FindObject("NofTracksInSB"))->Fill(counterSB);
			
			
			
			
		} // end if SE Analysis
		
		if(fmixing) { // analysis on Mixed Events
			if (pool->IsReady()|| pool->GetCurrentNEvents()>=5){ // check if the pool is ready
		
				pool->PrintInfo();
	
				Int_t multbinflag = pool->MultBinIndex();
				Int_t zvtxflag = pool->ZvtxBinIndex();
 
				if(isInPeak){ cout << "check 1" << endl;
					((TH2I*)fOutput->FindObject("EventMixingCheck"))->Fill(multbinflag,zvtxflag); 
					cout << "filling" << endl;}

				TObjArray* mixedtracks = 0x0;
				
				for (Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++) {//loop over the events in the pool
					mixedtracks = pool->GetEvent(jMix);
					if(!mixedtracks) cout << "No Mixed tracks " << endl;
					Int_t jMax = mixedtracks->GetEntriesFast();

					for(Int_t iMix =0; iMix<jMax; iMix++){ //loop on the tracks of the event
						AliVParticle *redpart = (AliVParticle*)mixedtracks->At(iMix);
						Double_t phiHad=redpart->Phi();
						Double_t etaHad=redpart->Eta();
						Double_t ptHad=redpart->Pt();
						Int_t label = redpart->GetLabel();
						
						
						if(fmontecarlo && isDStarMCtag){ // check correlations of MC tagged DStars in MonteCarlo
							
							Int_t PartSource = fCuts2->IsMCpartFromHF(label,fmcArray); // check source of associated particle (hadron/kaon/K0)
							
							FillMCTagCorrelations(ptDStar,phiDStar,etaDStar,ptHad,phiHad,etaHad,PartSource);
						}

						if(isInPeak) {
							FillCorrelations(ptDStar,phiDStar,etaDStar,phiHad,etaHad);// function for correlations

								if (phiDStar > 1.5*pi) phiDStar = phiDStar - 2*pi;
								if (phiDStar < -0.5*pi) phiDStar = phiDStar + 2*pi;
								
								((TH1F*)fOutput->FindObject("PhiTrigger"))->Fill(phiDStar);
								if (phiHad > 1.5*pi) phiHad = phiHad - 2*pi;
								if (phiHad < -0.5*pi) phiHad = phiHad + 2*pi;
								((TH1F*)fOutput->FindObject("PhiPart"))->Fill(phiHad);
							

							
						}
						
						if(isInSideBand) FillSideBandCorrelations(ptDStar,phiDStar,etaDStar,phiHad,etaHad); // function for sidebands

						} // end loop on tracks
				}// end loop on events in pool
			} // end if pool is ready
			
		} // end ME analysis
		
	}// end loop on D* candidates		
	
	if(fmixing) { // update the pool for Event Mixing
		if(fselect==1 || fselect==2)pool->UpdatePool(AcceptAndReduceTracks(aodEvent)); // updating the pool for hadrons
		if(fselect==3) pool->UpdatePool(AcceptAndReduceKZero(aodEvent,0,0)); // updating the pool for K0s
	}
	//cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> END OF THE EVENT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
	
	PostData(1,fOutput); // set the outputs
	PostData(3,fCounter); // set the outputs
		
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
TObjArray*  AliAnalysisTaskDStarCorrelations::AcceptAndReduceTracks(AliAODEvent* inputEvent){
	
	Int_t nTracks = inputEvent->GetNTracks();
	AliAODVertex * vtx = inputEvent->GetPrimaryVertex();
	Double_t Bz = inputEvent->GetMagneticField();

	
	TObjArray* tracksClone = new TObjArray;
	tracksClone->SetOwner(kTRUE);
	for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) {
		AliAODTrack* track = inputEvent->GetTrack(iTrack);
		if (! track) continue;

		if(!fCuts2->IsHadronSelected(track,vtx,Bz)) continue; // apply selection of hadrons
			
		
		if(fselect ==2){	
			if(!fCuts2->CheckKaonCompatibility(track,fmontecarlo,fmcArray)) continue; // check if it is a Kaon - data and MC
			}
      
		AliVParticle * part = (AliVParticle*)track;
		tracksClone->Add(new AliReducedParticle(part->Eta(), part->Phi(), part->Pt(),track->GetLabel(),track->GetID()));
		
	}
	return tracksClone;
}

//_____________________________________________________
TObjArray*  AliAnalysisTaskDStarCorrelations::AcceptAndReduceKZero(AliAODEvent* inputEvent, Int_t loopindex, Int_t plotassociation){
	
	Int_t nOfVZeros = inputEvent->GetNumberOfV0s();
	TObjArray* KZeroClone = new TObjArray;
	AliAODVertex *vertex1 = (AliAODVertex*)inputEvent->GetPrimaryVertex();
	Int_t v0label = -1;
	Int_t pdgDgK0toPipi[2] = {211,211};
	Double_t mPDGK0=TDatabasePDG::Instance()->GetParticle(310)->Mass();
	const Int_t size = inputEvent->GetNumberOfV0s();
	Int_t prevnegID[size];
	Int_t prevposID[size];
	for(Int_t iv0 =0; iv0< nOfVZeros; iv0++){// loop on all v0 candidates
		AliAODv0 *v0 = (static_cast<AliAODEvent*> (inputEvent))->GetV0(iv0);
		if(!v0) {cout << "is not a v0 at step " << iv0 << endl; continue;}
		
		if(!fCuts2->IsKZeroSelected(v0,vertex1)) continue; // check if it is a k0
	    
		// checks to avoid double counting
		Int_t negID = -999;
		Int_t posID = -998;
		//Int_t a = 0;
		prevnegID[iv0] = -997;
		prevposID[iv0]= -996;
		negID = v0->GetNegID();
		posID = v0->GetPosID();
		Bool_t DoubleCounting = kFALSE;
		
		for(Int_t k=0; k<iv0; k++){
			if(((negID==prevnegID[k])&&(posID==prevposID[k]))||((negID==prevposID[k])&&(posID==prevnegID[k]))){
				DoubleCounting = kTRUE;
				//a=k;
				break;
			}//end if
		} // end for
		
		if(DoubleCounting) {cout << "SKIPPING DOUBLE COUNTING" << endl;continue;}
		else{ prevposID[iv0] = posID; prevnegID[iv0] = negID;}
		
		if(fmontecarlo)	v0label = v0->MatchToMC(310,fmcArray, 0, pdgDgK0toPipi); //match a K0 short
		Double_t k0pt = v0->Pt();
		Double_t k0eta = v0->Eta();
		Double_t k0Phi = v0->Phi();
		Double_t k0InvMass = v0->MassK0Short();	

		if (loopindex == 0) {
			if(!plotassociation) ((TH2F*)fOutput->FindObject("KZeroSpectra"))->Fill(k0InvMass,k0pt); // spectra for all k0
			if(plotassociation) ((TH2F*)fOutput->FindObject("KZeroSpectraifHF"))->Fill(k0InvMass,k0pt); // spectra for k0 in association with a D*
		}
		// if there are more D* candidates per event, loopindex == 0 makes sure you fill the mass spectra only once!
		
		if(TMath::Abs(k0InvMass-mPDGK0)>3*0.004) continue; // select candidates within 3 sigma
		KZeroClone->Add(new AliReducedParticle(k0eta , k0Phi, k0pt,v0label,0));
		
	}
	
	return KZeroClone;
}

//_____________________________________________________
void AliAnalysisTaskDStarCorrelations::DefineHistoForAnalysis(){
	
	Double_t Pi = TMath::Pi();
	Int_t nbinscorr = 32;
	Double_t lowcorrbin = -0.5*Pi - Pi/32; // shift the bin by half the width so that at 0 is it the bin center
	Double_t upcorrbin = 1.5*Pi - Pi/32;	
	
	// ========================= histograms for both Data and MonteCarlo
	
	
	TH1D * NofEvents = new TH1D("NofEvents","NofEvents",1,0,1);
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
	
	TH1F * PhiTrigger = new TH1F("PhiTrigger","#phi distribution of the trigger particle",36,-0.5*Pi,1.5*Pi);
	fOutput->Add(PhiTrigger);
	
	TH1F * PhiSideBand = new TH1F("PhiSideBand","#phi distribution of the sideband particle",36,-0.5*Pi,1.5*Pi);
	fOutput->Add(PhiSideBand);
	
	TH1F * PhiPart = new TH1F("PhiPart","#phi distribution of the associated particle",36,-0.5*Pi,1.5*Pi);
	fOutput->Add(PhiPart);


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


//____________________________  Function for correlations ___________________________________________________
void AliAnalysisTaskDStarCorrelations::FillCorrelations(Double_t ptTrig, Double_t phiTrig, Double_t etaTrig, Double_t phiTrack, Double_t etaTrack){
	Double_t pi = TMath::Pi();
	Double_t deltaPhi, deltaEta;
	deltaPhi = phiTrig - phiTrack;
	deltaEta = etaTrig - etaTrack;
	// set correct Delta Phi range
	if (deltaPhi > 1.5*pi -pi/32) deltaPhi = deltaPhi - 2*pi;
	if (deltaPhi < -0.5*pi -pi/32) deltaPhi = deltaPhi + 2*pi;
	cout << "CRASH CHECK 1 " << endl;
		if(fselect==1) ((TH3D*)fOutput->FindObject("DPhiDStarHadron"))->Fill(deltaPhi,ptTrig,deltaEta);
		if(fselect==2) ((TH3D*)fOutput->FindObject("DPhiDStarKaon"))->Fill(deltaPhi,ptTrig,deltaEta);
		if(fselect==3) ((TH3D*)fOutput->FindObject("DPhiDStarKZero"))->Fill(deltaPhi,ptTrig,deltaEta);
	cout << "CRASH CHECK 2 " << endl;
	
	return;
}

//____________________________  Function for sidebands ___________________________________________________
void AliAnalysisTaskDStarCorrelations::FillSideBandCorrelations(Double_t ptTrig, Double_t phiTrig, Double_t etaTrig, Double_t phiTrack, Double_t etaTrack){
	
	Double_t pi = TMath::Pi();
	Double_t deltaPhi, deltaEta;
	deltaPhi = phiTrig - phiTrack;
	deltaEta = etaTrig - etaTrack;
	// set correct Delta Phi range
	if (deltaPhi > 1.5*pi -pi/32) deltaPhi = deltaPhi - 2*pi;
	if (deltaPhi < -0.5*pi -pi/32) deltaPhi = deltaPhi + 2*pi;
	cout << "CRASH CHECK bkg 1 " << endl;
		if(fselect==1) ((TH3D*)fOutput->FindObject("bkgDPhiDStarHadron"))->Fill(deltaPhi,ptTrig,deltaEta);
		if(fselect==2) ((TH3D*)fOutput->FindObject("bkgDPhiDStarKaon"))->Fill(deltaPhi,ptTrig,deltaEta);
		if(fselect==3) ((TH3D*)fOutput->FindObject("bkgDPhiDStarKZero"))->Fill(deltaPhi,ptTrig,deltaEta);
	cout << "CRASH CHECK bkg 2 " << endl;
	
	return;
	
	
}

//____________________________  Function for sidebands ___________________________________________________
void AliAnalysisTaskDStarCorrelations::FillMCTagCorrelations(Double_t ptTrig, Double_t phiTrig,  Double_t etaTrig, Double_t ptTrack, Double_t phiTrack, Double_t etaTrack, Int_t mcSource){
Double_t deltaPhi = phiTrig - phiTrack;
Double_t deltaEta = etaTrig - etaTrack;
// set correct Delta Phi range
	
	Double_t pi = TMath::Pi();
	
if (deltaPhi > 1.5*pi -pi/32) deltaPhi = deltaPhi - 2*pi;
if (deltaPhi < -0.5*pi -pi/32) deltaPhi = deltaPhi + 2*pi;
cout << "fill1" << endl;
	if(fselect==1) ((TH3D*)fOutput->FindObject("MCTagDPhiDStarHadron"))->Fill(deltaPhi,ptTrig,deltaEta);
	if(fselect==2 && ptTrack <1.5) ((TH3D*)fOutput->FindObject("MCTagDPhiDStarKaon"))->Fill(deltaPhi,ptTrig,deltaEta);
	if(fselect==3) ((TH3D*)fOutput->FindObject("MCTagDPhiDStarKZero"))->Fill(deltaPhi,ptTrig,deltaEta);



((TH1F*)fOutput->FindObject("MCSources"))->Fill(0);
cout << "fill2" << endl;
if (mcSource==44){ // is from charm ->D
	if(fselect==1) ((TH3D*)fOutput->FindObject("CharmDOriginDPhiDStarHadronMC"))->Fill(deltaPhi,ptTrig,deltaEta);
	if(fselect==2 && ptTrack <1.5) ((TH3D*)fOutput->FindObject("CharmDOriginDPhiDStarKaonMC"))->Fill(deltaPhi,ptTrig,deltaEta);
	if(fselect==3) ((TH3D*)fOutput->FindObject("CharmDOriginDPhiDStarKZeroMC"))->Fill(deltaPhi,ptTrig,deltaEta);
	
	
	((TH1F*)fOutput->FindObject("MCSources"))->Fill(1);
	((TH1F*)fOutput->FindObject("MCSources"))->Fill(2);
	}
cout << "fill3" << endl;
if (mcSource==54){ // is from beauty -> D
	if(fselect==1) ((TH3D*)fOutput->FindObject("BeautyDOriginDPhiDStarHadronMC"))->Fill(deltaPhi,ptTrig,deltaEta);
	if(fselect==2 && ptTrack <1.5) ((TH3D*)fOutput->FindObject("BeautyDOriginDPhiDStarKaonMC"))->Fill(deltaPhi,ptTrig,deltaEta);
	if(fselect==3) ((TH3D*)fOutput->FindObject("BeautyDOriginDPhiDStarKZeroMC"))->Fill(deltaPhi,ptTrig,deltaEta);
	if(fselect==3) ((TH1F*)fOutput->FindObject("MCSources"))->Fill(1);
	if(fselect==3) ((TH1F*)fOutput->FindObject("MCSources"))->Fill(3);
	}
cout << "fill4" << endl;
if (mcSource==55){ // is from beauty ->B
	if(fselect==1) ((TH3D*)fOutput->FindObject("BeautyBOriginDPhiDStarHadronMC"))->Fill(deltaPhi,ptTrig,deltaEta);
	if(fselect==2 && ptTrack <1.5) ((TH3D*)fOutput->FindObject("BeautyBOriginDPhiDStarKaonMC"))->Fill(deltaPhi,ptTrig,deltaEta);
	if(fselect==3) ((TH3D*)fOutput->FindObject("BeautyOriginBDPhiDStarKZeroMC"))->Fill(deltaPhi,ptTrig,deltaEta);
	if(fselect==3) ((TH1F*)fOutput->FindObject("MCSources"))->Fill(1);
	if(fselect==3) ((TH1F*)fOutput->FindObject("MCSources"))->Fill(4);
	}
	return;
}





