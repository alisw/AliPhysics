/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------
// Example of task (running locally, on AliEn and CAF),
// which provides standard way of calculating acceptance and efficiency
// between different steps of the procedure.
// The ouptut of the task is a AliCFContainer from which the efficiencies
// can be calculated
//-----------------------------------------------------------------------
// Author : R. Vernet, Consorzio Cometa - Catania (it)
//-----------------------------------------------------------------------
// Modification done by X. Lopez - LPC Clermont (fr)
//-----------------------------------------------------------------------
// Modification done by S. Ahn - LPC Clermont (fr), Konkuk University (kr)
//-----------------------------------------------------------------------



#ifndef ALICFMUONRESUPSILON_CXX
#define ALICFMUONRESUPSILON_CXX

#include "TH1.h"
#include "TParticle.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TCanvas.h"

#include "AliHeader.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"

#include "AliVEvent.h"
#include "AliPhysicsSelection.h"

#include "AliCFMuonResUpsilon.h"
#include "AliCFManager.h"
#include "AliCFCutBase.h"
#include "AliCFContainer.h"

#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliESDtrack.h"
#include "AliESDInputHandler.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"

ClassImp(AliCFMuonResUpsilon)

//__________________________________________________________________________
AliCFMuonResUpsilon::AliCFMuonResUpsilon() :
	AliAnalysisTaskSE(""),
	fReadAODData(kFALSE),
	fReadMCInfo(kFALSE),
  fCFManager(0x0),
	fnevts(0x0),
	fIsPhysSelMB(kFALSE),
	fIsPhysSelMUON(kFALSE),
  fQAHistList(0x0),
	fPDG(0),
	fPtMin(0),
	fPtMax(0),
	fYMin(0),
	fYMax(0),
	fTrigClassMuon(""),
	fTrigClassInteraction(""),
	fDistinguishTrigClass(kTRUE)
{

  //Default ctor

}
//___________________________________________________________________________
AliCFMuonResUpsilon::AliCFMuonResUpsilon(const Char_t* name) :
  AliAnalysisTaskSE(name),
	fReadAODData(kFALSE),
	fReadMCInfo(kFALSE),
  fCFManager(0x0),
	fnevts(0x0),
	fIsPhysSelMB(kFALSE),
	fIsPhysSelMUON(kFALSE),
  fQAHistList(0x0),
	fPDG(0),
	fPtMin(0),
	fPtMax(0),
	fYMin(0),
	fYMax(0),
	fTrigClassMuon(""),
	fTrigClassInteraction(""),
	fDistinguishTrigClass(kTRUE)
{

  // Constructor. Initialization of Inputs and Outputs

  Info("AliCFMuonResUpsilon","Calling Constructor");

 	fnevts	= new TH1D("fnevts","fnevts",5,0,5);							// nevent CINT1B,CMUS1B + PhysSel

	TString nameside[3]={"AC","B","E"};

  SetTrigClassMuonName();
  SetTrigClassInteracName();
  SetTrigClassSideName(nameside);

	DefineInput(0, TChain::Class());
	DefineOutput(1,TH1D::Class());
  DefineOutput(2,AliCFContainer::Class());
}

//___________________________________________________________________________
AliCFMuonResUpsilon& AliCFMuonResUpsilon::operator=(const AliCFMuonResUpsilon& c) 
{

  // Assignment operator
	
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c);
    fReadAODData = c.fReadAODData;
		fReadMCInfo = c.fReadMCInfo;
    fCFManager  = c.fCFManager;
    fQAHistList = c.fQAHistList;
		fnevts = c.fnevts;
		fPDG = c.fPDG;
		fPtMin = c.fPtMin;
		fPtMax = c.fPtMax;
		fYMin = c.fYMin;
		fYMax = c.fYMax;
  }
  return *this;
}

//___________________________________________________________________________
AliCFMuonResUpsilon::AliCFMuonResUpsilon(const AliCFMuonResUpsilon& c) :
  AliAnalysisTaskSE(c),
  fReadAODData(c.fReadAODData),
	fReadMCInfo(c.fReadMCInfo),
  fCFManager(c.fCFManager),
	fnevts(c.fnevts),
	fIsPhysSelMB(kFALSE),
	fIsPhysSelMUON(kFALSE),
  fQAHistList(c.fQAHistList),
	fPDG(c.fPDG),
	fPtMin(c.fPtMin),
	fPtMax(c.fPtMax),
	fYMin(c.fYMin),
	fYMax(c.fYMax),
  fTrigClassMuon(c.fTrigClassMuon),
  fTrigClassInteraction(c.fTrigClassInteraction),
  fDistinguishTrigClass(c.fDistinguishTrigClass)
{

  // Copy Constructor

}

//___________________________________________________________________________
AliCFMuonResUpsilon::~AliCFMuonResUpsilon() {
  
  //destructor
 
  Info("~AliCFMuonResUpsilon","Calling Destructor");
  if (fCFManager)           delete fCFManager ;
	if (fnevts)								delete fnevts ;
  if (fQAHistList) 					{fQAHistList->Clear(); delete fQAHistList;}
}
//___________________________________________________________________________
void AliCFMuonResUpsilon::UserCreateOutputObjects() {
	// UserCreateOutputObjects

 	//TH1D *fnevts	= new TH1D("fnevts","fnevts",5,0,5);							// nevent CINT1B,CMUS1B + PhysSel

  PostData(1,fnevts) ;
  PostData(2,fCFManager->GetParticleContainer()) ;
}

//_________________________________________________
void AliCFMuonResUpsilon::UserExec(Option_t *)
{
  
  // Main loop function
   
 	// variables
	fIsPhysSelMB=kFALSE;			// physics selection : MB
	fIsPhysSelMUON=kFALSE;		// physics selection : MUON

  Double_t containerInput[14] = {0,} ;

  fnevts->Fill(0.5);
	containerInput[0] = 0.5;
 
	if(!fReadAODData) { 	// ESD-based ANALYSIS

		Bool_t flag=kFALSE;

		if(fReadMCInfo) {
  		if (!fMCEvent) {
    		Error("UserExec","NO MC EVENT FOUND!");
    		return;
  		}

			// MC part ----------------------------------------------------------------------------------
			fCFManager->SetMCEventInfo(fMCEvent);  
  		AliStack *stack = fMCEvent->Stack();
			Int_t npart=fMCEvent->GetNumberOfTracks();

			//printf("ESD: npart=%d\n", npart);

  		for (Int_t ipart=0; ipart<npart; ipart++) { 

				AliMCParticle *mcPart  = (AliMCParticle*)fMCEvent->GetTrack(ipart);

    		// Mother kinematics
  	 		TParticle *part = mcPart->Particle(); 
    		Double_t e = part->Energy();
    		Double_t pz = part->Pz();           
    		Double_t rapmc = Rap(e,pz);

				// Selection of the resonance
    		if(!fCFManager->CheckParticleCuts(AliCFManager::kPartGenCuts,mcPart)) continue;
				if(part->Pt()<0 || part->Pt()>1000) continue;
				if(rapmc<-4 || rapmc>-2.4) continue;

				// Decays kinematics
    		Int_t p0 = part->GetDaughter(0);
    		TParticle *part0 = stack->Particle(p0); 
				AliMCParticle *mcpart0 = new AliMCParticle(part0);
    		Int_t pdg0 = part0->GetPdgCode();

    		Int_t p1 = part->GetDaughter(1);
    		TParticle *part1 = stack->Particle(p1);
				AliMCParticle *mcpart1 = new AliMCParticle(part1);
    		Int_t pdg1 = part1->GetPdgCode();

				Double_t e0 = part0->Energy();
				Double_t pz0 = part0->Pz();
				Double_t py0 = part0->Py();
				Double_t px0 = part0->Px();
				Double_t charge0 = part0->GetPDG()->Charge()/3;
				Double_t pt0 = TMath::Sqrt(px0*px0+py0*py0);
				Double_t theta0 = (180./TMath::Pi())*TMath::ATan2(TMath::Sqrt(px0*px0+py0*py0),pz0);
				Double_t eta0 = part0->Eta();
				Double_t vz0 = part0->Vz();
				Double_t mom0 = part0->P();

    		if(!fCFManager->CheckParticleCuts(AliCFManager::kPartAccCuts,mcpart0)) continue;
				if(pt0<0.5) continue;
				if(theta0<171. || theta0>178.) continue;

				Double_t e1 = part1->Energy();
				Double_t pz1 = part1->Pz();
				Double_t py1 = part1->Py();
				Double_t px1 = part1->Px();
				Double_t charge1 = part1->GetPDG()->Charge()/3;
				Double_t pt1 = TMath::Sqrt(px1*px1+py1*py1);
				Double_t theta1 = (180./TMath::Pi())*TMath::ATan2(TMath::Sqrt(px1*px1+py1*py1),pz1);
				Double_t eta1 = part1->Eta();
				Double_t vz1 = part1->Vz();
				Double_t mom1 = part1->P();

    		if(!fCFManager->CheckParticleCuts(AliCFManager::kPartAccCuts,mcpart1)) continue;
				if(pt1<0.5) continue;
				if(theta1<171. || theta1>178.) continue;

		    if(pdg0==-13 || pdg1==-13) { 

					Double_t ptmc = TMath::Sqrt((px0+px1)*(px0+px1)+(py0+py1)*(py0+py1));
					Double_t imassmc = Imass(e0,px0,py0,pz0,e1,px1,py1,pz1);
					if(!imassmc) continue;

					containerInput[1] = rapmc ;   
					containerInput[2] = ptmc ;
					containerInput[3] = imassmc ;   
					containerInput[4] = 1. ;
					containerInput[5] = pt0 ;   
					containerInput[5] = pt1 ;   
					containerInput[6] = mom0 ;   
					containerInput[6] = mom1 ;   
					containerInput[7] = 1. ;
					containerInput[8] = 0. ;
					containerInput[8] = 0. ;
					containerInput[9] = charge0+charge1 ;
					containerInput[10] = eta0;
					containerInput[10] = eta1;
					containerInput[11] = theta0;
					containerInput[11] = theta1;
					containerInput[12] = vz0;
					containerInput[12] = vz1;
					containerInput[13] = 0;
					containerInput[13] = 0;

					// fill the container at the first step
					fCFManager->GetParticleContainer()->Fill(containerInput,0);		// MC container
					flag=kTRUE;
    		}	// second mu loop
			}	// first mu loop
		} // end fReadMCInfo

		if(fReadMCInfo && !flag) return;
		// RECO : ESD part ----------------------------------------------------------------------------------

		// ESD handler
  	AliESDEvent *fESD = 0x0;

		AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
		if( ! esdH) {
			AliError("Cannot get input event handler");
			return;
		}
  	fESD = esdH->GetEvent();

	
		Int_t ntrk=fESD->GetNumberOfMuonTracks();

		Int_t trigfired=-1;
		Int_t trigside=-1;

		if(!fReadMCInfo) {
			fIsPhysSelMB=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kMB));
			fIsPhysSelMUON=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kMUON));

			if(fDistinguishTrigClass) {
				TString trigclass = fESD->GetFiredTriggerClasses();
				if(trigclass.Contains(fTrigClassMuon)) trigfired = 1;
				else if(trigclass.Contains(fTrigClassInteraction)) trigfired = 0;
				for(Int_t i=0;i<3;i++) {
					if(trigfired==1 && trigclass.Contains(fTrigClassSide[i])) {trigside=i;}
					if(trigfired==0 && trigclass.Contains(fTrigClassSide[i])) {trigside=i;}
				}
			}
	
  		if(trigside==1) {	// Beam-Beam event
				if(trigfired) {	// CMUS1-B
					containerInput[0]=3.5;
					fnevts->Fill(3.5);
					if(fIsPhysSelMUON) { 
						containerInput[0]=4.5;
						fnevts->Fill(4.5);	// PhysSel
					}
				} else {	// CINT1-B
					containerInput[0]=1.5;
					fnevts->Fill(1.5);
					if(fIsPhysSelMB) {
						containerInput[0]=2.5;
						fnevts->Fill(2.5);	// PhysSel
					}
				}
			}
		}
		else trigside = 0;	// for MC

		for(Int_t j=0; j<ntrk-1; j++) {
			for(Int_t jj=j+1; jj<ntrk; jj++) {
				AliESDMuonTrack *mu1 = new AliESDMuonTrack(*(fESD->GetMuonTrack(j)));
				AliESDMuonTrack *mu2 = new AliESDMuonTrack(*(fESD->GetMuonTrack(jj)));

				Double_t trigCondition = 0;
				if (mu1->GetMatchTrigger()>=mu2->GetMatchTrigger()) trigCondition = mu1->GetMatchTrigger()+0.1*mu2->GetMatchTrigger();
				else trigCondition = mu2->GetMatchTrigger()+0.1*mu1->GetMatchTrigger();

				// mu1
				Double_t zr1 = mu1->Charge();
				Double_t pxr1 = mu1->Px();
				Double_t pyr1 = mu1->Py();
				Double_t pzr1 = mu1->Pz();
     		Double_t ptr1 = TMath::Sqrt(pxr1*pxr1+pyr1*pyr1);
	 			Double_t er1 = mu1->E();
				Double_t momr1 = mu1->P();
				Double_t vzr1 = mu1->GetZ();
				Double_t dcar1 = mu1->GetDCA();
				Double_t rabsr1 = mu1->GetRAtAbsorberEnd();
				Double_t thetaabsr1 = 180.*(1.-TMath::ATan(rabsr1/505.)/TMath::Pi());  // [deg]
				Double_t thetaabsr12 = (thetaabsr1*TMath::Pi()/180.)/2.;
				Double_t tanthetaabsr1 = TMath::Tan(thetaabsr12);
				Double_t etar1 = -TMath::Log(tanthetaabsr1);

				if(Rap(er1,pzr1)<-4 || Rap(er1,pzr1)>-2.5) continue;

				// mu2
				Double_t zr2 = mu2->Charge();
				Double_t pxr2 = mu2->Px();
				Double_t pyr2 = mu2->Py();
				Double_t pzr2 = mu2->Pz();
     		Double_t ptr2 = TMath::Sqrt(pxr2*pxr2+pyr2*pyr2);
	 			Double_t er2 = mu2->E();
				Double_t momr2 = mu2->P();
				Double_t vzr2 = mu2->GetZ();
				Double_t dcar2 = mu2->GetDCA();
				Double_t rabsr2 = mu2->GetRAtAbsorberEnd();
				Double_t thetaabsr2 = 180.*(1.-TMath::ATan(rabsr2/505.)/TMath::Pi());  // [deg]
				Double_t thetaabsr22 = (thetaabsr2*TMath::Pi()/180.)/2.;
				Double_t tanthetaabsr2 = TMath::Tan(thetaabsr22);
				Double_t etar2 = -TMath::Log(tanthetaabsr2);

				if(Rap(er2,pzr2)<-4 || Rap(er2,pzr2)>-2.5) continue;

				// for MC
				if(fReadMCInfo) {
					//mu1
					rabsr1 = 0;
					Double_t thetar1 = (180./TMath::Pi())*TMath::ATan2(TMath::Sqrt(pxr1*pxr1+pyr1*pyr1),pzr1);
					thetaabsr1 = thetar1;
					etar1 = mu1->Eta();
					//mu2
					rabsr2 = 0;
					Double_t thetar2 = (180./TMath::Pi())*TMath::ATan2(TMath::Sqrt(pxr2*pxr2+pyr2*pyr2),pzr2);
					thetaabsr2 = thetar2;
					etar2 = mu2->Eta();
				}
					
				if(TMath::Abs(etar1) > 8 || TMath::Abs(etar2) > 8) continue;

				// dimuon
				Double_t ptrec = TMath::Sqrt((pxr1+pxr2)*(pxr1+pxr2)+(pyr1+pyr2)*(pyr1+pyr2));
				Double_t raprec= Rap((er1+er2),(pzr1+pzr2));
				Double_t imassrec = Imass(er1,pxr1,pyr1,pzr1,er2,pxr2,pyr2,pzr2);
				if(!imassrec) continue;

				containerInput[1] = raprec ;   
				containerInput[2] = ptrec ;
				containerInput[3] = imassrec ;   
				containerInput[4] = trigCondition+0.05 ;
				containerInput[5] = ptr1 ;   
				containerInput[5] = ptr2 ;   
				containerInput[6] = momr1;
				containerInput[6] = momr2;
				containerInput[7] = trigside;
				containerInput[8] = rabsr1;
				containerInput[8] = rabsr2;
				containerInput[9] = zr1 + zr2;
				containerInput[10] = etar1;
				containerInput[10] = etar2;
				containerInput[11] = thetaabsr1;
				containerInput[11] = thetaabsr2;
				containerInput[12] = vzr1;
				containerInput[12] = vzr2;
				containerInput[13] = dcar1;
				containerInput[13] = dcar2;
				
	  		if (trigfired==1 && trigside==1) fCFManager->GetParticleContainer()->Fill(containerInput,4);
	  		if (trigfired==0 && trigside==1) fCFManager->GetParticleContainer()->Fill(containerInput,1);
				if(fReadMCInfo) fCFManager->GetParticleContainer()->Fill(containerInput,1); // Rec container
			}	// second mu loop
		}	// first mu loop
	} // end ESD-based ANALYSIS

	else { 	// AOD-based ANALYSIS

		Bool_t flag=kTRUE;

		// AOD handler
		AliAODEvent *fAOD = 0x0;

  	AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
		if( ! aodH) {
			AliError("Cannot get input event handler");
			return;
		}
		fAOD = aodH->GetEvent();

		// MC part -----------------------------------------------------------------------------------

		if(fReadMCInfo) {
			TClonesArray *mcArr = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
			if( ! mcArr) {
			  AliError("Cannot get MC innf in AOD MC branch");
			  return;
			}
			Int_t npart=mcArr->GetEntries();

			for(Int_t i=0; i<npart; i++) {
				AliAODMCParticle *mctrack = (AliAODMCParticle*) mcArr->At(i);
				// select resonance
				if(mctrack->GetPdgCode()!=fPDG) continue;
				// cuts on resonance
				if(!(mctrack->Pt()>0 && mctrack->Pt()<1000) || !(Rap(mctrack->E(),mctrack->Pz())>-4 && Rap(mctrack->E(),mctrack->Pz())<-2.4)) continue;
				Int_t daug0 = mctrack->GetDaughter(0);
				Int_t daug1 = mctrack->GetDaughter(1);
				// daughter1
				AliAODMCParticle *mcdaug0 = (AliAODMCParticle*) mcArr->At(daug0);
				Double_t pt0 = mcdaug0->Pt();
				Double_t mom0 = mcdaug0->P();
				Double_t theta0 = (180./TMath::Pi())*mcdaug0->Theta();
				Double_t eta0 = mcdaug0->Eta();
				Int_t charge0 = (Int_t) mcdaug0->Charge()/3;
				Double_t vz0 = mcdaug0->Zv();
				if(!(pt0>0.5) || !(theta0>171. && theta0<178.)) continue;
				// daughter2
				AliAODMCParticle *mcdaug1 = (AliAODMCParticle*) mcArr->At(daug1);
				Double_t pt1 = mcdaug1->Pt();
				Double_t mom1 = mcdaug1->P();
				Double_t theta1 = (180./TMath::Pi())*mcdaug1->Theta();
				Double_t eta1 = mcdaug1->Eta();
				Int_t charge1 = (Int_t) mcdaug1->Charge()/3;
				Double_t vz1 = mcdaug1->Zv();
				if(!(pt1>0.5) || !(theta1>171. && theta1<178.)) continue;

				if (TMath::Abs(mcdaug0->GetPdgCode())!=13 || TMath::Abs(mcdaug1->GetPdgCode())!=13) continue;

				Double_t rapmc = Rap(mctrack->E(), mctrack->Pz());
				Double_t ptmc = mctrack->Pt();
				Double_t imassmc = mctrack->M();
				if(!imassmc) continue;

				containerInput[1] = rapmc ;   
				containerInput[2] = ptmc ;
				containerInput[3] = imassmc ;   
				containerInput[4] = 1. ;
				containerInput[5] = pt0 ;   
				containerInput[5] = pt1 ;   
				containerInput[6] = mom0 ;   
				containerInput[6] = mom1 ;   
				containerInput[7] = 1. ;
				containerInput[8] = 0. ;
				containerInput[8] = 0. ;
				containerInput[9] = charge0+charge1 ;
				containerInput[10] = eta0;
				containerInput[10] = eta1;
				containerInput[11] = theta0;
				containerInput[11] = theta1;
				containerInput[12] = vz0;
				containerInput[12] = vz1;
				containerInput[13] = 0;
				containerInput[13] = 0;

	  		fCFManager->GetParticleContainer()->Fill(containerInput,0);	// MC container
				flag=kTRUE;
			}
		}

if(fReadMCInfo && !flag) return;

		// RECO : AOD part ----------------------------------------------------------------------------------


		Int_t trigfired=-1;
		Int_t trigside=-1;

		Int_t ntrk = fAOD->GetNumberOfTracks();

		if(!fReadMCInfo) {
			fIsPhysSelMB=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kMB));
			fIsPhysSelMUON=(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kMUON));

			if(fDistinguishTrigClass) {
				TString trigclass = fAOD->GetFiredTriggerClasses();
				if(trigclass.Contains(fTrigClassMuon)) trigfired = 1;
				else if(trigclass.Contains(fTrigClassInteraction)) trigfired = 0;
				for(Int_t i=0;i<3;i++) {
					if(trigfired==1 && trigclass.Contains(fTrigClassSide[i].Data())) {trigside=i;} 
					if(trigfired==0 && trigclass.Contains(fTrigClassSide[i].Data())) {trigside=i;}
				}
			}
	
  		if(trigside==1) {	// Beam-Beam event
				if(trigfired) {	// CMUS1-B
					containerInput[0]=3.5;
					fnevts->Fill(3.5);
					if(fIsPhysSelMUON) { 
						containerInput[0]=4.5;
						fnevts->Fill(4.5);	// PhysSel
					}
				} else {	// CINT1-B
					containerInput[0]=1.5;
					fnevts->Fill(1.5);
					if(fIsPhysSelMB) {
						containerInput[0]=2.5;
						fnevts->Fill(2.5);	// PhysSel
					}
				}
			}
		}
		else trigside = 0;	// for MC


		for(Int_t j=0; j<ntrk; j++) {
			for(Int_t jj=j+1; jj<ntrk; jj++) {
				AliAODTrack *mu1 = fAOD->GetTrack(j);
				if(!mu1->IsMuonTrack() || !(mu1->Y()>-4 && mu1->Y()<-2.5)) continue;
				AliAODTrack *mu2 = fAOD->GetTrack(jj);
				if(!mu2->IsMuonTrack() || !(mu2->Y()>-4 && mu2->Y()<-2.5)) continue;

	     	Double_t trigCondition=0;
	     	if (mu1->GetMatchTrigger()>=mu2->GetMatchTrigger()) trigCondition = mu1->GetMatchTrigger()+0.1*mu2->GetMatchTrigger();
				else trigCondition = mu2->GetMatchTrigger()+0.1*mu1->GetMatchTrigger();	    

				// mu1
				Double_t zr1 = mu1->Charge();
				Double_t pxr1 = mu1->Px();
				Double_t pyr1 = mu1->Py();
				Double_t pzr1 = mu1->Pz();
				Double_t ptr1 = TMath::Sqrt(pxr1*pxr1+pyr1*pyr1);
				Double_t er1 = mu1->E();
				Double_t momr1 = mu1->P();
				Double_t vzr1 = mu1->Zv();
				Double_t dcar1 = mu1->DCA();
				Double_t rabsr1 = mu1->GetRAtAbsorberEnd();
				Double_t thetaabsr1 = 180.*(1.-TMath::ATan(rabsr1/505.)/TMath::Pi());	// [deg]
				Double_t thetaabsr12 = (thetaabsr1*TMath::Pi()/180.)/2.;
				Double_t tanthetaabsr1 = TMath::Tan(thetaabsr12);
				Double_t etar1 = -TMath::Log(tanthetaabsr1);

				// mu2
				Double_t zr2 = mu2->Charge();
				Double_t pxr2 = mu2->Px();
				Double_t pyr2 = mu2->Py();
				Double_t pzr2 = mu2->Pz();
				Double_t ptr2 = TMath::Sqrt(pxr2*pxr2+pyr2*pyr2);
				Double_t er2 = mu2->E();
				Double_t momr2 = mu2->P();
				Double_t vzr2 = mu2->Zv();
				Double_t dcar2 = mu2->DCA();
				Double_t rabsr2 = mu1->GetRAtAbsorberEnd();
				Double_t thetaabsr2 = 180.*(1.-TMath::ATan(rabsr2/505.)/TMath::Pi());	// [deg]
				Double_t thetaabsr22 = (thetaabsr2*TMath::Pi()/180.)/2.;
				Double_t tanthetaabsr2 = TMath::Tan(thetaabsr22);
				Double_t etar2 = -TMath::Log(tanthetaabsr2);

				// for MC
				if(fReadMCInfo) {
					//mu1
					rabsr1 = 0;
					Double_t thetar1 = (180./TMath::Pi())*TMath::ATan2(TMath::Sqrt(pxr1*pxr1+pyr1*pyr1),pzr1);
					thetaabsr1 = thetar1;
					etar1 = mu1->Eta();
					//mu2
					rabsr2 = 0;
					Double_t thetar2 = (180./TMath::Pi())*TMath::ATan2(TMath::Sqrt(pxr2*pxr2+pyr2*pyr2),pzr2);
					thetaabsr2 = thetar2;
					etar2 = mu2->Eta();
				}
	
				if(TMath::Abs(etar1) > 8 || TMath::Abs(etar2) > 8) continue;

				// dimuon
				Double_t ptrec = TMath::Sqrt((pxr1+pxr2)*(pxr1+pxr2)+(pyr1+pyr2)*(pyr1+pyr2));
				Double_t raprec= Rap((er1+er2),(pzr1+pzr2));
				Double_t imassrec = Imass(er1,pxr1,pyr1,pzr1,er2,pxr2,pyr2,pzr2);
				if(!imassrec) continue;

				containerInput[1] = raprec ;   
				containerInput[2] = ptrec ;
				containerInput[3] = imassrec ;   
				containerInput[4] = trigCondition+0.05 ;
				containerInput[5] = ptr1 ;   
				containerInput[5] = ptr2 ;   
				containerInput[6] = momr1;
				containerInput[6] = momr2;
				containerInput[7] = trigside;
				containerInput[8] = rabsr1;
				containerInput[8] = rabsr2;
				containerInput[9] = zr1 + zr2;
				containerInput[10] = etar1;
				containerInput[10] = etar2;
				containerInput[11] = thetaabsr1;
				containerInput[11] = thetaabsr2;
				containerInput[12] = vzr1;
				containerInput[12] = vzr2;
				containerInput[13] = dcar1;
				containerInput[13] = dcar2;

	  		if (trigfired==1 && trigside==1) fCFManager->GetParticleContainer()->Fill(containerInput,4);
	  		if (trigfired==0 && trigside==1) fCFManager->GetParticleContainer()->Fill(containerInput,1);
				if(fReadMCInfo) fCFManager->GetParticleContainer()->Fill(containerInput,1); // Rec container

			} // second mu loop
		} // first mu loop
	} // end AOD-based ANALYSIS

	// end User analysis loop
}
//________________________________________________________________________
Double_t AliCFMuonResUpsilon::Imass(Float_t e1, Float_t px1, Float_t py1, Float_t pz1,
				   Float_t e2, Float_t px2, Float_t py2, Float_t pz2) const
{
// invariant mass calculation
    Float_t imass = (e1+e2)*(e1+e2)-((px1+px2)*(px1+px2)+(py1+py2)*(py1+py2)+(pz1+pz2)*(pz1+pz2));
		if(imass<0) return 0;
    Float_t imassrec = TMath::Sqrt((e1+e2)*(e1+e2)-((px1+px2)*(px1+px2)+(py1+py2)*(py1+py2)+(pz1+pz2)*(pz1+pz2)));
    return imassrec;
}
//________________________________________________________________________
Double_t AliCFMuonResUpsilon::Rap(Float_t e, Float_t pz) const
{
// calculate rapidity
    Float_t rap;
    if(e!=pz){
	rap = 0.5*TMath::Log((e+pz)/(e-pz));
	return rap;
    }
    else{
	rap = -200;
	return rap;
    }
}
//________________________________________________________________________
Double_t AliCFMuonResUpsilon::CostCS(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2,
Double_t Energy)
{
	// CS angle
  TLorentzVector pMu1CM, pMu2CM, pProjCM, pTargCM, pDimuCM; // In the CM. frame
  TLorentzVector pMu1Dimu, pMu2Dimu, pProjDimu, pTargDimu; // In the dimuon rest frame
  TVector3 beta,zaxisCS;
  Double_t mp=0.93827231;
  //
  // --- Fill the Lorentz vector for projectile and target in the CM frame
  //
  pProjCM.SetPxPyPzE(0.,0.,-Energy,TMath::Sqrt(Energy*Energy+mp*mp)); 
  pTargCM.SetPxPyPzE(0.,0.,Energy,TMath::Sqrt(Energy*Energy+mp*mp)); 
  //
  // --- Get the muons parameters in the CM frame 
  //
  pMu1CM.SetPxPyPzE(px1,py1,pz1,e1);
  pMu2CM.SetPxPyPzE(px2,py2,pz2,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  pDimuCM=pMu1CM+pMu2CM;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  beta=(-1./pDimuCM.E())*pDimuCM.Vect();
  pMu1Dimu=pMu1CM;
  pMu2Dimu=pMu2CM;
  pProjDimu=pProjCM;
  pTargDimu=pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the CS angle 
  //
  zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  //				     
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  Double_t cost;
  
  if(charge1>0) cost = zaxisCS.Dot((pMu1Dimu.Vect()).Unit());
  else cost = zaxisCS.Dot((pMu2Dimu.Vect()).Unit());
  return cost;
}
//________________________________________________________________________

//________________________________________________________________________
Double_t AliCFMuonResUpsilon::CostHE(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2)
{
	// Helicity 
  TLorentzVector pMu1CM, pMu2CM, pDimuCM; // In the CM frame 
  TLorentzVector pMu1Dimu, pMu2Dimu; // In the dimuon rest frame
  TVector3 beta,zaxisCS;
  //
  // --- Get the muons parameters in the CM frame
  //
  pMu1CM.SetPxPyPzE(px1,py1,pz1,e1);
  pMu2CM.SetPxPyPzE(px2,py2,pz2,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  pDimuCM=pMu1CM+pMu2CM;
  //
  // --- Translate the muon parameters in the dimuon rest frame
  //
  beta=(-1./pDimuCM.E())*pDimuCM.Vect();
  pMu1Dimu=pMu1CM;
  pMu2Dimu=pMu2CM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  //
  // --- Determine the z axis for the calculation of the polarization angle (i.e. the direction of the dimuon in the CM system)
  //
  TVector3 zaxis;
  zaxis=(pDimuCM.Vect()).Unit();
  
  // --- Calculation of the polarization angle (angle between mu+ and the z axis defined above)
  Double_t cost;
  if(charge1>0) cost = zaxis.Dot((pMu1Dimu.Vect()).Unit()); 
  else cost = zaxis.Dot((pMu2Dimu.Vect()).Unit()); 
  
  return cost;
}
//________________________________________________________________________

//________________________________________________________________________
Double_t AliCFMuonResUpsilon::PhiCS(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2,
Double_t Energy)
{
	// CS phi
   TLorentzVector pMu1CM, pMu2CM, pProjCM, pTargCM, pDimuCM; // In the CM frame
   TLorentzVector pMu1Dimu, pMu2Dimu, pProjDimu, pTargDimu; // In the dimuon rest frame
   TVector3 beta,yaxisCS, xaxisCS, zaxisCS;
   Double_t mp=0.93827231;
   //
   // --- Fill the Lorentz vector for projectile and target in the CM frame
   //
   pProjCM.SetPxPyPzE(0.,0.,-Energy,TMath::Sqrt(Energy*Energy+mp*mp)); 
   pTargCM.SetPxPyPzE(0.,0.,Energy,TMath::Sqrt(Energy*Energy+mp*mp)); 
   //
   // --- Get the muons parameters in the CM frame 
   //
   pMu1CM.SetPxPyPzE(px1,py1,pz1,e1);
   pMu2CM.SetPxPyPzE(px2,py2,pz2,e2);
   //
   // --- Obtain the dimuon parameters in the CM frame
   //
   pDimuCM=pMu1CM+pMu2CM;
   //
   // --- Translate the dimuon parameters in the dimuon rest frame
   //
   beta=(-1./pDimuCM.E())*pDimuCM.Vect();
   pMu1Dimu=pMu1CM;
   pMu2Dimu=pMu2CM;
   pProjDimu=pProjCM;
   pTargDimu=pTargCM;
   pMu1Dimu.Boost(beta);
   pMu2Dimu.Boost(beta);
   pProjDimu.Boost(beta);
   pTargDimu.Boost(beta);
   //
   // --- Determine the z axis for the CS angle 
   //
   zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
   yaxisCS=(((pProjDimu.Vect()).Unit()).Cross((pTargDimu.Vect()).Unit())).Unit();
   xaxisCS=(yaxisCS.Cross(zaxisCS)).Unit();
 
   Double_t phi;
   if(charge1>0) phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxisCS),((pMu1Dimu.Vect()).Dot(xaxisCS)));
   else phi = TMath::ATan2((pMu2Dimu.Vect()).Dot(yaxisCS),((pMu2Dimu.Vect()).Dot(xaxisCS)));
     
   return phi;
}
//________________________________________________________________________

//________________________________________________________________________
Double_t AliCFMuonResUpsilon::PhiHE(Double_t px1, Double_t py1, Double_t pz1, Double_t e1,
Double_t charge1, Double_t px2, Double_t py2, Double_t pz2, Double_t e2, Double_t Energy)
{
	// Helicity phi
   TLorentzVector pMu1CM, pMu2CM, pProjCM, pTargCM, pDimuCM; // In the CM frame 
   TLorentzVector pMu1Dimu, pMu2Dimu, pProjDimu, pTargDimu; // In the dimuon rest frame
   TVector3 beta,xaxis,yaxis,zaxis;
   Double_t mp=0.93827231;
 
   //
   // --- Get the muons parameters in the CM frame
   //
   pMu1CM.SetPxPyPzE(px1,py1,pz1,e1);
   pMu2CM.SetPxPyPzE(px2,py2,pz2,e2);
   //
   // --- Obtain the dimuon parameters in the CM frame
   //
   pDimuCM=pMu1CM+pMu2CM;
   //
   // --- Translate the muon parameters in the dimuon rest frame
   // 
   zaxis=(pDimuCM.Vect()).Unit();
 
   beta=(-1./pDimuCM.E())*pDimuCM.Vect();
 
   pProjCM.SetPxPyPzE(0.,0.,-Energy,TMath::Sqrt(Energy*Energy+mp*mp));
   pTargCM.SetPxPyPzE(0.,0.,Energy,TMath::Sqrt(Energy*Energy+mp*mp)); 
 
   pProjDimu=pProjCM;
   pTargDimu=pTargCM;
 
   pProjDimu.Boost(beta);
   pTargDimu.Boost(beta);
   
   yaxis=((pProjDimu.Vect()).Cross(pTargDimu.Vect())).Unit();
   xaxis=(yaxis.Cross(zaxis)).Unit();
   
   pMu1Dimu=pMu1CM;
   pMu2Dimu=pMu2CM;
   pMu1Dimu.Boost(beta);
   pMu2Dimu.Boost(beta);
 
   Double_t phi;
   if(charge1>0) phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxis),(pMu1Dimu.Vect()).Dot(xaxis));
   else phi = TMath::ATan2((pMu2Dimu.Vect()).Dot(yaxis),(pMu2Dimu.Vect()).Dot(xaxis));
   
   return phi;
}

//________________________________________________________________________
void AliCFMuonResUpsilon::Terminate(Option_t *) 
{
  // draw result of the Invariant mass MC and ESD

/*
		TH1D *h1 = dynamic_cast<TH1D*>(GetOutputData(1));
    AliCFContainer *cont = dynamic_cast<AliCFContainer*> (GetOutputData(2));   

    TH1D *kmass = cont->ShowProjection(3,0);
    TH1D *rmass = cont->ShowProjection(3,1);
		TH1D *mmass = cont->ShowProjection(3,4);

    TCanvas *c1 = new TCanvas("AliCFMuonResUpsilon","UPSILON Container",0,0,800,800);
    c1->Divide(2,2);
    c1->cd(1);
    kmass->Draw("HIST");
    c1->cd(2);
    rmass->Draw("HIST");
		c1->cd(3);
		mmass->Draw("HIST");
		c1->cd(4);
		h1->Draw("HIST");
		*/
}
//________________________________________________________________________

#endif
