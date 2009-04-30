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

//-----------------------------------------------------------------------
// Class for HF corrections as a function of many variables
// 6 Steps introduced: MC, MC Acc, Reco, Reco Acc, Reco Acc + ITS Cl, 
// Reco Acc + ITS Cl + PPR cuts
// 11 variables used: pt, y, cosThetaStar, ptPi, ptK, ct,
// dca, d0Pi, d0K, d0Pixd0K, cosPointingAngle
//
//-----------------------------------------------------------------------
// Author : C. Zampolli, CERN
//-----------------------------------------------------------------------

#include <TCanvas.h>
#include <TParticle.h>
#include <TH1I.h>
#include <TStyle.h>

#include "AliCFHeavyFlavourTaskMultiVarMultiStep.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliCFManager.h"
#include "AliCFContainer.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODMCParticle.h"

//__________________________________________________________________________
AliCFHeavyFlavourTaskMultiVarMultiStep::AliCFHeavyFlavourTaskMultiVarMultiStep() :
	AliAnalysisTaskSE(),
	fPDG(0),
	fCFManager(0x0),
	fHistEventsProcessed(0x0),
	fCountMC(0),
	fCountAcc(0),
	fCountReco(0),
	fCountRecoAcc(0),
	fCountRecoITSClusters(0),
	fCountRecoPPR(0),
	fEvents(0),
	fFillFromGenerated(kFALSE),
	fMinITSClusters(5)
{
	//
	//Default ctor
	//
}
//___________________________________________________________________________
AliCFHeavyFlavourTaskMultiVarMultiStep::AliCFHeavyFlavourTaskMultiVarMultiStep(const Char_t* name) :
	AliAnalysisTaskSE(name),
	fPDG(0),
	fCFManager(0x0),
	fHistEventsProcessed(0x0),
	fCountMC(0),
	fCountAcc(0),
	fCountReco(0),
	fCountRecoAcc(0),
	fCountRecoITSClusters(0),
	fCountRecoPPR(0),
	fEvents(0),
	fFillFromGenerated(kFALSE),
	fMinITSClusters(5)
{
	//
	// Constructor. Initialization of Inputs and Outputs
	//
	Info("AliCFHeavyFlavourTaskMultiVarMultiStep","Calling Constructor");
	/*
	  DefineInput(0) and DefineOutput(0)
	  are taken care of by AliAnalysisTaskSE constructor
	*/
	DefineOutput(1,TH1I::Class());
	DefineOutput(2,AliCFContainer::Class());
}

//___________________________________________________________________________
AliCFHeavyFlavourTaskMultiVarMultiStep& AliCFHeavyFlavourTaskMultiVarMultiStep::operator=(const AliCFHeavyFlavourTaskMultiVarMultiStep& c) 
{
	//
	// Assignment operator
	//
	if (this!=&c) {
		AliAnalysisTaskSE::operator=(c) ;
		fPDG      = c.fPDG;
		fCFManager  = c.fCFManager;
		fHistEventsProcessed = c.fHistEventsProcessed;
	}
	return *this;
}

//___________________________________________________________________________
AliCFHeavyFlavourTaskMultiVarMultiStep::AliCFHeavyFlavourTaskMultiVarMultiStep(const AliCFHeavyFlavourTaskMultiVarMultiStep& c) :
	AliAnalysisTaskSE(c),
	fPDG(c.fPDG),
	fCFManager(c.fCFManager),
	fHistEventsProcessed(c.fHistEventsProcessed),
	fCountMC(c.fCountMC),
	fCountAcc(c.fCountAcc),
	fCountReco(c.fCountReco),
	fCountRecoAcc(c.fCountRecoAcc),
	fCountRecoITSClusters(c.fCountRecoITSClusters),
	fCountRecoPPR(c.fCountRecoPPR),
	fEvents(c.fEvents),
	fFillFromGenerated(c.fFillFromGenerated),
	fMinITSClusters(c.fMinITSClusters)
{
	//
	// Copy Constructor
	//
}

//___________________________________________________________________________
AliCFHeavyFlavourTaskMultiVarMultiStep::~AliCFHeavyFlavourTaskMultiVarMultiStep() {
	//
	//destructor
	//
	if (fCFManager)           delete fCFManager ;
	if (fHistEventsProcessed) delete fHistEventsProcessed ;
}

//_________________________________________________
void AliCFHeavyFlavourTaskMultiVarMultiStep::UserExec(Option_t *)
{
	//
	// Main loop function
	//
	
	if (fFillFromGenerated){
		AliWarning("Flag to fill container with generated value ON ---> dca, d0pi, d0K, d0xd0, cosPointingAngle will be set as dummy!");
	}

	if (!fInputEvent) {
		Error("UserExec","NO EVENT FOUND!");
		return;
	}
	
	fEvents++;
	if (fEvents%10000 ==0) AliInfo(Form("Event %d",fEvents));
	AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
	fCFManager->SetEventInfo(aodEvent);
	
	// MC-event selection
	Double_t containerInput[11] ;
        
	//loop on the MC event
	
	TClonesArray* mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	if (!mcArray) AliError("Could not find Monte-Carlo in AOD");
	Int_t icountMC = 0;
	Int_t icountAcc = 0;
	Int_t icountReco = 0;
	Int_t icountRecoAcc = 0;
	Int_t icountRecoITSClusters = 0;
	Int_t icountRecoPPR = 0;
	
	for (Int_t iPart=0; iPart<mcArray->GetEntriesFast(); iPart++) { 
		AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(iPart));
		if (!mcPart) {
			AliWarning("Particle not found in tree, skipping"); 
			continue;
		} 
		
		// check the MC-level cuts
		if (!fCFManager->CheckParticleCuts(0, mcPart)) continue;  // 0 stands for MC level
		
		// fill the container for Gen-level selection
		Double_t vectorMC[6] = {9999.,9999.,9999.,9999.,9999.,9999.};
		if (GetGeneratedValuesFromMCParticle(mcPart, mcArray, vectorMC)){
			containerInput[0] = vectorMC[0];
			containerInput[1] = vectorMC[1] ;
			containerInput[2] = vectorMC[2] ;
			containerInput[3] = vectorMC[3] ;
			containerInput[4] = vectorMC[4] ;
			containerInput[5] = vectorMC[5] ;  // in micron
			containerInput[6] = 0.;    // dummy value, meaningless in MC, in micron
			containerInput[7] = 0.;   // dummy value, meaningless in MC, in micron
			containerInput[8] = 0.;   // dummy value, meaningless in MC, in micron
			containerInput[9] = -100000.; // dummy value, meaningless in MC, in micron^2
			containerInput[10] = 1.01;    // dummy value, meaningless in MC
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepGenerated);
			icountMC++;

			// check the MC-Acceptance level cuts
			// since standard CF functions are not applicable, using Kine Cuts on daughters
			
			Int_t daughter0 = mcPart->GetDaughter(0);
			Int_t daughter1 = mcPart->GetDaughter(1);
			AliDebug(2, Form("daughter0 = %d and daughter1 = %d",daughter0,daughter1));
			if (daughter0 == 0 || daughter1 == 0) {
				AliDebug(2, "Error! the D0 MC doesn't have correct daughters!! But it should have, this check was already done...");
			}
			if (TMath::Abs(daughter1 - daughter0 != 1)) {
				AliDebug(2, "The D0 MC doesn't come from a 2-prong decay, but it should be, this check was already done...");
			}
			AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter0));
			AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter1));
			if (!mcPartDaughter0 || !mcPartDaughter1) {
				AliWarning("At least one Daughter Particle not found in tree, but it should be, this check was already done..."); 
			}
			Double_t eta0 = mcPartDaughter0->Eta();
			Double_t eta1 = mcPartDaughter1->Eta();
			Double_t y0 = mcPartDaughter0->Y();
			Double_t y1 = mcPartDaughter1->Y();
			Double_t pt0 = mcPartDaughter0->Pt();
			Double_t pt1 = mcPartDaughter1->Pt();
			AliDebug(2, Form("Daughter 0: eta = %f, y = %f, pt = %f", eta0, y0, pt0));
			AliDebug(2, Form("Daughter 1: eta = %f, y = %f, pt = %f", eta1, y1, pt1));
			Bool_t daught0inAcceptance = (TMath::Abs(eta0) < 0.9 && pt0 > 0.1); 
			Bool_t daught1inAcceptance = (TMath::Abs(eta1) < 0.9 && pt1 > 0.1); 
			if (daught0inAcceptance && daught1inAcceptance) {
				// checking whether the cuts implemented in the CF are equivalent - simply a cross-check
				AliDebug(2, "Daughter particles in acceptance");
				if (!fCFManager->CheckParticleCuts(1, mcPartDaughter0)) {
					AliInfo("Inconsistency with CF cut for daughter 0!");
				} 
				if (!fCFManager->CheckParticleCuts(1, mcPartDaughter1)) {
					AliInfo("Inconsistency with CF cut for daughter 1!");
				} 
				fCFManager->GetParticleContainer()->Fill(containerInput,kStepAcceptance);
				icountAcc++;
			}
		}
		else {
			AliDebug(3,"Problems in filling the container");
			continue;
		}
 	}    
	
	AliDebug(2, Form("Found %i MC particles that are D0!!",icountMC));
	AliDebug(2, Form("Found %i MC particles that are D0 and satisfy Acc cuts!!",icountAcc));

	// Now go to rec level
	fCountMC += icountMC;
	fCountAcc += icountAcc;

	// AOD primary vertex
	AliAODVertex *vtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();

	// load heavy flavour vertices
	TClonesArray *arrayD0toKpi = (TClonesArray*)((aodEvent->GetList())->FindObject("D0toKpi")); 	
	if (!arrayD0toKpi) AliError("Could not find array of HF vertices");
	AliDebug(2, Form("Found %d vertices",arrayD0toKpi->GetEntriesFast()));
	
	for (Int_t iD0toKpi = 0; iD0toKpi<arrayD0toKpi->GetEntriesFast(); iD0toKpi++) {
		
		AliAODRecoDecayHF2Prong* d0tokpi = (AliAODRecoDecayHF2Prong*)arrayD0toKpi->At(iD0toKpi);
		Bool_t unsetvtx=kFALSE;
		if(!d0tokpi->GetOwnPrimaryVtx()) {
		  d0tokpi->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
		  unsetvtx=kTRUE;
		}

		// find associated MC particle
		Int_t mcLabel = d0tokpi->MatchToMC(421, mcArray) ;
		if (mcLabel == -1) 
			{
				AliDebug(2,"No MC particle found");
			}
		else {
			AliAODMCParticle* mcVtxHF = (AliAODMCParticle*)mcArray->At(mcLabel);
			if (!mcVtxHF) {
				AliWarning("Could not find associated MC in AOD MC tree");
				continue;
			}
					
			// check if associated MC v0 passes the cuts
			if (!fCFManager->CheckParticleCuts(0 ,mcVtxHF)) {  // 0 stands for MC
			        AliDebug(2, "Skipping the particles due to cuts");
				continue; 
			}


			// fill the container...
			Double_t pt = d0tokpi->Pt();
			Double_t rapidity = d0tokpi->YD0();
			
			Double_t cosThetaStar = 9999.;
			Double_t pTpi = 0.;
			Double_t pTK = 0.;
			Double_t dca = d0tokpi->GetDCA();
			Double_t d0pi = 0.;
			Double_t d0K = 0.;
			Double_t d0xd0 = d0tokpi->Prodd0d0();
			Double_t cosPointingAngle = d0tokpi->CosPointingAngle();
			Int_t pdgCode = mcVtxHF->GetPdgCode();
			if (pdgCode > 0){
				cosThetaStar = d0tokpi->CosThetaStarD0();
				pTpi = d0tokpi->PtProng(0);
				pTK = d0tokpi->PtProng(1);
				d0pi = d0tokpi->Getd0Prong(0);
				d0K = d0tokpi->Getd0Prong(1);
			}
			else {
				cosThetaStar = d0tokpi->CosThetaStarD0bar();
				pTpi = d0tokpi->PtProng(1);
				pTK = d0tokpi->PtProng(0);
				d0pi = d0tokpi->Getd0Prong(1);
				d0K = d0tokpi->Getd0Prong(0);
			}

			Double_t cT = d0tokpi->CtD0();

			if (!fFillFromGenerated){
				// ...either with reconstructed values....
				containerInput[0] = pt;
				containerInput[1] = rapidity;
				containerInput[2] = cosThetaStar;
				containerInput[3] = pTpi;
				containerInput[4] = pTK;
				containerInput[5] = cT*1.E4;  // in micron
				containerInput[6] = dca*1.E4;  // in micron
				containerInput[7] = d0pi*1.E4;  // in micron
				containerInput[8] = d0K*1.E4;  // in micron
				containerInput[9] = d0xd0*1.E8;  // in micron^2
				containerInput[10] = cosPointingAngle;  // in micron
			}
			else {
				// ... or with generated values				
				Double_t vectorMC[6] = {9999.,9999.,9999.,9999.,9999.,9999.};
				if (GetGeneratedValuesFromMCParticle(mcVtxHF, mcArray, vectorMC)){
					containerInput[0] = vectorMC[0];
					containerInput[1] = vectorMC[1] ;
					containerInput[2] = vectorMC[2] ;
					containerInput[3] = vectorMC[3] ;
					containerInput[4] = vectorMC[4] ;
					containerInput[5] = vectorMC[5] ;  // in micron
					containerInput[6] = 0.;    // dummy value, meaningless in MC, in micron
					containerInput[7] = 0.;   // dummy value, meaningless in MC, in micron
					containerInput[8] = 0.;   // dummy value, meaningless in MC, in micron
					containerInput[9] = 100000.; // dummy value, meaningless in MC, in micron^2
					containerInput[10] = 1.01;    // dummy value, meaningless in MC
				}
				else {
					AliDebug(3,"Problems in filling the container");
					continue;
				}
			}
			AliDebug(2, Form("Filling the container with pt = %f, rapidity = %f, cosThetaStar = %f, pTpi = %f, pTK = %f, cT = %f", containerInput[0], containerInput[1], containerInput[2], containerInput[3], containerInput[4], containerInput[5]));
			icountReco++;
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructed) ;   

			// cut in acceptance
			Bool_t acceptanceProng0 = (TMath::Abs(d0tokpi->EtaProng(0))< 0.9 && d0tokpi->PtProng(0) > 0.1);
			Bool_t acceptanceProng1 = (TMath::Abs(d0tokpi->EtaProng(1))< 0.9 && d0tokpi->PtProng(1) > 0.1);
			if (acceptanceProng0 && acceptanceProng1) {
				AliDebug(2,"D0 reco daughters in acceptance");
				fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoAcceptance) ;
				icountRecoAcc++;   

				// cut on the min n. of clusters in ITS
				Int_t ncls0=0;
				for(Int_t l=0;l<6;l++) if(TESTBIT(d0tokpi->GetITSClusterMap(),l)) ncls0++;
				AliDebug(2, Form("n clusters = %d", ncls0));
				if (ncls0 >= fMinITSClusters){
					fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoITSClusters) ;
					icountRecoITSClusters++;   
					AliDebug(2,Form("pT = %f, dca = %f, cosThetaStar = %f, pTpi = %f, pTK = %f, d0pi = %f, d0K = %f, d0xd0 = %f, cosPointingAngle = %f", pt, dca, cosThetaStar,pTpi, pTK, d0pi*1E4, d0K*1E4, d0xd0*1E8, cosPointingAngle));

					// PPR cuts 
					Double_t cuts[6];
					if (pt <= 1){
						cuts[0] = 400;
						cuts[1] = 0.8;
						cuts[2] = 0.5;
						cuts[3] = 500;
						cuts[4] = -20000;
						cuts[5] = 0.5;
					}
					else if (pt > 1 && pt <= 2){
						cuts[0] = 300;
						cuts[1] = 0.8;
						cuts[2] = 0.6;
						cuts[3] = 500;
						cuts[4] = -20000;
						cuts[5] = 0.6;
					}
					else if (pt > 2 && pt <= 3){
						cuts[0] = 200;
						cuts[1] = 0.8;
						cuts[2] = 0.7;
						cuts[3] = 500;
						cuts[4] = -20000;
						cuts[5] = 0.8;
					}
					else if (pt > 3 && pt <= 5){
						cuts[0] = 200;
						cuts[1] = 0.8;
						cuts[2] = 0.7;
						cuts[3] = 500;
						cuts[4] = -10000;
						cuts[5] = 0.8;
					}
					else if (pt > 5){
						cuts[0] = 200;
						cuts[1] = 0.8;
						cuts[2] = 0.7;
						cuts[3] = 500;
						cuts[4] = -5000;
						cuts[5] = 0.8;
					}
					if (dca*1E4 < cuts[0] 
					    && TMath::Abs(cosThetaStar) < cuts[1]  
					    && pTpi > cuts[2] 
					    && pTK > cuts[2]  
					    && TMath::Abs(d0pi*1E4) < cuts[3] 
					    && TMath::Abs(d0K*1E4) < cuts[3]  
					    && d0xd0*1E8 < cuts[4] 
					    && cosPointingAngle > cuts[5]
					    ){

						AliDebug(2,"Particle passed PPR cuts");
						fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoPPR) ;   
						icountRecoPPR++;
					}
					else{
						AliDebug(2,"Particle skipped due to PPR cuts");
						if (dca*1E4 > cuts[0]){
							AliDebug(2,"Problems with dca");
						}
						if ( TMath::Abs(cosThetaStar) > cuts[1]){
							AliDebug(2,"Problems with cosThetaStar");
						}
						if (pTpi < cuts[2]){
							AliDebug(2,"Problems with pTpi");
						}
						if (pTK < cuts[2]){
							AliDebug(2,"Problems with pTK");
						}
						if (TMath::Abs(d0pi*1E4) > cuts[3]){
							AliDebug(2,"Problems with d0pi");
						}
						if (TMath::Abs(d0K*1E4) > cuts[3]){
							AliDebug(2,"Problems with d0K");
						}
						if (d0xd0*1E8 > cuts[4]){
							AliDebug(2,"Problems with d0xd0");
						}
						if (cosPointingAngle < cuts[5]){
							AliDebug(2,"Problems with cosPointingAngle");
						}
					}
				}
			}
		}
		if(unsetvtx) d0tokpi->UnsetOwnPrimaryVtx();
	} // end loop on D0->Kpi

	AliDebug(2, Form("Found %i Reco particles that are D0!!",icountReco));

	fCountReco+= icountReco;
	fCountRecoAcc+= icountRecoAcc;
	fCountRecoITSClusters+= icountRecoITSClusters;
	fCountRecoPPR+= icountRecoPPR;
	
	fHistEventsProcessed->Fill(0);
	/* PostData(0) is taken care of by AliAnalysisTaskSE */
	PostData(1,fHistEventsProcessed) ;
	PostData(2,fCFManager->GetParticleContainer()) ;
}


//___________________________________________________________________________
void AliCFHeavyFlavourTaskMultiVarMultiStep::Terminate(Option_t*)
{
	// The Terminate() function is the last function to be called during
	// a query. It always runs on the client, it can be used to present
	// the results graphically or save the results to file.
	
	Info("Terminate","");
	AliAnalysisTaskSE::Terminate();
	
	AliInfo(Form("Found %i MC particles that are D0 in MC, in %d events",fCountMC,fEvents));
	AliInfo(Form("Found %i MC particles that are D0 in MC and satisfy Acc cuts, in %d events",fCountAcc,fEvents));
	AliInfo(Form("Found %i reco D0 that are decaying in K+pi, in %d events",fCountReco,fEvents));
	AliInfo(Form("Among the above, found %i reco D0 that are decaying in K+pi and are in the requested acceptance, in %d events",fCountRecoAcc,fEvents));
	AliInfo(Form("Among the above, found %i reco D0 that are decaying in K+pi and have at least %d clusters in ITS, in %d events",fCountRecoITSClusters,fMinITSClusters,fEvents));
	AliInfo(Form("Among the above, found %i reco D0 that are decaying in K+pi and satisfy PPR cuts, in %d events",fCountRecoPPR,fEvents));
	
	// draw some example plots....
	
	AliCFContainer *cont= dynamic_cast<AliCFContainer*> (GetOutputData(2));
	
	// projecting the containers to obtain histograms
	// first argument = variable, second argument = step

	// MC-level
	TH1D* h00 =   cont->ShowProjection(0,0) ;   // pt
	TH1D* h10 =   cont->ShowProjection(1,0) ;   // rapidity
	TH1D* h20 =   cont->ShowProjection(2,0) ;   // cosThetaStar
	TH1D* h30 =   cont->ShowProjection(3,0) ;   // pTpi
	TH1D* h40 =   cont->ShowProjection(4,0) ;   // pTK
	TH1D* h50 =   cont->ShowProjection(5,0) ;   // cT
	TH1D* h60 =   cont->ShowProjection(6,0) ;   // dca
	TH1D* h70 =   cont->ShowProjection(7,0) ;   // d0pi
	TH1D* h80 =   cont->ShowProjection(8,0) ;   // d0K
	TH1D* h90 =   cont->ShowProjection(9,0) ;   // d0xd0
	TH1D* h100 =   cont->ShowProjection(10,0) ;   // cosPointingAngle
	
	// MC-Acceptance level
	TH1D* h01 =   cont->ShowProjection(0,1) ;   // pt
	TH1D* h11 =   cont->ShowProjection(1,1) ;   // rapidity
	TH1D* h21 =   cont->ShowProjection(2,1) ;   // cosThetaStar
	TH1D* h31 =   cont->ShowProjection(3,1) ;   // pTpi
	TH1D* h41 =   cont->ShowProjection(4,1) ;   // pTK
	TH1D* h51 =   cont->ShowProjection(5,1) ;   // cT
	TH1D* h61 =   cont->ShowProjection(6,1) ;   // dca
	TH1D* h71 =   cont->ShowProjection(7,1) ;   // d0pi
	TH1D* h81 =   cont->ShowProjection(8,1) ;   // d0K
	TH1D* h91 =   cont->ShowProjection(9,1) ;   // d0xd0
	TH1D* h101 =   cont->ShowProjection(10,1) ;   // cosPointingAngle

	// Reco-level
	TH1D* h02 =   cont->ShowProjection(0,2) ;   // pt
	TH1D* h12 =   cont->ShowProjection(1,2) ;   // rapidity
	TH1D* h22 =   cont->ShowProjection(2,2) ;   // cosThetaStar
	TH1D* h32 =   cont->ShowProjection(3,2) ;   // pTpi
	TH1D* h42 =   cont->ShowProjection(4,2) ;   // pTK
	TH1D* h52 =   cont->ShowProjection(5,2) ;   // cT
	TH1D* h62 =   cont->ShowProjection(6,2) ;   // dca
	TH1D* h72 =   cont->ShowProjection(7,2) ;   // d0pi
	TH1D* h82 =   cont->ShowProjection(8,2) ;   // d0K
	TH1D* h92 =   cont->ShowProjection(9,2) ;   // d0xd0
	TH1D* h102 =   cont->ShowProjection(10,2) ;   // cosPointingAngle
	
	h00->SetTitle("pT_D0 (GeV/c)");
	h10->SetTitle("rapidity");
	h20->SetTitle("cosThetaStar");
	h30->SetTitle("pT_pi (GeV/c)");
	h40->SetTitle("pT_K (Gev/c)");
	h50->SetTitle("cT (#mum)");
	h60->SetTitle("dca (#mum)");
	h70->SetTitle("d0_pi (#mum)");
	h80->SetTitle("d0_K (#mum)");
	h90->SetTitle("d0xd0 (#mum^2)");
	h100->SetTitle("cosPointingAngle");

	h00->GetXaxis()->SetTitle("pT_D0 (GeV/c)");
	h10->GetXaxis()->SetTitle("rapidity");
	h20->GetXaxis()->SetTitle("cosThetaStar");
	h30->GetXaxis()->SetTitle("pT_pi (GeV/c)");
	h40->GetXaxis()->SetTitle("pT_K (Gev/c)");
	h50->GetXaxis()->SetTitle("cT (#mum)");
	h60->GetXaxis()->SetTitle("dca (#mum)");
	h70->GetXaxis()->SetTitle("d0_pi (#mum)");
	h80->GetXaxis()->SetTitle("d0_K (#mum)");
	h90->GetXaxis()->SetTitle("d0xd0 (#mum^2)");
	h100->GetXaxis()->SetTitle("cosPointingAngle");

	h01->SetTitle("pT_D0 (GeV/c)");
	h11->SetTitle("rapidity");
	h21->SetTitle("cosThetaStar");
	h31->SetTitle("pT_pi (GeV/c)");
	h41->SetTitle("pT_K (Gev/c)");
	h51->SetTitle("cT (#mum)");
	h61->SetTitle("dca (#mum)");
	h71->SetTitle("d0_pi (#mum)");
	h81->SetTitle("d0_K (#mum)");
	h91->SetTitle("d0xd0 (#mum^2)");
	h101->SetTitle("cosPointingAngle");

	h01->GetXaxis()->SetTitle("pT_D0 (GeV/c)");
	h11->GetXaxis()->SetTitle("rapidity");
	h21->GetXaxis()->SetTitle("cosThetaStar");
	h31->GetXaxis()->SetTitle("pT_pi (GeV/c)");
	h41->GetXaxis()->SetTitle("pT_K (Gev/c)");
	h51->GetXaxis()->SetTitle("cT (#mum)");
	h61->GetXaxis()->SetTitle("dca (#mum)");
	h71->GetXaxis()->SetTitle("d0_pi (#mum)");
	h81->GetXaxis()->SetTitle("d0_K (#mum)");
	h91->GetXaxis()->SetTitle("d0xd0 (#mum^2)");
	h101->GetXaxis()->SetTitle("cosPointingAngle");

	h02->SetTitle("pT_D0 (GeV/c)");
	h12->SetTitle("rapidity");
	h22->SetTitle("cosThetaStar");
	h32->SetTitle("pT_pi (GeV/c)");
	h42->SetTitle("pT_K (Gev/c)");
	h52->SetTitle("cT (#mum)");
	h62->SetTitle("dca (#mum)");
	h72->SetTitle("d0_pi (#mum)");
	h82->SetTitle("d0_K (#mum)");
	h92->SetTitle("d0xd0 (#mum^2)");
	h102->SetTitle("cosPointingAngle");

	h02->GetXaxis()->SetTitle("pT_D0 (GeV/c)");
	h12->GetXaxis()->SetTitle("rapidity");
	h22->GetXaxis()->SetTitle("cosThetaStar");
	h32->GetXaxis()->SetTitle("pT_pi (GeV/c)");
	h42->GetXaxis()->SetTitle("pT_K (Gev/c)");
	h52->GetXaxis()->SetTitle("cT (#mum)");
	h62->GetXaxis()->SetTitle("dca (#mum)");
	h72->GetXaxis()->SetTitle("d0_pi (#mum)");
	h82->GetXaxis()->SetTitle("d0_K (#mum)");
	h92->GetXaxis()->SetTitle("d0xd0 (#mum^2)");
	h102->GetXaxis()->SetTitle("cosPointingAngle");

	Double_t max0 = h00->GetMaximum();
	Double_t max1 = h10->GetMaximum();
	Double_t max2 = h20->GetMaximum();
	Double_t max3 = h30->GetMaximum();
	Double_t max4 = h40->GetMaximum();
	Double_t max5 = h50->GetMaximum();
	Double_t max6 = h60->GetMaximum();
	Double_t max7 = h70->GetMaximum();
	Double_t max8 = h80->GetMaximum();
	Double_t max9 = h90->GetMaximum();
	Double_t max10 = h100->GetMaximum();
	
	h00->GetYaxis()->SetRangeUser(0,max0*1.2);
	h10->GetYaxis()->SetRangeUser(0,max1*1.2);
	h20->GetYaxis()->SetRangeUser(0,max2*1.2);
	h30->GetYaxis()->SetRangeUser(0,max3*1.2);
	h40->GetYaxis()->SetRangeUser(0,max4*1.2);
	h50->GetYaxis()->SetRangeUser(0,max5*1.2);
	h60->GetYaxis()->SetRangeUser(0,max6*1.2);
	h70->GetYaxis()->SetRangeUser(0,max7*1.2);
	h80->GetYaxis()->SetRangeUser(0,max8*1.2);
	h90->GetYaxis()->SetRangeUser(0,max9*1.2);
	h100->GetYaxis()->SetRangeUser(0,max10*1.2);
	
	h01->GetYaxis()->SetRangeUser(0,max0*1.2);
	h11->GetYaxis()->SetRangeUser(0,max1*1.2);
	h21->GetYaxis()->SetRangeUser(0,max2*1.2);
	h31->GetYaxis()->SetRangeUser(0,max3*1.2);
	h41->GetYaxis()->SetRangeUser(0,max4*1.2);
	h51->GetYaxis()->SetRangeUser(0,max5*1.2);
	h61->GetYaxis()->SetRangeUser(0,max6*1.2);
	h71->GetYaxis()->SetRangeUser(0,max7*1.2);
	h81->GetYaxis()->SetRangeUser(0,max8*1.2);
	h91->GetYaxis()->SetRangeUser(0,max9*1.2);
	h101->GetYaxis()->SetRangeUser(0,max10*1.2);
	
	h00->SetMarkerStyle(20);
	h10->SetMarkerStyle(24);
	h20->SetMarkerStyle(21);
	h30->SetMarkerStyle(25);
	h40->SetMarkerStyle(27);
	h50->SetMarkerStyle(28);
	h60->SetMarkerStyle(20);
	h70->SetMarkerStyle(24);
	h80->SetMarkerStyle(21);
	h90->SetMarkerStyle(25);
	h100->SetMarkerStyle(27);

	h00->SetMarkerColor(2);
	h10->SetMarkerColor(2);
	h20->SetMarkerColor(2);
	h30->SetMarkerColor(2);
	h40->SetMarkerColor(2);
	h50->SetMarkerColor(2);
 	h60->SetMarkerColor(2);
	h70->SetMarkerColor(2);
	h80->SetMarkerColor(2);
	h90->SetMarkerColor(2);
	h100->SetMarkerColor(2);

	h01->SetMarkerStyle(20) ;
	h11->SetMarkerStyle(24) ;
	h21->SetMarkerStyle(21) ;
	h31->SetMarkerStyle(25) ;
	h41->SetMarkerStyle(27) ;
	h51->SetMarkerStyle(28) ;
	h61->SetMarkerStyle(20);
	h71->SetMarkerStyle(24);
	h81->SetMarkerStyle(21);
	h91->SetMarkerStyle(25);
	h101->SetMarkerStyle(27);

	h01->SetMarkerColor(8);
	h11->SetMarkerColor(8);
	h21->SetMarkerColor(8);
	h31->SetMarkerColor(8);
	h41->SetMarkerColor(8);
	h51->SetMarkerColor(8);
 	h61->SetMarkerColor(8);
	h71->SetMarkerColor(8);
	h81->SetMarkerColor(8);
	h91->SetMarkerColor(8);
	h101->SetMarkerColor(8);

	h02->SetMarkerStyle(20) ;
	h12->SetMarkerStyle(24) ;
	h22->SetMarkerStyle(21) ;
	h32->SetMarkerStyle(25) ;
	h42->SetMarkerStyle(27) ;
	h52->SetMarkerStyle(28) ;
	h62->SetMarkerStyle(20);
	h72->SetMarkerStyle(24);
	h82->SetMarkerStyle(21);
	h92->SetMarkerStyle(25);
	h102->SetMarkerStyle(27);

	h02->SetMarkerColor(4);
	h12->SetMarkerColor(4);
	h22->SetMarkerColor(4);
	h32->SetMarkerColor(4);
	h42->SetMarkerColor(4);
	h52->SetMarkerColor(4);
	h62->SetMarkerColor(4);
	h72->SetMarkerColor(4);
	h82->SetMarkerColor(4);
	h92->SetMarkerColor(4);
	h102->SetMarkerColor(4);
	
	gStyle->SetCanvasColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetStatColor(0);

	// drawing in 2 separate canvas for a matter of clearity
	TCanvas * c1 =new TCanvas("c1","pT, rapidiy, cosThetaStar",1100,1600);
	c1->Divide(3,3);
	
	c1->cd(1);
	h00->Draw("p");
	c1->cd(1);
	c1->cd(2);
	h01->Draw("p");
	c1->cd(2);
	c1->cd(3);
	h02->Draw("p");
	c1->cd(3);
	c1->cd(4);
	h10->Draw("p");
	c1->cd(4);
	c1->cd(5);
	h11->Draw("p");
	c1->cd(5);
	c1->cd(6);
	h12->Draw("p");
	c1->cd(6);
	c1->cd(7);
	h20->Draw("p");
	c1->cd(7);
	c1->cd(8);
	h21->Draw("p");
	c1->cd(8);
	c1->cd(9);
	h22->Draw("p");
	c1->cd(9);
	c1->cd();

	TCanvas * c2 =new TCanvas("c2","pTpi, pTK, cT",1100,1600);
	c2->Divide(3,3);
	c2->cd(1);
	h30->Draw("p");
	c2->cd(1);
	c2->cd(2);
	h31->Draw("p");
	c2->cd(2);
	c2->cd(3);
	h32->Draw("p");
	c2->cd(3);
	c2->cd(4);
	h40->Draw("p");
	c2->cd(4);
	c2->cd(5);
	h41->Draw("p");
	c2->cd(5);
	c2->cd(6);
	h42->Draw("p");
	c2->cd(6);
	c2->cd(7);
	h50->Draw("p");
	c2->cd(7);
	c2->cd(8);
	h51->Draw("p");
	c2->cd(8);
	c2->cd(9);
	h52->Draw("p");
	c2->cd(9);
	c2->cd();

	TCanvas * c3 =new TCanvas("c3","dca, d0pi, d0K",1100,1600);
	c3->Divide(3,3);
	c3->cd(1);
	h60->Draw("p");
	c3->cd(1);
	c3->cd(2);
	h61->Draw("p");
	c3->cd(2);
	c3->cd(3);
	h62->Draw("p");
	c3->cd(3);
	c3->cd(4);
	h70->Draw("p");
	c3->cd(4);
	c3->cd(5);
	h71->Draw("p");
	c3->cd(5);
	c3->cd(6);
	h72->Draw("p");
	c3->cd(6);
	c3->cd(7);
	h80->Draw("p");
	c3->cd(7);
	c3->cd(8);
	h81->Draw("p");
	c3->cd(8);
	c3->cd(9);
	h82->Draw("p");
	c3->cd(9);
	c3->cd();

	TCanvas * c4 =new TCanvas("c4","d0xd0, cosPointingAngle",1100,770);
	c4->Divide(3,2);
	c4->cd(1);
	h90->Draw("p");
	c4->cd(1);
	c4->cd(2);
	h91->Draw("p");
	c4->cd(2);
	c4->cd(3);
	h92->Draw("p");
	c4->cd(3);
	c4->cd(4);
	h100->Draw("p");
	c4->cd(4);
	c4->cd(5);
	h101->Draw("p");
	c4->cd(5);
	c4->cd(6);
	h102->Draw("p");
	c4->cd(6);
	c4->cd();

	/*
	c1->SaveAs("Plots/pT_rapidity_cosThetaStar.eps");
	c2->SaveAs("Plots/pTpi_pTK_cT.eps");
	c3->SaveAs("Plots/dca_d0pi_d0TK.eps");
	c4->SaveAs("Plots/d0xd0_cosPointingAngle.eps");

	c1->SaveAs("Plots/pT_rapidity_cosThetaStar.gif");
	c2->SaveAs("Plots/pTpi_pTK_cT.gif");
	c3->SaveAs("Plots/dca_d0pi_d0TK.gif");
	c4->SaveAs("Plots/d0xd0_cosPointingAngle.gif");
	*/
}

//___________________________________________________________________________
void AliCFHeavyFlavourTaskMultiVarMultiStep::UserCreateOutputObjects() {
	//HERE ONE CAN CREATE OUTPUT OBJECTS, IN PARTICULAR IF THE OBJECT PARAMETERS DON'T NEED
	//TO BE SET BEFORE THE EXECUTION OF THE TASK
	//
	Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
	
	//slot #1
	OpenFile(1);
	fHistEventsProcessed = new TH1I("fHistEventsProcessed","",1,0,1) ;
}

//___________________________________________________________________________
Double_t AliCFHeavyFlavourTaskMultiVarMultiStep::CosThetaStar(AliAODMCParticle* mcPart, AliAODMCParticle* mcPartDaughter0, AliAODMCParticle* mcPartDaughter1) const {

	//
	// to calculate cos(ThetaStar) for generated particle
	// using the K, since mcPartDaughter0 and mcPartDaughter1 always correspond to K and pi respectively 
	// (see where the function is called)
	//

	Int_t pdgvtx = mcPart->GetPdgCode();
	/*	if (pdgvtx > 0) { // setting as the first daughter always the kaon, to be used to calculate cos(ThetaStar)
		Int_t pdgprong0 = TMath::Abs(mcPartDaughter0->GetPdgCode());
		Int_t pdgprong1 = TMath::Abs(mcPartDaughter1->GetPdgCode());
		AliInfo(Form("D0, with pdgprong0 = %d, pdgprong1 = %d",pdgprong0,pdgprong1));
		AliDebug(2,"This is a D0");
		AliAODMCParticle* mcPartdummy = mcPartDaughter0;
		mcPartDaughter0 = mcPartDaughter1;
		mcPartDaughter1 = mcPartdummy;
	} 
	else{
		AliInfo("D0bar");
	}
	*/
	Int_t pdgprong0 = TMath::Abs(mcPartDaughter0->GetPdgCode());
	Int_t pdgprong1 = TMath::Abs(mcPartDaughter1->GetPdgCode());
	if (pdgvtx > 0) { // setting as the first daughter always the kaon, to be used to calculate cos(ThetaStar)
		AliDebug(2,"D0");
	}
	else{
		AliDebug(2,"D0bar");
	}
	if (pdgprong0 == 211){
		AliDebug(2,Form("pdgprong0 = %d, pdgprong1 = %d, switching...",pdgprong0,pdgprong1));
		AliAODMCParticle* mcPartdummy = mcPartDaughter0;
		mcPartDaughter0 = mcPartDaughter1;
		mcPartDaughter1 = mcPartdummy;
		pdgprong0 = TMath::Abs(mcPartDaughter0->GetPdgCode());
		pdgprong1 = TMath::Abs(mcPartDaughter1->GetPdgCode());
	} 

	AliDebug(2,Form("After checking, pdgprong0 = %d, pdgprong1 = %d",pdgprong0,pdgprong1));
	Double_t massvtx = TDatabasePDG::Instance()->GetParticle(TMath::Abs(pdgvtx))->Mass();
	Double_t massp[2];
	massp[0] = TDatabasePDG::Instance()->GetParticle(pdgprong0)->Mass();
	massp[1] = TDatabasePDG::Instance()->GetParticle(pdgprong1)->Mass();

	Double_t pStar = TMath::Sqrt(TMath::Power(massvtx*massvtx-massp[0]*massp[0]-massp[1]*massp[1],2.)-4.*massp[0]*massp[0]*massp[1]*massp[1])/(2.*massvtx);

	Double_t px = mcPartDaughter0->Px()+mcPartDaughter1->Px();
	Double_t py = mcPartDaughter0->Py()+mcPartDaughter1->Py();
	Double_t pz = mcPartDaughter0->Pz()+mcPartDaughter1->Pz();
	Double_t p = TMath::Sqrt(px*px+py*py+pz*pz);
	Double_t e =  TMath::Sqrt(massvtx*massvtx+p*p);
	Double_t beta = p/e;
	Double_t gamma = e/massvtx;
	TVector3 mom(mcPartDaughter0->Px(),mcPartDaughter0->Py(),mcPartDaughter0->Pz());
	TVector3 momTot(mcPartDaughter0->Px()+mcPartDaughter1->Px(),mcPartDaughter0->Py()+mcPartDaughter1->Py(),mcPartDaughter0->Pz()+mcPartDaughter1->Pz());

	Double_t qlprong = mom.Dot(momTot)/momTot.Mag();  // analog to AliAODRecoDecay::QlProng(0)
	
	AliDebug(2,Form("pStar = %f, beta = %f, gamma = %f, qlprong = %f, massp[0] = %f", pStar, beta, gamma, qlprong, massp[0]));
	Double_t cts = (qlprong/gamma-beta*TMath::Sqrt(pStar*pStar+massp[0]*massp[0]))/pStar;
	AliDebug(2,Form("cts = %f", cts));
	return cts;
}
//___________________________________________________________________________
Double_t AliCFHeavyFlavourTaskMultiVarMultiStep::CT(AliAODMCParticle* mcPart, AliAODMCParticle* mcPartDaughter0, AliAODMCParticle* mcPartDaughter1) const {

	//
	// to calculate cT for generated particle
	//

	Double_t xmcPart[3] = {0,0,0};
	Double_t xdaughter0[3] = {0,0,0};
	Double_t xdaughter1[3] = {0,0,0};
	mcPart->XvYvZv(xmcPart);  // cm
	mcPartDaughter0->XvYvZv(xdaughter0);  // cm
	mcPartDaughter1->XvYvZv(xdaughter1);  //cm
	Double_t prodVtxD0 = TMath::Sqrt(xmcPart[0]*xmcPart[0]+
					 xmcPart[1]*xmcPart[1]+
					 xmcPart[2]*xmcPart[2]);
	Double_t prodVtxDaughter0 = TMath::Sqrt(xdaughter0[0]*xdaughter0[0]+
						xdaughter0[1]*xdaughter0[1]+
						xdaughter0[2]*xdaughter0[2]);
	Double_t prodVtxDaughter1 = TMath::Sqrt(xdaughter1[0]*xdaughter1[0]+
						xdaughter1[1]*xdaughter1[1]+
						xdaughter1[2]*xdaughter1[2]);
 
	AliDebug(2, Form("D0:        x = %f, y = %f, z = %f, production vertex distance = %f (cm), %f (micron)", xmcPart[0], xmcPart[1], xmcPart[2], prodVtxD0, prodVtxD0*1.E4));
	AliDebug(2, Form("Daughter0: x = %f, y = %f, z = %f, production vertex distance = %f (cm) %f (micron)", xdaughter0[0], xdaughter0[1], xdaughter0[2], prodVtxDaughter0, prodVtxDaughter0*1E4));
	AliDebug(2, Form("Daughter1: x = %f, y = %f, z = %f, production vertex distance = %f (cm) %f (micron)", xdaughter1[0], xdaughter1[1], xdaughter1[2], prodVtxDaughter1, prodVtxDaughter1*1.E4));

	Double_t cT0 = TMath::Sqrt((xdaughter0[0]-xmcPart[0])*(xdaughter0[0]-xmcPart[0])+
				     (xdaughter0[1]-xmcPart[1])*(xdaughter0[1]-xmcPart[1])+
				     (xdaughter0[2]-xmcPart[2])*(xdaughter0[2]-xmcPart[2]));

	Double_t cT1 = TMath::Sqrt((xdaughter1[0]-xmcPart[0])*(xdaughter1[0]-xmcPart[0])+
				     (xdaughter1[1]-xmcPart[1])*(xdaughter1[1]-xmcPart[1])+
				     (xdaughter1[2]-xmcPart[2])*(xdaughter1[2]-xmcPart[2]));

	if (cT0 != cT1) {
		AliWarning("cT from daughter 0 different from cT from daughter 1! Using cT from daughter 0, but PLEASE, CHECK....");
	}

	// calculating cT from cT0

	Double_t mass = TDatabasePDG::Instance()->GetParticle(mcPart->GetPdgCode())->Mass(); // mcPart->GetPdgCode() should return +/- 421 for the D0/D0bar
	Double_t p = mcPart-> P();
	Double_t cT = cT0*mass/p;
	AliDebug(2, Form("cT from daughter 0 = %f (micron)", cT0*1E4)); 
	AliDebug(2, Form("cT from daughter 1 = %f (micron)", cT1*1E4)); 
	AliDebug(2, Form("cT (with mass = %f and p = %f) = %f (micron)", mass, p, cT*1E4));
	return cT;
}
//________________________________________________________________________________
Bool_t AliCFHeavyFlavourTaskMultiVarMultiStep::GetGeneratedValuesFromMCParticle(AliAODMCParticle* mcPart, TClonesArray* mcArray, Double_t* vectorMC) const {

	// 
	// collecting all the necessary info (pt, y, cosThetaStar, ptPi, ptKa, cT) from MC particle
	//

	Bool_t isOk = kFALSE;

	// check whether the D0 decays in pi+K
	// to be added!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// could use a cut in the CF, but not clear how to define a TDecayChannel
	// implemented for the time being as a cut on the number of daughters - see later when 
	// getting the daughters

	// getting the daughters
	Int_t daughter0 = mcPart->GetDaughter(0);
	Int_t daughter1 = mcPart->GetDaughter(1);
	AliDebug(2, Form("daughter0 = %d and daughter1 = %d",daughter0,daughter1));
	if (daughter0 == 0 || daughter1 == 0) {
		AliDebug(2, "Error! the D0 MC doesn't have correct daughters!!");
		return isOk;  
	}
	if (TMath::Abs(daughter1 - daughter0 != 1)) {
		AliDebug(2, "The D0 MC doesn't come from a 2-prong decay, skipping!!");
		return isOk;  
	}
	AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter0));
	AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter1));
	if (!mcPartDaughter0 || !mcPartDaughter1) {
		AliWarning("At least one Daughter Particle not found in tree, skipping"); 
		return isOk;  
	}
	
	Double_t vtx1[3] = {0,0,0};   // primary vertex		
	Double_t vtx2daughter0[3] = {0,0,0};   // secondary vertex from daughter 0
	Double_t vtx2daughter1[3] = {0,0,0};   // secondary vertex from daughter 1
	mcPart->XvYvZv(vtx1);  // cm
	// getting vertex from daughters
	mcPartDaughter0->XvYvZv(vtx2daughter0);  // cm
	mcPartDaughter1->XvYvZv(vtx2daughter1);  //cm
	if (vtx2daughter0[0] != vtx2daughter1[0] && vtx2daughter0[1] != vtx2daughter1[1] && vtx2daughter0[2] != vtx2daughter1[2]) {
		AliError("Daughters have different secondary vertex, skipping the track");
		return isOk;
	}
	Int_t nprongs = 2;
	Short_t charge = 0;
	// always instantiate the AliAODRecoDecay with the positive daughter first, the negative second
	AliAODMCParticle* positiveDaugh = mcPartDaughter0;
	AliAODMCParticle* negativeDaugh = mcPartDaughter1;
	if (mcPartDaughter0->GetPdgCode()<0 && mcPartDaughter1->GetPdgCode()>0){
		// inverting in case the positive daughter is the second one
		positiveDaugh = mcPartDaughter1;
		negativeDaugh = mcPartDaughter0;
	}
	// getting the momentum from the daughters
	Double_t px[2] = {positiveDaugh->Px(), negativeDaugh->Px()};		
	Double_t py[2] = {positiveDaugh->Py(), negativeDaugh->Py()};		
	Double_t pz[2] = {positiveDaugh->Pz(), negativeDaugh->Pz()};

	Double_t d0[2] = {0.,0.};		

	AliAODRecoDecayHF* decay = new AliAODRecoDecayHF(vtx1,vtx2daughter0,nprongs,charge,px,py,pz,d0);

	Double_t cosThetaStar = 0.;
	Double_t cosThetaStarD0 = 0.;
	Double_t cosThetaStarD0bar = 0.;
	cosThetaStarD0 = decay->CosThetaStar(1,421,211,321);
	cosThetaStarD0bar = decay->CosThetaStar(0,421,321,211);
	if (mcPart->GetPdgCode() == 421){  // D0
		AliDebug(3, Form("D0, with pdgprong0 = %d, pdgprong1 = %d",mcPartDaughter0->GetPdgCode(),mcPartDaughter1->GetPdgCode()));
		cosThetaStar = cosThetaStarD0;
	}
	else if (mcPart->GetPdgCode() == -421){  // D0bar{
		AliDebug(3, Form("D0bar, with pdgprong0 = %d, pdgprong1 = %d",mcPartDaughter0->GetPdgCode(),mcPartDaughter1->GetPdgCode()));
		cosThetaStar = cosThetaStarD0bar;
	}
	else{
		AliWarning("There are problems!! particle was expected to be either a D0 or a D0bar, check...");
		return vectorMC;
	}
	if (cosThetaStar < -1 || cosThetaStar > 1) {
		AliWarning("Invalid value for cosine Theta star, returning");
		return isOk;
	}

	// calculate cos(Theta*) according to the method implemented herein

	Double_t cts = 9999.;
	cts = CosThetaStar(mcPart, mcPartDaughter0, mcPartDaughter1);
	if (cts < -1 || cts > 1) {
		AliWarning("Invalid value for cosine Theta star from AliCFHeavyFlavourTaskMultiVarMultiStep method");
	}
	if (TMath::Abs(cts - cosThetaStar)>0.001) {
		AliError(Form("cosThetaStar automatically calculated different from that manually calculated!!! cosThetaStar = %f, cosThetaStar = %f", cosThetaStar,cts));
	}
	
	Double_t cT = decay->Ct(421);

	// calculate cT from the method implemented herein
	Double_t cT1 = 0.;
	cT1 = CT(mcPart, mcPartDaughter0, mcPartDaughter1);

	if (TMath::Abs(cT1 - cT)>0.001) {
		AliError(Form("cT automatically calculated different from that manually calculated!!! cT = %f, cT1 = %f",cT,cT1));
	}
	
	// get the pT of the daughters
	
	Double_t pTpi = 0.;
	Double_t pTK = 0.;
	
	if (TMath::Abs(mcPartDaughter0->GetPdgCode()) == 211) {
		pTpi = mcPartDaughter0->Pt();
		pTK = mcPartDaughter1->Pt();
	}
	else {
		pTpi = mcPartDaughter1->Pt();
		pTK = mcPartDaughter0->Pt();
	}

	vectorMC[0] = mcPart->Pt();
	vectorMC[1] = mcPart->Y() ;
	vectorMC[2] = cosThetaStar ;
	vectorMC[3] = pTpi ;
	vectorMC[4] = pTK ;
	vectorMC[5] = cT*1.E4 ;  // in micron
	isOk = kTRUE;
	return isOk;
}

