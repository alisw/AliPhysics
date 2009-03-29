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
// Example of task
// which provides standard way of calculating acceptance and efficiency
// between different steps of the procedure.
// The ouptut of the task is a AliCFContainer from which the efficiencies
// can be calculated
// Applied on AliAODRecoDecayHF2Prong particles (D0->Kpi)
//-----------------------------------------------------------------------
// Author : C. Zampolli starting from an example from R. Vernet
//-----------------------------------------------------------------------

#include <TH1I.h>
#include <TCanvas.h>
#include <TStyle.h>

#include "AliCFHeavyFlavourTask.h"
#include "AliCFContainer.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODMCParticle.h"
#include "AliCFManager.h"

//__________________________________________________________________________
AliCFHeavyFlavourTask::AliCFHeavyFlavourTask() :
	AliAnalysisTaskSE(),
	fPDG(0),
	fCFManager(0x0),
	fHistEventsProcessed(0x0),
	fCountMC(0),
	fEvents(0)
{
	//
	//Default ctor
	//
}
//___________________________________________________________________________
AliCFHeavyFlavourTask::AliCFHeavyFlavourTask(const Char_t* name) :
	AliAnalysisTaskSE(name),
	fPDG(0),
	fCFManager(0x0),
	fHistEventsProcessed(0x0),
	fCountMC(0),
	fEvents(0)
{
	//
	// Constructor. Initialization of Inputs and Outputs
	//
	/*
	  DefineInput(0) and DefineOutput(0)
	  are taken care of by AliAnalysisTaskSE constructor
	*/
	DefineOutput(1,TH1I::Class());
	DefineOutput(2,AliCFContainer::Class());
}

//___________________________________________________________________________
AliCFHeavyFlavourTask& AliCFHeavyFlavourTask::operator=(const AliCFHeavyFlavourTask& c) 
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
AliCFHeavyFlavourTask::AliCFHeavyFlavourTask(const AliCFHeavyFlavourTask& c) :
	AliAnalysisTaskSE(c),
	fPDG(c.fPDG),
	fCFManager(c.fCFManager),
	fHistEventsProcessed(c.fHistEventsProcessed),
	fCountMC(c.fCountMC),
	fEvents(c.fEvents)
{
	//
	// Copy Constructor
	//
}

//___________________________________________________________________________
AliCFHeavyFlavourTask::~AliCFHeavyFlavourTask() {
	//
	//destructor
	//
	if (fCFManager)           delete fCFManager ;
	if (fHistEventsProcessed) delete fHistEventsProcessed ;
}

//_________________________________________________
void AliCFHeavyFlavourTask::UserExec(Option_t *)
{
	//
	// Main loop function
	//
	
	if (!fInputEvent) {
		Error("UserExec","NO EVENT FOUND!");
		return;
	}
	
	fEvents++;
	if (fEvents%10000 ==0) AliInfo(Form("Event %d",fEvents));
	AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
	fCFManager->SetEventInfo(aodEvent);
	
	// MC-event selection
	Double_t containerInput[2] ;
        
	// loop on the MC event
	TClonesArray* mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	if (!mcArray) AliError("Could not find Monte-Carlo in AOD");
	Int_t icountMC = 0;
	
	for (Int_t iPart=0; iPart<mcArray->GetEntriesFast(); iPart++) { 
		AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(iPart));
		if (!mcPart) AliWarning("Particle not found in tree");
		
		// check the MC-level cuts
		if (!fCFManager->CheckParticleCuts(0,mcPart)) continue; // 0 stands for MC level

		// check whether the D0 decays in pi+K
		// to be added!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// could use a cut in the CF, but not clear how to define a TDecayChannel
		// implemented for the time being as a cut on the number of daughters - see later when 
		// getting the daughters
		
		// getting the daughters: 
		// AliMCParticle::GetDaughter(0) returns the index of the 1st daughter
		// AliMCParticle::GetDaughter(1) returns the index of the last daughter
		// so when a particle has 2 daughters, the difference must be 1
		Int_t daughter0 = mcPart->GetDaughter(0);
		Int_t daughter1 = mcPart->GetDaughter(1);
		if (daughter0 == 0 || daughter1 == 0) {
			AliDebug(2, "Error! the D0 MC doesn't have correct daughters!!");
			continue;  
		}
		if (TMath::Abs(daughter1 - daughter0 != 1)) {
			AliDebug(2, "The D0 MC doesn't come from a 2-prong decay, skipping!!");
			continue;  
		}
		AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter0));
		AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter1));
		if (!mcPartDaughter0 || !mcPartDaughter1) {
			AliWarning("At least one Daughter Particle not found in tree, skipping"); 
			continue;  
		}

		// fill the container for Gen-level selection
		containerInput[0] = mcPart->Pt();
		containerInput[1] = mcPart->Y();
		fCFManager->GetParticleContainer()->Fill(containerInput,kStepGenerated);
		icountMC++;
	}    
	
	AliDebug(2,Form("Found %i MC particles that are D0!!",icountMC));
	// Now go to rec level
	fCountMC += icountMC;
	// load D0->Kpi candidates
	TClonesArray *arrayD0toKpi = (TClonesArray*)aodEvent->GetList()->FindObject("D0toKpi"); 
	
	if (!arrayD0toKpi) AliError("Could not find array of D0->Kpi candidates");
	
	AliDebug(2, Form("Found %d vertices",arrayD0toKpi->GetEntriesFast()));
	
	for (Int_t iD0 = 0; iD0<arrayD0toKpi->GetEntriesFast(); iD0++) {
		
		AliAODRecoDecayHF2Prong* rd = (AliAODRecoDecayHF2Prong*)arrayD0toKpi->At(iD0);
		
		// cuts can't be applied to RecoDecays particles
		// if (!fCFManager->CheckParticleCuts(1, rd)) continue;  // 1 stands for AOD level
		
		// check if associated MC v0 passes the cuts
		// first get the label
		Int_t mcLabel = rd->MatchToMC(421,mcArray);
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
					
			AliDebug(2, Form("pdg code from MC = %d",TMath::Abs(mcVtxHF->GetPdgCode())));
			if (!fCFManager->CheckParticleCuts(0, mcVtxHF)) {  // 0 stands for MC level
				continue; 
			}
			
			//fill the container
			Double_t pt = rd->Pt();
			Double_t rapidity = rd->YD0();
			
			AliDebug(2, Form("Filling the container with pt = %f and rapidity = %f", pt, rapidity));
			containerInput[0] = pt ;
			containerInput[1] = rapidity ;
			fCFManager->GetParticleContainer()->Fill(containerInput,1) ;   
		}
	}
	
	fHistEventsProcessed->Fill(0);
	// PostData(0) is taken care of by AliAnalysisTaskSE 
	PostData(1,fHistEventsProcessed) ;
	PostData(2,fCFManager->GetParticleContainer()) ;
}


//___________________________________________________________________________
void AliCFHeavyFlavourTask::Terminate(Option_t*)
{
	// The Terminate() function is the last function to be called during
	// a query. It always runs on the client, it can be used to present
	// the results graphically or save the results to file.
	
	AliAnalysisTaskSE::Terminate();
	
	AliInfo(Form("Found %i MC particles that are D0 in MC in %d events",fCountMC,fEvents));
	
	//draw some example plots....
	
	AliCFContainer *cont= dynamic_cast<AliCFContainer*> (GetOutputData(2));
	
	// projecting the containers to obtain histograms
	// first argument = variable, second argument = step
	TH1D* h00 =   cont->ShowProjection(0,0) ;  // pt, MC
	TH1D* h01 =   cont->ShowProjection(0,1) ;  // pt, Reco
	
	TH1D* h10 =   cont->ShowProjection(1,0) ;  // rapidity, MC
	TH1D* h11 =   cont->ShowProjection(1,1) ;  // rapidity, Reco
	
	Double_t max1 = h00->GetMaximum();
	Double_t max2 = h10->GetMaximum();
	
	// MC histos
	h00->SetTitle("p_{T} (GeV/c)");
	h10->SetTitle("rapidity");
	h00->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h10->GetXaxis()->SetTitle("rapidity");
	h00->GetYaxis()->SetRangeUser(0,max1*1.2);
	h10->GetYaxis()->SetRangeUser(0,max2*1.2);
	h00->SetMarkerStyle(20) ;
	h10->SetMarkerStyle(24) ;
	h00->SetMarkerColor(2);
	h10->SetMarkerColor(2);

	// Reco histos
	h01->SetTitle("p_{T} (GeV/c)");
	h11->SetTitle("rapidity");
	h01->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	h11->GetXaxis()->SetTitle("rapidity");
	h01->GetYaxis()->SetRangeUser(0,max1*1.2);
	h11->GetYaxis()->SetRangeUser(0,max2*1.2);
	h01->SetMarkerStyle(20) ;
	h11->SetMarkerStyle(24) ;
	h01->SetMarkerColor(4);
	h11->SetMarkerColor(4);

	gStyle->SetCanvasColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetStatColor(0);

	TCanvas * c =new TCanvas("c","",1000,600);
	c->Divide(2,2);
	
	c->cd(1);
	h00->Draw("p");
	c->cd(2);
	h01->Draw("p");
	c->cd(3);
	h10->Draw("p");
	c->cd(4);
	h11->Draw("p");
	c->cd();

	c->SaveAs("plots.eps");
}

//___________________________________________________________________________
void AliCFHeavyFlavourTask::UserCreateOutputObjects() {
	//HERE ONE CAN CREATE OUTPUT OBJECTS, IN PARTICULAR IF THE OBJECT PARAMETERS DON'T NEED
	//TO BE SET BEFORE THE EXECUTION OF THE TASK
	//
	Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
	
	//slot #1
	OpenFile(1);
	fHistEventsProcessed = new TH1I("fHistEventsProcessed","",1,0,1) ;
}

