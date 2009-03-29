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
//
// Example of task running on AliEn
// which provides standard way of calculating acceptance and efficiency
// between different steps of the procedure.
// The ouptut of the task is a AliCFContainer from which the efficiencies
// can be calculated
//-----------------------------------------------------------------------
// Author : C. Zampolli, CERN, on the basis of R. Vernet's example
//-----------------------------------------------------------------------

#include <TCanvas.h>
#include <TParticle.h>
#include <TH1I.h>
#include <TStyle.h>

#include "AliCFHeavyFlavourTaskMultiVar.h"
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
AliCFHeavyFlavourTaskMultiVar::AliCFHeavyFlavourTaskMultiVar() :
	AliAnalysisTaskSE(),
	fPDG(0),
	fCFManager(0x0),
	fHistEventsProcessed(0x0),
	fCountMC(0),
	fEvents(0),
	fFillFromGenerated(kFALSE)
{
	//
	//Default ctor
	//
}
//___________________________________________________________________________
AliCFHeavyFlavourTaskMultiVar::AliCFHeavyFlavourTaskMultiVar(const Char_t* name) :
	AliAnalysisTaskSE(name),
	fPDG(0),
	fCFManager(0x0),
	fHistEventsProcessed(0x0),
	fCountMC(0),
	fEvents(0),
	fFillFromGenerated(kFALSE)
{
	//
	// Constructor. Initialization of Inputs and Outputs
	//
	Info("AliCFHeavyFlavourTaskMultiVar","Calling Constructor");
	/*
	  DefineInput(0) and DefineOutput(0)
	  are taken care of by AliAnalysisTaskSE constructor
	*/
	DefineOutput(1,TH1I::Class());
	DefineOutput(2,AliCFContainer::Class());
}

//___________________________________________________________________________
AliCFHeavyFlavourTaskMultiVar& AliCFHeavyFlavourTaskMultiVar::operator=(const AliCFHeavyFlavourTaskMultiVar& c) 
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
AliCFHeavyFlavourTaskMultiVar::AliCFHeavyFlavourTaskMultiVar(const AliCFHeavyFlavourTaskMultiVar& c) :
	AliAnalysisTaskSE(c),
	fPDG(c.fPDG),
	fCFManager(c.fCFManager),
	fHistEventsProcessed(c.fHistEventsProcessed),
	fCountMC(c.fCountMC),
	fEvents(c.fEvents),
	fFillFromGenerated(c.fFillFromGenerated)
{
	//
	// Copy Constructor
	//
}

//___________________________________________________________________________
AliCFHeavyFlavourTaskMultiVar::~AliCFHeavyFlavourTaskMultiVar() {
	//
	//destructor
	//
	Info("~AliCFHeavyFlavourTaskMultiVar","Calling Destructor");
	if (fCFManager)           delete fCFManager ;
	if (fHistEventsProcessed) delete fHistEventsProcessed ;
}

//_________________________________________________
void AliCFHeavyFlavourTaskMultiVar::UserExec(Option_t *)
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
	Double_t containerInput[6] ;
        
	//loop on the MC event
	
	TClonesArray* mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	if (!mcArray) AliError("Could not find Monte-Carlo in AOD");
	Int_t icountMC = 0;
	
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
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepGenerated);
			icountMC++;
		}
		else {
			AliDebug(3,"Problems in filling the container");
			continue;
		}
 	}    
	
	AliDebug(2, Form("Found %i MC particles that are D0!!",icountMC));
	AliDebug(2, Form("aodEvent = %p ",aodEvent));

	// Now go to rec level
	fCountMC += icountMC;
	// load heavy flavour vertices

	TClonesArray *arrayVerticesHF = (TClonesArray*)((aodEvent->GetList())->FindObject("D0toKpi")); 	
	if (!arrayVerticesHF) AliError("Could not find array of HF vertices");
	AliDebug(2, Form("Found %d vertices",arrayVerticesHF->GetEntriesFast()));
	
	for (Int_t iVertex = 0; iVertex<arrayVerticesHF->GetEntriesFast(); iVertex++) {
		
		AliAODRecoDecayHF2Prong* vtx = (AliAODRecoDecayHF2Prong*)arrayVerticesHF->At(iVertex);
		
		// cuts can't be applied to RecoDecays particles
		// if (!fCFManager->CheckParticleCuts(1  , vtx)) continue;  // 1 stands for AOD level
		
		// find associated MC particle
		Int_t mcLabel = vtx->MatchToMC(421,mcArray) ;
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

			// fill the container
			// either with reconstructed values....
			Double_t pt = vtx->Pt();
			Double_t rapidity = vtx->YD0();
			
			Double_t cosThetaStar = 9999.;
			Double_t pTpi = 0.;
			Double_t pTK = 0.;
			Int_t pdgCode = mcVtxHF->GetPdgCode();
			if (pdgCode > 0){
				cosThetaStar = vtx->CosThetaStarD0();
				pTpi = vtx->PtProng(0);
				pTK = vtx->PtProng(1);
			}
			else {
				cosThetaStar = vtx->CosThetaStarD0bar();
				pTpi = vtx->PtProng(1);
				pTK = vtx->PtProng(0);
			}

			Double_t cT = vtx->CtD0();
			AliDebug(2, Form("cT from reconstructed vertex = %f (micron)",cT*1E4));
			AliDebug(2, Form("pdg code from MC = %d",TMath::Abs(mcVtxHF->GetPdgCode())));

			// ... or with generated values

			if (!fFillFromGenerated){
				containerInput[0] = pt;
				containerInput[1] = rapidity;
				containerInput[2] = cosThetaStar;
				containerInput[3] = pTpi;
				containerInput[4] = pTK;
				containerInput[5] = cT*1.E4;  // in micron
			}
			else {
				Double_t vectorMC[6] = {9999.,9999.,9999.,9999.,9999.,9999.};
				if (GetGeneratedValuesFromMCParticle(mcVtxHF, mcArray, vectorMC)){
					containerInput[0] = vectorMC[0];
					containerInput[1] = vectorMC[1] ;
					containerInput[2] = vectorMC[2] ;
					containerInput[3] = vectorMC[3] ;
					containerInput[4] = vectorMC[4] ;
					containerInput[5] = vectorMC[5] ;  // in micron
				}
				else {
					AliDebug(3,"Problems in filling the container");
					continue;
				}
			}
			AliDebug(2, Form("Filling the container with pt = %f, rapidity = %f, cosThetaStar = %f, pTpi = %f, pTK = %f, cT = %f", containerInput[0], containerInput[1], containerInput[2], containerInput[3], containerInput[4], containerInput[5]));
			fCFManager->GetParticleContainer()->Fill(containerInput,1) ;   
		}
		/*
		// for debugging only, filling container with dummy values
		Double_t pt = 1.5;
		Double_t rapidity = 0.5;
		Double_t cosThetaStar = 0.7;
		Double_t pTpi = 2.5;
		Double_t pTK = 1;
		Double_t cT = 3;
		containerInput[0] = pt;
		containerInput[1] = rapidity;
		containerInput[2] = cosThetaStar;
		containerInput[3] = pTpi;
		containerInput[4] = pTK;
		containerInput[5] = cT;
		fCFManager->GetParticleContainer()->Fill(containerInput,1) ;   
		*/

	}
	
	fHistEventsProcessed->Fill(0);
	/* PostData(0) is taken care of by AliAnalysisTaskSE */
	PostData(1,fHistEventsProcessed) ;
	PostData(2,fCFManager->GetParticleContainer()) ;
}


//___________________________________________________________________________
void AliCFHeavyFlavourTaskMultiVar::Terminate(Option_t*)
{
	// The Terminate() function is the last function to be called during
	// a query. It always runs on the client, it can be used to present
	// the results graphically or save the results to file.
	
	Info("Terminate","");
	AliAnalysisTaskSE::Terminate();
	
	AliInfo(Form("Found %i MC particles that are D0 in MC in %d events",fCountMC,fEvents));
	
	// draw some example plots....
	
	AliCFContainer *cont= dynamic_cast<AliCFContainer*> (GetOutputData(2));
	
	// projecting the containers to obtain histograms
	// first argument = variable, second argument = step
	TH1D* h00 =   cont->ShowProjection(0,0) ;   // pt
	TH1D* h10 =   cont->ShowProjection(1,0) ;   // rapidity
	TH1D* h20 =   cont->ShowProjection(2,0) ;   // cosThetaStar
	TH1D* h30 =   cont->ShowProjection(3,0) ;   // pTpi
	TH1D* h40 =   cont->ShowProjection(4,0) ;   // pTK
	TH1D* h50 =   cont->ShowProjection(5,0) ;   // cT
	
	TH1D* h01 =   cont->ShowProjection(0,1) ;   // pt
	TH1D* h11 =   cont->ShowProjection(1,1) ;   // rapidity
	TH1D* h21 =   cont->ShowProjection(2,1) ;   // cosThetaStar
	TH1D* h31 =   cont->ShowProjection(3,1) ;   // pTpi
	TH1D* h41 =   cont->ShowProjection(4,1) ;   // pTK
	TH1D* h51 =   cont->ShowProjection(5,1) ;   // cT
	
	h00->SetTitle("pT_D0 (GeV/c)");
	h10->SetTitle("rapidity");
	h20->SetTitle("cosThetaStar");
	h30->SetTitle("pT_pi (GeV/c)");
	h40->SetTitle("pT_K (Gev/c)");
	h50->SetTitle("cT (#mum)");

	h00->GetXaxis()->SetTitle("pT_D0 (GeV/c)");
	h10->GetXaxis()->SetTitle("rapidity");
	h20->GetXaxis()->SetTitle("cosThetaStar");
	h30->GetXaxis()->SetTitle("pT_pi (GeV/c)");
	h40->GetXaxis()->SetTitle("pT_K (Gev/c)");
	h50->GetXaxis()->SetTitle("cT (#mum)");

	h01->SetTitle("pT_D0 (GeV/c)");
	h11->SetTitle("rapidity");
	h21->SetTitle("cosThetaStar");
	h31->SetTitle("pT_pi (GeV/c)");
	h41->SetTitle("pT_K (Gev/c)");
	h51->SetTitle("cT (#mum)");

	h01->GetXaxis()->SetTitle("pT_D0 (GeV/c)");
	h11->GetXaxis()->SetTitle("rapidity");
	h21->GetXaxis()->SetTitle("cosThetaStar");
	h31->GetXaxis()->SetTitle("pT_pi (GeV/c)");
	h41->GetXaxis()->SetTitle("pT_K (Gev/c)");
	h51->GetXaxis()->SetTitle("cT (#mum)");

	Double_t max0 = h00->GetMaximum();
	Double_t max1 = h10->GetMaximum();
	Double_t max2 = h20->GetMaximum();
	Double_t max3 = h30->GetMaximum();
	Double_t max4 = h40->GetMaximum();
	Double_t max5 = h50->GetMaximum();
	
	h00->GetYaxis()->SetRangeUser(0,max0*1.2);
	h10->GetYaxis()->SetRangeUser(0,max1*1.2);
	h20->GetYaxis()->SetRangeUser(0,max2*1.2);
	h30->GetYaxis()->SetRangeUser(0,max3*1.2);
	h40->GetYaxis()->SetRangeUser(0,max4*1.2);
	h50->GetYaxis()->SetRangeUser(0,max5*1.2);
	
	h01->GetYaxis()->SetRangeUser(0,max0*1.2);
	h11->GetYaxis()->SetRangeUser(0,max1*1.2);
	h21->GetYaxis()->SetRangeUser(0,max2*1.2);
	h31->GetYaxis()->SetRangeUser(0,max3*1.2);
	h41->GetYaxis()->SetRangeUser(0,max4*1.2);
	h51->GetYaxis()->SetRangeUser(0,max5*1.2);
	
	h00->SetMarkerStyle(20);
	h10->SetMarkerStyle(24);
	h20->SetMarkerStyle(21);
	h30->SetMarkerStyle(25);
	h40->SetMarkerStyle(27);
	h50->SetMarkerStyle(28);

	h00->SetMarkerColor(2);
	h10->SetMarkerColor(2);
	h20->SetMarkerColor(2);
	h30->SetMarkerColor(2);
	h40->SetMarkerColor(2);
	h50->SetMarkerColor(2);

	h01->SetMarkerStyle(20) ;
	h11->SetMarkerStyle(24) ;
	h21->SetMarkerStyle(21) ;
	h31->SetMarkerStyle(25) ;
	h41->SetMarkerStyle(27) ;
	h51->SetMarkerStyle(28) ;

	h01->SetMarkerColor(4);
	h11->SetMarkerColor(4);
	h21->SetMarkerColor(4);
	h31->SetMarkerColor(4);
	h41->SetMarkerColor(4);
	h51->SetMarkerColor(4);
	
	gStyle->SetCanvasColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetStatColor(0);

	// drawing in 2 separate canvas for a matter of clearity
	TCanvas * c1 =new TCanvas("c1","pT, rapidiy, cosThetaStar",800,1600);
	TCanvas * c2 =new TCanvas("c2","pTpi, pTK, cT",800,1600);
	c1->Divide(2,3);
	c2->Divide(2,3);
	
	c1->cd(1);
	h00->Draw("p");
	c1->cd(2);
	h01->Draw("p");
	c1->cd(3);
	h10->Draw("p");
	c1->cd(4);
	h11->Draw("p");
	c1->cd(5);
	h20->Draw("p");
	c1->cd(6);
	h21->Draw("p");
	c1->cd();

	c2->cd(1);
	h30->Draw("p");
	c2->cd(2);
	h31->Draw("p");
	c2->cd(3);
	h40->Draw("p");
	c2->cd(4);
	h41->Draw("p");
	c2->cd(5);
	h50->Draw("p");
	c2->cd(6);
	h51->Draw("p");
	c2->cd();

	c1->SaveAs("Variables/pT_rapidity_cosThetaStar.eps");
	c2->SaveAs("Variables/pTpi_pTK_cT.eps");
	
}

//___________________________________________________________________________
void AliCFHeavyFlavourTaskMultiVar::UserCreateOutputObjects() {
	//HERE ONE CAN CREATE OUTPUT OBJECTS, IN PARTICULAR IF THE OBJECT PARAMETERS DON'T NEED
	//TO BE SET BEFORE THE EXECUTION OF THE TASK
	//
	Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
	
	//slot #1
	OpenFile(1);
	fHistEventsProcessed = new TH1I("fHistEventsProcessed","",1,0,1) ;
}

//___________________________________________________________________________
Double_t AliCFHeavyFlavourTaskMultiVar::CosThetaStar(AliAODMCParticle* mcPart, AliAODMCParticle* mcPartDaughter0, AliAODMCParticle* mcPartDaughter1) const {

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
Double_t AliCFHeavyFlavourTaskMultiVar::CT(AliAODMCParticle* mcPart, AliAODMCParticle* mcPartDaughter0, AliAODMCParticle* mcPartDaughter1) const {

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
Bool_t AliCFHeavyFlavourTaskMultiVar::GetGeneratedValuesFromMCParticle(AliAODMCParticle* mcPart, TClonesArray* mcArray, Double_t* vectorMC) const {

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
		AliWarning("Invalid value for cosine Theta star from AliCFHeavyFlavourTaskMultiVar method");
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

