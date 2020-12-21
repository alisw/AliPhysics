/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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
 * provided "as is" without express or implied warranty.                  * **************************************************************************/


#include "AliAnalysisTaskQAHighPtDeDx.h"

// ROOT includes
#include <TList.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TProfile.h>
#include <TParticle.h>
#include <TFile.h>

// AliRoot includes
#include <AliAnalysisManager.h>
#include <AliAnalysisFilter.h>
#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliESDVertex.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrackCuts.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>

#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>
#include "AliMultSelection.h"

#include <TTreeStream.h>

#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>

#include "AliCentrality.h" 
#include <AliESDv0.h>
#include <AliKFVertex.h>
#include <AliAODVertex.h>

#include <AliAODTrack.h> 
#include <AliAODPid.h> 
#include <AliAODMCHeader.h> 

// STL includes
#include <iostream>
using namespace std;


// Responsible:
// Antonio Ortiz (UNAM), antonio.ortiz@nucleares.unam.mx
// Peter Christiansen (Lund University)


const Double_t AliAnalysisTaskQAHighPtDeDx::fgkClight = 2.99792458e-2;
namespace {
	Float_t magf = -1;
	TF1 cutLow ("StandardPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 50);
	TF1 cutHigh("StandardPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 50);
	Double_t DeDxMIPMin  = 30;
	Double_t DeDxMIPMax  = 65;
	const Int_t nHists = 9;
	Float_t centralityGlobal = -10;
	Int_t etaLow[nHists]  = {-8, -8, -6, -4, -2, 0, 2, 4, 6};
	Int_t etaHigh[nHists] = { 8, -6, -4, -2,  0, 2, 4, 6, 8};

	Int_t nDeltaPiBins   = 80;
	Double_t deltaPiLow  = 20;
	Double_t deltaPiHigh = 100;
	const Char_t *Pid[7]={"Ch","Pion","Kaon","Proton","Electron","Muon","Oher"};
}

ClassImp(AliAnalysisTaskQAHighPtDeDx)

//_______________________________________________________
AliAnalysisTaskQAHighPtDeDx::AliAnalysisTaskQAHighPtDeDx():
		AliAnalysisTaskSE(),
		fESD(0x0),
		fAOD(0x0),
		fMC(0x0),
		fMCStack(0x0),
		fMCArray(0x0),
		fTrackFilterGolden(0x0),
		fTrackFilterTPC(0x0),
		fCentEst("V0M"),
		fAnalysisType("ESD"),
		fAnalysisMC(kFALSE),
		fAnalysisPbPb(kFALSE),
		ftrigBit(0x0),
		fRandom(0x0),
		fPileUpRej(kFALSE),
		fVtxCut(10.0),  
		fEtaCut(0.9),  
		fMinCent(0.0),
		fMaxCent(100.0),
		fStoreMcIn(kFALSE),//
		fMcProcessType(-999),
		fTriggeredEventMB(-999),
		fVtxStatus(-999),
		fZvtx(-999),
		fZvtxMC(-999),
		fRun(-999),
		fEventId(-999),
		fListOfObjects(0), 
		fEvents(0x0), fVtx(0x0), fVtxMC(0x0), fVtxBeforeCuts(0x0), fVtxAfterCuts(0x0),
		fn1(0x0),
		fcent(0x0),
		hMIPVsEta(0x0),
		pMIPVsEta(0x0),
		hMIPVsEtaV0s(0x0),
		pMIPVsEtaV0s(0x0),
		hPlateauVsEta(0x0),
		pPlateauVsEta(0x0),
		hPhi(0x0)


{
	//default constructor
	for(Int_t i=0;i<9;++i){

		hMIPVsNch[i]=0;//TH2D, MIP vs Nch for different eta intervals
		pMIPVsNch[i]=0;//TProfile, MIP vs Nch for different eta intervals
		hMIPVsPhi[i]=0;//TH2D, MIP vs phi for different eta intervals
		pMIPVsPhi[i]=0;//TProfile, MIP vs phi for different eta intervals
		hPlateauVsPhi[i]=0;//TH2D, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c
		pPlateauVsPhi[i]=0;//TProfile, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c
		histPiV0[i]=0;//TH2D, dE/dx vs p, pi id by V0s
		histpPiV0[i]=0;//TH1D, pi id by V0s
		histPV0[i]=0;// TH2D, dE/dx vs p, p id by V0s
		histpPV0[i]=0;// TH1D, p id by V0s
		histAllCh[i]=0;//TH2D, dE/dx vs p for all charged particles
		histPiTof[i]=0;//TH2D, dE/dx vs p for a "clean" sample of pions, beta>1
		histpPiTof[i]=0;//TH1D, for a "clean" sample of pions, beta>1
		histEV0[i]=0;

		for(Int_t pid=0;pid<7;++pid){
			hMcIn[pid][i]=0;
			hMcOut[pid][i]=0;
		}

	}



}


AliAnalysisTaskQAHighPtDeDx::AliAnalysisTaskQAHighPtDeDx(const char *name):
	AliAnalysisTaskSE(name),
	fESD(0x0),
	fAOD(0x0),
	fMC(0x0),
	fMCStack(0x0),
	fMCArray(0x0),
	fTrackFilterGolden(0x0),
	fTrackFilterTPC(0x0),
	fCentEst("V0M"),
	fAnalysisType("ESD"),
	fAnalysisMC(kFALSE),
	fAnalysisPbPb(kFALSE),
	ftrigBit(0x0),
	fRandom(0x0),
	fPileUpRej(kFALSE),
	fVtxCut(10.0),  
	fEtaCut(0.9),  
	fMinCent(0.0),
	fMaxCent(100.0),
	fStoreMcIn(kFALSE),//
	fMcProcessType(-999),
	fTriggeredEventMB(-999),
	fVtxStatus(-999),
	fZvtx(-999),
	fZvtxMC(-999),
	fRun(-999),
	fEventId(-999),
	fListOfObjects(0), 
	fEvents(0x0), fVtx(0x0), fVtxMC(0x0), fVtxBeforeCuts(0x0), fVtxAfterCuts(0x0),
	fn1(0x0),
	fcent(0x0),
	hMIPVsEta(0x0),
	pMIPVsEta(0x0),
	hMIPVsEtaV0s(0x0),
	pMIPVsEtaV0s(0x0),
	hPlateauVsEta(0x0),
	pPlateauVsEta(0x0),
	hPhi(0x0)


{
	// Default constructor (should not be used)
	for(Int_t i=0;i<9;++i){

		hMIPVsNch[i]=0;//TH2D, MIP vs Nch for different eta intervals
		pMIPVsNch[i]=0;//TProfile, MIP vs Nch for different eta intervals
		hMIPVsPhi[i]=0;//TH2D, MIP vs phi for different eta intervals
		pMIPVsPhi[i]=0;//TProfile, MIP vs phi for different eta intervals
		hPlateauVsPhi[i]=0;//TH2D, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c
		pPlateauVsPhi[i]=0;//TProfile, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c
		histPiV0[i]=0;//TH2D, dE/dx vs p, pi id by V0s
		histpPiV0[i]=0;//TH1D, pi id by V0s
		histPV0[i]=0;// TH2D, dE/dx vs p, p id by V0s
		histpPV0[i]=0;// TH1D, p id by V0s
		histAllCh[i]=0;//TH2D, dE/dx vs p for all charged particles
		histPiTof[i]=0;//TH2D, dE/dx vs p for a "clean" sample of pions, beta>1
		histpPiTof[i]=0;//TH1D, for a "clean" sample of pions, beta>1
		histEV0[i]=0;

		for(Int_t pid=0;pid<7;++pid){
			hMcIn[pid][i]=0;
			hMcOut[pid][i]=0;
		}

	}
	DefineOutput(1, TList::Class());//esto es nuevo
}




AliAnalysisTaskQAHighPtDeDx::~AliAnalysisTaskQAHighPtDeDx() {
	//
	// Destructor
	//

}

//______________________________________________________________________________
void AliAnalysisTaskQAHighPtDeDx::UserCreateOutputObjects()
{ 
	// This method is called once per worker node
	// Here we define the output: histograms and debug tree if requested 
	// We also create the random generator here so it might get different seeds...
	fRandom = new TRandom(0); // 0 means random seed



	//OpenFile(1);
	fListOfObjects = new TList();
	fListOfObjects->SetOwner();

	//
	// Histograms
	//  
	fEvents = new TH1I("fEvents","Number of analyzed events; Events; Counts", 3, 0, 3);
	fListOfObjects->Add(fEvents);

	fn1=new TH1F("fn1","fn1",11,-1,10);
	fListOfObjects->Add(fn1);

	fcent=new TH1F("fcent","fcent",104,-2,102);
	fListOfObjects->Add(fcent);

	fVtx = new TH1I("fVtx","Vtx info (0=no, 1=yes); Vtx; Counts", 2, -0.5, 1.5);
	fListOfObjects->Add(fVtx);

	fVtxBeforeCuts = new TH1F("fVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 120, -30, 30);
	fListOfObjects->Add(fVtxBeforeCuts);

	fVtxAfterCuts = new TH1F("fVtxAfterCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 120, -30, 30);
	fListOfObjects->Add(fVtxAfterCuts);


	const Int_t nPtBinsV0s = 25;
	Double_t ptBinsV0s[nPtBinsV0s+1] = { 0.0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 ,
		1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.5 , 3.0 , 3.5 , 4.0 , 5.0 , 7.0 ,
		9.0 , 12.0, 15.0, 20.0 };




	const Char_t* ending[nHists] = {"", "86", "64", "42", "20", "02", "24", "46", "68"};

	const Char_t* LatexEta[nHists] = {
		"|#eta|<0.8", "-0.8<#eta<-0.6", "-0.6<#eta<-0.4", "-0.4<#eta<-0.2", "-0.2<#eta<0", 
		"0<#eta<0.2", "0.2<#eta<0.4", "0.4<#eta<0.6", "0.6<#eta<0.8" 
	};
	hPhi = new TH2D("histPhi", "pt; #phi'", nPtBinsV0s, ptBinsV0s, 90, -0.05, 0.4); 
	//dE/dx vs phi, pions at the MIP
	fListOfObjects->Add(hPhi);




	Int_t nPhiBins = 36;

	for(Int_t i = 0; i < nHists; i++) {


		hMIPVsPhi[i]      = new TH2D(Form("hMIPVsPhi%s", ending[i]), Form("%s; #phi (rad); dE/dx MIP",LatexEta[i]), nPhiBins, 0, 2*TMath::Pi(), 
				DeDxMIPMax-DeDxMIPMin, DeDxMIPMin, DeDxMIPMax);
		hMIPVsPhi[i]->Sumw2();

		pMIPVsPhi[i]     = new TProfile(Form("pMIPVsPhi%s", ending[i]), Form("%s; #phi (rad); dE/dx MIP",LatexEta[i]),  nPhiBins, 0, 2*TMath::Pi(),
				DeDxMIPMin, DeDxMIPMax);
		pMIPVsPhi[i]->SetMarkerStyle(20);
		pMIPVsPhi[i]->SetMarkerColor(1);
		pMIPVsPhi[i]->SetLineColor(1);
		pMIPVsPhi[i]->Sumw2();

		fListOfObjects->Add(hMIPVsPhi[i]);
		fListOfObjects->Add(pMIPVsPhi[i]);

		hPlateauVsPhi[i]  = new TH2D(Form("hPlateauVsPhi%s", ending[i]), Form("%s; #phi (rad); dE/dx Plateau",LatexEta[i]),  nPhiBins, 0, 2*TMath::Pi(),
				95-DeDxMIPMax, DeDxMIPMax, 95);
		hPlateauVsPhi[i]->Sumw2();

		pPlateauVsPhi[i] = new TProfile(Form("pPlateauVsPhi%s", ending[i]), Form("%s; #phi (rad); dE/dx Plateau",LatexEta[i]), nPhiBins, 0, 2*TMath::Pi(),
				DeDxMIPMax, 95);
		pPlateauVsPhi[i]->SetMarkerStyle(20);
		pPlateauVsPhi[i]->SetMarkerColor(1);
		pPlateauVsPhi[i]->SetLineColor(1);
		pPlateauVsPhi[i]->Sumw2();

		fListOfObjects->Add(hPlateauVsPhi[i]);
		fListOfObjects->Add(pPlateauVsPhi[i]);


		hMIPVsNch[i]      = new TH2D(Form("hMIPVsNch%s", ending[i]), Form("%s; TPC track mult. |#eta|<0.8; dE/dx MIP",LatexEta[i]), 400, 1, 2001, 
				DeDxMIPMax-DeDxMIPMin, DeDxMIPMin, DeDxMIPMax);
		hMIPVsNch[i]->Sumw2();

		pMIPVsNch[i] = new TProfile(Form("pMIPVsNch%s", ending[i]), Form("%s; TPC track mult. |#eta|<0.8; dE/dx MIP",LatexEta[i]), 400, 1, 2001, DeDxMIPMin, DeDxMIPMax);
		pMIPVsNch[i]->SetMarkerStyle(20);
		pMIPVsNch[i]->SetMarkerColor(1);
		pMIPVsNch[i]->SetLineColor(1);
		pMIPVsNch[i]->Sumw2();

		fListOfObjects->Add(hMIPVsNch[i]);
		fListOfObjects->Add(pMIPVsNch[i]);

		//two dimmesional distributions dE/dx vs p for secondary pions
		histPiV0[i]  = new TH2D(Form("histPiV0%s", ending[i]), "Pions id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
		histPiV0[i]->Sumw2();
		fListOfObjects->Add(histPiV0[i]);

		histpPiV0[i]  = new TH1D(Form("histPiV01D%s", ending[i]), "Pions id by V0; #it{p} (GeV/#it{c}); counts", 200, 0, 20);
		histpPiV0[i]->Sumw2();
		fListOfObjects->Add(histpPiV0[i]);

		//two dimmesional distributions dE/dx vs p for secondary protons
		histPV0[i]   = new TH2D(Form("histPV0%s", ending[i]), "Protons id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
		histPV0[i]->Sumw2();
		fListOfObjects->Add(histPV0[i]);

		histpPV0[i]  = new TH1D(Form("histPV01D%s", ending[i]), "Protons id by V0; #it{p} (GeV/#it{c}); counts", 200, 0, 20);
		histpPV0[i]->Sumw2();
		fListOfObjects->Add(histpPV0[i]);

		//two dimmesional distributions dE/dx vs p for primary pions
		histPiTof[i] = new TH2D(Form("histPiTof%s", ending[i]), "all charged", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
		histPiTof[i]->Sumw2();

		histpPiTof[i]  = new TH1D(Form("histPiTof1D%s", ending[i]), "Protons id by V0; #it{p} (GeV/#it{c}); counts", 200, 0, 20);
		histpPiTof[i]->Sumw2();
		fListOfObjects->Add(histpPiTof[i]);


		histAllCh[i] = new TH2D(Form("histAllCh%s", ending[i]), "Pions id by TOF", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
		histAllCh[i]->Sumw2();

		fListOfObjects->Add(histPiTof[i]);
		fListOfObjects->Add(histAllCh[i]);

		histEV0[i]   = new TH2D(Form("histEV0%s", ending[i]), "Electrons id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
		histEV0[i]->Sumw2();
		fListOfObjects->Add(histEV0[i]);



	}


	hMIPVsEta = new TH2D("hMIPVsEta","; #eta; dE/dx_{MIP, primary tracks}",16,-0.8,0.8,DeDxMIPMax-DeDxMIPMin, DeDxMIPMin, DeDxMIPMax);
	pMIPVsEta = new TProfile("pMIPVsEta","; #eta; #LT dE/dx #GT_{MIP, primary tracks}",16,-0.8,0.8, DeDxMIPMin, DeDxMIPMax);
	hMIPVsEtaV0s = new TH2D("hMIPVsEtaV0s","; #eta; dE/dx_{MIP, secondary tracks}",16,-0.8,0.8,DeDxMIPMax-DeDxMIPMin, DeDxMIPMin, DeDxMIPMax);
	pMIPVsEtaV0s = new TProfile("pMIPVsEtaV0s","; #eta; #LT dE/dx #GT_{MIP, secondary tracks}",16,-0.8,0.8,DeDxMIPMin, DeDxMIPMax);

	hPlateauVsEta = new TH2D("hPlateauVsEta","; #eta; dE/dx_{Plateau, primary tracks}",16,-0.8,0.8,95-DeDxMIPMax, DeDxMIPMax, 95);
	pPlateauVsEta = new TProfile("pPlateauVsEta","; #eta; #LT dE/dx #GT_{Plateau, primary tracks}",16,-0.8,0.8, DeDxMIPMax, 95);

	fListOfObjects->Add(hMIPVsEta);
	fListOfObjects->Add(pMIPVsEta);
	fListOfObjects->Add(hMIPVsEtaV0s);
	fListOfObjects->Add(pMIPVsEtaV0s);
	fListOfObjects->Add(hPlateauVsEta);
	fListOfObjects->Add(pPlateauVsEta);


	if (fAnalysisMC) {    
		for(Int_t i = 0; i < nHists; i++) {
			for(Int_t pid = 0; pid < 7; pid++) {

				hMcIn[pid][i] = new TH1D(Form("hIn%s%s", Pid[pid],ending[i]), Form("MC in (pid %s)", Pid[pid]), 
						nPtBinsV0s, ptBinsV0s);
				fListOfObjects->Add(hMcIn[pid][i]);

				hMcOut[pid][i] = new TH1D(Form("hMcOut%s%s", Pid[pid],ending[i]), Form("MC out (pid %s)", Pid[pid]), 
						nPtBinsV0s, ptBinsV0s);
				fListOfObjects->Add(hMcOut[pid][i]);


			}
		}

		fVtxMC = new TH1F("fVtxMC","mc vtx", 120, -30, 30);
		fListOfObjects->Add(fVtxMC);

	}

	// Post output data.
	PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskQAHighPtDeDx::UserExec(Option_t *) 
{
	// Main loop

	//
	// First we make sure that we have valid input(s)!
	//


	AliVEvent *event = InputEvent();
	if (!event) {
		Error("UserExec", "Could not retrieve event");
		return;
	}



	if (fAnalysisType == "ESD"){
		fESD = dynamic_cast<AliESDEvent*>(event);
		if(!fESD){
			Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}    
	} else {
		fAOD = dynamic_cast<AliAODEvent*>(event);
		if(!fAOD){
			Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}    
	}


	if (fAnalysisMC) {

		if (fAnalysisType == "ESD"){
			fMC = dynamic_cast<AliMCEvent*>(MCEvent());
			if(!fMC){
				Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
				this->Dump();
				return;
			}    

			fMCStack = fMC->Stack();

			if(!fMCStack){
				Printf("%s:%d MCStack not found in Input Manager",(char*)__FILE__,__LINE__);
				this->Dump();
				return;
			}    
		} else { // AOD

			fMC = dynamic_cast<AliMCEvent*>(MCEvent());
			if(fMC)
				fMC->Dump();

			fMCArray = (TClonesArray*)fAOD->FindListObject("mcparticles");
			if(!fMCArray){
				Printf("%s:%d AOD MC array not found in Input Manager",(char*)__FILE__,__LINE__);
				this->Dump();
				return;
			}    
		}
	}


	// Get trigger decision
	fTriggeredEventMB = 0; //init

	fn1->Fill(0);

	if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
			->IsEventSelected() & ftrigBit ){
		fTriggeredEventMB = 1;  //event triggered as minimum bias
	}

	// Get process type for MC
	fMcProcessType = 0; // -1=invalid, 0=data, 1=ND, 2=SD, 3=DD

	// real data that are not triggered we skip
	if(!fAnalysisMC && !fTriggeredEventMB)
		return; 

	fn1->Fill(1);


	if (fAnalysisMC) {



		if (fAnalysisType == "ESD"){



			AliHeader* headerMC = fMC->Header();
			if (headerMC) {

				AliGenEventHeader* genHeader = headerMC->GenEventHeader();
				TArrayF vtxMC(3); // primary vertex  MC 
				vtxMC[0]=9999; vtxMC[1]=9999;  vtxMC[2]=9999; //initialize with dummy
				if (genHeader) {
					genHeader->PrimaryVertex(vtxMC);
				}
				fZvtxMC = vtxMC[2];

				// PYTHIA:
				AliGenPythiaEventHeader* pythiaGenHeader =
					dynamic_cast<AliGenPythiaEventHeader*>(headerMC->GenEventHeader());
				if (pythiaGenHeader) {  //works only for pythia
					fMcProcessType =  GetPythiaEventProcessType(pythiaGenHeader->ProcessType());
				}
				// PHOJET:
				AliGenDPMjetEventHeader* dpmJetGenHeader =
					dynamic_cast<AliGenDPMjetEventHeader*>(headerMC->GenEventHeader());
				if (dpmJetGenHeader) {
					fMcProcessType = GetDPMjetEventProcessType(dpmJetGenHeader->ProcessType());
				}
			}
		} else { // AOD



			AliAODMCHeader* mcHeader = dynamic_cast<AliAODMCHeader*>(fAOD->FindListObject("mcHeader")); 


			if(mcHeader) {
				fZvtxMC = mcHeader->GetVtxZ();



				if(strstr(mcHeader->GetGeneratorName(), "Pythia")) {
					fMcProcessType =  GetPythiaEventProcessType(mcHeader->GetEventType());
				} else {
					fMcProcessType =  GetDPMjetEventProcessType(mcHeader->GetEventType());
				}
			}
		}


	}



	if (fAnalysisType == "ESD"){

		const AliESDVertex *vtxESD = fESD->GetPrimaryVertexTracks();
		if(vtxESD->GetNContributors()<1) {
			// SPD vertex
			vtxESD = fESD->GetPrimaryVertexSPD();
			/* quality checks on SPD-vertex */
			if (vtxESD->IsFromVertexerZ() && (vtxESD->GetDispersion() > 0.04 || vtxESD->GetZRes() > 0.25))  
				fZvtx  = -1599; //vertex = 0x0; //
			else if (vtxESD->GetNContributors()<1) 
				fZvtx  = -999; //vertex = 0x0; //
			else
				fZvtx = vtxESD->GetZ();
		}  
		else
			fZvtx = vtxESD->GetZ();

	}
	else // AOD
		fZvtx = GetVertex(fAOD);

	fVtxBeforeCuts->Fill(fZvtx);

	//cut on the z position of vertex
	if (TMath::Abs(fZvtx) > fVtxCut) {	
		return;
	}
	fn1->Fill(2);

	Float_t centrality = -10;


	// only analyze triggered events
	if(fTriggeredEventMB) {

		if (fAnalysisType == "ESD"){
			AliMultSelection *MultSelection = (AliMultSelection*) fESD -> FindListObject("MultSelection");
			if(fCentEst == "V0M")
				centrality = MultSelection->GetMultiplicityPercentile("V0M");
			if(fCentEst == "V0A")
				centrality = MultSelection->GetMultiplicityPercentile("V0A");
			centralityGlobal = centrality;
			if((centrality>fMaxCent)||(centrality<fMinCent))return;
			fcent->Fill(centrality);
			fn1->Fill(3);
			if(fAnalysisMC){
				if(TMath::Abs(fZvtxMC)<fVtxCut){
					ProcessMCTruthESD();
					fVtxMC->Fill(fZvtxMC);
				}
			}
			AnalyzeESD(fESD);
		} else { // AOD
			AliMultSelection *MultSelection = (AliMultSelection*) fAOD -> FindListObject("MultSelection");
			if(fCentEst == "V0M")
				centrality = MultSelection->GetMultiplicityPercentile("V0M");
			if(fCentEst == "V0A")
				centrality = MultSelection->GetMultiplicityPercentile("V0A");
			if((centrality>fMaxCent)||(centrality<fMinCent))return;
			fcent->Fill(centrality);
			fn1->Fill(3);
			if(fAnalysisMC){
				if(TMath::Abs(fZvtxMC)<fVtxCut){

					ProcessMCTruthAOD();
					fVtxMC->Fill(fZvtxMC);
				}
			}
			AnalyzeAOD(fAOD);
		}
	}

	fVtxAfterCuts->Fill(fZvtx);




	// Post output data.
	PostData(1, fListOfObjects);
}

//________________________________________________________________________
void AliAnalysisTaskQAHighPtDeDx::AnalyzeESD(AliESDEvent* esdEvent)
{
	fRun  = esdEvent->GetRunNumber();
	fEventId = 0;
	if(esdEvent->GetHeader())
		fEventId = GetEventIdAsLong(esdEvent->GetHeader());



	Bool_t isPileup = esdEvent->IsPileupFromSPD();
	if(fPileUpRej)
		if(isPileup)
			return;
	fn1->Fill(4);

	magf      = esdEvent->GetMagneticField();





	if(fTriggeredEventMB) {// Only MC case can we have not triggered events

		// accepted event
		fEvents->Fill(0);


		ProduceArrayTrksESD( esdEvent );//produce array with global track parameters
		ProduceArrayV0ESD( esdEvent );//v0's


		fEvents->Fill(1);




	} // end if triggered


}

//________________________________________________________________________
void AliAnalysisTaskQAHighPtDeDx::AnalyzeAOD(AliAODEvent* aodEvent)
{
	fRun  = aodEvent->GetRunNumber();
	fEventId = 0;
	if(aodEvent->GetHeader())
		fEventId = GetEventIdAsLong(aodEvent->GetHeader());

	magf      = aodEvent->GetMagneticField();


	Bool_t isPileup = aodEvent->IsPileupFromSPD();
	if(fPileUpRej)
		if(isPileup)
			return;
	fn1->Fill(4);



	if(fTriggeredEventMB) {// Only MC case can we have not triggered events

		// accepted event
		fEvents->Fill(0);


		ProduceArrayTrksAOD( aodEvent );
		ProduceArrayV0AOD( aodEvent );

		fEvents->Fill(1);




	} // end if triggered

}

//_____________________________________________________________________________
Float_t AliAnalysisTaskQAHighPtDeDx::GetVertex(const AliVEvent* event) const
{
	Float_t zvtx = -999;

	const AliVVertex* primaryVertex = event->GetPrimaryVertex(); 

	if(primaryVertex->GetNContributors()>0)
		zvtx = primaryVertex->GetZ();

	return zvtx;
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskQAHighPtDeDx::GetPidCode(Int_t pdgCode) const 
{
	// return our internal code for pions, kaons, and protons

	Short_t pidCode = 6;

	switch (TMath::Abs(pdgCode)) {
		case 211:
			pidCode = 1; // pion
			break;
		case 321:
			pidCode = 2; // kaon
			break;
		case 2212:
			pidCode = 3; // proton
			break;
		case 11:
			pidCode = 4; // electron
			break;
		case 13:
			pidCode = 5; // muon
			break;
		default:
			pidCode = 6;  // something else?
	};

	return pidCode;
}

//_____________________________________________________________________________
void AliAnalysisTaskQAHighPtDeDx::ProcessMCTruthESD() 
{
	// Fill the special MC histogram with the MC truth info

	const Int_t nTracksMC = fMCStack->GetNtrack();

	for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {

		//Cuts
		if(!(fMCStack->IsPhysicalPrimary(iTracks)))
			continue;

		TParticle* trackMC = fMCStack->Particle(iTracks);

		TParticlePDG* pdgPart = trackMC ->GetPDG();
		Double_t chargeMC = pdgPart->Charge();

		if(chargeMC==0)
			continue;

		if (TMath::Abs(trackMC->Eta()) > fEtaCut )
			continue;

		Double_t etaMC = trackMC->Eta();
		Int_t pdgCode = trackMC->GetPdgCode();
		Short_t pidCodeMC = 0;
		pidCodeMC = GetPidCode(pdgCode);


		for(Int_t nh = 0; nh < 9; nh++) {

			if( etaMC > etaHigh[nh]/10.0 || etaMC < etaLow[nh]/10.0 )
				continue;

			hMcIn[0][nh]->Fill(trackMC->Pt());
			hMcIn[pidCodeMC][nh]->Fill(trackMC->Pt());


		}

	}//MC track loop



}

//_____________________________________________________________________________
void AliAnalysisTaskQAHighPtDeDx::ProcessMCTruthAOD() 
{
	// Fill the special MC histogram with the MC truth info


	const Int_t nTracksMC = fMCArray->GetEntriesFast();

	for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {

		AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));

		if(!trackMC){
			AliError("Cannot get MC particle");
			continue;
		}   

		//Cuts
		if(!(trackMC->IsPhysicalPrimary()))
			continue;


		Double_t chargeMC = trackMC->Charge();
		if(chargeMC==0)
			continue;


		if (TMath::Abs(trackMC->Eta()) > fEtaCut )
			continue;

		Double_t etaMC = trackMC->Eta();
		Int_t pdgCode = trackMC->GetPdgCode();
		Short_t pidCodeMC = 0;
		pidCodeMC = GetPidCode(pdgCode);

		for(Int_t nh = 0; nh < 9; nh++) {

			if( etaMC > etaHigh[nh]/10.0 || etaMC < etaLow[nh]/10.0 )
				continue;

			hMcIn[0][nh]->Fill(trackMC->Pt());
			hMcIn[pidCodeMC][nh]->Fill(trackMC->Pt());


		}

	}//MC track loop


}


//_____________________________________________________________________________
Short_t AliAnalysisTaskQAHighPtDeDx::GetPythiaEventProcessType(Int_t pythiaType) {
	//
	// Get the process type of the event.  PYTHIA
	//
	// source PWG0   dNdpt 

	Short_t globalType = -1; //init

	if(pythiaType==92||pythiaType==93){
		globalType = 2; //single diffractive
	}
	else if(pythiaType==94){
		globalType = 3; //double diffractive
	}
	else {
		globalType = 1;  //non diffractive
	}
	return globalType;
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskQAHighPtDeDx::GetDPMjetEventProcessType(Int_t dpmJetType) {
	//
	// get the process type of the event.  PHOJET
	//
	//source PWG0   dNdpt 
	// can only read pythia headers, either directly or from cocktalil header
	Short_t globalType = -1;

	if (dpmJetType == 1 || dpmJetType == 4) { // explicitly inelastic plus central diffraction
		globalType = 1;
	}
	else if (dpmJetType==5 || dpmJetType==6) {
		globalType = 2;
	}
	else if (dpmJetType==7) {
		globalType = 3;
	}
	return globalType;
}

//_____________________________________________________________________________
ULong64_t AliAnalysisTaskQAHighPtDeDx::GetEventIdAsLong(AliVHeader* header) const
{
	// To have a unique id for each event in a run!
	// Modified from AliRawReader.h
	return ((ULong64_t)header->GetBunchCrossNumber()+
			(ULong64_t)header->GetOrbitNumber()*3564+
			(ULong64_t)header->GetPeriodNumber()*16777215*3564);
}


//____________________________________________________________________
TParticle* AliAnalysisTaskQAHighPtDeDx::FindPrimaryMother(AliStack* stack, Int_t label)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// Taken from AliPWG0Helper class
	//

	Int_t motherLabel = FindPrimaryMotherLabel(stack, label);
	if (motherLabel < 0)
		return 0;

	return stack->Particle(motherLabel);
}

//____________________________________________________________________
Int_t AliAnalysisTaskQAHighPtDeDx::FindPrimaryMotherLabel(AliStack* stack, Int_t label)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// returns its label
	//
	// Taken from AliPWG0Helper class
	//
	const Int_t nPrim  = stack->GetNprimary();

	while (label >= nPrim) {


		TParticle* particle = stack->Particle(label);
		if (!particle) {

			AliDebugGeneral("FindPrimaryMotherLabel", AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack.", label));
			return -1;
		}

		// find mother
		if (particle->GetMother(0) < 0) {

			AliDebugGeneral("FindPrimaryMotherLabel", AliLog::kError, Form("UNEXPECTED: Could not find mother of secondary particle %d.", label));
			return -1;
		}

		label = particle->GetMother(0);
	}

	return label;
}

//____________________________________________________________________
AliAODMCParticle* AliAnalysisTaskQAHighPtDeDx::FindPrimaryMotherAOD(AliAODMCParticle* startParticle)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// Taken from AliPWG0Helper class
	//

	AliAODMCParticle* mcPart = startParticle;

	while (mcPart)
	{

		if(mcPart->IsPrimary())
			return mcPart;

		Int_t mother = mcPart->GetMother();

		mcPart = dynamic_cast<AliAODMCParticle*>(fMCArray->At(mother));
	}

	return 0;
}
//____________________________________________________________________
TParticle* AliAnalysisTaskQAHighPtDeDx::FindPrimaryMotherV0(AliStack* stack, Int_t label)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// Taken from AliPWG0Helper class
	//

	Int_t nSteps = 0;

	Int_t motherLabel = FindPrimaryMotherLabelV0(stack, label, nSteps);
	if (motherLabel < 0)
		return 0;

	return stack->Particle(motherLabel);
}

//____________________________________________________________________
Int_t AliAnalysisTaskQAHighPtDeDx::FindPrimaryMotherLabelV0(AliStack* stack, Int_t label, Int_t& nSteps)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// returns its label
	//
	// Taken from AliPWG0Helper class
	//
	nSteps = 0;
	const Int_t nPrim  = stack->GetNprimary();

	while (label >= nPrim) {

		nSteps++; // 1 level down

		TParticle* particle = stack->Particle(label);
		if (!particle) {

			AliDebugGeneral("FindPrimaryMotherLabelV0", AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack.", label));
			return -1;
		}

		// find mother
		if (particle->GetMother(0) < 0) {

			AliDebugGeneral("FindPrimaryMotherLabelV0", AliLog::kError, Form("UNEXPECTED: Could not find mother of secondary particle %d.", label));
			return -1;
		}

		label = particle->GetMother(0);
	}

	return label;
}

//____________________________________________________________________
AliAODMCParticle* AliAnalysisTaskQAHighPtDeDx::FindPrimaryMotherAODV0(AliAODMCParticle* startParticle, Int_t& nSteps)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// Taken from AliPWG0Helper class
	//

	nSteps = 0;

	AliAODMCParticle* mcPart = startParticle;

	while (mcPart)
	{

		if(mcPart->IsPrimary())
			return mcPart;

		Int_t mother = mcPart->GetMother();

		mcPart = dynamic_cast<AliAODMCParticle*>(fMCArray->At(mother));
		nSteps++; // 1 level down
	}

	return 0;
}



//__________________________________________________________________
void AliAnalysisTaskQAHighPtDeDx::ProduceArrayTrksESD( AliESDEvent *ESDevent ){

	const Int_t nESDTracks = ESDevent->GetNumberOfTracks();
	//Int_t trackmult=0;


	Int_t multTPC = 0;

	//get multiplicity tpc only track cuts
	for(Int_t iT = 0; iT < nESDTracks; iT++) {

		AliESDtrack* esdTrack = ESDevent->GetTrack(iT);


		if(TMath::Abs(esdTrack->Eta()) > fEtaCut)
			continue;

		//only golden track cuts
		UInt_t selectDebug = 0;
		if (fTrackFilterTPC) {
			selectDebug = fTrackFilterTPC->IsSelected(esdTrack);
			if (!selectDebug) {
				continue;
			}
		}

		multTPC++;

	}



	for(Int_t iT = 0; iT < nESDTracks; iT++) {

		AliESDtrack* esdTrack = ESDevent->GetTrack(iT);


		if(TMath::Abs(esdTrack->Eta()) > fEtaCut)
			continue;

		//only golden track cuts
		UInt_t selectDebug = 0;
		if (fTrackFilterGolden) {
			selectDebug = fTrackFilterGolden->IsSelected(esdTrack);
			if (!selectDebug) {
				continue;
			}
		}

		Short_t ncl     = esdTrack->GetTPCsignalN();


		if(ncl<70)
			continue;
		Double_t eta  = esdTrack->Eta();
		Double_t phi  = esdTrack->Phi();
		Double_t momentum = esdTrack->P();
		Float_t  dedx    = esdTrack->GetTPCsignal();

		if(!PhiCut(esdTrack->Pt(), phi, esdTrack->Charge(), magf, &cutLow, &cutHigh))
			continue;


		//TOF

		Bool_t IsTOFout=kFALSE;
		if ((esdTrack->GetStatus()&AliESDtrack::kTOFout)==0)
			IsTOFout=kTRUE;
		Float_t lengthtrack=esdTrack->GetIntegratedLength();
		Float_t timeTOF=esdTrack->GetTOFsignal();
		Double_t inttime[5]={0,0,0,0,0}; 
		esdTrack->GetIntegratedTimes(inttime);// Returns the array with integrated times for each particle hypothesis
		Float_t beta = -99;
		if ( !IsTOFout ){
			if ( ( lengthtrack != 0 ) && ( timeTOF != 0) )
				beta = inttime[0] / timeTOF;
		}

		Short_t pidCode     = 0; 

		if(fAnalysisMC) {

			const Int_t label = TMath::Abs(esdTrack->GetLabel());
			TParticle* mcTrack = fMCStack->Particle(label);	    
			if (mcTrack){

				Int_t pdgCode = mcTrack->GetPdgCode();
				pidCode = GetPidCode(pdgCode);

			}

		}


		if( momentum <= 0.6 && momentum >= 0.4  ){//only p:0.4-0.6 GeV, pion MIP

			if( dedx < DeDxMIPMax && dedx > DeDxMIPMin ){	  
				if(momentum<0.6&&momentum>0.4){
					hMIPVsEta->Fill(eta,dedx);
					pMIPVsEta->Fill(eta,dedx);
				}
			}
			if( dedx > DeDxMIPMax+1 && dedx < 95 ){
				if(TMath::Abs(beta-1)<0.1){
					hPlateauVsEta->Fill(eta,dedx);
					pPlateauVsEta->Fill(eta,dedx); 
				}
			}
		}


		for(Int_t nh = 0; nh < 9; nh++) {

			if( eta > etaHigh[nh]/10.0 || eta < etaLow[nh]/10.0 )
				continue;


			if(fAnalysisMC){
				hMcOut[0][nh]->Fill(esdTrack->Pt());
				hMcOut[pidCode][nh]->Fill(esdTrack->Pt());
			}

			histAllCh[nh]->Fill(momentum, dedx);
			if(beta>1){
				histPiTof[nh]->Fill(momentum, dedx);
				histpPiTof[nh]->Fill(momentum);
			}

			if( momentum <= 0.6 && momentum >= 0.4  ){//only p:0.4-0.6 GeV, pion MIP
				//Fill  pion MIP, before calibration
				if( dedx < DeDxMIPMax && dedx > DeDxMIPMin ){	  
					hMIPVsPhi[nh]->Fill(phi,dedx);
					pMIPVsPhi[nh]->Fill(phi,dedx);


					hMIPVsNch[nh]->Fill(multTPC,dedx);
					pMIPVsNch[nh]->Fill(multTPC,dedx);

				}

				//Fill electrons, before calibration
				if( dedx > DeDxMIPMax+1 && dedx < 95 ){
					if(TMath::Abs(beta-1)<0.1){
						hPlateauVsPhi[nh]->Fill(phi,dedx);
						pPlateauVsPhi[nh]->Fill(phi,dedx);
					}
				}

			}


		}//end loop over eta intervals





	}//end of track loop




}
//__________________________________________________________________
void AliAnalysisTaskQAHighPtDeDx::ProduceArrayTrksAOD( AliAODEvent *AODevent ){


	Int_t nAODTracks = AODevent->GetNumberOfTracks();
	Int_t multTPC = 0;

	//get multiplicity tpc only track cuts
	for(Int_t iT = 0; iT < nAODTracks; iT++) {

		AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(AODevent->GetTrack(iT));
		if(!aodTrack) AliFatal("Not a standard AOD");

		if(TMath::Abs(aodTrack->Eta()) > fEtaCut)
			continue;


		if (fTrackFilterTPC) {
			// TPC only cuts is bit 1
			if(!aodTrack->TestFilterBit(1))
				continue;
		}

		multTPC++;

	}


	for(Int_t iT = 0; iT < nAODTracks; iT++) {

		AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(AODevent->GetTrack(iT));
		if(!aodTrack) AliFatal("Not a standard AOD");

		if (fTrackFilterGolden) {     
			// "Global track RAA analysis QM2011 + Chi2ITS<36"; bit 1024
			if(!aodTrack->TestFilterBit(1024))
				continue;
		}


		if(TMath::Abs(aodTrack->Eta()) > fEtaCut)
			continue;


		Double_t eta  = aodTrack->Eta();
		Double_t phi  = aodTrack->Phi();
		Double_t momentum = aodTrack->P();


		if(!PhiCut(aodTrack->Pt(), phi, aodTrack->Charge(), magf, &cutLow, &cutHigh))
			continue;



		AliAODPid* aodPid = aodTrack->GetDetPid();
		Short_t ncl     = -10;
		Float_t dedx    = -10;

		//TOF    
		Float_t beta = -99;


		if(aodPid) {
			ncl     = aodPid->GetTPCsignalN();
			dedx    = aodPid->GetTPCsignal();
			//TOF
			Bool_t IsTOFout=kFALSE;
			Float_t lengthtrack=-999;//in aod we do not have that information, beta must be: beta=inttime/timeTOF 
			Float_t timeTOF=-999;

			if ((aodTrack->GetStatus()&AliESDtrack::kTOFout)==0)
				IsTOFout=kTRUE;

			lengthtrack=-999;//in aod we do not have that information, beta must be: beta=inttime/timeTOF 

			timeTOF=aodPid->GetTOFsignal();

			Double_t inttime[5]={0,0,0,0,0}; 
			aodPid->GetIntegratedTimes(inttime);// Returns the array with integrated times for each particle hypothesis


			if ( !IsTOFout ){
				if ( ( lengthtrack != 0 ) && ( timeTOF != 0) )
					beta = inttime[0] / timeTOF;
			}

		}


		if(ncl<70)
			continue;


		Short_t pidCode     = 0; 

		if(fAnalysisMC) {

			const Int_t label = TMath::Abs(aodTrack->GetLabel());
			AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(fMCArray->At(label));

			if (mcTrack){

				Int_t pdgCode = mcTrack->GetPdgCode();
				pidCode = GetPidCode(pdgCode);

			}

		}

		if( momentum <= 0.6 && momentum >= 0.4  ){//only p:0.4-0.6 GeV, pion MIP
			if( dedx < DeDxMIPMax && dedx > DeDxMIPMin ){	  
				if(momentum<0.6&&momentum>0.4){
					hMIPVsEta->Fill(eta,dedx);
					pMIPVsEta->Fill(eta,dedx);
				}
			}
			if( dedx > DeDxMIPMax+1 && dedx < 95 ){
				if(TMath::Abs(beta-1)<0.1){
					hPlateauVsEta->Fill(eta,dedx);
					pPlateauVsEta->Fill(eta,dedx); 
				}
			}
		}


		for(Int_t nh = 0; nh < 9; nh++) {

			if( eta > etaHigh[nh]/10.0 || eta < etaLow[nh]/10.0 )
				continue;


			if(fAnalysisMC){
				hMcOut[0][nh]->Fill(aodTrack->Pt());
				hMcOut[pidCode][nh]->Fill(aodTrack->Pt());
			}

			histAllCh[nh]->Fill(momentum, dedx);
			if(beta>1){
				histPiTof[nh]->Fill(momentum, dedx);
				histpPiTof[nh]->Fill(momentum);
			}

			if( momentum <= 0.6 && momentum >= 0.4  ){//only p:0.4-0.6 GeV, pion MIP
				//Fill  pion MIP, before calibration
				if( dedx < DeDxMIPMax && dedx > DeDxMIPMin ){	  
					hMIPVsPhi[nh]->Fill(phi,dedx);
					pMIPVsPhi[nh]->Fill(phi,dedx);


					hMIPVsNch[nh]->Fill(multTPC,dedx);
					pMIPVsNch[nh]->Fill(multTPC,dedx);

				}

				//Fill electrons, before calibration
				if( dedx > DeDxMIPMax+1 && dedx < 95 ){
					if(TMath::Abs(beta-1)<0.1){
						hPlateauVsPhi[nh]->Fill(phi,dedx);
						pPlateauVsPhi[nh]->Fill(phi,dedx);
					}
				}

			}

		}//end loop over eta intervals





	}//end of track loop







}
//__________________________________________________
void AliAnalysisTaskQAHighPtDeDx::ProduceArrayV0ESD( AliESDEvent *ESDevent ){
	Int_t nv0s = ESDevent->GetNumberOfV0s();

	const AliESDVertex *myBestPrimaryVertex = ESDevent->GetPrimaryVertex();
	if (!myBestPrimaryVertex) return;
	if (!(myBestPrimaryVertex->GetStatus())) return;
	Double_t  lPrimaryVtxPosition[3];
	myBestPrimaryVertex->GetXYZ(lPrimaryVtxPosition);

	Double_t  lPrimaryVtxCov[6];
	myBestPrimaryVertex->GetCovMatrix(lPrimaryVtxCov);
	Double_t  lPrimaryVtxChi2 = myBestPrimaryVertex->GetChi2toNDF();

	AliAODVertex* myPrimaryVertex = new AliAODVertex(lPrimaryVtxPosition, lPrimaryVtxCov, lPrimaryVtxChi2, NULL, -1, AliAODVertex::kPrimary);



	//
	// LOOP OVER V0s, K0s, L, AL
	//


	for (Int_t iV0 = 0; iV0 < nv0s; iV0++) {

		// This is the begining of the V0 loop  
		AliESDv0 *esdV0 = ESDevent->GetV0(iV0);
		if (!esdV0) continue;

		//check onfly status
		if(esdV0->GetOnFlyStatus()!=0)
			continue;


		// AliESDTrack (V0 Daughters)
		UInt_t lKeyPos = (UInt_t)TMath::Abs(esdV0->GetPindex());
		UInt_t lKeyNeg = (UInt_t)TMath::Abs(esdV0->GetNindex());

		AliESDtrack *pTrack = ESDevent->GetTrack(lKeyPos);
		AliESDtrack *nTrack = ESDevent->GetTrack(lKeyNeg);
		if (!pTrack || !nTrack) {
			Printf("ERROR: Could not retreive one of the daughter track");
			continue;
		}

		// Remove like-sign
		if (pTrack->GetSign() == nTrack->GetSign()) {
			//cout<< "like sign, continue"<< endl;
			continue;
		} 

		// Eta cut on decay products
		if(TMath::Abs(pTrack->Eta()) > fEtaCut || TMath::Abs(nTrack->Eta()) > fEtaCut)
			continue;


		// Check if switch does anything!
		Bool_t isSwitched = kFALSE;
		if (pTrack->GetSign() < 0) { // switch

			isSwitched = kTRUE;
			AliESDtrack* helpTrack = nTrack;
			nTrack = pTrack;
			pTrack = helpTrack;
		}	

		AliKFVertex primaryVtxKF( *myPrimaryVertex );
		AliKFParticle::SetField(ESDevent->GetMagneticField());

		// Also implement switch here!!!!!!
		AliKFParticle* negEKF  = 0; // e-
		AliKFParticle* posEKF  = 0; // e+
		AliKFParticle* negPiKF = 0; // pi -
		AliKFParticle* posPiKF = 0; // pi +
		AliKFParticle* posPKF  = 0; // p
		AliKFParticle* negAPKF = 0; // p-bar

		if(!isSwitched) {
			negEKF  = new AliKFParticle( *(esdV0->GetParamN()) , 11);
			posEKF  = new AliKFParticle( *(esdV0->GetParamP()) ,-11);
			negPiKF = new AliKFParticle( *(esdV0->GetParamN()) ,-211);
			posPiKF = new AliKFParticle( *(esdV0->GetParamP()) , 211);
			posPKF  = new AliKFParticle( *(esdV0->GetParamP()) , 2212);
			negAPKF = new AliKFParticle( *(esdV0->GetParamN()) ,-2212);
		} else { // switch + and - 
			negEKF  = new AliKFParticle( *(esdV0->GetParamP()) , 11);
			posEKF  = new AliKFParticle( *(esdV0->GetParamN()) ,-11);
			negPiKF = new AliKFParticle( *(esdV0->GetParamP()) ,-211);
			posPiKF = new AliKFParticle( *(esdV0->GetParamN()) , 211);
			posPKF  = new AliKFParticle( *(esdV0->GetParamN()) , 2212);
			negAPKF = new AliKFParticle( *(esdV0->GetParamP()) ,-2212);
		}

		AliKFParticle v0GKF;  // Gamma e.g. from pi0
		v0GKF+=(*negEKF);
		v0GKF+=(*posEKF);
		v0GKF.SetProductionVertex(primaryVtxKF);

		AliKFParticle v0K0sKF; // K0 short
		v0K0sKF+=(*negPiKF);
		v0K0sKF+=(*posPiKF);
		v0K0sKF.SetProductionVertex(primaryVtxKF);

		AliKFParticle v0LambdaKF; // Lambda
		v0LambdaKF+=(*negPiKF);
		v0LambdaKF+=(*posPKF);	
		v0LambdaKF.SetProductionVertex(primaryVtxKF);

		AliKFParticle v0AntiLambdaKF; // Lambda-bar
		v0AntiLambdaKF+=(*posPiKF);
		v0AntiLambdaKF+=(*negAPKF);
		v0AntiLambdaKF.SetProductionVertex(primaryVtxKF);

		Double_t dmassG     = v0GKF.GetMass();
		Double_t dmassK     = v0K0sKF.GetMass()-0.498;
		Double_t dmassL     = v0LambdaKF.GetMass()-1.116;
		Double_t dmassAL    = v0AntiLambdaKF.GetMass()-1.116;


		for( Int_t case_v0 = 0; case_v0 < 2; ++case_v0 ){


			switch(case_v0){
				case 0:{    

					       Bool_t fillPos = kFALSE;
					       Bool_t fillNeg = kFALSE;

					       if(dmassG < 0.1)
						       continue;

					       if(dmassK>0.01 && dmassL>0.01 && dmassAL>0.01){
						       continue;
					       }

					       if(dmassL<0.01){
						       fillPos = kTRUE;
					       }
					       if(dmassAL<0.01) {
						       if(fillPos)
							       continue;
						       fillNeg = kTRUE;
					       }

					       if(dmassK<0.01) {
						       if(fillPos||fillNeg)
							       continue;
						       fillPos = kTRUE;
						       fillNeg = kTRUE;
					       }


					       for(Int_t j = 0; j < 2; j++) {

						       AliESDtrack* track = 0;

						       if(j==0) {

							       if(fillNeg)
								       track = nTrack;
							       else
								       continue;
						       } else {

							       if(fillPos)
								       track = pTrack;
							       else
								       continue;
						       }

						       if(track->GetTPCsignalN()<=70)continue;
						       Double_t phi     = track->Phi();

						       if(!PhiCut(track->Pt(), phi, track->Charge(), magf, &cutLow, &cutHigh))
							       continue;


						       Double_t eta  = track->Eta();
						       Double_t momentum = track->Pt();
						       Double_t dedx = track->GetTPCsignal();

						       if(fillPos&&fillNeg){


							       if( dedx < DeDxMIPMax && dedx > DeDxMIPMin ){	  
								       if(momentum<0.6&&momentum>0.4){
									       hMIPVsEtaV0s->Fill(eta,dedx);
									       pMIPVsEtaV0s->Fill(eta,dedx);
								       }
							       }


						       }

						       for(Int_t nh = 0; nh < nHists; nh++) {



							       if( eta < etaLow[nh]/10.0 || eta > etaHigh[nh]/10.0 )
								       continue;

							       if(fillPos&&fillNeg){

								       histPiV0[nh]->Fill(momentum, dedx);	    
								       histpPiV0[nh]->Fill(momentum);

							       }
							       else{

								       histPV0[nh]->Fill(momentum, dedx);
								       histpPV0[nh]->Fill(momentum);

							       }

						       }

					       }//end loop over two tracks

				       };
				       break;

				case 1:{//gammas    

					       Bool_t fillPos = kFALSE;
					       Bool_t fillNeg = kFALSE;



					       if(dmassK>0.01 && dmassL>0.01 && dmassAL>0.01) {
						       if(dmassG<0.01 && dmassG>0.0001) {


							       if( TMath::Abs(nTrack->GetTPCsignal() - 85.0) < 5)
								       fillPos = kTRUE;
							       if( TMath::Abs(pTrack->GetTPCsignal() - 85.0) < 5)
								       fillNeg = kTRUE;

						       } else {
							       continue;
						       }
					       }

					       //cout<<"fillPos="<<fillPos<<endl;

					       if(fillPos == kTRUE && fillNeg == kTRUE)      
						       continue;


					       AliESDtrack* track = 0;
					       if(fillNeg)
						       track = nTrack;
					       else if(fillPos)
						       track = pTrack;
					       else
						       continue;

					       Double_t dedx  = track->GetTPCsignal();
					       Double_t eta  = track->Eta();
					       Double_t phi  = track->Phi();
					       Double_t momentum = track->P();

					       if(track->GetTPCsignalN()<=70)continue;

					       if(!PhiCut(track->Pt(), phi, track->Charge(), magf, &cutLow, &cutHigh))
						       continue;

					       for(Int_t nh = 0; nh < nHists; nh++) {

						       if( eta < etaLow[nh]/10.0 || eta > etaHigh[nh]/10.0 )
							       continue;

						       histEV0[nh]->Fill(momentum, dedx);

					       }

				       };
				       break;


			}//end switch

		}//end loop over case V0


		// clean up loop over v0

		delete negPiKF;
		delete posPiKF;
		delete posPKF;
		delete negAPKF;



	}


	delete myPrimaryVertex;


}
//__________________________________________________________________________
void AliAnalysisTaskQAHighPtDeDx::ProduceArrayV0AOD( AliAODEvent *AODevent ){
	Int_t nv0s = AODevent->GetNumberOfV0s();

	AliAODVertex *myBestPrimaryVertex = AODevent->GetPrimaryVertex();
	if (!myBestPrimaryVertex) return;



	// ################################
	// #### BEGINNING OF V0 CODE ######
	// ################################
	// This is the begining of the V0 loop  
	for (Int_t iV0 = 0; iV0 < nv0s; iV0++) {
		AliAODv0 *aodV0 = AODevent->GetV0(iV0);
		if (!aodV0) continue;


		//check onfly status
		if(aodV0->GetOnFlyStatus()!=0)
			continue;

		// AliAODTrack (V0 Daughters)
		AliAODVertex* vertex = aodV0->GetSecondaryVtx();
		if (!vertex) {
			Printf("ERROR: Could not retrieve vertex");
			continue;
		}

		AliAODTrack *pTrack = (AliAODTrack*)vertex->GetDaughter(0);
		AliAODTrack *nTrack = (AliAODTrack*)vertex->GetDaughter(1);
		if (!pTrack || !nTrack) {
			Printf("ERROR: Could not retrieve one of the daughter track");
			continue;
		}

		// Remove like-sign
		if (pTrack->Charge() == nTrack->Charge()) {
			//cout<< "like sign, continue"<< endl;
			continue;
		} 

		// Make sure charge ordering is ok
		if (pTrack->Charge() < 0) {
			AliAODTrack* helpTrack = pTrack;
			pTrack = nTrack;
			nTrack = helpTrack;
		} 

		// Eta cut on decay products
		if(TMath::Abs(pTrack->Eta()) > fEtaCut || TMath::Abs(nTrack->Eta()) > fEtaCut)
			continue;


		Double_t dmassG  = aodV0->InvMass2Prongs(0,1,11,11);
		Double_t dmassK  = aodV0->MassK0Short()-0.498;
		Double_t dmassL  = aodV0->MassLambda()-1.116;
		Double_t dmassAL = aodV0->MassAntiLambda()-1.116;

		for( Int_t case_v0 = 0; case_v0 < 2; ++case_v0 ){


			switch(case_v0){
				case 0:{   
					       Bool_t fillPos = kFALSE;
					       Bool_t fillNeg = kFALSE;


					       if(dmassG < 0.1)
						       continue;

					       if(dmassK>0.01 && dmassL>0.01 && dmassAL>0.01){
						       continue;
					       }

					       if(dmassL<0.01){
						       fillPos = kTRUE;
					       }
					       if(dmassAL<0.01) {
						       if(fillPos)
							       continue;
						       fillNeg = kTRUE;
					       }

					       if(dmassK<0.01) {
						       if(fillPos||fillNeg)
							       continue;
						       fillPos = kTRUE;
						       fillNeg = kTRUE;
					       }



					       for(Int_t j = 0; j < 2; j++) {

						       AliAODTrack* track = 0;

						       if(j==0) {

							       if(fillNeg)
								       track = nTrack;
							       else
								       continue;
						       } else {

							       if(fillPos)
								       track = pTrack;
							       else
								       continue;
						       }

						       if(track->GetTPCsignalN()<=70)continue;

						       Double_t phi     = track->Phi();

						       if(!PhiCut(track->Pt(), phi, track->Charge(), magf, &cutLow, &cutHigh))
							       continue;


						       //if(!PhiCut(pt, phi, charge, magf, &cutLow, &cutHigh, hPhi))
						       //	continue;

						       Double_t eta  = track->Eta();
						       Double_t momentum = track->Pt();
						       Double_t dedx = track->GetTPCsignal();

						       if(fillPos&&fillNeg){


							       if( dedx < DeDxMIPMax && dedx > DeDxMIPMin ){	  
								       if(momentum<0.6&&momentum>0.4){
									       hMIPVsEtaV0s->Fill(eta,dedx);
									       pMIPVsEtaV0s->Fill(eta,dedx);
								       }
							       }


						       }

						       for(Int_t nh = 0; nh < nHists; nh++) {



							       if( eta < etaLow[nh]/10.0 || eta > etaHigh[nh]/10.0 )
								       continue;

							       if(fillPos&&fillNeg){

								       histPiV0[nh]->Fill(momentum, dedx);	    
								       histpPiV0[nh]->Fill(momentum);	    

							       }
							       else{

								       histPV0[nh]->Fill(momentum, dedx);
								       histpPV0[nh]->Fill(momentum);

							       }

						       }


					       }//end loop over two tracks
				       };
				       break;

				case 1:{//gammas    

					       Bool_t fillPos = kFALSE;
					       Bool_t fillNeg = kFALSE;

					       if(dmassK>0.01 && dmassL>0.01 && dmassAL>0.01) {
						       if(dmassG<0.01 && dmassG>0.0001) {

							       if( TMath::Abs(nTrack->GetTPCsignal() - 85.0) < 5)
								       fillPos = kTRUE;
							       if( TMath::Abs(pTrack->GetTPCsignal() - 85.0) < 5)
								       fillNeg = kTRUE;

						       } else {
							       continue;
						       }
					       }


					       if(fillPos == kTRUE && fillNeg == kTRUE)      
						       continue;


					       AliAODTrack* track = 0;
					       if(fillNeg)
						       track = nTrack;
					       else if(fillPos)
						       track = pTrack;
					       else
						       continue;

					       Double_t dedx  = track->GetTPCsignal();
					       Double_t eta  = track->Eta();
					       Double_t phi  = track->Phi();
					       Double_t momentum = track->P();

					       if(track->GetTPCsignalN()<=70)continue;

					       if(!PhiCut(track->Pt(), phi, track->Charge(), magf, &cutLow, &cutHigh))
						       continue;

					       for(Int_t nh = 0; nh < nHists; nh++) {

						       if( eta < etaLow[nh]/10.0 || eta > etaHigh[nh]/10.0 )
							       continue;

						       histEV0[nh]->Fill(momentum, dedx);

					       }

				       };
				       break;


			}//end switch
		}//end loop over V0s cases

	}//end loop over v0's




}
Bool_t AliAnalysisTaskQAHighPtDeDx::PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t   mag, TF1* phiCutLow, TF1* phiCutHigh)
{
	if(pt < 2.0)
		return kTRUE;

	//Double_t phi = track->Phi();
	if(mag < 0)    // for negatve polarity field
		phi = TMath::TwoPi() - phi;
	if(q < 0) // for negatve charge
		phi = TMath::TwoPi()-phi;

	phi += TMath::Pi()/18.0; // to center gap in the middle
	phi = fmod(phi, TMath::Pi()/9.0);

	if(phi<phiCutHigh->Eval(pt) 
			&& phi>phiCutLow->Eval(pt))
		return kFALSE; // reject track

	hPhi->Fill(pt, phi);

	return kTRUE;
}
