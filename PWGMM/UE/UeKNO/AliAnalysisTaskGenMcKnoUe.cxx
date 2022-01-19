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
 * Author: Ahsan Mehmood Khan(ahsan.mehmood.khan@cern.ch)                 * 
 *         Feng Fan (Feng.Fan@cern.ch)	                                  *
 *         Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)                *
 *         Last modification: 25/11/2021                                  * 
 **************************************************************************/
//_____ ROOT headers
#include <TList.h>
#include <TTree.h>
#include <TMath.h>
#include <TFile.h>
#include <TF1.h>
#include <TRandom.h>
#include <TTreeStream.h>
#include "TChain.h"
#include <THnSparse.h>
#include <TDatabasePDG.h>
#include "TObjArray.h"
#include <TClonesArray.h>
#include "TH1.h"
#include "TH2.h"
#include "TH3D.h"
//_____ ALIROOT headers
#include "AliAnalysisTask.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliVParticle.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
//_____ Additional includes
#include "AliVEvent.h"
#include "AliGenEventHeader.h"

//_____ AnalysisTask headers
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskGenMcKnoUe.h"
//_____ STL includes
#include <iostream>
using namespace std;

TF1 * frho;
const Int_t nRegDphi = 4;
const Int_t nRegRho  = 3;
const Char_t * NameOfRegion[3]={"NS","AS","TS"};
const Char_t * NameOfRegionDPhi[nRegDphi]={"NS","AS","TS","TS2"};
const Char_t * NameOfRegionRho[nRegRho]={"allrho","lowrho","highrho"};
const Int_t NchNBins = 2000;
const Int_t NchNBinsRho = 600;
const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;
const Double_t  ptamin = 0.15;
const Double_t  ptamax = 20.0;
Int_t Nmpi = -1;

ClassImp(AliAnalysisTaskGenMcKnoUe)

	AliAnalysisTaskGenMcKnoUe::AliAnalysisTaskGenMcKnoUe(): 
		AliAnalysisTaskSE(),
		fMC(0x0),
		fMcHandler(0x0),
		fMCStack(0x0),
		fFirstPart(kTRUE),
		fGenerator(0),
		fEtaCut(0.8),
		fEtaCutRho(4.0),
		fnEtaBinsRho(1),
		fnPhiBinsRho(1),
		fIsPP(kTRUE),
		fPtMin(0.15),
		fGenLeadPhi(0x0),
		fGenLeadEta(0x0),
		fGenLeadPt(0x0),
		fGenLeadIn(0x0),
		hPtLeadingTrue(0x0),
		hPtLVsV0A(0x0),
		hnchmpi(0x0),
		hnchmpirho(0x0),
		hnchrho(0x0),
		hmpirho(0x0),
		hphiKNO(0x0),
		hphiKNO1(0x0),
		hphiKNO2(0x0),
		hNchTforKNOana(0x0),
		hNchTforKNOanaMin(0x0),
		hNchTforKNOanaMax(0x0),
		fOutputList(0)
{
	for(Int_t i=0;i<3;++i) 
		hPtVsPtLeadingTrue[i]=0;
	for(Int_t i_rho=0;i_rho<nRegRho;++i_rho){
		hPtLeadingRho[i_rho]=0;// 0: all, 1: low rho, 2: high rho
		hNchRho[i_rho]=0;
		hEtaLeadingRho[i_rho]=0;
		hDetaDphiRho[i_rho]=0;
		hDetaDphiRhoWideEta[i_rho]=0;
		for(Int_t i_dphi=0;i_dphi<nRegDphi;++i_dphi)
			hNchPtPidRho[i_dphi][i_rho]=0; // region, rho 
	}

}
//_____________________________________________________________________________
AliAnalysisTaskGenMcKnoUe::AliAnalysisTaskGenMcKnoUe(const char* name): 
	AliAnalysisTaskSE(name),
	fMC(0x0),
	fMcHandler(0x0),
	fMCStack(0x0),
	fFirstPart(kTRUE),
	fGenerator(0),
	fEtaCut(0.8),
	fEtaCutRho(4.0),
	fnEtaBinsRho(1),
	fnPhiBinsRho(1),
	fIsPP(kTRUE),
	fPtMin(0.15),
	fGenLeadPhi(0x0),
	fGenLeadEta(0x0),
	fGenLeadPt(0x0),
	fGenLeadIn(0x0),
	hPtLeadingTrue(0x0),
	hPtLVsV0A(0x0),
	hnchmpi(0x0),
	hnchmpirho(0x0),
	hnchrho(0x0),
	hmpirho(0x0),
	hphiKNO(0x0),
	hphiKNO1(0x0),
	hphiKNO2(0x0),
	hNchTforKNOana(0x0),
	hNchTforKNOanaMin(0x0),
	hNchTforKNOanaMax(0x0),
	fOutputList(0)
{
	for(Int_t i=0;i<3;++i)
		hPtVsPtLeadingTrue[i]=0;
	for(Int_t i_rho=0;i_rho<nRegRho;++i_rho){
		hPtLeadingRho[i_rho]=0;// 0: all, 1: low rho, 2: high rho
		hNchRho[i_rho]=0;
		hEtaLeadingRho[i_rho]=0;
		hDetaDphiRho[i_rho]=0;
		hDetaDphiRhoWideEta[i_rho]=0;
		for(Int_t i_dphi=0;i_dphi<nRegDphi;++i_dphi)
			hNchPtPidRho[i_dphi][i_rho]=0; // region, rho 
	}

	// constructor
	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskGenMcKnoUe::~AliAnalysisTaskGenMcKnoUe()
{
	// destructor
	if(fOutputList) {
		delete fOutputList;
		fOutputList = 0x0;
	}

}
//_____________________________________________________________________________
void AliAnalysisTaskGenMcKnoUe::UserCreateOutputObjects()
{
	// ### Analysis output
	fOutputList = new TList();
	fOutputList->SetOwner(kTRUE);
	// create output objects
	Double_t NchBins[NchNBins+1];
	for(Int_t i_nch=0; i_nch<NchNBins+1; ++i_nch){
		NchBins[i_nch]=0;
		if(i_nch<NchNBins)
			NchBins[i_nch]=i_nch*1.0-0.5;
		else
			NchBins[i_nch]=NchNBins*1.0+0.5;
	}



	const Int_t pTNBins = 36;
	Double_t pTNBins1[pTNBins+1] = {
		0.0,  0.1,  0.15,  0.2,  0.25,  0.3,  0.35,  0.4,  0.45,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5, 5.0, 6.0,    7.0,  8.0,  9.0,  10.0,  12.0,  14.0,  16.0,  18.0,  20.0,  25.0,  30.0,  40.0,  50.0
	};

	const Int_t pTNBinsL = 24;
	Double_t pTNBins1L[pTNBinsL+1] = {
		0.15, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50,
		5.00, 6.00, 7.00, 8.00, 9.00, 10.0, 12.0, 14.0, 16.0, 18.0,
		20.0, 25.0, 30.0, 40.0, 50.0
	};

	// Define array nch
	Double_t EtaBins[fnEtaBinsRho+1];
	Double_t deltaEta = (2.0*fEtaCutRho)/(1.0*fnEtaBinsRho);
	for(Int_t i_eta=0;i_eta<fnEtaBinsRho+1;++i_eta){
		EtaBins[i_eta]=0;
		if(i_eta<fnEtaBinsRho)
			EtaBins[i_eta]=i_eta*deltaEta - 1.0*fEtaCutRho;
		else
			EtaBins[i_eta]= 1.0*fEtaCutRho;
	}
	Double_t PhiBins[fnPhiBinsRho+1];
	Double_t deltaPhi = (2.0*pi)/(1.0*fnPhiBinsRho);
	for(Int_t i_phi=0;i_phi<fnPhiBinsRho+1;++i_phi){
		PhiBins[i_phi]=0;
		if(i_phi<fnPhiBinsRho)
			PhiBins[i_phi]=i_phi*deltaPhi;
		else
			PhiBins[i_phi]= 2.0*pi;
	}


	if(fFirstPart){
		hPtLVsV0A = 0;
		hPtLVsV0A = new TH2D("hPtLVsV0A","",pTNBinsL,pTNBins1L,NchNBins,NchBins);
		fOutputList->Add(hPtLVsV0A);

		hPtLeadingTrue = new TH1D("hPtLeadingTrue","",pTNBinsL,pTNBins1L);
		fOutputList->Add(hPtLeadingTrue);
		if(fIsPP){

			hnchmpi = 0;
			hnchmpi = new TH2D("hnchmpi","; #it{N}_{ch}; #it{N}_{mpi}",600,-0.5,599.5,500,-0.5,4999.5);
			fOutputList->Add(hnchmpi);

			hnchmpirho = 0;
			hnchmpirho = new TH3F("hnchmpirho","; #it{N}_{ch}; #it{N}_{mpi}; #rho",600,-0.5,599.5,500,-0.5,4999.5,500,0.0,5.0);
			fOutputList->Add(hnchmpirho);

			hnchrho = 0;
			hnchrho = new TH2D("hnchrho","; #it{N}_{ch}; #rho",600,-0.5,599.5,1000,0,5.0);
			fOutputList->Add(hnchrho);

			hmpirho = 0;
			hmpirho = new TH2D("hmpirho","; #it{N}_{mpi}; #rho",500,-0.5,4999.5,1000,0,5.0);
			fOutputList->Add(hmpirho);

		}

		hphiKNO = 0;
		hphiKNO = new TH1D("hphiKNO","",64,-pi/2.0,3*pi/2.0);
		fOutputList->Add(hphiKNO);
		hphiKNO1 = 0;
		hphiKNO1 = new TH1D("hphiKNO1","",64,-pi/2.0,3*pi/2.0);
		fOutputList->Add(hphiKNO1);
		hphiKNO2 = 0;
		hphiKNO2 = new TH1D("hphiKNO2","",64,-pi/2.0,3*pi/2.0);
		fOutputList->Add(hphiKNO2);
		hNchTforKNOana = 0;
		hNchTforKNOana = new TH1D("hNchTforKNOana","",100,-0.5,99.5);
		fOutputList->Add(hNchTforKNOana);
		hNchTforKNOanaMin = 0;
		hNchTforKNOanaMin = new TH1D("hNchTforKNOanaMin","",100,-0.5,99.5);
		fOutputList->Add(hNchTforKNOanaMin);
		hNchTforKNOanaMax = 0;
		hNchTforKNOanaMax = new TH1D("hNchTforKNOanaMax","",100,-0.5,99.5);
		fOutputList->Add(hNchTforKNOanaMax);




		// UE analysis
		for(Int_t i=0;i<3;++i){

			hPtVsPtLeadingTrue[i] = new TH3D(Form("hPtVsPtLeadingTrue_%s",NameOfRegion[i]),"",pTNBinsL,pTNBins1L,pTNBins,pTNBins1,NchNBins,NchBins);
			fOutputList->Add(hPtVsPtLeadingTrue[i]);

		}
	}else{
		Double_t NchBinsRho[NchNBinsRho+1];
		for(Int_t i_nch=0; i_nch<NchNBinsRho+1; ++i_nch){
			NchBinsRho[i_nch]=0;
			if(i_nch<NchNBinsRho)
				NchBinsRho[i_nch]=i_nch*1.0-0.5;
			else
				NchBinsRho[i_nch]=NchNBinsRho*1.0+0.5;
		}
		const Int_t NPid = 3;
		Double_t Pid[NPid+1];
		for(Int_t i_pid=0;i_pid<NPid+1;++i_pid)
			Pid[i_pid]=i_pid;

		for(Int_t i_rho=0;i_rho<nRegRho;++i_rho){

			hPtLeadingRho[i_rho]=0;// 0: all, 1: low rho, 2: high rho
			hPtLeadingRho[i_rho]=new TH1D(Form("hPtLeading_%s",NameOfRegionRho[i_rho]),"",pTNBinsL,pTNBins1L);
			fOutputList->Add(hPtLeadingRho[i_rho]);

			hNchRho[i_rho]=0;
			hNchRho[i_rho]=new TH2D(Form("hNch_%s",NameOfRegionRho[i_rho]),"",NchNBinsRho,NchBinsRho,pTNBinsL,pTNBins1L);
			fOutputList->Add(hNchRho[i_rho]);

			hEtaLeadingRho[i_rho]=0;
			hEtaLeadingRho[i_rho]=new TH1D(Form("hEtaLeading_%s",NameOfRegionRho[i_rho]),"",80,-4,4);
			fOutputList->Add(hEtaLeadingRho[i_rho]);

			hDetaDphiRho[i_rho]=0;
			hDetaDphiRho[i_rho]=new TH2D(Form("hDetaDphi_%s",NameOfRegionRho[i_rho]),"",100,-5,5,64,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
			fOutputList->Add(hDetaDphiRho[i_rho]);

			hDetaDphiRhoWideEta[i_rho]=0;
			hDetaDphiRhoWideEta[i_rho]=new TH2D(Form("hDetaDphiWideEta_%s",NameOfRegionRho[i_rho]),"",100,-5,5,64,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
			fOutputList->Add(hDetaDphiRhoWideEta[i_rho]);


			for(Int_t i_dphi=0;i_dphi<nRegDphi;++i_dphi){
				hNchPtPidRho[i_dphi][i_rho]=0; // region, rho
				hNchPtPidRho[i_dphi][i_rho]=new TH3D(Form("hNchPtPid_%s_%s",NameOfRegionRho[i_rho],NameOfRegionDPhi[i_dphi]),"",NchNBinsRho,NchBinsRho,pTNBins,pTNBins1,NPid,Pid);
				fOutputList->Add(hNchPtPidRho[i_dphi][i_rho]);
			} 
		}
	}
	PostData(1, fOutputList);

}
//_____________________________________________________________________________
void AliAnalysisTaskGenMcKnoUe::UserExec(Option_t *)
{
	// ### Initialize
	if(!fFirstPart){
		frho = new TF1("frho","[0]*exp([1]*x)+[2]*x^[3]",100.0,2000.0);
		SetParametersRho(fEtaCutRho);
	}

	fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());

	// ### MC handler


	if(fMcHandler)
		fMC = fMcHandler->MCEvent();

	else { if(fDebug > 1) printf("AliAnalysisTaskGenUeNchTS::Handler() fMcHandler = NULL\n"); return; }

	// ### MC event

	if(!fMC){
		Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
		this->Dump();
		return;
	}
	fMCStack = ((AliMCEvent*)fMC)->Stack();

	if (!fMCStack) {
		Printf("ERROR: Could not retrieve MC stack \n");
		cout << "Name of the file with pb :" << endl; 
		return;
	}

	// ### MC event selection

	Bool_t isEventMCSelected = IsMCEventSelected(fMC);
	if( !isEventMCSelected )
		return;

	AliHeader* headerMC = fMC->Header();

	Bool_t isGoodVtxPosMC = kFALSE;

	AliGenEventHeader* genHeader = headerMC->GenEventHeader();
	TArrayF vtxMC(3); // primary vertex  MC 
	vtxMC[0]=9999; vtxMC[1]=9999;  vtxMC[2]=9999; //initialize with dummy
	if (genHeader) {
		genHeader->PrimaryVertex(vtxMC);
	}
	if(TMath::Abs(vtxMC[2])<=10)
		isGoodVtxPosMC = kTRUE;
	Nmpi = genHeader->NProduced();
	// Before trigger selection
	GetGenLeadingObject();// leading particle at gen level

	if(fFirstPart){
		if(fIsPP)
			MakeALICE3Analysis();
		if(isGoodVtxPosMC)
			if(fGenLeadPt>=fPtMin)
				GetGenUEObservables();
	}else{
		if(fIsPP)
			MakeALICE3AnalysisP2();
	}

	PostData(1, fOutputList);
	return;
}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskGenMcKnoUe::IsMCEventSelected(TObject* obj){

	Bool_t isSelected = kTRUE;

	AliMCEvent *event = 0x0;
	event = dynamic_cast<AliMCEvent*>(obj);
	if( !event )
		isSelected = kFALSE;

	return isSelected;

}
//-------------------------------------------------
void AliAnalysisTaskGenMcKnoUe::GetGenLeadingObject() {

	Double_t flPt = 0;// leading pT
	Double_t flPhi = 0;
	Int_t flIndex = 0;
	Double_t flEta = 0;

	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {

		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;
		if (!fMC->IsPhysicalPrimary(i)) continue;
		if (particle->Charge() == 0) continue;
		if(fFirstPart){
			if ( TMath::Abs(particle->Eta()) > fEtaCut )continue;}
		else{
			if ( TMath::Abs(particle->Eta()) > fEtaCutRho )continue;}
		if( particle->Pt() < fPtMin)continue;

		if (flPt<particle->Pt()){
			flPt = particle->Pt();
			flPhi = particle->Phi();
			flEta = particle->Eta();
			flIndex = i;

		}
	}

	fGenLeadPhi = flPhi;
	fGenLeadPt  = flPt;
	fGenLeadEta = flEta;
	fGenLeadIn  = flIndex;
}
//----------------------
void AliAnalysisTaskGenMcKnoUe::GetGenUEObservables(){

	Int_t multV0Aeta = 0;
	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {
		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;
		if (!fMC->IsPhysicalPrimary(i)) continue;
		if (particle->Charge() == 0) continue;
		if(particle->Pt()<=0)continue;
		if ( particle->Eta() > 2.8 && particle->Eta()<5.1 )
			multV0Aeta++;
		else
			continue;
	}
	hPtLVsV0A->Fill(fGenLeadPt,multV0Aeta*1.0);
	Int_t multForKNOana = 0;
	Int_t multForKNOana1 = 0;
	Int_t multForKNOana2 = 0;
	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {

		if(i==fGenLeadIn)
			continue;

		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;

		if (!fMC->IsPhysicalPrimary(i)) continue;
		if (particle->Charge() == 0) continue;
		if ( TMath::Abs(particle->Eta()) > fEtaCut )continue;
		if( particle->Pt() < fPtMin)continue;

		Double_t DPhi = DeltaPhi(particle->Phi(), fGenLeadPhi);

		// definition of the topological regions
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			hPtVsPtLeadingTrue[0]->Fill(fGenLeadPt,particle->Pt(),multV0Aeta*1.0);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			hPtVsPtLeadingTrue[1]->Fill(fGenLeadPt,particle->Pt(),multV0Aeta*1.0);
		}
		else{// transverse side
			hPtVsPtLeadingTrue[2]->Fill(fGenLeadPt,particle->Pt(),multV0Aeta*1.0);
			if( particle->Pt() >= 0.5 ){
				multForKNOana++;
				hphiKNO->Fill(DPhi);
				if(DPhi-pi>=pi/3.0 || DPhi<-pi/3.0){
					multForKNOana2++;
					hphiKNO2->Fill(DPhi);
				}
				else{
					multForKNOana1++;
					hphiKNO1->Fill(DPhi);
				}
			}


		}

	}
	hPtLeadingTrue->Fill(fGenLeadPt);
	if(fGenLeadPt>=5.0 && fGenLeadPt<40.0){
		hNchTforKNOana->Fill(multForKNOana);
		if(multForKNOana2>multForKNOana1){// multForKNOana1: max, multForKNOana1: min
			hNchTforKNOanaMin->Fill(multForKNOana1);// hist min
			hNchTforKNOanaMax->Fill(multForKNOana2);//hist max
		}
		else{
			hNchTforKNOanaMin->Fill(multForKNOana2);// hist min
			hNchTforKNOanaMax->Fill(multForKNOana1);//hist max
		}
	}

}
//_______________________________________________________
void AliAnalysisTaskGenMcKnoUe::MakeALICE3AnalysisP2(){

	Double_t EtaBins[fnEtaBinsRho+1];
	Double_t deltaEta = (2.0*fEtaCutRho)/(1.0*fnEtaBinsRho);
	for(Int_t i_eta=0;i_eta<fnEtaBinsRho+1;++i_eta){
		EtaBins[i_eta]=0;
		if(i_eta<fnEtaBinsRho)
			EtaBins[i_eta]=i_eta*deltaEta - 1.0*fEtaCutRho;
		else
			EtaBins[i_eta]= 1.0*fEtaCutRho;
	}
	Double_t PhiBins[fnPhiBinsRho+1];
	Double_t deltaPhi = (2.0*pi)/(1.0*fnPhiBinsRho);
	for(Int_t i_phi=0;i_phi<fnPhiBinsRho+1;++i_phi){
		PhiBins[i_phi]=0;
		if(i_phi<fnPhiBinsRho)
			PhiBins[i_phi]=i_phi*deltaPhi;
		else
			PhiBins[i_phi]= 2.0*pi;
	}

	Int_t NchLattice[fnEtaBinsRho][fnPhiBinsRho];
	for(Int_t i_eta=0;i_eta<fnEtaBinsRho;++i_eta)
		for(Int_t i_phi=0;i_phi<fnPhiBinsRho;++i_phi)
			NchLattice[i_eta][i_phi]=0;
	Double_t MpTLattice[fnEtaBinsRho][fnPhiBinsRho];
	for(Int_t i_eta=0;i_eta<fnEtaBinsRho;++i_eta)
		for(Int_t i_phi=0;i_phi<fnPhiBinsRho;++i_phi)
			MpTLattice[i_eta][i_phi]=0;

	Double_t totalpt=0;
	Int_t totalnch_forpt=0;
	Int_t nchtotal=0;

	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {

		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;
		Float_t eta_a = particle->Eta();
		Float_t phi_a = particle->Phi();
		Float_t pt_a  = particle->Pt();

		if (!fMC->IsPhysicalPrimary(i)) continue;
		if (particle->Charge() == 0) continue;
		if ( TMath::Abs(eta_a) > fEtaCutRho )continue;
		if ( pt_a <= 0 )continue;
		nchtotal++;
		if(pt_a<ptamin||pt_a>ptamax)continue;
		totalpt+=pt_a;
		totalnch_forpt++;
		// loop over all eta and phi intervals
		for(Int_t i_eta=0;i_eta<fnEtaBinsRho;++i_eta){
			for(Int_t i_phi=0;i_phi<fnPhiBinsRho;++i_phi){
				if(eta_a>=EtaBins[i_eta]&&eta_a<EtaBins[i_eta+1]&&phi_a>=PhiBins[i_phi]&&phi_a<PhiBins[i_phi+1]){
					NchLattice[i_eta][i_phi]++;
					MpTLattice[i_eta][i_phi]+=pt_a;
				}
			}
		}
	}
	// analyzing array
	for(Int_t i_eta=0;i_eta<fnEtaBinsRho;++i_eta){
		for(Int_t i_phi=0;i_phi<fnPhiBinsRho;++i_phi){
			if(NchLattice[i_eta][i_phi]>0)
				MpTLattice[i_eta][i_phi]/=(1.0*NchLattice[i_eta][i_phi]);
			else
				MpTLattice[i_eta][i_phi]=0.0;
		}
	}
	Double_t mNch=0;
	Double_t mMpT=0;
	for(Int_t i_eta=0;i_eta<fnEtaBinsRho;++i_eta){
		for(Int_t i_phi=0;i_phi<fnPhiBinsRho;++i_phi){
			mMpT+=MpTLattice[i_eta][i_phi];
			mNch+=1.0*NchLattice[i_eta][i_phi];
		}
	}
	// average activity per cell
	mMpT/=(1.0*fnEtaBinsRho*fnPhiBinsRho);
	mNch/=(1.0*fnEtaBinsRho*fnPhiBinsRho);
	// get sigma
	Double_t sNch_tmp=0;
	Double_t sMpT_tmp=0;
	for(Int_t i_eta=0;i_eta<fnEtaBinsRho;++i_eta){
		for(Int_t i_phi=0;i_phi<fnPhiBinsRho;++i_phi){
			sMpT_tmp+=TMath::Power(MpTLattice[i_eta][i_phi]-mMpT,2);
			sNch_tmp+=TMath::Power(1.0*NchLattice[i_eta][i_phi]-mNch,2);
		}
	}
	sMpT_tmp/=(1.0*fnEtaBinsRho*fnPhiBinsRho);
	sNch_tmp/=(1.0*fnEtaBinsRho*fnPhiBinsRho);
	Double_t sMpT=TMath::Sqrt(sMpT_tmp);
	Int_t minmult = 21;// twice the average Nch (|eta|<0.8) in MB
	if(fEtaCutRho>3.0)
		minmult = 100;// twice the average Nch (|eta|<4) in MB
	//cout<<"minmult="<<minmult<<endl;
	if( nchtotal > minmult && TMath::Abs(fGenLeadEta) < fEtaCut ){
		Double_t rho = sMpT/mMpT;
		Double_t meanrho = frho->Eval(1.0*nchtotal);
		Int_t indexrho = -1;
		if((rho/meanrho)<0.85)
			indexrho = 1;
		else if((rho/meanrho)>1.15)
			indexrho = 2;
		hPtLeadingRho[0]->Fill(fGenLeadPt);
		hNchRho[0]->Fill(1.0*nchtotal,fGenLeadPt);
		hEtaLeadingRho[0]->Fill(fGenLeadEta);
		if(indexrho>0){
			hPtLeadingRho[indexrho]->Fill(fGenLeadPt);
			hNchRho[indexrho]->Fill(1.0*nchtotal,fGenLeadPt);
			hEtaLeadingRho[indexrho]->Fill(fGenLeadEta);
		}
		for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {

			if(i==fGenLeadIn)
				continue;

			AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
			if (!particle) continue;

			if (!fMC->IsPhysicalPrimary(i)) continue;
			if (particle->Charge() == 0) continue;
			if( particle->Pt() < ptamin)continue;
			Double_t Deta = fGenLeadEta-particle->Eta();
			Double_t DPhi = DeltaPhi(particle->Phi(), fGenLeadPhi);
			hDetaDphiRhoWideEta[0]->Fill(Deta,DPhi);
			if(indexrho>0){
				hDetaDphiRhoWideEta[indexrho]->Fill(Deta,DPhi);
			}
			if ( TMath::Abs(particle->Eta()) > fEtaCutRho )continue;
			Int_t pid = GetPidCode(particle->PdgCode());
			hDetaDphiRho[0]->Fill(Deta,DPhi);
			if(indexrho>0){
				hDetaDphiRho[indexrho]->Fill(Deta,DPhi);
			}
			// definition of the topological regions
			if(TMath::Abs(DPhi)<pi/3.0){// near side
				hNchPtPidRho[0][0]->Fill(1.0*nchtotal,particle->Pt(),pid);
				if(indexrho>0)
					hNchPtPidRho[0][indexrho]->Fill(1.0*nchtotal,particle->Pt(),pid);
			}
			else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
				hNchPtPidRho[1][0]->Fill(1.0*nchtotal,particle->Pt(),pid);
				if(indexrho>0)
					hNchPtPidRho[1][indexrho]->Fill(1.0*nchtotal,particle->Pt(),pid);
			}
			else{// transverse side
				hNchPtPidRho[2][0]->Fill(1.0*nchtotal,particle->Pt(),pid);
				if(indexrho>0)
					hNchPtPidRho[2][indexrho]->Fill(1.0*nchtotal,particle->Pt(),pid);
				if(TMath::Abs(Deta)>2.0){
					hNchPtPidRho[3][0]->Fill(1.0*nchtotal,particle->Pt(),pid);
					if(indexrho>0)
						hNchPtPidRho[3][indexrho]->Fill(1.0*nchtotal,particle->Pt(),pid);
				}
			}

		}

	}

}

//_______________________________________________________
void AliAnalysisTaskGenMcKnoUe::MakeALICE3Analysis(){

	Double_t EtaBins[fnEtaBinsRho+1];
	Double_t deltaEta = (2.0*fEtaCutRho)/(1.0*fnEtaBinsRho);
	for(Int_t i_eta=0;i_eta<fnEtaBinsRho+1;++i_eta){
		EtaBins[i_eta]=0;
		if(i_eta<fnEtaBinsRho)
			EtaBins[i_eta]=i_eta*deltaEta - 1.0*fEtaCutRho;
		else
			EtaBins[i_eta]= 1.0*fEtaCutRho;
	}
	Double_t PhiBins[fnPhiBinsRho+1];
	Double_t deltaPhi = (2.0*pi)/(1.0*fnPhiBinsRho);
	for(Int_t i_phi=0;i_phi<fnPhiBinsRho+1;++i_phi){
		PhiBins[i_phi]=0;
		if(i_phi<fnPhiBinsRho)
			PhiBins[i_phi]=i_phi*deltaPhi;
		else
			PhiBins[i_phi]= 2.0*pi;
	}
	Int_t NchLattice[fnEtaBinsRho][fnPhiBinsRho];
	for(Int_t i_eta=0;i_eta<fnEtaBinsRho;++i_eta)
		for(Int_t i_phi=0;i_phi<fnPhiBinsRho;++i_phi)
			NchLattice[i_eta][i_phi]=0;
	Double_t MpTLattice[fnEtaBinsRho][fnPhiBinsRho];
	for(Int_t i_eta=0;i_eta<fnEtaBinsRho;++i_eta)
		for(Int_t i_phi=0;i_phi<fnPhiBinsRho;++i_phi)
			MpTLattice[i_eta][i_phi]=0;

	Double_t totalpt=0;
	Int_t totalnch_forpt=0;
	Int_t nchtotal=0;

	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {

		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;
		Float_t eta_a = particle->Eta();
		Float_t phi_a = particle->Phi();
		Float_t pt_a  = particle->Pt();

		if (!fMC->IsPhysicalPrimary(i)) continue;
		if (particle->Charge() == 0) continue;
		if ( TMath::Abs(eta_a) > fEtaCutRho )continue;
		if ( pt_a <= 0 )continue;
		nchtotal++;
		if(pt_a<ptamin||pt_a>ptamax)continue;
		totalpt+=pt_a;
		totalnch_forpt++;
		// loop over all eta and phi intervals
		for(Int_t i_eta=0;i_eta<fnEtaBinsRho;++i_eta){
			for(Int_t i_phi=0;i_phi<fnPhiBinsRho;++i_phi){
				if(eta_a>=EtaBins[i_eta]&&eta_a<EtaBins[i_eta+1]&&phi_a>=PhiBins[i_phi]&&phi_a<PhiBins[i_phi+1]){
					NchLattice[i_eta][i_phi]++;
					MpTLattice[i_eta][i_phi]+=pt_a;
				}
			}
		}
	}
	// analyzing array
	for(Int_t i_eta=0;i_eta<fnEtaBinsRho;++i_eta){
		for(Int_t i_phi=0;i_phi<fnPhiBinsRho;++i_phi){
			if(NchLattice[i_eta][i_phi]>0)
				MpTLattice[i_eta][i_phi]/=(1.0*NchLattice[i_eta][i_phi]);
			else
				MpTLattice[i_eta][i_phi]=0.0;
		}
	}
	Double_t mNch=0;
	Double_t mMpT=0;
	for(Int_t i_eta=0;i_eta<fnEtaBinsRho;++i_eta){
		for(Int_t i_phi=0;i_phi<fnPhiBinsRho;++i_phi){
			mMpT+=MpTLattice[i_eta][i_phi];
			mNch+=1.0*NchLattice[i_eta][i_phi];
		}
	}
	// average activity per cell
	mMpT/=(1.0*fnEtaBinsRho*fnPhiBinsRho);
	mNch/=(1.0*fnEtaBinsRho*fnPhiBinsRho);
	// get sigma
	Double_t sNch_tmp=0;
	Double_t sMpT_tmp=0;
	for(Int_t i_eta=0;i_eta<fnEtaBinsRho;++i_eta){
		for(Int_t i_phi=0;i_phi<fnPhiBinsRho;++i_phi){
			sMpT_tmp+=TMath::Power(MpTLattice[i_eta][i_phi]-mMpT,2);
			sNch_tmp+=TMath::Power(1.0*NchLattice[i_eta][i_phi]-mNch,2);
		}
	}
	sMpT_tmp/=(1.0*fnEtaBinsRho*fnPhiBinsRho);
	sNch_tmp/=(1.0*fnEtaBinsRho*fnPhiBinsRho);
	Double_t sMpT=TMath::Sqrt(sMpT_tmp);
	hnchmpi->Fill(nchtotal,Nmpi);
	if(mMpT>0){
		hnchmpirho->Fill(1.0*nchtotal,1.0*Nmpi,sMpT/mMpT);
		hnchrho->Fill(1.0*nchtotal,sMpT/mMpT);
		hmpirho->Fill(1.0*Nmpi,sMpT/mMpT);
	}
}

//____________________________________________________________

Double_t AliAnalysisTaskGenMcKnoUe::DeltaPhi(Double_t phia, Double_t phib,
		Double_t rangeMin, Double_t rangeMax)
{
	Double_t dphi = -999;
	Double_t pi = TMath::Pi();

	if (phia < 0)         phia += 2*pi;
	else if (phia > 2*pi) phia -= 2*pi;
	if (phib < 0)         phib += 2*pi;
	else if (phib > 2*pi) phib -= 2*pi;
	dphi = phib - phia;
	if (dphi < rangeMin)      dphi += 2*pi;
	else if (dphi > rangeMax) dphi -= 2*pi;

	return dphi;
}
void AliAnalysisTaskGenMcKnoUe::SetParametersRho(Double_t etarange){

	if(etarange==0.8){
		//cout<<"seting the case eta 0.8"<<endl;
		switch(fGenerator){
			case 0: // Monash
				frho->SetParameter(0,-1.2796503);
				frho->SetParameter(1,-0.21406574);
				frho->SetParameter(2,7.9576486);
				frho->SetParameter(3,-0.50547998);
				break;
			case 1: // Monash NoCR
				frho->SetParameter(0,-0.62615697);
				frho->SetParameter(1,-0.12906934);
				frho->SetParameter(2,9.4532345);
				frho->SetParameter(3,-0.57046898);
				break;
			case 2: // Monash Ropes
				frho->SetParameter(0,-1.3036594);
				frho->SetParameter(1,-0.18819465);
				frho->SetParameter(2,8.2080068);
				frho->SetParameter(3,-0.51231662);
				break;
			case 3: // Epos LHC
				frho->SetParameter(0,-0.89774856);
				frho->SetParameter(1,-0.15704294);
				frho->SetParameter(2,8.7380225);
				frho->SetParameter(3,-0.52647214);
				break;
			case 4: // Herwig
				frho->SetParameter(0,-1.5178880);
				frho->SetParameter(1,-0.094629819);
				frho->SetParameter(2,17.327223);
				frho->SetParameter(3,-0.75409746);
				break;
			case 5: // AMPT (to be checked)
				frho->SetParameter(0,-1.8832819);
				frho->SetParameter(1,-0.18063863);
				frho->SetParameter(2,10.062015);
				frho->SetParameter(3,-0.57335253);
				break;
			case 6: // AMPT no string melting
				frho->SetParameter(0,-1.8299821);
				frho->SetParameter(1,-0.17885455);
				frho->SetParameter(2,10.072160);
				frho->SetParameter(3,-0.57360867);
				break;
			default:
				frho->SetParameter(0,-2034448.6);
				frho->SetParameter(1,-0.99051656);
				frho->SetParameter(2,14.116976);
				frho->SetParameter(3,-0.52448073);
		}
	}
	else{ // default eta < 4
		switch(fGenerator){
			case 0:
				frho->SetParameter(0,-2034448.6);
				frho->SetParameter(1,-0.99051656);
				frho->SetParameter(2,14.116976);
				frho->SetParameter(3,-0.52448073);
				break;
			case 1:
				frho->SetParameter(0,-7.8071100);
				frho->SetParameter(1,-1.9595275);
				frho->SetParameter(2,18.226513);
				frho->SetParameter(3,-0.58994765);
				break;
			case 2:
				frho->SetParameter(0,-7.8071100);
				frho->SetParameter(1,-1.1168320);
				frho->SetParameter(2,14.480357);
				frho->SetParameter(3,-0.52709865);
				break;
			case 3:
				frho->SetParameter(0,-1891978.1);
				frho->SetParameter(1,-1.8212843);
				frho->SetParameter(2,13.573444);
				frho->SetParameter(3,-0.51139756);
				break;
			case 4:
				frho->SetParameter(0,-7.8071100);
				frho->SetParameter(1,-2.0788016);
				frho->SetParameter(2,28.243893);
				frho->SetParameter(3,-0.70047215);
				break;
			case 5: // AMPT string melting (to be checked)
				frho->SetParameter(0,-7.8071100);
				frho->SetParameter(1,-1.0183988);
				frho->SetParameter(2,21.348248);
				frho->SetParameter(3,-0.61011648);
				break;
			case 6: // AMPT no string melting (to be checked)
				frho->SetParameter(0,-7.8071100);
				frho->SetParameter(1,-1.0183988);
				frho->SetParameter(2,21.348248);
				frho->SetParameter(3,-0.61011648);
				break;
			default:
				frho->SetParameter(0,-2034448.6);
				frho->SetParameter(1,-0.99051656);
				frho->SetParameter(2,14.116976);
				frho->SetParameter(3,-0.52448073);
		}
	}
}
//_______________________________________________________
Int_t AliAnalysisTaskGenMcKnoUe::GetPidCode(Int_t pdgCode)  {

	Int_t pidCode = 3;

	switch (TMath::Abs(pdgCode)) {
		case 211:
			pidCode = 0; // pion
			break;
		case 321:
			pidCode = 1; // kaon
			break;
		case 2212:
			pidCode = 2; // proton
			break;
			//default:
			//	3;
	};

	return pidCode;
}
//______________________________________________________________________________
void AliAnalysisTaskGenMcKnoUe::Terminate(Option_t *)
{
	fOutputList = dynamic_cast<TList*> (GetOutputData(1));
	if (!fOutputList) { Printf("ERROR: Output list not available"); return; }

	return;
}

