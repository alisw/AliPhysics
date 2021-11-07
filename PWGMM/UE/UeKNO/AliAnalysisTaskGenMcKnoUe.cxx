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
 *         Last modification: 03/11/2021                                  * 
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
#include "TH3.h"
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
// #include "AliAnalysisUtils.h"

//_____ AnalysisTask headers
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskGenMcKnoUe.h"
//_____ STL includes
#include <iostream>
using namespace std;

//using std::cout;
//using std::endl;
//using namespace std;

const Char_t * NameOfRegion[3]={"NS","AS","TS"};
const Int_t NchNBins = 2000;
const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;
const Double_t  ptamin = 0.15;
const Double_t  ptamax = 20.0;
const Double_t  etamax = 4.0;
const Int_t nEtaBins = 10;
const Int_t nPhiBins = 10;
Double_t Nmpi = -1;

ClassImp(AliAnalysisTaskGenMcKnoUe)

AliAnalysisTaskGenMcKnoUe::AliAnalysisTaskGenMcKnoUe() : AliAnalysisTaskSE(),
	fMC(0x0),fMcHandler(0x0),fMCStack(0x0),fEtaCut(0.8),fIsPP(kTRUE),fPtMin(0.5),fGenLeadPhi(0x0),fGenLeadPt(0x0),fGenLeadIn(0x0),hPtLeadingTrue(0x0),hPtLVsV0A(0x0),hetaphi(0x0),hnchmpi(0x0),hnchmpirho(0x0),hnchrho(0x0),hmpirho(0x0),fOutputList(0)
{
	for(Int_t i=0;i<3;++i) 
		hPtVsPtLeadingTrue[i]=0;

}
//_____________________________________________________________________________
AliAnalysisTaskGenMcKnoUe::AliAnalysisTaskGenMcKnoUe(const char* name) : AliAnalysisTaskSE(name),
	fMC(0x0),fMcHandler(0x0),fMCStack(0x0),fEtaCut(0.8),fIsPP(kTRUE),fPtMin(0.5),fGenLeadPhi(0x0),fGenLeadPt(0x0),fGenLeadIn(0x0),hPtLeadingTrue(0x0),hPtLVsV0A(0x0),hetaphi(0x0),hnchmpi(0x0),hnchmpirho(0x0),hnchrho(0x0),hmpirho(0x0),fOutputList(0)
{
	for(Int_t i=0;i<3;++i)
		hPtVsPtLeadingTrue[i]=0;

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
	Double_t EtaBins[nEtaBins+1];
	Double_t deltaEta = (2.0*etamax)/(1.0*nEtaBins);
	for(Int_t i_eta=0;i_eta<nEtaBins+1;++i_eta){
		EtaBins[i_eta]=0;
		if(i_eta<nEtaBins)
			EtaBins[i_eta]=i_eta*deltaEta - 1.0*etamax;
		else
			EtaBins[i_eta]= 1.0*etamax;
	}
	Double_t PhiBins[nPhiBins+1];
	Double_t deltaPhi = (2.0*pi)/(1.0*nPhiBins);
	for(Int_t i_phi=0;i_phi<nPhiBins+1;++i_phi){
		PhiBins[i_phi]=0;
		if(i_phi<nPhiBins)
			PhiBins[i_phi]=i_phi*deltaPhi - 1.0*pi;
		else
			PhiBins[i_phi]= 1.0*pi;
	}

	//OpenFile(1);
	fOutputList = new TList();
	fOutputList->SetOwner(kTRUE);

	hPtLVsV0A = 0;
	hPtLVsV0A = new TH2D("hPtLVsV0A","",pTNBinsL,pTNBins1L,NchNBins,NchBins);
	fOutputList->Add(hPtLVsV0A);

	hPtLeadingTrue = new TH1D("hPtLeadingTrue","",pTNBinsL,pTNBins1L);
	fOutputList->Add(hPtLeadingTrue);

	hetaphi = 0;
	hetaphi = new TH2D("hetaphi","; #eta; #phi (rad)",nEtaBins,EtaBins,nPhiBins,PhiBins);
	fOutputList->Add(hetaphi);

	hnchmpi = 0;
	hnchmpi = new TH2D("hnchmpi","; #it{N}_{ch}; #it{N}_{mpi}",600,-0.5,599.5,600,-0.5,599.5);
	fOutputList->Add(hnchmpi);


	hnchmpirho = 0;
	hnchmpirho = new TH3D("hnchmpirho","; #it{N}_{ch}; #it{N}_{mpi}; #rho",600,-0.5,599.5,600,-0.5,599.5,1000,0.0,5.0);
	fOutputList->Add(hnchmpirho);

	hnchrho = 0;
	hnchrho = new TH2D("hnchrho","; #it{N}_{ch}; #rho",600,-0.5,599.5,1000,0,5.0);
	fOutputList->Add(hnchrho);

	hmpirho = 0;
	hmpirho = new TH2D("hmpirho","; #it{N}_{mpi}; #rho",600,-0.5,599.5,1000,0,5.0);
	fOutputList->Add(hmpirho);

	// UE analysis
	for(Int_t i=0;i<3;++i){

		hPtVsPtLeadingTrue[i] = new TH3D(Form("hPtVsPtLeadingTrue_%s",NameOfRegion[i]),"",pTNBinsL,pTNBins1L,pTNBins,pTNBins1,NchNBins,NchBins);
		fOutputList->Add(hPtVsPtLeadingTrue[i]);

	}

	PostData(1, fOutputList);

}
//_____________________________________________________________________________
void AliAnalysisTaskGenMcKnoUe::UserExec(Option_t *)
{

	// ### Initialize

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

	if(fIsPP)
		MakeALICE3Analysis();

	if(isGoodVtxPosMC)
		if(fGenLeadPt>=fPtMin)
			GetGenUEObservables();

	PostData(1, fOutputList);
	return;
}

//______________________________________________________________________________
void AliAnalysisTaskGenMcKnoUe::Terminate(Option_t *)
{

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

	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {

		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;
		if (!fMC->IsPhysicalPrimary(i)) continue;
		if (particle->Charge() == 0) continue;
		if ( TMath::Abs(particle->Eta()) > fEtaCut )continue;
		if( particle->Pt() < fPtMin)continue;

		if (flPt<particle->Pt()){
			flPt = particle->Pt();
			flPhi = particle->Phi();
			flIndex = i;

		}
	}

	fGenLeadPhi = flPhi;
	fGenLeadPt  = flPt;
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
		}

	}
	hPtLeadingTrue->Fill(fGenLeadPt);
}
//_______________________________________________________
void AliAnalysisTaskGenMcKnoUe::MakeALICE3Analysis(){

	Double_t EtaBins[nEtaBins+1];
	Double_t deltaEta = (2.0*etamax)/(1.0*nEtaBins);
	for(Int_t i_eta=0;i_eta<nEtaBins+1;++i_eta){
		EtaBins[i_eta]=0;
		if(i_eta<nEtaBins)
			EtaBins[i_eta]=i_eta*deltaEta - 1.0*etamax;
		else
			EtaBins[i_eta]= 1.0*etamax;
	}
	Double_t PhiBins[nPhiBins+1];
	Double_t deltaPhi = (2.0*pi)/(1.0*nPhiBins);
	for(Int_t i_phi=0;i_phi<nPhiBins+1;++i_phi){
		PhiBins[i_phi]=0;
		if(i_phi<nPhiBins)
			PhiBins[i_phi]=i_phi*deltaPhi - 1.0*pi;
		else
			PhiBins[i_phi]= 1.0*pi;
	}

	Int_t NchLattice[nEtaBins][nPhiBins];
	for(Int_t i_eta=0;i_eta<nEtaBins;++i_eta)
		for(Int_t i_phi=0;i_phi<nPhiBins;++i_phi)
			NchLattice[i_eta][i_phi]=0;
	Double_t MpTLattice[nEtaBins][nPhiBins];
	for(Int_t i_eta=0;i_eta<nEtaBins;++i_eta)
		for(Int_t i_phi=0;i_phi<nPhiBins;++i_phi)
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
		if ( TMath::Abs(eta_a) > etamax )continue;
		if ( pt_a <= 0 )continue;
		nchtotal++;
		if(pt_a<ptamin||pt_a>ptamax)continue;
		totalpt+=pt_a;
		totalnch_forpt++;
		// loop over all eta and phi intervals
		for(Int_t i_eta=0;i_eta<nEtaBins;++i_eta){
			for(Int_t i_phi=0;i_phi<nPhiBins;++i_phi){
				if(eta_a>=EtaBins[i_eta]&&eta_a<EtaBins[i_eta+1]&&phi_a>=PhiBins[i_phi]&&phi_a<PhiBins[i_phi+1]){
					NchLattice[i_eta][i_phi]++;
					MpTLattice[i_eta][i_phi]+=pt_a;
				}
			}
		}
		hetaphi->Fill(eta_a,phi_a);
	}
	// analyzing array
	for(Int_t i_eta=0;i_eta<nEtaBins;++i_eta){
		for(Int_t i_phi=0;i_phi<nPhiBins;++i_phi){
			if(NchLattice[i_eta][i_phi]>0)
				MpTLattice[i_eta][i_phi]/=(1.0*NchLattice[i_eta][i_phi]);
			else
				MpTLattice[i_eta][i_phi]=0.0;
		}
	}
	Double_t mNch=0;
	Double_t mMpT=0;
	for(Int_t i_eta=0;i_eta<nEtaBins;++i_eta){
		for(Int_t i_phi=0;i_phi<nPhiBins;++i_phi){
			mMpT+=MpTLattice[i_eta][i_phi];
			mNch+=1.0*NchLattice[i_eta][i_phi];
		}
	}
	// average activity per cell
	mMpT/=(1.0*nEtaBins*nPhiBins);
	mNch/=(1.0*nEtaBins*nPhiBins);
	// get sigma
	Double_t sNch_tmp=0;
	Double_t sMpT_tmp=0;
	for(Int_t i_eta=0;i_eta<nEtaBins;++i_eta){
		for(Int_t i_phi=0;i_phi<nPhiBins;++i_phi){
			sMpT_tmp+=TMath::Power(MpTLattice[i_eta][i_phi]-mMpT,2);
			sNch_tmp+=TMath::Power(1.0*NchLattice[i_eta][i_phi]-mNch,2);
		}
	}
	sMpT_tmp/=(1.0*nEtaBins*nPhiBins);
	sNch_tmp/=(1.0*nEtaBins*nPhiBins);
	//Double_t sNch=TMath::Sqrt(sNch_tmp);
	Double_t sMpT=TMath::Sqrt(sMpT_tmp);

	hnchmpi->Fill(nchtotal,Nmpi);
	hnchmpirho->Fill(nchtotal,Nmpi,sMpT/mMpT);
	hnchrho->Fill(nchtotal,sMpT/mMpT);
	hmpirho->Fill(Nmpi,sMpT/mMpT);

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

