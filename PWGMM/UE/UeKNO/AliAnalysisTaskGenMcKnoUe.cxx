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

#include <Riostream.h>
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TList.h"

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVTrack.h"
#include "AliAnalysisTask.h"
#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include <AliHeader.h>
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliEventCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include <TTree.h>
#include <TMath.h>
#include <TRandom.h>
#include <TDirectory.h>
#include <TBits.h>
#include <AliAnalysisFilter.h>
#include "AliAnalysisManager.h"
using std::cout;
using std::endl;

#include "AliAnalysisTaskGenMcKnoUe.h"


const Char_t * NameOfRegion[3]={"NS","AS","TS"};
const Int_t NchNBins = 2000;

const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;
class AliAnalysisTaskGenMcKnoUe;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskGenMcKnoUe) // classimp: necessary for root

	AliAnalysisTaskGenMcKnoUe::AliAnalysisTaskGenMcKnoUe() : AliAnalysisTaskSE(),

	fMC(0x0),fMcHandler(0x0),fMCStack(0),fEtaCut(0.8),fPtMin(0.5),fOutputList(0),fGenLeadPhi(0),fGenLeadPt(0),fGenLeadIn(0), hPtLeadingTrue(0),hPtLVsV0A(0)

{
	for(Int_t i=0;i<3;++i) 
		hPtVsPtLeadingTrue[i]=0;

}
//_____________________________________________________________________________
AliAnalysisTaskGenMcKnoUe::AliAnalysisTaskGenMcKnoUe(const char* name) : AliAnalysisTaskSE(name),

	fMC(0x0),fMcHandler(0x0),fMCStack(0),fEtaCut(0.8),fPtMin(0.5),fOutputList(0),fGenLeadPhi(0),fGenLeadPt(0),fGenLeadIn(0), hPtLeadingTrue(0),hPtLVsV0A(0)

{
	for(Int_t i=0;i<3;++i)
		hPtVsPtLeadingTrue[i]=0;

	// constructor
	DefineInput(0, TChain::Class());    // define the input of the analysis: in this case you take a 'chain' of events
	// this chain is created by the analysis manager, so no need to worry about it, does its work automatically
	DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms

}
//_____________________________________________________________________________
AliAnalysisTaskGenMcKnoUe::~AliAnalysisTaskGenMcKnoUe()
{
	// destructor
	if(fOutputList) {
		delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
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

	OpenFile(1);
	fOutputList = new TList();          // this is a list which will contain all of your histograms
	// at the end of the analysis, the contents of this list are written  to the output file
	fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects and will delete them if requested

	hPtLVsV0A = 0;
	hPtLVsV0A = new TH2D("hPtLVsV0A","",pTNBinsL,pTNBins1L,NchNBins,NchBins);
	fOutputList->Add(hPtLVsV0A);

	hPtLeadingTrue = new TH1D("hPtLeadingTrue","",pTNBinsL,pTNBins1L);
	fOutputList->Add(hPtLeadingTrue);

	// UE analysis
	for(Int_t i=0;i<3;++i){


		hPtVsPtLeadingTrue[i] = new TH3D(Form("hPtVsPtLeadingTrue_%s",NameOfRegion[i]),"",pTNBinsL,pTNBins1L,pTNBins,pTNBins1,NchNBins,NchBins);
		fOutputList->Add(hPtVsPtLeadingTrue[i]);

	}

	PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the

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
		cout << "Name of the file with pb :" << endl; //fInputHandler->GetTree()->GetCurrentFile()->GetName() << 
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

	// Before trigger selection

	GetGenLeadingObject();// leading particle at gen level
	
	if(isGoodVtxPosMC){
		// UE analysis
		if(fGenLeadPt>=fPtMin){
			GetGenUEObservables();
		}

	}

	PostData(1, fOutputList); // stream the result of this event to the output manager which will write it to a file

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
		if (!fMC->IsPhysicalPrimary(i)) continue;  //  Particles produced including products of strong and electromagnetic decays but excluding feed-down from weak decays of strange particles like Ks,Lambda etc)
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

