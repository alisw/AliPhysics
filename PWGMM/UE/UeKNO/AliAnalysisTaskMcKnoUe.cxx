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
 * Author: Prabi, Feng Fan, and Paolo Bartalini, CCNU 2016               *
 **************************************************************************/

/* AliAnaysisTaskMcKnoUe source code
 *
 * simple task which can serve as a starting point for building an MPI UE analysis
 *
 */

class TTree;

class AliPPVsMultUtils;
class AliESDtrackCuts;


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
//#include "AliV0vertexer.h"


#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDVZERO.h"

#include "AliESDInputHandler.h"
#include "AliMultiplicity.h"
#include "AliPPVsMultUtils.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskESDfilter.h"

#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"

#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include <AliHeader.h>

#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultVariable.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"

#include "AliESDUtils.h"

#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
//#include "AliGenPythiaEventHeader.h"

#include "AliEventCuts.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliPPVsMultUtils.h"
#include <AliESDVertex.h>
#include <AliMultiplicity.h>
#include <TTree.h>
#include <TMath.h>
#include <TRandom.h>
#include <TDirectory.h>
#include <TBits.h>
#include <AliAnalysisFilter.h>

using std::cout;
using std::endl;

#include "AliAnalysisTaskMcKnoUe.h"
const Int_t ptNbins = 36;
Double_t ptbins1[ptNbins+1] = {
	0.0,  0.1,  0.15,  0.2,  0.25,  0.3,  0.35,  0.4,  0.45,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5, 5.0, 6.0,    7.0,  8.0,  9.0,  10.0,  12.0,  14.0,  16.0,  18.0,  20.0,  25.0,  30.0,  40.0,  50.0
};
const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;
class AliAnalysisTaskMcKnoUe;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskMcKnoUe) // classimp: necessary for root

AliAnalysisTaskMcKnoUe::AliAnalysisTaskMcKnoUe() : AliAnalysisTaskSE(),
	fESD(0), fEventCuts(0x0), fStack(0), fMC(0), fUseMC(kTRUE), fLeadingTrackFilter(0x0), fTrackFilter(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5), fLeadPtCutMin(5.0), fLeadPtCutMax(40.0), fGenLeadPhi(0), fGenLeadPt(0), fGenLeadIn(0), fRecLeadPhi(0), fRecLeadPt(0), fRecLeadIn(0), hNchTSGen(0), hNchTSGenTest(0), hNchTSRec(0), hNchTSRecTest(0), hNchResponse(0)  

{
	for(Int_t i=0;i<3;++i){ 
		hPhiGen[i]= 0;
		hPhiRec[i]= 0;
	}
	// default constructor, don't allocate memory here!  this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskMcKnoUe::AliAnalysisTaskMcKnoUe(const char* name) : AliAnalysisTaskSE(name),
	fESD(0), fEventCuts(0x0), fStack(0), fMC(0), fUseMC(kTRUE), fLeadingTrackFilter(0x0), fTrackFilter(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5), fLeadPtCutMin(5.0), fLeadPtCutMax(40.0), fGenLeadPhi(0), fGenLeadPt(0), fGenLeadIn(0), fRecLeadPhi(0), fRecLeadPt(0), fRecLeadIn(0), hNchTSGen(0), hNchTSGenTest(0), hNchTSRec(0), hNchTSRecTest(0), hNchResponse(0)
{
	for(Int_t i=0;i<3;++i){
		hPhiGen[i] = 0;
		hPhiRec[i] = 0;
	}
	// constructor
	DefineInput(0, TChain::Class());    // define the input of the analysis: in this case you take a 'chain' of events
	// this chain is created by the analysis manager, so no need to worry about it, does its work automatically
	DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms

}
//_____________________________________________________________________________
AliAnalysisTaskMcKnoUe::~AliAnalysisTaskMcKnoUe()
{
	// destructor
	if(fOutputList) {
		delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
		fOutputList = 0x0;
	}

}
//_____________________________________________________________________________
void AliAnalysisTaskMcKnoUe::UserCreateOutputObjects()
{
	// fCuts *** leading particle ***
	if(!fLeadingTrackFilter){
		fLeadingTrackFilter = new AliAnalysisFilter("trackFilter2015");
		AliESDtrackCuts * fCuts1 = new AliESDtrackCuts();
		fCuts1->SetMaxFractionSharedTPCClusters(0.4);//
		fCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);//
		fCuts1->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);//
		fCuts1->SetMaxChi2PerClusterTPC(4);//
		fCuts1->SetAcceptKinkDaughters(kFALSE);//
		fCuts1->SetRequireTPCRefit(kTRUE);//
		fCuts1->SetRequireITSRefit(kTRUE);//
		fCuts1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
				AliESDtrackCuts::kAny);//
		fCuts1->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");//
		fCuts1->SetMaxChi2TPCConstrainedGlobal(36);//
		fCuts1->SetMaxDCAToVertexZ(2);//
		fCuts1->SetDCAToVertex2D(kFALSE);//
		fCuts1->SetRequireSigmaToVertex(kFALSE);//
		fCuts1->SetMaxChi2PerClusterITS(36);//
		fLeadingTrackFilter->AddCuts(fCuts1);
	}

	///track quality =====
	// TPC ***  multiplicity in transverse side ***
	if(!fTrackFilter){
		fTrackFilter = new AliAnalysisFilter("trackFilterTPConly");
		AliESDtrackCuts * fCuts2 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
		fCuts2->SetRequireTPCRefit(kTRUE);
		fCuts2->SetRequireITSRefit(kTRUE);
		fCuts2->SetEtaRange(-0.8,0.8);
		fTrackFilter->AddCuts(fCuts2);
	}

	// create output objects

	OpenFile(1);
	fOutputList = new TList();          // this is a list which will contain all of your histograms
	// at the end of the analysis, the contents of this list are written  to the output file
	fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects and will delete them if requested


	const Char_t * nameReg[3]={"NS","AS","TS"};
	hNchTSGen = new TH1D("hNchTSGen","",100,-0.5,99.5);
	fOutputList->Add(hNchTSGen);
	for(Int_t i=0;i<3;++i){
		hPhiGen[i]= new TH1D(Form("hPhiGen_%s",nameReg[i]),"",64,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
		fOutputList->Add(hPhiGen[i]);
	}

	hNchTSRec = new TH1D("hNchTSRec","",100,-0.5,99.5);
	fOutputList->Add(hNchTSRec);
	for(Int_t i=0;i<3;++i){
		hPhiRec[i]= new TH1D(Form("hPhiRec_%s",nameReg[i]),"",64,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
		fOutputList->Add(hPhiRec[i]);
	}
	hNchResponse = new TH2D("hNchResponse","Detector response; rec mult; gen mult",100,-0.5,99.5,100,-0.5,99.5);
	fOutputList->Add(hNchResponse);


	fEventCuts.AddQAplotsToList(fOutputList);
	PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the

}
//_____________________________________________________________________________
void AliAnalysisTaskMcKnoUe::UserExec(Option_t *)
{

	AliVEvent *event = InputEvent();
	if (!event) {
		Error("UserExec", "Could not retrieve event");
		return;
	}

	fESD = dynamic_cast<AliESDEvent*>(event);

	if(!fESD){
		Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
		this->Dump();
		return;
	}

	if (fUseMC) {

		//      E S D
		fMC = dynamic_cast<AliMCEvent*>(MCEvent());
		if(!fMC){
			Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}

		/*
		   fMCStack = fMC->Stack();

		   if(!fMCStack){
		   Printf("%s:%d MCStack not found in Input Manager",(char*)__FILE__,__LINE__);
		   this->Dump();
		   return;
		   }
		 */
	}


	//AliMultSelection *MultSelection = (AliMultSelection*) fESD -> FindListObject("MultSelection");
	// Cuts at event level
	UInt_t fSelectMask= fInputHandler->IsEventSelected();
	Bool_t isINT7selected = fSelectMask&AliVEvent::kINT7;
	if(!isINT7selected)
		return;


/*
	if (!fEventCuts.AcceptEvent(event)) {
		PostData(1, fListOfObjects);
		return;
	}
*/


	//if(!IsEventSelected(fESD) ) {
//		PostData(1, fOutputList);   /// Event isn't selected, post output data, done here
	//	return;
//	}


	if (fUseMC){
		GetLeadingObject(kTRUE);// leading particle at gen level
	}

	GetLeadingObject(kFALSE);// leading particle at rec level

	// detector response (50% of full stat)
	if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax )){

		Double_t randomUE = gRandom->Uniform(0.0,1.0);
		if(randomUE<0.5)
			GetDetectorResponse();
		else{


		}

	}


	PostData(1, fOutputList); // stream the result of this event to the output manager which will write it to a file

}


//______________________________________________________________________________
void AliAnalysisTaskMcKnoUe::Terminate(Option_t *)
{

}


//_____________________________________________________________________________
void AliAnalysisTaskMcKnoUe::GetLeadingObject(Bool_t isMC) {

	Double_t flPt = 0;// leading pT
	Double_t flPhi = 0;
	Int_t flIndex = 0;

	if(isMC){
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
	else{

		Int_t iTracks(fESD->GetNumberOfTracks());           // see how many tracks there are in the event
		for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks

			AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

			if(!track) continue;

			if(!fLeadingTrackFilter->IsSelected(track))
				continue;

			if(TMath::Abs(track->Eta()) > fEtaCut)
				continue;

			if( track->Pt() < fPtMin)continue;

			if (flPt<track->Pt()){
				flPt  = track->Pt();
				flPhi = track->Phi();
				flIndex = i;
			}

		} 
		fRecLeadPhi = flPhi;
		fRecLeadPt  = flPt;
		fRecLeadIn  = flIndex;
	}


}

void AliAnalysisTaskMcKnoUe::GetDetectorResponse() {

	Int_t multTSgen=0;
	Int_t multTSrec=0;

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
			hPhiGen[0]->Fill(DPhi);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			hPhiGen[1]->Fill(DPhi);
		}
		else{// transverse side
			multTSgen++;
			hPhiGen[2]->Fill(DPhi);
		}


	}
	hNchTSGen->Fill(multTSgen);

	Int_t iTracks(fESD->GetNumberOfTracks());           // see how many tracks there are in the event
	for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks

		if(i==fRecLeadIn)
			continue;

		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

		if(!track) continue;

		if(!fTrackFilter->IsSelected(track))
			continue;

		if(TMath::Abs(track->Eta()) > fEtaCut)
			continue;

		if( track->Pt() < fPtMin)continue;

		Double_t DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		// definition of the topological regions
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			hPhiRec[0]->Fill(DPhi);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			hPhiRec[1]->Fill(DPhi);
		}
		else{// transverse side
			multTSrec++;
			hPhiRec[2]->Fill(DPhi);
		}

	}
	hNchTSRec->Fill(multTSrec); 

	hNchResponse->Fill(multTSrec,multTSgen);


}
void AliAnalysisTaskMcKnoUe::GetMultiplicityDistributions(){

	Int_t multTSgen=0;
	Int_t multTSrec=0;

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
			continue;
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			continue;
		}
		else{// transverse side
			multTSgen++;
		}


	}
	hNchTSGenTest->Fill(multTSgen);

	Int_t iTracks(fESD->GetNumberOfTracks());           // see how many tracks there are in the event
	for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks

		if(i==fRecLeadIn)
			continue;

		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

		if(!track) continue;

		if(!fTrackFilter->IsSelected(track))
			continue;

		if(TMath::Abs(track->Eta()) > fEtaCut)
			continue;

		if( track->Pt() < fPtMin)continue;

		Double_t DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		// definition of the topological regions
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			continue;
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			continue;
		}
		else{// transverse side
			multTSrec++;
		}

	}
	hNchTSRecTest->Fill(multTSrec); 


}

Double_t AliAnalysisTaskMcKnoUe::DeltaPhi(Double_t phia, Double_t phib,
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

