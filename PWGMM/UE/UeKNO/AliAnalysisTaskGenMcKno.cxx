/**************************************************************************
 * Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
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
/*  
 *  
 *  AliAnalysisTaskGenMcKno.cxx
 *  
 *
 *  Author: Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)
 *  First version:      April 23, 2020
 * 
 */

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
#include "AliAnalysisTaskGenMcKno.h"

//_____ STL includes
#include <iostream>
using namespace std;

const Int_t nTSBins=700;
const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;
const Char_t * estNames[5] = {"TS","Eta08","V0M","V0A","Eta10"};
const Char_t * PidNames[12] = { "Pion", "Kaon", "Proton", "K0Short", "Lambda", "Xi", "Omega", "Phi", "KStar", "KStarPM", "SigmaZero", "Charged" };
Bool_t iSprimary[11] = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE };
const Int_t npTBins = 36;
Double_t pTBins[npTBins+1] = {
	0.0,  0.1,  0.15,  0.2,  0.25,  0.3 ,  0.35,  0.4,  0.45,  0.5,  
	0.6,  0.7,  0.8 ,  0.9,  1.0 ,  1.25,  1.5 ,  2.0,  2.5 ,  3.0,
	3.5,  4.0,  4.5,   5.0,  6.0 ,  7.0 ,  8.0 ,  9.0,  10.0,  12.0,
	14.0,  16.0,  18.0,  20.0,  30.0,  40.0,  50.0
};


ClassImp( AliAnalysisTaskGenMcKno )

	//_____________________________________________________________________________

	AliAnalysisTaskGenMcKno::AliAnalysisTaskGenMcKno():
		AliAnalysisTaskSE(),
		fMcEvent(0x0),
		fMcHandler(0x0),
		fStack(0),
		fIndexLeading(-1),
		fMinPtLeading(5.0),
		fMaxPtLeading(5.0),
		fDphiNS(0x0),
		fDphiAS(0x0),
		fDphiTS(0x0),
		fMultTS(0x0),
		fptL(0x0),
		fListOfObjects(0)
{

	for(Int_t i_est=0; i_est<5; ++i_est){
		for(Int_t i_pid=0; i_pid<12; ++i_pid){

			fHistPtVsNchNS[i_pid][i_est] = 0;
			fHistPtVsNchAS[i_pid][i_est] = 0;
			fHistPtVsNchTS[i_pid][i_est] = 0;


		}
	}
	for(Int_t i_est=0; i_est<5; ++i_est){// loop over mult estimators

		fMult[i_est] = 0;

	}
	for(Int_t i_est=0; i_est<5; ++i_est){// loop over mult estimators

		fMult2[i_est] = 0;

	}




	// Default constructor (should not be used)
}

//______________________________________________________________________________

AliAnalysisTaskGenMcKno::AliAnalysisTaskGenMcKno(const char *name):
	AliAnalysisTaskSE(name),
	fMcEvent(0x0),
	fMcHandler(0x0),
	fStack(0),
	fIndexLeading(-1),
	fMinPtLeading(5.0),
	fMaxPtLeading(5.0),
	fDphiNS(0x0),
	fDphiAS(0x0),
	fDphiTS(0x0),
	fMultTS(0x0),
	fptL(0x0),
	fListOfObjects(0)
{
	for(Int_t i_est=0; i_est<5; ++i_est){
		for(Int_t i_pid=0; i_pid<12; ++i_pid){

			fHistPtVsNchNS[i_pid][i_est] = 0;
			fHistPtVsNchAS[i_pid][i_est] = 0;
			fHistPtVsNchTS[i_pid][i_est] = 0;


		}
	}
	for(Int_t i_est=0; i_est<5; ++i_est){// loop over mult estimators

		fMult[i_est] = 0;

	}
	for(Int_t i_est=0; i_est<5; ++i_est){// loop over mult estimators

		fMult2[i_est] = 0;

	}



	DefineInput( 0, TChain::Class());
	DefineOutput(1, TList::Class() ); // Basic output slot 
}


//_____________________________________________________________________________

AliAnalysisTaskGenMcKno::~AliAnalysisTaskGenMcKno(){
	// Destructor
	// histograms are in the output list and deleted when the output
	// list is deleted by the TSelector dtor
	if (fListOfObjects) { delete fListOfObjects; fListOfObjects=0x0; }
}

//______________________________________________________________________________

void AliAnalysisTaskGenMcKno::UserCreateOutputObjects(){


	// ### Analysis output
	fListOfObjects = new TList();
	fListOfObjects->SetOwner(kTRUE);

	// ### Create histograms

	fDphiNS = 0;
	fDphiNS = new TH1D("hDphiNS","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(fDphiNS);

	fDphiAS = 0;
	fDphiAS = new TH1D("hDphiAS","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(fDphiAS);

	fDphiTS = 0;
	fDphiTS = new TH1D("hDphiTS","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(fDphiTS);

	fMultTS = 0;
	fMultTS = new TH1D("fMultTS","",700,-0.5,699.5);
	fListOfObjects->Add(fMultTS);

	fptL = 0;
	fptL = new TH1D("fptL","",npTBins,pTBins);
	fListOfObjects->Add(fptL);

	Double_t TSBins[nTSBins+1]={0x0};
	for(Int_t i=0;i<nTSBins;++i){
		TSBins[i]=i*1.0-0.5;
	}
	TSBins[nTSBins]=1.0*nTSBins-0.5;

	for(Int_t i_est=0; i_est<5; ++i_est){
		for(Int_t i_pid=0; i_pid<12; ++i_pid){

			fHistPtVsNchNS[i_pid][i_est]=0;
			fHistPtVsNchNS[i_pid][i_est]=new TH2D(Form("fHistPtVsNchNS%s_%s",estNames[i_est],PidNames[i_pid]), "Generated #it{p}_{T} distribution NS",nTSBins, TSBins, npTBins,pTBins);
			fListOfObjects->Add(fHistPtVsNchNS[i_pid][i_est]);

			fHistPtVsNchAS[i_pid][i_est]=0;
			fHistPtVsNchAS[i_pid][i_est]=new TH2D(Form("fHistPtVsNchAS%s_%s",estNames[i_est],PidNames[i_pid]), "Generated #it{p}_{T} distribution AS",nTSBins, TSBins, npTBins,pTBins);
			fListOfObjects->Add(fHistPtVsNchAS[i_pid][i_est]);

			fHistPtVsNchTS[i_pid][i_est]=0;
			fHistPtVsNchTS[i_pid][i_est]=new TH2D(Form("fHistPtVsNchTS%s_%s",estNames[i_est],PidNames[i_pid]), "Generated #it{p}_{T} distribution TS",nTSBins, TSBins, npTBins,pTBins);
			fListOfObjects->Add(fHistPtVsNchTS[i_pid][i_est]);


		}
	}
	for(Int_t i=0; i<5; ++i){

		fMult[i] = 0;
		fMult[i] = new TH1D(Form("fMult_%s",estNames[i]),Form("%s; Nch; #eta",estNames[i]),700,-0.5,699.5); 
		fListOfObjects->Add(fMult[i]);

	}
	for(Int_t i=0; i<5; ++i){

		fMult2[i] = 0;
		fMult2[i] = new TH1D(Form("fMult2_%s",estNames[i]),Form("%s for leading pT sel; Nch; #eta",estNames[i]),700,-0.5,699.5); 
		fListOfObjects->Add(fMult2[i]);

	}


	// ### List of outputs
	PostData(1, fListOfObjects);

}
//______________________________________________________________________________
void AliAnalysisTaskGenMcKno::Init(){
	//
	fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());

}
//______________________________________________________________________________
void AliAnalysisTaskGenMcKno::UserExec(Option_t *){

	// ### Initialize
	Init();

	// ### MC handler
	if(fMcHandler)
		fMcEvent = fMcHandler->MCEvent();
	else { if(fDebug > 1) printf("AliAnalysisTaskGenMcKno::Handler() fMcHandler = NULL\n"); return; }

	// ### MC event
	if( !fMcEvent ) { if(fDebug > 1) printf("AliAnalysisTaskGenMcKno::UserExec() fMcEvent = NULL \n"); return; }

	fStack = ((AliMCEvent*)fMcEvent)->Stack();
	if (!fStack) {
		Printf("ERROR: Could not retrieve MC stack \n");
		cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
		return;
	}

	// ### MC event selection
	Bool_t isEventMCSelected = IsMCEventSelected(fMcEvent);
	if( !isEventMCSelected ) return;


	// GENERATOR LEVEL
	vector<Int_t> mult_estimators;
	// Get the leasding particle index
	fIndexLeading = -1;
	fIndexLeading = GetIndexLeading();

	GetMultipliciy( mult_estimators );


	// mult_estimators[0]: mult in transverse side
	// mult_estimators[1]: mult in |eta|<0.8
	// mult_estimators[2]: mult in V0M
	// mult_estimators[3]: mult in V0A
	// mult_estimators[4]: mult in |eta|<1
	Bool_t fIsInel0 = kFALSE;
	if(mult_estimators[4]>0)
		fIsInel0 = kTRUE;// is INEL>0

	if(fIsInel0){
		for(Int_t i=0;i<5;++i)
			fMult[i]->Fill(mult_estimators[i]);
	}

	// RT analysis, here a cut on pTleading is applied
	TParticle* mcPartLeadingGen         = 0x0;
	if(fIndexLeading>=0){
		mcPartLeadingGen                    = (TParticle *)fMcEvent->Particle(fIndexLeading);
		if( mcPartLeadingGen->Pt()>=fMinPtLeading &&mcPartLeadingGen->Pt()<=fMaxPtLeading  ){
			fptL->Fill(mcPartLeadingGen->Pt());
			for(Int_t i=0;i<5;++i)
				fMult2[i]->Fill(mult_estimators[i]);
			MakeRTAnalysis(mult_estimators);
		}
	}

	//MemInfo_t memInfo;
	//Int_t memUsage = 0;
	//gSystem->GetMemInfo(&memInfo);
	//memUsage = memInfo.fMemUsed;
	//cout<<"mem usage="<<memUsage<<endl;
	// ### Post data for all output slots

	mult_estimators.clear();

	PostData(1, fListOfObjects);

	return;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskGenMcKno::IsMCEventSelected(TObject* obj){

	Bool_t isSelected = kTRUE;

	AliMCEvent *event = 0x0;
	event = dynamic_cast<AliMCEvent*>(obj);
	if( !event ) 
		isSelected = kFALSE;

	return isSelected;
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskGenMcKno::GetIndexLeading(){

	Double_t ptleading = 0;
	Int_t index_leading = -1;
	Bool_t isPhysPrim = kFALSE;
	Double_t qPart = 0;
	Double_t etaPart = -10;

	// ### particle loop
	for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {

		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)fMcEvent->Particle(ipart);
		if (!mcPart) continue;
		//selection of primary charged particles
		if(!(mcPart->GetPDG())) continue;
		qPart = mcPart->GetPDG()->Charge()/3.;
		if(TMath::Abs(qPart)<0.001) continue;
		isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
		if(!isPhysPrim)
			continue;

		etaPart = mcPart -> Eta();
		if(TMath::Abs(etaPart) > 0.8)continue;
		if(mcPart -> Pt()<0.15)continue;
		if(mcPart -> Pt()>ptleading){
			ptleading = mcPart->Pt();
			index_leading = ipart;
		}


	} // particle loop
	cout<<"\n";
	return index_leading;
}
//______________________________________________________________________
Int_t AliAnalysisTaskGenMcKno::GetMultipliciy(vector<Int_t> &multArray){

	// Properties leading particle
	TParticle* mcPartTmp         = 0x0;
	Double_t phiL = 0;
	// Multiplicity transverse side
	if(fIndexLeading>=0){
		mcPartTmp                    = (TParticle *)fMcEvent->Particle(fIndexLeading);
		phiL = mcPartTmp->Phi();
	}

	multArray.clear();

	Bool_t isPhysPrim = kFALSE;

	Int_t mult_TS   = 0;
	Int_t mult_Eta8   = 0;
	Int_t mult_Eta1   = 0;
	Int_t mult_VZEROM = 0;
	Int_t mult_VZEROA = 0;
	Double_t qPart = 0;
	Double_t etaPart = -10;
	Double_t phiPart = -10.0;

	// ### particle loop
	for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {

		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)fMcEvent->Particle(ipart);
		if (!mcPart) continue;
		//selection of primary charged particles
		if(!(mcPart->GetPDG())) continue;
		qPart = mcPart->GetPDG()->Charge()/3.;
		if(TMath::Abs(qPart)<0.001) continue;
		isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
		if(!isPhysPrim)
			continue;

		etaPart = mcPart -> Eta();
		if( TMath::Abs(etaPart) < 0.8 ){
			mult_Eta8++;
			// only events with well defined leading particle
			if(fIndexLeading>=0){
				phiPart = mcPart -> Phi();
				Double_t DPhi = DeltaPhi(phiL,phiPart);

				if(TMath::Abs(DPhi)<pi/3.0){
					if(ipart!=fIndexLeading)
						fDphiNS->Fill(DPhi);
				}
				// away side
				else if(TMath::Abs(DPhi-pi)<pi/3.0){
					fDphiAS->Fill(DPhi);
				}
				// transverse side
				else{
					if(mcPart -> Pt()>=0.5)
						mult_TS++;
					fDphiTS->Fill(DPhi);
				}

			}
		}
		if( (2.8 < etaPart && etaPart < 5.1) || (-3.7 < etaPart && etaPart <-1.7) ) mult_VZEROM++;
		if( 2.8 < etaPart && etaPart < 5.1 ) mult_VZEROA++;

		if( TMath::Abs(etaPart) < 1 )
			if(mcPart -> Pt()>0)
				mult_Eta1++;// for INEL>0n

	} // particle loop
	multArray.push_back(mult_TS);
	multArray.push_back(mult_Eta8);
	multArray.push_back(mult_VZEROM);
	multArray.push_back(mult_VZEROA);
	multArray.push_back(mult_Eta1);

	return 1;
}
//______________________________________________________________________________
void AliAnalysisTaskGenMcKno::MakeRTAnalysis(vector<Int_t> &mult){

	// Properties leading particle
	TParticle* mcPartTmp         = 0x0;
	mcPartTmp                    = (TParticle *)fMcEvent->Particle(fIndexLeading);

	Double_t phiL = mcPartTmp->Phi();

	// Get multiplicity in transverse side
	Int_t pidCodeMC = 0;
	Double_t ipt = 0.;
	Double_t phiPart = -10.0;
	Bool_t isPhysPrim = kFALSE;
	Double_t qPart = 0;
	Int_t pPDG = -10;

	fMultTS->Fill(mult[0]);


	// ### particle loop
	for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {

		if(ipart == fIndexLeading)continue;

		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)fMcEvent->Particle(ipart);
		if (!mcPart) continue;

		if(!mcPart->GetPDG())continue;
		isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
		qPart = mcPart->GetPDG()->Charge()/3.;
		// only primary charged particles
		pPDG = TMath::Abs(mcPart->GetPdgCode());
		pidCodeMC = GetPidCode(pPDG);
		ipt = mcPart->Pt();
		if(ipt<0.15)continue;
		phiPart = mcPart -> Phi();
		Double_t DPhi = DeltaPhi(phiL,phiPart);
		if(TMath::Abs(mcPart->Eta()) >= 0.8)continue;

		for(Int_t j=0; j<5; ++j){// mult estimators
			for(Int_t i=0; i<11; ++i)// loop over particle species
			{
				if( pidCodeMC == i )
				{
					if( iSprimary[i] == kTRUE && isPhysPrim == kFALSE ) 
						continue;


					if(TMath::Abs(DPhi)<pi/3.0){// near side
						fHistPtVsNchNS[i][j]->Fill(1.0*mult[j],ipt);
					}       
					else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
						fHistPtVsNchAS[i][j]->Fill(1.0*mult[j],ipt);
					}       
					else{// transverse side
						fHistPtVsNchTS[i][j]->Fill(1.0*mult[j],ipt);
					}       


				}
			}
			if(isPhysPrim){// primary
				if(TMath::Abs(qPart)>0.001){// charged
					if(TMath::Abs(DPhi)<pi/3.0){// near side
						fHistPtVsNchNS[11][j]->Fill(1.0*mult[j],ipt);
					}       
					else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
						fHistPtVsNchAS[11][j]->Fill(1.0*mult[j],ipt);
					}       
					else{// transverse side
						fHistPtVsNchTS[11][j]->Fill(1.0*mult[j],ipt);
					}
				}
			}
		}

	} // particle loop

}

//______________________________________________________________________________

void AliAnalysisTaskGenMcKno::Terminate(Option_t*){

	fListOfObjects = dynamic_cast<TList*> (GetOutputData(1));
	if (!fListOfObjects) { Printf("ERROR: Output list not available"); return; }

	return;
}

//_____________________________________________________________________________

Int_t AliAnalysisTaskGenMcKno::GetPidCode(Int_t pdgCode) const  {

	Int_t pidCode = 999;

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
		case 310:
			pidCode = 3; // K0s
			break;
		case 3122:
			pidCode = 4; // Lambda
			break;
		case 3312:
			pidCode = 5; // Xi-
			break;
		case 3334:
			pidCode = 6; // Omega-
			break;
		case 333:
			pidCode = 7; // phi(1020)
			break;
		case 313:
			pidCode = 8; // K*(892)0
			break;
		case 323:
			pidCode = 9; // K*(892) +-
			break;
		case 3212:
			pidCode = 10; // Sigma 0
			break;    
		default:
			break;
	};

	return pidCode;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskGenMcKno::DeltaPhi(Double_t phia, Double_t phib,
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
