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
 *  AliAnalysisTaskGenRT.cxx
 *  
 *  Minimal analysis task meant to be used for on-the-fly simulations in LEGO trains (for generating light falvor particle pt spectra)
 *
 *  Gyula BENCEDI  <Gyula.Bencedi@cern.ch>, WIGNER RCP
 * 
 */

//_____ ROOT headers
#include <TList.h>
#include <TTree.h>
#include <TMath.h>
#include <TFile.h>
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
#include "AliAnalysisTaskGenRT.h"

//_____ STL includes
#include <iostream>
using namespace std;

const Char_t *estNames[6] = {"TS","Eta08","V0M","V0A","V0C","Eta10"};
const Char_t *pidNames[11] = { "Pion", "Kaon", "Proton", "K0Short", "Lambda", "Xi", "Omega", "Phi", "KStar", "KStarPM", "SigmaZero" };
Bool_t isPrimary[11] = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE };

ClassImp( AliAnalysisTaskGenRT )

	//_____________________________________________________________________________

	AliAnalysisTaskGenRT::AliAnalysisTaskGenRT():
		AliAnalysisTaskSE(),
		fMcEvent(0x0),
		fMcHandler(0x0),
		fStack(0),
		AllowXi(kFALSE),
		fY(0.5),
		fIndexLeadingGen(-1),
		fMinPtLeading(5.0),
		fMaxPtLeading(40.0),
		fHistEvt(0x0),
		fHistPart(0x0),
		hPtLeadingGen(0x0),
		hMultTSvsPtLeading(0x0),
		fDphiNS(0x0),
		fDphiAS(0x0),
		fDphiTS(0x0),
		fMultTS(0x0),
		fListOfObjects(0)
{

	for(Int_t i_pid=0; i_pid<11; ++i_pid){

		fHistPtVsNchNS[i_pid]=0;
		fHistPtVsNchAS[i_pid]=0;
		fHistPtVsNchTS[i_pid]=0;

	}

	for(Int_t i_est=0; i_est<6; ++i_est){// loop over mult estimators

		fMult[i_est] = 0;
		fMult2[i_est] = 0;

	}

	// Default constructor (should not be used)
}

//______________________________________________________________________________

AliAnalysisTaskGenRT::AliAnalysisTaskGenRT(const char *name):
	AliAnalysisTaskSE(name),
	fMcEvent(0x0),
	fMcHandler(0x0),
	fStack(0),
	AllowXi(kFALSE),
	fY(0.5),
	fIndexLeadingGen(-1),
	fMinPtLeading(5.0),
	fMaxPtLeading(40.0),
	fHistEvt(0x0),
	fHistPart(0x0),
	hPtLeadingGen(0x0),
	hMultTSvsPtLeading(0x0),
	fDphiNS(0x0),
	fDphiAS(0x0),
	fDphiTS(0x0),
	fMultTS(0x0),
	fListOfObjects(0)
{

	for(Int_t i_pid=0; i_pid<11; ++i_pid){

		fHistPtVsNchNS[i_pid]=0;
		fHistPtVsNchAS[i_pid]=0;
		fHistPtVsNchTS[i_pid]=0;

	}

	for(Int_t i_est=0; i_est<6; ++i_est){// loop over mult estimators

		fMult[i_est] = 0;
		fMult2[i_est] = 0;

	}

	DefineInput( 0, TChain::Class());
	DefineOutput(1, TList::Class() ); // Basic output slot 
}


//_____________________________________________________________________________

AliAnalysisTaskGenRT::~AliAnalysisTaskGenRT(){
	// Destructor
	// histograms are in the output list and deleted when the output
	// list is deleted by the TSelector dtor

	if (fHistEvt) { delete fHistEvt; fHistEvt=0x0; }
	if (fHistPart) { delete fHistPart; fHistPart=0x0; }
	if (fListOfObjects) { delete fListOfObjects; fListOfObjects=0x0; }
}

//______________________________________________________________________________

void AliAnalysisTaskGenRT::UserCreateOutputObjects(){

	// ### Analysis output
	fListOfObjects = new TList();
	fListOfObjects->SetOwner(kTRUE);

	////	TString pidNames[11] = { "Pion", "Kaon", "Proton", "K0Short", "Lambda", "Xi", "Omega", "Phi", "KStar", "KStarPM", "SigmaZero" };

	// ### Create histograms
	fHistEvt = new TH1I("fHistEvt","fHistEvt",2,0,2) ;
	fHistEvt->GetXaxis()->SetBinLabel(1,"All events");
	fHistEvt->GetXaxis()->SetBinLabel(2,"All particles");
	fHistEvt->Sumw2();
	fListOfObjects->Add(fHistEvt);

	InitHisto<TH1F>("fHistEta","Eta Distr.", 200, -1., 1., "#eta", "N_{part}");
	InitHisto<TH1F>("fHistY", "Y Distr.", 200, -1., 1., "#it{y}", "N_{part}");

	const Int_t nTSBins = 100;
	Double_t TSBins[nTSBins+1]={0x0};
	for(Int_t i=0;i<nTSBins;++i){
		TSBins[i]=i*1.0-0.5;
	}
	TSBins[nTSBins]=99.5;

	const Int_t nPtBins = 45;
	Double_t PtBins[nPtBins+1] = {
		0.0, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
		0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
		0.90, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
		1.80, 1.90, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
		3.60, 3.80, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10.0};

	const Int_t nPtBinsLeading = 24;
	Double_t PtBinsLeading[nPtBinsLeading+1] = {
		0.0 , 0.15 , 1.0 , 1.5  , 2.0  , 2.5  , 3.0  , 3.5  , 4.0  , 4.5  , 5.0  , 6.0, 
		7.0 , 8.0 , 9.0 , 10.0 , 12.0 , 14.0 , 16.0 , 18.0 , 20.0 , 25.0 , 30.0 , 40.0 , 100.0 };

	hPtLeadingGen = 0;
	hPtLeadingGen = new TH1F("hPtLeadingGen","; #it{p}_{T} (GeV/#it{c}) ; Entries", nPtBinsLeading, PtBinsLeading);
	fListOfObjects->Add(hPtLeadingGen);

	hMultTSvsPtLeading = 0;
	hMultTSvsPtLeading = new TH2F("hMultTSvsPtLeading","; NT ; #it{p}_{T} (GeV/#it{c})", nTSBins, TSBins, nPtBinsLeading, PtBinsLeading);
	fListOfObjects->Add(hMultTSvsPtLeading);

	fDphiNS = 0;
	fDphiNS = new TH1F("hDphiNS","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(fDphiNS);

	fDphiAS = 0;
	fDphiAS = new TH1F("hDphiAS","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(fDphiAS);

	fDphiTS = 0;
	fDphiTS = new TH1F("hDphiTS","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(fDphiTS);

	fMultTS = 0;
	fMultTS = new TH1F("fMultTS","; NT ; Entries", nTSBins, TSBins);
	fListOfObjects->Add(fMultTS);

	for(Int_t i=0; i<6; ++i){

		fMult[i] = 0;
		fMult[i] = new TH1F(Form("fMult_%s",estNames[i]),Form("%s; Nch; #eta",estNames[i]),3000,-0.5,2999.5); 
		fListOfObjects->Add(fMult[i]);

		fMult2[i] = 0;
		fMult2[i] = new TH1F(Form("fMult2_%s",estNames[i]),Form("%s; Nch; #eta",estNames[i]),3000,-0.5,2999.5); 
		fListOfObjects->Add(fMult2[i]);

	}

	for(Int_t i_pid=0; i_pid<11; ++i_pid){

		fHistPtVsNchNS[i_pid]=0;
		fHistPtVsNchNS[i_pid]=new TH2F(Form("fHistPtVsNchNS_%s",pidNames[i_pid]), "Generated #it{p}_{T} distribution NS; NT ; #it{p}_{T} (GeV/#it{c})",nTSBins, TSBins, nPtBins,PtBins);
		fListOfObjects->Add(fHistPtVsNchNS[i_pid]);

		fHistPtVsNchAS[i_pid]=0;
		fHistPtVsNchAS[i_pid]=new TH2F(Form("fHistPtVsNchAS_%s",pidNames[i_pid]), "Generated #it{p}_{T} distribution AS; NT ; #it{p}_{T} (GeV/#it{c})",nTSBins, TSBins, nPtBins,PtBins);
		fListOfObjects->Add(fHistPtVsNchAS[i_pid]);

		fHistPtVsNchTS[i_pid]=0;
		fHistPtVsNchTS[i_pid]=new TH2F(Form("fHistPtVsNchTS_%s",pidNames[i_pid]), "Generated #it{p}_{T} distribution TS; NT ; #it{p}_{T} (GeV/#it{c})",nTSBins, TSBins, nPtBins,PtBins);
		fListOfObjects->Add(fHistPtVsNchTS[i_pid]);

	}

	// ### List of outputs
	PostData(1, fListOfObjects);

}

//______________________________________________________________________________

inline void AliAnalysisTaskGenRT::FillHisto(const char* objkey, Double_t x)
{
	TH1* hTmp = static_cast<TH1*>(fListOfObjects->FindObject(objkey));
	if(!hTmp){
		AliError(Form("Cannot find histogram: %s",objkey)) ;
		return;
	}
	hTmp->Fill(x);
}

//______________________________________________________________________________

template <class T> T* AliAnalysisTaskGenRT::InitHisto(const char* hname, const char* htitle, Int_t nxbins, Double_t xmin, Double_t xmax, const char* xtitle, const char* ytitle)
{
	T* hTmp = new T(hname, htitle, nxbins, xmin, xmax);
	hTmp->GetXaxis()->SetTitle(xtitle);
	hTmp->GetYaxis()->SetTitle(ytitle);
	fListOfObjects->Add(hTmp);

	return hTmp;
}

//______________________________________________________________________________

void AliAnalysisTaskGenRT::Init(){
	//
	fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
}

//______________________________________________________________________________

void AliAnalysisTaskGenRT::UserExec(Option_t *){

	// ### Initialize
	Init();

	// ### MC handler
	if(fMcHandler)
		fMcEvent = fMcHandler->MCEvent();
	else { if(fDebug > 1) printf("AliAnalysisTaskGenRT::Handler() fMcHandler = NULL\n"); return; }

	// ### MC event
	if( !fMcEvent ) { if(fDebug > 1) printf("AliAnalysisTaskGenRT::UserExec() fMcEvent = NULL \n"); return; }

	fStack = ((AliMCEvent*)fMcEvent)->Stack();
	if (!fStack) {
		Printf("ERROR: Could not retrieve MC stack \n");
		cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
		return;
	}

	// ### MC event selection
	Bool_t isEventMCSelected = IsMCEventSelected(fMcEvent);
	if( !isEventMCSelected ) return;

	//! Get index of the Leading particle
	fIndexLeadingGen = -1;
	fIndexLeadingGen = GetIndexLeading();

	//! Multiplicity estimators
	vector<Int_t> mult_estimators;
	GetMultipliciy( mult_estimators );

	// mult_estimators[0]: mult in transverse side
	// mult_estimators[1]: mult in |eta|<0.8
	// mult_estimators[2]: mult in V0M
	// mult_estimators[3]: mult in V0A
	// mult_estimators[4]: mult in V0C
	// mult_estimators[5]: mult in |eta|<1
	Bool_t fIsInel0 = kFALSE;
	if(mult_estimators[5]>0)
		fIsInel0 = kTRUE;// is INEL>0

	if(fIsInel0){
		for(Int_t i=0;i<6;++i)
			fMult[i]->Fill(mult_estimators[i]);
	}

	//! Apply cut on the leading particle
	TParticle* mcPartLeadingGen = 0x0;
	if( fIndexLeadingGen >= 0 ){
		mcPartLeadingGen  = (TParticle *)fMcEvent->Particle( fIndexLeadingGen );
		hPtLeadingGen->Fill( mcPartLeadingGen->Pt() );
		if( ( mcPartLeadingGen->Pt() >= fMinPtLeading ) && ( mcPartLeadingGen->Pt() < fMaxPtLeading ) ){

			//! Fill Multiplicity with leading particle >= 5 GeV/c
			for(Int_t i=0;i<6;++i) { 
				fMult2[i]->Fill(mult_estimators[i]);
			}
			//! PID as a function of RT
			MakeRTAnalysis();
		}
	}

	// ### Post data for all output slots
	PostData(1, fListOfObjects);

	return;
}

//______________________________________________________________________________

Bool_t AliAnalysisTaskGenRT::IsMCEventSelected(TObject* obj){

	Bool_t isSelected = kTRUE;

	AliMCEvent *event = 0x0;
	event = dynamic_cast<AliMCEvent*>(obj);
	if( !event ) 
		isSelected = kFALSE;

	if( isSelected ) 
		FillHisto("fHistEvt",0.5);

	return isSelected;
}

//______________________________________________________________________________

Int_t AliAnalysisTaskGenRT::GetIndexLeading(){

	Double_t ptleading = 0.0;
	Int_t index_leading = -1;
	Bool_t isPhysPrim = kFALSE;
	Double_t qPart = 0.0;
	Double_t etaPart = -10.0;
	Double_t yPart = -10.0;
	Int_t pidCodeMC = 0;
	Int_t pPDG = -10;

	// ### particle loop
	for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {

		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)fMcEvent->Particle(ipart);
		if (!mcPart) continue;
		//! selection of primary charged particles
		if(!(mcPart->GetPDG())) continue;
		qPart = mcPart->GetPDG()->Charge()/3.;
		if(TMath::Abs(qPart)<0.001) continue;
		isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
		if(!isPhysPrim) continue;

		etaPart = mcPart->Eta();
		yPart = mcPart->Y();

		FillHisto("fHistEvt",1.5);
		FillHisto("fHistEta",etaPart);
		FillHisto("fHistY",yPart);

		pPDG = TMath::Abs(mcPart->GetPdgCode());
		pidCodeMC = GetPidCode(pPDG);

		//! Allow the leading to be a Xi
		if(!AllowXi) {
			if( pidCodeMC == 5 ) { continue; }
		}

		if(TMath::Abs(etaPart) > 0.8) { continue; }
		if(mcPart->Pt() < 0.15) { continue; }
		if(mcPart->Pt()>ptleading){
			ptleading = mcPart->Pt();
			index_leading = ipart;
		}

	} // particle loop

	return index_leading;

}

//______________________________________________________________________________

void AliAnalysisTaskGenRT::MakeRTAnalysis(){

	// Properties leading particle
	TParticle* mcPartTmp         = 0x0;
	mcPartTmp                    = (TParticle *)fMcEvent->Particle(fIndexLeadingGen);

	Double_t phiL = mcPartTmp->Phi();
	// Multiplicity transverse side
	Int_t multTS = 0;

	// Get multiplicity in transverse side
	Int_t pidCodeMC = 0;
	Double_t ipt = 0.;
	Double_t etaPart = -10.0;
	Double_t phiPart = -10.0;
	Bool_t isPhysPrim = kFALSE;
	Double_t qPart = 0;
	Int_t pPDG = -10;
	Double_t pi = TMath::Pi();

	// ### particle loop
	for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {

		if(ipart == fIndexLeadingGen) { continue; }

		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)fMcEvent->Particle(ipart);
		if (!mcPart) continue;
		//selection of primary charged particles
		if(!(mcPart->GetPDG())) continue;
		qPart = mcPart->GetPDG()->Charge()/3.;
		if(TMath::Abs(qPart)<0.001) continue;
		isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
		if(!isPhysPrim) continue;
		etaPart = mcPart->Eta();
		if(TMath::Abs(etaPart) > 0.8) continue;
		ipt = mcPart->Pt();
		if(ipt < 0.15) continue;
		phiPart = mcPart -> Phi();
		Double_t DPhi = DeltaPhi(phiL,phiPart);

		if(TMath::Abs(DPhi)<pi/3.0){
			fDphiNS->Fill(DPhi);
		}
		// away side
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			fDphiAS->Fill(DPhi);
		}
		// transverse side
		else{
			multTS++;
			fDphiTS->Fill(DPhi);
		}
	}

	fMultTS->Fill(multTS);
	hMultTSvsPtLeading->Fill(multTS,mcPartTmp->Pt());

	// selecting topological regions
	pidCodeMC = 0;
	ipt = 0.;
	etaPart = -10.0;
	phiPart = -10.0;
	isPhysPrim = kFALSE;
	qPart = 0;
	pPDG = -10;

	// ### particle loop
	for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {

		if(ipart == fIndexLeadingGen)continue;

		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)fMcEvent->Particle(ipart);
		if (!mcPart) continue;

		if(!mcPart->GetPDG())continue;
		isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
		qPart = mcPart->GetPDG()->Charge()/3.;
		pPDG = TMath::Abs(mcPart->GetPdgCode());
		pidCodeMC = GetPidCode(pPDG);
		Bool_t isSelectedPart = kTRUE;
		for(Int_t i=0; i<11; ++i){ 
			if( pidCodeMC == i ){
				isSelectedPart = kFALSE;
			}
		}

		if(isSelectedPart) continue;
		ipt = mcPart->Pt();
		if(ipt < 0.15) continue;
		phiPart = mcPart -> Phi();
		Double_t DPhi = DeltaPhi(phiL,phiPart);
		if(TMath::Abs(mcPart->Eta()) > 0.8) continue;

		for(Int_t i=0; i<11; ++i)
		{
			if( pidCodeMC == i )
			{
				if( isPrimary[i] == kTRUE && isPhysPrim == kFALSE ) 
					continue;

				if(TMath::Abs(DPhi)<pi/3.0){// near side
					fHistPtVsNchNS[i]->Fill(1.0*multTS,ipt);
				}       
				else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
					fHistPtVsNchAS[i]->Fill(1.0*multTS,ipt);
				}       
				else{// transverse side
					fHistPtVsNchTS[i]->Fill(1.0*multTS,ipt);
				}       
			}
		}

	} // particle loop


}

//______________________________________________________________________________

Int_t AliAnalysisTaskGenRT::GetMultipliciy(vector<Int_t> &multArray){

	// Properties leading particle
	TParticle* mcPartTmp         = 0x0;
	Double_t phiL = 0;
	// Multiplicity transverse side
	if(fIndexLeadingGen>=0){
		mcPartTmp                    = (TParticle *)fMcEvent->Particle(fIndexLeadingGen);
		phiL = mcPartTmp->Phi();
	}

	multArray.clear();

	Bool_t isPhysPrim = kFALSE;

	Int_t mult_TS   = 0;
	Int_t mult_Eta8   = 0;
	Int_t mult_Eta1   = 0;
	Int_t mult_VZEROM = 0;
	Int_t mult_VZEROA = 0;
	Int_t mult_VZEROC = 0;
	Double_t qPart = 0;
	Double_t etaPart = -10;
	Double_t phiPart = -10.0;
	Double_t pi = TMath::Pi();

	// ### particle loop
	for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {

		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)fMcEvent->Particle(ipart);
		if (!mcPart) continue;
		//! selection of primary charged particles
		if(!(mcPart->GetPDG())) continue;
		qPart = mcPart->GetPDG()->Charge()/3.;
		if(TMath::Abs(qPart)<0.001) continue;
		isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
		if(!isPhysPrim) continue;

		etaPart = mcPart->Eta();
		if( TMath::Abs(etaPart) < 0.8 ){
			mult_Eta8++;
			// only events with well defined leading particle
			if(fIndexLeadingGen>=0){
				phiPart = mcPart -> Phi();
				Double_t DPhi = DeltaPhi(phiL,phiPart);

				//! Near side
				if(TMath::Abs(DPhi)<pi/3.0){
					continue;
				}
				//! Away side
				else if(TMath::Abs(DPhi-pi)<pi/3.0){
					continue;
				}
				//! Transverse side
				else{
					if(mcPart -> Pt()>=0.15){
						mult_TS++;
					}
				}

			}
		}

		if( (2.8 < etaPart && etaPart < 5.1) || (-3.7 < etaPart && etaPart <-1.7) ) mult_VZEROM++;
		if( 2.8 < etaPart && etaPart < 5.1 ) mult_VZEROA++;
		if( -3.7 < etaPart && etaPart <-1.7 ) mult_VZEROC++;

		if( TMath::Abs(etaPart) < 1.0 ){
			if(mcPart -> Pt()>0.0){
				mult_Eta1++;// for INEL>0n
			}
		}

	} // particle loop

	multArray.push_back(mult_TS);
	multArray.push_back(mult_Eta8);
	multArray.push_back(mult_VZEROM);
	multArray.push_back(mult_VZEROA);
	multArray.push_back(mult_VZEROC);
	multArray.push_back(mult_Eta1);

	return 1;
}

//______________________________________________________________________________

Double_t AliAnalysisTaskGenRT::DeltaPhi(Double_t phia, Double_t phib,
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

//______________________________________________________________________________

void AliAnalysisTaskGenRT::Terminate(Option_t*){

	fListOfObjects = dynamic_cast<TList*> (GetOutputData(1));
	if (!fListOfObjects) { Printf("ERROR: Output list not available"); return; }

	return;
}

//_____________________________________________________________________________

Int_t AliAnalysisTaskGenRT::GetPidCode(Int_t pdgCode) const  {

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
