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
 *  AliAnalysisTaskGenUeSpherocity.cxx
 *  
 *
 *  Author: Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)
 *  Original analysistask provided by Gyula BENCEDI  <Gyula.Bencedi@cern.ch>, WIGNER RCP
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
#include "AliAnalysisTaskGenUeSpherocity.h"

//_____ STL includes
#include <iostream>
using namespace std;

const Char_t * estimators[3]={"Mid05","Mid08","V0M"};

ClassImp( AliAnalysisTaskGenUeSpherocity )

	//_____________________________________________________________________________

	AliAnalysisTaskGenUeSpherocity::AliAnalysisTaskGenUeSpherocity():
		AliAnalysisTaskSE(),
		fMcEvent(0x0),
		fMcHandler(0x0),
		fStack(0),
		fSizeStep(0.1),
		fNrec(-1),
		fY(0.5),
		fHistEvt(0x0),
		fHistPart(0x0),
		fListOfObjects(0)
{
	// Default constructor (should not be used)
}

//______________________________________________________________________________

AliAnalysisTaskGenUeSpherocity::AliAnalysisTaskGenUeSpherocity(const char *name):
	AliAnalysisTaskSE(name),
	fMcEvent(0x0),
	fMcHandler(0x0),
	fStack(0),
	fSizeStep(0.1),
	fNrec(-1),
	fY(0.5),
	fHistEvt(0x0),
	fHistPart(0x0),
	fListOfObjects(0)
{
	DefineInput( 0, TChain::Class());
	DefineOutput(1, TList::Class() ); // Basic output slot 
}


//_____________________________________________________________________________

AliAnalysisTaskGenUeSpherocity::~AliAnalysisTaskGenUeSpherocity(){
	// Destructor
	// histograms are in the output list and deleted when the output
	// list is deleted by the TSelector dtor

	if (fHistEvt) { delete fHistEvt; fHistEvt=0x0; }
	if (fHistPart) { delete fHistPart; fHistPart=0x0; }
	if (fListOfObjects) { delete fListOfObjects; fListOfObjects=0x0; }
}

//______________________________________________________________________________

void AliAnalysisTaskGenUeSpherocity::UserCreateOutputObjects(){

	// ### Analysis output
	fListOfObjects = new TList();
	fListOfObjects->SetOwner(kTRUE);

	TString pidNames[11] = { "Pion", "Kaon", "Proton", "K0Short", "Lambda", "Xi", "Omega", "Phi", "KStar", "KStarPM", "SigmaZero" };

	// ### Create histograms
	fHistEvt = new TH1I("fHistEvt","fHistEvt",2,0,2) ;
	fHistEvt->GetXaxis()->SetBinLabel(1,"All events");
	fHistEvt->GetXaxis()->SetBinLabel(2,"All particles");
	fHistEvt->Sumw2();
	fListOfObjects->Add(fHistEvt);

	InitHisto<TH1F>("fHistMultPrimary","Multiplicity Primary", 100, 0., 2000., "N_{prim.}", "Entries");
	InitHisto<TH1F>("fHistEta","Eta Distr.", 200, -1., 1., "#eta", "N_{part}");
	InitHisto<TH1F>("fHistY", "Y Distr.", 200, -1., 1., "#it{y}", "N_{part}");

	for(Int_t i=0; i<11; i++)
		InitHisto<TH1F>(Form("fHistPt_%s",pidNames[i].Data()), "Generated #it{p}_{T} distribution",2000,0.,20., "#it{p}_{T} (GeV/#it{c})", "Entries");

	for(Int_t i=0; i<3; i++){
		InitHisto<TH2D>(Form("fMult_%s",estimators[i]),Form("Selection bias %s",estimators[i]),300,-0.5,299.5,100,-5,5,"Nch","#eta");
		InitHisto<TH2D>(Form("fSoVsNch_%s",estimators[i]),Form("So vs Nch %s",estimators[i]),300,-0.5,299.5,200,0.0,1.0,"Nch","Spherocity");
		InitHisto<TH2D>(Form("fSoWeighedVsNch_%s",estimators[i]),Form("So vs Nch %s",estimators[i]),300,-0.5,299.5,200,0.0,1.0,"Nch","Spherocity (pT weighted)");

		InitHisto<TH1F>(Form("fNch_%s",estimators[i]),Form("Nch %s",estimators[i]),300,-0.5,299.5,"Nch","counts");
		InitHisto<TH1F>(Form("fNchSoSel_%s",estimators[i]),Form("Nch %s (after So cuts)",estimators[i]),300,-0.5,299.5,"Nch","counts");
	}
	// ### List of outputs
	PostData(1, fListOfObjects);

}

//______________________________________________________________________________

inline void AliAnalysisTaskGenUeSpherocity::FillHisto(const char* objkey, Double_t x)
{
	TH1* hTmp = static_cast<TH1*>(fListOfObjects->FindObject(objkey));
	if(!hTmp){
		AliError(Form("Cannot find histogram: %s",objkey)) ;
		return;
	}
	hTmp->Fill(x);
}

inline void AliAnalysisTaskGenUeSpherocity::FillHisto(const char* objkey, Double_t x, Double_t y)
{
	TH2* hTmp = static_cast<TH2*>(fListOfObjects->FindObject(objkey));
	if(!hTmp){
		AliError(Form("Cannot find histogram: %s",objkey)) ;
		return;
	}
	hTmp->Fill(x,y);
}

//______________________________________________________________________________

template <class T> T* AliAnalysisTaskGenUeSpherocity::InitHisto(const char* hname, const char* htitle, Int_t nxbins, Double_t xmin, Double_t xmax, const char* xtitle, const char* ytitle)
{
	T* hTmp = new T(hname, htitle, nxbins, xmin, xmax);
	hTmp->GetXaxis()->SetTitle(xtitle);
	hTmp->GetYaxis()->SetTitle(ytitle);
	hTmp->SetMarkerStyle(kFullCircle);
	hTmp->Sumw2();
	fListOfObjects->Add(hTmp);

	return hTmp;
}
template <class T> T* AliAnalysisTaskGenUeSpherocity::InitHisto(const char* hname, const char* htitle, Int_t nxbins, Double_t xmin, Double_t xmax, Int_t nybins, Double_t ymin, Double_t ymax, const char* xtitle, const char* ytitle)
{
	T* hTmp = new T(hname, htitle, nxbins, xmin, xmax, nybins, ymin, ymax);
	hTmp->GetXaxis()->SetTitle(xtitle);
	hTmp->GetYaxis()->SetTitle(ytitle);
	hTmp->Sumw2();
	fListOfObjects->Add(hTmp);

	return hTmp;
}

//______________________________________________________________________________

void AliAnalysisTaskGenUeSpherocity::Init(){
	//
	fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
}

//______________________________________________________________________________

void AliAnalysisTaskGenUeSpherocity::UserExec(Option_t *){

	// ### Initialize
	Init();

	// ### MC handler
	if(fMcHandler)
		fMcEvent = fMcHandler->MCEvent();
	else { if(fDebug > 1) printf("AliAnalysisTaskGenUeSpherocity::Handler() fMcHandler = NULL\n"); return; }

	// ### MC event
	if( !fMcEvent ) { if(fDebug > 1) printf("AliAnalysisTaskGenUeSpherocity::UserExec() fMcEvent = NULL \n"); return; }

	fStack = ((AliMCEvent*)fMcEvent)->Stack();
	if (!fStack) {
		Printf("ERROR: Could not retrieve MC stack \n");
		cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
		return;
	}

	// ### MC event selection
	Bool_t isEventMCSelected = IsMCEventSelected(fMcEvent);
	if( !isEventMCSelected ) return;

	vector<Int_t> mult_estimators;
	vector<Float_t> pt_so;
	vector<Float_t> eta_so;
	vector<Float_t> phi_so;

	fNrec = -1;
	fNrec = GetMultipliciy( mult_estimators, pt_so, eta_so, phi_so, fMcEvent );

	Float_t spherocity_ptWeighted=-0.5;
	Float_t spherocity=-0.5;

	if(fNrec>2){
		spherocity_ptWeighted = GetSpherocity(pt_so, eta_so, phi_so, kTRUE);
		spherocity = GetSpherocity(pt_so, eta_so, phi_so, kFALSE);
	}
	// ### Event and particle selection
	EventSel( fMcEvent );
	for(Int_t i=0;i<3;++i)
		FillHisto(Form("fNch_%s",estimators[i]),1.0*mult_estimators[i]);

	if(spherocity>=0&&spherocity<=1){

		for(Int_t i=0;i<3;++i){
			FillHisto(Form("fNchSoSel_%s",estimators[i]),1.0*mult_estimators[i]);
			FillHisto(Form("fSoVsNch_%s",estimators[i]),1.0*mult_estimators[i],spherocity);
		}

	}
	if(spherocity_ptWeighted>=0&&spherocity_ptWeighted<=1){

		for(Int_t i=0;i<3;++i){
			FillHisto(Form("fSoWeighedVsNch_%s",estimators[i]),1.0*mult_estimators[i],spherocity_ptWeighted);
		}
		ParticleSel( mult_estimators, fMcEvent );
	}



	// ### Post data for all output slots
	PostData(1, fListOfObjects);

	return;
}

//______________________________________________________________________________

Bool_t AliAnalysisTaskGenUeSpherocity::IsMCEventSelected(TObject* obj){

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

void AliAnalysisTaskGenUeSpherocity::EventSel(TObject* obj){

	if( !obj ) return;
	AliMCEvent *event = dynamic_cast<AliMCEvent*>(obj);
	if( !event ) return;

	Int_t multPrimary = fStack->GetNprimary();
	FillHisto("fHistMultPrimary",multPrimary);

}
//_____________________________________________________________________________

Int_t AliAnalysisTaskGenUeSpherocity::GetMultipliciy(vector<Int_t> &multArray, vector<Float_t> &ptArray,  vector<Float_t> &etaArray, vector<Float_t> &phiArray, TObject* obj){

	Int_t mult_so=0;

	multArray.clear();
	ptArray.clear();
	etaArray.clear();
	phiArray.clear();


	if ( !obj ) return -1;

	AliMCEvent *event = dynamic_cast<AliMCEvent*>(obj);
	if ( !event ) return -1;

	Bool_t isPhysPrim = kFALSE;

	Int_t mult_Eta5   = 0;
	Int_t mult_Eta8   = 0;
	Int_t mult_VZEROM = 0;

	// ### particle loop
	for (Int_t ipart = 0; ipart < event->GetNumberOfTracks(); ipart++) {

		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)event->Particle(ipart);
		if (!mcPart) continue;
		//selection of primary charged particles
		if(!mcPart->GetPDG()) continue;
		Double_t qPart = mcPart->GetPDG()->Charge()/3.;
		if(TMath::Abs(qPart)<0.001) continue;
		isPhysPrim = event->IsPhysicalPrimary(ipart);
		if(!isPhysPrim)
			continue;

		Double_t etaPart = mcPart -> Eta();
		if( TMath::Abs(etaPart) < 0.5 ) mult_Eta5++;
		if( (2.8 < etaPart && etaPart < 5.1) || (-3.7 < etaPart && etaPart <-1.7) ) mult_VZEROM++;
		if( TMath::Abs(etaPart) < 0.8 ){ 
			mult_Eta8++;

			if(mcPart -> Pt()>0.15){

				ptArray.push_back(mcPart->Pt());
				etaArray.push_back(mcPart->Eta());
				phiArray.push_back(mcPart->Phi());
				mult_so++;
			}

		}

	} // particle loop

	multArray.push_back(mult_Eta5);
	multArray.push_back(mult_Eta8);
	multArray.push_back(mult_VZEROM);

	return mult_so;
}

//______________________________________________________________________________

void AliAnalysisTaskGenUeSpherocity::ParticleSel(const vector<Int_t> &mult, TObject* obj){



	TString pidNames[11] = { "Pion", "Kaon", "Proton", "K0Short", "Lambda", "Xi", "Omega", "Phi", "KStar", "KStarPM", "SigmaZero" };

	if ( !obj ) return;

	AliMCEvent *event = dynamic_cast<AliMCEvent*>(obj);
	if ( !event ) return;

	Bool_t isPrimary[11] = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE };

	Int_t pidCodeMC = 0;
	Double_t ipt = 0.;

	Bool_t isPhysPrim = kFALSE;


	// ### particle loop
	for (Int_t ipart = 0; ipart < event->GetNumberOfTracks(); ipart++) {

		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)event->Particle(ipart);
		if (!mcPart) continue;

		FillHisto("fHistEvt",1.5);

		// only primary charged particles
		if( event->IsPhysicalPrimary(ipart) ){
			Double_t qPart = mcPart->GetPDG()->Charge()/3.;
			if(TMath::Abs(qPart)>0.001)
				for(Int_t i=0;i<3;i++)
					FillHisto(Form("fMult_%s",estimators[i]),1.0*mult[i],mcPart->Eta());
		}

		Int_t pPDG = TMath::Abs(mcPart->GetPdgCode());
		pidCodeMC = GetPidCode(pPDG);

		Bool_t isSelectedPart = kTRUE;
		for(Int_t i=0; i<11; i++) 
			if( pidCodeMC == i ) 
				isSelectedPart = kFALSE;
		if ( isSelectedPart ) continue;

		FillHisto("fHistEta",mcPart->Eta());


		if (!(TMath::Abs(mcPart->Energy()-mcPart->Pz())>0.)) continue;
		Double_t myY = (mcPart->Energy()+mcPart->Pz())/(mcPart->Energy()-mcPart->Pz());
		if( myY <= 0 ) continue;


		Double_t y = 0.5*TMath::Log(myY);

		ipt = mcPart->Pt();

		isPhysPrim = event->IsPhysicalPrimary(ipart);

		for(Int_t i=0; i<11; ++i)
		{
			if( pidCodeMC == i && TMath::Abs(y) < fY)
			{
				if( isPrimary[i] == kTRUE && isPhysPrim == kFALSE ) 
					continue;

				if(!i) FillHisto("fHistY",y);
				FillHisto(Form("fHistPt_%s",pidNames[i].Data()),ipt);


			}
		}

	} // particle loop

}

Float_t AliAnalysisTaskGenUeSpherocity::GetSpherocity( const vector<Float_t> &pt, const vector<Float_t> &eta, const vector<Float_t> &phi, const Bool_t isPtWeighted ){


	Float_t spherocity = -10.0;
	Float_t pFull = 0;
	Float_t Spherocity = 2;

	//computing total pt
	Float_t sumapt = 0;
	if(isPtWeighted)
		for(Int_t i1 = 0; i1 < fNrec; ++i1){
			sumapt += pt[i1];
		}
	else
		sumapt = 1.0*fNrec;
	//Getting thrust
	for(Int_t i = 0; i < 360/(fSizeStep); ++i){
		Float_t numerador = 0;
		Float_t phiparam  = 0;
		Float_t nx = 0;
		Float_t ny = 0;
		phiparam=( (TMath::Pi()) * i * fSizeStep ) / 180; // parametrization of the angle
		nx = TMath::Cos(phiparam);            // x component of an unitary vector n
		ny = TMath::Sin(phiparam);            // y component of an unitary vector n
		for(Int_t i1 = 0; i1 < fNrec; ++i1){

			Float_t pxA = 0;
			Float_t pyA = 0;
			if(isPtWeighted){
				pxA = pt[i1] * TMath::Cos( phi[i1] );
				pyA = pt[i1] * TMath::Sin( phi[i1] );
			}
			else{
				pxA = 1.0 * TMath::Cos( phi[i1] );
				pyA = 1.0 * TMath::Sin( phi[i1] );
			}



			numerador += TMath::Abs( ny * pxA - nx * pyA );//product between p  proyection in XY plane and the unitary vector
		}
		pFull=TMath::Power( (numerador / sumapt),2 );
		if(pFull < Spherocity)//maximization of pFull
		{
			Spherocity = pFull;
		}
	}

	spherocity=((Spherocity)*TMath::Pi()*TMath::Pi())/4.0;


	return spherocity;

}


//______________________________________________________________________________

void AliAnalysisTaskGenUeSpherocity::Terminate(Option_t*){

	fListOfObjects = dynamic_cast<TList*> (GetOutputData(1));
	if (!fListOfObjects) { Printf("ERROR: Output list not available"); return; }

	return;
}

//_____________________________________________________________________________

Int_t AliAnalysisTaskGenUeSpherocity::GetPidCode(Int_t pdgCode) const  {

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


