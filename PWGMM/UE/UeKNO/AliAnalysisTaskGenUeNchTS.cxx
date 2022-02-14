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
 *  AliAnalysisTaskGenUeNchTS.cxx
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
#include "AliAnalysisTaskGenUeNchTS.h"

//_____ STL includes
#include <iostream>
using namespace std;

const Char_t * Estimators[3]={"Mid05","Mid08","V0M"};
const Int_t NchPercBin=7;
const Int_t NchBin08=122;// multiplicity |eta|<0.8
const Int_t nTSBins=100;
Double_t NchBin_gen0[NchPercBin+1]={0x0};
Double_t NchBin_gen1[NchPercBin+1]={0x0};
Double_t NchBin_gen2[NchPercBin+1]={0x0};
Double_t NchBin_rec0[NchPercBin+1]={0x0};
Double_t NchBin_rec1[NchPercBin+1]={0x0};
Double_t NchBin_rec2[NchPercBin+1]={0x0};

//TF1* ch_Eff; // efficiency for pions

const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;
Bool_t IsPrimary[11] = { kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kFALSE, kFALSE, kFALSE };
const Char_t * NameReg[3]={"NS","AS","TS"};
const Int_t NchNbins = 100;
Double_t Nchbins[NchNbins+1]={-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5,29.5,30.5,31.5,32.5,33.5,34.5,35.5,36.5,37.5,38.5,39.5,40.5,41.5,42.5,43.5,44.5,45.5,46.5,47.5,48.5,49.5,50.5,51.5,52.5,53.5,54.5,55.5,56.5,57.5,58.5,59.5,60.5,61.5,62.5,63.5,64.5,65.5,66.5,67.5,68.5,69.5,70.5,71.5,72.5,73.5,74.5,75.5,76.5,77.5,78.5,79.5,80.5,81.5,82.5,83.5,84.5,85.5,86.5,87.5,88.5,89.5,90.5,91.5,92.5,93.5,94.5,95.5,96.5,97.5,98.5,99.5};
const Int_t NPtBins = 36;
Double_t PTBins[NPtBins+1] = {
	0.0,  0.1,  0.15,  0.2,  0.25,  0.3,  0.35,  0.4,  0.45,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5, 5.0, 6.0,7.0,  8.0,  9.0,  10.0,  12.0,  14.0,  16.0,  18.0,  20.0,  25.0,  30.0,  40.0,  50.0
};

ClassImp( AliAnalysisTaskGenUeNchTS )

	//_____________________________________________________________________________

	AliAnalysisTaskGenUeNchTS::AliAnalysisTaskGenUeNchTS():
		AliAnalysisTaskSE(),
		fMcEvent(0x0),
		fMcHandler(0x0),
		fStack(0),
		fGenerator("Pythia8"),
		fIndexLeadingGen(-1),
		fIndexLeadingRec(-1),
		fMinPtLeading(5.0),
		fSizeStep(0.1),
		//fNso_gen(3),
		//fNso_rec(3),
		fbinPerc0_gen(0),
		fbinPerc1_gen(0),
		fbinPerc2_gen(0),
		fbinPerc0_rec(0),
		fbinPerc1_rec(0),
		fbinPerc2_rec(0),
		fY(0.5),
		fHistEvt(0x0),
		fHistEta(0x0),
		fHistPart(0x0),
		fDphiNS(0x0),
		fDphiAS(0x0),
		fDphiTS(0x0),
		fMultTS(0x0),
		fDphiNSRec(0x0),
		fDphiASRec(0x0),
		fDphiTSRec(0x0),
		fMultTSRec(0x0),
		fListOfObjects(0)
{

	fHistPtVsUENS=0;
	fHistPtVsUEAS=0;
	fHistPtVsUETS=0;

	fHistPtVsUENSRec=0;
	fHistPtVsUEASRec=0;
	fHistPtVsUETSRec=0;



	for(Int_t i_pid=0; i_pid<11; ++i_pid){

		fHistPt[i_pid]=0;
		fHistPtRec[i_pid]=0;

	}
	for(Int_t i=0; i<3; ++i){// loop over mult estimators

		fMult[i] = 0;
		fMultRec[i] = 0;

		fNch[i]=0;
		fNchRec[i]=0;

	}


	// Default constructor (should not be used)
}

//______________________________________________________________________________

AliAnalysisTaskGenUeNchTS::AliAnalysisTaskGenUeNchTS(const char *name):
	AliAnalysisTaskSE(name),
	fMcEvent(0x0),
	fMcHandler(0x0),
	fStack(0),
	fGenerator("Pythia8"),
	fIndexLeadingGen(-1),
	fIndexLeadingRec(-1),
	fMinPtLeading(5.0),
	fSizeStep(0.1),
	//fNso_gen(3),
	//fNso_rec(3),
	fbinPerc0_gen(0),
	fbinPerc1_gen(0),
	fbinPerc2_gen(0),
	fbinPerc0_rec(0),
	fbinPerc1_rec(0),
	fbinPerc2_rec(0),
	fY(0.5),
	fHistEvt(0x0),
	fHistEta(0x0),
	fHistPart(0x0),
	fDphiNS(0x0),
	fDphiAS(0x0),
	fDphiTS(0x0),
	fMultTS(0x0),
	fDphiNSRec(0x0),
	fDphiASRec(0x0),
	fDphiTSRec(0x0),
	fMultTSRec(0x0),
	fListOfObjects(0)
{


	fHistPtVsUENS=0;
	fHistPtVsUEAS=0;
	fHistPtVsUETS=0;

	fHistPtVsUENSRec=0;
	fHistPtVsUEASRec=0;
	fHistPtVsUETSRec=0;


	for(Int_t i=0; i<3; ++i){
		fMult[i] = 0;
		fMultRec[i] = 0;

		fNch[i]=0;
		fNchRec[i]=0;

	}
	DefineInput( 0, TChain::Class());
	DefineOutput(1, TList::Class() ); // Basic output slot 
}


//_____________________________________________________________________________

AliAnalysisTaskGenUeNchTS::~AliAnalysisTaskGenUeNchTS(){
	// Destructor
	// histograms are in the output list and deleted when the output
	// list is deleted by the TSelector dtor
	if (fListOfObjects) { delete fListOfObjects; fListOfObjects=0x0; }
}

//______________________________________________________________________________

void AliAnalysisTaskGenUeNchTS::UserCreateOutputObjects(){

	// ### Analysis output
	fListOfObjects = new TList();
	fListOfObjects->SetOwner(kTRUE);

	// ### Create histograms
	fHistEvt = 0;	
	fHistEvt = new TH1I("fHistEvt","fHistEvt",2,0,2) ;
	fHistEvt->GetXaxis()->SetBinLabel(1,"All events");
	fHistEvt->GetXaxis()->SetBinLabel(2,"All particles");
	fHistEvt->Sumw2();
	fListOfObjects->Add(fHistEvt);

	fHistEta = 0;
	fHistEta = new TH1F("fHistEta","Eta Distr.; #eta; N_{part}", 200, -1., 1.);
	fListOfObjects->Add(fHistEta);

	InitHisto<TH1F>("fHistY", "Y Distr.", 200, -1., 1., "#it{y}", "N_{part}");
	InitHisto<TH1F>("fHistYRec", "Y Distr. rec", 200, -1., 1., "#it{y}", "N_{part}");

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
	fMultTS = new TH1D("fMultTS","",100,-0.5,99.5);
	fListOfObjects->Add(fMultTS);

	fDphiNSRec = 0;
	fDphiNSRec = new TH1D("hDphiNSRec","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(fDphiNSRec);

	fDphiASRec = 0;
	fDphiASRec = new TH1D("hDphiASRec","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(fDphiASRec);

	fDphiTSRec = 0;
	fDphiTSRec = new TH1D("hDphiTSRec","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(fDphiTSRec);

	fMultTSRec = 0;
	fMultTSRec = new TH1D("fMultTSRec","",100,-0.5,99.5);
	fListOfObjects->Add(fMultTSRec);


	Double_t TSBins[nTSBins+1]={0x0};
	for(Int_t i=0;i<nTSBins;++i){
		TSBins[i]=i*1.0-0.5;
	}
	TSBins[nTSBins]=99.5;

	
	fHistPtVsUENS = 0;
	fHistPtVsUENS=new TH2D("fHistPtVsUENS", "Generated #it{p}_{T} distribution NS",nTSBins, TSBins, NPtBins,PTBins);
	fListOfObjects->Add(fHistPtVsUENS);

	fHistPtVsUEAS=0;
	fHistPtVsUEAS=new TH2D("fHistPtVsUEAS", "Generated #it{p}_{T} distribution AS",nTSBins, TSBins, NPtBins,PTBins);
	fListOfObjects->Add(fHistPtVsUEAS);

	fHistPtVsUETS=0;
	fHistPtVsUETS=new TH2D("fHistPtVsUETS", "Generated #it{p}_{T} distribution TS",nTSBins, TSBins, NPtBins,PTBins);
	fListOfObjects->Add(fHistPtVsUETS);

	fHistPtVsUENSRec=0;
	fHistPtVsUENSRec=new TH2D("fHistPtVsUENSRec", "Rec #it{p}_{T} distribution NS",nTSBins, TSBins, NPtBins,PTBins);
	fListOfObjects->Add(fHistPtVsUENSRec);

	fHistPtVsUEASRec=0;
	fHistPtVsUEASRec=new TH2D("fHistPtVsUEASRec", "Rec #it{p}_{T} distribution AS",nTSBins, TSBins, NPtBins,PTBins);
	fListOfObjects->Add(fHistPtVsUEASRec);

	fHistPtVsUETSRec=0;
	fHistPtVsUETSRec=new TH2D("fHistPtVsUETSRec", "Rec #it{p}_{T} distribution TS",nTSBins, TSBins, NPtBins,PTBins);
	fListOfObjects->Add(fHistPtVsUETSRec);


	for(Int_t i=0; i<3; ++i){

		fMult[i] = 0;
		fMult[i] = new TH2D(Form("fMult_%s",Estimators[i]),Form("Selection bias %s; Nch; #eta",Estimators[i]),300,-0.5,299.5,100,-5,5); 
		fListOfObjects->Add(fMult[i]);

		fMultRec[i] = 0;
		fMultRec[i] = new TH2D(Form("fMultRec_%s",Estimators[i]),Form("Selection bias %s; Nch; #eta",Estimators[i]),300,-0.5,299.5,100,-5,5); 
		fListOfObjects->Add(fMultRec[i]);

		fNch[i]=0;
		fNch[i]=new TH1F(Form("fNch_%s",Estimators[i]),Form("Nch %s",Estimators[i]),300,-0.5,299.5);
		fListOfObjects->Add(fNch[i]);

		fNchRec[i]=0;
		fNchRec[i]=new TH1F(Form("fNchRec_%s",Estimators[i]),Form("Nch %s",Estimators[i]),300,-0.5,299.5);
		fListOfObjects->Add(fNchRec[i]);


	}


	// ### List of outputs
	PostData(1, fListOfObjects);

}

//______________________________________________________________________________

inline void AliAnalysisTaskGenUeNchTS::FillHisto(const char* objkey, Double_t x)
{
	TH1* hTmp = 0;
	hTmp = static_cast<TH1*>(fListOfObjects->FindObject(objkey));
	if(!hTmp){
		AliError(Form("Cannot find histogram: %s",objkey)) ;
		return;
	}
	hTmp->Fill(x);
}

inline void AliAnalysisTaskGenUeNchTS::FillHisto(const char* objkey, Double_t x, Double_t y)
{
	TH2* hTmp = 0;
	hTmp = static_cast<TH2*>(fListOfObjects->FindObject(objkey));
	if(!hTmp){
		AliError(Form("Cannot find histogram: %s",objkey)) ;
		return;
	}
	hTmp->Fill(x,y);
}

//______________________________________________________________________________

template <class T> T* AliAnalysisTaskGenUeNchTS::InitHisto(const char* hname, const char* htitle, Int_t nxbins, Double_t xmin, Double_t xmax, const char* xtitle, const char* ytitle)
{
	T* hTmp = 0;
	hTmp = new T(hname, htitle, nxbins, xmin, xmax);
	hTmp->GetXaxis()->SetTitle(xtitle);
	hTmp->GetYaxis()->SetTitle(ytitle);
	//	hTmp->SetMarkerStyle(kFullCircle);
	//	hTmp->Sumw2();
	fListOfObjects->Add(hTmp);

	return hTmp;
}
template <class T> T* AliAnalysisTaskGenUeNchTS::InitHisto(const char* hname, const char* htitle, Int_t nxbins, Double_t xmin, Double_t xmax, Int_t nybins, Double_t ymin, Double_t ymax, const char* xtitle, const char* ytitle)
{
	T* hTmp = 0;
	hTmp = new T(hname, htitle, nxbins, xmin, xmax, nybins, ymin, ymax);
	hTmp->GetXaxis()->SetTitle(xtitle);
	hTmp->GetYaxis()->SetTitle(ytitle);
	//	hTmp->Sumw2();
	fListOfObjects->Add(hTmp);

	return hTmp;
}

//______________________________________________________________________________

void AliAnalysisTaskGenUeNchTS::Init(){
	//
	fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());

	if(fGenerator == "Pythia6"){
		Double_t NchBin_gen0tmp[NchPercBin+1]={-0.5, 1.5, 3.5, 5.5, 7.5, 10.5, 15.5, 299.5};
		Double_t NchBin_gen1tmp[NchPercBin+1]={-0.5, 3.5, 6.5, 9.5, 13.5, 17.5, 25.5, 299.5};
		Double_t NchBin_gen2tmp[NchPercBin+1]={-0.5, 10.5, 18.5, 24.5, 31.5, 42.5, 57.5, 299.5};
		Double_t NchBin_rec0tmp[NchPercBin+1]={-0.5, 1.5, 2.5, 4.5, 5.5, 8.5, 11.5, 299.5};
		Double_t NchBin_rec1tmp[NchPercBin+1]={-0.5, 2.5, 4.5, 6.5, 9.5, 13.5, 18.5, 299.5};
		Double_t NchBin_rec2tmp[NchPercBin+1]={-0.5, 11.5, 19.5, 25.5, 32.5, 43.5, 57.5, 299.5};
		for(Int_t i=0;i<NchPercBin+1;++i){
			NchBin_gen0[i]=NchBin_gen0tmp[i];
			NchBin_gen1[i]=NchBin_gen1tmp[i];
			NchBin_gen2[i]=NchBin_gen2tmp[i];

			NchBin_rec0[i]=NchBin_rec0tmp[i];
			NchBin_rec1[i]=NchBin_rec1tmp[i];
			NchBin_rec2[i]=NchBin_rec2tmp[i];
		}
	}
	else if(fGenerator == "Epos"){
		Double_t NchBin_gen0tmp[NchPercBin+1]={-0.5, 1.5, 4.5, 5.5, 7.5, 10.5, 16.5, 299.5};
		Double_t NchBin_gen1tmp[NchPercBin+1]={-0.5, 3.5, 6.5, 9.5, 12.5, 17.5, 25.5, 299.5};
		Double_t NchBin_gen2tmp[NchPercBin+1]={-0.5, 12.5, 21.5, 27.5, 36.5, 49.5, 68.5, 299.5};
		Double_t NchBin_rec0tmp[NchPercBin+1]={-0.5, 1.5, 2.5, 3.5, 5.5, 7.5, 11.5, 299.5};
		Double_t NchBin_rec1tmp[NchPercBin+1]={-0.5, 2.5, 4.5, 6.5, 9.5, 12.5, 18.5, 299.5};
		Double_t NchBin_rec2tmp[NchPercBin+1]={-0.5, 13.5, 21.5, 28.5, 37.5, 50.5, 69.5, 299.5};
		for(Int_t i=0;i<NchPercBin+1;++i){
			NchBin_gen0[i]=NchBin_gen0tmp[i];
			NchBin_gen1[i]=NchBin_gen1tmp[i];
			NchBin_gen2[i]=NchBin_gen2tmp[i];
			NchBin_rec0[i]=NchBin_rec0tmp[i];
			NchBin_rec1[i]=NchBin_rec1tmp[i];
			NchBin_rec2[i]=NchBin_rec2tmp[i];
		}
	}
	else{
		Double_t NchBin_gen0tmp[NchPercBin+1]={-0.5, 1.5, 4.5, 5.5, 7.5, 10.5, 15.5, 299.5};
		Double_t NchBin_gen1tmp[NchPercBin+1]={-0.5, 3.5, 6.5, 9.5, 12.5, 17.5, 25.5, 299.5};
		Double_t NchBin_gen2tmp[NchPercBin+1]={-0.5, 12.5, 19.5, 26.5, 34.5, 46.5, 64.5, 299.5};
		Double_t NchBin_rec0tmp[NchPercBin+1]={-0.5, 1.5, 2.5, 3.5, 5.5, 7.5, 11.5, 299.5};
		Double_t NchBin_rec1tmp[NchPercBin+1]={-0.5, 2.5, 4.5, 6.5, 9.5, 12.5, 18.5, 299.5};
		Double_t NchBin_rec2tmp[NchPercBin+1]={-0.5, 12.5, 20.5, 26.5, 35.5, 47.5, 65.5, 299.5};
		for(Int_t i=0;i<NchPercBin+1;++i){
			NchBin_gen0[i]=NchBin_gen0tmp[i];
			NchBin_gen1[i]=NchBin_gen1tmp[i];
			NchBin_gen2[i]=NchBin_gen2tmp[i];
			NchBin_rec0[i]=NchBin_rec0tmp[i];
			NchBin_rec1[i]=NchBin_rec1tmp[i];
			NchBin_rec2[i]=NchBin_rec2tmp[i];
		}
	}
}

//______________________________________________________________________________

void AliAnalysisTaskGenUeNchTS::UserExec(Option_t *){

	// ### Initialize
	Init();

	// ### MC handler
	if(fMcHandler)
		fMcEvent = fMcHandler->MCEvent();
	else { if(fDebug > 1) printf("AliAnalysisTaskGenUeNchTS::Handler() fMcHandler = NULL\n"); return; }

	// ### MC event
	if( !fMcEvent ) { if(fDebug > 1) printf("AliAnalysisTaskGenUeNchTS::UserExec() fMcEvent = NULL \n"); return; }

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
	vector<Int_t> mult_estimators_gen;
	vector<Float_t> pt_so_gen;
	vector<Float_t> eta_so_gen;
	vector<Float_t> phi_so_gen;


	Int_t fNso_gen = -1;
	fNso_gen = GetMultipliciy( kFALSE, mult_estimators_gen, pt_so_gen, eta_so_gen, phi_so_gen );

	// PSEUDO REC LEVEL
	vector<Int_t> mult_estimators_rec;
	vector<Float_t> pt_so_rec;
	vector<Float_t> eta_so_rec;
	vector<Float_t> phi_so_rec;

	// RT analysis
	fIndexLeadingGen = -1;
	fIndexLeadingGen = GetIndexLeading(kFALSE);
	TParticle* mcPartLeadingGen         = 0x0;
	if(fIndexLeadingGen>=0){
		mcPartLeadingGen                    = (TParticle *)fMcEvent->Particle(fIndexLeadingGen);
		if(mcPartLeadingGen->Pt()>=fMinPtLeading){
			MakeRTAnalysis(kFALSE);
		}
	}

	fIndexLeadingRec = -1;
	fIndexLeadingRec = GetIndexLeading(kTRUE);
	TParticle* mcPartLeadingRec         = 0x0;
	if(fIndexLeadingRec>=0){
		mcPartLeadingRec                    = (TParticle *)fMcEvent->Particle(fIndexLeadingRec);
		if(mcPartLeadingRec->Pt()>=fMinPtLeading){
			MakeRTAnalysis(kTRUE);
		}
	}

	//MemInfo_t memInfo;
	//Int_t memUsage = 0;
	//gSystem->GetMemInfo(&memInfo);
	//memUsage = memInfo.fMemUsed;
	//cout<<"mem usage="<<memUsage<<endl;
	// ### Post data for all output slots
	mult_estimators_rec.clear();
	pt_so_rec.clear();
	eta_so_rec.clear();
	phi_so_rec.clear();

	mult_estimators_gen.clear();
	pt_so_gen.clear();
	eta_so_gen.clear();
	phi_so_gen.clear();


	PostData(1, fListOfObjects);

	return;
}

//______________________________________________________________________________

Bool_t AliAnalysisTaskGenUeNchTS::IsMCEventSelected(TObject* obj){

	Bool_t isSelected = kTRUE;

	AliMCEvent *event = 0x0;
	event = dynamic_cast<AliMCEvent*>(obj);
	if( !event ) 
		isSelected = kFALSE;

	if( isSelected ) 
		FillHisto("fHistEvt",0.5);

	return isSelected;
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskGenUeNchTS::GetIndexLeading(Bool_t fIsPseudoRec){

	Double_t ptleading = 0;
	Int_t index_leading = -1;
	Bool_t isPhysPrim = kFALSE;
	Double_t qPart = 0;
	Double_t etaPart = -10;
	Int_t pidCodeMC = 0;
	Int_t pPDG = -10;

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

		pPDG = TMath::Abs(mcPart->GetPdgCode());
		pidCodeMC = GetPidCode(pPDG);

		if(TMath::Abs(pidCodeMC)==5)
			continue;

		etaPart = mcPart -> Eta();
		if(fIsPseudoRec)
			if(!IsGoodTrack(pidCodeMC,mcPart->GetPdgCode(),mcPart->Pt())) continue;
		if(TMath::Abs(etaPart) > 0.8)continue;
		if(mcPart -> Pt()<1.0)continue;
		if(mcPart -> Pt()>ptleading){
			ptleading = mcPart->Pt();
			index_leading = ipart;
		}


	} // particle loop
	cout<<"\n";
	return index_leading;
}

Int_t AliAnalysisTaskGenUeNchTS::GetMultipliciy(Bool_t fIsPseudoRec, vector<Int_t> &multArray, vector<Float_t> &ptArray,  vector<Float_t> &etaArray, vector<Float_t> &phiArray){


	Int_t mult_so=0;
	multArray.clear();
	ptArray.clear();
	etaArray.clear();
	phiArray.clear();

	Bool_t isPhysPrim = kFALSE;

	Int_t mult_Eta5   = 0;
	Int_t mult_Eta8   = 0;
	Int_t mult_Eta1   = 0;
	Int_t mult_VZEROM = 0;
	Double_t qPart = 0;
	Double_t etaPart = -10;


	Int_t pidCodeMC = 0;
	Int_t pPDG = -10;

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

		pPDG = TMath::Abs(mcPart->GetPdgCode());
		pidCodeMC = GetPidCode(pPDG);

		etaPart = mcPart -> Eta();
		if( (2.8 < etaPart && etaPart < 5.1) || (-3.7 < etaPart && etaPart <-1.7) ) mult_VZEROM++;

		//	if(fIsPseudoRec)
		//	if(!IsGoodTrack(pidCodeMC,mcPart->GetPdgCode(),mcPart->Pt())) continue;

		if( TMath::Abs(etaPart) < 0.5 ) mult_Eta5++;
		if( TMath::Abs(etaPart) < 0.8 ){ 
			mult_Eta8++;

			if(mcPart -> Pt()>0.15){

				ptArray.push_back(mcPart->Pt());
				etaArray.push_back(mcPart->Eta());
				phiArray.push_back(mcPart->Phi());
				mult_so++;
			}

		}
		if( TMath::Abs(etaPart) < 1 )
			if(mcPart -> Pt()>0)
				mult_Eta1++;// for INEL>0n

	} // particle loop

	multArray.push_back(mult_Eta5);
	multArray.push_back(mult_Eta8);
	multArray.push_back(mult_VZEROM);
	multArray.push_back(mult_Eta1);

	return mult_so;
}

//______________________________________________________________________________
void AliAnalysisTaskGenUeNchTS::MakeRTAnalysis(Bool_t fIsPseudoRec){

	// Properties leading particle
	TParticle* mcPartTmp         = 0x0;
	if(fIsPseudoRec)
		mcPartTmp                    = (TParticle *)fMcEvent->Particle(fIndexLeadingRec);
	else
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
	// ### particle loop
	for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {

		if(fIsPseudoRec){
			if(ipart == fIndexLeadingRec)continue;
		}
		else{ 
			if(ipart == fIndexLeadingGen)continue;
		}

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
		pPDG = TMath::Abs(mcPart->GetPdgCode());
		pidCodeMC = GetPidCode(pPDG);
		etaPart = mcPart -> Eta();
		if(TMath::Abs(etaPart)>0.8)continue;
		ipt = mcPart->Pt();
		//if(fIsPseudoRec)
		//	if(!IsGoodTrack(pidCodeMC,mcPart->GetPdgCode(),ipt)) continue;
		if(ipt<0.15)continue;
		phiPart = mcPart -> Phi();
		Double_t DPhi = DeltaPhi(phiL,phiPart);

		if(fIsPseudoRec){
			if(TMath::Abs(DPhi)<pi/3.0){
				fDphiNSRec->Fill(DPhi);
			}
			// away side
			else if(TMath::Abs(DPhi-pi)<pi/3.0){
				fDphiASRec->Fill(DPhi);
			}
			// transverse side
			else{
				multTS++;
				fDphiTSRec->Fill(DPhi);
			}
		}
		else{
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


	}

	if(fIsPseudoRec)
		fMultTSRec->Fill(multTS);
	else
		fMultTS->Fill(multTS);


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

		if(fIsPseudoRec){
			if(ipart == fIndexLeadingRec)continue;
		}
		else{
			if(ipart == fIndexLeadingGen)continue;
		}


		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)fMcEvent->Particle(ipart);
		if (!mcPart) continue;

		if(!mcPart->GetPDG())continue;
		isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
		qPart = mcPart->GetPDG()->Charge()/3.;
		// only primary charged particles
		pPDG = TMath::Abs(mcPart->GetPdgCode());
		pidCodeMC = GetPidCode(pPDG);
		Bool_t isSelectedPart = kTRUE;
		for(Int_t i=0; i<11; ++i) 
			if( pidCodeMC == i ) 
				isSelectedPart = kFALSE;
		if ( isSelectedPart ) continue;
		ipt = mcPart->Pt();
		if(ipt<0.15)continue;
		phiPart = mcPart -> Phi();
		Double_t DPhi = DeltaPhi(phiL,phiPart);
		if(TMath::Abs(mcPart->Eta()) >= 0.8)continue;

		for(Int_t i=0; i<11; ++i)
		{
			if( pidCodeMC == i )
			{
				if( IsPrimary[i] == kTRUE && isPhysPrim == kFALSE ) 
					continue;

				if(fIsPseudoRec){
					if(!IsGoodTrack(pidCodeMC,mcPart->GetPdgCode(),ipt)) continue;
					if(TMath::Abs(DPhi)<pi/3.0){// near side
						fHistPtVsNchNSRec[i]->Fill(1.0*multTS,ipt);
					}
					else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
						fHistPtVsNchASRec[i]->Fill(1.0*multTS,ipt);
					}
					else{// transverse side
						fHistPtVsNchTSRec[i]->Fill(1.0*multTS,ipt);
					}


				}else{

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
		}

		if(fIsPseudoRec){
		  if(!IsGoodTrack(pidCodeMC,mcPart->GetPdgCode(),ipt)) continue;
		  if(TMath::Abs(DPhi)<pi/3.0){// near side
		    fHistPtVsUENSRec->Fill(1.0*multTS,ipt);
		  }
		  else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
		    fHistPtVsUEASRec->Fill(1.0*multTS,ipt);
		  }
		  else{// transverse side
		    fHistPtVsUETSRec->Fill(1.0*multTS,ipt);
		  }


		}else{

		  if(TMath::Abs(DPhi)<pi/3.0){// near side
		    fHistPtVsUENS->Fill(1.0*multTS,ipt);
		  }       
		  else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
		    fHistPtVsUEAS->Fill(1.0*multTS,ipt);
		  }       
		  else{// transverse side
		    fHistPtVsUETS->Fill(1.0*multTS,ipt);
		  }       
		}

	} // particle loop

}


void AliAnalysisTaskGenUeNchTS::ParticleSel(Bool_t fIsPseudoRec, const vector<Int_t> &mult){


	Int_t pidCodeMC = 0;
	Double_t ipt = 0.;

	Bool_t isPhysPrim = kFALSE;
	Double_t qPart = 0;
	Int_t pPDG = -10;
	Double_t y = -10;
	// ### particle loop
	for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {

		TParticle* mcPart         = 0x0;
		mcPart                    = (TParticle *)fMcEvent->Particle(ipart);
		if (!mcPart) continue;

		FillHisto("fHistEvt",1.5);
		if(!mcPart->GetPDG())continue;
		isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
		qPart = mcPart->GetPDG()->Charge()/3.;
		// only primary charged particles
		if( isPhysPrim ){
			if(TMath::Abs(qPart)>0.001)
				for(Int_t i=0;i<3;++i){
					if(fIsPseudoRec)
						fMultRec[i]->Fill(1.0*mult[i],mcPart->Eta());
					//FillHisto(Form("fMultRec_%s",estimators[i]),1.0*mult[i],mcPart->Eta());
					else
						fMult[i]->Fill(1.0*mult[i],mcPart->Eta());
					//	FillHisto(Form("fMult_%s",estimators[i]),1.0*mult[i],mcPart->Eta());
				}
		}

		pPDG = TMath::Abs(mcPart->GetPdgCode());
		pidCodeMC = GetPidCode(pPDG);
		Bool_t isSelectedPart = kTRUE;
		for(Int_t i=0; i<11; ++i) 
			if( pidCodeMC == i ) 
				isSelectedPart = kFALSE;
		if ( isSelectedPart ) continue;

		fHistEta->Fill(mcPart->Eta());

		if (!(TMath::Abs(mcPart->Energy()-mcPart->Pz())>0.)) continue;
		Double_t myY = (mcPart->Energy()+mcPart->Pz())/(mcPart->Energy()-mcPart->Pz());
		if( myY <= 0 ) continue;
		y = 0.5*TMath::Log(myY);
		//y = mcPart->Y(); 
		ipt = mcPart->Pt();


		for(Int_t i=0; i<11; ++i)
		{
			if( pidCodeMC == i && TMath::Abs(y) < fY)
			{
				if( IsPrimary[i] == kTRUE && isPhysPrim == kFALSE ) 
					continue;

				if(fIsPseudoRec){
					if(!IsGoodTrack(pidCodeMC,mcPart->GetPdgCode(),mcPart->Pt())) continue;
					if(i==0) FillHisto("fHistYRec",y);
					fHistPtRec[i]->Fill(ipt);
					


				}else{
					if(i==0) FillHisto("fHistY",y);
					fHistPt[i]->Fill(ipt);
				

				}

			}
		}

	} // particle loop

}


//______________________________________________________________________________

void AliAnalysisTaskGenUeNchTS::Terminate(Option_t*){

	fListOfObjects = dynamic_cast<TList*> (GetOutputData(1));
	if (!fListOfObjects) { Printf("ERROR: Output list not available"); return; }

	return;
}

//_____________________________________________________________________________

Int_t AliAnalysisTaskGenUeNchTS::GetPidCode(Int_t pdgCode) const  {

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

// Bool_t AliAnalysisTaskGenUeNchTS::IsGoodTrack(Int_t pid, Int_t pdgcode, Double_t pt){
// 	Double_t efficiency;
// 	if(pid==0)
// 		efficiency = ch_Eff->Eval(pt);
// 	else if(pid==1){
// 		if(pdgcode>0)
// 			efficiency = k_pos_Eff->Eval(pt);
// 		else
// 			efficiency = k_neg_Eff->Eval(pt);
// 	}
// 	else if(pid==2){
// 		if(pdgcode>0)
// 			efficiency = p_pos_Eff->Eval(pt);
// 		else
// 			efficiency = p_neg_Eff->Eval(pt);
// 	}
// 	else if(pid==3)
// 		efficiency = k0s_Eff->Eval(pt);
// 	else if(pid==4){
// 		if(pdgcode>0)
// 			efficiency = la_Eff->Eval(pt);
// 		else
// 			efficiency = labar_Eff->Eval(pt);
// 	}
// 	else if(pid==5)
// 		efficiency = xi_Eff->Eval(pt);
// 	else if(pid==7)
// 		efficiency = phi_Eff->Eval(pt);
// 	else if(pid==8)
// 		efficiency = k0star_Eff->Eval(pt);
// 	else
// 		efficiency = ch_Eff->Eval(pt);


// 	if(gRandom->Uniform(1.0)<efficiency)
// 		return kTRUE;
// 	else
// 		return kFALSE;
// }

Int_t AliAnalysisTaskGenUeNchTS::GetMultBin(Bool_t fIsPseudoRec,Int_t mult_int, Int_t mult_select){

	if(fIsPseudoRec){
		for(Int_t j=0;j<NchPercBin;++j){
			if(mult_select==0){
				if(mult_int>=NchBin_rec0[j] && mult_int<NchBin_rec0[j+1])
					return j;
			}
			if(mult_select==1){
				if(mult_int>=NchBin_rec1[j] && mult_int<NchBin_rec1[j+1])
					return j;
			}
			if(mult_select==2){
				if(mult_int>=NchBin_rec2[j] && mult_int<NchBin_rec2[j+1])
					return j;
			}
		}
	}
	else{
		for(Int_t j=0;j<NchPercBin;++j){
			if(mult_select==0){
				if(mult_int>=NchBin_gen0[j] && mult_int<NchBin_gen0[j+1])
					return j;
			}
			if(mult_select==1){
				if(mult_int>=NchBin_gen1[j] && mult_int<NchBin_gen1[j+1])
					return j;
			}
			if(mult_select==2){
				if(mult_int>=NchBin_gen2[j] && mult_int<NchBin_gen2[j+1])
					return j;
			}
		}
	}

	return -1;
}
Double_t AliAnalysisTaskGenUeNchTS::DeltaPhi(Double_t phia, Double_t phib,
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
