
/*

   This macro produces: leading particle dE/dx vs leading particle momentum vs Multiplicity (near side, away side and transverse side) 
   Antonio Ortiz, ICN-UNAM
   Please report bugs to: aortizve@cern.ch / antonio.ortiz@nucleares.unam.mx 
   First version: 20/02/2019
   Calibration factors of dE/dx were provided by Omar Vazquez (Lund University)

   Update 05/03/2019
   - dE/dx range was extended
   - Selction on momentum for pion MIP was added. dE/dx vs eta histos before recalibration were added 
   - Sum pT density has been added

 */

#include "AliAnalysisTaskUePid.h"

// ROOT includes
#include <TList.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TFile.h>


// AliRoot includes
#include <AliAnalysisManager.h>
#include <AliAnalysisFilter.h>
#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliEventCuts.h>
#include <AliESDVertex.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrackCuts.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>
#include <TTreeStream.h>
#include <AliHeader.h>
#include <AliAnalysisUtils.h>
#include <AliMultiplicity.h>
#include <AliMultSelection.h>
#include <AliAODInputHandler.h> 
#include <AliAODHandler.h> 
#include <AliAODVertex.h>
#include <AliAODTrack.h> 
#include <AliAODPid.h> 
#include <AliDataFile.h>

#include <iostream>
using namespace std;


const Char_t * legendsEta[5] = {"|#eta|<0.8","|#eta|<0.2","0.2<|#eta|<0.4","0.4<|#eta|<0.6","0.6<|#eta|<0.8"};
Float_t Magf                    = 1;
const Int_t nPeriod               = 2;
Int_t index_data = 0;
const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;
const Double_t aPos[nPeriod]      = {49.9799,49.9659};// period l, k...
const Double_t bPos[nPeriod]      = {2.99619,2.91366};
const Double_t cPos[nPeriod]      = {-45.718,-45.5994};
const Double_t dPos[nPeriod]      = {290.013,290.042};
const Double_t ePos[nPeriod]      = {-1018.42,-1014.49};
const Double_t fPos[nPeriod]      = {1948.68,1931.84};
const Double_t gPos[nPeriod]      = {-1864.06,-1839.36};
const Double_t hPos[nPeriod]      = {692.752,680.421};

const Double_t aNeg[nPeriod]      = {50.078,50.046};
const Double_t bNeg[nPeriod]      = {6.67199,6.79992};
const Double_t cNeg[nPeriod]      = {103.662,109.86};
const Double_t dNeg[nPeriod]      = {611.034,668.241};
const Double_t eNeg[nPeriod]      = {1695.63,1916.44};
const Double_t fNeg[nPeriod]      = {2395.88,2815.04};
const Double_t gNeg[nPeriod]      = {1669.22,2057.21};
const Double_t hNeg[nPeriod]      = {455.362,595.391};



ClassImp(AliAnalysisTaskUePid)

	//_____________________________________________________________________________
	AliAnalysisTaskUePid::AliAnalysisTaskUePid():
		AliAnalysisTaskSE(),
		fESD(0x0),
		fdata_set("k"),
		fEventCuts(0x0),
		fTrackFilter(0x0),
		fAnalysisType("ESD"),
		fcutLow(0x0),
		fcutHigh(0x0),
		fDeDxMIPMin(35),
		fDeDxMIPMax(95),
		fEtaCalibrationNeg(0x0),
		fEtaCalibrationPos(0x0),
		fNcl(70),
		fListOfObjects(0x0),
		fHistEventCounter(0x0),
		fEvents(0x0),
		hPhi(0x0),
		hMIPVsEtaBefore(0x0),
		pMIPVsEtaBefore(0x0),		
		hMIPVsEta(0x0),
		pMIPVsEta(0x0),
		hpT(0x0),
		hDphi(0x0),
		hDphiNNS(0x0),
		hDphiNS(0x0),
		hDphiAS(0x0),
		hDphiTS(0x0)

{
	// Default constructor (should not be used)
	for(Int_t i_eta = 0; i_eta < 5; ++i_eta){
		hPtL[i_eta] = 0;
		hDeDxL[i_eta] = 0;
		hEtaL[i_eta] = 0;
		hPhiL[i_eta] = 0;

		pNNS[i_eta] = 0;
		pNS[i_eta] = 0;
		pAS[i_eta] = 0;
		pTS[i_eta] = 0;

		hNNS[i_eta] = 0;
		hNS[i_eta] = 0;
		hAS[i_eta] = 0;
		hTS[i_eta] = 0;

		h3DNNS[i_eta] = 0;
		h3DNS[i_eta] = 0;
		h3DAS[i_eta] = 0;
		h3DTS[i_eta] = 0;

		pNNS2[i_eta] = 0;
		pNS2[i_eta] = 0;
		pAS2[i_eta] = 0;
		pTS2[i_eta] = 0;

		hNNS2[i_eta] = 0;
		hNS2[i_eta] = 0;
		hAS2[i_eta] = 0;
		hTS2[i_eta] = 0;

		h3DNNS2[i_eta] = 0;
		h3DNS2[i_eta] = 0;
		h3DAS2[i_eta] = 0;
		h3DTS2[i_eta] = 0;

		hptVSp[i_eta] = 0;
	}


}

//______________________________________________________________________________
AliAnalysisTaskUePid::AliAnalysisTaskUePid(const char *name):
	AliAnalysisTaskSE(name),
	fESD(0x0),
	fdata_set("k"),
	fEventCuts(0x0),
	fTrackFilter(0x0),
	fAnalysisType("ESD"),
	fcutLow(0x0),
	fcutHigh(0x0),
	fDeDxMIPMin(35),
	fDeDxMIPMax(95),
	fEtaCalibrationNeg(0x0),
	fEtaCalibrationPos(0x0),
	fNcl(70),
	fListOfObjects(0x0),
	fHistEventCounter(0x0), 
	fEvents(0x0),
	hPhi(0x0),
	hMIPVsEtaBefore(0x0),
	pMIPVsEtaBefore(0x0),	
	hMIPVsEta(0x0),
	pMIPVsEta(0x0),
	hpT(0x0),
	hDphi(0x0),
	hDphiNNS(0x0),
	hDphiNS(0x0),
	hDphiAS(0x0),
	hDphiTS(0x0)

{
	// Output slot #1 writes into a TList
	for(Int_t i_eta = 0; i_eta < 5; ++i_eta){
		hPtL[i_eta] = 0;
		hDeDxL[i_eta] = 0;
		hEtaL[i_eta] = 0;
		hPhiL[i_eta] = 0;

		pNNS[i_eta] = 0;
		pNS[i_eta] = 0;
		pAS[i_eta] = 0;
		pTS[i_eta] = 0;

		hNNS[i_eta] = 0;
		hNS[i_eta] = 0;
		hAS[i_eta] = 0;
		hTS[i_eta] = 0;

		h3DNNS[i_eta] = 0;
		h3DNS[i_eta] = 0;
		h3DAS[i_eta] = 0;
		h3DTS[i_eta] = 0;

		pNNS2[i_eta] = 0;
		pNS2[i_eta] = 0;
		pAS2[i_eta] = 0;
		pTS2[i_eta] = 0;

		hNNS2[i_eta] = 0;
		hNS2[i_eta] = 0;
		hAS2[i_eta] = 0;
		hTS2[i_eta] = 0;

		h3DNNS2[i_eta] = 0;
		h3DNS2[i_eta] = 0;
		h3DAS2[i_eta] = 0;
		h3DTS2[i_eta] = 0;

		hptVSp[i_eta] = 0;
	}

	DefineOutput(1, TList::Class());
}

//________________________________________________________________________

void AliAnalysisTaskUePid::Exit(const char *msg) {

	Printf("%s", msg);
	return;
}


//_____________________________________________________________________________
AliAnalysisTaskUePid::~AliAnalysisTaskUePid()
{
	// Destructor
	// histograms are in the output list and deleted when the output
	// list is deleted by the TSelector dtor
	if (fListOfObjects && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
		delete fListOfObjects;
		fListOfObjects = 0x0;
	}

}

//______________________________________________________________________________
void AliAnalysisTaskUePid::UserCreateOutputObjects()
{

	fcutLow = new TF1("StandardPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 50);
	fcutHigh = new TF1("StandardPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 50);

	fEtaCalibrationNeg = new TF1("fDeDxVsEtaNeg", "pol7", -1.0, 0.0);
	fEtaCalibrationPos    = new TF1("fDeDxVsEtaPos", "pol7", 0.0, 1.0);

	const Int_t nPtBins = 52;
	Double_t xBins[nPtBins+1] = {
		0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
		1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
		2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
		4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
		11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
		26.0, 28.0, 30.0 };

	const Int_t nNchBins = 21;
	Double_t NchBins[nNchBins+1]={
		-0.5, 0.5 , 1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 , 7.5 , 8.5 , 9.5 ,
		10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5
	};

	const Int_t nNchBins2 = 60;
	Double_t NchBins2[nNchBins2+1]={
		0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
		1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
		2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
		4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
		11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
		26.0, 28.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 80.0, 100.0, 200.0
	};



	const Int_t nDeDxBins = 60;
	Double_t DeDxBins[nDeDxBins+1]=
	{
		35,36,37,38,39,40,41,42,43,44,
		45,46,47,48,49,50,51,52,53,54,
		55,56,57,58,59,60,61,62,63,64,
		65,66,67,68,69,70,71,72,73,74,
		75,76,77,78,79,80,81,82,83,84,
		85,86,87,88,89,90,91,92,93,94,95
	};
	// This method is called once per worker node
	// Here we define the output: histograms and debug tree if requested 

	// Definition of trackcuts
	if(!fTrackFilter){	
		fTrackFilter = new AliAnalysisFilter("trackFilter2015");
		SetTrackCuts(fTrackFilter);
	}

	OpenFile(1);
	fListOfObjects = new TList();
	fListOfObjects->SetOwner();

	//
	// Histograms
	//  
	if(! fHistEventCounter ){
		fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",10,0,10);

		fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
		fHistEventCounter->GetXaxis()->SetBinLabel(2, "Trigger");//NotinVertexcut");
		fHistEventCounter->GetXaxis()->SetBinLabel(3, "Physics Selection"); //CINT7-B-NOPF-CENT");
		fListOfObjects->Add(fHistEventCounter);
	}

	fEvents = new TH1I("fEvents","Number of analyzed events; Events; Counts", 1, 0, 1);
	fListOfObjects->Add(fEvents);

	hPhi = new TH2D("histPhi", ";pt; #phi'", nPtBins,xBins, 90, -0.05, 0.4);
	fListOfObjects->Add(hPhi);

	hMIPVsEtaBefore = new TH2D("hMIPVsEtaBefore","; #eta; dE/dx_{MIP, primary tracks}",50,-0.8,0.8,fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
	fListOfObjects->Add(hMIPVsEtaBefore);

	pMIPVsEtaBefore = new TProfile("pMIPVsEtaBefore","; #eta; #LT dE/dx #GT_{MIP, primary tracks}",50,-0.8,0.8, fDeDxMIPMin, fDeDxMIPMax);
	fListOfObjects->Add(pMIPVsEtaBefore);


	hMIPVsEta = new TH2D("hMIPVsEta","; #eta; dE/dx_{MIP, primary tracks}",50,-0.8,0.8,fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
	fListOfObjects->Add(hMIPVsEta);

	pMIPVsEta = new TProfile("pMIPVsEta","; #eta; #LT dE/dx #GT_{MIP, primary tracks}",50,-0.8,0.8, fDeDxMIPMin, fDeDxMIPMax);
	fListOfObjects->Add(pMIPVsEta);

	hpT = new TH1D("hpT","",nPtBins,xBins);
	fListOfObjects->Add(hpT);

	hDphi = 0;
	hDphi = new TH1D("hDphi","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(hDphi);

	hDphiNNS = 0;
	hDphiNNS = new TH1D("hDphiSNN","new near side",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(hDphiNNS);

	hDphiNS = 0;
	hDphiNS = new TH1D("hDphiSN","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(hDphiNS);

	hDphiAS = 0;
	hDphiAS = new TH1D("hDphiAS","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(hDphiAS);

	hDphiTS = 0;
	hDphiTS = new TH1D("hDphiTS","",2*64,-2*TMath::Pi(),2*TMath::Pi());
	fListOfObjects->Add(hDphiTS);

	for(Int_t i_eta=0; i_eta<5; ++i_eta){
		hPtL[i_eta] = 0;
		hPtL[i_eta] = new TH1D(Form("hPtL_etaInt%d",i_eta),Form("%s; #it{p}_{T}^{leading} (GeV/#it{c});counts",legendsEta[i_eta]),nPtBins,xBins);
		fListOfObjects->Add( hPtL[i_eta] );

		hDeDxL[i_eta] = 0;
		hDeDxL[i_eta] = new TH2D(Form("hDeDxL_etaInt%d",i_eta),Form("%s; TPC-d#it{E}/d#it{x}^{leading};counts",legendsEta[i_eta]),nPtBins,xBins,nDeDxBins,DeDxBins);
		fListOfObjects->Add( hDeDxL[i_eta] );

		hEtaL[i_eta] = 0;
		hEtaL[i_eta] = new TH1D(Form("hEtaL_etaInt%d",i_eta),Form("%s; #eta^{leading};counts",legendsEta[i_eta]),20,-1,1);
		fListOfObjects->Add(hEtaL[i_eta]);

		hPhiL[i_eta] = 0;
		hPhiL[i_eta] = new TH1D(Form("hPhiL_etaInt%d",i_eta),Form("%s; #phi^{leading} (rad);counts",legendsEta[i_eta]),64,0,2.0*TMath::Pi());
		fListOfObjects->Add(hPhiL[i_eta]);

		// Number density 


		pNNS[i_eta] = 0;
		pNNS[i_eta] = new TProfile(Form("pNNS%d",i_eta),Form("R<0.5 (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,0,20);
		fListOfObjects->Add(pNNS[i_eta]);

		pNS[i_eta] = 0;
		pNS[i_eta] = new TProfile(Form("pNS%d",i_eta),Form("Near side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,0,20);
		fListOfObjects->Add(pNS[i_eta]);

		pAS[i_eta] = 0;
		pAS[i_eta] = new TProfile(Form("pAS%d",i_eta),Form("Away side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,0,20);
		fListOfObjects->Add(pAS[i_eta]);

		pTS[i_eta] = 0;
		pTS[i_eta] = new TProfile(Form("pTS%d",i_eta),Form("Transverse side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,0,20);
		fListOfObjects->Add(pTS[i_eta]);

		hNNS[i_eta] = 0;
		hNNS[i_eta] = new TH2D(Form("hNNS%d",i_eta),Form("R<0.5 (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,nNchBins,NchBins);
		fListOfObjects->Add(hNNS[i_eta]);

		hNS[i_eta] = 0;
		hNS[i_eta] = new TH2D(Form("hNS%d",i_eta),Form("Near side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,nNchBins,NchBins);
		fListOfObjects->Add(hNS[i_eta]);

		hAS[i_eta] = 0;
		hAS[i_eta] = new TH2D(Form("hAS%d",i_eta),Form("Away side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,nNchBins,NchBins);
		fListOfObjects->Add(hAS[i_eta]);

		hTS[i_eta] = 0;
		hTS[i_eta] = new TH2D(Form("hTS%d",i_eta),Form("Transverse side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,nNchBins,NchBins);
		fListOfObjects->Add(hTS[i_eta]);

		h3DNNS[i_eta] = 0;
		h3DNNS[i_eta] = new TH3D(Form("h3DNNS%d",i_eta),Form("#it{R}<0.5 (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c}); TPC-d#it{E}/d#it{x}^{leading}",legendsEta[i_eta]),nPtBins,xBins,nNchBins,NchBins,nDeDxBins,DeDxBins);
		fListOfObjects->Add(h3DNNS[i_eta]);

		h3DNS[i_eta] = 0;
		h3DNS[i_eta] = new TH3D(Form("h3DNS%d",i_eta),Form("Near side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c}); TPC-d#it{E}/d#it{x}^{leading}",legendsEta[i_eta]),nPtBins,xBins,nNchBins,NchBins,nDeDxBins,DeDxBins);
		fListOfObjects->Add(h3DNS[i_eta]);

		h3DAS[i_eta] = 0;
		h3DAS[i_eta] = new TH3D(Form("h3DAS%d",i_eta),Form("Away side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c}); TPC-d#it{E}/d#it{x}^{leading}",legendsEta[i_eta]),nPtBins,xBins,nNchBins,NchBins,nDeDxBins,DeDxBins);
		fListOfObjects->Add(h3DAS[i_eta]);

		h3DTS[i_eta] = 0;
		h3DTS[i_eta] = new TH3D(Form("h3DTS%d",i_eta),Form("Transverse side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c}); TPC-d#it{E}/d#it{x}^{leading}",legendsEta[i_eta]),nPtBins,xBins,nNchBins,NchBins,nDeDxBins,DeDxBins);
		fListOfObjects->Add(h3DTS[i_eta]);

		// Sum pT density [TODO: chanhe the sum pT intervals] 

		pNNS2[i_eta] = 0;
		pNNS2[i_eta] = new TProfile(Form("pNNS_pt_%d",i_eta),Form("R<0.5 (%s); #it{p}^{leading} (GeV/#it{c}); #sum #it{p}_{T} (#it{p}_{T}>0.5 GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,0,20);
		fListOfObjects->Add(pNNS2[i_eta]);

		pNS2[i_eta] = 0;
		pNS2[i_eta] = new TProfile(Form("pNS_pt_%d",i_eta),Form("Near side (%s); #it{p}^{leading} (GeV/#it{c}); #sum #it{p}_{T} (#it{p}_{T}>0.5 GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,0,20);
		fListOfObjects->Add(pNS2[i_eta]);

		pAS2[i_eta] = 0;
		pAS2[i_eta] = new TProfile(Form("pAS_pt_%d",i_eta),Form("Away side (%s); #it{p}^{leading} (GeV/#it{c}); #sum #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,0,20);
		fListOfObjects->Add(pAS2[i_eta]);

		pTS2[i_eta] = 0;
		pTS2[i_eta] = new TProfile(Form("pTS_pt_%d",i_eta),Form("Transverse side (%s); #it{p}^{leading} (GeV/#it{c}); #sum #it{p}_{T} (#it{p}_{T}>0.5 GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,0,20);
		fListOfObjects->Add(pTS2[i_eta]);

		hNNS2[i_eta] = 0;
		hNNS2[i_eta] = new TH2D(Form("hNNS_pt_%d",i_eta),Form("R<0.5 (%s); #it{p}^{leading} (GeV/#it{c}); #sum #it{p}_{T} (#it{p}_{T}>0.5 GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,nNchBins2,NchBins2);
		fListOfObjects->Add(hNNS2[i_eta]);

		hNS2[i_eta] = 0;
		hNS2[i_eta] = new TH2D(Form("hNS_pt_%d",i_eta),Form("Near side (%s); #it{p}^{leading} (GeV/#it{c}); #sum #it{p}_{T} (#it{p}_{T}>0.5 GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,nNchBins2,NchBins2);
		fListOfObjects->Add(hNS2[i_eta]);

		hAS2[i_eta] = 0;
		hAS2[i_eta] = new TH2D(Form("hAS_pt_%d",i_eta),Form("Away side (%s); #it{p}^{leading} (GeV/#it{c}); #sum #it{p}_{T} (#it{p}_{T}>0.5 GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,nNchBins2,NchBins2);
		fListOfObjects->Add(hAS2[i_eta]);

		hTS2[i_eta] = 0;
		hTS2[i_eta] = new TH2D(Form("hTS_pt_%d",i_eta),Form("Transverse side (%s); #it{p}^{leading} (GeV/#it{c}); #sum #it{p}_{T} (#it{p}_{T}>0.5 GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,nNchBins2,NchBins2);
		fListOfObjects->Add(hTS2[i_eta]);

		h3DNNS2[i_eta] = 0;
		h3DNNS2[i_eta] = new TH3D(Form("h3DNNS_pt_%d",i_eta),Form("#it{R}<0.5 (%s); #it{p}^{leading} (GeV/#it{c}); #sum #it{p}_{T} (#it{p}_{T}>0.5 GeV/#it{c}); TPC-d#it{E}/d#it{x}^{leading}",legendsEta[i_eta]),nPtBins,xBins,nNchBins2,NchBins2,nDeDxBins,DeDxBins);
		fListOfObjects->Add(h3DNNS2[i_eta]);

		h3DNS2[i_eta] = 0;
		h3DNS2[i_eta] = new TH3D(Form("h3DNS_pt_%d",i_eta),Form("Near side (%s); #it{p}^{leading} (GeV/#it{c}); #sum #it{p}_{T} (#it{p}_{T}>0.5 GeV/#it{c}); TPC-d#it{E}/d#it{x}^{leading}",legendsEta[i_eta]),nPtBins,xBins,nNchBins2,NchBins2,nDeDxBins,DeDxBins);
		fListOfObjects->Add(h3DNS2[i_eta]);

		h3DAS2[i_eta] = 0;
		h3DAS2[i_eta] = new TH3D(Form("h3DAS_pt_%d",i_eta),Form("Away side (%s); #it{p}^{leading} (GeV/#it{c}); #sum #it{p}_{T} (#it{p}_{T}>0.5 GeV/#it{c}); TPC-d#it{E}/d#it{x}^{leading}",legendsEta[i_eta]),nPtBins,xBins,nNchBins2,NchBins2,nDeDxBins,DeDxBins);
		fListOfObjects->Add(h3DAS2[i_eta]);

		h3DTS2[i_eta] = 0;
		h3DTS2[i_eta] = new TH3D(Form("h3DTS_pt_%d",i_eta),Form("Transverse side (%s); #it{p}^{leading} (GeV/#it{c}); #sum #it{p}_{T} (#it{p}_{T}>0.5 GeV/#it{c}); TPC-d#it{E}/d#it{x}^{leading}",legendsEta[i_eta]),nPtBins,xBins,nNchBins2,NchBins2,nDeDxBins,DeDxBins);
		fListOfObjects->Add(h3DTS2[i_eta]);

		// momentum vs pT


		hptVSp[i_eta] = 0;
		hptVSp[i_eta] = new TH2D(Form("hptVSp%d",i_eta),Form("%s; #it{p}_{T}^{leading} (GeV/#it{c}); #it{p}^{leading} (GeV/#it{c})",legendsEta[i_eta]),nPtBins,xBins,nPtBins,xBins);
		fListOfObjects->Add(hptVSp[i_eta]);


	}

	fEventCuts.AddQAplotsToList(fListOfObjects);

	PostData(1, fListOfObjects);

}

//______________________________________________________________________________
void AliAnalysisTaskUePid::UserExec(Option_t *)
{

	// -----------------------------------------------------
	//			 InputEvent
	// -----------------------------------------------------

	AliVEvent *event = InputEvent();
	if (!event) {
		Error("UserExec", "Could not retrieve event");
		return;
	}

	// -----------------------------------------------------
	//			 E S D
	// -----------------------------------------------------

	if (fAnalysisType == "ESD"){
		fESD = dynamic_cast<AliESDEvent*>(event);

		if(!fESD){
			Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}
	}

	// -----------------------------------------------------
	//			 M C
	// -----------------------------------------------------



	fHistEventCounter->Fill(0.5); 	//	All events

	AliAnalysisUtils * utils = new AliAnalysisUtils();
	if (!utils)
	{
		cout<<"------- No AnalysisUtils Object Found --------"<<utils<<endl;
		return;
	}

	//Bool_t fSPDCvsTCutStatus = kFALSE;

	//fSPDCvsTCutStatus = utils->IsSPDClusterVsTrackletBG(fESD);

	// Cuts at event level
	UInt_t fSelectMask= fInputHandler->IsEventSelected();
	Bool_t isINT7selected = fSelectMask&AliVEvent::kINT7;
	if(!isINT7selected)
		return;
	fHistEventCounter->Fill(1.5);

	if (!fEventCuts.AcceptEvent(event)) {
		PostData(1, fListOfObjects);
		return;
	}

	fHistEventCounter->Fill(2.5);
	MakeAnalysis(0.8);
	cout<<"hello!!!"<<endl;

	// Post output data.
	PostData(1, fListOfObjects);

}
//________________________________________________________________________
void AliAnalysisTaskUePid::MakeAnalysis( Double_t etaCut ){


	// selection on leading particle
	Double_t pt_leading    = 0;
	Double_t p_leading    = 0;
	Double_t eta_leading    = 0;
	Double_t phi_leading    = 0;
	Double_t dedx_leading    = 0;

	Int_t    i_leading = 0;

	index_data = GetIndexPeriod(fdata_set);

	for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

		AliESDtrack* esdTrack = fESD->GetTrack(i);

		Double_t eta      = esdTrack->Eta();
		Double_t phi      = esdTrack->Phi();
		Double_t momentum = esdTrack->P();
		Double_t pt       = esdTrack->Pt();
		Float_t  dedx     = esdTrack->GetTPCsignal();
		Float_t  dedxUnc  = esdTrack->GetTPCsignal();

		if(TMath::Abs(eta) > etaCut)
			continue;

		//quality cuts, standard 2015 track cuts
		if(!fTrackFilter->IsSelected(esdTrack))
			continue;

		if(pt<0.5)
			continue;

		Short_t ncl = esdTrack->GetTPCsignalN();
		if(ncl<fNcl)
			continue;



		// recalibration of dE/dx
		if(eta < 0){
			dedx *= 50/EtaCalibrationNeg(index_data,eta);
		}
		else{
			dedx *= 50/EtaCalibrationPos(index_data,eta);
		}

		// test using pions at MIP,  0.4<p<0.6 GeV/c
		if(momentum>0.4 && momentum<0.6){
			if( dedxUnc < fDeDxMIPMax && dedxUnc > fDeDxMIPMin ){
				hMIPVsEtaBefore->Fill(eta,dedxUnc);
				pMIPVsEtaBefore->Fill(eta,dedxUnc);

				hMIPVsEta->Fill(eta,dedx);
				pMIPVsEta->Fill(eta,dedx);
			}
		}
		// For high pT analysis
		if(!PhiCut(pt, phi, esdTrack->Charge(), Magf, fcutLow, fcutHigh))
			continue;

		if(pt>pt_leading){
			pt_leading      = pt;
			p_leading       = momentum;
			eta_leading     = eta;
			phi_leading     = phi;
			dedx_leading    = dedx;
			i_leading = i;
		}

		hpT->Fill(pt);

	}// end loop over tracks

	if(pt_leading<2.0)
		return;

	Int_t index_leading = GetIndexEta(eta_leading);
	hPtL[0]->Fill(pt_leading);
	hEtaL[0]->Fill(eta_leading);
	hPhiL[0]->Fill(phi_leading);
	hDeDxL[0]->Fill(p_leading,dedx_leading);
	if(index_leading>0){

		hPtL[index_leading]->Fill(pt_leading);
		hEtaL[index_leading]->Fill(eta_leading);
		hPhiL[index_leading]->Fill(phi_leading);
		hDeDxL[index_leading]->Fill(p_leading,dedx_leading);

	}

	// Next step: DeDx vs p vs Number density (NS, AS, TS)
	Double_t mult_ns = 0;
	Double_t mult_nns = 0;
	Double_t mult_as = 0;
	Double_t mult_ts = 0;
	// DeDx vs p vs Sum pT density (NS, AS, TS)
	Double_t pt_ns = 0;
	Double_t pt_nns = 0;
	Double_t pt_as = 0;
	Double_t pt_ts = 0;

	for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

		// exclude the auto-correlation
		if(i==i_leading)
			continue;

		AliESDtrack* esdTrack = fESD->GetTrack(i);

		Double_t eta      = esdTrack->Eta();
		Double_t phi      = esdTrack->Phi();
		Double_t pt       = esdTrack->Pt();

		if(TMath::Abs(eta) > etaCut)
			continue;

		//quality cuts, standard 2015 track cuts
		if(!fTrackFilter->IsSelected(esdTrack))
			continue;

		if(pt<0.5)// only above 500 GeV/c
			continue;

		Double_t DPhi = DeltaPhi( phi, phi_leading );
		Double_t DEta = TMath::Abs( eta -  eta_leading );
		Double_t R = TMath::Sqrt(DPhi*DPhi+DEta*DEta);

		hDphi->Fill(DPhi);

		// new near side
		if(R<0.2){
			mult_nns++;
			pt_nns+=pt;
			hDphiNNS->Fill(DPhi);
		}

		// near side
		if(TMath::Abs(DPhi)<pi/3.0){
			mult_ns++;
			pt_ns+=pt;
			hDphiNS->Fill(DPhi);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			mult_as++;
			pt_as+=pt;
			hDphiAS->Fill(DPhi);
		}
		else{
			mult_ts++;
			pt_ts+=pt;
			hDphiTS->Fill(DPhi);
		}


	}// end loop over tracks

	// areas
	Double_t total_a = 2*TMath::Pi()*1.6;
	Double_t ns_a = 2*(TMath::Pi()/3.0)*1.6;
	Double_t as_a = 2*(TMath::Pi()/3.0)*1.6;
	Double_t ts_a = total_a-(ns_a+as_a);

	// Plotting results for number density

	pNNS[0]->Fill(p_leading,mult_nns/(pi*0.5*0.5));
	pNS[0]->Fill(p_leading,mult_ns/ns_a);
	pAS[0]->Fill(p_leading,mult_as/as_a);
	pTS[0]->Fill(p_leading,mult_ts/ts_a);

	hNNS[0]->Fill(p_leading,1.0*mult_nns);
	hNS[0]->Fill(p_leading,1.0*mult_ns);
	hAS[0]->Fill(p_leading,1.0*mult_as);
	hTS[0]->Fill(p_leading,1.0*mult_ts);

	h3DNNS[0]->Fill(p_leading,mult_nns,dedx_leading);
	h3DNS[0]->Fill(p_leading,mult_ns,dedx_leading);
	h3DAS[0]->Fill(p_leading,mult_as,dedx_leading);
	h3DTS[0]->Fill(p_leading,mult_ts,dedx_leading);

	// Plotting results for sum pT density
	pNNS2[0]->Fill(p_leading,pt_nns/(pi*0.5*0.5));
	pNS2[0]->Fill(p_leading,pt_ns/ns_a);
	pAS2[0]->Fill(p_leading,pt_as/as_a);
	pTS2[0]->Fill(p_leading,pt_ts/ts_a);

	hNNS2[0]->Fill(p_leading,1.0*pt_nns);
	hNS2[0]->Fill(p_leading,1.0*pt_ns);
	hAS2[0]->Fill(p_leading,1.0*pt_as);
	hTS2[0]->Fill(p_leading,1.0*pt_ts);

	h3DNNS2[0]->Fill(p_leading,pt_nns,dedx_leading);
	h3DNS2[0]->Fill(p_leading,pt_ns,dedx_leading);
	h3DAS2[0]->Fill(p_leading,pt_as,dedx_leading);
	h3DTS2[0]->Fill(p_leading,pt_ts,dedx_leading);

	// momentum vs pt
	hptVSp[0]->Fill(pt_leading,p_leading);

	if(index_leading>0){

		// Plotting results for number density
		pNNS[index_leading]->Fill(p_leading,mult_nns/(pi*0.5*0.5));
		pNS[index_leading]->Fill(p_leading,mult_ns/ns_a);
		pAS[index_leading]->Fill(p_leading,mult_as/as_a);
		pTS[index_leading]->Fill(p_leading,mult_ts/ts_a);

		hNNS[index_leading]->Fill(p_leading,1.0*mult_nns);
		hNS[index_leading]->Fill(p_leading,1.0*mult_ns);
		hAS[index_leading]->Fill(p_leading,1.0*mult_as);
		hTS[index_leading]->Fill(p_leading,1.0*mult_ts);

		h3DNNS[index_leading]->Fill(p_leading,mult_nns,dedx_leading);
		h3DNS[index_leading]->Fill(p_leading,mult_ns,dedx_leading);
		h3DAS[index_leading]->Fill(p_leading,mult_as,dedx_leading);
		h3DTS[index_leading]->Fill(p_leading,mult_ts,dedx_leading);

		// Plotting results for sum pT density
		pNNS2[index_leading]->Fill(p_leading,pt_nns/(pi*0.5*0.5));
		pNS2[index_leading]->Fill(p_leading,pt_ns/ns_a);
		pAS2[index_leading]->Fill(p_leading,pt_as/as_a);
		pTS2[index_leading]->Fill(p_leading,pt_ts/ts_a);

		hNNS2[index_leading]->Fill(p_leading,1.0*pt_nns);
		hNS2[index_leading]->Fill(p_leading,1.0*pt_ns);
		hAS2[index_leading]->Fill(p_leading,1.0*pt_as);
		hTS2[index_leading]->Fill(p_leading,1.0*pt_ts);

		h3DNNS2[index_leading]->Fill(p_leading,pt_nns,dedx_leading);
		h3DNS2[index_leading]->Fill(p_leading,pt_ns,dedx_leading);
		h3DAS2[index_leading]->Fill(p_leading,pt_as,dedx_leading);
		h3DTS2[index_leading]->Fill(p_leading,pt_ts,dedx_leading);

		// momentum vs pt
		hptVSp[index_leading]->Fill(pt_leading,p_leading);

	}


}
//____________________________________________________________
void AliAnalysisTaskUePid::SetTrackCuts(AliAnalysisFilter* fTrackFilter){

	AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
	// TPC
	esdTrackCuts->SetMinNCrossedRowsTPC(70);
	esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);

	esdTrackCuts->SetMaxChi2PerClusterTPC(4);
	esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
	esdTrackCuts->SetRequireTPCRefit(kTRUE);
	// ITS
	esdTrackCuts->SetRequireITSRefit(kTRUE);
	esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
			AliESDtrackCuts::kAny);
	// 7*(0.0015+0.0050/pt^1.1)
	esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");

	esdTrackCuts->SetMaxDCAToVertexZ(2);
	esdTrackCuts->SetDCAToVertex2D(kFALSE);
	esdTrackCuts->SetRequireSigmaToVertex(kFALSE);

	esdTrackCuts->SetMaxChi2PerClusterITS(36);

	/*
	   AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts();
	   esdTrackCuts->SetMaxFractionSharedTPCClusters(0.4);//
	   esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);//
	   esdTrackCuts->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);//
	   esdTrackCuts->SetMaxChi2PerClusterTPC(4);//
	   esdTrackCuts->SetAcceptKinkDaughters(kFALSE);//
	   esdTrackCuts->SetRequireTPCRefit(kTRUE);//
	   esdTrackCuts->SetRequireITSRefit(kTRUE);//
	   esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
	   AliESDtrackCuts::kAny);//
	   esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");//
	   esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);//
	   esdTrackCuts->SetMaxDCAToVertexZ(2);//
	   esdTrackCuts->SetDCAToVertex2D(kFALSE);//
	   esdTrackCuts->SetRequireSigmaToVertex(kFALSE);//
	   esdTrackCuts->SetMaxChi2PerClusterITS(36);//
	 */
	fTrackFilter->AddCuts(esdTrackCuts);

}
//______________________________________________________________
Bool_t AliAnalysisTaskUePid::PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t   mag, TF1* phiCutLow, TF1* phiCutHigh)
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
// Calibration factors provided by Omar
Double_t AliAnalysisTaskUePid::EtaCalibrationNeg( const Int_t Cent, const Double_t eta){


	for(Int_t i=0; i<8; ++i)
		fEtaCalibrationNeg->SetParameter(i,0);

	fEtaCalibrationNeg->SetParameter(0,aNeg[Cent]);
	fEtaCalibrationNeg->SetParameter(1,bNeg[Cent]);
	fEtaCalibrationNeg->SetParameter(2,cNeg[Cent]);
	fEtaCalibrationNeg->SetParameter(3,dNeg[Cent]);
	fEtaCalibrationNeg->SetParameter(4,eNeg[Cent]);
	fEtaCalibrationNeg->SetParameter(5,fNeg[Cent]);
	fEtaCalibrationNeg->SetParameter(6,gNeg[Cent]);
	fEtaCalibrationNeg->SetParameter(7,hNeg[Cent]);

	return fEtaCalibrationNeg->Eval(eta);


}
Double_t AliAnalysisTaskUePid::EtaCalibrationPos( const Int_t Cent, const Double_t eta){


	for(Int_t i=0; i<8; ++i)
		fEtaCalibrationPos->SetParameter(i,0);

	fEtaCalibrationPos->SetParameter(0,aPos[Cent]);
	fEtaCalibrationPos->SetParameter(1,bPos[Cent]);
	fEtaCalibrationPos->SetParameter(2,cPos[Cent]);
	fEtaCalibrationPos->SetParameter(3,dPos[Cent]);
	fEtaCalibrationPos->SetParameter(4,ePos[Cent]);
	fEtaCalibrationPos->SetParameter(5,fPos[Cent]);
	fEtaCalibrationPos->SetParameter(6,gPos[Cent]);
	fEtaCalibrationPos->SetParameter(7,hPos[Cent]);

	return fEtaCalibrationPos->Eval(eta);

}
Int_t AliAnalysisTaskUePid::GetIndexEta(Double_t etain){

	Int_t index = -1;
	if(TMath::Abs(etain)<0.2)
		index = 1;
	if(TMath::Abs(etain)>=0.2&&TMath::Abs(etain)<0.4)
		index = 2;
	if(TMath::Abs(etain)>=0.4&&TMath::Abs(etain)<0.6)
		index = 3;
	if(TMath::Abs(etain)>=0.6&&TMath::Abs(etain)<=0.8)
		index = 4;

	return index;

}
Double_t AliAnalysisTaskUePid::DeltaPhi(Double_t phia, Double_t phib,
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

Int_t AliAnalysisTaskUePid::GetIndexPeriod(TString lPeriod)
{
	Int_t index_period = -1;
	if(lPeriod=="l")
		index_period = 0;
	if(lPeriod=="k")
		index_period = 1;
	return index_period;
}
