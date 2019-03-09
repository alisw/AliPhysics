/*   This macro produces: Delta phi correclation in different multiplicity classes and leading particle dE/dx vs leading particle momentum vs Multiplicity (near side, away side and transverse side) 
   Aditya Nath Mishra ICN-UNAM
   Please report bugs to: amishra@cern.ch / aditya.mishra@correo.nucleares.unam.mx 
   First version: 07/03/2019

   Calibration factors of dE/dx were provided by Omar Vazquez (Lund University)

 */

#include "AliAnalysisTaskUeMultDep.h"

// ROOT includes
#include <TList.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2D.h>
#include <TProfile.h>
#include <THnSparse.h>
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

const Char_t * legEta[5] = {"|#eta|<0.8","|#eta|<0.2","0.2<|#eta|<0.4","0.4<|#eta|<0.6","0.6<|#eta|<0.8"};
const Char_t * legpTL[10] = {"0.15<#it{p}^{leading}<200", "1.5<#it{p}^{leading}<2.0", "4.0<#it{p}^{leading}<4.5", "7.0<#it{p}^{leading}<8.0", "10<#it{p}^{leading}<20", "20<#it{p}^{leading}<30", "30<#it{p}^{leading}<50", "50<#it{p}^{leading}<70", "70<#it{p}^{leading}<100", "100<#it{p}^{leading}<300"};

Float_t MagF                    = 1;
const Int_t nPeriod               = 2;
Int_t index_data0 = 0;
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



ClassImp(AliAnalysisTaskUeMultDep)

	//_____________________________________________________________________________
	AliAnalysisTaskUeMultDep::AliAnalysisTaskUeMultDep():
		AliAnalysisTaskSE(),
		fESD(0x0),
		fdata_set("k"),
		fEventCuts(0x0),
		fTrackFilter(0x0),
		fAnalysisType("ESD"),
		fcutLow(0x0),
		fcutHigh(0x0),
		fDeDxMIPMin(40),
		fDeDxMIPMax(90),
		fEtaCalibrationNeg(0x0),
		fEtaCalibrationPos(0x0),
		fNcl(70),
		fListOfObjects(0x0),
		fHistEventCounter(0x0),
		fEvents(0x0),
		hPhi(0x0),
		hMIPVsEta(0x0),
		pMIPVsEta(0x0),
		hpT(0x0),
		hDphi(0x0),
		hDphiNNS(0x0),
		hDphiNS(0x0),
		hDphiAS(0x0),
	        hDphiTS(0x0),
	        hRefMult08(0x0),
		hV0Mmult(0x0),
		hpTvsNch(0x0),
		hpTLvsNch(0x0),
		hpTvsRefMult08(0x0),
		hpTLvsRefMult08(0x0),
		hpTvsV0Mmult(0x0),
		hpTLvsV0Mmult(0x0),
		hRefMultvsV0Mmult(0x0),
		hpTvsV0MmultvsRefMult08(0x0),
		hpTLvsV0MmultvsRefMult08(0x0),
                hpTLvsRefMult08vsDphi(0x0),
		hpTLvsV0MmultvsDphi(0x0),
                hptLvsv0MvsRefMultvsDphivsDeta(0x0),
                fMultSelection(0x0),
                ftrackmult08(-999),   
                fv0mpercentile(-999)

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

		hptVSp[i_eta] = 0;
	}


}

//______________________________________________________________________________
AliAnalysisTaskUeMultDep::AliAnalysisTaskUeMultDep(const char *name):
	AliAnalysisTaskSE(name),
	fESD(0x0),
	fdata_set("k"),
	fEventCuts(0x0),
	fTrackFilter(0x0),
	fAnalysisType("ESD"),
	fcutLow(0x0),
	fcutHigh(0x0),
	fDeDxMIPMin(40),
	fDeDxMIPMax(90),
	fEtaCalibrationNeg(0x0),
	fEtaCalibrationPos(0x0),
	fNcl(70),
	fListOfObjects(0x0),
	fHistEventCounter(0x0), 
	fEvents(0x0),
	hPhi(0x0),
	hMIPVsEta(0x0),
	pMIPVsEta(0x0),
	hpT(0x0),
	hDphi(0x0),
	hDphiNNS(0x0),
	hDphiNS(0x0),
	hDphiAS(0x0),
	hDphiTS(0x0),
	hRefMult08(0x0),
	hV0Mmult(0x0),
	hpTvsNch(0x0),
	hpTLvsNch(0x0),
	hpTvsRefMult08(0x0),
	hpTLvsRefMult08(0x0),
	hpTvsV0Mmult(0x0),
	hpTLvsV0Mmult(0x0),
	hRefMultvsV0Mmult(0x0),
	hpTvsV0MmultvsRefMult08(0x0),
	hpTLvsV0MmultvsRefMult08(0x0),
	hpTLvsRefMult08vsDphi(0x0),
	hpTLvsV0MmultvsDphi(0x0),
	hptLvsv0MvsRefMultvsDphivsDeta(0x0),
	fMultSelection(0x0),
	ftrackmult08(-999),   
	fv0mpercentile(-999)
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

		hptVSp[i_eta] = 0;
	}

	DefineOutput(1, TList::Class());
}

//________________________________________________________________________

void AliAnalysisTaskUeMultDep::Exit(const char *msg) {

	Printf("%s", msg);
	return;
}


//_____________________________________________________________________________
AliAnalysisTaskUeMultDep::~AliAnalysisTaskUeMultDep()
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
void AliAnalysisTaskUeMultDep::UserCreateOutputObjects()
{

	fcutLow = new TF1("StandardPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 50);
	fcutHigh = new TF1("StandardPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 50);

	fEtaCalibrationNeg = new TF1("fDeDxVsEtaNeg", "pol7", -1.0, 0.0);
	fEtaCalibrationPos    = new TF1("fDeDxVsEtaPos", "pol7", 0.0, 1.0);
	
	const Int_t nNchBins = 200;
	Double_t NchBins[nNchBins+1]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
				21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,
				39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,
				57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,
				75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,
				93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,
				109,110,111,112,113,114,115,116,117,118,119,120,121,122,
				123,124,125,126,127,128,129,130,131,132,133,134,135,136,
				137,138,139,140,141,142,143,144,145,146,147,148,149,150,
				151,152,153,154,155,156,157,158,159,160,161,162,163,164,
				165,166,167,168,169,170,171,172,173,174,175,176,177,178,
				179,180,181,182,183,184,185,186,187,188,189,190,191,192,
				193,194,195,196,197,198,199,200};

	const Int_t nPtBins      = 61;
	Double_t xBins[nPtBins+1] = {
                             0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,
			     1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.5,5,5.5,6,
			     6.5,7,8,9,10,11,12,13,14,15,16,18,20,22.0,24.0,26.0,28.0,32.0,36.0,
			     42.0,50.0,60.0,80.0,100.0,130.0,160.0,200.0};

	const Int_t nV0Mbins = 10;
	Double_t V0Mbins[nV0Mbins+1]={0.00, 1.00, 5.00, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100};

	const Int_t nDeDxBins = 50;
	Double_t DeDxBins[nDeDxBins+1]=
	{
		40,41,42,43,44,45,46,47,48,49,
		50,51,52,53,54,55,56,57,58,59,
		60,61,62,63,64,65,66,67,68,69,
		70,71,72,73,74,75,76,77,78,79,
		80,81,82,83,84,85,86,87,88,89,90
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

	hRefMult08 = 0;
	hRefMult08 = new TH1D("hRefMult08","Multiplicity (-0.8 < #eta < 0.8);N_{ch};count",nNchBins,NchBins);
	fListOfObjects->Add(hRefMult08);
	
	hV0Mmult = 0;
	hV0Mmult = new TH1D("hV0Mmult","V0M ;V0M percentile;count",200,0,200);
	fListOfObjects->Add(hV0Mmult);
	
	hpTvsNch = 0;
	hpTvsNch = new TH2D("hpTvsNch","p_{T} vs N_{ch}; #it{p}_{T} (GeV/#it{c}); N_{ch} (|#eta| < 0.8 & p_{T} > 0.5 GeV/c)",nPtBins,xBins,200,0,200);
	fListOfObjects->Add(hpTvsNch);
	
	hpTLvsNch = 0;
	hpTLvsNch = new TH2D("hpTLvsNch","p_{T}^{leading} vs N_{ch}; #it{p}_{T}^{leading} (GeV/#it{c});N_{ch} (|#eta| < 0.8 & p_{T} > 0.5 GeV/c)",nPtBins,xBins,200,0,200);
	fListOfObjects->Add(hpTLvsNch);
		
	hpTvsRefMult08 = 0;
	hpTvsRefMult08 = new TH2D("hpTvsRefMult08","p_{T} vs RefMult08; #it{p}_{T} (GeV/#it{c}); RefMult08 (|#eta| < 0.8 & p_{T} > 0.5 GeV/c)",nPtBins,xBins,200,0,200);
	fListOfObjects->Add(hpTvsRefMult08);
	
	hpTLvsRefMult08 = 0;
	hpTLvsRefMult08 = new TH2D("hpTLvsRefMult08","p_{T}^{leading} vs RefMult08; #it{p}_{T}^{leading} (GeV/#it{c}); RefMult08 (|#eta| < 0.8 & p_{T} > 0.5 GeV/c)",nPtBins,xBins,200,0,200);
	fListOfObjects->Add(hpTLvsRefMult08);
	
	hpTvsV0Mmult = 0;
	hpTvsV0Mmult = new TH2D("hpTvsV0Mmult","p_{T} vs V0Mmult; #it{p}_{T} (GeV/#it{c}); V0Mmult (|#eta| < 0.8 & p_{T} > 0.5 GeV/c)",nPtBins,xBins,200,0,200);
	fListOfObjects->Add(hpTvsV0Mmult);
	
	hpTLvsV0Mmult = 0;
	hpTLvsV0Mmult = new TH2D("hpTLvsV0Mmult","p_{T}^{leading} vs V0Mmult; #it{p}_{T}^{leading} (GeV/#it{c}); V0Mmult (|#eta| < 0.8 & p_{T} > 0.5 GeV/c)",nPtBins,xBins,200,0,200);
	fListOfObjects->Add(hpTLvsV0Mmult);
	
	hRefMultvsV0Mmult = 0;
	hRefMultvsV0Mmult = new TH2D("hRefMultvsV0Mmult","N_{ch} vs V0M percentile;N_{ch}; v0M percentile",nNchBins,NchBins,200,0,200);
	fListOfObjects->Add(hRefMultvsV0Mmult);

	hpTvsV0MmultvsRefMult08 = 0;
	hpTvsV0MmultvsRefMult08 = new TH3D("hpTvsV0MmultvsRefMult08","p_{T} vs v0M percentile vs RefMult08;#it{p}_{T} (GeV/c);V0M percentile;RefMult08",nPtBins,xBins,nV0Mbins,V0Mbins,nNchBins,NchBins);
	fListOfObjects->Add(hpTvsV0MmultvsRefMult08);

	hpTLvsV0MmultvsRefMult08 = 0;
	hpTLvsV0MmultvsRefMult08 = new TH3D("hpTLvsV0MmultvsRefMult08","p_{T}^{Leading} vs v0M percentile vs RefMult08;#it{p}_{T}^{Leading} (GeV/c);V0M percentile;RefMult08",nPtBins,xBins,nV0Mbins,V0Mbins,nNchBins,NchBins);
	fListOfObjects->Add(hpTLvsV0MmultvsRefMult08);
	
	hpTLvsRefMult08vsDphi = 0;
	hpTLvsRefMult08vsDphi = new TH3D("hpTLvsRefMult08vsDphi","p_{T}^{Leading} vs RefMult08 vs Dphi;#it{p}_{T}^{Leading} (GeV/c);RefMult08; #delta#phi (rad)",200,0,200,200,0,200,64,-(TMath::Pi())/2.0,3.0*(TMath::Pi())/2.0);
	fListOfObjects->Add(hpTLvsRefMult08vsDphi);

	hpTLvsV0MmultvsDphi = 0;
	hpTLvsV0MmultvsDphi = new TH3D("hpTLvsV0MmultvsDphi","p_{T}^{Leading} vs V0Mmult vs Dphi;#it{p}_{T}^{Leading} (GeV/c);V0Mmult; #delta#phi (rad)",200,0,200,200,0,200,64,-(TMath::Pi())/2.0,3.0*(TMath::Pi())/2.0);
	fListOfObjects->Add(hpTLvsV0MmultvsDphi);

	for(Int_t i_eta=0; i_eta<5; ++i_eta){
		hPtL[i_eta] = 0;
		hPtL[i_eta] = new TH1D(Form("hPtL_etaInt%d",i_eta),Form("%s; #it{p}_{T}^{leading} (GeV/#it{c});counts",legEta[i_eta]),nPtBins,xBins);
		fListOfObjects->Add( hPtL[i_eta] );

		hDeDxL[i_eta] = 0;
		hDeDxL[i_eta] = new TH2D(Form("hDeDxL_etaInt%d",i_eta),Form("%s; TPC-d#it{E}/d#it{x}^{leading};counts",legEta[i_eta]),nPtBins,xBins,nDeDxBins,DeDxBins);
		fListOfObjects->Add( hDeDxL[i_eta] );

		hEtaL[i_eta] = 0;
		hEtaL[i_eta] = new TH1D(Form("hEtaL_etaInt%d",i_eta),Form("%s; #eta^{leading};counts",legEta[i_eta]),20,-1,1);
		fListOfObjects->Add(hEtaL[i_eta]);

		hPhiL[i_eta] = 0;
		hPhiL[i_eta] = new TH1D(Form("hPhiL_etaInt%d",i_eta),Form("%s; #phi^{leading} (rad);counts",legEta[i_eta]),64,0,2.0*TMath::Pi());
		fListOfObjects->Add(hPhiL[i_eta]);

		pNNS[i_eta] = 0;
		pNNS[i_eta] = new TProfile(Form("pNNS%d",i_eta),Form("R<0.5 (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legEta[i_eta]),nPtBins,xBins,0,20);
		fListOfObjects->Add(pNNS[i_eta]);

		pNS[i_eta] = 0;
		pNS[i_eta] = new TProfile(Form("pNS%d",i_eta),Form("Near side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legEta[i_eta]),nPtBins,xBins,0,20);
		fListOfObjects->Add(pNS[i_eta]);

		pAS[i_eta] = 0;
		pAS[i_eta] = new TProfile(Form("pAS%d",i_eta),Form("Away side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legEta[i_eta]),nPtBins,xBins,0,20);
		fListOfObjects->Add(pAS[i_eta]);

		pTS[i_eta] = 0;
		pTS[i_eta] = new TProfile(Form("pTS%d",i_eta),Form("Transverse side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legEta[i_eta]),nPtBins,xBins,0,20);
		fListOfObjects->Add(pTS[i_eta]);

		hNNS[i_eta] = 0;
		hNNS[i_eta] = new TH2D(Form("hNNS%d",i_eta),Form("R<0.5 (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legEta[i_eta]),nPtBins,xBins,nNchBins,NchBins);
		fListOfObjects->Add(hNNS[i_eta]);

		hNS[i_eta] = 0;
		hNS[i_eta] = new TH2D(Form("hNS%d",i_eta),Form("Near side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legEta[i_eta]),nPtBins,xBins,nNchBins,NchBins);
		fListOfObjects->Add(hNS[i_eta]);

		hAS[i_eta] = 0;
		hAS[i_eta] = new TH2D(Form("hAS%d",i_eta),Form("Away side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legEta[i_eta]),nPtBins,xBins,nNchBins,NchBins);
		fListOfObjects->Add(hAS[i_eta]);

		hTS[i_eta] = 0;
		hTS[i_eta] = new TH2D(Form("hTS%d",i_eta),Form("Transverse side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c})",legEta[i_eta]),nPtBins,xBins,nNchBins,NchBins);
		fListOfObjects->Add(hTS[i_eta]);

		h3DNNS[i_eta] = 0;
		h3DNNS[i_eta] = new TH3D(Form("h3DNNS%d",i_eta),Form("#it{R}<0.5 (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c}); TPC-d#it{E}/d#it{x}^{leading}",legEta[i_eta]),nPtBins,xBins,nNchBins,NchBins,nDeDxBins,DeDxBins);
		fListOfObjects->Add(h3DNNS[i_eta]);

		h3DNS[i_eta] = 0;
		h3DNS[i_eta] = new TH3D(Form("h3DNS%d",i_eta),Form("Near side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c}); TPC-d#it{E}/d#it{x}^{leading}",legEta[i_eta]),nPtBins,xBins,nNchBins,NchBins,nDeDxBins,DeDxBins);
		fListOfObjects->Add(h3DNS[i_eta]);

		h3DAS[i_eta] = 0;
		h3DAS[i_eta] = new TH3D(Form("h3DAS%d",i_eta),Form("Away side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c}); TPC-d#it{E}/d#it{x}^{leading}",legEta[i_eta]),nPtBins,xBins,nNchBins,NchBins,nDeDxBins,DeDxBins);
		fListOfObjects->Add(h3DAS[i_eta]);

		h3DTS[i_eta] = 0;
		h3DTS[i_eta] = new TH3D(Form("h3DTS%d",i_eta),Form("Transverse side (%s); #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.5 GeV/#it{c}); TPC-d#it{E}/d#it{x}^{leading}",legEta[i_eta]),nPtBins,xBins,nNchBins,NchBins,nDeDxBins,DeDxBins);
		fListOfObjects->Add(h3DTS[i_eta]);


		hptVSp[i_eta] = 0;
		hptVSp[i_eta] = new TH2D(Form("hptVSp%d",i_eta),Form("%s; #it{p}_{T}^{leading} (GeV/#it{c}); #it{p}^{leading} (GeV/#it{c})",legEta[i_eta]),nPtBins,xBins,nPtBins,xBins);
		fListOfObjects->Add(hptVSp[i_eta]);


	}

		// defining THnSparseD....
	const Int_t dims = 5;
	
	Int_t bins[dims] = {100, 200, 200, 64, 20};
	Double_t xmin[dims] = {0., 0, 0, -(TMath::Pi())/2.0, -1.0};
	Double_t xmax[dims] = {100., 200, 200, 3.0*(TMath::Pi())/2.0, 1.0};
	
	hptLvsv0MvsRefMultvsDphivsDeta = new THnSparseD("hptLvsv0MvsRefMultvsDphivsDeta","", dims, bins, xmin, xmax);
	fListOfObjects->Add(hptLvsv0MvsRefMultvsDphivsDeta);

	fEventCuts.AddQAplotsToList(fListOfObjects);

	PostData(1, fListOfObjects);

}

//______________________________________________________________________________
void AliAnalysisTaskUeMultDep::UserExec(Option_t *)
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

	// -------------------------------------- multiplcity estimators section ------------------------------------------ //

	ftrackmult08 = -999;
	fv0mpercentile = -999;

	//ftrackmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);     //reference
	ftrackmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 0.8);     //tracklets
	//ftrackmult08 = AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD); //Combined estimator

	hRefMult08->Fill(ftrackmult08);

	fMultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection"); // Esto es para 13 TeV
	if (!fMultSelection)
		cout<<"------- No AliMultSelection Object Found --------"<<fMultSelection<<endl;
	else
	fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");
	hV0Mmult->Fill(fv0mpercentile);

	cout<<"------- V0M mult ==  "<<fv0mpercentile<<"--------"<<endl;
	
	hRefMultvsV0Mmult->Fill(ftrackmult08,fv0mpercentile);

	// ------------------------------------------ end of mult estimators -------------------------------------------------//
	
	MakeAnalysis(0.8);
	cout<<"hello!!!"<<endl;

	// Post output data.
	PostData(1, fListOfObjects);

}
//________________________________________________________________________
void AliAnalysisTaskUeMultDep::MakeAnalysis( Double_t etaCut ){


	// selection on leading particle
	Double_t pt_leading    = 0;
	Double_t p_leading    = 0;
	Double_t eta_leading    = 0;
	Double_t phi_leading    = 0;
	Double_t dedx_leading    = 0;

	Int_t    i_leading = 0;
	Double_t vals[5];

	index_data0 = GetIndexPeriod(fdata_set);

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
			dedx *= 50/EtaCalibrationNeg(index_data0,eta);
		}
		else{
			dedx *= 50/EtaCalibrationPos(index_data0,eta);
		}

		// test using pions at MIP
		if( dedxUnc < fDeDxMIPMax && dedxUnc > fDeDxMIPMin ){
			hMIPVsEta->Fill(eta,dedx);
			pMIPVsEta->Fill(eta,dedx);
		}
		// For high pT analysis
		if(!PhiCut(pt, phi, esdTrack->Charge(), MagF, fcutLow, fcutHigh))
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
	Double_t mult    = 0;

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
		mult++;
		
		Double_t DPhi = DeltaPhi( phi, phi_leading );
		Double_t DEta = TMath::Abs( eta -  eta_leading );
		Double_t R = TMath::Sqrt(DPhi*DPhi+DEta*DEta);
		Double_t DeltaEta = eta -  eta_leading;

		hDphi->Fill(DPhi);
		hpTvsNch->Fill(pt,mult);
		hpTvsRefMult08->Fill(pt,ftrackmult08);
		hpTvsV0Mmult->Fill(pt,fv0mpercentile);
		hpTvsV0MmultvsRefMult08->Fill(pt,fv0mpercentile,ftrackmult08);
		hpTLvsRefMult08vsDphi->Fill(pt_leading,ftrackmult08,DPhi);
		hpTLvsV0MmultvsDphi->Fill(pt_leading,fv0mpercentile,DPhi);

		vals[0]=pt_leading;
		vals[1]=fv0mpercentile;
		vals[2]=ftrackmult08;
		vals[3]=DPhi;
		vals[4]=DeltaEta;

		hptLvsv0MvsRefMultvsDphivsDeta->Fill(vals);

		// new near side
		if(R<0.5){
			mult_nns++;
			hDphiNNS->Fill(DPhi);
		}

		// near side
		if(TMath::Abs(DPhi)<pi/3.0){
			mult_ns++;
			hDphiNS->Fill(DPhi);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			mult_as++;
			hDphiAS->Fill(DPhi);
		}
		else{
			mult_ts++;
			hDphiTS->Fill(DPhi);
		}


	}// end loop over tracks

	hpTLvsNch->Fill(pt_leading,mult);
	hpTLvsRefMult08->Fill(pt_leading,ftrackmult08);
	hpTLvsV0Mmult->Fill(pt_leading,fv0mpercentile);
	hpTLvsV0MmultvsRefMult08->Fill(pt_leading,fv0mpercentile,ftrackmult08);
       
		
	// areas
	Double_t total_a = 2*TMath::Pi()*1.6;
	Double_t ns_a = 2*(TMath::Pi()/3.0)*1.6;
	Double_t as_a = 2*(TMath::Pi()/3.0)*1.6;
	Double_t ts_a = total_a-(ns_a+as_a);

	pNNS[0]->Fill(pt_leading,mult_nns/(pi*0.5*0.5));
	pNS[0]->Fill(pt_leading,mult_ns/ns_a);
	pAS[0]->Fill(pt_leading,mult_as/as_a);
	pTS[0]->Fill(pt_leading,mult_ts/ts_a);

	hNNS[0]->Fill(pt_leading,1.0*mult_nns);
	hNS[0]->Fill(pt_leading,1.0*mult_ns);
	hAS[0]->Fill(pt_leading,1.0*mult_as);
	hTS[0]->Fill(pt_leading,1.0*mult_ts);

	h3DNNS[0]->Fill(pt_leading,mult_nns,dedx_leading);
	h3DNS[0]->Fill(pt_leading,mult_ns,dedx_leading);
	h3DAS[0]->Fill(pt_leading,mult_as,dedx_leading);
	h3DTS[0]->Fill(pt_leading,mult_ts,dedx_leading);

	hptVSp[0]->Fill(pt_leading,p_leading);

	if(index_leading>0){
		pNNS[index_leading]->Fill(pt_leading,mult_nns/(pi*0.5*0.5));
		pNS[index_leading]->Fill(pt_leading,mult_ns/ns_a);
		pAS[index_leading]->Fill(pt_leading,mult_as/as_a);
		pTS[index_leading]->Fill(pt_leading,mult_ts/ts_a);

		hNNS[index_leading]->Fill(pt_leading,1.0*mult_nns);
		hNS[index_leading]->Fill(pt_leading,1.0*mult_ns);
		hAS[index_leading]->Fill(pt_leading,1.0*mult_as);
		hTS[index_leading]->Fill(pt_leading,1.0*mult_ts);

		h3DNNS[index_leading]->Fill(pt_leading,mult_nns,dedx_leading);
		h3DNS[index_leading]->Fill(pt_leading,mult_ns,dedx_leading);
		h3DAS[index_leading]->Fill(pt_leading,mult_as,dedx_leading);
		h3DTS[index_leading]->Fill(pt_leading,mult_ts,dedx_leading);

		hptVSp[index_leading]->Fill(pt_leading,p_leading);
	}


}
//____________________________________________________________
void AliAnalysisTaskUeMultDep::SetTrackCuts(AliAnalysisFilter* fTrackFilter){

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
Bool_t AliAnalysisTaskUeMultDep::PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t   mag, TF1* phiCutLow, TF1* phiCutHigh)
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
Double_t AliAnalysisTaskUeMultDep::EtaCalibrationNeg( const Int_t Cent, const Double_t eta){


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
Double_t AliAnalysisTaskUeMultDep::EtaCalibrationPos( const Int_t Cent, const Double_t eta){


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
Int_t AliAnalysisTaskUeMultDep::GetIndexEta(Double_t etain){

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
Double_t AliAnalysisTaskUeMultDep::DeltaPhi(Double_t phia, Double_t phib,
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

Int_t AliAnalysisTaskUeMultDep::GetIndexPeriod(TString lPeriod)
{
	Int_t index_period = -1;
	if(lPeriod=="l")
		index_period = 0;
	if(lPeriod=="k")
		index_period = 1;
	return index_period;
}
