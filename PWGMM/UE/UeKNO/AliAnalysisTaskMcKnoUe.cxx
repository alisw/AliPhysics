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
 *         Feng Fan (Feng.Fan@cern.ch)				                                  *
 *         Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)                *
 **************************************************************************/

/* AliAnaysisTaskMcKnoUe source code
 * The analysis task produce all the histos needed for MC closure test studies
 * Results include both KNO properties and traditional UE analysis
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

TF1* f_Eff;// efficiency for charged particles (2015 AA cuts) 
const Char_t * nameReg[3]={"NS","AS","TS"};
const Int_t nchNbins = 100;
Double_t nchbins[nchNbins+1]={-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5,29.5,30.5,31.5,32.5,33.5,34.5,35.5,36.5,37.5,38.5,39.5,40.5,41.5,42.5,43.5,44.5,45.5,46.5,47.5,48.5,49.5,50.5,51.5,52.5,53.5,54.5,55.5,56.5,57.5,58.5,59.5,60.5,61.5,62.5,63.5,64.5,65.5,66.5,67.5,68.5,69.5,70.5,71.5,72.5,73.5,74.5,75.5,76.5,77.5,78.5,79.5,80.5,81.5,82.5,83.5,84.5,85.5,86.5,87.5,88.5,89.5,90.5,91.5,92.5,93.5,94.5,95.5,96.5,97.5,98.5,99.5};

const Int_t ptNbins = 36;
Double_t ptbins1[ptNbins+1] = {
	0.0,  0.1,  0.15,  0.2,  0.25,  0.3,  0.35,  0.4,  0.45,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0,  1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5, 5.0, 6.0,    7.0,  8.0,  9.0,  10.0,  12.0,  14.0,  16.0,  18.0,  20.0,  25.0,  30.0,  40.0,  50.0
};

const Int_t ptNbinsL = 24;
Double_t ptbins1L[ptNbinsL+1] = {
        0.15, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50,
        5.00, 6.00, 7.00, 8.00, 9.00, 10.0, 12.0, 14.0, 16.0, 18.0,
        20.0, 25.0, 30.0, 40.0, 50.0
};


const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;
class AliAnalysisTaskMcKnoUe;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskMcKnoUe) // classimp: necessary for root

AliAnalysisTaskMcKnoUe::AliAnalysisTaskMcKnoUe() : AliAnalysisTaskSE(),
	fESD(0), fIsPythia(kTRUE), fEventCuts(0x0), fMCStack(0), fMC(0), fUseMC(kFALSE), fIsMCclosure(kFALSE), fLeadingTrackFilter(0x0), fTrackFilter(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5), fLeadPtCutMin(5.0), fLeadPtCutMax(40.0), fGenLeadPhi(0), fGenLeadPt(0), fGenLeadIn(0), fRecLeadPhi(0), fRecLeadPt(0), fRecLeadIn(0), hNchTSGen(0), hNchTSGenTest(0), hNchTSRec(0), hNchTSRecTest(0), hNchTSData(0), hNchResponse(0), hPtInPrim(0), hPtOut(0), hPtOutPrim(0), hPtOutSec(0), hCounter(0), hPtLeadingTrue(0), hPtLeadingMeasured(0), hPtLeadingData(0), hPtLeadingRecPS(0), hPtLeadingRecPSV(0), hPtLeadingRecGood(0), hPtLeadingGenPS(0), hPtLeadingGenPSV(0), hPtLeadingGenGood(0), hPtLeadingRecAll(0), hPtLeadingGenAll(0)   

{
	for(Int_t i=0;i<3;++i){ 
		hPhiGen[i]= 0;
		hPhiRec[i]= 0;
		hNumDenMC[i]=0;
		hSumPtMC[i]=0;
		hNumDenMCMatch[i]=0;
		hSumPtMCMatch[i]=0;
		// Data driven
		hNumDenMCDd[i]=0;
		hSumPtMCDd[i]=0;
		hNumDenMCMatchDd[i]=0;
		hSumPtMCMatchDd[i]=0;

		hPtVsPtLeadingMeasured[i]=0;
		hPtVsPtLeadingData[i]=0;// only for data
		hPtVsPtLeadingTrue[i]=0;

		pNumDenMeasured[i]=0;
		pSumPtMeasured[i]=0;
		pNumDenData[i]=0;// only for data
		pSumPtData[i]=0;// only for data


		pNumDenMeasuredAll[i]=0;
		pSumPtMeasuredAll[i]=0;
		pNumDenMeasuredPS[i]=0;
		pSumPtMeasuredPS[i]=0;
		pNumDenMeasuredPSV[i]=0;
		pSumPtMeasuredPSV[i]=0;
		pNumDenMeasuredGood[i]=0;
		pSumPtMeasuredGood[i]=0;


		pNumDenTrueAll[i]=0;
		pSumPtTrueAll[i]=0;
		pNumDenTruePS[i]=0;
		pSumPtTruePS[i]=0;
		pNumDenTruePSV[i]=0;
		pSumPtTruePSV[i]=0;
		pNumDenTrueGood[i]=0;
		pSumPtTrueGood[i]=0;

		pNumDenTrue[i]=0;
		pSumPtTrue[i]=0;
	}
	for(Int_t i=0;i<3;++i){
		hPtVsUEGenTest[i]=0;
		hPtVsUERecTest[i]=0;
		hPtVsUEData[i]=0;
	}


	// default constructor, don't allocate memory here!  this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskMcKnoUe::AliAnalysisTaskMcKnoUe(const char* name) : AliAnalysisTaskSE(name),
	fESD(0), fIsPythia(kTRUE), fEventCuts(0x0), fMCStack(0), fMC(0), fUseMC(kFALSE), fIsMCclosure(kFALSE), fLeadingTrackFilter(0x0), fTrackFilter(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5), fLeadPtCutMin(5.0), fLeadPtCutMax(40.0), fGenLeadPhi(0), fGenLeadPt(0), fGenLeadIn(0), fRecLeadPhi(0), fRecLeadPt(0), fRecLeadIn(0), hNchTSGen(0), hNchTSGenTest(0), hNchTSRec(0), hNchTSRecTest(0), hNchTSData(0), hNchResponse(0), hPtInPrim(0), hPtOut(0), hPtOutPrim(0), hPtOutSec(0), hCounter(0), hPtLeadingTrue(0), hPtLeadingMeasured(0), hPtLeadingData(0), hPtLeadingRecPS(0), hPtLeadingRecPSV(0), hPtLeadingRecGood(0), hPtLeadingGenPS(0), hPtLeadingGenPSV(0), hPtLeadingGenGood(0), hPtLeadingRecAll(0), hPtLeadingGenAll(0) 
{
	for(Int_t i=0;i<3;++i){

		hPhiGen[i] = 0;
		hPhiRec[i] = 0;
		hNumDenMC[i]=0;
		hSumPtMC[i]=0;
		hNumDenMCMatch[i]=0;
		hSumPtMCMatch[i]=0;
		// Data driven
		hNumDenMCDd[i]=0;
		hSumPtMCDd[i]=0;
		hNumDenMCMatchDd[i]=0;
		hSumPtMCMatchDd[i]=0;

		hPtVsPtLeadingMeasured[i]=0;
                hPtVsPtLeadingData[i]=0;// only for data
		hPtVsPtLeadingTrue[i]=0;		
		pNumDenMeasured[i]=0;
		pSumPtMeasured[i]=0;
		pNumDenData[i]=0;// only for data
		pSumPtData[i]=0;// only for data


		pNumDenMeasuredAll[i]=0;
		pSumPtMeasuredAll[i]=0;
		pNumDenMeasuredPS[i]=0;
		pSumPtMeasuredPS[i]=0;
		pNumDenMeasuredPSV[i]=0;
		pSumPtMeasuredPSV[i]=0;
		pNumDenMeasuredGood[i]=0;
		pSumPtMeasuredGood[i]=0;

		pNumDenTrueAll[i]=0;
		pSumPtTrueAll[i]=0;
		pNumDenTruePS[i]=0;
		pSumPtTruePS[i]=0;
		pNumDenTruePSV[i]=0;
		pSumPtTruePSV[i]=0;
		pNumDenTrueGood[i]=0;
		pSumPtTrueGood[i]=0;

		pNumDenTrue[i]=0;
		pSumPtTrue[i]=0;

	}
	for(Int_t i=0;i<3;++i){
		hPtVsUEGenTest[i]=0;
		hPtVsUERecTest[i]=0;
		hPtVsUEData[i]=0;
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

	// parametrization of efficiency
	f_Eff = 0;
	if(fIsPythia){
		f_Eff = new TF1("ch_Eff",
				"(x>=0.15&&x<[0])*([1]+x*[2])+(x>=[0]&&x<[3])*([4]+[5]*x*x+[6]*x*x*x+[7]*x)+(x>=[3])*([8])", 0.0, 1e2);
		f_Eff->SetParameters(9.00000e-01,9.30176e-01,-4.29864e-01,4.90000e+00,3.89778e-01,-5.81233e-02,5.41373e-03,2.20377e-01,7.10559e-01);
	}
	else{// epos lhc
		f_Eff = new TF1("ch_Eff",
				"(x>=0.15&&x<[0])*([1]+x*[2])+(x>=[0]&&x<[3])*([4]*x^[5]+exp([6]-[7]*x))+(x>=[3])*([8])", 0.0, 1e2);
		f_Eff->SetParameters(8.91265e-01,9.35736e-01,-4.45849e-01,9.50001e+00,-1.10333e-01,-1.02772e+00,-4.26261e-01,-7.56957e-03,7.13602e-01);
	}

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


	if(fUseMC)
	{
		hNchTSGen = new TH1D("hNchTSGen","",100,-0.5,99.5);
		fOutputList->Add(hNchTSGen);

		hNchTSGenTest = new TH1D("hNchTSGenTest","",100,-0.5,99.5); 
		fOutputList->Add(hNchTSGenTest);

		hNchTSRecTest = new TH1D("hNchTSRecTest","",100,-0.5,99.5); 
		fOutputList->Add(hNchTSRecTest);

		for(Int_t i=0;i<3;++i){
			hPhiGen[i]= new TH1D(Form("hPhiGen_%s",nameReg[i]),"",64,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
			fOutputList->Add(hPhiGen[i]);
		}

		hNchResponse = new TH2D("hNchResponse","Detector response; rec mult; gen mult",100,-0.5,99.5,100,-0.5,99.5);
		fOutputList->Add(hNchResponse);

		// UE analysis
		hPtInPrim = new TH1D("hPtInPrim","pT prim true; pT; Nch",ptNbins,ptbins1);
		fOutputList->Add(hPtInPrim);

		hPtOut = new TH1D("hPtOut","pT all rec; pT; Nch",ptNbins,ptbins1);
		fOutputList->Add(hPtOut);

		hPtOutPrim = new TH1D("hPtOutPrim","pT prim rec; pT; Nch",ptNbins,ptbins1);
		fOutputList->Add(hPtOutPrim);

		hPtOutSec = new TH1D("hPtOutSec","pT sec rec; pT; Nch",ptNbins,ptbins1);
		fOutputList->Add(hPtOutSec);

		for(Int_t i=0;i<3;++i){
			hNumDenMC[i]= new TH2D(Form("hNumDenMC_%s",nameReg[i]),"",ptNbinsL,ptbins1L,nchNbins,nchbins);
			fOutputList->Add(hNumDenMC[i]);
			hSumPtMC[i]= new TH2D(Form("hSumPtMC_%s",nameReg[i]),"",ptNbinsL,ptbins1L,ptNbins,ptbins1);
			fOutputList->Add(hSumPtMC[i]);
			hNumDenMCMatch[i]= new TH2D(Form("hNumDenMCMatch_%s",nameReg[i]),"",ptNbinsL,ptbins1L,nchNbins,nchbins);
			fOutputList->Add(hNumDenMCMatch[i]);
			hSumPtMCMatch[i]= new TH2D(Form("hSumPtMCMatch_%s",nameReg[i]),"",ptNbinsL,ptbins1L,ptNbins,ptbins1);
			fOutputList->Add(hSumPtMCMatch[i]);


			hPtVsPtLeadingTrue[i] = new TH2D(Form("hPtVsPtLeadingTrue_%s",nameReg[i]),"",ptNbinsL,ptbins1L,ptNbins,ptbins1);
			fOutputList->Add(hPtVsPtLeadingTrue[i]);

			hPtVsPtLeadingMeasured[i] = new TH2D(Form("hPtVsPtLeadingMeasured_%s",nameReg[i]),"",ptNbinsL,ptbins1L,ptNbins,ptbins1);
			fOutputList->Add(hPtVsPtLeadingMeasured[i]);


			pNumDenMeasured[i] = new TProfile(Form("pNumDenMeasured_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
			fOutputList->Add(pNumDenMeasured[i]);

			pSumPtMeasured[i] = new TProfile(Form("pSumPtMeasured_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
			fOutputList->Add(pSumPtMeasured[i]);

			// MC true
			// All no trigger sel
			pNumDenTrueAll[i] = new TProfile(Form("pNumDenTrueAll_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
			fOutputList->Add(pNumDenTrueAll[i]);

			pSumPtTrueAll[i] = new TProfile(Form("pSumPtTrueAll_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
			fOutputList->Add(pSumPtTrueAll[i]);
			// trigger sel
			pNumDenTruePS[i] = new TProfile(Form("pNumDenTruePS_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
			fOutputList->Add(pNumDenTruePS[i]);

			pSumPtTruePS[i] = new TProfile(Form("pSumPtTruePS_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
			fOutputList->Add(pSumPtTruePS[i]);
			// trigger sel + vtx rec
			pNumDenTruePSV[i] = new TProfile(Form("pNumDenTruePSV_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
			fOutputList->Add(pNumDenTruePSV[i]);

			pSumPtTruePSV[i] = new TProfile(Form("pSumPtTruePSV_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
			fOutputList->Add(pSumPtTruePSV[i]);
			// all sel criteria
			pNumDenTrueGood[i] = new TProfile(Form("pNumDenTrueGood_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
			fOutputList->Add(pNumDenTrueGood[i]);

			pSumPtTrueGood[i] = new TProfile(Form("pSumPtTrueGood_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
			fOutputList->Add(pSumPtTrueGood[i]);

			pNumDenTrue[i] = new TProfile(Form("pNumDenTrue_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
			fOutputList->Add(pNumDenTrue[i]);

			pSumPtTrue[i] = new TProfile(Form("pSumPtTrue_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
			fOutputList->Add(pSumPtTrue[i]);
		}

		hPtLeadingTrue = new TH1D("hPtLeadingTrue","",ptNbinsL,ptbins1L);
		fOutputList->Add(hPtLeadingTrue);

		for(Int_t i=0;i<3;++i){

			hPtVsUEGenTest[i] = new TH2D(Form("hPtVsUEGenTest_%s",nameReg[i]),"gen pT vs nch_transverse",nchNbins,nchbins,ptNbins,ptbins1);
			fOutputList->Add(hPtVsUEGenTest[i]);

			hPtVsUERecTest[i] = new TH2D(Form("hPtVsUERecTest_%s",nameReg[i]),"rec pT vs nch_transverse",nchNbins,nchbins,ptNbins,ptbins1);
			fOutputList->Add(hPtVsUERecTest[i]);

		}

		hPtLeadingMeasured = new TH1D("hPtLeadingMeasured","",ptNbinsL,ptbins1L);
		fOutputList->Add(hPtLeadingMeasured);

		hPtLeadingGenAll = new TH1D("hPtLeadingGenAll","gen pTleading before any selection",ptNbinsL,ptbins1L);
		fOutputList->Add(hPtLeadingGenAll);

		hPtLeadingGenPS = new TH1D("hPtLeadingGenPS","gen pTleading after physics selection",ptNbinsL,ptbins1L); 
		fOutputList->Add(hPtLeadingGenPS);

		hPtLeadingGenPSV = new TH1D("hPtLeadingGenPSV","gen pTleading after physics selection + vtx",ptNbinsL,ptbins1L); 
		fOutputList->Add(hPtLeadingGenPSV);

		hPtLeadingGenGood = new TH1D("hPtLeadingGenGood","gen pTleading after physics selection + vtx",ptNbinsL,ptbins1L); 
		fOutputList->Add(hPtLeadingGenGood);

	}

	hNchTSRec = new TH1D("hNchTSRec","",100,-0.5,99.5);
	fOutputList->Add(hNchTSRec);

	hNchTSData = new TH1D("hNchTSData","",100,-0.5,99.5); 
	fOutputList->Add(hNchTSData);

	for(Int_t i=0;i<3;++i){
		hPhiRec[i]= new TH1D(Form("hPhiRec_%s",nameReg[i]),"",64,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
		fOutputList->Add(hPhiRec[i]);
	}
	for(Int_t i=0;i<3;++i){

		hPtVsUEData[i] = new TH2D(Form("hPtVsUEData_%s",nameReg[i]),"rec pT vs nch_transverse",nchNbins,nchbins,ptNbins,ptbins1);
		fOutputList->Add(hPtVsUEData[i]);

	}

	hCounter = new TH1D("hCounter","Counter; sel; Nev",3,0,3);
	fOutputList->Add(hCounter);

	// Data driven
	for(Int_t i=0;i<3;++i){
		hNumDenMCDd[i]= new TH2D(Form("hNumDenMCDd_%s",nameReg[i]),"",ptNbinsL,ptbins1L,nchNbins,nchbins);
		fOutputList->Add(hNumDenMCDd[i]);
		hSumPtMCDd[i]= new TH2D(Form("hSumPtMCDd_%s",nameReg[i]),"",ptNbinsL,ptbins1L,ptNbins,ptbins1);
		fOutputList->Add(hSumPtMCDd[i]);
		hNumDenMCMatchDd[i]= new TH2D(Form("hNumDenMCMatchDd_%s",nameReg[i]),"",ptNbinsL,ptbins1L,nchNbins,nchbins);
		fOutputList->Add(hNumDenMCMatchDd[i]);
		hSumPtMCMatchDd[i]= new TH2D(Form("hSumPtMCMatchDd_%s",nameReg[i]),"",ptNbinsL,ptbins1L,ptNbins,ptbins1);
		fOutputList->Add(hSumPtMCMatchDd[i]);

	}


	for(Int_t i=0;i<3;++i){

		hPtVsPtLeadingData[i] = new TH2D(Form("hPtVsPtLeadingData_%s",nameReg[i]),"",ptNbinsL,ptbins1L,ptNbins,ptbins1);
		fOutputList->Add(hPtVsPtLeadingData[i]);


		pNumDenData[i] = new TProfile(Form("pNumDenData_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
		fOutputList->Add(pNumDenData[i]);

		pSumPtData[i] = new TProfile(Form("pSumPtData_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
		fOutputList->Add(pSumPtData[i]);

		// All no trigger sel
		pNumDenMeasuredAll[i] = new TProfile(Form("pNumDenMeasuredAll_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
		fOutputList->Add(pNumDenMeasuredAll[i]);

		pSumPtMeasuredAll[i] = new TProfile(Form("pSumPtMeasuredAll_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
		fOutputList->Add(pSumPtMeasuredAll[i]);
		// trigger sel
		pNumDenMeasuredPS[i] = new TProfile(Form("pNumDenMeasuredPS_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
		fOutputList->Add(pNumDenMeasuredPS[i]);

		pSumPtMeasuredPS[i] = new TProfile(Form("pSumPtMeasuredPS_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
		fOutputList->Add(pSumPtMeasuredPS[i]);
		// trigger sel + vtx rec
		pNumDenMeasuredPSV[i] = new TProfile(Form("pNumDenMeasuredPSV_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
		fOutputList->Add(pNumDenMeasuredPSV[i]);

		pSumPtMeasuredPSV[i] = new TProfile(Form("pSumPtMeasuredPSV_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
		fOutputList->Add(pSumPtMeasuredPSV[i]);
		// all sel criteria
		pNumDenMeasuredGood[i] = new TProfile(Form("pNumDenMeasuredGood_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
		fOutputList->Add(pNumDenMeasuredGood[i]);

		pSumPtMeasuredGood[i] = new TProfile(Form("pSumPtMeasuredGood_%s",nameReg[i]),"",ptNbinsL,ptbins1L);
		fOutputList->Add(pSumPtMeasuredGood[i]);




	}

	hPtLeadingData = new TH1D("hPtLeadingData","",ptNbinsL,ptbins1L);
	fOutputList->Add(hPtLeadingData);

	hPtLeadingRecAll = new TH1D("hPtLeadingRecAll","rec pTleading before any selection",ptNbinsL,ptbins1L); 
	fOutputList->Add(hPtLeadingRecAll);

	hPtLeadingRecPS = new TH1D("hPtLeadingRecPS","rec pTleading after physics selection",ptNbinsL,ptbins1L); 
	fOutputList->Add(hPtLeadingRecPS);

	hPtLeadingRecPSV = new TH1D("hPtLeadingRecPSV","rec pTleading after physics selection + vtx",ptNbinsL,ptbins1L); 
	fOutputList->Add(hPtLeadingRecPSV);

	hPtLeadingRecGood = new TH1D("hPtLeadingRecGood","rec pTleading after physics selection + vtx",ptNbinsL,ptbins1L); 
	fOutputList->Add(hPtLeadingRecGood);

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
		fMCStack = fMC->Stack();
	}


	hCounter->Fill(0);


	AliHeader* headerMC = fMC->Header();
	Bool_t isGoodVtxPosMC = kFALSE;
	if (fUseMC){
		AliGenEventHeader* genHeader = headerMC->GenEventHeader();
		TArrayF vtxMC(3); // primary vertex  MC 
		vtxMC[0]=9999; vtxMC[1]=9999;  vtxMC[2]=9999; //initialize with dummy
		if (genHeader) {
			genHeader->PrimaryVertex(vtxMC);
		}
		if(TMath::Abs(vtxMC[2])<=10)
			isGoodVtxPosMC = kTRUE;

		// Before trigger selection
		GetLeadingObject(kTRUE);// leading particle at gen level
	}
	// Before trigger selection
	GetLeadingObject(kFALSE);// leading particle at rec level
	// Now we get the leading particle pT
	hPtLeadingRecAll->Fill(fRecLeadPt);
	if(isGoodVtxPosMC)
		hPtLeadingGenAll->Fill(fGenLeadPt);

	vector<Double_t> ue_gen;// 0: nch_near, 1: nch_away, 2: nch_trans, 3: sumpt_near, ... 
	vector<Double_t> ue_rec;
	GetMeanUEObservables(ue_gen,ue_rec);

	for(Int_t i=0;i<3;++i){
		pNumDenMeasuredAll[i]->Fill(fRecLeadPt,ue_rec[i]);
		pSumPtMeasuredAll[i]->Fill(fRecLeadPt,ue_rec[i+3]);
		if(isGoodVtxPosMC){
			pNumDenTrueAll[i]->Fill(fGenLeadPt,ue_gen[i]);
			pSumPtTrueAll[i]->Fill(fGenLeadPt,ue_gen[i+3]);
		}
	}

	// Trigger selection
	UInt_t fSelectMask= fInputHandler->IsEventSelected();
	Bool_t isINT7selected = fSelectMask&AliVEvent::kINT7;
	if(!isINT7selected)
		return;
	hCounter->Fill(1);


	hPtLeadingRecPS->Fill(fRecLeadPt);
	if(isGoodVtxPosMC)
		hPtLeadingGenPS->Fill(fGenLeadPt);

	for(Int_t i=0;i<3;++i){
		pNumDenMeasuredPS[i]->Fill(fRecLeadPt,ue_rec[i]);
		pSumPtMeasuredPS[i]->Fill(fRecLeadPt,ue_rec[i+3]);
		if(isGoodVtxPosMC){
			pNumDenTruePS[i]->Fill(fGenLeadPt,ue_gen[i]);
			pSumPtTruePS[i]->Fill(fGenLeadPt,ue_gen[i+3]);
		}
	}

	// Good vertex
	Bool_t hasRecVertex = kFALSE;
	hasRecVertex=HasRecVertex();
	if(!hasRecVertex)return;

	hPtLeadingRecPSV->Fill(fRecLeadPt);
	if(isGoodVtxPosMC)
		hPtLeadingGenPSV->Fill(fGenLeadPt);

	for(Int_t i=0;i<3;++i){
		pNumDenMeasuredPSV[i]->Fill(fRecLeadPt,ue_rec[i]);
		pSumPtMeasuredPSV[i]->Fill(fRecLeadPt,ue_rec[i+3]);
		if(isGoodVtxPosMC){
			pNumDenTruePSV[i]->Fill(fGenLeadPt,ue_gen[i]);
			pSumPtTruePSV[i]->Fill(fGenLeadPt,ue_gen[i+3]);
		}
	}

	// Good events
	if (!fEventCuts.AcceptEvent(event)) {
		PostData(1, fOutputList);
		return;
	}

	hPtLeadingRecGood->Fill(fRecLeadPt);
	if(isGoodVtxPosMC)
		hPtLeadingGenGood->Fill(fGenLeadPt);

	for(Int_t i=0;i<3;++i){
		pNumDenMeasuredGood[i]->Fill(fRecLeadPt,ue_rec[i]);
		pSumPtMeasuredGood[i]->Fill(fRecLeadPt,ue_rec[i+3]);
		if(isGoodVtxPosMC){
			pNumDenTrueGood[i]->Fill(fGenLeadPt,ue_gen[i]);
			pSumPtTrueGood[i]->Fill(fGenLeadPt,ue_gen[i+3]);
		}
	}


	if(fIsMCclosure){
		Double_t randomUE = gRandom->Uniform(0.0,1.0);
		if(randomUE<0.5){// corrections (50% stat.)
			if(isGoodVtxPosMC){
				// KNO scaling
				if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax ))
					GetDetectorResponse();

				// UE analysis
				if(fGenLeadPt>=fPtMin){
					GetBinByBinCorrections();
					GetPtLeadingMisRecCorrection();
				}
			}
		}
		else{// for testing the method
			// KNO scaling
			if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax ))
				GetMultiplicityDistributions();

			// UE analysis
			if(fRecLeadPt>=fPtMin){
				GetUEObservables(); 
			}

		}
	}
	else{
		if(fUseMC){
			if(isGoodVtxPosMC){
				// KNO scaling
				if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax ))
					GetDetectorResponse();

				// UE analysis
				if(fGenLeadPt>=fPtMin){
					GetBinByBinCorrections();
					GetPtLeadingMisRecCorrection();
				}
			}
		}
		else{

			// KNO scaling
			if(( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax ))
				GetMultiplicityDistributionsData();

			// UE analysis
			if(fRecLeadPt>=fPtMin){
				GetPtLeadingMisRecCorrection();
				GetUEObservablesData();
			}

		}
	}



	ue_gen.clear();
	ue_rec.clear();

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
void AliAnalysisTaskMcKnoUe::GetPtLeadingMisRecCorrection(){

	Int_t nch_top[3];
	Double_t sumpt_top[3];
	for(Int_t i=0;i<3;++i){
		nch_top[i]=0;
		sumpt_top[i]=0;
	}

	vector<Float_t> ptArray;
	vector<Float_t> phiArray;
	vector<Int_t>   indexArray;


	Int_t iTracks(fESD->GetNumberOfTracks());           // see how many tracks there are in the event
	for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks

		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

		if(!track) continue;

		if(!fLeadingTrackFilter->IsSelected(track))
			continue;

		if(TMath::Abs(track->Eta()) > fEtaCut)
			continue;

		if( track->Pt() < fPtMin)continue;

		if(i!=fRecLeadIn){// here we exclude the auto correlation
			Double_t DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

			// definition of the topological regions
			if(TMath::Abs(DPhi)<pi/3.0){// near side
				nch_top[0]++; sumpt_top[0]+=track->Pt();	
			}
			else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
				nch_top[1]++; sumpt_top[1]+=track->Pt();
			}
			else{// transverse side
				nch_top[2]++; sumpt_top[2]+=track->Pt();
			}
		}
		// second track selection following the efficiency 
		if( f_Eff->Eval(track->Pt()) < gRandom->Uniform(0,1) )
			continue;
		ptArray.push_back(track->Pt());
		phiArray.push_back(track->Phi());
		indexArray.push_back(i);

	}

	if(fUseMC){
		AliESDtrack* ltrack = static_cast<AliESDtrack*>(fESD->GetTrack(fRecLeadIn));
		const Int_t label = TMath::Abs(ltrack->GetLabel());

		Double_t ptlrec = ltrack->Pt();
		for(Int_t i=0;i<3;++i){
			hNumDenMC[i]->Fill(ptlrec,nch_top[i]);
			hSumPtMC[i]->Fill(ptlrec,sumpt_top[i]);
			if(label==fGenLeadIn){
				hNumDenMCMatch[i]->Fill(ptlrec,nch_top[i]);
				hSumPtMCMatch[i]->Fill(ptlrec,sumpt_top[i]);
			}
		}
	}
	// Now the data driven approach
	Float_t flPt = 0;// leading pT
	Float_t flPhi = 0;
	Int_t flIndex = 0;
	Int_t ntrk = ptArray.size();

	for(Int_t i=0;i<ntrk;++i){


		if ( flPt < ptArray[i] ){
			flPt  = ptArray[i];
			flPhi = phiArray[i];
			flIndex = indexArray[i];
		}
	}

	Int_t nchm_top[3];
	Double_t sumptm_top[3];
	for(Int_t i=0;i<3;++i){
		nchm_top[i]=0;
		sumptm_top[i]=0;
	}
	for(Int_t i=0;i<ntrk;++i){

		if(indexArray[i]==flIndex)
			continue;

		Double_t DPhi = DeltaPhi(phiArray[i], flPhi);
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			nchm_top[0]++; sumptm_top[0]+=ptArray[i];
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			nchm_top[1]++; sumptm_top[1]+=ptArray[i];
		}
		else{// transverse side
			nchm_top[2]++; sumptm_top[2]+=ptArray[i];
		}
	}
	// Here I fill the histograms for the data driven (Dd) aprroach
	for(Int_t i=0;i<3;++i){
		hNumDenMCDd[i]->Fill(flPt,nchm_top[i]);
		hSumPtMCDd[i]->Fill(flPt,sumptm_top[i]);
		if(flIndex==fRecLeadIn){
			hNumDenMCMatchDd[i]->Fill(flPt,nchm_top[i]);
			hSumPtMCMatchDd[i]->Fill(flPt,sumptm_top[i]);
		}
	}

	ptArray.clear();
	phiArray.clear();
	indexArray.clear();

}


void AliAnalysisTaskMcKnoUe::GetBinByBinCorrections(){

	// Histos for efficiencyxacceptance
	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {

		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;

		if (!fMC->IsPhysicalPrimary(i)) continue; 
		if (particle->Charge() == 0) continue;
		if ( TMath::Abs(particle->Eta()) > fEtaCut )continue;
		if( particle->Pt() < fPtMin)continue;
		hPtInPrim->Fill(particle->Pt());// inital pT distribution (MC gen)

	}


	Int_t iTracks(fESD->GetNumberOfTracks());           // see how many tracks there are in the event
	for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks

		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

		if(!track) continue;

		if(!fLeadingTrackFilter->IsSelected(track))// for traditional UE analysis we consider TPCITS2015
			continue;

		if(TMath::Abs(track->Eta()) > fEtaCut)
			continue;

		if( track->Pt() < fPtMin)continue;

		const Int_t label = TMath::Abs(track->GetLabel());
		hPtOut->Fill(track->Pt());
		if( fMCStack->IsPhysicalPrimary(label) ){
			hPtOutPrim->Fill(track->Pt());
		}
		if( fMCStack->IsSecondaryFromWeakDecay(label) || fMCStack->IsSecondaryFromMaterial(label)){
			hPtOutSec->Fill(track->Pt());
		}


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

void AliAnalysisTaskMcKnoUe::GetMeanUEObservables(vector<Double_t> &genArray, vector<Double_t> &recArray){

	genArray.clear();
	recArray.clear();

	Int_t nch_top[3];
	Double_t sumpt_top[3];
	for(Int_t i=0;i<3;++i){
		nch_top[i]=0;
		sumpt_top[i]=0;
	}

	if (fUseMC){


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
				nch_top[0]++; sumpt_top[0]+=particle->Pt();
			}
			else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
				nch_top[1]++; sumpt_top[1]+=particle->Pt();
			}
			else{// transverse side
				nch_top[2]++; sumpt_top[2]+=particle->Pt();
			}



		}
	}
	Int_t nchm_top[3];
	Double_t sumptm_top[3];
	for(Int_t i=0;i<3;++i){
		nchm_top[i]=0;
		sumptm_top[i]=0;
	}


	Int_t iTracks(fESD->GetNumberOfTracks());           // see how many tracks there are in the event
	for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks

		if(i==fRecLeadIn)
			continue;

		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

		if(!track) continue;

		if(!fLeadingTrackFilter->IsSelected(track))
			continue;

		if(TMath::Abs(track->Eta()) > fEtaCut)
			continue;

		if( track->Pt() < fPtMin)continue;

		Double_t DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		// definition of the topological regions
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			nchm_top[0]++; sumptm_top[0]+=track->Pt();
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			nchm_top[1]++; sumptm_top[1]+=track->Pt();
		}
		else{// transverse side
			nchm_top[2]++; sumptm_top[2]+=track->Pt();
		}
	}

	if (fUseMC){
		for(Int_t i=0;i<3;++i)
			genArray.push_back(1.0*nch_top[i]);
		for(Int_t i=0;i<3;++i)
			genArray.push_back(sumpt_top[i]);
	}
	for(Int_t i=0;i<3;++i)
		recArray.push_back(1.0*nchm_top[i]);
	for(Int_t i=0;i<3;++i)
		recArray.push_back(sumptm_top[i]);


}

void AliAnalysisTaskMcKnoUe::GetUEObservablesData(){

	Int_t nchm_top[3];
	Double_t sumptm_top[3];
	for(Int_t i=0;i<3;++i){
		nchm_top[i]=0;
		sumptm_top[i]=0;
	}


	Int_t iTracks(fESD->GetNumberOfTracks());           // see how many tracks there are in the event
	for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks

		if(i==fRecLeadIn)
			continue;

		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

		if(!track) continue;

		if(!fLeadingTrackFilter->IsSelected(track))
			continue;

		if(TMath::Abs(track->Eta()) > fEtaCut)
			continue;

		if( track->Pt() < fPtMin)continue;

		Double_t DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		// definition of the topological regions
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			nchm_top[0]++; sumptm_top[0]+=track->Pt();
			hPtVsPtLeadingData[0]->Fill(fRecLeadPt,track->Pt());
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			nchm_top[1]++; sumptm_top[1]+=track->Pt();
			hPtVsPtLeadingData[1]->Fill(fRecLeadPt,track->Pt());
		}
		else{// transverse side
			nchm_top[2]++; sumptm_top[2]+=track->Pt();
			hPtVsPtLeadingData[2]->Fill(fRecLeadPt,track->Pt());
		}

	}

	for(Int_t i=0;i<3;++i){
		pNumDenData[i]->Fill(fRecLeadPt,nchm_top[i]);
		pSumPtData[i]->Fill(fRecLeadPt,sumptm_top[i]);
	}
	hPtLeadingData->Fill(fRecLeadPt);

}

void AliAnalysisTaskMcKnoUe::GetUEObservables(){

	Int_t nch_top[3];
	Double_t sumpt_top[3];
	for(Int_t i=0;i<3;++i){
		nch_top[i]=0;
		sumpt_top[i]=0;
	}


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
			nch_top[0]++; sumpt_top[0]+=particle->Pt();
			hPtVsPtLeadingTrue[0]->Fill(fGenLeadPt,particle->Pt());
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			nch_top[1]++; sumpt_top[1]+=particle->Pt();
			hPtVsPtLeadingTrue[1]->Fill(fGenLeadPt,particle->Pt());
		}
		else{// transverse side
			nch_top[2]++; sumpt_top[2]+=particle->Pt();
			hPtVsPtLeadingTrue[2]->Fill(fGenLeadPt,particle->Pt());
		}



	}

	Int_t nchm_top[3];
	Double_t sumptm_top[3];
	for(Int_t i=0;i<3;++i){
		nchm_top[i]=0;
		sumptm_top[i]=0;
	}


	Int_t iTracks(fESD->GetNumberOfTracks());           // see how many tracks there are in the event
	for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks

		if(i==fRecLeadIn)
			continue;

		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

		if(!track) continue;

		if(!fLeadingTrackFilter->IsSelected(track))
			continue;

		if(TMath::Abs(track->Eta()) > fEtaCut)
			continue;

		if( track->Pt() < fPtMin)continue;

		Double_t DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		// definition of the topological regions
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			nchm_top[0]++; sumptm_top[0]+=track->Pt();
			hPtVsPtLeadingMeasured[0]->Fill(fRecLeadPt,track->Pt());
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			nchm_top[1]++; sumptm_top[1]+=track->Pt();
			hPtVsPtLeadingMeasured[1]->Fill(fRecLeadPt,track->Pt());
		}
		else{// transverse side
			nchm_top[2]++; sumptm_top[2]+=track->Pt();
			hPtVsPtLeadingMeasured[2]->Fill(fRecLeadPt,track->Pt());
		}

	}

	for(Int_t i=0;i<3;++i){
		pNumDenMeasured[i]->Fill(fRecLeadPt,nchm_top[i]);
		pSumPtMeasured[i]->Fill(fRecLeadPt,sumptm_top[i]);

		pNumDenTrue[i]->Fill(fGenLeadPt,nch_top[i]);
		pSumPtTrue[i]->Fill(fGenLeadPt,sumpt_top[i]);
	}
	hPtLeadingTrue->Fill(fGenLeadPt);
	hPtLeadingMeasured->Fill(fRecLeadPt);



}
//______________________________________________________________
void AliAnalysisTaskMcKnoUe::GetMultiplicityDistributionsData(){

	Int_t multTSrec=0;

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
	hNchTSData->Fill(multTSrec); 

	// Filling rec pT vs UE (for pT I use 2015 track cuts, UE uses TPC-only)
	for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks

		if(i==fRecLeadIn)
			continue;

		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

		if(!track) continue;

		if(!fLeadingTrackFilter->IsSelected(track))
			continue;

		if(TMath::Abs(track->Eta()) > fEtaCut)
			continue;

		if( track->Pt() < fPtMin)continue;

		Double_t DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		// definition of the topological regions
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			hPtVsUEData[0]->Fill(multTSrec,track->Pt());
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			hPtVsUEData[1]->Fill(multTSrec,track->Pt());
		}
		else{// transverse side
			hPtVsUEData[2]->Fill(multTSrec,track->Pt());
		}

	}

}
//____________________________________________________________
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
	// Filling pT vs UE activity
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
			hPtVsUEGenTest[0]->Fill(multTSgen,particle->Pt());
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			hPtVsUEGenTest[1]->Fill(multTSgen,particle->Pt());
		}
		else{// transverse side
			hPtVsUEGenTest[2]->Fill(multTSgen,particle->Pt());
		}
	}


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

	// Filling rec pT vs UE (for pT I use 2015 track cuts, UE uses TPC-only)
	for(Int_t i=0; i < iTracks; i++) {                 // loop over all these tracks

		if(i==fRecLeadIn)
			continue;

		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  // get a track (type AliesdTrack)

		if(!track) continue;

		if(!fLeadingTrackFilter->IsSelected(track))
			continue;

		if(TMath::Abs(track->Eta()) > fEtaCut)
			continue;

		if( track->Pt() < fPtMin)continue;

		Double_t DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		// definition of the topological regions
		if(TMath::Abs(DPhi)<pi/3.0){// near side
			hPtVsUERecTest[0]->Fill(multTSrec,track->Pt());
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){// away side
			hPtVsUERecTest[1]->Fill(multTSrec,track->Pt());
		}
		else{// transverse side
			hPtVsUERecTest[2]->Fill(multTSrec,track->Pt());
		}

	}


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
Bool_t AliAnalysisTaskMcKnoUe::HasRecVertex(){


	float fMaxDeltaSpdTrackAbsolute = 0.5f;
	float fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
	float fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
	float fMaxResolutionSPDvertex = 0.25f;
	float fMaxDispersionSPDvertex = 1.e14f;

	Bool_t fRequireTrackVertex = true;
	unsigned long fFlag;
	fFlag =BIT(AliEventCuts::kNoCuts);

	const AliVVertex* vtTrc = fESD->GetPrimaryVertex();
	bool isTrackV = true;
	if(vtTrc->IsFromVertexer3D() || vtTrc->IsFromVertexerZ()) isTrackV=false;
	const AliVVertex* vtSPD = fESD->GetPrimaryVertexSPD();


	if (vtSPD->GetNContributors() > 0) fFlag |= BIT(AliEventCuts::kVertexSPD);

	if (vtTrc->GetNContributors() > 1 && isTrackV ) fFlag |= BIT(AliEventCuts::kVertexTracks);

	if (((fFlag & BIT(AliEventCuts::kVertexTracks)) ||  !fRequireTrackVertex) && (fFlag & BIT(AliEventCuts::kVertexSPD))) fFlag |= BIT(AliEventCuts::kVertex);

	const AliVVertex* &vtx = bool(fFlag & BIT(AliEventCuts::kVertexTracks)) ? vtTrc : vtSPD;
	AliVVertex   *fPrimaryVertex = const_cast<AliVVertex*>(vtx);
	if(!fPrimaryVertex)return kFALSE;

	/// Vertex quality cuts
	double covTrc[6],covSPD[6];
	vtTrc->GetCovarianceMatrix(covTrc);
	vtSPD->GetCovarianceMatrix(covSPD);
	double dz = bool(fFlag & AliEventCuts::kVertexSPD) && bool(fFlag & AliEventCuts::kVertexTracks) ? vtTrc->GetZ() - vtSPD->GetZ() : 0.; /// If one of the two vertices is not available this cut is always passed.
	double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
	double errTrc = bool(fFlag & AliEventCuts::kVertexTracks) ? TMath::Sqrt(covTrc[5]) : 1.;
	double nsigTot = TMath::Abs(dz) / errTot, nsigTrc = TMath::Abs(dz) / errTrc;
	/// vertex dispersion for run1, only for ESD, AOD code to be added here
	const AliESDVertex* vtSPDESD = dynamic_cast<const AliESDVertex*>(vtSPD);
	double vtSPDdispersion = vtSPDESD ? vtSPDESD->GetDispersion() : 0;
	if (
			(TMath::Abs(dz) <= fMaxDeltaSpdTrackAbsolute && nsigTot <= fMaxDeltaSpdTrackNsigmaSPD && nsigTrc <= fMaxDeltaSpdTrackNsigmaTrack) && // discrepancy track-SPD vertex
			(!vtSPD->IsFromVertexerZ() || TMath::Sqrt(covSPD[5]) <= fMaxResolutionSPDvertex) &&
			(!vtSPD->IsFromVertexerZ() || vtSPDdispersion <= fMaxDispersionSPDvertex) /// vertex dispersion cut for run1, only for ESD
	   ) // quality cut on vertexer SPD z
		fFlag |= BIT(AliEventCuts::kVertexQuality);  

	Bool_t hasVtx = (TESTBIT(fFlag,AliEventCuts::kVertex))&&(TESTBIT(fFlag,AliEventCuts::kVertexQuality));

	return hasVtx;

}

