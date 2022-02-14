/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *

 *                                                                        *
 * Author: Zhongbao Yin                                                  *
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

/* The task selects candidates for V0s and Cascade (as associated particles)// updated by (Donghai Liu & Mustafa Anaam) 2021
 * and calculates correlations with unidentified charged trigger particles in phi and eta. 
 * The task works with AOD events only and containes also mixing for acceptance corrections.
 */
#include <iostream>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODVertex.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliEventPoolManager.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskLambdaK0s.h"
#include "AliPhysicsSelectionTask.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliStack.h"
#include "AliAnalysisUtils.h"
#include "AliEventCuts.h"
#include "AliMultSelection.h"
#include <AliMultiInputEventHandler.h>
#include <AliMixInputEventHandler.h>



using namespace std;
	ClassImp(AliAnalysisTaskLambdaK0s)
ClassImp(AliV0XiParticleall)

	//_____________________________________________________________________
	AliAnalysisTaskLambdaK0s::AliAnalysisTaskLambdaK0s() : 
		AliAnalysisTaskSE(),
		fMixingTracks(200),
		fPoolSize(500),
		fPoolMgr(0x0),
		fEffList(0),
		fHistEffEtaPtK0s(0),
		fHistEffEtaPtLambda(0),
		fHistEffEtaPtAntiLambda(0),
		fHistEffEtaPtXiMinus(0),
		fOutput(0),
		fPIDResponse(0),
		fHistArmenterosPodolanski(0),
		fHistRCTrkPtEta(0),
		fHistRCPriTrkPtEta(0),
		fHistRCSecTrkPtEta(0),
		fHistMCtruthTrkPtEta(0),
		fHistRCK0sPt(0),
		fHistMCtruthK0sPt(0),
		fHistRCLambdaPt(0),
		fHistallRCLambdaPt(0),
		fHistRCXiMinusPt(0),
		fHistMCtruthXiMinusPt(0),
		fHistRCXiPt(0),
		HistRCXiPt(0),
		fHistMCXiMinusPt(0),
		fHistXiPtVSLambdaPt(0),
		fHistXiPtVSAntiLambdaPt(0),
		fRejectTrackPileUp(kTRUE),
		fIsPileUpCuts(kTRUE),
		fHistallDLLambdaPt(0),
		fHistallDCALambdaPt(0),
		fHistallRLambdaPt(0),
		fHistMCtruthLambdaPt(0),
		fHistRCAntiLambdaPt(0),
		fHistMCtruthAntiLambdaPt(0),
		fHistRCPrimLambdaDCAPt(0),
		fHistRCPrimLambdaDLPt(0),
		fHistRCPrimLambdaCTPt(0), 
		fHistRCPrimLambdaRPt(0),
		fHistRCSecLambdaDCAPt(0),
		fHistRCSecLambdaDLPt(0),
		fHistRCSecLambdaCTPt(0),
		fHistRCSecLambdaRPt(0),
		fHistRCPrimAntiLambdaDCAPt(0),
		fHistRCSecAntiLambdaDCAPt(0),
		fEventCuts(0x0),
		fEffCorr(kFALSE),
		fChXi(kFALSE),
		fChV0(kFALSE),
		fAnalysisMC(kFALSE),
		fEventMixing(kFALSE),
		fUseHybridGlobalTrack(kFALSE),
		fUtils(NULL),
		fMinCent(0),
		fMaxCent(0),
		fPrimVertexCut(0),
		fDCAToPrimVtx(0),
		fDCABetweenDaughters(0),
		fDCANegtoPrimVertexMinK0s(0),
		fDCAPostoPrimVertexMinK0s(0),
		fDCANegtoPrimVertexMinLambda(0),
		fDCAPostoPrimVertexMinLambda(0),
		fDCANegtoPrimVertexMinAntiLambda(0),
		fDCAPostoPrimVertexMinAntiLambda(0),
		fCPA(0),
		fPLtimeK0s(0),
		fPLtimeLambda(0),
		fLambdaCPA(0),
		fLambdaDCA(0),
		f2DFiducial(0),
		fNsigma(0), 
		fTPCNcls(0), 
		fMinCtau(0),
		fMaxCtau(0),
		fEtaCut(0),
		fV0Eta(0),
		fRapidity(0),
		fV0DaughterEtaCut(0),
		fV0DaughterPtMinCut(0),
		fPtAssoMin(0), 
		fPtAssoMax(0),
		fPtTrigMin(0),
		fPtTrigMax (0), 
		fMCArray(NULL),
		selectedTracks(NULL),
		selectedK0s(NULL),
		selectedLambda(NULL),
		selectedAntiLambda(NULL),
		selectedXiMinus(NULL),
		assocParticles(NULL)    

{
	fMassMean[0] = 0.497614; fMassMean[1] = 1.115683; 
	fMassRes[0] = 0.005; fMassRes[1] = 0.0025; 

	for(Int_t i = 0; i < kNVtxZ; i++){

		fHistRCTrk[i] = NULL;
		fHistRCPriTrk[i] = NULL;
		fHistRCSecTrk[i] = NULL;
		fHistMCtruthTrk[i] = NULL;
		fHistRCK0s[i] = NULL;
		fHistMCtruthK0s[i] = NULL;
		fHistRCLambda[i] = NULL;
		fHistMCtruthLambda[i] = NULL;
		fHistRCAntiLambda[i] = NULL;
		fHistMCtruthAntiLambda[i] = NULL;

	}


	for(Int_t i = 0; i < 3; i++){
		fBestPrimaryVtxPos[i] = -9999;
	}

}


//________________________________________________________________________
AliAnalysisTaskLambdaK0s::AliAnalysisTaskLambdaK0s(const char *name,Double_t centMin,	Double_t centMax,
		Bool_t effCorr)// All data members should be initialised here
	: AliAnalysisTaskSE(name),
	fMixingTracks(200),
	fPoolSize(500), 
	fPoolMgr(0x0),
	fEffList(0),
	fHistEffEtaPtK0s(0),
	fHistEffEtaPtLambda(0),
	fHistEffEtaPtAntiLambda(0),
	fHistEffEtaPtXiMinus(0),
	fOutput(0),
	fPIDResponse(0),
	fHistArmenterosPodolanski(0),     
	fHistRCSecTrkPtEta(0),
	fHistMCtruthTrkPtEta(0),
	fHistRCK0sPt(0),
	fHistMCtruthK0sPt(0),
	fHistRCLambdaPt(0),
	fHistMCtruthLambdaPt(0),
	fHistallRCLambdaPt(0),
	fHistallDLLambdaPt(0),
	fHistallDCALambdaPt(0),
	fHistallRLambdaPt(0),
	fHistRCAntiLambdaPt(0),
	fHistMCtruthAntiLambdaPt(0),
	fHistRCPrimLambdaDCAPt(0),
	fHistRCPrimLambdaDLPt(0),
	fHistRCPrimLambdaCTPt(0),
	fHistRCPrimLambdaRPt(0),
	fHistRCSecLambdaDCAPt(0),
	fHistRCSecLambdaDLPt(0),
	fHistRCSecLambdaCTPt(0),
	fHistRCSecLambdaRPt(0),
	fHistRCXiMinusPt(0),
	fHistMCtruthXiMinusPt(0),
	fHistRCXiPt(0),
	HistRCXiPt(0),
	fHistMCXiMinusPt(0),
	fHistXiPtVSLambdaPt(0),
	fHistXiPtVSAntiLambdaPt(0),
	fIsPileUpCuts(kTRUE),
	fRejectTrackPileUp(kTRUE),
	fHistRCPrimAntiLambdaDCAPt(0),
	fHistRCSecAntiLambdaDCAPt(0),
	fEventCuts(0x0), 
	fEffCorr(effCorr),
	fChXi(kFALSE),
	fChV0(kFALSE),
	fAnalysisMC(kFALSE),
	fEventMixing(kFALSE),
	fUseHybridGlobalTrack(kFALSE),
	fUtils(NULL),
	fMinCent(centMin),
	fMaxCent(centMax),
	fPrimVertexCut(0),
	fDCAToPrimVtx(0),
	fDCABetweenDaughters(0),
	fDCANegtoPrimVertexMinK0s(0),
	fDCAPostoPrimVertexMinK0s(0),
	fDCANegtoPrimVertexMinLambda(0),
	fDCAPostoPrimVertexMinLambda(0),
	fDCANegtoPrimVertexMinAntiLambda(0),
	fDCAPostoPrimVertexMinAntiLambda(0),
	fCPA(0),
	fPLtimeK0s(0),
	fPLtimeLambda(0),
	fLambdaCPA(0),
	fLambdaDCA(0),
	f2DFiducial(0),
	fNsigma(0), 
	fTPCNcls(0), 
	fMinCtau(0),
	fMaxCtau(0),
	fEtaCut(0),
	fV0Eta(0),
	fRapidity(0),
	fV0DaughterEtaCut(0),
	fV0DaughterPtMinCut(0),
	fPtAssoMin(0), 
	fPtAssoMax(0),
	fPtTrigMin(0),
	fPtTrigMax (0), 
	fMCArray(NULL),
	selectedTracks(NULL),
	selectedK0s(NULL),
	selectedLambda(NULL),
	selectedAntiLambda(NULL),
	selectedXiMinus(NULL),
	assocParticles(NULL)

{
	// Constructor
	// Define input and output slots here (never in the dummy constructor)
	// Input slot #0 works with a TChain - it is connected to the default input container
	// Output slot #1 writes into a TH1 container


	fMassMean[0] = 0.497614; fMassMean[1] = 1.115683; 
	fMassRes[0] = 0.005; fMassRes[1] = 0.0025;

	for(Int_t i = 0; i < kNVtxZ; i++){
		fHistRCPriTrk[i] = NULL;
		fHistRCSecTrk[i] = NULL;
		fHistMCtruthTrk[i] = NULL;
		fHistRCK0s[i] = NULL;
		fHistMCtruthK0s[i] = NULL;
		fHistRCLambda[i] = NULL;
		fHistMCtruthLambda[i] = NULL;
		fHistRCAntiLambda[i] = NULL;
		fHistMCtruthAntiLambda[i] = NULL;

	}


	for(Int_t i = 0; i < 3; i++){
		fBestPrimaryVtxPos[i] = -9999;
	}

	//DefineInput(0, TChain::Class());
	if(fEffCorr)
		DefineInput(1, TList::Class());
	DefineOutput(1, TList::Class());     // for output list
}

//________________________________________________________________________
AliAnalysisTaskLambdaK0s::~AliAnalysisTaskLambdaK0s()
{
	// Destructor. Clean-up the output list, but not the histograms that are put inside
	// (the list is owner and will clean-up these histograms). Protect in PROOF case.
	if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
		delete fOutput;
	}
	if(fEventCuts)
	{
		delete fEventCuts;
		fEventCuts = 0x0;
	}

	if(fUtils) delete fUtils;

	if(selectedTracks){delete selectedTracks;}

	if(selectedK0s){delete selectedK0s;}

	if(selectedLambda){delete selectedLambda;}

	if(selectedAntiLambda){delete selectedAntiLambda;}
	if(selectedXiMinus)   {	delete selectedXiMinus;}


}
//________________________________________________________________________
void AliAnalysisTaskLambdaK0s::UserCreateOutputObjects()
{

	if(fEffCorr){
		fEffList = dynamic_cast<TList*>(GetInputData(1));
		if(fEffList){
			fHistEffEtaPtK0s = (TH2F*)fEffList->FindObject("fHistEffEtaPtK0sCent0_10All");
			fHistEffEtaPtLambda = (TH2F*)fEffList->FindObject("fHistEffEtaPtLambdaCent0_10All");
			fHistEffEtaPtAntiLambda = (TH2F*)fEffList->FindObject("fHistEffEtaPtAntiLambdaCent0_10All");
			if(fChXi){
				fHistEffEtaPtXiMinus  = (TH2F*)fEffList->FindObject("fHistEffEtaPtXiMinusCent0_10All");
				if(!fHistEffEtaPtK0s || !fHistEffEtaPtLambda || !fHistEffEtaPtAntiLambda ||!fHistEffEtaPtXiMinus){
					std::cout<<"Efficiency histograms are not available!"<<std::endl;
				}
			} 
			else{   
				if(!fHistEffEtaPtK0s || !fHistEffEtaPtLambda || !fHistEffEtaPtAntiLambda){
					std::cout<<"Efficiency histograms are not available!"<<std::endl;
				}
			}

		}
	}
	fEventCuts = new AliEventCuts();
	// Create histograms for output
	fOutput = new TList();
	fOutput->SetOwner();  // IMPORTANT!
	TList *tEventCutQA = new TList();
	tEventCutQA->SetName("EventCutQA");
	tEventCutQA->SetOwner();
	fEventCuts->fUseVariablesCorrelationCuts = true;
	fEventCuts->AddQAplotsToList(tEventCutQA,true);
	fOutput->Add(tEventCutQA);

	AddQAEvent();
	AddAnalysisTrk();
	AddAnalysisK0s();
	AddAnalysisLambda();
	AddAnalysisAntiLambda();
	AddAnalysisXiMinus();
	AddQAV0Candidates();

	// defining bins for centrality
	const Int_t nCentralityBins  = 1;
	Double_t centBins[] = {fMinCent , fMaxCent};//centMin,centMax

	//-----Vtx bins 
	const Int_t nZvtxBins  = 6;
	Double_t vertexBins[] = {-8., -4., -2., 0., 2., 4., 8.};

	const Int_t nPtBinsV0Xi = 27;
	const Double_t PtBinsV0Xi[28] = {0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5, 5.5, 6, 7, 8, 9, 10};

	const Int_t nPtBinsHXi = 27;
	const Double_t PtBinsHXi[28] = {0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5, 5.5, 6, 7, 8, 9, 10};  

	// defining bins for mass distributions
	Int_t nBins = 160;
	Double_t mMassMin[] ={0.40, 1.04};
	Double_t mMassMax[] ={0.58, 1.21};

	// defining bins for dPhi distributions
	const Int_t ndPhiBins = 32;
	const Double_t kPi = TMath::Pi();
	Double_t PhiMin = -kPi/2.;
	Double_t PhiMax = -kPi/2. + 2.*kPi;

	// defining bins for dEta distributions
	const Int_t ndEtaBins = 32;
	Double_t EtaMin = -2.*fEtaCut;
	Double_t EtaMax = 2.*fEtaCut;


	const Int_t nEtaBins = 18;          
	Double_t EtaBins[nEtaBins+1] = {0.};
	EtaBins[0] = -fEtaCut;
	for (Int_t i=0; i<nEtaBins; i++) {
		EtaBins[i+1] = EtaBins[i] + 2.*fEtaCut/nEtaBins; }

	const Int_t nPhiBins = 36; 
	Double_t PhiBins[nPhiBins+1] = {0.};
	PhiBins[0] = 0;
	for (Int_t i=0; i<nPhiBins; i++) 
	{ PhiBins[i+1] = PhiBins[i] + 2.*TMath::Pi()/nPhiBins; }


	fHistArmenterosPodolanski  = new TH2F("fHistArmenterosPodolanski",
			"Armenteros-Podolanski phase space;#alpha;p_{t} arm",100,-1.0,1.0,50,0,0.5);

	if(fAnalysisMC){
		fHistRCTrkPtEta = new TH2F("fHistRCTrkPtEta",
				"pt and eta of reconstructed tracks; p_{T} [GeV/c]; #eta",
				30, 0, 15, nEtaBins, -fEtaCut, fEtaCut);

		fHistRCPriTrkPtEta = new TH2F("fHistRCPriTrkPtEta",
				"pt and eta of reconstructed prim. tracks; p_{T} [GeV/c]; #eta",
				30, 0, 15, nEtaBins, -fEtaCut, fEtaCut);

		fHistRCSecTrkPtEta = new TH2F("fHistRCSecTrkPtEta",
				"pt and eta of reconstructed sec. tracks; p_{T} [GeV/c]; #eta",
				30,0,15,nEtaBins, -fEtaCut, fEtaCut);

		fHistMCtruthTrkPtEta = new TH2F("fHistMCtruthTrkPtEta",
				"pt and eta of generated tracks; p_{T} [GeV/c]; #eta",
				30,0,15, nEtaBins, -fEtaCut, fEtaCut);

		fHistRCK0sPt = new TH2F("fHistRCK0sPt", "pt of reconstructed K0s; p_{T} [GeV/c]; #eta",
				nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins);

		fHistMCtruthK0sPt = new TH2F("fHistMCtruthK0sPt", "pt of generated K0s; p_{T} [GeV/c]; #eta",
				nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins);


		fHistRCLambdaPt = new TH2F("fHistRCLambdaPt", "pt of reconstructed Lambda; p_{T} [GeV/c]; #eta",
				nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins);

		fHistallRCLambdaPt = new TH2F("fHistallRCLambdaPt", "pt of reconstructed all Lambda;CT (cm); p_{T} [GeV/c] ", 80,0.,23.,nPtBinsV0Xi, PtBinsV0Xi);


		fHistallDLLambdaPt = new TH2F("fHistallDLLambdaPt", "pt of reconstructed all DL Lambda;DL(cm); p_{T} [GeV/c] ", 80,0.,23.,nPtBinsV0Xi, PtBinsV0Xi);
		fHistallDCALambdaPt = new TH2F("fHistallDCALambdaPt", "pt of reconstructed all DCA Lambda;CT (cm); p_{T} [GeV/c] ", 80,0.,1.,nPtBinsV0Xi, PtBinsV0Xi);
		fHistallRLambdaPt = new TH2F("fHistallRLambdaPt", "pt of reconstructed all R Lambda;R (cm); p_{T} [GeV/c] ", 80,5.,100.,nPtBinsV0Xi, PtBinsV0Xi);
		fHistRCAntiLambdaPt 
			= new TH2F("fHistRCAntiLambdaPt", "pt of reconstructed #bar{#Lambda}; p_{T} [GeV/c]; #eta",
					nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins);


		fHistMCtruthLambdaPt 
			= new TH2F("fHistMCtruthLambdaPt", "pt of generated #Lambda; p_{T} [GeV/c]; #eta",
					nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins);

		fHistMCtruthAntiLambdaPt = new TH2F("fHistMCtruthAntiLambdaPt", 
				"pt of generated #bar{#Lambda}; p_{T} [GeV/c]; #eta",
				nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins);

		fHistRCPrimLambdaDCAPt = new TH2F("fHistRCPrimLambdaDCAPt",  "DCA to PV of reco. prim. #Lambda; DCA (cm); p_{T} (GeV/c)", 
				80, 0., 1., nPtBinsV0Xi, PtBinsV0Xi);

		fHistRCPrimLambdaDLPt = new TH2F("fHistRCPrimLambdaDLPt",  "DL of reco. prim. #Lambda; DL (cm); p_{T} (GeV/c)", 
				80, 0., 23., nPtBinsV0Xi, PtBinsV0Xi);
		fHistRCPrimLambdaCTPt = new TH2F("fHistRCPrimLambdaCTPt",  "CT  of reco. prim. #Lambda;CT (cm) ; p_{T} (GeV/c)", 
				80, 0., 23., nPtBinsV0Xi, PtBinsV0Xi);
		fHistRCPrimLambdaRPt = new TH2F("fHistRCPrimLambdaRPt",  "fiducial volume of reco. prim. #Lambda; R(cm); p_{T} (GeV/c)", 
				80, 5., 100., nPtBinsV0Xi, PtBinsV0Xi);


		fHistRCSecLambdaDCAPt = new TH2F("fHistRCSecLambdaDCAPt", "DCA to PV of reco. sec. #Lambda; DCA (cm); p_{T} (GeV/c)",80, 0., 1., nPtBinsV0Xi, PtBinsV0Xi);

		fHistRCSecLambdaDLPt = new TH2F("fHistRCSecLambdaDLPt",  "DL of reco. Sec. #Lambda; DL (cm); p_{T} (GeV/c)", 
				80, 0., 23., nPtBinsV0Xi, PtBinsV0Xi);
		fHistRCSecLambdaCTPt = new TH2F("fHistRCSecLambdaCTPt",  "CT  of reco. Sec. #Lambda;CT (cm) ; p_{T} (GeV/c)", 
				80, 0., 23., nPtBinsV0Xi, PtBinsV0Xi);
		fHistRCSecLambdaRPt = new TH2F("fHistRCSecLambdaRPt",  "fiducial volume of reco. Sec. #Lambda; R(cm); p_{T} (GeV/c)", 
				80, 5., 100., nPtBinsV0Xi, PtBinsV0Xi);


		fHistRCPrimAntiLambdaDCAPt = new TH2F("fHistRCPrimAntiLambdaDCAPt",
				"DCA to PV of reco. prim. #bar{#Lambda}; DCA (cm); p_{T} (GeV/c)",
				80, 0., 4., nPtBinsV0Xi, PtBinsV0Xi);

		fHistRCSecAntiLambdaDCAPt = new TH2F("fHistRCSecAntiLambdaDCAPt",
				"DCA to PV of reco. sec. #bar{#Lambda}; DCA (cm); p_{T} (GeV/c)",
				80, 0., 4., nPtBinsV0Xi, PtBinsV0Xi);
		//----------------------------------------------------------- XI -------------------------------------------------------
		fHistRCXiMinusPt = new TH2F("fHistRCXiMinusPt", "pt of reconstructed #Xi^{-}; p_{T} of #Xi [GeV/c]]; #eta",nPtBinsHXi,PtBinsHXi, nEtaBins, EtaBins);

		fHistMCtruthXiMinusPt = new TH2F("fHistMCtruthXiMinusPt", "pt of generated #Xi^{-}; p_{T} of #Xi [GeV/c]; #eta",nPtBinsHXi,PtBinsHXi, nEtaBins, EtaBins); 
		//------------------------------ PT XI-------------------
		fHistRCXiPt = new TH1F("fHistRCXiPt", "pt of reconstructed #Xi^{-}; p_{T} of #Xi [GeV/c] ",nPtBinsHXi,PtBinsHXi); //mc xipt   
		TH1F*HistRCXiPt = new TH1F("HistRCXiPt", "Ppt of reconstructed #Xi^{-}; p_{T} of #Xi [GeV/c] ",nPtBinsHXi,PtBinsHXi);//mc lpt    


		fHistMCXiMinusPt = new TH1F("fHistMCXiMinusPt", "pt of generated #Xi^{-}; p_{T} of #Xi [GeV/c]",nPtBinsHXi,PtBinsHXi); 

		for(Int_t i = 0; i < kNVtxZ; i++){
			fHistRCTrk[i] = new TH3F(Form("fHistRCTrkVtx%d", i),
					"pt, eta, phi of reconstructed tracks; p_{T} [GeV/c]; #eta; #phi",
					30,0,15, nEtaBins, -fEtaCut, fEtaCut, nPhiBins, 0., 2.*TMath::Pi());

			fHistRCPriTrk[i] = new TH3F(Form("fHistRCPriTrkVtx%d", i),
					"pt, eta, phi of reconstructed prim. tracks; p_{T} [GeV/c]; #eta; #phi",
					30,0,15, nEtaBins, -fEtaCut, fEtaCut, nPhiBins, 0., 2.*TMath::Pi());

			fHistRCSecTrk[i] = new TH3F(Form("fHistRCSecTrkVtx%d", i),
					"pt, eta, phi of reconstructed sec. tracks; p_{T} [GeV/c]; #eta; #phi",
					30,0,15, nEtaBins, -fEtaCut, fEtaCut, nPhiBins, 0., 2.*TMath::Pi());

			fHistMCtruthTrk[i] = new TH3F(Form("fHistMCtruthTrkVtx%d", i),
					"pt, eta, phi of generated tracks; p_{T} [GeV/c]; #eta; #phi",
					30,0,15, nEtaBins, -fEtaCut, fEtaCut, 
					nPhiBins, 0., 2.*TMath::Pi());

			fHistRCK0s[i] = new TH3F(Form("fHistRCK0sVtx%d", i), 
					"pt, eta, phi of reconstructed K0s; p_{T} [GeV/c]; #eta; #phi",
					nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);

			fHistMCtruthK0s[i] = new TH3F(Form("fHistMCtruthK0sVtx%d", i), 
					"pt, eta and phi of generated K0s; p_{T} [GeV/c]; #eta; #phi",
					nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);


			fHistRCLambda[i] 
				= new TH3F(Form("fHistRCLambdaVtx%d", i), 
						"pt, eta and phi of reconstructed Lambda; p_{T} [GeV/c]; #eta; #phi",
						nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);

			fHistRCAntiLambda[i] 
				= new TH3F(Form("fHistRCAntiLambdaVtx%d", i), 
						"pt, eta and phi of reconstructed #bar{#Lambda}; p_{T} [GeV/c]; #eta; #phi",
						nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);


			fHistMCtruthLambda[i] 
				= new TH3F(Form("fHistMCtruthLambdaVtx%d", i), 
						"pt, eta and phi of generated #Lambda; p_{T} [GeV/c]; #eta; #phi",
						nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);

			fHistMCtruthAntiLambda[i] 
				= new TH3F(Form("fHistMCtruthAntiLambdaVtx%d", i),
						"pt, eta and phi of generated #bar{#Lambda}; p_{T} [GeV/c]; #eta; #phi",
						nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);
			//-----------------------------------------------  XI ----------------------------------------------------------------------

			fHistRCXiMinus[i] = new TH3F(Form("fHistRCXiMinusVtx%d", i), 
					"pt, eta , phi of reconstructed #Xi^{-}; p_{T} of #Xi [GeV/c]; #eta; #phi",
					nPtBinsHXi,PtBinsHXi,nEtaBins, EtaBins, nPhiBins, PhiBins);
			

			fHistMCtruthXiMinus[i] 
				= new TH3F(Form("fHistMCtruthXiMinusVtx%d", i), 
						"pt, eta and phi of generated #Xi^{-}; p_{T} of #Xi [GeV/c]; #eta; #phi",
						nPtBinsHXi,PtBinsHXi, nEtaBins, EtaBins,nPhiBins, PhiBins);
		}
	}

	fOutput->Add(fHistArmenterosPodolanski);

	if(fAnalysisMC){
		fOutput->Add(fHistRCTrkPtEta);
		fOutput->Add(fHistRCK0sPt);
		fOutput->Add(fHistRCLambdaPt);
		fOutput->Add(fHistallRCLambdaPt);
		fOutput->Add(fHistallDCALambdaPt);
		fOutput->Add(fHistallDLLambdaPt);
		fOutput->Add(fHistallRLambdaPt);
		fOutput->Add(fHistRCAntiLambdaPt);
		fOutput->Add(fHistRCPriTrkPtEta);
		fOutput->Add(fHistRCSecTrkPtEta);
		fOutput->Add(fHistMCtruthTrkPtEta);
		fOutput->Add(fHistMCtruthK0sPt);
		fOutput->Add(fHistMCtruthLambdaPt);
		fOutput->Add(fHistMCtruthAntiLambdaPt);
		fOutput->Add(fHistRCPrimLambdaDCAPt);
		fOutput->Add(fHistRCPrimLambdaDLPt);
		fOutput->Add(fHistRCPrimLambdaCTPt);
		fOutput->Add(fHistRCPrimLambdaRPt);
		fOutput->Add(fHistRCXiMinusPt);
		fOutput->Add(fHistMCtruthXiMinusPt);
		fOutput->Add(fHistRCXiPt);
		fOutput->Add(HistRCXiPt);
		fOutput->Add(fHistMCXiMinusPt);
		fOutput->Add(fHistRCSecLambdaDCAPt);
		fOutput->Add(fHistRCSecLambdaDLPt);
		fOutput->Add(fHistRCSecLambdaCTPt);
		fOutput->Add(fHistRCSecLambdaRPt);



		fOutput->Add(fHistRCPrimAntiLambdaDCAPt);
		fOutput->Add(fHistRCSecAntiLambdaDCAPt);

		for(Int_t i = 0 ; i < kNVtxZ; i++){

			fOutput->Add(fHistRCTrk[i]);
			fOutput->Add(fHistRCK0s[i]);
			fOutput->Add(fHistRCLambda[i]);
			fOutput->Add(fHistRCAntiLambda[i]);
			fOutput->Add(fHistRCPriTrk[i]);
			fOutput->Add(fHistRCSecTrk[i]);
			fOutput->Add(fHistMCtruthTrk[i]);
			fOutput->Add(fHistMCtruthK0s[i]);
			fOutput->Add(fHistMCtruthLambda[i]);
			fOutput->Add(fHistMCtruthAntiLambda[i]);
			fOutput->Add(fHistRCXiMinus[i]);
			fOutput->Add(fHistMCtruthXiMinus[i]);
		}

	}

	PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram


	// Settings for event mixing -------------------------------------
	if(fEventMixing){
	Int_t trackDepth = fMixingTracks;
	Int_t poolSize   = fPoolSize;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager

	fPoolMgr = new AliEventPoolManager(poolSize, trackDepth, nCentralityBins, centBins, nZvtxBins, vertexBins);
    }
//----------------------------------------------
}

void AliAnalysisTaskLambdaK0s::AddQAEvent()
{
	TList *tQAEvent;
	tQAEvent = new TList();
	tQAEvent->SetOwner();
	tQAEvent->SetName("EventInput");

	  
	TH1F *fhEventBf = new TH1F("fhEventBf", "Event Number; Counts; Number of Events",1, 0.,1);
	tQAEvent->Add(fhEventBf);

	TH1F *fhEventAllSel = new TH1F("fhEventAllSel", "Event Number; selections; Number of Events", 4, 0, 4);
	tQAEvent->Add(fhEventAllSel);

	TH1F *fHistVtxAf = new TH1F("fHistVtxAf", ";Z vertex", 8, -8., 8.);     
	tQAEvent->Add(fHistVtxAf);

	TH1F * fHistMultiMain = new TH1F("fHistMultiMain", "Multiplicity of main events", 2000, 0, 2000);
	tQAEvent->Add(fHistMultiMain);
	TH1F *fHistCentAf = new TH1F("fHistCentAf", "Centrality; centrality (%)", 100, 0., 100.);
	tQAEvent->Add(fHistCentAf);

	TH2F *fPrimayVtxGlobalvsSPD = new TH2F("fPrimayVtxGlobalvsSPD",";Z_{vtx,tr} (cm);Z_{SPD,tr} (cm)",200,-20,20,200,-20,20);
	tQAEvent->Add(fPrimayVtxGlobalvsSPD);

	fOutput->Add(tQAEvent);
}

//-------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskLambdaK0s::AddAnalysisTrk()
{
	TList *tAnalysisTrk;
	tAnalysisTrk = new TList();
	tAnalysisTrk->SetOwner();
	tAnalysisTrk->SetName("AnalysisTrk");

	// defining bins for centrality
	const Int_t nCentralityBins  = 1;//11for ppb

	Double_t centBins[] = {0.,  10.};
	const Int_t nZvtxBins  = 6;
	Double_t vertexBins[] = {-8., -4., -2., 0., 2., 4., 8.};

	const Int_t nEtaBins = 36;         
	Double_t EtaBins[nEtaBins+1] = {0.};
	EtaBins[0] = -fEtaCut;
	for (Int_t i=0; i<nEtaBins; i++) {
		EtaBins[i+1] = EtaBins[i] + 2.*fEtaCut/nEtaBins; }

	const Int_t nPhiBins = 36;
	Double_t PhiBins[nPhiBins+1] = {0.};
	PhiBins[0] = 0;
	for (Int_t i=0; i<nPhiBins; i++)
	{ PhiBins[i+1] = PhiBins[i] + 2.*TMath::Pi()/nPhiBins; }


	const Int_t nPtTrBins = 1;
	const Double_t PtTrBins[2] = {fPtTrigMin, fPtTrigMax}; 


	TH3F *fHistTrk = new TH3F("fHistTrk", 
			"pt, eta and phi of reconstructed tracks; p_{T} [GeV/c]; #eta; #phi", 
			30, 0, 15, nEtaBins, -fEtaCut, fEtaCut, nPhiBins, 0., 2.*TMath::Pi());
	tAnalysisTrk->Add(fHistTrk);

	TH2F *fHistTrkVzvsEta = new TH2F("fHistTrkVzvsEta",
			"Trk Vertex Z and eta;#eta;Vz [cm]",
			nEtaBins, -fEtaCut, fEtaCut, nZvtxBins,vertexBins);  
	tAnalysisTrk->Add(fHistTrkVzvsEta);

	TH2F *fHistTrigPtAll = new TH2F("fHistTrigPtAll", 
			"pt distribution for all triggered particles; p_{T} [GeV/c]; Z_{vtx} [cm]", 
			nPtTrBins, PtTrBins, nZvtxBins , vertexBins);
	tAnalysisTrk->Add(fHistTrigPtAll);

	TH2F *fTriggerEtaPhi = new TH2F("fTriggerEtaPhi","Trigger particle;#varphi (rad);#eta",nPhiBins,0.,2.*TMath::Pi(),100,-1.,1.);
	tAnalysisTrk->Add(fTriggerEtaPhi);

	TH1F *fTriggerEta = new TH1F("fTriggerEta","Trigger particles VS Eta;#eta",nEtaBins,-fEtaCut,fEtaCut);
	tAnalysisTrk->Add(fTriggerEta);

	TH1F *fTriggerPhi = new TH1F("fTriggerPhi","Trigger particles VS Phi;#phi",nPhiBins, 0., 2.*TMath::Pi());
	tAnalysisTrk->Add(fTriggerPhi);

	fOutput->Add(tAnalysisTrk);
} 
//________________________________________________________________________
void AliAnalysisTaskLambdaK0s::AddAnalysisK0s()
{
	TList *tAnalysisK0s;
	tAnalysisK0s = new TList();
	tAnalysisK0s->SetOwner();
	tAnalysisK0s->SetName("AnalysisK0s");

	// defining bins for centrality
	const Int_t nCentralityBins  = 1;//11 for ppb
	Double_t centBins[] = {0.,10.};

	const Int_t nZvtxBins  = 6;
	Double_t vertexBins[] = {-8., -4., -2., 0., 2., 4., 8.};

	
	const Int_t nPtBinsV0Xi = 11;
	const Double_t PtBinsV0Xi[12] = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5,4.0, 5.0, 6.0, 7.0, 8.0,10.0};
	const Int_t nPtAssoBins =nPtBinsV0Xi ;

	const Int_t nPtTrBins = 1;
	const Double_t PtTrBins[2] = {fPtTrigMin, fPtTrigMax}; 

	// defining bins for mass distributions
	Int_t nBins = 160;
	Double_t mMassks0[]={0.40, 0.58};
	Double_t mMassMin[] ={0.40, 1.04};
	Double_t mMassMax[] ={0.58, 1.21};

	Double_t massK[161]={0.};
	massK[0]=0.4;
	for(Int_t j=0;j<160;j++)
	{
		massK[j+1]= massK[j] +0.001125;
	}


	// defining bins for dPhi distributions
	const Int_t ndPhiBins = 32;
	const Double_t kPi = TMath::Pi();
	Double_t PhiMin = -kPi/2.;
	Double_t PhiMax = -kPi/2. + 2.*kPi;

	// defining bins for dEta distributions
	const Int_t ndEtaBins = 32;
	Double_t EtaMin = -2.*fEtaCut;
	Double_t EtaMax = 2.*fEtaCut;


	const Int_t nEtaBins = 18;          
	Double_t EtaBins[nEtaBins+1] = {0.};
	EtaBins[0] = -fEtaCut;
	for (Int_t i=0; i<nEtaBins; i++) {
		EtaBins[i+1] = EtaBins[i] + 2.*fEtaCut/nEtaBins; }

	const Int_t nPhiBins = 36;
	Double_t PhiBins[nPhiBins+1] = {0.};
	PhiBins[0] = 0;
	for (Int_t i=0; i<nPhiBins; i++)
	{ PhiBins[i+1] = PhiBins[i] + 2.*TMath::Pi()/nPhiBins; }

	TH3F *fHistK0s = new TH3F("fHistK0s",
			"pt, eta and phi of reconstructed K0s; p_{T} [GeV/c]; #eta; #phi",
			nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);
	tAnalysisK0s->Add(fHistK0s);

	TH1F *fHistK0sEta = new TH1F("fHistK0sEta","Eta distribution of K0s;#eta",nEtaBins,EtaBins);
	tAnalysisK0s->Add(fHistK0sEta);

	TH1F *fHistK0sPhi = new TH1F("fHistK0sPhi","Phi distribution of K0s;#phi",nPhiBins,PhiBins);
	tAnalysisK0s->Add(fHistK0sPhi);

	TH2F *fHistK0sVzvsEta = new TH2F("fHistK0sVzvsEta",
			"K0s Vertex Z and eta;#eta;Vz [cm]",
			nEtaBins, -fEtaCut, fEtaCut,50,-8.,8.);
	tAnalysisK0s->Add(fHistK0sVzvsEta);
	TH2F *fHistMassK0s = new TH2F("fHistMassK0s","invassmass k0s;mass (GeV/c^{2}); p_{T} (GeV/c)", nBins,massK,nPtBinsV0Xi,PtBinsV0Xi);

	tAnalysisK0s->Add(fHistMassK0s);
	fHistMassK0s->Sumw2();  


	//...........                  //0-dPhi, 1-dEta, 2-Pt trig, 3-Pt V0, 4-inv. mass,   5-Centrality, 6-Zvertex,
	const Int_t   corBinsK0s[7]    = {ndPhiBins,ndEtaBins,nPtTrBins,nPtAssoBins, nBins,        nCentralityBins, nZvtxBins};
	const Double_t corMinK0s[7]    = {PhiMin,   EtaMin,PtTrBins[0], fPtAssoMin,  mMassMin[0],  centBins[0],        vertexBins[0]};
	const Double_t corMaxK0s[7]    = {PhiMax,   EtaMax,PtTrBins[nPtTrBins],fPtAssoMax,mMassMax[0] , centBins[nCentralityBins],vertexBins[nZvtxBins]};

	const Int_t   corBinsK0sMix[7] = {ndPhiBins,ndEtaBins,nPtTrBins,nPtAssoBins, nBins, nCentralityBins, nZvtxBins };
	const Double_t corMinK0sMix[7] = {PhiMin,   EtaMin,PtTrBins[0], fPtAssoMin,  mMassMin[0],  centBins[0], vertexBins[0] };
	const Double_t corMaxK0sMix[7] = {PhiMax,   EtaMax,PtTrBins[nPtTrBins],fPtAssoMax,mMassMax[0] , centBins[nCentralityBins],vertexBins[nZvtxBins]}; 



	THnSparseF *fHistdPhidEtaSibK0s={NULL};
	THnSparseF *fHistdPhidEtaMixK0s={NULL};
	// Create correlation histograms
	fHistdPhidEtaSibK0s
		= new THnSparseF("fHistdPhidEtaSibK0s",
				";#Delta#phi; #Delta#eta;p_{T}^{Tri} (GeV/c); p_{T}^{As} (GeV/c); mass (GeV/c^{2})",
				7, corBinsK0s, corMinK0s, corMaxK0s);

	fHistdPhidEtaSibK0s->GetAxis(3)->Set(nPtAssoBins, PtBinsV0Xi);
	fHistdPhidEtaSibK0s->GetAxis(6)->Set(nZvtxBins,vertexBins);



	fHistdPhidEtaSibK0s->Sumw2();
	tAnalysisK0s->Add(fHistdPhidEtaSibK0s);


	if(fEventMixing){
		fHistdPhidEtaMixK0s
			= new THnSparseF("fHistdPhidEtaMixK0s",
					";#Delta#phi; #Delta#eta; p_{T}^{Tri} (GeV/c); p_{T}^{As} (GeV/c); mass (GeV/c^{2})",
					7, corBinsK0sMix, corMinK0sMix, corMaxK0sMix);
		fHistdPhidEtaMixK0s->GetAxis(3)->Set(nPtAssoBins, PtBinsV0Xi);
		fHistdPhidEtaMixK0s->GetAxis(6)->Set(nZvtxBins,vertexBins);
		tAnalysisK0s->Add(fHistdPhidEtaMixK0s);
		fHistdPhidEtaMixK0s->Sumw2();


	}
	fOutput->Add(tAnalysisK0s);

}

//________________________________________________________________________
void AliAnalysisTaskLambdaK0s::AddAnalysisLambda()
{
	TList *tAnalysisLambda;
	tAnalysisLambda = new TList();
	tAnalysisLambda->SetOwner();
	tAnalysisLambda->SetName("AnalysisLambda");

	// defining bins for centrality
	const Int_t nCentralityBins  = 1;//11
	Double_t centBins[] = {fMinCent , fMaxCent};

	const Int_t nZvtxBins  = 6;
	Double_t vertexBins[] = {-8., -4., -2., 0., 2., 4.,  8.};

	const Int_t nPtBinsV0Xi = 11;
	const Double_t PtBinsV0Xi[12] = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5,4.0, 5.0, 6.0, 7.0, 8.0,10.0};
	const Int_t nPtAssoBins =nPtBinsV0Xi ;


	// pt bins of associated particles (cascades) for the analysis
	const Int_t nPtBinsHXi = 27;
	const Double_t PtBinsHXi[28] = {0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5, 5.5, 6, 7, 8, 9, 10}; 


	const Int_t nPtTrBins = 1;
	const Double_t PtTrBins[2] = {fPtTrigMin, fPtTrigMax}; 

	// defining bins for mass distributions
	Int_t nBins = 160;
	Double_t mMasslambda[]={1.04,1.21 };
	Double_t mMassMin[] ={0.40, 1.04};
	Double_t mMassMax[] ={0.58, 1.21};

	Double_t massL[161]={0.};
	massL[0]=1.04;
	for(Int_t d=0;d<160;d++)
	{
		massL[d+1]= massL[d] +0.0010625;
	}

	// defining bins for dPhi distributions
	const Int_t ndPhiBins = 32;
	const Double_t kPi = TMath::Pi();
	Double_t PhiMin = -kPi/2.;
	Double_t PhiMax = -kPi/2. + 2.*kPi;

	// defining bins for dEta distributions
	const Int_t ndEtaBins = 32;
	Double_t EtaMin = -2.*fEtaCut;
	Double_t EtaMax = 2.*fEtaCut;


	const Int_t nEtaBins = 18;          
	Double_t EtaBins[nEtaBins+1] = {0.};
	EtaBins[0] = -fEtaCut;
	for (Int_t i=0; i<nEtaBins; i++) {
		EtaBins[i+1] = EtaBins[i] + 2.*fEtaCut/nEtaBins; }

	const Int_t nPhiBins = 36;
	Double_t PhiBins[nPhiBins+1] = {0.};
	PhiBins[0] = 0;
	for (Int_t i=0; i<nPhiBins; i++)
	{ PhiBins[i+1] = PhiBins[i] + 2.*TMath::Pi()/nPhiBins; }

	const Int_t nDCA = 80;

	Double_t DCAL[81]={0.};
	DCAL[0]=0.;
	for(Int_t i=0;i<80;i++)
	{
		DCAL[i+1]= DCAL[i] +0.0125;
	}


	TH2F *fHistDCALambda = new TH2F("fHistDCALambda","dca Lambda ;DCA(cm); p_{T} (GeV/c)",nPtBinsV0Xi,PtBinsV0Xi,nDCA,DCAL);
	TH2F *fHistMassLambda = new TH2F("fHistMassLambda","invassmass Lambda ;mass (GeV/c^{2}); p_{T} (GeV/c)",nBins,massL,nPtBinsV0Xi,PtBinsV0Xi);
	TH3F *fHistMassDCALambda = new TH3F("fHistMassDCALambda","invassmass DCA Lambda ;mass (GeV/c^{2}); p_{T} (GeV/c);DCA(cm)",nBins,massL,nPtBinsV0Xi,PtBinsV0Xi,nDCA,DCAL); 
	tAnalysisLambda->Add(fHistMassLambda);
	tAnalysisLambda->Add(fHistMassDCALambda);      
	tAnalysisLambda->Add(fHistDCALambda);
	fHistMassLambda->Sumw2();    
	fHistMassDCALambda->Sumw2();
	fHistDCALambda->Sumw2();
	TH2F *fHistXiPtVSLambdaPt = new TH2F("fHistXiPtVSLambdaPt", "Xi Pt VS Lambda Pt;p_{T} of #Xi [GeV/c];p_{T} of #Lambda from #Xi [GeV/c]",nPtBinsHXi,PtBinsHXi,nPtBinsV0Xi,PtBinsV0Xi);
	tAnalysisLambda->Add(fHistXiPtVSLambdaPt);

	TH3F *fHistLambda = new TH3F("fHistLambda",
			"pt, eta and phi of reconstructed Lambda; p_{T} [GeV/c]; #eta; #phi",
			nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);
	tAnalysisLambda->Add(fHistLambda);

	TH1F *fHistLambdaEta = new TH1F("fHistLambdaEta","Eta distribution of Lambda;#eta",nEtaBins,EtaBins);
	tAnalysisLambda->Add(fHistLambdaEta);

	TH1F *fHistLambdaPhi = new TH1F("fHistLambdaPhi","Phi distribution of Lambda;#Phi",nPhiBins,PhiBins);
	tAnalysisLambda->Add(fHistLambdaPhi);

	TH2F *fHistLambdaVzvsEta = new TH2F("fHistLambdaVzvsEta",
			"Lambda Vertex Z and eta;#eta;Vz [cm]",
			nEtaBins, -fEtaCut, fEtaCut,50,-8.,8.);
	tAnalysisLambda->Add(fHistLambdaVzvsEta);

	TH3F *fHistLambdaDphiDCAPtSig = new TH3F("fHistLambdaDphiDCAPtSig",
			";#Delta#phi; p_{T} [GeV/c]; DCA [cm]",
			ndPhiBins, PhiMin, PhiMax,
			nPtAssoBins, fPtAssoMin, fPtAssoMax,
			80, 0., 1.);//4
	TH3F *fHistLambdaDphiDLPtSig = new TH3F("fHistLambdaDphiDLPtSig",";#Delta#phi; p_{T} [GeV/c]; DL [cm]",ndPhiBins, PhiMin,                                                 PhiMax, nPtAssoBins, fPtAssoMin, fPtAssoMax,80,0.,23.);

	TH3F *fHistLambdaDphiCTPtSig = new TH3F("fHistLambdaDphiCTPtSig",";#Delta#phi; p_{T} [GeV/c];CT [cm]  ",ndPhiBins, PhiMin,                                                 PhiMax, nPtAssoBins, fPtAssoMin, fPtAssoMax,80,0.,23.);


	TH3F *fHistLambdaDphiRPtSig = new TH3F("fHistLambdaDphiRPtSig",";#Delta#phi; p_{T} [GeV/c]; R [cm]",ndPhiBins, PhiMin,                                                 PhiMax, nPtAssoBins, fPtAssoMin, fPtAssoMax,80,5.,100.);
	tAnalysisLambda->Add(fHistLambdaDphiDCAPtSig);
	tAnalysisLambda->Add(fHistLambdaDphiDLPtSig);
	tAnalysisLambda->Add(fHistLambdaDphiCTPtSig);
	tAnalysisLambda->Add(fHistLambdaDphiRPtSig);

	TH3F *fHistLambdaDphiDCAPtBkg = new TH3F("fHistLambdaDphiDCAPtBkg", ";#Delta#phi; p_{T} [GeV/c]; DCA [cm]",ndPhiBins, PhiMin,    PhiMax, nPtAssoBins, fPtAssoMin, fPtAssoMax, 80, 0., 1.);
	TH3F *fHistLambdaDphiDLPtBkg = new TH3F("fHistLambdaDphiDLPtBkg", ";#Delta#phi; p_{T} [GeV/c]; DL [cm]",ndPhiBins, PhiMin,    PhiMax, nPtAssoBins, fPtAssoMin, fPtAssoMax, 80, 0., 23.);
	TH3F *fHistLambdaDphiCTPtBkg = new TH3F("fHistLambdaDphiCTPtBkg", ";#Delta#phi; p_{T} [GeV/c];CT [cm] ",ndPhiBins, PhiMin,    PhiMax, nPtAssoBins, fPtAssoMin, fPtAssoMax, 80, 0., 23.);
	TH3F *fHistLambdaDphiRPtBkg = new TH3F("fHistLambdaDphiRPtBkg", ";#Delta#phi; p_{T} [GeV/c]; R [cm]",ndPhiBins, PhiMin,    PhiMax, nPtAssoBins, fPtAssoMin, fPtAssoMax, 80, 5., 100.);

	tAnalysisLambda->Add(fHistLambdaDphiDCAPtBkg);
	tAnalysisLambda->Add(fHistLambdaDphiDLPtBkg);
	tAnalysisLambda->Add(fHistLambdaDphiCTPtBkg);
	tAnalysisLambda->Add(fHistLambdaDphiRPtBkg);

	//...........                  //0-dPhi, 1-dEta, 2-Pt trig, 3-Pt V0, 4-inv. mass,   5-Centrality, 6-Zvertex,
	const Int_t   corBinsLambda[7]    = {ndPhiBins,ndEtaBins,nPtTrBins,nPtAssoBins, nBins,        nCentralityBins, nZvtxBins};
	const Double_t corMinLambda[7]    = {PhiMin,   EtaMin,PtTrBins[0], fPtAssoMin,  mMassMin[1],  centBins[0],        vertexBins[0]};
	const Double_t corMaxLambda[7]    = {PhiMax,   EtaMax,PtTrBins[nPtTrBins],fPtAssoMax,mMassMax[1] , centBins[nCentralityBins],vertexBins[nZvtxBins]};

	const Int_t   corBinsLambdaMix[7] = {ndPhiBins,ndEtaBins,nPtTrBins,nPtAssoBins, nBins,        nCentralityBins, nZvtxBins };
	const Double_t corMinLambdaMix[7] = {PhiMin,   EtaMin,PtTrBins[0], fPtAssoMin,  mMassMin[1],  centBins[0],        vertexBins[0] };
	const Double_t corMaxLambdaMix[7] = {PhiMax,   EtaMax,PtTrBins[nPtTrBins],fPtAssoMax,mMassMax[1] , centBins[nCentralityBins],vertexBins[nZvtxBins]};


	THnSparseF *fHistdPhidEtaSibLambda={NULL};//liu 5
	THnSparseF *fHistdPhidEtaMixLambda={NULL};//5
	// Create correlation histograms
	fHistdPhidEtaSibLambda
		= new THnSparseF("fHistdPhidEtaSibLambda",
				";#Delta#phi; #Delta#eta;p_{T}^{Tri} (GeV/c);  p_{T}^{As} (GeV/c); mass (GeV/c^{2})",
				7, corBinsLambda, corMinLambda, corMaxLambda);


	fHistdPhidEtaSibLambda->GetAxis(3)->Set(nPtAssoBins, PtBinsV0Xi);
	fHistdPhidEtaSibLambda->GetAxis(6)->Set(nZvtxBins,vertexBins);




	fHistdPhidEtaSibLambda->Sumw2();
	tAnalysisLambda->Add(fHistdPhidEtaSibLambda);

	if(fEventMixing){
		fHistdPhidEtaMixLambda
			= new THnSparseF("fHistdPhidEtaMixLambda",
					";#Delta#phi; #Delta#eta; p_{T}^{Tri} (GeV/c); p_{T}^{As} (GeV/c); mass (GeV/c^{2})",
					7, corBinsLambdaMix, corMinLambdaMix, corMaxLambdaMix);
		fHistdPhidEtaMixLambda->GetAxis(3)->Set(nPtAssoBins, PtBinsV0Xi);
		fHistdPhidEtaMixLambda->GetAxis(6)->Set(nZvtxBins,vertexBins);
		fHistdPhidEtaMixLambda->Sumw2();
		tAnalysisLambda->Add(fHistdPhidEtaMixLambda);
	}

	fOutput->Add(tAnalysisLambda);
}
//------------------------------------------------------------
void AliAnalysisTaskLambdaK0s::AddAnalysisAntiLambda()
{
	TList *tAnalysisAntiLambda;
	tAnalysisAntiLambda = new TList();
	tAnalysisAntiLambda->SetOwner();
	tAnalysisAntiLambda->SetName("AnalysisAntiLambda");

	// defining bins for centrality
	const Int_t nCentralityBins  = 1;//11
	Double_t centBins[] = {fMinCent , fMaxCent};


	const Int_t nZvtxBins  = 6;
	Double_t vertexBins[] = {-8.,-4., -2., 0., 2., 4., 8.};


	const Int_t nPtBinsV0Xi = 11;
	const Double_t PtBinsV0Xi[12] = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5,4.0, 5.0, 6.0, 7.0, 8.0,10.0};
	const Int_t nPtAssoBins =nPtBinsV0Xi ;


	// pt bins of associated particles (cascades) for the analysis
	const Int_t nPtBinsHXi = 27;
	const Double_t PtBinsHXi[28] = {0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5, 5.5, 6, 7, 8, 9, 10}; 


	const Int_t nPtTrBins = 1;
	const Double_t PtTrBins[2] = {fPtTrigMin, fPtTrigMax}; 

	// defining bins for mass distributions
	Int_t nBins = 160;
	Double_t mMassMin[] ={0.40, 1.04};
	Double_t mMassMax[] ={0.58, 1.21};

	// defining bins for dPhi distributions
	const Int_t ndPhiBins = 32;
	const Double_t kPi = TMath::Pi();
	Double_t PhiMin = -kPi/2.;
	Double_t PhiMax = -kPi/2. + 2.*kPi;

	// defining bins for dEta distributions
	const Int_t ndEtaBins = 32;
	Double_t EtaMin = -2.*fEtaCut;
	Double_t EtaMax = 2.*fEtaCut;


	const Int_t nEtaBins = 18;
	Double_t EtaBins[nEtaBins+1] = {0.};
	EtaBins[0] = -fEtaCut;
	for (Int_t i=0; i<nEtaBins; i++) {
		EtaBins[i+1] = EtaBins[i] + 2.*fEtaCut/nEtaBins; }

	const Int_t nPhiBins = 36;
	Double_t PhiBins[nPhiBins+1] = {0.};
	PhiBins[0] = 0;

	for (Int_t i=0; i<nPhiBins; i++)
	{ PhiBins[i+1] = PhiBins[i] + 2.*TMath::Pi()/nPhiBins; }

	// Xi-AntiLambda Histogram 
	TH2F *fHistXiPtVSAntiLambdaPt = new TH2F("fHistXiPtVSAntiLambdaPt","Xi Pt VS AntiLambda Pt;p_{T} of #Xi [GeV/c];p_{T} of #bar{#Lambda} from #Xi [GeV/c]",nPtBinsHXi,PtBinsHXi,nPtBinsV0Xi,PtBinsV0Xi);
	tAnalysisAntiLambda->Add(fHistXiPtVSAntiLambdaPt);

	//-------------
	TH3F *fHistAntiLambda = new TH3F("fHistAntiLambda",
			"pt, eta and phi of reconstructed AntiLambda; p_{T} [GeV/c]; #eta; #phi",
			nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);
	tAnalysisAntiLambda->Add(fHistAntiLambda);

	TH1F *fHistAntiLambdaEta = new TH1F("fHistAntiLambdaEta","Eta distribution of AntiLambda;#eta",nEtaBins,EtaBins);
	tAnalysisAntiLambda->Add(fHistAntiLambdaEta);

	TH1F *fHistAntiLambdaPhi = new TH1F("fHistAntiLambdaPhi","Phi distribution of AntiLambfa;#phi",nPhiBins,PhiBins);
	tAnalysisAntiLambda->Add(fHistAntiLambdaPhi);

	TH2F *fHistAntiLambdaVzvsEta = new TH2F("fHistAntiLambdaVzvsEta",
			"AntiLambda Vertex Z and eta;#eta;Vz [cm]",
			nEtaBins, -fEtaCut, fEtaCut,50,-8.,8.);
	tAnalysisAntiLambda->Add(fHistAntiLambdaVzvsEta);

	TH3F *fHistAntiLambdaDphiDCAPtSig = new TH3F("fHistAntiLambdaDphiDCAPtSig",
			";#Delta#phi; p_{T} [GeV/c]; DCA [cm]",
			ndPhiBins, PhiMin, PhiMax,
			nPtAssoBins, fPtAssoMin, fPtAssoMax,
			80, 0., 4.);
	tAnalysisAntiLambda->Add(fHistAntiLambdaDphiDCAPtSig);

	TH3F *fHistAntiLambdaDphiDCAPtBkg = new TH3F("fHistAntiLambdaDphiDCAPtBkg",
			";#Delta#phi; p_{T} [GeV/c]; DCA [cm]",
			ndPhiBins, PhiMin, PhiMax,
			nPtAssoBins, fPtAssoMin, fPtAssoMax,
			80, 0., 4.);
	tAnalysisAntiLambda->Add(fHistAntiLambdaDphiDCAPtBkg);


	//...........                  //0-dPhi, 1-dEta, 2-Pt trig, 3-Pt V0, 4-inv. mass,   5-Centrality, 6-Zvertex,
	const Int_t   corBinsAntiLambda[7]    = {ndPhiBins,ndEtaBins,nPtTrBins,nPtAssoBins, nBins,        nCentralityBins, nZvtxBins};
	const Double_t corMinAntiLambda[7]    = {PhiMin,   EtaMin,PtTrBins[0], fPtAssoMin,  mMassMin[1],  centBins[0],        vertexBins[0]};
	const Double_t corMaxAntiLambda[7]    = {PhiMax,   EtaMax,PtTrBins[nPtTrBins],fPtAssoMax,mMassMax[1] , centBins[nCentralityBins],vertexBins[nZvtxBins]};

	const Int_t   corBinsAntiLambdaMix[7] = {ndPhiBins,ndEtaBins,nPtTrBins,nPtAssoBins, nBins,        nCentralityBins, nZvtxBins };
	const Double_t corMinAntiLambdaMix[7] = {PhiMin,   EtaMin,PtTrBins[0], fPtAssoMin,  mMassMin[1],  centBins[0],        vertexBins[0] };
	const Double_t corMaxAntiLambdaMix[7] = {PhiMax,   EtaMax,PtTrBins[nPtTrBins],fPtAssoMax,mMassMax[1] , centBins[nCentralityBins],vertexBins[nZvtxBins]};

	THnSparseF *fHistdPhidEtaSibAntiLambda={NULL};
	THnSparseF *fHistdPhidEtaMixAntiLambda={NULL};
	// Create correlation histograms
	fHistdPhidEtaSibAntiLambda
		= new THnSparseF("fHistdPhidEtaSibAntiLambda",
				";#Delta#phi; #Delta#eta;p_{T}^{Tri} (GeV/c);  p_{T}^{As} (GeV/c); mass (GeV/c^{2})",
				7, corBinsAntiLambda, corMinAntiLambda, corMaxAntiLambda);

	fHistdPhidEtaSibAntiLambda->GetAxis(3)->Set(nPtAssoBins, PtBinsV0Xi);
	fHistdPhidEtaSibAntiLambda->GetAxis(6)->Set(nZvtxBins,vertexBins);


	fHistdPhidEtaSibAntiLambda->Sumw2();
	tAnalysisAntiLambda->Add(fHistdPhidEtaSibAntiLambda);

	if(fEventMixing){
		fHistdPhidEtaMixAntiLambda
			= new THnSparseF("fHistdPhidEtaMixAntiLambda",
					";#Delta#phi; #Delta#eta;p_{T}^{Tri} (GeV/c);  p_{T}^{As} (GeV/c); mass (GeV/c^{2})",
					7, corBinsAntiLambdaMix, corMinAntiLambdaMix, corMaxAntiLambdaMix);                                                                       
		fHistdPhidEtaMixAntiLambda->GetAxis(3)->Set(nPtAssoBins, PtBinsV0Xi);
		fHistdPhidEtaMixAntiLambda->GetAxis(6)->Set(nZvtxBins,vertexBins);
		fHistdPhidEtaMixAntiLambda->Sumw2();
		tAnalysisAntiLambda->Add(fHistdPhidEtaMixAntiLambda);
	}

	fOutput->Add(tAnalysisAntiLambda);
}

//========================================================== XiMinus ========================================================
void AliAnalysisTaskLambdaK0s::AddAnalysisXiMinus()
{

	TList *tAnalysisXiMinus;
	tAnalysisXiMinus = new TList();

	tAnalysisXiMinus->SetOwner();
	tAnalysisXiMinus->SetName("AnalysisXiMinus");

	// defining bins for centrality
	const Int_t nCentralityBins  = 1;
	Double_t centBins[] = {fMinCent , fMaxCent};
	const Double_t* centralityBins = centBins;

	// defining bins for Z vertex
	const Int_t nZvtxBins  = 6;
	Double_t vertexBins[] = {-8.,  -4., -2., 0., 2., 4., 8.};
	const Double_t* zvtxBins = vertexBins;

	// pt bins of associated particles (cascades) for the analysis
	const Int_t nPtBinsHXi = 27;
	const Double_t PtBinsHXi[28] = {0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5, 5.5, 6, 7, 8, 9, 10}; 


	// pt bins of for colrrlation
	const Int_t nPtAssoBinsXi= 5;
	const Double_t PtAssoBinsXi[6] = {1.,  2.,  3.,  4.,  5.,  7 }; 
	Double_t fPtAssoMinXi =1;
	Double_t fPtAssoMaxXi =7;


	const Int_t nPtTrBins = 1;
	const Double_t PtTrBins[2] = {fPtTrigMin, fPtTrigMax}; 

	// defining bins for mass distributions
	Int_t nBins = 160;
	Double_t mXiMin = 1.29;
	Double_t mXiMax = 1.37;
	Double_t masss[161]={0.};
	masss[0]=1.29;
	for(Int_t i=0;i<160;i++)
	{
		masss[i+1]= masss[i] +0.0005;
	}
	// defining bins for dPhi distributions
	const Int_t ndPhiBins = 32;
	const Double_t kPi = TMath::Pi();
	Double_t PhiMin = -kPi/2.;
	Double_t PhiMax = -kPi/2. + 2.*kPi;

	// defining bins for dEta distributions
	const Int_t ndEtaBins = 32;
	Double_t EtaMin = -1.6 ;
	Double_t EtaMax = 1.6 ;

	//defining bins of Eta distribution
	const Int_t nEtaBins = 16;         
	Double_t EtaBins[nEtaBins+1] = {0.};
	EtaBins[0] = -fEtaCut;
	for (Int_t i=0; i<nEtaBins; i++) {
		EtaBins[i+1] = EtaBins[i] + 2.*fEtaCut/nEtaBins; }
	//defining bins of Phi distribution
	const Int_t nPhiBins = 36;
	Double_t PhiBins[nPhiBins+1] = {0.};
	PhiBins[0] = 0;
	for (Int_t i=0; i<nPhiBins; i++)
	{ PhiBins[i+1] = PhiBins[i] + 2.*TMath::Pi()/nPhiBins; }

	//;;;;;;;;;;;;;;;;;;;;;

	//...........                  //0-dPhi, 1-dEta, 2-Pt trig, 3-Pt V0, 4-inv. mass,   5-Centrality, 6-Zvertex,
	const Int_t   corBinsXiMinus[7]    = {ndPhiBins,ndEtaBins,nPtTrBins,nPtAssoBinsXi, nBins,        nCentralityBins, nZvtxBins};
	const Double_t corMinXiMinus[7]    = {PhiMin,   EtaMin,PtTrBins[0], fPtAssoMinXi,  mXiMin,  centBins[0],        vertexBins[0]};
	const Double_t corMaxXiMinus[7]    = {PhiMax,   EtaMax,PtTrBins[nPtTrBins],fPtAssoMaxXi,mXiMax , centBins[nCentralityBins],vertexBins[nZvtxBins]};

	const Int_t   corBinsXiMinusMix[7] = {ndPhiBins,ndEtaBins,nPtTrBins,nPtAssoBinsXi, nBins,        nCentralityBins, nZvtxBins };
	const Double_t corMinXiMinusMix[7] = {PhiMin,   EtaMin,PtTrBins[0], fPtAssoMinXi,  mXiMin,  centBins[0],        vertexBins[0] };
	const Double_t corMaxXiMinusMix[7] = {PhiMax,   EtaMax,PtTrBins[nPtTrBins],fPtAssoMaxXi,mXiMax , centBins[nCentralityBins],vertexBins[nZvtxBins]};

	THnSparseF *fHistdPhidEtaSibXiMinus={NULL};//5
	THnSparseF *fHistdPhidEtaMixXiMinus={NULL};//
	// Create correlation histograms
	// THnSparseF *fHistdPhidEtaSibXiMinus[i]
	fHistdPhidEtaSibXiMinus
		= new THnSparseF("fHistdPhidEtaSibXiMinus",
				";#Delta#phi; #Delta#eta;p_{T}^{Tri} (GeV/c);  p_{T}^{As} (GeV/c); mass (GeV/c^{2})",
				7, corBinsXiMinus, corMinXiMinus, corMaxXiMinus);

	fHistdPhidEtaSibXiMinus->GetAxis(3)->Set(nPtAssoBinsXi, PtAssoBinsXi);
	fHistdPhidEtaSibXiMinus->GetAxis(6)->Set(nZvtxBins,vertexBins);


	fHistdPhidEtaSibXiMinus->Sumw2();
	tAnalysisXiMinus->Add(fHistdPhidEtaSibXiMinus);

	if(fEventMixing){
		fHistdPhidEtaMixXiMinus
			= new THnSparseF("fHistdPhidEtaMixXiMinus",
					";#Delta#phi; #Delta#eta;p_{T}^{Tri} (GeV/c);  p_{T}^{As} (GeV/c); mass (GeV/c^{2})",
					7, corBinsXiMinusMix, corMinXiMinusMix, corMaxXiMinusMix);                                                                       
		fHistdPhidEtaMixXiMinus->GetAxis(3)->Set(nPtAssoBinsXi, PtAssoBinsXi);
		fHistdPhidEtaMixXiMinus->GetAxis(6)->Set(nZvtxBins,vertexBins);
		fHistdPhidEtaMixXiMinus->Sumw2();
		tAnalysisXiMinus->Add(fHistdPhidEtaMixXiMinus);
	}


	const Int_t spBins[4] = {nBins, nPtBinsHXi, nCentralityBins, nZvtxBins };
	const Double_t spMinXi[4] = {mXiMin, PtBinsHXi[0], centralityBins[0], vertexBins[0]};
	const Double_t spMaxXi[4] = {mXiMax, PtBinsHXi[nPtBinsHXi], centralityBins[nCentralityBins], vertexBins[nZvtxBins]};

	TH3F *fHistMassXi = new TH3F("fHistMassXi","mass (GeV/c^{2}); p_{T} (GeV/c);Vz [cm]",nBins,masss,nPtBinsHXi,PtBinsHXi , nZvtxBins, vertexBins);

	tAnalysisXiMinus->Add(fHistMassXi);
	fHistMassXi->Sumw2();



	TH3F *fHistXiMinus = new TH3F("fHistXiMinus","pt, eta and phi of reconstructed XiMinus;p_{T} of #Xi [GeV/c]; #eta; #phi",nPtBinsHXi,PtBinsHXi,nEtaBins, EtaBins, nPhiBins,PhiBins);                      
	tAnalysisXiMinus->Add(fHistXiMinus);

	TH2F *fHistXiMinusVzvsEta = new TH2F("fHistXiMinusVzvsEta", "XiMinus Vertex Z and eta;#eta;Vz [cm]", nEtaBins, -fEtaCut, fEtaCut,80,-8.,8.);
	tAnalysisXiMinus->Add(fHistXiMinusVzvsEta);

	TH1F *fHistXiMinusEta = new TH1F("fHistXiMinusEta","Eta distribution of XiMinus;#eta",nEtaBins,EtaBins);
	tAnalysisXiMinus->Add(fHistXiMinusEta);

	TH1F *fHistXiMinusPhi = new TH1F("fHistXiMinusPhi","Phi distribution of XiMinus;#phi",nPhiBins,PhiBins);
	tAnalysisXiMinus->Add(fHistXiMinusPhi);
	//-------------------------------------------------
	TH1F *fHistXiMinusPT = new TH1F("fHistXiMinusPT","PT distribution of XiMinus; p_{T} of #Xi [GeV/c]; Number of Counts",nPtBinsHXi,PtBinsHXi);
	fHistXiMinusPT->Sumw2(); 
	tAnalysisXiMinus->Add(fHistXiMinusPT);  

	TH2F *fHistXiMinusPTMASS = new TH2F("fHistXiMinusPTMASS", "XiMinus PT distribution and Mass;p_{T} [GeV];Mass (GeV/^{2})",nPtBinsHXi,PtBinsHXi,
			nBins,1.29,1.37);
	tAnalysisXiMinus->Add(fHistXiMinusPTMASS);
	TH2F *fHistXiMinusPtMass = new TH2F("fHistXiMinusPtMass", "XiMinus PT distribution and Massnoteffcor;p_{T} [GeV];Mass (GeV/^{2})",nPtBinsHXi,PtBinsHXi,
			nBins,1.29,1.37);
	tAnalysisXiMinus->Add(fHistXiMinusPtMass);


	fOutput->Add(tAnalysisXiMinus);

}
//===============================================================================
void AliAnalysisTaskLambdaK0s::AddQAV0Candidates(){
	TList *tQACandidates;
	TH2F *tDL, *tDCA, *tDCA2PV, *tCTP, *tD0, *tDCA2PVAft;
	TH3F *tAPBef, *tAPAft;

	tQACandidates = new TList();
	tQACandidates->SetOwner();
	tQACandidates->SetName("V0Candidates");

	tDL   = new TH2F("BefDL", "DL;[cm];Pt [GeV]",  50, 0, 10,  24, 0, 12);
	tQACandidates->Add(tDL);

	tDCA  = new TH2F("BefDCA", "DCA;[cm];Pt [GeV]", 50,0,1.001, 24,0,12);
	tQACandidates->Add(tDCA);

	tDCA2PV  = new TH2F("BefDCA2PV", "DCA to Prim. Vertex;[cm]; Pt [GeV]", 80, 0, 4., 24,0,12);
	tQACandidates->Add(tDCA2PV);

	tCTP  = new TH2F("BefCTP", "CTP;;Pt [GeV]",   80,0.91,1.001, 24,0,12);
	tQACandidates->Add(tCTP);

	tD0   = new TH2F("BefD0",  "D0;[cm];Pt [GeV]",  50,-1,+1, 24,0,12);
	tQACandidates->Add(tD0);

	tAPBef   = new TH3F("BefAP",  "AP;#alpha;q_{t}[GeV];Pt [GeV]",
			80,-1,+1, 90, 0, 0.3, 24, 0, 12);
	tQACandidates->Add(tAPBef);

	tAPAft  = new TH3F("AftAP","AP;#alpha;q_{t}[GeV];Pt [GeV]",
			80,-1, +1, 90, 0, 0.3, 24, 0, 12);
	tQACandidates->Add(tAPAft);

	tDCA2PVAft  = new TH2F("AftDCA2PV", "DCA to Prim. Vertex;[cm]; Pt [GeV]", 80, 0, 4., 24,0,12);
	tQACandidates->Add(tDCA2PVAft);

	fOutput->Add(tQACandidates);  
}


//________________________________________________________________________
void AliAnalysisTaskLambdaK0s::Terminate(Option_t *)
{
	// Draw result to screen, or perform fitting, normalizations
	// Called once at the end of the query

	fOutput = dynamic_cast<TList*>(GetOutputData(1));
	if (!fOutput) { AliError("Could not retrieve TList fOutput"); return; }

	// NEW HISTO should be retrieved from the TList container in the above way,
	// so it is available to draw on a canvas such as below
}

//_________________________________________________________________________
void AliAnalysisTaskLambdaK0s::UserExec(Option_t *)
{
	

	AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
	AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());


	UInt_t fSelectMask= inputHandler->IsEventSelected();
	if (!fSelectMask) return;

	//Bool_t isINT7selected = fSelectMask& AliVEvent::kINT7;
	//if (!isINT7selected) return;

	AliAODEvent* aod = dynamic_cast<AliAODEvent*>(inputHandler->GetEvent());
	if(!aod) return;
	fPIDResponse = inputHandler->GetPIDResponse(); 
	((TH1F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fhEventAllSel"))->Fill(0);
	((TH1F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fhEventBf"))->Fill(0);
	((TH1F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fHistMultiMain"))->Fill(aod->GetNumberOfTracks());

	//remove pile-up events=============================================================== 
	if (!fEventCuts->AcceptEvent(aod))
	{	 
		PostData(1, fOutput);

		return;
	}

	((TH1F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fhEventAllSel"))->Fill(1);
	const AliVVertex * primVertex = fEventCuts->GetPrimaryVertex();
	if (! primVertex) return;
	if ( ( TMath::Abs(primVertex->GetZ()) ) >= fPrimVertexCut ) return ;
	Double_t lPVx = primVertex->GetX();
	Double_t lPVy = primVertex->GetY();
	Double_t lPVz = primVertex->GetZ();
	if (TMath::Abs(lPVx)<fVtxXMin && TMath::Abs(lPVy)<fVtxYMin&& TMath::Abs(lPVz)<fVtxZMin) return;
	Short_t binVertex = Short_t((lPVz+8.)/2.);      

	((TH1F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fhEventAllSel"))->Fill(2);


	// aod ensurance
	AliAODHeader *aodHeader = dynamic_cast<AliAODHeader*>(aod->GetHeader());
	if(!aodHeader) AliFatal("Not a standard AOD");
	if(!aodHeader) return;

	Float_t lPercentile = 0; 
	AliMultSelection *MultSelection = 0x0; 
	MultSelection = (AliMultSelection *) aod->FindListObject("MultSelection");
	if( !MultSelection) {
		//If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
		AliWarning("AliMultSelection object not found!");
	}else{
		lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
	}
	if ((lPercentile < fMinCent)||(lPercentile >= fMaxCent )) return;//liu change

	((TH1F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fhEventAllSel"))->Fill(3);
	((TH1F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fHistVtxAf"))->Fill(lPVz);
	((TH1F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fHistCentAf"))->Fill(lPercentile);



	//begin processing MC data=================================================================================================== 
	Int_t iMC = 0;
	if(fAnalysisMC){
		 if(fIsPileUpCuts){
		//// out of bunch Generated pile up
		AliAODMCHeader *aodMCheader = (AliAODMCHeader*)aod->FindListObject(AliAODMCHeader::StdBranchName());
		if(!aodMCheader) {
		printf("AliAnalysisTaskSEHFTreeCreator::UserExec: MC header branch not found!\n");
		return;
		}

		Bool_t isPileupInGeneratedEvent = kFALSE;
		isPileupInGeneratedEvent = AliAnalysisUtils::IsPileupInGeneratedEvent(aodMCheader,"Hijing");
		if(isPileupInGeneratedEvent) return;

		Float_t vzMC = aodMCheader->GetVtxZ();
		if (TMath::Abs(vzMC) >= fPrimVertexCut) return;   
		}
		//retreive MC particles from event
		fMCArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());

		if(!fMCArray){
			Printf("No MC particle branch found");
			return;
		}
		//begin MC truth--------------------------------------------------------------------------------- 


		Int_t nMCTracks = fMCArray->GetEntriesFast();
		for (iMC = 0; iMC < nMCTracks; iMC++){
			AliAODMCParticle *mcTrack = (AliAODMCParticle*)fMCArray->At(iMC);
			if (!mcTrack) {
				Error("ReadEventAODMC", "Could not receive particle %d", iMC);
				continue;
			}

			if ((mcTrack->GetStatus() == 21) 
					|| (mcTrack->GetPdgCode() == 443 && mcTrack->GetMother() == -1)) 
				break; 

			//track processing
			Double_t mcTrackEta = mcTrack->Eta();
			Double_t mcTrackPt = mcTrack->Pt();
			Double_t mcTrackPhi = mcTrack->Phi();
			Double_t mcTrackY = mcTrack->Y();


			if(!mcTrack->IsPhysicalPrimary()) continue;
			if(TMath::Abs(mcTrackEta) < fEtaCut && mcTrack->Charge() != 0){
				fHistMCtruthTrkPtEta->Fill(mcTrackPt,mcTrackEta);
				fHistMCtruthTrk[binVertex]->Fill(mcTrackPt,mcTrackEta, mcTrackPhi);
			}
			//V0 and cascade processing
              
			
			if(fRapidity){
			    if(TMath::Abs( mcTrackY )>= 0.5 &&  mcTrackPt <1) continue;
            }else {
				 if(TMath::Abs(mcTrackEta) >= fV0Eta && mcTrackPt < 1) continue;
			}

                     if(TMath::Abs(mcTrack->GetPdgCode()) == 310){
				     //K0S
				     fHistMCtruthK0sPt->Fill(mcTrackPt, mcTrackEta);
				      fHistMCtruthK0s[binVertex]->Fill(mcTrackPt, mcTrackEta, mcTrackPhi);
		         	}
			
			    	else if(mcTrack->GetPdgCode() == 3122){
				    	 //Lambda
				   		 fHistMCtruthLambdaPt->Fill(mcTrackPt, mcTrackEta);// mcTrackY
				   		 fHistMCtruthLambda[binVertex]->Fill(mcTrackPt, mcTrackEta, mcTrackPhi);
					}
					else if(mcTrack->GetPdgCode() == -3122){
						//Anti-Lambda
						fHistMCtruthAntiLambdaPt->Fill(mcTrackPt, mcTrackEta);
						fHistMCtruthAntiLambda[binVertex]->Fill(mcTrackPt, mcTrackEta, mcTrackPhi);
					}
						//------------------------------------------------------------------------------------------------------------
					else if (mcTrack->GetPdgCode() == 3312 || mcTrack->GetPdgCode() == -3312){
						//Xi-
						fHistMCtruthXiMinusPt->Fill(mcTrackPt, mcTrackEta);
						//add hist
						fHistMCtruthXiMinus[binVertex]->Fill(mcTrackPt, mcTrackEta, mcTrackPhi);
						fHistMCXiMinusPt->Fill(mcTrackPt);
					}  
				

			

			
		}//end MC truth-------------------------------------------------------------------------------------------


		// begin access the reconstructed data---------------------------------------------------------------------
		Int_t nRecTracks = aod->GetNumberOfTracks();
		for (Int_t i = 0; i < nRecTracks; i++){
			AliAODTrack* tr = (AliAODTrack*)aod->GetTrack(i);
			if(tr->Charge() == 0.) continue;
			if (!(IsGoodPrimaryTrack(tr))) continue;
			if(tr->GetLabel() ==-1) {
				continue;
			}

			if(TMath::Abs(tr->GetLabel()) > iMC ) {

				continue;
			}

			//out of Bunch rejection trk 

			if(fRejectTrackPileUp&&(!(tr->HasPointOnITSLayer(0) || tr->HasPointOnITSLayer(1)|| tr->GetTOFBunchCrossing()==0 ))) continue;

			fHistRCTrkPtEta->Fill(tr->Pt(), tr->Eta());
			fHistRCTrk[binVertex]->Fill(tr->Pt(), tr->Eta(), tr->Phi());

			if(! (static_cast<AliAODMCParticle*>(fMCArray->At(TMath::Abs(tr->GetLabel()))))->IsPhysicalPrimary()){
				fHistRCSecTrkPtEta->Fill(tr->Pt(), tr->Eta());
				fHistRCSecTrk[binVertex]->Fill(tr->Pt(), tr->Eta(), tr->Phi());
			}
			else{
				fHistRCPriTrkPtEta->Fill(tr->Pt(), tr->Eta());
				fHistRCPriTrk[binVertex]->Fill(tr->Pt(), tr->Eta(), tr->Phi());
			}
		}//end for MC reconstructed data(track)---------------------------
	}//end for MC data process(truth for all,reconstructed for track)-----



	// Track selection loop
	Int_t nTracks = aod->GetNumberOfTracks();
	// new tracks array
	TObjArray * selectedTracks = new TObjArray;
	selectedTracks->SetOwner(kTRUE);

	for (Int_t i = 0; i < nTracks; i++)
	{
		AliAODTrack* tr = (AliAODTrack*)aod->GetTrack(i);
		if (!(IsGoodPrimaryTrack(tr))) continue;
		if (tr->Charge()==0) continue;
		if (TMath::Abs(tr->Eta()) > fEtaCut) continue;

		((TH3F*)((TList*)fOutput->FindObject("AnalysisTrk"))->FindObject("fHistTrk"))->Fill(tr->Pt(), tr->Eta(), tr->Phi());
		((TH2F*)((TList*)fOutput->FindObject("AnalysisTrk"))->FindObject("fHistTrkVzvsEta"))->Fill(tr->Eta(),lPVz);

		if ( tr->Pt() < fPtTrigMin || tr->Pt() > fPtTrigMax) continue;

		if(fAnalysisMC) {
			if(tr->GetLabel() ==-1) continue;
			if(! (static_cast<AliAODMCParticle*>(fMCArray->At(TMath::Abs(tr->GetLabel()))))->IsPhysicalPrimary())
				continue;
			if(TMath::Abs( tr->GetLabel())  > iMC) continue;
		}
		selectedTracks->Add(tr);
	}
	//---------------------------------

	const AliAODVertex *primaryBestAODVtx = aod->GetPrimaryVertex();
	primaryBestAODVtx->GetXYZ(fBestPrimaryVtxPos);

	TObjArray * selectedK0s = new TObjArray();
	selectedK0s->SetOwner(kTRUE);

	TObjArray * selectedLambda = new TObjArray();
	selectedLambda->SetOwner(kTRUE);

	TObjArray * selectedAntiLambda = new TObjArray();
	selectedAntiLambda->SetOwner(kTRUE);

	TObjArray * assocParticles = new TObjArray();
	assocParticles->SetOwner(kTRUE);

	//loop for V0s
	Int_t nV0sTot = aod->GetNumberOfV0s();

	for (Int_t iV0 = 0; iV0 < nV0sTot; iV0++) {
		AliAODv0 *v0=aod->GetV0(iV0);
		if (!v0) continue;

		if (!IsGoodV0(v0)) continue;

		Double_t lEta  = v0->PseudoRapV0();
		Double_t lRapLabK0Short = v0->RapK0Short();
		Double_t lRapLabLambda = v0->RapLambda();                                                                                                     
        if(fRapidity){
			if ( -0.5 >= lRapLabK0Short || lRapLabK0Short >= 0.5 || -0.5>= lRapLabLambda || lRapLabLambda >= 0.5  )  continue; 
		}else{
            if (  TMath::Abs(lEta) >= fV0Eta ) continue;
			 
		}



		const AliAODTrack *ntrack=(AliAODTrack *)v0->GetDaughter(1);
		const AliAODTrack *ptrack=(AliAODTrack *)v0->GetDaughter(0);
		if (!ptrack || !ntrack) continue;


		Int_t iPos, iNeg;
		if(ptrack->Charge() > 0){
			iPos = 0; iNeg = 1;
		}
		else{
			iPos = 1; iNeg = 0;
		}

		// Decay vertex                                                             
		Double_t xyz[3]; v0->GetSecondaryVtx(xyz);
		Double_t dx=xyz[0]-lPVx, dy=xyz[1]-lPVy, dz=xyz[2]-lPVz; 

		// Momentum: 2D & 3D                                                        
		Double_t lPt=TMath::Sqrt(v0->Pt2V0());

		if(lPt < fPtAssoMin || lPt > fPtAssoMax) continue;

		// Decay length: 2D & 3D                                                    
		//Double_t lt=TMath::Sqrt(dx*dx + dy*dy);
		Float_t dl=TMath::Sqrt(dx*dx + dy*dy + dz*dz);                            

		Double_t dlK = 0.4977*dl/lPt;//ctau
		Double_t dlL = 1.1157*dl/lPt;

		Double_t CTV0 = v0->Ct(3122,fBestPrimaryVtxPos);                                                                                              

		// ctau                                                                             
		Bool_t ctK=kTRUE;  if (dlK >fPLtimeK0s || dlK < fMinCtau*2.68) ctK=kFALSE;
		Bool_t ctL=kTRUE;  if (dlL >fPLtimeLambda|| dlL < fMinCtau*7.89) ctL=kFALSE;


		Float_t xyn = v0->DcaNegToPrimVertex();
		Float_t xyp = v0->DcaPosToPrimVertex();
		//  ---- V0 candidate properties:      
		Double_t lAlphaV0      =  v0->AlphaV0();
		if(v0->ChargeProng(iPos) <0) lAlphaV0 = -lAlphaV0;
		Double_t lPtArmV0      =  v0->PtArmV0();
		fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0);
		//cosine pointing angle
		Double_t lCPA = v0->CosPointingAngle(fBestPrimaryVtxPos);
		//DCA to PV
		Double_t lDCA2PV = v0->DcaV0ToPrimVertex();
		// phi                                        
		Double_t lPhi  = v0->Phi();
		Double_t massK0s = v0->MassK0Short();
		Double_t massLambda = v0->MassLambda();
		Double_t massAntiLambda = v0->MassAntiLambda();

		Float_t nsigPosPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrack,AliPID::kPion));
		Float_t nsigPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrack,AliPID::kProton));

		Float_t nsigNegPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrack,AliPID::kPion));
		Float_t nsigNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrack,AliPID::kProton));

		//check whether it is K0s candidates
		if(ctK && lPtArmV0 > 0.2*TMath::Abs(lAlphaV0)&& nsigPosPion < fNsigma && nsigNegPion < fNsigma && xyn>fDCANegtoPrimVertexMinK0s &&xyp>fDCAPostoPrimVertexMinK0s){
			//	if(TMath::Abs(massLambda - fMassMean[1]) > 5.*fMassRes[1] 
			//			&& TMath::Abs(massAntiLambda - fMassMean[1]) > 5.*fMassRes[1]){
			selectedK0s->Add(v0);
			if( (massK0s > 0.40) && (     massK0s < 0.58))
			{
				((TH2F*)((TList*)fOutput->FindObject("AnalysisK0s"))->FindObject("fHistMassK0s"))->Fill(massK0s,lPt);//do the eff at the outside
			}    
			if(TMath::Abs(massK0s - fMassMean[0]) < 8.*fMassRes[0]){

				((TH3F*)((TList*)fOutput->FindObject("AnalysisK0s"))->FindObject("fHistK0s"))->Fill(lPt, lEta, lPhi);
				((TH1F*)((TList*)fOutput->FindObject("AnalysisK0s"))->FindObject("fHistK0sEta"))->Fill(lEta);	
				((TH1F*)((TList*)fOutput->FindObject("AnalysisK0s"))->FindObject("fHistK0sPhi"))->Fill(lPhi);
				((TH2F*)((TList*)fOutput->FindObject("AnalysisK0s"))->FindObject("fHistK0sVzvsEta"))->Fill(lEta,lPVz);  
				if(fAnalysisMC && IsMCV0Primary(v0, 0)){
					fHistRCK0sPt->Fill(lPt, lEta);
					fHistRCK0s[binVertex]->Fill(lPt, lEta, lPhi);
				}
				assocParticles->Add(new AliV0XiParticleall(lPt, lPhi, lEta, 
							v0->MassK0Short(), 0));}
				
		}

		// check whether it is Lambda candidates
		if(ctL && lCPA > fLambdaCPA && lAlphaV0 > 0 && lDCA2PV < fLambdaDCA
				&& nsigPosProton < fNsigma && nsigNegPion < fNsigma && xyn>fDCANegtoPrimVertexMinLambda &&xyp>fDCAPostoPrimVertexMinLambda ){
			selectedLambda->Add(v0);

			if( (massLambda > 1.04) && (massLambda      < 1.21))
			{ 
				((TH2F*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistMassLambda"))->Fill(massLambda,lPt);  
				((TH3F*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistMassDCALambda"))->Fill(massLambda,lPt,lDCA2PV);
			}  
			if(TMath::Abs(massLambda - fMassMean[1]) < 8.*fMassRes[1]){
				((TH2F*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistDCALambda"))->Fill(lPt,lDCA2PV);	  

				((TH3F*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistLambda"))->Fill(lPt, lEta, lPhi);
				((TH1F*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistLambdaEta"))->Fill(lEta);
				((TH1F*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistLambdaPhi"))->Fill(lPhi);
				((TH2F*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistLambdaVzvsEta"))->Fill(lEta,lPVz);

				Double_t r=TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);//fiducial volume
				if(fAnalysisMC){

					fHistallRCLambdaPt->Fill(v0->Ct(3122,fBestPrimaryVtxPos), lPt);
					fHistallDCALambdaPt->Fill(lDCA2PV, lPt);
					fHistallDLLambdaPt->Fill(v0->DecayLengthV0( fBestPrimaryVtxPos ), lPt);
					fHistallRLambdaPt->Fill(r, lPt);
					if(IsMCV0Primary(v0, 1)){
						fHistRCLambdaPt->Fill(lPt, lEta);
						fHistRCLambda[binVertex]->Fill(lPt, lEta, lPhi);

						fHistRCPrimLambdaDCAPt->Fill(lDCA2PV, lPt);
						fHistRCPrimLambdaDLPt->Fill(v0->DecayLengthV0( fBestPrimaryVtxPos ), lPt);
						fHistRCPrimLambdaCTPt->Fill(v0->Ct(3122,fBestPrimaryVtxPos), lPt);

						fHistRCPrimLambdaRPt->Fill(r, lPt);
					}
					else if(IsMCV0FromXi(v0, 1)){
						fHistRCSecLambdaDCAPt->Fill(lDCA2PV, lPt);	
						fHistRCSecLambdaDLPt->Fill(v0->DecayLengthV0( fBestPrimaryVtxPos ), lPt);
						fHistRCSecLambdaCTPt->Fill(v0->Ct(3122,fBestPrimaryVtxPos), lPt);
						fHistRCSecLambdaRPt->Fill(r, lPt);

						//----------------------------------------------------
						TClonesArray *mcArray = (TClonesArray*)aod->FindListObject(AliAODMCParticle::StdBranchName());
						if(!mcArray)return;

						Int_t myTrackPosLabel        = TMath::Abs(ptrack->GetLabel());
						Int_t myTrackNegLabel        = TMath::Abs(ntrack->GetLabel());

						AliAODMCParticle *mcPosTrack = (AliAODMCParticle*)mcArray->At(myTrackPosLabel);
						if(!mcPosTrack)continue;
						AliAODMCParticle *mcNegTrack = (AliAODMCParticle*)mcArray->At(myTrackNegLabel);
						if(!mcNegTrack)continue;

						Int_t PosDaughterPdg = mcPosTrack->GetPdgCode();
						Int_t NegDaughterPdg = mcNegTrack->GetPdgCode();

						Int_t myTrackPosMotherLabel = mcPosTrack->GetMother();
						Int_t myTrackNegMotherLabel = mcNegTrack->GetMother();

						if ((myTrackPosMotherLabel==-1)||(myTrackNegMotherLabel==-1)) continue;
						if (myTrackPosMotherLabel!=myTrackNegMotherLabel) continue;

						AliAODMCParticle *mcPosMother = (AliAODMCParticle*)mcArray->At(myTrackPosMotherLabel);
						if(!mcPosMother)continue;
						Int_t MotherPdg  = mcPosMother->GetPdgCode();
						Bool_t IsPrime   = mcPosMother->IsPhysicalPrimary();

						Int_t myGrandMotherLabel = mcPosMother->GetMother();
						if (myGrandMotherLabel != -1)
						{

							AliAODMCParticle *mcGrandMother = (AliAODMCParticle*)mcArray->At(myGrandMotherLabel);
							Int_t GrandMotherPdg     = mcGrandMother->GetPdgCode(); 

							if ((MotherPdg == 3122 && PosDaughterPdg== 2212 && NegDaughterPdg== -211 && (GrandMotherPdg==3322 ||GrandMotherPdg==3312))) 
							{

								((TH2F*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistXiPtVSLambdaPt"))->Fill(mcGrandMother->Pt(),lPt);

							}

						}		





					}
				}
				assocParticles->Add(new AliV0XiParticleall(lPt, lPhi, lEta,
							massLambda, 1));
			}
		}

		// check whether it is AntiLambda candidates   
		if(ctL && lCPA > fLambdaCPA && lAlphaV0 < 0 && lDCA2PV < fLambdaDCA
				&& nsigPosPion < fNsigma && nsigNegProton < fNsigma && xyn>fDCANegtoPrimVertexMinAntiLambda &&xyp>fDCAPostoPrimVertexMinAntiLambda){
			selectedAntiLambda->Add(v0);

			if(TMath::Abs(massAntiLambda - fMassMean[1]) < 8.*fMassRes[1]){
				((TH3F*)((TList*)fOutput->FindObject("AnalysisAntiLambda"))->FindObject("fHistAntiLambda"))->Fill(lPt, lEta, lPhi);
				((TH1F*)((TList*)fOutput->FindObject("AnalysisAntiLambda"))->FindObject("fHistAntiLambdaEta"))->Fill(lEta);
				((TH1F*)((TList*)fOutput->FindObject("AnalysisAntiLambda"))->FindObject("fHistAntiLambdaPhi"))->Fill(lPhi);
				((TH2F*)((TList*)fOutput->FindObject("AnalysisAntiLambda"))->FindObject("fHistAntiLambdaVzvsEta"))->Fill(lEta,lPVz);	 


				if(fAnalysisMC){
					if( IsMCV0Primary(v0, 2)){
						fHistRCAntiLambdaPt->Fill(lPt, lEta);
						fHistRCAntiLambda[binVertex]->Fill(lPt, lEta, lPhi);
						fHistRCPrimAntiLambdaDCAPt->Fill(lDCA2PV, lPt);
					}
					else if(IsMCV0FromXi(v0, 2)){
						fHistRCSecAntiLambdaDCAPt->Fill(lDCA2PV, lPt);
						//--------------------------------------------------------------------
						TClonesArray *mcArray = (TClonesArray*)aod->FindListObject(AliAODMCParticle::StdBranchName());
						if(!mcArray)return;

						Int_t myTrackPosLabel        = TMath::Abs(ptrack->GetLabel());
						Int_t myTrackNegLabel        = TMath::Abs(ntrack->GetLabel());

						AliAODMCParticle *mcPosTrack = (AliAODMCParticle*)mcArray->At(myTrackPosLabel);
						if(!mcPosTrack)continue;
						AliAODMCParticle *mcNegTrack = (AliAODMCParticle*)mcArray->At(myTrackNegLabel);
						if(!mcNegTrack)continue;

						Int_t PosDaughterPdg = mcPosTrack->GetPdgCode();
						Int_t NegDaughterPdg = mcNegTrack->GetPdgCode();

						Int_t myTrackPosMotherLabel = mcPosTrack->GetMother();
						Int_t myTrackNegMotherLabel = mcNegTrack->GetMother();

						if ((myTrackPosMotherLabel==-1)||(myTrackNegMotherLabel==-1)) continue;
						if (myTrackPosMotherLabel!=myTrackNegMotherLabel) continue;

						AliAODMCParticle *mcPosMother = (AliAODMCParticle*)mcArray->At(myTrackPosMotherLabel);
						if(!mcPosMother)continue;
						Int_t MotherPdg  = mcPosMother->GetPdgCode();
						Bool_t IsPrime   = mcPosMother->IsPhysicalPrimary();

						Int_t myGrandMotherLabel = mcPosMother->GetMother();

						if (myGrandMotherLabel != -1)
						{

							AliAODMCParticle *mcGrandMother = (AliAODMCParticle*)mcArray->At(myGrandMotherLabel);
							Int_t GrandMotherPdg     = mcGrandMother->GetPdgCode();
							if ((MotherPdg == -3122 && PosDaughterPdg== 211 && NegDaughterPdg== -2212 && (GrandMotherPdg== 3322 ||GrandMotherPdg==-3312))) 
							{
								((TH2F*)((TList*)fOutput->FindObject("AnalysisAntiLambda"))->FindObject("fHistXiPtVSAntiLambdaPt"))->Fill(mcGrandMother->Pt(),lPt);	
							}
						}
					}
				}
				assocParticles->Add(new AliV0XiParticleall(lPt, lPhi, lEta,
							massAntiLambda, 2));
			}
		}
	}
	//------------------------------------  Cascade  -------------------------------------------------------------------

	TObjArray * selectedXiMinus = new TObjArray();
	selectedXiMinus->SetOwner(kTRUE);

	Int_t nCascades = aod->GetNumberOfCascades();
	// This is the begining of the Cascade loop
	for(Int_t iXi = 0; iXi < nCascades; iXi++){
		Double_t invMassXi = 0.;
		Double_t invMassOmega = 0.;
		Double_t dcaXiDaughters = -1.;//dca between Xi daughters
		Double_t dcaPosToPrimaryVtxXi = -1.;//dca of V0 neg daughter to Prim. Vertex, in cascade
		Double_t dcaNegToPrimaryVtxXi = -1.;//dca of V0 pos daughter to Prim. Vertex, in cascade
		Double_t dcaBachToPrimaryVtxXi = -1.;// dca of bachelor to primary vertex, in cascade 
		Double_t dcaV0ToPrimaryVtxXi = -1.;//dca of V0 to primary vertex, in cascade 
		Double_t dcaV0DaughtersXi = -1.;// dca between V0 daughters, in cascade 
		Double_t XiCosOfPointingAngle = -1.;
		Double_t V0CosOfPointingAngleXi = -1.;//Cos of V0 Pointing Angle, in cascade
		Double_t posXi[3] = {-1000., -1000., -1000.};
		Double_t XiRadius = -1000.;
		Double_t posV0Xi[3] = {-1000., -1000., -1000.};//V0 pos daughters, in cascade
		Double_t V0RadiusXi = -1000.;//V0 decay radius, in cascade
		Double_t invMassLambdaAsCascDghter = 0.; 
		Double_t rapXi = -20.;
		Double_t etaXi = -20.;

		AliAODcascade *xi = aod->GetCascade(iXi);
		if (!xi) continue;
		dcaXiDaughters = xi->DcaXiDaughters();

		Bool_t isBachelorPionForTPC = kFALSE;
		Bool_t isNegPionForTPC = kFALSE;
		Bool_t isPosPionForTPC = kFALSE;
		Bool_t isNegProtonForTPC = kFALSE;
		Bool_t isPosProtonForTPC = kFALSE;

		posXi[0] = xi->DecayVertexXiX();
		posXi[1] = xi->DecayVertexXiY();
		posXi[2] = xi->DecayVertexXiZ();
		XiRadius = TMath::Sqrt(posXi[0]*posXi[0]
				+posXi[1]*posXi[1]
				+posXi[2]*posXi[2]);

		XiCosOfPointingAngle = xi->CosPointingAngleXi(fBestPrimaryVtxPos[0], fBestPrimaryVtxPos[1], fBestPrimaryVtxPos[2]);

		AliAODTrack *pTrkXi = dynamic_cast<AliAODTrack*>( xi->GetDaughter(0) );
		AliAODTrack *nTrkXi = dynamic_cast<AliAODTrack*>( xi->GetDaughter(1) );
		AliAODTrack *bTrkXi = dynamic_cast<AliAODTrack*>( xi->GetDecayVertexXi()->GetDaughter(0) ); 


		if(!pTrkXi || !nTrkXi || !bTrkXi) continue;

		UInt_t idxPosXi  = (UInt_t) TMath::Abs( pTrkXi->GetID() );
		UInt_t idxNegXi  = (UInt_t) TMath::Abs( nTrkXi->GetID() );
		UInt_t idxBach   = (UInt_t) TMath::Abs( bTrkXi->GetID() );

       // Int_t labelPos    = pTrackXi->GetLabel();//add
       // Int_t labelNeg    = nTrackXi->GetLabel();
       // Int_t labelBach = bachTrackXi->GetLabel();





		if(idxBach == idxNegXi || idxBach == idxPosXi) continue;

		if( !IsGoodDaughterTrack(pTrkXi) || !IsGoodDaughterTrack(nTrkXi)  || !IsGoodDaughterTrack(bTrkXi) ) continue;


		Double_t massxi = xi->MassXi();
		invMassXi = xi->MassXi();
		invMassOmega = xi->MassOmega();
		
		if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bTrkXi, AliPID::kPion))<4.)
			isBachelorPionForTPC = kTRUE;

		//Negative V0 daughter                                                               
		if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrkXi, AliPID::kPion))<4.)
			isNegPionForTPC = kTRUE;
		if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrkXi, AliPID::kProton))<4.)
			isNegProtonForTPC = kTRUE;

		//Positive V0 daughter                                                               
		if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrkXi, AliPID::kPion))<4.)
			isPosPionForTPC = kTRUE;
		if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrkXi, AliPID::kProton))<4.)
			isPosProtonForTPC = kTRUE; 

		//V0                                                                                 
		posV0Xi[0] = xi->DecayVertexV0X();
		posV0Xi[1] = xi->DecayVertexV0Y();
		posV0Xi[2] = xi->DecayVertexV0Z();
		V0RadiusXi = TMath::Sqrt(posV0Xi[0]*posV0Xi[0]
				+posV0Xi[1]*posV0Xi[1]
				+posV0Xi[2]*posV0Xi[2]);

		if(xi->ChargeXi() < 0)
			invMassLambdaAsCascDghter = xi->MassLambda();
		else
			invMassLambdaAsCascDghter = xi->MassAntiLambda();


		V0CosOfPointingAngleXi = xi->CosPointingAngle(fBestPrimaryVtxPos);

		dcaV0DaughtersXi = xi->DcaV0Daughters();

		dcaV0ToPrimaryVtxXi = xi->DcaV0ToPrimVertex();
		dcaBachToPrimaryVtxXi = xi->DcaBachToPrimVertex();


		dcaPosToPrimaryVtxXi = xi->DcaPosToPrimVertex();
		dcaNegToPrimaryVtxXi = xi->DcaNegToPrimVertex();
		//VO CUT
		if(TMath::Abs(invMassLambdaAsCascDghter-1.11568) > 0.008) continue;
		if(dcaV0DaughtersXi > 1.) continue;//2                                                
		if(dcaPosToPrimaryVtxXi < 0.2) continue; //0.03                                       
		if(dcaNegToPrimaryVtxXi < 0.2) continue; //0.04
		if(V0CosOfPointingAngleXi < 0.97) continue;
		if(V0RadiusXi < 3.0 ) continue;//1      

		// XI CUT 

		if(dcaXiDaughters > 1.3) continue;//1 
		if(dcaV0ToPrimaryVtxXi < 0.1) continue;//0.05                                           
		if(dcaBachToPrimaryVtxXi < 0.1) continue;//0.03
		if(XiCosOfPointingAngle < 0.97) continue; //0.98
		if(XiRadius < 1.2 ) continue; //0.4                                      


		// phi                                        
		Double_t lPhi  = xi->Phi();

		// Momentum: 2D & 3D                                                        
		Double_t lPt=TMath::Sqrt(xi->Pt2V0());

		if(lPt < 0.0 || lPt > 100.) continue;//100

		Double_t XiPx = 0., XiPy = 0., XiPz = 0.;                                                                                  
		Double_t XiPt = 0.;
		Double_t XiPtot = 0.;
		//-----------------------                                                            
		XiPx = xi->MomXiX();
		XiPy = xi->MomXiY();
		XiPz = xi->MomXiZ();
		XiPt = TMath::Sqrt(XiPx*XiPx + XiPy*XiPy);
		XiPtot= TMath::Sqrt(XiPx*XiPx + XiPy*XiPy + XiPz*XiPz);

		
		rapXi = xi->RapXi();
		etaXi = 0.5*TMath::Log((XiPtot+XiPz)/(XiPtot-XiPz+1.e-13));   
        //    lThetaXi = xi->Theta();//add
	    //   lEtaXi = -TMath::Log(TMath::Tan(lThetaXi/2));


		if (TMath::Abs(etaXi)>0.8) continue;
		if(lPt < 1 || lPt > 10.0) continue;


		if ((TMath::Abs(etaXi) < 0.7 && isBachelorPionForTPC) && ((xi->ChargeXi() < 0 && isPosProtonForTPC && isNegPionForTPC)|| (xi->ChargeXi() > 0 && isNegProtonForTPC && isPosPionForTPC) ) ){

			if( TMath::Abs(invMassOmega-1.67245) > 0.008 ){

				selectedXiMinus->Add(xi);

				((TH1F*)((TList*)fOutput->FindObject("AnalysisXiMinus"))->FindObject("fHistXiMinusEta"))->Fill(etaXi);//eta???
				((TH1F*)((TList*)fOutput->FindObject("AnalysisXiMinus"))->FindObject("fHistXiMinusPhi"))->Fill(lPhi);
				((TH3F*)((TList*)fOutput->FindObject("AnalysisXiMinus"))->FindObject("fHistXiMinus"))->Fill(lPt, etaXi, lPhi);
				((TH2F*)((TList*)fOutput->FindObject("AnalysisXiMinus"))->FindObject("fHistXiMinusVzvsEta"))->Fill(etaXi,lPVz);
				((TH1F*)((TList*)fOutput->FindObject("AnalysisXiMinus"))->FindObject("fHistXiMinusPT"))->Fill(lPt);
				((TH3F*)((TList*)fOutput->FindObject("AnalysisXiMinus"))->FindObject("fHistMassXi"))->Fill(invMassXi, lPt, lPVz);
				((TH2F*)((TList*)fOutput->FindObject("AnalysisXiMinus"))->FindObject("fHistXiMinusPTMASS"))->Fill(lPt,massxi);


				//MC to estimate the reconstruction eff. of XiMinus 
				if(fAnalysisMC ){
					if(IsMCXiPrimary(xi) ){
						fHistRCXiMinusPt->Fill(lPt, etaXi);
						fHistRCXiMinus[binVertex]->Fill(lPt, etaXi, lPhi); //rapXi


					}

				}

				assocParticles->Add(new AliV0XiParticleall(XiPt, lPhi, etaXi, massxi ,3 ));

			}
		}
	}//end for selecting Xi(-and +);
	//-------------------------------------------------------------------------------------  

	AliEventPool* pool = 0;
	if(fEventMixing){
		pool = fPoolMgr->GetEventPool(lPercentile, lPVz);
		if (!pool)
			AliFatal(Form("No pool found for centrality = %f, zVtx = %f", lPercentile, lPVz));
	}



	//build correlations
	for(Int_t iTrk = 0; iTrk < selectedTracks->GetEntries(); iTrk++){           
		AliAODTrack* tr = (AliAODTrack*) selectedTracks->At(iTrk);                
		if(!tr) continue;                

		((TH2F*)((TList*)fOutput->FindObject("AnalysisTrk"))->FindObject("fHistTrigPtAll"))->Fill(tr->Pt(),lPVz);
		((TH2F*)((TList*)fOutput->FindObject("AnalysisTrk"))->FindObject("fTriggerEtaPhi"))->Fill(tr->Phi(),tr->Eta());
		((TH1F*)((TList*)fOutput->FindObject("AnalysisTrk"))->FindObject("fTriggerEta"))->Fill(tr->Eta());
		((TH1F*)((TList*)fOutput->FindObject("AnalysisTrk"))->FindObject("fTriggerPhi"))->Fill(tr->Phi());

		Int_t trID = tr->GetID() >= 0 ? tr->GetID() : -1-tr->GetID();

		if(fChV0){




			//asso K0s.........   
			for(Int_t i = 0; i < selectedK0s->GetEntries(); i++){
				AliAODv0 * v0 = (AliAODv0*)selectedK0s->At(i);
				if(!v0) continue;
				AliAODTrack *pTrkXi = dynamic_cast<AliAODTrack*>( v0->GetDaughter(0) );
				AliAODTrack *nTrkXi = dynamic_cast<AliAODTrack*>( v0->GetDaughter(1) );

				Int_t pTrkID = pTrkXi->GetID()  >= 0 ? pTrkXi->GetID() : -1-pTrkXi->GetID();
				Int_t nTrkID = nTrkXi->GetID() >= 0 ? nTrkXi->GetID() : -1-nTrkXi->GetID();

				if (v0->Pt() >= tr->Pt()) continue; 
				//remove auto correlations
				if ( pTrkID == trID || nTrkID == trID) continue;
				// Correlation part                                    
				//=================================== 
				Double_t dEta = v0->Eta() - tr->Eta();
				Double_t dPhi = v0->Phi() - tr->Phi();
				if ( dPhi > 1.5*TMath::Pi() ) dPhi -= 2.0*TMath::Pi();
				if ( dPhi < -0.5*TMath::Pi() ) dPhi += 2.0*TMath::Pi();

				// Filling correlation histograms and histograms for triggers counting    
				Double_t spSig[7] = {dPhi, dEta,tr->Pt(), v0->Pt(), v0->MassK0Short(),lPercentile,lPVz }; 
				if(fEffCorr){
					Double_t weight = fHistEffEtaPtK0s/*[binVertex]*/->Interpolate(v0->Eta(), v0->Pt());
					if(weight == 0){
						continue;
					}
					((THnSparseF*)((TList*)fOutput->FindObject("AnalysisK0s"))->FindObject("fHistdPhidEtaSibK0s"))->Fill(spSig, 1./weight);
				}
				else
					((THnSparseF*)((TList*)fOutput->FindObject("AnalysisK0s"))->FindObject("fHistdPhidEtaSibK0s"))->Fill(spSig);

			}
			//asso Lambda...... 
			for(Int_t i = 0; i < selectedLambda->GetEntries(); i++){
				AliAODv0 * v0 = (AliAODv0*)selectedLambda->At(i);
				if(!v0) continue;
				AliAODTrack *pTrkXi = dynamic_cast<AliAODTrack*>( v0->GetDaughter(0) );
				AliAODTrack *nTrkXi = dynamic_cast<AliAODTrack*>( v0->GetDaughter(1) );

				Int_t pTrkID = pTrkXi->GetID()  >= 0 ? pTrkXi->GetID() : -1-pTrkXi->GetID();
				Int_t nTrkID = nTrkXi->GetID() >= 0 ? nTrkXi->GetID() : -1-nTrkXi->GetID();

				if (v0->Pt() >= tr->Pt()) continue; 

				//remove auto correlations                
				if ( pTrkID == trID || nTrkID == trID) continue;
				// Correlation part
				//===================================
				Double_t dEta = v0->Eta() - tr->Eta();
				Double_t dPhi = v0->Phi() - tr->Phi();
				if ( dPhi > 1.5*TMath::Pi() ) dPhi -= 2.0*TMath::Pi();
				if ( dPhi < -0.5*TMath::Pi() ) dPhi += 2.0*TMath::Pi();
				Double_t xyz[3]; v0->GetSecondaryVtx(xyz);
				Double_t r2 =TMath::Sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);

				// Filling correlation histograms and histograms for triggers counting
				if(TMath::Abs(v0->MassLambda()-fMassMean[1]) < 3.*fMassRes[1])
				{
					((TH3F*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistLambdaDphiDCAPtSig"))->Fill(dPhi, v0->Pt(),                                                                                                                v0->DcaV0ToPrimVertex());
					((TH3F*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistLambdaDphiDLPtSig"))->Fill(dPhi, v0->Pt(),                                                                                                 v0->DecayLengthV0( fBestPrimaryVtxPos ));
					((TH3F*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistLambdaDphiCTPtSig"))->Fill(dPhi, v0->Pt(),                                                                                              v0->Ct(3122,fBestPrimaryVtxPos));
					((TH3F*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistLambdaDphiRPtSig"))->Fill(dPhi, v0->Pt(),r2);
				}
				else if(TMath::Abs(v0->MassLambda()-fMassMean[1]) > 4.*fMassRes[1]
						&& TMath::Abs(v0->MassLambda()-fMassMean[1]) < 7.*fMassRes[1])
				{		 
					((TH3F*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistLambdaDphiDCAPtBkg"))->Fill(dPhi, v0->Pt(), v0->DcaV0ToPrimVertex() );
					((TH3F*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistLambdaDphiDLPtBkg"))->Fill(dPhi, v0->Pt(),  v0->DecayLengthV0( fBestPrimaryVtxPos ) );
					((TH3F*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistLambdaDphiCTPtBkg"))->Fill(dPhi, v0->Pt(), v0->Ct(3122,fBestPrimaryVtxPos));
					((TH3F*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistLambdaDphiRPtBkg"))->Fill(dPhi, v0->Pt(),r2);


				}
				Double_t spSig[7] = {dPhi, dEta,tr->Pt(), v0->Pt(), v0->MassLambda(),lPercentile,lPVz };



				if(fEffCorr){
					Double_t weight = fHistEffEtaPtLambda/*[binVertex]*/->Interpolate(v0->Eta(), v0->Pt());

					if(weight == 0){
						continue;
					}
					((THnSparseF*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistdPhidEtaSibLambda"))     ->Fill(spSig, 1./weight);


				}
				else
					((THnSparseF*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistdPhidEtaSibLambda"))->Fill(spSig);
			}
			//asso AntiLambda ...................
			for(Int_t i = 0; i < selectedAntiLambda->GetEntries(); i++){
				AliAODv0 * v0 = (AliAODv0*)selectedAntiLambda->At(i);
				if(!v0) continue;
				AliAODTrack *pTrkXi = dynamic_cast<AliAODTrack*>( v0->GetDaughter(0) );
				AliAODTrack *nTrkXi = dynamic_cast<AliAODTrack*>( v0->GetDaughter(1) );

				Int_t pTrkID = pTrkXi->GetID()  >= 0 ? pTrkXi->GetID() : -1-pTrkXi->GetID();
				Int_t nTrkID = nTrkXi->GetID() >= 0 ? nTrkXi->GetID() : -1-nTrkXi->GetID();

				if (v0->Pt() >= tr->Pt()) continue; 

				//remove auto correlations     
				if ( pTrkID == trID || nTrkID == trID) continue;
				// Correlation part
				//===================================
				Double_t dEta = v0->Eta() - tr->Eta();
				Double_t dPhi = v0->Phi() - tr->Phi();
				if ( dPhi > 1.5*TMath::Pi() ) dPhi -= 2.0*TMath::Pi();
				if ( dPhi < -0.5*TMath::Pi() ) dPhi += 2.0*TMath::Pi();

				// Filling correlation histograms and histograms for triggers counting
				if(TMath::Abs(v0->MassAntiLambda()-fMassMean[1]) < 3.*fMassRes[1])
					((TH3F*)((TList*)fOutput->FindObject("AnalysisAntiLambda"))->FindObject("fHistAntiLambdaDphiDCAPtSig"))->Fill(dPhi, v0->Pt(), v0->DcaV0ToPrimVertex());
				else if(TMath::Abs(v0->MassAntiLambda()-fMassMean[1]) > 4.*fMassRes[1] 
						&& TMath::Abs(v0->MassAntiLambda()-fMassMean[1]) < 7.*fMassRes[1])
					((TH3F*)((TList*)fOutput->FindObject("AnalysisAntiLambda"))->FindObject("fHistAntiLambdaDphiDCAPtBkg"))->Fill(dPhi, v0->Pt(), v0->DcaV0ToPrimVertex());
				Double_t spSig[7] = {dPhi, dEta,tr->Pt(), v0->Pt(), v0->MassAntiLambda(),lPercentile,lPVz };
				if(fEffCorr){
					Double_t weight = fHistEffEtaPtAntiLambda/*[binVertex]*/->Interpolate(v0->Eta(), v0->Pt());
					if(weight == 0){
						continue;
					}
					((THnSparseF*)((TList*)fOutput->FindObject("AnalysisAntiLambda"))->FindObject("fHistdPhidEtaSibAntiLambda"))->Fill(spSig, 1./weight); }
				else
					((THnSparseF*)((TList*)fOutput->FindObject("AnalysisAntiLambda"))->FindObject("fHistdPhidEtaSibAntiLambda"))->Fill(spSig);
			}//end for AntiLambda

		}

		if(fChXi){

			//asso XiMius..........
			for(Int_t i = 0; i < selectedXiMinus->GetEntries(); i++){
				AliAODcascade * xi = (AliAODcascade*)selectedXiMinus->At(i);
				if(!xi) continue;
				AliAODTrack *pTrkXi = dynamic_cast<AliAODTrack*>( xi->GetDaughter(0) );
				AliAODTrack *nTrkXi = dynamic_cast<AliAODTrack*>( xi->GetDaughter(1) );
				AliAODTrack *bTrkXi = dynamic_cast<AliAODTrack*>( xi->GetDecayVertexXi()->GetDaughter(0) );

				Int_t pTrkID = pTrkXi->GetID()  >= 0 ? pTrkXi->GetID() : -1-pTrkXi->GetID();
				Int_t nTrkID = nTrkXi->GetID() >= 0 ? nTrkXi->GetID() : -1-nTrkXi->GetID();
				Int_t bTrkID = bTrkXi->GetID() >=0 ? bTrkXi->GetID() : -1-bTrkXi->GetID();

				if (xi->Pt() >= tr->Pt()) continue; 
				//remove auto correlations
				if ( pTrkID == trID || nTrkID == trID || bTrkID == trID ) continue;
				// Correlation part                                    
				//=================================== 
				Double_t dEta = xi->Eta() - tr->Eta();
				Double_t dPhi = xi->Phi() - tr->Phi();
				if ( dPhi > 1.5*TMath::Pi() ) dPhi -= 2.0*TMath::Pi();
				if ( dPhi < -0.5*TMath::Pi() ) dPhi += 2.0*TMath::Pi();

				// Filling correlation histograms and histograms for triggers counting    
				Double_t spSig[7] = {dPhi, dEta,tr->Pt(), xi->Pt(), xi->MassXi(),lPercentile,lPVz }; 
				if(fEffCorr){

					Double_t weight = fHistEffEtaPtXiMinus->Interpolate(xi->Eta(), xi->Pt());
					if(weight == 0){
						continue;
					}
					((THnSparseF*)((TList*)fOutput->FindObject("AnalysisXiMinus"))->FindObject("fHistdPhidEtaSibXiMinus"))->Fill(spSig, 1./weight);
				}
				else
					((THnSparseF*)((TList*)fOutput->FindObject("AnalysisXiMinus"))->FindObject("fHistdPhidEtaSibXiMinus"))->Fill(spSig);

			}
		}
		// Mixing ==============================================
		if (fEventMixing &&(pool->IsReady()|| pool->NTracksInPool() > fMixingTracks || pool->GetCurrentNEvents() >= 5))//15



		{
			Int_t nMix = pool->GetCurrentNEvents();

			for (Int_t jMix = 0; jMix < nMix; jMix++){
				// loop through mixing events
				TObjArray* bgTracks = pool->GetEvent(jMix);

				for (Int_t j = 0; j < bgTracks->GetEntriesFast(); j++){
					// mixing tracks loop 
					AliV0XiParticleall* assoc = (AliV0XiParticleall*) bgTracks->At(j);
					if(!assoc) continue;
				
					Double_t aPt = assoc->Pt();
					Double_t aMass = assoc->M();

					Double_t dEtaMix = assoc->Eta() - tr->Eta();
					Double_t dPhiMix = assoc->Phi() - tr->Phi();
					if ( dPhiMix > 1.5*TMath::Pi() ) dPhiMix -= 2.0*TMath::Pi();
					if ( dPhiMix < -0.5*TMath::Pi() ) dPhiMix += 2.0*TMath::Pi();
					if (aPt >= tr->Pt()) continue; 
					Double_t spMix[7] = {dPhiMix, dEtaMix, tr->Pt(),  aPt, aMass,lPercentile,lPVz };     

					if(fChV0){

						if(assoc->WhichCandidate() == 0){
							if(fEffCorr){
								Double_t weight = fHistEffEtaPtK0s->Interpolate(assoc->Eta(), 
										assoc->Pt());
								if(weight == 0){
									continue;
								}
								((THnSparseF*)((TList*)fOutput->FindObject("AnalysisK0s"))->FindObject("fHistdPhidEtaMixK0s"))->Fill(spMix, 1./weight);

							}
							else
								((THnSparseF*)((TList*)fOutput->FindObject("AnalysisK0s"))->FindObject("fHistdPhidEtaMixK0s"))->Fill(spMix);
						}
						else if(assoc->WhichCandidate() == 1){
							if(fEffCorr){
								Double_t weight = fHistEffEtaPtLambda->Interpolate(assoc->Eta(), 
										assoc->Pt());
								if(weight == 0){
									continue;
								}
								((THnSparseF*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistdPhidEtaMixLambda"))->Fill(spMix, 1./weight);
							}
							else
								((THnSparseF*)((TList*)fOutput->FindObject("AnalysisLambda"))->FindObject("fHistdPhidEtaMixLambda"))->Fill(spMix);
						}
						else if(assoc->WhichCandidate() == 2){
							if(fEffCorr){
								
								Double_t weight = fHistEffEtaPtAntiLambda->Interpolate(assoc->Eta(), 
										assoc->Pt());
								if(weight == 0){
									continue;
								}
								((THnSparseF*)((TList*)fOutput->FindObject("AnalysisAntiLambda"))->FindObject("fHistdPhidEtaMixAntiLambda"))->Fill(spMix, 1./weight);}
							else{
								((THnSparseF*)((TList*)fOutput->FindObject("AnalysisAntiLambda"))->FindObject("fHistdPhidEtaMixAntiLambda"))->Fill(spMix);}
						}
					}


					if(fChXi){

						if(assoc->WhichCandidate() == 3){

							if(fEffCorr){

								Double_t weight = fHistEffEtaPtXiMinus->Interpolate(assoc->Eta(), 	assoc->Pt());
								if(weight == 0){
									continue;}
								((THnSparseF*)((TList*)fOutput->FindObject("AnalysisXiMinus"))->FindObject("fHistdPhidEtaMixXiMinus"))->Fill(spMix, 1./weight);}
							else
								((THnSparseF*)((TList*)fOutput->FindObject("AnalysisXiMinus"))->FindObject("fHistdPhidEtaMixXiMinus"))->Fill(spMix);




						}  
					}  
				}//end of loop on assoc particles
			}// end of loop of mixing events
		}//end if pool 
	}//build correlation

		if(fEventMixing){

	   TObjArray* tracksClone = (TObjArray*) assocParticles->Clone();
	   tracksClone->SetOwner(kTRUE);
	   pool->UpdatePool(tracksClone);
		}

	PostData(1, fOutput);
}

//___________________________________________
Bool_t AliAnalysisTaskLambdaK0s::IsGoodPrimaryTrack(const AliAODTrack *t)
{
	// Pseudorapidity cut
	if (TMath::Abs(t->Eta())>fEtaCut) return kFALSE;

	// Should correspond to set of cuts suitable for correlation analysis
	if(fUseHybridGlobalTrack){
		if(!t->IsHybridGlobalConstrainedGlobal()) return kFALSE;
	}
	else{
		if (!t->TestFilterBit(768)) return kFALSE;//768-Hybrid track 
	}

	//	if(t->GetTPCClusterInfo(2,1) < fTPCNcls) return kFALSE;

	return kTRUE;
}
//_____________________________________________
Bool_t AliAnalysisTaskLambdaK0s::IsGoodDaughterTrack(const AliAODTrack *t)
{
	// Pseudorapidity cut   
	if (TMath::Abs(t->Eta()) > fV0DaughterEtaCut) return kFALSE;

	//pt cut
	if(t->Pt() < fV0DaughterPtMinCut) return kFALSE;

	//Standard TPC only track
	//  if(!t->TestFilterBit(1<<0)) return kFALSE;//在v0粒子的框架中已经有这个设置了。

	// TPC refit
	if (!t->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;

	// Minimum number of clusters
	Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1);
	if (nCrossedRowsTPC < fTPCNcls) return kFALSE;

	//findable clusters
	Int_t findable=t->GetTPCNclsF();
	if (findable <= 0) return kFALSE;

	if (nCrossedRowsTPC/findable < 0.8) return kFALSE;

	return kTRUE;
}

//___________________________________________________________________________________________
Bool_t AliAnalysisTaskLambdaK0s::IsMCV0Primary(AliAODv0 *v0, 
		Int_t specie){
	Int_t motherLabel = IsMcV0(v0, specie);
	if(motherLabel == -1) return kFALSE;
	AliAODMCParticle * part
		=  dynamic_cast<AliAODMCParticle*>(fMCArray->At(TMath::Abs(motherLabel))) ;
	if(!part) return kFALSE;
	if(part->IsPhysicalPrimary()) return kTRUE;
	return kFALSE;
}

//_________________________________________________________________________________________
Bool_t AliAnalysisTaskLambdaK0s::IsMCV0FromXi(AliAODv0 *v0, Int_t specie){
	Int_t motherLabel = IsMcV0(v0, specie);
	if(motherLabel == -1) return kFALSE;

	AliAODMCParticle * part
		=  dynamic_cast<AliAODMCParticle*>(fMCArray->At(TMath::Abs(motherLabel))) ;

	if(!part) return kFALSE;
	if(part->IsPhysicalPrimary()) return kFALSE;

	Int_t grandMotherLabel = part->GetMother();
	if(grandMotherLabel == -1) return kFALSE;
	AliAODMCParticle * partXi
		=  dynamic_cast<AliAODMCParticle*>(fMCArray->At(TMath::Abs(grandMotherLabel))) ;
	if(!partXi) return kFALSE;

	switch(specie){
		case 1: 
			if(partXi->GetPdgCode() == 3322
					|| partXi->GetPdgCode() == 3312)
				return kTRUE;
			break;
		case 2:
			if(partXi->GetPdgCode() == -3322
					|| partXi->GetPdgCode() == -3312)
				return kTRUE;
			break;
		default:
			break;
	}

	return kFALSE;
}

//_________________________________________________________________________________________
Int_t AliAnalysisTaskLambdaK0s::IsMcV0(AliAODv0 *v0, 
		Int_t specie) const{
	//                                                                           
	// check if the passed V0 is associated to a MC one,                         
	//   and returns the corresponding geant label.                              
	// returns -1 if the V0 is fake (i.e. label<0).                              
	//       

	AliAODVertex *vtx=v0->GetSecondaryVtx();
	AliAODTrack *nTrack = (AliAODTrack*)vtx->GetDaughter(1);
	AliAODTrack *pTrack = (AliAODTrack*)vtx->GetDaughter(0);

	if (!nTrack || !pTrack) return -1 ;

	Int_t nlab  = nTrack->GetLabel() ;
	Int_t plab  = pTrack->GetLabel() ;

	if (nlab == -1 || plab == -1) return -1 ;

	return GetV0Label(nlab,plab, specie) ;
}

//_______________________________________________________________________________________
Int_t AliAnalysisTaskLambdaK0s::GetV0Label(Int_t lab1,
		Int_t lab2,
		Int_t specie) const{
	//                                                                            
	// returns the label of the V0, given the labels of the 2 daughter tracks    
	// returns -1 if the V0 is fake                                              
	//       

	AliAODMCParticle *mcPart1
		=  dynamic_cast<AliAODMCParticle*>(fMCArray->At(TMath::Abs(lab1))) ;
	AliAODMCParticle *mcPart2
		=  dynamic_cast<AliAODMCParticle*>(fMCArray->At(TMath::Abs(lab2))) ;

	Int_t part1MotherLab = mcPart1->GetMother();
	Int_t part2MotherLab = mcPart2->GetMother();

	if (part1MotherLab==-1 || part2MotherLab==-1) return -1 ;
	if (part1MotherLab != part2MotherLab )        return -1 ;

	AliAODMCParticle *mcMother
		=  dynamic_cast<AliAODMCParticle*>(fMCArray->At(TMath::Abs(part1MotherLab))) ;

	if(specie == 0){//K0S
		if(mcMother->GetPdgCode() == 310
				&& ( (mcPart1->GetPdgCode() == -211 && mcPart2->GetPdgCode() == 211)
					|| (mcPart1->GetPdgCode() == 211 && mcPart2->GetPdgCode() == -211)
				   )
		  )
			return part1MotherLab ;
	}
	else if(specie == 1){
		if(mcMother->GetPdgCode() == 3122
				&& ( (mcPart1->GetPdgCode() == -211
						&& mcPart2->GetPdgCode() == 2212)
					|| (mcPart1->GetPdgCode() == 2212
						&& mcPart2->GetPdgCode() == -211) )
		  )
			return part1MotherLab ;
	}
	else if(specie == 2){
		if ( mcMother->GetPdgCode() == -3122
				&& ( (mcPart1->GetPdgCode() == 211
						&& mcPart2->GetPdgCode() == -2212)
					|| (mcPart1->GetPdgCode() == -2212
						&& mcPart2->GetPdgCode() == 211) )
		   )
			return part1MotherLab;
	}

	return -1;
}
//______________________________________________  
Bool_t AliAnalysisTaskLambdaK0s::IsMCXiPrimary(AliAODcascade *xi){

	AliAODTrack *pTrkXi = dynamic_cast<AliAODTrack*>( xi->GetDaughter(0) );
	AliAODTrack *nTrkXi = dynamic_cast<AliAODTrack*>( xi->GetDaughter(1) );
	AliAODTrack *bTrkXi = dynamic_cast<AliAODTrack*>( xi->GetDecayVertexXi()->GetDaughter(0) );


	Int_t nlab = nTrkXi->GetLabel();
	Int_t plab = pTrkXi->GetLabel();
	Int_t blab = bTrkXi->GetLabel();

	if(nlab < 0 || plab < 0 || blab < 0) return kFALSE;


	AliAODMCParticle *mcPart1    =  dynamic_cast<AliAODMCParticle*>(fMCArray->At(nlab));


	AliAODMCParticle *mcPart2    =  dynamic_cast<AliAODMCParticle*>(fMCArray->At(plab));


	AliAODMCParticle *mcPart3    =  dynamic_cast<AliAODMCParticle*>(fMCArray->At(blab));


	Int_t part1MotherLab = mcPart1->GetMother();
	Int_t part2MotherLab = mcPart2->GetMother();
	Int_t part3MotherLab = mcPart3->GetMother();

	if (part1MotherLab==-1 || part2MotherLab==-1 || part3MotherLab ==-1) return kFALSE;
	if (part1MotherLab != part2MotherLab )        return kFALSE;


	AliAODMCParticle * mcV0  =  dynamic_cast<AliAODMCParticle*>(fMCArray->At(part1MotherLab));

	if(!mcV0) return kFALSE;
	if(TMath::Abs(mcV0->GetPdgCode()) != 3122) return kFALSE;

	if(mcV0->GetPdgCode() == 3122){ // Lambda 
		if((mcPart1->GetPdgCode() != -211 || mcPart2->GetPdgCode() != 2212)
				&& (mcPart1->GetPdgCode() != 2212 || mcPart2->GetPdgCode() != -211))
			return kFALSE;
	}
	else if(mcV0->GetPdgCode() == -3122){ // AntiLambda 
		if((mcPart1->GetPdgCode() != 211 || mcPart2->GetPdgCode() != -2212)
				&& (mcPart1->GetPdgCode() != -2212 || mcPart2->GetPdgCode() != 211))
			return kFALSE;
	}


	Int_t mcV0MotherLab = mcV0->GetMother();
	if(mcV0MotherLab != part3MotherLab) return kFALSE;

	AliAODMCParticle * mcXi =  dynamic_cast<AliAODMCParticle*>(fMCArray->At(mcV0MotherLab));


	if(!mcXi) return kFALSE; 
	// if(fSpecie == 0) {
	if (TMath::Abs(mcXi->GetPdgCode()) != 3312) return kFALSE;
	if (mcXi->GetPdgCode() == 3312) {
		if(mcV0->GetPdgCode() != 3122 || mcPart3->GetPdgCode() != -211) return kFALSE;
	}else if(mcXi->GetPdgCode() == -3312){
		if(mcV0->GetPdgCode() != -3122 || mcPart3->GetPdgCode() != 211) return kFALSE;
	}

	return kTRUE;

}


//______________________________________________
Bool_t AliAnalysisTaskLambdaK0s::IsGoodV0(AliAODv0* aodV0)
{
	if (!aodV0) {
		AliError(Form("ERROR: Could not retrieve aodV0"));
		return kFALSE;
	}

	// Offline reconstructed V0 only
	if (aodV0->GetOnFlyStatus()) return kFALSE;

	// Get daughters and check them     
	AliAODTrack *myTrackNegTest=dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(1));
	AliAODTrack *myTrackPosTest=dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(0));

	if (!myTrackPosTest || !myTrackNegTest) {
		Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
		return kFALSE;
	}

	// Track cuts for daugher tracks 
	if ( !(IsGoodDaughterTrack(myTrackPosTest)) || !(IsGoodDaughterTrack(myTrackNegTest)) )
		return kFALSE;

	// Unlike signs of daughters
	if (myTrackNegTest->Charge() == myTrackPosTest->Charge()) return kFALSE;

	// Decay length
	Double_t dDL = aodV0->DecayLengthV0( fBestPrimaryVtxPos );

	// DCA of daughter track to Primary Vertex
	Float_t xyn = aodV0->DcaNegToPrimVertex();
	Float_t xyp = aodV0->DcaPosToPrimVertex();

	// Cosinus of pointing angle      
	Double_t dCPA = aodV0->CosPointingAngle(fBestPrimaryVtxPos);

	// DCA of daughter tracks 
	Double_t dDCA = aodV0->DcaV0Daughters();

	// DCA or V0 to prim. vertex
	Double_t dDCAV0PV = aodV0->DcaV0ToPrimVertex();

	// 
	Double_t dQT=aodV0->PtArmV0();
	Double_t dALPHA=aodV0->AlphaV0(); 
	Double_t dPT=aodV0->Pt();
	// Double_t dETA=aodV0->Eta();

	((TH2F*)((TList*)fOutput->FindObject("V0Candidates"))->FindObject("BefDL"))  ->Fill(dDL,dPT);
	((TH2F*)((TList*)fOutput->FindObject("V0Candidates"))->FindObject("BefDCA")) ->Fill(dDCA,dPT);
	((TH2F*)((TList*)fOutput->FindObject("V0Candidates"))->FindObject("BefDCA2PV"))->Fill(dDCAV0PV, dPT);
	((TH2F*)((TList*)fOutput->FindObject("V0Candidates"))->FindObject("BefCTP")) ->Fill(dCPA,dPT);
	((TH2F*)((TList*)fOutput->FindObject("V0Candidates"))->FindObject("BefD0"))  ->Fill(xyn,dPT);
	((TH2F*)((TList*)fOutput->FindObject("V0Candidates"))->FindObject("BefD0"))  ->Fill(xyp,dPT);
	((TH3F*)((TList*)fOutput->FindObject("V0Candidates"))->FindObject("BefAP"))  ->Fill(dALPHA,dQT,dPT);

	// Fiducial volume 
	Double_t xyz[3]; aodV0->GetSecondaryVtx(xyz);
	Double_t r2=xyz[0]*xyz[0] + xyz[1]*xyz[1];

	// cut on DCA of daughter track to Primary Vertex
	if (TMath::Abs(xyn)<fDCAToPrimVtx || TMath::Abs(xyp)<fDCAToPrimVtx) return kFALSE;

	// cut on DCA of daughter tracks 
	if (dDCA > fDCABetweenDaughters) return kFALSE;

	// cut on Cosinus of pointing angle
	if (dCPA < fCPA) return kFALSE;

	// Fiducial volume cut
	if (r2 < f2DFiducial*f2DFiducial) return kFALSE;
	//	if (r2 > 100*100) return kFALSE;

	// c*tau cut - in main V0 loop - depends on particle hypothesis
	// rapidity cut - in main V0 loop - depends on particle hypothesis

	// Daughter PID cut - in main V0 loop - depends on particle hypothesis

	((TH3F*)((TList*)fOutput->FindObject("V0Candidates"))->FindObject("AftAP"))  ->Fill(dALPHA,dQT,dPT);
	((TH2F*)((TList*)fOutput->FindObject("V0Candidates"))->FindObject("AftDCA2PV"))->Fill(dDCAV0PV, dPT);

	return kTRUE;
}


