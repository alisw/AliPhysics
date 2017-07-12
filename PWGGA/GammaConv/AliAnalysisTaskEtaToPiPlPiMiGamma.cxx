/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock                                                *
 * Version 1                                                              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Analysis task for eta->pi+ +pi- gamma (pion Dalitz decay)

#include <vector>

#include "TParticle.h"
#include "TPDGCode.h"
#include "TMCProcess.h"
#include "TDatabasePDG.h"
#include "TList.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "THnSparse.h"
#include "TH2F.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliPID.h"
#include "AliLog.h"
#include "AliESDtrackCuts.h"
#include "AliESDpidCuts.h"
#include "AliMCEvent.h"
#include "AliESDv0.h"
#include "AliESDEvent.h"
#include "AliESDpid.h"
#include "AliKFParticle.h"
#include "AliMCEventHandler.h"
#include "AliKFVertex.h"
#include "AliTriggerAnalysis.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "AliAnalysisTaskEtaToPiPlPiMiGamma.h"
#include <vector>


ClassImp( AliAnalysisTaskEtaToPiPlPiMiGamma )

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskEtaToPiPlPiMiGamma::AliAnalysisTaskEtaToPiPlPiMiGamma():
	fV0Reader(NULL),
    fV0ReaderName("V0ReaderV1"),
	fPionSelector(NULL),
	fBGHandler(NULL),
	fESDEvent(NULL),
	fMCEvent(NULL),
	fCutFolder(NULL),
	fESDList(NULL),
	fBackList(NULL),
	fMotherList(NULL),
	fTrueList(NULL),
	fMCList(NULL),
	fOutputContainer(0),
	fReaderGammas(NULL),
	fSelectorNegPionIndex(0),
	fSelectorPosPionIndex(0),
	fGoodGammas(NULL),
	fGoodVirtualParticles(NULL),
	fEventCutArray(NULL),
	fGammaCutArray(NULL),
	fPionCutArray(NULL),
	fMesonCutArray(NULL),
	fEventCuts(NULL),
	fConversionCuts(NULL),
	fHistoConvGammaPt(NULL),
	fHistoConvGammaEta(NULL),
	fHistoNegPionPt(NULL),
	fHistoPosPionPt(NULL),
	fHistoNegPionPhi(NULL),
	fHistoPosPionPhi(NULL),
	fHistoNegPionEta(NULL),
	fHistoPosPionEta(NULL),
	fHistoNegPionClsTPC(NULL),
	fHistoPosPionClsTPC(NULL),
	fHistoPionDCAxy(NULL),
	fHistoPionDCAz(NULL),
	fHistoPionTPCdEdxNSigma(NULL),
	fHistoPionTPCdEdx(NULL),
	fHistoPionPionInvMassPt(NULL),
	fHistoMotherInvMassPt(NULL),
	fTHnSparseMotherInvMassPtZM(NULL),
	fHistoMotherBackInvMassPt(NULL),
	fTHnSparseMotherBackInvMassPtZM(NULL),
	fHistoMCAllGammaPt(NULL),
	fHistoMCConvGammaPt(NULL),
	fHistoMCAllPosPionsPt(NULL),
	fHistoMCAllNegPionsPt(NULL),
	fHistoMCGammaFromEtaPt(NULL),
	fHistoMCPosPionsFromEtaPt(NULL),
	fHistoMCNegPionsFromEtaPt(NULL),
	fHistoMCEtaPiPlPiMiGammaPt(NULL),
	fHistoMCEtaGGPt(NULL),
	fHistoMCEtaDalitzPt(NULL),
	fHistoMCEtaPiPlPiMiGammaInAccPt(NULL),
	fHistoTrueMotherPiPlPiMiGammaInvMassPt(NULL),
	fHistoTrueMotherGammaGammaInvMassPt(NULL),
	fHistoTrueMotherDalitzInvMassPt(NULL),
	fHistoTrueConvGammaPt(NULL),
	fHistoTrueConvGammaFromEtaPt(NULL),
	fHistoTruePosPionPt(NULL),
	fHistoTruePosPionFromEtaPt(NULL),
	fHistoTrueNegPionPt(NULL),
	fHistoTrueNegPionFromEtaPt(NULL),
	fHistoTruePionPionInvMassPt(NULL),
	fHistoTruePionPionFromEtaInvMassPt(NULL),
	fHistoDoubleCountTrueEtaInvMassPt(NULL),
	fHistoDoubleCountTrueConvGammaRPt(NULL),
	fVectorDoubleCountTrueEtas(0),
	fVectorDoubleCountTrueConvGammas(0),
	fHistoNEvents(NULL),
	fHistoNGoodESDTracks(NULL),
	fProfileEtaShift(NULL),
	fHistoSPDClusterTrackletBackground(NULL),
	fRandom(0),
	fnCuts(0),
	fiCut(0),
	fNumberOfESDTracks(0),
	fMoveParticleAccordingToVertex(kFALSE),
	fIsHeavyIon(kFALSE),
	fDoMesonAnalysis(kTRUE),
	fDoMesonQA(kFALSE),
	fIsFromMBHeader(kTRUE),
	fIsMC(kFALSE),
	fIsGammaEtaCand(kFALSE)
{

}

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskEtaToPiPlPiMiGamma::AliAnalysisTaskEtaToPiPlPiMiGamma( const char* name ):
	AliAnalysisTaskSE(name),
	fV0Reader(NULL),
    fV0ReaderName("V0ReaderV1"),
	fPionSelector(NULL),
	fBGHandler(NULL),
	fESDEvent(NULL),
	fMCEvent(NULL),
	fCutFolder(NULL),
	fESDList(NULL),
	fBackList(NULL),
	fMotherList(NULL),
	fTrueList(NULL),
	fMCList(NULL),
	fOutputContainer(0),
	fReaderGammas(NULL),
	fSelectorNegPionIndex(0),
	fSelectorPosPionIndex(0),
	fGoodGammas(NULL),
	fGoodVirtualParticles(NULL),
	fEventCutArray(NULL),
	fGammaCutArray(NULL),
	fPionCutArray(NULL),
	fMesonCutArray(NULL),
	fEventCuts(NULL),
	fConversionCuts(NULL),
	fHistoConvGammaPt(NULL),
	fHistoConvGammaEta(NULL),
	fHistoNegPionPt(NULL),
	fHistoPosPionPt(NULL),
	fHistoNegPionPhi(NULL),
	fHistoPosPionPhi(NULL),
	fHistoNegPionEta(NULL),
	fHistoPosPionEta(NULL),
	fHistoNegPionClsTPC(NULL),
	fHistoPosPionClsTPC(NULL),
	fHistoPionDCAxy(NULL),
	fHistoPionDCAz(NULL),
	fHistoPionTPCdEdxNSigma(NULL),
	fHistoPionTPCdEdx(NULL),
	fHistoPionPionInvMassPt(NULL),
	fHistoMotherInvMassPt(NULL),
	fTHnSparseMotherInvMassPtZM(NULL),
	fHistoMotherBackInvMassPt(NULL),
	fTHnSparseMotherBackInvMassPtZM(NULL),
	fHistoMCAllGammaPt(NULL),
	fHistoMCConvGammaPt(NULL),
	fHistoMCAllPosPionsPt(NULL),
	fHistoMCAllNegPionsPt(NULL),
	fHistoMCGammaFromEtaPt(NULL),
	fHistoMCPosPionsFromEtaPt(NULL),
	fHistoMCNegPionsFromEtaPt(NULL),
	fHistoMCEtaPiPlPiMiGammaPt(NULL),
	fHistoMCEtaGGPt(NULL),
	fHistoMCEtaDalitzPt(NULL),
	fHistoMCEtaPiPlPiMiGammaInAccPt(NULL),
	fHistoTrueMotherPiPlPiMiGammaInvMassPt(NULL),
	fHistoTrueMotherGammaGammaInvMassPt(NULL),
	fHistoTrueMotherDalitzInvMassPt(NULL),
	fHistoTrueConvGammaPt(NULL),
	fHistoTrueConvGammaFromEtaPt(NULL),
	fHistoTruePosPionPt(NULL),
	fHistoTruePosPionFromEtaPt(NULL),
	fHistoTrueNegPionPt(NULL),
	fHistoTrueNegPionFromEtaPt(NULL),
	fHistoTruePionPionInvMassPt(NULL),
	fHistoTruePionPionFromEtaInvMassPt(NULL),
	fHistoDoubleCountTrueEtaInvMassPt(NULL),
	fHistoDoubleCountTrueConvGammaRPt(NULL),
	fVectorDoubleCountTrueEtas(0),
	fVectorDoubleCountTrueConvGammas(0),
	fHistoNEvents(NULL),
	fHistoNGoodESDTracks(NULL),
	fProfileEtaShift(NULL),
	fHistoSPDClusterTrackletBackground(NULL),
	fRandom(0),
	fnCuts(0),
	fiCut(0),
	fNumberOfESDTracks(0),
	fMoveParticleAccordingToVertex(kFALSE),
	fIsHeavyIon(kFALSE),
	fDoMesonAnalysis(kTRUE),
	fDoMesonQA(kFALSE),
	fIsFromMBHeader(kTRUE),
	fIsMC(kFALSE),
	fIsGammaEtaCand(kFALSE)
{
   DefineOutput(1, TList::Class());
}

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskEtaToPiPlPiMiGamma::~AliAnalysisTaskEtaToPiPlPiMiGamma()
{
	//
	// virtual destructor
	//
	cout<<"Destructor"<<endl;
	if(fGoodGammas){
		delete fGoodGammas;
		fGoodGammas = 0x0;
	}
	if(fGoodVirtualParticles){
		delete fGoodVirtualParticles;
		fGoodGammas = 0x0;
	}
	if(fBGHandler){
		delete[] fBGHandler;
		fBGHandler = 0x0;
	}
}
//___________________________________________________________
void AliAnalysisTaskEtaToPiPlPiMiGamma::InitBack(){

	const Int_t nDim = 4;
	Int_t nBins[nDim] = {450,250,7,4};
	Double_t xMin[nDim] = {0.3,0, 0,0};
	Double_t xMax[nDim] = {0.75,25,7,4};
	
	fTHnSparseMotherInvMassPtZM = new THnSparseF*[fnCuts];
	fTHnSparseMotherBackInvMassPtZM = new THnSparseF*[fnCuts];

	fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		
		TString cutstringEvent		= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
		TString cutstringPion		= ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutNumber();
		TString cutstringMeson		= ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
		TString cutstringGamma		= ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutNumber();
		
		Int_t collisionSystem = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(0,1));
		Int_t centMin = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(1,1));
		Int_t centMax = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(2,1));
		
		if(collisionSystem == 1 || collisionSystem == 2 ||
			collisionSystem == 5 || collisionSystem == 8 ||
			collisionSystem == 9){
			centMin = centMin*10;
			centMax = centMax*10; 
		}
		else if(collisionSystem == 3 || collisionSystem == 6){
			centMin = centMin*5;
			centMax = centMax*5;
		}
		else if(collisionSystem == 4 || collisionSystem == 7){
			centMin = ((centMin*5)+45);
			centMax = ((centMax*5)+45);
		}


		fBackList[iCut] = new TList();
		fBackList[iCut]->SetName(Form("%s_%s_%s_%s Back histograms",cutstringEvent.Data(),cutstringGamma.Data(),cutstringPion.Data(),cutstringMeson.Data()));
		fBackList[iCut]->SetOwner(kTRUE);
		fCutFolder[iCut]->Add(fBackList[iCut]);

		fTHnSparseMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_m","Back_Back_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
		fBackList[iCut]->Add(fTHnSparseMotherBackInvMassPtZM[iCut]);

		fMotherList[iCut] = new TList();
		fMotherList[iCut]->SetName(Form("%s_%s_%s_%s Mother histograms",cutstringEvent.Data(),cutstringGamma.Data(),cutstringPion.Data(),cutstringMeson.Data()));
		fMotherList[iCut]->SetOwner(kTRUE);
		fCutFolder[iCut]->Add(fMotherList[iCut]);

		fTHnSparseMotherInvMassPtZM[iCut] = new THnSparseF("Back_Mother_InvMass_Pt_z_m","Back_Mother_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
		fMotherList[iCut]->Add(fTHnSparseMotherInvMassPtZM[iCut]);

		
		fBGHandler[iCut] = new AliGammaConversionAODBGHandler(	collisionSystem,centMin,centMax,
																((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfBGEvents(),
																((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity(),
																0,8,5);
		
	}
}

//______________________________________________________________________
void AliAnalysisTaskEtaToPiPlPiMiGamma::UserCreateOutputObjects()
{
	//
	// Create ouput objects
	//

	// Create the output container
	if(fOutputContainer != NULL){
		delete fOutputContainer;
		fOutputContainer = NULL;
	}
	if(fOutputContainer == NULL){
		fOutputContainer = new TList();
		fOutputContainer->SetOwner(kTRUE);
	}

	fGoodGammas = new TList();

	fGoodVirtualParticles = new TList();
	fGoodVirtualParticles->SetOwner(kTRUE);

	fCutFolder				= new TList*[fnCuts];
	fESDList				= new TList*[fnCuts];
	fBackList				= new TList*[fnCuts];
	fMotherList 			= new TList*[fnCuts];
	fHistoNEvents			= new TH1I*[fnCuts];
	fHistoNGoodESDTracks	= new TH1I*[fnCuts];
	fProfileEtaShift		= new TProfile*[fnCuts];
	fHistoSPDClusterTrackletBackground = new TH2F*[fnCuts];
	fHistoConvGammaPt		= new TH1F*[fnCuts];
	fHistoConvGammaEta		= new TH1F*[fnCuts];
	fHistoNegPionPt			= new TH1F*[fnCuts];
	fHistoPosPionPt			= new TH1F*[fnCuts];
	fHistoNegPionPhi		= new TH1F*[fnCuts];
	fHistoPosPionPhi		= new TH1F*[fnCuts];
	
	if( fDoMesonQA ) {			
		fHistoNegPionEta		= new TH1F*[fnCuts];
		fHistoPosPionEta		= new TH1F*[fnCuts];
		fHistoNegPionClsTPC		= new TH2F*[fnCuts];
		fHistoPosPionClsTPC		= new TH2F*[fnCuts];
		fHistoPionDCAxy			= new TH2F*[fnCuts];
		fHistoPionDCAz			= new TH2F*[fnCuts];
		fHistoPionTPCdEdxNSigma	= new TH2F*[fnCuts];
		fHistoPionTPCdEdx		= new TH2F*[fnCuts];
		fHistoPionPionInvMassPt	= new TH2F*[fnCuts];

	}
	
	fHistoMotherInvMassPt		= new TH2F*[fnCuts];
	fHistoMotherBackInvMassPt	= new TH2F*[fnCuts];

	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		TString cutstringEvent 	= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
		TString cutstringPion 	= ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutNumber();
		TString cutstringMeson	= ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
		TString cutstringGamma 	= ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutNumber();

		fCutFolder[iCut] = new TList();
		fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s_%s",cutstringEvent.Data(),cutstringGamma.Data(),cutstringPion.Data(),cutstringMeson.Data()));
		fCutFolder[iCut]->SetOwner(kTRUE);
		fOutputContainer->Add(fCutFolder[iCut]);

		fESDList[iCut] = new TList();
		fESDList[iCut]->SetName(Form("%s_%s_%s_%s ESD histograms",cutstringEvent.Data(),cutstringGamma.Data(),cutstringPion.Data(),cutstringMeson.Data()));
		fESDList[iCut]->SetOwner(kTRUE);

		fHistoNEvents[iCut] = new TH1I("NEvents","NEvents",14,-0.5,13.5);
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problems");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fESDList[iCut]->Add(fHistoNEvents[iCut]);

		if(fIsHeavyIon) fHistoNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",3000,0,3000);
			else fHistoNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",200,0,200);
		fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);

		fProfileEtaShift[iCut] = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
		fESDList[iCut]->Add(fProfileEtaShift[iCut]);
		fHistoSPDClusterTrackletBackground[iCut] = new TH2F("SPD tracklets vs SPD clusters","SPD tracklets vs SPD clusters",100,0,200,250,0,1000);
		fESDList[iCut]->Add(fHistoSPDClusterTrackletBackground[iCut]);
		fHistoConvGammaPt[iCut] = new TH1F("ESD_ConvGamma_Pt","ESD_ConvGamma_Pt",250,0,25);
		fESDList[iCut]->Add(fHistoConvGammaPt[iCut]);
		fHistoConvGammaEta[iCut] = new TH1F("ESD_ConvGamma_Eta","ESD_ConvGamma_Eta",600,-1.5,1.5);
		fESDList[iCut]->Add(fHistoConvGammaEta[iCut]);
		fHistoNegPionPt[iCut] = new TH1F("ESD_PrimaryNegPions_Pt","ESD_PrimaryNegPions_Pt",1000,0,25);
		fESDList[iCut]->Add(fHistoNegPionPt[iCut]);
		fHistoPosPionPt[iCut] = new TH1F("ESD_PrimaryPosPions_Pt","ESD_PrimaryPosPions_Pt",1000,0,25);
		fESDList[iCut]->Add(fHistoPosPionPt[iCut]);
		fHistoNegPionPhi[iCut] = new TH1F("ESD_PrimaryNegPions_Phi","ESD_PrimaryNegPions_Phi",360,0,2*TMath::Pi());
		fESDList[iCut]->Add(fHistoNegPionPhi[iCut]);
		fHistoPosPionPhi[iCut] = new TH1F("ESD_PrimaryPosPions_Phi","ESD_PrimaryPosPions_Phi",360,0,2*TMath::Pi());
		fESDList[iCut]->Add(fHistoPosPionPhi[iCut]);
		
		if ( fDoMesonQA ) {
			fHistoNegPionEta[iCut] = new TH1F("ESD_PrimaryNegPions_Eta","ESD_PrimaryNegPions_Eta",600,-1.5,1.5);
			fESDList[iCut]->Add(fHistoNegPionEta[iCut]);
			fHistoPosPionEta[iCut] = new TH1F("ESD_PrimaryPosPions_Eta","ESD_PrimaryPosPions_Eta",600,-1.5,1.5);
			fESDList[iCut]->Add(fHistoPosPionEta[iCut]);
			fHistoNegPionClsTPC[iCut]  = new TH2F("ESD_PrimaryNegPions_ClsTPC","ESD_PrimaryNegPions_ClsTPC",100,0,1,400,0.,10.);
			fESDList[iCut]->Add(fHistoNegPionClsTPC[iCut]);
			fHistoPosPionClsTPC[iCut]  = new TH2F("ESD_PrimaryPosPions_ClsTPC","ESD_PrimaryPosPions_ClsTPC",100,0,1,400,0.,10.);
			fESDList[iCut]->Add(fHistoPosPionClsTPC[iCut]);
			fHistoPionDCAxy[iCut] = new TH2F("ESD_PrimaryPions_DCAxy","ESD_PrimaryPions_DCAxy",800,-4.0,4.0,400,0.,10.);
			fESDList[iCut]->Add(fHistoPionDCAxy[iCut]);
			fHistoPionDCAz[iCut]  = new TH2F("ESD_PrimaryPions_DCAz","ESD_PrimaryPions_DCAz",800,-4.0,4.0,400,0.,10.);
			fESDList[iCut]->Add(fHistoPionDCAz[iCut]);
			fHistoPionTPCdEdxNSigma[iCut] = new TH2F("ESD_PrimaryPions_TPCdEdx","ESD_PrimaryPions_TPCdEdx",150,0.05,20,400,-10,10);
			fESDList[iCut]->Add(fHistoPionTPCdEdxNSigma[iCut]);
			fHistoPionTPCdEdx[iCut] =new TH2F("ESD_PrimaryPions_TPCdEdxSignal","ESD_PrimaryPions_TPCdEdxSignal" ,150,0.05,20.0,800,0.0,200);
			fESDList[iCut]->Add(fHistoPionTPCdEdx[iCut]);  			
			fHistoPionPionInvMassPt[iCut] = new TH2F("ESD_PiPlusPiNeg_InvMassPt","ESD_PiPlusPiNeg_InvMassPt",2000,0.,2.,200,0.,20.);
			fESDList[iCut]->Add(fHistoPionPionInvMassPt[iCut]);
		}

		fHistoMotherInvMassPt[iCut] = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",450,0.3,0.75,250,0,25);
		fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);
		fHistoMotherBackInvMassPt[iCut] = new TH2F("ESD_Background_InvMass_Pt","ESD_Background_InvMass_Pt",450,0.3,0.75,250,0,25);
		fESDList[iCut]->Add(fHistoMotherBackInvMassPt[iCut]);

		if ( fDoMesonQA ) {
			TAxis *AxisAfter = fHistoPionTPCdEdxNSigma[iCut]->GetXaxis(); 
			Int_t bins = AxisAfter->GetNbins();
			Double_t from = AxisAfter->GetXmin();
			Double_t to = AxisAfter->GetXmax();
			Double_t *newBins = new Double_t[bins+1];
			newBins[0] = from;
			Double_t factor = TMath::Power(to/from, 1./bins);
			for(Int_t i=1; i<=bins; ++i) newBins[i] = factor * newBins[i-1];

			AxisAfter->Set(bins, newBins);
			AxisAfter = fHistoPionTPCdEdx[iCut]->GetXaxis(); 
			AxisAfter->Set(bins, newBins);

			delete [] newBins;		
		}

		fCutFolder[iCut]->Add(fESDList[iCut]);

	}

	if( fIsMC ){
		// MC Histogramms
		fMCList = new TList*[fnCuts];
		// True Histogramms
		fTrueList = new TList*[fnCuts];
		fHistoTrueConvGammaPt = new TH1F*[fnCuts];
		fHistoDoubleCountTrueConvGammaRPt = new TH2F*[fnCuts];
		fHistoTrueConvGammaFromEtaPt = new TH1F*[fnCuts];
		fHistoTruePosPionPt  = new TH1F*[fnCuts];
		fHistoTrueNegPionPt  = new TH1F*[fnCuts];		
		fHistoTruePosPionFromEtaPt  = new TH1F*[fnCuts];
		fHistoTrueNegPionFromEtaPt  = new TH1F*[fnCuts];
		

		fHistoMCAllGammaPt  = new TH1F*[fnCuts];
		fHistoMCConvGammaPt = new TH1F*[fnCuts];
		fHistoMCAllPosPionsPt = new TH1F*[fnCuts];
		fHistoMCAllNegPionsPt = new TH1F*[fnCuts];
		fHistoMCGammaFromEtaPt  = new TH1F*[fnCuts];
		fHistoMCPosPionsFromEtaPt = new TH1F*[fnCuts];
		fHistoMCNegPionsFromEtaPt = new TH1F*[fnCuts];

// 		hMCPi0DalitzGammaPt    = new TH1F*[fnCuts];
// 		hMCPi0DalitzElectronPt = new TH1F*[fnCuts];
// 		hMCPi0DalitzPositronPt = new TH1F*[fnCuts];

		fHistoMCEtaPiPlPiMiGammaPt = new TH1F*[fnCuts];
		fHistoMCEtaDalitzPt = new TH1F*[fnCuts];
		fHistoMCEtaGGPt = new TH1F*[fnCuts];
		fHistoMCEtaPiPlPiMiGammaInAccPt = new TH1F*[fnCuts];
// 				
		fHistoDoubleCountTrueEtaInvMassPt = new TH2F*[fnCuts];
		fHistoTrueMotherPiPlPiMiGammaInvMassPt = new TH2F*[fnCuts];
		fHistoTrueMotherDalitzInvMassPt = new TH2F*[fnCuts];
		fHistoTrueMotherGammaGammaInvMassPt = new TH2F*[fnCuts];

		if (fDoMesonQA){
			fHistoTruePionPionInvMassPt = 			new TH2F*[fnCuts];
			fHistoTruePionPionFromEtaInvMassPt = 	new TH2F*[fnCuts];
		}
		
		for(Int_t iCut = 0; iCut<fnCuts;iCut++){
			TString cutstringEvent 	= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
			TString cutstringPion 	= ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutNumber();
			TString cutstringMeson	= ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
			TString cutstringGamma	= ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutNumber();

			fMCList[iCut] = new TList();
			fMCList[iCut]->SetName(Form("%s_%s_%s_%s MC histograms",cutstringEvent.Data(),cutstringGamma.Data(),cutstringPion.Data(),cutstringMeson.Data()));
			fMCList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fMCList[iCut]);

			fHistoMCAllGammaPt[iCut] = new TH1F("MC_AllGamma_Pt","MC_AllGamma_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCAllGammaPt[iCut]);			
			fHistoMCConvGammaPt[iCut] = new TH1F("MC_ConvGamma_Pt","MC_ConvGamma_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCConvGammaPt[iCut]);						
			fHistoMCAllPosPionsPt[iCut] = new TH1F("MC_AllPosPions_Pt","MC_AllPosPions_Pt",1000,0,25);
			fMCList[iCut]->Add(fHistoMCAllPosPionsPt[iCut]);
			fHistoMCAllNegPionsPt[iCut] = new TH1F("MC_AllNegPions_Pt","MC_AllNegPions_Pt",1000,0,25);
			fMCList[iCut]->Add(fHistoMCAllNegPionsPt[iCut]);
			fHistoMCGammaFromEtaPt[iCut] = new TH1F("MC_GammaFromEta_Pt","MC_GammaFromEta_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCGammaFromEtaPt[iCut]);	
			fHistoMCPosPionsFromEtaPt[iCut] = new TH1F("MC_PosPionsFromEta_Pt","MC_PosPionsFromEta_Pt",1000,0,25);
			fMCList[iCut]->Add(fHistoMCPosPionsFromEtaPt[iCut]);
			fHistoMCNegPionsFromEtaPt[iCut] = new TH1F("MC_NegPionsFromEta_Pt","MC_NegPionsFromEta_Pt",1000,0,25);
			fMCList[iCut]->Add(fHistoMCNegPionsFromEtaPt[iCut]);		

			fHistoMCEtaPiPlPiMiGammaPt[iCut] = new TH1F("MC_Eta_Pt","MC_Eta_Pt",250,0,25);
			fHistoMCEtaPiPlPiMiGammaPt[iCut]->Sumw2();
			fMCList[iCut]->Add(fHistoMCEtaPiPlPiMiGammaPt[iCut]);

			fHistoMCEtaDalitzPt[iCut] = new TH1F("MC_Eta_Dalitz_Pt","MC_Eta_Pt",250,0,25);
			fHistoMCEtaDalitzPt[iCut]->Sumw2();
			fMCList[iCut]->Add(fHistoMCEtaDalitzPt[iCut]);

			fHistoMCEtaGGPt[iCut] = new TH1F("MC_Eta_GG_Pt","MC_Eta_GG_Pt",250,0,25);
			fHistoMCEtaGGPt[iCut]->Sumw2();
			fMCList[iCut]->Add(fHistoMCEtaGGPt[iCut]);
			
			fHistoMCEtaPiPlPiMiGammaInAccPt[iCut] = new TH1F("MC_EtaInAcc_Pt","MC_EtaInAcc_Pt",250,0,25);
			fHistoMCEtaPiPlPiMiGammaInAccPt[iCut]->Sumw2();
			fMCList[iCut]->Add(fHistoMCEtaPiPlPiMiGammaInAccPt[iCut]);

			fTrueList[iCut] = new TList();
			fTrueList[iCut]->SetName(Form("%s_%s_%s_%s True histograms",cutstringEvent.Data(),cutstringGamma.Data(),cutstringPion.Data(),cutstringMeson.Data()));
			fTrueList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fTrueList[iCut]);

			fHistoTrueConvGammaPt[iCut] = new TH1F("ESD_TrueConvGamma_Pt","ESD_TrueConvGamma_Pt",250,0,25);
			fTrueList[iCut]->Add(fHistoTrueConvGammaPt[iCut]);
			fHistoDoubleCountTrueConvGammaRPt[iCut] = new TH2F("ESD_TrueDoubleCountConvGamma_R_Pt","ESD_TrueDoubleCountConvGamma_R_Pt",800,0,200,300,0,30);
			fTrueList[iCut]->Add(fHistoDoubleCountTrueConvGammaRPt[iCut]);
			fHistoTrueConvGammaFromEtaPt[iCut] = new TH1F("ESD_TrueConvGammaFromEta_Pt","ESD_TrueConvGammaFromEta_Pt",250,0,25);
			fTrueList[iCut]->Add(fHistoTrueConvGammaFromEtaPt[iCut]);
		
			fHistoTruePosPionPt[iCut] = new TH1F("ESD_TruePosPion_Pt","ESD_TruePosPion_Pt",1000,0,25);
			fTrueList[iCut]->Add(fHistoTruePosPionPt[iCut]);
			fHistoTrueNegPionPt[iCut] = new TH1F("ESD_TrueNegPion_Pt","ESD_TrueNegPion_Pt",1000,0,25);
			fTrueList[iCut]->Add(fHistoTrueNegPionPt[iCut]);	

			fHistoTrueNegPionFromEtaPt[iCut] = new TH1F("ESD_TrueNegPionFromEta_Pt","ESD_TrueNegPionFromEta_Pt",1000,0,25);
			fTrueList[iCut]->Add(fHistoTrueNegPionFromEtaPt[iCut]);
			fHistoTruePosPionFromEtaPt[iCut] = new TH1F("ESD_TruePosPionFromEta_Pt","ESD_TruePosPionFromEta_Pt",1000,0,25);
			fTrueList[iCut]->Add(fHistoTruePosPionFromEtaPt[iCut]);

			fHistoDoubleCountTrueEtaInvMassPt[iCut] = new TH2F("ESD_TrueDoubleCountEta_InvMass_Pt","ESD_TrueDoubleCountEta_InvMass_Pt",800,0,0.8,300,0,30);
			fTrueList[iCut]->Add(fHistoDoubleCountTrueEtaInvMassPt[iCut]);

			fHistoTrueMotherPiPlPiMiGammaInvMassPt[iCut] = new TH2F("ESD_TrueMotherPiPlPiMiGamma_InvMass_Pt","ESD_TrueMotherPiPlPiMiGamma_InvMass_Pt",450,0.3,0.75,250,0,25);
			fHistoTrueMotherPiPlPiMiGammaInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiGammaInvMassPt[iCut]);
		
			fHistoTrueMotherGammaGammaInvMassPt[iCut] = new TH2F("ESD_TrueMotherGG_InvMass_Pt","ESD_TrueMotherGG_InvMass_Pt",450,0.3,0.75,250,0,25);
			fHistoTrueMotherGammaGammaInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(fHistoTrueMotherGammaGammaInvMassPt[iCut]);
			
			fHistoTrueMotherDalitzInvMassPt[iCut] = new TH2F("ESD_TrueMotherDalitz_InvMass_Pt","ESD_TrueMotherDalitz_InvMass_Pt",450,0.3,0.75,250,0,25);
			fHistoTrueMotherDalitzInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(fHistoTrueMotherDalitzInvMassPt[iCut]);

			if (fDoMesonQA){
				fHistoTruePionPionInvMassPt[iCut] = new TH2F("ESD_TruePiPlusPiNeg_InvMassPt","ESD_TruePiPlusPiNeg_InvMassPt",2000,0.,2.,200,0.,20.);
				fTrueList[iCut]->Add(fHistoTruePionPionInvMassPt[iCut]);
				fHistoTruePionPionFromEtaInvMassPt[iCut] = new TH2F("ESD_TruePiPlusPiNegFromEta_InvMassPt","ESD_TruePiPlusPiNegFromEta_InvMassPt",2000,0.,2.,200,0.,20.);
				fTrueList[iCut]->Add(fHistoTruePionPionFromEtaInvMassPt[iCut]);
			}
		}
	}

	fVectorDoubleCountTrueEtas.clear();
	fVectorDoubleCountTrueConvGammas.clear();

	InitBack(); // Init Background Handler

    fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
	if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader
		
	if(fV0Reader)
		if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
			if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
				fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
		
		
		
	fPionSelector=(AliPrimaryPionSelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("PionSelector");
	if(!fPionSelector){printf("Error: No PionSelector");return;} // GetV0Reader
		
	if( fPionSelector ){
		if ( ((AliPrimaryPionCuts*)fPionSelector->GetPrimaryPionCuts())->GetCutHistograms() ){
			fOutputContainer->Add( ((AliPrimaryPionCuts*)fPionSelector->GetPrimaryPionCuts())->GetCutHistograms() );
		}
	}  

	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		if( fPionCutArray ){
			if( ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutHistograms() ) {
				fCutFolder[iCut]->Add( ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutHistograms() );
			}
		}
		if( fMesonCutArray  ) {
			if( ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms() ) {
				fCutFolder[iCut]->Add( ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
			}
		}
		if( fGammaCutArray ) {
			if( ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutHistograms() ) {
				fCutFolder[iCut]->Add( ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutHistograms()  );
			}
		}
	}

	PostData(1, fOutputContainer);

}

//______________________________________________________________________
void AliAnalysisTaskEtaToPiPlPiMiGamma::UserExec(Option_t *){

	//
	// Execute analysis for current event
	//

  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
	if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader

	Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
	if(InputEvent()->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
	if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete
		for(Int_t iCut = 0; iCut<fnCuts; iCut++){
			fHistoNEvents[iCut]->Fill(eventQuality);
		}
		return;
	}

	fPionSelector=(AliPrimaryPionSelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("PionSelector");
	if(!fPionSelector){printf("Error: No PionSelector");return;} // GetV0Reader

	if(fIsMC) fMCEvent     =  MCEvent();
	fESDEvent        = (AliESDEvent*)InputEvent();
	fReaderGammas    = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut
	fSelectorNegPionIndex = fPionSelector->GetReconstructedNegPionIndex(); // Electrons from default Cut
	fSelectorPosPionIndex = fPionSelector->GetReconstructedPosPionIndex(); // Positrons from default Cut

	fNumberOfESDTracks = fV0Reader->GetNumberOfPrimaryTracks();
	//AddTaskContainers(); //Add conatiner

	for(Int_t iCut = 0; iCut<fnCuts; iCut++){
		fiCut = iCut;
		Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);
		
		if(eventNotAccepted){
			// 			cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
			fHistoNEvents[iCut]->Fill(eventNotAccepted); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
			continue;
		}

		if(eventQuality != 0){// Event Not Accepted
			// 			cout << "event rejected due to: " <<eventQuality << endl;
			fHistoNEvents[iCut]->Fill(eventQuality);
			continue;
		}

		fHistoNEvents[iCut]->Fill(eventQuality);
		fHistoNGoodESDTracks[iCut]->Fill(fNumberOfESDTracks);
		fHistoSPDClusterTrackletBackground[iCut]->Fill(fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),(fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)));

		if(fMCEvent){ // Process MC Particle
			if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection() != 0){
				((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(), 
																						((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
																						fMCEvent);
			} 
			ProcessMCParticles();
		}

		fIsGammaEtaCand =kFALSE;
// 		cout << "new event" << endl;
		ProcessPhotonCandidates(); // Process this cuts gammas
		ProcessPionCandidates(); // Process this cuts gammas
		
	
		CalculateMesonCandidates();
		CalculateBackground();
		UpdateEventByEventData();

		fVectorDoubleCountTrueEtas.clear();
		fVectorDoubleCountTrueConvGammas.clear();
				
		fGoodGammas->Clear(); // delete this cuts good gammas
		fGoodVirtualParticles->Clear(); // delete this cuts good gammas
	}

	fSelectorNegPionIndex.clear();
	fSelectorPosPionIndex.clear();

	PostData( 1, fOutputContainer );
}

Bool_t AliAnalysisTaskEtaToPiPlPiMiGamma::Notify(){
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
     if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){        
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnum(fV0Reader->GetPeriodName());
    }  
   
		if( !((AliConvEventCuts*)fEventCutArray->At(iCut))->GetDoEtaShift() ){
			fProfileEtaShift[iCut]->Fill(0.,0.);
			continue; // No Eta Shift requested, continue
		}
		if( ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
			((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCorrectEtaShiftFromPeriod();
			((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once   
			((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->SetEtaShift( ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() );
			fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
			continue;
		} else {
			printf(" Eta t PiPlusPiMinus Gamma Task %s :: Eta Shift Manually Set to %f \n\n",
			(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber()).Data(),((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift());
			((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once   
			((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->SetEtaShift( ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() );
			fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
		}
	}
	return kTRUE;
}


void AliAnalysisTaskEtaToPiPlPiMiGamma::Terminate(const Option_t *){
///Grid
}

//________________________________________________________________________
void AliAnalysisTaskEtaToPiPlPiMiGamma::ProcessPhotonCandidates(){
	Int_t nV0 = 0;
	TList *GoodGammasStepOne = new TList();
	TList *GoodGammasStepTwo = new TList();
	// Loop over Photon Candidates allocated by ReaderV1
	
	for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
		AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
		if(!PhotonCandidate) continue;
		
		fIsFromMBHeader = kTRUE;
		
		if( fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0 ){		
			Int_t isPosFromMBHeader
                = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
			if(isPosFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
			Int_t isNegFromMBHeader
                = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent,fInputEvent);
			if(isNegFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
			if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
		}
		
		if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fESDEvent)) continue;

		if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut() &&
			!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){ // if no post reader loop is required add to events good gammas
			
			fGoodGammas->Add(PhotonCandidate);
		
			if(fIsFromMBHeader){
				fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
				fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
			}
		
			if(fMCEvent){
				ProcessTruePhotonCandidates(PhotonCandidate);
			}
		} else if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
			((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->FillElectonLabelArray(PhotonCandidate,nV0);
			nV0++;
			GoodGammasStepOne->Add(PhotonCandidate);
		} else if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut() &&
				((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
			GoodGammasStepTwo->Add(PhotonCandidate);
		}
	}
	
	
	if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut()){
		for(Int_t i = 0;i<GoodGammasStepOne->GetEntries();i++){
			AliAODConversionPhoton *PhotonCandidate= (AliAODConversionPhoton*) GoodGammasStepOne->At(i);
			if(!PhotonCandidate) continue;
			fIsFromMBHeader = kTRUE;
			if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
				Int_t isPosFromMBHeader
                = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent,fInputEvent);
				Int_t isNegFromMBHeader
                = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent,fInputEvent);
				if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
			}
			if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GoodGammasStepOne->GetEntries())) continue;
			if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
				fGoodGammas->Add(PhotonCandidate);
				if(fIsFromMBHeader){
					fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
					fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
				}
				if(fMCEvent){
					ProcessTruePhotonCandidates(PhotonCandidate);
				}
			}
			else GoodGammasStepTwo->Add(PhotonCandidate); // Close v0s cut enabled -> add to list two
		}
	}
	if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){
		for(Int_t i = 0;i<GoodGammasStepTwo->GetEntries();i++){
			AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GoodGammasStepTwo->At(i);
			if(!PhotonCandidate) continue;
			
			if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
				Int_t isPosFromMBHeader
                = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent,fInputEvent);
				Int_t isNegFromMBHeader
                = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent,fInputEvent);
				if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
			}
			
			if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GoodGammasStepTwo,i)) continue;
			fGoodGammas->Add(PhotonCandidate); // Add gamma to current cut TList
		
			if(fIsFromMBHeader){
				fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt()); // Differences to old V0Reader in p_t due to conversion KF->TLorentzVector
				fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
			}
		
			if(fMCEvent){
				ProcessTruePhotonCandidates(PhotonCandidate);
			}
		}
	}

	delete GoodGammasStepOne;
	GoodGammasStepOne = 0x0;
	delete GoodGammasStepTwo;
	GoodGammasStepTwo = 0x0;
}

//________________________________________________________________________
void AliAnalysisTaskEtaToPiPlPiMiGamma::ProcessTruePhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
	// Process True Photons
    TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(fMCEvent);
    TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCEvent);

	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();

	
	if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
	if(posDaughter->GetMother(0) != negDaughter->GetMother(0)){  // Not Same Mother == Combinatorial Bck
		return;
	}
	
	else if (posDaughter->GetMother(0) == -1){
		return;
	}
	
	if(TMath::Abs(posDaughter->GetPdgCode())!=11 || TMath::Abs(negDaughter->GetPdgCode())!=11) return; //One Particle is not electron
	if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge
	if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() !=5) return;// check if the daughters come from a conversion

    TParticle *Photon = TruePhotonCandidate->GetMCParticle(fMCEvent);
	if(Photon->GetPdgCode() != 22) return; // Mother is no Photon

	// True Photon

	if (CheckVectorForDoubleCount(fVectorDoubleCountTrueConvGammas,posDaughter->GetMother(0))) fHistoDoubleCountTrueConvGammaRPt[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),TruePhotonCandidate->Pt());
	
    Int_t labelGamma = TruePhotonCandidate->GetMCParticleLabel(fMCEvent);
    Bool_t gammaIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelGamma, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
	if( gammaIsPrimary ){
		if( fIsFromMBHeader ){
			fHistoTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
		}
	}
	
	if( IsEtaPiPlPiMiGammaDaughter(labelGamma) == kTRUE ) {
		if( gammaIsPrimary ) {
			if( fIsFromMBHeader ){
				fHistoTrueConvGammaFromEtaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				
// 				TParticle * gammaMC = (TParticle*)fMCEvent->Particle(labelGamma);
// 				Int_t gammaMotherLabel=gammaMC->GetFirstMother();
// 				for(Int_t index= ((TParticle*)fMCEvent->Particle(gammaMotherLabel))->GetFirstDaughter();index<= ((TParticle*)fMCEvent->Particle(gammaMotherLabel))->GetLastDaughter();index++){
// 					TParticle* temp = (TParticle*)fMCEvent->Particle( index );
// 					switch( temp->GetPdgCode() ) {
// 						case 211:
// 							cout << "pi- " << index << "\t" << temp->Pt() << "\t" << temp->Eta() << endl;
// 							break;
// 						case -211:
// 							cout << "pi+ " << index << "\t" << temp->Pt() << "\t" << temp->Eta() << endl;
// 							break;
// 						case ::kGamma:
// 							cout << "gamma " << index << "\t" << temp->Pt()<< "\t" << temp->Eta() << endl;
// 							break;
// 					}
// 				}						
				fIsGammaEtaCand = kTRUE;	
			}
		}
	} 
}

//________________________________________________________________________
void AliAnalysisTaskEtaToPiPlPiMiGamma::ProcessPionCandidates(){

	Double_t magField = fInputEvent->GetMagneticField();
	if( magField  < 0.0 ){
		magField =  1.0;
	} else {
		magField =  -1.0;
	}

	vector<Int_t> lGoodNegPionIndexPrev(0);
	vector<Int_t> lGoodPosPionIndexPrev(0);
	
    for(Long_t i = 0; i < fSelectorNegPionIndex.size(); i++){
		AliESDtrack* negPionCandidate = fESDEvent->GetTrack(fSelectorNegPionIndex[i]);
		if(! ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelected(negPionCandidate) ) continue;
		lGoodNegPionIndexPrev.push_back(   fSelectorNegPionIndex[i] );
		fHistoNegPionPt[fiCut]->Fill(negPionCandidate->Pt());
		fHistoNegPionPhi[fiCut]->Fill(negPionCandidate->Phi());
		if( fMCEvent ) {
			const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
			Double_t mcProdVtxX 	= primVtxMC->GetX();
			Double_t mcProdVtxY 	= primVtxMC->GetY();
			Double_t mcProdVtxZ 	= primVtxMC->GetZ();

			Int_t labelNegPion = TMath::Abs( negPionCandidate->GetLabel() );
            Bool_t negPionIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelNegPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
			if( negPionIsPrimary ){
                TParticle* negPion = fMCEvent->Particle(labelNegPion);
				if( negPion->GetPdgCode() ==  -211 ){
					fHistoTrueNegPionPt[fiCut]->Fill(negPionCandidate->Pt());    //primary negPion
					if( IsEtaPiPlPiMiGammaDaughter(labelNegPion) == kTRUE ) {
						fHistoTrueNegPionFromEtaPt[fiCut]->Fill(negPionCandidate->Pt());
// 							if (fIsGammaEtaCand) cout << "pi- rec" << labelNegPion << "\t" << negPionCandidate->Pt()<< endl;
					}	
				}
			}
		}
	}

    for(Long_t i = 0; i < fSelectorPosPionIndex.size(); i++){
		AliESDtrack* posPionCandidate = fESDEvent->GetTrack( fSelectorPosPionIndex[i] );
		if(! ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelected(posPionCandidate) ) continue;
		lGoodPosPionIndexPrev.push_back(   fSelectorPosPionIndex[i]  );
		fHistoPosPionPt[fiCut]->Fill( posPionCandidate->Pt() );
		fHistoPosPionPhi[fiCut]->Fill( posPionCandidate->Phi() );
		
		if( fMCEvent ) {
			const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
			Double_t mcProdVtxX 	= primVtxMC->GetX();
			Double_t mcProdVtxY 	= primVtxMC->GetY();
			Double_t mcProdVtxZ 	= primVtxMC->GetZ();

			Int_t labelPosPion = TMath::Abs( posPionCandidate->GetLabel() );
            Bool_t posPionIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelPosPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
			if( posPionIsPrimary ) {
                TParticle* posPion = fMCEvent->Particle(labelPosPion);
				if( posPion->GetPdgCode() ==  211 ){
					fHistoTruePosPionPt[fiCut]->Fill(posPionCandidate->Pt());
					if( IsEtaPiPlPiMiGammaDaughter(labelPosPion) == kTRUE ) {
						fHistoTruePosPionFromEtaPt[fiCut]->Fill(posPionCandidate->Pt());
// 							if (fIsGammaEtaCand) cout << "pi+ rec" << labelPosPion << "\t" << posPionCandidate->Pt()<< endl;
					}
				}
			}
		}
	}


    for(Long_t i = 0; i < lGoodNegPionIndexPrev.size(); i++){

		AliESDtrack *negPionCandidate = fESDEvent->GetTrack(lGoodNegPionIndexPrev[i]);
		AliKFParticle negPionCandidateKF( *negPionCandidate->GetConstrainedParam(), 211 );

        for(Long_t j = 0; j < lGoodPosPionIndexPrev.size(); j++){

			AliESDtrack *posPionCandidate = fESDEvent->GetTrack(lGoodPosPionIndexPrev[j]);
			AliKFParticle posPionCandidateKF( *posPionCandidate->GetConstrainedParam(), 211 );

			AliKFConversionPhoton* virtualPhoton = NULL;
			virtualPhoton = new AliKFConversionPhoton(negPionCandidateKF,posPionCandidateKF);
			AliKFVertex primaryVertexImproved(*fInputEvent->GetPrimaryVertex());
// 			primaryVertexImproved+=*virtualPhoton;
			virtualPhoton->SetProductionVertex(primaryVertexImproved);
			virtualPhoton->SetTrackLabels( lGoodPosPionIndexPrev[j], lGoodNegPionIndexPrev[i]);
			
			
			if( fMCEvent ) {
				Int_t labeln=TMath::Abs(negPionCandidate->GetLabel());
				Int_t labelp=TMath::Abs(posPionCandidate->GetLabel());
                TParticle *fNegativeMCParticle = fMCEvent->Particle(labeln);
                TParticle *fPositiveMCParticle = fMCEvent->Particle(labelp);
				
				if( fPositiveMCParticle && fNegativeMCParticle) {
					virtualPhoton->SetMCLabelPositive(labelp);
					virtualPhoton->SetMCLabelNegative(labeln);
				}
			}
			
			AliAODConversionPhoton *vParticle = new AliAODConversionPhoton(virtualPhoton); //To Apply PsiPairCut
			if (((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->DoMassCut()){
				if (vParticle->GetMass() < ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetMassCut()){
					fGoodVirtualParticles->Add(  vParticle );
				}
			} else {
				fGoodVirtualParticles->Add(  vParticle );
			}	
			delete virtualPhoton;
			virtualPhoton=NULL;
					
		}
	}
}

//_____________________________________________________________________________
void AliAnalysisTaskEtaToPiPlPiMiGamma::ProcessMCParticles(){

	// Loop over all primary MC particle

	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();
	
	
    for(Int_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
        if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){
            TParticle* particle = (TParticle *)fMCEvent->Particle(i);
			if (!particle) continue;

			Int_t isMCFromMBHeader = -1;
			if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
				isMCFromMBHeader
                    = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
				if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
			}

            if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent,fInputEvent)){
				
                if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kFALSE)){
					fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt()); // All MC Gamma
					if(particle->GetMother(0) >-1){
                        if (fMCEvent->Particle(particle->GetMother(0))->GetPdgCode() ==221 && fMCEvent->Particle(particle->GetMother(0))->GetNDaughters()==3 ) fHistoMCGammaFromEtaPt[fiCut]->Fill(particle->Pt()); // All pos from eta
					}	
				}
				
                if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kTRUE)){
					fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt());
				} // Converted MC Gamma
				
                if(((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(i,fMCEvent)){
					if( particle->GetPdgCode() == 211){
						fHistoMCAllPosPionsPt[fiCut]->Fill(particle->Pt()); // All pos pions
						if(particle->GetMother(0) >-1){
                            if (fMCEvent->Particle(particle->GetMother(0))->GetPdgCode() ==221) fHistoMCPosPionsFromEtaPt[fiCut]->Fill(particle->Pt()); // All pos from eta
						}	
					}	
					if( particle->GetPdgCode() == -211){
						fHistoMCAllNegPionsPt[fiCut]->Fill(particle->Pt()); // All neg pions
						if(particle->GetMother(0) >-1){
                            if (fMCEvent->Particle(particle->GetMother(0))->GetPdgCode() ==221) fHistoMCNegPionsFromEtaPt[fiCut]->Fill(particle->Pt()); // All pos from eta
						}	
					}
				}
				
				// \eta -> \gamma \gamma 
				
                if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMC( particle,fMCEvent,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift() ) ){
					Float_t weighted= 1;
					if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) { 
                        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
							if (particle->Pt()>0.005){
                                weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCEvent,fInputEvent);
							}
						}
					}
					if(particle->GetPdgCode() == 221)fHistoMCEtaGGPt[fiCut]->Fill( particle->Pt() , weighted); // All MC Eta GG decay
				}
				
				// \eta -> e+ e- \gamma 
				Int_t labelgamma 	  = -1;
				Int_t labelelectron = -1;
				Int_t labelpositron = -1;

                if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMCDalitz(particle,fMCEvent,labelelectron,labelpositron,labelgamma,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
					Float_t weighted= 1;
					if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) { 
                        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent,fInputEvent)){
							if (particle->Pt()>0.005){
                                weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCEvent,fInputEvent);
							}
						}
					}
					if(particle->GetPdgCode() == 221)fHistoMCEtaDalitzPt[fiCut]->Fill(particle->Pt(), weighted); // All MC Eta
				}		
				
				
				// \eta -> pi+ pi- \gamma 
				Int_t labelGamma3Body 	  = -1;
				Int_t labelNegPion = -1;
				Int_t labelPosPion = -1;

                if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMCEtaPiPlPiMiGamma(particle,fMCEvent,labelNegPion,labelPosPion,labelGamma3Body,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
					Float_t weighted= 1;
					if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) { 
                        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent,fInputEvent)){
							if (particle->Pt()>0.005){
                                weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCEvent,fInputEvent);
							}
						}
					}
					if(particle->GetPdgCode() == 221)fHistoMCEtaPiPlPiMiGammaPt[fiCut]->Fill(particle->Pt(), weighted); // All MC Eta
			
                    TParticle *gamma    = fMCEvent->Particle(labelGamma3Body);
                    Bool_t gammaIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelGamma3Body, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
                    Bool_t negPionIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelNegPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
                    Bool_t posPionIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelPosPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
					
					if( gammaIsPrimary && negPionIsPrimary && posPionIsPrimary &&
                        ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(gamma,fMCEvent,kFALSE) &&
                        ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelNegPion,fMCEvent) &&
                        ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelPosPion,fMCEvent) ) {
						if(particle->GetPdgCode() == 221)fHistoMCEtaPiPlPiMiGammaInAccPt[fiCut]->Fill(particle->Pt(), weighted ); // MC EtaDalitz with gamma and e+e- in acc
					}				
				}
			}
		}	
	}
}


//________________________________________________________________________
void AliAnalysisTaskEtaToPiPlPiMiGamma::CalculateMesonCandidates(){

	// Conversion Gammas
	if( fGoodGammas->GetEntries() > 0 && fGoodVirtualParticles->GetEntries() > 0 ){

		vector<Bool_t> lGoodVirtualParticle(fGoodVirtualParticles->GetEntries(), kFALSE);
		
		for(Int_t GammaIndex=0; GammaIndex<fGoodGammas->GetEntries(); GammaIndex++){

			AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fGoodGammas->At(GammaIndex));
			if (gamma==NULL) continue;
			for(Int_t virtualParticleIndex=0;virtualParticleIndex<fGoodVirtualParticles->GetEntries();virtualParticleIndex++){

				AliAODConversionPhoton *vParticle=dynamic_cast<AliAODConversionPhoton*>(fGoodVirtualParticles->At(virtualParticleIndex));
				if (vParticle==NULL) continue;
				//Check for same Electron ID
				if(gamma->GetTrackLabelPositive() == vParticle->GetTrackLabelPositive() ||
				gamma->GetTrackLabelNegative() == vParticle->GetTrackLabelNegative() ||
				gamma->GetTrackLabelNegative() == vParticle->GetTrackLabelPositive() ||
				gamma->GetTrackLabelPositive() == vParticle->GetTrackLabelNegative() ) continue;

				AliAODConversionMother *etacand = new AliAODConversionMother(gamma,vParticle);
				etacand->SetLabels(GammaIndex,virtualParticleIndex);
				
// 				if(fMCEvent){
// 					AliESDtrack *posPionVParticle = fESDEvent->GetTrack( vParticle->GetTrackLabelNegative() );
// 					AliESDtrack *negPionVParticle = fESDEvent->GetTrack( vParticle->GetTrackLabelPositive() );
// 
// 					Int_t labeln=TMath::Abs(negPionVParticle->GetLabel());
// 					Int_t labelp=TMath::Abs(posPionVParticle->GetLabel());
// 					
// 					cout << labeln << "\t" << labelp << endl;
// 				}	

				if( ( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(etacand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())) ){
			
// 					cout<< "Meson Accepted "<<endl;
					
					Int_t zbin= fBGHandler[fiCut]->GetZBinIndex(fESDEvent->GetPrimaryVertex()->GetZ());
					Int_t mbin = 0;
					if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
						mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
					} else {
						mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGoodGammas->GetEntries());
					}
					
					AliESDtrack *posPionVParticle = 0;
					AliESDtrack *negPionVParticle = 0;
					
					Double_t clsToFPos = -1.0;
					Double_t clsToFNeg = -1.0;
					
					Float_t dcaToVertexXYPos = -1.0;
					Float_t dcaToVertexZPos  = -1.0;
					Float_t dcaToVertexXYNeg = -1.0;
					Float_t dcaToVertexZNeg  = -1.0;
					
					
					if ( fDoMesonQA ) {
					
						fHistoPionPionInvMassPt[fiCut]->Fill( vParticle->GetMass(),vParticle->Pt());
						
						posPionVParticle = fESDEvent->GetTrack( vParticle->GetTrackLabelPositive() );
						negPionVParticle = fESDEvent->GetTrack( vParticle->GetTrackLabelNegative() );
						clsToFPos = ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetNFindableClustersTPC(posPionVParticle);
						clsToFNeg = ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetNFindableClustersTPC(negPionVParticle);
						
						Float_t bPos[2];
						Float_t bCovPos[3];
						posPionVParticle->GetImpactParameters(bPos,bCovPos);
						if (bCovPos[0]<=0 || bCovPos[2]<=0) {
							AliDebug(1, "Estimated b resolution lower or equal zero!");
							bCovPos[0]=0; bCovPos[2]=0;
						}
						
						Float_t bNeg[2];
						Float_t bCovNeg[3];
						posPionVParticle->GetImpactParameters(bNeg,bCovNeg);
						if (bCovNeg[0]<=0 || bCovNeg[2]<=0) {
							AliDebug(1, "Estimated b resolution lower or equal zero!");
							bCovNeg[0]=0; bCovNeg[2]=0;
						}
						
						dcaToVertexXYPos = bPos[0];
						dcaToVertexZPos  = bPos[1];
						dcaToVertexXYNeg = bNeg[0];
						dcaToVertexZNeg  = bNeg[1];
					}
					
					fHistoMotherInvMassPt[fiCut]->Fill(etacand->M(),etacand->Pt());
					Double_t sparesFill[4] = {etacand->M(),etacand->Pt(),(Double_t)zbin,(Double_t)mbin};
					fTHnSparseMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
					
					if ( fDoMesonQA ) {
						if( lGoodVirtualParticle[virtualParticleIndex] == kFALSE ) {
					
							fHistoNegPionEta[fiCut]->Fill( negPionVParticle->Eta() );
							fHistoPosPionEta[fiCut]->Fill( posPionVParticle->Eta() );
									
							fHistoNegPionClsTPC[fiCut]->Fill(clsToFNeg,negPionVParticle->Pt());
							fHistoPosPionClsTPC[fiCut]->Fill(clsToFPos,posPionVParticle->Pt());
							
							fHistoPionDCAxy[fiCut]->Fill(  dcaToVertexXYNeg, negPionVParticle->Pt() );
							fHistoPionDCAz[fiCut]->Fill(   dcaToVertexZNeg,  negPionVParticle->Pt() );
							fHistoPionDCAxy[fiCut]->Fill(  dcaToVertexXYPos, posPionVParticle->Pt() );
							fHistoPionDCAz[fiCut]->Fill(   dcaToVertexZPos,  posPionVParticle->Pt() );
							
							fHistoPionTPCdEdxNSigma[fiCut]->Fill( posPionVParticle->P(),((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetPIDResponse()->NumberOfSigmasTPC(posPionVParticle, AliPID::kPion) );
							fHistoPionTPCdEdxNSigma[fiCut]->Fill( negPionVParticle->P(),((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetPIDResponse()->NumberOfSigmasTPC(negPionVParticle, AliPID::kPion) );
							
							fHistoPionTPCdEdx[fiCut]->Fill( posPionVParticle->P(), TMath::Abs(posPionVParticle->GetTPCsignal()));
							fHistoPionTPCdEdx[fiCut]->Fill( negPionVParticle->P(), TMath::Abs(negPionVParticle->GetTPCsignal()));
							
							lGoodVirtualParticle[virtualParticleIndex] = kTRUE;
					
						}
					}
		
				
					if(fMCEvent){
						ProcessTrueMesonCandidates(etacand,gamma,vParticle);
					}
				}
				delete etacand;
				etacand=0x0;
			}
		}
	}
}

//________________________________________________________________________
void AliAnalysisTaskEtaToPiPlPiMiGamma::CalculateBackground(){

	Int_t zbin= fBGHandler[fiCut]->GetZBinIndex(fESDEvent->GetPrimaryVertex()->GetZ());
	Int_t mbin = 0;


	if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
		mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
	} else {
		mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGoodGammas->GetEntries());
	}

	Int_t method = 1;
	AliGammaConversionAODBGHandler::GammaConversionVertex *bgEventVertex = NULL;
	if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity() ) {
		for(Int_t nEventsInBG=0;nEventsInBG<fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
			AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
			if(fMoveParticleAccordingToVertex == kTRUE && method == 1){
				bgEventVertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
			}

            for(Int_t iCurrent=0;iCurrent<fGoodVirtualParticles->GetEntries();iCurrent++){
				AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGoodVirtualParticles->At(iCurrent));

                for(Int_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
					AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));

					if(fMoveParticleAccordingToVertex == kTRUE && method == 1 ){
						MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
					}

					AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
								

					if( ( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE, ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
						fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
						Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
						fTHnSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
					}
					delete backgroundCandidate;
					backgroundCandidate = 0x0;
				}
			}
		}
	} else {
		for(Int_t nEventsInBG=0;nEventsInBG <fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
			AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
			if(previousEventV0s){
				if(fMoveParticleAccordingToVertex == kTRUE && method == 1){
					bgEventVertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
				}
                for(Int_t iCurrent=0;iCurrent<fGoodVirtualParticles->GetEntries();iCurrent++){

					AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGoodVirtualParticles->At(iCurrent));

                    for(Int_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){

						AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));

						if(fMoveParticleAccordingToVertex == kTRUE && method ==1){
							MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
						}

						AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
								
						if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
							fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
							Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
							fTHnSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
						}
						delete backgroundCandidate;
						backgroundCandidate = 0x0;
					}
				}
			}
		}
	}
}

//______________________________________________________________________
void AliAnalysisTaskEtaToPiPlPiMiGamma::ProcessTrueMesonCandidates(AliAODConversionMother *EtaCandidate, AliAODConversionPhoton *TrueGammaCandidate, AliAODConversionPhoton *TrueVirtualParticleCandidate){

	// Process True Mesons

	if(	TrueGammaCandidate->GetV0Index()<fESDEvent->GetNumberOfV0s()	){
		
		Bool_t isTrueEta = kFALSE;
        Int_t gammaMCLabel = TrueGammaCandidate->GetMCParticleLabel(fMCEvent);
		Int_t gammaMotherLabel = -1;
// 		Bool_t gammaEtaCand = kFALSE;
		
		if(gammaMCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
			// Daughters Gamma 0
            TParticle * negativeMC = (TParticle*)TrueGammaCandidate->GetNegativeMCDaughter(fMCEvent);
            TParticle * positiveMC = (TParticle*)TrueGammaCandidate->GetPositiveMCDaughter(fMCEvent);
            TParticle * gammaMC = (TParticle*)fMCEvent->Particle(gammaMCLabel);

			if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
				if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
					if(gammaMC->GetPdgCode() == 22){ // ... with Gamma Mother
						gammaMotherLabel=gammaMC->GetFirstMother();
// 						if( ((TParticle*)fMCEvent->Particle(gammaMotherLabel))->GetNDaughters() == 3 && ((TParticle*)fMCEvent->Particle(gammaMotherLabel))->GetPdgCode() == 221 ) gammaEtaCand = kTRUE;
					}
				}
			}
			
			
		}

        Int_t virtualParticleMCLabel = TrueVirtualParticleCandidate->GetMCParticleLabel(fMCEvent);
		Int_t virtualParticleMotherLabel = -1;

		Bool_t isPiPiDecay = kFALSE;
		Bool_t isDalitz = kFALSE;
		Bool_t isRealGamma = kFALSE;
		
		if (fDoMesonQA){
            TParticle * negativeMC = (TParticle*)TrueVirtualParticleCandidate->GetNegativeMCDaughter(fMCEvent);
            TParticle * positiveMC = (TParticle*)TrueVirtualParticleCandidate->GetPositiveMCDaughter(fMCEvent);
// 			if (gammaEtaCand){
// 				cout << "neg Part: label - " <<  TrueVirtualParticleCandidate->GetMCLabelNegative() <<" pdg-code - " << negativeMC->GetPdgCode() << endl;
// 				cout << "pos Part: label - " <<  TrueVirtualParticleCandidate->GetMCLabelPositive() <<" pdg-code - " << positiveMC->GetPdgCode() << endl;			
// 			}
			if(TMath::Abs(negativeMC->GetPdgCode())==211 && TMath::Abs(positiveMC->GetPdgCode())==211){  // Pions ...
				fHistoTruePionPionInvMassPt[fiCut]->Fill(TrueVirtualParticleCandidate->GetMass(),TrueVirtualParticleCandidate->Pt());
			}
		}
		
		if(virtualParticleMCLabel != -1){ // if virtualParticleMCLabel==-1 particles don't have same mother 
            TParticle * negativeMC = (TParticle*)TrueVirtualParticleCandidate->GetNegativeMCDaughter(fMCEvent);
            TParticle * positiveMC = (TParticle*)TrueVirtualParticleCandidate->GetPositiveMCDaughter(fMCEvent);
            TParticle * virtualParticleMotherMC = (TParticle*)fMCEvent->Particle(virtualParticleMCLabel);
// 			cout << "pdg code same mother - " << virtualParticleMotherMC->GetPdgCode() << endl;
			
			if(TMath::Abs(negativeMC->GetPdgCode())==211 && TMath::Abs(positiveMC->GetPdgCode())==211){  // Pions ...
				virtualParticleMotherLabel=virtualParticleMCLabel;
				isPiPiDecay=kTRUE;
			} else if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
				if( virtualParticleMotherMC->GetPdgCode() != 22 ){
					virtualParticleMotherLabel=virtualParticleMCLabel;
					isDalitz = kTRUE;
				} else if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
					virtualParticleMotherLabel=virtualParticleMotherMC->GetFirstMother();
					isRealGamma = kTRUE; //no virtual gamma
				}
			}	
		}

		if(gammaMotherLabel >= 0 && ( gammaMotherLabel == virtualParticleMotherLabel) ){
            if(((TParticle*)fMCEvent->Particle(virtualParticleMotherLabel))->GetPdgCode() == 221){
				isTrueEta=kTRUE;
				if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,gammaMotherLabel)) fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(EtaCandidate->M(),EtaCandidate->Pt());
			}
		}

		if( isTrueEta ){ // True Eta
			if ( isPiPiDecay) { //real eta -> Pi+ Pi- Gamma
				Float_t weighted= 1;
				if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) { 
                    if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gammaMotherLabel, fMCEvent,fInputEvent)){
                        if (((TParticle*)fMCEvent->Particle(gammaMotherLabel))->Pt()>0.005){
                            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gammaMotherLabel,fMCEvent,fInputEvent);
						}
					}
				}
				fHistoTruePionPionFromEtaInvMassPt[fiCut]->Fill(TrueVirtualParticleCandidate->GetMass(),TrueVirtualParticleCandidate->Pt());
				fHistoTrueMotherPiPlPiMiGammaInvMassPt[fiCut]->Fill(EtaCandidate->M(),EtaCandidate->Pt(),weighted);
			} else if ( isRealGamma ){
				Float_t weighted= 1;
				if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) {
                    if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gammaMotherLabel, fMCEvent,fInputEvent)){
                        if (((TParticle*)fMCEvent->Particle(gammaMotherLabel))->Pt()>0.005){
                            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gammaMotherLabel,fMCEvent,fInputEvent);
						}
					}
				}

				fHistoTrueMotherGammaGammaInvMassPt[fiCut]->Fill(EtaCandidate->M(),EtaCandidate->Pt(),weighted); 
			} else if (isDalitz) {
				Float_t weighted= 1;
				if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) {
                    if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gammaMotherLabel, fMCEvent,fInputEvent)){
                        if (((TParticle*)fMCEvent->Particle(gammaMotherLabel))->Pt()>0.005){
                            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gammaMotherLabel,fMCEvent,fInputEvent);
						}
					}
				}
				fHistoTrueMotherDalitzInvMassPt[fiCut]->Fill(EtaCandidate->M(),EtaCandidate->Pt(),weighted);
			}
		}
	}
}


//________________________________________________________________________
void AliAnalysisTaskEtaToPiPlPiMiGamma::UpdateEventByEventData(){
   //see header file for documentation

	Int_t method = 1;
	if( method == 1 ) {
		if(fGoodGammas->GetEntries() >0 ){
			if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
				fBGHandler[fiCut]->AddEvent(fGoodGammas,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),0);
			} else{ // means we use #V0s for multiplicity
				fBGHandler[fiCut]->AddEvent(fGoodGammas,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGoodGammas->GetEntries(),0);
			}
		}
	} else if ( method == 2 ){
		if(fGoodVirtualParticles->GetEntries() > 0 ){
			if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
				fBGHandler[fiCut]->AddEvent(fGoodVirtualParticles,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),0);
			} else{ // means we use #V0s for multiplicity
				fBGHandler[fiCut]->AddEvent(fGoodVirtualParticles,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGoodVirtualParticles->GetEntries(),0);
			}
		}
	}
}

//________________________________________________________________________
void AliAnalysisTaskEtaToPiPlPiMiGamma::MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex){
   //see header file for documentation

   Double_t dx = vertex->fX - fInputEvent->GetPrimaryVertex()->GetX();
   Double_t dy = vertex->fY - fInputEvent->GetPrimaryVertex()->GetY();
   Double_t dz = vertex->fZ - fInputEvent->GetPrimaryVertex()->GetZ();

   Double_t movedPlace[3] = {particle->GetConversionX() - dx,particle->GetConversionY() - dy,particle->GetConversionZ() - dz};
   particle->SetConversionPoint(movedPlace);
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEtaToPiPlPiMiGamma::IsEtaPiPlPiMiGammaDaughter( Int_t label ) const {
//
// Returns true if the particle comes from eta -> pi+ pi- gamma
//
    Int_t motherLabel = fMCEvent->Particle( label )->GetMother(0);
    if( motherLabel < 0 || motherLabel >= fMCEvent->GetNumberOfTracks() ) return kFALSE;
	
    TParticle* mother = fMCEvent->Particle( motherLabel );
	if( mother->GetPdgCode() != 221 ) return kFALSE;
	if( IsPiPlPiMiGammaDecay( mother ) ) return kTRUE;	
	return kFALSE;       
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskEtaToPiPlPiMiGamma::IsPiPlPiMiGammaDecay(TParticle *fMCMother) const
{

	if( fMCMother->GetNDaughters() != 3 ) return kFALSE;
	if( fMCMother->GetPdgCode() != 221 ) return kFALSE;
	
	
	TParticle *posPion = 0x0;
	TParticle *negPion = 0x0;
	TParticle *gamma    = 0x0;
	
	for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){				
        TParticle* temp = (TParticle*)fMCEvent->Particle( index );
		
		switch( temp->GetPdgCode() ) {
		case 211:
			posPion =  temp;
			break;
		case -211:
			negPion =  temp;
			break;
		case ::kGamma:
			gamma    =  temp;
			break;
		}
	}
  
	if( posPion && negPion && gamma) return kTRUE;
	
	return kFALSE;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskEtaToPiPlPiMiGamma::CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked)
{
	if(tobechecked > -1)
	{
		vector<Int_t>::iterator it;
		it = find (vec.begin(), vec.end(), tobechecked);
		if (it != vec.end()) return true;
		else{
			vec.push_back(tobechecked);
			return false;
		}
	}
	return false;
}
