/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Friederike Bock 									              *
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

//

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
#include "AliStack.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliMCEvent.h"
#include "AliStack.h"
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
#include "AliAODEvent.h"
#include "AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero.h"


ClassImp( AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero )

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero():
	fV0Reader(NULL),
	fPionSelector(NULL),
	fBGHandler(NULL),
	fESDEvent(NULL),
	fMCEvent(NULL),
	fMCStack(NULL),
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
	fGoodConvGammas(NULL),
	fClusterCandidates(NULL),
	fNeutralPionCandidates(NULL),
	fGoodVirtualParticles(NULL),
	fEventCutArray(NULL),
	fGammaCutArray(NULL),
	fClusterCutArray(NULL),
	fPionCutArray(NULL),
	fNeutralPionMesonCutArray(NULL),
	fMesonCutArray(NULL),
	fEventCuts(NULL),
	fConversionCuts(NULL),
	fClusterCuts(NULL),
	fHistoConvGammaPt(NULL),
	fHistoConvGammaEta(NULL),
	fHistoClusterGammaPt(NULL),
	fHistoClusterGammaEta(NULL),
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
	fHistoGammaGammaInvMassPt(NULL),
	fHistoMotherInvMassPt(NULL),
	fTHnSparseMotherInvMassPtZM(NULL),
	fHistoMotherBackInvMassPt(NULL),
	fTHnSparseMotherBackInvMassPtZM(NULL),
	fHistoMCAllGammaPt(NULL),
	fHistoMCConvGammaPt(NULL),
	fHistoMCAllPosPionsPt(NULL),
	fHistoMCAllNegPionsPt(NULL),
	fHistoMCGammaFromNeutralMesonPt(NULL),
	fHistoMCPosPionsFromNeutralMesonPt(NULL),
	fHistoMCNegPionsFromNeutralMesonPt(NULL),
	fHistoMCEtaPiPlPiMiPiZeroPt(NULL),
	fHistoMCEtaPiPlPiMiPiZeroInAccPt(NULL),
	fHistoMCOmegaPiPlPiMiPiZeroPt(NULL),
	fHistoMCOmegaPiPlPiMiPiZeroInAccPt(NULL),
	fHistoTrueMotherPiPlPiMiPiZeroInvMassPt(NULL),
	fHistoTrueMotherGammaGammaInvMassPt(NULL),
	fHistoTrueMotherGammaGammaFromEtaInvMassPt(NULL),
	fHistoTrueMotherGammaGammaFromOmegaInvMassPt(NULL),
	fHistoTrueConvGammaPt(NULL),
	fHistoTrueConvGammaFromNeutralMesonPt(NULL),
	fHistoTrueClusterGammaPt(NULL),
	fHistoTrueClusterGammaFromNeutralMesonPt(NULL),
	fHistoTruePosPionPt(NULL),
	fHistoTruePosPionFromNeutralMesonPt(NULL),
	fHistoTrueNegPionPt(NULL),
	fHistoTrueNegPionFromNeutralMesonPt(NULL),
	fHistoTruePionPionInvMassPt(NULL),
	fHistoTruePionPionFromSameMotherInvMassPt(NULL),
	fHistoTruePionPionFromEtaInvMassPt(NULL),
	fHistoTruePionPionFromOmegaInvMassPt(NULL),
	fHistoNEvents(NULL),
	fHistoNGoodESDTracks(NULL),
	fProfileEtaShift(NULL),
	fRandom(0),
	fnCuts(0),
	fiCut(0),
	fNumberOfESDTracks(0),
	fMoveParticleAccordingToVertex(kFALSE),
	fIsHeavyIon(0),
	fDoMesonAnalysis(kTRUE),
	fDoMesonQA(kFALSE),
	fIsFromMBHeader(kTRUE),
	fIsMC(kFALSE),
	fNeutralPionMode(0)
{

}

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero( const char* name ):
	AliAnalysisTaskSE(name),
	fV0Reader(NULL),
	fPionSelector(NULL),
	fBGHandler(NULL),
	fESDEvent(NULL),
	fMCEvent(NULL),
	fMCStack(NULL),
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
	fGoodConvGammas(NULL),
	fClusterCandidates(NULL),
	fNeutralPionCandidates(NULL),	
	fGoodVirtualParticles(NULL),
	fEventCutArray(NULL),
	fGammaCutArray(NULL),
	fClusterCutArray(NULL),
	fPionCutArray(NULL),
	fNeutralPionMesonCutArray(NULL),
	fMesonCutArray(NULL),
	fEventCuts(NULL),
	fConversionCuts(NULL),
	fClusterCuts(NULL),
	fHistoConvGammaPt(NULL),
	fHistoConvGammaEta(NULL),
	fHistoClusterGammaPt(NULL),
	fHistoClusterGammaEta(NULL),
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
	fHistoGammaGammaInvMassPt(NULL),
	fHistoMotherInvMassPt(NULL),
	fTHnSparseMotherInvMassPtZM(NULL),
	fHistoMotherBackInvMassPt(NULL),
	fTHnSparseMotherBackInvMassPtZM(NULL),
	fHistoMCAllGammaPt(NULL),
	fHistoMCConvGammaPt(NULL),
	fHistoMCAllPosPionsPt(NULL),
	fHistoMCAllNegPionsPt(NULL),
	fHistoMCGammaFromNeutralMesonPt(NULL),
	fHistoMCPosPionsFromNeutralMesonPt(NULL),
	fHistoMCNegPionsFromNeutralMesonPt(NULL),
	fHistoMCEtaPiPlPiMiPiZeroPt(NULL),
	fHistoMCEtaPiPlPiMiPiZeroInAccPt(NULL),
	fHistoMCOmegaPiPlPiMiPiZeroPt(NULL),
	fHistoMCOmegaPiPlPiMiPiZeroInAccPt(NULL),
	fHistoTrueMotherPiPlPiMiPiZeroInvMassPt(NULL),
	fHistoTrueMotherGammaGammaInvMassPt(NULL),
	fHistoTrueMotherGammaGammaFromEtaInvMassPt(NULL),
	fHistoTrueMotherGammaGammaFromOmegaInvMassPt(NULL),
	fHistoTrueConvGammaPt(NULL),
	fHistoTrueConvGammaFromNeutralMesonPt(NULL),
	fHistoTrueClusterGammaPt(NULL),
	fHistoTrueClusterGammaFromNeutralMesonPt(NULL),
	fHistoTruePosPionPt(NULL),
	fHistoTruePosPionFromNeutralMesonPt(NULL),
	fHistoTrueNegPionPt(NULL),
	fHistoTrueNegPionFromNeutralMesonPt(NULL),
	fHistoTruePionPionInvMassPt(NULL),
	fHistoTruePionPionFromSameMotherInvMassPt(NULL),
	fHistoTruePionPionFromEtaInvMassPt(NULL),
	fHistoTruePionPionFromOmegaInvMassPt(NULL),
	fHistoNEvents(NULL),
	fHistoNGoodESDTracks(NULL),
	fProfileEtaShift(NULL),
	fRandom(0),
	fnCuts(0),
	fiCut(0),
	fNumberOfESDTracks(0),
	fMoveParticleAccordingToVertex(kFALSE),
	fIsHeavyIon(0),
	fDoMesonAnalysis(kTRUE),
	fDoMesonQA(kFALSE),
	fIsFromMBHeader(kTRUE),
	fIsMC(kFALSE),
	fNeutralPionMode(0)
{
   DefineOutput(1, TList::Class());
}

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::~AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero()
{
	//
	// virtual destructor
	//
	cout<<"Destructor"<<endl;
	if(fGoodConvGammas){
		delete fGoodConvGammas;
		fGoodConvGammas = 0x0;
	}
	if(fClusterCandidates){
		delete fClusterCandidates;
		fClusterCandidates = 0x0;
	}

	if(fNeutralPionCandidates){
		delete fNeutralPionCandidates;
		fNeutralPionCandidates = 0x0;
	}

	if(fGoodVirtualParticles){
		delete fGoodVirtualParticles;
		fGoodVirtualParticles = 0x0;
	}
	if(fBGHandler){
		delete[] fBGHandler;
		fBGHandler = 0x0;
	}
}
//___________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::InitBack(){

	const Int_t nDim = 4;
	Int_t nBins[nDim] = {500,250,7,4};
	Double_t xMin[nDim] = {0.4,0, 0,0};
	Double_t xMax[nDim] = {0.9,25,7,4};
	
	fTHnSparseMotherInvMassPtZM = new THnSparseF*[fnCuts];
	fTHnSparseMotherBackInvMassPtZM = new THnSparseF*[fnCuts];

	fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];

	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		
		TString cutstringEvent		= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
		TString cutstringPion		= ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutNumber();
		TString cutstringConvGamma = "";
		if (fNeutralPionMode < 2)  cutstringConvGamma = ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutNumber();
		TString cutstringCaloGamma = "";
		if (fNeutralPionMode > 0)  cutstringCaloGamma = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
		TString cutstringNeutralPion= ((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(iCut))->GetCutNumber();
		TString cutstringMeson		= ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
		
		TString fullCutString = "";
		if (fNeutralPionMode == 0) fullCutString = Form("%i_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
		else if (fNeutralPionMode == 1) fullCutString = Form("%i_%s_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringCaloGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
		else if (fNeutralPionMode == 2) fullCutString = Form("%i_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringCaloGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
		
		TString nameBackList = Form("%s Back histograms",fullCutString.Data());
		TString nameMotherList = Form("%s Mother histograms",fullCutString.Data());
		
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
		
		fBackList[iCut]->SetName(nameBackList.Data());
		fBackList[iCut]->SetOwner(kTRUE);
		fCutFolder[iCut]->Add(fBackList[iCut]);

		fTHnSparseMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_m","Back_Back_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
		fBackList[iCut]->Add(fTHnSparseMotherBackInvMassPtZM[iCut]);

		fMotherList[iCut] = new TList();
		fMotherList[iCut]->SetName(nameMotherList.Data());
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
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::UserCreateOutputObjects()
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

	fGoodConvGammas = new TList();
	fClusterCandidates = new TList();
	fNeutralPionCandidates = new TList();
	fGoodVirtualParticles = new TList();
	
	fCutFolder				= new TList*[fnCuts];
	fESDList				= new TList*[fnCuts];
	fBackList				= new TList*[fnCuts];
	fMotherList 			= new TList*[fnCuts];
	fHistoNEvents			= new TH1I*[fnCuts];
	fHistoNGoodESDTracks	= new TH1I*[fnCuts];
	fProfileEtaShift		= new TProfile*[fnCuts];
	if (fNeutralPionMode < 2){
		fHistoConvGammaPt		= new TH1F*[fnCuts];
		fHistoConvGammaEta		= new TH1F*[fnCuts];
	}
	if (fNeutralPionMode > 0){
		fHistoClusterGammaPt		= new TH1F*[fnCuts];
		fHistoClusterGammaEta		= new TH1F*[fnCuts];
	}
	fHistoNegPionPt			= new TH1F*[fnCuts];
	fHistoPosPionPt			= new TH1F*[fnCuts];
	fHistoNegPionPhi		= new TH1F*[fnCuts];
	fHistoPosPionPhi		= new TH1F*[fnCuts];
	fHistoPionPionInvMassPt	= new TH2F*[fnCuts];
	
	if( fDoMesonQA ) {			
		fHistoNegPionEta		= new TH1F*[fnCuts];
		fHistoPosPionEta		= new TH1F*[fnCuts];
		fHistoNegPionClsTPC		= new TH2F*[fnCuts];
		fHistoPosPionClsTPC		= new TH2F*[fnCuts];
		fHistoPionDCAxy			= new TH2F*[fnCuts];
		fHistoPionDCAz			= new TH2F*[fnCuts];
		fHistoPionTPCdEdxNSigma	= new TH2F*[fnCuts];
		fHistoPionTPCdEdx		= new TH2F*[fnCuts];
	}
	fHistoGammaGammaInvMassPt	= new TH2F*[fnCuts];
	fHistoMotherInvMassPt		= new TH2F*[fnCuts];
	fHistoMotherBackInvMassPt	= new TH2F*[fnCuts];

	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		TString cutstringEvent		= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
		TString cutstringPion		= ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutNumber();
		TString cutstringConvGamma = "";
		if (fNeutralPionMode < 2)  cutstringConvGamma = ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutNumber();
		TString cutstringCaloGamma = "";
		if (fNeutralPionMode > 0)  cutstringCaloGamma = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
		TString cutstringNeutralPion= ((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(iCut))->GetCutNumber();
		TString cutstringMeson		= ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
		
		TString fullCutString = "";
		if (fNeutralPionMode == 0) fullCutString = Form("%i_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
		else if (fNeutralPionMode == 1) fullCutString = Form("%i_%s_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringCaloGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
		else if (fNeutralPionMode == 2) fullCutString = Form("%i_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringCaloGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
		TString nameCutFolder = Form("Cut Number %s", fullCutString.Data());
		TString nameESDList = Form("%s ESD histograms", fullCutString.Data());
		
		
		fCutFolder[iCut] = new TList();
		fCutFolder[iCut]->SetName(nameCutFolder.Data());
		fCutFolder[iCut]->SetOwner(kTRUE);
		fOutputContainer->Add(fCutFolder[iCut]);

		fESDList[iCut] = new TList();
		fESDList[iCut]->SetName(nameESDList.Data());
		fESDList[iCut]->SetOwner(kTRUE);

		fHistoNEvents[iCut] = new TH1I("NEvents","NEvents",10,-0.5,9.5);
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Missing MC");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problem");
		fESDList[iCut]->Add(fHistoNEvents[iCut]);

		if(fIsHeavyIon>0) fHistoNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",3000,0,3000);
			else fHistoNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",200,0,200);
		fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);

		fProfileEtaShift[iCut] = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
		fESDList[iCut]->Add(fProfileEtaShift[iCut]);
		if (fNeutralPionMode < 2){
			fHistoConvGammaPt[iCut] = new TH1F("ESD_ConvGamma_Pt","ESD_ConvGamma_Pt",250,0,25);
			fESDList[iCut]->Add(fHistoConvGammaPt[iCut]);
			fHistoConvGammaEta[iCut] = new TH1F("ESD_ConvGamma_Eta","ESD_ConvGamma_Eta",600,-1.5,1.5);
			fESDList[iCut]->Add(fHistoConvGammaEta[iCut]);
		}
		if (fNeutralPionMode > 0){
			fHistoClusterGammaPt[iCut] = new TH1F("ESD_ClusterGamma_Pt","ESD_ClusterGamma_Pt",250,0,25);
			fESDList[iCut]->Add(fHistoClusterGammaPt[iCut]);
			fHistoClusterGammaEta[iCut] = new TH1F("ESD_ClusterGamma_Eta","ESD_ClusterGamma_Eta",600,-1.5,1.5);
			fESDList[iCut]->Add(fHistoClusterGammaEta[iCut]);
		}
		fHistoNegPionPt[iCut] = new TH1F("ESD_PrimaryNegPions_Pt","ESD_PrimaryNegPions_Pt",1000,0,25);
		fESDList[iCut]->Add(fHistoNegPionPt[iCut]);
		fHistoPosPionPt[iCut] = new TH1F("ESD_PrimaryPosPions_Pt","ESD_PrimaryPosPions_Pt",1000,0,25);
		fESDList[iCut]->Add(fHistoPosPionPt[iCut]);
		fHistoNegPionPhi[iCut] = new TH1F("ESD_PrimaryNegPions_Phi","ESD_PrimaryNegPions_Phi",360,0,2*TMath::Pi());
		fESDList[iCut]->Add(fHistoNegPionPhi[iCut]);
		fHistoPosPionPhi[iCut] = new TH1F("ESD_PrimaryPosPions_Phi","ESD_PrimaryPosPions_Phi",360,0,2*TMath::Pi());
		fESDList[iCut]->Add(fHistoPosPionPhi[iCut]);
		fHistoPionPionInvMassPt[iCut] = new TH2F("ESD_PiPlusPiNeg_InvMassPt","ESD_PiPlusPiNeg_InvMassPt",2000,0.,2.,200,0.,20.);
		fESDList[iCut]->Add(fHistoPionPionInvMassPt[iCut]);
		
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
		}
		fHistoGammaGammaInvMassPt[iCut] = new TH2F("ESD_GammaGamma_InvMass_Pt","ESD_GammaGamma_InvMass_Pt",450,0.,0.45,250,0,25);
		fESDList[iCut]->Add(fHistoGammaGammaInvMassPt[iCut]);
		fHistoMotherInvMassPt[iCut] = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",500,0.4,0.9,250,0,25);
		fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);
		fHistoMotherBackInvMassPt[iCut] = new TH2F("ESD_Background_InvMass_Pt","ESD_Background_InvMass_Pt",500,0.4,0.9,250,0,25);
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
		if (fNeutralPionMode < 2){
			fHistoTrueConvGammaPt = new TH1F*[fnCuts];
			fHistoTrueConvGammaFromNeutralMesonPt = new TH1F*[fnCuts];
		}	
		if (fNeutralPionMode > 0){
			fHistoTrueClusterGammaPt = new TH1F*[fnCuts];
			fHistoTrueClusterGammaFromNeutralMesonPt = new TH1F*[fnCuts];
		}	
		fHistoTruePosPionPt  = new TH1F*[fnCuts];
		fHistoTrueNegPionPt  = new TH1F*[fnCuts];		
		fHistoTruePosPionFromNeutralMesonPt  = new TH1F*[fnCuts];
		fHistoTrueNegPionFromNeutralMesonPt  = new TH1F*[fnCuts];
		

		fHistoMCAllGammaPt  = new TH1F*[fnCuts];
		if (fNeutralPionMode < 2){
			fHistoMCConvGammaPt = new TH1F*[fnCuts];
		}	
		fHistoMCAllPosPionsPt = new TH1F*[fnCuts];
		fHistoMCAllNegPionsPt = new TH1F*[fnCuts];
		fHistoMCGammaFromNeutralMesonPt  = new TH1F*[fnCuts];
		fHistoMCPosPionsFromNeutralMesonPt = new TH1F*[fnCuts];
		fHistoMCNegPionsFromNeutralMesonPt = new TH1F*[fnCuts];

// 		hMCPi0DalitzGammaPt    = new TH1F*[fnCuts];
// 		hMCPi0DalitzElectronPt = new TH1F*[fnCuts];
// 		hMCPi0DalitzPositronPt = new TH1F*[fnCuts];

		fHistoMCEtaPiPlPiMiPiZeroPt = new TH1F*[fnCuts];
		fHistoMCEtaPiPlPiMiPiZeroInAccPt = new TH1F*[fnCuts];
		fHistoMCOmegaPiPlPiMiPiZeroPt = new TH1F*[fnCuts];
		fHistoMCOmegaPiPlPiMiPiZeroInAccPt = new TH1F*[fnCuts];

		fHistoTrueMotherPiPlPiMiPiZeroInvMassPt = new TH2F*[fnCuts];
		fHistoTrueMotherGammaGammaInvMassPt = new TH2F*[fnCuts];
		fHistoTrueMotherGammaGammaFromEtaInvMassPt = new TH2F*[fnCuts];
		fHistoTrueMotherGammaGammaFromOmegaInvMassPt = new TH2F*[fnCuts];
// 		if (fDoMesonQA){
			fHistoTruePionPionInvMassPt = 			new TH2F*[fnCuts];
			fHistoTruePionPionFromSameMotherInvMassPt = 	new TH2F*[fnCuts];
			fHistoTruePionPionFromEtaInvMassPt = 	new TH2F*[fnCuts];
			fHistoTruePionPionFromOmegaInvMassPt = 	new TH2F*[fnCuts];
// 		}
		
		for(Int_t iCut = 0; iCut<fnCuts;iCut++){
			TString cutstringEvent		= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
			TString cutstringPion		= ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutNumber();
			TString cutstringConvGamma = "";
			if (fNeutralPionMode < 2)  cutstringConvGamma = ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutNumber();
			TString cutstringCaloGamma = "";
			if (fNeutralPionMode > 0)  cutstringCaloGamma = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
			TString cutstringNeutralPion= ((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(iCut))->GetCutNumber();
			TString cutstringMeson		= ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
			
			TString fullCutString = "";
			if (fNeutralPionMode == 0) fullCutString = Form("%i_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
			else if (fNeutralPionMode == 1) fullCutString = Form("%i_%s_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringConvGamma.Data(),cutstringCaloGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
			else if (fNeutralPionMode == 2) fullCutString = Form("%i_%s_%s_%s_%s_%s",fNeutralPionMode,cutstringEvent.Data(),cutstringCaloGamma.Data(),cutstringNeutralPion.Data(), cutstringPion.Data(),cutstringMeson.Data());
			TString nameMCList = Form("%s MC histograms", fullCutString.Data());
			TString nameTrueRecList = Form("%s True histograms", fullCutString.Data());

			fMCList[iCut] = new TList();
			fMCList[iCut]->SetName(nameMCList.Data());
			fMCList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fMCList[iCut]);

			fHistoMCAllGammaPt[iCut] = new TH1F("MC_AllGamma_Pt","MC_AllGamma_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCAllGammaPt[iCut]);			
			if (fNeutralPionMode < 2){
				fHistoMCConvGammaPt[iCut] = new TH1F("MC_ConvGamma_Pt","MC_ConvGamma_Pt",250,0,25);
				fMCList[iCut]->Add(fHistoMCConvGammaPt[iCut]);						
			}
			
			fHistoMCAllPosPionsPt[iCut] = new TH1F("MC_AllPosPions_Pt","MC_AllPosPions_Pt",1000,0,25);
			fMCList[iCut]->Add(fHistoMCAllPosPionsPt[iCut]);
			fHistoMCAllNegPionsPt[iCut] = new TH1F("MC_AllNegPions_Pt","MC_AllNegPions_Pt",1000,0,25);
			fMCList[iCut]->Add(fHistoMCAllNegPionsPt[iCut]);
			fHistoMCGammaFromNeutralMesonPt[iCut] = new TH1F("MC_GammaFromNeutralMeson_Pt","MC_GammaFromNeutralMeson_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCGammaFromNeutralMesonPt[iCut]);	
			fHistoMCPosPionsFromNeutralMesonPt[iCut] = new TH1F("MC_PosPionsFromNeutralMeson_Pt","MC_PosPionsFromNeutralMeson_Pt",1000,0,25);
			fMCList[iCut]->Add(fHistoMCPosPionsFromNeutralMesonPt[iCut]);
			fHistoMCNegPionsFromNeutralMesonPt[iCut] = new TH1F("MC_NegPionsFromNeutralMeson_Pt","MC_NegPionsFromNeutralMeson_Pt",1000,0,25);
			fMCList[iCut]->Add(fHistoMCNegPionsFromNeutralMesonPt[iCut]);		

			fHistoMCEtaPiPlPiMiPiZeroPt[iCut] = new TH1F("MC_Eta_Pt","MC_Eta_Pt",250,0,25);
			fHistoMCEtaPiPlPiMiPiZeroPt[iCut]->Sumw2();
			fMCList[iCut]->Add(fHistoMCEtaPiPlPiMiPiZeroPt[iCut]);
			
			fHistoMCEtaPiPlPiMiPiZeroInAccPt[iCut] = new TH1F("MC_EtaInAcc_Pt","MC_EtaInAcc_Pt",250,0,25);
			fHistoMCEtaPiPlPiMiPiZeroInAccPt[iCut]->Sumw2();
			fMCList[iCut]->Add(fHistoMCEtaPiPlPiMiPiZeroInAccPt[iCut]);

			fHistoMCOmegaPiPlPiMiPiZeroPt[iCut] = new TH1F("MC_Omega_Pt","MC_Omega_Pt",250,0,25);
			fHistoMCOmegaPiPlPiMiPiZeroPt[iCut]->Sumw2();
			fMCList[iCut]->Add(fHistoMCOmegaPiPlPiMiPiZeroPt[iCut]);
			
			fHistoMCOmegaPiPlPiMiPiZeroInAccPt[iCut] = new TH1F("MC_OmegaInAcc_Pt","MC_OmegaInAcc_Pt",250,0,25);
			fHistoMCOmegaPiPlPiMiPiZeroInAccPt[iCut]->Sumw2();
			fMCList[iCut]->Add(fHistoMCOmegaPiPlPiMiPiZeroInAccPt[iCut]);

			fTrueList[iCut] = new TList();
			fTrueList[iCut]->SetName(nameTrueRecList.Data());
			fTrueList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fTrueList[iCut]);

			if (fNeutralPionMode < 2){
				fHistoTrueConvGammaPt[iCut] = new TH1F("ESD_TrueConvGamma_Pt","ESD_TrueConvGamma_Pt",250,0,25);
				fTrueList[iCut]->Add(fHistoTrueConvGammaPt[iCut]);
				fHistoTrueConvGammaFromNeutralMesonPt[iCut] = new TH1F("ESD_TrueConvGammaFromNeutralMeson_Pt","ESD_TrueConvGammaFromNeutralMeson_Pt",250,0,25);
				fTrueList[iCut]->Add(fHistoTrueConvGammaFromNeutralMesonPt[iCut]);
			}
			if (fNeutralPionMode > 0){
				fHistoTrueClusterGammaPt[iCut] = new TH1F("ESD_TrueClusterGamma_Pt","ESD_TrueClusterGamma_Pt",250,0,25);
				fTrueList[iCut]->Add(fHistoTrueClusterGammaPt[iCut]);
				fHistoTrueClusterGammaFromNeutralMesonPt[iCut] = new TH1F("ESD_TrueClusterGammaFromNeutralMeson_Pt","ESD_TrueClusterGammaFromNeutralMeson_Pt",250,0,25);
				fTrueList[iCut]->Add(fHistoTrueClusterGammaFromNeutralMesonPt[iCut]);
			}
			fHistoTruePosPionPt[iCut] = new TH1F("ESD_TruePosPion_Pt","ESD_TruePosPion_Pt",1000,0,25);
			fTrueList[iCut]->Add(fHistoTruePosPionPt[iCut]);
			fHistoTrueNegPionPt[iCut] = new TH1F("ESD_TrueNegPion_Pt","ESD_TrueNegPion_Pt",1000,0,25);
			fTrueList[iCut]->Add(fHistoTrueNegPionPt[iCut]);	

			fHistoTrueNegPionFromNeutralMesonPt[iCut] = new TH1F("ESD_TrueNegPionFromNeutralMeson_Pt","ESD_TrueNegPionFromNeutralMeson_Pt",1000,0,25);
			fTrueList[iCut]->Add(fHistoTrueNegPionFromNeutralMesonPt[iCut]);
			fHistoTruePosPionFromNeutralMesonPt[iCut] = new TH1F("ESD_TruePosPionFromNeutralMeson_Pt","ESD_TruePosPionFromNeutralMeson_Pt",1000,0,25);
			fTrueList[iCut]->Add(fHistoTruePosPionFromNeutralMesonPt[iCut]);

			fHistoTrueMotherPiPlPiMiPiZeroInvMassPt[iCut] = new TH2F("ESD_TrueMotherPiPlPiMiPiZero_InvMass_Pt","ESD_TrueMotherPiPlPiMiPiZero_InvMass_Pt",500,0.4,0.9,250,0,25);
			fHistoTrueMotherPiPlPiMiPiZeroInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiPiZeroInvMassPt[iCut]);
		
			fHistoTrueMotherGammaGammaInvMassPt[iCut] = new TH2F("ESD_TrueMotherGG_InvMass_Pt","ESD_TrueMotherGG_InvMass_Pt",450,0.,0.45,250,0,25);
			fTrueList[iCut]->Add(fHistoTrueMotherGammaGammaInvMassPt[iCut]);
			fHistoTrueMotherGammaGammaFromEtaInvMassPt[iCut] = new TH2F("ESD_TrueMotherGGFromEta_InvMass_Pt","ESD_TrueMotherGGFromEta_InvMass_Pt",450,0.,0.45,250,0,25);
			fTrueList[iCut]->Add(fHistoTrueMotherGammaGammaFromEtaInvMassPt[iCut]);
			fHistoTrueMotherGammaGammaFromOmegaInvMassPt[iCut] = new TH2F("ESD_TrueMotherGGFromOmega_InvMass_Pt","ESD_TrueMotherGGFromOmega_InvMass_Pt",450,0.,0.45,250,0,25);
			fTrueList[iCut]->Add(fHistoTrueMotherGammaGammaFromOmegaInvMassPt[iCut]);

			
// 			if (fDoMesonQA){
				fHistoTruePionPionInvMassPt[iCut] = new TH2F("ESD_TruePiPlusPiNeg_InvMassPt","ESD_TruePiPlusPiNeg_InvMassPt",2000,0.,2.,200,0.,20.);
				fTrueList[iCut]->Add(fHistoTruePionPionInvMassPt[iCut]);
				fHistoTruePionPionFromSameMotherInvMassPt[iCut] = new TH2F("ESD_TruePiPlusPiNegFromSameMother_InvMassPt","ESD_TruePiPlusPiNegFromSameMother_InvMassPt",2000,0.,2.,200,0.,20.);
				fTrueList[iCut]->Add(fHistoTruePionPionFromSameMotherInvMassPt[iCut]);
				fHistoTruePionPionFromEtaInvMassPt[iCut] = new TH2F("ESD_TruePiPlusPiNegFromEta_InvMassPt","ESD_TruePiPlusPiNegFromEta_InvMassPt",2000,0.,2.,200,0.,20.);
				fTrueList[iCut]->Add(fHistoTruePionPionFromEtaInvMassPt[iCut]);
				fHistoTruePionPionFromOmegaInvMassPt[iCut] = new TH2F("ESD_TruePiPlusPiNegFromOmega_InvMassPt","ESD_TruePiPlusPiNegFromOmega_InvMassPt",2000,0.,2.,200,0.,20.);
				fTrueList[iCut]->Add(fHistoTruePionPionFromOmegaInvMassPt[iCut]);

				// 			}
		}
	}

	

	InitBack(); // Init Background Handler

	fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
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
		if (fNeutralPionMode < 2){
			if( fGammaCutArray ) {
				if( ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutHistograms() ) {
					fCutFolder[iCut]->Add( ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutHistograms()  );
				}
			}
		} 
		if (fNeutralPionMode > 0){
			if( fClusterCutArray ) {
				if( ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms() ) {
					fCutFolder[iCut]->Add( ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms()  );
				}
			}			
		}	
		if( fNeutralPionMesonCutArray  ) {
			if( ((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(iCut))->GetCutHistograms() ) {
				fCutFolder[iCut]->Add( ((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(iCut))->GetCutHistograms());
			}
		}
		if( fMesonCutArray  ) {
			if( ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms() ) {
				fCutFolder[iCut]->Add( ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
			}
		}
	}

	PostData(1, fOutputContainer);

}

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::UserExec(Option_t *){

	//
	// Execute analysis for current event
	//

	fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
	if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader

	Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();

	if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1
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
		
		Bool_t isRunningEMCALrelAna = kFALSE;
		if (fNeutralPionMode > 0){
			if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;
		}	
		
		Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon, isRunningEMCALrelAna);
		
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

		if(fMCEvent){ // Process MC Particle
			fMCStack = fMCEvent->Stack();			
			if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection() != 0){
				((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(), 
																						((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
																						fMCEvent);
			} 
			ProcessMCParticles();
		}

		if (fNeutralPionMode < 2){
			ProcessConversionPhotonCandidates(); // Process this cuts conversion gammas
		}
		if (fNeutralPionMode > 0){
			ProcessCaloPhotonCandidates(); // Process this cuts calo gammas
		}
		
		if (fNeutralPionMode == 0 ){
			ProcessNeutralPionCandidatesPureConversions(); // Process neutral pion candidates purely from conversions
		}
		if (fNeutralPionMode == 1){
			ProcessNeutralPionCandidatesMixedConvCalo(); // Process neutral pion candidates mixed conv and calo
		}	
		if (fNeutralPionMode == 2){
			ProcessNeutralPionCandidatesPureCalo(); // Process neutral pion candidates purely from calo
		}	
			
		ProcessPionCandidates(); // Process this cuts gammas
			
		CalculateMesonCandidates();
		CalculateBackground();
		UpdateEventByEventData();
				
		
		if (fGoodConvGammas->GetEntries()>0) fGoodConvGammas->Clear(); 
		if (fClusterCandidates->GetEntries()>0) fClusterCandidates->Clear();
		if (fNeutralPionCandidates->GetEntries()>0) fNeutralPionCandidates->Clear();
		fGoodVirtualParticles->Clear(); // delete this cuts good gammas
	}

	fSelectorNegPionIndex.clear();
	fSelectorPosPionIndex.clear();

	PostData( 1, fOutputContainer );
}

Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::Notify(){
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		if( !((AliConvEventCuts*)fEventCutArray->At(iCut))->GetDoEtaShift() ){
			fProfileEtaShift[iCut]->Fill(0.,0.);
			continue; // No Eta Shift requested, continue
		}
		if( ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
			((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCorrectEtaShiftFromPeriod(fV0Reader->GetPeriodName());
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


void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::Terminate(const Option_t *){
///Grid
}


//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessCaloPhotonCandidates()
{
	
	Int_t nclus = 0;
	nclus = fInputEvent->GetNumberOfCaloClusters();
	
// 	cout << nclus << endl;
	
	if(nclus == 0)	return;
	
	// vertex
	Double_t vertex[3] = {0};
	InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
	
	// Loop over EMCal clusters
	for(Long_t i = 0; i < nclus; i++){
		
		AliVCluster* clus = NULL;
		clus = fInputEvent->GetCaloCluster(i);		
		if (!clus) continue;
		if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fIsMC)) continue;
		// TLorentzvector with cluster
		TLorentzVector clusterVector;
		clus->GetMomentum(clusterVector,vertex);
		
		TLorentzVector* tmpvec = new TLorentzVector();
		tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());
		
		// convert to AODConversionPhoton
		AliAODConversionPhoton *PhotonCandidate=new AliAODConversionPhoton(tmpvec);
		if(!PhotonCandidate) continue;
		
		// Flag Photon as CaloPhoton
		PhotonCandidate->SetIsCaloPhoton();
		PhotonCandidate->SetCaloClusterRef(i);
		// get MC label
		if(fIsMC){
			Int_t* mclabelsCluster = clus->GetLabels();
			PhotonCandidate->SetNCaloPhotonMCLabels(clus->GetNLabels());
// 			cout << clus->GetNLabels() << endl;
			if (clus->GetNLabels()>0){
				for (Int_t k =0; k< (Int_t)clus->GetNLabels(); k++){
					if (k< 20)PhotonCandidate->SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
// 					Int_t pdgCode = fMCStack->Particle(mclabelsCluster[k])->GetPdgCode();
// 					cout << "label " << k << "\t" << mclabelsCluster[k] << " pdg code: " << pdgCode << endl;
				}	
			}
		}
		
		fIsFromMBHeader = kTRUE; 
		// test whether largest contribution to cluster orginates in added signals
		if (fIsMC && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetCaloPhotonMCLabel(0), fMCStack, fInputEvent) == 0) fIsFromMBHeader = kFALSE;
		
		if (fIsFromMBHeader){
			fHistoClusterGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
			fHistoClusterGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
		}	
		fClusterCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas
		
		if(fIsMC){
// 			if(fInputEvent->IsA()==AliESDEvent::Class()){
				ProcessTrueCaloPhotonCandidates(PhotonCandidate);
// 			} else {
// 				ProcessTrueClusterCandidatesAOD(PhotonCandidate);
// 			}	
		}
		
		delete tmpvec;
	}
	
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessTrueCaloPhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
	TParticle *Photon = NULL;
	if (!TruePhotonCandidate->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set task will abort");
// 	fHistoTrueNLabelsInClus[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels());

	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();

	if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0)Photon = fMCStack->Particle(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
		else return;
		
	if(Photon == NULL){
	//    cout << "no photon" << endl;
		return;
	}

// 	Int_t pdgCodeParticle = Photon->GetPdgCode();
	TruePhotonCandidate->SetCaloPhotonMCFlags(fMCStack);
	
	// True Photon
	if(fIsFromMBHeader){	
		Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, TruePhotonCandidate->GetCaloPhotonMCLabel(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
		if(isPrimary){
			if (TruePhotonCandidate->IsLargestComponentPhoton()){
				fHistoTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				if (GammaIsNeutralMesonPiPlPiMiPiZeroDaughter(TruePhotonCandidate->GetCaloPhotonMCLabel(0))){
					fHistoTrueClusterGammaFromNeutralMesonPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				}
			}	
			if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
					fHistoTrueClusterGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
					if (GammaIsNeutralMesonPiPlPiMiPiZeroDaughter(TruePhotonCandidate->GetCaloPhotonMCLabel(0))){
					fHistoTrueClusterGammaFromNeutralMesonPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				}
			}
		}	
	}
	return;
}



//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessConversionPhotonCandidates(){
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
				= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack, fInputEvent);
			if(isPosFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
			Int_t isNegFromMBHeader
				= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack,fInputEvent);
			if(isNegFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
			if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
		}
		
		if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fESDEvent)) continue;

		if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut() &&
			!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){ // if no post reader loop is required add to events good gammas
			
			fGoodConvGammas->Add(PhotonCandidate);
		
			if(fIsFromMBHeader){
				fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
				fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
			}
		
			if(fMCEvent){
				ProcessTrueConversionPhotonCandidates(PhotonCandidate);
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
				= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack,fInputEvent);
				Int_t isNegFromMBHeader
				= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack,fInputEvent);
				if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
			}
			if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GoodGammasStepOne->GetEntries())) continue;
			if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
				fGoodConvGammas->Add(PhotonCandidate);
				if(fIsFromMBHeader){
					fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
					fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
				}
				if(fMCEvent){
					ProcessTrueConversionPhotonCandidates(PhotonCandidate);
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
				= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack,fInputEvent);
				Int_t isNegFromMBHeader
				= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack,fInputEvent);
				if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
			}
			
			if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GoodGammasStepTwo,i)) continue;
			fGoodConvGammas->Add(PhotonCandidate); // Add gamma to current cut TList
		
			if(fIsFromMBHeader){
				fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt()); // Differences to old V0Reader in p_t due to conversion KF->TLorentzVector
				fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
			}
		
			if(fMCEvent){
				ProcessTrueConversionPhotonCandidates(PhotonCandidate);
			}
		}
	}

	delete GoodGammasStepOne;
	GoodGammasStepOne = 0x0;
	delete GoodGammasStepTwo;
	GoodGammasStepTwo = 0x0;
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessTrueConversionPhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
	// Process True Photons
	AliStack *MCStack = fMCEvent->Stack();
	TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(MCStack);
	TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(MCStack);

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

	TParticle *Photon = TruePhotonCandidate->GetMCParticle(MCStack);
	if(Photon->GetPdgCode() != 22) return; // Mother is no Photon

	// True Photon
	
	Int_t labelGamma = TruePhotonCandidate->GetMCParticleLabel(MCStack);
	Bool_t gammaIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, labelGamma, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
	if( gammaIsPrimary ){
		if( fIsFromMBHeader ){
			fHistoTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (GammaIsNeutralMesonPiPlPiMiPiZeroDaughter(labelGamma)){
				fHistoTrueConvGammaFromNeutralMesonPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			}	
		}
	}
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessNeutralPionCandidatesPureConversions(){
	// Conversion Gammas
	if(fGoodConvGammas->GetEntries()>1){
		for(Int_t firstGammaIndex=0;firstGammaIndex<fGoodConvGammas->GetEntries()-1;firstGammaIndex++){
			AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGoodConvGammas->At(firstGammaIndex));
			if (gamma0==NULL) continue;
			for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGoodConvGammas->GetEntries();secondGammaIndex++){
				AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGoodConvGammas->At(secondGammaIndex));
				//Check for same Electron ID
				if (gamma1==NULL) continue;
				if(gamma0->GetTrackLabelPositive() == gamma1->GetTrackLabelPositive() ||
				gamma0->GetTrackLabelNegative() == gamma1->GetTrackLabelNegative() ||
				gamma0->GetTrackLabelNegative() == gamma1->GetTrackLabelPositive() ||
				gamma0->GetTrackLabelPositive() == gamma1->GetTrackLabelNegative() ) continue;

				AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);
				pi0cand->SetLabels(firstGammaIndex,secondGammaIndex);

				pi0cand->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
				if((((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
					fHistoGammaGammaInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
					if (pi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->GetSelectionLow() && pi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->GetSelectionHigh()){
						fNeutralPionCandidates->Add(pi0cand);
// 						cout << "Pi0 candidate " << pi0cand->M() << "\t" << pi0cand->Pt() << endl;
					}	
				
					if(fIsMC){
						if(fInputEvent->IsA()==AliESDEvent::Class())
							ProcessTrueNeutralPionCandidatesPureConversions(pi0cand,gamma0,gamma1);
						if(fInputEvent->IsA()==AliAODEvent::Class())
							ProcessTrueNeutralPionCandidatesPureConversionsAOD(pi0cand,gamma0,gamma1);
					}
				}
// 				delete pi0cand;
// 				pi0cand=0x0;
			}
		}
	}
}


//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessNeutralPionCandidatesPureCalo(){
	
	// Conversion Gammas
	if(fClusterCandidates->GetEntries()>0){

		// vertex
		Double_t vertex[3] = {0};
		InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

		for(Int_t firstGammaIndex=0;firstGammaIndex<fClusterCandidates->GetEntries();firstGammaIndex++){
			AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(firstGammaIndex));
			if (gamma0==NULL) continue;
			
			for(Int_t secondGammaIndex=0;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
				if (firstGammaIndex == secondGammaIndex) continue;
				AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
				if (gamma1==NULL) continue;
				
				AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);
				pi0cand->SetLabels(firstGammaIndex,secondGammaIndex);

				if((((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
					fHistoGammaGammaInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
					if (pi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->GetSelectionLow() && pi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->GetSelectionHigh()){
						fNeutralPionCandidates->Add(pi0cand);
// 						cout << "Pi0 candidate " << pi0cand->M() << "\t" << pi0cand->Pt() << endl;
					}	
				
					if(fIsMC){
// 						if(fInputEvent->IsA()==AliESDEvent::Class())
							ProcessTrueNeutralPionCandidatesPureCalo(pi0cand,gamma0,gamma1);
// 						if(fInputEvent->IsA()==AliAODEvent::Class())
// 							ProcessTrueNeutralPionCandidatesPureConversionsAOD(pi0cand,gamma0,gamma1);
					}
				}	
			}
		}
	}
}	

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessTrueNeutralPionCandidatesPureCalo( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
	// Process True Mesons
	AliStack *MCStack = fMCEvent->Stack();
	
	Bool_t isTruePi0 = kFALSE;
	Int_t gamma0MCLabel = TrueGammaCandidate0->GetCaloPhotonMCLabel(0); 	// get most probable MC label
	Int_t gamma0MotherLabel = -1;
	Int_t motherRealLabel = -1;
	
	if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
		TParticle * gammaMC0 = (TParticle*)MCStack->Particle(gamma0MCLabel);
		if (TrueGammaCandidate0->IsLargestComponentPhoton() || TrueGammaCandidate0->IsLargestComponentElectron()){		// largest component is electro magnetic
			// get mother of interest (pi0 or eta)
			if (TrueGammaCandidate0->IsLargestComponentPhoton()){														// for photons its the direct mother 
				gamma0MotherLabel=gammaMC0->GetMother(0);
				motherRealLabel=gammaMC0->GetFirstMother();
			} else if (TrueGammaCandidate0->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
				if (TrueGammaCandidate0->IsConversion()){
					gamma0MotherLabel=MCStack->Particle(gammaMC0->GetMother(0))->GetMother(0);
					motherRealLabel=MCStack->Particle(gammaMC0->GetMother(0))->GetMother(0);
				} else {
					gamma0MotherLabel=gammaMC0->GetMother(0); 
					motherRealLabel=gammaMC0->GetMother(0); 
				}
			}
		}
	}
	
	if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
	
	Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); 	// get most probable MC label
	Int_t gamma1MotherLabel = -1;
	// check if 
	if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
		// Daughters Gamma 1
		TParticle * gammaMC1 = (TParticle*)MCStack->Particle(gamma1MCLabel);
		if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){		// largest component is electro magnetic
			// get mother of interest (pi0 or eta)
			if (TrueGammaCandidate1->IsLargestComponentPhoton()){														// for photons its the direct mother 
				gamma1MotherLabel=gammaMC1->GetMother(0);
			} else if (TrueGammaCandidate1->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
				if (TrueGammaCandidate1->IsConversion()) gamma1MotherLabel=MCStack->Particle(gammaMC1->GetMother(0))->GetMother(0);
				else gamma1MotherLabel=gammaMC1->GetMother(0); 
			}
		} 	
	}
			
	if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
		if(((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
			isTruePi0=kTRUE;
		}
	}
	
	if(isTruePi0){// True Pion
		Pi0Candidate->SetTrueMesonValue(1);
		Pi0Candidate->SetMCLabel(motherRealLabel);
		fHistoTrueMotherGammaGammaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt()); 		
		if( IsEtaPiPlPiMiPiZeroDaughter(motherRealLabel) ) { 
			fHistoTrueMotherGammaGammaFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
		}
		if( IsOmegaPiPlPiMiPiZeroDaughter(motherRealLabel) ) { 
			fHistoTrueMotherGammaGammaFromOmegaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
		}
	}
}



//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessTrueNeutralPionCandidatesPureConversions(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
	// Process True Mesons
	AliStack *MCStack = fMCEvent->Stack();
	if(TrueGammaCandidate0->GetV0Index()<fInputEvent->GetNumberOfV0s()){
		Bool_t isTruePi0 = kFALSE;
		Bool_t isTruePi0Dalitz = kFALSE;
		Bool_t gamma0DalitzCand = kFALSE;
		Bool_t gamma1DalitzCand = kFALSE;
		Int_t gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(MCStack);
		Int_t gamma0MotherLabel = -1;
		Int_t motherRealLabel = -1;
		if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
			// Daughters Gamma 0
			TParticle * negativeMC = (TParticle*)TrueGammaCandidate0->GetNegativeMCDaughter(MCStack);
			TParticle * positiveMC = (TParticle*)TrueGammaCandidate0->GetPositiveMCDaughter(MCStack);
			TParticle * gammaMC0 = (TParticle*)MCStack->Particle(gamma0MCLabel);
			if(abs(negativeMC->GetPdgCode())==11 && abs(positiveMC->GetPdgCode())==11){  // Electrons ...
				if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
					if(gammaMC0->GetPdgCode() == 22){ // ... with Gamma Mother
						gamma0MotherLabel=gammaMC0->GetFirstMother();
						motherRealLabel=gammaMC0->GetFirstMother();
					}
				}
				if(gammaMC0->GetPdgCode() ==111){ // Dalitz candidate
					gamma0DalitzCand = kTRUE;
					gamma0MotherLabel=-111;
					motherRealLabel=gamma0MCLabel;
				}
			}
		}
		if(TrueGammaCandidate1->GetV0Index()<fInputEvent->GetNumberOfV0s()){
			Int_t gamma1MCLabel = TrueGammaCandidate1->GetMCParticleLabel(MCStack);
			Int_t gamma1MotherLabel = -1;
			if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
				// Daughters Gamma 1
				TParticle * negativeMC = (TParticle*)TrueGammaCandidate1->GetNegativeMCDaughter(MCStack);
				TParticle * positiveMC = (TParticle*)TrueGammaCandidate1->GetPositiveMCDaughter(MCStack);
				TParticle * gammaMC1 = (TParticle*)MCStack->Particle(gamma1MCLabel);
				if(abs(negativeMC->GetPdgCode())==11 && abs(positiveMC->GetPdgCode())==11){  // Electrons ...
					if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
						if(gammaMC1->GetPdgCode() == 22){ // ... with Gamma Mother
							gamma1MotherLabel=gammaMC1->GetFirstMother();
						}
					}
					if(gammaMC1->GetPdgCode() ==111 ){ // Dalitz candidate
						gamma1DalitzCand = kTRUE;
						gamma1MotherLabel=-111;
					}
				}
			}
			if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
				if(((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
					isTruePi0=kTRUE;
				}
			}
			
			//Identify Dalitz candidate
			if (gamma1DalitzCand || gamma0DalitzCand){
				if (gamma0DalitzCand && gamma0MCLabel >=0 && gamma0MCLabel==gamma1MotherLabel){
					if (gamma0MotherLabel == -111) isTruePi0Dalitz = kTRUE;
				}   
				if (gamma1DalitzCand && gamma1MCLabel >=0 && gamma1MCLabel==gamma0MotherLabel){
					if (gamma1MotherLabel == -111) isTruePi0Dalitz = kTRUE;
				}
			}
			
			
			if(isTruePi0 || isTruePi0Dalitz){// True Pion 
				Pi0Candidate->SetTrueMesonValue(1);
				Pi0Candidate->SetMCLabel(motherRealLabel);
				fHistoTrueMotherGammaGammaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt()); 
				if( IsEtaPiPlPiMiPiZeroDaughter(motherRealLabel) ) { 
					fHistoTrueMotherGammaGammaFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				}
				if( IsOmegaPiPlPiMiPiZeroDaughter(motherRealLabel) ) { 
					fHistoTrueMotherGammaGammaFromOmegaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				}	
			}
		}	
	}
}

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessTrueNeutralPionCandidatesPureConversionsAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{

	// Process True Mesons
	TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	Bool_t isTruePi0 = kFALSE;
	Bool_t isTruePi0Dalitz = kFALSE;
	Bool_t gamma0DalitzCand = kFALSE;
	Bool_t gamma1DalitzCand = kFALSE;
	Int_t motherRealLabel = -1;
		
	if (AODMCTrackArray!=NULL && TrueGammaCandidate0 != NULL){
		AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
		AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));

		Int_t gamma0MCLabel = -1;
		Int_t gamma0MotherLabel = -1;
		if(!positiveMC||!negativeMC)
			return;
		
		if(positiveMC->GetMother()>-1&&(negativeMC->GetMother() == positiveMC->GetMother())){
			gamma0MCLabel = positiveMC->GetMother();
		}

		if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
			// Daughters Gamma 0
			AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
			if(abs(negativeMC->GetPdgCode())==11 && abs(positiveMC->GetPdgCode())==11){  // Electrons ...
				if(((positiveMC->GetMCProcessCode())) == 5 && ((negativeMC->GetMCProcessCode())) == 5){ // ... From Conversion ...     
					if(gammaMC0->GetPdgCode() == 22){ // ... with Gamma Mother
						gamma0MotherLabel=gammaMC0->GetMother();
						motherRealLabel=gammaMC0->GetMother();
					}
				}
				if(gammaMC0->GetPdgCode() ==111){ // Dalitz candidate
					gamma0DalitzCand = kTRUE;
					gamma0MotherLabel=-111;
					motherRealLabel=gamma0MCLabel;
				}
			}
		}
		positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelPositive()));
		negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelNegative()));
		
		Int_t gamma1MCLabel = -1;
		Int_t gamma1MotherLabel = -1;
		if(!positiveMC||!negativeMC)
			return;
		
		if(positiveMC->GetMother()>-1&&(negativeMC->GetMother() == positiveMC->GetMother())){
			gamma1MCLabel = positiveMC->GetMother();
		}
		if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
			// Daughters Gamma 1
			AliAODMCParticle * gammaMC1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MCLabel));
			if(abs(negativeMC->GetPdgCode())==11 && abs(positiveMC->GetPdgCode())==11){  // Electrons ...
				if(((positiveMC->GetMCProcessCode())) == 5 && ((negativeMC->GetMCProcessCode())) == 5){ // ... From Conversion ...     
					if(gammaMC1->GetPdgCode() == 22){ // ... with Gamma Mother
					gamma1MotherLabel=gammaMC1->GetMother();
					}
				}
				if(gammaMC1->GetPdgCode() ==111 ){ // Dalitz candidate
						gamma1DalitzCand = kTRUE;
						gamma1MotherLabel=-111;
				}
			}
		}
		if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
			if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
				isTruePi0=kTRUE;
			}
		}
		
		//Identify Dalitz candidate
		if (gamma1DalitzCand || gamma0DalitzCand){
			if (gamma0DalitzCand && gamma0MCLabel >=0 && gamma0MCLabel==gamma1MotherLabel){
				if (gamma0MotherLabel == -111) isTruePi0Dalitz = kTRUE;
			}   
			if (gamma1DalitzCand && gamma1MCLabel >=0 && gamma1MCLabel==gamma0MotherLabel){
				if (gamma1MotherLabel == -111) isTruePi0Dalitz = kTRUE;   
			}
		}
					
		if(isTruePi0 || isTruePi0Dalitz){// True Pion 
			Pi0Candidate->SetTrueMesonValue(1);
			Pi0Candidate->SetMCLabel(motherRealLabel);
			fHistoTrueMotherGammaGammaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			if( IsEtaPiPlPiMiPiZeroDaughter(motherRealLabel) ) { 
				fHistoTrueMotherGammaGammaFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
			if( IsOmegaPiPlPiMiPiZeroDaughter(motherRealLabel) ) { 
				fHistoTrueMotherGammaGammaFromOmegaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
		}	
	}
	return;
}


//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessNeutralPionCandidatesMixedConvCalo(){
	
	// Conversion Gammas
	if(fGoodConvGammas->GetEntries()>0){
		// vertex
		Double_t vertex[3] = {0};
		InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

		for(Int_t firstGammaIndex=0;firstGammaIndex<fGoodConvGammas->GetEntries();firstGammaIndex++){
			AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGoodConvGammas->At(firstGammaIndex));
			if (gamma0==NULL) continue;
			
			for(Int_t secondGammaIndex=0;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
				Bool_t matched = kFALSE;
				AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
				if (gamma1==NULL) continue;
				
				if (gamma1->GetIsCaloPhoton()){
					AliVCluster* cluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
					matched = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gamma0,cluster, fInputEvent );
				}	
				
				AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);
				pi0cand->SetLabels(firstGammaIndex,secondGammaIndex);
				
				if((((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
					if (!matched){
						fHistoGammaGammaInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
						if (pi0cand->M() > ((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->GetSelectionLow() && pi0cand->M() < ((AliConversionMesonCuts*)fNeutralPionMesonCutArray->At(fiCut))->GetSelectionHigh()){
							fNeutralPionCandidates->Add(pi0cand);
	// 						cout << "Pi0 candidate " << pi0cand->M() << "\t" << pi0cand->Pt() << endl;
						}	
					
						if(fIsMC){
// 							if(fInputEvent->IsA()==AliESDEvent::Class())
								ProcessTrueNeutralPionCandidatesMixedConvCalo(pi0cand,gamma0,gamma1);
// 							if(fInputEvent->IsA()==AliAODEvent::Class())
// 								ProcessTrueMesonCandidatesAOD(pi0cand,gamma0,gamma1, matched);
						}
					}
				}
			}
		}
	}
}

//______________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessTrueNeutralPionCandidatesMixedConvCalo( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
	// Process True Mesons
	AliStack *MCStack = fMCEvent->Stack();
	if(TrueGammaCandidate0->GetV0Index()<fInputEvent->GetNumberOfV0s()){
		Bool_t isTruePi0 = kFALSE;
		Bool_t isTruePi0Dalitz = kFALSE;
		Bool_t gamma0DalitzCand = kFALSE;
		
		Int_t gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(MCStack);
		Int_t gamma0MotherLabel = -1;
		Int_t motherRealLabel = -1;
		if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
			// Daughters Gamma 0
			TParticle * negativeMC = (TParticle*)TrueGammaCandidate0->GetNegativeMCDaughter(MCStack);
			TParticle * positiveMC = (TParticle*)TrueGammaCandidate0->GetPositiveMCDaughter(MCStack);
			TParticle * gammaMC0 = (TParticle*)MCStack->Particle(gamma0MCLabel);
			if(abs(negativeMC->GetPdgCode())==11 && abs(positiveMC->GetPdgCode())==11){  // Electrons ...
				if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
					if(gammaMC0->GetPdgCode() == 22){ // ... with Gamma Mother
						gamma0MotherLabel=gammaMC0->GetFirstMother();
						motherRealLabel=gammaMC0->GetFirstMother();
					}
				}
				if(gammaMC0->GetPdgCode() ==111){ // Dalitz candidate
					gamma0DalitzCand = kTRUE;
					gamma0MotherLabel=-111;
					motherRealLabel=gamma0MCLabel;
				}

			}
		}
		
		if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
		
		Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); 	// get most probable MC label
		Int_t gamma1MotherLabel = -1;
		// check if 

		if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
			// Daughters Gamma 1
			TParticle * gammaMC1 = (TParticle*)MCStack->Particle(gamma1MCLabel);
			if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){		// largest component is electro magnetic
				// get mother of interest (pi0 or eta)
				if (TrueGammaCandidate1->IsLargestComponentPhoton()){														// for photons its the direct mother 
					gamma1MotherLabel=gammaMC1->GetMother(0);
				} else if (TrueGammaCandidate1->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
					if (TrueGammaCandidate1->IsConversion()) gamma1MotherLabel=MCStack->Particle(gammaMC1->GetMother(0))->GetMother(0);
					else gamma1MotherLabel=gammaMC1->GetMother(0); 
				}
			} 	
		}
				
		if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
			if(((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
				isTruePi0=kTRUE;
			}
		}
		
		if (gamma0DalitzCand ){
			if (gamma0DalitzCand && gamma0MCLabel >=0 && gamma0MCLabel==gamma1MotherLabel){
				if (gamma0MotherLabel == -111) isTruePi0Dalitz = kTRUE;
			}	
		}
			
		if(isTruePi0 || isTruePi0Dalitz ){
			Pi0Candidate->SetTrueMesonValue(1);
			Pi0Candidate->SetMCLabel(motherRealLabel);
			fHistoTrueMotherGammaGammaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt()); 		
			if( IsEtaPiPlPiMiPiZeroDaughter(motherRealLabel) ) { 
				fHistoTrueMotherGammaGammaFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
			if( IsOmegaPiPlPiMiPiZeroDaughter(motherRealLabel) ) { 
				fHistoTrueMotherGammaGammaFromOmegaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
		}
	}
}



//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessPionCandidates(){

	Double_t magField = fInputEvent->GetMagneticField();
	if( magField  < 0.0 ){
		magField =  1.0;
	} else {
		magField =  -1.0;
	}

	vector<Int_t> lGoodNegPionIndexPrev(0);
	vector<Int_t> lGoodPosPionIndexPrev(0);
	
	for(UInt_t i = 0; i < fSelectorNegPionIndex.size(); i++){
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
			Bool_t negPionIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, labelNegPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
			if( labelNegPion < fMCStack->GetNtrack() ){
				TParticle* negPion = fMCStack->Particle(labelNegPion);
				if( negPion->GetPdgCode() ==  -211 ){
					if( negPionIsPrimary ){
						fHistoTrueNegPionPt[fiCut]->Fill(negPionCandidate->Pt());    //primary negPion
					}		
					if( IsEtaPiPlPiMiPiZeroDaughter(labelNegPion) || IsOmegaPiPlPiMiPiZeroDaughter(labelNegPion) ) {
						if( negPionIsPrimary ) {
							fHistoTrueNegPionFromNeutralMesonPt[fiCut]->Fill(negPionCandidate->Pt());
						} 
					}	
				}
			}
		}
	}

	for(UInt_t i = 0; i < fSelectorPosPionIndex.size(); i++){
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
			Bool_t posPionIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, labelPosPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
			if( labelPosPion < fMCStack->GetNtrack() ) {
				TParticle* posPion = fMCStack->Particle(labelPosPion);
				if( posPion->GetPdgCode() ==  211 ){
					if( posPionIsPrimary ){
						fHistoTruePosPionPt[fiCut]->Fill(posPionCandidate->Pt());
					} 
					if( IsEtaPiPlPiMiPiZeroDaughter(labelPosPion) || IsOmegaPiPlPiMiPiZeroDaughter(labelPosPion) ) {
						if(posPionIsPrimary){
							fHistoTruePosPionFromNeutralMesonPt[fiCut]->Fill(posPionCandidate->Pt());
						} 
					}
				}
			}
		}
	}


	for(UInt_t i = 0; i < lGoodNegPionIndexPrev.size(); i++){

		AliESDtrack *negPionCandidate = fESDEvent->GetTrack(lGoodNegPionIndexPrev[i]);
		AliKFParticle negPionCandidateKF( *negPionCandidate->GetConstrainedParam(), 211 );

		for(UInt_t j = 0; j < lGoodPosPionIndexPrev.size(); j++){

			AliESDtrack *posPionCandidate = fESDEvent->GetTrack(lGoodPosPionIndexPrev[j]);
			AliKFParticle posPionCandidateKF( *posPionCandidate->GetConstrainedParam(), 211 );

			AliKFConversionPhoton* virtualPhoton = NULL;
			virtualPhoton = new AliKFConversionPhoton(negPionCandidateKF,posPionCandidateKF);
			AliKFVertex primaryVertexImproved(*fInputEvent->GetPrimaryVertex());
// 			primaryVertexImproved+=*virtualPhoton;
			virtualPhoton->SetProductionVertex(primaryVertexImproved);
			virtualPhoton->SetTrackLabels( lGoodPosPionIndexPrev[j], lGoodNegPionIndexPrev[i]);
			
			Int_t labeln=0;
			Int_t labelp=0;
			Int_t motherlabelp = 0;
			Int_t motherlabeln = 0;
			TParticle *fNegativeMCParticle =NULL;
			TParticle *fPositiveMCParticle =NULL;
			if( fMCEvent ) {
				labeln=TMath::Abs(negPionCandidate->GetLabel());
				labelp=TMath::Abs(posPionCandidate->GetLabel());
				fNegativeMCParticle = fMCStack->Particle(labeln);
				fPositiveMCParticle = fMCStack->Particle(labelp);
				// check whether MC particles exist, else abort
				if (fNegativeMCParticle == NULL || fPositiveMCParticle == NULL) return;

				motherlabeln = fNegativeMCParticle->GetMother(0);
				motherlabelp = fPositiveMCParticle->GetMother(0);
				virtualPhoton->SetMCLabelPositive(labelp);
				virtualPhoton->SetMCLabelNegative(labeln);
			
			}
			
			AliAODConversionPhoton *vParticle = new AliAODConversionPhoton(virtualPhoton); //To apply mass 2 pion mass cut
			if (((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->DoMassCut()){
				if (vParticle->GetMass() < ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetMassCut()){
					fGoodVirtualParticles->Add(  vParticle );
					fHistoPionPionInvMassPt[fiCut]->Fill( vParticle->GetMass(),vParticle->Pt());
				}
			} else {
				fGoodVirtualParticles->Add(  vParticle );
				fHistoPionPionInvMassPt[fiCut]->Fill( vParticle->GetMass(),vParticle->Pt());
			}	

			Double_t clsToFPos = -1.0;
			Double_t clsToFNeg = -1.0;
			
			Float_t dcaToVertexXYPos = -1.0;
			Float_t dcaToVertexZPos  = -1.0;
			Float_t dcaToVertexXYNeg = -1.0;
			Float_t dcaToVertexZNeg  = -1.0;
			
			
			
			if ( fDoMesonQA ) {	
				clsToFPos = ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetNFindableClustersTPC(posPionCandidate);
				clsToFNeg = ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetNFindableClustersTPC(negPionCandidate);
				
				Float_t bPos[2];
				Float_t bCovPos[3];
				posPionCandidate->GetImpactParameters(bPos,bCovPos);
				if (bCovPos[0]<=0 || bCovPos[2]<=0) {
					AliDebug(1, "Estimated b resolution lower or equal zero!");
					bCovPos[0]=0; bCovPos[2]=0;
				}
				
				Float_t bNeg[2];
				Float_t bCovNeg[3];
				posPionCandidate->GetImpactParameters(bNeg,bCovNeg);
				if (bCovNeg[0]<=0 || bCovNeg[2]<=0) {
					AliDebug(1, "Estimated b resolution lower or equal zero!");
					bCovNeg[0]=0; bCovNeg[2]=0;
				}
				
				dcaToVertexXYPos = bPos[0];
				dcaToVertexZPos  = bPos[1];
				dcaToVertexXYNeg = bNeg[0];
				dcaToVertexZNeg  = bNeg[1];
				
				fHistoNegPionEta[fiCut]->Fill( negPionCandidate->Eta() );
				fHistoPosPionEta[fiCut]->Fill( posPionCandidate->Eta() );
						
				fHistoNegPionClsTPC[fiCut]->Fill(clsToFNeg,negPionCandidate->Pt());
				fHistoPosPionClsTPC[fiCut]->Fill(clsToFPos,posPionCandidate->Pt());
				
				fHistoPionDCAxy[fiCut]->Fill(  dcaToVertexXYNeg, negPionCandidate->Pt() );
				fHistoPionDCAz[fiCut]->Fill(   dcaToVertexZNeg,  negPionCandidate->Pt() );
				fHistoPionDCAxy[fiCut]->Fill(  dcaToVertexXYPos, posPionCandidate->Pt() );
				fHistoPionDCAz[fiCut]->Fill(   dcaToVertexZPos,  posPionCandidate->Pt() );
				
				fHistoPionTPCdEdxNSigma[fiCut]->Fill( posPionCandidate->P(),((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetPIDResponse()->NumberOfSigmasTPC(posPionCandidate, AliPID::kPion) );
				fHistoPionTPCdEdxNSigma[fiCut]->Fill( negPionCandidate->P(),((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetPIDResponse()->NumberOfSigmasTPC(negPionCandidate, AliPID::kPion) );
				
				fHistoPionTPCdEdx[fiCut]->Fill( posPionCandidate->P(), TMath::Abs(posPionCandidate->GetTPCsignal()));
				fHistoPionTPCdEdx[fiCut]->Fill( negPionCandidate->P(), TMath::Abs(negPionCandidate->GetTPCsignal()));
			}

			if (fMCEvent){
				if (fPositiveMCParticle && fNegativeMCParticle ) {
					if (((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->DoMassCut()){
						if (vParticle->GetMass() < ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetMassCut()){
							if(TMath::Abs(fNegativeMCParticle->GetPdgCode())==211 && TMath::Abs(fPositiveMCParticle->GetPdgCode())==211){  // Pions ...
								fHistoTruePionPionInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt());
								if (motherlabeln == motherlabelp){
									fHistoTruePionPionFromSameMotherInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt());
									if( IsEtaPiPlPiMiPiZeroDaughter(labeln) ) { //|| IsOmegaPiPlPiMiPiZeroDaughter(labeln) 
										fHistoTruePionPionFromEtaInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt());
									}
									if( IsOmegaPiPlPiMiPiZeroDaughter(labeln) ) { //||  
										fHistoTruePionPionFromOmegaInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt());
									}
								}
							}
						}
					} else { 
						if(TMath::Abs(fNegativeMCParticle->GetPdgCode())==211 && TMath::Abs(fPositiveMCParticle->GetPdgCode())==211){  // Pions ...
							fHistoTruePionPionInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt());
							if (motherlabeln == motherlabelp){
								fHistoTruePionPionFromSameMotherInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt());
								if( IsEtaPiPlPiMiPiZeroDaughter(labeln) ) { //|| IsOmegaPiPlPiMiPiZeroDaughter(labeln) 
									fHistoTruePionPionFromEtaInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt());
								}
								if( IsOmegaPiPlPiMiPiZeroDaughter(labeln) ) { //||  
									fHistoTruePionPionFromOmegaInvMassPt[fiCut]->Fill(vParticle->GetMass(),vParticle->Pt());
								}
							}
						}
					}	
				}	
			}
						
			delete virtualPhoton;
			virtualPhoton=NULL;
					
		}
	}
}

//_____________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessMCParticles(){

	// Loop over all primary MC particle
	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();

	for(Int_t i = 0; i < fMCStack->GetNtrack(); i++) {
		if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){ 
		
			TParticle* particle = (TParticle *)fMCStack->Particle(i);
			if (!particle) continue;

			Int_t isMCFromMBHeader = -1;
			if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
				isMCFromMBHeader
					= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent);
				if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
			}

			if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack,fInputEvent)){
				// find MC photons 
				if (fNeutralPionMode < 2){
					if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCStack,kFALSE)){
						fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt()); // All MC Gamma
						if(particle->GetMother(0) >-1){
							if (fMCStack->Particle(particle->GetMother(0))->GetPdgCode() ==111){
								if (fMCStack->Particle(particle->GetMother(0))->GetMother(0) > -1){
									if ( fMCStack->Particle((fMCStack->Particle(particle->GetMother(0)))->GetMother(0))->GetPdgCode() == 221 || 
										fMCStack->Particle((fMCStack->Particle(particle->GetMother(0)))->GetMother(0))->GetPdgCode() == 223 ){ 
										if ( fMCStack->Particle(particle->GetMother(0))->GetNDaughters()==3 ) 
											fHistoMCGammaFromNeutralMesonPt[fiCut]->Fill(particle->Pt()); // All photons from eta or omega via pi0 
									}		
								}		
							}		
						}	
					}
				} else if (fNeutralPionMode == 2){
					if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(particle,fMCStack)){
						fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt()); // All MC Gamma
						if(particle->GetMother(0) >-1){
							if (fMCStack->Particle(particle->GetMother(0))->GetPdgCode() ==111){
								if (fMCStack->Particle(particle->GetMother(0))->GetMother(0) > -1){
									if ( fMCStack->Particle((fMCStack->Particle(particle->GetMother(0)))->GetMother(0))->GetPdgCode() == 221 || 
										fMCStack->Particle((fMCStack->Particle(particle->GetMother(0)))->GetMother(0))->GetPdgCode() == 223 ){ 
										if ( fMCStack->Particle(particle->GetMother(0))->GetNDaughters()==3 ) 
											fHistoMCGammaFromNeutralMesonPt[fiCut]->Fill(particle->Pt()); // All photons from eta or omega via pi0 
									}		
								}		
							}		
						}	
					}
				}	
				if (fNeutralPionMode < 2){
					if (((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCStack,kTRUE)){
						fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt());
					} // Converted MC Gamma
				}
				if(((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(i,fMCStack)){
					if( particle->GetPdgCode() == 211){
						fHistoMCAllPosPionsPt[fiCut]->Fill(particle->Pt()); // All pos pions
						if(particle->GetMother(0) >-1){
							if (fMCStack->Particle(particle->GetMother(0))->GetPdgCode() ==221 || fMCStack->Particle(particle->GetMother(0))->GetPdgCode() ==223) 
								fHistoMCPosPionsFromNeutralMesonPt[fiCut]->Fill(particle->Pt()); // All pos from eta or omega
						}	
					}	
					if( particle->GetPdgCode() == -211){
						fHistoMCAllNegPionsPt[fiCut]->Fill(particle->Pt()); // All neg pions
						if(particle->GetMother(0) >-1){
							if (fMCStack->Particle(particle->GetMother(0))->GetPdgCode() ==221 || fMCStack->Particle(particle->GetMother(0))->GetPdgCode() ==223 )
								fHistoMCNegPionsFromNeutralMesonPt[fiCut]->Fill(particle->Pt()); // All pos from eta or omega
						}	
					}
				}
				
							
				// \eta -> pi+ pi- \gamma 
				Int_t labelNeutPion = -1;
				Int_t labelNegPion = -1;
				Int_t labelPosPion = -1;

				if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMCPiPlPiMiPiZero(particle,fMCStack,labelNegPion,labelPosPion,labelNeutPion,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){	
					Float_t weighted= 1;
					if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) { 
						if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack,fInputEvent)){
							if (particle->Pt()>0.005){
								weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, fMCStack,fInputEvent);
							}
						}
					}
					if(particle->GetPdgCode() == 221)fHistoMCEtaPiPlPiMiPiZeroPt[fiCut]->Fill(particle->Pt(), weighted); 						// All MC Eta in respective decay channel
					if(particle->GetPdgCode() == 223)fHistoMCOmegaPiPlPiMiPiZeroPt[fiCut]->Fill(particle->Pt(), weighted); 						// All MC Omega in respective decay channel
					
					TParticle *neutPion    = fMCStack->Particle(labelNeutPion);
					TParticle *gamma1 = fMCStack->Particle(neutPion->GetDaughter(0));
					TParticle *gamma2 = fMCStack->Particle(neutPion->GetDaughter(1));
					Bool_t kDaughter0IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, neutPion->GetDaughter(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
					Bool_t kDaughter1IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, neutPion->GetDaughter(1), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
					Bool_t kNegPionIsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, labelNegPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
					Bool_t kPosPionIsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, labelPosPion, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
					
					if (fNeutralPionMode < 2){					
						if( kDaughter0IsPrim && kDaughter1IsPrim && kNegPionIsPrim && kPosPionIsPrim &&
							((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(gamma1,fMCStack,kFALSE) &&					// test first daugther of pi0
							((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(gamma2,fMCStack,kFALSE) &&					// test second daughter of pi0
							((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelNegPion,fMCStack) &&								// test negative pion
							((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelPosPion,fMCStack) 								// test positive pion
						) {
								if(particle->GetPdgCode() == 221) fHistoMCEtaPiPlPiMiPiZeroInAccPt[fiCut]->Fill(particle->Pt(), weighted ); 		// MC Eta pi+ pi- pi0 with gamma's and e+e- in acc
								if(particle->GetPdgCode() == 223) fHistoMCOmegaPiPlPiMiPiZeroInAccPt[fiCut]->Fill(particle->Pt(), weighted ); 		// MC Omega pi+ pi- pi0 with gamma's and e+e- in acc
						}				
					} else if (fNeutralPionMode == 2){
						if( kDaughter0IsPrim && kDaughter1IsPrim && kNegPionIsPrim && kPosPionIsPrim &&
							((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma1,fMCStack) &&					// test first daugther of pi0
							((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma2,fMCStack) &&					// test second daughter of pi0
							((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelNegPion,fMCStack) &&								// test negative pion
							((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelPosPion,fMCStack) 								// test positive pion
						) {
								if(particle->GetPdgCode() == 221) fHistoMCEtaPiPlPiMiPiZeroInAccPt[fiCut]->Fill(particle->Pt(), weighted ); 		// MC Eta pi+ pi- pi0 with gamma's and e+e- in acc
								if(particle->GetPdgCode() == 223) fHistoMCOmegaPiPlPiMiPiZeroInAccPt[fiCut]->Fill(particle->Pt(), weighted ); 		// MC Omega pi+ pi- pi0 with gamma's and e+e- in acc
						}
					}	
						
				}
			}
		}	
	}
}


//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::CalculateMesonCandidates(){

	// Conversion Gammas
	if( fNeutralPionCandidates->GetEntries() > 0 && fGoodVirtualParticles->GetEntries() > 0 ){
		
		for(Int_t mesonIndex=0; mesonIndex<fNeutralPionCandidates->GetEntries(); mesonIndex++){
			AliAODConversionMother *neutralPion=dynamic_cast<AliAODConversionMother*>(fNeutralPionCandidates->At(mesonIndex));
			if (neutralPion==NULL) continue;
			for(Int_t virtualParticleIndex=0;virtualParticleIndex<fGoodVirtualParticles->GetEntries();virtualParticleIndex++){

				AliAODConversionPhoton *vParticle=dynamic_cast<AliAODConversionPhoton*>(fGoodVirtualParticles->At(virtualParticleIndex));
				if (vParticle==NULL) continue;
				//Check for same Electron ID

				AliAODConversionMother *mesoncand = new AliAODConversionMother(neutralPion,vParticle);
				mesoncand->SetLabels(mesonIndex,virtualParticleIndex);

				if( ( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(mesoncand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())) ){
			
// 					cout<< "Meson Accepted "<<endl;
					
					Int_t zbin= fBGHandler[fiCut]->GetZBinIndex(fESDEvent->GetPrimaryVertex()->GetZ());
					Int_t mbin = 0;
					if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
						mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
					} else {
						if (fNeutralPionMode < 2) mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGoodConvGammas->GetEntries());
						else mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fClusterCandidates->GetEntries());
					}
										
					fHistoMotherInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt());
					Double_t sparesFill[4] = {mesoncand->M(),mesoncand->Pt(),(Double_t)zbin,(Double_t)mbin};
					fTHnSparseMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
									
					if(fMCEvent){
						ProcessTrueMesonCandidates(mesoncand,neutralPion,vParticle);
					}
				}
				delete mesoncand;
				mesoncand=0x0;
			}
		}
	}
}

//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::CalculateBackground(){

	Int_t zbin= fBGHandler[fiCut]->GetZBinIndex(fESDEvent->GetPrimaryVertex()->GetZ());
	Int_t mbin = 0;


	if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
		mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
	} else {
		if (fNeutralPionMode < 2) mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGoodConvGammas->GetEntries());
		else mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fClusterCandidates->GetEntries());
	}

	Int_t method = 1;
	AliGammaConversionAODBGHandler::GammaConversionVertex *bgEventVertex = NULL;
	if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity() ) {
		for(Int_t nEventsInBG=0;nEventsInBG<fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
			AliGammaConversionMotherAODVector *previousEventMesons = fBGHandler[fiCut]->GetBGGoodMesons(zbin,mbin,nEventsInBG);
			if(fMoveParticleAccordingToVertex == kTRUE && method == 1){
				bgEventVertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
			}

			for(Int_t iCurrent=0;iCurrent<fGoodVirtualParticles->GetEntries();iCurrent++){
				AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGoodVirtualParticles->At(iCurrent));

				for(UInt_t iPrevious=0;iPrevious<previousEventMesons->size();iPrevious++){
					AliAODConversionMother previousGoodMeson = (AliAODConversionMother)(*(previousEventMesons->at(iPrevious)));

					if(fMoveParticleAccordingToVertex == kTRUE && method == 1 ){
						MoveParticleAccordingToVertex(&previousGoodMeson,bgEventVertex);
					}

					AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&previousGoodMeson,&currentEventGoodV0);
								

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
			AliGammaConversionMotherAODVector *previousEventMesons = fBGHandler[fiCut]->GetBGGoodMesons(zbin,mbin,nEventsInBG);
			if(previousEventMesons){
				if(fMoveParticleAccordingToVertex == kTRUE && method == 1){
					bgEventVertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
				}
				for(Int_t iCurrent=0;iCurrent<fGoodVirtualParticles->GetEntries();iCurrent++){

					AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGoodVirtualParticles->At(iCurrent));

					for(UInt_t iPrevious=0;iPrevious<previousEventMesons->size();iPrevious++){

						AliAODConversionMother previousGoodMeson = (AliAODConversionMother)(*(previousEventMesons->at(iPrevious)));

						if(fMoveParticleAccordingToVertex == kTRUE && method ==1){
							MoveParticleAccordingToVertex(&previousGoodMeson,bgEventVertex);
						}

						AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&previousGoodMeson,&currentEventGoodV0);
								
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
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::ProcessTrueMesonCandidates(AliAODConversionMother *mesoncand, AliAODConversionMother *TrueNeutralPionCandidate, AliAODConversionPhoton *TrueVirtualParticleCandidate){

	// Process True Mesons
	AliStack *MCStack = fMCEvent->Stack();
	
	Bool_t isTrueEta = kFALSE;
	Bool_t isTrueOmega = kFALSE;
	Int_t trueMesonFlag = TrueNeutralPionCandidate->GetTrueMesonValue();
	Int_t pi0MCLabel= TrueNeutralPionCandidate->GetMCLabel();

	
	if ( !(trueMesonFlag == 1 && pi0MCLabel != -1)) return;
// 	cout << trueMesonFlag << "\t" << pi0MCLabel << endl;

 
	Int_t virtualParticleMCLabel = TrueVirtualParticleCandidate->GetMCParticleLabel(MCStack);
	Int_t virtualParticleMotherLabel = -1;
	Bool_t isPiPiDecay = kFALSE;
	
// 	if (fDoMesonQA){
		TParticle * negativeMC = (TParticle*)TrueVirtualParticleCandidate->GetNegativeMCDaughter(MCStack);
		TParticle * positiveMC = (TParticle*)TrueVirtualParticleCandidate->GetPositiveMCDaughter(MCStack);
// 	}
	
	if(virtualParticleMCLabel != -1){ // if virtualParticleMCLabel==-1 particles don't have same mother 
// 		TParticle * negativeMC = (TParticle*)TrueVirtualParticleCandidate->GetNegativeMCDaughter(MCStack);
// 		TParticle * positiveMC = (TParticle*)TrueVirtualParticleCandidate->GetPositiveMCDaughter(MCStack);
// 		TParticle * virtualParticleMotherMC = (TParticle*)MCStack->Particle(virtualParticleMCLabel);
// 		cout << "pdg code same mother - " << virtualParticleMotherMC->GetPdgCode() << endl;
		
		if(TMath::Abs(negativeMC->GetPdgCode())==211 && TMath::Abs(positiveMC->GetPdgCode())==211){  // Pions ...
			virtualParticleMotherLabel=virtualParticleMCLabel;
			isPiPiDecay=kTRUE;
// 			} else if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
// 				if( virtualParticleMotherMC->GetPdgCode() != 22 ){
// 					virtualParticleMotherLabel=virtualParticleMCLabel;
// 					isDalitz = kTRUE;
// 				} else if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
// 					virtualParticleMotherLabel=virtualParticleMotherMC->GetFirstMother();
// 					isRealGamma = kTRUE; //no virtual gamma
// 				}
		}	
	}
	if (IsEtaPiPlPiMiPiZeroDaughter(pi0MCLabel) || IsOmegaPiPlPiMiPiZeroDaughter(pi0MCLabel)){
		Int_t pi0MotherMCLabel = ((TParticle*)MCStack->Particle(pi0MCLabel))->GetMother(0);		
		if(virtualParticleMCLabel != -1){
// 			cout << "pi+pi- mother: "<<  virtualParticleMCLabel << endl;
// 			cout << "pi0 mother: "<<  pi0MotherMCLabel << endl;

// 			TParticle * virtualParticleMotherMC = (TParticle*)MCStack->Particle(virtualParticleMCLabel);
// 			cout << "pdg code same mother - " << virtualParticleMotherMC->GetPdgCode() << endl;		
		}	
		if( pi0MotherMCLabel == virtualParticleMotherLabel ){
			if(((TParticle*)MCStack->Particle(virtualParticleMotherLabel))->GetPdgCode() == 221){
// 				cout << "found eta" << endl;
				isTrueEta=kTRUE;
			}
			if(((TParticle*)MCStack->Particle(virtualParticleMotherLabel))->GetPdgCode() == 223){
// 				cout << "found omega" << endl;
				isTrueOmega=kTRUE;
			}
		}
	}
	

	if( isTrueEta || isTrueOmega ){ // True Eta or Omega
		if ( isPiPiDecay) { //real eta -> Pi+ Pi- Pi0
			Float_t weighted= 1;
// 			if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) { 
// 				if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gammaMotherLabel, fMCStack,fInputEvent)){
// 					if (((TParticle*)MCStack->Particle(gammaMotherLabel))->Pt()>0.005){
// 						weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),gammaMotherLabel,fMCStack,fInputEvent);
// 					}
// 				}
// 			}
			fHistoTrueMotherPiPlPiMiPiZeroInvMassPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(),weighted);
		}	
	}

}


//________________________________________________________________________
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::UpdateEventByEventData(){
   //see header file for documentation

	Int_t method = 1;
	if( method == 1 ) {
		if(fNeutralPionCandidates->GetEntries() >0 ){
			if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
				fBGHandler[fiCut]->AddMesonEvent(fNeutralPionCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),0);
			} else { // means we use #V0s for multiplicity
				if (fNeutralPionMode < 2) fBGHandler[fiCut]->AddMesonEvent(fNeutralPionCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGoodConvGammas->GetEntries(),0);
				else fBGHandler[fiCut]->AddMesonEvent(fNeutralPionCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fClusterCandidates->GetEntries(),0);
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
void AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::MoveParticleAccordingToVertex(AliAODConversionMother* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex){
   //see header file for documentation

   Double_t dx = vertex->fX - fInputEvent->GetPrimaryVertex()->GetX();
   Double_t dy = vertex->fY - fInputEvent->GetPrimaryVertex()->GetY();
   Double_t dz = vertex->fZ - fInputEvent->GetPrimaryVertex()->GetZ();

   Double_t movedPlace[3] = {particle->GetProductionX() - dx,particle->GetProductionY() - dy,particle->GetProductionZ() - dz};
   particle->SetProductionPoint(movedPlace);
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::IsEtaPiPlPiMiPiZeroDaughter( Int_t label ) const {
//
// Returns true if the particle comes from eta -> pi+ pi- gamma
//
	Int_t motherLabel = fMCStack->Particle( label )->GetMother(0);
	if( motherLabel < 0 || motherLabel >= fMCStack->GetNtrack() ) return kFALSE;
	
	TParticle* mother = fMCStack->Particle( motherLabel );
// 	cout << "found eta? " << endl;
	if( mother->GetPdgCode() != 221 ) return kFALSE;
// 		else cout << "YES" << endl;
	if( IsPiPlPiMiPiZeroDecay( mother ) ) return kTRUE;	
	return kFALSE;       
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::IsOmegaPiPlPiMiPiZeroDaughter( Int_t label ) const {
//
// Returns true if the particle comes from eta -> pi+ pi- gamma
//
	Int_t motherLabel = fMCStack->Particle( label )->GetMother(0);
	if( motherLabel < 0 || motherLabel >= fMCStack->GetNtrack() ) return kFALSE;
	
	TParticle* mother = fMCStack->Particle( motherLabel );
// 	cout << "found omega? " << endl;
	if( mother->GetPdgCode() != 223 ) return kFALSE;
// 		else cout << "YES" << endl;
	if( IsPiPlPiMiPiZeroDecay( mother ) ) return kTRUE;	
	return kFALSE;       
}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::IsPiPlPiMiPiZeroDecay(TParticle *fMCMother) const
{
// 	cout << fMCMother->GetNDaughters() << endl;
	if( fMCMother->GetNDaughters() != 3 ) return kFALSE;
// 	cout << fMCMother->GetPdgCode() << endl;
	if( !(fMCMother->GetPdgCode() == 221 || fMCMother->GetPdgCode() == 223)  ) return kFALSE;
// 	cout << "made it til here" << endl;
	
	TParticle *posPion = 0x0;
	TParticle *negPion = 0x0;
	TParticle *neutPion    = 0x0;
	
	for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){				
		TParticle* temp = (TParticle*)fMCStack->Particle( index );
		
		switch( temp->GetPdgCode() ) {
		case 211:
			posPion =  temp;
			break;
		case -211:
			negPion =  temp;
			break;
		case 111:
			neutPion = temp;
			break;
		}
	}  
	if( posPion && negPion && neutPion) return kTRUE;
	
	return kFALSE;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskNeutralMesonToPiPlPiMiPiZero::GammaIsNeutralMesonPiPlPiMiPiZeroDaughter( Int_t label ) const {
//
// Returns true if the particle comes from eta -> pi+ pi- gamma
//
	Int_t motherLabel = fMCStack->Particle( label )->GetMother(0);
	if( motherLabel < 0 || motherLabel >= fMCStack->GetNtrack() ) return kFALSE;
	
	TParticle* mother = fMCStack->Particle( motherLabel );
// 	cout << "found omega? " << endl;
	if( mother->GetPdgCode() != 111 ) return kFALSE;
// 		else cout << "YES" << endl;
	Int_t grandMotherLabel = mother->GetMother(0);
	if( grandMotherLabel < 0 || grandMotherLabel >= fMCStack->GetNtrack() ) return kFALSE;
	TParticle* grandmother = fMCStack->Particle( grandMotherLabel );
	
	if( IsPiPlPiMiPiZeroDecay( grandmother ) ) return kTRUE;	
	return kFALSE;       
}
