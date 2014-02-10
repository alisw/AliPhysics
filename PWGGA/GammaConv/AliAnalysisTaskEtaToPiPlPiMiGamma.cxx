/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Pedro Gonz??lez, Pedro Ladr??n de Guevara, Ernesto L??pez Torres, *
 *         Eulogio Serradilla, Ana Marin, Friederike Bock                 *
 * Version 2                                                              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Analysis task for pi0->e+e-gamma (Dalitz decay)
// Analysis task for chic->JPsi+gamma

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
#include "AliAnalysisTaskEtaToPiPlPiMiGamma.h"


ClassImp( AliAnalysisTaskEtaToPiPlPiMiGamma )

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskEtaToPiPlPiMiGamma::AliAnalysisTaskEtaToPiPlPiMiGamma():
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
	fGoodGammas(NULL),
	fGoodVirtualParticles(NULL),
	fGammaCutArray(NULL),
	fPionCutArray(NULL),
	fMesonCutArray(NULL),
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
	fHistoNEvents(NULL),
	fHistoNGoodESDTracks(NULL),
	fProfileEtaShift(NULL),
	fRandom(0),
	fnCuts(0),
	fiCut(0),
	fNumberOfESDTracks(0),
	fMoveParticleAccordingToVertex(kFALSE),
	fIsHeavyIon(kFALSE),
	fDoMesonAnalysis(kTRUE),
	fDoMesonQA(kFALSE),
	fIsFromMBHeader(kTRUE),
	fIsMC(kFALSE)
{

}

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskEtaToPiPlPiMiGamma::AliAnalysisTaskEtaToPiPlPiMiGamma( const char* name ):
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
	fGoodGammas(NULL),
	fGoodVirtualParticles(NULL),
	fGammaCutArray(NULL),
	fPionCutArray(NULL),
	fMesonCutArray(NULL),
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
	fHistoNEvents(NULL),
	fHistoNGoodESDTracks(NULL),
	fProfileEtaShift(NULL),
	fRandom(0),
	fnCuts(0),
	fiCut(0),
	fNumberOfESDTracks(0),
	fMoveParticleAccordingToVertex(kFALSE),
	fIsHeavyIon(kFALSE),
	fDoMesonAnalysis(kTRUE),
	fDoMesonQA(kFALSE),
	fIsFromMBHeader(kTRUE),
	fIsMC(kFALSE)
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
	Int_t nBins[nDim] = {800,250,7,4};
	Double_t xMin[nDim] = {0,0, 0,0};
	Double_t xMax[nDim] = {0.8,25,7,4};
	
	fTHnSparseMotherInvMassPtZM = new THnSparseF*[fnCuts];
	fTHnSparseMotherBackInvMassPtZM = new THnSparseF*[fnCuts];

	fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		
		TString cutstringPion     =   ((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutNumber();
		TString cutstringMeson        =   ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
		TString cutstringGamma        =   ((AliConversionCuts*)fGammaCutArray->At(iCut))->GetCutNumber();
		
		Int_t collisionSystem = atoi((TString)(((AliConversionCuts*)fGammaCutArray->At(iCut))->GetCutNumber())(0,1));
		Int_t centMin = atoi((TString)(((AliConversionCuts*)fGammaCutArray->At(iCut))->GetCutNumber())(1,1));
		Int_t centMax = atoi((TString)(((AliConversionCuts*)fGammaCutArray->At(iCut))->GetCutNumber())(2,1));
		
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
		fBackList[iCut]->SetName(Form("%s_%s_%s Back histograms",cutstringGamma.Data(),cutstringPion.Data(),cutstringMeson.Data()));
		fBackList[iCut]->SetOwner(kTRUE);
		fCutFolder[iCut]->Add(fBackList[iCut]);

		fTHnSparseMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_m","Back_Back_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
		fBackList[iCut]->Add(fTHnSparseMotherBackInvMassPtZM[iCut]);

		fMotherList[iCut] = new TList();
		fMotherList[iCut]->SetName(Form("%s_%s_%s Mother histograms",cutstringGamma.Data(),cutstringPion.Data(),cutstringMeson.Data()));
		fMotherList[iCut]->SetOwner(kTRUE);
		fCutFolder[iCut]->Add(fMotherList[iCut]);

		fTHnSparseMotherInvMassPtZM[iCut] = new THnSparseF("Back_Mother_InvMass_Pt_z_m","Back_Mother_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
		fMotherList[iCut]->Add(fTHnSparseMotherInvMassPtZM[iCut]);

		
		fBGHandler[iCut] = new AliGammaConversionAODBGHandler(	collisionSystem,centMin,centMax,
																((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfBGEvents(),
																((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity());
		
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
	//fGoodGammas->SetOwner(kTRUE);

	fGoodVirtualParticles = new TList();
	//fGoodVirtualParticles->SetOwner(kTRUE);

	fCutFolder				= new TList*[fnCuts];
	fESDList				= new TList*[fnCuts];
	fBackList				= new TList*[fnCuts];
	fMotherList 			= new TList*[fnCuts];
	fHistoNEvents			= new TH1I*[fnCuts];
	fHistoNGoodESDTracks	= new TH1I*[fnCuts];
	fProfileEtaShift		= new TProfile*[fnCuts];
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
	}
	
	fHistoMotherInvMassPt		= new TH2F*[fnCuts];
	fHistoMotherBackInvMassPt	= new TH2F*[fnCuts];

	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		TString cutstringPion =((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutNumber();
		TString cutstringMeson= ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
		TString cutstringGamma = ((AliConversionCuts*)fGammaCutArray->At(iCut))->GetCutNumber();

		fCutFolder[iCut] = new TList();
		fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s",cutstringGamma.Data(),cutstringPion.Data(),cutstringMeson.Data()));
		fCutFolder[iCut]->SetOwner(kTRUE);
		fOutputContainer->Add(fCutFolder[iCut]);

		fESDList[iCut] = new TList();
		fESDList[iCut]->SetName(Form("%s_%s_%s ESD histograms",cutstringGamma.Data(),cutstringPion.Data(),cutstringMeson.Data()));
		fESDList[iCut]->SetOwner(kTRUE);

		fHistoNEvents[iCut] = new TH1I("NEvents","NEvents",9,-0.5,8.5);
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Missing MC");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
		fESDList[iCut]->Add(fHistoNEvents[iCut]);

		if(fIsHeavyIon) fHistoNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",3000,0,3000);
			else fHistoNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",200,0,200);
		fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);

		fProfileEtaShift[iCut] = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
		fESDList[iCut]->Add(fProfileEtaShift[iCut]);
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
		}

		fHistoMotherInvMassPt[iCut] = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",900,0.3,1.2,250,0,25);
		fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);
		fHistoMotherBackInvMassPt[iCut] = new TH2F("ESD_Background_InvMass_Pt","ESD_Background_InvMass_Pt",900,0.3,1.2,250,0,25);
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
		fHistoTrueMotherPiPlPiMiGammaInvMassPt = new TH2F*[fnCuts];
		fHistoTrueMotherDalitzInvMassPt = new TH2F*[fnCuts];
		fHistoTrueMotherGammaGammaInvMassPt = new TH2F*[fnCuts];

		for(Int_t iCut = 0; iCut<fnCuts;iCut++){
			TString cutstringPion =((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->GetCutNumber();
			TString cutstringMeson= ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
			TString cutstringGamma = ((AliConversionCuts*)fGammaCutArray->At(iCut))->GetCutNumber();

			fMCList[iCut] = new TList();
			fMCList[iCut]->SetName(Form("%s_%s_%s MC histograms",cutstringGamma.Data(),cutstringPion.Data(),cutstringMeson.Data()));
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
			
			fHistoMCEtaPiPlPiMiGammaInAccPt[iCut] = new TH1F("MC_EtaDalitzInAcc_Pt","MC_EtaDalitzInAcc_Pt",250,0,25);
			fHistoMCEtaPiPlPiMiGammaInAccPt[iCut]->Sumw2();
			fMCList[iCut]->Add(fHistoMCEtaPiPlPiMiGammaInAccPt[iCut]);

			fTrueList[iCut] = new TList();
			fTrueList[iCut]->SetName(Form("%s_%s_%s True histograms",cutstringGamma.Data(),cutstringPion.Data(),cutstringMeson.Data()));
			fTrueList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fTrueList[iCut]);

			fHistoTrueConvGammaPt[iCut] = new TH1F("ESD_TrueConvGamma_Pt","ESD_TrueConvGamma_Pt",250,0,25);
			fTrueList[iCut]->Add(fHistoTrueConvGammaPt[iCut]);
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

			fHistoTrueMotherPiPlPiMiGammaInvMassPt[iCut] = new TH2F("ESD_TrueMotherPiPlPiMiGamma_InvMass_Pt","ESD_TrueMotherPiPlPiMiGamma_InvMass_Pt",900,0.3,1.4,250,0,25);
			fHistoTrueMotherPiPlPiMiGammaInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(fHistoTrueMotherPiPlPiMiGammaInvMassPt[iCut]);
		
			fHistoTrueMotherGammaGammaInvMassPt[iCut] = new TH2F("ESD_TrueMotherGG_InvMass_Pt","ESD_TrueMotherGG_InvMass_Pt",900,0.3,1.4,250,0,25);
			fHistoTrueMotherGammaGammaInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(fHistoTrueMotherGammaGammaInvMassPt[iCut]);
			
			fHistoTrueMotherDalitzInvMassPt[iCut] = new TH2F("ESD_TrueMotherDalitz_InvMass_Pt","ESD_TrueMotherDalitz_InvMass_Pt",900,0.3,1.4,250,0,25);
			fHistoTrueMotherDalitzInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(fHistoTrueMotherDalitzInvMassPt[iCut]);
			
		}
	}

	

	InitBack(); // Init Background Handler

	fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
	if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader
		
	if(fV0Reader)
		if((AliConversionCuts*)fV0Reader->GetConversionCuts())
			if(((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
				fOutputContainer->Add(((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
		
		
		
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
			if( ((AliConversionCuts*)fGammaCutArray->At(iCut))->GetCutHistograms() ) {
				fCutFolder[iCut]->Add( ((AliConversionCuts*)fGammaCutArray->At(iCut))->GetCutHistograms()  );
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

	fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
	if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader

	Int_t eventQuality = ((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetEventQuality();

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
		Int_t eventNotAccepted =
			((AliConversionCuts*)fGammaCutArray->At(iCut))
			->IsEventAcceptedByConversionCut(fV0Reader->GetConversionCuts(),fInputEvent,fMCEvent,fIsHeavyIon);
		
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
			if(((AliConversionCuts*)fGammaCutArray->At(iCut))->GetSignalRejection() != 0){
				((AliConversionCuts*)fGammaCutArray->At(iCut))->GetNotRejectedParticles(((AliConversionCuts*)fGammaCutArray->At(iCut))->GetSignalRejection(), ((AliConversionCuts*)fGammaCutArray->At(iCut))->GetAcceptedHeader(),
													fMCEvent);
			} 
			ProcessMCParticles();
		}

// 		cout << "new event" << endl;
		ProcessPhotonCandidates(); // Process this cuts gammas
		ProcessPionCandidates(); // Process this cuts gammas
		
	
		CalculateMesonCandidates();
		CalculateBackground();
		UpdateEventByEventData();
				
		fGoodGammas->Clear(); // delete this cuts good gammas
		fGoodVirtualParticles->Clear(); // delete this cuts good gammas
	}

	fSelectorNegPionIndex.clear();
	fSelectorPosPionIndex.clear();

	PostData( 1, fOutputContainer );
}

Bool_t AliAnalysisTaskEtaToPiPlPiMiGamma::Notify(){
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		if( !((AliConversionCuts*)fGammaCutArray->At(iCut))->GetDoEtaShift() ){
			fProfileEtaShift[iCut]->Fill(0.,0.);
			continue; // No Eta Shift requested, continue
		}
		if( ((AliConversionCuts*)fGammaCutArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
			((AliConversionCuts*)fGammaCutArray->At(iCut))->GetCorrectEtaShiftFromPeriod(fV0Reader->GetPeriodName());
			((AliConversionCuts*)fGammaCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once   
			((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->SetEtaShift( ((AliConversionCuts*)fGammaCutArray->At(iCut))->GetEtaShift() );
			fProfileEtaShift[iCut]->Fill(0.,(((AliConversionCuts*)fGammaCutArray->At(iCut))->GetEtaShift()));
			continue;
		} else {
			printf(" Eta t PiPlusPiMinus Gamma Task %s :: Eta Shift Manually Set to %f \n\n",
			(((AliConversionCuts*)fGammaCutArray->At(iCut))->GetCutNumber()).Data(),((AliConversionCuts*)fGammaCutArray->At(iCut))->GetEtaShift());
			((AliConversionCuts*)fGammaCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once   
			((AliPrimaryPionCuts*)fPionCutArray->At(iCut))->SetEtaShift( ((AliConversionCuts*)fGammaCutArray->At(iCut))->GetEtaShift() );
			fProfileEtaShift[iCut]->Fill(0.,(((AliConversionCuts*)fGammaCutArray->At(iCut))->GetEtaShift()));
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
		
		if( fMCEvent && ((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetSignalRejection() != 0 ){		
			Int_t isPosFromMBHeader
				= ((AliConversionCuts*)fGammaCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack, fInputEvent);
			if(isPosFromMBHeader == 0 && ((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
			Int_t isNegFromMBHeader
				= ((AliConversionCuts*)fGammaCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack,fInputEvent);
			if(isNegFromMBHeader == 0 && ((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
			if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
		}
		
		if(!((AliConversionCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fESDEvent)) continue;

		if(!((AliConversionCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut() &&
			!((AliConversionCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){ // if no post reader loop is required add to events good gammas
			
			fGoodGammas->Add(PhotonCandidate);
		
			if(fIsFromMBHeader){
				fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
				fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
			}
		
			if(fMCEvent){
				ProcessTruePhotonCandidates(PhotonCandidate);
			}
		} else if(((AliConversionCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
			((AliConversionCuts*)fGammaCutArray->At(fiCut))->FillElectonLabelArray(PhotonCandidate,nV0);
			nV0++;
			GoodGammasStepOne->Add(PhotonCandidate);
		} else if(!((AliConversionCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut() &&
				((AliConversionCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
			GoodGammasStepTwo->Add(PhotonCandidate);
		}
	}
	
	
	if(((AliConversionCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut()){
		for(Int_t i = 0;i<GoodGammasStepOne->GetEntries();i++){
			AliAODConversionPhoton *PhotonCandidate= (AliAODConversionPhoton*) GoodGammasStepOne->At(i);
			if(!PhotonCandidate) continue;
			fIsFromMBHeader = kTRUE;
			if(fMCEvent && ((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetSignalRejection() != 0){
				Int_t isPosFromMBHeader
				= ((AliConversionCuts*)fGammaCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack,fInputEvent);
				Int_t isNegFromMBHeader
				= ((AliConversionCuts*)fGammaCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack,fInputEvent);
				if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
			}
			if(!((AliConversionCuts*)fGammaCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GoodGammasStepOne->GetEntries())) continue;
			if(!((AliConversionCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
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
	if(((AliConversionCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){
		for(Int_t i = 0;i<GoodGammasStepTwo->GetEntries();i++){
			AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GoodGammasStepTwo->At(i);
			if(!PhotonCandidate) continue;
			
			if(fMCEvent && ((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetSignalRejection() != 0){
				Int_t isPosFromMBHeader
				= ((AliConversionCuts*)fGammaCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack,fInputEvent);
				Int_t isNegFromMBHeader
				= ((AliConversionCuts*)fGammaCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack,fInputEvent);
				if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
			}
			
			if(!((AliConversionCuts*)fGammaCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GoodGammasStepTwo,i)) continue;
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
	AliStack *MCStack = fMCEvent->Stack();
	TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(MCStack);
	TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(MCStack);

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
	
	if( labelGamma < MCStack->GetNprimary() ){
		if( fIsFromMBHeader ){
		fHistoTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
		}
	}
	
	if( IsEtaPiPlPiMiGammaDaughter(labelGamma) == kTRUE ) {
		if( labelGamma < MCStack->GetNprimary() ) {
			if( fIsFromMBHeader ){
				fHistoTrueConvGammaFromEtaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
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
	
	for(UInt_t i = 0; i < fSelectorNegPionIndex.size(); i++){
		AliESDtrack* negPionCandidate = fESDEvent->GetTrack(fSelectorNegPionIndex[i]);
		if(! ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelected(negPionCandidate) ) continue;
		lGoodNegPionIndexPrev.push_back(   fSelectorNegPionIndex[i] );
		fHistoNegPionPt[fiCut]->Fill(negPionCandidate->Pt());
		fHistoNegPionPhi[fiCut]->Fill(negPionCandidate->Phi());
		if( fMCEvent ) {
			Int_t labelNegPion = TMath::Abs( negPionCandidate->GetLabel() );
			if( labelNegPion < fMCStack->GetNtrack() ){
				TParticle* negPion = fMCStack->Particle(labelNegPion);
				if( negPion->GetPdgCode() ==  -211 ){
					if( labelNegPion < fMCStack->GetNprimary() ){
						fHistoTrueNegPionPt[fiCut]->Fill(negPionCandidate->Pt());    //primary negPion
					}		
					if( IsEtaPiPlPiMiGammaDaughter(labelNegPion) == kTRUE ) {
						if( labelNegPion < fMCStack->GetNprimary() ) {
							fHistoTrueNegPionFromEtaPt[fiCut]->Fill(negPionCandidate->Pt());
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
			Int_t labelPosPion = TMath::Abs( posPionCandidate->GetLabel() );
			if( labelPosPion < fMCStack->GetNtrack() ) {
				TParticle* posPion = fMCStack->Particle(labelPosPion);
				if( posPion->GetPdgCode() ==  211 ){
					if( labelPosPion < fMCStack->GetNprimary() ){
						fHistoTruePosPionPt[fiCut]->Fill(posPionCandidate->Pt());
					} 
					if( IsEtaPiPlPiMiGammaDaughter(labelPosPion) == kTRUE ) {
						if( labelPosPion < fMCStack->GetNprimary() ){
							fHistoTruePosPionFromEtaPt[fiCut]->Fill(posPionCandidate->Pt());
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
			
			
			if( fMCEvent ) {
				Int_t labeln=TMath::Abs(negPionCandidate->GetLabel());
				Int_t labelp=TMath::Abs(posPionCandidate->GetLabel());
				TParticle *fNegativeMCParticle = fMCStack->Particle(labeln);
				TParticle *fPositiveMCParticle = fMCStack->Particle(labelp);
				
				if( fPositiveMCParticle && fNegativeMCParticle) {
					virtualPhoton->SetMCLabelPositive(labelp);
					virtualPhoton->SetMCLabelNegative(labeln);
				}
			}
			
			AliAODConversionPhoton *vphoton = new AliAODConversionPhoton(virtualPhoton); //To Apply PsiPairCut
			fGoodVirtualParticles->Add(  vphoton );
			delete virtualPhoton;
			virtualPhoton=NULL;
					
		}
	}
}

//_____________________________________________________________________________
void AliAnalysisTaskEtaToPiPlPiMiGamma::ProcessMCParticles(){

	// Loop over all primary MC particle

	for(Int_t i = 0; i < fMCStack->GetNprimary(); i++) {

		TParticle* particle = (TParticle *)fMCStack->Particle(i);
		if (!particle) continue;

		Int_t isMCFromMBHeader = -1;
		if(((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetSignalRejection() != 0){
			isMCFromMBHeader
				= ((AliConversionCuts*)fGammaCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent);
			if(isMCFromMBHeader == 0 && ((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
		}

		if(((AliConversionCuts*)fGammaCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack,fInputEvent)){
			
			if(((AliConversionCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCStack,kFALSE)){
				fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt()); // All MC Gamma
				if(particle->GetMother(0) >-1){
					if (fMCStack->Particle(particle->GetMother(0))->GetPdgCode() ==221 && fMCStack->Particle(particle->GetMother(0))->GetNDaughters()==3 ) fHistoMCGammaFromEtaPt[fiCut]->Fill(particle->Pt()); // All pos from eta
				}	
			}
			
			if(((AliConversionCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCStack,kTRUE)){
				fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt());
			} // Converted MC Gamma
			
			if(((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(i,fMCStack)){
				if( particle->GetPdgCode() == 211){
					fHistoMCAllPosPionsPt[fiCut]->Fill(particle->Pt()); // All pos pions
					if(particle->GetMother(0) >-1){
						if (fMCStack->Particle(particle->GetMother(0))->GetPdgCode() ==221) fHistoMCPosPionsFromEtaPt[fiCut]->Fill(particle->Pt()); // All pos from eta
					}	
				}	
				if( particle->GetPdgCode() == -211){
					fHistoMCAllNegPionsPt[fiCut]->Fill(particle->Pt()); // All neg pions
					if(particle->GetMother(0) >-1){
						if (fMCStack->Particle(particle->GetMother(0))->GetPdgCode() ==221) fHistoMCNegPionsFromEtaPt[fiCut]->Fill(particle->Pt()); // All pos from eta
					}	
				}
			}
			
			// \eta -> \gamma \gamma 
			
			if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMC( particle,fMCStack,((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetEtaShift() ) ){
				Float_t weighted= 1;
				if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) { 
					if(((AliConversionCuts*)fGammaCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent)){
						if (particle->Pt()>0.005){
							weighted= ((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, fMCStack,fInputEvent);
						}
					}
				}
				if(particle->GetPdgCode() == 221)fHistoMCEtaGGPt[fiCut]->Fill( particle->Pt() , weighted); // All MC Eta GG decay
			}
			
			// \eta -> e+ e- \gamma 
			Int_t labelgamma 	  = -1;
			Int_t labelelectron = -1;
			Int_t labelpositron = -1;

			if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMCDalitz(particle,fMCStack,labelelectron,labelpositron,labelgamma,((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetEtaShift())){
				Float_t weighted= 1;
				if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) { 
					if(((AliConversionCuts*)fGammaCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack,fInputEvent)){
						if (particle->Pt()>0.005){
							weighted= ((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, fMCStack,fInputEvent);
						}
					}
				}
				if(particle->GetPdgCode() == 221)fHistoMCEtaDalitzPt[fiCut]->Fill(particle->Pt(), weighted); // All MC Eta
			}		
			
			
			// \eta -> pi+ pi- \gamma 
			Int_t labelGamma3Body 	  = -1;
			Int_t labelNegPion = -1;
			Int_t labelPosPion = -1;

			if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMCEtaPiPlPiMiGamma(particle,fMCStack,labelNegPion,labelPosPion,labelGamma3Body,((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetEtaShift())){	
				Float_t weighted= 1;
				if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) { 
					if(((AliConversionCuts*)fGammaCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack,fInputEvent)){
						if (particle->Pt()>0.005){
							weighted= ((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, fMCStack,fInputEvent);
						}
					}
				}
				if(particle->GetPdgCode() == 221)fHistoMCEtaPiPlPiMiGammaPt[fiCut]->Fill(particle->Pt(), weighted); // All MC Eta
		
				TParticle *gamma    = fMCStack->Particle(labelGamma3Body);
				if(((AliConversionCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(gamma,fMCStack,kFALSE) &&
				((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelNegPion,fMCStack) &&
				((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->PionIsSelectedMC(labelPosPion,fMCStack) ) {
					if(particle->GetPdgCode() == 221)fHistoMCEtaPiPlPiMiGammaInAccPt[fiCut]->Fill(particle->Pt(), weighted ); // MC EtaDalitz with gamma and e+e- in acc
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

				AliAODConversionPhoton *Vgamma=dynamic_cast<AliAODConversionPhoton*>(fGoodVirtualParticles->At(virtualParticleIndex));
				if (Vgamma==NULL) continue;
				//Check for same Electron ID
				if(gamma->GetTrackLabelPositive() == Vgamma->GetTrackLabelPositive() ||
				gamma->GetTrackLabelNegative() == Vgamma->GetTrackLabelNegative() ||
				gamma->GetTrackLabelNegative() == Vgamma->GetTrackLabelPositive() ||
				gamma->GetTrackLabelPositive() == Vgamma->GetTrackLabelNegative() ) continue;

				AliAODConversionMother *etacand = new AliAODConversionMother(gamma,Vgamma);
				etacand->SetLabels(GammaIndex,virtualParticleIndex);
						

				if( ( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(etacand,kTRUE,((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetEtaShift())) ){
			
					//cout<< "Meson Accepted "<<endl;
					
					Int_t zbin= fBGHandler[fiCut]->GetZBinIndex(fESDEvent->GetPrimaryVertex()->GetZ());
					Int_t mbin = 0;
					if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
						mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
					} else {
						mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGoodGammas->GetEntries());
					}
					
					AliESDtrack *posPionVgamma = 0;
					AliESDtrack *negPionVgamma = 0;
					
					Double_t clsToFPos = -1.0;
					Double_t clsToFNeg = -1.0;
					
					Float_t dcaToVertexXYPos = -1.0;
					Float_t dcaToVertexZPos  = -1.0;
					Float_t dcaToVertexXYNeg = -1.0;
					Float_t dcaToVertexZNeg  = -1.0;
					
					
					if ( fDoMesonQA ) {
					
						posPionVgamma = fESDEvent->GetTrack( Vgamma->GetTrackLabelPositive() );
						negPionVgamma = fESDEvent->GetTrack( Vgamma->GetTrackLabelNegative() );
						clsToFPos = ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetNFindableClustersTPC(posPionVgamma);
						clsToFNeg = ((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetNFindableClustersTPC(negPionVgamma);
						
						Float_t bPos[2];
						Float_t bCovPos[3];
						posPionVgamma->GetImpactParameters(bPos,bCovPos);
						if (bCovPos[0]<=0 || bCovPos[2]<=0) {
							AliDebug(1, "Estimated b resolution lower or equal zero!");
							bCovPos[0]=0; bCovPos[2]=0;
						}
						
						Float_t bNeg[2];
						Float_t bCovNeg[3];
						posPionVgamma->GetImpactParameters(bNeg,bCovNeg);
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
					
							fHistoNegPionEta[fiCut]->Fill( negPionVgamma->Eta() );
							fHistoPosPionEta[fiCut]->Fill( posPionVgamma->Eta() );
									
							fHistoNegPionClsTPC[fiCut]->Fill(clsToFNeg,negPionVgamma->Pt());
							fHistoPosPionClsTPC[fiCut]->Fill(clsToFPos,posPionVgamma->Pt());
							
							fHistoPionDCAxy[fiCut]->Fill(  dcaToVertexXYNeg, negPionVgamma->Pt() );
							fHistoPionDCAz[fiCut]->Fill(   dcaToVertexZNeg,  negPionVgamma->Pt() );
							fHistoPionDCAxy[fiCut]->Fill(  dcaToVertexXYPos, posPionVgamma->Pt() );
							fHistoPionDCAz[fiCut]->Fill(   dcaToVertexZPos,  posPionVgamma->Pt() );
							
							fHistoPionTPCdEdxNSigma[fiCut]->Fill( posPionVgamma->P(),((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetPIDResponse()->NumberOfSigmasTPC(posPionVgamma, AliPID::kPion) );
							fHistoPionTPCdEdxNSigma[fiCut]->Fill( negPionVgamma->P(),((AliPrimaryPionCuts*)fPionCutArray->At(fiCut))->GetPIDResponse()->NumberOfSigmasTPC(negPionVgamma, AliPID::kPion) );
							
							fHistoPionTPCdEdx[fiCut]->Fill( posPionVgamma->P(), TMath::Abs(posPionVgamma->GetTPCsignal()));
							fHistoPionTPCdEdx[fiCut]->Fill( negPionVgamma->P(), TMath::Abs(negPionVgamma->GetTPCsignal()));
							
							lGoodVirtualParticle[virtualParticleIndex] = kTRUE;
					
						}
					}
		
				
					if(fMCEvent){
						ProcessTrueMesonCandidates(etacand,gamma,Vgamma);
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

				for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
					AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));

					if(fMoveParticleAccordingToVertex == kTRUE && method == 1 ){
						MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
					}

					AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
								

					if( ( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE, ((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetEtaShift()))){
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

					for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){

						AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));

						if(fMoveParticleAccordingToVertex == kTRUE && method ==1){
							MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
						}

						AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
								
						if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE,((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetEtaShift()))){
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

	AliStack *MCStack = fMCEvent->Stack();

	if(	TrueGammaCandidate->GetV0Index()<fESDEvent->GetNumberOfV0s()	){
		
		Bool_t isTrueEta = kFALSE;
		Int_t gammaMCLabel = TrueGammaCandidate->GetMCParticleLabel(MCStack);
		Int_t gammaMotherLabel = -1;
		
		if(gammaMCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
			// Daughters Gamma 0
			TParticle * negativeMC = (TParticle*)TrueGammaCandidate->GetNegativeMCDaughter(MCStack);
			TParticle * positiveMC = (TParticle*)TrueGammaCandidate->GetPositiveMCDaughter(MCStack);
			TParticle * gammaMC = (TParticle*)MCStack->Particle(gammaMCLabel);

			if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
				if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
					if(gammaMC->GetPdgCode() == 22){ // ... with Gamma Mother
						gammaMotherLabel=gammaMC->GetFirstMother();
					}
				}
			}
		}

		Int_t virtualParticleMCLabel = TrueVirtualParticleCandidate->GetMCParticleLabel(MCStack);
		Int_t virtualParticleMotherLabel = -1;

		Bool_t isPiPiDecay = kFALSE;
		Bool_t isDalitz = kFALSE;
		Bool_t isRealGamma = kFALSE;
		
		if(virtualParticleMCLabel != -1){ // if virtualParticleMCLabel==-1 particles don't have same mother 
			TParticle * negativeMC = (TParticle*)TrueVirtualParticleCandidate->GetNegativeMCDaughter(MCStack);
			TParticle * positiveMC = (TParticle*)TrueVirtualParticleCandidate->GetPositiveMCDaughter(MCStack);
// 			cout << "neg Part: label - " <<  TrueVirtualParticleCandidate->GetMCLabelNegative() <<" pdg-code - " << negativeMC->GetPdgCode() << endl;
// 			cout << "pos Part: label - " <<  TrueVirtualParticleCandidate->GetMCLabelPositive() <<" pdg-code - " << positiveMC->GetPdgCode() << endl;
			
			TParticle * virtualParticleMotherMC = (TParticle*)MCStack->Particle(virtualParticleMCLabel);
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
			if(((TParticle*)MCStack->Particle(virtualParticleMotherLabel))->GetPdgCode() == 221){
				isTrueEta=kTRUE;
			}
		}

		if( isTrueEta ){ // True Eta
			if ( isPiPiDecay) { //real eta -> Pi+ Pi- Gamma
				Float_t weighted= 1;
				if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) { 
					if(((AliConversionCuts*)fGammaCutArray->At(fiCut))->IsParticleFromBGEvent(gammaMotherLabel, fMCStack,fInputEvent)){
						if (((TParticle*)MCStack->Particle(gammaMotherLabel))->Pt()>0.005){
							weighted= ((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),gammaMotherLabel,fMCStack,fInputEvent);
						}
					}
				}

				fHistoTrueMotherPiPlPiMiGammaInvMassPt[fiCut]->Fill(EtaCandidate->M(),EtaCandidate->Pt(),weighted);
			} else if ( isRealGamma ){
				Float_t weighted= 1;
				if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) {
					if(((AliConversionCuts*)fGammaCutArray->At(fiCut))->IsParticleFromBGEvent(gammaMotherLabel, fMCStack,fInputEvent)){
						if (((TParticle*)MCStack->Particle(gammaMotherLabel))->Pt()>0.005){
							weighted= ((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),gammaMotherLabel,fMCStack,fInputEvent);
						}
					}
				}

				fHistoTrueMotherGammaGammaInvMassPt[fiCut]->Fill(EtaCandidate->M(),EtaCandidate->Pt(),weighted); 
			} else if (isDalitz) {
				Float_t weighted= 1;
				if( ((AliPrimaryPionCuts*) fPionCutArray->At(fiCut))->DoWeights() ) {
					if(((AliConversionCuts*)fGammaCutArray->At(fiCut))->IsParticleFromBGEvent(gammaMotherLabel, fMCStack,fInputEvent)){
						if (((TParticle*)MCStack->Particle(gammaMotherLabel))->Pt()>0.005){
							weighted= ((AliConversionCuts*)fGammaCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),gammaMotherLabel,fMCStack,fInputEvent);
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
	Int_t motherLabel = fMCStack->Particle( label )->GetMother(0);
	if( motherLabel < 0 || motherLabel >= fMCStack->GetNtrack() ) return kFALSE;
	
	TParticle* mother = fMCStack->Particle( motherLabel );
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
		TParticle* temp = (TParticle*)fMCStack->Particle( index );
		
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
