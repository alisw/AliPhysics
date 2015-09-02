/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *									 									  *
 * Author: Friederike Bock				 								  *
 * Version 1.0								 							  *
 *									 									  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its	  *
 * documentation strictly for non-commercial purposes is hereby granted	  *
 * without fee, provided that the above copyright notice appears in all	  *
 * copies and that both the copyright notice and this permission notice	  *
 * appear in the supporting documentation. The authors make no claims	  *
 * about the suitability of this software for any purpose. It is		  *
 * provided "as is" without express or implied warranty.       			  *
 **************************************************************************/

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class used to do analysis on conversion photons + calo photons
//----------------------------------------------------------------
//////////////////////////////////////////////////////////////////
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "AliESDpid.h"
#include "AliAnalysisTaskGammaCaloMerged.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliKFVertex.h"
#include "AliV0ReaderV1.h"
#include "AliGenCocktailEventHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliEventplane.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"
#include <vector>
#include <map>

ClassImp(AliAnalysisTaskGammaCaloMerged)

//________________________________________________________________________
AliAnalysisTaskGammaCaloMerged::AliAnalysisTaskGammaCaloMerged(): AliAnalysisTaskSE(),
	fV0Reader(NULL),
	fInputEvent(NULL),
	fMCEvent(NULL),
	fMCStack(NULL),
	fCutFolder(NULL),
	fESDList(NULL),
	fTrueList(NULL),
	fMCList(NULL),
	fHeaderNameList(NULL),
	fOutputContainer(NULL),
	fNClusterCandidates(0),
	fNClusterMergedCandidates(0),
	fEventCutArray(NULL),
	fEventCuts(NULL),
	fClusterCutArray(NULL),
	fClusterMergedCutArray(NULL),
	fMesonCutArray(NULL),
	fMesonCuts(NULL),
	fHistoMotherInvMassPt(NULL),
	fHistoMotherPi0PtY(NULL),
	fHistoMotherEtaPtY(NULL),
	fHistoMotherPi0PtAlpha(NULL),
	fHistoMotherEtaPtAlpha(NULL),
	fHistoMotherPi0PtOpenAngle(NULL),
	fHistoMotherEtaPtOpenAngle(NULL),
	fHistoClusGammaPt(NULL),
	fHistoClusOverlapHeadersGammaPt(NULL),
	fHistoClusMergedPtvsM02(NULL),
	fHistoMCHeaders(NULL),
	fHistoMCPi0Pt(NULL),
	fHistoMCPi0WOWeightPt(NULL),
	fHistoMCPi0WOEvtWeightPt(NULL),
	fHistoMCEtaPt(NULL),
	fHistoMCEtaWOWeightPt(NULL),
	fHistoMCEtaWOEvtWeightPt(NULL),
	fHistoMCPi0InAccPt(NULL),
	fHistoMCEtaInAccPt(NULL),
	fHistoMCPi0PtY(NULL),
	fHistoMCEtaPtY(NULL),
	fHistoMCPi0PtAlpha(NULL),
	fHistoMCEtaPtAlpha(NULL),
	fHistoMCPi0PtJetPt(NULL),
	fHistoMCEtaPtJetPt(NULL),
	fHistoTrueClusMergedPtvsM02(NULL),	
	fHistoTrueClusPi0PtvsM02(NULL),
	fHistoTrueClusPrimPi0PtvsM02(NULL),
	fHistoTrueClusSecPi0PtvsM02(NULL),
	fHistoTrueClusSecPi0FromK0sPtvsM02(NULL),
	fHistoTrueClusSecPi0FromLambdaPtvsM02(NULL),
	fHistoTrueClusEtaPtvsM02(NULL),
	fHistoTrueClusMergedPartConvPtvsM02(NULL),
	fHistoTrueClusMergedPartConvELeadPtvsM02(NULL),
	fHistoTrueClusPartConvPi0PtvsM02(NULL),
	fHistoTrueClusPartConvPrimPi0PtvsM02(NULL),
	fHistoTrueClusPartConvSecPi0PtvsM02(NULL),
	fHistoTrueClusPartConvSecPi0FromK0sPtvsM02(NULL),
	fHistoTrueClusPartConvSecPi0FromLambdaPtvsM02(NULL),
	fHistoTrueClusPartConvEtaPtvsM02(NULL),
	fHistoTrueClusBGPtvsM02(NULL),
	fHistoTrueClusGammaPtvsM02(NULL),
	fHistoTrueClusMergedInvMassvsPt(NULL),	
	fHistoTrueClusPi0InvMassvsPt(NULL),
	fHistoTrueClusPrimPi0InvMassvsPt(NULL),
	fHistoTrueClusSecPi0InvMassvsPt(NULL),
	fHistoTrueClusSecPi0FromK0sInvMassvsPt(NULL),
	fHistoTrueClusSecPi0FromLambdaInvMassvsPt(NULL),
	fHistoTrueClusEtaInvMassvsPt(NULL),
	fHistoTrueClusMergedPartConvInvMassvsPt(NULL),
	fHistoTrueClusMergedPartConvELeadInvMassvsPt(NULL),
	fHistoTrueClusPartConvPi0InvMassvsPt(NULL),
	fHistoTrueClusPartConvPrimPi0InvMassvsPt(NULL),
	fHistoTrueClusPartConvSecPi0InvMassvsPt(NULL),
	fHistoTrueClusPartConvSecPi0FromK0sInvMassvsPt(NULL),
	fHistoTrueClusPartConvSecPi0FromLambdaInvMassvsPt(NULL),
	fHistoTrueClusPartConvEtaInvMassvsPt(NULL),
	fHistoTrueClusBGInvMassvsPt(NULL),
	fHistoTrueClusGammaInvMassvsPt(NULL),	
	fHistoTrueClusBGPtvsSource(NULL),
	fHistoTruePi0PtY(NULL),
	fHistoTrueEtaPtY(NULL),
	fHistoTruePi0PtAlpha(NULL),
	fHistoTrueEtaPtAlpha(NULL),
	fHistoTruePi0PtOpenAngle(NULL),
	fHistoTrueEtaPtOpenAngle(NULL),
	fHistoDoubleCountTruePi0InvMassPt(NULL),
	fHistoDoubleCountTrueEtaInvMassPt(NULL),
	fVectorDoubleCountTruePi0s(0),
	fVectorDoubleCountTrueEtas(0),
	fHistoNEvents(NULL),
	fHistoNEventsWOWeight(NULL),
	fHistoNGoodESDTracks(NULL),
	fHistoVertexZ(NULL),
	fHistoNClusterCandidates(NULL),
	fHistoNClusterMergedCandidates(NULL),
	fHistoNGoodESDTracksVsNClusterCandidates(NULL),
	fHistoSPDClusterTrackletBackground(NULL),
	fHistoNV0Tracks(NULL),
	fProfileEtaShift(NULL),
	fProfileJetJetXSection(NULL),
	fHistoJetJetNTrials(NULL),
	fRandom(0),
	fnCuts(0),
	fiCut(0),
	fIsHeavyIon(0),
	fDoMesonQA(0),
	fDoClusterQA(0),
	fIsFromMBHeader(kTRUE),
	fIsOverlappingWithOtherHeader(kFALSE),
	fIsMC(0),
	fSetPlotHistsExtQA(kFALSE),
	fWeightJetJetMC(1)
{
  
}

//________________________________________________________________________
AliAnalysisTaskGammaCaloMerged::AliAnalysisTaskGammaCaloMerged(const char *name):
	AliAnalysisTaskSE(name),
	fV0Reader(NULL),
	fInputEvent(NULL),
	fMCEvent(NULL),
	fMCStack(NULL),
	fCutFolder(NULL),
	fESDList(NULL),
	fTrueList(NULL),
	fMCList(NULL),
	fHeaderNameList(NULL),
	fOutputContainer(NULL),
	fNClusterCandidates(0),
	fNClusterMergedCandidates(0),
	fEventCutArray(NULL),
	fEventCuts(NULL),
	fClusterCutArray(NULL),
	fClusterMergedCutArray(NULL),
	fMesonCutArray(NULL),
	fMesonCuts(NULL),
	fHistoMotherInvMassPt(NULL),
	fHistoMotherPi0PtY(NULL),
	fHistoMotherEtaPtY(NULL),
	fHistoMotherPi0PtAlpha(NULL),
	fHistoMotherEtaPtAlpha(NULL),
	fHistoMotherPi0PtOpenAngle(NULL),
	fHistoMotherEtaPtOpenAngle(NULL),
	fHistoClusGammaPt(NULL),
	fHistoClusOverlapHeadersGammaPt(NULL),
	fHistoClusMergedPtvsM02(NULL),
	fHistoMCHeaders(NULL),
	fHistoMCPi0Pt(NULL),
	fHistoMCPi0WOWeightPt(NULL),
	fHistoMCPi0WOEvtWeightPt(NULL),
	fHistoMCEtaPt(NULL),
	fHistoMCEtaWOWeightPt(NULL),
	fHistoMCEtaWOEvtWeightPt(NULL),
	fHistoMCPi0InAccPt(NULL),
	fHistoMCEtaInAccPt(NULL),
	fHistoMCPi0PtY(NULL),
	fHistoMCEtaPtY(NULL),
	fHistoMCPi0PtAlpha(NULL),
	fHistoMCEtaPtAlpha(NULL),
	fHistoMCPi0PtJetPt(NULL),
	fHistoMCEtaPtJetPt(NULL),
	fHistoTrueClusMergedPtvsM02(NULL),	
	fHistoTrueClusPi0PtvsM02(NULL),
	fHistoTrueClusPrimPi0PtvsM02(NULL),
	fHistoTrueClusSecPi0PtvsM02(NULL),
	fHistoTrueClusSecPi0FromK0sPtvsM02(NULL),
	fHistoTrueClusSecPi0FromLambdaPtvsM02(NULL),
	fHistoTrueClusEtaPtvsM02(NULL),
	fHistoTrueClusMergedPartConvPtvsM02(NULL),
	fHistoTrueClusMergedPartConvELeadPtvsM02(NULL),
	fHistoTrueClusPartConvPi0PtvsM02(NULL),
	fHistoTrueClusPartConvPrimPi0PtvsM02(NULL),
	fHistoTrueClusPartConvSecPi0PtvsM02(NULL),
	fHistoTrueClusPartConvSecPi0FromK0sPtvsM02(NULL),
	fHistoTrueClusPartConvSecPi0FromLambdaPtvsM02(NULL),
	fHistoTrueClusPartConvEtaPtvsM02(NULL),
	fHistoTrueClusBGPtvsM02(NULL),
	fHistoTrueClusGammaPtvsM02(NULL),
	fHistoTrueClusMergedInvMassvsPt(NULL),	
	fHistoTrueClusPi0InvMassvsPt(NULL),
	fHistoTrueClusPrimPi0InvMassvsPt(NULL),
	fHistoTrueClusSecPi0InvMassvsPt(NULL),
	fHistoTrueClusSecPi0FromK0sInvMassvsPt(NULL),
	fHistoTrueClusSecPi0FromLambdaInvMassvsPt(NULL),
	fHistoTrueClusEtaInvMassvsPt(NULL),
	fHistoTrueClusMergedPartConvInvMassvsPt(NULL),
	fHistoTrueClusMergedPartConvELeadInvMassvsPt(NULL),
	fHistoTrueClusPartConvPi0InvMassvsPt(NULL),
	fHistoTrueClusPartConvPrimPi0InvMassvsPt(NULL),
	fHistoTrueClusPartConvSecPi0InvMassvsPt(NULL),
	fHistoTrueClusPartConvSecPi0FromK0sInvMassvsPt(NULL),
	fHistoTrueClusPartConvSecPi0FromLambdaInvMassvsPt(NULL),
	fHistoTrueClusPartConvEtaInvMassvsPt(NULL),
	fHistoTrueClusBGInvMassvsPt(NULL),
	fHistoTrueClusGammaInvMassvsPt(NULL),
	fHistoTrueClusBGPtvsSource(NULL),
	fHistoTruePi0PtY(NULL),
	fHistoTrueEtaPtY(NULL),
	fHistoTruePi0PtAlpha(NULL),
	fHistoTrueEtaPtAlpha(NULL),
	fHistoTruePi0PtOpenAngle(NULL),
	fHistoTrueEtaPtOpenAngle(NULL),
	fHistoDoubleCountTruePi0InvMassPt(NULL),
	fHistoDoubleCountTrueEtaInvMassPt(NULL),
	fVectorDoubleCountTruePi0s(0),
	fVectorDoubleCountTrueEtas(0),
	fHistoNEvents(NULL),
	fHistoNEventsWOWeight(NULL),
	fHistoNGoodESDTracks(NULL),
	fHistoVertexZ(NULL),
	fHistoNClusterCandidates(NULL),
	fHistoNClusterMergedCandidates(NULL),
	fHistoNGoodESDTracksVsNClusterCandidates(NULL),
	fHistoSPDClusterTrackletBackground(NULL),
	fHistoNV0Tracks(NULL),
	fProfileEtaShift(NULL),
	fProfileJetJetXSection(NULL),
	fHistoJetJetNTrials(NULL),
	fRandom(0),
	fnCuts(0),
	fiCut(0),
	fIsHeavyIon(0),
	fDoMesonQA(0),
	fDoClusterQA(0),
	fIsFromMBHeader(kTRUE),
	fIsOverlappingWithOtherHeader(kFALSE),
	fIsMC(0),
	fSetPlotHistsExtQA(kFALSE),
	fWeightJetJetMC(1)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaCaloMerged::~AliAnalysisTaskGammaCaloMerged()
{

}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::UserCreateOutputObjects(){
  
	if (fIsMC == 2){
		fDoClusterQA = 0;
	}

	// Create histograms
	if(fOutputContainer != NULL){
		delete fOutputContainer;
		fOutputContainer = NULL;
	}
	if(fOutputContainer == NULL){
		fOutputContainer = new TList();
		fOutputContainer->SetOwner(kTRUE);
	}
  	
	fCutFolder = new TList*[fnCuts];
	fESDList = new TList*[fnCuts];

	fHistoNEvents = new TH1F*[fnCuts];
	if(fIsMC == 2){
		fHistoNEventsWOWeight = new TH1F*[fnCuts];
		fProfileJetJetXSection = new TProfile*[fnCuts];
		fHistoJetJetNTrials = new TH1F*[fnCuts];
	}
	
	fHistoNGoodESDTracks = new TH1F*[fnCuts];
	fHistoVertexZ = new TH1F*[fnCuts];
	fHistoNClusterCandidates = new TH1F*[fnCuts];
	fHistoNClusterMergedCandidates = new TH1F*[fnCuts];
	fHistoNGoodESDTracksVsNClusterCandidates = new TH2F*[fnCuts];
	fHistoSPDClusterTrackletBackground = new TH2F*[fnCuts];
	fHistoNV0Tracks = new TH1F*[fnCuts];
	fProfileEtaShift = new TProfile*[fnCuts];
  	
	fHistoMotherInvMassPt = new TH2F*[fnCuts];
	if (fDoMesonQA > 0 && fDoMesonQA < 3){
		fHistoMotherPi0PtY =  new TH2F*[fnCuts];
		fHistoMotherEtaPtY =  new TH2F*[fnCuts];
		fHistoMotherPi0PtAlpha =  new TH2F*[fnCuts];
		fHistoMotherEtaPtAlpha =  new TH2F*[fnCuts];
		fHistoMotherPi0PtOpenAngle =  new TH2F*[fnCuts];
		fHistoMotherEtaPtOpenAngle =  new TH2F*[fnCuts];
	}
		
	fHistoClusGammaPt = new TH1F*[fnCuts];
	fHistoClusOverlapHeadersGammaPt = new TH1F*[fnCuts];
	fHistoClusMergedPtvsM02  =  new TH2F*[fnCuts];
	
	
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		TString cutstringEvent 			= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
		TString cutstringCalo 			= ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
		TString cutstringCaloMerged 	= ((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))->GetCutNumber();
		TString cutstringMeson 			= ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
    
		Int_t nLMCut 	= ((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))->GetMinNLMCut();
		
		fCutFolder[iCut] = new TList();
		fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s_%s",cutstringEvent.Data(), cutstringCalo.Data(), cutstringCaloMerged.Data(), cutstringMeson.Data()));
		fCutFolder[iCut]->SetOwner(kTRUE);
		fOutputContainer->Add(fCutFolder[iCut]);
		fESDList[iCut] = new TList();
		fESDList[iCut]->SetName(Form("%s_%s_%s_%s ESD histograms",cutstringEvent.Data(), cutstringCalo.Data(), cutstringCaloMerged.Data(), cutstringMeson.Data()));
		fESDList[iCut]->SetOwner(kTRUE);
		fCutFolder[iCut]->Add(fESDList[iCut]);
    
		fHistoNEvents[iCut] = new TH1F("NEvents","NEvents",11,-0.5,10.5);
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Missing MC");
		if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){ 
			TString TriggerNames = "Not Trigger: ";
			TriggerNames = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
			fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
		} else {
			fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
		}
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problems");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
		fESDList[iCut]->Add(fHistoNEvents[iCut]);
		
		if (fIsMC == 2){
			fHistoNEventsWOWeight[iCut] = new TH1F("NEventsWOWeight","NEventsWOWeight",11,-0.5,10.5);
			fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
			fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
			fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(3,"Missing MC");
			if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
				TString TriggerNames = "Not Trigger: ";
				TriggerNames = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
				fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
			} else {
				fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
			}
			fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
			fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
			fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
			fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
			fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
			fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problem");
			fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
			fESDList[iCut]->Add(fHistoNEventsWOWeight[iCut]);

			fProfileJetJetXSection[iCut] = new TProfile("XSection","XSection",1,-0.5,0.5);
			fESDList[iCut]->Add(fProfileJetJetXSection[iCut]);
			fHistoJetJetNTrials[iCut] = new TH1F("NTrials","#sum{NTrials}",1,0,1);
			fHistoJetJetNTrials[iCut]->GetXaxis()->SetBinLabel(1,"#sum{NTrials}");
			fESDList[iCut]->Add(fHistoJetJetNTrials[iCut]);
		}	

		if(fIsHeavyIon == 1) fHistoNGoodESDTracks[iCut] = new TH1F("GoodESDTracks","GoodESDTracks",4000,0,4000);
		else if(fIsHeavyIon == 2) fHistoNGoodESDTracks[iCut] = new TH1F("GoodESDTracks","GoodESDTracks",400,0,400);
		else fHistoNGoodESDTracks[iCut] = new TH1F("GoodESDTracks","GoodESDTracks",200,0,200);
		fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);

		fHistoVertexZ[iCut] = new TH1F("VertexZ","VertexZ",1000,-50,50);
		fESDList[iCut]->Add(fHistoVertexZ[iCut]);

		if(fIsHeavyIon == 1) fHistoNClusterCandidates[iCut] = new TH1F("GammaCandidates","GammaCandidates",600,0,600);
		else if(fIsHeavyIon == 2) fHistoNClusterCandidates[iCut] = new TH1F("GammaCandidates","GammaCandidates",400,0,400);
		else fHistoNClusterCandidates[iCut] = new TH1F("GammaCandidates","GammaCandidates",100,0,100);
		fESDList[iCut]->Add(fHistoNClusterCandidates[iCut]);

		if(fIsHeavyIon == 1) fHistoNClusterMergedCandidates[iCut] = new TH1F("MergedCandidates","MergedCandidates",600,0,600);
		else if(fIsHeavyIon == 2) fHistoNClusterMergedCandidates[iCut] = new TH1F("MergedCandidates","MergedCandidates",400,0,400);
		else fHistoNClusterMergedCandidates[iCut] = new TH1F("MergedCandidates","MergedCandidates",100,0,100);
		fESDList[iCut]->Add(fHistoNClusterMergedCandidates[iCut]);

		
		if(fIsHeavyIon == 1) fHistoNGoodESDTracksVsNClusterCandidates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",4000,0,4000,600,0,600);
		else if(fIsHeavyIon == 2) fHistoNGoodESDTracksVsNClusterCandidates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",400,0,400,400,0,400);
		else fHistoNGoodESDTracksVsNClusterCandidates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",200,0,200,100,0,100);
		fESDList[iCut]->Add(fHistoNGoodESDTracksVsNClusterCandidates[iCut]);

		fHistoSPDClusterTrackletBackground[iCut] = new TH2F("SPD tracklets vs SPD clusters","SPD tracklets vs SPD clusters",100,0,200,250,0,1000);
		fESDList[iCut]->Add(fHistoSPDClusterTrackletBackground[iCut]);
		
		if(fIsHeavyIon == 1) fHistoNV0Tracks[iCut] = new TH1F("V0 Multiplicity","V0 Multiplicity",30000,0,30000);
		else if(fIsHeavyIon == 2) fHistoNV0Tracks[iCut] = new TH1F("V0 Multiplicity","V0 Multiplicity",2500,0,2500);
		else fHistoNV0Tracks[iCut] = new TH1F("V0 Multiplicity","V0 Multiplicity",1500,0,1500);
		fESDList[iCut]->Add(fHistoNV0Tracks[iCut]);
		fProfileEtaShift[iCut] = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
		fESDList[iCut]->Add(fProfileEtaShift[iCut]);

		if (fIsMC == 2){
			fHistoNEvents[iCut]->Sumw2();
			fHistoNGoodESDTracks[iCut]->Sumw2();
			fHistoVertexZ[iCut]->Sumw2();
			fHistoNClusterCandidates[iCut]->Sumw2();
			fHistoNGoodESDTracksVsNClusterCandidates[iCut]->Sumw2();
			fHistoSPDClusterTrackletBackground[iCut]->Sumw2();
			fHistoNV0Tracks[iCut]->Sumw2();
			fProfileEtaShift[iCut]->Sumw2();
		}

   		fHistoClusGammaPt[iCut] = new TH1F("ClusGamma_Pt","ClusGamma_Pt",500, 0, 50);
		fESDList[iCut]->Add(fHistoClusGammaPt[iCut]);
		fHistoClusOverlapHeadersGammaPt[iCut] = new TH1F("ClusGammaOverlapHeaders_Pt","ClusGammaOverlapHeaders_Pt",500, 0, 50);
		fESDList[iCut]->Add(fHistoClusOverlapHeadersGammaPt[iCut]);
		fHistoClusMergedPtvsM02[iCut]  = new TH2F("ClusMerged_Pt_M02","ClusMerged_Pt_M02",500,0,50,500,0,5);
		fESDList[iCut]->Add(fHistoClusMergedPtvsM02[iCut]);

		if (fIsMC == 2){
			fHistoClusGammaPt[iCut]->Sumw2();
			fHistoClusOverlapHeadersGammaPt[iCut]->Sumw2();
			fHistoClusMergedPtvsM02[iCut]->Sumw2();
		}
		
	
		fHistoMotherInvMassPt[iCut] = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",800,0,0.8,500, 0, 50);
		fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);

		if (fIsMC == 2){
			fHistoMotherInvMassPt[iCut]->Sumw2();
		}
		
		if (fDoMesonQA > 0 && fDoMesonQA < 3 ){
			fHistoMotherPi0PtY[iCut] = new TH2F("ESD_MotherPi0_Pt_Y","ESD_MotherPi0_Pt_Y",350,0.03,35.,150,-1.5,1.5);
			SetLogBinningXTH2(fHistoMotherPi0PtY[iCut]);
			fESDList[iCut]->Add(fHistoMotherPi0PtY[iCut]);
			fHistoMotherEtaPtY[iCut] = new TH2F("ESD_MotherEta_Pt_Y","ESD_MotherEta_Pt_Y",350,0.03,35.,150,-1.5,1.5);
			SetLogBinningXTH2(fHistoMotherEtaPtY[iCut]);
			fESDList[iCut]->Add(fHistoMotherEtaPtY[iCut]);
			fHistoMotherPi0PtAlpha[iCut] = new TH2F("ESD_MotherPi0_Pt_Alpha","ESD_MotherPi0_Pt_Alpha",350,0.03,35.,100,0,1);
			SetLogBinningXTH2(fHistoMotherPi0PtAlpha[iCut]);
			fESDList[iCut]->Add(fHistoMotherPi0PtAlpha[iCut]);
			fHistoMotherEtaPtAlpha[iCut] = new TH2F("ESD_MotherEta_Pt_Alpha","ESD_MotherEta_Pt_Alpha",350,0.03,35.,100,0,1);
			SetLogBinningXTH2(fHistoMotherEtaPtAlpha[iCut]);
			fESDList[iCut]->Add(fHistoMotherEtaPtAlpha[iCut]);
			fHistoMotherPi0PtOpenAngle[iCut] = new TH2F("ESD_MotherPi0_Pt_OpenAngle","ESD_MotherPi0_Pt_OpenAngle",350,0.03,35.,100,0, 0.5);
			SetLogBinningXTH2(fHistoMotherPi0PtOpenAngle[iCut]);
			fESDList[iCut]->Add(fHistoMotherPi0PtOpenAngle[iCut]);
			fHistoMotherEtaPtOpenAngle[iCut] = new TH2F("ESD_MotherEta_Pt_OpenAngle","ESD_MotherEta_Pt_OpenAngle",350,0.03,35.,180,0, 1.8);
			SetLogBinningXTH2(fHistoMotherEtaPtOpenAngle[iCut]);
			fESDList[iCut]->Add(fHistoMotherEtaPtOpenAngle[iCut]);
		}
		
		if (fIsMC == 2){
			fHistoMotherPi0PtY[iCut]->Sumw2();
			fHistoMotherEtaPtY[iCut]->Sumw2();
			fHistoMotherPi0PtAlpha[iCut]->Sumw2();
			fHistoMotherEtaPtAlpha[iCut]->Sumw2();
			fHistoMotherPi0PtOpenAngle[iCut]->Sumw2();
			fHistoMotherEtaPtOpenAngle[iCut]->Sumw2();
		}

	}
  
	if (fIsMC > 0){
		fMCList = new TList*[fnCuts];
		fTrueList = new TList*[fnCuts];

		fHistoMCPi0Pt	 		= new TH1F*[fnCuts];
		fHistoMCPi0WOWeightPt 	= new TH1F*[fnCuts];
		fHistoMCEtaPt 			= new TH1F*[fnCuts];
		fHistoMCEtaWOWeightPt	= new TH1F*[fnCuts];
		fHistoMCPi0InAccPt 		= new TH1F*[fnCuts];
		fHistoMCEtaInAccPt 		= new TH1F*[fnCuts];
		if (fIsMC == 2){
			fHistoMCPi0WOEvtWeightPt = new TH1F*[fnCuts];
			fHistoMCEtaWOEvtWeightPt = new TH1F*[fnCuts];
			if (fDoMesonQA > 0){
				fHistoMCPi0PtJetPt = new TH2F*[fnCuts];
				fHistoMCEtaPtJetPt = new TH2F*[fnCuts];
			}
		}

		if (fDoMesonQA > 0 && fDoMesonQA < 3){
			fHistoMCPi0PtY 		= new TH2F*[fnCuts];
			fHistoMCEtaPtY 		= new TH2F*[fnCuts];
			fHistoMCPi0PtAlpha 	= new TH2F*[fnCuts];
			fHistoMCEtaPtAlpha 	= new TH2F*[fnCuts];
		}

		fHistoTrueClusMergedPtvsM02 			= new TH2F*[fnCuts];
		fHistoTrueClusPi0PtvsM02 				= new TH2F*[fnCuts];
		fHistoTrueClusPrimPi0PtvsM02 			= new TH2F*[fnCuts];
		fHistoTrueClusSecPi0PtvsM02 			= new TH2F*[fnCuts];
		fHistoTrueClusSecPi0FromK0sPtvsM02 		= new TH2F*[fnCuts];
		fHistoTrueClusSecPi0FromLambdaPtvsM02 	= new TH2F*[fnCuts];
		fHistoTrueClusEtaPtvsM02 				= new TH2F*[fnCuts];
		fHistoTrueClusMergedPartConvPtvsM02 	= new TH2F*[fnCuts];
		fHistoTrueClusMergedPartConvELeadPtvsM02 		= new TH2F*[fnCuts];
		fHistoTrueClusPartConvPi0PtvsM02 		= new TH2F*[fnCuts];
		fHistoTrueClusPartConvPrimPi0PtvsM02 	= new TH2F*[fnCuts];
		fHistoTrueClusPartConvSecPi0PtvsM02 	= new TH2F*[fnCuts];
		fHistoTrueClusPartConvSecPi0FromK0sPtvsM02 		= new TH2F*[fnCuts];
		fHistoTrueClusPartConvSecPi0FromLambdaPtvsM02 	= new TH2F*[fnCuts];
		fHistoTrueClusPartConvEtaPtvsM02 		= new TH2F*[fnCuts];
		fHistoTrueClusBGPtvsM02 				= new TH2F*[fnCuts];
		fHistoTrueClusGammaPtvsM02 				= new TH2F*[fnCuts];

		fHistoTrueClusMergedInvMassvsPt 			= new TH2F*[fnCuts];
		fHistoTrueClusPi0InvMassvsPt 				= new TH2F*[fnCuts];
		fHistoTrueClusPrimPi0InvMassvsPt 			= new TH2F*[fnCuts];
		fHistoTrueClusSecPi0InvMassvsPt 			= new TH2F*[fnCuts];
		fHistoTrueClusSecPi0FromK0sInvMassvsPt 		= new TH2F*[fnCuts];
		fHistoTrueClusSecPi0FromLambdaInvMassvsPt 	= new TH2F*[fnCuts];
		fHistoTrueClusEtaInvMassvsPt 				= new TH2F*[fnCuts];
		fHistoTrueClusMergedPartConvInvMassvsPt 	= new TH2F*[fnCuts];
		fHistoTrueClusMergedPartConvELeadInvMassvsPt 		= new TH2F*[fnCuts];
		fHistoTrueClusPartConvPi0InvMassvsPt 		= new TH2F*[fnCuts];
		fHistoTrueClusPartConvPrimPi0InvMassvsPt 	= new TH2F*[fnCuts];
		fHistoTrueClusPartConvSecPi0InvMassvsPt 	= new TH2F*[fnCuts];
		fHistoTrueClusPartConvSecPi0FromK0sInvMassvsPt 		= new TH2F*[fnCuts];
		fHistoTrueClusPartConvSecPi0FromLambdaInvMassvsPt 	= new TH2F*[fnCuts];
		fHistoTrueClusPartConvEtaInvMassvsPt 		= new TH2F*[fnCuts];
		fHistoTrueClusBGInvMassvsPt 				= new TH2F*[fnCuts];
		fHistoTrueClusGammaInvMassvsPt 				= new TH2F*[fnCuts];
		fHistoTrueClusBGPtvsSource					= new TH2F*[fnCuts];
		
		if (fDoMesonQA > 0 && fDoMesonQA < 3){
			fHistoTruePi0PtY 							= new TH2F*[fnCuts];
			fHistoTrueEtaPtY 							= new TH2F*[fnCuts];
			fHistoTruePi0PtAlpha 						= new TH2F*[fnCuts];
			fHistoTrueEtaPtAlpha 						= new TH2F*[fnCuts];
			fHistoTruePi0PtOpenAngle 					= new TH2F*[fnCuts];
			fHistoTrueEtaPtOpenAngle 					= new TH2F*[fnCuts];
		}
		
		for(Int_t iCut = 0; iCut<fnCuts;iCut++){
			TString cutstringEvent 	= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
			TString cutstringCalo 	= ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
			TString cutstringCaloMerged 	= ((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))->GetCutNumber();
			TString cutstringMeson 	= ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
			
      
			fMCList[iCut] = new TList();
			fMCList[iCut]->SetName(Form("%s_%s_%s_%s MC histograms",cutstringEvent.Data(),cutstringCalo.Data(),cutstringCaloMerged.Data(), cutstringMeson.Data()));
			fMCList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fMCList[iCut]);

			fHistoMCPi0Pt[iCut] = new TH1F("MC_Pi0_Pt","MC_Pi0_Pt",500, 0, 50);
			fHistoMCPi0Pt[iCut]->Sumw2();
			fMCList[iCut]->Add(fHistoMCPi0Pt[iCut]);
			fHistoMCPi0WOWeightPt[iCut] = new TH1F("MC_Pi0_WOWeights_Pt","MC_Pi0_WOWeights_Pt",500, 0, 50);
			fHistoMCPi0WOWeightPt[iCut]->Sumw2();
			fMCList[iCut]->Add(fHistoMCPi0WOWeightPt[iCut]);
			
			fHistoMCEtaPt[iCut] = new TH1F("MC_Eta_Pt","MC_Eta_Pt",500, 0, 50);
			fHistoMCEtaPt[iCut]->Sumw2();
			fMCList[iCut]->Add(fHistoMCEtaPt[iCut]);
			fHistoMCEtaWOWeightPt[iCut] = new TH1F("MC_Eta_WOWeights_Pt","MC_Eta_WOWeights_Pt",500, 0, 50);
			fHistoMCEtaWOWeightPt[iCut]->Sumw2();
			fMCList[iCut]->Add(fHistoMCEtaWOWeightPt[iCut]);
			fHistoMCPi0InAccPt[iCut] = new TH1F("MC_Pi0InAcc_Pt","MC_Pi0InAcc_Pt",500, 0, 50);
			fHistoMCPi0InAccPt[iCut]->Sumw2();
			fMCList[iCut]->Add(fHistoMCPi0InAccPt[iCut]);
			fHistoMCEtaInAccPt[iCut] = new TH1F("MC_EtaInAcc_Pt","MC_EtaInAcc_Pt",500, 0, 50);
			fHistoMCEtaInAccPt[iCut]->Sumw2();
			fMCList[iCut]->Add(fHistoMCEtaInAccPt[iCut]);
			if (fIsMC == 2){
				fHistoMCPi0WOWeightPt[iCut]->Sumw2();
				fHistoMCEtaWOWeightPt[iCut]->Sumw2();
				fHistoMCPi0WOEvtWeightPt[iCut] = new TH1F("MC_Pi0_WOEventWeights_Pt","MC_Pi0_WOEventWeights_Pt",500, 0, 50);
				fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightPt[iCut]);
				fHistoMCEtaWOEvtWeightPt[iCut] = new TH1F("MC_Eta_WOEventWeights_Pt","MC_Eta_WOEventWeights_Pt",500, 0, 50);
				fMCList[iCut]->Add(fHistoMCEtaWOEvtWeightPt[iCut]);
				if (fDoMesonQA > 0){
					fHistoMCPi0PtJetPt[iCut] = new TH2F("MC_Pi0_Pt_JetPt","MC_Pi0_Pt_JetPt",350,0.03,35.,200,0,200);
					fHistoMCPi0PtJetPt[iCut]->Sumw2();
					SetLogBinningXTH2(fHistoMCPi0PtJetPt[iCut]);
					fMCList[iCut]->Add(fHistoMCPi0PtJetPt[iCut]);
					fHistoMCEtaPtJetPt[iCut] = new TH2F("MC_Eta_Pt_JetPt","MC_Eta_Pt_JetPt",350,0.03,35.,200,0,200);
					fHistoMCEtaPtJetPt[iCut]->Sumw2();
					SetLogBinningXTH2(fHistoMCEtaPtJetPt[iCut]);
					fMCList[iCut]->Add(fHistoMCEtaPtJetPt[iCut]);
				}
			}

			if (fDoMesonQA > 0 && fDoMesonQA < 3){
				fHistoMCPi0PtY[iCut] = new TH2F("MC_Pi0_Pt_Y","MC_Pi0_Pt_Y",350,0.03,35.,150,-1.5,1.5);
				fHistoMCPi0PtY[iCut]->Sumw2();
				SetLogBinningXTH2(fHistoMCPi0PtY[iCut]);
				fMCList[iCut]->Add(fHistoMCPi0PtY[iCut]);
				fHistoMCEtaPtY[iCut] = new TH2F("MC_Eta_Pt_Y","MC_Eta_Pt_Y",350,0.03,35.,150,-1.5,1.5);
				fHistoMCEtaPtY[iCut]->Sumw2();
				SetLogBinningXTH2(fHistoMCEtaPtY[iCut]);
				fMCList[iCut]->Add(fHistoMCEtaPtY[iCut]);
				fHistoMCPi0PtAlpha[iCut] = new TH2F("MC_Pi0_Pt_Alpha","MC_Pi0_Pt_Alpha",350,0.03,35.,100,0,1);
				SetLogBinningXTH2(fHistoMCPi0PtAlpha[iCut]);
				fMCList[iCut]->Add(fHistoMCPi0PtAlpha[iCut]);
				fHistoMCEtaPtAlpha[iCut] = new TH2F("MC_Eta_Pt_Alpha","MC_Eta_Pt_Alpha",350,0.03,35.,100,0,1);
				SetLogBinningXTH2(fHistoMCEtaPtAlpha[iCut]);
				fMCList[iCut]->Add(fHistoMCEtaPtAlpha[iCut]);

				if (fIsMC == 2){
					fHistoMCPi0PtAlpha[iCut]->Sumw2();
					fHistoMCEtaPtAlpha[iCut]->Sumw2();
				}		
			}
        
			fTrueList[iCut] = new TList();
			fTrueList[iCut]->SetName(Form("%s_%s_%s_%s True histograms",cutstringEvent.Data(),cutstringCalo.Data(),cutstringCaloMerged.Data(), cutstringMeson.Data()));
			fTrueList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fTrueList[iCut]);

			fHistoTrueClusMergedPtvsM02[iCut] = new TH2F("ESD_TrueClusMerged_Pt_M02","ESD_TrueClusMerged_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusMergedPtvsM02[iCut]); 
			fHistoTrueClusPi0PtvsM02 [iCut] = new TH2F("ESD_TrueClusFromPi0_Pt_M02","ESD_TrueClusFromPi0_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusPi0PtvsM02[iCut]);
			fHistoTrueClusPrimPi0PtvsM02[iCut] = new TH2F("ESD_TrueClusFromPrimPi0_Pt_M02","ESD_TrueClusFromPrimPi0_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusPrimPi0PtvsM02[iCut]);
			fHistoTrueClusSecPi0PtvsM02[iCut] = new TH2F("ESD_TrueClusFromSecPi0_Pt_M02","ESD_TrueClusFromSecPi0_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusSecPi0PtvsM02[iCut]);
			fHistoTrueClusSecPi0FromK0sPtvsM02[iCut]  = new TH2F("ESD_TrueClusFromSecPi0FromK0s_Pt_M02","ESD_TrueClusFromSecPi0FromK0s_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusSecPi0FromK0sPtvsM02[iCut]);
			fHistoTrueClusSecPi0FromLambdaPtvsM02[iCut]  = new TH2F("ESD_TrueClusFromSecPi0FromLambda_Pt_M02","ESD_TrueClusFromSecPi0FromLambda_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusSecPi0FromLambdaPtvsM02[iCut]);
			fHistoTrueClusEtaPtvsM02[iCut] = new TH2F("ESD_TrueClusFromEta_Pt_M02","ESD_TrueClusFromEta_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusEtaPtvsM02[iCut]);			
			fHistoTrueClusMergedPartConvPtvsM02[iCut]  = new TH2F("ESD_TrueClusMergedPartConv_Pt_M02","ESD_TrueClusMergedPartConv_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusMergedPartConvPtvsM02[iCut]);
			fHistoTrueClusMergedPartConvELeadPtvsM02[iCut]  = new TH2F("ESD_TrueClusMergedPartConvLeadE_Pt_M02","ESD_TrueClusMergedPartConvLeadE_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusMergedPartConvELeadPtvsM02[iCut]);
			fHistoTrueClusPartConvPi0PtvsM02[iCut] 	= new TH2F("ESD_TrueClusPartConvFromPi0_Pt_M02","ESD_TrueClusPartConvFromPi0_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusPartConvPi0PtvsM02[iCut]);
			fHistoTrueClusPartConvPrimPi0PtvsM02[iCut] 	= new TH2F("ESD_TrueClusPartConvFromPrimPi0_Pt_M02","ESD_TrueClusPartConvFromPrimPi0_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusPartConvPrimPi0PtvsM02[iCut]);
			fHistoTrueClusPartConvSecPi0PtvsM02[iCut] 	= new TH2F("ESD_TrueClusPartConvFromSecPi0_Pt_M02","ESD_TrueClusPartConvFromSecPi0_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusPartConvSecPi0PtvsM02[iCut]);
			fHistoTrueClusPartConvSecPi0FromK0sPtvsM02[iCut] 	= new TH2F("ESD_TrueClusPartConvFromSecPi0FromK0s_Pt_M02","ESD_TrueClusPartConvFromSecPi0FromK0s_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusPartConvSecPi0FromK0sPtvsM02[iCut]);
			fHistoTrueClusPartConvSecPi0FromLambdaPtvsM02[iCut] = new TH2F("ESD_TrueClusPartConvFromSecPi0FromLamdba_Pt_M02","ESD_TrueClusPartConvFromSecPi0FromLamdba_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusPartConvSecPi0FromLambdaPtvsM02[iCut]);
			fHistoTrueClusPartConvEtaPtvsM02[iCut] 	= new TH2F("ESD_TrueClusPartConvFromEta_Pt_M02","ESD_TrueClusPartConvFromEta_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusPartConvEtaPtvsM02[iCut]);
			fHistoTrueClusBGPtvsM02[iCut]  = new TH2F("ESD_TrueClusBG_Pt_M02","ESD_TrueClusBG_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusBGPtvsM02[iCut]);
			fHistoTrueClusGammaPtvsM02[iCut] = new TH2F("ESD_TrueClusGamma_Pt_M02","ESD_TrueClusGamma_Pt_M02",500,0,50,500,0,5);
			fTrueList[iCut]->Add(fHistoTrueClusGammaPtvsM02[iCut]);			

			fHistoTrueClusMergedInvMassvsPt[iCut] = new TH2F("ESD_TrueClusMerged_InvMass_Pt","ESD_TrueClusMerged_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusMergedInvMassvsPt[iCut]); 
			fHistoTrueClusPi0InvMassvsPt [iCut] = new TH2F("ESD_TrueClusFromPi0_InvMass_Pt","ESD_TrueClusFromPi0_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusPi0InvMassvsPt[iCut]);
			fHistoTrueClusPrimPi0InvMassvsPt[iCut] = new TH2F("ESD_TrueClusFromPrimPi0_InvMass_Pt","ESD_TrueClusFromPrimPi0_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusPrimPi0InvMassvsPt[iCut]);
			fHistoTrueClusSecPi0InvMassvsPt[iCut] = new TH2F("ESD_TrueClusFromSecPi0_InvMass_Pt","ESD_TrueClusFromSecPi0_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusSecPi0InvMassvsPt[iCut]);
			fHistoTrueClusSecPi0FromK0sInvMassvsPt[iCut]  = new TH2F("ESD_TrueClusFromSecPi0FromK0s_InvMass_Pt","ESD_TrueClusFromSecPi0FromK0s_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusSecPi0FromK0sInvMassvsPt[iCut]);
			fHistoTrueClusSecPi0FromLambdaInvMassvsPt[iCut]  = new TH2F("ESD_TrueClusFromSecPi0FromLambda_InvMass_Pt","ESD_TrueClusFromSecPi0FromLambda_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusSecPi0FromLambdaInvMassvsPt[iCut]);
			fHistoTrueClusEtaInvMassvsPt[iCut] = new TH2F("ESD_TrueClusFromEta_InvMass_Pt","ESD_TrueClusFromEta_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusEtaInvMassvsPt[iCut]);			
			fHistoTrueClusMergedPartConvInvMassvsPt[iCut]  = new TH2F("ESD_TrueClusMergedPartConv_InvMass_Pt","ESD_TrueClusMergedPartConv_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusMergedPartConvInvMassvsPt[iCut]);
			fHistoTrueClusMergedPartConvELeadInvMassvsPt[iCut]  = new TH2F("ESD_TrueClusMergedPartConvLeadE_InvMass_Pt","ESD_TrueClusMergedPartConvLeadE_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusMergedPartConvELeadInvMassvsPt[iCut]);
			fHistoTrueClusPartConvPi0InvMassvsPt[iCut] 	= new TH2F("ESD_TrueClusPartConvFromPi0_InvMass_Pt","ESD_TrueClusPartConvFromPi0_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusPartConvPi0InvMassvsPt[iCut]);
			fHistoTrueClusPartConvPrimPi0InvMassvsPt[iCut] 	= new TH2F("ESD_TrueClusPartConvFromPrimPi0_InvMass_Pt","ESD_TrueClusPartConvFromPrimPi0_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusPartConvPrimPi0InvMassvsPt[iCut]);
			fHistoTrueClusPartConvSecPi0InvMassvsPt[iCut] 	= new TH2F("ESD_TrueClusPartConvFromSecPi0_InvMass_Pt","ESD_TrueClusPartConvFromSecPi0_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusPartConvSecPi0InvMassvsPt[iCut]);
			fHistoTrueClusPartConvSecPi0FromK0sInvMassvsPt[iCut] 	= new TH2F("ESD_TrueClusPartConvFromSecPi0FromK0s_InvMass_Pt","ESD_TrueClusPartConvFromSecPi0FromK0s_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusPartConvSecPi0FromK0sInvMassvsPt[iCut]);
			fHistoTrueClusPartConvSecPi0FromLambdaInvMassvsPt[iCut] = new TH2F("ESD_TrueClusPartConvFromSecPi0FromLamdba_InvMass_Pt","ESD_TrueClusPartConvFromSecPi0FromLamdba_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusPartConvSecPi0FromLambdaInvMassvsPt[iCut]);
			fHistoTrueClusPartConvEtaInvMassvsPt[iCut] 	= new TH2F("ESD_TrueClusPartConvFromEta_InvMass_Pt","ESD_TrueClusPartConvFromEta_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusPartConvEtaInvMassvsPt[iCut]);
			fHistoTrueClusBGInvMassvsPt[iCut]  = new TH2F("ESD_TrueClusBG_InvMass_Pt","ESD_TrueClusBG_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusBGInvMassvsPt[iCut]);
			fHistoTrueClusGammaInvMassvsPt[iCut] = new TH2F("ESD_TrueClusGamma_InvMass_Pt","ESD_TrueClusGamma_InvMass_Pt",800,0,0.8,500, 0, 50);
			fTrueList[iCut]->Add(fHistoTrueClusGammaInvMassvsPt[iCut]);			

			fHistoTrueClusBGPtvsSource[iCut] = new TH2F("ESD_TrueClusBG_Pt_Source","ESD_TrueClusBG_Pt_Source",500, 0, 50,10, 0, 10);
			fTrueList[iCut]->Add(fHistoTrueClusBGPtvsSource[iCut]);			
			
			if (fDoMesonQA > 0 && fDoMesonQA < 3){
				fHistoTruePi0PtY[iCut] = new TH2F("ESD_TruePi0_Pt_Y","ESD_TruePi0_Pt_Y",350,0.03,35.,150,-1.5,1.5);
				SetLogBinningXTH2(fHistoTruePi0PtY[iCut]);
				fTrueList[iCut]->Add(fHistoTruePi0PtY[iCut]);
				fHistoTrueEtaPtY[iCut] = new TH2F("ESD_TrueEta_Pt_Y","ESD_TrueEta_Pt_Y",350,0.03,35.,150,-1.5,1.5);
				SetLogBinningXTH2(fHistoTrueEtaPtY[iCut]);
				fTrueList[iCut]->Add(fHistoTrueEtaPtY[iCut]);
				fHistoTruePi0PtAlpha[iCut] = new TH2F("ESD_TruePi0_Pt_Alpha","ESD_TruePi0_Pt_Alpha",350,0.03,35.,100,0,1);
				SetLogBinningXTH2(fHistoTruePi0PtAlpha[iCut]);
				fTrueList[iCut]->Add(fHistoTruePi0PtAlpha[iCut]);
				fHistoTrueEtaPtAlpha[iCut] = new TH2F("ESD_TrueEta_Pt_Alpha","ESD_TrueEta_Pt_Alpha",350,0.03,35.,100,0,1);
				SetLogBinningXTH2(fHistoTrueEtaPtAlpha[iCut]);
				fTrueList[iCut]->Add(fHistoTrueEtaPtAlpha[iCut]);			
				fHistoTruePi0PtOpenAngle[iCut] = new TH2F("ESD_TruePi0_Pt_OpenAngle","ESD_TruePi0_Pt_OpenAngle",350,0.03,35.,100,0,0.5);
				SetLogBinningXTH2(fHistoTruePi0PtOpenAngle[iCut]);
				fTrueList[iCut]->Add(fHistoTruePi0PtOpenAngle[iCut]);
				fHistoTrueEtaPtOpenAngle[iCut] = new TH2F("ESD_TrueEta_Pt_OpenAngle","ESD_TrueEta_Pt_OpenAngle",350,0.03,35.,180,0,1.8);
				SetLogBinningXTH2(fHistoTrueEtaPtOpenAngle[iCut]);
				fTrueList[iCut]->Add(fHistoTrueEtaPtOpenAngle[iCut]);
			}
			
			if (fIsMC == 2){
				fHistoTrueClusMergedPtvsM02[iCut]->Sumw2();
				fHistoTrueClusPi0PtvsM02[iCut]->Sumw2();
				fHistoTrueClusPrimPi0PtvsM02[iCut]->Sumw2();
				fHistoTrueClusSecPi0PtvsM02[iCut]->Sumw2();
				fHistoTrueClusSecPi0FromK0sPtvsM02[iCut]->Sumw2();
				fHistoTrueClusSecPi0FromLambdaPtvsM02[iCut]->Sumw2();
				fHistoTrueClusEtaPtvsM02[iCut]->Sumw2();
				fHistoTrueClusMergedPartConvPtvsM02[iCut]->Sumw2();
				fHistoTrueClusMergedPartConvELeadPtvsM02[iCut]->Sumw2();
				fHistoTrueClusPartConvPi0PtvsM02[iCut]->Sumw2();
				fHistoTrueClusPartConvPrimPi0PtvsM02[iCut]->Sumw2();
				fHistoTrueClusPartConvSecPi0PtvsM02[iCut]->Sumw2();
				fHistoTrueClusPartConvSecPi0FromK0sPtvsM02[iCut]->Sumw2();
				fHistoTrueClusPartConvSecPi0FromLambdaPtvsM02[iCut]->Sumw2();
				fHistoTrueClusPartConvEtaPtvsM02[iCut]->Sumw2();
				fHistoTrueClusBGPtvsM02[iCut]->Sumw2();
				fHistoTrueClusGammaPtvsM02[iCut]->Sumw2();

				fHistoTrueClusMergedInvMassvsPt[iCut]->Sumw2();
				fHistoTrueClusPi0InvMassvsPt[iCut]->Sumw2();
				fHistoTrueClusPrimPi0InvMassvsPt[iCut]->Sumw2();
				fHistoTrueClusSecPi0InvMassvsPt[iCut]->Sumw2();
				fHistoTrueClusSecPi0FromK0sInvMassvsPt[iCut]->Sumw2();
				fHistoTrueClusSecPi0FromLambdaInvMassvsPt[iCut]->Sumw2();
				fHistoTrueClusEtaInvMassvsPt[iCut]->Sumw2();
				fHistoTrueClusMergedPartConvInvMassvsPt[iCut]->Sumw2();
				fHistoTrueClusMergedPartConvELeadInvMassvsPt[iCut]->Sumw2();
				fHistoTrueClusPartConvPi0InvMassvsPt[iCut]->Sumw2();
				fHistoTrueClusPartConvPrimPi0InvMassvsPt[iCut]->Sumw2();
				fHistoTrueClusPartConvSecPi0InvMassvsPt[iCut]->Sumw2();
				fHistoTrueClusPartConvSecPi0FromK0sInvMassvsPt[iCut]->Sumw2();
				fHistoTrueClusPartConvSecPi0FromLambdaInvMassvsPt[iCut]->Sumw2();
				fHistoTrueClusPartConvEtaInvMassvsPt[iCut]->Sumw2();
				fHistoTrueClusBGInvMassvsPt[iCut]->Sumw2();
				fHistoTrueClusGammaInvMassvsPt[iCut]->Sumw2();

				if (fDoMesonQA > 0 && fDoMesonQA < 3){
					fHistoTruePi0PtY[iCut]->Sumw2();
					fHistoTrueEtaPtY[iCut]->Sumw2();
					fHistoTruePi0PtAlpha[iCut]->Sumw2();
					fHistoTrueEtaPtAlpha[iCut]->Sumw2();
					fHistoTruePi0PtOpenAngle[iCut]->Sumw2();
					fHistoTrueEtaPtOpenAngle[iCut]->Sumw2();
				}
			}
		}
		
	}	
  

	fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
	if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader
  
	if(fV0Reader)
		if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
			if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
				fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
	if(fV0Reader)
		if((AliConvEventCuts*)fV0Reader->GetEventCuts())
			if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
				fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());

			
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
		if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
			fCutFolder[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
		}
		if(!((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))) continue;
		if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms()){
			fCutFolder[iCut]->Add(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms());
		}
		if(!((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))) continue;
		if(((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))->GetCutHistograms()){
			fCutFolder[iCut]->Add(((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))->GetCutHistograms());
		}

		if(fSetPlotHistsExtQA){
			if(!((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))) continue;
			if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetExtQAHistograms()){
				fCutFolder[iCut]->Add(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetExtQAHistograms());
			}
		}
		if(!((AliConversionMesonCuts*)fMesonCutArray->At(iCut))) continue;
		if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms()){
			fCutFolder[iCut]->Add(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
		}
	}
	PostData(1, fOutputContainer);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaCaloMerged::Notify()
{
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		if(!((AliConvEventCuts*)fEventCutArray->At(iCut))->GetDoEtaShift()){
			fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
			continue; // No Eta Shift requested, continue
		}
		if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
			((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCorrectEtaShiftFromPeriod(fV0Reader->GetPeriodName());
			fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
			((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
			continue;
		}
		else{
			printf(" Gamma Conversion Task %s :: Eta Shift Manually Set to %f \n\n",
					(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber()).Data(),((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift());
			fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
			((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
		}
	}
	return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::UserExec(Option_t *)
{
	//
	// Called for each event
	//
	Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
	if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1
		for(Int_t iCut = 0; iCut<fnCuts; iCut++){
			fHistoNEvents[iCut]->Fill(eventQuality);
			if (fIsMC==2) 
				fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
		}
		return;
	}
	
	if(fIsMC> 0) fMCEvent = MCEvent();
	if(fMCEvent == NULL) fIsMC = 0;
	
	fInputEvent = InputEvent();
	
	if(fIsMC>0 && fInputEvent->IsA()==AliESDEvent::Class()){
		fMCStack = fMCEvent->Stack();
		if(fMCStack == NULL) fIsMC = 0;
	}
	
	// ------------------- BeginEvent ----------------------------
		
	for(Int_t iCut = 0; iCut<fnCuts; iCut++){	
		fiCut = iCut;
		
		Bool_t isRunningEMCALrelAna = kFALSE;
		if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;
		
		Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon, isRunningEMCALrelAna);
		
		fWeightJetJetMC = 1;
		// 		cout << fMCEvent << endl;
		Bool_t isMCJet = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC );
		if (!isMCJet){
			fHistoNEvents[iCut]->Fill(10,
									  fWeightJetJetMC);
			if (fIsMC==2) 
				fHistoNEventsWOWeight[iCut]->Fill(10);
			continue;
		}

		Bool_t triggered = kTRUE;
		if(eventNotAccepted){
		// cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
			fHistoNEvents[iCut]->Fill(eventNotAccepted, fWeightJetJetMC); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
			if (fIsMC==2) fHistoNEventsWOWeight[iCut]->Fill(eventNotAccepted);
			if (eventNotAccepted==3 && fIsMC>0){
				triggered = kFALSE;
			} else {	
				continue;
			}	
		}

		if(eventQuality != 0){// Event Not Accepted
			//cout << "event rejected due to: " <<eventQuality << endl;
			fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC);
			if (fIsMC==2) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);

			continue;
		}
		if (triggered == kTRUE) {
			fHistoNEvents[iCut]->Fill(	eventQuality, 
										fWeightJetJetMC); // Should be 0 here
			if (fIsMC==2){
				fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here
				Float_t xsection = -1.; Float_t ntrials = -1.;
				((AliConvEventCuts*)fEventCutArray->At(iCut))->GetXSectionAndNTrials(fMCEvent,xsection,ntrials);
				if((xsection==-1.) || (ntrials==-1.)) AliFatal("ERROR: GetXSectionAndNTrials returned invalid xsection/ntrials");
				fProfileJetJetXSection[iCut]->Fill(0.,xsection);
				fHistoJetJetNTrials[iCut]->Fill("#sum{NTrials}",ntrials);
			}

			fHistoNGoodESDTracks[iCut]->Fill(	fV0Reader->GetNumberOfPrimaryTracks(), 
												fWeightJetJetMC);
			fHistoVertexZ[iCut]->Fill(	fInputEvent->GetPrimaryVertex()->GetZ(), 
										fWeightJetJetMC);
			fHistoSPDClusterTrackletBackground[iCut]->Fill(	fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),
															(fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)), 
															fWeightJetJetMC);
			if(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsHeavyIon() == 2)
				fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A(),
											fWeightJetJetMC);
			else 
				fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C(),
											fWeightJetJetMC);
		}
		if(fIsMC> 0){
			// Process MC Particle
			if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection() != 0){
				if(fInputEvent->IsA()==AliESDEvent::Class()){
				((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(),
																					((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
																					fMCEvent);
				}
				else if(fInputEvent->IsA()==AliAODEvent::Class()){
				((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(),
																					((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
																					fInputEvent);
				}

				if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader()){
					for(Int_t i = 0;i<(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader())->GetEntries();i++){
						TString nameBin= fHistoMCHeaders[iCut]->GetXaxis()->GetBinLabel(i+1);
						if (nameBin.CompareTo("")== 0){
							TString nameHeader = ((TObjString*)((TList*)((AliConvEventCuts*)fEventCutArray->At(iCut))
																->GetAcceptedHeader())->At(i))->GetString();
							fHistoMCHeaders[iCut]->GetXaxis()->SetBinLabel(i+1,nameHeader.Data());
						}
					}
				}
			}
		}
		if(fIsMC> 0)
			ProcessMCParticles();
		
		if (triggered==kFALSE) continue;
		
		// it is in the loop to have the same conversion cut string (used also for MC stuff that should be same for V0 and Cluster)
		ProcessClusters();  					// process calo clusters

		fHistoNClusterCandidates[iCut]->Fill(	fNClusterCandidates, 
												fWeightJetJetMC);
		fHistoNClusterMergedCandidates[iCut]->Fill(	fNClusterMergedCandidates, 
												fWeightJetJetMC);

		fHistoNGoodESDTracksVsNClusterCandidates[iCut]->Fill(	fV0Reader->GetNumberOfPrimaryTracks(),
																fNClusterCandidates,
																fWeightJetJetMC);
// 	
// 			fVectorDoubleCountTruePi0s.clear();
// 			fVectorDoubleCountTrueEtas.clear();

	}
	
	PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::ProcessClusters(){
	
	Int_t nclus = 0;
	nclus = fInputEvent->GetNumberOfCaloClusters();
	
// 	cout << nclus << endl;
	
	if(nclus == 0)	return;
	
	fNClusterCandidates = 0;
	fNClusterMergedCandidates = 0;
	
	// plotting histograms on cell/tower level, only if extendedMatchAndQA > 1
	((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillHistogramsExtendedQA(fInputEvent);

	// vertex
	Double_t vertex[3] = {0};
	InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
	
	// Loop over EMCal clusters
	for(Long_t i = 0; i < nclus; i++){
	
		// select clusters with open cuts for normalizations histos
		AliVCluster* clus = NULL;
		if(fInputEvent->IsA()==AliESDEvent::Class()) clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
		else if(fInputEvent->IsA()==AliAODEvent::Class()) clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));

		if (!clus){
			delete clus;
			continue;
		}
		
		// if open cluster cuts are not fullfilled I can abort
		if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fIsMC)){
			delete clus;
			continue;
		}
		
		// TLorentzvector with cluster	
		TLorentzVector clusterVector;
		clus->GetMomentum(clusterVector,vertex);
			
		TLorentzVector* tmpvec = new TLorentzVector();
		tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());
	
		// convert to AODConversionPhoton
		AliAODConversionPhoton *PhotonCandidate=new AliAODConversionPhoton(tmpvec);
		if(!PhotonCandidate){
			delete clus;
			delete tmpvec;
			if (PhotonCandidate) delete PhotonCandidate;
			continue;
		}
		
		// Flag Photon as CaloPhoton
		PhotonCandidate->SetIsCaloPhoton();
		PhotonCandidate->SetCaloClusterRef(i);
		
		// get MC label
		if(fIsMC> 0){
			Int_t* mclabelsCluster = clus->GetLabels();
			PhotonCandidate->SetNCaloPhotonMCLabels(clus->GetNLabels());

			if (clus->GetNLabels()>0){
				for (Int_t k =0; k< (Int_t)clus->GetNLabels(); k++){
					if (k< 20)PhotonCandidate->SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
				}	
			}
		}
		
		fIsFromMBHeader = kTRUE; 
		fIsOverlappingWithOtherHeader = kFALSE;

		// test whether largest contribution to cluster orginates in added signals
		if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() > 0){
			if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetCaloPhotonMCLabel(0), fMCStack, fInputEvent) == 0) fIsFromMBHeader = kFALSE;
			if (clus->GetNLabels()>1){
				Int_t* mclabelsCluster = clus->GetLabels();
				for (Int_t l = 1; l < (Int_t)clus->GetNLabels(); l++ ){
					if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(mclabelsCluster[l], fMCStack, fInputEvent) == 0) fIsOverlappingWithOtherHeader = kTRUE;
				}	
			}
		}
		
		// check whether photon is from correct header 
		if (fIsFromMBHeader && fIsOverlappingWithOtherHeader)
			fHistoClusOverlapHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), 
														 fWeightJetJetMC);
			
		// only cluster candidates from correct header will be processed fully	
		if (fIsFromMBHeader && !fIsOverlappingWithOtherHeader){
			fHistoClusGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), 
										   fWeightJetJetMC);	
			fNClusterCandidates++;
		} else {
			delete clus;
			delete tmpvec;
			delete PhotonCandidate;
			continue;
		}
		
		
		// check whether photon fullfill merged cluster cuts as well
		if(!((AliCaloPhotonCuts*)fClusterMergedCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fIsMC)){
			delete clus;
			delete tmpvec;
			delete PhotonCandidate;
			continue;
		}	else {
			fNClusterMergedCandidates++; 	
		}	

		AliAODCaloCluster* clusSub1 = new AliAODCaloCluster();
		AliAODCaloCluster* clusSub2 = new AliAODCaloCluster();
		
		
		// split clusters according to their shares in the cluster (NLM == 1) needs to be treated differently
		if (((AliCaloPhotonCuts*)fClusterMergedCutArray->At(fiCut))->GetMinNLMCut() == 1 && 
			((AliCaloPhotonCuts*)fClusterMergedCutArray->At(fiCut))->GetMaxNLMCut() == 1 ){
			
			Int_t absCellIdFirst 	= ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FindLargestCellInCluster(clus, fInputEvent);
			Int_t absCellIdSecond 	= ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FindSecondLargestCellInCluster(clus, fInputEvent);
			
			((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->SplitEnergy(absCellIdFirst, absCellIdSecond,
																		   clus, fInputEvent, fIsMC, clusSub1, clusSub2);
			
			
		} else if (((AliCaloPhotonCuts*)fClusterMergedCutArray->At(fiCut))->GetMinNLMCut() > 1 ){
			const Int_t   nc = clus->GetNCells();		
			Int_t   absCellIdList[nc]; 
			Float_t maxEList[nc]; 		
			Int_t nlm = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetNumberOfLocalMaxima(clus, fInputEvent, absCellIdList, maxEList);
			((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->SplitEnergy(absCellIdList[0], absCellIdList[1],
																		clus, fInputEvent, fIsMC, clusSub1, clusSub2);
		}
			// 			cout << clusSub1->E() << "\t" <<clusSub2->E() << endl; 
			
		fHistoClusMergedPtvsM02[fiCut]->Fill( PhotonCandidate->Pt(),
													clus->GetM02(),
													fWeightJetJetMC);	
			
		// TLorentzvector with sub cluster 1
		TLorentzVector clusterVector1;
		clusSub1->GetMomentum(clusterVector1,vertex);
		TLorentzVector* tmpvec1 = new TLorentzVector();
		tmpvec1->SetPxPyPzE(clusterVector1.Px(),clusterVector1.Py(),clusterVector1.Pz(),clusterVector1.E());
		// convert to AODConversionPhoton
		AliAODConversionPhoton *PhotonCandidate1=new AliAODConversionPhoton(tmpvec1);
		if(!PhotonCandidate){
			delete clusSub1;
			delete tmpvec1;
			if (PhotonCandidate1) delete PhotonCandidate1;
			continue;
		}
		// Flag Photon as CaloPhoton
		PhotonCandidate1->SetIsCaloPhoton();
		// get MC label
		if(fIsMC> 0){
			Int_t* mclabelsCluster = clusSub1->GetLabels();
			PhotonCandidate1->SetNCaloPhotonMCLabels(clusSub1->GetNLabels());

			if (clusSub1->GetNLabels()>0){
				for (Int_t k =0; k< (Int_t)clusSub1->GetNLabels(); k++){
					if (k< 20)PhotonCandidate1->SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
				}	
			}
		}
		
		// TLorentzvector with sub cluster 2
		TLorentzVector clusterVector2;
		clusSub2->GetMomentum(clusterVector2,vertex);
		TLorentzVector* tmpvec2 = new TLorentzVector();
		tmpvec2->SetPxPyPzE(clusterVector2.Px(),clusterVector2.Py(),clusterVector2.Pz(),clusterVector2.E());
		// convert to AODConversionPhoton
		AliAODConversionPhoton *PhotonCandidate2=new AliAODConversionPhoton(tmpvec2);
		if(!PhotonCandidate){
			delete clusSub2;
			delete tmpvec2;
			if (PhotonCandidate2) delete PhotonCandidate2;
			continue;
		}
		// Flag Photon as CaloPhoton
		PhotonCandidate2->SetIsCaloPhoton();
		// get MC label
		if(fIsMC> 0){
			Int_t* mclabelsCluster = clusSub2->GetLabels();
			PhotonCandidate2->SetNCaloPhotonMCLabels(clusSub2->GetNLabels());

			if (clusSub2->GetNLabels()>0){
				for (Int_t k =0; k< (Int_t)clusSub2->GetNLabels(); k++){
					if (k< 20)PhotonCandidate2->SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
				}	
			}
		}
		
		CalculatePi0Candidate( PhotonCandidate1, PhotonCandidate2);
		
		if(fIsMC> 0 && PhotonCandidate && PhotonCandidate1 && PhotonCandidate2){
			ProcessTrueClusterCandidates(PhotonCandidate, clus->GetM02(), PhotonCandidate1, PhotonCandidate2);
		}
					
		delete clusSub1;
		delete tmpvec1;
		delete PhotonCandidate1;
		delete clusSub2;			
		delete tmpvec2;		
		delete PhotonCandidate2;		
			
		delete clus;
		delete tmpvec;
		delete PhotonCandidate;
		

	}
	
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::ProcessTrueClusterCandidates(AliAODConversionPhoton *TrueClusterCandidate, Float_t m02,
															      AliAODConversionPhoton *TrueSubClusterCandidate1,
															      AliAODConversionPhoton *TrueSubClusterCandidate2)
{

	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();

	TParticle *Photon = NULL;
	if (!TrueClusterCandidate->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set task will abort");
	
	if (TrueClusterCandidate->GetNCaloPhotonMCLabels()>0) Photon = fMCStack->Particle(TrueClusterCandidate->GetCaloPhotonMCLabel(0));
		else return;
		
	if(Photon == NULL){
	//    cout << "no photon" << endl;
		return;
	}

	AliAODConversionMother *mesoncand = new AliAODConversionMother(TrueSubClusterCandidate1,TrueSubClusterCandidate2);
	Bool_t mesonIsSelected = (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(mesoncand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()));
	
	Int_t pdgCodeParticle = Photon->GetPdgCode();
	TrueClusterCandidate->SetCaloPhotonMCFlags(fMCStack);
		

	Bool_t filled = kFALSE;
	
	if(fIsFromMBHeader && !fIsOverlappingWithOtherHeader){	
		Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, TrueClusterCandidate->GetCaloPhotonMCLabel(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
		if (TrueClusterCandidate->IsMerged() && !TrueClusterCandidate->IsMergedPartConv()){
			filled = kTRUE;
			fHistoTrueClusMergedPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
			if (mesonIsSelected) fHistoTrueClusMergedInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
			Int_t motherLab = Photon->GetMother(0);
			if (motherLab > -1){ 
				if (fMCStack->Particle(motherLab)->GetPdgCode() == 111){
					fHistoTrueClusPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
					if (mesonIsSelected){
						fHistoTrueClusPi0InvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
						if ( mesoncand->M() > 0.0 && mesoncand->M() < 0.4 && fDoMesonQA > 0 && fDoMesonQA < 3){
							fHistoTruePi0PtY[fiCut]->Fill(mesoncand->Pt(),mesoncand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
							fHistoTruePi0PtAlpha[fiCut]->Fill(mesoncand->Pt(),abs(mesoncand->GetAlpha()), fWeightJetJetMC);
							fHistoTruePi0PtOpenAngle[fiCut]->Fill(mesoncand->Pt(),mesoncand->GetOpeningAngle(), fWeightJetJetMC);
						}
					}	

					if (isPrimary) {
						fHistoTrueClusPrimPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
						if (mesonIsSelected) fHistoTrueClusPrimPi0InvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
					} else {
						fHistoTrueClusSecPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
						if (mesonIsSelected) fHistoTrueClusSecPi0InvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
						Int_t grandMaLab = fMCStack->Particle(motherLab)->GetMother(0);
						if (grandMaLab > -1){
							if (abs(fMCStack->Particle(grandMaLab)->GetPdgCode()) == 310){
								fHistoTrueClusSecPi0FromK0sPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
								if (mesonIsSelected) fHistoTrueClusSecPi0FromK0sInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
								
							} else if (abs(fMCStack->Particle(grandMaLab)->GetPdgCode()) == 3122){
								fHistoTrueClusSecPi0FromLambdaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
								if (mesonIsSelected) fHistoTrueClusSecPi0FromLambdaInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
								
							} 	
						}	
					}	
				}
				if 	(fMCStack->Particle(motherLab)->GetPdgCode() == 221){
					fHistoTrueClusEtaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
					if (mesonIsSelected){
						fHistoTrueClusEtaInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
						if ( mesoncand->M() > 0.4 && mesoncand->M() < 0.8 && fDoMesonQA > 0 && fDoMesonQA < 3){
							fHistoTrueEtaPtY[fiCut]->Fill(mesoncand->Pt(),mesoncand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
							fHistoTrueEtaPtAlpha[fiCut]->Fill(mesoncand->Pt(),abs(mesoncand->GetAlpha()), fWeightJetJetMC);
							fHistoTrueEtaPtOpenAngle[fiCut]->Fill(mesoncand->Pt(),mesoncand->GetOpeningAngle(), fWeightJetJetMC);						
						}
					}
				}	
			}	
		}
		if (TrueClusterCandidate->IsMergedPartConv() && !filled){
			filled = kTRUE;
			fHistoTrueClusMergedPartConvPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
			if (mesonIsSelected) fHistoTrueClusMergedPartConvInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
			
			if (TrueClusterCandidate->IsLargestComponentElectron()){
				fHistoTrueClusMergedPartConvELeadPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
				if (mesonIsSelected) fHistoTrueClusMergedPartConvELeadInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
			}
			Int_t motherLab = Photon->GetMother(0);
			if (motherLab > -1){ 
				if (fMCStack->Particle(motherLab)->GetPdgCode() == 111){
					fHistoTrueClusPartConvPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
					if (mesonIsSelected){
						fHistoTrueClusPartConvPi0InvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
						if ( mesoncand->M() > 0.0 && mesoncand->M() < 0.4 && fDoMesonQA > 0 && fDoMesonQA < 3){
							fHistoTruePi0PtY[fiCut]->Fill(mesoncand->Pt(),mesoncand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
							fHistoTruePi0PtAlpha[fiCut]->Fill(mesoncand->Pt(),abs(mesoncand->GetAlpha()), fWeightJetJetMC);
							fHistoTruePi0PtOpenAngle[fiCut]->Fill(mesoncand->Pt(),mesoncand->GetOpeningAngle(), fWeightJetJetMC);						
						}	
					}	
					if (isPrimary) {
						fHistoTrueClusPartConvPrimPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
						if (mesonIsSelected) fHistoTrueClusPartConvPrimPi0InvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
					} else {
						fHistoTrueClusPartConvSecPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
						if (mesonIsSelected) fHistoTrueClusPartConvSecPi0InvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
						Int_t grandMaLab = fMCStack->Particle(motherLab)->GetMother(0);
						if (grandMaLab > -1){
							if (abs(fMCStack->Particle(grandMaLab)->GetPdgCode()) == 310){
								fHistoTrueClusPartConvSecPi0FromK0sPtvsM02[fiCut]->Fill( TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
								if (mesonIsSelected) fHistoTrueClusPartConvSecPi0FromK0sInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
							} else if (abs(fMCStack->Particle(grandMaLab)->GetPdgCode()) == 3122){
								fHistoTrueClusPartConvSecPi0FromLambdaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
								if (mesonIsSelected) fHistoTrueClusPartConvSecPi0FromLambdaInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
							} 	
						}	
					}	
				}
				if 	(fMCStack->Particle(motherLab)->GetPdgCode() == 221){
					fHistoTrueClusPartConvEtaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
					if (mesonIsSelected){
						fHistoTrueClusPartConvEtaInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
						if ( mesoncand->M() > 0.4 && mesoncand->M() < 0.8 && fDoMesonQA > 0 && fDoMesonQA < 3){
							fHistoTrueEtaPtY[fiCut]->Fill(mesoncand->Pt(),mesoncand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
							fHistoTrueEtaPtAlpha[fiCut]->Fill(mesoncand->Pt(),abs(mesoncand->GetAlpha()), fWeightJetJetMC);
							fHistoTrueEtaPtOpenAngle[fiCut]->Fill(mesoncand->Pt(),mesoncand->GetOpeningAngle(), fWeightJetJetMC);						
						}
					}	
				}	
			}	
		}	
		
		if (!(TrueClusterCandidate->IsLargestComponentElectron() || TrueClusterCandidate->IsLargestComponentPhoton())){
			fHistoTrueClusBGPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
			if (mesonIsSelected){
				fHistoTrueClusBGInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
				if (pdgCodeParticle== abs(211)) fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 0.5, fWeightJetJetMC); // pi+/-
				else if (pdgCodeParticle== abs(2212)) fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 1.5, fWeightJetJetMC); // p
				else if (pdgCodeParticle== abs(321)) fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 2.5, fWeightJetJetMC); // K+-
				else if (pdgCodeParticle== abs(2112)) fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 3.5, fWeightJetJetMC); // n
				else if (pdgCodeParticle== abs(310)) fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 4.5, fWeightJetJetMC); // K0s
				else if (pdgCodeParticle== abs(3122)) fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 5.5, fWeightJetJetMC); // Lambda
				else if (pdgCodeParticle== abs(13)) fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 6.5, fWeightJetJetMC); // mu+/-
				else if (pdgCodeParticle== abs(130)) fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 7.5, fWeightJetJetMC); // K0l
				else fHistoTrueClusBGPtvsSource[fiCut]->Fill(mesoncand->Pt(), 8.5, fWeightJetJetMC); // Rest
				
			}	
		}	
			
		if (TrueClusterCandidate->IsLargestComponentPhoton() && !TrueClusterCandidate->IsMerged() && !filled){
			fHistoTrueClusGammaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, fWeightJetJetMC);
			if (mesonIsSelected) fHistoTrueClusGammaInvMassvsPt[fiCut]->Fill(mesoncand->M(),mesoncand->Pt(), fWeightJetJetMC);
		}
	}
	return;
}



//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::ProcessMCParticles()
{
	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();
	
	// Loop over all primary MC particle
	for(UInt_t i = 0; i < fMCStack->GetNtrack(); i++) {
		if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){ 

			TParticle* particle = (TParticle *)fMCStack->Particle(i);
			if (!particle) continue;
			
			Int_t isMCFromMBHeader = -1;
			if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
				isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent);
				if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
			}
			
			if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
				->MesonIsSelectedMC(particle,fMCStack,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
				TParticle* daughter0 = (TParticle*)fMCStack->Particle(particle->GetFirstDaughter());
				TParticle* daughter1 = (TParticle*)fMCStack->Particle(particle->GetLastDaughter());
				
				Float_t weighted= 1;
				if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent)){
					if (particle->Pt()>0.005){
						weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, fMCStack, fInputEvent);
					}
				}
				Double_t mesonY = 10.;
				if(particle->Energy() - particle->Pz() == 0 || particle->Energy() + particle->Pz() == 0){
					mesonY=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
				} else{
					mesonY = 0.5*(TMath::Log((particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
				}
	
				Double_t alpha = -1;
				if (particle->GetPdgCode() == 111 || particle->GetPdgCode() == 221){
					alpha = TMath::Abs((daughter0->Energy() - daughter1->Energy()))/(daughter0->Energy() + daughter1->Energy());
				}
	
				if(particle->GetPdgCode() == 111){
					fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // All MC Pi0
					fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
					if (fIsMC==2)fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
					if (fDoMesonQA > 0 && fDoMesonQA < 3){
						fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted* fWeightJetJetMC); 
						fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha, fWeightJetJetMC); 
						if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
					}
				} else if(particle->GetPdgCode() == 221){
					fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // All MC Eta
					fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt(), fWeightJetJetMC);
					if (fIsMC==2)fHistoMCEtaWOEvtWeightPt[fiCut]->Fill(particle->Pt());
					if (fDoMesonQA > 0 && fDoMesonQA < 3){
						fHistoMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted* fWeightJetJetMC); 
						fHistoMCEtaPtAlpha[fiCut]->Fill(particle->Pt(),alpha, fWeightJetJetMC); 
						if (fIsMC == 2) fHistoMCEtaPtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),fWeightJetJetMC);
					}
				}
			
				// Check the acceptance for both gammas & whether they are counted as primaries as well
				Bool_t kDaughter0IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, particle->GetFirstDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
				Bool_t kDaughter1IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, particle->GetLastDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
				
				if( kDaughter0IsPrim && kDaughter1IsPrim &&
					((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCStack) &&
					((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCStack) ){					
					if(particle->GetPdgCode() == 111){
						fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // MC Pi0 with gamma in acc
					} else if(particle->GetPdgCode() == 221){
						fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted* fWeightJetJetMC); // MC Eta with gamma in acc
					}
				}
			}
		}
	}
  

}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::CalculatePi0Candidate(AliAODConversionPhoton* photon1, AliAODConversionPhoton* photon2){
	
	AliAODConversionMother *pi0cand = new AliAODConversionMother(photon1,photon2);
					
	if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
		fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(), fWeightJetJetMC);
		
		if (fDoMesonQA > 0 && fDoMesonQA < 3){
			if ( pi0cand->M() > 0.0 && pi0cand->M() < 0.4){
				fHistoMotherPi0PtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
				fHistoMotherPi0PtAlpha[fiCut]->Fill(pi0cand->Pt(),abs(pi0cand->GetAlpha()), fWeightJetJetMC);
				fHistoMotherPi0PtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(), fWeightJetJetMC);
			}
			if ( pi0cand->M() > 0.4 && pi0cand->M() < 0.8){
				fHistoMotherEtaPtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
				fHistoMotherEtaPtAlpha[fiCut]->Fill(pi0cand->Pt(),abs(pi0cand->GetAlpha()), fWeightJetJetMC);
				fHistoMotherEtaPtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(),  fWeightJetJetMC);
			}
		}
	}
	
	delete pi0cand;
	pi0cand=0x0;
}



//_________________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::SetLogBinningXTH2(TH2* histoRebin){
	TAxis *axisafter = histoRebin->GetXaxis();
	Int_t bins = axisafter->GetNbins();
	Double_t from = axisafter->GetXmin();
	Double_t to = axisafter->GetXmax();
	Double_t *newbins = new Double_t[bins+1];
	newbins[0] = from;
	Double_t factor = TMath::Power(to/from, 1./bins);
	for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
	axisafter->Set(bins, newbins);
	delete [] newbins;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskGammaCaloMerged::CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked)
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

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::Terminate(const Option_t *)
{
  
  //fOutputContainer->Print(); // Will crash on GRID
}


//_________________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::FillMultipleCountMap(map<Int_t,Int_t> &ma, Int_t tobechecked){
	if( ma.find(tobechecked) != ma.end() ) ma[tobechecked] += 1;
	else ma[tobechecked] = 2;
	return;
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::FillMultipleCountHistoAndClear(map<Int_t,Int_t> &ma, TH1F* hist){
	map<Int_t, Int_t>::iterator it;
	for (it = ma.begin(); it != ma.end(); it++){
		hist->Fill(it->second, fWeightJetJetMC);
	}
	ma.clear();
	return;
}
