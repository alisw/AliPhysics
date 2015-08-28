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
	fClusterCandidates(NULL),
	fEventCutArray(NULL),
	fEventCuts(NULL),
	fClusterCutArray(NULL),
	fCaloPhotonCuts(NULL),
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
	fHistoTruePi0InvMassPt(NULL),
	fHistoTrueEtaInvMassPt(NULL),
	fHistoTruePrimaryPi0InvMassPt(NULL),
	fHistoTruePrimaryEtaInvMassPt(NULL),
	fHistoTruePrimaryPi0W0WeightingInvMassPt(NULL),
	fHistoTruePrimaryEtaW0WeightingInvMassPt(NULL),
	fProfileTruePrimaryPi0WeightsInvMassPt(NULL),
	fProfileTruePrimaryEtaWeightsInvMassPt(NULL),
	fHistoTruePrimaryPi0MCPtResolPt(NULL),
	fHistoTruePrimaryEtaMCPtResolPt(NULL),
	fHistoTrueSecondaryPi0InvMassPt(NULL),
	fHistoTrueSecondaryPi0FromK0sInvMassPt(NULL),
	fHistoTrueK0sWithPi0DaughterMCPt(NULL),
	fHistoTrueSecondaryPi0FromLambdaInvMassPt(NULL),
	fHistoTrueLambdaWithPi0DaughterMCPt(NULL),
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
	fHistoNGoodESDTracksVsNClusterCandidates(NULL),
	fHistoSPDClusterTrackletBackground(NULL),
	fHistoNV0Tracks(NULL),
	fProfileEtaShift(NULL),
	fRandom(0),
	fnCuts(0),
	fiCut(0),
	fIsHeavyIon(0),
	fDoMesonAnalysis(kTRUE),
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
	fClusterCandidates(NULL),
	fEventCutArray(NULL),
	fEventCuts(NULL),
	fClusterCutArray(NULL),
	fCaloPhotonCuts(NULL),
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
	fHistoTruePi0InvMassPt(NULL),
	fHistoTrueEtaInvMassPt(NULL),
	fHistoTruePrimaryPi0InvMassPt(NULL),
	fHistoTruePrimaryEtaInvMassPt(NULL),
	fHistoTruePrimaryPi0W0WeightingInvMassPt(NULL),
	fHistoTruePrimaryEtaW0WeightingInvMassPt(NULL),
	fProfileTruePrimaryPi0WeightsInvMassPt(NULL),
	fProfileTruePrimaryEtaWeightsInvMassPt(NULL),
	fHistoTruePrimaryPi0MCPtResolPt(NULL),
	fHistoTruePrimaryEtaMCPtResolPt(NULL),
	fHistoTrueSecondaryPi0InvMassPt(NULL),
	fHistoTrueSecondaryPi0FromK0sInvMassPt(NULL),
	fHistoTrueK0sWithPi0DaughterMCPt(NULL),
	fHistoTrueSecondaryPi0FromLambdaInvMassPt(NULL),
	fHistoTrueLambdaWithPi0DaughterMCPt(NULL),
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
	fHistoNGoodESDTracksVsNClusterCandidates(NULL),
	fHistoSPDClusterTrackletBackground(NULL),
	fHistoNV0Tracks(NULL),
	fProfileEtaShift(NULL),
	fRandom(0),
	fnCuts(0),
	fiCut(0),
	fIsHeavyIon(0),
	fDoMesonAnalysis(kTRUE),
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
	if(fClusterCandidates){
		delete fClusterCandidates;
		fClusterCandidates = 0x0;
	}
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
  
	// Array of current cut's gammas
	fClusterCandidates = new TList();
	fClusterCandidates->SetOwner(kTRUE);
  
	fCutFolder = new TList*[fnCuts];
	fESDList = new TList*[fnCuts];

	fHistoNEvents = new TH1F*[fnCuts];
	if(fIsMC == 2){
		fHistoNEventsWOWeight = new TH1F*[fnCuts];
	}
	
	fHistoNGoodESDTracks = new TH1F*[fnCuts];
	fHistoVertexZ = new TH1F*[fnCuts];
	fHistoNClusterCandidates = new TH1F*[fnCuts];
	fHistoNGoodESDTracksVsNClusterCandidates = new TH2F*[fnCuts];
	fHistoSPDClusterTrackletBackground = new TH2F*[fnCuts];
	fHistoNV0Tracks = new TH1F*[fnCuts];
	fProfileEtaShift = new TProfile*[fnCuts];
  	
	if(fDoMesonAnalysis){
		fHistoMotherInvMassPt = new TH2F*[fnCuts];
		if (fDoMesonQA > 0 && fDoMesonQA < 3){
			fHistoMotherPi0PtY =  new TH2F*[fnCuts];
			fHistoMotherEtaPtY =  new TH2F*[fnCuts];
			fHistoMotherPi0PtAlpha =  new TH2F*[fnCuts];
			fHistoMotherEtaPtAlpha =  new TH2F*[fnCuts];
			fHistoMotherPi0PtOpenAngle =  new TH2F*[fnCuts];
			fHistoMotherEtaPtOpenAngle =  new TH2F*[fnCuts];
		}
	}
		
	fHistoClusGammaPt = new TH1F*[fnCuts];
	fHistoClusOverlapHeadersGammaPt = new TH1F*[fnCuts];
	
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		TString cutstringEvent 	= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
		TString cutstringCalo 	= ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
		TString cutstringMeson 	= "NoMesonCut";
		if(fDoMesonAnalysis)cutstringMeson = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
    
		fCutFolder[iCut] = new TList();
		fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s",cutstringEvent.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
		fCutFolder[iCut]->SetOwner(kTRUE);
		fOutputContainer->Add(fCutFolder[iCut]);
		fESDList[iCut] = new TList();
		fESDList[iCut]->SetName(Form("%s_%s_%s ESD histograms",cutstringEvent.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
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

   		fHistoClusGammaPt[iCut] = new TH1F("ClusGamma_Pt","ClusGamma_Pt",350,0,35);
		fESDList[iCut]->Add(fHistoClusGammaPt[iCut]);
		fHistoClusOverlapHeadersGammaPt[iCut] = new TH1F("ClusGammaOverlapHeaders_Pt","ClusGammaOverlapHeaders_Pt",350,0,35);
		fESDList[iCut]->Add(fHistoClusOverlapHeadersGammaPt[iCut]);
		
		if (fIsMC == 2){
			fHistoClusGammaPt[iCut]->Sumw2();
			fHistoClusOverlapHeadersGammaPt[iCut]->Sumw2();
		}
		
		if(fDoMesonAnalysis){
			fHistoMotherInvMassPt[iCut] = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",800,0,0.8,350,0,35);
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
	}
  
	if(fIsMC> 0){
		// MC Histogramms
		fMCList 	= new TList*[fnCuts];
		// True Histogramms
		fTrueList 	= new TList*[fnCuts];
		// Selected Header List

		fHeaderNameList 					= new TList*[fnCuts];
		fHistoMCHeaders 					= new TH1I*[fnCuts];
	    
		if(fDoMesonAnalysis){
			fHistoMCPi0Pt 					= new TH1F*[fnCuts];
			fHistoMCPi0WOWeightPt 			= new TH1F*[fnCuts];
			fHistoMCEtaPt 					= new TH1F*[fnCuts];
			fHistoMCEtaWOWeightPt 			= new TH1F*[fnCuts];
			fHistoMCPi0InAccPt 				= new TH1F*[fnCuts];
			fHistoMCEtaInAccPt 				= new TH1F*[fnCuts];
      
			if (fIsMC == 2){
				fHistoMCPi0WOEvtWeightPt 			= new TH1F*[fnCuts];
				fHistoMCEtaWOEvtWeightPt 			= new TH1F*[fnCuts];
			}	
			
			fHistoTruePi0InvMassPt 					= new TH2F*[fnCuts];
			fHistoTrueEtaInvMassPt 					= new TH2F*[fnCuts];
			fHistoDoubleCountTruePi0InvMassPt			= new TH2F*[fnCuts];
			fHistoDoubleCountTrueEtaInvMassPt			= new TH2F*[fnCuts];
			fHistoTruePrimaryPi0InvMassPt 			= new TH2F*[fnCuts];
			fHistoTruePrimaryEtaInvMassPt 			= new TH2F*[fnCuts];
			fHistoTruePrimaryPi0W0WeightingInvMassPt = new TH2F*[fnCuts];
			fHistoTruePrimaryEtaW0WeightingInvMassPt = new TH2F*[fnCuts];
			fProfileTruePrimaryPi0WeightsInvMassPt 	= new TProfile2D*[fnCuts];
			fProfileTruePrimaryEtaWeightsInvMassPt 	= new TProfile2D*[fnCuts];
			fHistoTrueSecondaryPi0InvMassPt 			= new TH2F*[fnCuts];
			fHistoTrueSecondaryPi0FromK0sInvMassPt 	= new TH2F*[fnCuts];
			fHistoTrueSecondaryPi0FromLambdaInvMassPt = new TH2F*[fnCuts];

			if (fDoMesonQA > 0 && fDoMesonQA < 3 ){
				fHistoMCPi0PtY 								= new TH2F*[fnCuts];
				fHistoMCEtaPtY 								= new TH2F*[fnCuts];
				fHistoMCPi0PtAlpha 							= new TH2F*[fnCuts];
				fHistoMCEtaPtAlpha 							= new TH2F*[fnCuts];
				if (fIsMC == 2){
					fHistoMCPi0PtJetPt 						= new TH2F*[fnCuts];
					fHistoMCEtaPtJetPt 						= new TH2F*[fnCuts];				
				}	

				if (fIsMC != 2){
					fHistoTruePrimaryPi0MCPtResolPt 			= new TH2F*[fnCuts];
					fHistoTruePrimaryEtaMCPtResolPt 			= new TH2F*[fnCuts];
					fHistoTrueK0sWithPi0DaughterMCPt 			= new TH1F*[fnCuts];
					fHistoTrueLambdaWithPi0DaughterMCPt 		= new TH1F*[fnCuts];
				}	
				fHistoTruePi0PtY 							= new TH2F*[fnCuts];
				fHistoTrueEtaPtY 							= new TH2F*[fnCuts];
				fHistoTruePi0PtAlpha 						= new TH2F*[fnCuts];
				fHistoTrueEtaPtAlpha 						= new TH2F*[fnCuts];
				fHistoTruePi0PtOpenAngle 					= new TH2F*[fnCuts];
				fHistoTrueEtaPtOpenAngle 					= new TH2F*[fnCuts];
			}
		}
    
    
    
		for(Int_t iCut = 0; iCut<fnCuts;iCut++){
			TString cutstringEvent 	= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
			TString cutstringCalo 	= ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
			TString cutstringMeson 	= "NoMesonCut";
			if(fDoMesonAnalysis)cutstringMeson = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
      
			fMCList[iCut] = new TList();
			fMCList[iCut]->SetName(Form("%s_%s_%s MC histograms",cutstringEvent.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
			fMCList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fMCList[iCut]);

			fHistoMCHeaders[iCut] = new TH1I("MC_Headers","MC_Headers",20,0,20);
			fMCList[iCut]->Add(fHistoMCHeaders[iCut]);
			
			if(fDoMesonAnalysis){
				fHistoMCPi0Pt[iCut] = new TH1F("MC_Pi0_Pt","MC_Pi0_Pt",350,0,35);
				fHistoMCPi0Pt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCPi0Pt[iCut]);
				fHistoMCPi0WOWeightPt[iCut] = new TH1F("MC_Pi0_WOWeights_Pt","MC_Pi0_WOWeights_Pt",350,0,35);
				fHistoMCPi0WOWeightPt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCPi0WOWeightPt[iCut]);
				
				fHistoMCEtaPt[iCut] = new TH1F("MC_Eta_Pt","MC_Eta_Pt",350,0,35);
				fHistoMCEtaPt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCEtaPt[iCut]);
				fHistoMCEtaWOWeightPt[iCut] = new TH1F("MC_Eta_WOWeights_Pt","MC_Eta_WOWeights_Pt",350,0,35);
				fHistoMCEtaWOWeightPt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCEtaWOWeightPt[iCut]);
				fHistoMCPi0InAccPt[iCut] = new TH1F("MC_Pi0InAcc_Pt","MC_Pi0InAcc_Pt",350,0,35);
				fHistoMCPi0InAccPt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCPi0InAccPt[iCut]);
				fHistoMCEtaInAccPt[iCut] = new TH1F("MC_EtaInAcc_Pt","MC_EtaInAcc_Pt",350,0,35);
				fHistoMCEtaInAccPt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCEtaInAccPt[iCut]);
				if (fIsMC == 2){
					fHistoMCPi0WOWeightPt[iCut]->Sumw2();
					fHistoMCEtaWOWeightPt[iCut]->Sumw2();
					fHistoMCPi0WOEvtWeightPt[iCut] = new TH1F("MC_Pi0_WOEventWeights_Pt","MC_Pi0_WOEventWeights_Pt",350,0,35);
					fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightPt[iCut]);
					fHistoMCEtaWOEvtWeightPt[iCut] = new TH1F("MC_Eta_WOEventWeights_Pt","MC_Eta_WOEventWeights_Pt",350,0,35);
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

					fHistoMCPi0PtAlpha[iCut]->Sumw2();
					fHistoMCEtaPtAlpha[iCut]->Sumw2();
				}
        
			}
			fTrueList[iCut] = new TList();
			fTrueList[iCut]->SetName(Form("%s_%s_%s True histograms",cutstringEvent.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
			fTrueList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fTrueList[iCut]);
                  
			if(fDoMesonAnalysis){
				fHistoTruePi0InvMassPt[iCut] = new TH2F("ESD_TruePi0_InvMass_Pt","ESD_TruePi0_InvMass_Pt",800,0,0.8,350,0,35);
				fTrueList[iCut]->Add(fHistoTruePi0InvMassPt[iCut]);
				fHistoTrueEtaInvMassPt[iCut] = new TH2F("ESD_TrueEta_InvMass_Pt","ESD_TrueEta_InvMass_Pt",800,0,0.8,350,0,35);
				fTrueList[iCut]->Add(fHistoTrueEtaInvMassPt[iCut]);

				fHistoDoubleCountTruePi0InvMassPt[iCut] = new TH2F("ESD_TrueDoubleCountPi0_InvMass_Pt","ESD_TrueDoubleCountPi0_InvMass_Pt",800,0,0.8,350,0,35);
				fTrueList[iCut]->Add(fHistoDoubleCountTruePi0InvMassPt[iCut]);
				fHistoDoubleCountTrueEtaInvMassPt[iCut] = new TH2F("ESD_TrueDoubleCountEta_InvMass_Pt","ESD_TrueDoubleCountEta_InvMass_Pt",800,0,0.8,350,0,35);
				fTrueList[iCut]->Add(fHistoDoubleCountTrueEtaInvMassPt[iCut]);

				fHistoTruePrimaryPi0InvMassPt[iCut] = new TH2F("ESD_TruePrimaryPi0_InvMass_Pt", "ESD_TruePrimaryPi0_InvMass_Pt", 800,0,0.8,350,0,35);
				fHistoTruePrimaryPi0InvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTruePrimaryPi0InvMassPt[iCut]);
				fHistoTruePrimaryEtaInvMassPt[iCut] = new TH2F("ESD_TruePrimaryEta_InvMass_Pt", "ESD_TruePrimaryEta_InvMass_Pt", 800,0,0.8,350,0,35);
				fHistoTruePrimaryEtaInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTruePrimaryEtaInvMassPt[iCut]);

				fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut] = new TH2F("ESD_TruePrimaryPi0W0Weights_InvMass_Pt", "ESD_TruePrimaryPi0W0Weights_InvMass_Pt", 800,0,0.8,350,0,35);
				fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]);
				fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut] = new TH2F("ESD_TruePrimaryEtaW0Weights_InvMass_Pt", "ESD_TruePrimaryEtaW0Weights_InvMass_Pt", 800,0,0.8,350,0,35);
				fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]);

				fProfileTruePrimaryPi0WeightsInvMassPt[iCut] = new TProfile2D("ESD_TruePrimaryPi0Weights_InvMass_Pt", "ESD_TruePrimaryPi0Weights_InvMass_Pt", 800,0,0.8,350,0,35);
				fProfileTruePrimaryPi0WeightsInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fProfileTruePrimaryPi0WeightsInvMassPt[iCut]);
				fProfileTruePrimaryEtaWeightsInvMassPt[iCut] = new TProfile2D("ESD_TruePrimaryEtaWeights_InvMass_Pt", "ESD_TruePrimaryEtaWeights_InvMass_Pt", 800,0,0.8,350,0,35);
				fProfileTruePrimaryEtaWeightsInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fProfileTruePrimaryEtaWeightsInvMassPt[iCut]);

				fHistoTrueSecondaryPi0InvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0_InvMass_Pt", "ESD_TrueSecondaryPi0_InvMass_Pt", 800,0,0.8,350,0,35);
				fHistoTrueSecondaryPi0InvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTrueSecondaryPi0InvMassPt[iCut]);
				
				fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromK0s_InvMass_Pt","ESD_TrueSecondaryPi0FromK0s_InvMass_Pt",800,0,0.8,350,0,35);
				fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]);
				fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromLambda_InvMass_Pt","ESD_TrueSecondaryPi0FromLambda_InvMass_Pt",800,0,0.8,350,0,35);
				fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]);

				if (fIsMC == 2){
					fHistoTruePi0InvMassPt[iCut]->Sumw2();
					fHistoTrueEtaInvMassPt[iCut]->Sumw2();
					fHistoDoubleCountTruePi0InvMassPt[iCut]->Sumw2();
					fHistoDoubleCountTrueEtaInvMassPt[iCut]->Sumw2();
					fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->Sumw2();
					fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->Sumw2();
					fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]->Sumw2();
				}


				if (fDoMesonQA > 0 && fDoMesonQA < 3){
					if (fIsMC != 2){						
						fHistoTruePrimaryPi0MCPtResolPt[iCut] = new TH2F("ESD_TruePrimaryPi0_MCPt_ResolPt","ESD_TruePrimaryPi0_ResolPt_MCPt",500,0.03,35,1000,-1.,1.);
						fHistoTruePrimaryPi0MCPtResolPt[iCut]->Sumw2();
						SetLogBinningXTH2(fHistoTruePrimaryPi0MCPtResolPt[iCut]);
						fTrueList[iCut]->Add(fHistoTruePrimaryPi0MCPtResolPt[iCut]);
						fHistoTruePrimaryEtaMCPtResolPt[iCut]  = new TH2F("ESD_TruePrimaryEta_MCPt_ResolPt","ESD_TruePrimaryEta_ResolPt_MCPt",500,0.03,35,1000,-1.,1.);
						fHistoTruePrimaryEtaMCPtResolPt[iCut]->Sumw2();
						SetLogBinningXTH2(fHistoTruePrimaryEtaMCPtResolPt[iCut]);
						fTrueList[iCut]->Add(fHistoTruePrimaryEtaMCPtResolPt[iCut]);
						fHistoTrueK0sWithPi0DaughterMCPt[iCut] = new TH1F("ESD_TrueK0sWithPi0Daughter_MCPt","ESD_TrueK0sWithPi0Daughter_MCPt",350,0,35);
						fTrueList[iCut]->Add(fHistoTrueK0sWithPi0DaughterMCPt[iCut]);
						fHistoTrueLambdaWithPi0DaughterMCPt[iCut] = new TH1F("ESD_TrueLambdaWithPi0Daughter_MCPt","ESD_TrueLambdaWithPi0Daughter_MCPt",350,0,35);
						fTrueList[iCut]->Add(fHistoTrueLambdaWithPi0DaughterMCPt[iCut]);
					}
					
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
					
					if (fIsMC==2){
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
	}  
    
	fVectorDoubleCountTruePi0s.clear();
	fVectorDoubleCountTrueEtas.clear();

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
		if(fSetPlotHistsExtQA){
			if(!((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))) continue;
			if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetExtQAHistograms()){
				fCutFolder[iCut]->Add(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetExtQAHistograms());
			}
		}
		if(fDoMesonAnalysis){
			if(!((AliConversionMesonCuts*)fMesonCutArray->At(iCut))) continue;
			if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms()){
				fCutFolder[iCut]->Add(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
			}
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
			if (fIsMC==2) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
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
			fHistoNEvents[iCut]->Fill(10,fWeightJetJetMC);
			if (fIsMC==2) fHistoNEventsWOWeight[iCut]->Fill(10);
			continue;
		}

		Bool_t triggered = 1;
		if(eventNotAccepted){
		// cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
			fHistoNEvents[iCut]->Fill(eventNotAccepted, fWeightJetJetMC); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
			if (fIsMC==2) fHistoNEventsWOWeight[iCut]->Fill(eventNotAccepted);
			if (eventNotAccepted==3 && fIsMC>0){
				triggered = 0;
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
		if (triggered == 1) {
			fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC); // Should be 0 here
			if (fIsMC==2) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here

			fHistoNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(), fWeightJetJetMC);
			fHistoVertexZ[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetZ(), fWeightJetJetMC);
			fHistoSPDClusterTrackletBackground[iCut]->Fill(fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),(fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)), fWeightJetJetMC);
			if(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsHeavyIon() == 2)	fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A(), fWeightJetJetMC);
				else fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C(), fWeightJetJetMC);
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
		
		if (triggered==0) continue;
		
		// it is in the loop to have the same conversion cut string (used also for MC stuff that should be same for V0 and Cluster)
		ProcessClusters();  					// process calo clusters

		fHistoNClusterCandidates[iCut]->Fill(fClusterCandidates->GetEntries(), fWeightJetJetMC);
		fHistoNGoodESDTracksVsNClusterCandidates[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),fClusterCandidates->GetEntries(), fWeightJetJetMC);
		if(fDoMesonAnalysis){ // Meson Analysis
			
			CalculatePi0Candidates(); // Combine Gammas from conversion and from calo
	
			fVectorDoubleCountTruePi0s.clear();
			fVectorDoubleCountTrueEtas.clear();
		}

		fClusterCandidates->Clear(); // delete cluster candidates
	}
	
	PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::ProcessClusters(){
	
	Int_t nclus = 0;
	nclus = fInputEvent->GetNumberOfCaloClusters();
	
// 	cout << nclus << endl;
	
	if(nclus == 0)	return;
	
	// plotting histograms on cell/tower level, only if extendedMatchAndQA > 1
	((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillHistogramsExtendedQA(fInputEvent);

	// vertex
	Double_t vertex[3] = {0};
	InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
	
	// Loop over EMCal clusters
	for(Long_t i = 0; i < nclus; i++){
		
		AliVCluster* clus = NULL;
		if(fInputEvent->IsA()==AliESDEvent::Class()) clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
		else if(fInputEvent->IsA()==AliAODEvent::Class()) clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));

		if (!clus){
			delete clus;
			continue;
		}	
		if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fIsMC)) continue;
		// TLorentzvector with cluster
		
		const Int_t   nc = clus->GetNCells();		
		Int_t   absCellIdList[nc]; 
		Float_t maxEList[nc]; 

		Int_t nlm = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetNumberOfLocalMaxima(clus, fInputEvent, absCellIdList, maxEList);
		if (nlm == 2){
			AliAODCaloCluster* clusSub1 = new AliAODCaloCluster();
			AliAODCaloCluster* clusSub2 = new AliAODCaloCluster();
			
			((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->SplitEnergy(absCellIdList[0], absCellIdList[1],
																		   clus, fInputEvent, fIsMC, clusSub1, clusSub2);
			cout << clusSub1->E() << "\t" <<clusSub2->E() << endl; 
			
			delete clusSub1;
			delete clusSub2;
		}	
			
		
		TLorentzVector clusterVector;
		clus->GetMomentum(clusterVector,vertex);
		
		
		TLorentzVector* tmpvec = new TLorentzVector();
		tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());
	
		
		
		// convert to AODConversionPhoton
		AliAODConversionPhoton *PhotonCandidate=new AliAODConversionPhoton(tmpvec);
		if(!PhotonCandidate){
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

// 		if(fIsMC> 0){
// 			ProcessTrueClusterCandidates(PhotonCandidate);
// 		}

		if (fIsFromMBHeader && fIsOverlappingWithOtherHeader) fHistoClusOverlapHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), fWeightJetJetMC);
		if (fIsFromMBHeader && !fIsOverlappingWithOtherHeader){
			fHistoClusGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), fWeightJetJetMC);	
			fClusterCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas
		} else{
			delete PhotonCandidate;
		}
		
		delete clus;
		delete tmpvec;
	}
	
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::ProcessTrueClusterCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
		
	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();

	TParticle *Photon = NULL;
	if (!TruePhotonCandidate->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set task will abort");
	
	if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0) Photon = fMCStack->Particle(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
		else return;
		
	if(Photon == NULL){
	//    cout << "no photon" << endl;
		return;
	}

	Int_t pdgCodeParticle = Photon->GetPdgCode();
	TruePhotonCandidate->SetCaloPhotonMCFlags(fMCStack);
	
	// True Photon
	if(fIsFromMBHeader && !fIsOverlappingWithOtherHeader){
		
		Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, TruePhotonCandidate->GetCaloPhotonMCLabel(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
		
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
			
			if(fDoMesonAnalysis){
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
  

}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::CalculatePi0Candidates(){
	
	// Conversion Gammas
	if(fClusterCandidates->GetEntries()>0){
		
		for(Int_t firstGammaIndex=0;firstGammaIndex<fClusterCandidates->GetEntries();firstGammaIndex++){
			AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(firstGammaIndex));
			if (gamma0==NULL) continue;
			
			for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
				AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
				if (gamma1==NULL) continue;
				
				AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);
				pi0cand->SetLabels(firstGammaIndex,secondGammaIndex);
								
				if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
					fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt(), fWeightJetJetMC);
					
					if (fDoMesonQA > 0 && fDoMesonQA < 3){
						if ( pi0cand->M() > 0.05 && pi0cand->M() < 0.17){
							fHistoMotherPi0PtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
							fHistoMotherPi0PtAlpha[fiCut]->Fill(pi0cand->Pt(),abs(pi0cand->GetAlpha()), fWeightJetJetMC);
							fHistoMotherPi0PtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(), fWeightJetJetMC);
						}
						if ( pi0cand->M() > 0.45 && pi0cand->M() < 0.65){
							fHistoMotherEtaPtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
							fHistoMotherEtaPtAlpha[fiCut]->Fill(pi0cand->Pt(),abs(pi0cand->GetAlpha()), fWeightJetJetMC);
							fHistoMotherEtaPtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle(),  fWeightJetJetMC);
						}
					}
				
// 					if(fIsMC> 0){
// 						ProcessTrueMesonCandidates(pi0cand,gamma0,gamma1);
// 					}

				}
				
				delete pi0cand;
				pi0cand=0x0;
			}
		}
	}
}
//______________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
	// Process True Mesons
	AliStack *MCStack = fMCEvent->Stack();
	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();
	
	Bool_t isTruePi0 				= kFALSE;
	Bool_t isTrueEta 				= kFALSE;
	Bool_t isTrueGamma				= kFALSE;
	Bool_t isSameConvertedGamma 	= kFALSE;
	Int_t convertedPhotonLabel0		= -1;
	Int_t convertedPhotonLabel1		= -1;
	
	Int_t gamma0MCLabel = TrueGammaCandidate0->GetCaloPhotonMCLabel(0); 	// get most probable MC label
	Int_t gamma0MotherLabel = -1;

	TParticle * gammaMC0 = 0x0;
	if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
		gammaMC0 = (TParticle*)MCStack->Particle(gamma0MCLabel);
		if (TrueGammaCandidate0->IsLargestComponentPhoton() || TrueGammaCandidate0->IsLargestComponentElectron()){		// largest component is electro magnetic
			// get mother of interest (pi0 or eta)
			if (TrueGammaCandidate0->IsLargestComponentPhoton()){														// for photons its the direct mother 
				gamma0MotherLabel=gammaMC0->GetMother(0);
			} else if (TrueGammaCandidate0->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
				if (TrueGammaCandidate0->IsConversion()){
					convertedPhotonLabel0 = gammaMC0->GetMother(0);
					gamma0MotherLabel=MCStack->Particle(gammaMC0->GetMother(0))->GetMother(0);
				} else {
					gamma0MotherLabel=gammaMC0->GetMother(0); 
				}	
			}
		} 	
	}
	if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
	
	Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); 	// get most probable MC label
	Int_t gamma1MotherLabel = -1;
	// check if 

	TParticle * gammaMC1 = 0x0;
	if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
		// Daughters Gamma 1
		gammaMC1 = (TParticle*)MCStack->Particle(gamma1MCLabel);
		if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){		// largest component is electro magnetic
			// get mother of interest (pi0 or eta)
			if (TrueGammaCandidate1->IsLargestComponentPhoton()){														// for photons its the direct mother 
				gamma1MotherLabel			= gammaMC1->GetMother(0);
			} else if (TrueGammaCandidate1->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
				if (TrueGammaCandidate1->IsConversion()){
					convertedPhotonLabel1 	= gammaMC1->GetMother(0);
					gamma1MotherLabel		= MCStack->Particle(gammaMC1->GetMother(0))->GetMother(0);
				} else {
					gamma1MotherLabel		= gammaMC1->GetMother(0);
				}
			}
		}
	}
	
	if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
		if(((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
			isTruePi0=kTRUE;
		}
		if(((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetPdgCode() == 221){
			isTrueEta=kTRUE;
		}
	}

	
	if (convertedPhotonLabel0 > -1 && convertedPhotonLabel1 > -1){
		if (convertedPhotonLabel0==convertedPhotonLabel1){
			isSameConvertedGamma = kTRUE;
		}	
	}
	
	if(isTruePi0 || isTrueEta){// True Pion or Eta
		if (isTruePi0){
			fHistoTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  fWeightJetJetMC);
		}	
		if (isTrueEta){
			fHistoTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  fWeightJetJetMC);
		}
		
		if (fDoMesonQA > 0 && fDoMesonQA < 3){
			if (isTruePi0){
				if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
					fHistoTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
					fHistoTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),abs(Pi0Candidate->GetAlpha()), fWeightJetJetMC);
					fHistoTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle(), fWeightJetJetMC);
				}
			} else if (isTrueEta){
				if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
					fHistoTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), fWeightJetJetMC);
					fHistoTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),abs(Pi0Candidate->GetAlpha()), fWeightJetJetMC);
					fHistoTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle(), fWeightJetJetMC);
				}
			}
		}
		Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, gamma0MotherLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
		if(!isPrimary){ // Secondary Meson
			// filling secondary histograms
			Int_t secMotherLabel = ((TParticle*)MCStack->Particle(gamma0MotherLabel))->GetMother(0);
			Float_t weightedSec= 1;
			if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, fMCStack, fInputEvent) && MCStack->Particle(secMotherLabel)->GetPdgCode()==310){
				weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),secMotherLabel, fMCStack, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
				//cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
			}
			if (isTruePi0) fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* fWeightJetJetMC);
				
			
		} else { // Only primary pi0 for efficiency calculation
			// filling primary histograms
			Float_t weighted= 1;
			if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, fMCStack, fInputEvent)){
				if (((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt()>0.005){
					weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),gamma1MotherLabel, fMCStack, fInputEvent);
					//                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
				}
			}
			if (isTruePi0){
				fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* fWeightJetJetMC);
				fHistoTruePrimaryPi0W0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
				fProfileTruePrimaryPi0WeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* fWeightJetJetMC);
				if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), weighted*fWeightJetJetMC);
			} else if (isTrueEta){
				fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* fWeightJetJetMC);
				fHistoTruePrimaryEtaW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), fWeightJetJetMC);
				fProfileTruePrimaryEtaWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* fWeightJetJetMC);
				if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,gamma0MotherLabel)) fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), weighted*fWeightJetJetMC);
			}	
				
			if (fDoMesonQA > 0 && fDoMesonQA < 3 && fIsMC<2){
				if(isTruePi0){ // Only primary pi0 for resolution
					fHistoTruePrimaryPi0MCPtResolPt[fiCut]->Fill(((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt())/((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt(),weighted* fWeightJetJetMC);
				}
				if (isTrueEta){ // Only primary eta for resolution
					fHistoTruePrimaryEtaMCPtResolPt[fiCut]->Fill(((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt())/((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt(),weighted* fWeightJetJetMC);
				}
			}
		}	
	}
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
