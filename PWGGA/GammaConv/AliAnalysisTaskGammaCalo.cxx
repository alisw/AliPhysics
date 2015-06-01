/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *									 									  *
 * Author: Baldo Sahlmueller, Friederike Bock				       		  *
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
#include "AliAnalysisTaskGammaCalo.h"
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

ClassImp(AliAnalysisTaskGammaCalo)

//________________________________________________________________________
AliAnalysisTaskGammaCalo::AliAnalysisTaskGammaCalo(): AliAnalysisTaskSE(),
	fV0Reader(NULL),
	fBGHandler(NULL),
	fInputEvent(NULL),
	fMCEvent(NULL),
	fMCStack(NULL),
	fCutFolder(NULL),
	fESDList(NULL),
	fBackList(NULL),
	fMotherList(NULL),
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
	fHistoMotherInvMass3ClusterPt(NULL),
	fSparseMotherInvMassPtZM(NULL),
	fHistoMotherBackInvMassPt(NULL),
	fSparseMotherBackInvMassPtZM(NULL),
	fHistoMotherInvMassEalpha(NULL),
	fHistoMotherPi0PtY(NULL),
	fHistoMotherEtaPtY(NULL),
	fHistoMotherPi0PtAlpha(NULL),
	fHistoMotherEtaPtAlpha(NULL),
	fHistoMotherPi0PtOpenAngle(NULL),
	fHistoMotherEtaPtOpenAngle(NULL),
	fHistoMotherInvMassECalib(NULL),
	fHistoMotherInvMassECalibalpha(NULL),
	fHistoClusGammaPt(NULL),
	fHistoClusOverlapHeadersGammaPt(NULL),
	fHistoMCHeaders(NULL),
	fHistoMCAllGammaPt(NULL),
	fHistoMCDecayGammaPi0Pt(NULL),
	fHistoMCDecayGammaRhoPt(NULL),
	fHistoMCDecayGammaEtaPt(NULL),
	fHistoMCDecayGammaOmegaPt(NULL),
	fHistoMCDecayGammaEtapPt(NULL),
	fHistoMCDecayGammaPhiPt(NULL),
	fHistoMCDecayGammaSigmaPt(NULL),
	fHistoMCPi0Pt(NULL),
	fHistoMCPi0WOWeightPt(NULL),
	fHistoMCEtaPt(NULL),
	fHistoMCEtaWOWeightPt(NULL),
	fHistoMCPi0InAccPt(NULL),
	fHistoMCEtaInAccPt(NULL),
	fHistoMCPi0PtY(NULL),
	fHistoMCEtaPtY(NULL),
	fHistoMCPi0PtAlpha(NULL),
	fHistoMCEtaPtAlpha(NULL),
	fHistoMCK0sPt(NULL),
	fHistoMCK0sWOWeightPt(NULL),
	fHistoMCK0sPtY(NULL),
	fHistoMCSecPi0PtvsSource(NULL),
	fHistoMCSecPi0Source(NULL),
	fHistoMCSecEtaPt(NULL),
	fHistoMCSecEtaSource(NULL),
	fHistoTruePi0InvMassPt(NULL),
	fHistoTrueEtaInvMassPt(NULL),
	fHistoTruePi0CaloPhotonInvMassPt(NULL),
	fHistoTrueEtaCaloPhotonInvMassPt(NULL),
	fHistoTruePi0CaloConvertedPhotonInvMassPt(NULL),
	fHistoTrueEtaCaloConvertedPhotonInvMassPt(NULL),
	fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt(NULL),
	fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt(NULL),
	fHistoTruePi0CaloElectronInvMassPt(NULL),
	fHistoTrueEtaCaloElectronInvMassPt(NULL),
	fHistoTruePi0CaloMergedClusterInvMassPt(NULL),
	fHistoTrueEtaCaloMergedClusterInvMassPt(NULL),
	fHistoTruePi0CaloMergedClusterPartConvInvMassPt(NULL),
	fHistoTrueEtaCaloMergedClusterPartConvInvMassPt(NULL),
	fHistoTruePi0NonMergedElectronPhotonInvMassPt(NULL),
	fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt(NULL),
	fHistoTruePi0Category1(NULL),
	fHistoTrueEtaCategory1(NULL),
	fHistoTruePi0Category2(NULL),
	fHistoTrueEtaCategory2(NULL),
	fHistoTruePi0Category3(NULL),
	fHistoTrueEtaCategory3(NULL),
	fHistoTruePi0Category4_6(NULL),
	fHistoTrueEtaCategory4_6(NULL),
	fHistoTruePi0Category5(NULL),
	fHistoTrueEtaCategory5(NULL),
	fHistoTruePi0Category7(NULL),
	fHistoTrueEtaCategory7(NULL),
	fHistoTruePi0Category8(NULL),
	fHistoTrueEtaCategory8(NULL),
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
	fHistoTrueSecondaryPi0FromEtaInvMassPt(NULL),
	fHistoTrueEtaWithPi0DaughterMCPt(NULL),
	fHistoTrueSecondaryPi0FromLambdaInvMassPt(NULL),
	fHistoTrueLambdaWithPi0DaughterMCPt(NULL),
	fHistoTrueBckGGInvMassPt(NULL),
	fHistoTrueBckContInvMassPt(NULL),
	fHistoTruePi0PtY(NULL),
	fHistoTrueEtaPtY(NULL),
	fHistoTruePi0PtAlpha(NULL),
	fHistoTrueEtaPtAlpha(NULL),
	fHistoTruePi0PtOpenAngle(NULL),
	fHistoTrueEtaPtOpenAngle(NULL),
	fHistoClusPhotonBGPt(NULL),
	fHistoClusPhotonPlusConvBGPt(NULL),
	fHistoTrueClusGammaPt(NULL),
	fHistoTrueClusUnConvGammaPt(NULL),
	fHistoTrueClusUnConvGammaMCPt(NULL),
	fHistoTrueClusElectronPt(NULL),
	fHistoTrueClusConvGammaPt(NULL),
	fHistoTrueClusConvGammaMCPt(NULL),
	fHistoTrueClusConvGammaFullyPt(NULL),
	fHistoTrueClusMergedGammaPt(NULL),
	fHistoTrueClusMergedPartConvGammaPt(NULL),
	fHistoTrueClusDalitzPt(NULL),
	fHistoTrueClusDalitzMergedPt(NULL),
	fHistoTrueClusPhotonFromElecMotherPt(NULL),
	fHistoTrueClusShowerPt(NULL),
	fHistoTrueClusSubLeadingPt(NULL),
	fHistoTrueClusNParticles(NULL),
	fHistoTrueClusEMNonLeadingPt(NULL),
	fHistoTrueNLabelsInClus(NULL),
	fHistoTruePrimaryClusGammaPt(NULL),
	fHistoTruePrimaryClusGammaESDPtMCPt(NULL),
	fHistoTruePrimaryClusConvGammaPt(NULL),
	fHistoTruePrimaryClusConvGammaESDPtMCPt(NULL),
	fHistoTrueSecondaryClusGammaPt(NULL),
	fHistoTrueSecondaryClusConvGammaPt(NULL),
	fHistoTrueSecondaryClusGammaFromXFromK0sPt(NULL),
	fHistoTrueSecondaryClusConvGammaFromXFromK0sPt(NULL),
	fHistoTrueSecondaryClusGammaFromXFromLambdaPt(NULL),
	fHistoTrueSecondaryClusConvGammaFromXFromLambdaPt(NULL),
	fHistoTrueSecondaryClusGammaFromXFromEtasPt(NULL),
	fHistoTrueSecondaryClusConvGammaFromXFromEtasPt(NULL),
	fHistoDoubleCountTruePi0InvMassPt(NULL),
	fHistoDoubleCountTrueEtaInvMassPt(NULL),
	fVectorDoubleCountTruePi0s(0),
	fVectorDoubleCountTrueEtas(0),
	fHistoNEvents(NULL),
	fHistoNGoodESDTracks(NULL),
	fHistoVertexZ(NULL),
	fHistoNGammaCandidates(NULL),
	fHistoNGoodESDTracksVsNGammaCanditates(NULL),
	fHistoNV0Tracks(NULL),
	fProfileEtaShift(NULL),
	fEventPlaneAngle(-100),
	fRandom(0),
	fnCuts(0),
	fiCut(0),
	fIsHeavyIon(0),
	fDoMesonAnalysis(kTRUE),
	fDoMesonQA(0),
	fDoClusterQA(0),
	fIsFromMBHeader(kTRUE),
	fIsOverlappingWithOtherHeader(kFALSE),
	fIsMC(kFALSE),
	fDoTHnSparse(kTRUE)
{
  
}

//________________________________________________________________________
AliAnalysisTaskGammaCalo::AliAnalysisTaskGammaCalo(const char *name):
	AliAnalysisTaskSE(name),
	fV0Reader(NULL),
	fBGHandler(NULL),
	fInputEvent(NULL),
	fMCEvent(NULL),
	fMCStack(NULL),
	fCutFolder(NULL),
	fESDList(NULL),
	fBackList(NULL),
	fMotherList(NULL),
	fTrueList(NULL),
	fMCList(NULL),
	fHeaderNameList(NULL),
	fOutputContainer(0),
	fClusterCandidates(NULL),
	fEventCutArray(NULL),
	fEventCuts(NULL),
	fClusterCutArray(NULL),
	fCaloPhotonCuts(NULL),
	fMesonCutArray(NULL),
	fMesonCuts(NULL),
	fHistoMotherInvMassPt(NULL),
	fHistoMotherInvMass3ClusterPt(NULL),
	fSparseMotherInvMassPtZM(NULL),
	fHistoMotherBackInvMassPt(NULL),
	fSparseMotherBackInvMassPtZM(NULL),
	fHistoMotherInvMassEalpha(NULL),
	fHistoMotherPi0PtY(NULL),
	fHistoMotherEtaPtY(NULL),
	fHistoMotherPi0PtAlpha(NULL),
	fHistoMotherEtaPtAlpha(NULL),
	fHistoMotherPi0PtOpenAngle(NULL),
	fHistoMotherEtaPtOpenAngle(NULL),
	fHistoMotherInvMassECalib(NULL),
	fHistoMotherInvMassECalibalpha(NULL),
	fHistoClusGammaPt(NULL),
	fHistoClusOverlapHeadersGammaPt(NULL),
	fHistoMCHeaders(NULL),
	fHistoMCAllGammaPt(NULL),
	fHistoMCDecayGammaPi0Pt(NULL),
	fHistoMCDecayGammaRhoPt(NULL),
	fHistoMCDecayGammaEtaPt(NULL),
	fHistoMCDecayGammaOmegaPt(NULL),
	fHistoMCDecayGammaEtapPt(NULL),
	fHistoMCDecayGammaPhiPt(NULL),
	fHistoMCDecayGammaSigmaPt(NULL),
	fHistoMCPi0Pt(NULL),
	fHistoMCPi0WOWeightPt(NULL),
	fHistoMCEtaPt(NULL),
	fHistoMCEtaWOWeightPt(NULL),
	fHistoMCPi0InAccPt(NULL),
	fHistoMCEtaInAccPt(NULL),
	fHistoMCPi0PtY(NULL),
	fHistoMCEtaPtY(NULL),
	fHistoMCPi0PtAlpha(NULL),
	fHistoMCEtaPtAlpha(NULL),
	fHistoMCK0sPt(NULL),
	fHistoMCK0sWOWeightPt(NULL),
	fHistoMCK0sPtY(NULL),
	fHistoMCSecPi0PtvsSource(NULL),
	fHistoMCSecPi0Source(NULL),
	fHistoMCSecEtaPt(NULL),
	fHistoMCSecEtaSource(NULL),
	fHistoTruePi0InvMassPt(NULL),
	fHistoTrueEtaInvMassPt(NULL),
	fHistoTruePi0CaloPhotonInvMassPt(NULL),
	fHistoTrueEtaCaloPhotonInvMassPt(NULL),
	fHistoTruePi0CaloConvertedPhotonInvMassPt(NULL),
	fHistoTrueEtaCaloConvertedPhotonInvMassPt(NULL),
	fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt(NULL),
	fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt(NULL),
	fHistoTruePi0CaloElectronInvMassPt(NULL),
	fHistoTrueEtaCaloElectronInvMassPt(NULL),
	fHistoTruePi0CaloMergedClusterInvMassPt(NULL),
	fHistoTrueEtaCaloMergedClusterInvMassPt(NULL),
	fHistoTruePi0CaloMergedClusterPartConvInvMassPt(NULL),
	fHistoTrueEtaCaloMergedClusterPartConvInvMassPt(NULL),
	fHistoTruePi0NonMergedElectronPhotonInvMassPt(NULL),
	fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt(NULL),
	fHistoTruePi0Category1(NULL),
	fHistoTrueEtaCategory1(NULL),
	fHistoTruePi0Category2(NULL),
	fHistoTrueEtaCategory2(NULL),
	fHistoTruePi0Category3(NULL),
	fHistoTrueEtaCategory3(NULL),
	fHistoTruePi0Category4_6(NULL),
	fHistoTrueEtaCategory4_6(NULL),
	fHistoTruePi0Category5(NULL),
	fHistoTrueEtaCategory5(NULL),
	fHistoTruePi0Category7(NULL),
	fHistoTrueEtaCategory7(NULL),
	fHistoTruePi0Category8(NULL),
	fHistoTrueEtaCategory8(NULL),
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
	fHistoTrueSecondaryPi0FromEtaInvMassPt(NULL),
	fHistoTrueEtaWithPi0DaughterMCPt(NULL),
	fHistoTrueSecondaryPi0FromLambdaInvMassPt(NULL),
	fHistoTrueLambdaWithPi0DaughterMCPt(NULL),
	fHistoTrueBckGGInvMassPt(NULL),
	fHistoTrueBckContInvMassPt(NULL),
	fHistoTruePi0PtY(NULL),
	fHistoTrueEtaPtY(NULL),
	fHistoTruePi0PtAlpha(NULL),
	fHistoTrueEtaPtAlpha(NULL),
	fHistoTruePi0PtOpenAngle(NULL),
	fHistoTrueEtaPtOpenAngle(NULL),
	fHistoClusPhotonBGPt(NULL),
	fHistoClusPhotonPlusConvBGPt(NULL),
	fHistoTrueClusGammaPt(NULL),
	fHistoTrueClusUnConvGammaPt(NULL),
	fHistoTrueClusUnConvGammaMCPt(NULL),
	fHistoTrueClusElectronPt(NULL),
	fHistoTrueClusConvGammaPt(NULL),
	fHistoTrueClusConvGammaMCPt(NULL),
	fHistoTrueClusConvGammaFullyPt(NULL),
	fHistoTrueClusMergedGammaPt(NULL),
	fHistoTrueClusMergedPartConvGammaPt(NULL),
	fHistoTrueClusDalitzPt(NULL),
	fHistoTrueClusDalitzMergedPt(NULL),
	fHistoTrueClusPhotonFromElecMotherPt(NULL),
	fHistoTrueClusShowerPt(NULL),
	fHistoTrueClusSubLeadingPt(NULL),
	fHistoTrueClusNParticles(NULL),
	fHistoTrueClusEMNonLeadingPt(NULL),
	fHistoTrueNLabelsInClus(NULL),
	fHistoTruePrimaryClusGammaPt(NULL),
	fHistoTruePrimaryClusGammaESDPtMCPt(NULL),
	fHistoTruePrimaryClusConvGammaPt(NULL),
	fHistoTruePrimaryClusConvGammaESDPtMCPt(NULL),
	fHistoTrueSecondaryClusGammaPt(NULL),
	fHistoTrueSecondaryClusConvGammaPt(NULL),
	fHistoTrueSecondaryClusGammaFromXFromK0sPt(NULL),
	fHistoTrueSecondaryClusConvGammaFromXFromK0sPt(NULL),
	fHistoTrueSecondaryClusGammaFromXFromLambdaPt(NULL),
	fHistoTrueSecondaryClusConvGammaFromXFromLambdaPt(NULL),
	fHistoTrueSecondaryClusGammaFromXFromEtasPt(NULL),
	fHistoTrueSecondaryClusConvGammaFromXFromEtasPt(NULL),
	fHistoDoubleCountTruePi0InvMassPt(NULL),
	fHistoDoubleCountTrueEtaInvMassPt(NULL),
	fVectorDoubleCountTruePi0s(0),
	fVectorDoubleCountTrueEtas(0),
	fHistoNEvents(NULL),
	fHistoNGoodESDTracks(NULL),
	fHistoVertexZ(NULL),
	fHistoNGammaCandidates(NULL),
	fHistoNGoodESDTracksVsNGammaCanditates(NULL),
	fHistoNV0Tracks(NULL),
	fProfileEtaShift(NULL),
	fEventPlaneAngle(-100),
	fRandom(0),
	fnCuts(0),
	fiCut(0),
	fIsHeavyIon(0),
	fDoMesonAnalysis(kTRUE),
	fDoMesonQA(0),
	fDoClusterQA(0),
	fIsFromMBHeader(kTRUE),
	fIsOverlappingWithOtherHeader(kFALSE),
	fIsMC(kFALSE),
	fDoTHnSparse(kTRUE)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaCalo::~AliAnalysisTaskGammaCalo()
{
	if(fClusterCandidates){
		delete fClusterCandidates;
		fClusterCandidates = 0x0;
	}
	if(fBGHandler){
		delete[] fBGHandler;
		fBGHandler = 0x0;
	}
}
//___________________________________________________________
void AliAnalysisTaskGammaCalo::InitBack(){
	
	const Int_t nDim = 4;
	Int_t nBins[nDim] = {800,350,7,6};
	Double_t xMin[nDim] = {0,0, 0,0};
	Double_t xMax[nDim] = {0.8,35,7,6};
	
    if(fDoTHnSparse){
        fSparseMotherInvMassPtZM = new THnSparseF*[fnCuts];
        fSparseMotherBackInvMassPtZM = new THnSparseF*[fnCuts];
    }

	fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];

	
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		if (((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
			TString cutstringEvent 	= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
			TString cutstringCalo 	= ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
			TString cutstringMeson 	= ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
			
			Int_t collisionSystem = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(0,1));
			Int_t centMin = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(1,1));
			Int_t centMax = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(2,1));
			
			if(collisionSystem == 1 || collisionSystem == 2 ||
				collisionSystem == 5 || collisionSystem == 8 ||
				collisionSystem == 9){
				centMin = centMin*10;
				centMax = centMax*10;
				if(centMax ==0 && centMax!=centMin) centMax=100;
			} else if(collisionSystem == 3 || collisionSystem == 6){
				centMin = centMin*5;
				centMax = centMax*5;
			} else if(collisionSystem == 4 || collisionSystem == 7){
				centMin = ((centMin*5)+45);
				centMax = ((centMax*5)+45);
			}
			
			if(fDoTHnSparse){
				fBackList[iCut] = new TList();
				fBackList[iCut]->SetName(Form("%s_%s_%s Back histograms",cutstringEvent.Data(),cutstringCalo.Data(), cutstringMeson.Data()));
				fBackList[iCut]->SetOwner(kTRUE);
				fCutFolder[iCut]->Add(fBackList[iCut]);

				fSparseMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_m","Back_Back_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
				fBackList[iCut]->Add(fSparseMotherBackInvMassPtZM[iCut]);

				fMotherList[iCut] = new TList();
				fMotherList[iCut]->SetName(Form("%s_%s_%s Mother histograms",cutstringEvent.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
				fMotherList[iCut]->SetOwner(kTRUE);
				fCutFolder[iCut]->Add(fMotherList[iCut]);

				fSparseMotherInvMassPtZM[iCut] = new THnSparseF("Back_Mother_InvMass_Pt_z_m","Back_Mother_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
				fMotherList[iCut]->Add(fSparseMotherInvMassPtZM[iCut]);
			}

			if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
				fBGHandler[iCut] = new AliGammaConversionAODBGHandler(
																	collisionSystem,centMin,centMax,
																	((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
																	((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity(),
																	4,8,7);
			}
		}
	}
}
//________________________________________________________________________
void AliAnalysisTaskGammaCalo::UserCreateOutputObjects(){
  
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
  
	fCutFolder = new TList*[fnCuts];
	fESDList = new TList*[fnCuts];
    if(fDoTHnSparse){
        fBackList = new TList*[fnCuts];
        fMotherList = new TList*[fnCuts];
    }
	fHistoNEvents = new TH1I*[fnCuts];
	fHistoNGoodESDTracks = new TH1I*[fnCuts];
	fHistoVertexZ = new TH1F*[fnCuts];
	fHistoNGammaCandidates = new TH1I*[fnCuts];
	fHistoNGoodESDTracksVsNGammaCanditates = new TH2F*[fnCuts];
	fHistoNV0Tracks = new TH1I*[fnCuts];
	fProfileEtaShift = new TProfile*[fnCuts];
  	
	if(fDoMesonAnalysis){
		fHistoMotherInvMassPt = new TH2F*[fnCuts];
		fHistoMotherInvMass3ClusterPt = new TH2F*[fnCuts];
		fHistoMotherBackInvMassPt = new TH2F*[fnCuts];
		fHistoMotherInvMassEalpha = new TH2F*[fnCuts];
		if (fDoMesonQA > 0){
			fHistoMotherPi0PtY =  new TH2F*[fnCuts];
			fHistoMotherEtaPtY =  new TH2F*[fnCuts];
			fHistoMotherPi0PtAlpha =  new TH2F*[fnCuts];
			fHistoMotherEtaPtAlpha =  new TH2F*[fnCuts];
			fHistoMotherPi0PtOpenAngle =  new TH2F*[fnCuts];
			fHistoMotherEtaPtOpenAngle =  new TH2F*[fnCuts];
		}
		if(fDoMesonQA == 1){
			fHistoMotherInvMassECalib = new TH2F*[fnCuts];
			fHistoMotherInvMassECalibalpha = new TH2F*[fnCuts];
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
    
		fHistoNEvents[iCut] = new TH1I("NEvents","NEvents",10,-0.5,9.5);
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
		fESDList[iCut]->Add(fHistoNEvents[iCut]);
		
		if(fIsHeavyIon == 1) fHistoNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",4000,0,4000);
		else if(fIsHeavyIon == 2) fHistoNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",400,0,400);
		else fHistoNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",200,0,200);
		fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);

		fHistoVertexZ[iCut] = new TH1F("VertexZ","VertexZ",1000,-50,50);
		fESDList[iCut]->Add(fHistoVertexZ[iCut]);

		if(fIsHeavyIon == 1) fHistoNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",600,0,600);
		else if(fIsHeavyIon == 2) fHistoNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",400,0,400);
		else fHistoNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",100,0,100);
		fESDList[iCut]->Add(fHistoNGammaCandidates[iCut]);
		if(fIsHeavyIon == 1) fHistoNGoodESDTracksVsNGammaCanditates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",4000,0,4000,600,0,600);
		else if(fIsHeavyIon == 2) fHistoNGoodESDTracksVsNGammaCanditates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",400,0,400,400,0,400);
		else fHistoNGoodESDTracksVsNGammaCanditates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",200,0,200,100,0,100);
		fESDList[iCut]->Add(fHistoNGoodESDTracksVsNGammaCanditates[iCut]);
    
		
		if(fIsHeavyIon == 1) fHistoNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",30000,0,30000);
		else if(fIsHeavyIon == 2) fHistoNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",2500,0,2500);
		else fHistoNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",1500,0,1500);
		fESDList[iCut]->Add(fHistoNV0Tracks[iCut]);
		fProfileEtaShift[iCut] = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
		fESDList[iCut]->Add(fProfileEtaShift[iCut]);
    
   		fHistoClusGammaPt[iCut] = new TH1F("ClusGamma_Pt","ClusGamma_Pt",350,0,35);
		fESDList[iCut]->Add(fHistoClusGammaPt[iCut]);
		fHistoClusOverlapHeadersGammaPt[iCut] = new TH1F("ClusGammaOverlapHeaders_Pt","ClusGammaOverlapHeaders_Pt",350,0,35);
		fESDList[iCut]->Add(fHistoClusOverlapHeadersGammaPt[iCut]);
		
		if(fDoMesonAnalysis){
			fHistoMotherInvMassPt[iCut] = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",800,0,0.8,350,0,35);
			fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);
			fHistoMotherInvMass3ClusterPt[iCut] = new TH2F("ESD_Mother_InvMass3Cluster_Pt","ESD_Mother_InvMass3Cluster_Pt",800,0,0.8,350,0,35);
			fESDList[iCut]->Add(fHistoMotherInvMass3ClusterPt[iCut]);
			fHistoMotherBackInvMassPt[iCut] = new TH2F("ESD_Background_InvMass_Pt","ESD_Background_InvMass_Pt",800,0,0.8,350,0,35);
			fESDList[iCut]->Add(fHistoMotherBackInvMassPt[iCut]);
			fHistoMotherInvMassEalpha[iCut] = new TH2F("ESD_Mother_InvMass_vs_E_alpha","ESD_Mother_InvMass_vs_E_alpha",800,0,0.8,350,0,35);
			fESDList[iCut]->Add(fHistoMotherInvMassEalpha[iCut]);
			if(fDoMesonQA == 1){
				fHistoMotherInvMassECalib[iCut] = new TH2F("ESD_Mother_InvMass_Pt_Calib","ESD_Mother_InvMass_Pt_Calib",800,0,0.8,350,0,35);
				fESDList[iCut]->Add(fHistoMotherInvMassECalib[iCut]);
				fHistoMotherInvMassECalibalpha[iCut] = new TH2F("ESD_Mother_InvMass_vs_E_Calib_alpha","ESD_Mother_InvMass_vs_E_Calib_alpha",800,0,0.8,350,0,35);
				fESDList[iCut]->Add(fHistoMotherInvMassECalibalpha[iCut]);
			}

			if (fDoMesonQA > 0 ){
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
				fHistoMotherPi0PtOpenAngle[iCut] = new TH2F("ESD_MotherPi0_Pt_OpenAngle","ESD_MotherPi0_Pt_OpenAngle",350,0.03,35.,100,0,TMath::Pi());
				SetLogBinningXTH2(fHistoMotherPi0PtOpenAngle[iCut]);
				fESDList[iCut]->Add(fHistoMotherPi0PtOpenAngle[iCut]);
				fHistoMotherEtaPtOpenAngle[iCut] = new TH2F("ESD_MotherEta_Pt_OpenAngle","ESD_MotherEta_Pt_OpenAngle",350,0.03,35.,100,0,TMath::Pi());
				SetLogBinningXTH2(fHistoMotherEtaPtOpenAngle[iCut]);
				fESDList[iCut]->Add(fHistoMotherEtaPtOpenAngle[iCut]);
			}
		}    
	}
	if(fDoMesonAnalysis){
		InitBack(); // Init Background Handler
	}
  
	if(fIsMC){
		// MC Histogramms
		fMCList 	= new TList*[fnCuts];
		// True Histogramms
		fTrueList 	= new TList*[fnCuts];
		// Selected Header List
		fHeaderNameList 					= new TList*[fnCuts];
		fHistoMCHeaders 					= new TH1I*[fnCuts];
		fHistoMCAllGammaPt 					= new TH1F*[fnCuts];
		fHistoMCDecayGammaPi0Pt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaRhoPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaEtaPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaOmegaPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaEtapPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaPhiPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaSigmaPt 			= new TH1F*[fnCuts];
	
		fHistoClusPhotonBGPt 							= new TH2F*[fnCuts];
		fHistoClusPhotonPlusConvBGPt					= new TH2F*[fnCuts];
	
		fHistoTrueClusGammaPt 								= new TH1F*[fnCuts];
		fHistoTruePrimaryClusGammaPt 						= new TH1F*[fnCuts];
		fHistoTruePrimaryClusGammaESDPtMCPt 				= new TH2F*[fnCuts];
		fHistoTruePrimaryClusConvGammaPt 					= new TH1F*[fnCuts];
		fHistoTruePrimaryClusConvGammaESDPtMCPt 			= new TH2F*[fnCuts];
		fHistoTrueSecondaryClusGammaPt 						= new TH1F*[fnCuts];
		fHistoTrueSecondaryClusConvGammaPt 					= new TH1F*[fnCuts];
		fHistoTrueSecondaryClusGammaFromXFromK0sPt 			= new TH1F*[fnCuts];
		fHistoTrueSecondaryClusConvGammaFromXFromK0sPt 		= new TH1F*[fnCuts];
		fHistoTrueSecondaryClusGammaFromXFromLambdaPt 		= new TH1F*[fnCuts];
		fHistoTrueSecondaryClusConvGammaFromXFromLambdaPt 	= new TH1F*[fnCuts];
		fHistoTrueSecondaryClusGammaFromXFromEtasPt 		= new TH1F*[fnCuts];
		fHistoTrueSecondaryClusConvGammaFromXFromEtasPt 	= new TH1F*[fnCuts];
    
		if (fDoClusterQA > 0){	
			fHistoTrueClusUnConvGammaPt 		= new TH1F*[fnCuts];
			fHistoTrueClusUnConvGammaMCPt 		= new TH1F*[fnCuts];
			fHistoTrueClusElectronPt 			= new TH1F*[fnCuts];
			fHistoTrueClusConvGammaPt 			= new TH1F*[fnCuts];
			fHistoTrueClusConvGammaMCPt 			= new TH1F*[fnCuts];
			fHistoTrueClusConvGammaFullyPt 		= new TH1F*[fnCuts];
			fHistoTrueClusMergedGammaPt 		= new TH1F*[fnCuts];
			fHistoTrueClusMergedPartConvGammaPt = new TH1F*[fnCuts];
			fHistoTrueClusDalitzPt			 	= new TH1F*[fnCuts];
			fHistoTrueClusDalitzMergedPt	 	= new TH1F*[fnCuts];
			fHistoTrueClusPhotonFromElecMotherPt= new TH1F*[fnCuts];
			fHistoTrueClusShowerPt				= new TH1F*[fnCuts];
			fHistoTrueClusSubLeadingPt			= new TH1F*[fnCuts];
			fHistoTrueClusNParticles			= new TH1I*[fnCuts];
			fHistoTrueClusEMNonLeadingPt		= new TH1F*[fnCuts];
			fHistoTrueNLabelsInClus 			= new TH1F*[fnCuts];			
		}
    
		if(fDoMesonAnalysis){
			fHistoMCPi0Pt 					= new TH1F*[fnCuts];
			fHistoMCPi0WOWeightPt 			= new TH1F*[fnCuts];
			fHistoMCEtaPt 					= new TH1F*[fnCuts];
			fHistoMCEtaWOWeightPt 			= new TH1F*[fnCuts];
			fHistoMCPi0InAccPt 				= new TH1F*[fnCuts];
			fHistoMCEtaInAccPt 				= new TH1F*[fnCuts];
      
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
			fHistoTrueSecondaryPi0FromEtaInvMassPt 	= new TH2F*[fnCuts];
			fHistoTrueSecondaryPi0FromLambdaInvMassPt = new TH2F*[fnCuts];
			if (fDoMesonQA > 0){
				fHistoMCPi0PtY 								= new TH2F*[fnCuts];
				fHistoMCEtaPtY 								= new TH2F*[fnCuts];
				fHistoMCPi0PtAlpha 							= new TH2F*[fnCuts];
				fHistoMCEtaPtAlpha 							= new TH2F*[fnCuts];
				fHistoMCK0sPt 								= new TH1F*[fnCuts];
				fHistoMCK0sWOWeightPt 						= new TH1F*[fnCuts];
				fHistoMCK0sPtY	 							= new TH2F*[fnCuts];
				fHistoMCSecPi0PtvsSource 					= new TH2F*[fnCuts];
				fHistoMCSecPi0Source 						= new TH1F*[fnCuts];
				fHistoMCSecEtaPt 							= new TH1F*[fnCuts];
				fHistoMCSecEtaSource 						= new TH1F*[fnCuts];
				fHistoTruePi0CaloPhotonInvMassPt			= new TH2F*[fnCuts];
				fHistoTrueEtaCaloPhotonInvMassPt			= new TH2F*[fnCuts];
				fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt= new TH2F*[fnCuts];
				fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt= new TH2F*[fnCuts];
				fHistoTruePi0CaloConvertedPhotonInvMassPt 	= new TH2F*[fnCuts];
				fHistoTrueEtaCaloConvertedPhotonInvMassPt	= new TH2F*[fnCuts];
				fHistoTruePi0CaloElectronInvMassPt		= new TH2F*[fnCuts];
				fHistoTrueEtaCaloElectronInvMassPt		= new TH2F*[fnCuts];
				fHistoTruePi0CaloMergedClusterInvMassPt 	= new TH2F*[fnCuts];
				fHistoTrueEtaCaloMergedClusterInvMassPt 	= new TH2F*[fnCuts];
				fHistoTruePi0CaloMergedClusterPartConvInvMassPt 	= new TH2F*[fnCuts];
				fHistoTrueEtaCaloMergedClusterPartConvInvMassPt 	= new TH2F*[fnCuts];
				fHistoTruePi0NonMergedElectronPhotonInvMassPt 	= new TH2F*[fnCuts];
				fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt 	= new TH2F*[fnCuts];
				fHistoTruePrimaryPi0MCPtResolPt 			= new TH2F*[fnCuts];
				fHistoTruePrimaryEtaMCPtResolPt 			= new TH2F*[fnCuts];
				fHistoTrueK0sWithPi0DaughterMCPt 			= new TH1F*[fnCuts];
				fHistoTrueEtaWithPi0DaughterMCPt 			= new TH1F*[fnCuts];
				fHistoTrueLambdaWithPi0DaughterMCPt 		= new TH1F*[fnCuts];
				fHistoTrueBckGGInvMassPt 					= new TH2F*[fnCuts];
				fHistoTrueBckContInvMassPt 					= new TH2F*[fnCuts];
				fHistoTruePi0PtY 							= new TH2F*[fnCuts];
				fHistoTrueEtaPtY 							= new TH2F*[fnCuts];
				fHistoTruePi0PtAlpha 						= new TH2F*[fnCuts];
				fHistoTrueEtaPtAlpha 						= new TH2F*[fnCuts];
				fHistoTruePi0PtOpenAngle 					= new TH2F*[fnCuts];
				fHistoTrueEtaPtOpenAngle 					= new TH2F*[fnCuts];
			}
			if (fDoMesonQA==2){
				fHistoTruePi0Category1						= new TH2F*[fnCuts];
				fHistoTrueEtaCategory1						= new TH2F*[fnCuts];				
				fHistoTruePi0Category2						= new TH2F*[fnCuts];
				fHistoTrueEtaCategory2						= new TH2F*[fnCuts];				
				fHistoTruePi0Category3						= new TH2F*[fnCuts];
				fHistoTrueEtaCategory3						= new TH2F*[fnCuts];				
				fHistoTruePi0Category4_6					= new TH2F*[fnCuts];
				fHistoTrueEtaCategory4_6					= new TH2F*[fnCuts];				
				fHistoTruePi0Category5						= new TH2F*[fnCuts];
				fHistoTrueEtaCategory5						= new TH2F*[fnCuts];				
				fHistoTruePi0Category7						= new TH2F*[fnCuts];
				fHistoTrueEtaCategory7						= new TH2F*[fnCuts];				
				fHistoTruePi0Category8						= new TH2F*[fnCuts];
				fHistoTrueEtaCategory8						= new TH2F*[fnCuts];				
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
			fHistoMCAllGammaPt[iCut] = new TH1F("MC_AllGamma_Pt","MC_AllGamma_Pt",350,0,35);
			fMCList[iCut]->Add(fHistoMCAllGammaPt[iCut]);
			fHistoMCDecayGammaPi0Pt[iCut] = new TH1F("MC_DecayGammaPi0_Pt","MC_DecayGammaPi0_Pt",350,0,35);
			fMCList[iCut]->Add(fHistoMCDecayGammaPi0Pt[iCut]);
			fHistoMCDecayGammaRhoPt[iCut] = new TH1F("MC_DecayGammaRho_Pt","MC_DecayGammaRho_Pt",350,0,35);
			fMCList[iCut]->Add(fHistoMCDecayGammaRhoPt[iCut]);
			fHistoMCDecayGammaEtaPt[iCut] = new TH1F("MC_DecayGammaEta_Pt","MC_DecayGammaEta_Pt",350,0,35);
			fMCList[iCut]->Add(fHistoMCDecayGammaEtaPt[iCut]);
			fHistoMCDecayGammaOmegaPt[iCut] = new TH1F("MC_DecayGammaOmega_Pt","MC_DecayGammaOmmega_Pt",350,0,35);
			fMCList[iCut]->Add(fHistoMCDecayGammaOmegaPt[iCut]);
			fHistoMCDecayGammaEtapPt[iCut] = new TH1F("MC_DecayGammaEtap_Pt","MC_DecayGammaEtap_Pt",350,0,35);
			fMCList[iCut]->Add(fHistoMCDecayGammaEtapPt[iCut]);
			fHistoMCDecayGammaPhiPt[iCut] = new TH1F("MC_DecayGammaPhi_Pt","MC_DecayGammaPhi_Pt",350,0,35);
			fMCList[iCut]->Add(fHistoMCDecayGammaPhiPt[iCut]);
			fHistoMCDecayGammaSigmaPt[iCut] = new TH1F("MC_DecayGammaSigma_Pt","MC_DecayGammaSigma_Pt",350,0,35);
			fMCList[iCut]->Add(fHistoMCDecayGammaSigmaPt[iCut]);
			
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
				if (fDoMesonQA > 0){
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

					fHistoMCK0sPt[iCut] = new TH1F("MC_K0s_Pt","MC_K0s_Pt",350,0,35);
					fHistoMCK0sPt[iCut]->Sumw2();
					fMCList[iCut]->Add(fHistoMCK0sPt[iCut]);
					fHistoMCK0sWOWeightPt[iCut] = new TH1F("MC_K0s_WOWeights_Pt","MC_K0s_WOWeights_Pt",350,0,35);
					fHistoMCK0sWOWeightPt[iCut]->Sumw2();
					fMCList[iCut]->Add(fHistoMCK0sWOWeightPt[iCut]);
					fHistoMCK0sPtY[iCut] = new TH2F("MC_K0s_Pt_Y","MC_K0s_Pt_Y",350,0.03,35.,150,-1.5,1.5);
					fHistoMCK0sPtY[iCut]->Sumw2();
					SetLogBinningXTH2(fHistoMCK0sPtY[iCut]);
					fMCList[iCut]->Add(fHistoMCK0sPtY[iCut]);
					
					fHistoMCSecPi0Source[iCut] = new TH1F("MC_SecPi0_Source","MC_SecPi0_Source",5000,0.,5000);
					fMCList[iCut]->Add(fHistoMCSecPi0Source[iCut]);
					fHistoMCSecEtaSource[iCut] = new TH1F("MC_SecEta_Source","MC_SecEta_Source",5000,0,5000);
					fMCList[iCut]->Add(fHistoMCSecEtaSource[iCut]);
					fHistoMCSecPi0PtvsSource[iCut] = new TH2F("MC_SecPi0_Pt_Source","MC_SecPi0_Pt_Source",350,0.0,35.,16,-0.5,15.5);
					fHistoMCSecPi0PtvsSource[iCut]->Sumw2();
					fMCList[iCut]->Add(fHistoMCSecPi0PtvsSource[iCut]);
					fHistoMCSecEtaPt[iCut] = new TH1F("MC_SecEta_Pt","MC_SecEta_Pt",350,0,35);
					fHistoMCSecEtaPt[iCut]->Sumw2();
					fMCList[iCut]->Add(fHistoMCSecEtaPt[iCut]);
				}
        
			}
			fTrueList[iCut] = new TList();
			fTrueList[iCut]->SetName(Form("%s_%s_%s True histograms",cutstringEvent.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
			fTrueList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fTrueList[iCut]);
                  
			fHistoClusPhotonBGPt[iCut] = new TH2F("ESD_TrueClusPhotonBG_Pt","ESD_TrueClusPhotonBG_Pt",350,0,35,9,-0.5,8.5);
			fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 1,"Elec");
			fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 2,"Pion");
			fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 3,"Proton");
			fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 4,"Kaon");
			fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 5,"Neutron");
			fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 6,"K0s");
			fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 7,"Lambda");
			fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 8,"Muon");
			fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel(9,"Rest");
			fTrueList[iCut]->Add(fHistoClusPhotonBGPt[iCut]);
			fHistoClusPhotonPlusConvBGPt[iCut] = new TH2F("ESD_TrueClusPhotonPlusConvBG_Pt","ESD_TrueClusPhotonPlusConvBG_Pt",350,0,35,9,-0.5,8.5);
			fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 1,"Elec");
			fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 2,"Pion");
			fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 3,"Proton");
			fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 4,"Kaon");
			fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 5,"Neutron");
			fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 6,"K0s");
			fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 7,"Lambda");
			fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 8,"Muon");
			fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel(9,"Rest");
			fTrueList[iCut]->Add(fHistoClusPhotonPlusConvBGPt[iCut]);
		
			fHistoTrueClusGammaPt[iCut] = new TH1F("TrueClusGamma_Pt","ESD_TrueClusGamma_Pt",350,0,35);
			fTrueList[iCut]->Add(fHistoTrueClusGammaPt[iCut]);
			fHistoTruePrimaryClusGammaPt[iCut] = new TH1F("TruePrimaryClusGamma_Pt","ESD_TruePrimaryClusGamma_Pt",350,0,35);
			fTrueList[iCut]->Add(fHistoTruePrimaryClusGammaPt[iCut]);
			fHistoTruePrimaryClusGammaESDPtMCPt[iCut] = new TH2F("TruePrimaryClusGamma_Pt_MCPt","ESD_TruePrimaryClusGamma_MCPt",350,0,35,350,0,35);
			fTrueList[iCut]->Add(fHistoTruePrimaryClusGammaESDPtMCPt[iCut]);
			fHistoTruePrimaryClusConvGammaPt[iCut] = new TH1F("TruePrimaryClusConvGamma_Pt","ESD_TruePrimaryClusConvGamma_Pt",350,0,35);
			fTrueList[iCut]->Add(fHistoTruePrimaryClusConvGammaPt[iCut]);
			fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut] = new TH2F("TruePrimaryClusConvGamma_Pt_MCPt","ESD_TruePrimaryClusConvGamma_MCPt",350,0,35,350,0,35);
			fTrueList[iCut]->Add(fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut]);
			fHistoTrueSecondaryClusGammaPt[iCut] = new TH1F("ESD_TrueSecondaryClusGamma_Pt","ESD_TrueSecondaryClusGamma_Pt",350,0,35);
			fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaPt[iCut]);      
			fHistoTrueSecondaryClusConvGammaPt[iCut] = new TH1F("ESD_TrueSecondaryClusConvGamma_Pt","ESD_TrueSecondaryClusConvGamma_Pt",350,0,35);
			fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaPt[iCut]);      

			fHistoTrueSecondaryClusGammaFromXFromK0sPt[iCut] = new TH1F("ESD_TrueSecondaryClusGammaFromXFromK0s_Pt", "ESD_TrueSecondaryClusGammaFromXFromK0s_Pt",350,0,35);
			fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaFromXFromK0sPt[iCut]);
			fHistoTrueSecondaryClusConvGammaFromXFromK0sPt[iCut] = new TH1F("ESD_TrueSecondaryClusConvGammaFromXFromK0s_Pt", "ESD_TrueSecondaryClusConvGammaFromXFromK0s_Pt",350,0,35);
			fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaFromXFromK0sPt[iCut]);
			fHistoTrueSecondaryClusGammaFromXFromLambdaPt[iCut] = new TH1F("ESD_TrueSecondaryClusGammaFromXFromLambda_Pt", "ESD_TrueSecondaryClusGammaFromXFromLambda_Pt",350,0,35);
			fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaFromXFromLambdaPt[iCut]);
			fHistoTrueSecondaryClusConvGammaFromXFromLambdaPt[iCut] = new TH1F("ESD_TrueSecondaryClusConvGammaFromXFromLambda_Pt", "ESD_TrueSecondaryClusConvGammaFromXFromLambda_Pt",350,0,35);
			fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaFromXFromLambdaPt[iCut]);
			fHistoTrueSecondaryClusGammaFromXFromEtasPt[iCut] = new TH1F("ESD_TrueSecondaryClusGammaFromXFromEtas_Pt", "ESD_TrueSecondaryClusGammaFromXFromEtas_Pt",350,0,35);
			fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaFromXFromEtasPt[iCut]);
			fHistoTrueSecondaryClusConvGammaFromXFromEtasPt[iCut] = new TH1F("ESD_TrueSecondaryClusConvGammaFromXFromEtas_Pt", "ESD_TrueSecondaryClusConvGammaFromXFromEtas_Pt",350,0,35);
			fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaFromXFromEtasPt[iCut]);

			
			if (fDoClusterQA > 0){	
				fHistoTrueClusUnConvGammaPt[iCut] = new TH1F("TrueClusUnConvGamma_Pt","TrueClusUnConvGamma_Pt",350,0,35);
				fTrueList[iCut]->Add(fHistoTrueClusUnConvGammaPt[iCut]);
				fHistoTrueClusUnConvGammaMCPt[iCut] = new TH1F("TrueClusUnConvGamma_MCPt","TrueClusUnConvGamma_MCPt",350,0,35);
				fTrueList[iCut]->Add(fHistoTrueClusUnConvGammaMCPt[iCut]);
				fHistoTrueClusElectronPt[iCut] = new TH1F("TrueClusElectron_Pt","TrueElectronGamma_Pt",350,0,35);
				fTrueList[iCut]->Add(fHistoTrueClusElectronPt[iCut]);
				fHistoTrueClusConvGammaPt[iCut] = new TH1F("TrueClusConvGamma_Pt","TrueClusConvGamma_Pt",350,0,35);
				fTrueList[iCut]->Add(fHistoTrueClusConvGammaPt[iCut]);
				fHistoTrueClusConvGammaMCPt[iCut] = new TH1F("TrueClusConvGamma_MCPt","TrueClusConvGamma_MCPt",350,0,35);
				fTrueList[iCut]->Add(fHistoTrueClusConvGammaMCPt[iCut]);
				fHistoTrueClusConvGammaFullyPt[iCut] = new TH1F("TrueClusConvGammaFullyContained_Pt","TrueClusConvGammaFullyContained_Pt",350,0,35);
				fTrueList[iCut]->Add(fHistoTrueClusConvGammaFullyPt[iCut]);
				fHistoTrueClusMergedGammaPt[iCut] = new TH1F("TrueClusMergedGamma_Pt","TrueClusMergedGamma_Pt",350,0,35);
				fTrueList[iCut]->Add(fHistoTrueClusMergedGammaPt[iCut]);
				fHistoTrueClusMergedPartConvGammaPt[iCut] = new TH1F("TrueClusMergedPartConvGamma_Pt","TrueClusMergedPartConvGamma_Pt",350,0,35);
				fTrueList[iCut]->Add(fHistoTrueClusMergedPartConvGammaPt[iCut]);
				fHistoTrueClusDalitzPt[iCut] = new TH1F("TrueClusDalitz_Pt","TrueClusDalitz_Pt",350,0,35);
				fTrueList[iCut]->Add(fHistoTrueClusDalitzPt[iCut]);
				fHistoTrueClusDalitzMergedPt[iCut] = new TH1F("TrueClusDalitzMerged_Pt","TrueClusDalitzMerged_Pt",350,0,35);
				fTrueList[iCut]->Add(fHistoTrueClusDalitzMergedPt[iCut]);
				fHistoTrueClusPhotonFromElecMotherPt[iCut] = new TH1F("TrueClusPhotonFromElecMother_Pt","TrueClusPhotonFromElecMother_Pt",350,0,35);
				fTrueList[iCut]->Add(fHistoTrueClusPhotonFromElecMotherPt[iCut]);
				fHistoTrueClusShowerPt[iCut] = new TH1F("TrueClusShower_Pt","TrueClusShower_Pt",350,0,35);
				fTrueList[iCut]->Add(fHistoTrueClusShowerPt[iCut]);
				fHistoTrueClusSubLeadingPt[iCut] = new TH1F("TrueClusSubleading_Pt","TrueClusSubleading_Pt",350,0,35);
				fTrueList[iCut]->Add(fHistoTrueClusSubLeadingPt[iCut]);
				fHistoTrueClusNParticles[iCut] = new TH1I("TrueClusNParticles","TrueClusNParticles",20,0,20);
				fTrueList[iCut]->Add(fHistoTrueClusNParticles[iCut]);
				fHistoTrueClusEMNonLeadingPt[iCut] = new TH1F("TrueClusEMNonLeading_Pt","TrueClusEMNonLeading_Pt",350,0,35);
				fTrueList[iCut]->Add(fHistoTrueClusEMNonLeadingPt[iCut]);
				fHistoTrueNLabelsInClus[iCut] = new TH1F("TrueNLabelsInClus","TrueNLabelsInClus",100,-0.5,99.5);
				fTrueList[iCut]->Add(fHistoTrueNLabelsInClus[iCut]);	
			}	

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
				fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromEta_InvMass_Pt","ESD_TrueSecondaryPi0FromEta_InvMass_Pt",800,0,0.8,350,0,35);
				fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]);
				fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromLambda_InvMass_Pt","ESD_TrueSecondaryPi0FromLambda_InvMass_Pt",800,0,0.8,350,0,35);
				fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]);
				if (fDoMesonQA > 0){
					fHistoTruePi0CaloPhotonInvMassPt[iCut] = new TH2F("ESD_TruePi0CaloPhoton_InvMass_Pt","ESD_TruePi0CaloPhoton_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTruePi0CaloPhotonInvMassPt[iCut]);
					fHistoTrueEtaCaloPhotonInvMassPt[iCut] = new TH2F("ESD_TrueEtaCaloPhoton_InvMass_Pt","ESD_TrueEtaCaloPhoton_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTrueEtaCaloPhotonInvMassPt[iCut]);

					fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt[iCut] = new TH2F("ESD_TruePi0CaloMixedPhotonConvertedPhoton_InvMass_Pt","ESD_TruePi0CaloMixedPhotonConvertedPhoton_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt[iCut]);
					fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt[iCut] = new TH2F("ESD_TrueEtaCaloMixedPhotonConvertedPhoton_InvMass_Pt","ESD_TrueEtaCaloMixedPhotonConvertedPhoton_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt[iCut]);

					fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut] = new TH2F("ESD_TruePi0CaloConvertedPhoton_InvMass_Pt","ESD_TruePi0CaloConvertedPhoton_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut]);
					fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut] = new TH2F("ESD_TrueEtaCaloConvertedPhoton_InvMass_Pt","ESD_TrueEtaCaloConvertedPhoton_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut]);

					fHistoTruePi0CaloElectronInvMassPt[iCut] = new TH2F("ESD_TruePi0CaloElectron_InvMass_Pt","ESD_TruePi0CaloElectron_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTruePi0CaloElectronInvMassPt[iCut]);
					fHistoTrueEtaCaloElectronInvMassPt[iCut] = new TH2F("ESD_TrueEtaCaloElectron_InvMass_Pt","ESD_TrueEtaCaloElectron_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTrueEtaCaloElectronInvMassPt[iCut]);

					fHistoTruePi0CaloMergedClusterInvMassPt[iCut] = new TH2F("ESD_TruePi0CaloMergedCluster_InvMass_Pt","ESD_TruePi0CaloMergedCluster_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTruePi0CaloMergedClusterInvMassPt[iCut]);
					fHistoTrueEtaCaloMergedClusterInvMassPt[iCut] = new TH2F("ESD_TrueEtaCaloMergedCluster_InvMass_Pt","ESD_TrueEtaCaloMergedCluster_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTrueEtaCaloMergedClusterInvMassPt[iCut]);

					fHistoTruePi0CaloMergedClusterPartConvInvMassPt[iCut] = new TH2F("ESD_TruePi0CaloMergedClusterPartConv_InvMass_Pt","ESD_TruePi0CaloMergedClusterPartConv_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTruePi0CaloMergedClusterPartConvInvMassPt[iCut]);
					fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[iCut] = new TH2F("ESD_TrueEtaCaloMergedClusterPartConv_InvMass_Pt","ESD_TrueEtaCaloMergedClusterPartConv_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[iCut]);

					fHistoTruePi0NonMergedElectronPhotonInvMassPt[iCut] = new TH2F("ESD_TruePi0NonMergedElectronPhoton_InvMass_Pt","ESD_TruePi0NonMergedElectronPhoton_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTruePi0NonMergedElectronPhotonInvMassPt[iCut]);
					fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt[iCut] = new TH2F("ESD_TruePi0NonMergedElectronMergedPhoton_InvMass_Pt","ESD_TruePi0NonMergedElectronMergedPhoton_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt[iCut]);
					
					fHistoTruePrimaryPi0MCPtResolPt[iCut] = new TH2F("ESD_TruePrimaryPi0_MCPt_ResolPt","ESD_TruePrimaryPi0_ResolPt_MCPt",500,0.03,35,1000,-1.,1.);
					fHistoTruePrimaryPi0MCPtResolPt[iCut]->Sumw2();
					SetLogBinningXTH2(fHistoTruePrimaryPi0MCPtResolPt[iCut]);
					fTrueList[iCut]->Add(fHistoTruePrimaryPi0MCPtResolPt[iCut]);
					fHistoTruePrimaryEtaMCPtResolPt[iCut]  = new TH2F("ESD_TruePrimaryEta_MCPt_ResolPt","ESD_TruePrimaryEta_ResolPt_MCPt",500,0.03,35,1000,-1.,1.);
					fHistoTruePrimaryEtaMCPtResolPt[iCut]->Sumw2();
					SetLogBinningXTH2(fHistoTruePrimaryEtaMCPtResolPt[iCut]);
					fTrueList[iCut]->Add(fHistoTruePrimaryEtaMCPtResolPt[iCut]);
					fHistoTrueBckGGInvMassPt[iCut] = new TH2F("ESD_TrueBckGG_InvMass_Pt","ESD_TrueBckGG_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTrueBckGGInvMassPt[iCut]);
					fHistoTrueBckContInvMassPt[iCut] = new TH2F("ESD_TrueBckCont_InvMass_Pt","ESD_TrueBckCont_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTrueBckContInvMassPt[iCut]);
					fHistoTrueK0sWithPi0DaughterMCPt[iCut] = new TH1F("ESD_TrueK0sWithPi0Daughter_MCPt","ESD_TrueK0sWithPi0Daughter_MCPt",350,0,35);
					fTrueList[iCut]->Add(fHistoTrueK0sWithPi0DaughterMCPt[iCut]);
					fHistoTrueEtaWithPi0DaughterMCPt[iCut] = new TH1F("ESD_TrueEtaWithPi0Daughter_MCPt","ESD_TrueEtaWithPi0Daughter_MCPt",350,0,35);
					fTrueList[iCut]->Add(fHistoTrueEtaWithPi0DaughterMCPt[iCut]);
					fHistoTrueLambdaWithPi0DaughterMCPt[iCut] = new TH1F("ESD_TrueLambdaWithPi0Daughter_MCPt","ESD_TrueLambdaWithPi0Daughter_MCPt",350,0,35);
					fTrueList[iCut]->Add(fHistoTrueLambdaWithPi0DaughterMCPt[iCut]);
          
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
					
					fHistoTruePi0PtOpenAngle[iCut] = new TH2F("ESD_TruePi0_Pt_OpenAngle","ESD_TruePi0_Pt_OpenAngle",350,0.03,35.,200,0,2*TMath::Pi());
					SetLogBinningXTH2(fHistoTruePi0PtOpenAngle[iCut]);
					fTrueList[iCut]->Add(fHistoTruePi0PtOpenAngle[iCut]);
					fHistoTrueEtaPtOpenAngle[iCut] = new TH2F("ESD_TrueEta_Pt_OpenAngle","ESD_TrueEta_Pt_OpenAngle",350,0.03,35.,200,0,2*TMath::Pi());
					SetLogBinningXTH2(fHistoTrueEtaPtOpenAngle[iCut]);
					fTrueList[iCut]->Add(fHistoTrueEtaPtOpenAngle[iCut]);
				}
				
				if (fDoMesonQA == 2){
					fHistoTruePi0Category1[iCut] = new TH2F("ESD_TruePi0Category1_InvMass_Pt","ESD_TruePi0Category1_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTruePi0Category1[iCut]);
					fHistoTrueEtaCategory1[iCut] = new TH2F("ESD_TrueEtaCategory1_InvMass_Pt","ESD_TrueEtaCategory1_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTrueEtaCategory1[iCut]);					
					fHistoTruePi0Category2[iCut] = new TH2F("ESD_TruePi0Category2_InvMass_Pt","ESD_TruePi0Category2_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTruePi0Category2[iCut]);
					fHistoTrueEtaCategory2[iCut] = new TH2F("ESD_TrueEtaCategory2_InvMass_Pt","ESD_TrueEtaCategory2_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTrueEtaCategory2[iCut]);					
					fHistoTruePi0Category3[iCut] = new TH2F("ESD_TruePi0Category3_InvMass_Pt","ESD_TruePi0Category3_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTruePi0Category3[iCut]);
					fHistoTrueEtaCategory3[iCut] = new TH2F("ESD_TrueEtaCategory3_InvMass_Pt","ESD_TrueEtaCategory3_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTrueEtaCategory3[iCut]);					
					fHistoTruePi0Category4_6[iCut] = new TH2F("ESD_TruePi0Category4_6_InvMass_Pt","ESD_TruePi0Category4_6_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTruePi0Category4_6[iCut]);
					fHistoTrueEtaCategory4_6[iCut] = new TH2F("ESD_TrueEtaCategory4_6_InvMass_Pt","ESD_TrueEtaCategory4_6_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTrueEtaCategory4_6[iCut]);					
					fHistoTruePi0Category5[iCut] = new TH2F("ESD_TruePi0Category5_InvMass_Pt","ESD_TruePi0Category5_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTruePi0Category5[iCut]);
					fHistoTrueEtaCategory5[iCut] = new TH2F("ESD_TrueEtaCategory5_InvMass_Pt","ESD_TrueEtaCategory5_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTrueEtaCategory5[iCut]);					
					fHistoTruePi0Category7[iCut] = new TH2F("ESD_TruePi0Category7_InvMass_Pt","ESD_TruePi0Category7_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTruePi0Category7[iCut]);
					fHistoTrueEtaCategory7[iCut] = new TH2F("ESD_TrueEtaCategory7_InvMass_Pt","ESD_TrueEtaCategory7_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTrueEtaCategory7[iCut]);					
					fHistoTruePi0Category8[iCut] = new TH2F("ESD_TruePi0Category8_InvMass_Pt","ESD_TruePi0Category8_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTruePi0Category8[iCut]);
					fHistoTrueEtaCategory8[iCut] = new TH2F("ESD_TrueEtaCategory8_InvMass_Pt","ESD_TrueEtaCategory8_InvMass_Pt",800,0,0.8,350,0,35);
					fTrueList[iCut]->Add(fHistoTrueEtaCategory8[iCut]);					
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
Bool_t AliAnalysisTaskGammaCalo::Notify()
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
void AliAnalysisTaskGammaCalo::UserExec(Option_t *)
{
	//
	// Called for each event
	//
	Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
	if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1
		for(Int_t iCut = 0; iCut<fnCuts; iCut++){
		fHistoNEvents[iCut]->Fill(eventQuality);
		}
		return;
	}
	
	if(fIsMC) fMCEvent = MCEvent();
	if(fMCEvent == NULL) fIsMC = kFALSE;
	
	fInputEvent = InputEvent();
	
	if(fIsMC && fInputEvent->IsA()==AliESDEvent::Class()){
		fMCStack = fMCEvent->Stack();
		if(fMCStack == NULL) fIsMC = kFALSE;
	}
	
	// ------------------- BeginEvent ----------------------------
	
	AliEventplane *EventPlane = fInputEvent->GetEventplane();
	if(fIsHeavyIon ==1)fEventPlaneAngle = EventPlane->GetEventplane("V0",fInputEvent,2);
	else fEventPlaneAngle=0.0;
	
	for(Int_t iCut = 0; iCut<fnCuts; iCut++){
		
		fiCut = iCut;
		
		Bool_t isRunningEMCALrelAna = kFALSE;
		if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;
		
		Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon, isRunningEMCALrelAna);
		
		if(eventNotAccepted){
		// cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
			fHistoNEvents[iCut]->Fill(eventNotAccepted); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
			continue;
		}

		if(eventQuality != 0){// Event Not Accepted
			//cout << "event rejected due to: " <<eventQuality << endl;
			fHistoNEvents[iCut]->Fill(eventQuality);
			continue;
		}

		fHistoNEvents[iCut]->Fill(eventQuality); // Should be 0 here
		fHistoNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks());
		fHistoVertexZ[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetZ());
		if(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsHeavyIon() == 2)	fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A());
			else fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C());

		if(fIsMC){
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
		if(fIsMC){
		if(fInputEvent->IsA()==AliESDEvent::Class())
			ProcessMCParticles();
		if(fInputEvent->IsA()==AliAODEvent::Class())
			ProcessAODMCParticles();
		}
		
		// it is in the loop to have the same conversion cut string (used also for MC stuff that should be same for V0 and Cluster)
		ProcessClusters();  					// process calo clusters

		fHistoNGammaCandidates[iCut]->Fill(fClusterCandidates->GetEntries());
		fHistoNGoodESDTracksVsNGammaCanditates[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),fClusterCandidates->GetEntries());
		if(fDoMesonAnalysis){ // Meson Analysis
			
			CalculatePi0Candidates(); // Combine Gammas from conversion and from calo
			if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
				if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
					
					CalculateBackground(); // Combinatorial Background
					UpdateEventByEventData(); // Store Event for mixed Events
				}
				
			}
			fVectorDoubleCountTruePi0s.clear();
			fVectorDoubleCountTrueEtas.clear();
		}

		fClusterCandidates->Clear(); // delete cluster candidates
	}
	
	PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessClusters()
{
	
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
		fIsOverlappingWithOtherHeader = kFALSE;
		// test whether largest contribution to cluster orginates in added signals
		if (fIsMC && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() > 0){
			if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetCaloPhotonMCLabel(0), fMCStack, fInputEvent) == 0) fIsFromMBHeader = kFALSE;
			if (clus->GetNLabels()>1){
				Int_t* mclabelsCluster = clus->GetLabels();
				for (Int_t l = 1; l < (Int_t)clus->GetNLabels(); l++ ){
					if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(mclabelsCluster[l], fMCStack, fInputEvent) == 0) fIsOverlappingWithOtherHeader = kTRUE;
				}	
			}	
				
		}
		if (fIsFromMBHeader && !fIsOverlappingWithOtherHeader){
			fHistoClusGammaPt[fiCut]->Fill(PhotonCandidate->Pt());	
			fClusterCandidates->Add(PhotonCandidate);
		}
		if (fIsFromMBHeader && fIsOverlappingWithOtherHeader) fHistoClusOverlapHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
		 // if no second loop is required add to events good gammas
		
		if(fIsMC){
			if(fInputEvent->IsA()==AliESDEvent::Class()){
				ProcessTrueClusterCandidates(PhotonCandidate);
			} else {
				ProcessTrueClusterCandidatesAOD(PhotonCandidate);
			}	
		}
		
		delete tmpvec;
	}
	
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessTrueClusterCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
		
	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();

	TParticle *Photon = NULL;
	if (!TruePhotonCandidate->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set task will abort");
	if (fDoClusterQA > 0) fHistoTrueNLabelsInClus[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels());
	
	if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0)Photon = fMCStack->Particle(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
		else return;
		
	if(Photon == NULL){
	//    cout << "no photon" << endl;
		return;
	}

	Int_t pdgCodeParticle = Photon->GetPdgCode();
	TruePhotonCandidate->SetCaloPhotonMCFlags(fMCStack);
	
	// True Photon
	if(fIsFromMBHeader && !fIsOverlappingWithOtherHeader){
		if (TruePhotonCandidate->IsLargestComponentPhoton() || TruePhotonCandidate->IsLargestComponentElectron() )fHistoTrueClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			else fHistoTrueClusEMNonLeadingPt[fiCut]->Fill(TruePhotonCandidate->Pt());
		if (fDoClusterQA > 0){
			if (TruePhotonCandidate->IsLargestComponentPhoton()){ 
				fHistoTrueClusUnConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				fHistoTrueClusUnConvGammaMCPt[fiCut]->Fill(Photon->Pt());
			}	
			if (TruePhotonCandidate->IsLargestComponentElectron()) 
				fHistoTrueClusElectronPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
				fHistoTrueClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				fHistoTrueClusConvGammaMCPt[fiCut]->Fill(((TParticle*)fMCStack->Particle(Photon->GetMother(0)))->Pt());
			}	
			if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion() && TruePhotonCandidate->IsConversionFullyContained()) 
				fHistoTrueClusConvGammaFullyPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsMerged() || TruePhotonCandidate->IsMergedPartConv() || TruePhotonCandidate->IsDalitzMerged())
				fHistoTrueClusMergedGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsMergedPartConv())
				fHistoTrueClusMergedPartConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsDalitz()) 
				fHistoTrueClusDalitzPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsDalitzMerged()) 
				fHistoTrueClusDalitzMergedPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsPhotonWithElecMother()) 
				fHistoTrueClusPhotonFromElecMotherPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsShower()) 
				fHistoTrueClusShowerPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsSubLeadingEM())
				fHistoTrueClusSubLeadingPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			fHistoTrueClusNParticles[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMotherMCLabels());
			if (!TruePhotonCandidate->IsLargestComponentPhoton())
				FillPhotonBackgroundHist(TruePhotonCandidate,pdgCodeParticle);
			if (!(TruePhotonCandidate->IsLargestComponentPhoton() || (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())) )
				FillPhotonPlusConversionBackgroundHist(TruePhotonCandidate,pdgCodeParticle);
		}
		
		Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, TruePhotonCandidate->GetCaloPhotonMCLabel(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
		
		if(isPrimary){
			// filling primary histograms
			if (TruePhotonCandidate->IsLargestComponentPhoton()){
				fHistoTruePrimaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				fHistoTruePrimaryClusGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt()); // Allways Filled
			}
			if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
				fHistoTruePrimaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				fHistoTruePrimaryClusConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt()); // Allways Filled
			}
			
		} else {
			// filling secondary histograms
			if (TruePhotonCandidate->IsLargestComponentPhoton())
				fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
				fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if(fMCStack->Particle(Photon->GetMother(0))->GetMother(0) > -1){
				if(fMCStack->Particle(fMCStack->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 3122){
					if (TruePhotonCandidate->IsLargestComponentPhoton())
						fHistoTrueSecondaryClusGammaFromXFromLambdaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
					if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
						fHistoTrueSecondaryClusConvGammaFromXFromLambdaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				}
				if(fMCStack->Particle(fMCStack->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
					if (TruePhotonCandidate->IsLargestComponentPhoton())
						fHistoTrueSecondaryClusGammaFromXFromK0sPt[fiCut]->Fill(TruePhotonCandidate->Pt());
					if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
						fHistoTrueSecondaryClusConvGammaFromXFromK0sPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				}
				if(fMCStack->Particle(fMCStack->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 221){
					if (TruePhotonCandidate->IsLargestComponentPhoton())
						fHistoTrueSecondaryClusGammaFromXFromEtasPt[fiCut]->Fill(TruePhotonCandidate->Pt());
					if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
						fHistoTrueSecondaryClusConvGammaFromXFromEtasPt[fiCut]->Fill(TruePhotonCandidate->Pt());
						
				}	
			}
		}
	}
	return;
}


//________________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessTrueClusterCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate)
{
	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();

	AliAODMCParticle *Photon = NULL;
	TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	if (fDoClusterQA > 0) fHistoTrueNLabelsInClus[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels());
	if (AODMCTrackArray){
		if (!TruePhotonCandidate->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set task will abort");
		if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0) Photon = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
			else return;
	} else {
		AliInfo("AODMCTrackArray could not be loaded");
		return;
	}

	if(Photon == NULL){
	//	cout << "no photon" << endl;
		return;
	}
	Int_t pdgCodeParticle = Photon->GetPdgCode();
	TruePhotonCandidate->SetCaloPhotonMCFlagsAOD(fInputEvent);
	
	// True Photon
	if(fIsFromMBHeader && !fIsOverlappingWithOtherHeader){
		if (TruePhotonCandidate->IsLargestComponentPhoton() || TruePhotonCandidate->IsLargestComponentElectron() )fHistoTrueClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			else fHistoTrueClusEMNonLeadingPt[fiCut]->Fill(TruePhotonCandidate->Pt());
		if (fDoClusterQA > 0){
			if (TruePhotonCandidate->IsLargestComponentPhoton()) {
				fHistoTrueClusUnConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				fHistoTrueClusUnConvGammaMCPt[fiCut]->Fill(Photon->Pt());
			}	
			if (TruePhotonCandidate->IsLargestComponentElectron()) 
				fHistoTrueClusElectronPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()) {
				fHistoTrueClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				AliAODMCParticle *Mother = (AliAODMCParticle*) AODMCTrackArray->At(Photon->GetMother());
				fHistoTrueClusConvGammaMCPt[fiCut]->Fill(Mother->Pt());
			}	
			if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion() && TruePhotonCandidate->IsConversionFullyContained()) 
				fHistoTrueClusConvGammaFullyPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsMerged() || TruePhotonCandidate->IsMergedPartConv() || TruePhotonCandidate->IsDalitzMerged())
				fHistoTrueClusMergedGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsMergedPartConv())
				fHistoTrueClusMergedPartConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsDalitz()) 
				fHistoTrueClusDalitzPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsDalitzMerged()) 
				fHistoTrueClusDalitzMergedPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsPhotonWithElecMother()) 
				fHistoTrueClusPhotonFromElecMotherPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsShower()) 
				fHistoTrueClusShowerPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsSubLeadingEM())
				fHistoTrueClusSubLeadingPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			fHistoTrueClusNParticles[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMotherMCLabels());
			
			if (!TruePhotonCandidate->IsLargestComponentPhoton())
				FillPhotonBackgroundHist(TruePhotonCandidate,pdgCodeParticle);
			if (!(TruePhotonCandidate->IsLargestComponentPhoton() || (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())) )
				FillPhotonPlusConversionBackgroundHist(TruePhotonCandidate,pdgCodeParticle);
		}
		
		Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
		if(isPrimary){
			if (TruePhotonCandidate->IsLargestComponentPhoton()){
				fHistoTruePrimaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				fHistoTruePrimaryClusGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt()); // Allways Filled
			}
			if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
				fHistoTruePrimaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				fHistoTruePrimaryClusConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt()); // Allways Filled
			}
			
		} else {
			if (TruePhotonCandidate->IsLargestComponentPhoton())
				fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
				fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother() > -1){
				if(((AliAODMCParticle*)AODMCTrackArray->At(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 3122){
					if (TruePhotonCandidate->IsLargestComponentPhoton())
						fHistoTrueSecondaryClusGammaFromXFromLambdaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
					if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
						fHistoTrueSecondaryClusConvGammaFromXFromLambdaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				}
				if(((AliAODMCParticle*)AODMCTrackArray->At(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 310){
					if (TruePhotonCandidate->IsLargestComponentPhoton())
						fHistoTrueSecondaryClusGammaFromXFromK0sPt[fiCut]->Fill(TruePhotonCandidate->Pt());
					if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
						fHistoTrueSecondaryClusConvGammaFromXFromK0sPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				}
				if(((AliAODMCParticle*)AODMCTrackArray->At(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 221){
					if (TruePhotonCandidate->IsLargestComponentPhoton())
						fHistoTrueSecondaryClusGammaFromXFromEtasPt[fiCut]->Fill(TruePhotonCandidate->Pt());
					if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
						fHistoTrueSecondaryClusConvGammaFromXFromEtasPt[fiCut]->Fill(TruePhotonCandidate->Pt());
						
				}	
			}
		}	
	}
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessAODMCParticles()
{
	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();

	TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	if (AODMCTrackArray == NULL) return;
	
	// Loop over all primary MC particle
	for(Int_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {
		
		AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
		if (!particle) continue;

		Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, particle, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
		if (!isPrimary) continue;
		
		Int_t isMCFromMBHeader = -1;
		if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
			isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent);
			if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
		}
		
		if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(particle,AODMCTrackArray)){
			fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt()); // All MC Gamma
			if(particle->GetMother() >-1){ // Meson Decay Gamma
				switch((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother())))->GetPdgCode()){
				case 111: // Pi0
					fHistoMCDecayGammaPi0Pt[fiCut]->Fill(particle->Pt());
					break;
				case 113: // Rho0
					fHistoMCDecayGammaRhoPt[fiCut]->Fill(particle->Pt());
					break;
				case 221: // Eta
					fHistoMCDecayGammaEtaPt[fiCut]->Fill(particle->Pt());
					break;
				case 223: // Omega
					fHistoMCDecayGammaOmegaPt[fiCut]->Fill(particle->Pt());
					break;
				case 331: // Eta'
					fHistoMCDecayGammaEtapPt[fiCut]->Fill(particle->Pt());
					break;
				case 333: // Phi
					fHistoMCDecayGammaPhiPt[fiCut]->Fill(particle->Pt());
					break;
				case 3212: // Sigma
					fHistoMCDecayGammaSigmaPt[fiCut]->Fill(particle->Pt());
					break;
				}
			}
		}
		// Converted MC Gamma
		if(fDoMesonAnalysis){
			if(particle->GetPdgCode() == 310 && fDoMesonQA > 0){
				Double_t mesonY = 10.;
				if(particle->E() - particle->Pz() == 0 || particle->E() + particle->Pz() == 0){
					mesonY=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
				} else {
					mesonY = 0.5*(TMath::Log((particle->E()+particle->Pz()) / (particle->E()-particle->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
				}
				Float_t weightedK0s= 1;
				if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent)){
					if (particle->Pt()>0.005){
						weightedK0s= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, 0x0, fInputEvent);
						//cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
					}
				}
				fHistoMCK0sPt[fiCut]->Fill(particle->Pt(),weightedK0s);
				fHistoMCK0sWOWeightPt[fiCut]->Fill(particle->Pt());
				fHistoMCK0sPtY[fiCut]->Fill(particle->Pt(),mesonY,weightedK0s);
			}
			if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
				->MesonIsSelectedAODMC(particle,AODMCTrackArray,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
				AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetDaughter(0)));
				AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetDaughter(1)));
				Float_t weighted= 1;
				if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent)){
					if (particle->Pt()>0.005){
						weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, 0x0, fInputEvent);
					}
				}
				
				Double_t mesonY = 10.;
				if(particle->E() - particle->Pz() == 0 || particle->E() + particle->Pz() == 0){
					mesonY=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
				} else{
					mesonY = 0.5*(TMath::Log((particle->E()+particle->Pz()) / (particle->E()-particle->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
				}
				Double_t alpha = -1;
				if (particle->GetPdgCode() == 111 || particle->GetPdgCode() == 221){
					alpha = TMath::Abs((daughter0->E() - daughter1->E()))/(daughter0->E() + daughter1->E());
				}
				
				if(particle->GetPdgCode() == 111){
					fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted); // All MC Pi0
					fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt());
					if (fDoMesonQA > 0){
						fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
						fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha);
					}
				} else if(particle->GetPdgCode() == 221){
					fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted); // All MC Eta
					fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt());
					if (fDoMesonQA > 0){
						fHistoMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
						fHistoMCEtaPtAlpha[fiCut]->Fill(particle->Pt(),alpha);
					}
				}
				
				// Check the acceptance for both gammas
				if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter0,AODMCTrackArray) &&
				((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter1,AODMCTrackArray) ){
					
					if(particle->GetPdgCode() == 111){
						fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted); // MC Pi0 with gamma in acc
					} else if(particle->GetPdgCode() == 221){
						fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted); // MC Eta with gamma in acc
					}
				}
			}
		}
	}
	
}
//________________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessMCParticles()
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
			
			if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(particle,fMCStack)){
				fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt()); // All MC Gamma			
				if(particle->GetMother(0) >-1){ // Meson Decay Gamma
					switch(fMCStack->Particle(particle->GetMother(0))->GetPdgCode()){
					case 111: // Pi0
						fHistoMCDecayGammaPi0Pt[fiCut]->Fill(particle->Pt());
						break;
					case 113: // Rho0
						fHistoMCDecayGammaRhoPt[fiCut]->Fill(particle->Pt());
						break;
					case 221: // Eta
						fHistoMCDecayGammaEtaPt[fiCut]->Fill(particle->Pt());
						break;
					case 223: // Omega
						fHistoMCDecayGammaOmegaPt[fiCut]->Fill(particle->Pt());
						break;
					case 331: // Eta'
						fHistoMCDecayGammaEtapPt[fiCut]->Fill(particle->Pt());
						break;
					case 333: // Phi
						fHistoMCDecayGammaPhiPt[fiCut]->Fill(particle->Pt());
						break;
					case 3212: // Sigma
						fHistoMCDecayGammaSigmaPt[fiCut]->Fill(particle->Pt());
						break;
					}
				}
			}
			if(fDoMesonAnalysis){
				if(particle->GetPdgCode() == 310 && fDoMesonQA > 0){
					Double_t mesonY = 10.;
					if(particle->Energy() - particle->Pz() == 0 || particle->Energy() + particle->Pz() == 0){
						mesonY=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
					} else{
						mesonY = 0.5*(TMath::Log((particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
					}
					Float_t weightedK0s= 1;
					if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent)){
						if (particle->Pt()>0.005){
							weightedK0s= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, fMCStack, fInputEvent);
							//cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
						}
					}
					if (fMCStack->IsPhysicalPrimary(i)){
						fHistoMCK0sPt[fiCut]->Fill(particle->Pt(),weightedK0s);
						fHistoMCK0sWOWeightPt[fiCut]->Fill(particle->Pt());
						fHistoMCK0sPtY[fiCut]->Fill(particle->Pt(),mesonY,weightedK0s);
					}
				}
				if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
					->MesonIsSelectedMC(particle,fMCStack,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
					TParticle* daughter0 = (TParticle*)fMCStack->Particle(particle->GetFirstDaughter());
					TParticle* daughter1 = (TParticle*)fMCStack->Particle(particle->GetLastDaughter());
					
					Float_t weighted= 1;
					if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent)){
						if (particle->Pt()>0.005){
							weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, fMCStack, fInputEvent);
							//                   if(particle->GetPdgCode() == 221){
							//                      cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
							//                   }
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
						fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted); // All MC Pi0
						fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt());
						if (fDoMesonQA > 0){
							fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
							fHistoMCPi0PtAlpha[fiCut]->Fill(particle->Pt(),alpha); // All MC Pi0
						}
					} else if(particle->GetPdgCode() == 221){
						fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted); // All MC Eta
						fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt());
						if (fDoMesonQA > 0){
							fHistoMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
							fHistoMCEtaPtAlpha[fiCut]->Fill(particle->Pt(),alpha); // All MC Pi0
						}
					}
				
					// Check the acceptance for both gammas & whether they are counted as primaries as well
					Bool_t kDaughter0IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, particle->GetFirstDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
					Bool_t kDaughter1IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, particle->GetLastDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
					
					if( kDaughter0IsPrim && kDaughter1IsPrim &&
						((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCStack) &&
						((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCStack) ){					
						if(particle->GetPdgCode() == 111){
							fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted); // MC Pi0 with gamma in acc
						} else if(particle->GetPdgCode() == 221){
							fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted); // MC Eta with gamma in acc
						}
					}
				}
			}
		} else {
			if (fDoMesonQA){
				// fill secondary histograms
				TParticle* particle = (TParticle *)fMCStack->Particle(i);
				if (!particle) continue;
		
				Int_t isMCFromMBHeader = -1;
				if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
					isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent);
					if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
				}
		
				if(fDoMesonAnalysis){
					if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMC(particle,fMCStack,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
						Float_t weighted= 1;
						if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent)){
							if (particle->Pt()>0.005){
								weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, fMCStack, fInputEvent);
				//                   if(particle->GetPdgCode() == 221){
				//                      cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
				//                   }
							}
						}
						
						if(particle->GetPdgCode() == 111){
							Int_t pdgCode = ((TParticle*)fMCStack->Particle( particle->GetFirstMother() ))->GetPdgCode();
							Int_t source = GetSourceClassification(111,pdgCode);
							fHistoMCSecPi0PtvsSource[fiCut]->Fill(particle->Pt(),source,weighted); // All MC Pi0
							fHistoMCSecPi0Source[fiCut]->Fill(pdgCode);
						} else if(particle->GetPdgCode() == 221){
							Int_t pdgCode = ((TParticle*)fMCStack->Particle( particle->GetFirstMother() ))->GetPdgCode();
							fHistoMCSecEtaPt[fiCut]->Fill(particle->Pt(),weighted); // All MC Pi0
							fHistoMCSecEtaSource[fiCut]->Fill(pdgCode);
						}
					}
				}
			}			
		}	
	}
  

}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::CalculatePi0Candidates(){
	
	// Conversion Gammas
	if(fClusterCandidates->GetEntries()>0){

		// vertex
		Double_t vertex[3] = {0};
		InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

		for(Int_t firstGammaIndex=0;firstGammaIndex<fClusterCandidates->GetEntries();firstGammaIndex++){
			AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(firstGammaIndex));
			if (gamma0==NULL) continue;
			
			for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
				AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
				if (gamma1==NULL) continue;
				
				AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);
				pi0cand->SetLabels(firstGammaIndex,secondGammaIndex);
								
				if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
					fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
					// fill new histograms
					if(abs(pi0cand->GetAlpha())<0.1)
						fHistoMotherInvMassEalpha[fiCut]->Fill(pi0cand->M(),pi0cand->E());
					
					if (fDoMesonQA > 0){
						if ( pi0cand->M() > 0.05 && pi0cand->M() < 0.17){
							fHistoMotherPi0PtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
							fHistoMotherPi0PtAlpha[fiCut]->Fill(pi0cand->Pt(),abs(pi0cand->GetAlpha()));
							fHistoMotherPi0PtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle());	
						}
						if ( pi0cand->M() > 0.45 && pi0cand->M() < 0.65){
							fHistoMotherEtaPtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
							fHistoMotherEtaPtAlpha[fiCut]->Fill(pi0cand->Pt(),abs(pi0cand->GetAlpha()));
							fHistoMotherEtaPtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle());
						}
					}
                    if(fDoTHnSparse && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoBGCalculation()){
						Int_t zbin = 0;
						Int_t mbin = 0;
						
						if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->BackgroundHandlerType() == 0){
							zbin = fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
							if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
								mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
							} else {
								mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fClusterCandidates->GetEntries());
							}
						} 
						Double_t sparesFill[4] = {pi0cand->M(),pi0cand->Pt(),(Double_t)zbin,(Double_t)mbin};
						fSparseMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
					}
				
					if(fIsMC){
						if(fInputEvent->IsA()==AliESDEvent::Class())
							ProcessTrueMesonCandidates(pi0cand,gamma0,gamma1);
						if(fInputEvent->IsA()==AliAODEvent::Class())
							ProcessTrueMesonCandidatesAOD(pi0cand,gamma0,gamma1);
					}
					
					if (fDoMesonQA == 1){
						fHistoMotherInvMassECalib[fiCut]->Fill(pi0cand->M(),gamma1->E());
						if(abs(pi0cand->GetAlpha())<0.1)
						fHistoMotherInvMassECalibalpha[fiCut]->Fill(pi0cand->M(),gamma1->E());            
					}
					
				}
				for (Int_t thirdGammaIndex=0;thirdGammaIndex<fClusterCandidates->GetEntries();thirdGammaIndex++){
					if (firstGammaIndex == thirdGammaIndex || secondGammaIndex == thirdGammaIndex ) continue;
					AliAODConversionPhoton *gamma2=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(thirdGammaIndex));
					if (gamma2==NULL) continue;
					
					AliAODConversionMother *pi0cand2 = new AliAODConversionMother(pi0cand,gamma2);
					if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand2,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
						fHistoMotherInvMass3ClusterPt[fiCut]->Fill(pi0cand2->M(),pi0cand2->Pt());
					}
					delete pi0cand2;
					pi0cand2=0x0;
					
				}
				
				delete pi0cand;
				pi0cand=0x0;
			}
		}
	}
}
//______________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
	// Process True Mesons
	AliStack *MCStack = fMCEvent->Stack();
	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();

	
	Bool_t isTruePi0 				= kFALSE;
	Bool_t isTrueEta 				= kFALSE;
	Bool_t isSameConvertedGamma 	= kFALSE;
	Int_t convertedPhotonLabel0		= -1;
	Int_t convertedPhotonLabel1		= -1;
	
	Int_t gamma0MCLabel = TrueGammaCandidate0->GetCaloPhotonMCLabel(0); 	// get most probable MC label
	Int_t gamma0MotherLabel = -1;

	if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
		TParticle * gammaMC0 = (TParticle*)MCStack->Particle(gamma0MCLabel);
		if (TrueGammaCandidate0->IsLargestComponentPhoton() || TrueGammaCandidate0->IsLargestComponentElectron()){		// largest component is electro magnetic
			// get mother of interest (pi0 or eta)
			if (TrueGammaCandidate0->IsLargestComponentPhoton()){														// for photons its the direct mother 
				gamma0MotherLabel=gammaMC0->GetMother(0);
			} else if (TrueGammaCandidate0->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
				convertedPhotonLabel0 = gammaMC0->GetMother(0);
				if (TrueGammaCandidate0->IsConversion()) gamma0MotherLabel=MCStack->Particle(gammaMC0->GetMother(0))->GetMother(0);
				else gamma0MotherLabel=gammaMC0->GetMother(0); 
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
				convertedPhotonLabel1 = gammaMC1->GetMother(0);
				if (TrueGammaCandidate1->IsConversion()) gamma1MotherLabel=MCStack->Particle(gammaMC1->GetMother(0))->GetMother(0);
				else gamma1MotherLabel=gammaMC1->GetMother(0); 
			}
		} 	
	}
			
	if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
		if(((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
			isTruePi0=kTRUE;
			if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
		}
		if(((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetPdgCode() == 221){
			isTrueEta=kTRUE;
			if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,gamma0MotherLabel)) fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
		}
	}
	
	if (convertedPhotonLabel0 > -1 && convertedPhotonLabel1 > -1){
		if (convertedPhotonLabel0==convertedPhotonLabel1) isSameConvertedGamma = kTRUE;
	}
	
	if(isTruePi0 || isTrueEta){// True Pion or Eta
		if (isTruePi0) 	fHistoTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
		if (isTrueEta) 	fHistoTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
		if (fDoMesonQA > 0){
			// both gammas are real gammas
			if (TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate1->IsLargestComponentPhoton()) {
				if (isTruePi0) fHistoTruePi0CaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta) fHistoTrueEtaCaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	
			// both particles are electrons
			if (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate1->IsLargestComponentElectron() ) {
				if (isTruePi0) fHistoTruePi0CaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta) fHistoTrueEtaCaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	
			// both particles are converted electrons
			if ((TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ){
				if (isTruePi0 )fHistoTruePi0CaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta )fHistoTrueEtaCaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
			// 1 gamma is converted the other one is normals
			if ( (TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ||
				 (TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion())
			) {
				if (isTruePi0) fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta) fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
			
			if ( (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion() && TrueGammaCandidate1->IsLargestComponentPhoton() && !TrueGammaCandidate1->IsMerged()) ||
				 (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion() && TrueGammaCandidate0->IsLargestComponentPhoton() && !TrueGammaCandidate0->IsMerged()) 
			) {
				if (isTruePi0) fHistoTruePi0NonMergedElectronPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	

			if ( (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion() && TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate1->IsMerged()) ||
				 (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion() && TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate0->IsMerged()) 
			) {
				if (isTruePi0) fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	
			
			// at least one of the photon is merged
			if (TrueGammaCandidate0->IsMerged() || TrueGammaCandidate0->IsMergedPartConv() || TrueGammaCandidate0->IsDalitzMerged() || TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate1->IsDalitzMerged() ){
				if (isTruePi0) fHistoTruePi0CaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta) fHistoTruePi0CaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	
			// at least one of the photon is merged and part conv
			if (TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate0->IsMergedPartConv()) {	
				if (isTruePi0) fHistoTruePi0CaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta) fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	
		}
	
		if (fDoMesonQA == 2){
			// category 1: 2 real photons unmerged
			if (TrueGammaCandidate0->IsLargestComponentPhoton() && !TrueGammaCandidate0->IsMerged() && TrueGammaCandidate1->IsLargestComponentPhoton() && !TrueGammaCandidate1->IsMerged()) {
				if (isTruePi0) fHistoTruePi0Category1[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta) fHistoTrueEtaCategory1[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
			// category 2, 3: 1 real photons unmerged,  1 electron (category 2 merged, category 3 unmerged )
			// -> photon 0 is unconverted
			if ( (TrueGammaCandidate0->IsLargestComponentPhoton() && !TrueGammaCandidate0->IsMerged()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion())) {
				if (isTruePi0){
					if (TrueGammaCandidate1->IsMergedPartConv())	fHistoTruePi0Category2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (!TrueGammaCandidate1->IsMergedPartConv()){
						fHistoTruePi0Category3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					}	
				}
				if (isTrueEta){
					if (TrueGammaCandidate1->IsMergedPartConv())	fHistoTrueEtaCategory2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (!TrueGammaCandidate1->IsMergedPartConv())	fHistoTrueEtaCategory3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				}	
			}
			// -> photon 1 is unconverted
			if ( ( TrueGammaCandidate1->IsLargestComponentPhoton() && !TrueGammaCandidate1->IsMerged()) && (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion())) {
				if (isTruePi0){
					if (TrueGammaCandidate0->IsMergedPartConv())	fHistoTruePi0Category2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (!TrueGammaCandidate0->IsMergedPartConv()){
						fHistoTruePi0Category3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					}	
				}
				if (isTrueEta){
					if (TrueGammaCandidate0->IsMergedPartConv())	fHistoTrueEtaCategory2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (!TrueGammaCandidate0->IsMergedPartConv())	fHistoTrueEtaCategory3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				}	
			}
			
			// category 4 & 6, 5, 7, 8
			if ( (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ){
				if (isTruePi0){
					// category 4: both electrons are from same conversion 
					if (isSameConvertedGamma && !TrueGammaCandidate0->IsMergedPartConv() && !TrueGammaCandidate1->IsMergedPartConv() ) fHistoTruePi0Category4_6[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (!isSameConvertedGamma ){
						if (!TrueGammaCandidate0->IsMergedPartConv() && !TrueGammaCandidate1->IsMergedPartConv()){ 		// category 5: both electrons from different converted photons, electrons not merged
							fHistoTruePi0Category5[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
						} else if (TrueGammaCandidate0->IsMergedPartConv() && TrueGammaCandidate1->IsMergedPartConv()){ // category 8: both electrons from different converted photons, both electrons merged		
							fHistoTruePi0Category8[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());  
						} else { 																		// category 7: both electrons from different converted photons, 1 electrons not merged, 1 electron merged		
							fHistoTruePi0Category7[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
						}	
					}	
				}
				if (isTrueEta){
					// category 4: both electrons are from same conversion 
					if (isSameConvertedGamma && !TrueGammaCandidate0->IsMergedPartConv() && !TrueGammaCandidate1->IsMergedPartConv()) fHistoTrueEtaCategory4_6[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (!isSameConvertedGamma ){
						if (!TrueGammaCandidate0->IsMergedPartConv() && !TrueGammaCandidate1->IsMergedPartConv()){ 		// category 5: both electrons from different converted photons, electrons not merged
							fHistoTrueEtaCategory5[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
						} else if (TrueGammaCandidate0->IsMergedPartConv() && TrueGammaCandidate1->IsMergedPartConv()){ // category 8: both electrons from different converted photons, both electrons merged		
							fHistoTrueEtaCategory8[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());  
						} else { 																		// category 7: both electrons from different converted photons, 1 electrons not merged, 1 electron merged		
							fHistoTrueEtaCategory7[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
						}	
					}	
				}
			}
		}
	
		if (fDoMesonQA > 0){
			if (isTruePi0){
				if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
					fHistoTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
					fHistoTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),abs(Pi0Candidate->GetAlpha()));
					fHistoTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
				}
			} else if (isTrueEta){
				if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
					fHistoTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
					fHistoTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),abs(Pi0Candidate->GetAlpha()));
					fHistoTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
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
			if (isTruePi0) fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
			if (secMotherLabel >-1){
				if(MCStack->Particle(secMotherLabel)->GetPdgCode()==310 && isTruePi0){
					fHistoTrueSecondaryPi0FromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
					if (fDoMesonQA > 0)fHistoTrueK0sWithPi0DaughterMCPt[fiCut]->Fill(MCStack->Particle(secMotherLabel)->Pt());
				}
				if(MCStack->Particle(secMotherLabel)->GetPdgCode()==221 && isTruePi0){
					fHistoTrueSecondaryPi0FromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
					if (fDoMesonQA > 0)fHistoTrueEtaWithPi0DaughterMCPt[fiCut]->Fill(MCStack->Particle(secMotherLabel)->Pt());
				}
				if(MCStack->Particle(secMotherLabel)->GetPdgCode()==3122 && isTruePi0){
					fHistoTrueSecondaryPi0FromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
					if (fDoMesonQA > 0)fHistoTrueLambdaWithPi0DaughterMCPt[fiCut]->Fill(MCStack->Particle(secMotherLabel)->Pt());
				}
			}
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
				fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
				fHistoTruePrimaryPi0W0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				fProfileTruePrimaryPi0WeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
			} else if (isTrueEta){
				fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
				fHistoTruePrimaryEtaW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				fProfileTruePrimaryEtaWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
			}	
				
			if (fDoMesonQA > 0){
				if(isTruePi0){ // Only primary pi0 for resolution
					fHistoTruePrimaryPi0MCPtResolPt[fiCut]->Fill(((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt())/((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt(),weighted);
				}
				if (isTrueEta){ // Only primary eta for resolution
					fHistoTruePrimaryEtaMCPtResolPt[fiCut]->Fill(((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt())/((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt(),weighted);
				}
			}
		}	
	} else if(!isTruePi0 && !isTrueEta){ // Background
		if (fDoMesonQA > 0){
			if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
				fHistoTrueBckGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			} else { // No photon or without mother
				fHistoTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
		}
	}

}
//______________________________________________________________________
void AliAnalysisTaskGammaCalo::ProcessTrueMesonCandidatesAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();

	// Process True Mesons
	TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	if (AODMCTrackArray == NULL) return;
	
	Bool_t isTruePi0 				= kFALSE;
	Bool_t isTrueEta 				= kFALSE;
	Bool_t isSameConvertedGamma 	= kFALSE;
	Int_t convertedPhotonLabel0		= -1;
	Int_t convertedPhotonLabel1		= -1;
		
	Int_t gamma0MCLabel 			= TrueGammaCandidate0->GetCaloPhotonMCLabel(0); 	// get most probable MC label
	Int_t gamma0MotherLabel 		= -1;
	
	// check if 
	if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
		// Daughters Gamma 0
		AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
		if (TrueGammaCandidate0->IsLargestComponentPhoton() || TrueGammaCandidate0->IsLargestComponentElectron()){		// largest component is electro magnetic
			// get mother of interest (pi0 or eta)
			if (TrueGammaCandidate0->IsLargestComponentPhoton()){														// for photons its the direct mother 
				gamma0MotherLabel=gammaMC0->GetMother();
			} else if (TrueGammaCandidate0->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
				if (TrueGammaCandidate0->IsConversion()){
					convertedPhotonLabel0 = gammaMC0->GetMother();
					AliAODMCParticle * gammaGrandMotherMC0 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC0->GetMother()));
					gamma0MotherLabel=gammaGrandMotherMC0->GetMother();
				} else gamma0MotherLabel=gammaMC0->GetMother(); 
			}
		}	
	}

	Int_t gamma1MCLabel 			= TrueGammaCandidate1->GetCaloPhotonMCLabel(0); 	// get most probable MC label
	Int_t gamma1MotherLabel 		= -1;
	
	// check if 
	if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
		// Daughters Gamma 1
		AliAODMCParticle * gammaMC1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MCLabel));
		if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){		// largest component is electro magnetic
			// get mother of interest (pi0 or eta)
			if (TrueGammaCandidate1->IsLargestComponentPhoton()){														// for photons its the direct mother 
				gamma1MotherLabel=gammaMC1->GetMother();
			} else if (TrueGammaCandidate1->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
				if (TrueGammaCandidate1->IsConversion()){
					convertedPhotonLabel1 = gammaMC1->GetMother();
					AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC1->GetMother()));
					gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
				} else gamma1MotherLabel=gammaMC1->GetMother(); 
			}
		}	
	}
			
	if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
		if(((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
			isTruePi0=kTRUE;
			if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
		}
		if(((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 221){
			isTrueEta=kTRUE;
			if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,gamma0MotherLabel)) fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
		}
	}

	if (convertedPhotonLabel0 > -1 && convertedPhotonLabel1 > 1){
		if (convertedPhotonLabel0==convertedPhotonLabel1) isSameConvertedGamma = kTRUE;
	}

	
	if(isTruePi0 || isTrueEta){// True Pion or Eta
		if (isTruePi0)fHistoTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
		if (isTrueEta)fHistoTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
		if (fDoMesonQA > 0){
			// both gammas are real gammas
			if (TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate1->IsLargestComponentPhoton()) {
				if (isTruePi0)fHistoTruePi0CaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta)fHistoTrueEtaCaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	
			// both particles are electrons
			if (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate1->IsLargestComponentElectron() ) {
				if (isTruePi0) fHistoTruePi0CaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta) fHistoTrueEtaCaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	
			// both particles are converted electrons
			if ((TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ){
				if (isTruePi0 )fHistoTruePi0CaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta )fHistoTrueEtaCaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
			// 1 gamma is converted the other one is normals
			if ( (TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ||
				 (TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion())
			) {
				if (isTruePi0) fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta) fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
			
			if ( (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion() && TrueGammaCandidate1->IsLargestComponentPhoton() && !TrueGammaCandidate1->IsMerged()) ||
				 (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion() && TrueGammaCandidate0->IsLargestComponentPhoton() && !TrueGammaCandidate0->IsMerged()) 
			) {
				if (isTruePi0) fHistoTruePi0NonMergedElectronPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	

			if ( (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion() && TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate1->IsMerged()) ||
				 (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion() && TrueGammaCandidate0->IsLargestComponentPhoton() && TrueGammaCandidate0->IsMerged()) 
			) {
				if (isTruePi0) fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	

			// at least one of the photon is merged
			if (TrueGammaCandidate0->IsMerged() || TrueGammaCandidate0->IsMergedPartConv() || TrueGammaCandidate0->IsDalitzMerged() || TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate1->IsDalitzMerged() ){
				if (isTruePi0) fHistoTruePi0CaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta) fHistoTrueEtaCaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
			// at least one of the photon is merged and part conv
			if (TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate0->IsMergedPartConv()) {
				if (isTruePi0)fHistoTruePi0CaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta)fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	
		}

		if (fDoMesonQA == 2){
			// category 1: 2 real photons unmerged
			if (TrueGammaCandidate0->IsLargestComponentPhoton() && !TrueGammaCandidate0->IsMerged() && TrueGammaCandidate1->IsLargestComponentPhoton() && !TrueGammaCandidate1->IsMerged()) {
				if (isTruePi0) fHistoTruePi0Category1[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta) fHistoTrueEtaCategory1[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
			// category 2, 3: 1 real photons unmerged,  1 electron (category 2 merged, category 3 unmerged )
			// -> photon 0 is unconverted
			if ( (TrueGammaCandidate0->IsLargestComponentPhoton() && !TrueGammaCandidate0->IsMerged()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion())) {
				if (isTruePi0){
					if (TrueGammaCandidate1->IsMergedPartConv())	fHistoTruePi0Category2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (!TrueGammaCandidate1->IsMergedPartConv())	fHistoTruePi0Category3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				}
				if (isTrueEta){
					if (TrueGammaCandidate1->IsMergedPartConv())	fHistoTrueEtaCategory2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (!TrueGammaCandidate1->IsMergedPartConv())	fHistoTrueEtaCategory3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				}	
			}
			// -> photon 1 is unconverted
			if ( ( TrueGammaCandidate1->IsLargestComponentPhoton() && !TrueGammaCandidate1->IsMerged()) && (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion())) {
				if (isTruePi0){
					if (TrueGammaCandidate0->IsMergedPartConv())	fHistoTruePi0Category2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (!TrueGammaCandidate0->IsMergedPartConv())	fHistoTruePi0Category3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				}
				if (isTrueEta){
					if (TrueGammaCandidate0->IsMergedPartConv())	fHistoTrueEtaCategory2[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (!TrueGammaCandidate0->IsMergedPartConv())	fHistoTrueEtaCategory3[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				}	
			}
			
			// category 4 & 6, 5, 7, 8
			if ( (TrueGammaCandidate0->IsLargestComponentElectron() && TrueGammaCandidate0->IsConversion()) && (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()) ){
				if (isTruePi0){
					// category 4: both electrons are from same conversion 
					if (isSameConvertedGamma && !TrueGammaCandidate0->IsMergedPartConv() && !TrueGammaCandidate1->IsMergedPartConv() ) fHistoTruePi0Category4_6[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (!isSameConvertedGamma ){
						if (!TrueGammaCandidate0->IsMergedPartConv() && !TrueGammaCandidate1->IsMerged()){ 		// category 5: both electrons from different converted photons, electrons not merged
							fHistoTruePi0Category5[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
						} else if (TrueGammaCandidate0->IsMergedPartConv() && TrueGammaCandidate1->IsMerged()){ // category 8: both electrons from different converted photons, both electrons merged		
							fHistoTruePi0Category8[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());  
						} else { 																		// category 7: both electrons from different converted photons, 1 electrons not merged, 1 electron merged		
							fHistoTruePi0Category7[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
						}	
					}	
				}
				if (isTrueEta){
					// category 4: both electrons are from same conversion 
					if (isSameConvertedGamma && !TrueGammaCandidate0->IsMergedPartConv() && !TrueGammaCandidate1->IsMergedPartConv()) fHistoTrueEtaCategory4_6[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (!isSameConvertedGamma ){
						if (!TrueGammaCandidate0->IsMergedPartConv() && !TrueGammaCandidate1->IsMergedPartConv()){ 		// category 5: both electrons from different converted photons, electrons not merged
							fHistoTrueEtaCategory5[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
						} else if (TrueGammaCandidate0->IsMergedPartConv() && TrueGammaCandidate1->IsMergedPartConv()){ // category 8: both electrons from different converted photons, both electrons merged		
							fHistoTrueEtaCategory8[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());  
						} else { 																		// category 7: both electrons from different converted photons, 1 electrons not merged, 1 electron merged		
							fHistoTrueEtaCategory7[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
						}	
					}	
				}
			}
		}
		
		if (fDoMesonQA > 0){
			if (isTruePi0){
				if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
				fHistoTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
				fHistoTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),abs(Pi0Candidate->GetAlpha()));
				fHistoTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
				}
			} else if (isTrueEta){
				if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
				fHistoTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
				fHistoTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),abs(Pi0Candidate->GetAlpha()));
				fHistoTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
				}
			}
		}

		Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel)), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
		if(!isPrimary){ // Secondary Meson
			Int_t secMotherLabel = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->GetMother();
			Float_t weightedSec= 1;
			if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, 0x0, fInputEvent) && static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
				weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),secMotherLabel, 0x0, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
				//cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
			}
			if (isTruePi0) fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
			if (secMotherLabel >-1){
				if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310 && isTruePi0 ){
					fHistoTrueSecondaryPi0FromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
					if (fDoMesonQA > 0)fHistoTrueK0sWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->Pt());
				}
				if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==221 && isTruePi0){
					fHistoTrueSecondaryPi0FromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
					if (fDoMesonQA > 0)fHistoTrueEtaWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->Pt());
				}
				if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==3122 && isTruePi0){
					fHistoTrueSecondaryPi0FromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
					if (fDoMesonQA > 0)fHistoTrueLambdaWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->Pt());
				}
			}	
		} else{ // Only primary pi0 for efficiency calculation
			Float_t weighted= 1;
			if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, 0x0, fInputEvent)){
				if (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt()>0.005){
				weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),gamma1MotherLabel, 0x0, fInputEvent);
				//                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
				}
			}
			if (isTruePi0){
				fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
				fHistoTruePrimaryPi0W0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				fProfileTruePrimaryPi0WeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
			} else if (isTrueEta){
				fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
				fHistoTruePrimaryEtaW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				fProfileTruePrimaryEtaWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);	
			}	
			if (fDoMesonQA > 0){
				if(isTruePi0){ // Only primary pi0 for resolution
					fHistoTruePrimaryPi0MCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),
															(Pi0Candidate->Pt()-static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted);
				
				}
				if (isTrueEta){ // Only primary eta for resolution
					fHistoTruePrimaryEtaMCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),
															(Pi0Candidate->Pt()-static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted);
				}
			}
		}
	} else if(!isTruePi0 && !isTrueEta) { // Background
		if (fDoMesonQA > 0){
			if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
				fHistoTrueBckGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			} else { // No photon or without mother
				fHistoTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
		}
	}
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::CalculateBackground(){
	
	Int_t zbin= fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
	Int_t mbin = 0;
	
	if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
		mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
	} else {
		mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fClusterCandidates->GetEntries());
	}
	
	if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
		for(Int_t nEventsInBG=0;nEventsInBG<fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
			AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
			for(Int_t iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
				AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
				for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
					AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
					AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
					backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
				
					if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
						->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
						fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
                        if(fDoTHnSparse){
                            Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
                            fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
                        }
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
				for(Int_t iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
					AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
					for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
				
						AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
						AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
						backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
				
						if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
							fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
                            if(fDoTHnSparse){
                                Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
                                fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
                            }
						}
						delete backgroundCandidate;
						backgroundCandidate = 0x0;
					}
				}
			}
		}
	}
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::RotateParticle(AliAODConversionPhoton *gamma){
	Int_t fNDegreesPMBackground= ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->NDegreesRotation();
	Double_t nRadiansPM = fNDegreesPMBackground*TMath::Pi()/180;
	Double_t rotationValue = fRandom.Rndm()*2*nRadiansPM + TMath::Pi()-nRadiansPM;
	gamma->RotateZ(rotationValue);
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::UpdateEventByEventData(){
	//see header file for documentation
	if(fClusterCandidates->GetEntries() >0 ){
		if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
			fBGHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
		} else { // means we use #V0s for multiplicity
			fBGHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fClusterCandidates->GetEntries(),fEventPlaneAngle);
		}
	}
}


//________________________________________________________________________
void AliAnalysisTaskGammaCalo::FillPhotonBackgroundHist(AliAODConversionPhoton *TruePhotonCandidate, Int_t pdgCode)
{
	// Bck = 0 e+-, 1 pi+-, 2 p+-, 3 K+-, 4 n, 5 K0s, 6 Lambda, 7 mu+-, 8 rest
	if(fIsFromMBHeader){
		if(abs(pdgCode) == 11) 			fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0);
		else if( abs(pdgCode)==211) 	fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1);
		else if( abs(pdgCode)==2212) 	fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2);
		else if( abs(pdgCode)==321) 	fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3);
		else if( abs(pdgCode)==2112) 	fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4);
		else if( abs(pdgCode)==310) 	fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),5);   
		else if( abs(pdgCode)==3122) 	fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),6);
		else if( abs(pdgCode)==13) 		fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),7);
		else 							fHistoClusPhotonBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),8);
	}	
}

//________________________________________________________________________
void AliAnalysisTaskGammaCalo::FillPhotonPlusConversionBackgroundHist(AliAODConversionPhoton *TruePhotonCandidate, Int_t pdgCode)
{
	// Bck = 0 e+-, 1 pi+-, 2 p+-, 3 K+-, 4 n, 5 K0s, 6 Lambda, 7 mu+-, 8 rest
	if(fIsFromMBHeader){
		if(abs(pdgCode) == 11) 			fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0);
		else if( abs(pdgCode)==211) 	fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1);
		else if( abs(pdgCode)==2212) 	fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2);
		else if( abs(pdgCode)==321) 	fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3);
		else if( abs(pdgCode)==2112) 	fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4);
		else if( abs(pdgCode)==310) 	fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),5);   
		else if( abs(pdgCode)==3122) 	fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),6);
		else if( abs(pdgCode)==13) 		fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),7);
		else 							fHistoClusPhotonPlusConvBGPt[fiCut]->Fill(TruePhotonCandidate->Pt(),8);
	}	
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaCalo::SetLogBinningXTH2(TH2* histoRebin){
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
Bool_t AliAnalysisTaskGammaCalo::CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked)
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
void AliAnalysisTaskGammaCalo::Terminate(const Option_t *)
{
  
  //fOutputContainer->Print(); // Will crash on GRID
}

//________________________________________________________________________
Int_t AliAnalysisTaskGammaCalo::GetSourceClassification(Int_t daughter, Int_t pdgCode){
  
	if (daughter == 111) {
		if (abs(pdgCode) == 310) return 1; // k0s
		else if (abs(pdgCode) == 3122) return 2; // Lambda
		else if (abs(pdgCode) == 130) return 3; // K0L
		else if (abs(pdgCode) == 2212) return 4; // proton
		else if (abs(pdgCode) == 2112) return 5; // neutron
		else if (abs(pdgCode) == 211) return 6; // pion
		else if (abs(pdgCode) == 321) return 7; // kaon
		else if (abs(pdgCode) == 113 || abs(pdgCode) == 213 ) return 8; // rho 0,+,-
		else if (abs(pdgCode) == 3222 || abs(pdgCode) == 3212 || abs(pdgCode) == 3112  ) return 9; // Sigma
		else if (abs(pdgCode) == 2224 || abs(pdgCode) == 2214 || abs(pdgCode) == 2114 || abs(pdgCode) == 1114  ) return 10; // Delta
		else if (abs(pdgCode) == 313 || abs(pdgCode) == 323   ) return 11; // K*
		else return 15;
	}
	return 15;
  
}
