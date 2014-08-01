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
#include "AliAnalysisTaskGammaConvCalo.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliKFVertex.h"
#include "AliV0ReaderV1.h"
#include "AliGenCocktailEventHeader.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliEventplane.h"
#include "AliAnalysisTaskEMCALClusterizeFast.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"

ClassImp(AliAnalysisTaskGammaConvCalo)

//________________________________________________________________________
AliAnalysisTaskGammaConvCalo::AliAnalysisTaskGammaConvCalo(): AliAnalysisTaskSE(),
	fV0Reader(NULL),
	fBGHandler(NULL),
	fBGHandlerRP(NULL),
	fBGClusHandler(NULL),
	fBGClusHandlerRP(NULL),
	fInputEvent(NULL),
	fMCEvent(NULL),
	fMCStack(NULL),
	fCutFolder(NULL),
	fESDList(NULL),
	fBackList(NULL),
	fMotherList(NULL),
	fPhotonDCAList(NULL),
	fMesonDCAList(NULL),
	fTrueList(NULL),
	fMCList(NULL),
	fHeaderNameList(NULL),
	fTagOutputList(NULL),
	fOutputContainer(NULL),
	fReaderGammas(NULL),
	fGammaCandidates(NULL),
	fClusterCandidates(NULL),
	fEventCutArray(NULL),
	fEventCuts(NULL),
	fCutArray(NULL),
	fConversionCuts(NULL),
	fClusterCutArray(NULL),
	fCaloPhotonCuts(NULL),
	fMesonCutArray(NULL),
	fMesonCuts(NULL),
	fHistoConvGammaPt(NULL),
	fHistoConvGammaR(NULL),
	fHistoConvGammaEta(NULL),
	fTreeConvGammaPtDcazCat(NULL),
	fPtGamma(0),
	fDCAzPhoton(0),
	fRConvPhoton(0),
	fEtaPhoton(0),
	fCharCatPhoton(0),
	fCharPhotonMCInfo(0),
	fHistoMotherInvMassPt(NULL),
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
	fTreeMesonsInvMassPtDcazMinDcazMaxFlag(NULL),
	fInvMass(0),
	fPt(0),
	fDCAzGammaMin(0),
	fDCAzGammaMax(0),
	fCharFlag(0),
	fCharMesonMCInfo(0),
	fHistoConvGammaUntagged(NULL),
	fHistoConvGammaTagged(NULL),
	fHistoConvGammaPi0Tagged(NULL),
	fHistoConvGammaEtaTagged(NULL),
	fHistoPhotonPairAll(NULL),
	fHistoPhotonPairAllGam(NULL),
	fHistoClusGammaPt(NULL),
	fHistoMCHeaders(NULL),
	fHistoMCAllGammaPt(NULL),
	fHistoMCAllGammaEMCALAccPt(NULL),
	fHistoMCDecayGammaPi0Pt(NULL),
	fHistoMCDecayGammaRhoPt(NULL),
	fHistoMCDecayGammaEtaPt(NULL),
	fHistoMCDecayGammaOmegaPt(NULL),
	fHistoMCDecayGammaEtapPt(NULL),
	fHistoMCDecayGammaPhiPt(NULL),
	fHistoMCDecayGammaSigmaPt(NULL),
	fHistoMCConvGammaPt(NULL),
	fHistoMCConvGammaR(NULL),
	fHistoMCConvGammaEta(NULL),
	fHistoMCPi0Pt(NULL),
	fHistoMCPi0WOWeightPt(NULL),
	fHistoMCEtaPt(NULL),
	fHistoMCEtaWOWeightPt(NULL),
	fHistoMCPi0InAccPt(NULL),
	fHistoMCEtaInAccPt(NULL),
	fHistoMCPi0PtY(NULL),
	fHistoMCEtaPtY(NULL),
	fHistoMCK0sPt(NULL),
	fHistoMCK0sWOWeightPt(NULL),
	fHistoMCK0sPtY(NULL),
	fHistoMCSecPi0PtvsSource(NULL),
	fHistoMCSecPi0Source(NULL),
	fHistoMCSecEtaPt(NULL),
	fHistoMCSecEtaSource(NULL),
	fHistoTrueMotherInvMassPt(NULL),
	fHistoTrueMotherCaloPhotonInvMassPt(NULL),
	fHistoTrueMotherCaloConvertedPhotonInvMassPt(NULL),
	fHistoTruePi0CaloConvertedPhotonInvMassPt(NULL),
	fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt(NULL),
	fHistoTrueEtaCaloConvertedPhotonInvMassPt(NULL),
	fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt(NULL),
	fHistoTrueMotherCaloElectronInvMassPt(NULL),
	fHistoTrueMotherCaloMergedClusterInvMassPt(NULL),
	fHistoTrueMotherCaloEMNonLeadingInvMassPt(NULL),
	fHistoTrueMotherCaloMergedClusterPartConvInvMassPt(NULL),
	fHistoTruePrimaryMotherInvMassPt(NULL),
	fHistoTruePrimaryMotherW0WeightingInvMassPt(NULL),
	fProfileTruePrimaryMotherWeightsInvMassPt(NULL),
	fHistoTruePrimaryPi0MCPtResolPt(NULL),
	fHistoTruePrimaryEtaMCPtResolPt(NULL),
	fHistoTrueSecondaryMotherInvMassPt(NULL),
	fHistoTrueSecondaryMotherFromK0sInvMassPt(NULL),
	fHistoTrueK0sWithPi0DaughterMCPt(NULL),
	fHistoTrueSecondaryMotherFromEtaInvMassPt(NULL),
	fHistoTrueEtaWithPi0DaughterMCPt(NULL),
	fHistoTrueSecondaryMotherFromLambdaInvMassPt(NULL),
	fHistoTrueLambdaWithPi0DaughterMCPt(NULL),
	fHistoTrueBckGGInvMassPt(NULL),
	fHistoTrueBckContInvMassPt(NULL),
	fHistoTruePi0PtY(NULL),
	fHistoTrueEtaPtY(NULL),
	fHistoTruePi0PtAlpha(NULL),
	fHistoTrueEtaPtAlpha(NULL),
	fHistoTruePi0PtOpenAngle(NULL),
	fHistoTrueEtaPtOpenAngle(NULL),
	fHistoTrueConvGammaPt(NULL),
	fHistoTrueConvPi0GammaPt(NULL),
	fHistoTrueConvGammaEta(NULL),
	fHistoCombinatorialPt(NULL),
	fHistoTruePrimaryConvGammaPt(NULL),
	fHistoTruePrimaryConvGammaESDPtMCPt(NULL),
	fHistoTrueSecondaryConvGammaPt(NULL),
	fHistoTrueSecondaryConvGammaFromXFromK0sPt(NULL),
	fHistoTrueSecondaryConvGammaFromXFromLambdaPt(NULL),
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
	fHistoNEvents(NULL),
	fHistoNGoodESDTracks(NULL),
	fHistoNGammaCandidates(NULL),
	fHistoNGoodESDTracksVsNGammaCanditates(NULL),
	fHistoNV0Tracks(NULL),
	fProfileEtaShift(NULL),
	fEventPlaneAngle(-100),
	fRandom(0),
	fNGammaCandidates(0),
	fUnsmearedPx(NULL),
	fUnsmearedPy(NULL),
	fUnsmearedPz(NULL),
	fUnsmearedE(NULL),
	fMCStackPos(NULL),
	fMCStackNeg(NULL),
	fESDArrayPos(NULL),
	fESDArrayNeg(NULL),
	fnCuts(0),
	fiCut(0),
	fMoveParticleAccordingToVertex(kTRUE),
	fIsHeavyIon(0),
	fDoMesonAnalysis(kTRUE),
	fDoMesonQA(0),
	fDoPhotonQA(0),
	fDoClusterQA(0),
	fIsFromMBHeader(kTRUE),
	fIsMC(kFALSE),
	fMinE(0.1),
	fNminCells(2),
	fEMCm02cut(0.5)
{
  
}

//________________________________________________________________________
AliAnalysisTaskGammaConvCalo::AliAnalysisTaskGammaConvCalo(const char *name):
	AliAnalysisTaskSE(name),
	fV0Reader(NULL),
	fBGHandler(NULL),
	fBGHandlerRP(NULL),
	fBGClusHandler(NULL),
	fBGClusHandlerRP(NULL),
	fInputEvent(NULL),
	fMCEvent(NULL),
	fMCStack(NULL),
	fCutFolder(NULL),
	fESDList(NULL),
	fBackList(NULL),
	fMotherList(NULL),
	fPhotonDCAList(NULL),
	fMesonDCAList(NULL),
	fTrueList(NULL),
	fMCList(NULL),
	fHeaderNameList(NULL),
	fTagOutputList(NULL),
	fOutputContainer(0),
	fReaderGammas(NULL),
	fGammaCandidates(NULL),
	fClusterCandidates(NULL),
	fEventCutArray(NULL),
	fEventCuts(NULL),
	fCutArray(NULL),
	fConversionCuts(NULL),
	fClusterCutArray(NULL),
	fCaloPhotonCuts(NULL),
	fMesonCutArray(NULL),
	fMesonCuts(NULL),
	fHistoConvGammaPt(NULL),
	fHistoConvGammaR(NULL),
	fHistoConvGammaEta(NULL),
	fTreeConvGammaPtDcazCat(NULL),
	fPtGamma(0),
	fDCAzPhoton(0),
	fRConvPhoton(0),
	fEtaPhoton(0),
	fCharCatPhoton(0),
	fCharPhotonMCInfo(0),
	fHistoMotherInvMassPt(NULL),
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
	fTreeMesonsInvMassPtDcazMinDcazMaxFlag(NULL),
	fInvMass(0),
	fPt(0),
	fDCAzGammaMin(0),
	fDCAzGammaMax(0),
	fCharFlag(0),
	fCharMesonMCInfo(0),
	fHistoConvGammaUntagged(NULL),
	fHistoConvGammaTagged(NULL),
	fHistoConvGammaPi0Tagged(NULL),
	fHistoConvGammaEtaTagged(NULL),
	fHistoPhotonPairAll(NULL),
	fHistoPhotonPairAllGam(NULL),
	fHistoClusGammaPt(NULL),
	fHistoMCHeaders(NULL),
	fHistoMCAllGammaPt(NULL),
	fHistoMCAllGammaEMCALAccPt(NULL),
	fHistoMCDecayGammaPi0Pt(NULL),
	fHistoMCDecayGammaRhoPt(NULL),
	fHistoMCDecayGammaEtaPt(NULL),
	fHistoMCDecayGammaOmegaPt(NULL),
	fHistoMCDecayGammaEtapPt(NULL),
	fHistoMCDecayGammaPhiPt(NULL),
	fHistoMCDecayGammaSigmaPt(NULL),
	fHistoMCConvGammaPt(NULL),
	fHistoMCConvGammaR(NULL),
	fHistoMCConvGammaEta(NULL),
	fHistoMCPi0Pt(NULL),
	fHistoMCPi0WOWeightPt(NULL),
	fHistoMCEtaPt(NULL),
	fHistoMCEtaWOWeightPt(NULL),
	fHistoMCPi0InAccPt(NULL),
	fHistoMCEtaInAccPt(NULL),
	fHistoMCPi0PtY(NULL),
	fHistoMCEtaPtY(NULL),
	fHistoMCK0sPt(NULL),
	fHistoMCK0sWOWeightPt(NULL),
	fHistoMCK0sPtY(NULL),
	fHistoMCSecPi0PtvsSource(NULL),
	fHistoMCSecPi0Source(NULL),
	fHistoMCSecEtaPt(NULL),
	fHistoMCSecEtaSource(NULL),
	fHistoTrueMotherInvMassPt(NULL),
	fHistoTrueMotherCaloPhotonInvMassPt(NULL),
	fHistoTrueMotherCaloConvertedPhotonInvMassPt(NULL),
	fHistoTruePi0CaloConvertedPhotonInvMassPt(NULL),
	fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt(NULL),
	fHistoTrueEtaCaloConvertedPhotonInvMassPt(NULL),
	fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt(NULL),
	fHistoTrueMotherCaloElectronInvMassPt(NULL),
	fHistoTrueMotherCaloMergedClusterInvMassPt(NULL),
	fHistoTrueMotherCaloEMNonLeadingInvMassPt(NULL),
	fHistoTrueMotherCaloMergedClusterPartConvInvMassPt(NULL),
	fHistoTruePrimaryMotherInvMassPt(NULL),
	fHistoTruePrimaryMotherW0WeightingInvMassPt(NULL),
	fProfileTruePrimaryMotherWeightsInvMassPt(NULL),
	fHistoTruePrimaryPi0MCPtResolPt(NULL),
	fHistoTruePrimaryEtaMCPtResolPt(NULL),
	fHistoTrueSecondaryMotherInvMassPt(NULL),
	fHistoTrueSecondaryMotherFromK0sInvMassPt(NULL),
	fHistoTrueK0sWithPi0DaughterMCPt(NULL),
	fHistoTrueSecondaryMotherFromEtaInvMassPt(NULL),
	fHistoTrueEtaWithPi0DaughterMCPt(NULL),
	fHistoTrueSecondaryMotherFromLambdaInvMassPt(NULL),
	fHistoTrueLambdaWithPi0DaughterMCPt(NULL),
	fHistoTrueBckGGInvMassPt(NULL),
	fHistoTrueBckContInvMassPt(NULL),
	fHistoTruePi0PtY(NULL),
	fHistoTrueEtaPtY(NULL),
	fHistoTruePi0PtAlpha(NULL),
	fHistoTrueEtaPtAlpha(NULL),
	fHistoTruePi0PtOpenAngle(NULL),
	fHistoTrueEtaPtOpenAngle(NULL),
	fHistoTrueConvGammaPt(NULL),
	fHistoTrueConvPi0GammaPt(NULL),
	fHistoTrueConvGammaEta(NULL),
	fHistoCombinatorialPt(NULL),
	fHistoTruePrimaryConvGammaPt(NULL),
	fHistoTruePrimaryConvGammaESDPtMCPt(NULL),
	fHistoTrueSecondaryConvGammaPt(NULL),
	fHistoTrueSecondaryConvGammaFromXFromK0sPt(NULL),
	fHistoTrueSecondaryConvGammaFromXFromLambdaPt(NULL),
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
	fHistoNEvents(NULL),
	fHistoNGoodESDTracks(NULL),
	fHistoNGammaCandidates(NULL),
	fHistoNGoodESDTracksVsNGammaCanditates(NULL),
	fHistoNV0Tracks(NULL),
	fProfileEtaShift(NULL),
	fEventPlaneAngle(-100),
	fRandom(0),
	fNGammaCandidates(0),
	fUnsmearedPx(NULL),
	fUnsmearedPy(NULL),
	fUnsmearedPz(NULL),
	fUnsmearedE(NULL),
	fMCStackPos(NULL),
	fMCStackNeg(NULL),
	fESDArrayPos(NULL),
	fESDArrayNeg(NULL),
	fnCuts(0),
	fiCut(0),
	fMoveParticleAccordingToVertex(kTRUE),
	fIsHeavyIon(0),
	fDoMesonAnalysis(kTRUE),
	fDoMesonQA(0),
	fDoPhotonQA(0),
	fDoClusterQA(0),
	fIsFromMBHeader(kTRUE),
	fIsMC(kFALSE),
	fMinE(0.1),
	fNminCells(2),
	fEMCm02cut(0.5)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaConvCalo::~AliAnalysisTaskGammaConvCalo()
{
	if(fGammaCandidates){
		delete fGammaCandidates;
		fGammaCandidates = 0x0;
	}
	if(fClusterCandidates){
		delete fClusterCandidates;
		fClusterCandidates = 0x0;
	}
	if(fBGHandler){
		delete[] fBGHandler;
		fBGHandler = 0x0;
	}
	if(fBGHandlerRP){
		delete[] fBGHandlerRP;
		fBGHandlerRP = 0x0;
	}
	if(fBGClusHandler){
		delete[] fBGClusHandler;
		fBGClusHandler = 0x0;
	}
	if(fBGClusHandlerRP){
		delete[] fBGClusHandlerRP;
		fBGClusHandlerRP = 0x0;
	}
}
//___________________________________________________________
void AliAnalysisTaskGammaConvCalo::InitBack(){
	
	const Int_t nDim = 4;
	Int_t nBins[nDim] = {800,250,7,4};
	Double_t xMin[nDim] = {0,0, 0,0};
	Double_t xMax[nDim] = {0.8,25,7,4};
	
	fSparseMotherInvMassPtZM = new THnSparseF*[fnCuts];
	fSparseMotherBackInvMassPtZM = new THnSparseF*[fnCuts];
	
	fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];
	fBGHandlerRP = new AliConversionAODBGHandlerRP*[fnCuts];

	fBGClusHandler = new AliGammaConversionAODBGHandler*[fnCuts];
	fBGClusHandlerRP = new AliConversionAODBGHandlerRP*[fnCuts];
	
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		if (((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
			TString cutstringEvent 	= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
			TString cutstringPhoton	= ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
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
			
			fBackList[iCut] = new TList();
			fBackList[iCut]->SetName(Form("%s_%s_%s_%s Back histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(), cutstringMeson.Data()));
			fBackList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fBackList[iCut]);
			
			fSparseMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_m","Back_Back_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
			fBackList[iCut]->Add(fSparseMotherBackInvMassPtZM[iCut]);
			
			fMotherList[iCut] = new TList();
			fMotherList[iCut]->SetName(Form("%s_%s_%s_%s Mother histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
			fMotherList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fMotherList[iCut]);
			
			fSparseMotherInvMassPtZM[iCut] = new THnSparseF("Back_Mother_InvMass_Pt_z_m","Back_Mother_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
			fMotherList[iCut]->Add(fSparseMotherInvMassPtZM[iCut]);
			
			if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
				fBGHandler[iCut] = new AliGammaConversionAODBGHandler(
																	collisionSystem,centMin,centMax,
																	((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
																	((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity());
				fBGClusHandler[iCut] = new AliGammaConversionAODBGHandler(
																	collisionSystem,centMin,centMax,
																	((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
																	((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity());
				fBGHandlerRP[iCut] = NULL;
			} else{
				fBGHandlerRP[iCut] = new AliConversionAODBGHandlerRP(
																	((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsHeavyIon(),
																	((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity(),
																	((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents());
				fBGClusHandlerRP[iCut] = new AliConversionAODBGHandlerRP(
																	((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsHeavyIon(),
																	((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity(),
																	((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents());
				fBGHandler[iCut] = NULL;
			}
		}
	}
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::UserCreateOutputObjects(){
  
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
	fGammaCandidates = new TList();
	fClusterCandidates = new TList();
  
	fCutFolder = new TList*[fnCuts];
	fESDList = new TList*[fnCuts];
	fBackList = new TList*[fnCuts];
	fMotherList = new TList*[fnCuts];
	fHistoNEvents = new TH1I*[fnCuts];
	fHistoNGoodESDTracks = new TH1I*[fnCuts];
	fHistoNGammaCandidates = new TH1I*[fnCuts];
	fHistoNGoodESDTracksVsNGammaCanditates = new TH2F*[fnCuts];
	fHistoNV0Tracks = new TH1I*[fnCuts];
	fProfileEtaShift = new TProfile*[fnCuts];
	fHistoConvGammaPt = new TH1F*[fnCuts];
  
	if (fDoPhotonQA == 2){
		fPhotonDCAList = new TList*[fnCuts];
		fTreeConvGammaPtDcazCat = new TTree*[fnCuts];
	}
	if (fDoPhotonQA > 0){
		fHistoConvGammaR = new TH1F*[fnCuts];
		fHistoConvGammaEta = new TH1F*[fnCuts];
	}
	
	if(fDoMesonAnalysis){
		fHistoMotherInvMassPt = new TH2F*[fnCuts];
		fHistoMotherBackInvMassPt = new TH2F*[fnCuts];
		fHistoMotherInvMassEalpha = new TH2F*[fnCuts];
		if (fDoMesonQA == 2){
			fMesonDCAList = new TList*[fnCuts];
			fTreeMesonsInvMassPtDcazMinDcazMaxFlag = new TTree*[fnCuts];
		}
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
	fTagOutputList = new TList*[fnCuts];
	
	fHistoConvGammaUntagged = new TH1F*[fnCuts];
	fHistoConvGammaTagged = new TH1F*[fnCuts];
	fHistoConvGammaPi0Tagged = new TH1F*[fnCuts];
	fHistoConvGammaEtaTagged = new TH1F*[fnCuts];
	fHistoPhotonPairAll = new TH2F*[fnCuts];
	fHistoPhotonPairAllGam = new TH2F*[fnCuts];
	
	fHistoClusGammaPt = new TH1F*[fnCuts];
	
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		TString cutstringEvent 	= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
		TString cutstringPhoton = ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
		TString cutstringCalo 	= ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
		TString cutstringMeson 	= "NoMesonCut";
		if(fDoMesonAnalysis)cutstringMeson = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
    
		fCutFolder[iCut] = new TList();
		fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s_%s",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
		fCutFolder[iCut]->SetOwner(kTRUE);
		fOutputContainer->Add(fCutFolder[iCut]);
		fESDList[iCut] = new TList();
		fESDList[iCut]->SetName(Form("%s_%s_%s_%s ESD histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
		fESDList[iCut]->SetOwner(kTRUE);
		fCutFolder[iCut]->Add(fESDList[iCut]);
    
		fHistoNEvents[iCut] = new TH1I("NEvents","NEvents",9,-0.5,8.5);
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Missing MC");
		if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() == 4 ){
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
		fESDList[iCut]->Add(fHistoNEvents[iCut]);
		
		if(fIsHeavyIon == 1) fHistoNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",4000,0,4000);
		else if(fIsHeavyIon == 2) fHistoNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",400,0,400);
		else fHistoNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",200,0,200);
		fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);
		if(fIsHeavyIon == 1) fHistoNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",100,0,100);
		else if(fIsHeavyIon == 2) fHistoNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",50,0,50);
		else fHistoNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",50,0,50);
		fESDList[iCut]->Add(fHistoNGammaCandidates[iCut]);
		if(fIsHeavyIon == 1) fHistoNGoodESDTracksVsNGammaCanditates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",4000,0,4000,100,0,100);
		else if(fIsHeavyIon == 2) fHistoNGoodESDTracksVsNGammaCanditates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",400,0,400,50,0,50);
		else fHistoNGoodESDTracksVsNGammaCanditates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",200,0,200,50,0,50);
		fESDList[iCut]->Add(fHistoNGoodESDTracksVsNGammaCanditates[iCut]);
    
		
		if(fIsHeavyIon == 1) fHistoNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",30000,0,30000);
		else if(fIsHeavyIon == 2) fHistoNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",2500,0,2500);
		else fHistoNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",1500,0,1500);
		fESDList[iCut]->Add(fHistoNV0Tracks[iCut]);
		fProfileEtaShift[iCut] = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
		fESDList[iCut]->Add(fProfileEtaShift[iCut]);
		fHistoConvGammaPt[iCut] = new TH1F("ESD_ConvGamma_Pt","ESD_ConvGamma_Pt",250,0,25);
		fESDList[iCut]->Add(fHistoConvGammaPt[iCut]);
    
		if (fDoPhotonQA == 2){
			fPhotonDCAList[iCut] = new TList();
			fPhotonDCAList[iCut]->SetName(Form("%s_%s_%s_%s Photon DCA tree",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
			fPhotonDCAList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fPhotonDCAList[iCut]);
      
			fTreeConvGammaPtDcazCat[iCut] = new TTree("ESD_ConvGamma_Pt_Dcaz_R_Eta","ESD_ConvGamma_Pt_Dcaz_R_Eta_Cat");
			fTreeConvGammaPtDcazCat[iCut]->Branch("Pt",&fPtGamma,"fPtGamma/F");
			fTreeConvGammaPtDcazCat[iCut]->Branch("DcaZPhoton",&fDCAzPhoton,"fDCAzPhoton/F");
      //          fTreeConvGammaPtDcazCat[iCut]->Branch("R",&fRConvPhoton,"fRConvPhoton/F");
      //          fTreeConvGammaPtDcazCat[iCut]->Branch("Eta",&fEtaPhoton,"fEtaPhoton/F");
			
			fTreeConvGammaPtDcazCat[iCut]->Branch("cat",&fCharCatPhoton,"fCharCatPhoton/b");
			if(fIsMC){
				fTreeConvGammaPtDcazCat[iCut]->Branch("photonMCInfo",&fCharPhotonMCInfo,"fCharPhotonMCInfo/b");
			}
			fPhotonDCAList[iCut]->Add(fTreeConvGammaPtDcazCat[iCut]);
		}
    
		if (fDoPhotonQA > 0){
			fHistoConvGammaR[iCut] = new TH1F("ESD_ConvGamma_R","ESD_ConvGamma_R",800,0,200);
			fESDList[iCut]->Add(fHistoConvGammaR[iCut]);
			fHistoConvGammaEta[iCut] = new TH1F("ESD_ConvGamma_Eta","ESD_ConvGamma_Eta",2000,-2,2);
			fESDList[iCut]->Add(fHistoConvGammaEta[iCut]);
		}

		fTagOutputList[iCut] = new TList();
		fTagOutputList[iCut]->SetName(Form("%s_%s_%s_%s Tagging Output",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
		fTagOutputList[iCut]->SetOwner(1);
		fCutFolder[iCut]->Add(fTagOutputList[iCut]);
		
		const Int_t nptbins = 200;
		const Double_t ptmin = 0.;
		const Double_t ptmax = 20.;
		
		const Int_t nmbins = 180;
		const Double_t mmin = 0.;
		const Double_t mmax = 0.9;
		
		// photon candidates
		// this is maybe not necessary ...
			
		fHistoConvGammaUntagged[iCut] = new TH1F("ConvGammaUntagged","",nptbins,ptmin,ptmax);
		fHistoConvGammaUntagged[iCut]->SetXTitle("p_{T} (GeV/c)");
		fTagOutputList[iCut]->Add(fHistoConvGammaUntagged[iCut]);
		
		fHistoConvGammaTagged[iCut] = new TH1F("ConvGammaTagged","",nptbins,ptmin,ptmax);
		fHistoConvGammaTagged[iCut]->SetXTitle("p_{T} (GeV/c)");
		fTagOutputList[iCut]->Add(fHistoConvGammaTagged[iCut]);
		
		fHistoConvGammaPi0Tagged[iCut] = new TH1F("ConvGammaPi0Tagged","",nptbins,ptmin,ptmax);
		fHistoConvGammaPi0Tagged[iCut]->SetXTitle("p_{T} (GeV/c)");
		fTagOutputList[iCut]->Add(fHistoConvGammaPi0Tagged[iCut]);
		
		fHistoConvGammaEtaTagged[iCut] = new TH1F("ConvGammaEtaTagged","",nptbins,ptmin,ptmax);
		fHistoConvGammaEtaTagged[iCut]->SetXTitle("p_{T} (GeV/c)");
		fTagOutputList[iCut]->Add(fHistoConvGammaEtaTagged[iCut]);
		
		// pairs
		fHistoPhotonPairAll[iCut] = new TH2F("PhotonPairAll","",nmbins,mmin,mmax,nptbins,ptmin,ptmax);
		fHistoPhotonPairAll[iCut]->SetXTitle("M_{inv} (GeV/cc)");
		fHistoPhotonPairAll[iCut]->SetYTitle("p_{T} (GeV/c)");
		fTagOutputList[iCut]->Add(fHistoPhotonPairAll[iCut]);
		
		fHistoPhotonPairAllGam[iCut] = new TH2F("PhotonPairAllGammaConvPt","",nmbins,mmin,mmax,nptbins,ptmin,ptmax);
		fHistoPhotonPairAllGam[iCut]->SetXTitle("M_{inv} (GeV/cc)");
		fHistoPhotonPairAllGam[iCut]->SetYTitle("#gamma^{conv} p_{T} (GeV/c)");
		fTagOutputList[iCut]->Add(fHistoPhotonPairAllGam[iCut]);
	
		fHistoClusGammaPt[iCut] = new TH1F("ClusGamma_Pt","ClusGamma_Pt",250,0,25);
		fTagOutputList[iCut]->Add(fHistoClusGammaPt[iCut]);

		
		if(fDoMesonAnalysis){
			fHistoMotherInvMassPt[iCut] = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",800,0,0.8,250,0,25);
			fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);
			fHistoMotherBackInvMassPt[iCut] = new TH2F("ESD_Background_InvMass_Pt","ESD_Background_InvMass_Pt",800,0,0.8,250,0,25);
			fESDList[iCut]->Add(fHistoMotherBackInvMassPt[iCut]);
			fHistoMotherInvMassEalpha[iCut] = new TH2F("ESD_Mother_InvMass_vs_E_alpha","ESD_Mother_InvMass_vs_E_alpha",800,0,0.8,250,0,25);
			fESDList[iCut]->Add(fHistoMotherInvMassEalpha[iCut]);
			if (fDoMesonQA == 2){
				fMesonDCAList[iCut] = new TList();
				fMesonDCAList[iCut]->SetName(Form("%s_%s_%s_%s Meson DCA tree",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
				fMesonDCAList[iCut]->SetOwner(kTRUE);
				fCutFolder[iCut]->Add(fMesonDCAList[iCut]);
				
				fTreeMesonsInvMassPtDcazMinDcazMaxFlag[iCut] = new TTree("ESD_Mesons_InvMass_Pt_DcazMin_DcazMax_Flag","ESD_Mesons_InvMass_Pt_DcazMin_DcazMax_Flag");
				fTreeMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("InvMass",&fInvMass,"fInvMass/F");
				fTreeMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("Pt",&fPt,"fPt/F");
				fTreeMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("DcaZMin",&fDCAzGammaMin,"fDCAzGammaMin/F");
				fTreeMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("DcaZMax",&fDCAzGammaMax,"fDCAzGammaMax/F");
				fTreeMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("kind",&fCharFlag,"fCharFlag/b");
				if(fIsMC){
					fTreeMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("mesonMCInfo",&fCharMesonMCInfo,"fCharMesonMCInfo/b");
				}
				fMesonDCAList[iCut]->Add(fTreeMesonsInvMassPtDcazMinDcazMaxFlag[iCut]);
				
			}
			if(fDoMesonQA == 1){
				fHistoMotherInvMassECalib[iCut] = new TH2F("ESD_Mother_InvMass_Pt_Calib","ESD_Mother_InvMass_Pt_Calib",800,0,0.8,250,0,25);
				fESDList[iCut]->Add(fHistoMotherInvMassECalib[iCut]);
				fHistoMotherInvMassECalibalpha[iCut] = new TH2F("ESD_Mother_InvMass_vs_E_Calib_alpha","ESD_Mother_InvMass_vs_E_Calib_alpha",800,0,0.8,250,0,25);
				fESDList[iCut]->Add(fHistoMotherInvMassECalibalpha[iCut]);
			}

			if (fDoMesonQA > 0 ){
				fHistoMotherPi0PtY[iCut] = new TH2F("ESD_MotherPi0_Pt_Y","ESD_MotherPi0_Pt_Y",150,0.03,15.,150,-1.5,1.5);
				SetLogBinningXTH2(fHistoMotherPi0PtY[iCut]);
				fESDList[iCut]->Add(fHistoMotherPi0PtY[iCut]);
				fHistoMotherEtaPtY[iCut] = new TH2F("ESD_MotherEta_Pt_Y","ESD_MotherEta_Pt_Y",150,0.03,15.,150,-1.5,1.5);
				SetLogBinningXTH2(fHistoMotherEtaPtY[iCut]);
				fESDList[iCut]->Add(fHistoMotherEtaPtY[iCut]);
				fHistoMotherPi0PtAlpha[iCut] = new TH2F("ESD_MotherPi0_Pt_Alpha","ESD_MotherPi0_Pt_Alpha",150,0.03,15.,100,0,1);
				SetLogBinningXTH2(fHistoMotherPi0PtAlpha[iCut]);
				fESDList[iCut]->Add(fHistoMotherPi0PtAlpha[iCut]);
				fHistoMotherEtaPtAlpha[iCut] = new TH2F("ESD_MotherEta_Pt_Alpha","ESD_MotherEta_Pt_Alpha",150,0.03,15.,100,0,1);
				SetLogBinningXTH2(fHistoMotherEtaPtAlpha[iCut]);
				fESDList[iCut]->Add(fHistoMotherEtaPtAlpha[iCut]);
				fHistoMotherPi0PtOpenAngle[iCut] = new TH2F("ESD_MotherPi0_Pt_OpenAngle","ESD_MotherPi0_Pt_OpenAngle",150,0.03,15.,200,0,2*TMath::Pi());
				SetLogBinningXTH2(fHistoMotherPi0PtOpenAngle[iCut]);
				fESDList[iCut]->Add(fHistoMotherPi0PtOpenAngle[iCut]);
				fHistoMotherEtaPtOpenAngle[iCut] = new TH2F("ESD_MotherEta_Pt_OpenAngle","ESD_MotherEta_Pt_OpenAngle",150,0.03,15.,200,0,2*TMath::Pi());
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
		fHistoMCAllGammaEMCALAccPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaPi0Pt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaRhoPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaEtaPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaOmegaPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaEtapPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaPhiPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaSigmaPt 			= new TH1F*[fnCuts];
		fHistoMCConvGammaPt 				= new TH1F*[fnCuts];
		fHistoTrueConvGammaPt 				= new TH1F*[fnCuts];
		fHistoTrueConvPi0GammaPt 			= new TH1F*[fnCuts];
    
		fHistoCombinatorialPt 							= new TH2F*[fnCuts];
		fHistoTruePrimaryConvGammaPt 					= new TH1F*[fnCuts];
		fHistoTruePrimaryConvGammaESDPtMCPt 			= new TH2F*[fnCuts];
		fHistoTrueSecondaryConvGammaPt 					= new TH1F*[fnCuts];
		fHistoTrueSecondaryConvGammaFromXFromK0sPt 		= new TH1F*[fnCuts];
		fHistoTrueSecondaryConvGammaFromXFromLambdaPt 	= new TH1F*[fnCuts];
    
		fHistoTrueClusGammaPt 				= new TH1F*[fnCuts];
		fHistoTruePrimaryClusGammaPt 		= new TH1F*[fnCuts];
		fHistoTruePrimaryClusGammaESDPtMCPt = new TH2F*[fnCuts];

		if (fDoPhotonQA > 0){
			fHistoMCConvGammaR 					= new TH1F*[fnCuts];
			fHistoMCConvGammaEta 				= new TH1F*[fnCuts];
			fHistoTrueConvGammaEta 				= new TH1F*[fnCuts];
		}
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
      
			fHistoTrueMotherInvMassPt 					= new TH2F*[fnCuts];
			fHistoTruePrimaryMotherInvMassPt 			= new TH2F*[fnCuts];
			fHistoTruePrimaryMotherW0WeightingInvMassPt = new TH2F*[fnCuts];
			fProfileTruePrimaryMotherWeightsInvMassPt 	= new TProfile2D*[fnCuts];
			fHistoTrueSecondaryMotherInvMassPt 			= new TH2F*[fnCuts];
			fHistoTrueSecondaryMotherFromK0sInvMassPt 	= new TH2F*[fnCuts];
			fHistoTrueSecondaryMotherFromEtaInvMassPt 	= new TH2F*[fnCuts];
			fHistoTrueSecondaryMotherFromLambdaInvMassPt = new TH2F*[fnCuts];
			if (fDoMesonQA > 0){
				fHistoMCPi0PtY 								= new TH2F*[fnCuts];
				fHistoMCEtaPtY 								= new TH2F*[fnCuts];
				fHistoMCK0sPt 								= new TH1F*[fnCuts];
				fHistoMCK0sWOWeightPt 						= new TH1F*[fnCuts];
				fHistoMCK0sPtY	 							= new TH2F*[fnCuts];
				fHistoMCSecPi0PtvsSource 					= new TH2F*[fnCuts];
				fHistoMCSecPi0Source 						= new TH1F*[fnCuts];
				fHistoMCSecEtaPt 							= new TH1F*[fnCuts];
				fHistoMCSecEtaSource 						= new TH1F*[fnCuts];
				fHistoTrueMotherCaloPhotonInvMassPt			= new TH2F*[fnCuts];
				fHistoTrueMotherCaloConvertedPhotonInvMassPt= new TH2F*[fnCuts];
				fHistoTruePi0CaloConvertedPhotonInvMassPt 	= new TH2F*[fnCuts];
				fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt = new TH2F*[fnCuts];
				fHistoTrueEtaCaloConvertedPhotonInvMassPt	= new TH2F*[fnCuts];
				fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt = new TH2F*[fnCuts];
				fHistoTrueMotherCaloElectronInvMassPt		= new TH2F*[fnCuts];
				fHistoTrueMotherCaloMergedClusterInvMassPt 	= new TH2F*[fnCuts];
				fHistoTrueMotherCaloMergedClusterPartConvInvMassPt 	= new TH2F*[fnCuts];
				fHistoTrueMotherCaloEMNonLeadingInvMassPt 	= new TH2F*[fnCuts];
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
		}
    
    
    
		for(Int_t iCut = 0; iCut<fnCuts;iCut++){
			TString cutstringEvent 	= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
			TString cutstringPhoton = ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
			TString cutstringCalo 	= ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
			TString cutstringMeson 	= "NoMesonCut";
			if(fDoMesonAnalysis)cutstringMeson = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
      
			fMCList[iCut] = new TList();
			fMCList[iCut]->SetName(Form("%s_%s_%s_%s MC histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
			fMCList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fMCList[iCut]);
			fHistoMCHeaders[iCut] = new TH1I("MC_Headers","MC_Headers",20,0,20);
			fMCList[iCut]->Add(fHistoMCHeaders[iCut]);
			fHistoMCAllGammaPt[iCut] = new TH1F("MC_AllGamma_Pt","MC_AllGamma_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCAllGammaPt[iCut]);
			fHistoMCAllGammaEMCALAccPt[iCut] = new TH1F("MC_AllGammaEMCALAcc_Pt","MC_AllGamma_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCAllGammaEMCALAccPt[iCut]);
			fHistoMCDecayGammaPi0Pt[iCut] = new TH1F("MC_DecayGammaPi0_Pt","MC_DecayGammaPi0_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCDecayGammaPi0Pt[iCut]);
			fHistoMCDecayGammaRhoPt[iCut] = new TH1F("MC_DecayGammaRho_Pt","MC_DecayGammaRho_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCDecayGammaRhoPt[iCut]);
			fHistoMCDecayGammaEtaPt[iCut] = new TH1F("MC_DecayGammaEta_Pt","MC_DecayGammaEta_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCDecayGammaEtaPt[iCut]);
			fHistoMCDecayGammaOmegaPt[iCut] = new TH1F("MC_DecayGammaOmega_Pt","MC_DecayGammaOmmega_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCDecayGammaOmegaPt[iCut]);
			fHistoMCDecayGammaEtapPt[iCut] = new TH1F("MC_DecayGammaEtap_Pt","MC_DecayGammaEtap_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCDecayGammaEtapPt[iCut]);
			fHistoMCDecayGammaPhiPt[iCut] = new TH1F("MC_DecayGammaPhi_Pt","MC_DecayGammaPhi_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCDecayGammaPhiPt[iCut]);
			fHistoMCDecayGammaSigmaPt[iCut] = new TH1F("MC_DecayGammaSigma_Pt","MC_DecayGammaSigma_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCDecayGammaSigmaPt[iCut]);
			fHistoMCConvGammaPt[iCut] = new TH1F("MC_ConvGamma_Pt","MC_ConvGamma_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCConvGammaPt[iCut]);
      
			if (fDoPhotonQA > 0){
				fHistoMCConvGammaR[iCut] = new TH1F("MC_ConvGamma_R","MC_ConvGamma_R",800,0,200);
				fMCList[iCut]->Add(fHistoMCConvGammaR[iCut]);
				fHistoMCConvGammaEta[iCut] = new TH1F("MC_ConvGamma_Eta","MC_ConvGamma_Eta",2000,-2,2);
				fMCList[iCut]->Add(fHistoMCConvGammaEta[iCut]);
			}
      
			if(fDoMesonAnalysis){
				fHistoMCPi0Pt[iCut] = new TH1F("MC_Pi0_Pt","MC_Pi0_Pt",250,0,25);
				fHistoMCPi0Pt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCPi0Pt[iCut]);
				fHistoMCPi0WOWeightPt[iCut] = new TH1F("MC_Pi0_WOWeights_Pt","MC_Pi0_WOWeights_Pt",250,0,25);
				fHistoMCPi0WOWeightPt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCPi0WOWeightPt[iCut]);
				
				fHistoMCEtaPt[iCut] = new TH1F("MC_Eta_Pt","MC_Eta_Pt",250,0,25);
				fHistoMCEtaPt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCEtaPt[iCut]);
				fHistoMCEtaWOWeightPt[iCut] = new TH1F("MC_Eta_WOWeights_Pt","MC_Eta_WOWeights_Pt",250,0,25);
				fHistoMCEtaWOWeightPt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCEtaWOWeightPt[iCut]);
				fHistoMCPi0InAccPt[iCut] = new TH1F("MC_Pi0InAcc_Pt","MC_Pi0InAcc_Pt",250,0,25);
				fHistoMCPi0InAccPt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCPi0InAccPt[iCut]);
				fHistoMCEtaInAccPt[iCut] = new TH1F("MC_EtaInAcc_Pt","MC_EtaInAcc_Pt",250,0,25);
				fHistoMCEtaInAccPt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCEtaInAccPt[iCut]);
				if (fDoMesonQA > 0){
					fHistoMCPi0PtY[iCut] = new TH2F("MC_Pi0_Pt_Y","MC_Pi0_Pt_Y",150,0.03,15.,150,-1.5,1.5);
					fHistoMCPi0PtY[iCut]->Sumw2();
					SetLogBinningXTH2(fHistoMCPi0PtY[iCut]);
					fMCList[iCut]->Add(fHistoMCPi0PtY[iCut]);
					fHistoMCEtaPtY[iCut] = new TH2F("MC_Eta_Pt_Y","MC_Eta_Pt_Y",150,0.03,15.,150,-1.5,1.5);
					fHistoMCEtaPtY[iCut]->Sumw2();
					SetLogBinningXTH2(fHistoMCEtaPtY[iCut]);
					fMCList[iCut]->Add(fHistoMCEtaPtY[iCut]);
					fHistoMCK0sPt[iCut] = new TH1F("MC_K0s_Pt","MC_K0s_Pt",150,0,15);
					fHistoMCK0sPt[iCut]->Sumw2();
					fMCList[iCut]->Add(fHistoMCK0sPt[iCut]);
					fHistoMCK0sWOWeightPt[iCut] = new TH1F("MC_K0s_WOWeights_Pt","MC_K0s_WOWeights_Pt",150,0,15);
					fHistoMCK0sWOWeightPt[iCut]->Sumw2();
					fMCList[iCut]->Add(fHistoMCK0sWOWeightPt[iCut]);
					fHistoMCK0sPtY[iCut] = new TH2F("MC_K0s_Pt_Y","MC_K0s_Pt_Y",150,0.03,15.,150,-1.5,1.5);
					fHistoMCK0sPtY[iCut]->Sumw2();
					SetLogBinningXTH2(fHistoMCK0sPtY[iCut]);
					fMCList[iCut]->Add(fHistoMCK0sPtY[iCut]);
					
					fHistoMCSecPi0Source[iCut] = new TH1F("MC_SecPi0_Source","MC_SecPi0_Source",5000,0.,5000);
					fMCList[iCut]->Add(fHistoMCSecPi0Source[iCut]);
					fHistoMCSecEtaSource[iCut] = new TH1F("MC_SecEta_Source","MC_SecEta_Source",5000,0,5000);
					fMCList[iCut]->Add(fHistoMCSecEtaSource[iCut]);
					fHistoMCSecPi0PtvsSource[iCut] = new TH2F("MC_SecPi0_Pt_Source","MC_SecPi0_Pt_Source",250,0.0,25.,16,-0.5,15.5);
					fHistoMCSecPi0PtvsSource[iCut]->Sumw2();
					fMCList[iCut]->Add(fHistoMCSecPi0PtvsSource[iCut]);
					fHistoMCSecEtaPt[iCut] = new TH1F("MC_SecEta_Pt","MC_SecEta_Pt",250,0,25);
					fHistoMCSecEtaPt[iCut]->Sumw2();
					fMCList[iCut]->Add(fHistoMCSecEtaPt[iCut]);
				}
        
			}
			fTrueList[iCut] = new TList();
			fTrueList[iCut]->SetName(Form("%s_%s_%s_%s True histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
			fTrueList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fTrueList[iCut]);
      
			fHistoTrueConvGammaPt[iCut] = new TH1F("ESD_TrueConvGamma_Pt","ESD_TrueConvGamma_Pt",250,0,25);
			fTrueList[iCut]->Add(fHistoTrueConvGammaPt[iCut]);
      
			fHistoTrueConvPi0GammaPt[iCut] = new TH1F("ESD_TrueConvGamma_Pt","ESD_TrueConvGamma_Pt",250,0,25);
			fTrueList[iCut]->Add(fHistoTrueConvPi0GammaPt[iCut]);
      
			fHistoCombinatorialPt[iCut] = new TH2F("ESD_TrueCombinatorial_Pt","ESD_TrueCombinatorial_Pt",250,0,25,16,-0.5,15.5);
			fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 1,"Elec+Elec");
			fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 2,"Elec+Pion");
			fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 3,"Elec+Kaon");
			fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 4,"Elec+Proton");
			fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 5,"Elec+Muon");
			fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 6,"Pion+Pion");
			fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 7,"Pion+Kaon");
			fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 8,"Pion+Proton");
			fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 9,"Pion+Muon");
			fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(10,"Kaon+Kaon");
			fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(11,"Kaon+Proton");
			fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(12,"Kaon+Muon");
			fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(13,"Proton+Proton");
			fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(14,"Proton+Muon");
			fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(15,"Muon+Muon");
			fHistoCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(16,"Rest");
			fTrueList[iCut]->Add(fHistoCombinatorialPt[iCut]);
			fHistoTruePrimaryConvGammaPt[iCut] = new TH1F("ESD_TruePrimaryConvGamma_Pt","ESD_TruePrimaryConvGamma_Pt",250,0,25);
			fTrueList[iCut]->Add(fHistoTruePrimaryConvGammaPt[iCut]);
			fHistoTrueSecondaryConvGammaPt[iCut] = new TH1F("ESD_TrueSecondaryConvGamma_Pt","ESD_TrueSecondaryConvGamma_Pt",250,0,25);
			fTrueList[iCut]->Add(fHistoTrueSecondaryConvGammaPt[iCut]);
      
			fHistoTrueSecondaryConvGammaFromXFromK0sPt[iCut] = new TH1F("ESD_TrueSecondaryConvGammaFromXFromK0s_Pt", "ESD_TrueSecondaryConvGammaFromXFromK0s_Pt",250,0,25);
			fTrueList[iCut]->Add(fHistoTrueSecondaryConvGammaFromXFromK0sPt[iCut]);
			fHistoTrueSecondaryConvGammaFromXFromLambdaPt[iCut] = new TH1F("ESD_TrueSecondaryConvGammaFromXFromLambda_Pt", "ESD_TrueSecondaryConvGammaFromXFromLambda_Pt",250,0,25);
			fTrueList[iCut]->Add(fHistoTrueSecondaryConvGammaFromXFromLambdaPt[iCut]);
      			      
			fHistoTruePrimaryConvGammaESDPtMCPt[iCut] = new TH2F("ESD_TruePrimaryConvGammaESD_PtMCPt", "ESD_TruePrimaryConvGammaESD_PtMCPt",250,0,25,250,0,25);
			fTrueList[iCut]->Add(fHistoTruePrimaryConvGammaESDPtMCPt[iCut]);
		
			fHistoTrueClusGammaPt[iCut] = new TH1F("TrueClusGamma_Pt","ESD_TrueClusGamma_Pt",250,0,25);
			fTagOutputList[iCut]->Add(fHistoTrueClusGammaPt[iCut]);
			fHistoTruePrimaryClusGammaPt[iCut] = new TH1F("TruePrimaryClusGamma_Pt","ESD_TruePrimaryClusGamma_Pt",250,0,25);
			fTagOutputList[iCut]->Add(fHistoTruePrimaryClusGammaPt[iCut]);
			fHistoTruePrimaryClusGammaESDPtMCPt[iCut] = new TH2F("TruePrimaryClusGamma_Pt_MCPt","ESD_TruePrimaryClusGamma_MCPt",250,0,25,250,0,25);
			fTagOutputList[iCut]->Add(fHistoTruePrimaryClusGammaESDPtMCPt[iCut]);

			if (fDoPhotonQA > 0){
				fHistoTrueConvGammaEta[iCut] = new TH1F("ESD_TrueConvGamma_Eta","ESD_TrueConvGamma_Eta",2000,-2,2);
				fTrueList[iCut]->Add(fHistoTrueConvGammaEta[iCut]);		
			}	
			if (fDoClusterQA > 0){	
				fHistoTrueClusUnConvGammaPt[iCut] = new TH1F("TrueClusUnConvGamma_Pt","TrueClusUnConvGamma_Pt",250,0,25);
				fTagOutputList[iCut]->Add(fHistoTrueClusUnConvGammaPt[iCut]);
				fHistoTrueClusUnConvGammaMCPt[iCut] = new TH1F("TrueClusUnConvGamma_MCPt","TrueClusUnConvGamma_MCPt",250,0,25);
				fTagOutputList[iCut]->Add(fHistoTrueClusUnConvGammaMCPt[iCut]);
				fHistoTrueClusElectronPt[iCut] = new TH1F("TrueClusElectron_Pt","TrueElectronGamma_Pt",250,0,25);
				fTagOutputList[iCut]->Add(fHistoTrueClusElectronPt[iCut]);
				fHistoTrueClusConvGammaPt[iCut] = new TH1F("TrueClusConvGamma_Pt","TrueClusConvGamma_Pt",250,0,25);
				fTagOutputList[iCut]->Add(fHistoTrueClusConvGammaPt[iCut]);
				fHistoTrueClusConvGammaMCPt[iCut] = new TH1F("TrueClusConvGamma_MCPt","TrueClusConvGamma_MCPt",250,0,25);
				fTagOutputList[iCut]->Add(fHistoTrueClusConvGammaMCPt[iCut]);
				fHistoTrueClusConvGammaFullyPt[iCut] = new TH1F("TrueClusConvGammaFullyContained_Pt","TrueClusConvGammaFullyContained_Pt",250,0,25);
				fTagOutputList[iCut]->Add(fHistoTrueClusConvGammaFullyPt[iCut]);
				fHistoTrueClusMergedGammaPt[iCut] = new TH1F("TrueClusMergedGamma_Pt","TrueClusMergedGamma_Pt",250,0,25);
				fTagOutputList[iCut]->Add(fHistoTrueClusMergedGammaPt[iCut]);
				fHistoTrueClusMergedPartConvGammaPt[iCut] = new TH1F("TrueClusMergedPartConvGamma_Pt","TrueClusMergedPartConvGamma_Pt",250,0,25);
				fTagOutputList[iCut]->Add(fHistoTrueClusMergedPartConvGammaPt[iCut]);
				fHistoTrueClusDalitzPt[iCut] = new TH1F("TrueClusDalitz_Pt","TrueClusDalitz_Pt",250,0,25);
				fTagOutputList[iCut]->Add(fHistoTrueClusDalitzPt[iCut]);
				fHistoTrueClusDalitzMergedPt[iCut] = new TH1F("TrueClusDalitzMerged_Pt","TrueClusDalitzMerged_Pt",250,0,25);
				fTagOutputList[iCut]->Add(fHistoTrueClusDalitzMergedPt[iCut]);
				fHistoTrueClusPhotonFromElecMotherPt[iCut] = new TH1F("TrueClusPhotonFromElecMother_Pt","TrueClusPhotonFromElecMother_Pt",250,0,25);
				fTagOutputList[iCut]->Add(fHistoTrueClusPhotonFromElecMotherPt[iCut]);
				fHistoTrueClusShowerPt[iCut] = new TH1F("TrueClusShower_Pt","TrueClusShower_Pt",250,0,25);
				fTagOutputList[iCut]->Add(fHistoTrueClusShowerPt[iCut]);
				fHistoTrueClusSubLeadingPt[iCut] = new TH1F("TrueClusSubleading_Pt","TrueClusSubleading_Pt",250,0,25);
				fTagOutputList[iCut]->Add(fHistoTrueClusSubLeadingPt[iCut]);
				fHistoTrueClusNParticles[iCut] = new TH1I("TrueClusNParticles","TrueClusNParticles",20,0,20);
				fTagOutputList[iCut]->Add(fHistoTrueClusNParticles[iCut]);
				fHistoTrueClusEMNonLeadingPt[iCut] = new TH1F("TrueClusEMNonLeading_Pt","TrueClusEMNonLeading_Pt",250,0,25);
				fTagOutputList[iCut]->Add(fHistoTrueClusEMNonLeadingPt[iCut]);
				fHistoTrueNLabelsInClus[iCut] = new TH1F("TrueNLabelsInClus","TrueNLabelsInClus",100,-0.5,99.5);
				fTagOutputList[iCut]->Add(fHistoTrueNLabelsInClus[iCut]);	
			}	

			if(fDoMesonAnalysis){
				fHistoTrueMotherInvMassPt[iCut] = new TH2F("ESD_TrueMother_InvMass_Pt","ESD_TrueMother_InvMass_Pt",800,0,0.8,250,0,25);
				fTrueList[iCut]->Add(fHistoTrueMotherInvMassPt[iCut]);
				fHistoTruePrimaryMotherInvMassPt[iCut] = new TH2F("ESD_TruePrimaryMother_InvMass_Pt", "ESD_TruePrimaryMother_InvMass_Pt", 800,0,0.8,250,0,25);
				fHistoTruePrimaryMotherInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTruePrimaryMotherInvMassPt[iCut]);
				fHistoTruePrimaryMotherW0WeightingInvMassPt[iCut] = new TH2F("ESD_TruePrimaryMotherW0Weights_InvMass_Pt", "ESD_TruePrimaryMotherW0Weights_InvMass_Pt", 800,0,0.8,250,0,25);
				fHistoTruePrimaryMotherW0WeightingInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTruePrimaryMotherW0WeightingInvMassPt[iCut]);
				fProfileTruePrimaryMotherWeightsInvMassPt[iCut] = new TProfile2D("ESD_TruePrimaryMotherWeights_InvMass_Pt", "ESD_TruePrimaryMotherWeights_InvMass_Pt", 800,0,0.8,250,0,25);
				fProfileTruePrimaryMotherWeightsInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fProfileTruePrimaryMotherWeightsInvMassPt[iCut]);
				fHistoTrueSecondaryMotherInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryMother_InvMass_Pt", "ESD_TrueSecondaryMother_InvMass_Pt", 800,0,0.8,250,0,25);
				fHistoTrueSecondaryMotherInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTrueSecondaryMotherInvMassPt[iCut]);
				fHistoTrueSecondaryMotherFromK0sInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryMotherFromK0s_InvMass_Pt","ESD_TrueSecondaryMotherFromK0s_InvMass_Pt",800,0,0.8,250,0,25);
				fHistoTrueSecondaryMotherFromK0sInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTrueSecondaryMotherFromK0sInvMassPt[iCut]);
				fHistoTrueSecondaryMotherFromEtaInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryMotherFromEta_InvMass_Pt","ESD_TrueSecondaryMotherFromEta_InvMass_Pt",800,0,0.8,250,0,25);
				fTrueList[iCut]->Add(fHistoTrueSecondaryMotherFromEtaInvMassPt[iCut]);
				fHistoTrueSecondaryMotherFromLambdaInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryMotherFromLambda_InvMass_Pt","ESD_TrueSecondaryMotherFromLambda_InvMass_Pt",800,0,0.8,250,0,25);
				fTrueList[iCut]->Add(fHistoTrueSecondaryMotherFromLambdaInvMassPt[iCut]);
				if (fDoMesonQA > 0){
					fHistoTrueMotherCaloPhotonInvMassPt[iCut] = new TH2F("ESD_TrueMotherCaloPhoton_InvMass_Pt","ESD_TrueMotherCaloPhoton_InvMass_Pt",800,0,0.8,250,0,25);
					fTrueList[iCut]->Add(fHistoTrueMotherCaloPhotonInvMassPt[iCut]);
					fHistoTrueMotherCaloConvertedPhotonInvMassPt[iCut] = new TH2F("ESD_TrueMotherCaloConvertedPhoton_InvMass_Pt","ESD_TrueMotherCaloConvertedPhoton_InvMass_Pt",800,0,0.8,250,0,25);
					fTrueList[iCut]->Add(fHistoTrueMotherCaloConvertedPhotonInvMassPt[iCut]);

					fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut] = new TH2F("ESD_TruePi0CaloConvertedPhoton_InvMass_Pt","ESD_TruePi0CaloConvertedPhoton_InvMass_Pt",800,0,0.8,250,0,25);
					fTrueList[iCut]->Add(fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut]);
					fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt[iCut] = new TH2F("ESD_TruePi0CaloConvertedPhotonMatched_InvMass_Pt","ESD_TruePi0CaloConvertedPhotonMatched_InvMass_Pt",800,0,0.8,250,0,25);
					fTrueList[iCut]->Add(fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt[iCut]);

					fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut] = new TH2F("ESD_TrueEtaCaloConvertedPhoton_InvMass_Pt","ESD_TrueEtaCaloConvertedPhoton_InvMass_Pt",800,0,0.8,250,0,25);
					fTrueList[iCut]->Add(fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut]);
					fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt[iCut] = new TH2F("ESD_TrueEtaCaloConvertedPhotonMatched_InvMass_Pt","ESD_TrueEtaCaloConvertedPhotonMatched_InvMass_Pt",800,0,0.8,250,0,25);
					fTrueList[iCut]->Add(fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt[iCut]);

					fHistoTrueMotherCaloElectronInvMassPt[iCut] = new TH2F("ESD_TrueMotherCaloElectron_InvMass_Pt","ESD_TrueMotherCaloElectron_InvMass_Pt",800,0,0.8,250,0,25);
					fTrueList[iCut]->Add(fHistoTrueMotherCaloElectronInvMassPt[iCut]);
					fHistoTrueMotherCaloMergedClusterInvMassPt[iCut] = new TH2F("ESD_TrueMotherCaloMergedCluster_InvMass_Pt","ESD_TrueMotherCaloMergedCluster_InvMass_Pt",800,0,0.8,250,0,25);
					fTrueList[iCut]->Add(fHistoTrueMotherCaloMergedClusterInvMassPt[iCut]);
					fHistoTrueMotherCaloEMNonLeadingInvMassPt[iCut] = new TH2F("ESD_TrueMotherCaloEMNonLeading_InvMass_Pt","ESD_TrueMotherCaloEMNonLeading_InvMass_Pt",800,0,0.8,250,0,25);
					fTrueList[iCut]->Add(fHistoTrueMotherCaloEMNonLeadingInvMassPt[iCut]);
					fHistoTrueMotherCaloMergedClusterPartConvInvMassPt[iCut] = new TH2F("ESD_TrueMotherCaloMergedClusterPartConv_InvMass_Pt","ESD_TrueMotherCaloMergedClusterPartConv_InvMass_Pt",800,0,0.8,250,0,25);
					fTrueList[iCut]->Add(fHistoTrueMotherCaloMergedClusterPartConvInvMassPt[iCut]);
					
					fHistoTruePrimaryPi0MCPtResolPt[iCut] = new TH2F("ESD_TruePrimaryPi0_MCPt_ResolPt","ESD_TruePrimaryPi0_ResolPt_MCPt",500,0.03,25,1000,-1.,1.);
					fHistoTruePrimaryPi0MCPtResolPt[iCut]->Sumw2();
					SetLogBinningXTH2(fHistoTruePrimaryPi0MCPtResolPt[iCut]);
					fTrueList[iCut]->Add(fHistoTruePrimaryPi0MCPtResolPt[iCut]);
					fHistoTruePrimaryEtaMCPtResolPt[iCut]  = new TH2F("ESD_TruePrimaryEta_MCPt_ResolPt","ESD_TruePrimaryEta_ResolPt_MCPt",500,0.03,25,1000,-1.,1.);
					fHistoTruePrimaryEtaMCPtResolPt[iCut]->Sumw2();
					SetLogBinningXTH2(fHistoTruePrimaryEtaMCPtResolPt[iCut]);
					fTrueList[iCut]->Add(fHistoTruePrimaryEtaMCPtResolPt[iCut]);
					fHistoTrueBckGGInvMassPt[iCut] = new TH2F("ESD_TrueBckGG_InvMass_Pt","ESD_TrueBckGG_InvMass_Pt",800,0,0.8,250,0,25);
					fTrueList[iCut]->Add(fHistoTrueBckGGInvMassPt[iCut]);
					fHistoTrueBckContInvMassPt[iCut] = new TH2F("ESD_TrueBckCont_InvMass_Pt","ESD_TrueBckCont_InvMass_Pt",800,0,0.8,250,0,25);
					fTrueList[iCut]->Add(fHistoTrueBckContInvMassPt[iCut]);
					fHistoTrueK0sWithPi0DaughterMCPt[iCut] = new TH1F("ESD_TrueK0sWithPi0Daughter_MCPt","ESD_TrueK0sWithPi0Daughter_MCPt",250,0,25);
					fTrueList[iCut]->Add(fHistoTrueK0sWithPi0DaughterMCPt[iCut]);
					fHistoTrueEtaWithPi0DaughterMCPt[iCut] = new TH1F("ESD_TrueEtaWithPi0Daughter_MCPt","ESD_TrueEtaWithPi0Daughter_MCPt",250,0,25);
					fTrueList[iCut]->Add(fHistoTrueEtaWithPi0DaughterMCPt[iCut]);
					fHistoTrueLambdaWithPi0DaughterMCPt[iCut] = new TH1F("ESD_TrueLambdaWithPi0Daughter_MCPt","ESD_TrueLambdaWithPi0Daughter_MCPt",250,0,25);
					fTrueList[iCut]->Add(fHistoTrueLambdaWithPi0DaughterMCPt[iCut]);
          
					fHistoTruePi0PtY[iCut] = new TH2F("ESD_TruePi0_Pt_Y","ESD_TruePi0_Pt_Y",150,0.03,15.,150,-1.5,1.5);
					SetLogBinningXTH2(fHistoTruePi0PtY[iCut]);
					fTrueList[iCut]->Add(fHistoTruePi0PtY[iCut]);
					fHistoTrueEtaPtY[iCut] = new TH2F("ESD_TrueEta_Pt_Y","ESD_TrueEta_Pt_Y",150,0.03,15.,150,-1.5,1.5);
					SetLogBinningXTH2(fHistoTrueEtaPtY[iCut]);
					fTrueList[iCut]->Add(fHistoTrueEtaPtY[iCut]);
					fHistoTruePi0PtAlpha[iCut] = new TH2F("ESD_TruePi0_Pt_Alpha","ESD_TruePi0_Pt_Alpha",150,0.03,15.,100,0,1);
					SetLogBinningXTH2(fHistoTruePi0PtAlpha[iCut]);
					fTrueList[iCut]->Add(fHistoTruePi0PtAlpha[iCut]);
					fHistoTrueEtaPtAlpha[iCut] = new TH2F("ESD_TrueEta_Pt_Alpha","ESD_TrueEta_Pt_Alpha",150,0.03,15.,100,0,1);
					SetLogBinningXTH2(fHistoTrueEtaPtAlpha[iCut]);
					fTrueList[iCut]->Add(fHistoTrueEtaPtAlpha[iCut]);
					
					fHistoTruePi0PtOpenAngle[iCut] = new TH2F("ESD_TruePi0_Pt_OpenAngle","ESD_TruePi0_Pt_OpenAngle",150,0.03,15.,200,0,2*TMath::Pi());
					SetLogBinningXTH2(fHistoTruePi0PtOpenAngle[iCut]);
					fTrueList[iCut]->Add(fHistoTruePi0PtOpenAngle[iCut]);
					fHistoTrueEtaPtOpenAngle[iCut] = new TH2F("ESD_TrueEta_Pt_OpenAngle","ESD_TrueEta_Pt_OpenAngle",150,0.03,15.,200,0,2*TMath::Pi());
					SetLogBinningXTH2(fHistoTrueEtaPtOpenAngle[iCut]);
					fTrueList[iCut]->Add(fHistoTrueEtaPtOpenAngle[iCut]);
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
		if(!((AliConversionPhotonCuts*)fCutArray->At(iCut))) continue;
		if(((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutHistograms()){
			fCutFolder[iCut]->Add(((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutHistograms());
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
Bool_t AliAnalysisTaskGammaConvCalo::Notify()
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
void AliAnalysisTaskGammaConvCalo::UserExec(Option_t *)
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
	
	fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut

	// ------------------- BeginEvent ----------------------------
	
	AliEventplane *EventPlane = fInputEvent->GetEventplane();
	if(fIsHeavyIon ==1)fEventPlaneAngle = EventPlane->GetEventplane("V0",fInputEvent,2);
	else fEventPlaneAngle=0.0;
	
	if(fIsMC && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
		RelabelAODPhotonCandidates(kTRUE);    // In case of AODMC relabeling MC
		fV0Reader->RelabelAODs(kTRUE);
	}
	

	for(Int_t iCut = 0; iCut<fnCuts; iCut++){
		
		fiCut = iCut;
		Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon);
		
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
		ProcessPhotonCandidates(); 				// Process this cuts gammas

		fHistoNGammaCandidates[iCut]->Fill(fGammaCandidates->GetEntries());
		fHistoNGoodESDTracksVsNGammaCanditates[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),fGammaCandidates->GetEntries());
		if(fDoMesonAnalysis){ // Meson Analysis
		if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseMCPSmearing() && fIsMC){
			fUnsmearedPx = new Double_t[fGammaCandidates->GetEntries()]; // Store unsmeared Momenta
			fUnsmearedPy = new Double_t[fGammaCandidates->GetEntries()];
			fUnsmearedPz = new Double_t[fGammaCandidates->GetEntries()];
			fUnsmearedE =  new Double_t[fGammaCandidates->GetEntries()];
			
			for(Int_t gamma=0;gamma<fGammaCandidates->GetEntries();gamma++){ // Smear the AODPhotons in MC
			fUnsmearedPx[gamma] = ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->Px();
			fUnsmearedPy[gamma] = ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->Py();
			fUnsmearedPz[gamma] = ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->Pz();
			fUnsmearedE[gamma] =  ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->E();
			((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->SmearParticle(dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(gamma)));
			}
		}

		PhotonTagging(); // tag PCM photons with calorimeter

		CalculatePi0Candidates(); // Combine Gammas from conversion and from calo

		if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
			if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
			CalculateBackground(); // Combinatorial Background
			UpdateEventByEventData(); // Store Event for mixed Events
			}
			else{
			CalculateBackgroundRP(); // Combinatorial Background
			fBGHandlerRP[iCut]->AddEvent(fGammaCandidates,fInputEvent); // Store Event for mixed Events
			fBGClusHandlerRP[iCut]->AddEvent(fClusterCandidates,fInputEvent); // Store Event for mixed Events
			}
		}
		if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseMCPSmearing() && fIsMC){
			for(Int_t gamma=0;gamma<fGammaCandidates->GetEntries();gamma++){ // Smear the AODPhotons in MC
			((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetPx(fUnsmearedPx[gamma]); // Reset Unsmeared Momenta
			((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetPy(fUnsmearedPy[gamma]);
			((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetPz(fUnsmearedPz[gamma]);
			((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetE(fUnsmearedE[gamma]);
			}
			delete[] fUnsmearedPx; fUnsmearedPx = 0x0;
			delete[] fUnsmearedPy; fUnsmearedPy = 0x0;
			delete[] fUnsmearedPz; fUnsmearedPz = 0x0;
			delete[] fUnsmearedE;  fUnsmearedE  = 0x0;
		}
		}

		fGammaCandidates->Clear(); // delete this cuts good gammas
		fClusterCandidates->Clear(); // delete cluster candidates
	}
	
	if(fIsMC && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
		RelabelAODPhotonCandidates(kFALSE); // Back to ESDMC Label
		fV0Reader->RelabelAODs(kFALSE);
	}
	
	PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessClusters()
{
	
	Int_t nclus = 0;
	nclus = fInputEvent->GetNumberOfCaloClusters();
	
// 	cout << nclus << endl;
	
	if(nclus == 0)	return;
	
	// vertex
	Double_t vertex[3] = {0};
	InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
	
	// Loop over EMCal clusters
	for(Int_t i = 0; i < nclus; i++){
		
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
		
		fHistoClusGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
		fClusterCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas
		
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
void AliAnalysisTaskGammaConvCalo::ProcessTrueClusterCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
		
	TParticle *Photon = NULL;
	if (!TruePhotonCandidate->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set task will abort");
	fHistoTrueNLabelsInClus[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels());
	
	if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0)Photon = fMCStack->Particle(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
		else return;
		
	if(Photon == NULL){
	//    cout << "no photon" << endl;
		return;
	}

	TruePhotonCandidate->SetCaloPhotonMCFlags(fMCStack);
	
	// True Photon
	if(fIsFromMBHeader){
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
		}
	}

	if(Photon->GetMother(0) <= fMCStack->GetNprimary()){
		// Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
		if(fIsFromMBHeader){
			fHistoTruePrimaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			fHistoTruePrimaryClusGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt()); // Allways Filled
		}
	}	
	return;
}


//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessTrueClusterCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate)
{
	AliAODMCParticle *Photon = NULL;
	TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
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
	TruePhotonCandidate->SetCaloPhotonMCFlagsAOD(fInputEvent);
	
	// True Photon
	if(fIsFromMBHeader){
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
		}
	}

	// True Photon
	if(fIsFromMBHeader){
		fHistoTrueClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
// 		if (fDoPhotonQA > 0) fHistoTrueConvGammaEta[fiCut]->Fill(TruePhotonCandidate->Eta());
	}

	if(Photon->IsPrimary()){
		// Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
		if(fIsFromMBHeader){
			fHistoTruePrimaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			fHistoTruePrimaryClusGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt()); // Allways Filled
		}
	}	
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessPhotonCandidates()
{
	Int_t nV0 = 0;
	TList *GammaCandidatesStepOne = new TList();
	TList *GammaCandidatesStepTwo = new TList();
	// Loop over Photon Candidates allocated by ReaderV1
	for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
		AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
		if(!PhotonCandidate) continue;
		fIsFromMBHeader = kTRUE;
		if(fIsMC && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
			Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack, fInputEvent);
			if(isPosFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
			Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack, fInputEvent);
			if(isNegFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
			if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
		}
		
		if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
		if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(PhotonCandidate->GetPhotonPhi(),fEventPlaneAngle)) continue;
		if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
		!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
			fGammaCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas
			
			if(fIsFromMBHeader){
				fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
				if (fDoPhotonQA > 0){
				fHistoConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
				fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
				}
			}
			if(fIsMC){
				if(fInputEvent->IsA()==AliESDEvent::Class())
				ProcessTruePhotonCandidates(PhotonCandidate);
				if(fInputEvent->IsA()==AliAODEvent::Class())
				ProcessTruePhotonCandidatesAOD(PhotonCandidate);
			}
			if (fIsFromMBHeader && fDoPhotonQA == 2){
				if (fIsHeavyIon == 1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 12.){
					fPtGamma = PhotonCandidate->Pt();
					fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
					fRConvPhoton = PhotonCandidate->GetConversionRadius();
					fEtaPhoton = PhotonCandidate->GetPhotonEta();
					fCharCatPhoton = PhotonCandidate->GetPhotonQuality();
					fTreeConvGammaPtDcazCat[fiCut]->Fill();
				} else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 16.){
					fPtGamma = PhotonCandidate->Pt();
					fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
					fRConvPhoton = PhotonCandidate->GetConversionRadius();
					fEtaPhoton = PhotonCandidate->GetPhotonEta();
					fCharCatPhoton = PhotonCandidate->GetPhotonQuality();
					fTreeConvGammaPtDcazCat[fiCut]->Fill();
				}
			}
		} else if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
			((AliConversionPhotonCuts*)fCutArray->At(fiCut))->FillElectonLabelArray(PhotonCandidate,nV0);
			nV0++;
			GammaCandidatesStepOne->Add(PhotonCandidate);
		} else if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
				((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
			GammaCandidatesStepTwo->Add(PhotonCandidate);
		}
	}
	if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){
		for(Int_t i = 0;i<GammaCandidatesStepOne->GetEntries();i++){
			AliAODConversionPhoton *PhotonCandidate= (AliAODConversionPhoton*) GammaCandidatesStepOne->At(i);
			if(!PhotonCandidate) continue;
			fIsFromMBHeader = kTRUE;
			if(fMCStack && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
				Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack, fInputEvent);
				Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack, fInputEvent);
				if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
			}
			if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GammaCandidatesStepOne->GetEntries())) continue;
			if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
				fGammaCandidates->Add(PhotonCandidate);
				if(fIsFromMBHeader){
					fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
					if (fDoPhotonQA > 0){
						fHistoConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
						fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
					}
				}
			}
			if(fIsMC){
				if(fInputEvent->IsA()==AliESDEvent::Class())
					ProcessTruePhotonCandidates(PhotonCandidate);
				if(fInputEvent->IsA()==AliAODEvent::Class())
					ProcessTruePhotonCandidatesAOD(PhotonCandidate);
			} else GammaCandidatesStepTwo->Add(PhotonCandidate); // Close v0s cut enabled -> add to list two
			
			if (fIsFromMBHeader && fDoPhotonQA == 2){
				if (fIsHeavyIon ==1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 12.){
					fPtGamma = PhotonCandidate->Pt();
					fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
					fRConvPhoton = PhotonCandidate->GetConversionRadius();
					fEtaPhoton = PhotonCandidate->GetPhotonEta();
					fCharCatPhoton = PhotonCandidate->GetPhotonQuality();
					fTreeConvGammaPtDcazCat[fiCut]->Fill();
				} else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 16.){
					fPtGamma = PhotonCandidate->Pt();
					fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
					fRConvPhoton = PhotonCandidate->GetConversionRadius();
					fEtaPhoton = PhotonCandidate->GetPhotonEta();
					fCharCatPhoton = PhotonCandidate->GetPhotonQuality();
					fTreeConvGammaPtDcazCat[fiCut]->Fill();
				}
			}
		}
	}
	if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
		for(Int_t i = 0;i<GammaCandidatesStepTwo->GetEntries();i++){
			AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GammaCandidatesStepTwo->At(i);
			if(!PhotonCandidate) continue;
			fIsFromMBHeader = kTRUE;
			if(fMCStack && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
				Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack, fInputEvent);
				Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack, fInputEvent);
				if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
			}
			if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GammaCandidatesStepTwo,i)) continue;
			fGammaCandidates->Add(PhotonCandidate); // Add gamma to current cut TList
			if(fIsFromMBHeader){
				fHistoConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
				if (fDoPhotonQA > 0){
					fHistoConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
					fHistoConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
				}
			}
			if(fIsMC){
				if(fInputEvent->IsA()==AliESDEvent::Class())
					ProcessTruePhotonCandidates(PhotonCandidate);
				if(fInputEvent->IsA()==AliAODEvent::Class())
					ProcessTruePhotonCandidatesAOD(PhotonCandidate);
			}
			if (fIsFromMBHeader){
				if (fIsHeavyIon == 1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 12.){
					fPtGamma = PhotonCandidate->Pt();
					fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
					fRConvPhoton = PhotonCandidate->GetConversionRadius();
					fEtaPhoton = PhotonCandidate->GetPhotonEta();
					fCharCatPhoton = PhotonCandidate->GetPhotonQuality();
					fTreeConvGammaPtDcazCat[fiCut]->Fill();
				} else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 16.){
					fPtGamma = PhotonCandidate->Pt();
					fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
					fRConvPhoton = PhotonCandidate->GetConversionRadius();
					fEtaPhoton = PhotonCandidate->GetPhotonEta();
					fCharCatPhoton = PhotonCandidate->GetPhotonQuality();
					fTreeConvGammaPtDcazCat[fiCut]->Fill();
				}
			}
		}
	}
	
	delete GammaCandidatesStepOne;
	GammaCandidatesStepOne = 0x0;
	delete GammaCandidatesStepTwo;
	GammaCandidatesStepTwo = 0x0;
  
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessTruePhotonCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate)
{
	
	Double_t magField = fInputEvent->GetMagneticField();
	if( magField  < 0.0 ){
		magField =  1.0;
	}
	else {
		magField =  -1.0;
	}
	
	TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	AliAODMCParticle *posDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelPositive());
	AliAODMCParticle *negDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelNegative());
	fCharPhotonMCInfo = 0;
	
	if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
	Int_t pdgCode[2] = {abs(posDaughter->GetPdgCode()),abs(negDaughter->GetPdgCode())};
	
	if(posDaughter->GetMother() != negDaughter->GetMother()){
		FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode);
		fCharPhotonMCInfo = 1;
		return;
	}
	else if(posDaughter->GetMother() == -1){
		FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode);
		fCharPhotonMCInfo = 1;
		return;
	}
	
	if(pdgCode[0]!=11 || pdgCode[1]!=11){
		fCharPhotonMCInfo = 1;
		return; //One Particle is not a electron
	}
	
	if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()){
		fCharPhotonMCInfo = 1;
		return; // Same Charge
	}
	
	AliAODMCParticle *Photon = (AliAODMCParticle*) AODMCTrackArray->At(posDaughter->GetMother());	
	if(Photon->GetPdgCode() != 22){
		fCharPhotonMCInfo = 1;
		return; // Mother is no Photon
	}
	
	if(((posDaughter->GetMCProcessCode())) != 5 || ((negDaughter->GetMCProcessCode())) != 5){
		fCharPhotonMCInfo = 1;
		return;// check if the daughters come from a conversion
	}
	// STILL A BUG IN ALIROOT >>8 HAS TPO BE REMOVED AFTER FIX
	
	
	
	// True Photon
	if(fIsFromMBHeader){
		fHistoTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
		if (fDoPhotonQA > 0) fHistoTrueConvGammaEta[fiCut]->Fill(TruePhotonCandidate->Eta());
	}
	if(Photon->IsPrimary()){
		// Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
		if(fIsFromMBHeader){
			fCharPhotonMCInfo = 6;
			fHistoTruePrimaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			fHistoTruePrimaryConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt()); // Allways Filled
		}
		// (Not Filled for i6, Extra Signal Gamma (parambox) are secondary)
	} else {
		if(fIsFromMBHeader){
			fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			fCharPhotonMCInfo = 2;
			if(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother() > -1 &&
				((AliAODMCParticle*)AODMCTrackArray->At(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 3122){
				fCharPhotonMCInfo = 5;
				fHistoTrueSecondaryConvGammaFromXFromLambdaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			}
			if(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother() > -1 &&
				((AliAODMCParticle*)AODMCTrackArray->At(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 310){
				fCharPhotonMCInfo = 4;
				fHistoTrueSecondaryConvGammaFromXFromK0sPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			}
			if(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother() > -1 &&
				((AliAODMCParticle*)AODMCTrackArray->At(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 221){
				fCharPhotonMCInfo = 3;
			}
		}
	}
  
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessTruePhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{

	Double_t magField = fInputEvent->GetMagneticField();
	if( magField  < 0.0 ){
		magField =  1.0;
	}
	else {
		magField =  -1.0;
	}

	// Process True Photons
	TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(fMCStack);
	TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCStack);	
	fCharPhotonMCInfo = 0;
	
	if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
	Int_t pdgCode[2] = {abs(posDaughter->GetPdgCode()),abs(negDaughter->GetPdgCode())};
	fCharPhotonMCInfo = 1;
	if(posDaughter->GetMother(0) != negDaughter->GetMother(0)){
		FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode);
		return;
	}
	else if(posDaughter->GetMother(0) == -1){
		FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode);
		return;
	}
	
	if(pdgCode[0]!=11 || pdgCode[1]!=11) return; //One Particle is not a electron
	
	if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge

	TParticle *Photon = TruePhotonCandidate->GetMCParticle(fMCStack);
	
	if(Photon->GetPdgCode() != 22){
		return; // Mother is no Photon
	}
	
	if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() !=5) return;// check if the daughters come from a conversion
	
	// True Photon
	if(fIsFromMBHeader){
		fHistoTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
		if (fDoPhotonQA > 0) fHistoTrueConvGammaEta[fiCut]->Fill(TruePhotonCandidate->Eta());
	}
	if(posDaughter->GetMother(0) <= fMCStack->GetNprimary()){
		// Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
		if(fIsFromMBHeader){
			fCharPhotonMCInfo = 6;
			fHistoTruePrimaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			fHistoTruePrimaryConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt()); // Allways Filled
		
		}
		// (Not Filled for i6, Extra Signal Gamma (parambox) are secondary)
	}
	else{
		if(fIsFromMBHeader){
			fCharPhotonMCInfo = 2;
			fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if(fMCStack->Particle(Photon->GetMother(0))->GetMother(0) > -1 &&
				fMCStack->Particle(fMCStack->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 3122){
				fHistoTrueSecondaryConvGammaFromXFromLambdaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				fCharPhotonMCInfo = 5;
			}
			if(fMCStack->Particle(Photon->GetMother(0))->GetMother(0) > -1 &&
				fMCStack->Particle(fMCStack->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
				fHistoTrueSecondaryConvGammaFromXFromK0sPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				fCharPhotonMCInfo = 4;
			}
			if(fMCStack->Particle(Photon->GetMother(0))->GetMother(0) > -1 &&
				fMCStack->Particle(fMCStack->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 221){
				fCharPhotonMCInfo = 3;
			}
		}
	}

	// pi0 photon
	//Bool_t bpi0 = 0;
	Int_t imother = Photon->GetMother(0);
	if(imother > -1){
		AliMCParticle* McMother = static_cast<AliMCParticle*>(fMCEvent->GetTrack(imother));
		//cout << fMCEvent->GetRunNumber() << " " << imother << " " << fMCEvent->GetNumberOfTracks() << endl;
		if(McMother->PdgCode() == 111) fHistoTrueConvPi0GammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
	}
	return;
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessAODMCParticles()
{
	
	TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	
	// Loop over all primary MC particle
	for(Int_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {
		
		AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
		if (!particle) continue;
		if (!particle->IsPrimary()) continue;
		
		Int_t isMCFromMBHeader = -1;
		if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
			isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent);
			if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
		}
		
		if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) continue;
		if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,AODMCTrackArray,kFALSE)){
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
			if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,AODMCTrackArray,kTRUE)){
				Double_t rConv = 0;
				for(Int_t daughterIndex=particle->GetDaughter(0);daughterIndex<=particle->GetDaughter(1);daughterIndex++){
					AliAODMCParticle *tmpDaughter = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(daughterIndex));
					if(!tmpDaughter) continue;
					if(abs(tmpDaughter->GetPdgCode()) == 11){
						rConv = sqrt( (tmpDaughter->Xv()*tmpDaughter->Xv()) + (tmpDaughter->Yv()*tmpDaughter->Yv()) );
					}
				}
				fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt());
				if (fDoPhotonQA > 0){
					fHistoMCConvGammaR[fiCut]->Fill(rConv);
					fHistoMCConvGammaEta[fiCut]->Fill(particle->Eta());
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
						//                   if(particle->GetPdgCode() == 221){
						//                      cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
						//                   }
					}
				}
				Double_t mesonY = 10.;
				if(particle->E() - particle->Pz() == 0 || particle->E() + particle->Pz() == 0){
					mesonY=10.-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
				} else{
					mesonY = 0.5*(TMath::Log((particle->E()+particle->Pz()) / (particle->E()-particle->Pz())))-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
				}
				
				if(particle->GetPdgCode() == 111){
					fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted); // All MC Pi0
					fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt());
					if (fDoMesonQA > 0) fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
				} else if(particle->GetPdgCode() == 221){
					fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted); // All MC Eta
					fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt());
					if (fDoMesonQA > 0) fHistoMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
				}
				
				// Check the acceptance for both gammas
				if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter0,AODMCTrackArray,kFALSE) &&
				((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter1,AODMCTrackArray,kFALSE)  &&
				((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
				((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
					
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
void AliAnalysisTaskGammaConvCalo::ProcessMCParticles()
{
	// Loop over all primary MC particle
	for(Int_t i = 0; i < fMCStack->GetNprimary(); i++) {
		TParticle* particle = (TParticle *)fMCStack->Particle(i);
		if (!particle) continue;
		
		Int_t isMCFromMBHeader = -1;
		if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
			isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent);
			if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
		}
		
		if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) continue;
		if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCStack,kFALSE)){
			fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt()); // All MC Gamma
			if (abs(particle->Eta()) < 0.66 ){
				if (particle->Phi() > 1.39626 && particle->Phi() < 3.125) fHistoMCAllGammaEMCALAccPt[fiCut]->Fill(particle->Pt());
			}	
			
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
		if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCStack,kTRUE)){
			fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt());
			if (fDoPhotonQA > 0){
				fHistoMCConvGammaR[fiCut]->Fill(((TParticle*)fMCStack->Particle(particle->GetFirstDaughter()))->R());
				fHistoMCConvGammaEta[fiCut]->Fill(particle->Eta());
			}
		} // Converted MC Gamma
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
				
				if(particle->GetPdgCode() == 111){
					fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted); // All MC Pi0
					fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt());
					if (fDoMesonQA > 0) fHistoMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
				} else if(particle->GetPdgCode() == 221){
					fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted); // All MC Eta
					fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt());
					if (fDoMesonQA > 0) fHistoMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
				}
				
				// Check the acceptance for both gammas
				if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter0,fMCStack,kFALSE) &&
				((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter1,fMCStack,kFALSE)  &&
				((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
				((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
					
					if(particle->GetPdgCode() == 111){
						fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted); // MC Pi0 with gamma in acc
					} else if(particle->GetPdgCode() == 221){
						fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted); // MC Eta with gamma in acc
					}
				}
			}
		}
	}
  
	if (fDoMesonQA){
		for(Int_t i = fMCStack->GetNprimary(); i < fMCStack->GetNtrack(); i++) {
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

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::PhotonTagging(){
	
	// Conversion Gammas
	if(fGammaCandidates->GetEntries()>0){

		for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries()-1;firstGammaIndex++){
		
			// get conversion photon
			AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
			if (gamma0==NULL) continue;
			
			TLorentzVector photonVector;
			photonVector.SetPxPyPzE(gamma0->GetPx(),gamma0->GetPy(),gamma0->GetPz(),gamma0->GetPhotonP());
			
			Bool_t btagpi0 = 0;
			Bool_t btageta = 0;
			
			// loop over clusters
			for(Int_t secondGammaIndex = 0; secondGammaIndex<fClusterCandidates->GetEntries(); ++secondGammaIndex) {
				
				AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
				if (gamma1==NULL) continue;
				
				TLorentzVector clusterVector;
				clusterVector.SetPxPyPzE(gamma1->GetPx(),gamma1->GetPy(),gamma1->GetPz(),gamma1->GetPhotonP());
				
				// do the tagging
				TLorentzVector pairVector = photonVector+clusterVector;
				
				// see if pi0?
				if((pairVector.M() > 0.11 && pairVector.M() < 0.15)){
					btagpi0 = 1;	
				}
				// or eta
				if((pairVector.M() > 0.50 && pairVector.M() < 0.6)){
					btageta = 1;
				}
			}// end loop over clusters
			
			if(btagpi0 && btageta)
				fHistoConvGammaTagged[fiCut]->Fill(photonVector.Pt());
			else if(btagpi0 && !btageta)
				fHistoConvGammaPi0Tagged[fiCut]->Fill(photonVector.Pt());
			else if(btageta && !btagpi0)
				fHistoConvGammaEtaTagged[fiCut]->Fill(photonVector.Pt());
			else
				fHistoConvGammaUntagged[fiCut]->Fill(photonVector.Pt());

		}// end loop over gammas
	}// end if
	return;
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::CalculatePi0Candidates(){
	
	// Conversion Gammas
	if(fGammaCandidates->GetEntries()>0){

		// vertex
		Double_t vertex[3] = {0};
		InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

		for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
			AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
			if (gamma0==NULL) continue;
			
			for(Int_t secondGammaIndex=0;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){

				AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
				if (gamma1==NULL) continue;
				
				AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);
				pi0cand->SetLabels(firstGammaIndex,secondGammaIndex);
				
				if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
					fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
					
					// fill new histograms
					fHistoPhotonPairAll[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
					fHistoPhotonPairAllGam[fiCut]->Fill(pi0cand->M(),gamma0->Pt());
					
					if(pi0cand->GetAlpha()<0.1)
						fHistoMotherInvMassEalpha[fiCut]->Fill(pi0cand->M(),pi0cand->E());
					
					if (fDoMesonQA > 0){
						if ( pi0cand->M() > 0.05 && pi0cand->M() < 0.17){
							fHistoMotherPi0PtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
							fHistoMotherPi0PtAlpha[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetAlpha());
							fHistoMotherPi0PtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle());	
						}
						if ( pi0cand->M() > 0.45 && pi0cand->M() < 0.65){
							fHistoMotherEtaPtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
							fHistoMotherEtaPtAlpha[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetAlpha());
							fHistoMotherEtaPtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle());
						}
					}
					if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoBGCalculation()){
						Int_t zbin = 0;
						Int_t mbin = 0;
						
						if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->BackgroundHandlerType() == 0){
							zbin = fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
							if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
								mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
							} else {
								mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
							}
						} else{
							zbin = fBGHandlerRP[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
							if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
								mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
							} else {
								mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
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
					if (fDoMesonQA == 2){
						fInvMass = pi0cand->M();
						fPt  = pi0cand->Pt();
						if (abs(gamma0->GetDCAzToPrimVtx()) < abs(gamma1->GetDCAzToPrimVtx())){
							fDCAzGammaMin = gamma0->GetDCAzToPrimVtx();
							fDCAzGammaMax = gamma1->GetDCAzToPrimVtx();
						} else {
							fDCAzGammaMin = gamma1->GetDCAzToPrimVtx();
							fDCAzGammaMax = gamma0->GetDCAzToPrimVtx();
						}
						fCharFlag = pi0cand->GetMesonQuality();
						//                   cout << "gamma 0: " << gamma0->GetV0Index()<< "\t" << gamma0->GetPx() << "\t" << gamma0->GetPy() << "\t" <<  gamma0->GetPz() << "\t" << endl;
						//                   cout << "gamma 1: " << gamma1->GetV0Index()<< "\t"<< gamma1->GetPx() << "\t" << gamma1->GetPy() << "\t" <<  gamma1->GetPz() << "\t" << endl;
						//                    cout << "pi0: "<<fInvMass << "\t" << fPt <<"\t" << fDCAzGammaMin << "\t" << fDCAzGammaMax << "\t" << (Int_t)fCharFlag << "\t" << (Int_t)fCharMesonMCInfo <<endl;
						if (fIsHeavyIon == 1 && fPt > 0.399 && fPt < 20. ) {
							if (fInvMass > 0.08 && fInvMass < 0.2) fTreeMesonsInvMassPtDcazMinDcazMaxFlag[fiCut]->Fill();
							if ((fInvMass > 0.45 && fInvMass < 0.6) &&  (fPt > 0.999 && fPt < 20.) )fTreeMesonsInvMassPtDcazMinDcazMaxFlag[fiCut]->Fill();
						} else if (fPt > 0.299 && fPt < 20. )  {
							if ( (fInvMass > 0.08 && fInvMass < 0.2) || (fInvMass > 0.45 && fInvMass < 0.6)) fTreeMesonsInvMassPtDcazMinDcazMaxFlag[fiCut]->Fill();
						}
					}
					if (fDoMesonQA == 1){
						fHistoMotherInvMassECalib[fiCut]->Fill(pi0cand->M(),gamma1->E());
						if(pi0cand->GetAlpha()<0.1)
						fHistoMotherInvMassECalibalpha[fiCut]->Fill(pi0cand->M(),gamma1->E());            
					}

				}
				delete pi0cand;
				pi0cand=0x0;
			}
		}
	}
}
//______________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
	// Process True Mesons
	AliStack *MCStack = fMCEvent->Stack();
	fCharMesonMCInfo = 0;
	if(TrueGammaCandidate0->GetV0Index()<fInputEvent->GetNumberOfV0s()){
		Bool_t isTruePi0 = kFALSE;
		Bool_t isTrueEta = kFALSE;
		Int_t gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(MCStack);
		Int_t gamma0MotherLabel = -1;
		if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
			// Daughters Gamma 0
			TParticle * negativeMC = (TParticle*)TrueGammaCandidate0->GetNegativeMCDaughter(MCStack);
			TParticle * positiveMC = (TParticle*)TrueGammaCandidate0->GetPositiveMCDaughter(MCStack);
			TParticle * gammaMC0 = (TParticle*)MCStack->Particle(gamma0MCLabel);
			if(abs(negativeMC->GetPdgCode())==11 && abs(positiveMC->GetPdgCode())==11){  // Electrons ...
				if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
					if(gammaMC0->GetPdgCode() == 22){ // ... with Gamma Mother
						gamma0MotherLabel=gammaMC0->GetFirstMother();
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
			} else {
				if (fDoMesonQA > 0) fHistoTrueMotherCaloEMNonLeadingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
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
		
		if(isTruePi0 || isTrueEta){// True Pion or Eta
			fHistoTrueMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			if (fDoMesonQA > 0){
				if (TrueGammaCandidate1->IsLargestComponentPhoton()) 
					fHistoTrueMotherCaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (TrueGammaCandidate1->IsLargestComponentElectron()) 
					fHistoTrueMotherCaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()){
					fHistoTrueMotherCaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (isTruePi0)fHistoTruePi0CaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (isTrueEta)fHistoTrueEtaCaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if ((TrueGammaCandidate0->GetMCLabelPositive() == gamma1MCLabel || TrueGammaCandidate0->GetMCLabelNegative() == gamma1MCLabel) && isTruePi0)
						fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if ((TrueGammaCandidate0->GetMCLabelPositive() == gamma1MCLabel || TrueGammaCandidate0->GetMCLabelNegative() == gamma1MCLabel) && isTrueEta)
						fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				}	
				if (TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate1->IsDalitzMerged() )
					fHistoTrueMotherCaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (TrueGammaCandidate1->IsMergedPartConv()) 
					fHistoTrueMotherCaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
			if (fDoMesonQA > 0){
				if (isTruePi0){
					if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
						fHistoTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
						fHistoTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetAlpha());
						fHistoTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
					}
				} else if (isTrueEta){
					if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
						fHistoTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
						fHistoTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetAlpha());
						fHistoTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
					}
				}
			}
			if(gamma0MotherLabel >= MCStack->GetNprimary()){ // Secondary Meson
				Int_t secMotherLabel = ((TParticle*)MCStack->Particle(gamma0MotherLabel))->GetMother(0);
				Float_t weightedSec= 1;
				if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, fMCStack, fInputEvent) && MCStack->Particle(secMotherLabel)->GetPdgCode()==310){
					weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),secMotherLabel, fMCStack, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
					//cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
				}
				fHistoTrueSecondaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
				fCharMesonMCInfo = 2;
				if (secMotherLabel >-1){
					if(MCStack->Particle(secMotherLabel)->GetPdgCode()==310){
						fCharMesonMCInfo = 4;
						fHistoTrueSecondaryMotherFromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
						if (fDoMesonQA > 0)fHistoTrueK0sWithPi0DaughterMCPt[fiCut]->Fill(MCStack->Particle(secMotherLabel)->Pt());
					}
					if(MCStack->Particle(secMotherLabel)->GetPdgCode()==221){
						fCharMesonMCInfo = 3;
						fHistoTrueSecondaryMotherFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
						if (fDoMesonQA > 0)fHistoTrueEtaWithPi0DaughterMCPt[fiCut]->Fill(MCStack->Particle(secMotherLabel)->Pt());
					}
					if(MCStack->Particle(secMotherLabel)->GetPdgCode()==3122){
						fCharMesonMCInfo = 7;
						fHistoTrueSecondaryMotherFromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
						if (fDoMesonQA > 0)fHistoTrueLambdaWithPi0DaughterMCPt[fiCut]->Fill(MCStack->Particle(secMotherLabel)->Pt());
					}
				}
			} else { // Only primary pi0 for efficiency calculation
				fCharMesonMCInfo = 6;
				Float_t weighted= 1;
				if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, fMCStack, fInputEvent)){
					if (((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt()>0.005){
						weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),gamma1MotherLabel, fMCStack, fInputEvent);
						//                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
					}
				}
				fHistoTruePrimaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
				fHistoTruePrimaryMotherW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				fProfileTruePrimaryMotherWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
					
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
					fCharMesonMCInfo = 1;
				} else { // No photon or without mother
					fHistoTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				}
			}
		}
	}
}
//______________________________________________________________________
void AliAnalysisTaskGammaConvCalo::ProcessTrueMesonCandidatesAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
	
	// Process True Mesons
	TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	Bool_t isTruePi0 = kFALSE;
	Bool_t isTrueEta = kFALSE;
	
	AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
	AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));
	
	fCharMesonMCInfo = 0;
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
				}
			}
		}
	}	

	Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); 	// get most probable MC label
	Int_t gamma1MotherLabel = -1;
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
					AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gammaMC1->GetMother()));
					gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
				} else gamma1MotherLabel=gammaMC1->GetMother(); 
			}
		} else {
			if (fDoMesonQA > 0) fHistoTrueMotherCaloEMNonLeadingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
		}	
	}
			
	if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
		if(((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
			isTruePi0=kTRUE;
		}
		if(((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 221){
			isTrueEta=kTRUE;
		}
	}
	
	if(isTruePi0 || isTrueEta){// True Pion or Eta
		fHistoTrueMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
		if (fDoMesonQA > 0){
			if (TrueGammaCandidate1->IsLargestComponentPhoton()) 
				fHistoTrueMotherCaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			if (TrueGammaCandidate1->IsLargestComponentElectron()) 
				fHistoTrueMotherCaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			if (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()){
				fHistoTrueMotherCaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTruePi0)fHistoTruePi0CaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta)fHistoTrueEtaCaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if ((TrueGammaCandidate0->GetMCLabelPositive() == gamma1MCLabel || TrueGammaCandidate0->GetMCLabelNegative() == gamma1MCLabel) && isTruePi0)
					fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if ((TrueGammaCandidate0->GetMCLabelPositive() == gamma1MCLabel || TrueGammaCandidate0->GetMCLabelNegative() == gamma1MCLabel) && isTrueEta)
					fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	
			if (TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate1->IsDalitzMerged() )
				fHistoTrueMotherCaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			if (TrueGammaCandidate1->IsMergedPartConv()) 
				fHistoTrueMotherCaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
		}

		if (fDoMesonQA > 0){
			if (isTruePi0){
				if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
				fHistoTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
				fHistoTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetAlpha());
				fHistoTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
				}
			} else if (isTrueEta){
				if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
				fHistoTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
				fHistoTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetAlpha());
				fHistoTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
				}
			}
		}
		if(!(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel))->IsPrimary())){ // Secondary Meson
			Int_t secMotherLabel = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->GetMother();
			Float_t weightedSec= 1;
			if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, 0x0, fInputEvent) && static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
				weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),secMotherLabel, 0x0, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
				//cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
			}
			fHistoTrueSecondaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
			fCharMesonMCInfo = 2;
			if (secMotherLabel >-1){
				if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
					fCharMesonMCInfo = 4;
					fHistoTrueSecondaryMotherFromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
					if (fDoMesonQA > 0)fHistoTrueK0sWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->Pt());
				}
				if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==221){
					fCharMesonMCInfo = 3;
					fHistoTrueSecondaryMotherFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
					if (fDoMesonQA > 0)fHistoTrueEtaWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->Pt());
				}
				if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==3122){
					fCharMesonMCInfo = 7;
					fHistoTrueSecondaryMotherFromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
					if (fDoMesonQA > 0)fHistoTrueLambdaWithPi0DaughterMCPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->Pt());
				}
			}
		}else{ // Only primary pi0 for efficiency calculation
			Float_t weighted= 1;
			fCharMesonMCInfo = 6;
			if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, 0x0, fInputEvent)){
				if (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt()>0.005){
				weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),gamma1MotherLabel, 0x0, fInputEvent);
				//                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
				}
			}
			fHistoTruePrimaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
			fHistoTruePrimaryMotherW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			fProfileTruePrimaryMotherWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
			
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
				fCharMesonMCInfo = 1;
			} else { // No photon or without mother
				fHistoTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
		}
	}
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::CalculateBackground(){
	
	Int_t zbin= fBGClusHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
	Int_t mbin = 0;
	
	if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
		mbin = fBGClusHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
	} else {
		mbin = fBGClusHandler[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
	}
	
	AliGammaConversionAODBGHandler::GammaConversionVertex *bgEventVertex = NULL;    
	if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
		for(Int_t nEventsInBG=0;nEventsInBG<fBGClusHandler[fiCut]->GetNBGEvents();nEventsInBG++){
			AliGammaConversionAODVector *previousEventV0s = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
			if(fMoveParticleAccordingToVertex == kTRUE){
				bgEventVertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
			}
			
			for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){
				AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));
				for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
					AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
					if(fMoveParticleAccordingToVertex == kTRUE){
						MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
					}
					if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
						RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
					}
					
					AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
					backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
					if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
						->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
						fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
						Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
						fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
					}
					delete backgroundCandidate;
					backgroundCandidate = 0x0;
				}
			}
		}
	} else {
		for(Int_t nEventsInBG=0;nEventsInBG <fBGClusHandler[fiCut]->GetNBGEvents();nEventsInBG++){
			AliGammaConversionAODVector *previousEventV0s = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
			if(previousEventV0s){
				if(fMoveParticleAccordingToVertex == kTRUE){
					bgEventVertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
				}
				for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){
					AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));
					for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
				
						AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
					
						if(fMoveParticleAccordingToVertex == kTRUE){
							MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
						}
						if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
							RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
						}
					
						AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
						backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
						if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
							fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
							Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
							fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
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
void AliAnalysisTaskGammaConvCalo::CalculateBackgroundRP(){
	
	Int_t zbin= fBGHandlerRP[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
	Int_t mbin = 0;
	if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
		mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
	} else {
		mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
	}
	
	
	//Rotation Method
	if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseRotationMethod()){
		// Correct for the number of rotations
		// BG is for rotation the same, except for factor NRotations
		Double_t weight=1./Double_t(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfBGEvents());
		
		for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
			
			AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
			if (gamma0==NULL) continue;
			for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGammaCandidates->GetEntries();secondGammaIndex++){
				AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(secondGammaIndex));
				if (gamma1 == NULL) continue;
				if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelected(gamma1,fInputEvent))continue;
				for(Int_t nRandom=0;nRandom<((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfBGEvents();nRandom++){
				
					RotateParticle(gamma1);					
					AliAODConversionMother backgroundCandidate(gamma0,gamma1);
					backgroundCandidate.CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
					if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
						->MesonIsSelected(&backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
						fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(),backgroundCandidate.Pt());
						Double_t sparesFill[4] = {backgroundCandidate.M(),backgroundCandidate.Pt(),(Double_t)zbin,(Double_t)mbin};
						fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,weight);
					}
				}
			}
		}
	} else {
		// Do Event Mixing
		for(Int_t nEventsInBG=0;nEventsInBG <fBGHandlerRP[fiCut]->GetNBGEvents(fGammaCandidates,fInputEvent);nEventsInBG++){
		
			AliGammaConversionPhotonVector *previousEventGammas = fBGHandlerRP[fiCut]->GetBGGoodGammas(fGammaCandidates,fInputEvent,nEventsInBG);
			
			if(previousEventGammas){
				// test weighted background
				Double_t weight=1.0;
				// Correct for the number of eventmixing:
				// N gammas -> (N-1) + (N-2) +(N-3) ...+ (N-(N-1))  using sum formula sum(i)=N*(N-1)/2  -> N*(N-1)/2
				// real combinations (since you cannot combine a photon with its own)
				// but BG leads to N_{a}*N_{b} combinations
				weight*=0.5*(Double_t(fGammaCandidates->GetEntries()-1))/Double_t(previousEventGammas->size());
				
				for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){
					AliAODConversionPhoton *gamma0 = (AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));					
					for(UInt_t iPrevious=0;iPrevious<previousEventGammas->size();iPrevious++){
						
						AliAODConversionPhoton *gamma1 = (AliAODConversionPhoton*)(previousEventGammas->at(iPrevious));
						
						AliAODConversionMother backgroundCandidate(gamma0,gamma1);
						backgroundCandidate.CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
						if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
						->MesonIsSelected(&backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
						fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(),backgroundCandidate.Pt());
						Double_t sparesFill[4] = {backgroundCandidate.M(),backgroundCandidate.Pt(),(Double_t)zbin,(Double_t)mbin};
						fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,weight);
						}
					}
				}
			}
		}
	}
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::RotateParticle(AliAODConversionPhoton *gamma){
	Int_t fNDegreesPMBackground= ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->NDegreesRotation();
	Double_t nRadiansPM = fNDegreesPMBackground*TMath::Pi()/180;
	Double_t rotationValue = fRandom.Rndm()*2*nRadiansPM + TMath::Pi()-nRadiansPM;
	gamma->RotateZ(rotationValue);
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::RotateParticleAccordingToEP(AliAODConversionPhoton *gamma, Double_t previousEventEP, Double_t thisEventEP){
	
	previousEventEP=previousEventEP+TMath::Pi();
	thisEventEP=thisEventEP+TMath::Pi();
	Double_t rotationValue= thisEventEP-previousEventEP;
	gamma->RotateZ(rotationValue);
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex){
	//see header file for documentation
	
	Double_t dx = vertex->fX - fInputEvent->GetPrimaryVertex()->GetX();
	Double_t dy = vertex->fY - fInputEvent->GetPrimaryVertex()->GetY();
	Double_t dz = vertex->fZ - fInputEvent->GetPrimaryVertex()->GetZ();
	
	Double_t movedPlace[3] = {particle->GetConversionX() - dx,particle->GetConversionY() - dy,particle->GetConversionZ() - dz};
	particle->SetConversionPoint(movedPlace);
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::UpdateEventByEventData(){
	//see header file for documentation
	if(fGammaCandidates->GetEntries() >0 ){
		if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
			fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
			fBGClusHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
		} else { // means we use #V0s for multiplicity
			fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGammaCandidates->GetEntries(),fEventPlaneAngle);
			fBGClusHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGammaCandidates->GetEntries(),fEventPlaneAngle);
		}
	}
}


//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::FillPhotonCombinatorialBackgroundHist(AliAODConversionPhoton *TruePhotonCandidate, Int_t pdgCode[])
{
	// Combinatorial Bck = 0 ee, 1 ep,i 2 ek, 3 ep, 4 emu, 5 pipi, 6 pik, 7 pip, 8 pimu, 9 kk, 10 kp, 11 kmu, 12 pp, 13 pmu, 14 mumu, 15 Rest
	if(pdgCode[0]==11   && pdgCode[1]==11){if(fIsFromMBHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0);}
	else if( (pdgCode[0]==11   && pdgCode[1]==211) || (pdgCode[0]==211  && pdgCode[1]==11) )
	{if(fIsFromMBHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1);}
	else if( (pdgCode[0]==11   && pdgCode[1]==321) || (pdgCode[0]==321  && pdgCode[1]==11) )
	{if(fIsFromMBHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2);}
	else if( (pdgCode[0]==11   && pdgCode[1]==2212) || (pdgCode[0]==2212 && pdgCode[1]==11) )
	{if(fIsFromMBHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3);}
	else if( (pdgCode[0]==11   && pdgCode[1]==13) || (pdgCode[0]==13   && pdgCode[1]==11) )
	{if(fIsFromMBHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4);}
	else if(  pdgCode[0]==211  && pdgCode[1]==211 ){if(fIsFromMBHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),5);}
	else if( (pdgCode[0]==211  && pdgCode[1]==321) || (pdgCode[0]==321  && pdgCode[1]==211) )
	{if(fIsFromMBHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),6);}
	else if( (pdgCode[0]==211  && pdgCode[1]==2212) || (pdgCode[0]==2212 && pdgCode[1]==211) )
	{if(fIsFromMBHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),7);}
	else if( (pdgCode[0]==211  && pdgCode[1]==13) || (pdgCode[0]==13   && pdgCode[1]==211) )
	{if(fIsFromMBHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),8);}
	else if(  pdgCode[0]==321  && pdgCode[1]==321 ){if(fIsFromMBHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),9);}
	else if( (pdgCode[0]==321  && pdgCode[1]==2212) || (pdgCode[0]==2212 && pdgCode[1]==321) )
	{if(fIsFromMBHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),10);}
	else if( (pdgCode[0]==321  && pdgCode[1]==13) || (pdgCode[0]==13   && pdgCode[1]==321) )
	{if(fIsFromMBHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),11);}
	else if(  pdgCode[0]==2212   && pdgCode[1]==2212  ){if(fIsFromMBHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),12);}
	else if( (pdgCode[0]==2212  && pdgCode[1]==13) || (pdgCode[0]==13   && pdgCode[1]==2212) )
	{if(fIsFromMBHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),13);}
	else if(  pdgCode[0]==13   && pdgCode[1]==13  ){if(fIsFromMBHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),14);}
	else {if(fIsFromMBHeader)fHistoCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),15);}
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::RelabelAODPhotonCandidates(Bool_t mode){
	
	// Relabeling For AOD Event
	// ESDiD -> AODiD
	// MCLabel -> AODMCLabel
	
	if(mode){
		fMCStackPos = new Int_t[fReaderGammas->GetEntries()];
		fMCStackNeg = new Int_t[fReaderGammas->GetEntries()];
		fESDArrayPos = new Int_t[fReaderGammas->GetEntries()];
		fESDArrayNeg = new Int_t[fReaderGammas->GetEntries()];
	}
	
	for(Int_t iGamma = 0;iGamma<fReaderGammas->GetEntries();iGamma++){
		AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(iGamma);
		if(!PhotonCandidate) continue;
		if(!mode){// Back to ESD Labels
		PhotonCandidate->SetMCLabelPositive(fMCStackPos[iGamma]);
		PhotonCandidate->SetMCLabelNegative(fMCStackNeg[iGamma]);
		PhotonCandidate->SetLabelPositive(fESDArrayPos[iGamma]);
		PhotonCandidate->SetLabelNegative(fESDArrayNeg[iGamma]);
		continue;
		}
		fMCStackPos[iGamma] =  PhotonCandidate->GetMCLabelPositive();
		fMCStackNeg[iGamma] =  PhotonCandidate->GetMCLabelNegative();
		fESDArrayPos[iGamma] = PhotonCandidate->GetTrackLabelPositive();
		fESDArrayNeg[iGamma] = PhotonCandidate->GetTrackLabelNegative();
		
		Bool_t AODLabelPos = kFALSE;
		Bool_t AODLabelNeg = kFALSE;
		
		for(Int_t i = 0; i<fInputEvent->GetNumberOfTracks();i++){
		AliAODTrack *tempDaughter = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
		if(!AODLabelPos){
			if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelPositive() ){
			PhotonCandidate->SetMCLabelPositive(abs(tempDaughter->GetLabel()));
			PhotonCandidate->SetLabelPositive(i);
			AODLabelPos = kTRUE;
			}
		}
		if(!AODLabelNeg){
			if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelNegative()){
			PhotonCandidate->SetMCLabelNegative(abs(tempDaughter->GetLabel()));
			PhotonCandidate->SetLabelNegative(i);
			AODLabelNeg = kTRUE;
			}
		}
		if(AODLabelNeg && AODLabelPos){
			break;
		}
		}
		if(!AODLabelPos || !AODLabelNeg){
		cout<<"WARNING!!! AOD TRACKS NOT FOUND FOR"<<endl;
		}
	}
	
	
	if(!mode){
		delete[] fMCStackPos;
		delete[] fMCStackNeg;
		delete[] fESDArrayPos;
		delete[] fESDArrayNeg;
	}
}

void AliAnalysisTaskGammaConvCalo::SetLogBinningXTH2(TH2* histoRebin){
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

//________________________________________________________________________
void AliAnalysisTaskGammaConvCalo::Terminate(const Option_t *)
{
  
  //fOutputContainer->Print(); // Will crash on GRID
}

//________________________________________________________________________
Int_t AliAnalysisTaskGammaConvCalo::GetSourceClassification(Int_t daughter, Int_t pdgCode){
  
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

// //________________________________________________________________________
// Double_t AliAnalysisTaskGammaConvCalo::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
// {
//   // Get maximum energy of attached cell.
//   
//   id = -1;
//   Double_t maxe = 0;
//   Int_t ncells = cluster->GetNCells();
//   if (fEsdCells) {
//     for (Int_t i=0; i<ncells; i++) {
//       Double_t e = fEsdCells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
//       if (e>maxe) {
//         maxe = e;
//         id   = cluster->GetCellAbsId(i);
//       }
//     }
//   } else {
//     for (Int_t i=0; i<ncells; i++) {
//       Double_t e = fAodCells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
//       if (e>maxe)
//         maxe = e;
//       id   = cluster->GetCellAbsId(i);
//     }
//   }
//   return maxe;
// }
