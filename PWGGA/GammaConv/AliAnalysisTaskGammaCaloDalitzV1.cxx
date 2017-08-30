/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *									 									  *
 * Author: Pedro Gonzalez, Baldo Sahlmueller, Friederike Bock				       		  *
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
#include "TDatabasePDG.h"
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
#include "AliAnalysisTaskGammaCaloDalitzV1.h"
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
#include "AliCaloTrackMatcher.h"
#include <vector>

ClassImp(AliAnalysisTaskGammaCaloDalitzV1)

//________________________________________________________________________
AliAnalysisTaskGammaCaloDalitzV1::AliAnalysisTaskGammaCaloDalitzV1(): AliAnalysisTaskSE(),
	fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
	fElecSelector(NULL),
	fBGClusHandler(NULL),
	fInputEvent(NULL),
	fMCEvent(NULL),
	fCutFolder(NULL),
	fESDList(NULL),
	fBackList(NULL),
	fMotherList(NULL),
	fPhotonDCAList(NULL),
	fTrueList(NULL),
	fMCList(NULL),
	fHeaderNameList(NULL),
	fClusterOutputList(NULL),
	fOutputContainer(NULL),
	fReaderGammas(NULL),
	fSelectorElectronIndex(0),
	fSelectorPositronIndex(0),
	fGammaCandidates(NULL),
	fClusterCandidates(NULL),
	fVirtualGammaCandidates(NULL),
	fEventCutArray(NULL),
	fEventCuts(NULL),
	fGammaCutArray(NULL),
	fElectronCutArray(NULL),
	fConversionCuts(NULL),
	fClusterCutArray(NULL),
	fCaloPhotonCuts(NULL),
	fMesonCutArray(NULL),
	fMesonCuts(NULL),
	fHistoConvGammaPt(NULL),
	fHistoConvGammaR(NULL),
	fHistoConvGammaEta(NULL),
	fHistoDalitzElectronPt(NULL),
	fHistoDalitzPositronPt(NULL),
	fHistoDalitzElectronPhi(NULL),
	fHistoDalitzPositronPhi(NULL),
	fPtGamma(0),
	fDCAzPhoton(0),
	fRConvPhoton(0),
	fEtaPhoton(0),
	fCharCatPhoton(0),
	fCharPhotonMCInfo(0),
	fHistoMotherInvMassPt(NULL),
	fHistoMotherInvMassOpeningAngleGammaElectron(NULL),
	fHistoMotherMatchedInvMassPt(NULL),
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
	fHistoMotherPi0ConvPhotonEtaPhi(NULL),
	fHistoMotherEtaConvPhotonEtaPhi(NULL),
	fHistoMotherInvMassECalib(NULL),
	fHistoMotherInvMassECalibalpha(NULL),
	fHistoClusGammaPt(NULL),
	fHistoClusOverlapHeadersGammaPt(NULL),
	fHistoMCHeaders(NULL),
	fHistoMCAllGammaPt(NULL),
	fHistoMCAllGammaPi0Pt(NULL),
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
	fHistoMCAllPositronsPt(NULL),
	fHistoMCDecayPositronPi0Pt(NULL),
	fHistoMCAllElectronsPt(NULL),
	fHistoMCDecayElectronPi0Pt(NULL),
	fHistoMCDecayNoPrimElectronPi0DalitzR(NULL),
	fHistoMCDecayNoPrimPositronPi0DalitzR(NULL),
	fHistoMCDecayNoPrimElectronPi0DalitzID(NULL),
	fHistoMCDecayNoPrimPositronPi0DalitzID(NULL),
	fHistoMCPi0GGPt(NULL),
	fHistoMCPi0GGWOWeightPt(NULL),
	fHistoMCPi0Pt(NULL),
	fHistoMCPi0WOWeightPt(NULL),
	fHistoMCEtaGGPt(NULL),
	fHistoMCEtaGGWOWeightPt(NULL),
	fHistoMCEtaPt(NULL),
	fHistoMCEtaWOWeightPt(NULL),
	fHistoMCPi0InAccPt(NULL),
	fHistoMCPi0InAccOpeningAngleGammaElectron(NULL),
	fHistoMCEtaInAccPt(NULL),
	fHistoMCPi0PtY(NULL),
	fHistoMCEtaPtY(NULL),
	fHistoMCPi0PtAlpha(NULL),
	fHistoMCEtaPtAlpha(NULL),
	fHistoMCK0sPt(NULL),
	fHistoMCK0sWOWeightPt(NULL),
	fHistoMCK0sPtY(NULL),
	fHistoTruePi0InvMassPt(NULL),
	fHistoTrueEtaInvMassPt(NULL),
	fHistoTruePi0ShowerInvMassPt(NULL),
	fHistoTrueEtaShowerInvMassPt(NULL),
	fHistoTruePi0NoShowerInvMassPt(NULL),
	fHistoTrueEtaNoShowerInvMassPt(NULL),
	fHistoTruePi0OpeningAngleGammaElectron(NULL),
	fHistoTruePi0GGInvMassPt(NULL),
	fHistoTrueEtaGGInvMassPt(NULL),
	fHistoTruePi0CaloPhotonInvMassPt(NULL),
	fHistoTrueEtaCaloPhotonInvMassPt(NULL),
	fHistoTruePi0CaloConvertedPhotonInvMassPt(NULL),
	fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt(NULL),
	fHistoTrueEtaCaloConvertedPhotonInvMassPt(NULL),
	fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt(NULL),
	fHistoTruePi0CaloElectronInvMassPt(NULL),
	fHistoTrueEtaCaloElectronInvMassPt(NULL),
	fHistoTruePi0CaloMergedClusterInvMassPt(NULL),
	fHistoTrueEtaCaloMergedClusterInvMassPt(NULL),
	fHistoTrueMotherCaloEMNonLeadingInvMassPt(NULL),
	fHistoTruePi0CaloMergedClusterPartConvInvMassPt(NULL),
	fHistoTrueEtaCaloMergedClusterPartConvInvMassPt(NULL),
	fHistoTruePrimaryPi0InvMassPt(NULL),
	fHistoTruePrimaryPi0GGInvMassPt(NULL),
	fHistoTruePrimaryEtaInvMassPt(NULL),
	fHistoTruePrimaryEtaGGInvMassPt(NULL),
	fHistoTruePrimaryPi0W0WeightingInvMassPt(NULL),
	fHistoTruePrimaryEtaW0WeightingInvMassPt(NULL),
	fProfileTruePrimaryPi0WeightsInvMassPt(NULL),
	fProfileTruePrimaryEtaWeightsInvMassPt(NULL),
	fHistoTruePrimaryPi0MCPtResolPt(NULL),
	fHistoTruePrimaryEtaMCPtResolPt(NULL),
	fHistoTrueMotherPi0ConvPhotonEtaPhi(NULL),
	fHistoTrueMotherEtaConvPhotonEtaPhi(NULL),
	fHistoTrueSecondaryPi0InvMassPt(NULL),
	fHistoTrueSecondaryPi0GGInvMassPt(NULL),
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
	fHistoTrueConvGammaPt(NULL),
	fHistoTrueConvPi0GammaPt(NULL),
	fHistoTrueConvGammaEta(NULL),
	fHistoTruePositronPt(NULL),
	fHistoTrueElectronPt(NULL),
	fHistoTrueSecPositronPt(NULL),
	fHistoTrueSecElectronPt(NULL),								
	fHistoTruePi0DalitzPositronPt(NULL),
	fHistoTruePi0DalitzElectronPt(NULL),
	fHistoTruePi0DalitzSecPositronPt(NULL),
	fHistoTruePi0DalitzSecElectronPt(NULL),
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
	fHistoTruePi0DalitzClusGammaPt(NULL),
	fHistoTruePi0DalitzAllClusGammaPt(NULL),
	fHistoTruePi0DalitzClusGammaMCPt(NULL),
	fHistoTruePrimaryPi0PhotonPairPtconv(NULL),
	fHistoTruePrimaryPi0DCPtconv(NULL),
	fHistoTruePrimaryPi0MissingPtconv(NULL),
	fHistoTruePrimaryEtaPhotonPairPtconv(NULL),
	fHistoTruePrimaryEtaDCPtconv(NULL),
	fHistoTruePrimaryEtaMissingPtconv(NULL),
	fHistoTrueSecondaryPi0PhotonPairPtconv(NULL),
	fHistoTrueSecondaryPi0DCPtconv(NULL),
	fHistoTrueSecondaryPi0MissingPtconv(NULL),
	fStringRecTruePi0s(NULL),
	fStringRecTrueEtas(NULL),
	fHistoDoubleCountTruePi0InvMassPt(NULL),
	fHistoDoubleCountTrueEtaInvMassPt(NULL),
	fHistoDoubleCountTrueConvGammaRPt(NULL),
	fVectorDoubleCountTruePi0s(0),
	fVectorDoubleCountTrueEtas(0),
	fVectorDoubleCountTrueConvGammas(0),
	fHistoNEvents(NULL),
	fHistoNGoodESDTracks(NULL),
	fHistoNGammaCandidates(NULL),
	fHistoNGoodESDTracksVsNGammaCanditates(NULL),
	fHistoSPDClusterTrackletBackground(NULL),
	fHistoNV0Tracks(NULL),
	fProfileEtaShift(NULL),
	fEventPlaneAngle(-100),
	fRandom(0),
	fNGammaCandidates(0),
	fUnsmearedPx(NULL),
	fUnsmearedPy(NULL),
	fUnsmearedPz(NULL),
	fUnsmearedE(NULL),
    fMCEventPos(NULL),
    fMCEventNeg(NULL),
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
	fIsOverlappingWithOtherHeader(kFALSE),
	fIsMC(kFALSE),
	fDoTHnSparse(kTRUE)
{
  
}

//________________________________________________________________________
AliAnalysisTaskGammaCaloDalitzV1::AliAnalysisTaskGammaCaloDalitzV1(const char *name):
	AliAnalysisTaskSE(name),
	fV0Reader(NULL),
    fV0ReaderName("V0ReaderV1"),
	fElecSelector(NULL),
	fBGClusHandler(NULL),
	fInputEvent(NULL),
	fMCEvent(NULL),
	fCutFolder(NULL),
	fESDList(NULL),
	fBackList(NULL),
	fMotherList(NULL),
	fPhotonDCAList(NULL),
	fTrueList(NULL),
	fMCList(NULL),
	fHeaderNameList(NULL),
	fClusterOutputList(NULL),
	fOutputContainer(0),
	fReaderGammas(NULL),
	fSelectorElectronIndex(0),
	fSelectorPositronIndex(0),
	fGammaCandidates(NULL),
	fClusterCandidates(NULL),
	fVirtualGammaCandidates(NULL),
	fEventCutArray(NULL),
	fEventCuts(NULL),
	fGammaCutArray(NULL),
	fElectronCutArray(NULL),
	fConversionCuts(NULL),
	fClusterCutArray(NULL),
	fCaloPhotonCuts(NULL),
	fMesonCutArray(NULL),
	fMesonCuts(NULL),
	fHistoConvGammaPt(NULL),
	fHistoConvGammaR(NULL),
	fHistoConvGammaEta(NULL),
	fHistoDalitzElectronPt(NULL),
	fHistoDalitzPositronPt(NULL),
	fHistoDalitzElectronPhi(NULL),
	fHistoDalitzPositronPhi(NULL),
	fPtGamma(0),
	fDCAzPhoton(0),
	fRConvPhoton(0),
	fEtaPhoton(0),
	fCharCatPhoton(0),
	fCharPhotonMCInfo(0),
	fHistoMotherInvMassPt(NULL),
	fHistoMotherInvMassOpeningAngleGammaElectron(NULL),
	fHistoMotherMatchedInvMassPt(NULL),
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
	fHistoMotherPi0ConvPhotonEtaPhi(NULL),
	fHistoMotherEtaConvPhotonEtaPhi(NULL),
	fHistoMotherInvMassECalib(NULL),
	fHistoMotherInvMassECalibalpha(NULL),
	fHistoClusGammaPt(NULL),
	fHistoClusOverlapHeadersGammaPt(NULL),
	fHistoMCHeaders(NULL),
	fHistoMCAllGammaPt(NULL),
	fHistoMCAllGammaPi0Pt(NULL),
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
	fHistoMCAllPositronsPt(NULL),
	fHistoMCDecayPositronPi0Pt(NULL),
	fHistoMCAllElectronsPt(NULL),
	fHistoMCDecayElectronPi0Pt(NULL),
	fHistoMCDecayNoPrimElectronPi0DalitzR(NULL),
	fHistoMCDecayNoPrimPositronPi0DalitzR(NULL),
	fHistoMCDecayNoPrimElectronPi0DalitzID(NULL),
	fHistoMCDecayNoPrimPositronPi0DalitzID(NULL),
	fHistoMCPi0GGPt(NULL),
	fHistoMCPi0GGWOWeightPt(NULL),
	fHistoMCPi0Pt(NULL),
	fHistoMCPi0WOWeightPt(NULL),
	fHistoMCEtaGGPt(NULL),
	fHistoMCEtaGGWOWeightPt(NULL),
	fHistoMCEtaPt(NULL),
	fHistoMCEtaWOWeightPt(NULL),
	fHistoMCPi0InAccPt(NULL),
	fHistoMCPi0InAccOpeningAngleGammaElectron(NULL),
	fHistoMCEtaInAccPt(NULL),
	fHistoMCPi0PtY(NULL),
	fHistoMCEtaPtY(NULL),
	fHistoMCPi0PtAlpha(NULL),
	fHistoMCEtaPtAlpha(NULL),
	fHistoMCK0sPt(NULL),
	fHistoMCK0sWOWeightPt(NULL),
	fHistoMCK0sPtY(NULL),
	fHistoTruePi0InvMassPt(NULL),
	fHistoTrueEtaInvMassPt(NULL),
	fHistoTruePi0ShowerInvMassPt(NULL),
	fHistoTrueEtaShowerInvMassPt(NULL),
	fHistoTruePi0NoShowerInvMassPt(NULL),
	fHistoTrueEtaNoShowerInvMassPt(NULL),
	fHistoTruePi0OpeningAngleGammaElectron(NULL),
	fHistoTruePi0GGInvMassPt(NULL),
	fHistoTrueEtaGGInvMassPt(NULL),
	fHistoTruePi0CaloPhotonInvMassPt(NULL),
	fHistoTrueEtaCaloPhotonInvMassPt(NULL),
	fHistoTruePi0CaloConvertedPhotonInvMassPt(NULL),
	fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt(NULL),
	fHistoTrueEtaCaloConvertedPhotonInvMassPt(NULL),
	fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt(NULL),
	fHistoTruePi0CaloElectronInvMassPt(NULL),
	fHistoTrueEtaCaloElectronInvMassPt(NULL),
	fHistoTruePi0CaloMergedClusterInvMassPt(NULL),
	fHistoTrueEtaCaloMergedClusterInvMassPt(NULL),
	fHistoTrueMotherCaloEMNonLeadingInvMassPt(NULL),
	fHistoTruePi0CaloMergedClusterPartConvInvMassPt(NULL),
	fHistoTrueEtaCaloMergedClusterPartConvInvMassPt(NULL),
	fHistoTruePrimaryPi0InvMassPt(NULL),
	fHistoTruePrimaryPi0GGInvMassPt(NULL),
	fHistoTruePrimaryEtaInvMassPt(NULL),
	fHistoTruePrimaryEtaGGInvMassPt(NULL),
	fHistoTruePrimaryPi0W0WeightingInvMassPt(NULL),
	fHistoTruePrimaryEtaW0WeightingInvMassPt(NULL),
	fProfileTruePrimaryPi0WeightsInvMassPt(NULL),
	fProfileTruePrimaryEtaWeightsInvMassPt(NULL),
	fHistoTruePrimaryPi0MCPtResolPt(NULL),
	fHistoTruePrimaryEtaMCPtResolPt(NULL),
	fHistoTrueMotherPi0ConvPhotonEtaPhi(NULL),
	fHistoTrueMotherEtaConvPhotonEtaPhi(NULL),
	fHistoTrueSecondaryPi0InvMassPt(NULL),
	fHistoTrueSecondaryPi0GGInvMassPt(NULL),
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
	fHistoTrueConvGammaPt(NULL),
	fHistoTrueConvPi0GammaPt(NULL),
	fHistoTrueConvGammaEta(NULL),
	fHistoTruePositronPt(NULL),
	fHistoTrueElectronPt(NULL),
	fHistoTrueSecPositronPt(NULL),
	fHistoTrueSecElectronPt(NULL),
	fHistoTruePi0DalitzPositronPt(NULL),
	fHistoTruePi0DalitzElectronPt(NULL),
	fHistoTruePi0DalitzSecPositronPt(NULL),
	fHistoTruePi0DalitzSecElectronPt(NULL),
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
	fHistoTruePi0DalitzClusGammaPt(NULL),
	fHistoTruePi0DalitzAllClusGammaPt(NULL),
	fHistoTruePi0DalitzClusGammaMCPt(NULL),
	fHistoTruePrimaryPi0PhotonPairPtconv(NULL),
	fHistoTruePrimaryPi0DCPtconv(NULL),
	fHistoTruePrimaryPi0MissingPtconv(NULL),
	fHistoTruePrimaryEtaPhotonPairPtconv(NULL),
	fHistoTruePrimaryEtaDCPtconv(NULL),
	fHistoTruePrimaryEtaMissingPtconv(NULL),
	fHistoTrueSecondaryPi0PhotonPairPtconv(NULL),
	fHistoTrueSecondaryPi0DCPtconv(NULL),
	fHistoTrueSecondaryPi0MissingPtconv(NULL),
	fStringRecTruePi0s(NULL),
	fStringRecTrueEtas(NULL),
	fHistoDoubleCountTruePi0InvMassPt(NULL),
	fHistoDoubleCountTrueEtaInvMassPt(NULL),
	fHistoDoubleCountTrueConvGammaRPt(NULL),
	fVectorDoubleCountTruePi0s(0),
	fVectorDoubleCountTrueEtas(0),
	fVectorDoubleCountTrueConvGammas(0),
	fHistoNEvents(NULL),
	fHistoNGoodESDTracks(NULL),
	fHistoNGammaCandidates(NULL),
	fHistoNGoodESDTracksVsNGammaCanditates(NULL),
	fHistoSPDClusterTrackletBackground(NULL),
	fHistoNV0Tracks(NULL),
	fProfileEtaShift(NULL),
	fEventPlaneAngle(-100),
	fRandom(0),
	fNGammaCandidates(0),
	fUnsmearedPx(NULL),
	fUnsmearedPy(NULL),
	fUnsmearedPz(NULL),
	fUnsmearedE(NULL),
    fMCEventPos(NULL),
    fMCEventNeg(NULL),
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
	fIsOverlappingWithOtherHeader(kFALSE),
	fIsMC(kFALSE),
	fDoTHnSparse(kTRUE)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaCaloDalitzV1::~AliAnalysisTaskGammaCaloDalitzV1()
{
	if(fGammaCandidates){
		delete fGammaCandidates;
		fGammaCandidates = 0x0;
	}
	if(fClusterCandidates){
		delete fClusterCandidates;
		fClusterCandidates = 0x0;
	}
	if(fVirtualGammaCandidates){
		delete fVirtualGammaCandidates;
		fVirtualGammaCandidates = 0x0;
	}
		
	if(fBGClusHandler){
		delete[] fBGClusHandler;
		fBGClusHandler = 0x0;
	}
	
}
//___________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::InitBack(){
	
	const Int_t nDim = 4;
	Int_t nBins[nDim] = {800,250,7,4};
	Double_t xMin[nDim] = {0,0, 0,0};
	Double_t xMax[nDim] = {0.8,25,7,4};
	
	if( fDoTHnSparse ) {
	
	fSparseMotherInvMassPtZM = new THnSparseF*[fnCuts];
	fSparseMotherBackInvMassPtZM = new THnSparseF*[fnCuts];
	
	}
	
	
	fBGClusHandler = new AliGammaConversionAODBGHandler*[fnCuts];
	
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
	  
		if (((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
			TString cutstringEvent 	  = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
			TString cutstringPhoton	  = ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutNumber();
			TString cutstringElectron = ((AliDalitzElectronCuts*)fElectronCutArray->At(iCut))->GetCutNumber();
			TString cutstringCalo 	  = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
			TString cutstringMeson 	  = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
			
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
			
			if( fDoTHnSparse ) {
			
			fBackList[iCut] = new TList();				    
			fBackList[iCut]->SetName(Form("%s_%s_%s_%s_%s Back histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringElectron.Data(),cutstringCalo.Data(), cutstringMeson.Data()));
			fBackList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fBackList[iCut]);
			
			fSparseMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_m","Back_Back_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
			fBackList[iCut]->Add(fSparseMotherBackInvMassPtZM[iCut]);
			
			fMotherList[iCut] = new TList();
			fMotherList[iCut]->SetName(Form("%s_%s_%s_%s_%s Mother histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringElectron.Data(),cutstringCalo.Data(), cutstringMeson.Data()));
			fMotherList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fMotherList[iCut]);
			
			fSparseMotherInvMassPtZM[iCut] = new THnSparseF("Back_Mother_InvMass_Pt_z_m","Back_Mother_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
			fMotherList[iCut]->Add(fSparseMotherInvMassPtZM[iCut]);
			
			}
			
			if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
			  
			  
				fBGClusHandler[iCut] = new AliGammaConversionAODBGHandler(
																	collisionSystem,centMin,centMax,
																	((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
																	((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity(),
																	2,8,7);
				
			} 
		}
	}
}
//________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::UserCreateOutputObjects(){
  
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
	fClusterCandidates->SetOwner(kTRUE);
	fVirtualGammaCandidates = new TList();
	fVirtualGammaCandidates->SetOwner(kTRUE);
  
	fCutFolder 					= new TList*[fnCuts];
	fESDList 					= new TList*[fnCuts];
	fBackList 					= new TList*[fnCuts];
	fMotherList 					= new TList*[fnCuts];
	fHistoNEvents 					= new TH1I*[fnCuts];
	fHistoNGoodESDTracks 				= new TH1I*[fnCuts];
	fHistoNGammaCandidates 				= new TH1I*[fnCuts];
	fHistoNGoodESDTracksVsNGammaCanditates 		= new TH2F*[fnCuts];
	fHistoSPDClusterTrackletBackground = new TH2F*[fnCuts];
	fHistoNV0Tracks 				= new TH1I*[fnCuts];
	fProfileEtaShift 				= new TProfile*[fnCuts];
	fHistoConvGammaPt 				= new TH1F*[fnCuts];
  
	if (fDoPhotonQA == 2){
		fPhotonDCAList = new TList*[fnCuts];
		
	}
	if (fDoPhotonQA > 0){
		fHistoConvGammaR = new TH1F*[fnCuts];
		fHistoConvGammaEta = new TH1F*[fnCuts];
	}
	
	fHistoDalitzElectronPt		= new TH1F*[fnCuts];
	fHistoDalitzPositronPt		= new TH1F*[fnCuts];
	fHistoDalitzElectronPhi		= new TH1F*[fnCuts];
	fHistoDalitzPositronPhi		= new TH1F*[fnCuts];
	
	
	if(fDoMesonAnalysis){
		fHistoMotherInvMassPt 			= new TH2F*[fnCuts];
		fHistoMotherMatchedInvMassPt 		= new TH2F*[fnCuts];
		fHistoMotherBackInvMassPt 		= new TH2F*[fnCuts];
		fHistoMotherInvMassEalpha 		= new TH2F*[fnCuts];
		if (fDoMesonQA > 0){
			fHistoMotherPi0PtY 				=  new TH2F*[fnCuts];
			fHistoMotherEtaPtY 				=  new TH2F*[fnCuts];
			fHistoMotherPi0PtAlpha 				=  new TH2F*[fnCuts];
			fHistoMotherEtaPtAlpha 				=  new TH2F*[fnCuts];
			fHistoMotherPi0PtOpenAngle			=  new TH2F*[fnCuts];
			fHistoMotherEtaPtOpenAngle 			=  new TH2F*[fnCuts];
			fHistoMotherPi0ConvPhotonEtaPhi 		=  new TH2F*[fnCuts];
			fHistoMotherEtaConvPhotonEtaPhi 		=  new TH2F*[fnCuts];
			fHistoMotherInvMassOpeningAngleGammaElectron 	=  new TH1F*[fnCuts];
		}
		if(fDoMesonQA == 1){
			fHistoMotherInvMassECalib	 	= new TH2F*[fnCuts];
			fHistoMotherInvMassECalibalpha 		= new TH2F*[fnCuts];
		}
	}
	
	fClusterOutputList = new TList*[fnCuts];
	fHistoClusGammaPt  = new TH1F*[fnCuts];
	fHistoClusOverlapHeadersGammaPt = new TH1F*[fnCuts];

	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		TString cutstringEvent 	  = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
		TString cutstringPhoton   = ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutNumber();
		TString cutstringElectron = ((AliDalitzElectronCuts*)fElectronCutArray->At(iCut))->GetCutNumber();
		TString cutstringCalo 	  = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
		TString cutstringMeson 	  = "NoMesonCut";
		if(fDoMesonAnalysis)cutstringMeson = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
    
		fCutFolder[iCut] = new TList();
		fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s_%s_%s",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringElectron.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
		fCutFolder[iCut]->SetOwner(kTRUE);
		fOutputContainer->Add(fCutFolder[iCut]);
		fESDList[iCut] = new TList();
		fESDList[iCut]->SetName(Form("%s_%s_%s_%s_%s ESD histograms", cutstringEvent.Data(),cutstringPhoton.Data(),cutstringElectron.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
		fESDList[iCut]->SetOwner(kTRUE);
		fCutFolder[iCut]->Add(fESDList[iCut]);
    
		fHistoNEvents[iCut] = new TH1I("NEvents","NEvents",14,-0.5,13.5);
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
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problem");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
		fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
		fESDList[iCut]->Add(fHistoNEvents[iCut]);
		
		if(fIsHeavyIon == 1) fHistoNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",4000,0,4000);
		else if(fIsHeavyIon == 2) fHistoNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",400,0,400);
		else fHistoNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",200,0,200);
		fHistoNGoodESDTracks[iCut]->SetXTitle("# TPC tracks");
		fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);
		
		if(fIsHeavyIon == 1) fHistoNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",100,0,100);
		else if(fIsHeavyIon == 2) fHistoNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",50,0,50);
		else fHistoNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",50,0,50);
		fHistoNGammaCandidates[iCut]->SetXTitle("# accepted $#gamma_{conv}");
		fESDList[iCut]->Add(fHistoNGammaCandidates[iCut]);
		
		if(fIsHeavyIon == 1) fHistoNGoodESDTracksVsNGammaCanditates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",4000,0,4000,100,0,100);
		else if(fIsHeavyIon == 2) fHistoNGoodESDTracksVsNGammaCanditates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",400,0,400,50,0,50);
		else fHistoNGoodESDTracksVsNGammaCanditates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",200,0,200,50,0,50);
		fHistoNGoodESDTracksVsNGammaCanditates[iCut]->SetXTitle("# TPC tracks");
		fHistoNGoodESDTracksVsNGammaCanditates[iCut]->SetYTitle("# accepted $#gamma_{calo}");
		fESDList[iCut]->Add(fHistoNGoodESDTracksVsNGammaCanditates[iCut]);

		fHistoSPDClusterTrackletBackground[iCut] = new TH2F("SPD tracklets vs SPD clusters","SPD tracklets vs SPD clusters",100,0,200,250,0,1000);
		fESDList[iCut]->Add(fHistoSPDClusterTrackletBackground[iCut]);
		
		if(fIsHeavyIon == 1) fHistoNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",30000,0,30000);
		else if(fIsHeavyIon == 2) fHistoNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",2500,0,2500);
		else fHistoNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",1500,0,1500);
		fHistoNV0Tracks[iCut]->SetXTitle("VZERO amp [arb. units]");
		fESDList[iCut]->Add(fHistoNV0Tracks[iCut]);
		fProfileEtaShift[iCut] = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
		fESDList[iCut]->Add(fProfileEtaShift[iCut]);
		fHistoConvGammaPt[iCut] = new TH1F("ESD_ConvGamma_Pt","ESD_ConvGamma_Pt",250,0,25);
		fHistoConvGammaPt[iCut]->SetXTitle("p_{T,conv} (GeV/c)");
		fESDList[iCut]->Add(fHistoConvGammaPt[iCut]);
    
		if (fDoPhotonQA == 2){
			fPhotonDCAList[iCut] = new TList();
			fPhotonDCAList[iCut]->SetName(Form("%s_%s_%s_%s Photon DCA tree",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
			fPhotonDCAList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fPhotonDCAList[iCut]);
      
			
		}
    
		if (fDoPhotonQA > 0){
			fHistoConvGammaR[iCut] = new TH1F("ESD_ConvGamma_R","ESD_ConvGamma_R",800,0,200);
			fHistoConvGammaR[iCut]->SetXTitle("R_{conv} (cm)");
			fESDList[iCut]->Add(fHistoConvGammaR[iCut]);
			fHistoConvGammaEta[iCut] = new TH1F("ESD_ConvGamma_Eta","ESD_ConvGamma_Eta",2000,-2,2);
			fHistoConvGammaEta[iCut]->SetXTitle("#eta_{conv}");
			fESDList[iCut]->Add(fHistoConvGammaEta[iCut]);
		}
		
		
		
		
		
		fHistoDalitzElectronPt[iCut] = new TH1F("ESD_DalitzElectron_Pt","ESD_DalitzElectron_Pt",1000,0,25);
		fESDList[iCut]->Add(fHistoDalitzElectronPt[iCut]);

		fHistoDalitzPositronPt[iCut] = new TH1F("ESD_DalitzPositron_Pt","ESD_DalitzPositron_Pt",1000,0,25);
		fESDList[iCut]->Add(fHistoDalitzPositronPt[iCut]);
		
		fHistoDalitzElectronPhi[iCut] = new TH1F("ESD_DalitzElectron_Phi","ESD_DalitzElectron_Phi",360,0,2*TMath::Pi());
		fESDList[iCut]->Add(fHistoDalitzElectronPhi[iCut]);

		fHistoDalitzPositronPhi[iCut] = new TH1F("ESD_DalitzPositron_Phi","ESD_DalitzPositron_Phi",360,0,2*TMath::Pi());
		fESDList[iCut]->Add(fHistoDalitzPositronPhi[iCut]);
		
		
		
		
		
		

		fClusterOutputList[iCut] = new TList();
		fClusterOutputList[iCut]->SetName(Form("%s_%s_%s_%s_%s Cluster Output",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringElectron.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
		fClusterOutputList[iCut]->SetOwner(1);
		fCutFolder[iCut]->Add(fClusterOutputList[iCut]);
		
		fHistoClusGammaPt[iCut] = new TH1F("ClusGamma_Pt","ClusGamma_Pt",250,0,25);
		fHistoClusGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
		fClusterOutputList[iCut]->Add(fHistoClusGammaPt[iCut]);
		fHistoClusOverlapHeadersGammaPt[iCut] = new TH1F("ClusGammaOverlapHeaders_Pt","ClusGammaOverlapHeaders_Pt",250,0,25);
		fHistoClusOverlapHeadersGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c), selected header w/ overlap");
		fClusterOutputList[iCut]->Add(fHistoClusOverlapHeadersGammaPt[iCut]);

		if(fDoMesonAnalysis){
		  
			fHistoMotherInvMassPt[iCut] = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",800,0,0.8,250,0,25);
			fHistoMotherInvMassPt[iCut]->SetXTitle("M_{inv} (GeV/c^{2})");
			fHistoMotherInvMassPt[iCut]->SetYTitle("p_{T,pair} (GeV/c)");
			fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);
			
			fHistoMotherMatchedInvMassPt[iCut] = new TH2F("ESD_MotherMatched_InvMass_Pt","ESD_MotherMatched_InvMass_Pt",800,0,0.8,250,0,25);
			fHistoMotherMatchedInvMassPt[iCut]->SetXTitle("M_{inv} (GeV/c^{2}) matched conv e^{+/-} to cluster");
			fHistoMotherMatchedInvMassPt[iCut]->SetYTitle("p_{T,pair} (GeV/c)");
			fESDList[iCut]->Add(fHistoMotherMatchedInvMassPt[iCut]);
			
			fHistoMotherBackInvMassPt[iCut] = new TH2F("ESD_Background_InvMass_Pt","ESD_Background_InvMass_Pt",800,0,0.8,250,0,25);
			fHistoMotherBackInvMassPt[iCut]->SetXTitle("M_{inv, mxed} (GeV/c^{2})");
			fHistoMotherBackInvMassPt[iCut]->SetYTitle("p_{T,BG pair} (GeV/c)");
			fESDList[iCut]->Add(fHistoMotherBackInvMassPt[iCut]);
			
			fHistoMotherInvMassEalpha[iCut] = new TH2F("ESD_Mother_InvMass_vs_E_alpha","ESD_Mother_InvMass_vs_E_alpha",800,0,0.8,250,0,25);
			fHistoMotherInvMassEalpha[iCut]->SetXTitle("M_{inv} (GeV/c^{2})");
			fHistoMotherInvMassEalpha[iCut]->SetYTitle("p_{T,pair} (GeV/c)");
			fESDList[iCut]->Add(fHistoMotherInvMassEalpha[iCut]);
			
			if( fDoMesonQA > 0 ) {
			    fHistoMotherInvMassOpeningAngleGammaElectron[iCut] = new TH1F("ESD_MotherInvMass_OpeningAngle_GammaElectron", "ESD_MotherInvMass_OpeningAngle_GammaElectron",100,0.,TMath::Pi());
			    fESDList[iCut]->Add(fHistoMotherInvMassOpeningAngleGammaElectron[iCut]);
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
		fHeaderNameList 				= new TList*[fnCuts];
		fHistoMCHeaders 				= new TH1I*[fnCuts];
		fHistoMCAllGammaPt 				= new TH1F*[fnCuts];
		fHistoMCAllGammaPi0Pt                           = new TH1F*[fnCuts];
		fHistoMCAllGammaEMCALAccPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaPi0Pt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaRhoPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaEtaPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaOmegaPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaEtapPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaPhiPt 			= new TH1F*[fnCuts];
		fHistoMCDecayGammaSigmaPt 			= new TH1F*[fnCuts];
		fHistoMCConvGammaPt 				= new TH1F*[fnCuts];
		fHistoMCAllPositronsPt                          = new TH1F*[fnCuts];
		fHistoMCDecayPositronPi0Pt                      = new TH1F*[fnCuts];
		fHistoMCAllElectronsPt                          = new TH1F*[fnCuts];
		fHistoMCDecayElectronPi0Pt                      = new TH1F*[fnCuts];
		fHistoMCDecayNoPrimElectronPi0DalitzR		= new TH1F*[fnCuts];
		fHistoMCDecayNoPrimPositronPi0DalitzR		= new TH1F*[fnCuts];
		fHistoMCDecayNoPrimElectronPi0DalitzID		= new TH1F*[fnCuts];
		fHistoMCDecayNoPrimPositronPi0DalitzID		= new TH1F*[fnCuts];
	
		
		fHistoTrueConvGammaPt 				= new TH1F*[fnCuts];
		fHistoDoubleCountTrueConvGammaRPt		= new TH2F*[fnCuts];
		fHistoTrueConvPi0GammaPt 			= new TH1F*[fnCuts];
		fHistoTruePositronPt				= new TH1F*[fnCuts];
		fHistoTrueElectronPt				= new TH1F*[fnCuts];
		fHistoTrueSecPositronPt				= new TH1F*[fnCuts];
		fHistoTrueSecElectronPt				= new TH1F*[fnCuts];		
		fHistoTruePi0DalitzElectronPt			= new TH1F*[fnCuts];
		fHistoTruePi0DalitzPositronPt			= new TH1F*[fnCuts];
		fHistoTruePi0DalitzSecElectronPt		= new TH1F*[fnCuts];
		fHistoTruePi0DalitzSecPositronPt		= new TH1F*[fnCuts];
		
		
    
		
		fHistoTruePrimaryConvGammaPt 			= new TH1F*[fnCuts];
		fHistoTruePrimaryConvGammaESDPtMCPt 		= new TH2F*[fnCuts];
		fHistoTrueSecondaryConvGammaPt 			= new TH1F*[fnCuts];
		fHistoTrueSecondaryConvGammaFromXFromK0sPt 	= new TH1F*[fnCuts];
		fHistoTrueSecondaryConvGammaFromXFromLambdaPt 	= new TH1F*[fnCuts];
    
		fHistoTrueClusGammaPt 				= new TH1F*[fnCuts];
		fHistoTruePrimaryClusGammaPt 			= new TH1F*[fnCuts];
		fHistoTruePrimaryClusGammaESDPtMCPt 		= new TH2F*[fnCuts];
		
		fHistoTruePi0DalitzClusGammaPt			= new TH1F*[fnCuts];
		fHistoTruePi0DalitzClusGammaMCPt		= new TH1F*[fnCuts];
		fHistoTruePi0DalitzAllClusGammaPt               = new TH1F*[fnCuts];

		if (fDoPhotonQA > 0){
		  
			fHistoMCConvGammaR 			= new TH1F*[fnCuts];
			fHistoMCConvGammaEta 			= new TH1F*[fnCuts];
			fHistoTrueConvGammaEta 			= new TH1F*[fnCuts];
			
		}
		
		
		if (fDoClusterQA > 0){
		  
			fHistoTrueClusUnConvGammaPt 		= new TH1F*[fnCuts];
			fHistoTrueClusUnConvGammaMCPt 		= new TH1F*[fnCuts];
			fHistoTrueClusElectronPt 		= new TH1F*[fnCuts];
			fHistoTrueClusConvGammaPt 		= new TH1F*[fnCuts];
			fHistoTrueClusConvGammaMCPt 		= new TH1F*[fnCuts];
			fHistoTrueClusConvGammaFullyPt 		= new TH1F*[fnCuts];
			fHistoTrueClusMergedGammaPt 		= new TH1F*[fnCuts];
			fHistoTrueClusMergedPartConvGammaPt 	= new TH1F*[fnCuts];
			fHistoTrueClusDalitzPt			= new TH1F*[fnCuts];
			fHistoTrueClusDalitzMergedPt	 	= new TH1F*[fnCuts];
			fHistoTrueClusPhotonFromElecMotherPt	= new TH1F*[fnCuts];
			fHistoTrueClusShowerPt			= new TH1F*[fnCuts];
			fHistoTrueClusSubLeadingPt		= new TH1F*[fnCuts];
			fHistoTrueClusNParticles		= new TH1I*[fnCuts];
			fHistoTrueClusEMNonLeadingPt		= new TH1F*[fnCuts];
			fHistoTrueNLabelsInClus 		= new TH1F*[fnCuts];			
		}
		
		
		if(fDoMesonAnalysis){
		  
			fHistoMCPi0GGPt 				= new TH1F*[fnCuts];
			fHistoMCPi0GGWOWeightPt 			= new TH1F*[fnCuts];
			fHistoMCPi0Pt					= new TH1F*[fnCuts];
			fHistoMCPi0WOWeightPt				= new TH1F*[fnCuts];
			
			fHistoMCEtaGGPt 				= new TH1F*[fnCuts];
			fHistoMCEtaGGWOWeightPt 			= new TH1F*[fnCuts];
			fHistoMCEtaPt 					= new TH1F*[fnCuts];
			fHistoMCEtaWOWeightPt 				= new TH1F*[fnCuts];
			
			fHistoMCPi0InAccPt 				= new TH1F*[fnCuts];
			
			fHistoMCEtaInAccPt 				= new TH1F*[fnCuts];
      
			fHistoTruePi0InvMassPt 				= new TH2F*[fnCuts];
			fHistoTrueEtaInvMassPt 				= new TH2F*[fnCuts];
			fHistoDoubleCountTruePi0InvMassPt			= new TH2F*[fnCuts];
			fHistoDoubleCountTrueEtaInvMassPt			= new TH2F*[fnCuts];
			fHistoTruePi0GGInvMassPt			= new TH2F*[fnCuts];
			fHistoTrueEtaGGInvMassPt			= new TH2F*[fnCuts];
			fHistoTruePrimaryPi0InvMassPt 			= new TH2F*[fnCuts];
			fHistoTruePrimaryPi0GGInvMassPt			= new TH2F*[fnCuts];
			fHistoTruePrimaryEtaInvMassPt 			= new TH2F*[fnCuts];
			fHistoTruePrimaryEtaGGInvMassPt			= new TH2F*[fnCuts];
			fHistoTruePrimaryPi0W0WeightingInvMassPt 	= new TH2F*[fnCuts];
			fHistoTruePrimaryEtaW0WeightingInvMassPt	= new TH2F*[fnCuts];
			fProfileTruePrimaryPi0WeightsInvMassPt 		= new TProfile2D*[fnCuts];
			fProfileTruePrimaryEtaWeightsInvMassPt 		= new TProfile2D*[fnCuts];
			fHistoTrueSecondaryPi0InvMassPt 		= new TH2F*[fnCuts];
			fHistoTrueSecondaryPi0GGInvMassPt               = new TH2F*[fnCuts];
			fHistoTrueSecondaryPi0FromK0sInvMassPt 		= new TH2F*[fnCuts];
			fHistoTrueSecondaryPi0FromEtaInvMassPt 		= new TH2F*[fnCuts];
			fHistoTrueSecondaryPi0FromLambdaInvMassPt 	= new TH2F*[fnCuts];
			fHistoTruePrimaryPi0PhotonPairPtconv 		= new TH2F*[fnCuts];
			fHistoTrueSecondaryPi0PhotonPairPtconv 		= new TH2F*[fnCuts];
			fHistoTruePrimaryEtaPhotonPairPtconv 		= new TH2F*[fnCuts];
			fHistoTruePrimaryPi0DCPtconv 			= new TH1F*[fnCuts];
			fHistoTrueSecondaryPi0DCPtconv 			= new TH1F*[fnCuts];
			fHistoTruePrimaryEtaDCPtconv 			= new TH1F*[fnCuts];
			fHistoTruePrimaryPi0MissingPtconv 		= new TH1F*[fnCuts];
			fHistoTrueSecondaryPi0MissingPtconv 		= new TH1F*[fnCuts];
			fHistoTruePrimaryEtaMissingPtconv 		= new TH1F*[fnCuts];
			fStringRecTruePi0s				= new TString[fnCuts];
			fStringRecTrueEtas				= new TString[fnCuts];
			
			
			if (fDoMesonQA > 0){
			  
				fHistoMCPi0PtY 						= new TH2F*[fnCuts];
				fHistoMCEtaPtY 						= new TH2F*[fnCuts];
				fHistoMCPi0PtAlpha 					= new TH2F*[fnCuts];
				fHistoMCEtaPtAlpha 					= new TH2F*[fnCuts];
				fHistoMCK0sPt 						= new TH1F*[fnCuts];
				fHistoMCK0sWOWeightPt 					= new TH1F*[fnCuts];
				fHistoMCK0sPtY	 					= new TH2F*[fnCuts];
				fHistoTruePi0CaloPhotonInvMassPt			= new TH2F*[fnCuts];
				fHistoTrueEtaCaloPhotonInvMassPt			= new TH2F*[fnCuts];
				fHistoTruePi0CaloConvertedPhotonInvMassPt 		= new TH2F*[fnCuts];
				fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt 	= new TH2F*[fnCuts];
				fHistoTrueEtaCaloConvertedPhotonInvMassPt		= new TH2F*[fnCuts];
				fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt 	= new TH2F*[fnCuts];
				fHistoTruePi0CaloElectronInvMassPt			= new TH2F*[fnCuts];
				fHistoTrueEtaCaloElectronInvMassPt			= new TH2F*[fnCuts];
				fHistoTruePi0CaloMergedClusterInvMassPt 		= new TH2F*[fnCuts];
				fHistoTrueEtaCaloMergedClusterInvMassPt 		= new TH2F*[fnCuts];
				fHistoTruePi0CaloMergedClusterPartConvInvMassPt 	= new TH2F*[fnCuts];
				fHistoTrueEtaCaloMergedClusterPartConvInvMassPt 	= new TH2F*[fnCuts];
				fHistoTrueMotherCaloEMNonLeadingInvMassPt 		= new TH2F*[fnCuts];
				fHistoTruePrimaryPi0MCPtResolPt 			= new TH2F*[fnCuts];
				fHistoTruePrimaryEtaMCPtResolPt 			= new TH2F*[fnCuts];
				fHistoTrueK0sWithPi0DaughterMCPt 			= new TH1F*[fnCuts];
				fHistoTrueEtaWithPi0DaughterMCPt 			= new TH1F*[fnCuts];
				fHistoTrueLambdaWithPi0DaughterMCPt 			= new TH1F*[fnCuts];
				fHistoTrueBckGGInvMassPt 				= new TH2F*[fnCuts];
				fHistoTrueBckContInvMassPt 				= new TH2F*[fnCuts];
				fHistoTruePi0PtY 					= new TH2F*[fnCuts];
				fHistoTrueEtaPtY 					= new TH2F*[fnCuts];
				fHistoTruePi0PtAlpha 					= new TH2F*[fnCuts];
				fHistoTrueEtaPtAlpha 					= new TH2F*[fnCuts];
				fHistoTruePi0PtOpenAngle 				= new TH2F*[fnCuts];
				fHistoTrueEtaPtOpenAngle 				= new TH2F*[fnCuts];
				fHistoTrueMotherPi0ConvPhotonEtaPhi         		= new TH2F*[fnCuts];
				fHistoTrueMotherEtaConvPhotonEtaPhi			= new TH2F*[fnCuts];
				fHistoMCPi0InAccOpeningAngleGammaElectron       	= new TH1F*[fnCuts];
				fHistoTruePi0ShowerInvMassPt				= new TH2F*[fnCuts];
				fHistoTrueEtaShowerInvMassPt				= new TH2F*[fnCuts];
				fHistoTruePi0NoShowerInvMassPt				= new TH2F*[fnCuts];
				fHistoTrueEtaNoShowerInvMassPt				= new TH2F*[fnCuts];
				fHistoTruePi0OpeningAngleGammaElectron                  = new TH1F*[fnCuts];
				
			}
		}
     
		for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		  
			TString cutstringEvent 	  = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
			TString cutstringPhoton   = ((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutNumber();
			TString cutstringElectron = ((AliDalitzElectronCuts*)fElectronCutArray->At(iCut))->GetCutNumber();
			TString cutstringCalo 	  = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
			TString cutstringMeson 	  = "NoMesonCut";
			if(fDoMesonAnalysis)cutstringMeson = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
      
			fMCList[iCut] = new TList();
			fMCList[iCut]->SetName(Form("%s_%s_%s_%s_%s MC histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringElectron.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
			fMCList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fMCList[iCut]);
			fHistoMCHeaders[iCut] = new TH1I("MC_Headers","MC_Headers",20,0,20);
			fMCList[iCut]->Add(fHistoMCHeaders[iCut]);
			fHistoMCAllGammaPt[iCut] = new TH1F("MC_AllGamma_Pt","MC_AllGamma_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCAllGammaPt[iCut]);
			
			fHistoMCAllGammaPi0Pt[iCut] = new TH1F("MC_AllGammaPi0_Pt","MC_AllGammaPi0_Pt",250,0,25);
			fMCList[iCut]->Add(fHistoMCAllGammaPi0Pt[iCut]);
			
			fHistoMCAllGammaEMCALAccPt[iCut] = new TH1F("MC_AllGammaEMCALAcc_Pt","MC_AllGammaEMCALAcc_Pt",250,0,25);
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
			
			fHistoMCAllPositronsPt[iCut] = new TH1F("MC_AllPositrons_Pt","MC_AllPositrons_Pt",1000,0,25);
			fMCList[iCut]->Add(fHistoMCAllPositronsPt[iCut]);
			
			
			fHistoMCDecayPositronPi0Pt[iCut] = new TH1F("MC_DecayPositronPi0_Pt","MC_DecayPositronPi0_Pt",1000,0,25);
			fMCList[iCut]->Add(fHistoMCDecayPositronPi0Pt[iCut]);
			
			fHistoMCAllElectronsPt[iCut] = new TH1F("MC_AllElectrons_Pt","MC_AllElectrons_Pt",1000,0,25);
			fMCList[iCut]->Add(fHistoMCAllElectronsPt[iCut]);
			
			fHistoMCDecayElectronPi0Pt[iCut] = new TH1F("MC_DecayElectronPi0_Pt","MC_DecayElectronPi0_Pt",1000,0,25);
			fMCList[iCut]->Add(fHistoMCDecayElectronPi0Pt[iCut]);
			
			fHistoMCDecayNoPrimElectronPi0DalitzR[iCut] = new TH1F("MC_DecayNoPrimElectronPi0Dalitz_R","MC_DecayNoPrimElectronPi0Dalitz_R",200,0,200);
			fMCList[iCut]->Add(fHistoMCDecayNoPrimElectronPi0DalitzR[iCut]);
			
			fHistoMCDecayNoPrimPositronPi0DalitzR[iCut]  = new TH1F("MC_DecayNoPrimPositronPi0Dalitz_R","MC_DecayNoPrimPositronPi0Dalitz_R",200,0,200);
			fMCList[iCut]->Add(fHistoMCDecayNoPrimPositronPi0DalitzR[iCut]);
			
			fHistoMCDecayNoPrimElectronPi0DalitzID[iCut] = new TH1F("MC_DecayNoPrimElectronPi0Dalitz_ID","MC_DecayNoPrimElectronPi0Dalitz_ID",51,-1,50);
			fMCList[iCut]->Add(fHistoMCDecayNoPrimElectronPi0DalitzID[iCut]);
			
			fHistoMCDecayNoPrimPositronPi0DalitzID[iCut] = new TH1F("MC_DecayNoPrimPositronPi0Dalitz_ID","MC_DecayNoPrimPositronPi0Dalitz_ID",51,-1,50);
			fMCList[iCut]->Add(fHistoMCDecayNoPrimPositronPi0DalitzID[iCut]);
					
			
      
			if (fDoPhotonQA > 0){
			  
				fHistoMCConvGammaR[iCut] = new TH1F("MC_ConvGamma_R","MC_ConvGamma_R",800,0,200);
				fMCList[iCut]->Add(fHistoMCConvGammaR[iCut]);
				fHistoMCConvGammaEta[iCut] = new TH1F("MC_ConvGamma_Eta","MC_ConvGamma_Eta",2000,-2,2);
				fMCList[iCut]->Add(fHistoMCConvGammaEta[iCut]);
				
			}
      
			if(fDoMesonAnalysis){
			  
				fHistoMCPi0GGPt[iCut] = new TH1F("MC_Pi0_GG_Pt","MC_Pi0_GG_Pt",250,0,25);
				fHistoMCPi0GGPt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCPi0GGPt[iCut]);
				fHistoMCPi0GGWOWeightPt[iCut] = new TH1F("MC_Pi0_GG_WOWeights_Pt","MC_Pi0_GG_WOWeights_Pt",250,0,25);
				fHistoMCPi0GGWOWeightPt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCPi0GGWOWeightPt[iCut]);
				
				fHistoMCPi0Pt[iCut] = new TH1F("MC_Pi0_Pt","MC_Pi0_Pt",250,0,25);
				fHistoMCPi0Pt[iCut]->Sumw2();
				
				fMCList[iCut]->Add(fHistoMCPi0Pt[iCut]);
				fHistoMCPi0WOWeightPt[iCut] = new TH1F("MC_Pi0_WOWeights_Pt","MC_Pi0_WOWeights_Pt",250,0,25);
				fHistoMCPi0WOWeightPt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCPi0WOWeightPt[iCut]);
				
				
				fHistoMCEtaGGPt[iCut] = new TH1F("MC_Eta_GG_Pt","MC_Eta_GG_Pt",250,0,25);
				fHistoMCEtaGGPt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCEtaGGPt[iCut]);
				fHistoMCEtaGGWOWeightPt[iCut] = new TH1F("MC_Eta_GG_WOWeights_Pt","MC_Eta_GG_WOWeights_Pt",250,0,25);
				fHistoMCEtaGGWOWeightPt[iCut]->Sumw2();
				fMCList[iCut]->Add(fHistoMCEtaGGWOWeightPt[iCut]);
				
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
					fHistoMCPi0PtAlpha[iCut] = new TH2F("MC_Pi0_Pt_Alpha","MC_Pi0_Pt_Alpha",150,0.03,15.,100,0,1);
					SetLogBinningXTH2(fHistoMCPi0PtAlpha[iCut]);
					fMCList[iCut]->Add(fHistoMCPi0PtAlpha[iCut]);
					fHistoMCEtaPtAlpha[iCut] = new TH2F("MC_Eta_Pt_Alpha","MC_Eta_Pt_Alpha",150,0.03,15.,100,0,1);
					SetLogBinningXTH2(fHistoMCEtaPtAlpha[iCut]);
					fMCList[iCut]->Add(fHistoMCEtaPtAlpha[iCut]);

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
					
					fHistoMCPi0InAccOpeningAngleGammaElectron[iCut] = new TH1F("MC_Pi0InAcc_OpeningAngle_GammaElectron","MC_Pi0InAcc_OpeningAngle_GammaElectron",100,0,TMath::Pi());
					fMCList[iCut]->Add(fHistoMCPi0InAccOpeningAngleGammaElectron[iCut]);
					
					
				}
        
			}
			fTrueList[iCut] = new TList();
			fTrueList[iCut]->SetName(Form("%s_%s_%s_%s_%s True histograms",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringElectron.Data(),cutstringCalo.Data(),cutstringMeson.Data()));
			fTrueList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fTrueList[iCut]);
      
			fHistoTrueConvGammaPt[iCut] = new TH1F("ESD_TrueConvGamma_Pt","ESD_TrueConvGamma_Pt",250,0,25);
			fTrueList[iCut]->Add(fHistoTrueConvGammaPt[iCut]);

			fHistoDoubleCountTrueConvGammaRPt[iCut] = new TH2F("ESD_TrueDoubleCountConvGamma_R_Pt","ESD_TrueDoubleCountConvGamma_R_Pt",800,0,200,300,0,30);
			fTrueList[iCut]->Add(fHistoDoubleCountTrueConvGammaRPt[iCut]);

			fHistoTrueConvPi0GammaPt[iCut] = new TH1F("ESD_TrueConvPi0Gamma_Pt","ESD_TrueConvPi0Gamma_Pt",250,0,25);
			fTrueList[iCut]->Add(fHistoTrueConvPi0GammaPt[iCut]);
			
			
			fHistoTruePositronPt[iCut] = new TH1F("ESD_TruePositron_Pt","ESD_TruePositron_Pt",1000,0,25);
			fTrueList[iCut]->Add(fHistoTruePositronPt[iCut]);

			fHistoTrueElectronPt[iCut] = new TH1F("ESD_TrueElectron_Pt","ESD_TrueElectron_Pt",1000,0,25);
			fTrueList[iCut]->Add(fHistoTrueElectronPt[iCut]);
		
			fHistoTrueSecPositronPt[iCut] = new TH1F("ESD_TrueSecPositron_Pt","ESD_TrueSecPositron_Pt",1000,0,25);
			fTrueList[iCut]->Add(fHistoTrueSecPositronPt[iCut]);

			fHistoTrueSecElectronPt[iCut] = new TH1F("ESD_TrueSecElectron_Pt","ESD_TrueSecElectron_Pt",1000,0,25);
			fTrueList[iCut]->Add(fHistoTrueSecElectronPt[iCut]); 
													
			fHistoTruePi0DalitzElectronPt[iCut] = new TH1F("ESD_TruePi0DalitzElectron_Pt","ESD_TruePi0DalitzElectron_Pt",1000,0,25);
			fTrueList[iCut]->Add(fHistoTruePi0DalitzElectronPt[iCut]);

			fHistoTruePi0DalitzPositronPt[iCut] = new TH1F("ESD_TruePi0DalitzPositron_Pt","ESD_TruePi0DalitzPositron_Pt",1000,0,25);
			fTrueList[iCut]->Add(fHistoTruePi0DalitzPositronPt[iCut]);
		
			fHistoTruePi0DalitzSecElectronPt[iCut] = new TH1F("ESD_TruePi0DalitzSecElectron_Pt","ESD_TruePi0DalitzSecElectron_Pt",1000,0,25);
			fTrueList[iCut]->Add(fHistoTruePi0DalitzSecElectronPt[iCut]);

			fHistoTruePi0DalitzSecPositronPt[iCut] = new TH1F("ESD_TruePi0DalitzSecPositron_Pt","ESD_TruePi0DalitzSecPositron_Pt",1000,0,25);
			fTrueList[iCut]->Add(fHistoTruePi0DalitzSecPositronPt[iCut]);

				
			
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
			fClusterOutputList[iCut]->Add(fHistoTrueClusGammaPt[iCut]);
			fHistoTruePrimaryClusGammaPt[iCut] = new TH1F("TruePrimaryClusGamma_Pt","ESD_TruePrimaryClusGamma_Pt",250,0,25);
			fClusterOutputList[iCut]->Add(fHistoTruePrimaryClusGammaPt[iCut]);
			fHistoTruePrimaryClusGammaESDPtMCPt[iCut] = new TH2F("TruePrimaryClusGamma_Pt_MCPt","ESD_TruePrimaryClusGamma_MCPt",250,0,25,250,0,25);
			fClusterOutputList[iCut]->Add(fHistoTruePrimaryClusGammaESDPtMCPt[iCut]);
			
			 fHistoTruePi0DalitzClusGammaPt[iCut] = new TH1F("TruePi0DalitzClusGamma_Pt","fHistoTruePi0DalitzClusGamma_Pt",250,0,25);
			 fClusterOutputList[iCut]->Add(fHistoTruePi0DalitzClusGammaPt[iCut]);
			 
	       	         fHistoTruePi0DalitzClusGammaMCPt[iCut]= new TH1F("TruePi0DalitzClusGamma_MCPt","fHistoTruePi0DalitzClusGamma_MCPt",250,0,25);
			 fClusterOutputList[iCut]->Add(fHistoTruePi0DalitzClusGammaMCPt[iCut]);
			 
			
			 fHistoTruePi0DalitzAllClusGammaPt[iCut] = new TH1F("TruePi0DalitzAllClusGamma_Pt","TruePi0DalitzAllClusGamma_Pt",250,0,25);
			 fClusterOutputList[iCut]->Add(fHistoTruePi0DalitzAllClusGammaPt[iCut]);
			
			
			

			if (fDoPhotonQA > 0){
				fHistoTrueConvGammaEta[iCut] = new TH1F("ESD_TrueConvGamma_Eta","ESD_TrueConvGamma_Eta",2000,-2,2);
				fTrueList[iCut]->Add(fHistoTrueConvGammaEta[iCut]);		
			}	
			if (fDoClusterQA > 0){	
				fHistoTrueClusUnConvGammaPt[iCut] = new TH1F("TrueClusUnConvGamma_Pt","TrueClusUnConvGamma_Pt",250,0,25);
				fClusterOutputList[iCut]->Add(fHistoTrueClusUnConvGammaPt[iCut]);
				fHistoTrueClusUnConvGammaMCPt[iCut] = new TH1F("TrueClusUnConvGamma_MCPt","TrueClusUnConvGamma_MCPt",250,0,25);
				fClusterOutputList[iCut]->Add(fHistoTrueClusUnConvGammaMCPt[iCut]);
				fHistoTrueClusElectronPt[iCut] = new TH1F("TrueClusElectron_Pt","TrueElectronGamma_Pt",250,0,25);
				fClusterOutputList[iCut]->Add(fHistoTrueClusElectronPt[iCut]);
				fHistoTrueClusConvGammaPt[iCut] = new TH1F("TrueClusConvGamma_Pt","TrueClusConvGamma_Pt",250,0,25);
				fClusterOutputList[iCut]->Add(fHistoTrueClusConvGammaPt[iCut]);
				fHistoTrueClusConvGammaMCPt[iCut] = new TH1F("TrueClusConvGamma_MCPt","TrueClusConvGamma_MCPt",250,0,25);
				fClusterOutputList[iCut]->Add(fHistoTrueClusConvGammaMCPt[iCut]);
				fHistoTrueClusConvGammaFullyPt[iCut] = new TH1F("TrueClusConvGammaFullyContained_Pt","TrueClusConvGammaFullyContained_Pt",250,0,25);
				fClusterOutputList[iCut]->Add(fHistoTrueClusConvGammaFullyPt[iCut]);
				fHistoTrueClusMergedGammaPt[iCut] = new TH1F("TrueClusMergedGamma_Pt","TrueClusMergedGamma_Pt",250,0,25);
				fClusterOutputList[iCut]->Add(fHistoTrueClusMergedGammaPt[iCut]);
				fHistoTrueClusMergedPartConvGammaPt[iCut] = new TH1F("TrueClusMergedPartConvGamma_Pt","TrueClusMergedPartConvGamma_Pt",250,0,25);
				fClusterOutputList[iCut]->Add(fHistoTrueClusMergedPartConvGammaPt[iCut]);
				fHistoTrueClusDalitzPt[iCut] = new TH1F("TrueClusDalitz_Pt","TrueClusDalitz_Pt",250,0,25);
				fClusterOutputList[iCut]->Add(fHistoTrueClusDalitzPt[iCut]);
				fHistoTrueClusDalitzMergedPt[iCut] = new TH1F("TrueClusDalitzMerged_Pt","TrueClusDalitzMerged_Pt",250,0,25);
				fClusterOutputList[iCut]->Add(fHistoTrueClusDalitzMergedPt[iCut]);
				fHistoTrueClusPhotonFromElecMotherPt[iCut] = new TH1F("TrueClusPhotonFromElecMother_Pt","TrueClusPhotonFromElecMother_Pt",250,0,25);
				fClusterOutputList[iCut]->Add(fHistoTrueClusPhotonFromElecMotherPt[iCut]);
				fHistoTrueClusShowerPt[iCut] = new TH1F("TrueClusShower_Pt","TrueClusShower_Pt",250,0,25);
				fClusterOutputList[iCut]->Add(fHistoTrueClusShowerPt[iCut]);
				fHistoTrueClusSubLeadingPt[iCut] = new TH1F("TrueClusSubleading_Pt","TrueClusSubleading_Pt",250,0,25);
				fClusterOutputList[iCut]->Add(fHistoTrueClusSubLeadingPt[iCut]);
				fHistoTrueClusNParticles[iCut] = new TH1I("TrueClusNParticles","TrueClusNParticles",20,0,20);
				fClusterOutputList[iCut]->Add(fHistoTrueClusNParticles[iCut]);
				fHistoTrueClusEMNonLeadingPt[iCut] = new TH1F("TrueClusEMNonLeading_Pt","TrueClusEMNonLeading_Pt",250,0,25);
				fClusterOutputList[iCut]->Add(fHistoTrueClusEMNonLeadingPt[iCut]);
				fHistoTrueNLabelsInClus[iCut] = new TH1F("TrueNLabelsInClus","TrueNLabelsInClus",100,-0.5,99.5);
				fClusterOutputList[iCut]->Add(fHistoTrueNLabelsInClus[iCut]);	
			}	

			if(fDoMesonAnalysis){
				fHistoTruePi0InvMassPt[iCut] = new TH2F("ESD_TruePi0_InvMass_Pt","ESD_TruePi0_InvMass_Pt",800,0,0.8,250,0,25);
				fHistoTruePi0InvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}} (GeV/c^{2})");
				fHistoTruePi0InvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
				fTrueList[iCut]->Add(fHistoTruePi0InvMassPt[iCut]);
				
				fHistoTrueEtaInvMassPt[iCut] = new TH2F("ESD_TrueEta_InvMass_Pt","ESD_TrueEta_InvMass_Pt",800,0,0.8,250,0,25);
				fHistoTrueEtaInvMassPt[iCut]->SetXTitle("M_{inv,#eta} (GeV/c^{2})");
				fHistoTrueEtaInvMassPt[iCut]->SetYTitle("#eta p_{T} (GeV/c)");
				fTrueList[iCut]->Add(fHistoTrueEtaInvMassPt[iCut]);
				
				fHistoDoubleCountTruePi0InvMassPt[iCut] = new TH2F("ESD_TrueDoubleCountPi0_InvMass_Pt","ESD_TrueDoubleCountPi0_InvMass_Pt",800,0,0.8,300,0,30);
				fTrueList[iCut]->Add(fHistoDoubleCountTruePi0InvMassPt[iCut]);
				fHistoDoubleCountTrueEtaInvMassPt[iCut] = new TH2F("ESD_TrueDoubleCountEta_InvMass_Pt","ESD_TrueDoubleCountEta_InvMass_Pt",800,0,0.8,300,0,30);
				fTrueList[iCut]->Add(fHistoDoubleCountTrueEtaInvMassPt[iCut]);

				fHistoTruePi0GGInvMassPt[iCut] = new TH2F("ESD_TruePi0GG_InvMass_Pt","ESD_TruePi0GG_InvMass_Pt",800,0,0.8,250,0,25);
				fHistoTruePi0GGInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}} (GeV/c^{2})");
				fHistoTruePi0GGInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
				fTrueList[iCut]->Add(fHistoTruePi0GGInvMassPt[iCut]);
				
				fHistoTrueEtaGGInvMassPt[iCut] = new TH2F("ESD_TrueEtaGG_InvMass_Pt","ESD_TrueEtaGG_InvMass_Pt",800,0,0.8,250,0,25);
				fHistoTrueEtaGGInvMassPt[iCut]->SetXTitle("M_{inv,#eta} (GeV/c^{2})");
				fHistoTrueEtaGGInvMassPt[iCut]->SetYTitle("#eta p_{T} (GeV/c)");
				fTrueList[iCut]->Add(fHistoTrueEtaGGInvMassPt[iCut]);
				
				

				fHistoTruePrimaryPi0InvMassPt[iCut] = new TH2F("ESD_TruePrimaryPi0_InvMass_Pt", "ESD_TruePrimaryPi0_InvMass_Pt", 800,0,0.8,250,0,25);
				fHistoTruePrimaryPi0InvMassPt[iCut]->SetXTitle("M_{inv,prim #pi^{0}} (GeV/c^{2})");
				fHistoTruePrimaryPi0InvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
				fHistoTruePrimaryPi0InvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTruePrimaryPi0InvMassPt[iCut]);
				
				
				fHistoTruePrimaryPi0GGInvMassPt[iCut] = new TH2F("ESD_TruePrimaryPi0GG_InvMass_Pt", "ESD_TruePrimaryPi0GG_InvMass_Pt", 800,0,0.8,250,0,25);
				fHistoTruePrimaryPi0GGInvMassPt[iCut]->SetXTitle("M_{inv,prim #pi^{0}} (GeV/c^{2})");
				fHistoTruePrimaryPi0GGInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
				fHistoTruePrimaryPi0GGInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTruePrimaryPi0GGInvMassPt[iCut]);
				
				
				fHistoTruePrimaryEtaInvMassPt[iCut] = new TH2F("ESD_TruePrimaryEta_InvMass_Pt", "ESD_TruePrimaryEta_InvMass_Pt", 800,0,0.8,250,0,25);
				fHistoTruePrimaryEtaInvMassPt[iCut]->SetXTitle("M_{inv,prim #eta} (GeV/c^{2})");
				fHistoTruePrimaryEtaInvMassPt[iCut]->SetYTitle("#eta p_{T} (GeV/c)");
				fHistoTruePrimaryEtaInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTruePrimaryEtaInvMassPt[iCut]);
				
				fHistoTruePrimaryEtaGGInvMassPt[iCut] = new TH2F("ESD_TruePrimaryEtaGG_InvMass_Pt", "ESD_TruePrimaryEtaGG_InvMass_Pt", 800,0,0.8,250,0,25);
				fHistoTruePrimaryEtaGGInvMassPt[iCut]->SetXTitle("M_{inv,prim #eta} (GeV/c^{2})");
				fHistoTruePrimaryEtaGGInvMassPt[iCut]->SetYTitle("#eta p_{T} (GeV/c)");
				fHistoTruePrimaryEtaGGInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTruePrimaryEtaGGInvMassPt[iCut]);
				

				fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut] = new TH2F("ESD_TruePrimaryPi0W0Weights_InvMass_Pt", "ESD_TruePrimaryPi0W0Weights_InvMass_Pt", 800,0,0.8,250,0,25);
				fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->SetXTitle("M_{inv,prim #pi^{0}} (GeV/c^{2})");
				fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
				fTrueList[iCut]->Add(fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]);
				
				fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut] = new TH2F("ESD_TruePrimaryEtaW0Weights_InvMass_Pt", "ESD_TruePrimaryEtaW0Weights_InvMass_Pt", 800,0,0.8,250,0,25);
				fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->SetXTitle("M_{inv,prim #eta} (GeV/c^{2})");
				fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->SetYTitle("#eta p_{T} (GeV/c)");
				fTrueList[iCut]->Add(fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]);

				fProfileTruePrimaryPi0WeightsInvMassPt[iCut] = new TProfile2D("ESD_TruePrimaryPi0Weights_InvMass_Pt", "ESD_TruePrimaryPi0Weights_InvMass_Pt", 800,0,0.8,250,0,25);
				fProfileTruePrimaryPi0WeightsInvMassPt[iCut]->SetXTitle("M_{inv,prim #pi^{0}} (GeV/c^{2})");
				fProfileTruePrimaryPi0WeightsInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
				fProfileTruePrimaryPi0WeightsInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fProfileTruePrimaryPi0WeightsInvMassPt[iCut]);
				fProfileTruePrimaryEtaWeightsInvMassPt[iCut] = new TProfile2D("ESD_TruePrimaryEtaWeights_InvMass_Pt", "ESD_TruePrimaryEtaWeights_InvMass_Pt", 800,0,0.8,250,0,25);
				fProfileTruePrimaryEtaWeightsInvMassPt[iCut]->SetXTitle("M_{inv,prim #eta} (GeV/c^{2})");
				fProfileTruePrimaryEtaWeightsInvMassPt[iCut]->SetYTitle("#eta p_{T} (GeV/c)");
				fProfileTruePrimaryEtaWeightsInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fProfileTruePrimaryEtaWeightsInvMassPt[iCut]);

				fHistoTrueSecondaryPi0InvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0_InvMass_Pt", "ESD_TrueSecondaryPi0_InvMass_Pt", 800,0,0.8,250,0,25);
				fHistoTrueSecondaryPi0InvMassPt[iCut]->SetXTitle("M_{inv,sec #pi^{0}} (GeV/c^{2})");
				fHistoTrueSecondaryPi0InvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
				fHistoTrueSecondaryPi0InvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTrueSecondaryPi0InvMassPt[iCut]);
				
				
				fHistoTrueSecondaryPi0GGInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0GG_InvMass_Pt", "ESD_TrueSecondaryPi0GG_InvMass_Pt", 800,0,0.8,250,0,25);
				fHistoTrueSecondaryPi0GGInvMassPt[iCut]->SetXTitle("M_{inv,sec #pi^{0}} (GeV/c^{2})");
				fHistoTrueSecondaryPi0GGInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
				fHistoTrueSecondaryPi0GGInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTrueSecondaryPi0GGInvMassPt[iCut]);
				
				
				

				fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromK0s_InvMass_Pt","ESD_TrueSecondaryPi0FromK0s_InvMass_Pt",800,0,0.8,250,0,25);
				fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0} from K^{0}_{S}} (GeV/c^{2})");
				fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
				fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]);
				fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromEta_InvMass_Pt","ESD_TrueSecondaryPi0FromEta_InvMass_Pt",800,0,0.8,250,0,25);
				fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0} from #eta} (GeV/c^{2})");
				fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]->SetYTitle("#pi^{0}  p_{T} (GeV/c)");
				fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]);
				fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromLambda_InvMass_Pt","ESD_TrueSecondaryPi0FromLambda_InvMass_Pt",800,0,0.8,250,0,25);
				fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0} from #Lambda} (GeV/c^{2})");
				fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]->SetYTitle("#pi^{0}  p_{T} (GeV/c)");
				fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]);

				fHistoTruePrimaryPi0PhotonPairPtconv[iCut] = new TH2F("ESD_TruePrimaryPi0_InvMass_PtConv","",800,0,0.8,250,0,25);
				fHistoTruePrimaryPi0PhotonPairPtconv[iCut]->SetXTitle("M_{inv,#pi^{0}} (GeV/c^{2})");
				fHistoTruePrimaryPi0PhotonPairPtconv[iCut]->SetYTitle("#gamma^{conv} p_{T} (GeV/c)");
				fHistoTruePrimaryPi0PhotonPairPtconv[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTruePrimaryPi0PhotonPairPtconv[iCut]);

				fHistoTrueSecondaryPi0PhotonPairPtconv[iCut] = new TH2F("ESD_TrueSecondaryPi0_InvMass_PtConv","",800,0,0.8,250,0,25);
				fHistoTrueSecondaryPi0PhotonPairPtconv[iCut]->SetXTitle("M_{inv,#pi^{0}} (GeV/c^{2})");
				fHistoTrueSecondaryPi0PhotonPairPtconv[iCut]->SetYTitle("#gamma^{conv} p_{T} (GeV/c)");
				fHistoTrueSecondaryPi0PhotonPairPtconv[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTrueSecondaryPi0PhotonPairPtconv[iCut]);

				fHistoTruePrimaryEtaPhotonPairPtconv[iCut] = new TH2F("ESD_TruePrimaryEta_InvMass_PtConv","",800,0,0.8,250,0,25);
				fHistoTruePrimaryEtaPhotonPairPtconv[iCut]->SetXTitle("M_{inv,#eta} (GeV/c^{2})");
				fHistoTruePrimaryEtaPhotonPairPtconv[iCut]->SetYTitle("#gamma^{conv} p_{T} (GeV/c)");
				fHistoTruePrimaryEtaPhotonPairPtconv[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTruePrimaryEtaPhotonPairPtconv[iCut]);

				fHistoTruePrimaryPi0DCPtconv[iCut] = new TH1F("ESD_TruePrimaryPi0DC_PtConv","",250,0,25);
				fHistoTruePrimaryPi0DCPtconv[iCut]->SetXTitle("#gamma^{conv} p_{T} (GeV/c)");
				fHistoTruePrimaryPi0DCPtconv[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTruePrimaryPi0DCPtconv[iCut]);

				fHistoTrueSecondaryPi0DCPtconv[iCut] = new TH1F("ESD_TrueSecondaryPi0DC_PtConv","",250,0,25);
				fHistoTrueSecondaryPi0DCPtconv[iCut]->SetXTitle("#gamma^{conv} p_{T} (GeV/c)");
				fHistoTrueSecondaryPi0DCPtconv[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTrueSecondaryPi0DCPtconv[iCut]);

				fHistoTruePrimaryEtaDCPtconv[iCut] = new TH1F("ESD_TruePrimaryEtaDC_PtConv","",250,0,25);
				fHistoTruePrimaryEtaDCPtconv[iCut]->SetXTitle("#gamma^{conv} p_{T} (GeV/c)");
				fHistoTruePrimaryEtaDCPtconv[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTruePrimaryEtaDCPtconv[iCut]);

				fHistoTruePrimaryPi0MissingPtconv[iCut] = new TH1F("ESD_TruePrimaryPi0Missing_PtConv","",250,0,25);
				fHistoTruePrimaryPi0MissingPtconv[iCut]->SetXTitle("#gamma^{conv} p_{T} (GeV/c)");
				fHistoTruePrimaryPi0MissingPtconv[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTruePrimaryPi0MissingPtconv[iCut]);

				fHistoTrueSecondaryPi0MissingPtconv[iCut] = new TH1F("ESD_TrueSecondaryPi0Missing_PtConv","",250,0,25);
				fHistoTrueSecondaryPi0MissingPtconv[iCut]->SetXTitle("#gamma^{conv} p_{T} (GeV/c)");
				fHistoTrueSecondaryPi0MissingPtconv[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTrueSecondaryPi0MissingPtconv[iCut]);

				fHistoTruePrimaryEtaMissingPtconv[iCut] = new TH1F("ESD_TruePrimaryEtaMissing_PtConv","",250,0,25);
				fHistoTruePrimaryEtaMissingPtconv[iCut]->SetXTitle("#gamma^{conv} p_{T} (GeV/c)");
				fHistoTruePrimaryEtaMissingPtconv[iCut]->Sumw2();
				fTrueList[iCut]->Add(fHistoTruePrimaryEtaMissingPtconv[iCut]);
				
				
				if (fDoMesonQA > 0){
					fHistoTruePi0CaloPhotonInvMassPt[iCut] = new TH2F("ESD_TruePi0CaloPhoton_InvMass_Pt","ESD_TruePi0CaloPhoton_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTruePi0CaloPhotonInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}} (GeV/c^{2}) #gamma #gamma");
					fHistoTruePi0CaloPhotonInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTruePi0CaloPhotonInvMassPt[iCut]);
					
					fHistoTrueEtaCaloPhotonInvMassPt[iCut] = new TH2F("ESD_TrueEtaCaloPhoton_InvMass_Pt","ESD_TrueEtaCaloPhoton_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTrueEtaCaloPhotonInvMassPt[iCut]->SetXTitle("M_{inv,#eta} (GeV/c^{2}) #gamma #gamma");
					fHistoTrueEtaCaloPhotonInvMassPt[iCut]->SetYTitle("#eta p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTrueEtaCaloPhotonInvMassPt[iCut]);
					
					fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut] = new TH2F("ESD_TruePi0CaloConvertedPhoton_InvMass_Pt","ESD_TruePi0CaloConvertedPhoton_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}} (GeV/c^{2}) #gamma #gamma_{conv}");
					fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTruePi0CaloConvertedPhotonInvMassPt[iCut]);
					
					fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt[iCut] = new TH2F("ESD_TruePi0CaloConvertedPhotonMatched_InvMass_Pt","ESD_TruePi0CaloConvertedPhotonMatched_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}} (GeV/c^{2}) #gamma #gamma_{conv,matched}");
					fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt[iCut]);
					
					fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut] = new TH2F("ESD_TrueEtaCaloConvertedPhoton_InvMass_Pt","ESD_TrueEtaCaloConvertedPhoton_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut]->SetXTitle("M_{inv,#eta} (GeV/c^{2}) #gamma #gamma_{conv}");
					fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut]->SetYTitle("#eta p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTrueEtaCaloConvertedPhotonInvMassPt[iCut]);
					
					fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt[iCut] = new TH2F("ESD_TrueEtaCaloConvertedPhotonMatched_InvMass_Pt","ESD_TrueEtaCaloConvertedPhotonMatched_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt[iCut]->SetXTitle("M_{inv,#eta} (GeV/c^{2}) #gamma #gamma_{conv,matched}");
					fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt[iCut]->SetYTitle("#eta p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt[iCut]);

					fHistoTruePi0CaloElectronInvMassPt[iCut] = new TH2F("ESD_TruePi0CaloElectron_InvMass_Pt","ESD_TruePi0CaloElectron_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTruePi0CaloElectronInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}} (GeV/c^{2}) #gamma e^{#pm}");
					fHistoTruePi0CaloElectronInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTruePi0CaloElectronInvMassPt[iCut]);
					fHistoTrueEtaCaloElectronInvMassPt[iCut] = new TH2F("ESD_TrueEtaCaloElectron_InvMass_Pt","ESD_TrueEtaCaloElectron_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTrueEtaCaloElectronInvMassPt[iCut]->SetXTitle("M_{inv,#eta} (GeV/c^{2}) #gamma e^{#pm}");
					fHistoTrueEtaCaloElectronInvMassPt[iCut]->SetYTitle("#eta p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTrueEtaCaloElectronInvMassPt[iCut]);

					fHistoTruePi0CaloMergedClusterInvMassPt[iCut] = new TH2F("ESD_TruePi0CaloMergedCluster_InvMass_Pt","ESD_TruePi0CaloMergedCluster_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTruePi0CaloMergedClusterInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}} (GeV/c^{2}) #gamma merged cluster");
					fHistoTruePi0CaloMergedClusterInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTruePi0CaloMergedClusterInvMassPt[iCut]);
					fHistoTrueEtaCaloMergedClusterInvMassPt[iCut] = new TH2F("ESD_TrueEtaCaloMergedCluster_InvMass_Pt","ESD_TrueEtaCaloMergedCluster_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTrueEtaCaloMergedClusterInvMassPt[iCut]->SetXTitle("M_{inv,#eta} (GeV/c^{2}) #gamma merged cluster");
					fHistoTrueEtaCaloMergedClusterInvMassPt[iCut]->SetYTitle("#eta p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTrueEtaCaloMergedClusterInvMassPt[iCut]);

					fHistoTrueMotherCaloEMNonLeadingInvMassPt[iCut] = new TH2F("ESD_TrueMotherCaloEMNonLeading_InvMass_Pt","ESD_TrueMotherCaloEMNonLeading_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTrueMotherCaloEMNonLeadingInvMassPt[iCut]->SetXTitle("M_{inv} (GeV/c^{2}) #gamma cluster no leading EM");
					fHistoTrueMotherCaloEMNonLeadingInvMassPt[iCut]->SetYTitle("#pair p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTrueMotherCaloEMNonLeadingInvMassPt[iCut]);
					fHistoTruePi0CaloMergedClusterPartConvInvMassPt[iCut] = new TH2F("ESD_TruePi0CaloMergedClusterPartConv_InvMass_Pt","ESD_TruePi0CaloMergedClusterPartConv_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTruePi0CaloMergedClusterPartConvInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}} (GeV/c^{2}) #gamma merged cluster, part conv");
					fHistoTruePi0CaloMergedClusterPartConvInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTruePi0CaloMergedClusterPartConvInvMassPt[iCut]);
					fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[iCut] = new TH2F("ESD_TrueEtaCaloMergedClusterPartConv_InvMass_Pt","ESD_TrueEtaCaloMergedClusterPartConv_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[iCut]->SetXTitle("M_{inv,#eta} (GeV/c^{2}) #gamma merged cluster, part conv");
					fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[iCut]->SetYTitle("#eta p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[iCut]);
					
					fHistoTruePrimaryPi0MCPtResolPt[iCut] = new TH2F("ESD_TruePrimaryPi0_MCPt_ResolPt","ESD_TruePrimaryPi0_ResolPt_MCPt",500,0.03,25,1000,-1.,1.);
					fHistoTruePrimaryPi0MCPtResolPt[iCut]->SetXTitle("#pi^{0} p_{T,MC} (GeV/c)");
					fHistoTruePrimaryPi0MCPtResolPt[iCut]->SetYTitle("#pi^{0} (p_{T,rec}-p_{T,MC})/p_{T,MC} ()");	
					fHistoTruePrimaryPi0MCPtResolPt[iCut]->Sumw2();
					SetLogBinningXTH2(fHistoTruePrimaryPi0MCPtResolPt[iCut]);
					fTrueList[iCut]->Add(fHistoTruePrimaryPi0MCPtResolPt[iCut]);
					
					fHistoTruePrimaryEtaMCPtResolPt[iCut]  = new TH2F("ESD_TruePrimaryEta_MCPt_ResolPt","ESD_TruePrimaryEta_ResolPt_MCPt",500,0.03,25,1000,-1.,1.);
					fHistoTruePrimaryEtaMCPtResolPt[iCut]->SetXTitle("#eta p_{T,MC} (GeV/c)");
					fHistoTruePrimaryEtaMCPtResolPt[iCut]->SetYTitle("#eta (p_{T,rec}-p_{T,MC})/p_{T,MC} ()");	
					fHistoTruePrimaryEtaMCPtResolPt[iCut]->Sumw2();
					SetLogBinningXTH2(fHistoTruePrimaryEtaMCPtResolPt[iCut]);
					fTrueList[iCut]->Add(fHistoTruePrimaryEtaMCPtResolPt[iCut]);
					
					fHistoTrueBckGGInvMassPt[iCut] = new TH2F("ESD_TrueBckGG_InvMass_Pt","ESD_TrueBckGG_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTrueBckGGInvMassPt[iCut]->SetXTitle("M_{inv} (GeV/c^{2}) #gamma #gamma no signal");
					fHistoTrueBckGGInvMassPt[iCut]->SetYTitle("#pair p_{T} (GeV/c)");					
					fTrueList[iCut]->Add(fHistoTrueBckGGInvMassPt[iCut]);
					fHistoTrueBckContInvMassPt[iCut] = new TH2F("ESD_TrueBckCont_InvMass_Pt","ESD_TrueBckCont_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTrueBckContInvMassPt[iCut]->SetXTitle("M_{inv} (GeV/c^{2}) contamination");
					fHistoTrueBckContInvMassPt[iCut]->SetYTitle("#pair p_{T} (GeV/c)");					
					fTrueList[iCut]->Add(fHistoTrueBckContInvMassPt[iCut]);
					fHistoTrueK0sWithPi0DaughterMCPt[iCut] = new TH1F("ESD_TrueK0sWithPi0Daughter_MCPt","ESD_TrueK0sWithPi0Daughter_MCPt",250,0,25);
					fHistoTrueK0sWithPi0DaughterMCPt[iCut]->SetXTitle("K^{0}_{s} p_{MC,T} (GeV/c) for K^{0}_{s} where #pi^{0} rec ");			
					fTrueList[iCut]->Add(fHistoTrueK0sWithPi0DaughterMCPt[iCut]);
					fHistoTrueEtaWithPi0DaughterMCPt[iCut] = new TH1F("ESD_TrueEtaWithPi0Daughter_MCPt","ESD_TrueEtaWithPi0Daughter_MCPt",250,0,25);
					fHistoTrueEtaWithPi0DaughterMCPt[iCut]->SetXTitle("#eta p_{MC,T} (GeV/c) for #eta where #pi^{0} rec ");			
					fTrueList[iCut]->Add(fHistoTrueEtaWithPi0DaughterMCPt[iCut]);
					fHistoTrueLambdaWithPi0DaughterMCPt[iCut] = new TH1F("ESD_TrueLambdaWithPi0Daughter_MCPt","ESD_TrueLambdaWithPi0Daughter_MCPt",250,0,25);
					fHistoTrueLambdaWithPi0DaughterMCPt[iCut]->SetXTitle("#Lambda p_{MC,T} (GeV/c) for #Lambda where #pi^{0} rec ");			
					fTrueList[iCut]->Add(fHistoTrueLambdaWithPi0DaughterMCPt[iCut]);
          
					fHistoTruePi0PtY[iCut] = new TH2F("ESD_TruePi0_Pt_Y","ESD_TruePi0_Pt_Y",150,0.03,15.,150,-1.5,1.5);
					fHistoTruePi0PtY[iCut]->SetYTitle("Y_{#pi^{0}}");
					fHistoTruePi0PtY[iCut]->SetXTitle("#pi^{0} p_{T} (GeV/c)");					
					SetLogBinningXTH2(fHistoTruePi0PtY[iCut]);
					fTrueList[iCut]->Add(fHistoTruePi0PtY[iCut]);
					fHistoTrueEtaPtY[iCut] = new TH2F("ESD_TrueEta_Pt_Y","ESD_TrueEta_Pt_Y",150,0.03,15.,150,-1.5,1.5);
					fHistoTrueEtaPtY[iCut]->SetYTitle("Y_{#eta}");
					fHistoTrueEtaPtY[iCut]->SetXTitle("#eta p_{T} (GeV/c)");					
					SetLogBinningXTH2(fHistoTrueEtaPtY[iCut]);
					fTrueList[iCut]->Add(fHistoTrueEtaPtY[iCut]);
					fHistoTruePi0PtAlpha[iCut] = new TH2F("ESD_TruePi0_Pt_Alpha","ESD_TruePi0_Pt_Alpha",150,0.03,15.,100,0,1);
					fHistoTruePi0PtAlpha[iCut]->SetYTitle("#alpha_{#pi^{0}}");
					fHistoTruePi0PtAlpha[iCut]->SetXTitle("#pi^{0} p_{T} (GeV/c)");		
					SetLogBinningXTH2(fHistoTruePi0PtAlpha[iCut]);
					fTrueList[iCut]->Add(fHistoTruePi0PtAlpha[iCut]);
					fHistoTrueEtaPtAlpha[iCut] = new TH2F("ESD_TrueEta_Pt_Alpha","ESD_TrueEta_Pt_Alpha",150,0.03,15.,100,0,1);
					fHistoTrueEtaPtAlpha[iCut]->SetYTitle("#alpha_{#eta}");
					fHistoTrueEtaPtAlpha[iCut]->SetXTitle("#eta p_{T} (GeV/c)");		
					SetLogBinningXTH2(fHistoTrueEtaPtAlpha[iCut]);
					fTrueList[iCut]->Add(fHistoTrueEtaPtAlpha[iCut]);
					
					fHistoTruePi0PtOpenAngle[iCut] = new TH2F("ESD_TruePi0_Pt_OpenAngle","ESD_TruePi0_Pt_OpenAngle",150,0.03,15.,200,0,2*TMath::Pi());
					fHistoTruePi0PtOpenAngle[iCut]->SetYTitle("#theta_{#pi^{0}}");
					fHistoTruePi0PtOpenAngle[iCut]->SetXTitle("#pi^{0} p_{T} (GeV/c)");
					SetLogBinningXTH2(fHistoTruePi0PtOpenAngle[iCut]);
					fTrueList[iCut]->Add(fHistoTruePi0PtOpenAngle[iCut]);
					fHistoTrueEtaPtOpenAngle[iCut] = new TH2F("ESD_TrueEta_Pt_OpenAngle","ESD_TrueEta_Pt_OpenAngle",150,0.03,15.,200,0,2*TMath::Pi());
					fHistoTrueEtaPtOpenAngle[iCut]->SetYTitle("#theta_{#eta}");
					fHistoTrueEtaPtOpenAngle[iCut]->SetXTitle("#eta p_{T} (GeV/c)");
					SetLogBinningXTH2(fHistoTrueEtaPtOpenAngle[iCut]);
					fTrueList[iCut]->Add(fHistoTrueEtaPtOpenAngle[iCut]);
					
					fHistoTrueMotherPi0ConvPhotonEtaPhi[iCut] = new TH2F("ESD_TrueMotherPi0ConvPhoton_Eta_Phi","conv photons for true #pi^{0}",600,0,2*TMath::Pi(),200,-1,1);
					fHistoTrueMotherPi0ConvPhotonEtaPhi[iCut]->SetXTitle("#phi_{#gamma_{conv}} (rad)");
					fHistoTrueMotherPi0ConvPhotonEtaPhi[iCut]->SetYTitle("#eta_{#gamma_{conv}}");
					fTrueList[iCut]->Add(fHistoTrueMotherPi0ConvPhotonEtaPhi[iCut]);
					fHistoTrueMotherEtaConvPhotonEtaPhi[iCut] = new TH2F("ESD_TrueMotherEtaConvPhoton_Eta_Phi","conv photons for true #eta",600,0,2*TMath::Pi(),200,-1,1);
					fHistoTrueMotherEtaConvPhotonEtaPhi[iCut]->SetXTitle("#phi_{#gamma_{conv}} (rad)");
					fHistoTrueMotherEtaConvPhotonEtaPhi[iCut]->SetYTitle("#eta_{#gamma_{conv}}");
					fTrueList[iCut]->Add(fHistoTrueMotherEtaConvPhotonEtaPhi[iCut]);
					
					
					fHistoTruePi0ShowerInvMassPt[iCut] = new TH2F("ESD_TruePi0Shower_InvMass_Pt","ESD_TruePi0Shower_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTruePi0ShowerInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}} (GeV/c^{2})");
					fHistoTruePi0ShowerInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTruePi0ShowerInvMassPt[iCut]);
				
					fHistoTrueEtaShowerInvMassPt[iCut] = new TH2F("ESD_TrueEtaShower_InvMass_Pt","ESD_TrueEtaShower_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTrueEtaShowerInvMassPt[iCut]->SetXTitle("M_{inv,#eta} (GeV/c^{2})");
					fHistoTrueEtaShowerInvMassPt[iCut]->SetYTitle("#eta p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTrueEtaShowerInvMassPt[iCut]);
					
					fHistoTruePi0NoShowerInvMassPt[iCut] = new TH2F("ESD_TruePi0NoShower_InvMass_Pt","ESD_TruePi0NoShower_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTruePi0NoShowerInvMassPt[iCut]->SetXTitle("M_{inv,#pi^{0}} (GeV/c^{2})");
					fHistoTruePi0NoShowerInvMassPt[iCut]->SetYTitle("#pi^{0} p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTruePi0NoShowerInvMassPt[iCut]);
				
					fHistoTrueEtaNoShowerInvMassPt[iCut] = new TH2F("ESD_TrueEtaNoShower_InvMass_Pt","ESD_TrueEtaNoShower_InvMass_Pt",800,0,0.8,250,0,25);
					fHistoTrueEtaNoShowerInvMassPt[iCut]->SetXTitle("M_{inv,#eta} (GeV/c^{2})");
					fHistoTrueEtaNoShowerInvMassPt[iCut]->SetYTitle("#eta p_{T} (GeV/c)");
					fTrueList[iCut]->Add(fHistoTrueEtaNoShowerInvMassPt[iCut]);
				
					fHistoTruePi0OpeningAngleGammaElectron[iCut] = new TH1F("ESD_TruePi0_OpeningAngle_GammaElectron", "ESD_TruePi0_OpeningAngle_GammaElectron",100,0.,TMath::Pi());
					fTrueList[iCut]->Add(fHistoTruePi0OpeningAngleGammaElectron[iCut]);
					
					

				}
			}
		}
	}  

	fVectorDoubleCountTruePi0s.clear();
	fVectorDoubleCountTrueEtas.clear();
	fVectorDoubleCountTrueConvGammas.clear();

    fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
	if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader
  
	if(fV0Reader)
		if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
			if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
				fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
	if(fV0Reader)
		if((AliConvEventCuts*)fV0Reader->GetEventCuts())
			if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
				fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());

    for(Int_t iMatcherTask = 0; iMatcherTask < 3; iMatcherTask++){
      AliCaloTrackMatcher* temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcher_%i",iMatcherTask)));
      if(temp) fOutputContainer->Add(temp->GetCaloTrackMatcherHistograms());
    }
			
	fElecSelector=(AliDalitzElectronSelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("ElectronSelector");
	if(!fElecSelector){printf("Error: No ElectronSelector");return;} // GetV0Reader
		
	if( fElecSelector ){
		if ( ((AliDalitzElectronCuts*)fElecSelector->GetDalitzElectronCuts())->GetCutHistograms() ){
				fOutputContainer->Add( ((AliDalitzElectronCuts*)fElecSelector->GetDalitzElectronCuts())->GetCutHistograms() );
		}
	}  

			
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
	  
		if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
		if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
			fCutFolder[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
		}
		if(!((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))) continue;
		if(((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutHistograms()){
			fCutFolder[iCut]->Add(((AliConversionPhotonCuts*)fGammaCutArray->At(iCut))->GetCutHistograms());
		}
		if(!((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))) continue;
		if(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms()){
			fCutFolder[iCut]->Add(((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutHistograms());
		}
		if(!((AliDalitzElectronCuts*)fElectronCutArray->At(iCut))) continue;
		if( ((AliDalitzElectronCuts*)fElectronCutArray->At(iCut))->GetCutHistograms() ) {
			fCutFolder[iCut]->Add( ((AliDalitzElectronCuts*)fElectronCutArray->At(iCut))->GetCutHistograms() );
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
Bool_t AliAnalysisTaskGammaCaloDalitzV1::Notify()
{
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){        
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnum(fV0Reader->GetPeriodName());
    }  
    
    if(!((AliConvEventCuts*)fEventCutArray->At(iCut))->GetDoEtaShift()){
      fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
      continue; // No Eta Shift requested, continue
    }
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCorrectEtaShiftFromPeriod();
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
void AliAnalysisTaskGammaCaloDalitzV1::UserExec(Option_t *)
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
	
	fInputEvent = InputEvent();
	
	if(fInputEvent->IsA()==AliAODEvent::Class()){
		fInputEvent->InitMagneticField();
	}
	
	fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut
	fSelectorElectronIndex = fElecSelector->GetReconstructedElectronsIndex(); // Electrons from default Cut
	fSelectorPositronIndex = fElecSelector->GetReconstructedPositronsIndex(); // Positrons from default Cut

	

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
		
		Bool_t isRunningEMCALrelAna = kFALSE;
		
		if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;

		Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,isRunningEMCALrelAna);
		
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
		fHistoSPDClusterTrackletBackground[iCut]->Fill(fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),(fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)));
		
		if(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsHeavyIon() == 2)	fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A());
			else fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C());

		if(fIsMC){
		  
		  
			// Process MC Particle
			fStringRecTruePi0s[iCut] = "";
			fStringRecTrueEtas[iCut] = "";
			
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
			
			if(fInputEvent->IsA()==AliESDEvent::Class())
				ProcessMCParticles();
			if(fInputEvent->IsA()==AliAODEvent::Class())
				ProcessAODMCParticles();
			
		}
		
				
		// it is in the loop to have the same conversion cut string (used also for MC stuff that should be same for V0 and Cluster)
		ProcessClusters();  					// process calo clusters
		ProcessPhotonCandidates(); 				// Process this cuts gammas
		ProcessElectronCandidates();

		fHistoNGammaCandidates[iCut]->Fill(fClusterCandidates->GetEntries());
		fHistoNGoodESDTracksVsNGammaCanditates[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),fClusterCandidates->GetEntries());
		
		if(fDoMesonAnalysis){ // Meson Analysis
			CalculatePi0DalitzCandidates(); //Combine gamma candidates

			if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){

				if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
					CalculateBackground(); // Combinatorial Background
					UpdateEventByEventData(); // Store Event for mixed Events
				}

			}

			if(fIsMC && fInputEvent->IsA()==AliAODEvent::Class()){
				ProcessConversionPhotonsForMissingTagsAOD(); //Count missing tags
			} else if (fIsMC && fInputEvent->IsA()==AliESDEvent::Class()){
				ProcessConversionPhotonsForMissingTags(); //Count missing tags
			}
			fVectorDoubleCountTruePi0s.clear();
			fVectorDoubleCountTrueEtas.clear();
		}

		fVectorDoubleCountTrueConvGammas.clear();

		fGammaCandidates->Clear(); // delete this cuts good gammas
		fClusterCandidates->Clear(); // delete cluster candidates
		fVirtualGammaCandidates->Clear(); // delete this cuts good  virtual gammas
	}
	
	if(fIsMC && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
		RelabelAODPhotonCandidates(kFALSE); // Back to ESDMC Label
		fV0Reader->RelabelAODs(kFALSE);
	}
	
	fSelectorElectronIndex.clear();
	fSelectorPositronIndex.clear();
	
	PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::ProcessClusters()
{
	
	Int_t nclus = 0;
	nclus = fInputEvent->GetNumberOfCaloClusters();
	
	
	if(nclus == 0)	return;
	
	
	// vertex
	Double_t vertex[3] = {0};
	InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
	
	
	// Loop over EMCal clusters
	for(Long_t i = 0; i < nclus; i++){
		
		AliVCluster* clus = NULL;
		if(fInputEvent->IsA()==AliESDEvent::Class()) clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
		else if(fInputEvent->IsA()==AliAODEvent::Class()) clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));

		if ( !clus ) continue;
		if ( !((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC,1.,i) ){ delete clus; continue;}
		
		// TLorentzvector with cluster
		TLorentzVector clusterVector;
		clus->GetMomentum(clusterVector,vertex);
		
		TLorentzVector* tmpvec = new TLorentzVector();
		tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());
		
		// convert to AODConversionPhoton
		AliAODConversionPhoton *PhotonCandidate=new AliAODConversionPhoton(tmpvec);
		if(!PhotonCandidate){ delete clus; delete tmpvec; continue;}
		
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
					if (k< 50)PhotonCandidate->SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
// 					Int_t pdgCode = fMCEvent->Particle(mclabelsCluster[k])->GetPdgCode();
// 					cout << "label " << k << "\t" << mclabelsCluster[k] << " pdg code: " << pdgCode << endl;
				}	
			}
		}
		
		fIsFromMBHeader = kTRUE; 
		fIsOverlappingWithOtherHeader = kFALSE;
		// test whether largest contribution to cluster orginates in added signals
		
		if ( fIsMC && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() > 0 ){
		
            if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 0 ) fIsFromMBHeader = kFALSE;
		
		  
			if (clus->GetNLabels()>1){
				Int_t* mclabelsCluster = clus->GetLabels();
				for (Int_t l = 1; l < (Int_t)clus->GetNLabels(); l++ ){
                    if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(mclabelsCluster[l], fMCEvent, fInputEvent) == 0) fIsOverlappingWithOtherHeader = kTRUE;
				}	
			}	
		}

		if(fIsMC){
			if(fInputEvent->IsA()==AliESDEvent::Class()){
				ProcessTrueClusterCandidates(PhotonCandidate);
			} else {
				ProcessTrueClusterCandidatesAOD(PhotonCandidate);
			}
		}

		if (fIsFromMBHeader && fIsOverlappingWithOtherHeader) fHistoClusOverlapHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
		if (fIsFromMBHeader && !fIsOverlappingWithOtherHeader){
			fHistoClusGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
			fClusterCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas
		} else{
			delete PhotonCandidate;
		}
		
		delete clus;
		delete tmpvec;
	}
	
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::ProcessTrueClusterCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
		
	TParticle *Photon = NULL;
	if (!TruePhotonCandidate->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set task will abort");
	fHistoTrueNLabelsInClus[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels());
	
    if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0)Photon = fMCEvent->Particle(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
		else return;
		
	if(Photon == NULL){
	//    cout << "no photon" << endl;
		return;
	}

    TruePhotonCandidate->SetCaloPhotonMCFlags(fMCEvent,kFALSE);
	
	// True Photon
	if(fIsFromMBHeader && !fIsOverlappingWithOtherHeader){
	  
		if (TruePhotonCandidate->IsLargestComponentPhoton() || TruePhotonCandidate->IsLargestComponentElectron() ) {
		  
		  fHistoTrueClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
		  
		  
		  
		  const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
		  Double_t mcProdVtxX 	= primVtxMC->GetX();
		  Double_t mcProdVtxY 	= primVtxMC->GetY();
		  Double_t mcProdVtxZ 	= primVtxMC->GetZ();	
          Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, TruePhotonCandidate->GetCaloPhotonMCLabel(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);

	
	
		  if(isPrimary){
		    
		      // Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
		      
			fHistoTruePrimaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			fHistoTruePrimaryClusGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt()); // Allways Filled
			
			
			Int_t gammaMotherLabel = -1;
			
			if (TruePhotonCandidate->IsLargestComponentPhoton()){ // for photons its the direct mother 

				gammaMotherLabel=Photon->GetMother(0);
				
					
			} else if ( TruePhotonCandidate->IsLargestComponentElectron() ){  // for electrons its either the direct mother or for conversions the grandmother

                    if (TruePhotonCandidate->IsConversion()) gammaMotherLabel=fMCEvent->Particle(Photon->GetMother(0))->GetMother(0);
					else gammaMotherLabel=Photon->GetMother(0); 
					
			} 
			
			if( gammaMotherLabel > -1  ){
			  
			      Int_t motherFromShower = -1;
			  
			         
                   TParticle* motherGamma = fMCEvent->Particle( gammaMotherLabel );
			       
			       //if( IsDalitz ( motherGamma ) == kTRUE ){
				 if( motherGamma->GetPdgCode() == 111 ) {
				    
				    if( IsDalitz( motherGamma ) == kTRUE ){
				 
				      fHistoTruePi0DalitzClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				      fHistoTruePi0DalitzClusGammaMCPt[fiCut]->Fill(Photon->Pt());
				    }
				  
				 
			       } else if ( ( motherFromShower = FindMotherOfPhoton( gammaMotherLabel ) ) > -1 ){
				 
				 
				 
				 
                   TParticle* motherGammaShower = fMCEvent->Particle( motherFromShower );
				   
				   
				   	   
				   
				       if( motherGammaShower->GetPdgCode() == 111 ) {
					    if( IsDalitz ( motherGammaShower ) ) {
					  					  
					    fHistoTruePi0DalitzClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
					    fHistoTruePi0DalitzAllClusGammaPt[fiCut]->Fill( TruePhotonCandidate->Pt() );
					    fHistoTruePi0DalitzClusGammaMCPt[fiCut]->Fill(Photon->Pt());
				  
					  
					  }
				 
				       }	    
			       }
			}
			
		     }	
		  
		}		  
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
                fHistoTrueClusConvGammaMCPt[fiCut]->Fill(((TParticle*)fMCEvent->Particle(Photon->GetMother(0)))->Pt());
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
			fHistoTrueClusNParticles[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels());
		}
	}

	
	return;
}


//________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::ProcessTrueClusterCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate)
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
	TruePhotonCandidate->SetCaloPhotonMCFlagsAOD(fInputEvent, kFALSE);
	fHistoTrueNLabelsInClus[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels());	
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
			fHistoTrueClusNParticles[fiCut]->Fill(TruePhotonCandidate->GetNCaloPhotonMCLabels());
		}
	}

	// True Photon
	if(fIsFromMBHeader){
		fHistoTrueClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
 		if (fDoPhotonQA > 0) fHistoTrueConvGammaEta[fiCut]->Fill(TruePhotonCandidate->Eta());
	}

	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();
	Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);

	if(isPrimary){
		// Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
		if(fIsFromMBHeader){
			fHistoTruePrimaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			fHistoTruePrimaryClusGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt()); // Allways Filled
		}
	}	
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::ProcessAODMCParticles()
{
	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();

	TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	
	// Loop over all primary MC particle
	for(Int_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {
		
		AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
		if (!particle) continue;

		Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, particle, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
		if (!isPrimary) continue;
		
		Int_t isMCFromMBHeader = -1;
		if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
            isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
			if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
		}
		
		if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) continue;
		if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,AODMCTrackArray,kFALSE)){
			fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt()); // All MC Gamma
            if (TMath::Abs(particle->Eta()) < 0.66 ){
				if (particle->Phi() > 1.39626 && particle->Phi() < 3.125) fHistoMCAllGammaEMCALAccPt[fiCut]->Fill(particle->Pt());
			}
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
		if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,AODMCTrackArray,kTRUE)){
			Double_t rConv = 0;
			for(Int_t daughterIndex=particle->GetDaughter(0);daughterIndex<=particle->GetDaughter(1);daughterIndex++){
				AliAODMCParticle *tmpDaughter = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(daughterIndex));
				if(!tmpDaughter) continue;
                if(TMath::Abs(tmpDaughter->GetPdgCode()) == 11){
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
                if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
					if (particle->Pt()>0.005){
						weightedK0s= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, 0x0, fInputEvent);
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
                if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
					if (particle->Pt()>0.005){
						weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, 0x0, fInputEvent);
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
				
				Double_t alpha = -1;
				if (particle->GetPdgCode() == 111 || particle->GetPdgCode() == 221){
					alpha = TMath::Abs((daughter0->E() - daughter1->E()))/(daughter0->E() + daughter1->E());
				}

				
				if(particle->GetPdgCode() == 111){
					fHistoMCPi0GGPt[fiCut]->Fill(particle->Pt(),weighted); // All MC Pi0
					fHistoMCPi0GGWOWeightPt[fiCut]->Fill(particle->Pt());
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
				
				// Check the acceptance for both gammas
				if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter0,AODMCTrackArray,kFALSE) &&
				((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter1,AODMCTrackArray,kFALSE)  &&
				((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
				((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){
					
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
void AliAnalysisTaskGammaCaloDalitzV1::ProcessMCParticles()
{
	
	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 		= primVtxMC->GetX();
	Double_t mcProdVtxY 		= primVtxMC->GetY();
	Double_t mcProdVtxZ 		= primVtxMC->GetZ();
    
  
  
	// Loop over all primary MC particle
    for(Int_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
	  
	  
	 
      if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){
		
        TParticle* particle = (TParticle *)fMCEvent->Particle(i);
		if (!particle) continue;
		
		
		
		Int_t isMCFromMBHeader = -1;
		if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
		  
            isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
			if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
			
		}
		
		
			
				
		if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) continue;
		
		
        if (  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(particle,fMCEvent) == kTRUE ){
		  
			fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt()); // All MC Gamma
			
			
			if( IsPi0DalitzDaughter ( i ) == kTRUE ) {
			  
			       fHistoMCAllGammaPi0Pt[fiCut]->Fill(particle->Pt());
			}
		  
		}
		
		
        if(((AliDalitzElectronCuts*)fElectronCutArray->At(fiCut))->ElectronIsSelectedMC(i,fMCEvent) == kTRUE){
		  
		Double_t deltaX = particle->Vx() - mcProdVtxX;
		Double_t deltaY = particle->Vy() - mcProdVtxY;
		Double_t deltaZ = particle->Vz() - mcProdVtxZ;

		Double_t realRadius3D = TMath::Sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ);
		  
		  
		  	if( particle->GetPdgCode() == -11) {
				fHistoMCAllPositronsPt[fiCut]->Fill(particle->Pt()); // All positrons
				
								  
				
				if ( particle->GetMother(0) > -1 ) {
				    
				      if( IsPi0DalitzDaughter( i ) == kTRUE ){
					fHistoMCDecayPositronPi0Pt[fiCut]->Fill(particle->Pt()); //Positrons from Pi0->Dalitz
					
                      if( i >= fMCEvent->GetNumberOfPrimaries() ){
				  
						fHistoMCDecayNoPrimPositronPi0DalitzR[fiCut]->Fill(realRadius3D);
						fHistoMCDecayNoPrimPositronPi0DalitzID[fiCut]->Fill( particle->GetUniqueID());
				      
					  }
					
				      }
				  
				}
			}
			if( particle->GetPdgCode() ==  11){
			  
				fHistoMCAllElectronsPt[fiCut]->Fill(particle->Pt()); // All electrons
				
				if (  particle->GetMother(0) > -1 ) {
				    
				      if( IsPi0DalitzDaughter( i ) == kTRUE ){
					fHistoMCDecayElectronPi0Pt[fiCut]->Fill(particle->Pt()); //Electrons from Pi0->Dalitz
				      }
				      
                      if( i >= fMCEvent->GetNumberOfPrimaries() ){
				  
						fHistoMCDecayNoPrimElectronPi0DalitzR[fiCut]->Fill(realRadius3D);
						fHistoMCDecayNoPrimElectronPi0DalitzID[fiCut]->Fill( particle->GetUniqueID());
				  
					  }
				  
				}

			}		
		}
		
		
        if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kFALSE)){
		  
			
			
            if ( TMath::Abs(particle->Eta()) < 0.66 ){
				if (particle->Phi() > 1.39626 && particle->Phi() < 3.125) fHistoMCAllGammaEMCALAccPt[fiCut]->Fill(particle->Pt());
			}	
			
			if(particle->GetMother(0) >-1){ // Meson Decay Gamma
                switch(fMCEvent->Particle(particle->GetMother(0))->GetPdgCode()){
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
        if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kTRUE)){
			fHistoMCConvGammaPt[fiCut]->Fill(particle->Pt());
			if (fDoPhotonQA > 0){
                fHistoMCConvGammaR[fiCut]->Fill(((TParticle*)fMCEvent->Particle(particle->GetFirstDaughter()))->R());
				fHistoMCConvGammaEta[fiCut]->Fill(particle->Eta());
			}
		} // Converted MC Gamma
		
		if(fDoMesonAnalysis){
			
			
			if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
                ->MesonIsSelectedMC(particle,fMCEvent,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
				
				Float_t weighted= 1;
			
				if( ((AliDalitzElectronCuts*) fElectronCutArray->At(fiCut))->DoWeights() ) {
				  
                  if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
				  
					if (particle->Pt()>0.005){
                        weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCEvent, fInputEvent);
						
					}
				  }
				}
				
				if(particle->GetPdgCode() == 111){
				  
					fHistoMCPi0GGPt[fiCut]->Fill(particle->Pt(),weighted); // All MC Pi0
					fHistoMCPi0GGWOWeightPt[fiCut]->Fill(particle->Pt());
					
					
				} else if(particle->GetPdgCode() == 221){
				  
					fHistoMCEtaGGPt[fiCut]->Fill(particle->Pt(),weighted); // All MC Eta
					fHistoMCEtaGGWOWeightPt[fiCut]->Fill(particle->Pt());
					
				}
								
				
			}
			
			
			Int_t labelgamma	= -1;
			Int_t labelelectron 	= -1;
			Int_t labelpositron 	= -1;


            if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMCDalitz(particle,fMCEvent,labelelectron,labelpositron,labelgamma,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
			
			      Float_t weighted= 1;
			      
			      if( ((AliDalitzElectronCuts*) fElectronCutArray->At(fiCut))->DoWeights() ) { 
                  if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent,fInputEvent) ) {
					  if (particle->Pt()>0.005){
                        weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCEvent,fInputEvent);
					  }
				  }
			      }
			    
			      if(particle->GetPdgCode() == 111){
				    fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(), weighted); // All MC Pi0 Dalitz
				    fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt() );  //All MC Pi0 Dalitz without weights
			      }
			      if(particle->GetPdgCode() == 221){
				    fHistoMCEtaPt[fiCut]->Fill(particle->Pt(), weighted); // All MC Eta
				    fHistoMCEtaWOWeightPt[fiCut]->Fill( particle->Pt() ); //All MC Eta Dalitz without weights
			      }
			
			      // Check the acceptance for gamma and electrons
			
                  TParticle *gamma    = fMCEvent->Particle(labelgamma);
                  TParticle *electron = fMCEvent->Particle(labelelectron);
                  TParticle *positron = fMCEvent->Particle(labelpositron);
			
	  
			      
			      
			      
                  Bool_t kDaughGammaIsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelgamma,     mcProdVtxX, mcProdVtxY, mcProdVtxZ);
                  Bool_t kDaughElectIsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelelectron,  mcProdVtxX, mcProdVtxY, mcProdVtxZ);
                  Bool_t kDaughPositIsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelpositron,  mcProdVtxX, mcProdVtxY, mcProdVtxZ);
					
			
			
			      if( kDaughElectIsPrim && kDaughPositIsPrim && kDaughGammaIsPrim &&				
                 ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma,fMCEvent) &&
                 ((AliDalitzElectronCuts*)fElectronCutArray->At(fiCut))->ElectronIsSelectedMC(labelelectron,fMCEvent) &&
                 ((AliDalitzElectronCuts*)fElectronCutArray->At(fiCut))->ElectronIsSelectedMC(labelpositron,fMCEvent) ) {
				
				
				        TLorentzVector TLVEpos,TLVEneg,TLVDalitz;
					
					Double_t electronMass = TDatabasePDG::Instance()->GetParticle(  ::kElectron   )->Mass();					
									        
				        TLVEpos.SetXYZM(positron->Px(),positron->Py(),positron->Pz(),electronMass);
					TLVEneg.SetXYZM(electron->Px(),electron->Py(),electron->Pz(),electronMass);
					TVector3 V3gamma(gamma->Px(),gamma->Py(),gamma->Pz());
					Double_t angleGammaEpos = V3gamma.Angle(TLVEpos.Vect());
					Double_t angleGammaEneg = V3gamma.Angle(TLVEneg.Vect());
				
				
				
				
				
				
				  if( particle->GetPdgCode() == 111 ){ 
				    
				    
					fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt() , weighted); // MC Pi0Dalitz with gamma and e+e- in acc
					
					if( fDoMesonQA > 0 ) {
					  
					    fHistoMCPi0InAccOpeningAngleGammaElectron[fiCut]->Fill(angleGammaEpos);
					    fHistoMCPi0InAccOpeningAngleGammaElectron[fiCut]->Fill(angleGammaEneg);
					    
					}
					
				
				  }
				  if( particle->GetPdgCode() == 221 ){ 
				    
					fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),  weighted); // MC EtaDalitz with gamma and e+e- in acc
					
				  }
			  }	
		  }
	    }
	}
   }
}


//________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::CalculatePi0DalitzCandidates(){

	// Conversion Gammas
	
	if( fClusterCandidates->GetEntries() > 0 && fVirtualGammaCandidates->GetEntries() > 0 ){
	  
		vector<Bool_t> lGoodVirtualGamma(fVirtualGammaCandidates->GetEntries(), kFALSE);
		
		for(Int_t GammaIndex=0; GammaIndex<fClusterCandidates->GetEntries(); GammaIndex++){
		  
		    AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(GammaIndex));
		    
		      if (gamma==NULL) continue;
		      
		      
			for(Int_t virtualGammaIndex=0;virtualGammaIndex<fVirtualGammaCandidates->GetEntries();virtualGammaIndex++){
			  
			        Bool_t matched = kFALSE; 
				
				AliAODConversionPhoton *Vgamma=dynamic_cast<AliAODConversionPhoton*>(fVirtualGammaCandidates->At(virtualGammaIndex));
				if (Vgamma==NULL) continue;
				
				//Check for same Electron ID
			     	if ( gamma->GetIsCaloPhoton() ){
				  
					AliVCluster* cluster = fInputEvent->GetCaloCluster(gamma->GetCaloClusterRef());
					matched = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(Vgamma,cluster,fInputEvent );
					
				}	
				
				
				AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma,Vgamma);
				pi0cand->SetLabels(GammaIndex,virtualGammaIndex);
				
				
				
				
				if( (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())) ){
				  
				Int_t zbin = 0;
				Int_t mbin = 0;
				
				Double_t angleGammaPositron = -999;
				Double_t angleGammaElectron = -999;
				
				
				 if ( fDoMesonQA > 0 ) {
				
					AliESDtrack* positronVgamma = (AliESDtrack*)fInputEvent->GetTrack( Vgamma->GetTrackLabelPositive() );
					Double_t momPositron[3];
					positronVgamma->GetConstrainedPxPyPz(momPositron);
					
							
					AliESDtrack* electronVgamma = (AliESDtrack*)fInputEvent->GetTrack( Vgamma->GetTrackLabelNegative() );
					Double_t momElectron[3];
					electronVgamma->GetConstrainedPxPyPz(momElectron);
					
							
					TVector3 vGamma(gamma->GetPx(),gamma->GetPy(),gamma->GetPz());;
					TVector3 vPositron(momPositron[0],momPositron[1],momPositron[2]);
					TVector3 vElectron(momElectron[0],momElectron[1],momElectron[2]);
					
					angleGammaPositron = vGamma.Angle(vPositron);
					angleGammaElectron = vGamma.Angle(vElectron);
					
				 }
					
				
				  
				    if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoBGCalculation()){
				      
				      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->BackgroundHandlerType() == 0){
					
					    zbin = fBGClusHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
					  
					  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
					    mbin = fBGClusHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
					  } else {
					    mbin = fBGClusHandler[fiCut]->GetMultiplicityBinIndex(fClusterCandidates->GetEntries());
					  }
				      } 
				    }
				  				  
				    if(  ((AliDalitzElectronCuts*)fElectronCutArray->At(fiCut))->DoMassCut() == kTRUE ) {
				      if( ((AliDalitzElectronCuts*) fElectronCutArray->At(fiCut))->MassCut( pi0cand->Pt() , Vgamma->M() ) == kTRUE ){
				  				    
					      if (matched){
						
						fHistoMotherMatchedInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
						
					      } else {
						
						fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
						
                            if(TMath::Abs(pi0cand->GetAlpha())<0.1)
							fHistoMotherInvMassEalpha[fiCut]->Fill(pi0cand->M(),pi0cand->E());
							
							
							if( fDoMesonQA > 0 ) {  
					  
							fHistoMotherInvMassOpeningAngleGammaElectron[fiCut]->Fill( angleGammaPositron );
							fHistoMotherInvMassOpeningAngleGammaElectron[fiCut]->Fill( angleGammaElectron );
							
							
							}
							  
						
							if( fDoTHnSparse ) {
							  Double_t sparesFill[4] = {pi0cand->M(),pi0cand->Pt(),(Double_t)zbin,(Double_t)mbin};
							  fSparseMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
							}
						
						if(fMCEvent){
							  ProcessTrueMesonCandidates(pi0cand,Vgamma,gamma,matched);
						}
						
					      }
				    	}
				    } else {
				      
					      if (matched){
						fHistoMotherMatchedInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
					      } else {
						
						    fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
						
                            if(TMath::Abs( pi0cand->GetAlpha()) < 0.1 )
						    fHistoMotherInvMassEalpha[fiCut]->Fill(pi0cand->M(),pi0cand->E());
					      
						    if ( fDoTHnSparse ) {
							Double_t sparesFill[4] = {pi0cand->M(),pi0cand->Pt(),(Double_t)zbin,(Double_t)mbin};
							fSparseMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
						    }
						    
						    if( fDoMesonQA > 0 ) {  
					  
							fHistoMotherInvMassOpeningAngleGammaElectron[fiCut]->Fill( angleGammaPositron );
							fHistoMotherInvMassOpeningAngleGammaElectron[fiCut]->Fill( angleGammaElectron );
							
						    }
						    
						
						    if(fMCEvent){
							  ProcessTrueMesonCandidates(pi0cand,Vgamma,gamma,matched);
						    }
					      }				      
				      
				    }
				}
				delete pi0cand;
				pi0cand=0x0;
			}
		}
	}
}

//______________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueVirtualGammaCandidate, AliAODConversionPhoton *TrueGammaCandidate1, Bool_t matched)
{
	// Process True Mesons
	
		const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
		Double_t mcProdVtxX 		= primVtxMC->GetX();
		Double_t mcProdVtxY 		= primVtxMC->GetY();
		Double_t mcProdVtxZ 		= primVtxMC->GetZ();	
		
		
		Bool_t isTruePi0 = kFALSE;
		Bool_t isTrueEta = kFALSE;
		
		
        Int_t virtualGammaMCLabel = TrueVirtualGammaCandidate->GetMCParticleLabel(fMCEvent);
		Int_t virtualGammaMotherLabel = -1;
		Int_t virtualGamma = -1;
		
		Bool_t motherFromShower  = kFALSE;
		
		Double_t angleGammaPositron = -999;
		Double_t angleGammaElectron = -999;
		
		
	       if( fDoMesonQA > 0 ) {
		
		AliESDtrack* positronVgamma = (AliESDtrack*)fInputEvent->GetTrack( TrueVirtualGammaCandidate->GetTrackLabelPositive() );
		Double_t momPositron[3];
		positronVgamma->GetConstrainedPxPyPz(momPositron);
					
							
		AliESDtrack* electronVgamma = (AliESDtrack*)fInputEvent->GetTrack( TrueVirtualGammaCandidate->GetTrackLabelNegative() );
		Double_t momElectron[3];
		electronVgamma->GetConstrainedPxPyPz(momElectron);
					
							
		TVector3 vGamma(TrueGammaCandidate1->GetPx(),TrueGammaCandidate1->GetPy(),TrueGammaCandidate1->GetPz());;
		TVector3 vPositron(momPositron[0],momPositron[1],momPositron[2]);
		TVector3 vElectron(momElectron[0],momElectron[1],momElectron[2]);
		
		angleGammaPositron = vGamma.Angle(vPositron);
		angleGammaElectron = vGamma.Angle(vElectron);
		
		}
		
		
	
		if( virtualGammaMCLabel != -1 ){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
			// Daughters Gamma 1
            TParticle * negativeMC = (TParticle*)TrueVirtualGammaCandidate->GetNegativeMCDaughter(fMCEvent);
            TParticle * positiveMC = (TParticle*)TrueVirtualGammaCandidate->GetPositiveMCDaughter(fMCEvent);
            TParticle * virtualGammaMotherMC = (TParticle*)fMCEvent->Particle(virtualGammaMCLabel);

			if( TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11) {  // Electrons ...

				if( virtualGammaMotherMC->GetPdgCode() != 22 ){
					virtualGammaMotherLabel=virtualGammaMCLabel;
					virtualGamma = 1;
					
				} else if (negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
					virtualGammaMotherLabel=virtualGammaMotherMC->GetMother(0);
					virtualGamma = 0; //no virtual gamma
				}
			}
		}
		
		if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
		
		Int_t gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0); 	// get most probable MC label
		Int_t gamma1MotherLabel = -1;
		// check if 

		if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother

			// Daughters Gamma 1
            TParticle * gammaMC1 = (TParticle*)fMCEvent->Particle(gamma1MCLabel);
			
			if ( TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron() ){		// largest component is electro magnetic
				// get mother of interest (pi0 or eta)
				if (TrueGammaCandidate1->IsLargestComponentPhoton()){														// for photons its the direct mother 
					gamma1MotherLabel=gammaMC1->GetMother(0);
				} else if (TrueGammaCandidate1->IsLargestComponentElectron()){ 												// for electrons its either the direct mother or for conversions the grandmother
                    if (TrueGammaCandidate1->IsConversion()) gamma1MotherLabel=fMCEvent->Particle(gammaMC1->GetMother(0))->GetMother(0);
					else gamma1MotherLabel=gammaMC1->GetMother(0); 
				}
				
				
				if( gamma1MotherLabel > -1  ){
				  
                    Int_t pdgCodeMother = ((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode();
				 
				    if( pdgCodeMother != 111 && pdgCodeMother != 221 ){ //The mother is not and eta and MC
				
					  gamma1MotherLabel = FindMotherOfPhoton( gamma1MotherLabel );
					  motherFromShower = kTRUE;
									      
				    }
				}
				
				
				
			} else {
				if (fDoMesonQA > 0) fHistoTrueMotherCaloEMNonLeadingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	
		}
				
		if(virtualGammaMotherLabel >= 0 && virtualGammaMotherLabel==gamma1MotherLabel){
		  
            if(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
				isTruePi0=kTRUE;
				if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma1MotherLabel)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
			
            if(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->GetPdgCode() == 221){
				isTrueEta=kTRUE;
				if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,gamma1MotherLabel)) fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
			
		}
		
		if(isTruePi0 || isTrueEta){// True Pion or Eta

						  
			  if( virtualGamma == 1 ){
			  			  
				if( !matched ) {
			    
				    if( isTruePi0 ) {
				      
					    fHistoTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					    if( fDoMesonQA > 0 ) {
					      fHistoTruePi0OpeningAngleGammaElectron[fiCut]->Fill(angleGammaPositron);
					      fHistoTruePi0OpeningAngleGammaElectron[fiCut]->Fill(angleGammaElectron);
					    }
				    }
				    if (isTrueEta)fHistoTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				
				}
							
								
				if (fDoMesonQA > 0){
				  
				  
				      if( motherFromShower == kTRUE  && !matched ){
					
					  if (isTruePi0)fHistoTruePi0ShowerInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					  if (isTrueEta)fHistoTrueEtaShowerInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					
				      } else if ( motherFromShower == kFALSE  && !matched ){
					
					  if (isTruePi0)fHistoTruePi0NoShowerInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					  if (isTrueEta)fHistoTrueEtaNoShowerInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					 					
				      }
				  
				      if (TrueGammaCandidate1->IsLargestComponentPhoton() && !matched){
					if(isTruePi0) fHistoTruePi0CaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if(isTrueEta) fHistoTrueEtaCaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				      }	
				      if (TrueGammaCandidate1->IsLargestComponentElectron() && !matched){
					if(isTruePi0) fHistoTruePi0CaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if(isTrueEta) fHistoTrueEtaCaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				      }	
				      if (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion() ){
					if (isTruePi0 && !matched)fHistoTruePi0CaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (isTrueEta && !matched)fHistoTrueEtaCaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if ((TrueVirtualGammaCandidate->GetMCLabelPositive() == gamma1MCLabel || TrueVirtualGammaCandidate->GetMCLabelNegative() == gamma1MCLabel) && isTruePi0){
						fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					}	
					if ((TrueVirtualGammaCandidate->GetMCLabelPositive() == gamma1MCLabel || TrueVirtualGammaCandidate->GetMCLabelNegative() == gamma1MCLabel) && isTrueEta){
						fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					}	
				      }	
				      if ((TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate1->IsDalitzMerged()) && !matched ){
					if (isTruePi0 )fHistoTruePi0CaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (isTrueEta )fHistoTrueEtaCaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				      }	
				      if (TrueGammaCandidate1->IsMergedPartConv() && !matched){
					if (isTruePi0) fHistoTruePi0CaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					if (isTrueEta) fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				      }	
				}
				
				if ( !matched ) {
				  
                    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, gamma1MotherLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
				
			
				    if(!isPrimary){ // Secondary Meson

                    Int_t secMotherLabel = ((TParticle*)fMCEvent->Particle(virtualGammaMotherLabel))->GetMother(0);
					Float_t weightedSec = 1;
					
                    if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, fMCEvent, fInputEvent) && fMCEvent->Particle(secMotherLabel)->GetPdgCode()==310){
                        weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, fMCEvent, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
						
					}
					
					if (isTruePi0){
						fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
						fHistoTrueSecondaryPi0PhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate1->Pt(),weightedSec);
						if(CheckIfContainedInString(fStringRecTruePi0s[fiCut],virtualGammaMotherLabel)){
							fHistoTrueSecondaryPi0DCPtconv[fiCut]->Fill(TrueGammaCandidate1->Pt(),weightedSec);
						}	
					}	
					
					if (secMotherLabel >-1){
                        if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==310 && isTruePi0 ){
							fHistoTrueSecondaryPi0FromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
                            if (fDoMesonQA > 0)fHistoTrueK0sWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt());
						}
                        if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==221 && isTruePi0){
							fHistoTrueSecondaryPi0FromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
                            if (fDoMesonQA > 0)fHistoTrueEtaWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt());
						}
                        if(fMCEvent->Particle(secMotherLabel)->GetPdgCode()==3122 && isTruePi0){
							fHistoTrueSecondaryPi0FromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
                            if (fDoMesonQA > 0)fHistoTrueLambdaWithPi0DaughterMCPt[fiCut]->Fill(fMCEvent->Particle(secMotherLabel)->Pt());
						}
					}
				} else { // Only primary pi0 for efficiency calculation
					Float_t weighted= 1;
                    if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, fMCEvent, fInputEvent)){
                        if (((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt()>0.005){
                            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, fMCEvent, fInputEvent);
							//                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
						}
					}
					if (isTruePi0){
						fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
						fHistoTruePrimaryPi0PhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate1->Pt(),weighted);
						fHistoTruePrimaryPi0W0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
						fProfileTruePrimaryPi0WeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
						if(CheckIfContainedInString(fStringRecTruePi0s[fiCut],virtualGammaMotherLabel)){
							fHistoTruePrimaryPi0DCPtconv[fiCut]->Fill(TrueGammaCandidate1->Pt(),weighted);
						}	

					} else if (isTrueEta) {
						fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
						fHistoTruePrimaryEtaPhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate1->Pt(),weighted);
						fHistoTruePrimaryEtaW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
						fProfileTruePrimaryEtaWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
						if(CheckIfContainedInString(fStringRecTrueEtas[fiCut],virtualGammaMotherLabel)){
							fHistoTruePrimaryEtaDCPtconv[fiCut]->Fill(TrueGammaCandidate1->Pt(),weighted);
						}							
					}	
						
					if (fDoMesonQA > 0){
						if(isTruePi0){ // Only primary pi0 for resolution
                            fHistoTruePrimaryPi0MCPtResolPt[fiCut]->Fill(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt())/((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),weighted);
						}
						if (isTrueEta){ // Only primary eta for resolution
                            fHistoTruePrimaryEtaMCPtResolPt[fiCut]->Fill(((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt())/((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt(),weighted);
						}
					}
				    } 
				 }
			    } else if (	virtualGamma == 0 ) {
			      
			      	if( !matched ) {   
			      
				  if (isTruePi0)fHistoTruePi0GGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				  if (isTrueEta)fHistoTrueEtaGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			
                if(virtualGammaMotherLabel >= fMCEvent->GetNumberOfPrimaries()){ // Secondary Meson

                    Int_t secMotherLabel = ((TParticle*)fMCEvent->Particle(virtualGammaMotherLabel))->GetMother(0);
					Float_t weightedSec= 1;
					
                    if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, fMCEvent, fInputEvent) && fMCEvent->Particle(secMotherLabel)->GetPdgCode()==310){
                        weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, fMCEvent, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
						//cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
					}
					if (isTruePi0){
						fHistoTrueSecondaryPi0GGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
					}	
					
				} else { // Only primary pi0 for efficiency calculation
					Float_t weighted= 1;
                    if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, fMCEvent, fInputEvent)){
                        if (((TParticle*)fMCEvent->Particle(gamma1MotherLabel))->Pt()>0.005){
                            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, fMCEvent, fInputEvent);
							//                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
						}
					}
					if (isTruePi0){
						fHistoTruePrimaryPi0GGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
						
					} else if (isTrueEta) {
						fHistoTruePrimaryEtaGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
					}	
						
				}
			      }
			    
			    }  
			      
			} else if( !isTruePi0 && !isTrueEta ){ // Background

			    if (fDoMesonQA > 0){
				  if(virtualGammaMotherLabel>-1 && gamma1MotherLabel>-1 && virtualGamma == 0){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
					fHistoTrueBckGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				  } else { // No photon or without mother
					fHistoTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				  }
			    }
			}
		
		
		if (isTrueEta && !matched){
			fStringRecTrueEtas[fiCut].Append(Form("%i,",virtualGammaMotherLabel));
		}
		if (isTruePi0 && !matched){
			fStringRecTruePi0s[fiCut].Append(Form("%i,",virtualGammaMotherLabel));
		}	

}
//______________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::ProcessTrueMesonCandidatesAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1, Bool_t matched)
{
	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();

	// Process True Mesons
	TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	Bool_t isTruePi0 = kFALSE;
	Bool_t isTrueEta = kFALSE;
	
	AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
	AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));
	
	Int_t gamma0MCLabel = -1;
	Int_t gamma0MotherLabel = -1;
	if(!positiveMC||!negativeMC)
		return;
	
	if (TrueGammaCandidate0->IsTrueConvertedPhoton()){
		gamma0MCLabel = positiveMC->GetMother();
		AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
		gamma0MotherLabel=gammaMC0->GetMother();
	}	

	if (!TrueGammaCandidate1->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set. Aborting");
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
			if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma1MotherLabel)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
		}
		if(((AliAODMCParticle*)AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 221){
			isTrueEta=kTRUE;
			if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,gamma1MotherLabel)) fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
		}
	}
	
	if(isTruePi0 || isTrueEta){// True Pion or Eta
		if (!matched){
			if (isTruePi0)fHistoTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			if (isTrueEta)fHistoTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
		}	
		if (fDoMesonQA > 0){
			if (TrueGammaCandidate1->IsLargestComponentPhoton() && !matched){
				if (isTruePi0) fHistoTruePi0CaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta) fHistoTrueEtaCaloPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	
			if (TrueGammaCandidate1->IsLargestComponentElectron() && !matched) {
				if (isTruePi0) fHistoTruePi0CaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta) fHistoTrueEtaCaloElectronInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	
			if (TrueGammaCandidate1->IsLargestComponentElectron() && TrueGammaCandidate1->IsConversion()){
			  
				if (isTruePi0 && !matched)fHistoTruePi0CaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta && !matched)fHistoTrueEtaCaloConvertedPhotonInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if ((TrueGammaCandidate0->GetMCLabelPositive() == gamma1MCLabel || TrueGammaCandidate0->GetMCLabelNegative() == gamma1MCLabel) && isTruePi0)
					fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if ((TrueGammaCandidate0->GetMCLabelPositive() == gamma1MCLabel || TrueGammaCandidate0->GetMCLabelNegative() == gamma1MCLabel) && isTrueEta)
					fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	
			if ((TrueGammaCandidate1->IsMerged() || TrueGammaCandidate1->IsMergedPartConv() || TrueGammaCandidate1->IsDalitzMerged()) && !matched ){
				if (isTruePi0) fHistoTruePi0CaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta) fHistoTrueEtaCaloMergedClusterInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	
			if (TrueGammaCandidate1->IsMergedPartConv() && !matched) {
				if (isTruePi0) fHistoTruePi0CaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				if (isTrueEta) fHistoTrueEtaCaloMergedClusterPartConvInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}	
		}

		if ( !matched){
			if (fDoMesonQA > 0){
				if (isTruePi0){
					if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
						fHistoTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
                        fHistoTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()));
						fHistoTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
						fHistoTrueMotherPi0ConvPhotonEtaPhi[fiCut]->Fill(TrueGammaCandidate0->GetPhotonPhi(), TrueGammaCandidate0->GetPhotonEta());
					}
				} else if (isTrueEta){
					if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
						fHistoTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift());
                        fHistoTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()));
						fHistoTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle());
						fHistoTrueMotherEtaConvPhotonEtaPhi[fiCut]->Fill(TrueGammaCandidate0->GetPhotonPhi(), TrueGammaCandidate0->GetPhotonEta());
					}
				}
			}
			Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel)), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
			if(isPrimary){ // Secondary Meson
				Int_t secMotherLabel = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->GetMother();
				Float_t weightedSec= 1;
				if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, 0x0, fInputEvent) && static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
					weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, 0x0, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
					//cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
				}
				if (isTruePi0){
					fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
					fHistoTrueSecondaryPi0PhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->Pt(),weightedSec);
					if(CheckIfContainedInString(fStringRecTruePi0s[fiCut],gamma0MotherLabel)){
						fHistoTrueSecondaryPi0DCPtconv[fiCut]->Fill(TrueGammaCandidate0->Pt(),weightedSec);
					}						
				}	
				if (secMotherLabel >-1){
					if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310 && isTruePi0){
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
					weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, 0x0, fInputEvent);
					//                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
					}
				}
				if (isTruePi0){
					fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
					fHistoTruePrimaryPi0W0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					fHistoTruePrimaryPi0PhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->Pt(),weighted);
					fProfileTruePrimaryPi0WeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
					if(CheckIfContainedInString(fStringRecTruePi0s[fiCut],gamma0MotherLabel)){
						fHistoTruePrimaryPi0DCPtconv[fiCut]->Fill(TrueGammaCandidate0->Pt(),weighted);
					}						
				} else if (isTrueEta){
					fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
					fHistoTruePrimaryEtaW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					fHistoTruePrimaryEtaPhotonPairPtconv[fiCut]->Fill(Pi0Candidate->M(),TrueGammaCandidate0->Pt(),weighted);
					fProfileTruePrimaryEtaWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);					
					if(CheckIfContainedInString(fStringRecTrueEtas[fiCut],gamma0MotherLabel)){
						fHistoTruePrimaryEtaDCPtconv[fiCut]->Fill(TrueGammaCandidate0->Pt(),weighted);
					}												
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
	
	if (isTrueEta && !matched){
		fStringRecTrueEtas[fiCut].Append(Form("%i,",gamma0MotherLabel));
	}
	if (isTruePi0 && !matched){
		fStringRecTruePi0s[fiCut].Append(Form("%i,",gamma0MotherLabel));
	}	


	
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::ProcessPhotonCandidates()
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
		  
            Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
			if(isPosFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
            Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
			if(isNegFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
			if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
		}
		
		if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
		if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->InPlaneOutOfPlaneCut(PhotonCandidate->GetPhotonPhi(),fEventPlaneAngle)) continue;
		if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut() &&
		!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){
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
			
		} else if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
			((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->FillElectonLabelArray(PhotonCandidate,nV0);
			nV0++;
			GammaCandidatesStepOne->Add(PhotonCandidate);
		} else if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut() &&
				((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
			GammaCandidatesStepTwo->Add(PhotonCandidate);
		}
	}
	
	if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseElecSharingCut()){
		for(Int_t i = 0;i<GammaCandidatesStepOne->GetEntries();i++){
			AliAODConversionPhoton *PhotonCandidate= (AliAODConversionPhoton*) GammaCandidatesStepOne->At(i);
			if(!PhotonCandidate) continue;
			fIsFromMBHeader = kTRUE;
            if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
                Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
                Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
				if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
			}
			if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GammaCandidatesStepOne->GetEntries())) continue;
			if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
				fGammaCandidates->Add(PhotonCandidate);
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
				
			}
			else GammaCandidatesStepTwo->Add(PhotonCandidate); // Close v0s cut enabled -> add to list two
			
			
		}
	}
	if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->UseToCloseV0sCut()){
		for(Int_t i = 0;i<GammaCandidatesStepTwo->GetEntries();i++){
			AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GammaCandidatesStepTwo->At(i);
			if(!PhotonCandidate) continue;
			fIsFromMBHeader = kTRUE;
            if(fMCEvent && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
                Int_t isPosFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
                Int_t isNegFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent, fInputEvent);
				if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
			}
			if(!((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GammaCandidatesStepTwo,i)) continue;
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
			
		}
	}
	
	delete GammaCandidatesStepOne;
	GammaCandidatesStepOne = 0x0;
	delete GammaCandidatesStepTwo;
	GammaCandidatesStepTwo = 0x0;
  
}
//________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::ProcessTruePhotonCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate)
{
	
	Double_t magField = fInputEvent->GetMagneticField();
	if( magField  < 0.0 ){
		magField =  1.0;
	}
	else {
		magField =  -1.0;
	}
	
	TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	if (AODMCTrackArray == NULL) return;
	AliAODMCParticle *posDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelPositive());
	AliAODMCParticle *negDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelNegative());
	fCharPhotonMCInfo = 0;
	
	if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
    Int_t pdgCode[2] = {TMath::Abs(posDaughter->GetPdgCode()),TMath::Abs(negDaughter->GetPdgCode())};
	
	if(posDaughter->GetMother() != negDaughter->GetMother()){
		
		fCharPhotonMCInfo = 1;
		return;
	}
	else if(posDaughter->GetMother() == -1){
		
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
		if (CheckVectorForDoubleCount(fVectorDoubleCountTrueConvGammas,posDaughter->GetMother())) fHistoDoubleCountTrueConvGammaRPt[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),TruePhotonCandidate->Pt());
	}

	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();
	Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);

	if(isPrimary){
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
	TruePhotonCandidate->SetIsTrueConvertedPhoton();
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::ProcessTruePhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();	
	
	Double_t magField = fInputEvent->GetMagneticField();
	if( magField  < 0.0 ){
		magField =  1.0;
	}
	else {
		magField =  -1.0;
	}

	// Process True Photons
    TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(fMCEvent);
    TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCEvent);
	fCharPhotonMCInfo = 0;
	
	if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
    Int_t pdgCode[2] = {TMath::Abs(posDaughter->GetPdgCode()),TMath::Abs(negDaughter->GetPdgCode())};
	fCharPhotonMCInfo = 1;
	if(posDaughter->GetMother(0) != negDaughter->GetMother(0)){
		
		return;
	}
	else if(posDaughter->GetMother(0) == -1){
		
		return;
	}
	
	if(pdgCode[0]!=11 || pdgCode[1]!=11) return; //One Particle is not a electron
	
	if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge

    TParticle *Photon = TruePhotonCandidate->GetMCParticle(fMCEvent);
	
	if(Photon->GetPdgCode() != 22){
		return; // Mother is no Photon
	}
	
	if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() !=5) return;// check if the daughters come from a conversion
	
	// True Photon
	if(fIsFromMBHeader){
		fHistoTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
		if (fDoPhotonQA > 0) fHistoTrueConvGammaEta[fiCut]->Fill(TruePhotonCandidate->Eta());
		if (CheckVectorForDoubleCount(fVectorDoubleCountTrueConvGammas,posDaughter->GetMother(0))) fHistoDoubleCountTrueConvGammaRPt[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),TruePhotonCandidate->Pt());
	}
    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, posDaughter->GetMother(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
	if(isPrimary){
		// Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
		if(fIsFromMBHeader){
			fCharPhotonMCInfo = 6;
			fHistoTruePrimaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			fHistoTruePrimaryConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt()); // Allways Filled
		}
		// (Not Filled for i6, Extra Signal Gamma (parambox) are secondary)
	} else {
		// fill secondary histograms
		if(fIsFromMBHeader){
			fCharPhotonMCInfo = 2;
			fHistoTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
            if(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0) > -1 &&
                fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 3122){
				fHistoTrueSecondaryConvGammaFromXFromLambdaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				fCharPhotonMCInfo = 5;
			}
            if(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0) > -1 &&
                fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
				fHistoTrueSecondaryConvGammaFromXFromK0sPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				fCharPhotonMCInfo = 4;
			}
            if(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0) > -1 &&
                fMCEvent->Particle(fMCEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 221){
				fCharPhotonMCInfo = 3;
			}
		}
	}
	TruePhotonCandidate->SetIsTrueConvertedPhoton();
	return;
}
//________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::CalculateBackground(){
	
	Int_t zbin= fBGClusHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
	Int_t mbin = 0;
	
	
	if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
		mbin = fBGClusHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
	} else {
		mbin = fBGClusHandler[fiCut]->GetMultiplicityBinIndex(fClusterCandidates->GetEntries());
	}
	
	
	AliGammaConversionAODBGHandler::GammaConversionVertex *bgEventVertex = NULL;    
	
	if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
		
		for(Int_t nEventsInBG=0;nEventsInBG<fBGClusHandler[fiCut]->GetNBGEvents();nEventsInBG++){
		  
			AliGammaConversionAODVector *previousEventV0s = fBGClusHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
			if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
				bgEventVertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
			}
			
			for(Int_t iCurrent=0;iCurrent<fVirtualGammaCandidates->GetEntries();iCurrent++){
			  
				AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fVirtualGammaCandidates->At(iCurrent));
				
				for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
					AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
					
					if(fMoveParticleAccordingToVertex == kTRUE){
						if (bgEventVertex){
							MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
						}	
					}
					if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
						if (bgEventVertex){
							RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
						}	
					}
					
					AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
					backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
					if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
						->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
					  
					  
					       if(  ((AliDalitzElectronCuts*)fElectronCutArray->At(fiCut))->DoMassCut() == kTRUE ) {
						 
							if( ((AliDalitzElectronCuts*) fElectronCutArray->At(fiCut))->MassCut( backgroundCandidate->Pt() , currentEventGoodV0.M() ) == kTRUE ){
					  
					  					  
							    fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
							    
							    if(fDoTHnSparse){
								Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
								fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
							    }
							}
					       }
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
				if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0 ){
					bgEventVertex = fBGClusHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
				}
				for(Int_t iCurrent=0;iCurrent<fVirtualGammaCandidates->GetEntries();iCurrent++){
					AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fVirtualGammaCandidates->At(iCurrent));
					for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
				
						AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
					
						if(fMoveParticleAccordingToVertex == kTRUE){
							if (bgEventVertex){
								MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
							}	
						}
						if(((AliConversionPhotonCuts*)fGammaCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
							if (bgEventVertex){
								RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
							}	
						}
					
						AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
						backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
						if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
						  
						   if(  ((AliDalitzElectronCuts*)fElectronCutArray->At(fiCut))->DoMassCut() == kTRUE ) {
						 
							if( ((AliDalitzElectronCuts*) fElectronCutArray->At(fiCut))->MassCut( backgroundCandidate->Pt() , currentEventGoodV0.M() ) == kTRUE ){
					  
						  
							      fHistoMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
							      							
							      if(fDoTHnSparse){
								Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
								fSparseMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
							      }
							}
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
void AliAnalysisTaskGammaCaloDalitzV1::RotateParticle(AliAODConversionPhoton *gamma){
	Int_t fNDegreesPMBackground= ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->NDegreesRotation();
	Double_t nRadiansPM = fNDegreesPMBackground*TMath::Pi()/180;
	Double_t rotationValue = fRandom.Rndm()*2*nRadiansPM + TMath::Pi()-nRadiansPM;
	gamma->RotateZ(rotationValue);
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::RotateParticleAccordingToEP(AliAODConversionPhoton *gamma, Double_t previousEventEP, Double_t thisEventEP){
	
	previousEventEP=previousEventEP+TMath::Pi();
	thisEventEP=thisEventEP+TMath::Pi();
	Double_t rotationValue= thisEventEP-previousEventEP;
	gamma->RotateZ(rotationValue);
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex){
	//see header file for documentation
	
	Double_t dx = vertex->fX - fInputEvent->GetPrimaryVertex()->GetX();
	Double_t dy = vertex->fY - fInputEvent->GetPrimaryVertex()->GetY();
	Double_t dz = vertex->fZ - fInputEvent->GetPrimaryVertex()->GetZ();
	
	Double_t movedPlace[3] = {particle->GetConversionX() - dx,particle->GetConversionY() - dy,particle->GetConversionZ() - dz};
	particle->SetConversionPoint(movedPlace);
}
//________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::UpdateEventByEventData(){
	//see header file for documentation
	if( fClusterCandidates->GetEntries() > 0  ){
		if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
			fBGClusHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
		} else { // means we use #V0s for multiplicity
			fBGClusHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fClusterCandidates->GetEntries(),fEventPlaneAngle);
		}
	}
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::RelabelAODPhotonCandidates(Bool_t mode){
	
	// Relabeling For AOD Event
	// ESDiD -> AODiD
	// MCLabel -> AODMCLabel
	
	if(mode){
        fMCEventPos = new Int_t[fReaderGammas->GetEntries()];
        fMCEventNeg = new Int_t[fReaderGammas->GetEntries()];
		fESDArrayPos = new Int_t[fReaderGammas->GetEntries()];
		fESDArrayNeg = new Int_t[fReaderGammas->GetEntries()];
	}
	
	for(Int_t iGamma = 0;iGamma<fReaderGammas->GetEntries();iGamma++){
		AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(iGamma);
		if(!PhotonCandidate) continue;
		if(!mode){// Back to ESD Labels
        PhotonCandidate->SetMCLabelPositive(fMCEventPos[iGamma]);
        PhotonCandidate->SetMCLabelNegative(fMCEventNeg[iGamma]);
		PhotonCandidate->SetLabelPositive(fESDArrayPos[iGamma]);
		PhotonCandidate->SetLabelNegative(fESDArrayNeg[iGamma]);
		continue;
		}
        fMCEventPos[iGamma] =  PhotonCandidate->GetMCLabelPositive();
        fMCEventNeg[iGamma] =  PhotonCandidate->GetMCLabelNegative();
		fESDArrayPos[iGamma] = PhotonCandidate->GetTrackLabelPositive();
		fESDArrayNeg[iGamma] = PhotonCandidate->GetTrackLabelNegative();
		
		Bool_t AODLabelPos = kFALSE;
		Bool_t AODLabelNeg = kFALSE;
		
		for(Int_t i = 0; i<fInputEvent->GetNumberOfTracks();i++){
		AliAODTrack *tempDaughter = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
		if(!AODLabelPos){
			if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelPositive() ){
            PhotonCandidate->SetMCLabelPositive(TMath::Abs(tempDaughter->GetLabel()));
			PhotonCandidate->SetLabelPositive(i);
			AODLabelPos = kTRUE;
			}
		}
		if(!AODLabelNeg){
			if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelNegative()){
            PhotonCandidate->SetMCLabelNegative(TMath::Abs(tempDaughter->GetLabel()));
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
        delete[] fMCEventPos;
        delete[] fMCEventNeg;
		delete[] fESDArrayPos;
		delete[] fESDArrayNeg;
	}
}

void AliAnalysisTaskGammaCaloDalitzV1::SetLogBinningXTH2(TH2* histoRebin){
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
void AliAnalysisTaskGammaCaloDalitzV1::Terminate(const Option_t *)
{
  
  //fOutputContainer->Print(); // Will crash on GRID
}

//________________________________________________________________________
Int_t AliAnalysisTaskGammaCaloDalitzV1::GetSourceClassification(Int_t daughter, Int_t pdgCode){
  
	if (daughter == 111) {
        if (TMath::Abs(pdgCode) == 310) return 1; // k0s
        else if (TMath::Abs(pdgCode) == 3122) return 2; // Lambda
        else if (TMath::Abs(pdgCode) == 130) return 3; // K0L
        else if (TMath::Abs(pdgCode) == 2212) return 4; // proton
        else if (TMath::Abs(pdgCode) == 2112) return 5; // neutron
        else if (TMath::Abs(pdgCode) == 211) return 6; // pion
        else if (TMath::Abs(pdgCode) == 321) return 7; // kaon
        else if (TMath::Abs(pdgCode) == 113 || TMath::Abs(pdgCode) == 213 ) return 8; // rho 0,+,-
        else if (TMath::Abs(pdgCode) == 3222 || TMath::Abs(pdgCode) == 3212 || TMath::Abs(pdgCode) == 3112  ) return 9; // Sigma
        else if (TMath::Abs(pdgCode) == 2224 || TMath::Abs(pdgCode) == 2214 || TMath::Abs(pdgCode) == 2114 || TMath::Abs(pdgCode) == 1114  ) return 10; // Delta
        else if (TMath::Abs(pdgCode) == 313 || TMath::Abs(pdgCode) == 323   ) return 11; // K*
		else return 15;
	}
	return 15;
  
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskGammaCaloDalitzV1::CheckIfContainedInString(TString input, Int_t tobechecked){
	TObjArray *arr = input.Tokenize(",");
	for (Int_t i = 0; i < arr->GetEntriesFast();i++){
		TString tempStr = ((TObjString*)arr->At(i))->GetString();
		if (tempStr.Atoi() == tobechecked) return kTRUE;
	}	
	return kFALSE;
}

//_________________________________________________________________________________
Bool_t AliAnalysisTaskGammaCaloDalitzV1::CheckIfContainedInStringAndAppend(TString &input, Int_t tobechecked){
	TObjArray *arr = input.Tokenize(",");
	Bool_t isContained = kFALSE;
	for (Int_t i = 0; i < arr->GetEntriesFast();i++){
		TString tempStr = ((TObjString*)arr->At(i))->GetString();
		if (tempStr.Atoi() == tobechecked) isContained= kTRUE;
	}	
	if (!isContained)input.Append(Form("%i,",tobechecked));	
	return isContained;
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::ProcessConversionPhotonsForMissingTags (){

    if (!fMCEvent) return;
	
	const AliVVertex* 	primVtxMC = fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();

	for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
	  
		AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
		
		if (gamma0->IsTrueConvertedPhoton()){
			Int_t gamma0MotherLabel = -1;
            Int_t gamma0MCLabel = gamma0->GetMCParticleLabel(fMCEvent);
			if(gamma0MCLabel != -1){ 
                TParticle * gammaMC0 = (TParticle*)fMCEvent->Particle(gamma0MCLabel);
				gamma0MotherLabel=gammaMC0->GetFirstMother();
				if (gamma0MotherLabel>-1){
                    if(((TParticle*)fMCEvent->Particle(gamma0MotherLabel))->GetPdgCode() == 111){
						if (!CheckIfContainedInStringAndAppend(fStringRecTruePi0s[fiCut],gamma0MotherLabel)){ 
						  
                                Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, gamma0MotherLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
													      
							if (!isPrimary){
                                Int_t secMotherLabel = ((TParticle*)fMCEvent->Particle(gamma0MotherLabel))->GetMother(0);
								Float_t weightedSec= 1;
                                if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, fMCEvent, fInputEvent) && fMCEvent->Particle(secMotherLabel)->GetPdgCode()==310){
                                    weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, fMCEvent, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
								}
								fHistoTrueSecondaryPi0MissingPtconv[fiCut]->Fill(gamma0->Pt(),weightedSec);
							} else {
								Float_t weighted= 1;
                                if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma0MotherLabel, fMCEvent, fInputEvent)){
                                    if (((TParticle*)fMCEvent->Particle(gamma0MotherLabel))->Pt()>0.005){
                                        weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma0MotherLabel, fMCEvent, fInputEvent);
									}
								}
								fHistoTruePrimaryPi0MissingPtconv[fiCut]->Fill(gamma0->Pt(),weighted);		
							}	
						}
                    } else if(((TParticle*)fMCEvent->Particle(gamma0MotherLabel))->GetPdgCode() == 221){
						if (!CheckIfContainedInStringAndAppend(fStringRecTrueEtas[fiCut],gamma0MotherLabel)){
							Float_t weighted= 1;
                            if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma0MotherLabel, fMCEvent, fInputEvent)){
                                if (((TParticle*)fMCEvent->Particle(gamma0MotherLabel))->Pt()>0.005){
                                    weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma0MotherLabel, fMCEvent, fInputEvent);
								}
							}
							fHistoTruePrimaryEtaMissingPtconv[fiCut]->Fill(gamma0->Pt(),weighted);
						}	
					}	
				}
			} 	
		}	
	}	
	
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::ProcessConversionPhotonsForMissingTagsAOD (){

	const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
	Double_t mcProdVtxX 	= primVtxMC->GetX();
	Double_t mcProdVtxY 	= primVtxMC->GetY();
	Double_t mcProdVtxZ 	= primVtxMC->GetZ();

	TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
		
	for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){
		AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));

		if (gamma0->IsTrueConvertedPhoton()){
			AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0->GetMCLabelPositive()));
			AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0->GetMCLabelNegative()));
			
			Int_t gamma0MCLabel = -1;
			Int_t gamma0MotherLabel = -1;
			if(!positiveMC||!negativeMC)
				return;
			
			if (gamma0->IsTrueConvertedPhoton()){
				gamma0MCLabel = positiveMC->GetMother();
				AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
				gamma0MotherLabel=gammaMC0->GetMother();

				if (gamma0MotherLabel>-1){
					if(((AliAODMCParticle*)AODMCTrackArray->At(gamma0MotherLabel))->GetPdgCode() == 111){
						if (!CheckIfContainedInStringAndAppend(fStringRecTruePi0s[fiCut],gamma0MotherLabel)){ 
							Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel)), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
							if (isPrimary){
								Int_t secMotherLabel = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel))->GetMother();
								Float_t weightedSec= 1;
								if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, 0x0, fInputEvent) && static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
									weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, 0x0, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
								}
								fHistoTrueSecondaryPi0MissingPtconv[fiCut]->Fill(gamma0->Pt(),weightedSec);
							} else {
								Float_t weighted= 1;
								if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma0MotherLabel, 0x0, fInputEvent)){
									if (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel))->Pt()>0.005){
										weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma0MotherLabel, 0x0, fInputEvent);
									}
								}
								fHistoTruePrimaryPi0MissingPtconv[fiCut]->Fill(gamma0->Pt(),weighted);		
							}	
						}
					} else if(((AliAODMCParticle*)AODMCTrackArray->At(gamma0MotherLabel))->GetPdgCode() == 221){
						if (!CheckIfContainedInStringAndAppend(fStringRecTrueEtas[fiCut],gamma0MotherLabel)){
							Float_t weighted= 1;
							if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma0MotherLabel, 0x0, fInputEvent)){
								if (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel))->Pt()>0.005){
									weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma0MotherLabel, 0x0, fInputEvent);
								}
							}
							fHistoTruePrimaryEtaMissingPtconv[fiCut]->Fill(gamma0->Pt(),weighted);
						}	
					}	
				}
			} 	
		}	
	}	
	
	return;
}	

//________________________________________________________________________
void AliAnalysisTaskGammaCaloDalitzV1::ProcessElectronCandidates(){

	 const AliVVertex* primVtxMC 	= 0;
	 Double_t mcProdVtxX 	= 0;
	 Double_t mcProdVtxY 	= 0;
	 Double_t mcProdVtxZ 	= 0;	
	
  
  
  
	if ( fMCEvent ) {
  
	  primVtxMC 	= fMCEvent->GetPrimaryVertex();
	  mcProdVtxX 	= primVtxMC->GetX();
	  mcProdVtxY 	= primVtxMC->GetY();
	  mcProdVtxZ 	= primVtxMC->GetZ();	
	
	}
  
  
   
	Double_t magField = fInputEvent->GetMagneticField();
	

	if( magField  < 0.0 ){
	  
		magField =   1.0;
		
	} else {
	  
		magField =  -1.0;
	}
	
		 
	vector<Int_t> lGoodElectronIndexPrev(0);
	vector<Int_t> lGoodPositronIndexPrev(0);
	
	
	
	for(UInt_t i = 0; i < fSelectorElectronIndex.size(); i++){
		
		AliESDtrack* electronCandidate = (AliESDtrack*)fInputEvent->GetTrack(fSelectorElectronIndex[i]);
		
		if( fMCEvent ) {
		  
			Int_t labelelectron = TMath::Abs( electronCandidate->GetLabel() );
			Int_t isMCFromMBHeader = -1;
			
			if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0) {
                isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(labelelectron,fMCEvent,fInputEvent);
				if( isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut) )->GetSignalRejection() != 3) continue;
			}    
			
		}
			
			
		if(! ((AliDalitzElectronCuts*)fElectronCutArray->At(fiCut))->ElectronIsSelected(electronCandidate) ) continue;
		
		lGoodElectronIndexPrev.push_back(   fSelectorElectronIndex[i] );
		fHistoDalitzElectronPt[fiCut]->Fill(electronCandidate->Pt());
		fHistoDalitzElectronPhi[fiCut]->Fill(electronCandidate->Phi());
		
		
		if( fMCEvent ) {
			
			Int_t labelelectron = TMath::Abs( electronCandidate->GetLabel() );
			
            Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelelectron, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
				
			
            if( labelelectron >=0 && labelelectron < fMCEvent->GetNumberOfTracks() ){
                TParticle* electron = fMCEvent->Particle(labelelectron);
				if( electron->GetPdgCode() ==  11 ){
					if( isPrimary ){
						fHistoTrueElectronPt[fiCut]->Fill(electronCandidate->Pt());    //primary electron
					} else {
						fHistoTrueSecElectronPt[fiCut]->Fill(electronCandidate->Pt()); //secondary electron
					}
					if( IsPi0DalitzDaughter(labelelectron) == kTRUE ) {
						if( isPrimary ) {
							fHistoTruePi0DalitzElectronPt[fiCut]->Fill(electronCandidate->Pt());
						} else {
							fHistoTruePi0DalitzSecElectronPt[fiCut]->Fill(electronCandidate->Pt());
						}
					}
				}
			}
		}
		
	}
	
	

	for(UInt_t i = 0; i < fSelectorPositronIndex.size(); i++){
	  

		AliESDtrack* positronCandidate = (AliESDtrack*)fInputEvent->GetTrack( fSelectorPositronIndex[i] );
		
		if( fMCEvent ) {
		
			Int_t labelpositron = TMath::Abs( positronCandidate->GetLabel() );
			Int_t isMCFromMBHeader = -1;
			
			if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0) {
                isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(labelpositron,fMCEvent,fInputEvent);
				if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
			}  
		}
		
	
		if(! ((AliDalitzElectronCuts*)fElectronCutArray->At(fiCut))->ElectronIsSelected(positronCandidate) ) continue;
		
		lGoodPositronIndexPrev.push_back( fSelectorPositronIndex[i] );
		fHistoDalitzPositronPt[fiCut]->Fill( positronCandidate->Pt() );
		fHistoDalitzPositronPhi[fiCut]->Fill( positronCandidate->Phi() );
			
		if( fMCEvent ) {
			Int_t labelpositron = TMath::Abs( positronCandidate->GetLabel() );
            Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, labelpositron, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
			
			
            if( labelpositron >=0 && labelpositron < fMCEvent->GetNumberOfTracks() ) {
                TParticle* positron = fMCEvent->Particle(labelpositron);
				if( positron->GetPdgCode() ==  -11 ){
					if( isPrimary ){
						fHistoTruePositronPt[fiCut]->Fill(positronCandidate->Pt());
					} else {
						fHistoTrueSecPositronPt[fiCut]->Fill(positronCandidate->Pt());
					}
					if( IsPi0DalitzDaughter(labelpositron) == kTRUE ) {
						if( isPrimary ){
							fHistoTruePi0DalitzPositronPt[fiCut]->Fill(positronCandidate->Pt());
						} else {
							fHistoTruePi0DalitzSecPositronPt[fiCut]->Fill(positronCandidate->Pt());
						}
					}
				}
			}
		}
	}

	vector<Bool_t> lElectronPsiIndex(lGoodElectronIndexPrev.size(), kTRUE);
	vector<Bool_t> lPositronPsiIndex(lGoodPositronIndexPrev.size(), kTRUE);

	if( ((AliDalitzElectronCuts*)fElectronCutArray->At(fiCut))->DoPsiPairCut() == kTRUE ){
		for( UInt_t i = 0; i < lGoodElectronIndexPrev.size(); i++ ) {
			AliESDtrack *electronCandidate = (AliESDtrack*) fInputEvent->GetTrack(lGoodElectronIndexPrev[i]);
			for(UInt_t j = 0; j <  lGoodPositronIndexPrev.size(); j++){
			  
				AliESDtrack *positronCandidate = (AliESDtrack*) fInputEvent->GetTrack(lGoodPositronIndexPrev[j]);
				Double_t psiPair = GetPsiPair(positronCandidate,electronCandidate);
				Double_t deltaPhi = magField * TVector2::Phi_mpi_pi( electronCandidate->GetConstrainedParam()->Phi()-positronCandidate->GetConstrainedParam()->Phi());

				if( ((AliDalitzElectronCuts*)fElectronCutArray->At(fiCut))->IsFromGammaConversion(psiPair,deltaPhi ) ){
					lElectronPsiIndex[i] = kFALSE;
					lPositronPsiIndex[j] = kFALSE;
				}
			}
		}
	}

	vector<Int_t> lGoodElectronIndex(0);
	vector<Int_t> lGoodPositronIndex(0);
	
	for( UInt_t i = 0; i < lGoodElectronIndexPrev.size(); i++ ) {
		if(  lElectronPsiIndex[i] == kTRUE )
		lGoodElectronIndex.push_back(   lGoodElectronIndexPrev[i]  );
	}

	for( UInt_t i = 0; i < lGoodPositronIndexPrev.size(); i++ ) {
		if(  lPositronPsiIndex[i] == kTRUE )
		lGoodPositronIndex.push_back(   lGoodPositronIndexPrev[i]  );
	}
	

	for(UInt_t i = 0; i < lGoodElectronIndex.size(); i++){
	  
		AliESDtrack *electronCandidate = (AliESDtrack*)fInputEvent->GetTrack(lGoodElectronIndex[i]);
		AliKFParticle electronCandidateKF( *electronCandidate->GetConstrainedParam(), ::kElectron );
		
		TLorentzVector electronCandidateTLV;
		
		electronCandidateTLV.SetXYZM(electronCandidate->GetConstrainedParam()->Px(),electronCandidate->GetConstrainedParam()->Py(),electronCandidate->GetConstrainedParam()->Pz(),TDatabasePDG::Instance()->GetParticle(  ::kElectron   )->Mass());
		

		for(UInt_t j = 0; j < lGoodPositronIndex.size(); j++){

			AliESDtrack *positronCandidate = (AliESDtrack*) fInputEvent->GetTrack(lGoodPositronIndex[j]);
			AliKFParticle positronCandidateKF( *positronCandidate->GetConstrainedParam(), ::kPositron );
			TLorentzVector positronCandidateTLV;
		
			positronCandidateTLV.SetXYZM(positronCandidate->GetConstrainedParam()->Px(),positronCandidate->GetConstrainedParam()->Py(),positronCandidate->GetConstrainedParam()->Pz(),TDatabasePDG::Instance()->GetParticle(  ::kPositron   )->Mass());
			
			TLorentzVector *virtualPhotonTLV = 0;
			AliKFConversionPhoton* virtualPhoton = NULL;
			AliAODConversionPhoton *vphoton;
			
			if( ((AliDalitzElectronCuts*)fElectronCutArray->At(fiCut))->GetUseElectronMCSmearing() && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseMCPSmearing() && fMCEvent){
				TLorentzVector smearelectronCandidateTLV = ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->SmearElectron(electronCandidateTLV);
				TLorentzVector smearpositronCandidateTLV = ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->SmearElectron(positronCandidateTLV);
				virtualPhotonTLV = new TLorentzVector( smearelectronCandidateTLV + smearpositronCandidateTLV );
				vphoton= new AliAODConversionPhoton(virtualPhotonTLV); 
				vphoton->SetMass(virtualPhotonTLV->M());
			} else {
				virtualPhoton = new AliKFConversionPhoton(electronCandidateKF,positronCandidateKF);
				
					AliKFVertex primaryVertexImproved(*fInputEvent->GetPrimaryVertex());
					primaryVertexImproved+=*virtualPhoton;
					virtualPhoton->SetProductionVertex(primaryVertexImproved);			
				
				        vphoton= new AliAODConversionPhoton(virtualPhoton);
			}
			vphoton->SetTrackLabels( lGoodPositronIndex[j], lGoodElectronIndex[i]);
			
			if( fMCEvent ) {
				Int_t labeln=TMath::Abs(electronCandidate->GetLabel());
				Int_t labelp=TMath::Abs(positronCandidate->GetLabel());
                TParticle *fNegativeMCParticle = fMCEvent->Particle(labeln);
                TParticle *fPositiveMCParticle = fMCEvent->Particle(labelp);
				
				if( fPositiveMCParticle && fNegativeMCParticle) {		
					vphoton->SetMCLabelPositive(labelp);
					vphoton->SetMCLabelNegative(labeln);
						
				}		
			}
			
					
			fVirtualGammaCandidates->Add(  vphoton );
			
			if( virtualPhoton ){
				delete virtualPhoton;
				virtualPhoton=NULL;
			}
			if ( virtualPhotonTLV ){
				delete virtualPhotonTLV;
				virtualPhotonTLV=NULL;
			}
					
		}
	}
	
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskGammaCaloDalitzV1::IsPi0DalitzDaughter( Int_t label ) const
{
//
// Returns true if the particle comes from Pi0 -> e+ e- gamma
//        
    Int_t motherLabel = fMCEvent->Particle( label )->GetMother(0);
	
    if( motherLabel < 0 || motherLabel >= fMCEvent->GetNumberOfTracks() ) return kFALSE;
	
    TParticle* mother = fMCEvent->Particle( motherLabel );

	if( mother->GetPdgCode() != 111 ) return kFALSE;
	
	if( IsDalitz( mother ) ) return kTRUE;
		
	return kFALSE;
   
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskGammaCaloDalitzV1::GetPsiPair( const AliESDtrack *trackPos, const AliESDtrack *trackNeg ) const
{
	//
	// This angle is a measure for the contribution of the opening in polar
	// direction ?0 to the opening angle ? Pair
	//
	// Ref. Measurement of photons via conversion pairs with the PHENIX experiment at RHIC
	//      Master Thesis. Thorsten Dahms. 2005
	// https://twiki.cern.ch/twiki/pub/ALICE/GammaPhysicsPublications/tdahms_thesis.pdf
	//
	Double_t momPos[3];
	Double_t momNeg[3];
	if( trackPos->GetConstrainedPxPyPz(momPos) == 0 ) trackPos->GetPxPyPz( momPos );
	if( trackNeg->GetConstrainedPxPyPz(momNeg) == 0 ) trackNeg->GetPxPyPz( momNeg );

	TVector3 posDaughter;
	TVector3 negDaughter;

	posDaughter.SetXYZ( momPos[0], momPos[1], momPos[2] );
	negDaughter.SetXYZ( momNeg[0], momNeg[1], momNeg[2] );

	Double_t deltaTheta = negDaughter.Theta() - posDaughter.Theta();
	Double_t openingAngle =  posDaughter.Angle( negDaughter );  //TMath::ACos( posDaughter.Dot(negDaughter)/(negDaughter.Mag()*posDaughter.Mag()) );

	if( openingAngle < 1e-20 ) return 0.;

	Double_t psiAngle = TMath::ASin( deltaTheta/openingAngle );

	return psiAngle;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaCaloDalitzV1::IsDalitz(TParticle *fMCMother) const
{

	if( fMCMother->GetNDaughters() != 3 ) return kFALSE;
	if( fMCMother->GetPdgCode() != 111 && fMCMother->GetPdgCode() != 221 ) return kFALSE;
	
	
	TParticle *positron = 0x0;
	TParticle *electron = 0x0;
	TParticle *gamma    = 0x0;
	
	for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){				
        TParticle* temp = (TParticle*)fMCEvent->Particle( index );
		switch( temp->GetPdgCode() ) {
			case ::kPositron:
				positron =  temp;
				break;
			case ::kElectron:
				electron =  temp;
				break;
			case ::kGamma:
				gamma    =  temp;
				break;
		}
	}
  
	if( positron && electron && gamma) return kTRUE;
	
	return kFALSE;
}

//_________________________________________________________________
Int_t AliAnalysisTaskGammaCaloDalitzV1::FindMotherOfPhoton(Int_t particleLabel ){
  
    
  
	    if( particleLabel < 0 ) return -1;
	    
        TParticle* particle = (TParticle*)fMCEvent->Particle( particleLabel );
	    
	    if( particle->GetMother(0) < 0 ) return -1;
	    
        else if ( ((TParticle*)fMCEvent->Particle( particle->GetMother(0) ) )->GetPdgCode() == 111 ||
              ((TParticle*)fMCEvent->Particle( particle->GetMother(0) ) )->GetPdgCode() == 221 )  {
			
			return particle->GetMother(0);
		      }
	     
	    else return FindMotherOfPhoton( particle->GetMother(0) );
	   
         
}
  
//_________________________________________________________________________________
Bool_t AliAnalysisTaskGammaCaloDalitzV1::CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked)
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

