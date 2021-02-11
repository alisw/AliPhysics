/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                               *
 * Author: Friederike Bock                                                         *
 * Version 1.0                                                                           *
 *                                                                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its    *
 * documentation strictly for non-commercial purposes is hereby granted    *
 * without fee, provided that the above copyright notice appears in all    *
 * copies and that both the copyright notice and this permission notice    *
 * appear in the supporting documentation. The authors make no claims    *
 * about the suitability of this software for any purpose. It is      *
 * provided "as is" without express or implied warranty.                           *
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
#include "AliCaloTrackMatcher.h"
#include <vector>
#include <map>


ClassImp(AliAnalysisTaskGammaCaloMerged)

//________________________________________________________________________
AliAnalysisTaskGammaCaloMerged::AliAnalysisTaskGammaCaloMerged(): AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fCorrTaskSetting(""),
  fDoLightOutput(kFALSE),
  fInputEvent(NULL),
  fMCEvent(NULL),
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
  fOutlierJetReader(NULL),
  fConvJetReader(NULL),
  fDoJetAnalysis(kFALSE),
  fDoJetQA(kFALSE),
  fDoOutOfJet(0),
  fJetHistograms(NULL),
  fAODMCTrackArray(NULL),
  farrClustersProcess(NULL),
  fMapNeutralPionOverlap(NULL),
  fMesonOLFromCluster(0),
  fHistoMotherInvMassPt(NULL),
  fHistoMotherPtY(NULL),
  fHistoMotherPtAlpha(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusGammaE(NULL),
  fHistoClusOverlapHeadersGammaPt(NULL),
  fHistoClusNLMPt(NULL),
  fHistoClusMergedPtvsM02(NULL),
  fHistoClusMergedEvsM02(NULL),
  fHistoClusMergedPtvsM02Accepted(NULL),
  fHistoClusMergedEvsM02Accepted(NULL),
  fHistoClusMergedEvsM20(NULL),
  fHistoClusNCellsPt(NULL),
  fHistoClusMergedNCellsPt(NULL),
  fHistoClusMergedNParticlePt(NULL),
  fHistoClusMergedNCellsAroundPt(NULL),
  fHistoClusMergedNCellsAroundAndInPt(NULL),
  fHistoClusMergedEAroundE(NULL),
  fHistoPtJet(NULL),
  fHistoJetEta(NULL),
  fHistoJetPhi(NULL),
  fHistoJetArea(NULL),
  fClusterEtaPhiJets(NULL),
  fHistoNJets(NULL),
  fHistoEventwJets(NULL),
  fHistoTruevsRecJetPt(NULL),
  fHistoClusMergedPtvsRJetAccepted(NULL),
  fHistoJetFragmFunc(NULL),
  fHistoMCHeaders(NULL),
  fHistoMCPi0Pt(NULL),
  fHistoMCPi0ReducedPt(NULL),
  fHistoMCPi0WOWeightPt(NULL),
  fHistoMCPi0WOEvtWeightPt(NULL),
  fHistoMCEtaPt(NULL),
  fHistoMCEtaWOWeightPt(NULL),
  fHistoMCEtaWOEvtWeightPt(NULL),
  fHistoMCPi0DalitzPt(NULL),
  fHistoMCPi0DalitzWOWeightPt(NULL),
  fHistoMCPi0DalitzWOEvtWeightPt(NULL),
  fHistoMCEtaDalitzPt(NULL),
  fHistoMCEtaDalitzWOWeightPt(NULL),
  fHistoMCEtaDalitzWOEvtWeightPt(NULL),
  fHistoMCPi0InAccPt(NULL),
  fHistoMCPi0ReducedInAccPt(NULL),
  fHistoMCEtaInAccPt(NULL),
  fHistoMCPi0WOEvtWeightInAccPt(NULL),
  fHistoMCEtaWOEvtWeightInAccPt(NULL),
  fHistoMCPi0DalitzInAccPt(NULL),
  fHistoMCEtaDalitzInAccPt(NULL),
  fHistoMCPi0DalitzWOEvtWeightInAccPt(NULL),
  fHistoMCEtaDalitzWOEvtWeightInAccPt(NULL),
  fHistoMCSecPi0PtvsSource(NULL),
  fHistoMCSecPi0InAccPtvsSource(NULL),
  fHistoMCPi0PtJetPt(NULL),
  fHistoMCEtaPtJetPt(NULL),
  fHistoMCPrimaryPtvsSource(NULL),
  fHistoMCPrimaryYvsSource(NULL),
  fHistoMCDecayGammaPt(NULL),
  fHistoMCAllGammaPt(NULL),
  fHistoTrueClusEFracFirstLabel(NULL),
  fHistoTrueClusEFracLeadingPi0(NULL),
  fHistoTrueClusNeutralContamination(NULL),
  fHistoTrueClusPhotonContamination(NULL),
  fHistoTrueClusElectronContamination(NULL),
  fHistoTrueClusOtherContamination(NULL),
  fHistoTrueClusMergedPtvsM02(NULL),
  fHistoTrueClusPi0PtvsM02(NULL),
  fHistoTrueClusMultiplePi0PtvsM02(NULL),
  fHistoTrueClusPi0DalitzPtvsM02(NULL),
  fHistoTrueClusPrimPi0PtvsM02(NULL),
  fHistoTrueClusPrimPi0PtMCPt(NULL),
  fHistoTruePureMergedClusPrimPi0PtvsM02(NULL),
  fHistoTruePartConvClusPrimPi0PtvsM02(NULL),
  fHistoTrueClusMergedPartConvPi0EVsM20(NULL),
  fHistoTrueClusPrimPi0PureMergedPtMCPt(NULL),
  fHistoTrueClusPrimPi0PartConvMergedPtMCPt(NULL),
  fHistoTrueClusPrimPi01GammaMergedPtMCPt(NULL),
  fHistoTrueClusPrimPi01ElectronMergedPtMCPt(NULL),
  fHistoTrueClusMultiplePrimPi0PtvsM02(NULL),
  fHistoTrueClusSecPi0PtvsM02(NULL),
  fHistoTrueClusSecPi0FromK0sPtvsM02(NULL),
  fHistoTrueClusSecPi0FromK0lPtvsM02(NULL),
  fHistoTrueClusSecPi0FromLambdaPtvsM02(NULL),
  fHistoTrueClusEtaPtvsM02(NULL),
  fHistoTrueClusEtaDalitzPtvsM02(NULL),
  fHistoTrueClusMergedPureFromPi0PtvsM02(NULL),
  fHistoTrueClusMergedPureFromEtaPtvsM02(NULL),
  fHistoTrueClusMergedPartConvFromPi0PtvsM02(NULL),
  fHistoTrueClusMergedPartConvFromEtaPtvsM02(NULL),
  fHistoNTrueMultiplePi0vsPt(NULL),
  fHistoTrueClusGammaFromPi0PtvsM02(NULL),
  fHistoTrueClusGammaFromEtaPtvsM02(NULL),
  fHistoTrueClusElectronFromPi0PtvsM02(NULL),
  fHistoTrueClusElectronFromEtaPtvsM02(NULL),
  fHistoTrueSecPi0PtvsDiffReco(NULL),
  fHistoTrueClusBGPtvsM02(NULL),
  fHistoTrueClusGammaPtvsM02(NULL),
  fHistoTrueClusGammaPartConvPtvsM02(NULL),
  fHistoTrueClusElectronPtvsM02(NULL),
  fHistoTrueClusElectronFromGammaPtvsM02(NULL),
  fHistoTrueClusMergedInvMassvsPt(NULL),
  fHistoTrueClusPi0InvMassvsPt(NULL),
  fHistoTrueClusPrimPi0InvMassvsPt(NULL),
  fHistoTrueClusEtaInvMassvsPt(NULL),
  fHistoTrueClusBGInvMassvsPt(NULL),
  fHistoTrueClusGammaInvMassvsPt(NULL),
  fHistoTrueClusElectronInvMassvsPt(NULL),
  fHistoTrueClusBGPtvsSource(NULL),
  fHistoTrueClusGammaPtvsSource(NULL),
  fHistoTrueClusElectronPtvsSource(NULL),
//  fHistoTrueClusElectronPtvsMotherID(NULL),
//  fHistoTrueClusElectronPtvsTopMotherID(NULL),
//  fHistoTrueClusElectronPtvsConvPhotonTopMotherID(NULL),
  fHistoTrueMergedMissedPDG(NULL),
  fHistoTrueClusMergedPtvsRJet(NULL),
  fHistoTrueClusPi0PtvsRJet(NULL),
  fHistoTrueClusEtaPtvsRJet(NULL),
  fHistoTruePi0PtY(NULL),
  fHistoTrueEtaPtY(NULL),
  fHistoTruePi0PtAlpha(NULL),
  fHistoTrueEtaPtAlpha(NULL),
  fHistoTrueClusGammaEM02(NULL),
  fHistoTrueClusElectronEM02(NULL),
  fHistoTrueClusPi0EM02(NULL),
  fHistoTrueClusEtaEM02(NULL),
  fHistoTrueClusBGEvsM02(NULL),
  fHistoTrueClusGammaEvsM20(NULL),
  fHistoTrueClusElectronEM20(NULL),
  fHistoTrueClusMergedPi0EVsM20(NULL),
  fHistoTrueClusMergedPi0EVsM02(NULL),
  fHistoTrueClusMergedEtaEVsM20(NULL),
  fHistoTrueClusMergedEtaEVsM02(NULL),
  fHistoTrueClusBGEvsM20(NULL),
  fHistoTruePrimaryPi0MCPtResolPt(NULL),
  fHistoTruePrimaryPi0PureMergedMCPtResolPt(NULL),
  fHistoTruePrimaryPi0MergedPartConvMCPtResolPt(NULL),
  fHistoTruePrimaryPi01GammaMCPtResolPt(NULL),
  fHistoTruePrimaryPi01ElectronMCPtResolPt(NULL),
  fHistoTruePrimaryEtaMCPtResolPt(NULL),
  fHistoTrueSecondaryPi0MCPtResolPt(NULL),
  fHistoDoubleCountTruePi0PtvsM02(NULL),
  fHistoDoubleCountTrueSecPi0Pt(NULL),
  fHistoDoubleCountTrueMultiplePi0PtvsM02(NULL),
  fHistoDoubleCountTrueMultipleSecPi0Pt(NULL),
  fHistoDoubleCountTrueEtaPtvsM02(NULL),
  fVectorLabelsLeadingPi0(0),
  fVectorLabelsMultiplePi0(0),
  fVectorLabelsMultiplePi0Reduced(0),
  fVectorDoubleCountTruePi0s(0),
  fVectorDoubleCountTrueMultilePi0s(0),
  fVectorDoubleCountTrueEtas(0),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoVertexZ(NULL),
  fHistoOverlapsPi0All(NULL),
  fHistoOverlapsPi0Accepted(NULL),
  fHistoNClusterCandidates(NULL),
  fHistoNClusterMergedCandidates(NULL),
  fHistoNGoodESDTracksVsNClusterCandidates(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoNV0Tracks(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoJetJetNTrials(NULL),
  fVectorJetPt(0),
  fVectorJetPx(0),
  fVectorJetPy(0),
  fVectorJetPz(0),
  fVectorJetEta(0),
  fVectorJetPhi(0),
  fVectorJetArea(0),
  fTrueVectorJetPt(0),
  fTrueVectorJetPx(0),
  fTrueVectorJetPy(0),
  fTrueVectorJetPz(0),
  fTrueVectorJetEta(0),
  fTrueVectorJetPhi(0),
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
  fWeightJetJetMC(1),
  fSelectedMesonID(0),
  fEnableDetailedPrintOut(kFALSE),
  fEnableSortForClusMC(kFALSE),
  tBrokenFiles(NULL),
  fFileNameBroken(NULL),
  fDoDetailedM02(kFALSE),
  fTrackMatcherRunningMode(0),
  fMinAllowedPi0OverlapsMC(-1),
  fMaxAllowedPi0OverlapsMC(-1),
  fHistoPi0EvsGammaOverlapE(NULL)
{

}

//________________________________________________________________________
AliAnalysisTaskGammaCaloMerged::AliAnalysisTaskGammaCaloMerged(const char *name):
  AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fCorrTaskSetting(""),
  fDoLightOutput(kFALSE),
  fInputEvent(NULL),
  fMCEvent(NULL),
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
  fOutlierJetReader(NULL),
  fConvJetReader(NULL),
  fDoJetAnalysis(kFALSE),
  fDoJetQA(kFALSE),
  fDoOutOfJet(0),
  fJetHistograms(NULL),
  fAODMCTrackArray(NULL),
  farrClustersProcess(NULL),
  fMapNeutralPionOverlap(NULL),
  fMesonOLFromCluster(0),
  fHistoMotherInvMassPt(NULL),
  fHistoMotherPtY(NULL),
  fHistoMotherPtAlpha(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusGammaE(NULL),
  fHistoClusOverlapHeadersGammaPt(NULL),
  fHistoClusNLMPt(NULL),
  fHistoClusMergedPtvsM02(NULL),
  fHistoClusMergedEvsM02(NULL),
  fHistoClusMergedPtvsM02Accepted(NULL),
  fHistoClusMergedEvsM02Accepted(NULL),
  fHistoClusMergedEvsM20(NULL),
  fHistoClusNCellsPt(NULL),
  fHistoClusMergedNCellsPt(NULL),
  fHistoClusMergedNParticlePt(NULL),
  fHistoClusMergedNCellsAroundPt(NULL),
  fHistoClusMergedNCellsAroundAndInPt(NULL),
  fHistoClusMergedEAroundE(NULL),
  fHistoPtJet(NULL),
  fHistoJetEta(NULL),
  fHistoJetPhi(NULL),
  fHistoJetArea(NULL),
  fClusterEtaPhiJets(NULL),
  fHistoNJets(NULL),
  fHistoEventwJets(NULL),
  fHistoTruevsRecJetPt(NULL),
  fHistoClusMergedPtvsRJetAccepted(NULL),
  fHistoJetFragmFunc(NULL),
  fHistoMCHeaders(NULL),
  fHistoMCPi0Pt(NULL),
  fHistoMCPi0ReducedPt(NULL),
  fHistoMCPi0WOWeightPt(NULL),
  fHistoMCPi0WOEvtWeightPt(NULL),
  fHistoMCEtaPt(NULL),
  fHistoMCEtaWOWeightPt(NULL),
  fHistoMCEtaWOEvtWeightPt(NULL),
  fHistoMCPi0DalitzPt(NULL),
  fHistoMCPi0DalitzWOWeightPt(NULL),
  fHistoMCPi0DalitzWOEvtWeightPt(NULL),
  fHistoMCEtaDalitzPt(NULL),
  fHistoMCEtaDalitzWOWeightPt(NULL),
  fHistoMCEtaDalitzWOEvtWeightPt(NULL),
  fHistoMCPi0InAccPt(NULL),
  fHistoMCPi0ReducedInAccPt(NULL),
  fHistoMCEtaInAccPt(NULL),
  fHistoMCPi0WOEvtWeightInAccPt(NULL),
  fHistoMCEtaWOEvtWeightInAccPt(NULL),
  fHistoMCPi0DalitzInAccPt(NULL),
  fHistoMCEtaDalitzInAccPt(NULL),
  fHistoMCPi0DalitzWOEvtWeightInAccPt(NULL),
  fHistoMCEtaDalitzWOEvtWeightInAccPt(NULL),
  fHistoMCSecPi0PtvsSource(NULL),
  fHistoMCSecPi0InAccPtvsSource(NULL),
  fHistoMCPi0PtJetPt(NULL),
  fHistoMCEtaPtJetPt(NULL),
  fHistoMCPrimaryPtvsSource(NULL),
  fHistoMCPrimaryYvsSource(NULL),
  fHistoMCDecayGammaPt(NULL),
  fHistoMCAllGammaPt(NULL),
  fHistoTrueClusEFracFirstLabel(NULL),
  fHistoTrueClusEFracLeadingPi0(NULL),
  fHistoTrueClusNeutralContamination(NULL),
  fHistoTrueClusPhotonContamination(NULL),
  fHistoTrueClusElectronContamination(NULL),
  fHistoTrueClusOtherContamination(NULL),
  fHistoTrueClusMergedPtvsM02(NULL),
  fHistoTrueClusPi0PtvsM02(NULL),
  fHistoTrueClusMultiplePi0PtvsM02(NULL),
  fHistoTrueClusPi0DalitzPtvsM02(NULL),
  fHistoTrueClusPrimPi0PtvsM02(NULL),
  fHistoTrueClusPrimPi0PtMCPt(NULL),
  fHistoTruePureMergedClusPrimPi0PtvsM02(NULL),
  fHistoTruePartConvClusPrimPi0PtvsM02(NULL),
  fHistoTrueClusMergedPartConvPi0EVsM20(NULL),
  fHistoTrueClusPrimPi0PureMergedPtMCPt(NULL),
  fHistoTrueClusPrimPi0PartConvMergedPtMCPt(NULL),
  fHistoTrueClusPrimPi01GammaMergedPtMCPt(NULL),
  fHistoTrueClusPrimPi01ElectronMergedPtMCPt(NULL),
  fHistoTrueClusMultiplePrimPi0PtvsM02(NULL),
  fHistoTrueClusSecPi0PtvsM02(NULL),
  fHistoTrueClusSecPi0FromK0sPtvsM02(NULL),
  fHistoTrueClusSecPi0FromK0lPtvsM02(NULL),
  fHistoTrueClusSecPi0FromLambdaPtvsM02(NULL),
  fHistoTrueClusEtaPtvsM02(NULL),
  fHistoTrueClusEtaDalitzPtvsM02(NULL),
  fHistoTrueClusMergedPureFromPi0PtvsM02(NULL),
  fHistoTrueClusMergedPureFromEtaPtvsM02(NULL),
  fHistoTrueClusMergedPartConvFromPi0PtvsM02(NULL),
  fHistoTrueClusMergedPartConvFromEtaPtvsM02(NULL),
  fHistoNTrueMultiplePi0vsPt(NULL),
  fHistoTrueClusGammaFromPi0PtvsM02(NULL),
  fHistoTrueClusGammaFromEtaPtvsM02(NULL),
  fHistoTrueClusElectronFromPi0PtvsM02(NULL),
  fHistoTrueClusElectronFromEtaPtvsM02(NULL),
  fHistoTrueSecPi0PtvsDiffReco(NULL),
  fHistoTrueClusBGPtvsM02(NULL),
  fHistoTrueClusGammaPtvsM02(NULL),
  fHistoTrueClusGammaPartConvPtvsM02(NULL),
  fHistoTrueClusElectronPtvsM02(NULL),
  fHistoTrueClusElectronFromGammaPtvsM02(NULL),
  fHistoTrueClusMergedInvMassvsPt(NULL),
  fHistoTrueClusPi0InvMassvsPt(NULL),
  fHistoTrueClusPrimPi0InvMassvsPt(NULL),
  fHistoTrueClusEtaInvMassvsPt(NULL),
  fHistoTrueClusBGInvMassvsPt(NULL),
  fHistoTrueClusGammaInvMassvsPt(NULL),
  fHistoTrueClusElectronInvMassvsPt(NULL),
  fHistoTrueClusBGPtvsSource(NULL),
  fHistoTrueClusGammaPtvsSource(NULL),
  fHistoTrueClusElectronPtvsSource(NULL),
//  fHistoTrueClusElectronPtvsMotherID(NULL),
//  fHistoTrueClusElectronPtvsTopMotherID(NULL),
//  fHistoTrueClusElectronPtvsConvPhotonTopMotherID(NULL),
  fHistoTrueMergedMissedPDG(NULL),
  fHistoTrueClusMergedPtvsRJet(NULL),
  fHistoTrueClusPi0PtvsRJet(NULL),
  fHistoTrueClusEtaPtvsRJet(NULL),
  fHistoTruePi0PtY(NULL),
  fHistoTrueEtaPtY(NULL),
  fHistoTruePi0PtAlpha(NULL),
  fHistoTrueEtaPtAlpha(NULL),
  fHistoTrueClusGammaEM02(NULL),
  fHistoTrueClusElectronEM02(NULL),
  fHistoTrueClusPi0EM02(NULL),
  fHistoTrueClusEtaEM02(NULL),
  fHistoTrueClusBGEvsM02(NULL),
  fHistoTrueClusGammaEvsM20(NULL),
  fHistoTrueClusElectronEM20(NULL),
  fHistoTrueClusMergedPi0EVsM20(NULL),
  fHistoTrueClusMergedPi0EVsM02(NULL),
  fHistoTrueClusMergedEtaEVsM20(NULL),
  fHistoTrueClusMergedEtaEVsM02(NULL),
  fHistoTrueClusBGEvsM20(NULL),
  fHistoTruePrimaryPi0MCPtResolPt(NULL),
  fHistoTruePrimaryPi0PureMergedMCPtResolPt(NULL),
  fHistoTruePrimaryPi0MergedPartConvMCPtResolPt(NULL),
  fHistoTruePrimaryPi01GammaMCPtResolPt(NULL),
  fHistoTruePrimaryPi01ElectronMCPtResolPt(NULL),
  fHistoTruePrimaryEtaMCPtResolPt(NULL),
  fHistoTrueSecondaryPi0MCPtResolPt(NULL),
  fHistoDoubleCountTruePi0PtvsM02(NULL),
  fHistoDoubleCountTrueSecPi0Pt(NULL),
  fHistoDoubleCountTrueMultiplePi0PtvsM02(NULL),
  fHistoDoubleCountTrueMultipleSecPi0Pt(NULL),
  fHistoDoubleCountTrueEtaPtvsM02(NULL),
  fVectorLabelsLeadingPi0(0),
  fVectorLabelsMultiplePi0(0),
  fVectorLabelsMultiplePi0Reduced(0),
  fVectorDoubleCountTruePi0s(0),
  fVectorDoubleCountTrueMultilePi0s(0),
  fVectorDoubleCountTrueEtas(0),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoVertexZ(NULL),
  fHistoOverlapsPi0All(NULL),
  fHistoOverlapsPi0Accepted(NULL),
  fHistoNClusterCandidates(NULL),
  fHistoNClusterMergedCandidates(NULL),
  fHistoNGoodESDTracksVsNClusterCandidates(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoNV0Tracks(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoJetJetNTrials(NULL),
  fVectorJetPt(0),
  fVectorJetPx(0),
  fVectorJetPy(0),
  fVectorJetPz(0),
  fVectorJetEta(0),
  fVectorJetPhi(0),
  fVectorJetArea(0),
  fTrueVectorJetPt(0),
  fTrueVectorJetPx(0),
  fTrueVectorJetPy(0),
  fTrueVectorJetPz(0),
  fTrueVectorJetEta(0),
  fTrueVectorJetPhi(0),
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
  fWeightJetJetMC(1),
  fSelectedMesonID(0),
  fEnableDetailedPrintOut(kFALSE),
  fEnableSortForClusMC(kFALSE),
  tBrokenFiles(NULL),
  fFileNameBroken(NULL),
  fDoDetailedM02(kFALSE),
  fTrackMatcherRunningMode(0),
  fMinAllowedPi0OverlapsMC(-1),
  fMaxAllowedPi0OverlapsMC(-1),
  fHistoPi0EvsGammaOverlapE(NULL)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaCaloMerged::~AliAnalysisTaskGammaCaloMerged()
{

}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::UserCreateOutputObjects(){

  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader

  if(((AliConvEventCuts*)fEventCutArray->At(0))->GetUseJetFinderForOutliers()){
    fOutlierJetReader=(AliAnalysisTaskJetOutlierRemoval*)AliAnalysisManager::GetAnalysisManager()->GetTask("AliAnalysisTaskJetOutlierRemoval");
    if(!fOutlierJetReader){AliFatal("Error: No AliAnalysisTaskJetOutlierRemoval");} // GetV0Reader
    else{printf("Found AliAnalysisTaskJetOutlierRemoval used for outlier removal!\n");}
  }

  if( ((AliConversionMesonCuts*)fMesonCutArray->At(0))->DoJetAnalysis())  fDoJetAnalysis = kTRUE;
  if( ((AliConversionMesonCuts*)fMesonCutArray->At(0))->DoJetQA())        fDoJetQA       = kTRUE;
  if( ((AliConversionMesonCuts*)fMesonCutArray->At(0))->DoOutOfJet() > 0 )fDoOutOfJet    = ((AliConversionMesonCuts*)fMesonCutArray->At(0))->DoOutOfJet();

  if(fDoJetAnalysis){
    fConvJetReader=(AliAnalysisTaskConvJet*)AliAnalysisManager::GetAnalysisManager()->GetTask("AliAnalysisTaskConvJet");
    if(!fConvJetReader){printf("Error: No AliAnalysisTaskConvJet");return;}
  }

  fMapNeutralPionOverlap = new std::map<Int_t, Int_t>[fnCuts];

  Int_t invMassBins                           = 800;
  Float_t startMass                           = 0;
  Float_t endMass                             = 0.8;
  Double_t *arrPtBinning                      = new Double_t[1000];

  if (GetSelectedMesonID() == 1) {
    invMassBins                               = 400;
    startMass                                 = 0;
    endMass                                   = 0.4;
  } else if (GetSelectedMesonID() == 2) {
    invMassBins                               = 800;
    startMass                                 = 0.;
    endMass                                   = 0.8;
  }

  Int_t ptBinsDefClus                         = 500;
  Float_t startPtDefClus                      = 0;
  Float_t endPtDefClus                        = 50;

  Int_t ptBins                                = 400;
  Float_t startPt                             = 10;
  Float_t endPt                               = 50;

  Int_t ptBinsLog                             = 200;
  Float_t startPtLog                          = 10;
  Float_t endPtLog                            = 50;
  if (GetSelectedMesonID() == 2 ) {
    ptBins                                    = 500;
    startPt                                   = 20;
    endPt                                     = 70;
    for(Int_t i=0; i<ptBins+1;i++){
      arrPtBinning[i]                         = 10.+((endPt-startPt)/ptBins)*i;
    }
    ptBinsLog                                 = 250;
    startPtLog                                = 20;
    endPtLog                                  = 70;
  } else if (((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k5TeV ||
             ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k7TeV ||
             ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k8TeV ||
             ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb8TeV ||
             ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb5TeVR2  ||
             ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kpPb5TeV ||
             ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kPbPb5TeV){
    ptBins                                    = 120;
    startPt                                   = 0;
    endPt                                     = 200;
    // pT dependent binning for very high pT analyses
    for(Int_t i=0; i<ptBins+1;i++){
      if(i<100)      arrPtBinning[i]          = 1.0*i;
      else if(i<120) arrPtBinning[i]          = 100.+5*(i-100);
      else           arrPtBinning[i]          = endPt;
    }
    ptBinsLog                                 = 475;
    startPtLog                                = 10;
    endPtLog                                  = 200;
    ptBinsDefClus                             = 1000;
    startPtDefClus                            = 0;
    endPtDefClus                              = 200;
  } else if (((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k13TeVLowB  ){
    ptBins                                    = 900;
    startPt                                   = 10;
    endPt                                     = 100;
    for(Int_t i=0; i<ptBins+1;i++){
      arrPtBinning[i]                         = 10.+((endPt-startPt)/ptBins)*i;
    }
    ptBinsLog                                 = 450;
    startPtLog                                = 10;
    endPtLog                                  = 100;
    ptBinsDefClus                             = 1000;
    startPtDefClus                            = 0;
    endPtDefClus                              = 100;
  } else if (((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::k13TeV){
    ptBins                                    = 160;
    startPt                                   = 0;
    endPt                                     = 400;
    // pT dependent binning for very high pT (400GeV) analyses
    for(Int_t i=0; i<ptBins+1;i++){
      if(i<100)      arrPtBinning[i]          = 1.0*i;
      else if(i<160) arrPtBinning[i]          = 100.+5*(i-100);
      else           arrPtBinning[i]          = endPt;
    }
    ptBinsLog                                 = 780;
    startPtLog                                = 10;
    endPtLog                                  = 400;
    ptBinsDefClus                             = 800;
    startPtDefClus                            = 0;
    endPtDefClus                              = 400;
  }else {
    for(Int_t i=0; i<ptBins+1;i++){
      arrPtBinning[i]                          = 10.+((endPt-startPt)/ptBins)*i;
    }
  }

  Int_t showerShapeBins                       = 500;
  Float_t startShowerShape                    = 0;
  Float_t endShowerShape                      = 5;

  if (fDoDetailedM02) showerShapeBins         = 5000;

  // Create histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer                          = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer                          = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

  fCutFolder                                  = new TList*[fnCuts];
  fESDList                                    = new TList*[fnCuts];

  fHistoNEvents                               = new TH1F*[fnCuts];
  if(fIsMC == 2){
    fHistoNEventsWOWeight                     = new TH1F*[fnCuts];
    fProfileJetJetXSection                    = new TProfile*[fnCuts];
    fHistoJetJetNTrials                       = new TH1F*[fnCuts];
  }

  fHistoNGoodESDTracks                        = new TH1F*[fnCuts];
  fHistoVertexZ                               = new TH1F*[fnCuts];
  fHistoNClusterCandidates                    = new TH1F*[fnCuts];
  fHistoNClusterMergedCandidates              = new TH1F*[fnCuts];
  fHistoNGoodESDTracksVsNClusterCandidates    = new TH2F*[fnCuts];
  fHistoSPDClusterTrackletBackground          = new TH2F*[fnCuts];
  fHistoNV0Tracks                             = new TH1F*[fnCuts];
  fProfileEtaShift                            = new TProfile*[fnCuts];

  fHistoMotherInvMassPt                       = new TH2F*[fnCuts];
  if (fDoMesonQA > 0 ){
    fHistoMotherPtY                           = new TH2F*[fnCuts];
    fHistoMotherPtAlpha                       = new TH2F*[fnCuts];
  }

  fHistoClusGammaPt                           = new TH1F*[fnCuts];
  fHistoClusGammaE                            = new TH1F*[fnCuts];
  fHistoClusOverlapHeadersGammaPt             = new TH1F*[fnCuts];
  fHistoClusMergedPtvsM02                     = new TH2F*[fnCuts];
  fHistoClusMergedPtvsM02Accepted             = new TH2F*[fnCuts];
  fHistoClusMergedEvsM02Accepted              = new TH2F*[fnCuts];
  fHistoClusNLMPt                             = new TH2F*[fnCuts];
  if(fDoMesonQA > 1){
    fHistoClusMergedEvsM02                    = new TH2F*[fnCuts];
    fHistoClusMergedEvsM20                    = new TH2F*[fnCuts];
  }
  if (fDoClusterQA > 0){
    fHistoClusNCellsPt                        = new TH2F*[fnCuts];
    fHistoClusMergedNCellsPt                  = new TH2F*[fnCuts];
    fHistoClusMergedNCellsAroundPt            = new TH2F*[fnCuts];
    fHistoClusMergedNCellsAroundAndInPt       = new TH2F*[fnCuts];
    fHistoClusMergedEAroundE                  = new TH2F*[fnCuts];

    if (fIsMC > 0){
      fHistoClusMergedNParticlePt             = new TH2F*[fnCuts];
      fHistoOverlapsPi0Accepted               = new TH2F*[fnCuts];
      fHistoOverlapsPi0All                    = new TH2F*[fnCuts];
    }

  }

  if(fDoJetAnalysis){
    fJetHistograms            = new TList*[fnCuts];

    fHistoPtJet               = new TH1F*[fnCuts];
    fHistoJetEta              = new TH1F*[fnCuts];
    fHistoJetPhi              = new TH1F*[fnCuts];
    fHistoJetArea             = new TH1F*[fnCuts];
    fClusterEtaPhiJets        = new TH2F*[fnCuts];
    fHistoNJets               = new TH1F*[fnCuts];
    if(fDoJetQA){
      fHistoEventwJets          = new TH1F*[fnCuts];
    }
    if(!fDoLightOutput){
      fHistoTruevsRecJetPt             = new TH2F*[fnCuts];
      if(fDoOutOfJet != 1 ){
        fHistoClusMergedPtvsRJetAccepted = new TH2F*[fnCuts];
        fHistoJetFragmFunc               = new TH2F*[fnCuts];
      }
    }
  }

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    TString cutstringEvent                        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringCalo                         = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
    TString cutstringCaloMerged                   = ((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))->GetCutNumber();
    TString cutstringMeson                        = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

    //Int_t nLMCut   = ((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))->GetMinNLMCut();

    fCutFolder[iCut]                              = new TList();
    fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s_%s",cutstringEvent.Data(), cutstringCalo.Data(), cutstringCaloMerged.Data(), cutstringMeson.Data()));
    fCutFolder[iCut]->SetOwner(kTRUE);
    fOutputContainer->Add(fCutFolder[iCut]);
    fESDList[iCut]                                = new TList();
    fESDList[iCut]->SetName(Form("%s_%s_%s_%s ESD histograms",cutstringEvent.Data(), cutstringCalo.Data(), cutstringCaloMerged.Data(), cutstringMeson.Data()));
    fESDList[iCut]->SetOwner(kTRUE);
    fCutFolder[iCut]->Add(fESDList[iCut]);

    fHistoNEvents[iCut]                           = new TH1F("NEvents","NEvents",14,-0.5,13.5);
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
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
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL/TPC problems");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fESDList[iCut]->Add(fHistoNEvents[iCut]);

    if (fIsMC == 2){
      fHistoNEventsWOWeight[iCut]                 = new TH1F("NEventsWOWeight","NEventsWOWeight",14,-0.5,13.5);
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
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
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
      fESDList[iCut]->Add(fHistoNEventsWOWeight[iCut]);

      fProfileJetJetXSection[iCut]                = new TProfile("XSection","XSection",1,-0.5,0.5);
      fESDList[iCut]->Add(fProfileJetJetXSection[iCut]);
      fHistoJetJetNTrials[iCut]                   = new TH1F("NTrials","#sum{NTrials}",1,0,1);
      fHistoJetJetNTrials[iCut]->GetXaxis()->SetBinLabel(1,"#sum{NTrials}");
      fESDList[iCut]->Add(fHistoJetJetNTrials[iCut]);
    }

    if(fIsHeavyIon == 1)
      fHistoNGoodESDTracks[iCut]                  = new TH1F("GoodESDTracks","GoodESDTracks",4000,-0.5,3999.5);
    else if(fIsHeavyIon == 2)
      fHistoNGoodESDTracks[iCut]                  = new TH1F("GoodESDTracks","GoodESDTracks",400,-0.5,399.5);
    else
      fHistoNGoodESDTracks[iCut]                  = new TH1F("GoodESDTracks","GoodESDTracks",200,-0.5,199.5);
    fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);

    fHistoVertexZ[iCut]                           = new TH1F("VertexZ","VertexZ",1000,-50,50);
    fESDList[iCut]->Add(fHistoVertexZ[iCut]);

    if(fIsHeavyIon == 1)
      fHistoNClusterCandidates[iCut]              = new TH1F("GammaCandidates","GammaCandidates",600,-0.5,599.5);
    else if(fIsHeavyIon == 2)
      fHistoNClusterCandidates[iCut]              = new TH1F("GammaCandidates","GammaCandidates",400,-0.5,399.5);
    else
      fHistoNClusterCandidates[iCut]              = new TH1F("GammaCandidates","GammaCandidates",100,-0.5,99.5);
    fESDList[iCut]->Add(fHistoNClusterCandidates[iCut]);

    if(fIsHeavyIon == 1)
      fHistoNClusterMergedCandidates[iCut]        = new TH1F("MergedCandidates","MergedCandidates",300,-0.5,299.5);
    else if(fIsHeavyIon == 2)
      fHistoNClusterMergedCandidates[iCut]        = new TH1F("MergedCandidates","MergedCandidates",200,-0.5,199.5);
    else
      fHistoNClusterMergedCandidates[iCut]        = new TH1F("MergedCandidates","MergedCandidates",50,-0.5,49.5);
    fESDList[iCut]->Add(fHistoNClusterMergedCandidates[iCut]);


    if(fIsHeavyIon == 1)
      fHistoNGoodESDTracksVsNClusterCandidates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",4000,-0.5,3999.5,600,-0.5,599.5);
    else if(fIsHeavyIon == 2)
      fHistoNGoodESDTracksVsNClusterCandidates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",400,-0.5,399.5,400,-0.5,399.5);
    else
      fHistoNGoodESDTracksVsNClusterCandidates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",200,-0.5,199.5,100,-0.5,99.5);
    fESDList[iCut]->Add(fHistoNGoodESDTracksVsNClusterCandidates[iCut]);

    fHistoSPDClusterTrackletBackground[iCut]      = new TH2F("SPD tracklets vs SPD clusters","SPD tracklets vs SPD clusters",100,0,200,250,0,1000);
    fESDList[iCut]->Add(fHistoSPDClusterTrackletBackground[iCut]);

    if(fIsHeavyIon == 1)
      fHistoNV0Tracks[iCut]                       = new TH1F("V0 Multiplicity","V0 Multiplicity",30000,-0.5,29999.5);
    else if(fIsHeavyIon == 2)
      fHistoNV0Tracks[iCut]                       = new TH1F("V0 Multiplicity","V0 Multiplicity",2500,-0.5,2499.5);
    else
      fHistoNV0Tracks[iCut]                       = new TH1F("V0 Multiplicity","V0 Multiplicity",1500,-0.5,1499.5);
    fESDList[iCut]->Add(fHistoNV0Tracks[iCut]);
    fProfileEtaShift[iCut]                        = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
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

    fHistoClusGammaPt[iCut]                       = new TH1F("ClusGamma_Pt","ClusGamma_Pt",ptBinsDefClus, startPtDefClus, endPtDefClus);
    fESDList[iCut]->Add(fHistoClusGammaPt[iCut]);
    fHistoClusGammaE[iCut]                        = new TH1F("ClusGamma_E","ClusGamma_E",ptBinsDefClus, startPtDefClus, endPtDefClus);
    fESDList[iCut]->Add(fHistoClusGammaE[iCut]);
    fHistoClusOverlapHeadersGammaPt[iCut]         = new TH1F("ClusGammaOverlapHeaders_Pt","ClusGammaOverlapHeaders_Pt",ptBins, arrPtBinning);
    fESDList[iCut]->Add(fHistoClusOverlapHeadersGammaPt[iCut]);
    fHistoClusMergedPtvsM02[iCut]                 = new TH2F("ClusMerged_Pt_M02","ClusMerged_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
    fESDList[iCut]->Add(fHistoClusMergedPtvsM02[iCut]);
    fHistoClusMergedPtvsM02Accepted[iCut]         = new TH2F("ClusMerged_Pt_M02_AcceptedMeson","ClusMerged_Pt_M02_AcceptedMeson",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
    fESDList[iCut]->Add(fHistoClusMergedPtvsM02Accepted[iCut]);
    fHistoClusMergedEvsM02Accepted[iCut]          = new TH2F("ClusMerged_E_M02_AcceptedMeson","ClusMerged_E_M02_AcceptedMeson",ptBins, arrPtBinning, showerShapeBins, startShowerShape, endShowerShape);
    fESDList[iCut]->Add(fHistoClusMergedEvsM02Accepted[iCut]);
    fHistoClusNLMPt[iCut]                         = new TH2F("ClusMerged_NLM_Pt_AcceptedMeson","ClusMerged_NLM_Pt_AcceptedMeson",12, -0.5, 11.5, ptBins, arrPtBinning);
    fESDList[iCut]->Add(fHistoClusNLMPt[iCut]);

    if(fDoMesonQA > 1){
      fHistoClusMergedEvsM02[iCut]                 = new TH2F("ClusMerged_E_M02","ClusMerged_E_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fESDList[iCut]->Add(fHistoClusMergedEvsM02[iCut]);
      fHistoClusMergedEvsM20[iCut]                 = new TH2F("ClusMerged_E_M20","ClusMerged_E_M20",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fESDList[iCut]->Add(fHistoClusMergedEvsM20[iCut]);
    }


    if (fIsMC == 2){
      fHistoClusGammaPt[iCut]->Sumw2();
      fHistoClusGammaE[iCut]->Sumw2();
      fHistoClusOverlapHeadersGammaPt[iCut]->Sumw2();
      fHistoClusMergedPtvsM02[iCut]->Sumw2();
      fHistoClusMergedPtvsM02Accepted[iCut]->Sumw2();
      fHistoClusMergedEvsM02Accepted[iCut]->Sumw2();
      fHistoClusNLMPt[iCut]->Sumw2();
      if(fDoMesonQA > 1){
        fHistoClusMergedEvsM02[iCut]->Sumw2();
        fHistoClusMergedEvsM20[iCut]->Sumw2();
      }
    }

    if (fDoClusterQA > 0){
      fHistoClusNCellsPt[iCut]                      = new TH2F("Clus_NCells_Pt","Clus_NCells_Pt",100,-0.5,99.5,ptBins, arrPtBinning);
      fESDList[iCut]->Add(fHistoClusNCellsPt[iCut]);
      fHistoClusMergedNCellsPt[iCut]                = new TH2F("ClusMerged_NCells_Pt","ClusMerged_NCells_Pt",100,-0.5,99.5,ptBins, arrPtBinning);
      fESDList[iCut]->Add(fHistoClusMergedNCellsPt[iCut]);
      fHistoClusMergedNCellsAroundPt[iCut]          = new TH2F("ClusMerged_NCellsAroundClus_Pt","ClusMerged_NCellsAroundClus_Pt",100,-0.5,99.5,ptBins, arrPtBinning);
      fESDList[iCut]->Add(fHistoClusMergedNCellsAroundPt[iCut]);
      fHistoClusMergedNCellsAroundAndInPt[iCut]     = new TH2F("ClusMerged_NCellsAroundAndInClus_Pt","ClusMerged_NCellsAroundAndInClus_Pt",100,-0.5,99.5,ptBins, arrPtBinning);
      fESDList[iCut]->Add(fHistoClusMergedNCellsAroundAndInPt[iCut]);
      fHistoClusMergedEAroundE[iCut]                = new TH2F("ClusMerged_EAroundClus_E","ClusMerged_EAroundClus_E",ptBinsDefClus, startPtDefClus, endPtDefClus, ptBins, arrPtBinning);
      fESDList[iCut]->Add(fHistoClusMergedEAroundE[iCut]);

      if (fIsMC > 0){
        fHistoClusMergedNParticlePt[iCut]           = new TH2F("ClusMerged_NPart_Pt","ClusMerged_NPart_Pt",100,-0.5,99.5,ptBins, arrPtBinning);
        fESDList[iCut]->Add(fHistoClusMergedNParticlePt[iCut]);
        fHistoOverlapsPi0Accepted[iCut]           = new TH2F("Pi0_NOverlapsAccepted_GenPt","Pi0_NOverlapsAccepted_GenPt",50,-0.5,49.5,ptBins, arrPtBinning);
        fESDList[iCut]->Add(fHistoOverlapsPi0Accepted[iCut]);
        fHistoOverlapsPi0All[iCut]           = new TH2F("Pi0_NOverlapsAll_GenPt","Pi0_NOverlapsAll_GenPt",50,-0.5,49.5,ptBins, arrPtBinning);
        fESDList[iCut]->Add(fHistoOverlapsPi0All[iCut]);
        if (fIsMC == 2){
          fHistoClusNCellsPt[iCut]->Sumw2();
          fHistoClusMergedNCellsPt[iCut]->Sumw2();
          fHistoClusMergedNCellsAroundPt[iCut]->Sumw2();
          fHistoClusMergedNCellsAroundAndInPt[iCut]->Sumw2();
          fHistoClusMergedEAroundE[iCut]->Sumw2();
          fHistoClusMergedNParticlePt[iCut]->Sumw2();
          fHistoClusMergedEAroundE[iCut]->Sumw2();
          fHistoOverlapsPi0Accepted[iCut]->Sumw2();
          fHistoOverlapsPi0All[iCut]->Sumw2();
        }
      }
    }

    fHistoMotherInvMassPt[iCut]                   = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",invMassBins, startMass, endMass,ptBins, arrPtBinning);
    fESDList[iCut]->Add(fHistoMotherInvMassPt[iCut]);

    if (fIsMC == 2){
      fHistoMotherInvMassPt[iCut]->Sumw2();
    }

    if (fDoMesonQA > 0 ){
      fHistoMotherPtY[iCut]                       = new TH2F("ESD_Mother_Pt_Y","ESD_Mother_Pt_Y", ptBinsLog, startPtLog, endPtLog, 150, -1.5, 1.5);
      SetLogBinningXTH2(fHistoMotherPtY[iCut]);
      fESDList[iCut]->Add(fHistoMotherPtY[iCut]);
      fHistoMotherPtAlpha[iCut]                   = new TH2F("ESD_Mother_Pt_Alpha","ESD_Mother_Pt_Alpha", ptBinsLog, startPtLog, endPtLog, 100, 0, 1);
      SetLogBinningXTH2(fHistoMotherPtAlpha[iCut]);
      fESDList[iCut]->Add(fHistoMotherPtAlpha[iCut]);
      if (fIsMC == 2){
        fHistoMotherPtY[iCut]->Sumw2();
        fHistoMotherPtAlpha[iCut]->Sumw2();
      }
    }

    if(fDoJetAnalysis){

      fJetHistograms[iCut] = new TList();
      fJetHistograms[iCut]->SetOwner(kTRUE);
      fJetHistograms[iCut]->SetName(Form("%s_%s_%s Jet histograms", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));

      fHistoPtJet[iCut] = new TH1F("JetPt", "JetPt", 150, 0, 150);
      fJetHistograms[iCut]->Add(fHistoPtJet[iCut]);
      fHistoJetEta[iCut] = new TH1F("JetEta", "JetEta", 100, -1, 1);
      fJetHistograms[iCut]->Add(fHistoJetEta[iCut]);
      fHistoJetPhi[iCut] = new TH1F("JetPhi", "JetPhi", 70, 0, 7);
      fJetHistograms[iCut]->Add(fHistoJetPhi[iCut]);
      fHistoJetArea[iCut] = new TH1F("JetArea", "JetArea", 50, 0, 1);
      fJetHistograms[iCut]->Add(fHistoJetArea[iCut]);
      fHistoNJets[iCut] = new TH1F("NJets", "NJets", 10, 0, 10);
      fJetHistograms[iCut]->Add(fHistoNJets[iCut]);

      if(fIsMC > 0){
        if(fDoJetQA){
          fHistoEventwJets[iCut] = new TH1F("NEvents_with_Jets", "NEvents_with_Jets", 3, 0, 3);
          fHistoEventwJets[iCut]->GetXaxis()->SetBinLabel(1,"rec. > 10GeV");
          fHistoEventwJets[iCut]->GetXaxis()->SetBinLabel(2,"rec. < 10GeV, true > 10GeV");
          fHistoEventwJets[iCut]->GetXaxis()->SetBinLabel(3,"rec. > 10GeV, true < 10GeV");
          fJetHistograms[iCut]->Add(fHistoEventwJets[iCut]);
        }
      }

      if(!fDoLightOutput){
        fClusterEtaPhiJets[iCut] = new TH2F("JetEtaPhiMap", "JetEtaPhiMap", 110, -0.7, 0.7, 462, 0, 2*TMath::Pi());
        fJetHistograms[iCut]->Add(fClusterEtaPhiJets[iCut]);
        if(fDoOutOfJet != 1 ){
          fHistoClusMergedPtvsRJetAccepted[iCut] = new TH2F("ESD_ClusMerged_Pt_RJet","ESD_ClusMerged_Pt_RJet",ptBins, arrPtBinning,50, 0, 1.);
          fJetHistograms[iCut]->Add(fHistoClusMergedPtvsRJetAccepted[iCut]);
          fHistoJetFragmFunc[iCut] = new TH2F("ESD_Pi0inJetPt_FragmentationFunc","ESD_Pi0inJetPt_FragmentationFunc",ptBins, arrPtBinning,50, 0, 1);
          fJetHistograms[iCut]->Add(fHistoJetFragmFunc[iCut]);
        }
        if(fIsMC > 0){
          fHistoTruevsRecJetPt[iCut] = new TH2F("True_JetPt_vs_Rec_JetPt", "True_JetPt_vs_Rec_JetPt", ptBins, arrPtBinning, ptBins, arrPtBinning);
          fJetHistograms[iCut]->Add(fHistoTruevsRecJetPt[iCut]);
        }
      }

    if (fIsMC == 2){
      fHistoPtJet[iCut]->Sumw2();
      fHistoJetEta[iCut]->Sumw2();
      fHistoJetPhi[iCut]->Sumw2();
      fHistoJetArea[iCut]->Sumw2();
      fHistoNJets[iCut]->Sumw2();
      if(fDoJetQA){
        fHistoEventwJets[iCut]->Sumw2();
      }
      if(!fDoLightOutput){
        fClusterEtaPhiJets[iCut]->Sumw2();
        fHistoTruevsRecJetPt[iCut]->Sumw2();
        if(fDoOutOfJet != 1 ){
          fHistoClusMergedPtvsRJetAccepted[iCut]->Sumw2();
          fHistoJetFragmFunc[iCut]->Sumw2();
        }
      }
    }
  }
}

  if (fIsMC > 0){
    fMCList = new TList*[fnCuts];
    fTrueList = new TList*[fnCuts];

    fHistoMCHeaders                               = new TH1I*[fnCuts];
    fHistoMCPi0Pt                                 = new TH1F*[fnCuts];
    fHistoMCPi0ReducedPt                          = new TH1F*[fnCuts];
    fHistoMCEtaPt                                 = new TH1F*[fnCuts];
    fHistoMCPrimaryPtvsSource                     = new TH2F*[fnCuts];
    fHistoMCPrimaryYvsSource                      = new TH2F*[fnCuts];
    fHistoMCDecayGammaPt                          = new TH1F*[fnCuts];
    fHistoMCAllGammaPt                            = new TH1F*[fnCuts];
    if (GetSelectedMesonID() != 2){
      fHistoMCPi0WOWeightPt                       = new TH1F*[fnCuts];
      fHistoMCPi0InAccPt                          = new TH1F*[fnCuts];
      fHistoMCPi0ReducedInAccPt                   = new TH1F*[fnCuts];
      fHistoMCPi0WOEvtWeightInAccPt               = new TH1F*[fnCuts];
      fHistoMCSecPi0PtvsSource                    = new TH2F*[fnCuts];
      fHistoMCSecPi0InAccPtvsSource               = new TH2F*[fnCuts];
      fHistoMCPi0DalitzPt                         = new TH1F*[fnCuts];
      fHistoMCPi0DalitzWOWeightPt                 = new TH1F*[fnCuts];
      fHistoMCPi0DalitzInAccPt                    = new TH1F*[fnCuts];
      fHistoMCPi0DalitzWOEvtWeightInAccPt         = new TH1F*[fnCuts];
            fHistoPi0EvsGammaOverlapE             = new TH2F*[fnCuts];

    }
    if (GetSelectedMesonID() != 1){
      fHistoMCEtaWOWeightPt                       = new TH1F*[fnCuts];
      fHistoMCEtaInAccPt                          = new TH1F*[fnCuts];
      fHistoMCEtaWOEvtWeightInAccPt               = new TH1F*[fnCuts];
      fHistoMCEtaDalitzPt                         = new TH1F*[fnCuts];
      fHistoMCEtaDalitzWOWeightPt                 = new TH1F*[fnCuts];
      fHistoMCEtaDalitzInAccPt                    = new TH1F*[fnCuts];
      fHistoMCEtaDalitzWOEvtWeightInAccPt         = new TH1F*[fnCuts];
    }
    if (fIsMC == 2){
      if (GetSelectedMesonID() != 2){
        fHistoMCPi0WOEvtWeightPt                  = new TH1F*[fnCuts];
        fHistoMCPi0DalitzWOEvtWeightPt            = new TH1F*[fnCuts];
      }
      if (GetSelectedMesonID() != 1){
        fHistoMCEtaWOEvtWeightPt                  = new TH1F*[fnCuts];
        fHistoMCEtaDalitzWOEvtWeightPt            = new TH1F*[fnCuts];
      }
      if (fDoMesonQA > 0){
        if (GetSelectedMesonID() != 2)
          fHistoMCPi0PtJetPt                      = new TH2F*[fnCuts];
        if (GetSelectedMesonID() != 1)
          fHistoMCEtaPtJetPt                      = new TH2F*[fnCuts];
      }
    }

    fHistoTrueClusEFracFirstLabel                 = new TH2F*[fnCuts];
    fHistoTrueClusEFracLeadingPi0                 = new TH2F*[fnCuts];
    fHistoTrueClusMergedPtvsM02                   = new TH2F*[fnCuts];
    fHistoTrueClusPi0PtvsM02                      = new TH2F*[fnCuts];
    fHistoTrueClusMultiplePi0PtvsM02              = new TH2F*[fnCuts];
    fHistoTrueClusPi0DalitzPtvsM02                = new TH2F*[fnCuts];
    if (GetSelectedMesonID() < 2){
      fHistoTrueClusPrimPi0PtvsM02                = new TH2F*[fnCuts];
      fHistoTrueClusPrimPi0PtMCPt                 = new TH2F*[fnCuts];
      fHistoTruePureMergedClusPrimPi0PtvsM02                 = new TH2F*[fnCuts];
      fHistoTruePartConvClusPrimPi0PtvsM02                 = new TH2F*[fnCuts];
      if (fDoMesonQA > 1){
        fHistoTrueClusMergedPartConvPi0EVsM20     = new TH2F*[fnCuts];
        fHistoTrueClusPrimPi0PureMergedPtMCPt     = new TH2F*[fnCuts];
        fHistoTrueClusPrimPi0PartConvMergedPtMCPt = new TH2F*[fnCuts];
        fHistoTrueClusPrimPi01GammaMergedPtMCPt   = new TH2F*[fnCuts];
        fHistoTrueClusPrimPi01ElectronMergedPtMCPt= new TH2F*[fnCuts];
        fHistoTrueClusNeutralContamination            = new TH2F*[fnCuts];
        fHistoTrueClusPhotonContamination             = new TH2F*[fnCuts];
        fHistoTrueClusElectronContamination           = new TH2F*[fnCuts];
        fHistoTrueClusOtherContamination              = new TH2F*[fnCuts];
      }
      fHistoTrueClusMultiplePrimPi0PtvsM02        = new TH2F*[fnCuts];
      fHistoTrueClusSecPi0PtvsM02                 = new TH2F*[fnCuts];
      fHistoTrueClusSecPi0FromK0sPtvsM02          = new TH2F*[fnCuts];
      fHistoTrueClusSecPi0FromK0lPtvsM02          = new TH2F*[fnCuts];
      fHistoTrueClusSecPi0FromLambdaPtvsM02       = new TH2F*[fnCuts];
      fHistoTrueSecPi0PtvsDiffReco                = new TH2F*[fnCuts];
    }
    fHistoTrueClusEtaPtvsM02                      = new TH2F*[fnCuts];
    fHistoTrueClusEtaDalitzPtvsM02                = new TH2F*[fnCuts];
    fHistoTrueClusMergedPureFromPi0PtvsM02        = new TH2F*[fnCuts];
    fHistoTrueClusMergedPureFromEtaPtvsM02        = new TH2F*[fnCuts];
    fHistoTrueClusMergedPartConvFromPi0PtvsM02    = new TH2F*[fnCuts];
    fHistoTrueClusMergedPartConvFromEtaPtvsM02    = new TH2F*[fnCuts];
    fHistoNTrueMultiplePi0vsPt                    = new TH2F*[fnCuts];
    fHistoTrueClusBGPtvsM02                       = new TH2F*[fnCuts];
    fHistoTrueClusGammaPtvsM02                    = new TH2F*[fnCuts];
    fHistoTrueClusGammaFromPi0PtvsM02             = new TH2F*[fnCuts];
    fHistoTrueClusGammaFromEtaPtvsM02             = new TH2F*[fnCuts];
    fHistoTrueClusElectronPtvsM02                 = new TH2F*[fnCuts];
    fHistoTrueClusElectronFromPi0PtvsM02          = new TH2F*[fnCuts];
    fHistoTrueClusElectronFromEtaPtvsM02          = new TH2F*[fnCuts];
    fHistoTrueClusElectronFromGammaPtvsM02        = new TH2F*[fnCuts];

    fHistoTrueClusBGPtvsSource                    = new TH2F*[fnCuts];
    fHistoTrueClusGammaPtvsSource                 = new TH2F*[fnCuts];
    fHistoTrueClusElectronPtvsSource              = new TH2F*[fnCuts];
//    fHistoTrueClusElectronPtvsMotherID            = new TH2F*[fnCuts];
//    fHistoTrueClusElectronPtvsTopMotherID         = new TH2F*[fnCuts];
//    fHistoTrueClusElectronPtvsConvPhotonTopMotherID = new TH2F*[fnCuts];
    fHistoTrueMergedMissedPDG                     = new TH1F*[fnCuts];

    fHistoDoubleCountTruePi0PtvsM02               = new TH2F*[fnCuts];
    fHistoDoubleCountTrueSecPi0Pt                 = new TH1F*[fnCuts];
    fHistoDoubleCountTrueMultiplePi0PtvsM02       = new TH2F*[fnCuts];
    fHistoDoubleCountTrueMultipleSecPi0Pt         = new TH1F*[fnCuts];
    fHistoDoubleCountTrueEtaPtvsM02               = new TH2F*[fnCuts];

    if (fDoMesonQA > 1){
      fHistoTrueClusMergedInvMassvsPt               = new TH2F*[fnCuts];
      fHistoTrueClusPi0InvMassvsPt                  = new TH2F*[fnCuts];
      fHistoTrueClusEtaInvMassvsPt                  = new TH2F*[fnCuts];
      fHistoTrueClusBGInvMassvsPt                   = new TH2F*[fnCuts];
      fHistoTrueClusGammaInvMassvsPt                = new TH2F*[fnCuts];
      fHistoTrueClusElectronInvMassvsPt             = new TH2F*[fnCuts];
      if (GetSelectedMesonID() < 2){
        fHistoTrueClusPrimPi0InvMassvsPt                  = new TH2F*[fnCuts];
      }
      fHistoTrueClusGammaEvsM20                     = new TH2F*[fnCuts];
      fHistoTrueClusBGEvsM02                        = new TH2F*[fnCuts];
      fHistoTrueClusBGEvsM20                        = new TH2F*[fnCuts];
      fHistoTrueClusElectronEM20                    = new TH2F*[fnCuts];
      fHistoTrueClusMergedPi0EVsM20                 = new TH2F*[fnCuts];
      fHistoTrueClusMergedPi0EVsM02                 = new TH2F*[fnCuts];
      fHistoTrueClusMergedEtaEVsM20                 = new TH2F*[fnCuts];
      fHistoTrueClusMergedEtaEVsM02                 = new TH2F*[fnCuts];
    }

    if (fDoMesonQA > 0){
      fHistoTrueClusGammaEM02                     = new TH2F*[fnCuts];
      fHistoTrueClusElectronEM02                     = new TH2F*[fnCuts];
      fHistoTrueClusPi0EM02                       = new TH2F*[fnCuts];
      fHistoTrueClusEtaEM02                       = new TH2F*[fnCuts];

      if (GetSelectedMesonID() != 2){
        fHistoTruePi0PtY                          = new TH2F*[fnCuts];
        fHistoTruePi0PtAlpha                      = new TH2F*[fnCuts];
        fHistoTruePrimaryPi0PureMergedMCPtResolPt     = new TH2F*[fnCuts];
        fHistoTruePrimaryPi0MCPtResolPt               = new TH2F*[fnCuts];
        fHistoTruePrimaryPi0MergedPartConvMCPtResolPt = new TH2F*[fnCuts];
        fHistoTruePrimaryPi01GammaMCPtResolPt         = new TH2F*[fnCuts];
        fHistoTruePrimaryPi01ElectronMCPtResolPt      = new TH2F*[fnCuts];
        fHistoTrueSecondaryPi0MCPtResolPt             = new TH2F*[fnCuts];
      }
      if (GetSelectedMesonID() != 1){
        fHistoTrueEtaPtY                          = new TH2F*[fnCuts];
        fHistoTrueEtaPtAlpha                      = new TH2F*[fnCuts];
        fHistoTruePrimaryEtaMCPtResolPt           = new TH2F*[fnCuts];
      }
    }
    if(fDoJetAnalysis){
      if(!fDoLightOutput){
        fHistoTrueClusMergedPtvsRJet                 = new TH2F*[fnCuts];
        fHistoTrueClusPi0PtvsRJet                    = new TH2F*[fnCuts];
        fHistoTrueClusEtaPtvsRJet                    = new TH2F*[fnCuts];
      }
    }



    for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      TString cutstringEvent                        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringCalo                         = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
      TString cutstringCaloMerged                   = ((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson                        = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();


      fMCList[iCut]                                 = new TList();
      fMCList[iCut]->SetName(Form("%s_%s_%s_%s MC histograms",cutstringEvent.Data(),cutstringCalo.Data(),cutstringCaloMerged.Data(), cutstringMeson.Data()));
      fMCList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fMCList[iCut]);
      fHistoMCHeaders[iCut]             = new TH1I("MC_Headers", "MC_Headers", 20, 0, 20);
      fMCList[iCut]->Add(fHistoMCHeaders[iCut]);
      fHistoMCPi0Pt[iCut]                         = new TH1F("MC_Pi0_Pt","MC_Pi0_Pt",ptBins, arrPtBinning);
      fHistoMCPi0Pt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCPi0Pt[iCut]);
      fHistoMCPi0ReducedPt[iCut]                         = new TH1F("MC_Pi0MultipleSubtraction_Pt","MC_Pi0MultipleSubtraction_Pt",ptBins, arrPtBinning);
      fHistoMCPi0ReducedPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCPi0ReducedPt[iCut]);
      fHistoMCEtaPt[iCut]                         = new TH1F("MC_Eta_Pt","MC_Eta_Pt",ptBins, arrPtBinning);
      fHistoMCEtaPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCEtaPt[iCut]);

      fHistoMCPrimaryPtvsSource[iCut]   = new TH2F("MC_Primary_Pt_Source","MC_Primary_Pt_Source",ptBins, arrPtBinning, 7, -0.5, 6.5);
      fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(1,"Pi+");
      fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(2,"Pi-");
      fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(3,"K+");
      fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(4,"K-");
      fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(5,"K0s");
      fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(6,"K0l");
      fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(7,"Lambda");
      fHistoMCPrimaryPtvsSource[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCPrimaryPtvsSource[iCut]);

      fHistoMCPrimaryYvsSource[iCut]   = new TH2F("MC_Primary_Y_Source","MC_Primary_Y_Source",400, -20, 20, 7, -0.5, 6.5);
      fHistoMCPrimaryYvsSource[iCut]->GetYaxis()->SetBinLabel(1,"Pi+");
      fHistoMCPrimaryYvsSource[iCut]->GetYaxis()->SetBinLabel(2,"Pi-");
      fHistoMCPrimaryYvsSource[iCut]->GetYaxis()->SetBinLabel(3,"K+");
      fHistoMCPrimaryYvsSource[iCut]->GetYaxis()->SetBinLabel(4,"K-");
      fHistoMCPrimaryYvsSource[iCut]->GetYaxis()->SetBinLabel(5,"K0s");
      fHistoMCPrimaryYvsSource[iCut]->GetYaxis()->SetBinLabel(6,"K0l");
      fHistoMCPrimaryYvsSource[iCut]->GetYaxis()->SetBinLabel(7,"Lambda");
      fHistoMCPrimaryYvsSource[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCPrimaryYvsSource[iCut]);

      fHistoMCDecayGammaPt[iCut]                  = new TH1F("MC_DecayGamma_Pt","MC_DecayGamma_Pt",ptBins, arrPtBinning);
      fHistoMCDecayGammaPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCDecayGammaPt[iCut]);
      fHistoMCAllGammaPt[iCut]                    = new TH1F("MC_AllGamma_Pt","MC_AllGamma_Pt",ptBins, arrPtBinning);
      fHistoMCAllGammaPt[iCut]->Sumw2();
      fMCList[iCut]->Add(fHistoMCAllGammaPt[iCut]);

      if (GetSelectedMesonID() != 2){
        fHistoMCPi0WOWeightPt[iCut]                 = new TH1F("MC_Pi0_WOWeights_Pt","MC_Pi0_WOWeights_Pt",ptBins, arrPtBinning);
        fMCList[iCut]->Add(fHistoMCPi0WOWeightPt[iCut]);
        fHistoMCPi0InAccPt[iCut]                    = new TH1F("MC_Pi0InAcc_Pt","MC_Pi0InAcc_Pt",ptBins, arrPtBinning);
        fHistoMCPi0InAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0InAccPt[iCut]);
        fHistoMCPi0ReducedInAccPt[iCut]                    = new TH1F("MC_Pi0InAccMulipleSubtraction_Pt","MC_Pi0InAccMulipleSubtraction_Pt",ptBins, arrPtBinning);
        fHistoMCPi0ReducedInAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0ReducedInAccPt[iCut]);
        fHistoMCPi0DalitzPt[iCut]                   = new TH1F("MC_Pi0Dalitz_Pt","MC_Pi0Dalitz_Pt",ptBins, arrPtBinning);
        fHistoMCPi0DalitzPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0DalitzPt[iCut]);
        fHistoMCPi0DalitzWOWeightPt[iCut]           = new TH1F("MC_Pi0Dalitz_WOWeights_Pt","MC_Pi0Dalitz_WOWeights_Pt",ptBins, arrPtBinning);
        fMCList[iCut]->Add(fHistoMCPi0DalitzWOWeightPt[iCut]);
        fHistoMCPi0DalitzInAccPt[iCut]              = new TH1F("MC_Pi0DalitzInAcc_Pt","MC_Pi0DalitzInAcc_Pt",ptBins, arrPtBinning);
        fHistoMCPi0DalitzInAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0DalitzInAccPt[iCut]);

        // initialization for secondary histograms
        fHistoMCSecPi0PtvsSource[iCut]              = new TH2F("MC_SecPi0_Pt_Source","MC_SecPi0_Pt_Source", ptBins, arrPtBinning, 16, -0.5, 15.5);
        fHistoMCSecPi0PtvsSource[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCSecPi0PtvsSource[iCut]);
        fHistoMCSecPi0InAccPtvsSource[iCut]         = new TH2F("MC_SecPi0InAcc_Pt_Source","MC_SecPi0InAcc_Pt_Source", ptBins, arrPtBinning, 16, -0.5, 15.5);
        fHistoMCSecPi0InAccPtvsSource[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCSecPi0InAccPtvsSource[iCut]);

        fHistoPi0EvsGammaOverlapE[iCut]           = new TH2F("MC_Pi0_E_OverlapGamma_E","MC_Pi0_E_OverlapGamma_E",200,0,200,200,0,200);
        fMCList[iCut]->Add(fHistoPi0EvsGammaOverlapE[iCut]);
      }

      if (GetSelectedMesonID() != 1){
        fHistoMCEtaWOWeightPt[iCut]                 = new TH1F("MC_Eta_WOWeights_Pt","MC_Eta_WOWeights_Pt",ptBins, arrPtBinning);
        fMCList[iCut]->Add(fHistoMCEtaWOWeightPt[iCut]);
        fHistoMCEtaInAccPt[iCut]                    = new TH1F("MC_EtaInAcc_Pt","MC_EtaInAcc_Pt",ptBins, arrPtBinning);
        fHistoMCEtaInAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaInAccPt[iCut]);
        fHistoMCEtaDalitzPt[iCut]                   = new TH1F("MC_EtaDalitz_Pt","MC_EtaDalitz_Pt",ptBins, arrPtBinning);
        fHistoMCEtaDalitzPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaDalitzPt[iCut]);
        fHistoMCEtaDalitzWOWeightPt[iCut]           = new TH1F("MC_EtaDalitz_WOWeights_Pt","MC_EtaDalitz_WOWeights_Pt",ptBins, arrPtBinning);
        fMCList[iCut]->Add(fHistoMCEtaDalitzWOWeightPt[iCut]);
        fHistoMCEtaDalitzInAccPt[iCut]              = new TH1F("MC_EtaDalitzInAcc_Pt","MC_EtaDalitzInAcc_Pt",ptBins, arrPtBinning);
        fHistoMCEtaDalitzInAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCEtaDalitzInAccPt[iCut]);
      }

      if (fIsMC == 2){
        if (GetSelectedMesonID() != 2){
          fHistoMCPi0WOWeightPt[iCut]->Sumw2();
          fHistoMCPi0DalitzWOWeightPt[iCut]->Sumw2();
          fHistoMCPi0WOEvtWeightPt[iCut]            = new TH1F("MC_Pi0_WOEventWeights_Pt","MC_Pi0_WOEventWeights_Pt",ptBins, arrPtBinning);
          fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightPt[iCut]);
          fHistoMCPi0DalitzWOEvtWeightPt[iCut]      = new TH1F("MC_Pi0Dalitz_WOEventWeights_Pt","MC_Pi0Dalitz_WOEventWeights_Pt",ptBins, arrPtBinning);
          fMCList[iCut]->Add(fHistoMCPi0DalitzWOEvtWeightPt[iCut]);
          fHistoMCPi0DalitzWOEvtWeightInAccPt[iCut]   = new TH1F("MC_Pi0DalitzWOEvtWeightInAcc_Pt","MC_Pi0DalitzWOEvtWeightInAcc_Pt",ptBins, arrPtBinning);
          fMCList[iCut]->Add(fHistoMCPi0DalitzWOEvtWeightInAccPt[iCut]);
          fHistoMCPi0WOEvtWeightInAccPt[iCut]         = new TH1F("MC_Pi0WOEvtWeightInAcc_Pt","MC_Pi0WOEvtWeightInAcc_Pt",ptBins, arrPtBinning);
          fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightInAccPt[iCut]);


        }
        if (GetSelectedMesonID() != 1){
          fHistoMCEtaWOWeightPt[iCut]->Sumw2();
          fHistoMCEtaDalitzWOWeightPt[iCut]->Sumw2();
          fHistoMCEtaWOEvtWeightPt[iCut]            = new TH1F("MC_Eta_WOEventWeights_Pt","MC_Eta_WOEventWeights_Pt",ptBins, arrPtBinning);
          fMCList[iCut]->Add(fHistoMCEtaWOEvtWeightPt[iCut]);
          fHistoMCEtaDalitzWOEvtWeightPt[iCut]      = new TH1F("MC_EtaDalitz_WOEventWeights_Pt","MC_EtaDalitz_WOEventWeights_Pt",ptBins, arrPtBinning);
          fMCList[iCut]->Add(fHistoMCEtaDalitzWOEvtWeightPt[iCut]);
          fHistoMCEtaWOEvtWeightInAccPt[iCut]         = new TH1F("MC_EtaWOEvtWeightInAcc_Pt","MC_EtaWOEvtWeightInAcc_Pt",ptBins, arrPtBinning);
          fMCList[iCut]->Add(fHistoMCEtaWOEvtWeightInAccPt[iCut]);
          fHistoMCEtaDalitzWOEvtWeightInAccPt[iCut]   = new TH1F("MC_EtaDalitzWOEvtWeightInAcc_Pt","MC_EtaDalitzWOEvtWeightInAcc_Pt",ptBins, arrPtBinning);
          fMCList[iCut]->Add(fHistoMCEtaDalitzWOEvtWeightInAccPt[iCut]);
        }
        if (fDoMesonQA > 0){
          if (GetSelectedMesonID() != 2){
            fHistoMCPi0PtJetPt[iCut]                = new TH2F("MC_Pi0_Pt_JetPt","MC_Pi0_Pt_JetPt", ptBinsLog, startPtLog, endPtLog, 200, -0.5, 199.5);
            fHistoMCPi0PtJetPt[iCut]->Sumw2();
            SetLogBinningXTH2(fHistoMCPi0PtJetPt[iCut]);
            fMCList[iCut]->Add(fHistoMCPi0PtJetPt[iCut]);
          }
          if (GetSelectedMesonID() != 1){
            fHistoMCEtaPtJetPt[iCut]                = new TH2F("MC_Eta_Pt_JetPt","MC_Eta_Pt_JetPt", ptBinsLog, startPtLog, endPtLog, 200, -0.5, 199.5);
            fHistoMCEtaPtJetPt[iCut]->Sumw2();
            SetLogBinningXTH2(fHistoMCEtaPtJetPt[iCut]);
            fMCList[iCut]->Add(fHistoMCEtaPtJetPt[iCut]);
          }
        }
      }

      fTrueList[iCut]                               = new TList();
      fTrueList[iCut]->SetName(Form("%s_%s_%s_%s True histograms",cutstringEvent.Data(),cutstringCalo.Data(),cutstringCaloMerged.Data(), cutstringMeson.Data()));
      fTrueList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fTrueList[iCut]);

      fHistoTrueClusEFracFirstLabel[iCut]             = new TH2F("ESD_TrueClusEFracFirstLabel_E_Frac","ESD_TrueClusEFracFirstLabel_E_Frac",ptBins, arrPtBinning,50, -0.01, 1.01);
      fTrueList[iCut]->Add(fHistoTrueClusEFracFirstLabel[iCut]);
      fHistoTrueClusEFracLeadingPi0[iCut]             = new TH2F("ESD_TrueClusEFracLeadingPi0_E_Frac","ESD_TrueClusEFracLeadingPi0_E_Frac",ptBins, arrPtBinning,50, -0.01, 1.01);
      fTrueList[iCut]->Add(fHistoTrueClusEFracLeadingPi0[iCut]);
      fHistoTrueClusMergedPtvsM02[iCut]             = new TH2F("ESD_TrueClusMerged_Pt_M02","ESD_TrueClusMerged_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusMergedPtvsM02[iCut]);
      fHistoTrueClusPi0PtvsM02 [iCut]               = new TH2F("ESD_TrueClusFromPi0_Pt_M02","ESD_TrueClusFromPi0_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusPi0PtvsM02[iCut]);
      fHistoTrueClusMultiplePi0PtvsM02 [iCut]       = new TH2F("ESD_TrueClusFromMultiplePi0_Pt_M02","ESD_TrueClusFromMultiplePi0_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusMultiplePi0PtvsM02[iCut]);
      fHistoTrueClusPi0DalitzPtvsM02 [iCut]         = new TH2F("ESD_TrueClusFromPi0Dalitz_Pt_M02","ESD_TrueClusFromPi0Dalitz_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusPi0DalitzPtvsM02[iCut]);
      fHistoTrueClusEtaPtvsM02[iCut]                = new TH2F("ESD_TrueClusFromEta_Pt_M02","ESD_TrueClusFromEta_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusEtaPtvsM02[iCut]);
      fHistoTrueClusEtaDalitzPtvsM02[iCut]          = new TH2F("ESD_TrueClusFromEtaDalitz_Pt_M02","ESD_TrueClusFromEtaDalitz_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusEtaDalitzPtvsM02[iCut]);
      fHistoTrueClusMergedPureFromPi0PtvsM02[iCut]     = new TH2F("ESD_TrueClusMergedPureFromPi0_Pt_M02","ESD_TrueClusMergedPureFromPi0_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusMergedPureFromPi0PtvsM02[iCut]);
      fHistoTrueClusMergedPureFromEtaPtvsM02[iCut]     = new TH2F("ESD_TrueClusMergedPureFromEta_Pt_M02","ESD_TrueClusMergedPureFromEta_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusMergedPureFromEtaPtvsM02[iCut]);
      fHistoTrueClusMergedPartConvFromPi0PtvsM02[iCut]     = new TH2F("ESD_TrueClusMergedPartConvFromPi0_Pt_M02","ESD_TrueClusMergedPartConvFromPi0_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusMergedPartConvFromPi0PtvsM02[iCut]);
      fHistoTrueClusMergedPartConvFromEtaPtvsM02[iCut]     = new TH2F("ESD_TrueClusMergedPartConvFromEta_Pt_M02","ESD_TrueClusMergedPartConvFromEta_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusMergedPartConvFromEtaPtvsM02[iCut]);
      fHistoNTrueMultiplePi0vsPt [iCut]       = new TH2F("ESD_TrueNumberOfMultiplePi0_Pt","ESD_TrueNumberOfMultiplePi0_Pt",ptBins, arrPtBinning,20, 0, 20);
      fTrueList[iCut]->Add(fHistoNTrueMultiplePi0vsPt[iCut]);

      fHistoTrueClusBGPtvsM02[iCut]                 = new TH2F("ESD_TrueClusBG_Pt_M02","ESD_TrueClusBG_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusBGPtvsM02[iCut]);
      fHistoTrueClusGammaPtvsM02[iCut]              = new TH2F("ESD_TrueClusGamma_Pt_M02","ESD_TrueClusGamma_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusGammaPtvsM02[iCut]);
      fHistoTrueClusGammaFromPi0PtvsM02[iCut]       = new TH2F("ESD_TrueClusGamma_FromPi0_Pt_M02","ESD_TrueClusGamma_FromPi0_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusGammaFromPi0PtvsM02[iCut]);
      fHistoTrueClusGammaFromEtaPtvsM02[iCut]       = new TH2F("ESD_TrueClusGamma_FromEta_Pt_M02","ESD_TrueClusGamma_FromEta_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusGammaFromEtaPtvsM02[iCut]);
      fHistoTrueClusElectronPtvsM02[iCut]           = new TH2F("ESD_TrueClusElectron_Pt_M02","ESD_TrueClusElectron_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusElectronPtvsM02[iCut]);
      fHistoTrueClusElectronFromPi0PtvsM02[iCut]    = new TH2F("ESD_TrueClusElectron_FromPi0_Pt_M02","ESD_TrueClusElectron_FromPi0_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusElectronFromPi0PtvsM02[iCut]);
      fHistoTrueClusElectronFromEtaPtvsM02[iCut]    = new TH2F("ESD_TrueClusElectron_FromEta_Pt_M02","ESD_TrueClusElectron_FromEta_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusElectronFromEtaPtvsM02[iCut]);
      fHistoTrueClusElectronFromGammaPtvsM02[iCut]  = new TH2F("ESD_TrueClusElectron_FromGamma_Pt_M02","ESD_TrueClusElectron_FromGamma_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoTrueClusElectronFromGammaPtvsM02[iCut]);

      if (GetSelectedMesonID() < 2){
        fHistoTrueClusPrimPi0PtvsM02[iCut]                  = new TH2F("ESD_TrueClusFromPrimPi0_Pt_M02","ESD_TrueClusFromPrimPi0_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusPrimPi0PtvsM02[iCut]);
        fHistoTrueClusPrimPi0PtMCPt[iCut]                  = new TH2F("ESD_TrueClusFromPrimPi0_Pt_MCPt","ESD_TrueClusFromPrimPi0_Pt_MCPt",ptBins, arrPtBinning,ptBins, arrPtBinning);
        fTrueList[iCut]->Add(fHistoTrueClusPrimPi0PtMCPt[iCut]);
        fHistoTruePureMergedClusPrimPi0PtvsM02[iCut]      = new TH2F("ESD_TruePureMergedClusFromPrimPi0_Pt_M02","ESD_TruePureMergedClusFromPrimPi0_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTruePureMergedClusPrimPi0PtvsM02[iCut]);
        fHistoTruePartConvClusPrimPi0PtvsM02[iCut]        = new TH2F("ESD_TruePartConvClusFromPrimPi0_Pt_M02","ESD_TruePartConvClusFromPrimPi0_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTruePartConvClusPrimPi0PtvsM02[iCut]);
        if (fDoMesonQA > 1){
          fHistoTrueClusMergedPartConvPi0EVsM20[iCut]        = new TH2F("ESD_TruePartConvClusFromPrimPi0_E_M20","ESD_TruePartConvClusFromPrimPi0_E_M20",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
          fTrueList[iCut]->Add(fHistoTrueClusMergedPartConvPi0EVsM20[iCut]);
          fHistoTrueClusPrimPi0PureMergedPtMCPt[iCut]        = new TH2F("ESD_TrueClusFromPrimPi0PureMerged_Pt_MCPt","ESD_TrueClusFromPrimPi0PureMerged_Pt_MCPt",ptBins, arrPtBinning,ptBins, arrPtBinning);
          fTrueList[iCut]->Add(fHistoTrueClusPrimPi0PureMergedPtMCPt[iCut]);
          fHistoTrueClusPrimPi0PartConvMergedPtMCPt[iCut]    = new TH2F("ESD_TrueClusFromPrimPi0PartConvMerged_Pt_MCPt","ESD_TrueClusFromPrimPi0PartConvMerged_Pt_MCPt",ptBins, arrPtBinning,ptBins, arrPtBinning);
          fTrueList[iCut]->Add(fHistoTrueClusPrimPi0PartConvMergedPtMCPt[iCut]);
          fHistoTrueClusPrimPi01GammaMergedPtMCPt[iCut]      = new TH2F("ESD_TrueClusFromPrimPi01GammaMerged_Pt_MCPt","ESD_TrueClusFromPrimPi01GammaMerged_Pt_MCPt",ptBins, arrPtBinning,ptBins, arrPtBinning);
          fTrueList[iCut]->Add(fHistoTrueClusPrimPi01GammaMergedPtMCPt[iCut]);
          fHistoTrueClusPrimPi01ElectronMergedPtMCPt[iCut]   = new TH2F("ESD_TrueClusFromPrimPi01ElectronMerged_Pt_MCPt","ESD_TrueClusFromPrimPi01ElectronMerged_Pt_MCPt",ptBins, arrPtBinning,ptBins, arrPtBinning);
          fTrueList[iCut]->Add(fHistoTrueClusPrimPi01ElectronMergedPtMCPt[iCut]);

          fHistoTrueClusNeutralContamination[iCut]             = new TH2F("ESD_TrueClusE_RelContaminationNeutralPart","ESD_TrueClusE_RelContaminationNeutralPart",ptBins, arrPtBinning,100, -1.0, 1.0);
          fTrueList[iCut]->Add(fHistoTrueClusNeutralContamination[iCut]);
          fHistoTrueClusPhotonContamination[iCut]             = new TH2F("ESD_TrueClusE_RelContaminationPhotons","ESD_TrueClusE_RelContaminationPhotons",ptBins, arrPtBinning,100, -1.0, 1.0);
          fTrueList[iCut]->Add(fHistoTrueClusPhotonContamination[iCut]);
          fHistoTrueClusElectronContamination[iCut]             = new TH2F("ESD_TrueClusE_RelContaminationElectrons","ESD_TrueClusE_RelContaminationElectrons",ptBins, arrPtBinning,100, -1.0, 1.0);
          fTrueList[iCut]->Add(fHistoTrueClusElectronContamination[iCut]);
          fHistoTrueClusOtherContamination[iCut]             = new TH2F("ESD_TrueClusE_RelContaminationOthers","ESD_TrueClusE_RelContaminationOthers",ptBins, arrPtBinning,100, -1.0, 1.0);
          fTrueList[iCut]->Add(fHistoTrueClusOtherContamination[iCut]);
        }
        fHistoTrueClusMultiplePrimPi0PtvsM02[iCut]          = new TH2F("ESD_TrueClusFromMultiplePrimPi0_Pt_M02","ESD_TrueClusFromMultiplePrimPi0_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusMultiplePrimPi0PtvsM02[iCut]);
        fHistoTrueClusSecPi0PtvsM02[iCut]                   = new TH2F("ESD_TrueClusFromSecPi0_Pt_M02","ESD_TrueClusFromSecPi0_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusSecPi0PtvsM02[iCut]);
        fHistoTrueClusSecPi0FromK0sPtvsM02[iCut]            = new TH2F("ESD_TrueClusFromSecPi0FromK0s_Pt_M02","ESD_TrueClusFromSecPi0FromK0s_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusSecPi0FromK0sPtvsM02[iCut]);
        fHistoTrueClusSecPi0FromK0lPtvsM02[iCut]            = new TH2F("ESD_TrueClusFromSecPi0FromK0l_Pt_M02","ESD_TrueClusFromSecPi0FromK0l_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusSecPi0FromK0lPtvsM02[iCut]);
        fHistoTrueClusSecPi0FromLambdaPtvsM02[iCut]         = new TH2F("ESD_TrueClusFromSecPi0FromLambda_Pt_M02","ESD_TrueClusFromSecPi0FromLambda_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusSecPi0FromLambdaPtvsM02[iCut]);
        fHistoTrueSecPi0PtvsDiffReco[iCut]                  = new TH2F("ESD_TrueClusFromSecPi0_Pt_RecoMethod","ESD_TrueClusFromSecPi0_Pt_RecoMethod",ptBins, arrPtBinning, 4, -0.5, 3.5);
        fTrueList[iCut]->Add(fHistoTrueSecPi0PtvsDiffReco[iCut]);
      }


      fHistoTrueClusBGPtvsSource[iCut]                      = new TH2F("ESD_TrueClusBG_Pt_Source","ESD_TrueClusBG_Pt_Source",ptBins, arrPtBinning, 10, 0, 10);
      fTrueList[iCut]->Add(fHistoTrueClusBGPtvsSource[iCut]);
      fHistoTrueClusGammaPtvsSource[iCut]                   = new TH2F("ESD_TrueClusGamma_Pt_Source","ESD_TrueClusGamma_Pt_Source",ptBins, arrPtBinning, 8, 0, 8);
      fTrueList[iCut]->Add(fHistoTrueClusGammaPtvsSource[iCut]);
      fHistoTrueClusElectronPtvsSource[iCut]                = new TH2F("ESD_TrueClusElectron_Pt_Source","ESD_TrueClusElectron_Pt_Source",ptBins, arrPtBinning, 9, 0, 9);
      fTrueList[iCut]->Add(fHistoTrueClusElectronPtvsSource[iCut]);
//      fHistoTrueClusElectronPtvsMotherID[iCut]              = new TH2F("ESD_TrueClusElectron_Pt_MotherID","ESD_TrueClusElectron_Pt_MotherID",ptBins, arrPtBinning, 570, -1.5, 568.5);
//      fTrueList[iCut]->Add(fHistoTrueClusElectronPtvsMotherID[iCut]);
//      fHistoTrueClusElectronPtvsTopMotherID[iCut]           = new TH2F("ESD_TrueClusElectron_Pt_TopMotherID","ESD_TrueClusElectron_Pt_TopMotherID",ptBins, arrPtBinning, 570, -1.5, 568.5);
//      fTrueList[iCut]->Add(fHistoTrueClusElectronPtvsTopMotherID[iCut]);
//      fHistoTrueClusElectronPtvsConvPhotonTopMotherID[iCut] = new TH2F("ESD_TrueClusElectron_Pt_ConvPhotonTopMotherID","ESD_TrueClusElectron_Pt_ConvPhotonTopMotherID",ptBins, arrPtBinning, 570, -1.5, 568.5);
//      fTrueList[iCut]->Add(fHistoTrueClusElectronPtvsConvPhotonTopMotherID[iCut]);
      fHistoTrueMergedMissedPDG[iCut]                       = new TH1F("ESD_TrueMergedMissed_PDG","ESD_TrueMergedMissed_PDG",10000, -1.5, 9998.5);
      fTrueList[iCut]->Add(fHistoTrueMergedMissedPDG[iCut]);

      fHistoDoubleCountTruePi0PtvsM02[iCut]         = new TH2F("ESD_TrueDoubleCountPi0_Pt_M02","ESD_TrueDoubleCountPi0_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoDoubleCountTruePi0PtvsM02[iCut]);
      fHistoDoubleCountTrueSecPi0Pt[iCut]           = new TH1F("ESD_TrueDoubleCountSecPi0_Pt","ESD_TrueDoubleCountSecPi0_Pt",ptBins, arrPtBinning);
      fTrueList[iCut]->Add(fHistoDoubleCountTrueSecPi0Pt[iCut]);
      fHistoDoubleCountTrueMultiplePi0PtvsM02[iCut]         = new TH2F("ESD_TrueDoubleCountMultiplePi0_Pt_M02","ESD_TrueDoubleCountMultiplePi0_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoDoubleCountTrueMultiplePi0PtvsM02[iCut]);
      fHistoDoubleCountTrueMultipleSecPi0Pt[iCut]           = new TH1F("ESD_TrueDoubleCountMultipleSecPi0_Pt","ESD_TrueDoubleCountMultipleSecPi0_Pt",ptBins, arrPtBinning);
      fTrueList[iCut]->Add(fHistoDoubleCountTrueMultipleSecPi0Pt[iCut]);
      fHistoDoubleCountTrueEtaPtvsM02[iCut]         = new TH2F("ESD_TrueDoubleCountEta_Pt_M02","ESD_TrueDoubleCountEta_Pt_M02",ptBins, arrPtBinning,showerShapeBins, startShowerShape, endShowerShape);
      fTrueList[iCut]->Add(fHistoDoubleCountTrueEtaPtvsM02[iCut]);


      if (fDoMesonQA > 1){
        fHistoTrueClusMergedInvMassvsPt[iCut]                 = new TH2F("ESD_TrueClusMerged_InvMass_Pt","ESD_TrueClusMerged_InvMass_Pt",invMassBins, startMass, endMass,ptBins, arrPtBinning);
        fTrueList[iCut]->Add(fHistoTrueClusMergedInvMassvsPt[iCut]);
        fHistoTrueClusPi0InvMassvsPt [iCut]                   = new TH2F("ESD_TrueClusFromPi0_InvMass_Pt","ESD_TrueClusFromPi0_InvMass_Pt",invMassBins, startMass, endMass,ptBins, arrPtBinning);
        fTrueList[iCut]->Add(fHistoTrueClusPi0InvMassvsPt[iCut]);
        fHistoTrueClusEtaInvMassvsPt[iCut]                    = new TH2F("ESD_TrueClusFromEta_InvMass_Pt","ESD_TrueClusFromEta_InvMass_Pt",invMassBins, startMass, endMass,ptBins, arrPtBinning);
        fTrueList[iCut]->Add(fHistoTrueClusEtaInvMassvsPt[iCut]);
        fHistoTrueClusBGInvMassvsPt[iCut]                     = new TH2F("ESD_TrueClusBG_InvMass_Pt","ESD_TrueClusBG_InvMass_Pt",invMassBins, startMass, endMass,ptBins, arrPtBinning);
        fTrueList[iCut]->Add(fHistoTrueClusBGInvMassvsPt[iCut]);
        fHistoTrueClusGammaInvMassvsPt[iCut]                  = new TH2F("ESD_TrueClusGamma_InvMass_Pt","ESD_TrueClusGamma_InvMass_Pt",invMassBins, startMass, endMass,ptBins, arrPtBinning);
        fTrueList[iCut]->Add(fHistoTrueClusGammaInvMassvsPt[iCut]);
        fHistoTrueClusElectronInvMassvsPt[iCut]               = new TH2F("ESD_TrueClusElectron_InvMass_Pt","ESD_TrueClusElectron_InvMass_Pt",invMassBins, startMass, endMass,ptBins, arrPtBinning);
        fTrueList[iCut]->Add(fHistoTrueClusElectronInvMassvsPt[iCut]);
        if (GetSelectedMesonID() < 2){
          fHistoTrueClusPrimPi0InvMassvsPt[iCut]                    = new TH2F("ESD_TrueClusFromPrimPi0_InvMass_Pt","ESD_TrueClusFromPrimPi0_InvMass_Pt",invMassBins, startMass, endMass,ptBins, arrPtBinning);
          fTrueList[iCut]->Add(fHistoTrueClusPrimPi0InvMassvsPt[iCut]);
        }
        fHistoTrueClusGammaEvsM20[iCut]                       = new TH2F("ESD_TrueClusGammaEM20","ESD_TrueClusGammaEM20",ptBins, arrPtBinning, showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusGammaEvsM20[iCut]);
        fHistoTrueClusBGEvsM02[iCut]                          = new TH2F("ESD_TrueClusBGEvsM02","ESD_TrueClusBGEvsM02",ptBins, arrPtBinning, showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusBGEvsM02[iCut]);
        fHistoTrueClusBGEvsM20[iCut]                          = new TH2F("ESD_TrueClusBGEvsM20","ESD_TrueClusBGEvsM20",ptBins, arrPtBinning, showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusBGEvsM20[iCut]);
        fHistoTrueClusElectronEM20[iCut]                       = new TH2F("ESD_TrueClusElectronEM20","ESD_TrueClusElectronEM20",ptBins, arrPtBinning, showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusElectronEM20[iCut]);
        fHistoTrueClusMergedPi0EVsM20[iCut]                            = new TH2F("ESD_TrueMergedClusPi0EM20","ESD_TrueMergedClusPi0EM20",ptBins, arrPtBinning, showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusMergedPi0EVsM20[iCut]);
        fHistoTrueClusMergedPi0EVsM02[iCut]                            = new TH2F("ESD_TrueMergedClusPi0EM02","ESD_TrueMergedClusPi0EM02",ptBins, arrPtBinning, showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusMergedPi0EVsM02[iCut]);
        fHistoTrueClusMergedEtaEVsM20[iCut]                            = new TH2F("ESD_TrueMergedClusEtaEM20","ESD_TrueMergedClusEtaEM20",ptBins, arrPtBinning, showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusMergedEtaEVsM20[iCut]);
        fHistoTrueClusMergedEtaEVsM02[iCut]                            = new TH2F("ESD_TrueMergedClusEtaEM02","ESD_TrueMergedClusEtaEM02",ptBins, arrPtBinning, showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusMergedEtaEVsM02[iCut]);
      }

      if (fDoMesonQA > 0){
        fHistoTrueClusGammaEM02[iCut]                       = new TH2F("TrueClusGammaEM02","TrueClusGammaEM02",ptBins, arrPtBinning, showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusGammaEM02[iCut]);
        fHistoTrueClusElectronEM02[iCut]                       = new TH2F("TrueClusElectronEM02","TrueClusElectronEM02",ptBins, arrPtBinning, showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusElectronEM02[iCut]);
        fHistoTrueClusPi0EM02[iCut]                         = new TH2F("TrueClusPi0EM02","TrueClusPi0EM02",ptBins, arrPtBinning, showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusPi0EM02[iCut]);
        fHistoTrueClusEtaEM02[iCut]                         = new TH2F("TrueClusEtaEM02","TrueClusEtaEM02",ptBins, arrPtBinning, showerShapeBins, startShowerShape, endShowerShape);
        fTrueList[iCut]->Add(fHistoTrueClusEtaEM02[iCut]);
        if (GetSelectedMesonID() != 2){
          fHistoTruePi0PtY[iCut]                            = new TH2F("ESD_TruePi0_Pt_Y","ESD_TruePi0_Pt_Y",ptBinsLog, startPtLog, endPtLog,150,-1.5,1.5);
          SetLogBinningXTH2(fHistoTruePi0PtY[iCut]);
          fTrueList[iCut]->Add(fHistoTruePi0PtY[iCut]);
          fHistoTruePi0PtAlpha[iCut]                        = new TH2F("ESD_TruePi0_Pt_Alpha","ESD_TruePi0_Pt_Alpha",ptBinsLog, startPtLog, endPtLog,100,0,1);
          SetLogBinningXTH2(fHistoTruePi0PtAlpha[iCut]);
          fTrueList[iCut]->Add(fHistoTruePi0PtAlpha[iCut]);
          fHistoTruePrimaryPi0MCPtResolPt[iCut]       = new TH2F("ESD_TruePrimaryPi0_MCPt_ResolPt","ESD_TruePrimaryPi0_ResolPt_MCPt",ptBinsLog, startPtLog, endPtLog,1000,-1.,1.);
          SetLogBinningXTH2(fHistoTruePrimaryPi0MCPtResolPt[iCut]);
          fTrueList[iCut]->Add(fHistoTruePrimaryPi0MCPtResolPt[iCut]);
          fHistoTruePrimaryPi0PureMergedMCPtResolPt[iCut]       = new TH2F("ESD_TruePrimaryPi0PureMerged_MCPt_ResolPt","ESD_TruePrimaryPi0PureMerged_ResolPt_MCPt",ptBinsLog, startPtLog, endPtLog,1000,-1.,1.);
          SetLogBinningXTH2(fHistoTruePrimaryPi0PureMergedMCPtResolPt[iCut]);
          fTrueList[iCut]->Add(fHistoTruePrimaryPi0PureMergedMCPtResolPt[iCut]);
          fHistoTruePrimaryPi0MergedPartConvMCPtResolPt[iCut]   = new TH2F("ESD_TruePrimaryPi0MergedPartConv_MCPt_ResolPt","ESD_TruePrimaryPi0MergedPartConv_ResolPt_MCPt",ptBinsLog, startPtLog, endPtLog,1000,-1.,1.);
          SetLogBinningXTH2(fHistoTruePrimaryPi0MergedPartConvMCPtResolPt[iCut]);
          fTrueList[iCut]->Add(fHistoTruePrimaryPi0MergedPartConvMCPtResolPt[iCut]);
          fHistoTruePrimaryPi01GammaMCPtResolPt[iCut]           = new TH2F("ESD_TruePrimaryPi01Gamma_MCPt_ResolPt","ESD_TruePrimaryPi01Gamma_ResolPt_MCPt",ptBinsLog, startPtLog, endPtLog,1000,-1.,1.);
          SetLogBinningXTH2(fHistoTruePrimaryPi01GammaMCPtResolPt[iCut]);
          fTrueList[iCut]->Add(fHistoTruePrimaryPi01GammaMCPtResolPt[iCut]);
          fHistoTruePrimaryPi01ElectronMCPtResolPt[iCut]        = new TH2F("ESD_TruePrimaryPi01Electron_MCPt_ResolPt","ESD_TruePrimaryPi01Electron_ResolPt_MCPt",ptBinsLog, startPtLog, endPtLog,1000,-1.,1.);
          SetLogBinningXTH2(fHistoTruePrimaryPi01ElectronMCPtResolPt[iCut]);
          fTrueList[iCut]->Add(fHistoTruePrimaryPi01ElectronMCPtResolPt[iCut]);
          fHistoTrueSecondaryPi0MCPtResolPt[iCut]               = new TH2F("ESD_TrueSecondaryPi0_MCPt_ResolPt","ESD_TrueSecondaryPi0_MCPt_ResolPt",ptBinsLog, startPtLog, endPtLog,2000,-10.,10.);
          SetLogBinningXTH2(fHistoTrueSecondaryPi0MCPtResolPt[iCut]);
          fTrueList[iCut]->Add(fHistoTrueSecondaryPi0MCPtResolPt[iCut]);
        }
        if (GetSelectedMesonID() != 1){
          fHistoTrueEtaPtY[iCut]                            = new TH2F("ESD_TrueEta_Pt_Y","ESD_TrueEta_Pt_Y",ptBinsLog, startPtLog, endPtLog,150,-1.5,1.5);
          SetLogBinningXTH2(fHistoTrueEtaPtY[iCut]);
          fTrueList[iCut]->Add(fHistoTrueEtaPtY[iCut]);
          fHistoTrueEtaPtAlpha[iCut]                        = new TH2F("ESD_TrueEta_Pt_Alpha","ESD_TrueEta_Pt_Alpha",ptBinsLog, startPtLog, endPtLog,100,0,0.5);
          SetLogBinningXTH2(fHistoTrueEtaPtAlpha[iCut]);
          fTrueList[iCut]->Add(fHistoTrueEtaPtAlpha[iCut]);
          fHistoTruePrimaryEtaMCPtResolPt[iCut]             = new TH2F("ESD_TruePrimaryEta_MCPt_ResolPt","ESD_TruePrimaryEta_ResolPt_MCPt",ptBinsLog, startPtLog, endPtLog,1000,-1.,1.);
          SetLogBinningXTH2(fHistoTruePrimaryEtaMCPtResolPt[iCut]);
          fTrueList[iCut]->Add(fHistoTruePrimaryEtaMCPtResolPt[iCut]);
        }
      }

      if(fDoJetAnalysis){
        if(!fDoLightOutput){
          fHistoTrueClusMergedPtvsRJet[iCut]          = new TH2F("TrueClusPtVsR","TrueClusPtVsR",ptBins, arrPtBinning, 20, 0, 0.4);
          fTrueList[iCut]->Add(fHistoTrueClusMergedPtvsRJet[iCut]);
          fHistoTrueClusPi0PtvsRJet[iCut]             = new TH2F("TrueClusPi0PtVsR","TrueClusPi0PtVsR",ptBins, arrPtBinning, 20, 0, 0.4);
          fTrueList[iCut]->Add(fHistoTrueClusPi0PtvsRJet[iCut]);
          fHistoTrueClusEtaPtvsRJet[iCut]             = new TH2F("TrueClusEtaPtVsR","TrueClusEtaPtVsR",ptBins, arrPtBinning, 20, 0, 0.4);
          fTrueList[iCut]->Add(fHistoTrueClusEtaPtvsRJet[iCut]);
        }
      }

      if (fIsMC == 2){
        fHistoTrueClusEFracFirstLabel[iCut]->Sumw2();
        fHistoTrueClusEFracLeadingPi0[iCut]->Sumw2();
        fHistoTrueClusMergedPtvsM02[iCut]->Sumw2();
        fHistoTrueClusPi0PtvsM02[iCut]->Sumw2();
        fHistoTrueClusMultiplePi0PtvsM02[iCut]->Sumw2();
        fHistoNTrueMultiplePi0vsPt[iCut]->Sumw2();
        fHistoTrueClusPi0DalitzPtvsM02[iCut]->Sumw2();
        fHistoTrueClusEtaPtvsM02[iCut]->Sumw2();
        fHistoTrueClusEtaDalitzPtvsM02[iCut]->Sumw2();
        fHistoTrueClusMergedPureFromPi0PtvsM02[iCut]->Sumw2();
        fHistoTrueClusMergedPureFromEtaPtvsM02[iCut]->Sumw2();
        fHistoTrueClusMergedPartConvFromPi0PtvsM02[iCut]->Sumw2();
        fHistoTrueClusMergedPartConvFromEtaPtvsM02[iCut]->Sumw2();
        fHistoTrueMergedMissedPDG[iCut]->Sumw2();
        fHistoTrueClusBGPtvsM02[iCut]->Sumw2();
        fHistoTrueClusGammaPtvsM02[iCut]->Sumw2();
        fHistoTrueClusGammaFromPi0PtvsM02[iCut]->Sumw2();
        fHistoTrueClusGammaFromEtaPtvsM02[iCut]->Sumw2();
        fHistoTrueClusElectronPtvsM02[iCut]->Sumw2();
        fHistoTrueClusElectronFromPi0PtvsM02[iCut]->Sumw2();
        fHistoTrueClusElectronFromEtaPtvsM02[iCut]->Sumw2();
        fHistoTrueClusElectronFromGammaPtvsM02[iCut]->Sumw2();
        fHistoTrueClusBGPtvsSource[iCut]->Sumw2();
        fHistoTrueClusGammaPtvsSource[iCut]->Sumw2();
        fHistoTrueClusElectronPtvsSource[iCut]->Sumw2();
//        fHistoTrueClusElectronPtvsMotherID[iCut]->Sumw2();
//        fHistoTrueClusElectronPtvsTopMotherID[iCut]->Sumw2();
//        fHistoTrueClusElectronPtvsConvPhotonTopMotherID[iCut]->Sumw2();
        fHistoDoubleCountTruePi0PtvsM02[iCut]->Sumw2();
        fHistoDoubleCountTrueSecPi0Pt[iCut]->Sumw2();
        fHistoDoubleCountTrueMultiplePi0PtvsM02[iCut]->Sumw2();
        fHistoDoubleCountTrueMultipleSecPi0Pt[iCut]->Sumw2();
        fHistoDoubleCountTrueEtaPtvsM02[iCut]->Sumw2();

        if (GetSelectedMesonID() < 2){
          fHistoTrueClusPrimPi0PtvsM02[iCut]->Sumw2();
          fHistoTrueClusPrimPi0PtMCPt[iCut]->Sumw2();
          fHistoTruePureMergedClusPrimPi0PtvsM02[iCut]->Sumw2();
          fHistoTruePartConvClusPrimPi0PtvsM02[iCut]->Sumw2();
          if (fDoMesonQA > 1){
            fHistoTrueClusMergedPartConvPi0EVsM20[iCut]->Sumw2();
            fHistoTrueClusPrimPi0PureMergedPtMCPt[iCut]->Sumw2();
            fHistoTrueClusPrimPi0PartConvMergedPtMCPt[iCut]->Sumw2();
            fHistoTrueClusPrimPi01GammaMergedPtMCPt[iCut]->Sumw2();
            fHistoTrueClusPrimPi01ElectronMergedPtMCPt[iCut]->Sumw2();

            fHistoTrueClusNeutralContamination[iCut]->Sumw2();
            fHistoTrueClusPhotonContamination[iCut]->Sumw2();
            fHistoTrueClusElectronContamination[iCut]->Sumw2();
            fHistoTrueClusOtherContamination[iCut]->Sumw2();
          }
          fHistoTrueClusMultiplePrimPi0PtvsM02[iCut]->Sumw2();
          fHistoTrueClusSecPi0PtvsM02[iCut]->Sumw2();
          fHistoTrueClusSecPi0FromK0sPtvsM02[iCut]->Sumw2();
          fHistoTrueClusSecPi0FromK0lPtvsM02[iCut]->Sumw2();
          fHistoTrueClusSecPi0FromLambdaPtvsM02[iCut]->Sumw2();
          fHistoTrueSecPi0PtvsDiffReco[iCut]->Sumw2();
        }

        if (fDoMesonQA > 1){
          fHistoTrueClusMergedInvMassvsPt[iCut]->Sumw2();
          fHistoTrueClusPi0InvMassvsPt[iCut]->Sumw2();
          fHistoTrueClusEtaInvMassvsPt[iCut]->Sumw2();
          fHistoTrueClusBGInvMassvsPt[iCut]->Sumw2();
          fHistoTrueClusGammaInvMassvsPt[iCut]->Sumw2();
          fHistoTrueClusElectronInvMassvsPt[iCut]->Sumw2();
          if (GetSelectedMesonID() < 2){
            fHistoTrueClusPrimPi0InvMassvsPt[iCut]->Sumw2();
          }
          fHistoTrueClusGammaEvsM20[iCut]->Sumw2();
          fHistoTrueClusBGEvsM02[iCut]->Sumw2();
          fHistoTrueClusBGEvsM20[iCut]->Sumw2();
          fHistoTrueClusElectronEM20[iCut]->Sumw2();
          fHistoTrueClusMergedPi0EVsM20[iCut]->Sumw2();
          fHistoTrueClusMergedPi0EVsM02[iCut]->Sumw2();
          fHistoTrueClusMergedEtaEVsM20[iCut]->Sumw2();
          fHistoTrueClusMergedEtaEVsM02[iCut]->Sumw2();
        }

        if (fDoMesonQA > 0 ){
          if (GetSelectedMesonID() != 2){
            fHistoTruePi0PtY[iCut]->Sumw2();
            fHistoTruePi0PtAlpha[iCut]->Sumw2();
            fHistoTruePrimaryPi0MCPtResolPt[iCut]->Sumw2();
            fHistoTruePrimaryPi0PureMergedMCPtResolPt[iCut]->Sumw2();
            fHistoTruePrimaryPi0MergedPartConvMCPtResolPt[iCut]->Sumw2();
            fHistoTruePrimaryPi01GammaMCPtResolPt[iCut]->Sumw2();
            fHistoTruePrimaryPi01ElectronMCPtResolPt[iCut]->Sumw2();
            fHistoTrueSecondaryPi0MCPtResolPt[iCut]->Sumw2();
            fHistoTrueClusPi0EM02[iCut]->Sumw2();
            fHistoTrueClusEtaEM02[iCut]->Sumw2();
            fHistoTrueClusGammaEM02[iCut]->Sumw2();
            fHistoTrueClusElectronEM02[iCut]->Sumw2();
          }
          if (GetSelectedMesonID() != 1){
            fHistoTrueEtaPtY[iCut]->Sumw2();
            fHistoTrueEtaPtAlpha[iCut]->Sumw2();
            fHistoTruePrimaryEtaMCPtResolPt[iCut]->Sumw2();
          }
        }
        if(fDoJetAnalysis){
          if(!fDoLightOutput){
            fHistoTrueClusMergedPtvsRJet[iCut]->Sumw2();
            fHistoTrueClusPi0PtvsRJet[iCut]->Sumw2();
            fHistoTrueClusEtaPtvsRJet[iCut]->Sumw2();
          }
        }
      }
    }

  }

  if(fV0Reader)
    if((AliConvEventCuts*)fV0Reader->GetEventCuts())
      if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());

  for(Int_t iMatcherTask = 0; iMatcherTask < 5; iMatcherTask++){
    AliCaloTrackMatcher* temp = 0x0;
    if(!fCorrTaskSetting.CompareTo("")){
      temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcher_%i_%i",iMatcherTask,fTrackMatcherRunningMode)));
    } else {
      temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcher_%i_%i_%s",iMatcherTask,fTrackMatcherRunningMode,fCorrTaskSetting.Data())));
    }
    if(temp) fOutputContainer->Add(temp->GetCaloTrackMatcherHistograms());
  }


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
    if(fDoJetAnalysis){
      fCutFolder[iCut]->Add(fJetHistograms[iCut]);
    }
  }

  if (fIsMC > 0){
    tBrokenFiles = new TTree("BrokenFiles","BrokenFiles");
    tBrokenFiles->Branch("fileName",&fFileNameBroken);
    fOutputContainer->Add(tBrokenFiles);
  }


  PostData(1, fOutputContainer);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaCaloMerged::Notify()
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
    } else {
      printf(" Gamma Conversion Task %s :: Eta Shift Manually Set to %f \n\n",
          (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber()).Data(),((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift());
      fProfileEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
    }

    // Check whether non linearity correction has been applied already
    Bool_t doNonLinCorr   = kTRUE;
    if (((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetNonLinearity() == ((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))->GetNonLinearity()){
      doNonLinCorr   = kFALSE;
    } else if (((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetNonLinearity() == 0 && ((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))->GetNonLinearity() > 0){
      doNonLinCorr   = kTRUE;
    } else if (((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetNonLinearity() > 0 && ((AliCaloPhotonCuts*)fClusterMergedCutArray->At(iCut))->GetNonLinearity() == 0){
      doNonLinCorr   = kTRUE;
    } else {
      doNonLinCorr   = kFALSE;
      cout << "ERROR: something went wrong in the configuration, different non linearity corrections have been chosen in the standard and merged cluster, which are incompatible" << endl;
      cout << "INFO: switching off the non lin corr for merged cluster" << endl;
    }
    // check whether photon fullfill merged cluster cuts as well
    if (doNonLinCorr) ((AliCaloPhotonCuts*)fClusterMergedCutArray->At(fiCut))->SetUseNonLinearitySwitch(doNonLinCorr);


  }
  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::UserExec(Option_t *)
{
  //
  // Called for each event
  //
  fInputEvent = InputEvent();

  if(fIsMC> 0) fMCEvent = MCEvent();

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(fInputEvent->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event

  // Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete => abort processing of this file/event
  if(eventQuality == 2 || eventQuality == 3){
    // write out name of broken file for first event
    if (fIsMC > 0){
      if (fInputEvent->IsA()==AliESDEvent::Class()){
        if (((AliESDEvent*)fInputEvent)->GetEventNumberInFile() == 0){
          fFileNameBroken = new TObjString(Form("%s",((TString)fV0Reader->GetCurrentFileName()).Data()));
          if (tBrokenFiles) tBrokenFiles->Fill();
          delete fFileNameBroken;
        }
      }
    }

    for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      fHistoNEvents[iCut]->Fill(eventQuality);
      if (fIsMC==2)
        fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
    }
    return;
  }

  // ------------------- BeginEvent ----------------------------

  for(Int_t iCut = 0; iCut<fnCuts; iCut++){
    fiCut = iCut;

    Bool_t isRunningEMCALrelAna = kFALSE;
    if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;

    Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon, isRunningEMCALrelAna);

    fWeightJetJetMC = 1;
    //     cout << fMCEvent << endl;
    Float_t maxjetpt      = -1.;
    Float_t pthard = -1;
    if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetUseJetFinderForOutliers()) maxjetpt = fOutlierJetReader->GetMaxJetPt();
    Bool_t isMCJet        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC ,pthard, fInputEvent, maxjetpt);
    if (!isMCJet){
      fHistoNEvents[iCut]->Fill(10,fWeightJetJetMC);
      if (fIsMC==2) fHistoNEventsWOWeight[iCut]->Fill(10);
      continue;
    }

    if(fIsMC==2){
      Float_t xsection = -1.; Float_t ntrials = -1.;
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetXSectionAndNTrials(fMCEvent,xsection,ntrials,fInputEvent);
      if((xsection==-1.) || (ntrials==-1.)) AliFatal(Form("ERROR: GetXSectionAndNTrials returned invalid xsection/ntrials, periodName from V0Reader: '%s'",fV0Reader->GetPeriodName().Data()));
      fProfileJetJetXSection[iCut]->Fill(0.,xsection);
      fHistoJetJetNTrials[iCut]->Fill("#sum{NTrials}",ntrials);
    }

    Bool_t triggered = kTRUE;
    if(eventNotAccepted!= 0){
      fHistoNEvents[iCut]->Fill(eventNotAccepted, fWeightJetJetMC); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventNotAccepted);
      // cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
      if (eventNotAccepted==3 && fIsMC > 0){
        triggered = kFALSE;
      }else {
        continue;
      }
    }

    if(eventQuality != 0 && triggered== kTRUE){// Event Not Accepted
      //cout << "event rejected due to: " <<eventQuality << endl;
      fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC);
      if (fIsMC==2) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);

      continue;
    }

    // reject events if there are more overlaps of pi0 in R<0.05 than allowed
    if(fIsMC> 0 && NumberOfMCEventNeutralPionOverlapInEMCal(fMCEvent)){
      fHistoNEvents[iCut]->Fill(3, fWeightJetJetMC);
      if (fIsMC==2) fHistoNEventsWOWeight[iCut]->Fill(3);
      continue;
    }

    if (triggered == kTRUE) {
      fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC); // Should be 0 here
      if (fIsMC==2) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here

      fHistoNGoodESDTracks[iCut]->Fill(  fV0Reader->GetNumberOfPrimaryTracks(),
                        fWeightJetJetMC);
      fHistoVertexZ[iCut]->Fill(  fInputEvent->GetPrimaryVertex()->GetZ(),
                    fWeightJetJetMC);
      fHistoSPDClusterTrackletBackground[iCut]->Fill(  fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),
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

    if(fDoJetAnalysis)   ProcessJets();

    if(fIsMC>0){
      // ProcessNeutralOverlapsMC(fMCEvent);
      // create new vector that only contains those pi0s that would be lost due to overlap in the clusters
      if(fVectorLabelsLeadingPi0.size() > 0){
        for (std::vector<Int_t>::const_iterator i = fVectorLabelsMultiplePi0.begin(); i != fVectorLabelsMultiplePi0.end(); ++i){
          // check that the labels from fVectorLabelsMultiplePi0 are not contained in fVectorLabelsLeadingPi0 and fill fVectorLabelsMultiplePi0Reduced
          CheckVector1ForEntryAndFillVector2(fVectorLabelsLeadingPi0, fVectorLabelsMultiplePi0Reduced, (Int_t)*i);
        }
      }

      if(fInputEvent->IsA()==AliESDEvent::Class())
        ProcessMCParticles();
      if(fInputEvent->IsA()==AliAODEvent::Class())
        ProcessAODMCParticles();
    }

    if (triggered==kFALSE) continue;

    // it is in the loop to have the same conversion cut string (used also for MC stuff that should be same for V0 and Cluster)
    ProcessClusters();                      // process calo clusters


    fHistoNClusterCandidates[iCut]->Fill(  fNClusterCandidates,
                        fWeightJetJetMC);
    fHistoNClusterMergedCandidates[iCut]->Fill(  fNClusterMergedCandidates,
                        fWeightJetJetMC);

    fHistoNGoodESDTracksVsNClusterCandidates[iCut]->Fill(  fV0Reader->GetNumberOfPrimaryTracks(),
                                fNClusterCandidates,
                                fWeightJetJetMC);

    fVectorLabelsLeadingPi0.clear();
    fVectorLabelsMultiplePi0.clear();
    fVectorLabelsMultiplePi0Reduced.clear();
    fVectorDoubleCountTruePi0s.clear();
    fVectorDoubleCountTrueMultilePi0s.clear();
    fVectorDoubleCountTrueEtas.clear();

  }

  PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::ProcessClusters(){

  Int_t nclus = 0;
  Bool_t IsClusAcceptedByJet = kFALSE;
  if(!fCorrTaskSetting.CompareTo("")){
    nclus = fInputEvent->GetNumberOfCaloClusters();
  } else {
    if(!farrClustersProcess) farrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!farrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaCaloMerged! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclus = farrClustersProcess->GetEntries();
  }

  if(nclus == 0)  return;

  fNClusterCandidates       = 0;
  fNClusterMergedCandidates = 0;

  // plotting histograms on cell/tower level, only if extendedMatchAndQA > 1
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillHistogramsExtendedQA(fInputEvent,fIsMC);

  // match tracks to clusters
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchTracksToClusters(fInputEvent, fWeightJetJetMC, kTRUE, fMCEvent);
  // match tracks to clusters also for mergedCutArray
  ((AliCaloPhotonCuts*)fClusterMergedCutArray->At(fiCut))->MatchTracksToClusters(fInputEvent, fWeightJetJetMC, kTRUE, fMCEvent);

  // vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  // Loop over EMCal clusters
  for(Long_t i = 0; i < nclus; i++){

    Double_t tempClusterWeight        = fWeightJetJetMC;
    Double_t tempPhotonWeight         = fWeightJetJetMC;
    // select clusters with open cuts for normalizations histos
    AliVCluster* clus = NULL;
    if(fInputEvent->IsA()==AliESDEvent::Class()){
      if(farrClustersProcess)
        clus = new AliESDCaloCluster(*(AliESDCaloCluster*)farrClustersProcess->At(i));
      else
        clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
    } else if(fInputEvent->IsA()==AliAODEvent::Class()){
      if(farrClustersProcess)
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)farrClustersProcess->At(i));
      else
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
    }

    if (!clus) continue;



    // Set the jetjet weight to 1 in case the cluster orignated from the minimum bias header
    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
      if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(clus->GetLabelAt(0), fMCEvent, fInputEvent) == 2)
        tempClusterWeight = 1;
    }

    // if open cluster cuts are not fullfilled I can abort
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC,tempClusterWeight, i)){
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
      continue;
    }

    // Flag Photon as CaloPhoton
    PhotonCandidate->SetIsCaloPhoton(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType());
    PhotonCandidate->SetCaloClusterRef(i);

    // get MC label
    if(fIsMC> 0){
      Int_t* mclabelsCluster  = clus->GetLabels();
      // cout << "cluster: " << i << endl;
      Int_t nValidClusters    = 0;

      if (clus->GetNLabels()>0){
        for (Int_t k =0; k< (Int_t)clus->GetNLabels(); k++){
          // TParticle *dummy    = NULL;
          if (mclabelsCluster[k]>0){
            // dummy             = fMCEvent->Particle(mclabelsCluster[k]);
            // if (dummy->R() < 407.0){
              PhotonCandidate->SetCaloPhotonMCLabel(nValidClusters,mclabelsCluster[k]);
              nValidClusters++;
            // }
          }
        }
      }
      PhotonCandidate->SetNCaloPhotonMCLabels(nValidClusters);

    }

    fIsFromMBHeader = kTRUE;
    fIsOverlappingWithOtherHeader = kFALSE;

    // test whether largest contribution to cluster orginates in added signals
    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() > 0){
      // Set the jetjet weight to 1 in case the photon candidate orignated from the minimum bias header
      if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempPhotonWeight = 1;
      if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 0) fIsFromMBHeader = kFALSE;
      if (clus->GetNLabels()>1){
        Int_t* mclabelsCluster = clus->GetLabels();
        for (Int_t l = 1; l < (Int_t)clus->GetNLabels(); l++ ){
          if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(mclabelsCluster[l], fMCEvent, fInputEvent) == 0) fIsOverlappingWithOtherHeader = kTRUE;
        }
      }
    }

    // test whether neutral pion overlaps with other neutral pions
    if (fIsMC>0 && fMesonOLFromCluster){
      AliAODMCParticle* dummyPart        = (AliAODMCParticle*) fAODMCTrackArray->At(clus->GetLabelAt(0));
      AliAODMCParticle* dummyMother        = (AliAODMCParticle*) fAODMCTrackArray->At(dummyPart->GetMother());
      if( ( (fMinAllowedPi0OverlapsMC>-1) ? (fMapNeutralPionOverlap[fiCut][dummyMother->GetLabel()] < fMinAllowedPi0OverlapsMC ) : kFALSE ) || ( (fMaxAllowedPi0OverlapsMC>-1) ? (fMapNeutralPionOverlap[fiCut][dummyMother->GetLabel()] > fMaxAllowedPi0OverlapsMC ) : kFALSE ) ){
        delete clus;
        delete tmpvec;
        delete PhotonCandidate;
        continue;
      }
    }

    // check whether photon is from correct header
    if (fIsFromMBHeader && fIsOverlappingWithOtherHeader){
      fHistoClusOverlapHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), tempPhotonWeight);
    }
    // only cluster candidates from correct header will be processed fully
    if (fIsFromMBHeader && !fIsOverlappingWithOtherHeader){
      fHistoClusGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), tempPhotonWeight);
      fHistoClusGammaE[fiCut]->Fill(PhotonCandidate->E(), tempPhotonWeight);
      if (fDoClusterQA > 0){
        fHistoClusNCellsPt[fiCut]->Fill(clus->GetNCells(), PhotonCandidate->Pt(), tempPhotonWeight);
      }
      fNClusterCandidates++;
    } else {
      delete clus;
      delete tmpvec;
      delete PhotonCandidate;
      continue;
    }

    Int_t matchedJet = 0;
    if(fDoJetAnalysis){
      if(fDoOutOfJet == 1) IsClusAcceptedByJet = kTRUE;
      else IsClusAcceptedByJet = kFALSE;
      Float_t clusPos[3]={0,0,0};
      clus->GetPosition(clusPos);
      TVector3 clusterVectorJets(clusPos[0],clusPos[1],clusPos[2]);
      Double_t etaCluster = clusterVectorJets.Eta();
      Double_t phiCluster = (clusterVectorJets.Phi() > 0) ? clusterVectorJets.Phi() : clusterVectorJets.Phi() + 2*TMath::Pi();
      if(fConvJetReader->GetNJets()>0){
        fVectorJetEta = fConvJetReader->GetVectorJetEta();
        fVectorJetPhi = fConvJetReader->GetVectorJetPhi();
        for(Int_t ij=0; ij<fConvJetReader->GetNJets(); ij++){
          Double_t DeltaEta = fVectorJetEta.at(ij)-etaCluster;
          Double_t DeltaPhi = abs(fVectorJetPhi.at(ij)-phiCluster);
          if(DeltaPhi > TMath::Pi()) {
            DeltaPhi = 2*TMath::Pi() - DeltaPhi;
          }
          if(fDoOutOfJet == 2){ // check if on opposite side of jet (DeltaEta/Phi = 0 if directly opposite)
            DeltaEta = fVectorJetEta.at(ij) + etaCluster;
            DeltaPhi = abs(TMath::Pi() - DeltaPhi);
          }
          Double_t RJetPi0Cand = TMath::Sqrt(DeltaEta*DeltaEta+DeltaPhi*DeltaPhi);
          if(fConvJetReader->Get_Jet_Radius() > 0 ){
            if(fDoOutOfJet == 0){ // in jet
              if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                IsClusAcceptedByJet = kTRUE;
                matchedJet = ij;
                break;
              }
            } else if(fDoOutOfJet == 1){ // out of jet
              if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                IsClusAcceptedByJet = kFALSE;
                matchedJet = ij;
                break;
              }
            } else if(fDoOutOfJet == 2){ // out of jet on away side
              if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                IsClusAcceptedByJet = kTRUE;
                matchedJet = ij;
                break;
              }
            } else if(fDoOutOfJet == 3){ // out of jet in interval [R, R+0.2]
              if((RJetPi0Cand > fConvJetReader->Get_Jet_Radius()) && (RJetPi0Cand < fConvJetReader->Get_Jet_Radius() + 0.2)){
                IsClusAcceptedByJet = kTRUE;
                matchedJet = ij;
                break;
              }
            }
          }
        }
        fVectorJetEta.clear();
        fVectorJetPhi.clear();
      }
      if(!IsClusAcceptedByJet){
        delete clus;
        delete tmpvec;
        delete PhotonCandidate;
        continue;
      }
      if(!fDoLightOutput)fClusterEtaPhiJets[fiCut]->Fill(etaCluster, phiCluster);
    }

    if(!((AliCaloPhotonCuts*)fClusterMergedCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC,tempPhotonWeight,i)){
      delete clus;
      delete tmpvec;
      delete PhotonCandidate;
      continue;
    }else {
      fNClusterMergedCandidates++;
    }

    AliAODCaloCluster* clusSub1 = new AliAODCaloCluster();
    AliAODCaloCluster* clusSub2 = new AliAODCaloCluster();

    const Int_t   nc = clus->GetNCells();
    Int_t     absCellIdList[nc];
    Float_t   cellEnergyList[nc];
    Int_t nLM = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetNumberOfLocalMaxima(clus, fInputEvent, absCellIdList, cellEnergyList);

    // split clusters according to their shares in the cluster (NLM == 1) needs to be treated differently
    if ( (((AliCaloPhotonCuts*)fClusterMergedCutArray->At(fiCut))->GetMinNLMCut() == 1 &&
      ((AliCaloPhotonCuts*)fClusterMergedCutArray->At(fiCut))->GetMaxNLMCut() == 1 ) ||
      nLM == 1){
      Int_t absCellIdFirst    = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FindLargestCellInCluster(clus, fInputEvent);
      Int_t absCellIdSecond   = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FindSecondLargestCellInCluster(clus, fInputEvent);

      ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->SplitEnergy(absCellIdFirst, absCellIdSecond,
                                      clus, fInputEvent, fIsMC, clusSub1, clusSub2);


    } else {
      ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->SplitEnergy(absCellIdList[0], absCellIdList[1],
                                                                     clus, fInputEvent, fIsMC, clusSub1, clusSub2);
    }
    fHistoClusMergedPtvsM02[fiCut]->Fill( PhotonCandidate->Pt(),
                          clus->GetM02(),
                          tempPhotonWeight);
    if(fDoMesonQA > 1){
      fHistoClusMergedEvsM02[fiCut]->Fill( PhotonCandidate->E(),
                            clus->GetM02(),
                            tempPhotonWeight);
      fHistoClusMergedEvsM20[fiCut]->Fill( PhotonCandidate->E(),
                            clus->GetM20(),
                            tempPhotonWeight);
    }

    // TLorentzvector with sub cluster 1
    TLorentzVector clusterVector1;
    clusSub1->GetMomentum(clusterVector1,vertex);
    TLorentzVector* tmpvec1 = new TLorentzVector();
    tmpvec1->SetPxPyPzE(clusterVector1.Px(),clusterVector1.Py(),clusterVector1.Pz(),clusterVector1.E());

    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate1=new AliAODConversionPhoton(tmpvec1);
    if(!PhotonCandidate1){
      delete clusSub1;
      delete tmpvec1;
      continue;
    }
    // Flag Photon as CaloPhoton
    PhotonCandidate1->SetIsCaloPhoton(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType());

    // TLorentzvector with sub cluster 2
    TLorentzVector clusterVector2;
    clusSub2->GetMomentum(clusterVector2,vertex);
    TLorentzVector* tmpvec2 = new TLorentzVector();
    tmpvec2->SetPxPyPzE(clusterVector2.Px(),clusterVector2.Py(),clusterVector2.Pz(),clusterVector2.E());
    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate2=new AliAODConversionPhoton(tmpvec2);
    if(!PhotonCandidate2){
      delete clusSub2;
      delete tmpvec2;
      continue;
    }
    // Flag Photon as CaloPhoton
    PhotonCandidate1->SetIsCaloPhoton(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType());
    // create pi0 candidate
    AliAODConversionMother *pi0cand = new AliAODConversionMother(PhotonCandidate1,PhotonCandidate2);

    if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()))){
      fHistoMotherInvMassPt[fiCut]->Fill(pi0cand->M(),PhotonCandidate->Pt(), tempPhotonWeight);
      fHistoClusMergedPtvsM02Accepted[fiCut]->Fill( PhotonCandidate->Pt(), clus->GetM02(), tempPhotonWeight);
      fHistoClusMergedEvsM02Accepted[fiCut]->Fill( PhotonCandidate->E(), clus->GetM02(), tempPhotonWeight);

      if(fDoJetAnalysis && !fDoLightOutput && (fDoOutOfJet != 1 )){
        fVectorJetPt  = fConvJetReader->GetVectorJetPt();
        fVectorJetPx  = fConvJetReader->GetVectorJetPx();
        fVectorJetPy  = fConvJetReader->GetVectorJetPy();
        fVectorJetPz  = fConvJetReader->GetVectorJetPz();
        fVectorJetEta = fConvJetReader->GetVectorJetEta();
        fVectorJetPhi = fConvJetReader->GetVectorJetPhi();
        Double_t DeltaEta = fVectorJetEta.at(matchedJet)-PhotonCandidate->Eta();
        Double_t DeltaPhi = abs(fVectorJetPhi.at(matchedJet)-PhotonCandidate->Phi());
        if(DeltaPhi > M_PI) {
            DeltaPhi = 2*M_PI - DeltaPhi;
        }
        Double_t RJetPi0Cand = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
        if(fDoOutOfJet == 2){
          RJetPi0Cand = abs(TMath::Pi() - RJetPi0Cand);
        }
        if(fConvJetReader->Get_Jet_Radius() > 0){
          if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius() || (fDoOutOfJet == 3 && RJetPi0Cand > fConvJetReader->Get_Jet_Radius()) ){
            Double_t dotproduct = fVectorJetPx.at(matchedJet)*pi0cand->Px() + fVectorJetPy.at(matchedJet)*pi0cand->Py() + fVectorJetPz.at(matchedJet)*pi0cand->Pz();
            Double_t magn = pow(fVectorJetPx.at(matchedJet),2) + pow(fVectorJetPy.at(matchedJet),2) + pow(fVectorJetPz.at(matchedJet),2);
            Double_t z = dotproduct/magn;
            fHistoJetFragmFunc[fiCut]->Fill(PhotonCandidate->Pt(), z, tempPhotonWeight);
            fHistoClusMergedPtvsRJetAccepted[fiCut]->Fill(PhotonCandidate->Pt(), RJetPi0Cand, tempPhotonWeight);
          }
        }
      }
      Int_t nlm = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetNumberOfLocalMaxima(clus, fInputEvent);
      if ( nlm < 11 ){
        fHistoClusNLMPt[fiCut]->Fill( nlm, PhotonCandidate->Pt());
      } else {
        fHistoClusNLMPt[fiCut]->Fill( 11., PhotonCandidate->Pt());
      }

      if (fDoClusterQA > 0){
        fHistoClusMergedNCellsPt[fiCut]->Fill(clus->GetNCells(), PhotonCandidate->Pt(), tempPhotonWeight);
        Double_t energyAround   = 0;
        Double_t nCellsAround   = 0;
        for (Int_t j = 0; j < nclus; j++){
          if (j == i) continue;
          AliVCluster* clusTemp   = NULL;
          if(fInputEvent->IsA()==AliESDEvent::Class()){
            if(farrClustersProcess)
              clusTemp = new AliESDCaloCluster(*(AliESDCaloCluster*)farrClustersProcess->At(j));
            else
              clusTemp = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(j));
          } else if(fInputEvent->IsA()==AliAODEvent::Class()){
            if(farrClustersProcess)
              clusTemp = new AliAODCaloCluster(*(AliAODCaloCluster*)farrClustersProcess->At(j));
            else
              clusTemp = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(j));
          }

          if (!clusTemp) continue;

          Double_t R = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDistanceBetweenClusters(clus, clusTemp);

          if (R < 0.15){
            nCellsAround        = nCellsAround+clusTemp->GetNCells();
            energyAround        = energyAround+clusTemp->E();
          }
          delete clusTemp;
        }
        fHistoClusMergedNCellsAroundPt[fiCut]->Fill(nCellsAround, PhotonCandidate->Pt(), tempPhotonWeight);
        fHistoClusMergedNCellsAroundAndInPt[fiCut]->Fill(nCellsAround+clus->GetNCells(), PhotonCandidate->Pt(), tempPhotonWeight);
        fHistoClusMergedEAroundE[fiCut]->Fill(energyAround, clus->E(), tempPhotonWeight);
        if (fIsMC > 0){
          fHistoClusMergedNParticlePt[fiCut]->Fill(clus->GetNLabels(), PhotonCandidate->Pt(), tempPhotonWeight);
        }
      }
      if (fDoMesonQA > 0 ){
        fHistoMotherPtY[fiCut]->Fill(PhotonCandidate->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), tempPhotonWeight);
        fHistoMotherPtAlpha[fiCut]->Fill(PhotonCandidate->Pt(),TMath::Abs(pi0cand->GetAlpha()), tempPhotonWeight);
      }
      if(fIsMC> 0 && PhotonCandidate && PhotonCandidate1 && PhotonCandidate2){

        // PrintFractionsOfMCLabels(clus,fInputEvent,fMCEvent);

        if(fInputEvent->IsA()==AliESDEvent::Class())
          ProcessTrueClusterCandidates(PhotonCandidate, clus, PhotonCandidate1, PhotonCandidate2);
        if(fInputEvent->IsA()==AliAODEvent::Class())
          ProcessTrueClusterCandidatesAOD(PhotonCandidate, clus, PhotonCandidate1, PhotonCandidate2);
      }
    } // meson is selected



    delete pi0cand;
    delete clusSub1;
    delete tmpvec1;
    delete PhotonCandidate1;
    delete clusSub2;
    delete tmpvec2;
    delete PhotonCandidate2;

    delete clus;
    delete tmpvec;
    delete PhotonCandidate;
  } // cluster loop
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::ProcessJets()
{
  fHistoNJets[fiCut]->Fill(fConvJetReader->GetNJets());
  if(fConvJetReader->GetNJets()>0){
    fVectorJetPt  = fConvJetReader->GetVectorJetPt();
    fVectorJetEta = fConvJetReader->GetVectorJetEta();
    fVectorJetPhi = fConvJetReader->GetVectorJetPhi();
    fVectorJetArea = fConvJetReader->GetVectorJetArea();
    if(fIsMC > 0 && fConvJetReader->GetTrueNJets()>0){
      fTrueVectorJetPt = fConvJetReader->GetTrueVectorJetPt();
      fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
      fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
    }
    if(fVectorJetPt.size() == fConvJetReader->GetNJets() && fVectorJetEta.size() == fConvJetReader->GetNJets() && fVectorJetPhi.size() == fConvJetReader->GetNJets() && fVectorJetArea.size() == fConvJetReader->GetNJets()){

      for(Int_t i=0; i<fConvJetReader->GetNJets(); i++){
        fHistoPtJet[fiCut]->Fill(fVectorJetPt.at(i));
        fHistoJetEta[fiCut]->Fill(fVectorJetEta.at(i));
        fHistoJetPhi[fiCut]->Fill(fVectorJetPhi.at(i));
        fHistoJetArea[fiCut]->Fill(fVectorJetArea.at(i));
        if(fIsMC > 0 && fConvJetReader->GetNJets()>0 && fConvJetReader->GetTrueNJets()>0){
          Double_t min = 100;
          Int_t match = 0;
          for(Int_t j = 0; j<fConvJetReader->GetTrueNJets(); j++){
            Double_t R_jetjet;
            Double_t DeltaEta = fVectorJetEta.at(i)-fTrueVectorJetEta.at(j);
            Double_t DeltaPhi = abs(fVectorJetPhi.at(i)-fTrueVectorJetPhi.at(j));
            if(DeltaPhi > TMath::Pi()) {
              DeltaPhi = 2*TMath::Pi() - DeltaPhi;
            }
            R_jetjet = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
            if(R_jetjet < min){
              min = R_jetjet;
              match = j;
            }
          }
          if(!fDoLightOutput) fHistoTruevsRecJetPt[fiCut]->Fill(fVectorJetPt.at(i), fTrueVectorJetPt.at(match));
          if(fDoJetQA){
            if(fVectorJetPt.at(i) >= 10) fHistoEventwJets[fiCut]->Fill(0.5);
            if(fVectorJetPt.at(i) < 10 && fTrueVectorJetPt.at(match) >= 10) fHistoEventwJets[fiCut]->Fill(1.5);
            if(fVectorJetPt.at(i) >= 10 && fTrueVectorJetPt.at(match) < 10) fHistoEventwJets[fiCut]->Fill(2.5);
          }
        }
      }
    }

    fVectorJetPt.clear();
    fVectorJetEta.clear();
    fVectorJetPhi.clear();
    fVectorJetArea.clear();
    if(fIsMC > 0 && fConvJetReader->GetTrueNJets()>0){
      fTrueVectorJetPt.clear();
      fTrueVectorJetEta.clear();
      fTrueVectorJetPhi.clear();
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::ProcessTrueClusterCandidates(AliAODConversionPhoton *TrueClusterCandidate, AliVCluster* cluster,
                                    AliAODConversionPhoton *TrueSubClusterCandidate1,
                                    AliAODConversionPhoton *TrueSubClusterCandidate2)
{
  Float_t m02 = cluster->GetM02();
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  Double_t tempClusterWeight       = fWeightJetJetMC;
  TParticle *Photon = NULL;
  if (TrueClusterCandidate->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
  if (TrueClusterCandidate->GetCaloPhotonMCLabel(0) < 0) return;

  AliAODConversionMother *mesoncand = new AliAODConversionMother(TrueSubClusterCandidate1,TrueSubClusterCandidate2);
  Bool_t mesonIsSelected            = (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(mesoncand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()));
  if (!mesonIsSelected){
    delete mesoncand;
    return;
  }

  // Setting all MC Flags (IsMerged, etc)
  TrueClusterCandidate->SetCaloPhotonMCFlags(fMCEvent, fEnableSortForClusMC, kTRUE,cluster);

  if (TrueClusterCandidate->GetNCaloPhotonMCLabels()>0 && fEnableDetailedPrintOut){
    cout << endl << endl << "Cluster energy: " << TrueClusterCandidate->E() << endl;;
    TrueClusterCandidate->PrintCaloMCLabelsAndInfo(fMCEvent);
    TrueClusterCandidate->PrintCaloMCFlags();
  }


  if (TrueClusterCandidate->GetNCaloPhotonMCLabels()>0){
    fHistoTrueClusEFracFirstLabel[fiCut]->Fill(cluster->E(), cluster->GetClusterMCEdepFraction(0), tempClusterWeight);
    if(TrueClusterCandidate->GetNNeutralPionMCLabels()>0)
      fHistoTrueClusEFracLeadingPi0[fiCut]->Fill(cluster->E(), TrueClusterCandidate->GetNeutralPionEnergyFraction(TrueClusterCandidate->GetLeadingNeutralPionIndex()), tempClusterWeight);
    // check if leading pi0 comes not from label 0 in cluster
    // for this do:
    // -> check if neutral pions were found in cluster
    // -> if the leading daughter index is not 0
    // -> the leading neutral pion has a larger cluster energy fraction than the cluster label 0
    if( TrueClusterCandidate->GetNNeutralPionMCLabels()>0 && TrueClusterCandidate->GetLeadingNeutralPionDaughterIndex()!=0 && TrueClusterCandidate->GetNeutralPionEnergyFraction(TrueClusterCandidate->GetLeadingNeutralPionIndex())>cluster->GetClusterMCEdepFraction(0)) {
      // load particle corresponding to largest daughter of leading pi0
      Photon         = fMCEvent->Particle(TrueClusterCandidate->GetCaloPhotonMCLabel(TrueClusterCandidate->GetLeadingNeutralPionDaughterIndex()));
      // fill vector with leading pi0 MC label
      CheckVectorForDoubleCount(fVectorLabelsLeadingPi0,TrueClusterCandidate->GetNeutralPionMCLabel(TrueClusterCandidate->GetLeadingNeutralPionIndex()));
      // fill vector will ALL pi0 MC labels that are found in the cluster
      for(Int_t npion = 0; npion < TrueClusterCandidate->GetNNeutralPionMCLabels(); npion++){
        CheckVectorForDoubleCount(fVectorLabelsMultiplePi0,TrueClusterCandidate->GetNeutralPionMCLabel(npion));
      }
    } else {
      // load particle corresponding to MC label 0 in cluster
      Photon         = fMCEvent->Particle(TrueClusterCandidate->GetCaloPhotonMCLabel(0));
      // important check that leading daughter corresponds to the cluster MC label 0
      if(TrueClusterCandidate->GetLeadingNeutralPionDaughterIndex()==0){
        // fill vector with leading pi0 MC label
        CheckVectorForDoubleCount(fVectorLabelsLeadingPi0,TrueClusterCandidate->GetNeutralPionMCLabel(TrueClusterCandidate->GetLeadingNeutralPionIndex()));
        // fill vector will ALL pi0 MC labels that are found in the cluster
        for(Int_t npion = 0; npion < TrueClusterCandidate->GetNNeutralPionMCLabels(); npion++){
          CheckVectorForDoubleCount(fVectorLabelsMultiplePi0,TrueClusterCandidate->GetNeutralPionMCLabel(npion));
        }
      }
    }
  } else {
    // return if there are no MC labels in the cluster
    return;
  }

  if(Photon == NULL){
    return;
  }

  // Get Distance to Jet
  Double_t RJetPi0Cand = 0.;
  if(fDoJetAnalysis && !fDoLightOutput){
    for(Int_t j=0; j<fConvJetReader->GetNJets(); j++){
      Double_t DeltaEta = fVectorJetEta.at(j)-Photon->Eta();
      Double_t DeltaPhi = abs(fVectorJetPhi.at(j)-Photon->Phi());
      if(DeltaPhi > TMath::Pi()) {
        DeltaPhi = 2*TMath::Pi() - DeltaPhi;
      }
      RJetPi0Cand = TMath::Sqrt(DeltaEta*DeltaEta+DeltaPhi*DeltaPhi);
      if(fConvJetReader->Get_Jet_Radius() > 0 ){
        if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
          break;
        }
      }
    }
  }


  Int_t pdgCodeParticle             = Photon->GetPdgCode();

  if(fIsFromMBHeader && !fIsOverlappingWithOtherHeader){
    Int_t clusterClass    = 0;
    Bool_t isPrimary      = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, TrueClusterCandidate->GetCaloPhotonMCLabel(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);

    // cluster classification:
    // 1    - nice merged cluster (2 gamma | contributions from 2 gamma) from pi0/eta
    // 2    - contribution from only 1 partner (1 gamma, 1 fully coverted gamma) from pi0/eta
    // 3    - contribution from part of 1 partner (1 electron) from pi0/eta
    Long_t motherLab = -1;
    if (TrueClusterCandidate->IsMerged() || TrueClusterCandidate->IsMergedPartConv()){
        clusterClass    = 1;
        motherLab       = TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0);
    } else if (TrueClusterCandidate->GetNCaloPhotonMotherMCLabels()> 0){
      // cout << TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0) << endl;
      if (TrueClusterCandidate->IsLargestComponentElectron() || TrueClusterCandidate->IsLargestComponentPhoton()){
        if (TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0) > -1 && (fMCEvent->Particle(TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0))->GetPdgCode() == 111 || fMCEvent->Particle(TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0))->GetPdgCode() == 221) ){
          if ( TrueClusterCandidate->IsConversion() && !TrueClusterCandidate->IsConversionFullyContained() ){
            clusterClass  = 3;
            motherLab       = TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0);
          } else {
            clusterClass  = 2;
            motherLab       = TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0);
          }
        }
      } else if (TrueClusterCandidate->IsSubLeadingEM()){
        if (TrueClusterCandidate->GetNCaloPhotonMotherMCLabels()> 1){
          if (fEnableDetailedPrintOut) cout << "Is Subleading EM: "<<  TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1) << endl;
          if ( TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1) > -1){
            if (fMCEvent->Particle(TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1))->GetPdgCode() == 111 || fMCEvent->Particle(TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1))->GetPdgCode() == 221 ){
              clusterClass  = 2;
              motherLab       = TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1);
            }
          }
        }
      } else {
        motherLab       = TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0);
      }
    }

    // Get Mother particle
    TParticle *mother = NULL;
    Int_t motherPDG   = -1;
    if (motherLab > -1){
       mother           = fMCEvent->Particle(motherLab);
    }
    if (fEnableDetailedPrintOut) cout << "cluster class: " << clusterClass << "\t mother lab: "<< motherLab ;
    if (mother){
        motherPDG = TMath::Abs(mother->GetPdgCode());
        if (fEnableDetailedPrintOut) cout  << "\t mother pdg: " << motherPDG << endl;
    } else {
      if (fEnableDetailedPrintOut) cout << endl;
    }
    // Set the jetjet weight to 1 in case the photon orignated from the minimum bias header
    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
      if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(motherLab, fMCEvent, fInputEvent) == 2) tempClusterWeight = 1;
    }

    //
    if (clusterClass == 1 || clusterClass == 2 || clusterClass == 3 ){
      fHistoTrueClusMergedPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
      if (fDoMesonQA > 1)fHistoTrueClusMergedInvMassvsPt[fiCut]->Fill(mesoncand->M(),TrueClusterCandidate->Pt(), tempClusterWeight);
      if(fDoJetAnalysis && !fDoLightOutput) fHistoTrueClusMergedPtvsRJet[fiCut]->Fill(TrueClusterCandidate->Pt(), RJetPi0Cand, tempClusterWeight);

      // separate different components
      if (clusterClass == 1 && TrueClusterCandidate->IsMerged()){
        if (motherPDG == 111){
          fHistoTrueClusMergedPureFromPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
          if (GetSelectedMesonID() < 2 && !isPrimary && m02 >= 0 && m02 <= 4.8 ){
            fHistoTrueSecPi0PtvsDiffReco[fiCut]->Fill(TrueClusterCandidate->Pt(), 0., tempClusterWeight);
          }
        }
        if (motherPDG == 221)
          fHistoTrueClusMergedPureFromEtaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
      } else if (clusterClass == 1 && TrueClusterCandidate->IsMergedPartConv()){
        if (motherPDG == 111){
          fHistoTrueClusMergedPartConvFromPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
          if (GetSelectedMesonID() < 2 && !isPrimary && m02 >= 0 && m02 <= 4.8 ){
            fHistoTrueSecPi0PtvsDiffReco[fiCut]->Fill(TrueClusterCandidate->Pt(), 1., tempClusterWeight);
          }
        }
        if (motherPDG == 221)
          fHistoTrueClusMergedPartConvFromEtaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
      } else if (clusterClass == 2){
        if (motherPDG == 111){
          fHistoTrueClusGammaFromPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
          if (GetSelectedMesonID() < 2 && !isPrimary && m02 >= 0 && m02 <= 4.8 ){
            fHistoTrueSecPi0PtvsDiffReco[fiCut]->Fill(TrueClusterCandidate->Pt(), 2., tempClusterWeight);
          }
        }if (motherPDG == 221)
          fHistoTrueClusGammaFromEtaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
      } else if (clusterClass == 3){
        if (motherPDG == 111) {
          fHistoTrueClusElectronFromPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
          if (GetSelectedMesonID() < 2 && !isPrimary && m02 >= 0 && m02 <= 4.8 ){
            fHistoTrueSecPi0PtvsDiffReco[fiCut]->Fill(TrueClusterCandidate->Pt(), 3., tempClusterWeight);
          }
        }
        if (motherPDG == 221)
          fHistoTrueClusElectronFromEtaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
      }

      // deal with pi0 only
      if (motherPDG == 111){
        fHistoTrueClusPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
        if(fDoJetAnalysis && !fDoLightOutput) fHistoTrueClusPi0PtvsRJet[fiCut]->Fill(TrueClusterCandidate->Pt(), RJetPi0Cand, tempClusterWeight);
        if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,motherLab)){
          fHistoDoubleCountTruePi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
          if (!isPrimary)
            fHistoDoubleCountTrueSecPi0Pt[fiCut]->Fill(TrueClusterCandidate->Pt(), tempClusterWeight);
        }

        // Treatment of multiple pi0 in cluster (filling of true candidates, amount of pions and double counting)
        if(TrueClusterCandidate->GetNNeutralPionMCLabels()>1){
          for(Int_t npion = 0; npion < TrueClusterCandidate->GetNNeutralPionMCLabels(); npion++){
            // fill multiple true pions in same cluster
            if(npion>0) fHistoTrueClusMultiplePi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);

            // check MC labels for double counting (each particle should only be reconstructed once)
            if (CheckVectorForDoubleCount(fVectorDoubleCountTrueMultilePi0s,TrueClusterCandidate->GetNeutralPionMCLabel(npion))){
              fHistoDoubleCountTrueMultiplePi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);

              // check if particle is not a primary -> secondary particle
              if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, TrueClusterCandidate->GetNeutralPionMCLabel(npion), mcProdVtxX, mcProdVtxY, mcProdVtxZ)) // check if pion is primary
                fHistoDoubleCountTrueMultipleSecPi0Pt[fiCut]->Fill(TrueClusterCandidate->Pt(), tempClusterWeight);
            }
          }
        }
        fHistoNTrueMultiplePi0vsPt[fiCut]->Fill(TrueClusterCandidate->Pt(),TrueClusterCandidate->GetNNeutralPionMCLabels(), tempClusterWeight);


        if (TrueClusterCandidate->IsDalitz()){
          fHistoTrueClusPi0DalitzPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
        }
        if (fDoMesonQA > 0){
          fHistoTrueClusPi0EM02[fiCut]->Fill(TrueClusterCandidate->E(), m02, tempClusterWeight);
          if (fDoMesonQA > 1){
            fHistoTrueClusPi0InvMassvsPt[fiCut]->Fill(mesoncand->M(),TrueClusterCandidate->Pt(), tempClusterWeight);
          }
        }
        if (fDoMesonQA > 0 && GetSelectedMesonID() != 2){
          fHistoTruePi0PtY[fiCut]->Fill(TrueClusterCandidate->Pt(),mesoncand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), tempClusterWeight);
          fHistoTruePi0PtAlpha[fiCut]->Fill(TrueClusterCandidate->Pt(),TMath::Abs(mesoncand->GetAlpha()), tempClusterWeight);
        }

        if (GetSelectedMesonID() < 2) {
          Float_t weighted= 1;
          if (mother->Pt()>0.005){
            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(motherLab, fMCEvent, fInputEvent);
          }
          if (isPrimary) {
            fHistoTrueClusPrimPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight*weighted);
            if (clusterClass == 1 && TrueClusterCandidate->IsMerged()){
              fHistoTruePureMergedClusPrimPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight*weighted);
              if (fDoMesonQA > 1){
                fHistoTrueClusMergedPi0EVsM20[fiCut]->Fill(TrueClusterCandidate->E(), cluster->GetM20(), tempClusterWeight*weighted);
                fHistoTrueClusMergedPi0EVsM02[fiCut]->Fill(TrueClusterCandidate->E(), m02, tempClusterWeight*weighted);
              }
            }else if (clusterClass == 1 && TrueClusterCandidate->IsMergedPartConv()){
              fHistoTruePartConvClusPrimPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight*weighted);
              if (fDoMesonQA > 1){
                fHistoTrueClusMergedPartConvPi0EVsM20[fiCut]->Fill(TrueClusterCandidate->E(), cluster->GetM20(), tempClusterWeight);
              }
            }
            if(TrueClusterCandidate->GetNNeutralPionMCLabels()>1){;
              for(Int_t npion = 1; npion < TrueClusterCandidate->GetNNeutralPionMCLabels(); npion++){
                if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, TrueClusterCandidate->GetNeutralPionMCLabel(npion), mcProdVtxX, mcProdVtxY, mcProdVtxZ)) // check if pion is primary
                  fHistoTrueClusMultiplePrimPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
              }
            }
            if (fDoMesonQA > 1) fHistoTrueClusPrimPi0InvMassvsPt[fiCut]->Fill(mesoncand->M(),TrueClusterCandidate->Pt(), tempClusterWeight);
            if (fDoMesonQA > 0 && mother->Pt()>0){
              fHistoTruePrimaryPi0MCPtResolPt[fiCut]->Fill(mother->Pt(),(TrueClusterCandidate->Pt()-mother->Pt())/mother->Pt(),tempClusterWeight);
                if (fDoMesonQA > 1) fHistoTrueClusPrimPi0PtMCPt[fiCut]->Fill(TrueClusterCandidate->Pt(), mother->Pt(),tempClusterWeight);
              if (clusterClass == 1 && TrueClusterCandidate->IsMerged()){
                fHistoTruePrimaryPi0PureMergedMCPtResolPt[fiCut]->Fill(mother->Pt(),(TrueClusterCandidate->Pt()-mother->Pt())/mother->Pt(),tempClusterWeight);
                if (fDoMesonQA > 1) fHistoTrueClusPrimPi0PureMergedPtMCPt[fiCut]->Fill(TrueClusterCandidate->Pt(), mother->Pt(),tempClusterWeight);
              }else if (clusterClass == 1 && TrueClusterCandidate->IsMergedPartConv()){
                fHistoTruePrimaryPi0MergedPartConvMCPtResolPt[fiCut]->Fill(mother->Pt(),(TrueClusterCandidate->Pt()-mother->Pt())/mother->Pt(),tempClusterWeight);
                if (fDoMesonQA > 1) fHistoTrueClusPrimPi0PartConvMergedPtMCPt[fiCut]->Fill(TrueClusterCandidate->Pt(), mother->Pt(),tempClusterWeight);
              }else if (clusterClass == 2){
                fHistoTruePrimaryPi01GammaMCPtResolPt[fiCut]->Fill(mother->Pt(),(TrueClusterCandidate->Pt()-mother->Pt())/mother->Pt(),tempClusterWeight);
                if (fDoMesonQA > 1) fHistoTrueClusPrimPi01GammaMergedPtMCPt[fiCut]->Fill(TrueClusterCandidate->Pt(), mother->Pt(),tempClusterWeight);
              }else if (clusterClass == 3){
                fHistoTruePrimaryPi01ElectronMCPtResolPt[fiCut]->Fill(mother->Pt(),(TrueClusterCandidate->Pt()-mother->Pt())/mother->Pt(),tempClusterWeight);
                if (fDoMesonQA > 1) fHistoTrueClusPrimPi01ElectronMergedPtMCPt[fiCut]->Fill(TrueClusterCandidate->Pt(), mother->Pt(),tempClusterWeight);
              }
            }
          } else {
            fHistoTrueClusSecPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
            if (fDoMesonQA > 0 && mother->Pt()>0){
              fHistoTrueSecondaryPi0MCPtResolPt[fiCut]->Fill(mother->Pt(),(TrueClusterCandidate->Pt()-mother->Pt())/mother->Pt(),tempClusterWeight);
            }
            Int_t grandMaLab = mother->GetMother(0);
            if (grandMaLab > -1){
              if (TMath::Abs(fMCEvent->Particle(grandMaLab)->GetPdgCode()) == 310){
                fHistoTrueClusSecPi0FromK0sPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
              } else if (TMath::Abs(fMCEvent->Particle(grandMaLab)->GetPdgCode()) == 130){
                fHistoTrueClusSecPi0FromK0lPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
              } else if (TMath::Abs(fMCEvent->Particle(grandMaLab)->GetPdgCode()) == 3122){
                fHistoTrueClusSecPi0FromLambdaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
              }
            }
          }
        }
      // deal with eta only
      } else if (motherPDG == 221){
          fHistoTrueClusEtaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
          if(fDoJetAnalysis && !fDoLightOutput) fHistoTrueClusEtaPtvsRJet[fiCut]->Fill(TrueClusterCandidate->Pt(), RJetPi0Cand, tempClusterWeight);
          if (TrueClusterCandidate->IsDalitz()){
            fHistoTrueClusEtaDalitzPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
          }
          if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,motherLab)) fHistoDoubleCountTrueEtaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
          if (fDoMesonQA > 0){
            fHistoTrueClusEtaEM02[fiCut]->Fill(TrueClusterCandidate->E(), m02, tempClusterWeight);
            if (fDoMesonQA > 1){
              fHistoTrueClusEtaInvMassvsPt[fiCut]->Fill(mesoncand->M(),TrueClusterCandidate->Pt(), tempClusterWeight);
              if(TrueClusterCandidate->IsMerged()){
                fHistoTrueClusMergedEtaEVsM02[fiCut]->Fill(TrueClusterCandidate->E(), m02, tempClusterWeight);
                fHistoTrueClusMergedEtaEVsM20[fiCut]->Fill(TrueClusterCandidate->E(), cluster->GetM20(), tempClusterWeight);
              }
            }
          }
          if ( fDoMesonQA > 0 && GetSelectedMesonID() != 1 ){
            fHistoTrueEtaPtY[fiCut]->Fill(TrueClusterCandidate->Pt(),mesoncand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), tempClusterWeight);
            fHistoTrueEtaPtAlpha[fiCut]->Fill(TrueClusterCandidate->Pt(),TMath::Abs(mesoncand->GetAlpha()), tempClusterWeight);
            if( mother->Pt()>0)
              fHistoTruePrimaryEtaMCPtResolPt[fiCut]->Fill(mother->Pt(),(TrueClusterCandidate->Pt()-mother->Pt())/mother->Pt(),tempClusterWeight);
          }
      } else {
        fHistoTrueMergedMissedPDG[fiCut]->Fill(motherPDG, tempClusterWeight);
      }

    // leading particle is a photon or the conversion is fully contained and its not from pi0 || eta
    } else if (TrueClusterCandidate->IsLargestComponentPhoton() || TrueClusterCandidate->IsConversionFullyContained()
               || TrueClusterCandidate->IsElectronFromFragPhoton()){
      if (fEnableDetailedPrintOut) cout << "photon" << endl;
      fHistoTrueClusGammaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);

      if (fDoMesonQA > 0){
        fHistoTrueClusGammaEM02[fiCut]->Fill(TrueClusterCandidate->E(), m02, tempClusterWeight);
        if (fDoMesonQA > 1){
          fHistoTrueClusGammaEvsM20[fiCut]->Fill(TrueClusterCandidate->E(), cluster->GetM20(), tempClusterWeight);
          fHistoTrueClusGammaInvMassvsPt[fiCut]->Fill(mesoncand->M(),TrueClusterCandidate->Pt(), tempClusterWeight);
        }
      }

      if (motherLab == -1){
        fHistoTrueClusGammaPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 0.5, tempClusterWeight); // direct photon
      } else {
        if (motherPDG == 111)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 1.5, tempClusterWeight); // pi0
        else if (motherPDG == 221)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 2.5, tempClusterWeight); // eta
        else if (motherPDG == 331)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 3.5, tempClusterWeight); // eta'
        else if (motherPDG == 223)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 4.5, tempClusterWeight); // omega
        else if (motherPDG == 333)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 5.5, tempClusterWeight); // phi
        else if (motherPDG == 3122)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 6.5, tempClusterWeight); // Lambda
        else
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 7.5, tempClusterWeight); // rest
      }

    // leading particle is an electron and its not from pi0 || eta and no electron from fragmentation photon conversion
    } else if (TrueClusterCandidate->IsLargestComponentElectron() &&
               !TrueClusterCandidate->IsElectronFromFragPhoton()){
      if (fEnableDetailedPrintOut) cout << "electron" << endl;
      fHistoTrueClusElectronPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
      if (fDoMesonQA > 0){
        fHistoTrueClusElectronEM02[fiCut]->Fill(TrueClusterCandidate->E(), m02, tempClusterWeight);
        if (fDoMesonQA > 1){
          fHistoTrueClusElectronEM20[fiCut]->Fill(TrueClusterCandidate->E(), cluster->GetM20(), tempClusterWeight);
          fHistoTrueClusElectronInvMassvsPt[fiCut]->Fill(mesoncand->M(),TrueClusterCandidate->Pt(), tempClusterWeight);
        }
      }

      Int_t motherLab = Photon->GetMother(0);
      mother = fMCEvent->Particle(motherLab);
      if (mother){
        TrueClusterCandidate->GetLabel(0);
        motherPDG = TMath::Abs(mother->GetPdgCode());
      }

      if (motherLab == -1){
        fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 0.5, tempClusterWeight); // direct electron
      } else {
        if (motherPDG == 22){
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 1.5, tempClusterWeight); // gamma
          fHistoTrueClusElectronFromGammaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
        } else if (motherPDG == 111){
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 2.5, tempClusterWeight); // pi0
        } else if (motherPDG == 221){
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 3.5, tempClusterWeight); // eta
        } else if ( int(motherPDG/100.)==5 || int(motherPDG/1000.)==5 ){
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 5.5, tempClusterWeight); // b
        } else if ( int(motherPDG/100.)==4 || int(motherPDG/1000.)==4 ){
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 4.5, tempClusterWeight); // c
        } else if (motherPDG == 23 || motherPDG == 24){
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 6.5, tempClusterWeight); // W/Z
        } else if (motherPDG == 15) {
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 7.5, tempClusterWeight); // tau
        } else {
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 8.5, tempClusterWeight); // rest
        }
      }

    // leading particle is a hadron
    } else {
      if (fEnableDetailedPrintOut) cout << "BG" << endl;
      fHistoTrueClusBGPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
      if(fDoMesonQA > 1){
        fHistoTrueClusBGEvsM02[fiCut]->Fill(TrueClusterCandidate->E(), m02, tempClusterWeight);
        fHistoTrueClusBGEvsM20[fiCut]->Fill(TrueClusterCandidate->E(), cluster->GetM20(), tempClusterWeight);
        fHistoTrueClusBGInvMassvsPt[fiCut]->Fill(mesoncand->M(),TrueClusterCandidate->Pt(), tempClusterWeight);
      }
      Double_t maxM02 = 4.8;
      if (m02 >= 0 && m02 < maxM02){
        if (TMath::Abs(pdgCodeParticle) == 211) fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 0.5, tempClusterWeight); // pi+/-
        else if (TMath::Abs(pdgCodeParticle) == 2212) fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 1.5, tempClusterWeight); // p
        else if (TMath::Abs(pdgCodeParticle) == 321) fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 2.5, tempClusterWeight); // K+-
        else if (TMath::Abs(pdgCodeParticle) == 2112) fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 3.5, tempClusterWeight); // n
        else if (TMath::Abs(pdgCodeParticle) == 310) fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 4.5, tempClusterWeight); // K0s
        else if (TMath::Abs(pdgCodeParticle) == 3122) fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 5.5, tempClusterWeight); // Lambda
        else if (TMath::Abs(pdgCodeParticle) == 13) fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 6.5, tempClusterWeight); // mu+/-
        else if (TMath::Abs(pdgCodeParticle) == 130) fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 7.5, tempClusterWeight); // K0l
        else fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 8.5, tempClusterWeight); // Rest
      }
    }
  }
  delete mesoncand;
  return;
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::ProcessTrueClusterCandidatesAOD(AliAODConversionPhoton *TrueClusterCandidate, AliVCluster* cluster,
                                    AliAODConversionPhoton *TrueSubClusterCandidate1,
                                    AliAODConversionPhoton *TrueSubClusterCandidate2)
{
  Float_t m02 = cluster->GetM02();
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  Double_t tempClusterWeight       = fWeightJetJetMC;
  AliAODMCParticle *Photon = NULL;
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray){
    if (TrueClusterCandidate->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
    if (TrueClusterCandidate->GetCaloPhotonMCLabel(0) < 0) return;
  } else {
    AliInfo("fAODMCTrackArray could not be loaded");
    return;
  }

  AliAODConversionMother *mesoncand = new AliAODConversionMother(TrueSubClusterCandidate1,TrueSubClusterCandidate2);
  Bool_t mesonIsSelected            = (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(mesoncand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()));
  if (!mesonIsSelected){
    delete mesoncand;
    return;
  }

  // Setting all MC Flags (IsMerged, etc)
  TrueClusterCandidate->SetCaloPhotonMCFlagsAOD(fAODMCTrackArray, fEnableSortForClusMC, kTRUE,cluster);

  if (TrueClusterCandidate->GetNCaloPhotonMCLabels()>0 && fEnableDetailedPrintOut && fDoMesonQA > 1){
    // cout << endl << "Cluster energy: " << TrueClusterCandidate->E() << endl;;
    // TrueClusterCandidate->PrintCaloMCLabelsAndInfo(fMCEvent);
    // TrueClusterCandidate->PrintCaloMCLabelsAndInfoAOD(fInputEvent);
      PrintCaloMCLabelsAndInfoAOD(fInputEvent,TrueClusterCandidate,cluster);
    // TrueClusterCandidate->PrintCaloMCFlags();
  }

  if (TrueClusterCandidate->GetNCaloPhotonMCLabels()>0){
    fHistoTrueClusEFracFirstLabel[fiCut]->Fill(cluster->E(), cluster->GetClusterMCEdepFraction(0), tempClusterWeight);
    if(TrueClusterCandidate->GetNNeutralPionMCLabels()>0)
      fHistoTrueClusEFracLeadingPi0[fiCut]->Fill(cluster->E(), TrueClusterCandidate->GetNeutralPionEnergyFraction(TrueClusterCandidate->GetLeadingNeutralPionIndex()), tempClusterWeight);
    // check if leading pi0 comes not from label 0 in cluster
    // for this do:
    // -> check if neutral pions were found in cluster
    // -> if the leading daughter index is not 0
    // -> the leading neutral pion has a larger cluster energy fraction than the cluster label 0
    if( TrueClusterCandidate->GetNNeutralPionMCLabels()>0 && TrueClusterCandidate->GetLeadingNeutralPionDaughterIndex()!=0 && TrueClusterCandidate->GetNeutralPionEnergyFraction(TrueClusterCandidate->GetLeadingNeutralPionIndex())>cluster->GetClusterMCEdepFraction(0)) {
      // load particle corresponding to largest daughter of leading pi0
      Photon         = (AliAODMCParticle*) fAODMCTrackArray->At(TrueClusterCandidate->GetCaloPhotonMCLabel(TrueClusterCandidate->GetLeadingNeutralPionDaughterIndex()));
      // fill vector with leading pi0 MC label
      CheckVectorForDoubleCount(fVectorLabelsLeadingPi0,TrueClusterCandidate->GetNeutralPionMCLabel(TrueClusterCandidate->GetLeadingNeutralPionIndex()));
      // fill vector will ALL pi0 MC labels that are found in the cluster
      for(Int_t npion = 0; npion < TrueClusterCandidate->GetNNeutralPionMCLabels(); npion++){
        CheckVectorForDoubleCount(fVectorLabelsMultiplePi0,TrueClusterCandidate->GetNeutralPionMCLabel(npion));
      }
    } else {
      // load particle corresponding to MC label 0 in cluster
      Photon        = (AliAODMCParticle*) fAODMCTrackArray->At(TrueClusterCandidate->GetCaloPhotonMCLabel(0));
      // important check that leading daughter corresponds to the cluster MC label 0
      if(TrueClusterCandidate->GetLeadingNeutralPionDaughterIndex()==0){
        // fill vector with leading pi0 MC label
        CheckVectorForDoubleCount(fVectorLabelsLeadingPi0,TrueClusterCandidate->GetNeutralPionMCLabel(TrueClusterCandidate->GetLeadingNeutralPionIndex()));
        // fill vector will ALL pi0 MC labels that are found in the cluster
        for(Int_t npion = 0; npion < TrueClusterCandidate->GetNNeutralPionMCLabels(); npion++){
          CheckVectorForDoubleCount(fVectorLabelsMultiplePi0,TrueClusterCandidate->GetNeutralPionMCLabel(npion));
        }
      }
    }
  } else {
    // return if there are no MC labels in the cluster
    return;
  }

  if(Photon == NULL){
    return;
  }

  // Get Distance to Jet
  Double_t RJetPi0Cand = 0.;
  if(fDoJetAnalysis && !fDoLightOutput){
    fVectorJetEta = fConvJetReader->GetVectorJetEta();
    fVectorJetPhi = fConvJetReader->GetVectorJetPhi();
    for(Int_t j=0; j<fConvJetReader->GetNJets(); j++){
      Double_t DeltaEta = fVectorJetEta.at(j)-Photon->Eta();
      Double_t DeltaPhi = abs(fVectorJetPhi.at(j)-Photon->Phi());
      if(DeltaPhi > TMath::Pi()) {
        DeltaPhi = 2*TMath::Pi() - DeltaPhi;
      }
      RJetPi0Cand = TMath::Sqrt(DeltaEta*DeltaEta+DeltaPhi*DeltaPhi);
      if(fConvJetReader->Get_Jet_Radius() > 0 ){
        if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
          break;
        }
      }
    }
  }

  Int_t pdgCodeParticle             = Photon->GetPdgCode();

  if(fIsFromMBHeader && !fIsOverlappingWithOtherHeader){
    Int_t clusterClass    = 0;
    Bool_t isPrimary      = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);

    // cluster classification:
    // 1    - nice merged cluster (2 gamma | contributions from 2 gamma) from pi0/eta
    // 2    - contribution from only 1 partner (1 gamma, 1 fully coverted gamma) from pi0/eta
    // 3    - contribution from part of 1 partner (1 electron) from pi0/eta
    Long_t motherLab = -1;
    if (TrueClusterCandidate->IsMerged() || TrueClusterCandidate->IsMergedPartConv()){
        clusterClass    = 1;
        motherLab       = TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0);
    } else if (TrueClusterCandidate->GetNCaloPhotonMotherMCLabels()> 0){
      if (TrueClusterCandidate->IsLargestComponentElectron() || TrueClusterCandidate->IsLargestComponentPhoton()){
        if (TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0) > -1 && (((AliAODMCParticle*) fAODMCTrackArray->At(TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0)))->GetPdgCode() == 111 || ((AliAODMCParticle*) fAODMCTrackArray->At(TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0)))->GetPdgCode() == 221) ){
          if ( TrueClusterCandidate->IsConversion() && !TrueClusterCandidate->IsConversionFullyContained() ){
            clusterClass  = 3;
            motherLab       = TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0);
          } else {
            clusterClass  = 2;
            motherLab       = TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0);
          }
        }
      } else if (TrueClusterCandidate->IsSubLeadingEM()){
        if (TrueClusterCandidate->GetNCaloPhotonMotherMCLabels()> 1){
          if (fEnableDetailedPrintOut) cout << "Is Subleading EM: "<<  TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1) << endl;
          if ( TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1) > -1){
            if (TMath::Abs(((AliAODMCParticle*) fAODMCTrackArray->At(TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1)))->GetPdgCode()) == 111 || TMath::Abs(((AliAODMCParticle*) fAODMCTrackArray->At(TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1)))->GetPdgCode()) == 221 ){
              clusterClass  = 2;
              motherLab       = TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1);
            }
          }
        }
      } else {
        motherLab       = TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0);
      }
    }

    // Get Mother particle
    AliAODMCParticle *mother = NULL;
    Int_t motherPDG   = -1;
    if (motherLab > -1){
       mother           = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(motherLab));
    }
    if (fEnableDetailedPrintOut) cout << "cluster class: " << clusterClass << "\t mother lab: "<< motherLab ;
    if (mother){
        motherPDG = TMath::Abs(mother->GetPdgCode());
        if (fEnableDetailedPrintOut) cout  << "\t mother pdg: " << motherPDG << endl;
    } else {
      if (fEnableDetailedPrintOut) cout << endl;
    }

    // Set the jetjet weight to 1 in case the cluster orignated from the minimum bias header
    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
      if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(motherLab, fMCEvent, fInputEvent) == 2) tempClusterWeight = 1;
    }

    //
    if (clusterClass == 1 || clusterClass == 2 || clusterClass == 3 ){
      fHistoTrueClusMergedPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
      if(fDoJetAnalysis && !fDoLightOutput) fHistoTrueClusMergedPtvsRJet[fiCut]->Fill(TrueClusterCandidate->Pt(), RJetPi0Cand, tempClusterWeight);
      if (fDoMesonQA > 1)fHistoTrueClusMergedInvMassvsPt[fiCut]->Fill(mesoncand->M(),TrueClusterCandidate->Pt(), tempClusterWeight);

      // separate different components
      if (clusterClass == 1 && TrueClusterCandidate->IsMerged()){
        if (motherPDG == 111){
          fHistoTrueClusMergedPureFromPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
          if (GetSelectedMesonID() < 2 && !isPrimary && m02 >= 0 && m02 <= 4.8 ){
            fHistoTrueSecPi0PtvsDiffReco[fiCut]->Fill(TrueClusterCandidate->Pt(), 0., tempClusterWeight);
          }
        }
        if (motherPDG == 221)
          fHistoTrueClusMergedPureFromEtaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
      } else if (clusterClass == 1 && TrueClusterCandidate->IsMergedPartConv()){
        if (motherPDG == 111){
          fHistoTrueClusMergedPartConvFromPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
          if (GetSelectedMesonID() < 2 && !isPrimary && m02 >= 0 && m02 <= 4.8 ){
            fHistoTrueSecPi0PtvsDiffReco[fiCut]->Fill(TrueClusterCandidate->Pt(), 1., tempClusterWeight);
          }
        }
        if (motherPDG == 221)
          fHistoTrueClusMergedPartConvFromEtaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
      } else if (clusterClass == 2){
        if (motherPDG == 111){
          fHistoTrueClusGammaFromPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
          if (GetSelectedMesonID() < 2 && !isPrimary && m02 >= 0 && m02 <= 4.8 ){
            fHistoTrueSecPi0PtvsDiffReco[fiCut]->Fill(TrueClusterCandidate->Pt(), 2., tempClusterWeight);
          }
        }if (motherPDG == 221)
          fHistoTrueClusGammaFromEtaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
      } else if (clusterClass == 3){
        if (motherPDG == 111) {
          fHistoTrueClusElectronFromPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
          if (GetSelectedMesonID() < 2 && !isPrimary && m02 >= 0 && m02 <= 4.8 ){
            fHistoTrueSecPi0PtvsDiffReco[fiCut]->Fill(TrueClusterCandidate->Pt(), 3., tempClusterWeight);
          }
        }
        if (motherPDG == 221)
          fHistoTrueClusElectronFromEtaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
      }
      // deal with pi0 only
      if (motherPDG == 111){
        fHistoTrueClusPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
        if(fDoJetAnalysis && !fDoLightOutput) fHistoTrueClusPi0PtvsRJet[fiCut]->Fill(TrueClusterCandidate->Pt(), RJetPi0Cand, tempClusterWeight);
        if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,motherLab)){
          fHistoDoubleCountTruePi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
          if (!isPrimary)
            fHistoDoubleCountTrueSecPi0Pt[fiCut]->Fill(TrueClusterCandidate->Pt(), tempClusterWeight);
        }
        // Treatment of multiple pi0 in cluster (filling of true candidates, amount of pions and double counting)
        if(TrueClusterCandidate->GetNNeutralPionMCLabels()>1){
          for(Int_t npion = 0; npion < TrueClusterCandidate->GetNNeutralPionMCLabels(); npion++){
            // fill multiple true pions in same cluster
            if(npion>0) fHistoTrueClusMultiplePi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);

            // check MC labels for double counting (each particle should only be reconstructed once)
            if (CheckVectorForDoubleCount(fVectorDoubleCountTrueMultilePi0s,TrueClusterCandidate->GetNeutralPionMCLabel(npion))){
              fHistoDoubleCountTrueMultiplePi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);

              AliAODMCParticle *PhotonDummyPrimary =  (AliAODMCParticle*) fAODMCTrackArray->At(TrueClusterCandidate->GetNeutralPionMCLabel(npion));
              // check if particle is not a primary -> secondary particle
              if(!((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, PhotonDummyPrimary, mcProdVtxX, mcProdVtxY, mcProdVtxZ))
                fHistoDoubleCountTrueMultipleSecPi0Pt[fiCut]->Fill(TrueClusterCandidate->Pt(), tempClusterWeight);
            }
          }
        }
        fHistoNTrueMultiplePi0vsPt[fiCut]->Fill(TrueClusterCandidate->Pt(),TrueClusterCandidate->GetNNeutralPionMCLabels(), tempClusterWeight);

        if (TrueClusterCandidate->IsDalitz()){
          fHistoTrueClusPi0DalitzPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
        }
        if (fDoMesonQA > 0){
          fHistoTrueClusPi0EM02[fiCut]->Fill(TrueClusterCandidate->E(), m02, tempClusterWeight);
          if (fDoMesonQA > 1){
            fHistoTrueClusPi0InvMassvsPt[fiCut]->Fill(mesoncand->M(),TrueClusterCandidate->Pt(), tempClusterWeight);
          }
        }
        if (fDoMesonQA > 0 && GetSelectedMesonID() != 2){
          fHistoTruePi0PtY[fiCut]->Fill(TrueClusterCandidate->Pt(),mesoncand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), tempClusterWeight);
          fHistoTruePi0PtAlpha[fiCut]->Fill(TrueClusterCandidate->Pt(),TMath::Abs(mesoncand->GetAlpha()), tempClusterWeight);
        }

        if (GetSelectedMesonID() < 2) {
          Float_t weighted= 1;
          if (mother->Pt()>0.005){
            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(motherLab, fMCEvent, fInputEvent);
          }
          if (isPrimary) {
            fHistoTrueClusPrimPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight*weighted);
            if (clusterClass == 1 && TrueClusterCandidate->IsMerged()){
              fHistoTruePureMergedClusPrimPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight*weighted);
              if (fDoMesonQA > 1){
                fHistoTrueClusMergedPi0EVsM20[fiCut]->Fill(TrueClusterCandidate->E(), cluster->GetM20(), tempClusterWeight*weighted);
                fHistoTrueClusMergedPi0EVsM02[fiCut]->Fill(TrueClusterCandidate->E(), m02, tempClusterWeight*weighted);
              }
            }else if (clusterClass == 1 && TrueClusterCandidate->IsMergedPartConv()){
              fHistoTruePartConvClusPrimPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight*weighted);
              if (fDoMesonQA > 1){
                fHistoTrueClusMergedPartConvPi0EVsM20[fiCut]->Fill(TrueClusterCandidate->E(), cluster->GetM20(), tempClusterWeight);
              }
            }
            if(TrueClusterCandidate->GetNNeutralPionMCLabels()>1){
              for(Int_t npion = 1; npion < TrueClusterCandidate->GetNNeutralPionMCLabels(); npion++){
                AliAODMCParticle *PhotonDummyPrimary =  (AliAODMCParticle*) fAODMCTrackArray->At(TrueClusterCandidate->GetNeutralPionMCLabel(npion));
                if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, PhotonDummyPrimary, mcProdVtxX, mcProdVtxY, mcProdVtxZ)) // check if pion is primary
                  fHistoTrueClusMultiplePrimPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
              }
            }
            if (fDoMesonQA > 1) fHistoTrueClusPrimPi0InvMassvsPt[fiCut]->Fill(mesoncand->M(),TrueClusterCandidate->Pt(), tempClusterWeight);
            if (fDoMesonQA > 0 && mother->Pt()>0){
              fHistoTruePrimaryPi0MCPtResolPt[fiCut]->Fill(mother->Pt(),(TrueClusterCandidate->Pt()-mother->Pt())/mother->Pt(),tempClusterWeight);
                if (fDoMesonQA > 1) fHistoTrueClusPrimPi0PtMCPt[fiCut]->Fill(TrueClusterCandidate->Pt(), mother->Pt(),tempClusterWeight);
              if (clusterClass == 1 && TrueClusterCandidate->IsMerged()){
                fHistoTruePrimaryPi0PureMergedMCPtResolPt[fiCut]->Fill(mother->Pt(),(TrueClusterCandidate->Pt()-mother->Pt())/mother->Pt(),tempClusterWeight);
                if (fDoMesonQA > 1) fHistoTrueClusPrimPi0PureMergedPtMCPt[fiCut]->Fill(TrueClusterCandidate->Pt(), mother->Pt(),tempClusterWeight);
              }else if (clusterClass == 1 && TrueClusterCandidate->IsMergedPartConv()){
                fHistoTruePrimaryPi0MergedPartConvMCPtResolPt[fiCut]->Fill(mother->Pt(),(TrueClusterCandidate->Pt()-mother->Pt())/mother->Pt(),tempClusterWeight);
                if (fDoMesonQA > 1) fHistoTrueClusPrimPi0PartConvMergedPtMCPt[fiCut]->Fill(TrueClusterCandidate->Pt(), mother->Pt(),tempClusterWeight);
              }else if (clusterClass == 2){
                fHistoTruePrimaryPi01GammaMCPtResolPt[fiCut]->Fill(mother->Pt(),(TrueClusterCandidate->Pt()-mother->Pt())/mother->Pt(),tempClusterWeight);
                if (fDoMesonQA > 1) fHistoTrueClusPrimPi01GammaMergedPtMCPt[fiCut]->Fill(TrueClusterCandidate->Pt(), mother->Pt(),tempClusterWeight);
              }else if (clusterClass == 3){
                fHistoTruePrimaryPi01ElectronMCPtResolPt[fiCut]->Fill(mother->Pt(),(TrueClusterCandidate->Pt()-mother->Pt())/mother->Pt(),tempClusterWeight);
                if (fDoMesonQA > 1) fHistoTrueClusPrimPi01ElectronMergedPtMCPt[fiCut]->Fill(TrueClusterCandidate->Pt(), mother->Pt(),tempClusterWeight);
              }
            }
          } else {
            fHistoTrueClusSecPi0PtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
            if (fDoMesonQA > 0 && mother->Pt()>0){
              fHistoTrueSecondaryPi0MCPtResolPt[fiCut]->Fill(mother->Pt(),(TrueClusterCandidate->Pt()-mother->Pt())/mother->Pt(),tempClusterWeight);
            }
            Int_t grandMaLab = mother->GetMother();
            if (grandMaLab > -1){
              if (TMath::Abs(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(grandMaLab))->GetPdgCode()) == 310){
                fHistoTrueClusSecPi0FromK0sPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
              } else if (TMath::Abs(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(grandMaLab))->GetPdgCode()) == 130){
                fHistoTrueClusSecPi0FromK0lPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
              } else if (TMath::Abs(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(grandMaLab))->GetPdgCode()) == 3122){
                fHistoTrueClusSecPi0FromLambdaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
              }
            }
          }
        }
      // deal with eta only
      } else if (motherPDG == 221){
          fHistoTrueClusEtaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
          if(fDoJetAnalysis && !fDoLightOutput) fHistoTrueClusEtaPtvsRJet[fiCut]->Fill(TrueClusterCandidate->Pt(), RJetPi0Cand, tempClusterWeight);
          if (TrueClusterCandidate->IsDalitz()){
            fHistoTrueClusEtaDalitzPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
          }
          if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,motherLab)) fHistoDoubleCountTrueEtaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
          if (fDoMesonQA > 0){
            fHistoTrueClusEtaEM02[fiCut]->Fill(TrueClusterCandidate->E(), m02, tempClusterWeight);
            if (fDoMesonQA > 1){
              fHistoTrueClusEtaInvMassvsPt[fiCut]->Fill(mesoncand->M(),TrueClusterCandidate->Pt(), tempClusterWeight);
              if(TrueClusterCandidate->IsMerged()){
                fHistoTrueClusMergedEtaEVsM02[fiCut]->Fill(TrueClusterCandidate->E(), m02, tempClusterWeight);
                fHistoTrueClusMergedEtaEVsM20[fiCut]->Fill(TrueClusterCandidate->E(), cluster->GetM20(), tempClusterWeight);
              }
            }
          }
          if ( fDoMesonQA > 0 && GetSelectedMesonID() != 1 ){
            fHistoTrueEtaPtY[fiCut]->Fill(TrueClusterCandidate->Pt(),mesoncand->Rapidity()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), tempClusterWeight);
            fHistoTrueEtaPtAlpha[fiCut]->Fill(TrueClusterCandidate->Pt(),TMath::Abs(mesoncand->GetAlpha()), tempClusterWeight);
            if(mother->Pt()>0)
              fHistoTruePrimaryEtaMCPtResolPt[fiCut]->Fill(mother->Pt(),(TrueClusterCandidate->Pt()-mother->Pt())/mother->Pt(),tempClusterWeight);
          }
      } else {
        fHistoTrueMergedMissedPDG[fiCut]->Fill(motherPDG, tempClusterWeight);
      }
    // leading particle is a photon or the conversion is fully contained and its not from pi0 || eta
    } else if (TrueClusterCandidate->IsLargestComponentPhoton() || TrueClusterCandidate->IsConversionFullyContained()){
      if (fEnableDetailedPrintOut) cout << "photon" << endl;
      fHistoTrueClusGammaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
      if (fDoMesonQA > 0){
        fHistoTrueClusGammaEM02[fiCut]->Fill(TrueClusterCandidate->E(), m02, tempClusterWeight);
        if (fDoMesonQA > 1) {
          fHistoTrueClusGammaEvsM20[fiCut]->Fill(TrueClusterCandidate->E(), cluster->GetM20(), tempClusterWeight);
          fHistoTrueClusGammaInvMassvsPt[fiCut]->Fill(mesoncand->M(),TrueClusterCandidate->Pt(), tempClusterWeight);
        }
    }
      if (motherLab == -1){
        fHistoTrueClusGammaPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 0.5, tempClusterWeight); // direct photon
      } else {
        if (motherPDG == 111)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 1.5, tempClusterWeight); // pi0
        else if (motherPDG == 221)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 2.5, tempClusterWeight); // eta
        else if (motherPDG == 331)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 3.5, tempClusterWeight); // eta'
        else if (motherPDG == 223)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 4.5, tempClusterWeight); // omega
        else if (motherPDG == 333)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 5.5, tempClusterWeight); // phi
        else if (motherPDG == 3122)
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 6.5, tempClusterWeight); // Lambda
        else
          fHistoTrueClusGammaPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 7.5, tempClusterWeight); // rest
      }
    // leading particle is an electron and its not from pi0 || eta
    } else if (TrueClusterCandidate->IsLargestComponentElectron()){
      if (fEnableDetailedPrintOut) cout << "electron" << endl;
      fHistoTrueClusElectronPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
      if (fDoMesonQA > 0){
        fHistoTrueClusElectronEM02[fiCut]->Fill(TrueClusterCandidate->E(), m02, tempClusterWeight);
        if (fDoMesonQA > 1){
          fHistoTrueClusElectronEM20[fiCut]->Fill(TrueClusterCandidate->E(), cluster->GetM20(), tempClusterWeight);
          fHistoTrueClusElectronInvMassvsPt[fiCut]->Fill(mesoncand->M(),TrueClusterCandidate->Pt(), tempClusterWeight);
        }
      }

      Int_t motherLab = Photon->GetMother();
      if (motherLab == -1){
        fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 0.5, tempClusterWeight); // direct electron
      } else {
        if (motherPDG == 22){
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 1.5, tempClusterWeight); // gamma
          fHistoTrueClusElectronFromGammaPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
        } else if (motherPDG == 111){
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 2.5, tempClusterWeight); // pi0
        } else if (motherPDG == 221){
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 3.5, tempClusterWeight); // eta
        } else if ( int(motherPDG/100.)==5 || int(motherPDG/1000.)==5 ){
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 5.5, tempClusterWeight); // b
        } else if ( int(motherPDG/100.)==4 || int(motherPDG/1000.)==4 ){
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 4.5, tempClusterWeight); // c
        } else if (motherPDG == 23 || motherPDG == 24){
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 6.5, tempClusterWeight); // W/Z
        } else if (motherPDG == 15) {
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 7.5, tempClusterWeight); // tau
        } else {
          fHistoTrueClusElectronPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 8.5, tempClusterWeight); // rest
        }
      }
    // leading particle is a hadron
    } else {
      if (fEnableDetailedPrintOut) cout << "BG" << endl;
      fHistoTrueClusBGPtvsM02[fiCut]->Fill(TrueClusterCandidate->Pt(), m02, tempClusterWeight);
      if(fDoMesonQA > 1){
        fHistoTrueClusBGEvsM02[fiCut]->Fill(TrueClusterCandidate->E(), m02, tempClusterWeight);
        fHistoTrueClusBGEvsM20[fiCut]->Fill(TrueClusterCandidate->E(), cluster->GetM20(), tempClusterWeight);
        fHistoTrueClusBGInvMassvsPt[fiCut]->Fill(mesoncand->M(),TrueClusterCandidate->Pt(), tempClusterWeight);
      }

      Double_t maxM02 = 4.8;
      if (m02 >= 0 && m02 < maxM02){
        if (TMath::Abs(pdgCodeParticle) == 211) fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 0.5, tempClusterWeight); // pi+/-
        else if (TMath::Abs(pdgCodeParticle) == 2212) fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 1.5, tempClusterWeight); // p
        else if (TMath::Abs(pdgCodeParticle) == 321) fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 2.5, tempClusterWeight); // K+-
        else if (TMath::Abs(pdgCodeParticle) == 2112) fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 3.5, tempClusterWeight); // n
        else if (TMath::Abs(pdgCodeParticle) == 310) fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 4.5, tempClusterWeight); // K0s
        else if (TMath::Abs(pdgCodeParticle) == 3122) fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 5.5, tempClusterWeight); // Lambda
        else if (TMath::Abs(pdgCodeParticle) == 13) fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 6.5, tempClusterWeight); // mu+/-
        else if (TMath::Abs(pdgCodeParticle) == 130) fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 7.5, tempClusterWeight); // K0l
        else fHistoTrueClusBGPtvsSource[fiCut]->Fill(TrueClusterCandidate->Pt(), 8.5, tempClusterWeight); // Rest
      }
    }
  }
  delete mesoncand;
  return;
}


//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::ProcessMCParticles()
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  Bool_t particleInJet = kFALSE;
  if(fDoJetAnalysis){
    if(fConvJetReader->GetTrueNJets()>0){
      fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
      fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
    }
  }

  // Loop over all primary MC particle
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    Double_t tempParticleWeight       = fWeightJetJetMC;
    if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){

      TParticle* particle = (TParticle *)fMCEvent->Particle(i);
      if (!particle) continue;

      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
        // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
        if(isMCFromMBHeader == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempParticleWeight = 1;
      }

      // check if particle is in Jet
      if(fDoJetAnalysis){
        if(fDoOutOfJet == 1) particleInJet = kTRUE;
        else particleInJet = kFALSE;
        for(Int_t j=0; j<fConvJetReader->GetTrueNJets(); j++){
          Double_t DeltaEta = fVectorJetEta.at(j)-particle->Eta();
          Double_t DeltaPhi = abs(fVectorJetPhi.at(j)-particle->Phi());
          if(fDoOutOfJet == 2){ // check if on opposite side of jet (DeltaEta/Phi = 0 if directly opposite)
            DeltaEta = fVectorJetEta.at(j) + particle->Eta();
            DeltaPhi = abs(TMath::Pi() - DeltaPhi);
          }
          if(DeltaPhi > TMath::Pi()) {
            DeltaPhi = 2*TMath::Pi() - DeltaPhi;
          }
          Double_t RJetPi0Cand = TMath::Sqrt(DeltaEta*DeltaEta+DeltaPhi*DeltaPhi);
          if(fConvJetReader->Get_Jet_Radius() > 0 ){
            if(fDoOutOfJet == 0){ // in jet
              if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                particleInJet = kTRUE;
                break;
              }
            } else if(fDoOutOfJet == 1){ // out of jet
              if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                particleInJet = kFALSE;
                break;
              }
            } else if(fDoOutOfJet == 2){ // out of jet on away side
              if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                particleInJet = kTRUE;
                break;
              }
            } else if(fDoOutOfJet == 3){ // out of jet in interval [R, R+0.2]
              if((RJetPi0Cand > fConvJetReader->Get_Jet_Radius()) && (RJetPi0Cand < fConvJetReader->Get_Jet_Radius() + 0.2)){
                particleInJet = kTRUE;
                break;
              }
            }
          }
        }
      }

      Double_t mesonY = 1.e30;
      Double_t ratio  = 0;
      if (particle->Energy() != TMath::Abs(particle->Pz())){
        ratio         = (particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz());
      }
      if( !(ratio <= 0) ){
        mesonY = particle->Y()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
      }

      // fill Primary Y hist
      if ( particle->GetPdgCode() == 211 ){  // positve pions
        fHistoMCPrimaryYvsSource[fiCut]->Fill(mesonY, 0., tempParticleWeight);
      } else if ( particle->GetPdgCode() == -211 ){  // negative pions
        fHistoMCPrimaryYvsSource[fiCut]->Fill(mesonY, 1., tempParticleWeight);
      } else if ( particle->GetPdgCode() == 321 ){  // positve kaons
        fHistoMCPrimaryYvsSource[fiCut]->Fill(mesonY, 2., tempParticleWeight);
      } else if ( particle->GetPdgCode() == -321 ){  // negative kaons
        fHistoMCPrimaryYvsSource[fiCut]->Fill(mesonY, 3., tempParticleWeight);
      } else if ( TMath::Abs(particle->GetPdgCode()) == 310 ){  // K0s
        fHistoMCPrimaryYvsSource[fiCut]->Fill(mesonY, 4., tempParticleWeight);
      } else if ( TMath::Abs(particle->GetPdgCode()) == 130 ){  // K0l
        fHistoMCPrimaryYvsSource[fiCut]->Fill(mesonY, 5., tempParticleWeight);
      } else if ( TMath::Abs(particle->GetPdgCode()) == 3122 ){  // Lambda/ AntiLambda
        fHistoMCPrimaryYvsSource[fiCut]->Fill(mesonY, 6., tempParticleWeight);
      }

      if ((mesonY > ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValueMin()) && (mesonY < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValueMax())){
        if ( particle->GetPdgCode() == 211 ){  // positve pions
          fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 0., tempParticleWeight);
        } else if ( particle->GetPdgCode() == -211 ){  // negative pions
          fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 1., tempParticleWeight);
        } else if ( particle->GetPdgCode() == 321 ){  // positve kaons
          fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 2., tempParticleWeight);
        } else if ( particle->GetPdgCode() == -321 ){  // negative kaons
          fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 3., tempParticleWeight);
        } else if ( TMath::Abs(particle->GetPdgCode()) == 310 ){  // K0s
          fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 4., tempParticleWeight);
        } else if ( TMath::Abs(particle->GetPdgCode()) == 130 ){  // K0l
          fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 5., tempParticleWeight);
        } else if ( TMath::Abs(particle->GetPdgCode()) == 3122 ){  // Lambda/ AntiLambda
          fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 6., tempParticleWeight);
        } else if ( particle->GetPdgCode() == 22 ){  // photons
          fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight); // direct photons
          if(particle->GetMother(0) > -1){
            TParticle* mother = (TParticle*)fMCEvent->Particle(particle->GetMother(0));
            if (  TMath::Abs(mother->GetPdgCode()) == 111  ||
                  TMath::Abs(mother->GetPdgCode()) == 113  ||
                  TMath::Abs(mother->GetPdgCode()) == 221  ||
                  TMath::Abs(mother->GetPdgCode()) == 223  ||
                  TMath::Abs(mother->GetPdgCode()) == 331  ||
                  TMath::Abs(mother->GetPdgCode()) == 333  ||
                  TMath::Abs(mother->GetPdgCode()) == 3212 ||
                  TMath::Abs(mother->GetPdgCode()) == 213
              ){
              fHistoMCDecayGammaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight); // decay photons
            }
          }
        }
      }
      // continue if in Jet analysis and cluster is not within Jet radius
      if(fDoJetAnalysis && !particleInJet) continue;

      // check if particle is pi0/eta from di-photon decay
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
        ->MesonIsSelectedMC(particle,fMCEvent,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
        TParticle* daughter0 = (TParticle*)fMCEvent->Particle(particle->GetFirstDaughter());
        TParticle* daughter1 = (TParticle*)fMCEvent->Particle(particle->GetLastDaughter());

        Float_t weighted= 1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
          if (particle->Pt()>0.005){
            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCEvent, fInputEvent);
          }
        }

        if( TMath::Abs(particle->GetPdgCode()) == 111 ){ // neutral pions
          fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // All MC Pi0
          if(std::find(fVectorLabelsMultiplePi0Reduced.begin(), fVectorLabelsMultiplePi0Reduced.end(), i) != fVectorLabelsMultiplePi0Reduced.end()) {
            fHistoMCPi0ReducedPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc, ONLY overlapping pi0 in cluster
          }
          if (GetSelectedMesonID() != 2){
            fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
            if (fIsMC==2)fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0 ){
              if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),tempParticleWeight);
            }
          }
        } else if( TMath::Abs(particle->GetPdgCode()) == 221 ){ // eta mesons
          fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // All MC Eta
          if (GetSelectedMesonID() != 1){
            fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
            if (fIsMC==2)fHistoMCEtaWOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0 ){
              if (fIsMC == 2) fHistoMCEtaPtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),tempParticleWeight);
            }
          }
        }
        // Check the acceptance for both gammas & whether they are counted as primaries as well
        Bool_t kDaughter0IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, particle->GetFirstDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        Bool_t kDaughter1IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, particle->GetLastDaughter(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);

        if( kDaughter0IsPrim && kDaughter1IsPrim &&
           (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCEvent) ||
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCEvent) ) ){
          if(particle->GetPdgCode() == 111 && GetSelectedMesonID() != 2){
            fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc
            if(std::find(fVectorLabelsMultiplePi0Reduced.begin(), fVectorLabelsMultiplePi0Reduced.end(), i) != fVectorLabelsMultiplePi0Reduced.end()) {
              fHistoMCPi0ReducedInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc, ONLY overlapping pi0 in cluster
            }
            if (fIsMC == 2)fHistoMCPi0WOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0 w/o event weights with gamma in acc
          } else if(particle->GetPdgCode() == 221 && GetSelectedMesonID() != 1){
            fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Eta with gamma in acc
            if (fIsMC == 2)fHistoMCEtaWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Eta w/o event weights with gamma in acc
          }
        }
      }
      Int_t gammaLabel    = -1;
      Int_t electronLabel = -1;
      Int_t positronLabel = -1;
      // check if particle is pi0/eta from Dalitz decay
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
        ->MesonIsSelectedMCDalitz(particle,fMCEvent, electronLabel, positronLabel, gammaLabel, ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
       if( (gammaLabel > -1) && (electronLabel > -1) && (positronLabel > -1) ){
        TParticle* gamma    = (TParticle*)fMCEvent->Particle(gammaLabel);
        TParticle* electron = (TParticle*)fMCEvent->Particle(electronLabel);
        TParticle* positron = (TParticle*)fMCEvent->Particle(positronLabel);

        Float_t weighted= 1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
          if (particle->Pt()>0.005){
            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCEvent, fInputEvent);
          }
        }

        if(particle->GetPdgCode() == 111 && GetSelectedMesonID() != 2){
          fHistoMCPi0DalitzPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // All MC Pi0
          fHistoMCPi0DalitzWOWeightPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
          if (fIsMC==2)fHistoMCPi0DalitzWOEvtWeightPt[fiCut]->Fill(particle->Pt());
          if (fDoMesonQA > 0 ){
            if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),tempParticleWeight);
          }
        } else if(particle->GetPdgCode() == 221 && GetSelectedMesonID() != 1){
          fHistoMCEtaDalitzPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // All MC Eta
          fHistoMCEtaDalitzWOWeightPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
          if (fIsMC==2)fHistoMCEtaDalitzWOEvtWeightPt[fiCut]->Fill(particle->Pt());
          if (fDoMesonQA > 0 ){
            if (fIsMC == 2) fHistoMCEtaPtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),tempParticleWeight);
          }
        }

        // Check the acceptance for both gammas & whether they are counted as primaries as well
        Bool_t kGammaIsPrim     = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, gammaLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        Bool_t kElectronIsPrim  = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, electronLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        Bool_t kPositronIsPrim  = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, positronLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        if( kGammaIsPrim && kElectronIsPrim && kPositronIsPrim &&
            (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma,fMCEvent) ||
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedElecMC(electron,fMCEvent) ||
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedElecMC(positron,fMCEvent) )
          ){
          if(particle->GetPdgCode() == 111 && GetSelectedMesonID() != 2){
            fHistoMCPi0DalitzInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc
            if (fIsMC == 2) fHistoMCPi0DalitzWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc
          } else if(particle->GetPdgCode() == 221 && GetSelectedMesonID() != 1){
            fHistoMCEtaDalitzInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Eta with gamma in acc
            if (fIsMC == 2) fHistoMCEtaDalitzWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc
          }
        }
       }
      }
    // End of primary threatment, now secondary treatment
    } else {

      TParticle* particle = (TParticle *)fMCEvent->Particle(i);
      if (!particle) continue;

      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
        // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
        if(isMCFromMBHeader == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempParticleWeight = 1;
      }

      // check if particle is pi0 from di-photon decay
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
        ->MesonIsSelectedMC(particle,fMCEvent,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
        TParticle* daughter0  = (TParticle*)fMCEvent->Particle(particle->GetFirstDaughter());
        TParticle* daughter1  = (TParticle*)fMCEvent->Particle(particle->GetLastDaughter());
        TParticle* mother     = NULL;
        Int_t motherPDG       = -1000000;
        if (particle->GetMother(0) > -1){
            mother            = (TParticle*)fMCEvent->Particle(particle->GetMother(0));
            if (mother)
              motherPDG       = TMath::Abs(mother->GetPdgCode());
        }

        if( TMath::Abs(particle->GetPdgCode()) == 111 ){ // neutral pions
          Int_t source        = GetSourceClassification(111,motherPDG);
          fHistoMCSecPi0PtvsSource[fiCut]->Fill(particle->Pt(),source, tempParticleWeight); // All secondary MC Pi0
        }

        // check whether pi0 landet in acceptance
        if( (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCEvent) ||
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCEvent) ) ){
          if( TMath::Abs(particle->GetPdgCode()) == 111 ){ // neutral pions
            Int_t source      = GetSourceClassification(111,motherPDG);
            fHistoMCSecPi0InAccPtvsSource[fiCut]->Fill(particle->Pt(),source, tempParticleWeight); // All secondary MC Pi0
          }
        }
      }

      Int_t gammaLabel    = -1;
      Int_t electronLabel = -1;
      Int_t positronLabel = -1;
      // check if particle is pi0/eta from Dalitz decay
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
        ->MesonIsSelectedMCDalitz(particle,fMCEvent, electronLabel, positronLabel, gammaLabel, ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
       if( (gammaLabel > -1) && (electronLabel > -1) && (positronLabel > -1) ){
        TParticle* gamma    = (TParticle*)fMCEvent->Particle(gammaLabel);
        TParticle* electron = (TParticle*)fMCEvent->Particle(electronLabel);
        TParticle* positron = (TParticle*)fMCEvent->Particle(positronLabel);

        TParticle* mother     = NULL;
        Int_t motherPDG       = -1000000;
        if (particle->GetMother(0) > -1){
            mother            = (TParticle*)fMCEvent->Particle(particle->GetMother(0));
            if (mother)
              motherPDG       = TMath::Abs(mother->GetPdgCode());
        }

        if( TMath::Abs(particle->GetPdgCode()) == 111 ){ // neutral pions
          Int_t source        = GetSourceClassification(111,motherPDG);
          fHistoMCSecPi0PtvsSource[fiCut]->Fill(particle->Pt(),source, tempParticleWeight); // All secondary MC Pi0
        }

        // check whether pi0 landet in acceptance
        if( (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(gamma,fMCEvent) ||
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedElecMC(electron,fMCEvent) ||
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedElecMC(positron,fMCEvent) ) ){
          if( TMath::Abs(particle->GetPdgCode()) == 111 ){ // neutral pions
            Int_t source      = GetSourceClassification(111,motherPDG);
            fHistoMCSecPi0InAccPtvsSource[fiCut]->Fill(particle->Pt(),source, tempParticleWeight); // All secondary MC Pi0
          }
        }
       }
      }
    }
  } // end of particle loop
}

//________________________________________________________________________
void AliAnalysisTaskGammaCaloMerged::ProcessAODMCParticles()
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray == NULL) return;

  // Load Jet eta phi vectors
  if(fDoJetAnalysis){
    if(fConvJetReader->GetTrueNJets()>0){
      fTrueVectorJetEta = fConvJetReader->GetTrueVectorJetEta();
      fTrueVectorJetPhi = fConvJetReader->GetTrueVectorJetPhi();
    }
  }

  // Loop over all primary MC particle
  for(Long_t i = 0; i < fAODMCTrackArray->GetEntriesFast(); i++) {
    Double_t tempParticleWeight       = fWeightJetJetMC;
    AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(i));
    if (!particle) continue;

    Bool_t isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, particle, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if (isPrimary) {

      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
        // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
        if(isMCFromMBHeader == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempParticleWeight = 1;
      }
      // check if particle is in Jet
      Bool_t particleInJet = kFALSE;
      if(fDoJetAnalysis){
        if(fDoOutOfJet == 1) particleInJet = kTRUE;
        else particleInJet = kFALSE;
        for(Int_t j=0; j<fConvJetReader->GetTrueNJets(); j++){
          Double_t DeltaEta = fTrueVectorJetEta.at(j)-particle->Eta();
          Double_t DeltaPhi = abs(fTrueVectorJetPhi.at(j)-particle->Phi());
          if(fDoOutOfJet == 2){ // check if on opposite side of jet (DeltaEta/Phi = 0 if directly opposite)
            DeltaEta = fTrueVectorJetEta.at(j) + particle->Eta();
            DeltaPhi = abs(TMath::Pi() - DeltaPhi);
          }
          if(DeltaPhi > TMath::Pi()) {
            DeltaPhi = 2*TMath::Pi() - DeltaPhi;
          }
          Double_t RJetPi0Cand = TMath::Sqrt(DeltaEta*DeltaEta+DeltaPhi*DeltaPhi);
          if(fConvJetReader->Get_Jet_Radius() > 0 ){
            if(fDoOutOfJet == 0){ // in jet
              if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                particleInJet = kTRUE;
                break;
              }
            } else if(fDoOutOfJet == 1){ // out of jet
              if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                particleInJet = kFALSE;
                break;
              }
            } else if(fDoOutOfJet == 2){ // out of jet on away side
              if(RJetPi0Cand < fConvJetReader->Get_Jet_Radius()){
                particleInJet = kTRUE;
                break;
              }
            } else if(fDoOutOfJet == 3){ // out of jet in interval [R, R+0.2]
              if((RJetPi0Cand > fConvJetReader->Get_Jet_Radius()) && (RJetPi0Cand < fConvJetReader->Get_Jet_Radius() + 0.2)){
                particleInJet = kTRUE;
                break;
              }
            }
          }
        }
      }
      Double_t mesonY = 1.e30;
      Double_t ratio  = 0;
      if (particle->E() != TMath::Abs(particle->Pz())){
        ratio         = (particle->E()+particle->Pz()) / (particle->E()-particle->Pz());
      }
      if( !(ratio <= 0) ){
        mesonY = particle->Y()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
      }

      // fill Primary Y hist
      if ( particle->GetPdgCode() == 211 ){  // positve pions
        fHistoMCPrimaryYvsSource[fiCut]->Fill(mesonY, 0., tempParticleWeight);
      } else if ( particle->GetPdgCode() == -211 ){  // negative pions
        fHistoMCPrimaryYvsSource[fiCut]->Fill(mesonY, 1., tempParticleWeight);
      } else if ( particle->GetPdgCode() == 321 ){  // positve kaons
        fHistoMCPrimaryYvsSource[fiCut]->Fill(mesonY, 2., tempParticleWeight);
      } else if ( particle->GetPdgCode() == -321 ){  // negative kaons
        fHistoMCPrimaryYvsSource[fiCut]->Fill(mesonY, 3., tempParticleWeight);
      } else if ( TMath::Abs(particle->GetPdgCode()) == 310 ){  // K0s
        fHistoMCPrimaryYvsSource[fiCut]->Fill(mesonY, 4., tempParticleWeight);
      } else if ( TMath::Abs(particle->GetPdgCode()) == 130 ){  // K0l
        fHistoMCPrimaryYvsSource[fiCut]->Fill(mesonY, 5., tempParticleWeight);
      } else if ( TMath::Abs(particle->GetPdgCode()) == 3122 ){  // Lambda/ AntiLambda
        fHistoMCPrimaryYvsSource[fiCut]->Fill(mesonY, 6., tempParticleWeight);
      }

      if ((mesonY > ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValueMin()) && (mesonY < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValueMax())){
        if ( particle->GetPdgCode() == 211 ){  // positve pions
          fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 0., tempParticleWeight);
        } else if ( particle->GetPdgCode() == -211 ){  // negative pions
          fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 1., tempParticleWeight);
        } else if ( particle->GetPdgCode() == 321 ){  // positve kaons
          fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 2., tempParticleWeight);
        } else if ( particle->GetPdgCode() == -321 ){  // negative kaons
          fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 3., tempParticleWeight);
        } else if ( TMath::Abs(particle->GetPdgCode()) == 310 ){  // K0s
          fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 4., tempParticleWeight);
        } else if ( TMath::Abs(particle->GetPdgCode()) == 130 ){  // K0l
          fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 5., tempParticleWeight);
        } else if ( TMath::Abs(particle->GetPdgCode()) == 3122 ){  // Lambda/ AntiLambda
          fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 6., tempParticleWeight);
        } else if ( particle->GetPdgCode() == 22 ){  // photons
          fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight); // direct photons
          if(particle->GetMother() > -1){
            AliAODMCParticle *mother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetMother()));
            if (  TMath::Abs(mother->GetPdgCode()) == 111  ||
                  TMath::Abs(mother->GetPdgCode()) == 113  ||
                  TMath::Abs(mother->GetPdgCode()) == 221  ||
                  TMath::Abs(mother->GetPdgCode()) == 223  ||
                  TMath::Abs(mother->GetPdgCode()) == 331  ||
                  TMath::Abs(mother->GetPdgCode()) == 333  ||
                  TMath::Abs(mother->GetPdgCode()) == 3212 ||
                  TMath::Abs(mother->GetPdgCode()) == 213
              ){
              fHistoMCDecayGammaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight); // decay photons
            }
          }
        }
      }

      if(fDoJetAnalysis && !particleInJet) continue;

      // check if particle is pi0/eta from di-photon decay
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
        ->MesonIsSelectedAODMC(particle,fAODMCTrackArray,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
        AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(0)));
        AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(1)));

        Float_t weighted= 1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
          if (particle->Pt()>0.005){
            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, 0x0, fInputEvent);
          }
        }
        if( TMath::Abs(particle->GetPdgCode()) == 111 ){ // neutral pions
          fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // All MC Pi0
          if(std::find(fVectorLabelsMultiplePi0Reduced.begin(), fVectorLabelsMultiplePi0Reduced.end(), i) != fVectorLabelsMultiplePi0Reduced.end()) {
            fHistoMCPi0ReducedPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc, ONLY overlapping pi0 in cluster
          }
          if (GetSelectedMesonID() != 2){
            fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
            if (fIsMC==2)fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0 ){
              if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),tempParticleWeight);
            }
          }
        } else if( TMath::Abs(particle->GetPdgCode()) == 221 ){ // eta mesons
          fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // All MC Eta
          if (GetSelectedMesonID() != 1){
            fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
            if (fIsMC==2)fHistoMCEtaWOEvtWeightPt[fiCut]->Fill(particle->Pt());
            if (fDoMesonQA > 0 ){
              if (fIsMC == 2) fHistoMCEtaPtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),tempParticleWeight);
            }
          }
        }
        // Check the acceptance for both gammas & whether they are counted as primaries as well
        Bool_t kDaughter0IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, daughter0, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        Bool_t kDaughter1IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, daughter1, mcProdVtxX, mcProdVtxY, mcProdVtxZ);

        if( kDaughter0IsPrim && kDaughter1IsPrim &&
           (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter0,fAODMCTrackArray) ||
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter1,fAODMCTrackArray) ) ){
          if(particle->GetPdgCode() == 111 && GetSelectedMesonID() != 2){
            fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc
            if(std::find(fVectorLabelsMultiplePi0Reduced.begin(), fVectorLabelsMultiplePi0Reduced.end(), i) != fVectorLabelsMultiplePi0Reduced.end()) {
              fHistoMCPi0ReducedInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc, ONLY overlapping pi0 in cluster
            }
            if (fIsMC == 2)fHistoMCPi0WOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0 w/o event weights with gamma in acc
          } else if(particle->GetPdgCode() == 221 && GetSelectedMesonID() != 1){
            fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Eta with gamma in acc
            if (fIsMC == 2)fHistoMCEtaWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Eta w/o event weights with gamma in acc
          }
        }
      }
      Int_t gammaLabel    = -1;
      Int_t electronLabel = -1;
      Int_t positronLabel = -1;
      // check if particle is pi0/eta from Dalitz decay
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
        ->MesonIsSelectedAODMCDalitz(particle,fAODMCTrackArray, electronLabel, positronLabel, gammaLabel, ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
       if( (gammaLabel > -1) && (electronLabel > -1) && (positronLabel > -1) ){
         AliAODMCParticle* gamma = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaLabel));
         AliAODMCParticle* electron = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(electronLabel));
         AliAODMCParticle* positron = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(positronLabel));

        Float_t weighted= 1;
        if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
          if (particle->Pt()>0.005){
            weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCEvent, fInputEvent);
          }
        }

        if(particle->GetPdgCode() == 111 && GetSelectedMesonID() != 2){
          fHistoMCPi0DalitzPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // All MC Pi0
          fHistoMCPi0DalitzWOWeightPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
          if (fIsMC==2)fHistoMCPi0DalitzWOEvtWeightPt[fiCut]->Fill(particle->Pt());
          if (fDoMesonQA > 0 ){
            if (fIsMC == 2) fHistoMCPi0PtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),tempParticleWeight);
          }
        } else if(particle->GetPdgCode() == 221 && GetSelectedMesonID() != 1){
          fHistoMCEtaDalitzPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // All MC Eta
          fHistoMCEtaDalitzWOWeightPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
          if (fIsMC==2)fHistoMCEtaDalitzWOEvtWeightPt[fiCut]->Fill(particle->Pt());
          if (fDoMesonQA > 0 ){
            if (fIsMC == 2) fHistoMCEtaPtJetPt[fiCut]->Fill(particle->Pt(),((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetMaxPtJet(),tempParticleWeight);
          }
        }

        // Check the acceptance for both gammas & whether they are counted as primaries as well
        Bool_t kGammaIsPrim     = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, gamma, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        Bool_t kElectronIsPrim  = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, electron, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        Bool_t kPositronIsPrim  = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD( fInputEvent, positron, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
        if( kGammaIsPrim && kElectronIsPrim && kPositronIsPrim &&
            (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(gamma,fAODMCTrackArray) ||
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedElecAODMC(electron,fAODMCTrackArray) ||
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedElecAODMC(positron,fAODMCTrackArray) )
          ){
          if(particle->GetPdgCode() == 111 && GetSelectedMesonID() != 2){
            fHistoMCPi0DalitzInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc
            if (fIsMC == 2) fHistoMCPi0DalitzWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc
          } else if(particle->GetPdgCode() == 221 && GetSelectedMesonID() != 1){
            fHistoMCEtaDalitzInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Eta with gamma in acc
            if (fIsMC == 2) fHistoMCEtaDalitzWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc
          }
        }
       }
      }
    // End of primary threatment, now secondary treatment
    } else {

      AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(i));
      if (!particle) continue;

      Int_t isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
        if(isMCFromMBHeader == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempParticleWeight = 1;
      }

      // check if particle is pi0 from di-photon decay
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
        ->MesonIsSelectedAODMC(particle,fAODMCTrackArray,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
        AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(0)));
        AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(1)));
        AliAODMCParticle* mother     = NULL;
        Int_t motherPDG       = -1000000;
        if (particle->GetMother() > -1){
            mother            = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetMother()));
            if (mother)
              motherPDG       = TMath::Abs(mother->GetPdgCode());
        }

        if( TMath::Abs(particle->GetPdgCode()) == 111 ){ // neutral pions
          Int_t source        = GetSourceClassification(111,motherPDG);
          fHistoMCSecPi0PtvsSource[fiCut]->Fill(particle->Pt(),source, tempParticleWeight); // All secondary MC Pi0
        }

        // check whether pi0 landet in acceptance
        if( (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter0,fAODMCTrackArray) ||
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter1,fAODMCTrackArray) ) ){
          if( TMath::Abs(particle->GetPdgCode()) == 111 ){ // neutral pions
            Int_t source      = GetSourceClassification(111,motherPDG);
            fHistoMCSecPi0InAccPtvsSource[fiCut]->Fill(particle->Pt(),source, tempParticleWeight); // All secondary MC Pi0
          }
        }
      }

      Int_t gammaLabel    = -1;
      Int_t electronLabel = -1;
      Int_t positronLabel = -1;
      // check if particle is pi0/eta from Dalitz decay
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
        ->MesonIsSelectedAODMCDalitz(particle,fAODMCTrackArray, electronLabel, positronLabel, gammaLabel, ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift())){
       if( (gammaLabel > -1) && (electronLabel > -1) && (positronLabel > -1) ){
         AliAODMCParticle* gamma = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaLabel));
         AliAODMCParticle* electron = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(electronLabel));
         AliAODMCParticle* positron = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(positronLabel));

        AliAODMCParticle* mother     = NULL;
        Int_t motherPDG       = -1000000;
        if (particle->GetMother() > -1){
            mother            = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetMother()));
            if (mother)
              motherPDG       = TMath::Abs(mother->GetPdgCode());
        }

        if( TMath::Abs(particle->GetPdgCode()) == 111 ){ // neutral pions
          Int_t source        = GetSourceClassification(111,motherPDG);
          fHistoMCSecPi0PtvsSource[fiCut]->Fill(particle->Pt(),source, tempParticleWeight); // All secondary MC Pi0
        }

        // check whether pi0 landet in acceptance
        if( (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(gamma,fAODMCTrackArray) ||
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedElecAODMC(electron,fAODMCTrackArray) ||
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedElecAODMC(positron,fAODMCTrackArray) ) ){
          if( TMath::Abs(particle->GetPdgCode()) == 111 ){ // neutral pions
            Int_t source      = GetSourceClassification(111,motherPDG);
            fHistoMCSecPi0InAccPtvsSource[fiCut]->Fill(particle->Pt(),source, tempParticleWeight); // All secondary MC Pi0
          }
        }
       }
      }
    }
  } // end of particle loop
  return;
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

//_________________________________________________________________________________
Bool_t AliAnalysisTaskGammaCaloMerged::CheckVector1ForEntryAndFillVector2(vector<Int_t> vecCheck, vector<Int_t> &vecFill, Int_t tobechecked)
{
  if(tobechecked > -1)
  {
    vector<Int_t>::iterator itCheck;
    vector<Int_t>::iterator itFill;
    // search for entry in vecCheck
    itCheck = find (vecCheck.begin(), vecCheck.end(), tobechecked);
    if (itCheck != vecCheck.end()){
      // entry in vecCheck -> may not be filled in vecFill -> return false
      return false;
    } else {
      // entry NOT in vecCheck -> entry is allowed to be filled in vecFill
      // check vecFill if entry already exists
      itFill = find (vecFill.begin(), vecFill.end(), tobechecked);
      if (itFill != vecFill.end()){
        // entry already in vecFill -> return false
      } else {
        // entry does not yet exist in vecFill -> fill
        vecFill.push_back(tobechecked);
       return false;
      }
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


//________________________________________________________________________
Int_t AliAnalysisTaskGammaCaloMerged::GetSourceClassification(Int_t daughter, Int_t pdgCode){

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


Bool_t AliAnalysisTaskGammaCaloMerged::NumberOfMCEventNeutralPionOverlapInEMCal(AliMCEvent *mcEvent){
  if (mcEvent && (fMaxAllowedPi0OverlapsMC>-1 || fMinAllowedPi0OverlapsMC > -1)){
    fMapNeutralPionOverlap[fiCut].clear();
    Bool_t failedOverlapCriterium = kFALSE;
    for(Long_t i = 0; i < mcEvent->GetNumberOfPrimaries(); i++) {
      AliMCParticle* particle1 = (AliMCParticle*) mcEvent->GetTrack(i);
      if (!particle1) continue;
      if (TMath::Abs(particle1->PdgCode()) == 111){
        if(TMath::Abs(particle1->Eta()) > 0.67) continue;

        if ( (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) && !(particle1->Phi() > 1.396 && particle1->Phi() < 3.28)) continue;
        if ( (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 3) && !(particle1->Phi() > 4.55 && particle1->Phi() < 5.70)) continue;
        if ( (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 4) && !( (particle1->Phi() > 1.396 && particle1->Phi() < 3.28) || (particle1->Phi() > 4.55 && particle1->Phi() < 5.70) )) continue;

        Int_t nOverlapsFound = 0;
        for(Long_t j = i+1; j < mcEvent->GetNumberOfPrimaries(); j++) {
          AliMCParticle* particle2 = (AliMCParticle*) mcEvent->GetTrack(j);
          if (!particle2) continue;
          if(i==j) continue;
          if((particle1->GetDaughterFirst() == j) || (particle1->GetDaughterLast() == j)) continue;

          if(TMath::Abs(particle2->Eta()) > 0.7) continue;
          if(TMath::Abs(particle2->PdgCode()) == 111 || TMath::Abs(particle2->PdgCode()) == 22){
            Double_t DeltaEta = particle1->Eta()-particle2->Eta();
            Double_t DeltaPhi = abs(particle1->Phi()-particle2->Phi());
            Double_t RneutralParts = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
            if(RneutralParts<0.05){
              nOverlapsFound++;
            }
          }
        }
        if(!fMesonOLFromCluster){
          if( (fMinAllowedPi0OverlapsMC>-1 ? nOverlapsFound<fMinAllowedPi0OverlapsMC : kFALSE) || (fMaxAllowedPi0OverlapsMC>-1 ? nOverlapsFound>fMaxAllowedPi0OverlapsMC : kFALSE)){
            // returns kTRUE as soon as one pi0 with enough overlaps is found
            // if QA is running, then it will loop first over all generated particles and then return failedOverlapCriterium
            if(fDoClusterQA==0) return kTRUE;
            else failedOverlapCriterium = kTRUE;
          }
        } else {
          fMapNeutralPionOverlap[fiCut].insert({particle1->GetLabel(), nOverlapsFound});
        }
        if(fDoClusterQA > 0){
          fHistoOverlapsPi0All[fiCut]->Fill(nOverlapsFound,particle1->Pt(),fWeightJetJetMC);
          if(!failedOverlapCriterium) fHistoOverlapsPi0Accepted[fiCut]->Fill(nOverlapsFound,particle1->Pt(),fWeightJetJetMC);
        }
      }
    }
    return failedOverlapCriterium;
  }
  return kFALSE;
}

void AliAnalysisTaskGammaCaloMerged::ProcessNeutralOverlapsMC(AliMCEvent *mcEvent){
  // cout << "---------------------------------------------------------------------------------------------" << endl;
  if (mcEvent){
    for(Long_t i = 0; i < mcEvent->GetNumberOfPrimaries(); i++) {
      AliMCParticle* particle1 = (AliMCParticle*) mcEvent->GetTrack(i);
      if (!particle1) continue;
      if (TMath::Abs(particle1->PdgCode()) == 111){
        if(particle1->Pt() < 10) continue;
        AliMCParticle* particle1D1 = (AliMCParticle*) mcEvent->GetTrack(particle1->GetDaughterFirst());
        AliMCParticle* particle1D2 = (AliMCParticle*) mcEvent->GetTrack(particle1->GetDaughterLast());
        if(TMath::Abs(particle1->Eta()) > 0.67 || TMath::Abs(particle1D1->Eta()) > 0.67 || TMath::Abs(particle1D2->Eta()) > 0.67) continue;

        // if ( (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) && !(particle1->Phi() > 1.396 && particle1->Phi() < 3.28)) continue;
        // if ( (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 4) && !( (particle1->Phi() > 1.396 && particle1->Phi() < 3.28) || (particle1->Phi() > 4.55 && particle1->Phi() < 5.70) )) continue;

        if ( (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) && !(particle1D1->Phi() > 1.396 && particle1D1->Phi() < 3.28)) continue;
        if ( (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 4) && !( (particle1D1->Phi() > 1.396 && particle1D1->Phi() < 3.28) || (particle1D1->Phi() > 4.55 && particle1D1->Phi() < 5.70) )) continue;
        if ( (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) && !(particle1D2->Phi() > 1.396 && particle1D2->Phi() < 3.28)) continue;
        if ( (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 4) && !( (particle1D2->Phi() > 1.396 && particle1D2->Phi() < 3.28) || (particle1D2->Phi() > 4.55 && particle1D2->Phi() < 5.70) )) continue;
        // both photons must be in a R<0.05 cone
        if(TMath::Sqrt(pow((particle1D1->Eta()-particle1D2->Eta()),2)+pow((abs(particle1D1->Phi()-particle1D2->Phi())),2))>0.05) continue;


        Double_t gammaEnergyInCone = 0;
        Int_t nOverlapsFound = 0;
        // cout << "E_pi0: " << particle1->E() << "\tID:" << i << "\tFirstDaughter:" << particle1->GetDaughterFirst()<< "\tLastDaughter:" << particle1->GetDaughterLast()<< endl;
        // cout << "  pT_g1: " << particle1D1->Pt() << "\tEta: " << particle1D1->Eta()<< "\tPhi: " << particle1D1->Phi() << "\tdR = " << TMath::Sqrt(pow((particle1->Eta()-particle1D1->Eta()),2)+pow((abs(particle1->Phi()-particle1D1->Phi())),2))<< endl;
        // cout << "  pT_g2: " << particle1D2->Pt() << "\tEta: " << particle1D2->Eta()<< "\tPhi: " << particle1D2->Phi() << "\tdR = " << TMath::Sqrt(pow((particle1->Eta()-particle1D2->Eta()),2)+pow((abs(particle1->Phi()-particle1D2->Phi())),2))<< endl;

        Bool_t foundHigherEParticleInCone = kFALSE;
        for(Long_t j = i+1; j < mcEvent->GetNumberOfPrimaries(); j++) {
          AliMCParticle* particle2 = (AliMCParticle*) mcEvent->GetTrack(j);
          if (!particle2) continue;
          if(i==j) continue;
          if(TMath::Abs(particle2->Eta()) > 0.7) continue;
          if(TMath::Abs(particle2->PdgCode()) == 22 && (particle1->GetDaughterFirst() != j) && (particle1->GetDaughterLast() != j)){
            Double_t DeltaEta = particle1->Eta()-particle2->Eta();
            Double_t DeltaPhi = abs(particle1->Phi()-particle2->Phi());
            Double_t RneutralParts = TMath::Sqrt(pow((DeltaEta),2)+pow((DeltaPhi),2));
            if(RneutralParts<0.05){
              if(particle2->Pt()>(particle1D1->Pt()+particle1D2->Pt()) ) foundHigherEParticleInCone = kTRUE;
              nOverlapsFound++;
              gammaEnergyInCone+=particle2->Pt();
              // cout << "\tE_gamma_" << j << ": " << particle2->Pt() << "\tdEta: " << DeltaEta<< "\tdPhi: " << DeltaPhi<< "\tdR = " << RneutralParts << "\tMother: " << particle2->GetMother() << "\tMotherPDG: " << ((AliMCParticle*) mcEvent->GetTrack(particle2->GetMother()))->PdgCode()<< "\tMotherE: " << ((AliMCParticle*) mcEvent->GetTrack(particle2->GetMother()))->E() << endl;
            }
          }
        }
        if(nOverlapsFound>0 && !foundHigherEParticleInCone){
          // cout << "\t\tE_gamma surr: " << gammaEnergyInCone << " from " << nOverlapsFound << " gammas" << endl;
          fHistoPi0EvsGammaOverlapE[fiCut]->Fill(particle1->E(),gammaEnergyInCone);
        }
      }
    }
  }
}


void AliAnalysisTaskGammaCaloMerged::PrintCaloMCLabelsAndInfoAOD(AliVEvent* event, AliAODConversionPhoton *TrueClusterCandidate, AliVCluster* cluster){
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!fAODMCTrackArray) return;
  // cout << endl << endl << "particles contributing: " << endl;
  Double_t energyFromPhotons = 0;
  Double_t energyFromElectrons = 0;
  Double_t energyFromOthers = 0;
  Double_t energyFromLeadingPionParticles = 0;
  for (Int_t i =0 ; i < TrueClusterCandidate->GetNCaloPhotonMCLabels(); i++ ){
    AliAODMCParticle *dummy    = NULL;
    if (TrueClusterCandidate->GetCaloPhotonMCLabel(i)>0){
      dummy        = (AliAODMCParticle*) fAODMCTrackArray->At(TrueClusterCandidate->GetCaloPhotonMCLabel(i));
      // cout << i << "\t ID: "<< TrueClusterCandidate->GetCaloPhotonMCLabel(i) << "\t pdg: " <<dummy->GetPdgCode() << "\t R: "<< TMath::Sqrt( (dummy->Xv()*dummy->Xv()) + (dummy->Yv()*dummy->Yv()) ) << "\t E: " << dummy->E() << "\t Erec: " << cluster->E()* cluster->GetClusterMCEdepFraction(i) << " ( " <<  cluster->GetClusterMCEdepFraction(i)*100 << " %)\t ClsLbl: " << cluster->GetLabelAt(i) ;
      if(TMath::Abs(dummy->GetPdgCode())==22) energyFromPhotons+=cluster->E()* cluster->GetClusterMCEdepFraction(i);
      else if(TMath::Abs(dummy->GetPdgCode())==11) energyFromElectrons+=cluster->E()* cluster->GetClusterMCEdepFraction(i);
      else energyFromOthers+=cluster->E()* cluster->GetClusterMCEdepFraction(i);
      if (dummy->GetMother() > -1){
          AliAODMCParticle* dummyMother   = (AliAODMCParticle*) fAODMCTrackArray->At(dummy->GetMother());
          // cout << "\t MID: " << dummy->GetMother() << "\t Mpdg: " << dummyMother->GetPdgCode() << "\t E: " << dummyMother->E();
          if (dummyMother->GetMother() > -1){
            // AliAODMCParticle* dummyGrandMother   = (AliAODMCParticle*) fAODMCTrackArray->At(dummyMother->GetMother());
            // cout << "\t GMID: " << dummyMother->GetMother() << "\t GMcode: " << dummyGrandMother->GetPdgCode() << "\t E: " << dummyGrandMother->E() << endl;
          } else {
              // cout << endl;
          }
      } else {
        // cout << endl;
      }
    }
  }
  // cout << "leading pion: " << endl;
  if(TrueClusterCandidate->GetNNeutralPionMCLabels()>0){
      // cout << "ID: " << TrueClusterCandidate->GetNeutralPionMCLabel(TrueClusterCandidate->GetLeadingNeutralPionIndex()) << "\t fracE: " << TrueClusterCandidate->GetNeutralPionEnergyFraction(TrueClusterCandidate->GetLeadingNeutralPionIndex()) << endl;
      // cout << "\t leading daughter: ID: " << TrueClusterCandidate->GetCaloPhotonMCLabel(TrueClusterCandidate->GetLeadingNeutralPionDaughterIndex()) << "\t fracE: " << TrueClusterCandidate->GetNeutralPionEnergyFraction(TrueClusterCandidate->GetLeadingNeutralPionDaughterIndex()) << endl;
      energyFromLeadingPionParticles = TrueClusterCandidate->GetNeutralPionEnergyFraction(TrueClusterCandidate->GetLeadingNeutralPionIndex())*cluster->E();
      if(TrueClusterCandidate->GetNeutralPionEnergyFraction(TrueClusterCandidate->GetLeadingNeutralPionIndex())>cluster->GetClusterMCEdepFraction(0)){ // only fill histo if this cluster would be selected for a pi0
      // cout << "frac pi0: " << TrueClusterCandidate->GetNeutralPionEnergyFraction(TrueClusterCandidate->GetLeadingNeutralPionIndex()) << "\t filling " << ((energyFromPhotons+energyFromElectrons)-(energyFromLeadingPionParticles))/(energyFromLeadingPionParticles) << endl;
        if((energyFromPhotons+energyFromElectrons)>0)fHistoTrueClusNeutralContamination[fiCut]->Fill(cluster->E(), ((energyFromPhotons+energyFromElectrons)-(energyFromLeadingPionParticles))/(cluster->E()), fWeightJetJetMC);
        if(energyFromPhotons>0)fHistoTrueClusPhotonContamination[fiCut]->Fill(cluster->E(), ((energyFromPhotons)-(energyFromLeadingPionParticles))/(cluster->E()), fWeightJetJetMC);
        if(energyFromElectrons>0)fHistoTrueClusElectronContamination[fiCut]->Fill(cluster->E(), ((energyFromElectrons)-(energyFromLeadingPionParticles))/(cluster->E()), fWeightJetJetMC);
        if(energyFromOthers>0)fHistoTrueClusOtherContamination[fiCut]->Fill(cluster->E(), ((energyFromOthers)-(energyFromLeadingPionParticles))/(cluster->E()), fWeightJetJetMC);
        // if((energyFromPhotons+energyFromElectrons)>0)fHistoTrueClusNeutralContamination[fiCut]->Fill(cluster->E(), ((energyFromPhotons+energyFromElectrons)-(energyFromLeadingPionParticles))/(energyFromLeadingPionParticles), fWeightJetJetMC);
        // if(energyFromPhotons>0)fHistoTrueClusPhotonContamination[fiCut]->Fill(cluster->E(), ((energyFromPhotons)-(energyFromLeadingPionParticles))/(energyFromLeadingPionParticles), fWeightJetJetMC);
        // if(energyFromElectrons>0)fHistoTrueClusElectronContamination[fiCut]->Fill(cluster->E(), ((energyFromElectrons)-(energyFromLeadingPionParticles))/(energyFromLeadingPionParticles), fWeightJetJetMC);
        // if(energyFromOthers>0)fHistoTrueClusOtherContamination[fiCut]->Fill(cluster->E(), ((energyFromOthers)-(energyFromLeadingPionParticles))/(energyFromLeadingPionParticles), fWeightJetJetMC);
      }
  }

}
