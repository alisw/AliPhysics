/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                       *
 * Author: Baldo Sahlmueller, Friederike Bock                     *
 * Version 1.0                                 *
 *                                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its    *
 * documentation strictly for non-commercial purposes is hereby granted    *
 * without fee, provided that the above copyright notice appears in all    *
 * copies and that both the copyright notice and this permission notice    *
 * appear in the supporting documentation. The authors make no claims    *
 * about the suitability of this software for any purpose. It is      *
 * provided "as is" without express or implied warranty.               *
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
#include "AliAnalysisTaskPi0v2Calo.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliGAKFVertex.h"
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
#include "AliCaloTriggerMimicHelper.h"
#include <vector>
#include <map>
#include <fstream>

ClassImp(AliAnalysisTaskPi0v2Calo)

//________________________________________________________________________
AliAnalysisTaskPi0v2Calo::AliAnalysisTaskPi0v2Calo(): AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fCaloTriggerHelperName(""),
  fCorrTaskSetting(""),
  fPeriod(""),
  fBGHandler(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fQAList(NULL),
  fBackList(NULL),
  fMotherList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fOutputContainer(NULL),
  fReaderGammas(NULL),
  fGammaCandidates(NULL),
  fClusterCandidates(NULL),
  fEventCutArray(NULL),
  fEventCuts(NULL),
  fClusterCutArray(NULL),
  fCaloPhotonCuts(NULL),
  fMesonCutArray(NULL),
  fMesonCuts(NULL),
  fOutlierJetReader(NULL),
  fConversionCuts(NULL),
  fCaloTriggerMimicHelper(NULL),
  fSetEventCutsOutputlist(),
  fHistoMotherInvMassPtPhiV0A(NULL),
  fHistoMotherBackInvMassPtPhiV0A(NULL),
  fHistoMotherInvMassPtPhiV0C(NULL),
  fHistoMotherBackInvMassPtPhiV0C(NULL),
  fHistoMotherInvMassPtV0CInPlane(NULL),
  fHistoMotherInvMassPtV0AInPlane(NULL),
  fHistoMotherInvMassPtV0COutPlane(NULL),
  fHistoMotherInvMassPtV0AOutPlane(NULL),
  fHistoMotherBackInvMassPtV0CInPlane(NULL),
  fHistoMotherBackInvMassPtV0AInPlane(NULL),
  fHistoMotherBackInvMassPtV0COutPlane(NULL),
  fHistoMotherBackInvMassPtV0AOutPlane(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusGammaE(NULL),
  fHistoClusOverlapHeadersGammaPt(NULL),
  fHistoClusOverlapMBHeaderGammaPt(NULL),
  fHistoClusAllHeadersGammaPt(NULL),
  fHistoClusRejectedHeadersGammaPt(NULL),
  fHistoClusGammaPtM02(NULL),
  fHistoMCHeaders(NULL),
  fHistoMCEventsTrigg(NULL),
  fHistoMCGammaPtNotTriggered(NULL),
  fHistoMCGammaPtNoVertex(NULL),
  fHistoMCAllGammaPt(NULL),
  fHistoMCAllSecondaryGammaPt(NULL),
  fHistoMCDecayGammaPi0Pt(NULL),
  fHistoMCDecayGammaRhoPt(NULL),
  fHistoMCDecayGammaEtaPt(NULL),
  fHistoMCDecayGammaOmegaPt(NULL),
  fHistoMCDecayGammaEtapPt(NULL),
  fHistoMCDecayGammaPhiPt(NULL),
  fHistoMCDecayGammaSigmaPt(NULL),
  fHistoMCPi0Pt(NULL),
  fHistoMCPi0PtNotTriggered(NULL),
  fHistoMCPi0PtNoVertex(NULL),
  fHistoMCPi0WOWeightPt(NULL),
  fHistoMCPi0WOEvtWeightPt(NULL),
  fHistoMCEtaPt(NULL),
  fHistoMCEtaPtNotTriggered(NULL),
  fHistoMCEtaPtNoVertex(NULL),
  fHistoMCEtaWOWeightPt(NULL),
  fHistoMCEtaWOEvtWeightPt(NULL),
  fHistoMCPi0InAccPt(NULL),
  fHistoMCEtaInAccPt(NULL),
  fHistoMCPi0InAccPtNotTriggered(NULL),
  fHistoMCEtaInAccPtNotTriggered(NULL),
  fHistoMCPi0WOEvtWeightInAccPt(NULL),
  fHistoMCEtaWOEvtWeightInAccPt(NULL),
  fHistoMCPrimaryPtvsSource(NULL),
  fHistoMCSecPi0PtvsSource(NULL),
  fHistoMCSecPi0Source(NULL),
  fHistoMCSecPi0InAccPtvsSource(NULL),
  fHistoMCSecEtaPt(NULL),
  fHistoMCSecEtaSource(NULL),
  fHistoMCPi0PtJetPt(NULL),
  fHistoMCEtaPtJetPt(NULL),
  fHistoTruePi0InvMassPt(NULL),
  fHistoTruePi0InvMassPtAdditional(NULL),
  fHistoTrueEtaInvMassPt(NULL),
  fHistoTrueEtaInvMassPtAdditional(NULL),
  fHistoTruePrimaryPi0InvMassPt(NULL),
  fHistoTruePrimaryEtaInvMassPt(NULL),
  fHistoTruePrimaryPi0W0WeightingInvMassPt(NULL),
  fHistoTruePrimaryEtaW0WeightingInvMassPt(NULL),
  fProfileTruePrimaryPi0WeightsInvMassPt(NULL),
  fProfileTruePrimaryEtaWeightsInvMassPt(NULL),
  fHistoTrueSecondaryPi0InvMassPt(NULL),
  fHistoTrueSecondaryPi0FromK0sInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromK0lInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromEtaInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromLambdaInvMassPt(NULL),
  fHistoClusPhotonBGPt(NULL),
  fHistoClusPhotonPlusConvBGPt(NULL),
  fHistoTrueClusGammaPt(NULL),
  fHistoTrueClusGammaEResE(NULL),
  fHistoTrueClusPhotonGammaEResE(NULL),
  fHistoTrueNLabelsInClus(NULL),
  fHistoTruePrimaryClusGammaPt(NULL),
  fHistoTruePrimaryClusGammaESDPtMCPt(NULL),
  fHistoTruePrimaryClusConvGammaPt(NULL),
  fHistoTruePrimaryClusConvGammaESDPtMCPt(NULL),
  fHistoTrueSecondaryClusGammaPt(NULL),
  fHistoTrueSecondaryClusConvGammaPt(NULL),
  fHistoTrueSecondaryClusGammaMCPt(NULL),
  fHistoTrueSecondaryClusConvGammaMCPt(NULL),
  fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt(NULL),
  fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt(NULL),
  fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt(NULL),
  fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt(NULL),
  fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt(NULL),
  fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt(NULL),
  fHistoDoubleCountTruePi0InvMassPt(NULL),
  fHistoDoubleCountTrueEtaInvMassPt(NULL),
  fHistoDoubleCountTrueClusterGammaPt(NULL),
  fVectorDoubleCountTruePi0s(0),
  fVectorDoubleCountTrueEtas(0),
  fVectorDoubleCountTrueClusterGammas(0),
  fHistoMultipleCountTrueClusterGamma(NULL),
  fMapMultipleCountTrueClusterGammas(),
  fHistoTruePi0InvMassPtAlpha(NULL),
  fHistoTruePi0PureGammaInvMassPtAlpha(NULL),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoVertexZ(NULL),
  fHistoNGammaCandidates(NULL),
  fHistoNGammaCandidatesBasic(NULL),
  fHistoNGoodESDTracksVsNGammaCandidates(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoNV0Tracks(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoJetJetNTrials(NULL),
  fHistoPtHardJJWeight(NULL),
  fHistoPionSpectrum(NULL),
  fHistoProtonSpectrum(NULL),
  fHistoKaonSpectrum(NULL),
  fHistoNPionSpectrum(NULL),
  fHistoEtaSpectrum(NULL),
  fHistoDMesonSpectrum(NULL),
  fHistoMotherInvMassPtV0CCos2phi(NULL),
  fHistoMotherInvMassPtV0ACos2phi(NULL),
  fHistoMotherBackInvMassPtV0CCos2phi(NULL),
  fHistoMotherBackInvMassPtV0ACos2phi(NULL),
  fRunNumber(0),
  fOldRunNumber(0),
  IsVZEROCalibOn(kTRUE),
  IsQAVZERO(kTRUE),
  fUseInOutPlane(kFALSE),
  fListVZEROCalib(NULL),
  fVZEROCalibFile(NULL),
  fPsi2V0C(-999),
  fPsi2V0A(-999),
  fHist2DPsi2V0CCent(NULL),
  fHist2DPsi2V0ACent(NULL),
  hMultV0(NULL),
  contMult(NULL),
  contQxncm(NULL),
  contQyncm(NULL),
  contQxnam(NULL),
  contQynam(NULL),
  fHCorrectV0ChWeghts(NULL),
  fProfileV0CQxCentGE(NULL),
  fProfileV0CQyCentGE(NULL),
  fProfileV0CQxVtxGE(NULL),
  fProfileV0CQyVtxGE(NULL),
  fHist2CalibPsi2V0CCentGE(NULL),
  fProfileV0AQxCentGE(NULL),
  fProfileV0AQyCentGE(NULL),
  fProfileV0AQxVtxGE(NULL),
  fProfileV0AQyVtxGE(NULL),
  fHist2CalibPsi2V0ACentGE(NULL),
  fProfileV0CQxCentRC(NULL),
  fProfileV0CQyCentRC(NULL),
  fProfileV0CQxVtxRC(NULL),
  fProfileV0CQyVtxRC(NULL),
  fHist2CalibPsi2V0CCentRC(NULL),
  fProfileV0AQxCentRC(NULL),
  fProfileV0AQyCentRC(NULL),
  fProfileV0AQxVtxRC(NULL),
  fProfileV0AQyVtxRC(NULL),
  fHist2CalibPsi2V0ACentRC(NULL),
  fHist2V0Res(NULL),
  fCent(-1),
  fEventPlaneAngle(-100),
  fRandom(0),
  fnCuts(0),
  fiCut(0),
  fIsHeavyIon(0),
  fDoLightOutput(0),
  fDoMesonAnalysis(kTRUE),
  fIsFromDesiredHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fIsOverlapWithMBHeader(kFALSE),
  fIsMC(0),
  fWeightJetJetMC(1),
  fDoInOutTimingCluster(kFALSE),
  fMinTimingCluster(0),
  fMaxTimingCluster(0),
  fEnableSortForClusMC(kFALSE),
  fDoPrimaryTrackMatching(kTRUE),
  tBrokenFiles(NULL),
  fFileNameBroken(NULL),
  fFileNameTrigger(NULL),
  fGenPhaseSpace(),
  fAODMCTrackArray(NULL),
  fLocalDebugFlag(0),
  fAllowOverlapHeaders(kTRUE),
  fNCurrentClusterBasic(0),
  fTrackMatcherRunningMode(0),
  fDoPi0Only(kFALSE)
{
  for (int i = 0; i < 2; ++i)
    {
      hQx2mV0[i] = NULL;
      hQy2mV0[i] = NULL;
    }
}

//________________________________________________________________________
AliAnalysisTaskPi0v2Calo::AliAnalysisTaskPi0v2Calo(const char *name):
  AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fCaloTriggerHelperName(""),
  fCorrTaskSetting(""),
  fPeriod(""),
  fBGHandler(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fQAList(NULL),
  fBackList(NULL),
  fMotherList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fOutputContainer(0),
  fReaderGammas(NULL),
  fGammaCandidates(NULL),
  fClusterCandidates(NULL),
  fEventCutArray(NULL),
  fEventCuts(NULL),
  fClusterCutArray(NULL),
  fCaloPhotonCuts(NULL),
  fMesonCutArray(NULL),
  fMesonCuts(NULL),
  fOutlierJetReader(NULL),
  fConversionCuts(NULL),
  fCaloTriggerMimicHelper(NULL),
  fSetEventCutsOutputlist(),
  fHistoMotherInvMassPtPhiV0A(NULL),
  fHistoMotherBackInvMassPtPhiV0A(NULL),
  fHistoMotherInvMassPtPhiV0C(NULL),
  fHistoMotherBackInvMassPtPhiV0C(NULL),
  fHistoMotherInvMassPtV0CInPlane(NULL),
  fHistoMotherInvMassPtV0AInPlane(NULL),
  fHistoMotherInvMassPtV0COutPlane(NULL),
  fHistoMotherInvMassPtV0AOutPlane(NULL),
  fHistoMotherBackInvMassPtV0CInPlane(NULL),
  fHistoMotherBackInvMassPtV0AInPlane(NULL),
  fHistoMotherBackInvMassPtV0COutPlane(NULL),
  fHistoMotherBackInvMassPtV0AOutPlane(NULL),
  fHistoClusGammaPt(NULL),
  fHistoClusGammaE(NULL),
  fHistoClusOverlapHeadersGammaPt(NULL),
  fHistoClusOverlapMBHeaderGammaPt(NULL),
  fHistoClusAllHeadersGammaPt(NULL),
  fHistoClusRejectedHeadersGammaPt(NULL),
  fHistoClusGammaPtM02(NULL),
  fHistoMCHeaders(NULL),
  fHistoMCEventsTrigg(NULL),
  fHistoMCGammaPtNotTriggered(NULL),
  fHistoMCGammaPtNoVertex(NULL),
  fHistoMCAllGammaPt(NULL),
  fHistoMCAllSecondaryGammaPt(NULL),
  fHistoMCDecayGammaPi0Pt(NULL),
  fHistoMCDecayGammaRhoPt(NULL),
  fHistoMCDecayGammaEtaPt(NULL),
  fHistoMCDecayGammaOmegaPt(NULL),
  fHistoMCDecayGammaEtapPt(NULL),
  fHistoMCDecayGammaPhiPt(NULL),
  fHistoMCDecayGammaSigmaPt(NULL),
  fHistoMCPi0Pt(NULL),
  fHistoMCPi0PtNotTriggered(NULL),
  fHistoMCPi0PtNoVertex(NULL),
  fHistoMCPi0WOWeightPt(NULL),
  fHistoMCPi0WOEvtWeightPt(NULL),
  fHistoMCEtaPt(NULL),
  fHistoMCEtaPtNotTriggered(NULL),
  fHistoMCEtaPtNoVertex(NULL),
  fHistoMCEtaWOWeightPt(NULL),
  fHistoMCEtaWOEvtWeightPt(NULL),
  fHistoMCPi0InAccPt(NULL),
  fHistoMCEtaInAccPt(NULL),
  fHistoMCPi0InAccPtNotTriggered(NULL),
  fHistoMCEtaInAccPtNotTriggered(NULL),
  fHistoMCPi0WOEvtWeightInAccPt(NULL),
  fHistoMCEtaWOEvtWeightInAccPt(NULL),
  fHistoMCPrimaryPtvsSource(NULL),
  fHistoMCSecPi0PtvsSource(NULL),
  fHistoMCSecPi0Source(NULL),
  fHistoMCSecPi0InAccPtvsSource(NULL),
  fHistoMCSecEtaPt(NULL),
  fHistoMCSecEtaSource(NULL),
  fHistoMCPi0PtJetPt(NULL),
  fHistoMCEtaPtJetPt(NULL),
  fHistoTruePi0InvMassPt(NULL),
  fHistoTruePi0InvMassPtAdditional(NULL),
  fHistoTrueEtaInvMassPt(NULL),
  fHistoTrueEtaInvMassPtAdditional(NULL),
  fHistoTruePrimaryPi0InvMassPt(NULL),
  fHistoTruePrimaryEtaInvMassPt(NULL),
  fHistoTruePrimaryPi0W0WeightingInvMassPt(NULL),
  fHistoTruePrimaryEtaW0WeightingInvMassPt(NULL),
  fProfileTruePrimaryPi0WeightsInvMassPt(NULL),
  fProfileTruePrimaryEtaWeightsInvMassPt(NULL),
  fHistoTrueSecondaryPi0InvMassPt(NULL),
  fHistoTrueSecondaryPi0FromK0sInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromK0lInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromEtaInvMassPt(NULL),
  fHistoTrueSecondaryPi0FromLambdaInvMassPt(NULL),
  fHistoClusPhotonBGPt(NULL),
  fHistoClusPhotonPlusConvBGPt(NULL),
  fHistoTrueClusGammaPt(NULL),
  fHistoTrueClusGammaEResE(NULL),
  fHistoTrueClusPhotonGammaEResE(NULL),
  fHistoTrueNLabelsInClus(NULL),
  fHistoTruePrimaryClusGammaPt(NULL),
  fHistoTruePrimaryClusGammaESDPtMCPt(NULL),
  fHistoTruePrimaryClusConvGammaPt(NULL),
  fHistoTruePrimaryClusConvGammaESDPtMCPt(NULL),
  fHistoTrueSecondaryClusGammaPt(NULL),
  fHistoTrueSecondaryClusConvGammaPt(NULL),
  fHistoTrueSecondaryClusGammaMCPt(NULL),
  fHistoTrueSecondaryClusConvGammaMCPt(NULL),
  fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt(NULL),
  fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt(NULL),
  fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt(NULL),
  fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt(NULL),
  fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt(NULL),
  fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt(NULL),
  fHistoDoubleCountTruePi0InvMassPt(NULL),
  fHistoDoubleCountTrueEtaInvMassPt(NULL),
  fHistoDoubleCountTrueClusterGammaPt(NULL),
  fVectorDoubleCountTruePi0s(0),
  fVectorDoubleCountTrueEtas(0),
  fVectorDoubleCountTrueClusterGammas(0),
  fHistoMultipleCountTrueClusterGamma(NULL),
  fMapMultipleCountTrueClusterGammas(),
  fHistoTruePi0InvMassPtAlpha(NULL),
  fHistoTruePi0PureGammaInvMassPtAlpha(NULL),
  fHistoNEvents(NULL),
  fHistoNEventsWOWeight(NULL),
  fHistoNGoodESDTracks(NULL),
  fHistoVertexZ(NULL),
  fHistoNGammaCandidates(NULL),
  fHistoNGammaCandidatesBasic(NULL),
  fHistoNGoodESDTracksVsNGammaCandidates(NULL),
  fHistoSPDClusterTrackletBackground(NULL),
  fHistoNV0Tracks(NULL),
  fProfileEtaShift(NULL),
  fProfileJetJetXSection(NULL),
  fHistoJetJetNTrials(NULL),
  fHistoPtHardJJWeight(NULL),
  fHistoPionSpectrum(NULL),
  fHistoProtonSpectrum(NULL),
  fHistoKaonSpectrum(NULL),
  fHistoNPionSpectrum(NULL),
  fHistoEtaSpectrum(NULL),
  fHistoDMesonSpectrum(NULL),
  fHistoMotherInvMassPtV0CCos2phi(NULL),
  fHistoMotherInvMassPtV0ACos2phi(NULL),
  fHistoMotherBackInvMassPtV0CCos2phi(NULL),
  fHistoMotherBackInvMassPtV0ACos2phi(NULL),
  fRunNumber(0),
  fOldRunNumber(0),
  IsVZEROCalibOn(kTRUE),
  IsQAVZERO(kTRUE),
  fUseInOutPlane(kFALSE),
  fListVZEROCalib(NULL),
  fVZEROCalibFile(NULL),
  fPsi2V0C(-999),
  fPsi2V0A(-999),
  fHist2DPsi2V0CCent(NULL),
  fHist2DPsi2V0ACent(NULL),
  hMultV0(NULL),
  contMult(NULL),
  contQxncm(NULL),
  contQyncm(NULL),
  contQxnam(NULL),
  contQynam(NULL),
  fHCorrectV0ChWeghts(NULL),
  fProfileV0CQxCentGE(NULL),
  fProfileV0CQyCentGE(NULL),
  fProfileV0CQxVtxGE(NULL),
  fProfileV0CQyVtxGE(NULL),
  fHist2CalibPsi2V0CCentGE(NULL),
  fProfileV0AQxCentGE(NULL),
  fProfileV0AQyCentGE(NULL),
  fProfileV0AQxVtxGE(NULL),
  fProfileV0AQyVtxGE(NULL),
  fHist2CalibPsi2V0ACentGE(NULL),
  fProfileV0CQxCentRC(NULL),
  fProfileV0CQyCentRC(NULL),
  fProfileV0CQxVtxRC(NULL),
  fProfileV0CQyVtxRC(NULL),
  fHist2CalibPsi2V0CCentRC(NULL),
  fProfileV0AQxCentRC(NULL),
  fProfileV0AQyCentRC(NULL),
  fProfileV0AQxVtxRC(NULL),
  fProfileV0AQyVtxRC(NULL),
  fHist2CalibPsi2V0ACentRC(NULL),
  fHist2V0Res(NULL),
  fCent(-1),
  fEventPlaneAngle(-100),
  fRandom(0),
  fnCuts(0),
  fiCut(0),
  fIsHeavyIon(0),
  fDoLightOutput(0),
  fDoMesonAnalysis(kTRUE),
  fIsFromDesiredHeader(kTRUE),
  fIsOverlappingWithOtherHeader(kFALSE),
  fIsOverlapWithMBHeader(kFALSE),
  fIsMC(0),
  fWeightJetJetMC(1),
  fDoInOutTimingCluster(kFALSE),
  fMinTimingCluster(0),
  fMaxTimingCluster(0),
  fEnableSortForClusMC(kFALSE),
  fDoPrimaryTrackMatching(kTRUE),
  tBrokenFiles(NULL),
  fFileNameBroken(NULL),
  fFileNameTrigger(NULL),
  fGenPhaseSpace(),
  fAODMCTrackArray(NULL),
  fLocalDebugFlag(0),
  fAllowOverlapHeaders(kTRUE),
  fNCurrentClusterBasic(0),
  fTrackMatcherRunningMode(0),
  fDoPi0Only(kFALSE)
{
  for (int i = 0; i < 2; ++i)
  {
    hQx2mV0[i] = NULL;
    hQy2mV0[i] = NULL;
  }
  // Define output slots here
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
  DefineOutput(5, TTree::Class());
  DefineOutput(6, TTree::Class());
  DefineOutput(7, TTree::Class());
  DefineOutput(8, TTree::Class());
  DefineOutput(9, TTree::Class());}

AliAnalysisTaskPi0v2Calo::~AliAnalysisTaskPi0v2Calo()
{
  if(fClusterCandidates){
    delete fClusterCandidates;
    fClusterCandidates = 0x0;
  }
  if(fGammaCandidates){
    delete fGammaCandidates;
    fGammaCandidates = 0x0;
  }
  if(fBGHandler){
    delete[] fBGHandler;
    fBGHandler = 0x0;
  }
}
//___________________________________________________________
void AliAnalysisTaskPi0v2Calo::InitBack(){

  fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];


  for(int iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
      TString cutstringEvent  = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringCalo   = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson  = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

      int collisionSystem   = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(0,1));
      int centMin           = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(1,1));
      int centMax           = atoi((TString)(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber())(2,1));

      if(collisionSystem == 1 || collisionSystem == 2 ||
        collisionSystem == 5 || collisionSystem == 8 ||
        collisionSystem == 9){
        centMin   = centMin*10;
        centMax   = centMax*10;
        if(centMax ==0 && centMax!=centMin) centMax=100;
      } else if(collisionSystem == 3 || collisionSystem == 6){
        centMin   = centMin*5;
        centMax   = centMax*5;
      } else if(collisionSystem == 4 || collisionSystem == 7){
        centMin   = ((centMin*5)+45);
        centMax   = ((centMax*5)+45);
      }

      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
        if( ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoSectorMixing() ){
          fBGHandler[iCut] = new AliGammaConversionAODBGHandler(
                                    collisionSystem,centMin,centMax,
                                    ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
                                    ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity(),
                                    4,9,7);
        } else {
          fBGHandler[iCut] = new AliGammaConversionAODBGHandler(
                                    collisionSystem,centMin,centMax,
                                    ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
                                    ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity(),
                                    4,8,7);
        }
      }
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskPi0v2Calo::UserCreateOutputObjects(){

  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader
  if(((AliConvEventCuts*)fEventCutArray->At(0))->GetUseJetFinderForOutliers()){
    fOutlierJetReader=(AliAnalysisTaskJetOutlierRemoval*)AliAnalysisManager::GetAnalysisManager()->GetTask("AliAnalysisTaskJetOutlierRemoval");
    if(!fOutlierJetReader){AliFatal("Error: No AliAnalysisTaskJetOutlierRemoval");} // GetV0Reader
    else{printf("Found AliAnalysisTaskJetOutlierRemoval used for outlier removal!\n");}
  }
  if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoSecondaryTrackMatching()){
    fConversionCuts = new AliConversionPhotonCuts();
    fConversionCuts->SetV0ReaderName(fV0ReaderName.Data());
    fConversionCuts->InitializeCutsFromCutString("00200009327000008250400000"); //Use standard cuts
    fConversionCuts->SetIsHeavyIon(fIsHeavyIon);
    fConversionCuts->SetFillCutHistograms("",kFALSE);
  }

  // set common binning in pT for mesons and photons
  double epsilon              = 1.e-6;
  float binWidthPt          = 0.1;
  int nBinsPt               = 250;
  float minPt               = 0;
  float maxPt               = 25;
  int nBinsClusterPt        = 500;
  float minClusterPt        = 0;
  float maxClusterPt        = 50;
  int nBinsMinv             = 800;
  float maxMinv             = 0.8;
  float minRes                = -5.f;
  float maxRes                = +1.f; 
  double *arrPtBinning      = new double[nBinsPt+1];
  double *arrClusPtBinning  = new double[nBinsClusterPt+1];
  std::vector<double> arrResBinning(201, 0.);
  if( fDoPi0Only ){
    nBinsMinv                 = 150;
    maxMinv                  = 0.3;
  }

  if (((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEnergyEnum() == AliConvEventCuts::kPbPb5TeV  ){
    nBinsMinv                 = 400;
    nBinsPt                   = 49;
    minPt                     = 0;
    maxPt                     = 40;
    for(int i=0; i<nBinsPt+1;i++){
      if(i<=20)           arrPtBinning[i]  = 0.0   + 0.2*i;             // 0.2 GeV bin width until 4 GeV
      else if(i<=30)      arrPtBinning[i]  = 4.0   + 0.5*(i-20);        // 0.5 GeV                 9 GeV
      else if(i<=41)      arrPtBinning[i]  = 9.0   + 1.0*(i-30);        // 1.0 GeV                 15 GeV
      else if(i<=49) arrPtBinning[i]       = 20.   + 2.5*(i-41);        // 2.5 GeV                 40 GeV
      else arrPtBinning[i]                    = maxPt;
    }
    nBinsClusterPt            = 148;
    minClusterPt              = 0;
    maxClusterPt              = 40;
    for(int i=0; i<nBinsClusterPt+1;i++){
      if (i < 1) arrClusPtBinning[i]          = 0.3*i;
      else if(i<98) arrClusPtBinning[i]       = 0.3+0.1*(i-1);
      else if(i<123) arrClusPtBinning[i]      = 10.+0.2*(i-98);
      else if(i<148) arrClusPtBinning[i]      = 15.+1.0*(i-123);
      else arrClusPtBinning[i]                = maxClusterPt;
    }
  } else {
    for(int i=0; i<nBinsPt+1;i++){
      arrPtBinning[i]         = ((maxPt-minPt)/nBinsPt)*i;
    }
    for(int i=0; i<nBinsClusterPt+1;i++){
      arrClusPtBinning[i]     = ((maxClusterPt-minClusterPt)/nBinsClusterPt)*i;
    }
  }
  double minPhi = 0.;
  double maxPhi = TMath::Pi();
  int nPhiBins = 6;
  double arrPhiBin[nPhiBins + 1];
  for (int i = 0; i < nPhiBins + 1; i++) {
    arrPhiBin[i] = minPhi + ((maxPhi - minPhi) / nPhiBins) * i;
  }
  double arrMinvBin[nBinsMinv + 1];
  for (int i = 0; i < nBinsMinv + 1; i++) {
    arrMinvBin[i] = (maxMinv / nBinsMinv) * i;
  }

  double valRes = minRes;
  for (size_t i = 0; i < arrResBinning.size(); ++i) {
    arrResBinning.at(i) = valRes;
    if (valRes < -0.5 - epsilon)
      valRes += 0.05;
    else if (valRes < 0.5 - epsilon)
      valRes += 0.01;
    else if (valRes < maxRes - epsilon)
      valRes += 0.05;
    else
      break;
  }

  // Create histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer  = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer  = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

  // Array of current cut's gammas
  if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoSecondaryTrackMatching()){
    fGammaCandidates    = new TList();
  }
  fClusterCandidates  = new TList();
  fClusterCandidates->SetOwner(kTRUE);
  fCutFolder          = new TList*[fnCuts];
  fESDList            = new TList*[fnCuts];
  fQAList             = new TList*[fnCuts];
  fHistoNEvents       = new TH1F*[fnCuts];
  if(fIsMC > 1){
    fHistoNEventsWOWeight   = new TH1F*[fnCuts];
  }
  if(fIsMC == 2){
    fProfileJetJetXSection  = new TProfile*[fnCuts];
    fHistoJetJetNTrials     = new TH1F*[fnCuts];
    fHistoPtHardJJWeight     = new TH2F*[fnCuts];
  }

  fHistoNGoodESDTracks        = new TH1F*[fnCuts];
  fHistoVertexZ               = new TH1F*[fnCuts];
  fHistoNGammaCandidates      = new TH1F*[fnCuts];
  fHistoNGammaCandidatesBasic = new TH1F*[fnCuts];
  if(!fDoLightOutput){
    fHistoNGoodESDTracksVsNGammaCandidates  = new TH2F*[fnCuts];
    fHistoSPDClusterTrackletBackground      = new TH2F*[fnCuts];
    fHistoNV0Tracks                         = new TH1F*[fnCuts];
  }
  if(fIsHeavyIon==2) fProfileEtaShift          = new TProfile*[fnCuts];
  if(fDoMesonAnalysis){
    if(!fUseInOutPlane){
      fHistoMotherInvMassPtPhiV0A     = new TH3F*[fnCuts];
      fHistoMotherBackInvMassPtPhiV0A = new TH3F*[fnCuts];
      fHistoMotherInvMassPtPhiV0C     = new TH3F*[fnCuts];
      fHistoMotherBackInvMassPtPhiV0C = new TH3F*[fnCuts];
    } else{
      fHistoMotherInvMassPtV0CInPlane = new TH2F *[fnCuts];
      fHistoMotherInvMassPtV0AInPlane = new TH2F *[fnCuts];
      fHistoMotherInvMassPtV0COutPlane = new TH2F *[fnCuts];
      fHistoMotherInvMassPtV0AOutPlane = new TH2F *[fnCuts];
      fHistoMotherBackInvMassPtV0CInPlane = new TH2F *[fnCuts];
      fHistoMotherBackInvMassPtV0AInPlane = new TH2F *[fnCuts];
      fHistoMotherBackInvMassPtV0COutPlane = new TH2F *[fnCuts];
      fHistoMotherBackInvMassPtV0AOutPlane = new TH2F *[fnCuts];
    }
    fHistoMotherInvMassPtV0CCos2phi = new TProfile2D *[fnCuts];
    fHistoMotherInvMassPtV0ACos2phi = new TProfile2D *[fnCuts];

    fHistoMotherInvMassPtV0CCos2phi = new TProfile2D *[fnCuts];
    fHistoMotherInvMassPtV0ACos2phi = new TProfile2D *[fnCuts];
    fHistoMotherBackInvMassPtV0CCos2phi = new TProfile2D *[fnCuts];
    fHistoMotherBackInvMassPtV0ACos2phi = new TProfile2D *[fnCuts];
    fHist2DPsi2V0CCent = new TH2D *[fnCuts];
    fHist2DPsi2V0ACent = new TH2D *[fnCuts];
    fProfileV0CQxCentGE = new TProfile *[fnCuts];
    fProfileV0CQyCentGE = new TProfile *[fnCuts];
    fProfileV0CQxVtxGE = new TProfile *[fnCuts];
    fProfileV0CQyVtxGE = new TProfile *[fnCuts];
    fHist2CalibPsi2V0CCentGE = new TH2D *[fnCuts];
    fProfileV0AQxCentGE = new TProfile *[fnCuts];
    fProfileV0AQyCentGE = new TProfile *[fnCuts];
    fProfileV0AQxVtxGE = new TProfile *[fnCuts];
    fProfileV0AQyVtxGE = new TProfile *[fnCuts];
    fHist2CalibPsi2V0ACentGE = new TH2D *[fnCuts];
    fProfileV0CQxCentRC = new TProfile *[fnCuts];
    fProfileV0CQyCentRC = new TProfile *[fnCuts];
    fProfileV0CQxVtxRC = new TProfile *[fnCuts];
    fProfileV0CQyVtxRC = new TProfile *[fnCuts];
    fHist2CalibPsi2V0CCentRC = new TH2D *[fnCuts];
    fProfileV0AQxCentRC = new TProfile *[fnCuts];
    fProfileV0AQyCentRC = new TProfile *[fnCuts];
    fProfileV0AQxVtxRC = new TProfile *[fnCuts];
    fProfileV0AQyVtxRC = new TProfile *[fnCuts];
    fHist2CalibPsi2V0ACentRC = new TH2D *[fnCuts];
    fHist2V0Res = new TProfile *[fnCuts];
  }
  fHistoClusGammaPt                 = new TH1F*[fnCuts];
  fHistoClusGammaE                  = new TH1F*[fnCuts];

  fHistoClusOverlapHeadersGammaPt   = new TH1F*[fnCuts];
  fHistoClusOverlapMBHeaderGammaPt  = new TH1F*[fnCuts]; 
  fHistoClusAllHeadersGammaPt       = new TH1F*[fnCuts];
  fHistoClusRejectedHeadersGammaPt  = new TH1F*[fnCuts];


  for(int iCut = 0; iCut<fnCuts;iCut++){
    TString cutstringEvent    = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringCalo     = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
    TString cutstringMeson    = "NoMesonCut";
    if(fDoMesonAnalysis)
      cutstringMeson          = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

    fCutFolder[iCut]        = new TList();
    fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
    fCutFolder[iCut]->SetOwner(kTRUE);
    fOutputContainer->Add(fCutFolder[iCut]);
    fESDList[iCut]          = new TList();
    fESDList[iCut]->SetName(Form("%s_%s_%s ESD histograms", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
    fESDList[iCut]->SetOwner(kTRUE);
    fCutFolder[iCut]->Add(fESDList[iCut]);
    fQAList[iCut]          = new TList();
    fQAList[iCut]->SetName(Form("%s_%s_%s QA histograms", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
    fQAList[iCut]->SetOwner(kTRUE);
    fCutFolder[iCut]->Add(fQAList[iCut]);

    fHistoNEvents[iCut]     = new TH1F("NEvents", "NEvents", 15, -0.5, 14.5);
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
      TString TriggerNames  = "Not Trigger: ";
      TriggerNames          = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
      fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
    } else {
      fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
    }
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(7,"SPD Pile-Up");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL/TPC problems");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fHistoNEvents[iCut]->GetXaxis()->SetBinLabel(15,"Sphericity");
    fESDList[iCut]->Add(fHistoNEvents[iCut]);

    if (fIsMC > 1){
      fHistoNEventsWOWeight[iCut] = new TH1F("NEventsWOWeight", "NEventsWOWeight", 15, -0.5, 14.5);
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
      if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
        TString TriggerNames  = "Not Trigger: ";
        TriggerNames          = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
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
      fHistoNEventsWOWeight[iCut]->GetXaxis()->SetBinLabel(15,"Sphericity");
      fESDList[iCut]->Add(fHistoNEventsWOWeight[iCut]);
    }
    if (fIsMC == 2){
      fProfileJetJetXSection[iCut]  = new TProfile("XSection", "XSection", 1, -0.5, 0.5);
      fESDList[iCut]->Add(fProfileJetJetXSection[iCut]);
      fHistoJetJetNTrials[iCut]     = new TH1F("NTrials", "#sum{NTrials}", 1, 0, 1);
      fHistoJetJetNTrials[iCut]->GetXaxis()->SetBinLabel(1,"#sum{NTrials}");
      fESDList[iCut]->Add(fHistoJetJetNTrials[iCut]);
      fHistoPtHardJJWeight[iCut]     = new TH2F("fHistoPtHardJJWeight", "fHistoPtHardJJWeight", 400, 0, 200, 60, 0, 30);
      fESDList[iCut]->Add(fHistoPtHardJJWeight[iCut]);
    }

    if(fIsHeavyIon == 1)
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks", "GoodESDTracks", 4000, 0, 4000);
    else if(fIsHeavyIon == 2)
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks", "GoodESDTracks", 400, 0, 400);
    else
      fHistoNGoodESDTracks[iCut]    = new TH1F("GoodESDTracks", "GoodESDTracks", 200, 0, 200);
    fHistoNGoodESDTracks[iCut]->GetXaxis()->SetTitle("#primary tracks");
    fESDList[iCut]->Add(fHistoNGoodESDTracks[iCut]);

    fHistoVertexZ[iCut]             = new TH1F("VertexZ", "VertexZ", 200, -10, 10);
    fHistoVertexZ[iCut]->GetXaxis()->SetTitle("Z_{vtx} (cm)");
    fESDList[iCut]->Add(fHistoVertexZ[iCut]);

    if(fIsHeavyIon == 1)
      fHistoNGammaCandidatesBasic[iCut]  = new TH1F("GammaCandidatesBasic", "GammaCandidatesBasic", 600, 0, 600);
    else if(fIsHeavyIon == 2)
      fHistoNGammaCandidatesBasic[iCut]  = new TH1F("GammaCandidatesBasic", "GammaCandidatesBasic", 400, 0, 400);
    else
      fHistoNGammaCandidatesBasic[iCut]  = new TH1F("GammaCandidatesBasic", "GammaCandidatesBasic", 100, 0, 100);
    fHistoNGammaCandidatesBasic[iCut]->GetXaxis()->SetTitle("#cluster candidates basic");
    fESDList[iCut]->Add(fHistoNGammaCandidatesBasic[iCut]);


    if(fIsHeavyIon == 1)
      fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates", "GammaCandidates", 200, 0, 200);
    else if(fIsHeavyIon == 2)
      fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates", "GammaCandidates", 100, 0, 100);
    else
      fHistoNGammaCandidates[iCut]  = new TH1F("GammaCandidates", "GammaCandidates", 50, 0, 50);
    fHistoNGammaCandidates[iCut]->GetXaxis()->SetTitle("#cluster candidates with current cut");
    fESDList[iCut]->Add(fHistoNGammaCandidates[iCut]);

    if(!fDoLightOutput){
      if(fIsHeavyIon == 1)
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates", "GoodESDTracksVsGammaCandidates", 4000, 0, 4000, 200, 0, 200);
      else if(fIsHeavyIon == 2)
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates", "GoodESDTracksVsGammaCandidates", 400, 0, 400, 100, 0, 100);
      else
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]  = new TH2F("GoodESDTracksVsGammaCandidates", "GoodESDTracksVsGammaCandidates", 200, 0, 200, 50, 0, 50);
      fHistoNGoodESDTracksVsNGammaCandidates[iCut]->SetXTitle("#good tracks");
      fHistoNGoodESDTracksVsNGammaCandidates[iCut]->SetYTitle("#cluster candidates");
      fESDList[iCut]->Add(fHistoNGoodESDTracksVsNGammaCandidates[iCut]);

      fHistoSPDClusterTrackletBackground[iCut]        = new TH2F("SPD tracklets vs SPD clusters", "SPD tracklets vs SPD clusters", 100, 0, 200, 250, 0, 1000);
      fESDList[iCut]->Add(fHistoSPDClusterTrackletBackground[iCut]);

      if(fIsHeavyIon == 1)
        fHistoNV0Tracks[iCut]       = new TH1F("V0 Multiplicity", "V0 Multiplicity", 20000, 0, 40000);
      else if(fIsHeavyIon == 2)
        fHistoNV0Tracks[iCut]       = new TH1F("V0 Multiplicity", "V0 Multiplicity", 2500, 0, 2500);
      else
        fHistoNV0Tracks[iCut]       = new TH1F("V0 Multiplicity", "V0 Multiplicity", 1500, 0, 1500);
      fHistoNV0Tracks[iCut]->SetXTitle("V0 amplitude");
      fESDList[iCut]->Add(fHistoNV0Tracks[iCut]);
    }

    if(fIsHeavyIon==2) {
      fProfileEtaShift[iCut]        = new TProfile("Eta Shift", "Eta Shift", 1, -0.5, 0.5);
      fProfileEtaShift[iCut]->SetXTitle("#eta shift");
      fESDList[iCut]->Add(fProfileEtaShift[iCut]);
    }

    if (fIsMC > 1){
      fHistoNEvents[iCut]->Sumw2();
      fHistoNGoodESDTracks[iCut]->Sumw2();
      fHistoVertexZ[iCut]->Sumw2();
      fHistoNGammaCandidates[iCut]->Sumw2();
      fHistoNGammaCandidatesBasic[iCut]->Sumw2();
      if(!fDoLightOutput){
        fHistoNGoodESDTracksVsNGammaCandidates[iCut]->Sumw2();
        fHistoSPDClusterTrackletBackground[iCut]->Sumw2();
        fHistoNV0Tracks[iCut]->Sumw2();
      }
      if(fIsHeavyIon==2) fProfileEtaShift[iCut]->Sumw2();
    }

    fHistoClusGammaPt[iCut]               = new TH1F("ClusGamma_Pt", "ClusGamma_Pt", nBinsClusterPt, arrClusPtBinning);
    fHistoClusGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c)");
    fESDList[iCut]->Add(fHistoClusGammaPt[iCut]);
    fHistoClusGammaE[iCut]               = new TH1F("ClusGamma_E", "ClusGamma_E", nBinsClusterPt, arrClusPtBinning);
    fHistoClusGammaE[iCut]->SetXTitle("E_{clus} (GeV/c)");
    fESDList[iCut]->Add(fHistoClusGammaE[iCut]);
    fHistoClusOverlapHeadersGammaPt[iCut]   = new TH1F("ClusGammaOverlapHeaders_Pt", "ClusGammaOverlapHeaders_Pt", nBinsClusterPt, arrClusPtBinning);
    fHistoClusOverlapHeadersGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c), selected header w/ overlap");
    fESDList[iCut]->Add(fHistoClusOverlapHeadersGammaPt[iCut]);
    fHistoClusOverlapMBHeaderGammaPt[iCut]   = new TH1F("ClusGammaOverlapMBHeader_Pt", "ClusGammaOverlapMBHeader_Pt", nBinsClusterPt, arrClusPtBinning);
    fHistoClusOverlapMBHeaderGammaPt[iCut]->SetXTitle("#it{p}_{T,clus} (GeV/#it{c}), selected header w/ MB overlap");
    fESDList[iCut]->Add(fHistoClusOverlapMBHeaderGammaPt[iCut]);
    fHistoClusAllHeadersGammaPt[iCut]       = new TH1F("ClusGammaAllHeaders_Pt", "ClusGammaAllHeaders_Pt", nBinsClusterPt, arrClusPtBinning);
    fHistoClusAllHeadersGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c), all headers");
    fESDList[iCut]->Add(fHistoClusAllHeadersGammaPt[iCut]);
    fHistoClusRejectedHeadersGammaPt[iCut]  = new TH1F("ClusGammaRejectedHeaders_Pt", "ClusGammaRejectedHeaders_Pt", nBinsClusterPt, arrClusPtBinning);
    fHistoClusRejectedHeadersGammaPt[iCut]->SetXTitle("p_{T,clus} (GeV/c), rejected headers");
    fESDList[iCut]->Add(fHistoClusRejectedHeadersGammaPt[iCut]);

    if (fIsMC > 1){
      fHistoClusGammaPt[iCut]->Sumw2();
      fHistoClusGammaE[iCut]->Sumw2();
      fHistoClusOverlapHeadersGammaPt[iCut]->Sumw2();
      fHistoClusOverlapMBHeaderGammaPt[iCut]->Sumw2();
      fHistoClusAllHeadersGammaPt[iCut]->Sumw2();
      fHistoClusRejectedHeadersGammaPt[iCut]->Sumw2();
    }
    if(fDoMesonAnalysis){
      if(!fUseInOutPlane){
        fHistoMotherInvMassPtPhiV0A[iCut] = new TH3F("ESD_Mother_InvMass_Pt_PhiV0A", "ESD_Mother_InvMass_Pt_PhiV0A", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning, nPhiBins, arrPhiBin);
        fHistoMotherInvMassPtPhiV0A[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherInvMassPtPhiV0A[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoMotherInvMassPtPhiV0A[iCut]->SetZTitle("phi");
        fESDList[iCut]->Add(fHistoMotherInvMassPtPhiV0A[iCut]);
        fHistoMotherInvMassPtPhiV0C[iCut] = new TH3F("ESD_Mother_InvMass_Pt_PhiV0C", "ESD_Mother_InvMass_Pt_PhiV0C", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning, nPhiBins, arrPhiBin);
        fHistoMotherInvMassPtPhiV0C[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherInvMassPtPhiV0C[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoMotherInvMassPtPhiV0C[iCut]->SetZTitle("phi");
        fESDList[iCut]->Add(fHistoMotherInvMassPtPhiV0C[iCut]);
        fHistoMotherBackInvMassPtPhiV0A[iCut] = new TH3F("ESD_Background_InvMass_Pt_PhiV0A", "ESD_Background_InvMass_Pt_PhiV0A", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning, nPhiBins, arrPhiBin);
        fHistoMotherBackInvMassPtPhiV0A[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherBackInvMassPtPhiV0A[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoMotherBackInvMassPtPhiV0A[iCut]->SetZTitle("phi");
        fESDList[iCut]->Add(fHistoMotherBackInvMassPtPhiV0A[iCut]);
        fHistoMotherBackInvMassPtPhiV0C[iCut] = new TH3F("ESD_Background_InvMass_Pt_PhiV0C", "ESD_Background_InvMass_Pt_PhiV0C", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning, nPhiBins, arrPhiBin);
        fHistoMotherBackInvMassPtPhiV0C[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherBackInvMassPtPhiV0C[iCut]->SetYTitle("p_{T} (GeV/c)");
        fHistoMotherBackInvMassPtPhiV0C[iCut]->SetZTitle("phi");
        fESDList[iCut]->Add(fHistoMotherBackInvMassPtPhiV0C[iCut]);
      } else{
        fHistoMotherInvMassPtV0CInPlane[iCut] = new TH2F("ESD_Mother_InvMass_Pt_V0C_InPlane", "ESD_Mother_InvMass_Pt_V0C_InPlane", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherInvMassPtV0CInPlane[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherInvMassPtV0CInPlane[iCut]->SetYTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherInvMassPtV0CInPlane[iCut]);
        fHistoMotherInvMassPtV0AInPlane[iCut] = new TH2F("ESD_Mother_InvMass_Pt_V0A_InPlane", "ESD_Mother_InvMass_Pt_V0A_InPlane", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherInvMassPtV0AInPlane[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherInvMassPtV0AInPlane[iCut]->SetYTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherInvMassPtV0AInPlane[iCut]);
        fHistoMotherInvMassPtV0COutPlane[iCut] = new TH2F("ESD_Mother_InvMass_Pt_V0C_OutPlane", "ESD_Mother_InvMass_Pt_V0C_OutPlane", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherInvMassPtV0COutPlane[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherInvMassPtV0COutPlane[iCut]->SetYTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherInvMassPtV0COutPlane[iCut]);
        fHistoMotherInvMassPtV0AOutPlane[iCut] = new TH2F("ESD_Mother_InvMass_Pt_V0A_OutPlane", "ESD_Mother_InvMass_Pt_V0A_OutPlane", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherInvMassPtV0AOutPlane[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherInvMassPtV0AOutPlane[iCut]->SetYTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherInvMassPtV0AOutPlane[iCut]);
        fHistoMotherBackInvMassPtV0CInPlane[iCut] = new TH2F("ESD_Background_InvMass_Pt_V0C_InPlane", "ESD_Background_InvMass_Pt_V0C_InPlane", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherBackInvMassPtV0CInPlane[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherBackInvMassPtV0CInPlane[iCut]->SetYTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherBackInvMassPtV0CInPlane[iCut]);
        fHistoMotherBackInvMassPtV0AInPlane[iCut] = new TH2F("ESD_Background_InvMass_Pt_V0A_InPlane", "ESD_Background_InvMass_Pt_V0A_InPlane", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherBackInvMassPtV0AInPlane[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherBackInvMassPtV0AInPlane[iCut]->SetYTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherBackInvMassPtV0AInPlane[iCut]);
        fHistoMotherBackInvMassPtV0COutPlane[iCut] = new TH2F("ESD_Background_InvMass_Pt_V0C_OutPlane", "ESD_Background_InvMass_Pt_V0C_OutPlane", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherBackInvMassPtV0COutPlane[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherBackInvMassPtV0COutPlane[iCut]->SetYTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherBackInvMassPtV0COutPlane[iCut]);
        fHistoMotherBackInvMassPtV0AOutPlane[iCut] = new TH2F("ESD_Background_InvMass_Pt_V0A_OutPlane", "ESD_Background_InvMass_Pt_V0A_OutPlane", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
        fHistoMotherBackInvMassPtV0AOutPlane[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoMotherBackInvMassPtV0AOutPlane[iCut]->SetYTitle("p_{T} (GeV/c)");
        fESDList[iCut]->Add(fHistoMotherBackInvMassPtV0AOutPlane[iCut]);
      }
      fHistoMotherInvMassPtV0CCos2phi[iCut] = new TProfile2D("InvMassPtV0CCos2phi", "InvMassPtV0CCos2phi", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
      fHistoMotherInvMassPtV0ACos2phi[iCut] = new TProfile2D("InvMassPtV0ACos2phi", "InvMassPtV0ACos2phi", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
      fHistoMotherBackInvMassPtV0CCos2phi[iCut] = new TProfile2D("BackInvMassPtV0CCos2phi", "BackInvMassPtV0CCos2phi", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
      fHistoMotherBackInvMassPtV0ACos2phi[iCut] = new TProfile2D("BackInvMassPtV0ACos2phi", "BackInvMassPtV0ACos2phi", nBinsMinv, arrMinvBin, nBinsPt, arrPtBinning);
      fESDList[iCut]->Add(fHistoMotherInvMassPtV0CCos2phi[iCut]);
      fESDList[iCut]->Add(fHistoMotherInvMassPtV0ACos2phi[iCut]);
      fESDList[iCut]->Add(fHistoMotherBackInvMassPtV0CCos2phi[iCut]);
      fESDList[iCut]->Add(fHistoMotherBackInvMassPtV0ACos2phi[iCut]);
      // VZERO
      fHist2DPsi2V0CCent[iCut] = new TH2D("fHist2DPsi2V0CCent", "", 20, 0, 100, 100, -2 * TMath::Pi(), 2 * TMath::Pi());
      fHist2DPsi2V0ACent[iCut] = new TH2D("fHist2DPsi2V0ACent", "", 20, 0, 100, 100, -2 * TMath::Pi(), 2 * TMath::Pi());
      fHist2V0Res[iCut] = new TProfile("fHist2V0Res", "", 4, 0, 4);
      fESDList[iCut]->Add(fHist2DPsi2V0CCent[iCut]);
      fESDList[iCut]->Add(fHist2DPsi2V0ACent[iCut]);
      fESDList[iCut]->Add(fHist2V0Res[iCut]);
      if(IsQAVZERO){
        fProfileV0CQxCentGE[iCut] = new TProfile("fProfileV0CQxCentGE", "", 100, 0, 100.);
        fProfileV0CQyCentGE[iCut] = new TProfile("fProfileV0CQyCentGE", "", 100, 0, 100.);
        fProfileV0CQxVtxGE[iCut] = new TProfile("fProfileV0CQxVzGE", "", 20, -10, 10);
        fProfileV0CQyVtxGE[iCut] = new TProfile("fProfileV0CQyVzGE", "", 20, -10, 10);
        fHist2CalibPsi2V0CCentGE[iCut] = new TH2D("fHist2CalibPsi2V0CCentGE", "", 20, 0, 100, 50, 0, TMath::Pi());
        fProfileV0AQxCentGE[iCut] = new TProfile("fProfileV0AQxCentGE", "", 100, 0, 100.);
        fProfileV0AQyCentGE[iCut] = new TProfile("fProfileV0AQyCentGE", "", 100, 0, 100.);
        fProfileV0AQxVtxGE[iCut] = new TProfile("fProfileV0AQxVzGE", "", 20, -10, 10);
        fProfileV0AQyVtxGE[iCut] = new TProfile("fProfileV0AQyVzGE", "", 20, -10, 10);
        fHist2CalibPsi2V0ACentGE[iCut] = new TH2D("fHist2CalibPsi2V0ACentGE", "", 20, 0, 100, 50, 0, TMath::Pi());
        fProfileV0CQxCentRC[iCut] = new TProfile("fProfileV0CQxCentRC", "", 100, 0, 100.);
        fProfileV0CQyCentRC[iCut] = new TProfile("fProfileV0CQyCentRC", "", 100, 0, 100.);
        fProfileV0CQxVtxRC[iCut] = new TProfile("fProfileV0CQxVzRC", "", 20, -10, 10);
        fProfileV0CQyVtxRC[iCut] = new TProfile("fProfileV0CQyVzRC", "", 20, -10, 10);
        fHist2CalibPsi2V0CCentRC[iCut] = new TH2D("fHist2CalibPsi2V0CCentRC", "", 20, 0, 100, 50, 0, TMath::Pi());
        fProfileV0AQxCentRC[iCut] = new TProfile("fProfileV0AQxCentRC", "", 100, 0, 100.);
        fProfileV0AQyCentRC[iCut] = new TProfile("fProfileV0AQyCentRC", "", 100, 0, 100.);
        fProfileV0AQxVtxRC[iCut] = new TProfile("fProfileV0AQxVzRC", "", 20, -10, 10);
        fProfileV0AQyVtxRC[iCut] = new TProfile("fProfileV0AQyVzRC", "", 20, -10, 10);
        fHist2CalibPsi2V0ACentRC[iCut] = new TH2D("fHist2CalibPsi2V0ACentRC", "", 20, 0, 100, 50, 0, TMath::Pi());
        fQAList[iCut]->Add(fProfileV0CQxCentGE[iCut]);
        fQAList[iCut]->Add(fProfileV0CQyCentGE[iCut]);
        fQAList[iCut]->Add(fProfileV0CQxVtxGE[iCut]);
        fQAList[iCut]->Add(fProfileV0CQyVtxGE[iCut]);
        fQAList[iCut]->Add(fHist2CalibPsi2V0CCentGE[iCut]);
        fQAList[iCut]->Add(fProfileV0AQxCentGE[iCut]);
        fQAList[iCut]->Add(fProfileV0AQyCentGE[iCut]);
        fQAList[iCut]->Add(fProfileV0AQxVtxGE[iCut]);
        fQAList[iCut]->Add(fProfileV0AQyVtxGE[iCut]);
        fQAList[iCut]->Add(fHist2CalibPsi2V0ACentGE[iCut]);
        fQAList[iCut]->Add(fProfileV0CQxCentRC[iCut]);
        fQAList[iCut]->Add(fProfileV0CQyCentRC[iCut]);
        fQAList[iCut]->Add(fProfileV0CQxVtxRC[iCut]);
        fQAList[iCut]->Add(fProfileV0CQyVtxRC[iCut]);
        fQAList[iCut]->Add(fHist2CalibPsi2V0CCentRC[iCut]);
        fQAList[iCut]->Add(fProfileV0AQxCentRC[iCut]);
        fQAList[iCut]->Add(fProfileV0AQyCentRC[iCut]);
        fQAList[iCut]->Add(fProfileV0AQxVtxRC[iCut]);
        fQAList[iCut]->Add(fProfileV0AQyVtxRC[iCut]);
        fQAList[iCut]->Add(fHist2CalibPsi2V0ACentRC[iCut]);
      }

      if (fIsMC > 1){
        if(!fUseInOutPlane){
          fHistoMotherInvMassPtPhiV0A[iCut]->Sumw2();
          fHistoMotherInvMassPtPhiV0C[iCut]->Sumw2();
          fHistoMotherBackInvMassPtPhiV0A[iCut]->Sumw2();
          fHistoMotherBackInvMassPtPhiV0C[iCut]->Sumw2();
        } else{
          fHistoMotherInvMassPtV0CInPlane[iCut]->Sumw2();
          fHistoMotherInvMassPtV0AInPlane[iCut]->Sumw2();
          fHistoMotherInvMassPtV0COutPlane[iCut]->Sumw2();
          fHistoMotherInvMassPtV0AOutPlane[iCut]->Sumw2();
          fHistoMotherBackInvMassPtV0CInPlane[iCut]->Sumw2();
          fHistoMotherBackInvMassPtV0AInPlane[iCut]->Sumw2();
          fHistoMotherBackInvMassPtV0COutPlane[iCut]->Sumw2();
          fHistoMotherBackInvMassPtV0AOutPlane[iCut]->Sumw2();
        }
        fHistoMotherInvMassPtV0CCos2phi[iCut]->Sumw2();
        fHistoMotherInvMassPtV0ACos2phi[iCut]->Sumw2();
        fHistoMotherBackInvMassPtV0CCos2phi[iCut]->Sumw2();
        fHistoMotherBackInvMassPtV0ACos2phi[iCut]->Sumw2();
      }
    }

    if( ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoSectorMixing() ){
      fV0Reader->SetCalcSector(kTRUE);
    }
  }
  if(fDoMesonAnalysis){
    InitBack(); // Init Background Handler
  }


  if(fIsMC> 0){
    // MC Histogramms
    fMCList     = new TList*[fnCuts];
    // True Histogramms
    fTrueList   = new TList*[fnCuts];
    // Selected Header List

    fHistoMCAllGammaPt              = new TH1F*[fnCuts];
    fHistoMCGammaPtNotTriggered     = new TH1F*[fnCuts];
    fHistoMCGammaPtNoVertex         = new TH1F*[fnCuts];
    
    fHistoMCEventsTrigg             = new TH1D*[fnCuts];

    if(!fDoLightOutput){
      fHistoMCHeaders                 = new TH1I*[fnCuts];
      fHistoMCAllSecondaryGammaPt     = new TH2F*[fnCuts];
      fHistoMCDecayGammaPi0Pt         = new TH1F*[fnCuts];
      fHistoMCDecayGammaRhoPt         = new TH1F*[fnCuts];
      fHistoMCDecayGammaEtaPt         = new TH1F*[fnCuts];
      fHistoMCDecayGammaOmegaPt       = new TH1F*[fnCuts];
      fHistoMCDecayGammaEtapPt        = new TH1F*[fnCuts];
      fHistoMCDecayGammaPhiPt         = new TH1F*[fnCuts];
      fHistoMCDecayGammaSigmaPt       = new TH1F*[fnCuts];
      fHistoClusPhotonBGPt            = new TH2F*[fnCuts];
      fHistoClusPhotonPlusConvBGPt    = new TH2F*[fnCuts];
    }

    fHistoTrueClusGammaPt                                       = new TH1F*[fnCuts];
    if(!fDoLightOutput){
      fHistoTruePrimaryClusGammaPt                              = new TH1F*[fnCuts];
      fHistoTruePrimaryClusGammaESDPtMCPt                       = new TH2F*[fnCuts];
      fHistoTruePrimaryClusConvGammaPt                          = new TH1F*[fnCuts];
      fHistoTruePrimaryClusConvGammaESDPtMCPt                   = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusGammaPt                            = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusConvGammaPt                        = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusGammaMCPt                          = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusConvGammaMCPt                      = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt         = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt     = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt         = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt     = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt      = new TH2F*[fnCuts];
      fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt  = new TH2F*[fnCuts];
      fHistoTrueNLabelsInClus                                   = new TH1F*[fnCuts];
    }
    fHistoDoubleCountTrueClusterGammaPt                         = new TH2F*[fnCuts];
    fHistoMultipleCountTrueClusterGamma                         = new TH1F*[fnCuts];
    fHistoTrueClusGammaEResE                                    = new TH2F*[fnCuts];
    fHistoTrueClusPhotonGammaEResE                              = new TH2F*[fnCuts];

    if(fDoMesonAnalysis){
      fHistoMCPi0Pt             = new TH1F*[fnCuts];
      fHistoMCPi0PtNotTriggered = new TH1F*[fnCuts];
      fHistoMCPi0PtNoVertex     = new TH1F*[fnCuts];
      fHistoMCPi0WOWeightPt     = new TH1F*[fnCuts];
      fHistoMCPi0InAccPt        = new TH1F*[fnCuts];
      fHistoMCPi0InAccPtNotTriggered   = new TH1F*[fnCuts];

      if (fIsMC > 1){
        fHistoMCPi0WOEvtWeightPt       = new TH1F*[fnCuts];
        fHistoMCPi0WOEvtWeightInAccPt  = new TH1F*[fnCuts];
        if( !fDoPi0Only ){
          fHistoMCEtaWOEvtWeightPt       = new TH1F*[fnCuts];
          fHistoMCEtaWOEvtWeightInAccPt  = new TH1F*[fnCuts];
        }
      }

      fHistoTruePi0InvMassPt                    = new TH2F*[fnCuts];
      if (!fDoLightOutput){
          fHistoTruePi0InvMassPtAdditional            = new TH2F*[fnCuts];
      }
      fHistoDoubleCountTruePi0InvMassPt         = new TH2F*[fnCuts];
      fHistoTruePrimaryPi0InvMassPt             = new TH2F*[fnCuts];
      fHistoTruePrimaryPi0W0WeightingInvMassPt  = new TH2F*[fnCuts];
      fProfileTruePrimaryPi0WeightsInvMassPt    = new TProfile2D*[fnCuts];
      fHistoTrueSecondaryPi0InvMassPt           = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromK0sInvMassPt    = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromK0lInvMassPt    = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromLambdaInvMassPt = new TH2F*[fnCuts];
      fHistoTrueSecondaryPi0FromEtaInvMassPt    = new TH2F*[fnCuts];
      if( !fDoPi0Only ){
        fHistoMCEtaPt                             = new TH1F*[fnCuts];
        fHistoMCEtaPtNotTriggered                 = new TH1F*[fnCuts];
        fHistoMCEtaPtNoVertex                     = new TH1F*[fnCuts];
        fHistoMCEtaWOWeightPt                     = new TH1F*[fnCuts];
        fHistoMCEtaInAccPt                        = new TH1F*[fnCuts];
        fHistoMCEtaInAccPtNotTriggered            = new TH1F*[fnCuts];
        fHistoTrueEtaInvMassPt                    = new TH2F*[fnCuts];
        if (!fDoLightOutput){
            fHistoTrueEtaInvMassPtAdditional            = new TH2F*[fnCuts];
        }
        fHistoDoubleCountTrueEtaInvMassPt         = new TH2F*[fnCuts];
        fHistoTruePrimaryEtaInvMassPt             = new TH2F*[fnCuts];
        fHistoTruePrimaryEtaW0WeightingInvMassPt  = new TH2F*[fnCuts];
        fProfileTruePrimaryEtaWeightsInvMassPt    = new TProfile2D*[fnCuts];
        fHistoMCSecEtaPt                          = new TH1F*[fnCuts];
        fHistoMCSecEtaSource                      = new TH1F*[fnCuts];
      }
      if(!fDoLightOutput){
        fHistoTruePi0InvMassPtAlpha               = new TH2F*[fnCuts];
        fHistoTruePi0PureGammaInvMassPtAlpha      = new TH2F*[fnCuts];
      }
      fHistoMCPrimaryPtvsSource                   = new TH2F*[fnCuts];
      fHistoMCSecPi0PtvsSource                    = new TH2F*[fnCuts];
      fHistoMCSecPi0InAccPtvsSource               = new TH2F*[fnCuts];
      fHistoMCSecPi0Source                        = new TH1F*[fnCuts];
    }


    for(int iCut = 0; iCut<fnCuts;iCut++){
      TString cutstringEvent    = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
      TString cutstringCalo     = ((AliCaloPhotonCuts*)fClusterCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson    = "NoMesonCut";
      if(fDoMesonAnalysis)
        cutstringMeson          = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

      fMCList[iCut]                   = new TList();
      fMCList[iCut]->SetName(Form("%s_%s_%s MC histograms", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
      fMCList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fMCList[iCut]);

      fHistoMCAllGammaPt[iCut]          = new TH1F("MC_AllGamma_Pt", "MC_AllGamma_Pt", nBinsClusterPt, arrClusPtBinning);
      fHistoMCAllGammaPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fMCList[iCut]->Add(fHistoMCAllGammaPt[iCut]);
      fHistoMCGammaPtNotTriggered[iCut]          = new TH1F("MC_AllGammaNotTriggered_Pt", "MC_AllGammaNotTriggered_Pt", nBinsClusterPt, arrClusPtBinning);
      fHistoMCGammaPtNotTriggered[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fMCList[iCut]->Add(fHistoMCGammaPtNotTriggered[iCut]);
      fHistoMCGammaPtNoVertex[iCut]          = new TH1F("MC_AllGammaNoVertex_Pt", "MC_AllGammaNoVertex_Pt", nBinsClusterPt, arrClusPtBinning);
      fHistoMCGammaPtNoVertex[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fMCList[iCut]->Add(fHistoMCGammaPtNoVertex[iCut]);

      fHistoMCEventsTrigg[iCut]          = new TH1D("MC_NEvents", "MC_NEvents", 4, -0.5, 3.5);
      fHistoMCEventsTrigg[iCut]->GetXaxis()->SetBinLabel(1, "accepted, trig");
      fHistoMCEventsTrigg[iCut]->GetXaxis()->SetBinLabel(2, "accepted, not trig");
      fHistoMCEventsTrigg[iCut]->GetXaxis()->SetBinLabel(3, "rejected, trig");
      fHistoMCEventsTrigg[iCut]->GetXaxis()->SetBinLabel(4, "rejected, not trig");
      fMCList[iCut]->Add(fHistoMCEventsTrigg[iCut]);

      if(!fDoLightOutput){
        fHistoMCHeaders[iCut]             = new TH1I("MC_Headers", "MC_Headers", 20, 0, 20);
        fHistoMCHeaders[iCut]->SetXTitle("accepted headers");
        fMCList[iCut]->Add(fHistoMCHeaders[iCut]);
        fHistoMCAllSecondaryGammaPt[iCut] = new TH2F("MC_AllSecondaryGamma_Pt", "MC_AllSecondaryGamma_Pt", nBinsClusterPt, arrClusPtBinning, 5, -0.5, 4.5);
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(1,"K0s");
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(2,"K0l");
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(3,"Lambda");
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(4,"Eta");
        fHistoMCAllSecondaryGammaPt[iCut]->GetYaxis()->SetBinLabel(5,"rest");
        fHistoMCAllSecondaryGammaPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoMCAllSecondaryGammaPt[iCut]->SetYTitle("sec. particle");
        fMCList[iCut]->Add(fHistoMCAllSecondaryGammaPt[iCut]);
        fHistoMCDecayGammaPi0Pt[iCut]     = new TH1F("MC_DecayGammaPi0_Pt", "MC_DecayGammaPi0_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoMCDecayGammaPi0Pt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fMCList[iCut]->Add(fHistoMCDecayGammaPi0Pt[iCut]);
        fHistoMCDecayGammaRhoPt[iCut]     = new TH1F("MC_DecayGammaRho_Pt", "MC_DecayGammaRho_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoMCDecayGammaRhoPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fMCList[iCut]->Add(fHistoMCDecayGammaRhoPt[iCut]);
        fHistoMCDecayGammaEtaPt[iCut]     = new TH1F("MC_DecayGammaEta_Pt", "MC_DecayGammaEta_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoMCDecayGammaEtaPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fMCList[iCut]->Add(fHistoMCDecayGammaEtaPt[iCut]);
        fHistoMCDecayGammaOmegaPt[iCut]   = new TH1F("MC_DecayGammaOmega_Pt", "MC_DecayGammaOmmega_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoMCDecayGammaOmegaPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fMCList[iCut]->Add(fHistoMCDecayGammaOmegaPt[iCut]);
        fHistoMCDecayGammaEtapPt[iCut]    = new TH1F("MC_DecayGammaEtap_Pt", "MC_DecayGammaEtap_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoMCDecayGammaEtapPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fMCList[iCut]->Add(fHistoMCDecayGammaEtapPt[iCut]);
        fHistoMCDecayGammaPhiPt[iCut]     = new TH1F("MC_DecayGammaPhi_Pt", "MC_DecayGammaPhi_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoMCDecayGammaPhiPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fMCList[iCut]->Add(fHistoMCDecayGammaPhiPt[iCut]);
        fHistoMCDecayGammaSigmaPt[iCut]   = new TH1F("MC_DecayGammaSigma_Pt", "MC_DecayGammaSigma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoMCDecayGammaSigmaPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fMCList[iCut]->Add(fHistoMCDecayGammaSigmaPt[iCut]);

        if (fIsMC > 1){
          fHistoMCAllSecondaryGammaPt[iCut]->Sumw2();
          fHistoMCDecayGammaPi0Pt[iCut]->Sumw2();
          fHistoMCDecayGammaRhoPt[iCut]->Sumw2();
          fHistoMCDecayGammaEtaPt[iCut]->Sumw2();
          fHistoMCDecayGammaOmegaPt[iCut]->Sumw2();
          fHistoMCDecayGammaEtapPt[iCut]->Sumw2();
          fHistoMCDecayGammaPhiPt[iCut]->Sumw2();
          fHistoMCDecayGammaSigmaPt[iCut]->Sumw2();
        }
      }
      if (fIsMC > 1){
        fHistoMCAllGammaPt[iCut]->Sumw2();
        fHistoMCGammaPtNotTriggered[iCut]->Sumw2();
        fHistoMCGammaPtNoVertex[iCut]->Sumw2();
        fHistoMCEventsTrigg[iCut]->Sumw2();
      }
      if(fDoMesonAnalysis){
        fHistoMCPi0Pt[iCut]           = new TH1F("MC_Pi0_Pt", "MC_Pi0_Pt", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoMCPi0Pt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoMCPi0Pt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0Pt[iCut]);
        fHistoMCPi0PtNotTriggered[iCut]           = new TH1F("MC_Pi0_Pt_NotTriggered", "MC_Pi0_Pt_NotTriggered", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoMCPi0PtNotTriggered[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoMCPi0PtNotTriggered[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0PtNotTriggered[iCut]);
        fHistoMCPi0PtNoVertex[iCut]           = new TH1F("MC_Pi0_Pt_NoVertex", "MC_Pi0_Pt_NoVertex", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoMCPi0PtNoVertex[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoMCPi0PtNoVertex[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0PtNoVertex[iCut]);
        fHistoMCPi0WOWeightPt[iCut]   = new TH1F("MC_Pi0_WOWeights_Pt", "MC_Pi0_WOWeights_Pt", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoMCPi0WOWeightPt[iCut]->Sumw2();
        fHistoMCPi0WOWeightPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fMCList[iCut]->Add(fHistoMCPi0WOWeightPt[iCut]);
        fHistoMCPi0InAccPt[iCut]      = new TH1F("MC_Pi0InAcc_Pt", "MC_Pi0InAcc_Pt", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoMCPi0InAccPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoMCPi0InAccPt[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0InAccPt[iCut]);
        fHistoMCPi0InAccPtNotTriggered[iCut]      = new TH1F("MC_Pi0InAcc_Pt_NotTriggered", "MC_Pi0InAcc_Pt_NotTriggered", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
        fHistoMCPi0InAccPtNotTriggered[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoMCPi0InAccPtNotTriggered[iCut]->Sumw2();
        fMCList[iCut]->Add(fHistoMCPi0InAccPtNotTriggered[iCut]);
        if( !fDoPi0Only ){
          fHistoMCEtaPt[iCut]           = new TH1F("MC_Eta_Pt", "MC_Eta_Pt", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fHistoMCEtaPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fHistoMCEtaPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCEtaPt[iCut]);
          fHistoMCEtaPtNotTriggered[iCut]           = new TH1F("MC_Eta_Pt_NotTriggered", "MC_Eta_Pt_NotTriggered", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fHistoMCEtaPtNotTriggered[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fHistoMCEtaPtNotTriggered[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCEtaPtNotTriggered[iCut]);
          fHistoMCEtaPtNoVertex[iCut]           = new TH1F("MC_Eta_Pt_NoVertex", "MC_Eta_Pt_NoVertex", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fHistoMCEtaPtNoVertex[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fHistoMCEtaPtNoVertex[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCEtaPtNoVertex[iCut]);
          fHistoMCEtaWOWeightPt[iCut]   = new TH1F("MC_Eta_WOWeights_Pt", "MC_Eta_WOWeights_Pt", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fHistoMCEtaWOWeightPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fHistoMCEtaWOWeightPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCEtaWOWeightPt[iCut]);
          fHistoMCEtaInAccPt[iCut]      = new TH1F("MC_EtaInAcc_Pt", "MC_EtaInAcc_Pt", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fHistoMCEtaInAccPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fHistoMCEtaInAccPt[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCEtaInAccPt[iCut]);
          fHistoMCEtaInAccPtNotTriggered[iCut]      = new TH1F("MC_EtaInAcc_Pt_NotTriggered", "MC_EtaInAcc_Pt_NotTriggered", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fHistoMCEtaInAccPtNotTriggered[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fHistoMCEtaInAccPtNotTriggered[iCut]->Sumw2();
          fMCList[iCut]->Add(fHistoMCEtaInAccPtNotTriggered[iCut]);
        }
        if (fIsMC > 1){
          fHistoMCPi0WOWeightPt[iCut]->Sumw2();
          fHistoMCPi0WOEvtWeightPt[iCut] = new TH1F("MC_Pi0_WOEventWeights_Pt", "MC_Pi0_WOEventWeights_Pt", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fHistoMCPi0WOEvtWeightPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightPt[iCut]);
          fHistoMCPi0WOEvtWeightInAccPt[iCut] = new TH1F("MC_Pi0WOEvtWeightInAcc_Pt", "MC_Pi0WOEvtWeightInAcc_Pt", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fHistoMCPi0WOEvtWeightInAccPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fMCList[iCut]->Add(fHistoMCPi0WOEvtWeightInAccPt[iCut]);
          if( !fDoPi0Only ){
            fHistoMCEtaWOWeightPt[iCut]->Sumw2();
            fHistoMCEtaWOEvtWeightPt[iCut] = new TH1F("MC_Eta_WOEventWeights_Pt", "MC_Eta_WOEventWeights_Pt", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
            fHistoMCEtaWOEvtWeightPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fMCList[iCut]->Add(fHistoMCEtaWOEvtWeightPt[iCut]);
            fHistoMCEtaWOEvtWeightInAccPt[iCut] = new TH1F("MC_EtaWOEvtWeightInAcc_Pt", "MC_EtaWOEvtWeightInAcc_Pt", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
            fHistoMCEtaWOEvtWeightInAccPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
            fMCList[iCut]->Add(fHistoMCEtaWOEvtWeightInAccPt[iCut]);
          }

        }
        fHistoMCPrimaryPtvsSource[iCut]   = new TH2F("MC_Primary_Pt_Source", "MC_Primary_Pt_Source", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt, 7, -0.5, 6.5);
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(1,"Pi+");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(2,"Pi-");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(3,"K+");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(4,"K-");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(5,"K0s");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(6,"K0l");
        fHistoMCPrimaryPtvsSource[iCut]->GetYaxis()->SetBinLabel(7,"Lambda");
        fHistoMCPrimaryPtvsSource[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoMCPrimaryPtvsSource[iCut]->SetYTitle("particle");
        fMCList[iCut]->Add(fHistoMCPrimaryPtvsSource[iCut]);

        fHistoMCSecPi0Source[iCut]      = new TH1F("MC_SecPi0_Source", "MC_SecPi0_Source", 5000, 0., 5000);
        fHistoMCSecPi0Source[iCut]->SetYTitle("source PDG");
        fMCList[iCut]->Add(fHistoMCSecPi0Source[iCut]);
        if( !fDoPi0Only ){
          fHistoMCSecEtaSource[iCut]      = new TH1F("MC_SecEta_Source", "MC_SecEta_Source", 5000, 0, 5000);
          fHistoMCSecEtaSource[iCut]->SetYTitle("source PDG");
          fMCList[iCut]->Add(fHistoMCSecEtaSource[iCut]);
        }
        fHistoMCSecPi0PtvsSource[iCut]  = new TH2F("MC_SecPi0_Pt_Source", "MC_SecPi0_Pt_Source", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt, 16, -0.5, 15.5);
        fHistoMCSecPi0PtvsSource[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoMCSecPi0PtvsSource[iCut]->SetYTitle("source");
        fMCList[iCut]->Add(fHistoMCSecPi0PtvsSource[iCut]);
        fHistoMCSecPi0InAccPtvsSource[iCut]  = new TH2F("MC_SecPi0InAcc_Pt_Source", "MC_SecPi0InAcc_Pt_Source", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt,  16, -0.5, 15.5);
        fHistoMCSecPi0InAccPtvsSource[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoMCSecPi0InAccPtvsSource[iCut]->SetYTitle("source");
        fMCList[iCut]->Add(fHistoMCSecPi0InAccPtvsSource[iCut]);
        if( !fDoPi0Only ){
          fHistoMCSecEtaPt[iCut]          = new TH1F("MC_SecEta_Pt", "MC_SecEta_Pt", (int)((maxPt-minPt)/binWidthPt), minPt, maxPt);
          fHistoMCSecEtaPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
          fMCList[iCut]->Add(fHistoMCSecEtaPt[iCut]);
        }
        if (fIsMC == 2) {
          fHistoMCPrimaryPtvsSource[iCut]->Sumw2();
          fHistoMCSecPi0PtvsSource[iCut]->Sumw2();
          fHistoMCSecPi0InAccPtvsSource[iCut]->Sumw2();
          if( !fDoPi0Only ) fHistoMCSecEtaPt[iCut]->Sumw2();
        }
      }
      fTrueList[iCut]                 = new TList();
      fTrueList[iCut]->SetName(Form("%s_%s_%s True histograms", cutstringEvent.Data(), cutstringCalo.Data(), cutstringMeson.Data()));
      fTrueList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fTrueList[iCut]);

      if(!fDoLightOutput){
        fHistoClusPhotonBGPt[iCut]          = new TH2F("ESD_TrueClusPhotonBG_Pt", "ESD_TrueClusPhotonBG_Pt", nBinsClusterPt, arrClusPtBinning,10,-0.5,9.5);
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 1,"Elec");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 2,"Pion");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 3,"Proton");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 4,"Kaon");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 5,"Neutron");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 6,"K0s");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 7,"Lambda");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 8,"Muon");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 9,"K0l");
        fHistoClusPhotonBGPt[iCut]->GetYaxis()->SetBinLabel( 10,"Rest");
        fHistoClusPhotonBGPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoClusPhotonBGPt[iCut]->SetYTitle("source");
        fTrueList[iCut]->Add(fHistoClusPhotonBGPt[iCut]);
        fHistoClusPhotonPlusConvBGPt[iCut]  = new TH2F("ESD_TrueClusPhotonPlusConvBG_Pt", "ESD_TrueClusPhotonPlusConvBG_Pt", nBinsClusterPt, arrClusPtBinning,10,-0.5,9.5);
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 1,"Elec");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 2,"Pion");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 3,"Proton");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 4,"Kaon");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 5,"Neutron");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 6,"K0s");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 7,"Lambda");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 8,"Muon");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel( 9,"K0l");
        fHistoClusPhotonPlusConvBGPt[iCut]->GetYaxis()->SetBinLabel(10,"Rest");
        fHistoClusPhotonPlusConvBGPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoClusPhotonPlusConvBGPt[iCut]->SetYTitle("source");
        fTrueList[iCut]->Add(fHistoClusPhotonPlusConvBGPt[iCut]);

      }
      fHistoTrueClusGammaPt[iCut]                   = new TH1F("TrueClusGamma_Pt", "ESD_TrueClusGamma_Pt", nBinsClusterPt, arrClusPtBinning);
      fHistoTrueClusGammaPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fTrueList[iCut]->Add(fHistoTrueClusGammaPt[iCut]);
      if(!fDoLightOutput){
        fHistoTruePrimaryClusGammaPt[iCut]            = new TH1F("TruePrimaryClusGamma_Pt", "ESD_TruePrimaryClusGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTruePrimaryClusGammaPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fTrueList[iCut]->Add(fHistoTruePrimaryClusGammaPt[iCut]);
        fHistoTruePrimaryClusGammaESDPtMCPt[iCut]     = new TH2F("TruePrimaryClusGamma_Pt_MCPt", "ESD_TruePrimaryClusGamma_MCPt", nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTruePrimaryClusGammaESDPtMCPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoTruePrimaryClusGammaESDPtMCPt[iCut]->SetYTitle("p_{T, MC} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTruePrimaryClusGammaESDPtMCPt[iCut]);
        fHistoTruePrimaryClusConvGammaPt[iCut]        = new TH1F("TruePrimaryClusConvGamma_Pt", "ESD_TruePrimaryClusConvGamma_Pt", nBinsClusterPt, arrClusPtBinning);
        fHistoTruePrimaryClusConvGammaPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fTrueList[iCut]->Add(fHistoTruePrimaryClusConvGammaPt[iCut]);
        fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut] = new TH2F("TruePrimaryClusConvGamma_Pt_MCPt", "ESD_TruePrimaryClusConvGamma_MCPt", nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut]->SetYTitle("p_{T, MC} (GeV/c)");
        fTrueList[iCut]->Add(fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut]);
        fHistoTrueSecondaryClusGammaPt[iCut]          = new TH2F("ESD_TrueSecondaryClusGamma_Pt", "ESD_TrueSecondaryClusGamma_Pt", nBinsClusterPt, arrClusPtBinning, 5, -0.5, 4.5);
        fHistoTrueSecondaryClusGammaPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoTrueSecondaryClusGammaPt[iCut]->SetYTitle("source");
        fHistoTrueSecondaryClusGammaPt[iCut]->GetYaxis()->SetBinLabel( 1,"K0s");
        fHistoTrueSecondaryClusGammaPt[iCut]->GetYaxis()->SetBinLabel( 2,"K0l");
        fHistoTrueSecondaryClusGammaPt[iCut]->GetYaxis()->SetBinLabel( 3,"Lambda");
        fHistoTrueSecondaryClusGammaPt[iCut]->GetYaxis()->SetBinLabel( 4,"Eta");
        fHistoTrueSecondaryClusGammaPt[iCut]->GetYaxis()->SetBinLabel( 5,"rest");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaPt[iCut]);
        fHistoTrueSecondaryClusConvGammaPt[iCut]      = new TH2F("ESD_TrueSecondaryClusConvGamma_Pt", "ESD_TrueSecondaryClusConvGamma_Pt", nBinsClusterPt, arrClusPtBinning, 5, -0.5, 4.5);
        fHistoTrueSecondaryClusConvGammaPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoTrueSecondaryClusConvGammaPt[iCut]->SetYTitle("source");
        fHistoTrueSecondaryClusConvGammaPt[iCut]->GetYaxis()->SetBinLabel( 1,"K0s");
        fHistoTrueSecondaryClusConvGammaPt[iCut]->GetYaxis()->SetBinLabel( 2,"K0l");
        fHistoTrueSecondaryClusConvGammaPt[iCut]->GetYaxis()->SetBinLabel( 3,"Lambda");
        fHistoTrueSecondaryClusConvGammaPt[iCut]->GetYaxis()->SetBinLabel( 4,"Eta");
        fHistoTrueSecondaryClusConvGammaPt[iCut]->GetYaxis()->SetBinLabel( 5,"rest");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaPt[iCut]);
        fHistoTrueSecondaryClusGammaMCPt[iCut]          = new TH2F("ESD_TrueSecondaryClusGamma_MCPt", "ESD_TrueSecondaryClusGamma_MCPt", nBinsClusterPt, arrClusPtBinning, 5, -0.5, 4.5);
        fHistoTrueSecondaryClusGammaMCPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fHistoTrueSecondaryClusGammaMCPt[iCut]->SetYTitle("source");
        fHistoTrueSecondaryClusGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 1,"K0s");
        fHistoTrueSecondaryClusGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 2,"K0l");
        fHistoTrueSecondaryClusGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 3,"Lambda");
        fHistoTrueSecondaryClusGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 4,"Eta");
        fHistoTrueSecondaryClusGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 5,"rest");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaMCPt[iCut]);
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]      = new TH2F("ESD_TrueSecondaryClusConvGamma_MCPt", "ESD_TrueSecondaryClusConvGamma_MCPt", nBinsClusterPt, arrClusPtBinning, 5, -0.5, 4.5);
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->SetYTitle("source");
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 1,"K0s");
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 2,"K0l");
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 3,"Lambda");
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 4,"Eta");
        fHistoTrueSecondaryClusConvGammaMCPt[iCut]->GetYaxis()->SetBinLabel( 5,"rest");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaMCPt[iCut]);

        fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusGammaFromXFromK0s_MCPt_Pt", "ESD_TrueSecondaryClusGammaFromXFromK0s_MCPt_Pt",
                                                                           nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[iCut]);
        fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusConvGammaFromXFromK0s_MCPt_Pt", "ESD_TrueSecondaryClusConvGammaFromXFromK0s_MCPt_Pt",
                                                                               nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[iCut]);
        fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusGammaFromXFromK0l_MCPt_Pt", "ESD_TrueSecondaryClusGammaFromXFromK0l_MCPt_Pt",
                                                                           nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[iCut]);
        fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusConvGammaFromXFromK0l_MCPt_Pt", "ESD_TrueSecondaryClusConvGammaFromXFromK0l_MCPt_Pt",
                                                                               nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[iCut]);
        fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusGammaFromXFromLambda_MCPt_Pt", "ESD_TrueSecondaryClusGammaFromXFromLambda_MCPt_Pt",
                                                                              nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[iCut]);
        fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[iCut] = new TH2F("ESD_TrueSecondaryClusConvGammaFromXFromLambda_MCPt_Pt", "ESD_TrueSecondaryClusConvGammaFromXFromLambda_MCPt_Pt",
                                                                                  nBinsClusterPt, arrClusPtBinning, nBinsClusterPt, arrClusPtBinning);
        fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[iCut]->SetXTitle("p_{T, MC} (GeV/c)");
        fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
        fTrueList[iCut]->Add(fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[iCut]);

        fHistoTrueNLabelsInClus[iCut]                 = new TH1F("TrueNLabelsInClus", "TrueNLabelsInClus", 100, -0.5, 99.5);
        fHistoTrueNLabelsInClus[iCut]->SetXTitle("# labels");
        fTrueList[iCut]->Add(fHistoTrueNLabelsInClus[iCut]);
      }
      fHistoDoubleCountTrueClusterGammaPt[iCut]     = new TH2F("TrueDoubleCountClusterGamma_Pt", "TrueDoubleCountClusterGamma_Pt", nBinsClusterPt, arrClusPtBinning, 2, 0, 2);
      fHistoDoubleCountTrueClusterGammaPt[iCut]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
      fTrueList[iCut]->Add(fHistoDoubleCountTrueClusterGammaPt[iCut]);
      fHistoMultipleCountTrueClusterGamma[iCut]     = new TH1F("TrueMultipleCountClusterGamma", "TrueMultipleCountClusterGamma", 10, 1, 11);
      fHistoMultipleCountTrueClusterGamma[iCut]->SetXTitle("# multiple");
      fTrueList[iCut]->Add(fHistoMultipleCountTrueClusterGamma[iCut]);
      fHistoTrueClusGammaEResE[iCut]          = new TH2F("TrueClusGammaERes_E", "TrueClusGammaERes_E", nBinsClusterPt, arrClusPtBinning,  arrResBinning.size()-1, arrResBinning.data());
      fHistoTrueClusGammaEResE[iCut]->SetXTitle("#it{E}_{rec} (GeV)");
      fHistoTrueClusGammaEResE[iCut]->SetYTitle("(#it{E}_{rec}-#it{E}_{true})/#it{E}_{rec}");
      fTrueList[iCut]->Add(fHistoTrueClusGammaEResE[iCut]);
      fHistoTrueClusPhotonGammaEResE[iCut]          = new TH2F("TrueClusPhotonGammaERes_E", "TrueClusPhotonGammaERes_E", nBinsClusterPt, arrClusPtBinning, arrResBinning.size()-1, arrResBinning.data());
      fHistoTrueClusPhotonGammaEResE[iCut]->SetXTitle("#it{E}_{rec} (GeV)");
      fHistoTrueClusPhotonGammaEResE[iCut]->SetYTitle("(#it{E}_{rec}-#it{E}_{true})/#it{E}_{rec}");
      fTrueList[iCut]->Add(fHistoTrueClusPhotonGammaEResE[iCut]);

      if (fIsMC > 1){
        fHistoTrueClusGammaPt[iCut]->Sumw2();
        fHistoDoubleCountTrueClusterGammaPt[iCut]->Sumw2();
        fHistoMultipleCountTrueClusterGamma[iCut]->Sumw2();
        if(!fDoLightOutput){
          fHistoTrueNLabelsInClus[iCut]->Sumw2();
          fHistoTruePrimaryClusGammaPt[iCut]->Sumw2();
          fHistoTruePrimaryClusGammaESDPtMCPt[iCut]->Sumw2();
          fHistoTruePrimaryClusConvGammaPt[iCut]->Sumw2();
          fHistoTruePrimaryClusConvGammaESDPtMCPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusGammaPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusConvGammaPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusGammaMCPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusConvGammaMCPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[iCut]->Sumw2();
          fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[iCut]->Sumw2();
        }
      }

      if(fDoMesonAnalysis){
        fHistoTruePi0InvMassPt[iCut]                    = new TH2F("ESD_TruePi0_InvMass_Pt", "ESD_TruePi0_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoTruePi0InvMassPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoTruePi0InvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fTrueList[iCut]->Add(fHistoTruePi0InvMassPt[iCut]);
        if (!fDoLightOutput){
            fHistoTruePi0InvMassPtAdditional[iCut]                = new TH2F("ESD_TruePi0_InvMass_Pt_Additional", "ESD_TruePi0_InvMass_Pt_Additional", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
            fHistoTruePi0InvMassPtAdditional[iCut]->SetXTitle("M_{inv,#pi^{0}}(GeV/c^{2})");
            fHistoTruePi0InvMassPtAdditional[iCut]->SetYTitle("#pi^{0}p_{T}(GeV/c)");
            fTrueList[iCut]->Add(fHistoTruePi0InvMassPtAdditional[iCut]);
        }
        fHistoDoubleCountTruePi0InvMassPt[iCut]         = new TH2F("ESD_TrueDoubleCountPi0_InvMass_Pt", "ESD_TrueDoubleCountPi0_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoDoubleCountTruePi0InvMassPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoDoubleCountTruePi0InvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fTrueList[iCut]->Add(fHistoDoubleCountTruePi0InvMassPt[iCut]);
        if( !fDoPi0Only ){
          fHistoTrueEtaInvMassPt[iCut]                    = new TH2F("ESD_TrueEta_InvMass_Pt", "ESD_TrueEta_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTrueEtaInvMassPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
          fHistoTrueEtaInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fTrueList[iCut]->Add(fHistoTrueEtaInvMassPt[iCut]);
          if (!fDoLightOutput){
              fHistoTrueEtaInvMassPtAdditional[iCut]                =  new TH2F("ESD_TrueEta_InvMass_Pt_Additional", "ESD_TrueEta_InvMass_PtAdditional", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
              fHistoTrueEtaInvMassPtAdditional[iCut]->SetXTitle("M_{inv,#eta}(GeV/c^{2})");
              fHistoTrueEtaInvMassPtAdditional[iCut]->SetYTitle("#eta p_{T}(GeV/c)");
              fTrueList[iCut]->Add(fHistoTrueEtaInvMassPtAdditional[iCut]);
          }
          fHistoDoubleCountTrueEtaInvMassPt[iCut]         = new TH2F("ESD_TrueDoubleCountEta_InvMass_Pt", "ESD_TrueDoubleCountEta_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoDoubleCountTrueEtaInvMassPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
          fHistoDoubleCountTrueEtaInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fTrueList[iCut]->Add(fHistoDoubleCountTrueEtaInvMassPt[iCut]);
          fHistoTruePrimaryEtaInvMassPt[iCut]             = new TH2F("ESD_TruePrimaryEta_InvMass_Pt", "ESD_TruePrimaryEta_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTruePrimaryEtaInvMassPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
          fHistoTruePrimaryEtaInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fHistoTruePrimaryEtaInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePrimaryEtaInvMassPt[iCut]);
          fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]  = new TH2F("ESD_TruePrimaryEtaW0Weights_InvMass_Pt", "ESD_TruePrimaryEtaW0Weights_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
          fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]);
          fProfileTruePrimaryEtaWeightsInvMassPt[iCut]    = new TProfile2D("ESD_TruePrimaryEtaWeights_InvMass_Pt", "ESD_TruePrimaryEtaWeights_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fProfileTruePrimaryEtaWeightsInvMassPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
          fProfileTruePrimaryEtaWeightsInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fProfileTruePrimaryEtaWeightsInvMassPt[iCut]->Sumw2();
          fTrueList[iCut]->Add(fProfileTruePrimaryEtaWeightsInvMassPt[iCut]);
        }
        fHistoTruePrimaryPi0InvMassPt[iCut]             = new TH2F("ESD_TruePrimaryPi0_InvMass_Pt", "ESD_TruePrimaryPi0_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoTruePrimaryPi0InvMassPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoTruePrimaryPi0InvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoTruePrimaryPi0InvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTruePrimaryPi0InvMassPt[iCut]);
        fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]  = new TH2F("ESD_TruePrimaryPi0W0Weights_InvMass_Pt", "ESD_TruePrimaryPi0W0Weights_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]);
        fProfileTruePrimaryPi0WeightsInvMassPt[iCut]    = new TProfile2D("ESD_TruePrimaryPi0Weights_InvMass_Pt", "ESD_TruePrimaryPi0Weights_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fProfileTruePrimaryPi0WeightsInvMassPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
        fProfileTruePrimaryPi0WeightsInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fProfileTruePrimaryPi0WeightsInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fProfileTruePrimaryPi0WeightsInvMassPt[iCut]);
        fHistoTrueSecondaryPi0InvMassPt[iCut]           = new TH2F("ESD_TrueSecondaryPi0_InvMass_Pt", "ESD_TrueSecondaryPi0_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryPi0InvMassPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoTrueSecondaryPi0InvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoTrueSecondaryPi0InvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0InvMassPt[iCut]);

        if(!fDoLightOutput){
          fHistoTruePi0InvMassPtAlpha[iCut]               = new TH2F("ESD_TruePi0_InvMass_vs_Pt_Alpha", "ESD_TruePi0_InvMass_vs_Pt_Alpha", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTruePi0InvMassPtAlpha[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
          fHistoTruePi0InvMassPtAlpha[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fHistoTruePi0InvMassPtAlpha[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoTruePi0InvMassPtAlpha[iCut]);
          fHistoTruePi0PureGammaInvMassPtAlpha[iCut]      = new TH2F("ESD_TruePi0PureGamma_InvMass_vs_Pt_Alpha", "ESD_TruePi0PureGamma_InvMass_vs_Pt_Alpha", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
          fHistoTruePi0PureGammaInvMassPtAlpha[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
          fHistoTruePi0PureGammaInvMassPtAlpha[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
          fHistoTruePi0PureGammaInvMassPtAlpha[iCut]->Sumw2();
          fESDList[iCut]->Add(fHistoTruePi0PureGammaInvMassPtAlpha[iCut]);
        }
        fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]    = new TH2F("ESD_TrueSecondaryPi0FromK0s_InvMass_Pt", "ESD_TrueSecondaryPi0FromK0s_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromK0sInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]    = new TH2F("ESD_TrueSecondaryPi0FromK0l_InvMass_Pt", "ESD_TrueSecondaryPi0FromK0l_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]->Sumw2();
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromK0lInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]    = new TH2F("ESD_TrueSecondaryPi0FromEta_InvMass_Pt", "ESD_TrueSecondaryPi0FromEta_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]);
        fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryPi0FromLambda_InvMass_Pt", "ESD_TrueSecondaryPi0FromLambda_InvMass_Pt", nBinsMinv, 0, maxMinv, nBinsPt, arrPtBinning);
        fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]->SetYTitle("#it{p}_{T} (GeV/#it{c})");
        fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
        fTrueList[iCut]->Add(fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]);
        if (fIsMC > 1){
          fHistoTruePi0InvMassPt[iCut]->Sumw2();
          fHistoDoubleCountTruePi0InvMassPt[iCut]->Sumw2();
          fHistoTruePrimaryPi0W0WeightingInvMassPt[iCut]->Sumw2();
          fHistoTrueSecondaryPi0FromEtaInvMassPt[iCut]->Sumw2();
          fHistoTrueSecondaryPi0FromLambdaInvMassPt[iCut]->Sumw2();
          if (!fDoLightOutput){
              fHistoTruePi0InvMassPtAdditional[iCut]->Sumw2();
          }
          if( !fDoPi0Only ){
            fHistoTrueEtaInvMassPt[iCut]->Sumw2();
            if (!fDoLightOutput){
                fHistoTrueEtaInvMassPtAdditional[iCut]->Sumw2();
            }
            fHistoDoubleCountTrueEtaInvMassPt[iCut]->Sumw2();
            fHistoTruePrimaryEtaW0WeightingInvMassPt[iCut]->Sumw2();
          }
        }
      }
    }
  }
  fVectorDoubleCountTruePi0s.clear();
  fVectorDoubleCountTrueEtas.clear();
  fVectorDoubleCountTrueClusterGammas.clear();

  fMapMultipleCountTrueClusterGammas.clear();

  if(fV0Reader)
    if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
      if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
  if(fV0Reader)
    if((AliConvEventCuts*)fV0Reader->GetEventCuts())
      if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
        fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());

  for(int iMatcherTask = 0; iMatcherTask < 5; iMatcherTask++){
    AliCaloTrackMatcher* temp = 0x0;
    if(!fCorrTaskSetting.CompareTo("")){
      temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcher_%i_%i",iMatcherTask,fTrackMatcherRunningMode)));
    } else {
      temp = (AliCaloTrackMatcher*) (AliAnalysisManager::GetAnalysisManager()->GetTask(Form("CaloTrackMatcher_%i_%i_%s",iMatcherTask,fTrackMatcherRunningMode,fCorrTaskSetting.Data())));
    }
    if(temp) fOutputContainer->Add(temp->GetCaloTrackMatcherHistograms());
  }

  for(int iCut = 0; iCut<fnCuts;iCut++){
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

  if (fIsMC > 0 ){
    tBrokenFiles = new TTree("BrokenFiles", "BrokenFiles");
    tBrokenFiles->Branch("fileName",&fFileNameBroken);
    fOutputContainer->Add(tBrokenFiles);
  }

  if(fLocalDebugFlag > 0){
    fstream fOutputLocalDebug;
    fOutputLocalDebug.open("debugOutput.txt",ios::out);
    fOutputLocalDebug.close();
  }

  ////////////////////////
  // VZERO
  ////////////////////////
  LoadVZEROCalibration();

  OpenFile(1);
  PostData(1, fOutputContainer);
}
//_____________________________________________________________________________
bool AliAnalysisTaskPi0v2Calo::Notify()
{
  for(int iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){
        ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnum(fV0Reader->GetPeriodName());
    }
    if(fIsHeavyIon==2) {
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
  }

  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskPi0v2Calo::UserExec(Option_t *)
{
  //
  // Called for each event
  //
  fInputEvent           = InputEvent();
  if(fIsMC> 0) fMCEvent = MCEvent();
  int eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();

  if(fInputEvent->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
  if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1 or because it is incomplete
    // write out name of broken file for first event
    if (fIsMC > 0){
      if (fInputEvent->IsA()==AliESDEvent::Class()){
        if (((AliESDEvent*)fInputEvent)->GetEventNumberInFile() == 0){
          fFileNameBroken = new TObjString(Form("%s", ((TString)fV0Reader->GetCurrentFileName()).Data()));
          if (tBrokenFiles) tBrokenFiles->Fill();
          delete fFileNameBroken;
        }
      }
    }

    for(int iCut = 0; iCut<fnCuts; iCut++){
      fHistoNEvents[iCut]->Fill(eventQuality);
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality);
    }
    return;
  }
  if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoSecondaryTrackMatching())fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut

  // ------------------- BeginEvent ----------------------------

  AliEventplane *EventPlane = fInputEvent->GetEventplane();
  if(fIsHeavyIon ==1)fEventPlaneAngle = EventPlane->GetEventplane("V0",fInputEvent,2);
  else fEventPlaneAngle=0.0;
  for(int iCut = 0; iCut<fnCuts; iCut++){

    fiCut = iCut;
    fNCurrentClusterBasic       = 0;
    bool isRunningEMCALrelAna = kFALSE;
    if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType() == 1) isRunningEMCALrelAna = kTRUE;
    fCent = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetCentrality(fInputEvent);

    // VZERO Plane
    fRunNumber = fInputEvent->GetRunNumber();
    if (IsVZEROCalibOn)
    {
      if (fRunNumber != fOldRunNumber) {
        if (!LoadCalibHistForThisRun()) {
          AliInfo("=====faild to LoadCalibHistForThisRun=====");
          return;
        }
        fOldRunNumber = fRunNumber;
      }
      if (!GetVZEROPlane()) {
        AliInfo("==============faild to GetVZEROPlane===========");
        return;
      }
    }
    int centBin = 999; //  depend on GammaCalo trainconfig
    if (fCent >= 0 && fCent < 10) {
      centBin = 0;
    } else if (fCent >= 10 && fCent < 30) {
      centBin = 1;
    } else if (fCent >= 30 && fCent < 50) {
      centBin = 2;
    } else if (fCent >= 50 && fCent < 90) {
      centBin = 3;
    }
    fHist2V0Res[iCut]->Fill(centBin + 0.5, cos(2 * (fPsi2V0A - fPsi2V0C)));
    fHist2DPsi2V0CCent[iCut]->Fill(fCent, fPsi2V0C);
    fHist2DPsi2V0ACent[iCut]->Fill(fCent, fPsi2V0A);

    int eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon, isRunningEMCALrelAna);

    if(fIsMC==2){
      float xsection      = -1.;
      float ntrials       = -1.;
      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetXSectionAndNTrials(fMCEvent,xsection,ntrials, fInputEvent );
      if((xsection==-1.) || (ntrials==-1.)) AliFatal(Form("ERROR: GetXSectionAndNTrials returned invalid xsection/ntrials, periodName from V0Reader: '%s'",fV0Reader->GetPeriodName().Data()));
      fProfileJetJetXSection[iCut]->Fill(0.,xsection);
      fHistoJetJetNTrials[iCut]->Fill("#sum{NTrials}", ntrials);
    }

    if (fIsMC > 0){
      fWeightJetJetMC       = 1;
      float maxjetpt      = -1.;
      float pthard = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetUseJetFinderForOutliers()) maxjetpt = fOutlierJetReader->GetMaxJetPt();
      bool isMCJet        = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsJetJetMCEventAccepted( fMCEvent, fWeightJetJetMC ,pthard, fInputEvent, maxjetpt);
      if(isMCJet && (fIsMC==2))           fHistoPtHardJJWeight[iCut]->Fill(pthard,fWeightJetJetMC);
      if (fIsMC == 3){
        double weightMult   = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetWeightForMultiplicity(fV0Reader->GetNumberOfPrimaryTracks());
        fWeightJetJetMC       = fWeightJetJetMC*weightMult;
      }

      if (!isMCJet){
        fHistoNEvents[iCut]->Fill(10,fWeightJetJetMC);
        if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(10);
        continue;
      }
    }
    if(eventNotAccepted!= 0){
      fHistoNEvents[iCut]->Fill(eventNotAccepted, fWeightJetJetMC); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventNotAccepted);

      if(fIsMC > 0 && (eventQuality == 0 || eventQuality == 3 || eventQuality == 5)){
        if (eventNotAccepted==3){ // wrong trigger selected. However, we still want to count the MC particles fot these events! If MC particles should be rejected in addition, use IsMCTriggerSelected function
          if(fInputEvent->IsA()==AliESDEvent::Class())
            ProcessMCParticles(1);
          if(fInputEvent->IsA()==AliAODEvent::Class())
            ProcessAODMCParticles(1);
        }
      }
      continue;
    }
    if(eventQuality != 0){// Event Not Accepted
      fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC);
      if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here

      if(fIsMC > 0){
        if(eventQuality == 3 || eventQuality == 5){ // 3 = trigger, 5 = contr. to vertex
          if(fInputEvent->IsA()==AliESDEvent::Class())
            ProcessMCParticles(2);
          if(fInputEvent->IsA()==AliAODEvent::Class())
            ProcessAODMCParticles(2);
        }
      }
      continue;
    }

    fHistoNEvents[iCut]->Fill(eventQuality, fWeightJetJetMC); // Should be 0 here
    if (fIsMC>1) fHistoNEventsWOWeight[iCut]->Fill(eventQuality); // Should be 0 here
    fHistoNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(), fWeightJetJetMC);
    fHistoVertexZ[iCut]->Fill(fInputEvent->GetPrimaryVertex()->GetZ(), fWeightJetJetMC);
    if(!fDoLightOutput){
      fHistoSPDClusterTrackletBackground[iCut]->Fill(fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),(fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)), fWeightJetJetMC);
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsHeavyIon() == 2)  fHistoNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A(), fWeightJetJetMC);
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
        if(!fDoLightOutput){
          if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader()){
            for(int i = 0;i<(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader())->GetEntries();i++){
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
    }

    if(fIsMC> 0){
    if(fInputEvent->IsA()==AliESDEvent::Class())
      ProcessMCParticles(0);
    if(fInputEvent->IsA()==AliAODEvent::Class())
      ProcessAODMCParticles(0);
    }

    // it is in the loop to have the same conversion cut string (used also for MC stuff that should be same for V0 and Cluster)
    ProcessClusters();            // process calo clusters
    if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoSecondaryTrackMatching()) ProcessConversionCandidates(); // process conversion candidates for secondary track matching

    fHistoNGammaCandidatesBasic[iCut]->Fill(fNCurrentClusterBasic, fWeightJetJetMC);
    fHistoNGammaCandidates[iCut]->Fill(fClusterCandidates->GetEntries(), fWeightJetJetMC);
    if(!fDoLightOutput) fHistoNGoodESDTracksVsNGammaCandidates[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),fClusterCandidates->GetEntries(), fWeightJetJetMC);
    if(fDoMesonAnalysis){ // Meson Analysis
      CalculatePi0Candidates(); // Combine Gammas from conversion and from calo
      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
          CalculateBackground(); // Combinatorial Background
          UpdateEventByEventData(); // Store Event for mixed Events
        } else if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 2){
          CalculateBackgroundSwapp(); // Combinatorial Background
        }

      }
      fVectorDoubleCountTruePi0s.clear();
      fVectorDoubleCountTrueEtas.clear();
    }
    if(fIsMC> 0){
      fVectorDoubleCountTrueClusterGammas.clear();
      FillMultipleCountHistoAndClear(fMapMultipleCountTrueClusterGammas,fHistoMultipleCountTrueClusterGamma[iCut]);
    }
    if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoSecondaryTrackMatching()) fGammaCandidates->Clear();
    fClusterCandidates->Clear(); // delete cluster candidates
  }
  PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2Calo::ProcessClusters()
{

  int nclus                       = 0;
  TClonesArray * arrClustersProcess = NULL;
  fNCurrentClusterBasic             = 0;
  if(!fCorrTaskSetting.CompareTo("")){
    nclus = fInputEvent->GetNumberOfCaloClusters();
  } else {
    arrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!arrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskPi0v2Calo! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclus = arrClustersProcess->GetEntries();
  }

  // energy correction for neutral overlap!
  if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoEnergyCorrectionForOverlap() > 0 && ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoEnergyCorrectionForOverlap() <= 2){
    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->SetPoissonParamCentFunction(fIsMC);
    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->SetNMatchedTracksFunc(fCent);
  }

  vector<AliAODConversionPhoton*>         vectorCurrentClusters;
  vector<int>                           vectorRejectCluster;
  vector<double>                        vectorPhotonWeight;
  vector<double>                        vectorClusterM02;
  vector<bool>                          vectorIsFromDesiredHeader;
  vector<int>                           vectorCurrentClusters_DDL;

  if(nclus == 0)  return;
  // plotting histograms on cell/tower level, only if extendedMatchAndQA > 1
  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillHistogramsExtendedQA(fInputEvent,fIsMC);

  // match tracks to clusters
  if(fDoPrimaryTrackMatching) ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kTRUE, fMCEvent);

  // vertex
  double vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  int totalActiveCells = 0;
  double totalCellsinClusters = 0;
  double totalUnclusteredE = 0;
  if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoFlatEnergySubtraction()){
    double totalClusterEnergy = 0;
    totalActiveCells = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetNactiveEmcalCells();
    totalCellsinClusters = 0;
    for(long i = 0; i < nclus; i++){
      AliVCluster* clus = NULL;
      if(fInputEvent->IsA()==AliESDEvent::Class()){
        if(arrClustersProcess)
          clus = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(i));
        else
          clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
      } else if(fInputEvent->IsA()==AliAODEvent::Class()){
        if(arrClustersProcess)
          clus = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(i));
        else
          clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
      }
      if(!clus) continue;
      totalClusterEnergy += clus->E();
      totalCellsinClusters += clus->GetNCells();
      delete clus;
    }
    double totalEDeposit = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetTotalEnergyDeposit(fInputEvent);
    totalUnclusteredE = totalEDeposit - totalClusterEnergy;
  }

  map<long,int> mapIsClusterAccepted;
  map<long,int> mapIsClusterAcceptedWithoutTrackMatch;
  // Loop over EMCal clusters
  for(long i = 0; i < nclus; i++){
    double tempClusterWeight        = fWeightJetJetMC;
    double tempPhotonWeight         = fWeightJetJetMC;
    AliVCluster* clus = NULL;
    if(fInputEvent->IsA()==AliESDEvent::Class()){
      if(arrClustersProcess)
        clus = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(i));
      else
        clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
    } else if(fInputEvent->IsA()==AliAODEvent::Class()){
      if(arrClustersProcess)
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(i));
      else
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
    }
    if(!clus) continue;

    // energy correction for neutral overlap!
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoFlatEnergySubtraction()){
      if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoEnergyCorrectionForOverlap() > 0){
        // Neutral Overlap correction via mean number of matched primary tracks per cluster
        if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoEnergyCorrectionForOverlap() <= 2){
          clus->SetE(clus->E() - ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CorrectEnergyForOverlap(fCent));
        } else if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoEnergyCorrectionForOverlap() == 3){
          // Neutral Overlap correction via mean number of charged particles per cell
          clus->SetE(clus->E() - 12.f * ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CorrectEnergyForOverlap(fCent));
        } else if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoEnergyCorrectionForOverlap() >= 4){
          // Neutral Overlap correction via NonLin like correction
          clus->SetE(clus->E() * ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CorrectEnergyForOverlap(fCent, clus->E()));
        }
      }
    }

    if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoFlatEnergySubtraction()){
      if(totalUnclusteredE!=0 && totalActiveCells!=0 && totalCellsinClusters!=0){
        clus->SetE(clus->E()-(clus->GetNCells()*(totalUnclusteredE)/(totalActiveCells-totalCellsinClusters)));
      }
    }

    // Set the jetjet weight to 1 in case the cluster orignated from the minimum bias header
    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
      if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(clus->GetLabelAt(0), fMCEvent, fInputEvent) == 2)
        tempClusterWeight = 1;
    }

    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
      if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetIsAcceptedForBasicCounting())fNCurrentClusterBasic++;
      delete clus;
      continue;
    }

    fNCurrentClusterBasic++;

    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    clus->GetMomentum(clusterVector,vertex);

    TLorentzVector* tmpvec = new TLorentzVector();
    tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());

    // convert to AODConversionPhoton
    AliAODConversionPhoton *PhotonCandidate=new AliAODConversionPhoton(tmpvec);
    if(!PhotonCandidate){ delete clus; delete tmpvec; continue;}

    // Flag Photon as CaloPhoton
    PhotonCandidate->SetIsCaloPhoton(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetClusterType());
    PhotonCandidate->SetCaloClusterRef(i);
    PhotonCandidate->SetLeadingCellID(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FindLargestCellInCluster(clus,fInputEvent));
    // get MC label
    if(fIsMC> 0){
      int* mclabelsCluster = clus->GetLabels();
      PhotonCandidate->SetNCaloPhotonMCLabels(clus->GetNLabels());
     // cout << clus->GetNLabels() << endl;
      if (clus->GetNLabels()>0){
        for (int k =0; k< (int)clus->GetNLabels(); k++){
          PhotonCandidate->SetCaloPhotonMCLabel(k,mclabelsCluster[k]);
          // int pdgCode = fMCEvent->GetTrack(mclabelsCluster[k])->PdgCode();
          // cout << "label " << k << "\t" << mclabelsCluster[k] << " pdg code: " << pdgCode << endl;
        } // end of label loop
      }
    }
    fIsFromDesiredHeader          = kTRUE;
    fIsOverlappingWithOtherHeader = kFALSE;
    fIsOverlapWithMBHeader        = kFALSE;
    int particleFromBGEvent = 0;
    int labelFromBGEvent = 0;
    // test whether largest contribution to cluster orginates in added signals
    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() > 0){
      // Set the jetjet weight to 1 in case the photon candidate orignated from the minimum bias header
      particleFromBGEvent = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent);
      if ( particleFromBGEvent == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempPhotonWeight = 1;
      if ( particleFromBGEvent == 0) fIsFromDesiredHeader = kFALSE;
      if (clus->GetNLabels()>1){
        int* mclabelsCluster = clus->GetLabels();
        if (fLocalDebugFlag > 1)   cout << "testing if other labels in cluster belong to different header, need to test " << (int)clus->GetNLabels()-1 << " additional labels" << endl;
        for (int l = 1; l < (int)clus->GetNLabels(); l++ ){
          labelFromBGEvent = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(mclabelsCluster[l], fMCEvent, fInputEvent, fLocalDebugFlag);
          if (labelFromBGEvent == 0) fIsOverlappingWithOtherHeader = kTRUE;
          if (labelFromBGEvent == 2) fIsOverlapWithMBHeader = kTRUE;
        }
        if (fLocalDebugFlag > 1 && fIsOverlappingWithOtherHeader) cout << "found overlapping header: " << endl;
      }
    }

    fHistoClusAllHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), tempPhotonWeight);
    if (!fIsFromDesiredHeader) fHistoClusRejectedHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), tempPhotonWeight);
    if (fIsFromDesiredHeader && fIsOverlapWithMBHeader) fHistoClusOverlapMBHeaderGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), tempPhotonWeight);
    if (fIsFromDesiredHeader && fIsOverlappingWithOtherHeader) fHistoClusOverlapHeadersGammaPt[fiCut]->Fill(PhotonCandidate->Pt(), tempPhotonWeight);


    if ( (fIsFromDesiredHeader && !fIsOverlappingWithOtherHeader && !fAllowOverlapHeaders) || (fIsFromDesiredHeader && fAllowOverlapHeaders) ){
      vectorCurrentClusters.push_back(PhotonCandidate);
      vectorPhotonWeight.push_back(tempPhotonWeight);
      vectorClusterM02.push_back(clus->GetM02());
      vectorIsFromDesiredHeader.push_back(fIsFromDesiredHeader);
    } else{
      delete PhotonCandidate;
    }
    delete clus;
    delete tmpvec;
  }
  bool rejected = kFALSE;
  // run conversion recovery in addition
  if (((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetIsConversionRecovery()){
    rejected = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckForReconstructedConversionPairs(vectorCurrentClusters,vectorRejectCluster);
    if (fLocalDebugFlag > 1 && rejected) cout << "found clusters to be rejected" << endl;
  }

  for (int iter = 0; iter < (int)vectorCurrentClusters.size();iter++){

    if (!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckVectorForIndexAndAdd(vectorRejectCluster, iter,kFALSE)){
      fIsFromDesiredHeader = vectorIsFromDesiredHeader.at(iter);
      fHistoClusGammaPt[fiCut]->Fill(vectorCurrentClusters.at(iter)->Pt(), vectorPhotonWeight.at(iter));
      fHistoClusGammaE[fiCut]->Fill(vectorCurrentClusters.at(iter)->E(), vectorPhotonWeight.at(iter));
      if(fIsMC> 0){
        if(fInputEvent->IsA()==AliESDEvent::Class()){
          ProcessTrueClusterCandidates(vectorCurrentClusters.at(iter),vectorClusterM02.at(iter));
        } else {
          ProcessTrueClusterCandidatesAOD(vectorCurrentClusters.at(iter),vectorClusterM02.at(iter));
        }
      }
      fClusterCandidates->Add(vectorCurrentClusters.at(iter));
    } else {
      if (fLocalDebugFlag > 1) cout << "removed cluster: " << iter << endl;
    }
  }
  if (fLocalDebugFlag > 1 && rejected) cout << "reduced from: " << vectorCurrentClusters.size() << " to " << fClusterCandidates->GetEntries() << endl;
  vectorRejectCluster.clear();
  vectorPhotonWeight.clear();
  vectorClusterM02.clear();
  vectorIsFromDesiredHeader.clear();

  if(fLocalDebugFlag == 2) EventDebugMethod();

  return;
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2Calo::ProcessConversionCandidates(){

  // Loop over Photon Candidates allocated by ReaderV1
  for(int i = 0; i < fReaderGammas->GetEntriesFast(); i++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
    if(!PhotonCandidate) continue;
    if(!fConversionCuts->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
    if(!fConversionCuts->InPlaneOutOfPlaneCut(PhotonCandidate->GetPhotonPhi(),fEventPlaneAngle)) continue;
    if(!fConversionCuts->UseElecSharingCut() && !fConversionCuts->UseToCloseV0sCut()){
      fGammaCandidates->Add(PhotonCandidate);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2Calo::ProcessTrueClusterCandidates(AliAODConversionPhoton *TruePhotonCandidate, double clusterM02)
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  double mcProdVtxX   = primVtxMC->GetX();
  double mcProdVtxY   = primVtxMC->GetY();
  double mcProdVtxZ   = primVtxMC->GetZ();

  double tempPhotonWeight       = fWeightJetJetMC;
  AliMCParticle *Photon = NULL;
  if (TruePhotonCandidate->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");

  if (TruePhotonCandidate->GetCaloPhotonMCLabel(0) < 0) return;
  if (TruePhotonCandidate->GetNCaloPhotonMCLabels() > 0) Photon = (AliMCParticle*) fMCEvent->GetTrack(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
    else return;

  if(Photon == NULL){
  //    cout << "no photon" << endl;
    return;
  }
  // Set the jetjet weight to 1 in case the photon orignated from the minimum bias header
  if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
    if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TruePhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2) tempPhotonWeight = 1;
  }

  TruePhotonCandidate->SetCaloPhotonMCFlags(fMCEvent, fEnableSortForClusMC);

  // True Photon

  if (TruePhotonCandidate->IsLargestComponentPhoton() || (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())) {
    fHistoTrueClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    fHistoTrueClusGammaEResE[fiCut]->Fill(TruePhotonCandidate->E(), (TruePhotonCandidate->E()-Photon->E())/TruePhotonCandidate->E(), tempPhotonWeight);
  }
  if (TruePhotonCandidate->IsLargestComponentPhoton()) fHistoTrueClusPhotonGammaEResE[fiCut]->Fill(TruePhotonCandidate->E(), (TruePhotonCandidate->E()-Photon->E())/TruePhotonCandidate->E(), tempPhotonWeight);
  if(!fDoLightOutput){
    bool isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, TruePhotonCandidate->GetCaloPhotonMCLabel(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
      isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, Photon->GetMother(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);

    if(isPrimary){
      // filling primary histograms
      if (TruePhotonCandidate->IsLargestComponentPhoton()){
        fHistoTruePrimaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
        fHistoTruePrimaryClusGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(), tempPhotonWeight); // Allways Filled
      }
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        fHistoTruePrimaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
        fHistoTruePrimaryClusConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(), tempPhotonWeight); // Allways Filled
      }

    } else {
      // filling secondary histograms
      int secondaryClass    = -1;
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
        secondaryClass        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->SecondaryClassificationPhoton( Photon, fMCEvent, kTRUE);
      else
        secondaryClass        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->SecondaryClassificationPhoton( Photon, fMCEvent, kFALSE);

      // all secondaries
      if (TruePhotonCandidate->IsLargestComponentPhoton()){
        if (secondaryClass == 2) {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 5) {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 3) {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 4) {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,tempPhotonWeight);
        } else {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),4.,tempPhotonWeight);
        }
      }
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        if (secondaryClass == 2) {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 5) {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 3) {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 4) {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,tempPhotonWeight);
        } else {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),4.,tempPhotonWeight);
        }
      }
    }
  }
  return;
}


//________________________________________________________________________
void AliAnalysisTaskPi0v2Calo::ProcessTrueClusterCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate, double clusterM02)
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  double mcProdVtxX   = primVtxMC->GetX();
  double mcProdVtxY   = primVtxMC->GetY();
  double mcProdVtxZ   = primVtxMC->GetZ();

  double tempPhotonWeight       = fWeightJetJetMC;
  AliAODMCParticle *Photon = NULL;
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray){
    if (TruePhotonCandidate->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set task will abort");
    if (TruePhotonCandidate->GetNCaloPhotonMCLabels()>0) Photon = (AliAODMCParticle*) fAODMCTrackArray->At(TruePhotonCandidate->GetCaloPhotonMCLabel(0));
      else return;
  } else {
    AliInfo("AODMCTrackArray could not be loaded");
    return;
  }

  if(Photon == NULL){
  //  cout << "no photon" << endl;
    return;
  }
  TruePhotonCandidate->SetCaloPhotonMCFlagsAOD(fAODMCTrackArray, fEnableSortForClusMC);

  // Set the jetjet weight to 1 in case the cluster orignated from the minimum bias header
  if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
    if ( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TruePhotonCandidate->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2) tempPhotonWeight = 1;
  }

  // True Photon
  if (TruePhotonCandidate->IsLargestComponentPhoton() || (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())) {
    fHistoTrueClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
    fHistoTrueClusGammaEResE[fiCut]->Fill(TruePhotonCandidate->E(), (TruePhotonCandidate->E()-Photon->E())/TruePhotonCandidate->E(), tempPhotonWeight);
  }
  if (TruePhotonCandidate->IsLargestComponentPhoton()) fHistoTrueClusPhotonGammaEResE[fiCut]->Fill(TruePhotonCandidate->E(), (TruePhotonCandidate->E()-Photon->E())/TruePhotonCandidate->E(), tempPhotonWeight);
  if(!fDoLightOutput){
    bool isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
      if (Photon->GetMother()> -1){
        AliAODMCParticle *Mother  = (AliAODMCParticle*) fAODMCTrackArray->At(Photon->GetMother());
        isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, Mother, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
      }
    }
    if(isPrimary){
      if (TruePhotonCandidate->IsLargestComponentPhoton()){
        fHistoTruePrimaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
        fHistoTruePrimaryClusGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(), tempPhotonWeight); // Allways Filled
      }
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        fHistoTruePrimaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(), tempPhotonWeight);
        fHistoTruePrimaryClusConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt(), tempPhotonWeight); // Allways Filled
      }

    } else {
      // filling secondary histograms
      int secondaryClass    = -1;
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion())
        secondaryClass        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->SecondaryClassificationPhotonAOD( Photon, fAODMCTrackArray, kTRUE);
      else
        secondaryClass        = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->SecondaryClassificationPhotonAOD( Photon, fAODMCTrackArray, kFALSE);

      // all secondaries
      if (TruePhotonCandidate->IsLargestComponentPhoton()){
        if (secondaryClass == 2) {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 5) {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 3) {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 4) {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,tempPhotonWeight);
        } else {
          fHistoTrueSecondaryClusGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4.,tempPhotonWeight);
          fHistoTrueSecondaryClusGammaMCPt[fiCut]->Fill(Photon->Pt(),4.,tempPhotonWeight);
        }
      }
      if (TruePhotonCandidate->IsLargestComponentElectron() && TruePhotonCandidate->IsConversion()){
        if (secondaryClass == 2) {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),0.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 5) {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),1.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 3) {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),2.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt[fiCut]->Fill(Photon->Pt(),TruePhotonCandidate->Pt(),tempPhotonWeight);
        } else if (secondaryClass == 4) {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),3.,tempPhotonWeight);
        } else {
          fHistoTrueSecondaryClusConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4.,tempPhotonWeight);
          fHistoTrueSecondaryClusConvGammaMCPt[fiCut]->Fill(Photon->Pt(),4.,tempPhotonWeight);
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2Calo::ProcessAODMCParticles(int isCurrentEventSelected)
{
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  double mcProdVtxX   = primVtxMC->GetX();
  double mcProdVtxY   = primVtxMC->GetY();
  double mcProdVtxZ   = primVtxMC->GetZ();


  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray == NULL) return;

  // Check if MC generated particles should be filled for this event using the selected trigger
  if( !((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsMCTriggerSelected(fInputEvent, fMCEvent)){
    if(isCurrentEventSelected == 0) fHistoMCEventsTrigg[fiCut]->Fill(2., fWeightJetJetMC); // event is rejected and triggered
    else fHistoMCEventsTrigg[fiCut]->Fill(3., fWeightJetJetMC); // event is rejected but not triggered
    return;
  }
  if(isCurrentEventSelected == 0) fHistoMCEventsTrigg[fiCut]->Fill(0., fWeightJetJetMC); // event is accepted and triggered
  else fHistoMCEventsTrigg[fiCut]->Fill(1., fWeightJetJetMC); // event is accepted but not triggered

  // Loop over all primary MC particle
  for(long i = 0; i < fAODMCTrackArray->GetEntriesFast(); i++) {
    double tempParticleWeight       = fWeightJetJetMC;

    AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(i));
    if (!particle) continue;

    bool isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, particle, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if (isPrimary){

      int isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
        // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
        if(isMCFromMBHeader == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempParticleWeight = 1;
      }

      if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(particle,fAODMCTrackArray)){
        if (isCurrentEventSelected==1) {
          fHistoMCGammaPtNotTriggered[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Gamma
        } else if (isCurrentEventSelected==2) {
          fHistoMCGammaPtNoVertex[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Gamma
        }
        fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Gamma

        if(!fDoLightOutput){
          if(particle->GetMother() >-1){ // Meson Decay Gamma
            switch((static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetMother())))->GetPdgCode()){
            case 111: // Pi0
              fHistoMCDecayGammaPi0Pt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 113: // Rho0
              fHistoMCDecayGammaRhoPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 221: // Eta
              fHistoMCDecayGammaEtaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 223: // Omega
              fHistoMCDecayGammaOmegaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 331: // Eta'
              fHistoMCDecayGammaEtapPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 333: // Phi
              fHistoMCDecayGammaPhiPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 3212: // Sigma
              fHistoMCDecayGammaSigmaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            }
          }
        }
      }
      // Converted MC Gamma
      if(fDoMesonAnalysis){
        double mesonY = 1.e30;
        double ratio  = 0;
        if (particle->E() != TMath::Abs(particle->Pz())){
          ratio         = (particle->E()+particle->Pz()) / (particle->E()-particle->Pz());
        }
        if( !(ratio <= 0) ){
          mesonY = particle->Y()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
        }
        if(fHistoMCPrimaryPtvsSource[fiCut]){
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
            }
          }
        }

        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
          ->MesonIsSelectedAODMC(particle,fAODMCTrackArray,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()) && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonLeadTrackSelectionAODMC(fInputEvent, particle)){
          AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(0)));
          AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(1)));
          float weighted= 1;
          if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
            if (particle->Pt()>0.005){
              weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, 0x0, fInputEvent);
            }
          }

          if(particle->GetPdgCode() == 111){
            if (isCurrentEventSelected == 1) fHistoMCPi0PtNotTriggered[fiCut]->Fill(particle->Pt(),weighted*tempParticleWeight); // All MC Pi0 in not triggered collisions
            else if (isCurrentEventSelected == 2) fHistoMCPi0PtNoVertex[fiCut]->Fill(particle->Pt(),weighted*tempParticleWeight); // All MC Pi0 in not triggered collisions
            fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted*tempParticleWeight); // All MC Pi0
          
            fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
            if (fIsMC > 1) fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
          } else if(particle->GetPdgCode() == 221 && !fDoPi0Only){
            if (isCurrentEventSelected == 1) fHistoMCEtaPtNotTriggered[fiCut]->Fill(particle->Pt(),weighted*tempParticleWeight); // All MC Eta in not triggered collisions
            else if (isCurrentEventSelected == 2) fHistoMCEtaPtNoVertex[fiCut]->Fill(particle->Pt(),weighted*tempParticleWeight); // All MC Eta in not triggered collisions
            fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted*tempParticleWeight); // All MC Eta
            fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
            if (fIsMC > 1) fHistoMCEtaWOEvtWeightPt[fiCut]->Fill(particle->Pt());
          }

          // Check the acceptance for both gammas
          if( ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter0,fAODMCTrackArray) &&
              ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter1,fAODMCTrackArray) ){

            if(particle->GetPdgCode() == 111){
              fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc
              if (isCurrentEventSelected==1) fHistoMCPi0InAccPtNotTriggered[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc for events which are not triggered
              if(fIsMC > 1) fHistoMCPi0WOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0 with gamma in acc
            } else if(particle->GetPdgCode() == 221 && !fDoPi0Only){
              fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Eta with gamma in acc
              if(isCurrentEventSelected == 1) fHistoMCEtaInAccPtNotTriggered[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Eta with gamma in acc for not triggered events
              if(fIsMC > 1) fHistoMCEtaWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Eta with gamma in acc
            }
          }
        }
      }
    // fill secondaries
    } else {
      int isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
        if(isMCFromMBHeader == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempParticleWeight = 1;
      }

      if(!fDoLightOutput) {
        if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(particle,fAODMCTrackArray)){
          if(particle->GetMother() >-1){
            AliAODMCParticle *tmpMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetMother()));
            if(tmpMother->GetMother() >-1){
              AliAODMCParticle *tmpGrandMother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(tmpMother->GetMother()));
              if(tmpGrandMother->GetPdgCode() == 310) {
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),0.,tempParticleWeight);
              } else if (tmpGrandMother->GetPdgCode() == 130) {
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),1.,tempParticleWeight);
              } else if (tmpGrandMother->GetPdgCode() == 3122) {
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),2.,tempParticleWeight);
              } else if (tmpGrandMother->GetPdgCode() == 221) {
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,tempParticleWeight);
              } else {
                if( !(TMath::Abs(tmpMother->GetPdgCode()) == 11 && tmpGrandMother->GetPdgCode() == 22) )
                  fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),4.,tempParticleWeight);
              }
            } else {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),4.,tempParticleWeight);
            }
          } else {
            fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),4.,tempParticleWeight);
          }
        }
      }

      if(fDoMesonAnalysis){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedAODMC(particle,fAODMCTrackArray,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()) && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonLeadTrackSelectionAODMC(fInputEvent, particle)){
          AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(0)));
          AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetDaughterLabel(1)));
          AliAODMCParticle* mother = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particle->GetMother()));
          int pdgCode = mother->GetPdgCode();
          if(particle->GetPdgCode() == 111){
            int source = GetSourceClassification(111,pdgCode);
            fHistoMCSecPi0PtvsSource[fiCut]->Fill(particle->Pt(),source,tempParticleWeight); // All MC Pi0
            fHistoMCSecPi0Source[fiCut]->Fill(pdgCode);
          }else if(particle->GetPdgCode() == 221 && !fDoPi0Only){

            fHistoMCSecEtaPt[fiCut]->Fill(particle->Pt(),tempParticleWeight); // All MC Pi0
            fHistoMCSecEtaSource[fiCut]->Fill(pdgCode);
          }

          // check if conversion where within acceptance
          if( ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter0,fAODMCTrackArray) &&
              ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedAODMC(daughter1,fAODMCTrackArray)){
            if(particle->GetPdgCode() == 111){
              int source = GetSourceClassification(111,pdgCode);
              fHistoMCSecPi0InAccPtvsSource[fiCut]->Fill(particle->Pt(),source,tempParticleWeight); // All MC Pi0
            }
          }
        }
      }
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskPi0v2Calo::ProcessMCParticles(int isCurrentEventSelected){
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  double mcProdVtxX   = primVtxMC->GetX();
  double mcProdVtxY   = primVtxMC->GetY();
  double mcProdVtxZ   = primVtxMC->GetZ();

  // Check if MC generated particles should be filled for this event using the selected trigger
  if( !((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsMCTriggerSelected(fInputEvent, fMCEvent)){
    if(isCurrentEventSelected == 0) fHistoMCEventsTrigg[fiCut]->Fill(2., fWeightJetJetMC); // event is rejected and triggered
    else fHistoMCEventsTrigg[fiCut]->Fill(3., fWeightJetJetMC); // event is rejected but not triggered
    return;
  }
  if(isCurrentEventSelected == 0) fHistoMCEventsTrigg[fiCut]->Fill(0., fWeightJetJetMC); // event is accepted and triggered
  else fHistoMCEventsTrigg[fiCut]->Fill(1., fWeightJetJetMC); // event is accepted but not triggered

  // Loop over all primary MC particle
  for(long i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    if (((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, i, mcProdVtxX, mcProdVtxY, mcProdVtxZ)){
      double tempParticleWeight       = fWeightJetJetMC;

      AliMCParticle* particle = (AliMCParticle *)fMCEvent->GetTrack(i);
      if (!particle) continue;

      int isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
        // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
        if(isMCFromMBHeader == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempParticleWeight = 1;
      }
      if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(particle,fMCEvent)){
        if (isCurrentEventSelected==1) {
          fHistoMCGammaPtNotTriggered[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Gamma
        } else if (isCurrentEventSelected==2) {
          fHistoMCGammaPtNoVertex[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Gamma
        }
        fHistoMCAllGammaPt[fiCut]->Fill(particle->Pt(),fWeightJetJetMC); // All MC Gamma

        if(!fDoLightOutput){
          if(particle->GetMother() >-1){ // Meson Decay Gamma
            switch(fMCEvent->GetTrack(particle->GetMother())->PdgCode()){
            case 111: // Pi0
              fHistoMCDecayGammaPi0Pt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 113: // Rho0
              fHistoMCDecayGammaRhoPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 221: // Eta
              fHistoMCDecayGammaEtaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 223: // Omega
              fHistoMCDecayGammaOmegaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 331: // Eta'
              fHistoMCDecayGammaEtapPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 333: // Phi
              fHistoMCDecayGammaPhiPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            case 3212: // Sigma
              fHistoMCDecayGammaSigmaPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
              break;
            }
          }
        }
      }
      if(fDoMesonAnalysis){

        double mesonY = 1.e30;
        double ratio  = 0;
        if (particle->E() != TMath::Abs(particle->Pz())){
          ratio         = (particle->E()+particle->Pz()) / (particle->E()-particle->Pz());
        }
        if( !(ratio <= 0) ){
          mesonY = particle->Y()-((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift();
        }
        if(fHistoMCPrimaryPtvsSource[fiCut]){
          if ((mesonY > ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValueMin()) && (mesonY < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetRapidityCutValueMax())){
            if ( particle->PdgCode() == 211 ){  // positve pions
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 0., tempParticleWeight);
            } else if ( particle->PdgCode() == -211 ){  // negative pions
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 1., tempParticleWeight);
            } else if ( particle->PdgCode() == 321 ){  // positve kaons
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 2., tempParticleWeight);
            } else if ( particle->PdgCode() == -321 ){  // negative kaons
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 3., tempParticleWeight);
            } else if ( TMath::Abs(particle->PdgCode()) == 310 ){  // K0s
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 4., tempParticleWeight);
            } else if ( TMath::Abs(particle->PdgCode()) == 130 ){  // K0l
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 5., tempParticleWeight);
            } else if ( TMath::Abs(particle->PdgCode()) == 3122 ){  // Lambda/ AntiLambda
              fHistoMCPrimaryPtvsSource[fiCut]->Fill(particle->Pt(), 6., tempParticleWeight);
            }
          }
        }
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
          ->MesonIsSelectedMC(particle,fMCEvent,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()) && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonLeadTrackSelectionMC(fInputEvent, particle)){
          AliMCParticle* daughter0 = (AliMCParticle*)fMCEvent->GetTrack(particle->GetDaughterFirst());
          AliMCParticle* daughter1 = (AliMCParticle*)fMCEvent->GetTrack(particle->GetDaughterLast());

          float weighted= 1;
          if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
            if (particle->Pt()>0.005){
              weighted = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(i, fMCEvent, fInputEvent);
            }
          }

          if(particle->PdgCode() == 111){
            if (isCurrentEventSelected == 1) fHistoMCPi0PtNotTriggered[fiCut]->Fill(particle->Pt(),weighted*tempParticleWeight); // All MC Pi0 in not triggered collisions
            else if (isCurrentEventSelected == 2) fHistoMCPi0PtNoVertex[fiCut]->Fill(particle->Pt(),weighted*tempParticleWeight); // All MC Pi0 in not triggered collisions
            fHistoMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted*tempParticleWeight); // All MC Pi0
            fHistoMCPi0WOWeightPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
            if (fIsMC > 1)fHistoMCPi0WOEvtWeightPt[fiCut]->Fill(particle->Pt());
          } else if(particle->PdgCode() == 221 && !fDoPi0Only){
            if (isCurrentEventSelected == 1) fHistoMCEtaPtNotTriggered[fiCut]->Fill(particle->Pt(),weighted*tempParticleWeight); // All MC Eta in not triggered collisions
            else if (isCurrentEventSelected == 2) fHistoMCEtaPtNoVertex[fiCut]->Fill(particle->Pt(),weighted*tempParticleWeight); // All MC Eta in not triggered collisions
            fHistoMCEtaPt[fiCut]->Fill(particle->Pt(),weighted*tempParticleWeight); // All MC Eta
            fHistoMCEtaWOWeightPt[fiCut]->Fill(particle->Pt(), tempParticleWeight);
            if (fIsMC > 1)fHistoMCEtaWOEvtWeightPt[fiCut]->Fill(particle->Pt());
          }
          // Check the acceptance for both gammas & whether they are counted as primaries as well
          bool kDaughter0IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, particle->GetDaughterFirst(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
          bool kDaughter1IsPrim = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, particle->GetDaughterLast(), mcProdVtxX, mcProdVtxY, mcProdVtxZ);

          if( kDaughter0IsPrim && kDaughter1IsPrim &&
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCEvent) &&
            ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCEvent) ){
            if(particle->PdgCode() == 111){
              fHistoMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc
              if (isCurrentEventSelected==1) fHistoMCPi0InAccPtNotTriggered[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Pi0 with gamma in acc for events which are not triggered
              if(fIsMC > 1) fHistoMCPi0WOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0 with gamma in acc
            } else if(particle->PdgCode() == 221 && !fDoPi0Only){
              fHistoMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Eta with gamma in acc
              if(isCurrentEventSelected == 1) fHistoMCEtaInAccPtNotTriggered[fiCut]->Fill(particle->Pt(),weighted* tempParticleWeight); // MC Eta with gamma in acc for not triggered events
              if(fIsMC > 1) fHistoMCEtaWOEvtWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Eta with gamma in acc
            }
          }
        }
      }
    // fill secondary histograms
    } else {
      AliMCParticle* particle = (AliMCParticle *)fMCEvent->GetTrack(i);
      if (!particle) continue;
      double tempParticleWeight       = fWeightJetJetMC;
      int isMCFromMBHeader = -1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
        isMCFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent);
        if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
        // Set the jetjet weight to 1 in case the particle orignated from the minimum bias header
        if(isMCFromMBHeader == 2 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4) tempParticleWeight = 1;
      }

      if(!fDoLightOutput) {
        if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(particle,fMCEvent)){
          if (particle->GetMother() > -1 && fMCEvent->GetTrack(particle->GetMother())->GetMother() > -1) {
            if (fMCEvent->GetTrack(fMCEvent->GetTrack(particle->GetMother())->GetMother())->PdgCode() == 310){
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),0.,tempParticleWeight);
            } else if (fMCEvent->GetTrack(fMCEvent->GetTrack(particle->GetMother())->GetMother())->PdgCode() == 130) {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),1.,tempParticleWeight);
            } else if (fMCEvent->GetTrack(fMCEvent->GetTrack(particle->GetMother())->GetMother())->PdgCode() == 3122) {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),2.,tempParticleWeight);
            } else if (fMCEvent->GetTrack(fMCEvent->GetTrack(particle->GetMother())->GetMother())->PdgCode() == 221) {
              fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),3.,tempParticleWeight);
            } else {
              if ( !(TMath::Abs(fMCEvent->GetTrack(particle->GetMother())->PdgCode()) == 11 && fMCEvent->GetTrack(fMCEvent->GetTrack(particle->GetMother())->GetMother())->PdgCode() == 22) )
                fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),4.,tempParticleWeight);
            }
          } else {
            fHistoMCAllSecondaryGammaPt[fiCut]->Fill(particle->Pt(),4.,tempParticleWeight);
          }
        }
      }

      if(fDoMesonAnalysis){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMC(particle,fMCEvent,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()) && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonLeadTrackSelectionMC(fInputEvent, particle)){
          AliMCParticle* daughter0  = (AliMCParticle*)fMCEvent->GetTrack(particle->GetDaughterFirst());
          AliMCParticle* daughter1  = (AliMCParticle*)fMCEvent->GetTrack(particle->GetDaughterLast());
          int pdgCode = -1;
          if(particle->GetMother() > -1) pdgCode = ((AliMCParticle*)fMCEvent->GetTrack(particle->GetMother()))->PdgCode();

          if(particle->PdgCode() == 111){
            int source = GetSourceClassification(111,pdgCode);
            fHistoMCSecPi0PtvsSource[fiCut]->Fill(particle->Pt(),source,tempParticleWeight);
            fHistoMCSecPi0Source[fiCut]->Fill(pdgCode);
          } else if(particle->PdgCode() == 221 && !fDoPi0Only){
            fHistoMCSecEtaPt[fiCut]->Fill(particle->Pt(),tempParticleWeight);
            fHistoMCSecEtaSource[fiCut]->Fill(pdgCode);
          }

          // check if photons where within acceptance
          if( ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter0,fMCEvent) &&
                ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelectedMC(daughter1,fMCEvent)){
            if(particle->PdgCode() == 111){
              int source = GetSourceClassification(111,pdgCode);
              fHistoMCSecPi0InAccPtvsSource[fiCut]->Fill(particle->Pt(),source,tempParticleWeight);
            }
          }
        }
      }

    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2Calo::CalculatePi0Candidates(){

  // Conversion Gammas
  if(fClusterCandidates->GetEntries()>0){

    for(int firstGammaIndex=0;firstGammaIndex<fClusterCandidates->GetEntries();firstGammaIndex++){
      AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(firstGammaIndex));
      if (gamma0==NULL) continue;
      if ( fDoInOutTimingCluster ){
        double tof = fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef())->GetTOF();
        if ( tof < fMinTimingCluster || tof > fMaxTimingCluster ) continue;
      }
      for(int secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fClusterCandidates->GetEntries();secondGammaIndex++){
        AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fClusterCandidates->At(secondGammaIndex));
        if (gamma1==NULL) continue;
        if ( fDoInOutTimingCluster ){
          double tof = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef())->GetTOF();
          if ( tof > fMinTimingCluster && tof < fMaxTimingCluster ) continue;
        }

        double tempPi0CandWeight       = fWeightJetJetMC;
        // Set the pi0 candidate jetjet weight to 1 in case both photons orignated from the minimum bias header
        if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
          if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma0->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2 &&
              ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2)
            tempPi0CandWeight = 1;
        }

        if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaMinEnergyCut() ){
          float minDaughterEnergy = ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetSingleDaughterMinE();
          if( (gamma0->E() < minDaughterEnergy)  && (gamma1->E() < minDaughterEnergy)) { // at least one over threshold
            continue;
          }
        }
        AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);
        pi0cand->SetLabels(firstGammaIndex,secondGammaIndex);

        if(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetDoSecondaryTrackMatching()){
          TClonesArray * arrClustersProcess = NULL;
          arrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
          AliVCluster* Cluster0 = NULL;
          AliVCluster* Cluster1 = NULL;
          if (gamma0->GetIsCaloPhoton() > 0 && gamma1->GetIsCaloPhoton() > 0){
            if(fInputEvent->IsA()==AliESDEvent::Class()){
              if(arrClustersProcess){
                Cluster0 = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(gamma0->GetCaloClusterRef()));
                Cluster1 = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(gamma1->GetCaloClusterRef()));
              }else{
                Cluster0 = fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef());
                Cluster1 = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
              }
            } else if(fInputEvent->IsA()==AliAODEvent::Class()){
              if(arrClustersProcess){
                Cluster0 = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(gamma0->GetCaloClusterRef()));
                Cluster1 = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(gamma1->GetCaloClusterRef()));
              } else{
                Cluster0 = fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef());
                Cluster1 = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef());
              }
            }
          }
          bool ClusterMatched = kFALSE;
          for(int ConversionIndex=0;ConversionIndex<fGammaCandidates->GetEntries();ConversionIndex++){
            AliAODConversionPhoton *gammaConversion=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(ConversionIndex));
            if (gammaConversion==NULL) continue;
            bool matchedGamma0 =  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gammaConversion, Cluster0, fInputEvent, tempPi0CandWeight);
            bool matchedGamma1 =  ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->MatchConvPhotonToCluster(gammaConversion, Cluster1, fInputEvent, tempPi0CandWeight);
            if(matchedGamma0 || matchedGamma1) {
              ClusterMatched = kTRUE;
            }
          }
          if(arrClustersProcess){
              delete Cluster0;
              delete Cluster1;
          }
          if(ClusterMatched) continue;
        }
        if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),gamma0->GetLeadingCellID(),gamma1->GetLeadingCellID(), gamma0->GetIsCaloPhoton(), gamma1->GetIsCaloPhoton())) && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonLeadTrackSelection(fInputEvent, pi0cand)){
          if(fLocalDebugFlag == 1) DebugMethodPrint1(pi0cand,gamma0,gamma1);
          FillFlowHistograms(pi0cand);

          if(fIsMC> 0){
            if(fInputEvent->IsA()==AliESDEvent::Class())
              ProcessTrueMesonCandidates(pi0cand,gamma0,gamma1);
            if(fInputEvent->IsA()==AliAODEvent::Class())
              ProcessTrueMesonCandidatesAOD(pi0cand,gamma0,gamma1);
          }
        }
        delete pi0cand;
        pi0cand=0x0;
      }
    }
  }
}

//______________________________________________________________________
void AliAnalysisTaskPi0v2Calo::ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1){
  // Process True Mesons
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  double mcProdVtxX   = primVtxMC->GetX();
  double mcProdVtxY   = primVtxMC->GetY();
  double mcProdVtxZ   = primVtxMC->GetZ();

  double tempTruePi0CandWeight       = fWeightJetJetMC;

  bool isTruePi0         = kFALSE;
  bool isTrueEta         = kFALSE;
  int convertedPhotonLabel1    = -1;

  int gamma0MCLabel = TrueGammaCandidate0->GetCaloPhotonMCLabel(0);   // get most probable MC label
  int gamma0MotherLabel = -1;
  int tmpGammaMotherlabel = -1;

  AliMCParticle * gammaMC0 = 0x0;
  if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    gammaMC0 = (AliMCParticle*)fMCEvent->GetTrack(gamma0MCLabel);
    if (TrueGammaCandidate0->IsLargestComponentPhoton() || TrueGammaCandidate0->IsLargestComponentElectron()){    // largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      tmpGammaMotherlabel=gammaMC0->GetMother();
      if (TrueGammaCandidate0->IsLargestComponentPhoton()){                            // for photons its the direct mother
        gamma0MotherLabel=gammaMC0->GetMother();
      } else if (TrueGammaCandidate0->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate0->IsConversion() && (gammaMC0->GetMother() > -1)){
          gamma0MotherLabel=fMCEvent->GetTrack(gammaMC0->GetMother())->GetMother();
        } else {
          gamma0MotherLabel=gammaMC0->GetMother();
        }
      }
    }
  }
  if (TrueGammaCandidate1->GetIsCaloPhoton() == 0) AliFatal("CaloPhotonFlag has not been set. Aborting");

  bool previouslyNotFoundTrueMesons = kFALSE;
  int SaftyLoopCounter = 0;
  while (tmpGammaMotherlabel > 0 && SaftyLoopCounter < 100) {
      SaftyLoopCounter++;
      if(((AliMCParticle*)fMCEvent->GetTrack(tmpGammaMotherlabel))->PdgCode() != 111 && ((AliMCParticle*)fMCEvent->GetTrack(tmpGammaMotherlabel))->PdgCode() != 221) {
          tmpGammaMotherlabel = ((AliMCParticle*)fMCEvent->GetTrack(tmpGammaMotherlabel))->GetMother();
      } else {
          if (tmpGammaMotherlabel != gamma0MotherLabel) {
              previouslyNotFoundTrueMesons = kTRUE;
          }
          gamma0MotherLabel = tmpGammaMotherlabel;
          break;
      }
  }
  int gamma1MCLabel = TrueGammaCandidate1->GetCaloPhotonMCLabel(0);   // get most probable MC label
  int gamma1MotherLabel = -1;
  tmpGammaMotherlabel = -1;
  // check if

  AliMCParticle * gammaMC1 = 0x0;
  if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    // Daughters Gamma 1
    gammaMC1 = (AliMCParticle*)fMCEvent->GetTrack(gamma1MCLabel);
    if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){    // largest component is electro magnetic
      // get mother of interest (pi0 or eta)
      tmpGammaMotherlabel = gammaMC1->GetMother();
      if (TrueGammaCandidate1->IsLargestComponentPhoton()){                            // for photons its the direct mother
        gamma1MotherLabel = gammaMC1->GetMother();
      } else if (TrueGammaCandidate1->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate1->IsConversion()){
          convertedPhotonLabel1 = gammaMC1->GetMother();
          if(convertedPhotonLabel1 > -1) gamma1MotherLabel = fMCEvent->GetTrack(convertedPhotonLabel1)->GetMother();
        } else {
          gamma1MotherLabel = gammaMC1->GetMother();
        }
      }
    }
  }
  SaftyLoopCounter = 0;
  while (tmpGammaMotherlabel > 0 && SaftyLoopCounter < 100) {
      SaftyLoopCounter++;
      if(((AliMCParticle*)fMCEvent->GetTrack(tmpGammaMotherlabel))->PdgCode() != 111 && ((AliMCParticle*)fMCEvent->GetTrack(tmpGammaMotherlabel))->PdgCode() != 221) {
          tmpGammaMotherlabel = ((AliMCParticle*)fMCEvent->GetTrack(tmpGammaMotherlabel))->GetMother();
      } else {
          if (tmpGammaMotherlabel != gamma1MotherLabel) {
              previouslyNotFoundTrueMesons = kTRUE;
          }
          gamma1MotherLabel = tmpGammaMotherlabel;
          break;
      }
  }

  if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
    if(((AliMCParticle*)fMCEvent->GetTrack(gamma1MotherLabel))->PdgCode() == 111){
      isTruePi0=kTRUE;
    }
    if(((AliMCParticle*)fMCEvent->GetTrack(gamma1MotherLabel))->PdgCode() == 221){
      isTrueEta=kTRUE;
    }
  }

  // Set the pi0 candidate jetjet weight to 1 in case both photons orignated from the minimum bias header
  if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
    if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TrueGammaCandidate0->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2 &&
        ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TrueGammaCandidate1->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2)
      tempTruePi0CandWeight = 1;
  }

  if(isTruePi0 || isTrueEta){// True Pion or Eta
    if (isTruePi0){
      fHistoTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  tempTruePi0CandWeight);
      if(previouslyNotFoundTrueMesons && !fDoLightOutput && fHistoTruePi0InvMassPtAdditional[fiCut]) fHistoTruePi0InvMassPtAdditional[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),tempTruePi0CandWeight);
      
      if (!fDoLightOutput && TMath::Abs(Pi0Candidate->GetAlpha())< 0.1){
        fHistoTruePi0InvMassPtAlpha[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  tempTruePi0CandWeight);
        if (TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate0->IsLargestComponentPhoton())
          fHistoTruePi0PureGammaInvMassPtAlpha[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  tempTruePi0CandWeight);
      }
    }
    if (isTrueEta && !fDoPi0Only){
      fHistoTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  tempTruePi0CandWeight);
      if(previouslyNotFoundTrueMesons && !fDoLightOutput && fHistoTrueEtaInvMassPtAdditional[fiCut]) fHistoTrueEtaInvMassPtAdditional[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),tempTruePi0CandWeight);
    }

    bool isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCEvent, gamma0MotherLabel, mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if(!isPrimary){ // Secondary Meson
      // filling secondary histograms
      long secMotherLabel = -1;
      if(gamma0MotherLabel > -1) secMotherLabel = ((AliMCParticle*)fMCEvent->GetTrack(gamma0MotherLabel))->GetMother();

      float weightedSec= 1;
      if((secMotherLabel > -1) && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, fMCEvent, fInputEvent) && fMCEvent->GetTrack(secMotherLabel)->PdgCode()==310){
        weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, fMCEvent, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
        //cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
      }
      if (isTruePi0) {
        fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
      }
      if (secMotherLabel >-1){
        if(fMCEvent->GetTrack(secMotherLabel)->PdgCode() == 310 && isTruePi0){
          fHistoTrueSecondaryPi0FromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
        } else if(fMCEvent->GetTrack(secMotherLabel)->PdgCode()==130 && isTruePi0){
          fHistoTrueSecondaryPi0FromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
        } else if(fMCEvent->GetTrack(secMotherLabel)->PdgCode()==221 && isTruePi0){
          fHistoTrueSecondaryPi0FromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
        } else if(fMCEvent->GetTrack(secMotherLabel)->PdgCode()==3122 && isTruePi0){
          fHistoTrueSecondaryPi0FromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
        } else if (isTruePi0){
        } else if (isTrueEta && !fDoPi0Only){
        }
      }
    } else { // Only primary pi0 for efficiency calculation
      // filling primary histograms
      float weighted= 1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, fMCEvent, fInputEvent)){
        if (((AliMCParticle*)fMCEvent->GetTrack(gamma1MotherLabel))->Pt()>0.005){
          weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, fMCEvent, fInputEvent);
          // cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
        }
      }
      if (isTruePi0){
        fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
        fHistoTruePrimaryPi0W0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
        fProfileTruePrimaryPi0WeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
        if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), weighted*tempTruePi0CandWeight);
      } else if (isTrueEta && !fDoPi0Only){
        fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
        fHistoTruePrimaryEtaW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
        fProfileTruePrimaryEtaWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
        if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,gamma0MotherLabel)) fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), weighted*tempTruePi0CandWeight);
      }

    }
  }

}
//______________________________________________________________________
void AliAnalysisTaskPi0v2Calo::ProcessTrueMesonCandidatesAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1){
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  double mcProdVtxX   = primVtxMC->GetX();
  double mcProdVtxY   = primVtxMC->GetY();
  double mcProdVtxZ   = primVtxMC->GetZ();

  double tempTruePi0CandWeight       = fWeightJetJetMC;

  // Process True Mesons
  if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (fAODMCTrackArray == NULL) return;

  bool isTruePi0         = kFALSE;
  bool isTrueEta         = kFALSE;

  int gamma0MCLabel       = TrueGammaCandidate0->GetCaloPhotonMCLabel(0);   // get most probable MC label
  int gamma0MotherLabel     = -1;
  int tmpGammaMotherlabel = -1;

  // check if
  AliAODMCParticle * gammaMC0 = 0x0;
  if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    // Daughters Gamma 0
    gammaMC0 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MCLabel));
    if (TrueGammaCandidate0->IsLargestComponentPhoton() || TrueGammaCandidate0->IsLargestComponentElectron()){    // largest component is electro magnetic
      tmpGammaMotherlabel = gammaMC0->GetMother();
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate0->IsLargestComponentPhoton()){                            // for photons its the direct mother
        gamma0MotherLabel=gammaMC0->GetMother();
      } else if (TrueGammaCandidate0->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate0->IsConversion()){
          AliAODMCParticle * gammaGrandMotherMC0 =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMC0->GetMother()));
          gamma0MotherLabel=gammaGrandMotherMC0->GetMother();
        } else gamma0MotherLabel=gammaMC0->GetMother();
      }
    }
  }
  bool previouslyNotFoundTrueMesons = kFALSE;
  int SaftyLoopCounter = 0;
  while (tmpGammaMotherlabel > 0 && SaftyLoopCounter < 100) {
      SaftyLoopCounter++;
      if(((AliAODMCParticle*)fAODMCTrackArray->At(tmpGammaMotherlabel))->GetPdgCode() != 111 && ((AliAODMCParticle*)fAODMCTrackArray->At(tmpGammaMotherlabel))->GetPdgCode() != 221) {
          tmpGammaMotherlabel = ((AliAODMCParticle*)fAODMCTrackArray->At(tmpGammaMotherlabel))->GetMother();
      } else {
          if (tmpGammaMotherlabel != gamma0MotherLabel) {
              previouslyNotFoundTrueMesons = kTRUE;
          }
          gamma0MotherLabel = tmpGammaMotherlabel;
          break;
      }
  }

  int gamma1MCLabel         = TrueGammaCandidate1->GetCaloPhotonMCLabel(0);   // get most probable MC label
  int gamma1MotherLabel     = -1;
  tmpGammaMotherlabel = -1;

  // check if
  AliAODMCParticle *gammaMC1  = 0x0;
  if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
    // Daughters Gamma 1
    gammaMC1 = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MCLabel));
    if (TrueGammaCandidate1->IsLargestComponentPhoton() || TrueGammaCandidate1->IsLargestComponentElectron()){    // largest component is electro magnetic
      tmpGammaMotherlabel = gammaMC1->GetMother();
      // get mother of interest (pi0 or eta)
      if (TrueGammaCandidate1->IsLargestComponentPhoton()){                            // for photons its the direct mother
        gamma1MotherLabel=gammaMC1->GetMother();
      } else if (TrueGammaCandidate1->IsLargestComponentElectron()){                         // for electrons its either the direct mother or for conversions the grandmother
        if (TrueGammaCandidate1->IsConversion()){
          AliAODMCParticle * gammaGrandMotherMC1 =  static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gammaMC1->GetMother()));
          gamma1MotherLabel=gammaGrandMotherMC1->GetMother();
        } else gamma1MotherLabel=gammaMC1->GetMother();
      }
    }
  }

  SaftyLoopCounter = 0;
  while (tmpGammaMotherlabel > 0 && SaftyLoopCounter < 100) {
      SaftyLoopCounter++;
      if(((AliAODMCParticle*)fAODMCTrackArray->At(tmpGammaMotherlabel))->GetPdgCode() != 111 && ((AliAODMCParticle*)fAODMCTrackArray->At(tmpGammaMotherlabel))->GetPdgCode() != 221) {
          tmpGammaMotherlabel = ((AliAODMCParticle*)fAODMCTrackArray->At(tmpGammaMotherlabel))->GetMother();
      } else {
          if (tmpGammaMotherlabel != gamma1MotherLabel) {
              previouslyNotFoundTrueMesons = kTRUE;
          }
          gamma1MotherLabel = tmpGammaMotherlabel;
          break;
      }
  }

  if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
    if(((AliAODMCParticle*)fAODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
      isTruePi0=kTRUE;
    }
    if(((AliAODMCParticle*)fAODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 221){
      isTrueEta=kTRUE;
    }
  }

  // Set the pi0 candidate jetjet weight to 1 in case both photons orignated from the minimum bias header
  if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
    if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TrueGammaCandidate0->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2 &&
        ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TrueGammaCandidate1->GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2)
      tempTruePi0CandWeight = 1;
  }

  if(isTruePi0 || isTrueEta){// True Pion or Eta
    if(isTruePi0){
      fHistoTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
      if(previouslyNotFoundTrueMesons && !fDoLightOutput && fHistoTruePi0InvMassPtAdditional[fiCut]) fHistoTruePi0InvMassPtAdditional[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),tempTruePi0CandWeight);
      if (!fDoLightOutput && TMath::Abs(Pi0Candidate->GetAlpha())< 0.1){
        fHistoTruePi0InvMassPtAlpha[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  tempTruePi0CandWeight);
        if (TrueGammaCandidate1->IsLargestComponentPhoton() && TrueGammaCandidate0->IsLargestComponentPhoton())
          fHistoTruePi0PureGammaInvMassPtAlpha[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),  tempTruePi0CandWeight);
      }
    }
    if(isTrueEta && !fDoPi0Only){
      fHistoTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
      if(previouslyNotFoundTrueMesons && !fDoLightOutput && fHistoTrueEtaInvMassPtAdditional[fiCut]) fHistoTrueEtaInvMassPtAdditional[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),tempTruePi0CandWeight);
    }

    bool isPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryAOD(fInputEvent, static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma0MotherLabel)), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if(!isPrimary){ // Secondary Meson
      long secMotherLabel = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->GetMother();
      float weightedSec= 1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, 0x0, fInputEvent) && static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
        weightedSec= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(secMotherLabel, 0x0, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
        //cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
      }
      if (isTruePi0) {
        fHistoTrueSecondaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
      }
      if (secMotherLabel >-1){
        if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310 && isTruePi0 ){
          fHistoTrueSecondaryPi0FromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
        }
        if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==130 && isTruePi0 ){
          fHistoTrueSecondaryPi0FromK0lInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
        }
        if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==221 && isTruePi0){
          fHistoTrueSecondaryPi0FromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
        }
        if(static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(secMotherLabel))->GetPdgCode()==3122 && isTruePi0){
          fHistoTrueSecondaryPi0FromLambdaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec* tempTruePi0CandWeight);
        }
      }
    } else{ // Only primary pi0 for efficiency calculation
      float weighted= 1;
      if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, 0x0, fInputEvent)){
        if (static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(gamma1MotherLabel))->Pt()>0.005){
        weighted= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetWeightForMeson(gamma1MotherLabel, 0x0, fInputEvent);
        //                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
        }
      }
      if (isTruePi0){
        fHistoTruePrimaryPi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
        fHistoTruePrimaryPi0W0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
        fProfileTruePrimaryPi0WeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
        if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gamma0MotherLabel)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), weighted*tempTruePi0CandWeight);
      } else if (isTrueEta && !fDoPi0Only){
        fHistoTruePrimaryEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
        fHistoTruePrimaryEtaW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), tempTruePi0CandWeight);
        fProfileTruePrimaryEtaWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted* tempTruePi0CandWeight);
        if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,gamma0MotherLabel)) fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(), weighted*tempTruePi0CandWeight);
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2Calo::CalculateBackground(){

  int zbin= 0;
  if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoSectorMixing() ) {
    zbin = fBGHandler[fiCut]->GetZBinIndex(fV0Reader->GetPtMaxSector());
  } else {
    zbin = fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
  }
  int mbin = 0;

  double tempBGCandidateWeight       = fWeightJetJetMC;

  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
    mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
  } else {
    mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fClusterCandidates->GetEntries());
  }

  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
    for(int nEventsInBG=0;nEventsInBG<fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
      AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
      for(int iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
        AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
        for(unsigned int iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
          AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
          AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
          backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());

          // Set the BG candidate jetjet weight to 1 in case both photons orignated from the minimum bias header
          if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
            if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(previousGoodV0.GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2 &&
                ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(currentEventGoodV0.GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2)
              tempBGCandidateWeight = 1;
          }

          if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
            ->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()),currentEventGoodV0.GetLeadingCellID(),previousGoodV0.GetLeadingCellID(),currentEventGoodV0.GetIsCaloPhoton(), previousGoodV0.GetIsCaloPhoton() )  && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonLeadTrackSelection(fInputEvent, backgroundCandidate)){
            FillFlowBackHistograms(backgroundCandidate, tempBGCandidateWeight);
          }
          delete backgroundCandidate;
          backgroundCandidate = 0x0;
        }
      }
    }
  } else if ( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UsePtmaxMethod() ) {

    double currentPtMax   = 0;  double previousPtMax  = 0;
    double currentAvePt   = 0;  double previousAvePt  = 0;
    double currentAveEta  = 0;  double previousAveEta = 0;
    double currentAvePhi  = 0;  double previousAvePhi = 0;
    bool acceptedPtMax    = kFALSE;

    for(int nEventsInBG=0;nEventsInBG <fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
      AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
      if(previousEventV0s){
        acceptedPtMax = kFALSE;
        currentPtMax = 0; previousPtMax = 0;
        currentAvePt = 0; previousAvePt = 0; currentAveEta = 0; previousAveEta = 0; currentAvePhi = 0; previousAvePhi = 0;
        for(int iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
          AliAODConversionPhoton *currentV0 = (AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
          currentAvePt += currentV0->GetPhotonPt();
          currentAveEta += currentV0->GetPhotonPt()*currentV0->GetPhotonEta();
          currentAvePhi += currentV0->GetPhotonPt()*currentV0->GetPhotonPhi();
          if(currentV0->GetPhotonPt() > currentPtMax){ currentPtMax = currentV0->GetPhotonPt(); }
        }
        currentAveEta /= currentAvePt;
        currentAvePhi /= currentAvePt;
        currentAvePt /= fClusterCandidates->GetEntries();
        for(unsigned int iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
            AliAODConversionPhoton *previousV0 = (AliAODConversionPhoton*)(previousEventV0s->at(iPrevious));
            previousAvePt += previousV0->GetPhotonPt();
            previousAveEta += previousV0->GetPhotonPt()*previousV0->GetPhotonEta();
            previousAvePhi += previousV0->GetPhotonPt()*previousV0->GetPhotonPhi();
            if(previousV0->GetPhotonPt() > previousPtMax){ previousPtMax = previousV0->GetPhotonPt(); }
        }
        previousAveEta /= previousAvePt;
        previousAvePhi /= previousAvePt;
        previousAvePt /= previousEventV0s->size();
        if(currentPtMax > 0. && previousPtMax > 0.){
         //if(TMath::Sqrt(pow((currentEta-previousEta),2)+pow((currentPhi-previousPhi),2)) < 0.2) acceptedPtMax = kTRUE;
         if(TMath::Abs(previousAveEta-currentAveEta)<0.4 && TMath::Abs(previousAvePhi-currentAvePhi)<0.6 && (previousAvePt/currentAvePt)<4. && (previousAvePt/currentAvePt)>0.25) acceptedPtMax = kTRUE;
        }
        if(acceptedPtMax){
          for(int iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
            AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
            for(unsigned int iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
              AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
              AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
              backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());

              // Set the BG candidate jetjet weight to 1 in case both photons orignated from the minimum bias header
              if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
                if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(previousGoodV0.GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2 &&
                    ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(currentEventGoodV0.GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2)
                  tempBGCandidateWeight = 1;
              }

              if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
                ->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift()),currentEventGoodV0.GetLeadingCellID(),previousGoodV0.GetLeadingCellID(), currentEventGoodV0.GetIsCaloPhoton(), previousGoodV0.GetIsCaloPhoton())  && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonLeadTrackSelection(fInputEvent, backgroundCandidate)){
                FillFlowBackHistograms(backgroundCandidate, tempBGCandidateWeight);
              }
              delete backgroundCandidate;
              backgroundCandidate = 0x0;
            }
          }
        }
      }
    }
  } else if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoSectorMixing() ) {
    if(fClusterCandidates->GetEntries()>0){
      for(int nEventsInBG=0;nEventsInBG <fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
        AliGammaConversionAODVector *previousEventV0s = NULL;
        previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
        if(previousEventV0s && previousEventV0s->size()>0){
              for(int iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
                AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
                for(unsigned int iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
                  AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
                  AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
                  backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());

                  if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),currentEventGoodV0.GetLeadingCellID(),previousGoodV0.GetLeadingCellID(), currentEventGoodV0.GetIsCaloPhoton(), previousGoodV0.GetIsCaloPhoton()))  && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonLeadTrackSelection(fInputEvent, backgroundCandidate)){
                    // Set the BG candidate jetjet weight to 1 in case both photons orignated from the minimum bias header
                    if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
                      if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(previousGoodV0.GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2 &&
                          ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(currentEventGoodV0.GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2)
                        tempBGCandidateWeight = 1;
                    }
                    FillFlowBackHistograms(backgroundCandidate, tempBGCandidateWeight);
                  }
                  delete backgroundCandidate;
                  backgroundCandidate = 0x0;
                }
              }
        }
      }
    }
  } else {
    for(int nEventsInBG=0;nEventsInBG <fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
      AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
      if(previousEventV0s){
        for(int iCurrent=0;iCurrent<fClusterCandidates->GetEntries();iCurrent++){
          AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent));
          for(unsigned int iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){

            AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
            std::unique_ptr<AliAODConversionMother> backgroundCandidate (new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0));
            backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());

            if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(),currentEventGoodV0.GetLeadingCellID(),previousGoodV0.GetLeadingCellID(),currentEventGoodV0.GetIsCaloPhoton(),  previousGoodV0.GetIsCaloPhoton() )) && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonLeadTrackSelection(fInputEvent, backgroundCandidate.get())){
              // Set the BG candidate jetjet weight to 1 in case both photons orignated from the minimum bias header
              if (fIsMC>0 && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() == 4){
                if( ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(previousGoodV0.GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2 &&
                    ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(currentEventGoodV0.GetCaloPhotonMCLabel(0), fMCEvent, fInputEvent) == 2)
                  tempBGCandidateWeight = 1;
              }
              FillFlowBackHistograms(backgroundCandidate.get(), tempBGCandidateWeight);
            }
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2Calo::CalculateBackgroundSwapp(){

  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoGammaSwappForBg()) {

    double rotationAngle = TMath::Pi()/2.0; //0.78539816339; // rotaion angle 90°

    TLorentzVector lvRotationPhoton1;   // photon candidates which get rotated
    TLorentzVector lvRotationPhoton2;   // photon candidates which get rotated
    TVector3 lvRotationPion;            // reconstructed mother particle from the two photons
    // Needed for TGenPhaseSpace
    TVector3 tvEtaPhigamma1, tvEtaPhigamma2, tvEtaPhigamma1Decay, tvEtaPhigamma2Decay, tvNormBeforeDecay, tvNormAfterDecay;
    float asymBeforeDecay = 0.;
    float asymAfterDecay = 0.;
    double massGamma[2] = {0,0};

    int cellIDRotatedPhoton1 = -1; // cell ID of the cluster after rotation
    int cellIDRotatedPhoton2 = -1; // cell ID of the cluster after rotation

    std::vector<std::array<double, 2>> vSwappingInvMassPT;
    std::vector<std::array<double, 2>> vSwappingInvMassEAlphaCut;
    vSwappingInvMassPT.clear();
    vSwappingInvMassEAlphaCut.clear();
    vSwappingInvMassPT.resize(0);
    vSwappingInvMassEAlphaCut.resize(0);

    // curcial requierment is that the event has at least 3 cluster candidates
    if(fClusterCandidates->GetEntries() > 2 ){

      for(int iCurrent1=0;iCurrent1<fClusterCandidates->GetEntries();iCurrent1++){
        AliAODConversionPhoton* currentEventGoodV0Temp1 = (AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent1));

        for(int iCurrent2=iCurrent1+1;iCurrent2<fClusterCandidates->GetEntries();iCurrent2++){
          AliAODConversionPhoton* currentEventGoodV0Temp2 = (AliAODConversionPhoton*)(fClusterCandidates->At(iCurrent2));

          for(int iSwapp = 0; iSwapp < ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfSwappsForBg(); ++iSwapp){

            lvRotationPhoton1.SetX(currentEventGoodV0Temp1->Px());
            lvRotationPhoton1.SetY(currentEventGoodV0Temp1->Py());
            lvRotationPhoton1.SetZ(currentEventGoodV0Temp1->Pz());
            lvRotationPhoton1.SetE(currentEventGoodV0Temp1->E());

            lvRotationPhoton2.SetX(currentEventGoodV0Temp2->Px());
            lvRotationPhoton2.SetY(currentEventGoodV0Temp2->Py());
            lvRotationPhoton2.SetZ(currentEventGoodV0Temp2->Pz());
            lvRotationPhoton2.SetE(currentEventGoodV0Temp2->E());

            lvRotationPion = (lvRotationPhoton1 + lvRotationPhoton2).Vect();

            // rotate both photons around the momentum vector of their hypothetical mother particle
            if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 0 || ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 1)){
              if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 0) rotationAngle = TMath::Pi()/2.0; // rotate by 90 degree
              else if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 1){  // rotate by random angle between
                 double temp = (fRandom.Rndm() < 0.5) ? 0 : TMath::Pi();
                 rotationAngle = temp + TMath::Pi()/3.0 + fRandom.Rndm()*TMath::Pi()/3.0;
              }
              lvRotationPhoton1.Rotate(rotationAngle, lvRotationPion);
              lvRotationPhoton2.Rotate(rotationAngle, lvRotationPion);
            } else if (((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() >=10){ // generate new decay with TGenPhaseSpace
              if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 11){
                tvEtaPhigamma1 = lvRotationPhoton1.Vect();
                tvEtaPhigamma2 = lvRotationPhoton2.Vect();
                tvNormBeforeDecay = tvEtaPhigamma1.Cross(tvEtaPhigamma2);
                asymBeforeDecay = fabs((lvRotationPhoton1.E()-lvRotationPhoton2.E())/(lvRotationPhoton1.E()+lvRotationPhoton2.E()));
              }

              TLorentzVector lvRotationMother = lvRotationPhoton1 + lvRotationPhoton2;
              fGenPhaseSpace.SetDecay(lvRotationMother, 2, massGamma);
              fGenPhaseSpace.Generate();
              lvRotationPhoton1 = *fGenPhaseSpace.GetDecay(0);
              lvRotationPhoton2 = *fGenPhaseSpace.GetDecay(1);

              if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GammaSwappMethodBg() == 11){
                tvEtaPhigamma1Decay = lvRotationPhoton1.Vect();
                tvEtaPhigamma2Decay = lvRotationPhoton2.Vect();
                tvNormAfterDecay = tvEtaPhigamma1Decay.Cross(tvEtaPhigamma2Decay);  // norm vector to decay plane
                asymAfterDecay = fabs((lvRotationPhoton1.E()-lvRotationPhoton2.E())/(lvRotationPhoton1.E()+lvRotationPhoton2.E()));
                // check if decay is nearly the same as original decay: if yes continue with next decay
                if((tvNormAfterDecay.Angle(tvNormBeforeDecay) < 20*TMath::Pi()/180. || tvNormAfterDecay.Angle(tvNormBeforeDecay) > 340*TMath::Pi()/180.) && ( fabs(asymBeforeDecay - asymAfterDecay) < 0.05 )   ) continue;
              }

            }


            cellIDRotatedPhoton1 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton1.Eta(), static_cast<double>((lvRotationPhoton1.Phi()<0) ? lvRotationPhoton1.Phi() + TMath::Pi()*2. : lvRotationPhoton1.Phi()));
            cellIDRotatedPhoton2 = ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetCaloCellIdFromEtaPhi(lvRotationPhoton2.Eta(), static_cast<double>((lvRotationPhoton2.Phi()<0) ? lvRotationPhoton2.Phi() + TMath::Pi()*2. : lvRotationPhoton2.Phi()));

            if(!fDoLightOutput){
              if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton1, fInputEvent, ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetDistanceToBorderForBg()))){
                ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillEtaPhiMapForClusterInBg(lvRotationPhoton1.Eta(), static_cast<double>((lvRotationPhoton1.Phi()<0) ? lvRotationPhoton1.Phi() + TMath::Pi()*2. : lvRotationPhoton1.Phi()), fWeightJetJetMC);
              }
              if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton2, fInputEvent, ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetDistanceToBorderForBg()))){
                ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->FillEtaPhiMapForClusterInBg(lvRotationPhoton2.Eta(), static_cast<double>((lvRotationPhoton2.Phi()<0) ? lvRotationPhoton2.Phi() + TMath::Pi()*2. : lvRotationPhoton2.Phi()), fWeightJetJetMC);
              }
            }

            std::unique_ptr<AliAODConversionPhoton> currentEventGoodV0Rotation1 (new AliAODConversionPhoton(&lvRotationPhoton1));
            std::unique_ptr<AliAODConversionPhoton> currentEventGoodV0Rotation2 (new AliAODConversionPhoton(&lvRotationPhoton2));

            for(auto const& kCurrentClusterCandidates  : *fClusterCandidates){
              if(currentEventGoodV0Temp1 == ((AliAODConversionPhoton*) kCurrentClusterCandidates) || currentEventGoodV0Temp2 == ((AliAODConversionPhoton*) kCurrentClusterCandidates)){ continue;}

              std::unique_ptr<AliAODConversionMother> backgroundCandidate1(new AliAODConversionMother(currentEventGoodV0Rotation1.get(), ((AliAODConversionPhoton*) kCurrentClusterCandidates)));
              std::unique_ptr<AliAODConversionMother> backgroundCandidate2(new AliAODConversionMother(currentEventGoodV0Rotation2.get(), ((AliAODConversionPhoton*) kCurrentClusterCandidates)));

              if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton1, fInputEvent, ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetDistanceToBorderForBg())) && lvRotationPhoton1.E() > ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetMinClusterEnergy())
              {
                if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate1.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), cellIDRotatedPhoton1, ((AliAODConversionPhoton*) kCurrentClusterCandidates)->GetLeadingCellID())  && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonLeadTrackSelection(fInputEvent, backgroundCandidate1.get()))
                {
                  vSwappingInvMassPT.push_back({backgroundCandidate1->M(),backgroundCandidate1->Pt()});
                }
              }
              if(!(((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->CheckDistanceToBadChannelSwapping(cellIDRotatedPhoton2, fInputEvent, ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetDistanceToBorderForBg())) && lvRotationPhoton2.E() > ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetMinClusterEnergy())
              {
                if(((AliConversionMesonCuts*) fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate2.get(),kFALSE,((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetEtaShift(), cellIDRotatedPhoton2, ((AliAODConversionPhoton*) kCurrentClusterCandidates)->GetLeadingCellID()) && ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonLeadTrackSelection(fInputEvent, backgroundCandidate2.get()))
                {
                  vSwappingInvMassPT.push_back({backgroundCandidate2->M(),backgroundCandidate2->Pt()});
                }
              }
            }
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskPi0v2Calo::UpdateEventByEventData(){
  //see header file for documentation
  if(fClusterCandidates->GetEntries() >1 ){
    if( ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoSectorMixing() ){
      fBGHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fV0Reader->GetPtMaxSector(),fClusterCandidates->GetEntries(),fEventPlaneAngle);
    } else if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
      fBGHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
    } else { // means we use #V0s for multiplicity
      fBGHandler[fiCut]->AddEvent(fClusterCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fClusterCandidates->GetEntries(),fEventPlaneAngle);
    }
  }
}

//_________________________________________________________________________________
bool AliAnalysisTaskPi0v2Calo::CheckVectorForDoubleCount(vector<int> &vec, int tobechecked)
{
  if(tobechecked > -1)
  {
    vector<int>::iterator it;
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
void AliAnalysisTaskPi0v2Calo::Terminate(const Option_t *)
{

  //fOutputContainer->Print(); // Will crash on GRID
}

//________________________________________________________________________
int AliAnalysisTaskPi0v2Calo::GetSourceClassification(int daughter, int pdgCode){

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
void AliAnalysisTaskPi0v2Calo::FillMultipleCountHistoAndClear(map<int,int> &ma, TH1F* hist){
  map<int, int>::iterator it;
  for (it = ma.begin(); it != ma.end(); it++){
    hist->Fill(it->second, fWeightJetJetMC);
  }
  ma.clear();
  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskPi0v2Calo::DebugMethodPrint1(AliAODConversionMother *pi0cand, AliAODConversionPhoton *gamma0, AliAODConversionPhoton *gamma1){
  if(pi0cand->GetOpeningAngle() < 0.02) DebugMethod(pi0cand,gamma0,gamma1);
  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskPi0v2Calo::DebugMethod(AliAODConversionMother *pi0cand, AliAODConversionPhoton *gamma0, AliAODConversionPhoton *gamma1){
  if(fLocalDebugFlag != 1) return;
  fstream fOutputLocalDebug;
  fOutputLocalDebug.open("debugOutput.txt",ios::out|ios::app);
  if(!fOutputLocalDebug.is_open()) return;

  AliVCaloCells* cells = fInputEvent->GetEMCALCells();
  fOutputLocalDebug << "--pi0cand--" << endl;
  fOutputLocalDebug << "openingAngle " << pi0cand->GetOpeningAngle() << endl;
  fOutputLocalDebug << "pT " << pi0cand->Pt() << endl;
  fOutputLocalDebug << "invMass " << pi0cand->M() << endl;
  fOutputLocalDebug << "--cluster1--" << endl;
  float clusPos[3]={0,0,0};
  fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef())->GetPosition(clusPos);
  TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
  double etaCluster = clusterVector.Eta();
  double phiCluster = clusterVector.Phi();

  int nCellCluster = fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef())->GetNCells();
  for(int iCell=0;iCell<nCellCluster;iCell++){
    int nSupMod, nModule, nIphi, nIeta;
    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetGeomEMCAL()->GetCellIndex(
          fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef())->GetCellAbsId(iCell),
          nSupMod,
          nModule,
          nIphi,
          nIeta);
    int nphi, neta;
    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetGeomEMCAL()->GetCellPhiEtaIndexInSModule(
          nSupMod,
          nModule,
          nIphi,
          nIeta,
          nphi,
          neta);
    fOutputLocalDebug << fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef())->GetCellAbsId(iCell) << " " << nSupMod << " " << nphi << " " << neta << " " << cells->GetCellAmplitude(fInputEvent->GetCaloCluster(gamma0->GetCaloClusterRef())->GetCellAbsId(iCell)) << endl;
  }
  fOutputLocalDebug << "phi " << phiCluster << endl;
  fOutputLocalDebug << "eta " << etaCluster << endl;

  fOutputLocalDebug << "--cluster2--" << endl;
  fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef())->GetPosition(clusPos);
  TVector3 clusterVector2(clusPos[0],clusPos[1],clusPos[2]);
  double etaCluster2 = clusterVector2.Eta();
  double phiCluster2 = clusterVector2.Phi();
  nCellCluster = fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef())->GetNCells();
  for(int iCell=0;iCell<nCellCluster;iCell++){
    int nSupMod, nModule, nIphi, nIeta;
    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetGeomEMCAL()->GetCellIndex(
          fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef())->GetCellAbsId(iCell),
          nSupMod,
          nModule,
          nIphi,
          nIeta);
    int nphi, neta;
    ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetGeomEMCAL()->GetCellPhiEtaIndexInSModule(
          nSupMod,
          nModule,
          nIphi,
          nIeta,
          nphi,
          neta);
    fOutputLocalDebug << fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef())->GetCellAbsId(iCell) << " " << nSupMod << " " << nphi << " " << neta << " " << cells->GetCellAmplitude(fInputEvent->GetCaloCluster(gamma1->GetCaloClusterRef())->GetCellAbsId(iCell)) << endl;
  }
  fOutputLocalDebug << "phi " << phiCluster2 << endl;
  fOutputLocalDebug << "eta " << etaCluster2 << endl;
  fOutputLocalDebug << "---------------------------------------" << endl;
  fOutputLocalDebug.close();

  return;
}

//_________________________________________________________________________________
void AliAnalysisTaskPi0v2Calo::EventDebugMethod(){
  if(fLocalDebugFlag != 2) return;
  fstream fOutputLocalDebug;
  fOutputLocalDebug.open("debugOutput.txt",ios::out|ios::app);
  if(!fOutputLocalDebug.is_open()) return;

  int nclus = 0;
  TClonesArray * arrClustersDebug = NULL;
  if(!fCorrTaskSetting.CompareTo("")){
    nclus = fInputEvent->GetNumberOfCaloClusters();
  } else {
    arrClustersDebug = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!arrClustersDebug)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskPi0v2Calo! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclus = arrClustersDebug->GetEntries();
  }

  AliVCaloCells* cells = fInputEvent->GetEMCALCells();
  fOutputLocalDebug << "--event--" << endl;
  fOutputLocalDebug << "nclusters " << nclus << endl;

  // vertex
  double vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  for(long i = 0; i < nclus; i++){
    AliVCluster* clus = NULL;
    if(fInputEvent->IsA()==AliESDEvent::Class()){
      if(arrClustersDebug)
        clus = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersDebug->At(i));
      else
        clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
    } else if(fInputEvent->IsA()==AliAODEvent::Class()){
      if(arrClustersDebug)
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersDebug->At(i));
      else
        clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
    }

    if(!clus) continue;
    if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC,fWeightJetJetMC,i)){
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
    if(!PhotonCandidate){ delete clus; delete tmpvec; continue;}

    fOutputLocalDebug << "--cluster" << i << "--" << endl;

    double etaCluster = clusterVector.Eta();
    double phiCluster = clusterVector.Phi();
    int nCellCluster = clus->GetNCells();
    for(int iCell=0;iCell<nCellCluster;iCell++){
      int nSupMod, nModule, nIphi, nIeta;
      ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetGeomEMCAL()->GetCellIndex(
            clus->GetCellAbsId(iCell),
            nSupMod,
            nModule,
            nIphi,
            nIeta);
      int nphi, neta;
      ((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->GetGeomEMCAL()->GetCellPhiEtaIndexInSModule(
            nSupMod,
            nModule,
            nIphi,
            nIeta,
            nphi,
            neta);
      fOutputLocalDebug << clus->GetCellAbsId(iCell) << " " << nSupMod << " " << nphi << " " << neta << " " << cells->GetCellAmplitude(clus->GetCellAbsId(iCell)) << endl;
    }
    fOutputLocalDebug << "phi " << phiCluster << endl;
    fOutputLocalDebug << "eta " << etaCluster << endl;

    delete clus;
    delete tmpvec;
  }
  fOutputLocalDebug << "---------------------------------------" << endl;
  fOutputLocalDebug.close();

  return;
}
// Function that searches iteratively for
// contributions of a mcparticle (and its daughters)
// in a cluster and returns the cluster id
int AliAnalysisTaskPi0v2Calo::CheckClustersForMCContribution(int mclabel, bool leading){
    int clusID = -1; // found contribution in ID cluster
    int nclus                       = 0;
    TClonesArray * arrClustersProcess = NULL;
    fNCurrentClusterBasic             = 0;
    AliMCParticle* particle = (AliMCParticle *)fMCEvent->GetTrack(mclabel);
    if(!fCorrTaskSetting.CompareTo("")){
      nclus = fInputEvent->GetNumberOfCaloClusters();
    } else {
      arrClustersProcess = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
      if(!arrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskPi0v2Calo! Check the correction framework settings!",fCorrTaskSetting.Data()));
      nclus = arrClustersProcess->GetEntries();
    }

    for(long i = 0; i < nclus; i++){
      AliVCluster* clus = NULL;
      if(fInputEvent->IsA()==AliESDEvent::Class()){
        if(arrClustersProcess)
          clus = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(i));
        else
          clus = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
      } else if(fInputEvent->IsA()==AliAODEvent::Class()){
        if(arrClustersProcess)
          clus = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(i));
        else
          clus = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
      }
      if(!clus) continue;

      if(!((AliCaloPhotonCuts*)fClusterCutArray->At(fiCut))->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, 1.,i)){
         delete clus;
         continue;
      }

      int *mclabelsCluster = clus->GetLabels();
      if (clus->GetNLabels() > 0)
      {
          if(!leading){
            for (int k = 0; k < (int)clus->GetNLabels(); k++)
            {
              if ((mclabelsCluster[k] == mclabel) ) clusID = i;
            }
          } else{
            if ((mclabelsCluster[0] == mclabel)) clusID = i;
          }
      }

      if(clusID != -1){ // do a check if its a conversion electron
          bool motherIsPhoton = kFALSE;
          int absPdg = TMath::Abs(particle->PdgCode());
          if(absPdg == 11){
              int mothLabel = particle->GetMother();
              if(mothLabel!= -1){
                 AliMCParticle* mother = (AliMCParticle *)fMCEvent->GetTrack(mothLabel);
                 if(mother->PdgCode() == 22){
                        motherIsPhoton = kTRUE;
                 }
                 if(motherIsPhoton){
                   // dont use conv electrons that have smaller than 0.5 of mother photon
                   if(particle->E()/mother->E() < 0.5) clusID = -1;
                 }
              }
          }

      }
      delete clus;
    }
    // If i did not find a contribution for this mc label, check iteratively for daughters
    // and break as soon as i found a contribution of a daughter in a cluster
    if(clusID == -1 && (particle->GetNDaughters()>0)){
       for (int daughter = particle->GetDaughterFirst(); daughter <= particle->GetDaughterLast(); daughter++)
       {
         clusID = CheckClustersForMCContribution(daughter, leading);
         if(clusID!=-1){
           break; // break when you found a contribution
         }
       }

    }
    return clusID;
}

// Function to load the calibration containers
void AliAnalysisTaskPi0v2Calo::LoadCalibrationContainers() {
  if (fPeriod.EqualTo("LHC15o")) {
    contQxncm = dynamic_cast<AliOADBContainer *>(fListVZEROCalib->FindObject(Form("fqxc%im", 2)));  // V0C Qx Mean
    contQyncm = dynamic_cast<AliOADBContainer *>(fListVZEROCalib->FindObject(Form("fqyc%im", 2)));  // V0C Qy Mean
    contQxnam = dynamic_cast<AliOADBContainer *>(fListVZEROCalib->FindObject(Form("fqxa%im", 2)));  // V0A Qx Mean
    contQynam = dynamic_cast<AliOADBContainer *>(fListVZEROCalib->FindObject(Form("fqya%im", 2)));  // V0A Qy Mean
    contMult = dynamic_cast<AliOADBContainer *>(fListVZEROCalib->FindObject("hMultV0BefCorPfpx"));
  } else if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    contQxncm = static_cast<AliOADBContainer *>(fVZEROCalibFile->GetObjectChecked("fqxc2m", "AliOADBContainer"));  // V0C Qx Mean
    contQyncm = static_cast<AliOADBContainer *>(fVZEROCalibFile->GetObjectChecked("fqyc2m", "AliOADBContainer"));  // V0C Qy Mean
    contQxnam = static_cast<AliOADBContainer *>(fVZEROCalibFile->GetObjectChecked("fqxa2m", "AliOADBContainer"));  // V0A Qx Mean
    contQynam = static_cast<AliOADBContainer *>(fVZEROCalibFile->GetObjectChecked("fqya2m", "AliOADBContainer"));  // V0A Qy Mean
  }
  if (!contQxncm || !contQyncm || !contQxnam || !contQynam) {
      AliInfo("One or more VZERO calibration containers not found.");
  }
}

bool AliAnalysisTaskPi0v2Calo::LoadVZEROCalibration(){
  // Connect to the alien grid if not already connected
  if (!gGrid) {
    TGrid::Connect("alien://");
  }


  // If VZERO calibration is enabled
  if (IsVZEROCalibOn) {
    // Load calibration for LHC15o period
    if (fPeriod.EqualTo("LHC15o")) {
      fVZEROCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC15o/VZEROCalibFile15o.root", "READ");
      if (fVZEROCalibFile) {
        fListVZEROCalib = dynamic_cast<TList *>(fVZEROCalibFile->Get("VZEROCalibList"));
        if (fListVZEROCalib) {
          LoadCalibrationContainers();
        } else {
          AliInfo("VZERO calibration list not found for LHC15o.");
          return false;
        }
      } else {
        AliInfo("VZERO calibration file not found for LHC15o.");
        return false;
      }
    }
    // Load calibration for LHC18q and LHC18r periods
    else if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
      fVZEROCalibFile = TFile::Open("alien:///alice/cern.ch/user/j/jwan/CalibFile/calibq2V0C18qrP3.root", "READ");
      if (fVZEROCalibFile) {
        LoadCalibrationContainers();
      } else {
        AliInfo("VZERO calibration file not found for LHC18q/r.");
        return false;
      }
    }
  }
  return true;
}

//________________________________________________________________________
bool AliAnalysisTaskPi0v2Calo::LoadCalibHistForThisRun() {
  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    // Check if calibration containers are available
    if (!contQxncm || !contQyncm || !contQxnam || !contQynam) {
      AliInfo(Form("Calibration containers not found for run %i", fRunNumber));
      return false;
    }

    if (fPeriod.EqualTo("LHC15o")) {
      // Load histograms for LHC15o period
      hMultV0 = static_cast<TH1D *>(contMult->GetObject(fRunNumber));
      hQx2mV0[0] = static_cast<TH1D *>(contQxncm->GetObject(fRunNumber));
      hQy2mV0[0] = static_cast<TH1D *>(contQyncm->GetObject(fRunNumber));
      hQx2mV0[1] = static_cast<TH1D *>(contQxnam->GetObject(fRunNumber));
      hQy2mV0[1] = static_cast<TH1D *>(contQynam->GetObject(fRunNumber));

      if (!hMultV0 || !hQx2mV0[0] || !hQy2mV0[0] || !hQx2mV0[1] || !hQy2mV0[1]) {
        AliInfo(Form("Histograms not found for LHC15o run %i", fRunNumber));
        return false;
      }
    } else if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
      // Load default objects for LHC18q/r period
      hQx2mV0[0] = static_cast<TH1D *>(contQxncm->GetDefaultObject(Form("hV0QxMeanCRun%d", fRunNumber)));
      hQy2mV0[0] = static_cast<TH1D *>(contQyncm->GetDefaultObject(Form("hV0QyMeanCRun%d", fRunNumber)));
      hQx2mV0[1] = static_cast<TH1D *>(contQxnam->GetDefaultObject(Form("hV0QxMeanARun%d", fRunNumber)));
      hQy2mV0[1] = static_cast<TH1D *>(contQynam->GetDefaultObject(Form("hV0QyMeanARun%d", fRunNumber)));

      if (!hQx2mV0[0] || !hQy2mV0[0] || !hQx2mV0[1] || !hQy2mV0[1]) {
        AliInfo(Form("Histograms not found for LHC18q/r run %i", fRunNumber));
        return false;
      }

      // Load channel correction weights
      fHCorrectV0ChWeghts = static_cast<TH2F *>(fVZEROCalibFile->GetObjectChecked(Form("hWgtV0ChannelsvsVzRun%d", fRunNumber), "TH2F"));
      if (!fHCorrectV0ChWeghts) {
        AliInfo(Form("Channel weights not found for run %i", fRunNumber));
        return false;
      }
    }
  } else {
    AliInfo(TString("Unknown period: " + fPeriod).Data());
    return false;
  }
  return true;
}


//________________________________________________________________________
bool AliAnalysisTaskPi0v2Calo::GetVZEROPlane() {
  // Initialize arrays
  double multV0Ch[64] = {0};
  double V0XMean[2] = {-999, -999}, V0YMean[2] = {-999, -999};
  double qxGE[2] = {0, 0}, qyGE[2] = {0, 0};
  double qxRC[2] = {0, 0}, qyRC[2] = {0, 0};
  double multRingGE[2] = {0, 0};
  double psi2GE[2] = {0, 0}, psi2RC[2] = {0, 0};

  // Get the primary vertex
  double vertex[3] = {0, 0, 0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  // Load V0 channel multiplicities if in LHC15o period
  if (fPeriod.EqualTo("LHC15o")) {
    for (int iCh = 0; iCh < 64; ++iCh) {
      multV0Ch[iCh] = hMultV0->GetBinContent(iCh + 1);
    }
  }

  // Get centrality and check valid range
  int iCentSPD = static_cast<int>(fCent);
  if (iCentSPD >= 90) return false;

  // Load Qx/Qy mean values for V0C and V0A
  for (int i = 0; i < 2; ++i) {
    V0XMean[i] = hQx2mV0[i]->GetBinContent(iCentSPD + 1);
    V0YMean[i] = hQy2mV0[i]->GetBinContent(iCentSPD + 1);
  }

  // Loop over V0 channels and calculate q-vectors for V0C (0-31) and V0A (32-63)
  AliAODVZERO *aodV0 = static_cast<AliAODVZERO *>(fInputEvent->GetVZEROData());

  for (int iCh = 0; iCh < 64; ++iCh) {
    double phi = TMath::Pi() / 8. + TMath::Pi() / 4. * (iCh % 8);
    double multCh = aodV0->GetMultiplicity(iCh);
    double multChGE = -1;

    if (fPeriod.EqualTo("LHC15o")) {
      int binGroup = (iCh < 32) ? iCh / 8 * 8 : (iCh - 32) / 8 * 8 + 32;
      multChGE = multCh / multV0Ch[iCh] * multV0Ch[binGroup];
    } else if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
      int ibinV0 = fHCorrectV0ChWeghts->FindBin(vertex[2], iCh);
      multChGE = multCh * fHCorrectV0ChWeghts->GetBinContent(ibinV0);
    }

    if (multChGE < 0) continue;

    int ringIndex = (iCh < 32) ? 0 : 1;  // 0 for V0C, 1 for V0A
    qxGE[ringIndex] += multChGE * TMath::Cos(2 * phi);
    qyGE[ringIndex] += multChGE * TMath::Sin(2 * phi);
    multRingGE[ringIndex] += multChGE;
  }

  // Check for valid q-vectors
  if (multRingGE[0] < 1.e-6 || multRingGE[1] < 1.e-6) return false;

  // Calculate VZERO GE Plane angles
  for (int i = 0; i < 2; ++i) {
    psi2GE[i] = GetEventPlane(qxGE[i], qyGE[i], 2.);
    if (TMath::IsNaN(psi2GE[i])) return false;
  }

  // Recenter VZERO Q-vectors
  for (int i = 0; i < 2; ++i) {
    double qxMean = V0XMean[i], qyMean = V0YMean[i];
    if (TMath::IsNaN(qxMean) || TMath::IsNaN(qyMean) || qxMean < -900 || qyMean < -900) continue;

    qxRC[i] = qxGE[i] - qxMean;
    qyRC[i] = qyGE[i] - qyMean;
    psi2RC[i] = GetEventPlane(qxRC[i], qyRC[i], 2.);
    if (TMath::IsNaN(psi2RC[i])) return false;
  }

  // Quality Assurance (QA) for VZERO if enabled
  if (IsQAVZERO) {
    for (int i = 0; i < 2; ++i) {
      // Profiles for V0C (0) and V0A (1)
      fProfileV0CQxCentGE[fiCut]->Fill(fCent, qxGE[i]);
      fProfileV0CQyCentGE[fiCut]->Fill(fCent, qyGE[i]);
      fProfileV0CQxVtxGE[fiCut]->Fill(vertex[2], qxGE[i]);
      fProfileV0CQyVtxGE[fiCut]->Fill(vertex[2], qyGE[i]);
      fHist2CalibPsi2V0CCentGE[fiCut]->Fill(fCent, psi2GE[i]);

      fProfileV0CQxCentRC[fiCut]->Fill(fCent, qxRC[i]);
      fProfileV0CQyCentRC[fiCut]->Fill(fCent, qyRC[i]);
      fProfileV0CQxVtxRC[fiCut]->Fill(vertex[2], qxRC[i]);
      fProfileV0CQyVtxRC[fiCut]->Fill(vertex[2], qyRC[i]);
      fHist2CalibPsi2V0CCentRC[fiCut]->Fill(fCent, psi2RC[i]);
    }
  }

  // Set Psi2 values for V0C and V0A
  fPsi2V0C = psi2RC[0];
  fPsi2V0A = psi2RC[1];
  return true;
}

//________________________________________________________________________
double AliAnalysisTaskPi0v2Calo::GetEventPlane(double qx, double qy, double harmonic){
  double psi = (1. / harmonic) * TMath::ATan2(qy, qx);
  if (psi < 0)
    return psi += TMath::TwoPi() / harmonic;
  else
    return psi;
}

void AliAnalysisTaskPi0v2Calo::FillFlowHistograms(AliAODConversionMother *pi0cand){
  fHistoMotherInvMassPtV0CCos2phi[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), cos(2 * (pi0cand->Phi() - fPsi2V0C)));
  fHistoMotherInvMassPtV0ACos2phi[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), cos(2 * (pi0cand->Phi() - fPsi2V0A)));
  double dphiV0A = TVector2::Phi_0_2pi(pi0cand->Phi() - fPsi2V0A);
  if (dphiV0A > TMath::Pi()) {
    dphiV0A -= TMath::Pi();
  }
  double dphiV0C = TVector2::Phi_0_2pi(pi0cand->Phi() - fPsi2V0C);
  if (dphiV0C > TMath::Pi()) {
    dphiV0C -= TMath::Pi();
  }
  if (!fUseInOutPlane) {
    fHistoMotherInvMassPtPhiV0C[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), dphiV0C, pi0cand->GetWeight());
    fHistoMotherInvMassPtPhiV0A[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), dphiV0A, pi0cand->GetWeight());
  } else
  {
    // V0A Plane
    if (dphiV0A < 0.25 * TMath::Pi() || (dphiV0A > 0.75 * TMath::Pi() && dphiV0A < 1.25 * TMath::Pi()) || dphiV0A > 1.75 * TMath::Pi()) {
      fHistoMotherInvMassPtV0AInPlane[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), pi0cand->GetWeight());
    } else {
      fHistoMotherInvMassPtV0AOutPlane[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), pi0cand->GetWeight());
    }
    // V0C Plane
    if (dphiV0C < 0.25 * TMath::Pi() || (dphiV0C > 0.75 * TMath::Pi() && dphiV0C < 1.25 * TMath::Pi()) || dphiV0C > 1.75 * TMath::Pi()) {
      fHistoMotherInvMassPtV0CInPlane[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), pi0cand->GetWeight());
    } else {
      fHistoMotherInvMassPtV0COutPlane[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), pi0cand->GetWeight());
    }
  }
}

void AliAnalysisTaskPi0v2Calo::FillFlowBackHistograms(AliAODConversionMother *pi0cand, double weight){
  fHistoMotherBackInvMassPtV0CCos2phi[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), cos(2 * (pi0cand->Phi() - fPsi2V0C)), weight);
  fHistoMotherBackInvMassPtV0ACos2phi[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), cos(2 * (pi0cand->Phi() - fPsi2V0A)), weight);
  double dphiV0A = TVector2::Phi_0_2pi(pi0cand->Phi() - fPsi2V0A);
  if (dphiV0A > TMath::Pi()) {
    dphiV0A -= TMath::Pi();
  }
  double dphiV0C = TVector2::Phi_0_2pi(pi0cand->Phi() - fPsi2V0C);
  if (dphiV0C > TMath::Pi()) {
    dphiV0C -= TMath::Pi();
  }
  if (!fUseInOutPlane) {
    fHistoMotherBackInvMassPtPhiV0C[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), dphiV0C, pi0cand->GetWeight() * weight);
    fHistoMotherBackInvMassPtPhiV0A[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), dphiV0A, pi0cand->GetWeight() * weight);
  } else
  {
    // V0A Plane
    if (dphiV0A < 0.25 * TMath::Pi() || (dphiV0A > 0.75 * TMath::Pi() && dphiV0A < 1.25 * TMath::Pi()) || dphiV0A > 1.75 * TMath::Pi()) {
      fHistoMotherBackInvMassPtV0AInPlane[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), pi0cand->GetWeight() * weight);
    } else {
      fHistoMotherBackInvMassPtV0AOutPlane[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), pi0cand->GetWeight() * weight);
    }
    // V0C Plane
    if (dphiV0C < 0.25 * TMath::Pi() || (dphiV0C > 0.75 * TMath::Pi() && dphiV0C < 1.25 * TMath::Pi()) || dphiV0C > 1.75 * TMath::Pi()) {
      fHistoMotherBackInvMassPtV0CInPlane[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), pi0cand->GetWeight() * weight);
    } else {
      fHistoMotherBackInvMassPtV0COutPlane[fiCut]->Fill(pi0cand->M(), pi0cand->Pt(), pi0cand->GetWeight() * weight);
    }
  }
}
