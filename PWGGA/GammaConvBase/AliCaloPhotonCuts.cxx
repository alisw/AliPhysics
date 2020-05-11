/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Friederike Bock, Daniel Muehlheim                             *
 * Version 1.0                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class handling all kinds of selection cuts for
// Photon from EMCAL clusters
//---------------------------------------------
////////////////////////////////////////////////

#include "AliCaloPhotonCuts.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODHandler.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "AliMCEvent.h"
#include "AliAODConversionMother.h"
#include "AliAODConversionPhoton.h"
#include "TObjString.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "TList.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliV0ReaderV1.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliPicoTrack.h"
#include "AliPHOSGeoUtils.h"
#include "AliTrackerBase.h"
#include "AliVCaloCells.h"
#include "AliVCluster.h"
#include "AliEmcalCorrectionTask.h"
#include "AliEmcalCorrectionComponent.h"
#include "AliTender.h"
#include "AliTenderSupply.h"
#include "AliEMCALTenderSupply.h"
#include "AliEmcalTenderTask.h"
#include "AliPHOSTenderSupply.h"
#include "AliPHOSTenderTask.h"
#include "AliOADBContainer.h"
#include "AliESDtrackCuts.h"
#include "AliCaloTrackMatcher.h"
#include "AliPhotonIsolation.h"
#include <memory>
#include <vector>

class iostream;

using namespace std;

ClassImp(AliCaloPhotonCuts)


const char* AliCaloPhotonCuts::fgkCutNames[AliCaloPhotonCuts::kNCuts] = {
  "ClusterType",          //0    0: all,    1: EMCAL,   2: PHOS
  "EtaMin",               //1    0: -10,    1: -0.6687, 2: -0,5, 3: -2
  "EtaMax",               //2    0: 10,     1: 0.66465, 2: 0.5,  3: 2
  "PhiMin",               //3    0: -10000, 1: 1.39626
  "PhiMax",               //4    0: 10000, 1: 3.125
  "NonLinearity1"         //5
  "NonLinearity2"         //6
  "DistanceToBadChannel", //7    0: 0,      1: 5
  "Timing",               //8    0: no cut
  "TrackMatching",        //9    0: 0,      1: 5
  "ExoticCluster",        //10   0: no cut
  "MinEnergy",            //11   0: no cut, 1: 0.05,    2: 0.1,  3: 0.15, 4: 0.2, 5: 0.3, 6: 0.5, 7: 0.75, 8: 1, 9: 1.25 (all GeV)
  "MinNCells",            //12   0: no cut, 1: 1,       2: 2,    3: 3,    4: 4,   5: 5,   6: 6
  "MinM02",               //13
  "MaxM02",               //14
  "MinMaxM20",            //15
  "RecConv",              //16
  "MaximumDispersion",    //17
  "NLM"                   //18
};


//________________________________________________________________________
AliCaloPhotonCuts::AliCaloPhotonCuts(Int_t isMC, const char *name,const char *title) :
  AliAnalysisCuts(name,title),
  fHistograms(NULL),
  fHistExtQA(NULL),
  fCaloTrackMatcher(NULL),
  fCaloIsolation(NULL),
  fGeomEMCAL(NULL),
  fEMCALRecUtils(NULL),
  fEMCALInitialized(kFALSE),
  fGeomPHOS(NULL),
  fPHOSGeoUtils(NULL),
  fPHOSInitialized(kFALSE),
  fPHOSCurrentRun(-1),
  fEMCALBadChannelsMap(NULL),
  fEMCALBadChannelsMap1D(NULL),
  fPHOSBadChannelsMap(NULL),
  fBadChannels(NULL),
  fNMaxEMCalModules(12),
  fNMaxPHOSModules(5),
  fHistoModifyAcc(NULL),
  fAODMCTrackArray(NULL),
  fDoLightOutput(0),
  fIsMC(0),
  fIsCurrentClusterAcceptedBeforeTM(kFALSE),
  fUseEtaPhiMapForBackCand(kFALSE),
  fV0ReaderName("V0ReaderV1"),
  fCorrTaskSetting(""),
  fCaloTrackMatcherName("CaloTrackMatcher_1_0"),
  fCaloIsolationName("PhotonIsolation"),
  fPeriodName(""),
  fCurrentMC(kNoMC),
  fClusterType(0),
  fIsolationRadius(0.4),
  fMomPercentage(0.1),
  fUsePhotonIsolation(kFALSE),
  fMinEtaCut(-10),
  fMinEtaInnerEdge(0),
  fMaxEtaCut(10),
  fMaxEtaInnerEdge(0),
  fUseEtaCut(0),
  fMinPhiCut(-10000),
  fMaxPhiCut(-10000),
  fMinPhiCutDMC(-10000),
  fMaxPhiCutDMC(-10000),
  fUsePhiCut(0),
  fMinDistanceToBadChannel(0),
  fUseDistanceToBadChannel(0),
  fMaxTimeDiff(10e10),
  fMinTimeDiff(-10e10),
  fMaxTimeDiffHighPt(10e10),
  fMinTimeDiffHighPt(-10e10),
  fUseTimeDiff(0),
  fMaxDistTrackToClusterEta(0),
  fMinDistTrackToClusterPhi(0),
  fMaxDistTrackToClusterPhi(0),
  fUseDistTrackToCluster(0),
  fUsePtDepTrackToCluster(0),
  fFuncPtDepEta(0),
  fFuncPtDepPhi(0),
  fRandom(0),
  fUseTimingEfficiencyMCSimCluster(0),
  fFuncTimingEfficiencyMCSimCluster(0),
  fFuncTimingEfficiencyMCSimClusterHighPt(0),
  fTimingEfficiencyMCSimClusterLowPtEnd(4.),
  fTimingEfficiencyMCSimClusterHighPtStart(6.),
  fFuncNCellCutEfficiencyEMCal(0),
  fMinTMDistSigma(10),
  fUseEOverPVetoTM(0),
  fEOverPMax(0.),
  fUseTMMIPsubtraction(0),
  fUseElectronClusterCalibration(0),
  fExtendedMatchAndQA(0),
  fExoticEnergyFracCluster(0),
  fExoticMinEnergyTCard(0),
  fExoticMinEnergyCell(1),
  fUseExoticCluster(0),
  fDoExoticsQA(kFALSE),
  fMinEnergy(0),
  fDoFlatEnergySubtraction(0),
  fSeedEnergy(0.1),
  fLocMaxCutEDiff(0.03),
  fUseMinEnergy(0),
  fMinNCells(0),
  fUseNCells(0),
  fMaxM02(1000),
  fMinM02(0),
  fUseM02(0),
  fMaxM02CutNr(0),
  fMinM02CutNr(0),
  fMaxM20(1000),
  fMinM20(0),
  fUseM20(0),
  fMaxMGGRecConv(0),
  fUseRecConv(0),
  fMaxDispersion(1000),
  fUseDispersion(0),
  fMinNLM(0),
  fMaxNLM(1000),
  fUseNLM(0),
  fNonLinearity1(0),
  fNonLinearity2(0),
  fSwitchNonLinearity(0),
  fUseNonLinearity(kFALSE),
  fIsPureCalo(0),
  fNactiveEmcalCells(0),
  fDoSecondaryTrackMatching(kFALSE),
  fVectorMatchedClusterIDs(0),
  fCutString(NULL),
  fCutStringRead(""),
  fHistCutIndex(NULL),
  fHistAcceptanceCuts(NULL),
  fHistClusterIdentificationCuts(NULL),
  fHistClusterEtavsPhiBeforeAcc(NULL),
  fHistClusterEtavsPhiAfterAcc(NULL),
  fHistClusterEtavsPhiAfterQA(NULL),
  fHistClusterEtavsPhiBackground(NULL),
  fHistClusterTimevsEBeforeQA(NULL),
  fHistClusterTimevsEAfterQA(NULL),
  fHistEnergyOfClusterBeforeNL(NULL),
  fHistEnergyOfClusterAfterNL(NULL),
  fHistEnergyOfClusterBeforeQA(NULL),
  fHistEnergyOfClusterAfterQA(NULL),
  fHistNCellsBeforeQA(NULL),
  fHistNCellsAfterQA(NULL),
  fHistM02BeforeQA(NULL),
  fHistM02AfterQA(NULL),
  fHistM20BeforeQA(NULL),
  fHistM20AfterQA(NULL),
  fHistDispersionBeforeQA(NULL),
  fHistDispersionAfterQA(NULL),
  fHistNLMBeforeQA(NULL),
  fHistNLMAfterQA(NULL),
  fHistNLMVsNCellsAfterQA(NULL),
  fHistNLMVsEAfterQA(NULL),
//   fHistNLMAvsNLMBBeforeQA(NULL),
  fHistClusterEnergyvsMod(NULL),
  fHistNCellsBigger100MeVvsMod(NULL),
  fHistNCellsBigger1500MeVvsMod(NULL),
  fHistEnergyOfModvsMod(NULL),
  fHistClusterEnergyvsNCellsBeforeQA(NULL),
  fHistClusterEnergyvsNCellsAfterQA(NULL),
  fHistCellEnergyvsCellID(NULL),
  fHistCellTimevsCellID(NULL),
  fHistClusterEM02BeforeQA(NULL),
  fHistClusterEM02AfterQA(NULL),
  fHistClusterIncludedCellsBeforeQA(NULL),
  fHistClusterIncludedCellsAfterQA(NULL),
  fHistClusterEnergyFracCellsBeforeQA(NULL),
  fHistClusterEnergyFracCellsAfterQA(NULL),
  fHistClusterIncludedCellsTimingAfterQA(NULL),
  fHistClusterIncludedCellsTimingEnergyAfterQA(NULL),
  fHistClusterDistanceInTimeCut(NULL),
  fHistClusterDistanceOutTimeCut(NULL),
  fHistClusterDistance1DInTimeCut(NULL),
  fHistClusterRBeforeQA(NULL),
  fHistClusterRAfterQA(NULL),
  fHistClusterdEtadPhiBeforeQA(NULL),
  fHistClusterdEtadPhiAfterQA(NULL),
  fHistDistanceTrackToClusterBeforeQA(NULL),
  fHistDistanceTrackToClusterAfterQA(NULL),
  fHistClusterdEtadPhiPosTracksBeforeQA(NULL),
  fHistClusterdEtadPhiNegTracksBeforeQA(NULL),
  fHistClusterdEtadPhiPosTracksAfterQA(NULL),
  fHistClusterdEtadPhiNegTracksAfterQA(NULL),
  fHistClusterdEtadPhiPosTracksP_000_075BeforeQA(NULL),
  fHistClusterdEtadPhiPosTracksP_075_125BeforeQA(NULL),
  fHistClusterdEtadPhiPosTracksP_125_999BeforeQA(NULL),
  fHistClusterdEtadPhiNegTracksP_000_075BeforeQA(NULL),
  fHistClusterdEtadPhiNegTracksP_075_125BeforeQA(NULL),
  fHistClusterdEtadPhiNegTracksP_125_999BeforeQA(NULL),
  fHistClusterdEtadPtBeforeQA(NULL),
  fHistClusterdEtadPtAfterQA(NULL),
  fHistClusterdEtadPtTrueMatched(NULL),
  fHistClusterdPhidPtPosTracksBeforeQA(NULL),
  fHistClusterdPhidPtNegTracksBeforeQA(NULL),
  fHistClusterdPhidPtAfterQA(NULL),
  fHistClusterdPhidPtPosTracksTrueMatched(NULL),
  fHistClusterdPhidPtNegTracksTrueMatched(NULL),
  fHistClusterM02M20BeforeQA(NULL),
  fHistClusterM02M20AfterQA(NULL),
  fHistClusterEtavsPhiExotics(NULL),
  fHistClusterEM02Exotics(NULL),
  fHistClusterEnergyvsNCellsExotics(NULL),
  fHistClusterEEstarExotics(NULL),
  fHistClusterTMEffiInput(NULL),
  fHistClusterTrueElecEtaPhiBeforeTM_30_00(NULL),
  fHistClusterTrueElecEtaPhiAfterTM_30_00(NULL),
  fHistClusterEvsTrackECharged(NULL),
  fHistClusterEvsTrackEChargedLead(NULL),
  fHistClusterEvsTrackENeutral(NULL),
  fHistClusterEvsTrackENeutralSubCharged(NULL),
  fHistClusterEvsTrackEGamma(NULL),
  fHistClusterEvsTrackEGammaSubCharged(NULL),
  fHistClusterEvsTrackEConv(NULL),
  fHistClusterENMatchesNeutral(NULL),
  fHistClusterENMatchesCharged(NULL),
  fHistClusterEvsTrackEPrimaryButNoElec(NULL),
  fHistClusterEvsTrackSumEPrimaryButNoElec(NULL),
  fHistClusETruePi0_BeforeTM(NULL),
  fHistClusETruePi0_Matched(NULL),
  fHistMatchedTrackPClusE(NULL),
  fHistMatchedTrackPClusEAfterEOverPVeto(NULL),
  fHistMatchedTrackPClusETruePi0Clus(NULL),
  fHistElectronPositronClusterMatch(NULL),
  fHistElectronPositronClusterMatchSub(NULL),
  fHistElectronClusterMatch(NULL),
  fHistPositronClusterMatch(NULL),
  fHistTrueElectronPositronClusterMatch(NULL),
  fHistTrueNoElectronPositronClusterMatch(NULL),
  fHistElectronClusterMatchTruePID(NULL),
  fHistInvMassDiCluster(NULL),
  fHistInvMassConvFlagging(NULL),
  fNMaxDCalModules(8),
  fgkDCALCols(32),
  fIsAcceptedForBasic(kFALSE)
{
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=0;}
  fCutString=new TObjString((GetCutNumber()).Data());

  fIsMC = isMC;
}

//________________________________________________________________________
AliCaloPhotonCuts::AliCaloPhotonCuts(const AliCaloPhotonCuts &ref) :
  AliAnalysisCuts(ref),
  fHistograms(NULL),
  fHistExtQA(NULL),
  fCaloTrackMatcher(NULL),
  fCaloIsolation(NULL),
  fGeomEMCAL(NULL),
  fEMCALRecUtils(NULL),
  fEMCALInitialized(kFALSE),
  fGeomPHOS(NULL),
  fPHOSGeoUtils(NULL),
  fPHOSInitialized(kFALSE),
  fPHOSCurrentRun(-1),
  fEMCALBadChannelsMap(NULL),
  fEMCALBadChannelsMap1D(NULL),
  fPHOSBadChannelsMap(NULL),
  fBadChannels(NULL),
  fNMaxEMCalModules(ref.fNMaxEMCalModules),
  fNMaxPHOSModules(ref.fNMaxPHOSModules),
  fHistoModifyAcc(NULL),
  fAODMCTrackArray(NULL),
  fDoLightOutput(ref.fDoLightOutput),
  fIsMC(ref.fIsMC),
  fIsCurrentClusterAcceptedBeforeTM(kFALSE),
  fUseEtaPhiMapForBackCand(ref.fUseEtaPhiMapForBackCand),
  fV0ReaderName(ref.fV0ReaderName),
  fCorrTaskSetting(ref.fCorrTaskSetting),
  fCaloTrackMatcherName(ref.fCaloTrackMatcherName),
  fCaloIsolationName("PhotonIsolation"),
  fPeriodName(ref.fPeriodName),
  fCurrentMC(ref.fCurrentMC),
  fClusterType(ref.fClusterType),
  fIsolationRadius(0.4),
  fMomPercentage(0.1),
  fUsePhotonIsolation(kFALSE),
  fMinEtaCut(ref.fMinEtaCut),
  fMinEtaInnerEdge(ref.fMinEtaInnerEdge),
  fMaxEtaCut(ref.fMaxEtaCut),
  fMaxEtaInnerEdge(ref.fMaxEtaInnerEdge),
  fUseEtaCut(ref.fUseEtaCut),
  fMinPhiCut(ref.fMinPhiCut),
  fMaxPhiCut(ref.fMaxPhiCut),
  fMinPhiCutDMC(ref.fMinPhiCutDMC),
  fMaxPhiCutDMC(ref.fMaxPhiCutDMC),
  fUsePhiCut(ref.fUsePhiCut),
  fMinDistanceToBadChannel(ref.fMinDistanceToBadChannel),
  fUseDistanceToBadChannel(ref.fUseDistanceToBadChannel),
  fMaxTimeDiff(ref.fMaxTimeDiff),
  fMinTimeDiff(ref.fMinTimeDiff),
  fMaxTimeDiffHighPt(ref.fMaxTimeDiffHighPt),
  fMinTimeDiffHighPt(ref.fMinTimeDiffHighPt),
  fUseTimeDiff(ref.fUseTimeDiff),
  fMaxDistTrackToClusterEta(ref.fMaxDistTrackToClusterEta),
  fMinDistTrackToClusterPhi(ref.fMinDistTrackToClusterPhi),
  fMaxDistTrackToClusterPhi(ref.fMaxDistTrackToClusterPhi),
  fUseDistTrackToCluster(ref.fUseDistTrackToCluster),
  fUsePtDepTrackToCluster(ref.fUsePtDepTrackToCluster),
  fFuncPtDepEta(ref.fFuncPtDepEta),
  fFuncPtDepPhi(ref.fFuncPtDepPhi),
  fRandom(ref.fRandom),
  fUseTimingEfficiencyMCSimCluster(ref.fUseTimingEfficiencyMCSimCluster),
  fFuncTimingEfficiencyMCSimCluster(ref.fFuncTimingEfficiencyMCSimCluster),
  fFuncTimingEfficiencyMCSimClusterHighPt(ref.fFuncTimingEfficiencyMCSimClusterHighPt),
  fTimingEfficiencyMCSimClusterLowPtEnd(ref.fTimingEfficiencyMCSimClusterLowPtEnd),
  fTimingEfficiencyMCSimClusterHighPtStart(ref.fTimingEfficiencyMCSimClusterHighPtStart),
  fFuncNCellCutEfficiencyEMCal(ref.fFuncNCellCutEfficiencyEMCal),
  fMinTMDistSigma(ref.fMinTMDistSigma),
  fUseEOverPVetoTM(ref.fUseEOverPVetoTM),
  fEOverPMax(ref.fEOverPMax),
  fUseTMMIPsubtraction(ref.fUseTMMIPsubtraction),
  fUseElectronClusterCalibration(ref.fUseElectronClusterCalibration),
  fExtendedMatchAndQA(ref.fExtendedMatchAndQA),
  fExoticEnergyFracCluster(ref.fExoticEnergyFracCluster),
  fExoticMinEnergyTCard(ref.fExoticMinEnergyTCard),
  fExoticMinEnergyCell(ref.fExoticMinEnergyCell),
  fUseExoticCluster(ref.fUseExoticCluster),
  fDoExoticsQA(ref.fDoExoticsQA),
  fMinEnergy(ref.fMinEnergy),
  fDoFlatEnergySubtraction(ref.fDoFlatEnergySubtraction),
  fSeedEnergy(ref.fSeedEnergy),
  fLocMaxCutEDiff(ref.fLocMaxCutEDiff),
  fUseMinEnergy(ref.fUseMinEnergy),
  fMinNCells(ref.fMinNCells),
  fUseNCells(ref.fUseNCells),
  fMaxM02(ref.fMaxM02),
  fMinM02(ref.fMinM02),
  fUseM02(ref.fUseM02),
  fMaxM02CutNr(ref.fMaxM02CutNr),
  fMinM02CutNr(ref.fMinM02CutNr),
  fMaxM20(ref.fMaxM20),
  fMinM20(ref.fMinM20),
  fUseM20(ref.fUseM20),
  fMaxMGGRecConv(ref.fMaxMGGRecConv),
  fUseRecConv(ref.fUseRecConv),
  fMaxDispersion(ref.fMaxDispersion),
  fUseDispersion(ref.fUseDispersion),
  fMinNLM(ref.fMinNLM),
  fMaxNLM(ref.fMaxNLM),
  fUseNLM(ref.fUseNLM),
  fNonLinearity1(ref.fNonLinearity1),
  fNonLinearity2(ref.fNonLinearity2),
  fSwitchNonLinearity(ref.fSwitchNonLinearity),
  fUseNonLinearity(ref.fUseNonLinearity),
  fIsPureCalo(ref.fIsPureCalo),
  fNactiveEmcalCells(ref.fNactiveEmcalCells),
  fDoSecondaryTrackMatching(ref.fDoSecondaryTrackMatching),
  fVectorMatchedClusterIDs(0),
  fCutString(NULL),
  fCutStringRead(""),
  fHistCutIndex(NULL),
  fHistAcceptanceCuts(NULL),
  fHistClusterIdentificationCuts(NULL),
  fHistClusterEtavsPhiBeforeAcc(NULL),
  fHistClusterEtavsPhiAfterAcc(NULL),
  fHistClusterEtavsPhiAfterQA(NULL),
  fHistClusterEtavsPhiBackground(NULL),
  fHistClusterTimevsEBeforeQA(NULL),
  fHistClusterTimevsEAfterQA(NULL),
  fHistEnergyOfClusterBeforeNL(NULL),
  fHistEnergyOfClusterAfterNL(NULL),
  fHistEnergyOfClusterBeforeQA(NULL),
  fHistEnergyOfClusterAfterQA(NULL),
  fHistNCellsBeforeQA(NULL),
  fHistNCellsAfterQA(NULL),
  fHistM02BeforeQA(NULL),
  fHistM02AfterQA(NULL),
  fHistM20BeforeQA(NULL),
  fHistM20AfterQA(NULL),
  fHistDispersionBeforeQA(NULL),
  fHistDispersionAfterQA(NULL),
  fHistNLMBeforeQA(NULL),
  fHistNLMAfterQA(NULL),
//   fHistNLMAvsNLMBBeforeQA(NULL),
  fHistNLMVsNCellsAfterQA(NULL),
  fHistNLMVsEAfterQA(NULL),
  fHistClusterEnergyvsMod(NULL),
  fHistNCellsBigger100MeVvsMod(NULL),
  fHistNCellsBigger1500MeVvsMod(NULL),
  fHistEnergyOfModvsMod(NULL),
  fHistClusterEnergyvsNCellsBeforeQA(NULL),
  fHistClusterEnergyvsNCellsAfterQA(NULL),
  fHistCellEnergyvsCellID(NULL),
  fHistCellTimevsCellID(NULL),
  fHistClusterEM02BeforeQA(NULL),
  fHistClusterEM02AfterQA(NULL),
  fHistClusterIncludedCellsBeforeQA(NULL),
  fHistClusterIncludedCellsAfterQA(NULL),
  fHistClusterEnergyFracCellsBeforeQA(NULL),
  fHistClusterEnergyFracCellsAfterQA(NULL),
  fHistClusterIncludedCellsTimingAfterQA(NULL),
  fHistClusterIncludedCellsTimingEnergyAfterQA(NULL),
  fHistClusterDistanceInTimeCut(NULL),
  fHistClusterDistanceOutTimeCut(NULL),
  fHistClusterDistance1DInTimeCut(NULL),
  fHistClusterRBeforeQA(NULL),
  fHistClusterRAfterQA(NULL),
  fHistClusterdEtadPhiBeforeQA(NULL),
  fHistClusterdEtadPhiAfterQA(NULL),
  fHistDistanceTrackToClusterBeforeQA(NULL),
  fHistDistanceTrackToClusterAfterQA(NULL),
  fHistClusterdEtadPhiPosTracksBeforeQA(NULL),
  fHistClusterdEtadPhiNegTracksBeforeQA(NULL),
  fHistClusterdEtadPhiPosTracksAfterQA(NULL),
  fHistClusterdEtadPhiNegTracksAfterQA(NULL),
  fHistClusterdEtadPhiPosTracksP_000_075BeforeQA(NULL),
  fHistClusterdEtadPhiPosTracksP_075_125BeforeQA(NULL),
  fHistClusterdEtadPhiPosTracksP_125_999BeforeQA(NULL),
  fHistClusterdEtadPhiNegTracksP_000_075BeforeQA(NULL),
  fHistClusterdEtadPhiNegTracksP_075_125BeforeQA(NULL),
  fHistClusterdEtadPhiNegTracksP_125_999BeforeQA(NULL),
  fHistClusterdEtadPtBeforeQA(NULL),
  fHistClusterdEtadPtAfterQA(NULL),
  fHistClusterdEtadPtTrueMatched(NULL),
  fHistClusterdPhidPtPosTracksBeforeQA(NULL),
  fHistClusterdPhidPtNegTracksBeforeQA(NULL),
  fHistClusterdPhidPtAfterQA(NULL),
  fHistClusterdPhidPtPosTracksTrueMatched(NULL),
  fHistClusterdPhidPtNegTracksTrueMatched(NULL),
  fHistClusterM02M20BeforeQA(NULL),
  fHistClusterM02M20AfterQA(NULL),
  fHistClusterEtavsPhiExotics(NULL),
  fHistClusterEM02Exotics(NULL),
  fHistClusterEnergyvsNCellsExotics(NULL),
  fHistClusterEEstarExotics(NULL),
  fHistClusterTMEffiInput(NULL),
  fHistClusterTrueElecEtaPhiBeforeTM_30_00(NULL),
  fHistClusterTrueElecEtaPhiAfterTM_30_00(NULL),
  fHistClusterEvsTrackECharged(NULL),
  fHistClusterEvsTrackEChargedLead(NULL),
  fHistClusterEvsTrackENeutral(NULL),
  fHistClusterEvsTrackENeutralSubCharged(NULL),
  fHistClusterEvsTrackEGamma(NULL),
  fHistClusterEvsTrackEGammaSubCharged(NULL),
  fHistClusterEvsTrackEConv(NULL),
  fHistClusterENMatchesNeutral(NULL),
  fHistClusterENMatchesCharged(NULL),
  fHistClusterEvsTrackEPrimaryButNoElec(NULL),
  fHistClusterEvsTrackSumEPrimaryButNoElec(NULL),
  fHistClusETruePi0_BeforeTM(NULL),
  fHistClusETruePi0_Matched(NULL),
  fHistMatchedTrackPClusE(NULL),
  fHistMatchedTrackPClusEAfterEOverPVeto(NULL),
  fHistMatchedTrackPClusETruePi0Clus(NULL),
  fHistElectronPositronClusterMatch(NULL),
  fHistElectronPositronClusterMatchSub(NULL),
  fHistElectronClusterMatch(NULL),
  fHistPositronClusterMatch(NULL),
  fHistTrueElectronPositronClusterMatch(NULL),
  fHistTrueNoElectronPositronClusterMatch(NULL),
  fHistElectronClusterMatchTruePID(NULL),
  fHistInvMassDiCluster(NULL),
  fHistInvMassConvFlagging(NULL),
  fNMaxDCalModules(ref.fNMaxDCalModules),
  fgkDCALCols(ref.fgkDCALCols),
  fIsAcceptedForBasic(ref.fIsAcceptedForBasic)
{
  // Copy Constructor
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=ref.fCuts[jj];}
  fCutString=new TObjString((GetCutNumber()).Data());

}


//________________________________________________________________________
AliCaloPhotonCuts::~AliCaloPhotonCuts() {
  // Destructor
  if(fCutString != NULL){
    delete fCutString;
    fCutString = NULL;
  }

  if(fPHOSBadChannelsMap != NULL){
    delete[] fPHOSBadChannelsMap;
    fPHOSBadChannelsMap = NULL;
  }

  if(fAODMCTrackArray){
    delete[] fAODMCTrackArray;
    fAODMCTrackArray = 0x0;
  }

  if(fFuncPtDepEta) delete fFuncPtDepEta;
  if(fFuncPtDepPhi) delete fFuncPtDepPhi;
  if(fFuncTimingEfficiencyMCSimCluster) delete fFuncTimingEfficiencyMCSimCluster;
  if(fFuncTimingEfficiencyMCSimClusterHighPt) delete fFuncTimingEfficiencyMCSimClusterHighPt;
  if(fFuncNCellCutEfficiencyEMCal) delete fFuncNCellCutEfficiencyEMCal;
}

//________________________________________________________________________
void AliCaloPhotonCuts::InitCutHistograms(TString name){

  // Initialize Cut Histograms for QA (only initialized and filled if function is called)
  TH1::AddDirectory(kFALSE);



  if (fDoLightOutput)
    fDoExoticsQA    = kFALSE;

  if(fHistograms != NULL){
    delete fHistograms;
    fHistograms=NULL;
  }
  if(fDoLightOutput==2) {
      AliInfo("Minimal output chosen");
      return;
  }

  if(fHistograms==NULL){
    fHistograms     = new TList();
    fHistograms->SetOwner(kTRUE);
    if(name=="")fHistograms->SetName(Form("CaloCuts_%s",GetCutNumber().Data()));
    else fHistograms->SetName(Form("%s_%s",name.Data(),GetCutNumber().Data()));
  }

  if(fHistExtQA != NULL){
    delete fHistExtQA;
    fHistExtQA=NULL;
  }
  if(fHistExtQA==NULL){
    fHistExtQA      = new TList();
    fHistExtQA->SetOwner(kTRUE);
    if(name=="")fHistExtQA->SetName(Form("CaloExtQA_%s",GetCutNumber().Data()));
    else fHistExtQA->SetName(Form("%s_ExtQA_%s",name.Data(),GetCutNumber().Data()));
  }

  Double_t *arrClusEBinning        = new Double_t[400];
  Double_t *arrClusEBinningCoarse  = new Double_t[400];
  Double_t *arrClusEBinningOnlyHighPt  = new Double_t[400];
  Int_t nBinsClusterE              = 235;
  Int_t nBinsClusterEMod           = 249;
  Int_t nBinsClusterECell          = 119;
  Int_t nBinsClusterECellCoarse    = 109;
  Int_t nBinsClusterEOnlyHighPt    = 134;
  Double_t maxClusterE             = 400;
  for(Int_t i=0; i<nBinsClusterEMod+1;i++){
    if (i < 1) arrClusEBinning[i]          = 0.3*i;
    else if(i<55) arrClusEBinning[i]       = 0.3+0.05*(i-1);
    else if(i<105) arrClusEBinning[i]      = 3.+0.1*(i-55);
    else if(i<140) arrClusEBinning[i]      = 8.+0.2*(i-105);
    else if(i<170) arrClusEBinning[i]      = 15.+0.5*(i-140);
    else if(i<190) arrClusEBinning[i]      = 30.+1.0*(i-170);
    else if(i<215) arrClusEBinning[i]      = 50.+2.0*(i-190);
    else if(i<235) arrClusEBinning[i]      = 100.+5.0*(i-215);
    else if(i<245) arrClusEBinning[i]      = 200.+10.0*(i-235);
    else if(i<249) arrClusEBinning[i]      = 300.+25.0*(i-245);
    else arrClusEBinning[i]                = maxClusterE;
  }
  for(Int_t i=0; i<nBinsClusterECell+1;i++){
    if(i<20) arrClusEBinningCoarse[i]             = 0.05*(i);
    else if(i<50) arrClusEBinningCoarse[i]        = 1.+0.1*(i-20);
    else if(i<70) arrClusEBinningCoarse[i]        = 4.+0.2*(i-50);
    else if(i<74) arrClusEBinningCoarse[i]        = 8.+0.5*(i-70);
    else if(i<90) arrClusEBinningCoarse[i]        = 10.+1.0*(i-74);
    else if(i<97) arrClusEBinningCoarse[i]        = 26.+2.0*(i-90);
    else if(i<109) arrClusEBinningCoarse[i]       = 40.+5.0*(i-97);
    else if(i<119) arrClusEBinningCoarse[i]       = 100.+10.0*(i-109);
    else arrClusEBinningCoarse[i]                 = 200;
  }
  for(Int_t i=0; i<nBinsClusterEOnlyHighPt+1;i++){
    if(i<25) arrClusEBinningOnlyHighPt[i]         = 10.+0.2*i;
    else if(i<55) arrClusEBinningOnlyHighPt[i]    = 15.+0.5*(i-25);
    else if(i<75) arrClusEBinningOnlyHighPt[i]    = 30.+1.0*(i-55);
    else if(i<100) arrClusEBinningOnlyHighPt[i]   = 50.+2.0*(i-75);
    else if(i<120) arrClusEBinningOnlyHighPt[i]   = 100.+5.0*(i-100);
    else if(i<130) arrClusEBinningOnlyHighPt[i]   = 200.+10.0*(i-120);
    else if(i<134) arrClusEBinningOnlyHighPt[i]   = 300.+25.0*(i-130);
    else arrClusEBinningOnlyHighPt[i]             = maxClusterE;
  }

  // IsPhotonSelected
  fHistCutIndex                   = new TH1F(Form("IsPhotonSelected %s",GetCutNumber().Data()),"IsPhotonSelected",5,-0.5,4.5);
  fHistCutIndex->GetXaxis()->SetBinLabel(kPhotonIn+1,"in");
  fHistCutIndex->GetXaxis()->SetBinLabel(kDetector+1,"detector");
  fHistCutIndex->GetXaxis()->SetBinLabel(kAcceptance+1,"acceptance");
  fHistCutIndex->GetXaxis()->SetBinLabel(kClusterQuality+1,"cluster QA");
  fHistCutIndex->GetXaxis()->SetBinLabel(kPhotonOut+1,"out");
  fHistograms->Add(fHistCutIndex);

  // Acceptance Cuts
  fHistAcceptanceCuts             = new TH1F(Form("AcceptanceCuts %s",GetCutNumber().Data()),"AcceptanceCuts",5,-0.5,4.5);
  fHistAcceptanceCuts->GetXaxis()->SetBinLabel(1,"in");
  fHistAcceptanceCuts->GetXaxis()->SetBinLabel(2,"eta");
  fHistAcceptanceCuts->GetXaxis()->SetBinLabel(3,"phi");
  fHistAcceptanceCuts->GetXaxis()->SetBinLabel(4,"dist bad channel/acc map");
  fHistAcceptanceCuts->GetXaxis()->SetBinLabel(5,"out");
  fHistograms->Add(fHistAcceptanceCuts);

  // Cluster Cuts
  if (GetIsConversionRecovery() == 0)
    fHistClusterIdentificationCuts  = new TH2F(Form("ClusterQualityCuts vs E %s",GetCutNumber().Data()),"ClusterQualityCuts",11,-0.5,10.5,nBinsClusterE, arrClusEBinning);
  else
    fHistClusterIdentificationCuts  = new TH2F(Form("ClusterQualityCuts vs E %s",GetCutNumber().Data()),"ClusterQualityCuts",12,-0.5,11.5,nBinsClusterE, arrClusEBinning);
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(1,"in");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(2,"timing");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(3,"Exotics");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(4,"minimum NCells");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(5,"NLM");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(6,"M02");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(7,"M20");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(8,"dispersion");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(9,"track matching");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(10,"minimum energy");
  fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(11,"out");
  if (GetIsConversionRecovery() > 0)
    fHistClusterIdentificationCuts->GetXaxis()->SetBinLabel(12,"rejected by conv Cl. flag");
  fHistograms->Add(fHistClusterIdentificationCuts);

  // Acceptance related histogramms
  if( fClusterType == 1 ){ //EMCAL
    const Int_t nEmcalEtaBins             = 96;
    const Int_t nEmcalPhiBins             = 124;
    Float_t EmcalEtaBins[nEmcalEtaBins+1] = {-0.66687,-0.653,-0.63913,-0.62526,-0.61139,-0.59752,-0.58365,-0.56978,-0.55591,-0.54204,-0.52817,-0.5143,-0.50043,-0.48656,-0.47269,-0.45882,-0.44495,-0.43108,-0.41721,-0.40334,-0.38947,-0.3756,-0.36173,-0.34786,-0.33399,-0.32012,-0.30625,-0.29238,-0.27851,-0.26464,-0.25077,-0.2369,-0.22303,-0.20916,-0.19529,-0.18142,-0.16755,-0.15368,-0.13981,-0.12594,-0.11207,-0.0982,-0.08433,-0.07046,-0.05659,-0.04272,-0.02885,-0.01498,-0.00111,0.01276,0.02663,0.0405,0.05437,0.06824,0.08211,0.09598,0.10985,0.12372,0.13759,0.15146,0.16533,0.1792,0.19307,0.20694,0.22081,0.23468,0.24855,0.26242,0.27629,0.29016,0.30403,0.3179,0.33177,0.34564,0.35951,0.37338,0.38725,0.40112,0.41499,0.42886,0.44273,0.4566,0.47047,0.48434,0.49821,0.51208,0.52595,0.53982,0.55369,0.56756,0.58143,0.5953,0.60917,0.62304,0.63691,0.65078,0.66465};
    Float_t EmcalPhiBins[nEmcalPhiBins+1] = {1.408,1.4215,1.435,1.4485,1.462,1.4755,1.489,1.5025,1.516,1.5295,1.543,1.5565,1.57,1.5835,1.597,1.6105,1.624,1.6375,1.651,1.6645,1.678,1.6915,1.705,1.7185,1.732,1.758,1.7715,1.785,1.7985,1.812,1.8255,1.839,1.8525,1.866,1.8795,1.893,1.9065,1.92,1.9335,1.947,1.9605,1.974,1.9875,2.001,2.0145,2.028,2.0415,2.055,2.0685,2.082,2.108,2.1215,2.135,2.1485,2.162,2.1755,2.189,2.2025,2.216,2.2295,2.243,2.2565,2.27,2.2835,2.297,2.3105,2.324,2.3375,2.351,2.3645,2.378,2.3915,2.405,2.4185,2.432,2.456,2.4695,2.483,2.4965,2.51,2.5235,2.537,2.5505,2.564,2.5775,2.591,2.6045,2.618,2.6315,2.645,2.6585,2.672,2.6855,2.699,2.7125,2.726,2.7395,2.753,2.7665,2.78,2.804,2.8175,2.831,2.8445,2.858,2.8715,2.885,2.8985,2.912,2.9255,2.939,2.9525,2.966,2.9795,2.993,3.0065,3.02,3.0335,3.047,3.0605,3.074,3.0875,3.101,3.1145,3.128};

    if(!fDoLightOutput){
      fHistClusterEtavsPhiBeforeAcc = new TH2F(Form("EtaPhi_beforeAcceptance %s",GetCutNumber().Data()),"EtaPhi_beforeAcceptance",nEmcalPhiBins,EmcalPhiBins,nEmcalEtaBins,EmcalEtaBins);
      fHistClusterEtavsPhiBeforeAcc->GetXaxis()->SetTitle("#varphi (rad)");
      fHistClusterEtavsPhiBeforeAcc->GetYaxis()->SetTitle("#eta");
      fHistograms->Add(fHistClusterEtavsPhiBeforeAcc);
      fHistClusterEtavsPhiAfterAcc  = new TH2F(Form("EtaPhi_afterAcceptance %s",GetCutNumber().Data()),"EtaPhi_afterAcceptance",nEmcalPhiBins,EmcalPhiBins,nEmcalEtaBins,EmcalEtaBins);
      fHistClusterEtavsPhiAfterAcc->GetXaxis()->SetTitle("#varphi (rad)");
      fHistClusterEtavsPhiAfterAcc->GetYaxis()->SetTitle("#eta");
      fHistograms->Add(fHistClusterEtavsPhiAfterAcc);
      if(fUseEtaPhiMapForBackCand){
        fHistClusterEtavsPhiBackground  = new TH2F(Form("EtaPhi_Background %s",GetCutNumber().Data()),"EtaPhi_Background",nEmcalPhiBins,EmcalPhiBins,nEmcalEtaBins,EmcalEtaBins);
        fHistClusterEtavsPhiBackground->GetXaxis()->SetTitle("#varphi (rad)");
        fHistClusterEtavsPhiBackground->GetYaxis()->SetTitle("#eta");
        fHistograms->Add(fHistClusterEtavsPhiBackground);
      }
    }
    fHistClusterEtavsPhiAfterQA     = new TH2F(Form("EtaPhi_afterClusterQA %s",GetCutNumber().Data()),"EtaPhi_afterClusterQA",nEmcalPhiBins,EmcalPhiBins,nEmcalEtaBins,EmcalEtaBins);
    fHistClusterEtavsPhiAfterQA->GetXaxis()->SetTitle("#varphi (rad)");
    fHistClusterEtavsPhiAfterQA->GetYaxis()->SetTitle("#eta");
    fHistograms->Add(fHistClusterEtavsPhiAfterQA);
  } else if( fClusterType == 2 ){ //PHOS
    const Int_t nPhosEtaBins        = 56;
    const Int_t nPhosPhiBins        = 256;
    const Float_t PhosEtaRange[2]   = {-0.16, 0.16};
    const Float_t PhosPhiRange[2]   = {4.355, 5.6};

    if(!fDoLightOutput){
      fHistClusterEtavsPhiBeforeAcc = new TH2F(Form("EtaPhi_beforeAcceptance %s",GetCutNumber().Data()),"EtaPhi_beforeAcceptance",nPhosPhiBins,PhosPhiRange[0],PhosPhiRange[1],nPhosEtaBins,PhosEtaRange[0],PhosEtaRange[1]);
      fHistClusterEtavsPhiBeforeAcc->GetXaxis()->SetTitle("#varphi (rad)");
      fHistClusterEtavsPhiBeforeAcc->GetYaxis()->SetTitle("#eta");
      fHistograms->Add(fHistClusterEtavsPhiBeforeAcc);
      fHistClusterEtavsPhiAfterAcc  = new TH2F(Form("EtaPhi_afterAcceptance %s",GetCutNumber().Data()),"EtaPhi_afterAcceptance",nPhosPhiBins,PhosPhiRange[0],PhosPhiRange[1],nPhosEtaBins,PhosEtaRange[0],PhosEtaRange[1]);
      fHistClusterEtavsPhiAfterAcc->GetXaxis()->SetTitle("#varphi (rad)");
      fHistClusterEtavsPhiAfterAcc->GetYaxis()->SetTitle("#eta");
      fHistograms->Add(fHistClusterEtavsPhiAfterAcc);
      if(fUseEtaPhiMapForBackCand){
        fHistClusterEtavsPhiBackground  = new TH2F(Form("EtaPhi_Background %s",GetCutNumber().Data()),"EtaPhi_Background",nPhosPhiBins,PhosPhiRange[0],PhosPhiRange[1],nPhosEtaBins,PhosEtaRange[0],PhosEtaRange[1]);
        fHistClusterEtavsPhiBackground->GetXaxis()->SetTitle("#varphi (rad)");
        fHistClusterEtavsPhiBackground->GetYaxis()->SetTitle("#eta");
        fHistograms->Add(fHistClusterEtavsPhiBackground);
      }
    }
    fHistClusterEtavsPhiAfterQA     = new TH2F(Form("EtaPhi_afterClusterQA %s",GetCutNumber().Data()),"EtaPhi_afterClusterQA",nPhosPhiBins,PhosPhiRange[0],PhosPhiRange[1],nPhosEtaBins,PhosEtaRange[0],PhosEtaRange[1]);
    fHistClusterEtavsPhiAfterQA->GetXaxis()->SetTitle("#varphi (rad)");
    fHistClusterEtavsPhiAfterQA->GetYaxis()->SetTitle("#eta");
    fHistograms->Add(fHistClusterEtavsPhiAfterQA);
  } else if( fClusterType == 3){ //DCAL
    const Int_t nDcalEtaBins             = 96;
    const Int_t nDcalPhiBins             = 124;
    if(!fDoLightOutput){
      fHistClusterEtavsPhiBeforeAcc = new TH2F(Form("EtaPhi_beforeAcceptance %s",GetCutNumber().Data()),"EtaPhi_beforeAcceptance",nDcalPhiBins,4.5,5.7,nDcalEtaBins,-0.66687,0.66465);
      fHistClusterEtavsPhiBeforeAcc->GetXaxis()->SetTitle("#varphi (rad)");
      fHistClusterEtavsPhiBeforeAcc->GetYaxis()->SetTitle("#eta");
      fHistograms->Add(fHistClusterEtavsPhiBeforeAcc);
      fHistClusterEtavsPhiAfterAcc  = new TH2F(Form("EtaPhi_afterAcceptance %s",GetCutNumber().Data()),"EtaPhi_afterAcceptance",nDcalPhiBins,4.5,5.7,nDcalEtaBins,-0.66687,0.66465);
      fHistClusterEtavsPhiAfterAcc->GetXaxis()->SetTitle("#varphi (rad)");
      fHistClusterEtavsPhiAfterAcc->GetYaxis()->SetTitle("#eta");
      fHistograms->Add(fHistClusterEtavsPhiAfterAcc);
      if(fUseEtaPhiMapForBackCand){
        fHistClusterEtavsPhiBackground  = new TH2F(Form("EtaPhi_Background %s",GetCutNumber().Data()),"EtaPhi_Background",nDcalPhiBins,4.5,5.7,nDcalEtaBins,-0.66687,0.66465);
        fHistClusterEtavsPhiBackground->GetXaxis()->SetTitle("#varphi (rad)");
        fHistClusterEtavsPhiBackground->GetYaxis()->SetTitle("#eta");
        fHistograms->Add(fHistClusterEtavsPhiBackground);
      }
    }
    fHistClusterEtavsPhiAfterQA     = new TH2F(Form("EtaPhi_afterClusterQA %s",GetCutNumber().Data()),"EtaPhi_afterClusterQA",nDcalPhiBins,4.5,5.7,nDcalEtaBins,-0.66687,0.66465);
    fHistClusterEtavsPhiAfterQA->GetXaxis()->SetTitle("#varphi (rad)");
    fHistClusterEtavsPhiAfterQA->GetYaxis()->SetTitle("#eta");
    fHistograms->Add(fHistClusterEtavsPhiAfterQA);
  } else if( fClusterType == 4){ //EMCAL+DCAL
    if(!fDoLightOutput){
      fHistClusterEtavsPhiBeforeAcc = new TH2F(Form("EtaPhi_beforeAcceptance %s",GetCutNumber().Data()),"EtaPhi_beforeAcceptance",462,0,2*TMath::Pi(),110,-0.7,0.7);
      fHistClusterEtavsPhiBeforeAcc->GetXaxis()->SetTitle("#varphi (rad)");
      fHistClusterEtavsPhiBeforeAcc->GetYaxis()->SetTitle("#eta");
      fHistograms->Add(fHistClusterEtavsPhiBeforeAcc);
      fHistClusterEtavsPhiAfterAcc  = new TH2F(Form("EtaPhi_afterAcceptance %s",GetCutNumber().Data()),"EtaPhi_afterAcceptance",462,0,2*TMath::Pi(),110,-0.7,0.7);
      fHistClusterEtavsPhiAfterAcc->GetXaxis()->SetTitle("#varphi (rad)");
      fHistClusterEtavsPhiAfterAcc->GetYaxis()->SetTitle("#eta");
      fHistograms->Add(fHistClusterEtavsPhiAfterAcc);
      if(fUseEtaPhiMapForBackCand){
        fHistClusterEtavsPhiBackground  = new TH2F(Form("EtaPhi_Background %s",GetCutNumber().Data()),"EtaPhi_Background",462,0,2*TMath::Pi(),110,-0.7,0.7);
        fHistClusterEtavsPhiBackground->GetXaxis()->SetTitle("#varphi (rad)");
        fHistClusterEtavsPhiBackground->GetYaxis()->SetTitle("#eta");
        fHistograms->Add(fHistClusterEtavsPhiBackground);
      }
    }
    fHistClusterEtavsPhiAfterQA     = new TH2F(Form("EtaPhi_afterClusterQA %s",GetCutNumber().Data()),"EtaPhi_afterClusterQA",462,0,2*TMath::Pi(),110,-0.7,0.7);
    fHistClusterEtavsPhiAfterQA->GetXaxis()->SetTitle("#varphi (rad)");
    fHistClusterEtavsPhiAfterQA->GetYaxis()->SetTitle("#eta");
    fHistograms->Add(fHistClusterEtavsPhiAfterQA);
  } else if( fClusterType == 0 ){ //all
    if(!fDoLightOutput){
      fHistClusterEtavsPhiBeforeAcc = new TH2F(Form("EtaPhi_beforeAcceptance %s",GetCutNumber().Data()),"EtaPhi_beforeAcceptance",462,0,2*TMath::Pi(),110,-0.7,0.7);
      fHistClusterEtavsPhiBeforeAcc->GetXaxis()->SetTitle("#varphi (rad)");
      fHistClusterEtavsPhiBeforeAcc->GetYaxis()->SetTitle("#eta");
      fHistograms->Add(fHistClusterEtavsPhiBeforeAcc);
      fHistClusterEtavsPhiAfterAcc  = new TH2F(Form("EtaPhi_afterAcceptance %s",GetCutNumber().Data()),"EtaPhi_afterAcceptance",462,0,2*TMath::Pi(),110,-0.7,0.7);
      fHistClusterEtavsPhiAfterAcc->GetXaxis()->SetTitle("#varphi (rad)");
      fHistClusterEtavsPhiAfterAcc->GetYaxis()->SetTitle("#eta");
      fHistograms->Add(fHistClusterEtavsPhiAfterAcc);
      if(fUseEtaPhiMapForBackCand){
        fHistClusterEtavsPhiBackground  = new TH2F(Form("EtaPhi_Background %s",GetCutNumber().Data()),"EtaPhi_Background",462,0,2*TMath::Pi(),110,-0.7,0.7);
        fHistClusterEtavsPhiBackground->GetXaxis()->SetTitle("#varphi (rad)");
        fHistClusterEtavsPhiBackground->GetYaxis()->SetTitle("#eta");
        fHistograms->Add(fHistClusterEtavsPhiBackground);
      }
    }
    fHistClusterEtavsPhiAfterQA     = new TH2F(Form("EtaPhi_afterClusterQA %s",GetCutNumber().Data()),"EtaPhi_afterClusterQA",462,0,2*TMath::Pi(),110,-0.7,0.7);
    fHistClusterEtavsPhiAfterQA->GetXaxis()->SetTitle("#varphi (rad)");
    fHistClusterEtavsPhiAfterQA->GetYaxis()->SetTitle("#eta");
    fHistograms->Add(fHistClusterEtavsPhiAfterQA);
  } else{AliError(Form("Cluster Type is not EMCAL nor PHOS nor all: %i",fClusterType));}

  if(fIsMC > 1){
    if(!fDoLightOutput){
      fHistClusterEtavsPhiBeforeAcc->Sumw2();
      fHistClusterEtavsPhiAfterAcc->Sumw2();
      if(fUseEtaPhiMapForBackCand){
        fHistClusterEtavsPhiBackground->Sumw2();
      }
    }
    fHistClusterEtavsPhiAfterQA->Sumw2();
  }

  // Cluster quality related histograms
  Double_t timeMin                  = -2e-6;
  Double_t timeMax                  = 8e-6;
  if( fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
    timeMin                         = -2e-7;
    timeMax                         = 12e-7;
  }

  if(!fDoLightOutput){
    fHistClusterTimevsEBeforeQA     = new TH2F(Form("ClusterTimeVsE_beforeClusterQA %s",GetCutNumber().Data()),"ClusterTimeVsE_beforeClusterQA",800, timeMin, timeMax,
                                               nBinsClusterE, arrClusEBinning);
    fHistClusterTimevsEBeforeQA->GetXaxis()->SetTitle("t_{cl} (s)");
    fHistClusterTimevsEBeforeQA->GetYaxis()->SetTitle("E_{cl} (GeV)");
    fHistograms->Add(fHistClusterTimevsEBeforeQA);
  }
  fHistClusterTimevsEAfterQA        = new TH2F(Form("ClusterTimeVsE_afterClusterQA %s",GetCutNumber().Data()),"ClusterTimeVsE_afterClusterQA",800, timeMin, timeMax,
                                               nBinsClusterE, arrClusEBinning);
  fHistClusterTimevsEAfterQA->GetXaxis()->SetTitle("t_{cl} (s)");
  fHistClusterTimevsEAfterQA->GetYaxis()->SetTitle("E_{cl} (GeV)");
  fHistograms->Add(fHistClusterTimevsEAfterQA);
  if(fUseNonLinearity && !fDoLightOutput){
    fHistEnergyOfClusterBeforeNL    = new TH1F(Form("EnergyOfCluster_beforeNonLinearity %s",GetCutNumber().Data()), "EnergyOfCluster_beforeNonLinearity",
                                               nBinsClusterE, arrClusEBinning);
    fHistEnergyOfClusterBeforeNL->GetXaxis()->SetTitle("E_{cl} (GeV)");
    fHistograms->Add(fHistEnergyOfClusterBeforeNL);
    fHistEnergyOfClusterAfterNL     = new TH1F(Form("EnergyOfCluster_afterNonLinearity %s",GetCutNumber().Data()), "EnergyOfCluster_afterNonLinearity",
                                               nBinsClusterE, arrClusEBinning);
    fHistEnergyOfClusterAfterNL->GetXaxis()->SetTitle("E_{cl} (GeV)");
    fHistograms->Add(fHistEnergyOfClusterAfterNL);
  }
  fHistEnergyOfClusterBeforeQA      = new TH1F(Form("EnergyOfCluster_beforeClusterQA %s",GetCutNumber().Data()),"EnergyOfCluster_beforeClusterQA", nBinsClusterE, arrClusEBinning);
  fHistEnergyOfClusterBeforeQA->GetXaxis()->SetTitle("E_{cl} (GeV)");
  fHistograms->Add(fHistEnergyOfClusterBeforeQA);
  fHistEnergyOfClusterAfterQA       = new TH1F(Form("EnergyOfCluster_afterClusterQA %s",GetCutNumber().Data()),"EnergyOfCluster_afterClusterQA",nBinsClusterE, arrClusEBinning);
  fHistEnergyOfClusterAfterQA->GetXaxis()->SetTitle("E_{cl} (GeV)");
  fHistograms->Add(fHistEnergyOfClusterAfterQA);
  if(!fDoLightOutput){
    fHistNCellsBeforeQA             = new TH1F(Form("NCellPerCluster_beforeClusterQA %s",GetCutNumber().Data()),"NCellPerCluster_beforeClusterQA",50,0,50);
    fHistNCellsBeforeQA->GetXaxis()->SetTitle("N_{cells} in cluster");
    fHistograms->Add(fHistNCellsBeforeQA);
    fHistNCellsAfterQA              = new TH1F(Form("NCellPerCluster_afterClusterQA %s",GetCutNumber().Data()),"NCellPerCluster_afterClusterQA",50,0,50);
    fHistNCellsAfterQA->GetXaxis()->SetTitle("N_{cells} in cluster");
    fHistograms->Add(fHistNCellsAfterQA);
    fHistM02BeforeQA                = new TH1F(Form("M02_beforeClusterQA %s",GetCutNumber().Data()),"M02_beforeClusterQA",400,0,5);
    fHistM02BeforeQA->GetXaxis()->SetTitle("#sigma_{long}^2");
    fHistograms->Add(fHistM02BeforeQA);
    fHistM02AfterQA                 = new TH1F(Form("M02_afterClusterQA %s",GetCutNumber().Data()),"M02_afterClusterQA",400,0,5);
    fHistM02AfterQA->GetXaxis()->SetTitle("#sigma_{long}^2");
    fHistograms->Add(fHistM02AfterQA);
    fHistM20BeforeQA                = new TH1F(Form("M20_beforeClusterQA %s",GetCutNumber().Data()),"M20_beforeClusterQA",400,0,2.5);
    fHistM20BeforeQA->GetXaxis()->SetTitle("#sigma_{short}^2");
    fHistograms->Add(fHistM20BeforeQA);
    fHistM20AfterQA                 = new TH1F(Form("M20_afterClusterQA %s",GetCutNumber().Data()),"M20_afterClusterQA",400,0,2.5);
    fHistM20AfterQA->GetXaxis()->SetTitle("#sigma_{short}^2");
    fHistograms->Add(fHistM20AfterQA);
    fHistDispersionBeforeQA         = new TH1F(Form("Dispersion_beforeClusterQA %s",GetCutNumber().Data()),"Dispersion_beforeClusterQA",200,0,20);
    fHistDispersionBeforeQA->GetXaxis()->SetTitle("dispersion");
    fHistograms->Add(fHistDispersionBeforeQA);
    fHistDispersionAfterQA          = new TH1F(Form("Dispersion_afterClusterQA %s",GetCutNumber().Data()),"Dispersion_afterClusterQA",200,0,20);
    fHistDispersionAfterQA->GetXaxis()->SetTitle("dispersion");
    fHistograms->Add(fHistDispersionAfterQA);
    fHistNLMBeforeQA                = new TH1F(Form("NLM_beforeClusterQA %s",GetCutNumber().Data()),"NLM_beforeClusterQA",10,0,10);
    fHistNLMBeforeQA->GetXaxis()->SetTitle("N_{LM} in cluster");
    fHistograms->Add(fHistNLMBeforeQA);
    fHistNLMAfterQA                 = new TH1F(Form("NLM_afterClusterQA %s",GetCutNumber().Data()),"NLM_afterClusterQA",10,0,10);
    fHistNLMAfterQA->GetXaxis()->SetTitle("N_{LM} in cluster");
    fHistograms->Add(fHistNLMAfterQA);
//     fHistNLMAvsNLMBBeforeQA         = new TH2F(Form("NLMAvsNLMB_beforeClusterQA %s",GetCutNumber().Data()),"NLMAvsNLMB_beforeClusterQA",10,0,10,10,0,10);
//     fHistograms->Add(fHistNLMAvsNLMBBeforeQA);
    fHistNLMVsNCellsAfterQA         = new TH2F(Form("NLM_NCells_afterClusterQA %s",GetCutNumber().Data()),"NLM_NCells_afterClusterQA",10,0,10,50,0,50);
    fHistNLMVsNCellsAfterQA->GetXaxis()->SetTitle("N_{LM} in cluster");
    fHistNLMVsNCellsAfterQA->GetYaxis()->SetTitle("N_{cells} in cluster");
    fHistograms->Add(fHistNLMVsNCellsAfterQA);
    fHistNLMVsEAfterQA              = new TH2F(Form("NLM_E_afterClusterQA %s",GetCutNumber().Data()),"NLM_E_afterClusterQA",10,0,10,nBinsClusterE, arrClusEBinning);
    fHistNLMVsEAfterQA->GetXaxis()->SetTitle("N_{LM} in cluster");
    fHistNLMVsEAfterQA->GetYaxis()->SetTitle("E_{cl} (GeV)");
    fHistograms->Add(fHistNLMVsEAfterQA);

    if (GetIsConversionRecovery() > 0){
      fHistInvMassDiCluster      = new TH2F(Form("InvMass_ClusterPair %s",GetCutNumber().Data()),"InvMass_ClusterPair",100,0,0.1,100,0,20);
      fHistInvMassDiCluster->GetXaxis()->SetTitle("M_{cl,cl} (GeV/c^{2})");
      fHistInvMassDiCluster->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      fHistograms->Add(fHistInvMassDiCluster);
      fHistInvMassConvFlagging      = new TH2F(Form("InvMass_Rejected_ClusterPair %s",GetCutNumber().Data()),"InvMass_Rejected_ClusterPair",100,0,0.1,100,0,20);
      fHistInvMassConvFlagging->GetXaxis()->SetTitle("M_{cl,cl} (GeV/c^{2})");
      fHistInvMassConvFlagging->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      fHistograms->Add(fHistInvMassConvFlagging);
    }

    if(fExtendedMatchAndQA > 0 || fIsPureCalo > 0){
      fHistClusterEM02AfterQA       = new TH2F(Form("EVsM02_afterClusterQA %s",GetCutNumber().Data()),"EVsM02_afterClusterQA",nBinsClusterE, arrClusEBinning,400,0,5);
      fHistClusterEM02AfterQA->GetYaxis()->SetTitle("#sigma_{long}^2");
      fHistClusterEM02AfterQA->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistograms->Add(fHistClusterEM02AfterQA);
      fHistClusterEM02BeforeQA      = new TH2F(Form("EVsM02_beforeClusterQA %s",GetCutNumber().Data()),"EVsM02_beforeClusterQA",nBinsClusterE, arrClusEBinning,400,0,5);
      fHistClusterEM02BeforeQA->GetYaxis()->SetTitle("#sigma_{long}^2");
      fHistClusterEM02BeforeQA->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistograms->Add(fHistClusterEM02BeforeQA);
      // QA histos for shower shape correlations
      fHistClusterM02M20BeforeQA                      = new TH2F(Form("M02VsM20_beforeClusterQA %s",GetCutNumber().Data()),"M02VsM20_beforeClusterQA",400,0,5,200,0,2.5);
      fHistClusterM02M20BeforeQA->GetXaxis()->SetTitle("#sigma_{short}^2");
      fHistClusterM02M20BeforeQA->GetYaxis()->SetTitle("#sigma_{long}^2");
      fHistograms->Add(fHistClusterM02M20BeforeQA);
      fHistClusterM02M20AfterQA                       = new TH2F(Form("M02VsM20_afterClusterQA %s",GetCutNumber().Data()),"M02VsM20_afterClusterQA",400,0,5,200,0,2.5);
      fHistClusterM02M20AfterQA->GetXaxis()->SetTitle("#sigma_{short}^2");
      fHistClusterM02M20AfterQA->GetYaxis()->SetTitle("#sigma_{long}^2");
      fHistograms->Add(fHistClusterM02M20AfterQA);
    }
  }
  //----------------
  if(fIsMC > 1){
    if(!fDoLightOutput) fHistClusterTimevsEBeforeQA->Sumw2();
    fHistClusterTimevsEAfterQA->Sumw2();
    if(fUseNonLinearity && !fDoLightOutput){
      fHistEnergyOfClusterBeforeNL->Sumw2();
      fHistEnergyOfClusterAfterNL->Sumw2();
    }
    fHistEnergyOfClusterBeforeQA->Sumw2();
    fHistEnergyOfClusterAfterQA->Sumw2();
    if(!fDoLightOutput){
      fHistNCellsBeforeQA->Sumw2();
      fHistNCellsAfterQA->Sumw2();
      fHistM02BeforeQA->Sumw2();
      fHistM02AfterQA->Sumw2();
      fHistM20BeforeQA->Sumw2();
      fHistM20AfterQA->Sumw2();
      fHistDispersionBeforeQA->Sumw2();
      fHistDispersionAfterQA->Sumw2();
      fHistNLMBeforeQA->Sumw2();
      fHistNLMAfterQA->Sumw2();
//       fHistNLMAvsNLMBBeforeQA->Sumw2();
      fHistNLMVsNCellsAfterQA->Sumw2();
      fHistNLMVsEAfterQA->Sumw2();
      if(fExtendedMatchAndQA > 0 || fIsPureCalo > 0){
        fHistClusterEM02AfterQA->Sumw2();
        fHistClusterEM02BeforeQA->Sumw2();
        fHistClusterM02M20BeforeQA->Sumw2();
        fHistClusterM02M20AfterQA->Sumw2();
      }
    }
  }
//----------------
  if(!fDoLightOutput){
    TString namePeriod = ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->GetPeriodName();
    if( fClusterType == 1 ){
      Int_t nMaxCellsEMCAL      = fNMaxEMCalModules*48*24;
      fBadChannels              = new TProfile("EMCal - Bad Channels","EMCal - Bad Channels",nMaxCellsEMCAL,-0.5,nMaxCellsEMCAL-0.5);
      fBadChannels->GetXaxis()->SetTitle("channel ID");
      fHistExtQA->Add(fBadChannels);
    } else if( fClusterType == 2 ){
      if (namePeriod.Contains("LHC10")) fNMaxPHOSModules=3;
      Int_t nMaxCellsPHOS       = fNMaxPHOSModules*56*64;
      fBadChannels              = new TProfile("PHOS - Bad Channels","PHOS - Bad Channels",nMaxCellsPHOS,-0.5,nMaxCellsPHOS-0.5);
      fBadChannels->GetXaxis()->SetTitle("channel ID");
      fHistExtQA->Add(fBadChannels);
    } else if( fClusterType == 3 ){
      Int_t nStartCellDCAL      = 12288;
      Int_t nMaxCellsDCAL       = fNMaxDCalModules*32*24;
      fBadChannels              = new TProfile("DCAL - Bad Channels","DCAL - Bad Channels",nMaxCellsDCAL,nStartCellDCAL-0.5,nStartCellDCAL+nMaxCellsDCAL-0.5);
      fBadChannels->GetXaxis()->SetTitle("channel ID");
      fHistExtQA->Add(fBadChannels);
    } else if( fClusterType == 4 ){
      Int_t nMaxCellsEMCAL      = fNMaxEMCalModules*48*24;
      Int_t nMaxCellsDCAL       = fNMaxDCalModules*32*24;
      fBadChannels              = new TProfile("EMCAL+DCAL - Bad Channels","EMCAL+DCAL - Bad Channels",nMaxCellsEMCAL+nMaxCellsDCAL,-0.5,nMaxCellsEMCAL+nMaxCellsDCAL-0.5);
      fBadChannels->GetXaxis()->SetTitle("channel ID");
      fHistExtQA->Add(fBadChannels);
    }
  }

  if(fExtendedMatchAndQA > 1){
    if( fClusterType == 1 ){ //EMCAL
      // detailed cluster QA histos for EMCAL
      Int_t nMaxCellsEMCAL            = fNMaxEMCalModules*48*24;
      fHistClusterEnergyvsMod         = new TH2F(Form("ClusterEnergyVsModule_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyVsModule_afterClusterQA",
                                                 nBinsClusterE, arrClusEBinning, fNMaxEMCalModules, 0, fNMaxEMCalModules);
      fHistClusterEnergyvsMod->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistClusterEnergyvsMod->GetYaxis()->SetTitle("module ID");
      fHistExtQA->Add(fHistClusterEnergyvsMod);
      fHistNCellsBigger100MeVvsMod    = new TH2F(Form("NCellsAbove100VsModule %s",GetCutNumber().Data()),"NCellsAbove100VsModule",200,0,200,fNMaxEMCalModules,0,fNMaxEMCalModules);
      fHistNCellsBigger100MeVvsMod->GetXaxis()->SetTitle("N_{cells} with E_{cell} > 0.1 GeV");
      fHistNCellsBigger100MeVvsMod->GetYaxis()->SetTitle("module ID");
      fHistExtQA->Add(fHistNCellsBigger100MeVvsMod);
      fHistNCellsBigger1500MeVvsMod   = new TH2F(Form("NCellsAbove1500VsModule %s",GetCutNumber().Data()),"NCellsAbove1500VsModule",100,0,100,fNMaxEMCalModules,0,fNMaxEMCalModules);
      fHistNCellsBigger1500MeVvsMod->GetXaxis()->SetTitle("N_{cells} with E_{cell} > 1.5 GeV");
      fHistNCellsBigger1500MeVvsMod->GetYaxis()->SetTitle("module ID");
      fHistExtQA->Add(fHistNCellsBigger1500MeVvsMod);
      fHistEnergyOfModvsMod           = new TH2F(Form("ModuleEnergyVsModule %s",GetCutNumber().Data()),"ModuleEnergyVsModule", nBinsClusterEMod, arrClusEBinning,
                                                 fNMaxEMCalModules,0,fNMaxEMCalModules);
      fHistEnergyOfModvsMod->GetXaxis()->SetTitle("E_{mod} (GeV)");
      fHistEnergyOfModvsMod->GetYaxis()->SetTitle("module ID");
      fHistExtQA->Add(fHistEnergyOfModvsMod);
      fHistClusterEnergyvsNCellsBeforeQA      = new TH2F(Form("ClusterEnergyVsNCells_beforeQA %s",GetCutNumber().Data()),"ClusterEnergyVsNCells_beforeQA",nBinsClusterE, arrClusEBinning,
                                                          50, 0, 50);
      fHistClusterEnergyvsNCellsBeforeQA->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistClusterEnergyvsNCellsBeforeQA->GetYaxis()->SetTitle("N_{cells}");
      fHistExtQA->Add(fHistClusterEnergyvsNCellsBeforeQA);
      fHistClusterEnergyvsNCellsAfterQA      = new TH2F(Form("ClusterEnergyVsNCells_afterQA %s",GetCutNumber().Data()),"ClusterEnergyVsNCells_afterQA",nBinsClusterE, arrClusEBinning,
                                                 50, 0, 50);
      fHistClusterEnergyvsNCellsAfterQA->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistClusterEnergyvsNCellsAfterQA->GetYaxis()->SetTitle("N_{cells}");
      fHistExtQA->Add(fHistClusterEnergyvsNCellsAfterQA);
      fHistClusterDistanceInTimeCut   = new TH2F(Form("ClusterDistanceTo_withinTimingCut %s",GetCutNumber().Data()),"ClusterDistanceTo_withinTimingCut",20,-10,10,20,-10,10);
      fHistClusterDistanceInTimeCut->GetXaxis()->SetTitle("R_{cl,row} within time cut (cell)");
      fHistClusterDistanceInTimeCut->GetYaxis()->SetTitle("R_{cl,col} within time cut (cell)");

      fHistExtQA->Add(fHistClusterDistanceInTimeCut);
      fHistClusterDistanceOutTimeCut  = new TH2F(Form("ClusterDistanceTo_outsideTimingCut %s",GetCutNumber().Data()),"ClusterDistanceTo_outsideTimingCut",20,-10,10,20,-10,10);
      fHistClusterDistanceOutTimeCut->GetXaxis()->SetTitle("R_{cl,row} outside time cut (cell)");
      fHistClusterDistanceOutTimeCut->GetYaxis()->SetTitle("R_{cl,col} outside time cut (cell)");

      fHistExtQA->Add(fHistClusterDistanceOutTimeCut);
      fHistClusterDistance1DInTimeCut   = new TH1F(Form("Cluster1D_DistanceTo_withinTimingCut %s",GetCutNumber().Data()),"Cluster1D_DistanceTo_withinTimingCut",200,0.,0.5);
      fHistClusterDistance1DInTimeCut->GetXaxis()->SetTitle("R_{cl,1D} within time cut (cell)");

      fHistExtQA->Add(fHistClusterDistance1DInTimeCut);

      // detailed cell QA histos for EMCAL
      if(fExtendedMatchAndQA > 3){
        fHistCellEnergyvsCellID                 = new TH2F(Form("CellEnergyVsCellID %s",GetCutNumber().Data()),"CellEnergyVsCellID",nBinsClusterECellCoarse,  arrClusEBinningCoarse,
                                                           nMaxCellsEMCAL, 0, nMaxCellsEMCAL);
        fHistCellEnergyvsCellID->GetXaxis()->SetTitle("E_{cell} (GeV)");
        fHistCellEnergyvsCellID->GetYaxis()->SetTitle("Cell ID");
        fHistExtQA->Add(fHistCellEnergyvsCellID);
        fHistCellTimevsCellID                   = new TH2F(Form("CellTimeVsCellID %s",GetCutNumber().Data()),"CellTimeVsCellID",600,-timeMax,timeMax,nMaxCellsEMCAL,0,nMaxCellsEMCAL);
        fHistCellTimevsCellID->GetXaxis()->SetTitle("t_{cell} (GeV)");
        fHistCellTimevsCellID->GetYaxis()->SetTitle("Cell ID");
        fHistExtQA->Add(fHistCellTimevsCellID);
        fHistClusterIncludedCellsBeforeQA       = new TH1F(Form("ClusterIncludedCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_beforeClusterQA",
                                                           nMaxCellsEMCAL,0,nMaxCellsEMCAL);
        fHistClusterIncludedCellsBeforeQA->GetXaxis()->SetTitle("Cell ID");
        fHistClusterIncludedCellsBeforeQA->GetYaxis()->SetTitle("# included in cluster");
        fHistExtQA->Add(fHistClusterIncludedCellsBeforeQA);
        fHistClusterIncludedCellsAfterQA        = new TH1F(Form("ClusterIncludedCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_afterClusterQA",
                                                           nMaxCellsEMCAL,0,nMaxCellsEMCAL);
        fHistClusterIncludedCellsAfterQA->GetXaxis()->SetTitle("Cell ID");
        fHistClusterIncludedCellsAfterQA->GetYaxis()->SetTitle("# included in cluster");
        fHistExtQA->Add(fHistClusterIncludedCellsAfterQA);
        fHistClusterEnergyFracCellsBeforeQA     = new TH1F(Form("ClusterEnergyFracCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_beforeClusterQA",
                                                           nMaxCellsEMCAL,0,nMaxCellsEMCAL);
        fHistClusterEnergyFracCellsBeforeQA->GetXaxis()->SetTitle("Cell ID");
        fHistClusterEnergyFracCellsBeforeQA->GetYaxis()->SetTitle("fraction of energy cell is contrib. to cl.");
        fHistClusterEnergyFracCellsBeforeQA->Sumw2();
        fHistExtQA->Add(fHistClusterEnergyFracCellsBeforeQA);
        fHistClusterEnergyFracCellsAfterQA      = new TH1F(Form("ClusterEnergyFracCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_afterClusterQA",
                                                           nMaxCellsEMCAL,0,nMaxCellsEMCAL);
        fHistClusterEnergyFracCellsAfterQA->GetXaxis()->SetTitle("Cell ID");
        fHistClusterEnergyFracCellsAfterQA->GetYaxis()->SetTitle("fraction of energy cell is contrib. to cl.");
        fHistClusterEnergyFracCellsAfterQA->Sumw2();
        fHistExtQA->Add(fHistClusterEnergyFracCellsAfterQA);
        fHistClusterIncludedCellsTimingAfterQA  = new TH2F(Form("ClusterIncludedCellsTiming_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCellsTiming_afterClusterQA",
                                                           nBinsClusterE, arrClusEBinning, 200, -500, 500);
        fHistClusterIncludedCellsTimingAfterQA->GetXaxis()->SetTitle("E_{cl} (GeV)");
        fHistClusterIncludedCellsTimingAfterQA->GetYaxis()->SetTitle("t_{cell} (ns)");
        fHistExtQA->Add(fHistClusterIncludedCellsTimingAfterQA);
        fHistClusterIncludedCellsTimingEnergyAfterQA  = new TH2F(Form("ClusterIncludedCellsTimingEnergy_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCellsTimingEnergy_afterClusterQA",
                                                                 nBinsClusterECell,  arrClusEBinningCoarse, 200, -500, 500);
        fHistClusterIncludedCellsTimingEnergyAfterQA->GetXaxis()->SetTitle("E_{cell} (GeV)");
        fHistClusterIncludedCellsTimingEnergyAfterQA->GetYaxis()->SetTitle("t_{cell} (ns)");
        fHistExtQA->Add(fHistClusterIncludedCellsTimingEnergyAfterQA);
      }

    } else if( fClusterType == 2 ){ //PHOS
      // detailed cluster QA histos for PHOS
      Int_t nMaxCellsPHOS             = fNMaxPHOSModules*56*64;
      fHistClusterEnergyvsMod         = new TH2F(Form("ClusterEnergyVsModule_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyVsModule_afterClusterQA",nBinsClusterE, arrClusEBinning,
                                                 fNMaxPHOSModules, 0, fNMaxPHOSModules);
      fHistClusterEnergyvsMod->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistClusterEnergyvsMod->GetYaxis()->SetTitle("module ID");
      fHistExtQA->Add(fHistClusterEnergyvsMod);
      fHistNCellsBigger100MeVvsMod    = new TH2F(Form("NCellsAbove100VsModule %s",GetCutNumber().Data()),"NCellsAbove100VsModule",200,0,200,fNMaxPHOSModules,0,fNMaxPHOSModules);
      fHistNCellsBigger100MeVvsMod->GetXaxis()->SetTitle("N_{cells} with E_{cell} > 0.1 GeV");
      fHistNCellsBigger100MeVvsMod->GetYaxis()->SetTitle("module ID");
      fHistExtQA->Add(fHistNCellsBigger100MeVvsMod);
      fHistNCellsBigger1500MeVvsMod   = new TH2F(Form("NCellsAbove1500VsModule %s",GetCutNumber().Data()),"NCellsAbove1500VsModule",100,0,100,fNMaxPHOSModules,0,fNMaxPHOSModules);
      fHistExtQA->Add(fHistNCellsBigger1500MeVvsMod);
      fHistNCellsBigger1500MeVvsMod->GetXaxis()->SetTitle("N_{cells} with E_{cell} > 1.5 GeV");
      fHistNCellsBigger1500MeVvsMod->GetYaxis()->SetTitle("module ID");
      fHistEnergyOfModvsMod           = new TH2F(Form("ModuleEnergyVsModule %s",GetCutNumber().Data()),"ModuleEnergyVsModule",nBinsClusterEMod, arrClusEBinning,
                                                 fNMaxPHOSModules, 0, fNMaxPHOSModules);
      fHistEnergyOfModvsMod->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistEnergyOfModvsMod->GetYaxis()->SetTitle("N_{cells}");
      fHistExtQA->Add(fHistEnergyOfModvsMod);
      fHistClusterEnergyvsNCellsBeforeQA      = new TH2F(Form("ClusterEnergyVsNCells_beforeQA %s",GetCutNumber().Data()),"ClusterEnergyVsNCells_beforeQA",nBinsClusterE, arrClusEBinning, 50, 0, 50);
      fHistClusterEnergyvsNCellsBeforeQA->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistClusterEnergyvsNCellsBeforeQA->GetYaxis()->SetTitle("N_{cells}");
      fHistExtQA->Add(fHistClusterEnergyvsNCellsBeforeQA);
      fHistClusterEnergyvsNCellsAfterQA      = new TH2F(Form("ClusterEnergyVsNCells_afterQA %s",GetCutNumber().Data()),"ClusterEnergyVsNCells_afterQA",nBinsClusterE, arrClusEBinning, 50, 0, 50);
      fHistClusterEnergyvsNCellsAfterQA->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistClusterEnergyvsNCellsAfterQA->GetYaxis()->SetTitle("N_{cells}");
      fHistExtQA->Add(fHistClusterEnergyvsNCellsAfterQA);
      fHistClusterDistanceInTimeCut   = new TH2F(Form("ClusterDistanceTo_withinTimingCut %s",GetCutNumber().Data()),"ClusterDistanceTo_withinTimingCut",20,-10,10,20,-10,10);
      fHistClusterDistanceInTimeCut->GetXaxis()->SetTitle("R_{cl,row} within time cut (cell)");
      fHistClusterDistanceInTimeCut->GetYaxis()->SetTitle("R_{cl,col} within time cut (cell)");
      fHistExtQA->Add(fHistClusterDistanceInTimeCut);
      fHistClusterDistanceOutTimeCut  = new TH2F(Form("ClusterDistanceTo_outsideTimingCut %s",GetCutNumber().Data()),"ClusterDistanceTo_outsideTimingCut",20,-10,10,20,-10,10);
      fHistClusterDistanceOutTimeCut->GetXaxis()->SetTitle("R_{cl,row} outside time cut (cell)");
      fHistClusterDistanceOutTimeCut->GetYaxis()->SetTitle("R_{cl,col} outside time cut (cell)");
      fHistExtQA->Add(fHistClusterDistanceOutTimeCut);
      fHistClusterDistance1DInTimeCut   = new TH1F(Form("Cluster1D_DistanceTo_withinTimingCut %s",GetCutNumber().Data()),"Cluster1D_DistanceTo_withinTimingCut",200,0.,0.5);
      fHistClusterDistance1DInTimeCut->GetXaxis()->SetTitle("R_{cl,1D} within time cut (cell)");
      fHistExtQA->Add(fHistClusterDistance1DInTimeCut);

      // detailed cell QA histos for PHOS
      if(fExtendedMatchAndQA > 3){
        fHistCellEnergyvsCellID         = new TH2F(Form("CellEnergyVsCellID %s",GetCutNumber().Data()),"CellEnergyVsCellID",nBinsClusterECellCoarse,  arrClusEBinningCoarse,
                                                   nMaxCellsPHOS,0,nMaxCellsPHOS);
        fHistCellEnergyvsCellID->GetXaxis()->SetTitle("E_{cell} (GeV)");
        fHistCellEnergyvsCellID->GetYaxis()->SetTitle("Cell ID");
        fHistExtQA->Add(fHistCellEnergyvsCellID);
        fHistCellTimevsCellID           = new TH2F(Form("CellTimeVsCellID %s",GetCutNumber().Data()),"CellTimeVsCellID",600,-timeMax,timeMax,nMaxCellsPHOS,0,nMaxCellsPHOS);
        fHistCellTimevsCellID->GetXaxis()->SetTitle("t_{cell} (GeV)");
        fHistCellTimevsCellID->GetYaxis()->SetTitle("Cell ID");
        fHistExtQA->Add(fHistCellTimevsCellID);
        fHistClusterIncludedCellsBeforeQA   = new TH1F(Form("ClusterIncludedCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_beforeClusterQA",
                                                       nMaxCellsPHOS, 0, nMaxCellsPHOS);
        fHistClusterIncludedCellsBeforeQA->GetXaxis()->SetTitle("Cell ID");
        fHistClusterIncludedCellsBeforeQA->GetYaxis()->SetTitle("# included in cluster");
        fHistExtQA->Add(fHistClusterIncludedCellsBeforeQA);
        fHistClusterIncludedCellsAfterQA    = new TH1F(Form("ClusterIncludedCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_afterClusterQA",
                                                       nMaxCellsPHOS, 0, nMaxCellsPHOS);
        fHistClusterIncludedCellsAfterQA->GetXaxis()->SetTitle("Cell ID");
        fHistClusterIncludedCellsAfterQA->GetYaxis()->SetTitle("# included in cluster");
        fHistExtQA->Add(fHistClusterIncludedCellsAfterQA);
        fHistClusterEnergyFracCellsBeforeQA = new TH1F(Form("ClusterEnergyFracCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_beforeClusterQA",
                                                       nMaxCellsPHOS,0,nMaxCellsPHOS);
        fHistClusterEnergyFracCellsBeforeQA->GetXaxis()->SetTitle("Cell ID");
        fHistClusterEnergyFracCellsBeforeQA->GetYaxis()->SetTitle("fraction of energy cell is contrib. to cl.");
        fHistClusterEnergyFracCellsBeforeQA->Sumw2();
        fHistExtQA->Add(fHistClusterEnergyFracCellsBeforeQA);
        fHistClusterEnergyFracCellsAfterQA  = new TH1F(Form("ClusterEnergyFracCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_afterClusterQA",
                                                       nMaxCellsPHOS, 0, nMaxCellsPHOS);
        fHistClusterEnergyFracCellsAfterQA->GetXaxis()->SetTitle("Cell ID");
        fHistClusterEnergyFracCellsAfterQA->GetYaxis()->SetTitle("fraction of energy cell is contrib. to cl.");
        fHistClusterEnergyFracCellsAfterQA->Sumw2();
        fHistExtQA->Add(fHistClusterEnergyFracCellsAfterQA);
        fHistClusterIncludedCellsTimingAfterQA        = new TH2F(Form("ClusterIncludedCellsTiming_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCellsTiming_afterClusterQA",
                                                                 nBinsClusterE, arrClusEBinning, 200, -500, 500);
        fHistClusterIncludedCellsTimingAfterQA->GetXaxis()->SetTitle("E_{cl} (GeV)");
        fHistClusterIncludedCellsTimingAfterQA->GetYaxis()->SetTitle("t_{cell} (ns)");
        fHistExtQA->Add(fHistClusterIncludedCellsTimingAfterQA);
        fHistClusterIncludedCellsTimingEnergyAfterQA  = new TH2F(Form("ClusterIncludedCellsTimingEnergy_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCellsTimingEnergy_afterClusterQA",
                                                                 nBinsClusterECell,  arrClusEBinningCoarse, 200, -500, 500);
        fHistClusterIncludedCellsTimingEnergyAfterQA->GetXaxis()->SetTitle("E_{cell} (GeV)");
        fHistClusterIncludedCellsTimingEnergyAfterQA->GetYaxis()->SetTitle("t_{cell} (ns)");
        fHistExtQA->Add(fHistClusterIncludedCellsTimingEnergyAfterQA);
      }
    } else if( fClusterType == 3 ){ //DCAL
      Int_t nModulesStart = 12;
      Int_t nCellsStart = 12288;
      Int_t nMaxCellsDCAL            = fNMaxDCalModules*32*24;

      fHistClusterEnergyvsMod         = new TH2F(Form("ClusterEnergyVsModule_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyVsModule_afterClusterQA",
                                                 nBinsClusterE, arrClusEBinning, fNMaxDCalModules, nModulesStart, fNMaxDCalModules+nModulesStart);
      fHistClusterEnergyvsMod->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistClusterEnergyvsMod->GetYaxis()->SetTitle("module ID");
      fHistExtQA->Add(fHistClusterEnergyvsMod);
      fHistNCellsBigger100MeVvsMod    = new TH2F(Form("NCellsAbove100VsModule %s",GetCutNumber().Data()),"NCellsAbove100VsModule",200 , 0, 200,
                                                 fNMaxDCalModules, nModulesStart, fNMaxDCalModules+nModulesStart);
      fHistNCellsBigger100MeVvsMod->GetXaxis()->SetTitle("N_{cells} with E_{cell} > 0.1 GeV");
      fHistNCellsBigger100MeVvsMod->GetYaxis()->SetTitle("module ID");
      fHistExtQA->Add(fHistNCellsBigger100MeVvsMod);
      fHistNCellsBigger1500MeVvsMod   = new TH2F(Form("NCellsAbove1500VsModule %s",GetCutNumber().Data()),"NCellsAbove1500VsModule", 100, 0, 100,
                                                 fNMaxDCalModules,nModulesStart,nModulesStart+fNMaxDCalModules);
      fHistNCellsBigger1500MeVvsMod->GetXaxis()->SetTitle("N_{cells} with E_{cell} > 1.5 GeV");
      fHistNCellsBigger1500MeVvsMod->GetYaxis()->SetTitle("module ID");
      fHistExtQA->Add(fHistNCellsBigger1500MeVvsMod);
      fHistEnergyOfModvsMod           = new TH2F(Form("ModuleEnergyVsModule %s",GetCutNumber().Data()),"ModuleEnergyVsModule", nBinsClusterEMod, arrClusEBinning,
                                                 fNMaxDCalModules, nModulesStart, fNMaxDCalModules+nModulesStart);
      fHistEnergyOfModvsMod->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistEnergyOfModvsMod->GetYaxis()->SetTitle("N_{cells}");
      fHistExtQA->Add(fHistEnergyOfModvsMod);
      fHistClusterEnergyvsNCellsBeforeQA      = new TH2F(Form("ClusterEnergyVsNCells_beforeQA %s",GetCutNumber().Data()),"ClusterEnergyVsNCells_beforeQA",nBinsClusterE, arrClusEBinning,
                                                50, 0, 50);
      fHistClusterEnergyvsNCellsBeforeQA->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistClusterEnergyvsNCellsBeforeQA->GetYaxis()->SetTitle("N_{cells}");
      fHistExtQA->Add(fHistClusterEnergyvsNCellsBeforeQA);
      fHistClusterEnergyvsNCellsAfterQA      = new TH2F(Form("ClusterEnergyVsNCells_afterQA %s",GetCutNumber().Data()),"ClusterEnergyVsNCells_afterQA",nBinsClusterE, arrClusEBinning,
                                                 50, 0, 50);
      fHistClusterEnergyvsNCellsAfterQA->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistClusterEnergyvsNCellsAfterQA->GetYaxis()->SetTitle("N_{cells}");
      fHistExtQA->Add(fHistClusterEnergyvsNCellsAfterQA);
      fHistClusterDistanceInTimeCut   = new TH2F(Form("ClusterDistanceTo_withinTimingCut %s",GetCutNumber().Data()),"ClusterDistanceTo_withinTimingCut",20,-10,10,20,-10,10);
      fHistClusterDistanceInTimeCut->GetXaxis()->SetTitle("R_{cl,row} within time cut (cell)");
      fHistClusterDistanceInTimeCut->GetYaxis()->SetTitle("R_{cl,col} within time cut (cell)");
      fHistExtQA->Add(fHistClusterDistanceInTimeCut);
      fHistClusterDistanceOutTimeCut  = new TH2F(Form("ClusterDistanceTo_outsideTimingCut %s",GetCutNumber().Data()),"ClusterDistanceTo_outsideTimingCut",20,-10,10,20,-10,10);
      fHistClusterDistanceOutTimeCut->GetXaxis()->SetTitle("R_{cl,row} outside time cut (cell)");
      fHistClusterDistanceOutTimeCut->GetYaxis()->SetTitle("R_{cl,col} outside time cut (cell)");
      fHistExtQA->Add(fHistClusterDistanceOutTimeCut);
      fHistClusterDistance1DInTimeCut   = new TH1F(Form("Cluster1D_DistanceTo_withinTimingCut %s",GetCutNumber().Data()),"Cluster1D_DistanceTo_withinTimingCut",200,0.,0.5);
      fHistClusterDistance1DInTimeCut->GetXaxis()->SetTitle("R_{cl,1D} within time cut (cell)");
      fHistExtQA->Add(fHistClusterDistance1DInTimeCut);

      // detailed cell QA histos for DCAL
      if(fExtendedMatchAndQA > 3){
        fHistCellEnergyvsCellID         = new TH2F(Form("CellEnergyVsCellID %s",GetCutNumber().Data()),"CellEnergyVsCellID",nBinsClusterECellCoarse,  arrClusEBinningCoarse,
                                                   nMaxCellsDCAL, nCellsStart, nMaxCellsDCAL+nCellsStart);
        fHistCellEnergyvsCellID->GetXaxis()->SetTitle("E_{cell} (GeV)");
        fHistCellEnergyvsCellID->GetYaxis()->SetTitle("Cell ID");
        fHistExtQA->Add(fHistCellEnergyvsCellID);
        fHistCellTimevsCellID           = new TH2F(Form("CellTimeVsCellID %s",GetCutNumber().Data()),"CellTimeVsCellID",600,-timeMax,timeMax, nMaxCellsDCAL, nCellsStart,
                                                   nMaxCellsDCAL+nCellsStart);
        fHistCellTimevsCellID->GetXaxis()->SetTitle("t_{cell} (GeV)");
        fHistCellTimevsCellID->GetYaxis()->SetTitle("Cell ID");
        fHistExtQA->Add(fHistCellTimevsCellID);
        fHistClusterIncludedCellsBeforeQA       = new TH1F(Form("ClusterIncludedCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_beforeClusterQA",
                                                           nMaxCellsDCAL, nCellsStart, nMaxCellsDCAL+nCellsStart);
        fHistClusterIncludedCellsBeforeQA->GetXaxis()->SetTitle("Cell ID");
        fHistClusterIncludedCellsBeforeQA->GetYaxis()->SetTitle("# included in cluster");
        fHistExtQA->Add(fHistClusterIncludedCellsBeforeQA);
        fHistClusterIncludedCellsAfterQA        = new TH1F(Form("ClusterIncludedCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_afterClusterQA",
                                                           nMaxCellsDCAL, nCellsStart, nMaxCellsDCAL+nCellsStart);
        fHistClusterIncludedCellsAfterQA->GetXaxis()->SetTitle("Cell ID");
        fHistClusterIncludedCellsAfterQA->GetYaxis()->SetTitle("# included in cluster");
        fHistExtQA->Add(fHistClusterIncludedCellsAfterQA);
        fHistClusterEnergyFracCellsBeforeQA     = new TH1F(Form("ClusterEnergyFracCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_beforeClusterQA",
                                                           nMaxCellsDCAL, nCellsStart, nMaxCellsDCAL+nCellsStart);
        fHistClusterEnergyFracCellsBeforeQA->GetXaxis()->SetTitle("Cell ID");
        fHistClusterEnergyFracCellsBeforeQA->GetYaxis()->SetTitle("fraction of energy cell is contrib. to cl.");
        fHistClusterEnergyFracCellsBeforeQA->Sumw2();
        fHistExtQA->Add(fHistClusterEnergyFracCellsBeforeQA);
        fHistClusterEnergyFracCellsAfterQA      = new TH1F(Form("ClusterEnergyFracCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_afterClusterQA",
                                                           nMaxCellsDCAL, nCellsStart, nMaxCellsDCAL+nCellsStart);
        fHistClusterEnergyFracCellsAfterQA->GetXaxis()->SetTitle("Cell ID");
        fHistClusterEnergyFracCellsAfterQA->GetYaxis()->SetTitle("fraction of energy cell is contrib. to cl.");
        fHistClusterEnergyFracCellsAfterQA->Sumw2();
        fHistExtQA->Add(fHistClusterEnergyFracCellsAfterQA);
        fHistClusterIncludedCellsTimingAfterQA  = new TH2F(Form("ClusterIncludedCellsTiming_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCellsTiming_afterClusterQA",
                                                           nBinsClusterE, arrClusEBinning, 200, -500, 500);
        fHistClusterIncludedCellsTimingAfterQA->GetXaxis()->SetTitle("E_{cl} (GeV)");
        fHistClusterIncludedCellsTimingAfterQA->GetYaxis()->SetTitle("t_{cell} (ns)");
        fHistExtQA->Add(fHistClusterIncludedCellsTimingAfterQA);
        fHistClusterIncludedCellsTimingEnergyAfterQA  = new TH2F(Form("ClusterIncludedCellsTimingEnergy_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCellsTimingEnergy_afterClusterQA",
                                                                 nBinsClusterECell,  arrClusEBinningCoarse, 200, -500, 500);
        fHistClusterIncludedCellsTimingEnergyAfterQA->GetXaxis()->SetTitle("E_{cell} (GeV)");
        fHistClusterIncludedCellsTimingEnergyAfterQA->GetYaxis()->SetTitle("t_{cell} (ns)");
        fHistExtQA->Add(fHistClusterIncludedCellsTimingEnergyAfterQA);
      }
    } else if( fClusterType == 4 ){ //EMCAL+DCAL
      Int_t nMaxCellsEMCAL           = fNMaxEMCalModules*48*24;
      Int_t nMaxCellsDCAL            = fNMaxDCalModules*32*24;

      fHistClusterEnergyvsMod         = new TH2F(Form("ClusterEnergyVsModule_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyVsModule_afterClusterQA",
                                                 nBinsClusterE, arrClusEBinning, fNMaxDCalModules+fNMaxEMCalModules, 0, fNMaxDCalModules+fNMaxEMCalModules);
      fHistClusterEnergyvsMod->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistClusterEnergyvsMod->GetYaxis()->SetTitle("module ID");
      fHistExtQA->Add(fHistClusterEnergyvsMod);
      fHistNCellsBigger100MeVvsMod    = new TH2F(Form("NCellsAbove100VsModule %s",GetCutNumber().Data()),"NCellsAbove100VsModule",200 , 0, 200,
                                                 fNMaxDCalModules+fNMaxEMCalModules, 0, fNMaxDCalModules+fNMaxEMCalModules);
      fHistNCellsBigger100MeVvsMod->GetXaxis()->SetTitle("N_{cells} with E_{cell} > 0.1 GeV");
      fHistNCellsBigger100MeVvsMod->GetYaxis()->SetTitle("module ID");
      fHistExtQA->Add(fHistNCellsBigger100MeVvsMod);
      fHistNCellsBigger1500MeVvsMod   = new TH2F(Form("NCellsAbove1500VsModule %s",GetCutNumber().Data()),"NCellsAbove1500VsModule", 100, 0, 100,
                                                 fNMaxDCalModules+fNMaxEMCalModules, 0, fNMaxDCalModules+fNMaxEMCalModules);
      fHistNCellsBigger1500MeVvsMod->GetXaxis()->SetTitle("N_{cells} with E_{cell} > 1.5 GeV");
      fHistNCellsBigger1500MeVvsMod->GetYaxis()->SetTitle("module ID");
      fHistExtQA->Add(fHistNCellsBigger1500MeVvsMod);
      fHistEnergyOfModvsMod           = new TH2F(Form("ModuleEnergyVsModule %s",GetCutNumber().Data()),"ModuleEnergyVsModule", nBinsClusterEMod, arrClusEBinning,
                                                 fNMaxDCalModules+fNMaxEMCalModules, 0, fNMaxDCalModules+fNMaxEMCalModules);
      fHistEnergyOfModvsMod->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistEnergyOfModvsMod->GetYaxis()->SetTitle("N_{cells}");
      fHistExtQA->Add(fHistEnergyOfModvsMod);
      fHistClusterEnergyvsNCellsBeforeQA      = new TH2F(Form("ClusterEnergyVsNCells_beforeQA %s",GetCutNumber().Data()),"ClusterEnergyVsNCells_beforeQA",nBinsClusterE, arrClusEBinning,
                                                50, 0, 50);
      fHistClusterEnergyvsNCellsBeforeQA->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistClusterEnergyvsNCellsBeforeQA->GetYaxis()->SetTitle("N_{cells}");
      fHistExtQA->Add(fHistClusterEnergyvsNCellsBeforeQA);
      fHistClusterEnergyvsNCellsAfterQA      = new TH2F(Form("ClusterEnergyVsNCells_afterQA %s",GetCutNumber().Data()),"ClusterEnergyVsNCells_afterQA",nBinsClusterE, arrClusEBinning,
                                                 50, 0, 50);
      fHistClusterEnergyvsNCellsAfterQA->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistClusterEnergyvsNCellsAfterQA->GetYaxis()->SetTitle("N_{cells}");
      fHistExtQA->Add(fHistClusterEnergyvsNCellsAfterQA);
      fHistClusterDistanceInTimeCut   = new TH2F(Form("ClusterDistanceTo_withinTimingCut %s",GetCutNumber().Data()),"ClusterDistanceTo_withinTimingCut",20,-10,10,20,-10,10);
      fHistClusterDistanceInTimeCut->GetXaxis()->SetTitle("R_{cl,row} within time cut (cell)");
      fHistClusterDistanceInTimeCut->GetYaxis()->SetTitle("R_{cl,col} within time cut (cell)");
      fHistExtQA->Add(fHistClusterDistanceInTimeCut);
      fHistClusterDistanceOutTimeCut  = new TH2F(Form("ClusterDistanceTo_outsideTimingCut %s",GetCutNumber().Data()),"ClusterDistanceTo_outsideTimingCut",20,-10,10,20,-10,10);
      fHistClusterDistanceOutTimeCut->GetXaxis()->SetTitle("R_{cl,row} outside time cut (cell)");
      fHistClusterDistanceOutTimeCut->GetYaxis()->SetTitle("R_{cl,col} outside time cut (cell)");
      fHistExtQA->Add(fHistClusterDistanceOutTimeCut);
      fHistClusterDistance1DInTimeCut   = new TH1F(Form("Cluster1D_DistanceTo_withinTimingCut %s",GetCutNumber().Data()),"Cluster1D_DistanceTo_withinTimingCut",200,0.,0.5);
      fHistClusterDistance1DInTimeCut->GetXaxis()->SetTitle("R_{cl,1D} within time cut (cell)");
      fHistExtQA->Add(fHistClusterDistance1DInTimeCut);

      // detailed cell QA histos for DCAL
      if(fExtendedMatchAndQA > 3){
        fHistCellEnergyvsCellID         = new TH2F(Form("CellEnergyVsCellID %s",GetCutNumber().Data()),"CellEnergyVsCellID",nBinsClusterECellCoarse,  arrClusEBinningCoarse,
                                                   nMaxCellsEMCAL+nMaxCellsDCAL, 0, nMaxCellsEMCAL+nMaxCellsDCAL);
        fHistCellEnergyvsCellID->GetXaxis()->SetTitle("E_{cell} (GeV)");
        fHistCellEnergyvsCellID->GetYaxis()->SetTitle("Cell ID");
        fHistExtQA->Add(fHistCellEnergyvsCellID);
        fHistCellTimevsCellID           = new TH2F(Form("CellTimeVsCellID %s",GetCutNumber().Data()),"CellTimeVsCellID",600,-timeMax,timeMax,
                                                   nMaxCellsEMCAL+nMaxCellsDCAL, 0, nMaxCellsEMCAL+nMaxCellsDCAL);
        fHistCellTimevsCellID->GetXaxis()->SetTitle("t_{cell} (GeV)");
        fHistCellTimevsCellID->GetYaxis()->SetTitle("Cell ID");
        fHistExtQA->Add(fHistCellTimevsCellID);
        fHistClusterIncludedCellsBeforeQA       = new TH1F(Form("ClusterIncludedCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_beforeClusterQA",
                                                           nMaxCellsEMCAL+nMaxCellsDCAL, 0, nMaxCellsEMCAL+nMaxCellsDCAL);
        fHistClusterIncludedCellsBeforeQA->GetXaxis()->SetTitle("Cell ID");
        fHistClusterIncludedCellsBeforeQA->GetYaxis()->SetTitle("# included in cluster");
        fHistExtQA->Add(fHistClusterIncludedCellsBeforeQA);
        fHistClusterIncludedCellsAfterQA        = new TH1F(Form("ClusterIncludedCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCells_afterClusterQA",
                                                           nMaxCellsEMCAL+nMaxCellsDCAL, 0, nMaxCellsEMCAL+nMaxCellsDCAL);
        fHistClusterIncludedCellsAfterQA->GetXaxis()->SetTitle("Cell ID");
        fHistClusterIncludedCellsAfterQA->GetYaxis()->SetTitle("# included in cluster");
        fHistExtQA->Add(fHistClusterIncludedCellsAfterQA);
        fHistClusterEnergyFracCellsBeforeQA     = new TH1F(Form("ClusterEnergyFracCells_beforeClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_beforeClusterQA",
                                                           nMaxCellsEMCAL+nMaxCellsDCAL, 0, nMaxCellsEMCAL+nMaxCellsDCAL);
        fHistClusterEnergyFracCellsBeforeQA->GetXaxis()->SetTitle("Cell ID");
        fHistClusterEnergyFracCellsBeforeQA->GetYaxis()->SetTitle("fraction of energy cell is contrib. to cl.");
        fHistClusterEnergyFracCellsBeforeQA->Sumw2();
        fHistExtQA->Add(fHistClusterEnergyFracCellsBeforeQA);
        fHistClusterEnergyFracCellsAfterQA      = new TH1F(Form("ClusterEnergyFracCells_afterClusterQA %s",GetCutNumber().Data()),"ClusterEnergyFracCells_afterClusterQA",
                                                           nMaxCellsEMCAL+nMaxCellsDCAL, 0, nMaxCellsEMCAL+nMaxCellsDCAL);
        fHistClusterEnergyFracCellsAfterQA->GetXaxis()->SetTitle("Cell ID");
        fHistClusterEnergyFracCellsAfterQA->GetYaxis()->SetTitle("fraction of energy cell is contrib. to cl.");
        fHistClusterEnergyFracCellsAfterQA->Sumw2();
        fHistExtQA->Add(fHistClusterEnergyFracCellsAfterQA);
        fHistClusterIncludedCellsTimingAfterQA  = new TH2F(Form("ClusterIncludedCellsTiming_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCellsTiming_afterClusterQA",
                                                           nBinsClusterE, arrClusEBinning, 200, -500, 500);
        fHistClusterIncludedCellsTimingAfterQA->GetXaxis()->SetTitle("E_{cl} (GeV)");
        fHistClusterIncludedCellsTimingAfterQA->GetYaxis()->SetTitle("t_{cell} (ns)");
        fHistExtQA->Add(fHistClusterIncludedCellsTimingAfterQA);
        fHistClusterIncludedCellsTimingEnergyAfterQA  = new TH2F(Form("ClusterIncludedCellsTimingEnergy_afterClusterQA %s",GetCutNumber().Data()),"ClusterIncludedCellsTimingEnergy_afterClusterQA",
                                                                 nBinsClusterECell,  arrClusEBinningCoarse, 200, -500, 500);
        fHistClusterIncludedCellsTimingEnergyAfterQA->GetXaxis()->SetTitle("E_{cell} (GeV)");
        fHistClusterIncludedCellsTimingEnergyAfterQA->GetYaxis()->SetTitle("t_{cell} (ns)");
        fHistExtQA->Add(fHistClusterIncludedCellsTimingEnergyAfterQA);
      }
    } else{AliError(Form("fExtendedMatchAndQA (%i) not (yet) defined for cluster type (%i)",fExtendedMatchAndQA,fClusterType));}
  }

  //TrackMatching histograms
  if(fUseDistTrackToCluster && !fDoLightOutput) {
    const Int_t nEtaBins        = 300;
    const Int_t nPhiBins        = 300;
    const Float_t EtaRange[2]   = {-0.3, 0.3};
    const Float_t PhiRange[2]   = {-0.3, 0.3};

    fHistClusterdEtadPhiBeforeQA    = new TH2F(Form("dEtaVsdPhi_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_beforeClusterQA", nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
    fHistograms->Add(fHistClusterdEtadPhiBeforeQA);
    fHistClusterdEtadPhiAfterQA     = new TH2F(Form("dEtaVsdPhi_afterClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_afterClusterQA", nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
    fHistograms->Add(fHistClusterdEtadPhiAfterQA);
    //----------------
    if(fIsMC > 1){
      fHistClusterdEtadPhiBeforeQA->Sumw2();
      fHistClusterdEtadPhiAfterQA->Sumw2();
    }
    //----------------
    if(fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5 ){

      // QA histograms for track matching
      fHistClusterRBeforeQA                           = new TH1F(Form("R_Cluster_beforeClusterQA %s",GetCutNumber().Data()),"R of cluster",200,400,500);
      fHistClusterRBeforeQA->GetXaxis()->SetTitle("R_{prim vtx} (m)");
      fHistograms->Add(fHistClusterRBeforeQA);
      fHistClusterRAfterQA                            = new TH1F(Form("R_Cluster_afterClusterQA %s",GetCutNumber().Data()),"R of cluster_matched",200,400,500);
      fHistClusterRAfterQA->GetXaxis()->SetTitle("R_{prim vtx} (m)");
      fHistograms->Add(fHistClusterRAfterQA);
      fHistDistanceTrackToClusterBeforeQA             = new TH1F(Form("DistanceToTrack_beforeClusterQA %s",GetCutNumber().Data()),"DistanceToTrack_beforeClusterQA",200,0,2);
      fHistDistanceTrackToClusterBeforeQA->GetXaxis()->SetTitle("R_{cl-track}");
      fHistograms->Add(fHistDistanceTrackToClusterBeforeQA);
      fHistDistanceTrackToClusterAfterQA              = new TH1F(Form("DistanceToTrack_afterClusterQA %s",GetCutNumber().Data()),"DistanceToTrack_afterClusterQA",200,0,2);
      fHistDistanceTrackToClusterAfterQA->GetXaxis()->SetTitle("R_{cl-track}");
      fHistograms->Add(fHistDistanceTrackToClusterAfterQA);
      fHistClusterdEtadPhiPosTracksBeforeQA           = new TH2F(Form("dEtaVsdPhi_posTracks_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_posTracks_beforeClusterQA",
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistDistanceTrackToClusterAfterQA->GetXaxis()->SetTitle("d#eta_{cl-track}");
      fHistDistanceTrackToClusterAfterQA->GetYaxis()->SetTitle("d#varphi_{cl-track} (rad)");
      fHistograms->Add(fHistClusterdEtadPhiPosTracksBeforeQA);
      fHistClusterdEtadPhiNegTracksBeforeQA           = new TH2F(Form("dEtaVsdPhi_negTracks_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_negTracks_beforeClusterQA",
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistClusterdEtadPhiNegTracksBeforeQA->GetXaxis()->SetTitle("d#eta_{cl-track}");
      fHistClusterdEtadPhiNegTracksBeforeQA->GetYaxis()->SetTitle("d#varphi_{cl-track} (rad)");
      fHistograms->Add(fHistClusterdEtadPhiNegTracksBeforeQA);
      fHistClusterdEtadPhiPosTracksAfterQA            = new TH2F(Form("dEtaVsdPhi_posTracks_afterClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_posTracks_afterClusterQA",
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistClusterdEtadPhiPosTracksAfterQA->GetXaxis()->SetTitle("d#eta_{cl-track}");
      fHistClusterdEtadPhiPosTracksAfterQA->GetYaxis()->SetTitle("d#varphi_{cl-track} (rad)");
      fHistograms->Add(fHistClusterdEtadPhiPosTracksAfterQA);
      fHistClusterdEtadPhiNegTracksAfterQA            = new TH2F(Form("dEtaVsdPhi_negTracks_afterClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_negTracks_afterClusterQA",
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistClusterdEtadPhiNegTracksAfterQA->GetXaxis()->SetTitle("d#eta_{cl-track}");
      fHistClusterdEtadPhiNegTracksAfterQA->GetYaxis()->SetTitle("d#varphi_{cl-track} (rad)");
      fHistograms->Add(fHistClusterdEtadPhiNegTracksAfterQA);
      fHistClusterdEtadPhiPosTracksP_000_075BeforeQA  = new TH2F(Form("dEtaVsdPhi_posTracks_P<0.75_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_posTracks_P<0.75_beforeClusterQA",
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistClusterdEtadPhiPosTracksP_000_075BeforeQA->GetXaxis()->SetTitle("d#eta_{cl-track}");
      fHistClusterdEtadPhiPosTracksP_000_075BeforeQA->GetYaxis()->SetTitle("d#varphi_{cl-track} (rad)");
      fHistograms->Add(fHistClusterdEtadPhiPosTracksP_000_075BeforeQA);
      fHistClusterdEtadPhiPosTracksP_075_125BeforeQA  = new TH2F(Form("dEtaVsdPhi_posTracks_0.75<P<1.25_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_posTracks_0.75<P<1.25_beforeClusterQA",
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistClusterdEtadPhiPosTracksP_075_125BeforeQA->GetXaxis()->SetTitle("d#eta_{cl-track}");
      fHistClusterdEtadPhiPosTracksP_075_125BeforeQA->GetYaxis()->SetTitle("d#varphi_{cl-track} (rad)");
      fHistograms->Add(fHistClusterdEtadPhiPosTracksP_075_125BeforeQA);
      fHistClusterdEtadPhiPosTracksP_125_999BeforeQA  = new TH2F(Form("dEtaVsdPhi_posTracks_P>1.25_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_posTracks_P>1.25_beforeClusterQA",
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistClusterdEtadPhiPosTracksP_125_999BeforeQA->GetXaxis()->SetTitle("d#eta_{cl-track}");
      fHistClusterdEtadPhiPosTracksP_125_999BeforeQA->GetYaxis()->SetTitle("d#varphi_{cl-track} (rad)");
      fHistograms->Add(fHistClusterdEtadPhiPosTracksP_125_999BeforeQA);
      fHistClusterdEtadPhiNegTracksP_000_075BeforeQA  = new TH2F(Form("dEtaVsdPhi_negTracks_P<0.75_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_negTracks_P<0.75_beforeClusterQA",
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistClusterdEtadPhiNegTracksP_000_075BeforeQA->GetXaxis()->SetTitle("d#eta_{cl-track}");
      fHistClusterdEtadPhiNegTracksP_000_075BeforeQA->GetYaxis()->SetTitle("d#varphi_{cl-track} (rad)");
      fHistograms->Add(fHistClusterdEtadPhiNegTracksP_000_075BeforeQA);
      fHistClusterdEtadPhiNegTracksP_075_125BeforeQA  = new TH2F(Form("dEtaVsdPhi_negTracks_0.75<P<1.25_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_negTracks_0.75<P<1.25_beforeClusterQA",
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistClusterdEtadPhiNegTracksP_075_125BeforeQA->GetXaxis()->SetTitle("d#eta_{cl-track}");
      fHistClusterdEtadPhiNegTracksP_075_125BeforeQA->GetYaxis()->SetTitle("d#varphi_{cl-track} (rad)");
      fHistograms->Add(fHistClusterdEtadPhiNegTracksP_075_125BeforeQA);
      fHistClusterdEtadPhiNegTracksP_125_999BeforeQA  = new TH2F(Form("dEtaVsdPhi_negTracks_P>1.25_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsdPhi_negTracks_P>1.25_beforeClusterQA",
                                                                 nEtaBins, EtaRange[0], EtaRange[1], nPhiBins, PhiRange[0], PhiRange[1]);
      fHistClusterdEtadPhiNegTracksP_125_999BeforeQA->GetXaxis()->SetTitle("d#eta_{cl-track}");
      fHistClusterdEtadPhiNegTracksP_125_999BeforeQA->GetYaxis()->SetTitle("d#varphi_{cl-track} (rad)");
      fHistograms->Add(fHistClusterdEtadPhiNegTracksP_125_999BeforeQA);
      fHistClusterdEtadPtBeforeQA                     = new TH2F(Form("dEtaVsPt_beforeClusterQA %s",GetCutNumber().Data()),"dEtaVsPt_beforeClusterQA",nEtaBins,EtaRange[0],EtaRange[1],
                                                                 nBinsClusterEMod, arrClusEBinning);
      fHistClusterdEtadPtBeforeQA->GetXaxis()->SetTitle("d#eta_{cl-track}");
      fHistClusterdEtadPtBeforeQA->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      fHistograms->Add(fHistClusterdEtadPtBeforeQA);
      fHistClusterdEtadPtAfterQA                      = new TH2F(Form("dEtaVsPt_afterClusterQA %s",GetCutNumber().Data()),"dEtaVsPt_afterClusterQA",nEtaBins,EtaRange[0],EtaRange[1],
                                                                 nBinsClusterEMod, arrClusEBinning);
      fHistClusterdEtadPtAfterQA->GetXaxis()->SetTitle("d#eta_{cl-track}");
      fHistClusterdEtadPtAfterQA->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      fHistograms->Add(fHistClusterdEtadPtAfterQA);
      fHistClusterdPhidPtPosTracksBeforeQA            = new TH2F(Form("dPhiVsPt_posTracks_beforeClusterQA %s",GetCutNumber().Data()),"dPhiVsPt_posTracks_beforeClusterQA",
                                                                 2*nPhiBins, 2*PhiRange[0], 2*PhiRange[1], nBinsClusterEMod, arrClusEBinning);
      fHistClusterdPhidPtPosTracksBeforeQA->GetXaxis()->SetTitle("d#varphi_{cl-track} (rad)");
      fHistClusterdPhidPtPosTracksBeforeQA->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      fHistograms->Add(fHistClusterdPhidPtPosTracksBeforeQA);
      fHistClusterdPhidPtNegTracksBeforeQA            = new TH2F(Form("dPhiVsPt_negTracks_beforeClusterQA %s",GetCutNumber().Data()),"dPhiVsPt_negTracks_beforeClusterQA",
                                                                 2*nPhiBins, 2*PhiRange[0], 2*PhiRange[1], nBinsClusterEMod, arrClusEBinning);
      fHistClusterdPhidPtNegTracksBeforeQA->GetXaxis()->SetTitle("d#varphi_{cl-track} (rad)");
      fHistClusterdPhidPtNegTracksBeforeQA->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      fHistograms->Add(fHistClusterdPhidPtNegTracksBeforeQA);
      fHistClusterdPhidPtAfterQA                      = new TH2F(Form("dPhiVsPt_afterClusterQA %s",GetCutNumber().Data()),"dPhiVsPt_afterClusterQA",2*nPhiBins,2*PhiRange[0],2*PhiRange[1],
                                                                 nBinsClusterEMod, arrClusEBinning);
      fHistClusterdPhidPtAfterQA->GetXaxis()->SetTitle("d#varphi_{cl-track} (rad)");
      fHistClusterdPhidPtAfterQA->GetYaxis()->SetTitle("p_{T} (GeV/c)");
      fHistograms->Add(fHistClusterdPhidPtAfterQA);

      if(fIsMC > 0 && fIsPureCalo == 0){ // these histograms are so far only used in conjunction with PCM, namely in MatchConvPhotonToCluster
        fHistClusterdEtadPtTrueMatched                = new TH2F(Form("dEtaVsPt_TrueMatched %s",GetCutNumber().Data()),"dEtaVsPt_TrueMatched",nEtaBins,EtaRange[0],EtaRange[1],
                                                                  nBinsClusterEMod, arrClusEBinning);
        fHistClusterdEtadPtTrueMatched->GetXaxis()->SetTitle("d#eta_{cl-track}");
        fHistClusterdEtadPtTrueMatched->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fHistograms->Add(fHistClusterdEtadPtTrueMatched);
        fHistClusterdPhidPtPosTracksTrueMatched       = new TH2F(Form("dPhiVsPt_posTracks_TrueMatched %s",GetCutNumber().Data()),"dPhiVsPt_posTracks_TrueMatched",
                                                                 2*nPhiBins,2*PhiRange[0],2*PhiRange[1], nBinsClusterEMod, arrClusEBinning);
        fHistClusterdPhidPtPosTracksTrueMatched->GetXaxis()->SetTitle("d#varphi_{cl-track} (rad)");
        fHistClusterdPhidPtPosTracksTrueMatched->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fHistograms->Add(fHistClusterdPhidPtPosTracksTrueMatched);
        fHistClusterdPhidPtNegTracksTrueMatched       = new TH2F(Form("dPhiVsPt_negTracks_TrueMatched %s",GetCutNumber().Data()),"dPhiVsPt_negTracks_TrueMatched",
                                                                 2*nPhiBins,2*PhiRange[0],2*PhiRange[1], nBinsClusterEMod, arrClusEBinning);
        fHistClusterdPhidPtNegTracksTrueMatched->GetXaxis()->SetTitle("d#varphi_{cl-track} (rad)");
        fHistClusterdPhidPtNegTracksTrueMatched->GetYaxis()->SetTitle("p_{T} (GeV/c)");
        fHistograms->Add(fHistClusterdPhidPtNegTracksTrueMatched);
      }

      if(fUseEOverPVetoTM){
        // plot trackP vs. clusterE in case of a match
        fHistMatchedTrackPClusE  = new TH2F(Form("MatchedTrackPClusE %s",GetCutNumber().Data()), "Matched tracks",
                                            nBinsClusterEOnlyHighPt, arrClusEBinningOnlyHighPt, nBinsClusterEOnlyHighPt, arrClusEBinningOnlyHighPt);
        fHistMatchedTrackPClusE->GetXaxis()->SetTitle("E_{cl} (GeV)");
        fHistMatchedTrackPClusE->GetYaxis()->SetTitle("P_{track} (GeV/c)");
        fHistograms->Add(fHistMatchedTrackPClusE);

        // plot trackP vs. clusterE in case of a match AFTER EOverP veto has been applied
        fHistMatchedTrackPClusEAfterEOverPVeto = new TH2F(Form("MatchedTrackPClusEAfterEOverPVeto  %s",GetCutNumber().Data()), "Matched tracks after EOverP veto",
                                                        nBinsClusterEOnlyHighPt, arrClusEBinningOnlyHighPt, nBinsClusterEOnlyHighPt, arrClusEBinningOnlyHighPt);
        fHistMatchedTrackPClusEAfterEOverPVeto ->GetXaxis()->SetTitle("E_{cl} (GeV)");
        fHistMatchedTrackPClusEAfterEOverPVeto ->GetYaxis()->SetTitle("P_{track} (GeV/c)");
        fHistograms->Add(fHistMatchedTrackPClusEAfterEOverPVeto);
      }
      //----------------
      if(fIsMC > 1){
        fHistClusterRBeforeQA->Sumw2();
        fHistClusterRAfterQA->Sumw2();
        fHistDistanceTrackToClusterBeforeQA->Sumw2();
        fHistDistanceTrackToClusterAfterQA->Sumw2();
        fHistClusterdEtadPhiPosTracksBeforeQA->Sumw2();
        fHistClusterdEtadPhiNegTracksBeforeQA->Sumw2();
        fHistClusterdEtadPhiPosTracksAfterQA->Sumw2();
        fHistClusterdEtadPhiNegTracksAfterQA->Sumw2();
        fHistClusterdEtadPhiPosTracksP_000_075BeforeQA->Sumw2();
        fHistClusterdEtadPhiPosTracksP_075_125BeforeQA->Sumw2();
        fHistClusterdEtadPhiPosTracksP_125_999BeforeQA->Sumw2();
        fHistClusterdEtadPhiNegTracksP_000_075BeforeQA->Sumw2();
        fHistClusterdEtadPhiNegTracksP_075_125BeforeQA->Sumw2();
        fHistClusterdEtadPhiNegTracksP_125_999BeforeQA->Sumw2();
        fHistClusterdEtadPtBeforeQA->Sumw2();
        fHistClusterdEtadPtAfterQA->Sumw2();
        fHistClusterdPhidPtPosTracksBeforeQA->Sumw2();
        fHistClusterdPhidPtNegTracksBeforeQA->Sumw2();
        fHistClusterdPhidPtAfterQA->Sumw2();
        if(fIsPureCalo == 0){
          fHistClusterdEtadPtTrueMatched->Sumw2();
          fHistClusterdPhidPtPosTracksTrueMatched->Sumw2();
          fHistClusterdPhidPtNegTracksTrueMatched->Sumw2();
        }

        if(fUseEOverPVetoTM){
          fHistMatchedTrackPClusE->Sumw2();
          fHistMatchedTrackPClusEAfterEOverPVeto->Sumw2();
        }
      }
      //----------------
    }
  }
  if( fUseDistTrackToCluster && fIsMC && (fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5 )){
    // TM efficiency histograms
    const Int_t nEmcalEtaBins             = 96;
    const Int_t nEmcalPhiBins             = 124;
    Float_t EmcalEtaBins[nEmcalEtaBins+1] = {-0.66687,-0.653,-0.63913,-0.62526,-0.61139,-0.59752,-0.58365,-0.56978,-0.55591,-0.54204,-0.52817,-0.5143,-0.50043,-0.48656,-0.47269,-0.45882,-0.44495,-0.43108,-0.41721,-0.40334,-0.38947,-0.3756,-0.36173,-0.34786,-0.33399,-0.32012,-0.30625,-0.29238,-0.27851,-0.26464,-0.25077,-0.2369,-0.22303,-0.20916,-0.19529,-0.18142,-0.16755,-0.15368,-0.13981,-0.12594,-0.11207,-0.0982,-0.08433,-0.07046,-0.05659,-0.04272,-0.02885,-0.01498,-0.00111,0.01276,0.02663,0.0405,0.05437,0.06824,0.08211,0.09598,0.10985,0.12372,0.13759,0.15146,0.16533,0.1792,0.19307,0.20694,0.22081,0.23468,0.24855,0.26242,0.27629,0.29016,0.30403,0.3179,0.33177,0.34564,0.35951,0.37338,0.38725,0.40112,0.41499,0.42886,0.44273,0.4566,0.47047,0.48434,0.49821,0.51208,0.52595,0.53982,0.55369,0.56756,0.58143,0.5953,0.60917,0.62304,0.63691,0.65078,0.66465};
    Float_t EmcalPhiBins[nEmcalPhiBins+1] = {1.408,1.4215,1.435,1.4485,1.462,1.4755,1.489,1.5025,1.516,1.5295,1.543,1.5565,1.57,1.5835,1.597,1.6105,1.624,1.6375,1.651,1.6645,1.678,1.6915,1.705,1.7185,1.732,1.758,1.7715,1.785,1.7985,1.812,1.8255,1.839,1.8525,1.866,1.8795,1.893,1.9065,1.92,1.9335,1.947,1.9605,1.974,1.9875,2.001,2.0145,2.028,2.0415,2.055,2.0685,2.082,2.108,2.1215,2.135,2.1485,2.162,2.1755,2.189,2.2025,2.216,2.2295,2.243,2.2565,2.27,2.2835,2.297,2.3105,2.324,2.3375,2.351,2.3645,2.378,2.3915,2.405,2.4185,2.432,2.456,2.4695,2.483,2.4965,2.51,2.5235,2.537,2.5505,2.564,2.5775,2.591,2.6045,2.618,2.6315,2.645,2.6585,2.672,2.6855,2.699,2.7125,2.726,2.7395,2.753,2.7665,2.78,2.804,2.8175,2.831,2.8445,2.858,2.8715,2.885,2.8985,2.912,2.9255,2.939,2.9525,2.966,2.9795,2.993,3.0065,3.02,3.0335,3.047,3.0605,3.074,3.0875,3.101,3.1145,3.128};

    fHistClusterTrueElecEtaPhiBeforeTM_30_00  = new TH2F(Form("ElecEtaPhiBeforeTM_clusterE>30_Histo %s",GetCutNumber().Data()),"electron clusters before track matching, E_{cl}>30GeV",nEmcalPhiBins,EmcalPhiBins,nEmcalEtaBins,EmcalEtaBins);
    fHistClusterTrueElecEtaPhiBeforeTM_30_00->GetXaxis()->SetTitle("#varphi (rad)");
    fHistClusterTrueElecEtaPhiBeforeTM_30_00->GetYaxis()->SetTitle("#eta");
    fHistograms->Add(fHistClusterTrueElecEtaPhiBeforeTM_30_00);


    fHistClusterTrueElecEtaPhiAfterTM_30_00  = new TH2F(Form("ElecEtaPhiAfterTM_clusterE>30_Histo %s",GetCutNumber().Data()),"electron clusters before track matching, E_{cl}>30GeV",nEmcalPhiBins,EmcalPhiBins,nEmcalEtaBins,EmcalEtaBins);
    fHistClusterTrueElecEtaPhiAfterTM_30_00->GetXaxis()->SetTitle("#varphi (rad)");
    fHistClusterTrueElecEtaPhiAfterTM_30_00->GetYaxis()->SetTitle("#eta");
    fHistograms->Add(fHistClusterTrueElecEtaPhiAfterTM_30_00);

    fHistClusterTMEffiInput                       = new TH2F(Form("TMEffiInputHisto %s",GetCutNumber().Data()),"TMEffiInputHisto",nBinsClusterE, arrClusEBinning, 22, -0.5, 21.5);
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(1,"All cl");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(2,"Ch cl");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(3,"Ne cl");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(4,"Ne cl sub ch");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(5,"Ga cl");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(6,"Ga cl sub ch");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(7,"conv cl");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(8,"Ch cl prim");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(9,"El cl");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(10,"All cl match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(11,"Ch cl match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(12,"Ch cl match w lead");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(13,"Ne cl match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(14,"Ne cl sub ch match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(15,"Ga cl match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(16,"Ga cl sub ch match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(17,"conv cl match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(18,"conv cl match w lead");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(19,"Ch cl match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(20,"Ch cl match w lead");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(21,"El cl match");
    fHistClusterTMEffiInput->GetYaxis()->SetBinLabel(22,"El cl match w lead");
    fHistClusterTMEffiInput->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistograms->Add(fHistClusterTMEffiInput);

    fHistClusterEvsTrackECharged                  = new TH2F(Form("ClusterE_TrackE_ChargedCluster %s",GetCutNumber().Data()),"ClusterE TrackE ChargedCluster",
                                                             nBinsClusterE, arrClusEBinning, nBinsClusterE, arrClusEBinning);
    fHistClusterEvsTrackECharged->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackECharged->GetYaxis()->SetTitle("#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackECharged);
    fHistClusterEvsTrackEChargedLead              = new TH2F(Form("ClusterE_TrackE_ChargedCluster_LeadMatched %s",GetCutNumber().Data()),"ClusterE TrackE ChargedCluster LeadMatched",
                                                             nBinsClusterE, arrClusEBinning, nBinsClusterE, arrClusEBinning);
    fHistClusterEvsTrackEChargedLead->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackEChargedLead->GetYaxis()->SetTitle("#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackEChargedLead);
    fHistClusterEvsTrackENeutral                  = new TH2F(Form("ClusterE_TrackE_NeutralCluster %s",GetCutNumber().Data()),"ClusterE TrackE NeutralCluster",
                                                             nBinsClusterE, arrClusEBinning, nBinsClusterE, arrClusEBinning);
    fHistClusterEvsTrackENeutral->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackENeutral->GetYaxis()->SetTitle("#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackENeutral);
    fHistClusterEvsTrackENeutralSubCharged        = new TH2F(Form("ClusterE_TrackE_NeutralClusterSubCharged %s",GetCutNumber().Data()),"ClusterE TrackE NeutralCluster SubCharged",
                                                             nBinsClusterE, arrClusEBinning, nBinsClusterE, arrClusEBinning);
    fHistClusterEvsTrackENeutralSubCharged->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackENeutralSubCharged->GetYaxis()->SetTitle("#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackENeutralSubCharged);
    fHistClusterEvsTrackEGamma                    = new TH2F(Form("ClusterE_TrackE_GammaCluster %s",GetCutNumber().Data()),"ClusterE TrackE GammaCluster",
                                                             nBinsClusterE, arrClusEBinning, nBinsClusterE, arrClusEBinning);
    fHistClusterEvsTrackEGamma->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackEGamma->GetYaxis()->SetTitle("#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackEGamma);
    fHistClusterEvsTrackEGammaSubCharged          = new TH2F(Form("ClusterE_TrackE_GammaClusterSubCharged %s",GetCutNumber().Data()),"ClusterE TrackE GammaCluster SubCharged",
                                                             nBinsClusterE, arrClusEBinning, nBinsClusterE, arrClusEBinning);
    fHistClusterEvsTrackEGammaSubCharged->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackEGammaSubCharged->GetYaxis()->SetTitle("#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackEGammaSubCharged);
    fHistClusterEvsTrackEConv                    = new TH2F(Form("ClusterE_TrackE_ConvCluster %s",GetCutNumber().Data()),"ClusterE TrackE ConvCluster",
                                                             nBinsClusterE, arrClusEBinning, nBinsClusterE, arrClusEBinning);
    fHistClusterEvsTrackEConv->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackEConv->GetYaxis()->SetTitle("#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackEConv);

    fHistClusterENMatchesNeutral                  = new TH2F(Form("ClusterE_NMatches_NeutralCluster %s",GetCutNumber().Data()),"ClusterE NMatches NeutralCluster",
                                                             nBinsClusterE, arrClusEBinning, 20, -0.5, 19.5);
    fHistClusterENMatchesNeutral->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterENMatchesNeutral->GetYaxis()->SetTitle("#it{N}_{matches}");
    fHistograms->Add(fHistClusterENMatchesNeutral);
    fHistClusterENMatchesCharged                  = new TH2F(Form("ClusterE_NMatches_ChargedCluster %s",GetCutNumber().Data()),"ClusterE NMatches ChargedCluster",
                                                             nBinsClusterE, arrClusEBinning, 20, -0.5, 19.5);
    fHistClusterENMatchesCharged->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterENMatchesCharged->GetYaxis()->SetTitle("#it{N}_{matches}");
    fHistograms->Add(fHistClusterENMatchesCharged);

    fHistClusterEvsTrackEPrimaryButNoElec        = new TH2F(Form("ClusterE_TrackE_ChargedClusterNoElec %s",GetCutNumber().Data()),"ClusterE TrackE ChargedClusterNoElec",
                                                             nBinsClusterE, arrClusEBinning, nBinsClusterE, arrClusEBinning);
    fHistClusterEvsTrackEPrimaryButNoElec->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackEPrimaryButNoElec->GetYaxis()->SetTitle("#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackEPrimaryButNoElec);
    fHistClusterEvsTrackSumEPrimaryButNoElec        = new TH2F(Form("ClusterE_TrackSumE_ChargedClusterNoElec %s",GetCutNumber().Data()),"ClusterE TrackSumE ChargedClusterNoElec",
                                                             nBinsClusterE, arrClusEBinning, nBinsClusterE, arrClusEBinning);
    fHistClusterEvsTrackSumEPrimaryButNoElec->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
    fHistClusterEvsTrackSumEPrimaryButNoElec->GetYaxis()->SetTitle("#sum#it{E}_{tr} (GeV)");
    fHistograms->Add(fHistClusterEvsTrackSumEPrimaryButNoElec);

    if(fUseEOverPVetoTM){
      fHistClusETruePi0_BeforeTM = new TH1F(Form("ClusETruePi0_BeforeTM %s",GetCutNumber().Data()), "true #pi^{0} clusters before TM",
                                            nBinsClusterEOnlyHighPt, arrClusEBinningOnlyHighPt);
      fHistClusETruePi0_BeforeTM->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
      fHistClusETruePi0_BeforeTM->GetYaxis()->SetTitle("weighted counts");
      fHistograms->Add(fHistClusETruePi0_BeforeTM);

      fHistClusETruePi0_Matched = new TH1F(Form("ClusETruePi0_Matched %s",GetCutNumber().Data()), "matched true #pi^{0} clusters",
                                           nBinsClusterEOnlyHighPt, arrClusEBinningOnlyHighPt);
      fHistClusETruePi0_Matched->GetXaxis()->SetTitle("#it{E}_{cl} (GeV)");
      fHistClusETruePi0_Matched->GetYaxis()->SetTitle("weighted counts");
      fHistograms->Add(fHistClusETruePi0_Matched);

      // plot trackP vs. clusterE in case of a match with pi0
      fHistMatchedTrackPClusETruePi0Clus  = new TH2F(Form("MatchedTrackPClusETruePi0Clus %s",GetCutNumber().Data()), "Matched tracks to #pi^{0} clusters",
                                                      nBinsClusterEOnlyHighPt, arrClusEBinningOnlyHighPt, nBinsClusterEOnlyHighPt, arrClusEBinningOnlyHighPt);
      fHistMatchedTrackPClusETruePi0Clus->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistMatchedTrackPClusETruePi0Clus->GetYaxis()->SetTitle("P_{track} (GeV/c)");
      fHistograms->Add(fHistMatchedTrackPClusETruePi0Clus);
    }


    if(fIsMC > 1){
      fHistClusterTMEffiInput->Sumw2();
      fHistClusterTrueElecEtaPhiBeforeTM_30_00->Sumw2();
      fHistClusterTrueElecEtaPhiAfterTM_30_00->Sumw2();
      fHistClusterEvsTrackECharged->Sumw2();
      fHistClusterEvsTrackEChargedLead->Sumw2();
      fHistClusterEvsTrackENeutral->Sumw2();
      fHistClusterEvsTrackENeutralSubCharged->Sumw2();
      fHistClusterEvsTrackEGamma->Sumw2();
      fHistClusterEvsTrackEGammaSubCharged->Sumw2();
      fHistClusterEvsTrackEConv->Sumw2();
      fHistClusterENMatchesNeutral->Sumw2();
      fHistClusterENMatchesCharged->Sumw2();
      fHistClusterEvsTrackEPrimaryButNoElec->Sumw2();
      fHistClusterEvsTrackSumEPrimaryButNoElec->Sumw2();
      if(fUseEOverPVetoTM){
        fHistMatchedTrackPClusETruePi0Clus->Sumw2();
        fHistClusETruePi0_BeforeTM->Sumw2();
        fHistClusETruePi0_Matched->Sumw2();
      }
    }
  }

  if(fUseElectronClusterCalibration){
    // Propagate electrons to EMCal and match tracks with clusters to compare their energy for calibration

    fHistElectronPositronClusterMatch = new TH2F(Form("MatchedElectronPositronTrackPClusE %s",GetCutNumber().Data()), "Matched Electron Positron tracks with P on EMC",
                                                      nBinsClusterE, arrClusEBinning, nBinsClusterE, arrClusEBinning);
    fHistElectronPositronClusterMatch->GetXaxis()->SetTitle("E_{cl} (GeV)");
    fHistElectronPositronClusterMatch->GetYaxis()->SetTitle("P_{track, EMC}  (GeV/c)");
    fHistograms->Add(fHistElectronPositronClusterMatch);

    fHistElectronPositronClusterMatchSub = new TH2F(Form("MatchedElectronPositronESubP %s",GetCutNumber().Data()), "Matched Electron Positron tracks with P on EMC E Sub P",
                                                      nBinsClusterE, arrClusEBinning, 1000, -5, 5);
    fHistElectronPositronClusterMatchSub->GetXaxis()->SetTitle("E_{cl} (GeV)");
    fHistElectronPositronClusterMatchSub->GetYaxis()->SetTitle("E_{cl} - P_{track, EMC}");
    fHistograms->Add(fHistElectronPositronClusterMatchSub);

    if(fExtendedMatchAndQA > 1 ){
      fHistElectronClusterMatch = new TH2F(Form("MatchedElectronTrackPClusE %s",GetCutNumber().Data()), "Matched Electron tracks with P on EMC",
                                                nBinsClusterE, arrClusEBinning, nBinsClusterE, arrClusEBinning);
      fHistElectronClusterMatch->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistElectronClusterMatch->GetYaxis()->SetTitle("P_{track, EMC} (GeV/c)");
      fHistograms->Add(fHistElectronClusterMatch);

      fHistPositronClusterMatch = new TH2F(Form("MatchedPositronTrackPClusE %s",GetCutNumber().Data()), "Matched Positron tracks with P on EMC",
                                                nBinsClusterE, arrClusEBinning, nBinsClusterE, arrClusEBinning);
      fHistPositronClusterMatch->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistPositronClusterMatch->GetYaxis()->SetTitle("P_{track, EMC} (GeV/c)");
      fHistograms->Add(fHistPositronClusterMatch);
    }


    if(fIsMC > 0){
      fHistTrueElectronPositronClusterMatch = new TH2F(Form("TrueMatchedElectronPositronTrackPClusE %s",GetCutNumber().Data()), "True Matched Electron Positron tracks",
                                                  nBinsClusterE, arrClusEBinning, nBinsClusterE, arrClusEBinning);
      fHistTrueElectronPositronClusterMatch->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistTrueElectronPositronClusterMatch->GetYaxis()->SetTitle("P_{track} (GeV/c)");
      fHistograms->Add(fHistTrueElectronPositronClusterMatch);

      fHistTrueNoElectronPositronClusterMatch = new TH2F(Form("TrueMatchedNoElectronPositronTrackPClusE %s",GetCutNumber().Data()), "True Matched No Electron Positron tracks",
                                                  nBinsClusterE, arrClusEBinning, nBinsClusterE, arrClusEBinning);
      fHistTrueNoElectronPositronClusterMatch->GetXaxis()->SetTitle("E_{cl} (GeV)");
      fHistTrueNoElectronPositronClusterMatch->GetYaxis()->SetTitle("P_{track} (GeV/c)");
      fHistograms->Add(fHistTrueNoElectronPositronClusterMatch);


      fHistElectronClusterMatchTruePID = new TH2F(Form("MatchedElectronPositronTrackPID %s",GetCutNumber().Data()),"Matched Electron Positron tracks true PID",
                                                  5,0,5, nBinsClusterE, arrClusEBinning);
      fHistElectronClusterMatchTruePID->GetXaxis()->SetTitle("true PID");
      fHistElectronClusterMatchTruePID->GetYaxis()->SetTitle("P_{track} (GeV/c)");
      fHistElectronClusterMatchTruePID->GetXaxis()->SetBinLabel(1, "e^{#pm}");
      fHistElectronClusterMatchTruePID->GetXaxis()->SetBinLabel(2, "#pi^{#pm}");
      fHistElectronClusterMatchTruePID->GetXaxis()->SetBinLabel(3, "proton");
      fHistElectronClusterMatchTruePID->GetXaxis()->SetBinLabel(4, "K^{#pm}");
      fHistElectronClusterMatchTruePID->GetXaxis()->SetBinLabel(5, "rest");
      fHistograms->Add(fHistElectronClusterMatchTruePID);
    }

    if(fIsMC > 1){
      fHistElectronPositronClusterMatch->Sumw2();
      fHistElectronPositronClusterMatchSub->Sumw2();
      fHistTrueElectronPositronClusterMatch->Sumw2();
      fHistTrueNoElectronPositronClusterMatch->Sumw2();
      fHistElectronClusterMatchTruePID->Sumw2();
      if(fExtendedMatchAndQA > 1 ){
        fHistElectronClusterMatch->Sumw2();
        fHistPositronClusterMatch->Sumw2();
      }
    }

  }


  if (fDoExoticsQA){
    if( fClusterType == 1 ){ //EMCAL
      const Int_t nEmcalEtaBins             = 96;
      const Int_t nEmcalPhiBins             = 124;
      Float_t EmcalEtaBins[nEmcalEtaBins+1] = {-0.66687,-0.653,-0.63913,-0.62526,-0.61139,-0.59752,-0.58365,-0.56978,-0.55591,-0.54204,-0.52817,-0.5143,-0.50043,-0.48656,-0.47269,-0.45882,-0.44495,-0.43108,-0.41721,-0.40334,-0.38947,-0.3756,-0.36173,-0.34786,-0.33399,-0.32012,-0.30625,-0.29238,-0.27851,-0.26464,-0.25077,-0.2369,-0.22303,-0.20916,-0.19529,-0.18142,-0.16755,-0.15368,-0.13981,-0.12594,-0.11207,-0.0982,-0.08433,-0.07046,-0.05659,-0.04272,-0.02885,-0.01498,-0.00111,0.01276,0.02663,0.0405,0.05437,0.06824,0.08211,0.09598,0.10985,0.12372,0.13759,0.15146,0.16533,0.1792,0.19307,0.20694,0.22081,0.23468,0.24855,0.26242,0.27629,0.29016,0.30403,0.3179,0.33177,0.34564,0.35951,0.37338,0.38725,0.40112,0.41499,0.42886,0.44273,0.4566,0.47047,0.48434,0.49821,0.51208,0.52595,0.53982,0.55369,0.56756,0.58143,0.5953,0.60917,0.62304,0.63691,0.65078,0.66465};
      Float_t EmcalPhiBins[nEmcalPhiBins+1] = {1.408,1.4215,1.435,1.4485,1.462,1.4755,1.489,1.5025,1.516,1.5295,1.543,1.5565,1.57,1.5835,1.597,1.6105,1.624,1.6375,1.651,1.6645,1.678,1.6915,1.705,1.7185,1.732,1.758,1.7715,1.785,1.7985,1.812,1.8255,1.839,1.8525,1.866,1.8795,1.893,1.9065,1.92,1.9335,1.947,1.9605,1.974,1.9875,2.001,2.0145,2.028,2.0415,2.055,2.0685,2.082,2.108,2.1215,2.135,2.1485,2.162,2.1755,2.189,2.2025,2.216,2.2295,2.243,2.2565,2.27,2.2835,2.297,2.3105,2.324,2.3375,2.351,2.3645,2.378,2.3915,2.405,2.4185,2.432,2.456,2.4695,2.483,2.4965,2.51,2.5235,2.537,2.5505,2.564,2.5775,2.591,2.6045,2.618,2.6315,2.645,2.6585,2.672,2.6855,2.699,2.7125,2.726,2.7395,2.753,2.7665,2.78,2.804,2.8175,2.831,2.8445,2.858,2.8715,2.885,2.8985,2.912,2.9255,2.939,2.9525,2.966,2.9795,2.993,3.0065,3.02,3.0335,3.047,3.0605,3.074,3.0875,3.101,3.1145,3.128};

      fHistClusterEtavsPhiExotics     = new TH2F(Form("EtaPhi_Exotics %s",GetCutNumber().Data()),"EtaPhi_Exotics",nEmcalPhiBins, EmcalPhiBins, nEmcalEtaBins, EmcalEtaBins);
      fHistClusterEtavsPhiExotics->GetXaxis()->SetTitle("#eta");
      fHistClusterEtavsPhiExotics->GetYaxis()->SetTitle("#varphi (rad)");
      fHistograms->Add(fHistClusterEtavsPhiExotics);
    } else if( fClusterType == 2 ){ //PHOS
      const Int_t nPhosEtaBins        = 56;
      const Int_t nPhosPhiBins        = 192;
      const Float_t PhosEtaRange[2]   = {-0.16, 0.16};
      const Float_t PhosPhiRange[2]   = {4.5, 5.6};

      fHistClusterEtavsPhiExotics     = new TH2F(Form("EtaPhi_Exotics %s",GetCutNumber().Data()),"EtaPhi_Exotics",nPhosPhiBins, PhosPhiRange[0], PhosPhiRange[1],
                                                 nPhosEtaBins, PhosEtaRange[0], PhosEtaRange[1]);
      fHistClusterEtavsPhiExotics->GetXaxis()->SetTitle("#eta");
      fHistClusterEtavsPhiExotics->GetYaxis()->SetTitle("#varphi (rad)");
      fHistograms->Add(fHistClusterEtavsPhiExotics);
    } else if( fClusterType == 3 ){ //DCAL
      const Int_t nDcalEtaBins             = 96;
      const Int_t nDcalPhiBins             = 124;

      fHistClusterEtavsPhiExotics     = new TH2F(Form("EtaPhi_Exotics %s",GetCutNumber().Data()),"EtaPhi_Exotics",nDcalPhiBins,4.5,5.7,nDcalEtaBins,-0.66687,0.66465);
      fHistClusterEtavsPhiExotics->GetXaxis()->SetTitle("#eta");
      fHistClusterEtavsPhiExotics->GetYaxis()->SetTitle("#varphi (rad)");
      fHistograms->Add(fHistClusterEtavsPhiExotics);
    } else if( fClusterType == 4 ){ //EMCAL+DCAL
      fHistClusterEtavsPhiExotics     = new TH2F(Form("EtaPhi_Exotics %s",GetCutNumber().Data()),"EtaPhi_Exotics",462,0,2*TMath::Pi(),110,-0.7,0.7);
      fHistClusterEtavsPhiExotics->GetXaxis()->SetTitle("#eta");
      fHistClusterEtavsPhiExotics->GetYaxis()->SetTitle("#varphi (rad)");
      fHistograms->Add(fHistClusterEtavsPhiExotics);
    }
    fHistClusterEM02Exotics           = new TH2F(Form("EVsM02_Exotics %s",GetCutNumber().Data()),"EVsM02_afterClusterQA",nBinsClusterE, arrClusEBinning,200,0,2);
    fHistClusterEM02Exotics->GetXaxis()->SetTitle("E_{cl,exotics} (GeV)");
    fHistClusterEM02Exotics->GetYaxis()->SetTitle("#sigma_{long}^2");
    fHistograms->Add(fHistClusterEM02Exotics);
    fHistClusterEnergyvsNCellsExotics = new TH2F(Form("ClusterEnergyVsNCells_Exotics %s",GetCutNumber().Data()),"ClusterEnergyVsNCells_Exotics",
                                                 nBinsClusterE, arrClusEBinning, 50, 0, 50);
    fHistClusterEnergyvsNCellsExotics->GetXaxis()->SetTitle("E_{cl,exotics} (GeV)");
    fHistClusterEnergyvsNCellsExotics->GetYaxis()->SetTitle("N_{LM}");
    fHistograms->Add(fHistClusterEnergyvsNCellsExotics);
    fHistClusterEEstarExotics         = new TH2F(Form("ClusterEnergyVsEnergystar_Exotics %s",GetCutNumber().Data()),"ClusterEnergyVsEnergystar_Exotics",
                                                 nBinsClusterE, arrClusEBinning, nBinsClusterE, arrClusEBinning);
    fHistClusterEEstarExotics->GetXaxis()->SetTitle("E_{cl,exotics} (GeV)");
    fHistClusterEEstarExotics->GetYaxis()->SetTitle("E_{cl,#star} (GeV)");
    fHistograms->Add(fHistClusterEEstarExotics);

    if (fIsMC > 1){
      fHistClusterEtavsPhiExotics->Sumw2();
      fHistClusterEM02Exotics->Sumw2();
      fHistClusterEnergyvsNCellsExotics->Sumw2();
      fHistClusterEEstarExotics->Sumw2();
    }
  }

  fVectorMatchedClusterIDs.clear();

  TH1::AddDirectory(kTRUE);
  return;
}

//________________________________________________________________________
void AliCaloPhotonCuts::InitializeEMCAL(AliVEvent *event){

  if (fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
    AliTender* alitender=0x0;
    AliEmcalTenderTask* emcaltender=0x0;
    AliEmcalCorrectionTask* emcalCorrTask=0x0;
    if(event->IsA()==AliESDEvent::Class()){
      alitender         = (AliTender*) AliAnalysisManager::GetAnalysisManager()->GetTopTasks()->FindObject("AliTender");
      if(!alitender){
        emcalCorrTask  = (AliEmcalCorrectionTask*) AliAnalysisManager::GetAnalysisManager()->GetTopTasks()->FindObject("AliEmcalCorrectionTask_defaultSetting");
      }
    } else if( event->IsA()==AliAODEvent::Class()){
      emcaltender       = (AliEmcalTenderTask*) AliAnalysisManager::GetAnalysisManager()->GetTopTasks()->FindObject("AliEmcalTenderTask");
      if(!emcaltender)
        emcalCorrTask  = (AliEmcalCorrectionTask*) AliAnalysisManager::GetAnalysisManager()->GetTopTasks()->FindObject("AliEmcalCorrectionTask_defaultSetting");
    }
    if(alitender){
      TIter next(alitender->GetSupplies());
      AliTenderSupply *supply;
      while ((supply=(AliTenderSupply*)next())) if(supply->IsA()==AliEMCALTenderSupply::Class()) break;
      fEMCALRecUtils        = ((AliEMCALTenderSupply*)supply)->GetRecoUtils();
      fEMCALBadChannelsMap  = fEMCALRecUtils->GetEMCALBadChannelStatusMapArray();
    } else if(emcaltender){
      fEMCALRecUtils        = ((AliEMCALTenderSupply*)emcaltender->GetEMCALTenderSupply())->GetRecoUtils();
      fEMCALBadChannelsMap  = fEMCALRecUtils->GetEMCALBadChannelStatusMapArray();
    } else if(emcalCorrTask){
      AliEmcalCorrectionComponent * emcalCorrComponent = emcalCorrTask->GetCorrectionComponent("AliEmcalCorrectionCellBadChannel_defaultSetting");
      if(emcalCorrComponent){
        fEMCALRecUtils        = emcalCorrComponent->GetRecoUtils();
        fEMCALBadChannelsMap  = fEMCALRecUtils->GetEMCALBadChannelStatusMapArray();
        fEMCALBadChannelsMap1D  = fEMCALRecUtils->GetEMCALChannelStatusMap1D();
        if(!fEMCALBadChannelsMap1D || fEMCALBadChannelsMap1D->GetNbinsX()<1e3)
          fEMCALBadChannelsMap1D = NULL;
      }
    }
    if (fEMCALRecUtils) fEMCALInitialized = kTRUE;

    fGeomEMCAL              = AliEMCALGeometry::GetInstance();
    if(!fGeomEMCAL){ AliFatal("EMCal geometry not initialized!");}


    //retrieve pointer to CaloIsolation Instance
    if(fUsePhotonIsolation) fCaloIsolation = (AliPhotonIsolation*)AliAnalysisManager::GetAnalysisManager()->GetTask(fCaloIsolationName.Data());
    if(!fCaloIsolation && fUsePhotonIsolation){ AliFatal("AliPhotonIsolation instance could not be initialized!");}

    //retrieve pointer to trackMatcher Instance
    if(fUseDistTrackToCluster || fUseElectronClusterCalibration) fCaloTrackMatcher = (AliCaloTrackMatcher*)AliAnalysisManager::GetAnalysisManager()->GetTask(fCaloTrackMatcherName.Data());
    if(!fCaloTrackMatcher && ( fUseDistTrackToCluster || fUseElectronClusterCalibration ) ){ AliFatal("CaloTrackMatcher instance could not be initialized!");}

    if(!fDoLightOutput || fDoFlatEnergySubtraction){
      Int_t nMaxCellsEMCAL  = fNMaxEMCalModules*48*24;
      Int_t nMinCellsDCAL = 12288;
      Int_t nMaxCellsDCAL = nMinCellsDCAL+fNMaxDCalModules*32*24;
      Int_t nMaxCells = 0;
      Int_t nMinCells = 0;
      if(fClusterType == 1){
        nMaxCells = nMaxCellsEMCAL;
        nMinCells = 0;
      } else if(fClusterType == 3){
        nMaxCells = nMaxCellsDCAL;
        nMinCells = nMinCellsDCAL;
      } else if(fClusterType == 4){
        nMaxCells = nMaxCellsEMCAL+nMaxCellsDCAL;
        nMinCells = 0;
      }


      Int_t imod = -1;Int_t iTower = -1, iIphi = -1, iIeta = -1;
      Int_t icol = -1;Int_t irow = -1;

      fNactiveEmcalCells = 0;
      if(fEMCALBadChannelsMap1D){
        for(Int_t iCell=nMinCells;iCell<nMaxCells;iCell++){
          Int_t iBadCell      = (Int_t) fEMCALBadChannelsMap1D->GetBinContent(iCell);
          if(iBadCell > 0) fBadChannels->Fill(iCell,1);
          else { fBadChannels->Fill(iCell,0); fNactiveEmcalCells++; }
        }
      } else {
        for(Int_t iCell=nMinCells;iCell<nMaxCells;iCell++){
          fGeomEMCAL->GetCellIndex(iCell,imod,iTower,iIphi,iIeta);
          if (fEMCALBadChannelsMap->GetEntries() <= imod) continue;
          fGeomEMCAL->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi,iIeta,irow,icol);
          Int_t iBadCell      = (Int_t) ((TH2I*)fEMCALBadChannelsMap->At(imod))->GetBinContent(icol,irow);
          if(iBadCell > 0) fBadChannels->Fill(iCell,1);
          else { fBadChannels->Fill(iCell,0); fNactiveEmcalCells++; }
        }
      }
    }
  }
  return;
}

//________________________________________________________________________
void AliCaloPhotonCuts::InitializePHOS (AliVEvent *event){

  if (fClusterType == 2){
    if(!fDoLightOutput){
      fGeomPHOS = AliPHOSGeometry::GetInstance();
      if(!fGeomPHOS) AliFatal("PHOS geometry not initialized!");
      Int_t nModules = fGeomPHOS->GetNModules();
      //cout << nModules << endl;

      fPHOSBadChannelsMap = new TH2I*[nModules];
      AliPHOSTenderTask* aliphostender = (AliPHOSTenderTask*) AliAnalysisManager::GetAnalysisManager()->GetTopTasks()->FindObject("PHOSTenderTask");
      AliPHOSTenderSupply *PHOSSupply  =((AliPHOSTenderSupply*) aliphostender->GetPHOSTenderSupply()) ;

      if(!PHOSSupply){
        AliError(Form("Can not find PHOSTenderSupply in run %d. No bad channel map could be found for QA!\n",event->GetRunNumber())) ;
        for(Int_t mod=0;mod<nModules;mod++) fPHOSBadChannelsMap[mod] = NULL;
      }else{
        AliInfo("Setting PHOS bad map from PHOSSupply \n") ;
        for(Int_t mod=0;mod<nModules;mod++){
          TH2I * h = (TH2I*)PHOSSupply->GetPHOSBadChannelStatusMap(mod);
          if(h){
              fPHOSBadChannelsMap[mod] = new TH2I(*h);
              AliInfo(Form("using bad map for module %d with nch=%f\n",mod,h->Integral()));
          }
          else fPHOSBadChannelsMap[mod] = NULL;
        }

        Int_t nMaxCellsPHOS = (fNMaxPHOSModules*56*64);
        Int_t relid[4];

        for(Int_t iCell=0;iCell<nMaxCellsPHOS;iCell++){
          fGeomPHOS->AbsToRelNumbering(iCell,relid);
          //cout << relid[0] << ", "  << relid[1] << ", "  << relid[2] << ", "  << relid[3] << ", " << __LINE__ << endl;
          if(relid[1]!=0) AliFatal("PHOS CPV in PHOS cell array?");
          if(relid[0]>=nModules || relid[0]<0 || !fPHOSBadChannelsMap[relid[0]]) continue;
          Int_t iBadCell = (Int_t) ((TH2I*)fPHOSBadChannelsMap[relid[0]])->GetBinContent(relid[2],relid[3]);
          if(iBadCell > 0) fBadChannels->Fill(iCell,1);
          else fBadChannels->Fill(iCell,0);
        }
      }
    }

    if(!fPHOSGeoUtils) fPHOSGeoUtils = new AliPHOSGeoUtils("IHEP","");
    //retrieve pointer to trackMatcher Instance
    if(fUseDistTrackToCluster || fUseElectronClusterCalibration) fCaloTrackMatcher = (AliCaloTrackMatcher*)AliAnalysisManager::GetAnalysisManager()->GetTask(fCaloTrackMatcherName.Data());
    if(!fCaloTrackMatcher && ( fUseDistTrackToCluster || fUseElectronClusterCalibration )){ AliFatal("CaloTrackMatcher instance could not be initialized!");}

    fPHOSInitialized = kTRUE;
    fPHOSCurrentRun = event->GetRunNumber();
  }
  return;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterIsSelectedMC(TParticle *particle,AliMCEvent *mcEvent){
   // MonteCarlo Photon Selection

  if(!mcEvent)return kFALSE;
  if(!particle) return kFALSE;

  if (particle->GetPdgCode() == 22){
    if ( fClusterType == 4){
      //pseudorapidty range same for EMC and DMC
      if ( particle->Eta() < fMinEtaCut || particle->Eta() > fMaxEtaCut ) return kFALSE;
      //check if outside of EMC and DMC phi acceptance
      if ( (particle->Phi() < fMinPhiCut || particle->Phi() > fMaxPhiCut) && (particle->Phi() < fMinPhiCutDMC || particle->Phi() > fMaxPhiCutDMC) ) return kFALSE;
      //if in DMC phi range, reject clusters in hole
      // if ( particle->Phi() > fMinPhiCutDMC && particle->Phi() < fMaxPhiCutDMC && particle->Eta() < fMaxEtaInnerEdge && particle->Eta() > fMinEtaInnerEdge ) return kFALSE;
    } else {
      if ( particle->Eta() < fMinEtaCut || particle->Eta() > fMaxEtaCut ) return kFALSE;
      if ( particle->Phi() < fMinPhiCut || particle->Phi() > fMaxPhiCut ) return kFALSE;
      if ( fClusterType == 3 && particle->Eta() < fMaxEtaInnerEdge && particle->Eta() > fMinEtaInnerEdge ) return kFALSE;
    }
    if(particle->GetMother(0) >-1 && mcEvent->Particle(particle->GetMother(0))->GetPdgCode() == 22){
      return kFALSE;// no photon as mothers!
    }
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterIsSelectedElecMC(TParticle *particle,AliMCEvent *mcEvent){
   // MonteCarlo Photon Selection

  if(!mcEvent)return kFALSE;
  if(!particle) return kFALSE;

  if (TMath::Abs(particle->GetPdgCode()) == 11){

    if ( fClusterType == 4){
      //pseudorapidty range same for EMC and DMC
      if ( particle->Eta() < fMinEtaCut || particle->Eta() > fMaxEtaCut ) return kFALSE;
      //check if outside of EMC and DMC phi acceptance
      if ( (particle->Phi() < fMinPhiCut || particle->Phi() > fMaxPhiCut) && (particle->Phi() < fMinPhiCutDMC || particle->Phi() > fMaxPhiCutDMC) ) return kFALSE;
      //if in DMC phi range, reject clusters in hole
      // if ( particle->Phi() > fMinPhiCutDMC && particle->Phi() < fMaxPhiCutDMC && particle->Eta() < fMaxEtaInnerEdge && particle->Eta() > fMinEtaInnerEdge ) return kFALSE;
    } else {
      if ( particle->Eta() < fMinEtaCut || particle->Eta() > fMaxEtaCut ) return kFALSE;
      if ( particle->Phi() < fMinPhiCut || particle->Phi() > fMaxPhiCut ) return kFALSE;
      if ( fClusterType == 3 && particle->Eta() < fMaxEtaInnerEdge && particle->Eta() > fMinEtaInnerEdge ) return kFALSE;
    }
    if(particle->GetMother(0) >-1 && mcEvent->Particle(particle->GetMother(0))->GetPdgCode() == 11){
      return kFALSE;// no photon as mothers!
    }
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterIsSelectedElecAODMC(AliAODMCParticle *particle,TClonesArray *aodmcArray){
   // MonteCarlo Photon Selection

  if(!aodmcArray)return kFALSE;
  if(!particle) return kFALSE;

  if (TMath::Abs(particle->GetPdgCode()) == 11){

    if ( fClusterType == 4){
      //pseudorapidty range same for EMC and DMC
      if ( particle->Eta() < fMinEtaCut || particle->Eta() > fMaxEtaCut ) return kFALSE;
      //check if outside of EMC and DMC phi acceptance
      if ( (particle->Phi() < fMinPhiCut || particle->Phi() > fMaxPhiCut) && (particle->Phi() < fMinPhiCutDMC || particle->Phi() > fMaxPhiCutDMC) ) return kFALSE;
      //if in DMC phi range, reject clusters in hole
      // if ( particle->Phi() > fMinPhiCutDMC && particle->Phi() < fMaxPhiCutDMC && particle->Eta() < fMaxEtaInnerEdge && particle->Eta() > fMinEtaInnerEdge ) return kFALSE;
    } else {
      if ( particle->Eta() < fMinEtaCut || particle->Eta() > fMaxEtaCut ) return kFALSE;
      if ( particle->Phi() < fMinPhiCut || particle->Phi() > fMaxPhiCut ) return kFALSE;
      if ( fClusterType == 3 && particle->Eta() < fMaxEtaInnerEdge && particle->Eta() > fMinEtaInnerEdge ) return kFALSE;
    }
    if(particle->GetMother() >-1 && (static_cast<AliAODMCParticle*>(aodmcArray->At(particle->GetMother())))->GetPdgCode() == 11){
      return kFALSE;// no photon as mothers!
    }
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterIsSelectedAODMC(AliAODMCParticle *particle,TClonesArray *aodmcArray){
  // MonteCarlo Photon Selection

  if(!aodmcArray)return kFALSE;
  if(!particle) return kFALSE;

  if (particle->GetPdgCode() == 22){
    if ( fClusterType == 4){
      //pseudorapidty range same for EMC and DMC
      if ( particle->Eta() < fMinEtaCut || particle->Eta() > fMaxEtaCut ) return kFALSE;
      //check if outside of EMC and DMC phi acceptance
      if ( (particle->Phi() < fMinPhiCut || particle->Phi() > fMaxPhiCut) && (particle->Phi() < fMinPhiCutDMC || particle->Phi() > fMaxPhiCutDMC) ) return kFALSE;
      //if in DMC phi range, reject clusters in hole
      // if ( particle->Phi() > fMinPhiCutDMC && particle->Phi() < fMaxPhiCutDMC && particle->Eta() < fMaxEtaInnerEdge && particle->Eta() > fMinEtaInnerEdge ) return kFALSE;
    } else {
      if ( particle->Eta() < fMinEtaCut || particle->Eta() > fMaxEtaCut ) return kFALSE;
      if ( particle->Phi() < fMinPhiCut || particle->Phi() > fMaxPhiCut ) return kFALSE;
      if ( fClusterType == 3 && particle->Eta() < fMaxEtaInnerEdge && particle->Eta() > fMinEtaInnerEdge ) return kFALSE;
    }
    if(particle->GetMother() > -1 && (static_cast<AliAODMCParticle*>(aodmcArray->At(particle->GetMother())))->GetPdgCode() == 22){
        return kFALSE;// no photon as mothers!
    }
    return kTRUE;// return in case of accepted gamma
  }
  return kFALSE;
}

//________________________________________________________________________
// This function selects the clusters based on their quality criteria
//________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterQualityCuts(AliVCluster* cluster, AliVEvent *event, AliMCEvent* mcEvent, Int_t isMC, Double_t weight, Long_t clusterID)
{   // Specific Photon Cuts

  // Initialize EMCAL rec utils if not initialized
  if(!fEMCALInitialized && (fClusterType == 1 || fClusterType == 3 || fClusterType == 4)) InitializeEMCAL(event);

  fIsCurrentClusterAcceptedBeforeTM = kFALSE;
  fIsAcceptedForBasic               = kFALSE;
  Int_t cutIndex = 0;
  if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex,cluster->E());
  cutIndex++;

  // cluster position defintion
  Float_t clusPos[3]={0,0,0};
  cluster->GetPosition(clusPos);
  TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
  Double_t etaCluster = clusterVector.Eta();
  Double_t phiCluster = clusterVector.Phi();
  if (phiCluster < 0) phiCluster += 2*TMath::Pi();


  Int_t nLM = GetNumberOfLocalMaxima(cluster, event);
//   Int_t nLMGustavo = fEMCALCaloUtils->GetNumberOfLocalMaxima(cluster, event->GetEMCALCells()) ;
//   cout << "mine: " << nLM << "\t Gustavo: " << nLMGustavo << endl;

  // Fill Histos before Cuts
  if(fHistClusterTimevsEBeforeQA) fHistClusterTimevsEBeforeQA->Fill(cluster->GetTOF(), cluster->E(), weight);
  if(fHistEnergyOfClusterBeforeQA) fHistEnergyOfClusterBeforeQA->Fill(cluster->E(), weight);
  if(fHistNCellsBeforeQA) fHistNCellsBeforeQA->Fill(cluster->GetNCells(), weight);
  if(fHistM02BeforeQA) fHistM02BeforeQA->Fill(cluster->GetM02(), weight);
  if(fHistM20BeforeQA) fHistM20BeforeQA->Fill(cluster->GetM20(), weight);
  if(fHistDispersionBeforeQA) fHistDispersionBeforeQA->Fill(cluster->GetDispersion(), weight);
  if(fHistNLMBeforeQA) fHistNLMBeforeQA->Fill(nLM, weight);
//   if(fHistNLMAvsNLMBBeforeQA) fHistNLMAvsNLMBBeforeQA->Fill(nLM, nLMGustavo, weight);
  if(fHistClusterEM02BeforeQA) fHistClusterEM02BeforeQA->Fill(cluster->E(),cluster->GetM02(), weight);

  AliVCaloCells* cells = NULL;
  if(fExtendedMatchAndQA > 1){

    if(fHistClusterEnergyvsNCellsBeforeQA) fHistClusterEnergyvsNCellsBeforeQA->Fill(cluster->E(),cluster->GetNCells());
    if(cluster->IsEMCAL()){ //EMCAL
      cells = event->GetEMCALCells();
    }else if(cluster->IsPHOS()){ //PHOS
      cells = event->GetPHOSCells();
    }
    if(fHistClusterIncludedCellsBeforeQA){
      Int_t nCellCluster = cluster->GetNCells();
      for(Int_t iCell=0;iCell<nCellCluster;iCell++){
        fHistClusterIncludedCellsBeforeQA->Fill(cluster->GetCellAbsId(iCell));
        if(cluster->E()>0.) fHistClusterEnergyFracCellsBeforeQA->Fill(cluster->GetCellAbsId(iCell),cells->GetCellAmplitude(cluster->GetCellAbsId(iCell))/cluster->E());
      }
    }
  }

  // Check wether timing is ok
  if (fUseTimeDiff){
    if(fUseTimingEfficiencyMCSimCluster==2){
      if ( cluster->E() < 5) {
        if( (cluster->GetTOF() < fMinTimeDiff || cluster->GetTOF() > fMaxTimeDiff) && !(isMC>0)){
          if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//1
          return kFALSE;
        }
      } else {
        if( (cluster->GetTOF() < fMinTimeDiffHighPt || cluster->GetTOF() > fMaxTimeDiffHighPt) && !(isMC>0)){
          if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//1
          return kFALSE;
        }
      }
    } else {
      if( (cluster->GetTOF() < fMinTimeDiff || cluster->GetTOF() > fMaxTimeDiff) && !(isMC>0)){
        if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//1
        return kFALSE;
      }
    }
    if( ((fUseTimingEfficiencyMCSimCluster==1) || (fUseTimingEfficiencyMCSimCluster==2)) && isMC && cluster->E() < fTimingEfficiencyMCSimClusterLowPtEnd && cluster->E() > fMinEnergy ){
      fRandom.SetSeed(0);
      if( fRandom.Uniform(1) > fFuncTimingEfficiencyMCSimCluster->Eval(cluster->E()) ){
        if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//1
        return kFALSE;
      }
    }
    if(cluster->IsPHOS() && fUseTimingEfficiencyMCSimCluster==1 && isMC && cluster->E() > fTimingEfficiencyMCSimClusterHighPtStart){
      fRandom.SetSeed(0);
      if( fRandom.Uniform(1) > fFuncTimingEfficiencyMCSimClusterHighPt->Eval(cluster->E()) ){
        if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//1
        return kFALSE;
      }
    }
  }
  cutIndex++;//2, next cut

  // exotic cluster cut
  Float_t energyStar      = 0;
  if(fUseExoticCluster && IsExoticCluster(cluster, event, energyStar)){
    if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//3
    if (fDoExoticsQA){
      // replay cuts
      Bool_t failed     = kFALSE;
      Bool_t failedM02  = kFALSE;
      Bool_t passedNCellSpecial  = kFALSE;
      if (fUseMinEnergy)
        if(cluster->E() < fMinEnergy)
          failed = kTRUE;
      if (fUseNCells == 1) {
          if (cluster->GetNCells() < fMinNCells)
            failed = kTRUE;
      // special case for PHOS: only apply Ncell cut for clusters with a minimum energy of 1 GeV
      } else if (fUseNCells == 2){
          if (cluster->GetNCells() < fMinNCells && cluster->E() > 1)
            failed = kTRUE;
      // special case for EMCal MC (allow passing of NCell<2 clusters depending on cut efficiency)
      } else if (fUseNCells == 3){
        if(isMC){
          fRandom.SetSeed(0);
          // evaluate effi function and compare to random number between 1 and 2
          // if function value greater than random number, reject cluster. otherwise let it pass
          // function is 1 for E>4 GeV -> will apply standard NCell cut then
          if( (cluster->GetNCells() < fMinNCells) && (fRandom.Uniform(0,1) > fFuncNCellCutEfficiencyEMCal->Eval(cluster->E() )) ){
            failed = kTRUE;
          } else {
            passedNCellSpecial = kTRUE;
          }
        } else {
          if (cluster->GetNCells() < fMinNCells)
            failed = kTRUE;
        }
      } else if (fUseNCells == 4){
        if (cluster->GetNCells() < fMinNCells){
            failed = kTRUE;
        } else if(isMC==0){
          fRandom.SetSeed(0);
          // evaluate effi function and compare to random number between 1 and 2
          // if function value greater than random number, reject cluster. otherwise let it pass
          // function is 1 for E>4 GeV -> will apply standard NCell cut then
          if( (fRandom.Uniform(1,2) < fFuncNCellCutEfficiencyEMCal->Eval(cluster->E() )) ){
            failed = kTRUE;
          }
        }
      }
      if (fUseNLM)
        if( nLM < fMinNLM || nLM > fMaxNLM )
          failed = kTRUE;
      if(!fUseNCells && cluster->GetNCells()<2 && cluster->E()<4){
        // no cut to be applied in this case on M20
        // as cluster needs at least 2 cells for M20 calculation
      } else {
        if (fUseM02 == 1 && !passedNCellSpecial){
          if( cluster->GetM02()< fMinM02 || cluster->GetM02() > fMaxM02 )
            failedM02  = kTRUE;
        } else if (fUseM02 ==2  && !passedNCellSpecial) {
          if( cluster->GetM02()< CalculateMinM02(fMinM02CutNr, cluster->E()) ||
              cluster->GetM02() > CalculateMaxM02(fMaxM02CutNr, cluster->E()) )
            failedM02  = kTRUE;
        // special case for PHOS: only apply M02 cut for clusters with a minimum energy of 1 GeV
        } else if (fUseM02 == 3){
            if( (cluster->GetM02()< fMinM02 || cluster->GetM02() > fMaxM02)  && cluster->E() > 1 )
              failedM02  = kTRUE;
        }
        if (fUseM20 && !passedNCellSpecial)
          if( cluster->GetM20()< fMinM20 || cluster->GetM20() > fMaxM20 )
            failed = kTRUE;
        if (fUseDispersion && !passedNCellSpecial)
          if( cluster->GetDispersion()> fMaxDispersion)
            failed = kTRUE;
      }
      if (fVectorMatchedClusterIDs.size()>0 && fUseDistTrackToCluster){
        if( CheckClusterForTrackMatch(cluster) )
          failed = kTRUE;
      } else if (fUseDistTrackToCluster && fUsePtDepTrackToCluster == 2){
        if (cluster->GetEmcCpvDistance() < fMinTMDistSigma )
          failed = kTRUE;
      }
      if ( !( failed || failedM02 ) ){
        if(fHistClusterEtavsPhiExotics) fHistClusterEtavsPhiExotics->Fill(phiCluster, etaCluster, weight);
        if(fHistClusterEnergyvsNCellsExotics) fHistClusterEnergyvsNCellsExotics->Fill(cluster->E(), cluster->GetNCells(), weight);
        if(fHistClusterEEstarExotics) fHistClusterEEstarExotics->Fill(cluster->E(),energyStar, weight);
      }
      if ( !failed ){
        if(fHistClusterEM02Exotics) fHistClusterEM02Exotics->Fill(cluster->E(), cluster->GetM02(), weight);
      }
    }
    return kFALSE;
  }
  cutIndex++;//3, next cut

  // minimum number of cells
  Bool_t passedSpecialNCell = kFALSE;
  if (fUseNCells == 1){
    if(cluster->GetNCells() < fMinNCells) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//5
      return kFALSE;
    }
  } else if (fUseNCells == 2){ // special case for PHOS: only apply Ncell cut for clusters with a minimum energy of 1 GeV
      if (cluster->GetNCells() < fMinNCells && cluster->E() > 1){
        if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//5
        return kFALSE;
      }
  // special case for EMCal MC (allow passing of NCell<2 clusters depending on cut efficiency)
  } else if (fUseNCells == 3){
    if(isMC>0){
      fRandom.SetSeed(0);
      if(  (cluster->GetNCells() < fMinNCells)){
        // evaluate effi function and compare to random number between 1 and 2
        // if function value greater than random number, reject cluster. otherwise let it pass
        // function is 1 for E>6 GeV -> will apply standard NCell cut then
        if((cluster->E()<6) && (fRandom.Uniform(0,1) < fFuncNCellCutEfficiencyEMCal->Eval(cluster->E()) ) ){
          passedSpecialNCell = kTRUE;
        } else {
          if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//5
          return kFALSE;
        }
      }
    } else {
      if (cluster->GetNCells() < fMinNCells){
        if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//5
        return kFALSE;
      }
    }
  // special case for EMCal Data (rejection of NCell<2 clusters depending on cut efficiency)
  } else if (fUseNCells == 4){
    if (cluster->GetNCells() < fMinNCells){
        if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//5
        return kFALSE;
    } else if(isMC==0){
      fRandom.SetSeed(0);
      // evaluate effi function and compare to random number between 1 and 2
      // if function value greater than random number, reject cluster. otherwise let it pass
      if((cluster->E()<6) && (fRandom.Uniform(1,2) < fFuncNCellCutEfficiencyEMCal->Eval(cluster->E()) ) ){
        if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//5
        return kFALSE;
      }
    }
  }
  cutIndex++;//4, next cut

  // NLM cut
  if (fUseNLM){
    if( nLM < fMinNLM || nLM > fMaxNLM ) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//9
      return kFALSE;
    }
  }
  cutIndex++;//5, next cut

  // M02 cut
  if(passedSpecialNCell || (!fUseNCells && cluster->GetNCells()<2 && cluster->E()<4)){
    // no cut to be applied in this case on M02
    // as cluster needs at least 2 cells for M02 calculation
  } else if (fUseM02 == 1){
    if( cluster->GetM02()< fMinM02 || cluster->GetM02() > fMaxM02 ) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//6
      return kFALSE;
    }
  } else if (fUseM02 ==2 ) {
    if(  cluster->GetM02()< CalculateMinM02(fMinM02CutNr, cluster->E()) ||
      cluster->GetM02() > CalculateMaxM02(fMaxM02CutNr, cluster->E())  ) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//6
      return kFALSE;
    }
  } else if (fUseM02 == 3 && cluster->GetNCells() > 1){ // special case for PHOS: only apply M02 cut for clusters with a minimum energy of 1 GeV
      if( (cluster->GetM02()< fMinM02 || cluster->GetM02() > fMaxM02)  && cluster->E() > 1 ){
        if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//6
        return kFALSE;
      }
  }
  cutIndex++;//6, next cut

  // M20 cut
  if(passedSpecialNCell || (!fUseNCells && cluster->GetNCells()<2 && cluster->E()<4)){
    // no cut to be applied in this case on M20
    // as cluster needs at least 2 cells for M20 calculation
  } else if (fUseM20){
    if( cluster->GetM20()< fMinM20 || cluster->GetM20() > fMaxM20 ) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//7
      return kFALSE;
    }
  }
  cutIndex++;//7, next cut

  // dispersion cut
  if (fUseDispersion && !passedSpecialNCell){
    if( cluster->GetDispersion()> fMaxDispersion) {
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//8
      return kFALSE;
    }
  }
  cutIndex++;//8, next cut


  // Classify events
  AliESDEvent *esdev  = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aodev  = 0;
  Bool_t isESD        = kTRUE;
  if (!esdev) {
    isESD             = kFALSE;
    aodev             = dynamic_cast<AliAODEvent*>(event);
    if (!aodev) {
      AliError("Task needs AOD or ESD event, returning");
      return kFALSE;
    }
  }

  fIsCurrentClusterAcceptedBeforeTM = kTRUE;
  // classification of clusters for TM efficiency
  // 0: Neutral cluster
  // 1: Neutral cluster sub charged
  // 2: Gamma cluster
  // 3: Gamma cluster sub charged
  // 4: Gamma conv cluster
  // 5: Charged cluster
  // 6: Electron
  // 7: prim charged cluster
  Int_t classification  = -1;
  Long_t leadMCLabel    = -1;
  if (isESD)
    leadMCLabel         = ((AliESDCaloCluster*)cluster)->GetLabel();
  else
    leadMCLabel         = ((AliAODCaloCluster*)cluster)->GetLabel();

  // TM efficiency histograms before TM
  if (fIsMC && isMC && (fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5) && fUseDistTrackToCluster  && !(fIsPureCalo > 0 && cluster->E() < 10.)
      && fUsePtDepTrackToCluster < 2){
    classification    = ClassifyClusterForTMEffi(cluster, event, mcEvent, isESD);
    if(IsClusterPi0(event, mcEvent, cluster))
      if(fHistClusETruePi0_BeforeTM) fHistClusETruePi0_BeforeTM->Fill(cluster->E(), weight);
    fHistClusterTMEffiInput->Fill(cluster->E(), 0., weight); //All cl
    if (classification == 5 )
      fHistClusterTMEffiInput->Fill(cluster->E(), 1., weight); //Ch cl
    if (classification == 7 )
      fHistClusterTMEffiInput->Fill(cluster->E(), 7., weight); //Ch cl
    if (classification == 4)
      fHistClusterTMEffiInput->Fill(cluster->E(), 6., weight); //conv electron cl
    if (classification == 6){
      fHistClusterTMEffiInput->Fill(cluster->E(), 8., weight); // electron cl
      if (cluster->E() > 30.)
        fHistClusterTrueElecEtaPhiBeforeTM_30_00->Fill(phiCluster, etaCluster, weight);
    }

    if (classification == 0 || classification == 1)
      fHistClusterTMEffiInput->Fill(cluster->E(), 2., weight); // Ne cl match
    if (classification == 1)
      fHistClusterTMEffiInput->Fill(cluster->E(), 3., weight); // Ne cl sub ch match
    if (classification == 2 || classification == 3)
      fHistClusterTMEffiInput->Fill(cluster->E(), 4., weight); // Ga cl match
    if ( classification == 3)
      fHistClusterTMEffiInput->Fill(cluster->E(), 5., weight); // Ga cl sub ch match

    Int_t nlabelsMatchedTracks      = 0;
    if (fUsePtDepTrackToCluster == 0)
      nlabelsMatchedTracks          = fCaloTrackMatcher->GetNMatchedTrackIDsForCluster(event, cluster->GetID(), fMaxDistTrackToClusterEta, -fMaxDistTrackToClusterEta,
                                                                                      fMaxDistTrackToClusterPhi, fMinDistTrackToClusterPhi);
    else if (fUsePtDepTrackToCluster == 1)
      nlabelsMatchedTracks          = fCaloTrackMatcher->GetNMatchedTrackIDsForCluster(event, cluster->GetID(), fFuncPtDepEta, fFuncPtDepPhi);

    if (classification < 4 && classification > -1)
      fHistClusterENMatchesNeutral->Fill(cluster->E(), nlabelsMatchedTracks);
    else
      fHistClusterENMatchesCharged->Fill(cluster->E(), nlabelsMatchedTracks);

    // plot electrons that survived the track matching
    if (!CheckClusterForTrackMatch(cluster) && classification == 6){ // electrons that survived the matching
      if (cluster->E() > 30.)
        fHistClusterTrueElecEtaPhiAfterTM_30_00->Fill(phiCluster, etaCluster, weight);
    }
  }

  if (fVectorMatchedClusterIDs.size()>0 && fUseDistTrackToCluster && fUsePtDepTrackToCluster < 2){
    if( CheckClusterForTrackMatch(cluster) ){
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//2
      // TM efficiency histos after TM
      if (fIsMC && isMC && (fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5)  && !(fIsPureCalo > 0 && cluster->E() < 10.)){ // ignore low energies for merged analysis
        if(IsClusterPi0(event, mcEvent, cluster))
          if(fHistClusETruePi0_Matched) fHistClusETruePi0_Matched->Fill(cluster->E(), weight);
        fHistClusterTMEffiInput->Fill(cluster->E(), 9., weight);
        if (classification == 5 )
          fHistClusterTMEffiInput->Fill(cluster->E(), 10., weight); //Ch cl match
        if (classification == 4)
          fHistClusterTMEffiInput->Fill(cluster->E(), 16., weight); //conv cl match
        if (classification == 0 || classification == 1)
          fHistClusterTMEffiInput->Fill(cluster->E(), 12., weight); // Ne cl match
        if ( classification == 1)
          fHistClusterTMEffiInput->Fill(cluster->E(), 13., weight); // Ne cl sub ch match
        if (classification == 2 || classification == 3)
          fHistClusterTMEffiInput->Fill(cluster->E(), 14., weight); // Ga cl match
        if ( classification == 3)
          fHistClusterTMEffiInput->Fill(cluster->E(), 15., weight); // Ga cl sub ch match
        if ( classification == 7)
          fHistClusterTMEffiInput->Fill(cluster->E(), 18., weight); // Ch prim cl match
        if ( classification == 6)
          fHistClusterTMEffiInput->Fill(cluster->E(), 20., weight); // El cl match

        vector<Int_t> labelsMatchedTracks;
        if (fUsePtDepTrackToCluster == 0)
          labelsMatchedTracks           = fCaloTrackMatcher->GetMatchedTrackIDsForCluster(event, cluster->GetID(), fMaxDistTrackToClusterEta, -fMaxDistTrackToClusterEta,
                                                                                          fMaxDistTrackToClusterPhi, fMinDistTrackToClusterPhi);
        else if  (fUsePtDepTrackToCluster == 1)
          labelsMatchedTracks           = fCaloTrackMatcher->GetMatchedTrackIDsForCluster(event, cluster->GetID(), fFuncPtDepEta, fFuncPtDepPhi);

        //Int_t idHighestPt = -1;
        Double_t ptMax    = -1;
        Double_t eMax     = -1;
        Double_t eSum     = 0;
        Bool_t foundLead  = kFALSE;
        Double_t eLead    = -1;
        //Int_t idLead      = -1;
        for (Int_t i = 0; i < (Int_t)labelsMatchedTracks.size(); i++){
          AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(labelsMatchedTracks.at(i)));
          eSum += currTrack->E();
          if (ptMax < currTrack->Pt()){
            ptMax               = currTrack->Pt();
            eMax                = currTrack->E();
            //idHighestPt         = labelsMatchedTracks.at(i);
          }
          if (classification == 4 || classification == 5 || classification == 6 || classification == 7){
            Long_t mcLabelTrack = -1;
            if (isESD)
              mcLabelTrack      = TMath::Abs(((AliESDtrack*)currTrack)->GetLabel());
            else
              mcLabelTrack      = TMath::Abs(((AliAODTrack*)currTrack)->GetLabel());
            if (mcLabelTrack!= -1 && mcLabelTrack == leadMCLabel){
              foundLead         = kTRUE;
              eLead             = currTrack->E();
              //idLead            = labelsMatchedTracks.at(i);
            }
          }
        }
        if (classification == 5 || classification == 7 || classification == 6){
          fHistClusterEvsTrackECharged->Fill(cluster->E(), eMax, weight);
          if (classification == 5 || classification == 7){
            fHistClusterEvsTrackEPrimaryButNoElec->Fill(cluster->E(), eMax, weight);
            fHistClusterEvsTrackSumEPrimaryButNoElec->Fill(cluster->E(), eSum, weight);
          }

          if (foundLead ){
            if (classification == 5)
              fHistClusterTMEffiInput->Fill(cluster->E(), 11., weight); //Ch cl match w lead
            if (classification == 7)
              fHistClusterTMEffiInput->Fill(cluster->E(), 19., weight); //Ch prim cl match w lead
            if (classification == 6)
              fHistClusterTMEffiInput->Fill(cluster->E(), 21., weight); //El cl match w lead
            fHistClusterEvsTrackEChargedLead->Fill(cluster->E(), eLead, weight);
          }
        }
        if (classification == 4){
          fHistClusterEvsTrackEConv->Fill(cluster->E(), eMax, weight);
          if (foundLead)
          fHistClusterTMEffiInput->Fill(cluster->E(), 17., weight); //conv cl match w lead
        }
        if (classification == 0 )
          fHistClusterEvsTrackENeutral->Fill(cluster->E(), eMax, weight);
        if (classification == 1)
          fHistClusterEvsTrackENeutralSubCharged->Fill(cluster->E(), eMax, weight);
        if (classification == 2)
          fHistClusterEvsTrackEGamma->Fill(cluster->E(), eMax, weight);
        if (classification == 3)
          fHistClusterEvsTrackEGammaSubCharged->Fill(cluster->E(), eMax, weight);

        labelsMatchedTracks.clear();
      }

      return kFALSE;
    }
  // special case for PHOS TM from tender
  } else if (fUseDistTrackToCluster && fUsePtDepTrackToCluster == 2){
      if( cluster->GetEmcCpvDistance() < fMinTMDistSigma ){
          if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//2
          return kFALSE;
      }
  }


  cutIndex++;//9, next cut

  if(GetClusterType() == 1 || GetClusterType() == 3) {
    if (cluster->E() > 0.5) fIsAcceptedForBasic = kTRUE;
  } else if (GetClusterType() == 2 ){
    if (cluster->E() > 0.3) fIsAcceptedForBasic = kTRUE;
  }

  // minimum cluster energy cut
  if (fUseMinEnergy){
    if(cluster->E() < fMinEnergy){
      if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//4
      return kFALSE;
    }
  }
  cutIndex++;//10, next cut


  // DONE with selecting photons
  if(fHistClusterIdentificationCuts)fHistClusterIdentificationCuts->Fill(cutIndex, cluster->E());//10

  // Histos after Cuts

  if(fHistClusterEtavsPhiAfterQA && cluster->E()>0. ) fHistClusterEtavsPhiAfterQA->Fill(phiCluster, etaCluster, weight);
  if(fHistClusterTimevsEAfterQA) fHistClusterTimevsEAfterQA->Fill(cluster->GetTOF(), cluster->E(), weight);
  if(fHistEnergyOfClusterAfterQA) fHistEnergyOfClusterAfterQA->Fill(cluster->E(), weight);
  if(fHistNCellsAfterQA) fHistNCellsAfterQA->Fill(cluster->GetNCells(), weight);
  if(fHistM02AfterQA) fHistM02AfterQA->Fill(cluster->GetM02(), weight);
  if(fHistM20AfterQA) fHistM20AfterQA->Fill(cluster->GetM20(), weight);
  if(fHistDispersionAfterQA) fHistDispersionAfterQA->Fill(cluster->GetDispersion(), weight);
  if(fHistNLMAfterQA) fHistNLMAfterQA->Fill(nLM, weight);
  if(fHistNLMVsNCellsAfterQA) fHistNLMVsNCellsAfterQA->Fill(nLM,cluster->GetNCells(), weight);
  if(fHistNLMVsEAfterQA) fHistNLMVsEAfterQA->Fill(nLM, cluster->E(), weight);
  if(fHistClusterEM02AfterQA) fHistClusterEM02AfterQA->Fill(cluster->E(), cluster->GetM02(), weight);

  if(fExtendedMatchAndQA > 1){
    if(fHistClusterIncludedCellsAfterQA){
      Int_t nCellCluster = cluster->GetNCells();
      for(Int_t iCell=0;iCell<nCellCluster;iCell++){
        Int_t cellID = cluster->GetCellAbsId(iCell);
        Double_t cellAmp = cells->GetCellAmplitude(cellID);
        Double_t cellTime = cells->GetCellTime(cellID);
        fHistClusterIncludedCellsAfterQA->Fill(cellID);
        if(cluster->E()>0.) fHistClusterEnergyFracCellsAfterQA->Fill(cellID,cellAmp/cluster->E());
        if(isMC==0){
          fHistClusterIncludedCellsTimingAfterQA->Fill(cluster->E(),cellTime*1E9);
          fHistClusterIncludedCellsTimingEnergyAfterQA->Fill(cellAmp,cellTime*1E9);
        }else{
          fHistClusterIncludedCellsTimingAfterQA->Fill(cluster->E(),cellTime*1E8);
          fHistClusterIncludedCellsTimingEnergyAfterQA->Fill(cellAmp,cellTime*1E8);
        }
      }
    }

    if(fHistClusterEnergyvsNCellsAfterQA) fHistClusterEnergyvsNCellsAfterQA->Fill(cluster->E(),cluster->GetNCells());
    if(cluster->IsEMCAL()){
      Int_t iSuperModule = -1;
      fGeomEMCAL = AliEMCALGeometry::GetInstance();
      if(!fGeomEMCAL){ AliFatal("EMCal geometry not initialized!");}
      if(fHistClusterEnergyvsMod && fGeomEMCAL->SuperModuleNumberFromEtaPhi(clusterVector.Eta(),clusterVector.Phi(),iSuperModule)){
        fHistClusterEnergyvsMod->Fill(cluster->E(),iSuperModule);
      }
    }else if(cluster->IsPHOS()){
      Int_t relId[4] = {0,0,0,0};
      fGeomPHOS = AliPHOSGeometry::GetInstance();
      if(!fGeomPHOS){ AliFatal("PHOS geometry not initialized!");}
      if(fHistClusterEnergyvsMod && fGeomPHOS->GlobalPos2RelId(clusterVector,relId)){
        fHistClusterEnergyvsMod->Fill(cluster->E(),relId[0]);
      }
    }
  }

  return kTRUE;
}



//________________________________________________________________________
void AliCaloPhotonCuts::FillHistogramsExtendedQA(AliVEvent *event, Int_t isMC)
{
  if(fExtendedMatchAndQA < 2) return;

  AliVCaloCells* cells = 0x0;

  Int_t nModules = 0;
  if( (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) && !fEMCALInitialized ) InitializeEMCAL(event);
  if( fClusterType == 2 && ( !fPHOSInitialized || (fPHOSCurrentRun != event->GetRunNumber()) ) ) InitializePHOS(event);

  Int_t nModulesStart = 0;
  if( fClusterType == 1 || fClusterType == 3 || fClusterType == 4){ //EMCAL & DCAL
    cells = event->GetEMCALCells();
    fGeomEMCAL = AliEMCALGeometry::GetInstance();
    if(!fGeomEMCAL) AliFatal("EMCal geometry not initialized!");
    if(!fEMCALBadChannelsMap && !fEMCALBadChannelsMap1D) AliFatal("EMCal bad channels map not initialized!");
    nModules = fGeomEMCAL->GetNumberOfSuperModules();
    if( fClusterType == 3) {nModules = 8; nModulesStart = 12;}
    if( fClusterType == 4) {nModules = 20;}
  } else if( fClusterType == 2 ){ //PHOS
    cells = event->GetPHOSCells();
    fGeomPHOS = AliPHOSGeometry::GetInstance();
    if(!fGeomPHOS) AliFatal("PHOS geometry not initialized!");
    nModules = fGeomPHOS->GetNModules();
  } else{
    AliError(Form("fExtendedMatchAndQA(%i):FillHistogramsExtendedMatchAndQA() not (yet) defined for cluster type (%i)",fExtendedMatchAndQA,fClusterType));
  }

  std::vector<Int_t> nCellsBigger100MeV(nModules);
  std::vector<Int_t> nCellsBigger1500MeV(nModules);
  std::vector<Double_t> EnergyOfMod(nModules);

  for(Int_t iModule=0;iModule<nModules;iModule++){nCellsBigger100MeV[iModule]=0;nCellsBigger1500MeV[iModule]=0;EnergyOfMod[iModule]=0;}


  for(Int_t iCell=0;iCell<cells->GetNumberOfCells();iCell++){
    Short_t cellNumber=0;
    Double_t cellAmplitude=0;
    Double_t cellTime=0;
    Double_t cellEFrac=0;
    Int_t cellMCLabel=0;
    Int_t nMod = -1;

    cells->GetCell(iCell,cellNumber,cellAmplitude,cellTime,cellMCLabel,cellEFrac);
    if( fClusterType == 3 && cellNumber < 12288){continue;}
    if( fClusterType == 2 && cellNumber < 0){continue;} //Scip CPV cells in PHOS case
    Int_t imod = -1;Int_t iTower = -1, iIphi = -1, iIeta = -1;
    Int_t icol = -1;Int_t irow = -1;
    Int_t relid[4];// for PHOS

    Bool_t doBadCell = kTRUE;
    if( fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
      nMod = fGeomEMCAL->GetSuperModuleNumber(cellNumber);
      fGeomEMCAL->GetCellIndex(cellNumber,imod,iTower,iIphi,iIeta);
      if (fEMCALBadChannelsMap->GetEntries() <= imod && !fEMCALBadChannelsMap1D) doBadCell=kFALSE;
      fGeomEMCAL->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi,iIeta,irow,icol);
    }else if( fClusterType == 2 ){
      fGeomPHOS->AbsToRelNumbering(cellNumber,relid);
      if(relid[1]!=0) AliFatal("PHOS CPV in PHOS cell array?");
      nMod = relid[0];//(Int_t) (1 + (cellNumber-1)/3584);
      if(nMod>=nModules || nMod<0 || !fPHOSBadChannelsMap[nMod]) doBadCell=kFALSE;
    }

    Int_t iBadCell = 0;
    if( (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) && doBadCell){
      if(fEMCALBadChannelsMap1D)
        iBadCell = (Int_t) fEMCALBadChannelsMap1D->GetBinContent(cellNumber);
      else
        iBadCell = (Int_t) ((TH2I*)fEMCALBadChannelsMap->At(imod))->GetBinContent(icol,irow);
    }else if( fClusterType == 2 && doBadCell){
      iBadCell = (Int_t) ((TH2I*)fPHOSBadChannelsMap[nMod])->GetBinContent(relid[2],relid[3]);
    }

    if(iBadCell > 0) continue;
	// nModulesStart == 0 for EMCAL and PHOS
    if(cellAmplitude > 0.1) nCellsBigger100MeV[nMod-nModulesStart]++;
    if(cellAmplitude > 1.5) nCellsBigger1500MeV[nMod-nModulesStart]++;
    if(cellAmplitude > 0.05) EnergyOfMod[nMod-nModulesStart]+=cellAmplitude;

    if(fExtendedMatchAndQA > 3){
      if(fHistCellEnergyvsCellID && (cellAmplitude > 0.05) && (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) )
          fHistCellEnergyvsCellID->Fill(cellAmplitude,cellNumber);
      if(fHistCellEnergyvsCellID && (cellAmplitude > 0.01) && fClusterType == 2 )
          fHistCellEnergyvsCellID->Fill(cellAmplitude,cellNumber);
      if(fHistCellTimevsCellID && (cellAmplitude > 0.01) && fClusterType == 2 )
          fHistCellTimevsCellID->Fill(cellTime,cellNumber);
      else if(fHistCellTimevsCellID && (cellAmplitude > 0.2)) fHistCellTimevsCellID->Fill(cellTime,cellNumber);
    }
  }

  for(Int_t iModule=0;iModule<nModules;iModule++){
    if(fHistNCellsBigger100MeVvsMod) fHistNCellsBigger100MeVvsMod->Fill(nCellsBigger100MeV[iModule],iModule+nModulesStart);
    if(fHistNCellsBigger1500MeVvsMod) fHistNCellsBigger1500MeVvsMod->Fill(nCellsBigger1500MeV[iModule],iModule+nModulesStart);
    if(fHistEnergyOfModvsMod) fHistEnergyOfModvsMod->Fill(EnergyOfMod[iModule],iModule+nModulesStart);
  }

  //fill distClusterTo_withinTiming/outsideTiming
  Int_t nclus = 0;
  TClonesArray * arrClustersExtQA = NULL;
  if(!fCorrTaskSetting.CompareTo("")){
    nclus = event->GetNumberOfCaloClusters();
  } else {
    arrClustersExtQA = dynamic_cast<TClonesArray*>(event->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!arrClustersExtQA)
      AliFatal(Form("%sClustersBranch was not found in AliCaloPhotonCuts::FillHistogramsExtendedQA! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclus = arrClustersExtQA->GetEntries();
  }
  std::unique_ptr<AliVCluster> cluster, clusterMatched;
  for(Int_t iClus=0; iClus<nclus ; iClus++){
    if(event->IsA()==AliESDEvent::Class()){
      if(arrClustersExtQA)
        cluster = std::unique_ptr<AliVCluster>(new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersExtQA->At(iClus)));
      else
        cluster = std::unique_ptr<AliVCluster>(new AliESDCaloCluster(*(AliESDCaloCluster*)event->GetCaloCluster(iClus)));
    } else if(event->IsA()==AliAODEvent::Class()){
      if(arrClustersExtQA)
        cluster = std::unique_ptr<AliVCluster>(new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersExtQA->At(iClus)));
      else
        cluster = std::unique_ptr<AliVCluster>(new AliAODCaloCluster(*(AliAODCaloCluster*)event->GetCaloCluster(iClus)));
    }

    if( (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) && !cluster->IsEMCAL()){continue;}
    if( fClusterType == 2 && cluster->GetType() !=AliVCluster::kPHOSNeutral){continue;}

    Float_t clusPos[3]={0,0,0};
    cluster->GetPosition(clusPos);
    TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
    Double_t etaCluster = clusterVector.Eta();
    Double_t phiCluster = clusterVector.Phi();
    if (phiCluster < 0) phiCluster += 2*TMath::Pi();
    Int_t nLM = GetNumberOfLocalMaxima(cluster.get(), event);

    //acceptance cuts
    if (fUseEtaCut && (etaCluster < fMinEtaCut || etaCluster > fMaxEtaCut)){continue;}
    if (fUseEtaCut && fClusterType == 3 && etaCluster < fMaxEtaInnerEdge && etaCluster > fMinEtaInnerEdge ) {continue;}
    if (fClusterType == 4){
      if (fUsePhiCut && (phiCluster < fMinPhiCut || phiCluster > fMaxPhiCut) && (phiCluster < fMinPhiCutDMC || phiCluster > fMaxPhiCutDMC)){continue;}
    } else {
      if (fUsePhiCut && (phiCluster < fMinPhiCut || phiCluster > fMaxPhiCut)){continue;}
    }
    if (fUseDistanceToBadChannel>0 && CheckDistanceToBadChannel(cluster.get(),event)){continue;}
    //cluster quality cuts
    if (fVectorMatchedClusterIDs.size()>0 && CheckClusterForTrackMatch(cluster.get())){continue;}
    if (fUseMinEnergy && (cluster->E() < fMinEnergy)){continue;}
    if (fUseNCells && (cluster->GetNCells() < fMinNCells)){continue;}
    if (fUseNLM && (nLM < fMinNLM || nLM > fMaxNLM)){continue;}
    if(!fUseNCells && cluster->GetNCells()<2 && cluster->E()<4){
      // no cut to be applied in this case on M20
      // as cluster needs at least 2 cells for M20 calculation
    } else {
      if (fUseM02 == 1 && (cluster->GetM02() < fMinM02 || cluster->GetM02() > fMaxM02)){continue;}
      if (fUseM02 == 2 && (cluster->GetM02() < CalculateMinM02(fMinM02CutNr, cluster->E()) || cluster->GetM02() > CalculateMaxM02(fMaxM02CutNr, cluster->E()))){continue;}
      if (fUseM20 && (cluster->GetM20() < fMinM20 || cluster->GetM20() > fMaxM20)){continue;}
      if (fUseDispersion && (cluster->GetDispersion() > fMaxDispersion)){continue;}
    }
    //cluster within timing cut
    if( fUseTimingEfficiencyMCSimCluster==2 ){
      if ( cluster->E() < 5 ) {
        if (!(isMC>0) && (cluster->GetTOF() < fMinTimeDiff || cluster->GetTOF() > fMaxTimeDiff)){continue;}
      } else {
        if (!(isMC>0) && (cluster->GetTOF() < fMinTimeDiffHighPt || cluster->GetTOF() > fMaxTimeDiffHighPt)){continue;}
      }
    } else {
      if (!(isMC>0) && (cluster->GetTOF() < fMinTimeDiff || cluster->GetTOF() > fMaxTimeDiff)){continue;}
    }
    Int_t largestCellicol = -1, largestCellirow = -1;
    Int_t largestCellID = FindLargestCellInCluster(cluster.get(),event);
    if(largestCellID==-1) AliFatal("FillHistogramsExtendedQA: FindLargestCellInCluster found cluster with NCells<1?");
    Int_t largestCelliMod = GetModuleNumberAndCellPosition(largestCellID, largestCellicol, largestCellirow);
    if(largestCelliMod < 0) AliFatal("FillHistogramsExtendedQA: GetModuleNumberAndCellPosition found SM with ID<0?");

    for(Int_t iClus2=iClus+1; iClus2<nclus; iClus2++){
      if(event->IsA()==AliESDEvent::Class()){
        if(arrClustersExtQA)
          clusterMatched = std::unique_ptr<AliVCluster>(new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersExtQA->At(iClus2)));
        else
          clusterMatched = std::unique_ptr<AliVCluster>(new AliESDCaloCluster(*(AliESDCaloCluster*)event->GetCaloCluster(iClus2)));
      } else if(event->IsA()==AliAODEvent::Class()){
        if(arrClustersExtQA)
          clusterMatched = std::unique_ptr<AliVCluster>(new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersExtQA->At(iClus2)));
        else
          clusterMatched = std::unique_ptr<AliVCluster>(new AliAODCaloCluster(*(AliAODCaloCluster*)event->GetCaloCluster(iClus2)));
      }

      if( (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) && !clusterMatched->IsEMCAL()){continue;}
      if( fClusterType == 2 && clusterMatched->GetType() !=AliVCluster::kPHOSNeutral){continue;}

      Float_t clusPos2[3]={0,0,0};
      clusterMatched->GetPosition(clusPos2);
      TVector3 clusterMatchedVector(clusPos2[0],clusPos2[1],clusPos2[2]);
      Double_t etaclusterMatched = clusterMatchedVector.Eta();
      Double_t phiclusterMatched = clusterMatchedVector.Phi();
      if (phiclusterMatched < 0) phiclusterMatched += 2*TMath::Pi();
      Int_t nLMMatched = GetNumberOfLocalMaxima(clusterMatched.get(), event);

      //acceptance cuts
      if (fUseEtaCut && (etaclusterMatched < fMinEtaCut || etaclusterMatched > fMaxEtaCut)){continue;}
      if (fUseEtaCut && fClusterType == 3 && etaclusterMatched < fMaxEtaInnerEdge && etaclusterMatched > fMinEtaInnerEdge ) {continue;}
      if (fClusterType == 4){
        if (fUsePhiCut && (phiclusterMatched < fMinPhiCut || phiclusterMatched > fMaxPhiCut) && (phiclusterMatched < fMinPhiCutDMC || phiclusterMatched > fMaxPhiCutDMC)){continue;}
      } else {
        if (fUsePhiCut && (phiclusterMatched < fMinPhiCut || phiclusterMatched > fMaxPhiCut)){continue;}
      }
      if (fUseDistanceToBadChannel>0 && CheckDistanceToBadChannel(clusterMatched.get(),event)){continue;}
      //cluster quality cuts
      if (fVectorMatchedClusterIDs.size()>0 && CheckClusterForTrackMatch(clusterMatched.get())){continue;}
      if (fUseMinEnergy && (clusterMatched->E() < fMinEnergy)){continue;}
      if (fUseNCells && (clusterMatched->GetNCells() < fMinNCells)){continue;}
      if (fUseNLM && (nLMMatched < fMinNLM || nLMMatched > fMaxNLM)){continue;}
      if(!fUseNCells && cluster->GetNCells()<2 && cluster->E()<4){
        // no cut to be applied in this case on M20
        // as cluster needs at least 2 cells for M20 calculation
      } else {
        if (fUseM02 == 1 && (clusterMatched->GetM02() < fMinM02 || clusterMatched->GetM02() > fMaxM02)){continue;}
        if (fUseM02 == 2 && (clusterMatched->GetM02() < CalculateMinM02(fMinM02CutNr, clusterMatched->E()) || cluster->GetM02() > CalculateMaxM02(fMaxM02CutNr, clusterMatched->E()))){continue;}
        if (fUseM20 && (clusterMatched->GetM20() < fMinM20 || clusterMatched->GetM20() > fMaxM20)){continue;}
        if (fUseDispersion && (clusterMatched->GetDispersion() > fMaxDispersion)){continue;}
      }
      // Get rowdiff and coldiff

      Int_t matched_largestCellicol = -1, matched_largestCellirow = -1;
      Int_t matched_largestCellID = FindLargestCellInCluster(clusterMatched.get(),event);
      if(matched_largestCellID==-1) AliFatal("FillHistogramsExtendedQA: FindLargestCellInCluster found cluster with NCells<1?");
      Int_t matched_largestCelliMod = GetModuleNumberAndCellPosition(matched_largestCellID, matched_largestCellicol, matched_largestCellirow);
      if(matched_largestCelliMod < 0) AliFatal("FillHistogramsExtendedQA: GetModuleNumberAndCellPosition found SM with ID<0?");

//      cout << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
//      cout << "Cluster: " << largestCelliMod << ", " << largestCellirow << ", " << largestCellicol << " , time: " << cluster->GetTOF() << endl;
//      cout << "Matched: " << matched_largestCelliMod << ", " << matched_largestCellirow << ", " << matched_largestCellicol << " , time: " << clusterMatched->GetTOF() << endl;

      Int_t rowdiff = -100;
      Int_t coldiff = -100;
      Bool_t calculatedDiff = kFALSE;

      Int_t ClusID = largestCelliMod/2;
      Int_t matchClusID = matched_largestCelliMod/2;

      if( matched_largestCelliMod == largestCelliMod){
        rowdiff = largestCellirow - matched_largestCellirow;
        coldiff = largestCellicol - matched_largestCellicol;
        calculatedDiff = kTRUE;
      }else if( TMath::Abs(matched_largestCelliMod - largestCelliMod) == 1 && (ClusID == matchClusID) ){
        if(matched_largestCelliMod%2){
          matched_largestCelliMod -= 1;
          matched_largestCellicol += AliEMCALGeoParams::fgkEMCALCols;
        }else{
          matched_largestCelliMod += 1;
          matched_largestCellicol -= AliEMCALGeoParams::fgkEMCALCols;
        }

        if( matched_largestCelliMod == largestCelliMod ){
          rowdiff = largestCellirow - matched_largestCellirow;
          coldiff = largestCellicol - matched_largestCellicol;
          calculatedDiff = kTRUE;
        }
      }
//      cout << "\t\t ROWDIFF: " << rowdiff << endl;
//      cout << "\t\t COLDIFF: " << coldiff << endl;
//      cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n" << endl;
      //cluster outside timing cut
      if( calculatedDiff ){
        Float_t dist1D = TMath::Sqrt(TMath::Power(etaCluster-etaclusterMatched,2)+TMath::Power(phiCluster-phiclusterMatched,2));
        if( !(isMC>0) ){
          if( (fUseTimingEfficiencyMCSimCluster==2) ){
            if ( cluster->E() < 5) {
              if( (clusterMatched->GetTOF() > fMinTimeDiff && clusterMatched->GetTOF() < fMaxTimeDiff) ){
                fHistClusterDistanceInTimeCut->Fill(rowdiff,coldiff);
                fHistClusterDistance1DInTimeCut->Fill(dist1D);
              }
              else fHistClusterDistanceOutTimeCut->Fill(rowdiff,coldiff);
            }
            else {
              if( (clusterMatched->GetTOF() > fMinTimeDiffHighPt && clusterMatched->GetTOF() < fMaxTimeDiffHighPt) ){
                fHistClusterDistanceInTimeCut->Fill(rowdiff,coldiff);
                fHistClusterDistance1DInTimeCut->Fill(dist1D);
              }
              else fHistClusterDistanceOutTimeCut->Fill(rowdiff,coldiff);
            }
          } else {
            if( (clusterMatched->GetTOF() > fMinTimeDiff && clusterMatched->GetTOF() < fMaxTimeDiff) ){
              fHistClusterDistanceInTimeCut->Fill(rowdiff,coldiff);
              fHistClusterDistance1DInTimeCut->Fill(dist1D);
            }
            else fHistClusterDistanceOutTimeCut->Fill(rowdiff,coldiff);
          }
        }else{
          fHistClusterDistanceInTimeCut->Fill(rowdiff,coldiff);
          fHistClusterDistance1DInTimeCut->Fill(dist1D);
        }
      }
    }
  }
  return;
}
//________________________________________________________________________
Double_t AliCaloPhotonCuts::GetTotalEnergyDeposit(AliVEvent *event)
{

  AliVCaloCells* cells = 0x0;
  Int_t nModules = 0;
  Double_t totalCellAmplitude = 0;

  if( (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) && !fEMCALInitialized ) InitializeEMCAL(event);
  if( fClusterType == 2 && ( !fPHOSInitialized || (fPHOSCurrentRun != event->GetRunNumber()) ) ) InitializePHOS(event);

  if( fClusterType == 1 || fClusterType == 3 || fClusterType == 4){ //EMCAL & DCAL
    cells = event->GetEMCALCells();
    fGeomEMCAL = AliEMCALGeometry::GetInstance();
    if(!fGeomEMCAL) AliFatal("EMCal geometry not initialized!");
    if(!fEMCALBadChannelsMap && !fEMCALBadChannelsMap1D) AliFatal("EMCal bad channels map not initialized!");
    nModules = fGeomEMCAL->GetNumberOfSuperModules();
  } else if( fClusterType == 2 ){ //PHOS
    cells = event->GetPHOSCells();
    fGeomPHOS = AliPHOSGeometry::GetInstance();
    if(!fGeomPHOS) AliFatal("PHOS geometry not initialized!");
    nModules = fGeomPHOS->GetNModules();
  } else{
    AliError(Form("GetTotalEnergyDeposit not (yet) defined for cluster type (%i)",fClusterType));
  }

  for(Int_t iCell=0;iCell<cells->GetNumberOfCells();iCell++){
    Short_t cellNumber=0;
    Double_t cellAmplitude=0;
    Double_t cellTime=0;
    Double_t cellEFrac=0;
    Int_t cellMCLabel=0;
    Int_t nMod = -1;

    cells->GetCell(iCell,cellNumber,cellAmplitude,cellTime,cellMCLabel,cellEFrac);
    if( fClusterType == 3 && cellNumber < 12288){continue;}
    if( fClusterType == 2 && cellNumber < 0){continue;} //Scip CPV cells in PHOS case
    Int_t imod = -1;Int_t iTower = -1, iIphi = -1, iIeta = -1;
    Int_t icol = -1;Int_t irow = -1;
    Int_t relid[4];// for PHOS

    Bool_t doBadCell = kTRUE;
    if( fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
      nMod = fGeomEMCAL->GetSuperModuleNumber(cellNumber);
      fGeomEMCAL->GetCellIndex(cellNumber,imod,iTower,iIphi,iIeta);
      if (fEMCALBadChannelsMap->GetEntries() <= imod && !fEMCALBadChannelsMap1D) doBadCell=kFALSE;
      fGeomEMCAL->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi,iIeta,irow,icol);
    }else if( fClusterType == 2 ){
      fGeomPHOS->AbsToRelNumbering(cellNumber,relid);
      if(relid[1]!=0) AliFatal("PHOS CPV in PHOS cell array?");
      nMod = relid[0];//(Int_t) (1 + (cellNumber-1)/3584);
      if(nMod>=nModules || nMod<0 || !fPHOSBadChannelsMap[nMod]) doBadCell=kFALSE;
    }

    Int_t iBadCell = 0;
    if( (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) && doBadCell){
      if(fEMCALBadChannelsMap1D)
        iBadCell = (Int_t) fEMCALBadChannelsMap1D->GetBinContent(cellNumber);
      else
        iBadCell = (Int_t) ((TH2I*)fEMCALBadChannelsMap->At(imod))->GetBinContent(icol,irow);
    }else if( fClusterType == 2 && doBadCell){
      iBadCell = (Int_t) ((TH2I*)fPHOSBadChannelsMap[nMod])->GetBinContent(relid[2],relid[3]);
    }

    if(iBadCell > 0) continue;
    totalCellAmplitude += cellAmplitude;
  }

  return totalCellAmplitude;
}

//________________________________________________________________________
//************** Find number of local maxima in cluster ******************
//* derived from G. Conesa Balbastre's AliCalorimeterUtils *******************
//************************************************************************
Int_t AliCaloPhotonCuts::GetNumberOfLocalMaxima(AliVCluster* cluster, AliVEvent * event){

  const Int_t   nc = cluster->GetNCells();

  Int_t   absCellIdList[nc];
  Float_t maxEList[nc];

  Int_t nMax = GetNumberOfLocalMaxima(cluster, event, absCellIdList, maxEList);

  return nMax;
}

//________________________________________________________________________
Int_t AliCaloPhotonCuts::FindSecondLargestCellInCluster(AliVCluster* cluster, AliVEvent* event){

  const Int_t nCells      = cluster->GetNCells();
  AliVCaloCells* cells    = NULL;

  if (fClusterType == 1 || fClusterType == 3 || fClusterType == 4)
    cells                 = event->GetEMCALCells();
  else if (fClusterType ==2 )
    cells                 = event->GetPHOSCells();

//   cout << "NCells: "<< nCells<< " cluster energy: " << cluster->E() << endl;
  Float_t eMax            = 0.;
  Int_t idMax             = -1;
  Int_t idMax2            = -1;
  Int_t iCellMax          = -1;

  if (nCells < 2) return idMax;
  for (Int_t iCell = 1;iCell < nCells;iCell++){
    if (cells->GetCellAmplitude(cluster->GetCellsAbsId()[iCell])> eMax){
      eMax                = cells->GetCellAmplitude(cluster->GetCellsAbsId()[iCell]);
      idMax               = cluster->GetCellsAbsId()[iCell];
      iCellMax            = iCell;
    }
  }

  eMax                    = 0.;
  for (Int_t iCell = 1;iCell < nCells;iCell++){
    if (iCell == iCellMax) continue;
    if (cells->GetCellAmplitude(cluster->GetCellsAbsId()[iCell])> eMax){
      eMax                = cells->GetCellAmplitude(cluster->GetCellsAbsId()[iCell]);
      idMax2              = cluster->GetCellsAbsId()[iCell];
    }
  }

  return idMax2;
}

//________________________________________________________________________
Int_t AliCaloPhotonCuts::FindLargestCellInCluster(AliVCluster* cluster, AliVEvent* event){

  const Int_t nCells      = cluster->GetNCells();
  AliVCaloCells* cells    = NULL;

  if (fClusterType == 1 || fClusterType == 3 || fClusterType == 4)
    cells                 = event->GetEMCALCells();
  else if (fClusterType ==2 )
    cells                 = event->GetPHOSCells();

//   cout << "NCells: "<< nCells<< " cluster energy: " << cluster->E() << endl;
  Float_t eMax            = 0.;
  Int_t idMax             = -1;

  if (nCells < 1) return idMax;
  for (Int_t iCell = 0;iCell < nCells;iCell++){
    Int_t cellAbsID       = cluster->GetCellsAbsId()[iCell];
    if (cells->GetCellAmplitude(cellAbsID)> eMax){
      eMax                = cells->GetCellAmplitude(cellAbsID);
      idMax               = cellAbsID;
    }
  }
  return idMax;

}


//________________________________________________________________________
//************** Find number of local maxima in cluster ******************
//* derived from G. Conesa Balbastre's AliCalorimeterUtils ***************
//************************************************************************
Int_t AliCaloPhotonCuts::GetNumberOfLocalMaxima(AliVCluster* cluster, AliVEvent * event, Int_t *absCellIdList, Float_t* maxEList){

  Int_t absCellId1        = -1;
  Int_t absCellId2        = -1;
  const Int_t nCells      = cluster->GetNCells();
  AliVCaloCells* cells    = NULL;

  if (fClusterType == 1 || fClusterType == 3 || fClusterType == 4)
    cells                 = event->GetEMCALCells();
  else if (fClusterType ==2 )
    cells                 = event->GetPHOSCells();

//   cout << "NCells: "<< nCells<< " cluster energy: " << cluster->E() << endl;
  Float_t eMax            = 0.;
  Int_t idMax             = -1;

  for (Int_t iCell = 0;iCell < nCells;iCell++){
    absCellIdList[iCell]  = cluster->GetCellsAbsId()[iCell];
//     Int_t imod = -1, icol = -1, irow = -1;
//     imod = GetModuleNumberAndCellPosition(absCellIdList[iCell], icol, irow);
//     cout << absCellIdList[iCell] <<"\t" << cells->GetCellAmplitude(absCellIdList[iCell]) << "\t"<< imod << "\t" << icol << "\t" << irow << endl;
    if (cells->GetCellAmplitude(absCellIdList[iCell])> eMax){
      eMax                = cells->GetCellAmplitude(absCellIdList[iCell]);
      idMax               = absCellIdList[iCell];
    }
  }

  // find the largest separated cells in cluster
  for (Int_t iCell = 0;iCell < nCells;iCell++){
    // check whether valid cell number is selected
    if (absCellIdList[iCell] >= 0){
      // store current energy and cell id
      absCellId1          = cluster->GetCellsAbsId()[iCell];
      Float_t en1         = cells->GetCellAmplitude(absCellId1);

      // loop over other cells in cluster
      for (Int_t iCellN = 0;iCellN < nCells;iCellN++){
        // jump out if array has changed in the meantime
        if (absCellIdList[iCell] == -1) continue;
        // get cell id & check whether its valid
        absCellId2        = cluster->GetCellsAbsId()[iCellN];

        // don't compare to yourself
        if (absCellId2 == -1) continue;
        if (absCellId1 == absCellId2) continue;

        // get cell energy of second cell
        Float_t en2       = cells->GetCellAmplitude(absCellId2);

        // check if cells are Neighbours
        if (AreNeighbours(absCellId1, absCellId2)){
          // determine which cell has larger energy, mask the other
//           cout << "found neighbour: " << absCellId1 << "\t" << absCellId2 << endl;
//           cout << "energies: " << en1 << "\t" << en2 << endl;
          if (en1 > en2 ){
            absCellIdList[iCellN]       = -1;
            if (en1 < en2 + fLocMaxCutEDiff)
                absCellIdList[iCell]    = -1;
          } else {
            absCellIdList[iCell]        = -1;
            if (en1 > en2 - fLocMaxCutEDiff)
                absCellIdList[iCellN]   = -1;
          }
        }
      }
    }
  }

  // shrink list of cells to only maxima
  Int_t nMaximaNew        = 0;
  for (Int_t iCell = 0;iCell < nCells;iCell++){
//     cout << iCell << "\t" << absCellIdList[iCell] << endl;
    if (absCellIdList[iCell] > -1){
      Float_t en          = cells->GetCellAmplitude(absCellIdList[iCell]);
      // check whether cell energy is larger than required seed
      if (en < fSeedEnergy) continue;
      absCellIdList[nMaximaNew]   = absCellIdList[iCell];
      maxEList[nMaximaNew]        = en;
      nMaximaNew++;
    }
  }

  // check whether a local maximum was found
  // if no maximum was found use highest cell as maximum
  if (nMaximaNew == 0){
    nMaximaNew            = 1;
    maxEList[0]           = eMax;
    absCellIdList[0]      = idMax;
  }

  return nMaximaNew;
}

//________________________________________________________________________
//************** Function to determine neighbours of cells ***************
//* derived from G. Conesa Balbastre's AliCalorimeterUtils ***************
//************************************************************************
Bool_t AliCaloPhotonCuts::AreNeighbours(Int_t absCellId1, Int_t absCellId2){
  Bool_t areNeighbours = kFALSE ;
  Int_t irow1 = -1, icol1 = -1;
  Int_t irow2 = -1, icol2 = -1;

  Int_t rowdiff =  0, coldiff =  0;

  Int_t nSupMod1          = GetModuleNumberAndCellPosition(absCellId1, icol1, irow1);
  Int_t nSupMod2          = GetModuleNumberAndCellPosition(absCellId2, icol2, irow2);
  // check if super modules are correct
  if (nSupMod1== -1 || nSupMod2 == -1) return areNeighbours;

  if((fClusterType==1 || fClusterType==4) && nSupMod1!=nSupMod2) {
    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2-1
    // C Side impair SM, nSupMod%2=1;A side pair SM nSupMod%2=0
    if(nSupMod1%2) icol1+=AliEMCALGeoParams::fgkEMCALCols;
    else           icol2+=AliEMCALGeoParams::fgkEMCALCols;
  }

  rowdiff = TMath::Abs( irow1 - irow2 ) ;
  coldiff = TMath::Abs( icol1 - icol2 ) ;

  if ((coldiff + rowdiff == 1 )){
    areNeighbours         = kTRUE ;
  }

  return areNeighbours;
}


//________________________________________________________________________
//************** Function to obtain module number, row and column ********
//* derived from G. Conesa Balbastre's AliCalorimeterUtils ***************
//************************************************************************
Int_t AliCaloPhotonCuts::GetModuleNumberAndCellPosition(Int_t absCellId, Int_t & icol, Int_t & irow){
  if( fClusterType == 1 || fClusterType == 3 || fClusterType == 4){ //EMCAL & DCAL
    fGeomEMCAL = AliEMCALGeometry::GetInstance();
    if(!fGeomEMCAL) AliFatal("EMCal geometry not initialized!");
  } else if( fClusterType == 2 ){ //PHOS
    fGeomPHOS = AliPHOSGeometry::GetInstance();
    if(!fGeomPHOS) AliFatal("PHOS geometry not initialized!");
  }

  Int_t imod = -1;Int_t iTower = -1, iIphi = -1, iIeta = -1;
  if( fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
    fGeomEMCAL->GetCellIndex(absCellId,imod,iTower,iIphi,iIeta);
    fGeomEMCAL->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi,iIeta,irow,icol);
  } else if ( fClusterType == 2 ){
    Int_t relId[4];
    fGeomPHOS->AbsToRelNumbering(absCellId,relId);
    irow                  = relId[2];
    icol                  = relId[3];
    imod                  = relId[0]-1;
  }
  return imod;
}

//___________________________________________________________________________
// Split energy of cluster between the 2 local maxima, sum energy on 3x3, and if the 2
// maxima are too close and have common cells, split the energy between the 2.
//* derived from G. Conesa Balbastre's AliCalorimeterUtils *******************
//___________________________________________________________________________
void AliCaloPhotonCuts::SplitEnergy(Int_t absCellId1, Int_t absCellId2,
                                    AliVCluster* cluster,
                                    AliVEvent* event,
                                    Int_t isMC,
                                    AliAODCaloCluster* cluster1,
                                    AliAODCaloCluster* cluster2){

  const Int_t ncells      = cluster->GetNCells();
  Int_t absCellIdList[ncells];

  AliVCaloCells* cells    = NULL;
  if (fClusterType == 1 || fClusterType == 3 || fClusterType == 4)
    cells                 = event->GetEMCALCells();
  else if (fClusterType ==2 )
    cells                 = event->GetPHOSCells();

  Float_t e1              = 0;
  Float_t e2              = 0;
  Float_t eCluster        = 0;

  for(Int_t iCell    = 0;iCell < ncells;iCell++ ) {
    absCellIdList[iCell]  = cluster->GetCellsAbsId()[iCell];
    Float_t ec            = cells->GetCellAmplitude(absCellIdList[iCell]);
    eCluster+=ec;
  }

  UShort_t absCellIdList1 [12];
  Double_t fracList1      [12];
  UShort_t absCellIdList2 [12];
  Double_t fracList2      [12];

  // Init counters and variables
  Int_t ncells1         = 1 ;
  absCellIdList1[0]     = absCellId1 ;
  fracList1 [0]         = 1. ;

  Float_t ecell1        = cells->GetCellAmplitude(absCellId1);
  e1                    = ecell1;

  Int_t ncells2         = 1 ;
  absCellIdList2[0]     = absCellId2 ;
  fracList2 [0]         = 1. ;

  Float_t ecell2        = cells->GetCellAmplitude(absCellId2);
  e2                    = ecell2;

//   cout << "Cluster: " << eCluster << "\t cell1: " << absCellId1 << "\t" << e1 << "\t cell2: " << absCellId2 << "\t" << e2 << endl;
  // Very rough way to share the cluster energy
  Float_t eRemain           = (eCluster-ecell1-ecell2)/2;
  Float_t shareFraction1    = (ecell1+eRemain)/eCluster;
  Float_t shareFraction2    = (ecell2+eRemain)/eCluster;

//   cout << eRemain << "\t" << shareFraction1<< "\t" << shareFraction2 << endl;

  for(Int_t iCell = 0;iCell < ncells;iCell++){

    Int_t absId         = absCellIdList[iCell];
    if ( absId==absCellId1 || absId==absCellId2 || absId < 0 ) continue;

    Float_t ecell = cells->GetCellAmplitude(absId);
    if(AreNeighbours(absCellId1,absId )){
      absCellIdList1[ncells1] = absId;
      if(AreNeighbours(absCellId2,absId )){
        fracList1[ncells1] = shareFraction1;
        e1 += ecell*shareFraction1;
      } else {
        fracList1[ncells1] = 1.;
        e1 += ecell;
      }
      ncells1++;
    } // neigbour to cell1

    if(AreNeighbours(absCellId2,absId )) {
      absCellIdList2[ncells2]= absId;

      if(AreNeighbours(absCellId1,absId )){
        fracList2[ncells2] = shareFraction2;
        e2 += ecell*shareFraction2;
      } else {
        fracList2[ncells2] = 1.;
        e2 += ecell;
      }
      ncells2++;
    } // neigbour to cell2
  }
//   cout << "Cluster: " << eCluster << "\t cell1: " << absCellId1 << "\t" << e1 << "\t cell2: " << absCellId2 << "\t" << e2 << endl;

  cluster1->SetE(e1);
  cluster2->SetE(e2);

  cluster1->SetNCells(ncells1);
  cluster2->SetNCells(ncells2);

  cluster1->SetCellsAbsId(absCellIdList1);
  cluster2->SetCellsAbsId(absCellIdList2);

  cluster1->SetCellsAmplitudeFraction(fracList1);
  cluster2->SetCellsAmplitudeFraction(fracList2);

  // Correct linearity
  ApplyNonLinearity(cluster1, isMC, event) ;
  ApplyNonLinearity(cluster2, isMC, event) ;
  if(isMC == 0){
    ApplySMWiseEnergyCorrection(cluster1, isMC, event);
    ApplySMWiseEnergyCorrection(cluster2, isMC, event);
  }

  // Initialize EMCAL rec utils if not initialized
  if(!fEMCALInitialized && (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) ) InitializeEMCAL(event);

  if(fEMCALInitialized && (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) ){
    fEMCALRecUtils->RecalculateClusterPosition(fGeomEMCAL, cells, cluster1);
    fEMCALRecUtils->RecalculateClusterPosition(fGeomEMCAL, cells, cluster2);
  }
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::CheckDistanceToBadChannel(AliVCluster* cluster, AliVEvent* event)
{
  if(fUseDistanceToBadChannel != 1 && fUseDistanceToBadChannel != 2) return kFALSE;

  if( (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) && !fEMCALInitialized ) InitializeEMCAL(event);
  if( fClusterType == 2 && ( !fPHOSInitialized || (fPHOSCurrentRun != event->GetRunNumber()) ) ) InitializePHOS(event);
  if( fClusterType == 2 ) fGeomPHOS = AliPHOSGeometry::GetInstance();

  Int_t largestCellID = FindLargestCellInCluster(cluster,event);
  if(largestCellID==-1) AliFatal("CheckDistanceToBadChannel: FindLargestCellInCluster found cluster with NCells<1?");

  Int_t largestCellicol = -1, largestCellirow = -1;
  Int_t rowdiff =  0, coldiff =  0;

  Int_t largestCelliMod = GetModuleNumberAndCellPosition(largestCellID, largestCellicol, largestCellirow);
  if(largestCelliMod < 0) AliFatal("CheckDistanceToBadChannel: GetModuleNumberAndCellPosition found SM with ID<0?");

  Int_t nMinRows = 0, nMaxRows = 0;
  Int_t nMinCols = 0, nMaxCols = 0;

  Bool_t checkNextSM = kFALSE;
  Int_t distanceForLoop = fMinDistanceToBadChannel+1;

  Bool_t isDCal = kFALSE;
  if(fClusterType == 4){
    Float_t clusPos[3]={0,0,0};
    cluster->GetPosition(clusPos);
    TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
    Double_t phiCluster = clusterVector.Phi();
    if (phiCluster < 0) phiCluster += 2*TMath::Pi();
    if (phiCluster > fMinPhiCut && phiCluster < fMaxPhiCut)
      isDCal = kFALSE;
    else if (phiCluster > fMinPhiCutDMC && phiCluster < fMaxPhiCutDMC)
      isDCal = kTRUE;
  }

  if( fClusterType == 1 || (fClusterType == 4 && !isDCal)){
    nMinRows = largestCellirow - distanceForLoop;
    nMaxRows = largestCellirow + distanceForLoop;
    if(nMinRows < 0) nMinRows = 0;
    if(nMaxRows > AliEMCALGeoParams::fgkEMCALRows) nMaxRows = AliEMCALGeoParams::fgkEMCALRows;

    nMinCols = largestCellicol - distanceForLoop;
    nMaxCols = largestCellicol + distanceForLoop;

    if(largestCelliMod%2){
      if(nMinCols < 0){
        nMinCols = 0;
        checkNextSM = kTRUE;
      }
      if(nMaxCols > AliEMCALGeoParams::fgkEMCALCols) nMaxCols = AliEMCALGeoParams::fgkEMCALCols;
    }else{
      if(nMinCols < 0) nMinCols = 0;
      if(nMaxCols > AliEMCALGeoParams::fgkEMCALCols){
        nMaxCols = AliEMCALGeoParams::fgkEMCALCols;
        checkNextSM = kTRUE;
      }
    }
  }else if( fClusterType == 3 || (fClusterType == 4 && isDCal)){
    nMinRows = largestCellirow - distanceForLoop;
    nMaxRows = largestCellirow + distanceForLoop;
    if(nMinRows < 0) nMinRows = 0;
    if(nMaxRows > AliEMCALGeoParams::fgkEMCALRows) nMaxRows = AliEMCALGeoParams::fgkEMCALRows; //AliEMCALGeoParams::fgkDCALRows; <- doesnt exist yet (DCAl = EMCAL here)

    nMinCols = largestCellicol - distanceForLoop;
    nMaxCols = largestCellicol + distanceForLoop;
    if(nMinCols < 0) nMinCols = 0;
    if(nMaxCols > fgkDCALCols) nMaxCols = fgkDCALCols; // AliEMCALGeoParams::fgkDCALCols; <- doesnt exist yet

  }else if( fClusterType == 2 ){
    nMinRows = largestCellirow - distanceForLoop;
    nMaxRows = largestCellirow + distanceForLoop;
    if (nMinRows < 0) nMinRows = 0;
    if (nMaxRows > fGeomPHOS->GetNPhi()) nMaxRows = fGeomPHOS->GetNPhi();

    nMinCols = largestCellicol - distanceForLoop;
    nMaxCols = largestCellicol + distanceForLoop;
    if(nMinCols < 0) nMinCols = 0;
    if(nMaxCols > fGeomPHOS->GetNZ()) nMaxCols = fGeomPHOS->GetNZ();
  }

//  cout << "Cluster: " << fClusterType << ",checkNextSM: " << checkNextSM << endl;
//  cout << "largestCell: " << largestCellID << ",mod: " << largestCelliMod << ",col: " << largestCellicol << ",row: " << largestCellirow << endl;
//  cout << "distanceForLoop: " << distanceForLoop << ",nMinRows: " << nMinRows << ",nMaxRows: " << nMaxRows << ",nMinCols: " << nMinCols << ",nMaxCols: " << nMaxCols << endl;

  //check bad cells within respective SM
  for (Int_t irow = nMinRows;irow < nMaxRows;irow++)
  {
    for (Int_t icol = nMinCols;icol < nMaxCols;icol++)
    {
      if(irow == largestCellirow && icol == largestCellicol) continue;

      Int_t iBadCell = 0;
      if( (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) && (largestCelliMod<fEMCALBadChannelsMap->GetEntries() || fEMCALBadChannelsMap1D) ){
        if(fEMCALBadChannelsMap1D)
          iBadCell = (Int_t) fEMCALBadChannelsMap1D->GetBinContent(fGeomEMCAL->GetAbsCellIdFromCellIndexes(largestCelliMod, icol, irow));
        else
          iBadCell = (Int_t) ((TH2I*)fEMCALBadChannelsMap->At(largestCelliMod))->GetBinContent(icol,irow);
      }else if( fClusterType == 2 && fPHOSBadChannelsMap[largestCelliMod+1]){
        iBadCell = (Int_t) ((TH2I*)fPHOSBadChannelsMap[largestCelliMod+1])->GetBinContent(icol,irow);
      }
      //cout << "largestCelliMod: " << largestCelliMod << ",iBadCell: " << iBadCell << ",icol: " << icol << ",irow: " << irow << endl;
      if(iBadCell==0) continue;

      rowdiff = TMath::Abs( largestCellirow - irow ) ;
      coldiff = TMath::Abs( largestCellicol - icol ) ;
      //cout << "rowdiff: " << rowdiff << ",coldiff: " << coldiff << endl;
      if(fUseDistanceToBadChannel==1){
        if ((coldiff + rowdiff <= fMinDistanceToBadChannel )) return kTRUE;
      }else if(fUseDistanceToBadChannel==2){
        if (( coldiff <= fMinDistanceToBadChannel )  && ( rowdiff <= fMinDistanceToBadChannel ) && (coldiff + rowdiff <= fMinDistanceToBadChannel*2)) return kTRUE;
      }
      //cout << "not within distanceToBadChannel!" << endl;
    }
  }

  //check bad cells in neighboring SM only if within chosen distanceToBadChannel from maxEnergyCell the next SM could be reached
  if(checkNextSM) {
    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2-1
    // C Side impair SM, nSupMod%2=1;A side pair SM nSupMod%2=0
    if( fClusterType == 1 || fClusterType == 4){
      if(largestCelliMod%2){
        nMinCols = largestCellicol - distanceForLoop + AliEMCALGeoParams::fgkEMCALCols;
        nMaxCols = AliEMCALGeoParams::fgkEMCALCols;

        largestCelliMod -= 1;
        largestCellicol += AliEMCALGeoParams::fgkEMCALCols;
      }else{
        nMinCols = 0;
        nMaxCols = largestCellicol + distanceForLoop - AliEMCALGeoParams::fgkEMCALCols;

        largestCelliMod += 1;
        largestCellicol -= AliEMCALGeoParams::fgkEMCALCols;
      }
    }else if( fClusterType == 2 ){
//      nMaxRows = 64;
//      nMaxCols = 56;
    }
    //cout << "largestCell: " << largestCellID << ",mod: " << largestCelliMod << ",col: " << largestCellicol << ",row: " << largestCellirow << endl;
    //cout << "distanceForLoop: " << distanceForLoop << ",nMinRows: " << nMinRows << ",nMaxRows: " << nMaxRows << ",nMinCols: " << nMinCols << ",nMaxCols: " << nMaxCols << endl;
    for (Int_t irow = nMinRows;irow < nMaxRows;irow++)
    {
      for (Int_t icol = nMinCols;icol < nMaxCols;icol++)
      {
        Int_t iBadCell = 0;
        if( (fClusterType == 1 || fClusterType == 4) && (largestCelliMod<fEMCALBadChannelsMap->GetEntries() || fEMCALBadChannelsMap1D)){
          if(fEMCALBadChannelsMap1D)
            iBadCell = (Int_t) fEMCALBadChannelsMap1D->GetBinContent(fGeomEMCAL->GetAbsCellIdFromCellIndexes(largestCelliMod, icol, irow));
          else
            iBadCell = (Int_t) ((TH2I*)fEMCALBadChannelsMap->At(largestCelliMod))->GetBinContent(icol,irow);
        }else if( fClusterType == 2 && fPHOSBadChannelsMap[largestCelliMod+1]){
          iBadCell = (Int_t) ((TH2I*)fPHOSBadChannelsMap[largestCelliMod+1])->GetBinContent(icol,irow);
        }
        //cout << "largestCelliMod: " << largestCelliMod << ",iBadCell: " << iBadCell << ",icol: " << icol << ",irow: " << irow << endl;
        if(iBadCell==0) continue;

        rowdiff = TMath::Abs( largestCellirow - irow ) ;
        coldiff = TMath::Abs( largestCellicol - icol ) ;
        //cout << "rowdiff: " << rowdiff << ",coldiff: " << coldiff << endl;
        if(fUseDistanceToBadChannel==1){
          if ((coldiff + rowdiff <= fMinDistanceToBadChannel )) return kTRUE;
        }else if(fUseDistanceToBadChannel==2){
          if (( coldiff <= fMinDistanceToBadChannel )  && ( rowdiff <= fMinDistanceToBadChannel ) && (coldiff + rowdiff <= fMinDistanceToBadChannel*2)) return kTRUE;
        }
        //cout << "not within distanceToBadChannel!" << endl;
      }
    }
  }

  return kFALSE;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::CheckDistanceToBadChannelSwapping(const Int_t CellID, Double_t phiCluster, AliVEvent* event)
{
  if(CellID < 0) return kTRUE;
  if( (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) && !fEMCALInitialized ) InitializeEMCAL(event);
  if( fClusterType == 2 && ( !fPHOSInitialized || (fPHOSCurrentRun != event->GetRunNumber()) ) ) InitializePHOS(event);
  if( fClusterType == 2 && !fGeomPHOS) fGeomPHOS = AliPHOSGeometry::GetInstance();

  Int_t iBadCell = -1;
  Int_t iSupMod, iMod, iPhi, iEta;
  if( (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) ){ // EMCal Case
    if(!fEMCALBadChannelsMap1D) fEMCALBadChannelsMap1D  = fEMCALRecUtils->GetEMCALChannelStatusMap1D();

    if(fEMCALBadChannelsMap1D){
        iBadCell      = (Int_t) fEMCALBadChannelsMap1D->GetBinContent(CellID);
    } else if(fEMCALBadChannelsMap){
        fGeomEMCAL->GetCellIndex(CellID, iSupMod, iMod, iPhi, iEta);
        if(iSupMod >= 0 && iSupMod < 20) iBadCell = (Int_t) ((TH2I*)fEMCALBadChannelsMap->At(iSupMod))->GetBinContent(iPhi,iEta);
    }

  } else if( fClusterType == 2){ // PHOS case
    Int_t relid[4];
    fGeomPHOS->AbsToRelNumbering(CellID,relid);
    if(relid[1]!=0) AliFatal("PHOS CPV in PHOS cell array (rotation background method)?");
    if(relid[0]>=fGeomPHOS->GetNModules() || relid[0]<0 || !fPHOSBadChannelsMap[relid[0]]) return kTRUE;
    iBadCell = (Int_t) ((TH2I*)fPHOSBadChannelsMap[relid[0]])->GetBinContent(relid[2],relid[3]);
  }


  if(iBadCell > 0) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Int_t  AliCaloPhotonCuts::GetCaloCellIdFromEtaPhi(const Double_t eta, const Double_t phi){
  Int_t cellId = -1;
  if(fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
    if(!fGeomEMCAL) fGeomEMCAL = AliEMCALGeometry::GetInstance();
    if(!fGeomEMCAL){ AliFatal("ERROR: EMCal geometry not initialized for cluster swapping method");}
    fGeomEMCAL->GetAbsCellIdFromEtaPhi(eta, phi, cellId);
  }
  else if(fClusterType == 2){
    if(!fPHOSGeoUtils){ AliFatal("PHOS geoUtils not initialized!");}
    Double_t tmpVtx[] = {0,0,0};
    Int_t modNr;
    Double_t x, z;
    if(fPHOSGeoUtils->ImpactOnEmc(tmpVtx, 2*atan(exp(-eta)), phi, modNr, z, x)){
      fPHOSGeoUtils->RelPosToAbsId(modNr, x, z, cellId);
    }
  }
  return cellId;

}


//________________________________________________________________________
Bool_t AliCaloPhotonCuts::ClusterIsSelected(AliVCluster *cluster, AliVEvent * event, AliMCEvent * mcEvent, Int_t isMC, Double_t weight, Long_t clusterID)
{
  //Selection of Reconstructed photon clusters with Calorimeters
  fIsAcceptedForBasic               = kFALSE;
  FillClusterCutIndex(kPhotonIn);

//  Double_t vertex[3] = {0,0,0};
//  event->GetPrimaryVertex()->GetXYZ(vertex);
    // TLorentzvector with cluster
//  TLorentzVector clusterVector;
//  cluster->GetMomentum(clusterVector,vertex);

  Float_t clusPos[3]={0,0,0};
  cluster->GetPosition(clusPos);
  TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
  Double_t etaCluster = clusterVector.Eta();
  Double_t phiCluster = clusterVector.Phi();
  if (phiCluster < 0) phiCluster += 2*TMath::Pi();

  // Histos before cuts
  if(fHistClusterEtavsPhiBeforeAcc && cluster->E()>0. ) fHistClusterEtavsPhiBeforeAcc->Fill(phiCluster,etaCluster,weight);

  // Cluster Selection - 0= accept any calo cluster
  if (fClusterType > 0){
    //Select EMCAL cluster
    if ( (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) && !cluster->IsEMCAL()){
      FillClusterCutIndex(kDetector);
      return kFALSE;
    }
    //Select PHOS cluster
    if (fClusterType == 2 && !( cluster->GetType() == AliVCluster::kPHOSNeutral)){
      FillClusterCutIndex(kDetector);
      return kFALSE;
    }
    // do NonLinearity if switched on
    if(fUseNonLinearity){
      if(fHistEnergyOfClusterBeforeNL) fHistEnergyOfClusterBeforeNL->Fill(cluster->E(),weight);
      ApplyNonLinearity(cluster,isMC,event);
      if(fHistEnergyOfClusterAfterNL) fHistEnergyOfClusterAfterNL->Fill(cluster->E(),weight);
    }
  }
  if(isMC == 0){
    ApplySMWiseEnergyCorrection(cluster, isMC, event);
  }

  // Acceptance Cuts
  if(!AcceptanceCuts(cluster,event,weight)){
    FillClusterCutIndex(kAcceptance);
    return kFALSE;
  }
  // Cluster Quality Cuts
  if(!ClusterQualityCuts(cluster,event,mcEvent,isMC,weight,clusterID)){
    FillClusterCutIndex(kClusterQuality);
    return kFALSE;
  }

  // Photon passed cuts
  FillClusterCutIndex(kPhotonOut);
  return kTRUE;
}


//________________________________________________________________________
Bool_t AliCaloPhotonCuts::AcceptanceCuts(AliVCluster *cluster, AliVEvent* event, Double_t weight)
{
   // Exclude certain areas for photon reconstruction

  Int_t cutIndex=0;
  if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);
  cutIndex++;

  Float_t clusPos[3]={0,0,0};
  cluster->GetPosition(clusPos);
  TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
  Double_t etaCluster = clusterVector.Eta();
  Double_t phiCluster = clusterVector.Phi();
  if (phiCluster < 0) phiCluster += 2*TMath::Pi();

  // check eta range
  if (fUseEtaCut){
    if(fClusterType == 4){
      // additional phi requirement needed for the DCal hole check in ClusterType 4 case (otherwise hole is also cutted in EMCal)
      if (etaCluster < fMinEtaCut || etaCluster > fMaxEtaCut ){
        if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);
        return kFALSE;
      }
    } else {
      if (etaCluster < fMinEtaCut || etaCluster > fMaxEtaCut || (fClusterType == 3 && etaCluster < fMaxEtaInnerEdge && etaCluster > fMinEtaInnerEdge ) ){
        if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);
        return kFALSE;
      }
    }
  }
  cutIndex++;

  // check phi range
  if (fUsePhiCut ){
    if(fClusterType == 4){
      if ( (phiCluster < fMinPhiCut || phiCluster > fMaxPhiCut) && (phiCluster < fMinPhiCutDMC || phiCluster > fMaxPhiCutDMC) ){
        if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);
        return kFALSE;
      }
    } else {
      if (phiCluster < fMinPhiCut || phiCluster > fMaxPhiCut){
        if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);
        return kFALSE;
      }
    }
  }
  cutIndex++;

  // check distance to bad channel
  if (fUseDistanceToBadChannel>0){
    if (CheckDistanceToBadChannel(cluster,event)){
      if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
    }
  }
  //alternatively check histogram fHistoModifyAcc if cluster should be rejected
  if(fHistoModifyAcc){
    if(fHistoModifyAcc->GetBinContent(FindLargestCellInCluster(cluster,event)) < 1){
      if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);
      return kFALSE;
    }
  }
  cutIndex++;
  if(fHistAcceptanceCuts)fHistAcceptanceCuts->Fill(cutIndex);

  // Histos after cuts
  if(fHistClusterEtavsPhiAfterAcc && cluster->E()>0. ) fHistClusterEtavsPhiAfterAcc->Fill(phiCluster,etaCluster,weight);

  return kTRUE;
}

Bool_t  AliCaloPhotonCuts::ClusterIsIsolated(Int_t clusterID, AliAODConversionPhoton *PhotonCandidate)
{

  if(fUsePhotonIsolation){
    Float_t ClusterPt = PhotonCandidate->Pt();
    Float_t pTCone = fMomPercentage*ClusterPt;
    //Float_t pTCone = 0.05*ClusterPt;
    //Float_t pTCone = 2;

    if(fCaloIsolation->GetIsolation(clusterID,fIsolationRadius,pTCone)){
      return kTRUE;
    }else{
      return kFALSE;
    }
  }else{
    return kTRUE;//if there's no isolation, all of them can be seen as isolated and pass the cut
  }

}


//________________________________________________________________________
Bool_t AliCaloPhotonCuts::MatchConvPhotonToCluster(AliAODConversionPhoton* convPhoton, AliVCluster* cluster, AliVEvent* event, Double_t weight){

  if (!fUseDistTrackToCluster || fUseElectronClusterCalibration) return kFALSE;
  if( (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) && !fEMCALInitialized ) InitializeEMCAL(event);
  if( fClusterType == 2 && ( !fPHOSInitialized || (fPHOSCurrentRun != event->GetRunNumber()) ) ) InitializePHOS(event);

  AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aodev = 0;
  if (!esdev) {
    aodev = dynamic_cast<AliAODEvent*>(event);
    if (!aodev) {
      AliError("Task needs AOD or ESD event, returning");
      return kFALSE;
    }
  }

  if(!cluster->IsEMCAL() && !cluster->IsPHOS()){AliError("Cluster is neither EMCAL nor PHOS, returning");return kFALSE;}

  Float_t clusterPosition[3] = {0,0,0};
  cluster->GetPosition(clusterPosition);
  Double_t clusterR = TMath::Sqrt( clusterPosition[0]*clusterPosition[0] + clusterPosition[1]*clusterPosition[1] );
  if(fHistClusterRBeforeQA) fHistClusterRBeforeQA->Fill(clusterR,weight);

//cout << "+++++++++ Cluster: x, y, z, R" << clusterPosition[0] << ", " << clusterPosition[1] << ", " << clusterPosition[2] << ", " << clusterR << "+++++++++" << endl;

  Bool_t matched = kFALSE;
  for (Int_t i = 0;i < 2;i++){
    Int_t tracklabel = convPhoton->GetLabel(i);
    AliVTrack *inTrack = 0x0;
    if(esdev) {
      if(tracklabel > event->GetNumberOfTracks() ) continue;
      inTrack = esdev->GetTrack(tracklabel);
    } else {
      if(((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->AreAODsRelabeled()){
        inTrack = dynamic_cast<AliVTrack*>(event->GetTrack(tracklabel));
      } else {
        for(Int_t ii=0;ii<event->GetNumberOfTracks();ii++) {
          inTrack = dynamic_cast<AliVTrack*>(event->GetTrack(ii));
          if(inTrack){
            if(inTrack->GetID() == tracklabel) {
              break;
            }
          }
        }
      }
    }

    Float_t dEta = 0;
    Float_t dPhi = 0;
    Bool_t propagated = fCaloTrackMatcher->PropagateV0TrackToClusterAndGetMatchingResidual(inTrack,cluster,event,dEta,dPhi);
    if (propagated){
      Float_t dR2 = dPhi*dPhi + dEta*dEta;
      if(fHistDistanceTrackToClusterBeforeQA)fHistDistanceTrackToClusterBeforeQA->Fill(TMath::Sqrt(dR2), weight);
      if(fHistClusterdEtadPhiBeforeQA) fHistClusterdEtadPhiBeforeQA->Fill(dEta, dPhi, weight);

      Float_t clusM02 = (Float_t) cluster->GetM02();
      Float_t clusM20 = (Float_t) cluster->GetM20();
      if(!fDoLightOutput && (fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5 )){
        if(inTrack->Charge() > 0) {
          fHistClusterdEtadPhiPosTracksBeforeQA->Fill(dEta, dPhi, weight);
          fHistClusterdPhidPtPosTracksBeforeQA->Fill(dPhi, inTrack->Pt(), weight);
          if(inTrack->P() < 0.75) fHistClusterdEtadPhiPosTracksP_000_075BeforeQA->Fill(dEta, dPhi, weight);
          else if(inTrack->P() < 1.25) fHistClusterdEtadPhiPosTracksP_075_125BeforeQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiPosTracksP_125_999BeforeQA->Fill(dEta, dPhi, weight);
        } else {
          fHistClusterdEtadPhiNegTracksBeforeQA->Fill(dEta, dPhi, weight);
          fHistClusterdPhidPtNegTracksBeforeQA->Fill(dPhi, inTrack->Pt(), weight);
          if(inTrack->P() < 0.75) fHistClusterdEtadPhiNegTracksP_000_075BeforeQA->Fill(dEta, dPhi, weight);
          else if(inTrack->P() < 1.25) fHistClusterdEtadPhiNegTracksP_075_125BeforeQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiNegTracksP_125_999BeforeQA->Fill(dEta, dPhi, weight);
        }
        fHistClusterdEtadPtBeforeQA->Fill(dEta, inTrack->Pt(), weight);
        fHistClusterM02M20BeforeQA->Fill(clusM02, clusM20, weight);
        if(fCurrentMC != kNoMC && fIsMC > 0){
          Int_t clusterMCLabel = cluster->GetLabel();
          Int_t convPhotonDaughterLabel = -1;
          if(inTrack->Charge() > 0) convPhotonDaughterLabel = convPhoton->GetMCLabelPositive();
          else convPhotonDaughterLabel = convPhoton->GetMCLabelNegative();
          if( (convPhotonDaughterLabel != -1) && (clusterMCLabel != -1) && (convPhotonDaughterLabel == clusterMCLabel)){ //real match
            fHistClusterdEtadPtTrueMatched->Fill(dEta, inTrack->Pt(), weight);
            if(inTrack->Charge() > 0) fHistClusterdPhidPtPosTracksTrueMatched->Fill(dPhi, inTrack->Pt(), weight);
            else fHistClusterdPhidPtNegTracksTrueMatched->Fill(dPhi, inTrack->Pt(), weight);
          }
        }
      }

      Bool_t match_dEta = (TMath::Abs(dEta) < fMaxDistTrackToClusterEta) ? kTRUE : kFALSE;
      Bool_t match_dPhi = kFALSE;
      if( (inTrack->Charge() > 0) && (dPhi > fMinDistTrackToClusterPhi) && (dPhi < fMaxDistTrackToClusterPhi) ) match_dPhi = kTRUE;
      else if( (inTrack->Charge() < 0) && (dPhi < -fMinDistTrackToClusterPhi) && (dPhi > -fMaxDistTrackToClusterPhi) ) match_dPhi = kTRUE;

      if(fUsePtDepTrackToCluster == 1){
        if( TMath::Abs(dEta) < fFuncPtDepEta->Eval(inTrack->Pt())) match_dEta = kTRUE;
        else match_dEta = kFALSE;

        if( TMath::Abs(dPhi) < fFuncPtDepPhi->Eval(inTrack->Pt())) match_dPhi = kTRUE;
        else match_dPhi = kFALSE;
      }
//
      if(match_dEta && match_dPhi){
            //if(dR2 < fMinDistTrackToCluster*fMinDistTrackToCluster){
        matched = kTRUE;
        if(fHistClusterdEtadPtAfterQA) fHistClusterdEtadPtAfterQA->Fill(dEta,inTrack->Pt());
        if(fHistClusterdPhidPtAfterQA) fHistClusterdPhidPtAfterQA->Fill(dPhi,inTrack->Pt());
      } else {
        if(fHistDistanceTrackToClusterAfterQA)fHistDistanceTrackToClusterAfterQA->Fill(TMath::Sqrt(dR2), weight);
        if(fHistClusterdEtadPhiAfterQA) fHistClusterdEtadPhiAfterQA->Fill(dEta, dPhi, weight);
        if(fHistClusterRAfterQA) fHistClusterRAfterQA->Fill(clusterR, weight);
        if(!fDoLightOutput && (fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5 )){
          if(inTrack->Charge() > 0) fHistClusterdEtadPhiPosTracksAfterQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiNegTracksAfterQA->Fill(dEta, dPhi, weight);
          fHistClusterM02M20AfterQA->Fill(clusM02, clusM20, weight);
        }
      }
    }
  }

  return matched;

}

//________________________________________________________________________
void AliCaloPhotonCuts::MatchElectronTracksToClusters(AliVEvent* event, AliMCEvent* MCevent, AliVCluster* cluster, Int_t isMC, vector<Int_t> vElectronTracks, Double_t weight){

  if (!cluster){
    return;
  }
  if( (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) && !fEMCALInitialized ) InitializeEMCAL(event);
  if( fClusterType == 2 && ( !fPHOSInitialized || (fPHOSCurrentRun != event->GetRunNumber()) ) ) InitializePHOS(event);

//  Int_t nClus = 0;
//  TClonesArray * arrClustersMatch = NULL;
//  if(!fCorrTaskSetting.CompareTo("")){
//    nClus = event->GetNumberOfCaloClusters();
//  } else {
//    arrClustersMatch = dynamic_cast<TClonesArray*>(event->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
//    if(!arrClustersMatch)
//      AliFatal(Form("%sClustersBranch was not found in AliCaloPhotonCuts::FillHistogramsExtendedQA! Check the correction framework settings!",fCorrTaskSetting.Data()));
//    nClus = arrClustersMatch->GetEntries();
//  }

  if(fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
    fGeomEMCAL = AliEMCALGeometry::GetInstance();
    if(!fGeomEMCAL){ AliFatal("EMCal geometry not initialized!");}
    //nModules = fGeomEMCAL->GetNumberOfSuperModules();
  }else if(fClusterType == 2){
    fGeomPHOS = AliPHOSGeometry::GetInstance();
    if(!fGeomPHOS){ AliFatal("PHOS geometry not initialized!");}
    //nModules = fGeomPHOS->GetNModules();
  }

  AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aodev = 0;
  if (!esdev) {
    aodev = dynamic_cast<AliAODEvent*>(event);
    if (!aodev) {
      AliError("Task needs AOD or ESD event, returning");
      return;
    }
  }

  //loop over all electron candidates
  for (UInt_t itr=0;itr<vElectronTracks.size();itr++){
    AliVTrack *inTrack = 0x0;
    if(esdev){
      inTrack = esdev->GetTrack(vElectronTracks[itr]);
      if(!inTrack) continue;
//      AliESDtrack *esdt = dynamic_cast<AliESDtrack*>(inTrack);
//      const AliExternalTrackParam *in = esdt->GetInnerParam();
//      if (!in){AliDebug(2, "Could not get InnerParam of Track, continue");continue;}
    } else if(aodev) {
      inTrack = dynamic_cast<AliVTrack*>(aodev->GetTrack(vElectronTracks[itr]));
      if(!inTrack) {cout<<"track not valid..."<<endl; continue;}
//      AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(inTrack);

    }

    Float_t dEta, dPhi;
    Float_t clsPos[3] = {0.,0.,0.};
    if(!fCaloTrackMatcher->GetTrackClusterMatchingResidual(inTrack->GetID(),cluster->GetID(),dEta,dPhi)){
      if(!fCaloTrackMatcher->PropagateV0TrackToClusterAndGetMatchingResidual(inTrack, cluster, event, dEta, dPhi)){
        continue;
      }
    }

    cluster->GetPosition(clsPos);
//    Float_t clusterR = TMath::Sqrt( clsPos[0]*clsPos[0] + clsPos[1]*clsPos[1] );
//    Float_t dR2 = dPhi*dPhi + dEta*dEta;

    Bool_t match_dEta = (TMath::Abs(dEta) < fMaxDistTrackToClusterEta) ? kTRUE : kFALSE;
    Bool_t match_dPhi = kFALSE;
//    Bool_t vetoEOverP = kFALSE;

    if( (inTrack->Charge() > 0) && (dPhi > fMinDistTrackToClusterPhi) && (dPhi < fMaxDistTrackToClusterPhi) ) match_dPhi = kTRUE;
    else if( (inTrack->Charge() < 0) && (dPhi < -fMinDistTrackToClusterPhi) && (dPhi > -fMaxDistTrackToClusterPhi) ) match_dPhi = kTRUE;

    if(fUsePtDepTrackToCluster == 1){
      if( TMath::Abs(dEta) < fFuncPtDepEta->Eval(inTrack->Pt())) match_dEta = kTRUE;
      else match_dEta = kFALSE;

      if( TMath::Abs(dPhi) < fFuncPtDepPhi->Eval(inTrack->Pt())) match_dPhi = kTRUE;
      else match_dPhi = kFALSE;
    }

    if(match_dEta && match_dPhi){
      if(fExtendedMatchAndQA > 1){
        if(inTrack->Charge() < 0){
          fHistElectronClusterMatch->Fill(cluster->E(), inTrack->GetTrackPOnEMCal(), weight);
        } else {
          fHistPositronClusterMatch->Fill(cluster->E(), inTrack->GetTrackPOnEMCal(), weight);
        }
      }
      fHistElectronPositronClusterMatch->Fill(cluster->E(), inTrack->GetTrackPOnEMCal(), weight);
      fHistElectronPositronClusterMatchSub->Fill(cluster->E(), cluster->E() - inTrack->GetTrackPOnEMCal(), weight);

      if(isMC){
        if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
        if (fAODMCTrackArray == NULL){
          AliError("No MC particle list available in AOD");
          return;
        }
        Int_t tmpLabel = (Int_t) ((AliAODTrack*)inTrack)->GetLabel();
        if(tmpLabel > 0){
          AliAODMCParticle* trackPart    = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(tmpLabel));
          if(!trackPart) continue;

          if(TMath::Abs(trackPart->GetPdgCode()) == 11){
            fHistElectronClusterMatchTruePID->Fill(0.5, trackPart->P(), weight);
            fHistTrueElectronPositronClusterMatch->Fill(cluster->E(), inTrack->GetTrackPOnEMCal(), weight);
          } else if(TMath::Abs(trackPart->GetPdgCode()) == 211){
            fHistElectronClusterMatchTruePID->Fill(1.5, trackPart->P(), weight);
            fHistTrueNoElectronPositronClusterMatch->Fill(cluster->E(), inTrack->GetTrackPOnEMCal(), weight);
          } else if(TMath::Abs(trackPart->GetPdgCode()) == 2212){
            fHistElectronClusterMatchTruePID->Fill(2.5, trackPart->P(), weight);
            fHistTrueNoElectronPositronClusterMatch->Fill(cluster->E(), inTrack->GetTrackPOnEMCal(), weight);
          } else if(TMath::Abs(trackPart->GetPdgCode()) == 321){
            fHistElectronClusterMatchTruePID->Fill(3.5, trackPart->P(), weight);
            fHistTrueNoElectronPositronClusterMatch->Fill(cluster->E(), inTrack->GetTrackPOnEMCal(), weight);
          } else {
            fHistElectronClusterMatchTruePID->Fill(4.5, trackPart->P(), weight);
            fHistTrueNoElectronPositronClusterMatch->Fill(cluster->E(), inTrack->GetTrackPOnEMCal(), weight);
          }
        }
      }
    }
  }
}

//________________________________________________________________________
void AliCaloPhotonCuts::MatchTracksToClusters(AliVEvent* event, Double_t weight, Bool_t isEMCalOnly, AliMCEvent* mcEvent){
  if( !fUseDistTrackToCluster || fUseElectronClusterCalibration) return;
  if( (fClusterType == 1 || fClusterType == 3 || fClusterType == 4) && !fEMCALInitialized ) InitializeEMCAL(event);
  if( fClusterType == 2 && ( !fPHOSInitialized || (fPHOSCurrentRun != event->GetRunNumber()) ) ) InitializePHOS(event);

  // not yet fully implemented + tested for PHOS
  // if( fClusterType != 1 && fClusterType != 3) return;

  fVectorMatchedClusterIDs.clear();

  Int_t nClus = 0;
  TClonesArray * arrClustersMatch = NULL;
  if(!fCorrTaskSetting.CompareTo("")){
    nClus = event->GetNumberOfCaloClusters();
  } else {
    arrClustersMatch = dynamic_cast<TClonesArray*>(event->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!arrClustersMatch)
      AliFatal(Form("%sClustersBranch was not found in AliCaloPhotonCuts::FillHistogramsExtendedQA! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nClus = arrClustersMatch->GetEntries();
  }
  //Int_t nModules = 0;

  if(fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
    fGeomEMCAL = AliEMCALGeometry::GetInstance();
    if(!fGeomEMCAL){ AliFatal("EMCal geometry not initialized!");}
    //nModules = fGeomEMCAL->GetNumberOfSuperModules();
  }else if(fClusterType == 2){
    fGeomPHOS = AliPHOSGeometry::GetInstance();
    if(!fGeomPHOS){ AliFatal("PHOS geometry not initialized!");}
    //nModules = fGeomPHOS->GetNModules();
  }

  AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aodev = 0;
  if (!esdev) {
    aodev = dynamic_cast<AliAODEvent*>(event);
    if (!aodev) {
      AliError("Task needs AOD or ESD event, returning");
      return;
    }
  }

  // if not EMCal only reconstruction (-> hybrid PCM+EMCal), use only primary tracks for basic track matching procedure
  std::unique_ptr<AliESDtrackCuts> EsdTrackCuts;
  if(!isEMCalOnly && esdev){
    // Using standard function for setting Cuts
    Int_t runNumber = event->GetRunNumber();
    // if LHC11a or earlier or if LHC13g or if LHC12a-i -> use 2010 cuts
    if( (runNumber<=146860) || (runNumber>=197470 && runNumber<=197692) || (runNumber>=172440 && runNumber<=193766) ){
      EsdTrackCuts = std::unique_ptr<AliESDtrackCuts>(AliESDtrackCuts::GetStandardITSTPCTrackCuts2010());
    // else if run2 data use 2015 PbPb cuts
    }else if (runNumber>=209122){
      // EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb();
      // hard coded track cuts for the moment, because AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb() gives spams warnings
      EsdTrackCuts = std::unique_ptr<AliESDtrackCuts>(new  AliESDtrackCuts());
      EsdTrackCuts->AliESDtrackCuts::SetMinNCrossedRowsTPC(70);
      EsdTrackCuts->AliESDtrackCuts::SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
      EsdTrackCuts->AliESDtrackCuts::SetCutOutDistortedRegionsTPC(kTRUE);
      EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterTPC(4);
      EsdTrackCuts->AliESDtrackCuts::SetAcceptKinkDaughters(kFALSE);
      EsdTrackCuts->AliESDtrackCuts::SetRequireTPCRefit(kTRUE);
      // ITS
      EsdTrackCuts->AliESDtrackCuts::SetRequireITSRefit(kTRUE);
      EsdTrackCuts->AliESDtrackCuts::SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                                              AliESDtrackCuts::kAny);
      EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
      EsdTrackCuts->AliESDtrackCuts::SetMaxChi2TPCConstrainedGlobal(36);
      EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexZ(2);
      EsdTrackCuts->AliESDtrackCuts::SetDCAToVertex2D(kFALSE);
      EsdTrackCuts->AliESDtrackCuts::SetRequireSigmaToVertex(kFALSE);
      EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterITS(36);
    // else use 2011 version of track cuts
    }else{
      EsdTrackCuts = std::unique_ptr<AliESDtrackCuts>(AliESDtrackCuts::GetStandardITSTPCTrackCuts2011());
    }
    EsdTrackCuts->SetMaxDCAToVertexZ(2);
    EsdTrackCuts->SetEtaRange(-0.8, 0.8);
    EsdTrackCuts->SetPtRange(0.15);
  }

//  cout << "MatchTracksToClusters: " << event->GetNumberOfTracks() << ", " << fIsPureCalo << ", " << fUseDistTrackToCluster << endl;

  for (Int_t itr=0;itr<event->GetNumberOfTracks();itr++){
    AliVTrack *inTrack = 0x0;
    if(esdev){
      inTrack = esdev->GetTrack(itr);
      if(!inTrack) continue;
      AliESDtrack *esdt = dynamic_cast<AliESDtrack*>(inTrack);
      if(!isEMCalOnly){ //match only primaries for hybrid reconstruction schemes
         if(!EsdTrackCuts->AcceptTrack(esdt)) continue;
      }
      const AliExternalTrackParam *in = esdt->GetInnerParam();
      if (!in){AliDebug(2, "Could not get InnerParam of Track, continue");continue;}
    } else if(aodev) {
      inTrack = dynamic_cast<AliVTrack*>(aodev->GetTrack(itr));
      if(!inTrack) continue;
      AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(inTrack);
      if(!isEMCalOnly){ //match only primaries for hybrid reconstruction schemes
        if(!aodt->IsHybridGlobalConstrainedGlobal()) continue;
        if(TMath::Abs(aodt->Eta())>0.8) continue;
        if(aodt->Pt()<0.15) continue;
      }

    }

    Float_t clsPos[3] = {0.,0.,0.};
    for(Int_t iclus=0;iclus < nClus;iclus++){
      AliVCluster * cluster = NULL;
      std::unique_ptr<AliVCluster> tmpcluster;
      if(arrClustersMatch){
        if(esdev){
          tmpcluster = std::unique_ptr<AliVCluster>(new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersMatch->At(iclus)));
          cluster = tmpcluster.get();
        }
        else if(aodev){
          tmpcluster = std::unique_ptr<AliVCluster>(new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersMatch->At(iclus)));
          cluster = tmpcluster.get();
        }
      } else {
        cluster = event->GetCaloCluster(iclus);
      }

      if (!cluster){
        continue;
      }
      Float_t dEta, dPhi;
      if(!fCaloTrackMatcher->GetTrackClusterMatchingResidual(inTrack->GetID(),cluster->GetID(),dEta,dPhi)){
        continue;
      }
      cluster->GetPosition(clsPos);
      Float_t clusterR = TMath::Sqrt( clsPos[0]*clsPos[0] + clsPos[1]*clsPos[1] );
      Float_t dR2 = dPhi*dPhi + dEta*dEta;

      if(isEMCalOnly && fHistDistanceTrackToClusterBeforeQA)fHistDistanceTrackToClusterBeforeQA->Fill(TMath::Sqrt(dR2), weight);
      if(isEMCalOnly && fHistClusterdEtadPhiBeforeQA) fHistClusterdEtadPhiBeforeQA->Fill(dEta, dPhi, weight);
      if(isEMCalOnly && fHistClusterRBeforeQA) fHistClusterRBeforeQA->Fill(clusterR, weight);

      Float_t clusM02 = (Float_t) cluster->GetM02();
      Float_t clusM20 = (Float_t) cluster->GetM20();
      if(isEMCalOnly && !fDoLightOutput && (fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5 )){
        if(inTrack->Charge() > 0) {
          fHistClusterdEtadPhiPosTracksBeforeQA->Fill(dEta, dPhi, weight);
          fHistClusterdPhidPtPosTracksBeforeQA->Fill(dPhi, inTrack->Pt(), weight);
          if(inTrack->P() < 0.75) fHistClusterdEtadPhiPosTracksP_000_075BeforeQA->Fill(dEta, dPhi, weight);
          else if(inTrack->P() < 1.25) fHistClusterdEtadPhiPosTracksP_075_125BeforeQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiPosTracksP_125_999BeforeQA->Fill(dEta, dPhi, weight);
        }
        else{
          fHistClusterdEtadPhiNegTracksBeforeQA->Fill(dEta, dPhi, weight);
          fHistClusterdPhidPtNegTracksBeforeQA->Fill(dPhi, inTrack->Pt(), weight);
          if(inTrack->P() < 0.75) fHistClusterdEtadPhiNegTracksP_000_075BeforeQA->Fill(dEta, dPhi, weight);
          else if(inTrack->P() < 1.25) fHistClusterdEtadPhiNegTracksP_075_125BeforeQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiNegTracksP_125_999BeforeQA->Fill(dEta, dPhi, weight);
        }
        fHistClusterdEtadPtBeforeQA->Fill(dEta, inTrack->Pt(), weight);
        fHistClusterM02M20BeforeQA->Fill(clusM02, clusM20, weight);
      }

      Bool_t match_dEta = (TMath::Abs(dEta) < fMaxDistTrackToClusterEta) ? kTRUE : kFALSE;
      Bool_t match_dPhi = kFALSE;
      Bool_t vetoEOverP = kFALSE;

      if( (inTrack->Charge() > 0) && (dPhi > fMinDistTrackToClusterPhi) && (dPhi < fMaxDistTrackToClusterPhi) ) match_dPhi = kTRUE;
      else if( (inTrack->Charge() < 0) && (dPhi < -fMinDistTrackToClusterPhi) && (dPhi > -fMaxDistTrackToClusterPhi) ) match_dPhi = kTRUE;

      if(fUsePtDepTrackToCluster == 1){
        if( TMath::Abs(dEta) < fFuncPtDepEta->Eval(inTrack->Pt())) match_dEta = kTRUE;
        else match_dEta = kFALSE;

        if( TMath::Abs(dPhi) < fFuncPtDepPhi->Eval(inTrack->Pt())) match_dPhi = kTRUE;
        else match_dPhi = kFALSE;
      }

      if(fUseEOverPVetoTM && cluster->E()/inTrack->P() > fEOverPMax)
        vetoEOverP = kTRUE;

      if(match_dEta && match_dPhi){
        if(fUseEOverPVetoTM){
          if(!vetoEOverP){
            fVectorMatchedClusterIDs.push_back(cluster->GetID());
            if(fHistMatchedTrackPClusEAfterEOverPVeto) fHistMatchedTrackPClusEAfterEOverPVeto->Fill(cluster->E(),inTrack->P());
          }
        }else if(fUseTMMIPsubtraction){
          //Subtracting the MIP energy is there is a match
          cluster->SetE(cluster->E()-0.290);
        }else{
          fVectorMatchedClusterIDs.push_back(cluster->GetID());
        }
        if(fHistMatchedTrackPClusE && fUseEOverPVetoTM) fHistMatchedTrackPClusE->Fill(cluster->E(),inTrack->P());
        if(fIsMC && (fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5 ) && fUseEOverPVetoTM){
          if(IsClusterPi0(event, mcEvent, cluster) && fHistMatchedTrackPClusETruePi0Clus)
            fHistMatchedTrackPClusETruePi0Clus->Fill(cluster->E(),inTrack->P());
        }

        if(isEMCalOnly){
          if(fHistClusterdEtadPtAfterQA) fHistClusterdEtadPtAfterQA->Fill(dEta,inTrack->Pt());
          if(fHistClusterdPhidPtAfterQA) fHistClusterdPhidPtAfterQA->Fill(dPhi,inTrack->Pt());
        }
        if(!fUseTMMIPsubtraction) break;
      } else if(isEMCalOnly){
        if(fHistDistanceTrackToClusterAfterQA)fHistDistanceTrackToClusterAfterQA->Fill(TMath::Sqrt(dR2), weight);
        if(fHistClusterdEtadPhiAfterQA) fHistClusterdEtadPhiAfterQA->Fill(dEta, dPhi, weight);
        if(fHistClusterRAfterQA) fHistClusterRAfterQA->Fill(clusterR, weight);
        if(!fDoLightOutput && (fExtendedMatchAndQA == 1 || fExtendedMatchAndQA == 3 || fExtendedMatchAndQA == 5 )){
          if(inTrack->Charge() > 0) fHistClusterdEtadPhiPosTracksAfterQA->Fill(dEta, dPhi, weight);
          else fHistClusterdEtadPhiNegTracksAfterQA->Fill(dEta, dPhi, weight);
          fHistClusterM02M20AfterQA->Fill(clusM02, clusM20, weight);
        }
        // cout << "no match" << endl;
      }
    }
  }
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::CheckClusterForTrackMatch(AliVCluster* cluster){
  vector<Int_t>::iterator it;
  it = find (fVectorMatchedClusterIDs.begin(), fVectorMatchedClusterIDs.end(), cluster->GetID());
  if (it != fVectorMatchedClusterIDs.end()) return kTRUE;
  else return kFALSE;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::UpdateCutString() {
   ///Update the cut string (if it has been created yet)

   if(fCutString && fCutString->GetString().Length() == kNCuts) {
      fCutString->SetString(GetCutNumber());
   } else {
      return kFALSE;
   }
   return kTRUE;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::InitializeCutsFromCutString(const TString analysisCutSelection ) {
  fCutStringRead = Form("%s",analysisCutSelection.Data());

  // Initialize Cuts from a given Cut string
  AliInfo(Form("Set CaloCut Number: %s",analysisCutSelection.Data()));
  if(analysisCutSelection.Length()!=kNCuts) {
    AliError(Form("Cut selection has the wrong length! size is %d, number of cuts is %d", analysisCutSelection.Length(), kNCuts));
    return kFALSE;
  }
  if(!analysisCutSelection.IsAlnum()){
    AliError("Cut selection is not alphanumeric");
    return kFALSE;
  }

  TString analysisCutSelectionLowerCase = Form("%s",analysisCutSelection.Data());
  analysisCutSelectionLowerCase.ToLower();
  const char *cutSelection = analysisCutSelectionLowerCase.Data();
  #define ASSIGNARRAY(i)  fCuts[i] = ((int)cutSelection[i]>=(int)'a') ? cutSelection[i]-'a'+10 : cutSelection[i]-'0'
  for(Int_t ii=0;ii<kNCuts;ii++){
    ASSIGNARRAY(ii);
  }

  // Set Individual Cuts
  for(Int_t ii=0;ii<kNCuts;ii++){
    if(!SetCut(cutIds(ii),fCuts[ii]))return kFALSE;
  }

  PrintCutsWithValues(analysisCutSelection);
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliCaloPhotonCuts::SetCut(cutIds cutID, const Int_t value) {
  ///Set individual cut ID

  switch (cutID) {

    case kClusterType:
      if( SetClusterTypeCut(value)) {
        fCuts[kClusterType] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kEtaMin:
      if( SetMinEtaCut(value)) {
        fCuts[kEtaMin] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kEtaMax:
      if( SetMaxEtaCut(value)) {
        fCuts[kEtaMax] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kPhiMin:
      if( SetMinPhiCut(value)) {
        fCuts[kPhiMin] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kPhiMax:
      if( SetMaxPhiCut(value)) {
        fCuts[kPhiMax] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kDistanceToBadChannel:
      if( SetDistanceToBadChannelCut(value)) {
        fCuts[kDistanceToBadChannel] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kTiming:
      if( SetTimingCut(value)) {
        fCuts[kTiming] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kTrackMatching:
      if( SetTrackMatchingCut(value)) {
        fCuts[kTrackMatching] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kExoticCluster:
      if( SetExoticClusterCut(value)) {
        fCuts[kExoticCluster] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kMinEnergy:
      if( SetMinEnergyCut(value)) {
        fCuts[kMinEnergy] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kNMinCells:
      if( SetMinNCellsCut(value)) {
        fCuts[kNMinCells] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kMinM02:
      if( SetMinM02(value)) {
        fCuts[kMinM02] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kMaxM02:
      if( SetMaxM02(value)) {
        fCuts[kMaxM02] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kMinMaxM20:
      if( SetMinMaxM20(value)) {
        fCuts[kMinMaxM20] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kRecConv:
      if( SetRecConv(value)) {
        fCuts[kRecConv] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kDispersion:
      if( SetDispersion(value)) {
        fCuts[kDispersion] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kNLM:
      if( SetNLM(value)) {
        fCuts[kNLM] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kNonLinearity1:
      if( SetNonLinearity1(value)) {
        fCuts[kNonLinearity1] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kNonLinearity2:
      if( SetNonLinearity2(value)) {
        fCuts[kNonLinearity2] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kNCuts:
      AliError("Cut id out of range");
      return kFALSE;
  }

  AliError("Cut id %d not recognized");
  return kFALSE;


}

//________________________________________________________________________
void AliCaloPhotonCuts::PrintCuts() {
   // Print out current Cut Selection
   for(Int_t ic = 0;ic < kNCuts;ic++) {
      printf("%-30s : %d \n", fgkCutNames[ic], fCuts[ic]);
   }
}

//________________________________________________________________________
void AliCaloPhotonCuts::PrintCutsWithValues(const TString analysisCutSelection) {
  // Print out current Cut Selection with value
  printf("\nCluster cutnumber \n %s", analysisCutSelection.Data());
  printf("\n\n");
  if (fIsPureCalo>0) printf("Merged cluster analysis was specified, mode: '%i'\n", fIsPureCalo);

  printf("Acceptance cuts: \n");
  if (fClusterType == 0) printf("\tall calorimeter clusters are used\n");
  if (fClusterType == 1) printf("\tEMCAL calorimeter clusters are used\n");
  if (fClusterType == 2) printf("\tPHOS calorimeter clusters are used\n");
  if (fClusterType == 3) printf("\tDCAL calorimeter clusters are used\n");
  if (fClusterType == 4) printf("\tEMCAL and DCAL calorimeter clusters are used together\n");
  if (fUseEtaCut) printf("\t%3.2f < eta_{cluster} < %3.2f\n", fMinEtaCut, fMaxEtaCut );
  if (fUsePhiCut) printf("\t%3.2f < phi_{cluster} < %3.2f\n", fMinPhiCut, fMaxPhiCut );
  if (fUsePhiCut && fClusterType == 4) printf("\t%3.2f < phi_{cluster}^{DCAL} < %3.2f\n", fMinPhiCutDMC, fMaxPhiCutDMC );
  if (fUseDistanceToBadChannel>0) printf("\tdistance to bad channel used in mode '%i', distance in cells: %f \n",fUseDistanceToBadChannel, fMinDistanceToBadChannel);

  if (fUsePhotonIsolation){
    printf("PhotonIsolation Cuts: \n");
    if (fClusterType == 1) printf("\tEMCAL calorimeter clusters are used\n");
    if (fUsePhotonIsolation) printf("\tPhotonIsolation is turned on\n");
    if (fIsolationRadius < 0.11 && fIsolationRadius > 0.09) printf("\tIsolation Radius = 0.1\n");
    if (fIsolationRadius < 0.21 && fIsolationRadius > 0.19) printf("\tIsolation Radius = 0.2\n");
    if (fIsolationRadius < 0.31 && fIsolationRadius > 0.29) printf("\tIsolation Radius = 0.3\n");
    if (fIsolationRadius < 0.41 && fIsolationRadius > 0.39) printf("\tIsolation Radius = 0.4\n");
  }

  printf("Cluster Quality cuts: \n");
  if (fUseTimeDiff) printf("\t %6.2f ns < time difference < %6.2f ns\n", fMinTimeDiff*1e9, fMaxTimeDiff*1e9 );
  if ((fUseTimeDiff)&&(fUseTimingEfficiencyMCSimCluster==2)) printf("\t %6.2f ns < time difference HighPt < %6.2f ns\n", fMinTimeDiffHighPt*1e9, fMaxTimeDiffHighPt*1e9 );
  if (fUseDistTrackToCluster) printf("\tmin distance to track in eta > %3.2f, min phi < %3.2f and max phi > %3.2f\n", fMaxDistTrackToClusterEta, fMinDistTrackToClusterPhi, fMaxDistTrackToClusterPhi );
  if (fUseExoticCluster && fUseExoticCluster != 3)printf("\t exotic cluster: %3.2f\n", fExoticEnergyFracCluster );
  if (fUseExoticCluster == 2)printf("\t exotic cluster above: %3.2f in same T-Card\n", fExoticMinEnergyTCard );
  if (fUseExoticCluster == 3)printf("\t exotic cluster rejection from correction framework\n" );
  if (fUseMinEnergy)printf("\t E_{cluster} > %3.2f\n", fMinEnergy );
  if (fUseNCells) printf("\t number of cells per cluster >= %d\n", fMinNCells );
  if (fUseM02 == 1) printf("\t %3.2f < M02 < %3.2f\n", fMinM02, fMaxM02 );
  if (fUseM02 == 2) printf("\t energy dependent M02 cut used with cutnumber min: %d  max: %d \n", fMinM02CutNr, fMaxM02CutNr );
  if (fUseM20) printf("\t %3.2f < M20 < %3.2f\n", fMinM20, fMaxM20 );
  if (fUseRecConv) printf("\t recovering conversions for  Mgg < %3.3f\n", fMaxMGGRecConv );
  if (fUseDispersion) printf("\t dispersion < %3.2f\n", fMaxDispersion );
  if (fUseNLM) printf("\t %d < NLM < %d\n", fMinNLM, fMaxNLM );
  printf("Correction Task Setting: %s \n",fCorrTaskSetting.Data());

  printf("NonLinearity Correction: \n");
  printf("VO Reader name: %s \n",fV0ReaderName.Data());
  TString namePeriod = ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->GetPeriodName();
  if (namePeriod.CompareTo("") != 0) fCurrentMC = FindEnumForMCSet(namePeriod);
  if (fUseNonLinearity) printf("\t Chose NonLinearity cut '%i', Period name: %s, MCSet: %i \n", fSwitchNonLinearity, namePeriod.Data(), fCurrentMC );
  else printf("\t No NonLinearity Correction on AnalysisTask level has been chosen\n");

}

// EMCAL acceptance 2011
// 1.39626, 3.125 (phi)
// -0.66687,,0.66465


//________________________________________________________________________
Bool_t AliCaloPhotonCuts::SetClusterTypeCut(Int_t clusterType)
{   // Set Cut
  switch(clusterType){
  case 0: // all clusters
    fClusterType=0;
    break;
  case 1: // EMCAL clusters
    fClusterType=1;
    break;
  case 2: // PHOS clusters
    fClusterType=2;
    break;
  case 3: // DCAL clusters
    fClusterType=3;
    break;
  case 4: // EMCAL and DCAL clusters
    fClusterType=4;
    break;
  case 11: //EMCAL clusters with isolation R=0.1 and pTCone=0.1*ET_Cluster, accessible via "b"
    fClusterType=1;
    fIsolationRadius=0.1;
    fMomPercentage=0.1;
    fUsePhotonIsolation=kTRUE;
    break;
  case 12: //EMCAL clusters with isolation R=0.2 and pTCone=0.1*ET_Cluster, accessible via "c"
    fClusterType=1;
    fIsolationRadius=0.2;
    fMomPercentage=0.1;
    fUsePhotonIsolation=kTRUE;
    break;
  case 13: //EMCAL clusters with isolation R=0.3 and pTCone=0.1*ET_Cluster, accessible via "d"
    fClusterType=1;
    fIsolationRadius=0.3;
    fMomPercentage=0.1;
    fUsePhotonIsolation=kTRUE;
    break;
  case 14: //EMCAL clusters with isolation R=0.4 and pTCone=0.1*ET_Cluster, accessible via "e"
    fClusterType=1;
    fIsolationRadius=0.4;
    fMomPercentage=0.1;
    fUsePhotonIsolation=kTRUE;
    break;
  default:
    AliError(Form("ClusterTypeCut not defined %d",clusterType));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinEtaCut(Int_t minEta)
{
  switch(minEta){
  case 0:
    if (!fUseEtaCut) fUseEtaCut=0;
    fMinEtaCut=-10.;
    break;
  case 1:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut=-0.6687;
    break;
  case 2:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut=-0.5;
    break;
  case 3:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut=-2;
    break;
  case 4:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut = -0.13;
    break;
  case 5:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut=-0.7;
    break;
  case 6:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut=-0.3;
    break;
  case 7:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut=-0.4;
    break;
  case 8:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut=-0.66112;
    fMinEtaInnerEdge=-0.227579;
    break;
  case 9:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut=-0.6687; // use EMCal cut also for DCal
    fMinEtaInnerEdge=-0.227579; // DCal hole
    break;
  case 10:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMinEtaCut=0.;
    break;
  default:
    AliError(Form("MinEta Cut not defined %d",minEta));
    return kFALSE;
  }
  return kTRUE;
}


//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMaxEtaCut(Int_t maxEta)
{
  switch(maxEta){
  case 0:
    if (!fUseEtaCut) fUseEtaCut=0;
    fMaxEtaCut=10;
    break;
  case 1:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut=0.66465;
    break;
  case 2:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut=0.5;
    break;
  case 3:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut=2;
    break;
  case 4:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut= 0.13;
    break;
  case 5:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut=0.7;
    break;
  case 6:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut=0.3;
    break;
  case 7:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut=0.4;
    break;
  case 8:
    if(!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut=0.66112;
    fMaxEtaInnerEdge=0.227579;
    break;
  case 9:
    if(!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut=0.66465; // use EMCal cut also for DCal
    fMaxEtaInnerEdge=0.227579; // DCal hole
    break;
  case 10:
    if (!fUseEtaCut) fUseEtaCut=1;
    fMaxEtaCut=0.0;
    break;
  default:
    AliError(Form("MaxEta Cut not defined %d",maxEta));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinPhiCut(Int_t minPhi)
{
  switch(minPhi){
  case 0:
    if (!fUsePhiCut) fUsePhiCut=0;
    fMinPhiCut=-10000;
    break;
  case 1: // min EMCAL
    if (!fUsePhiCut) fUsePhiCut=1;
    fMinPhiCut=1.39626;
    break;
  case 2: // min EMCAL with TRD 2012
    if (!fUsePhiCut) fUsePhiCut=1;
    fMinPhiCut=2.10;
    break;
  case 3: // min EMCAL with TRD 2011
    if (!fUsePhiCut) fUsePhiCut=1;
    fMinPhiCut=2.45;
    break;
  case 4:
    if( !fUsePhiCut ) fUsePhiCut=1;
    fMinPhiCut = 4.54;//PHOS acceptance
    break;
  case 5:
    if( !fUsePhiCut ) fUsePhiCut=1;
    fMinPhiCut = 4.5572;//DCal acceptance
    break;
  case 6:
    if( !fUsePhiCut ) fUsePhiCut=1;
    fMinPhiCut = 4.36;//PHOS acceptance RUN2
    break;
  case 7: // min EMCAL
    if (!fUsePhiCut) fUsePhiCut=1;
    fMinPhiCut=1.39626;
    fMinPhiCutDMC=4.5572;
    break;

  default:
    AliError(Form("MinPhi Cut not defined %d",minPhi));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMaxPhiCut(Int_t maxPhi)
{
  switch(maxPhi){
  case 0:
    if (!fUsePhiCut) fUsePhiCut=0;
    fMaxPhiCut=10000;
    break;
  case 1: // max EMCAL
    if (!fUsePhiCut) fUsePhiCut=1;
    fMaxPhiCut=3.15;
    break;
  case 2: // max EMCAL with TRD 2011
    if (!fUsePhiCut) fUsePhiCut=1;
    fMaxPhiCut=2.45;
    break;
  case 3: // max EMCAL with TRD 2012
    if (!fUsePhiCut) fUsePhiCut=1;
    fMaxPhiCut=2.10;
    break;
  case 4:
    if( !fUsePhiCut ) fUsePhiCut=1;
    fMaxPhiCut = 5.59;//PHOS acceptance
    break;
  case 5:
    if( !fUsePhiCut ) fUsePhiCut=1;
    fMaxPhiCut = 5.5658;//DCal acceptance
    break;
  case 6:
    if( !fUsePhiCut ) fUsePhiCut=1;
    fMaxPhiCut = 5.59;//PHOS acceptance RUN2
    break;
  case 7:
    if( !fUsePhiCut ) fUsePhiCut=1;
    fMaxPhiCut = 3.15;//EMCal acceptance Run2 w/o stripe
    fMaxPhiCutDMC = 5.5658;//DCal acceptance
    break;
  case 8:
    if( !fUsePhiCut ) fUsePhiCut=1;
    fMaxPhiCut = 3.28;//EMCal acceptance Run2 with stripe (w/o stripe 3.15)
    fMaxPhiCutDMC = 5.5658;//DCal acceptance
    break;
  case 9:
    if( !fUsePhiCut ) fUsePhiCut=1;
    fMaxPhiCut = 3.28;//EMCal acceptance Run2 with stripe
    fMaxPhiCutDMC = 5.70;//DCal acceptance with stripe
    break;
  case 10:
    if( !fUsePhiCut ) fUsePhiCut=1;
    fMaxPhiCut = 2.09;//EMCal acceptance 2010 (1.39626 + 40 degrees)
    break;
  default:
    AliError(Form("Max Phi Cut not defined %d",maxPhi));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetDistanceToBadChannelCut(Int_t distanceToBadChannel)
{
  switch(distanceToBadChannel){
  case 0:
    fUseDistanceToBadChannel=0;
    fMinDistanceToBadChannel=0;
    break;
  case 1:
    if(fUseDistanceToBadChannel!=1) fUseDistanceToBadChannel=1;
    fMinDistanceToBadChannel=1;
    break;
  case 2:
    if(fUseDistanceToBadChannel!=1) fUseDistanceToBadChannel=1;
    fMinDistanceToBadChannel=2;
    break;
  case 3:
    if(fUseDistanceToBadChannel!=1) fUseDistanceToBadChannel=1;
    fMinDistanceToBadChannel=3;
    break;
  case 4:
    if(fUseDistanceToBadChannel!=1) fUseDistanceToBadChannel=1;
    fMinDistanceToBadChannel=4;
    break;
  case 5:
    if(fUseDistanceToBadChannel!=2) fUseDistanceToBadChannel=2;
    fMinDistanceToBadChannel=1;
    break;
  case 6:
    if(fUseDistanceToBadChannel!=2) fUseDistanceToBadChannel=2;
    fMinDistanceToBadChannel=2;
    break;
  case 7:
    if(fUseDistanceToBadChannel!=2) fUseDistanceToBadChannel=2;
    fMinDistanceToBadChannel=3;
    break;
  case 8:
    if(fUseDistanceToBadChannel!=2) fUseDistanceToBadChannel=2;
    fMinDistanceToBadChannel=4;
    break;
  default:
    AliError(Form("minimum distance to bad channel Cut not defined %d",distanceToBadChannel));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetTimingCut(Int_t timing)
{
  switch(timing){
  case 0:
    fUseTimeDiff=0;
    fMinTimeDiff=-500;
    fMaxTimeDiff=500;
    break;
  case 1:
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-10e-7;
    fMaxTimeDiff=10e-7;//1000ns
    break;
  case 2:
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-50e-8;
    fMaxTimeDiff=50e-8;//500ns
    break;
  case 3:
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-20e-8;
    fMaxTimeDiff=20e-8;//200ns
    break;
  case 4:
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-10e-8;
    fMaxTimeDiff=10e-8;//100ns
    break;
  case 5:
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-50e-9;
    fMaxTimeDiff=50e-9;//50ns
    break;
  case 6:
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=35e-9;
    break;
  case 7:
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=30e-9;//30ns
    break;
  case 8:
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-20e-9;
    fMaxTimeDiff=30e-9;
    break;
  case 9:
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-20e-9;
    fMaxTimeDiff=25e-9;
    break;
  case 10: //a
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-12.5e-9;
    fMaxTimeDiff=13e-9;
    break;
  case 11: //b
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-130e-9;
    fMaxTimeDiff=130e-9;
    break;
  case 12: //c
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-110e-9;
    fMaxTimeDiff=110e-9;
    break;
  case 13: //d
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-120e-9;
    fMaxTimeDiff=120e-9;
    break;
  case 14: //e
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-90e-9;
    fMaxTimeDiff=90e-9;
    break;
  case 15: //f
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-80e-9;
    fMaxTimeDiff=80e-9;
    break;
  case 16: //g PHOS timing cut, applying timing cut efficiency in MC
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=30e-9;//30ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 4.0;
    fTimingEfficiencyMCSimClusterHighPtStart = 6.0;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "1 /([0]/([1]*(1./(1.+[2]*exp(-x/[3]))* 1./(1.+[4]*exp((x-[5])/[6])))))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(1.51165e+00,6.41558e-02,1.24776e+01,1.32035e-01,-1.15887e+00,3.89796e+02,2.02598e+03);
    fFuncTimingEfficiencyMCSimClusterHighPt = new TF1("FuncTimingEfficiencyMCSimClusterHighPt", "[0]*x+[1]");
    fFuncTimingEfficiencyMCSimClusterHighPt->SetParameters(-0.008333,1.05);
    break;
  case 17: //h PHOS timing cut, applying timing cut efficiency in MC
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-50e-9;
    fMaxTimeDiff=50e-9;//50ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 4.0;
    fTimingEfficiencyMCSimClusterHighPtStart = 6.0;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "1 /([0]/([1]*(1./(1.+[2]*exp(-x/[3]))* 1./(1.+[4]*exp((x-[5])/[6])))))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(8.36250e-01,1.00398e-01,1.43170e+01,1.04184e-01,-1.24269e+00,3.30702e+02,9.49252e+02);
    fFuncTimingEfficiencyMCSimClusterHighPt = new TF1("FuncTimingEfficiencyMCSimClusterHighPt", "[0]*x+[1]");
    fFuncTimingEfficiencyMCSimClusterHighPt->SetParameters(-0.008333,1.05);
    break;
  case 18: //i PHOS timing cut, applying timing cut efficiency in MC
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-12.5e-9;//-12.5ns
    fMaxTimeDiff=13e-9;//13ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 4.0;
    fTimingEfficiencyMCSimClusterHighPtStart = 6.0;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "1 /([0]/([1]*(1./(1.+[2]*exp(-x/[3]))* 1./(1.+[4]*exp((x-[5])/[6])))))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(1.29282e+00,6.50756e-02,9.57716e+00,2.44441e-01,-1.29253e+00,3.00901e+02,9.62463e+02);
    fFuncTimingEfficiencyMCSimClusterHighPt = new TF1("FuncTimingEfficiencyMCSimClusterHighPt", "[0]*x+[1]");
    fFuncTimingEfficiencyMCSimClusterHighPt->SetParameters(-0.008333,1.05);
    break;
  case 19: //j EMCal timing cut, applying timing cut efficiency in MC
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=30e-9;//12.5ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 4.0;
    fTimingEfficiencyMCSimClusterHighPtStart = 6.0;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "1 /([0]/([1]*(1./(1.+[2]*exp(-x/[3]))* 1./(1.+[4]*exp((x-[5])/[6])))))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(9.92263e-01,1.00921e+00,1.37239e-02,6.56658e-01,2.75953e-01,2.24610e+02,7.98508e+01);
    break;
  case 20: //k EMCal timing cut, applying timing cut efficiency in MC
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-50e-9;
    fMaxTimeDiff=50e-9;//50ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 4.0;
    fTimingEfficiencyMCSimClusterHighPtStart = 6.0;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "1 /([0]/([1]*(1./(1.+[2]*exp(-x/[3]))* 1./(1.+[4]*exp((x-[5])/[6])))))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(9.91952e-01,1.01113e+00,1.08929e-02,6.86663e-01,1.75156e-01,2.19392e+02,9.84879e+01);
    break;
  case 21: //l EMCal timing cut, applying timing cut efficiency in MC
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-12.5e-9;//-12.5ns
    fMaxTimeDiff=13e-9;//13ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 4.0;
    fTimingEfficiencyMCSimClusterHighPtStart = 6.0;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "1 /([0]/([1]*(1./(1.+[2]*exp(-x/[3]))* 1./(1.+[4]*exp((x-[5])/[6])))))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(1.61742e+00,1.69595e+00,1.07106e-01,4.83419e-01,1.42459e-01,1.97986e+02,1.82539e+02);
    break;
  case 22: //m PHOS timing cut, 13TeV MB 30ns, applying timing cut efficiency in MC
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=30e-9;//30ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 4.0;
    fTimingEfficiencyMCSimClusterHighPtStart = 6.0;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "1 /([0]/([1]*(1./(1.+[2]*exp(-x/[3]))* 1./(1.+[4]*exp((x-[5])/[6])))))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(1.01021e+00,1.00143e+00,1.36545e+01,1.49372e-01,-1.09826e-01,5.56485e+02,1.25420e+01);
    fFuncTimingEfficiencyMCSimClusterHighPt = new TF1("FuncTimingEfficiencyMCSimClusterHighPt", "[0]*(x-[1])");
    fFuncTimingEfficiencyMCSimClusterHighPt->SetParameters(-3.89726e-02, 6.00000e+00);
    break;
  case 23: //n PHOS timing cut, 13TeV Trigger 30ns, applying timing cut efficiency in MC
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=30e-9;//30ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 4.0;
    fTimingEfficiencyMCSimClusterHighPtStart = 6.0;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "1 /([0]/([1]*(1./(1.+[2]*exp(-x/[3]))* 1./(1.+[4]*exp((x-[5])/[6])))))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(1.01021e+00,1.00143e+00,1.36545e+01,1.49372e-01,-1.09826e-01,5.56485e+02,1.25420e+01);
    fFuncTimingEfficiencyMCSimClusterHighPt = new TF1("FuncTimingEfficiencyMCSimClusterHighPt", "(x<[3])*(((1.-[2])*exp(-([1]*(x-[0]))))+[2])+(x>[3])*((((1.-[2])*exp(-([1]*([3]-[0]))))+[2])+((x-[3])*[4]))");
    fFuncTimingEfficiencyMCSimClusterHighPt->SetParameters(6.00000e+00, 3.50809e-01, 6.96152e-01, 1.54686e+01, 2.55793e-03);
    break;
  case 24: //o PHOS timing cut, 13TeV MB 30ns, applying timing cut efficiency in MC
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=30e-9;//30ns
    fMinTimeDiffHighPt=-150e-9;
    fMaxTimeDiffHighPt=150e-9;//150ns
    fUseTimingEfficiencyMCSimCluster = 2;
    fTimingEfficiencyMCSimClusterLowPtEnd = 4.0;
    fTimingEfficiencyMCSimClusterHighPtStart = 6.0;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "1 /([0]/([1]*(1./(1.+[2]*exp(-x/[3]))* 1./(1.+[4]*exp((x-[5])/[6])))))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(1.01021e+00,1.00143e+00,1.36545e+01,1.49372e-01,-1.09826e-01,5.56485e+02,1.25420e+01);
    break;
  case 25: //p PHOS timing cut, 13TeV Trigger 30ns, applying timing cut efficiency in MC
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=30e-9;//30ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 4.0;
    fTimingEfficiencyMCSimClusterHighPtStart = 6.0;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "exp(([0]+[1]*x-[2]*x*x+x*x*x)/(1.-[3]*x+[4]*x*x+x*x*x))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(-8.26274e+01,6.47055e+01,2.18057e+01,1.36203e+02,1.06448e+03);
    fFuncTimingEfficiencyMCSimClusterHighPt = new TF1("FuncTimingEfficiencyMCSimClusterHighPt", "(x<[3])*(((1.-[2])*exp(-([1]*(x-[0]))))+[2])+(x>[3])*((((1.-[2])*exp(-([1]*([3]-[0]))))+[2])+((x-[3])*[4]))");
    fFuncTimingEfficiencyMCSimClusterHighPt->SetParameters(6.00000e+00, 3.24850e-01, 6.75747e-01, 1.58392e+01, 4.18167e-03);
    break;
  case 26: //q PHOS timing cut, 13TeV Trigger 30ns, applying timing cut efficiency in MC
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=30e-9;//30ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 4.0;
    fTimingEfficiencyMCSimClusterHighPtStart = 6.0;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "exp(([0]+[1]*x-[2]*x*x+x*x*x)/(1.-[3]*x+[4]*x*x+x*x*x))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(-8.44603e+01,5.31853e+01,1.77562e+01,1.07382e+02,9.87056e+02);
    fFuncTimingEfficiencyMCSimClusterHighPt = new TF1("FuncTimingEfficiencyMCSimClusterHighPt", "(x<[3])*(((1.-[2])*exp(-([1]*(x-[0]))))+[2])+(x>[3])*((((1.-[2])*exp(-([1]*([3]-[0]))))+[2])+((x-[3])*[4]))");
    fFuncTimingEfficiencyMCSimClusterHighPt->SetParameters(6.00000e+00, 3.50809e-01, 6.96152e-01, 1.54686e+01, 2.55793e-03);
    break;
  case 27: //r PHOS timing cut, 13TeV Trigger 30ns, by Signal Extraction; 2GeV<ETag<5.5GeV, |TimingTag|<30ns, |TimingProbe|<1000ns, LowPt from MB; HighPt from Trigger
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=30e-9;//30ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 4.0;
    fTimingEfficiencyMCSimClusterHighPtStart = 5.60995e+00;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "exp(([0]+[1]*x-[2]*x*x+x*x*x)/(1.-[3]*x+[4]*x*x+x*x*x))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(-9.60861e+00,  1.37754e+01, 7.23636e+00, -9.95086e+00, 1.07023e+02);
    fFuncTimingEfficiencyMCSimClusterHighPt = new TF1("FuncTimingEfficiencyMCSimClusterHighPt", "(x<[3])*(((1.-[2])*exp(-([1]*(x-[0]))))+[2])+(x>[3])*((((1.-[2])*exp(-([1]*([3]-[0]))))+[2])+((x-[3])*[4]))");
    fFuncTimingEfficiencyMCSimClusterHighPt->SetParameters(5.60995e+00, 1.74831e-01, 9.28195e-01, 8.92803e+00,  0.00000e+00);
    break;
  case 28: //s PHOS timing cut, 13TeV Trigger 30ns, by Signal Extraction; 2GeV<ETag<5.5GeV, |TimingTag|<30ns, |TimingProbe|<1000ns, LowPt from Trigger; HighPt from Trigger
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=30e-9;//30ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 4.0;
    fTimingEfficiencyMCSimClusterHighPtStart = 5.60995e+00;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "exp(([0]+[1]*x-[2]*x*x+x*x*x)/(1.-[3]*x+[4]*x*x+x*x*x))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(-1.50291e+01,  1.78481e+01, 8.18420e+00, -3.18997e+01, 1.66224e+02);
    fFuncTimingEfficiencyMCSimClusterHighPt = new TF1("FuncTimingEfficiencyMCSimClusterHighPt", "(x<[3])*(((1.-[2])*exp(-([1]*(x-[0]))))+[2])+(x>[3])*((((1.-[2])*exp(-([1]*([3]-[0]))))+[2])+((x-[3])*[4]))");
    fFuncTimingEfficiencyMCSimClusterHighPt->SetParameters(5.60995e+00, 1.74831e-01, 9.28195e-01, 8.92803e+00,  0.00000e+00);
    break;
  case 29: //t PHOS timing cut, 13TeV Trigger 25ns, by Signal Extraction; 2GeV<ETag<5.5GeV, |TimingTag|<25ns, |TimingProbe|<1000ns, LowPt from Trigger; HighPt from Trigger
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-25e-9;
    fMaxTimeDiff=25e-9;//30ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 4.0;
    fTimingEfficiencyMCSimClusterHighPtStart = 6.0;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "exp(([0]+[1]*x-[2]*x*x+x*x*x)/(1.-[3]*x+[4]*x*x+x*x*x))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(-1.92595e+01,  2.16880e+01, 9.09869e+00, -5.34259e+01, 9.34326e+01);
    fFuncTimingEfficiencyMCSimClusterHighPt = new TF1("FuncTimingEfficiencyMCSimClusterHighPt", "(x<[3])*(((1.-[2])*exp(-([1]*(x-[0]))))+[2])+(x>[3])*((((1.-[2])*exp(-([1]*([3]-[0]))))+[2])+((x-[3])*[4]))");
    fFuncTimingEfficiencyMCSimClusterHighPt->SetParameters(5.50000e+00, 1.75640e-01, 8.83017e-01, 9.15908e+00, 0.00000e+00);
    break;
  case 30: //u PHOS timing cut, 13TeV Trigger 50ns, by Signal Extraction; 2GeV<ETag<5.5GeV, |TimingTag|<50ns, |TimingProbe|<1000ns, LowPt from Trigger; HighPt from Trigger
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-50e-9;
    fMaxTimeDiff=50e-9;//30ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 5.5;
    fTimingEfficiencyMCSimClusterHighPtStart = 100e+00;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "exp(([0]+[1]*x-[2]*x*x+x*x*x)/(1.-[3]*x+[4]*x*x+x*x*x))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(-6.60829e+00, 1.08812e+01, 7.38930e+00, 9.63108e+01, 5.55664e+02);
    fFuncTimingEfficiencyMCSimClusterHighPt = new TF1("FuncTimingEfficiencyMCSimClusterHighPt", "(x<[3])*(((1.-[2])*exp(-([1]*(x-[0]))))+[2])+(x>[3])*((((1.-[2])*exp(-([1]*([3]-[0]))))+[2])+((x-[3])*[4]))");
    fFuncTimingEfficiencyMCSimClusterHighPt->SetParameters(5.50000e+00, 2.62837e-01, 9.84211e-01,  9.25566e+00, 0.00000e+00);
    break;
  case 31: //v PHOS timing cut, 13TeV Trigger 30ns by Signal Extraction, applying timing cut efficiency in MC; 2GeV<ETag<5.5GeV, |TimingTag|<30ns, |TimingProbe|<1000ns, LowPt from Tr; HighPt constant
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=30e-9;//30ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 5.5;
    fTimingEfficiencyMCSimClusterHighPtStart = 100e+00;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "exp(([0]+[1]*x-[2]*x*x+x*x*x)/(1.-[3]*x+[4]*x*x+x*x*x))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(-2.36670e+01, 2.27558e+01, 9.74890e+00,  2.59148e+01, 4.35144e+02);
    fFuncTimingEfficiencyMCSimClusterHighPt = new TF1("FuncTimingEfficiencyMCSimClusterHighPt", "[0]");
    fFuncTimingEfficiencyMCSimClusterHighPt->SetParameter(0, 1.);
    break;
  case 32: //w PHOS timing cut, 13TeV Trigger 30ns by Signal Extraction, applying timing cut efficiency in MC; 2GeV<ETag<5.5GeV, |TimingTag|<30ns, |TimingProbe|<1000ns, LowPt from MB; HighPt from Trigger
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=30e-9;//30ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 5.5;
    fTimingEfficiencyMCSimClusterHighPtStart = 5.71917e+00;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "exp(([0]+[1]*x-[2]*x*x+x*x*x)/(1.-[3]*x+[4]*x*x+x*x*x))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(-2.81491e+01, 2.71243e+01, 1.14455e+01, 1.12762e+02, 7.06898e+02);
    fFuncTimingEfficiencyMCSimClusterHighPt = new TF1("FuncTimingEfficiencyMCSimClusterHighPt", "(x<[3])*(((1.-[2])*exp(-([1]*(x-[0]))))+[2])+(x>[3])*((((1.-[2])*exp(-([1]*([3]-[0]))))+[2])+((x-[3])*[4]))");
    fFuncTimingEfficiencyMCSimClusterHighPt->SetParameters(5.71917e+00, 1.87419e-01, 5.80457e-01, 1.28769e+01, 5.37268e-03);
    break;
  case 33: //x PHOS timing cut, 13TeV Trigger 30ns by Signal Extraction, applying timing cut efficiency in MC; 2GeV<ETag<5.5GeV, |TimingTag|<30ns, |TimingProbe|<1000ns, LowPt from Tr; HighPt from Trigger
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=30e-9;//30ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 5.5;
    fTimingEfficiencyMCSimClusterHighPtStart = 5.71917e+00;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "exp(([0]+[1]*x-[2]*x*x+x*x*x)/(1.-[3]*x+[4]*x*x+x*x*x))");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(-2.36670e+01, 2.27558e+01, 9.74890e+00,  2.59148e+01, 4.35144e+02);
    fFuncTimingEfficiencyMCSimClusterHighPt = new TF1("FuncTimingEfficiencyMCSimClusterHighPt", "(x<[3])*(((1.-[2])*exp(-([1]*(x-[0]))))+[2])+(x>[3])*((((1.-[2])*exp(-([1]*([3]-[0]))))+[2])+((x-[3])*[4]))");
    fFuncTimingEfficiencyMCSimClusterHighPt->SetParameters(5.71917e+00, 1.87419e-01, 5.80457e-01, 1.28769e+01, 5.37268e-03);
    break;
  case 34: //y PHOS timing cut, pPb8TeV Trigger 30ns from Dmitri
    if (!fUseTimeDiff) fUseTimeDiff=1;
    fMinTimeDiff=-30e-9;
    fMaxTimeDiff=30e-9;//30ns
    fUseTimingEfficiencyMCSimCluster = 1;
    fTimingEfficiencyMCSimClusterLowPtEnd = 1e5;
    fTimingEfficiencyMCSimClusterHighPtStart = 1e6;
    fFuncTimingEfficiencyMCSimCluster = new TF1("FuncTimingEfficiencyMCSimCluster", "(x<2.5)*exp(([0]+[1]*x-[2]*x*x+x*x*x)/(1.-[3]*x+[4]*x*x+x*x*x))+(x>=2.5)*[5]");
    fFuncTimingEfficiencyMCSimCluster->SetParameters(-7.35340e+01, 7.14029e+01, 2.25335e+01,  4.99060e+01, 1.28905e+03, 0.9975);
    fFuncTimingEfficiencyMCSimClusterHighPt = new TF1("FuncTimingEfficiencyMCSimClusterHighPt", "[0]");
    fFuncTimingEfficiencyMCSimClusterHighPt->SetParameter(0, 1.0);
    break;
  default:
    AliError(Form("Timing Cut not defined %d",timing));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetTrackMatchingCut(Int_t trackMatching)
{
  // matching parameters for EMCal clusters
  if(fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
    switch(trackMatching){
    case 0:
      fUseDistTrackToCluster = kFALSE;
      fMaxDistTrackToClusterEta = 0;
      fMinDistTrackToClusterPhi = 0;
      fMaxDistTrackToClusterPhi = 0;
      break;
    case 1:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fMaxDistTrackToClusterEta = 0.008;//0.015;
      fMinDistTrackToClusterPhi = -0.03;//-0.01;
      fMaxDistTrackToClusterPhi = 0.03;//0.03;//0.04;
      break;
    case 2:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fMaxDistTrackToClusterEta = 0.012;//0.015;
      fMinDistTrackToClusterPhi = -0.05;//-0.01;
      fMaxDistTrackToClusterPhi = 0.04;//0.035;//0.05;
      break;
    case 3:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fMaxDistTrackToClusterEta = 0.016;//0.015;
      fMinDistTrackToClusterPhi = -0.09;//-0.015;
      fMaxDistTrackToClusterPhi = 0.06;//0.04;//0.1;
      break;
    case 4:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fMaxDistTrackToClusterEta = 0.018;//0.015;
      fMinDistTrackToClusterPhi = -0.11;//-0.015;
      fMaxDistTrackToClusterPhi = 0.07;//0.045;//0.13;
      break;
    case 5:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fMaxDistTrackToClusterEta = 0.02;//0.015;
      fMinDistTrackToClusterPhi = -0.13;//-0.02;
      fMaxDistTrackToClusterPhi = 0.08;//0.05;//0.15
      break;
    // pT dependent matching parameters
    case 6:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fUsePtDepTrackToCluster = 1;
      fFuncPtDepEta = new TF1("funcEta6", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.03, 0.010, 2.5);
      fFuncPtDepPhi = new TF1("funcPhi6", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.08, 0.015, 2.);
      break;
    case 7:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fUsePtDepTrackToCluster = 1;
      fFuncPtDepEta = new TF1("funcEta7", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.04, 0.010, 2.5);
      fFuncPtDepPhi = new TF1("funcPhi7", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.09, 0.015, 2.);
      break;
    case 8:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fUsePtDepTrackToCluster = 1;
      fFuncPtDepEta = new TF1("funcEta8", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.05, 0.010, 2.5);
      fFuncPtDepPhi = new TF1("funcPhi8", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.10, 0.015, 1.75);
      break;
    case 9:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fUsePtDepTrackToCluster = 1;
      fFuncPtDepEta = new TF1("funcEta9", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.06, 0.015, 2.5);
      fFuncPtDepPhi = new TF1("funcPhi9", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.12, 0.020, 1.75);
      break;
    case 10:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fUsePtDepTrackToCluster = 1;
      fFuncPtDepEta = new TF1("funcEta10", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.035, 0.010, 2.5);
      fFuncPtDepPhi = new TF1("funcPhi10", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.085, 0.015, 2.0);
      break;
    case 11:
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fUsePtDepTrackToCluster = 1;
      fFuncPtDepEta = new TF1("funcEta11", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.045, 0.010, 2.5);
      fFuncPtDepPhi = new TF1("funcPhi11", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.095, 0.015, 1.75);
      break;
    case 12: // here starts clusterE/trackP veto for TM, taking case 7 as usual TM cuts; cut char 'c'
      // note: this case does not apply the cut (fEOverPMax = 9e9), but sets fUseEOverPVetoTM = kTRUE, so that reference histos are still produced
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      if (!fUseEOverPVetoTM) fUseEOverPVetoTM=kTRUE;
      fUsePtDepTrackToCluster = 1;
      fFuncPtDepEta = new TF1("funcEta12", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.04, 0.010, 2.5);
      fFuncPtDepPhi = new TF1("funcPhi12", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.09, 0.015, 2.);

      fEOverPMax = 9e9;
      break;
    case 13: // loosest E/P cut; cut char 'd'
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      if (!fUseEOverPVetoTM) fUseEOverPVetoTM=kTRUE;
      fUsePtDepTrackToCluster = 1;
      fFuncPtDepEta = new TF1("funcEta13", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.04, 0.010, 2.5);
      fFuncPtDepPhi = new TF1("funcPhi13", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.09, 0.015, 2.);

      fEOverPMax = 3.0;
      break;
    case 14: // cut char 'e'
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      if (!fUseEOverPVetoTM) fUseEOverPVetoTM=kTRUE;
      fUsePtDepTrackToCluster = 1;
      fFuncPtDepEta = new TF1("funcEta14", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.04, 0.010, 2.5);
      fFuncPtDepPhi = new TF1("funcPhi14", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.09, 0.015, 2.);

      fEOverPMax = 2.0;
      break;
    case 15: // cut char 'f'
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      if (!fUseEOverPVetoTM) fUseEOverPVetoTM=kTRUE;
      fUsePtDepTrackToCluster = 1;
      fFuncPtDepEta = new TF1("funcEta15", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.04, 0.010, 2.5);
      fFuncPtDepPhi = new TF1("funcPhi15", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.09, 0.015, 2.);

      fEOverPMax = 1.75;
      break;
    case 16: // cut char 'g'
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      if (!fUseEOverPVetoTM) fUseEOverPVetoTM=kTRUE;
      fUsePtDepTrackToCluster = 1;
      fFuncPtDepEta = new TF1("funcEta16", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.04, 0.010, 2.5);
      fFuncPtDepPhi = new TF1("funcPhi16", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.09, 0.015, 2.);

      fEOverPMax = 1.5;
      break;

    case 17: // hardest E/p cut; cut char 'h'
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      if (!fUseEOverPVetoTM) fUseEOverPVetoTM=kTRUE;
      fUsePtDepTrackToCluster = 1;
      fFuncPtDepEta = new TF1("funcEta17", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.04, 0.010, 2.5);
      fFuncPtDepPhi = new TF1("funcPhi17", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.09, 0.015, 2.);

      fEOverPMax = 1.25;
      break;
    case 18: //i TM cut for PbPb EMC clusters - tight
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      if (!fUseTMMIPsubtraction) fUseTMMIPsubtraction=kTRUE;
      fMaxDistTrackToClusterEta = 0.010;
      fMinDistTrackToClusterPhi = -0.011;
      fMaxDistTrackToClusterPhi = 0.011;
      break;
    case 19: //j TM cut for PbPb EMC clusters - loose
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      if (!fUseTMMIPsubtraction) fUseTMMIPsubtraction=kTRUE;
      fMaxDistTrackToClusterEta = 0.015;
      fMinDistTrackToClusterPhi = -0.02;
      fMaxDistTrackToClusterPhi = 0.02;
      break;
    case 20: //k TM cut for PbPb EMC clusters
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      if (!fUseTMMIPsubtraction) fUseTMMIPsubtraction=kTRUE;
      fUsePtDepTrackToCluster = 1;
      fFuncPtDepEta = new TF1("funcEta20", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.04, 0.010, 2.5);
      fFuncPtDepPhi = new TF1("funcPhi20", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.09, 0.015, 2.);
      break;
    case 21: //l Standard TM cut + secondary trackmatching
      if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
      fUsePtDepTrackToCluster = 1;
      fFuncPtDepEta = new TF1("func", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepEta->SetParameters(0.04, 0.010, 2.5);
      fFuncPtDepPhi = new TF1("func", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
      fFuncPtDepPhi->SetParameters(0.09, 0.015, 2.);
      fDoSecondaryTrackMatching = kTRUE;
      break;

    default:
      AliError(Form("Track Matching Cut not defined %d",trackMatching));
      return kFALSE;
    }
    return kTRUE;
  // matching parameters for PHOS clusters
  }else if(fClusterType == 2) {
    switch(trackMatching){
      case 0:
        fUseDistTrackToCluster = kFALSE;
        fMaxDistTrackToClusterEta = 0;
        fMinDistTrackToClusterPhi = 0;
        fMaxDistTrackToClusterPhi = 0;
        break;
      case 1:
        if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
        fMaxDistTrackToClusterEta = 0.005;//0.015;
        fMinDistTrackToClusterPhi = -0.03;//-0.025;
        fMaxDistTrackToClusterPhi = 0.03;//0.06;//0.3;
        break;
      case 2:
        if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
        fMaxDistTrackToClusterEta = 0.01;//0.015;
        fMinDistTrackToClusterPhi = -0.09;//-0.025;
        fMaxDistTrackToClusterPhi = 0.07;//0.07;//0.4;
        break;
      case 3:
        if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
        fMaxDistTrackToClusterEta = 0.015;//0.02;
        fMinDistTrackToClusterPhi = -0.15;//-0.03;
        fMaxDistTrackToClusterPhi = 0.11;//0.1;//0.5;
        break;
      case 4: //pT dependent for PCM-PHOS "default" selection
        if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
        fUsePtDepTrackToCluster = 1;
        fFuncPtDepEta = new TF1("funcEta4PHOS", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
        fFuncPtDepEta->SetParameters(0.05, 0.005, 3.0);

        fFuncPtDepPhi = new TF1("funcPhi4PHOS", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
        fFuncPtDepPhi->SetParameters(0.33, 0.005, 2.3);
        break;
      case 5: //pT dependent for PCM-PHOS tight selection
        if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
        fUsePtDepTrackToCluster = 1;
        fFuncPtDepEta = new TF1("funcEta5PHOS", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
        fFuncPtDepEta->SetParameters(0.025, 0.002, 3.0);

        fFuncPtDepPhi = new TF1("funcPhi5PHOS", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
        fFuncPtDepPhi->SetParameters(0.17, 0.005, 2.5);
        break;
      case 6: //pT dependent for PCM-PHOS loose selection
        if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
        fUsePtDepTrackToCluster = 1;
        fFuncPtDepEta = new TF1("funcEta6PHOS", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
        fFuncPtDepEta->SetParameters(0.07, 0.003, 2.5);

        fFuncPtDepPhi = new TF1("funcPhi6PHOS", "[1] + 1 / pow(x + pow(1 / ([0] - [1]), 1 / [2]), [2])");
        fFuncPtDepPhi->SetParameters(0.45, 0.010, 2.0);
        break;
      case 7: //
        if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
        fUsePtDepTrackToCluster = 2;
        fMinTMDistSigma         = 2.;
        break;
      case 8: //
        if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
        fUsePtDepTrackToCluster = 2;
        fMinTMDistSigma         = 2.5;
        break;
      case 9: //
        if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
        fUsePtDepTrackToCluster = 2;
        fMinTMDistSigma         = 3;
        break;
      case 10: // a loose phi
        if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
        fMaxDistTrackToClusterEta = 0.02;
        fMinDistTrackToClusterPhi = -0.08;
        fMaxDistTrackToClusterPhi = 0.08;
      case 11: // b strict phi
        if (!fUseDistTrackToCluster) fUseDistTrackToCluster=kTRUE;
        fMaxDistTrackToClusterEta = 0.02;
        fMinDistTrackToClusterPhi = -0.03;
        fMaxDistTrackToClusterPhi = 0.03;
        break;
      default:
        AliError(Form("Track Matching Cut not defined %d",trackMatching));
        return kFALSE;
    }
    return kTRUE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetExoticClusterCut(Int_t exoticCell)
{
  switch(exoticCell){
  case 0:
    fUseExoticCluster         = 0;
    fExoticEnergyFracCluster  = 0;
    break;
  case 1:
    if (fUseExoticCluster != 1)
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.995;
    break;
  case 2:
    if (fUseExoticCluster != 1)
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.99;
    break;
  case 3:
    if (fUseExoticCluster != 1)
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.98;
    break;
  case 4:
    if (fUseExoticCluster != 1)
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.975;
    break;
  case 5:
    if (fUseExoticCluster != 1)
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.97;
    break;
  case 6:
    if (fUseExoticCluster != 1)
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.965;
    break;
  case 7:
    if (fUseExoticCluster != 1)
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.96;
    break;
  case 8:
    if (fUseExoticCluster != 1)
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.95;
    break;
  case 9:
    if (fUseExoticCluster != 1)
      fUseExoticCluster       = 1;
    fExoticEnergyFracCluster  = 0.94;
    break;

    // exotics rejection with cluster only in one TCard above threshold
  case 10: //a
    if (fUseExoticCluster != 2)
      fUseExoticCluster       = 2;
    fExoticEnergyFracCluster  = 0.95;
    fExoticMinEnergyTCard     = 40;
    break;
  case 11: //b
    if (fUseExoticCluster != 2)
      fUseExoticCluster       = 2;
    fExoticEnergyFracCluster  = 0.95;
    fExoticMinEnergyTCard     = 50;
    break;
  case 12: //c
    if (fUseExoticCluster != 2)
      fUseExoticCluster       = 2;
    fExoticEnergyFracCluster  = 0.95;
    fExoticMinEnergyTCard     = 60;
    break;
  case 13: //d
    if (fUseExoticCluster != 2)
      fUseExoticCluster       = 2;
    fExoticEnergyFracCluster  = 0.97;
    fExoticMinEnergyTCard     = 40;
    break;
  case 14: //e
    if (fUseExoticCluster != 2)
      fUseExoticCluster       = 2;
    fExoticEnergyFracCluster  = 0.97;
    fExoticMinEnergyTCard     = 50;
    break;
  case 15: //f
    if (fUseExoticCluster != 2)
      fUseExoticCluster       = 2;
    fExoticEnergyFracCluster  = 0.97;
    fExoticMinEnergyTCard     = 60;
    break;
  case 16: //g
    if (fUseExoticCluster != 3)
      fUseExoticCluster       = 3;
    fExoticEnergyFracCluster  = 0;
    break;
  default:
    AliError(Form("Exotic cell Cut not defined %d",exoticCell));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinEnergyCut(Int_t minEnergy)
{
  if (fIsPureCalo != 1){
    if (fClusterType!=2) {
      switch(minEnergy){
        case 0:
          if (!fUseMinEnergy) fUseMinEnergy=0;
          fMinEnergy=0.1;
          break;
        case 1:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.5;
          break;
        case 2:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.6;
          break;
        case 3:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.7;
          break;
        case 4:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.8;
          break;
        case 5:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.9;
          break;
        case 6:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=4.5;
          break;
        case 7:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=5.0;
          break;
        case 8:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=5.5;
          break;
        case 9:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=6.0;
          break;
        case 10: // a
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=1.5;
          break;
        case 11: // b
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=1.0;
          break;
        case 12: // c
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.65;
          break;
        case 13: // d
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.675;
          break;
        case 14: // e
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.625;
          break;
        case 15: // f
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.7;
          fDoFlatEnergySubtraction=kTRUE;
          break;
        case 16: // g
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.3;
          break;
        case 17: // h
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.4;
          break;
        case 18: // i
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.2;
          break;
        default:
          AliError(Form("Minimum Energy Cut not defined %d",minEnergy));
          return kFALSE;
      }
    } else {
      switch(minEnergy){
        case 0:
          if (!fUseMinEnergy) fUseMinEnergy=0;
          fMinEnergy=0.1;
          break;
        case 1:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.3;
          break;
        case 2:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.5;
          break;
        case 3:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.6;
          break;
        case 4:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.7;
          break;
        case 5:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.8;
          break;
        case 6:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.9;
          break;
        case 7:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.2;
          break;
        case 8:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.4;
          break;
        case 9:
          if (!fUseMinEnergy) fUseMinEnergy=1;
          fMinEnergy=0.1;
          break;
        default:
          AliError(Form("Minimum Energy Cut not defined %d",minEnergy));
          return kFALSE;
      }
    }
    return kTRUE;
  } else {
    switch(minEnergy){
    case 0:
      if (!fUseMinEnergy) fUseMinEnergy=0;
      fMinEnergy=0.1;
      break;
    case 1:
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=4.;
      break;
    case 2:
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=5.;
      break;
    case 3:
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=6.;
      break;
    case 4:
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=7.;
      break;
    case 5:
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=7.5;
      break;
    case 6:
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=8.;
      break;
    case 7:
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=8.5;
      break;
    case 8:
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=9.;
      break;
    case 9:
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=9.5;
      break;
    case 10:
      if (!fUseMinEnergy) fUseMinEnergy=1;
      fMinEnergy=1.5;
      break;
    default:
      AliError(Form("Minimum Energy Cut not defined %d",minEnergy));
      return kFALSE;
    }
    return kTRUE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinNCellsCut(Int_t minNCells)
{
  switch(minNCells){
  case 0:
    if (!fUseNCells) fUseNCells=0;
    fMinNCells=0;
    break;
  case 1:
    if (!fUseNCells) fUseNCells=1;
    fMinNCells=1;
    break;
  case 2:
    if (!fUseNCells) fUseNCells=1;
    fMinNCells=2;
    break;
  case 3:
    if (!fUseNCells) fUseNCells=1;
    fMinNCells=3;
    break;
  case 4:
    if (!fUseNCells) fUseNCells=1;
    fMinNCells=4;
    break;
  case 5:
    if (!fUseNCells) fUseNCells=1;
    fMinNCells=5;
    break;
  case 6:
    if (!fUseNCells) fUseNCells=1;
    fMinNCells=6;
    break;

  // special cases for PHOS: only use the Ncell cut for clusters with a minimal energy
  // if the first number is 1, the chosen cut will only be used if Ecluster > 1 GeV
  case 12: // c
    if (!fUseNCells) fUseNCells=2;
    fMinNCells=2;
    break;
  case 13: // d
    if (!fUseNCells) fUseNCells=2;
    fMinNCells=3;
    break;

  // special cases for EMCal: this will randomly evaluate the NCell cut efficiency for MC
  // and let clusters with NCell<2 pass if sucessful, for data the normal NCell cut is applied
  case 17: // h
    if (!fUseNCells) fUseNCells=3;
    fMinNCells=2;
    fFuncNCellCutEfficiencyEMCal = new TF1("fFuncNCellCutEfficiencyEMCal", "gaus(0)+gaus(3)+gaus(6)");
    fFuncNCellCutEfficiencyEMCal->SetParameters(2.55402e-01, 7.55563e-01, 4.26875e-01, -2.73506e-01, 8.06193e-01, 9.24605e-01, 3.30671e-01, 1.63428e+00, 5.30527e-01);
    break;
  case 18: // i
    if (!fUseNCells) fUseNCells=3;
    fMinNCells=3;
    fFuncNCellCutEfficiencyEMCal = new TF1("fFuncNCellCutEfficiencyEMCal", "gaus(0)+gaus(3)+gaus(6)");
    fFuncNCellCutEfficiencyEMCal->SetParameters(2.55402e-01, 7.55563e-01, 4.26875e-01, -2.73506e-01, 8.06193e-01, 9.24605e-01, 3.30671e-01, 1.63428e+00, 5.30527e-01);
    break;
  case 19: // j
    if (!fUseNCells) fUseNCells=4;
    fMinNCells=2;
    fFuncNCellCutEfficiencyEMCal = new TF1("fFuncNCellCutEfficiencyEMCal", "[0] + TMath::Exp([1]+[2]*x) + [3]*x");
    fFuncNCellCutEfficiencyEMCal->SetParameters(1.06254e+00, 3.81052e+00, -9.93005e+00, -1.95966e-02);
    break;
  case 20: // k
    if (!fUseNCells) fUseNCells=4;
    fMinNCells=2;
    fFuncNCellCutEfficiencyEMCal = new TF1("fFuncNCellCutEfficiencyEMCal", "[0] + TMath::Exp([1]+[2]*x)");
    fFuncNCellCutEfficiencyEMCal->SetParameters(9.98136e-01, -4.25858e-01, -9.53998e-01);
    break;
  default:
    AliError(Form("Min N cells Cut not defined %d",minNCells));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMaxM02(Int_t maxM02)
{
  fMaxM02CutNr = maxM02;
  if (fIsPureCalo == 1){
    fUseM02 = 2;
    return kTRUE;
  }

  switch(maxM02){
    case 0:
      if (!fUseM02) fUseM02=0;
      fMaxM02=100;
      break;
    case 1:
      if (!fUseM02) fUseM02=1;
      fMaxM02=1.;
      break;
    case 2:
      if (!fUseM02) fUseM02=1;
      fMaxM02=0.7;
      break;
    case 3:
      if (!fUseM02) fUseM02=1;
      fMaxM02=0.5;
      break;
    case 4:
      if (!fUseM02) fUseM02=1;
      fMaxM02=0.4;
      break;
    case 5:
      if (!fUseM02) fUseM02=1;
      fMaxM02=0.3;
      break;
    case 6:
      if (!fUseM02) fUseM02=1;
      fMaxM02=0.27;
      break;
    case 7: // PHOS cuts
      if (!fUseM02) fUseM02=1;
      fMaxM02=1.3;
      break;
    case 8: // PHOS cuts
      if (!fUseM02) fUseM02=1;
      fMaxM02=2.5;
      break;
    case 9:
      if (!fUseM02) fUseM02=1;
      fMaxM02=0.35;
      break;
    case 10: // a
      if (!fUseM02) fUseM02=1;
      fMaxM02=0.33;
      break;
    case 11: // b
      if (!fUseM02) fUseM02=1;
      fMaxM02=0.28;
      break;
    case 12: // c
      if (!fUseM02) fUseM02=1;
      fMaxM02=0.32;
      break;

    // E dependent M02 variations
    case 13:  // d
      //(0.27 + 0.0072 * TMath::Power(clusEnergy,2));
      fUseM02=2;
      fMinM02CutNr=9;
      fMaxM02=0.4;
      break;
    case 14:  // e
      //(0.31 + 0.0072 * TMath::Power(clusEnergy,2));
      fUseM02=2;
      fMinM02CutNr=9;
      fMaxM02=0.5;
      break;
    case 15:  // f
      //(0.36 + 0.0072 * TMath::Power(clusEnergy,2));
      fUseM02=2;
      fMinM02CutNr=9;
      fMaxM02=0.7;
      break;
    case 16:  // g
      //(0.37 + 0.0072 * TMath::Power(clusEnergy,2))
      fUseM02=2;
      fMinM02CutNr=9;
      fMaxM02=0.7;
      break;
    case 17:  // h
      //(0.30 + 0.0072 * TMath::Power(clusEnergy,2));
      fUseM02=2;
      fMinM02CutNr=9;
      fMaxM02=0.5;
      break;
    case 18:  // i
      // (0.35 + 0.0072 * TMath::Power(clusEnergy,2));
      fUseM02=2;
      fMinM02CutNr=9;
      fMaxM02=0.7;
      break;
    case 19:  // j
      // (0.25 + 0.0072 * TMath::Power(clusEnergy,2));
      fUseM02=2;
      fMinM02CutNr=9;
      fMaxM02=0.39;
      break;
    case 20:  //k
      //(0.27 + 0.0092 * TMath::Power(clusEnergy,2));
      fUseM02=2;
      fMinM02CutNr=9;
      fMaxM02=0.5;
      break;
    case 21:  // l
      //(0.32 + 0.0072 * TMath::Power(clusEnergy,2));
      fUseM02=2;
      fMinM02CutNr=9;
      fMaxM02=0.5;
      break;
    case 22:  // m
      // (0.32 + 0.0152 * TMath::Power(clusEnergy,2));
      fUseM02=2;
      fMinM02CutNr=9;
      fMaxM02=0.5;
      break;
    case 23:  // n
      // (0.32 + 0.0238 * TMath::Power(clusEnergy,2));
      fUseM02=2;
      fMinM02CutNr=9;
      fMaxM02=0.7;
      break;
    case 24:  // o
      // (0.27 + 0.0092 * TMath::Power(clusEnergy,2));
      fUseM02=2;
      fMinM02CutNr=9;
      fMaxM02=0.7;
      break;
    case 25:  // p
      // (0.32 + 0.0072 * TMath::Power(clusEnergy,2));
      fUseM02=2;
      fMinM02CutNr=9;
      fMaxM02=0.7;
      break;
    case 26:  // q
      // (0.34 + 0.0072 * TMath::Power(clusEnergy,2));
      fUseM02=2;
      fMinM02CutNr=9;
      fMaxM02=0.7;
      break;
    case 27:  // r
      // (0.25 + 0.0072 * TMath::Power(clusEnergy,2));
      fUseM02=2;
      fMinM02CutNr=9;
      fMaxM02=0.5;
      break;
    case 28:  // s
      // (0.32 + 0.0238 * TMath::Power(clusEnergy,2));
      fUseM02=2;
      fMinM02CutNr=9;
      fMaxM02=0.5;
      break;
    default:
      AliError(Form("Max M02 Cut not defined %d",maxM02));
      return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Float_t AliCaloPhotonCuts::CalculateMaxM02 (Int_t maxM02, Float_t clusEnergy){
  switch (maxM02){
    case 0:
      return 10;
    case 1:
      if (fMinNLM == 1 && fMaxNLM == 1 ){
        return FunctionM02(clusEnergy, 0.0662, -0.0201, -0.0955, 1.86e-3, 9.91 );
      } else if (fMinNLM == 2 && fMaxNLM == 2 ){
        return FunctionM02(clusEnergy, 0.353, -0.0264, -0.524, 5.59e-3, 21.9 );
      } else {
        return 10;
      }
    case 2:
      if (fMinNLM == 1 && fMaxNLM == 1 ){
        return FunctionM02(clusEnergy, 0.0662, -0.0201, -0.0, 1.86e-3, 9.91 );
      } else if (fMinNLM == 2 && fMaxNLM == 2 ){
        return FunctionM02(clusEnergy, 0.353, -0.0264, -0.424, 5.59e-3, 21.9 );
      } else {
        return 10;
      }
    case 3:
      if (fMinNLM == 1 && fMaxNLM == 1 ){
        return FunctionM02(clusEnergy, 0.0662, -0.0201, -0.2, 1.86e-3, 9.91 );
      } else if (fMinNLM == 2 && fMaxNLM == 2 ){
        return FunctionM02(clusEnergy, 0.353, -0.0264, -0.624, 5.59e-3, 21.9 );
      } else {
        return 10;
      }

    case 13:  // d
      if( (0.27 + 0.0072 * TMath::Power(clusEnergy,2)) >= 0.4) return 0.4;
      else return (0.27 + 0.0072 * TMath::Power(clusEnergy,2));
    case 14:  // e
      if( (0.31 + 0.0072 * TMath::Power(clusEnergy,2)) >= 0.5) return 0.5;
      else return (0.31 + 0.0072 * TMath::Power(clusEnergy,2));
    case 15:  // f
      if( (0.36 + 0.0072 * TMath::Power(clusEnergy,2)) >= 0.7) return 0.7;
      else return (0.36 + 0.0072 * TMath::Power(clusEnergy,2));
    case 16:  // g
      if( (0.37 + 0.0072 * TMath::Power(clusEnergy,2)) >= 0.7) return 0.7;
      else return (0.37 + 0.0072 * TMath::Power(clusEnergy,2));
    case 17:  // h
      if( (0.30 + 0.0072 * TMath::Power(clusEnergy,2)) >= 0.5) return 0.5;
      else return (0.30 + 0.0072 * TMath::Power(clusEnergy,2));
    case 18:  // i
      if( (0.35 + 0.0072 * TMath::Power(clusEnergy,2)) >= 0.7) return 0.7;
      else return (0.35 + 0.0072 * TMath::Power(clusEnergy,2));
    case 19:  // j
      if( (0.25 + 0.0072 * TMath::Power(clusEnergy,2)) >= 0.39) return 0.39;
      else return (0.25 + 0.0072 * TMath::Power(clusEnergy,2));
    case 20:  // k
      if( (0.27 + 0.0092 * TMath::Power(clusEnergy,2)) >= 0.5) return 0.5;
      else return (0.27 + 0.0092 * TMath::Power(clusEnergy,2));
    case 21:  // l
      if( (0.32 + 0.0072 * TMath::Power(clusEnergy,2)) >= 0.5) return 0.5;
      else return (0.32 + 0.0072 * TMath::Power(clusEnergy,2));
    case 22:  // m
      if( (0.32 + 0.0152 * TMath::Power(clusEnergy,2)) >= 0.5) return 0.5;
      else return (0.32 + 0.0152 * TMath::Power(clusEnergy,2));
    case 23:  // n
      if( (0.32 + 0.0238 * TMath::Power(clusEnergy,2)) >= 0.7) return 0.7;
      else return (0.32 + 0.0238 * TMath::Power(clusEnergy,2));
    case 24:  // o
      if( (0.27 + 0.0092 * TMath::Power(clusEnergy,2)) >= 0.7) return 0.7;
      else return (0.27 + 0.0092 * TMath::Power(clusEnergy,2));
    case 25:  // p - standard
      if( (0.32 + 0.0072 * TMath::Power(clusEnergy,2)) >= 0.7) return 0.7;
      else return (0.32 + 0.0072 * TMath::Power(clusEnergy,2));
    case 26:  // q
      if( (0.34 + 0.0072 * TMath::Power(clusEnergy,2)) >= 0.7) return 0.7;
      else return (0.34 + 0.0072 * TMath::Power(clusEnergy,2));
    case 27:  // r
      if( (0.25 + 0.0072 * TMath::Power(clusEnergy,2)) >= 0.5) return 0.5;
      else return (0.25 + 0.0072 * TMath::Power(clusEnergy,2));
    case 28:  // s
      if( (0.32 + 0.0238 * TMath::Power(clusEnergy,2)) >= 0.5) return 0.5;
      else return (0.32 + 0.0238 * TMath::Power(clusEnergy,2));

    default:
      AliError(Form("Max M02 for merged cluster Cut not defined %d",maxM02));
      return 10;
  }
  return 10;

}

//___________________________________________________________________
Float_t AliCaloPhotonCuts::CalculateMinM02 (Int_t minM02, Float_t clusEnergy){
  switch (minM02){
    case 0:
      return 0.;
    case 1:
      if (FunctionM02(clusEnergy, 2.135, -0.245, 0., 0., 0. ) > 0.3)
        return FunctionM02(clusEnergy, 2.135, -0.245, 0., 0., 0. );
      else
        return 0.3;
    case 2:
      if (FunctionM02(clusEnergy, 2.135, -0.245, 0., 0., 0. ) > 0.27)
        return FunctionM02(clusEnergy, 2.135, -0.245, 0., 0., 0. );
      else
        return 0.27;
    case 3:
      if (FunctionM02(clusEnergy, 2.135, -0.245, 0., 0., 0. ) > 0.25)
        return FunctionM02(clusEnergy, 2.135, -0.245, 0., 0., 0. );
      else
        return 0.25;
    case 4:
      if (FunctionM02(clusEnergy, 2.135, -0.245, 0.1, 0., 0. ) > 0.27)
        return FunctionM02(clusEnergy, 2.135, -0.245, 0.1, 0., 0. );
      else
        return 0.27;
    case 5:
      if (FunctionM02(clusEnergy, 2.135, -0.245, -0.1, 0., 0. ) > 0.27)
        return FunctionM02(clusEnergy, 2.135, -0.245, -0.1, 0., 0. );
      else
        return 0.27;
    case 6:
      return 0.3;
    case 7:
      return 0.27;
    case 8:
      return 0.25;
    case 9:
      return 0.1;
    case 10: //a
      return 0.26;
    case 11: //b
      return 0.28;
    case 12: //c
      return 0.29;
    case 13: //d
      return 0.33;
    case 14: //e
      return 0.36;
    case 15: //f
      return 0.39;
    case 16: //g
      return 0.24;
    case 17: //h
      return 0.23;
    case 18: //i
      return 0.22;
    case 19: //j
      return 0.21;
    case 20: //k
      return 0.20;
    case 21: //l
      return 0.19;
    case 22: //m
      return 0.18;

    default:
      AliError(Form("Min M02 for merged cluster Cut not defined %d",minM02));
      return -1;
  }
  return -1;
}


//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinM02(Int_t minM02)
{
  fMinM02CutNr = minM02;
  if (fIsPureCalo == 1){
    fUseM02 = 2;
    return kTRUE;
  }

  switch(minM02){
  case 0:
    if (!fUseM02) fUseM02=0;
    fMinM02=0;
    break;
  case 1:
    if (!fUseM02) fUseM02=1;
    fMinM02=0.002;
    break;
  case 2:
    if (!fUseM02) fUseM02=1;
    fMinM02=0.1;
    break;
  case 3:
    if (!fUseM02) fUseM02=1;
    fMinM02=0.2;
    break;

  // special PHOS cases: apply cut only if the cluster energy is bigger than 1 GeV
  case 11: // b
    if (!fUseM02) fUseM02=3;
    fMinM02 = 0.002;
    break;
  case 12: // c
    if (!fUseM02) fUseM02=3;
    fMinM02 = 0.1;
    break;
  case 13: // d
    if (!fUseM02) fUseM02=3;
    fMinM02 = 0.2;
    break;

  default:
    AliError(Form("Min M02 not defined %d",minM02));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetMinMaxM20(Int_t minM20)
{
  switch(minM20){
  case 0:
    if (!fUseM20) fUseM20=0;
    fMinM20=0;
    fMaxM20=100;
    break;
  case 1:
    if (!fUseM20) fUseM20=1;
    fMinM20=0.002;
    fMaxM20=100;
    break;
  case 2:
    if (!fUseM20) fUseM20=1;
    fMinM20=0;
    fMaxM20=0.5;
    break;
  case 3:
    if (!fUseM20) fUseM20=1;
    fMinM20=0.002;
    fMaxM20=0.5;
    break;
  default:
    AliError(Form("Min M20 Cut not defined %d",minM20));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetRecConv(Int_t recConv)
{
  switch(recConv){
    case 0:
      if (!fUseRecConv) fUseRecConv=0;
      fMaxMGGRecConv=0;
      break;
    case 1:
      if (!fUseRecConv) fUseRecConv=1;
      fMaxMGGRecConv=0.02;
      break;
    case 2:
      if (!fUseRecConv) fUseRecConv=1;
      fMaxMGGRecConv=0.025;
      break;
    case 3:
      if (!fUseRecConv) fUseRecConv=1;
      fMaxMGGRecConv=0.03;
      break;
    case 4:
      if (!fUseRecConv) fUseRecConv=1;
      fMaxMGGRecConv=0.035;
      break;
    case 5:
      if (!fUseRecConv) fUseRecConv=1;
      fMaxMGGRecConv=0.04;
      break;
    case 6:
      if (!fUseRecConv) fUseRecConv=1;
      fMaxMGGRecConv=0.045;
      break;
    case 7:
      if (!fUseRecConv) fUseRecConv=1;
      fMaxMGGRecConv=0.05;
      break;
    case 8:
      if (!fUseRecConv) fUseRecConv=1;
      fMaxMGGRecConv=0.01;
      break;
    case 9:
      if (!fUseRecConv) fUseRecConv=1;
      fMaxMGGRecConv=0.007;
      break;
    default:
      AliError(Form("Conversion Recovery Cut not defined %d",recConv));
      return kFALSE;
  }
  return kTRUE;
}


//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetDispersion(Int_t dispersion)
{
  switch(dispersion){
  case 0:
    if (!fUseDispersion) fUseDispersion=0;
    fMaxDispersion =100;
    break;
  case 1:
    if (!fUseDispersion) fUseDispersion=1;
    fMaxDispersion=2.;
    break;
  case 2:
    if (!fUseDispersion) fUseDispersion=1;
    fMaxDispersion=2.5*2.5;
    break;
  case 3:
    if (!fUseDispersion) fUseDispersion=1;
    fMaxDispersion=2*2;
    break;
  case 4:
    if (!fUseDispersion) fUseDispersion=1;
    fMaxDispersion=3*3;
    break;
  default:
    AliError(Form("Maximum Dispersion Cut not defined %d",dispersion));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetNLM(Int_t nlm)
{
  switch(nlm){
  case 0:
    if (!fUseNLM) fUseNLM=0;
    fMinNLM =0;
    fMaxNLM =100;
    break;
  case 1:
    if (!fUseNLM) fUseNLM=1;
    fMinNLM =1;
    fMaxNLM =1;
    break;
  case 2:
    if (!fUseNLM) fUseNLM=1;
    fMinNLM =2;
    fMaxNLM =2;
    break;
  case 3:
    if (!fUseNLM) fUseNLM=1;
    fMinNLM =1;
    fMaxNLM =2;
    break;
  default:
    AliError(Form("NLM Cut not defined %d",nlm));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetNonLinearity1(Int_t nl1)
{
  if( nl1 >= 0 && nl1 <=9){
    fNonLinearity1 = nl1;
  }
  else{
    AliError(Form("NonLinearity Correction (part1) not defined %d",nl1));
    return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________
Bool_t AliCaloPhotonCuts::SetNonLinearity2(Int_t nl2)
{
  if( nl2 >= 0 && nl2 <=9){
    fNonLinearity2 = nl2;
    if(nl2 == 0) fUseNonLinearity = kFALSE;
    else if(nl2 > 0) fUseNonLinearity = kTRUE;
    fSwitchNonLinearity = fNonLinearity1*10 + fNonLinearity2;
  }
  else{
    AliError(Form("NonLinearity Correction (part2) not defined %d",nl2));
    return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
void AliCaloPhotonCuts::ApplyNonLinearity(AliVCluster* cluster, Int_t isMC, AliVEvent *event)
{
  if(!fUseNonLinearity) return;

  if (!cluster) {
    AliInfo("Cluster pointer null!");
    return;
  }

  Float_t energy = cluster->E();

  Bool_t isDCal = kFALSE;
  Int_t clusterSMID = -1;
  if(fClusterType == 4 && event){
    Int_t largestCellIDcluster = FindLargestCellInCluster(cluster,event);
    if(largestCellIDcluster>-1){
      Int_t dummycol = -1, dummyrow = -1;
      clusterSMID = GetModuleNumberAndCellPosition(largestCellIDcluster, dummycol, dummyrow);
      if(clusterSMID>11)
        isDCal = kTRUE;
    }
  }

  if( fClusterType == 1 || fClusterType == 3|| fClusterType == 4){
    if (energy < 0.05) {
      // Clusters with less than 50 MeV or negative are not possible
      AliInfo(Form("Too Low Cluster energy!, E = %f < 0.05 GeV",energy));
      return;
    }
  } else {
    if (energy < 0.01) {
      // Clusters with less than 50 MeV or negative are not possible
      AliInfo(Form("Too Low Cluster energy!, E = %f < 0.01 GeV",energy));
      return;
    }
  }

  if(fCurrentMC==kNoMC){
    AliV0ReaderV1* V0Reader = (AliV0ReaderV1*) AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
    if( V0Reader == NULL ){
      AliFatal(Form("No V0Reader called '%s' could be found within AliCaloPhotonCuts::ApplyNonLinearity",fV0ReaderName.Data()));
      return;
    }
    fPeriodName = V0Reader->GetPeriodName();
    fCurrentMC = FindEnumForMCSet(fPeriodName);

    printf("AliCaloPhotonCuts:Period name has been set to %s, period-enum: %o\n",fPeriodName.Data(),fCurrentMC ) ;
  }


  Bool_t fPeriodNameAvailable = kTRUE;
  switch(fSwitchNonLinearity){

    // Standard NonLinearity -
    case 1:
      label_case_01:
      if( fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
        // TB parametrization from Nico on Martin 100MeV points (final version incl. fine tuning)
        if(isMC){
          energy /= FunctionNL_OfficialTB_100MeV_MC_V2(energy);
          // fine tuning based on gaussian fits on PCMEMC in pPb5TeV
          energy /= FunctionNL_kSDM(energy, 0.987912, -2.94105, -0.273207) ;
        } else {
          energy /= FunctionNL_OfficialTB_100MeV_Data_V2(energy);
        }
      } else if ( fClusterType == 2 ){
          // Nonlin from PHOS group only MC part
          if(isMC != 0) {
              if( fCurrentMC==k14j4 ){
                  energy *= FunctionNL_PHOSOnlyMC(energy, 1.008, 0.015, 0.4);
                  // for LHC13bc
              } else if( fCurrentMC==kPPb5T13P2DPMJet || fCurrentMC==kPPb5T13P4DPMJet || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c || fCurrentMC == kPPb5T13P2HIJAdd || fCurrentMC == kPPb5T13P4JJhigh || fCurrentMC == kPPb5T13P4JJlow){
                  energy *= FunctionNL_PHOSOnlyMC(energy, 1.0135, 0.018, 1.9);
              } else if(  // pp 5 TeV 2015
                  fCurrentMC == k16h3  || fCurrentMC == k16h8a || fCurrentMC == k16h8b || fCurrentMC == k16k3a  || fCurrentMC == k16k5a ||  fCurrentMC == k16k5b || fCurrentMC == k17e2 || fCurrentMC == k18j3 ||
                  // PbPb 5 TeV 2015
                  fCurrentMC == k16k3b ||
                  // pPb 5 TeV 2016
                  fCurrentMC == kPPb5T16EPOS || fCurrentMC == kPPb5T16DPMJet ||
                  // pPb 8 TeV 2016
                  fCurrentMC == k17f3a || fCurrentMC == k17f3b || fCurrentMC == k17f4a || fCurrentMC == k17f4b ||
                  // XeXe 5.44 TeV 2017
                  fCurrentMC == kXeXe5T17HIJING
              ){
                  energy *= FunctionNL_PHOSOnlyMC(energy, 1.012, -0.06, 0.7);
              }
          }
      }
      break;

    case 2:
      if( fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
        // TB parametrization from Nico on Martin 50MeV points
        if(isMC){
          energy /= FunctionNL_OfficialTB_50MeV_MC(energy);
        } else {
          energy /= FunctionNL_OfficialTB_50MeV_Data(energy);
        }
      }
      break;

    case 3:
      if( fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
        // TB parametrization from Nico on Martin 150MeV points
        if(isMC){
          energy /= FunctionNL_OfficialTB_150MeV_MC(energy);
        } else {
          energy /= FunctionNL_OfficialTB_150MeV_Data(energy);
        }
      }
      break;

    case 4:
      if( fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
        // TB parametrization from Nico on Martin 300MeV points
        if(isMC){
          energy /= FunctionNL_OfficialTB_300MeV_MC(energy);
        } else {
          energy /= FunctionNL_OfficialTB_300MeV_Data(energy);
        }
      }
      break;
    // kPi0MCv2 for MC and kTestBeamv2 for data
    case 5:
      if (fClusterType == 1|| fClusterType == 3 || fClusterType == 4){
        if(isMC == 0) energy *= FunctionNL_kTestBeamv2(energy);
        else energy *= FunctionNL_kPi0MCv2(energy);
      }
      break;

    // kPi0MCv3 for MC and kTestBeamv4 for data (same as case 02, but with new testbeam param. for > 100 gev)
    case 6:
      if (fClusterType == 1|| fClusterType == 3 || fClusterType == 4){
        if(isMC == 0) energy *= FunctionNL_kTestBeamv4(energy);
        else energy *= FunctionNL_kPi0MCv3(energy);
      }
      break;

    // kPi0MCv6 for MC and kSDMv6 for data
    case 8:
      if (fClusterType == 1|| fClusterType == 3 || fClusterType == 4){
        if(isMC == 0) energy *= FunctionNL_kSDMv6(energy);
        else energy *= FunctionNL_kPi0MCv6(energy);
      }
      break;
    // case 11 of the 8 TeV (LHC15h1 PYTHIA8) nonlinearity as a general case
    case 9:
      if(isMC>0){
        if (fClusterType == 1){
          energy /= FunctionNL_kSDM(energy, 0.96874*0.991*0.9958*0.999, -3.76064, -0.193181);
        }
      }
      break;



//----------------------------------------------------------------------------------------------------------

// *************** 10 + x **** default tender settings - pp

    // NonLinearity pp ConvCalo - only shifting MC - no timing cut
    case 11:
      label_case_11:
      if(isMC>0){
        // 8TeV LHC12x
        //pass1
        if( fCurrentMC==k14e2b ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.983251, -3.44339, -1.70998);

        //pass2
        } else if( fCurrentMC == kPP8T12P2Pyt8 ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.96874*0.991*0.9958*0.999, -3.76064, -0.193181);

        } else if( fCurrentMC == kPP8T12P2Pho ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.969703*0.989*0.9969*0.9991, -3.80387, -0.200546);

        } else if( fCurrentMC == kPP8T12P2JJ ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.974859*0.987*0.996, -3.85842, -0.405277);

        // 2.76TeV LHC11a/LHC13g
        } else if( fCurrentMC==k12f1a || fCurrentMC==k12i3 || fCurrentMC==k15g2 ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.984889*0.995*0.9970, -3.65456, -1.12744);

        } else if(fCurrentMC==k12f1b){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.984384*0.995*0.9970, -3.30287, -1.48516);

        } else if( fCurrentMC==kPP2T11P4JJ || fCurrentMC==k15g1b || fCurrentMC==kPP2T13P1JJ || fCurrentMC==k15a3b ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.981892*0.995*0.9970, -5.43438, -1.05468);

        // 7 TeV LHC10x
        } else if( fCurrentMC==k14j4 ){ //v3
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.974525*0.986*0.999, -4.00247, -0.453046) ;
            energy /= FunctionNL_kSDM(energy, 0.988038, -4.27667, -0.196969);
            energy /= FunctionNL_kSDM(energy, 0.997544, -4.5662, -0.459687);
          } else if(fClusterType==2){
            energy /= 0.9827788048; //const fit
          }
        // 7 TeV LHC11x (using MB and EMC7)
        } else if( fCurrentMC==k14b7 || fCurrentMC==k14k1ab ){
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.923536, -2.63133, -0.152515) ;
            energy /= FunctionNL_kSDM(energy, 0.991343, -4.26797, -0.156959) ;
          }
        // pp 5.02 TeV LHC15n
        // pass2
        } else if( fCurrentMC==k16h8a ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.969799, -4.11836, -0.293151);

        } else if( fCurrentMC==k16h8b ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.969944, -4.02916, -0.366743);

        } else if( fCurrentMC==k16h3 ||fCurrentMC==k16k5a ||  fCurrentMC==k17e2 || fCurrentMC == k18j3) {
          if(fClusterType==1 || (fClusterType == 4 && !isDCal)){
            energy /= FunctionNL_kSDM(energy, 0.959714, -3.77191, -1.41648);
          }
          if(fClusterType==3 || (fClusterType == 4 && isDCal)){
            if(fCurrentMC==k16k5a) {
              energy /= 0.9870110951;
              energy /= FunctionNL_kSDM(energy, 0.992345, -2.33772, -6.1127);
            }
            if(fCurrentMC==k17e2 || fCurrentMC == k18j3) {
              energy /= FunctionNL_kSDM(energy, 0.986513, 0.430032, -10.99999);
              energy /= 0.9908118231;
            }
          }
          // pp 5.02 TeV LHC17pq
          // pass1
        } else if(fCurrentMC==k17l3b || fCurrentMC==k18j2 || fCurrentMC==k17l4b || fCurrentMC==k18b8 || fCurrentMC==k18b10 || fCurrentMC==k18l2) {
          if(fClusterType==1 || (fClusterType == 4 && !isDCal)){
            energy /= FunctionNL_kSDM(energy, 0.95515, -3.19364, -0.936124);
          } else if (fClusterType==3 || (fClusterType == 4 && isDCal)){
            if(fCurrentMC==k17l3b || fCurrentMC==k18j2){
              energy /= FunctionNL_kSDM(energy, 0.982718, -4.93042, -0.48982);
              energy /= FunctionNL_kSDM(energy, 0.992193, -6.97473, 0.0997857);
              energy /= 0.9965307684;
            } else if(fCurrentMC==k17l4b){
              energy /= 0.9840216957*0.9919000094*0.9955886197;
            }
          }
        } else if( fCurrentMC==k16k5b ) {
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.974173, -4.07732, -0.570223);
          if(fClusterType==3) {
            energy /= 0.9872826260;
            energy /= 0.9930726691;
          }
        //pp 13 TeV LHC16 || LHC17 || LHC18
        } else if ( fCurrentMC==kPP13T16P1Pyt8 || fCurrentMC==kPP13T17P1Pyt8 || fCurrentMC==kPP13T18P1Pyt8 || fCurrentMC==kPP13T16P1JJ || fCurrentMC==kPP13T17P1JJ || fCurrentMC==kPP13T18P1JJ || fCurrentMC==kPP13T16P1JJTrigger || fCurrentMC==kPP13T17P1JJTrigger || fCurrentMC==kPP13T18P1JJTrigger){
          if(fClusterType==2) { //13 TeV PCM-PHOS Exponential function fitted
            // energy /= FunctionNL_kSDM(energy, 0.964058, -2.46552, -0.384301); //old
            energy /= FunctionNL_kSDM(energy, 0.966115, -2.7256, -1.02957, 1.0);
            energy /= 1.022224;
          }
          if(fClusterType==1 || fClusterType==3 || fClusterType==4){
            energy /= FunctionNL_OfficialTB_100MeV_MC_V2(energy);
            energy /= FunctionNL_kSDM(energy, 0.987912, -2.94105, -0.273207) ;
          }

        } else if (fCurrentMC==kPP13T16P1Pyt8LowB || fCurrentMC==kPP13T17P1Pyt8LowB || fCurrentMC==kPP13T18P1Pyt8LowB ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.922912, -2.97895, -0.132756);
          if(fClusterType==2) { //13 TeV Low-B PCM-PHOS Exponential function fitted
              energy /= FunctionNL_kSDM(energy, 1.00571, -2.03882, -2.12252);
          }
          if(fClusterType==4){
              energy /= FunctionNL_OfficialTB_100MeV_MC(energy);
              energy /= FunctionNL_kSDM(energy, 1.00892, -3.9119, -0.339742);
          }
        } else fPeriodNameAvailable = kFALSE;

      } else if (isMC == 0){  // Test Beam Non Lin applied on data
        if( fCurrentMC == k16pp13TeV || fCurrentMC == k17pp13TeV || fCurrentMC == k18pp13TeV ){
          if(fClusterType==1 || fClusterType==3 || fClusterType==4){
            energy /= FunctionNL_OfficialTB_100MeV_Data_V2(energy);
          } else if(fClusterType==2) {
            energy /= 1.022224;
          }
        } else if( fCurrentMC == k16pp13TeVLow || fCurrentMC == k17pp13TeVLow || fCurrentMC == k18pp13TeVLow ){
          if(fClusterType==1 || fClusterType==3 || fClusterType==4){
            energy /= FunctionNL_OfficialTB_100MeV_Data(energy);
          } else if(fClusterType==2) {
            energy /= 1.022224;
          }
        }
      }
      break;

    // NonLinearity pp Calo - only shifting MC - no timing cut
    case 12:
      label_case_12:
      if(isMC>0){
        // 8TeV LHC12x
        //pass1
        if( fCurrentMC==k14e2b ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(2.0*energy, 0.967301, -3.1683, -0.653058);

          //pass2
        } else if( fCurrentMC == kPP8T12P2Pyt8 ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.963379*0.9985*0.9992, -3.61217, -0.614043);

        } else if( fCurrentMC == kPP8T12P2Pho ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.96105*0.999*0.9996, -3.62239, -0.556256);

        } else if( fCurrentMC == kPP8T12P2JJ  ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.960596*0.999*0.999, -3.48444, -0.766862);

        // 2.76TeV LHC11a/LHC13g
        } else if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 || fCurrentMC==k15g2 ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(2.0*energy, 0.966151*0.995*0.9981, -2.97974, -0.29463);

        } else if( fCurrentMC==k12f1b ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(2.0*energy, 0.988814*0.995*0.9981, 0.335011, -4.30322);

        } else if( fCurrentMC==kPP2T11P4JJ || fCurrentMC==k15g1b || fCurrentMC==kPP2T13P1JJ || fCurrentMC==k15a3b ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(2.0*energy, 0.979994*0.995*0.9981, -3.24431, -0.760205);

        // 7TeV LHC10x
        } else if(  fCurrentMC==k14j4 ){ //v3
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.962095*0.9991*0.9993, -3.63967, -0.747825) ;
            energy /= FunctionNL_kSDM(energy, 0.988922, -4.47811, -0.132757);
            energy /= FunctionNL_kSDM(energy, 0.99738, -4.82724, -0.281305);
          } else if(fClusterType==2){ //const fit
            energy /= 0.988937947;
          }
        // 5 TeV LHC15n
        //pass2
        } else if( fCurrentMC==k16h8a ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.958994, -4.48233, -0.0314569);

        } else if( fCurrentMC==k16h8b ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.960074, -3.31954, -1.14748);

        //pass3
        } else if( fCurrentMC==k16h3 ||fCurrentMC==k16k5a ||  fCurrentMC==k17e2 || fCurrentMC == k18j3) {
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.958885, -3.78953, -0.575111);
          if(fClusterType==3) energy /= 0.9835764493;

        } else if(fCurrentMC==k17l3b || fCurrentMC==k18j2 || fCurrentMC==k17l4b || fCurrentMC==k18b8 || fCurrentMC==k18b10 || fCurrentMC==k18l2) {
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.95565, -3.39479, -0.510495);
            energy /= 0.9972974486;
          }
          if(fClusterType==3) {
            if(fCurrentMC==k17l3b || fCurrentMC==k18j2){
              energy /= FunctionNL_kSDM(energy, 0.985487, -4.98658, -1.867);
              energy /= 0.9934716565;
            } else if(fCurrentMC==k17l4b){
              energy /= FunctionNL_kSDM(energy, 0.982095, -2.10719, -3.83615);
              energy /= 0.9949694001;
            }
          }

        } else if( fCurrentMC==k16k5b ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.960074, -3.31954, -1.14748);
          if(fClusterType==3) energy /= FunctionNL_kSDM(energy, 0.981191, -1.93399, -2.60859);


        //pp 13 TeV LHC16 || LHC17 || LHC18
        } else if (fCurrentMC==kPP13T16P1Pyt8LowB || fCurrentMC==kPP13T17P1Pyt8LowB || fCurrentMC==kPP13T18P1Pyt8LowB ){
          if(fClusterType==1) energy /= FunctionNL_kSDM(energy, 0.957323, -3.55283, -0.608886);
          if(fClusterType==4){
              energy /= FunctionNL_OfficialTB_100MeV_MC(energy);
              energy /= FunctionNL_kSDM(energy, 1.00714,-3.65528,-0.199527);
          }
        } else if ( fCurrentMC==kPP13T16P1Pyt8 || fCurrentMC==kPP13T17P1Pyt8 || fCurrentMC==kPP13T18P1Pyt8 || fCurrentMC==kPP13T16P1JJ || fCurrentMC==kPP13T17P1JJ || fCurrentMC==kPP13T18P1JJ || fCurrentMC==kPP13T16P1JJTrigger || fCurrentMC==kPP13T17P1JJTrigger || fCurrentMC==kPP13T18P1JJTrigger){
          if(fClusterType==2) { //13 TeV PHOS-PHOS Exponential function fitted
              //energy /= FunctionNL_kSDM(energy, 0.967918, -2.81051, -1.04303, 1.0); //old
              energy /= FunctionNL_kSDM(energy, 0.972774, -2.77133, -1.39596, 1.0);
              energy /= 1.022224;
          }
          if(fClusterType==1 || fClusterType==3 || fClusterType==4){
            energy /= FunctionNL_OfficialTB_100MeV_MC_V2(energy);
            energy /= FunctionNL_kSDM(energy, 0.987912, -2.94105, -0.273207) ;
            energy /= 1.007;
          }

        } else fPeriodNameAvailable = kFALSE;
      } else if (isMC == 0){  // Test Beam Non Lin applied on data
        if( fCurrentMC == k16pp13TeV || fCurrentMC == k17pp13TeV || fCurrentMC == k18pp13TeV ){
          if(fClusterType==1 || fClusterType==3 || fClusterType==4){
            energy /= FunctionNL_OfficialTB_100MeV_Data_V2(energy);
          } else if(fClusterType==2) {
            energy /= 1.022224;
          }
        } else if( fCurrentMC == k16pp13TeVLow || fCurrentMC == k17pp13TeVLow || fCurrentMC == k18pp13TeVLow ){
          if(fClusterType==1 || fClusterType==3 || fClusterType==4){
            energy /= FunctionNL_OfficialTB_100MeV_Data(energy);
          } else if(fClusterType==2) {
            energy /= 1.022224;
          }
        }
      }
      break;

    // NonLinearity ConvCalo - kTestBeamv3 + shifting MC
    case 13:
      if (fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
        energy *= FunctionNL_kTestBeamv3(energy);
        goto label_case_11;// goto previous case for shifting MC
      }
      break;

    // NonLinearity Calo - kTestBeamv3 + shifting MC
    case 14:
      if (fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
        energy *= FunctionNL_kTestBeamv3(energy);
        goto label_case_12;// goto previous case for shifting MC
      }
      break;

    // NonLinearity ConvCalo - kPi0MC + kSDM
    case 15:
      if (fClusterType == 1 || fClusterType == 3){
        // 8TeV LHC12x
        if ( fCurrentMC==k14e2b || fCurrentMC == kPP8T12P2Pyt8 || fCurrentMC == kPP8T12P2Pho  || fCurrentMC == k12pp8TeV || fCurrentMC == kPP8T12P2JJ ){
          energy *= FunctionNL_kPi0MC(energy, 1.0, 0.04979, 1.3, 0.0967998, 219.381, 63.1604, 1.011);
          if(isMC == 0) energy *= FunctionNL_kSDM(energy, 0.9846, -3.319, -2.033);

        // 2.76TeV LHC11a/LHC13g
        } else if ( fCurrentMC == k12f1a || fCurrentMC == k12i3 || fCurrentMC == k15g2 || fCurrentMC == k12f1b ||
                    fCurrentMC == kPP2T11P4JJ || fCurrentMC == k15g1b || fCurrentMC == kPP2T13P1JJ || fCurrentMC == k15a3b ||
                    fCurrentMC == k11pp2760GeV || fCurrentMC == k13pp2760GeV
                  ) {
          energy *= FunctionNL_kPi0MC(energy, 1.0, 0.04123, 1.045, 0.0967998, 219.381, 63.1604, 1.014);
          if(isMC == 0) energy *= FunctionNL_kSDM(energy, 0.9807*0.995*0.9970, -3.377, -0.8535);
        }
        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity Calo - kPi0MC + kSDM
    case 16:
      if (fClusterType == 1 || fClusterType == 3){
        // 8TeV LHC12x
        if ( fCurrentMC==k14e2b  || fCurrentMC == kPP8T12P2Pyt8 || fCurrentMC == kPP8T12P2Pho  || fCurrentMC == k12pp8TeV || fCurrentMC == kPP8T12P2JJ ){
          energy *= FunctionNL_kPi0MC(energy, 1.0, 0.06539, 1.121, 0.0967998, 219.381, 63.1604, 1.011);
          if(isMC == 0) energy *= FunctionNL_kSDM(2.0*energy, 0.9676, -3.216, -0.6828);

        // 2.76TeV LHC11a/LHC13g
        } else if ( fCurrentMC == k12f1a || fCurrentMC == k12i3 || fCurrentMC == k15g2 || fCurrentMC == k12f1b ||
                    fCurrentMC == kPP2T11P4JJ || fCurrentMC == k15g1b || fCurrentMC == kPP2T13P1JJ || fCurrentMC == k15a3b ||
                    fCurrentMC == k11pp2760GeV || fCurrentMC == k13pp2760GeV
                  ) {
          energy *= FunctionNL_kPi0MC(energy, 1.0, 0.06115, 0.9535, 0.0967998, 219.381, 63.1604, 1.013);
          if(isMC == 0) energy *= FunctionNL_kSDM(2.0*energy, 0.9772*0.995*0.9981, -3.256, -0.4449);
        }
        else fPeriodNameAvailable = kFALSE;
      }
      break;

    // New PCM-EMC based SM-wise relative correction applied on data
    case 17:
      if(isMC==0){
        switch (clusterSMID){
          // values determined on LHC17pq
          case 0: energy/=0.996406; break;
          case 1: energy/=0.995672; break;
          case 2: energy/=1.000260; break;
          case 3: energy/=0.998378; break;
          case 4: energy/=0.998326; break;
          case 5: energy/=0.999704; break;
          case 6: energy/=1.001140; break;
          case 7: energy/=0.999142; break;
          case 8: energy/=1.002450; break;
          case 9: energy/=1.002400; break;
          case 10: energy/=1.006850; break;
          case 11: energy/=1.005900; break;
          case 12: energy/=0.999443; break;
          case 13: energy/=0.996047; break;
          case 14: energy/=0.997640; break;
          case 15: energy/=0.996848; break;
          case 16: energy/=1.005100; break;
          case 17: energy/=1.007150; break;
          case 18: energy/=1.006660; break;
          case 19: energy/=1.008700; break;
          default: energy/=1.0; break;
        }
      }
      else
        goto label_case_18;
      break;
    // PCM-EDC based nonlinearity kSDM
    case 18:
      label_case_18:
      if(isMC>0){
         // pp 5 TeV LHC15n+LHC17pq anchored MCs combined
         if( fCurrentMC==k17l3b || fCurrentMC==k18j2 ||  fCurrentMC==k17e2 || fCurrentMC == k18j3 || fCurrentMC == k16h3 || fCurrentMC == k18b8 || fCurrentMC == k18b10 || fCurrentMC==k18l2) {
          if(fClusterType==4){
            energy /= (FunctionNL_kSDM(energy, 0.94723, -3.44986, -0.483821));
          }
        } else fPeriodNameAvailable = kFALSE;
      }
      break;
    case 19:
      if(isMC>0){
         // pp 13 TeV PCM-PHOS NL with PHOS finetuning
         if ( fCurrentMC==kPP13T16P1Pyt8 || fCurrentMC==kPP13T17P1Pyt8 || fCurrentMC==kPP13T18P1Pyt8 || fCurrentMC==kPP13T16P1JJ || fCurrentMC==kPP13T17P1JJ || fCurrentMC==kPP13T18P1JJ || fCurrentMC==kPP13T16P1JJTrigger  || fCurrentMC==kPP13T17P1JJTrigger || fCurrentMC==kPP13T18P1JJTrigger || fCurrentMC==kPP13T16P1Pyt8LowB || fCurrentMC==kPP13T17P1Pyt8LowB || fCurrentMC==kPP13T18P1Pyt8LowB || fCurrentMC==kPP13T16P1JJLowB){
           if(fClusterType==2) { //13 TeV PCM-PHOS Exponential function fitted, corrected by PHOS
               energy /= FunctionNL_kSDM(energy, 0.966115, -2.7256, -1.02957, 1.0);
               energy /= 1.022224;
               energy /= FunctionNL_LinLogConst(energy,  0.374346, 2.08291, 1.12166, -0.33141, 0.00247156, -0.124062, -0.119848);
           }
        } else fPeriodNameAvailable = kFALSE;
      } else if (isMC==0){
          if( fCurrentMC == k16pp13TeV || fCurrentMC == k17pp13TeV || fCurrentMC == k18pp13TeV || fCurrentMC == k16pp13TeVLow || fCurrentMC == k17pp13TeVLow || fCurrentMC == k18pp13TeVLow ){
              if(fClusterType==2) {
                energy /= 1.022224;
              }
          }
      }
      break;

// *************** 20 + x **** modified tender Settings 1 - pp
    // NonLinearity pp ConvCalo - only shifting MC - no timing cut
    case 21:
      label_case_21:
      if(isMC>0){
        // 2.76TeV LHC11a/LHC13g
        if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0443938253, -0.0691830812, -0.1247555443, 1.1673716264, -0.1853095466, -0.0848801702) - 0.0055);
        } else if(fCurrentMC==k15g2){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.1716155406, -0.1962930603, -0.0193959829, 1.0336659741, -0.0467778485, -0.4407662248) - 0.0055);
        } else if(fCurrentMC==k12f1b){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0166321784, -0.0440799552, -0.2611899222, 1.0636538464, -0.0816662488, -0.2173961316) - 0.007);
        } else if( fCurrentMC==kPP2T11P4JJ || fCurrentMC==k15g1b ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.1100193881, -0.1389194936, -0.0800000242, 1.1673716264, -0.1853095466, -0.0848801702) - 0.017);
        } else if( fCurrentMC==kPP2T13P1JJ || fCurrentMC==k15a3b ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0520183153, -0.0806102847, -0.1450415920, 1.0336724056, -0.0467844121, -0.4406992764) - 0.016);
        // 8TeV
        } else if( fCurrentMC == kPP8T12P2Pyt8 ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0654169768, -0.0935785719, -0.1137883054, 1.1814766150, -0.1980098061, -0.0854569214) - 0.0138);
        } else if( fCurrentMC == kPP8T12P2Pho ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0652493513, -0.0929276101, -0.1113762695, 1.1837801885, -0.1999914832, -0.0854569214) - 0.0145);
        } else if( fCurrentMC == kPP8T12P2JJ ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0489259285, -0.0759079646, -0.1239772934, 1.1835846739, -0.1998987993, -0.0854186691) - 0.014);
        // 7 TeV
        } else if( fCurrentMC == k14j4 ){ //v3
          if(fClusterType==1){
            energy /= (FunctionNL_DPOW(energy, 1.1082846035, -0.1369968318, -0.0800000002, 1.1850179319, -0.1999999950, -0.0863054172) - 0.015);
            energy /= FunctionNL_kSDM(energy, 0.988248, -4.26369, -0.208921) ;
            energy /= FunctionNL_kSDM(energy, 0.997359, -4.51031, -0.460041) ;
          } else if(fClusterType==2){
            // ext PHOS + additional correction PCM-PHOS
            energy /= FunctionNL_kSDM(energy, 0.990156, -3.86446, -2.17927);
            goto label_case_01;
          }
        // 7 TeV LHC11x (with MB and EMC7)
        } else if( fCurrentMC==k14b7 || fCurrentMC==k14k1ab ){
          if(fClusterType==1){
            energy /= (FunctionNL_DExp(energy, 0.9898309757, 0.2470261543, -3.0333352573, 1.0932323138, 0.1173201428, -2.0088635359));
            energy /= (FunctionNL_DExp(energy, 1.0830004158, 0.1234703150, -2.0898601923, 1.1073840098, 0.1042803065, -1.8965494640));
          }
        // 5 TeV LHC15n
        //pass2
        } else if( fCurrentMC==k16h8a ){
          if(fClusterType==1) energy /= (FunctionNL_DExp(energy, 0.9831956962, 1.2383793944, -3.2676359751, 1.0121710221, 0.6588125132, -3.1578818630));
        } else if( fCurrentMC==k16h8b ){
          if(fClusterType==1) energy /= (FunctionNL_DExp(energy, 0.9912139474, 0.3721971884, -3.6440765835, 1.0141024579, 0.5574244401, -3.1894624833));
        //pass3
        } else if( fCurrentMC==k16h3 ||fCurrentMC==k16k5a ||  fCurrentMC==k17e2 || fCurrentMC == k18j3 || fCurrentMC==k17l3b || fCurrentMC==k18j2 || fCurrentMC==k17l4b || fCurrentMC==k18b8 || fCurrentMC==k18b10 || fCurrentMC==k18l2) {
          if(fClusterType==1){
            if(fCurrentMC==k17e2 || fCurrentMC == k18j3){
              energy /= (FunctionNL_DExp(energy, 0.9801638660, 0.5771548181, -2.8646883006, 1.0351523739, 0.4332338997, -2.5048158533));
            }else{
              energy /= (FunctionNL_DExp(energy, 0.9814238552, 0.4630567354, -2.9816023028, 1.0293861994, 0.5615679532, -2.3995137175));
            }
          } else if( fClusterType==3) {
            if( fCurrentMC==k16k5a ) {
              energy /= (FunctionNL_DPOW(energy, 0.9943969544,-0.0181151588,-0.4999998851,1.0288066416,-0.0367913727,-0.4995137932));
              energy /= FunctionNL_DPOW(energy, 1.0055560859, -0.0213391278, -0.4999999991, 1.1047136553, -0.1141567995, -0.1573142879);
              energy /= FunctionNL_DPOW(energy, 1.0275381918, -0.0400165029, -0.4999999995, 1.0703233524, -0.0855441426, -0.2099590700);
            } else if( fCurrentMC==k17e2 || fCurrentMC == k18j3 ) {
              energy /= (FunctionNL_DPOW(energy, 0.9943969544,-0.0181151588,-0.4999998851,1.0288066416,-0.0367913727,-0.4995137932));
              energy /= FunctionNL_DPOW(energy, 1.0055560859, -0.0213391278, -0.4999999991, 1.1047136553, -0.1141567995, -0.1573142879);
              energy /= FunctionNL_DPOW(energy, 1.0275381918, -0.0400165029, -0.4999999995, 1.0703233524, -0.0855441426, -0.2099590700);
            } else if( fCurrentMC==k17l3b || fCurrentMC==k18j2 || fCurrentMC==k18b8 || fCurrentMC==k18b10 || fCurrentMC==k18l2) {
              energy /= (FunctionNL_DPOW(energy, 1.1178961735, -0.1395676170, -0.0800000116, 1.1609614690, -0.1698895555, -0.0800001311));
              energy /= 0.9945;
              energy /= (FunctionNL_DPOW(energy, 1.1329092078, -0.1458542910, -0.0800001476, 1.1609614690, -0.1698895555, -0.0800001311));
              energy /= 0.9982018995;
             } else if(fCurrentMC==k17l4b ){
              energy /= (FunctionNL_DPOW(energy, 1.1196043752, -0.1442572200, -0.0800000673, 1.1609614690, -0.1698895555, -0.0800001311));
              energy /= (FunctionNL_DPOW(energy, 1.1580381376,  -0.1746731633,  -0.0800000020,  1.1609614690,  -0.1698895555,  -0.0800001311));
              energy /= 0.9957873506;
              }
          }
        } else if( fCurrentMC==k16k5b ){
          if(fClusterType==1) energy /= (FunctionNL_DExp(energy, 0.9842689920, 0.9150246921, -3.6796298486, 1.0113148506, 0.6876891951, -3.1672234730));
          if(fClusterType==3) {
            energy /= (FunctionNL_DPOW(energy, 1.1343351836,-0.1571288013,-0.0800000607,1.0288066416,-0.0367913727,-0.4995137932));
            energy /= FunctionNL_DPOW(energy, 1.1105555600, -0.1266067088, -0.0800000497, 1.1047136553, -0.1141567995, -0.1573142879);
          }
        //pp 13 TeV LHC16 || LHC17 || LHC18
        } else if ( fCurrentMC==kPP13T16P1Pyt8  || fCurrentMC==kPP13T17P1Pyt8 || fCurrentMC==kPP13T18P1Pyt8 || fCurrentMC==kPP13T16P1JJ || fCurrentMC==kPP13T17P1JJ || fCurrentMC==kPP13T18P1JJ  || fCurrentMC==kPP13T16P1JJTrigger  || fCurrentMC==kPP13T17P1JJTrigger || fCurrentMC==kPP13T18P1JJTrigger){
          if(fClusterType==2) energy /= (FunctionNL_DPOW(energy, 0.9893461252, 0.0541088219, -0.4999999904, 1.0204701327, 0.0010000000, 1.7769590236));
          if(fClusterType==1 || fClusterType==3 || fClusterType==4){
            energy /= FunctionNL_OfficialTB_100MeV_MC_V2(energy);
            energy /= FunctionNL_kSDM(energy, 0.987912, -2.94105, -0.273207) ;
            energy /= 1.00349;
          }

        } else if ( fCurrentMC==kPP13T16P1Pyt8LowB || fCurrentMC==kPP13T17P1Pyt8LowB  || fCurrentMC==kPP13T18P1Pyt8LowB){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0496452471, -0.1047424135, -0.2108759639, 1.1740021856, -0.2000000000, -0.1917378883));
          if(fClusterType==2) energy /= (FunctionNL_DPOW(energy, 1.0167588250, 0.0501002307, -0.8336787497, 0.9500009312, 0.0944118922, -0.1043983134));
          if(fClusterType==4){
              energy /= FunctionNL_OfficialTB_100MeV_MC(energy);
              energy /= (FunctionNL_DExp(energy, 1.0152044323, 1.1015035488, -3.5544630200, 1.0030278679, 0.6774469194, -3.2423591684));
          }
        } else fPeriodNameAvailable = kFALSE;

      } else if (isMC == 0){  // Test Beam Non Lin applied on data
        if( fCurrentMC == k16pp13TeV || fCurrentMC == k17pp13TeV || fCurrentMC == k18pp13TeV ){
          if(fClusterType==1 || fClusterType==3 || fClusterType==4){
            energy /= FunctionNL_OfficialTB_100MeV_Data_V2(energy);
          }
        } else if( fCurrentMC == k16pp13TeVLow || fCurrentMC == k17pp13TeVLow || fCurrentMC == k18pp13TeVLow ){
          if(fClusterType==1 || fClusterType==3 || fClusterType==4){
            energy /= FunctionNL_OfficialTB_100MeV_Data(energy);
          }
        }
      }
      break;

    // NonLinearity pp Calo - only shifting MC - no timing cut
    case 22:
      label_case_22:
      if(isMC>0){
        // 2.76TeV LHC11a/LHC13g
        if(  fCurrentMC==k12f1a || fCurrentMC==k12i3 ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 0.9980625418, -0.0564782662, -0.5, 1.0383412435, -0.0851830429, -0.4999999996) - 0.00175);
        } else if( fCurrentMC==k15g2 ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0795372569, -0.1347324732, -0.1630736190, 1.1614181498, -0.199995361, -0.1711378093) - 0.0035);
        } else if( fCurrentMC==k12f1b ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0232969083, -0.090409434, -0.3592406513, 1.0383412435, -0.0851830429, -0.4999999996) + 0.0007);
        } else if( fCurrentMC==kPP2T11P4JJ || fCurrentMC==k15g1b ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0106037132, -0.0748250591, -0.4999999996, 1.0383412435, -0.0851830429, -0.4999999996) - 0.014);
        } else if( fCurrentMC==kPP2T13P1JJ || fCurrentMC==k15a3b ) {
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0119417393, -0.0755250741, -0.4999999996, 1.1614181498, -0.1999995361, -0.1711378093) - 0.006);
        //8TeV
        } else if( fCurrentMC == kPP8T12P2Pyt8 ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.1389201636, -0.1999994717, -0.1622237979, 1.1603460704, -0.1999999989, -0.2194447313) - 0.0025);
        } else if( fCurrentMC == kPP8T12P2Pho ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0105301622, -0.0732424689, -0.5000000000, 1.0689250170, -0.1082682369, -0.4388156470) - 0.001);
        } else if( fCurrentMC == kPP8T12P2JJ ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 0.9922456908, -0.0551212559, -0.5000000000, 1.0513459039, -0.0894163252, -0.5000000000) + 0.002);
        // 7 TeV
        } else if( fCurrentMC == k14j4 ){ //v3
          if(fClusterType==1){
            energy /= (FunctionNL_DPOW(energy, 1.0074002842, -0.0682543971, -0.4509341085, 1.1224162203, -0.1586806096, -0.2458351112) - 0.003) ;
            energy /= FunctionNL_kSDM(energy, 0.99598, -5.03134, -0.269278) ;
            energy /= FunctionNL_kSDM(energy, 0.997738, -4.91921, -0.377381) ;
          }
        // 5 TeV LHC15n
        //pass2
        } else if( fCurrentMC==k16h8a ){
          if(fClusterType==1) energy /= (FunctionNL_DExp(energy, 0.9747084556, 1.3652950049, -1.7832191813, 1.0039014622, 1.3657547071, -1.7852900827));
        } else if( fCurrentMC==k16h8b ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0193460981, -0.0851635674, -0.4984580141, 1.0588985795, -0.0957023147, -0.4999999998));
        //pass3
        } else if( fCurrentMC==k16h3 ||fCurrentMC==k16k5a ||  fCurrentMC==k17e2 || fCurrentMC == k18j3 || fCurrentMC==k17l3b || fCurrentMC==k18j2 || fCurrentMC==k17l4b || fCurrentMC==k18b8 || fCurrentMC==k18b10 || fCurrentMC==k18l2) {
          if(fClusterType==1){
            if(fCurrentMC==k17e2 || fCurrentMC == k18j3){
              energy /= (FunctionNL_DExp(energy, 0.9781525413, 0.6119550658, -2.5791945651, 1.0124545460, 0.7684107762, -2.2569192409));
            }else{
              energy /= (FunctionNL_DExp(energy, 0.9671756224, 0.9580061524, -2.4592540166, 1.0144265411, 0.7007731928, -2.1689124045));
              energy /= 0.9973908612;
            }
          } else if(fClusterType==3){
            if(fCurrentMC==k17e2 || fCurrentMC == k18j3) energy /= 0.9825370234*0.9993152454;
            else if(fCurrentMC==k16k5a) energy /= 0.9825370234*0.9993152454;
            else if( fCurrentMC==k17l3b || fCurrentMC==k18j2 || fCurrentMC==k18b8 || fCurrentMC==k18b10 || fCurrentMC==k18l2) {
              energy /= (FunctionNL_DExp(energy, 0.9836617483, 0.6608383279, -2.2721654573, 1.0021574658, 0.6692791379, -2.1508057626));
              energy /= 0.9938508694;
            } else if(fCurrentMC==k17l4b){
              energy /= (FunctionNL_DExp(energy, 0.9645652272, 1.0183616681, -2.3081105248, 1.0021574658, 0.6692791379, -2.1508057626));
              energy /= 0.9960109849;
            }
          }
        } else if( fCurrentMC==k16k5b ){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0193460981, -0.0851635674, -0.4984580141, 1.0588985795, -0.0957023147, -0.4999999998));
          if(fClusterType==3) energy /= (FunctionNL_DPOW(energy, 0.9629798154, -0.0178058455, -0.4999999880, 1.1467423891, -0.1999980199, -0.1753999427));

        //pp 13 TeV LHC16 || LHC17 || LHC18
        } else if ( fCurrentMC==kPP13T16P1Pyt8 || fCurrentMC==kPP13T17P1Pyt8 || fCurrentMC==kPP13T18P1Pyt8 || fCurrentMC==kPP13T16P1JJ || fCurrentMC==kPP13T17P1JJ || fCurrentMC==kPP13T18P1JJ || fCurrentMC==kPP13T16P1JJTrigger  || fCurrentMC==kPP13T17P1JJTrigger || fCurrentMC==kPP13T18P1JJTrigger){
          if(fClusterType==1 || fClusterType==3 || fClusterType==4){
            energy /= FunctionNL_OfficialTB_100MeV_MC_V2(energy);
            energy /= FunctionNL_kSDM(energy, 0.987912, -2.94105, -0.273207) ;
            energy /= 1.002;
          }

        } else if (fCurrentMC==kPP13T16P1Pyt8LowB || fCurrentMC==kPP13T17P1Pyt8LowB  || fCurrentMC==kPP13T18P1Pyt8LowB){
          if(fClusterType==1) energy /= (FunctionNL_DPOW(energy, 1.0187401756, -0.0857332791, -0.5000000000, 1.1585209386, -0.1999999989, -0.2646540338));
          if(fClusterType==4){
              energy /= FunctionNL_OfficialTB_100MeV_MC(energy);
              energy /= (FunctionNL_DExp(energy, 1.0027689964, 1.4242237614, -2.3044616171, 0.9872275956, 0.8291951277, -2.5477399129));
          }

        } else fPeriodNameAvailable = kFALSE;
      } else if (isMC == 0){  // Test Beam Non Lin applied on data
        if( fCurrentMC == k16pp13TeV || fCurrentMC == k17pp13TeV || fCurrentMC == k18pp13TeV ){
          if(fClusterType==1 || fClusterType==3 || fClusterType==4){
            energy /= FunctionNL_OfficialTB_100MeV_Data_V2(energy);
          }
        } else if( fCurrentMC == k16pp13TeVLow || fCurrentMC == k17pp13TeVLow || fCurrentMC == k18pp13TeVLow ){
          if(fClusterType==1 || fClusterType==3 || fClusterType==4){
            energy /= FunctionNL_OfficialTB_100MeV_Data(energy);
          }
        }
      }
      break;
    // NonLinearity ConvCalo - kTestBeamv3 + shifting MC
    case 23:
      if (fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
        energy *= FunctionNL_kTestBeamv3(energy);
        goto label_case_21;// goto previous case for shifting MC
      }
      break;

    // NonLinearity Calo - kTestBeamv3 + shifting MC
    case 24:
      if (fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
        energy *= FunctionNL_kTestBeamv3(energy);
        goto label_case_22;// goto previous case for shifting MC
      }
      break;

    // New PCM-EMC nonlinearity with energy squared
    case 27:
      if(isMC>0){
         if( fCurrentMC==k16h3 ||fCurrentMC==k16k5a ||fCurrentMC==k16k5b ||  fCurrentMC==k17e2 || fCurrentMC == k18j3 ) {
          if(fClusterType==1){
            energy /= (FunctionNL_DPOW(energy, 1.1497456392, -0.1999999732, -0.0839303140, 1.1818406492, -0.1999998957, -0.1434322871) + 0.0055);
          }
        } else fPeriodNameAvailable = kFALSE;
      }
      break;
    // PCM-EDC based nonlinearity DExp or DPow
    case 28:
      if(isMC>0){
         // pp 5 TeV LHC15n+LHC17pq anchored MCs combined
         if( fCurrentMC==k17l3b || fCurrentMC==k18j2 ||  fCurrentMC==k17e2 || fCurrentMC == k18j3 ) {
          if(fClusterType==4){
            energy /= (FunctionNL_DExp(energy, 0.9857080359, 0.3587374604, -2.9204912856, 1.0362777077, 0.4371862377, -2.4503073858));
          }
        } else fPeriodNameAvailable = kFALSE;
      }
      break;
// *************** 30 + x **** modified tender Settings 2 - pp
    // PCM-EDC based nonlinearity kSDM
    case 31:
      // apply testbeam nonlinearity (same as case 1) and further fine tuning
      if(fClusterType==1 || fClusterType==3 || fClusterType==4){
        if(isMC){
          energy /= FunctionNL_OfficialTB_100MeV_MC(energy);
          energy /= FunctionNL_kSDM(energy, 0.987534, -3.87469, -0.128085) ;
          if(fCurrentMC==k14j4) energy /= 1.0073711044; // additional finetuning needed prob. due to diff nmb of SM
        } else {
          energy /= FunctionNL_OfficialTB_100MeV_Data(energy);
        }
      }
      break;
    // EDC based nonlinearity kSDM
    case 32:
      // new nonlin for shaper corrected cells
      if(fClusterType==1 || fClusterType==3 || fClusterType==4){
        if(isMC){
          energy /= FunctionNL_OfficialTB_100MeV_MC_V2(energy);
          energy /= FunctionNL_kSDM(energy, 0.987912, -2.94105, -0.273207) ;
          if(fCurrentMC==kPP8T12P2Pho || fCurrentMC==kPP8T12P2Pyt8 || fCurrentMC==kPP8T12P2JJ || fCurrentMC==kPP8T12P2GJLow || fCurrentMC==kPP8T12P2GJHigh) energy /= 0.9875; // additional finetuning needed for pp 8 TeV
        } else {
          energy /= FunctionNL_OfficialTB_100MeV_Data_V2(energy);
        }
      }
      break;
    //
    case 33:
      // apply testbeam nonlinearity - systematic variation 1 (const shift down)
      if(fClusterType==1 || fClusterType==3 || fClusterType==4){
        if(isMC){
          energy /= FunctionNL_OfficialTB_100MeV_MC_V2(energy);
          energy /= FunctionNL_kSDM(energy, 0.977966, -3.0556, -0.27946) ;
        } else {
          energy /= FunctionNL_OfficialTB_100MeV_Data_Sys1(energy);
        }
      }
      break;
    //
    case 34:
      // apply testbeam nonlinearity - systematic variation 2 (tilt low->high)
      if(fClusterType==1 || fClusterType==3 || fClusterType==4){
        if(isMC){
          energy /= FunctionNL_OfficialTB_100MeV_MC_V2(energy);
        } else {
          energy /= FunctionNL_OfficialTB_100MeV_Data_Sys2(energy);
        }
      }
      break;

    // PCM-EDC based nonlinearity for LHC16x,17x,18x pp 13TeV ******* shifting   MC
    case 35:
      if(isMC>0){
        //pp 13 TeV MCs for LHC16 || LHC17 || LHC18
        if ( fCurrentMC==kPP13T16P1Pyt8 || fCurrentMC==kPP13T17P1Pyt8 || fCurrentMC==kPP13T18P1Pyt8 || fCurrentMC==kPP13T16P1JJ || fCurrentMC==kPP13T17P1JJ || fCurrentMC==kPP13T18P1JJ){
          if(fClusterType==4 || fClusterType==1 || fClusterType==3){
            energy /= (FunctionNL_ExpExp(energy, 0.9872432434, 0.3665071019, -2.8842177373, 7.4181896132)/FunctionNL_ExpExp(energy, 1.0469170329, 0.2974710295, -2.4204052267, 7.2038176960));
            energy /= (FunctionNL_DExp(energy, 1.0069808932, 1.0667129308, -1.9891679083, 1.0064686109, 0.9237756300, -2.0610042221));
          }
        }
      }
      break;
    // PCM-EDC based nonlinearity for LHC16x,17x,18x pp 13TeV ******* shifting  data and MC
    case 36:
      if(isMC>0){
        //pp 13 TeV MCs for LHC16 || LHC17 || LHC18
        if ( fCurrentMC==kPP13T16P1Pyt8 || fCurrentMC==kPP13T17P1Pyt8 || fCurrentMC==kPP13T18P1Pyt8 || fCurrentMC==kPP13T16P1JJ || fCurrentMC==kPP13T17P1JJ || fCurrentMC==kPP13T18P1JJ){
          if(fClusterType==4 || fClusterType==1 || fClusterType==3){
            energy /= FunctionNL_ExpExp(energy, 0.9872432434, 0.3665071019, -2.8842177373, 7.4181896132);
          }
        }
      } else if (isMC == 0){
        //pp 13 TeV LHC16 || LHC17 || LHC18
        if( fCurrentMC == k16pp13TeV || fCurrentMC == k17pp13TeV || fCurrentMC == k18pp13TeV ){
          if(fClusterType==1 || fClusterType==3 || fClusterType==4){
            energy /= FunctionNL_ExpExp(energy, 1.0469170329, 0.2974710295, -2.4204052267, 7.2038176960);
          }
        }
      }
      break;
    // PCM-EDC based nonlinearity for LHC16x,17x,18x pp 13TeV ******* shifting   MC
    case 37:
      if( fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
        // TB parametrization from Nico on Martin 100MeV points (final version without. fine tuning)
        if(isMC){
          energy /= FunctionNL_OfficialTB_100MeV_MC_V2(energy);
        } else {
          energy /= FunctionNL_OfficialTB_100MeV_Data_V2(energy);
        }
      }
      break;
    case 38:
      // apply testbeam nonlinearity - systematic variation 0 (const shift up)
      if(fClusterType==1 || fClusterType==3 || fClusterType==4){
        if(isMC){
          energy /= FunctionNL_OfficialTB_100MeV_MC_V2(energy);
          energy /= FunctionNL_kSDM(energy, 0.996916, -2.77626, -0.282107) ;
        } else {
          energy /= FunctionNL_OfficialTB_100MeV_Data_Sys0(energy);
        }
      }
      break;
    //
    case 39:
      // apply testbeam nonlinearity - systematic variation 3 (tilt high->low)
      if(fClusterType==1 || fClusterType==3 || fClusterType==4){
        if(isMC){
          energy /= FunctionNL_OfficialTB_100MeV_MC_V2(energy);
          energy /= FunctionNL_kSDM(energy, 0.986212, -2.63601, -0.27951) ;
      } else {
          energy /= FunctionNL_OfficialTB_100MeV_Data_Sys3(energy);
        }
      }
      break;
// *************** 40 + x **** default tender Settings - pPb
    // NonLinearity LHC13 pPb ConvCalo  - only shifting MC
    case 41:
      label_case_41:
      if(isMC>0){
        if( fCurrentMC==kPPb5T13P2DPMJet || fCurrentMC==kPPb5T13P4DPMJet || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c ){
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.967546, -3.57657, -0.233837) ; // with TM pt dep
            energy /= FunctionNL_kSDM(energy, 0.987513, -4.34641, -0.522125) ;
            energy /= 0.9935;
          }
        } else if( fCurrentMC==kPPb5T13P4JJ || fCurrentMC == kPPb5T13P4JJhigh || fCurrentMC == kPPb5T13P4JJlow  ){
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.965112, -3.45484, -1.33685) ;
          }
        } else if( fCurrentMC==kPPb5T13P2HIJAdd ) {
          if(fClusterType==1){
            energy /=  FunctionNL_kSDM(energy, 0.968868, -3.38407, -0.318188) ;
            energy /=  (FunctionNL_kSDM(energy, 0.987931, -4.13218, -0.583746)*0.9953479301) ;//with TM pt dep
          }
        } else if( fCurrentMC==kPPb5T16EPOS ) {
          if(fClusterType==1 || fClusterType==3 ){
            energy /= FunctionNL_kSDM(energy, 0.951944, -3.38177, -0.597868) ; //  2018 02 20
            energy /= FunctionNL_DPOW(energy, 1.1682295822, -0.1999999973, -0.2088780018, 1.1653083892, -0.1999999998, -0.2697014136 ); //  2018 03 22
          } else if (fClusterType==2){
            energy /= FunctionNL_kSDM(energy, 0.977985, -2.97438, -0.598613) ;
            energy /= 1.02231;
          }
        } else if( fCurrentMC==kPPb5T16DPMJet || fCurrentMC==k17g8a) {
          if(fClusterType==1 || fClusterType==3 || fClusterType== 4){
            energy /= FunctionNL_kSDM(energy, 0.949519, -3.34992, -0.472659) ;   //  2019 08 26
          } else if (fClusterType==2){
            energy /= FunctionNL_kSDM(energy, 0.975133, -2.44519, -1.89808) ; // 2019 08 28
          }
        } else if( fCurrentMC==k18f3bc || fCurrentMC==k18b9b || fCurrentMC==k18b9c ) {
          if(fClusterType==1 ){
            energy /= FunctionNL_kSDM(energy, 0.94058, -3.01613, -0.633012) ;
            energy /= 0.9977804942 ;
          } else if (fClusterType==3){
            energy /= FunctionNL_kSDM(energy, 0.96232, -3.92883, -0.735246) ;
            energy /= FunctionNL_kSDM(energy, 0.968415, -3.71181, -0.169342) ;
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;

    // NonLinearity LHC13 pPb Calo  - only shifting MC
    case 42:
      label_case_42:
      if(isMC>0){
        if( fCurrentMC==kPPb5T13P2DPMJet || fCurrentMC==kPPb5T13P4DPMJet || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c ) {
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.973301, -3.66136, -1.20116) ; //with TM pt dep
            energy /= (FunctionNL_kSDM(energy, 0.987611, -4.14227, -0.282541) * 1.0036264536 );
            energy /= 0.9935;
          }
        } else if( fCurrentMC==kPPb5T13P4JJ || fCurrentMC == kPPb5T13P4JJhigh || fCurrentMC == kPPb5T13P4JJlow) {
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.967048, -2.30814, -2.0672) ;
          }
        } else if( fCurrentMC==kPPb5T13P2HIJAdd ) {
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.962047, -3.18433, -0.586904); //with TM pt dep
            energy /= FunctionNL_kSDM(energy, 0.990771, -4.29086, -0.27403);
          }
        } else if( fCurrentMC==kPPb5T16EPOS ) {
          if(fClusterType==1 || fClusterType==3 ){
            energy /= FunctionNL_kSDM(energy, 0.949402, -3.17052, -0.57999) ; // 2018 02 20
            energy /= FunctionNL_DPOW(energy, 1.1749094095, -0.1999999946, -0.2097130855, 1.1820209963, -0.1999999999, -0.1811167881 ); // 2018 03 22
          }
        } else if( fCurrentMC==kPPb5T16DPMJet || fCurrentMC==k17g8a ) {
          if(fClusterType==1 || fClusterType==3 ){
            energy /= FunctionNL_kSDM(energy, 0.950272, -3.25783, -0.48271) ; //  2019 01 23
          } else if (fClusterType==2) {
            energy /= FunctionNL_kSDM(energy, 0.982663, -2.48524, -0.615548, 1.);
            energy /= FunctionNL_DPOW(energy, 8.4639110152, -7.4514726252, -0.0004532696,  2.3679260731, -1.3567522221, -0.0044275381) ;
            energy /= FunctionNL_kSDM(0.978439, -3.73828, -0.452777, -1);
          }
        } else if( fCurrentMC==k18f3bc || fCurrentMC==k18b9b || fCurrentMC==k18b9c ) {
          if(fClusterType==1 ){
            energy /= FunctionNL_kSDM(energy, 0.959162, -4.58126, -0.495856) ;
            energy /= 0.998 ;
            energy /= FunctionNL_DPOW(energy, 1.0670830446, -0.1052003823, -0.4999999897, 1.0561213150, -0.0918775002, -0.4999999972 ) ;
          } else if (fClusterType==3){
            energy /= FunctionNL_kSDM(energy, 0.943033, -3.96729, -0.383147) ;
            energy /= 0.99858 ;
          }
        } else fPeriodNameAvailable = kFALSE;
      }
      break;

    // NonLinearity LHC13 pPb ConvCalo  - kTestBeamv3 + shifting MC
    case 43:
      if (fClusterType == 1 || fClusterType == 3){
        energy *= FunctionNL_kTestBeamv3(energy);
        goto label_case_41;// goto previous case for shifting MC
      }
      break;

    // NonLinearity LHC13 pPb Calo  - kTestBeamv3 + shifting MC
    case 44:
      if (fClusterType == 1 || fClusterType == 3){
        energy *= FunctionNL_kTestBeamv3(energy);
        goto label_case_42;// goto previous case for shifting MC
      }
      break;
    // NonLinearity LHC13 pPb ConvCalo  - applying f^2
    case 45:
      if(isMC>0){
        if( fCurrentMC==kPPb5T13P2DPMJet || fCurrentMC==kPPb5T13P4DPMJet || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c ){
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 0.967546, -3.57657, -0.233837) ;
            //apply again the same
            energy /= FunctionNL_kSDM(energy, 0.967546, -3.57657, -0.233837) ;
          }
        } else if( fCurrentMC==kPPb5T13P2HIJAdd ) {
          if(fClusterType==1){
            energy /=  FunctionNL_kSDM(energy, 0.968868, -3.38407, -0.318188) ;
            //apply again the same
            energy /=  FunctionNL_kSDM(energy, 0.968868, -3.38407, -0.318188) ;
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;
    case 47:
      if (fClusterType == 1|| fClusterType == 3 || fClusterType == 4){
        if(isMC == 0)
          energy *= FunctionNL_kTestBeamMod(energy, 0.972947, 0.986154, 0.214860, 0.717724, 0.069200, 155.497605, 48.868069);
        else {
          energy *= FunctionNL_kPi0MCMod(energy, 1.004055, 1.009121, 0.083153, 1.444362, 0.100294, 416.897753, 324.246101);
          if( fCurrentMC==k18f3bc || fCurrentMC==k18b9b || fCurrentMC==k18b9c ) {
            if(fClusterType==1 || (fClusterType == 4 && !isDCal)){
              energy /= FunctionNL_SPOW(energy, 0.9911279061, 0.0046281250, -5.1486761047) ;
            } else if (fClusterType==3 || (fClusterType == 4 && isDCal)){
              energy /= FunctionNL_kSDM(energy, 0.976447, -3.88631, -0.337078) ;
            }
          } else fPeriodNameAvailable = kFALSE;
        }
      }
      break;
    case 48:
      if (fClusterType == 1|| fClusterType == 3 || fClusterType == 4){
        if(isMC == 0)
          energy *= FunctionNL_kTestBeamMod(energy, 0.972947, 0.986154, 0.214860, 0.717724, 0.069200, 155.497605, 48.868069);
        else {
          energy *= FunctionNL_kPi0MCMod(energy, 1.004055, 1.009121, 0.083153, 1.444362, 0.100294, 416.897753, 324.246101);
          if( fCurrentMC==k18f3bc || fCurrentMC==k18b9b || fCurrentMC==k18b9c ) {
            if(fClusterType==1  || (fClusterType == 4 && !isDCal)){
              energy /= FunctionNL_SPOW(energy, 1.0010549936, 0.0155962913, -2.5907698895) ;
            } else if (fClusterType==3 || (fClusterType == 4 && isDCal)){
              energy /= FunctionNL_SPOW(energy, 0.9905473008, 0.0157137852, -2.3456519285) ;
            }
          } else fPeriodNameAvailable = kFALSE;
        }
      }
      break;

// *************** 50 + x **** modified tender Settings 1 - pPb
    // NonLinearity LHC13 pPb ConvCalo  - only shifting MC
    case 51:
      label_case_51:
      if(isMC>0){
        if( fCurrentMC==kPPb5T13P2DPMJet || fCurrentMC==kPPb5T13P4DPMJet || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c || fCurrentMC == kPPb5T13P4JJhigh || fCurrentMC == kPPb5T13P4JJlow || fCurrentMC == kPPb5T13P4JJ) {
          if(fClusterType==1){
            energy /= FunctionNL_DExp(energy, 0.9910691195, 0.4901455923, -3.6647921806, 1.0255088817, 0.3070452373, -2.9149185308); //with TM pt dep
            energy /= FunctionNL_kSDM(energy, 0.989111, -4.26219, -0.819192);
            energy /= 0.9935;
          } else if(fClusterType==2){
            energy /= ( 0.994914734 * 0.9964 * 0.9975659906); // additional factors
          }
        } else if( fCurrentMC==kPPb5T13P2HIJAdd ) {
          if(fClusterType==1){
            energy /= FunctionNL_DExp(energy, 0.9978241421, 0.2054669115, -3.7888984452, 1.0255088817, 0.3070452373, -2.9149185308) ; //with TM pt dep
            energy /= (FunctionNL_kSDM(energy, 0.986673, -4.14594, -0.450765)* 0.9953727823);
          } else if(fClusterType==2){
            energy /= ( 0.993485*0.9971126333 );
          }
        } else if( fCurrentMC==kPPb5T16EPOS ) {
          if(fClusterType==1 || fClusterType==3 ){
            energy /= FunctionNL_DExp(energy, 0.9772393830, 0.9651600903, -2.8485741777, 1.0436698408, 0.4584792411, -2.3634185342); // 2018 02 20
            energy /= FunctionNL_DPOW(energy, 1.1771657563, -0.1999999970, -0.1955919316, 1.1820209963, -0.1999999999, -0.1811167881 ) ; // 2018 03 22
          } else if(fClusterType==2){
            energy /= (0.949117*1.02231) ; //first iteration with constant
          }
        } else if( fCurrentMC==kPPb5T16DPMJet || fCurrentMC==k17g8a ) {
          if(fClusterType==1 || fClusterType==3 ){
            energy /= FunctionNL_DExp(energy, 0.9822883129, 0.4885659420, -2.8847511641, 1.0461791627, 0.4463927828, -2.2994619341); // 2018 02 20
            energy /= FunctionNL_kSDM(energy, 0.992184, -4.22097, -0.561982); // 2018 03 22
            energy /= 0.9907587828 ; // 2018 12 09
          } else if(fClusterType==2){
            energy /= ( 0.997*0.9965200155 ); // additional factors
          }
        } else if(fCurrentMC==k17e2 || fCurrentMC == k18j3 || fCurrentMC == k16h3) {
          if(fClusterType==2){
            energy /= FunctionNL_kSDM(energy, 0.991778, -2.60609, -1.63899);
          }
        } else if(fCurrentMC==k17l3b || fCurrentMC==k18j2 || fCurrentMC==k17l4b || fCurrentMC == k18b8 || fCurrentMC == k18b10 || fCurrentMC==k18l2) {
          if(fClusterType==2){
            energy /= FunctionNL_kSDM(energy, 0.991778, -2.60609, -1.63899);
          }
        } else if( fCurrentMC==k18f3bc || fCurrentMC==k18b9b || fCurrentMC==k18b9c ) {
          if(fClusterType==1 ){
            energy /= FunctionNL_DExp(energy, 0.9725126980, 0.6103856108, -3.0948482125, 1.0353908764, 0.5838986758, -2.3256942180 );
            energy /= FunctionNL_DExp(energy, 1.0275565185, 0.7150701150, -2.2907879895, 1.0356734923, 0.5692435467, -2.3326775814 );
          } else if (fClusterType==3){
            energy /= FunctionNL_DPOW(energy, 0.9984202064, -0.0234949772, -0.4999999701, 1.0436429747, -0.0402608501, -0.4999988777 );
            energy /= FunctionNL_DPOW(energy, 1.0206520740, -0.0284859767, -0.4999999308, 1.0418924549, -0.0385631705, -0.4999999717 );
          }
        } else if ( fCurrentMC==kPP13T16P1Pyt8  || fCurrentMC==kPP13T17P1Pyt8 || fCurrentMC==kPP13T18P1Pyt8 || fCurrentMC==kPP13T16P1JJ || fCurrentMC==kPP13T17P1JJ || fCurrentMC==kPP13T18P1JJ ){
            if(fClusterType==2){
              energy /= FunctionNL_kSDM(energy, 0.991778, -2.60609, -1.63899);
            }
        } else fPeriodNameAvailable = kFALSE;
      }
      break;


    // NonLinearity LHC13 pPb Calo  - only shifting MC
    case 52:
      label_case_52:
      if(isMC>0){
        if( fCurrentMC==kPPb5T13P2DPMJet || fCurrentMC==kPPb5T13P4DPMJet || fCurrentMC == k16c3a || fCurrentMC == k16c3b || fCurrentMC == k16c3c || fCurrentMC == kPPb5T13P4JJhigh || fCurrentMC == kPPb5T13P4JJlow || fCurrentMC == kPPb5T13P4JJ) {
          if(fClusterType==1){
            energy /= FunctionNL_DExp(energy, 0.9795532189, 0.8578583955, -2.3447892540, 1.0165873637, 0.6999387334, -2.1324782465) ;//with TM pt dep
            energy /= (FunctionNL_kSDM(energy, 0.990609, -4.37834, -0.304314) * 1.0040232773) ;
            energy /= 0.9935;
          } else if(fClusterType==2) {
            energy /= (FunctionNL_DExp(energy, 1.0154938040, 0.3062978125, -3.9089772679, 1.0061692542, 513.7621552761, -3566.4426936867 ) * 0.996512);
          }
        } else if( fCurrentMC==kPPb5T13P2HIJAdd ) {
          if(fClusterType==1){
            energy /= FunctionNL_DExp(energy, 0.9746342307, 0.9576270870, -2.5098585110, 1.0165871862, 0.6999571530, -2.1324658480) ; //with TM pt dep
            energy /= FunctionNL_kSDM(energy, 0.993562, -4.52817, -0.366368) ;
          } else if(fClusterType==2){
            energy /= (FunctionNL_DPOW(energy, 1.0154784875, -0.0161589457, -0.4999999976, 1.0086650887, -0.0010000001, -0.0800000139 ) * 0.9983468115 );
          }
       } else if( fCurrentMC==kPPb5T16EPOS ) {
          if(fClusterType==1|| fClusterType==3 ){
            energy /= FunctionNL_DExp(energy, 0.9707112053, 1.4050445187, -2.0357906356, 1.0241095707, 0.9217457498, -1.9020815528) ;//2018 02 20
            energy /= FunctionNL_DPOW(energy, 1.1753551048, -0.1999999981, -0.2060941701, 1.1822845437, -0.2, -0.1793151582 ) ; // 2018 03 22
          } else if(fClusterType==2) {
            energy /= (FunctionNL_DExp(energy, 1.0154938040, 0.3062978125, -3.9089772679, 1.0061692542, 513.7621552761, -3566.4426936867 ) * 0.996512);
          }
       } else if( fCurrentMC==kPPb5T16DPMJet || fCurrentMC==k17g8a ) {
          if(fClusterType==1|| fClusterType==3 ){
            energy /= FunctionNL_DExp(energy, 0.9706706146, 0.9781531357, -2.5633710383, 1.0355397924, 0.6750800461, -1.9817285526) ;//2018 02 20
            energy /= FunctionNL_DPOW(energy, 1.1717517490, -0.1999999942, -0.2126460833, 1.1820209963, -0.1999999999, -0.1811167881 ) ; // 2018 03 22
            energy /= 0.992867 ; // 2018 12 09
          } else if (fClusterType==2){
            energy /= FunctionNL_DExp(energy, 1.0135811994, 0.8812760922, -2.8245351546, 1.0202921434, 0.6804648476, -4.0012344773, -1., 1.);
            energy /= FunctionNL_DExp(energy, 1.0160625375, 0.7320554259, -4.8983050388, 1.0204802067, 0.6653009626, -3.9947369662, 1., 1.);
          } else if(fClusterType==2) {
            energy /= (FunctionNL_DExp(energy, 1.0154938040, 0.3062978125, -3.9089772679, 1.0061692542, 513.7621552761, -3566.4426936867 ) * 0.996512);
          }
       } else if(fCurrentMC==k17l3b || fCurrentMC==k18j2 || fCurrentMC==k17l4b || fCurrentMC == k18b8 || fCurrentMC == k18b10 || fCurrentMC==k18l2) {
          if(fClusterType==2){
            energy /= FunctionNL_kSDM(energy, 1.01008, -2.41975, -1.43326);
          }
       } else if( fCurrentMC==k18f3bc || fCurrentMC==k18b9b || fCurrentMC==k18b9c ) {
          if(fClusterType==1 ){
            energy /= FunctionNL_DPOW(energy, 1.0115133381, -0.0800376623, -0.4999999999, 1.0586512028, -0.0940942550, -0.4999999982 );
            energy /= 0.9963;
          } else if (fClusterType==3){
            energy /= FunctionNL_DExp(energy, 0.9615753216, 1.0083998720, -2.3788968272, 1.0114418980, 1.0697887890, -2.1228591720 );
            energy /= 0.993;
          }
        } else fPeriodNameAvailable = kFALSE;
      }
      break;
    // NonLinearity LHC13 pPb Calo  - kTestBeamv3 + shifting MC
    case 53:
      if (fClusterType == 1 || fClusterType == 3){
        energy *= FunctionNL_kTestBeamv3(energy);
        goto label_case_51;// goto previous case for shifting MC
      } else if (fClusterType == 2) { // PHOS case: shift data and MC to pi0 mass (Calo-Calo)
          if ( fCurrentMC == kPPb5T13P4DPMJet ){ // RUN1, MB MC
            energy /= FunctionNL_SPOW(energy, 1.00132, -0.00598632, -2.167);
          } else if ( fCurrentMC == kPPb5T13P4JJ || fCurrentMC == kPPb5T13P4JJhigh || fCurrentMC == kPPb5T13P4JJlow || fCurrentMC == k16c3a || fCurrentMC == k16c3b) { // RUN1, JJ MC
            energy /= FunctionNL_SExp(energy, 1.00355, 756.545, -20.3831, -1.);
          } else if ( fCurrentMC == k13pPb5023GeV ) { // RUN1, data
            energy /= FunctionNL_SPOW(energy, 0.997167, -0.000759949, -5.04513);
          } else if ( fCurrentMC == kPPb5T16DPMJet ) { // RUN2, MB MC
            energy /= FunctionNL_SExp(energy, 0.991632, 1.22134, -2.97826, -1.);
          } else if ( fCurrentMC == k17g8a ) { // RUN2, JJ MC
            energy /= 1.0145754170;
          } else if ( fCurrentMC == k16pPb5023GeV ) { // RUN2, data
            energy /= FunctionNL_SExp(energy, 1.02521, 0.59747, -4.17083);
          }
      }
      break;

    // NonLinearity LHC13 pPb Calo  - kTestBeamv3 + shifting MC
    case 54:
      if (fClusterType == 1 || fClusterType == 3){
        energy *= FunctionNL_kTestBeamv3(energy);
        goto label_case_52;// goto previous case for shifting MC
      } else if (fClusterType == 2) { // PHOS case: shift data and MC to pi0 mass (PCM-Calo)
          if ( fCurrentMC == kPPb5T13P4DPMJet ){ // RUN1, MB MC
            energy /= FunctionNL_SExp(energy, 1.00703, 16753.6, -23418.6, -1);
          } else if ( fCurrentMC == kPPb5T13P4JJ || fCurrentMC == kPPb5T13P4JJhigh || fCurrentMC == kPPb5T13P4JJlow || fCurrentMC == k16c3a || fCurrentMC == k16c3b) { // RUN1, JJ MC
            energy /= FunctionNL_SExp(energy, 0.996626, 1.99301, -2.33312, -1.);
          } else if ( fCurrentMC == k13pPb5023GeV ) { // RUN1, data
            energy /= FunctionNL_SExp(energy, 1.0128, 11.7606, 1.10537, -1.);
          } else if ( fCurrentMC == kPPb5T16DPMJet ) { // RUN2, MB MC
            energy /= FunctionNL_SPOW(energy, 0.965941, 0.0458828, -0.5);
          } else if ( fCurrentMC == k17g8a ) { // RUN2, JJ MC
            energy /= FunctionNL_SExp(energy, 0.983876, 0.851007, -3.41359, -1);
          } else if ( fCurrentMC == k16pPb5023GeV ) { // RUN2, data
            energy /= 1.0268;
          }
      }
      break;

    // NonLinearity LHC13 pPb Calo  - kTestBeamv3 /kMCTestBeam+shift from PCM-EMC
    case 55:
      if (fClusterType == 1 || fClusterType == 3){
        if(isMC == 0) {
          energy *= FunctionNL_kTestBeamv3(energy);
        } else{
          energy *= FunctionNL_kPi0MCv3(energy);
          if( fCurrentMC==kPPb5T16EPOS )          energy /= 0.9854013683;
          else if ( fCurrentMC==kPPb5T16DPMJet )  energy /= 0.9818891524;
        }

      }
      break;
    // NonLinearity LHC13 pPb Calo  - kTestBeamv3 /kMCTestBeam+shift from EMC
    case 56:
      if (fClusterType == 1 || fClusterType == 3){
        if(isMC == 0) {
          energy *= FunctionNL_kTestBeamv3(energy);
        } else{
          energy *= FunctionNL_kPi0MCv3(energy);
          if( fCurrentMC==kPPb5T16EPOS )          energy /= 0.9887044419;
          else if ( fCurrentMC==kPPb5T16DPMJet )  energy /= 0.9891917142;
        }
      }
      break;
    case 57:
      if (fClusterType == 1|| fClusterType == 3 || fClusterType == 4){
        if(isMC == 0)
          energy *= FunctionNL_kTestBeamMod(energy, 0.972947, 0.986154, 0.214860, 0.717724, 0.069200, 155.497605, 48.868069);
        else {
          energy *= FunctionNL_kPi0MCMod(energy, 1.004055, 1.009121, 0.083153, 1.444362, 0.100294, 416.897753, 324.246101);
          if( fCurrentMC==k18f3bc || fCurrentMC==k18b9b || fCurrentMC==k18b9c ) {
            if(fClusterType==1  || (fClusterType == 4 && !isDCal)){
              energy /= FunctionNL_DPOW(energy, 0.9952456542, 0.0138341175,  -1.2855539410, 1.0092843163, 0.0123126175, -4.1691957278);
              energy /= FunctionNL_DExp(energy, 1.0108067622, 3.7654440410, -13.6062772214, 1.0733438107, 0.0784398968, -2.5492190281) ;
            } else if (fClusterType==3 || (fClusterType == 4 && isDCal)){
              energy /= FunctionNL_DExp(energy, 0.9968158185, 15.1688406741, 3.6964381952, 1.0242323030, 0.1955714892, -3.8051558657) ;
              energy /= FunctionNL_DExp(energy, 1.0108067622, 3.7654440410, -13.6062772214, 1.0733438107, 0.0784398968, -2.5492190281) ;
            }
          } else fPeriodNameAvailable = kFALSE;
        }
      }
      break;
    case 58:
      if (fClusterType == 1|| fClusterType == 3 || fClusterType == 4){
        if(isMC == 0)
          energy *= FunctionNL_kTestBeamMod(energy, 0.972947, 0.986154, 0.214860, 0.717724, 0.069200, 155.497605, 48.868069);
        else {
          energy *= FunctionNL_kPi0MCMod(energy, 1.004055, 1.009121, 0.083153, 1.444362, 0.100294, 416.897753, 324.246101);
          if( fCurrentMC==k18f3bc || fCurrentMC==k18b9b || fCurrentMC==k18b9c ) {
            if(fClusterType==1 || (fClusterType == 4 && !isDCal) ){
              energy /= FunctionNL_DPOW(energy, 0.9975703755, -0.0069259279, -4.0620824681, 0.9968577540, 0.0085314562, -2.0417077818 );
            } else if (fClusterType==3 || (fClusterType == 4 && isDCal)){
              energy /= FunctionNL_DPOW(energy, 0.9883179003, -0.0022477965, -7.4912876804, 0.9983309708, 0.0123415948, -2.2360698158 );
            }
          } else fPeriodNameAvailable = kFALSE;
        }
      }
      break;
    case 59: //PHOS data shift + finetuning based on symmetric decays PHOS
      if (fClusterType == 2){
        if(isMC == 0){
          energy /= 1.012;
        } else {
          energy /= FunctionNL_DPOW(energy, 0.9566250057, 0.0443087138, -0.3580211849, 1.0135534518, -0.0003626298, -3.6314926044 );
          energy /= 1.012;
        }
      }
      break;

// *************** 60 + x **** modified tender Settings 2 - pPb
// PCM-EDC based nonlinearity kSDM
    case 61:
      // apply testbeam nonlinearity (same as case 1) but with resolution uncertainy
      if(isMC){
        energy /= FunctionNL_OfficialTB_100MeV_MC(energy);
      } else {
        energy /= FunctionNL_OfficialTB_100MeV_Data(energy);
      }
      break;
    // EDC based nonlinearity kSDM
    case 62:
      // apply testbeam nonlinearity (same as case 2) but with resolution uncertainy
      if(isMC){
        energy /= FunctionNL_OfficialTB_50MeV_MC(energy);
      } else {
        energy /= FunctionNL_OfficialTB_50MeV_Data(energy);
      }
      break;
    // PCM-EDC based nonlinearity DExp or DPow
    case 63:
      // apply testbeam nonlinearity (same as case 3) but with resolution uncertainy
      if(isMC){
        energy /= FunctionNL_OfficialTB_150MeV_MC(energy);
      } else {
        energy /= FunctionNL_OfficialTB_150MeV_Data(energy);
      }
      break;
    // EDC based nonlinearity DExp or DPow
    case 64:
      // apply testbeam nonlinearity (same as case 4) but with resolution uncertainy
      if (fClusterType == 2){
        if(isMC){
          energy *= 1.0245*(1.- 0.013*TMath::Exp(energy*-0.5));
        }
      }
      break;
    case 65: //50MeV TB update
      if (fClusterType == 1|| fClusterType == 3 || fClusterType == 4){
        if(isMC == 0)
          energy *= FunctionNL_kTestBeamMod(energy, 0.972947, 0.986154, 0.214860, 0.717724, 0.069200, 155.497605, 48.868069);
        else energy *= FunctionNL_kPi0MCMod(energy, 1.004055, 1.009121, 0.083153, 1.444362, 0.100294, 416.897753, 324.246101);
      }
      break;
    case 66: //100MeV TB update
      if (fClusterType == 1|| fClusterType == 3 || fClusterType == 4){
        if(isMC == 0) energy *= FunctionNL_kTestBeamMod(energy, 0.978526, 0.980207, 0.311926, 0.596070, 0.062173, 160.714850, 41.312966);
        else energy *= FunctionNL_kPi0MCMod(energy, 1.011841, 0.999134, 0.110756, 2.736180, 0.002017, 110.402252, -73.996864);
      }
      break;
    case 67: //EVI TB + EMC Corr update
      if (fClusterType == 1|| fClusterType == 3 || fClusterType == 4){
        if(isMC == 0) energy *= FunctionNL_kTestBeamv4(energy);
        else{
          energy *= FunctionNL_kPi0MCv3(energy);
          energy /= FunctionNL_DPOW(energy, -0.8802886739, 1.8764944987, -0.0020594487, 0.9891399006, 0.0139889085, -2.0388063034);
        }
      }
      break;
    case 68: //50MeV TB + EMC Corr update
      if (fClusterType == 1|| fClusterType == 3 || fClusterType == 4){
        if(isMC == 0)
          energy *= FunctionNL_kTestBeamMod(energy, 0.972947, 0.986154, 0.214860, 0.717724, 0.069200, 155.497605, 48.868069);
        else{
          energy *= FunctionNL_kPi0MCMod(energy, 0.898861, 1.109502, 0.083507, 1.910893, 0.285469, 469.340995, 3605.837798);
          energy /= FunctionNL_DPOW(energy, 0.9969964995, -0.0023796998, -6.1181405682, 0.9926523280, 0.0080737917, -3.6584464200);
        }
      }
      break;
    case 69: //PHOS data shift + finetuning based on PCM-PHOS
      if (fClusterType == 2){
        if(isMC == 0){
          energy /= 1.012;
        } else {
          energy /= FunctionNL_DPOW(energy, 0.9516188999, 0.0430079212, -0.4999998720, 1.0095759080, 0.0010000001, 0.0800000000 );
          energy /= 1.012;
        }
      }
      break;

// *************** 70 + x **** default tender Settings - PbPb

    // NonLinearity LHC11h - PbPb 2.76TeV - 0-10% centrality
    case 71:
      if(isMC>0){
        if( fCurrentMC==k14a1 ){
          if (fClusterType == 1 || fClusterType == 3){
            energy /= 0.972607; //
          }
        } else if( fCurrentMC==kXeXe5T17HIJING ){
          if (fClusterType == 1 ){
            energy /= FunctionNL_kSDM(energy, 0.973646, -0.901289, -4.32682) ; // based on peripheral XeXe
          } else if (fClusterType == 2 ){
            energy /= (0.9599764493*0.9873); // based on peripheral XeXe
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;
    // NonLinearity LHC11h - PbPb 2.76TeV - 20-50% centrality
    case 72:
      if(isMC>0){
        if( fCurrentMC==k14a1 ){
          if(fClusterType==1){
            energy /= FunctionNL_kSDM(energy, 1.00926, -2.42107, -1.60995);
          }
        } else if( fCurrentMC==kXeXe5T17HIJING ){
          if (fClusterType == 1 ){
            energy /= 0.9737701536; // based on peripheral XeXe
          } else if (fClusterType == 2 ){
            energy /= (0.9350697962*1.01); // based on semi-central XeXe
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;

    // *************** 80 + x **** modified tender Settings 1 - PbPb

    // NonLinearity LHC15o PbPb ConvCalo  - only shifting MC
    case 81:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= (FunctionNL_kSDM(energy, 0.91336, -2.57765, -0.334976) ); //updated 7 Dec 2018 - peripheral (60-80%) pcmcalo correction
          }
          if(fClusterType==2){
            energy /= FunctionNL_kSDM(energy, 1.02357, -2.19441, -3.10045); //updated on 2019 01 18 - PCMPHOS based on centrality 20-50%
            energy /= FunctionNL_kSDM(energy, 1.00081, -2.09128, -2.41587); //updated on 2019 02 18 - PCMPHOS based on centrality 20-50%
            energy /= 0.9835; // 2019 08 30 - PCMPHOS
          }
        } else if( fCurrentMC==kPbPb5T18HIJING ){
          if (fClusterType == 2 ){
            energy /= FunctionNL_kSDM(energy, 0.991778, -2.60609, -1.63899); // added on 2019 07 29 , based on pp 5TeV
            energy /= 1.0073; // 2019 08 30 - PCMPHOS
          }
        } else if( fCurrentMC==kXeXe5T17HIJING ){
          if (fClusterType == 1 ){
            energy /= (FunctionNL_DPOW(energy, 1.0547527663, -0.0927180446, -0.0800012482, 1.0254208020, -0.0345156682, -0.4999999199)); // based on peripheral XeXe
          } else if (fClusterType == 2 ){
            energy /= (0.9764119296*0.9794); // based on peripheral XeXe
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;

    // NonLinearity LHC15o PbPb Calo  - only shifting MC
    case 82:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= 0.96 ; //updated 7 Dec 2018 - peripheral (60-80%) calo correction
          }
        } else if( fCurrentMC==kXeXe5T17HIJING ){
          if (fClusterType == 1 ){
            energy /= (FunctionNL_DPOW(energy, 1.1223479533, -0.1999999659, -0.1954398178, 1.0373041075, -0.0913333671, -0.5000000000)); // based on peripheral XeXe
          } else if (fClusterType == 2 ){
            energy /= (FunctionNL_DExp(energy, 0.9840385879, 0.4801589926, -2.8482099501, 1.0214220397, 5.8987542970, -12.5701079799)*1.0148) ; // based on  semi-central XeXe
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;
      // NonLinearity LHC15o PbPb Calo  - only shifting MC - MB
    case 83:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1 || fClusterType==4){
            energy /= (FunctionNL_kSDM(energy, 0.91336, -2.57765, -0.334976) ); //updated 7 Dec 2018 - peripheral (60-80%) pcmcalo correction
            energy /= (FunctionNL_DPOW(energy, 1.1698155939, -0.1999998275, -0.2656304851, 1.1754720500, -0.1999999998, -0.1467153129) ); //additional calo correction 20190215
          }
        } else if( fCurrentMC==kPbPb5T18HIJING ){
          if(fClusterType==1 || fClusterType==4){
            energy /= 0.9633; // added on 2019 06 18 , based on PCM-EMC peri
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;
      // NonLinearity LHC15o PbPb Calo  - only shifting MC - 0-10%
    case 84:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= (0.954566) ;
            energy /= FunctionNL_DExp(energy, 1.0548582854, 1.5096237243, -1.6079078305, 1.0538380642, 124049.7, -38409.5) ;
            energy /= (0.9669) ;
            energy /= FunctionNL_DExp(energy, 1.0300358811, 0.8758778788, -2.4142808809, 1.0279722682, 0.5471598153, -3.5373985248) ;
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;
      // NonLinearity LHC15o PbPb Calo  - only shifting MC - 10-20%
    case 85:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= (0.9365) ;
            energy /= FunctionNL_DExp(energy, 1.0380275426, 0.7534354400, -2.2110408210, 1.0408879042, 0.4399353376, -2.9554918759) ;
            energy /= FunctionNL_DExp(energy, 1.0300358811, 0.8758778788, -2.4142808809, 1.0279722682, 0.5471598153, -3.5373985248) ;
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;
      // NonLinearity LHC15o PbPb Calo  - only shifting MC - 20-50%
    case 86:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= (0.948553) ;
            energy /= (0.998) ;
            energy /= FunctionNL_DExp(energy, 1.0300358811, 0.8758778788, -2.4142808809, 1.0279722682, 0.5471598153, -3.5373985248) ;
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;
      // NonLinearity LHC15o PbPb Calo  - only shifting MC - 50-90%
    case 87:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= (0.95306) ;
            energy /= FunctionNL_DExp(energy, 1.0300358811, 0.8758778788, -2.4142808809, 1.0279722682, 0.5471598153, -3.5373985248) ;
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;

// *************** 90 + x **** modified tender Settings 2 - PbPb

      // NonLinearity LHC15o PbPb ConvCalo  - only shifting MC
    case 91:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= FunctionNL_DExp(energy, 1.0026971373, -0.0320283624, -0.4999999953, 1.0750656618, -0.0855019990, -0.4571523301);
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;

      // NonLinearity LHC15o PbPb Calo  - only shifting MC
    case 92:
      if(isMC>0){
        if( fCurrentMC== kPbPb5T15HIJING){
          if(fClusterType==1){
            energy /= FunctionNL_DExp(energy, 1.0541217488, -0.1111428177, -0.4999999983, 1.0782958817, -0.0706389211, -0.4999999959);
          }
        } else {
          fPeriodNameAvailable = kFALSE;
        }
      }
      break;



//----------------------------------------------------------------------------------------------------------

    default:
      AliFatal(Form("NonLinearity correction not defined for cut: '%d' ! Returning...",fSwitchNonLinearity));
      return;

  }

  if(!fPeriodNameAvailable){
    AliFatal(Form("NonLinearity correction not defined for fPeriodName: '%s'! Please check cut number (%d) as well as function AliCaloPhotonCuts::ApplyNonLinearity. Correction failed, returning...",fPeriodName.Data(),fSwitchNonLinearity));
    return;
  }

  cluster->SetE(energy);

  return;
}

//________________________________________________________________________
void AliCaloPhotonCuts::ApplySMWiseEnergyCorrection(AliVCluster* cluster, Int_t isMC, AliVEvent *event)
{
  if (!cluster) {
    AliInfo("Cluster pointer null!");
    return;
  }
  Float_t energy = cluster->E();

  Int_t clusterSMID = -1;
  if(fClusterType == 4 && event){
    Int_t largestCellIDcluster = FindLargestCellInCluster(cluster,event);
    if(largestCellIDcluster>-1){
      Int_t dummycol = -1, dummyrow = -1;
      clusterSMID = GetModuleNumberAndCellPosition(largestCellIDcluster, dummycol, dummyrow);
    }
  }
  if( fClusterType == 1 || fClusterType == 3|| fClusterType == 4){
    if (energy < 0.05) {
      // Clusters with less than 50 MeV or negative are not possible
      AliInfo(Form("Too Low Cluster energy!, E = %f < 0.05 GeV",energy));
      return;
    }
  } else {
    if (energy < 0.01) {
      // Clusters with less than 10 MeV or negative are not possible
      AliInfo(Form("Too Low Cluster energy!, E = %f < 0.01 GeV",energy));
      return;
    }
  }
  if(fCurrentMC==kNoMC){
    AliV0ReaderV1* V0Reader = (AliV0ReaderV1*) AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
    if( V0Reader == NULL ){
      AliFatal(Form("No V0Reader called '%s' could be found within AliCaloPhotonCuts::ApplySMWiseEnergyCorrection",fV0ReaderName.Data()));
      return;
    }
    fPeriodName = V0Reader->GetPeriodName();
    fCurrentMC = FindEnumForMCSet(fPeriodName);

    printf("AliCaloPhotonCuts:Period name has been set to %s, period-enum: %o\n",fPeriodName.Data(),fCurrentMC ) ;
  }


  if(isMC == 0){         // data; SM wise correction
    if( fCurrentMC == k16pp13TeV || fCurrentMC == k17pp13TeV || fCurrentMC == k18pp13TeV ){
      switch (clusterSMID){
        // values determined on LHC16x & LHC17c
        case 0: energy/=0.994364; break;
        case 1: energy/=0.991352; break;
        case 2: energy/=1.000522; break;
        case 3: energy/=0.995918; break;
        case 4: energy/=0.995661; break;
        case 5: energy/=0.998285; break;
        case 6: energy/=1.000275; break;
        case 7: energy/=1.003544; break;
        case 8: energy/=1.007220; break;
        case 9: energy/=1.000911; break;
        case 10: energy/=1.012508; break;
        case 11: energy/=1.012867; break;
        case 12: energy/=1.001028; break;
        case 13: energy/=0.995514; break;
        case 14: energy/=0.994373; break;
        case 15: energy/=0.997765; break;
        case 16: energy/=1.009084; break;
        case 17: energy/=1.011123; break;
        case 18: energy/=1.007118; break;
        case 19: energy/=1.018894; break;
        default: energy/=1.0; break;
      }
    } else if(fCurrentMC == k17pp5TeV){
      switch (clusterSMID){
        // values determined on LHC17pq
        case 0: energy/=0.996406; break;
        case 1: energy/=0.995672; break;
        case 2: energy/=1.000260; break;
        case 3: energy/=0.998378; break;
        case 4: energy/=0.998326; break;
        case 5: energy/=0.999704; break;
        case 6: energy/=1.001140; break;
        case 7: energy/=0.999142; break;
        case 8: energy/=1.002450; break;
        case 9: energy/=1.002400; break;
        case 10: energy/=1.006850; break;
        case 11: energy/=1.005900; break;
        case 12: energy/=0.999443; break;
        case 13: energy/=0.996047; break;
        case 14: energy/=0.997640; break;
        case 15: energy/=0.996848; break;
        case 16: energy/=1.005100; break;
        case 17: energy/=1.007150; break;
        case 18: energy/=1.006660; break;
        case 19: energy/=1.008700; break;
        default: energy/=1.0; break;
      }
    }
  }

  cluster->SetE(energy);

  return;
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kPi0MC(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5, Float_t p6){
  return ( p6 / ( p0 * ( 1. / ( 1. + p1 * exp( -e / p2 ) ) * 1. / ( 1. + p3 * exp( ( e - p4 ) / p5 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kSDM(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3){
  return ( p0 + p3 * exp( p1 + ( p2 * e ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_DPOW(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5){
  Float_t ret = 1;
  if ((p3 +  p4 * TMath::Power(e,p5 ) ) != 0)
    ret = ( (p0 +  p1 * TMath::Power(e,p2 ) )/(p3 +  p4 * TMath::Power(e,p5 ) ) );
  if (ret != 0.)
    return ret;
  else
    return 1.;
}
//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_SPOW(Float_t e, Float_t p0, Float_t p1, Float_t p2){
  Float_t ret = 1;
  if ((p0 -  p1 * TMath::Power(e,p2 ) ) != 0)
    ret = (p0 -  p1 * TMath::Power(e,p2 ) );
  if (ret != 0.)
    return ret;
  else
    return 1.;
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_DExp(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5, Float_t p6, Float_t p7){
  Float_t ret = 1;
  if ( (p3 - p7*TMath::Exp(-p4*e+p5) ) != 0)
    ret = ( (p0 - p6*TMath::Exp(-p1*e+p2) )/(p3 - p7*TMath::Exp(-p4*e+p5) ) );
  if (ret != 0.)
    return ret;
  else
    return 1.;
}
Float_t AliCaloPhotonCuts::FunctionNL_SExp(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3){
  return ( p0 - p3 * TMath::Exp( - p1 * e + p2 ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_PHOSOnlyMC(Float_t e, Float_t p0, Float_t p1, Float_t p2){
  return p0*(1+p1/(1.+e*e/p2/p2)) ;
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_ExpExp(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3){
    // "[0] - TMath::Exp(-[1]*x+[2]) + TMath::Exp(-[3]*x)";
    Float_t ret = ( p0 - TMath::Exp(-p1*e+p2) + TMath::Exp(-p3*e));
    if (ret != 0.)
      return ret;
    else
      return 1.;
}


//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_LinLogConst(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t const1, Float_t const2){
    //Function splitted into multiple parts:
    //Constant Correction, whole energery range
    //Linear Function, from 0 to p0
    //Logarithmic Function from p0 to p1
    //Constant from p1
    //    p2+((e<p0)*(p3*e))+((e>=p0)*(e<p1)*( (p3*p0)+(TMath::Log(p4*(e-p0)+1))))+((e>=p1)*((p3*p0)+(TMath::Log(p4*(p1-p0)+1))));
    // => p2+((e<p0)*(p3*e))+((e>=p0)*(e<p1)*((const1)+(TMath::Log(p4*(e-p0)+1))))+((e>=p1)*const2);
    Float_t ret=p2;
    if (e<p0){
        ret+=e*p3;
    } else if ((e>=p0)&&(e<p1)) {
        ret+=const1+TMath::Log(p4*(e-p0)+1);
    } else {
        ret+=const2;
    }
    if (ret != 0.)
      return ret;
    else
      return 1.;
}


//************************************************************************
// predefined functions:
//________________________________________________________________________
// testbeam parametrizations by Martin and Nico for different aggregation thresholds (still work in progress)
// Data params
Float_t AliCaloPhotonCuts::FunctionNL_OfficialTB_50MeV_Data(Float_t e){
  Double_t funcParams[5] = {0.960211, 0.0142135, 0.0786752, 130.306, 65.9035};
  return ( 1.0585 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
}
Float_t AliCaloPhotonCuts::FunctionNL_OfficialTB_100MeV_Data(Float_t e){
  Double_t funcParams[5] = {0.941138, 0.0172153, 0.0783153, 130.869, 64.9742};
  return ( 1.0505 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
}
Float_t AliCaloPhotonCuts::FunctionNL_OfficialTB_100MeV_Data_V2(Float_t e){
  Double_t funcParams[5] = {1.91897, 0.0264988, 0.965663, -187.501, 2762.51};
  return ( 1.0505 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
}
Float_t AliCaloPhotonCuts::FunctionNL_OfficialTB_100MeV_Data_Sys0(Float_t e){ // refit on systematic 1sigma shift upwards of all TB points
  Double_t funcParams[5] = {2.864, 0.031267, 1.89089, -425.59, 8525.01};
  return ( 1.0505 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
}
Float_t AliCaloPhotonCuts::FunctionNL_OfficialTB_100MeV_Data_Sys1(Float_t e){ // refit on systematic 1sigma shift downwards of all TB points
  Double_t funcParams[5] = {0.98992, 0.015482, 0.0993751, 126.234, 293.826};
  return ( 1.0505 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
}
Float_t AliCaloPhotonCuts::FunctionNL_OfficialTB_100MeV_Data_Sys2(Float_t e){ // refit on systematic 1sigma tilt shift (-1sigma to +1sigma) of TB points
  Double_t funcParams[5] = {3.3084, 0.0561285, 2.16752, -822.236, 5185.83};
  return ( 1.0505 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
}
Float_t AliCaloPhotonCuts::FunctionNL_OfficialTB_100MeV_Data_Sys3(Float_t e){ // refit on systematic 1sigma tilt shift (+1sigma to -1sigma) of TB points
  Double_t funcParams[5] = {0.96096, 0.00786329, 0.0299616, 151.315, 94.2094};
  return ( 1.0505 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
}
Float_t AliCaloPhotonCuts::FunctionNL_OfficialTB_150MeV_Data(Float_t e){
  Double_t funcParams[5] = {0.921363, 0.0200311, 0.0776928, 132.598, 62.9008};
  return ( 1.0505 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
}
Float_t AliCaloPhotonCuts::FunctionNL_OfficialTB_300MeV_Data(Float_t e){
  Double_t funcParams[5] = {0.88448, 0.0240087, 0.0712406, 136.93, 55.1195};
  return ( 1.0505 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
}
// MC params
Float_t AliCaloPhotonCuts::FunctionNL_OfficialTB_50MeV_MC(Float_t e){
  Double_t funcParams[5] = {4.24777, 0.0383424, 3.00719, -536.41, 4329.67};
  return ( 1.00 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
}
Float_t AliCaloPhotonCuts::FunctionNL_OfficialTB_100MeV_MC(Float_t e){
  Double_t funcParams[5] = {4.37267, 0.0636557, 3.11756, -613.126, 3913.14};
  return ( 1.00 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
}
Float_t AliCaloPhotonCuts::FunctionNL_OfficialTB_100MeV_MC_V2(Float_t e){
  Double_t funcParams[5] = {1.09357, 0.0192266, 0.291993, 370.927, 694.656};
  return ( 1.00 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
}
Float_t AliCaloPhotonCuts::FunctionNL_OfficialTB_150MeV_MC(Float_t e){
  Double_t funcParams[5] = {4.68642, 0.0844255, 3.49109, -572.864, 3786.48};
  return ( 1.00 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
}
Float_t AliCaloPhotonCuts::FunctionNL_OfficialTB_300MeV_MC(Float_t e){
  Double_t funcParams[5] = {3.13707, 0.0675494, 2.16932, -483.678, 3730.81};
  return ( 1.00 * (funcParams[0] + funcParams[1] * TMath::Log(e) ) / ( 1 + ( funcParams[2] * TMath::Exp( ( e - funcParams[3] ) / funcParams[4] ) ) ) );
}

// other testbeam parametrization
Float_t AliCaloPhotonCuts::FunctionNL_kPi0MCv1(Float_t e){
  return ( 1.014 * exp( 0.03329 / e ) ) + ( ( -0.3853 / ( 0.5423 * 2. * TMath::Pi() ) * exp( -( e + 0.4335 ) * ( e + 0.4335 ) / (2. * 0.5423 * 0.5423 ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kPi0MCv2(Float_t e){
  return ( 0.311111 / TMath::Power( e - 0.571666, 0.567995 ) + 1 );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kPi0MCv3(Float_t e){
  return ( 1.0 / ( 0.981039 * ( 1. / ( 1. + 0.113508 * exp( -e / 1.00173 ) ) * 1. / ( 1. + 0.0967998 * exp( ( e - 219.381 ) / 63.1604 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kPi0MCv5(Float_t e){
  return ( 1.01286 / ( 1.0 * ( 1. / ( 1. + 0.0664778 * exp( -e / 1.57 ) ) * 1. / ( 1. + 0.0967998 * exp( ( e - 219.381 ) / 63.1604 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kPi0MCv6(Float_t e){
  return ( 1.00437 / ( 1.0 * ( 1. / ( 1. + 0.0797873 * exp( -e / 1.68322 ) ) * 1. / ( 1. + 0.0806098 * exp( ( e - 244.586 ) / 116.938 ) ) ) ) );
}

// only shifting data, to be used with kPi0MCv5 before
//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kSDMv5(Float_t e){
  return ( 0.964 + exp( -3.132 + ( -0.435 * 2.0 * e ) ) );
}

// be careful: different definition than kSDMv5
//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kSDMv6(Float_t e){
  return ( 0.987054 / ( 1.0 * ( 1. / ( 1. + 0.237767 * exp( -e / 0.651203 ) ) * 1. / ( 1. + 0.183741 * exp( ( e - 155.427 ) / 17.0335 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kTestBeamv2(Float_t e){
  return ( 0.968 / ( 0.983504 *( 1. / ( 1. + 0.210106 * exp( -e / 0.897274 ) ) * 1. / ( 1. + 0.0829064 * exp( ( e - 152.299 ) / 31.5028 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kTestBeamv3(Float_t e){
  return ( 0.9615 / ( 0.976941 *( 1. / ( 1. + 0.162310 * exp( -e / 1.08689 ) ) * 1. / ( 1. + 0.0819592 * exp( ( e - 152.338 ) / 30.9594 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kTestBeamv4(Float_t e){ // supplied by Evi Ganoti, cf. https://alice.its.cern.ch/jira/browse/EMCAL-190
  return ( 0.97 / ( 0.9892 *( 1. / ( 1. + 0.1976 * exp( -e / 0.865 ) ) * 1. / ( 1. + 0.06775 * exp( ( e - 156.6 ) / 47.18 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kTestBeamMod(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5, Float_t p6){
  return ( p0 / ( p1 *( 1. / ( 1. + p2 * exp( -e / p3 ) ) * 1. / ( 1. + p4 * exp( ( e - p5 ) / p6 ) ) ) ) );
}

//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionNL_kPi0MCMod(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5, Float_t p6){
  return ( p0 / ( p1 *( 1. / ( 1. + p2 * exp( -e / p3 ) ) * 1. / ( 1. + p4 * exp( ( e - p5 ) / p6 ) ) ) ) );
}
//************************************************************************
//________________________________________________________________________
Float_t AliCaloPhotonCuts::FunctionM02(Float_t E, Float_t a, Float_t b, Float_t c, Float_t d, Float_t e){
  return ( exp( a+ b*E ) + c + d*E + e/E);
}

//________________________________________________________________________
AliCaloPhotonCuts::MCSet AliCaloPhotonCuts::FindEnumForMCSet(TString namePeriod){
  // pp 7TeV MB MCs pass 4
  if(       namePeriod.CompareTo("LHC14j4b")==0 ||
            namePeriod.CompareTo("LHC14j4c")==0 ||
            namePeriod.CompareTo("LHC14j4d")==0 ||
            namePeriod.CompareTo("LHC14j4e")==0 ||
            namePeriod.CompareTo("LHC14j4f")==0)        return k14j4;

  // pp 7 TeV LHC11 anch MC pass 1
  else if(  namePeriod.CompareTo("LHC14b7")==0)        return k14b7;
  else if(  namePeriod.CompareTo("LHC14k1a")==0)        return k14k1ab;
  else if(  namePeriod.CompareTo("LHC14k1b")==0)        return k14k1ab;
  // pp 2.76 TeV LHC11a anch MC's pass 4
  else if(  namePeriod.CompareTo("LHC12f1a")==0)        return k12f1a;
  else if(  namePeriod.CompareTo("LHC12f1b")==0)        return k12f1b;
  else if(  namePeriod.CompareTo("LHC12i3")==0)         return k12i3;
  else if(  namePeriod.CompareTo("LHC15g1a")==0)        return kPP2T11P4JJ;
  else if(  namePeriod.CompareTo("LHC15g1b")==0)        return k15g1b;

  // PbPb 2.76TeV 2011 MC add sig
  else if(  namePeriod.CompareTo("LHC14a1a")==0 ||
            namePeriod.CompareTo("LHC14a1b")==0 ||
            namePeriod.CompareTo("LHC14a1c")==0)        return k14a1;

    // pp 8 TeV MC MCs
  // pass 1
  else if(  namePeriod.CompareTo("LHC14e2b")==0)        return k14e2b;
  // pass 2
  else if(  namePeriod.CompareTo("LHC12P2Pyt8")==0 ||
            namePeriod.CompareTo("LHC15h1")==0 ||
            namePeriod.CompareTo("LHC15h1a1")==0 ||
            namePeriod.CompareTo("LHC15h1b")==0 ||
            namePeriod.CompareTo("LHC15h1c")==0 ||
            namePeriod.CompareTo("LHC15h1d")==0 ||
            namePeriod.CompareTo("LHC15h1f")==0 ||
            namePeriod.CompareTo("LHC15h1g")==0 ||
            namePeriod.CompareTo("LHC15h1h")==0 ||
            namePeriod.CompareTo("LHC15h1i")==0)        return kPP8T12P2Pyt8;
  else if(  namePeriod.CompareTo("LHC12P2Pho")==0 ||
            namePeriod.CompareTo("LHC15h2")==0 ||
            namePeriod.CompareTo("LHC15h2a")==0 ||
            namePeriod.CompareTo("LHC15h2b")==0 ||
            namePeriod.CompareTo("LHC15h2c")==0 ||
            namePeriod.CompareTo("LHC15h2d")==0 ||
            namePeriod.CompareTo("LHC15h2f")==0 ||
            namePeriod.CompareTo("LHC15h2g")==0 ||
            namePeriod.CompareTo("LHC15h2h")==0 ||
            namePeriod.CompareTo("LHC15h2i")==0)        return kPP8T12P2Pho;
  // pp 8 TeV JJ MC pass 2
  else if ( namePeriod.CompareTo("LHC12P2JJ") == 0 ||
            namePeriod.CompareTo("LHC16c2") == 0 ||
            namePeriod.CompareTo("LHC16c2_plus") == 0 ) return kPP8T12P2JJ;
  else if ( namePeriod.CompareTo("LHC17g5b") == 0 ) return kPP8T12P2GJLow;
  else if ( namePeriod.CompareTo("LHC17g5c") == 0 ) return kPP8T12P2GJHigh;

  // pPb 5 TeV 2013 MC pass 2
  else if(  namePeriod.Contains("LHC13b2_efix"))        return kPPb5T13P2DPMJet;
  else if(  namePeriod.Contains("LHC13b4_fix") ||
            namePeriod.Contains("LHC13b4_plus"))        return kPPb5T13P4JJ;
  else if(  namePeriod.Contains("LHC13e7"))             return kPPb5T13P2HIJAdd;
  else if(  namePeriod.Contains("LHC18j5"))             return kPPb5T13P4DPMJet;
  else if(  namePeriod.Contains("LHC19a4"))             return kPPb5T13P4JJ;
  // pPb 5 TeV 2013 MC JJ EMC enhanced
  else if ( namePeriod.CompareTo("LHC16c3a") == 0 ||
            namePeriod.CompareTo("LHC16c3a2") == 0 )    return k16c3a;
  else if ( namePeriod.CompareTo("LHC16c3b") == 0 ||
            namePeriod.CompareTo("LHC16c3b2") == 0 )    return k16c3b;
  else if ( namePeriod.CompareTo("LHC17g6a2") == 0 )    return kPPb5T13P4JJlow;
  else if ( namePeriod.CompareTo("LHC17g6a3") == 0 )    return kPPb5T13P4JJhigh;
  // pPb 5 TeV 2013 MC GJ
  else if ( namePeriod.CompareTo("LHC16c3c") == 0 ||
            namePeriod.CompareTo("LHC16c3c2") == 0 )    return k16c3c;

  // pp 2.76 TeV LHC13g anch MC's pass 4
  else if(  namePeriod.CompareTo("LHC15g2")==0)         return k15g2;
  else if(  namePeriod.CompareTo("LHC15a3a")==0 ||
            namePeriod.CompareTo("LHC15a3a_plus")==0)   return kPP2T13P1JJ;
  else if(  namePeriod.CompareTo("LHC15a3b")==0)        return k15a3b;

  // pp 13 TeV 2015 MB pass 2
  else if ( namePeriod.CompareTo("LHC15P2Pyt8") == 0 ||
            namePeriod.CompareTo("LHC17i4") == 0 ||
            namePeriod.CompareTo("LHC17i4_2") == 0 ||
            namePeriod.CompareTo("LHC17g7") == 0 )      return kPP13T15P2Pyt8;
  else if ( namePeriod.CompareTo("LHC15P2EPos") == 0 ||
            namePeriod.CompareTo("LHC16d3") == 0 )      return kPP13T15P2EPOS;
  // pp 13 TeV 2015 HF prod
  else if ( namePeriod.CompareTo("LHC15k5a") == 0 ||
            namePeriod.CompareTo("LHC15k5b") == 0 ||
            namePeriod.CompareTo("LHC15k5c") == 0 ||
            namePeriod.CompareTo("LHC15k5a2") == 0 ||
            namePeriod.CompareTo("LHC15k5b2") == 0 ||
            namePeriod.CompareTo("LHC15k5c2") == 0 )    return k15k5;

  // pp 5 TeV 2015 MB MC pass 2
  else if ( namePeriod.CompareTo("LHC16h8a") == 0 )     return k16h8a;
  else if ( namePeriod.CompareTo("LHC16h8b") == 0 )     return k16h8b;
  else if ( namePeriod.CompareTo("LHC16k3a") == 0 ||                    // special pileup prods
            namePeriod.CompareTo("LHC16k3a2") == 0 )    return k16k3a;
  // pp 5 TeV 2015 MB MC pass 3
  else if ( namePeriod.CompareTo("LHC16k5a") == 0  )    return k16k5a;
  else if ( namePeriod.CompareTo("LHC16k5b") == 0  )    return k16k5b;
  // pp 5 TeV 2015 MB MC pass 4
  else if ( namePeriod.CompareTo("LHC17e2") == 0  )     return k17e2;
  else if ( namePeriod.CompareTo("LHC18j3") == 0  )     return k18j3;
  // pp 5 TeV 2015 JJ pass 3
  else if ( namePeriod.CompareTo("LHC16h3") == 0 )      return k16h3;
  // pp 5 TeV 2017 MB MC pass 1
  else if ( namePeriod.CompareTo("LHC17l4b") == 0 ||
            namePeriod.CompareTo("LHC17l4b_fast") == 0 ||
            namePeriod.CompareTo("LHC17l4b_cent") == 0 ||
            namePeriod.CompareTo("LHC17l4b_cent_woSDD") == 0)     return k17l4b;
  else if ( namePeriod.CompareTo("LHC17l3b") == 0 ||
            namePeriod.CompareTo("LHC17l3b_fast") == 0 ||
            namePeriod.CompareTo("LHC17l3b_cent") == 0 ||
            namePeriod.CompareTo("LHC17l3b_cent_woSDD") == 0)     return k17l3b;
  else if ( namePeriod.CompareTo("LHC18j2") == 0 ||
            namePeriod.CompareTo("LHC18j2_fast") == 0 ||
            namePeriod.CompareTo("LHC18j2_cent") == 0 ||
            namePeriod.CompareTo("LHC18j2_cent_woSDD") == 0)     return k18j2;
  // pp 5 TeV 2015 JJ pass 4
  else if ( namePeriod.Contains("LHC18b8") )      return k18b8;
  // pp 5 TeV 2015 GJ pass 4
  else if ( namePeriod.Contains("LHC18b10") )      return k18b10;
  else if ( namePeriod.Contains("LHC18l2") ||
            namePeriod.Contains("LHC18g7") )      return k18l2;
  // PbPb 5TeV 2015 MB prods
  else if ( namePeriod.Contains("LHC18e1") || namePeriod.CompareTo("LHC16h4") == 0 )    return kPbPb5T15HIJING;
  else if ( namePeriod.CompareTo("LHC16k3b") == 0 ||                    // special pileup prods
            namePeriod.CompareTo("LHC16k3b2") == 0 )    return k16k3b;
  // PbPb 5TeV 2018 MB prods
  else if ( namePeriod.CompareTo("LHC18l8a") == 0 ||
            namePeriod.CompareTo("LHC18l8b") == 0 ||
            namePeriod.CompareTo("LHC18l8c") == 0 ||
            namePeriod.CompareTo("LHC19h2a") == 0 ||
            namePeriod.CompareTo("LHC19h2b") == 0 ||
            namePeriod.CompareTo("LHC19h2c") == 0 ||
            namePeriod.CompareTo("LHC19h3")  == 0)   return kPbPb5T18HIJING;

  // pp 13 TeV 2016 MB prod
  else if ( namePeriod.CompareTo("LHC16P1Pyt8") == 0 ||
            namePeriod.CompareTo("LHC16P1Pyt8NomB") ==0 ||
            namePeriod.CompareTo("LHC17f6") == 0 ||
            namePeriod.CompareTo("LHC17d17") == 0 ||
            namePeriod.CompareTo("LHC17f5") == 0 ||
            namePeriod.CompareTo("LHC17d3") == 0 ||
            namePeriod.CompareTo("LHC17e5") == 0 ||
            namePeriod.CompareTo("LHC17d20a1") == 0 ||
            namePeriod.CompareTo("LHC17d20a1_extra") == 0 ||
            namePeriod.CompareTo("LHC17d20a2") == 0 ||
            namePeriod.CompareTo("LHC17d20a2_extra") == 0 ||
            namePeriod.CompareTo("LHC17d16") == 0 ||
            namePeriod.CompareTo("LHC17d18") == 0 ||
            namePeriod.CompareTo("LHC17f9") == 0 ||
            namePeriod.CompareTo("LHC17f9_test") == 0||
            namePeriod.CompareTo("LHC18f1") == 0 ||
            namePeriod.CompareTo("LHC18d8") == 0 ||
            namePeriod.CompareTo("LHC19g7a") == 0 )     return kPP13T16P1Pyt8;
  else if ( namePeriod.CompareTo("LHC16P1Pyt8LowB") == 0 ||
            namePeriod.CompareTo("LHC17d1") == 0 )      return kPP13T16P1Pyt8LowB;
  else if ( namePeriod.CompareTo("LHC16P1EPOS") == 0 ||
            namePeriod.CompareTo("LHC17d20b1") == 0 ||
            namePeriod.CompareTo("LHC17d20b2") == 0 )   return kPP13T16P1EPOS;
  // pp 13 TeV 2016 JJ prod
  else if ( namePeriod.CompareTo("LHC16P1JJ") == 0 ||
            namePeriod.CompareTo("LHC17f8a") == 0 ||
            namePeriod.CompareTo("LHC17f8c") == 0 ||
            namePeriod.CompareTo("LHC17f8d") == 0 ||
            namePeriod.CompareTo("LHC17f8e") == 0 ||
            namePeriod.CompareTo("LHC17f8f") == 0 ||
            namePeriod.CompareTo("LHC17f8g") == 0 ||
            namePeriod.CompareTo("LHC17f8h") == 0 ||
            namePeriod.CompareTo("LHC17f8i") == 0 ||
            namePeriod.CompareTo("LHC17f8j") == 0 ||
            namePeriod.CompareTo("LHC17f8k") == 0 )     return kPP13T16P1JJ;
  else if ( namePeriod.CompareTo("LHC16P1JJLowB") == 0 ||
            namePeriod.CompareTo("LHC17f8b") == 0  )    return kPP13T16P1JJLowB;
  //pp 13 TeV LHC16 JJ MC with high enrgy gamma in EMC, DMC, PHOS acc.
  else if ( namePeriod.CompareTo("LHC17i3b1") == 0 ||
            namePeriod.CompareTo("LHC17i3b2") == 0 ||
            namePeriod.CompareTo("LHC17i3c1") == 0 ||
            namePeriod.CompareTo("LHC17i3c2") == 0 ||
            namePeriod.CompareTo("LHC20b1b1") == 0 ||
            namePeriod.CompareTo("LHC20b1b2") == 0 ||
            namePeriod.CompareTo("LHC20b1c1") == 0 ||
            namePeriod.CompareTo("LHC20b1c2") == 0 ) return kPP13T16P1JJTrigger;
  // pp 13 TeV 2016 HF prods
  else if ( namePeriod.CompareTo("LHC17h8a") == 0 )     return k17h8a;
  else if ( namePeriod.CompareTo("LHC17h8b") == 0 )     return k17h8b;
  else if ( namePeriod.CompareTo("LHC17h8c") == 0 )     return k17h8c;
  else if ( namePeriod.CompareTo("LHC17c3b1") == 0 )    return k17c3b1;
  else if ( namePeriod.CompareTo("LHC17c3a1") == 0 )    return k17c3a1;
  else if ( namePeriod.CompareTo("LHC17c3b2") == 0 )    return k17c3b2;
  else if ( namePeriod.CompareTo("LHC17c3a2") == 0 )    return k17c3a2;
  // pp 13 TeV 2016 GJ prod
  else if ( namePeriod.CompareTo("LHC17i3a1") == 0 )    return k17i3a1;

  // pPb 5 TeV 2016 MB MC
  else if ( namePeriod.CompareTo("LHC17f2a") == 0  )    return kPPb5T16EPOS;
  else if ( namePeriod.CompareTo("LHC17f2b") == 0 ||
            namePeriod.CompareTo("LHC18f3")  == 0 )     return kPPb5T16DPMJet;
  // pPb 5 TeV 2016 JJ MC
  else if ( namePeriod.CompareTo("LHC17g8a") == 0  )    return k17g8a;
  // pPb 5 TeV 2016 HF prod
  else if ( namePeriod.CompareTo("LHC17d2a") == 0 )     return k17d2a;
  else if ( namePeriod.CompareTo("LHC17d2b") == 0 )     return k17d2b;

  // pPb 8 TeV 2016 MB MC
  else if ( namePeriod.CompareTo("LHC17f3a") == 0  )    return k17f3a;
  else if ( namePeriod.CompareTo("LHC17f3b") == 0  )    return k17f3b;
  else if ( namePeriod.CompareTo("LHC17f4a") == 0  )    return k17f4a;
  else if ( namePeriod.CompareTo("LHC17f4b") == 0  )    return k17f4b;
  else if ( namePeriod.Contains("LHC18f3b") ||
            namePeriod.Contains("LHC18f3c") )           return k18f3bc;
  // pPb 8 TeV 2016 decay gamma MC
  else if ( namePeriod.CompareTo("LHC17g6b2a") == 0  )    return k17g6b2a;
  else if ( namePeriod.CompareTo("LHC17g6b2b") == 0  )    return k17g6b2b;
  else if ( namePeriod.CompareTo("LHC17g6b3a") == 0  )    return k17g6b3a;
  else if ( namePeriod.CompareTo("LHC17g6b3b") == 0  )    return k17g6b3b;
  // pPb 8 TeV 2016 JJ MC
  else if ( namePeriod.CompareTo("LHC17g8b") == 0  )    return k17g8b;
  else if ( namePeriod.CompareTo("LHC17g8c") == 0  )    return k17g8c;
  else if ( namePeriod.CompareTo("LHC18b9b") == 0  )    return k18b9b;
  else if ( namePeriod.CompareTo("LHC18b9c") == 0  )    return k18b9c;

  //pp 13 TeV LHC17
  else if ( namePeriod.CompareTo("LHC17k1") ==0 )       return k17k1; // HF low B
  else if ( namePeriod.CompareTo("LHC17P1Pyt8NomB") ==0 ||
            namePeriod.CompareTo("LHC17k4") ==0 ||
            namePeriod.CompareTo("LHC17h11") ==0 ||
            namePeriod.CompareTo("LHC17h1") == 0 ||
            namePeriod.CompareTo("LHC17l5") == 0 ||
            namePeriod.CompareTo("LHC18c13") == 0 ||
            namePeriod.CompareTo("LHC18a8") == 0 ||
            namePeriod.CompareTo("LHC18a9") == 0 ||
            namePeriod.CompareTo("LHC18a1") == 0 ||
            namePeriod.CompareTo("LHC18c12") == 0 ||
            namePeriod.CompareTo("LHC17P1Pyt8") == 0 ||
            namePeriod.CompareTo("LHC19g7b") == 0 )     return kPP13T17P1Pyt8;
  else if ( namePeriod.CompareTo("LHC17h7b") ==0 )      return kPP13T17P1Pho;
  else if ( namePeriod.CompareTo("LHC17h7a") ==0 )      return kPP13T17P1Pyt6;

  else if ( namePeriod.CompareTo("LHC17j5a") ==0 ||
            namePeriod.CompareTo("LHC17j5b") ==0 ||
            namePeriod.CompareTo("LHC17j5c") ==0 ||
            namePeriod.CompareTo("LHC17j5d") ==0 ||
            namePeriod.CompareTo("LHC17j5e") ==0 )      return kPP13T17P1Pyt8Str;
  else if ( namePeriod.CompareTo("LHC17P1Pyt8LowB") ==0 ||
            namePeriod.CompareTo("LHC17h3") == 0 )      return kPP13T17P1Pyt8LowB;
  else if ( namePeriod.CompareTo("LHC17P1JJ") == 0 ||
            namePeriod.CompareTo("LHC18f5") == 0 )      return kPP13T17P1JJ;
  else if ( namePeriod.CompareTo("LHC18l6b1") == 0 ||
            namePeriod.CompareTo("LHC18l6b2") == 0 ||
            namePeriod.CompareTo("LHC18l6c1") == 0 ||
            namePeriod.CompareTo("LHC18l6c2") == 0 )    return kPP13T17P1JJTrigger;
  // XeXe 5.44 TeV 2017 MB MC
  else if ( namePeriod.CompareTo("LHC17j7") == 0 ||
            namePeriod.CompareTo("LHC18d2") == 0 ||
            namePeriod.CompareTo("LHC18d2_1") == 0 ||
            namePeriod.CompareTo("LHC18d2_2") == 0 ||
            namePeriod.CompareTo("LHC18d2_3") == 0 )    return kXeXe5T17HIJING;
  //pp 13 TeV LHC18
  else if ( namePeriod.CompareTo("LHC18P1Pyt8NomB") ==0 ||
            namePeriod.CompareTo("LHC18P1Pyt8") ==0 ||
            namePeriod.CompareTo("LHC18g4") ==0 ||
            namePeriod.CompareTo("LHC18g5") ==0 ||
            namePeriod.CompareTo("LHC18g6") == 0 ||
            namePeriod.CompareTo("LHC19g7c")==0 )      return kPP13T18P1Pyt8;
  else if ( namePeriod.CompareTo("LHC18P1Pyt8LowB") ==0 ||
            namePeriod.CompareTo("LHC18h1") ==0  )      return kPP13T18P1Pyt8LowB;
  //pp 13 TeV LHC18 JJ MCs
  else if ( namePeriod.CompareTo("LHC18P1JJ") == 0 ||
            namePeriod.Contains("LHC19d3") )            return kPP13T18P1JJ;
  //pp 13 TeV LHC18 JJ MC with high enrgy gamma in EMC, DMC, PHOS acc.
  else if ( namePeriod.CompareTo("LHC19i3b1") == 0 ||
            namePeriod.CompareTo("LHC19i3b2") == 0 ||
            namePeriod.CompareTo("LHC19i3c1") == 0 ||
            namePeriod.CompareTo("LHC19i3c2") == 0 ) return kPP13T18P1JJTrigger;
  // PbPb 5 TeV 2015 Gamma-Jet MC
  else if ( namePeriod.CompareTo("LHC18b11c") == 0 ) return  kLHC18b11c;

  // data starts here
  else if ( namePeriod.CompareTo("LHC10b") == 0 ||
            namePeriod.CompareTo("LHC10c") == 0 ||
            namePeriod.CompareTo("LHC10d") == 0 ||
            namePeriod.CompareTo("LHC10e") == 0 ||
            namePeriod.CompareTo("LHC10f") == 0 ||
            namePeriod.CompareTo("LHC10g") == 0 )       return k10pp7TeV;
  else if ( namePeriod.CompareTo("LHC10h") == 0 )       return k10PbPb2760GeV;
  else if ( namePeriod.CompareTo("LHC11a") == 0 )       return k11pp2760GeV;
  else if ( namePeriod.CompareTo("LHC11b") == 0 ||
            namePeriod.CompareTo("LHC11c") == 0 ||
            namePeriod.CompareTo("LHC11d") == 0 ||
            namePeriod.CompareTo("LHC11e") == 0 ||
            namePeriod.CompareTo("LHC11f") == 0 ||
            namePeriod.CompareTo("LHC11g") == 0 )       return k11pp7TeV;
  else if ( namePeriod.CompareTo("LHC11h") == 0 )       return k11PbPb2760GeV;
  else if ( namePeriod.CompareTo("LHC12a") == 0 ||
            namePeriod.CompareTo("LHC12b") == 0 ||
            namePeriod.CompareTo("LHC12c") == 0 ||
            namePeriod.CompareTo("LHC12d") == 0 ||
            namePeriod.CompareTo("LHC12e") == 0 ||
            namePeriod.CompareTo("LHC12f") == 0 ||
            namePeriod.CompareTo("LHC12g") == 0 ||
            namePeriod.CompareTo("LHC12h") == 0 ||
            namePeriod.CompareTo("LHC12i") == 0 )       return k12pp8TeV;
  else if ( namePeriod.CompareTo("LHC13b") == 0 ||
            namePeriod.CompareTo("LHC13c") == 0 ||
            namePeriod.CompareTo("LHC13d") == 0 ||
            namePeriod.CompareTo("LHC13e") == 0 ||
            namePeriod.CompareTo("LHC13f") == 0 )       return k13pPb5023GeV;
  else if ( namePeriod.CompareTo("LHC13g") == 0 )       return k13pp2760GeV;
  else if ( namePeriod.CompareTo("LHC15fm") == 0 ||
            namePeriod.CompareTo("LHC15f") == 0 ||
            namePeriod.CompareTo("LHC15g") == 0 ||
            namePeriod.CompareTo("LHC15h") == 0 ||
            namePeriod.CompareTo("LHC15i") == 0 ||
            namePeriod.CompareTo("LHC15j") == 0 ||
            namePeriod.CompareTo("LHC15k") == 0 ||
            namePeriod.CompareTo("LHC15l") == 0 ||
            namePeriod.CompareTo("LHC15m") == 0 )       return k15pp13TeV;
  else if ( namePeriod.CompareTo("LHC15n") == 0 )       return k15pp5TeV;
  else if ( namePeriod.CompareTo("LHC15o") == 0 )       return k15PbPb5TeV;
  else if ( namePeriod.CompareTo("LHC16f") == 0 )       return k16pp13TeVLow;
  else if ( namePeriod.CompareTo("LHC16dp") == 0 ||
            namePeriod.CompareTo("LHC16d") == 0 ||
            namePeriod.CompareTo("LHC16e") == 0 ||
            namePeriod.CompareTo("LHC16g") == 0 ||
            namePeriod.CompareTo("LHC16h") == 0 ||
            namePeriod.CompareTo("LHC16i") == 0 ||
            namePeriod.CompareTo("LHC16j") == 0 ||
            namePeriod.CompareTo("LHC16k") == 0 ||
            namePeriod.CompareTo("LHC16l") == 0 ||
            namePeriod.CompareTo("LHC16m") == 0 ||
            namePeriod.CompareTo("LHC16n") == 0 ||
            namePeriod.CompareTo("LHC16o") == 0 ||
            namePeriod.CompareTo("LHC16p") == 0 )       return k16pp13TeV;
  else if ( namePeriod.CompareTo("LHC16qt") == 0 ||
            namePeriod.CompareTo("LHC16q") == 0 ||
            namePeriod.CompareTo("LHC16t") == 0 )       return k16pPb5023GeV;
  else if ( namePeriod.CompareTo("LHC16rs") == 0 ||
            namePeriod.CompareTo("LHC16r") == 0 ||
            namePeriod.CompareTo("LHC16s") == 0 )       return k16pPb8TeV;
  else if ( namePeriod.CompareTo("LHC17cr") == 0 ||
            namePeriod.CompareTo("LHC17c") == 0 ||
            namePeriod.CompareTo("LHC17d") == 0 ||
            namePeriod.CompareTo("LHC17e") == 0 ||
            namePeriod.CompareTo("LHC17f") == 0 ||
            namePeriod.CompareTo("LHC17h") == 0 ||
            namePeriod.CompareTo("LHC17i") == 0 ||
            namePeriod.CompareTo("LHC17j") == 0 ||
            namePeriod.CompareTo("LHC17k") == 0 ||
            namePeriod.CompareTo("LHC17l") == 0 ||
            namePeriod.CompareTo("LHC17m") == 0 ||
            namePeriod.CompareTo("LHC17o") == 0 ||
            namePeriod.CompareTo("LHC17r") == 0 )       return k17pp13TeV;
  else if ( namePeriod.CompareTo("LHC17g") == 0 )       return k17pp13TeVLow;
  else if ( namePeriod.CompareTo("LHC17pq") == 0 ||
            namePeriod.CompareTo("LHC17p") == 0 ||
            namePeriod.CompareTo("LHC17q") == 0  )      return k17pp5TeV;
  else if ( namePeriod.CompareTo("LHC17n") == 0 )       return k17XeXe5440GeV;
  else if ( namePeriod.CompareTo("LHC18b") == 0 ||
            namePeriod.CompareTo("LHC18d") == 0 ||
            namePeriod.CompareTo("LHC18e") == 0 ||
            namePeriod.CompareTo("LHC18f") == 0 ||
            namePeriod.CompareTo("LHC18g") == 0 ||
            namePeriod.CompareTo("LHC18h") == 0 ||
            namePeriod.CompareTo("LHC18i") == 0 ||
            namePeriod.CompareTo("LHC18j") == 0 ||
            namePeriod.CompareTo("LHC18k") == 0 ||
            namePeriod.CompareTo("LHC18l") == 0 ||
            namePeriod.CompareTo("LHC18m") == 0 ||
            namePeriod.CompareTo("LHC18n") == 0 ||
            namePeriod.CompareTo("LHC18o") == 0 ||
            namePeriod.CompareTo("LHC18p") == 0 )       return k18pp13TeV;
  else if ( namePeriod.CompareTo("LHC18c") == 0 )       return k18pp13TeVLow;
  else if ( namePeriod.CompareTo("LHC18q") == 0 ||
            namePeriod.CompareTo("LHC18r") == 0 ||
            namePeriod.CompareTo("LHC18qr") == 0 )      return k18PbPb5TeV;
  else return kNoMC;
}

//________________________________________________________________________
void AliCaloPhotonCuts::SetLogBinningXTH1(TH1* histoRebin){
  TAxis *axisafter  = histoRebin->GetXaxis();
  Int_t bins        = axisafter->GetNbins();
  Double_t from     = axisafter->GetXmin();
  Double_t to       = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0]        = from;
  Double_t factor   = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
  return;
}


//________________________________________________________________________
void AliCaloPhotonCuts::SetLogBinningXTH2(TH2* histoRebin){
  TAxis *axisafter  = histoRebin->GetXaxis();
  Int_t bins        = axisafter->GetNbins();
  Double_t from     = axisafter->GetXmin();
  Double_t to       = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0]        = from;
  Double_t factor   = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
  return;
}

//________________________________________________________________________
void AliCaloPhotonCuts::SetLogBinningYTH2(TH2* histoRebin){
  TAxis *axisafter  = histoRebin->GetYaxis();
  Int_t bins        = axisafter->GetNbins();
  Double_t from     = axisafter->GetXmin();
  Double_t to       = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0]        = from;
  Double_t factor   = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
  return;
}

//________________________________________________________________________
Double_t AliCaloPhotonCuts::GetDistanceBetweenClusters(AliVCluster* cluster1, AliVCluster* cluster2){

  Float_t clusPos1[3]   = {0,0,0};
  cluster1->GetPosition(clusPos1);
  TVector3 clusterVector1(clusPos1[0],clusPos1[1],clusPos1[2]);
  Double_t etaCluster1  = clusterVector1.Eta();
  Double_t phiCluster1  = clusterVector1.Phi();
  if (phiCluster1 < 0) phiCluster1 += 2*TMath::Pi();

  Float_t clusPos2[3]   = {0,0,0};
  cluster2->GetPosition(clusPos2);
  TVector3 clusterVector2(clusPos2[0],clusPos2[1],clusPos2[2]);
  Double_t etaCluster2  = clusterVector2.Eta();
  Double_t phiCluster2  = clusterVector2.Phi();
  if (phiCluster2 < 0) phiCluster2 += 2*TMath::Pi();

  Double_t deltaEta     = TMath::Abs(etaCluster1-etaCluster2);
  Double_t deltaPhi     = 0;
  if (phiCluster1 > phiCluster2){
    deltaPhi            = phiCluster1-phiCluster2;
    if (deltaPhi > TMath::Pi())
      deltaPhi          = 2*TMath::Pi()-deltaPhi;
  } else {
    deltaPhi            = phiCluster2-phiCluster1;
    if (deltaPhi > TMath::Pi())
      deltaPhi          = 2*TMath::Pi()-deltaPhi;
  }
  Double_t r            = TMath::Sqrt(deltaPhi*deltaPhi+deltaEta*deltaEta);
  return r;
}

//________________________________________________________________________
TString AliCaloPhotonCuts::GetCutNumber(){
   // returns TString with current cut number
   return fCutStringRead;
}


//___________________________________________________________________
// Check if the cluster highest energy tower is exotic.
//___________________________________________________________________
Bool_t AliCaloPhotonCuts::IsExoticCluster( AliVCluster *cluster, AliVEvent *event, Float_t &energyStar ) {

  if (!cluster) {
    AliInfo("Cluster pointer null!");
    return kFALSE;
  }

  // case where exotic clusters are determined in correction framework
  if ( fUseExoticCluster == 3){
    if(cluster->GetIsExotic()){
      return kTRUE;
    } else {
      return kFALSE;
    }
  }

  energyStar              = 0;

  AliVCaloCells* cells    = NULL;
  if (fClusterType == 1 || fClusterType == 3 || fClusterType == 4)
    cells                 = event->GetEMCALCells();
  else if (fClusterType ==2 )
    return kFALSE;
//     cells                 = event->GetPHOSCells();


  Int_t largestCellID     = FindLargestCellInCluster(cluster,event);
  Float_t ecell1          = cells->GetCellAmplitude(largestCellID); ;
  Float_t eCross          = GetECross(largestCellID,cells);
  energyStar              = ecell1+eCross;

  if (ecell1 < fExoticMinEnergyCell)
    return kFALSE;

  if (1-eCross/ecell1 > fExoticEnergyFracCluster) {
    return kTRUE;
  } else if ( fUseExoticCluster == 2 && (cluster->E() > fExoticMinEnergyTCard) ){
    for(Int_t i = 1; i < cluster->GetNCells(); i++){  // check if cells of cluster are all in the same T-Card
      if( !IsAbsIDsFromTCard( cluster->GetCellAbsId(0), cluster->GetCellAbsId(i)) ){
        return kFALSE;
      }
    }
    return kTRUE;
  } else {
    return kFALSE;
  }
  return kFALSE;
}

//___________________________________________________________________________
// Calculate the energy in the cross around the energy of a given cell.
// Used in exotic clusters/cells rejection.
//___________________________________________________________________________
Float_t AliCaloPhotonCuts::GetECross( Int_t absID, AliVCaloCells* cells ){

  Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1;
  fGeomEMCAL->GetCellIndex(absID,imod,iTower,iIphi,iIeta);
  fGeomEMCAL->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);

  // Get close cells index, energy and time, not in corners
  Int_t absID1 = -1;
  Int_t absID2 = -1;

  if ( iphi < AliEMCALGeoParams::fgkEMCALRows-1) absID1 = fGeomEMCAL->GetAbsCellIdFromCellIndexes(imod, iphi+1, ieta);
  if ( iphi > 0 )                                absID2 = fGeomEMCAL->GetAbsCellIdFromCellIndexes(imod, iphi-1, ieta);

  // In case of cell in eta = 0 border, depending on SM shift the cross cell index
  Int_t absID3 = -1;
  Int_t absID4 = -1;
  if ( fClusterType == 1 && ieta == AliEMCALGeoParams::fgkEMCALCols-1 && !(imod%2) ){
    absID3 = fGeomEMCAL-> GetAbsCellIdFromCellIndexes(imod+1, iphi, 0);
    absID4 = fGeomEMCAL-> GetAbsCellIdFromCellIndexes(imod,   iphi, ieta-1);
  } else if ( fClusterType == 1 && ieta == 0 && imod%2 ) {
    absID3 = fGeomEMCAL-> GetAbsCellIdFromCellIndexes(imod,   iphi, ieta+1);
    absID4 = fGeomEMCAL-> GetAbsCellIdFromCellIndexes(imod-1, iphi, AliEMCALGeoParams::fgkEMCALCols-1);
  } else {
    if ( ieta < AliEMCALGeoParams::fgkEMCALCols-1 )
      absID3 = fGeomEMCAL-> GetAbsCellIdFromCellIndexes(imod, iphi, ieta+1);
    if ( ieta > 0 )
      absID4 = fGeomEMCAL-> GetAbsCellIdFromCellIndexes(imod, iphi, ieta-1);
  }

  Float_t  ecell1  = 0, ecell2  = 0, ecell3  = 0, ecell4  = 0;


  // Do not include bad channels found in analysis,
  if (AcceptCellByBadChannelMap(absID1)) ecell1 = cells->GetCellAmplitude(absID1); ;
  if (AcceptCellByBadChannelMap(absID2)) ecell2 = cells->GetCellAmplitude(absID2); ;
  if (AcceptCellByBadChannelMap(absID3)) ecell3 = cells->GetCellAmplitude(absID3); ;
  if (AcceptCellByBadChannelMap(absID4)) ecell4 = cells->GetCellAmplitude(absID4); ;

  return ecell1+ecell2+ecell3+ecell4;
}

//_______________________________________________________________________________
// Determine whether cell would be accepted according to bad channel map
//_______________________________________________________________________________
Bool_t AliCaloPhotonCuts::AcceptCellByBadChannelMap(Int_t absID ){
  if (fClusterType == 1 || fClusterType == 3 || fClusterType == 4){
    if(!fGeomEMCAL) {
      AliFatal("No instance of the geometry is available");
      return kFALSE;
    }
    if ( absID < 0 || absID >= 24*48*fGeomEMCAL->GetNumberOfSuperModules() )
      return kFALSE;

    Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1;
    fGeomEMCAL->GetCellIndex(absID,imod,iTower,iIphi,iIeta);
    fGeomEMCAL->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,iphi,ieta);

    // Do not include bad channels found in analysis,
    Int_t status = 0;
    if (fEMCALRecUtils->GetEMCALChannelStatus(imod, ieta, iphi, status) == 0 )
      return kTRUE;
    else
      return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________________________
/// Check if two cells are in the same T-Card.
/// Exotic cluster rejection for very high energy.
/// derived from AliAnalysisTaskEMCALPhotonIsolation
//________________________________________________________________________________________
Bool_t  AliCaloPhotonCuts::IsAbsIDsFromTCard(Int_t absId1, Int_t absId2) const {

  if(absId1 == absId2) return kFALSE;

    // Check if in same SM, if not for sure not same TCard
  Int_t sm1 = fGeomEMCAL->GetSuperModuleNumber(absId1);
  Int_t sm2 = fGeomEMCAL->GetSuperModuleNumber(absId2);
  if ( sm1 != sm2 ) return kFALSE ;

    // Get the column and row of each absId
  Int_t iTower = -1, iIphi = -1, iIeta = -1;

  Int_t col1, row1;
  fGeomEMCAL->GetCellIndex(absId1,sm1,iTower,iIphi,iIeta);
  fGeomEMCAL->GetCellPhiEtaIndexInSModule(sm1,iTower,iIphi, iIeta,row1,col1);

  Int_t col2, row2;
  fGeomEMCAL->GetCellIndex(absId2,sm2,iTower,iIphi,iIeta);
  fGeomEMCAL->GetCellPhiEtaIndexInSModule(sm2,iTower,iIphi, iIeta,row2,col2);

  Int_t row0 = Int_t(row1-row1%8);
  Int_t col0 = Int_t(col1-col1%2);

  Int_t rowDiff0 = row2-row0;
  Int_t colDiff0 = col2-col0;

    // TCard is of 2x8 towers
  if ( colDiff0 >=0 && colDiff0 < 2 && rowDiff0 >=0 && rowDiff0 < 8 )
    return kTRUE;
  else
    return kFALSE;
}

//_______________________________________________________________________________
// Classify clusters for TM effi
// 0: Neutral cluster
// 1: Neutral cluster sub charged
// 2: Gamma cluster
// 3: Gamma cluster sub charged
// 4: Gamma conv cluster
// 5: Charged cluster
// 6: electron
// 7: primary charged
//_______________________________________________________________________________
Int_t AliCaloPhotonCuts::ClassifyClusterForTMEffi(AliVCluster* cluster, AliVEvent* event, AliMCEvent* mcEvent, Bool_t isESD){
  Int_t* mclabelsCluster  = cluster->GetLabels();
  Int_t classification    = -1;

  const AliVVertex* primVtxMC   = mcEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  if (isESD){
    if (cluster->GetNLabels()>0){
      TParticle* particleLead   = (TParticle*)mcEvent->Particle(mclabelsCluster[0]);
      Double_t charge           = ((TParticlePDG*)particleLead->GetPDG())->Charge();
      if (charge == 0 || charge == -0){
        classification        = 0;
        if (particleLead->GetPdgCode() == 22)
          classification      = 2;
      } else {
        classification        = 5;
        if (particleLead->GetPdgCode() == 11 || particleLead->GetPdgCode() == -11){
          classification      = 6;
          if (particleLead->GetMother(0) > -1)
            if (((TParticle*)mcEvent->Particle(particleLead->GetMother(0)))->GetPdgCode() == 22)
              classification  = 4;
        } else {
            Double_t deltaX = particleLead->Vx() - mcProdVtxX;
            Double_t deltaY = particleLead->Vy() - mcProdVtxY;
            Double_t deltaZ = particleLead->Vz() - mcProdVtxZ;

            Double_t realRadius3D = TMath::Sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ);
            if (realRadius3D < 5.0)
              classification  = 7;
        }
      }
      if ((classification == 0 || classification == 2) && cluster->GetNLabels() > 1){
        Bool_t goOut = kFALSE;
        for (UInt_t i = 1; (i< cluster->GetNLabels() && !goOut); i++){
          TParticle* particleSub    = (TParticle*)mcEvent->Particle(mclabelsCluster[i]);
          Double_t charge           = ((TParticlePDG*)particleSub->GetPDG())->Charge();
          if (!(charge == 0 || charge == -0)){
            classification++;
            goOut = kTRUE;
          }
        }
      }
    }
  } else {
    if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
    if (fAODMCTrackArray == NULL){
      AliError("No MC particle list available in AOD");
      return -1;
    }

    if (cluster->GetNLabels()>0){
      AliAODMCParticle* particleLead    = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(mclabelsCluster[0]));
      Short_t charge                    = particleLead->Charge();
      if (charge == 0 ){
        classification        = 0;
        if (particleLead->GetPdgCode() == 22)
          classification      = 2;
      } else {
        classification        = 5;
        if (particleLead->GetPdgCode() == 11 || particleLead->GetPdgCode() == -11) {
          classification      = 6;
          if (particleLead->GetMother() > -1){
            AliAODMCParticle* mother    = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particleLead->GetMother()));
            if (mother->GetPdgCode() == 22)
              classification  = 4;
          }
        } else {
          Double_t deltaX = particleLead->Xv() - mcProdVtxX;
          Double_t deltaY = particleLead->Yv() - mcProdVtxY;
          Double_t deltaZ = particleLead->Zv() - mcProdVtxZ;

          Double_t realRadius3D = TMath::Sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ);
          if (realRadius3D < 5.0)
            classification    = 7;
        }
      }
      if ((classification == 0 || classification == 2) && cluster->GetNLabels() > 1){
        Bool_t goOut = kFALSE;
        for (UInt_t i = 1; (i< cluster->GetNLabels() && !goOut); i++){
          AliAODMCParticle* particleSub = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(mclabelsCluster[i]));
          Double_t charge               = particleSub->Charge();
          if ( charge != 0 ){
            classification++;
            goOut = kTRUE;
          }
        }
      }
    }
  }
  return classification;
}

//_______________________________________________________________________________
std::vector<Int_t> AliCaloPhotonCuts::GetVectorMatchedTracksToCluster(AliVEvent* event, AliVCluster* cluster){
  vector<Int_t> labelsMatched(0);
  if(!fUseDistTrackToCluster || fUseElectronClusterCalibration) return labelsMatched;

  if (fUsePtDepTrackToCluster == 0)
    labelsMatched = fCaloTrackMatcher->GetMatchedTrackIDsForCluster(event, cluster->GetID(), fMaxDistTrackToClusterEta, -fMaxDistTrackToClusterEta,
                                                                    fMaxDistTrackToClusterPhi, fMinDistTrackToClusterPhi);
  else if (fUsePtDepTrackToCluster == 1)
    labelsMatched = fCaloTrackMatcher->GetMatchedTrackIDsForCluster(event, cluster->GetID(), fFuncPtDepEta, fFuncPtDepPhi);

  return labelsMatched;
}

//_______________________________________________________________________________
Bool_t AliCaloPhotonCuts::GetClosestMatchedTrackToCluster(AliVEvent* event, AliVCluster* cluster, Int_t &trackLabel){
  if(!fUseDistTrackToCluster || fUseElectronClusterCalibration) return kFALSE;
  vector<Int_t> labelsMatched = GetVectorMatchedTracksToCluster(event,cluster);

  if((Int_t) labelsMatched.size()<1) return kFALSE;

  Float_t dEta = -100;
  Float_t dPhi = -100;
  Int_t smallestDistTrack = -1;
  Float_t smallestDist = 100;
  Bool_t isSucessful = kFALSE;
  for (Int_t i = 0; i < (Int_t)labelsMatched.size(); i++){
    AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(labelsMatched.at(i)));
    if(!fCaloTrackMatcher->GetTrackClusterMatchingResidual(currTrack->GetID(), cluster->GetID(), dEta, dPhi)) continue;

    Float_t tempDist = TMath::Sqrt((dEta*dEta)+(dPhi*dPhi));
    if(tempDist < smallestDist){
      smallestDist = tempDist;
      smallestDistTrack = labelsMatched.at(i);
      isSucessful = kTRUE;
    }
  }
  if(!isSucessful) return kFALSE;
  trackLabel = smallestDistTrack;
  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliCaloPhotonCuts::GetHighestPtMatchedTrackToCluster(AliVEvent* event, AliVCluster* cluster, Int_t &trackLabel){
  if(!fUseDistTrackToCluster || fUseElectronClusterCalibration) return kFALSE;
  vector<Int_t> labelsMatched = GetVectorMatchedTracksToCluster(event,cluster);

  if((Int_t) labelsMatched.size()<1) return kFALSE;

  Int_t highestPtTrack = -1;
  Float_t highestPt = 0;
  Bool_t isSucessful = kFALSE;
  for (Int_t i = 0; i < (Int_t)labelsMatched.size(); i++){
    AliVTrack* currTrack  = dynamic_cast<AliVTrack*>(event->GetTrack(labelsMatched.at(i)));
    if(currTrack->Pt() > highestPt){
      highestPt = currTrack->Pt();
      highestPtTrack = labelsMatched.at(i);
      isSucessful = kTRUE;
    }
  }
  if(!isSucessful) return kFALSE;
  trackLabel = highestPtTrack;
  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliCaloPhotonCuts::IsClusterPi0(AliVEvent *event,  AliMCEvent* mcEvent, AliVCluster *cluster){

  Int_t* mclabelsCluster  = cluster->GetLabels();
  // check if esd or aod
  AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent *aodev = 0;
  if (!esdev) {
    aodev = dynamic_cast<AliAODEvent*>(event);
    if (!aodev) {
      AliError("Task needs AOD or ESD event, returning");
      return kFALSE;
    }
  }

  if(esdev){
    // check if cluster is from pi0
    if (cluster->GetNLabels()>0){
      TParticle* particleLead   = (TParticle*)mcEvent->Particle(mclabelsCluster[0]);
      if(particleLead->GetPdgCode() == 22 || particleLead->GetPdgCode() == -11 || particleLead->GetPdgCode() == 11){
        if (particleLead->GetMother(0) > -1){
          TParticle* motherDummy = mcEvent->Particle(particleLead->GetMother(0));
//          printf("mother pdg = %d\n",motherDummy->GetPdgCode());
          if(motherDummy->GetPdgCode() == 111) return kTRUE;
          if(motherDummy->GetPdgCode() == 22){ // check also conversions
            UInt_t whileCounter = 0;// to be 100% sure against infinite while loop
            while(motherDummy->GetMother(0) > -1 && whileCounter++ < 5){
              motherDummy = mcEvent->Particle(motherDummy->GetMother(0));
//              printf("grandmother pdg = %d\n",motherDummy->GetPdgCode());
              if(motherDummy->GetPdgCode() == 111) return kTRUE;
            }
          }
        }
      }
    }
  }else if(aodev){ // same procedure for AODsesdt
    if(!fAODMCTrackArray) fAODMCTrackArray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
    if (fAODMCTrackArray == NULL){
      AliError("No MC particle list available in AOD");
      return kFALSE;
    }
    // check if cluster is from pi0
    if (cluster->GetNLabels()>0){
      AliAODMCParticle* particleLead = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(mclabelsCluster[0]));
      if(particleLead->GetPdgCode() == 22 || particleLead->GetPdgCode() == -11 || particleLead->GetPdgCode() == 11){
        if (particleLead->GetMother() > -1){
          AliAODMCParticle* motherDummy = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(particleLead->GetMother()));
//          printf("mother pdg = %d\n",motherDummy->GetPdgCode());
          if(motherDummy->GetPdgCode() == 111) return kTRUE;
          if(motherDummy->GetPdgCode() == 22){ // check also conversions
            UInt_t whileCounter = 0;// to be 100% sure against infinite while loop
            while(motherDummy->GetMother() > -1 && whileCounter++ < 5){
              motherDummy    = static_cast<AliAODMCParticle*>(fAODMCTrackArray->At(motherDummy->GetMother()));
//              printf("grandmother pdg = %d\n",motherDummy->GetPdgCode());
              if(motherDummy->GetPdgCode() == 111) return kTRUE;
            }
          }
        }
      }
    }
  }

  return kFALSE;
}

//_________________________________________________________________________________
// function to find possible conversion candidates and clean up cluster array
// fMaxMGGRecConv  determines cut off for inv Mass
// lower energetic photon is rejected to to likely worse energy resolution
//_________________________________________________________________________________
Bool_t AliCaloPhotonCuts::CheckForReconstructedConversionPairs( vector<AliAODConversionPhoton*> &vecPhotons,
                                                                vector<Int_t> &vecReject
                                                              ){

  Bool_t rejected   = kFALSE;
  if(vecPhotons.size()>0){
    for(Int_t firstGammaIndex=0;firstGammaIndex<(Int_t)vecPhotons.size();firstGammaIndex++){
      AliAODConversionPhoton *gamma1=vecPhotons.at(firstGammaIndex);
      if (gamma1==NULL){
        CheckVectorForIndexAndAdd(vecReject, firstGammaIndex,kTRUE);
        continue;
      }
      TLorentzVector photon1;
      photon1.SetPxPyPzE (gamma1->Px(), gamma1->Py(), gamma1->Pz(), gamma1->E() );

      for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<(Int_t)vecPhotons.size();secondGammaIndex++){
        AliAODConversionPhoton *gamma2=vecPhotons.at(secondGammaIndex);
        if (gamma2==NULL){
          CheckVectorForIndexAndAdd(vecReject, secondGammaIndex,kTRUE);
          continue;
        }
        TLorentzVector photon2;
        photon2.SetPxPyPzE (gamma2->Px(), gamma2->Py(), gamma2->Pz(), gamma2->E() );

        TLorentzVector mesonCand;
        mesonCand = photon1+photon2;
        fHistInvMassDiCluster->Fill(mesonCand.M(), mesonCand.Pt());
        if (mesonCand.M() < fMaxMGGRecConv){
          fHistInvMassConvFlagging->Fill(mesonCand.M(), mesonCand.Pt());
          if (CheckVectorForIndexAndAdd(vecReject, firstGammaIndex,kFALSE)) AliDebug(2,"1st gamma already rejected");
          if (CheckVectorForIndexAndAdd(vecReject, secondGammaIndex,kFALSE)) AliDebug(2,"2nd gamma already rejected");
          rejected = kTRUE;
          if (gamma1->E() < gamma2->E())
            CheckVectorForIndexAndAdd(vecReject, firstGammaIndex,kTRUE);
          else
            CheckVectorForIndexAndAdd(vecReject, secondGammaIndex,kTRUE);
        }
      }
    }
  }
  if (rejected){
    AliDebug(2,"================================================================================");
    AliDebug(2,"================================================================================");
    AliDebug(2,Form("array of rejected gamma from initial: %i/%i ", (int)vecReject.size(), (int)vecPhotons.size()));
    for (Int_t index=0;index<(Int_t)vecReject.size();index++)
       AliDebug(2,Form("index: %i", (int)vecReject.at(index)));
    AliDebug(2,"================================================================================");
  }
  return rejected;
}

//_________________________________________________________________________________
// function to check if integer is contained in current vector
// if not already contained will be aded if addIndex == kTRUE
//_________________________________________________________________________________
Bool_t AliCaloPhotonCuts::CheckVectorForIndexAndAdd(vector<Int_t> &vec, Int_t tobechecked, Bool_t addIndex )
{
  if(tobechecked > -1){
    vector<Int_t>::iterator it;
    it = find (vec.begin(), vec.end(), tobechecked);
    if (it != vec.end()) return true;
    else if (addIndex){
      vec.push_back(tobechecked);
      return kFALSE;
    }
  }
  return false;
}
