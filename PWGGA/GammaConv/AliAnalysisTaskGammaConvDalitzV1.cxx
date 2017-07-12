/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Pedro Gonzalez, Pedro Ladron de Guevara, Ernesto Lopez Torres, *
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
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliMCEvent.h"
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
#include "AliAnalysisTaskGammaConvDalitzV1.h"
#include "AliEventplane.h"
#include <vector>


ClassImp( AliAnalysisTaskGammaConvDalitzV1 )

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskGammaConvDalitzV1::AliAnalysisTaskGammaConvDalitzV1():
	fV0Reader(NULL),
    fV0ReaderName("V0ReaderV1"),
	fElecSelector(NULL),
	fBGHandler(NULL),
	fESDEvent(NULL),
	fMCEvent(NULL),
	fCutFolder(NULL),
	fESDList(NULL),
	fBackList(NULL),
	fMotherList(NULL),
	fTrueList(NULL),
	fMCList(NULL),
	fQAFolder(NULL),
	fOutputContainer(0),
	fReaderGammas(NULL),
	fSelectorElectronIndex(0),
	fSelectorPositronIndex(0),
	fGoodGammas(NULL),
	fGoodVirtualGammas(NULL),
	fGoodElectrons(NULL),
	fGoodPositrons(NULL),
	fCutEventArray(NULL),
	fCutGammaArray(NULL),
	fCutElectronArray(NULL),
	fCutMesonArray(NULL),
	fEventCuts(NULL),
	fConversionCuts(NULL),
	hESDConvGammaPt(NULL),
	hESDConvGammaEta(NULL),
	//hESDConvGammaZR(NULL),
	sESDConvGammaZR(NULL),
	sESDConvGammaXY(NULL),
	hESDDalitzElectronPt(NULL),
	hESDDalitzPositronPt(NULL),
	hESDDalitzElectronPhi(NULL),
	hESDDalitzPositronPhi(NULL),
	hESDDalitzElectronAfterPt(NULL),
	hESDDalitzPositronAfterPt(NULL),
	hESDDalitzElectronAfterEta(NULL),
	hESDDalitzElectronAfterEtaPCut(NULL),
	hESDDalitzPositronAfterEta(NULL),
	hESDDalitzPositronAfterEtaPCut(NULL),
	hESDDalitzElectronAfterPhi(NULL),
	hESDDalitzPositronAfterPhi(NULL),
	hESDDalitzElectronAfterNClsITS(NULL),
	hESDDalitzElectronAfterNClsITSPCut(NULL),
	hESDDalitzPositronAfterNClsITS(NULL),
	hESDDalitzPositronAfterNClsITSPCut(NULL),
	hESDDalitzElectronAfterNFindClsTPC(NULL),
	hESDDalitzElectronAfterNFindClsTPCPCut(NULL),
	hESDDalitzPositronAfterNFindClsTPC(NULL),
	hESDDalitzPositronAfterNFindClsTPCPCut(NULL),
	hESDDalitzElectronAfterNClsTPC(NULL),
	hESDDalitzElectronAfterNClsTPCPCut(NULL),
	hESDDalitzPositronAfterNClsTPC(NULL),
	hESDDalitzPositronAfterNClsTPCPCut(NULL),
	hESDDalitzElectronAfterNCrossedRowsTPC(NULL),
	hESDDalitzElectronAfterNCrossedRowsTPCPCut(NULL),
	hESDDalitzPositronAfterNCrossedRowsTPC(NULL),
	hESDDalitzPositronAfterNCrossedRowsTPCPCut(NULL),
	hESDDalitzPosEleAfterDCAxy(NULL),
	hESDDalitzPosEleAfterDCAz(NULL),
	hESDDalitzElectronAfterTPCdEdxVsP(NULL),
	hESDDalitzPositronAfterTPCdEdxVsP(NULL),
	hESDDalitzElectronAfterTPCdEdxSignalVsP(NULL),
	hESDDalitzPositronAfterTPCdEdxSignalVsP(NULL),
	hESDDalitzElectronAfterTPCdEdxVsEta(NULL),
	hESDDalitzPositronAfterTPCdEdxVsEta(NULL),
	hESDDalitzElectronAfterTPCdEdxVsPhi(NULL),
	hESDDalitzPositronAfterTPCdEdxVsPhi(NULL),
	hESDEposEnegPsiPairDPhi(NULL),
	hESDEposEnegPsiPairEta(NULL),
	hESDEposEnegDPhiEta(NULL),
	hESDEposEnegInvMassPt(NULL),
	hESDTruePi0EposEnegInvMassPi0Pt(NULL),
	hESDEposEnegLikeSignBackInvMassPt(NULL),
	hESDMotherInvMassPt(NULL),
	hESDPi0MotherInvMassPt(NULL),
	hESDPi0MotherDiffInvMassPt(NULL),
	hESDPi0MotherDiffLimInvMassPt(NULL),
	hESDMotherInvMassOpeningAngleGammaElectron(NULL),
	sESDMotherInvMassPtZM(NULL),
	hESDMotherBackInvMassPt(NULL),
	sESDMotherBackInvMassPtZM(NULL),
	hESDMotherPi0PtY(NULL),
	hESDMotherPi0PtAlpha(NULL),
	hESDMotherPi0PtOpenAngle(NULL),
	sESDMotherDalitzPlot(NULL),
	hMCAllGammaPt(NULL),
	hMCAllGammaPi0Pt(NULL),
	hMCConvGammaPt(NULL),
	hMCConvGammaPtR(NULL),
	hMCConvGammaRSPt(NULL),
	hMCConvGammaPi0Pt(NULL),
	hMCAllPositronsPt(NULL),
	hMCDecayPositronPi0Pt(NULL),
	hMCAllElectronsPt(NULL),
	hMCDecayElectronPi0Pt(NULL),
	hMCConvGammaEta(NULL),
	hMCConvGammaR(NULL),
	hMCAllPositronsEta(NULL),
	hMCAllElectronsEta(NULL),
	hMCPi0DalitzGammaPt(NULL),
	hMCPi0DalitzElectronPt(NULL),
	hMCPi0DalitzPositronPt(NULL),
	hMCPi0Pt(NULL),
	hMCPi0WOWeightPt(NULL),
	hMCPi0GGPt(NULL),
	hMCEtaPt(NULL),
	hMCEtaGGPt(NULL), 
	hMCPi0InAccPt(NULL),
	hMCPi0WOWeightInAccPt(NULL),
	hMCPi0InAccOpeningAngleGammaElectron(NULL),
	sMCPi0DalitzPlot(NULL),
	hMCEtaInAccPt(NULL),
	hMCChiCPt(NULL),
	hMCChiCInAccPt(NULL),
	hMCPi0EposEnegInvMassPt(NULL),
	hMCEtaEposEnegInvMassPt(NULL),
	hESDEposEnegTruePi0DalitzInvMassPt(NULL),
	hESDEposEnegTruePrimPi0DalitzInvMass(NULL),
	hESDEposEnegTruePi0DalitzPsiPairDPhi(NULL),
	hESDEposEnegTruePi0DalitzPsiPairMC(NULL),
	hESDEposEnegTruePi0DalitzPsiPairEta(NULL),
	hESDEposEnegTruePi0DalitzDPhiEta(NULL),
	hESDEposEnegTrueEtaDalitzInvMassPt(NULL),
	hESDEposEnegTruePrimEtaDalitzInvMass(NULL),
	hESDEposEnegTrueEtaDalitzPsiPairDPhi(NULL),
	hESDEposEnegTruePhotonInvMassPt(NULL),
	hESDEposEnegTrueInvMassPt(NULL),
	hESDEposEnegTrueMotherInvMassPt(NULL),
	hESDEposEnegTruePhotonPsiPairDPhi(NULL),
	hESDEposEnegTruePhotonPsiPairDPhiPtCut(NULL),
	hESDEposEnegTrueJPsiInvMassPt(NULL),
	hESDTrueMotherChiCInvMassPt(NULL),
	hESDTrueMotherChiCDiffInvMassPt(NULL),
	hESDTrueMotherInvMassPt(NULL),
	hESDTrueMotherW0WeightsInvMassPt(NULL),
	hESDTrueMotherDalitzInvMassPt(NULL),
	hESDTrueMotherPi0GGInvMassPt(NULL),
	hESDTrueMotherPi0GGW0WeightsInvMassPt(NULL),
	hESDTruePi0PtY(NULL),
	hESDTruePi0PtAlpha(NULL),
	hESDTruePi0PtOpenAngle(NULL),
	sESDTruePi0DalitzPlot(NULL),
	hESDTruePrimaryMotherPi0GGInvMassPt(NULL),
	hESDTrueSecondaryMotherPi0GGInvMassPt(NULL),
	hESDTruePrimaryMotherInvMassMCPt(NULL),
	hESDTruePrimaryMotherInvMassPt(NULL),
	hESDTruePrimaryMotherW0WeightingInvMassPt(NULL),
	hESDTruePrimaryPi0DalitzESDPtMCPt(NULL),
	hESDTrueSecondaryMotherInvMassPt(NULL),
	hESDTrueSecondaryMotherFromK0sInvMassPt(NULL),
	hESDTrueBckGGInvMassPt(NULL),
	hESDTrueBckContInvMassPt(NULL),
	hESDTrueMotherGGInvMassPt(NULL),
	hESDTrueConvGammaPt(NULL),
	hESDTrueConvGammaPtMC(NULL),
	hESDTrueConvGammaR(NULL),
	hESDTrueConvGammaRMC(NULL),
	hESDTruePositronPt(NULL),
	hESDTrueElectronPt(NULL),
	hESDTrueSecConvGammaPt(NULL),
	hESDTrueSecPositronPt(NULL),
	hESDTrueSecElectronPt(NULL),
	hESDTruePi0DalitzConvGammaPt(NULL),
	hESDTruePi0DalitzConvGammaR(NULL),
	hESDTruePi0DalitzPositronPt(NULL),
	hESDTruePi0DalitzPositronPtMB(NULL),
	hESDTruePi0DalitzElectronPt(NULL),
	hESDTruePi0DalitzElectronPtMB(NULL),
	hESDTruePi0DalitzSecConvGammaPt(NULL),
	hESDTruePi0DalitzSecPositronPt(NULL),
	hESDTruePi0DalitzSecElectronPt(NULL),
	hNEvents(NULL),
	hNGoodESDTracks(NULL),
	hNGoodESDTracksVsNGoodGammas(NULL),
	hNGoodESDTracksVsNGoodVGammas(NULL),
	fHistoSPDClusterTrackletBackground(NULL),
	hNV0Tracks(NULL),
	hEtaShift(NULL),
	fHistoDoubleCountTruePi0InvMassPt(NULL),
	fHistoDoubleCountTrueEtaInvMassPt(NULL),
	fHistoDoubleCountTrueConvGammaRPt(NULL),
	fVectorDoubleCountTruePi0s(0),
	fVectorDoubleCountTrueEtas(0),
	fVectorDoubleCountTrueConvGammas(0),
	fEventPlaneAngle(-100),
	fRandom(0),
	fUnsmearedPx(NULL),
	fUnsmearedPy(NULL),
	fUnsmearedPz(NULL),
	fUnsmearedE(NULL),
	fUnsmearedVPx(NULL),
	fUnsmearedVPy(NULL),
	fUnsmearedVPz(NULL),
	fUnsmearedVE(NULL),
	fnCuts(0),
	fiCut(0),
	fNumberOfESDTracks(0),
	fNumberOfESDTrackskBoth(0),
	fNVirtualGammas(0),
	fMoveParticleAccordingToVertex(kFALSE),
	fIsHeavyIon(0),
	fDoMesonAnalysis(kTRUE),
	fDoChicAnalysis(kFALSE),
	fDoMesonQA(0),
	fSetProductionVertextoVGamma(kTRUE),
	fIsFromMBHeader(kTRUE),
	fIsMC(kFALSE),
	fDoTHnSparse(kTRUE),
	fDoHistoDalitzMassLog(kFALSE),
    fDoMaterialBudgetWeightingOfGammasForTrueMesons(kFALSE)
{

}

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskGammaConvDalitzV1::AliAnalysisTaskGammaConvDalitzV1( const char* name ):
	AliAnalysisTaskSE(name),
	fV0Reader(NULL),
    fV0ReaderName("V0ReaderV1"),
	fElecSelector(NULL),
	fBGHandler(NULL),
	fESDEvent(NULL),
	fMCEvent(NULL),
	fCutFolder(NULL),
	fESDList(NULL),
	fBackList(NULL),
	fMotherList(NULL),
	fTrueList(NULL),
	fMCList(NULL),
	fQAFolder(NULL),
	fOutputContainer(0),
	fReaderGammas(NULL),
	fSelectorElectronIndex(0),
	fSelectorPositronIndex(0),
	fGoodGammas(NULL),
	fGoodVirtualGammas(NULL),
	fGoodElectrons(NULL),
	fGoodPositrons(NULL),
	fCutEventArray(NULL),
	fCutGammaArray(NULL),
	fCutElectronArray(NULL),
	fCutMesonArray(NULL),
	fEventCuts(NULL),
	fConversionCuts(NULL),
	hESDConvGammaPt(NULL),
	hESDConvGammaEta(NULL),
	//hESDConvGammaZR(NULL),
	sESDConvGammaZR(NULL),
	sESDConvGammaXY(NULL),
	hESDDalitzElectronPt(NULL),
	hESDDalitzPositronPt(NULL),
	hESDDalitzElectronPhi(NULL),
	hESDDalitzPositronPhi(NULL),
	hESDDalitzElectronAfterPt(NULL),
	hESDDalitzPositronAfterPt(NULL),
	hESDDalitzElectronAfterEta(NULL),
	hESDDalitzElectronAfterEtaPCut(NULL),
	hESDDalitzPositronAfterEta(NULL),
	hESDDalitzPositronAfterEtaPCut(NULL),
	hESDDalitzElectronAfterPhi(NULL),
	hESDDalitzPositronAfterPhi(NULL),
	hESDDalitzElectronAfterNClsITS(NULL),
	hESDDalitzElectronAfterNClsITSPCut(NULL),
	hESDDalitzPositronAfterNClsITS(NULL),
	hESDDalitzPositronAfterNClsITSPCut(NULL),
	hESDDalitzElectronAfterNFindClsTPC(NULL),
	hESDDalitzElectronAfterNFindClsTPCPCut(NULL),
	hESDDalitzPositronAfterNFindClsTPC(NULL),
	hESDDalitzPositronAfterNFindClsTPCPCut(NULL),
	hESDDalitzElectronAfterNClsTPC(NULL),
	hESDDalitzElectronAfterNClsTPCPCut(NULL),
	hESDDalitzPositronAfterNClsTPC(NULL),
	hESDDalitzPositronAfterNClsTPCPCut(NULL),
	hESDDalitzElectronAfterNCrossedRowsTPC(NULL),
	hESDDalitzElectronAfterNCrossedRowsTPCPCut(NULL),
	hESDDalitzPositronAfterNCrossedRowsTPC(NULL),
	hESDDalitzPositronAfterNCrossedRowsTPCPCut(NULL),
	hESDDalitzPosEleAfterDCAxy(NULL),
	hESDDalitzPosEleAfterDCAz(NULL),
	hESDDalitzElectronAfterTPCdEdxVsP(NULL),
	hESDDalitzPositronAfterTPCdEdxVsP(NULL),
	hESDDalitzElectronAfterTPCdEdxSignalVsP(NULL),
	hESDDalitzPositronAfterTPCdEdxSignalVsP(NULL),
	hESDDalitzElectronAfterTPCdEdxVsEta(NULL),
	hESDDalitzPositronAfterTPCdEdxVsEta(NULL),
	hESDDalitzElectronAfterTPCdEdxVsPhi(NULL),
	hESDDalitzPositronAfterTPCdEdxVsPhi(NULL),
	hESDEposEnegPsiPairDPhi(NULL),
	hESDEposEnegPsiPairEta(NULL),
	hESDEposEnegDPhiEta(NULL),
	hESDEposEnegInvMassPt(NULL),
	hESDTruePi0EposEnegInvMassPi0Pt(NULL),
	hESDEposEnegLikeSignBackInvMassPt(NULL),
	hESDMotherInvMassPt(NULL),
	hESDPi0MotherInvMassPt(NULL),
	hESDPi0MotherDiffInvMassPt(NULL),
	hESDPi0MotherDiffLimInvMassPt(NULL),
	hESDMotherInvMassOpeningAngleGammaElectron(NULL),
	sESDMotherInvMassPtZM(NULL),
	hESDMotherBackInvMassPt(NULL),
	sESDMotherBackInvMassPtZM(NULL),
	hESDMotherPi0PtY(NULL),
	hESDMotherPi0PtAlpha(NULL),
	hESDMotherPi0PtOpenAngle(NULL),
	sESDMotherDalitzPlot(NULL),
	hMCAllGammaPt(NULL),
	hMCAllGammaPi0Pt(NULL),
	hMCConvGammaPt(NULL),
	hMCConvGammaPtR(NULL),
	hMCConvGammaRSPt(NULL),
	hMCConvGammaPi0Pt(NULL),
	hMCAllPositronsPt(NULL),
	hMCDecayPositronPi0Pt(NULL),
	hMCAllElectronsPt(NULL),
	hMCDecayElectronPi0Pt(NULL),
	hMCConvGammaEta(NULL),
	hMCConvGammaR(NULL),
	hMCAllPositronsEta(NULL),
	hMCAllElectronsEta(NULL),
	hMCPi0DalitzGammaPt(NULL),
	hMCPi0DalitzElectronPt(NULL),
	hMCPi0DalitzPositronPt(NULL),
	hMCPi0Pt(NULL),
	hMCPi0WOWeightPt(NULL),
	hMCPi0GGPt(NULL),
	hMCEtaPt(NULL),
	hMCEtaGGPt(NULL),
	hMCPi0InAccPt(NULL),
	hMCPi0WOWeightInAccPt(NULL),
	hMCPi0InAccOpeningAngleGammaElectron(NULL),
	sMCPi0DalitzPlot(NULL),
	hMCEtaInAccPt(NULL),
	hMCChiCPt(NULL),
	hMCChiCInAccPt(NULL),
	hMCPi0EposEnegInvMassPt(NULL),
	hMCEtaEposEnegInvMassPt(NULL),
	hESDEposEnegTruePi0DalitzInvMassPt(NULL),
	hESDEposEnegTruePrimPi0DalitzInvMass(NULL),
	hESDEposEnegTruePi0DalitzPsiPairDPhi(NULL),
	hESDEposEnegTruePi0DalitzPsiPairMC(NULL),
	hESDEposEnegTruePi0DalitzPsiPairEta(NULL),
	hESDEposEnegTruePi0DalitzDPhiEta(NULL),
	hESDEposEnegTrueEtaDalitzInvMassPt(NULL),
	hESDEposEnegTruePrimEtaDalitzInvMass(NULL),
	hESDEposEnegTrueEtaDalitzPsiPairDPhi(NULL),
	hESDEposEnegTruePhotonInvMassPt(NULL),
	hESDEposEnegTrueInvMassPt(NULL),
	hESDEposEnegTrueMotherInvMassPt(NULL),
	hESDEposEnegTruePhotonPsiPairDPhi(NULL),
	hESDEposEnegTruePhotonPsiPairDPhiPtCut(NULL),
	hESDEposEnegTrueJPsiInvMassPt(NULL),
	hESDTrueMotherChiCInvMassPt(NULL),
	hESDTrueMotherChiCDiffInvMassPt(NULL),
	hESDTrueMotherInvMassPt(NULL),
	hESDTrueMotherW0WeightsInvMassPt(NULL),
	hESDTrueMotherDalitzInvMassPt(NULL),
	hESDTrueMotherPi0GGInvMassPt(NULL),
	hESDTrueMotherPi0GGW0WeightsInvMassPt(NULL),
	hESDTruePi0PtY(NULL),
	hESDTruePi0PtAlpha(NULL),
	hESDTruePi0PtOpenAngle(NULL),
	sESDTruePi0DalitzPlot(NULL),
	hESDTruePrimaryMotherPi0GGInvMassPt(NULL),
	hESDTrueSecondaryMotherPi0GGInvMassPt(NULL),
	hESDTruePrimaryMotherInvMassMCPt(NULL),
	hESDTruePrimaryMotherInvMassPt(NULL),
	hESDTruePrimaryMotherW0WeightingInvMassPt(NULL),
	hESDTruePrimaryPi0DalitzESDPtMCPt(NULL),
	hESDTrueSecondaryMotherInvMassPt(NULL),
	hESDTrueSecondaryMotherFromK0sInvMassPt(NULL),
	hESDTrueBckGGInvMassPt(NULL),
	hESDTrueBckContInvMassPt(NULL),
	hESDTrueMotherGGInvMassPt(NULL),
	hESDTrueConvGammaPt(NULL),
	hESDTrueConvGammaPtMC(NULL),
	hESDTrueConvGammaR(NULL),
	hESDTrueConvGammaRMC(NULL),
	hESDTruePositronPt(NULL),
	hESDTrueElectronPt(NULL),
	hESDTrueSecConvGammaPt(NULL),
	hESDTrueSecPositronPt(NULL),
	hESDTrueSecElectronPt(NULL),
	hESDTruePi0DalitzConvGammaPt(NULL),
	hESDTruePi0DalitzConvGammaR(NULL),
	hESDTruePi0DalitzPositronPt(NULL),
	hESDTruePi0DalitzPositronPtMB(NULL),
	hESDTruePi0DalitzElectronPt(NULL),
	hESDTruePi0DalitzElectronPtMB(NULL),
	hESDTruePi0DalitzSecConvGammaPt(NULL),
	hESDTruePi0DalitzSecPositronPt(NULL),
	hESDTruePi0DalitzSecElectronPt(NULL),
	hNEvents(NULL),
	hNGoodESDTracks(NULL),
	hNGoodESDTracksVsNGoodGammas(NULL),
	hNGoodESDTracksVsNGoodVGammas(NULL),
	fHistoSPDClusterTrackletBackground(NULL),
	hNV0Tracks(NULL),
	hEtaShift(NULL),
	fHistoDoubleCountTruePi0InvMassPt(NULL),
	fHistoDoubleCountTrueEtaInvMassPt(NULL),
	fHistoDoubleCountTrueConvGammaRPt(NULL),
	fVectorDoubleCountTruePi0s(0),
	fVectorDoubleCountTrueEtas(0),
	fVectorDoubleCountTrueConvGammas(0),
	fEventPlaneAngle(-100),
	fRandom(0),
	fUnsmearedPx(NULL),
	fUnsmearedPy(NULL),
	fUnsmearedPz(NULL),
	fUnsmearedE(NULL),
	fUnsmearedVPx(NULL),
	fUnsmearedVPy(NULL),
	fUnsmearedVPz(NULL),
	fUnsmearedVE(NULL),
	fnCuts(0),
	fiCut(0),
	fNumberOfESDTracks(0),
	fNumberOfESDTrackskBoth(0),
	fNVirtualGammas(0),
	fMoveParticleAccordingToVertex(kFALSE),
	fIsHeavyIon(0),
	fDoMesonAnalysis(kTRUE),
	fDoChicAnalysis(kFALSE),
	fDoMesonQA(0),
	fSetProductionVertextoVGamma(kTRUE),
	fIsFromMBHeader(kTRUE),
	fIsMC(kFALSE),
	fDoTHnSparse(kTRUE),
	fDoHistoDalitzMassLog(kFALSE),
    fDoMaterialBudgetWeightingOfGammasForTrueMesons(kFALSE)
{
	DefineOutput(1, TList::Class());
}

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskGammaConvDalitzV1::~AliAnalysisTaskGammaConvDalitzV1()
{
	//
	// virtual destructor
	//

	cout<<"Destructor"<<endl;

	if(fGoodGammas){
		delete fGoodGammas;
		fGoodGammas = 0x0;
	}
	if(fGoodVirtualGammas){
		delete fGoodVirtualGammas;
		fGoodVirtualGammas = 0x0;
	}
	if(fGoodElectrons){
		delete fGoodElectrons;
		fGoodElectrons = 0x0;
	}
	if(fGoodPositrons){
		delete fGoodPositrons;
		fGoodPositrons = 0x0;
	}
	if(fBGHandler){
		delete[] fBGHandler;
		fBGHandler = 0x0;
	}
	
}

//___________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::InitBack(){

	const Int_t nDim = 4;
	Int_t nBins[nDim] = {800,250,7,4};
	Double_t xMin[nDim] = {0,0, 0,0};
	Double_t xMax[nDim] = {0.8,25,7,4};
	
	if( fDoTHnSparse ) {
	sESDMotherInvMassPtZM = new THnSparseF*[fnCuts];
	sESDMotherBackInvMassPtZM = new THnSparseF*[fnCuts];
	}

	fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];
	
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
	
		TString cutstringEvent		= ((AliConvEventCuts*)fCutEventArray->At(iCut))->GetCutNumber();
		TString cutstringElectron	= ((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->GetCutNumber();
		TString cutstringMeson		= ((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->GetCutNumber();
		TString cutstringGamma		= ((AliConversionPhotonCuts*)fCutGammaArray->At(iCut))->GetCutNumber();
		
		Int_t collisionSystem = atoi((TString)(((AliConvEventCuts*)fCutEventArray->At(iCut))->GetCutNumber())(0,1));
		Int_t centMin = atoi((TString)(((AliConvEventCuts*)fCutEventArray->At(iCut))->GetCutNumber())(1,1));
		Int_t centMax = atoi((TString)(((AliConvEventCuts*)fCutEventArray->At(iCut))->GetCutNumber())(2,1));
		
		if(collisionSystem == 1 || collisionSystem == 2 ||
			collisionSystem == 5 || collisionSystem == 8 ||
			collisionSystem == 9){
			centMin = centMin*10;
			centMax = centMax*10; 
		        if(centMax ==0 && centMax!=centMin) centMax=100;
		}
		else if(collisionSystem == 3 || collisionSystem == 6){
			centMin = centMin*5;
			centMax = centMax*5;
		}
		else if(collisionSystem == 4 || collisionSystem == 7){
			centMin = ((centMin*5)+45);
			centMax = ((centMax*5)+45);
		}

		if( fDoTHnSparse ) {
		fBackList[iCut] = new TList();
		fBackList[iCut]->SetName(Form("%s_%s_%s_%s Back histograms",cutstringEvent.Data(),cutstringGamma.Data(),cutstringElectron.Data(),cutstringMeson.Data()));
		fBackList[iCut]->SetOwner(kTRUE);
		fCutFolder[iCut]->Add(fBackList[iCut]);

		sESDMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_m","Back_Back_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
		fBackList[iCut]->Add(sESDMotherBackInvMassPtZM[iCut]);

		fMotherList[iCut] = new TList();
		fMotherList[iCut]->SetName(Form("%s_%s_%s_%s Mother histograms",cutstringEvent.Data(),cutstringGamma.Data(),cutstringElectron.Data(),cutstringMeson.Data()));
		fMotherList[iCut]->SetOwner(kTRUE);
		fCutFolder[iCut]->Add(fMotherList[iCut]);

		sESDMotherInvMassPtZM[iCut] = new THnSparseF("Back_Mother_InvMass_Pt_z_m","Back_Mother_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
		fMotherList[iCut]->Add(sESDMotherInvMassPtZM[iCut]);
		}
		
		if(((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->BackgroundHandlerType() == 0){
			
		
		fBGHandler[iCut] = new AliGammaConversionAODBGHandler(
																collisionSystem,centMin,centMax,
																((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->GetNumberOfBGEvents(),
																((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->UseTrackMultiplicity(),
																1,8,5);
		}
		
		
			
	}
}

//______________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::UserCreateOutputObjects()
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

	fGoodVirtualGammas = new TList();
	fGoodVirtualGammas->SetOwner(kTRUE);

	
	fCutFolder				= new TList*[fnCuts];
	fESDList				= new TList*[fnCuts];
	
	if( fDoTHnSparse ){
	fBackList				= new TList*[fnCuts];
	fMotherList				= new TList*[fnCuts];
	}
	hNEvents				= new TH1I*[fnCuts];
	hNGoodESDTracks			= new TH1I*[fnCuts];
	hNV0Tracks				= new TH1I*[fnCuts];
	hEtaShift				= new TProfile*[fnCuts];
	hESDConvGammaPt			= new TH1F*[fnCuts];
	hESDConvGammaEta		= new TH1F*[fnCuts];
	
	hESDDalitzElectronPt	= new TH1F*[fnCuts];
	hESDDalitzPositronPt	= new TH1F*[fnCuts];
	hESDDalitzElectronPhi	= new TH1F*[fnCuts];
	hESDDalitzPositronPhi	= new TH1F*[fnCuts];
	
	const Int_t nBinsMassPair = 800;
	Double_t xMinMassPair = 0.0;
	Double_t xMaxMassPair = 0.8;
	
	if( fDoMesonQA > 0) {
		fQAFolder  					= new TList*[fnCuts];  
		hNGoodESDTracksVsNGoodGammas			= new TH2F*[fnCuts];
		hNGoodESDTracksVsNGoodVGammas			= new TH2F*[fnCuts]; 
		fHistoSPDClusterTrackletBackground		= new TH2F*[fnCuts];
		hESDDalitzElectronAfterPt			= new TH1F*[fnCuts];
		hESDDalitzPositronAfterPt			= new TH1F*[fnCuts];
		hESDDalitzElectronAfterEta			= new TH1F*[fnCuts];
		hESDDalitzElectronAfterEtaPCut			= new TH1F*[fnCuts];
		hESDDalitzPositronAfterEta			= new TH1F*[fnCuts];
		hESDDalitzPositronAfterEtaPCut			= new TH1F*[fnCuts];
		hESDDalitzElectronAfterPhi			= new TH1F*[fnCuts];
		hESDDalitzPositronAfterPhi			= new TH1F*[fnCuts];
		hESDDalitzElectronAfterNClsITS			= new TH1F*[fnCuts];
		hESDDalitzElectronAfterNClsITSPCut              = new TH1F*[fnCuts];
		hESDDalitzPositronAfterNClsITS			= new TH1F*[fnCuts];
		hESDDalitzPositronAfterNClsITSPCut		= new TH1F*[fnCuts];
		hESDDalitzElectronAfterNFindClsTPC		= new TH2F*[fnCuts];
		hESDDalitzElectronAfterNFindClsTPCPCut          = new TH1F*[fnCuts];
		hESDDalitzPositronAfterNFindClsTPC		= new TH2F*[fnCuts];
		hESDDalitzPositronAfterNFindClsTPCPCut		= new TH1F*[fnCuts];
		hESDDalitzElectronAfterNClsTPC			= new TH2F*[fnCuts];
		hESDDalitzElectronAfterNClsTPCPCut              = new TH1F*[fnCuts];
		hESDDalitzPositronAfterNClsTPC			= new TH2F*[fnCuts];
		hESDDalitzPositronAfterNClsTPCPCut		= new TH1F*[fnCuts];
		hESDDalitzElectronAfterNCrossedRowsTPC 		= new TH2F*[fnCuts];
		hESDDalitzElectronAfterNCrossedRowsTPCPCut      = new TH1F*[fnCuts];
		hESDDalitzPositronAfterNCrossedRowsTPC		= new TH2F*[fnCuts];
		hESDDalitzPositronAfterNCrossedRowsTPCPCut	= new TH1F*[fnCuts];
		hESDDalitzPosEleAfterDCAxy			= new TH2F*[fnCuts];
		hESDDalitzPosEleAfterDCAz			= new TH2F*[fnCuts];
		hESDDalitzElectronAfterTPCdEdxVsP		= new TH2F*[fnCuts];
		hESDDalitzPositronAfterTPCdEdxVsP		= new TH2F*[fnCuts];
		hESDDalitzElectronAfterTPCdEdxSignalVsP 	= new TH2F*[fnCuts];
		hESDDalitzPositronAfterTPCdEdxSignalVsP 	= new TH2F*[fnCuts];
		hESDDalitzElectronAfterTPCdEdxVsEta 		= new TH2F*[fnCuts];
		hESDDalitzPositronAfterTPCdEdxVsEta 		= new TH2F*[fnCuts];
		hESDDalitzElectronAfterTPCdEdxVsPhi 		= new TH2F*[fnCuts];
		hESDDalitzPositronAfterTPCdEdxVsPhi 		= new TH2F*[fnCuts];
		hESDEposEnegPsiPairDPhi				= new TH2F*[fnCuts];
		hESDEposEnegPsiPairEta				= new TH2F*[fnCuts];
		hESDEposEnegDPhiEta				= new TH2F*[fnCuts];
		hESDEposEnegDPhiEta				= new TH2F*[fnCuts];
		hESDEposEnegInvMassPt				= new TH2F*[fnCuts];
		hESDTruePi0EposEnegInvMassPi0Pt			= new TH2F*[fnCuts];
		hESDEposEnegLikeSignBackInvMassPt 		= new TH2F*[fnCuts];
		hESDMotherInvMassOpeningAngleGammaElectron 	= new TH1F*[fnCuts];
		hESDMotherPi0PtY 				= new TH2F*[fnCuts];
		hESDMotherPi0PtAlpha 				= new TH2F*[fnCuts];
		hESDMotherPi0PtOpenAngle 			= new TH2F*[fnCuts];
		
		
		if ( fDoMesonQA > 1 ){
		      sESDConvGammaZR 	    = new THnSparseF*[fnCuts];
		      sESDConvGammaXY  	    = new THnSparseF*[fnCuts];
		      sESDMotherDalitzPlot  = new THnSparseF*[fnCuts];
		}
		
		if( fDoHistoDalitzMassLog ){
			  
			  xMinMassPair = 0.001;
			  xMaxMassPair = 0.801;
			  
		}
		
				
	}
	
	
	
	hESDMotherInvMassPt				= new TH2F*[fnCuts];
	
	if(fDoChicAnalysis) {
		hESDPi0MotherInvMassPt 				 	= new TH2F*[fnCuts];
		hESDPi0MotherDiffInvMassPt				= new TH2F*[fnCuts];
		hESDPi0MotherDiffLimInvMassPt				= new TH2F*[fnCuts];
	}
	hESDMotherBackInvMassPt	= new TH2F*[fnCuts];


	for(Int_t iCut = 0; iCut<fnCuts;iCut++){

		TString cutstringEvent 		= ((AliConvEventCuts*)fCutEventArray->At(iCut))->GetCutNumber();
		TString cutstringElectron 	= ((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->GetCutNumber();
		TString cutstringMeson		= ((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->GetCutNumber();
		TString cutstringGamma 		= ((AliConversionPhotonCuts*)fCutGammaArray->At(iCut))->GetCutNumber();

		fCutFolder[iCut] = new TList();
		fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s_%s",cutstringEvent.Data(),cutstringGamma.Data(),cutstringElectron.Data(),cutstringMeson.Data()));
		fCutFolder[iCut]->SetOwner(kTRUE);
		fOutputContainer->Add(fCutFolder[iCut]);

		fESDList[iCut] = new TList();
		fESDList[iCut]->SetName(Form("%s_%s_%s_%s ESD histograms",cutstringEvent.Data(),cutstringGamma.Data(),cutstringElectron.Data(),cutstringMeson.Data()));
		fESDList[iCut]->SetOwner(kTRUE);
		
		hNEvents[iCut] = new TH1I("NEvents","NEvents",13,-0.5,12.5);
		hNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Missing MC");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problem");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
		fESDList[iCut]->Add(hNEvents[iCut]);

		if(fIsHeavyIon == 1 || fIsHeavyIon == 2) hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",3000,0,3000);
		else hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",200,0,200);
		fESDList[iCut]->Add(hNGoodESDTracks[iCut]);
		
		
		if(fIsHeavyIon == 1) hNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",30000,0,30000);
		else if(fIsHeavyIon == 2) hNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",2500,0,2500);
		else hNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",1500,0,1500);
		fESDList[iCut]->Add(hNV0Tracks[iCut]);
		
		hEtaShift[iCut] = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
		fESDList[iCut]->Add(hEtaShift[iCut]);

		hESDConvGammaPt[iCut] = new TH1F("ESD_ConvGamma_Pt","ESD_ConvGamma_Pt",250,0,25);
		fESDList[iCut]->Add(hESDConvGammaPt[iCut]);
		
		hESDConvGammaEta[iCut] = new TH1F("ESD_ConvGamma_Eta","ESD_ConvGamma_Eta",600,-1.5,1.5);
		fESDList[iCut]->Add(hESDConvGammaEta[iCut]);
		
		hESDDalitzElectronPt[iCut] = new TH1F("ESD_DalitzElectron_Pt","ESD_DalitzElectron_Pt",1000,0,25);
		fESDList[iCut]->Add(hESDDalitzElectronPt[iCut]);

		hESDDalitzPositronPt[iCut] = new TH1F("ESD_DalitzPositron_Pt","ESD_DalitzPositron_Pt",1000,0,25);
		fESDList[iCut]->Add(hESDDalitzPositronPt[iCut]);
		
		
		hESDDalitzElectronPhi[iCut] = new TH1F("ESD_DalitzElectron_Phi","ESD_DalitzElectron_Phi",360,0,2*TMath::Pi());
		fESDList[iCut]->Add(hESDDalitzElectronPhi[iCut]);

		hESDDalitzPositronPhi[iCut] = new TH1F("ESD_DalitzPositron_Phi","ESD_DalitzPositron_Phi",360,0,2*TMath::Pi());
		fESDList[iCut]->Add(hESDDalitzPositronPhi[iCut]);
		
		
		
		if ( fDoMesonQA > 0 ) {

			fQAFolder[iCut] = new TList();
			fQAFolder[iCut]->SetName(Form("%s_%s_%s_%s QA histograms",cutstringEvent.Data(),cutstringGamma.Data(),cutstringElectron.Data(),cutstringMeson.Data()));
			fQAFolder[iCut]->SetOwner(kTRUE);
				
			const Int_t kPtBins=110;
			Double_t binsPtDummy[kPtBins+1];
			const Int_t kPBins = 109;
			Double_t binsPDummy[kPBins+1];
			binsPtDummy[0]=0.0;
			binsPDummy[0]=0.05;
			
			
			
			
			
			for(Int_t i=1;i<kPtBins+1;i++){
				if(binsPtDummy[i-1]+0.05<1.01)
						binsPtDummy[i]=binsPtDummy[i-1]+0.05;
				else
						binsPtDummy[i]=binsPtDummy[i-1]+0.1;		
			}
			for(Int_t i=1; i<kPBins+1;i++){	
				if( binsPDummy[i-1]+0.05<1.01)
						binsPDummy[i] = binsPDummy[i-1]+0.05;
				else
					binsPDummy[i] = binsPDummy[i-1]+0.1;
			}
				
			hNGoodESDTracksVsNGoodGammas[iCut] = new TH2F("hNGoodESDTracksVsNGoodGammas","hNGoodESDTracksVsNGoodGammas",200,-0.5,199.5,100,-0.5,99.5);
			fQAFolder[iCut]->Add(hNGoodESDTracksVsNGoodGammas[iCut]);
			
			hNGoodESDTracksVsNGoodVGammas[iCut] = new TH2F("hNGoodESDTracksVsNVGoodVGammas","hNGoodESDTracksVsNGoodVGammas",200,-0.5,199.5,100,-0.5,99.5);
			fQAFolder[iCut]->Add(hNGoodESDTracksVsNGoodVGammas[iCut]);

			fHistoSPDClusterTrackletBackground[iCut] = new TH2F("SPD tracklets vs SPD clusters","SPD tracklets vs SPD clusters",100,0,200,250,0,1000);
			fQAFolder[iCut]->Add(fHistoSPDClusterTrackletBackground[iCut]);
				
			
			hESDDalitzElectronAfterPt[iCut] = new TH1F("ESD_DalitzElectron_After_Pt","ESD_DalitzElectron_After_Pt",1000,0,25);
			fQAFolder[iCut]->Add(hESDDalitzElectronAfterPt[iCut]);

			hESDDalitzPositronAfterPt[iCut] = new TH1F("ESD_DalitzPositron_After_Pt","ESD_DalitzPositron_After_Pt",1000,0,25);
			fQAFolder[iCut]->Add(hESDDalitzPositronAfterPt[iCut]);
				
			hESDDalitzElectronAfterEta[iCut] = new TH1F("ESD_DalitzElectron_After_Eta","ESD_DalitzElectron_After_Eta",600,-1.5,1.5);
			fQAFolder[iCut]->Add(hESDDalitzElectronAfterEta[iCut]);
			
			hESDDalitzElectronAfterEtaPCut[iCut] = new TH1F("ESD_DalitzElectron_After_Eta_PCut","ESD_DalitzElectron_After_Eta_PCut",600,-1.5,1.5);
			fQAFolder[iCut]->Add(hESDDalitzElectronAfterEtaPCut[iCut]);
			
			hESDDalitzPositronAfterEta[iCut] = new TH1F("ESD_DalitzPositron_After_Eta","ESD_DalitzElectron_After_Eta",600,-1.5,1.5);
			fQAFolder[iCut]->Add(hESDDalitzPositronAfterEta[iCut]);
			
			hESDDalitzPositronAfterEtaPCut[iCut] = new TH1F("ESD_DalitzPositron_After_Eta_PCut","ESD_DalitzElectron_After_Eta_PCut",600,-1.5,1.5);
			fQAFolder[iCut]->Add(hESDDalitzPositronAfterEtaPCut[iCut]);
			
			hESDDalitzElectronAfterPhi[iCut] = new TH1F("ESD_DalitzElectron_After_Phi","ESD_DalitzElectron_After_Phi",360,0,2*TMath::Pi());
			fQAFolder[iCut]->Add(hESDDalitzElectronAfterPhi[iCut]);

			hESDDalitzPositronAfterPhi[iCut] = new TH1F("ESD_DalitzPositron_After_Phi","ESD_DalitzPositron_After_Phi",360,0,2*TMath::Pi());
			fQAFolder[iCut]->Add(hESDDalitzPositronAfterPhi[iCut]);
			
			hESDDalitzElectronAfterNClsITS[iCut]  = new TH1F("ESD_DalitzElectron_After_NClsITS","ESD_DalitzElectron_After_NClsITS",7,-0.5,6.5);
			fQAFolder[iCut]->Add(hESDDalitzElectronAfterNClsITS[iCut]);
			
			hESDDalitzElectronAfterNClsITSPCut[iCut]  = new TH1F("ESD_DalitzElectron_After_NClsITS_PCut","ESD_DalitzElectron_After_NClsITS_PCut",7,-0.5,6.5);
			fQAFolder[iCut]->Add(hESDDalitzElectronAfterNClsITSPCut[iCut]);
			
			hESDDalitzPositronAfterNClsITS[iCut]  = new TH1F("ESD_DalitzPositron_After_NClsITS","ESD_DalitzPositron_After_NClsITS",7,-0.5,6.5);
			fQAFolder[iCut]->Add(hESDDalitzPositronAfterNClsITS[iCut]);
			
			hESDDalitzPositronAfterNClsITSPCut[iCut]  = new TH1F("ESD_DalitzPositron_After_NClsITS_PCut","ESD_DalitzPositron_After_NClsITS_PCut",7,-0.5,6.5);
			fQAFolder[iCut]->Add(hESDDalitzPositronAfterNClsITSPCut[iCut]);
						
			hESDDalitzElectronAfterNFindClsTPC[iCut]  = new TH2F("ESD_DalitzElectron_After_NFindClsTPC","ESD_DalitzElectron_After_NFindClsTPC",60,0,1.5,kPtBins,binsPtDummy);
			fQAFolder[iCut]->Add(hESDDalitzElectronAfterNFindClsTPC[iCut]);
			
			hESDDalitzElectronAfterNFindClsTPCPCut[iCut]  = new TH1F("ESD_DalitzElectron_After_NFindClsTPC_PCut","ESD_DalitzElectron_After_NFindClsTPC_PCut",60,0,1.5);
			fQAFolder[iCut]->Add(hESDDalitzElectronAfterNFindClsTPCPCut[iCut]);
			
			
			hESDDalitzPositronAfterNFindClsTPC[iCut]  = new TH2F("ESD_DalitzPositron_After_NFindClsTPC","ESD_DalitzPositron_After_NFindClsTPC",60,0,1.5,kPtBins,binsPtDummy);
			fQAFolder[iCut]->Add(hESDDalitzPositronAfterNFindClsTPC[iCut]);
			
			hESDDalitzPositronAfterNFindClsTPCPCut[iCut]  = new TH1F("ESD_DalitzPositron_After_NFindClsTPC_PCut","ESD_DalitzPositron_After_NFindClsTPC_PCut",60,0,1.5);
			fQAFolder[iCut]->Add(hESDDalitzPositronAfterNFindClsTPCPCut[iCut]);
			
			
			hESDDalitzElectronAfterNClsTPC[iCut]  = new TH2F("ESD_DalitzElectron_After_NClsTPC","ESD_DalitzElectron_After_NClsTPC",200,-0.5,199.5,kPtBins,binsPtDummy);
			fQAFolder[iCut]->Add(hESDDalitzElectronAfterNClsTPC[iCut]);
			
			hESDDalitzElectronAfterNClsTPCPCut[iCut]  = new TH1F("ESD_DalitzElectron_After_NClsTPC_PCut","ESD_DalitzElectron_After_NClsTPC_PCut",200,-0.5,199.5);
			fQAFolder[iCut]->Add(hESDDalitzElectronAfterNClsTPCPCut[iCut]);
			
			hESDDalitzPositronAfterNClsTPC[iCut]  = new TH2F("ESD_DalitzPositron_After_NClsTPC","ESD_DalitzPositron_After_NClsTPC",200,-0.5,199.5,kPtBins,binsPtDummy);
			fQAFolder[iCut]->Add(hESDDalitzPositronAfterNClsTPC[iCut]);
			
			hESDDalitzPositronAfterNClsTPCPCut[iCut]  = new TH1F("ESD_DalitzPositron_After_NClsTPC_PCut","ESD_DalitzPositron_After_NClsTPC_PCut",200,-0.5,199.5);
			fQAFolder[iCut]->Add(hESDDalitzPositronAfterNClsTPCPCut[iCut]);
					
			hESDDalitzPositronAfterNCrossedRowsTPC[iCut] =  new TH2F("ESD_DalitzPositron_After_NCrossedRowsTPC","ESD_DalitzPositron_After_NCrossedRowsTPC",165,-0.5,164.5,kPtBins,binsPtDummy);
			fQAFolder[iCut]->Add(hESDDalitzPositronAfterNCrossedRowsTPC[iCut]);
			
			hESDDalitzPositronAfterNCrossedRowsTPCPCut[iCut] =  new TH1F("ESD_DalitzPositron_After_NCrossedRowsTPC_PCut","ESD_DalitzPositron_After_NCrossedRowsTPC_PCut",165,-0.5,164.5);
			fQAFolder[iCut]->Add(hESDDalitzPositronAfterNCrossedRowsTPCPCut[iCut]);
						
			hESDDalitzElectronAfterNCrossedRowsTPC[iCut] =  new TH2F("ESD_DalitzElectron_After_NCrossedRowsTPC","ESD_DalitzElectron_After_NCrossedRowsTPC",165,-0.5,164.5,kPtBins,binsPtDummy);
			fQAFolder[iCut]->Add(hESDDalitzElectronAfterNCrossedRowsTPC[iCut]);
			
			hESDDalitzElectronAfterNCrossedRowsTPCPCut[iCut] =  new TH1F("ESD_DalitzElectron_After_NCrossedRowsTPC_PCut","ESD_DalitzElectron_After_NCrossedRowsTPC_PCut",165,-0.5,164.5);
			fQAFolder[iCut]->Add(hESDDalitzElectronAfterNCrossedRowsTPCPCut[iCut]);
					
			hESDDalitzPosEleAfterDCAxy[iCut] = new TH2F("ESD_DalitzPosEle_After_DCAxy","ESD_DalitzPosEle_After_DCAxy",124,-0.62,0.62,kPtBins,binsPtDummy);
			fQAFolder[iCut]->Add(hESDDalitzPosEleAfterDCAxy[iCut]);
			
			hESDDalitzPosEleAfterDCAz[iCut]  = new TH2F("ESD_DalitzPosEle_After_DCAz","ESD_DalitzPosEle_After_DCAz",200,-1.0,1.0,kPtBins,binsPtDummy);
			fQAFolder[iCut]->Add(hESDDalitzPosEleAfterDCAz[iCut]);
			
			hESDDalitzElectronAfterTPCdEdxVsP[iCut] = new TH2F("ESD_DalitzElectron_After_TPCdEdxVsP","ESD_DalitzElectron_After_TPCdEdxVsP_After_TPCdEdx",kPBins,binsPDummy,200,-10,10);
			SetLogBinningXTH2(hESDDalitzElectronAfterTPCdEdxVsP[iCut]);
			fQAFolder[iCut]->Add(hESDDalitzElectronAfterTPCdEdxVsP[iCut]);
			
			hESDDalitzPositronAfterTPCdEdxVsP[iCut] = new TH2F("ESD_DalitzPositron_After_TPCdEdxVsP","ESD_DalitzPositron_After_TPCdEdxVsP",kPBins,binsPDummy,200,-10,10);
			SetLogBinningXTH2(hESDDalitzPositronAfterTPCdEdxVsP[iCut]);
			fQAFolder[iCut]->Add(hESDDalitzPositronAfterTPCdEdxVsP[iCut]);
			
			hESDDalitzElectronAfterTPCdEdxSignalVsP[iCut] =new TH2F("ESD_DalitzElectron_After_TPCdEdxSignalVsP","ESD_DalitzElectron_After_TPCdEdxSignalVsP" ,kPBins,binsPDummy,200,0.0,200);
			SetLogBinningXTH2(hESDDalitzElectronAfterTPCdEdxSignalVsP[iCut]);
			fQAFolder[iCut]->Add(hESDDalitzElectronAfterTPCdEdxSignalVsP[iCut]);  
			
			hESDDalitzPositronAfterTPCdEdxSignalVsP[iCut] =new TH2F("ESD_DalitzPositron_After_TPCdEdxSignalVsP","ESD_DalitzPositron_After_TPCdEdxSignalVsP" ,kPBins,binsPDummy,200,0.0,200);
			SetLogBinningXTH2(hESDDalitzPositronAfterTPCdEdxSignalVsP[iCut]);
			fQAFolder[iCut]->Add(hESDDalitzPositronAfterTPCdEdxSignalVsP[iCut]); 
			
			hESDDalitzElectronAfterTPCdEdxVsEta[iCut] = new TH2F("ESD_DalitzElectron_After_TPCdEdxVsEta","ESD_DalitzElectron_After_TPCdEdxVsEta",140,-1.4,1.4,200,-10,10);
			fQAFolder[iCut]->Add(hESDDalitzElectronAfterTPCdEdxVsEta[iCut]);
			
			hESDDalitzPositronAfterTPCdEdxVsEta[iCut] = new  TH2F("ESD_DalitzPositron_After_TPCdEdxVsEta","ESD_DalitzPositron_After_TPCdEdxVsEta",140,-1.4,1.4,200,-10,10);
			fQAFolder[iCut]->Add(hESDDalitzPositronAfterTPCdEdxVsEta[iCut]);
			
			hESDDalitzElectronAfterTPCdEdxVsPhi[iCut] = new TH2F("ESD_DalitzElectron_After_TPCdEdxVsPhi","ESD_DalitzElectron_After_TPCdEdxVsPhi",180,0,2*TMath::Pi(),200,-10,10);
			fQAFolder[iCut]->Add(hESDDalitzElectronAfterTPCdEdxVsPhi[iCut]);
			
			hESDDalitzPositronAfterTPCdEdxVsPhi[iCut] = new TH2F("ESD_DalitzPositron_After_TPCdEdxVsPhi","ESD_DalitzPositron_After_TPCdEdxVsPhi",180,0,2*TMath::Pi(),200,-10,10);
			fQAFolder[iCut]->Add(hESDDalitzPositronAfterTPCdEdxVsPhi[iCut]);
					
			hESDEposEnegPsiPairDPhi[iCut] = new TH2F("ESD_EposEneg_PsiPair_DPhi","ESD_EposEneg_PsiPair_DPhi",100,-1.0,1.0,100,-1.0,1.0 );
			fQAFolder[iCut]->Add(hESDEposEnegPsiPairDPhi[iCut]);
			
			hESDEposEnegPsiPairEta[iCut] = new TH2F("ESD_EposEneg_PsiPair_Eta","ESD_EposEneg_PsiPair_Eta",100,-1.0,1.0,600,-1.5,1.5);
			fQAFolder[iCut]->Add(hESDEposEnegPsiPairEta[iCut]);
			
			hESDEposEnegDPhiEta[iCut] = new TH2F("ESD_EposEneg_DPhi_Eta","ESD_EposEneg_DPhi_Eta",100,-1.0,1.0,600,-1.5,1.5);
			fQAFolder[iCut]->Add(hESDEposEnegDPhiEta[iCut]);
			

			hESDEposEnegInvMassPt[iCut] = new TH2F("ESD_EposEneg_InvMassPt","ESD_EposEneg_InvMassPt",nBinsMassPair,xMinMassPair,xMaxMassPair,100,0.,10.);
			if(fDoHistoDalitzMassLog) SetLogBinningXTH2(hESDEposEnegInvMassPt[iCut]);
			fQAFolder[iCut]->Add(hESDEposEnegInvMassPt[iCut]);
			
				
			hESDEposEnegLikeSignBackInvMassPt[iCut]  = new TH2F("ESD_EposEneg_LikeSignBack_InvMassPt","ESD_EposEneg_LikeSignBack_InvMassPt",4000,0.0,2.,100,0.,10.);
			fQAFolder[iCut]->Add(hESDEposEnegLikeSignBackInvMassPt[iCut]);
			
			
			if( fDoMesonQA > 1 ) {
          
        const Int_t nDimRZ = 2;
        Int_t   nBinsRZ[nDimRZ] = { 1200, 800};
        Double_t xMinRZ[nDimRZ] = { -150, 0};
        Double_t xMaxRZ[nDimRZ] = {  150, 200};
        
        const Int_t nDimXY = 2;
        Int_t   nBinsXY[nDimXY] = { 1200, 1200};
        Double_t xMinXY[nDimXY] = { -150, -150};
        Double_t xMaxXY[nDimXY] = {  150,  150};
        
        sESDConvGammaZR[iCut] = new THnSparseF("ESD_ConvGamma_ZR","ESD_ConvGamma_ZR",nDimRZ,nBinsRZ,xMinRZ,xMaxRZ);
        fQAFolder[iCut]->Add(sESDConvGammaZR[iCut]);
        
        sESDConvGammaXY[iCut] = new THnSparseF("ESD_ConvGamma_XY","ESD_ConvGamma_XY",nDimXY,nBinsXY,xMinXY,xMaxXY);
        fQAFolder[iCut]->Add(sESDConvGammaXY[iCut]);
        
        
        const Int_t nDimDalPlot = 2;
        Int_t    nBinsDalPlot[nDimDalPlot] = {  4000, 4000};
        Double_t xMinDalPlot[nDimDalPlot]  = {    0,     0};
        Double_t xMaxDalPlot[nDimDalPlot]  = {  4.0,   4.0};
        
        
        sESDMotherDalitzPlot[iCut] = new THnSparseF("ESD_Mother_DalitzPlot","ESD_Mother_DalitzPlot",nDimDalPlot,nBinsDalPlot,xMinDalPlot,xMaxDalPlot);
        fQAFolder[iCut]->Add(sESDMotherDalitzPlot[iCut]);
        
			}
			fCutFolder[iCut]->Add(fQAFolder[iCut]);
		}
		
		hESDMotherInvMassPt[iCut] = new TH2F("ESD_DalitzMother_InvMass_Pt","ESD_DalitzMother_InvMass_Pt",800,0,0.8,250,0,25);
		fESDList[iCut]->Add(hESDMotherInvMassPt[iCut]);
		
		if( fDoMesonQA > 0 ) {
		  
		    hESDMotherInvMassOpeningAngleGammaElectron[iCut] = new TH1F("ESD_MotherInvMass_OpeningAngle_GammaElectron", "ESD_MotherInvMass_OpeningAngle_GammaElectron",100,0.,TMath::Pi());
		    fESDList[iCut]->Add(hESDMotherInvMassOpeningAngleGammaElectron[iCut]);
		    
		    hESDMotherPi0PtY[iCut] 		= new TH2F("ESD_MotherPi0_Pt_Y","ESD_MotherPi0_Pt_Y",150,0.03,15.,150,-1.5,1.5);            
		    SetLogBinningXTH2(hESDMotherPi0PtY[iCut]);
		    fESDList[iCut]->Add(hESDMotherPi0PtY[iCut]);
		    hESDMotherPi0PtAlpha[iCut] 		= new TH2F("ESD_MotherPi0_Pt_Alpha","ESD_MotherPi0_Pt_Alpha",150,0.03,15.,100,0,1);            
		    SetLogBinningXTH2(hESDMotherPi0PtAlpha[iCut]);
		    fESDList[iCut]->Add(hESDMotherPi0PtAlpha[iCut]);
		    hESDMotherPi0PtOpenAngle[iCut] 	= new TH2F("ESD_MotherPi0_Pt_OpenAngle","ESD_MotherPi0_Pt_OpenAngle",150,0.03,15.,100,0,TMath::Pi());
		    SetLogBinningXTH2(hESDMotherPi0PtOpenAngle[iCut]);
		    fESDList[iCut]->Add(hESDMotherPi0PtOpenAngle[iCut]);
                }
											
		
		if( fDoChicAnalysis) {
			hESDPi0MotherInvMassPt[iCut] = new TH2F("ESD_Pi0Mother_InvMass_Pt","ESD_Pi0Mother_InvMass_Pt",4000,0,4,250,0,25);
			fESDList[iCut]->Add(hESDPi0MotherInvMassPt[iCut]);
			hESDPi0MotherDiffInvMassPt[iCut] = new TH2F("ESD_Pi0Mother_DiffInvMass_Pt","ESD_Pi0Mother_DiffInvMass_Pt",2000,0,2,250,0,25);
			fESDList[iCut]->Add(hESDPi0MotherDiffInvMassPt[iCut]);
			hESDPi0MotherDiffLimInvMassPt[iCut] = new TH2F("ESD_Pi0Mother_DiffLimInvMass_Pt","ESD_Pi0Mother_DiffLimInvMass_Pt",2000,0,2,250,0,25);
			fESDList[iCut]->Add(hESDPi0MotherDiffLimInvMassPt[iCut]);
		}

		hESDMotherBackInvMassPt[iCut] = new TH2F("ESD_DalitzBackground_InvMass_Pt","ESD_DalitzBackground_InvMass_Pt",800,0,0.8,250,0,25);
		fESDList[iCut]->Add(hESDMotherBackInvMassPt[iCut]);

		fCutFolder[iCut]->Add(fESDList[iCut]);		
	}


	InitBack(); // Init Background Handler


	if( fIsMC ){
		// MC Histogramms
		fMCList = new TList*[fnCuts];
		// True Histogramms
		fTrueList = new TList*[fnCuts];
		hESDTrueConvGammaPt = new TH1F*[fnCuts];
		fHistoDoubleCountTrueConvGammaRPt = new TH2F*[fnCuts];
		hESDTruePositronPt  = new TH1F*[fnCuts];
		hESDTrueElectronPt  = new TH1F*[fnCuts];
		hESDTrueSecConvGammaPt = new TH1F*[fnCuts];
		hESDTrueSecPositronPt  = new TH1F*[fnCuts];
		hESDTrueSecElectronPt  = new TH1F*[fnCuts];
		hESDTruePi0DalitzConvGammaPt = new TH1F*[fnCuts];
		
		fHistoDoubleCountTruePi0InvMassPt = new TH2F*[fnCuts];
		fHistoDoubleCountTrueEtaInvMassPt = new TH2F*[fnCuts];
		hESDTruePi0DalitzPositronPt  = new TH1F*[fnCuts];
		hESDTruePi0DalitzPositronPtMB = new TH1F*[fnCuts];
		hESDTruePi0DalitzElectronPt  = new TH1F*[fnCuts];
		hESDTruePi0DalitzElectronPtMB = new TH1F*[fnCuts];
		hESDTruePi0DalitzSecConvGammaPt = new TH1F*[fnCuts];
		hESDTruePi0DalitzSecPositronPt  = new TH1F*[fnCuts];
		hESDTruePi0DalitzSecElectronPt  = new TH1F*[fnCuts];

		hMCAllGammaPt  		= new TH1F*[fnCuts];
		hMCAllGammaPi0Pt        = new TH1F*[fnCuts];
		hMCConvGammaPt 		= new TH1F*[fnCuts];
		hMCConvGammaRSPt 	= new TH1F*[fnCuts];
		hMCConvGammaPi0Pt       = new TH1F*[fnCuts];
		hMCAllPositronsPt 	= new TH1F*[fnCuts];
		hMCDecayPositronPi0Pt   = new TH1F*[fnCuts];
		hMCAllElectronsPt 	= new TH1F*[fnCuts];
		hMCDecayElectronPi0Pt   = new TH1F*[fnCuts];
		
		if( fDoMesonQA > 0 ) {

		hMCConvGammaEta    = new TH1F*[fnCuts];
		hMCConvGammaR 	   = new TH1F*[fnCuts];
		hMCAllPositronsEta = new TH1F*[fnCuts];
		hMCAllElectronsEta = new TH1F*[fnCuts];
		hMCConvGammaPtR    = new TH2F*[fnCuts];
		
		      		
		}
		hMCPi0DalitzGammaPt    = new TH1F*[fnCuts];
		hMCPi0DalitzElectronPt = new TH1F*[fnCuts];
		hMCPi0DalitzPositronPt = new TH1F*[fnCuts];
	
		hMCPi0Pt = new TH1F*[fnCuts];
		hMCPi0WOWeightPt = new TH1F*[fnCuts];
		hMCPi0GGPt =  new TH1F*[fnCuts];
		hMCEtaPt = new TH1F*[fnCuts];
		hMCEtaGGPt = new TH1F*[fnCuts];
		hMCPi0InAccPt = new TH1F*[fnCuts];
		hMCPi0WOWeightInAccPt = new TH1F*[fnCuts];
		hMCPi0InAccOpeningAngleGammaElectron = new TH1F*[fnCuts];
		hMCEtaInAccPt = new TH1F*[fnCuts];
		hMCChiCPt = new TH1F*[fnCuts];
		hMCChiCInAccPt = new TH1F*[fnCuts];
				
				
		if ( fDoMesonQA > 0 ) {
		
		hMCPi0EposEnegInvMassPt			     = new TH2F*[fnCuts];	
		hMCEtaEposEnegInvMassPt			     = new TH2F*[fnCuts];
		hESDEposEnegTruePi0DalitzInvMassPt           = new TH2F*[fnCuts];
		hESDEposEnegTruePrimPi0DalitzInvMass         = new TH1F*[fnCuts];
		hESDEposEnegTruePi0DalitzPsiPairDPhi         = new TH2F*[fnCuts];
		hESDEposEnegTruePi0DalitzPsiPairMC	     = new TH1F*[fnCuts];
		hESDEposEnegTruePi0DalitzPsiPairEta	     = new TH2F*[fnCuts];
		hESDEposEnegTruePi0DalitzDPhiEta	     = new TH2F*[fnCuts];
		hESDEposEnegTrueEtaDalitzInvMassPt           = new TH2F*[fnCuts];
		hESDEposEnegTruePrimEtaDalitzInvMass         = new TH1F*[fnCuts];
		hESDEposEnegTrueEtaDalitzPsiPairDPhi         = new TH2F*[fnCuts];
		hESDEposEnegTruePhotonInvMassPt              = new TH2F*[fnCuts];
		hESDEposEnegTrueInvMassPt                    = new TH2F*[fnCuts];
		hESDEposEnegTrueMotherInvMassPt		     = new TH2F*[fnCuts];
		hESDEposEnegTruePhotonPsiPairDPhi            = new TH2F*[fnCuts];
		hESDEposEnegTruePhotonPsiPairDPhiPtCut       = new TH2F*[fnCuts];
		hESDEposEnegTrueJPsiInvMassPt                = new TH2F*[fnCuts];
		hESDTrueConvGammaR                           = new TH1F*[fnCuts];
		hESDTrueConvGammaRMC                         = new TH1F*[fnCuts];
		hESDTrueConvGammaPtMC                        = new TH1F*[fnCuts];
		hESDTruePi0PtY 				     = new TH2F*[fnCuts];
		hESDTruePi0PtAlpha 			     = new TH2F*[fnCuts];
		hESDTruePi0PtOpenAngle 			     = new TH2F*[fnCuts];
		
		    if ( fDoMesonQA > 1 ) {
			  hESDTruePi0DalitzConvGammaR  	     = new TH1F*[fnCuts];
			  sESDTruePi0DalitzPlot  	     = new THnSparseF*[fnCuts];
			  sMCPi0DalitzPlot                   = new THnSparseF*[fnCuts];
		    }
		
		
		}
		
		
		if( fDoChicAnalysis ){
		hESDTrueMotherChiCInvMassPt = new TH2F*[fnCuts];
		hESDTrueMotherChiCDiffInvMassPt = new TH2F*[fnCuts];
		}
		
		
		hESDTrueMotherInvMassPt = new TH2F*[fnCuts];
		hESDTrueMotherW0WeightsInvMassPt = new TH2F*[fnCuts];
		hESDTrueMotherDalitzInvMassPt = new TH2F*[fnCuts];
		hESDTrueMotherPi0GGInvMassPt = new TH2F*[fnCuts];
		hESDTrueMotherPi0GGW0WeightsInvMassPt = new TH2F*[fnCuts];
		hESDTruePrimaryMotherPi0GGInvMassPt = new TH2F*[fnCuts];
		hESDTrueSecondaryMotherPi0GGInvMassPt = new TH2F*[fnCuts];
		hESDTruePrimaryPi0DalitzESDPtMCPt = new TH2F*[fnCuts];
		hESDTruePrimaryMotherInvMassMCPt = new TH2F*[fnCuts];
		hESDTruePrimaryMotherInvMassPt   = new TH2F*[fnCuts];
		hESDTruePrimaryMotherW0WeightingInvMassPt = new TH2F*[fnCuts];
		hESDTrueSecondaryMotherInvMassPt = new TH2F*[fnCuts];
		hESDTrueSecondaryMotherFromK0sInvMassPt = new TH2F*[fnCuts];
		hESDTrueBckGGInvMassPt = new TH2F*[fnCuts];
		hESDTrueBckContInvMassPt = new TH2F*[fnCuts];
		hESDTrueMotherGGInvMassPt = new TH2F*[fnCuts];
		//}

		for(Int_t iCut = 0; iCut<fnCuts;iCut++){
			TString cutstringEvent 		= ((AliConvEventCuts*)fCutEventArray->At(iCut))->GetCutNumber();
			TString cutstringElectron 	= ((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->GetCutNumber();
			TString cutstringMeson		= ((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->GetCutNumber();
			TString cutstringGamma	 	= ((AliConversionPhotonCuts*)fCutGammaArray->At(iCut))->GetCutNumber();

			fMCList[iCut] = new TList();
			fMCList[iCut]->SetName(Form("%s_%s_%s_%s MC histograms",cutstringEvent.Data(),cutstringGamma.Data(),cutstringElectron.Data(),cutstringMeson.Data()));
			fMCList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fMCList[iCut]);

			hMCAllGammaPt[iCut] = new TH1F("MC_AllGamma_Pt","MC_AllGamma_Pt",250,0,25);
			fMCList[iCut]->Add(hMCAllGammaPt[iCut]);
			
			hMCAllGammaPi0Pt[iCut] = new TH1F("MC_AllGammaPi0_Pt","MC_AllGammaPi0_Pt",250,0,25);
			fMCList[iCut]->Add(hMCAllGammaPi0Pt[iCut]);
			
			hMCConvGammaPt[iCut] = new TH1F("MC_ConvGamma_Pt","MC_ConvGamma_Pt",250,0,25);
			fMCList[iCut]->Add(hMCConvGammaPt[iCut]);
			
			hMCConvGammaRSPt[iCut] = new TH1F("MC_ConvGamma_RS_Pt","MC_ConvGamma_RS_Pt",250,0,25);
			fMCList[iCut]->Add(hMCConvGammaRSPt[iCut]);
			
			hMCConvGammaPi0Pt[iCut] = new TH1F("MC_ConvGammaPi0_Pt","MC_ConvGammaPi0_Pt",250,0,25);
			fMCList[iCut]->Add(hMCConvGammaPi0Pt[iCut]);
										
			hMCAllPositronsPt[iCut] = new TH1F("MC_AllPositrons_Pt","MC_AllPositrons_Pt",1000,0,25);
			fMCList[iCut]->Add(hMCAllPositronsPt[iCut]);
			
			hMCDecayPositronPi0Pt[iCut] = new TH1F("MC_DecayPositronPi0_Pt","MC_DecayPositronPi0_Pt",1000,0,25);
			fMCList[iCut]->Add(hMCDecayPositronPi0Pt[iCut]);
				
			hMCAllElectronsPt[iCut] = new TH1F("MC_AllElectrons_Pt","MC_AllElectrons_Pt",1000,0,25);
			fMCList[iCut]->Add(hMCAllElectronsPt[iCut]);
			
			hMCDecayElectronPi0Pt[iCut] = new TH1F("MC_DecayElectronPi0_Pt","MC_DecayElectronPi0_Pt",1000,0,25);
			fMCList[iCut]->Add(hMCDecayElectronPi0Pt[iCut]);
			
			hMCPi0DalitzGammaPt[iCut] = new TH1F("MC_Pi0DalitzGamma_Pt","MC_Pi0DalitzGamma_Pt",250,0,25);
			hMCPi0DalitzGammaPt[iCut]->Sumw2();
			fMCList[iCut]->Add(hMCPi0DalitzGammaPt[iCut]);
				
			if ( fDoMesonQA > 0 ){
				hMCConvGammaEta[iCut] = new TH1F("MC_ConvGamma_Eta","MC_ConvGamma_Eta",600,-1.5,1.5);
				fMCList[iCut]->Add(hMCConvGammaEta[iCut]);
				hMCConvGammaR[iCut] = new TH1F("MC_ConvGamma_R","MC_ConvGamma_R",800,0,200);
				fMCList[iCut]->Add(hMCConvGammaR[iCut]);
				hMCAllPositronsEta[iCut] = new TH1F("MC_AllPositrons_Eta","MC_AllPositrons_Eta",600,-1.5,1.5);
				fMCList[iCut]->Add(hMCAllPositronsEta[iCut]);
				hMCAllElectronsEta[iCut] = new TH1F("MC_AllElectrons_Eta","MC_AllElectrons_Eta",600,-1.5,1.5);
				fMCList[iCut]->Add(hMCAllElectronsEta[iCut]);
				hMCConvGammaPtR[iCut] = new TH2F("MC_ConvGamma_Pt_R","MC_ConvGamma_Pt_R",250,0,25,180,0.,180.0);
				fMCList[iCut]->Add(hMCConvGammaPtR[iCut]);

			}
				
			hMCPi0DalitzPositronPt[iCut] = new TH1F("MC_Pi0DalitzPositron_Pt","MC_Pi0DalitzPositron_Pt",1000,0,25);
			hMCPi0DalitzPositronPt[iCut]->Sumw2();
			fMCList[iCut]->Add(hMCPi0DalitzPositronPt[iCut]);
			
			hMCPi0DalitzElectronPt[iCut] = new TH1F("MC_Pi0DalitzElectron_Pt","MC_Pi0DalitzElectron_Pt",1000,0,25);
			hMCPi0DalitzElectronPt[iCut]->Sumw2();
			fMCList[iCut]->Add(hMCPi0DalitzElectronPt[iCut]);
			
			hMCPi0Pt[iCut] = new TH1F("MC_Pi0_Pt","MC_Pi0_Pt",250,0,25);
			hMCPi0Pt[iCut]->Sumw2();
			fMCList[iCut]->Add(hMCPi0Pt[iCut]);
			
			hMCPi0WOWeightPt[iCut] = new TH1F("MC_Pi0_WOWeights_Pt","MC_Pi0_WOWeights_Pt",250,0,25);
			hMCPi0WOWeightPt[iCut]->Sumw2();
			fMCList[iCut]->Add(hMCPi0WOWeightPt[iCut]);
			
			hMCPi0GGPt[iCut] = new TH1F("MC_Pi0_GG_Pt","MC_Pi0_GG_Pt",250,0,25);
			hMCPi0GGPt[iCut]->Sumw2();
			fMCList[iCut]->Add(hMCPi0GGPt[iCut]);
			
			hMCEtaPt[iCut] = new TH1F("MC_Eta_Pt","MC_Eta_Pt",250,0,25);
			hMCEtaPt[iCut]->Sumw2();
			fMCList[iCut]->Add(hMCEtaPt[iCut]);

			hMCEtaGGPt[iCut] = new TH1F("MC_Eta_GG_Pt","MC_Eta_GG_Pt",250,0,25);
			hMCEtaGGPt[iCut]->Sumw2();
			fMCList[iCut]->Add(hMCEtaGGPt[iCut]);
			
			hMCPi0InAccPt[iCut] = new TH1F("MC_Pi0DalitzInAcc_Pt","MC_Pi0DalitzInAcc_Pt",250,0,25);
			hMCPi0InAccPt[iCut]->Sumw2();
			fMCList[iCut]->Add(hMCPi0InAccPt[iCut]);
			
			hMCPi0WOWeightInAccPt[iCut] = new TH1F("MC_Pi0WOWeightInAcc_Pt","MC_Pi0WOWeightInAcc_Pt",250,0,25);
			hMCPi0WOWeightInAccPt[iCut]->Sumw2();
			fMCList[iCut]->Add(hMCPi0WOWeightInAccPt[iCut]);
			
			hMCPi0InAccOpeningAngleGammaElectron[iCut] = new TH1F("MC_Pi0InAcc_OpeningAngle_GammaElectron","MC_Pi0InAcc_OpeningAngle_GammaElectron",100,0,TMath::Pi());
			fMCList[iCut]->Add(hMCPi0InAccOpeningAngleGammaElectron[iCut]);
			
			

			hMCEtaInAccPt[iCut] = new TH1F("MC_EtaDalitzInAcc_Pt","MC_EtaDalitzInAcc_Pt",250,0,25);
			hMCEtaInAccPt[iCut]->Sumw2();
			fMCList[iCut]->Add(hMCEtaInAccPt[iCut]);
			
			hMCChiCPt[iCut] = new TH1F("MC_ChiC_Pt","MC_ChiC_Pt",250,0,25);
			fMCList[iCut]->Add(hMCChiCPt[iCut]);

			hMCChiCInAccPt[iCut] = new TH1F("MC_ChiCInAcc_Pt","MC_ChiCInAcc_Pt",250,0,25);
			fMCList[iCut]->Add(hMCChiCInAccPt[iCut]);
			
			if ( fDoMesonQA > 0 ) {
				hMCPi0EposEnegInvMassPt[iCut] = new TH2F("MC_Pi0EposEneg_InvMassPt","MC_Pi0EposEneg_InvMassPt",100,0.0,0.5,100,0.,10.);
				fMCList[iCut]->Add(hMCPi0EposEnegInvMassPt[iCut]);
				hMCEtaEposEnegInvMassPt[iCut] = new TH2F("MC_EtaEposEneg_InvMassPt","MC_EtaEposEneg_InvMassPt",140,0.,0.7,100,0.,10.);
				fMCList[iCut]->Add(hMCEtaEposEnegInvMassPt[iCut]);	
			}
				
			fTrueList[iCut] = new TList();
			fTrueList[iCut]->SetName(Form("%s_%s_%s_%s True histograms",cutstringEvent.Data(),cutstringGamma.Data(),cutstringElectron.Data(),cutstringMeson.Data()));
			fTrueList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fTrueList[iCut]);

			if ( fDoMesonQA > 0 ) {
				
			        hESDEposEnegTruePi0DalitzInvMassPt[iCut] = new TH2F("ESD_EposEneg_TruePi0Dalitz_InvMassPt","ESD_EposEneg_TruePi0Dalitz_InvMassPt",nBinsMassPair,xMinMassPair,xMaxMassPair,100,0.,10.);
				if(fDoHistoDalitzMassLog)SetLogBinningXTH2(hESDEposEnegTruePi0DalitzInvMassPt[iCut]);
				fTrueList[iCut]->Add(hESDEposEnegTruePi0DalitzInvMassPt[iCut]);
				
				hESDEposEnegTruePrimPi0DalitzInvMass[iCut] = new TH1F("ESD_EposEneg_TruePrimPi0Dalitz_InvMass","ESD_EposEneg_TruePrimPi0Dalitz_InvMass",nBinsMassPair,xMinMassPair,xMaxMassPair);
				fTrueList[iCut]->Add(hESDEposEnegTruePrimPi0DalitzInvMass[iCut]);
				
				hESDEposEnegTrueInvMassPt[iCut] = new TH2F("ESD_EposEneg_True_InvMassPt","ESD_EposEneg_True_InvMassPt",nBinsMassPair,xMinMassPair,xMaxMassPair,100,0.,10.);
				if(fDoHistoDalitzMassLog)SetLogBinningXTH2(hESDEposEnegTrueInvMassPt[iCut]);
				fTrueList[iCut]->Add(hESDEposEnegTrueInvMassPt[iCut]);
				
				hESDEposEnegTrueMotherInvMassPt[iCut] = new TH2F("ESD_EposEneg_TrueMother_InvMassPt","ESD_EposEneg_TrueMother_InvMassPt",nBinsMassPair,xMinMassPair,xMaxMassPair,100,0.,10.);
				if(fDoHistoDalitzMassLog)SetLogBinningXTH2(hESDEposEnegTrueMotherInvMassPt[iCut]);
				fTrueList[iCut]->Add(hESDEposEnegTrueMotherInvMassPt[iCut]);
				
				hESDEposEnegTrueEtaDalitzInvMassPt[iCut] = new TH2F("ESD_EposEneg_TrueEtaDalitz_InvMassPt","ESD_EposEneg_TrueEtaDalitz_InvMassPt",nBinsMassPair,xMinMassPair,xMaxMassPair,100,0.,10.);
				if(fDoHistoDalitzMassLog)SetLogBinningXTH2(hESDEposEnegTrueEtaDalitzInvMassPt[iCut]);
				fTrueList[iCut]->Add(hESDEposEnegTrueEtaDalitzInvMassPt[iCut]);
				
				hESDEposEnegTruePrimEtaDalitzInvMass[iCut] = new TH1F("ESD_EposEneg_TruePrimEtaDalitz_InvMass","ESD_EposEneg_TruePrimEtaDalitz_InvMass",nBinsMassPair,xMinMassPair,xMaxMassPair);			
				//if(fDoHistoDalitzMassLog)SetLogBinningXTH2(hESDEposEnegTruePrimEtaDalitzInvMass[iCut]);
				fTrueList[iCut]->Add(hESDEposEnegTruePrimEtaDalitzInvMass[iCut]);
				
				hESDEposEnegTruePhotonInvMassPt[iCut] = new TH2F("ESD_EposEneg_TruePhoton_InvMassPt","ESD_EposEneg_TruePhoton_InvMassPt",nBinsMassPair,xMinMassPair,xMaxMassPair,100,0.,10.);
				if(fDoHistoDalitzMassLog) SetLogBinningXTH2(hESDEposEnegTruePhotonInvMassPt[iCut]);
				fTrueList[iCut]->Add(hESDEposEnegTruePhotonInvMassPt[iCut]);
				
				hESDEposEnegTrueJPsiInvMassPt[iCut] = new TH2F("ESD_EposEneg_TrueJPsi_InvMassPt","ESD_EposEneg_TrueJPsi_InvMassPt",1000,0.,5.,100,0.,10.);
				fTrueList[iCut]->Add(hESDEposEnegTrueJPsiInvMassPt[iCut]);
				
				
				hESDEposEnegTruePi0DalitzPsiPairDPhi[iCut] = new TH2F("ESD_EposEneg_TruePi0Dalitz_PsiPair_DPhi","ESD_EposEneg_TruePi0Dalitz_PsiPair_DPhi", 100, -1.0,1.0,100,-1.0,1.0 );
				fTrueList[iCut]->Add(hESDEposEnegTruePi0DalitzPsiPairDPhi[iCut]);
				
				hESDEposEnegTruePi0DalitzPsiPairMC[iCut] = new TH1F("ESD_EposEneg_TruePi0Dalitz_PsiPairMC","ESD_EposEneg_TruePi0Dalitz_PsiPairMC", 100, -1.0,1.0);
				fTrueList[iCut]->Add(hESDEposEnegTruePi0DalitzPsiPairMC[iCut]); 
				
				hESDEposEnegTruePi0DalitzPsiPairEta[iCut] = new TH2F("ESD_EposEneg_TruePi0Dalitz_PsiPair_Eta","ESD_EposEneg_TruePi0Dalitz_PsiPair_Eta", 100, -1.0,1.0,600,-1.5,1.5 );
				fTrueList[iCut]->Add(hESDEposEnegTruePi0DalitzPsiPairEta[iCut]);
				
				hESDEposEnegTruePi0DalitzDPhiEta[iCut] = new TH2F("ESD_EposEneg_TruePi0Dalitz_DPhi_Eta","ESD_EposEneg_TruePi0Dalitz_DPhi_Eta", 100, -1.0,1.0,600,-1.5,1.5 );
				fTrueList[iCut]->Add(hESDEposEnegTruePi0DalitzDPhiEta[iCut]);
				
				
				hESDEposEnegTrueEtaDalitzPsiPairDPhi[iCut] = new TH2F("ESD_EposEneg_TrueEtaDalitz_PsiPair_DPhi","ESD_EposEneg_TrueEtaDalitz_PsiPair_DPhi", 100, -1.0,1.0,100,-1.0,1.0 );
				fTrueList[iCut]->Add(hESDEposEnegTrueEtaDalitzPsiPairDPhi[iCut]);
			
								
				hESDEposEnegTruePhotonPsiPairDPhi[iCut] = new TH2F("ESD_EposEneg_TruePhoton_PsiPair_DPhi","ESD_EposEneg_TruePhoton_PsiPair_DPhi", 100, -1.0,1.0,100,-1.0,1.0 );
				fTrueList[iCut]->Add(hESDEposEnegTruePhotonPsiPairDPhi[iCut]);
			
				hESDEposEnegTruePhotonPsiPairDPhiPtCut[iCut] = new TH2F("ESD_EposEneg_TruePhoton_PsiPair_DPhi_PtCut","ESD_EposEneg_TruePhoton_PsiPair_DPhi_PtCut", 100, -1.0,1.0,100,-1.0,1.0 );
				fTrueList[iCut]->Add(hESDEposEnegTruePhotonPsiPairDPhiPtCut[iCut]);
			
								
				hESDTrueConvGammaR[iCut] = new TH1F("ESD_TrueConvGamma_R","ESD_TrueConvGamma_R",800,0,200);
				fTrueList[iCut]->Add(hESDTrueConvGammaR[iCut]);
				
				hESDTrueConvGammaRMC[iCut] = new TH1F("ESD_TrueConvGamma_R_MC","ESD_TrueConvGamma_R_MC",800,0,200);
				fTrueList[iCut]->Add(hESDTrueConvGammaRMC[iCut]);
				
				hESDTrueConvGammaPtMC[iCut] = new TH1F("ESD_TrueConvGamma_Pt_MC","ESD_TrueConvGamma_Pt_MC",250,0,25);
				fTrueList[iCut]->Add(hESDTrueConvGammaPtMC[iCut]);
				
				hESDTruePi0EposEnegInvMassPi0Pt[iCut] = new TH2F("ESD_TruePi0EposEnegInvMass_Pi0Pt","ESD_TruePi0EposEnegInvMass_Pi0Pt",100,0.0,0.5,100,0.,10.);
				fTrueList[iCut]->Add(hESDTruePi0EposEnegInvMassPi0Pt[iCut]);
				
				
				hESDTruePi0PtY[iCut] = new TH2F("ESD_TruePi0_Pt_Y","ESD_TruePi0_Pt_Y",150,0.03,15.,150,-1.5,1.5);
				SetLogBinningXTH2(hESDTruePi0PtY[iCut]);
				fTrueList[iCut]->Add(hESDTruePi0PtY[iCut]);
				
				hESDTruePi0PtAlpha[iCut] = new TH2F("ESD_TruePi0_Pt_Alpha","ESD_TruePi0_Pt_Alpha",150,0.03,15.,100,0,1);
				SetLogBinningXTH2(hESDTruePi0PtAlpha[iCut]);
				fTrueList[iCut]->Add(hESDTruePi0PtAlpha[iCut]);
				
				hESDTruePi0PtOpenAngle[iCut] = new TH2F("ESD_TruePi0_Pt_OpenAngle","ESD_TruePi0_Pt_OpenAngle",150,0.03,15.,200,0,2*TMath::Pi());            
				SetLogBinningXTH2(hESDTruePi0PtOpenAngle[iCut]);
				fTrueList[iCut]->Add(hESDTruePi0PtOpenAngle[iCut]);
				
				if(fDoMesonQA > 1 ) {
				hESDTruePi0DalitzConvGammaR[iCut] = new TH1F("ESD_TruePi0DalitzConvGamma_R","hESDTruePi0DalitzConvGammaR",800,0,200);
				fTrueList[iCut]->Add(hESDTruePi0DalitzConvGammaR[iCut]);
				
				
				const Int_t nDimDalPlot = 2;
				Int_t    nBinsDalPlot[nDimDalPlot] = {  4000, 4000};
				Double_t xMinDalPlot[nDimDalPlot]  = {    0,     0};
				Double_t xMaxDalPlot[nDimDalPlot]  = {  4.0,   4.0};
			
			
				sESDTruePi0DalitzPlot[iCut] = new THnSparseF("ESD_TruePi0Dalitz_DalitzPlot","ESD_TruePi0Dalitz_DalitzPlot",nDimDalPlot,nBinsDalPlot,xMinDalPlot,xMaxDalPlot);
				fQAFolder[iCut]->Add(sESDTruePi0DalitzPlot[iCut]);
				
				sMCPi0DalitzPlot[iCut]      = new THnSparseF("MC_Pi0Dalitz_DalitzPlot","MC_Pi0Dalitz_DalitzPlot",nDimDalPlot,nBinsDalPlot,xMinDalPlot,xMaxDalPlot);
				fQAFolder[iCut]->Add(sMCPi0DalitzPlot[iCut]);
				
				}
				
			}

			hESDTruePositronPt[iCut] = new TH1F("ESD_TruePositron_Pt","ESD_TruePositron_Pt",1000,0,25);
			fTrueList[iCut]->Add(hESDTruePositronPt[iCut]);

			hESDTrueElectronPt[iCut] = new TH1F("ESD_TrueElectron_Pt","ESD_TrueElectron_Pt",1000,0,25);
			fTrueList[iCut]->Add(hESDTrueElectronPt[iCut]);
		
			hESDTrueSecPositronPt[iCut] = new TH1F("ESD_TrueSecPositron_Pt","ESD_TrueSecPositron_Pt",1000,0,25);
			fTrueList[iCut]->Add(hESDTrueSecPositronPt[iCut]);

			hESDTrueSecElectronPt[iCut] = new TH1F("ESD_TrueSecElectron_Pt","ESD_TrueSecElectron_Pt",1000,0,25);
			fTrueList[iCut]->Add(hESDTrueSecElectronPt[iCut]); 

			hESDTrueConvGammaPt[iCut] = new TH1F("ESD_TrueConvGamma_Pt","ESD_TrueConvGamma_Pt",250,0,25);
			fTrueList[iCut]->Add(hESDTrueConvGammaPt[iCut]);

			fHistoDoubleCountTrueConvGammaRPt[iCut] = new TH2F("ESD_TrueDoubleCountConvGamma_R_Pt","ESD_TrueDoubleCountConvGamma_R_Pt",800,0,200,300,0,30);
			fTrueList[iCut]->Add(fHistoDoubleCountTrueConvGammaRPt[iCut]);
		
			hESDTrueSecConvGammaPt[iCut] = new TH1F("ESD_TrueSecConvGamma_Pt","ESD_TrueSecConvGamma_Pt",250,0,25);
			fTrueList[iCut]->Add(hESDTrueSecConvGammaPt[iCut]);

			hESDTruePi0DalitzConvGammaPt[iCut] = new TH1F("ESD_TruePi0DalitzConvGamma_Pt","ESD_TruePi0DalitzConvGamma_Pt",250,0,25);
			fTrueList[iCut]->Add(hESDTruePi0DalitzConvGammaPt[iCut]);
			
				

			fHistoDoubleCountTruePi0InvMassPt[iCut] = new TH2F("ESD_TrueDoubleCountPi0_InvMass_Pt","ESD_TrueDoubleCountPi0_InvMass_Pt",800,0,0.8,300,0,30);
			fTrueList[iCut]->Add(fHistoDoubleCountTruePi0InvMassPt[iCut]);
			fHistoDoubleCountTrueEtaInvMassPt[iCut] = new TH2F("ESD_TrueDoubleCountEta_InvMass_Pt","ESD_TrueDoubleCountEta_InvMass_Pt",800,0,0.8,300,0,30);
			fTrueList[iCut]->Add(fHistoDoubleCountTrueEtaInvMassPt[iCut]);

			hESDTruePi0DalitzElectronPt[iCut] = new TH1F("ESD_TruePi0DalitzElectron_Pt","ESD_TruePi0DalitzElectron_Pt",1000,0,25);
			fTrueList[iCut]->Add(hESDTruePi0DalitzElectronPt[iCut]);
			
			hESDTruePi0DalitzElectronPtMB[iCut] = new TH1F("ESD_TruePi0DalitzElectron_Pt_MB","ESD_TruePi0DalitzElectron_Pt_MB",1000,0,25);
			fTrueList[iCut]->Add(hESDTruePi0DalitzElectronPtMB[iCut]);

			hESDTruePi0DalitzPositronPt[iCut] = new TH1F("ESD_TruePi0DalitzPositron_Pt","ESD_TruePi0DalitzPositron_Pt",1000,0,25);
			fTrueList[iCut]->Add(hESDTruePi0DalitzPositronPt[iCut]);
			
			hESDTruePi0DalitzPositronPtMB[iCut] = new TH1F("ESD_TruePi0DalitzPositron_Pt_MB","ESD_TruePi0DalitzPositron_Pt_MB",1000,0,25);
			fTrueList[iCut]->Add(hESDTruePi0DalitzPositronPtMB[iCut]);
			
		
			hESDTruePi0DalitzSecConvGammaPt[iCut] = new TH1F("ESD_TruePi0DalitzSecConvGamma_Pt","ESD_TruePi0DalitzSecConvGamma_Pt",250,0,25);
			fTrueList[iCut]->Add(hESDTruePi0DalitzSecConvGammaPt[iCut]);
		
			hESDTruePi0DalitzSecElectronPt[iCut] = new TH1F("ESD_TruePi0DalitzSecElectron_Pt","ESD_TruePi0DalitzSecElectron_Pt",1000,0,25);
			fTrueList[iCut]->Add(hESDTruePi0DalitzSecElectronPt[iCut]);

			hESDTruePi0DalitzSecPositronPt[iCut] = new TH1F("ESD_TruePi0DalitzSecPositron_Pt","ESD_TruePi0DalitzSecPositron_Pt",1000,0,25);
			fTrueList[iCut]->Add(hESDTruePi0DalitzSecPositronPt[iCut]);

			if( fDoChicAnalysis) { 
				hESDTrueMotherChiCInvMassPt[iCut] = new TH2F("ESD_TrueMotherChiC_InvMass_Pt","ESD_TrueMotherChiC_InvMass_Pt",4000,0,4,250,0,25);
				fTrueList[iCut]->Add(hESDTrueMotherChiCInvMassPt[iCut]);
				hESDTrueMotherChiCDiffInvMassPt[iCut] = new TH2F("ESD_TrueMotherChiCDiff_InvMass_Pt","ESD_TrueMotherChiCDiff_InvMass_Pt",2000,0,2,250,0,25);
				fTrueList[iCut]->Add(hESDTrueMotherChiCDiffInvMassPt[iCut]);
			}

			hESDTrueMotherInvMassPt[iCut] = new TH2F("ESD_TrueMother_InvMass_Pt","ESD_TrueMother_InvMass_Pt",800,0,0.8,250,0,25);
			hESDTrueMotherInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(hESDTrueMotherInvMassPt[iCut]);
			
			hESDTrueMotherW0WeightsInvMassPt[iCut] = new TH2F("ESD_TrueMotherW0Weights_InvMass_Pt","ESD_TrueMotherW0Weights_InvMass_Pt",800,0,0.8,250,0,25);
			hESDTrueMotherW0WeightsInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(hESDTrueMotherW0WeightsInvMassPt[iCut]);
			
			
		
			hESDTrueMotherDalitzInvMassPt[iCut] = new TH2F("ESD_TrueMother_Dalitz_InvMass_Pt","ESD_TrueMother_Dalitz_InvMass_Pt",800,0,0.8,250,0,25);
			hESDTrueMotherDalitzInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(hESDTrueMotherDalitzInvMassPt[iCut]);
		
		
		
		

			hESDTrueMotherPi0GGInvMassPt[iCut] = new TH2F("ESD_TrueMotherPi0GG_InvMass_Pt","ESD_TrueMotherPi0GG_InvMass_Pt",800,0,0.8,250,0,25);
			hESDTrueMotherPi0GGInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(hESDTrueMotherPi0GGInvMassPt[iCut]);
			
			hESDTrueMotherPi0GGW0WeightsInvMassPt[iCut] = new TH2F("ESD_TrueMotherPi0GGW0Weights_InvMass_Pt","ESD_TrueMotherPi0GGW0Weights_InvMass_Pt",800,0,0.8,250,0,25);
			hESDTrueMotherPi0GGW0WeightsInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(hESDTrueMotherPi0GGW0WeightsInvMassPt[iCut]);
			
			hESDTruePrimaryMotherPi0GGInvMassPt[iCut] = new TH2F("ESD_TruePrimaryMotherPi0GG_InvMass_Pt","ESD_TruePrimaryMotherPi0GG_InvMass_Pt",800,0,0.8,250,0,25);
			hESDTruePrimaryMotherPi0GGInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(hESDTruePrimaryMotherPi0GGInvMassPt[iCut]);

			hESDTrueSecondaryMotherPi0GGInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryMotherPi0GG_InvMass_Pt","ESD_TrueSecondaryMotherPi0GG_InvMass_Pt",800,0,0.8,250,0,25);
			hESDTrueSecondaryMotherPi0GGInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(hESDTrueSecondaryMotherPi0GGInvMassPt[iCut]);

			hESDTruePrimaryPi0DalitzESDPtMCPt[iCut] = new TH2F("ESD_TruePrimaryPi0Dalitz_ESDPt_MCPt","ESD_TruePrimaryPi0Dalitz_ESDPt_MCPt",250,0,25,250,0,25);
			hESDTruePrimaryPi0DalitzESDPtMCPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(hESDTruePrimaryPi0DalitzESDPtMCPt[iCut]);
			hESDTruePrimaryMotherInvMassMCPt[iCut] = new TH2F("ESD_TruePrimaryMother_InvMass_MCPt","ESD_TrueDalitzPrimaryMother_InvMass_MCPt",800,0,0.8,250,0,25);
			hESDTruePrimaryMotherInvMassMCPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(hESDTruePrimaryMotherInvMassMCPt[iCut]);
			hESDTruePrimaryMotherInvMassPt[iCut] = new TH2F("ESD_TruePrimaryMother_InvMass_Pt","ESD_TruePrimaryMother_InvMass_Pt",800,0,0.8,250,0,25);
			hESDTruePrimaryMotherInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(hESDTruePrimaryMotherInvMassPt[iCut]);
			hESDTruePrimaryMotherW0WeightingInvMassPt[iCut] = new TH2F("ESD_TruePrimaryMotherW0Weighting_InvMass_Pt","ESD_TruePrimaryMotherW0Weighting_InvMass_Pt",800,0,0.8,250,0,25);
			hESDTruePrimaryMotherW0WeightingInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(hESDTruePrimaryMotherW0WeightingInvMassPt[iCut]); 

			hESDTrueSecondaryMotherInvMassPt[iCut] = new TH2F("ESD_TrueDalitzSecondaryMother_InvMass_Pt","ESD_TrueDalitzSecondaryMother_InvMass_Pt",800,0,0.8,250,0,25);
			hESDTrueSecondaryMotherInvMassPt[iCut]->Sumw2();
			fTrueList[iCut]->Add(hESDTrueSecondaryMotherInvMassPt[iCut]);
			// 				hESDTrueSecondaryMotherFromK0sInvMassPt[iCut] = new TH2F("ESD_TrueDalitzSecondaryMotherFromK0s_InvMass_Pt","ESD_TrueDalitzSecondaryMotherFromK0s_InvMass_Pt",1000,0,1,250,0,25);
			// 				fTrueList[iCut]->Add(hESDTrueSecondaryMotherFromK0sInvMassPt[iCut]);
			hESDTrueBckGGInvMassPt[iCut] = new TH2F("ESD_TrueDalitzBckGG_InvMass_Pt","ESD_TrueDalitzBckGG_InvMass_Pt",800,0,0.8,250,0,25);
			fTrueList[iCut]->Add(hESDTrueBckGGInvMassPt[iCut]);
			hESDTrueBckContInvMassPt[iCut] = new TH2F("ESD_TrueDalitzBckCont_InvMass_Pt","ESD_TrueDalitzBckCont_InvMass_Pt",800,0,0.8,250,0,25);
			fTrueList[iCut]->Add(hESDTrueBckContInvMassPt[iCut]);
		// 				hESDTrueMotherGGInvMassPt[iCut] = new TH2F("ESD_TrueGammaGamma_InvMass_Pt","ESD_TrueGammaGamma_InvMass_Pt",1000,0,1,250,0,25);
		// 				fTrueList[iCut]->Add(hESDTrueMotherGGInvMassPt[iCut]);

		}
	}

	fVectorDoubleCountTruePi0s.clear();
	fVectorDoubleCountTrueEtas.clear();
	fVectorDoubleCountTrueConvGammas.clear();

    fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
	if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader


	if(fV0Reader)
		if((AliConvEventCuts*)fV0Reader->GetEventCuts())
			if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
				fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());
	
	if(fV0Reader)
		if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
			if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
				fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
			
	if(fV0Reader && fV0Reader->GetProduceV0FindingEfficiency())
                if (fV0Reader->GetV0FindingEfficiencyHistograms())
                        fOutputContainer->Add(fV0Reader->GetV0FindingEfficiencyHistograms());


	
	
	
	fElecSelector=(AliDalitzElectronSelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("ElectronSelector");
	if(!fElecSelector){printf("Error: No ElectronSelector");return;} // GetV0Reader
		
	if( fElecSelector ){
		if ( ((AliDalitzElectronCuts*)fElecSelector->GetDalitzElectronCuts())->GetCutHistograms() ){
				fOutputContainer->Add( ((AliDalitzElectronCuts*)fElecSelector->GetDalitzElectronCuts())->GetCutHistograms() );
		}
	}  

	for(Int_t iCut = 0; iCut<fnCuts;iCut++){

		if( fCutElectronArray ){
			if( ((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->GetCutHistograms() ) {
				fCutFolder[iCut]->Add( ((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->GetCutHistograms() );
			}
		}

		if( fCutMesonArray  ) {
			if( ((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->GetCutHistograms() ) {
				fCutFolder[iCut]->Add( ((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->GetCutHistograms());
			}
		}

		if( fCutGammaArray ) {
			if( ((AliConversionPhotonCuts*)fCutGammaArray->At(iCut))->GetCutHistograms() ) {
				fCutFolder[iCut]->Add( ((AliConversionPhotonCuts*)fCutGammaArray->At(iCut))->GetCutHistograms()  );
			}
		}
		
		if( fCutEventArray )
		if( ((AliConvEventCuts*)fCutEventArray->At(iCut))->GetCutHistograms()){
			fCutFolder[iCut]->Add(((AliConvEventCuts*)fCutEventArray->At(iCut))->GetCutHistograms());
		}
	}

	PostData(1, fOutputContainer);

}

//______________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::UserExec(Option_t *)
{

  //
  // Execute analysis for current event
  //

  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader


  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();

  // Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1
  if(eventQuality == 2 || eventQuality == 3){
    for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      hNEvents[iCut]->Fill(eventQuality);
    }
    return;
  }

  fElecSelector=(AliDalitzElectronSelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("ElectronSelector");
  if(!fElecSelector){printf("Error: No ElectronSelector");return;} // GetV0Reader


  if(fIsMC) 	      fMCEvent        =  MCEvent();
  fESDEvent        = (AliESDEvent*)InputEvent();
  fReaderGammas    = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut
  
  AliEventplane *EventPlane = fInputEvent->GetEventplane();
  if(fIsHeavyIon == 1)fEventPlaneAngle = EventPlane->GetEventplane("V0",fInputEvent,2);
  else fEventPlaneAngle=0.0;
  
  
  fSelectorElectronIndex = fElecSelector->GetReconstructedElectronsIndex(); // Electrons from default Cut
  fSelectorPositronIndex = fElecSelector->GetReconstructedPositronsIndex(); // Positrons from default Cut

  CountESDTracks(); // Estimate Event Multiplicity
  fNumberOfESDTracks = fV0Reader->GetNumberOfPrimaryTracks();
  //AddTaskContainers(); //Add conatiner

  for(Int_t iCut = 0; iCut<fnCuts; iCut++){		
    fiCut = iCut;
    fNVirtualGammas = 0;

    Int_t eventNotAccepted =
      ((AliConvEventCuts*)fCutEventArray->At(iCut))
      ->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);
    
    
    //Int_t eventNotAccepted = ((AliConvEventCuts*)fCutGammaArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);
    
    if(eventNotAccepted){
      // 			cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
      hNEvents[iCut]->Fill(eventNotAccepted); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
      continue;
    }

    if(eventQuality != 0){// Event Not Accepted
      // 			cout << "event rejected due to: " <<eventQuality << endl;
      hNEvents[iCut]->Fill(eventQuality);
      continue;
    }

    hNEvents[iCut]->Fill(eventQuality);
    hNGoodESDTracks[iCut]->Fill(fNumberOfESDTracks);
    
    if(((AliConvEventCuts*)fCutEventArray->At(iCut))->IsHeavyIon() == 2) hNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A());
    else hNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C());

    if(fMCEvent){ // Process MC Particle
      if(((AliConvEventCuts*)fCutEventArray->At(iCut))->GetSignalRejection() != 0){
          ((AliConvEventCuts*)fCutEventArray->At(iCut))->GetNotRejectedParticles( ((AliConvEventCuts*)fCutEventArray->At(iCut))->GetSignalRejection(),
                                              ((AliConvEventCuts*)fCutEventArray->At(iCut))->GetAcceptedHeader(),
                                              fMCEvent);
      } 
      ProcessMCParticles();
    }

    ProcessPhotonCandidates(); // Process this cuts gammas
    ProcessElectronCandidates(); // Process this cuts gammas
    
    if(((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->UseMCPSmearing() && fMCEvent){
    
      fUnsmearedPx = new Double_t[fGoodGammas->GetEntries()]; // Store unsmeared Momenta
      fUnsmearedPy = new Double_t[fGoodGammas->GetEntries()];
      fUnsmearedPz = new Double_t[fGoodGammas->GetEntries()];
      fUnsmearedE =  new Double_t[fGoodGammas->GetEntries()];

      for(Int_t gamma=0;gamma<fGoodGammas->GetEntries();gamma++){ // Smear the AODPhotons in MC
        fUnsmearedPx[gamma] = ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->Px();
        fUnsmearedPy[gamma] = ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->Py();
        fUnsmearedPz[gamma] = ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->Pz();
        fUnsmearedE[gamma] =  ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->E();
        ((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->SmearParticle(dynamic_cast<AliAODConversionPhoton*>(fGoodGammas->At(gamma)));
      }           
    }
    
    if( ((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->GetUseVPhotonMCPmearing() && ((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->UseMCPSmearing() && fMCEvent){
      // cout<<"Entro virtual photon smearing"<<endl;
      fUnsmearedVPx = new Double_t[fGoodVirtualGammas->GetEntries()]; // Store unsmeared Momenta
      fUnsmearedVPy = new Double_t[fGoodVirtualGammas->GetEntries()];
      fUnsmearedVPz = new Double_t[fGoodVirtualGammas->GetEntries()];
      fUnsmearedVE =  new Double_t[fGoodVirtualGammas->GetEntries()];

      for(Int_t Vgamma=0;Vgamma<fGoodVirtualGammas->GetEntries();Vgamma++){ // Smear the AODPhotons in MC
        fUnsmearedVPx[Vgamma] = ((AliAODConversionPhoton*)fGoodVirtualGammas->At(Vgamma))->Px();
        fUnsmearedVPy[Vgamma] = ((AliAODConversionPhoton*)fGoodVirtualGammas->At(Vgamma))->Py();
        fUnsmearedVPz[Vgamma] = ((AliAODConversionPhoton*)fGoodVirtualGammas->At(Vgamma))->Pz();
        fUnsmearedVE[Vgamma] =  ((AliAODConversionPhoton*)fGoodVirtualGammas->At(Vgamma))->E();
        ((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->SmearVirtualPhoton(dynamic_cast<AliAODConversionPhoton*>(fGoodVirtualGammas->At(Vgamma)));
      }           
    }
    
    ProcessVirtualGammasCandidates();
    CalculatePi0DalitzCandidates();
    
    if(((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->DoBGCalculation()){
        if(((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->BackgroundHandlerType() == 0){
      
    
    CalculateBackground();
    UpdateEventByEventData();
    
        }
    }
        
    
    if ( fDoMesonQA > 0) {
      hNGoodESDTracksVsNGoodGammas[iCut]->Fill(fNumberOfESDTrackskBoth,fGoodGammas->GetEntries());
      hNGoodESDTracksVsNGoodVGammas[iCut]->Fill(fNumberOfESDTrackskBoth,fNVirtualGammas);
      fHistoSPDClusterTrackletBackground[iCut]->Fill(fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),(fInputEvent->GetNumberOfITSClusters(0)+fInputEvent->GetNumberOfITSClusters(1)));
    }
    
    
    if(((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->UseMCPSmearing() && fMCEvent){
    
      for(Int_t gamma=0;gamma<fGoodGammas->GetEntries();gamma++){ // Smear the AODPhotons in MC
        ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->SetPx(fUnsmearedPx[gamma]); // Reset Unsmeared Momenta
        ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->SetPy(fUnsmearedPy[gamma]);
        ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->SetPz(fUnsmearedPz[gamma]);
        ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->SetE(fUnsmearedE[gamma]);
      }
      delete[] fUnsmearedPx; fUnsmearedPx = 0x0;
      delete[] fUnsmearedPy; fUnsmearedPy = 0x0;
      delete[] fUnsmearedPz; fUnsmearedPz = 0x0;
      delete[] fUnsmearedE;  fUnsmearedE  = 0x0;
    }
    
    if( ((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->GetUseVPhotonMCPmearing() && ((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->UseMCPSmearing() && fMCEvent){
    
      for(Int_t Vgamma=0;Vgamma<fGoodVirtualGammas->GetEntries();Vgamma++){ // Smear the AODPhotons in MC
        ((AliAODConversionPhoton*)fGoodVirtualGammas->At(Vgamma))->SetPx(fUnsmearedVPx[Vgamma]); // Reset Unsmeared Momenta
        ((AliAODConversionPhoton*)fGoodVirtualGammas->At(Vgamma))->SetPy(fUnsmearedVPy[Vgamma]);
        ((AliAODConversionPhoton*)fGoodVirtualGammas->At(Vgamma))->SetPz(fUnsmearedVPz[Vgamma]);
        ((AliAODConversionPhoton*)fGoodVirtualGammas->At(Vgamma))->SetE(fUnsmearedVE[Vgamma]);
      } 
      delete[] fUnsmearedVPx; fUnsmearedVPx = 0x0;
      delete[] fUnsmearedVPy; fUnsmearedVPy = 0x0;
      delete[] fUnsmearedVPz; fUnsmearedVPz = 0x0;
      delete[] fUnsmearedVE;  fUnsmearedVE  = 0x0;
    }
    fVectorDoubleCountTruePi0s.clear();
    fVectorDoubleCountTrueEtas.clear();
    fVectorDoubleCountTrueConvGammas.clear();
    fGoodGammas->Clear(); // delete this cuts good gammas
    fGoodVirtualGammas->Clear(); // delete this cuts good gammas
  }

  fSelectorElectronIndex.clear();
  fSelectorPositronIndex.clear();

  PostData( 1, fOutputContainer );
}

Bool_t AliAnalysisTaskGammaConvDalitzV1::Notify()
{
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if (((AliConvEventCuts*)fCutEventArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){        
        ((AliConvEventCuts*)fCutEventArray->At(iCut))->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (((AliConvEventCuts*)fCutEventArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      ((AliConvEventCuts*)fCutEventArray->At(iCut))->SetPeriodEnum(fV0Reader->GetPeriodName());
    }  
    
    
    if( !((AliConvEventCuts*)fCutEventArray->At(iCut))->GetDoEtaShift() ){

      hEtaShift[iCut]->Fill(0.,0.);
      continue; // No Eta Shift requested, continue
    }
    
    
    if( ((AliConvEventCuts*)fCutEventArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
      ((AliConvEventCuts*)fCutEventArray->At(iCut))->GetCorrectEtaShiftFromPeriod();
      ((AliConvEventCuts*)fCutEventArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once   
      hEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fCutEventArray->At(iCut))->GetEtaShift()));
      continue;
    } else {
      printf(" Gamma Conversion Dalitz Task %s :: Eta Shift Manually Set to %f \n\n",
      (((AliConvEventCuts*)fCutEventArray->At(iCut))->GetCutNumber()).Data(),((AliConvEventCuts*)fCutEventArray->At(iCut))->GetEtaShift());
      ((AliConvEventCuts*)fCutEventArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once   
      hEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fCutEventArray->At(iCut))->GetEtaShift()));
    }
  }
  
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::Terminate(const Option_t *)
{
	///Grid
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::ProcessPhotonCandidates()
{
	Int_t nV0 = 0;
	TList *GoodGammasStepOne = new TList();
	TList *GoodGammasStepTwo = new TList();
	// Loop over Photon Candidates allocated by ReaderV1
	
	for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
		AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
		if(!PhotonCandidate) continue;
		fIsFromMBHeader = kTRUE;
		if( fMCEvent && ((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetSignalRejection() != 0 ){
            Int_t isPosFromMBHeader = ((AliConvEventCuts*)fCutEventArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent, fInputEvent);
			if(isPosFromMBHeader == 0 && ((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetSignalRejection() != 3) continue;
            Int_t isNegFromMBHeader = ((AliConvEventCuts*)fCutEventArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent,fInputEvent);
			if(isNegFromMBHeader == 0 && ((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetSignalRejection() != 3) continue;
			if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
		}
		
		if(!((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fESDEvent)) continue;

		if(!((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->UseElecSharingCut() &&
			!((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->UseToCloseV0sCut()){ // if no post reader loop is required add to events good gammas
			
			fGoodGammas->Add(PhotonCandidate);
			if(fIsFromMBHeader){
				hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
				hESDConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
				if( fDoMesonQA > 1 ) {
				  
				  Double_t sparesFillZR[2] = {PhotonCandidate->GetConversionZ(),PhotonCandidate->GetConversionRadius()};										  
				  sESDConvGammaZR[fiCut]->Fill(sparesFillZR);
				  Double_t sparesFillXY[2] = {PhotonCandidate->GetConversionX(),PhotonCandidate->GetConversionY()};
				  sESDConvGammaXY[fiCut]->Fill(sparesFillXY);
				
				  
				  
				}
			}
			if(fMCEvent){
				ProcessTruePhotonCandidates(PhotonCandidate);
			}
		} else if(((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
			((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->FillElectonLabelArray(PhotonCandidate,nV0);
			nV0++;
			GoodGammasStepOne->Add(PhotonCandidate);
		} else if(!((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->UseElecSharingCut() &&
				((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
			GoodGammasStepTwo->Add(PhotonCandidate);
		}
	}
	
	
	if(((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->UseElecSharingCut()){
		for(Int_t i = 0;i<GoodGammasStepOne->GetEntries();i++){
			AliAODConversionPhoton *PhotonCandidate= (AliAODConversionPhoton*) GoodGammasStepOne->At(i);
			if(!PhotonCandidate) continue;	
			fIsFromMBHeader = kTRUE;
			if(fMCEvent && ((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetSignalRejection() != 0){
                Int_t isPosFromMBHeader = ((AliConvEventCuts*)fCutEventArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent,fInputEvent);
                Int_t isNegFromMBHeader = ((AliConvEventCuts*)fCutEventArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent,fInputEvent);
				if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
			}
			
			if(!((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GoodGammasStepOne->GetEntries())) continue;
			if(!((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
				fGoodGammas->Add(PhotonCandidate);
			
				if(fIsFromMBHeader){
					hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
					hESDConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
					if( fDoMesonQA > 1 ){
					    Double_t sparesFillZR[2] = {PhotonCandidate->GetConversionZ(),PhotonCandidate->GetConversionRadius()};										  
					    sESDConvGammaZR[fiCut]->Fill(sparesFillZR);
					    Double_t sparesFillXY[2] = {PhotonCandidate->GetConversionX(),PhotonCandidate->GetConversionY()};
					    sESDConvGammaXY[fiCut]->Fill(sparesFillXY);
					  }
				}
			
				if(fMCEvent){
					ProcessTruePhotonCandidates(PhotonCandidate);
				}
			} else GoodGammasStepTwo->Add(PhotonCandidate); // Close v0s cut enabled -> add to list two
		}
	}
	if(((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->UseToCloseV0sCut()){
		for(Int_t i = 0;i<GoodGammasStepTwo->GetEntries();i++){
			AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GoodGammasStepTwo->At(i);
			if(!PhotonCandidate) continue;
			if(fMCEvent && ((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetSignalRejection() != 0){
                Int_t isPosFromMBHeader = ((AliConvEventCuts*)fCutEventArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCEvent,fInputEvent);
                Int_t isNegFromMBHeader = ((AliConvEventCuts*)fCutEventArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCEvent,fInputEvent);
				if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
			}
			if(!((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GoodGammasStepTwo,i)) continue;
			fGoodGammas->Add(PhotonCandidate); // Add gamma to current cut TList
		
			if(fIsFromMBHeader){
				hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt()); // Differences to old V0Reader in p_t due to conversion KF->TLorentzVector
				hESDConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
				if(fDoMesonQA > 1 ){
				   Double_t sparesFillZR[2] = {PhotonCandidate->GetConversionZ(),PhotonCandidate->GetConversionRadius()};										  
				   sESDConvGammaZR[fiCut]->Fill(sparesFillZR);
				   Double_t sparesFillXY[2] = {PhotonCandidate->GetConversionX(),PhotonCandidate->GetConversionY()};
				   sESDConvGammaXY[fiCut]->Fill(sparesFillXY);
				}
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
void AliAnalysisTaskGammaConvDalitzV1::ProcessTruePhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
	// Process True Photons
    TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(fMCEvent);
    TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCEvent);

	if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
	if(posDaughter->GetMother(0) != negDaughter->GetMother(0)){  // Not Same Mother == Combinatorial Bck
		return;
	} else if (posDaughter->GetMother(0) == -1){
		return;
	}
	
	if(TMath::Abs(posDaughter->GetPdgCode())!=11 || TMath::Abs(negDaughter->GetPdgCode())!=11) return; //One Particle is not electron
	if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge
	if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() !=5) return;// check if the daughters come from a conversion

    TParticle *Photon = TruePhotonCandidate->GetMCParticle(fMCEvent);
	if(Photon->GetPdgCode() != 22) return; // Mother is no Photon

	// True Photon
	
    Int_t labelGamma = TruePhotonCandidate->GetMCParticleLabel(fMCEvent);
	
    if( labelGamma < fMCEvent->GetNumberOfPrimaries() ){
		if( fIsFromMBHeader ){
			hESDTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			if (CheckVectorForDoubleCount(fVectorDoubleCountTrueConvGammas,posDaughter->GetMother(0))) fHistoDoubleCountTrueConvGammaRPt[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius(),TruePhotonCandidate->Pt());
			if( fDoMesonQA > 0 ){
			  hESDTrueConvGammaR[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius());			 
			  hESDTrueConvGammaRMC[fiCut]->Fill( posDaughter->R() );
			  hESDTrueConvGammaPtMC[fiCut]->Fill( Photon->Pt());
			}
		}
	} else {
		if( fIsFromMBHeader){
			hESDTrueSecConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
		}
	}
	
	if( IsPi0DalitzDaughter(labelGamma) == kTRUE ) {
        if( labelGamma < fMCEvent->GetNumberOfPrimaries() ) {
			if( fIsFromMBHeader ){
				hESDTruePi0DalitzConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
				if( fDoMesonQA > 1 ){
				hESDTruePi0DalitzConvGammaR[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius());
				}
			}
		} else {
			if( fIsFromMBHeader ) {
				hESDTruePi0DalitzSecConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
			}
		}
	}
   
 
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::ProcessVirtualGammasCandidates(){
	
	Double_t magField = fInputEvent->GetMagneticField();


	if( magField  < 0.0 ){
		magField =  1.0;
	}
	else {
		magField =  -1.0;
	}
	
	
	for(Int_t virtualGammaIndex=0; virtualGammaIndex < fGoodVirtualGammas->GetEntries();  virtualGammaIndex++ ){
		AliAODConversionPhoton *Vgamma=dynamic_cast<AliAODConversionPhoton*>(fGoodVirtualGammas->At(virtualGammaIndex));
		
		AliESDtrack* positronVgamma = fESDEvent->GetTrack( Vgamma->GetTrackLabelPositive() );
		AliESDtrack* electronVgamma = fESDEvent->GetTrack( Vgamma->GetTrackLabelNegative() );
		
		Bool_t isPhoton         = kFALSE;
		Bool_t isPi0Dalitz      = kFALSE;
		Bool_t isEtaDalitz      = kFALSE;
		Bool_t isJPsi           = kFALSE;
		Bool_t isTrueEposENeg   = kFALSE;
		Bool_t isMotherPrimary  = kFALSE;

		
		if ( fDoMesonQA > 0 ) {
		  
			Double_t psiPair = GetPsiPair(positronVgamma,electronVgamma);
			Double_t deltaPhi = magField * TVector2::Phi_mpi_pi( electronVgamma->GetConstrainedParam()->Phi()-positronVgamma->GetConstrainedParam()->Phi());
		
		
			hESDEposEnegPsiPairDPhi[fiCut]->Fill(deltaPhi,psiPair);  
			hESDEposEnegPsiPairEta[fiCut]->Fill(psiPair,Vgamma->Eta());
			hESDEposEnegDPhiEta[fiCut]->Fill(deltaPhi,Vgamma->Eta());
			
			hESDEposEnegInvMassPt[fiCut]->Fill(Vgamma->M(),Vgamma->Pt());
			
			if( fMCEvent ) {

                Int_t virtualGammaMCLabel = Vgamma->GetMCParticleLabel(fMCEvent);
				
				TParticle * negativeMC = 0x0; 
				TParticle * positiveMC = 0x0;   
				TParticle * mcVgamma   = 0x0;
				
				if( virtualGammaMCLabel != -1 ){
                    mcVgamma = (TParticle*)fMCEvent->Particle(virtualGammaMCLabel);
				}
				
                negativeMC = (TParticle*)Vgamma->GetNegativeMCDaughter(fMCEvent);
                positiveMC = (TParticle*)Vgamma->GetPositiveMCDaughter(fMCEvent);
				
						
				if( negativeMC && positiveMC ) {
					if( positiveMC->GetPdgCode() == ::kPositron && negativeMC->GetPdgCode() == ::kElectron ) 
						isTrueEposENeg = kTRUE;
                    if( positiveMC->GetMother(0) > -1 && positiveMC->GetMother(0) < fMCEvent->GetNumberOfPrimaries() )
						isMotherPrimary = kTRUE;
				}
				if(mcVgamma){
					// Check if it is a true photon
					if(mcVgamma->GetPdgCode() == 22){
						isPhoton = kTRUE;
					} else if(mcVgamma->GetPdgCode() == 443){
						isJPsi = kTRUE;
					} else if( IsDalitz( mcVgamma ) ){
						if     ( mcVgamma->GetPdgCode() == 111 ) isPi0Dalitz = kTRUE;
						else if( mcVgamma->GetPdgCode() == 221 ) isEtaDalitz = kTRUE;
					}
				}

				if(isPhoton){
					hESDEposEnegTruePhotonInvMassPt[fiCut]->Fill(Vgamma->M(),Vgamma->Pt());
					hESDEposEnegTruePhotonPsiPairDPhi[fiCut]->Fill(deltaPhi,psiPair);
					if( Vgamma->Pt() > 1.5 ){
						hESDEposEnegTruePhotonPsiPairDPhiPtCut[fiCut]->Fill(deltaPhi,psiPair);
					}
				} else if(isJPsi){
					hESDEposEnegTrueJPsiInvMassPt[fiCut]->Fill(Vgamma->M(),Vgamma->Pt());
				} else if(isPi0Dalitz){
				  
				        Double_t psiPairMC = GetPsiPairMC(positiveMC, negativeMC);
					hESDEposEnegTruePi0DalitzPsiPairMC[fiCut]->Fill(psiPairMC);
					hESDEposEnegTruePi0DalitzInvMassPt[fiCut]->Fill(Vgamma->M(),Vgamma->Pt());
					hESDEposEnegTruePi0DalitzPsiPairDPhi[fiCut]->Fill(deltaPhi,psiPair);
					hESDEposEnegTruePi0DalitzPsiPairEta[fiCut]->Fill(psiPair,Vgamma->Eta());
					hESDEposEnegTruePi0DalitzDPhiEta[fiCut]->Fill(deltaPhi,Vgamma->Eta());
					
				        if( isMotherPrimary ) hESDEposEnegTruePrimPi0DalitzInvMass[fiCut]->Fill( Vgamma->M() );
				} else if(isEtaDalitz){
					hESDEposEnegTrueEtaDalitzInvMassPt[fiCut]->Fill(Vgamma->M(),Vgamma->Pt());
					hESDEposEnegTrueEtaDalitzPsiPairDPhi[fiCut]->Fill(deltaPhi,psiPair);
					if( isMotherPrimary ) hESDEposEnegTruePrimEtaDalitzInvMass[fiCut]->Fill( Vgamma->M() );
				} else if ( isTrueEposENeg && mcVgamma ){
					
					hESDEposEnegTrueMotherInvMassPt[fiCut]->Fill(Vgamma->M(), Vgamma->Pt());    
				}	
				  
				if ( isTrueEposENeg ) hESDEposEnegTrueInvMassPt[fiCut]->Fill(Vgamma->M(), Vgamma->Pt());	    
				
			}
		}
	}
}  
  
//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::ProcessElectronCandidates(){

	Double_t magField = fInputEvent->GetMagneticField();

	if( magField  < 0.0 ){
		magField =  1.0;
	} else {
		magField =  -1.0;
	}

	vector<Int_t> lGoodElectronIndexPrev(0);
	vector<Int_t> lGoodPositronIndexPrev(0);
	
	for(UInt_t i = 0; i < fSelectorElectronIndex.size(); i++){
		
		AliESDtrack* electronCandidate = fESDEvent->GetTrack(fSelectorElectronIndex[i]);
		
		Bool_t IsMCFromMBHeader = kTRUE;
		
		if( fMCEvent ) {
		  
			Int_t isMCFromMBHeader = -1;
			Int_t labelelectron = TMath::Abs( electronCandidate->GetLabel() );
			
			if(((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetSignalRejection() != 0) {
                isMCFromMBHeader = ((AliConvEventCuts*)fCutEventArray->At(fiCut))->IsParticleFromBGEvent(labelelectron,fMCEvent,fInputEvent);
				if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetSignalRejection() != 3) continue;
			}    
			
			if( isMCFromMBHeader != 2 ) IsMCFromMBHeader = kFALSE;
		}
			
		if(! ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelected(electronCandidate) ) continue;
		
		lGoodElectronIndexPrev.push_back(   fSelectorElectronIndex[i] );
		hESDDalitzElectronPt[fiCut]->Fill(electronCandidate->Pt());
		hESDDalitzElectronPhi[fiCut]->Fill(electronCandidate->Phi());
		
		
		if( fMCEvent ) {
			
		      Int_t labelelectron = TMath::Abs( electronCandidate->GetLabel() );
            if( labelelectron < fMCEvent->GetNumberOfTracks() ){
                TParticle* electron = fMCEvent->Particle(labelelectron);
				if( electron->GetPdgCode() ==  11 ){
                    if( labelelectron < fMCEvent->GetNumberOfPrimaries() ){
						hESDTrueElectronPt[fiCut]->Fill(electronCandidate->Pt());    //primary electron
					} else {
						hESDTrueSecElectronPt[fiCut]->Fill(electronCandidate->Pt()); //secondary electron
					}
                    if( IsPi0DalitzDaughter(labelelectron) == kTRUE && labelelectron < fMCEvent->GetNumberOfPrimaries()  ) {
                        if( electron->GetMother(0) < fMCEvent->GetNumberOfPrimaries() ) {
							hESDTruePi0DalitzElectronPt[fiCut]->Fill(electronCandidate->Pt());
						  if( IsMCFromMBHeader == kTRUE ) hESDTruePi0DalitzElectronPtMB[fiCut]->Fill(electronCandidate->Pt());
						} else {
							hESDTruePi0DalitzSecElectronPt[fiCut]->Fill(electronCandidate->Pt());
						}
					}
				}
			}
		}
	}

	for(UInt_t i = 0; i < fSelectorPositronIndex.size(); i++){

		AliESDtrack* positronCandidate = fESDEvent->GetTrack( fSelectorPositronIndex[i] );
		
		Bool_t IsMCFromMBHeader = kTRUE;
		
		if( fMCEvent ) {
		
			Int_t isMCFromMBHeader = -1;
			Int_t labelpositron = TMath::Abs( positronCandidate->GetLabel() );
			
			
			if(((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetSignalRejection() != 0) {
                isMCFromMBHeader = ((AliConvEventCuts*)fCutEventArray->At(fiCut))->IsParticleFromBGEvent(labelpositron,fMCEvent,fInputEvent);
				if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetSignalRejection() != 3) continue;
			}  
			
			if( isMCFromMBHeader != 2 ) IsMCFromMBHeader = kFALSE;
		}
	
		if(! ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelected(positronCandidate) ) continue;
		
		lGoodPositronIndexPrev.push_back( fSelectorPositronIndex[i] );
		hESDDalitzPositronPt[fiCut]->Fill( positronCandidate->Pt() );
		hESDDalitzPositronPhi[fiCut]->Fill( positronCandidate->Phi() );
			
		if( fMCEvent ) {
			Int_t labelpositron = TMath::Abs( positronCandidate->GetLabel() );
			
            if( labelpositron < fMCEvent->GetNumberOfTracks() ) {
                TParticle* positron = fMCEvent->Particle(labelpositron);
				if( positron->GetPdgCode() ==  -11 ){
                    if( labelpositron < fMCEvent->GetNumberOfPrimaries() ){
						hESDTruePositronPt[fiCut]->Fill(positronCandidate->Pt());
					} else {
						hESDTrueSecPositronPt[fiCut]->Fill(positronCandidate->Pt());
					}
                    if( IsPi0DalitzDaughter(labelpositron) == kTRUE && labelpositron <  fMCEvent->GetNumberOfPrimaries() ) {
                        if( positron->GetMother(0) < fMCEvent->GetNumberOfPrimaries() ){
							hESDTruePi0DalitzPositronPt[fiCut]->Fill(positronCandidate->Pt());
							if ( IsMCFromMBHeader == kTRUE )hESDTruePi0DalitzPositronPtMB[fiCut]->Fill(positronCandidate->Pt());
						} else {
							hESDTruePi0DalitzSecPositronPt[fiCut]->Fill(positronCandidate->Pt());
						}
					}
				}
			}
		}
	}

	vector<Bool_t> lElectronPsiIndex(lGoodElectronIndexPrev.size(), kTRUE);
	vector<Bool_t> lPositronPsiIndex(lGoodPositronIndexPrev.size(), kTRUE);

	if( ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->DoPsiPairCut() == kTRUE ){
		for( UInt_t i = 0; i < lGoodElectronIndexPrev.size(); i++ ) {
			AliESDtrack *electronCandidate = fESDEvent->GetTrack(lGoodElectronIndexPrev[i]);
			for(UInt_t j = 0; j <  lGoodPositronIndexPrev.size(); j++){
				AliESDtrack *positronCandidate = fESDEvent->GetTrack(lGoodPositronIndexPrev[j]);
				Double_t psiPair = GetPsiPair(positronCandidate,electronCandidate);
				Double_t deltaPhi = magField * TVector2::Phi_mpi_pi( electronCandidate->GetConstrainedParam()->Phi()-positronCandidate->GetConstrainedParam()->Phi());

				if( ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->IsFromGammaConversion(psiPair,deltaPhi ) ){
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
		AliESDtrack *electronCandidate = fESDEvent->GetTrack(lGoodElectronIndex[i]);
		AliKFParticle electronCandidateKF( *electronCandidate->GetConstrainedParam(), ::kElectron );
		TLorentzVector electronCandidateTLV;
		electronCandidateTLV.SetXYZM(electronCandidate->GetConstrainedParam()->Px(),electronCandidate->GetConstrainedParam()->Py(),electronCandidate->GetConstrainedParam()->Pz(),TDatabasePDG::Instance()->GetParticle(  ::kElectron   )->Mass());
		

		for(UInt_t j = 0; j < lGoodPositronIndex.size(); j++){

			AliESDtrack *positronCandidate = fESDEvent->GetTrack(lGoodPositronIndex[j]);
			AliKFParticle positronCandidateKF( *positronCandidate->GetConstrainedParam(), ::kPositron );
			TLorentzVector positronCandidateTLV;
		
			positronCandidateTLV.SetXYZM(positronCandidate->GetConstrainedParam()->Px(),positronCandidate->GetConstrainedParam()->Py(),positronCandidate->GetConstrainedParam()->Pz(),TDatabasePDG::Instance()->GetParticle(  ::kPositron   )->Mass());
			TLorentzVector *virtualPhotonTLV = 0;
			AliKFConversionPhoton* virtualPhoton = NULL;
			AliAODConversionPhoton *vphoton;
			
			if( ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->GetUseElectronMCSmearing() && ((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->UseMCPSmearing() && fMCEvent){
				TLorentzVector smearelectronCandidateTLV = ((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->SmearElectron(electronCandidateTLV);
				TLorentzVector smearpositronCandidateTLV = ((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->SmearElectron(positronCandidateTLV);
				virtualPhotonTLV = new TLorentzVector( smearelectronCandidateTLV + smearpositronCandidateTLV );
				//TLorentzVector temp = electronCandidateTLV + positronCandidateTLV;
				//cout<<"Mass:  "<<"NoSmearElectrons  "<<temp.M() <<"   SmearElectrons  "<<virtualPhotonTLV->M()<<endl;
				vphoton= new AliAODConversionPhoton(virtualPhotonTLV); 
				vphoton->SetMass(virtualPhotonTLV->M());
			} else {
				virtualPhoton = new AliKFConversionPhoton(electronCandidateKF,positronCandidateKF);
				if( fSetProductionVertextoVGamma == kTRUE ){
					AliKFVertex primaryVertexImproved(*fInputEvent->GetPrimaryVertex());
					primaryVertexImproved+=*virtualPhoton;
					virtualPhoton->SetProductionVertex(primaryVertexImproved);			
				}

				vphoton= new AliAODConversionPhoton(virtualPhoton);
			}
			//virtualPhoton->SetTrackLabels( lGoodPositronIndex[j], lGoodElectronIndex[i]);
			vphoton->SetTrackLabels( lGoodPositronIndex[j], lGoodElectronIndex[i]);
			
			if( fMCEvent ) {
				Int_t labeln=TMath::Abs(electronCandidate->GetLabel());
				Int_t labelp=TMath::Abs(positronCandidate->GetLabel());
                TParticle *fNegativeMCParticle = fMCEvent->Particle(labeln);
                TParticle *fPositiveMCParticle = fMCEvent->Particle(labelp);
				if( fPositiveMCParticle && fNegativeMCParticle) {		
					//virtualPhoton->SetMCLabelPositive(labelp);
					//virtualPhoton->SetMCLabelNegative(labeln);
					vphoton->SetMCLabelPositive(labelp);
					vphoton->SetMCLabelNegative(labeln);
						
				}		
			}
			
					
			fGoodVirtualGammas->Add(  vphoton );
			
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
	
	//Computing mixing event
	if(  fDoMesonQA > 0 ) {
		for(UInt_t i = 0; i < lGoodElectronIndex.size(); i++){
			AliESDtrack *electronCandidate1 = fESDEvent->GetTrack(lGoodElectronIndex[i]);
			AliKFParticle electronCandidate1KF( *electronCandidate1->GetConstrainedParam(), ::kElectron );
			
			for(UInt_t j = i+1; j < lGoodElectronIndex.size(); j++){
			
				AliESDtrack *electronCandidate2 = fESDEvent->GetTrack(lGoodElectronIndex[j]);
				AliKFParticle electronCandidate2KF( *electronCandidate2->GetConstrainedParam(), ::kElectron );
				AliKFConversionPhoton* virtualPhoton = new AliKFConversionPhoton(electronCandidate1KF,electronCandidate2KF);
				
				//AliKFVertex primaryVertexImproved(*fInputEvent->GetPrimaryVertex());
				//primaryVertexImproved+=*virtualPhoton;
				//virtualPhoton->SetProductionVertex(primaryVertexImproved);

				AliAODConversionPhoton *vphoton = new AliAODConversionPhoton(virtualPhoton); 
				hESDEposEnegLikeSignBackInvMassPt[fiCut]->Fill(vphoton->M(),vphoton->Pt());
				delete vphoton;
				delete virtualPhoton;
				vphoton = 0x0;
				virtualPhoton = 0x0;            
			}
		}   
	
	
		for(UInt_t i = 0; i < lGoodPositronIndex.size(); i++){
			AliESDtrack *positronCandidate1 = fESDEvent->GetTrack(lGoodPositronIndex[i]);
			AliKFParticle positronCandidate1KF( *positronCandidate1->GetConstrainedParam(), ::kPositron );
			
			for(UInt_t j = i+1; j < lGoodPositronIndex.size(); j++){
			
				AliESDtrack *positronCandidate2 = fESDEvent->GetTrack(lGoodPositronIndex[j]);
				AliKFParticle positronCandidate2KF( *positronCandidate2->GetConstrainedParam(), ::kPositron );
				AliKFConversionPhoton* virtualPhoton = new AliKFConversionPhoton(positronCandidate1KF,positronCandidate2KF);
				//AliKFVertex primaryVertexImproved(*fInputEvent->GetPrimaryVertex());
				//primaryVertexImproved+=*virtualPhoton;
				//virtualPhoton->SetProductionVertex(primaryVertexImproved);
				
				AliAODConversionPhoton *vphoton = new AliAODConversionPhoton(virtualPhoton); 
				hESDEposEnegLikeSignBackInvMassPt[fiCut]->Fill(vphoton->M(),vphoton->Pt());	
				delete vphoton;
				delete virtualPhoton;  
				vphoton = 0x0;
				virtualPhoton = 0x0;
		
			}
		}   
	}
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::CalculatePi0DalitzCandidates(){

	// Conversion Gammas
	if( fGoodGammas->GetEntries() > 0 && fGoodVirtualGammas->GetEntries() > 0 ){
		vector<Bool_t> lGoodVirtualGamma(fGoodVirtualGammas->GetEntries(), kFALSE);
		
		for(Int_t GammaIndex=0; GammaIndex<fGoodGammas->GetEntries(); GammaIndex++){
			AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fGoodGammas->At(GammaIndex));
			if (gamma==NULL) continue;
			for(Int_t virtualGammaIndex=0;virtualGammaIndex<fGoodVirtualGammas->GetEntries();virtualGammaIndex++){
				AliAODConversionPhoton *Vgamma=dynamic_cast<AliAODConversionPhoton*>(fGoodVirtualGammas->At(virtualGammaIndex));
				if (Vgamma==NULL) continue;
				//Check for same Electron ID
			     if(gamma->GetTrackLabelPositive()    == Vgamma->GetTrackLabelPositive() ||
				gamma->GetTrackLabelNegative()    == Vgamma->GetTrackLabelNegative() ||
				gamma->GetTrackLabelNegative()    == Vgamma->GetTrackLabelPositive() ||
				gamma->GetTrackLabelPositive()    == Vgamma->GetTrackLabelNegative() ) continue;

				AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma,Vgamma);
				pi0cand->SetLabels(GammaIndex,virtualGammaIndex);
						

				if( ( ((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetEtaShift())) ){
			
					//cout<< "Meson Accepted "<<endl;
						Int_t zbin = 0;
						Int_t mbin = 0;
					     if(fDoTHnSparse && ((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->DoBGCalculation()){
					        if(((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->BackgroundHandlerType() == 0){
						
						      zbin= fBGHandler[fiCut]->GetZBinIndex(fESDEvent->GetPrimaryVertex()->GetZ());
						      if( ((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->UseTrackMultiplicity() ){
							mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
						      } else {
							mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGoodGammas->GetEntries());
						      }
						}
					     }
					if ( fDoMesonQA > 0 ) {
					  
						      
						      
						        
						      if( fMCEvent ) {
							
                                Int_t virtualGammaMCLabel = Vgamma->GetMCParticleLabel(fMCEvent);
							    
							    if( virtualGammaMCLabel > -1 ) {
							    
                                TParticle * negativeMC = (TParticle*)Vgamma->GetNegativeMCDaughter(fMCEvent);
                                TParticle * positiveMC = (TParticle*)Vgamma->GetPositiveMCDaughter(fMCEvent);
                               // TParticle * virtualGammaMotherMC = (TParticle*)fMCEvent->Particle(virtualGammaMCLabel);
							    
							    if( negativeMC->GetPdgCode() == 11 && positiveMC->GetPdgCode() == -11) {  // Electrons ...

								  
								if( IsPi0DalitzDaughter( Vgamma->GetMCLabelPositive() ) ){
								  
								  hESDTruePi0EposEnegInvMassPi0Pt[fiCut]->Fill( Vgamma->M(), pi0cand->Pt() );
								
								
								}
							    							    
							      }

							    }
						      }
										
					  }
			
					if(  ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoMassCut() == kTRUE ) {
						if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->MassCut( pi0cand->Pt() , Vgamma->M() ) == kTRUE ){

							hESDMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
							
							if( fDoTHnSparse ){
							Double_t sparesFill[4] = {pi0cand->M(),pi0cand->Pt(),(Double_t)zbin,(Double_t)mbin};
							sESDMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
							}
														
							
							if(fMCEvent){
							  ProcessTrueMesonCandidates(pi0cand,gamma,Vgamma);
							}
						
							if ( fDoMesonQA > 0 ) {
							  
								
							  
											
								
								AliESDtrack* positronVgamma = fESDEvent->GetTrack( Vgamma->GetTrackLabelPositive() );
								Double_t momPositron[3];
								positronVgamma->GetConstrainedPxPyPz(momPositron);
							
								AliESDtrack* electronVgamma = fESDEvent->GetTrack( Vgamma->GetTrackLabelNegative() );
								Double_t momElectron[3];
								electronVgamma->GetConstrainedPxPyPz(momElectron);
							
								TVector3 vGamma(gamma->GetPx(),gamma->GetPy(),gamma->GetPz());;
								TVector3 vPositron(momPositron[0],momPositron[1],momPositron[2]);
								TVector3 vElectron(momElectron[0],momElectron[1],momElectron[2]);
	
							
								hESDMotherInvMassOpeningAngleGammaElectron[fiCut]->Fill( vGamma.Angle(vPositron) );
								hESDMotherInvMassOpeningAngleGammaElectron[fiCut]->Fill( vGamma.Angle(vElectron) );
								
								
								
								if ( pi0cand->M() > 0.05 && pi0cand->M() < 0.17){
								  
								    if( fDoMesonQA > 1 ){
								    TLorentzVector electronCandidateTLV(vElectron,TDatabasePDG::Instance()->GetParticle(  ::kElectron   )->Mass());
								    TLorentzVector positronCandidateTLV(vPositron,TDatabasePDG::Instance()->GetParticle(  ::kPositron   )->Mass());
								    TLorentzVector gammaTLV(vGamma,0);
								    TLorentzVector GammaElectronCandidateTLV = gammaTLV + electronCandidateTLV;
								    TLorentzVector GammaPositronCandidateTLV = gammaTLV + positronCandidateTLV;
								    Double_t sparesDalitzPlot[2] = {GammaElectronCandidateTLV.M()*GammaElectronCandidateTLV.M(),GammaPositronCandidateTLV.M()*GammaPositronCandidateTLV.M()};										  
								    sESDMotherDalitzPlot[fiCut]->Fill(sparesDalitzPlot);
								    
								    
								    }
								
								  
								    hESDMotherPi0PtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetEtaShift());  
                                    hESDMotherPi0PtAlpha[fiCut]->Fill(pi0cand->Pt(),TMath::Abs(pi0cand->GetAlpha()));
								    hESDMotherPi0PtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle()); 
								    
								} 
								
								
								if( lGoodVirtualGamma[virtualGammaIndex] == kFALSE ) {
						
								       
								        if( pi0cand->M() > 0.1 && pi0cand->M() < 0.145 ){
									FillElectronQAHistos(Vgamma);
									
									lGoodVirtualGamma[virtualGammaIndex] = kTRUE;
									fNVirtualGammas++;
									}
								}
							}
						}
					} else {
						hESDMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
						if( fDoTHnSparse ){
						Double_t sparesFill[4] = {pi0cand->M(),pi0cand->Pt(),(Double_t)zbin,(Double_t)mbin};
						sESDMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
						}
						
						if(fMCEvent){
							  ProcessTrueMesonCandidates(pi0cand,gamma,Vgamma);
						}
						
						if ( fDoMesonQA > 0 ) {
							
							
							if ( pi0cand->M() > 0.05 && pi0cand->M() < 0.17){
								  
								    hESDMotherPi0PtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetEtaShift());  
                                    hESDMotherPi0PtAlpha[fiCut]->Fill(pi0cand->Pt(),TMath::Abs(pi0cand->GetAlpha()));
								    hESDMotherPi0PtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle()); 
								    
								    
							} 
							
							if( lGoodVirtualGamma[virtualGammaIndex] == kFALSE ) {
							  
								if( pi0cand->M() > 0.1 && pi0cand->M() < 0.145 ){
									
								  FillElectronQAHistos(Vgamma);
								  lGoodVirtualGamma[virtualGammaIndex] = kTRUE;
								  fNVirtualGammas++;
								
								}
							}
						}
					}
				
					if( fDoChicAnalysis) {
					
						hESDPi0MotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());						
						Double_t diffMass = pi0cand->M() - Vgamma->M();
						hESDPi0MotherDiffInvMassPt[fiCut]->Fill( diffMass , pi0cand->Pt() );
				
						if( Vgamma->M() > 2.5 && Vgamma->M() < 3.4){
							hESDPi0MotherDiffLimInvMassPt[fiCut]->Fill( diffMass , pi0cand->Pt() );
						}
						
						if(fMCEvent){
						  ProcessTrueChicCandidates(pi0cand,gamma,Vgamma);
						}
						
						
					}
					
					
				}
				delete pi0cand;
				pi0cand=0x0;
			}
		}
	}
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::CalculateBackground(){

	Int_t zbin= fBGHandler[fiCut]->GetZBinIndex(fESDEvent->GetPrimaryVertex()->GetZ());
	Int_t mbin = 0;

		
	if(((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->UseTrackMultiplicity()){
        mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
	} else {
        mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGoodGammas->GetEntries());
	}
	

		AliGammaConversionAODBGHandler::GammaConversionVertex *bgEventVertex = NULL;

		  if(((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->UseTrackMultiplicity()){
			for(Int_t nEventsInBG=0;nEventsInBG<fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){

				AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
				if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0 ){
					bgEventVertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
				}

				for(Int_t iCurrent=0;iCurrent<fGoodVirtualGammas->GetEntries();iCurrent++){
					AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGoodVirtualGammas->At(iCurrent));

					for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
						AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));

						if(fMoveParticleAccordingToVertex == kTRUE ){
							MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
						}
						
						if(((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
							RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
						}

						AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
						backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
						
						if( ( ((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE, ((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetEtaShift()))){
							if(  ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoMassCut() == kTRUE ) {
							
								if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->MassCut( backgroundCandidate->Pt() , currentEventGoodV0.M() ) == kTRUE ){
							
									hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
									if(fDoTHnSparse){
									Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
									sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
									}
								}
							} else {
								hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
								if(fDoTHnSparse){
								  Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
								  sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
								}
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
					if(fMoveParticleAccordingToVertex == kTRUE || ((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0 ){
						bgEventVertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
					}
					for(Int_t iCurrent=0;iCurrent<fGoodVirtualGammas->GetEntries();iCurrent++){
						AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGoodVirtualGammas->At(iCurrent));
						for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
							AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
							if(fMoveParticleAccordingToVertex == kTRUE ){
								MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
							}
							
							if(((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
							    RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
							}

							AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
							backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
						
							if((((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE,((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetEtaShift()))){
								if(  ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoMassCut() == kTRUE ) {
									if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->MassCut( backgroundCandidate->Pt() , currentEventGoodV0.M() ) == kTRUE ){
										hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
										if(fDoTHnSparse){
										Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
										sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
										}
									}
								} else {
									hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
									if(fDoTHnSparse){
									Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
									sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
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
void AliAnalysisTaskGammaConvDalitzV1::UpdateEventByEventData(){
	//see header file for documentation

		if(fGoodGammas->GetEntries() > 0 ){
			if( ((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->UseTrackMultiplicity() ){
				fBGHandler[fiCut]->AddEvent(fGoodGammas,fESDEvent->GetPrimaryVertex()->GetX(),fESDEvent->GetPrimaryVertex()->GetY(),fESDEvent->GetPrimaryVertex()->GetZ(),fNumberOfESDTracks,fEventPlaneAngle);
			} else { // means we use #V0s for multiplicity
				fBGHandler[fiCut]->AddEvent(fGoodGammas,fESDEvent->GetPrimaryVertex()->GetX(),fESDEvent->GetPrimaryVertex()->GetY(),fESDEvent->GetPrimaryVertex()->GetZ(),fGoodGammas->GetEntries(),fEventPlaneAngle);
			}
		}
	
}

//______________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate, AliAODConversionPhoton *TrueVirtualGammaCandidate)
{
	// Process True Mesons
	
	if( TrueGammaCandidate->GetV0Index() < fESDEvent->GetNumberOfV0s() ){
		
	  

		Bool_t isTruePi0 = kFALSE;
		Bool_t isTrueEta = kFALSE;
		
            Int_t gammaMCLabel = TrueGammaCandidate->GetMCParticleLabel(fMCEvent);
		Int_t gammaMotherLabel = -1;
		
		

		
		
		//Checking if the gamma candidate is a real gamma
		
		
		if(gammaMCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
			// Daughters Gamma 0
            TParticle * negativeMC = (TParticle*)TrueGammaCandidate->GetNegativeMCDaughter(fMCEvent);
            TParticle * positiveMC = (TParticle*)TrueGammaCandidate->GetPositiveMCDaughter(fMCEvent);
            TParticle * gammaMC = (TParticle*)fMCEvent->Particle(gammaMCLabel);

			if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
				if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
					if(gammaMC->GetPdgCode() == 22){ // ... with Gamma Mother
						gammaMotherLabel=gammaMC->GetFirstMother();
					}
				}
			}
		}


		
		//Checking if the virtual gamma is a real virtual gamma
		
        Int_t virtualGammaMCLabel = TrueVirtualGammaCandidate->GetMCParticleLabel(fMCEvent);
		Int_t virtualGammaMotherLabel = -1;
		Int_t virtualGamma = -1;
		


		if( virtualGammaMCLabel != -1 ){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
			// Daughters Gamma 1
            TParticle * negativeMC = (TParticle*)TrueVirtualGammaCandidate->GetNegativeMCDaughter(fMCEvent);
            TParticle * positiveMC = (TParticle*)TrueVirtualGammaCandidate->GetPositiveMCDaughter(fMCEvent);
            TParticle * virtualGammaMotherMC = (TParticle*)fMCEvent->Particle(virtualGammaMCLabel);

			if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...

				if( virtualGammaMotherMC->GetPdgCode() != 22 ){
					virtualGammaMotherLabel=virtualGammaMCLabel;
					virtualGamma = 1;
					
				} else if (negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
					virtualGammaMotherLabel=virtualGammaMotherMC->GetFirstMother();
					virtualGamma = 0; //no virtual gamma
				}
			}
		}

		
		//Checking if both gamma and virtual gamma comming from Pi0 or Eta
		

		if( gammaMotherLabel >= 0 && ( gammaMotherLabel == virtualGammaMotherLabel) ){
            if(((TParticle*)fMCEvent->Particle(virtualGammaMotherLabel))->GetPdgCode() == 111){
				isTruePi0=kTRUE;
				if (CheckVectorForDoubleCount(fVectorDoubleCountTruePi0s,gammaMotherLabel)) fHistoDoubleCountTruePi0InvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
            if(((TParticle*)fMCEvent->Particle(virtualGammaMotherLabel))->GetPdgCode() == 221){
				isTrueEta=kTRUE;
				if (CheckVectorForDoubleCount(fVectorDoubleCountTrueEtas,gammaMotherLabel)) fHistoDoubleCountTrueEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
		}

		

		if(  isTruePi0 || isTrueEta ){ // True Pion or Eta
            
            
			if ( virtualGamma == 1 ) { //True Dalitz	
				Float_t weighted= 1;
				
                Float_t weightMatBudget = 1.;       
                if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->GetMaterialBudgetWeightsInitialized ()) {
                    weightMatBudget = ((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(TrueGammaCandidate);
                }
              
				
				if ( isTruePi0 && fDoMesonQA > 0 ) {
						
						if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
						  
						        
						   
							hESDTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetEtaShift()); 
                            hESDTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),TMath::Abs(Pi0Candidate->GetAlpha()));
							hESDTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle()); 
							
							if( fDoMesonQA > 1 ) {
								AliESDtrack* positronVgamma = fESDEvent->GetTrack( TrueVirtualGammaCandidate->GetTrackLabelPositive() );
								Double_t momPositron[3];
								positronVgamma->GetConstrainedPxPyPz(momPositron);
							
								AliESDtrack* electronVgamma = fESDEvent->GetTrack( TrueVirtualGammaCandidate->GetTrackLabelNegative() );
								Double_t momElectron[3];
								electronVgamma->GetConstrainedPxPyPz(momElectron);
							
								TVector3 vGamma(TrueGammaCandidate->GetPx(),TrueGammaCandidate->GetPy(),TrueGammaCandidate->GetPz());;
								TVector3 vPositron(momPositron[0],momPositron[1],momPositron[2]);
								TVector3 vElectron(momElectron[0],momElectron[1],momElectron[2]);
	                                                        
								TLorentzVector electronCandidateTLV(vElectron,TDatabasePDG::Instance()->GetParticle(  ::kElectron   )->Mass());
								TLorentzVector positronCandidateTLV(vPositron,TDatabasePDG::Instance()->GetParticle(  ::kPositron   )->Mass());
								TLorentzVector gammaTLV(vGamma,0);
								TLorentzVector GammaElectronCandidateTLV = gammaTLV + electronCandidateTLV;
								TLorentzVector GammaPositronCandidateTLV = gammaTLV + positronCandidateTLV;
								Double_t sparesDalitzPlot[2] = {GammaElectronCandidateTLV.M()*GammaElectronCandidateTLV.M(),GammaPositronCandidateTLV.M()*GammaPositronCandidateTLV.M()};										  
								sESDTruePi0DalitzPlot[fiCut]->Fill(sparesDalitzPlot);
							}
							
						}
				}
				
				
				
				if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoWeights() ) { 
                    if(((AliConvEventCuts*)fCutEventArray->At(fiCut))->IsParticleFromBGEvent(gammaMotherLabel, fMCEvent,fInputEvent)){
                        if (((TParticle*)fMCEvent->Particle(gammaMotherLabel))->Pt()>0.005){
                            weighted= ((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetWeightForMeson(gammaMotherLabel,fMCEvent,fInputEvent);
						}
					}
				}

				hESDTrueMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*weightMatBudget);
				hESDTrueMotherW0WeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				hESDTrueMotherDalitzInvMassPt[fiCut]->Fill( TrueVirtualGammaCandidate->M(),Pi0Candidate->Pt(),weighted*weightMatBudget);

                if(gammaMotherLabel < fMCEvent->GetNumberOfPrimaries()){ // Only primary pi0 for efficiency calculation
					hESDTruePrimaryMotherInvMassPt[fiCut]->Fill( Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*weightMatBudget);
					hESDTruePrimaryMotherW0WeightingInvMassPt[fiCut]->Fill( Pi0Candidate->M(), Pi0Candidate->Pt() );
				
                    hESDTruePrimaryMotherInvMassMCPt[fiCut]->Fill(Pi0Candidate->M(),((TParticle*)fMCEvent->Particle(virtualGammaMotherLabel))->Pt(),weighted*weightMatBudget);
					if(isTruePi0){ // Only primaries for unfolding
                        hESDTruePrimaryPi0DalitzESDPtMCPt[fiCut]->Fill(Pi0Candidate->Pt(),((TParticle*)fMCEvent->Particle(virtualGammaMotherLabel))->Pt(),weighted*weightMatBudget);
					}
				} else { // Secondary Meson
					Float_t weightedSec= 1;
					if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoWeights() ) { 
                        Int_t secMotherLabel = ((TParticle*)fMCEvent->Particle(gammaMotherLabel))->GetMother(0);
                        if(((AliConvEventCuts*)fCutEventArray->At(fiCut))->IsParticleFromBGEvent(gammaMotherLabel, fMCEvent, fInputEvent) && fMCEvent->Particle(gammaMotherLabel)->GetPdgCode()==310){
                            weightedSec= ((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetWeightForMeson(secMotherLabel, fMCEvent, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
						}
					}
					hESDTrueSecondaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*weightMatBudget); 
				}
			} else if ( virtualGamma == 0 ){
                
                Float_t weightMatBudget = 1.;       
                if (fDoMaterialBudgetWeightingOfGammasForTrueMesons && ((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->GetMaterialBudgetWeightsInitialized ()) {
                    weightMatBudget = ((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(TrueGammaCandidate) * ((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->GetMaterialBudgetCorrectingWeightForTrueGamma(TrueVirtualGammaCandidate);
                }
                
				Float_t weighted= 1;
				if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoWeights() ) {
                    if(((AliConvEventCuts*)fCutEventArray->At(fiCut))->IsParticleFromBGEvent(gammaMotherLabel, fMCEvent,fInputEvent)){
                        if (((TParticle*)fMCEvent->Particle(gammaMotherLabel))->Pt()>0.005){
                            weighted= ((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetWeightForMeson(gammaMotherLabel,fMCEvent,fInputEvent);
						}
					}
				}
				hESDTrueMotherPi0GGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*weightMatBudget); // Pi0 from GG
				hESDTrueMotherPi0GGW0WeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
				
                if( gammaMotherLabel < fMCEvent->GetNumberOfPrimaries() ){
					hESDTruePrimaryMotherPi0GGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted*weightMatBudget);
				} else {
					Float_t weightedSec= 1;
					if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoWeights() ) { 
                        Int_t secMotherLabel = ((TParticle*)fMCEvent->Particle(gammaMotherLabel))->GetMother(0);
                        if(((AliConvEventCuts*)fCutEventArray->At(fiCut))->IsParticleFromBGEvent(gammaMotherLabel, fMCEvent, fInputEvent) && fMCEvent->Particle(gammaMotherLabel)->GetPdgCode()==310){
                            weightedSec= ((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetWeightForMeson(secMotherLabel, fMCEvent, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
						}
					}
					hESDTrueSecondaryMotherPi0GGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec*weightMatBudget);
				}
			}
		}

		if( !isTruePi0 && !isTrueEta ){ // Background
			if( gammaMotherLabel > -1 && virtualGammaMotherLabel > -1 && virtualGamma == 0 ){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
				hESDTrueBckGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			} else { // No photon or without mother
				hESDTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
			}
		}
	}
}

//______________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::ProcessTrueChicCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate, AliAODConversionPhoton *TruejpsiCandidate)
{


	
	if(	TrueGammaCandidate->GetV0Index() < fESDEvent->GetNumberOfV0s()	){
		
	  
	  
		//Checking gamma
		
        Int_t gammaMCLabel = TrueGammaCandidate->GetMCParticleLabel(fMCEvent);
		Int_t gammaMotherLabel = -1;

		
		
		if( gammaMCLabel != -1){// Gamma is Combinatorial; MC Particles don't belong to the same Mother
			// Daughters Gamma 0
			
            TParticle * negativeMC = (TParticle*)TrueGammaCandidate->GetNegativeMCDaughter(fMCEvent);
            TParticle * positiveMC = (TParticle*)TrueGammaCandidate->GetPositiveMCDaughter(fMCEvent);
            TParticle * gammaMC = (TParticle*)fMCEvent->Particle(gammaMCLabel);

			if( TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11 ){  // Electrons ...
				if( negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() == 5 ){ // ... From Conversion ...
					if(gammaMC->GetPdgCode() == 22){ // ... with Gamma Mother
						gammaMotherLabel=gammaMC->GetFirstMother();
					}
				}
			}
		}


		
		//Checking jpsi
		
        Int_t jpsiMCLabel 		= TruejpsiCandidate->GetMCParticleLabel(fMCEvent);
		Int_t jpsiMotherLabel 		= -1;

		

		if( jpsiMCLabel != -1 ){ 
			
            TParticle * negativeMC = (TParticle*)TruejpsiCandidate->GetNegativeMCDaughter(fMCEvent);
            TParticle * positiveMC = (TParticle*)TruejpsiCandidate->GetPositiveMCDaughter(fMCEvent);
            TParticle * jpsiMC = (TParticle*)fMCEvent->Particle(jpsiMCLabel);

			if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...

									
				if(jpsiMC->GetPdgCode() == 443){
				  
						jpsiMotherLabel=jpsiMC->GetFirstMother();
						
				}
				 
								
			}
		}


		if( gammaMotherLabel>=0 && ( gammaMotherLabel == jpsiMotherLabel) ){
          if(			((TParticle*)fMCEvent->Particle(jpsiMotherLabel))->GetPdgCode() == 445   ||
                    ((TParticle*)fMCEvent->Particle(jpsiMotherLabel))->GetPdgCode() == 10443 ||
                    ((TParticle*)fMCEvent->Particle(jpsiMotherLabel))->GetPdgCode() == 20443 ){
				
					hESDTrueMotherChiCInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
					hESDTrueMotherChiCDiffInvMassPt[fiCut]->Fill(Pi0Candidate->M()-TruejpsiCandidate->M(),Pi0Candidate->Pt());
		  }
		}  
		
	}
}



//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::RotateParticleAccordingToEP(AliAODConversionPhoton *gamma, Double_t previousEventEP, Double_t thisEventEP){

	previousEventEP=previousEventEP+TMath::Pi();
	thisEventEP=thisEventEP+TMath::Pi();
	Double_t rotationValue= thisEventEP-previousEventEP;
	gamma->RotateZ(rotationValue);
}


//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex){
   //see header file for documentation

	Double_t dx = vertex->fX - fESDEvent->GetPrimaryVertex()->GetX();
	Double_t dy = vertex->fY - fESDEvent->GetPrimaryVertex()->GetY();
	Double_t dz = vertex->fZ - fESDEvent->GetPrimaryVertex()->GetZ();

	Double_t movedPlace[3] = {particle->GetConversionX() - dx,particle->GetConversionY() - dy,particle->GetConversionZ() - dz};
	particle->SetConversionPoint(movedPlace);
}


//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::CountESDTracks(){

	// Using standard function for setting Cuts
	Bool_t selectPrimaries=kTRUE;
	AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
	EsdTrackCuts->SetMaxDCAToVertexZ(2);
	EsdTrackCuts->SetEtaRange(-0.8, 0.8);
	EsdTrackCuts->SetPtRange(0.15);
	EsdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kBoth);

	fNumberOfESDTrackskBoth = 0;
	for(Int_t iTracks = 0; iTracks < fESDEvent->GetNumberOfTracks(); iTracks++){
		AliESDtrack* curTrack = fESDEvent->GetTrack(iTracks);
		if(!curTrack) continue;
		if(EsdTrackCuts->AcceptTrack(curTrack) ) fNumberOfESDTrackskBoth++;
	}
	delete EsdTrackCuts;
	EsdTrackCuts=0x0;

	return;
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::ProcessMCParticles()
{

// Loop over all primary MC particle

    for(Int_t i = 0; i < fMCEvent->GetNumberOfPrimaries(); i++) {
	
        TParticle* particle = (TParticle *)fMCEvent->Particle(i);
		if (!particle) continue;
	
		Bool_t mcIsFromMB = kTRUE;
		Int_t isMCFromMBHeader = -1;
			
		if(((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetSignalRejection() != 0) {
            isMCFromMBHeader = ((AliConvEventCuts*)fCutEventArray->At(fiCut))->IsParticleFromBGEvent(i,fMCEvent,fInputEvent);
			if(isMCFromMBHeader == 0 && ((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetSignalRejection() != 3) continue;
			if(isMCFromMBHeader != 2) mcIsFromMB = kFALSE;
		}    
		
        if(((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kFALSE)){
			hMCAllGammaPt[fiCut]->Fill(particle->Pt()); // All MC Gamma
			
			if( IsPi0DalitzDaughter( i ) == kTRUE ){
			   
			    hMCAllGammaPi0Pt[fiCut]->Fill(particle->Pt());
			   
			}
			
			
			
		}
		
        if(((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCEvent,kTRUE)){
			
			hMCConvGammaPt[fiCut]->Fill(particle->Pt());
			
			if(fDoMesonQA > 0 ) {
			  hMCConvGammaEta[fiCut]->Fill( particle->Eta());
			  if(mcIsFromMB) {
                hMCConvGammaR[fiCut]->Fill( ((TParticle*)fMCEvent->Particle(particle->GetFirstDaughter()))->R() );
                hMCConvGammaPtR[fiCut]->Fill( particle->Pt(), ((TParticle*)fMCEvent->Particle(particle->GetFirstDaughter()))->R() );
			  }
			}
		
			if(mcIsFromMB){
			  
			    hMCConvGammaRSPt[fiCut]->Fill(particle->Pt());
			    
			}
			
			 if( IsPi0DalitzDaughter( i ) == kTRUE ){
			   
			    hMCConvGammaPi0Pt[fiCut]->Fill(particle->Pt());
			   
			 }
			
		} // Converted MC Gamma
	
        if(((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelectedMC(i,fMCEvent)){
		  	if( particle->GetPdgCode() == -11) {
				hMCAllPositronsPt[fiCut]->Fill(particle->Pt()); // All positrons
				if(fDoMesonQA > 0 ) hMCAllPositronsEta[fiCut]->Fill( particle->Eta() );
				
                if ( particle->GetMother(0) > -1 && particle->GetMother(0) < fMCEvent->GetNumberOfPrimaries() ) {
				    
				      if( IsPi0DalitzDaughter( i ) == kTRUE ){
					hMCDecayPositronPi0Pt[fiCut]->Fill(particle->Pt()); //Positrons from Pi0->Dalitz
				      }
				  
				}
			}
			if( particle->GetPdgCode() ==  11){
				hMCAllElectronsPt[fiCut]->Fill(particle->Pt()); // All electrons
				if (fDoMesonQA > 0 )  hMCAllElectronsEta[fiCut]->Fill(particle->Eta()); // All electrons

                if ( particle->GetMother(0) > -1 && particle->GetMother(0) < fMCEvent->GetNumberOfPrimaries() ) {
				    
				      if( IsPi0DalitzDaughter( i ) == kTRUE ){
					hMCDecayElectronPi0Pt[fiCut]->Fill(particle->Pt()); //Electrons from Pi0->Dalitz
				      }
				  
				}

			}		
		}
	
        if(((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelectedMC( particle,fMCEvent,((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetEtaShift() ) ){
			Float_t weighted= 1;
			if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoWeights() ) { 
                if(((AliConvEventCuts*)fCutEventArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent, fInputEvent)){
					if (particle->Pt()>0.005){
                        weighted= ((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetWeightForMeson(i, fMCEvent,fInputEvent);
					}
				}
			}
			if(particle->GetPdgCode() == 111)hMCPi0GGPt[fiCut]->Fill( particle->Pt() , weighted); // All MC Pi0 GG decay
			if(particle->GetPdgCode() == 221)hMCEtaGGPt[fiCut]->Fill( particle->Pt() , weighted); // All MC Eta GG decay
		}
	
		
		Int_t labelgamma 	  = -1;
		Int_t labelelectron = -1;
		Int_t labelpositron = -1;


        if( ((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelectedMCDalitz(particle,fMCEvent,labelelectron,labelpositron,labelgamma,((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetEtaShift())){
			Float_t weighted= 1;
			if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoWeights() ) { 
                if(((AliConvEventCuts*)fCutEventArray->At(fiCut))->IsParticleFromBGEvent(i, fMCEvent,fInputEvent)){
					if (particle->Pt()>0.005){
                        weighted= ((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetWeightForMeson(i, fMCEvent,fInputEvent);
					}
				}
			}
			if(particle->GetPdgCode() == 111){
			  hMCPi0Pt[fiCut]->Fill(particle->Pt(), weighted); // All MC Pi0
			  hMCPi0WOWeightPt[fiCut]->Fill(particle->Pt());
						
			}
			if(particle->GetPdgCode() == 221)hMCEtaPt[fiCut]->Fill(particle->Pt(), weighted); // All MC Eta
			
			// Check the acceptance for gamma and electrons
			
            TParticle *gamma    = fMCEvent->Particle(labelgamma);
            TParticle *electron = fMCEvent->Particle(labelelectron);
            TParticle *positron = fMCEvent->Particle(labelpositron);
			
	
            if(((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->PhotonIsSelectedMC(gamma,fMCEvent,kFALSE) &&
            ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelectedMC(labelelectron,fMCEvent) &&
            ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelectedMC(labelpositron,fMCEvent) ) {
				
				Double_t massDalitz 	= -1;
				Double_t angleGammaEpos = -1;
				Double_t angleGammaEneg = -1;
				Double_t massGammaEpos  = -1;
				Double_t massGammaEneg  = -1;
				
				if( fDoMesonQA > 0 ){
				  
					TLorentzVector TLVEpos,TLVEneg,TLVgamma,TLVDalitz,TLVGammaEpos,TLVGammaEneg;
					Double_t electronMass = TDatabasePDG::Instance()->GetParticle(  ::kElectron   )->Mass();					
					TLVEpos.SetXYZM(positron->Px(),positron->Py(),positron->Pz(),electronMass);
					TLVEneg.SetXYZM(electron->Px(),electron->Py(),electron->Pz(),electronMass);
					TLVgamma.SetXYZM(gamma->Px(),gamma->Py(),gamma->Pz(),0);
					
					TVector3 V3gamma(gamma->Px(),gamma->Py(),gamma->Pz());
					angleGammaEpos = V3gamma.Angle(TLVEpos.Vect());
					angleGammaEneg = V3gamma.Angle(TLVEneg.Vect());
					
					TLVDalitz =  TLVEpos + TLVEneg ;
					massDalitz = TLVDalitz.M();
					TLVGammaEpos = TLVEpos + TLVgamma;
					TLVGammaEneg = TLVEneg + TLVgamma;
					
					massGammaEpos = TLVGammaEpos.M();
					massGammaEneg = TLVGammaEneg.M();
					
					
					
				}
			
				if(particle->GetPdgCode() == 111){ 
				  
					hMCPi0InAccPt[fiCut]->Fill(particle->Pt() , weighted); // MC Pi0Dalitz with gamma and e+e- in acc
					hMCPi0WOWeightInAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0 with gamma in acc NOT weighted
					hMCPi0DalitzGammaPt[fiCut]->Fill( gamma->Pt() );
					hMCPi0DalitzPositronPt[fiCut]->Fill( positron->Pt()  );
					hMCPi0DalitzElectronPt[fiCut]->Fill( electron->Pt()  );
					
					
					
					if ( fDoMesonQA > 0 ){
					  
					      hMCPi0EposEnegInvMassPt[fiCut]->Fill(massDalitz,particle->Pt());
					      hMCPi0InAccOpeningAngleGammaElectron[fiCut]->Fill(angleGammaEpos);
					      hMCPi0InAccOpeningAngleGammaElectron[fiCut]->Fill(angleGammaEneg);
					      
					      if ( fDoMesonQA > 1 ) {
						
						Double_t sMCDalitzPlot[2] = {massGammaEneg*massGammaEneg,massGammaEpos*massGammaEpos};
						
						sMCPi0DalitzPlot[fiCut]->Fill(sMCDalitzPlot);
						
					      }
							      
					}
				}
				if(particle->GetPdgCode() == 221){ 
					hMCEtaInAccPt[fiCut]->Fill(particle->Pt(), weighted ); // MC EtaDalitz with gamma and e+e- in acc
					if( fDoMesonQA  > 0 ) hMCEtaEposEnegInvMassPt[fiCut]->Fill( massDalitz, particle->Pt());
				}
			}	
		}
		
		Int_t labelgammaChiC=-1;
		Int_t labelpositronChiC=-1;
		Int_t labelelectronChiC=-1;

        if(((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelectedMCChiC(particle,fMCEvent,labelelectronChiC,labelpositronChiC,labelgammaChiC,((AliConvEventCuts*)fCutEventArray->At(fiCut))->GetEtaShift())){
		
			hMCChiCPt[fiCut]->Fill(particle->Pt()); // All MC ChiC
            TParticle * gammaChiC  =fMCEvent->Particle(labelgammaChiC);

            if( ((AliConversionPhotonCuts*)fCutGammaArray->At(fiCut))->PhotonIsSelectedMC( gammaChiC,fMCEvent,kFALSE) &&
                ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelectedMC(labelelectronChiC,fMCEvent) &&
                ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelectedMC(labelpositronChiC,fMCEvent) ){
				hMCChiCInAccPt[fiCut]->Fill(particle->Pt()); // All MC ChiC
			} 
		}
	}
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaConvDalitzV1::IsDalitz(TParticle *fMCMother) const
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

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskGammaConvDalitzV1::IsPi0DalitzDaughter( Int_t label ) const
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

void AliAnalysisTaskGammaConvDalitzV1::FillElectronQAHistos(AliAODConversionPhoton *Vgamma) const
{

	AliESDtrack *positronVgamma = 0;
	AliESDtrack *electronVgamma = 0;
	
	Double_t clsToFPos = -1.0;
	Double_t clsToFNeg = -1.0;
	
	Double_t NumClsITSPos = -1.0;
	Double_t NumClsITSNeg = -1.0;
	Double_t NumClsTPCPos = -1.0;
	Double_t NumClsTPCNeg = -1.0;
	Double_t nCrossedRowsTPCPos = -1.0;
	Double_t nCrossedRowsTPCNeg = -1.0;
	
	Float_t dcaToVertexXYPos = -1.0;
	Float_t dcaToVertexZPos  = -1.0;
	Float_t dcaToVertexXYNeg = -1.0;
	Float_t dcaToVertexZNeg  = -1.0;
	
	Double_t nSigmaPosTPC = -999.;
	Double_t nSigmaNegTPC = -999.;
	
	positronVgamma = fESDEvent->GetTrack( Vgamma->GetTrackLabelPositive() );
	electronVgamma = fESDEvent->GetTrack( Vgamma->GetTrackLabelNegative() );
	clsToFPos = ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->GetNFindableClustersTPC(positronVgamma);
	clsToFNeg = ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->GetNFindableClustersTPC(electronVgamma);
	
	nSigmaPosTPC = ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->GetPIDResponse()->NumberOfSigmasTPC(positronVgamma, AliPID::kElectron) ;
	nSigmaNegTPC = ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->GetPIDResponse()->NumberOfSigmasTPC(electronVgamma, AliPID::kElectron) ;	
	
	NumClsITSPos =  positronVgamma->GetNcls(0); //Get number of ITS clusters
	NumClsITSNeg =  electronVgamma->GetNcls(0);
	NumClsTPCPos =  positronVgamma->GetNcls(1);  //Get number of TPC clusters
	NumClsTPCNeg =  electronVgamma->GetNcls(1);
	
	nCrossedRowsTPCPos = positronVgamma->GetTPCCrossedRows();
	nCrossedRowsTPCNeg = electronVgamma->GetTPCCrossedRows();

	Float_t bPos[2];
	Float_t bCovPos[3];
	positronVgamma->GetImpactParameters(bPos,bCovPos);
	
	if (bCovPos[0]<=0 || bCovPos[2]<=0) {
		AliDebug(1, "Estimated b resolution lower or equal zero!");
		bCovPos[0]=0; bCovPos[2]=0;
	}
	
	Float_t bNeg[2];
	Float_t bCovNeg[3];
	positronVgamma->GetImpactParameters(bNeg,bCovNeg);
	
	if (bCovNeg[0]<=0 || bCovNeg[2]<=0) {
		AliDebug(1, "Estimated b resolution lower or equal zero!");
		bCovNeg[0]=0; bCovNeg[2]=0;
	}
	
	dcaToVertexXYPos = bPos[0];
	dcaToVertexZPos  = bPos[1];
	dcaToVertexXYNeg = bNeg[0];
	dcaToVertexZNeg  = bNeg[1];
		
	if( electronVgamma->P() > 0.3 &&  electronVgamma->P() < 0.45 ){ 
	  
		hESDDalitzElectronAfterEtaPCut[fiCut]->Fill( electronVgamma->Eta() );
		hESDDalitzElectronAfterNClsITSPCut[fiCut]->Fill( NumClsITSNeg );
		hESDDalitzElectronAfterNFindClsTPCPCut[fiCut]->Fill(clsToFNeg);
		hESDDalitzElectronAfterNClsTPCPCut[fiCut]->Fill( NumClsTPCNeg );
		hESDDalitzElectronAfterNCrossedRowsTPCPCut[fiCut]->Fill( nCrossedRowsTPCNeg );
	}
	if( positronVgamma->P() > 0.3 &&  positronVgamma->P() < 0.45 ){ 
	  
		hESDDalitzPositronAfterEtaPCut[fiCut]->Fill( positronVgamma->Eta() );
		hESDDalitzPositronAfterNClsITSPCut[fiCut]->Fill( NumClsITSPos );
		hESDDalitzPositronAfterNFindClsTPCPCut[fiCut]->Fill(clsToFPos);
		hESDDalitzPositronAfterNClsTPCPCut[fiCut]->Fill( NumClsTPCPos );
		hESDDalitzPositronAfterNCrossedRowsTPCPCut[fiCut]->Fill( nCrossedRowsTPCPos );
	}
	
	hESDDalitzElectronAfterPt[fiCut]->Fill(  electronVgamma->Pt()  );
	hESDDalitzPositronAfterPt[fiCut]->Fill(  positronVgamma->Pt()  );
			
	hESDDalitzElectronAfterEta[fiCut]->Fill( electronVgamma->Eta() );
	hESDDalitzPositronAfterEta[fiCut]->Fill( positronVgamma->Eta() );
	
	hESDDalitzElectronAfterPhi[fiCut]->Fill( electronVgamma->Phi() );
	hESDDalitzPositronAfterPhi[fiCut]->Fill( positronVgamma->Phi() );
		
	hESDDalitzElectronAfterNFindClsTPC[fiCut]->Fill(clsToFNeg,electronVgamma->Pt());
	hESDDalitzPositronAfterNFindClsTPC[fiCut]->Fill(clsToFPos,positronVgamma->Pt());
	
	hESDDalitzElectronAfterNClsTPC[fiCut]->Fill(  NumClsTPCNeg,electronVgamma->Pt());
	hESDDalitzPositronAfterNClsTPC[fiCut]->Fill(  NumClsTPCPos,positronVgamma->Pt());
	
	hESDDalitzElectronAfterNCrossedRowsTPC[fiCut]->Fill( nCrossedRowsTPCNeg, electronVgamma->Pt() );
	hESDDalitzPositronAfterNCrossedRowsTPC[fiCut]->Fill( nCrossedRowsTPCPos, positronVgamma->Pt() );
	
	hESDDalitzElectronAfterNClsITS[fiCut]->Fill( NumClsITSNeg);
	hESDDalitzPositronAfterNClsITS[fiCut]->Fill( NumClsITSPos);
	
	hESDDalitzPosEleAfterDCAxy[fiCut]->Fill(  dcaToVertexXYNeg, electronVgamma->Pt() );
	hESDDalitzPosEleAfterDCAz[fiCut]->Fill(   dcaToVertexZNeg,  electronVgamma->Pt() );
	hESDDalitzPosEleAfterDCAxy[fiCut]->Fill(  dcaToVertexXYPos, positronVgamma->Pt() );
	hESDDalitzPosEleAfterDCAz[fiCut]->Fill(   dcaToVertexZPos,  positronVgamma->Pt() );
		
	hESDDalitzElectronAfterTPCdEdxVsP[fiCut]->Fill( electronVgamma->P(),nSigmaNegTPC);
	hESDDalitzPositronAfterTPCdEdxVsP[fiCut]->Fill( positronVgamma->P(), nSigmaPosTPC);
	
	hESDDalitzElectronAfterTPCdEdxVsEta[fiCut]->Fill( electronVgamma->Eta(),nSigmaNegTPC);
	hESDDalitzPositronAfterTPCdEdxVsEta[fiCut]->Fill( positronVgamma->Eta(),nSigmaPosTPC);

	hESDDalitzElectronAfterTPCdEdxVsPhi[fiCut]->Fill( electronVgamma->Phi(),nSigmaNegTPC);
	hESDDalitzPositronAfterTPCdEdxVsPhi[fiCut]->Fill( positronVgamma->Phi(),nSigmaPosTPC);
				
	hESDDalitzElectronAfterTPCdEdxSignalVsP[fiCut]->Fill( electronVgamma->P(), TMath::Abs(electronVgamma->GetTPCsignal()));
	hESDDalitzPositronAfterTPCdEdxSignalVsP[fiCut]->Fill( positronVgamma->P(), TMath::Abs(positronVgamma->GetTPCsignal()));
	
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskGammaConvDalitzV1::GetPsiPairMC(const TParticle *fMCPosParticle, const TParticle *fMCNegParticle) const
{



        TVector3 posDaughter;
	TVector3 negDaughter;

	posDaughter.SetXYZ( fMCPosParticle->Px(), fMCPosParticle->Py(), fMCPosParticle->Pz() );
	negDaughter.SetXYZ( fMCNegParticle->Px(), fMCNegParticle->Py(), fMCNegParticle->Pz() );

	Double_t deltaTheta = negDaughter.Theta() - posDaughter.Theta();
	Double_t openingAngle =  posDaughter.Angle( negDaughter );  //TMath::ACos( posDaughter.Dot(negDaughter)/(negDaughter.Mag()*posDaughter.Mag()) );

	if( openingAngle < 1e-20 ) return 0.;

	Double_t psiAngle = TMath::ASin( deltaTheta/openingAngle );

	return psiAngle;


}



//_____________________________________________________________________________
Double_t AliAnalysisTaskGammaConvDalitzV1::GetPsiPair( const AliESDtrack *trackPos, const AliESDtrack *trackNeg ) const
{
	//
	// This angle is a measure for the contribution of the opening in polar
	// direction ?0 to the opening angle ? Pair
	//
	// Ref. Measurement of photons via conversion pairs with the PHENIX experiment at RHIC
	//      Mas   ter Thesis. Thorsten Dahms. 2005
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

//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::SetLogBinningXTH2(TH2* histoRebin){
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
Bool_t AliAnalysisTaskGammaConvDalitzV1::CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked)
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
