/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  Task for Hadron/LeadingParticle-HFE & HFE-HFE Correlation         //
//  Author: Florian Herrmann  &  Denise Godoy                         //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include <set>
#include <vector>
#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TMath.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"
#include "THashList.h"
#include "TRandom2.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDHandler.h"
#include  "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliVEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"
#include "AliVertexingHFUtils.h"

#include "AliAnalysisTaskHaHFECorrel.h"
#include "TGeoGlobalMagField.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"
#include "AliAODInputHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "AliPhysicsSelection.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"
#include "AliESDv0KineCuts.h"
#include "AliESDv0.h"
#include "AliEMCALTrack.h"
#include "AliMagF.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"

#include "AliEventplane.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"

#include "AliEventPoolManager.h"


ClassImp(AliAnalysisTaskHaHFECorrel)
ClassImp(AliBasicParticleHaHFE)

using std::cout;
using std::endl;

//________________________________________________________________________
AliAnalysisTaskHaHFECorrel::AliAnalysisTaskHaHFECorrel(const char *name)
: AliAnalysisTaskSE(name)
,fRunNumber(0)
,fUseTender(kFALSE)
,fWhichPeriod(2016)
,fUseEPOS(kFALSE)
,fUseKFforPhotonicPartner(kFALSE)
,fMaxPtEvent(999)
,fMinPtEvent(0)
,fMaxNTr(999)
,fMinNTr(0)
,fVarZVTXCut(0)
,fMaxElectronEta(0.8)
,fMinElectronEta(-0.8)
,fMaxHadronEta(0.8)
,fMinHadronEta(-0.8)
,fVarEleOpt(0)
,fElectronkAny(kFALSE)
,fElectronkFirst(kFALSE)
,fTPCnCut(100)
,fTPCndEdxCut(80)
,fITSnCut(3)
,fITSSharedClusterCut(999)
,fEleDCAr(1)
,fEleDCAz(2)
,fUseTRD(0)
,fUseITSsa(1)
,fSigmaITScut(20)
,fSigmaTOFcut(2.)
,fSigmaTPCcutLow(-1.)
,fSigmaTPCcutHigh(3)
,fVarPhotOpt(0)
,fPhotElecPtCut(0.0)
,fPhotElecSigmaTPCcut(3)
,fPhotElecTPCnCut(80)
,fPhotElecITSrefitCut(kTRUE)
,fPhotCorrCase(2)
,fAssNonEleTPCcut(-4)
,fVarHadOpt(0)
,fHTPCnCut(100)
,fHITSrefitCut(kTRUE)
,fHTPCrefitCut(kTRUE)
,fHadDCAr(0.25)
,fHadDCAz(1.)
,fHadkAny(kTRUE)
,fHadTOFmatch(kFALSE)
,fOpeningAngleCut(1000.)
,fInvmassCut(0.14)
,fChi2Cut(3.5)
,fDCAcut(999)
,fTRDQA(kFALSE)
,fMCTrueCorrelation(kTRUE)
,fUseEventWeights(kTRUE)
,fCorrHadron(kTRUE)
,fCorrLParticle(kTRUE)
,fMixedEvent(kTRUE)
,fPionEtaProduction(kFALSE)
,fRecEff(kTRUE)
,fTagEff(kTRUE)
,fHadCont(kTRUE)
,fOneTimeCheck(kFALSE)
,fLParticle(kTRUE)
,fESD(0)
,fesdTrackCuts(0)
,fAOD(0)
,fVevent(0)
,fpidResponse(0)
,fMultSelection(0)
,fCentrality(0)
,fPoolMgr(0)
,fPoolIsFilled(0)
,fMC(0)
,fStack(0)
,fMCparticle(0)
,fMCarray(0)
,fMCheader(0)
,PdgTable(0)
,PDGMap()
,fTracks_tender(0)
,fCaloClusters_tender(0)
,fEventCuts()
,fCuts(0)
,fIsMC(kFALSE)
,fIsAOD(kTRUE)
,fCFM(0)
,fPID(0)
,fPIDqa(0)
,fOutputList(0)
,fOutputListMain(0)
,fOutputListLP(0)
,fOutputListHadron(0)
,fOutputListQA(0)
,fNoEvents(0)
,fNoEventsNTr(0)
,fMCNoEvents(0)
,fHFENoEvents(0)
,fDiffractiveType(0)
,fV0ACTrueInel(0)
,fV0TrueMinInel(0)
,fV0TrueMinInelNTr(0)
,fV0ACTriggered(0)
,fV0MinTriggered(0)
,fV0MinTriggeredNTr(0)
,fTriggerWeight()
,fVtxEtaNTr(0)
,fVtxBeforeNTrAcc(0)
,fVtxAfterNTrAcc()
,fVtxRecBeforeNTr(0)
,fVtxRecAfterNTr(0)
,fVtxWeight()
,fTrkpt(0)
,fEtaVtxZ(0)
,fSPDVtxRes(0)              
,fDiffSPDPrimVtx(0)         
,fSPDnTrAcc(0)    
,fSPDnTrCorrMax(0)
,fSPDnTrGen(0)              
,fDiffSPDMCVtx(0)           
,fnTrAccMaxGen(0)
,fnTrAccGen(0)
,fnTrAccGenTrueInel(0)
,fnTrAccMinGen(0)           
,fnTrAccMeanGen(0)          
,fnTrAccMax(0)              
,fnTrAccMin(0)              
,fnTrAccMean(0)
,fMCThrustTagged(0)
,fMCSpherTagged(0)
,fRecLPTagged(0)
,fMultCorrTagged(0)
,fNHadTagged(0)
,fNHadTaggedA(0)
,fNHadTaggedB(0)
,fNHadTaggedC(0)
,fMeanPtTagged(0)
,fMeanPtTaggedA(0)
,fMeanPtTaggedB(0)
,fMeanPtTaggedC(0)
,fMCThrustNTagged(0)
,fMCSpherNTagged(0)
,fRecLPNTagged(0)
,fMultCorrNTagged(0)
,fNHadNTagged(0)
,fNHadNTaggedA(0)
,fNHadNTaggedB(0)
,fNHadNTaggedC(0)
,fMeanPtNTagged(0)
,fMeanPtNTaggedA(0)
,fMeanPtNTaggedB(0)
,fMeanPtNTaggedC(0)
,fPt2Tagged(0)
,fPt2NTagged(0)
,fMothMCThrustTagged(0)
,fMothMCSpherTagged(0)
,fMothRecLPTagged(0)
,fMothMultCorrTagged(0)
,fMothNHadTagged(0)
,fMothMeanPtTagged(0)
,fMothMCThrustNTagged(0)
,fMothMCSpherNTagged(0)
,fMothRecLPNTagged(0)
,fMothMultCorrNTagged(0)
,fMothNHadNTagged(0)
,fMothMeanPtNTagged(0)
,fMCThrustTaggedH(0)
,fMCSpherTaggedH(0)
,fRecLPTaggedH(0)
,fMultCorrTaggedH(0)
,fNHadTaggedH(0)
,fMeanPtTaggedH(0)
,fMCThrustNTaggedH(0)
,fMCSpherNTaggedH(0)
,fRecLPNTaggedH(0)
,fMultCorrNTaggedH(0)
,fNHadNTaggedH(0)
,fMeanPtNTaggedH(0)
,fMothMCThrustTaggedH(0)
,fMothMCSpherTaggedH(0)
,fMothRecLPTaggedH(0)
,fMothMultCorrTaggedH(0)
,fMothNHadTaggedH(0)
,fMothMeanPtTaggedH(0)
,fMothMCThrustNTaggedH(0)
,fMothMCSpherNTaggedH(0)
,fMothRecLPNTaggedH(0)
,fMothMultCorrNTaggedH(0)
,fMothNHadNTaggedH(0)
,fMothMeanPtNTaggedH(0)
,fMultiplicity(0)
,fSPDMultiplicity(0)
,fRunList(0)
,fElectronTrackCuts(0)
,fElectronTrackTPCChi2(0)
,fElectronTrackTPCCrossedRows(0)
,fElectronTrackTPCNcls(0)
,fElectronTrackTPCNclsdEdx(0)
,fElectronTrackTPCFrac(0)
,fElectronTrackITSNcls(0)
,fElectronTrackITSChi2(0)
,fElectronTrackITSLayer(0)
,fElectronTrackRefit(0)
,fElectronTrackDCA(0)
,fElectronTrackITSCuts(0)
,fPhotTrackITSCuts(0)
,fHadronTrackCuts(0)
,fHadronTrackTPCNcls(0)
,fHadronTrackRefit(0)
,fHadronTrackDCA(0)
,fHadronTrackDCA_woITSAny(0)
,fHadronTrackDCA_wITSAny(0)
,fHistITSnSig(0)
,fHistTOFnSig(0)
,fHistTPCnSig(0)
,fHistTPCnSigITScut(0)
,fHistTPCnSigTOFcut(0)
,fHistTPCnSigITSTOFcut(0)
,fHistITSnSigTOFTPCcut(0)
,fCheckNHadronScaling(0)
,fCheckNPhotHadScaling(0)
,fCheckTaggedEvent(0)
,fHadContPvsPt(0)
,fHadContEtaPhiPt(0)
,fHadContTPCEtaPhiPt(0)
,fHadContPPhiEtaTPC(0)
,fHadContamination(0)
,fHadContaminationPt(0)
,fHadContMC(0)
,fHadContMCPt(0)
,fInclElecPtEta(0)
,fInclElecPtEtaWRecEff(0)
,fInclElecP(0)
,fULSElecPt(0)
,fULSElecPtWRecEff(0)
,fLSElecPt(0)
,fLSElecPtWRecEff(0)
,fInvmassLS(0)
,fInvmassULS(0)
,fInvmassMCTrue(0)
,fRecMCInvMass(0)
,fPhotMixULS(0)
,fPhotMixLS(0)
,fPhotPt1PtMTag(0)
,fPhotPt1PtMNTag(0)
,fPhotPt1Pt2(0)
,fPhotPt1Pt2Only(0)
,fPhotPt1Pt2Corr(0)
,fPhotPt1Pt2MC(0)
,fPhotPt1RecPt2(0)
,fPhotPt1RecPt2Corr(0)
,fPhotPt1RecPt2Rec(0)
,fPhotPt1RecPt2RecCorr(0)
,fPhotPt1Pt2Rec(0)
,fPhotPt1Pt2RecCorr(0)
,fPhotPt2MCRec(0)
,fPhotPt1E(0)
,fPhotPt1Pt2E(0)
,fPhotPt1ECorr(0)
,fPhotPt1Pt2ECorr(0)
,fPhotPt1Mass(0)
,fPhotPt1Mom(0)
,fOpeningAngleLS(0)
,fOpeningAngleULS(0)
,fCheckLSULS(0)
,fTagEtaPt1Pt2(0)
,fTagEtaPhiPt(0)
,fTagEtaZvtxPt(0)
,fTagEtaPhiPtwW(0)
,fTagEtaZvtxPtwW(0)
,fNonTagEtaPt1Pt2(0)
,fNonTagEtaPhiPt(0)
,fNonTagEtaZvtxPt(0)
,fNonTagEtaPhiPtwW(0)
,fNonTagEtaZvtxPtwW(0)
,fTagMotherPt(0)
,fTagEffInclMult(0)
,fTagEffULSMult(0)
,fTagEffInclBGMult(0)
,fTagEffULSBGMult(0)
,fTagTruePairsMult(0)
,fTagEffInclMultWoW(0)
,fTagEffULSMultWoW(0)
,fTagEffInclBGMultWoW(0)
,fTagEffULSBGMultWoW(0)
,fTagTruePairsMultWoW(0)
,fTagEffInclMultWoWS(0)
,fTagEffULSMultWoWS(0)
,fTagEffInclBGMultWoWS(0)
,fTagEffULSBGMultWoWS(0)
,fTagTruePairsMultWoWS(0)
,fTagEffIncl(0)
,fTagEffLS(0)
,fTagEffULS(0)
,fTagTruePairs(0)
,fTagEffInclWoWeight(0)
,fTagEffLSWoWeight(0)
,fTagEffULSWoWeight(0)
,fTagTruePairsWoWeight(0)
,fCorrectPiontoData()
,fCorrectEtatoData()
,fBgWeight()
,fHadRecEff()
,fEleRecEff()
,fSPDConfig(0)
,fSPDConfigHist()
,fSPDConfigProfiles()
,fSPDnTrAvg(0)
,fNonTagCorr()
,fAssPtHad_Nbins(0)
,fAssPtHad_Xmin()
,fAssPtHad_Xmax()
,fAssPtElec_Nbins(0)
,fAssPtElec_Xmin()
,fAssPtElec_Xmax()
,fElecTrigger(0)
,fInclElecPhi(0)
,fInclElecEta(0)
,fInclElecPhiEta(0)
,fULSElecPhi(0)
,fLSElecPhi(0)
,fElecDphi(0)
,fULSElecDphi(0)
,fLSElecDphi(0)
,fULSElecDphiDiffMethod(0)
,fLSElecDphiDiffMethod(0)
,fNoPartnerNoT(0)
,fNoPartnerNoTPt2(0)
,fTPartnerNoT(0)
,fTPartnerNoTPt2(0)
,fElecHadTrigger(0)
,fElecHadTriggerLS(0)
,fElecHadTriggerULS(0)
,fElecHadTriggerLSNoP(0)
,fElecHadTriggerULSNoP(0)
,fElecHadTriggerLSNoPCorr(0)
,fElecHadTriggerULSNoPCorr(0)
,fElecHadTriggerULSNoPCorrTrue(0)
,fHadContTrigger(0)
,fHadElecTrigger(0)
,fNonElecHadTrigger(0)
,fHadNonElecTrigger(0)
,fInclElecHa(0)
,fLSElecHa(0)
,fULSElecHa(0)
,fULSElecHaTrue(0)
,fSignalElecHa(0)
,fBackgroundElecHa(0)
,fMCElecHaHadron(0)
,fElecHaHa(0)
,fElecHaLSNoPartner(0)
,fElecHaULSNoPartner(0)
,fElecHaLSNoPartnerCorrTrue(0)
,fElecHaULSNoPartnerCorrTrue(0)
,fElecHaLSNoPartnerCorr(0)
,fElecHaULSNoPartnerCorr(0)
,fMCElecHaTruePartner(0)
,fMCElecHaNoPartner(0)
,fMCElecHaRemovedPartner(0)
,fMCElecHaTruePartnerTrigger(0)
,fMCElecHaTruePartnerTriggerWW(0)
,fMCElecHaNoPartnerTrigger(0)
,fMCElecHaNoPartnerTriggerWW(0)
,fMCElecHaRemovedPartnerTrigger(0)
,fElecHaMixedEvent(0)
,fLSElecHaMixedEvent(0)
,fULSElecHaMixedEvent(0)
,fTagHaMixedEvent(0)
,fNonTagHaMixedEvent(0)
,fElecLPTrigger(0)
,fElecLPTriggerLS(0)
,fElecLPTriggerULS(0)
,fElecLPTriggerLSNoP(0)
,fElecLPTriggerULSNoP(0)
,fElecLPTriggerULSNoPCorr(0)
,fHadContLPTrigger(0)
,fLPElecTrigger(0)
,fLPNonElecTrigger(0)
,fNonElecLPTrigger(0)
,fInclElecLP(0)
,fLSElecLP(0)
,fULSElecLP(0)
,fMCElecLPHadron(0)
,fElecLPHa(0)
,fElecLPLSNoPartner(0)
,fElecLPULSNoPartner(0)
,fElecLPLSNoPartnerCorrTrue(0)
,fElecLPULSNoPartnerCorrTrue(0)
,fElecLPLSNoPartnerCorr(0)
,fElecLPULSNoPartnerCorr(0)
,fMCElecLPTruePartner(0)
,fMCElecLPNoPartner(0)
,fMCElecLPRemovedPartner(0)
,fMCElecLPTruePartnerTrigger(0)
,fMCElecLPNoPartnerTrigger(0)
,fMCElecLPRemovedPartnerTrigger(0)
,fElecLPMixedEvent(0)
,fLSElecLPMixedEvent(0)
,fULSElecLPMixedEvent(0)
,fTagLPMixedEvent(0)
,fNonTagLPMixedEvent(0)
,fCheckMCVertex(0)
,fCheckMCPtvsRecPtHad(0)
,fCheckMCEtavsRecEtaHad(0)
,fCheckMCPhivsRecPhiHad(0)
,fMCHadPtEtaPhiVtx(0)
,fRecHadMCSecondaryCont(0)
,fRecHadMCPtEtaPhiVtx(0)
,fRecHadPtEtaPhiVtx(0)
,fRecHadPtEtaPhiVtxWRecEff(0)
,fCheckMCPtvsRecPtEle(0)
,fRecHFE(0)
,fMCElecPtEtaPhiVtx(0)
,fRecElecMCSecondaryCont(0)
,fRecElecPtEtaPhiVtx(0)
,fRecElecPtEtaPhiVtxWRecEff(0)
,fRecElecMCPtEtaPhiVtx(0)
,fMCElecPDG(0)
,fMCElecPtEtaPhiStrictVtx(0)
,fMCPi0Prod(0)
,fMCEtaProd(0)
,fMCPiPlusProd(0)
,fMCPiPlusProdV2(0)
,fMCBGProd(0)
,fMCLeadingParticle(0)
,fCompareLPRecCheck(0)
,fMCTruePoolMgr(0)
,fTrueMCHadronEventCuts(0)   
,fTrueMCHadronEventCutsZvtx(0)   
,fTrueMCHadronEventCutsZvtxMEv(0)   
,fTrueMCHadron(0)   
,fTrueMCElecHaTriggerEventCuts(0) 
,fTrueMCElecHaTrigger(0) 
,fTrueMCLPEventCuts(0)   
,fTrueMCLPEventCutsZvtx(0)   
,fTrueMCLPEventCutsZvtxMEv(0)   
,fTrueMCLP(0)   
,fTrueMCElecLPTriggerEventCuts(0) 
,fTrueMCElecLPTrigger(0) 
,fTrueElectronEta(0) 
,fRecHFEEtaWRecEff(0) 
,fTrueLPinAcceptanceEta(0) 
,fTrueLPEta(0) 
,fRecLPEta(0) 
,fTrueHadronEta(0) 
,fRecHadronEtaWRecEff(0) 
,fCompareLP(0) 
,fV0cutsESD(0)
,fV0cuts(0)
,fV0electrons(0)
,fV0pions(0)
,fV0protons(0)
,fhArmenteros(0)
,fEventsPerRun(0)
,fTRDnTrackRun(0)
,fV0tags(0)
,fTRDEtaPhi(0)
,fTRDNTracklets(0)
,fTRDV0NTracklets(0)
,fTRDSpectra(0)
,fTRDV0Spectra(0)
,fTRDMCSpectra(0)
{
    //Named constructor
   
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2,TList::Class());
    DefineOutput(3,TList::Class());
    DefineOutput(4,TList::Class());
    DefineOutput(5,TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskHaHFECorrel::AliAnalysisTaskHaHFECorrel()
: AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisElecHadCorrel")
,fRunNumber(0)
,fUseTender(kFALSE)
,fWhichPeriod(2016)
,fUseEPOS(kFALSE)
,fUseKFforPhotonicPartner(kFALSE)
,fMaxPtEvent(999)
,fMinPtEvent(0)
,fMaxNTr(999)
,fMinNTr(0)
,fVarZVTXCut(0)
,fMaxElectronEta(0.8)
,fMinElectronEta(-0.8)
,fMaxHadronEta(0.8)
,fMinHadronEta(-0.8)
,fVarEleOpt(0)
,fElectronkAny(kFALSE)
,fElectronkFirst(kFALSE)
,fTPCnCut(100)
,fTPCndEdxCut(80)
,fITSnCut(3)
,fITSSharedClusterCut(999)
,fEleDCAr(1)
,fEleDCAz(2)
,fUseTRD(0)
,fUseITSsa(1)
,fSigmaITScut(20)
,fSigmaTOFcut(2.)
,fSigmaTPCcutLow(-1.)
,fSigmaTPCcutHigh(3)
,fVarPhotOpt(0)
,fPhotElecPtCut(0.0)
,fPhotElecSigmaTPCcut(3)
,fPhotElecTPCnCut(80)
,fPhotElecITSrefitCut(kTRUE)
,fPhotCorrCase(2)
,fAssNonEleTPCcut(-4)
,fVarHadOpt(0)
,fHTPCnCut(100)
,fHITSrefitCut(kTRUE)
,fHTPCrefitCut(kTRUE)
,fHadDCAr(0.25)
,fHadDCAz(1.)
,fHadkAny(kTRUE)
,fHadTOFmatch(kFALSE)
,fOpeningAngleCut(1000.)
,fInvmassCut(0.14)
,fChi2Cut(3.5)
,fDCAcut(999)
,fTRDQA(kFALSE)
,fMCTrueCorrelation(kTRUE)
,fUseEventWeights(kTRUE)
,fCorrHadron(kTRUE)
,fCorrLParticle(kTRUE)
,fMixedEvent(kTRUE)
,fPionEtaProduction(kFALSE)
,fRecEff(kTRUE)
,fTagEff(kTRUE)
,fHadCont(kTRUE)
,fOneTimeCheck(kFALSE)
,fLParticle(kTRUE)
,fESD(0)
,fesdTrackCuts(0)
,fAOD(0)
,fVevent(0)
,fpidResponse(0)
,fMultSelection(0)
,fCentrality(0)
,fPoolMgr(0)
,fPoolIsFilled(0)
,fMC(0)
,fStack(0)
,fMCparticle(0)
,fMCarray(0)
,fMCheader(0)
,PdgTable(0)
,PDGMap()
,fTracks_tender(0)
,fCaloClusters_tender(0)
,fEventCuts()
,fCuts(0)
,fIsMC(kFALSE)
,fIsAOD(kTRUE)
,fCFM(0)
,fPID(0)
,fPIDqa(0)
,fOutputList(0)
,fOutputListMain(0)
,fOutputListLP(0)
,fOutputListHadron(0)
,fOutputListQA(0)
,fNoEvents(0)
,fNoEventsNTr(0)
,fMCNoEvents(0)
,fHFENoEvents(0)
,fDiffractiveType(0)
,fV0ACTrueInel(0)
,fV0TrueMinInel(0)
,fV0TrueMinInelNTr(0)
,fV0ACTriggered(0)
,fV0MinTriggered(0)
,fV0MinTriggeredNTr(0)
,fTriggerWeight()
,fVtxEtaNTr(0)
,fVtxBeforeNTrAcc(0)
,fVtxAfterNTrAcc()
,fVtxRecBeforeNTr(0)
,fVtxRecAfterNTr(0)
,fVtxWeight()
,fTrkpt(0)
,fEtaVtxZ(0)
,fSPDVtxRes(0)              
,fDiffSPDPrimVtx(0)         
,fSPDnTrAcc(0)    
,fSPDnTrCorrMax(0)
,fSPDnTrGen(0)              
,fDiffSPDMCVtx(0)           
,fnTrAccMaxGen(0)
,fnTrAccGen(0)
,fnTrAccGenTrueInel(0)
,fnTrAccMinGen(0)           
,fnTrAccMeanGen(0)          
,fnTrAccMax(0)              
,fnTrAccMin(0)              
,fnTrAccMean(0)
,fMCThrustTagged(0)
,fMCSpherTagged(0)
,fRecLPTagged(0)
,fMultCorrTagged(0)
,fNHadTagged(0)
,fNHadTaggedA(0)
,fNHadTaggedB(0)
,fNHadTaggedC(0)
,fMeanPtTagged(0)
,fMeanPtTaggedA(0)
,fMeanPtTaggedB(0)
,fMeanPtTaggedC(0)
,fMCThrustNTagged(0)
,fMCSpherNTagged(0)
,fRecLPNTagged(0)
,fMultCorrNTagged(0)
,fNHadNTagged(0)
,fNHadNTaggedA(0)
,fNHadNTaggedB(0)
,fNHadNTaggedC(0)
,fMeanPtNTagged(0)
,fMeanPtNTaggedA(0)
,fMeanPtNTaggedB(0)
,fMeanPtNTaggedC(0)
,fPt2Tagged(0)
,fPt2NTagged(0)
,fMothMCThrustTagged(0)
,fMothMCSpherTagged(0)
,fMothRecLPTagged(0)
,fMothMultCorrTagged(0)
,fMothNHadTagged(0)
,fMothMeanPtTagged(0)
,fMothMCThrustNTagged(0)
,fMothMCSpherNTagged(0)
,fMothRecLPNTagged(0)
,fMothMultCorrNTagged(0)
,fMothNHadNTagged(0)
,fMothMeanPtNTagged(0)
,fMCThrustTaggedH(0)
,fMCSpherTaggedH(0)
,fRecLPTaggedH(0)
,fMultCorrTaggedH(0)
,fNHadTaggedH(0)
,fMeanPtTaggedH(0)
,fMCThrustNTaggedH(0)
,fMCSpherNTaggedH(0)
,fRecLPNTaggedH(0)
,fMultCorrNTaggedH(0)
,fNHadNTaggedH(0)
,fMeanPtNTaggedH(0)
,fMothMCThrustTaggedH(0)
,fMothMCSpherTaggedH(0)
,fMothRecLPTaggedH(0)
,fMothMultCorrTaggedH(0)
,fMothNHadTaggedH(0)
,fMothMeanPtTaggedH(0)
,fMothMCThrustNTaggedH(0)
,fMothMCSpherNTaggedH(0)
,fMothRecLPNTaggedH(0)
,fMothMultCorrNTaggedH(0)
,fMothNHadNTaggedH(0)
,fMothMeanPtNTaggedH(0)
,fMultiplicity(0)
,fSPDMultiplicity(0)
,fRunList(0)
,fElectronTrackCuts(0)
,fElectronTrackTPCChi2(0)
,fElectronTrackTPCCrossedRows(0)
,fElectronTrackTPCNcls(0)
,fElectronTrackTPCNclsdEdx(0)
,fElectronTrackTPCFrac(0)
,fElectronTrackITSNcls(0)
,fElectronTrackITSChi2(0)
,fElectronTrackITSLayer(0)
,fElectronTrackRefit(0)
,fElectronTrackDCA(0)
,fElectronTrackITSCuts(0)
,fPhotTrackITSCuts(0)
,fHadronTrackCuts(0)
,fHadronTrackTPCNcls(0)
,fHadronTrackRefit(0)
,fHadronTrackDCA(0)
,fHadronTrackDCA_woITSAny(0)
,fHadronTrackDCA_wITSAny(0)
,fHistITSnSig(0)
,fHistTOFnSig(0)
,fHistTPCnSig(0)
,fHistTPCnSigITScut(0)
,fHistTPCnSigTOFcut(0)
,fHistTPCnSigITSTOFcut(0)
,fHistITSnSigTOFTPCcut(0)
,fCheckNHadronScaling(0)
,fCheckNPhotHadScaling(0)
,fCheckTaggedEvent(0)
,fHadContPvsPt(0)
,fHadContEtaPhiPt(0)
,fHadContTPCEtaPhiPt(0)
,fHadContPPhiEtaTPC(0)
,fHadContamination(0)
,fHadContaminationPt(0)
,fHadContMC(0)
,fHadContMCPt(0)
,fInclElecPtEta(0)
,fInclElecPtEtaWRecEff(0)
,fInclElecP(0)
,fULSElecPt(0)
,fULSElecPtWRecEff(0)
,fLSElecPt(0)
,fLSElecPtWRecEff(0)
,fInvmassLS(0)
,fInvmassULS(0)
,fInvmassMCTrue(0)
,fRecMCInvMass(0)
,fPhotMixULS(0)
,fPhotMixLS(0)
,fPhotPt1PtMTag(0)
,fPhotPt1PtMNTag(0)
,fPhotPt1Pt2(0)
,fPhotPt1Pt2Only(0)
,fPhotPt1Pt2Corr(0)
,fPhotPt1Pt2MC(0)
,fPhotPt1RecPt2(0)
,fPhotPt1RecPt2Corr(0)
,fPhotPt1RecPt2Rec(0)
,fPhotPt1RecPt2RecCorr(0)
,fPhotPt1Pt2Rec(0)
,fPhotPt1Pt2RecCorr(0)
,fPhotPt2MCRec(0)
,fPhotPt1E(0)
,fPhotPt1Pt2E(0)
,fPhotPt1ECorr(0)
,fPhotPt1Pt2ECorr(0)
,fPhotPt1Mass(0)
,fPhotPt1Mom(0)
,fOpeningAngleLS(0)
,fOpeningAngleULS(0)
,fCheckLSULS(0)
,fTagEtaPt1Pt2(0)
,fTagEtaPhiPt(0)
,fTagEtaZvtxPt(0)
,fTagEtaPhiPtwW(0)
,fTagEtaZvtxPtwW(0)
,fNonTagEtaPt1Pt2(0)
,fNonTagEtaPhiPt(0)
,fNonTagEtaZvtxPt(0)
,fNonTagEtaPhiPtwW(0)
,fNonTagEtaZvtxPtwW(0)
,fTagMotherPt(0)
,fTagEffInclMult(0)
,fTagEffULSMult(0)
,fTagEffInclBGMult(0)
,fTagEffULSBGMult(0)
,fTagTruePairsMult(0)
,fTagEffInclMultWoW(0)
,fTagEffULSMultWoW(0)
,fTagEffInclBGMultWoW(0)
,fTagEffULSBGMultWoW(0)
,fTagTruePairsMultWoW(0)
,fTagEffInclMultWoWS(0)
,fTagEffULSMultWoWS(0)
,fTagEffInclBGMultWoWS(0)
,fTagEffULSBGMultWoWS(0)
,fTagTruePairsMultWoWS(0)
,fTagEffIncl(0)
,fTagEffLS(0)
,fTagEffULS(0)
,fTagTruePairs(0)
,fTagEffInclWoWeight(0)
,fTagEffLSWoWeight(0)
,fTagEffULSWoWeight(0)
,fTagTruePairsWoWeight(0)
,fCorrectPiontoData()
,fCorrectEtatoData()
,fBgWeight()
,fHadRecEff()
,fEleRecEff()
,fSPDConfig(0)
,fSPDConfigHist()
,fSPDConfigProfiles()
,fSPDnTrAvg(0)
,fNonTagCorr()
,fAssPtHad_Nbins(0)
,fAssPtHad_Xmin()
,fAssPtHad_Xmax()
,fAssPtElec_Nbins(0)
,fAssPtElec_Xmin()
,fAssPtElec_Xmax()
,fElecTrigger(0)
,fInclElecPhi(0)
,fInclElecEta(0)
,fInclElecPhiEta(0)
,fULSElecPhi(0)
,fLSElecPhi(0)
,fElecDphi(0)
,fULSElecDphi(0)
,fLSElecDphi(0)
,fULSElecDphiDiffMethod(0)
,fLSElecDphiDiffMethod(0)
,fNoPartnerNoT(0)
,fNoPartnerNoTPt2(0)
,fTPartnerNoT(0)
,fTPartnerNoTPt2(0)
,fElecHadTrigger(0)
,fElecHadTriggerLS(0)
,fElecHadTriggerULS(0)
,fElecHadTriggerLSNoP(0)
,fElecHadTriggerULSNoP(0)
,fElecHadTriggerLSNoPCorr(0)
,fElecHadTriggerULSNoPCorr(0)
,fElecHadTriggerULSNoPCorrTrue(0)
,fHadContTrigger(0)
,fHadElecTrigger(0)
,fNonElecHadTrigger(0)
,fHadNonElecTrigger(0)
,fInclElecHa(0)
,fLSElecHa(0)
,fULSElecHa(0)
,fULSElecHaTrue(0)
,fSignalElecHa(0)
,fBackgroundElecHa(0)
,fMCElecHaHadron(0)
,fElecHaHa(0)
,fElecHaLSNoPartner(0)
,fElecHaULSNoPartner(0)
,fElecHaLSNoPartnerCorrTrue(0)
,fElecHaULSNoPartnerCorrTrue(0)
,fElecHaLSNoPartnerCorr(0)
,fElecHaULSNoPartnerCorr(0)
,fMCElecHaTruePartner(0)
,fMCElecHaNoPartner(0)
,fMCElecHaRemovedPartner(0)
,fMCElecHaTruePartnerTrigger(0)
,fMCElecHaTruePartnerTriggerWW(0)
,fMCElecHaNoPartnerTrigger(0)
,fMCElecHaNoPartnerTriggerWW(0)
,fMCElecHaRemovedPartnerTrigger(0)
,fElecHaMixedEvent(0)
,fLSElecHaMixedEvent(0)
,fULSElecHaMixedEvent(0)
,fTagHaMixedEvent(0)
,fNonTagHaMixedEvent(0)
,fElecLPTrigger(0)
,fElecLPTriggerLS(0)
,fElecLPTriggerULS(0)
,fElecLPTriggerLSNoP(0)
,fElecLPTriggerULSNoP(0)
,fElecLPTriggerULSNoPCorr(0)
,fHadContLPTrigger(0)
,fLPElecTrigger(0)
,fLPNonElecTrigger(0)
,fNonElecLPTrigger(0)
,fInclElecLP(0)
,fLSElecLP(0)
,fULSElecLP(0)
,fMCElecLPHadron(0)
,fElecLPHa(0)
,fElecLPLSNoPartner(0)
,fElecLPULSNoPartner(0)
,fElecLPLSNoPartnerCorrTrue(0)
,fElecLPULSNoPartnerCorrTrue(0)
,fElecLPLSNoPartnerCorr(0)
,fElecLPULSNoPartnerCorr(0)
,fMCElecLPTruePartner(0)
,fMCElecLPNoPartner(0)
,fMCElecLPRemovedPartner(0)
,fMCElecLPTruePartnerTrigger(0)
,fMCElecLPNoPartnerTrigger(0)
,fMCElecLPRemovedPartnerTrigger(0)
,fElecLPMixedEvent(0)
,fLSElecLPMixedEvent(0)
,fULSElecLPMixedEvent(0)
,fTagLPMixedEvent(0)
,fNonTagLPMixedEvent(0)
,fCheckMCVertex(0)
,fCheckMCPtvsRecPtHad(0)
,fCheckMCEtavsRecEtaHad(0)
,fCheckMCPhivsRecPhiHad(0)
,fMCHadPtEtaPhiVtx(0)
,fRecHadMCSecondaryCont(0)
,fRecHadMCPtEtaPhiVtx(0)
,fRecHadPtEtaPhiVtx(0)
,fRecHadPtEtaPhiVtxWRecEff(0)
,fCheckMCPtvsRecPtEle(0)
,fRecHFE(0)
,fMCElecPtEtaPhiVtx(0)
,fRecElecMCSecondaryCont(0)
,fRecElecPtEtaPhiVtx(0)
,fRecElecPtEtaPhiVtxWRecEff(0)
,fRecElecMCPtEtaPhiVtx(0)
,fMCElecPDG(0)
,fMCElecPtEtaPhiStrictVtx(0)
,fMCPi0Prod(0)
,fMCEtaProd(0)
,fMCPiPlusProd(0)
,fMCPiPlusProdV2(0)
,fMCBGProd(0)
,fMCLeadingParticle(0)
,fCompareLPRecCheck(0)
,fMCTruePoolMgr(0)
,fTrueMCHadronEventCuts(0)   
,fTrueMCHadronEventCutsZvtx(0)   
,fTrueMCHadronEventCutsZvtxMEv(0)   
,fTrueMCHadron(0)   
,fTrueMCElecHaTriggerEventCuts(0) 
,fTrueMCElecHaTrigger(0) 
,fTrueMCLPEventCuts(0)   
,fTrueMCLPEventCutsZvtx(0)   
,fTrueMCLPEventCutsZvtxMEv(0)   
,fTrueMCLP(0)   
,fTrueMCElecLPTriggerEventCuts(0) 
,fTrueMCElecLPTrigger(0) 
,fTrueElectronEta(0) 
,fRecHFEEtaWRecEff(0) 
,fTrueLPinAcceptanceEta(0) 
,fTrueLPEta(0) 
,fRecLPEta(0) 
,fTrueHadronEta(0) 
,fRecHadronEtaWRecEff(0) 
,fCompareLP(0) 
,fV0cutsESD(0)
,fV0cuts(0)
,fV0electrons(0)
,fV0pions(0)
,fV0protons(0)
,fhArmenteros(0)
,fEventsPerRun(0)
,fTRDnTrackRun(0)
,fV0tags(0)
,fTRDEtaPhi(0)
,fTRDNTracklets(0)
,fTRDV0NTracklets(0)
,fTRDSpectra(0)
,fTRDV0Spectra(0)
,fTRDMCSpectra(0)
{

    // Default constructor
    
    DefineInput(0, TChain::Class());
    //DefineInput(1, TList::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2,TList::Class());
    DefineOutput(3,TList::Class());
    DefineOutput(4,TList::Class());
    DefineOutput(5,TList::Class());
}
//_________________________________________

AliAnalysisTaskHaHFECorrel::~AliAnalysisTaskHaHFECorrel()
{
  //Destructor
  delete fOutputList;
  delete fOutputListMain;
  delete fOutputListLP;
  delete fOutputListHadron;
  delete fOutputListQA;
  delete fPID;
  delete fCFM;
  delete fCuts;
  delete fPIDqa;
  delete fTracks_tender;
  delete fCaloClusters_tender;
  delete fPoolMgr;
  delete fV0cuts;
  delete fV0cutsESD;
  delete fRunList;
  //delete fTrackCuts;
  //delete fAssTrackCuts;
}
//_________________________________________

void AliAnalysisTaskHaHFECorrel::UserExec(Option_t*)
{
  //Main loop
  //Called for each event 



  // create pointer to event
  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  if(!fVevent){
    printf("ERROR: fVEvent not available\n");
    return;
  }
  
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!(fESD || fAOD))
    {
      printf("ERROR: fESD & fAOD not available\n");
      return;
    }
  
  // Tender
  if(fUseTender){
    fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AODFilterTracks"));
    fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));
    if (!fTracks_tender || !fCaloClusters_tender) return;
  }
 
  
  // Initialize MC
  if(fIsMC) {
    if (fIsAOD) {
      AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
      fMC = eventHandler->MCEvent();
      fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
      if (!fMCheader) AliError("fMCheader could not be initialised");
      }
    else {
      AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      fMC = eventHandler->MCEvent();
      fStack = fMC->Stack();
    }
    if (!fMC) {
	cout << "No fMC" << endl;
    }
  }

    // get MCNCharged particle and observabls for TagEffCrossCheck
   Int_t nTrMCAcc=0;
   Double_t SumMCHadronsPt[4]={0., 0., 0., 0.}; // >all, 0.5-2, 2-5, 5-10;
   Double_t AverageMCPt[4]={0., 0., 0., 0.};
   Int_t MCnHadrons[4]={0, 0, 0, 0};
   TObjArray* MCHadrons;
   if (fIsMC) {
     if (fOneTimeCheck)MCHadrons = new TObjArray(1000, 0);
     for (Int_t i=0; i<fMC->GetNumberOfTracks(); i++) {
       AliAODMCParticle *mcPart  = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(i));
       if (mcPart->Charge()==0) continue;
       if (!mcPart->IsPhysicalPrimary()) continue;
       if (TMath::Abs(mcPart->Eta())<1.) {
 	 nTrMCAcc++;
	 if (fOneTimeCheck) MCHadrons->Add(mcPart);
	 SumMCHadronsPt[0]+=mcPart->Pt();
	 MCnHadrons[0]++;
	 Double_t ptHadron = mcPart->Pt();
	 if (ptHadron>0.5 && ptHadron<=2) {
	   SumMCHadronsPt[1]+=mcPart->Pt();
	   MCnHadrons[1]++;
	 }
	 if (ptHadron>2 && ptHadron<=5) {
	   SumMCHadronsPt[2]+=mcPart->Pt();
	   MCnHadrons[2]++;
	 }
	 if (ptHadron>5 && ptHadron<=10) {
	   SumMCHadronsPt[3]+=mcPart->Pt();
	   MCnHadrons[3]++;
	 }
       }
     }	  
   }
   for (Int_t j=0; j<4; j++)  if (MCnHadrons[j]>0) AverageMCPt[j] = SumMCHadronsPt[j]/(1.*MCnHadrons[j]);



  Int_t LPinAccBeforeEventCuts=0, LPBeforeEventCuts=0;
  Int_t EventHasElectroninPtBin[fAssPtElec_Nbins];

  // MC Truth correlations and find HFEs (for EventBias)
  for (Int_t i=0; i<fAssPtElec_Nbins; i++) EventHasElectroninPtBin[i]=0;
  if (fIsMC) {
    TObjArray* MCTrueRedTracks=new TObjArray;
    MCTrueRedTracks->SetOwner(kTRUE);
    if (fMCTrueCorrelation) MCTruthCorrelation(MCTrueRedTracks, kFALSE, 0, 0, nTrMCAcc,  LPinAccBeforeEventCuts, LPBeforeEventCuts, 1.) ; // EventWeight are 1 
    if (fMCTrueCorrelation && MCTrueRedTracks->GetEntriesFast()>0) {
      for (Int_t i=0; i<MCTrueRedTracks->GetEntriesFast(); i++) {
	AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) MCTrueRedTracks->At(i);
	AliAODMCParticle* MCTrack = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(RedTrack->GetLabel()));
	if (abs(MCTrack->GetPdgCode())!=11) continue;
	Double_t pt = MCTrack->Pt();
	Double_t eta = MCTrack->Eta();
	if (eta<fMinElectronEta || eta>fMaxElectronEta) continue;
	for (Int_t j=0; j<fAssPtElec_Nbins; j++) {
	  if (pt>fAssPtElec_Xmin[j] && pt<fAssPtElec_Xmax[j])  EventHasElectroninPtBin[j]=1;
	}
      }
      delete MCTrueRedTracks;
    }
  }

 
   
  // Get Vertex 
  const AliVVertex *pVtx=0;  
  const AliVVertex *spdVtx=0;
  if (fAOD) {
    pVtx =   fAOD->GetPrimaryVertex();
    spdVtx = fAOD->GetPrimaryVertexSPD();
  }
  else if (fESD) {
    pVtx =   fESD->GetPrimaryVertex();
    spdVtx = fESD->GetPrimaryVertexSPD();
  }

    

  //  GetSPDTracklets
  AliAODTracklets* SPDtracklets= ((AliAODEvent*)fAOD)->GetTracklets(); 
  Int_t nTr=SPDtracklets->GetNumberOfTracklets();
  Int_t nTr15Acc = 0;
  Int_t nTrAcc=0;
  for(Int_t iTr=0; iTr<nTr; iTr++){
    Double_t theta=SPDtracklets->GetTheta(iTr);
    Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
    if(TMath::Abs(eta)<1.) nTrAcc++;
    if(TMath::Abs(eta)<1.5) nTr15Acc++;
  }

  

  // Corrected SPDTracklets
  Double_t RefMinSPD, RefMaxSPD, RefMeanSPD;
  //TString RefMaxSPDString =RefMaxSPD->GetTitle();
  if (fIsMC) { // weventweights
    RefMinSPD=8.32;
    RefMaxSPD=11.33;
    RefMeanSPD=11.07;

    if (fUseEPOS) {
      RefMaxSPD=11.69;
      RefMinSPD=8.11;
      RefMeanSPD=10.76;
    }
  }
  else {
    RefMinSPD=8.28;
    RefMaxSPD=11.33;
    RefMeanSPD=11.11;
  }
  //fSPDnTrackAvg - only temporary, adjust per run, period  mc


  fRunNumber = fVevent->GetRunNumber();
  Int_t SPDConfigBin=-1;
  // implement FindFixBin (Root6) for backward compatibility to Root5 (instead of FindBin were bin is added)
  //  cout << "fSPDConfigHist " << &fSPDConfigHist << endl;
  TObjString *obj = (TObjString*)  fSPDConfigHist.GetXaxis()->GetLabels()->FindObject(Form("%i", fRunNumber));
  if (obj) SPDConfigBin= (Int_t)obj->GetUniqueID();
  else return;  // exit event if unknow SPD configuration
  if (fSPDConfig!=(Int_t) fSPDConfigHist.GetBinContent(SPDConfigBin)) {
    // for (Int_t i=1; i<fSPDConfigHist.GetXaxis()->GetNbins(); i++) cout << i << "\t" << fSPDConfigHist.GetXaxis()->GetBinLabel(i) << endl;
    //    cout<< SPDConfigBin << endl;
    fSPDConfig = fSPDConfigHist.GetBinContent(SPDConfigBin);
    //cout << fSPDConfig << endl;
    Int_t SPDConfigProfBin;
    TObjString *obj = (TObjString*)  fSPDConfigProfiles.GetXaxis()->GetLabels()->FindObject(Form("%i", fSPDConfig));
    if (obj) SPDConfigProfBin= (Int_t)obj->GetUniqueID();
    else  SPDConfigProfBin=1;// select first bin if configuration does not exist
    // cout << "fSPDCOnfigProfiles " << &fSPDConfigProfiles << endl;
    // cout << SPDConfigProfBin << endl;
    fSPDConfigProfiles.GetXaxis()->SetRange(SPDConfigProfBin, SPDConfigProfBin);
    TH2F* Configuration = 0;
    Configuration = (TH2F*)fSPDConfigProfiles.Project3D("zy");
    Configuration->SetName(Form("ConfigHist2D_%s", GetName()));
    fSPDConfigProfiles.GetXaxis()->SetRange(0, 0);
    // cout << "Configuration " << Configuration << endl;
    //  cout<< "GetTaskName" << GetName() << endl;
    if (fSPDnTrAvg!=0) {
      //cout << "Delete " << fSPDnTrAvg << endl;
      delete fSPDnTrAvg;
      fSPDnTrAvg=0;
    }
    if (Configuration!=0) {
      fSPDnTrAvg = (TProfile*) Configuration->ProfileX(Form("Prof_%i_%s", fSPDConfig, GetName()), 2, 1000); // neglecting 0 bin
      RefMaxSPD = fSPDnTrAvg->GetMaximum();
      // cout << "RefMaxSPD " << RefMaxSPD << endl;
      delete Configuration;
    }
    else{
      cout << "ProfileHist not created" << endl;
      return;
    }
    Configuration=0;
    //  cout <<  fRunNumber << "\t" << fSPDConfig << "\t" << fSPDnTrAvg << endl;
    //   for (Int_t i=1; i<fSPDnTrAvg->GetXaxis()->GetNbins(); i++) cout << fSPDnTrAvg->GetBinContent(i) << endl;
  }
  Double_t nTrAccCorrMin, nTrAccCorrMax, nTrAccCorrMean;
  if (fSPDnTrAvg!=0) {
     nTrAccCorrMin=AliVertexingHFUtils::GetCorrectedNtracklets(fSPDnTrAvg,nTrAcc*1.,spdVtx->GetZ(),RefMinSPD); 
     nTrAccCorrMax=AliVertexingHFUtils::GetCorrectedNtracklets(fSPDnTrAvg,nTrAcc*1.,spdVtx->GetZ(),RefMaxSPD); 
     nTrAccCorrMean=AliVertexingHFUtils::GetCorrectedNtracklets(fSPDnTrAvg,nTrAcc*1.,spdVtx->GetZ(),RefMeanSPD);
  }
  else {
    nTrAccCorrMin= nTrAcc*1.;
    nTrAccCorrMax=nTrAcc*1.;
    nTrAccCorrMean=nTrAcc*1.;
  }


  
  // Multiplicity estimates
   // Old get Multitplicity
  Double_t fMultV0Per, fMultSPDPer, fMultV0Tot, fMultSPD;
  Float_t mult = nTrAccCorrMax*1.;
  fMultV0Per = 1.;
  fMultSPDPer =1.;
  fMultV0Tot = 1.;
  fMultSPD = 1;
  
  fMultSelection = (AliMultSelection * ) fVevent->FindListObject("MultSelection");
  if (fMultSelection) {
    fMultV0Per = fMultSelection->GetMultiplicityPercentile("V0M", kFALSE); // Method, Embed Event selection (kFALSE)
    fMultSPDPer = fMultSelection->GetMultiplicityPercentile("SPDTracklets", kFALSE); // Method, Embed Event selection (kFALSE)
  }

  AliVVZERO* AODV0 = fVevent->GetVZEROData();
  Float_t multV0A=AODV0->GetMTotV0A();
  Float_t multV0C=AODV0->GetMTotV0C();
  fMultV0Tot=multV0A+multV0C;
  Float_t minV0=TMath::Min(multV0A, multV0C);
  Double_t TriggerWeight = GetTriggerWeight(fRunNumber, minV0, 1.*nTrAcc);
  Double_t VtxWeight = GetVtxWeight(fRunNumber, 1.*nTrAcc);
  Double_t EventWeight = TriggerWeight*VtxWeight;


  if (kTRUE) { // remove runs with different SPD configuration (active staves) and too low statistic to correct for
    // for Mult estimate
    if (fRunNumber==254983 || fRunNumber==254984) return; //SPDConfig 33440011
    if (fRunNumber==257960) return; //SPDConfig 22441122
    if (fRunNumber==263331 || fRunNumber==263332) return; //SPDConfig 22442211
    if (fRunNumber==275149) return; //SPDConfig 4445511;
    if (fRunNumber==277194 || fRunNumber==282607) return; //SPDConfig 22551122
    if (fRunNumber==280405) return; //SPDConfig 22661111
    if (fRunNumber==281920) return; //SPDConfig 22551121

    // because of HFE efficiency (SPD same for 16d to 17f(wo17c), 17ghlmor, 17ijk
    if (fRunNumber==254983 || fRunNumber==254984) return; //16h  more (in)active modules in inner/outer layer (already removed by Mult cut)
    if (fRunNumber==257138 || fRunNumber==257960) return; // 16k inner/outer layer
    if (fRunNumber>=259747 && fRunNumber<=259860) return; // 16l inner layer
    if (fRunNumber==263331 || fRunNumber==263332) return; // 16o outer layer (see above)
    if (fRunNumber>=270581 && fRunNumber<=270667) return; // complete period 17c
    if (fRunNumber==271774) return; // 17g outer layer
    if (fRunNumber==272690 || (fRunNumber>=272747 && fRunNumber<=273100)); // cuts away most of 17h or separate efficiency?
    if (fRunNumber==275149) return; // 17k (see above)
    if (fRunNumber==277194) return; // 17l (see above)
    if (fRunNumber==281920) return; // 17o (see above)
    if (fRunNumber==282607) return; // 17r

    // because of V0 modified voltage
    if (fRunNumber==274594) return; // 17j
    
  }

  // cout << "fNoEvents " << fNoEvents << endl;
  
  fNoEvents->Fill(0);

  if (fIsMC && fOneTimeCheck) fDiffractiveType->Fill(0., 1.*fMCheader->GetEventType(), nTrMCAcc);
  Bool_t IsTrueInelastic=kTRUE;
  if (fIsMC) {
    IsTrueInelastic=fEventCuts.IsTrueINELgtZero(fVevent); // MC True INEL>0
    fHFENoEvents->Fill(nTrMCAcc, 0 );
    fMCNoEvents->Fill(nTrMCAcc, 0);
    Double_t mcVtx[3];
    fMCheader->GetVertex(mcVtx);
    if (IsTrueInelastic) {
      fNoEvents->Fill(1);
      fHFENoEvents->Fill(nTrMCAcc, 1);
      fMCNoEvents->Fill(nTrMCAcc, 1);
      if (fOneTimeCheck) {
	fDiffractiveType->Fill(1., 1.*fMCheader->GetEventType(), minV0);      
	fV0ACTrueInel->Fill(multV0A, multV0C);
      }
      fV0TrueMinInel->Fill(minV0, mcVtx[2]);
      fV0TrueMinInelNTr->Fill(minV0, nTrAcc, mcVtx[2]);
      if (TMath::Abs(mcVtx[2])<10) fnTrAccGenTrueInel->Fill(nTrAcc, nTrMCAcc);
    }
    if (IsTrueInelastic && EventHasElectroninPtBin[0]==1) fHFENoEvents->Fill(nTrMCAcc, 2);    
  } 


  fEventCuts.SetManualMode();
  fEventCuts.fRequireTrackVertex = false; // not in default pp cuts
  fEventCuts.fMinVtz = -10.f;
  fEventCuts.fMaxVtz = 10.f;
  if (fVarZVTXCut==1) fEventCuts.fMaxVtz = 0.52;
  if (fVarZVTXCut==2) fEventCuts.fMinVtz = 0.52;
  
  fEventCuts.fMaxDeltaSpdTrackAbsolute = 0.5f;
  fEventCuts.fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
  fEventCuts.fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
  fEventCuts.fMaxResolutionSPDvertex = 0.25f;
  fEventCuts.fTriggerMask = AliVEvent::kINT7;
  fEventCuts.fRejectDAQincomplete = true;
  fEventCuts.fRequiredSolenoidPolarity = 0; // chose polarity (0 cut off)

  fEventCuts.fUseMultiplicityDependentPileUpCuts= true;
  //fEventCuts.fSPDpileupMinContributors = 3;
  fEventCuts.fSPDpileupMinZdist = 0.8;
  fEventCuts.fSPDpileupNsigmaZdist = 3.;
  fEventCuts.fSPDpileupNsigmaDiamXY = 2.;
  fEventCuts.fSPDpileupNsigmaDiamZ = 5.;
  fEventCuts.fTrackletBGcut = true;

  fEventCuts.SelectOnlyInelGt0(false); // this cuts away all inel>0 with not reconstructed tracks
  
  fEventCuts.fUseVariablesCorrelationCuts = false; // maybe remove this cut, small effect, not recommended
  fEventCuts.fFB128vsTrklLinearCut[0] = 32.077;
  fEventCuts.fFB128vsTrklLinearCut[1] = 0.932;

  Bool_t EventIsAccepted = fEventCuts.AcceptEvent(fVevent);

  if (fEventCuts.CheckNormalisationMask(AliEventCuts::kTriggeredEvent)) { // AliPhysicsSelction: kINT7 + Pileup
    fNoEventsNTr->Fill(0.,  nTrMCAcc>0 ? 1. : 0.);
    if (fIsMC && IsTrueInelastic) fNoEventsNTr->Fill(1., nTrMCAcc>0 ? 1. : 0.);

    if (nTr>0.5) { // All events having at leas one SPD tracklet (not restricted to be eta<1)
      fNoEventsNTr->Fill(2., nTrMCAcc>0 ? 1. : 0.);
      if (fIsMC && IsTrueInelastic) fNoEventsNTr->Fill(3., nTrMCAcc>0 ? 1. : 0.);
      if (nTrAcc>0.5) {
	fNoEventsNTr->Fill(4., nTrMCAcc>0 ? 1. : 0.);
	if (fIsMC && IsTrueInelastic) fNoEventsNTr->Fill(5., nTrMCAcc>0 ? 1. : 0.);
      }
      
      if (fIsMC) {
	if (fOneTimeCheck) {
	  fDiffractiveType->Fill(2., 1.*fMCheader->GetEventType(), minV0);
	  if (nTrAcc>0.5)  fDiffractiveType->Fill(4., 1.*fMCheader->GetEventType(), minV0); // with rec INEL>0
	  if (IsTrueInelastic) {
	  fDiffractiveType->Fill(3., 1.*fMCheader->GetEventType(), minV0);
	  fDiffractiveType->Fill(5., 1.*fMCheader->GetEventType(), minV0, TriggerWeight);
	  fV0ACTriggered->Fill(multV0A, multV0C);
	  }
	}
	fMCNoEvents->Fill(nTrMCAcc,2, TriggerWeight);
	if (IsTrueInelastic) {
	  if (EventHasElectroninPtBin[0]==1) fHFENoEvents->Fill(nTrMCAcc, 3, TriggerWeight);
	  fMCNoEvents->Fill(nTrMCAcc, 3, TriggerWeight);
	  Double_t mcVtx[3];
	  fMCheader->GetVertex(mcVtx);
	  fV0MinTriggered->Fill(minV0, mcVtx[2]);
	  fV0MinTriggeredNTr->Fill(minV0, nTrAcc, mcVtx[2]);
	}
      }
      else {
	if (fOneTimeCheck) fV0ACTriggered->Fill(multV0A, multV0C);
	fV0MinTriggered->Fill(minV0, pVtx->GetZ());
	fV0MinTriggeredNTr->Fill(minV0, nTrAcc, pVtx->GetZ());
      }
      fNoEvents->Fill(3, TriggerWeight);  
    }
  }
  
  if (fEventCuts.CheckNormalisationMask(AliEventCuts::kPassesNonVertexRelatedSelections) && nTr>0.5) { //Trigger+nTr +, Pileup, DAQ, BField
    if (fIsMC) {
      if (EventHasElectroninPtBin[0]==1) fHFENoEvents->Fill(nTrMCAcc, 4, TriggerWeight);
      fMCNoEvents->Fill(nTrMCAcc, 4, TriggerWeight);
      Double_t mcVtx[3];
      fMCheader->GetVertex(mcVtx);
      fVtxBeforeNTrAcc->Fill(nTrAcc*1., mcVtx[2], TriggerWeight);
    }
    else{
      fVtxBeforeNTrAcc->Fill(nTrAcc*1., 0., TriggerWeight);
    }
    fNoEvents->Fill(4, TriggerWeight);
    //  fnTrAccGenTrueInelTrig->Fill(nTrAcc, nTrMCAcc);
  }

  if (fEventCuts.CheckNormalisationMask(AliEventCuts::kHasReconstructedVertex) && nTr>0.5) { // + VertexExistence, VertexQuality
    if (fIsMC) {// && IsTrueInelastic) {
      if (EventHasElectroninPtBin[0]==1) fHFENoEvents->Fill(nTrMCAcc, 5, EventWeight);
      fMCNoEvents->Fill(nTrMCAcc, 5, EventWeight);
    }
    fNoEvents->Fill(5, EventWeight);
    for(Int_t iTr=0; iTr<nTr; iTr++){
      Double_t theta=SPDtracklets->GetTheta(iTr);
      Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
      if (fOneTimeCheck) fVtxEtaNTr->Fill(pVtx->GetZ(), eta);
    }
    //    fnTrAccGenTrueInelVtxQA->Fill(nTrAcc, nTrMCAcc);
  }

  if (fEventCuts.CheckNormalisationMask(AliEventCuts::kPassesAllCuts) && nTr>0.5) { // + VertexPosition
    if (fIsMC){ // && IsTrueInelastic) {
      if (EventHasElectroninPtBin[0]==1) fHFENoEvents->Fill(nTrMCAcc, 6, EventWeight);
      fMCNoEvents->Fill(nTrMCAcc, 6, EventWeight);
    }
    fNoEvents->Fill(6, EventWeight);
    fVtxAfterNTrAcc->Fill(nTrAcc*1., pVtx->GetZ(), TriggerWeight);
    //fnTrAccGenTrueInelVtxEx->Fill(nTrAcc, nTrMCAcc);
  }
    
  // EventCuts
  if(!EventIsAccepted || !(nTr>0.5) ) {
    PostData(1, fOutputList);
    PostData(2, fOutputListMain);
    return;
  }
 
  Double_t pVtxZ = -999.;
  pVtxZ = pVtx->GetZ();
  if(TMath::Abs(pVtxZ)>10. || TMath::Abs(spdVtx->GetZ())>10. ){ // minor case in which the spd vertex is out of range, while track vertex is in range
    cout << pVtxZ << "\t spd vtx " << spdVtx->GetZ() << endl;
    PostData(2, fOutputListMain);
    return;
  }
  fNoEvents->Fill(7, TriggerWeight*VtxWeight);
  if (fIsMC && IsTrueInelastic) {
    if (EventHasElectroninPtBin[0]==1) fHFENoEvents->Fill(nTrMCAcc, 7, EventWeight);
    fMCNoEvents->Fill(nTrMCAcc, 7, EventWeight);
  }

  UInt_t fSelectMask = fInputHandler->IsEventSelected(); // should have no impact
  Bool_t isINT7selected = fSelectMask & AliVEvent::kINT7;
  if (!isINT7selected){
    PostData(2, fOutputListMain);
    printf("Event not selected \n");
    return;
  }
  fNoEvents->Fill(8, EventWeight);

  // Perform Event Bias  // Find MotherKinks
  Int_t *listofmotherkink=0;
  Int_t nMotherKink=0;
  if (fIsAOD) {
    Int_t nVertices = 1;
    nVertices = fAOD->GetNumberOfVertices();
    listofmotherkink = new Int_t[nVertices];
      for(Int_t ivertex=0; ivertex < nVertices; ivertex++) {
      AliAODVertex *vertex = fAOD->GetVertex(ivertex);
      if(!vertex) continue;
      if(vertex->GetType()==AliAODVertex::kKink) {
	AliAODTrack *mother = (AliAODTrack *) vertex->GetParent();
	if(!mother) continue;
	Int_t idmother = mother->GetID();
	listofmotherkink[nMotherKink] = idmother;
	nMotherKink++;
      }
    }
  }


  if (fMinPtEvent > 0.1 || fMaxPtEvent <100) { // currently no used (default range 1-999
    if (!PassEventBias(pVtx,nMotherKink,listofmotherkink, EventWeight)) {
      delete [] listofmotherkink;
      PostData(2, fOutputListMain);
      return;
    }
  }
  fNoEvents->Fill(9, EventWeight);

    //  cout << "nTrAcc " << nTrAcc << "\t" << nTrAccCorrMin << endl;
  if (nTrAccCorrMax<fMinNTr || nTrAccCorrMax>fMaxNTr) {
    delete [] listofmotherkink;
    PostData(1, fOutputList);
    PostData(2, fOutputListMain);
    PostData(5, fOutputListQA);
    return;
  }
  fNoEvents->Fill(10, EventWeight);
  

  // Initialize PID Resonse
  fpidResponse = fInputHandler->GetPIDResponse();
  if(!fpidResponse){
    AliDebug(1, "Using default PID Response");
    fpidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class());
  }  

  // Vertex Control hists
  Double_t covSPD[6];
  spdVtx->GetCovarianceMatrix(covSPD);
  fSPDVtxRes->Fill(spdVtx->GetNContributors(), TMath::Sqrt(covSPD[5]));
  Double_t DiffVertexZ = spdVtx->GetZ()-pVtx->GetZ();
  fDiffSPDPrimVtx->Fill(spdVtx->GetNContributors(), (DiffVertexZ));


  // Multiplicity Histogram
  if (nTr>0.5) {
    fSPDnTrAcc->Fill(spdVtx->GetZ(), nTrAcc, EventWeight);
    fSPDnTrCorrMax->Fill(spdVtx->GetZ(), 1.*nTrAccCorrMax, EventWeight);
  }
  
  Double_t fillSparse[3]={spdVtx->GetZ(), 1.*nTrAcc,1.};
  fillSparse[2]=1.*nTrAccCorrMax;
  if (nTr>0.5) fnTrAccMax->Fill(fillSparse, EventWeight);
  // fillSparse[2]=1.*nTrAccCorrMin;
  // fnTrAccMin->Fill(fillSparse);
  // fillSparse[2]=1.*nTrAccCorrMean;
  // fnTrAccMean->Fill(fillSparse);

 
  if (fIsMC) {
    Double_t mcVtx[3];
    fMCheader->GetVertex(mcVtx);
    DiffVertexZ = spdVtx->GetZ()-mcVtx[2];
    fDiffSPDMCVtx->Fill(spdVtx->GetNContributors(), DiffVertexZ);
    if (nTr>0.5) fSPDnTrGen->Fill(spdVtx->GetZ(), nTrMCAcc,EventWeight);
    
    fillSparse[1]=1.*nTrMCAcc;
    fillSparse[2]=1.*nTrAccCorrMax;
    if (nTr>0.5) fnTrAccMaxGen->Fill(fillSparse, EventWeight);
    fillSparse[2]=1.*nTrAcc;
    if (nTr>0.5) fnTrAccGen->Fill(fillSparse, EventWeight);

    // fillSparse[2]=1.*nTrAccCorrMin;
    // fnTrAccMinGen->Fill(fillSparse);
    // fillSparse[2]=1.*nTrAccCorrMean;
    // fnTrAccMeanGen->Fill(fillSparse);
  }
 
  /*
  AliVMultiplicity* AliMult = fVevent->GetMultiplicity();
  if (AliMult) fMultSPD=AliMult->GetNumberOfTracklets();  
  if (fIsAOD)  AliAODTracklets *AODtracklets = ((AliAODEvent*)fAOD)->GetTracklets();
  mult=fMultSPDPer;
  Double_t fillSparse[4]={fMultV0Per, fMultSPDPer, fMultV0Tot, fMultSPD};
  fMultiplicity->Fill(fillSparse); 
  fSPDMultiplicity->Fill(fMultSPDPer, fMultSPD, pVtxZ);
  */
    
  // Efficiency Corrections
  if(fIsMC) {
    if (fIsAOD) {
      MCEfficiencyCorrections(pVtx, EventWeight); //  Electron reconstruction, Hadron reconstruction
      // TList *lh=fMCheader->GetCocktailHeaders();
      // Int_t nh=lh->GetEntries();  
      /* for(Int_t i=0;i<nh;i++)	{
	  AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(i);
	  TString genname=gh->GetName();
	  Int_t npart=gh->NProduced();
	  cout << genname << "with " << npart << " Particles" << endl;
	  } */
    }
    else {
      //  MCEfficiencyCorrections(pVtx); //  Electron reconstruction, Hadron reconstruction
    }
  }

  // Cross check variables for event bias
  Double_t ThrustVar[2]={-1., -1.};
  Double_t SpherVar=-1.;
  if (fIsMC && fOneTimeCheck) {
    SpherVar = Sphericity(MCHadrons, 1, 0.5);
    MCHadrons->Clear();
    delete MCHadrons;
  }
 

  ///////////////////////
  // Preparational Tasks
  ///////////////////////



  // List of HFE for analysis and mixed event
  TObjArray* RedTracksHFE = new TObjArray;
  RedTracksHFE->SetOwner(kTRUE);

  // FirstLoop: Find leading particle (LP) and HFE (check for LS/ULS partner)
  AliVTrack* LPtrack;
  fLParticle=kFALSE;
  Bool_t EvContainsTaggedPhot=kFALSE;
  Bool_t EvContainsNonTaggedPhot=kFALSE;

  PhotULSLSElectronAcceptance(pVtx, mult,  nMotherKink, listofmotherkink, EventWeight);
  LPtrack=FindLPAndHFE(RedTracksHFE, pVtx, nMotherKink,listofmotherkink, mult, EvContainsTaggedPhot, EvContainsNonTaggedPhot, EventWeight);
  if (fLParticle) {
    //    cout << "LPTrackPt " << LPtrack->Pt() << endl;
    if (LPtrack->Pt()>=1000) return;
  }

  // Control hists for TagEff
  if (fIsMC) {
    Double_t MotherPt=0; // MC
    for (Int_t i=0; i<RedTracksHFE->GetEntriesFast(); i++) {
      AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(i);
      AliAODMCParticle* MCTrack = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(RedTrack->GetLabel())); 
      if (abs(MCTrack->GetPdgCode())!=11) continue;
      Int_t MCMotherLabel = abs(MCTrack->GetMother());
      Int_t MCGMotherLabel=-999, MCGGMotherLabel=-999;
      Int_t MCMotherPDG = -999, MCGMotherPDG = -999, MCGGMotherPDG=-999;
      Bool_t IsFeedDown=kFALSE;
      Int_t Case=-999; // Case 0 gamma pion, 1 gamma pion feed down, 2 gamma eta, 3 gammaeta feed down, 4 pion 5, pion feed deon , 6, eta, 7, eta feed down
      if (MCMotherLabel>0) {
	AliAODMCParticle *MCMother = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(MCMotherLabel));  
	MCMotherPDG = abs(MCMother->GetPdgCode());
	MCGMotherLabel = abs(MCMother->GetMother());
	MotherPt=MCMother->Pt();
      }
      if (MCGMotherLabel>0) {
	AliAODMCParticle *MCGMother = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(MCGMotherLabel));  
	MCGMotherPDG = abs(MCGMother->GetPdgCode());
	MCGGMotherLabel = abs(MCGMother->GetMother());
      }
      if (MCMotherPDG==22) {
	if (MCGMotherPDG==111) Case =0;
	else if (MCGMotherPDG==221) Case=2;
	if (MCGGMotherLabel>0) {
	  AliAODMCParticle *MCGGMother = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(MCGGMotherLabel));  
	  MCGGMotherPDG = abs(MCGGMother->GetPdgCode());
	}
	if (MCGGMotherPDG>100) Case++;
      }
      else if (MCMotherPDG==111) {
	Case=4;
	if (MCGMotherPDG>100) Case++;
      }
      else if (MCMotherPDG==221) {
	Case=6;
	if (MCGMotherPDG>100) Case++;
      }
   
      if (!RedTrack->IsPhotonic()) continue;
      if (RedTrack->TruePartner()) {
	if (fLParticle) fPt2Tagged->Fill(RedTrack->Pt(), RedTrack->TruePartnerMCPt(), LPtrack->Pt(), EventWeight);
	//	fMCThrustTagged->Fill(RedTrack->Pt(), ThrustVar[0], Case);
	fMCSpherTagged->Fill(RedTrack->Pt(), SpherVar, Case, EventWeight);
	if (fLParticle) if (abs(RedTrack->ID())!=abs(LPtrack->GetID())) fRecLPTagged->Fill(RedTrack->Pt(), LPtrack->Pt(), Case, EventWeight);
	//	fMultCorrTagged->Fill(RedTrack->Pt(), nTrAccCorrMax, Case);
	fNHadTagged->Fill(RedTrack->Pt(), MCnHadrons[0], Case, EventWeight);
	//	fNHadTaggedA->Fill(RedTrack->Pt(), MCnHadrons[1], Case);
	fNHadTaggedB->Fill(RedTrack->Pt(), MCnHadrons[2], Case, EventWeight);
	//fNHadTaggedC->Fill(RedTrack->Pt(), MCnHadrons[3], Case);
	fMeanPtTagged->Fill(RedTrack->Pt(), AverageMCPt[0], Case, EventWeight);
	//	fMeanPtTaggedA->Fill(RedTrack->Pt(), AverageMCPt[1], Case);
	fMeanPtTaggedB->Fill(RedTrack->Pt(), AverageMCPt[2], Case, EventWeight);
	//fMeanPtTaggedC->Fill(RedTrack->Pt(), AverageMCPt[3], Case);
	

	//	fMothMCThrustTagged->Fill(MotherPt, ThrustVar[0]);
	fMothMCSpherTagged->Fill(MotherPt, SpherVar, EventWeight);
	if (fLParticle) if (abs(RedTrack->ID())!=abs(LPtrack->GetID())) fMothRecLPTagged->Fill(MotherPt, LPtrack->Pt(), EventWeight);
	//	fMothMultCorrTagged->Fill(MotherPt, nTrAccCorrMax);
	fMothNHadTagged->Fill(MotherPt, MCnHadrons[0], EventWeight);
	fMothMeanPtTagged->Fill(MotherPt, AverageMCPt[0],EventWeight);
	if (fLParticle) {
	  if (LPtrack->Pt()>2) {
	    //  fMCThrustTaggedH->Fill(RedTrack->Pt(), ThrustVar[0]);
	    fMCSpherTaggedH->Fill(RedTrack->Pt(), SpherVar, EventWeight);
	    if (abs(RedTrack->ID())!=abs(LPtrack->GetID())) fRecLPTaggedH->Fill(RedTrack->Pt(), LPtrack->Pt(), EventWeight);
	    //	  fMultCorrTaggedH->Fill(RedTrack->Pt(), nTrAccCorrMax);
	    fNHadTaggedH->Fill(RedTrack->Pt(), MCnHadrons[0], EventWeight);
	    fMeanPtTaggedH->Fill(RedTrack->Pt(), AverageMCPt[0], EventWeight);
	    //  fMothMCThrustTaggedH->Fill(MotherPt, ThrustVar[0]);
	    fMothMCSpherTaggedH->Fill(MotherPt, SpherVar, EventWeight);
	    if (abs(RedTrack->ID())!=abs(LPtrack->GetID())) fMothRecLPTaggedH->Fill(MotherPt, LPtrack->Pt(), EventWeight);
	    // fMothMultCorrTaggedH->Fill(MotherPt, nTrAccCorrMax);
	    fMothNHadTaggedH->Fill(MotherPt, MCnHadrons[0], EventWeight);
	    fMothMeanPtTaggedH->Fill(MotherPt, AverageMCPt[0], EventWeight);
	  }
	}
      }
      else {
	
	if (fLParticle) fPt2NTagged->Fill(RedTrack->Pt(), RedTrack->TruePartnerMCPt(), LPtrack->Pt(), EventWeight);
	
	//	fMCThrustNTagged->Fill(RedTrack->Pt(), ThrustVar[0], Case);
	fMCSpherNTagged->Fill(RedTrack->Pt(), SpherVar, Case, EventWeight);
	if (fLParticle)	if (abs(RedTrack->ID())!=abs(LPtrack->GetID())) fRecLPNTagged->Fill(RedTrack->Pt(), LPtrack->Pt(), Case, EventWeight);
	//	fMultCorrNTagged->Fill(RedTrack->Pt(), nTrAccCorrMax, Case);
	fNHadNTagged->Fill(RedTrack->Pt(), MCnHadrons[0], Case, EventWeight);
	//	fNHadNTaggedA->Fill(RedTrack->Pt(), MCnHadrons[1], Case);
	fNHadNTaggedB->Fill(RedTrack->Pt(), MCnHadrons[2], Case, EventWeight);
	//	fNHadNTaggedC->Fill(RedTrack->Pt(), MCnHadrons[3], Case);
	fMeanPtNTagged->Fill(RedTrack->Pt(), AverageMCPt[0], Case, EventWeight);
	//	fMeanPtNTaggedA->Fill(RedTrack->Pt(), AverageMCPt[1], Case);
	fMeanPtNTaggedB->Fill(RedTrack->Pt(), AverageMCPt[2], Case, EventWeight);
	//	fMeanPtNTaggedC->Fill(RedTrack->Pt(), AverageMCPt[3], Case);


	//	fMothMCThrustNTagged->Fill(MotherPt, ThrustVar[0]);
	fMothMCSpherNTagged->Fill(MotherPt, SpherVar, EventWeight);
	if (fLParticle)	if (abs(RedTrack->ID())!=abs(LPtrack->GetID())) 	fMothRecLPNTagged->Fill(MotherPt, LPtrack->Pt(), EventWeight);
	//	fMothMultCorrNTagged->Fill(MotherPt, nTrAccCorrMax);
	fMothNHadNTagged->Fill(MotherPt, MCnHadrons[0], EventWeight);
	fMothMeanPtNTagged->Fill(MotherPt, AverageMCPt[0], EventWeight);
	if (fLParticle) {
	  if (LPtrack->Pt()>2) {
	    //  fMCThrustNTaggedH->Fill(RedTrack->Pt(), ThrustVar[0]);
	    fMCSpherNTaggedH->Fill(RedTrack->Pt(), SpherVar, EventWeight);
	    if (abs(RedTrack->ID())!=abs(LPtrack->GetID())) fRecLPNTaggedH->Fill(RedTrack->Pt(), LPtrack->Pt(), EventWeight);
	    //  fMultCorrNTaggedH->Fill(RedTrack->Pt(), nTrAccCorrMax);
	    fNHadNTaggedH->Fill(RedTrack->Pt(), MCnHadrons[0], EventWeight);
	    fMeanPtNTaggedH->Fill(RedTrack->Pt(), AverageMCPt[0], EventWeight);
	    //  fMothMCThrustNTaggedH->Fill(MotherPt, ThrustVar[0]);
	    fMothMCSpherNTaggedH->Fill(MotherPt, SpherVar, EventWeight);
	    if (abs(RedTrack->ID())!=abs(LPtrack->GetID())) fMothRecLPNTaggedH->Fill(MotherPt, LPtrack->Pt(),EventWeight);
	    //  fMothMultCorrNTaggedH->Fill(MotherPt, nTrAccCorrMax);
	    fMothNHadNTaggedH->Fill(MotherPt, MCnHadrons[0], EventWeight);
	    fMothMeanPtNTaggedH->Fill(MotherPt, AverageMCPt[0], EventWeight);
	  }
	}
      }
    }
  }


  /*

  if (EvContainsTaggedPhot & !EvContainsNonTaggedPhot)  fCheckTaggedEvent->Fill(1.,LPtrack->Pt(), mult);
  else if (EvContainsNonTaggedPhot & !EvContainsTaggedPhot) fCheckTaggedEvent->Fill(2., LPtrack->Pt(), mult);
  else if (EvContainsTaggedPhot && EvContainsNonTaggedPhot) fCheckTaggedEvent->Fill(0., LPtrack->Pt(), mult);
  */

  // MC Truth correlation after event Cuts
  Int_t LPinAccAfterEventCuts=-999, LPAfterEventCuts=-999;
  AliAODTrack * LPtrackAOD=0;
  if (fLParticle) LPtrackAOD = dynamic_cast<AliAODTrack*>(LPtrack);
  Int_t LPtrackLabel=0;
  if (LPtrackAOD) LPtrackLabel=abs(LPtrackAOD->GetLabel());
  if (fIsMC && fMCTrueCorrelation) {  
    TObjArray* MCTrueRedTracks = new TObjArray;
    MCTrueRedTracks->SetOwner(kTRUE);
    MCTruthCorrelation(MCTrueRedTracks, kTRUE, LPtrackLabel, pVtx->GetZ(), mult,  LPinAccAfterEventCuts, LPAfterEventCuts, EventWeight) ;
    //    delete MCTrueRedTracks;
    if (fIsAOD && fLParticle && RedTracksHFE->GetEntriesFast()>0 && LPinAccAfterEventCuts>=0) {
      AliAODMCParticle* LPinAcc = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(LPinAccAfterEventCuts));  
      AliAODMCParticle* LP = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(LPAfterEventCuts));  
      if (LPinAcc && LP) {
	fCompareLP->Fill(LPtrack->Pt(), LPinAcc->Pt(), LP->Pt());
      }
    }
  }

  /*
  if(fTRDQA) {
    Int_t RunNumber = fVevent->GetRunNumber();
    Int_t RunIndex=0;
    while (RunNumber!=fRunList[RunIndex]) {
      RunIndex++;
      if (RunIndex>199) {
	cout << "Run not found in run list" << endl;
	break;
      }
    }
    if (fIsAOD)  FindV0CandidatesAOD(fAOD);
    else FindV0CandidatesESD(fESD);
    fEventsPerRun->Fill(RunIndex);
    TRDQA(RunIndex,pVtx,nMotherKink,listofmotherkink);
    }*/

  //////////////////////////
  // Main analysis
  ///////////////////////////


  if (fCorrLParticle && fLParticle) {
    // LP - two different functions for Same Event and mixed Event
    if (RedTracksHFE->GetEntriesFast()>0) CorrelateLP(LPtrack,pVtx, nMotherKink, listofmotherkink, RedTracksHFE, EventWeight);
    if (fMixedEvent) CorrelateLPMixedEvent(LPtrack, mult, pVtx->GetZ(), LPtrack->Pt(),  EvContainsTaggedPhot, EvContainsNonTaggedPhot, EventWeight); // condition that electron track is in event has been removed!
  }

  // Hadron - only one function for both options, as only one loop over Hadron tracks
  // Mixed event is called inside this function
  if (fCorrHadron && fLParticle) {
    if (RedTracksHFE->GetEntriesFast()>0) CorrelateHadron(RedTracksHFE, pVtx, nMotherKink, listofmotherkink, mult, LPtrack->Pt(), EventWeight);
    if (fMixedEvent) CorrelateHadronMixedEvent( mult, pVtx, LPtrack->Pt(), nMotherKink, listofmotherkink,  EvContainsTaggedPhot, EvContainsNonTaggedPhot, EventWeight);
  }
  // Electrons - currently no mixed event (>1 as two electrons are required), work in progress
  // if (RedTracksHFE->GetEntriesFast()>1) CorrelateElectron(RedTracksHFE);
 

  /////////////////////////
  //Fill Mixed event pool//
  /////////////////////////
  if (RedTracksHFE->GetEntriesFast()>0 && fLParticle) {
    AliEventPool * HFEPool = fPoolMgr->GetEventPool(mult, pVtxZ, LPtrack->Pt()); // Get the buffer associated with the current centrality and z-vtx
    //HFEPool->SetDebug(kTRUE);
    // cout << Form("\nPool found for centrality = %f, zVtx = %f, leading pt = %f end.\n", mult, pVtx->GetZ(), LPtrack->Pt()) << endl;
    if (!HFEPool)
      {
	AliFatal(Form("No pool found for centrality = %f, zVtx = %f, leading pt = %f", mult, pVtx->GetZ(), LPtrack->Pt()));
	return;
      }
    HFEPool->UpdatePool(RedTracksHFE);
    fPoolIsFilled->Fill(mult, pVtxZ,LPtrack->Pt());
  }
  else {
    delete RedTracksHFE;
  }
  
  ClearV0PIDList();
  delete [] listofmotherkink;
  
  PostData(1,fOutputList);
  PostData(2,fOutputListMain);
  PostData(3,fOutputListLP);
  PostData(4,fOutputListHadron);
  PostData(5,fOutputListQA);
}

//_________________________________________
void AliAnalysisTaskHaHFECorrel::UserCreateOutputObjects()
{      
  
  AliLog::SetClassDebugLevel("AliAnalysisTaskHaHFECorrel", 0);
  AliDebug(1, Form("This is AliDebug level %i", 1));
  AliDebug(2, Form("This is AliDebug level %i", 2));
  AliDebug(3, Form("This is AliDebug level %i", 3));
   
  //---------Output Tlist
  fOutputList = new TList();
  fOutputList->SetOwner();

  fOutputListMain = new TList();
  fOutputListMain->SetOwner();
  fOutputListLP = new TList();
  fOutputListLP->SetOwner();
  fOutputListHadron = new TList();
  fOutputListHadron->SetOwner();
  fOutputListQA = new TList();
  fOutputListQA->SetOwner(); 


  if(fIsAOD) fesdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 

  // V0 Kine cuts 
  fV0cuts = new AliAODv0KineCuts();
  fV0cutsESD = new AliESDv0KineCuts();

  // V0 PID Obj arrays
  fV0electrons = new TObjArray;
  fV0pions     = new TObjArray;
  fV0protons   = new TObjArray;

  fNoEvents = new TH1F("fNoEvents","",11,-0.5,10.5);
  fOutputListMain->Add(fNoEvents);
  TString BinLabelsEvents[11]={"AllEvents", "+IsTrueIne", "+IsHFE", "+KINT7", "+PILEUP", "+VERTEX", "+VERTEXPOS", "+SPDzVtx", "+CheckKINT7", "+PTCut", "+MultCut"};
  
  if(fIsMC) {
    fMCNoEvents = new TH2F("fMCNoEvents", "MCAccEvents; NCharged; EventCut", 150, -0.5, 149.5, 11, -0.5, 10.5);
    fOutputListMain->Add(fMCNoEvents);
    
    fHFENoEvents = new TH2F("fHFENoEvents", "HFEMCAccEvents; NCharged; EventCut", 150, -0.5, 149.5, 11, -0.5, 10.5);
    fOutputListMain->Add(fHFENoEvents);
  }
  for (Int_t i=1; i<=11; i++) {
    fNoEvents->GetXaxis()->SetBinLabel(i, BinLabelsEvents[i-1]);
    if (fIsMC) {
      fMCNoEvents->GetYaxis()->SetBinLabel(i, BinLabelsEvents[i-1]);
      fHFENoEvents->GetYaxis()->SetBinLabel(i, BinLabelsEvents[i-1]);
    }
  }

  
  fNoEventsNTr = new TH2F("fNoEventsNtr","ControlHist for 0Bin; EventSelection; NCh > 0 ?", 6, -0.5, 5.5, 2, -0.5, 1.5);
  TString BinLabelsEventsNTr[6]={"Trigger", "Trigger (wTI)", "+Ntr", "+Ntr (wTI)", "+NtrAcc", "+NtrAcc (wTI)"}; //
  for (Int_t i=1; i<=6; i++) {
    fNoEventsNTr->GetXaxis()->SetBinLabel(i, BinLabelsEventsNTr[i-1]);
  }
  fOutputListMain->Add(fNoEventsNTr);

  if (fOneTimeCheck) {
    fV0ACTriggered = new TH2F("fV0ACTriggered", "V0 kINT7 + trueInel; V0A; V0C", 200, 0, 40, 200, 0, 40);
    fOutputListMain->Add(fV0ACTriggered);
  }
  
  fV0MinTriggered = new TH2F("fV0MinTriggered", "V0min kINT7 + trueInel; Min(V0A, V0C); zVtx", 100, 0, 20, 40, -20, 20);
  fOutputListMain->Add(fV0MinTriggered);

  fV0MinTriggeredNTr = new TH3F("fV0MinTriggeredNTr", "V0min kINT7 + trueInel; Min(V0A, V0C); nTrAcc SPD; zVtx", 100, 0, 20, 25, -0.5,24.5, 4, -20, 20);
  fOutputListMain->Add(fV0MinTriggeredNTr);
  

  if (fIsMC) {
  
    if (fOneTimeCheck) {
        TString BinLabelsDiff[6]={"AllEv", "TrueInel>0", "Triggered", "Triggered+TrueInel>0", "Triggered+RecInel>0", "Trig+TrueInel+Weight"};
	fDiffractiveType = new TH3F("fDiffType", "fDiffType; EventSelection; MCEventType; minV0", 6, -0.5, 5.5, 7, 99.5, 106.5, 50, 0, 20);
	fOutputListMain->Add(fDiffractiveType);
	for (Int_t i=1; i<=6; i++) fDiffractiveType->GetXaxis()->SetBinLabel(i, BinLabelsDiff[i-1]);
 
	fV0ACTrueInel = new TH2F("fV0ACTrueInel", "V0 true inelastic; V0A; V0C", 200, 0, 40, 200, 0, 40);
	fOutputListMain->Add(fV0ACTrueInel);
    }
    
    fV0TrueMinInel = new TH2F("fV0TrueMinInel", "V0min for true inelastic; Min(V0A, V0C); zVtx", 100, 0, 20, 40, -20, 20);
    fOutputListMain->Add(fV0TrueMinInel);

    fV0TrueMinInelNTr = new TH3F("fV0TrueMinInelNTr", "V0min for true inelastic; Min(V0A, V0C); nTrAcc SPD; zVtx", 100, 0, 20, 25, -0.5, 24.5, 4, -20, 20);
    fOutputListMain->Add(fV0TrueMinInelNTr);
  }
    
 
  
  

  fVtxBeforeNTrAcc = new TH2F("fVtxBeforeNTrAcc", "VtxBeforeNTrAcc; nTrAcc; zVtx (MC)", 50, -0.5, 49.5, 40, -10, 10);
  fOutputListMain->Add(fVtxBeforeNTrAcc);

  fVtxAfterNTrAcc = new TH2F("fVtxAfterNTrAcc", "VtxAfterNTrAcc; nTrAcc; zVtx", 50, -0.5, 49.5, 40, -10, 10);
  fOutputListMain->Add(fVtxAfterNTrAcc);

  if (fOneTimeCheck) {
    fVtxEtaNTr = new TH2F("fVtxEtaTr", "#of SPDtracklets for zVtx, eta; zTrackVtx; eta", 96, -12, 12, 88, -2.2, 2.2);
    fOutputListMain->Add(fVtxEtaNTr);
    
    fEtaVtxZ = new TH2F("fEtaVtxZ", "Eta vs VtxZ after hadron track cuts; #eta; zVtx", 90, -0.9, 0.9, 100, -10, 10);
    fOutputListMain->Add(fEtaVtxZ);
  }
  /*
  fVtxRecBeforeNTr= new TH1F("fVtxRecBeforeNTr", "VtxRecBefore (ntr eta unrestricted); ntr SPD", 50, -0.5, 49.5);
  fOutputListMain->Add(fVtxRecBeforeNTr);
  
  fVtxRecAfterNTr = new TH2F("fVtxRecAfterNTr", "VtxRecAfter (ntr eta unrestricted); ntr SPD;  zVtxSPD", 50, -0.5, 49.5, 40, -10, 10);
  fOutputListMain->Add(fVtxRecAfterNTr);
  */
  
 


  PdgTable =  new TDatabasePDG();
  PdgTable->ReadPDGTable();
  const THashList* PDGParticleList = PdgTable->ParticleList();
  // cout << "PDGSize " << PDGParticleList->GetEntries() << endl;
  Int_t PDGSize = PDGParticleList->GetEntries()+20;
  cout << "PDG Size" << PDGSize << endl;
  TIter next(PDGParticleList);
  TParticlePDG *p;
  Int_t Bin=0;
  std::vector<TString> PDGLabel;
  PDGLabel.push_back("Empty");

  if (!fOneTimeCheck) { // reduced memory consumption by only filling relevant bins
                              //e, t ,ga,pi0,r0 ,k0l,pi+,r+, eta,w ,k0s ,k0 ,k+,eta',phi,psi,Y 
    std::set<Int_t> SelectBins {11,15,22,111,113,130,211,213,221,223,310,311,321,331,333,443,553};
    Bin++;
    PDGLabel.push_back(Form("<"));
    while ((p = (TParticlePDG *)next())) {
      if (SelectBins.find(p->PdgCode())!=SelectBins.end()){
	Bin++;
	PDGMap[p->PdgCode()] = Bin;
	PDGLabel.push_back(Form("%i",p->PdgCode()));
	Bin++;
	PDGLabel.push_back(Form(">%i",p->PdgCode()));
      }
      else{
	PDGMap[p->PdgCode()] = Bin;
      }
      // cout << p->PdgCode() << "\t" << Bin << "\t" << PDGLabel[Bin] << endl;
    }
  }
  else {
    while ((p = (TParticlePDG *)next())) {
      Bin++;
      PDGMap[p->PdgCode()] = Bin;
      PDGLabel.push_back(Form("%i",p->PdgCode()));
      //  cout << p->PdgCode() << "\t" << Bin << "\t" << PDGLabel[Bin] << endl;
    }
  }
  PDGSize=Bin+1;
   

  // Int_t PDGLabelShort={11,13,15, 22, 111, 113, 130,  211, 213, 221, 223, 310, 311, 313, 323,  321, 331, 333, 

    

  // General  Binning 		    
  // fAssPtHad_Nbins = 15;
  const Int_t TmpHad_Nbins = 12;
  fAssPtHad_Nbins = TmpHad_Nbins;
  Float_t TmpArrayLow[12]= {0.25, 0.25,  0.5, 0.5, 1., 1., 2,  5, 5, 10, 10, 15};
  Float_t TmpArrayUp[12]=  {1.0 ,  2.0,  2.0, 5.0, 2., 5., 5, 10,15, 15, 30, 30};
  fAssPtHad_Xmin.Set(TmpHad_Nbins, TmpArrayLow);
  fAssPtHad_Xmax.Set(TmpHad_Nbins, TmpArrayUp);
  
  const Int_t TmpElec_Nbins =12;
  fAssPtElec_Nbins = TmpElec_Nbins;
  Float_t TmpArrayELow[12]={0 ,  0.5, 0.5, 1., 1., 1., 2,  2, 4,  4,  6,  8};
  Float_t TmpArrayEUp[12]={999, 2.0, 999, 2., 4.,999, 4,999, 6,999, 10,999};
  fAssPtElec_Xmin.Set(TmpElec_Nbins, TmpArrayELow);
  fAssPtElec_Xmax.Set(TmpElec_Nbins, TmpArrayEUp);
  
  // HFE Electron
  Int_t    NBinsElectron =30;
  Double_t XBinElectronArray[31];
  for (Int_t BinE=0; BinE<=NBinsElectron;BinE++) {
    XBinElectronArray[BinE]=0.25+BinE*0.125;
    cout << XBinElectronArray[BinE] << endl;
  }
  Double_t XminElectron=0.25;
  Double_t XmaxElectron=4.;
  const Int_t    NBinsElectronRed = 19;
  Double_t XBinsElectronRed[]={0.25,0.5,0.75, 1., 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4, 5, 6, 10, 100};

  const Int_t    NBinsElectronRecEff = 22;
  Double_t XBinsElectronRecEff[]={0.25,0.375, 0.5,0.625, 0.75,0.875, 1., 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4, 5, 6, 10, 100};


  
  Int_t     NBinsHadron=200 ;
  Double_t  XminHadron=0.0;
  Double_t  XmaxHadron=100;
  const Int_t   NBinsHadRed=8;
  Double_t  XBinsHadRed[]={0.25,  0.5, 1., 2., 5., 10, 15, 30, 50}; 

  const  Int_t     NMultBinsRed=4;
  Double_t    XMultBinsRed[]={-0.5,24.5,44.5,64.5,500};

  const Int_t NMultBins=7;
  Double_t XMultBins[]={0.5, 14.5, 24.5, 34.5, 44.5,54.5,64.5,200.5};

  
  const Int_t NVertexBins = 8;
  Double_t XVertexBins[]={-10,-7,-4,-2,0,2,4,7,10};  // Quantile
  
  // Double_t XVertexBins[]={-10,-4.32,-2.32,-0.82,0.52,1.87,3.36,5.31, 10};  // Quantile


  Int_t    NBinsPhi=64;
  Int_t    NBinsEta=90;
  Double_t XminEta=-0.9;
  Double_t XmaxEta=0.9;
  Int_t    NBinsDPhi=24;
  Double_t DEtaMax=fMaxElectronEta+fMaxHadronEta;
  Double_t DEtaMin=fMinElectronEta+fMinHadronEta;
  Int_t    NBinsDEta=(DEtaMax-DEtaMin+0.005)/0.2; // 0.005 to prevent round errors
  cout << DEtaMax << "\t" << DEtaMin << endl;
  cout << "NBinsDEta" <<  NBinsDEta << endl;


    
  // Multiplicity Sparse
  fSPDVtxRes = new TH2F("fSPDVtxRes", "fSPDVtxRes; NContributors; Resolution",  200, -0.5, 199.5, 250, 0., 5.);
  fOutputListQA->Add(fSPDVtxRes);
 
  fDiffSPDPrimVtx = new TH2F("fDiffSPDPrimVtx", "fDiffSPDPrimVtx;NContributors; zSPDVtx-zTrackVtx", 200, -0.5, 199.5, 300, -3., 3.);
  fOutputListQA->Add(fDiffSPDPrimVtx);

  fSPDnTrAcc = new TH2F("fSPDnTrAcc", "fSPDnTrAcc; zSPDVtx; SPDnTr (eta<1)", 220, -11, 11, 300, -0.5, 299.5);
  fOutputListQA->Add(fSPDnTrAcc);

  fSPDnTrCorrMax = new TH2F("fSPDnTrCorrMax", "fSPDnTrCorrMax; zSPDVtx; SPDnTr*PoisCorr(max) (eta<1)", 220, -11,11,300, -0.5, 299.5);
  fOutputListQA->Add(fSPDnTrCorrMax);

  Int_t    nBinsSPD[3]={10, 200, 200};
  Double_t xminSPD[3]={-10,-0.5,-0.5};
  Double_t xmaxSPD[3]={10,199.5, 199.5};

  if (fIsMC) {
    fSPDnTrGen = new TH2F("fSPDnTrGen", "fSPDnTrGen; zSPDVtx; MC_ChargedParticles (eta<1)", 220, -11, 11, 300, -0.5, 299.5);
    fOutputListQA->Add(fSPDnTrGen);

    
    fDiffSPDMCVtx = new TH2F("fDiffSPDMCVtx", "fDiffSPDMCVtx; NContributors; zSPDVtx-zMCVtx", 200, -0.5, 199.5, 300, -3., 3.);
    fOutputListQA->Add(fDiffSPDMCVtx);
     
    fCheckMCVertex = new TH2F("fCheckMCVertex", "TruthVsRec/track Vertex; MC Vtx; Rec Vtx", 200, -10, 10, 200, -10, 10);
    fOutputListQA->Add(fCheckMCVertex);
  
    
    fnTrAccMaxGen = new THnSparseF("fnTrAccMaxGen", "fnTrAccMaxGen; zSPDVtx; MC_ChargedParticles (eta<1);SPDnTR*PoisCorr(max) (eta<1);", 3, nBinsSPD, xminSPD, xmaxSPD);
    fOutputListQA->Add(fnTrAccMaxGen);

    fnTrAccGen = new THnSparseF("fnTrAccGen", "fnTrAccGen; zSPDVtx; MC_ChargedParticles (eta<1); MC_SPDnTr reconstructed (eta<1);", 3, nBinsSPD, xminSPD, xmaxSPD);
    fOutputListQA->Add(fnTrAccGen);

    fnTrAccGenTrueInel = new TH2F("fnTrAccGenTrueInel", "fnTrAccGenTrueInel;  MC_SPDnTr reconstructed (eta<1); MC_ChargedParticles(eta<1)", nBinsSPD[1], xminSPD[1], xmaxSPD[1], nBinsSPD[2], xminSPD[2], xmaxSPD[2]);
    fOutputListQA->Add(fnTrAccGenTrueInel);

    // fnTrAccGenTrueInelTrig = new TH2F("fnTrAccGenTrueInelTrig", "fnTrAccGenTrueInel;  MC_SPDnTr reconstructed (eta<1); MC_ChargedParticles(eta<1)", nBinsSPD[1], xminSPD[1], xmaxSPD[1], nBinsSPD[2], xminSPD[2], xmaxSPD[2]);
    // fOutputListQA->Add(fnTrAccGenTrueInelTrig);

    // fnTrAccGenTrueInelVtxQA = new TH2F("fnTrAccGenTrueInelVtxQA", "fnTrAccGenTrueInel;  MC_SPDnTr reconstructed (eta<1); MC_ChargedParticles(eta<1)", nBinsSPD[1], xminSPD[1], xmaxSPD[1], nBinsSPD[2], xminSPD[2], xmaxSPD[2]);
    // fOutputListQA->Add(fnTrAccGenTrueInelVtxQA);

    // fnTrAccGenTrueInelVtxEx = new TH2F("fnTrAccGenTrueInelVtxEx", "fnTrAccGenTrueInel;  MC_SPDnTr reconstructed (eta<1); MC_ChargedParticles(eta<1)", nBinsSPD[1], xminSPD[1], xmaxSPD[1], nBinsSPD[2], xminSPD[2], xmaxSPD[2]);
    // fOutputListQA->Add(fnTrAccGenTrueInelVtxEx);
    
    // fnTrAccMinGen = new THnSparseF("fnTrAccMinGen", "fnTrAccMinGen", 3, nBinsSPD, xminSPD, xmaxSPD);
    //  fOutputListQA->Add(fnTrAccMinGen);
  
    //  fnTrAccMeanGen = new THnSparseF("fnTrAccMenGen", "fnTrAccMeanGen",  3, nBinsSPD, xminSPD, xmaxSPD);
    //  fOutputListQA->Add(fnTrAccMeanGen);


    fPt2Tagged = new TH3F("fPt2Tagged", "Pt2Tagged; RecPtE; MCPt2; RecLP", 20, 0.5,3,60, 0.,3., 20,0.,10.);
    fOutputListQA->Add(fPt2Tagged);

    fPt2NTagged = new TH3F("fPt2NTagged", "Pt2NTagged; RecPtE; MCPt2; RecLP", 20, 0.5,3,60, 0.,3., 20,0.,10.);
    fOutputListQA->Add(fPt2NTagged);

    // fMCThrustTagged = new TH3F("fMCThrustTagged", "Thrust of Tagged; ptElectron; Thrust", 44, 0.5, 6., 100, 0., 1., 8, -0.5, 7.5);
    //fOutputListQA->Add(fMCThrustTagged);

    fMCSpherTagged = new TH3F("fMCSpherTagged", "Sphericity of Tagged; ptElectron; Thrust", 44, 0.5, 6., 100, 0., 1., 8, -0.5, 7.5);
    fOutputListQA->Add(fMCSpherTagged);
    
    fRecLPTagged = new TH3F("fRecLPTagged", "LP of Tagged; ptElectron; ptLP", 44, 0.5, 6., 200, 0, 20, 8, -0.5, 7.5);
    fOutputListQA->Add(fRecLPTagged);

    // fMultCorrTagged = new TH3F("fMultCorrTagged", "Corr SPD mult of tagged; ptElectron; nTr", 44, 0.5, 6., 150, -0.5, 149.5, 8, -0.5, 7.5);
    // fOutputListQA->Add(fMultCorrTagged);

    fNHadTagged = new TH3F("fNHadTagged", "NHadrons for tagged; ptElectron; nHad", 44, 0.5, 6., 150, -0.5, 149.5, 8, -0.5, 7.5);
    fOutputListQA->Add(fNHadTagged);

    //    fNHadTaggedA = new TH3F("fNHadTaggedA", "NHadrons for tagged; ptElectron; nHad", 44, 0.5, 6., 150, -0.5, 149.5, 8, -0.5, 7.5);
    //fOutputListQA->Add(fNHadTaggedA);

    fNHadTaggedB = new TH3F("fNHadTaggedB", "NHadrons for tagged; ptElectron; nHad", 44, 0.5, 6., 150, -0.5, 149.5, 8, -0.5, 7.5);
    fOutputListQA->Add(fNHadTaggedB);

    //fNHadTaggedC = new TH3F("fNHadTaggedC", "NHadrons for tagged; ptElectron; nHad", 44, 0.5, 6., 150, -0.5, 149.5, 8, -0.5, 7.5);
    //fOutputListQA->Add(fNHadTaggedC);

    fMeanPtTagged = new TH3F("fMeanPtTagged", "Mean pt for tagged; ptElectron; average pt", 44, 0.5, 6., 100, 0, 10, 8, -0.5, 7.5);
    fOutputListQA->Add(fMeanPtTagged);

    //   fMeanPtTaggedA = new TH3F("fMeanPtTaggedA", "Mean pt for tagged; ptElectron; average pt", 44, 0.5, 6., 100, 0, 10, 8, -0.5, 7.5);
    //   fOutputListQA->Add(fMeanPtTaggedA);

    fMeanPtTaggedB = new TH3F("fMeanPtTaggedB", "Mean pt for tagged; ptElectron; average pt; case", 44, 0.5, 6., 100, 0, 10, 8, -0.5, 7.5);
    fOutputListQA->Add(fMeanPtTaggedB);

    //fMeanPtTaggedC = new TH3F("fMeanPtTaggedC", "Mean pt for tagged; ptElectron; average pt", 44, 0.5, 6., 100, 0, 10, 8, -0.5, 7.5);
    // fOutputListQA->Add(fMeanPtTaggedC);

    // fMCThrustNTagged = new TH3F("fMCThrustNTagged", "Thrust of NTagged; ptElectron; Thrust", 44, 0.5, 6., 100, 0., 1., 8, -0.5, 7.5);
    // fOutputListQA->Add(fMCThrustNTagged);

    fMCSpherNTagged = new TH3F("fMCSpherNTagged", "Sphericity of NTagged; ptElectron; Thrust; case", 44, 0.5, 6., 100, 0., 1., 8, -0.5, 7.5);
    fOutputListQA->Add(fMCSpherNTagged);

    fRecLPNTagged = new TH3F("fRecLPNTagged", "LP of NTagged; ptElectron; ptLP; case", 44, 0.5, 6., 200, 0, 20, 8, -0.5, 7.5);
    fOutputListQA->Add(fRecLPNTagged);

    // fMultCorrNTagged = new TH3F("fMultCorrNTagged", "Corr SPD mult of tagged; ptElectron; nTr", 44, 0.5, 6., 150, -0.5, 149.5, 8, -0.5, 7.5);
    //fOutputListQA->Add(fMultCorrNTagged);

    fNHadNTagged = new TH3F("fNHadNTagged", "NHadrons for tagged; ptElectron; nHad; case", 44, 0.5, 6., 150, -0.5, 149.5, 8, -0.5, 7.5);
    fOutputListQA->Add(fNHadNTagged);

    //  fNHadNTaggedA = new TH3F("fNHadNTaggedA", "NHadrons for tagged; ptElectron; nHad", 44, 0.5, 6., 150, -0.5, 149.5, 8, -0.5, 7.5);
    // fOutputListQA->Add(fNHadNTaggedA);

    fNHadNTaggedB = new TH3F("fNHadNTaggedB", "NHadrons for tagged; ptElectron; nHad; case", 44, 0.5, 6., 150, -0.5, 149.5, 8, -0.5, 7.5);
    fOutputListQA->Add(fNHadNTaggedB);

    //fNHadNTaggedC = new TH3F("fNHadNTaggedC", "NHadrons for tagged; ptElectron; nHad", 44, 0.5, 6., 150, -0.5, 149.5, 8, -0.5, 7.5);
    //fOutputListQA->Add(fNHadNTaggedC);

    fMeanPtNTagged = new TH3F("fMeanPtNTagged", "Mean pt for tagged; ptElectron; average pt; case", 44, 0.5, 6., 100, 0, 10, 8, -0.5, 7.5);
    fOutputListQA->Add(fMeanPtNTagged);

    // fMeanPtNTaggedA = new TH3F("fMeanPtNTaggedA", "Mean pt for tagged; ptElectron; average pt", 44, 0.5, 6., 100, 0, 10, 8, -0.5, 7.5);
    //   fOutputListQA->Add(fMeanPtNTaggedA);

    fMeanPtNTaggedB = new TH3F("fMeanPtNTaggedB", "Mean pt for tagged; ptElectron; average pt; case", 44, 0.5, 6., 100, 0, 10, 8, -0.5, 7.5);
    fOutputListQA->Add(fMeanPtNTaggedB);

    //fMeanPtNTaggedC = new TH3F("fMeanPtNTaggedC", "Mean pt for tagged; ptElectron; average pt", 44, 0.5, 6., 100, 0, 10, 8, -0.5, 7.5);
    //fOutputListQA->Add(fMeanPtNTaggedC);



    // fMothMCThrustTagged = new TH2F("fMothMCThrustTagged", "Thrust of Tagged; ptMother; Thrust", 44, 0.5, 6., 100, 0., 1.);
    // fOutputListQA->Add(fMothMCThrustTagged);

    fMothMCSpherTagged = new TH2F("fMothMCSpherTagged", "Sphericity of Tagged; ptMother; Thrust", 44, 0.5, 6., 100, 0., 1.);
     fOutputListQA->Add(fMothMCSpherTagged);
    
    fMothRecLPTagged = new TH2F("fMothRecLPTagged", "LP of Tagged; ptMother; ptLP", 44, 0.5, 6., 200, 0, 20);
    fOutputListQA->Add(fMothRecLPTagged);

    //  fMothMultCorrTagged = new TH2F("fMothMultCorrTagged", "Corr SPD mult of tagged; ptMother; nTr", 44, 0.5, 6., 150, -0.5, 149.5);
    // fOutputListQA->Add(fMothMultCorrTagged);

    fMothNHadTagged = new TH2F("fMothNHadTagged", "NHadrons for tagged; ptMother; nHad", 44, 0.5, 6., 150, -0.5, 149.5);
    fOutputListQA->Add(fMothNHadTagged);

    fMothMeanPtTagged = new TH2F("fMothMeanPtTagged", "Mean pt for tagged; ptMother; average pt", 44, 0.5, 6., 100, 0, 10);
    fOutputListQA->Add(fMothMeanPtTagged);

    
    //fMothMCThrustNTagged = new TH2F("fMothMCThrustNTagged", "Thrust of NTagged; ptMother; Thrust", 44, 0.5, 6., 100, 0., 1.);
    // fOutputListQA->Add(fMothMCThrustNTagged);

    fMothMCSpherNTagged = new TH2F("fMothMCSpherNTagged", "Sphericity of NTagged; ptMother; Thrust", 44, 0.5, 6., 100, 0., 1.);
    fOutputListQA->Add(fMothMCSpherNTagged);

    fMothRecLPNTagged = new TH2F("fMothRecLPNTagged", "LP of NTagged; ptMother; ptLP", 44, 0.5, 6., 200, 0, 20);
    fOutputListQA->Add(fMothRecLPNTagged);

    // fMothMultCorrNTagged = new TH2F("fMothMultCorrNTagged", "Corr SPD mult of tagged; ptMother; nTr", 44, 0.5, 6., 150, -0.5, 149.5);
    // fOutputListQA->Add(fMothMultCorrNTagged);

    fMothNHadNTagged = new TH2F("fMothNHadNTagged", "NHadrons of tagged; ptMother; nHad", 44, 0.5, 6., 150, -0.5, 149.5);
    fOutputListQA->Add(fMothNHadNTagged);

    fMothMeanPtNTagged = new TH2F("fMothMeanPtNTagged", "Mean pt for tagged; ptMother; average pt", 44, 0.5, 6., 100, 0, 10);
    fOutputListQA->Add(fMothMeanPtNTagged);



    //fMCThrustTaggedH = new TH2F("fMCThrustTaggedH", "Thrust of TaggedH; ptElectron; Thrust", 44, 0.5, 6., 100, 0., 1.);
    //fOutputListQA->Add(fMCThrustTaggedH);

    fMCSpherTaggedH = new TH2F("fMCSpherTaggedH", "Sphericity of TaggedH; ptElectron; Thrust", 44, 0.5, 6., 100, 0., 1.);
    fOutputListQA->Add(fMCSpherTaggedH);
    
    fRecLPTaggedH = new TH2F("fRecLPTaggedH", "LP of TaggedH; ptElectron; ptLP", 44, 0.5, 6., 200, 0, 20);
    fOutputListQA->Add(fRecLPTaggedH);

    //  fMultCorrTaggedH = new TH2F("fMultCorrTaggedH", "Corr SPD mult of tagged; ptElectron; nTr", 44, 0.5, 6., 150, -0.5, 149.5);
    // fOutputListQA->Add(fMultCorrTaggedH);

    fNHadTaggedH = new TH2F("fNHadTaggedH", "NHadrons for tagged; ptElectron; nHad", 44, 0.5, 6., 150, -0.5, 149.5);
    fOutputListQA->Add(fNHadTaggedH);

    fMeanPtTaggedH = new TH2F("fMeanPtTaggedH", "Mean pt for tagged; ptElectron; average pt", 44, 0.5, 6., 100, 0, 10);
    fOutputListQA->Add(fMeanPtTaggedH);


    // fMCThrustNTaggedH = new TH2F("fMCThrustNTaggedH", "Thrust of NTaggedH; ptElectron; Thrust", 44, 0.5, 6., 100, 0., 1.);
    // fOutputListQA->Add(fMCThrustNTaggedH);

    fMCSpherNTaggedH = new TH2F("fMCSpherNTaggedH", "Sphericity of NTaggedH; ptElectron; Thrust", 44, 0.5, 6., 100, 0., 1.);
    fOutputListQA->Add(fMCSpherNTaggedH);

    fRecLPNTaggedH = new TH2F("fRecLPNTaggedH", "LP of NTaggedH; ptElectron; ptLP", 44, 0.5, 6., 200, 0, 20);
    fOutputListQA->Add(fRecLPNTaggedH);

    // fMultCorrNTaggedH = new TH2F("fMultCorrNTaggedH", "Corr SPD mult of tagged; ptElectron; nTr", 44, 0.5, 6., 150, -0.5, 149.5);
    // fOutputListQA->Add(fMultCorrNTaggedH);

    fNHadNTaggedH = new TH2F("fNHadNTaggedH", "NHadrons for tagged; ptElectron; nHad", 44, 0.5, 6., 150, -0.5, 149.5);
    fOutputListQA->Add(fNHadNTaggedH);

    fMeanPtNTaggedH = new TH2F("fMeanPtNTaggedH", "Mean pt for tagged; ptElectron; average pt", 44, 0.5, 6., 100, 0, 10);
    fOutputListQA->Add(fMeanPtNTaggedH);


    //fMothMCThrustTaggedH = new TH2F("fMothMCThrustTaggedH", "Thrust of TaggedH; ptMother; Thrust", 44, 0.5, 6., 100, 0., 1.);
    //fOutputListQA->Add(fMothMCThrustTaggedH);

    fMothMCSpherTaggedH = new TH2F("fMothMCSpherTaggedH", "Sphericity of TaggedH; ptMother; Thrust", 44, 0.5, 6., 100, 0., 1.);
    fOutputListQA->Add(fMothMCSpherTaggedH);
    
    fMothRecLPTaggedH = new TH2F("fMothRecLPTaggedH", "LP of TaggedH; ptMother; ptLP", 44, 0.5, 6., 200, 0, 20);
    fOutputListQA->Add(fMothRecLPTaggedH);

    //fMothMultCorrTaggedH = new TH2F("fMothMultCorrTaggedH", "Corr SPD mult of tagged; ptMother; nTr", 44, 0.5, 6., 150, -0.5, 149.5);
    //fOutputListQA->Add(fMothMultCorrTaggedH);

    fMothNHadTaggedH = new TH2F("fMothNHadTaggedH", "NHadrons for tagged; ptMother; nHad", 44, 0.5, 6., 150, -0.5, 149.5);
    fOutputListQA->Add(fMothNHadTaggedH);

    fMothMeanPtTaggedH = new TH2F("fMothMeanPtTaggedH", "Mean pt for tagged; ptMother; average pt", 44, 0.5, 6., 100, 0, 10);
    fOutputListQA->Add(fMothMeanPtTaggedH);


    //fMothMCThrustNTaggedH = new TH2F("fMothMCThrustNTaggedH", "Thrust of NTaggedH; ptMother; Thrust", 44, 0.5, 6., 100, 0., 1.);
    //fOutputListQA->Add(fMothMCThrustNTaggedH);

    fMothMCSpherNTaggedH = new TH2F("fMothMCSpherNTaggedH", "Sphericity of NTaggedH; ptMother; Thrust", 44, 0.5, 6., 100, 0., 1.);
    fOutputListQA->Add(fMothMCSpherNTaggedH);

    fMothRecLPNTaggedH = new TH2F("fMothRecLPNTaggedH", "LP of NTaggedH; ptMother; ptLP", 44, 0.5, 6., 200, 0, 20);
    fOutputListQA->Add(fMothRecLPNTaggedH);

    // fMothMultCorrNTaggedH = new TH2F("fMothMultCorrNTaggedH", "Corr SPD mult of tagged; ptMother; nTr", 44, 0.5, 6., 150, -0.5, 149.5);
    //fOutputListQA->Add(fMothMultCorrNTaggedH);

    fMothNHadNTaggedH = new TH2F("fMothNHadNTaggedH", "NHadrons of tagged; ptMother; nHad", 44, 0.5, 6., 150, -0.5, 149.5);
    fOutputListQA->Add(fMothNHadNTaggedH);

    fMothMeanPtNTaggedH = new TH2F("fMothMeanPtNTaggedH", "Mean pt for tagged; ptMother; average pt", 44, 0.5, 6., 100, 0, 10);
    fOutputListQA->Add(fMothMeanPtNTaggedH);

  }

  fnTrAccMax = new THnSparseF("fnTrAccMax", "fnTrAccMax; SPDzVtx; SPDnTr; SPDnTr+Pois(max)", 3, nBinsSPD, xminSPD, xmaxSPD);
  fOutputListQA->Add(fnTrAccMax);

  //fnTrAccMin = new THnSparseF("fnTrAccMin", "fnTrAccMin",  3, nBinsSPD, xminSPD, xmaxSPD);
  // fOutputListQA->Add(fnTrAccMin);

  // fnTrAccMean = new THnSparseF("fnTrAccMean", "fnTrAccMean",  3, nBinsSPD, xminSPD, xmaxSPD);
  // fOutputListQA->Add(fnTrAccMean);
  /*
  Int_t    nBinsMult[4]={20, 20, 50, 50};
  Double_t xminMult[4]={0,0,0,0};
  Double_t xmaxMult[4]={100,100, 3000, 300};
  */
  //  fMultiplicity = new THnSparseF("fMultiplicity", "Multiplicity: VOPerc, SPDPerc, TotV0, SPD", 4, nBinsMult, xminMult, xmaxMult);
  //fOutputListQA->Add(fMultiplicity);
  
  //  fSPDMultiplicity = new TH3F("fSPDMultiplicity", "Multiplicity: SPDPerc, SPDTracklets, Zvtx", 20, 0., 100.001, 301, -0.5, 300.5, 100, -10., 10);
  // fOutputListQA->Add(fSPDMultiplicity);


  // Track Cuts
  const char*  ElecTrackCutsLabels[15]= {"All", "FilterBit", "PtCut", "EtaCut", "TPCNcls", "TPCNclsdEdx", "TPCrefit", "TPCfrac", "KinkCut", "ITSNcls", "ITSrefit", "SPDAny", "SPDBoth", "DCACut", "ITSSCluster"};

 fElectronTrackCuts = new TH2F("fElectronTrackCuts", "fElectronTrackCuts", 20, 0, 10, 15, 0,15);
  for (Int_t i=0; i<15; i++) fElectronTrackCuts->GetYaxis()->SetBinLabel(i+1, ElecTrackCutsLabels[i]);
  fOutputListQA->Add(fElectronTrackCuts);



  
  fElectronTrackTPCChi2 = new TH2F("fElectronTrackTPCChi2", "fElectronTrackTPCChi2", 20, 0, 10,50, 0, 10);
  fOutputListQA->Add(fElectronTrackTPCChi2);
  
  fElectronTrackTPCCrossedRows = new TH2F("fElectronTrackTPCCrossedRows", "fElectronTrackTPCCrossedRows", 20, 0, 10,170, 0, 170);
  fOutputListQA->Add(fElectronTrackTPCCrossedRows);
  
  fElectronTrackTPCNcls = new TH2F("fElectronTrackTPCNcls", "fElectronTrackTPCNcls", 20, 0, 10,170, 0, 170);
  fOutputListQA->Add(fElectronTrackTPCNcls);

  fElectronTrackTPCNclsdEdx= new TH2F("fElectronTrackTPCNclsdEdx", "fElectronTrackTPCNclsdEdx",20, 0, 10, 170, 0, 170);
  fOutputListQA->Add(fElectronTrackTPCNclsdEdx);

  fElectronTrackTPCFrac = new TH2F("fElectronTrackTPCFrac", "fElectronTrackTPCFrac", 20, 0, 10,110, 0, 1.1);
  fOutputListQA->Add(fElectronTrackTPCFrac);
  
  fElectronTrackITSNcls = new TH2F("fElectronTrackITSNcls", "fElectronTrackITSNcls", 20, 0, 10,10, 0, 10);
  fOutputListQA->Add(fElectronTrackITSNcls);
  
  fElectronTrackITSChi2 = new TH2F("fElectronTrackITSChi2", "fElectronTrackITSChi2", 20, 0, 10,100, 0, 40);
  fOutputListQA->Add(fElectronTrackITSChi2);

  fElectronTrackITSLayer = new TH3F("fElectronTrackITSLayer", "fElectronTrackITSLayer", 20, 0, 10,2, -0.5, 1.5, 2, -0.5, 1.5);
  fOutputListQA->Add(fElectronTrackITSLayer);
 
  fElectronTrackRefit = new TH3F("fElectronTrackRefit", "fElectronTrackReft", 20, 0, 10,2, -0.5, 1.5, 2, -0.5, 1.5);
  fOutputListQA->Add(fElectronTrackRefit);

  fElectronTrackDCA = new TH3F("fElectronTrackDCA", "fElectronTrackDCA: pt, r, z", 10, 0, 5, 60, -3.0, 3.0, 60, -3.0, 3.0);
  fOutputListQA->Add(fElectronTrackDCA);

 
  Int_t BinsITSCuts[5]={9, 7, 7, 20, 8};
  Double_t BinsITSMin[5]={0.0,-0.5,-0.5,0.,0.};
  Double_t BinsITSMax[5]={4.5, 6.5,6.5,20.,399.};

  fElectronTrackITSCuts = new THnSparseF("fElectronTrackITSCuts", "fElectronTrackITSCuts; Pt; NCls; SharedCls; X2/NCls; ParticleID;", 5, BinsITSCuts, BinsITSMin, BinsITSMax);
  Double_t ParticleIDBins[9]={0.,4.5, 5.5,20., 110.,220.,222.,224.,399.};
  fElectronTrackITSCuts->GetAxis(4)->Set(BinsITSCuts[4], ParticleIDBins);
  fOutputListQA->Add(fElectronTrackITSCuts);

  /*
  Int_t BinsPITSCuts[5]={30, 7, 7, 20, 2};
  Double_t BinsPITSMin[5]={0.0,-0.5,-0.5,0.,-0.5};
  Double_t BinsPITSMax[5]={1.5, 6.5,6.5,20.,1.5};

  fPhotTrackITSCuts = new THnSparseF("fPhotTrackITSCuts", "fPhotTrackITSCuts; Pt; NCls; SharedCls; X2/NCls; IsPhotonic", 5, BinsPITSCuts, BinsPITSMin, BinsPITSMax);
  fOutputListQA->Add(fPhotTrackITSCuts);
  */
 

  const char*  HadTrackCutsLabels[11]= {"All", "FilterBit", "PtCut", "EtaCut",  "KinkCut", "TPCNcls",  "TPCrefit", "ITSrefit", "DCACut", "SPDAny", "TOFmatch"};
  fHadronTrackCuts = new TH2F("fHadronTrackCuts", "fHadronTrackCuts", 20, 0, 10, 11, 0,11);
  for (Int_t i=0; i<11; i++) fHadronTrackCuts->GetYaxis()->SetBinLabel(i+1, HadTrackCutsLabels[i]);
  fOutputListQA->Add(fHadronTrackCuts);

  fHadronTrackTPCNcls = new TH2F("fHadronTrackTPCNcls", "fHadronTrackTPCNcls", 20, 0, 10,170, 0, 170);
  fOutputListQA->Add(fHadronTrackTPCNcls);

  fHadronTrackRefit = new TH3F("fHadronTrackRefit", "fHadronTrackReft", 20, 0, 10,2, -0.5, 1.5, 2, -0.5, 1.5);
  fOutputListQA->Add(fHadronTrackRefit);

  fHadronTrackDCA = new TH3F("fHadronTrackDCA", "fHadronTrackDCA: pt, r, z", 20, 0, 10,  60, -3.0, 3.0, 60, -3.0, 3.0);
  fOutputListQA->Add(fHadronTrackDCA);

  
  fHadronTrackDCA_woITSAny = new TH3F("fHadronTrackDCA_woITSAny", "fHadronTrackDCA_woITSAny: pt, r, z", 20, 0, 10,  60, -3.0, 3.0, 60, -3.0, 3.0);
  fOutputListQA->Add(fHadronTrackDCA_woITSAny);

  fHadronTrackDCA_wITSAny = new TH3F("fHadronTrackDCA_wITSAny", "fHadronTrackDCA_wITSAny: pt, r, z", 20, 0, 10,  60, -3.0, 3.0, 60, -3.0, 3.0);
  fOutputListQA->Add(fHadronTrackDCA_wITSAny);


  fHistITSnSig = new TH2F("fHistITSnSig","fHistITSnSig",50,0.3,5,200,-10,10);
  BinLogX(fHistITSnSig->GetXaxis());
  fOutputListQA->Add(fHistITSnSig);
  
  fHistTOFnSig = new TH2F("fHistTOFnSig","fHistTOFnSig",100,0.3,10,200,-10,10);
  BinLogX(fHistTOFnSig->GetXaxis());
  fOutputListQA->Add(fHistTOFnSig);
  
  fHistTPCnSig = new TH2F("fHistTPCnSig","fHistTPCnSig",150,0.3,15,200,-10,10);
  BinLogX(fHistTPCnSig->GetXaxis());
  fOutputListQA->Add(fHistTPCnSig);
    
  fHistTPCnSigITScut = new TH2F("fHistTPCnSigITScut","fHistTPCnSigITScut",150,0.3,15,200,-10,10);
  BinLogX(fHistTPCnSigITScut->GetXaxis());
  fOutputListQA->Add(fHistTPCnSigITScut);
    
  fHistTPCnSigTOFcut = new TH2F("fHistTPCnSigTOFcut","fHistTPCnSigTOFcut",150,0.3,15,200,-10,10);
  BinLogX(fHistTPCnSigTOFcut->GetXaxis());
  fOutputListQA->Add(fHistTPCnSigTOFcut);
    
  fHistTPCnSigITSTOFcut = new TH2F("fHistTPCnSigITSTOFcut","fHistTPCnSigITSTOFcut",150,0.3,15,200,-10,10);
  BinLogX(fHistTPCnSigITSTOFcut->GetXaxis());
  fOutputListQA->Add(fHistTPCnSigITSTOFcut);

  fHistITSnSigTOFTPCcut = new TH2F("fHistITSnSigTOFTPCcut","fHistITSnSigTOFTPCcut",150,0.3,15,200,-10,10);
  BinLogX(fHistITSnSigTOFTPCcut->GetXaxis());
  fOutputListQA->Add(fHistITSnSigTOFTPCcut);

  

  // 
  Int_t    binHadScaling[5]={25,25, 25, 25, 10};
  Double_t xminHadScaling[5]={-0.5,-0.5,-0.5,-0.5, 0};
  Double_t xmaxHadScaling[5]={99.5,49.5, 49.5, 49.5, 100.01};

  //  fCheckNHadronScaling = new THnSparseF("fCheckNHadronScaling", "NHadScaling: NHadron, NElectron, NNonElectron, HWrongElectron, Mult", 5, binHadScaling, xminHadScaling, xmaxHadScaling);
  //  fOutputList->Add(fCheckNHadronScaling);

  //  fCheckNPhotHadScaling = new THnSparseF("fCheckNPhotHadScaling", "NHadScaling: NHadron, NElectron, NTagged, NNotTagged, Mult", 5, binHadScaling, xminHadScaling, xmaxHadScaling);
  // fOutputList->Add(fCheckNPhotHadScaling);

  // fCheckTaggedEvent = new TH3F("fChecTaggedEvent", "CheckTaggedEvent: Case (0-both, 1tagged, 2nontagged), LPpt, Mult for ptphot>1", 3, -0.5, 2.5, 40, 0., 20., 20, -0.5, 99.5);
  // fOutputList->Add(fCheckTaggedEvent);

  if (fHadCont) {

    Int_t    binHC2[4] =  {NBinsElectron  ,10   , 10 ,100}; //p, ITS, TOF, TPC
    Double_t xminHC2[4] = {XminElectron   ,-5  ,-5 ,-10};
    Double_t xmaxHC2[4] = {XmaxElectron   ,5   , 5 ,10};  
    fHadContamination = new THnSparseF("fHadContamination", "HadCont; P; ITS; TOF; TPC;", 4, binHC2, xminHC2, xmaxHC2);
    fOutputListMain->Add(fHadContamination);

    
    if (fOneTimeCheck) {
      fHadContPvsPt = new TH2F("fHadContPvsPt", "PvsPt; P;Pt", 100, 0, 10, 100, 0, 10);
      fOutputListMain->Add(fHadContPvsPt);  
    
      Int_t    binHC[4] =  {NBinsElectron  ,NBinsPhi/8      ,NBinsEta/2    ,200}; //p, Phi, Eta, TPC
      Double_t xminHC[4] = {XminElectron   ,0             ,-0.9  ,-10};
      Double_t xmaxHC[4] = {XmaxElectron   ,TMath::TwoPi(), 0.9  ,10};
      fHadContPPhiEtaTPC = new THnSparseF("fHadContPPhiEtaTPC", "HadCont; P; Phi; Eta; TPC;", 4, binHC, xminHC, xmaxHC);
      fOutputListMain->Add(fHadContPPhiEtaTPC); // select electrons with had cont

      fHadContTPCEtaPhiPt = new TH3F("fHadContTPCEtaPhiPt", "MCHadContTPC; Eta; Phi; Pt", 36, -0.9, 0.9, 32, 0, TMath::TwoPi(), NBinsElectron, XminElectron, XmaxElectron);
      fOutputListMain->Add(fHadContTPCEtaPhiPt); // select hadrons by intend
    }
 
   
    if (fIsMC) {
      if(fOneTimeCheck) fHadContEtaPhiPt = new TH3F("fHadContEtaPhiPt", "MCHadCont; Eta; Phi; Pt", 36, -0.9, 0.9, 32, 0, TMath::TwoPi(), NBinsElectron, XminElectron, XmaxElectron);
      fOutputListMain->Add(fHadContEtaPhiPt);   
      
      Int_t    binHMC[4] = {NBinsElectron  , 7,  10, 100}; // p PDG ITS, TPC
      Double_t xminHMC[4] = {XminElectron   , 0, -5, -10};
      Double_t xmaxHMC[4] = {XmaxElectron   , 7,  5,  10};
      fHadContMC = new THnSparseF("fHadContMC", "HadContMC; P; PDG; ITS; TPC;", 4, binHMC, xminHMC, xmaxHMC);
      fOutputListMain->Add(fHadContMC); // had cont after tof cut
    }
  }



  fTrkpt = new TH2F("fTrkpt","track pt (0),after EleTrack (1), after ElePID(2); pt ",200,0,20,3,-0.5,2.5);
  fOutputListMain->Add(fTrkpt);
  
  fInclElecPtEta = new TH3F("fInclElePtEta", "fInclElePtEta; #it{p}_{T}; #eta",NBinsElectron,XminElectron, XmaxElectron, NBinsEta, -0.9, 0.9, NMultBins, 0, NMultBins);
  fInclElecPtEta->GetZaxis()->Set(NMultBins,XMultBins);
  fOutputListMain->Add(fInclElecPtEta);

  fInclElecPtEtaWRecEff = new TH3F("fInclElePtEtaWRecEff", "fInclElePtEtaWRecEff; #it{p}_{T}; #eta",NBinsElectron,XminElectron, XmaxElectron,  NBinsEta, -0.9, 0.9 , NMultBins, 0, NMultBins);
  fInclElecPtEtaWRecEff->GetZaxis()->Set(NMultBins,XMultBins);
  fOutputListMain->Add(fInclElecPtEtaWRecEff);

  if (fOneTimeCheck) {
    fInclElecPhi = new TH2F("fInclElePhi", "fInclElePhi",NBinsElectron, XminElectron, XmaxElectron,NBinsPhi,0,TMath::TwoPi());
    fOutputListMain->Add(fInclElecPhi);
    
    fInclElecPhiEta = new TH2F("fInclElePhiEta", "fInclElePhiEta", NBinsEta,XminEta,XmaxEta, NBinsPhi,0,TMath::TwoPi());
    fOutputListMain->Add(fInclElecPhiEta);
  }

  fULSElecPt = new TH2F("fULSElePt", "fULSElePt; #it{p}_{T}",NBinsElectron,XminElectron, XmaxElectron, NMultBins, 0, NMultBins);
  fULSElecPt->GetYaxis()->Set(NMultBins,XMultBins);
  fOutputListMain->Add(fULSElecPt);
    
  fULSElecPtWRecEff = new TH2F("fULSElePtWRecEff", "fULSElePtWRecEff; #it{p}_{T}",NBinsElectron,XminElectron, XmaxElectron, NMultBins, 0, NMultBins);
  fULSElecPtWRecEff->GetYaxis()->Set(NMultBins,XMultBins);
  fOutputListMain->Add(fULSElecPtWRecEff);
    
  fLSElecPt = new TH2F("fLSElePt", "fLSElePt; #it{p}_{T}",NBinsElectron,XminElectron, XmaxElectron, NMultBins, 0, NMultBins);
  fLSElecPt->GetYaxis()->Set(NMultBins,XMultBins);
  fOutputListMain->Add(fLSElecPt);

  fLSElecPtWRecEff = new TH2F("fLSElePtWRecEff", "fLSElePtWRecEff; #it{p}_{T}",NBinsElectron,XminElectron, XmaxElectron, NMultBins, 0, NMultBins);
  fLSElecPtWRecEff->GetYaxis()->Set(NMultBins,XMultBins);
  fOutputListMain->Add(fLSElecPtWRecEff);
    
  if (fOneTimeCheck) {
    fULSElecPhi= new TH2F("fULSElePhi", "fULSElePhi", NBinsElectron, XminElectron, XmaxElectron ,NBinsPhi,0,TMath::TwoPi());
    fOutputListMain->Add(fULSElecPhi);
    
    fLSElecPhi= new TH2F("fLSElePhi", "fLSElePhi",NBinsElectron, XminElectron, XmaxElectron,NBinsPhi,0,TMath::TwoPi());
    fOutputListMain->Add(fLSElecPhi);
  }

  fInvmassULS = new TH2F("fInvmassULS", "fInvmassULS; m_{ee}; pt", 500,0,0.5,NBinsElectronRed, XBinsElectronRed);
  fOutputListMain->Add(fInvmassULS);
    
  fInvmassLS = new TH2F("fInvmassLS", "fInvmassLS; m_{ee}; pt", 500,0,0.5,NBinsElectronRed, XBinsElectronRed);
  fOutputListMain->Add(fInvmassLS);


  fPhotPt1Pt2Only = new TH2F("fPhotPt1Pt2Only", "fPhotPt1Pt2 for tagged; Pt1Rec; Pt2Rec", 48, 0., 6., 48, 0., 6.);
  fOutputListMain->Add(fPhotPt1Pt2Only);
    
  if (fOneTimeCheck) {
    fPhotMixULS = new TH1F("fPhotMixULS", "fPhotMixULS; #it{p}_{T}", NBinsElectron,XminElectron, XmaxElectron);
    fOutputListMain->Add(fPhotMixULS);
    
    fPhotMixLS = new TH2F("fPhotMixLS", "fPhotMixcLS; #it{p}_{T}; Charge", NBinsElectron,XminElectron, XmaxElectron, 3, -1.5, 1.5);
    fOutputListMain->Add(fPhotMixLS);
    
 
  }

  
  if (fIsMC) {

    //    fRecMCInvMass = new TH3F("fRecMCInvMass", "fRecMCInvMass; MCPt; rec. inv. mass; MC inv. mass", 6, 0., 3., 14, 0, 0.14, 18, 0, 0.18);
    // fOutputListMain->Add(fRecMCInvMass);

    fInvmassMCTrue = new TH3F("fInvmassMCTrue", "fInvmassMCTrue; mee; pt; pdg", 100,0,0.5,24, 0, 6, PDGSize, 0.5, PDGSize+0.5);
    SetPDGAxis(fInvmassMCTrue->GetZaxis(), PDGLabel);
    fOutputListMain->Add(fInvmassMCTrue);

    Int_t    PhotBins[5]={14 , 240, 48 ,  4, 2};
    Double_t PhotXMin[5]={0.5,-3., 0.0, -0.5, -0.5};
    Double_t PhotXMax[5]={4. , 3., 12.,  3.5, 1.5};

    fPhotPt1PtMTag = new TH2F("fPhotPt1PtMTag", "fPhotPt1PtMTag; MCPt_e; MCPt_mother", 50, 0.5,3., 120, 0., 12.);
    fOutputListMain->Add(fPhotPt1PtMTag);
    fPhotPt1PtMNTag = new TH2F("fPhotPt1PtMNTag", "fPhotPt1PtMNTag; MCPt_e; MCPt_mother", 50, 0.5,3., 120, 0., 12.);
    fOutputListMain->Add(fPhotPt1PtMNTag);
    
    Double_t Pt2CorrBins[13]={0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 2., 10.};
    fPhotPt2MCRec = new TH2F("fPhotPt2MCRec", "fPhotPt2MCRec; Pt2Rec; Pt2MC", 12, Pt2CorrBins, 12, Pt2CorrBins); //48, 0., 6., 48, 0., 6.);
    fOutputListMain->Add(fPhotPt2MCRec);

    fPhotPt1Pt2 = new THnSparseF("fPhotPt1Pt2",         "fPhotPt1Pt2; Pt1Rec; Pt2MC; PtMotherMC; IdMother; PrimSec;",  5, PhotBins, PhotXMin, PhotXMax);
    fOutputListMain->Add(fPhotPt1Pt2);
    fPhotPt1Pt2Corr = new THnSparseF("fPhotPt1Pt2Corr", "fPhotPt1Pt2Corr; Pt1Rec; Pt2MC; PtMotherMC; IdMother; PrimSec;",  5, PhotBins, PhotXMin, PhotXMax);
    fOutputListMain->Add(fPhotPt1Pt2Corr);
    fPhotPt1Pt2MC = new THnSparseF("fPhotPt1Pt2MC",     "fPhotPt1Pt2MC-woPhotWeights; Pt1Rec; Pt2MC; PtMotherMC; IdMother; PrimSec;",  5, PhotBins, PhotXMin, PhotXMax);
    fOutputListMain->Add(fPhotPt1Pt2MC);
  
    /*
    fPhotPt1RecPt2 = new THnSparseF("fPhotPt1RecPt2", "fPhotPt1RecPt2; Pt1Rec; Pt2; PtMother; IdMother; PrimSec",  5, PhotBins, PhotXMin, PhotXMax);
    fOutputListMain->Add(fPhotPt1RecPt2);
    fPhotPt1RecPt2Corr = new THnSparseF("fPhotPt1RecPt2Corr", "fPhotPt1RecPt2COrr; Pt1Rec; Pt2; PtMother; IdMother; PrimSec",  5, PhotBins, PhotXMin, PhotXMax);
    fOutputListMain->Add(fPhotPt1RecPt2Corr);

    */

    //  fPhotPt1Pt2Rec = new THnSparseF("fPhotPt1Pt2Rec", "fPhotPt1Pt2Rec; Pt1Rec; Pt2Rec; PtMotherMC; IdMother; PrimSec",  5, PhotBins, PhotXMin, PhotXMax);
    // fOutputListMain->Add(fPhotPt1Pt2Rec);

    /*  
    PhotXMin[1]=-0.24;
    PhotXMax[1]=0.24;
    PhotBins[2]=24;
    PhotBins[3]=2;
    PhotBins[4]=1;
    // fPhotPt1E = new THnSparseF("fPhotPt1E", "fPhotPt1; Pt1Rec;InvMass;PtMothMC;IDMoth",  5, PhotBins, PhotXMin, PhotXMax);
    //  fOutputListMain->Add(fPhotPt1E);
    //  fPhotPt1ECorr = new THnSparseF("fPhotPt1ECorr", "fPhotPt1E; Pt1Rec; InvMass; PtMotherMC; IdMother",  5, PhotBins, PhotXMin, PhotXMax);
    //  fOutputListMain->Add(fPhotPt1ECorr);

    PhotBins[2]=60;
    PhotXMin[2]=0.;
    PhotXMax[2]=3.;
    fPhotPt1Pt2E = new THnSparseF("fPhotPt1Pt2EDef", "fPhotPt1Pt2E; Pt1Rec;InvMass;Pt2MC;IDMoth",  5, PhotBins, PhotXMin, PhotXMax);
    fOutputListMain->Add(fPhotPt1Pt2E);
    fPhotPt1Pt2ECorr = new THnSparseF("fPhotPt1Pt2ECorr", "fPhotPt1Pt2E; Pt1; InvMass; Pt2MC; IdMother",  5, PhotBins, PhotXMin, PhotXMax);
    fOutputListMain->Add(fPhotPt1Pt2ECorr);
    */
  }

  if (fOneTimeCheck){
    fOpeningAngleULS = new TH2F("fOpeningAngleULS","fOpeningAngleULS;angle;pt",100,0,1,NBinsElectron,XminElectron,XmaxElectron);
    fOutputListMain->Add(fOpeningAngleULS);
    
    fOpeningAngleLS = new TH2F("fOpeningAngleLS","fOpeningAngleLS;angle;pt",100,0,1,NBinsElectron,XminElectron,XmaxElectron);
    fOutputListMain->Add(fOpeningAngleLS);
  }
  
  fCheckLSULS = new TH2F("fCheckLSULS", "LSULS; #LS; #ULS",5,0,5,5,0,5);
  fOutputListMain->Add(fCheckLSULS);
   
  //  fTagEtaPt1Pt2 = new TH3F("fTagEtaPt1Pt2", "Tagged; Eta; Pt1; Pt2", 36, -0.9, 0.9, NBinsElectron/2, XminElectron, XmaxElectron, 30, 0, 1.5);
  // fOutputList->Add(fTagEtaPt1Pt2);

  /*
  fTagEtaPhiPt = new TH3F("fTagEtaPhiPt", "Tagged: Eta, Phi, Pt", 72, -0.9, 0.9, 64, 0, TMath::TwoPi(), NBinsElectron, XminElectron, XmaxElectron);
  fOutputList->Add(fTagEtaPhiPt);

  fTagEtaZvtxPt = new TH3F("fTagEtaZvtxPt", "Tagged: Eta, Phi, Pt", 72, -0.9, 0.9,  32, -10,10,  NBinsElectron, XminElectron, XmaxElectron);
  fOutputList->Add(fTagEtaZvtxPt);

  fTagEtaPhiPtwW = new TH3F("fTagEtaPhiPtwW", "Tagged: Eta, Phi, Pt", 72, -0.9, 0.9,64, 0, TMath::TwoPi(), NBinsElectron, XminElectron, XmaxElectron);
  fOutputList->Add(fTagEtaPhiPtwW);

  fTagEtaZvtxPtwW= new TH3F("fTagEtaZvtxPtwW", "Tagged: Eta, Phi, Pt", 72, -0.9, 0.9, 32, -10,10, NBinsElectron, XminElectron, XmaxElectron);
  fOutputList->Add(fTagEtaZvtxPtwW);
  */

  if (fIsMC) {
    Int_t NPDGBG = 7;
    Double_t XPDGBG[]={-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5};
    // 0 others, KC 1, KOL 2  K0S 3, w 4, rho 5,  phi 6
    
    fTagEffInclMult = new TH2F("fTagEffInclMult", "fTagEffInclMult; pt; mult", NBinsElectron, XminElectron, XmaxElectron, NMultBins, XMultBins);
    fOutputListMain->Add(fTagEffInclMult);

    fTagEffULSMult = new TH2F("fTagEffULSMult", "fTagEffULSMult; pt; mult", NBinsElectron, XminElectron, XmaxElectron, NMultBins, XMultBins);
    fOutputListMain->Add(fTagEffULSMult);
    
    fTagTruePairsMult = new TH2F("fTagTruePairsMult", "fTagTruePairsMult; pt; mult", NBinsElectron, XminElectron, XmaxElectron, NMultBins, XMultBins);
    fOutputListMain->Add(fTagTruePairsMult);

    fTagEffInclBGMult = new TH3F("fTagEffInclBGMult", "fTagEffInclBGMult; pt; mult;pdg", NBinsElectron, XBinElectronArray, NMultBins, XMultBins, NPDGBG, XPDGBG);
    fOutputListMain->Add(fTagEffInclBGMult);

    fTagEffULSBGMult = new TH3F("fTagEffULSBGMult", "fTagEffULSBGMult; pt; mult; pdg",  NBinsElectron, XBinElectronArray, NMultBins, XMultBins, NPDGBG, XPDGBG);
    fOutputListMain->Add(fTagEffULSBGMult);

    // woW for MCClosure
    fTagEffInclMultWoW = new TH2F("fTagEffInclMultWoW", "fTagEffInclMultWoW; pt; mult", NBinsElectron, XminElectron, XmaxElectron, NMultBins, XMultBins);
    fOutputListMain->Add(fTagEffInclMultWoW);

    fTagEffULSMultWoW = new TH2F("fTagEffULSMultWoW", "fTagEffULSMultWoW; pt; mult", NBinsElectron, XminElectron, XmaxElectron, NMultBins, XMultBins);
    fOutputListMain->Add(fTagEffULSMultWoW);
    
    fTagTruePairsMultWoW = new TH2F("fTagTruePairsMultWoW", "fTagTruePairsMultWoW; pt; mult", NBinsElectron, XminElectron, XmaxElectron, NMultBins, XMultBins);
    fOutputListMain->Add(fTagTruePairsMultWoW);

    fTagEffInclBGMultWoW = new TH3F("fTagEffInclBGMultWoW", "fTagEffInclBGMultWoW; pt; mult",  NBinsElectron, XBinElectronArray, NMultBins, XMultBins, NPDGBG, XPDGBG);
    fOutputListMain->Add(fTagEffInclBGMultWoW);

    fTagEffULSBGMultWoW = new TH3F("fTagEffULSBGMultWoW", "fTagEffULSBGMultWoW; pt; mult",  NBinsElectron, XBinElectronArray, NMultBins, XMultBins, NPDGBG, XPDGBG);
    fOutputListMain->Add(fTagEffULSBGMultWoW);


      // woW for MCClosure
    fTagEffInclMultWoWS = new TH2F("fTagEffInclMultWoWS", "fTagEffInclMultWoWS; pt; mult", NBinsElectron, XminElectron, XmaxElectron, NMultBins, XMultBins);
    fOutputListMain->Add(fTagEffInclMultWoWS);

    fTagEffULSMultWoWS = new TH2F("fTagEffULSMultWoWS", "fTagEffULSMultWoWS; pt; mult", NBinsElectron, XminElectron, XmaxElectron, NMultBins, XMultBins);
    fOutputListMain->Add(fTagEffULSMultWoWS);
    
    fTagTruePairsMultWoWS = new TH2F("fTagTruePairsMultWoWS", "fTagTruePairsMultWoWS; pt; mult", NBinsElectron, XminElectron, XmaxElectron, NMultBins, XMultBins);
    fOutputListMain->Add(fTagTruePairsMultWoWS);

    fTagEffInclBGMultWoWS = new TH3F("fTagEffInclBGMultWoWS", "fTagEffInclBGMultWoWS; pt; mult",  NBinsElectron, XBinElectronArray, NMultBins, XMultBins, NPDGBG, XPDGBG);
    fOutputListMain->Add(fTagEffInclBGMultWoWS);

    fTagEffULSBGMultWoWS = new TH3F("fTagEffULSBGMultWoWS", "fTagEffULSBGMultWoWS; pt; mult", NBinsElectron, XBinElectronArray, NMultBins, XMultBins, NPDGBG, XPDGBG);
    fOutputListMain->Add(fTagEffULSBGMultWoWS);
    
  }

  
  if (fIsMC && fTagEff) {

    // fNonTagEtaPt1Pt2 = new TH3F("fNonTagEtaPt1Pt2", "NonTagged; Eta; Pt1; Pt2", 36, -0.9, 0.9, NBinsElectron/2, XminElectron, XmaxElectron, 30, 0, 1.5);
    // fOutputList->Add(fNonTagEtaPt1Pt2);
    /*
    fNonTagEtaPhiPt = new TH3F("fNonTagEtaPhiPt", "Tagged: Eta, Phi, Pt", 72, -0.9, 0.9, 64, 0, TMath::TwoPi(), NBinsElectron, XminElectron, XmaxElectron);
    fOutputList->Add(fNonTagEtaPhiPt);

    fNonTagEtaZvtxPt = new TH3F("fNonTagEtaZvtxPt", "Tagged: Eta, Phi, Pt", 72, -0.9, 0.9, 32, -10, 10,NBinsElectron, XminElectron, XmaxElectron);
    fOutputList->Add(fNonTagEtaZvtxPt);

    fNonTagEtaPhiPtwW = new TH3F("fNonTagEtaPhiPtwW", "Tagged: Eta, Phi, Pt", 72, -0.9, 0.9,  64, 0, TMath::TwoPi(), NBinsElectron, XminElectron, XmaxElectron);
    fOutputList->Add(fNonTagEtaPhiPtwW);

    fNonTagEtaZvtxPtwW = new TH3F("fNonTagEtaZvtxPtwW", "Tagged: Eta, Phi, Pt", 72, -0.9, 0.9, 32, -10, 10, NBinsElectron, XminElectron, XmaxElectron);
    fOutputList->Add(fNonTagEtaZvtxPtwW);
    */

    /*
    Int_t binMothPt[4]= {NBinsElectron, 50, 10001, 10001};
    Double_t xminMothPt[4]={XminElectron, 0, -0.5, -0.5};
    Double_t xmaxMothPt[4]={XmaxElectron, 25, 9999.5, 9999.5};
    
    fTagMotherPt = new THnSparseF("fTagMotherPt", "Incl: ptElectron, ptMother, Mother, Grandmother",4, binMothPt, xminMothPt, xmaxMothPt);
    fOutputList->Add(fTagMotherPt);
*/
    Int_t binMothPt[4]= {NBinsElectron, 50, PDGSize, PDGSize};
    Double_t xminMothPt[4]={XminElectron, 0, 0.5, 0.5};
    Double_t xmaxMothPt[4]={XmaxElectron, 25, PDGSize+0.5, PDGSize+0.5};
    
    fTagMotherPt = new THnSparseF("fTagMotherPt", "Incl; ptElectron; ptMother; PDGMother; PDGGrandmother;",4, binMothPt, xminMothPt, xmaxMothPt);
    SetPDGAxis(fTagMotherPt->GetAxis(2), PDGLabel);
    SetPDGAxis(fTagMotherPt->GetAxis(3), PDGLabel);
    fOutputListMain->Add(fTagMotherPt);
    /*
    Int_t    binTagEff[5] =  {NBinsElectron   ,10000 ,10001, 10001, 10001}; //p, pdg, pdgmother
    Double_t xminTagEff[5] = {XminElectron    ,-0.5   , -1.5,-1.5, -1.5};
    Double_t xmaxTagEff[5] = {XmaxElectron    ,9999.5, 9999.5,9999.5, 9999.5};  
    */
    // PtElectron, ElectronOrHadron, ElectronMother (Gamma, Pion, Eta, Other), ElectronGrandMother( Eta-Gamma, Pion-Gamma, Other), HFGGMother, isprimary).
    Int_t    binTagEff[6] =  {NBinsElectron   ,2 ,PDGSize, PDGSize, 2, 2}; //p, pdg, pdgmother
    
    Double_t xminTagEff[6] = {XminElectron    ,-0.5   , 0.5,0.5, -0.5, -0.5};
    Double_t xmaxTagEff[6] = {XmaxElectron    ,1.5, PDGSize+0.5,PDGSize+0.5, 1.5, 1.5};  

   

    fTagEffIncl = new THnSparseF("fTagEffIncl", "Incl tag.eff.; Pt; IsElectron; PDGMother; PDGGMother; PDGGGMIsNonHF;IsPrimary;", 6, binTagEff, xminTagEff, xmaxTagEff);
    // SetPDGAxis(fTagEffIncl->GetAxis(1), PDGLabel);
    SetPDGAxis(fTagEffIncl->GetAxis(2), PDGLabel);
    SetPDGAxis(fTagEffIncl->GetAxis(3), PDGLabel);
    // SetPDGAxis(fTagEffIncl->GetAxis(4), PDGLabel);cd
    fOutputListMain->Add(fTagEffIncl);
   
    
    fTagEffULS = new THnSparseF("fTagEffULS", "ULS-LS tag.eff.;  Pt; IsElectron; PDGMother; PDGGMother; PDGGGMIsNonHF;IsPrimary;", 6, binTagEff, xminTagEff, xmaxTagEff);
    //    SetPDGAxis(fTagEffULS->GetAxis(1), PDGLabel);
    SetPDGAxis(fTagEffULS->GetAxis(2), PDGLabel);
    SetPDGAxis(fTagEffULS->GetAxis(3), PDGLabel);
    //    SetPDGAxis(fTagEffULS->GetAxis(4), PDGLabel);
    fOutputListMain->Add(fTagEffULS);

    fTagTruePairs = new THnSparseF("fTagTruePairs", "ULS true pairs tag.eff.;  Pt; IsElectron; PDGMother; PDGGMother; PDGGGMIsNonHF;IsPrimary;", 6, binTagEff, xminTagEff, xmaxTagEff);
    //    SetPDGAxis(fTagTruePairs->GetAxis(1), PDGLabel);
    SetPDGAxis(fTagTruePairs->GetAxis(2), PDGLabel);
    SetPDGAxis(fTagTruePairs->GetAxis(3), PDGLabel);
    //    SetPDGAxis(fTagTruePairs->GetAxis(4), PDGLabel);
    fOutputListMain->Add(fTagTruePairs);

    if (fOneTimeCheck) {
      // fTagEffInclWoWeight = new THnSparseF("fTagEffInclWoWeight", "Incl tag.eff. woW;  Pt; IsElectron; PDGMother; PDGGMother; PDGGGMIsNonHF;IsPrimary;", 6, binTagEff, xminTagEff, xmaxTagEff);
      // //    SetPDGAxis(fTagEffInclWoWeight->GetAxis(1), PDGLabel);
      // SetPDGAxis(fTagEffInclWoWeight->GetAxis(2), PDGLabel);
      // SetPDGAxis(fTagEffInclWoWeight->GetAxis(3), PDGLabel);
      // //    SetPDGAxis(fTagEffInclWoWeight->GetAxis(4), PDGLabel);
      // fOutputListMain->Add(fTagEffInclWoWeight);
  
      // fTagEffULSWoWeight = new THnSparseF("fTagEffULSWoWeight", "ULS-LS tag.eff. woW;  Pt; IsElectron; PDGMother; PDGGMother; PDGGGMIsNonHF;IsPrimary;", 6, binTagEff, xminTagEff, xmaxTagEff);
      // //    SetPDGAxis(fTagEffULSWoWeight->GetAxis(1), PDGLabel);
      // SetPDGAxis(fTagEffULSWoWeight->GetAxis(2), PDGLabel);
      // SetPDGAxis(fTagEffULSWoWeight->GetAxis(3), PDGLabel);
      // //    SetPDGAxis(fTagEffULSWoWeight->GetAxis(4), PDGLabel);
      // fOutputListMain->Add(fTagEffULSWoWeight);

      // fTagTruePairsWoWeight = new THnSparseF("fTagTruePairsWoWeight", "ULS true pairs tag.eff. woW:  Pt; IsElectron; PDGMother; PDGGMother; PDGGGMIsNonHF;IsPrimary;", 6, binTagEff, xminTagEff, xmaxTagEff);
      // //    SetPDGAxis(fTagTruePairsWoWeight->GetAxis(1), PDGLabel);
      // SetPDGAxis(fTagTruePairsWoWeight->GetAxis(2), PDGLabel);
      // SetPDGAxis(fTagTruePairsWoWeight->GetAxis(3), PDGLabel);
      // //    SetPDGAxis(fTagTruePairsWoWeight->GetAxis(4), PDGLabel);
      // fOutputListMain->Add(fTagTruePairsWoWeight);
    }
  }
  

  

  Int_t     bin[5] = {NBinsHadRed   ,NBinsElectron, NBinsDPhi, NBinsDEta, NVertexBins}; //ptH, ptE, Dphi, Deta
  Double_t  xmin[5] = {XminHadron   ,XminElectron ,-TMath::Pi()/2    ,DEtaMin, -10};
  Double_t  xmax[5] = {XmaxHadron   ,XmaxElectron ,(3*TMath::Pi())/2 ,DEtaMax, 10};

  if (fCorrHadron) {

    fNoPartnerNoT = new TH1F("fNoPartnerNoT", "fNoParnterNoT; pt", bin[1], xmin[1], xmax[1]);
    fOutputListHadron->Add(fNoPartnerNoT);

    Int_t TriggerBin[4]={5, 30, 30, fAssPtHad_Nbins};
    Double_t TrigMinBin[4]={0.5, 0., 0., -0.5};
    Double_t TrigMaxBin[4]={3., 1.5, 1.5, fAssPtHad_Nbins-0.5};
    fNoPartnerNoTPt2 = new THnSparseF("fNoPartnerNoTPt2", "fNoParnterNoTPt2; Pt1; Pt2Rec; Pt2MCm; AssHadBin;", 4, TriggerBin, TrigMinBin, TrigMaxBin);
    SetTriggerAxis(fNoPartnerNoTPt2->GetAxis(3));
    fOutputListHadron->Add(fNoPartnerNoTPt2);
  
    fTPartnerNoT = new TH1F("fTPartnerNoT", "fTPartnerNoT; pt", bin[1], xmin[1], xmax[1]);
    fOutputListHadron->Add(fTPartnerNoT);

    fTPartnerNoTPt2 =new THnSparseF("fTPartnerNoTPt2", "fNoParnterNoTPt2; Pt1; Pt2Rec; Pt2MCm; AssHadBin;", 4, TriggerBin, TrigMinBin, TrigMaxBin);
    SetTriggerAxis(fTPartnerNoTPt2->GetAxis(3));
    fOutputListHadron->Add(fTPartnerNoTPt2);

    fElecHadTrigger = new TH3F("fElecHadTrigger", "fElecHadTrigger; ptE; AssPtHad; EtaE", bin[1], xmin[1], xmax[1], fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5, 90, 0, 0.9);
    SetTriggerAxis(fElecHadTrigger->GetYaxis());
    fOutputListHadron->Add(fElecHadTrigger);
    
    fElecHadTriggerULS = new TH2F("fElecHadTriggerULS", "fElecHadTriggerULS-LS; ptE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    SetTriggerAxis(fElecHadTriggerULS->GetYaxis());
    fOutputListHadron->Add(fElecHadTriggerULS);		

    fElecHadTriggerULSNoP = new TH2F("fElecHadTriggerULSNoP", "fElecHadTriggerULSNoP ULS-LS; ptE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    SetTriggerAxis(fElecHadTriggerULSNoP->GetYaxis());
    fOutputListHadron->Add(fElecHadTriggerULSNoP);						      
    fElecHadTriggerULSNoPCorr = new TH2F("fElecHadTriggerULSNoPCorr", "fElecHadTriggerULSNoPCorr ULS-LS; ptE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    SetTriggerAxis(fElecHadTriggerULSNoPCorr->GetYaxis());
    fOutputListHadron->Add(fElecHadTriggerULSNoPCorr);				      		      

    fElecHadTriggerULSNoPCorrTrue = new TH2F("fElecHadTriggerULSNoPCorrTrue", "fElecHadTriggerULSNoPCorrTrue ULS-LS; ptE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    SetTriggerAxis(fElecHadTriggerULSNoPCorrTrue->GetYaxis());
    fOutputListHadron->Add(fElecHadTriggerULSNoPCorrTrue);	

    fHadContTrigger = new TH2F("fHadContTrigger", "fHadContTrigger; PtE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    SetTriggerAxis(fHadContTrigger->GetYaxis());
    fOutputListHadron->Add(fHadContTrigger);						      

    fNonElecHadTrigger = new TH2F("fNonElecHadTrigger", "fNonElecHadTrigger; PtE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    SetTriggerAxis(fNonElecHadTrigger->GetYaxis());
    fOutputListHadron->Add(fNonElecHadTrigger);						      

    fMCElecHaTruePartnerTrigger = new TH2F("fMCElecHaTruePartnerTrigger", "fMCElecHaTruePartnerTrigger; PtE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    SetTriggerAxis(fMCElecHaTruePartnerTrigger->GetYaxis());
    fOutputListHadron->Add(fMCElecHaTruePartnerTrigger);        

    fMCElecHaTruePartnerTriggerWW = new TH2F("fMCElecHaTruePartnerTriggerWW", "fMCElecHaTruePartnerTriggerWW; PtE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    SetTriggerAxis(fMCElecHaTruePartnerTriggerWW->GetYaxis());
    fOutputListHadron->Add(fMCElecHaTruePartnerTriggerWW);      	

    fMCElecHaNoPartnerTrigger = new TH2F("fMCElecHaNoPartnerTrigger", "fMCElecHaNoPartnerTrigger; PtE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    SetTriggerAxis(fMCElecHaNoPartnerTrigger->GetYaxis());
    fOutputListHadron->Add(fMCElecHaNoPartnerTrigger);

    fMCElecHaNoPartnerTriggerWW = new TH2F("fMCElecHaNoPartnerTriggerWW", "fMCElecHaNoPartnerTriggerWW; PtE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    SetTriggerAxis(fMCElecHaNoPartnerTriggerWW->GetYaxis());
    fOutputListHadron->Add(fMCElecHaNoPartnerTriggerWW);

    fHadElecTrigger = new TH2F("fHadElecTrigger", "fHadElecTrigger; PtH; AssPtE", bin[0], xmin[0], xmax[0],  fAssPtElec_Nbins, -0.5, fAssPtElec_Nbins-0.5);
    fHadElecTrigger->GetXaxis()->Set(NBinsHadRed, XBinsHadRed);
    SetTriggerAxis(fHadElecTrigger->GetYaxis());
    fOutputListHadron->Add(fHadElecTrigger);						      

    fHadNonElecTrigger = new TH2F("fHadNonElecTrigger", "fHadNonElecTrigger; PtH; AssPtE", bin[0], xmin[0], xmax[0],  fAssPtElec_Nbins, -0.5, fAssPtElec_Nbins-0.5);
    fHadNonElecTrigger->GetXaxis()->Set(NBinsHadRed, XBinsHadRed);
    SetTriggerAxis(fHadNonElecTrigger->GetYaxis());
    fOutputListHadron->Add(fHadNonElecTrigger);	

    fInclElecHa = new THnSparseF("fEleHaIncl", "Sparse for Ele-Had : PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
    fInclElecHa->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fInclElecHa->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fInclElecHa);

    fULSElecHa = new THnSparseF("fEleHaULS", "Sparse for ULS-LS Ele-Had; PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
    fULSElecHa->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fULSElecHa->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fULSElecHa);
 
    if (fIsMC) { 
      fSignalElecHa = new THnSparseF("fSignalEleHa", "Sparse for ULSEle-Had; PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
      fSignalElecHa->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fSignalElecHa->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListHadron->Add(fSignalElecHa);
     
      fBackgroundElecHa = new THnSparseF("fBackgroundEleHa", "Sparse for ULSEle-Had; PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
      fBackgroundElecHa->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fBackgroundElecHa->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListHadron->Add(fBackgroundElecHa);
          
      fULSElecHaTrue = new THnSparseF("fEleHaULSTrue", "Sparse for ULS-LS True Ele-Had; PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
      fULSElecHaTrue->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fULSElecHaTrue->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListHadron->Add(fULSElecHaTrue);

      fMCElecHaHadron = new THnSparseF("fMCEleHaHadron", "Sparse for Ele-Had; PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
      fMCElecHaHadron->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fMCElecHaHadron->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListHadron->Add(fMCElecHaHadron);
    }
   
    fElecHaULSNoPartner = new THnSparseF("fEleHaULSNoPartner", "Sparse for ULS-LS Ele-Had with no Partner; PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
    fElecHaULSNoPartner->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fElecHaULSNoPartner->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fElecHaULSNoPartner);

    fElecHaULSNoPartnerCorr = new THnSparseF("fEleHaULSNoPartnerCorr", "Sparse for ULS-LS Ele-Had with no Partner; PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
    fElecHaULSNoPartnerCorr->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fElecHaULSNoPartnerCorr->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fElecHaULSNoPartnerCorr);

    fElecHaULSNoPartnerCorrTrue = new THnSparseF("fEleHaULSNoPartnerCorrTrue", "Sparse for ULSEle-Had with no Partner; PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
    fElecHaULSNoPartnerCorrTrue->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fElecHaULSNoPartnerCorrTrue->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fElecHaULSNoPartnerCorrTrue);
    

    if (fIsMC) {
      fMCElecHaTruePartner = new THnSparseF("fMCElecHaTruePartner", "Sparse for MC true photonics with no Partner; PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
      fMCElecHaTruePartner->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fMCElecHaTruePartner->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListHadron->Add(fMCElecHaTruePartner);
    
      fMCElecHaNoPartner = new THnSparseF("fMCElecHaNoPartner", "Sparse for MC true photonics with no Partner: PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
      fMCElecHaNoPartner->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fMCElecHaNoPartner->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListHadron->Add(fMCElecHaNoPartner);	
    }

    fElecHaHa = new THnSparseF("fEleHaHa", "Sparse for Hadron Contamination; PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
    fElecHaHa->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fElecHaHa->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fElecHaHa);

    fElecHaMixedEvent = new THnSparseF("fEleHaMixedEv", "Sparse for Ele-Had MixEvent; PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
    fElecHaMixedEvent->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fElecHaMixedEvent->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fElecHaMixedEvent);
  
    //    fLSElecHaMixedEvent = new THnSparseF("fEleHaLSMixedEv", "Sparse for LSEle-Had MixEvent; PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
    //fLSElecHaMixedEvent->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    //fLSElecHaMixedEvent->GetAxis(4)->Set(NVertexBins, XVertexBins);
    //fOutputListHadron->Add(fLSElecHaMixedEvent);
  
    fULSElecHaMixedEvent = new THnSparseF("fEleHaULSMixedEv", "Sparse for ULS-LS Ele-Had MixEvent; PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
    fULSElecHaMixedEvent->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fULSElecHaMixedEvent->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fULSElecHaMixedEvent);

    if (fIsMC) {
      fTagHaMixedEvent = new THnSparseF("fTagHaMixedEv", "Sparse for Tag Ele-Had MixEvent; PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
      fTagHaMixedEvent->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fTagHaMixedEvent->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListHadron->Add(fTagHaMixedEvent);
      
      fNonTagHaMixedEvent = new THnSparseF("fNonTagHaMixedEv", "Sparse for NonTag Ele-Had MixEvent; PtH; PtE; Dphi; Deta; zVtx;", 5, bin, xmin, xmax);
      fNonTagHaMixedEvent->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fNonTagHaMixedEvent->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListHadron->Add(fNonTagHaMixedEvent);
    }
  }
  // LP HFE   

  if (fCorrLParticle) {


    
    if (fIsMC) {
      Int_t LPBins[4]={NBinsHadRed,NBinsHadRed, PDGSize, PDGSize};
      Double_t LPBinsXmin[4]={XminHadron,XminHadron, 0.5, 0.5};
      Double_t LPBinsXmax[4]={XmaxHadron,XmaxHadron, PDGSize+0.5, PDGSize+0.5};
      fMCLeadingParticle = new THnSparseF("fMCLeadingParticle", "fMCLeadingParticle; pt; mcPt; pdg; pdgmother", 4, LPBins, LPBinsXmin, LPBinsXmax);
      fMCLeadingParticle->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fMCLeadingParticle->GetAxis(1)->Set(NBinsHadRed, XBinsHadRed);
      SetPDGAxis(fMCLeadingParticle->GetAxis(2), PDGLabel);
      SetPDGAxis(fMCLeadingParticle->GetAxis(3), PDGLabel);
      fOutputListLP->Add(fMCLeadingParticle);
      
      fCompareLPRecCheck = new TH3F("fCompareLPRecCheck", "fMCLPRecCheck; #Delta pt; #Delta phi; #Delta eta", 50, -1, 4, 36, -TMath::Pi()/2, (3.*TMath::Pi())/2, 40, -4, 4);
      fOutputListLP->Add(fCompareLPRecCheck);
    }

    fElecLPTrigger = new TH3F("fElecLPTrigger", "fElecLPTrigger", bin[1], xmin[1], xmax[1], fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5, 90, 0, 0.9);
    fOutputListLP->Add(fElecLPTrigger);	
						      
    //    fElecLPTriggerLS = new TH2F("fElecLPTriggerLS", "fElecLPTriggerLS", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    // fOutputListLP->Add(fElecLPTriggerLS);						      

    fElecLPTriggerULS = new TH2F("fElecLPTriggerULS", "fElecLPTriggerULS-LS", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListLP->Add(fElecLPTriggerULS);	
    
    //    fElecLPTriggerLSNoP = new TH2F("fElecLPTriggerLSNoP", "fElecLPTriggerLSNoP", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    //    fOutputListLP->Add(fElecLPTriggerLSNoP);						      

    fElecLPTriggerULSNoP = new TH2F("fElecLPTriggerULSNoP", "fElecLPTriggerULS-LSNoP", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListLP->Add(fElecLPTriggerULSNoP);

    fElecLPTriggerULSNoPCorr = new TH2F("fElecLPTriggerULSNoPCorr", "fElecLPTriggerULSNoPCorr ULS-LS; ptE; AssPtLP", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListLP->Add(fElecLPTriggerULSNoPCorr);				      		      

    //  fElecLPTriggerULSNoPCorrTrue = new TH2F("fElecLPTriggerULSNoPCorrTrue", "fElecLPTriggerULSNoPCorrTrue ULS-LS; ptE; AssPtLP", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    // fOutputListLP->Add(fElecLPTriggerULSNoPCorrTrue);	

    

    fHadContLPTrigger = new TH2F("fHadContLPTrigger", "fHadContLPTrigger", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListLP->Add(fHadContLPTrigger);		

    fNonElecLPTrigger = new TH2F("fNonElecLPTrigger", "fNonElecLPTrigger", bin[1], xmin[1], xmax[1], fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListLP->Add(fNonElecLPTrigger);					      
					      
    fMCElecLPTruePartnerTrigger = new TH2F("fMCElecLPTruePartnerTrigger", "fMCElecLPTruePartnerTrigger", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListLP->Add(fMCElecLPTruePartnerTrigger);	

    fMCElecLPNoPartnerTrigger = new TH2F("fMCElecLPNoPartnerTrigger", "fMCElecLPNoPartnerTrigger", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListLP->Add(fMCElecLPNoPartnerTrigger);

    fLPElecTrigger = new TH2F("fLPElecTrigger", "fLPElecTrigger", bin[0], xmin[0], xmax[0],  fAssPtElec_Nbins, -0.5, fAssPtElec_Nbins-0.5);
    fLPElecTrigger->GetXaxis()->Set(NBinsHadRed, XBinsHadRed);
    fOutputListLP->Add(fLPElecTrigger);						      

    fLPNonElecTrigger = new TH2F("fLPNonElecTrigger", "fLPNonElecTrigger", bin[0], xmin[0], xmax[0], fAssPtElec_Nbins, -0.5, fAssPtElec_Nbins-0.5);
    fLPNonElecTrigger->GetXaxis()->Set(NBinsHadRed, XBinsHadRed);
    fOutputListLP->Add(fLPNonElecTrigger);						      


    fInclElecLP = new THnSparseF("fEleLPIncl", "Sparse for Ele-LP ; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fInclElecLP->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fInclElecLP->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListLP->Add(fInclElecLP);
    
    //fLSElecLP = new THnSparseF("fEleLPLS", "Sparse for LSEle-LP ; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    //fLSElecLP->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    //fLSElecLP->GetAxis(4)->Set(NVertexBins, XVertexBins);
    //fOutputListLP->Add(fLSElecLP);
    
    fULSElecLP = new THnSparseF("fEleLPULS", "Sparse for ULS-LS Ele-LP ; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fULSElecLP->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fULSElecLP->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListLP->Add(fULSElecLP);

    fMCElecLPHadron = new THnSparseF("fMCEleLPHadron", "Sparse for Ele-Had ; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fMCElecLPHadron->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fMCElecLPHadron->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListLP->Add(fMCElecLPHadron);
    
    //    fElecLPLSNoPartner = new THnSparseF("fEleLPLSNoPartner", "Sparse for LSEle-Had with no Partner; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    //fElecLPLSNoPartner->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    //fElecLPLSNoPartner->GetAxis(4)->Set(NVertexBins, XVertexBins);
    //fOutputListLP->Add(fElecLPLSNoPartner);
    
    fElecLPULSNoPartner = new THnSparseF("fEleLPULSNoPartner", "Sparse for ULSEle-Had with no Partner; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fElecLPULSNoPartner->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fElecLPULSNoPartner->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListLP->Add(fElecLPULSNoPartner);


    fElecLPULSNoPartnerCorr = new THnSparseF("fEleLPULSNoPartnerCorr", "Sparse for ULS-LS Ele-LP with no Partner; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fElecLPULSNoPartnerCorr->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fElecLPULSNoPartnerCorr->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListLP->Add(fElecLPULSNoPartnerCorr);

    fElecLPULSNoPartnerCorrTrue = new THnSparseF("fEleLPULSNoPartnerCorrTrue", "Sparse for ULSEle-LP with no Partner; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fElecLPULSNoPartnerCorrTrue->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fElecLPULSNoPartnerCorrTrue->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListLP->Add(fElecLPULSNoPartnerCorrTrue);
    
    
    fElecLPHa = new THnSparseF("fEleLPHa", "Sparse for Hadron Contamination; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fElecLPHa->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fElecLPHa->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListLP->Add(fElecLPHa);

    if (fIsMC) {
      fMCElecLPTruePartner = new THnSparseF("fMCElecLPTruePartner", "Sparse for MC true photonics with true Partner; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
      fMCElecLPTruePartner->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fMCElecLPTruePartner->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListLP->Add(fMCElecLPTruePartner);

      fMCElecLPNoPartner = new THnSparseF("fMCElecLPNoPartner", "Sparse for MC true photonics with no Partner; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
      fMCElecLPNoPartner->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fMCElecLPNoPartner->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListLP->Add(fMCElecLPNoPartner);	
    }

    fElecLPMixedEvent = new THnSparseF("fEleLPMixedEv", "Sparse for Ele-LP MixEvent; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fElecLPMixedEvent->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fElecLPMixedEvent->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListLP->Add(fElecLPMixedEvent);
    
    //    fLSElecLPMixedEvent = new THnSparseF("fEleLPLSMixedEv", "Sparse for LSEle-LP MixEvent; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    //fLSElecLPMixedEvent->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    //fLSElecLPMixedEvent->GetAxis(4)->Set(NVertexBins, XVertexBins);
    //fOutputListLP->Add(fLSElecLPMixedEvent);
    
    fULSElecLPMixedEvent = new THnSparseF("fEleLPULSMixedEv", "Sparse for ULS-LS Ele-LP MixEvent; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fULSElecLPMixedEvent->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fULSElecLPMixedEvent->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListLP->Add(fULSElecLPMixedEvent);
    
    if (fIsMC) {
      fTagLPMixedEvent = new THnSparseF("fTagLPMixedEv", "Sparse for Tag-LP MixEvent; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
      fTagLPMixedEvent->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fTagLPMixedEvent->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListLP->Add(fTagLPMixedEvent);
      
      fNonTagLPMixedEvent = new THnSparseF("fNonTagLPMixedEv", "Sparse for NonTag-LP MixEvent; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
      fNonTagLPMixedEvent->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fNonTagLPMixedEvent->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListLP->Add(fNonTagLPMixedEvent);

    }


  }





  //////////// Hadron and Electron RecEff //////////////////////

  Int_t    EffHBins[4]={NBinsHadRed, 16, 16, NVertexBins};
  Double_t EffHXmin[4]={XminHadron, -0.8, 0, -10};
  Double_t EffHXmax[4]={XmaxHadron, 0.8, TMath::TwoPi(), 10};

  Int_t    EffEBins[4]={NBinsElectronRecEff, 16, 16, NVertexBins};
  Double_t EffEXmin[4]={XminElectron, -0.8, 0, -10};
  Double_t EffEXmax[4]={100, 0.8, TMath::TwoPi(), 10};
    
  if (!fOneTimeCheck) { // maybe reduce axis further
    EffHBins[1]=9; // hadrons for pt, eta, zVtx
    EffHBins[2]=1; // hadrons for pt, eta, zVtx
    EffEBins[1]=2; // electrons only pt, zVtx
    EffEBins[2]=1;
  }

 
  // each time count corrected hadrons and electrons (check for right efficiency)
  fRecHadPtEtaPhiVtxWRecEff=new THnSparseF("fRecHadPtEtaPhiVtxWRecEff", "Rec hadrons wRecEff; rec #it{p}_{T}; rec #eta; rec #varphi; rec zVtx;", 4, EffHBins, EffHXmin, EffHXmax);
  fRecHadPtEtaPhiVtxWRecEff->GetAxis(3)->Set(NVertexBins, XVertexBins);
  fRecHadPtEtaPhiVtxWRecEff->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
  fOutputListMain->Add(fRecHadPtEtaPhiVtxWRecEff);   

 

  if (fIsMC) {

    //fRecHFE = new TH1F("fRecHFE", "fRecHFE; #it{p}_{T}", NBinsElectronRed, XminElectron, XmaxElectron);
    // fOutputListMain->Add(fRecHFE);
     
    //fRecHFEEtaWRecEff = new TH2F("fRecHFEEtaWRecEff", "fRecHFEEtaWRecEff; Pt; Eta", NBinsElectron , XminElectron, XmaxElectron, 250, -2.5, 2.5); 
    //fOutputList->Add(fRecHFEEtaWRecEff);

    // accordingly check truth generated hadrons and electrons
    fMCHadPtEtaPhiVtx=new THnSparseF("fMCHadPtEtaPhiVtx", "MC truth gen. hadrons; MC #it{p}_{T}; MC #eta; MC #varphi; MC zVtx;", 4, EffHBins, EffHXmin, EffHXmax);
    fMCHadPtEtaPhiVtx->GetAxis(3)->Set(NVertexBins, XVertexBins);
    fMCHadPtEtaPhiVtx->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fOutputListMain->Add(fMCHadPtEtaPhiVtx);

    fRecElecPtEtaPhiVtxWRecEff=new THnSparseF("fRecElePtEtaPhiWRecEff", "Rec electrons w. rec.Eff; rec #it{p}_{T}; rec #eta; rec #varphi; rec zVtx;", 4, EffEBins, EffEXmin, EffEXmax);
    fRecElecPtEtaPhiVtxWRecEff->GetAxis(3)->Set(NVertexBins, XVertexBins);
    fRecElecPtEtaPhiVtxWRecEff->GetAxis(0)->Set(NBinsElectronRecEff, XBinsElectronRecEff);
    fOutputListMain->Add(fRecElecPtEtaPhiVtxWRecEff);
    
    fMCElecPtEtaPhiVtx=new THnSparseF("fMCElePtEtaPhiVtx", "MC truth gen. electrons;  MC #it{p}_{T}; MC #eta; MC #varphi; MC zVtx;", 4, EffEBins, EffEXmin, EffEXmax);
    fMCElecPtEtaPhiVtx->GetAxis(3)->Set(NVertexBins, XVertexBins);
    fMCElecPtEtaPhiVtx->GetAxis(0)->Set(NBinsElectronRecEff, XBinsElectronRecEff);
    fOutputListMain->Add(fMCElecPtEtaPhiVtx);
   
    if (fRecEff) { // for rec eff count rec. hadrons and electrons wo. receff correction
      fRecHadPtEtaPhiVtx=new THnSparseF("fRecHadPtEtaPhiVtx", "Rec hadrons; rec #it{p}_{T}; rec #eta; rec #varphi; rec zVtx;", 4, EffHBins, EffHXmin, EffHXmax);
      fRecHadPtEtaPhiVtx->GetAxis(3)->Set(NVertexBins, XVertexBins);
      fRecHadPtEtaPhiVtx->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fOutputListMain->Add(fRecHadPtEtaPhiVtx);
      
      fRecElecPtEtaPhiVtx=new THnSparseF("fRecElePtEtaPhi", "Rec electrons; rec #it{p}_{T}; rec #eta; rec #varphi; rec zVtx;", 4, EffEBins, EffEXmin, EffEXmax);
      fRecElecPtEtaPhiVtx->GetAxis(3)->Set(NVertexBins, XVertexBins);
      fRecElecPtEtaPhiVtx->GetAxis(0)->Set(NBinsElectronRecEff, XBinsElectronRecEff);
      fOutputListMain->Add(fRecElecPtEtaPhiVtx);

      if (fOneTimeCheck) { // check mc pt .. for rec. hadrons and electrons wo receff correction      
	fRecHadMCPtEtaPhiVtx=new THnSparseF("fRecHadMCPtEtaPhiVtx", "MC rec hadrons w. MC ; MC #it{p}_{T}; MC #eta; MC #varphi; MC zVtx;", 4, EffHBins, EffHXmin, EffHXmax);
	fRecHadMCPtEtaPhiVtx->GetAxis(3)->Set(NVertexBins, XVertexBins);
	fRecHadMCPtEtaPhiVtx->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
	fOutputListMain->Add(fRecHadMCPtEtaPhiVtx);    
	
	fRecElecMCPtEtaPhiVtx=new THnSparseF("fRecEleMCPtEtaPhi", "MC rec electrons w. MC; MC #it{p}_{T}; MC #eta; MC #varphi; MC zVtx;", 4, EffEBins, EffEXmin, EffEXmax);
	fRecElecMCPtEtaPhiVtx->GetAxis(3)->Set(NVertexBins, XVertexBins);
	fRecElecMCPtEtaPhiVtx->GetAxis(0)->Set(NBinsElectronRecEff, XBinsElectronRecEff);
	fOutputListMain->Add(fRecElecMCPtEtaPhiVtx);

	fCheckMCPtvsRecPtHad = new TH2F("fCheckMCPtvsRecPtHad", "MCvsRec Pt; MC pt; Rec pt", NBinsHadron, XminHadron, XmaxHadron, NBinsHadron, XminHadron, XmaxHadron); 
	fOutputListMain->Add(fCheckMCPtvsRecPtHad);
	
	fCheckMCEtavsRecEtaHad = new TH2F("fCheckMCEtavsRecEtaHad", "MC vs. Rec; rec. #eta; MC #eta", 36, -0.9, 0.9, 36, -0.9, 0.9);
	fOutputListMain->Add(fCheckMCEtavsRecEtaHad);
	
	fCheckMCPhivsRecPhiHad = new TH2F("fCheckMCPhivsRecPhiHad", "MC vs. Rec; rec. #eta; MC #eta", NBinsPhi/2, 0, TMath::TwoPi(), NBinsPhi/2, 0, TMath::TwoPi());
	fOutputListMain->Add(fCheckMCPhivsRecPhiHad);
	
	fCheckMCPtvsRecPtEle = new TH2F("fCheckMCPtvsRecPtEle", "MCvsRec Pt", NBinsElectron, XminElectron, XmaxElectron, NBinsElectron, XminElectron, XmaxElectron); //200, 0, 20, 200, 0, 20);
	fOutputListMain->Add(fCheckMCPtvsRecPtEle);
	
	fMCElecPDG=new TH1F("fMCElePDG", "MC truth mother of heavy electrons", 10000, -0.5, 9999.5);
	fOutputListMain->Add(fMCElecPDG);
	
	// fMCElecPtEtaPhiStrictVtx=new THnSparseF("fMCElePtEtaPhiStrictVtx", "MC truth electrons  pt, eta, phi, vtx", 4, EffEBins, EffEXmin, EffEXmax);
	// fMCElecPtEtaPhiStrictVtx->GetAxis(3)->Set(NVertexBins, XVertexBins);
	// fMCElecPtEtaPhiStrictVtx->GetAxis(0)->Set(NBinsElectronRecEff, XBinsElectronRecEff);
	// fOutputListMain->Add(fMCElecPtEtaPhiStrictVtx);
      
      }
    }
    fRecHadMCSecondaryCont = new TH2F("fRecHadMCSecondaryCont", "Cont. from secondaries; #it{p}_{T}; Prim or Sec", NBinsHadRed, XBinsHadRed, 2, -0.5, 1.5);
    fRecHadMCSecondaryCont->GetYaxis()->SetBinLabel(1, "Secondary");
    fRecHadMCSecondaryCont->GetYaxis()->SetBinLabel(2, "Primary");
    fOutputListMain->Add(fRecHadMCSecondaryCont);

    fRecElecMCSecondaryCont = new TH2F("fRecElecMCSecondaryCont", "Cont. from secondaries; #it{p}_{T}; Prim or Sec", NBinsElectronRecEff, XBinsElectronRecEff, 2, -0.5, 1.5);
    fRecElecMCSecondaryCont->GetYaxis()->SetBinLabel(1, "Secondary");
    fRecElecMCSecondaryCont->GetYaxis()->SetBinLabel(2, "Primary");
    fOutputListMain->Add(fRecElecMCSecondaryCont);
  
    
  }


  Int_t Pi0EtaBins[5] ={100, 4, 40, PDGSize,5}; // pt,eta,y, 
  Double_t Pi0EtaXmin[5]={0.1, -2,-2, 0.5, 0};
  Double_t Pi0EtaXmax[5]={25, 2,2,  PDGSize+0.5,5};

  if (fIsMC && fPionEtaProduction) {
    fMCPi0Prod = new THnSparseF("fMCPi0Prod", "fMCPi0Prod; Pt; Eta; Y; PDGMoth; Cocktail", 5, Pi0EtaBins, Pi0EtaXmin, Pi0EtaXmax);
    BinLogX(fMCPi0Prod->GetAxis(0));
    SetPDGAxis(fMCPi0Prod->GetAxis(3), PDGLabel);
    fOutputListMain->Add(fMCPi0Prod);

    fMCEtaProd = new THnSparseF("fMCEtaProd", "fMCEtaProd; Pt; Eta; Y; PDGMoth; Cocktail", 5, Pi0EtaBins, Pi0EtaXmin, Pi0EtaXmax);
    BinLogX(fMCEtaProd->GetAxis(0));
    SetPDGAxis(fMCEtaProd->GetAxis(3), PDGLabel);
    fOutputListMain->Add(fMCEtaProd);

    fMCPiPlusProd = new THnSparseF("fMCPiPlusProd", "fMCPiPlusProd; Pt; Eta; Y; PDGMoth; Cocktail;", 5, Pi0EtaBins, Pi0EtaXmin, Pi0EtaXmax);
    BinLogX(fMCPiPlusProd->GetAxis(0));
    SetPDGAxis(fMCPiPlusProd->GetAxis(3), PDGLabel);
    fOutputListMain->Add(fMCPiPlusProd);

    fMCPiPlusProdV2 = new THnSparseF("fMCPiPlusProdV2", "fMCPiPlusProdV2; Pt; Eta; Y; PDGMoth; Cocktail;", 5, Pi0EtaBins, Pi0EtaXmin, Pi0EtaXmax);
    BinLogX(fMCPiPlusProdV2->GetAxis(0));
    SetPDGAxis(fMCPiPlusProdV2->GetAxis(3), PDGLabel);
    fOutputListMain->Add(fMCPiPlusProdV2);

    Pi0EtaBins[4] = Pi0EtaBins[3];
    Pi0EtaXmin[4] = Pi0EtaXmin[3];
    Pi0EtaXmax[4] = Pi0EtaXmax[3];
    fMCBGProd = new THnSparseF("fMCBGProd", "fMCBGProd; Pt; Eta; Y; PDG; PDGMoth;", 5, Pi0EtaBins, Pi0EtaXmin, Pi0EtaXmax);
    BinLogX(fMCBGProd->GetAxis(0));
    SetPDGAxis(fMCBGProd->GetAxis(3), PDGLabel);
    SetPDGAxis(fMCBGProd->GetAxis(4), PDGLabel);
    fOutputListMain->Add(fMCBGProd);
  }




  if (fMCTrueCorrelation) {
    Int_t    poolSize = 100;
    Int_t    trackDepth = 2000; 
    const Int_t    nMaxPtBins=5;
    Double_t maxPtBins[]={0, 1000, 1001, 1002, 1003, 9999}; //{0.5, 2., 5., 10, 9999};

    //  8 Vertex bins * 3 MultBins, * 5 PtBins
    fMCTruePoolMgr = new AliEventPoolManager(poolSize, trackDepth, NMultBinsRed, (Double_t*) XMultBinsRed, NVertexBins, (Double_t*) XVertexBins, nMaxPtBins,(Double_t*)maxPtBins);
    fMCTruePoolMgr->Validate();  


    Int_t MCTrueBins[5]= {NBinsHadRed   ,NBinsElectron, NBinsDPhi, NBinsDEta, 10}; 
    Double_t MCTrueXmin[5] = {XminHadron   ,XminElectron ,-TMath::Pi()/2    ,DEtaMin, 0.5};
    Double_t MCTrueXmax[5] ={XmaxHadron   ,XmaxElectron ,(3*TMath::Pi())/2 ,DEtaMax, 10.5};

    Int_t MCTrueMEvBins[5]= {NBinsHadRed   ,NBinsElectron, NBinsDPhi, NBinsDEta, NVertexBins}; 
    Double_t MCTrueMEvXmin[5] = {XminHadron   ,XminElectron ,-TMath::Pi()/2    ,DEtaMin,-10};
    Double_t MCTrueMEvXmax[5] ={XmaxHadron   ,XmaxElectron ,(3*TMath::Pi())/2 ,DEtaMax, 10};

    
    fCompareLP = new TH3F("fCompareLP", "fComparelP; RecLP; LPinAcc; LP", NBinsHadRed, XBinsHadRed, NBinsHadRed, XBinsHadRed,NBinsHadRed, XBinsHadRed);
    fOutputList->Add(fCompareLP);

   
    // fRecHadronEtaWRecEff = new TH2F("fRecHadronEtaWRecEff", "fRecHadronEtaWRecEff; Pt; Eta", NBinsHadron , XminHadron, XmaxHadron, 250, -2.5, 2.5); 
    //fOutputList->Add(fRecHadronEtaWRecEff);

    fRecLPEta = new TH2F("fRecLPEta", "fRecLPEta; Pt; Eta", NBinsHadRed , XminHadron, XmaxHadron, 250, -2.5, 2.5); 
    fRecLPEta->GetXaxis()->Set(NBinsHadRed, XBinsHadRed);
    fOutputList->Add(fRecLPEta);

    if (fIsMC) {
      fTrueElectronEta = new TH3F("fTrueElectronEta", "fTrueElectronEta; Pt; Eta; mult", NBinsElectron , XminElectron, XmaxElectron, 250, -2.5, 2.5, 200, -0.5, 199.5); 
      fOutputList->Add(fTrueElectronEta);
	
    }


    if (fIsMC && (fCorrHadron || fCorrLParticle)) {
      fTrueMCHadronEventCutsZvtx = new THnSparseF("fTrueMCHadronEventCutsZvtx", "fTrueMCHadronEventCutsZvtx: ptH, ptE, dphi, deta, zVtx", 5, MCTrueMEvBins, MCTrueMEvXmin, MCTrueMEvXmax);
      fTrueMCHadronEventCutsZvtx->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fTrueMCHadronEventCutsZvtx->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputList->Add(fTrueMCHadronEventCutsZvtx);

      fTrueMCHadronEventCutsZvtxMEv = new THnSparseF("fTrueMCHadronEventCutsZvtMEvx", "fTrueMCHadronEventCutsZvtxMEv: ptH, ptE, dphi, deta, zVtx", 5, MCTrueMEvBins, MCTrueMEvXmin, MCTrueMEvXmax);
      fTrueMCHadronEventCutsZvtxMEv->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fTrueMCHadronEventCutsZvtxMEv->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputList->Add(fTrueMCHadronEventCutsZvtxMEv);

      fTrueMCElecHaTriggerEventCuts = new TH3F("fMCTrueEleHaTriggerEvCuts", "fMCTrueEleHaTriggerEvCuts: pt, case, assbin",MCTrueBins[1], MCTrueXmin[1], MCTrueXmax[1], 10, 0.5, 10.5,  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5); 
      fOutputList->Add(fTrueMCElecHaTriggerEventCuts);

      fTrueMCLPEventCutsZvtx = new THnSparseF("fTrueMCLPEventCutsZvtx", "fTrueMCLPEventCutsZvtx: ptH, ptE, dphi, deta, zVtx", 5, MCTrueMEvBins, MCTrueMEvXmin, MCTrueMEvXmax);
      fTrueMCLPEventCutsZvtx->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fTrueMCLPEventCutsZvtx->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputList->Add(fTrueMCLPEventCutsZvtx);

      fTrueMCLPEventCutsZvtxMEv = new THnSparseF("fTrueMCLPEventCutsZvtMEvx", "fTrueMCLPEventCutsZvtxMEv: ptH, ptE, dphi, deta, zVtx", 5, MCTrueMEvBins, MCTrueMEvXmin, MCTrueMEvXmax);
      fTrueMCLPEventCutsZvtxMEv->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fTrueMCLPEventCutsZvtxMEv->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputList->Add(fTrueMCLPEventCutsZvtxMEv);

      fTrueMCElecLPTriggerEventCuts = new TH3F("fMCTrueEleLPTriggerEvCuts", "fMCTrueEleLPTriggerEvCuts: pt, case, assbin",MCTrueBins[1], MCTrueXmin[1], MCTrueXmax[1], 10, 0.5, 10.5,  fAssPtHad_Nbins,-0.5, fAssPtHad_Nbins-0.5); 
      fOutputList->Add(fTrueMCElecLPTriggerEventCuts);
     
      if (fOneTimeCheck) {

	fTrueHadronEta = new TH2F("fTrueHadronEta", "fTrueHadronEta; Pt; Eta", NBinsHadron , XminHadron, XmaxHadron, 250, -2.5, 2.5); 
	fOutputList->Add(fTrueHadronEta);
	
	fTrueLPinAcceptanceEta= new TH2F("fTrueLPinAcceptanceEta", "fTrueLPinAcceptanceEta; Pt; Eta", NBinsHadRed , XminHadron, XmaxHadron, 250, -2.5, 2.5); 
	fTrueLPinAcceptanceEta->GetXaxis()->Set(NBinsHadRed, XBinsHadRed);
	fOutputList->Add(fTrueLPinAcceptanceEta);
	
	fTrueLPEta = new TH2F("fTrueLPEta", "fTrueLPEta; Pt; Eta", NBinsHadRed , XminHadron, XmaxHadron, 250, -2.5, 2.5); 
	fTrueLPEta->GetXaxis()->Set(NBinsHadRed, XBinsHadRed);
	fOutputList->Add(fTrueLPEta);

	fTrueMCHadronEventCuts = new THnSparseF("fMCTrueHadronEvCuts", "fMCTrueHadronEvCuts; ptH; ptE; dphi; deta; case", 5, MCTrueBins, MCTrueXmin, MCTrueXmax);
	fTrueMCHadronEventCuts->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
	fOutputList->Add(fTrueMCHadronEventCuts);

	fTrueMCHadron  = new THnSparseF("fMCTrueHadron", "fMCTrueHadron; ptH; ptE; dphi; deta; case", 5, MCTrueBins, MCTrueXmin, MCTrueXmax);
	fTrueMCHadron->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
	fOutputList->Add(fTrueMCHadron);

	fTrueMCElecHaTrigger = new TH3F("fMCTrueEleHaTrigger", "fMCTrueEleHaTrigger: pt, case, assbin",MCTrueBins[1], MCTrueXmin[1], MCTrueXmax[1], 10, 0.5, 10.5, fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
	fOutputList->Add(fTrueMCElecHaTrigger);

	fTrueMCLPEventCuts = new THnSparseF("fMCTrueLPEvCuts", "fMCTrueLPEvCuts; ptH; ptE; dphi; deta; case", 5, MCTrueBins, MCTrueXmin, MCTrueXmax);
	fTrueMCLPEventCuts->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
	fOutputList->Add(fTrueMCLPEventCuts);

	fTrueMCLP  = new THnSparseF("fMCTrueLP", "fMCTrueLP; ptH; ptE; dphi; deta; case", 5, MCTrueBins, MCTrueXmin, MCTrueXmax);
	fTrueMCLP->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
	fOutputList->Add(fTrueMCLP);

	fTrueMCElecLPTrigger = new TH3F("fMCTrueEleLPTrigger", "fMCTrueEleLPTrigger: pt, case, assbin",MCTrueBins[1], MCTrueXmin[1], MCTrueXmax[1], 10, 0.5, 10.5,  fAssPtHad_Nbins,-0.5, fAssPtHad_Nbins-0.5);
	fOutputList->Add(fTrueMCElecLPTrigger);

      }
    }
  }




  if (fTRDQA) {
    fhArmenteros  = new TH2F("fhArmenteros","Armenteros plot",200,-1.,1.,200,0.,0.4);
    fOutputListQA->Add(fhArmenteros);

    fEventsPerRun = new TH1F("fEventsPerRun", "EventsPerRun",500, 0.5, 500.5);
    fOutputListQA->Add(fEventsPerRun);
    
    fTRDnTrackRun = new TH2F("fTRDnTrackRun", "TrackletCharge, TrackletPID", 6, 0.5, 6.5, 6, 0.5, 6.5);
    fOutputListQA->Add(fTRDnTrackRun);
    
    Int_t TRDBins[5]={11, 30, 54, 7, 3};
    Double_t TRDBMin[5]={0.5, XminEta, 0, -0.5, -1.5};
    Double_t TRDBMax[5]={XmaxElectron, XmaxEta, TMath::TwoPi(), 6.5, 1.5};

    fTRDEtaPhi = new THnSparseF("fTRDEtaPhi", "fTRDEtaPhi: Pt, Eta, Phi, Layer, Charge", 6, TRDBins, TRDBMin, TRDBMax);
    fOutputListQA->Add(fTRDEtaPhi);

    Int_t TRDBins3[7]=   {11             ,   9  ,      3,    7,     7,    3, 500};
    Double_t TRD3BMin[7]={0.5         , XminEta,   0.5,   -3,  -0.5, -1.5, 0.5};
    Double_t TRD3BMax[7]={XmaxElectron, XmaxEta,   3.5,    4,   6.5,  1.5, 500.5};

    fTRDV0NTracklets = new THnSparseF("fTRDV0NTracklets", "fTRDnTracklets: Pt, Eta, PID, TPC, NTracklets, Charge", 7, TRDBins3, TRD3BMin, TRD3BMax);
    fOutputListQA->Add(fTRDV0NTracklets);
  
    fTRDNTracklets = new THnSparseF("fTRDNTracklets", "fTRDnTracklets: Pt, Eta, PID, TPC, NTracklets, Charge",  7, TRDBins3, TRD3BMin, TRD3BMax);
    fOutputList->Add(fTRDNTracklets);

    Int_t TRDBins2[5]={NBinsElectron/2,    7,   3,   7, 3};
    Double_t TRDB2Min[5]={XminElectron, -0.5, 0.5,  -3, -1.5};
    Double_t TRDB2Max[5]={XmaxElectron,  6.5, 3.5,   4,  1.5};   

    fTRDV0Spectra = new THnSparseF("fTRDV0Spectra", "fTRDV0Spectra: Pt, NTracklets, PID, TPC, Charge", 5, TRDBins2, TRDB2Min, TRDB2Max);
    fOutputList->Add(fTRDV0Spectra);

    fTRDSpectra = new THnSparseF("fTRDSpectra", "fTRDSpectra: Pt, NTracklets, PID, TPC, Charge", 5, TRDBins2, TRDB2Min, TRDB2Max); 
    fOutputList->Add(fTRDSpectra);

    fTRDMCSpectra= new THnSparseF("fTRDMCSpectra","fTRDMCSpectra: Pt, NTracklets, PID, TPC, Charge", 5, TRDBins2, TRDB2Min, TRDB2Max);  
    fOutputList->Add(fTRDMCSpectra);
  }

  Int_t    poolSize = 100;
  Int_t    trackDepth = 2000; 
  const Int_t    nMaxPtBins=5;
  Double_t maxPtBins[]={0,1000, 1001, 1002, 1003, 9999}; //{ 0.5, 2., 5., 10, 9999};

  //  8 Vertex bins * 3 MultBins, * 5 PtBins
  fPoolMgr = new AliEventPoolManager(poolSize, trackDepth, NMultBinsRed, (Double_t*) XMultBinsRed, NVertexBins, (Double_t*) XVertexBins, nMaxPtBins,(Double_t*)maxPtBins);
  fPoolMgr->Validate();
    
  fPoolIsFilled = new TH3F("fPoolIsFilled", "fPoolIsFilled: mult, vertex, LP pt", NMultBinsRed, (Double_t*) XMultBins, NVertexBins, (Double_t*) XVertexBins, nMaxPtBins,(Double_t*)maxPtBins);
  fOutputListQA->Add(fPoolIsFilled);
  
  for (Int_t i=0; i < fOutputList->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutputList->At(i));
    if (h1) h1->Sumw2();
    THnSparse *hSparse = dynamic_cast<THnSparse*>(fOutputList->At(i));
    if (hSparse)  hSparse->Sumw2();
  }
  for (Int_t i=0; i < fOutputListMain->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutputListMain->At(i));
    if (h1) h1->Sumw2();
    THnSparse *hSparse = dynamic_cast<THnSparse*>(fOutputListMain->At(i));
    if (hSparse)  hSparse->Sumw2();
  }
  for (Int_t i=0; i < fOutputListLP->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutputListLP->At(i));
    if (h1) h1->Sumw2();
    THnSparse *hSparse = dynamic_cast<THnSparse*>(fOutputListLP->At(i));
    if (hSparse)  hSparse->Sumw2();
  }
  for (Int_t i=0; i < fOutputListHadron->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutputListHadron->At(i));
    if (h1) h1->Sumw2();
    THnSparse *hSparse = dynamic_cast<THnSparse*>(fOutputListHadron->At(i));
    if (hSparse)  hSparse->Sumw2();
  }
  for (Int_t i=0; i < fOutputListQA->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutputListQA->At(i));
    if (h1) h1->Sumw2();
    THnSparse *hSparse = dynamic_cast<THnSparse*>(fOutputListQA->At(i));
    if (hSparse)  hSparse->Sumw2();
  }
     





  fEventCuts.AddQAplotsToList(fOutputListQA, kTRUE);
 
  PostData(1,fOutputList);
  PostData(2,fOutputListMain);
  PostData(3,fOutputListLP);
  PostData(4,fOutputListHadron);
  PostData(5,fOutputListQA);

}

//________________________________________________________________________
void AliAnalysisTaskHaHFECorrel::Terminate(Option_t *)
{
  // fPoolMgr->ClearPools();
  // Info("Terminate");
  AliAnalysisTaskSE::Terminate();
}

//_________________________________________
Double_t AliAnalysisTaskHaHFECorrel::GetDeltaPhi(Double_t phiA,Double_t phiB) const
{
  //Get phi_tag-psi_assoc
  Double_t dPhi = phiA - phiB;
  if (dPhi < -0.5*TMath::Pi()) dPhi = dPhi + TMath::TwoPi();
  if (dPhi >  1.5*TMath::Pi()) dPhi = dPhi - TMath::TwoPi();
  return dPhi;
}


//_________________________________________
Double_t AliAnalysisTaskHaHFECorrel::GetDeltaEta(Double_t etaA,Double_t etaB) const
{
  Double_t dEta=etaB-etaA;  
  return dEta;
}

//_________________________________________
AliVTrack*  AliAnalysisTaskHaHFECorrel::FindLPAndHFE( TObjArray* RedTracks, const AliVVertex *pVtx, Int_t nMother, Int_t listMother[], Double_t mult, Bool_t & EvContTP, Bool_t & EvContNTP, Double_t EventWeight)
{
  AliVTrack* LPtrack=0;
  AliVTrack* LPRecChecktrack=0;
  fLParticle=kFALSE;
  Double_t fLParticleRecCheck=kFALSE;
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();
  
  // Leading Particle
  Double_t ptH=-999; // tmp highest pt particle in acceptance
  Double_t ptHRecCheck = -999; // tmp highest pt patricle in acceptance with random rejection base on rec eff

  // Check NHadron Correlation
  Double_t NHadrons=0, NElectrons=0, NNonElectrons=0, NHadCont=0;
  Double_t NPhotElectronsUntagged=0, NPhotElectronsTagged=0;

  // Loop over all tracks to find LP and HFE
  for(Int_t jTracks = 0; jTracks < ntracks; jTracks++){
          
    AliVParticle* VHtrack = 0x0;
    if(!fUseTender) VHtrack  = fVevent->GetTrack(jTracks);
    if(fUseTender) VHtrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jTracks));
    if (!VHtrack) {
      printf("ERROR: Could not receive associated track %d\n", jTracks);
      continue;
    }
    AliVTrack *Vtrack = dynamic_cast<AliVTrack*>(VHtrack);
    if(!Vtrack) continue;
    AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(VHtrack);
    if (fIsAOD && AODtrack==0) continue;

    Double_t p=-9.,pt=-9.,eta =-9.,phi=-9., recEffH=-9., recEffE=-9.;
    pt = Vtrack->Pt();
    p = Vtrack->P();
    phi = Vtrack->Phi();
    eta = Vtrack->Eta();
    fTrkpt->Fill(pt,0., EventWeight);

    // track cuts
    Bool_t passHadTrackCut=kFALSE;
    passHadTrackCut = ChargedHadronTrackCuts(pVtx, Vtrack, nMother, listMother, EventWeight, kTRUE); // last argument = fillHist
   
    Bool_t passHadPIDCut=kFALSE;
    passHadPIDCut = ChargedHadronPIDCuts(Vtrack,EventWeight); // currently empty
    
    if (passHadTrackCut && passHadPIDCut) { // calculate Hadron reconstruction efficiency
      recEffH = GetHadronRecEff(fVevent->GetRunNumber(), pt, phi, eta, pVtx->GetZ());
      if (recEffH>0) {
	NHadrons+=(1./recEffH);
	//	if (fRecEff) fRecHadronEtaWRecEff->Fill(pt, eta, EventWeight/recEffH);
      }
    }


    // find hadron with the largest pT -> leading particle
    if (passHadTrackCut && passHadPIDCut) {

      // check for LP
      if(Vtrack->Pt()>ptH) {
	ptH = Vtrack->Pt();
	LPtrack=Vtrack;
	fLParticle=kTRUE;
      }
      // check for LP reconstruction (data driven)
      TRandom2 RandGen(1);
      if ( (RandGen.Rndm())>recEffH) {
	if (Vtrack->Pt()>ptHRecCheck) {
	  ptHRecCheck = Vtrack->Pt();
	  LPRecChecktrack=Vtrack;
	  fLParticleRecCheck=kTRUE;
	}
      }


      if (fOneTimeCheck) fEtaVtxZ->Fill(eta, pVtx->GetZ(), EventWeight);
     
      // fill rec hadron histograms for RecEff 
      Double_t fillSparse[4]={-999,-999,-999,-999};
      fillSparse[0]=pt;
      fillSparse[1]=eta;
      fillSparse[2]=phi;
      fillSparse[3]=pVtx->GetZ();
      if (recEffH>0) fRecHadPtEtaPhiVtxWRecEff->Fill(fillSparse, EventWeight/recEffH);

      if (fRecEff && fIsMC && fIsAOD) {
	fRecHadPtEtaPhiVtx->Fill(fillSparse, EventWeight);
	Int_t MClabel=AODtrack->GetLabel();
       	AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(MClabel));
	fRecHadMCSecondaryCont->Fill(pt, MCParticle->IsPhysicalPrimary(), EventWeight);
	                                           
	if (fOneTimeCheck) {  
	  Double_t mcVtx[3];
	  fMCheader->GetVertex(mcVtx);
	  fillSparse[0]=MCParticle->Pt();
	  fillSparse[1]=MCParticle->Eta();
	  fillSparse[2]=MCParticle->Phi();
	  fillSparse[3]=mcVtx[2];
	  fRecHadMCPtEtaPhiVtx->Fill(fillSparse, EventWeight);

	  // controll plot for MC vs Rec pt;
	  fCheckMCPtvsRecPtHad->Fill(MCParticle->Pt(), pt, EventWeight);
	  fCheckMCEtavsRecEtaHad->Fill(MCParticle->Eta(), eta, EventWeight);
	  fCheckMCPhivsRecPhiHad->Fill(MCParticle->Phi(), phi, EventWeight);
	}
      }
    }  
 

    // Select HFE candidates
    Bool_t passHFETrackCut=kFALSE;
    passHFETrackCut= InclElecTrackCuts(pVtx,Vtrack,nMother,listMother, EventWeight, kTRUE); // last argument = fill hists
    if (passHFETrackCut)  fTrkpt->Fill(pt,1, EventWeight); // after track cuts
 
    Bool_t passHFEPIDCut=kFALSE;
    if (passHFETrackCut) passHFEPIDCut= InclElecPIDCuts(Vtrack,  EventWeight, kTRUE); // second argument = fill hists
    if (passHFETrackCut && passHFEPIDCut) {
      recEffE = GetElectronRecEff(fVevent->GetRunNumber(), pt, phi, eta, pVtx->GetZ());
      fTrkpt->Fill(pt,2, EventWeight);
    }
    
    Int_t lsPartner=0, ulsPartner=0;
    Int_t lsPartnerID[20], ulsPartnerID[20];
    Float_t lsPartnerWeight[20], ulsPartnerWeight[20];
    Double_t recTruePairMass=-999.;
    Bool_t trueULSPartner = kFALSE;
    Bool_t isPhotonic = kFALSE;
    Bool_t isHadron = kFALSE;
    Float_t RecPartnerPt=-999., MCPartnerPt=-999.;
    if (passHFETrackCut && passHFEPIDCut && recEffE>0) { // if HFE is found, look for ls and uls partner
      NElectrons+=(1./recEffE);
      FindPhotonicPartner(jTracks, Vtrack, pVtx, nMother, listMother, lsPartner, ulsPartner, lsPartnerID, ulsPartnerID, lsPartnerWeight, ulsPartnerWeight, trueULSPartner, isPhotonic, MCPartnerPt, RecPartnerPt, EventWeight, mult);
      if (fIsMC) CheckPhotonicPartner(Vtrack, trueULSPartner, MCPartnerPt, RecPartnerPt, EventWeight);
      if (fIsMC && fTagEff) {
	if (isPhotonic) {
	  if (trueULSPartner) {
	    NPhotElectronsTagged+=(1./recEffE);
	    // fTagEtaPt1Pt2->Fill(Vtrack->Eta(), Vtrack->Pt(), RecPartnerPt, EventWeight); // EvW not req but consistent
	    // fTagEtaPhiPt->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt());
	    // fTagEtaZvtxPt->Fill(Vtrack->Eta(), pVtx->GetZ(), Vtrack->Pt());
	    // fTagEtaPhiPtwW->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt(), 1./recEffE);
	    // fTagEtaZvtxPtwW->Fill(Vtrack->Eta(), pVtx->GetZ(), Vtrack->Pt(), 1./recEffE);
	    if (Vtrack->Pt()>1.0) EvContTP=kTRUE;
	  }
	  else  {
	    NPhotElectronsUntagged+=(1./recEffE);
	    // fNonTagEtaPt1Pt2->Fill(Vtrack->Eta(), Vtrack->Pt(), RecPartnerPt, EventWeight); // EvW not req but consistent
	    // fNonTagEtaPhiPt->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt());
	    // fNonTagEtaZvtxPt->Fill(Vtrack->Eta(), pVtx->GetZ(), Vtrack->Pt());
	    // fNonTagEtaPhiPtwW->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt(), 1./recEffE);
	    // fNonTagEtaZvtxPtwW->Fill(Vtrack->Eta(), pVtx->GetZ(), Vtrack->Pt(), 1./recEffE);
	    if (Vtrack->Pt()>1.0) EvContNTP=kTRUE;
	  }
	}
	EvaluateTaggingEfficiency(Vtrack, lsPartner, ulsPartner, trueULSPartner, EventWeight, mult, recEffE);
      }
      else if (fTagEff) {
	/*
	fTagEtaPhiPt->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt(), ulsPartner);
	fTagEtaPhiPt->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt(), -lsPartner);
	fTagEtaZvtxPt->Fill(Vtrack->Eta(), pVtx->GetZ(), Vtrack->Pt(), ulsPartner);
	fTagEtaZvtxPt->Fill(Vtrack->Eta(), pVtx->GetZ(), Vtrack->Pt(), -lsPartner);
	fTagEtaPhiPtwW->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt(), ulsPartner/recEffE);
	fTagEtaPhiPtwW->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt(), -lsPartner/recEffE);
	fTagEtaZvtxPtwW->Fill(Vtrack->Eta(), pVtx->GetZ(), Vtrack->Pt(), ulsPartner/recEffE);
	fTagEtaZvtxPtwW->Fill(Vtrack->Eta(), pVtx->GetZ(), Vtrack->Pt(), -lsPartner/recEffE);
	*/

      }
      fInclElecPtEta->Fill(pt, eta, mult,  EventWeight);
      fInclElecPtEtaWRecEff->Fill(pt, eta, mult, EventWeight/recEffE);
      //      fInclElecP->Fill(p);
      if (fOneTimeCheck) {
	fInclElecPhi->Fill(pt,phi, EventWeight); 
	fInclElecPhiEta->Fill(eta, phi, EventWeight);
      }
      
      // count reconstructed heavy electrons for recEff
      if (fIsMC & fIsAOD) {
	Int_t MClabel=AODtrack->GetLabel();
	AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(MClabel));     
	Int_t PDGCode = abs(MCParticle->GetPdgCode());
	Double_t mcVtx[3];
	fMCheader->GetVertex(mcVtx);

	if (PDGCode ==11) {
	  AliAODMCParticle* mcPartMother=(AliAODMCParticle*)fMC->GetTrack(MCParticle->GetMother());
	  Int_t mcMotherPDG = abs(mcPartMother->GetPdgCode());

	  if (mcMotherPDG==11 || mcMotherPDG==15) { // to includ HF->e->e HF->tau->e
	    AliAODMCParticle* mcPartGMother=(AliAODMCParticle*)fMC->GetTrack(mcPartMother->GetMother());
	    mcMotherPDG = abs(mcPartGMother->GetPdgCode());
	  }

	  
	  Int_t MIsHeavy  = Int_t (mcMotherPDG / TMath::Power(10, Int_t(TMath::Log10(mcMotherPDG))));
	  fRecElecMCSecondaryCont->Fill(pt, MCParticle->IsPhysicalPrimary(), EventWeight);

	  if (MIsHeavy>3 && MIsHeavy<6) {
	    // if (MCParticle->IsPhysicalPrimary()) cout << "IsPrimary" << endl;
	    // else cout << "Is not primary" << endl;

	    
	    Double_t fillSparse[4];
	    fillSparse[0]=pt;
	    fillSparse[1]=eta;
	    fillSparse[2]=phi;
	    fillSparse[3]=pVtx->GetZ();
	    fRecElecPtEtaPhiVtxWRecEff->Fill(fillSparse, EventWeight/recEffE); // for checks if right efficiency has been used
	    
	    //fRecHFE->Fill(pt, EventWeight);
	    // fRecHFEEtaWRecEff->Fill(pt, eta, EventWeight/recEffE);
	    if (fRecEff) {
	      // Fill Reconstructed Histogram
	      fRecElecPtEtaPhiVtx->Fill(fillSparse, EventWeight);
	      fillSparse[0]=MCParticle->Pt();
	      fillSparse[1]=MCParticle->Eta();
	      fillSparse[2]=MCParticle->Phi();
	      fillSparse[3]=mcVtx[2];
	      if (fOneTimeCheck) {
		fRecElecMCPtEtaPhiVtx->Fill(fillSparse, EventWeight);
		fCheckMCPtvsRecPtEle->Fill(MCParticle->Pt(), pt, EventWeight);
	      }
	    }
	  }
	}

	if (PDGCode != 11) {
	  NHadCont+=(1./recEffE);
	  isHadron=kTRUE;
	  if (fHadCont && fOneTimeCheck) fHadContEtaPhiPt->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt(), EventWeight);
	}
      }
        
      // store only essential informations of these electrons for later correlations and mixed event
      CloneAndReduceTrackList(RedTracks, Vtrack, lsPartner, ulsPartner, lsPartnerID, ulsPartnerID,  lsPartnerWeight,  ulsPartnerWeight, trueULSPartner, MCPartnerPt, RecPartnerPt, isPhotonic, isHadron);
    }

    Bool_t passNonElecPIDCut=kFALSE;
    passNonElecPIDCut=AssoHadronPIDCuts(Vtrack, EventWeight);
    if (passHFETrackCut && passNonElecPIDCut && recEffE>0) {
      NNonElectrons+=(1./recEffE);
      if (fOneTimeCheck) fHadContTPCEtaPhiPt->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt(),EventWeight);
    }
  }

  // CheckNoOfElectrons/Event
  Double_t fillSparse[5];
  fillSparse[0]=NHadrons;
  fillSparse[1]=NElectrons;
  fillSparse[2]=NNonElectrons;
  fillSparse[4]=mult;
  if (fIsMC) fillSparse[3] = NHadCont;
  else fillSparse[3]=0;
  //  fCheckNHadronScaling->Fill(fillSparse);

  if (fIsMC) { 
    fillSparse[2]= NPhotElectronsTagged;
    fillSparse[3]= NPhotElectronsUntagged;
    //fCheckNPhotHadScaling->Fill(fillSparse);
  }

  // if no leading particle is found, mark this in flag fLParticle;
  if (!LPtrack) {
    fLParticle=kFALSE;
  }
  else if (fIsMC && fIsAOD && RedTracks->GetEntriesFast()>0) {
    AliAODTrack *LPAODtrack = dynamic_cast<AliAODTrack*>(LPtrack);   
    if (!LPAODtrack) {
      fLParticle=kFALSE;
      return 0;
    }
    Int_t MClabel=LPAODtrack->GetLabel();
    AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(MClabel));     
    Int_t PDGCode = abs(MCParticle->GetPdgCode());
    Int_t PDGCodeMother=0;
    if (MCParticle->GetMother()>=0) {
      AliAODMCParticle* MCParticleMother = (AliAODMCParticle*) fMC->GetTrack(MCParticle->GetMother());
      PDGCodeMother = abs(MCParticleMother->GetPdgCode());
    }
    else {
      PDGCodeMother=-1;
    }
    Double_t fillSparse[3];
    fillSparse[0]=LPAODtrack->Pt();
    fillSparse[1]=PDGMap.find(PDGCode)->second;
    fillSparse[2]=PDGMap.find(PDGCodeMother)->second;
    if (fCorrLParticle)  {
      fMCLeadingParticle->Fill(fillSparse, EventWeight);

	if (LPRecChecktrack) {
	  AliAODTrack *LPRecCheckTrack = dynamic_cast<AliAODTrack*>(LPRecChecktrack);   
	  Double_t LPDeltaPhi = GetDeltaPhi(LPAODtrack->Phi(), LPRecCheckTrack->Phi());
	  Double_t LPDeltaEta = GetDeltaEta(LPAODtrack->Eta(), LPRecCheckTrack->Eta());
	  Double_t LPDeltaPt  = LPAODtrack->Pt()  -LPRecCheckTrack->Pt();
	  fCompareLPRecCheck->Fill(LPDeltaPt, LPDeltaPhi, LPDeltaEta, EventWeight); // EvW not req but consistent
	}
    }




  }
  //  if (fLParticle && RedTracks->GetEntriesFast()>0) cout<< "RecLP " << LPtrack->GetLabel() <<  "\t" << LPtrack->Pt() << endl;

  return LPtrack;
}

//_________________________________________
void AliAnalysisTaskHaHFECorrel::FindPhotonicPartner(Int_t iTracks, AliVTrack* Vtrack, const AliVVertex *pVtx, Int_t nMother, Int_t listMother[], Int_t &lsPartner, Int_t &ulsPartner, Int_t *lsPartnerID, Int_t *ulsPartnerID,  Float_t *lsPartnerWeight, Float_t *ulsPartnerWeight, Bool_t &trueULSPartner, Bool_t &isPhotonic, Float_t &MCPartnerPt, Float_t &RecPartnerPt, Double_t EventWeight, Double_t mult) {
  

  AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);
  if(!AODtrack) return;
  lsPartner=0;
  ulsPartner=0;
  trueULSPartner=kFALSE; // only for MC use
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();



  // tagged particle
  Double_t pt=-9.,phi=-9.,eta =-9.,recEff=-9.;
  Int_t charge = 0;      
  pt = Vtrack->Pt();
  phi = Vtrack->Phi();
  eta = Vtrack->Eta();
  charge = Vtrack->Charge();
  recEff = GetElectronRecEff(fVevent->GetRunNumber(), pt,phi,eta,pVtx->GetZ());
  if (recEff<0) return;
        
    
  // Loop over all tracks to find photonic partner
  for(Int_t jTracks = 0; jTracks < ntracks; jTracks++){
        
    if(jTracks==iTracks) continue;
    
    AliVParticle* Vassotrack = 0x0;
    if(!fUseTender) Vassotrack  = fVevent->GetTrack(jTracks);
    if(fUseTender)  Vassotrack  = dynamic_cast<AliVTrack*>(fTracks_tender->At(jTracks));
    if (!Vassotrack) {
      printf("ERROR: Could not receive associated track %d\n", jTracks);
      continue;
    }
    AliVTrack *VtrackAsso = dynamic_cast<AliVTrack*>(Vassotrack);
    if(!VtrackAsso) continue;
    AliAODTrack *AODtrackAsso = dynamic_cast<AliAODTrack*>(VtrackAsso);   
    if(!AODtrackAsso) continue;
  
  
    // associated particle variables
    Double_t phiAsso=-9.;
    Int_t chargeAsso = 0;
    phiAsso = VtrackAsso->Phi();
    chargeAsso = VtrackAsso->Charge();

    double dphi = -9;

    // Cuts for associated electrons for HFE-HFE correlations
    Bool_t passAssoTrackCutIncl = kFALSE;
    //  passAssoTrackCutIncl = InclElecTrackCuts(pVtx,VtrackAsso,nMother,listMother, EventWeight);

    Bool_t passAssoPIDCutIncl = kFALSE;
    //passAssoPIDCutIncl = InclElecPIDCuts(VtrackAsso, kFALSE, EventWeight);

    if (passAssoTrackCutIncl && passAssoPIDCutIncl) {
      dphi = GetDeltaPhi(phi,phiAsso); // electron-electron dphi
      //  fElecDphi->Fill(pt,dphi);        // HFEHFE
    }

    // looser Track cuts for associated photonic eletron
    Bool_t passAssoTrackCutPhot = kFALSE;
    passAssoTrackCutPhot = PhotElecTrackCuts(pVtx,VtrackAsso,nMother,listMother, EventWeight);
    if(!passAssoTrackCutPhot) continue;
   
    Bool_t passAssoPIDCutPhot = kFALSE;
    passAssoPIDCutPhot = PhotElecPIDCuts(VtrackAsso, EventWeight); 
    if (!passAssoPIDCutPhot) continue;         
   
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
    Double_t openingAngle = -999., mass=999., width = -999;
    Double_t energy=-999., transenergy=-999., mom=-999., tansmom=-999.;
        
    Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
    if(charge>0) fPDGe1 = -11;
    if(chargeAsso>0) fPDGe2 = -11;
        
    if(charge == chargeAsso) fFlagLS = kTRUE;
    if(charge != chargeAsso) fFlagULS = kTRUE;
    
    // Two different Alogrithms in use (check difference) A. KalmanFilter (old), B. DCA
    if (fUseKFforPhotonicPartner) {        
      AliKFParticle::SetField(fVevent->GetMagneticField());
      AliKFParticle ge1(*Vtrack, fPDGe1);
      AliKFParticle ge2(*VtrackAsso, fPDGe2);
      AliKFParticle recg(ge1, ge2);
        
      if(recg.GetNDF()<1) continue;
      Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
      if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
      
      openingAngle = ge1.GetAngle(ge2);
      //if(openingAngle > fOpeningAngleCut) continue;
      
      recg.GetMass(mass,width);     
    }
    else {
      //Variables
      Double_t p1[3], p2[3];
      Double_t xt1, xt2; //radial position track 1 and 2  at the DCA point
      Double_t dca12; //DCA 1-2
      Bool_t hasdcaT1, hasdcaT2;
      Double_t bfield = fVevent->GetMagneticField();
            
      AliExternalTrackParam extTrackParam1;
      extTrackParam1.CopyFromVTrack(Vtrack);
      AliExternalTrackParam extTrackParam2;
      extTrackParam2.CopyFromVTrack(VtrackAsso);

      //DCA track1-track2
      dca12 = extTrackParam2.GetDCA(&extTrackParam1,bfield,xt2,xt1);
                
      //Momento of the track extrapolated to DCA track-track
      //Track1
      hasdcaT1 = extTrackParam1.GetPxPyPzAt(xt1,bfield,p1);
      //Track2
      hasdcaT2 = extTrackParam2.GetPxPyPzAt(xt2,bfield,p2);
      
      if(!hasdcaT1 || !hasdcaT2) AliWarning("There could be a problem in the extrapolation");
      //track1-track2 Invariant Mass
      Double_t eMass = 0.000510998910; //Electron mass in GeV
      Double_t pP1 = sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]); //Track 1 momentum
      Double_t pP2 = sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]); //Track 1 momentum

      TLorentzVector v1(p1[0],p1[1],p1[2],sqrt(eMass*eMass+pP1*pP1));
      TLorentzVector v2(p2[0],p2[1],p2[2],sqrt(eMass*eMass+pP2*pP2));
      mass = (v1+v2).M(); //Invariant Mass
      openingAngle = v1.Angle(v2.Vect()); //Opening Angle (Total Angle)
    }

    if(fFlagLS){
      if (fOneTimeCheck) fOpeningAngleLS->Fill(openingAngle,pt, EventWeight);
      fInvmassLS->Fill(mass,pt, EventWeight);
    }
    if(fFlagULS){
      if (fOneTimeCheck)  fOpeningAngleULS->Fill(openingAngle,pt, EventWeight);
      fInvmassULS->Fill(mass,pt, EventWeight);
    }

    if (fIsAOD && fIsMC) {
      AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(AODtrack->GetLabel()));  
      AliAODMCParticle* MCParticleAsso = (AliAODMCParticle*) fMC->GetTrack(abs(AODtrackAsso->GetLabel())); 
      Int_t Mother=(MCParticle->GetMother());
      AliAODMCParticle* MCMother=(AliAODMCParticle*) fMC->GetTrack(abs(Mother));
      Int_t PdgCode = abs(MCMother->GetPdgCode()); 
      if (abs(MCParticle->GetPdgCode())==11) {
	if (fFlagULS && abs(MCParticleAsso->GetPdgCode())==11) {
	  if ( HaveSameMother(AODtrack->GetLabel(), AODtrackAsso->GetLabel())) {
	    fInvmassMCTrue->Fill(mass,pt,PDGMap.find(PdgCode)->second, EventWeight);
	    fInvmassMCTrue->GetZaxis()->SetBinLabel(PDGMap.find(PdgCode)->second, Form("%i", PdgCode));
	    RecPartnerPt=AODtrackAsso->Pt();
	    MCPartnerPt=MCParticleAsso->Pt();
	    if (IsPhotonicElectron(AODtrack->GetLabel())&& mass<fInvmassCut){
	      fPhotPt2MCRec->Fill(AODtrackAsso->Pt(), MCParticleAsso->Pt(), EventWeight);
	    }
	  }
	}
      }    
    }
   
        
    if(mass<fInvmassCut && fFlagULS) {
      if (fIsAOD && fIsMC && !trueULSPartner) {
	trueULSPartner = HaveSameMother(AODtrack->GetLabel(), AODtrackAsso->GetLabel());
	if (trueULSPartner) {
	  //  cout << "TruePartnerRecPt" << AODtrackAsso->Pt() << endl;
	  //	  cout << "RecPt " << RecPartnerPt << "\tMCPt " << MCPartnerPt << endl;
	}
      }
      /*
      Int_t nclustersOuter=0;
      for(Int_t ily = 2; ily < 6; ily++) {
	if(AODtrackAsso->HasPointOnITSLayer(ily)) nclustersOuter++;
      }
      if (!trueULSPartner && AODtrackAsso->TestFilterMask(AliAODTrack::kTrkITSsa)) cout << "tagged, but not true:" << nclustersOuter << endl;
      if (trueULSPartner && nclustersOuter<2) cout << "true tagged:" << nclustersOuter << endl;
      */



      if (fOneTimeCheck) fULSElecPhi->Fill(pt,phi, EventWeight);
      fULSElecPt->Fill(pt, mult, EventWeight);
      fULSElecPtWRecEff->Fill(pt, mult,EventWeight/recEff);
      dphi = GetDeltaPhi(phi,phiAsso); // electron-electron dphi
      // fULSElecDphi->Fill(pt,dphi);
      ulsPartnerID[ulsPartner]=VtrackAsso->GetID();
      ulsPartnerWeight[ulsPartner]=GetNonTagCorr(pt, VtrackAsso->Pt());
      ulsPartner++;
      // if (fIsAOD && fIsMC && trueULSPartner) recTruePairMass=mass;
      if (fIsAOD)   fPhotPt1Pt2Only->Fill(AODtrack->Pt(), AODtrackAsso->Pt(), EventWeight);
    }
        
    if(mass<fInvmassCut && fFlagLS){
      if (fOneTimeCheck) fLSElecPhi->Fill(pt,phi, EventWeight);
      fLSElecPt->Fill(pt,mult, EventWeight);
      fLSElecPtWRecEff->Fill(pt, mult, EventWeight/recEff);
      dphi = GetDeltaPhi(phi,phiAsso); // electron-electron dphi
      // fLSElecDphi->Fill(pt,dphi);
      lsPartnerID[lsPartner]=VtrackAsso->GetID();
      lsPartnerWeight[lsPartner]=GetNonTagCorr(pt, VtrackAsso->Pt());
      lsPartner++;
      if (fIsAOD)   fPhotPt1Pt2Only->Fill(AODtrack->Pt(), AODtrackAsso->Pt(), -EventWeight);
    }
  
  }//track loop
  if (fIsMC) { 
    isPhotonic=IsPhotonicElectron(AODtrack->GetLabel());
  } 
}



void AliAnalysisTaskHaHFECorrel::CheckPhotonicPartner(AliVTrack* Vtrack, Bool_t Tagged, Float_t& MCPartnerPt, Float_t RecPartnerPt, Double_t EventWeight) {


  AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);
  if(!AODtrack) return;
  Bool_t isPhotonic=IsPhotonicElectron(AODtrack->GetLabel());
  AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(AODtrack->GetLabel()));  
  if (abs(MCParticle->GetPdgCode())!=11 || !isPhotonic) return;
  Int_t Mother=(MCParticle->GetMother());
 
  Bool_t foundPartner=kFALSE;
  Int_t motherPartner=-999;
  Double_t mass = -999., energy=-999., transenergy=-999., mom=-999., transmom=-999., ptPartner=-999.;
  for (Int_t i=0; i<fMC->GetNumberOfTracks(); i++) {
    if (i==abs(AODtrack->GetLabel())) continue;
    AliAODMCParticle *mcPart  = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(i));
    if (abs(mcPart->GetPdgCode())!=11) continue;
    if (mcPart->GetMother()!=Mother) continue;
    // cout << "partner " << i << endl;
    if (foundPartner) cout << "p0->eeee? " << endl; // but in p0->eeee
    Double_t p1[3], p2[3];
    MCParticle->PxPyPz(p1);
    mcPart->PxPyPz(p2);
    TLorentzVector v1(p1[0],p1[1],p1[2], MCParticle->E());
    TLorentzVector v2(p2[0],p2[1],p2[2], mcPart->E());
    TLorentzVector Sum=v1+v2;
    mass = Sum.M(); //Invariant Mass
    energy = Sum.E();
    transenergy = Sum.Et();
    mom = Sum.P();
    transmom = Sum.Pt();
    // cout << "InvMass " << mass << "\tEnergy " << energy << "\tEt " << transenergy << "\tmom " << mom << "\t transmom " << transmom << endl;
    ptPartner=mcPart->Pt();
    motherPartner=mcPart->GetMother();
    foundPartner=kTRUE;
  }
  
  if (!foundPartner) {  return; }
  if (Tagged && TMath::Abs(ptPartner-MCPartnerPt)>0.001) {
    cout << "Error in Partner pt" << endl;
    cout << "ptPartner " << ptPartner << "\tMCPartnerPt " << MCPartnerPt << " \t diff " << TMath::Abs(ptPartner-MCPartnerPt) << endl;
  }
  else MCPartnerPt=ptPartner;
  // compare invariant mass
  //  if (recTruePairMass>0) fRecMCInvMass->Fill(MCParticle->Pt(), mass, recTruePairMass);

  AliAODMCParticle *MCMother  = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(abs(Mother)));
  Int_t PdgCodeMother = abs(MCMother->GetPdgCode());
  Double_t MotherPt=  MCMother->Pt();;
  Int_t GMother = abs(MCMother->GetMother());
  AliAODMCParticle *MCGMother;
  Int_t PdgCodeGrandMother = -999;
  Int_t PdgCodeGGM = -999;
  Double_t PtGrandMother=-999.;
  if (GMother>=0) {
    AliAODMCParticle* MCGMother = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(abs(MCMother->GetMother())));
    PdgCodeGrandMother = abs(MCGMother->GetPdgCode());
    PtGrandMother=MCGMother->Pt();
    AliAODMCParticle* MCGGMother = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(abs(MCGMother->GetMother())));
    PdgCodeGGM = abs(MCGGMother->GetPdgCode());
  }

  Double_t factor = -999;
  if (Tagged) factor = 1.;
  else factor= -1.;


  Double_t PtMotherWeight=1.;
  if (PdgCodeMother==22) {
    if (PdgCodeGrandMother == 221) {
      PtMotherWeight=GetEtaWeight(PtGrandMother); 
    }
    else if (PdgCodeGrandMother == 111) {
      PtMotherWeight=GetPionWeight(PtGrandMother); 
    }
  }
  else  {
    if (PdgCodeMother == 221) {
      PtMotherWeight=GetEtaWeight(MotherPt); 
    }
    else if (PdgCodeMother == 111) {
      PtMotherWeight=GetPionWeight(MotherPt); 
    }
  }
    
	    
  Double_t fillSparse[] = {AODtrack->Pt(), factor*ptPartner, MCMother->Pt(), 1., 1.};

  if (abs(PdgCodeMother)==22) {
    if (PdgCodeGrandMother == 111) fillSparse[3]=0;
    else if (PdgCodeGrandMother ==221) fillSparse[3]=1;
    if (PdgCodeGGM>100) fillSparse[4]=1;
    else fillSparse[4]=0;
  }
  else{
    if (abs(PdgCodeMother)==111)   fillSparse[3]=2.;
    else if (abs(PdgCodeMother)==221) fillSparse[3]=3.;
    if (PdgCodeGrandMother>100) fillSparse[4]=1;
    else fillSparse[4]=0;
  }
  if (factor>0) fPhotPt1PtMTag->Fill(MCParticle->Pt(), MCMother->Pt(), EventWeight);
  else fPhotPt1PtMNTag->Fill(MCParticle->Pt(), MCMother->Pt(), EventWeight);

  fPhotPt1Pt2->Fill(fillSparse, PtMotherWeight*EventWeight);
  if (factor>0) fPhotPt1Pt2Corr->Fill(fillSparse, EventWeight*PtMotherWeight*GetNonTagCorr(AODtrack->Pt(), RecPartnerPt)); 
  else  fPhotPt1Pt2Corr->Fill(fillSparse, EventWeight*PtMotherWeight);
  fPhotPt1Pt2MC->Fill(fillSparse, EventWeight);

  // fillSparse[1]=factor*RecPartnerPt;
  // if (factor<0) fillSparse[1]=factor*ptPartner;
  //fPhotPt1Pt2Rec->Fill(fillSparse, PtMotherWeight);
 
  // only to check the invariant mass distribution
  // fillSparse[1]=mass*factor;
  /*
  fPhotPt1E->Fill(fillSparse, PtMotherWeight);
  if (factor>0)   fPhotPt1ECorr->Fill(fillSparse, PtMotherWeight*GetNonTagCorr(AODtrack->Pt(), RecPartnerPt));
  else   fPhotPt1ECorr->Fill(fillSparse, PtMotherWeight);
  */
  // fillSparse[2]=ptPartner;
  // fPhotPt1Pt2E->Fill(fillSparse, PtMotherWeight*EventWeight);
  // if (factor>0)   fPhotPt1Pt2ECorr->Fill(fillSparse, EventWeight*PtMotherWeight*GetNonTagCorr(AODtrack->Pt(), RecPartnerPt));
  // else   fPhotPt1Pt2ECorr->Fill(fillSparse, EventWeight*PtMotherWeight);
  // 


  return;



}
  

void AliAnalysisTaskHaHFECorrel::CorrelateElectron(TObjArray* RedTracksHFE) {

  // work in progress - only correlate electrons which fullfill inclusive cuts

  for (Int_t j=0; j<RedTracksHFE->GetEntriesFast(); j++) {

    AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(j);

    // tagged particle
    Double_t pt=-9.,eta =-9.,phi=-9.;
    Int_t charge = 0;
    pt = RedTrack->Pt();
    phi = RedTrack->Phi();
    eta = RedTrack->Eta();
    charge = RedTrack->Charge();

    /*
    // fill trigger electron
    fElecTrigger->Fill(pt);

    for (Int_t k=0; k<RedTracksHFE->GetEntriesFast(); k++) {
      if (k==j) continue;
      AliBasicParticleHaHFE *trackAsso = (AliBasicParticleHaHFE*) RedTracksHFE->At(k);
      Double_t ptAsso=-9.,etaAsso =-9.,phiAsso=-9.;
      Int_t lsAsso=0, ulsAsso=0;
        
      ptAsso = trackAsso->Pt();
      phiAsso = trackAsso->Phi();
      etaAsso = trackAsso->Eta();
      lsAsso = trackAsso->LS();
      ulsAsso = trackAsso->ULS();
      
      // Fill Incl-Incl (already done above)
      Double_t dphi = -99;
      dphi = GetDeltaPhi(phi, phiAsso);
    
      for (Int_t l=0; l<lsAsso;l++)   fLSElecDphiDiffMethod->Fill(pt, dphi);
      for (Int_t l=0; l<ulsAsso; l++) fULSElecDphiDiffMethod->Fill(pt, dphi);
      
    }
    */
  }
}


/// Maybe include LP case
void AliAnalysisTaskHaHFECorrel::CorrelateHadron(TObjArray* RedTracksHFE,  const AliVVertex* pVtx, Int_t nMother, Int_t listMother[], Float_t mult, Float_t maxPt, Double_t EventWeight) {
    
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();

  Bool_t *HadronTrigger=new Bool_t[fAssPtElec_Nbins];
  Bool_t **ElectronIsTrigger = new Bool_t*[RedTracksHFE->GetEntriesFast()];
  Bool_t **ElectronIsTriggerNoP = new Bool_t*[RedTracksHFE->GetEntriesFast()];
  Bool_t **ElectronIsTriggerNoPTrue = new Bool_t*[RedTracksHFE->GetEntriesFast()];
  Bool_t **HadContIsTrigger=new Bool_t*[RedTracksHFE->GetEntriesFast()];
  Bool_t **PhotElecWPartnerTrigger=new Bool_t*[RedTracksHFE->GetEntriesFast()];
  Bool_t **PhotElecWoPartnerTrigger=new Bool_t*[RedTracksHFE->GetEntriesFast()];
  for (Int_t j=0; j<RedTracksHFE->GetEntriesFast(); j++)
    {
      ElectronIsTrigger[j]=new Bool_t[fAssPtHad_Nbins];
      ElectronIsTriggerNoP[j]=new Bool_t[fAssPtHad_Nbins];
      ElectronIsTriggerNoPTrue[j] = new Bool_t[fAssPtHad_Nbins];
      HadContIsTrigger[j]=new Bool_t[fAssPtHad_Nbins];
      PhotElecWPartnerTrigger[j] = new Bool_t[fAssPtHad_Nbins];
      PhotElecWoPartnerTrigger[j] = new Bool_t[fAssPtHad_Nbins];
      for (Int_t m=0; m<fAssPtHad_Nbins; m++) {
	ElectronIsTrigger[j][m]=kFALSE;
	ElectronIsTriggerNoP[j][m]=kFALSE;
	ElectronIsTriggerNoPTrue[j][m]=kFALSE;
	HadContIsTrigger[j][m]=kFALSE;
	PhotElecWPartnerTrigger[j][m]=kFALSE;
	PhotElecWoPartnerTrigger[j][m]=kFALSE;
      }
    }
  Double_t *ElectronIsTriggerPt = new Double_t[RedTracksHFE->GetEntriesFast()];
  Double_t *MotherWeight = new Double_t[RedTracksHFE->GetEntriesFast()];
  Bool_t **NonElectronIsTrigger = new Bool_t*[ntracks];
  for (Int_t j=0; j<ntracks; j++) { 
    NonElectronIsTrigger[j]=new Bool_t[fAssPtHad_Nbins]; 
    for (Int_t m=0; m<fAssPtHad_Nbins; m++) {
      NonElectronIsTrigger[j][m]=kFALSE;
    }
  }
  Double_t *NonElectronIsTriggerPt= new Double_t[ntracks];
  Double_t *NonElectronIsTriggerWeight = new Double_t[ntracks];
  

  // ------------------------  // Track loop for hadron correlations
   for(Int_t iTracks = 0; iTracks < ntracks; iTracks++) {   
    for (Int_t l=0; l<fAssPtElec_Nbins; l++) HadronTrigger[l]=kFALSE;

    AliVParticle* Vtrack = 0x0;
    if(!fUseTender) Vtrack  = fVevent->GetTrack(iTracks);
    if(fUseTender)  Vtrack  = dynamic_cast<AliVTrack*>(fTracks_tender->At(iTracks));
    if (!Vtrack) {
      printf("ERROR: Could not receive tagged  track %d\n", iTracks);
      continue;
    }
    AliVTrack   *track = dynamic_cast<AliVTrack*>(Vtrack);
    if (!track) continue;
    AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);
    AliESDtrack *ESDtrack = dynamic_cast<AliESDtrack*>(Vtrack);
    if (fIsAOD) {  if(!AODtrack) continue; }
    else {if (!ESDtrack) continue;}
    
    // track cuts
    Bool_t passTrackCut=kFALSE;
    passTrackCut = ChargedHadronTrackCuts(pVtx, track, nMother, listMother, EventWeight);
    if (!passTrackCut) continue;

    // pid cuts
    Bool_t passPIDCut=kFALSE;
    passPIDCut = ChargedHadronPIDCuts(track, EventWeight);
    if (!passPIDCut) continue;
   
    // associated particle
    Double_t ptH=-9.,etaH =-9.,phiH=-9.,recEffH=-9.;
    Int_t idH=-9;
    ptH = track->Pt();
    phiH = track->Phi();
    etaH = track->Eta();
    idH = track->GetID();
    recEffH = GetHadronRecEff(fVevent->GetRunNumber(), ptH,phiH,etaH,pVtx->GetZ());
    if (recEffH<0) continue;
 
    // Correlate hadrons with Hadron (kTRUE, kFALSE);
    CorrelateWithHadrons(track, pVtx, nMother, listMother, kTRUE, kFALSE, NonElectronIsTrigger, NonElectronIsTriggerPt, NonElectronIsTriggerWeight, RedTracksHFE->GetEntriesFast(), EventWeight); 

    // loop over all electrons
    for (Int_t k=0; k<RedTracksHFE->GetEntriesFast(); k++) {

      AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(k);

      if (RedTrack->ID()==idH) continue;
      // tagged particle
      Double_t pt=-9.,eta =-9.,phi=-9., recEffE=-9.;
      Int_t charge = 0, ls=0, uls=0;    
      pt = RedTrack->Pt();
      phi = RedTrack->Phi();
      eta = RedTrack->Eta();
      charge = RedTrack->Charge();
      ls = RedTrack->LS();
      uls = RedTrack->ULS();
      // cout << "ls " <<  ls << "\t" << uls << endl;
      recEffE = GetElectronRecEff(fVevent->GetRunNumber(), pt,phi,eta,pVtx->GetZ());
      if (recEffE<0) continue;

      Double_t dphi = -99, deta=-99;
      dphi = GetDeltaPhi(phi, phiH);
      deta = GetDeltaEta(eta, etaH);

      Double_t fillSparse[5]={-999,-999,-999,-999,-999};
      fillSparse[0]=ptH;
      fillSparse[1]=pt;
      fillSparse[2]=dphi;
      fillSparse[3]=deta;
      fillSparse[4]=pVtx->GetZ();
      
      CheckElectronIsTrigger(ptH, ElectronIsTrigger[k]);
      ElectronIsTriggerPt[k]=pt;
      CheckHadronIsTrigger(pt, HadronTrigger);  
      fInclElecHa->Fill(fillSparse, EventWeight/(recEffH*recEffE)); 
      //   fLSElecHa->Fill(fillSparse, ls/(recEffH*recEffE)); // not needed anymore
      for (int j=0; j<uls; j++) fULSElecHa->Fill(fillSparse, EventWeight/(recEffH*recEffE)); // TriggerConditionProvedAtnd
      for (int j=0; j<ls; j++)  fULSElecHa->Fill(fillSparse, -1.*EventWeight/(recEffH*recEffE)); // TriggerConditionProvedAtnd
      if (RedTrack->IsPhotonic()) {
	for (int j=0; j<uls; j++) fULSElecHaTrue->Fill(fillSparse, EventWeight/(recEffH*recEffE));
	for (int j=0; j<ls; j++) fULSElecHaTrue->Fill(fillSparse, -1.*EventWeight/(recEffH*recEffE));
      }

      
      // Fill Histograms with missing partner
      Bool_t HadIsULSPartner=kFALSE;
      Bool_t HadIsLSPartner=kFALSE;
      for (int j=0; j<uls; j++) {
	if (idH==RedTrack->ULSPartner(j)) {
	  HadIsULSPartner =kTRUE;
	}
      }
      for (int j=0; j<ls; j++) {
	if (idH==RedTrack->LSPartner(j)) {
	  HadIsLSPartner=kTRUE;
	}
      }
      


      if (!HadIsULSPartner && !HadIsLSPartner) {
	CheckElectronIsTrigger(ptH, ElectronIsTriggerNoP[k]);
	
	if (uls>0 || ls>0) {
	  for (int j=0; j<uls; j++) fElecHaULSNoPartner->Fill(fillSparse, EventWeight/(recEffH*recEffE));
	  for (int j=0; j<ls; j++) fElecHaULSNoPartner->Fill(fillSparse, -1.*EventWeight/(recEffH*recEffE));
	  
	  for (int j=0; j<ls; j++) {
	    fElecHaULSNoPartnerCorr->Fill(fillSparse, -RedTrack->LSPartnerWeight(j)*EventWeight/(recEffH*recEffE));
	  }
	  for (int j=0; j<uls; j++) {
	    fElecHaULSNoPartnerCorr->Fill(fillSparse, RedTrack->ULSPartnerWeight(j)*EventWeight/(recEffH*recEffE));
	  }

	  if (RedTrack->IsPhotonic()) {
	    CheckElectronIsTrigger(ptH, ElectronIsTriggerNoPTrue[k]);
	    for (int j=0; j<uls; j++) fElecHaULSNoPartnerCorrTrue->Fill(fillSparse, EventWeight/(recEffH*recEffE));
	    for (int j=0; j<ls; j++) fElecHaULSNoPartnerCorrTrue->Fill(fillSparse, -1.*EventWeight/(recEffH*recEffE));
	  }
	}
      }
   
      

      if (fIsMC) {
	if (RedTrack->IsHadron()){
	  CheckElectronIsTrigger(ptH, HadContIsTrigger[k]);
	  fMCElecHaHadron->Fill(fillSparse, EventWeight/(recEffH*recEffE));
	}
	if  (RedTrack->IsPhotonic()) {
	  if (RedTrack->TruePartner()) {
	    CheckElectronIsTrigger(ptH, PhotElecWPartnerTrigger[k]);
	    fMCElecHaTruePartner->Fill(fillSparse, EventWeight/(recEffH*recEffE));
	  }
	  else {
	    CheckElectronIsTrigger(ptH,  PhotElecWoPartnerTrigger[k]);
	    fMCElecHaNoPartner->Fill(fillSparse, EventWeight/(recEffH*recEffE));
	  }
	}
	if (fIsAOD) {
	  AliAODMCParticle* MCParticleAOD = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(abs(RedTrack->GetLabel())));  
	  if (!MCParticleAOD) continue;
	  Int_t PDGCode = abs(MCParticleAOD->GetPdgCode());
	  if (PDGCode!=11) {
	    // fBackgroundElecHa->Fill(fillSparse, EventWeight/(recEffH*recEffE));
	  }
	  else if (PDGCode == 11) {
	    Int_t Mother = MCParticleAOD->GetMother();
	    AliAODMCParticle* MCMotherAOD = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(abs(Mother)));
	    Int_t MotherPDG = abs(MCMotherAOD->GetPdgCode());
	    Double_t MotherPt = MCMotherAOD->Pt();
	    Int_t GMother = MCMotherAOD->GetMother();
	    Int_t GMotherPDG = -999;
	    Double_t GMotherPt = -999.;
	    if (abs(GMother)>0) {
	      AliAODMCParticle* MCGMotherAOD = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(abs(GMother)));
	      GMotherPt = MCGMotherAOD->Pt();
	      GMotherPDG =abs(MCGMotherAOD->GetPdgCode());
	    }
	    if (MotherPDG==11 || MotherPDG==15) {
	      Int_t IsHeavyGM =  Int_t (GMotherPDG / TMath::Power(10, Int_t(TMath::Log10(GMotherPDG))));
	      if (IsHeavyGM>3 && IsHeavyGM<6)  fSignalElecHa->Fill(fillSparse, EventWeight/(recEffH*recEffE));
	      else if (!RedTrack->IsPhotonic()){
		Double_t BGWeight = GetBackgroundWeight(MotherPDG, MCParticleAOD->Pt());
		fBackgroundElecHa->Fill(fillSparse, BGWeight*EventWeight/(recEffH*recEffE));
	      }
	    }
	    else {
	      Int_t  IsHeavy =  Int_t (MotherPDG / TMath::Power(10, Int_t(TMath::Log10(MotherPDG))));
	      if (IsHeavy>3 && IsHeavy<6)  fSignalElecHa->Fill(fillSparse, EventWeight/(recEffH*recEffE));
	      else if (!RedTrack->IsPhotonic()) {
		Double_t BGWeight = GetBackgroundWeight(MotherPDG, MCParticleAOD->Pt());
		fBackgroundElecHa->Fill(fillSparse, BGWeight*EventWeight/(recEffH*recEffE));
	      }
	    }
	    if  (RedTrack->IsPhotonic()) {
	      Double_t PtMotherWeight=1.;
	      if (MotherPDG==22) {
		if (GMotherPDG == 221)  PtMotherWeight=GetEtaWeight(GMotherPt); 
		else if (GMotherPDG == 111)  PtMotherWeight=GetPionWeight(GMotherPt); 
	      }
	      else  {
		if (MotherPDG == 221) PtMotherWeight=GetEtaWeight(MotherPt); 
		else if (MotherPDG == 111)  PtMotherWeight=GetPionWeight(MotherPt); 
	      }
	      if (TMath::Abs(PtMotherWeight-1.)<1E-5) {
		//cout << "PDG Info " << PDGCode << " " << MotherPDG <<  " " <<GMotherPDG << endl;
		//cout << MCParticleAOD->Pt() << " " << MotherPt << " " << GMotherPt << endl;
	      }
	      MotherWeight[k]=PtMotherWeight;
	      //   cout << "k " << k << " weight " << MotherWeight[k] << endl;
	    }
	  }  
	}
      }
    }

    // Fill Trigger Hadrons
    for (Int_t AssPtBin=0; AssPtBin<fAssPtElec_Nbins; AssPtBin++) {
      if (HadronTrigger[AssPtBin]) fHadElecTrigger->Fill(ptH, AssPtBin, EventWeight/recEffH);
    }

  } //end of track loop

  // Fill Trigger Histograms
  for (Int_t l=0; l<RedTracksHFE->GetEntriesFast(); l++) {
    AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(l);
    Double_t recEff = GetElectronRecEff(fVevent->GetRunNumber(), RedTrack->Pt(), RedTrack->Phi(), RedTrack->Eta(), pVtx->GetZ());
    if (recEff<0) continue;
    for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
      if (ElectronIsTrigger[l][AssPtBin]) {
	fElecHadTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin, RedTrack->Eta(), 1./recEff);

	for (int j=0; j<RedTrack->LS(); j++) fElecHadTriggerULS->Fill(ElectronIsTriggerPt[l], AssPtBin,  -1.*EventWeight/recEff);
	for (int j=0; j<RedTrack->ULS(); j++) fElecHadTriggerULS->Fill(ElectronIsTriggerPt[l], AssPtBin,  EventWeight/recEff);
      }
      if (ElectronIsTriggerNoP[l][AssPtBin]) {
	for (int j=0; j<RedTrack->LS(); j++) fElecHadTriggerULSNoP->Fill(ElectronIsTriggerPt[l], AssPtBin,   -1*EventWeight/recEff);
	for (int j=0; j<RedTrack->ULS(); j++) fElecHadTriggerULSNoP->Fill(ElectronIsTriggerPt[l], AssPtBin,  EventWeight/recEff);
     	
	for (int j=0; j<RedTrack->LS(); j++) fElecHadTriggerULSNoPCorr->Fill(ElectronIsTriggerPt[l], AssPtBin,   -RedTrack->LSPartnerWeight(j)*EventWeight/recEff);
	for (int j=0; j<RedTrack->ULS(); j++) fElecHadTriggerULSNoPCorr->Fill(ElectronIsTriggerPt[l], AssPtBin,   RedTrack->ULSPartnerWeight(j)*EventWeight/recEff);
      }
      if (ElectronIsTriggerNoPTrue[l][AssPtBin]) {
	for (int j=0; j<RedTrack->LS(); j++) fElecHadTriggerULSNoPCorrTrue->Fill(ElectronIsTriggerPt[l], AssPtBin,  -1.*EventWeight/recEff);
	for (int j=0; j<RedTrack->ULS(); j++) fElecHadTriggerULSNoPCorrTrue->Fill(ElectronIsTriggerPt[l], AssPtBin,  EventWeight/recEff);
      }
      if (PhotElecWPartnerTrigger[l][AssPtBin]) {
	fMCElecHaTruePartnerTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin, EventWeight/recEff);
	fMCElecHaTruePartnerTriggerWW->Fill(ElectronIsTriggerPt[l], AssPtBin, MotherWeight[l]*EventWeight/recEff);
      }
      if (PhotElecWoPartnerTrigger[l][AssPtBin]) {
	fMCElecHaNoPartnerTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin, EventWeight/recEff);
	fMCElecHaNoPartnerTriggerWW->Fill(ElectronIsTriggerPt[l], AssPtBin, MotherWeight[l]*EventWeight/recEff);
      }
      if (HadContIsTrigger[l][AssPtBin]) fHadContTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin, EventWeight/recEff);
    }
  }


  for (Int_t l=0; l<RedTracksHFE->GetEntriesFast(); l++) {
    AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(l);
    Double_t recEff = GetElectronRecEff(fVevent->GetRunNumber(), RedTrack->Pt(), RedTrack->Phi(), RedTrack->Eta(), pVtx->GetZ());
    if (recEff<0) continue;
    //  for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
    Bool_t NoPNoT = kTRUE;
    Bool_t TPNoT = kTRUE;
    if (RedTrack->IsPhotonic()) { // missleading notation - only trigger particles are stored
      for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
	Double_t fillSparse[4]={RedTrack->Pt(), RedTrack->TruePartnerRecPt(), RedTrack->TruePartnerMCPt(), 1.0*AssPtBin};
	if (PhotElecWPartnerTrigger[l][AssPtBin])  fTPartnerNoTPt2->Fill(fillSparse, EventWeight*MotherWeight[l]);
	if (PhotElecWoPartnerTrigger[l][AssPtBin]) fNoPartnerNoTPt2->Fill(fillSparse, EventWeight*MotherWeight[l]);

	if (PhotElecWPartnerTrigger[l][AssPtBin] && TPNoT) TPNoT=kFALSE;
	if (PhotElecWoPartnerTrigger[l][AssPtBin] && NoPNoT) NoPNoT = kFALSE;
      }
      if (RedTrack->TruePartner() && TPNoT) {
	fTPartnerNoT->Fill(RedTrack->Pt(), EventWeight/recEff); //cout << "TruePartnerNoTrigger" << endl;
      }
      if (!RedTrack->TruePartner() && NoPNoT) {
	fNoPartnerNoT->Fill(RedTrack->Pt(), EventWeight/recEff); // cout << "NoPartnerNoTrigger" << endl;
      }
    }
  }
  
  for (int l=0; l<ntracks; l++) {
    for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
      if (NonElectronIsTrigger[l][AssPtBin]) fNonElecHadTrigger->Fill(NonElectronIsTriggerPt[l], AssPtBin, EventWeight/NonElectronIsTriggerWeight[l]);
    }
  }
  
  // Clear Trigger Arrays
  for (Int_t j=0; j < RedTracksHFE->GetEntriesFast(); j++) {
    delete [] ElectronIsTrigger[j];
    delete [] ElectronIsTriggerNoP[j];
    delete [] ElectronIsTriggerNoPTrue[j];
    delete [] HadContIsTrigger[j];
    delete [] PhotElecWPartnerTrigger[j];
    delete [] PhotElecWoPartnerTrigger[j];
  }
  delete [] ElectronIsTrigger;
  delete [] ElectronIsTriggerNoP;
  delete [] ElectronIsTriggerNoPTrue;
  delete [] HadContIsTrigger;
  delete [] PhotElecWPartnerTrigger;
  delete [] PhotElecWoPartnerTrigger;
  delete [] ElectronIsTriggerPt;
  delete [] MotherWeight;
  for (Int_t j=0; j<ntracks; j++) delete[] NonElectronIsTrigger[j];
  delete [] NonElectronIsTrigger;
  delete [] NonElectronIsTriggerPt;
  delete [] NonElectronIsTriggerWeight;
  delete [] HadronTrigger;
  
}


void AliAnalysisTaskHaHFECorrel::CorrelateHadronMixedEvent(Float_t mult, const AliVVertex* pVtx, Float_t maxPt, Int_t nMother, Int_t listMother[], Bool_t EvContTP, Bool_t EvContNTP, Double_t EventWeight) {
  AliEventPool* fPool;
  fPool = fPoolMgr->GetEventPool(mult, pVtx->GetZ(), maxPt); // Get the buffer associated with the current centrality and z-vtx
  //  fPool->SetDebug(kTRUE);
  if (!fPool)
  {
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f, maxPt = %f", mult, pVtx->GetZ(), maxPt));
    return;
  }
  
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();

  // Track loop for hadron correlations
  for(Int_t iTracks = 0; iTracks < ntracks; iTracks++) {     
  
    AliVParticle* Vtrack = 0x0;
    if(!fUseTender) Vtrack  = fVevent->GetTrack(iTracks);
    if(fUseTender)  Vtrack  = dynamic_cast<AliVTrack*>(fTracks_tender->At(iTracks));
    if (!Vtrack) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }
    AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
    AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);
    AliESDtrack *ESDtrack = dynamic_cast<AliESDtrack*>(Vtrack);
    if(!track) continue;
    if (fIsAOD && !AODtrack) continue;
    if (!fIsAOD && !ESDtrack) continue;
        
    // track cuts
    Bool_t passTrackCut=kFALSE;
    passTrackCut = ChargedHadronTrackCuts(pVtx, track, nMother, listMother, EventWeight);
    if (!passTrackCut) continue;

    // pid cuts
    Bool_t passPIDCut=kFALSE;
    passPIDCut = ChargedHadronPIDCuts(track, EventWeight);
    if (!passPIDCut) continue;
   
    Double_t ptH=-9.,etaH =-9.,phiH=-9., recEffH=-9.;
    Int_t idH=-9;
    ptH = track->Pt();
    phiH = track->Phi();
    etaH = track->Eta();
    idH = track->GetID();
    recEffH = GetHadronRecEff(fVevent->GetRunNumber(), ptH, phiH, etaH, pVtx->GetZ());
    if (recEffH<0) continue;
 
    if (fPool->GetCurrentNEvents() >= 3){// start mixing when 3 events are in the buffer
      Int_t nMix = fPool->GetCurrentNEvents();   
      for (Int_t jMix=0; jMix<nMix; jMix++){  // mix with each event in the buffer
	TObjArray* mixTracks = fPool->GetEvent(jMix);

       	for (Int_t i=0;i<mixTracks->GetEntriesFast(); i++) {	
	  AliBasicParticleHaHFE* mixtrk = (AliBasicParticleHaHFE*) mixTracks->At(i);
	  if (!mixtrk) {
	    printf("ERROR: Could not receive mix pool track %d\n",i);
	    continue;
	  }
	  
	  Double_t ptMix=-9., phiMix=-9., etaMix=-9, recEff=-9;
	  Int_t ls=-9, uls=-9;
	  ptMix=mixtrk->Pt();
	  phiMix=mixtrk->Phi();
	  etaMix=mixtrk->Eta();
	  ls=mixtrk->LS();
	  uls=mixtrk->ULS();
	  recEff=GetElectronRecEff(fVevent->GetRunNumber(), ptMix,phiMix,etaMix,pVtx->GetZ()); // assumption same vtx as original event (same Vtx binning)
	  if (recEff<0) continue;
	  
	  Double_t dphi=-99, deta=-99;
	  dphi=GetDeltaPhi(phiMix, phiH);
	  deta=GetDeltaEta(etaMix, etaH);
		
	  Double_t fillSparse[5]={-999,-999,-999,-999,-999};
	  fillSparse[0]=ptH;
	  fillSparse[1]=ptMix;
	  fillSparse[2]=dphi;
	  fillSparse[3]=deta;
	  fillSparse[4]=pVtx->GetZ();
	  
	  fElecHaMixedEvent->Fill(fillSparse, 1./(recEff*recEffH));
	  for (Int_t k=0; k<ls; k++)  fULSElecHaMixedEvent->Fill(fillSparse, -1./(recEff*recEffH)); // Remark: wo EventWeights
	  for (Int_t k=0; k<uls; k++) fULSElecHaMixedEvent->Fill(fillSparse, 1./(recEff*recEffH)); // Remark: wo EventWeights
	  
	  if (fIsMC) {
	    if (mixtrk->IsPhotonic()) {
	      if (mixtrk->TruePartner()) {
		fTagHaMixedEvent->Fill(fillSparse, 1./(recEff*recEffH));
	      }
	      else if (!mixtrk->TruePartner()) {
		fNonTagHaMixedEvent->Fill(fillSparse, 1./(recEff*recEffH));
	      }
	    }
	  }
	}
      }
    }
  }
  return;
}

/////////////--------------

void AliAnalysisTaskHaHFECorrel::CorrelateWithHadrons(AliVTrack* Vtrack, const AliVVertex* pVtx, Int_t nMother, Int_t listMother[], Bool_t FillHadron, Bool_t FillLP, Bool_t **NonElecIsTrigger, Double_t *NonElecIsTriggerPt, Double_t * NonElecIsTriggerWeight, Int_t NumElectronsInEvent, Double_t EventWeight) {
  
  // Trigger Hadron
  Double_t ptH=-9.,etaH =-9.,phiH=-9., recEffH=-9;
  Int_t idH=-9;
  ptH = Vtrack->Pt();
  phiH = Vtrack->Phi();
  etaH = Vtrack->Eta();
  idH = Vtrack->GetID();
  recEffH = GetHadronRecEff(fVevent->GetRunNumber(), ptH,phiH,etaH, pVtx->GetZ());
  if (recEffH<0) return;

  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();

  Bool_t Trigger[fAssPtElec_Nbins];
  for (Int_t AssPt=0; AssPt<fAssPtElec_Nbins; AssPt++) Trigger[AssPt]=kFALSE;           // flag to store information, if used as trigger particle
  Int_t CountFakeElectrons=0;

  // Track loop for wrong electron candidates (TPCs<-4)
  for(Int_t iTracks = 0; iTracks < ntracks; iTracks++) {     
    
    AliVParticle* Vtrack = 0x0;
    if(!fUseTender) Vtrack  = fVevent->GetTrack(iTracks);
    if(fUseTender)  Vtrack  = dynamic_cast<AliVTrack*>(fTracks_tender->At(iTracks));
    if (!Vtrack) {
      printf("ERROR: Could not receive tagged  track %d\n", iTracks);
      continue;
    }
    AliVTrack *VtrackAssoHad = dynamic_cast<AliVTrack*>(Vtrack);
    if(!VtrackAssoHad) continue;
        
    if(VtrackAssoHad->GetID()==idH) continue;

    // track cuts
    Bool_t passTrackCut=kFALSE;
    passTrackCut = InclElecTrackCuts(pVtx, VtrackAssoHad, nMother, listMother, EventWeight);
    if (!passTrackCut) continue;

    // pid cuts
    Bool_t passPIDCut=kFALSE;
    passPIDCut = AssoHadronPIDCuts(VtrackAssoHad, EventWeight);
    if (!passPIDCut) continue;
   
    // associated particle
    Double_t ptAsso=-9.,etaAsso=-9.,phiAsso=-9.,recEffAsso=-9;
    ptAsso = VtrackAssoHad->Pt();
    phiAsso = VtrackAssoHad->Phi();
    etaAsso = VtrackAssoHad->Eta();
    recEffAsso = GetElectronRecEff(fVevent->GetRunNumber(), ptAsso, phiAsso, etaAsso, pVtx->GetZ());
    if (recEffAsso<0) continue;

    CountFakeElectrons++;
    


    
    CheckHadronIsTrigger(ptAsso, Trigger); // for hadrons check that ptWrongElectron=ptAsso is trigger
    CheckElectronIsTrigger(ptH, NonElecIsTrigger[iTracks]); // for 
    NonElecIsTriggerPt[iTracks]=ptAsso;
    NonElecIsTriggerWeight[iTracks]=recEffAsso;
    if(FillLP) {
      for (Int_t ptAssHad=0; ptAssHad<fAssPtHad_Nbins; ptAssHad++) {
	if (NonElecIsTrigger[iTracks][ptAssHad]) fNonElecLPTrigger->Fill(ptAsso, ptAssHad, EventWeight/recEffAsso);
      }
    }
  
    Double_t dphi = -99, deta=-99;
    dphi = GetDeltaPhi(phiAsso, phiH);
    deta = GetDeltaEta(etaAsso, etaH);

    // Fill Sparse
    Double_t fillSparse[5]={-999,-999,-999,-999,-999};
    fillSparse[0]=ptH;
    fillSparse[1]=ptAsso;
    fillSparse[2]=dphi;
    fillSparse[3]=deta;
    fillSparse[4]=pVtx->GetZ();
    if (FillHadron) fElecHaHa->Fill(fillSparse, EventWeight/(recEffH*recEffAsso));
    if (FillLP)     fElecLPHa->Fill(fillSparse, EventWeight/(recEffAsso));

    if (CountFakeElectrons>NumElectronsInEvent) iTracks=ntracks;
  }

  if (FillHadron) {
    for (Int_t ptAssElec=0; ptAssElec<fAssPtElec_Nbins; ptAssElec++) {
      if (Trigger[ptAssElec]) fHadNonElecTrigger->Fill(ptH, ptAssElec, EventWeight/(recEffH));
    }
  }
  if (FillLP){
    for (Int_t ptAssElec=0; ptAssElec<fAssPtElec_Nbins; ptAssElec++) {
      if (Trigger[ptAssElec]) fLPNonElecTrigger->Fill(ptH, ptAssElec, EventWeight/recEffH);
    }
  }
}



//_______________________

void AliAnalysisTaskHaHFECorrel::CorrelateLP(AliVTrack* LPtrack,  const AliVVertex* pVtx, Int_t nMother, Int_t listMother[], TObjArray *RedTracksHFE, Double_t EventWeight) 
{ 

  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();

  Bool_t *LPTrigger= new Bool_t[fAssPtElec_Nbins];
  for (Int_t l=0; l<fAssPtElec_Nbins; l++) LPTrigger[l]=kFALSE;
  Bool_t **ElectronIsTrigger = new Bool_t*[RedTracksHFE->GetEntriesFast()];
  Bool_t **ElectronIsTriggerNoP = new Bool_t*[RedTracksHFE->GetEntriesFast()];
  Bool_t **HadContIsTrigger=new Bool_t*[RedTracksHFE->GetEntriesFast()];
  Bool_t **PhotElecWPartnerTrigger=new Bool_t*[RedTracksHFE->GetEntriesFast()];
  Bool_t **PhotElecWoPartnerTrigger=new Bool_t*[RedTracksHFE->GetEntriesFast()];
  for (Int_t j=0; j<RedTracksHFE->GetEntriesFast(); j++)
    {
      ElectronIsTrigger[j]=new Bool_t[fAssPtHad_Nbins];
      ElectronIsTriggerNoP[j]=new Bool_t[fAssPtHad_Nbins];
      HadContIsTrigger[j]=new Bool_t[fAssPtHad_Nbins];
      PhotElecWPartnerTrigger[j] = new Bool_t[fAssPtHad_Nbins];
      PhotElecWoPartnerTrigger[j] = new Bool_t[fAssPtHad_Nbins];
      for (Int_t m=0; m<fAssPtHad_Nbins; m++) {
	ElectronIsTrigger[j][m]=kFALSE;
	ElectronIsTriggerNoP[j][m]=kFALSE;
	HadContIsTrigger[j][m]=kFALSE;
	PhotElecWPartnerTrigger[j][m]=kFALSE;
	PhotElecWoPartnerTrigger[j][m]=kFALSE;
      }
    }
  Double_t *ElectronIsTriggerPt = new Double_t[RedTracksHFE->GetEntriesFast()];
  Bool_t **NonElectronIsTrigger = new Bool_t*[ntracks];
  for (Int_t j=0; j<ntracks; j++) { 
    NonElectronIsTrigger[j]=new Bool_t[fAssPtHad_Nbins]; 
    for (Int_t m=0; m<fAssPtHad_Nbins; m++) {
      NonElectronIsTrigger[j][m]=kFALSE;
    }
  }
  Double_t *NonElectronIsTriggerPt = new Double_t[ntracks];
  Double_t *NonElectronIsTriggerWeight = new Double_t[ntracks];

  // leading Particle neq 
  Double_t ptH=-9.,etaH =-9.,phiH=-9.,recEffH=-9.;
  Int_t idH=-9;
  ptH = LPtrack->Pt();
  phiH = LPtrack->Phi();
  etaH = LPtrack->Eta();
  idH = LPtrack->GetID();
  recEffH = GetHadronRecEff(fVevent->GetRunNumber(), ptH, phiH, etaH, pVtx->GetZ());
  if (recEffH<0) return;

  CorrelateWithHadrons(LPtrack, pVtx, nMother, listMother, kFALSE, kTRUE, NonElectronIsTrigger, NonElectronIsTriggerPt, NonElectronIsTriggerWeight, RedTracksHFE->GetEntriesFast(), EventWeight); // correlate LPHadron (kFALSE, kTRUE);
		       
  // Only loop over electrons in event
  for (Int_t k=0; k<RedTracksHFE->GetEntriesFast(); k++) {
    AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(k);
  

    if (RedTrack->ID()==idH) continue; // no self-correlation
    Double_t pt=-9.,eta =-9.,phi=-9.,recEffE=-9.;
    Int_t ls=0, uls=0;    
    pt = RedTrack->Pt();
    phi = RedTrack->Phi();
    eta = RedTrack->Eta();
    ls = RedTrack->LS();
    uls = RedTrack->ULS();
    recEffE = GetElectronRecEff(fVevent->GetRunNumber(), pt,phi,eta,pVtx->GetZ());
    if (recEffE<0) continue;

    Double_t dphi = -99, deta = -99;
    dphi = GetDeltaPhi(phi, phiH);
    deta = GetDeltaEta(eta, etaH);

    Double_t fillSparse[5]={-999,-999,-999,-999,-999};
    fillSparse[0]=ptH;
    fillSparse[1]=pt;
    fillSparse[2]=dphi;
    fillSparse[3]=deta;
    fillSparse[4]=pVtx->GetZ();


    CheckElectronIsTrigger(ptH, ElectronIsTrigger[k]);
    ElectronIsTriggerPt[k]=pt;
    CheckHadronIsTrigger(pt, LPTrigger);
    fInclElecLP->Fill(fillSparse,EventWeight/(recEffE));
    for (Int_t j=0; j<ls; j++)  fULSElecLP->Fill(fillSparse, -EventWeight/(recEffE));
    for (Int_t j=0; j<uls; j++) fULSElecLP->Fill(fillSparse, EventWeight/(recEffE));


    Bool_t HadIsULSPartner=kFALSE;
    Bool_t HadIsLSPartner=kFALSE;
    for (int j=0; j<uls; j++) {
      if (idH==RedTrack->ULSPartner(j)) HadIsULSPartner =kTRUE;
    }
    for (int j=0; j<ls; j++) {
      if (idH==RedTrack->LSPartner(j)) HadIsLSPartner=kTRUE;
    }
    if (!HadIsULSPartner && !HadIsLSPartner) {
      CheckElectronIsTrigger(ptH, ElectronIsTriggerNoP[k]);
      if (uls>0 || ls>0) {
	for (int j=0; j<ls; j++) fElecLPULSNoPartner->Fill(fillSparse, -1.*EventWeight/(recEffE));
	for (int j=0; j<uls; j++) fElecLPULSNoPartner->Fill(fillSparse, EventWeight/(recEffE));
    	for (int j=0; j<ls; j++) {
	  fElecLPULSNoPartnerCorr->Fill(fillSparse, -RedTrack->LSPartnerWeight(j)*EventWeight/(recEffE));
	}
	for (int j=0; j<uls; j++) {
	  fElecLPULSNoPartnerCorr->Fill(fillSparse, RedTrack->ULSPartnerWeight(j)*EventWeight/(recEffE));
	}
	if (RedTrack->IsPhotonic()) {
	  //  CheckElectronIsTrigger(ptH, ElectronIsTriggerNoPTrue[k]);
	  //  fElecLPULSNoPartnerCorrTrue->Fill(fillSparse, uls*EventWeight/(recEffE));
	  // fElecLPULSNoPartnerCorrTrue->Fill(fillSparse, -ls*EventWeight/(recEffE));
	}
      }
    }


    

    if (fIsMC) {
      if (RedTrack->IsHadron()){
	CheckElectronIsTrigger(ptH, HadContIsTrigger[k]);
	fMCElecLPHadron->Fill(fillSparse, EventWeight/(recEffE));
      }
      if (RedTrack->IsPhotonic()) {
	if (RedTrack->TruePartner()) {
	  CheckElectronIsTrigger(ptH, PhotElecWPartnerTrigger[k]);
	  fMCElecLPTruePartner->Fill(fillSparse, EventWeight/(recEffE));
	}
	else {
	  CheckElectronIsTrigger(ptH,  PhotElecWoPartnerTrigger[k]);
	  fMCElecLPNoPartner->Fill(fillSparse, EventWeight/(recEffE));
	}
      }
    }
  }
   
  //save LP as Trigger (existens of LP and HFE are satisfied in function call)
  for (Int_t PtAssElec=0; PtAssElec<fAssPtElec_Nbins; PtAssElec++) {
    if (LPTrigger[PtAssElec]) fLPElecTrigger->Fill(ptH, PtAssElec, EventWeight);
  }


  // Fill Trigger Histograms
  for (Int_t l=0; l<RedTracksHFE->GetEntriesFast(); l++) {
    AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(l);
    Double_t recEffE = GetElectronRecEff(fVevent->GetRunNumber(), RedTrack->Pt(), RedTrack->Phi(), RedTrack->Eta(), pVtx->GetZ());
    if (recEffE<0) continue;
    for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
      if (ElectronIsTrigger[l][AssPtBin]) {
	fElecLPTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin, RedTrack->Eta(), 1.0/recEffE);
	for (int j=0; j<RedTrack->LS(); j++) fElecLPTriggerULS->Fill(ElectronIsTriggerPt[l], AssPtBin,   -EventWeight*1./recEffE);
	for (int j=0; j<RedTrack->ULS(); j++) fElecLPTriggerULS->Fill(ElectronIsTriggerPt[l], AssPtBin,  EventWeight/recEffE);
      }
      if (ElectronIsTriggerNoP[l][AssPtBin]) {
	for (int j=0; j<RedTrack->LS(); j++)  fElecLPTriggerULSNoP->Fill(ElectronIsTriggerPt[l], AssPtBin, -EventWeight*1./recEffE);
	for (int j=0; j<RedTrack->ULS(); j++) fElecLPTriggerULSNoP->Fill(ElectronIsTriggerPt[l], AssPtBin,  EventWeight/recEffE);
	for (int j=0; j<RedTrack->LS(); j++) fElecLPTriggerULSNoPCorr->Fill(ElectronIsTriggerPt[l], AssPtBin,   -RedTrack->LSPartnerWeight(j)*EventWeight/recEffE);
	for (int j=0; j<RedTrack->ULS(); j++) fElecLPTriggerULSNoPCorr->Fill(ElectronIsTriggerPt[l], AssPtBin,   RedTrack->ULSPartnerWeight(j)*EventWeight/recEffE);	
      }
      if (PhotElecWPartnerTrigger[l][AssPtBin]) fMCElecLPTruePartnerTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin, EventWeight/recEffE);
      if (PhotElecWoPartnerTrigger[l][AssPtBin]) fMCElecLPNoPartnerTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin, EventWeight/recEffE);
      if (HadContIsTrigger[l][AssPtBin]) fHadContLPTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin, EventWeight/recEffE);
    }
  }

  for (int l=0; l<ntracks; l++) {
    for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
      if (NonElectronIsTrigger[l][AssPtBin]) fNonElecLPTrigger->Fill(NonElectronIsTriggerPt[l], AssPtBin, EventWeight/NonElectronIsTriggerWeight[l]);
    }
  }


  // Clear Trigger Arrays
  for (Int_t j=0; j < RedTracksHFE->GetEntriesFast(); j++) {
    delete [] ElectronIsTrigger[j];
    delete [] ElectronIsTriggerNoP[j];
    delete [] HadContIsTrigger[j];
    delete [] PhotElecWPartnerTrigger[j];
    delete [] PhotElecWoPartnerTrigger[j];
  }
  delete [] ElectronIsTrigger;
  delete [] ElectronIsTriggerNoP;
  delete [] HadContIsTrigger;
  delete [] PhotElecWPartnerTrigger;
  delete [] PhotElecWoPartnerTrigger;
  delete [] LPTrigger;
  delete [] ElectronIsTriggerPt;
  for (Int_t j=0; j<ntracks; j++) delete[] NonElectronIsTrigger[j];
  delete [] NonElectronIsTrigger;
  delete [] NonElectronIsTriggerPt;
  delete [] NonElectronIsTriggerWeight;

}



//_______________________

void AliAnalysisTaskHaHFECorrel::CorrelateLPMixedEvent(AliVTrack* LPtrack, Float_t mult, Float_t zVtx, Float_t maxPt, Bool_t EvContTP, Bool_t EvContNTP, Double_t EventWeight) 
{ 
  // leading Particle 
  Double_t ptH=-9.,etaH =-9.,phiH=-9.,recEffH=-9.;
  ptH = LPtrack->Pt();
  phiH = LPtrack->Phi();
  etaH = LPtrack->Eta();
  recEffH = GetHadronRecEff(fVevent->GetRunNumber(), ptH, phiH, etaH, zVtx);
  if (recEffH<0) return;

  AliEventPool* fPool;
  fPool = fPoolMgr->GetEventPool(mult, zVtx, maxPt); // Get the buffer associated with the current centrality and z-vtx
  //  fPool->SetDebug(kTRUE);
  if (!fPool)
  {
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f, maxPt = %f", mult, zVtx, maxPt));
    return;
  }
  
  if (fPool->GetCurrentNEvents() >= 3) // start mixing when 3 events are in the buffer
  {
    Int_t nMix = fPool->GetCurrentNEvents();   
    for (Int_t jMix=0; jMix<nMix; jMix++)  // mix with each event in the buffer
    {
      TObjArray* mixTracks = fPool->GetEvent(jMix);
      for (Int_t i=0;i<mixTracks->GetEntriesFast(); i++)
      {	
        AliBasicParticleHaHFE* mixtrk = (AliBasicParticleHaHFE*) mixTracks->At(i);
	if (!mixtrk) {
          printf("ERROR: Could not receive mix pool track %d\n",i);
          continue;
        }
	
	Double_t ptMix=-9., phiMix=-9., etaMix=-9, recEffE=-9.;
	Int_t ls=-9, uls=-9;
	ptMix=mixtrk->Pt();
	phiMix=mixtrk->Phi();
	etaMix=mixtrk->Eta();
	ls=mixtrk->LS();
	uls=mixtrk->ULS();
	recEffE = GetElectronRecEff(fVevent->GetRunNumber(), ptMix, phiMix, etaMix, zVtx);
	if (recEffE<0) continue;

	Double_t dphi=-99, deta=-99;
	dphi=GetDeltaPhi(phiMix, phiH);
	deta=GetDeltaEta(etaMix, etaH);

	Double_t fillSparse[5]={-999,-999,-999,-999,-999};
	fillSparse[0]=ptH;
	fillSparse[1]=ptMix;
	fillSparse[2]=dphi;
	fillSparse[3]=deta;
	fillSparse[4]=zVtx;

	fElecLPMixedEvent->Fill(fillSparse, 1./(recEffE));
	for (Int_t k=0; k<ls; k++)  fULSElecLPMixedEvent->Fill(fillSparse, -1./(recEffE)); // Remark NoEvWeights
	for (Int_t k=0; k<uls; k++) fULSElecLPMixedEvent->Fill(fillSparse, 1./(recEffE)); // Remark NoEvWeights

	if (fIsMC) {
	  if (mixtrk->IsPhotonic()) {
	    if (mixtrk->TruePartner()) {
	      fTagLPMixedEvent->Fill(fillSparse, 1./(recEffE)); // NoEvWeights
	    }
	    else if (!mixtrk->TruePartner()) {
	      fNonTagLPMixedEvent->Fill(fillSparse, 1./(recEffE)); // NoEvWeights
	    }
	  }
       	}
      }
    }
  }
}


//_________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::InclElecTrackCuts(const AliVVertex *pVtx,AliVTrack *Vtrack,Int_t nMother, Int_t listMother[], Double_t EventWeight, Bool_t fillHists )
{
  AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);   
  AliESDtrack *ESDtrack = dynamic_cast<AliESDtrack*>(Vtrack);   
  if (fIsAOD && !AODtrack) return kFALSE;
  if (!fIsAOD && !ESDtrack) return kFALSE;

  // track cuts for inclusive electrons
 

  /*

  if (fIsAOD) { 
    cout << "TPCFraction" << endl;
    cout << AODtrack->GetTPCFoundFraction() << endl;
    cout << AODtrack->GetTPCClusterInfo(2,0) << endl;  // not always same variable
    fElectronTrackTPCFrac->Fill(AODtrack->Pt(), AODtrack->GetTPCFoundFraction()); 
  }
  //  if(ietrack->PropagateToDCA(pVietx, fVevent->GetMagneticField(), 20., d0z0, cov))
  //fElectronTrackDCA->Fill(TMath::Abs(d0z0[0]), TMath::Abs(d0z0[1]));
  */

  // First Acceptance Cuts and Filter bit
  if (fillHists) fElectronTrackCuts->Fill(Vtrack->Pt(), 0., EventWeight);

  if (fIsAOD) { 
    if(!AODtrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE;
  }
  else if(!ESDkTrkGlobalNoDCA(Vtrack)) return kFALSE;

  

  

  if (fillHists)  fElectronTrackCuts->Fill(Vtrack->Pt(), 1., EventWeight);

  if(Vtrack->Pt()<0.25) return kFALSE;
   if (fillHists)  fElectronTrackCuts->Fill(Vtrack->Pt(), 2., EventWeight);

  if(Vtrack->Eta()<fMinElectronEta) return kFALSE;
  if(Vtrack->Eta()>fMaxElectronEta) return kFALSE;
  if (fillHists) fElectronTrackCuts->Fill(Vtrack->Pt(), 3., EventWeight);


  // Store QA informations



  Double_t d0z0[2]={-999,-999}, cov[3];
  Float_t DCAr=99, DCAz=99;
  if(fIsAOD) {
    AODtrack->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov);
    DCAr=d0z0[0];
    DCAz=d0z0[1];
  }
  else ESDtrack->GetImpactParameters(DCAr, DCAz);

  Int_t ITSNcls;
  if (fIsAOD) ITSNcls = AODtrack->GetITSNcls();
  else ITSNcls = ESDtrack->GetNcls(0);

  Bool_t ITSPointOnLayer0=kFALSE;
  Bool_t ITSPointOnLayer1=kFALSE;
  ITSPointOnLayer0=Vtrack->HasPointOnITSLayer(0);
  ITSPointOnLayer1=Vtrack->HasPointOnITSLayer(1);
  
  //* ITS shared cluster cut - not useful for single electron studies
  Int_t ITSNclsShared = 0;
  for(int i=0; i<6; i++){
    if( Vtrack->HasSharedPointOnITSLayer(i) ) ITSNclsShared++;
  }
  Float_t ITSFracSharedNCls = ITSNclsShared/(ITSNcls*1.);
  Float_t ITSRedChi2 = Vtrack->GetITSchi2()/(ITSNcls*1);


  if (fillHists){
    if (fIsAOD) fElectronTrackTPCChi2->Fill(Vtrack->Pt(), AODtrack->GetTPCchi2perCluster(), EventWeight);
    fElectronTrackTPCCrossedRows->Fill(Vtrack->Pt(), Vtrack->GetTPCCrossedRows(), EventWeight);
    fElectronTrackTPCNcls->Fill(Vtrack->Pt(), Vtrack->GetTPCNcls(), EventWeight);
    fElectronTrackTPCNclsdEdx->Fill(Vtrack->Pt(), Vtrack->GetTPCsignalN(), EventWeight);
    fElectronTrackITSNcls->Fill(Vtrack->Pt(),ITSNcls, EventWeight);
    fElectronTrackITSChi2->Fill(Vtrack->Pt(), ITSRedChi2, EventWeight);
    fElectronTrackITSLayer->Fill(Vtrack->Pt(),ITSPointOnLayer0, ITSPointOnLayer1, EventWeight);
    fElectronTrackRefit->Fill(Vtrack->Pt(), Vtrack->IsOn(AliAODTrack::kTPCrefit), Vtrack->IsOn(AliAODTrack::kITSrefit), EventWeight);
    fElectronTrackDCA->Fill(Vtrack->Pt(), DCAr, DCAz, EventWeight);
  }

  //MC study of shared cluster cut

  
  Int_t ParticleID=-99;
  if (fIsMC && fIsAOD && (ITSPointOnLayer0 && ITSPointOnLayer1)) {
    Int_t MClabel=AODtrack->GetLabel();
    AliAODMCParticle* MCParticleAOD = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(abs(MClabel)));  
    if (!MCParticleAOD) return kFALSE;
    Int_t PDGCode = abs(MCParticleAOD->PdgCode());  
    if (PDGCode==11) {
      Int_t Mother = abs(MCParticleAOD->GetMother());
      AliAODMCParticle* MCMotherAOD =  dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(Mother));  
      if (!MCMotherAOD) return kFALSE;;
      PDGCode = abs(MCMotherAOD->PdgCode());  
      if (PDGCode==22 || PDGCode==111 || PDGCode==221 || PDGCode==223 || PDGCode==333) ParticleID = PDGCode;
      if  (Int_t (PDGCode / TMath::Power(10, Int_t(TMath::Log10(PDGCode)))) ==4) ParticleID = 4;
      if  (Int_t (PDGCode / TMath::Power(10, Int_t(TMath::Log10(PDGCode)))) ==5) ParticleID = 5;
    }
    Double_t fillSparse[5]={Vtrack->Pt(), 1.*ITSNcls, 1.*ITSNclsShared, 1.*ITSRedChi2, 1.*ParticleID};
    //  fElectronTrackITSCuts->Fill(fillSparse);
  }








  if(Vtrack->GetTPCNcls() < fTPCnCut) return kFALSE;
  if (fillHists) fElectronTrackCuts->Fill(Vtrack->Pt(), 4, EventWeight);
  if(Vtrack->GetTPCsignalN() < fTPCndEdxCut) return kFALSE ;
  if (fillHists) fElectronTrackCuts->Fill(Vtrack->Pt(), 5, EventWeight);
						
  if (!Vtrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
  if (fillHists) fElectronTrackCuts->Fill(Vtrack->Pt(), 6, EventWeight);
  //if(ietrack->GetTPCFoundFraction() < 0.6) return kFALSE;   
  if (fillHists) fElectronTrackCuts->Fill(Vtrack->Pt(), 7, EventWeight);

  if (fIsAOD) {
    Bool_t kinkmotherpass = kTRUE;
    for(Int_t kinkmother = 0; kinkmother < nMother; kinkmother++) {
      if(Vtrack->GetID() == listMother[kinkmother]) {
	kinkmotherpass = kFALSE;
	continue;
      }
    }
    if(!kinkmotherpass) return kFALSE;
  }
  else {
    AliESDtrack *Htrack = dynamic_cast<AliESDtrack*>(Vtrack);   
    if(!Htrack) return kFALSE;
    if(Htrack->GetKinkIndex(0) != 0) return kFALSE;
  }
  if (fillHists) fElectronTrackCuts->Fill(Vtrack->Pt(), 8, EventWeight);

  if(ITSNcls < fITSnCut) return kFALSE;
  if (fillHists) fElectronTrackCuts->Fill(Vtrack->Pt(), 9, EventWeight);

  if (!Vtrack->IsOn(AliAODTrack::kITSrefit)) return kFALSE ;   
  if (fillHists) fElectronTrackCuts->Fill(Vtrack->Pt(), 10, EventWeight);

  if(!(ITSPointOnLayer0 || ITSPointOnLayer1)) return kFALSE; // now kAny
  if (fillHists)fElectronTrackCuts->Fill(Vtrack->Pt(), 11, EventWeight);

  if(!(ITSPointOnLayer0) && fElectronkFirst) return kFALSE;

  if(!(ITSPointOnLayer0 && ITSPointOnLayer1) && !(fElectronkAny || fElectronkFirst)) return kFALSE; // now kBoth
  if (fillHists) fElectronTrackCuts->Fill(Vtrack->Pt(), 12, EventWeight);
  

  if(TMath::Abs(DCAr) > fEleDCAr || TMath::Abs(DCAz) > fEleDCAz) return kFALSE;
  if (fillHists) fElectronTrackCuts->Fill(Vtrack->Pt(), 13, EventWeight); 


  if (fITSSharedClusterCut>0.9) {
    //    cout << "ITSShardeClusterCUt" << fITSSharedClusterCut << "\t" << ITSNclsShared << endl;
    if (ITSNclsShared>fITSSharedClusterCut) return kFALSE;
  }
  else {
    if (ITSFracSharedNCls>fITSSharedClusterCut) return kFALSE;
  }
  // cout << "Cut " << fITSSharedClusterCut << "\t" << ITSNclsShared << endl;
  if (fillHists) fElectronTrackCuts->Fill(Vtrack->Pt(), 14, EventWeight); 
 


 if (fUseTRD) {
    
    if (fIsAOD && AODtrack->GetTRDntrackletsPID()<4) return kFALSE;
    else if (!fIsAOD && ESDtrack->GetTRDntrackletsPID()<4) return kFALSE;

    if (fIsAOD && AODtrack->GetTRDncls()/ AODtrack->GetTRDntrackletsPID()<17) return kFALSE;
    else if (!fIsAOD && ESDtrack->GetTRDncls()/ESDtrack->GetTRDntrackletsPID()<17) return kFALSE;
  }


 


  return kTRUE;

 }

 //________________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::InclElecPIDCuts(AliVTrack* Vtrack,Double_t EventWeight,  Bool_t fillHists ) {

  AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack); 
  if (fIsAOD && !AODtrack) return kFALSE;
  AliESDtrack *ESDtrack = dynamic_cast<AliESDtrack*>(Vtrack); 
  if (!fIsAOD && !ESDtrack) return kFALSE;

  // PID cuts
  Double_t fITSnSigma=-9.,fTOFnSigma=-9.,fTPCnSigma=-9.;
  Bool_t PassITSCut=kTRUE, PassTPCCut=kTRUE, PassTOFCut=kTRUE;

  fITSnSigma = fpidResponse->NumberOfSigmasITS(Vtrack, AliPID::kElectron);
  fTPCnSigma = fpidResponse->NumberOfSigmasTPC(Vtrack, AliPID::kElectron);
  fTOFnSigma = fpidResponse->NumberOfSigmasTOF(Vtrack, AliPID::kElectron);
  

  // if (TMath::Abs(fITSnSigma)>10) cout << "ITSNSigma" << fITSnSigma << endl;
  
  if (fITSnSigma>fSigmaITScut) PassITSCut=kFALSE; //(fITSnSigma < (-1.*fSigmaITScut) ) || (fITSnSigma > fSigmaITScut)) PassITSCut=kFALSE;
  if ((fTOFnSigma < (-1.*fSigmaTOFcut) ) || (fTOFnSigma > fSigmaTOFcut)) PassTOFCut=kFALSE;
  if ((fTPCnSigma<fSigmaTPCcutLow)  || (fTPCnSigma>fSigmaTPCcutHigh)) PassTPCCut=kFALSE; 


  if (fillHists) { //Fill PID histograms only in first round
    fHistITSnSig->Fill(Vtrack->P(),fITSnSigma, EventWeight);
    fHistTOFnSig->Fill(Vtrack->P(),fTOFnSigma, EventWeight);
    fHistTPCnSig->Fill(Vtrack->P(),fTPCnSigma, EventWeight);
    
    if (PassITSCut) fHistTPCnSigITScut->Fill(Vtrack->P(),fTPCnSigma, EventWeight);
    if (PassTOFCut) fHistTPCnSigTOFcut->Fill(Vtrack->P(),fTPCnSigma, EventWeight);
    if (PassITSCut && PassTOFCut) fHistTPCnSigITSTOFcut->Fill(Vtrack->P(),fTPCnSigma, EventWeight);
    if (PassTOFCut && PassTPCCut) fHistITSnSigTOFTPCcut->Fill(Vtrack->P(),fITSnSigma, EventWeight);

    if (PassTOFCut && fHadCont && fOneTimeCheck) fHadContPvsPt->Fill(Vtrack->P(), Vtrack->Pt(), EventWeight);

    Double_t fillSparse[4];   
    fillSparse[0]=Vtrack->P();
    fillSparse[1]=Vtrack->Phi();
    fillSparse[2]=Vtrack->Eta();
    fillSparse[3]=fTPCnSigma;

    if (PassTOFCut && fHadCont && fOneTimeCheck) fHadContPPhiEtaTPC->Fill(fillSparse, EventWeight);

    fillSparse[1]=fITSnSigma;
    fillSparse[2]=fTOFnSigma;
    if (fHadCont) fHadContamination->Fill(fillSparse, EventWeight);

    fillSparse[0]=Vtrack->Pt();
    //    if (fHadCont) fHadContaminationPt->Fill(fillSparse, EventWeight);

    if (fIsAOD && PassTOFCut && fIsMC && fHadCont) {

      Int_t MClabel=AODtrack->GetLabel();
      AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(MClabel));     
      Int_t PDGCode = abs(MCParticle->GetPdgCode());
      fillSparse[0]=AODtrack->P();
      if (PDGCode==11) fillSparse[1]=0; // electron
      else if (PDGCode==13) fillSparse[1]=1; // muon
      else if (PDGCode==111 || PDGCode==211) fillSparse[1]=2; // p0 and p+
      else if (PDGCode==130 || PDGCode==310 || PDGCode==311 || PDGCode==321) fillSparse[1]=3; // K0L KOS K0 K0+
      else if (PDGCode==2212) fillSparse[1]=4; // proton
      else if (PDGCode==1000010020) fillSparse[1]=5;
      else fillSparse[1]=6;
      fillSparse[2]=fITSnSigma;
      fHadContMC->Fill(fillSparse, EventWeight);
	   // fillSparse[0]=AODtrack->Pt();
      //      fHadContMCPt->Fill(fillSparse);
    }
   }
 
  if (fSigmaITScut>5) PassITSCut = kTRUE; // currently no ITS cut applied (pp)
  if (!PassITSCut || !PassTOFCut || !PassTPCCut) return kFALSE;
  return kTRUE;
}




// _________________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::PhotElecPIDCuts(AliVTrack* Vtrack, Double_t EventWeight) {
  // associated particle variables
  Double_t fITSnSigmaAsso=-9.,fTPCnSigmaAsso=-9.;
  
  if (fIsAOD) {
    AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);   
    Bool_t FBTPCOnly =  AODtrack->TestFilterMask(AliAODTrack::kTrkTPCOnly);
    Bool_t FBITSsa = AODtrack->TestFilterMask(AliAODTrack::kTrkITSsa);
    Bool_t FBITSConst =  AODtrack->TestFilterMask(AliAODTrack::kTrkITSConstrained);


   // looser PID cuts
   fITSnSigmaAsso = fpidResponse->NumberOfSigmasITS(Vtrack, AliPID::kElectron);
   fTPCnSigmaAsso = fpidResponse->NumberOfSigmasTPC(Vtrack, AliPID::kElectron);

   
   if(FBTPCOnly && TMath::Abs(fTPCnSigmaAsso)>fPhotElecSigmaTPCcut) return kFALSE;
   if(FBITSsa && TMath::Abs(fITSnSigmaAsso)>3) return kFALSE;
  }
  return kTRUE;
  
}


//_________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::PhotElecTrackCuts(const AliVVertex *pVtx,AliVTrack *Vtrack,Int_t nMother, Int_t listMother[], Double_t EventWeight)
{
  if (fIsAOD){
    AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);   
    // quality track cuts for associate tracks of photonic electrons
    Bool_t FBTPCOnly  = AODtrack->TestFilterMask(AliAODTrack::kTrkTPCOnly);
    Bool_t FBITSsa    = AODtrack->TestFilterMask(AliAODTrack::kTrkITSsa);
    Bool_t FBITSConst = AODtrack->TestFilterMask(AliAODTrack::kTrkITSConstrained);
    //  if (FBTPCOnly && FBITSsa) cout << "ERRRRRRRRRRRRRORR" << endl;
     // if (!AODtrack->IsOn(AliAODTrack::kITSrefit) && FBITSsa) cout << "No REFIT" << endl;

    if(TMath::Abs(Vtrack->Eta())>0.9) return kFALSE;
    if(Vtrack->Pt() < fPhotElecPtCut) return kFALSE;
    // if (fPhotElecITSrefitCut && !(Vtrack->GetStatus()&AliESDtrack::kITSrefit)) return kFALSE;
    // if (fPhotElecITSrefitCut && !(Vtrack->GetStatus()&AliESDtrack::kTPCrefit)) return kFALSE;


    Int_t ITSNcls;
    ITSNcls = AODtrack->GetITSNcls();
   
    //* ITS shared cluster cut - not useful for single electron studies
    Int_t ITSNclsShared = 0;
    for(int i=0; i<6; i++){
      if( Vtrack->HasSharedPointOnITSLayer(i) ) ITSNclsShared++;
    }
    Float_t ITSFracSharedNCls = ITSNclsShared/(ITSNcls*1.);
    Float_t ITSRedChi2 = Vtrack->GetITSchi2()/(ITSNcls*1);
    


    //  Int_t ParticleID=-99;
    if (fIsMC && FBITSsa) {
      Double_t fillSparse[5]={Vtrack->Pt(), 1.*ITSNcls, 1.*ITSNclsShared, 1.*ITSRedChi2, 1.*IsPhotonicElectron(AODtrack->GetLabel())};
      //fPhotTrackITSCuts->Fill(fillSparse, EventWeight);
    }

    if(!(FBTPCOnly || (FBITSsa && fUseITSsa)) ) return kFALSE;
    if(FBITSsa) {
      if (AODtrack->GetITSNcls()<4) return kFALSE;
      if (fPhotElecITSrefitCut && !(Vtrack->GetStatus()&AliESDtrack::kITSrefit)) return kFALSE;
    }
    if(FBTPCOnly) {
      if (Vtrack->GetTPCNcls() < fPhotElecTPCnCut) return kFALSE;
      if (fPhotElecITSrefitCut && !(Vtrack->GetStatus()&AliESDtrack::kITSrefit)) return kFALSE;
      if (fPhotElecITSrefitCut && !(Vtrack->GetStatus()&AliESDtrack::kTPCrefit)) return kFALSE;
    }

    //    if (FBITSsa) cout << "ITSStandalone" << endl;
    
    return kTRUE;
  }
  return kFALSE; // if not AOD
}

 //_________________________________________

Bool_t AliAnalysisTaskHaHFECorrel::ChargedHadronTrackCuts(const AliVVertex *pVtx,AliVTrack *Vtrack,Int_t nMother, Int_t listMother[], Double_t EventWeight, Bool_t fillHists)
 {

   if (fillHists) fHadronTrackCuts->Fill(Vtrack->Pt(), 0., EventWeight);
   AliAODTrack *AODtrack=0;
   AliESDtrack *ESDtrack=0;
   if (fIsAOD) {
     AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);   
     if (!AODtrack) return kFALSE;
     if(!AODtrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) return kFALSE; 
   }
   else {
     ESDtrack = dynamic_cast<AliESDtrack*>(Vtrack);   
     if(!ESDtrack) return kFALSE;
     if(!ESDkTrkGlobalNoDCA(Vtrack)) return kFALSE;
   }
   if (fillHists) fHadronTrackCuts->Fill(Vtrack->Pt(), 1, EventWeight);

   if(Vtrack->Pt()<0.25) return kFALSE;
   if (fillHists)  fHadronTrackCuts->Fill(Vtrack->Pt(), 2, EventWeight);
    
   if(Vtrack->Eta()>fMaxHadronEta) return kFALSE;
   if(Vtrack->Eta()<fMinHadronEta) return kFALSE;
   if (fillHists) fHadronTrackCuts->Fill(Vtrack->Pt(), 3, EventWeight);

   Bool_t ITSPointOnLayer0=kFALSE;
   Bool_t ITSPointOnLayer1=kFALSE;
   ITSPointOnLayer0=Vtrack->HasPointOnITSLayer(0);
   ITSPointOnLayer1=Vtrack->HasPointOnITSLayer(1);
  

    // quality track cuts for charged hadrons

   Double_t d0z0[2]={-999,-999}, cov[3];
   Float_t DCAr=99, DCAz=99;
   if (fIsAOD) {
     Bool_t kinkmotherpass = kTRUE;
     for(Int_t kinkmother = 0; kinkmother < nMother; kinkmother++) {
       if(AODtrack->GetID() == listMother[kinkmother]) {
	 kinkmotherpass = kFALSE;
	 continue;
       }
     }
     if(!kinkmotherpass) return kFALSE;

     AODtrack->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov);
     DCAr=d0z0[0];
     DCAz=d0z0[1];
   }
   else {
     if(ESDtrack->GetKinkIndex(0) != 0) return kFALSE;
     ESDtrack->GetImpactParameters(DCAr, DCAz);
   }
   if (fillHists) fHadronTrackCuts->Fill(Vtrack->Pt(), 4, EventWeight);
   //   cout << DCAr << "\t" << DCAz << endl;

   if (fillHists) fHadronTrackTPCNcls->Fill(Vtrack->Pt(), Vtrack->GetTPCNcls(), EventWeight);
   if (fillHists) fHadronTrackRefit->Fill(Vtrack->Pt(), Vtrack->IsOn(AliAODTrack::kTPCrefit), Vtrack->IsOn(AliAODTrack::kITSrefit), EventWeight);
   if (fillHists) fHadronTrackDCA->Fill(Vtrack->Pt(), DCAr, DCAz, EventWeight);


   if(Vtrack->GetTPCNcls() < fHTPCnCut) return kFALSE;
   if (fillHists) fHadronTrackCuts->Fill(Vtrack->Pt(), 5, EventWeight);
   
   if(fHTPCrefitCut && !(Vtrack->GetStatus()&AliESDtrack::kTPCrefit)) return kFALSE;
   if (fillHists) fHadronTrackCuts->Fill(Vtrack->Pt(), 6, EventWeight);
   
   if(fHITSrefitCut && !(Vtrack->GetStatus()&AliESDtrack::kITSrefit)) return kFALSE;
   if (fillHists) fHadronTrackCuts->Fill(Vtrack->Pt(), 7, EventWeight);
   
   if(TMath::Abs(DCAr) > fHadDCAr || TMath::Abs(DCAz) > fHadDCAz) {
     //   cout << "Hadron fails DCA cut " << Vtrack->Pt() << endl;
     return kFALSE;
   }
   if (fillHists) fHadronTrackCuts->Fill(Vtrack->Pt(), 8, EventWeight);

   if (fillHists) fHadronTrackDCA_woITSAny->Fill(Vtrack->Pt(), DCAr, DCAz, EventWeight);
   
   if (fHadkAny && !(ITSPointOnLayer0 || ITSPointOnLayer1)) return kFALSE; // now kAny
   if (fillHists) fHadronTrackCuts->Fill(Vtrack->Pt(), 9, EventWeight);

   if (fillHists) fHadronTrackDCA_wITSAny->Fill(Vtrack->Pt(), DCAr, DCAz, EventWeight);

   if (fHadTOFmatch & !(Vtrack->GetTOFBunchCrossing()==0)) return kFALSE;
   if (fillHists) fHadronTrackCuts->Fill(Vtrack->Pt(), 10, EventWeight);
   
   return kTRUE;   
}


//_________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::ChargedHadronPIDCuts(AliVTrack *Vtrack, Double_t EventWeight)
{   
  return kTRUE;   
}


Bool_t AliAnalysisTaskHaHFECorrel::AssoHadronPIDCuts(AliVTrack *Vtrack, Double_t EventWeight) 
{
  // associated particle variables
  Double_t fTPCnSigmaAsso=-9.;
  Double_t fTOFnSigmaAsso=-9.;
  


  fTOFnSigmaAsso = fpidResponse->NumberOfSigmasTOF(Vtrack, AliPID::kElectron);
   if ((fTOFnSigmaAsso < (-1.*fSigmaTOFcut) ) || (fTOFnSigmaAsso > fSigmaTOFcut)) return kFALSE;
  // looser PID cuts
  fTPCnSigmaAsso = fpidResponse->NumberOfSigmasTPC(Vtrack, AliPID::kElectron);
  if(TMath::Abs(fTPCnSigmaAsso)<3.5 || TMath::Abs(fTPCnSigmaAsso)>4 ) return kFALSE; //fAssNonEleTPCcut) return kFALSE;
   
  return kTRUE;
}

void AliAnalysisTaskHaHFECorrel::MCEfficiencyCorrections(const AliVVertex * RecVertex, Double_t EventWeight) {
  for (Int_t i=0; i<fMC->GetNumberOfTracks(); i++) {
    AliAODMCParticle *mcPart  = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(i));
    Double_t mcPt, mcPhi, mcEta, mcVtx[3];
    Int_t mcPDG, mcMotherID;
    mcPt=mcPart->Pt();
    mcPhi=mcPart->Phi();
    mcEta=mcPart->Eta();
    mcPDG=abs(mcPart->GetPdgCode());
    mcMotherID=mcPart->GetMother();
    
    fMCheader->GetVertex(mcVtx);
    fCheckMCVertex->Fill(mcVtx[2], RecVertex->GetZ()); // TH2D

    if (fPionEtaProduction) {

      // PI0 and Eta
      if (mcPart->IsPrimary()) {
	TString genname;
	Bool_t GetCocktailHeader=fMC->GetCocktailGenerator(i,genname);
	if(!GetCocktailHeader) Printf("no cocktail header list was found for this event");
	Double_t fillSparse[5];
	if(GetCocktailHeader) {
	  //  Printf("cocktail header name is %s for particle %i", genname.Data(), mcPDG); 
	  if(genname.Contains("pi0") || genname.Contains("eta")) {
	    fillSparse[4]=1;
	  } 
	  else fillSparse[4]=0;
	}
	else fillSparse[4]=2;

	// Pi0 = 111.      eta             pi+            r0         r+                eta'          w             phi0  
	if (mcPDG==111 || mcPDG==221 || mcPDG== 211 || mcPDG==113 || mcPDG ==213 || mcPDG==331 || mcPDG==223 || mcPDG==333 || mcPDG==130|| mcPDG==310 || mcPDG==311 || mcPDG==321||  mcPDG==443 || mcPDG==553) {
	  if (mcMotherID<=0) {
	    //  fillSparse[4]=1;
	    fillSparse[3]=PDGMap.find(-1)->second;
	  }
	  else {
	    AliAODMCParticle* mcPartMother=(AliAODMCParticle*)fMC->GetTrack(mcMotherID);
	    Int_t mcMotherPDG = abs(mcPartMother->GetPdgCode());
	    // fillSparse[4]=0;
	    fillSparse[3]=PDGMap.find(mcMotherPDG)->second;
	  }
	  fillSparse[0]=mcPt;
	  fillSparse[1]=mcEta;
	  fillSparse[2]=Eta2y(mcPt, mcPart->M(), mcEta);
	
	  if (mcPDG ==111) fMCPi0Prod->Fill(fillSparse);
	  else if (mcPDG==221) fMCEtaProd->Fill(fillSparse);
	  else if (mcPDG==211) fMCPiPlusProd->Fill(fillSparse);
	  else{
	    // if (mcPart->IsSecondaryFromWeakDecay()) cout<< "SecFromWeakDecay" << mcPDG << "\t" << endl; // mcMotherPDG << endl;
	    // else cout << "primary" << endl;
	    fillSparse[4]=fillSparse[3]; //Mother
	    fillSparse[3]=PDGMap.find(mcPDG)->second;
	    fMCBGProd->Fill(fillSparse);
	  }
	}
      }


      if (mcPart->IsPhysicalPrimary() && mcPDG==211) {
	if (!mcPart->IsSecondaryFromMaterial() && ! mcPart->IsSecondaryFromWeakDecay()) {
	  Double_t fillSparse[5];
	  TString genname;
	  Bool_t GetCocktailHeader=fMC->GetCocktailGenerator(i,genname);
	  if(!GetCocktailHeader) Printf("no cocktail header list was found for this event");
	  if(GetCocktailHeader) {
	    //    Printf("cocktail header name is %s for particle %i", genname.Data(), mcPDG); 
	    if(genname.Contains("pi0") || genname.Contains("eta")) {
	      fillSparse[4]=1;
	    } 
	    else fillSparse[4]=0;
	  }
	  else fillSparse[4]=2;

	  if (mcMotherID<=0) {
	    //   fillSparse[4]=1;
	    fillSparse[3]=PDGMap.find(-1)->second;
	  }
	  else {
	    AliAODMCParticle* mcPartMother=(AliAODMCParticle*)fMC->GetTrack(mcMotherID);
	    Int_t mcMotherPDG = abs(mcPartMother->GetPdgCode());
	    //  fillSparse[4]=0;
	    fillSparse[3]=PDGMap.find(mcMotherPDG)->second;
	  }
	  fillSparse[0]=mcPt;
	  fillSparse[1]=mcEta;
	  fillSparse[2]=Eta2y(mcPt, mcPart->M(), mcEta);
	  fMCPiPlusProdV2->Fill(fillSparse);
	}
      }
    }
    
   
    // select true generate hadrons and electrons for rec efficiency
    if (mcPart->Charge()==0) continue;    
    
    if (mcPart->IsPhysicalPrimary()) {
      Double_t fillSparse[4];
      fillSparse[0]=mcPt;
      fillSparse[1]=mcEta;
      fillSparse[2]=mcPhi;
      fillSparse[3]=mcVtx[2];
      
      if (mcEta> fMinHadronEta && mcEta < fMaxHadronEta) {
	fMCHadPtEtaPhiVtx->Fill(fillSparse, EventWeight);
      }
      
      if (abs(mcPDG)==11) {
	if (mcMotherID>0) {
	  AliAODMCParticle* mcPartMother=(AliAODMCParticle*)fMC->GetTrack(mcMotherID);
	  Int_t mcMotherPDG = abs(mcPartMother->GetPdgCode()); 

	  if (mcMotherPDG==11 || mcMotherPDG==15) { // to includ HF->e->e HF->tau->e
	    AliAODMCParticle* mcPartGMother=(AliAODMCParticle*)fMC->GetTrack(mcPartMother->GetMother());
	    mcMotherPDG = abs(mcPartGMother->GetPdgCode());
	  }




	  Int_t MIsHeavy  = Int_t (mcMotherPDG / TMath::Power(10, Int_t(TMath::Log10(mcMotherPDG))));
	  if (MIsHeavy>3 && MIsHeavy<6 && mcEta> fMinElectronEta && mcEta < fMaxElectronEta) {
	    fillSparse[2]=mcPhi;
	    fMCElecPtEtaPhiVtx->Fill(fillSparse, EventWeight);
	    if (fOneTimeCheck && fRecEff) {
	      fMCElecPDG->Fill(mcMotherPDG, EventWeight);
	      if (mcMotherPDG==411 || mcMotherPDG==421 || mcMotherPDG==431 || mcMotherPDG==511 || mcMotherPDG ==521 || mcMotherPDG ==531 || mcMotherPDG==541) {
		//	fMCElecPtEtaPhiStrictVtx->Fill(fillSparse);
	      }
	    } 
	  }
	}     
      }
    }
  }
}



void AliAnalysisTaskHaHFECorrel::EvaluateTaggingEfficiency(AliVTrack * Vtrack, Int_t LSPartner, Int_t ULSPartner, Bool_t  trueULSPartner, Double_t EventWeight, Double_t mult, Double_t recEffE) {
  if (fIsAOD) {
    AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);   
    if (!AODtrack) return;
    Int_t MClabel=AODtrack->GetLabel();
    AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(MClabel));     
    Int_t PDGCode = abs(MCParticle->GetPdgCode());
    Int_t PDGCodeMother=0;
    Bool_t MotherIsPrimary=0;
    Int_t PDGCodeGrandMother=0;
    Bool_t GrandMotherIsPrimary=0;
    Int_t PDGCodeGGMother=0;
    Double_t PtMother=0; 
    Double_t PtGrandMother=0;


    //  cout << MClabel << "\t" << AODtrack->Pt() << "\t" << recEffE << endl;
    
    // get ancestors (eta-pi0-e, eta-gamma-e, pi0-gamma-e, eta-e, pi0-e, ...) - remember pt of pi0/Eta for correction
    if (MCParticle->GetMother()>=0) {
      AliAODMCParticle* MCParticleMother = (AliAODMCParticle*) fMC->GetTrack(MCParticle->GetMother());
      PDGCodeMother = abs(MCParticleMother->GetPdgCode());
      MotherIsPrimary = MCParticleMother->IsPrimary();
      PtMother = MCParticleMother->Pt();
      if (MCParticleMother->GetMother()>=0) {
	AliAODMCParticle* MCParticleGrandMother=(AliAODMCParticle*) fMC->GetTrack(MCParticleMother->GetMother());
	PDGCodeGrandMother=abs(MCParticleGrandMother->GetPdgCode());
	GrandMotherIsPrimary = MCParticleGrandMother->IsPrimary();
	PtGrandMother = MCParticleGrandMother->Pt();
	if (MCParticleGrandMother->GetMother()>=0) {
	  AliAODMCParticle* MCParticleGGMother = (AliAODMCParticle*) fMC->GetTrack(MCParticleGrandMother->GetMother());
	  PDGCodeGGMother=abs(MCParticleGGMother->GetPdgCode());
	}
	else{
	  //	cout << "No Mother? " << MCParticleGrandMother->GetMother() << endl;
	  PDGCodeGGMother=-1;
	}
      }
      else {
	PDGCodeGrandMother=-1;
	PDGCodeGGMother=-1;
      } 
    }
    else{
      PDGCodeMother=-1;
      PDGCodeGrandMother=-1;
      PDGCodeGGMother=-1;
    } 
  
 

    Double_t fillSparse[6]={-999,-999,-999, -999, -999, -999};
    fillSparse[0]=Vtrack->Pt(); // sparse will be filled with rec pt 
    if (PDGCode==11) fillSparse[1]=1.; //PDGMap.find(PDGCode)->second;
    else (fillSparse[1]=0.);
    fillSparse[2]=PDGMap.find(PDGCodeMother)->second;
    fillSparse[3]=PDGMap.find(PDGCodeGrandMother)->second; 
    if (PDGCodeGGMother<400) fillSparse[4]=1;
    else fillSparse[4]=0.; //PDGMap.find(PDGCodeGGMother)->second;
    Bool_t IsPrimary=kFALSE;
 
    Double_t PtMotherWeight=1.; 
    // Double_t ExcludePion[9] = {130,211, 310, 321, 2212, 3122,3222,3322,5122}; // K0L, PI+, K0S, K+, p, Lambd0, Sigma+, Xci0, lambda b0
    // Double_t ExcludeEta[0] = {};

    Int_t MotherIsHeavy =  Int_t (PDGCodeMother / TMath::Power(10, Int_t(TMath::Log10(PDGCodeMother))));
    Int_t GrandMotherIsHeavy =  Int_t (PDGCodeGrandMother / TMath::Power(10, Int_t(TMath::Log10(PDGCodeGrandMother))));
    Int_t GGMotherIsHeavy =  Int_t (PDGCodeGGMother / TMath::Power(10, Int_t(TMath::Log10(PDGCodeGGMother))));
    
    if (PDGCode==11) { // only correct contributions which have been used for correction (primaries
      if (PDGCodeMother==22 && GrandMotherIsPrimary) {
	IsPrimary=kTRUE;
	if (PDGCodeGrandMother == 221) {
	  PtMotherWeight=GetEtaWeight(PtGrandMother); 
	}
	else if (PDGCodeGrandMother == 111 ) {
	  PtMotherWeight=GetPionWeight(PtGrandMother); 
	  //	for (Int_t j=0; j<9; j++) if (PDGCodeGGMother==ExcludePion[j]) PtMotherWeight=1.;
	}
	else if (GrandMotherIsHeavy<4 &&  GGMotherIsHeavy<4) { //  exclude HF but include quarkonia
	  PtMotherWeight=GetBackgroundWeight(PDGCodeGrandMother, PtGrandMother);
	  //	  CheckParticleOrigin(MClabel);
	}
      }
      else if (MotherIsPrimary) {
	IsPrimary=kTRUE;
	if (PDGCodeMother == 221) {
	  PtMotherWeight=GetEtaWeight(PtMother); 
	}
	else if (PDGCodeMother == 111) {
	  PtMotherWeight=GetPionWeight(PtMother); 
	  //for (Int_t j=0; j<9; j++) if (PDGCodeGGMother==ExcludePion[j]) PtMotherWeight=1.;
	}
	else if (MotherIsHeavy<4 &&  GrandMotherIsHeavy<4 && GGMotherIsHeavy<4) {
	  PtMotherWeight=GetBackgroundWeight(PDGCodeMother, PtMother);
	  //	  CheckParticleOrigin(MClabel);
	}
      }
    }
    fillSparse[5]=IsPrimary;
    //    if (!IsPrimary) cout << mult << "\t" << PDGCode << "\t" << PDGCodeMother << "\t" << PDGCodeGrandMother << "\t" << PDGCodeGGMother << endl;


    fTagEffIncl->Fill(fillSparse, PtMotherWeight*EventWeight);
    for (int j=0; j<LSPartner; j++)  fTagEffULS->Fill(fillSparse, -1*PtMotherWeight*EventWeight);
    for (int j=0; j<ULSPartner; j++) fTagEffULS->Fill(fillSparse, PtMotherWeight*EventWeight);
    if (trueULSPartner) fTagTruePairs->Fill(fillSparse, PtMotherWeight*EventWeight);

    if (fOneTimeCheck) { // check tag eff wo MC reweighting
      // fTagEffInclWoWeight->Fill(fillSparse, EventWeight);
      // for (int j=0; j<LSPartner; j++)  fTagEffULSWoWeight->Fill(fillSparse, -1.*EventWeight);
      // for (int j=0; j<ULSPartner; j++) fTagEffULSWoWeight->Fill(fillSparse, EventWeight);
      // if (trueULSPartner) fTagTruePairsWoWeight->Fill(fillSparse, EventWeight);
    }

    // Reduced TagEff for MultMeasurement
    Double_t PDGBG=0;
    if (PDGCodeMother ==22) {
      if (PDGCodeGrandMother == 113) PDGBG=5; // rho
      else if (PDGCodeGrandMother == 223) PDGBG=4; // omega
      else if (PDGCodeGrandMother == 333) PDGBG=6; // phi
      else if (PDGCodeGrandMother == 130) PDGBG=2; // K0L
      else if (PDGCodeGrandMother == 310) PDGBG=3; // K0S
      else if (PDGCodeGrandMother == 321) PDGBG=1; // K+-
    }
    else {
      if (PDGCodeMother == 113) PDGBG=5; // rho
      else if (PDGCodeMother == 223) PDGBG=4; // omega
      else if (PDGCodeMother == 333) PDGBG=6; // phi
      else if (PDGCodeMother == 130) PDGBG=2; // K0L
      else if (PDGCodeMother == 310) PDGBG=3; // K0S
      else if (PDGCodeMother == 321) PDGBG=1; // K+-
    }








    Double_t pt=Vtrack->Pt(); 
    if (PDGCode==11 && IsPrimary) {     
      if (PDGCodeMother==111 || PDGCodeMother==221) {
	fTagEffInclMult->Fill(pt, mult, PtMotherWeight*EventWeight);
        for (int j=0; j<LSPartner; j++)   fTagEffULSMult->Fill(pt, mult, -1.*PtMotherWeight*EventWeight);
	for (int j=0; j<ULSPartner; j++)   fTagEffULSMult->Fill(pt, mult, PtMotherWeight*EventWeight);
	if (trueULSPartner) fTagTruePairsMult->Fill(pt, mult,  PtMotherWeight*EventWeight);
	if (trueULSPartner && !(LSPartner>0 || ULSPartner>0)) cout << "ERROR ULSLS" << endl;
      }
      else if (PDGCodeMother==22 &&( PDGCodeGrandMother==111 || PDGCodeGrandMother==221)) {
	fTagEffInclMult->Fill(pt, mult, PtMotherWeight*EventWeight);
	for (int j=0; j<LSPartner; j++)  fTagEffULSMult->Fill(pt,mult, -1.*PtMotherWeight*EventWeight);
	for (int j=0; j<ULSPartner; j++) fTagEffULSMult->Fill(pt,mult, PtMotherWeight*EventWeight);
	if (trueULSPartner) fTagTruePairsMult->Fill(pt,mult, PtMotherWeight*EventWeight);
	if (trueULSPartner && !(LSPartner>0 || ULSPartner>0)) cout << "ERROR ULSLS" << endl;
      }
      else if ((MotherIsHeavy<4 &&  GrandMotherIsHeavy<4 && GGMotherIsHeavy<4)  && ((PDGCodeMother>100) || (PDGCodeMother==22 && PDGCodeGrandMother>100))) {
	fTagEffInclBGMult->Fill(pt, mult, PDGBG, PtMotherWeight*EventWeight/recEffE);
	for (int j=0; j<LSPartner; j++)  fTagEffULSBGMult->Fill(pt,mult, PDGBG, -1.*PtMotherWeight*EventWeight/recEffE);
	for (int j=0; j<ULSPartner; j++) fTagEffULSBGMult->Fill(pt,mult, PDGBG, PtMotherWeight*EventWeight/recEffE);
      }
      else {
	//	cout << "primary, but neither photonic or default bg" << PDGCodeMother <<  endl; // usually heavy particles
      }
    }
    else {
      //   cout << "not primary" <<endl;
    }
	

    if (kTRUE) { // same as above but without pt mother weights for MC closure     
      if (PDGCode==11 && IsPrimary) {     
	if (PDGCodeMother==111 || PDGCodeMother==221) {
	  fTagEffInclMultWoW->Fill(pt, mult, EventWeight);
	  for (int j=0; j<LSPartner; j++)   fTagEffULSMultWoW->Fill(pt, mult, -1.*EventWeight);
	  for (int j=0; j<ULSPartner; j++)   fTagEffULSMultWoW->Fill(pt, mult, EventWeight);
	  if (trueULSPartner) fTagTruePairsMultWoW->Fill(pt, mult,  EventWeight);
	  if (trueULSPartner && !(LSPartner>0 || ULSPartner>0)) cout << "ERROR ULSLS" << endl;
	}
	else if (PDGCodeMother==22 &&( PDGCodeGrandMother==111 || PDGCodeGrandMother==221)) {
	  fTagEffInclMultWoW->Fill(pt, mult, EventWeight);
	  for (int j=0; j<LSPartner; j++)  fTagEffULSMultWoW->Fill(pt,mult, -1.*EventWeight);
	  for (int j=0; j<ULSPartner; j++) fTagEffULSMultWoW->Fill(pt,mult, EventWeight);
	  if (trueULSPartner) fTagTruePairsMultWoW->Fill(pt,mult, EventWeight);
	  if (trueULSPartner && !(LSPartner>0 || ULSPartner>0)) cout << "ERROR ULSLS" << endl;
	}
	else if ((MotherIsHeavy<4 &&  GrandMotherIsHeavy<4 && GGMotherIsHeavy<4)  && ((PDGCodeMother>100) || (PDGCodeMother==22 && PDGCodeGrandMother>100))) {
	  fTagEffInclBGMultWoW->Fill(pt, mult, PDGBG, EventWeight/recEffE);
	  for (int j=0; j<LSPartner; j++)  fTagEffULSBGMultWoW->Fill(pt,mult, PDGBG, -1.*EventWeight/recEffE);
	  for (int j=0; j<ULSPartner; j++) fTagEffULSBGMultWoW->Fill(pt,mult, PDGBG, EventWeight/recEffE);
	}
	else {
	  //	cout << "primary, but neither photonic or default bg" << PDGCodeMother <<  endl; // usually heavy particles
	}
      }
      else {
	//   cout << "not primary" <<endl;
      }



      if (PDGCode==11 && !IsPrimary) {     
	if (PDGCodeMother==111 || PDGCodeMother==221) {
	  fTagEffInclMultWoWS->Fill(pt, mult, EventWeight);
	  for (int j=0; j<LSPartner; j++)   fTagEffULSMultWoWS->Fill(pt, mult, -1.*EventWeight);
	  for (int j=0; j<ULSPartner; j++)   fTagEffULSMultWoWS->Fill(pt, mult, EventWeight);
	  if (trueULSPartner) fTagTruePairsMultWoWS->Fill(pt, mult,  EventWeight);
	  if (trueULSPartner && !(LSPartner>0 || ULSPartner>0)) cout << "ERROR ULSLS" << endl;
	}
	else if (PDGCodeMother==22 &&( PDGCodeGrandMother==111 || PDGCodeGrandMother==221)) {
	  fTagEffInclMultWoWS->Fill(pt, mult, EventWeight);
	  for (int j=0; j<LSPartner; j++)  fTagEffULSMultWoWS->Fill(pt,mult, -1.*EventWeight);
	  for (int j=0; j<ULSPartner; j++) fTagEffULSMultWoWS->Fill(pt,mult, EventWeight);
	  if (trueULSPartner) fTagTruePairsMultWoWS->Fill(pt,mult, EventWeight);
	  if (trueULSPartner && !(LSPartner>0 || ULSPartner>0)) cout << "ERROR ULSLS" << endl;
	}
	else if ((MotherIsHeavy<4 &&  GrandMotherIsHeavy<4 && GGMotherIsHeavy<4)  && ((PDGCodeMother>100) || (PDGCodeMother==22 && PDGCodeGrandMother>100))) {
	  fTagEffInclBGMultWoWS->Fill(pt, mult, PDGBG, EventWeight/recEffE);
	  for (int j=0; j<LSPartner; j++)  fTagEffULSBGMultWoWS->Fill(pt,mult, PDGBG, -1.*EventWeight/recEffE);
	  for (int j=0; j<ULSPartner; j++) fTagEffULSBGMultWoWS->Fill(pt,mult, PDGBG, EventWeight/recEffE);
	}
	else {
	  //	cout << "primary, but neither photonic or default bg" << PDGCodeMother <<  endl; // usually heavy particles
	}
      }
      else {
	//   cout << "not primary" <<endl;
      }
    }



    
    // extra sparse to control pt mother, pt electron, ratio of primary to secondary
    if (PDGCode==11) {
      if (PDGCodeMother==22) {
	fillSparse[1]=PtGrandMother;
	fillSparse[2]=PDGMap.find(PDGCodeGrandMother)->second;
	fillSparse[3]=PDGMap.find(PDGCodeGGMother)->second;
      }
      else {
	fillSparse[1]=PtMother;
	fillSparse[2]=PDGMap.find(PDGCodeMother)->second;
	fillSparse[3]=PDGMap.find(PDGCodeGrandMother)->second;
      }
      fTagMotherPt->Fill(fillSparse, PtMotherWeight*EventWeight);
    }
  }
}




//_________________________________
Bool_t  AliAnalysisTaskHaHFECorrel::CloneAndReduceTrackList(TObjArray* RedTracks, AliVTrack* Vtrack, Int_t LSPartner, Int_t ULSPartner, Int_t* LSPartnerID, Int_t* ULSPartnerID,  Float_t* LSPartnerWeight, Float_t* ULSPartnerWeight, Bool_t trueULSPartner, Float_t MCPartnerPt, Float_t RecPartnerPt, Bool_t isPhotonic, Bool_t isHadron) {
  
  // Copied form AliAnalysisTaksPhiCorrelations and following instructions on
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PAGCorrelationsFAQ#In_Memory_Event_Mixing
  // Main difference is that track selection is in main to reduce computation time


  fCheckLSULS->Fill(LSPartner, ULSPartner);

  //  cout << "ParticleAdded" << endl;

  // clones a track list by using AliBasicParticle which uses much less memory (used for event mixing) -- modified class for my use
  // Clone only a certain pt bin on demand
  AliBasicParticleHaHFE * bparticle  = 0;
  Int_t label=-999;
  if (fIsMC && fIsAOD) {
    AliAODTrack * AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);
    if (AODtrack) label=AODtrack->GetLabel();
  }
  else if (fIsMC && !fIsAOD) {
      AliESDtrack* ESDtrack = dynamic_cast<AliESDtrack*>(Vtrack);
    if (ESDtrack) label=ESDtrack->GetLabel();
  }


 AliExternalTrackParam extTrackParam;
 extTrackParam.CopyFromVTrack(Vtrack);

  bparticle = new AliBasicParticleHaHFE(Vtrack->GetID(), Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt(), Vtrack->Charge(), LSPartner, ULSPartner, LSPartnerID, ULSPartnerID, LSPartnerWeight, ULSPartnerWeight, trueULSPartner, MCPartnerPt, RecPartnerPt, isPhotonic, isHadron, abs(label), extTrackParam);
  if (!bparticle) return kFALSE;
  
  RedTracks->Add(bparticle);
  return kTRUE;

}

Double_t AliAnalysisTaskHaHFECorrel::Eta2y(Double_t pt, Double_t m, Double_t eta) const
{
  // convert eta to y
  Double_t mt = TMath::Sqrt(m * m + pt * pt);
  return TMath::ASinH(pt / mt * TMath::SinH(eta));
}

void AliAnalysisTaskHaHFECorrel::BinLogX(TAxis *axis)
{
  // Method for the correct logarithmic binning of histograms
  const Int_t bins    = axis->GetNbins();
  const Double_t from = axis->GetXmin();
  const Double_t to   = axis->GetXmax();
    
  if (from<10e-9){
    //printf("BinLogX warning xmin < epsilon! nothing done, axis not set. %e\n", from);  exit(1);
    return;
  }
  Double_t * new_bins = new Double_t[bins + 1];
  new_bins[0] = from;
  const Double_t factor = pow(to/from, 1./bins);
  for (int i = 1; i <= bins; i++) {
    new_bins[i] = factor * new_bins[i-1];
  }
  axis->Set(bins, new_bins);
  delete [] new_bins;
}

void AliAnalysisTaskHaHFECorrel::SetPDGAxis(TAxis *axis, std::vector<TString> PDGLabel)
{
  for (Int_t i=1; i<axis->GetNbins(); i++) {
    if (i<PDGLabel.size()) {
      axis->SetBinLabel(i,PDGLabel[i].Data());
      // cout << PDGLabel[i] << endl;
    }
    else axis->SetBinLabel(i, "-");
  }
  //  cout << "Axis Label" << axis->GetBinLabel(22) << endl;
}

void AliAnalysisTaskHaHFECorrel::SetTriggerAxis(TAxis *axis)
{
  for (Int_t i=1; i<=axis->GetNbins(); i++) {
    TString AxisLabel = Form("%4.2f_%4.2f", fAssPtHad_Xmin[i-1], fAssPtHad_Xmax[i-1]); 
    axis->SetBinLabel(i, AxisLabel.Data());
   }
}

Bool_t AliAnalysisTaskHaHFECorrel::HaveSameMother(Int_t Label1, Int_t Label2) const
{
  AliAODMCParticle* MCParticle1 = (AliAODMCParticle*) fMC->GetTrack(abs(Label1));  
  AliAODMCParticle* MCParticle2 = (AliAODMCParticle*) fMC->GetTrack(abs(Label2));
  if (MCParticle1->GetMother() ==  MCParticle2->GetMother()) return kTRUE;
  else {
    Int_t Mother1=(MCParticle1->GetMother());
    Int_t Mother2=(MCParticle2->GetMother());
    AliAODMCParticle* GMother1=(AliAODMCParticle*) fMC->GetTrack(abs(Mother1));
    AliAODMCParticle* GMother2=(AliAODMCParticle*) fMC->GetTrack(abs(Mother2));

    if (Mother1==11 && (GMother1->GetMother() == Mother2)) return kTRUE;
    if (Mother2==11 && (GMother2->GetMother() == Mother1)) return kTRUE;
    if (Mother1==11 && Mother2==11 && (GMother1->GetMother() == GMother2->GetMother())) return kTRUE;
    return kFALSE;
  }
}

Bool_t AliAnalysisTaskHaHFECorrel::IsPhotonicElectron(Int_t Label1) const
{
  Int_t PDGCodeMother=-999, PDGCodeGrandMother=-999, PDGCodeGGMother=-999;
  AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(Label1));
  if (abs(MCParticle->GetPdgCode())!=11) return kFALSE; // false if no electron  
  if (MCParticle->GetMother()>=0) {
     AliAODMCParticle* MCParticleMother = (AliAODMCParticle*) fMC->GetTrack(MCParticle->GetMother());
     PDGCodeMother = abs(MCParticleMother->GetPdgCode());
     if (PDGCodeMother == 111 || PDGCodeMother == 221) { // if mother is pion or eta
       if (MCParticleMother->GetMother()>=0) { // check grandmother<399
	 AliAODMCParticle* MCParticleGrandMother=(AliAODMCParticle*) fMC->GetTrack(MCParticleMother->GetMother());
	 PDGCodeGrandMother=abs(MCParticleGrandMother->GetPdgCode());
	 if (PDGCodeGrandMother<399) return kTRUE; // ignore HF feed-down return kTRUE; // true if mother is eta or pion
	 else return kFALSE;
       }
       else return kTRUE;
     }
     else if (PDGCodeMother == 22 && MCParticleMother->GetMother()>=0) { // if mother is gamma
       AliAODMCParticle* MCParticleGrandMother=(AliAODMCParticle*) fMC->GetTrack(MCParticleMother->GetMother());
       PDGCodeGrandMother=abs(MCParticleGrandMother->GetPdgCode());
       if (PDGCodeGrandMother == 111 || PDGCodeGrandMother == 221) { // if grandmother is pion or eta
	 if (MCParticleGrandMother->GetMother()>=0) {// check if ggmother is<400
	   AliAODMCParticle* MCParticleGGMother=(AliAODMCParticle*) fMC->GetTrack(MCParticleGrandMother->GetMother());
	   PDGCodeGGMother=abs(MCParticleGGMother->GetPdgCode());
	   if (PDGCodeGGMother<399) return kTRUE; // ignore HF feed-down
	   else return kFALSE;
	 }
	 else return kTRUE;
       }
       else return kFALSE;
     }
     else return kFALSE;
  }
  else return kFALSE;

  return kFALSE;
}




void AliAnalysisTaskHaHFECorrel::FindV0CandidatesESD(AliESDEvent *event) 
{

  fV0cutsESD->SetMode(AliESDv0KineCuts::kEffGamma,AliESDv0KineCuts::kPP);
  fV0cutsESD->SetEvent(event);

  const Int_t numTracks = event->GetNumberOfTracks();
  fV0tags = new Int_t[numTracks];
  for (Int_t i = 0; i < numTracks; i++) fV0tags[i] = 0;

  // loop over V0 particles
  for(Int_t iv0=0; iv0<event->GetNumberOfV0s();iv0++){
    AliESDv0 *v0 = (AliESDv0 *) event->GetV0(iv0);
 
    if(!v0) continue;
    if(v0->GetOnFlyStatus()) continue; 
    // Get the particle selection 
    Bool_t foundV0 = kFALSE;
    Int_t pdgV0, pdgP, pdgN;
    foundV0 = fV0cutsESD->ProcessV0(v0, pdgV0, pdgP, pdgN);
    if(!foundV0) continue;
    Int_t iTrackP = v0->GetPindex();  // positive track
    Int_t iTrackN = v0->GetNindex();  // negative track

    // v0 Armenteros plot (QA)
    Float_t armVar[2] = {0.0,0.0};
    fV0cutsESD->Armenteros(v0, armVar);
    fhArmenteros->Fill(armVar[0],armVar[1]);
      
    if( pdgP == -11){
      fV0electrons->Add((AliVTrack*)event->GetTrack(iTrackP));
      // fV0tags[iTrackP] = 11;
    }
    else if( pdgP == 211){
      fV0pions->Add((AliVTrack*)event->GetTrack(iTrackP));
      // fV0tags[iTrackP] = 211;
    }
    else if( pdgP == 2212){
      fV0protons->Add((AliVTrack*)event->GetTrack(iTrackP));
      //fV0tags[iTrackP] = 2212;
    }

    // negative particles
    if( pdgN == 11){
      fV0electrons->Add((AliVTrack*)event->GetTrack(iTrackN));
      // fV0tags[iTrackN] = -11;
    }
    else if( pdgN == -211){
      fV0pions->Add((AliVTrack*)event->GetTrack(iTrackN));
      //        fV0tags[iTrackN] = -211;
    }
    else if( pdgN == -2212){
      fV0protons->Add((AliVTrack*)event->GetTrack(iTrackN));
      //        fV0tags[iTrackN] = -2212;
    }
  }
}



void AliAnalysisTaskHaHFECorrel::FindV0CandidatesAOD(AliAODEvent *event) 
{
  
  fV0cuts->SetMode(AliAODv0KineCuts::kEffGamma,AliAODv0KineCuts::kPP); // set mode (maxEfficiency or Purity kPurity);
  fV0cuts->SetEvent(event);

  const Int_t numTracks = event->GetNumberOfTracks();
  fV0tags = new Int_t[numTracks];
  for (Int_t i = 0; i < numTracks; i++) fV0tags[i] = 0;

  // loop over V0 particles
  for(Int_t iv0=0; iv0<event->GetNumberOfV0s();iv0++){

    AliAODv0 *v0 = (AliAODv0 *) event->GetV0(iv0);
      
    if(!v0) continue;
    if(v0->GetOnFlyStatus()) continue; 
    // Get the particle selection 
    Bool_t foundV0 = kFALSE;
    Int_t pdgV0, pdgP, pdgN;
    foundV0 = fV0cuts->ProcessV0(v0, pdgV0, pdgP, pdgN);
    if(!foundV0) continue;
    Int_t iTrackP = v0->GetPosID();  // positive track
    Int_t iTrackN = v0->GetNegID();  // negative track
      
    if (iTrackP<0 || iTrackN <0) {
      cout << "negative ID" << endl;
      continue;
    }
      
    fhArmenteros->Fill(v0->Alpha(),v0->QtProng());

    // AliAODTrack (V0 Daughters)
    AliAODVertex* vertex = v0->GetSecondaryVtx();
    if (!vertex) {
      Printf("ERROR: Could not retrieve vertex");
      continue;
    }
      
    AliAODTrack *pTrack=0, *nTrack=0;
    Int_t NDaughters = vertex->GetNDaughters(); // should be only two from construction
    for (Int_t i = 0; i<NDaughters; i++) {
      if (((AliAODTrack*)vertex->GetDaughter(i))->GetID()==iTrackP) pTrack= (AliAODTrack*)vertex->GetDaughter(i);
      else if  (((AliAODTrack*)vertex->GetDaughter(i))->GetID()==iTrackN) nTrack= (AliAODTrack*)vertex->GetDaughter(i);
    }
    if (pTrack==0 || nTrack ==0) continue;
    // if (iTrackP<0 || iTrackN<0) continue;
      
    // fill the Object arrays
    // positive particles
    if( pdgP == -11){
      fV0electrons->Add(pTrack);
      //     fV0tags[iTrackP] = 11;
    }
    else if( pdgP == 211){
      fV0pions->Add(pTrack);
      // fV0tags[iTrackP] = 211;
    }
    else if( pdgP == 2212){
      fV0protons->Add(pTrack);
      //	fV0tags[iTrackP] = 2212;
    }
      
    // negative particles
    if( pdgN == 11){
      fV0electrons->Add(nTrack);
      // fV0tags[iTrackN] = -11;
    }
    else if( pdgN == -211){
      fV0pions->Add(nTrack);
      // fV0tags[iTrackN] = -211;
    }
    else if( pdgN == -2212){
      fV0protons->Add(nTrack);
      //fV0tags[iTrackN] = -2212;
    }
  }
}


void AliAnalysisTaskHaHFECorrel::ClearV0PIDList(){

  // Clear the PID object arrays
  fV0electrons->Clear();
  fV0pions->Clear();
  fV0protons->Clear();
  if (fV0tags!=0) delete[] fV0tags;
  fV0tags = 0;
}

void AliAnalysisTaskHaHFECorrel::TRDQA(Int_t RunNumber, const AliVVertex *pVtx, Int_t nMother, Int_t listMother[], Double_t EventWeight) {

  // Fill full Histograms
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();
 
  // Loop over all tracks 
  for(Int_t jTracks = 0; jTracks < ntracks; jTracks++){
          
    AliVParticle* VHtrack = 0x0;
    if(!fUseTender) VHtrack  = fVevent->GetTrack(jTracks);
    if(fUseTender) VHtrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jTracks));
    if (!VHtrack) {
      printf("ERROR: Could not receive associated track %d\n", jTracks);
      continue;
    }
    AliVTrack *Vtrack = dynamic_cast<AliVTrack*>(VHtrack);
    if(!Vtrack) continue;
    AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);   
    if (fIsAOD && !AODtrack) continue;
    AliESDtrack *ESDtrack = dynamic_cast<AliESDtrack*>(Vtrack);   
    if (!fIsAOD && !ESDtrack) continue;


    if (fIsAOD) {
      if (!AODtrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;
    }
    else {
      if(!ESDkTrkGlobalNoDCA(Vtrack)) continue;
    }

    Bool_t passHFETrackCut=kFALSE;
    passHFETrackCut= InclElecTrackCuts(pVtx,Vtrack,nMother,listMother, EventWeight);

    Bool_t passHFEPIDCut=kFALSE;
    if (passHFETrackCut) passHFEPIDCut= InclElecPIDCuts(Vtrack, EventWeight, kFALSE);  // do not fill hists

    Double_t fillSparse[7];
    fillSparse[0]=Vtrack->Pt();
    fillSparse[1]=Vtrack->Eta();
    fillSparse[2]=Vtrack->Phi();
    fillSparse[4]=Vtrack->Charge();
    // fillSparse[5]=RunNumber;
    
    Int_t nTracklets=0;
    for (Int_t i=0; i<6; i++) {
      if (Vtrack->GetTRDmomentum(i)>0.01) { // TRDMomentumCriterion
	if (!fIsAOD && !ESDtrack->IsTRDtrkltChmbGood(i)) {
	  cout << "Chamber not good" << i << "\t" <<  Vtrack->Eta() << "\t" << Vtrack->Phi() << endl;
	  continue;
	}
	nTracklets++;
	fillSparse[3]=i+1; // Layer
	fTRDEtaPhi->Fill(fillSparse);
      }
    }
       
    if (passHFETrackCut)  fillSparse[2]=1; //PID - pass HFEcut
    else fillSparse[2]=2; // PID
    fillSparse[3]=fpidResponse->NumberOfSigmasTPC(Vtrack, AliPID::kElectron);
    fillSparse[4]=nTracklets;
    fillSparse[5]=Vtrack->Charge();
    fillSparse[6]=RunNumber;
    fTRDNTracklets->Fill(fillSparse);
    
    if (Vtrack->Eta()>0.8 || Vtrack->Eta()<-0.8) continue;
    fTRDnTrackRun->Fill(nTracklets,  (Int_t)Vtrack->GetTRDntrackletsPID(), RunNumber);

    fillSparse[1]=nTracklets;
    if (passHFETrackCut)  fillSparse[2]=1;
    else fillSparse[2]=2;
    fillSparse[3]=fpidResponse->NumberOfSigmasTPC(Vtrack, AliPID::kElectron);
    fillSparse[4]=Vtrack->Charge();
    fTRDSpectra->Fill(fillSparse);

    if (fIsMC) {
      Int_t MClabel=-999;
      Int_t PDGCode=-990; 
      if (fIsAOD) {
	if (!AODtrack) continue;
	MClabel=AODtrack->GetLabel();
	AliAODMCParticle* MCParticleAOD = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(abs(MClabel)));  
	if (!MCParticleAOD) continue;
	PDGCode = abs(MCParticleAOD->GetPdgCode());
      }
      else {
	if (!ESDtrack) continue;
	MClabel=ESDtrack->GetLabel();
	fMC->GetTrack(abs(MClabel));
	//cout << MClabel << endl;
	//cout << fMC->GetNumberOfTracks() << endl;
	AliMCParticle* MCParticleESD = dynamic_cast<AliMCParticle*>(fMC->GetTrack(abs(MClabel)));  
	if (!MCParticleESD) continue;
	PDGCode = abs(MCParticleESD->PdgCode());  
	//cout << MClabel << " TrackPt " << Vtrack->Pt() << " MCpt " << MCParticleESD->Pt() << endl;
      }

    
      if (PDGCode == 11) {
	fillSparse[2]=1;
	fTRDMCSpectra->Fill(fillSparse);
      }
      else if (PDGCode == 211) {
	fillSparse[2]=3;
	fTRDMCSpectra->Fill(fillSparse);
      }
      else if (PDGCode == 2212) {
	fillSparse[2]=2;
	fTRDMCSpectra->Fill(fillSparse);
      }
    }
  }



  for (Int_t i=0; i<fV0electrons->GetEntriesFast(); i++) {
     AliAODTrack *track = (AliAODTrack*) fV0electrons->At(i);
     FillV0Histograms(track, 1, RunNumber);
  }
  for (Int_t i=0; i<fV0protons->GetEntriesFast(); i++) {
     AliAODTrack *track = (AliAODTrack*) fV0protons->At(i);
     FillV0Histograms(track, 2, RunNumber);
  }
  for (Int_t i=0; i<fV0pions->GetEntriesFast(); i++) {
     AliAODTrack *track = (AliAODTrack*) fV0pions->At(i);
     FillV0Histograms(track, 3, RunNumber);

  }
}



void AliAnalysisTaskHaHFECorrel::FillV0Histograms(AliVTrack* Vtrack, Int_t Species, Int_t RunNumber) {
   
  AliESDtrack *ESDtrack = dynamic_cast<AliESDtrack*>(Vtrack);
  AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);
  if (fIsAOD) {
    if (!AODtrack) return;
    if (!AODtrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return;
  }
  else {
    if (!ESDtrack) return;
    if(!ESDkTrkGlobalNoDCA(Vtrack)) return;
  }

  Double_t fillSparse[7];
  fillSparse[0]=Vtrack->Pt();
  fillSparse[1]=Vtrack->Eta();
  fillSparse[2]=Species;
  fillSparse[3]=fpidResponse->NumberOfSigmasTPC(Vtrack, AliPID::kElectron);
  fillSparse[5]=Vtrack->Charge();
  fillSparse[6]=RunNumber;
  
  Int_t nTracklets=0;
  for (Int_t i=0; i<7; i++) {
    if (Vtrack->GetTRDmomentum(i)>0.01) {
      if (!fIsAOD && !ESDtrack->IsTRDtrkltChmbGood(i)) {
	continue;
      }
      nTracklets++;
    }
  }
  
  fillSparse[4]=nTracklets;
  fTRDV0NTracklets->Fill(fillSparse);
  
  if (Vtrack->Eta()>0.8 || Vtrack->Eta()<-0.8) return;
  
  fillSparse[1]= nTracklets;
  fillSparse[2]= Species;
  fillSparse[3]= fpidResponse->NumberOfSigmasTPC(Vtrack, AliPID::kElectron);
  fillSparse[4]= Vtrack->Charge();
  fTRDV0Spectra->Fill(fillSparse);
    
}


void AliAnalysisTaskHaHFECorrel::CheckElectronIsTrigger(Double_t ptH, Bool_t *ElectronIsTrigger) {
  for (Int_t i=0; i<fAssPtHad_Nbins; i++) {
       if (ptH>fAssPtHad_Xmin[i] && ptH< fAssPtHad_Xmax[i]) ElectronIsTrigger[i]=kTRUE;
  }
  return;
}
  
void AliAnalysisTaskHaHFECorrel::CheckHadronIsTrigger(Double_t ptE, Bool_t* HadronIsTrigger) {
  for (Int_t i=0; i<fAssPtElec_Nbins; i++) {
    if (ptE>fAssPtElec_Xmin[i] && ptE<fAssPtElec_Xmax[i]) HadronIsTrigger[i]=kTRUE;
  }
  return;
}
 
Double_t AliAnalysisTaskHaHFECorrel::GetEtaWeight(Double_t pt) {
  if (fCorrectEtatoData.IsZombie()) return 1.;
  Int_t Bin = fCorrectEtatoData.FindBin(pt);
  if (fCorrectEtatoData.IsBinUnderflow(Bin) || fCorrectEtatoData.IsBinOverflow(Bin)) {
    return 0;
  }
  return  fCorrectEtatoData.GetBinContent(Bin);
}

Double_t AliAnalysisTaskHaHFECorrel::GetPionWeight(Double_t pt) {
  if (fCorrectPiontoData.IsZombie()) return 1.;
  Int_t Bin = fCorrectPiontoData.FindBin(pt);
  if (fCorrectPiontoData.IsBinUnderflow(Bin) || fCorrectPiontoData.IsBinOverflow(Bin)) {
    return 0;
  }
  return  fCorrectPiontoData.GetBinContent(Bin);
}

Double_t AliAnalysisTaskHaHFECorrel::GetBackgroundWeight(Int_t PDGMother, Double_t pt) {

  if (fBgWeight.IsZombie()) {
    // cout<< "Zombie Histogram" << endl;
    return 1.;
  }
  Int_t Bin = fBgWeight.FindBin(1.*PDGMother, pt);
  if (fBgWeight.IsBinUnderflow(Bin) || fBgWeight.IsBinOverflow(Bin)){
    //  cout << "BGWeightDef1 " << PDGMother << "\t" << 1. << endl;
 
    return 1.;
  }
  //  cout << "BGWeight " << PDGMother << "\t" << fBgWeight.GetBinContent(Bin) << endl;  
  return fBgWeight.GetBinContent(Bin);


}


Double_t AliAnalysisTaskHaHFECorrel::GetHadronRecEff(Int_t run, Double_t pt, Double_t phi, Double_t eta, Double_t zVtx) {
  if (pt<0.25) return -1.;
  
  Int_t Bin = fHadRecEff.FindBin(1.*run,zVtx, pt);
  if (fHadRecEff.IsBinUnderflow(Bin) || fHadRecEff.IsBinOverflow(Bin) ) {
    return -1.;
  }
  Double_t RecEff = fHadRecEff.GetBinContent(Bin);

  //  cout << "HR " << run << ", " << pt <<  ", " << RecEff << endl;

  if (RecEff>0.05) return RecEff;
  else {
    return -1.;
  }
}

Double_t AliAnalysisTaskHaHFECorrel::GetElectronRecEff(Int_t run, Double_t pt, Double_t phi, Double_t eta, Double_t zVtx) {

  //  return 1.;

  if (pt<0.5) return -1;
  Int_t Bin = fEleRecEff.FindBin(1.*run,zVtx, pt);
  //  cout << "EleBin " << Bin << endl;
  if (fEleRecEff.IsBinUnderflow(Bin) || fEleRecEff.IsBinOverflow(Bin)) {
    // cout <<  "ElecRecEff: " << pt << "\t" << eta << "\t" << zVtx << endl;
    return -1.;
  }
  Double_t RecEff= fEleRecEff.GetBinContent(Bin);
  //  cout << "ER " << run << ", " << pt <<  ", " << RecEff << endl;
  if (RecEff>0.05) return RecEff;
  else {
    //cout <<  "ElecRecEff: " << pt << "\t" << eta << "\t" << zVtx << endl;
    return -1.;
  }
}

Double_t AliAnalysisTaskHaHFECorrel::GetTriggerWeight(Int_t run, Double_t minV0, Double_t nTrAcc) {

  if (!fUseEventWeights) return 1.;
  if (minV0>20) return 1.;
  // cout << "XAxis " << fTriggerWeight.GetXaxis()->GetNbins() << endl;
  // cout << "Yaxis" << fTriggerWeight.GetYaxis()->GetNbins() << endl;
  // cout << "Zaxis" << fTriggerWeight.GetZaxis()->GetNbins() << endl;

  
  Int_t BinX= 0;
  BinX = fTriggerWeight.GetXaxis()->FindFixBin(run);
  if (BinX>fTriggerWeight.GetXaxis()->GetNbins()) cout << "TriggerWeight  - Run out of range" << endl;
  Int_t BinY=0;
  BinY = fTriggerWeight.GetYaxis()->FindFixBin(minV0);
  if (BinY>fTriggerWeight.GetYaxis()->GetNbins()) return 1.;

  Int_t BinZ= fTriggerWeight.GetZaxis()->FindFixBin(1.*nTrAcc);
  if (BinZ>fTriggerWeight.GetZaxis()->GetNbins()) BinZ=fTriggerWeight.GetZaxis()->GetLast();
  Double_t TriggerWeight = fTriggerWeight.GetBinContent(BinX, BinY, BinZ);
  // Int_t Bin = fTriggerWeight.FindBin(1.*run, minV0, 1.*nTrAcc);
  // cout << "TriggerWeightNew" << fTriggerWeight.GetBinContent(Bin) << endl;
  return 1./TriggerWeight;
}

Double_t AliAnalysisTaskHaHFECorrel::GetVtxWeight(Int_t run, Double_t nTrAcc) {
  if (!fUseEventWeights) return 1.;
  Int_t BinX= fVtxWeight.GetXaxis()->FindFixBin(run);
  Int_t BinY= fVtxWeight.GetYaxis()->FindFixBin(nTrAcc);
  if (BinY>=fVtxWeight.GetYaxis()->GetNbins()) BinY = fVtxWeight.GetXaxis()->GetLast();
  Double_t VtxWeight = fVtxWeight.GetBinContent(BinX, BinY);
  if (VtxWeight<0.001) return 1.;
  return 1./VtxWeight;
}

Double_t AliAnalysisTaskHaHFECorrel::GetNonTagCorr(Double_t ptTrack, Double_t ptAsso) {
  Double_t Pt2Max = 5.;
  if (fNonTagCorr.IsZombie()) return -1;
  Int_t Bin = fNonTagCorr.FindBin(ptTrack);
  if (fNonTagCorr.IsBinUnderflow(Bin) || fNonTagCorr.IsBinOverflow(Bin)) {
    return -1.;
  }
  Pt2Max = fNonTagCorr.GetBinContent(Bin);
 
  /*
  if (invMass>0.05) CorrFactor = 1.105;
  else {
    Int_t Bin = fNonTagCorr.FindBin(invMass);
    if (fNonTagCorr.IsBinUnderflow(Bin) || fNonTagCorr.IsBinOverflow(Bin)) {
      // cout <<  "ElecRecEff: " << pt << "\t" << eta << "\t" << zVtx << endl;
      return -1.;
    }
    CorrFactor= fNonTagCorr.GetBinContent(Bin);
  }
  //  cout << "NonTagCorr " << CorrFactor << endl;
  if (CorrFactor>0.000001) return CorrFactor;
  else return -1.;
  */

  if (fPhotCorrCase==2) {
    if (ptAsso<Pt2Max) return 1.;
    else return 0;
  }
  else if (fPhotCorrCase==3) {
    if (ptAsso<Pt2Max*1.1) return 1.;
    else return 0;
  }
  else if (fPhotCorrCase==4) {
    if (ptAsso<Pt2Max*0.9) return 1.;
    else return 0;
  }
  else return 1.;
	 
	 /*
  // Case 0 MCtoMc Case !=0 MctoMcRecCorr
  if (ptAsso<0.2) CorrFactor = 2.95+1.65*ptTrack; 
  else if (ptAsso>0.2 && ptAsso<0.35) {
    if (fPhotCorrCase==0) CorrFactor =  0.055+0.133*ptTrack;
    else CorrFactor=-0.277-0.042*ptTrack;
  }
  else if (ptAsso>0.35 && ptAsso<2.0) {
    if (fPhotCorrCase==0) CorrFactor= 0.052+0.087*ptTrack; 
    else CorrFactor=0.041+0.078*ptTrack;
  }
  else if (ptAsso>2.0) {
    if (fPhotCorrCase==0) CorrFactor=0.060+0.117*ptTrack;
    else CorrFactor = 0.061+0.122*ptTrack;
  }
  return CorrFactor;
	 */
}


Bool_t AliAnalysisTaskHaHFECorrel::PassEventBias( const AliVVertex *pVtx, Int_t nMother, Int_t *listMother, Double_t EventWeight) {
 
  
  // Leading Particle
  Double_t ptMax=-999;
  Double_t pt =-999;
  
  if (fIsMC) {
    for (Int_t i=0; i<fMC->GetNumberOfTracks(); i++) {
      AliVParticle *MCVParticle  = dynamic_cast<AliVParticle*>(fMC->GetTrack(i));
      AliAODMCParticle *MCAODParticle = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(i));
      AliMCParticle *MCESDParticle = dynamic_cast<AliMCParticle*>(fMC->GetTrack(i));
      if (fIsAOD)  {
	if (!MCAODParticle) continue;
	if (!MCAODParticle->IsPhysicalPrimary()) continue;
      }
      else {
	if (!fStack || !MCESDParticle) continue;
	if (!fStack->IsPhysicalPrimary(MCESDParticle->Label())) continue;
      }

      if(MCVParticle->Eta()>fMaxHadronEta) continue;
      if(MCVParticle->Eta()<fMinHadronEta) continue;
      pt=MCVParticle->Pt();
      if (pt>ptMax) ptMax=pt;
    }
  }
  else {
    Int_t ntracks = -999;
    if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
    if(fUseTender)  ntracks = fTracks_tender->GetEntries();

    // Loop over all tracks to find LP and HFE
    for(Int_t jTracks = 0; jTracks < ntracks; jTracks++){
          
      AliVParticle* VHtrack = 0x0;
      if(!fUseTender) VHtrack  = fVevent->GetTrack(jTracks);
      if(fUseTender) VHtrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jTracks));
        
      if (!VHtrack) {
	printf("ERROR: Could not receive associated track %d\n", jTracks);
	continue;
      }
      AliVTrack *Vtrack = dynamic_cast<AliVTrack*>(VHtrack);
      if(!Vtrack) {
	//	cout << "noVTrack" << endl;
	continue;
      }
    
      pt = Vtrack->Pt();

      // track cuts
      Bool_t passHadTrackCut=kFALSE;
      passHadTrackCut = ChargedHadronTrackCuts(pVtx, Vtrack, nMother, listMother, EventWeight);
   
      Bool_t passHadPIDCut=kFALSE;
      passHadPIDCut = ChargedHadronPIDCuts(Vtrack, EventWeight); // currently empty
    
      // find hadron with the largest pT -> leading particle
      if (passHadTrackCut && passHadPIDCut) {
	if (pt>ptMax) ptMax=pt;
      }
    }
  }
  
    if (ptMax>fMinPtEvent && ptMax<fMaxPtEvent) return kTRUE;
    else return kFALSE;
}

Bool_t AliAnalysisTaskHaHFECorrel::ESDkTrkGlobalNoDCA(AliVTrack* Vtrack) {
  AliESDtrack *ESDtrack = dynamic_cast<AliESDtrack*>(Vtrack);   
  if (!fIsAOD && !ESDtrack && !fesdTrackCuts) return kFALSE;
  fesdTrackCuts->SetMaxDCAToVertexXY(2.4);
  fesdTrackCuts->SetMaxDCAToVertexZ(3.2);
  fesdTrackCuts->SetDCAToVertex2D(kTRUE);
  if (fesdTrackCuts->AcceptTrack(ESDtrack)) return kTRUE;
  else return kFALSE;
}

void AliAnalysisTaskHaHFECorrel::MCTruthCorrelation(TObjArray* MCTrueRedTracks, Bool_t AfterEventCuts, Int_t RecLPLabel,  Float_t pVtxZ, Float_t mult,  Int_t &LeadingParticleInAcceptance, Int_t &LeadingParticle, Double_t EventWeight) {
  Int_t PDGCode=-999, Mother=-999, MotherPDG=-999, MotherIsHeavy=-999; // GrandMother=-999, GrandMotherPDG=-999, GrandMotherIsHeavy=-999;
  Bool_t ElectronInAcceptanceCut=kFALSE;
  Bool_t HadronInAcceptanceCut=kFALSE;
  LeadingParticleInAcceptance=-999, LeadingParticle=-999;
  Double_t  LeadingParticlePtInAcceptance=-99, LeadingParticlePt=-99;

  Bool_t **ElectronIsTrigger = new Bool_t*[11]; // Check if electron is trigger for all cases
  for (Int_t i=0; i<11; i++) { // case
    ElectronIsTrigger[i]= new Bool_t[fAssPtHad_Nbins]; // associated hadron bin
    for (Int_t j=0; j<fAssPtHad_Nbins; j++) ElectronIsTrigger[i][j]=kFALSE;
  }

 
 
  if (fIsAOD) {
    // Find HFE Eelcton
    for(Int_t iMCElectron = 1; iMCElectron < (fMC->GetNumberOfTracks()); iMCElectron++) {
      AliAODMCParticle* MCElectron = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(iMCElectron));  
      if (!MCElectron) continue;
      if (!MCElectron->IsPhysicalPrimary() || MCElectron->Charge()==0) continue; // strong and electronweak decays safe for particles not reaching the detector
      PDGCode = abs(MCElectron->GetPdgCode());
   
      if (PDGCode==11) { // should be an electron
	for (Int_t i=0; i<11; i++) { 
	  for (Int_t j=0; j<fAssPtHad_Nbins; j++) ElectronIsTrigger[i][j]=kFALSE; // reset trigger 
	}
	Mother = MCElectron->GetMother();
	if (Mother>=0) { // Mother exists
	  AliAODMCParticle* MCElectronMother= dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(Mother));
	  MotherPDG = abs(MCElectronMother->GetPdgCode()); 
	  if (MotherPDG ==11 || MotherPDG==15) { // e from e e from tau
	    Mother = MCElectronMother->GetMother();
	    if (Mother>=0) {
	      AliAODMCParticle* MCElectronGrandMother= dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(Mother));
	      MotherPDG = abs(MCElectronGrandMother->GetPdgCode()); 
	      // GrandMotherIsHeavy  = Int_t (GrandMotherPDG / TMath::Power(10, Int_t(TMath::Log10(GrandMotherPDG))));
	      // if (GrandMotherIsHeavy>3) cout << "GrandMotherIsHeavy: " << GrandMotherIsHeavy << endl;
	    }
	  }
	  MotherIsHeavy  = Int_t (MotherPDG / TMath::Power(10, Int_t(TMath::Log10(MotherPDG)))); // only mother of first particle with rough selection
	 
	  if (MotherIsHeavy>3 && MotherIsHeavy<6) {
	    MotherIsHeavy = HFEisCharmOrBeauty(iMCElectron); // looking in history for B-D return 4,5 or 0
	  }  

	  ElectronInAcceptanceCut = kFALSE;
	  if ((MCElectron->Eta() < fMaxElectronEta) && (MCElectron->Eta() > fMinElectronEta)) ElectronInAcceptanceCut = kTRUE;

	  if (MotherIsHeavy>3 && MotherIsHeavy<6) { // start Hadron loop
	    if (AfterEventCuts) {
	      // fTrueElectronEta->Fill(1., 1., 1., 0.5);
	      if (!MCElectron->IsPhysicalPrimary()) cout << "TruEle not primary" << endl;
	      fTrueElectronEta->Fill(MCElectron->Pt(), MCElectron->Eta(), 1.*mult,  EventWeight);
	    }
	    
	    if (ElectronInAcceptanceCut) {
	      AliBasicParticleHaHFE * ElectronParticle  = 0;
	      AliExternalTrackParam ExtTrackParam;
	      ElectronParticle = new AliBasicParticleHaHFE(MotherIsHeavy, MCElectron->Eta(), MCElectron->Phi(), MCElectron->Pt(), MCElectron->Charge(), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, MCElectron->Label(),ExtTrackParam);  // Original mixed event particle class: id has been filled with c or b flag
	      if (ElectronParticle) MCTrueRedTracks->Add(ElectronParticle);
	    }

	    Int_t CharmOrBeauty = MotherIsHeavy*5-20;
	    LeadingParticleInAcceptance = -999;
	    LeadingParticlePtInAcceptance = -999;
	    LeadingParticle = -999;
	    LeadingParticlePt = -999;
	    if (fCorrHadron || fCorrLParticle) {
	      for (Int_t iMCHadron = 1; iMCHadron< fMC->GetNumberOfTracks(); iMCHadron++) {
		AliAODMCParticle* MCHadron = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(iMCHadron));  
		if (!MCHadron) continue;
		if ((!MCHadron->IsPhysicalPrimary()) || ( MCHadron->Charge()==0)) continue; 
		if (AfterEventCuts && fOneTimeCheck) fTrueHadronEta->Fill(MCHadron->Pt(), MCHadron->Eta(), EventWeight);
		
		HadronInAcceptanceCut = kFALSE;
		if ((MCHadron->Eta() < fMaxHadronEta) && (MCHadron->Eta() > fMinHadronEta)) HadronInAcceptanceCut = kTRUE;
		
		if (HadronInAcceptanceCut  && (MCHadron->Pt() > LeadingParticlePtInAcceptance)) {
		  LeadingParticleInAcceptance = abs(MCHadron->Label());
		  LeadingParticlePtInAcceptance = MCHadron->Pt();
		}
		if (MCHadron->Pt() > LeadingParticlePt) {
		  LeadingParticle = abs(MCHadron->Label());
		  LeadingParticlePt = MCHadron->Pt();
		}
		
		if (MCElectron->Label() == MCHadron->Label()) continue; // self correlation
	      
		Double_t FillSparse[5];
		FillSparse[0]=MCHadron->Pt();
		FillSparse[1]=MCElectron->Pt();
		FillSparse[2]=GetDeltaPhi(MCElectron->Phi(), MCHadron->Phi());
		FillSparse[3]=GetDeltaEta(MCElectron->Eta(), MCHadron->Eta());
	      
		if (ElectronInAcceptanceCut && HadronInAcceptanceCut) {
		  FillSparse[4]=CharmOrBeauty+1;
		  CheckElectronIsTrigger(MCHadron->Pt(), ElectronIsTrigger[CharmOrBeauty+1]);
		  if (fOneTimeCheck) {
		    if (AfterEventCuts) fTrueMCHadronEventCuts->Fill(FillSparse, EventWeight);
		    else fTrueMCHadron->Fill(FillSparse, EventWeight);
		  }
		  if (AfterEventCuts) {
		    FillSparse[4]=pVtxZ;
		    fTrueMCHadronEventCutsZvtx->Fill(FillSparse, EventWeight);
		  }
		}
		if (ElectronInAcceptanceCut) {
		  FillSparse[4]=CharmOrBeauty+2;
		  CheckElectronIsTrigger(MCHadron->Pt(), ElectronIsTrigger[CharmOrBeauty+2]);
		  if (fOneTimeCheck) {
		    if (AfterEventCuts) fTrueMCHadronEventCuts->Fill(FillSparse, EventWeight);
		    else fTrueMCHadron->Fill(FillSparse, EventWeight);
		  }
		}
		if (HadronInAcceptanceCut) {
		  FillSparse[4]=CharmOrBeauty+3;
		  CheckElectronIsTrigger(MCHadron->Pt(), ElectronIsTrigger[CharmOrBeauty+3]);
		  if (fOneTimeCheck) {
		    if (AfterEventCuts) fTrueMCHadronEventCuts->Fill(FillSparse, EventWeight);
		    else fTrueMCHadron->Fill(FillSparse, EventWeight);
		  }
		}
	      }
	   
	    
	      // fill trigger histogram for elec - hadron
	      for (Int_t i=0; i<11; i++) { // case
		for (Int_t j=0; j<fAssPtHad_Nbins; j++) {
		  if (ElectronIsTrigger[i][j]) {
		    if (AfterEventCuts) {
		      fTrueMCElecHaTriggerEventCuts->Fill(MCElectron->Pt(), i, j, EventWeight);
		    }
		    else if (fOneTimeCheck) {
		      fTrueMCElecHaTrigger->Fill(MCElectron->Pt(), i, j, EventWeight);
		    }
		  }
		  ElectronIsTrigger[i][j]=kFALSE;
		}
	      }
	      
	  
	      // Fill LeadingPartilce
	      if (LeadingParticleInAcceptance!=-999 && LeadingParticle !=-999) {
		
		Double_t FillSparseLP[5];
		
		if (MCElectron->Label() != LeadingParticleInAcceptance) {
		AliAODMCParticle* MCLPinAcceptance = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(LeadingParticleInAcceptance));  
		if (AfterEventCuts && fOneTimeCheck) fTrueLPinAcceptanceEta->Fill(MCLPinAcceptance->Pt(), MCLPinAcceptance->Eta(), EventWeight);

	
		FillSparseLP[0]=MCLPinAcceptance->Pt();
		FillSparseLP[1]=MCElectron->Pt();
		FillSparseLP[2]=GetDeltaPhi(MCElectron->Phi(), MCLPinAcceptance->Phi());
		FillSparseLP[3]=GetDeltaEta(MCElectron->Eta(), MCLPinAcceptance->Eta());
		
		if (ElectronInAcceptanceCut) { //electron in acceptancen and LP within particles in acceptance
		  FillSparseLP[4]=CharmOrBeauty+1;
		  CheckElectronIsTrigger(MCLPinAcceptance->Pt(), ElectronIsTrigger[CharmOrBeauty+1]);
		  if (fOneTimeCheck) {		 
		    if (AfterEventCuts) fTrueMCLPEventCuts->Fill(FillSparseLP, EventWeight);
		    else fTrueMCLP->Fill(FillSparseLP, EventWeight);
		  }
		}
		FillSparseLP[4]=CharmOrBeauty+3; // all eclectron with LP from acceptance
		CheckElectronIsTrigger(MCLPinAcceptance->Pt(), ElectronIsTrigger[CharmOrBeauty+3]);
		if (fOneTimeCheck) {
		  if (AfterEventCuts) fTrueMCLPEventCuts->Fill(FillSparseLP, EventWeight);
		  else fTrueMCLP->Fill(FillSparseLP, EventWeight);
		}
		}	 

		if (MCElectron->Label() != LeadingParticle) {
		  AliAODMCParticle* MCLP = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(LeadingParticle));  
		  if(AfterEventCuts && fOneTimeCheck) fTrueLPEta->Fill(MCLP->Pt(), MCLP->Eta(), EventWeight);
		  
		  FillSparseLP[0]=MCLP->Pt();
		  FillSparseLP[2]=GetDeltaPhi(MCElectron->Phi(), MCLP->Phi());
		  FillSparseLP[3]=GetDeltaEta(MCElectron->Eta(), MCLP->Eta());
		  if (ElectronInAcceptanceCut) { // electron in acceptance with true LP
		    FillSparseLP[4]=CharmOrBeauty+2;
		    CheckElectronIsTrigger(MCLP->Pt(), ElectronIsTrigger[CharmOrBeauty+2]);
		    if (fOneTimeCheck) {
		      if (AfterEventCuts) fTrueMCLPEventCuts->Fill(FillSparseLP, EventWeight);
		      else fTrueMCLP->Fill(FillSparseLP,EventWeight);
		    }
		  }
		  
		  FillSparseLP[4]=CharmOrBeauty+4;
		  CheckElectronIsTrigger(MCLP->Pt(), ElectronIsTrigger[CharmOrBeauty+4]);
		  if (fOneTimeCheck) { // neighter in acceptance
		    if (AfterEventCuts) fTrueMCLPEventCuts->Fill(FillSparseLP,EventWeight);
		    else fTrueMCLP->Fill(FillSparseLP, EventWeight);
		  }
		}
		
		if (RecLPLabel>0 && ( abs(MCElectron->Label()) != abs(RecLPLabel))) {
		  AliAODMCParticle* MCLPRec = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(RecLPLabel));  
		  if(AfterEventCuts) fRecLPEta->Fill(MCLPRec->Pt(), MCLPRec->Eta(), EventWeight);
		  
		  FillSparseLP[0]=MCLPRec->Pt();
		  FillSparseLP[2]=GetDeltaPhi(MCElectron->Phi(), MCLPRec->Phi());
		  FillSparseLP[3]=GetDeltaEta(MCElectron->Eta(), MCLPRec->Eta());
		  if (ElectronInAcceptanceCut) { // electron in acceptance with rec LP
		    FillSparseLP[4]=CharmOrBeauty+5;
		    CheckElectronIsTrigger(MCLPRec->Pt(), ElectronIsTrigger[CharmOrBeauty+5]);
		    if (fOneTimeCheck) {
		      if (AfterEventCuts) fTrueMCLPEventCuts->Fill(FillSparseLP, EventWeight);
		      else fTrueMCLP->Fill(FillSparseLP, EventWeight);
		    }
		    
		    if (AfterEventCuts ) {
		      FillSparseLP[4]=pVtxZ;
		      fTrueMCLPEventCutsZvtx->Fill(FillSparseLP, EventWeight);
		    }
		  }
		}
	      }
	      
	      
	      // fill trigger histogram for electron leading particle
	      for (Int_t i=0; i<11; i++) { // case
		for (Int_t j=0; j<fAssPtHad_Nbins; j++) {
		  if (ElectronIsTrigger[i][j]) {
		    if (AfterEventCuts) {
		      fTrueMCElecLPTriggerEventCuts->Fill(MCElectron->Pt(), i, j, EventWeight);;
		    }
		    else if (fOneTimeCheck) {
		      fTrueMCElecLPTrigger->Fill(MCElectron->Pt(), i, j, EventWeight);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  for (Int_t i=0; i<11; i++) { // case
    delete [] ElectronIsTrigger[i];
  }
  delete [] ElectronIsTrigger;



  if (fCorrHadron || fCorrLParticle) {
    if (AfterEventCuts && MCTrueRedTracks->GetEntriesFast()>0 && RecLPLabel>0) {
      AliAODMCParticle* MCLPRec = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(abs(RecLPLabel)));  
      
      AliEventPool * HFEPool = fMCTruePoolMgr->GetEventPool(mult, pVtxZ, MCLPRec->Pt()); // Get the buffer associated with the current centrality and z-vtx
      
      if (!HFEPool) {
	AliFatal(Form("No pool found for centrality = %f, zVtx = %f, leading pt = %f", mult, pVtxZ, MCLPRec->Pt()));
	return;
      }
      
      if (HFEPool->GetCurrentNEvents() >= 3) { // start mixing when 3 events are in the buffer
	Int_t nMix = HFEPool->GetCurrentNEvents();   
	
	
	for (Int_t iMCHadron = 1; iMCHadron< fMC->GetNumberOfTracks(); iMCHadron++) {
	  AliAODMCParticle* MCHadron = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(iMCHadron));  
	  if (!MCHadron) continue;
	  if ((!MCHadron->IsPhysicalPrimary()) || ( MCHadron->Charge()==0)) continue; 
	  if ((MCHadron->Eta() > fMaxHadronEta) || (MCHadron->Eta() < fMinHadronEta)) continue;
	  
	  
	  for (Int_t jMix=0; jMix<nMix; jMix++) { // mix with each event in the buffer
	    TObjArray* mixTracks = HFEPool->GetEvent(jMix);
	    for (Int_t i=0;i<mixTracks->GetEntriesFast(); i++)  {	
	      AliBasicParticleHaHFE* mixtrk = (AliBasicParticleHaHFE*) mixTracks->At(i);
	      if (!mixtrk) {
		printf("ERROR: Could not receive mix pool track %d\n",i);
		continue;
	      }
	      
	      Double_t FillSparse[5];
	      FillSparse[0]=MCHadron->Pt();
	      FillSparse[1]=mixtrk->Pt();
	      FillSparse[2]=GetDeltaPhi(mixtrk->Phi(), MCHadron->Phi());
	      FillSparse[3]=GetDeltaEta(mixtrk->Eta(), MCHadron->Eta());
	      FillSparse[4]=pVtxZ;
	      fTrueMCHadronEventCutsZvtxMEv->Fill(FillSparse, 1.); // no EventWeight in MEv
	      if (abs(MCHadron->Label()) == abs(RecLPLabel))  fTrueMCLPEventCutsZvtxMEv->Fill(FillSparse, 1.); // no EventWeight in MEv
	    }
	  }
	}   
      }     
      HFEPool->UpdatePool(MCTrueRedTracks);
    }
  }
  //else delete MCTrueRedTracks;
  // cout << "e" << endl;
}




Int_t AliAnalysisTaskHaHFECorrel::HFEisCharmOrBeauty(Int_t ElectronIndex) {
  AliAODMCParticle* HFParticle = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(ElectronIndex));  
  if (!HFParticle) return 0;
  if (!HFParticle->IsPhysicalPrimary() || HFParticle->Charge()==0) return 0; // strong and e
  Int_t HFDaughter = ElectronIndex;
  Int_t HFMother = HFParticle->GetMother();
  Int_t PDGCode = abs(HFParticle->GetPdgCode());
  //  cout<< "Initial PDG " << PDGCode << endl;
  while ( PDGCode>100 || PDGCode ==11) {
    HFDaughter = HFMother;
    HFMother = HFParticle->GetMother();
    if (HFMother>=0) {
      HFParticle = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(HFMother));  
      if (!HFParticle) break;   
      PDGCode = abs(HFParticle->GetPdgCode());
      //   cout << PDGCode << endl;
    }
    else break;
  }
  if ( PDGCode==4) return 4; // charm quark
  else if (PDGCode==5) return 5; // bottom quark
  else if (PDGCode<100) { // neither charm or bottom - look at first Meson
    AliAODMCParticle* HFParticle = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(HFDaughter));  
    if (!HFParticle) return 0;
    Int_t FirstMesonPDG = abs(HFParticle->GetPdgCode());
    return Int_t (FirstMesonPDG / TMath::Power(10, Int_t(TMath::Log10(FirstMesonPDG))));
  }
  else if (PDGCode>100) return 0;

  return 0;
}

Bool_t AliAnalysisTaskHaHFECorrel::Thrust(const TObjArray* tracks, Double_t t[2], Double_t MaxEta, Double_t MinPt)
{
  // compute thrust value[0] and direction[1]

  memset(t, 0, 2*sizeof(Double_t));
  Int_t ntracks(0);
  if(!(ntracks=tracks->GetEntries())) return kFALSE;

  Double_t deltaPhi=0.05*TMath::Pi()/180., //the resolution of the thrust orientation
           absPtSum(0.);

  AliVParticle *track(NULL);
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    if(!(track = dynamic_cast<AliVParticle*>(tracks->At(iTracks)))) continue;
    if (TMath::Abs(track->Eta())>MaxEta) continue;
    if (track->Pt()<MinPt) continue;
    absPtSum += track->Pt();
  }

  for (Double_t nphi(0.); nphi<=TMath::TwoPi(); nphi+=deltaPhi){//nphi loop
    Double_t thru(0.);
    for(Int_t j=0;j<ntracks;j++){
      if(!(track = dynamic_cast<AliVParticle*> (tracks->At(j)))) continue;
    if (TMath::Abs(track->Eta())>MaxEta) continue;
    if (track->Pt()<MinPt) continue;
      thru+=track->Pt()*TMath::Abs(TMath::Cos(track->Phi()-nphi));
    }
    if(thru>t[0]){
      t[0] = thru;
      t[1] = nphi;
    }
  }//nphi loop

  if(absPtSum<kAlmost0) return kFALSE;
  t[0] /= absPtSum;
  return kTRUE;
}

Double_t AliAnalysisTaskHaHFECorrel::Sphericity(const TObjArray* tracks, Double_t MaxEta, Double_t MinPt)
{
    // compute sphericity
  
  Int_t ntracks(0);
  if(!(ntracks=tracks->GetEntries())) return -1.;

  Double_t a(0.), b(0.), c(0.), d(0.);
  AliAODMCParticle* track;
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    if(!(track =  dynamic_cast<AliAODMCParticle*> (tracks->At(iTracks)))) continue;
    if(TMath::Abs(track->Eta()) > MaxEta) continue;
    if(track->Pt()<MinPt) continue;
    a+=track->Px()*track->Px();
    b+=track->Px()*track->Py();
    d+=track->Py()*track->Py();
  }
  c = b;
  Double_t ad=a+d;
  Double_t delta=TMath::Sqrt((ad)*(ad)-4.*(a*d-b*c));
  Double_t lambda1=(ad-delta)/2.;
  //lambda2=(ad+delta)/2.;
  return ad>kAlmost0?(2.*lambda1/ad):-1;
}


void AliAnalysisTaskHaHFECorrel::PhotULSLSElectronAcceptance(const AliVVertex *pVtx, Float_t mult,  Int_t nMother, Int_t listMother[], Double_t EventWeight) {
  
  AliEventPool* fPool;
  fPool = fPoolMgr->GetEventPool(mult, pVtx->GetZ(), 1.0); // Get the buffer associated with the current centrality and z-vtx
  //  fPool->SetDebug(kTRUE);
  if (!fPool)
  {
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f, maxPt = %f", mult, pVtx->GetZ(), 1.));
    return;
  }

  if (fPool->GetCurrentNEvents() >= 3){// start mixing when 3 events are in the buffer
    Int_t nMix = fPool->GetCurrentNEvents();   
    for (Int_t jMix=0; jMix<nMix; jMix++){  // mix with each event in the buffer
      TObjArray* mixTracks = fPool->GetEvent(jMix);

      for (Int_t i=0;i<mixTracks->GetEntriesFast(); i++) {	
	AliBasicParticleHaHFE* mixtrk = (AliBasicParticleHaHFE*) mixTracks->At(i);
	if (!mixtrk) {
	  printf("ERROR: Could not receive mix pool track %d\n",i);
	  continue;
	}
	  
	Double_t ptMix=-9., phiMix=-9., etaMix=-9, recEff=-9;
	Int_t ls=-9, uls=-9;
	Short_t chargeMix=0;;
	ptMix=mixtrk->Pt();
	phiMix=mixtrk->Phi();
	etaMix=mixtrk->Eta();
	chargeMix = mixtrk->Charge();
	
	Int_t ntracks = -999;
	if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
	if(fUseTender)  ntracks = fTracks_tender->GetEntries();
	
	// Loop over all tracks to find photonic partner
	for(Int_t jTracks = 0; jTracks < ntracks; jTracks++){
	  
	 
	  AliVParticle* Vassotrack = 0x0;
	  if(!fUseTender) Vassotrack  = fVevent->GetTrack(jTracks);
	  if(fUseTender)  Vassotrack  = dynamic_cast<AliVTrack*>(fTracks_tender->At(jTracks));
	  if (!Vassotrack) {
	    printf("ERROR: Could not receive associated track %d\n", jTracks);
	    continue;
	  }
	  AliVTrack *VtrackAsso = dynamic_cast<AliVTrack*>(Vassotrack);
	  if(!VtrackAsso) continue;
	  AliAODTrack *AODtrackAsso = dynamic_cast<AliAODTrack*>(VtrackAsso);   
	  if(!AODtrackAsso) continue;
	  
	  // associated particle variables
	  Double_t phiAsso=-9.;
	  Int_t chargeAsso = 0;
	  phiAsso = VtrackAsso->Phi();
	  chargeAsso = VtrackAsso->Charge();

	  // looser Track cuts for associated photonic eletron
	  Bool_t passAssoTrackCutPhot = kFALSE;
	  passAssoTrackCutPhot = PhotElecTrackCuts(pVtx,VtrackAsso,nMother,listMother, EventWeight);
	  if(!passAssoTrackCutPhot) continue;
   
	  Bool_t passAssoPIDCutPhot = kFALSE;
	  passAssoPIDCutPhot = PhotElecPIDCuts(VtrackAsso, EventWeight); 
	  if (!passAssoPIDCutPhot) continue;         
   
	  Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
	  Double_t openingAngle = -999., mass=999., width = -999;
	  Double_t energy=-999., transenergy=-999., mom=-999., tansmom=-999.;
        
	  Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
	  if(chargeMix>0) fPDGe1 = -11;
	  if(chargeAsso>0) fPDGe2 = -11;
        
	  if(chargeMix == chargeAsso) fFlagLS = kTRUE;
	  if(chargeMix != chargeAsso) fFlagULS = kTRUE;
    
	  //Variables
	  Double_t p1[3], p2[3];
	  Double_t xt1, xt2; //radial position track 1 and 2  at the DCA point
	  Double_t dca12; //DCA 1-2
	  Bool_t hasdcaT1, hasdcaT2;
	  Double_t bfield = fVevent->GetMagneticField();
	  
	  AliExternalTrackParam extTrackParam1;
	  extTrackParam1=mixtrk->GetExtTrackParam();
	  AliExternalTrackParam extTrackParam2;
	  extTrackParam2.CopyFromVTrack(VtrackAsso);

	  //DCA track1-track2
	  dca12 = extTrackParam2.GetDCA(&extTrackParam1,bfield,xt2,xt1);
	  
	  //Momento of the track extrapolated to DCA track-track
	  //Track1
	  hasdcaT1 = extTrackParam1.GetPxPyPzAt(xt1,bfield,p1);
	  //Track2
	  hasdcaT2 = extTrackParam2.GetPxPyPzAt(xt2,bfield,p2);
	  
	  if(!hasdcaT1 || !hasdcaT2) AliWarning("There could be a problem in the extrapolation");
	  //track1-track2 Invariant Mass
	  Double_t eMass = 0.000510998910; //Electron mass in GeV
	  Double_t pP1 = sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]); //Track 1 momentum
	  Double_t pP2 = sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]); //Track 1 momentum

	  TLorentzVector v1(p1[0],p1[1],p1[2],sqrt(eMass*eMass+pP1*pP1));
	  TLorentzVector v2(p2[0],p2[1],p2[2],sqrt(eMass*eMass+pP2*pP2));
	  mass = (v1+v2).M(); //Invariant Mass
	  openingAngle = v1.Angle(v2.Vect()); //Opening Angle (Total Angle)
    
	  if(fFlagLS){
	    //   fInvmassLS->Fill(mass,ptMix);
	  }
	  if(fFlagULS){
	    //  fInvmassULS->Fill(mass,ptMix);
	  }
      
	  if(mass<fInvmassCut && fFlagULS) {
	    if (fOneTimeCheck) fPhotMixULS->Fill(ptMix); // EventWeight are not import, events are mixed within the same mult class
	  }
        
	  if(mass<fInvmassCut && fFlagLS){
	    if (fOneTimeCheck) fPhotMixLS->Fill(ptMix, chargeMix); // EventWeight are not important
	  }
  
	}//track loop
      } // mixtracks
    }// mix events
  }	  
}


void AliAnalysisTaskHaHFECorrel::SetEleVarOpt(Int_t VarOpt) {
  fVarEleOpt=VarOpt;
  switch (fVarEleOpt) {
  case 0:
    return;
    break;
  case 1:
    fMinElectronEta=-0.7; // defaul 0.8
    fMaxElectronEta=0.7;
    break;
  case 2:
    fElectronkAny=kTRUE;
    break;
  case 3:
    fITSnCut=4; // default 3
    break;
  case 4:
    fITSnCut=5;
    break;
  case 5:
    fEleDCAr=0.75; // default 1. and 2.
    fEleDCAz=1.5;
    break;
  case 6:
    fEleDCAr=1.25;
    fEleDCAz=2.5;
    break;
  case 7:
    fTPCnCut=90; // default 100
    break;
  case 8:
    fTPCnCut=110;
    break;
  case 9:
    fTPCndEdxCut=70;
    break;
  case 10:
    fTPCndEdxCut=90;
    break;
  case 11:
    fSigmaTPCcutLow=0; // default -1;
    break;
  case 12:
    fSigmaTPCcutHigh=2;
    break;
  case 13:
    fSigmaTOFcut=2;
    break;
  case 14:
    fSigmaITScut=2;
    break;
  case 15:
    fSigmaTPCcutLow=-1.5; // default -1;
    break;
  case 16:
    fSigmaTOFcut=2.5;
    break;
  case 17:
    fMinElectronEta=-0.6; // default 0.8
    fMaxElectronEta=0.6;
    break;
  case 18:
    fMinElectronEta=0; // default 0.8
    fMaxElectronEta=0.8;
    break;
  case 19:
    fMinElectronEta=-0.8; // default 0.8
    fMaxElectronEta=0;
    break;
  case 20:
    fElectronkFirst=kTRUE;
    break;
  }


}



void AliAnalysisTaskHaHFECorrel::SetHadVarOpt(Int_t VarOpt) {
  fVarHadOpt=VarOpt;
  switch (fVarHadOpt) {
  case 0:
    return;
    break;
  case 1:
    fMinHadronEta=-0.7;
    fMaxHadronEta=0.7;
    break;
  case 2:
    fMinHadronEta=-0.6;
    fMaxHadronEta=0.6;
    break;
  case 3:
    fHTPCnCut=80;
    break;
  case 4:
    fHTPCnCut=120;
    break;
  case 5:
    fHadDCAr=0.5;
    fHadDCAz=1.5;
    break;
  case 6:
    fHadDCAr=0.75;
    fHadDCAz=2;
    break;
  case 7:
    fHadDCAr=0.5;
    break;
  case 8:
    fHadkAny=kFALSE;
    break;
  case 9:
    fHadTOFmatch=kTRUE;
    break;
  }
      
}

void AliAnalysisTaskHaHFECorrel::SetPhotVarOpt(Int_t VarOpt) {
  fVarPhotOpt=VarOpt;
  switch (fVarPhotOpt) {
  case 0:
    return;
    break;
  case 1:
    fUseITSsa = kFALSE;
    break;
  case 2:
    fPhotElecTPCnCut = 60;
    break;
  case 3:
    fPhotElecTPCnCut = 100;
    break;
  case 4:
    fPhotElecSigmaTPCcut = 2;
    break;
  case 5:
    fInvmassCut=0.08;
    break;
  case 6:
    fInvmassCut=0.10;
    break;
  case 7:
    fInvmassCut=0.12;
    break;
  case 8:
    fInvmassCut=0.16;
    break;
  case 9:
    fPhotElecPtCut=0.1;
    break;
  case 10:
    fPhotElecPtCut=0.2;
    break;
  case 11:
    fPhotElecPtCut=0.3;
    break;
  case 12:
    fPhotCorrCase = 3;
    break;
  case 13:
    fPhotCorrCase = 4;
    break;
  case 14: // Tilted up weights
    return;
    break;
  case 15: // Tilted down weights
    return;
    break;
  }
}

Int_t AliAnalysisTaskHaHFECorrel::CheckParticleOrigin(Int_t Label) {
  AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(Label));    
  Int_t PDGCode = abs(MCParticle->GetPdgCode());
  cout << "Track with label: " << Label << "\t is a " << PDGCode << endl;
  while( MCParticle->GetMother()>=0) {
    MCParticle = (AliAODMCParticle*) fMC->GetTrack(MCParticle->GetMother());
    PDGCode = abs(MCParticle->GetPdgCode());
    cout << PDGCode << "\t";
  }
  cout << endl;
  return PDGCode;

}
