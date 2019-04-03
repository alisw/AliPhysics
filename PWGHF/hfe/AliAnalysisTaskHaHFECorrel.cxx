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
,fUseTender(kFALSE)
,fWhichPeriod(2016)
,fUseKFforPhotonicPartner(kFALSE)
,fMaxPtEvent(999)
,fMinPtEvent(0)
,fMaxElectronEta(0.8)
,fMinElectronEta(-0.8)
,fMaxHadronEta(0.9)
,fMinHadronEta(-0.9)
,fTPCnCut(100)
,fTPCndEdxCut(80)
,fITSnCut(3)
,fITSSharedClusterCut(999)
,fUseTRD(0)
,fUseITS(0)
,fSigmaITScut(2.)
,fSigmaTOFcut(2.)
,fSigmaTPCcut(-1.)
,fPhotElecPtCut(0.125)
,fPhotElecSigmaTPCcut(3)
,fPhotElecTPCnCut(80)
,fPhotElecITSrefitCut(kTRUE)
,fPhotCorrCase(2)
,fAssNonEleTPCcut(-4)
,fHTPCnCut(100)
,fHITSrefitCut(kTRUE)
,fHTPCrefitCut(kTRUE)
,fOpeningAngleCut(1000.)
,fInvmassCut(0.14)
,fChi2Cut(3.5)
,fDCAcut(999)
,fTRDQA(kFALSE)
,fMCTrueCorrelation(kTRUE)
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
,fTrkpt(0)
,fEtaVtxZ(0)
,fSPDVtxRes(0)              
,fDiffSPDPrimVtx(0)         
,fSPDnTrAcc(0)    
,fSPDnTrCorrMax(0)
,fSPDnTrGen(0)              
,fDiffSPDMCVtx(0)           
,fnTrAccMaxGen(0)           
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
,fElectronTrackTPCNcls(0)
,fElectronTrackTPCNclsdEdx(0)
,fElectronTrackTPCFrac(0)
,fElectronTrackITSNcls(0)
,fElectronTrackITSLayer(0)
,fElectronTrackRefit(0)
,fElectronTrackDCA(0)
,fElectronTrackITSCuts(0)
,fPhotTrackITSCuts(0)
,fHadronTrackTPCNcls(0)
,fHadronTrackRefit(0)
,fHadronTrackDCA(0)
,fHistITSnSig(0)
,fHistTOFnSig(0)
,fHistTPCnSig(0)
,fHistTPCnSigITScut(0)
,fHistTPCnSigTOFcut(0)
,fHistTPCnSigITSTOFcut(0)
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
,fHadRecEff()
,fEleRecEff()
,fSPDnTrAvg()
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
,fRecHadMCPtEtaPhiVtx(0)
,fRecHadPtEtaPhiVtx(0)
,fRecHadPtEtaPhiVtxWRecEff(0)
,fCheckMCPtvsRecPtEle(0)
,fRecHFE(0)
,fMCElecPtEtaPhiVtx(0)
,fRecElecPtEtaPhiVtx(0)
,fRecElecPtEtaPhiVtxWRecEff(0)
,fRecElecMCPtEtaPhiVtx(0)
,fMCElecPDG(0)
,fMCElecPtEtaPhiStrictVtx(0)
,fMCPi0Prod(0)
,fMCEtaProd(0)
,fMCPiPlusProd(0)
,fMCPiPlusProdV2(0)
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
,fUseTender(kFALSE)
,fWhichPeriod(2016)
,fUseKFforPhotonicPartner(kFALSE)
,fMaxPtEvent(999)
,fMinPtEvent(0)
,fMaxElectronEta(0.8)
,fMinElectronEta(-0.8)
,fMaxHadronEta(0.9)
,fMinHadronEta(-0.9)
,fTPCnCut(100)
,fTPCndEdxCut(80)
,fITSnCut(3)
,fITSSharedClusterCut(999)
,fUseTRD(0)
,fUseITS(0)
,fSigmaITScut(2.)
,fSigmaTOFcut(2.)
,fSigmaTPCcut(-1.)
,fPhotElecPtCut(0.125)
,fPhotElecSigmaTPCcut(3)
,fPhotElecTPCnCut(80)
,fPhotElecITSrefitCut(kTRUE)
,fPhotCorrCase(1)
,fAssNonEleTPCcut(-4)
,fHTPCnCut(100)
,fHITSrefitCut(kTRUE)
,fHTPCrefitCut(kTRUE)
,fOpeningAngleCut(1000.)
,fInvmassCut(0.14)
,fChi2Cut(3.5)
,fDCAcut(999)
,fTRDQA(kFALSE)
,fMCTrueCorrelation(kTRUE)
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
,fTrkpt(0)
,fEtaVtxZ(0)
,fSPDVtxRes(0)              
,fDiffSPDPrimVtx(0)         
,fSPDnTrAcc(0)  
,fSPDnTrCorrMax(0)            
,fSPDnTrGen(0)              
,fDiffSPDMCVtx(0)           
,fnTrAccMaxGen(0)           
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
,fElectronTrackTPCNcls(0)
,fElectronTrackTPCNclsdEdx(0)
,fElectronTrackTPCFrac(0)
,fElectronTrackITSNcls(0)
,fElectronTrackITSLayer(0)
,fElectronTrackRefit(0)
,fElectronTrackDCA(0)
,fElectronTrackITSCuts(0)
,fPhotTrackITSCuts(0)
,fHadronTrackTPCNcls(0)
,fHadronTrackRefit(0)
,fHadronTrackDCA(0)
,fHistITSnSig(0)
,fHistTOFnSig(0)
,fHistTPCnSig(0)
,fHistTPCnSigITScut(0)
,fHistTPCnSigTOFcut(0)
,fHistTPCnSigITSTOFcut(0)
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
,fHadRecEff()
,fEleRecEff()
,fSPDnTrAvg()
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
,fTPartnerNoTPt2(0),fElecHadTrigger(0)
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
,fRecHadMCPtEtaPhiVtx(0)
,fRecHadPtEtaPhiVtx(0)
,fRecHadPtEtaPhiVtxWRecEff(0)
,fCheckMCPtvsRecPtEle(0)
,fRecHFE(0)
,fMCElecPtEtaPhiVtx(0)
,fRecElecPtEtaPhiVtx(0)
,fRecElecPtEtaPhiVtxWRecEff(0)
,fRecElecMCPtEtaPhiVtx(0)
,fMCElecPDG(0)
,fMCElecPtEtaPhiStrictVtx(0)
,fMCPi0Prod(0)
,fMCEtaProd(0)
,fMCPiPlusProd(0)
,fMCPiPlusProdV2(0)
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

  Int_t LPinAccBeforeEventCuts=0, LPBeforeEventCuts=0;
  if (fIsMC && fMCTrueCorrelation) MCTruthCorrelation(kFALSE, 0, 0, 0,  LPinAccBeforeEventCuts, LPBeforeEventCuts) ;


  // Get Vertex and cut > 10cm and min NumberOfTracks ( suggested by DPG to remove outliers)
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
    
  fNoEvents->Fill(0);
 

  fEventCuts.SetManualMode();
  fEventCuts.fRequireTrackVertex = true; // not in default pp cuts
  fEventCuts.fMinVtz = -10.f;
  fEventCuts.fMaxVtz = 10.f;
  fEventCuts.fMaxDeltaSpdTrackAbsolute = 0.5f;
  fEventCuts.fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
  fEventCuts.fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
  fEventCuts.fMaxResolutionSPDvertex = 0.25f;
  fEventCuts.fTriggerMask = AliVEvent::kINT7;
  fEventCuts.fRejectDAQincomplete = true;
  fEventCuts.fRequiredSolenoidPolarity = 0;

  fEventCuts.fUseMultiplicityDependentPileUpCuts= true;
  //fEventCuts.fSPDpileupMinContributors = 3;
  fEventCuts.fSPDpileupMinZdist = 0.8;
  fEventCuts.fSPDpileupNsigmaZdist = 3.;
  fEventCuts.fSPDpileupNsigmaDiamXY = 2.;
  fEventCuts.fSPDpileupNsigmaDiamZ = 5.;
  fEventCuts.fTrackletBGcut = true;

  fEventCuts.fUseVariablesCorrelationCuts = true; // maybe remove this cut, small effect
  fEventCuts.fFB128vsTrklLinearCut[0] = 32.077;
  fEventCuts.fFB128vsTrklLinearCut[1] = 0.932;

  // EventCuts
  if(!fEventCuts.AcceptEvent(fVevent)) {
    PostData(1, fOutputList);
    PostData(2, fOutputListMain);
    return;
  }
 
  fNoEvents->Fill(1);

  // old cuts which should not have any relevance
  Int_t fNOtrks = fVevent->GetNumberOfTracks();
  if (fNOtrks<2) {
    PostData(2, fOutputListMain);
    return;
  }
  fNoEvents->Fill(2);

  Double_t pVtxZ = -999.;
  pVtxZ = pVtx->GetZ();
  if(TMath::Abs(pVtxZ)>10. || TMath::Abs(spdVtx->GetZ())>10. ){
    cout << pVtxZ << "\t spd vtx " << spdVtx->GetZ() << endl;
    PostData(2, fOutputListMain);
    return;
  }
  fNoEvents->Fill(3);

  UInt_t fSelectMask = fInputHandler->IsEventSelected();
  Bool_t isINT7selected = fSelectMask & AliVEvent::kINT7;
  if (!isINT7selected){
    PostData(2, fOutputListMain);
    printf("Event not selected \n");
    return;
  }
  fNoEvents->Fill(4);

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



  if (fMinPtEvent > 0.1 || fMaxPtEvent <100) {
    if (!PassEventBias(pVtx,nMotherKink,listofmotherkink)) {
      delete [] listofmotherkink;
      PostData(2, fOutputListMain);
      return;
    }
  }
  fNoEvents->Fill(5);


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
			     

  //  SPDTracklets
  AliAODTracklets* SPDtracklets= ((AliAODEvent*)fAOD)->GetTracklets(); 
  Int_t nTr=SPDtracklets->GetNumberOfTracklets();
  Int_t nTrAcc=0;
  for(Int_t iTr=0; iTr<nTr; iTr++){
    Double_t theta=SPDtracklets->GetTheta(iTr);
    Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
    if(TMath::Abs(eta)<1.) nTrAcc++;
  }
  fSPDnTrAcc->Fill(spdVtx->GetZ(), nTrAcc);
  Int_t nTrMCAcc=0;
  Double_t SumMCHadronsPt[4]={0., 0., 0., 0.}; // >all, 0.5-2, 2-5, 5-10;
  Double_t AverageMCPt[4]={0., 0., 0., 0.};
  Int_t MCnHadrons[4]={0, 0, 0, 0};
  TObjArray* MCHadrons;
  if (fIsMC) {
    MCHadrons = new TObjArray(1000, 0);
    Double_t mcVtx[3];
    fMCheader->GetVertex(mcVtx);
    for (Int_t i=0; i<fMC->GetNumberOfTracks(); i++) {
      AliAODMCParticle *mcPart  = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(i));
      if (mcPart->Charge()==0) continue;
      if (!mcPart->IsPhysicalPrimary()) continue;
      if (TMath::Abs(mcPart->Eta())<1.) {
	nTrMCAcc++;
	MCHadrons->Add(mcPart);
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
    DiffVertexZ = spdVtx->GetZ()-mcVtx[2];
    fDiffSPDMCVtx->Fill(spdVtx->GetNContributors(), DiffVertexZ);
    fSPDnTrGen->Fill(spdVtx->GetZ(), nTrMCAcc);
  }
  for (Int_t j=0; j<4; j++)  if (MCnHadrons[j]>0) AverageMCPt[j] = SumMCHadronsPt[j]/(1.*MCnHadrons[j]);



  // Corrected SPDTracklets
  Double_t RefMinSPD, RefMaxSPD, RefMeanSPD;
  if (fIsMC) {
    RefMinSPD=8.32;
    RefMaxSPD=11.94;
    RefMeanSPD=11.07;
  }
  else {
    RefMinSPD=8.28;
    RefMaxSPD=12.02;
    RefMeanSPD=11.11;

  }
  //fSPDnTrackAvg - only temporary, adjust per run, period  mc 
  Double_t nTrAccCorrMin=AliVertexingHFUtils::GetCorrectedNtracklets(&fSPDnTrAvg,nTrAcc*1.,spdVtx->GetZ(),RefMinSPD); 
  Double_t nTrAccCorrMax=AliVertexingHFUtils::GetCorrectedNtracklets(&fSPDnTrAvg,nTrAcc*1.,spdVtx->GetZ(),RefMaxSPD); 
  Double_t nTrAccCorrMean=AliVertexingHFUtils::GetCorrectedNtracklets(&fSPDnTrAvg,nTrAcc*1.,spdVtx->GetZ(),RefMeanSPD);

    //  cout << "nTrAcc " << nTrAcc << "\t" << nTrAccCorrMin << endl;

  fSPDnTrCorrMax->Fill(spdVtx->GetZ(), 1.*nTrAccCorrMax);

  Double_t fillSparse[3]={spdVtx->GetZ(), 1.*nTrAcc,1.};
  fillSparse[2]=1.*nTrAccCorrMax;
  fnTrAccMax->Fill(fillSparse);
  fillSparse[2]=1.*nTrAccCorrMin;
  // fnTrAccMin->Fill(fillSparse);
  fillSparse[2]=1.*nTrAccCorrMean;
  // fnTrAccMean->Fill(fillSparse);


  if (fIsMC) {
    fillSparse[1]=1.*nTrMCAcc;
    fillSparse[2]=1.*nTrAccCorrMax;
    fnTrAccMaxGen->Fill(fillSparse);
    fillSparse[2]=1.*nTrAccCorrMin;
    // fnTrAccMinGen->Fill(fillSparse);
    fillSparse[2]=1.*nTrAccCorrMean;
    //  fnTrAccMeanGen->Fill(fillSparse);
  }
 

  
  // Old get Multitplicity
  Double_t fMultV0Per, fMultSPDPer, fMultV0Tot, fMultSPD;
  Float_t mult = 1;
  fMultV0Per = 1.;
  fMultSPDPer =1.;
  fMultV0Tot = 1.;
  fMultSPD = 1;
  
  fMultSelection = (AliMultSelection * ) fVevent->FindListObject("MultSelection");
  if (fMultSelection) {
    fMultV0Per = fMultSelection->GetMultiplicityPercentile("V0M", kFALSE); // Method, Embed Event selection (kFALSE)
    fMultSPDPer = fMultSelection->GetMultiplicityPercentile("SPDTracklets", kFALSE); // Method, Embed Event selection (kFALSE)
  }

  // Multiplicity estimates
  AliVVZERO* AODV0 = fVevent->GetVZEROData();
  Float_t multV0A=AODV0->GetMTotV0A();
  Float_t multV0C=AODV0->GetMTotV0C();
  fMultV0Tot=multV0A+multV0C;
  
  AliVMultiplicity* AliMult = fVevent->GetMultiplicity();
  if (AliMult) fMultSPD=AliMult->GetNumberOfTracklets();  
  // if (fIsAOD)  AliAODTracklets *AODtracklets = ((AliAODEvent*)fAOD)->GetTracklets();
  mult=fMultSPDPer;
  //Double_t fillSparse[4]={fMultV0Per, fMultSPDPer, fMultV0Tot, fMultSPD};
  //  fMultiplicity->Fill(fillSparse); 
  // fSPDMultiplicity->Fill(fMultSPDPer, fMultSPD, pVtxZ);
    
  // Efficiency Corrections
  if(fIsMC) {
    if (fIsAOD) {
      MCEfficiencyCorrections(pVtx); //  Electron reconstruction, Hadron reconstruction
      // TList *lh=fMCheader->GetCocktailHeaders();
      //      Int_t nh=lh->GetEntries();  
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

  Double_t ThrustVar[2]={-1., -1.};
  Double_t SpherVar=-1.;
  // if (Thrust(MCHadrons, ThrustVar, 1, 0.5)); // cout << ThrustVar[0] << "\t" << ThrustVar[1];
  if (fIsMC) {
    SpherVar = Sphericity(MCHadrons, 1, 0.5);
    MCHadrons->Clear();
    delete MCHadrons;
  }
  //  cout << "Sphericity " << SpherVar << endl;




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

  PhotULSLSElectronAcceptance(pVtx, mult,  nMotherKink, listofmotherkink);
  LPtrack=FindLPAndHFE(RedTracksHFE, pVtx,nMotherKink,listofmotherkink, mult, EvContainsTaggedPhot, EvContainsNonTaggedPhot);
  if (fLParticle) if (LPtrack->Pt()>=1000) return;


  if (fIsMC) {
    Double_t MotherPt=0; // MC
    for (Int_t i=0; i<RedTracksHFE->GetEntriesFast(); i++) {
    
      AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(i);
      //  cout << "TrackLabel " << RedTrack->GetLabel() << endl;
      AliAODMCParticle* MCTrack = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(RedTrack->GetLabel())); 
      if (abs(MCTrack->GetPdgCode())!=11) continue;
      Int_t MCMotherLabel = abs(MCTrack->GetMother());
      //   cout << "Mother " << MCMotherLabel << endl;
      Int_t MCGMotherLabel=-999, MCGGMotherLabel=-999;
      Int_t MCMotherPDG = -999, MCGMotherPDG = -999, MCGGMotherPDG=-999;
      Bool_t IsFeedDown=kFALSE;
      Int_t Case=-999; // Case 0 gamma pion, 1 gamma pion feed down, 2 gamma eta, 3 gammaeta feed down, 4 pion 5, pion feed deon , 6, eta, 7, eta feed down
      if (MCMotherLabel>0) {
	AliAODMCParticle *MCMother = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(MCMotherLabel));  
	MCMotherPDG = abs(MCMother->GetPdgCode());
	MCGMotherLabel = abs(MCMother->GetMother());
	MotherPt=MCMother->Pt();
	//     cout << "GMother " << MCGMotherLabel << endl;
      }
      if (MCGMotherLabel>0) {
	//     cout << "GMother " << MCGMotherLabel << endl;
	AliAODMCParticle *MCGMother = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(MCGMotherLabel));  
	MCGMotherPDG = abs(MCGMother->GetPdgCode());
	MCGGMotherLabel = abs(MCGMother->GetMother());
	//   cout << "GGMother " << MCGGMotherLabel << endl;
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
	fPt2Tagged->Fill(RedTrack->Pt(), RedTrack->TruePartnerMCPt(), LPtrack->Pt());
	
	//	fMCThrustTagged->Fill(RedTrack->Pt(), ThrustVar[0], Case);
	fMCSpherTagged->Fill(RedTrack->Pt(), SpherVar, Case);
	if (abs(RedTrack->ID())!=abs(LPtrack->GetID())) fRecLPTagged->Fill(RedTrack->Pt(), LPtrack->Pt(), Case);
	//	fMultCorrTagged->Fill(RedTrack->Pt(), nTrAccCorrMax, Case);
	fNHadTagged->Fill(RedTrack->Pt(), MCnHadrons[0], Case);
	//	fNHadTaggedA->Fill(RedTrack->Pt(), MCnHadrons[1], Case);
	fNHadTaggedB->Fill(RedTrack->Pt(), MCnHadrons[2], Case);
	//fNHadTaggedC->Fill(RedTrack->Pt(), MCnHadrons[3], Case);
	fMeanPtTagged->Fill(RedTrack->Pt(), AverageMCPt[0], Case);
	//	fMeanPtTaggedA->Fill(RedTrack->Pt(), AverageMCPt[1], Case);
	fMeanPtTaggedB->Fill(RedTrack->Pt(), AverageMCPt[2], Case);
	//fMeanPtTaggedC->Fill(RedTrack->Pt(), AverageMCPt[3], Case);
	



	//	fMothMCThrustTagged->Fill(MotherPt, ThrustVar[0]);
	fMothMCSpherTagged->Fill(MotherPt, SpherVar);
	if (abs(RedTrack->ID())!=abs(LPtrack->GetID())) fMothRecLPTagged->Fill(MotherPt, LPtrack->Pt());
	//	fMothMultCorrTagged->Fill(MotherPt, nTrAccCorrMax);
	fMothNHadTagged->Fill(MotherPt, MCnHadrons[0]);
	fMothMeanPtTagged->Fill(MotherPt, AverageMCPt[0]);
	if (LPtrack->Pt()>2) {
	  //  fMCThrustTaggedH->Fill(RedTrack->Pt(), ThrustVar[0]);
	  fMCSpherTaggedH->Fill(RedTrack->Pt(), SpherVar);
	  if (abs(RedTrack->ID())!=abs(LPtrack->GetID())) fRecLPTaggedH->Fill(RedTrack->Pt(), LPtrack->Pt());
	  //	  fMultCorrTaggedH->Fill(RedTrack->Pt(), nTrAccCorrMax);
	  fNHadTaggedH->Fill(RedTrack->Pt(), MCnHadrons[0]);
	  fMeanPtTaggedH->Fill(RedTrack->Pt(), AverageMCPt[0]);
	  //  fMothMCThrustTaggedH->Fill(MotherPt, ThrustVar[0]);
	  fMothMCSpherTaggedH->Fill(MotherPt, SpherVar);
	  if (abs(RedTrack->ID())!=abs(LPtrack->GetID())) fMothRecLPTaggedH->Fill(MotherPt, LPtrack->Pt());
	  // fMothMultCorrTaggedH->Fill(MotherPt, nTrAccCorrMax);
	  fMothNHadTaggedH->Fill(MotherPt, MCnHadrons[0]);
	  fMothMeanPtTaggedH->Fill(MotherPt, AverageMCPt[0]);
	
	}
      }
      else {

	fPt2NTagged->Fill(RedTrack->Pt(), RedTrack->TruePartnerMCPt(), LPtrack->Pt());
	
	//	fMCThrustNTagged->Fill(RedTrack->Pt(), ThrustVar[0], Case);
	fMCSpherNTagged->Fill(RedTrack->Pt(), SpherVar, Case);
	if (abs(RedTrack->ID())!=abs(LPtrack->GetID())) fRecLPNTagged->Fill(RedTrack->Pt(), LPtrack->Pt(), Case);
	//	fMultCorrNTagged->Fill(RedTrack->Pt(), nTrAccCorrMax, Case);
	fNHadNTagged->Fill(RedTrack->Pt(), MCnHadrons[0], Case);
	//	fNHadNTaggedA->Fill(RedTrack->Pt(), MCnHadrons[1], Case);
	fNHadNTaggedB->Fill(RedTrack->Pt(), MCnHadrons[2], Case);
	//	fNHadNTaggedC->Fill(RedTrack->Pt(), MCnHadrons[3], Case);
	fMeanPtNTagged->Fill(RedTrack->Pt(), AverageMCPt[0], Case);
	//	fMeanPtNTaggedA->Fill(RedTrack->Pt(), AverageMCPt[1], Case);
	fMeanPtNTaggedB->Fill(RedTrack->Pt(), AverageMCPt[2], Case);
	//	fMeanPtNTaggedC->Fill(RedTrack->Pt(), AverageMCPt[3], Case);




	//	fMothMCThrustNTagged->Fill(MotherPt, ThrustVar[0]);
	fMothMCSpherNTagged->Fill(MotherPt, SpherVar);
	if (abs(RedTrack->ID())!=abs(LPtrack->GetID())) 	fMothRecLPNTagged->Fill(MotherPt, LPtrack->Pt());
	//	fMothMultCorrNTagged->Fill(MotherPt, nTrAccCorrMax);
	fMothNHadNTagged->Fill(MotherPt, MCnHadrons[0]);
	fMothMeanPtNTagged->Fill(MotherPt, AverageMCPt[0]);
	if (LPtrack->Pt()>2) {
	  //  fMCThrustNTaggedH->Fill(RedTrack->Pt(), ThrustVar[0]);
	  fMCSpherNTaggedH->Fill(RedTrack->Pt(), SpherVar);
	  if (abs(RedTrack->ID())!=abs(LPtrack->GetID())) fRecLPNTaggedH->Fill(RedTrack->Pt(), LPtrack->Pt());
	  //  fMultCorrNTaggedH->Fill(RedTrack->Pt(), nTrAccCorrMax);
	  fNHadNTaggedH->Fill(RedTrack->Pt(), MCnHadrons[0]);
	  fMeanPtNTaggedH->Fill(RedTrack->Pt(), AverageMCPt[0]);
	  //  fMothMCThrustNTaggedH->Fill(MotherPt, ThrustVar[0]);
	  fMothMCSpherNTaggedH->Fill(MotherPt, SpherVar);
	  if (abs(RedTrack->ID())!=abs(LPtrack->GetID())) fMothRecLPNTaggedH->Fill(MotherPt, LPtrack->Pt());
	  //  fMothMultCorrNTaggedH->Fill(MotherPt, nTrAccCorrMax);
	  fMothNHadNTaggedH->Fill(MotherPt, MCnHadrons[0]);
	  fMothMeanPtNTaggedH->Fill(MotherPt, AverageMCPt[0]);
	}
      }
    }
  }


  /*

  if (EvContainsTaggedPhot & !EvContainsNonTaggedPhot)  fCheckTaggedEvent->Fill(1.,LPtrack->Pt(), mult);
  else if (EvContainsNonTaggedPhot & !EvContainsTaggedPhot) fCheckTaggedEvent->Fill(2., LPtrack->Pt(), mult);
  else if (EvContainsTaggedPhot && EvContainsNonTaggedPhot) fCheckTaggedEvent->Fill(0., LPtrack->Pt(), mult);
  */

  Int_t LPinAccAfterEventCuts=-999, LPAfterEventCuts=-999;
  AliAODTrack * LPtrackAOD = dynamic_cast<AliAODTrack*>(LPtrack);
  Int_t LPtrackLabel=0;
  if (LPtrackAOD) LPtrackLabel=abs(LPtrackAOD->GetLabel());
  if (fIsMC && fMCTrueCorrelation) {  
    MCTruthCorrelation(kTRUE, LPtrackLabel, pVtx->GetZ(), mult,  LPinAccAfterEventCuts, LPAfterEventCuts) ;
  
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
    if (RedTracksHFE->GetEntriesFast()>0) CorrelateLP(LPtrack,pVtx, nMotherKink, listofmotherkink, RedTracksHFE);
    if (fMixedEvent) CorrelateLPMixedEvent(LPtrack, mult, pVtx->GetZ(), LPtrack->Pt(),  EvContainsTaggedPhot, EvContainsNonTaggedPhot); // condition that electron track is in event has been removed!
  }

  // Hadron - only one function for both options, as only one loop over Hadron tracks
  // Mixed event is called inside this function
  if (fCorrHadron && fLParticle) {
    if (RedTracksHFE->GetEntriesFast()>0) CorrelateHadron(RedTracksHFE, pVtx, nMotherKink, listofmotherkink, mult, LPtrack->Pt());
    if (fMixedEvent) CorrelateHadronMixedEvent( mult, pVtx, LPtrack->Pt(), nMotherKink, listofmotherkink,  EvContainsTaggedPhot, EvContainsNonTaggedPhot);
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


  fesdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 

  // V0 Kine cuts 
  fV0cuts = new AliAODv0KineCuts();
  fV0cutsESD = new AliESDv0KineCuts();

  // V0 PID Obj arrays
  fV0electrons = new TObjArray;
  fV0pions     = new TObjArray;
  fV0protons   = new TObjArray;

  fNoEvents = new TH1F("fNoEvents","",8,-0.5,7.5);
  fOutputListMain->Add(fNoEvents);
    
  fTrkpt = new TH2F("fTrkpt","track pt,after EleTrack, after ElePID",200,0,20,3,-0.5,2.5);
  fOutputListMain->Add(fTrkpt);

  fEtaVtxZ = new TH2F("fEtaVtxZ", "Eta vs VtxZ after hadron track cuts", 90, -0.9, 0.9, 100, -10, 10);
  fOutputListMain->Add(fEtaVtxZ);

  PdgTable =  new TDatabasePDG();
  PdgTable->ReadPDGTable();
  const THashList* PDGParticleList = PdgTable->ParticleList();
  // cout << "PDGSize " << PDGParticleList->GetEntries() << endl;
  Int_t PDGSize = PDGParticleList->GetEntries()+20;

  TIter next(PDGParticleList);
  TParticlePDG *p;
  Int_t Bin=0;
  std::vector<TString> PDGLabel;
  PDGLabel.push_back("Empty");
  while ((p = (TParticlePDG *)next())) {
    Bin++;
    PDGMap[p->PdgCode()] = Bin;
    PDGLabel.push_back(Form("%i",p->PdgCode()));
    //  cout << p->PdgCode() << "\t" << Bin << "\t" << PDGLabel[Bin] << endl;
  }
  //  cout << PDGMap.find(-4)->second << endl;
 


    

  // General  Binning 		    
  // fAssPtHad_Nbins = 15;
  const Int_t TmpHad_Nbins = 15;
  fAssPtHad_Nbins = TmpHad_Nbins;
  Float_t TmpArrayLow[15]= {0,  0.5, 0.5, 0.5, 1., 1., 1., 2,  2, 5, 5,  5, 10, 10, 15};
  Float_t TmpArrayUp[15]={ 999, 2.0, 5.0, 999, 2., 5.,999, 5,999,10,15,999, 15,999,999};
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
  Double_t XminElectron=0.25;
  Double_t XmaxElectron=4.;
  const Int_t    NBinsElectronRed = 19;
  Double_t XBinsElectronRed[]={0.25,0.5,0.75, 1., 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4, 5, 6, 10, 100};

  Int_t     NBinsHadron=200 ;
  Double_t  XminHadron=0.0;
  Double_t  XmaxHadron=100;
  const Int_t   NBinsHadRed=8;
  Double_t  XBinsHadRed[]={0., 0.5, 1., 2., 5., 10, 15, 50, 100}; 

  const  Int_t     NMultBins=4;
  Double_t    XMultBins[]={0,25,50,75,100};
 
  const Int_t NVertexBins = 8;
  Double_t XVertexBins[]={-10,-4.32,-2.32,-0.82,0.52,1.87,3.36,5.31, 10};  // Quantile


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
  fSPDVtxRes = new TH2F("fSPDVtxRes", "fSPDVtxRes",  300, -0.5, 299.5, 500, 0., 10.);
  fOutputListQA->Add(fSPDVtxRes);
 
  fDiffSPDPrimVtx = new TH2F("fDiffSPDPrimVtx", "fDiffSPDPrimVtx", 300, -0.5, 299.5, 500, -5., 5.);
  fOutputListQA->Add(fDiffSPDPrimVtx);

  fSPDnTrAcc = new TH2F("fSPDnTrAcc", "fSPDnTrAcc", 220, -11, 11, 300, -0.5, 299.5);
  fOutputListQA->Add(fSPDnTrAcc);

  fSPDnTrCorrMax = new TH2F("fSPDnTrCorrMax", "fSPDnTrCorrMax", 220, -11,11,300, -0.5, 299.5);
  fOutputListQA->Add(fSPDnTrCorrMax);

  Int_t    nBinsSPD[3]={10, 200, 200};
  Double_t xminSPD[3]={-10,-0.5,-0.5};
  Double_t xmaxSPD[3]={10,199.5, 199.5};

  if (fIsMC) {
    fSPDnTrGen = new TH2F("fSPDnTrGen", "fSPDnTrGen", 220, -11, 11, 300, -0.5, 299.5);
    fOutputListQA->Add(fSPDnTrGen);

    fDiffSPDMCVtx = new TH2F("fDiffSPDMCVtx", "fDiffSPDMCVtx", 200, -0.5, 199.5, 500, -5., 5.);
    fOutputListQA->Add(fDiffSPDMCVtx);

    fnTrAccMaxGen = new THnSparseF("fnTrAccMaxGen", "fnTrAccMaxGen", 3, nBinsSPD, xminSPD, xmaxSPD);
    fOutputListQA->Add(fnTrAccMaxGen);


    // fnTrAccMinGen = new THnSparseF("fnTrAccMinGen", "fnTrAccMinGen", 3, nBinsSPD, xminSPD, xmaxSPD);
    //  fOutputListQA->Add(fnTrAccMinGen);
  
    //  fnTrAccMeanGen = new THnSparseF("fnTrAccMenGen", "fnTrAccMeanGen",  3, nBinsSPD, xminSPD, xmaxSPD);
    //  fOutputListQA->Add(fnTrAccMeanGen);


    fPt2Tagged = new TH3F("fPt2Tagged", "Pt2Tagged; RecPtE; MCPt2; RecLP", 20, 0.5,3,60, 0.,3., 20,0.,10.);
    fOutputListMain->Add(fPt2Tagged);

    fPt2NTagged = new TH3F("fPt2NTagged", "Pt2NTagged; RecPtE; MCPt2; RecLP", 20, 0.5,3,60, 0.,3., 20,0.,10.);
    fOutputListMain->Add(fPt2NTagged);

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

    fMeanPtTaggedB = new TH3F("fMeanPtTaggedB", "Mean pt for tagged; ptElectron; average pt", 44, 0.5, 6., 100, 0, 10, 8, -0.5, 7.5);
    fOutputListQA->Add(fMeanPtTaggedB);

    //fMeanPtTaggedC = new TH3F("fMeanPtTaggedC", "Mean pt for tagged; ptElectron; average pt", 44, 0.5, 6., 100, 0, 10, 8, -0.5, 7.5);
    // fOutputListQA->Add(fMeanPtTaggedC);

    // fMCThrustNTagged = new TH3F("fMCThrustNTagged", "Thrust of NTagged; ptElectron; Thrust", 44, 0.5, 6., 100, 0., 1., 8, -0.5, 7.5);
    // fOutputListQA->Add(fMCThrustNTagged);

    fMCSpherNTagged = new TH3F("fMCSpherNTagged", "Sphericity of NTagged; ptElectron; Thrust", 44, 0.5, 6., 100, 0., 1., 8, -0.5, 7.5);
    fOutputListQA->Add(fMCSpherNTagged);

    fRecLPNTagged = new TH3F("fRecLPNTagged", "LP of NTagged; ptElectron; ptLP", 44, 0.5, 6., 200, 0, 20, 8, -0.5, 7.5);
    fOutputListQA->Add(fRecLPNTagged);

    // fMultCorrNTagged = new TH3F("fMultCorrNTagged", "Corr SPD mult of tagged; ptElectron; nTr", 44, 0.5, 6., 150, -0.5, 149.5, 8, -0.5, 7.5);
    //fOutputListQA->Add(fMultCorrNTagged);

    fNHadNTagged = new TH3F("fNHadNTagged", "NHadrons for tagged; ptElectron; nHad", 44, 0.5, 6., 150, -0.5, 149.5, 8, -0.5, 7.5);
    fOutputListQA->Add(fNHadNTagged);

    //  fNHadNTaggedA = new TH3F("fNHadNTaggedA", "NHadrons for tagged; ptElectron; nHad", 44, 0.5, 6., 150, -0.5, 149.5, 8, -0.5, 7.5);
    // fOutputListQA->Add(fNHadNTaggedA);

    fNHadNTaggedB = new TH3F("fNHadNTaggedB", "NHadrons for tagged; ptElectron; nHad", 44, 0.5, 6., 150, -0.5, 149.5, 8, -0.5, 7.5);
    fOutputListQA->Add(fNHadNTaggedB);

    //fNHadNTaggedC = new TH3F("fNHadNTaggedC", "NHadrons for tagged; ptElectron; nHad", 44, 0.5, 6., 150, -0.5, 149.5, 8, -0.5, 7.5);
    //fOutputListQA->Add(fNHadNTaggedC);

    fMeanPtNTagged = new TH3F("fMeanPtNTagged", "Mean pt for tagged; ptElectron; average pt", 44, 0.5, 6., 100, 0, 10, 8, -0.5, 7.5);
    fOutputListQA->Add(fMeanPtNTagged);

    // fMeanPtNTaggedA = new TH3F("fMeanPtNTaggedA", "Mean pt for tagged; ptElectron; average pt", 44, 0.5, 6., 100, 0, 10, 8, -0.5, 7.5);
    //   fOutputListQA->Add(fMeanPtNTaggedA);

    fMeanPtNTaggedB = new TH3F("fMeanPtNTaggedB", "Mean pt for tagged; ptElectron; average pt", 44, 0.5, 6., 100, 0, 10, 8, -0.5, 7.5);
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

  fnTrAccMax = new THnSparseF("fnTrAccMax", "fnTrAccMax", 3, nBinsSPD, xminSPD, xmaxSPD);
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

  fElectronTrackTPCNcls = new TH2F("fElectronTrackTPCNcls", "fElectronTrackTPCNcls", 20, 0, 10,170, 0, 170);
  fOutputListQA->Add(fElectronTrackTPCNcls);

  fElectronTrackTPCNclsdEdx= new TH2F("fElectronTrackTPCNclsdEdx", "fElectronTrackTPCNclsdEdx",20, 0, 10, 170, 0, 170);
  fOutputListQA->Add(fElectronTrackTPCNclsdEdx);

  fElectronTrackTPCFrac = new TH2F("fElectronTrackTPCFrac", "fElectronTrackTPCFrac", 20, 0, 10,110, 0, 1.1);
  fOutputListQA->Add(fElectronTrackTPCFrac);
  
  fElectronTrackITSNcls = new TH2F("fElectronTrackITSNcls", "fElectronTrackITSNcls", 20, 0, 10,10, 0, 10);
  fOutputListQA->Add(fElectronTrackITSNcls);

  fElectronTrackITSLayer = new TH3F("fElectronTrackITSLayer", "fElectronTrackITSLayer", 20, 0, 10,2, -0.5, 1.5, 2, -0.5, 1.5);
  fOutputListQA->Add(fElectronTrackITSLayer);
 
  fElectronTrackRefit = new TH3F("fElectronTrackRefit", "fElectronTrackReft", 20, 0, 10,2, -0.5, 1.5, 2, -0.5, 1.5);
  fOutputListQA->Add(fElectronTrackRefit);

  fElectronTrackDCA = new TH3F("fElectronTrackDCA", "fElectronTrackDCA: pt, r, z", 10, 0, 5, 60, -3.0, 3.0, 60, -3.0, 3.0);
  fOutputListQA->Add(fElectronTrackDCA);

  /*
  Int_t BinsITSCuts[5]={9, 7, 7, 20, 8};
  Double_t BinsITSMin[5]={0.0,-0.5,-0.5,0.,0.};
  Double_t BinsITSMax[5]={4.5, 6.5,6.5,20.,399.};

  fElectronTrackITSCuts = new THnSparseF("fElectronTrackITSCuts", "fElectronTrackITSCuts; Pt; NCls; SharedCls; X2/NCls; ParticleID", 5, BinsITSCuts, BinsITSMin, BinsITSMax);
  Double_t ParticleIDBins[9]={0.,4.5, 5.5,20., 110.,220.,222.,224.,399.};
  fElectronTrackITSCuts->GetAxis(4)->Set(BinsITSCuts[4], ParticleIDBins);
  fOutputListQA->Add(fElectronTrackITSCuts);

 
  Int_t BinsPITSCuts[5]={30, 7, 7, 20, 2};
  Double_t BinsPITSMin[5]={0.0,-0.5,-0.5,0.,-0.5};
  Double_t BinsPITSMax[5]={1.5, 6.5,6.5,20.,1.5};

  fPhotTrackITSCuts = new THnSparseF("fPhotTrackITSCuts", "fPhotTrackITSCuts; Pt; NCls; SharedCls; X2/NCls; IsPhotonic", 5, BinsPITSCuts, BinsPITSMin, BinsPITSMax);
  fOutputListQA->Add(fPhotTrackITSCuts);
  */
 

  fHadronTrackTPCNcls = new TH2F("fHadronTrackTPCNcls", "fHadronTrackTPCNcls", 20, 0, 10,170, 0, 170);
  fOutputListQA->Add(fHadronTrackTPCNcls);

  fHadronTrackRefit = new TH3F("fHadronTrackRefit", "fHadronTrackReft", 20, 0, 10,2, -0.5, 1.5, 2, -0.5, 1.5);
  fOutputListQA->Add(fHadronTrackRefit);

  fHadronTrackDCA = new TH3F("fHadronTrackDCA", "fHadronTrackDCA: pt, r, z", 20, 0, 10,  60, -3.0, 3.0, 60, -3.0, 3.0);
  fOutputListQA->Add(fHadronTrackDCA);

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

    // HadContSparse
    Int_t    binHC[4] =  {NBinsElectron  ,NBinsPhi/8      ,NBinsEta/2    ,200}; //p, Phi, Eta, TPC
    Double_t xminHC[4] = {XminElectron   ,0             ,-0.9  ,-10};
    Double_t xmaxHC[4] = {XmaxElectron   ,TMath::TwoPi(), 0.9  ,10};
    
    fHadContPvsPt = new TH2F("fHadContPvsPt", "PvsPt; P;Pt", 100, 0, 10, 100, 0, 10);
    fOutputListMain->Add(fHadContPvsPt);  
    
    fHadContPPhiEtaTPC = new THnSparseF("fHadContPPhiEtaTPC", "HadCont; P; Phi; Eta; TPC", 4, binHC, xminHC, xmaxHC);
    fOutputListMain->Add(fHadContPPhiEtaTPC);
    
    Int_t    binHC2[4] =  {NBinsElectron  ,10   , 20 ,100}; //p, ITS, TOF, TPC
    Double_t xminHC2[4] = {XminElectron   ,-5  ,-10 ,-10};
    Double_t xmaxHC2[4] = {XmaxElectron   ,5   , 10 ,10};  
    fHadContamination = new THnSparseF("fHadContamination", "HadCont; P; ITS; TOF; TPC", 4, binHC2, xminHC2, xmaxHC2);
    fOutputListMain->Add(fHadContamination);

    //    fHadContaminationPt = new THnSparseF("fHadContaminationPt", "HadCont: Pt, ITS, TOF, TPC", 4, binHC2, xminHC2, xmaxHC2);
    //  fOutputListMain->Add(fHadContaminationPt);
    
    fHadContTPCEtaPhiPt = new TH3F("fHadContTPCEtaPhiPt", "MCHadContTPC; Eta; Phi; Pt", 36, -0.9, 0.9, 32, 0, TMath::TwoPi(), NBinsElectron, XminElectron, XmaxElectron);
    fOutputListMain->Add(fHadContTPCEtaPhiPt);
  
    if (fIsMC) {
      fHadContEtaPhiPt = new TH3F("fHadContEtaPhiPt", "MCHadCont; Eta; Phi; Pt", 36, -0.9, 0.9, 32, 0, TMath::TwoPi(), NBinsElectron, XminElectron, XmaxElectron);
      fOutputListMain->Add(fHadContEtaPhiPt);   
      
      Int_t    binHMC[4] = {NBinsElectron  , 7,  20, 100}; // p PDG ITS, TPC
      Double_t xminHMC[4] = {XminElectron   , 0, -10, -10};
      Double_t xmaxHMC[4] = {XmaxElectron   , 7,  10,  10};
      fHadContMC = new THnSparseF("fHadContMC", "HadContMC; P; PDG; ITS; TPC", 4, binHMC, xminHMC, xmaxHMC);
      fOutputListMain->Add(fHadContMC);
      // fHadContMCPt = new THnSparseF("fHadContMCPt", "HadCont: Pt, PDG, ITS, TPC", 4, binHMC, xminHMC, xmaxHMC);
      // fOutputListMain->Add(fHadContMCPt);
    }    
  }

  fInclElecPtEta = new TH2F("fInclElePtEta", "fInclElePtEta; #it{p}_{T}; #eta",NBinsElectron,XminElectron, XmaxElectron, NBinsEta, -0.9, 0.9 );
  fOutputListMain->Add(fInclElecPtEta);

  fInclElecPtEtaWRecEff = new TH2F("fInclElePtEtaWRecEff", "fInclElePtEtaWRecEff; #it{p}_{T}; #eta",NBinsElectron,XminElectron, XmaxElectron,  NBinsEta, -0.9, 0.9 );
  fOutputListMain->Add(fInclElecPtEtaWRecEff);

  // fInclElecP = new TH1F("fInclEleP", "fInclEleP",NBinsElectronRed, XBinsElectronRed);
  // fOutputListMain->Add(fInclElecP);

  fULSElecPt = new TH1F("fULSElePt", "fULSElePt; #it{p}_{T}",NBinsElectron,XminElectron, XmaxElectron);
  fOutputListMain->Add(fULSElecPt);
    
  fULSElecPtWRecEff = new TH1F("fULSElePtWRecEff", "fULSElePtWRecEff; #it{p}_{T}",NBinsElectron,XminElectron, XmaxElectron);
  fOutputListMain->Add(fULSElecPtWRecEff);
    
  fLSElecPt = new TH1F("fLSElePt", "fLSElePt; #it{p}_{T}",NBinsElectron,XminElectron, XmaxElectron);
  fOutputListMain->Add(fLSElecPt);

  fLSElecPtWRecEff = new TH1F("fLSElePtWRecEff", "fLSElePtWRecEff; #it{p}_{T}",NBinsElectron,XminElectron, XmaxElectron);
  fOutputListMain->Add(fLSElecPtWRecEff);
    
    
  fInvmassULS = new TH2F("fInvmassULS", "fInvmassULS", 500,0,0.5,NBinsElectronRed, XBinsElectronRed);
  fOutputListMain->Add(fInvmassULS);
    
  fInvmassLS = new TH2F("fInvmassLS", "fInvmassLS", 500,0,0.5,NBinsElectronRed, XBinsElectronRed);
  fOutputListMain->Add(fInvmassLS);


  fPhotMixULS = new TH1F("fPhotMixULS", "fPhotMixULS; #it{p}_{T}", NBinsElectron,XminElectron, XmaxElectron);
  fOutputListMain->Add(fPhotMixULS);

  fPhotMixLS = new TH2F("fPhotMixLS", "fPhotMixcLS; #it{p}_{T}; Charge", NBinsElectron,XminElectron, XmaxElectron, 3, -1.5, 1.5);
  fOutputListMain->Add(fPhotMixLS);




  fPhotPt1Pt2Only = new TH2F("fPhotPt1Pt2Only", "fPhotPt1Pt2 for tagged; Pt1Rec; Pt2Rec", 48, 0., 6., 48, 0., 6.);
  fOutputListMain->Add(fPhotPt1Pt2Only);

  if (fIsMC) {

    //    fRecMCInvMass = new TH3F("fRecMCInvMass", "fRecMCInvMass; MCPt; rec. inv. mass; MC inv. mass", 6, 0., 3., 14, 0, 0.14, 18, 0, 0.18);
    // fOutputListMain->Add(fRecMCInvMass);

    fInvmassMCTrue = new TH3F("fInvmassMCTrue", "fInvmassMCTrue; mee; pt; pdg", 100,0,0.5,24, 0, 6, PDGSize, 0.5, PDGSize+0.5);
    SetPDGAxis(fInvmassMCTrue->GetZaxis(), PDGLabel);
    fOutputListMain->Add(fInvmassMCTrue);

    Int_t    PhotBins[5]={10 , 120, 48 ,  4, 2};
    Double_t PhotXMin[5]={0.5,-3., 0.0, -0.5, -0.5};
    Double_t PhotXMax[5]={3. , 3., 12.,  3.5, 1.5};

    fPhotPt1PtMTag = new TH2F("fPhotPt1PtMTag", "fPhotPt1PtMTag; MCPt_e; MCPt_mother", 50, 0.5,3., 120, 0., 12.);
    fOutputListMain->Add(fPhotPt1PtMTag);
    fPhotPt1PtMNTag = new TH2F("fPhotPt1PtMNTag", "fPhotPt1PtMTag; MCPt_e; MCPt_mother", 50, 0.5,3., 120, 0., 12.);
    fOutputListMain->Add(fPhotPt1PtMNTag);
    
    Double_t Pt2CorrBins[5]={0., 0.2, 0.35, 2., 10.};
    fPhotPt2MCRec = new TH2F("fPhotPt2MCRec", "fPhotPt2MCRec; Pt2Rec; Pt2MC", 4, Pt2CorrBins, 4, Pt2CorrBins); //48, 0., 6., 48, 0., 6.);
    fOutputListMain->Add(fPhotPt2MCRec);

    fPhotPt1Pt2 = new THnSparseF("fPhotPt1Pt2", "fPhotPt1Pt2; Pt1Rec; Pt2MC; PtMotherMC; IdMother; PrimSec",  5, PhotBins, PhotXMin, PhotXMax);
    fOutputListMain->Add(fPhotPt1Pt2);
    fPhotPt1Pt2Corr = new THnSparseF("fPhotPt1Pt2Corr", "fPhotPt1Pt2COrr; Pt1Rec; Pt2MC; PtMotherMC; IdMother; PrimSec",  5, PhotBins, PhotXMin, PhotXMax);
    fOutputListMain->Add(fPhotPt1Pt2Corr);
    fPhotPt1Pt2MC = new THnSparseF("fPhotPt1Pt2MC", "fPhotPt1Pt2MC; Pt1Rec; Pt2MC; PtMotherMC; IdMother; PrimSec",  5, PhotBins, PhotXMin, PhotXMax);
    fOutputListMain->Add(fPhotPt1Pt2MC);
  
    /*
    fPhotPt1RecPt2 = new THnSparseF("fPhotPt1RecPt2", "fPhotPt1RecPt2; Pt1Rec; Pt2; PtMother; IdMother; PrimSec",  5, PhotBins, PhotXMin, PhotXMax);
    fOutputListMain->Add(fPhotPt1RecPt2);
    fPhotPt1RecPt2Corr = new THnSparseF("fPhotPt1RecPt2Corr", "fPhotPt1RecPt2COrr; Pt1Rec; Pt2; PtMother; IdMother; PrimSec",  5, PhotBins, PhotXMin, PhotXMax);
    fOutputListMain->Add(fPhotPt1RecPt2Corr);

    */

    //  fPhotPt1Pt2Rec = new THnSparseF("fPhotPt1Pt2Rec", "fPhotPt1Pt2Rec; Pt1Rec; Pt2Rec; PtMotherMC; IdMother; PrimSec",  5, PhotBins, PhotXMin, PhotXMax);
    // fOutputListMain->Add(fPhotPt1Pt2Rec);

    
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
   
  }

  fOpeningAngleULS = new TH2F("fOpeningAngleULS","fOpeningAngleULS",100,0,1,NBinsElectron,XminElectron,XmaxElectron);
  fOutputListMain->Add(fOpeningAngleULS);

  fOpeningAngleLS = new TH2F("fOpeningAngleLS","fOpeningAngleLS",100,0,1,NBinsElectron,XminElectron,XmaxElectron);
  fOutputListMain->Add(fOpeningAngleLS);
    
  fCheckLSULS = new TH2F("fCheckLSULS", "LSULS",5,0,5,5,0,5);
  fOutputListMain->Add(fCheckLSULS);
   
  fTagEtaPt1Pt2 = new TH3F("fTagEtaPt1Pt2", "Tagged; Eta; Pt1; Pt2", 36, -0.9, 0.9, NBinsElectron/2, XminElectron, XmaxElectron, 30, 0, 1.5);
  fOutputList->Add(fTagEtaPt1Pt2);

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

  if (fIsMC && fTagEff) {

    fNonTagEtaPt1Pt2 = new TH3F("fNonTagEtaPt1Pt2", "NonTagged; Eta; Pt1; Pt2", 36, -0.9, 0.9, NBinsElectron/2, XminElectron, XmaxElectron, 30, 0, 1.5);
    fOutputList->Add(fNonTagEtaPt1Pt2);
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
    
    fTagMotherPt = new THnSparseF("fTagMotherPt", "Incl; ptElectron; ptMother; PDGMother; PDGGrandmother",4, binMothPt, xminMothPt, xmaxMothPt);
    SetPDGAxis(fTagMotherPt->GetAxis(2), PDGLabel);
    SetPDGAxis(fTagMotherPt->GetAxis(3), PDGLabel);
    fOutputListMain->Add(fTagMotherPt);
    /*
    Int_t    binTagEff[5] =  {NBinsElectron   ,10000 ,10001, 10001, 10001}; //p, pdg, pdgmother
    Double_t xminTagEff[5] = {XminElectron    ,-0.5   , -1.5,-1.5, -1.5};
    Double_t xmaxTagEff[5] = {XmaxElectron    ,9999.5, 9999.5,9999.5, 9999.5};  
    */
    // PtElectron, ElectronOrHadron, ElectronMother (Gamma, Pion, Eta, Other), ElectronGrandMother( Eta-Gamma, Pion-Gamma, Other), HFGGMother).
    Int_t    binTagEff[5] =  {NBinsElectron   ,2 ,PDGSize, PDGSize, 2}; //p, pdg, pdgmother
    Double_t xminTagEff[5] = {XminElectron    ,-0.5   , 0.5,0.5, -0.5};
    Double_t xmaxTagEff[5] = {XmaxElectron    ,1.5, PDGSize+0.5,PDGSize+0.5, 1.5};  

    fTagEffIncl = new THnSparseF("fTagEffIncl", "Incl tag.eff.; Pt; IsElectron; PDGMother; PDGGMother; PDGGGMIsNonHF", 5, binTagEff, xminTagEff, xmaxTagEff);
    //    SetPDGAxis(fTagEffIncl->GetAxis(1), PDGLabel);
    SetPDGAxis(fTagEffIncl->GetAxis(2), PDGLabel);
    SetPDGAxis(fTagEffIncl->GetAxis(3), PDGLabel);
    // SetPDGAxis(fTagEffIncl->GetAxis(4), PDGLabel);
    fOutputListMain->Add(fTagEffIncl);
    
    /*
    fTagEffLS = new THnSparseF("fTagEffLS", "LS tag.eff.:  Pt; IsElectron; PDGMother; PDGGMother; PDGGGMIsNonHF", 5, binTagEff, xminTagEff, xmaxTagEff);
    // SetPDGAxis(fTagEffLS->GetAxis(1), PDGLabel);
    SetPDGAxis(fTagEffLS->GetAxis(2), PDGLabel);
    SetPDGAxis(fTagEffLS->GetAxis(3), PDGLabel);
    // SetPDGAxis(fTagEffLS->GetAxis(4), PDGLabel);
    fOutputListMain->Add(fTagEffLS);
    */

    fTagEffULS = new THnSparseF("fTagEffULS", "ULS-LS tag.eff.;  Pt; IsElectron; PDGMother; PDGGMother; PDGGGMIsNonHF", 5, binTagEff, xminTagEff, xmaxTagEff);
    //    SetPDGAxis(fTagEffULS->GetAxis(1), PDGLabel);
    SetPDGAxis(fTagEffULS->GetAxis(2), PDGLabel);
    SetPDGAxis(fTagEffULS->GetAxis(3), PDGLabel);
    //    SetPDGAxis(fTagEffULS->GetAxis(4), PDGLabel);
    fOutputListMain->Add(fTagEffULS);

    fTagTruePairs = new THnSparseF("fTagTruePairs", "ULS true pairs tag.eff.;  Pt; IsElectron; PDGMother; PDGGMother; PDGGGMIsNonHF", 5, binTagEff, xminTagEff, xmaxTagEff);
    //    SetPDGAxis(fTagTruePairs->GetAxis(1), PDGLabel);
    SetPDGAxis(fTagTruePairs->GetAxis(2), PDGLabel);
    SetPDGAxis(fTagTruePairs->GetAxis(3), PDGLabel);
    //    SetPDGAxis(fTagTruePairs->GetAxis(4), PDGLabel);
    fOutputListMain->Add(fTagTruePairs);

    fTagEffInclWoWeight = new THnSparseF("fTagEffInclWoWeight", "Incl tag.eff. woW;  Pt; IsElectron; PDGMother; PDGGMother; PDGGGMIsNonHF", 5, binTagEff, xminTagEff, xmaxTagEff);
    //    SetPDGAxis(fTagEffInclWoWeight->GetAxis(1), PDGLabel);
    SetPDGAxis(fTagEffInclWoWeight->GetAxis(2), PDGLabel);
    SetPDGAxis(fTagEffInclWoWeight->GetAxis(3), PDGLabel);
    //    SetPDGAxis(fTagEffInclWoWeight->GetAxis(4), PDGLabel);
    fOutputListMain->Add(fTagEffInclWoWeight);
  
    /*
    fTagEffLSWoWeight = new THnSparseF("fTagEffLSWoWeight", "LS tag.eff. woW;  Pt; IsElectron; PDGMother; PDGGMother; PDGGGMIsNonHF", 5, binTagEff, xminTagEff, xsfmaxTagEff);
    //    SetPDGAxis(fTagEffLSWoWeight->GetAxis(1), PDGLabel);
    SetPDGAxis(fTagEffLSWoWeight->GetAxis(2), PDGLabel);
    SetPDGAxis(fTagEffLSWoWeight->GetAxis(3), PDGLabel);
    //    SetPDGAxis(fTagEffLSWoWeight->GetAxis(4), PDGLabel);
    fOutputListMain->Add(fTagEffLSWoWeight);
    */


    fTagEffULSWoWeight = new THnSparseF("fTagEffULSWoWeight", "ULS-LS tag.eff. woW;  Pt; IsElectron; PDGMother; PDGGMother; PDGGGMIsNonHF", 5, binTagEff, xminTagEff, xmaxTagEff);
    //    SetPDGAxis(fTagEffULSWoWeight->GetAxis(1), PDGLabel);
    SetPDGAxis(fTagEffULSWoWeight->GetAxis(2), PDGLabel);
    SetPDGAxis(fTagEffULSWoWeight->GetAxis(3), PDGLabel);
    //    SetPDGAxis(fTagEffULSWoWeight->GetAxis(4), PDGLabel);
    fOutputListMain->Add(fTagEffULSWoWeight);


    fTagTruePairsWoWeight = new THnSparseF("fTagTruePairsWoWeight", "ULS true pairs tag.eff. woW:  Pt; IsElectron; PDGMother; PDGGMother; PDGGGMIsNonHF", 5, binTagEff, xminTagEff, xmaxTagEff);
    //    SetPDGAxis(fTagTruePairsWoWeight->GetAxis(1), PDGLabel);
    SetPDGAxis(fTagTruePairsWoWeight->GetAxis(2), PDGLabel);
    SetPDGAxis(fTagTruePairsWoWeight->GetAxis(3), PDGLabel);
    //    SetPDGAxis(fTagTruePairsWoWeight->GetAxis(4), PDGLabel);
    fOutputListMain->Add(fTagTruePairsWoWeight);

  }
  
  // HFE-HFE Correlation

  fElecTrigger = new TH1F("fEleTrigger","fEleTrigger", NBinsElectron, XminElectron, XmaxElectron);
  fOutputListMain->Add(fElecTrigger);
    
  fInclElecPhi = new TH2F("fInclElePhi", "fInclElePhi",NBinsElectron, XminElectron, XmaxElectron,NBinsPhi,0,TMath::TwoPi());
  fOutputListMain->Add(fInclElecPhi);
  
  fInclElecEta = new TH2F("fInclEleEta", "fInclEleEta",NBinsElectron, XminElectron, XmaxElectron,NBinsEta,XminEta,XmaxEta);
  fOutputListMain->Add(fInclElecEta);

  fInclElecPhiEta = new TH2F("fInclElePhiEta", "fInclElePhiEta", NBinsEta,XminEta,XmaxEta, NBinsPhi,0,TMath::TwoPi());
  fOutputListMain->Add(fInclElecPhiEta);
  
  fULSElecPhi= new TH2F("fULSElePhi", "fULSElePhi", NBinsElectron, XminElectron, XmaxElectron ,NBinsPhi,0,TMath::TwoPi());
  fOutputListMain->Add(fULSElecPhi);
  
  fLSElecPhi= new TH2F("fLSElePhi", "fLSElePhi",NBinsElectron, XminElectron, XmaxElectron,NBinsPhi,0,TMath::TwoPi());
  fOutputListMain->Add(fLSElecPhi);
  
  /*
  fElecDphi = new TH2F("fEleDphi", "fEleDphi",NBinsElectron, XminElectron, XmaxElectron,NBinsDPhi,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fOutputListMain->Add(fElecDphi);
  
  fULSElecDphi = new TH2F("fULSEleDphi", "fULSEleDphi",NBinsElectron, XminElectron, XmaxElectron,NBinsDPhi,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fOutputListMain->Add(fULSElecDphi);
  
  fLSElecDphi = new TH2F("fLSEleDphi", "fLSEleDphi",NBinsElectron, XminElectron, XmaxElectron,NBinsDPhi,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fOutputListMain->Add(fLSElecDphi);
  
  fULSElecDphiDiffMethod = new TH2F("fULSEleDphiDiffMethod", "fULSEleDphiDiffMethod",NBinsElectron, XminElectron, XmaxElectron,NBinsDPhi,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fOutputListMain->Add(fULSElecDphiDiffMethod);
  
  fLSElecDphiDiffMethod = new TH2F("fLSEleDphiDiffMethod", "fLSEleDphiDiffMethod",NBinsElectron, XminElectron, XmaxElectron,NBinsDPhi,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fOutputListMain->Add(fLSElecDphiDiffMethod);
  */

  Int_t     bin[5] = {NBinsHadRed   ,NBinsElectron, NBinsDPhi, NBinsDEta, NVertexBins}; //ptH, ptE, Dphi, Deta
  Double_t  xmin[5] = {XminHadron   ,XminElectron ,-TMath::Pi()/2    ,DEtaMin, -10};
  Double_t  xmax[5] = {XmaxHadron   ,XmaxElectron ,(3*TMath::Pi())/2 ,DEtaMax, 10};

  if (fCorrHadron) {

    fNoPartnerNoT = new TH1F("fNoPartnerNoT", "fNoParnterNoT", bin[1], xmin[1], xmax[1]);
    fOutputListHadron->Add(fNoPartnerNoT);

    Int_t TriggerBin[4]={5, 30, 30, fAssPtHad_Nbins};
    Double_t TrigMinBin[4]={0.5, 0., 0., -0.5};
    Double_t TrigMaxBin[4]={3., 1.5, 1.5, fAssPtHad_Nbins-0.5};
    fNoPartnerNoTPt2 = new THnSparseF("fNoPartnerNoTPt2", "fNoParnterNoTPt2; Pt1; Pt2Rec; Pt2MCm; AssHadBin", 4, TriggerBin, TrigMinBin, TrigMaxBin);
    fOutputListHadron->Add(fNoPartnerNoTPt2);
  
    fTPartnerNoT = new TH1F("fTPartnerNoT", "fTPartnerNoT", bin[1], xmin[1], xmax[1]);
    fOutputListHadron->Add(fTPartnerNoT);

    fTPartnerNoTPt2 =new THnSparseF("fTPartnerNoTPt2", "fNoParnterNoTPt2; Pt1; Pt2Rec; Pt2MCm; AssHadBin", 4, TriggerBin, TrigMinBin, TrigMaxBin);
    fOutputListHadron->Add(fTPartnerNoTPt2);

    fElecHadTrigger = new TH3F("fElecHadTrigger", "fElecHadTrigger; ptE; AssPtHad; EtaE", bin[1], xmin[1], xmax[1], fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5, 90, 0, 0.9);
    fOutputListHadron->Add(fElecHadTrigger);						      

    // fElecHadTriggerLS = new TH2F("fElecHadTriggerLS", "fElecHadTriggerLS", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    // fOutputListHadron->Add(fElecHadTriggerLS);						      

    fElecHadTriggerULS = new TH2F("fElecHadTriggerULS", "fElecHadTriggerULS-LS; ptE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListHadron->Add(fElecHadTriggerULS);		

    //    fElecHadTriggerLSNoP = new TH2F("fElecHadTriggerLSNoP", "fElecHadTriggerLSNoP", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    //    fOutputListHadron->Add(fElecHadTriggerLSNoP);						      

    fElecHadTriggerULSNoP = new TH2F("fElecHadTriggerULSNoP", "fElecHadTriggerULSNoP ULS-LS; ptE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListHadron->Add(fElecHadTriggerULSNoP);						      

    //    fElecHadTriggerLSNoPCorr = new TH2F("fElecHadTriggerLSNoPCorr", "fElecHadTriggerLSNoPCorr", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    //    fOutputListHadron->Add(fElecHadTriggerLSNoPCorr);						      

    fElecHadTriggerULSNoPCorr = new TH2F("fElecHadTriggerULSNoPCorr", "fElecHadTriggerULSNoPCorr ULS-LS; ptE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListHadron->Add(fElecHadTriggerULSNoPCorr);				      		      

    fElecHadTriggerULSNoPCorrTrue = new TH2F("fElecHadTriggerULSNoPCorrTrue", "fElecHadTriggerULSNoPCorrTrue ULS-LS; ptE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListHadron->Add(fElecHadTriggerULSNoPCorrTrue);	

    fHadContTrigger = new TH2F("fHadContTrigger", "fHadContTrigger; PtE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListHadron->Add(fHadContTrigger);						      

    fNonElecHadTrigger = new TH2F("fNonElecHadTrigger", "fNonElecHadTrigger; PtE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListHadron->Add(fNonElecHadTrigger);						      

    fMCElecHaTruePartnerTrigger = new TH2F("fMCElecHaTruePartnerTrigger", "fMCElecHaTruePartnerTrigger; PtE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListHadron->Add(fMCElecHaTruePartnerTrigger);        

    fMCElecHaTruePartnerTriggerWW = new TH2F("fMCElecHaTruePartnerTriggerWW", "fMCElecHaTruePartnerTriggerWW; PtE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListHadron->Add(fMCElecHaTruePartnerTriggerWW);      	

    fMCElecHaNoPartnerTrigger = new TH2F("fMCElecHaNoPartnerTrigger", "fMCElecHaNoPartnerTrigger; PtE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListHadron->Add(fMCElecHaNoPartnerTrigger);

    fMCElecHaNoPartnerTriggerWW = new TH2F("fMCElecHaNoPartnerTriggerWW", "fMCElecHaNoPartnerTriggerWW; PtE; AssPtHad", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
    fOutputListHadron->Add(fMCElecHaNoPartnerTriggerWW);

    fHadElecTrigger = new TH2F("fHadElecTrigger", "fHadElecTrigger; PtH; AssPtE", bin[0], xmin[0], xmax[0],  fAssPtElec_Nbins, -0.5, fAssPtElec_Nbins-0.5);
    fHadElecTrigger->GetXaxis()->Set(NBinsHadRed, XBinsHadRed);
    fOutputListHadron->Add(fHadElecTrigger);						      

    fHadNonElecTrigger = new TH2F("fHadNonElecTrigger", "fHadNonElecTrigger; PtH; AssPtE", bin[0], xmin[0], xmax[0],  fAssPtElec_Nbins, -0.5, fAssPtElec_Nbins-0.5);
    fHadNonElecTrigger->GetXaxis()->Set(NBinsHadRed, XBinsHadRed);
    fOutputListHadron->Add(fHadNonElecTrigger);	

    fInclElecHa = new THnSparseF("fEleHaIncl", "Sparse for Ele-Had : PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fInclElecHa->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fInclElecHa->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fInclElecHa);

    //    fLSElecHa = new THnSparseF("fEleHaLS", "Sparse for LSEle-Had : PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    // fLSElecHa->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    // fLSElecHa->GetAxis(4)->Set(NVertexBins, XVertexBins);
    //fOutputListHadron->Add(fLSElecHa);

    fULSElecHa = new THnSparseF("fEleHaULS", "Sparse for ULS-LS Ele-Had; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fULSElecHa->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fULSElecHa->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fULSElecHa);
 
    if (fIsMC) {
      
      fSignalElecHa = new THnSparseF("fSignalEleHa", "Sparse for ULSEle-Had : PtH, PtE, Dphi, Deta", 5, bin, xmin, xmax);
      fSignalElecHa->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fSignalElecHa->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListHadron->Add(fSignalElecHa);
     
      /*
      fBackgroundElecHa = new THnSparseF("fBackgroundEleHa", "Sparse for ULSEle-Had : PtH, PtE, Dphi, Deta", 5, bin, xmin, xmax);
      fBackgroundElecHa->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fBackgroundElecHa->GetAxis(4)->Set(NVertexBins, XVertexBins);
      //      fOutputListHadron->Add(fBackgroundElecHa);
      */ 
      fULSElecHaTrue = new THnSparseF("fEleHaULSTrue", "Sparse for ULS-LS True Ele-Had; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
      fULSElecHaTrue->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fULSElecHaTrue->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListHadron->Add(fULSElecHaTrue);
 
      fMCElecHaHadron = new THnSparseF("fMCEleHaHadron", "Sparse for Ele-Had; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
      fMCElecHaHadron->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fMCElecHaHadron->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListHadron->Add(fMCElecHaHadron);
      
    }

    //    fElecHaLSNoPartner = new THnSparseF("fEleHaLSNoPartner", "Sparse for LSEle-Had with no Partner: PtH, PtE, Dphi, Deta", 5, bin, xmin, xmax);
    // fElecHaLSNoPartner->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    //fElecHaLSNoPartner->GetAxis(4)->Set(NVertexBins, XVertexBins);
    //fOutputListHadron->Add(fElecHaLSNoPartner);

    fElecHaULSNoPartner = new THnSparseF("fEleHaULSNoPartner", "Sparse for ULS-LS Ele-Had with no Partner; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fElecHaULSNoPartner->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fElecHaULSNoPartner->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fElecHaULSNoPartner);

    
    //    fElecHaLSNoPartnerCorr = new THnSparseF("fEleHaLSNoPartnerCorr", "Sparse for LSEle-Had with no Partner: PtH, PtE, Dphi, Deta", 5, bin, xmin, xmax);
    // fElecHaLSNoPartnerCorr->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    //fElecHaLSNoPartnerCorr->GetAxis(4)->Set(NVertexBins, XVertexBins);
    //fOutputListHadron->Add(fElecHaLSNoPartnerCorr);

    fElecHaULSNoPartnerCorr = new THnSparseF("fEleHaULSNoPartnerCorr", "Sparse for ULS-LS Ele-Had with no Partner; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fElecHaULSNoPartnerCorr->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fElecHaULSNoPartnerCorr->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fElecHaULSNoPartnerCorr);

    /*
    fElecHaLSNoPartnerCorrTrue = new THnSparseF("fEleHaLSNoPartnerCorrTrue", "Sparse for LSEle-Had with no Partner: PtH, PtE, Dphi, Deta", 5, bin, xmin, xmax);
    fElecHaLSNoPartnerCorrTrue->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fElecHaLSNoPartnerCorrTrue->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fElecHaLSNoPartnerCorrTrue);
    */

    fElecHaULSNoPartnerCorrTrue = new THnSparseF("fEleHaULSNoPartnerCorrTrue", "Sparse for ULSEle-Had with no Partner; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fElecHaULSNoPartnerCorrTrue->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fElecHaULSNoPartnerCorrTrue->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fElecHaULSNoPartnerCorrTrue);
    

    if (fIsMC) {
      fMCElecHaTruePartner = new THnSparseF("fMCElecHaTruePartner", "Sparse for MC true photonics with no Partner; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
      fMCElecHaTruePartner->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fMCElecHaTruePartner->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListHadron->Add(fMCElecHaTruePartner);
    
      fMCElecHaNoPartner = new THnSparseF("fMCElecHaNoPartner", "Sparse for MC true photonics with no Partner: PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
      fMCElecHaNoPartner->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fMCElecHaNoPartner->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListHadron->Add(fMCElecHaNoPartner);	
    }

    fElecHaHa = new THnSparseF("fEleHaHa", "Sparse for Hadron Contamination; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fElecHaHa->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fElecHaHa->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fElecHaHa);

    fElecHaMixedEvent = new THnSparseF("fEleHaMixedEv", "Sparse for Ele-Had MixEvent; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fElecHaMixedEvent->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fElecHaMixedEvent->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fElecHaMixedEvent);
  
    //    fLSElecHaMixedEvent = new THnSparseF("fEleHaLSMixedEv", "Sparse for LSEle-Had MixEvent; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    //fLSElecHaMixedEvent->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    //fLSElecHaMixedEvent->GetAxis(4)->Set(NVertexBins, XVertexBins);
    //fOutputListHadron->Add(fLSElecHaMixedEvent);
  
    fULSElecHaMixedEvent = new THnSparseF("fEleHaULSMixedEv", "Sparse for ULS-LS Ele-Had MixEvent; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
    fULSElecHaMixedEvent->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fULSElecHaMixedEvent->GetAxis(4)->Set(NVertexBins, XVertexBins);
    fOutputListHadron->Add(fULSElecHaMixedEvent);

    if (fIsMC) {
      fTagHaMixedEvent = new THnSparseF("fTagHaMixedEv", "Sparse for Tag Ele-Had MixEvent; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
      fTagHaMixedEvent->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fTagHaMixedEvent->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListHadron->Add(fTagHaMixedEvent);
      
      fNonTagHaMixedEvent = new THnSparseF("fNonTagHaMixedEv", "Sparse for NonTag Ele-Had MixEvent; PtH; PtE; Dphi; Deta", 5, bin, xmin, xmax);
      fNonTagHaMixedEvent->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fNonTagHaMixedEvent->GetAxis(4)->Set(NVertexBins, XVertexBins);
      fOutputListHadron->Add(fNonTagHaMixedEvent);
    }
  }
  // LP HFE   

  if (fCorrLParticle) {

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


  if (fIsMC) {
    fCheckMCVertex = new TH2F("fCheckMCVertex", "TruthVsRec Vertex", 200, -10, 10, 200, -10, 10);
    fOutputListQA->Add(fCheckMCVertex);
  }



  //////////// Hadron and Electron RecEff //////////////////////

  Int_t    EffHBins[4]={NBinsHadRed, 18, 16, NVertexBins};
  Double_t EffHXmin[4]={XminHadron, -0.9, 0, -10};
  Double_t EffHXmax[4]={XmaxHadron, 0.9, TMath::TwoPi(), 10};

  Int_t    EffEBins[4]={NBinsElectronRed, 18, 16, NVertexBins};
  Double_t EffEXmin[4]={XminElectron, -0.9, 0, -10};
  Double_t EffEXmax[4]={100, 0.9, TMath::TwoPi(), 10};
    
  if (!fOneTimeCheck) { // maybe reduce axis further
    EffHBins[2]=1; // hadrons for pt, eta, zVtx
    EffEBins[1]=1; // electrons only pt, zVtx
    EffEBins[2]=1;
  }

  fRecHFE = new TH1F("fRecHFE", "fRecHFE; #it{p}_{T}", NBinsElectronRed, XminElectron, XmaxElectron);
  fOutputListMain->Add(fRecHFE);

  fRecHFEEtaWRecEff = new TH2F("fRecHFEEtaWRecEff", "fRecHFEEtaWRecEff; Pt; Eta", NBinsElectron , XminElectron, XmaxElectron, 250, -2.5, 2.5); 
  fOutputList->Add(fRecHFEEtaWRecEff);


  // each time count corrected hadrons and electrons (check for right efficiency)
  fRecHadPtEtaPhiVtxWRecEff=new THnSparseF("fRecHadPtEtaPhiVtxWRecEff", "Rec hadrons wRecEff; rec #it{p}_{T}; rec #eta; rec #varphi; rec zVtx", 4, EffHBins, EffHXmin, EffHXmax);
  fRecHadPtEtaPhiVtxWRecEff->GetAxis(3)->Set(NVertexBins, XVertexBins);
  fRecHadPtEtaPhiVtxWRecEff->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
  fOutputListMain->Add(fRecHadPtEtaPhiVtxWRecEff);   

  fRecElecPtEtaPhiVtxWRecEff=new THnSparseF("fRecElePtEtaPhiWRecEff", "Rec electrons w. rec.Eff; rec #it{p}_{T}; rec #eta; rec #varphi; rec zVtx", 4, EffEBins, EffEXmin, EffEXmax);
  fRecElecPtEtaPhiVtxWRecEff->GetAxis(3)->Set(NVertexBins, XVertexBins);
  fRecElecPtEtaPhiVtxWRecEff->GetAxis(0)->Set(NBinsElectronRed, XBinsElectronRed);
  fOutputListMain->Add(fRecElecPtEtaPhiVtxWRecEff);

  if (fIsMC) {
    // accordingly check truth generated hadrons and electrons
    fMCHadPtEtaPhiVtx=new THnSparseF("fMCHadPtEtaPhiVtx", "MC truth gen. hadrons; MC #it{p}_{T}; MC #eta; MC #varphi; MC zVtx", 4, EffHBins, EffHXmin, EffHXmax);
    fMCHadPtEtaPhiVtx->GetAxis(3)->Set(NVertexBins, XVertexBins);
    fMCHadPtEtaPhiVtx->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fOutputListMain->Add(fMCHadPtEtaPhiVtx);
    
    fMCElecPtEtaPhiVtx=new THnSparseF("fMCElePtEtaPhiVtx", "MC truth gen. electrons;  MC #it{p}_{T}; MC #eta; MC #varphi; MC zVt", 4, EffEBins, EffEXmin, EffEXmax);
    fMCElecPtEtaPhiVtx->GetAxis(3)->Set(NVertexBins, XVertexBins);
    fMCElecPtEtaPhiVtx->GetAxis(0)->Set(NBinsElectronRed, XBinsElectronRed);
    fOutputListMain->Add(fMCElecPtEtaPhiVtx);
  }
  if (fRecEff) { // for rec eff count rec. hadrons and electrons wo. receff correction
    fRecHadPtEtaPhiVtx=new THnSparseF("fRecHadPtEtaPhiVtx", "Rec hadrons; rec #it{p}_{T}; rec #eta; rec #varphi; rec zVtx", 4, EffHBins, EffHXmin, EffHXmax);
    fRecHadPtEtaPhiVtx->GetAxis(3)->Set(NVertexBins, XVertexBins);
    fRecHadPtEtaPhiVtx->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
    fOutputListMain->Add(fRecHadPtEtaPhiVtx);
  
    fRecElecPtEtaPhiVtx=new THnSparseF("fRecElePtEtaPhi", "Rec electrons; rec #it{p}_{T}; rec #eta; rec #varphi; rec zVtx", 4, EffEBins, EffEXmin, EffEXmax);
    fRecElecPtEtaPhiVtx->GetAxis(3)->Set(NVertexBins, XVertexBins);
    fRecElecPtEtaPhiVtx->GetAxis(0)->Set(NBinsElectronRed, XBinsElectronRed);
    fOutputListMain->Add(fRecElecPtEtaPhiVtx);

    if (fIsMC & fOneTimeCheck) { // check mc pt .. for rec. hadrons and electrons wo receff correction
       
      fRecHadMCPtEtaPhiVtx=new THnSparseF("fRecHadMCPtEtaPhiVtx", "MC rec hadrons w. MC ; MC #it{p}_{T}; MC #eta; MC #varphi; MC zVtx", 4, EffHBins, EffHXmin, EffHXmax);
      fRecHadMCPtEtaPhiVtx->GetAxis(3)->Set(NVertexBins, XVertexBins);
      fRecHadMCPtEtaPhiVtx->GetAxis(0)->Set(NBinsHadRed, XBinsHadRed);
      fOutputListMain->Add(fRecHadMCPtEtaPhiVtx);    


      fRecElecMCPtEtaPhiVtx=new THnSparseF("fRecEleMCPtEtaPhi", "MC rec electrons w. MC; MC #it{p}_{T}; MC #eta; MC #varphi; MC zVtx", 4, EffEBins, EffEXmin, EffEXmax);
      fRecElecMCPtEtaPhiVtx->GetAxis(3)->Set(NVertexBins, XVertexBins);
      fRecElecMCPtEtaPhiVtx->GetAxis(0)->Set(NBinsElectronRed, XBinsElectronRed);
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
      // fMCElecPtEtaPhiStrictVtx->GetAxis(0)->Set(NBinsElectronRed, XBinsElectronRed);
      // fOutputListMain->Add(fMCElecPtEtaPhiStrictVtx);
      
    }
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

    fMCPiPlusProd = new THnSparseF("fMCPiPlusProd", "fMCPiPlusProd; Pt; Eta; Y; PDGMoth; Cocktail", 5, Pi0EtaBins, Pi0EtaXmin, Pi0EtaXmax);
    BinLogX(fMCPiPlusProd->GetAxis(0));
    SetPDGAxis(fMCPiPlusProd->GetAxis(3), PDGLabel);
    fOutputListMain->Add(fMCPiPlusProd);

    fMCPiPlusProdV2 = new THnSparseF("fMCPiPlusProdV2", "fMCPiPlusProdV2; Pt; Eta; Y; PDGMoth; Cocktail", 5, Pi0EtaBins, Pi0EtaXmin, Pi0EtaXmax);
    BinLogX(fMCPiPlusProdV2->GetAxis(0));
    SetPDGAxis(fMCPiPlusProdV2->GetAxis(3), PDGLabel);
    fOutputListMain->Add(fMCPiPlusProdV2);
    
   
  }

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
  }

  fCompareLPRecCheck = new TH3F("fCompareLPRecCheck", "fMCLPRecCheck; #Delta pt; #Delta phi; #Delta eta", 50, -1, 4, 36, -TMath::Pi()/2, (3.*TMath::Pi())/2, 40, -4, 4);
  fOutputListLP->Add(fCompareLPRecCheck);


  if (fMCTrueCorrelation) {
    Int_t    poolSize = 1000;
    Int_t    trackDepth = 2000; 
    const Int_t    nMaxPtBins=5;
    Double_t maxPtBins[]={0, 1000, 1001, 1002, 1003, 9999}; //{0.5, 2., 5., 10, 9999};

    //  8 Vertex bins * 3 MultBins, * 5 PtBins
    fMCTruePoolMgr = new AliEventPoolManager(poolSize, trackDepth, NMultBins, (Double_t*) XMultBins, NVertexBins, (Double_t*) XVertexBins, nMaxPtBins,(Double_t*)maxPtBins);
    fMCTruePoolMgr->Validate();  


    Int_t MCTrueBins[5]= {NBinsHadRed   ,NBinsElectron, NBinsDPhi, NBinsDEta, 10}; 
    Double_t MCTrueXmin[5] = {XminHadron   ,XminElectron ,-TMath::Pi()/2    ,DEtaMin, 0.5};
    Double_t MCTrueXmax[5] ={XmaxHadron   ,XmaxElectron ,(3*TMath::Pi())/2 ,DEtaMax, 10.5};

    Int_t MCTrueMEvBins[5]= {NBinsHadRed   ,NBinsElectron, NBinsDPhi, NBinsDEta, NVertexBins}; 
    Double_t MCTrueMEvXmin[5] = {XminHadron   ,XminElectron ,-TMath::Pi()/2    ,DEtaMin,-10};
    Double_t MCTrueMEvXmax[5] ={XmaxHadron   ,XmaxElectron ,(3*TMath::Pi())/2 ,DEtaMax, 10};

    
    fCompareLP = new TH3F("fCompareLP", "fComparelP; RecLP; LPinAcc; LP", NBinsHadRed, XBinsHadRed, NBinsHadRed, XBinsHadRed,NBinsHadRed, XBinsHadRed);
    fOutputList->Add(fCompareLP);

   
    fRecHadronEtaWRecEff = new TH2F("fRecHadronEtaWRecEff", "fRecHadronEtaWRecEff; Pt; Eta", NBinsHadron , XminHadron, XmaxHadron, 250, -2.5, 2.5); 
    fOutputList->Add(fRecHadronEtaWRecEff);

    fRecLPEta = new TH2F("fRecLPEta", "fRecLPEta; Pt; Eta", NBinsHadRed , XminHadron, XmaxHadron, 250, -2.5, 2.5); 
    fRecLPEta->GetXaxis()->Set(NBinsHadRed, XBinsHadRed);
    fOutputList->Add(fRecLPEta);


    if (fIsMC) {
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
	fTrueElectronEta = new TH2F("fTrueElectronEta", "fTrueElectronEta; Pt; Eta", NBinsElectron , XminElectron, XmaxElectron, 250, -2.5, 2.5); 
	fOutputList->Add(fTrueElectronEta);
	
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

  Int_t    poolSize = 1000;
  Int_t    trackDepth = 2000; 
  const Int_t    nMaxPtBins=5;
  Double_t maxPtBins[]={0,1000, 1001, 1002, 1003, 9999}; //{ 0.5, 2., 5., 10, 9999};

  //  8 Vertex bins * 3 MultBins, * 5 PtBins
  fPoolMgr = new AliEventPoolManager(poolSize, trackDepth, NMultBins, (Double_t*) XMultBins, NVertexBins, (Double_t*) XVertexBins, nMaxPtBins,(Double_t*)maxPtBins);
  fPoolMgr->Validate();
    
  fPoolIsFilled = new TH3F("fPoolIsFilled", "fPoolIsFilled: mult, vertex, LP pt", NMultBins, (Double_t*) XMultBins, NVertexBins, (Double_t*) XVertexBins, nMaxPtBins,(Double_t*)maxPtBins);
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
AliVTrack*  AliAnalysisTaskHaHFECorrel::FindLPAndHFE( TObjArray* RedTracks, const AliVVertex *pVtx, Int_t nMother, Int_t listMother[], Double_t mult, Bool_t & EvContTP, Bool_t & EvContNTP)
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
    fTrkpt->Fill(pt,0);

    // track cuts
    Bool_t passHadTrackCut=kFALSE;
    passHadTrackCut = ChargedHadronTrackCuts(pVtx, Vtrack, nMother, listMother);
   
    Bool_t passHadPIDCut=kFALSE;
    passHadPIDCut = ChargedHadronPIDCuts(Vtrack); // currently empty
    
    if (passHadTrackCut && passHadPIDCut) { // calculate Hadron reconstruction efficiency
      recEffH = GetHadronRecEff(pt, phi, eta, pVtx->GetZ());
      if (recEffH>0) {
	NHadrons+=(1./recEffH);
	if (fRecEff) fRecHadronEtaWRecEff->Fill(pt, eta, 1./recEffH);
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


      fEtaVtxZ->Fill(eta, pVtx->GetZ());
     
      // fill rec hadron histograms for RecEff 
      Double_t fillSparse[4]={-999,-999,-999,-999};
      fillSparse[0]=pt;
      fillSparse[1]=eta;
      fillSparse[2]=phi;
      fillSparse[3]=pVtx->GetZ();
      if (recEffH>0) fRecHadPtEtaPhiVtxWRecEff->Fill(fillSparse, 1./recEffH);
      if (fRecEff && fIsMC)  fRecHadPtEtaPhiVtx->Fill(fillSparse);
      
      if (fRecEff && fOneTimeCheck) {
	if (fIsAOD && fIsMC) {
	  Int_t MClabel=AODtrack->GetLabel();
	  AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(MClabel));  
	  Double_t mcVtx[3];
	  fMCheader->GetVertex(mcVtx);
	  fillSparse[0]=MCParticle->Pt();
	  fillSparse[1]=MCParticle->Eta();
	  fillSparse[2]=MCParticle->Phi();
	  fillSparse[3]=mcVtx[2];
	  fRecHadMCPtEtaPhiVtx->Fill(fillSparse);

	  // controll plot for MC vs Rec pt;
	  fCheckMCPtvsRecPtHad->Fill(MCParticle->Pt(), pt);
	  fCheckMCEtavsRecEtaHad->Fill(MCParticle->Eta(), eta);
	  fCheckMCPhivsRecPhiHad->Fill(MCParticle->Phi(), phi);
	}
      }
    }  
 

    // Select HFE candidates
    Bool_t passHFETrackCut=kFALSE;
    passHFETrackCut= InclElecTrackCuts(pVtx,Vtrack,nMother,listMother);
    if (passHFETrackCut)  fTrkpt->Fill(pt,1); // after track cuts
 
    Bool_t passHFEPIDCut=kFALSE;
    if (passHFETrackCut) passHFEPIDCut= InclElecPIDCuts(Vtrack, kTRUE);
    if (passHFETrackCut && passHFEPIDCut) {
      recEffE = GetElectronRecEff(pt, phi, eta, pVtx->GetZ());
      fTrkpt->Fill(pt,2);
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
      FindPhotonicPartner(jTracks, Vtrack, pVtx, nMother, listMother, lsPartner, ulsPartner, lsPartnerID, ulsPartnerID, lsPartnerWeight, ulsPartnerWeight, trueULSPartner, isPhotonic, MCPartnerPt, RecPartnerPt);
      if (fIsMC) CheckPhotonicPartner(Vtrack, trueULSPartner, MCPartnerPt, RecPartnerPt);
      if (fIsMC && fTagEff) {
	if (isPhotonic) {
	  if (trueULSPartner) {
	    NPhotElectronsTagged+=(1./recEffE);
	    fTagEtaPt1Pt2->Fill(Vtrack->Eta(), Vtrack->Pt(), RecPartnerPt);
	    // fTagEtaPhiPt->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt());
	    // fTagEtaZvtxPt->Fill(Vtrack->Eta(), pVtx->GetZ(), Vtrack->Pt());
	    // fTagEtaPhiPtwW->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt(), 1./recEffE);
	    // fTagEtaZvtxPtwW->Fill(Vtrack->Eta(), pVtx->GetZ(), Vtrack->Pt(), 1./recEffE);
	    if (Vtrack->Pt()>1.0) EvContTP=kTRUE;
	  }
	  else  {
	    NPhotElectronsUntagged+=(1./recEffE);
	    fNonTagEtaPt1Pt2->Fill(Vtrack->Eta(), Vtrack->Pt(), RecPartnerPt);
	    // fNonTagEtaPhiPt->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt());
	    // fNonTagEtaZvtxPt->Fill(Vtrack->Eta(), pVtx->GetZ(), Vtrack->Pt());
	    // fNonTagEtaPhiPtwW->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt(), 1./recEffE);
	    // fNonTagEtaZvtxPtwW->Fill(Vtrack->Eta(), pVtx->GetZ(), Vtrack->Pt(), 1./recEffE);
	    if (Vtrack->Pt()>1.0) EvContNTP=kTRUE;
	  }
	}
	EvaluateTaggingEfficiency(Vtrack, lsPartner, ulsPartner, trueULSPartner);
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
      fInclElecPtEta->Fill(pt, eta);
      fInclElecPtEtaWRecEff->Fill(pt, eta, 1./recEffE);
      //      fInclElecP->Fill(p);
      fInclElecPhi->Fill(pt,phi); 
      fInclElecEta->Fill(pt,eta); 
      fInclElecPhiEta->Fill(eta, phi);

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
	  Int_t MIsHeavy  = Int_t (mcMotherPDG / TMath::Power(10, Int_t(TMath::Log10(mcMotherPDG))));
	  if (MIsHeavy>3 && MIsHeavy<6) {
	    Double_t fillSparse[4];
	    fillSparse[0]=pt;
	    fillSparse[1]=eta;
	    fillSparse[2]=phi;
	    fillSparse[3]=pVtx->GetZ();
	    fRecElecPtEtaPhiVtxWRecEff->Fill(fillSparse, 1./recEffE); // for checks if right efficiency has been used
	    
	    fRecHFE->Fill(pt);
	    fRecHFEEtaWRecEff->Fill(pt, eta, 1./recEffE);
	    if (fRecEff) {
	      // Fill Reconstructed Histogram
	      fRecElecPtEtaPhiVtx->Fill(fillSparse);                 
	      fillSparse[0]=MCParticle->Pt();
	      fillSparse[1]=MCParticle->Eta();
	      fillSparse[2]=MCParticle->Phi();
	      fillSparse[3]=mcVtx[2];
	      if (fOneTimeCheck) {
		fRecElecMCPtEtaPhiVtx->Fill(fillSparse);
		fCheckMCPtvsRecPtEle->Fill(MCParticle->Pt(), pt);
	      }
	    }
	  }
	}

	if (PDGCode != 11) {
	  NHadCont+=(1./recEffE);
	  isHadron=kTRUE;
	  if (fHadCont) fHadContEtaPhiPt->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt());
	}
      }
        
      // store only essential informations of these electrons for later correlations and mixed event
      CloneAndReduceTrackList(RedTracks, Vtrack, lsPartner, ulsPartner, lsPartnerID, ulsPartnerID,  lsPartnerWeight,  ulsPartnerWeight, trueULSPartner, MCPartnerPt, RecPartnerPt, isPhotonic, isHadron);
    }

    Bool_t passNonElecPIDCut=kFALSE;
    passNonElecPIDCut=AssoHadronPIDCuts(Vtrack);
    if (passHFETrackCut && passNonElecPIDCut && recEffE>0) {
      NNonElectrons+=(1./recEffE);
      fHadContTPCEtaPhiPt->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt());
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
    fMCLeadingParticle->Fill(fillSparse);

    if (LPRecChecktrack) {
      AliAODTrack *LPRecCheckTrack = dynamic_cast<AliAODTrack*>(LPRecChecktrack);   
      Double_t LPDeltaPhi = GetDeltaPhi(LPAODtrack->Phi(), LPRecCheckTrack->Phi());
      Double_t LPDeltaEta = GetDeltaEta(LPAODtrack->Eta(), LPRecCheckTrack->Eta());
      Double_t LPDeltaPt  = LPAODtrack->Pt()  -LPRecCheckTrack->Pt();
      fCompareLPRecCheck->Fill(LPDeltaPt, LPDeltaPhi, LPDeltaEta);
    }




  }
  //  if (fLParticle && RedTracks->GetEntriesFast()>0) cout<< "RecLP " << LPtrack->GetLabel() <<  "\t" << LPtrack->Pt() << endl;

  return LPtrack;
}

//_________________________________________
void AliAnalysisTaskHaHFECorrel::FindPhotonicPartner(Int_t iTracks, AliVTrack* Vtrack, const AliVVertex *pVtx, Int_t nMother, Int_t listMother[], Int_t &lsPartner, Int_t &ulsPartner, Int_t *lsPartnerID, Int_t *ulsPartnerID,  Float_t *lsPartnerWeight, Float_t *ulsPartnerWeight, Bool_t &trueULSPartner, Bool_t &isPhotonic, Float_t &MCPartnerPt, Float_t &RecPartnerPt) {
  

  AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);
  if(!AODtrack) return;

  lsPartner=0;
  ulsPartner=0;
  trueULSPartner=kFALSE; // only for MC use
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();
    
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
  
    // tagged particle
    Double_t pt=-9.,phi=-9.,eta =-9.,recEff=-9.;
    Int_t charge = 0;      
    pt = Vtrack->Pt();
    phi = Vtrack->Phi();
    eta = Vtrack->Eta();
    charge = Vtrack->Charge();
    recEff = GetElectronRecEff(pt,phi,eta,pVtx->GetZ());
    if (recEff<0) continue;
        
    // associated particle variables
    Double_t phiAsso=-9.;
    Int_t chargeAsso = 0;
    phiAsso = VtrackAsso->Phi();
    chargeAsso = VtrackAsso->Charge();

    double dphi = -9;

    // Cuts for associated electrons for HFE-HFE correlations
    Bool_t passAssoTrackCutIncl = kFALSE;
    passAssoTrackCutIncl = InclElecTrackCuts(pVtx,VtrackAsso,nMother,listMother);

    Bool_t passAssoPIDCutIncl = kFALSE;
    passAssoPIDCutIncl = InclElecPIDCuts(VtrackAsso, kFALSE);

    if (passAssoTrackCutIncl && passAssoPIDCutIncl) {
      dphi = GetDeltaPhi(phi,phiAsso); // electron-electron dphi
      //  fElecDphi->Fill(pt,dphi);        // HFEHFE
    }

    // looser Track cuts for associated photonic eletron
    Bool_t passAssoTrackCutPhot = kFALSE;
    passAssoTrackCutPhot = PhotElecTrackCuts(pVtx,VtrackAsso,nMother,listMother);
    if(!passAssoTrackCutPhot) continue;
   
    Bool_t passAssoPIDCutPhot = kFALSE;
    passAssoPIDCutPhot = PhotElecPIDCuts(VtrackAsso); 
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
      fOpeningAngleLS->Fill(openingAngle,pt);
      fInvmassLS->Fill(mass,pt);
    }
    if(fFlagULS){
      fOpeningAngleULS->Fill(openingAngle,pt);
      fInvmassULS->Fill(mass,pt);
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
	    fInvmassMCTrue->Fill(mass,pt,PDGMap.find(PdgCode)->second);
	    fInvmassMCTrue->GetZaxis()->SetBinLabel(PDGMap.find(PdgCode)->second, Form("%i", PdgCode));
	    RecPartnerPt=AODtrackAsso->Pt();
	    MCPartnerPt=MCParticleAsso->Pt();
	    if (IsPhotonicElectron(AODtrack->GetLabel())&& mass<fInvmassCut){
	      fPhotPt2MCRec->Fill(AODtrackAsso->Pt(), MCParticleAsso->Pt());
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



      fULSElecPhi->Fill(pt,phi);
      fULSElecPt->Fill(pt);
      fULSElecPtWRecEff->Fill(pt, 1./recEff);
      dphi = GetDeltaPhi(phi,phiAsso); // electron-electron dphi
      // fULSElecDphi->Fill(pt,dphi);
      ulsPartnerID[ulsPartner]=VtrackAsso->GetID();
      ulsPartnerWeight[ulsPartner]=GetNonTagCorr(pt, VtrackAsso->Pt());
      ulsPartner++;
      // if (fIsAOD && fIsMC && trueULSPartner) recTruePairMass=mass;
      if (fIsAOD)   fPhotPt1Pt2Only->Fill(AODtrack->Pt(), AODtrackAsso->Pt(), 1.);
    }
        
    if(mass<fInvmassCut && fFlagLS){
      fLSElecPhi->Fill(pt,phi);
      fLSElecPt->Fill(pt);
      //  cout << recEff << endl;
      fLSElecPtWRecEff->Fill(pt, 1./recEff);
      dphi = GetDeltaPhi(phi,phiAsso); // electron-electron dphi
      // fLSElecDphi->Fill(pt,dphi);
      lsPartnerID[lsPartner]=VtrackAsso->GetID();
      lsPartnerWeight[lsPartner]=GetNonTagCorr(pt, VtrackAsso->Pt());
      lsPartner++;
      if (fIsAOD)   fPhotPt1Pt2Only->Fill(AODtrack->Pt(), AODtrackAsso->Pt(), -1.);
    }
  
  }//track loop
  if (fIsMC) { 
    isPhotonic=IsPhotonicElectron(AODtrack->GetLabel());
  } 
}



void AliAnalysisTaskHaHFECorrel::CheckPhotonicPartner(AliVTrack* Vtrack, Bool_t Tagged, Float_t& MCPartnerPt, Float_t RecPartnerPt) {


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
    if (foundPartner) cout << "should never happen " << endl; // but in p0->eeee
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
  if (factor>0) fPhotPt1PtMTag->Fill(MCParticle->Pt(), MCMother->Pt());
  else fPhotPt1PtMNTag->Fill(MCParticle->Pt(), MCMother->Pt());

  fPhotPt1Pt2->Fill(fillSparse, PtMotherWeight);
  if (factor>0) fPhotPt1Pt2Corr->Fill(fillSparse, PtMotherWeight*GetNonTagCorr(AODtrack->Pt(), RecPartnerPt)); 
  else  fPhotPt1Pt2Corr->Fill(fillSparse, PtMotherWeight);
  fPhotPt1Pt2MC->Fill(fillSparse);

  fillSparse[1]=factor*RecPartnerPt;
  if (factor<0) fillSparse[1]=factor*ptPartner;
  //fPhotPt1Pt2Rec->Fill(fillSparse, PtMotherWeight);
 
  // only to check the invariant mass distribution
  fillSparse[1]=mass*factor;
  /*
  fPhotPt1E->Fill(fillSparse, PtMotherWeight);
  if (factor>0)   fPhotPt1ECorr->Fill(fillSparse, PtMotherWeight*GetNonTagCorr(AODtrack->Pt(), RecPartnerPt));
  else   fPhotPt1ECorr->Fill(fillSparse, PtMotherWeight);
  */
  fillSparse[2]=ptPartner;
  fPhotPt1Pt2E->Fill(fillSparse, PtMotherWeight);
  if (factor>0)   fPhotPt1Pt2ECorr->Fill(fillSparse, PtMotherWeight*GetNonTagCorr(AODtrack->Pt(), RecPartnerPt));
  else   fPhotPt1Pt2ECorr->Fill(fillSparse, PtMotherWeight);
  


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
void AliAnalysisTaskHaHFECorrel::CorrelateHadron(TObjArray* RedTracksHFE,  const AliVVertex* pVtx, Int_t nMother, Int_t listMother[], Float_t mult, Float_t maxPt) {
    
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
    passTrackCut = ChargedHadronTrackCuts(pVtx, track, nMother, listMother);
    if (!passTrackCut) continue;

    // pid cuts
    Bool_t passPIDCut=kFALSE;
    passPIDCut = ChargedHadronPIDCuts(track);
    if (!passPIDCut) continue;
   
    // associated particle
    Double_t ptH=-9.,etaH =-9.,phiH=-9.,recEffH=-9.;
    Int_t idH=-9;
    ptH = track->Pt();
    phiH = track->Phi();
    etaH = track->Eta();
    idH = track->GetID();
    recEffH = GetHadronRecEff(ptH,phiH,etaH,pVtx->GetZ());
    if (recEffH<0) continue;
 
    // Correlate hadrons with Hadron (kTRUE, kFALSE);
    CorrelateWithHadrons(track, pVtx, nMother, listMother, kTRUE, kFALSE, NonElectronIsTrigger, NonElectronIsTriggerPt, NonElectronIsTriggerWeight, RedTracksHFE->GetEntriesFast()); 

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
      recEffE = GetElectronRecEff(pt,phi,eta,pVtx->GetZ());
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
      fInclElecHa->Fill(fillSparse, 1./(recEffH*recEffE)); 
      //   fLSElecHa->Fill(fillSparse, ls/(recEffH*recEffE)); // not needed anymore
      fULSElecHa->Fill(fillSparse, uls*1./(recEffH*recEffE)); // TriggerConditionProvedAtnd
      fULSElecHa->Fill(fillSparse, -ls*1./(recEffH*recEffE)); // TriggerConditionProvedAtnd
      if (RedTrack->IsPhotonic()) {
	fULSElecHaTrue->Fill(fillSparse, uls*1./(recEffH*recEffE));
	fULSElecHaTrue->Fill(fillSparse, -ls*1./(recEffH*recEffE));
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
	  fElecHaULSNoPartner->Fill(fillSparse, uls*1./(recEffH*recEffE));
	  fElecHaULSNoPartner->Fill(fillSparse, -ls*1./(recEffH*recEffE));
	  //fElecHaLSNoPartner->Fill(fillSparse, ls/(recEffH*recEffE));

	  for (int j=0; j<ls; j++) {
	    fElecHaULSNoPartnerCorr->Fill(fillSparse, -RedTrack->LSPartnerWeight(j)/(recEffH*recEffE));
	  }
	  for (int j=0; j<uls; j++) {
	    fElecHaULSNoPartnerCorr->Fill(fillSparse, RedTrack->ULSPartnerWeight(j)/(recEffH*recEffE));
	  }

	  if (RedTrack->IsPhotonic()) {
	    CheckElectronIsTrigger(ptH, ElectronIsTriggerNoPTrue[k]);
	    fElecHaULSNoPartnerCorrTrue->Fill(fillSparse, uls*1./(recEffH*recEffE));
	    fElecHaULSNoPartnerCorrTrue->Fill(fillSparse, -ls*1./(recEffH*recEffE));
	  }
	}
      }
   
      

      if (fIsMC) {
	if (RedTrack->IsHadron()){
	  CheckElectronIsTrigger(ptH, HadContIsTrigger[k]);
	  fMCElecHaHadron->Fill(fillSparse, 1./(recEffH*recEffE));
	}
	if  (RedTrack->IsPhotonic()) {
	  if (RedTrack->TruePartner()) {
	    CheckElectronIsTrigger(ptH, PhotElecWPartnerTrigger[k]);
	    fMCElecHaTruePartner->Fill(fillSparse, 1./(recEffH*recEffE));
	  }
	  else {
	    CheckElectronIsTrigger(ptH,  PhotElecWoPartnerTrigger[k]);
	    fMCElecHaNoPartner->Fill(fillSparse, 1./(recEffH*recEffE));
	  }
	}
	if (fIsAOD) {
	  AliAODMCParticle* MCParticleAOD = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(abs(RedTrack->GetLabel())));  
	  if (!MCParticleAOD) continue;
	  Int_t PDGCode = abs(MCParticleAOD->GetPdgCode());
	  if (PDGCode!=11) {
	    // fBackgroundElecHa->Fill(fillSparse, 1./(recEffH*recEffE));
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
	    if (MotherPDG==11) {
	      Int_t IsHeavyGM =  Int_t (GMotherPDG / TMath::Power(10, Int_t(TMath::Log10(GMotherPDG))));
	      if (IsHeavyGM>3 && IsHeavyGM<6)  fSignalElecHa->Fill(fillSparse, 1./(recEffH*recEffE));
	      //else fBackgroundElecHa->Fill(fillSparse, 1./(recEffH*recEffE));
	    }
	    else {
	      Int_t  IsHeavy =  Int_t (MotherPDG / TMath::Power(10, Int_t(TMath::Log10(MotherPDG))));
	      if (IsHeavy>3 && IsHeavy<6)  fSignalElecHa->Fill(fillSparse, 1./(recEffH*recEffE));
	      //else  fBackgroundElecHa->Fill(fillSparse, 1./(recEffH*recEffE));
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
      if (HadronTrigger[AssPtBin]) fHadElecTrigger->Fill(ptH, AssPtBin, 1./recEffH);
    }

  } //end of track loop

  // Fill Trigger Histograms
  for (Int_t l=0; l<RedTracksHFE->GetEntriesFast(); l++) {
    AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(l);
    Double_t recEff = GetElectronRecEff(RedTrack->Pt(), RedTrack->Phi(), RedTrack->Eta(), pVtx->GetZ());
    if (recEff<0) continue;
    for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
      if (ElectronIsTrigger[l][AssPtBin]) {
	fElecHadTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin, RedTrack->Eta(), 1./recEff);
	if( RedTrack->LS()>0) fElecHadTriggerULS->Fill(ElectronIsTriggerPt[l], AssPtBin,  -RedTrack->LS()/recEff);
	if( RedTrack->ULS()>0) fElecHadTriggerULS->Fill(ElectronIsTriggerPt[l], AssPtBin,  RedTrack->ULS()/recEff);
      }
      if (ElectronIsTriggerNoP[l][AssPtBin]) {
	if( RedTrack->LS()>0) fElecHadTriggerULSNoP->Fill(ElectronIsTriggerPt[l], AssPtBin,   -RedTrack->LS()/recEff);
	if( RedTrack->ULS()>0) fElecHadTriggerULSNoP->Fill(ElectronIsTriggerPt[l], AssPtBin,  RedTrack->ULS()/recEff);
     	
	for (int j=0; j<RedTrack->LS(); j++) fElecHadTriggerULSNoPCorr->Fill(ElectronIsTriggerPt[l], AssPtBin,   -RedTrack->LSPartnerWeight(j)/recEff);
	for (int j=0; j<RedTrack->ULS(); j++) fElecHadTriggerULSNoPCorr->Fill(ElectronIsTriggerPt[l], AssPtBin,   RedTrack->ULSPartnerWeight(j)/recEff);
      }
      if (ElectronIsTriggerNoPTrue[l][AssPtBin]) {
	if( RedTrack->LS()>0)  fElecHadTriggerULSNoPCorrTrue->Fill(ElectronIsTriggerPt[l], AssPtBin,  -RedTrack->LS()/recEff);
	if( RedTrack->ULS()>0) fElecHadTriggerULSNoPCorrTrue->Fill(ElectronIsTriggerPt[l], AssPtBin,  RedTrack->ULS()/recEff);
      }
      if (PhotElecWPartnerTrigger[l][AssPtBin]) {
	fMCElecHaTruePartnerTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin, 1./recEff);
	fMCElecHaTruePartnerTriggerWW->Fill(ElectronIsTriggerPt[l], AssPtBin, MotherWeight[l]*1./recEff);
      }
      if (PhotElecWoPartnerTrigger[l][AssPtBin]) {
	fMCElecHaNoPartnerTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin, 1./recEff);
	fMCElecHaNoPartnerTriggerWW->Fill(ElectronIsTriggerPt[l], AssPtBin, MotherWeight[l]*1./recEff);
      }
      if (HadContIsTrigger[l][AssPtBin]) fHadContTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin, 1./recEff);
    }
  }


  for (Int_t l=0; l<RedTracksHFE->GetEntriesFast(); l++) {
    AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(l);
    Double_t recEff = GetElectronRecEff(RedTrack->Pt(), RedTrack->Phi(), RedTrack->Eta(), pVtx->GetZ());
    if (recEff<0) continue;
    //  for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
    Bool_t NoPNoT = kTRUE;
    Bool_t TPNoT = kTRUE;
    if (RedTrack->IsPhotonic()) {
      for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
	Double_t fillSparse[4]={RedTrack->Pt(), RedTrack->TruePartnerRecPt(), RedTrack->TruePartnerMCPt(), 1.0*AssPtBin};
	if (PhotElecWPartnerTrigger[l][AssPtBin])  fTPartnerNoTPt2->Fill(fillSparse, MotherWeight[l]);
	if (PhotElecWoPartnerTrigger[l][AssPtBin]) fNoPartnerNoTPt2->Fill(fillSparse, MotherWeight[l]);

	if (PhotElecWPartnerTrigger[l][AssPtBin] && TPNoT) TPNoT=kFALSE;
	if (PhotElecWoPartnerTrigger[l][AssPtBin] && NoPNoT) NoPNoT = kFALSE;
	}
      if (RedTrack->TruePartner() && TPNoT) {
	fTPartnerNoT->Fill(RedTrack->Pt(), 1./recEff); //cout << "TruePartnerNoTrigger" << endl;
      }
      if (!RedTrack->TruePartner() && NoPNoT) {
	fNoPartnerNoT->Fill(RedTrack->Pt(), 1./recEff); // cout << "NoPartnerNoTrigger" << endl;
      }
    }
  }
  
  for (int l=0; l<ntracks; l++) {
    for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
      if (NonElectronIsTrigger[l][AssPtBin]) fNonElecHadTrigger->Fill(NonElectronIsTriggerPt[l], AssPtBin, 1./NonElectronIsTriggerWeight[l]);
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


void AliAnalysisTaskHaHFECorrel::CorrelateHadronMixedEvent(Float_t mult, const AliVVertex* pVtx, Float_t maxPt, Int_t nMother, Int_t listMother[], Bool_t EvContTP, Bool_t EvContNTP) {
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
    passTrackCut = ChargedHadronTrackCuts(pVtx, track, nMother, listMother);
    if (!passTrackCut) continue;

    // pid cuts
    Bool_t passPIDCut=kFALSE;
    passPIDCut = ChargedHadronPIDCuts(track);
    if (!passPIDCut) continue;
   
    Double_t ptH=-9.,etaH =-9.,phiH=-9., recEffH=-9.;
    Int_t idH=-9;
    ptH = track->Pt();
    phiH = track->Phi();
    etaH = track->Eta();
    idH = track->GetID();
    recEffH = GetHadronRecEff(ptH, phiH, etaH, pVtx->GetZ());
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
	  recEff=GetElectronRecEff(ptMix,phiMix,etaMix,pVtx->GetZ()); // assumption same vtx as original event (same Vtx binning)
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
	  for (Int_t k=0; k<ls; k++)  fULSElecHaMixedEvent->Fill(fillSparse, -1./(recEff*recEffH));
	  for (Int_t k=0; k<uls; k++) fULSElecHaMixedEvent->Fill(fillSparse, 1./(recEff*recEffH));
	  
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

 void AliAnalysisTaskHaHFECorrel::CorrelateWithHadrons(AliVTrack* Vtrack, const AliVVertex* pVtx, Int_t nMother, Int_t listMother[], Bool_t FillHadron, Bool_t FillLP, Bool_t **NonElecIsTrigger, Double_t *NonElecIsTriggerPt, Double_t * NonElecIsTriggerWeight, Int_t NumElectronsInEvent) {
  
  // Trigger Hadron
  Double_t ptH=-9.,etaH =-9.,phiH=-9., recEffH=-9;
  Int_t idH=-9;
  ptH = Vtrack->Pt();
  phiH = Vtrack->Phi();
  etaH = Vtrack->Eta();
  idH = Vtrack->GetID();
  recEffH = GetHadronRecEff(ptH,phiH,etaH, pVtx->GetZ());
  if (recEffH<0) return;

  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();

  Bool_t Trigger[fAssPtElec_Nbins];
  for (Int_t AssPt=0; AssPt<fAssPtElec_Nbins; AssPt++) Trigger[AssPt]=kFALSE;           // flag to store information, if used as trigger particle
 

  // Track loop for hadron correlations
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
    passTrackCut = InclElecTrackCuts(pVtx, VtrackAssoHad, nMother, listMother);
    if (!passTrackCut) continue;

    // pid cuts
    Bool_t passPIDCut=kFALSE;
    passPIDCut = AssoHadronPIDCuts(VtrackAssoHad);
    if (!passPIDCut) continue;
   
    // associated particle
    Double_t ptAsso=-9.,etaAsso=-9.,phiAsso=-9.,recEffAsso=-9;
    ptAsso = VtrackAssoHad->Pt();
    phiAsso = VtrackAssoHad->Phi();
    etaAsso = VtrackAssoHad->Eta();
    recEffAsso = GetElectronRecEff(ptAsso, phiAsso, etaAsso, pVtx->GetZ());
    if (recEffAsso<0) continue;
    
    CheckHadronIsTrigger(ptAsso, Trigger);
    CheckElectronIsTrigger(ptH, NonElecIsTrigger[iTracks]);
    NonElecIsTriggerPt[iTracks]=ptAsso;
    NonElecIsTriggerWeight[iTracks]=recEffAsso;
    if(FillLP) {
      for (Int_t ptAssHad=0; ptAssHad<fAssPtHad_Nbins; ptAssHad++) {
	if (NonElecIsTrigger[iTracks][ptAssHad]) fNonElecLPTrigger->Fill(ptAsso, ptAssHad, 1./recEffAsso);
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
    if (FillHadron) fElecHaHa->Fill(fillSparse, 1./(recEffH*recEffAsso));
    if (FillLP)     fElecLPHa->Fill(fillSparse, 1./(recEffAsso));
  }

  if (FillHadron) {
    for (Int_t ptAssElec=0; ptAssElec<fAssPtElec_Nbins; ptAssElec++) {
      if (Trigger[ptAssElec]) fHadNonElecTrigger->Fill(ptH, ptAssElec, 1./(recEffH));
    }
  }
  if (FillLP){
    for (Int_t ptAssElec=0; ptAssElec<fAssPtElec_Nbins; ptAssElec++) {
      if (Trigger[ptAssElec]) fLPNonElecTrigger->Fill(ptH, ptAssElec, 1./recEffH);
    }
  }
}



//_______________________

void AliAnalysisTaskHaHFECorrel::CorrelateLP(AliVTrack* LPtrack,  const AliVVertex* pVtx, Int_t nMother, Int_t listMother[], TObjArray *RedTracksHFE) 
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
  recEffH = GetHadronRecEff(ptH, phiH, etaH, pVtx->GetZ());
  if (recEffH<0) return;

  CorrelateWithHadrons(LPtrack, pVtx, nMother, listMother, kFALSE, kTRUE, NonElectronIsTrigger, NonElectronIsTriggerPt, NonElectronIsTriggerWeight, RedTracksHFE->GetEntriesFast()); // correlate LPHadron (kFALSE, kTRUE);
		       
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
    recEffE = GetElectronRecEff(pt,phi,eta,pVtx->GetZ());
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
    fInclElecLP->Fill(fillSparse,1./(recEffE));
    for (Int_t j=0; j<ls; j++)  fULSElecLP->Fill(fillSparse, -1./(recEffE));
    for (Int_t j=0; j<uls; j++) fULSElecLP->Fill(fillSparse, 1./(recEffE));


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
      for (Int_t j=0; j<ls; j++) fElecLPULSNoPartner->Fill(fillSparse, -1./(recEffE));
      for (Int_t j=0; j<uls; j++) fElecLPULSNoPartner->Fill(fillSparse, 1./(recEffE));
    }

    if (fIsMC) {
      if (RedTrack->IsHadron()){
	CheckElectronIsTrigger(ptH, HadContIsTrigger[k]);
	fMCElecLPHadron->Fill(fillSparse, 1./(recEffE));
      }
      if (RedTrack->IsPhotonic()) {
	if (RedTrack->TruePartner()) {
	  CheckElectronIsTrigger(ptH, PhotElecWPartnerTrigger[k]);
	  fMCElecLPTruePartner->Fill(fillSparse, 1./(recEffE));
	}
	else {
	  CheckElectronIsTrigger(ptH,  PhotElecWoPartnerTrigger[k]);
	  fMCElecLPNoPartner->Fill(fillSparse, 1./(recEffE));
	}
      }
    }
  }
   
  //save LP as Trigger (existens of LP and HFE are satisfied in function call)
  for (Int_t PtAssElec=0; PtAssElec<fAssPtElec_Nbins; PtAssElec++) {
    if (LPTrigger[PtAssElec]) fLPElecTrigger->Fill(ptH, PtAssElec, 1.);
  }


  // Fill Trigger Histograms
  for (Int_t l=0; l<RedTracksHFE->GetEntriesFast(); l++) {
    AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(l);
    Double_t recEffE = GetElectronRecEff(RedTrack->Pt(), RedTrack->Phi(), RedTrack->Eta(), pVtx->GetZ());
    if (recEffE<0) continue;
    for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
      if (ElectronIsTrigger[l][AssPtBin]) {
	fElecLPTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin, RedTrack->Eta(), 1.0/recEffE);
	if( RedTrack->LS()>0)  fElecLPTriggerULS->Fill(ElectronIsTriggerPt[l], AssPtBin,   -RedTrack->LS()/recEffE);
	if( RedTrack->ULS()>0) fElecLPTriggerULS->Fill(ElectronIsTriggerPt[l], AssPtBin,  RedTrack->ULS()/recEffE);
      }
      if (ElectronIsTriggerNoP[l][AssPtBin]) {
	if( RedTrack->LS()>0)  fElecLPTriggerULSNoP->Fill(ElectronIsTriggerPt[l], AssPtBin, -RedTrack->LS()/recEffE);
	if( RedTrack->ULS()>0) fElecLPTriggerULSNoP->Fill(ElectronIsTriggerPt[l], AssPtBin,  RedTrack->ULS()/recEffE);
      }
      if (PhotElecWPartnerTrigger[l][AssPtBin]) fMCElecLPTruePartnerTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin, 1./recEffE);
      if (PhotElecWoPartnerTrigger[l][AssPtBin]) fMCElecLPNoPartnerTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin, 1./recEffE);
      if (HadContIsTrigger[l][AssPtBin]) fHadContLPTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin, 1./recEffE);
    }
  }

  for (int l=0; l<ntracks; l++) {
    for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
      if (NonElectronIsTrigger[l][AssPtBin]) fNonElecLPTrigger->Fill(NonElectronIsTriggerPt[l], AssPtBin, 1./NonElectronIsTriggerWeight[l]);
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

void AliAnalysisTaskHaHFECorrel::CorrelateLPMixedEvent(AliVTrack* LPtrack, Float_t mult, Float_t zVtx, Float_t maxPt, Bool_t EvContTP, Bool_t EvContNTP) 
{ 
  // leading Particle 
  Double_t ptH=-9.,etaH =-9.,phiH=-9.,recEffH=-9.;
  ptH = LPtrack->Pt();
  phiH = LPtrack->Phi();
  etaH = LPtrack->Eta();
  recEffH = GetHadronRecEff(ptH, phiH, etaH, zVtx);
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
	recEffE = GetElectronRecEff(ptMix, phiMix, etaMix, zVtx);
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
	for (Int_t k=0; k<ls; k++)  fULSElecLPMixedEvent->Fill(fillSparse, -1./(recEffE));
	for (Int_t k=0; k<uls; k++) fULSElecLPMixedEvent->Fill(fillSparse, 1./(recEffE));

	if (fIsMC) {
	  if (mixtrk->IsPhotonic()) {
	    if (mixtrk->TruePartner()) {
	      fTagLPMixedEvent->Fill(fillSparse, 1./(recEffE));
	    }
	    else if (!mixtrk->TruePartner()) {
	      fNonTagLPMixedEvent->Fill(fillSparse, 1./(recEffE));
	    }
	  }
       	}
      }
    }
  }
}


//_________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::InclElecTrackCuts(const AliVVertex *pVtx,AliVTrack *Vtrack,Int_t nMother, Int_t listMother[])
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
  fElectronTrackCuts->Fill(Vtrack->Pt(), 0);

  if (fIsAOD) { 
    if(!AODtrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE;
  }
  else if(!ESDkTrkGlobalNoDCA(Vtrack)) return kFALSE;

  fElectronTrackCuts->Fill(Vtrack->Pt(), 1);

  if(Vtrack->Pt()<0.25) return kFALSE;
  fElectronTrackCuts->Fill(Vtrack->Pt(), 2);

  if(Vtrack->Eta()<fMinElectronEta) return kFALSE;
  if(Vtrack->Eta()>fMaxElectronEta) return kFALSE;
  fElectronTrackCuts->Fill(Vtrack->Pt(), 3);


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


  fElectronTrackTPCNcls->Fill(Vtrack->Pt(), Vtrack->GetTPCNcls());
  fElectronTrackTPCNclsdEdx->Fill(Vtrack->Pt(), Vtrack->GetTPCsignalN());
  fElectronTrackITSNcls->Fill(Vtrack->Pt(),ITSNcls);
  fElectronTrackITSLayer->Fill(Vtrack->Pt(),ITSPointOnLayer0, ITSPointOnLayer1);
  fElectronTrackRefit->Fill(Vtrack->Pt(), Vtrack->IsOn(AliAODTrack::kTPCrefit), Vtrack->IsOn(AliAODTrack::kITSrefit));
  fElectronTrackDCA->Fill(Vtrack->Pt(), DCAr, DCAz);


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
  fElectronTrackCuts->Fill(Vtrack->Pt(), 4);
  if(Vtrack->GetTPCsignalN() < fTPCndEdxCut) return kFALSE ;
  fElectronTrackCuts->Fill(Vtrack->Pt(), 5);
						
  if (!Vtrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
  fElectronTrackCuts->Fill(Vtrack->Pt(), 6);
  //if(ietrack->GetTPCFoundFraction() < 0.6) return kFALSE;   
  fElectronTrackCuts->Fill(Vtrack->Pt(), 7);

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
  fElectronTrackCuts->Fill(Vtrack->Pt(), 8);
  
  if(ITSNcls < fITSnCut) return kFALSE;
  fElectronTrackCuts->Fill(Vtrack->Pt(), 9);

  if (!Vtrack->IsOn(AliAODTrack::kITSrefit)) return kFALSE ;   
  fElectronTrackCuts->Fill(Vtrack->Pt(), 10);

  if(!(ITSPointOnLayer0 || ITSPointOnLayer1)) return kFALSE; // now kAny
  fElectronTrackCuts->Fill(Vtrack->Pt(), 11);

  if(!(ITSPointOnLayer0 && ITSPointOnLayer1)) return kFALSE; // now kBoth
  fElectronTrackCuts->Fill(Vtrack->Pt(), 12);
  

  if(TMath::Abs(DCAr) > 1 || TMath::Abs(DCAz) > 2) return kFALSE;
  fElectronTrackCuts->Fill(Vtrack->Pt(), 13); 


  if (fITSSharedClusterCut>0.9) {
    if (ITSNclsShared>fITSSharedClusterCut) return kFALSE;
  }
  else {
    if (ITSFracSharedNCls>fITSSharedClusterCut) return kFALSE;
  }
  // cout << "Cut " << fITSSharedClusterCut << "\t" << ITSNclsShared << endl;
  fElectronTrackCuts->Fill(Vtrack->Pt(), 14); 
 


 if (fUseTRD) {
    
    if (fIsAOD && AODtrack->GetTRDntrackletsPID()<4) return kFALSE;
    else if (!fIsAOD && ESDtrack->GetTRDntrackletsPID()<4) return kFALSE;

    if (fIsAOD && AODtrack->GetTRDncls()/ AODtrack->GetTRDntrackletsPID()<17) return kFALSE;
    else if (!fIsAOD && ESDtrack->GetTRDncls()/ESDtrack->GetTRDntrackletsPID()<17) return kFALSE;
  }


 


  return kTRUE;

 }

 //________________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::InclElecPIDCuts(AliVTrack* Vtrack, Bool_t IsPrimary) {

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
  
  
  if ((fITSnSigma < (-1.*fSigmaITScut) ) || (fITSnSigma > fSigmaITScut)) PassITSCut=kFALSE;
  if ((fTOFnSigma < (-1.*fSigmaTOFcut) ) || (fTOFnSigma > fSigmaTOFcut)) PassTOFCut=kFALSE;
  if ((fTPCnSigma<fSigmaTPCcut)  || (fTPCnSigma>3)) PassTPCCut=kFALSE; 


  if (IsPrimary) { //Fill PID histograms only in first round
    fHistITSnSig->Fill(Vtrack->P(),fITSnSigma);
    fHistTOFnSig->Fill(Vtrack->P(),fTOFnSigma);
    fHistTPCnSig->Fill(Vtrack->P(),fTPCnSigma);
    
    if (PassITSCut) fHistTPCnSigITScut->Fill(Vtrack->P(),fTPCnSigma);
    if (PassTOFCut) fHistTPCnSigTOFcut->Fill(Vtrack->P(),fTPCnSigma);
    if (PassITSCut && PassTOFCut) fHistTPCnSigITSTOFcut->Fill(Vtrack->P(),fTPCnSigma);

    if (PassTOFCut && fHadCont) fHadContPvsPt->Fill(Vtrack->P(), Vtrack->Pt());

    Double_t fillSparse[4];   
    fillSparse[0]=Vtrack->P();
    fillSparse[1]=Vtrack->Phi();
    fillSparse[2]=Vtrack->Eta();
    fillSparse[3]=fTPCnSigma;

    if (PassTOFCut && fHadCont) fHadContPPhiEtaTPC->Fill(fillSparse);

    fillSparse[1]=fITSnSigma;
    fillSparse[2]=fTOFnSigma;
    if (fHadCont) fHadContamination->Fill(fillSparse);

    fillSparse[0]=Vtrack->Pt();
    //    if (fHadCont) fHadContaminationPt->Fill(fillSparse);

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
      fHadContMC->Fill(fillSparse);
	   // fillSparse[0]=AODtrack->Pt();
      //      fHadContMCPt->Fill(fillSparse);
    }
   }
 
  if (!fUseITS) PassITSCut = kTRUE;
  if (!PassITSCut || !PassTOFCut || !PassTPCCut) return kFALSE;
  return kTRUE;
}




// _________________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::PhotElecPIDCuts(AliVTrack* Vtrack) {
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
Bool_t AliAnalysisTaskHaHFECorrel::PhotElecTrackCuts(const AliVVertex *pVtx,AliVTrack *Vtrack,Int_t nMother, Int_t listMother[])
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
    //if(Vtrack->Pt() < fPhotElecPtCut) return kFALSE;
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
      //fPhotTrackITSCuts->Fill(fillSparse);
    }

    if(!(FBTPCOnly || FBITSsa) ) return kFALSE;
    if(FBITSsa) {
      if (AODtrack->GetITSNcls()<5) return kFALSE;
      if (fPhotElecITSrefitCut && !(Vtrack->GetStatus()&AliESDtrack::kITSrefit)) return kFALSE;
    }
    if(FBTPCOnly) {
      if (Vtrack->GetTPCNcls() < fPhotElecTPCnCut) return kFALSE;
      if (fPhotElecITSrefitCut && !(Vtrack->GetStatus()&AliESDtrack::kITSrefit)) return kFALSE;
      if (fPhotElecITSrefitCut && !(Vtrack->GetStatus()&AliESDtrack::kTPCrefit)) return kFALSE;
    }

    return kTRUE;
  }
  return kFALSE; // if not AOD
}

 //_________________________________________

Bool_t AliAnalysisTaskHaHFECorrel::ChargedHadronTrackCuts(const AliVVertex *pVtx,AliVTrack *Vtrack,Int_t nMother, Int_t listMother[])
 {
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

   if(Vtrack->Pt()<0.25) return kFALSE;
   if(Vtrack->Eta()>fMaxHadronEta) return kFALSE;
   if(Vtrack->Eta()<fMinHadronEta) return kFALSE;

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
   //   cout << DCAr << "\t" << DCAz << endl;

   fHadronTrackTPCNcls->Fill(Vtrack->Pt(), Vtrack->GetTPCNcls());
   fHadronTrackRefit->Fill(Vtrack->Pt(), Vtrack->IsOn(AliAODTrack::kTPCrefit), Vtrack->IsOn(AliAODTrack::kITSrefit));
   fHadronTrackDCA->Fill(Vtrack->Pt(), DCAr, DCAz);


   if(fHITSrefitCut && !(Vtrack->GetStatus()&AliESDtrack::kITSrefit)) return kFALSE;
   if(fHTPCrefitCut && !(Vtrack->GetStatus()&AliESDtrack::kTPCrefit)) return kFALSE;
   
   if(Vtrack->GetTPCNcls() < fHTPCnCut) return kFALSE;
   
   if(TMath::Abs(DCAr) > 0.25 || TMath::Abs(DCAz) > 1) {
     //   cout << "Hadron fails DCA cut " << Vtrack->Pt() << endl;
     return kFALSE;
   }  

   return kTRUE;   
}


//_________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::ChargedHadronPIDCuts(AliVTrack *Vtrack)
{   
  return kTRUE;   
}


Bool_t AliAnalysisTaskHaHFECorrel::AssoHadronPIDCuts(AliVTrack *Vtrack) 
{
  // associated particle variables
  Double_t fTPCnSigmaAsso=-9.;
  Double_t fTOFnSigmaAsso=-9.;
  


  fTOFnSigmaAsso = fpidResponse->NumberOfSigmasTOF(Vtrack, AliPID::kElectron);
   if ((fTOFnSigmaAsso < (-1.*fSigmaTOFcut) ) || (fTOFnSigmaAsso > fSigmaTOFcut)) return kFALSE;
  // looser PID cuts
  fTPCnSigmaAsso = fpidResponse->NumberOfSigmasTPC(Vtrack, AliPID::kElectron);
  if(fTPCnSigmaAsso>fAssNonEleTPCcut) return kFALSE;
   
  return kTRUE;
}

void AliAnalysisTaskHaHFECorrel::MCEfficiencyCorrections(const AliVVertex * RecVertex) {
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

	// Pi0 = 111. Pi+ 211. Eta 221
	if (mcPDG==111 || mcPDG==221 || mcPDG== 211) {
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
	    if(genname.Contains("pi0") || genname.Contains("etagit")) {
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
	fMCHadPtEtaPhiVtx->Fill(fillSparse);
      }
      
      if (abs(mcPDG)==11) {
	if (mcMotherID>0) {
	  AliAODMCParticle* mcPartMother=(AliAODMCParticle*)fMC->GetTrack(mcMotherID);
	  Int_t mcMotherPDG = abs(mcPartMother->GetPdgCode()); 
	  Int_t MIsHeavy  = Int_t (mcMotherPDG / TMath::Power(10, Int_t(TMath::Log10(mcMotherPDG))));
	  if (MIsHeavy>3 && mcEta> fMinElectronEta && mcEta < fMaxElectronEta) {
	    fillSparse[2]=mcPhi;
	    fMCElecPtEtaPhiVtx->Fill(fillSparse);
	    if (fOneTimeCheck && fRecEff) {
	      fMCElecPDG->Fill(mcMotherPDG);
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



void AliAnalysisTaskHaHFECorrel::EvaluateTaggingEfficiency(AliVTrack * Vtrack, Int_t LSPartner, Int_t ULSPartner, Bool_t  trueULSPartner) {
  if (fIsAOD) {
    AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);   
    if (!AODtrack) return;
    Int_t MClabel=AODtrack->GetLabel();
    AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(MClabel));     
    Int_t PDGCode = abs(MCParticle->GetPdgCode());
    Int_t PDGCodeMother=0;
    Int_t PDGCodeGrandMother=0;
    Int_t PDGCodeGGMother=0;
    Double_t PtMother=0; 
    Double_t PtGrandMother=0;

    // get ancestors (eta-pi0-e, eta-gamma-e, pi0-gamma-e, eta-e, pi0-e, ...) - remember pt of pi0/Eta for correction
    if (MCParticle->GetMother()>=0) {
      AliAODMCParticle* MCParticleMother = (AliAODMCParticle*) fMC->GetTrack(MCParticle->GetMother());
      PDGCodeMother = abs(MCParticleMother->GetPdgCode());
      PtMother = MCParticleMother->Pt();
      if (MCParticleMother->GetMother()>=0) {
	AliAODMCParticle* MCParticleGrandMother=(AliAODMCParticle*) fMC->GetTrack(MCParticleMother->GetMother());
	PDGCodeGrandMother=abs(MCParticleGrandMother->GetPdgCode());
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
  
 

    Double_t fillSparse[5]={-999,-999,-999, -999, -999};
    fillSparse[0]=Vtrack->Pt(); // sparse will be filled with rec pt 
    if (PDGCode==11) fillSparse[1]=1.; //PDGMap.find(PDGCode)->second;
    else (fillSparse[1]=0.);
    fillSparse[2]=PDGMap.find(PDGCodeMother)->second;
    fillSparse[3]=PDGMap.find(PDGCodeGrandMother)->second; 
    if (PDGCodeGGMother<400) fillSparse[4]=1;
    else fillSparse[4]=0.; //PDGMap.find(PDGCodeGGMother)->second;
  
 
    Double_t PtMotherWeight=1.; 
    // Double_t ExcludePion[9] = {130,211, 310, 321, 2212, 3122,3222,3322,5122}; // K0L, PI+, K0S, K+, p, Lambd0, Sigma+, Xci0, lambda b0
    // Double_t ExcludeEta[0] = {};
    if (PDGCode==11) {
      if (PDGCodeMother==22) {
	if (PDGCodeGrandMother == 221) {
	  PtMotherWeight=GetEtaWeight(PtGrandMother); 
	}
	else if (PDGCodeGrandMother == 111) {
	  PtMotherWeight=GetPionWeight(PtGrandMother); 
	  //	for (Int_t j=0; j<9; j++) if (PDGCodeGGMother==ExcludePion[j]) PtMotherWeight=1.;
	}
      }
      else  {
	if (PDGCodeMother == 221) {
	  PtMotherWeight=GetEtaWeight(PtMother); 
	}
	else if (PDGCodeMother == 111) {
	  PtMotherWeight=GetPionWeight(PtMother); 
	  //for (Int_t j=0; j<9; j++) if (PDGCodeGGMother==ExcludePion[j]) PtMotherWeight=1.;
	}
      }
    }
  
    // real tag corr sparse
    fTagEffIncl->Fill(fillSparse, PtMotherWeight);
    if (LSPartner>0)  fTagEffULS->Fill(fillSparse, -LSPartner*PtMotherWeight);
    if (ULSPartner>0)  fTagEffULS->Fill(fillSparse, ULSPartner*PtMotherWeight);
    if (trueULSPartner) fTagTruePairs->Fill(fillSparse, PtMotherWeight);

    fTagEffInclWoWeight->Fill(fillSparse);
    if (LSPartner>0)  fTagEffULSWoWeight->Fill(fillSparse, -LSPartner);
    if (ULSPartner>0)  fTagEffULSWoWeight->Fill(fillSparse, ULSPartner);
    if (trueULSPartner) fTagTruePairsWoWeight->Fill(fillSparse);

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
      fTagMotherPt->Fill(fillSparse, PtMotherWeight);
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
     if (PDGCodeMother == 111 || PDGCodeMother == 221) return kTRUE; // true if mother is eta or pion
     else if (PDGCodeMother == 22 && MCParticleMother->GetMother()>=0) {
       AliAODMCParticle* MCParticleGrandMother=(AliAODMCParticle*) fMC->GetTrack(MCParticleMother->GetMother());
       PDGCodeGrandMother=abs(MCParticleGrandMother->GetPdgCode());
       if (PDGCodeGrandMother == 111 || PDGCodeGrandMother == 221) {
	 AliAODMCParticle* MCParticleGGMother=(AliAODMCParticle*) fMC->GetTrack(MCParticleMother->GetMother());
	 PDGCodeGGMother=abs(MCParticleGGMother->GetPdgCode());
	 if (PDGCodeGGMother<399) return kTRUE; // ignore HF feed-down 
	 else return kFALSE;
       }
       else return kFALSE;
     }
     else return kFALSE;
  }
  else return kFALSE;
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

void AliAnalysisTaskHaHFECorrel::TRDQA(Int_t RunNumber, const AliVVertex *pVtx, Int_t nMother, Int_t listMother[]) {

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
    passHFETrackCut= InclElecTrackCuts(pVtx,Vtrack,nMother,listMother);

    Bool_t passHFEPIDCut=kFALSE;
    if (passHFETrackCut) passHFEPIDCut= InclElecPIDCuts(Vtrack, kFALSE);

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
  Int_t Bin = fCorrectEtatoData.FindBin(pt);
  if (fCorrectEtatoData.IsBinUnderflow(Bin) || fCorrectEtatoData.IsBinOverflow(Bin)) {
    return 0;
  }
  return  fCorrectEtatoData.GetBinContent(Bin);
}

Double_t AliAnalysisTaskHaHFECorrel::GetPionWeight(Double_t pt) {
  Int_t Bin = fCorrectPiontoData.FindBin(pt);
  if (fCorrectPiontoData.IsBinUnderflow(Bin) || fCorrectPiontoData.IsBinOverflow(Bin)) {
    return 0;
  }
  return  fCorrectPiontoData.GetBinContent(Bin);
}


Double_t AliAnalysisTaskHaHFECorrel::GetHadronRecEff(Double_t pt, Double_t phi, Double_t eta, Double_t zVtx) {
  if (pt<0.5) return -1.;
  
  Int_t Bin = fHadRecEff.FindBin(pt,eta,zVtx);
  if (fHadRecEff.IsBinUnderflow(Bin) || fHadRecEff.IsBinOverflow(Bin) ) {
    return -1.;
  }
  Double_t RecEff = fHadRecEff.GetBinContent(Bin);
  if (RecEff>0.05) return RecEff;
  else {
    return -1.;
  }
}

Double_t AliAnalysisTaskHaHFECorrel::GetElectronRecEff(Double_t pt, Double_t phi, Double_t eta, Double_t zVtx) {

  //  return 1.;

  // Apply electron rec eff offline!

  if (pt<0.5) return -1;
  Int_t Bin = fEleRecEff.FindBin(zVtx, pt);
  //  cout << "EleBin " << Bin << endl;
  if (fEleRecEff.IsBinUnderflow(Bin) || fEleRecEff.IsBinOverflow(Bin)) {
    // cout <<  "ElecRecEff: " << pt << "\t" << eta << "\t" << zVtx << endl;
    return -1.;
  }
  Double_t RecEff= fEleRecEff.GetBinContent(Bin);
  if (RecEff>0.05) return RecEff;
  else {
    //cout <<  "ElecRecEff: " << pt << "\t" << eta << "\t" << zVtx << endl;
    return -1.;
  }
}

Double_t AliAnalysisTaskHaHFECorrel::GetNonTagCorr(Double_t ptTrack, Double_t ptAsso) {
  Double_t CorrFactor = -1.;
  // Int_t Bin = fNonTagCorr.FindBin(ptAsso);
  //  if (fNonTagCorr.IsBinUnderflow(Bin) || fNonTagCorr.IsBinOverflow(Bin)) {
  //    return -1.;
  //  }
  // CorrFactor = fNonTagCorr.GetBinContent(Bin);
 
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
    if (ptAsso<0.15) return 1.;
    else return 0;

  }


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
}


Bool_t AliAnalysisTaskHaHFECorrel::PassEventBias( const AliVVertex *pVtx, Int_t nMother, Int_t *listMother) {
 
  
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
      passHadTrackCut = ChargedHadronTrackCuts(pVtx, Vtrack, nMother, listMother);
   
      Bool_t passHadPIDCut=kFALSE;
      passHadPIDCut = ChargedHadronPIDCuts(Vtrack); // currently empty
    
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
  if (!fIsAOD && !ESDtrack) return kFALSE;
  fesdTrackCuts->SetMaxDCAToVertexXY(2.4);
  fesdTrackCuts->SetMaxDCAToVertexZ(3.2);
  fesdTrackCuts->SetDCAToVertex2D(kTRUE);
  if (fesdTrackCuts->AcceptTrack(ESDtrack)) return kTRUE;
  else return kFALSE;
}

void AliAnalysisTaskHaHFECorrel::MCTruthCorrelation(Bool_t AfterEventCuts, Int_t RecLPLabel,  Float_t pVtxZ, Float_t mult,  Int_t &LeadingParticleInAcceptance, Int_t &LeadingParticle) {
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

  TObjArray* MCTrueRedTracks = new TObjArray;
  MCTrueRedTracks->SetOwner(kTRUE);
 
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
	  if (MotherPDG ==11) { // e from e
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
	    if (AfterEventCuts && fOneTimeCheck) fTrueElectronEta->Fill(MCElectron->Pt(), MCElectron->Eta());

	    if (AfterEventCuts && ElectronInAcceptanceCut) {
	      AliBasicParticleHaHFE * ElectronParticle  = 0;
	      AliExternalTrackParam ExtTrackParam;
	      ElectronParticle = new AliBasicParticleHaHFE(MotherIsHeavy, MCElectron->Eta(), MCElectron->Phi(), MCElectron->Pt(), MCElectron->Charge(), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,ExtTrackParam);  // Original mixed event particle class: id has been filled with c or b flag
	      if (ElectronParticle) MCTrueRedTracks->Add(ElectronParticle);
	    }

	    Int_t CharmOrBeauty = MotherIsHeavy*5-20;
	    LeadingParticleInAcceptance = -999;
	    LeadingParticlePtInAcceptance = -999;
	    LeadingParticle = -999;
	    LeadingParticlePt = -999;
	    for (Int_t iMCHadron = 1; iMCHadron< fMC->GetNumberOfTracks(); iMCHadron++) {
	      AliAODMCParticle* MCHadron = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(iMCHadron));  
	      if (!MCHadron) continue;
	      if ((!MCHadron->IsPhysicalPrimary()) || ( MCHadron->Charge()==0)) continue; 
	      if (AfterEventCuts && fOneTimeCheck) fTrueHadronEta->Fill(MCHadron->Pt(), MCHadron->Eta());

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
		  if (AfterEventCuts) fTrueMCHadronEventCuts->Fill(FillSparse);
		  else fTrueMCHadron->Fill(FillSparse);
		}
		if (AfterEventCuts) {
		  FillSparse[4]=pVtxZ;
		  fTrueMCHadronEventCutsZvtx->Fill(FillSparse);
		}
	      }
	      if (ElectronInAcceptanceCut) {
		FillSparse[4]=CharmOrBeauty+2;
		CheckElectronIsTrigger(MCHadron->Pt(), ElectronIsTrigger[CharmOrBeauty+2]);
		if (fOneTimeCheck) {
		  if (AfterEventCuts) fTrueMCHadronEventCuts->Fill(FillSparse);
		  else fTrueMCHadron->Fill(FillSparse);
		}
	      }
	      if (HadronInAcceptanceCut) {
		FillSparse[4]=CharmOrBeauty+3;
		CheckElectronIsTrigger(MCHadron->Pt(), ElectronIsTrigger[CharmOrBeauty+3]);
		if (fOneTimeCheck) {
		  if (AfterEventCuts) fTrueMCHadronEventCuts->Fill(FillSparse);
		  else fTrueMCHadron->Fill(FillSparse);
		}
	      }
	    }
	   
	    
	    // fill trigger histogram for elec - hadron
	    for (Int_t i=0; i<11; i++) { // case
	      for (Int_t j=0; j<fAssPtHad_Nbins; j++) {
		if (ElectronIsTrigger[i][j]) {
		  if (AfterEventCuts) {
		    fTrueMCElecHaTriggerEventCuts->Fill(MCElectron->Pt(), i, j);
		  }
		  else if (fOneTimeCheck) {
		    fTrueMCElecHaTrigger->Fill(MCElectron->Pt(), i, j);
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
		if (AfterEventCuts && fOneTimeCheck) fTrueLPinAcceptanceEta->Fill(MCLPinAcceptance->Pt(), MCLPinAcceptance->Eta());

	
		FillSparseLP[0]=MCLPinAcceptance->Pt();
		FillSparseLP[1]=MCElectron->Pt();
		FillSparseLP[2]=GetDeltaPhi(MCElectron->Phi(), MCLPinAcceptance->Phi());
		FillSparseLP[3]=GetDeltaEta(MCElectron->Eta(), MCLPinAcceptance->Eta());
	      
		if (ElectronInAcceptanceCut) { //electron in acceptancen and LP within particles in acceptance
		  FillSparseLP[4]=CharmOrBeauty+1;
		  CheckElectronIsTrigger(MCLPinAcceptance->Pt(), ElectronIsTrigger[CharmOrBeauty+1]);
		  if (fOneTimeCheck) {		 
		    if (AfterEventCuts) fTrueMCLPEventCuts->Fill(FillSparseLP);
		    else fTrueMCLP->Fill(FillSparseLP);
		  }
		}
		FillSparseLP[4]=CharmOrBeauty+3; // all eclectron with LP from acceptance
		CheckElectronIsTrigger(MCLPinAcceptance->Pt(), ElectronIsTrigger[CharmOrBeauty+3]);
		if (fOneTimeCheck) {
		  if (AfterEventCuts) fTrueMCLPEventCuts->Fill(FillSparseLP);
		  else fTrueMCLP->Fill(FillSparseLP);
		}
	      }	 

	      if (MCElectron->Label() != LeadingParticle) {
		AliAODMCParticle* MCLP = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(LeadingParticle));  
		if(AfterEventCuts && fOneTimeCheck) fTrueLPEta->Fill(MCLP->Pt(), MCLP->Eta());
		
		FillSparseLP[0]=MCLP->Pt();
		FillSparseLP[2]=GetDeltaPhi(MCElectron->Phi(), MCLP->Phi());
		FillSparseLP[3]=GetDeltaEta(MCElectron->Eta(), MCLP->Eta());
		if (ElectronInAcceptanceCut) { // electron in acceptance with true LP
		  FillSparseLP[4]=CharmOrBeauty+2;
		  CheckElectronIsTrigger(MCLP->Pt(), ElectronIsTrigger[CharmOrBeauty+2]);
		  if (fOneTimeCheck) {
		    if (AfterEventCuts) fTrueMCLPEventCuts->Fill(FillSparseLP);
		    else fTrueMCLP->Fill(FillSparseLP);
		  }
		}
		
		FillSparseLP[4]=CharmOrBeauty+4;
		CheckElectronIsTrigger(MCLP->Pt(), ElectronIsTrigger[CharmOrBeauty+4]);
		if (fOneTimeCheck) { // neighter in acceptance
		  if (AfterEventCuts) fTrueMCLPEventCuts->Fill(FillSparseLP);
		  else fTrueMCLP->Fill(FillSparseLP);
		}
	      }

	      if (RecLPLabel>0 && ( abs(MCElectron->Label()) != abs(RecLPLabel))) {
		AliAODMCParticle* MCLPRec = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(RecLPLabel));  
		if(AfterEventCuts) fRecLPEta->Fill(MCLPRec->Pt(), MCLPRec->Eta());
		
		FillSparseLP[0]=MCLPRec->Pt();
		FillSparseLP[2]=GetDeltaPhi(MCElectron->Phi(), MCLPRec->Phi());
		FillSparseLP[3]=GetDeltaEta(MCElectron->Eta(), MCLPRec->Eta());
		if (ElectronInAcceptanceCut) { // electron in acceptance with rec LP
		  FillSparseLP[4]=CharmOrBeauty+5;
		  CheckElectronIsTrigger(MCLPRec->Pt(), ElectronIsTrigger[CharmOrBeauty+5]);
		  if (fOneTimeCheck) {
		    if (AfterEventCuts) fTrueMCLPEventCuts->Fill(FillSparseLP);
		    else fTrueMCLP->Fill(FillSparseLP);
		  }
	
		  if (AfterEventCuts ) {
		    FillSparseLP[4]=pVtxZ;
		    fTrueMCLPEventCutsZvtx->Fill(FillSparseLP);
		  }
		}
	      }
	    }


	    // fill trigger histogram for electron leading particle
	    for (Int_t i=0; i<11; i++) { // case
	      for (Int_t j=0; j<fAssPtHad_Nbins; j++) {
		if (ElectronIsTrigger[i][j]) {
		  if (AfterEventCuts) {
		    fTrueMCElecLPTriggerEventCuts->Fill(MCElectron->Pt(), i, j);
		  }
		  else if (fOneTimeCheck) {
		    fTrueMCElecLPTrigger->Fill(MCElectron->Pt(), i, j);
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
	    fTrueMCHadronEventCutsZvtxMEv->Fill(FillSparse);
	    if (abs(MCHadron->Label()) == abs(RecLPLabel))  fTrueMCLPEventCutsZvtxMEv->Fill(FillSparse);
	  }
	}
      }   
    }     
    HFEPool->UpdatePool(MCTrueRedTracks);
  }
  else delete MCTrueRedTracks;
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


void AliAnalysisTaskHaHFECorrel::PhotULSLSElectronAcceptance(const AliVVertex *pVtx, Float_t mult,  Int_t nMother, Int_t listMother[]) {
  
  AliEventPool* fPool;
  fPool = fPoolMgr->GetEventPool(mult, pVtx->GetZ(), 0.); // Get the buffer associated with the current centrality and z-vtx
  //  fPool->SetDebug(kTRUE);
  if (!fPool)
  {
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f, maxPt = %f", mult, pVtx->GetZ(), 0.));
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
	  passAssoTrackCutPhot = PhotElecTrackCuts(pVtx,VtrackAsso,nMother,listMother);
	  if(!passAssoTrackCutPhot) continue;
   
	  Bool_t passAssoPIDCutPhot = kFALSE;
	  passAssoPIDCutPhot = PhotElecPIDCuts(VtrackAsso); 
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
	    fPhotMixULS->Fill(ptMix);
	  }
        
	  if(mass<fInvmassCut && fFlagLS){
	    fPhotMixLS->Fill(ptMix, chargeMix);
	  }
  
	}//track loop
      } // mixtracks
    }// mix events
  }	  
}
