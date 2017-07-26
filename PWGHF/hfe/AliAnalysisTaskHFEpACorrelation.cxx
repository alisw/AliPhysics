
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
//      Task for Heavy-flavour electron analysis in pPb collisions    //
//      (+ Electron-Hadron Jetlike Azimuthal Correlation)             //
//																	  //
//		version: Fev 6, 2016.	     						          //
//                                                                    //
//	    Authors 							                          //
//		Elienos Pereira de Oliveira Filho (epereira@cern.ch)	      //
//		Cristiane Jahnke		(cristiane.jahnke@cern.ch)            //
//      Henrique Zanoli (h.zanoli@cern.ch)                            //
//      Alexis Mas (aleximas@if.usp.br)                               //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliPhysicsSelection.h"
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"
#include "AliLog.h"
#include "AliStack.h"
#include "AliCentrality.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
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
#include "AliSelectNonHFE.h"
#include "AliHFEpidTPC.h"
#include "AliAnalysisTaskHFEpACorrelation.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"
#include "AliESDHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "TParticle.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TVector.h"
#include "stdio.h"
#include "iostream"
#include "fstream"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliEventPoolManager.h"
#include "TObjArray.h"
#include "AliMultSelection.h"
//include to use reader as Lucile does
#include "AliAODHeader.h"
#include "AliAnalysisUtils.h"
#include "AliGenEventHeader.h"



// --- ANALYSIS system ---
#include "AliCalorimeterUtils.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMixedEvent.h"
#include "AliAODCaloCluster.h"
#include "AliOADBContainer.h"
#include "AliAnalysisManager.h"

// --- Detector ---

//______________________________________________________________________

//______________________________________________________________________
ClassImp(AliAnalysisTaskHFEpACorrelation)

//______________________________________________________________________
AliAnalysisTaskHFEpACorrelation::AliAnalysisTaskHFEpACorrelation(const char *name)
: AliAnalysisTaskSE(name)
,fCorrelationFlag(0)
,fUseKF(kFALSE)
,fIspp(kFALSE)
,fIsMC(0)
,fUseAlternativeBinnig(kFALSE)
,fAssocWithSPD(kFALSE)
,fIsHFE1(kFALSE)
,fIsHFE2(kFALSE)
,fIsNonHFE(kFALSE)
,fIsFromD(kFALSE)
,fIsFromB(kFALSE)
,fIsFromPi0(kFALSE)
,fIsFromEta(kFALSE)
,fIsFromGamma(kFALSE)
,fpTBins(0)
,fESD(0)
,fAOD(0)
,fVevent(0)
,fPartnerCuts(0)
,fOutputList(0)
,fPidResponse(0)
,fNonHFE(0)
,fIsAOD(kTRUE)
,fCentrality(0)
,fCentralityMin(0)
,fCentralityMax(100)
,fHasCentralitySelection(kFALSE)
,fCentralityHist(0)
,fCentralityHistPass(0)
,fZvtx(0)
,fEstimator(0)
,fUseDCACutforHadrons(kTRUE)
,fEffHadron(0)
,fNevent(0)
,fNevent2(0)
,fPtElec_Inc(0)
,fPtElec_ULS(0)
,fPtElec_LS(0)
,fPtElec_ULS_NoPid(0)
,fPtElec_LS_NoPid(0)
,fPtElec_ULS_weight(0)
,fPtElec_LS_weight(0)
,fPtElec_ULS2_weight(0)
,fPtElec_LS2_weight(0)
,fTOF01(0)
,fTOF02(0)
,fTOF03(0)
,fpid(0)
,fTPC_pt(0)
,fTPC_p(0)
,fTPC_momentum(0)
,fTPC_eta(0)
,fTPC_momentum1(0)
,fTPC_eta1(0)
,fTPCnsigma_pt(0)
,fTOFTPCnsigma_pt(0)
,fTPCnsigma_p(0)
,fTPCnsigma_eta(0)
,fTPCnsigma_phi(0)
,fVtxZ(0)
,fVtxZ_new1(0)
,fVtxZ_new2(0)
,fVtxZ_new3(0)
,fVtxZ_new4(0)
,fzRes1(0)
,fzRes2(0)
,fSPD_track_vtx1(0)
,fSPD_track_vtx2(0)
,fEtad(0)
,fNTracks(0)
,fTrack_Multi(0)
,fNTracks_pt(0)
,fNTracks_eta(0)
,fNTracks_phi(0)
,fTPCNcls_pid(0)
,fCEtaPhi_Inc(0)
,fCEtaPhi_ULS_Weight(0)
,fCEtaPhi_LS_Weight(0)
,fCEtaPhi_ULS_NoP_Weight(0)
,fCEtaPhi_LS_NoP_Weight(0)
,fInvMassULS(0)
,fInvMassLS(0)
,fDCA(0)
,fDCABack(0)
,fOpAngle(0)
,fOpAngleBack(0)
,fInvMass2(0)
,fInvMass2_weight(0)
,fInvMassBack2(0)
,fInvMassBack2_weight(0)
,fInvMass_pT(0)
,fInvMassBack_pT(0)
,fDCA2(0)
,fDCABack2(0)
,fOpAngle2(0)
,fOpAngleBack2(0)
,fMassCut(0.1)
,fEtaCutMin(-0.8)
,fEtaCutMax(0.8)
,fMinpTElec(0.5)
,fMaxpTElec(6.)
,fChi2Cut(3.5)
,fDCAcutzHadron(0.25)
,fDCAcutrHadron(1)
,fDCAcut(999)//dca between two tracks
,fDCAcutr(999)//dca to vertex
,fDCAcutz(999)//dca to vertex
,fMassCutFlag(kTRUE)
,fAngleCutFlag(kFALSE)
,fAngleCut(999)
,fChi2CutFlag(kFALSE)
,fDCAcutFlag(kFALSE)
,fAssHadronPtMin(0.3)
,fAssHadronPtMax(2.0)
,fPtBackgroundBeforeReco(0)
,fPtMinAsso(0.0)
,fTpcNclsAsso(80)
,fCuts(0)
,fCFM(0)
,fPID()
//Lucile
,fPIDqa(0)
,fMCstack(0)
,fRejectKinkMother(kFALSE)
,fMCtrack(0)
,fMCtrackMother(0)
,fMCtrackGMother(0)
,fMCtrackGGMother(0)
,fMCtrackGGGMother(0)
,fMCarray(0)
,fMCheader(0)
,fMCparticle(0)
,fMCparticle2(0)
,fMCparticleMother(0)
,fMCparticleGMother(0)
,fMCparticleGGMother(0)
,fMCparticleGGGMother(0)
,fEventHandler(0)
,fMCevent(0)
,fPoolMgr(0)
,fPool(0)
,fTracksClone(0)
,fTracks(0)
,fCEtaPhi_Inc_EM(0)
,fCEtaPhi_ULS_Weight_EM(0)
,fCEtaPhi_LS_Weight_EM(0)
,fPoolNevents(0)
,fEventMixingFlag(0)
,fCEtaPhi_Inc_DiHadron(0)
,fPtTrigger_Inc(0)
,fAnalysisUtils(0)

//pT,eta,zvtx
,fNoEtaCutElectronGeneratedSignalPtEtaZvtx(0)
,fEtaCutElectronGeneratedSignalPtEtaZvtx(0)
//,fEtaCutElectronInclusiveGeneratedPtEtaZvtx(0)
//,fEtaCutElectronBKGGeneratedPtEtaZvtx(0)
//,fEtaCutElectronHFeGeneratedPtEtaZvtx(0)
//,fEtaCutElectronOtherGeneratedPtEtaZvtx(0)
//,fEtaCutElectronNoMotherGeneratedPtEtaZvtx(0)
,fEtaCutElectronInclusiveRecoPtEtaZvtx(0)
,fEtaCutElectronBKWithLabelULS(0)
,fEtaCutElectronBKWithLabelLS(0)
,fEtaCutElectronBKNoTag(0)
,fEtaCutElectronRecoHFEMC(0)
,fEtaCutElectronRecoOtherMC(0)
,fMissIDElectronsReco(0)
,fHadronsReco(0)
,fHadronsRecoPP(0)
,fHadronsGenerated(0)

//DPhi MC
,fCEtaPhiNoEtaCutInclusive(0)
,fCEtaPhiNoEtaCutBKG(0)
,fCEtaPhiNoEtaCutHFe(0)
,fCEtaPhiNoEtaCutHFeNoDCA(0)
,fCEtaPhiNoEtaCutOther(0)
,fCEtaPhiNoEtaCutNoMother(0)
,fCEtaPhiCutInclusive(0)
,fCEtaPhiCutBKG(0)
,fCEtaPhiCutHFe(0)
,fCEtaPhiCutOther(0)
,fCEtaPhiCutNoMother(0)
,fCEtaPhi_Back_MC_Tag(0)
,fCEtaPhi_Other_MC_Tag(0)
,fCEtaPhi_HFe_MC_Tag(0)
,fCEtaPhi_MC_NoMother_Tag(0)
,fDCAElectronXY(0)
,fDCAElectronZ(0)
,fElectronNoLabel(0)
,fElectronNoLabelULS(0)
,fElectronNoLabelLS(0)
,fEtaCutElectronBKULSMainSources(0)
,fEtaCutElectronBKLSMainSources(0)
,fEtaCutElectronBKULSOtherSources(0)
,fEtaCutElectronBKLSOtherSources(0)
,fEtaCutElectronHFEULS(0)
,fEtaCutElectronHFELS(0)
,fEtaCutElectronMissIDULS(0)
,fEtaCutElectronMissIDLS(0)
//Test NHFE weight
,fEtaCutElectronBKULSMainSources_NW(0)
,fEtaCutElectronBKLSMainSources_NW(0)
//Test Background weight (as Cris)
,fEtaCutElectronBKNoTag_WithMotherW(0)
,fEtaCutElectronBKULSMainSources_WithMotherW(0)
,fEtaCutElectronBKULSMainSources_WithMotherW_NW(0)
,fEtaCutElectronBKLSMainSources_WithMotherW(0)
,fEtaCutElectronBKLSMainSources_WithMotherW_NW(0)
,fUseGlobalTracksHadron(kTRUE)
,fCentralityValue(-1)
,fPtMCpi0_PureHijing(0)
,fPtMCpi0_NoMother(0)
,fPtMCEta_NoMother(0)
,fPtMCEta_PureHijing(0)
,fBkgPi0Weight(0)
,fBkgEtaWeight(0)
,fElectronBKGNoEnhULS(0)
,fElectronBKGNoEnhLS(0)
,fElectronBKGNoEnhTotalNumber(0)
,fElectronBKGWToDataTotal(0)
,fElectronBKGWToDataULS(0)
,fElectronBKGWToDataLS(0)
,fBkgPi0WeightToData(0)
,fBkgEtaWeightToData(0)
,fElectronBKGNoEnhULS_WithW(0)
,fElectronBKGNoEnhLS_WithW(0)
,fElectronBKGNoEnhTotalNumber_WithW(0)

{
    //Named constructor
    // Define input and output slots here
    // Input slot #0 works with a TChain
    
    //exotic
    //fEMCALRecoUtils  = new AliEMCALRecoUtils();
    fPID = new AliHFEpid("hfePid");
    
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    // DefineOutput(1, TH1I::Class());
    DefineOutput(1, TList::Class());
    //  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskHFEpACorrelation::AliAnalysisTaskHFEpACorrelation()
: AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTaskHFEpACorrelation")
,fCorrelationFlag(0)
,fUseKF(kFALSE)
,fIspp(kFALSE)
,fIsMC(0)
,fUseAlternativeBinnig(kFALSE)
,fAssocWithSPD(kFALSE)
,fIsHFE1(kFALSE)
,fIsHFE2(kFALSE)
,fIsNonHFE(kFALSE)
,fIsFromD(kFALSE)
,fIsFromB(kFALSE)
,fIsFromPi0(kFALSE)
,fIsFromEta(kFALSE)
,fIsFromGamma(kFALSE)
,fpTBins(0)
,fESD(0)
,fAOD(0)
,fVevent(0)
,fPartnerCuts(0)
,fOutputList(0)
,fPidResponse(0)
,fNonHFE(0)
,fIsAOD(kTRUE)
,fCentrality(0)
,fCentralityMin(0)
,fCentralityMax(100)
,fHasCentralitySelection(kFALSE)
,fCentralityHist(0)
,fCentralityHistPass(0)
,fZvtx(0)
,fEstimator(0)
,fUseDCACutforHadrons(kTRUE)
,fEffHadron(0)
,fNevent(0)
,fNevent2(0)
,fPtElec_Inc(0)
,fPtElec_ULS(0)
,fPtElec_LS(0)
,fPtElec_ULS_NoPid(0)
,fPtElec_LS_NoPid(0)
,fPtElec_ULS_weight(0)
,fPtElec_LS_weight(0)
,fPtElec_ULS2_weight(0)
,fPtElec_LS2_weight(0)
,fTOF01(0)
,fTOF02(0)
,fTOF03(0)
,fpid(0)
,fTPC_pt(0)
,fTPC_p(0)
,fTPC_momentum(0)
,fTPC_eta(0)
,fTPC_momentum1(0)
,fTPC_eta1(0)
,fTPCnsigma_pt(0)
,fTOFTPCnsigma_pt(0)
,fTPCnsigma_p(0)
,fTPCnsigma_eta(0)
,fTPCnsigma_phi(0)
,fVtxZ(0)
,fVtxZ_new1(0)
,fVtxZ_new2(0)
,fVtxZ_new3(0)
,fVtxZ_new4(0)
,fzRes1(0)
,fzRes2(0)
,fSPD_track_vtx1(0)
,fSPD_track_vtx2(0)
,fEtad(0)
,fNTracks(0)
,fTrack_Multi(0)
,fNTracks_pt(0)
,fNTracks_eta(0)
,fNTracks_phi(0)
,fTPCNcls_pid(0)
,fCEtaPhi_Inc(0)
,fCEtaPhi_ULS_Weight(0)
,fCEtaPhi_LS_Weight(0)
,fCEtaPhi_ULS_NoP_Weight(0)
,fCEtaPhi_LS_NoP_Weight(0)
,fInvMassULS(0)
,fInvMassLS(0)
,fDCA(0)
,fDCABack(0)
,fOpAngle(0)
,fOpAngleBack(0)
,fInvMass2(0)
,fInvMass2_weight(0)
,fInvMassBack2(0)
,fInvMassBack2_weight(0)
,fInvMass_pT(0)
,fInvMassBack_pT(0)
,fDCA2(0)
,fDCABack2(0)
,fOpAngle2(0)
,fOpAngleBack2(0)
,fMassCut(0.1)
,fEtaCutMin(-0.8)
,fEtaCutMax(0.8)
,fMinpTElec(0.5)
,fMaxpTElec(6.)
,fChi2Cut(3.5)
,fDCAcutzHadron(0.25)
,fDCAcutrHadron(1)
,fDCAcut(999)//dca between two tracks
,fDCAcutr(999)//dca to vertex
,fDCAcutz(999)//dca to vertex
,fMassCutFlag(kTRUE)
,fAngleCutFlag(kFALSE)
,fAngleCut(999)
,fChi2CutFlag(kFALSE)
,fDCAcutFlag(kFALSE)
,fAssHadronPtMin(0.3)
,fAssHadronPtMax(2.0)
,fPtBackgroundBeforeReco(0)
,fPtMinAsso(0.0)
,fTpcNclsAsso(80)
,fCuts(0)
,fCFM(0)
,fPID()
//Lucile
,fPIDqa(0)
,fMCstack(0)
,fRejectKinkMother(kFALSE)
,fMCtrack(0)
,fMCtrackMother(0)
,fMCtrackGMother(0)
,fMCtrackGGMother(0)
,fMCtrackGGGMother(0)
,fMCarray(0)
,fMCheader(0)
,fMCparticle(0)
,fMCparticle2(0)
,fMCparticleMother(0)
,fMCparticleGMother(0)
,fMCparticleGGMother(0)
,fMCparticleGGGMother(0)
,fEventHandler(0)
,fMCevent(0)
,fPoolMgr(0)
,fPool(0)
,fTracksClone(0)
,fTracks(0)
,fCEtaPhi_Inc_EM(0)
,fCEtaPhi_ULS_Weight_EM(0)
,fCEtaPhi_LS_Weight_EM(0)
,fPoolNevents(0)
,fEventMixingFlag(0)
,fCEtaPhi_Inc_DiHadron(0)
,fPtTrigger_Inc(0)
,fAnalysisUtils(0)

//pT,eta,zvtx
,fNoEtaCutElectronGeneratedSignalPtEtaZvtx(0)
,fEtaCutElectronGeneratedSignalPtEtaZvtx(0)
,fEtaCutElectronInclusiveRecoPtEtaZvtx(0)
,fEtaCutElectronBKWithLabelULS(0)
,fEtaCutElectronBKWithLabelLS(0)
,fEtaCutElectronBKNoTag(0)
,fEtaCutElectronRecoHFEMC(0)
,fEtaCutElectronRecoOtherMC(0)
,fMissIDElectronsReco(0)
,fHadronsReco(0)
,fHadronsRecoPP(0)
,fHadronsGenerated(0)

//DPhi MC
,fCEtaPhiNoEtaCutInclusive(0)
,fCEtaPhiNoEtaCutBKG(0)
,fCEtaPhiNoEtaCutHFe(0)
,fCEtaPhiNoEtaCutHFeNoDCA(0)
,fCEtaPhiNoEtaCutOther(0)
,fCEtaPhiNoEtaCutNoMother(0)
,fCEtaPhiCutInclusive(0)
,fCEtaPhiCutBKG(0)
,fCEtaPhiCutHFe(0)
,fCEtaPhiCutOther(0)
,fCEtaPhiCutNoMother(0)
,fCEtaPhi_Back_MC_Tag(0)
,fCEtaPhi_Other_MC_Tag(0)
,fCEtaPhi_HFe_MC_Tag(0)
,fCEtaPhi_MC_NoMother_Tag(0)
,fDCAElectronXY(0)
,fDCAElectronZ(0)
,fElectronNoLabel(0)
,fElectronNoLabelULS(0)
,fElectronNoLabelLS(0)
,fEtaCutElectronBKULSMainSources(0)
,fEtaCutElectronBKLSMainSources(0)
,fEtaCutElectronBKULSOtherSources(0)
,fEtaCutElectronBKLSOtherSources(0)
,fEtaCutElectronHFEULS(0)
,fEtaCutElectronHFELS(0)
,fEtaCutElectronMissIDULS(0)
,fEtaCutElectronMissIDLS(0)
//Test NHFE weight
,fEtaCutElectronBKULSMainSources_NW(0)
,fEtaCutElectronBKLSMainSources_NW(0)
//Test Background weight (as Cris)
,fEtaCutElectronBKNoTag_WithMotherW(0)
,fEtaCutElectronBKULSMainSources_WithMotherW(0)
,fEtaCutElectronBKULSMainSources_WithMotherW_NW(0)
,fEtaCutElectronBKLSMainSources_WithMotherW(0)
,fEtaCutElectronBKLSMainSources_WithMotherW_NW(0)
,fUseGlobalTracksHadron(kTRUE)
,fCentralityValue(-1)
,fPtMCpi0_PureHijing(0)
,fPtMCpi0_NoMother(0)
,fPtMCEta_NoMother(0)
,fPtMCEta_PureHijing(0)
,fBkgPi0Weight(0)
,fBkgEtaWeight(0)
,fElectronBKGNoEnhULS(0)
,fElectronBKGNoEnhLS(0)
,fElectronBKGNoEnhTotalNumber(0)
,fElectronBKGWToDataTotal(0)
,fElectronBKGWToDataULS(0)
,fElectronBKGWToDataLS(0)
,fBkgPi0WeightToData(0)
,fBkgEtaWeightToData(0)
,fElectronBKGNoEnhULS_WithW(0)
,fElectronBKGNoEnhLS_WithW(0)
,fElectronBKGNoEnhTotalNumber_WithW(0)

{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

//______________________________________________________________________
AliAnalysisTaskHFEpACorrelation::~AliAnalysisTaskHFEpACorrelation()
{
    //Destructor
    delete fOutputList;
    delete fPID;
    delete fCFM;
    delete fPIDqa;
    delete fNonHFE;
    delete fPartnerCuts;
    delete fAnalysisUtils;
}

//______________________________________________________________________
//Create Output Objects
//Here we can define the histograms and others output files
//Called once
void AliAnalysisTaskHFEpACorrelation::UserCreateOutputObjects()
{
    //______________________________________________________________________
    //Initialize PID
    
    fAnalysisUtils = new AliAnalysisUtils;
    //
    fPartnerCuts = new AliESDtrackCuts();
    
    
    if(!fPID->GetNumberOfPIDdetectors())
    {
        fPID->AddDetector("TPC", 0);
    }
    
    fPID->SortDetectors();
    
    fPIDqa = new AliHFEpidQAmanager();
    fPIDqa->Initialize(fPID);
    
    fNonHFE = new AliSelectNonHFE();
    //______________________________________________________________________
    
    
    //______________________________________________________________________
    //Initialize correction Framework and Cuts
    fCFM = new AliCFManager;
    const Int_t kNcutSteps = AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack + AliHFEcuts::kNcutStepsDETrack;
    fCFM->SetNStepParticle(kNcutSteps);
    for(Int_t istep = 0; istep < kNcutSteps; istep++) fCFM->SetParticleCutsList(istep, NULL);
    
    if(!fCuts)
    {
        AliWarning("================ \n Cuts not available. Default cuts will be used ================ \n");
        AliWarning("================ \n You should really consider to abort the analysis ================ \n");
        fCuts = new AliHFEcuts;
        fCuts->CreateStandardCuts();
    }
    
    //from Andrea Dubla
    //fCuts->SetAOD();
    
    fCuts->Initialize(fCFM);
    
    
    ///Output Tlist
    //Create TList
    fOutputList = new TList();
    fOutputList->SetOwner();
    
    
    
    //PIDqa
    //fOutputList->Add(fPIDqa->MakeList("PIDQA"));
    
    //Store the number of events
    //Define the histo
    fNevent = new TH1F("fNevent","Number of Events",30,0,30);
    fNevent2 = new TH1F("fNevent2","Number of Events 2",30,0,30);
    fOutputList->Add(fNevent);
    fOutputList->Add(fNevent2);
    //And then, add to the output list
    
    
    fpid = new TH1F("fpid","PID flag",5,0,5);
    fOutputList->Add(fpid);
    
    //pt Distribution
    fPtElec_Inc = new TH1F("fPtElec_Inc","Inclusive Electrons; p_{T} (GeV/c); Count",110,0.5,6);
    fPtElec_ULS = new TH1F("fPtElec_ULS","ULS; p_{T} (GeV/c); Count",110,0.5,6);
    fPtElec_LS = new TH1F("fPtElec_LS","LS; p_{T} (GeV/c); Count",110,0.5,6);
    
    if (fIsMC)
    {
        fPtElec_ULS_NoPid = new TH1F("fPtElec_ULS_NoPid","ULS; p_{T} (GeV/c); Count",110,0.5,6);
        fPtElec_LS_NoPid = new TH1F("fPtElec_LS_NoPid","LS; p_{T} (GeV/c); Count",110,0.5,6);
        fOutputList->Add(fPtElec_ULS_NoPid);
        fOutputList->Add(fPtElec_LS_NoPid);
    }
    
    fOutputList->Add(fPtElec_Inc);
    fOutputList->Add(fPtElec_ULS);
    fOutputList->Add(fPtElec_LS);
    
    //fTOF01 = new TH2F("fTOF01","",100,-20,20,100,-20,20);
    //fTOF02 = new TH2F("fTOF02","",100,-20,20,100,-20,20);
    //fTOF03 = new TH2F("fTOF03","",100,-20,20,100,-20,20);
    
    fPtTrigger_Inc = new TH1F("fPtTrigger_Inc","pT dist for Hadron Contamination; p_{t} (GeV/c); Count",110,0.5,6);
    
    fOutputList->Add(fPtTrigger_Inc);
    //fOutputList->Add(fTOF01);
    //fOutputList->Add(fTOF02);
    //fOutputList->Add(fTOF03);
    
    fVtxZ_new1= new  TH1F("fVtxZ_new1","fVtxZ_new1",100, -50,50);
    fVtxZ_new2= new  TH1F("fVtxZ_new2","fVtxZ_new2",100, -50,50);
    fVtxZ_new3= new  TH1F("fVtxZ_new3","fVtxZ_new3",100, -50,50);
    fVtxZ_new4= new  TH1F("fVtxZ_new4","fVtxZ_new4",100, -50,50);
    
    fOutputList->Add(fVtxZ_new1);
    fOutputList->Add(fVtxZ_new2);
    fOutputList->Add(fVtxZ_new3);
    fOutputList->Add(fVtxZ_new4);
    
    fzRes1= new  TH1F("fzRes1","fzRes1",100, 0,1);
    fzRes2= new  TH1F("fzRes2","fzRes2",100, 0,1);
    fSPD_track_vtx1= new  TH1F("fSPD_track_vtx1","fSPD_track_vtx1",100, -5,5);
    fSPD_track_vtx2= new  TH1F("fSPD_track_vtx2","fSPD_track_vtx2",100, -5,5);
    
    fOutputList->Add(fzRes1);
    fOutputList->Add(fzRes2);
    fOutputList->Add(fSPD_track_vtx1);
    fOutputList->Add(fSPD_track_vtx2);
    
    
    //General Histograms
    
    //Steps
    //Step 1: Before Track cuts
    //Step 2: Before PID
    //Step 3: After PID
    
    fTPC_p = new TH2F *[3];
    fTPCnsigma_p = new TH2F *[3];
    fVtxZ= new  TH1F *[3];
    fEtad= new  TH1F *[3];
    fNTracks= new  TH1F *[3];
    
    fNTracks_pt= new  TH2F *[3];
    fNTracks_eta= new  TH2F *[3];
    fNTracks_phi= new  TH2F *[3];
    fTPCNcls_pid = new  TH2F *[3];
    
    
    fTPC_momentum = new TH2F("fTPC_momentum",";p (GeV/c);TPC dE/dx (a. u.)",100,0,10,100,-20,20);
    fTPC_eta = new TH2F("fTPC_eta",";eta (GeV/c);TPC dE/dx (a. u.)",80,-2,2,220,-20,200);
    fOutputList->Add(fTPC_momentum);
    fOutputList->Add(fTPC_eta);
    
    fTPC_momentum1 = new TH2F("fTPC_momentum1",";p (GeV/c);TPC dE/dx (a. u.)",100,0,10,100,-20,20);
    fTPC_eta1 = new TH2F("fTPC_eta1",";eta (GeV/c);TPC dE/dx (a. u.)",80,-2,2,220,-20,200);
    fOutputList->Add(fTPC_momentum1);
    fOutputList->Add(fTPC_eta1);
    
    for(Int_t i = 0; i < 3; i++)
    {
        fTPC_p[i] = new TH2F(Form("fTPC_p%d",i),";p (GeV/c);TPC dE/dx (a. u.)",1000,0.3,15,1000,-20,200);
        fTPCnsigma_p[i] = new TH2F(Form("fTPCnsigma_p%d",i),";p (GeV/c);TPC Electron N#sigma",1000,0.3,15,1000,-15,10);
        
        fOutputList->Add(fTPC_p[i]);
        fOutputList->Add(fTPCnsigma_p[i]);
        
        fVtxZ[i]= new  TH1F(Form("fVtxZ%d",i),"VtxZ",100, -50,50);
        fEtad[i]= new  TH1F(Form("fEtad%d",i),"Eta distribution",100, -1.2,1.2);
        fNTracks[i]= new  TH1F(Form("fNTracks%d",i),"NTracks",1000, 0,5000);
        
        fOutputList->Add(fVtxZ[i]);
        fOutputList->Add(fEtad[i]);
        fOutputList->Add(fNTracks[i]);
        
        fNTracks_pt[i]= new  TH2F(Form("fNTracks_pt%d",i),"NTracks vs. p_{T}; NTracks; p_{T}",100, 0,5000, 100, 0, 10);
        fNTracks_eta[i]= new  TH2F(Form("fNTracks_eta%d",i),"NTracks vs. #eta; NTracks; #eta",100, 0,5000, 100, -1.0, 1.0);
        fNTracks_phi[i]= new  TH2F(Form("fNTracks_phi%d",i),"NTracks vs. #varph; NTracks; #varphi",100, 0,5000, 100, 0, TMath::Pi());
        
        
        fOutputList->Add(fNTracks_pt[i]);
        fOutputList->Add(fNTracks_eta[i]);
        fOutputList->Add(fNTracks_phi[i]);
        
    }
    
    for(Int_t i = 0; i < 4; i++)
    {
        fTPCNcls_pid[i]= new TH2F(Form("fTPCNcls_pid%d",i),"fTPCNcls_pid;NCls;NCls for PID",200,0,200,200,0,200);
        fOutputList->Add(fTPCNcls_pid[i]);
    }
    //pt bin
    
    
    fTPC_pt = new TH1F *[fpTBins.GetSize()];
    
    
    fTPCnsigma_pt = new TH1F *[fpTBins.GetSize()];
    fTOFTPCnsigma_pt = new TH2F *[fpTBins.GetSize()];
    
    
    fDCAElectronXY = new TH1F *[2];
    fDCAElectronZ = new TH1F *[2];
    
    for (Int_t i = 0; i < 2; i++)
    {
        fDCAElectronXY[i] = new TH1F(Form("fDCAElectronXY%d",i),"",40,0,4);
        fDCAElectronZ[i] = new TH1F(Form("fDCAElectronZ%d",i),"",40,0,4);
        fOutputList->Add(fDCAElectronXY[i]);
        fOutputList->Add(fDCAElectronZ[i]);
    }
    
    if (fIsMC)
    {
        //binnig for unfolding
        
        Double_t EtaBins[] = {-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}; //16 bins
        Int_t NEtaBins = 16;
        
        Double_t ZVtxBins[] = {-10, -7, -5, -3, -1, 1, 3, 5, 7, 10.01}; //9 bins
        Int_t nZvtxBinsMC = 9;
        
        Double_t pTBinsH[] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.8,2.0,2.5,3,3.5,4,5,6}; //21 bins
        Int_t NpTbinsH = 21;
        
        Int_t NumberofBins = fpTBins.GetSize();
        
        fCEtaPhiNoEtaCutHFe = new TH2F *[NumberofBins];
        fCEtaPhiCutHFe = new TH2F *[NumberofBins];
        
        fCEtaPhi_Back_MC_Tag = new TH2F *[NumberofBins];
        fCEtaPhi_Other_MC_Tag  = new TH2F *[NumberofBins];
        fCEtaPhi_HFe_MC_Tag = new TH2F *[NumberofBins];
        fCEtaPhi_MC_NoMother_Tag = new TH2F *[NumberofBins];
        
        fNoEtaCutElectronGeneratedSignalPtEtaZvtx = new TH1F("fNoEtaCutElectronGeneratedSignalPtEtaZvtx", "",110,0.5,6);
        fOutputList->Add(fNoEtaCutElectronGeneratedSignalPtEtaZvtx);
        
        fEtaCutElectronGeneratedSignalPtEtaZvtx = new TH1F("fEtaCutElectronGeneratedSignalPtEtaZvtx", "",110,0.5,6);
        fOutputList->Add(fEtaCutElectronGeneratedSignalPtEtaZvtx);
        
        
        
        fEtaCutElectronBKWithLabelULS = new TH1F("fEtaCutElectronBKWithLabelULS", "" ,110,0.5,6);
        
        fEtaCutElectronBKWithLabelLS = new TH1F("fEtaCutElectronBKWithLabelLS", "", 110,0.5,6);
        
        
        fEtaCutElectronInclusiveRecoPtEtaZvtx = new TH1F("fEtaCutElectronInclusiveRecoPtEtaZvtx", "", 110,0.5,6);
        //Backgound for data with MC info
        fEtaCutElectronBKNoTag = new TH1F("fEtaCutElectronBKNoTag", "", 110,0.5,6);
        
        //Tagged Background
        
        //HFe with MC information
        fEtaCutElectronRecoHFEMC = new TH1F("fEtaCutElectronRecoHFEMC", "", 110,0.5,6);
        //Others with MC information
        fEtaCutElectronRecoOtherMC = new TH1F("fEtaCutElectronRecoOtherMC", "", 110,0.5,6);
        
        //Miss ID
        fMissIDElectronsReco = new TH1F("fMissIDElectronsReco", "", 110,0.5,6);
        
        fElectronNoLabel = new TH1F("fElectronNoLabel", "", 110,0.5,6);
        fElectronNoLabelULS = new TH1F("fElectronNoLabelULS", "", 110,0.5,6);
        fElectronNoLabelLS = new TH1F("fElectronNoLabelLS", "", 110,0.5,6);
        fEtaCutElectronBKULSMainSources = new TH1F("fEtaCutElectronBKULSMainSources", "", 110,0.5,6);
        fEtaCutElectronBKLSMainSources = new TH1F("fEtaCutElectronBKLSMainSources", "", 110,0.5,6);
        fEtaCutElectronBKULSOtherSources = new TH1F("fEtaCutElectronBKULSOtherSources", "", 110,0.5,6);
        fEtaCutElectronBKLSOtherSources = new TH1F("fEtaCutElectronBKLSOtherSources", "", 110,0.5,6);
        fEtaCutElectronHFEULS = new TH1F("fEtaCutElectronHFEULS", "", 110,0.5,6);
        fEtaCutElectronHFELS = new TH1F("fEtaCutElectronHFELS", "", 110,0.5,6);
        fEtaCutElectronMissIDULS = new TH1F("fEtaCutElectronMissIDULS", "", 110,0.5,6);
        fEtaCutElectronMissIDLS = new TH1F("fEtaCutElectronMissIDLS", "", 110,0.5,6);
        
        //Test NHFE weight
        fEtaCutElectronBKULSMainSources_NW = new TH1F("fEtaCutElectronBKULSMainSources_NW", "", 110,0.5,6);
        fEtaCutElectronBKLSMainSources_NW = new TH1F("fEtaCutElectronBKLSMainSources_NW", "", 110,0.5,6);
        
        //Test Background weight (as Cris)
        fEtaCutElectronBKNoTag_WithMotherW = new TH1F("fEtaCutElectronBKNoTag_WithMotherW", "", 110,0.5,6);
        fEtaCutElectronBKULSMainSources_WithMotherW = new TH1F("fEtaCutElectronBKULSMainSources_WithMotherW", "", 110,0.5,6);
        fEtaCutElectronBKULSMainSources_WithMotherW_NW = new TH1F("fEtaCutElectronBKULSMainSources_WithMotherW_NW", "", 110,0.5,6);
        fEtaCutElectronBKLSMainSources_WithMotherW = new TH1F("fEtaCutElectronBKLSMainSources_WithMotherW", "", 110,0.5,6);
        fEtaCutElectronBKLSMainSources_WithMotherW_NW = new TH1F("fEtaCutElectronBKLSMainSources_WithMotherW_NW", "", 110,0.5,6);
        
        fElectronBKGNoEnhTotalNumber = new TH1F("fElectronBKGNoEnhTotalNumber", "", 110,0.5,6);
        fElectronBKGNoEnhULS = new TH1F("fElectronBKGNoEnhULS", "", 110,0.5,6);
        fElectronBKGNoEnhLS = new TH1F("fElectronBKGNoEnhLS", "", 110,0.5,6);
        
        fOutputList->Add(fElectronBKGNoEnhTotalNumber);
        fOutputList->Add(fElectronBKGNoEnhULS);
        fOutputList->Add(fElectronBKGNoEnhLS);
        
        fElectronBKGNoEnhTotalNumber_WithW = new TH1F("fElectronBKGNoEnhTotalNumber_WithW", "", 110,0.5,6);
        fElectronBKGNoEnhULS_WithW = new TH1F("fElectronBKGNoEnhULS_WithW", "", 110,0.5,6);
        fElectronBKGNoEnhLS_WithW = new TH1F("fElectronBKGNoEnhLS_WithW", "", 110,0.5,6);
        
        fOutputList->Add(fElectronBKGNoEnhTotalNumber_WithW);
        fOutputList->Add(fElectronBKGNoEnhULS_WithW);
        fOutputList->Add(fElectronBKGNoEnhLS_WithW);
        
        
        fOutputList->Add(fEtaCutElectronBKULSMainSources_NW);
        fOutputList->Add(fEtaCutElectronBKLSMainSources_NW);
        
        fOutputList->Add(fEtaCutElectronBKNoTag_WithMotherW);
        fOutputList->Add(fEtaCutElectronBKULSMainSources_WithMotherW);
        fOutputList->Add(fEtaCutElectronBKULSMainSources_WithMotherW_NW);
        fOutputList->Add(fEtaCutElectronBKLSMainSources_WithMotherW);
        fOutputList->Add(fEtaCutElectronBKLSMainSources_WithMotherW_NW);
        
        
        fOutputList->Add(fEtaCutElectronInclusiveRecoPtEtaZvtx);
        fOutputList->Add(fEtaCutElectronBKNoTag);
        fOutputList->Add(fEtaCutElectronBKWithLabelLS);
        fOutputList->Add(fEtaCutElectronBKWithLabelULS);
        fOutputList->Add(fEtaCutElectronRecoHFEMC);
        fOutputList->Add(fEtaCutElectronRecoOtherMC);
        fOutputList->Add(fMissIDElectronsReco);
        
        fOutputList->Add(fElectronNoLabel);
        fOutputList->Add(fElectronNoLabelULS);
        fOutputList->Add(fElectronNoLabelLS);
        fOutputList->Add(fEtaCutElectronBKULSMainSources);
        fOutputList->Add(fEtaCutElectronBKLSMainSources);
        fOutputList->Add(fEtaCutElectronBKULSOtherSources);
        fOutputList->Add(fEtaCutElectronBKLSOtherSources);
        fOutputList->Add(fEtaCutElectronHFEULS);
        fOutputList->Add(fEtaCutElectronHFELS);
        fOutputList->Add(fEtaCutElectronMissIDULS);
        fOutputList->Add(fEtaCutElectronMissIDLS);
        
        //Hadrons
        fHadronsReco = new TH3F("fHadronsReco",  "p_{T} x #eta x  ZVtx; p_{T}; #eta; Zvtx", NpTbinsH, pTBinsH, NEtaBins, EtaBins, nZvtxBinsMC, ZVtxBins);
        fHadronsRecoPP = new TH3F("fHadronsRecoPP",  "p_{T} x #eta x  ZVtx; p_{T}; #eta; Zvtx", NpTbinsH, pTBinsH, NEtaBins, EtaBins, nZvtxBinsMC, ZVtxBins);
        fHadronsGenerated = new TH3F("fHadronsGenerated", "p_{T} x #eta x  ZVtx; p_{T}; #eta; Zvtx", NpTbinsH, pTBinsH, NEtaBins, EtaBins, nZvtxBinsMC, ZVtxBins);
        
        fOutputList->Add(fHadronsReco);
        fOutputList->Add(fHadronsRecoPP);
        fOutputList->Add(fHadronsGenerated);
        
        
        //Add MC pT of Pi0, eta
        fPtMCpi0_NoMother = new TH1F("fPtMCpi0_NoMother","#pi^{0} distribution from MC with No Mother;p_{t} (GeV/c);Count",2000,0,100);
        fPtMCpi0_PureHijing = new TH1F("fPtMCpi0_PureHijing","#pi^{0} distribution from MC with no Enh. ;p_{t} (GeV/c);Count",2000,0,100);
        
        fPtMCEta_NoMother = new TH1F("fPtMCEta_NoMother","#eta distribution from MC with No Mother;p_{t} (GeV/c);Count",2000,0,100);
        fPtMCEta_PureHijing = new TH1F("fPtMCEta_PureHijing","#eta distribution from MC with no Enh. ;p_{t} (GeV/c);Count",2000,0,100);
        
        fElectronBKGWToDataTotal = new TH1F("fElectronBKGWToDataTotal","Total Bkg weighted to data;p_{t} (GeV/c);Count",110,0.5,6);
        fElectronBKGWToDataULS = new TH1F("fElectronBKGWToDataULS","Bkg ULS weighted to data;p_{t} (GeV/c);Count",110,0.5,6);
        fElectronBKGWToDataLS = new TH1F("fElectronBKGWToDataLS","Bkg LS weighted to data;p_{t} (GeV/c);Count",110,0.5,6);
        
        fOutputList->Add(fElectronBKGWToDataTotal);
        fOutputList->Add(fElectronBKGWToDataULS);
        fOutputList->Add(fElectronBKGWToDataLS);
        
        fOutputList->Add(fPtMCpi0_NoMother);
        fOutputList->Add(fPtMCpi0_PureHijing);
        
        fOutputList->Add(fPtMCEta_NoMother);
        fOutputList->Add(fPtMCEta_PureHijing);
        
        
    }
    
    
    //fInvMass = new TH1F("fInvMass","",100,0,0.3);
    //fInvMassBack = new TH1F("fInvMassBack","",100,0,0.3);
    
    //Inv Mass in pT Bins
    fInvMassULS = new TH1F *[fpTBins.GetSize()];
    fInvMassLS = new TH1F *[fpTBins.GetSize()];
    
    for (Int_t i = 0 ; i < fpTBins.GetSize()-1 ; i++)
    {
        fInvMassULS[i] = new TH1F(Form("fInvMassULS%d",i), Form("ULS Inv Mass distribution for %1.2f < p_{T}^{e} < %1.2f",fpTBins.At(i),fpTBins.At(i+1)),100,0,0.3);
        fInvMassLS[i] = new TH1F(Form("fInvMassLS%d",i), Form("LS Inv Mass distribution for %1.2f < p_{T}^{e} < %1.2f",fpTBins.At(i),fpTBins.At(i+1)),100,0,0.3);
        fOutputList->Add(fInvMassULS[i]);
        fOutputList->Add(fInvMassLS[i]);
        
    }
    
    fDCA = new TH1F("fDCA","",100,0,1);
    fDCABack = new TH1F("fDCABack","",100,0,1);
    fOpAngle = new TH1F("fOpAngle","",100,0,0.5);
    fOpAngleBack = new TH1F("fOpAngleBack","",100,0,0.5);
    
    if(fCorrelationFlag)
    {
        Int_t NumberBins = fpTBins.GetSize();
        
        fCEtaPhi_Inc = new TH2F *[NumberBins];
        fCEtaPhi_Inc_DiHadron = new TH2F *[NumberBins];
        fCEtaPhi_ULS_Weight = new TH2F *[NumberBins];
        fCEtaPhi_LS_Weight = new TH2F *[NumberBins];
        fCEtaPhi_ULS_NoP_Weight = new TH2F *[NumberBins];
        fCEtaPhi_LS_NoP_Weight = new TH2F *[NumberBins];
        
        fCEtaPhi_Inc_EM = new TH2F *[NumberBins];
        fCEtaPhi_ULS_Weight_EM = new TH2F *[NumberBins];
        fCEtaPhi_LS_Weight_EM = new TH2F *[NumberBins];
        
        fOutputList->Add(fDCA);
        fOutputList->Add(fDCABack);
        fOutputList->Add(fOpAngle);
        fOutputList->Add(fOpAngleBack);
    }
    
    
    
    
    for(Int_t i = 0; i < fpTBins.GetSize()-1; i++)
    {
        fTPC_pt[i] = new TH1F(Form("fTPC_pt%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;TPC Electron N#sigma;Count",fpTBins.At(i),fpTBins.At(i+1)),100,20,200);
        fTPCnsigma_pt[i] = new TH1F(Form("fTPCnsigma_pt%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;TPC Electron N#sigma;Count",fpTBins.At(i),fpTBins.At(i+1)),100,-15,10);
        fTOFTPCnsigma_pt[i] = new TH2F(Form("fTOFTPCnsigma_pt%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c; TOF Electron N#sigma; TPC Electron N#sigma; Count",fpTBins.At(i),fpTBins.At(i+1)), 150,-10,20,180,-12,8  );
        
        fOutputList->Add(fTPC_pt[i]);
        fOutputList->Add(fTPCnsigma_pt[i]);
        fOutputList->Add(fTOFTPCnsigma_pt[i]);
        
        if(fCorrelationFlag)
        {
            
            Int_t BinSize = 40;
            fCEtaPhi_Inc[i] = new TH2F(Form("fCEtaPhi_Inc%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
            fCEtaPhi_Inc_DiHadron[i] = new TH2F(Form("fCEtaPhi_Inc_DiHadron%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
            
            
            fCEtaPhi_ULS_Weight[i] = new TH2F(Form("fCEtaPhi_ULS_Weight%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
            fCEtaPhi_LS_Weight[i] = new TH2F(Form("fCEtaPhi_LS_Weight%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
            fCEtaPhi_ULS_NoP_Weight[i] = new TH2F(Form("fCEtaPhi_ULS_NoP_Weight%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
            fCEtaPhi_LS_NoP_Weight[i] = new TH2F(Form("fCEtaPhi_LS_NoP_Weight%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
            
            
            if (fIsMC)
            {
                
                fCEtaPhiNoEtaCutHFe[i] = new TH2F(Form("fCEtaPhiNoEtaCutHFe%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                
                fOutputList->Add(fCEtaPhiNoEtaCutHFe[i]);
                
                fCEtaPhiCutHFe[i] = new TH2F(Form("fCEtaPhiCutHFe%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                
                fOutputList->Add(fCEtaPhiCutHFe[i]);
                
                fCEtaPhi_Back_MC_Tag[i] = new TH2F(Form("fCEtaPhi_Back_MC_Tag%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                fCEtaPhi_Other_MC_Tag[i] = new TH2F(Form("fCEtaPhi_Other_MC_Tag%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                fCEtaPhi_HFe_MC_Tag[i] = new TH2F(Form("fCEtaPhi_HFe_MC_Tag%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                fCEtaPhi_MC_NoMother_Tag[i] = new TH2F(Form("fCEtaPhi_MC_NoMother%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                
                fOutputList->Add(fCEtaPhi_Back_MC_Tag[i]);
                fOutputList->Add(fCEtaPhi_Other_MC_Tag[i]);
                fOutputList->Add(fCEtaPhi_HFe_MC_Tag[i]);
                fOutputList->Add(fCEtaPhi_MC_NoMother_Tag[i]);
                
                
            }
            
            fOutputList->Add(fCEtaPhi_Inc[i]);
            fOutputList->Add(fCEtaPhi_Inc_DiHadron[i]);
            
            fOutputList->Add(fCEtaPhi_ULS_Weight[i]);
            fOutputList->Add(fCEtaPhi_LS_Weight[i]);
            fOutputList->Add(fCEtaPhi_ULS_NoP_Weight[i]);
            fOutputList->Add(fCEtaPhi_LS_NoP_Weight[i]);
            
            
            
            if(fEventMixingFlag)
            {
                Int_t BinSize = 40;
                fCEtaPhi_Inc_EM[i] = new TH2F(Form("fCEtaPhi_Inc_EM%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                
                fCEtaPhi_ULS_Weight_EM[i] = new TH2F(Form("fCEtaPhi_ULS_Weight_EM%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                fCEtaPhi_LS_Weight_EM[i] = new TH2F(Form("fCEtaPhi_LS_Weight_EM%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                
                fOutputList->Add(fCEtaPhi_Inc_EM[i]);
                fOutputList->Add(fCEtaPhi_ULS_Weight_EM[i]);
                fOutputList->Add(fCEtaPhi_LS_Weight_EM[i]);
            }
        }
    }
    
    //pt integrated
    fTPCnsigma_eta = new TH2F("fTPCnsigma_eta",";Pseudorapidity #eta; TPC signal - <TPC signal>_{elec} (#sigma)",100,-0.9,0.9,100,-15,15);
    fTPCnsigma_phi = new TH2F("fTPCnsigma_phi",";Azimuthal Angle #phi; TPC signal - <TPC signal>_{elec} (#sigma)",100,0,2*TMath::Pi(),100,-15,15);
    
    fOutputList->Add(fTPCnsigma_eta);
    fOutputList->Add(fTPCnsigma_phi);
    
    
    
    
    //__________________________________________________________________
    //Efficiency studies
    if(fIsMC)
    {
        fPtBackgroundBeforeReco = new TH1F("fPtBackgroundBeforeReco",";p_{T} (GeV/c);Count",110,0.5,6);
        fOutputList->Add(fPtBackgroundBeforeReco);
    }
    
    if(!fIspp){
        fCentralityHist = new TH1F("fCentralityHist",";Centrality (%); Count",100,0,100);
        fCentralityHistPass = new TH1F("fCentralityHistPass",";Centrality (%); Count",100,0,100);
        fOutputList->Add(fCentralityHist);
        fOutputList->Add(fCentralityHistPass);
        
    }
    
    for (Int_t i=0; i < fOutputList->GetEntries(); ++i)
    {
        TH1 *h1 = dynamic_cast<TH1*>(fOutputList->At(i));
        if (h1)
        {
            h1->Sumw2();
        }
        else printf("Failed to evaluate Sumw2() \n");
    }
    
    
    //______________________________________________________________________
    //Mixed event analysis -- it was removed from is pp because
    if(fEventMixingFlag && fCorrelationFlag)
    {
        fPoolNevents = new TH1F("fPoolNevents","Event Mixing Statistics; Number of events; Count",1000,0,1000);
        fOutputList->Add(fPoolNevents);
        
        Int_t trackDepth = 100000; // number of objects (tracks) kept per event buffer bin. Once the number of stored objects (tracks) is above that limit, the oldest ones are removed.
        Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
        
        Int_t nCentralityBins  = 15;
        Double_t centralityBins[] = { 0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1 };
        
        Int_t nZvtxBins  = 9;
        Double_t vertexBins[] = {-10, -7, -5, -3, -1, 1, 3, 5, 7, 10.01};
        
        fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, (Double_t*) centralityBins, nZvtxBins, (Double_t*) vertexBins);
    }
    
    //______________________________________________________________________
    
    PostData(1, fOutputList);
    
    ///______________________________________________________________________
}

//______________________________________________________________________
//Main loop
//Called for each event
void AliAnalysisTaskHFEpACorrelation::UserExec(Option_t *)
{
    //Check Event
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    
    if(!(fESD || fAOD))
    {
        printf("ERROR: fESD & fAOD not available\n");
        return;
    }
    
    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    
    if(!fVevent)
    {
        printf("ERROR: fVEvent not available\n");
        return;
    }
    
    //Check Cuts
    if(!fCuts)
    {
        AliError("HFE cuts not available");
        return;
    }
    //Check PID
    if(!fPID->IsInitialized())
    {
        // Initialize PID with the given run number
        AliWarning("PID not initialised, get from Run no");
        
        if(fIsAOD)
        {
            fPID->InitializePID(fAOD->GetRunNumber());
        }
        else
        {
            fPID->InitializePID(fESD->GetRunNumber());
        }
    }
    
    //PID response
    fPidResponse = fInputHandler->GetPIDResponse();
    
    
    //Check PID response
    if(!fPidResponse)
    {
        AliDebug(1, "Using default PID Response");
        fPidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class());
    }
    
    fPID->SetPIDResponse(fPidResponse);
    
    fCFM->SetRecEventInfo(fVevent);
    
    Double_t *fListOfmotherkink = 0;
    Int_t fNumberOfVertices = 0;
    Int_t fNumberOfMotherkink = 0;
    
    
    //total event before event selection
    fNevent->Fill(1);
    
    //Physics Selection + Pileup rejection from Task
    
    if (fAOD->GetRunNumber()>200000 && !fIsMC)
    {
        UInt_t fSelectMask = fInputHandler->IsEventSelected();
        Bool_t isINT7selected = fSelectMask & AliVEvent::kINT7;
        
        if (!isINT7selected){
            fNevent2->Fill(25);
            return;
        }
    }
    
    //______________________________________________________________________
    //Vertex Selection
    if(!fIspp){
        
        fNevent2->Fill(8);
        
        if(fIsAOD)
        {
            const AliAODVertex* trkVtx = fAOD->GetPrimaryVertex();
            
            
            Float_t zvtx = -100;
            zvtx=trkVtx->GetZ();
            fZvtx = zvtx;
            //x
            Float_t xvtx = -100;
            xvtx=trkVtx->GetX();
            
            //y
            Float_t yvtx = -100;
            yvtx=trkVtx->GetY();
            
            
            //events without vertex
            if(zvtx==0 && xvtx==0 && yvtx==0) fNevent2->Fill(0);
            
            
            if(!trkVtx || trkVtx->GetNContributors()<=0) return;
            TString vtxTtl = trkVtx->GetTitle();
            if(!vtxTtl.Contains("VertexerTracks")) return;
            //Float_t zvtx = trkVtx->GetZ();
            
            trkVtx->GetZ();
            fZvtx = zvtx;
            
            
            fVtxZ_new1->Fill(fZvtx);
            
            const AliAODVertex* spdVtx = fAOD->GetPrimaryVertexSPD();
            if(spdVtx->GetNContributors()<=0) return;
            TString vtxTyp = spdVtx->GetTitle();
            Double_t cov[6]={0};
            spdVtx->GetCovarianceMatrix(cov);
            Double_t zRes = TMath::Sqrt(cov[5]);
            
            fzRes1->Fill(zRes);
            //Yvonne e-mail from 12 June 2015 says it has a bug on "vertexer:Z".
            //if(vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return;
            
            //new line:
            if (spdVtx->IsFromVertexerZ() && (zRes>0.25)) return;
            
            fzRes2->Fill(zRes);
            
            fSPD_track_vtx1->Fill(spdVtx->GetZ() - trkVtx->GetZ());
            if(TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return;
            fSPD_track_vtx2->Fill(spdVtx->GetZ() - trkVtx->GetZ());
            
            
            
            if(TMath::Abs(zvtx) > 10){
                fNevent2->Fill(2);
                fVtxZ_new2->Fill(fZvtx);
                return;
            }
            
            //if(fabs(zvtx>10.0))return;
            
            fVtxZ_new3->Fill(fZvtx);
            fNevent2->Fill(4);
            
            //Look for kink mother for AOD
            
            fNumberOfVertices = 0;
            fNumberOfMotherkink = 0;
            
            if(fIsAOD)
            {
                fNumberOfVertices = fAOD->GetNumberOfVertices();
                
                fListOfmotherkink = new Double_t[fNumberOfVertices];
                
                for(Int_t ivertex=0; ivertex < fNumberOfVertices; ivertex++)
                {
                    AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
                    if(!aodvertex) continue;
                    if(aodvertex->GetType()==AliAODVertex::kKink)
                    {
                        AliAODTrack *mother1 = (AliAODTrack *) aodvertex->GetParent();
                        if(!mother1) continue;
                        Int_t idmother = mother1->GetID();
                        fListOfmotherkink[fNumberOfMotherkink] = idmother;
                        fNumberOfMotherkink++;
                    }
                }
            }
        }
        
        else
        {
            
            
            
            /// ESD
            const AliESDVertex* trkVtx = fESD->GetPrimaryVertex();
            if(!trkVtx || trkVtx->GetNContributors()<=0) return;
            TString vtxTtl = trkVtx->GetTitle();
            if(!vtxTtl.Contains("VertexerTracks")) return;
            Float_t zvtx = -100;
            zvtx=trkVtx->GetZ();
            
            
            const AliESDVertex* spdVtx = fESD->GetPrimaryVertexSPD();
            if(spdVtx->GetNContributors()<=0) return;
            TString vtxTyp = spdVtx->GetTitle();
            Double_t cov[6]={0};
            spdVtx->GetCovarianceMatrix(cov);
            Double_t zRes = TMath::Sqrt(cov[5]);
            if(vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return;
            if(TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return;
            if(TMath::Abs(zvtx) > 10) return;
        }
        
        
        
        fNevent2->Fill(12);
        
        //check pA pileup cut as it is done in the hfe package
        if(fAnalysisUtils->IsPileUpEvent(fVevent)){
            fNevent2->Fill(14);
            return;
        }
        fNevent2->Fill(16);
        
    }//close !Ispp flag
    
    //========================== vertex selection for pp
    
    if(fIspp)
    {
        if(fIsAOD)
        {
            const AliAODVertex* trkVtx = fAOD->GetPrimaryVertex();
            if(!trkVtx || trkVtx->GetNContributors()<=0) return;
            Float_t zvtx = -100;
            zvtx=trkVtx->GetZ();
            fZvtx = zvtx;
            //fVtxZ_new1->Fill(fZvtx);
            
            if(TMath::Abs(zvtx) > 10) return;
            //fVtxZ_new2->Fill(fZvtx);
            
            //Look for kink mother for AOD
            
            fNumberOfVertices = 0;
            fNumberOfMotherkink = 0;
            
            if(fIsAOD)
            {
                fNumberOfVertices = fAOD->GetNumberOfVertices();
                
                fListOfmotherkink = new Double_t[fNumberOfVertices];
                
                for(Int_t ivertex=0; ivertex < fNumberOfVertices; ivertex++)
                {
                    AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
                    if(!aodvertex) continue;
                    if(aodvertex->GetType()==AliAODVertex::kKink)
                    {
                        AliAODTrack *mother1 = (AliAODTrack *) aodvertex->GetParent();
                        if(!mother1) continue;
                        Int_t idmother = mother1->GetID();
                        fListOfmotherkink[fNumberOfMotherkink] = idmother;
                        fNumberOfMotherkink++;
                    }
                }
            }
        }
        else
        {
            
            
            
            /// ESD
            const AliESDVertex* trkVtx = fESD->GetPrimaryVertex();
            if(!trkVtx || trkVtx->GetNContributors()<=0) return;
            Float_t zvtx = -100;
            zvtx=trkVtx->GetZ();
            if(TMath::Abs(zvtx) > 10) return;
        }
    }
    //______________________________________________________________________
    //after vertex selection
    fNevent->Fill(10);
    
    //Only events with at least 2 tracks are accepted
    Int_t fNOtrks =  fVevent->GetNumberOfTracks();
    if(fNOtrks<2) return;
    
    fNevent->Fill(11);
    
    if(fIsAOD){
        Int_t fNOtrks2 =  fAOD->GetNumberOfTracks();
        if(fNOtrks2<2) return;
    }
    
    fNevent->Fill(12);
    
    
    //______________________________________________________________________
    //______________________________________________________________________
    //______________________________________________________________________
    //Centrality Selection
    if(!fIspp){
        if(fHasCentralitySelection)
        {
            fCentralityValue = -1;
            //Run 2 (:
            if (fAOD->GetRunNumber()>200000)
            {
                AliMultSelection *MultSelection = 0x0;
                MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
                if(!MultSelection)
                {
                    AliWarning("AliMultSelection object not found!");
                }
                else
                {
                    if (fEstimator==0)
                    {
                        fCentralityValue =  MultSelection->GetMultiplicityPercentile("ZNA");
                    }
                    else
                    {
                        fCentralityValue = MultSelection->GetMultiplicityPercentile("V0A");
                    }
                    
                }
            }
            else
            {
                //Run 1 centrality
                fCentrality = ((AliAODHeader*)fAOD->GetHeader())->GetCentralityP();
                if(fEstimator==0)
                    fCentralityValue = fCentrality->GetCentralityPercentile("ZNA");
                else
                    fCentralityValue = fCentrality->GetCentralityPercentile("V0A");
            }
            
            fCentralityHist->Fill(fCentralityValue);
            
            if(fCentralityValue<fCentralityMin || fCentralityValue>fCentralityMax) return;
            
            fCentralityHistPass->Fill(fCentralityValue);
        }
    }
    //______________________________________________________________________
    
    
    //______________________________________________________________________
    //threshold selection was here
    //______________________________________________________________________
    //all events selected
    
    fNevent->Fill(0);
    
    
    //______________________________________________________________________
    //events in the threshold
    
    fVtxZ_new4->Fill(fZvtx);
    
    
    
    
    if (fIsMC)
    {
        //Inicialzing MC Array
        fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
        fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
        
        Bool_t isNHFe = kFALSE;
        Bool_t isHFe = kFALSE;
        Bool_t isFromOtherLightMeson = kFALSE;
        Bool_t isOther = kFALSE;
        Bool_t HasNoMother = kFALSE;
        Bool_t IsOnDCAcut = kFALSE;
        
        //Vertex from MC
        Double_t vtxMC[3];
        fMCheader->GetVertex(vtxMC);
        
        ComputeWeightInEnhancedSample();
        
        for(Int_t iMC = 0; iMC < fMCarray->GetEntriesFast(); iMC++)
        {
            isNHFe = kFALSE;
            isHFe = kFALSE;
            isOther = kFALSE;
            HasNoMother = kFALSE;
            isFromOtherLightMeson = kFALSE;
            IsOnDCAcut = kFALSE;
            
            
            fMCparticle = (AliAODMCParticle*) fMCarray->At(iMC);
            Int_t ParticlePDG = TMath::Abs(fMCparticle->GetPdgCode());
            Double_t EtaMC = fMCparticle->Eta();
            
            //Check Pi0 and eta spectra in MC
            
            if (fMCparticle->Charge() == 0) continue;
            //Save the pT of all Charged hadrons in the acceptance (This is the denominator of the efficiency)
            if (fMCparticle->IsPhysicalPrimary())
                fHadronsGenerated->Fill(fMCparticle->Pt(),EtaMC, fZvtx);
            
            if(fMCparticle->Pt() < fMinpTElec || fMCparticle->Pt() > fMaxpTElec) continue;
            
            
            Int_t MotherPDG = -999;
            
            //electron correlation
            if (ParticlePDG == 11)
            {
                if (fMCparticle->IsPhysicalPrimary())
                {
                    if(fMCparticle->GetMother()>0)
                    {
                        fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                        MotherPDG = TMath::Abs(fMCparticleMother->GetPdgCode());
                        
                        Int_t MotherPDGHeavy  = Int_t (MotherPDG / TMath::Power(10, Int_t(TMath::Log10(MotherPDG))));
                        
                        if (MotherPDGHeavy > 3)
                        {
                            fNoEtaCutElectronGeneratedSignalPtEtaZvtx->Fill(fMCparticle->Pt());
                            if (EtaMC >= fEtaCutMin && EtaMC <= fEtaCutMax)
                            {
                                fEtaCutElectronGeneratedSignalPtEtaZvtx->Fill(fMCparticle->Pt());
                                isHFe = kTRUE;
                            }
                            
                        }
                    }
                    
                    if (!fCorrelationFlag) continue;
                    
                    if (!isHFe) continue;
                    
                    for (Int_t iHadron = 0 ; iHadron < fMCarray->GetEntriesFast(); iHadron++)
                    {
                        
                        fMCparticle2 = (AliAODMCParticle*) fMCarray->At(iHadron);
                        
                        if (iHadron == iMC) continue;
                        
                        if (!fMCparticle2->IsPhysicalPrimary()) continue;
                        
                        if (fMCparticle2->Charge() == 0) continue;
                        
                        if(fMCparticle2->Pt() <fAssHadronPtMin || fMCparticle2->Pt() > fAssHadronPtMax) continue;
                        Double_t EtaHadron = fMCparticle2->Eta();
                        Double_t dEta =  EtaMC - EtaHadron;
                        Double_t dPhi = fMCparticle->Phi() - fMCparticle2->Phi();
                        
                        if (dPhi > 3*TMath::Pi()/2) dPhi = dPhi - 2*TMath::Pi();
                        if (dPhi < -TMath::Pi()/2)  dPhi = dPhi + 2*TMath::Pi();
                        
                        for (Int_t i = 0 ; i < fpTBins.GetSize()-1; i++)
                        {
                            if(fMCparticle->Pt()>=fpTBins.At(i) && fMCparticle->Pt()<fpTBins.At(i+1))
                            {
                                fCEtaPhiNoEtaCutHFe[i]->Fill(dPhi,dEta);
                                
                                if ( EtaMC>= fEtaCutMin &&  EtaMC <= fEtaCutMax && EtaHadron >= fEtaCutMin &&  EtaHadron <= fEtaCutMax )
                                {
                                    fCEtaPhiCutHFe[i]->Fill(dPhi,dEta);
                                    
                                }
                                
                            }
                        }
                        
                    }
                    
                    
                    
                    
                }
                
            }
            
            
            
        }
    }
    
    //______________________________________________________________________
    
    ///_____________________________________________________________________
    ///Track loop
    Int_t NTracks=0;
    
    NTracks=fVevent->GetNumberOfTracks();
    
    for(Int_t iTracks = 0; iTracks < NTracks; iTracks++)
    {
        //AliVParticle* Vtrack = fVevent->GetTrack(iTracks);
        AliVParticle* Vtrack = 0x0;
        Vtrack  = fVevent->GetTrack(iTracks);
        
        
        if (!Vtrack)
        {
            printf("ERROR: Could not receive track %d\n", iTracks);
            continue;
        }
        
        
        
        AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
        AliESDtrack *etrack = dynamic_cast<AliESDtrack*>(Vtrack);
        AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);
        
        
        Double_t fTPCnSigma0 = -999;
        Double_t fTPCnSigma = -999;
        Double_t fTOFnSigma = -999;
        Double_t fTPCsignal = -999;
        Double_t fPt = -999;
        Double_t fP = -999;
        Double_t EtaTrig = track->Eta();
        
        
        if (fIsMC)
        {
            Bool_t PassAODFilterBitSelection = kFALSE;
            Bool_t IsOnHadronDCAcut = kTRUE;
            
            if (fUseGlobalTracksHadron)
            {
                if (atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA))
                    PassAODFilterBitSelection = kTRUE;
            }
            else
            {
                if (atrack->TestFilterMask(AliAODTrack::kTrkTPCOnly))
                    PassAODFilterBitSelection = kTRUE;
            }
            if(EtaTrig>=fEtaCutMin && EtaTrig<=fEtaCutMax && PassAODFilterBitSelection && atrack->GetStatus()&AliESDtrack::kITSrefit && atrack->GetStatus()&AliESDtrack::kTPCrefit && atrack->GetTPCNcls() >= 80)
            {
                
                if(fUseDCACutforHadrons)
                {
                    Double_t d0z0[2], cov[3];
                    AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
                    track->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 20., d0z0, cov);
                    Double_t DCAxy = d0z0[0];
                    Double_t DCAz = d0z0[1];
                    if(TMath::Abs(DCAxy) >= fDCAcutrHadron || TMath::Abs(DCAz)>=fDCAcutzHadron)
                        IsOnHadronDCAcut = kFALSE;
                }
                
                if (IsOnHadronDCAcut)
                {
                    if (track->GetLabel() >= 0)
                    {
                        fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
                        if (fMCparticle->IsPhysicalPrimary())
                            fHadronsRecoPP->Fill(fMCparticle->Pt(),fMCparticle->Eta(), fZvtx);
                        
                    }
                    
                    fHadronsReco->Fill(track->Pt(),EtaTrig, fZvtx);
                }
            }
        }
        
        
        
        //Etacut test on the begging
        fEtad[0]->Fill(EtaTrig);
        
        ///_____________________________________________________________________________
        ///Fill QA plots without track selection
        fPt = track->Pt();
        fP = TMath::Sqrt((track->Pt())*(track->Pt()) + (track->Pz())*(track->Pz()));
        
        fTPCsignal = track->GetTPCsignal();
        fTPCnSigma = fPidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        fTOFnSigma = fPidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
        
        
        fTPC_momentum->Fill(fP,fTPCsignal);
        fTPC_eta->Fill(EtaTrig,fTPCsignal);
        
        fTPC_p[0]->Fill(fP,fTPCsignal);
        fTPCnsigma_p[0]->Fill(fP,fTPCnSigma);
        Float_t TPCNcls = track->GetTPCNcls();
        //TPC Ncls for pid
        Float_t TPCNcls_pid = track->GetTPCsignalN();
        //______________________________________________________________
        // Vertex
        
        fVtxZ[0]->Fill(fZvtx);
        if(iTracks == 0)fNTracks[0]->Fill(NTracks);
        fNTracks_pt[0]->Fill(NTracks, fPt);
        fNTracks_eta[0]->Fill(NTracks, EtaTrig);
        fNTracks_phi[0]->Fill(NTracks, track->Phi());
        
        
        fTPCNcls_pid[0]->Fill(TPCNcls, TPCNcls_pid);
        //______________________________________________________________
        
        ///Fill QA plots without track selection
        ///_________________________________________________________________________
        //__________________________________________________________________________
        //Track Selection Cuts
        
        if(EtaTrig<fEtaCutMin || EtaTrig>fEtaCutMax) continue;
        
        
        //AOD (Test Filter Bit)
        if(fIsAOD)
        {
            // standard cuts with very loose DCA - BIT(4)
            // Description:
            /*
             GetStandardITSTPCTrackCuts2011(kFALSE)
             SetMaxChi2PerClusterTPC(4);
             SetAcceptKinkDaughters(kFALSE);
             SetRequireTPCRefit(kTRUE);
             SetRequireITSRefit(kTRUE);
             SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
             SetMaxDCAToVertexZ(2);
             SetMaxDCAToVertex2D(kFALSE);
             SetRequireSigmaToVertex(kFALSE);
             SetMaxChi2PerClusterITS(36);
             SetMaxDCAToVertexXY(2.4)
             SetMaxDCAToVertexZ(3.2)
             SetDCaToVertex2D(kTRUE)
             */
            
            if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;
        }
        
        
        
        
        //RecKine: ITSTPC cuts
        if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
        
        //RecKink
        if(fRejectKinkMother)
        {
            if(fIsAOD)
            {
                Bool_t kinkmotherpass = kTRUE;
                for(Int_t kinkmother = 0; kinkmother < fNumberOfMotherkink; kinkmother++)
                {
                    if(track->GetID() == fListOfmotherkink[kinkmother])
                    {
                        kinkmotherpass = kFALSE;
                        continue;
                    }
                }
                if(!kinkmotherpass) continue;
            }
            else
            {
                if(etrack->GetKinkIndex(0) != 0) continue;
            }
        }
        
        
        if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;
        
        //HFEcuts: ITS layers cuts
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;
        
        //HFE cuts: TPC PID cleanup
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;
        
        fEtad[1]->Fill(EtaTrig);
        fTPC_p[1]->Fill(fP,fTPCsignal);
        fTPCnsigma_p[1]->Fill(fP,fTPCnSigma);
        TPCNcls = track->GetTPCNcls();
        
        //______________________________________________________________
        // Vertex
        
        fVtxZ[1]->Fill(fZvtx);
        if(iTracks == 0)fNTracks[1]->Fill(NTracks);
        fNTracks_pt[1]->Fill(NTracks, fPt);
        fNTracks_eta[1]->Fill(NTracks, EtaTrig);
        fNTracks_phi[1]->Fill(NTracks, track->Phi());
        fTPCNcls_pid[1]->Fill(TPCNcls, TPCNcls_pid);
        
        //______________________________________________________________
        
        ///______________________________________________________________________
        ///Histograms for PID Studies
        
        for(Int_t i = 0; i < fpTBins.GetSize()-1; i++)
        {
            if(fPt>=fpTBins.At(i) && fPt<fpTBins.At(i+1))
            {
                fTPC_pt[i]->Fill(fTPCsignal);
                fTPCnsigma_pt[i]->Fill(fTPCnSigma);
                fTOFTPCnsigma_pt[i]->Fill(fTOFnSigma,fTPCnSigma);
                
            }
        }
        
        //_______________________________________________________
        //Correlation Analysis - DiHadron
        if(fTPCnSigma < 3.5 && fCorrelationFlag)
        {
            fPtTrigger_Inc->Fill(fPt);
            DiHadronCorrelation(track, iTracks);
        }
        
        //Add TOF-Only PID cuts to study contamination as function of p
        if (TMath::Abs(fTOFnSigma) <= 3)
        {
            fTPC_p[2]->Fill(fP,fTPCsignal);
            fTPCnsigma_p[2]->Fill(fP,fTPCnSigma);
        }
        
        
        ///________________________________________________________________________
        ///PID
        ///Here the PID cuts defined in the file "ConfigEMCalHFEpA.C" is applied
        
        
        Int_t pidpassed = 1;
        AliHFEpidObject hfetrack;
        hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
        hfetrack.SetRecTrack(track);
        hfetrack.SetPP();	//proton-proton analysis
        if(!fPID->IsSelected(&hfetrack, NULL, "", fPIDqa)) pidpassed = 0;
        fpid->Fill(pidpassed);
        
        if(pidpassed==0) continue;
        
        
        //______________________________________________________________
        // Vertex
        
        fVtxZ[2]->Fill(fZvtx);
        if(iTracks == 0)fNTracks[2]->Fill(NTracks);
        fNTracks_pt[2]->Fill(NTracks, fPt);
        fNTracks_eta[2]->Fill(NTracks, EtaTrig);
        fNTracks_phi[2]->Fill(NTracks, track->Phi());
        fTPCNcls_pid[2]->Fill(TPCNcls, TPCNcls_pid);
        
        //______________________________________________________________
        
        //_______________________________________________________
        //Correlation Analysis
        
        
        fPtElec_Inc->Fill(fPt);
        
        Double_t d0z0[2], cov[3];
        
        AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
        Double_t DCAxy = 100;
        Double_t DCAz = 100;
        if (track->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 20., d0z0, cov))
        {
            DCAxy = TMath::Abs(d0z0[0]);
            DCAz =  TMath::Abs(d0z0[1]);
        }
        
        fDCAElectronXY[0]->Fill(DCAxy);
        fDCAElectronZ[0]->Fill(DCAz);
        
        fEtad[2]->Fill(EtaTrig);
        
        if(track->Pt()< fMinpTElec || track->Pt() > fMaxpTElec) continue;
        
        ElectronHadronCorrelation(track, iTracks, Vtrack);
        
        //_______________________________________________________
        
        ///________________________________________________________________________
    }
    
    //__________________________________________________________________
    //Event Mixing Analysis
    //Filling pool
    if(fEventMixingFlag && fCorrelationFlag)
    {
        if(fIspp)
        {
            fPool = fPoolMgr->GetEventPool(1.5, fZvtx); // Get the buffer associated with the current centrality and z-vtx
            
            //if(!fPool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fCentrality->GetCentralityPercentile("V0A"), fZvtx));
            if(!fPool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f", 1.5, fZvtx));
            
            TObjArray *fArrayTracksMixed = SelectedHadrons();
            fArrayTracksMixed->SetOwner(kTRUE);
            
            fPool->UpdatePool(fArrayTracksMixed);
        }
        else
        {
            fPool = fPoolMgr->GetEventPool(fCentralityValue, fZvtx); // Get the buffer associated with the current centrality and z-vtx
            
            if(!fPool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fCentralityValue, fZvtx));
            
            fPool->UpdatePool(SelectedHadrons()); // fill the tracks in the event buffer. The ownership is handed over to the event mixing class. We are not allowed to delete tracksClone anymore!
        }
        
    }
    
    //__________________________________________________________________
    
    delete fListOfmotherkink;
    PostData(1, fOutputList);
}


//______________________________________________________________________
void AliAnalysisTaskHFEpACorrelation::Terminate(Option_t *)
{
    //Draw result to the screen
    //Called once at the end of the query
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    
    if(!fOutputList)
    {
        printf("ERROR: Output list not available\n");
        return;
    }
}

//______________________________________________________________________
Bool_t AliAnalysisTaskHFEpACorrelation::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
    //Check single track cuts for a given cut step
    //Note this function is called inside the UserExec function
    const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
    if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
    return kTRUE;
}




//______________________________________________________________________

//______________________________________________________________________
void AliAnalysisTaskHFEpACorrelation::ElectronHadronCorrelation(AliVTrack *track, Int_t trackIndex, AliVParticle *vtrack)
{
    
    ///_________________________________________________________________
    ///MC analysis
    Bool_t lIsNHFe = kFALSE;
    Bool_t lIsHFe = kFALSE;
    Bool_t lIsOther = kFALSE;
    Bool_t lHasMother = kFALSE;
    
    //Electron Information
    Double_t fPhiE = -999;
    Double_t fEtaE = -999;
    Double_t fPhiH = -999;
    Double_t fEtaH = -999;
    Double_t fDphi = -999;
    Double_t fDeta = -999;
    Double_t fPtE = -999;
    Double_t fPtH = -999;
    
    Double_t pi = TMath::Pi();
    
    fPhiE = track->Phi();
    fEtaE = track->Eta();
    fPtE = track->Pt();
    
    
    ///_________________________________________________________________
    
    //________________________________________________
    //Associated particle cut
    fPartnerCuts->SetAcceptKinkDaughters(kFALSE);
    fPartnerCuts->SetRequireITSRefit(kTRUE);
    fPartnerCuts->SetRequireTPCRefit(kTRUE);
    fPartnerCuts->SetEtaRange(fEtaCutMin,fEtaCutMax);
    fPartnerCuts->SetMaxChi2PerClusterTPC(4.0);
    fPartnerCuts->SetMinNClustersTPC(80);
    fPartnerCuts->SetPtRange(fPtMinAsso,1e10);
    //fPartnerCuts->SetRequireSigmaToVertex(kTRUE);
    //fPartnerCuts->SetMaxDCAToVertexXY(1);
    //fPartnerCuts->SetMaxDCAToVertexZ(3);
    //_________________________________________________
    
    ///#################################################################
    //Non-HFE reconstruction
    //fNonHFE = new AliSelectNonHFE(); I do not need to create if again
    fNonHFE->SetAODanalysis(fIsAOD);
    
    if(fMassCutFlag) fNonHFE->SetInvariantMassCut(fMassCut);
    if(fAngleCutFlag) fNonHFE->SetOpeningAngleCut(fAngleCut);
    if(fChi2CutFlag) fNonHFE->SetChi2OverNDFCut(fChi2Cut);
    if(fDCAcutFlag) fNonHFE->SetDCACut(fDCAcut);
    fNonHFE->SetAlgorithm("DCA"); //KF
    
    //Remove strong cuts
    //fNonHFE->SetUseGlobalTracks();
    //fNonHFE->SetNClustITS(2); Remove ITS cut from the partner. It generates aditional problems for NonID Background subtraction.
    //fNonHFE->SetDCAPartnerCuts(fDCAcutr, fDCAcutz);
    //fNonHFE->SetEtaCutForPart(fEtaCutMin, fEtaCutMax);
    //fNonHFE->SetTPCNclsForPID(60);
    fNonHFE->SetUseITSTPCRefit(kFALSE);
    
    
    fNonHFE->SetPIDresponse(fPidResponse);
    fNonHFE->SetTrackCuts(-3.0,3.0,fPartnerCuts);
    fNonHFE->SetAdditionalCuts(fPtMinAsso,fTpcNclsAsso);
    
    fNonHFE->SetHistAngleBack(fOpAngleBack);
    fNonHFE->SetHistAngle(fOpAngle);
    fNonHFE->SetHistDCABack(fDCABack);
    fNonHFE->SetHistDCA(fDCA);
    
    //"SetHistMassBack" sets the LS histogram in the invariant mass
    //"SetHistMass" sets the ULS histogram in the invariant mass
    
    for (Int_t pTbin = 0; pTbin < fpTBins.GetSize()-1; pTbin++ )
    {
        if(fPtE>=fpTBins.At(pTbin) && fPtE<fpTBins.At(pTbin+1))
        {
            fNonHFE->SetHistMass(fInvMassULS[pTbin]);
            fNonHFE->SetHistMassBack(fInvMassLS[pTbin]);
        }
    }
    
    fNonHFE->FindNonHFE(trackIndex,vtrack,fVevent);
    
    Int_t *fUlsPartner = fNonHFE->GetPartnersULS();
    Int_t *fLsPartner = fNonHFE->GetPartnersLS();
    Bool_t fUlsIsPartner = kFALSE;
    Bool_t fLsIsPartner = kFALSE;
    ///#################################################################
    
    if (fIsMC)
        TaggingEfficiencyCalculation(track,&lIsNHFe,&lIsHFe,&lIsOther,&lHasMother);
    
    
    if(fNonHFE->IsULS()) fPtElec_ULS->Fill(fPtE,fNonHFE->GetNULS());
    if(fNonHFE->IsLS()) fPtElec_LS->Fill(fPtE,fNonHFE->GetNLS());
    
    if (!fCorrelationFlag) return;
    
    if(fEventMixingFlag)
    {
        //hadling pp in the same task
        if (fIspp)
        {
            fPool = fPoolMgr->GetEventPool(1.5, fZvtx);
            
            if(!fPool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f",1.5, fZvtx));
        }
        else
        {
            fPool = fPoolMgr->GetEventPool(fCentralityValue, fZvtx); // Get the buffer associated with the current centrality and z-vtx
            
            if(!fPool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f",fCentralityValue, fZvtx));
        }
        if(fPool->GetCurrentNEvents() >= 5) // start mixing when 5 events are in the buffer
        {
            fPoolNevents->Fill(fPool->GetCurrentNEvents());
            
            for (Int_t jMix = 0; jMix < fPool->GetCurrentNEvents(); jMix++)  // mix with each event in the buffer
            {
                TObjArray* bgTracks = fPool->GetEvent(jMix);
                
                for (Int_t kMix = 0; kMix < bgTracks->GetEntriesFast(); kMix++)  // mix with each track in the event
                {
                    const AliHFEHCParticle* MixedTrack(dynamic_cast<AliHFEHCParticle*>(bgTracks->At(kMix)));
                    if (NULL == MixedTrack) continue;
                    
                    fPhiH = MixedTrack->Phi();
                    fEtaH = MixedTrack->Eta();
                    fPtH = MixedTrack->Pt();
                    
                    if(fPtH<fAssHadronPtMin || fPtH>fAssHadronPtMax) continue;
                    
                    fDphi = fPhiE - fPhiH;
                    
                    if (fDphi > 3*pi/2) fDphi = fDphi - 2*pi;
                    if (fDphi < -pi/2)  fDphi = fDphi + 2*pi;
                    
                    fDeta = fEtaE - fEtaH;
                    
                    Double_t WeightHadron = 1./GetHadronEfficiency(fPtH,fEtaH,fZvtx);
                    
                    for(Int_t i = 0; i < fpTBins.GetSize()-1; i++)
                    {
                        if(fPtE>=fpTBins.At(i) && fPtE<fpTBins.At(i+1))
                        {
                            fCEtaPhi_Inc_EM[i]->Fill(fDphi,fDeta,WeightHadron);
                            
                            if(fNonHFE->IsULS()) fCEtaPhi_ULS_Weight_EM[i]->Fill(fDphi,fDeta,fNonHFE->GetNULS()*WeightHadron);
                            if(fNonHFE->IsLS()) fCEtaPhi_LS_Weight_EM[i]->Fill(fDphi,fDeta,fNonHFE->GetNLS()*WeightHadron);
                            
                            
                        }
                    }
                }
            }
        }
        
    }
    //__________________________________________________________________
    
    //__________________________________________________________________
    //Same Event Analysis - Hadron Loop
    for(Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++)
    {
        if(trackIndex==iTracks) continue;
        
        AliVParticle* Vtrack2 = fVevent->GetTrack(iTracks);
        if (!Vtrack2)
        {
            printf("ERROR: Could not receive track %d\n", iTracks);
            continue;
        }
        
        AliVTrack *track2 = dynamic_cast<AliVTrack*>(Vtrack2);
        
        fPhiH = track2->Phi();
        fEtaH = track2->Eta();
        fPtH = track2->Pt();
        
        if(fEtaH<fEtaCutMin || fEtaH>fEtaCutMax ) continue;
        
        if(fIsAOD)
        {
            AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(Vtrack2);
            
            if (fUseGlobalTracksHadron)
            {
                if (!atrack2->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA))
                {
                    continue;
                }
            }
            else
                if(!atrack2->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
            
            if((!(atrack2->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack2->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
            if(atrack2->GetTPCNcls() < 80) continue;
            if(fAssocWithSPD && ((!(atrack2->HasPointOnITSLayer(0))) && (!(atrack2->HasPointOnITSLayer(1))))) continue;
        }
        else
        {
            AliESDtrack *etrack2 = dynamic_cast<AliESDtrack*>(Vtrack2);
            if(!fPartnerCuts->AcceptTrack(etrack2)) continue;
        }
        
        //DCA cut for hadrons
        if(fIsAOD && fUseDCACutforHadrons)
        {
            Double_t d0z0[2], cov[3];
            AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
            track2->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 20., d0z0, cov);
            Double_t DCAxy = d0z0[0];
            Double_t DCAz = d0z0[1];
            if(TMath::Abs(DCAxy) >= fDCAcutrHadron || TMath::Abs(DCAz)>=fDCAcutzHadron) continue;
        }
        
        
        
        if(fPtH<fAssHadronPtMin || fPtH>fAssHadronPtMax) continue;
        
        fDphi = fPhiE - fPhiH;
        
        if (fDphi > 3*pi/2) fDphi = fDphi - 2*pi;
        if (fDphi < -pi/2)  fDphi = fDphi + 2*pi;
        
        fDeta = fEtaE - fEtaH;
        
        
        
        //______________________________________________________________
        //Check if this track is a Non-HFE partner
        fUlsIsPartner = kFALSE;
        fLsIsPartner = kFALSE;
        for(Int_t i = 0; i < fNonHFE->GetNULS(); i++)
        {
            if(fUlsPartner[i]==iTracks) fUlsIsPartner=kTRUE;
        }
        for(Int_t i = 0; i < fNonHFE->GetNLS(); i++)
        {
            if(fLsPartner[i]==iTracks) fLsIsPartner=kTRUE;
        }
        //______________________________________________________________
        
        //Double_t WeightInclusive = 1./(GetInclusiveEfficiency(fPtH,fEtaH,fZvtx) * GetHadronEfficiency(fPtH,fEtaH,fZvtx));
        //Double_t WeightBKG = 1./(GetBackgroundEfficiency(fPtH,fEtaH,fZvtx) * GetHadronEfficiency(fPtH,fEtaH,fZvtx));
        
        Double_t WeightHadron = 1./GetHadronEfficiency(fPtH,fEtaH,fZvtx);
        
        for(Int_t i = 0; i < fpTBins.GetSize()-1; i++)
        {
            if(fPtE>=fpTBins.At(i) && fPtE<fpTBins.At(i+1))
            {
                
                //Filling histograms: Only hadron weight for now.
                fCEtaPhi_Inc[i]->Fill(fDphi,fDeta,WeightHadron);
                
                if(fNonHFE->IsULS())
                    fCEtaPhi_ULS_Weight[i]->Fill(fDphi,fDeta,fNonHFE->GetNULS()*WeightHadron);
                
                if(fNonHFE->IsLS())
                    fCEtaPhi_LS_Weight[i]->Fill(fDphi,fDeta,fNonHFE->GetNLS()*WeightHadron);
                
                if(fNonHFE->IsULS() && !fUlsIsPartner && !fLsIsPartner)
                    fCEtaPhi_ULS_NoP_Weight[i]->Fill(fDphi,fDeta,fNonHFE->GetNULS()*WeightHadron);
                
                if(fNonHFE->IsLS() && !fUlsIsPartner && !fLsIsPartner)
                    fCEtaPhi_LS_NoP_Weight[i]->Fill(fDphi,fDeta,fNonHFE->GetNLS()*WeightHadron);
                
                
                if (fIsMC)
                {
                    if (lHasMother)
                    {
                        if(lIsNHFe)
                            fCEtaPhi_Back_MC_Tag[i]->Fill(fDphi,fDeta,WeightHadron);
                        if(lIsOther)
                            fCEtaPhi_Other_MC_Tag[i]->Fill(fDphi,fDeta,WeightHadron);
                        if(lIsHFe)
                            fCEtaPhi_HFe_MC_Tag[i]->Fill(fDphi,fDeta,WeightHadron);
                        
                    }
                    else
                        fCEtaPhi_MC_NoMother_Tag[i]->Fill(fDphi,fDeta,WeightHadron);
                    
                }
                
            }
        }
    }
    
    
}

//____________________________________________________________________________________________________________
//Create a TObjArray with selected hadrons, for the mixed event analysis
TObjArray* AliAnalysisTaskHFEpACorrelation::SelectedHadrons()
{
    fTracksClone = new TObjArray;
    fTracksClone->SetOwner(kTRUE);
    
    for(Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++)
    {
        AliVParticle* Vtrack2 = fVevent->GetTrack(iTracks);
        if (!Vtrack2)
        {
            printf("ERROR: Could not receive track %d\n", iTracks);
            continue;
        }
        
        
        AliVTrack *track2 = dynamic_cast<AliVTrack*>(Vtrack2);
        
        if(track2->Pt() < fAssHadronPtMin || track2->Pt() > fAssHadronPtMax) continue;
        
        Double_t HadronEta = track2->Eta();
        
        if(HadronEta<fEtaCutMin || HadronEta>fEtaCutMax ) continue;
        
        if(fIsAOD)
        {
            AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(Vtrack2);
            
            if (fUseGlobalTracksHadron)
            {
                if (!atrack2->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA))
                {
                    continue;
                }
            }
            else
                if(!atrack2->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
            
            if((!(atrack2->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack2->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
            if(atrack2->GetTPCNcls() < 80) continue;
            if(fAssocWithSPD && ((!(atrack2->HasPointOnITSLayer(0))) && (!(atrack2->HasPointOnITSLayer(1))))) continue;
        }
        else
        {
            AliESDtrack *etrack2 = dynamic_cast<AliESDtrack*>(Vtrack2);
            if(!fPartnerCuts->AcceptTrack(etrack2)) continue;
        }
        
        //DCA cut for hadrons
        if(fIsAOD && fUseDCACutforHadrons)
        {
            Double_t d0z0[2], cov[3];
            AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
            track2->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 20., d0z0, cov);
            Double_t DCAxy = d0z0[0];
            Double_t DCAz = d0z0[1];
            if(TMath::Abs(DCAxy) >= fDCAcutrHadron || TMath::Abs(DCAz)>=fDCAcutzHadron) continue;
        }
        
        fTracksClone->Add(new AliHFEHCParticle(HadronEta, track2->Phi(), track2->Pt(), 1./(GetHadronEfficiency(track2->Pt(),HadronEta, fZvtx)) ) );
    }
    return fTracksClone;
}
//____________________________________________________________________________________________________________

//______________________________________________________________________
void AliAnalysisTaskHFEpACorrelation::DiHadronCorrelation(AliVTrack *track, Int_t trackIndex)
{
    //________________________________________________
    //Associated particle cut
    fPartnerCuts->SetAcceptKinkDaughters(kFALSE);
    fPartnerCuts->SetRequireITSRefit(kTRUE);
    fPartnerCuts->SetRequireTPCRefit(kTRUE);
    fPartnerCuts->SetEtaRange(-0.9,0.9);
    fPartnerCuts->SetMaxChi2PerClusterTPC(4.0);
    fPartnerCuts->SetMinNClustersTPC(80);
    fPartnerCuts->SetPtRange(0.3,1e10);
    //fPartnerCuts->SetRequireSigmaToVertex(kTRUE);
    //fPartnerCuts->SetMaxDCAToVertexXY(1);
    //fPartnerCuts->SetMaxDCAToVertexZ(3);
    //_________________________________________________
    
    //Electron Information
    Double_t fPhiE = -999;
    Double_t fEtaE = -999;
    Double_t fPhiH = -999;
    Double_t fEtaH = -999;
    Double_t fDphi = -999;
    Double_t fDeta = -999;
    Double_t fPtE = -999;
    Double_t fPtH = -999;
    
    Double_t pi = TMath::Pi();
    
    fPhiE = track->Phi();
    fEtaE = track->Eta();
    fPtE = track->Pt();
    
    //__________________________________________________________________
    //Same Event Analysis - Hadron Loop
    for(Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++)
    {
        if(trackIndex==iTracks) continue;
        
        AliVParticle* Vtrack2 = fVevent->GetTrack(iTracks);
        if (!Vtrack2)
        {
            printf("ERROR: Could not receive track %d\n", iTracks);
            continue;
        }
        
        AliVTrack *track2 = dynamic_cast<AliVTrack*>(Vtrack2);
        
        fPhiH = track2->Phi();
        fEtaH = track2->Eta();
        fPtH = track2->Pt();
        
        if(fEtaH<fEtaCutMin || fEtaH>fEtaCutMax ) continue;
        
        
        
        if(fIsAOD)
        {
            AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(Vtrack2);
            
            if (fUseGlobalTracksHadron)
            {
                if (!atrack2->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA))
                {
                    continue;
                }
            }
            else
                if(!atrack2->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
            
            if((!(atrack2->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack2->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
            if(atrack2->GetTPCNcls() < 80) continue;
            if(fAssocWithSPD && ((!(atrack2->HasPointOnITSLayer(0))) && (!(atrack2->HasPointOnITSLayer(1))))) continue;
        }
        else
        {
            AliESDtrack *etrack2 = dynamic_cast<AliESDtrack*>(Vtrack2);
            if(!fPartnerCuts->AcceptTrack(etrack2)) continue;
        }
        
        //DCA cut for hadrons
        if(fIsAOD && fUseDCACutforHadrons)
        {
            Double_t d0z0[2], cov[3];
            AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
            track2->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 20., d0z0, cov);
            Double_t DCAxy = d0z0[0];
            Double_t DCAz = d0z0[1];
            if(TMath::Abs(DCAxy) >= fDCAcutrHadron || TMath::Abs(DCAz)>=fDCAcutzHadron) continue;
        }
        
        
        
        
        if(fPtH<fAssHadronPtMin || fPtH>fAssHadronPtMax) continue;
        
        fDphi = fPhiE - fPhiH;
        
        if (fDphi > 3*pi/2) fDphi = fDphi - 2*pi;
        if (fDphi < -pi/2)  fDphi = fDphi + 2*pi;
        
        fDeta = fEtaE - fEtaH;
        
        
        for(Int_t i = 0; i<fpTBins.GetSize()-1; i++)
        {
            if(fPtE>=fpTBins.At(i) && fPtE<fpTBins.At(i+1))
            {
                fCEtaPhi_Inc_DiHadron[i]->Fill(fDphi,fDeta);
            }
        }
    }
}
//____________________________________________________________________________________________________________

//______________________________________________________________________


Double_t AliAnalysisTaskHFEpACorrelation::GetHadronEfficiency(Double_t pT, Double_t eta, Double_t zvtx)
{
    if (!fEffHadron)
        return 1.;
    Int_t bin = fEffHadron->FindBin(pT,eta,zvtx);
    if ( fEffHadron->IsBinUnderflow(bin) || fEffHadron->IsBinOverflow(bin) ) {return 1.;}
    //printf ("Setting hadron Eff as %f\n",fEffHadron->GetBinContent(bin));
    return fEffHadron->GetBinContent(bin);
    
}

Double_t AliAnalysisTaskHFEpACorrelation::CalculateWeight(Int_t pdg_particle, Double_t x)
{
    //weight for d3 based on MinJung parametrization //sent by Jan
    Double_t weight=1.;
    //MB
    
    if(pdg_particle==111){
        if(x>= 0.100000 &&  x < 0.112797 ) weight=1.262120;
        if(x>= 0.112797 &&  x < 0.127231 ) weight=1.277765;
        if(x>= 0.127231 &&  x < 0.143512 ) weight=1.295605;
        if(x>= 0.143512 &&  x < 0.161877 ) weight=1.318155;
        if(x>= 0.161877 &&  x < 0.182592 ) weight=1.348693;
        if(x>= 0.182592 &&  x < 0.205957 ) weight=1.388636;
        if(x>= 0.205957 &&  x < 0.232313 ) weight=1.439122;
        if(x>= 0.232313 &&  x < 0.262041 ) weight=1.497452;
        if(x>= 0.262041 &&  x < 0.295573 ) weight=1.559409;
        if(x>= 0.295573 &&  x < 0.333397 ) weight=1.615169;
        if(x>= 0.333397 &&  x < 0.376060 ) weight=1.654954;
        if(x>= 0.376060 &&  x < 0.424183 ) weight=1.668753;
        if(x>= 0.424183 &&  x < 0.478465 ) weight=1.652225;
        if(x>= 0.478465 &&  x < 0.539692 ) weight=1.603119;
        if(x>= 0.539692 &&  x < 0.608754 ) weight=1.526049;
        if(x>= 0.608754 &&  x < 0.686654 ) weight=1.426724;
        if(x>= 0.686654 &&  x < 0.774523 ) weight=1.312684;
        if(x>= 0.774523 &&  x < 0.873636 ) weight=1.195395;
        if(x>= 0.873636 &&  x < 0.985432 ) weight=1.086264;
        if(x>= 0.985432 &&  x < 1.111534 ) weight=0.993666;
        if(x>= 1.111534 &&  x < 1.253773 ) weight=0.922587;
        if(x>= 1.253773 &&  x < 1.414214 ) weight=0.875739;
        if(x>= 1.414214 &&  x < 1.595185 ) weight=0.852181;
        if(x>= 1.595185 &&  x < 1.799315 ) weight=0.847828;
        if(x>= 1.799315 &&  x < 2.029567 ) weight=0.863875;
        if(x>= 2.029567 &&  x < 2.289283 ) weight=0.899112;
        if(x>= 2.289283 &&  x < 2.582235 ) weight=0.955194;
        if(x>= 2.582235 &&  x < 2.912674 ) weight=1.033824;
        if(x>= 2.912674 &&  x < 3.285398 ) weight=1.133714;
        if(x>= 3.285398 &&  x < 3.705818 ) weight=1.259471;
        if(x>= 3.705818 &&  x < 4.180038 ) weight=1.406883;
        if(x>= 4.180038 &&  x < 4.714942 ) weight=1.578923;
        if(x>= 4.714942 &&  x < 5.318296 ) weight=1.778513;
        if(x>= 5.318296 &&  x < 5.998859 ) weight=2.001171;
        if(x>= 5.998859 &&  x < 6.766511 ) weight=2.223161;
        if(x>= 6.766511 &&  x < 7.632396 ) weight=2.449445;
        if(x>= 7.632396 &&  x < 8.609086 ) weight=2.661734;
        if(x>= 8.609086 &&  x < 9.710759 ) weight=2.851935;
        if(x>= 9.710759 &&  x < 10.953409 ) weight=2.974319;
        if(x>= 10.953409 &&  x < 12.355077 ) weight=3.106314;
        if(x>= 12.355077 &&  x < 13.936111 ) weight=3.126815;
        if(x>= 13.936111 &&  x < 15.719464 ) weight=3.150053;
        if(x>= 15.719464 &&  x < 17.731026 ) weight=3.218509;
        if(x>= 17.731026 &&  x < 20.000000 ) weight=3.252141;
        
    }
    else if(pdg_particle==221){
        if(x>= 0.100000 &&  x < 0.112797 ) weight=2.159358;
        if(x>= 0.112797 &&  x < 0.127231 ) weight=2.145546;
        if(x>= 0.127231 &&  x < 0.143512 ) weight=2.132390;
        if(x>= 0.143512 &&  x < 0.161877 ) weight=2.109918;
        if(x>= 0.161877 &&  x < 0.182592 ) weight=2.084920;
        if(x>= 0.182592 &&  x < 0.205957 ) weight=2.054302;
        if(x>= 0.205957 &&  x < 0.232313 ) weight=2.015202;
        if(x>= 0.232313 &&  x < 0.262041 ) weight=1.966068;
        if(x>= 0.262041 &&  x < 0.295573 ) weight=1.912255;
        if(x>= 0.295573 &&  x < 0.333397 ) weight=1.844087;
        if(x>= 0.333397 &&  x < 0.376060 ) weight=1.767913;
        if(x>= 0.376060 &&  x < 0.424183 ) weight=1.680366;
        if(x>= 0.424183 &&  x < 0.478465 ) weight=1.583468;
        if(x>= 0.478465 &&  x < 0.539692 ) weight=1.475131;
        if(x>= 0.539692 &&  x < 0.608754 ) weight=1.361436;
        if(x>= 0.608754 &&  x < 0.686654 ) weight=1.244388;
        if(x>= 0.686654 &&  x < 0.774523 ) weight=1.125973;
        if(x>= 0.774523 &&  x < 0.873636 ) weight=1.015769;
        if(x>= 0.873636 &&  x < 0.985432 ) weight=0.914579;
        if(x>= 0.985432 &&  x < 1.111534 ) weight=0.830529;
        if(x>= 1.111534 &&  x < 1.253773 ) weight=0.766397;
        if(x>= 1.253773 &&  x < 1.414214 ) weight=0.723663;
        if(x>= 1.414214 &&  x < 1.595185 ) weight=0.701236;
        if(x>= 1.595185 &&  x < 1.799315 ) weight=0.695605;
        if(x>= 1.799315 &&  x < 2.029567 ) weight=0.707578;
        if(x>= 2.029567 &&  x < 2.289283 ) weight=0.735194;
        if(x>= 2.289283 &&  x < 2.582235 ) weight=0.781052;
        if(x>= 2.582235 &&  x < 2.912674 ) weight=0.842350;
        if(x>= 2.912674 &&  x < 3.285398 ) weight=0.923676;
        if(x>= 3.285398 &&  x < 3.705818 ) weight=1.028317;
        if(x>= 3.705818 &&  x < 4.180038 ) weight=1.154029;
        if(x>= 4.180038 &&  x < 4.714942 ) weight=1.296915;
        if(x>= 4.714942 &&  x < 5.318296 ) weight=1.463674;
        if(x>= 5.318296 &&  x < 5.998859 ) weight=1.651985;
        if(x>= 5.998859 &&  x < 6.766511 ) weight=1.847912;
        if(x>= 6.766511 &&  x < 7.632396 ) weight=2.066284;
        if(x>= 7.632396 &&  x < 8.609086 ) weight=2.202231;
        if(x>= 8.609086 &&  x < 9.710759 ) weight=2.399942;
        if(x>= 9.710759 &&  x < 10.953409 ) weight=2.555106;
        if(x>= 10.953409 &&  x < 12.355077 ) weight=2.632377;
        if(x>= 12.355077 &&  x < 13.936111 ) weight=2.609660;
        if(x>= 13.936111 &&  x < 15.719464 ) weight=2.656343;
        if(x>= 15.719464 &&  x < 17.731026 ) weight=2.615438;
        if(x>= 17.731026 &&  x < 20.000000 ) weight=2.576269;
        
    }
    else weight=1.;
    
    
    return weight/40000.;
    
}

Double_t AliAnalysisTaskHFEpACorrelation::CalculateWeightRun2(Int_t pdg_particle, Double_t pT)
{
    if (TMath::Abs(pdg_particle) == 111)
    {
        if (!fBkgPi0Weight)
            return 1.0;
        
        Int_t bin = fBkgPi0Weight->FindBin(pT);
        return fBkgPi0Weight->GetBinContent(bin);
    }
    else if (TMath::Abs(pdg_particle) == 221)
    {
        if (!fBkgEtaWeight)
            return 1.0;
        
        Int_t bin = fBkgEtaWeight->FindBin(pT);
        return fBkgEtaWeight->GetBinContent(bin);
    }
    
    return 1.0;
    
}

void AliAnalysisTaskHFEpACorrelation::ComputeWeightInEnhancedSample()
{
    for (Int_t MCIndex = 0; MCIndex < fMCarray->GetEntriesFast() ; MCIndex++)
    {
        AliAODMCParticle* particle = (AliAODMCParticle*) fMCarray->At(MCIndex);
        Int_t PDGCode = TMath::Abs(particle->GetPdgCode());
        
        Bool_t IsEnhancedPi0Eta = kFALSE;
        Bool_t IsEnhancedHF = kFALSE;
        
        Double_t Eta = particle->Eta();
        if (Eta<-0.8 || Eta > 0.8)
            continue;
        
        CocktailType_t Type = FindTrackGenerator(MCIndex, fMCheader,fMCarray);
        
        if ( Type == kBackgroundEnhanced || Type == kHFEnhanced)
            continue;
        
        if (particle->IsPrimary() && particle->GetMother()<0)
        {
            if (PDGCode == 111)
                fPtMCpi0_NoMother->Fill(particle->Pt());
            else if (PDGCode == 221)
                fPtMCEta_NoMother->Fill(particle->Pt());
        }

    }
}

void AliAnalysisTaskHFEpACorrelation::TaggingEfficiencyCalculation(AliVTrack *track,Bool_t* lIsNHFe,Bool_t* lIsHFe,Bool_t* lIsOther,Bool_t* lHasMother)
{
    if (fAOD->GetRunNumber()< 200000)
        TaggingEfficiencyCalculationRun1(track,lIsNHFe,lIsHFe,lIsOther,lHasMother);
    else
        TaggingEfficiencyCalculationRun2(track,lIsNHFe,lIsHFe,lIsOther,lHasMother);
}

void AliAnalysisTaskHFEpACorrelation::TaggingEfficiencyCalculationRun1(AliVTrack *track,Bool_t *lIsNHFe, Bool_t *lIsHFe,Bool_t *lIsOther,Bool_t *lHasMother)
{
    if(fIsMC)
    {
        if(track->GetLabel() < 0)
        {
            fElectronNoLabel->Fill(track->Pt());
            if (fNonHFE->IsULS(),fNonHFE->GetNULS())
                fElectronNoLabelULS->Fill(track->Pt());
            if (fNonHFE->IsLS())
                fElectronNoLabelLS->Fill(track->Pt(),fNonHFE->GetNLS());
        }
        else
        {
            if (fNonHFE->IsULS())
                fEtaCutElectronBKWithLabelULS->Fill(track->Pt(),fNonHFE->GetNULS());
            if (fNonHFE->IsLS())
                fEtaCutElectronBKWithLabelLS->Fill(track->Pt(),fNonHFE->GetNLS());
        }
        
        if(fIsAOD)
        {
            fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
            
            Int_t ElectronPDG = TMath::Abs(fMCparticle->GetPdgCode());
            
            if (ElectronPDG == 11)
            {
                fEtaCutElectronInclusiveRecoPtEtaZvtx->Fill(track->Pt());
                
                if(fMCparticle->GetMother()>=0)
                {
                    *lHasMother = kTRUE;
                    fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                    
                    Int_t MotherPDGAfterReco = TMath::Abs(fMCparticleMother->GetPdgCode());
                    Int_t MotherPDGHeavy  = Int_t (MotherPDGAfterReco / TMath::Power(10, Int_t(TMath::Log10(MotherPDGAfterReco))));
                    
                    //NHFE
                    if( MotherPDGAfterReco==22 || MotherPDGAfterReco ==111 || MotherPDGAfterReco ==221)
                    {
                        *lIsNHFe = kTRUE;
                        fEtaCutElectronBKNoTag->Fill(track->Pt());
                        
                        if (fNonHFE->IsULS())
                        {
                            fEtaCutElectronBKULSMainSources->Fill(track->Pt(),fNonHFE->GetNULS());
                            fEtaCutElectronBKULSMainSources_NW->Fill(track->Pt());
                        }
                        if (fNonHFE->IsLS())
                        {
                            fEtaCutElectronBKLSMainSources->Fill(track->Pt(),fNonHFE->GetNLS());
                            fEtaCutElectronBKLSMainSources_NW->Fill(track->Pt());
                            
                        }
                        
                        //UseWeights for background
                        
                        if ( MotherPDGAfterReco ==111 || MotherPDGAfterReco ==221 )
                        {
                            Double_t mPt=fMCparticleMother->Pt();
                            Double_t mweight=1;
                            //________________________________________________________________
                            //correction for d3 based on data
                            mweight=CalculateWeight(MotherPDGAfterReco, mPt);
                            
                            fEtaCutElectronBKNoTag_WithMotherW->Fill(track->Pt(), 1./mweight );
                            
                            if (fNonHFE->IsULS())
                            {
                                fEtaCutElectronBKULSMainSources_WithMotherW->Fill(track->Pt(),fNonHFE->GetNULS() * 1./mweight );
                                fEtaCutElectronBKULSMainSources_WithMotherW_NW->Fill(track->Pt(), 1./mweight );
                            }
                            
                            if (fNonHFE->IsLS())
                            {
                                fEtaCutElectronBKLSMainSources_WithMotherW->Fill(track->Pt(),fNonHFE->GetNULS() * 1./mweight );
                                fEtaCutElectronBKLSMainSources_WithMotherW_NW->Fill(track->Pt(), 1./mweight);
                            }
                        }
                        else if (fMCparticleMother->GetMother() >=0)
                        {
                            //correcting also gamma from pi0 and eta
                            fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
                            Int_t GMotherPDG = TMath::Abs(fMCparticleGMother->GetPdgCode());
                            
                            if (GMotherPDG == 111 || GMotherPDG == 221)
                            {
                                Double_t GMotherPt = fMCparticleGMother->Pt();
                                Double_t Gweight = 1;
                                
                                Gweight = CalculateWeight(GMotherPDG, GMotherPt);
                                
                                fEtaCutElectronBKNoTag_WithMotherW->Fill(track->Pt(), 1./Gweight );
                                
                                if (fNonHFE->IsULS())
                                {
                                    fEtaCutElectronBKULSMainSources_WithMotherW->Fill(track->Pt(),fNonHFE->GetNULS() * 1./Gweight );
                                    fEtaCutElectronBKULSMainSources_WithMotherW_NW->Fill(track->Pt(), 1./Gweight );
                                }
                                
                                if (fNonHFE->IsLS())
                                {
                                    fEtaCutElectronBKLSMainSources_WithMotherW->Fill(track->Pt(),fNonHFE->GetNULS() * 1./Gweight );
                                    fEtaCutElectronBKLSMainSources_WithMotherW_NW->Fill(track->Pt(), 1./Gweight);
                                }
                                
                                
                                
                                
                            }
                        }
                        
                        
                    }
                    else if(MotherPDGHeavy<4) //NHFE
                    {
                        fEtaCutElectronRecoOtherMC->Fill(track->Pt());
                        *lIsOther = kTRUE;
                        
                        if (fNonHFE->IsULS())
                            fEtaCutElectronBKULSOtherSources->Fill(track->Pt(),fNonHFE->GetNULS());
                        if (fNonHFE->IsLS())
                            fEtaCutElectronBKLSOtherSources->Fill(track->Pt(),fNonHFE->GetNLS());
                        
                        
                    }
                    else
                    {
                        fEtaCutElectronRecoHFEMC->Fill(track->Pt());
                        *lIsHFe = kTRUE;
                        if (fNonHFE->IsULS())
                            fEtaCutElectronHFEULS->Fill(track->Pt(),fNonHFE->GetNULS());
                        if (fNonHFE->IsLS())
                            fEtaCutElectronHFELS->Fill(track->Pt(),fNonHFE->GetNLS());
                        
                    }
                    
                }
            }
            else
            {
                fMissIDElectronsReco->Fill(track->Pt());
                if (fNonHFE->IsULS())
                    fEtaCutElectronMissIDULS->Fill(track->Pt(),fNonHFE->GetNULS());
                if (fNonHFE->IsLS())
                    fEtaCutElectronMissIDLS->Fill(track->Pt(),fNonHFE->GetNLS());
                
            }
            
        }
    }
}

void AliAnalysisTaskHFEpACorrelation::TaggingEfficiencyCalculationRun2(AliVTrack *track, Bool_t *lIsNHFe,Bool_t *lIsHFe,Bool_t *lIsOther, Bool_t *lHasMother)
{
    //Run 2 tagging efficiency calculation: use only electrons from the enhanced sample to evaluate the efficiencies
    //Require: pi0/eta <- e, pi0 and eta with No mother
    //Use the weight from the first generated pi0/eta
    //Flags (CT) not implemented now
    
    Int_t LabelMC = TMath::Abs(track->GetLabel());
    CocktailType_t Type = FindTrackGenerator(LabelMC, fMCheader, fMCarray);
    
    //Calculate using only the NON-enhanced sample (HF and eta/Pi0)
    
    if ( Type != kBackgroundEnhanced &&  Type != kHFEnhanced)
    {
        
        AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMCarray->At(LabelMC);
        if (MCParticle->GetMother()>=0)
        {
            AliAODMCParticle* MCMother = (AliAODMCParticle*) fMCarray->At(MCParticle->GetMother());
            Int_t MotherPDGAfterReco = TMath::Abs(MCMother->GetPdgCode());
            
            if( MotherPDGAfterReco==22 || MotherPDGAfterReco ==111 || MotherPDGAfterReco ==221)
            {
                //No weight
                
                fElectronBKGNoEnhTotalNumber->Fill(track->Pt());
                
                if (fNonHFE->IsULS())
                    fElectronBKGNoEnhULS->Fill(track->Pt(),fNonHFE->GetNULS());
                
                if (fNonHFE->IsLS())
                    fElectronBKGNoEnhLS->Fill(track->Pt(),fNonHFE->GetNLS());
                
                //Add Weight for Non-Enhanced sample
                
                if (MotherPDGAfterReco ==111 || MotherPDGAfterReco ==221)
                {
                    //Pi0/Eta should be primary and have no mother (same condition used to calculate the weight)

                    if (!MCMother->IsPrimary() || MCMother->GetMother()>0)
                        return;
                    
                    Double_t Weight = CalculateWeightRun2ToData(MotherPDGAfterReco,MCMother->Pt());
                    fElectronBKGNoEnhTotalNumber_WithW->Fill(track->Pt(), Weight);
                    
                    if (fNonHFE->IsULS())
                        fElectronBKGNoEnhULS_WithW->Fill(track->Pt(),fNonHFE->GetNULS()*Weight);
                    
                    if (fNonHFE->IsLS())
                        fElectronBKGNoEnhLS_WithW->Fill(track->Pt(),fNonHFE->GetNLS()*Weight);
                    
                }
                else if (MotherPDGAfterReco==22)
                {
                    if (MCMother->GetMother()>=0)
                    {
                        AliAODMCParticle* MCGMother = (AliAODMCParticle*) fMCarray->At(MCMother->GetMother());
                        Int_t GMotherPDGAfterReco = TMath::Abs(MCGMother->GetPdgCode());
                        
                        //Pi0/Eta should be primary and have no mother (same condition used to calculate the weight)
                        if (!MCGMother->IsPrimary() || MCGMother->GetMother()>0)
                            return;
                        
                        if (GMotherPDGAfterReco ==111 || GMotherPDGAfterReco ==221)
                        {
                            Double_t Weight = CalculateWeightRun2ToData(GMotherPDGAfterReco,MCGMother->Pt());
                            
                            fElectronBKGNoEnhTotalNumber_WithW->Fill(track->Pt(), Weight);
                            
                            if (fNonHFE->IsULS())
                                fElectronBKGNoEnhULS_WithW->Fill(track->Pt(),fNonHFE->GetNULS()*Weight);
                            
                            if (fNonHFE->IsLS())
                                fElectronBKGNoEnhLS_WithW->Fill(track->Pt(),fNonHFE->GetNLS()*Weight);
                            
                        }

                        
                    }
                    
                }

            }
            
            

        }
    }
    
    //Calculate tagging efficiency using all the produced Pi0 and eta weighted to data (that were obtained using mT scaling)
    
    AliAODMCParticle* MCParticleWtoData = (AliAODMCParticle*) fMCarray->At(LabelMC);
    
    if (MCParticleWtoData->GetMother()>=0)
    {
        AliAODMCParticle* MCMotheWtoData = (AliAODMCParticle*) fMCarray->At(MCParticleWtoData->GetMother());
        Int_t MotherPDGAfterReco = TMath::Abs(MCMotheWtoData->GetPdgCode());
        
        if( MotherPDGAfterReco==22 || MotherPDGAfterReco ==111 || MotherPDGAfterReco ==221)
        {
            //Eta-> No decay from Pi0 or gamma
            if (MotherPDGAfterReco ==221)
                FillHistBkgWtoData(MCMotheWtoData, track);
            
            //pi0
            if (MotherPDGAfterReco ==111)
            {
                if (MCMotheWtoData->GetMother() >=0 )
                {
                    AliAODMCParticle* MCGMotheWtoData = (AliAODMCParticle*) fMCarray->At(MCMotheWtoData->GetMother());
                    //Check if the mother is a Eta. If not, use the normal one
                    if (TMath::Abs(MCGMotheWtoData->GetMother()) == 221)
                        FillHistBkgWtoData(MCGMotheWtoData, track);
                    else
                        FillHistBkgWtoData(MCMotheWtoData, track);
                    
                }
                else
                     FillHistBkgWtoData(MCMotheWtoData, track);
            }
            
            //Photon
            if (MotherPDGAfterReco==22)
            {
                //Check if photon comes from pi0 or eta
                if (MCMotheWtoData->GetMother() >=0 )
                {
                    AliAODMCParticle* MCGMotheWtoData = (AliAODMCParticle*) fMCarray->At(MCMotheWtoData->GetMother());
                    
                    if (TMath::Abs(MCGMotheWtoData->GetMother()) == 111)
                    {
                        //Check if GGmother is an eta
                        if (MCGMotheWtoData->GetMother()>=0)
                        {
                            AliAODMCParticle* MCGGMotheWtoData = (AliAODMCParticle*) fMCarray->At(MCGMotheWtoData->GetMother());
                            if (TMath::Abs(MCGGMotheWtoData->GetMother()) == 221)
                                FillHistBkgWtoData(MCGGMotheWtoData, track);
                            else
                                FillHistBkgWtoData(MCGMotheWtoData, track);
                            
                        }
                        else
                             FillHistBkgWtoData(MCGMotheWtoData, track);
                        
                    }
                    else if (TMath::Abs(MCGMotheWtoData->GetMother()) == 221)
                    {
                        FillHistBkgWtoData(MCGMotheWtoData, track);
                    }
                    
                }
                
            }
            
            
        }
        
    }
    
    //Calculate HFe efficiency in any case of Cocktail type
    
    AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMCarray->At(LabelMC);
    
    if (MCParticle->GetMother()>=0)
    {
        AliAODMCParticle* MCMother = (AliAODMCParticle*) fMCarray->At(MCParticle->GetMother());
        Int_t MotherPDGAfterReco = TMath::Abs(MCMother->GetPdgCode());
        Int_t MotherPDGHeavy  = Int_t (MotherPDGAfterReco / TMath::Power(10, Int_t(TMath::Log10(MotherPDGAfterReco))));
        
        //HFe
        if (MotherPDGHeavy >3)
            fEtaCutElectronRecoHFEMC->Fill(track->Pt());

    }
    
}


void AliAnalysisTaskHFEpACorrelation::FillHistBkgWtoData(AliAODMCParticle* ParticleToUseForW, AliVTrack *track)
{
    Double_t Weight = CalculateWeightRun2ToData(TMath::Abs(ParticleToUseForW->GetPdgCode()), track->Pt());
    
    fElectronBKGWToDataTotal->Fill(track->Pt(),Weight);
    
    if (fNonHFE->IsULS())
        fElectronBKGWToDataULS->Fill(track->Pt(),fNonHFE->GetNULS()*Weight);
    if (fNonHFE->IsLS())
        fElectronBKGWToDataLS->Fill(track->Pt(),fNonHFE->GetNULS()*Weight);
}


Double_t AliAnalysisTaskHFEpACorrelation::CalculateWeightRun2ToData(Int_t pdg_particle, Double_t pT)
{
    
    if (TMath::Abs(pdg_particle) == 111)
    {
        if (!fBkgPi0WeightToData)
            return 1.0;
        
        Int_t bin = fBkgPi0WeightToData->FindBin(pT);
        return 1./fBkgPi0WeightToData->GetBinContent(bin);
    }
    else if (TMath::Abs(pdg_particle) == 221)
    {
        if (!fBkgEtaWeightToData)
            return 1.0;
        
        Int_t bin = fBkgEtaWeightToData->FindBin(pT);
        return 1./fBkgEtaWeightToData->GetBinContent(bin);
    }
    
    return 1.0;
    
}

//Methods to get the generator type and check if the particle is enhanced (adapted from AliVertexingHFUtils)

TString AliAnalysisTaskHFEpACorrelation::GetGenerator(Int_t label, AliAODMCHeader* header){
    Int_t nsumpart=0;
    TList *lh=header->GetCocktailHeaders();
    Int_t nh=lh->GetEntries();
    
    for(Int_t i=0;i<nh;i++)
    {
        AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(i);
        TString genname=gh->GetName();
        Int_t npart=gh->NProduced();
        if(label>=nsumpart && label<(nsumpart+npart))
            return genname;
        nsumpart+=npart;
    }
    TString empty="";
    return empty;
}

void AliAnalysisTaskHFEpACorrelation::GetTrackPrimaryGenerator(Int_t lab, AliAODMCHeader *header,TClonesArray *arrayMC,TString &nameGen){
    
    lab = TMath::Abs(lab);
    
    nameGen=GetGenerator(lab,header);
    
    
    while(nameGen.IsWhitespace()){
        AliAODMCParticle *mcpart= (AliAODMCParticle*)arrayMC->At(lab);
        if(!mcpart){
            break;
        }
        Int_t mother = mcpart->GetMother();
        if(mother<0){
            break;
        }
        lab=mother;
        nameGen=GetGenerator(mother,header);
    }
    
    return;
}


AliAnalysisTaskHFEpACorrelation::CocktailType_t AliAnalysisTaskHFEpACorrelation::FindTrackGenerator(Int_t label,AliAODMCHeader *header,TClonesArray *arrayMC){
    TString nameGen;
    GetTrackPrimaryGenerator(label,header,arrayMC,nameGen);
    
    if (nameGen.IsWhitespace())
        return kNoCoktail;
    if (nameGen.Contains("Hijing"))
        return kHijing;
    if (nameGen.Contains("pi0") || nameGen.Contains("eta"))
        return kBackgroundEnhanced;
    if (nameGen.Contains("bele") || nameGen.Contains("cele"))
        return kHFEnhanced;
    return kUndefined;
}






