
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
//		version: April  26, 2016.							          //
//                                                                    //
//	    Authors 							                          //
//		Elienos Pereira de Oliveira Filho (epereira@cern.ch)	      //
//		Cristiane Jahnke		(cristiane.jahnke@cern.ch)            //
//                                                                    //
//      Updates                                                       //
//      Henrique Zanoli (h.zanoli@cern.ch)                            //
//      Alexis Mas (aleximas@if.usp.br)                               //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TTree.h"
#include "TNtuple.h"
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
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliEMCALTrack.h"
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
#include "THnSparse.h"
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
#include "TRefArray.h"
#include "TVector.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliEventPoolManager.h"
#include "TObjArray.h"
//include to use reader as Lucile does
#include "AliCaloTrackAODReader.h"
#include "AliCaloTrackReader.h"
#include "AliEMCALRecoUtils.h" //to remove exotics
#include "AliAODHeader.h"
#include "AliEMCALGeometry.h"

#include "AliAnalysisUtils.h"



// --- ANALYSIS system ---
#include "AliCalorimeterUtils.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAODPWG4Particle.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliMixedEvent.h"
#include "AliAODCaloCluster.h"
#include "AliOADBContainer.h"
#include "AliAnalysisManager.h"

// --- Detector ---
#include "AliEMCALGeometry.h"
#include "AliPHOSGeoUtils.h"

//______________________________________________________________________

//______________________________________________________________________
ClassImp(AliAnalysisTaskHFEpACorrelation)

//______________________________________________________________________
AliAnalysisTaskHFEpACorrelation::AliAnalysisTaskHFEpACorrelation(const char *name)
: AliAnalysisTaskSE(name)
,fCorrelationFlag(0)
,fIspp(kFALSE)
,fIsMC(0)
,fUseEMCal(kFALSE)
,fUseTrigger(kFALSE)
,fUseTender(kFALSE)
,fUseShowerShapeCut(kFALSE)
,fCalibrateTPC(0)
,fCalibrateTPC_mean(0)
,fCalibrateTPC_sigma(1)
,fFillBackground(kFALSE)
,fEoverPnsigma(kFALSE)
,fAssocWithSPD(kFALSE)
,fEMCEG1(kFALSE)
,fEMCEG2(kFALSE)
,fIsHFE1(kFALSE)
,fIsHFE2(kFALSE)
,fIsNonHFE(kFALSE)
,fIsFromD(kFALSE)
,fIsFromB(kFALSE)
,fIsFromPi0(kFALSE)
,fIsFromEta(kFALSE)
,fIsFromGamma(kFALSE)
,fESD(0)
,fAOD(0)
,fVevent(0)
,fPartnerCuts(new AliESDtrackCuts())
,fOutputList(0)
,fPidResponse(0)
,fNonHFE(new AliSelectNonHFE())
,fIsAOD(kFALSE)
,fCentrality(0)
,fCentralityMin(0)
,fCentralityMax(100)
,fHasCentralitySelection(kFALSE)
,fUseDCACutforHadrons(kFALSE)
,fUseAlternativeBinnig(kFALSE)
,fCentralityHist(0)
,fCentralityHistPass(0)
,fZvtx(0)
,fEstimator(0)
,fClus(0)
//,fClusESD(0)
,fNevent(0)
,fNevent2(0)
,fPtElec_Inc(0)
,fPtPrim(0)
,fPtSec(0)
,fPtPrim2(0)
,fPtSec2(0)
,fCharge_n(0)
,fCharge_p(0)
,fTime(0)
,fTime2(0)
,ftimingEle(0)
,ftimingEle2(0)
,fPtElec_ULS(0)
,fPtElec_LS(0)

,fPtElec_ULS_NoPid(0)
,fPtElec_LS_NoPid(0)
,fPtElec_ULS_MC(0)
,fPtElec_ULS_MC_weight(0)
,fPtElec_ULS2(0)
,fPtElec_LS2(0)

,fPtElec_ULS_mc_closure(0)
,fPtElec_LS_mc_closure(0)
,fPtElec_ULS2_mc_closure(0)
,fPtElec_LS2_mc_closure(0)

,fPtElec_ULS_weight(0)
,fPtElec_LS_weight(0)
,fPtElec_ULS2_weight(0)
,fPtElec_LS2_weight(0)
,fTOF01(0)
,fTOF02(0)
,fTOF03(0)
,fpid(0)
,fEoverP_pt_true_electrons(0)
,fEoverP_pt_true_electrons_weight(0)
,fEoverP_pt_true_HFE(0)
,fEoverP_pt_not_HFE(0)

,fEoverP_ntracks_matched(0)

,fEoverP_ncells(0)

,fEmc_Ereco_gamma0(0)
,fEmc_Ereco_gamma_ratio0(0)
,fEmc_Ereco_ele0(0)
,fEmc_Ereco_ele_ratio0(0)

,fEmc_Ereco_gamma_all(0)
,fEmc_Ereco_gamma_ratio_all(0)
,fEmc_Ereco_ele_all(0)
,fEmc_Ereco_ele_ratio_all(0)

,fEmc_Ereco_gamma(0)
,fEmc_Ereco_gamma_ratio(0)
,fEmc(0)
,fEmc_Ereco_ele(0)
,fEmc_Ereco_ele_ratio(0)

,fEmc_Ereco_ratio_large_EoverP0(0)
,fEmc_Ereco_ratio_small_EoverP0(0)

,fEmc_Ereco_ratio_large_EoverP(0)
,fEmc_Ereco_ratio_small_EoverP(0)

,fEoverP_pt_true_hadrons(0)
,fEoverP_pt_true_electrons0(0)
,fEoverP_pt_true_hadrons0(0)
,fEoverP_pt(0)
,fEoverP_pt_highE0(0)
,fEoverP_pt_highE1(0)
,fEoverP_tpc(0)
,fEoverP_tpc_p_trigger(0)
,fEoverP_tpc_pt_trigger(0)
,fTPC_pt(0)
,fTPC_p(0)

,fTPC_momentum(0)
,fTPC_eta(0)
,fTPC_momentum1(0)
,fTPC_eta1(0)

,fTPCnsigma_pt(0)
,fTPCnsigma_p(0)

,fTPCnsigma_p_eta1(0)
,fTPCnsigma_p_eta2(0)
,fTPCnsigma_p_eta3(0)
,fTPCnsigma_p_eta4(0)

,fTPCnsigma_p_TPC(0)
,fTPCnsigma_p_TPC_on_EMCal_acc(0)
,fTPCnsigma_p_TPC_EoverP_cut(0)
,fTPCnsigma_pt_2D(0)
,fTPCnsigma_pt_2D0(0)
,fTPCnsigma_pt_2D1(0)
,fTPCnsigma_pt_2D2(0)
,fTPCnsigma_pt_2D3(0)
,fTPCnsigma_pt_2D4(0)
,fTPCnsigma_pt_2D5(0)
,fShowerShapeCut(0)
,fShowerShapeM02_EoverP(0)
,fShowerShapeM20_EoverP(0)

,fShowerShapeM02_EoverP_Ecut12_MC(0)
,fShowerShapeM20_EoverP_Ecut12_MC(0)
,fShowerShapeM02_EoverP_Ecut8_MC(0)
,fShowerShapeM20_EoverP_Ecut8_MC(0)

,fShowerShape_ha(0)
,fShowerShape_ele(0)
,fTPCnsigma_eta(0)
,fTPCnsigma_phi(0)
,fECluster(0)
,fECluster_pure(0)
,fECluster_highpT0(0)
,fECluster_highpT1(0)
,fECluster_not_exotic(0)
,fECluster_not_exotic1(0)
,fECluster_not_exotic2(0)
,fECluster_exotic(0)

,fNcells_Energy_highE0(0)
,fEtaPhi_highE0(0)
,fEta_highE0(0)
,fPhi_highE0(0)
,fR_highE0(0)
,fNCluster_highE0(0)
,fNcells_Energy_highE1(0)
,fEtaPhi_highE1(0)
,fEta_highE1(0)
,fPhi_highE1(0)
,fR_highE1(0)
,fNCluster_highE1(0)

,fNCluster_pure(0)
,fNCluster_pure_aod(0)
,fNCluster_ECluster(0)
,fNcells_energy(0)
,fNcells_energy_elec_selected(0)
,fNcells_energy_not_exotic(0)

,fEtaPhi(0)
,fEtaPhi_large_EoverP(0)
,fEtaPhi_small_EoverP(0)
,fEtaPhi_num(0)
,fEtaPhi_den(0)
,fEtaPhi_data(0)
,fpt_reco_pt_MC_num(0)
,fpt_reco_pt_MC_den(0)
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
,fNClusters(0)
,fTPCNcls_EoverP(0)
,fTPCNcls_pid(0)
,fEta(0)
,fPhi(0)
,fR(0)
,fR_EoverP(0)
,fNcells(0)
,fNcells_EoverP(0)
,fNcells_electrons(0)
,fNcells_hadrons(0)
,fECluster_ptbins(0)
,fEoverP_ptbins(0)
,fEoverP_wSSCut(0)
,fM02_EoverP(0)
,fM20_EoverP(0)
,fTPCnsigma_eta_electrons(0)
,fTPCnsigma_eta_hadrons(0)
,fEoverP_pt_pions(0)
,ftpc_p_EoverPcut(0)
,fnsigma_p_EoverPcut(0)

,fEoverP_pt_pions2(0)
,fEoverP_pt_pions3(0)
,fEoverP_pt_pions2_highE0(0)
,fEoverP_pt_pions2_highE1(0)

,fNcells_pt(0)
,fEoverP_pt_hadrons(0)
,fCEtaPhi_Inc(0)
,fCEtaPhi_ULS(0)
,fCEtaPhi_LS(0)
,fCEtaPhi_ULS_NoP(0)
,fCEtaPhi_LS_NoP(0)
,fCEtaPhi_ULS_Weight(0)
,fCEtaPhi_LS_Weight(0)
,fCEtaPhi_ULS_NoP_Weight(0)
,fCEtaPhi_LS_NoP_Weight(0)
,fCetaPhi_MC_with_partner_greater(0) //new 15-09
,fCetaPhi_MC_with_partner_below(0) // new 15-09
,fCetaPhi_MC_NHFE_1partner_reco(0) // new 17-09
,fpT_MC_with_partner_greater(0)
,fpT_MC_with_partner_below(0)
,fMCEffPID_beforePID(0)
,fMCEffPID_afterPID(0)
,fpTMCGeneratedBKG(0)
,fInvMass(0)
,fInvMassBack(0)
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
,fEtaCutMin(-0.9)
,fEtaCutMax(0.9)
,fTPCcal_CutMin(-1)
,fTPCcal_CutMax(3)
,fdPhiCut(0.05)
,fdEtaCut(0.05)
,fEoverPCutMin(0.8)
,fEoverPCutMax(1.2)
,fM20CutMin(0.0)
,fM20CutMax(10)
,fM02CutMin(0.0)
,fM02CutMax(10)
,fAngleCut(999)
,fChi2Cut(3.5)
,fDCAcut(999)//dca between two tracks
,fDCAcutr(999)//dca to vertex
,fDCAcutz(999)//dca to vertex
,fMassCutFlag(kTRUE)
,fAngleCutFlag(kFALSE)
,fChi2CutFlag(kFALSE)
,fDCAcutFlag(kFALSE)
,fAssHadronPtMin(0.5)
,fAssHadronPtMax(2.0)
,fPtBackgroundBeforeReco(0)
,fPtBackgroundBeforeReco2(0)
,fPtBackgroundBeforeReco_weight(0)
,fPtBackgroundBeforeReco2_weight(0)
,fpT_m_electron(0)
,fpT_gm_electron(0)
,fPtBackgroundAfterReco(0)
,fPtMinAsso(0.3)
,fTpcNclsAsso(80)
,fPtMCparticleAll(0)
,fPtMCparticleAll_nonPrimary(0)
,fPtMCparticleAlle_nonPrimary(0)
,fPtMCparticleAlle_Primary(0)
,fPtMCparticleReco(0)
,fPtMCparticleReco_nonPrimary(0)
,fPtMCparticleAllHfe1(0)
,fPtMCparticleRecoHfe1(0)
,fPtMCparticleAllHfe2(0)
,fPtMCparticleRecoHfe2(0)
,fPtMCelectronAfterAll(0)
,fPtMCelectronAfterAll_unfolding(0)
,fPtMCelectronAfterAll_nonPrimary(0)
,fPtMCelectronAfterAll_Primary(0)
,fPtMCpi0(0)
,fPtMCeta(0)
,fPtMCpi02(0)
,fPtMCeta2(0)
,fPtMCpi03(0)
,fPtMCeta3(0)
,fPtMC_EMCal_All(0)
,fPtMC_EMCal_Selected(0)
,fPtMC_TPC_All(0)
,fPtMC_TPC_Selected(0)
,fPt_track_match_den(0)
,fPt_track_match_num(0)
,fPtMCWithLabel(0)
,fPtMCWithoutLabel(0)
,fPtIsPhysicaPrimary(0)
,fCuts(0)
//,reader(0)
,fCFM(0)
,fPID(new AliHFEpid("hfePid"))
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
,fCEtaPhi_ULS_EM(0)
,fCEtaPhi_LS_EM(0)
,fCEtaPhi_ULS_Weight_EM(0)
,fCEtaPhi_LS_Weight_EM(0)
,fPoolNevents(0)
,fEventMixingFlag(0)
,fCEtaPhi_Inc_DiHadron(0)
,fPtTrigger_Inc(0)
//,fEMCALRecoUtils(new AliEMCALRecoUtils)
//,fEMCALGeo(0x0)
//,fCaloUtils(0x0)

,fBitEGA(0)
,fAnalysisUtils(0)
,fCEtaPhi_Inc_PH(0)
,fCEtaPhi_ULS_Weight_PH(0)
,fCEtaPhi_LS_Weight_PH(0)
,fCEtaPhi_ULS_NoP_Weight_PH(0)
,fCEtaPhi_LS_NoP_Weight_PH(0)
,fpT_MC_HFE_RECO_MC(0)
,fCetaPhi_MC_HFE_RECO_MC_PhyPrimH(0)
,fCetaPhi_Data_Data_RECO_MC_PhyPrimH(0)
,fpT_Data_HFE_RECO_Data(0)
,fCetaPhi_MC_HFE_pTofReco(0)
,fCetaPhi_Reco_with_MChadrons(0)
,fpTReco_vs_MC(0)
,fpTShifHadronsMC(0)
,fpTShiftHadronsReco(0)
,fpTEtaEffHadronsReco(0)
,fpTPhiEffHadronsReco(0)
,fpTEtaEffHadronsMC(0)
,fpTPhiEffHadronsMC(0)
,fCEtaPhi_DataHFE_with_onlyPhysPriHadron(0)
,fDCAcutrHadron(999)
,fDCAcutzHadron(999)
,fpTPhiHadron(0)
,fpTEtaHadron(0)
,fpTPhiSecPartNoLabel(0)
,fpTEtaSecPartNoLabel(0)
,fpTPhiSecPartPhysPri(0)
,fpTEtaSecPartPhysPri(0)
,fpTPhiSecPartNonPhysPri(0)
,fpTEtaSecPartNonPhysPri(0)
//,fEMCALRecoUtils(0)//exotic
{
    //Named constructor
    // Define input and output slots here
    // Input slot #0 works with a TChain
    
    //exotic
    //fEMCALRecoUtils  = new AliEMCALRecoUtils();
    
    fAnalysisUtils = new AliAnalysisUtils;
    
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
,fIspp(kFALSE)
,fIsMC(0)
,fUseEMCal(kFALSE)
,fUseTrigger(kFALSE)
,fUseTender(kFALSE)
,fUseShowerShapeCut(kFALSE)
,fCalibrateTPC(0)
,fCalibrateTPC_mean(0)
,fCalibrateTPC_sigma(1)
,fFillBackground(kFALSE)
,fEoverPnsigma(kFALSE)
,fAssocWithSPD(kFALSE)
,fEMCEG1(kFALSE)
,fEMCEG2(kFALSE)
,fIsHFE1(kFALSE)
,fIsHFE2(kFALSE)
,fIsNonHFE(kFALSE)
,fIsFromD(kFALSE)
,fIsFromB(kFALSE)
,fIsFromPi0(kFALSE)
,fIsFromEta(kFALSE)
,fIsFromGamma(kFALSE)
,fESD(0)
,fAOD(0)
,fVevent(0)
,fPartnerCuts(new AliESDtrackCuts())
,fOutputList(0)
,fPidResponse(0)
,fNonHFE(new AliSelectNonHFE())
,fIsAOD(kFALSE)
,fCentrality(0)
,fCentralityMin(0)
,fCentralityMax(100)
,fHasCentralitySelection(kFALSE)
,fUseDCACutforHadrons(kFALSE)
,fUseAlternativeBinnig(kFALSE)
,fCentralityHist(0)
,fCentralityHistPass(0)
,fZvtx(0)
,fEstimator(0)
,fClus(0)
//,fClusESD(0)
,fNevent(0)
,fNevent2(0)
,fPtElec_Inc(0)
,fPtPrim(0)
,fPtSec(0)
,fPtPrim2(0)
,fPtSec2(0)
,fCharge_n(0)
,fCharge_p(0)
,fTime(0)
,fTime2(0)
,ftimingEle(0)
,ftimingEle2(0)
,fPtElec_ULS(0)
,fPtElec_LS(0)
,fPtElec_ULS_NoPid(0)
,fPtElec_LS_NoPid(0)
,fPtElec_ULS_MC(0)
,fPtElec_ULS_MC_weight(0)
,fPtElec_ULS2(0)
,fPtElec_LS2(0)

,fPtElec_ULS_mc_closure(0)
,fPtElec_LS_mc_closure(0)
,fPtElec_ULS2_mc_closure(0)
,fPtElec_LS2_mc_closure(0)

,fPtElec_ULS_weight(0)
,fPtElec_LS_weight(0)
,fPtElec_ULS2_weight(0)
,fPtElec_LS2_weight(0)
,fTOF01(0)
,fTOF02(0)
,fTOF03(0)
,fpid(0)
,fEoverP_pt_true_electrons(0)
,fEoverP_pt_true_electrons_weight(0)
,fEoverP_pt_true_HFE(0)
,fEoverP_pt_not_HFE(0)

,fEoverP_ntracks_matched(0)
,fEoverP_ncells(0)

,fEmc_Ereco_gamma0(0)
,fEmc_Ereco_gamma_ratio0(0)

,fEmc_Ereco_ele0(0)
,fEmc_Ereco_ele_ratio0(0)

,fEmc_Ereco_gamma_all(0)
,fEmc_Ereco_gamma_ratio_all(0)
,fEmc_Ereco_ele_all(0)
,fEmc_Ereco_ele_ratio_all(0)


,fEmc_Ereco_gamma(0)
,fEmc_Ereco_gamma_ratio(0)
,fEmc(0)
,fEmc_Ereco_ele(0)
,fEmc_Ereco_ele_ratio(0)

,fEmc_Ereco_ratio_large_EoverP0(0)
,fEmc_Ereco_ratio_small_EoverP0(0)

,fEmc_Ereco_ratio_large_EoverP(0)
,fEmc_Ereco_ratio_small_EoverP(0)

,fEoverP_pt_true_hadrons(0)
,fEoverP_pt_true_electrons0(0)

,fEoverP_pt_true_hadrons0(0)
,fEoverP_pt(0)
,fEoverP_pt_highE0(0)
,fEoverP_pt_highE1(0)
,fEoverP_tpc(0)
,fEoverP_tpc_p_trigger(0)
,fEoverP_tpc_pt_trigger(0)
,fTPC_pt(0)
,fTPC_p(0)
,fTPC_momentum(0)
,fTPC_eta(0)
,fTPC_momentum1(0)
,fTPC_eta1(0)
,fTPCnsigma_pt(0)

,fTPCnsigma_p(0)
,fTPCnsigma_p_eta1(0)
,fTPCnsigma_p_eta2(0)
,fTPCnsigma_p_eta3(0)
,fTPCnsigma_p_eta4(0)

,fTPCnsigma_p_TPC(0)
,fTPCnsigma_p_TPC_on_EMCal_acc(0)
,fTPCnsigma_p_TPC_EoverP_cut(0)
,fTPCnsigma_pt_2D(0)
,fTPCnsigma_pt_2D0(0)
,fTPCnsigma_pt_2D1(0)
,fTPCnsigma_pt_2D2(0)
,fTPCnsigma_pt_2D3(0)
,fTPCnsigma_pt_2D4(0)
,fTPCnsigma_pt_2D5(0)
,fShowerShapeCut(0)
,fShowerShapeM02_EoverP(0)
,fShowerShapeM20_EoverP(0)

,fShowerShapeM02_EoverP_Ecut12_MC(0)
,fShowerShapeM20_EoverP_Ecut12_MC(0)
,fShowerShapeM02_EoverP_Ecut8_MC(0)
,fShowerShapeM20_EoverP_Ecut8_MC(0)

,fShowerShape_ha(0)
,fShowerShape_ele(0)
,fTPCnsigma_eta(0)
,fTPCnsigma_phi(0)
,fECluster(0)
,fECluster_pure(0)
,fECluster_highpT0(0)
,fECluster_highpT1(0)
,fECluster_not_exotic(0)
,fECluster_not_exotic1(0)
,fECluster_not_exotic2(0)
,fECluster_exotic(0)

,fNcells_Energy_highE0(0)
,fEtaPhi_highE0(0)
,fEta_highE0(0)
,fPhi_highE0(0)
,fR_highE0(0)
,fNCluster_highE0(0)
,fNcells_Energy_highE1(0)
,fEtaPhi_highE1(0)
,fEta_highE1(0)
,fPhi_highE1(0)
,fR_highE1(0)
,fNCluster_highE1(0)

,fNCluster_pure(0)
,fNCluster_pure_aod(0)
,fNCluster_ECluster(0)
,fNcells_energy(0)
,fNcells_energy_elec_selected(0)
,fNcells_energy_not_exotic(0)
,fEtaPhi(0)
,fEtaPhi_large_EoverP(0)
,fEtaPhi_small_EoverP(0)
,fEtaPhi_num(0)
,fEtaPhi_den(0)
,fEtaPhi_data(0)
,fpt_reco_pt_MC_num(0)
,fpt_reco_pt_MC_den(0)
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
,fNClusters(0)
,fTPCNcls_EoverP(0)
,fTPCNcls_pid(0)
,fEta(0)
,fPhi(0)
,fR(0)
,fR_EoverP(0)
,fNcells(0)
,fNcells_EoverP(0)
,fNcells_electrons(0)
,fNcells_hadrons(0)
,fECluster_ptbins(0)
,fEoverP_ptbins(0)
,fEoverP_wSSCut(0)
,fM02_EoverP(0)
,fM20_EoverP(0)
,fTPCnsigma_eta_electrons(0)
,fTPCnsigma_eta_hadrons(0)
,fEoverP_pt_pions(0)
,ftpc_p_EoverPcut(0)
,fnsigma_p_EoverPcut(0)
,fEoverP_pt_pions2(0)
,fEoverP_pt_pions3(0)
,fEoverP_pt_pions2_highE0(0)
,fEoverP_pt_pions2_highE1(0)
,fNcells_pt(0)
,fEoverP_pt_hadrons(0)
,fCEtaPhi_Inc(0)
,fCEtaPhi_ULS(0)
,fCEtaPhi_LS(0)
,fCEtaPhi_ULS_NoP(0)
,fCEtaPhi_LS_NoP(0)
,fCEtaPhi_ULS_Weight(0)
,fCEtaPhi_LS_Weight(0)
,fCEtaPhi_ULS_NoP_Weight(0)
,fCEtaPhi_LS_NoP_Weight(0)
,fCetaPhi_MC_with_partner_greater(0) //new 15-09
,fCetaPhi_MC_with_partner_below(0) // new 15-09
,fCetaPhi_MC_NHFE_1partner_reco(0) // new 17-09
,fpT_MC_with_partner_greater(0)
,fpT_MC_with_partner_below(0)
,fMCEffPID_beforePID(0)
,fMCEffPID_afterPID(0)
,fpTMCGeneratedBKG(0)
,fInvMass(0)
,fInvMassBack(0)
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
,fEtaCutMin(-0.9)
,fEtaCutMax(0.9)
,fTPCcal_CutMin(-1)
,fTPCcal_CutMax(3)
,fdPhiCut(0.05)
,fdEtaCut(0.05)
,fEoverPCutMin(0.8)
,fEoverPCutMax(1.2)
,fM20CutMin(0)
,fM20CutMax(10)
,fM02CutMin(0)
,fM02CutMax(10)
,fAngleCut(999)
,fChi2Cut(3.5)
,fDCAcut(999)//dca between two tracks
,fDCAcutr(999)//dca to vertex
,fDCAcutz(999)//dca to vertex
,fMassCutFlag(kTRUE)
,fAngleCutFlag(kFALSE)
,fChi2CutFlag(kFALSE)
,fDCAcutFlag(kFALSE)
,fAssHadronPtMin(0.5)
,fAssHadronPtMax(2.0)
,fPtBackgroundBeforeReco(0)
,fPtBackgroundBeforeReco2(0)
,fPtBackgroundBeforeReco_weight(0)
,fPtBackgroundBeforeReco2_weight(0)
,fpT_m_electron(0)
,fpT_gm_electron(0)
,fPtBackgroundAfterReco(0)
,fPtMinAsso(0.3)
,fTpcNclsAsso(80)
,fPtMCparticleAll(0)
,fPtMCparticleAll_nonPrimary(0)
,fPtMCparticleAlle_nonPrimary(0)
,fPtMCparticleAlle_Primary(0)
,fPtMCparticleReco(0)
,fPtMCparticleReco_nonPrimary(0)
,fPtMCparticleAllHfe1(0)
,fPtMCparticleRecoHfe1(0)
,fPtMCparticleAllHfe2(0)
,fPtMCparticleRecoHfe2(0)
,fPtMCelectronAfterAll(0)
,fPtMCelectronAfterAll_unfolding(0)
,fPtMCelectronAfterAll_nonPrimary(0)
,fPtMCelectronAfterAll_Primary(0)
,fPtMCpi0(0)
,fPtMCeta(0)
,fPtMCpi02(0)
,fPtMCeta2(0)
,fPtMCpi03(0)
,fPtMCeta3(0)
,fPtMC_EMCal_All(0)
,fPtMC_EMCal_Selected(0)
,fPtMC_TPC_All(0)
,fPtMC_TPC_Selected(0)
,fPt_track_match_den(0)
,fPt_track_match_num(0)
,fPtMCWithLabel(0)
,fPtMCWithoutLabel(0)
,fPtIsPhysicaPrimary(0)
,fCuts(0)
//Lucile
//,reader(0)
,fCFM(0)
,fPID(new AliHFEpid("hfePid"))
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
,fCEtaPhi_ULS_EM(0)
,fCEtaPhi_LS_EM(0)
,fCEtaPhi_ULS_Weight_EM(0)
,fCEtaPhi_LS_Weight_EM(0)
,fPoolNevents(0)
,fEventMixingFlag(0)
,fCEtaPhi_Inc_DiHadron(0)
,fPtTrigger_Inc(0)
//,fEMCALRecoUtils(new AliEMCALRecoUtils)
//,fEMCALGeo(0x0)
//,fCaloUtils(0x0)
,fBitEGA(0)
,fAnalysisUtils(0)
,fCEtaPhi_Inc_PH(0)
,fCEtaPhi_ULS_Weight_PH(0)
,fCEtaPhi_LS_Weight_PH(0)
,fCEtaPhi_ULS_NoP_Weight_PH(0)
,fCEtaPhi_LS_NoP_Weight_PH(0)
,fpT_MC_HFE_RECO_MC(0)
,fCetaPhi_MC_HFE_RECO_MC_PhyPrimH(0)
,fCetaPhi_Data_Data_RECO_MC_PhyPrimH(0)
,fpT_Data_HFE_RECO_Data(0)
,fCetaPhi_MC_HFE_pTofReco(0)
,fCetaPhi_Reco_with_MChadrons(0)
,fpTReco_vs_MC(0)
,fpTShifHadronsMC(0)
,fpTShiftHadronsReco(0)
,fpTEtaEffHadronsReco(0)
,fpTPhiEffHadronsReco(0)
,fpTEtaEffHadronsMC(0)
,fpTPhiEffHadronsMC(0)
,fCEtaPhi_DataHFE_with_onlyPhysPriHadron(0)
,fDCAcutrHadron(999)
,fDCAcutzHadron(999)
,fpTPhiHadron(0)
,fpTEtaHadron(0)
,fpTPhiSecPartNoLabel(0)
,fpTEtaSecPartNoLabel(0)
,fpTPhiSecPartPhysPri(0)
,fpTEtaSecPartPhysPri(0)
,fpTPhiSecPartNonPhysPri(0)
,fpTEtaSecPartNonPhysPri(0)
//,fEMCALRecoUtils(0)//exotic
{
    // Constructor
    // Define input and output slots here
    // Input slot #0 works with a TChain
    
    //exotic
    // fEMCALRecoUtils  = new AliEMCALRecoUtils();
    
    fAnalysisUtils = new AliAnalysisUtils;
    
    
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    // DefineOutput(1, TH1I::Class());
    DefineOutput(1, TList::Class());
    //DefineOutput(3, TTree::Class());
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
    if(fAnalysisUtils) delete fAnalysisUtils;
}

//______________________________________________________________________
//Create Output Objects
//Here we can define the histograms and others output files
//Called once
void AliAnalysisTaskHFEpACorrelation::UserCreateOutputObjects()
{
    //______________________________________________________________________
    //Initialize PID
    if(!fPID->GetNumberOfPIDdetectors())
    {
        fPID->AddDetector("TPC", 0);
    }
    
    fPID->SortDetectors();
    
    fPIDqa = new AliHFEpidQAmanager();
    fPIDqa->Initialize(fPID);
    //______________________________________________________________________
    
    //______________________________________________________________________
    //Initialize correction Framework and Cuts
    fCFM = new AliCFManager;
    const Int_t kNcutSteps = AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack + AliHFEcuts::kNcutStepsDETrack;
    fCFM->SetNStepParticle(kNcutSteps);
    for(Int_t istep = 0; istep < kNcutSteps; istep++) fCFM->SetParticleCutsList(istep, NULL);
    
    if(!fCuts)
    {
        AliWarning("Cuts not available. Default cuts will be used");
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
    fOutputList->Add(fPIDqa->MakeList("PIDQA"));
    
    //Store the number of events
    //Define the histo
    fNevent = new TH1F("fNevent","Number of Events",30,0,30);
    fNevent2 = new TH1F("fNevent2","Number of Events 2",30,0,30);
    //And then, add to the output list
    fOutputList->Add(fNevent);
    fOutputList->Add(fNevent2);
    
    fpid = new TH1F("fpid","PID flag",5,0,5);
    fOutputList->Add(fpid);
    
    //pt Distribution
    fPtElec_Inc = new TH1F("fPtElec_Inc","Inclusive Electrons; p_{T} (GeV/c); Count",300,0,30);
    
    
    
    
    fPtElec_ULS = new TH1F("fPtElec_ULS","ULS; p_{T} (GeV/c); Count",300,0,30);
    fPtElec_LS = new TH1F("fPtElec_LS","LS; p_{T} (GeV/c); Count",300,0,30);
    
    
    fTOF01 = new TH2F("fTOF01","",200,-20,20,200,-20,20);
    fTOF02 = new TH2F("fTOF02","",200,-20,20,200,-20,20);
    fTOF03 = new TH2F("fTOF03","",200,-20,20,200,-20,20);
    
    if(fIsMC){
        
        //pT shift
        fpTShifHadronsMC = new TH1F("fpTShifHadronsMC","Hadrons generated pT; p_{T} (GeV/c); Count",300,0,30);
        fOutputList->Add(fpTShifHadronsMC);
        fpTShiftHadronsReco = new TH1F("fpTShiftHadronsReco","Hadrons generated pT; p_{T} (GeV/c); Count",300,0,30);
        fOutputList->Add(fpTShiftHadronsReco);
        
        //Efficiency of Hadrons
        fpTEtaEffHadronsMC = new TH2F("fpTEtaEffHadronsMC ","Hadrons generated pT; p_{T} (GeV/c);  #eta ; Count",300,0,30,200,-1,1);
        fOutputList->Add(fpTEtaEffHadronsMC);
        
        fpTPhiEffHadronsMC = new TH2F("fpTPhiEffHadronsMC","Hadrons generated pT; p_{T} (GeV/c); #varphi; Count",300,0,30,200,0,2*TMath::Pi());
        fOutputList->Add(fpTPhiEffHadronsMC);
        
        
        fpTEtaEffHadronsReco = new TH2F("fpTEtaEffHadronsReco","Hadrons generated pT; p_{T} (GeV/c); #eta; Count",300,0,30,200,-1,1);
        fOutputList->Add(fpTEtaEffHadronsReco);
        
        fpTPhiEffHadronsReco = new TH2F("fpTPhiEffHadronsReco","Hadrons generated pT; p_{T} (GeV/c); #varphi; Count",300,0,30,200,0,2*TMath::Pi());
        fOutputList->Add(fpTPhiEffHadronsReco);
        
        //Secundary particles
        fpTPhiHadron = new TH2F("fpTPhiHadron","All hadrons pT; p_{T} (GeV/c); #varphi; Count",300,0,30,200,0,2*TMath::Pi());
        fOutputList->Add(fpTPhiHadron);

        fpTEtaHadron = new TH2F("fpTEtaHadron","All hadrons pT; p_{T} (GeV/c); #eta; Count",300,0,30,200,-1,1);
        fOutputList->Add(fpTEtaHadron);

        fpTPhiSecPartNoLabel = new TH2F("fpTPhiSecPartNoLabel","No Label Hadrons pT; p_{T} (GeV/c); #varphi; Count",300,0,30,200,0,2*TMath::Pi());
        fOutputList->Add(fpTPhiSecPartNoLabel);
        //parei aqui
        fpTEtaSecPartNoLabel = new TH2F("fpTEtaSecPartNoLabel","No Label pT; p_{T} (GeV/c); #eta; Count",300,0,30,200,-1,1);
        fOutputList->Add(fpTEtaSecPartNoLabel);

        fpTPhiSecPartPhysPri = new TH2F("fpTPhiSecPartPhysPri","Phys Pri hadrons pT; p_{T} (GeV/c); #varphi; Count",300,0,30,200,0,2*TMath::Pi());
        fOutputList->Add(fpTPhiSecPartPhysPri);

        fpTEtaSecPartPhysPri = new TH2F("fpTEtaSecPartPhysPri","Phys Pri pT; p_{T} (GeV/c); #eta; Count",300,0,30,200,-1,1);
        fOutputList->Add(fpTEtaSecPartPhysPri);

        fpTPhiSecPartNonPhysPri = new TH2F("fpTPhiSecPartNonPhysPri","Non Phys Pri hadrons pT; p_{T} (GeV/c); #varphi; Count",300,0,30,200,0,2*TMath::Pi());
        fOutputList->Add(fpTPhiSecPartNonPhysPri);

        fpTEtaSecPartNonPhysPri = new TH2F("fpTEtaSecPartNonPhysPri","Non Phys Pri pT; p_{T} (GeV/c); #eta; Count",300,0,30,200,-1,1);
        fOutputList->Add(fpTEtaSecPartNonPhysPri);

        
        fCEtaPhi_DataHFE_with_onlyPhysPriHadron = new TH2F *[6];
        
        fCetaPhi_MC_HFE_pTofReco = new TH2F *[6];
        fCetaPhi_Reco_with_MChadrons = new TH2F *[6];
        
        fpTReco_vs_MC = new TH2F("fpTReco_vs_MC","pT Reco vs MC;pT Reco (Gev/c); pT MC (GeV/c)",1000,0,10,1000,0,10);
        fOutputList->Add(fpTReco_vs_MC);
        fpT_MC_HFE_RECO_MC = new TH1F("fpT_MC_HFE_RECO_MC","HFE generated pT; p_{T} (GeV/c); Count",300,0,30);
        fOutputList->Add(fpT_MC_HFE_RECO_MC);
        fpT_Data_HFE_RECO_Data =new TH1F("fpT_Data_HFE_RECO_Data","HFE reco pT; p_{T} (GeV/c); Count",300,0,30);
        fOutputList->Add(fpT_Data_HFE_RECO_Data);
        
        fCetaPhi_MC_HFE_RECO_MC_PhyPrimH = new TH2F *[6];
        fCetaPhi_Data_Data_RECO_MC_PhyPrimH = new TH2F *[6];
        
        Double_t fPtBin[7];
        if(fUseAlternativeBinnig)
        {
            fPtBin[0] = 0.5;
            fPtBin[1] = 1.0;
            fPtBin[2] = 1.5;
            fPtBin[3] = 2.0;
            fPtBin[4] = 3.0;
            fPtBin[5] = 4.0;
            fPtBin[6] = 6.0;
        }
        else
        {
            fPtBin[0] = 1.0;
            fPtBin[1] = 2.0;
            fPtBin[2] = 4.0;
            fPtBin[3] = 6.0;
            fPtBin[4] = 8.0;
            fPtBin[5] = 10.0;
            fPtBin[6] = 15.0;
        }
        
        
        for(Int_t i = 0; i < 6; i++)
        {
            fCEtaPhi_DataHFE_with_onlyPhysPriHadron[i] = new TH2F(Form("fCEtaPhi_DataHFE_with_onlyPhysPriHadron%d",i),Form("Full MC WITH PT Reco: (reco) HFE-h %1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
            fOutputList->Add(fCEtaPhi_DataHFE_with_onlyPhysPriHadron[i]);
            
            fCetaPhi_MC_HFE_pTofReco[i] = new TH2F(Form("fCetaPhi_MC_HFE_pTofReco%d",i),Form("Full MC WITH PT Reco: (reco) HFE-h %1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
            fOutputList->Add(fCetaPhi_MC_HFE_pTofReco[i]);
            
            fCetaPhi_Reco_with_MChadrons[i] = new TH2F(Form("fCetaPhi_Reco_with_MChadrons%d",i),Form("HFE(reco)-h(MC) %1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
            fOutputList->Add(fCetaPhi_Reco_with_MChadrons[i]);
            fCetaPhi_MC_HFE_RECO_MC_PhyPrimH[i] = new TH2F(Form("fCetaPhi_MC_HFE_RECO_MC_PhyPrimH%d",i),Form("Full MC: (reco) HFE-h %1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
            fCetaPhi_Data_Data_RECO_MC_PhyPrimH[i] = new TH2F(Form("CetaPhi_Data_Data_RECO_MC_PhyPrimH%d",i),Form("Full Data: (reco) HFE-h %1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
            
            fOutputList->Add(fCetaPhi_MC_HFE_RECO_MC_PhyPrimH[i]);
            fOutputList->Add(fCetaPhi_Data_Data_RECO_MC_PhyPrimH[i]);
        }
        
        fPtPrim = new TH1F("fPtPrim","Primary Electrons aod track; p_{T} (GeV/c); Count",300,0,30);
        fPtSec = new TH1F("fPtSec","Secundary Electrons aod track; p_{T} (GeV/c); Count",300,0,30);
        fPtPrim2 = new TH1F("fPtPrim2","Primary Electrons vtrack; p_{T} (GeV/c); Count",300,0,30);
        fPtSec2 = new TH1F("fPtSec2","Secundary Electrons vtrack; p_{T} (GeV/c); Count",300,0,30);
        fPtElec_ULS_NoPid = new TH1F("fPtElec_ULS_NoPid","ULS; p_{T} (GeV/c); Count",300,0,30);
        fPtElec_LS_NoPid = new TH1F("fPtElec_LS_NoPid","LS; p_{T} (GeV/c); Count",300,0,30);
        fPtElec_ULS_MC = new TH1F("fPtElec_ULS_MC","ULS; p_{T} (GeV/c); Count",300,0,30);
        fPtElec_ULS_MC_weight = new TH1F("fPtElec_ULS_MC_weight","ULS; p_{T} (GeV/c); Count",300,0,30);
        
        fPtElec_ULS_weight = new TH1F("fPtElec_ULS_weight","ULS; p_{T} (GeV/c); Count",300,0,30);
        fPtElec_LS_weight = new TH1F("fPtElec_LS_weight","LS; p_{T} (GeV/c); Count",300,0,30);
        
        fPtElec_ULS2 = new TH1F("fPtElec_ULS2","ULS; p_{T} (GeV/c); Count",300,0,30);
        fPtElec_LS2 = new TH1F("fPtElec_LS2","LS; p_{T} (GeV/c); Count",300,0,30);
        
        fPtElec_ULS2_weight = new TH1F("fPtElec_ULS2_weight","ULS; p_{T} (GeV/c); Count",300,0,30);
        fPtElec_LS2_weight = new TH1F("fPtElec_LS2_weight","LS; p_{T} (GeV/c); Count",300,0,30);
        
        //mc closure
        fPtElec_ULS_mc_closure = new TH1F("fPtElec_ULS_mc_closure","ULS; p_{T} (GeV/c); Count",300,0,30);
        fPtElec_LS_mc_closure = new TH1F("fPtElec_LS_mc_closure","LS; p_{T} (GeV/c); Count",300,0,30);
        fPtElec_ULS2_mc_closure = new TH1F("fPtElec_ULS2_mc_closure","ULS; p_{T} (GeV/c); Count",300,0,30);
        fPtElec_LS2_mc_closure = new TH1F("fPtElec_LS2_mc_closure","LS; p_{T} (GeV/c); Count",300,0,30);
        
        
    }
    
    fPtTrigger_Inc = new TH1F("fPtTrigger_Inc","pT dist for Hadron Contamination; p_{t} (GeV/c); Count",300,0,30);
    fTPCnsigma_pt_2D = new TH2F("fTPCnsigma_pt_2D",";pt (GeV/c);TPC Electron N#sigma",1000,0.3,30,1000,-15,10);
    
    fTPCnsigma_pt_2D0 = new TH2F("fTPCnsigma_pt_2D0",";pt (GeV/c);TPC Electron N#sigma",1000,0.3,30,1000,-15,10);
    fTPCnsigma_pt_2D1 = new TH2F("fTPCnsigma_pt_2D1",";pt (GeV/c);TPC Electron N#sigma",1000,0.3,30,1000,-15,10);
    fTPCnsigma_pt_2D2 = new TH2F("fTPCnsigma_pt_2D2",";pt (GeV/c);TPC Electron N#sigma",1000,0.3,30,1000,-15,10);
    fTPCnsigma_pt_2D3 = new TH2F("fTPCnsigma_pt_2D3",";pt (GeV/c);TPC Electron N#sigma",1000,0.3,30,1000,-15,10);
    fTPCnsigma_pt_2D4 = new TH2F("fTPCnsigma_pt_2D4",";pt (GeV/c);TPC Electron N#sigma",1000,0.3,30,1000,-15,10);
    fTPCnsigma_pt_2D5 = new TH2F("fTPCnsigma_pt_2D5",";pt (GeV/c);TPC Electron N#sigma",1000,0.3,30,1000,-15,10);
    
    
    
    //new histos for TPC signal -> Can be used for any p range
    fTPCnsigma_p_TPC = new TH2F("fTPCnsigma_p_TPC",";p (GeV/c);TPC Electron N#sigma",3000,0,30,1000,-15,10);
    fTPCnsigma_p_TPC_on_EMCal_acc = new TH2F("fTPCnsigma_p_TPC_on_EMCal_acc",";p (GeV/c);TPC Electron N#sigma",3000,0,30,1000,-15,10);
    fTPCnsigma_p_TPC_EoverP_cut = new TH2F("fTPCnsigma_p_TPC_EoverP_cut",";p (GeV/c);TPC Electron N#sigma",3000,0,30,1000,-15,10);
    
    fShowerShapeCut = new TH2F("fShowerShapeCut","Shower Shape;M02;M20",500,0,1.8,500,0,1.8);
    
    fOutputList->Add(fShowerShapeCut);
    
    fCharge_n = new TH1F("fCharge_n","Inclusive Electrons (Negative Charge); p_{t} (GeV/c); Count",200,0,30);
    fCharge_p = new TH1F("fCharge_p","Inclusive Positrons (Positive Charge); p_{t} (GeV/c); Count",200,0,30);
    
    
    fOutputList->Add(fCharge_n);
    fOutputList->Add(fCharge_p);
    
    if(fIsMC)
    {
        
        fEtaPhi_num=new TH2F("fEtaPhi_num","#eta x #phi track;#phi;#eta",200,0.,5,50,-1.,1.);
        fEtaPhi_den=new TH2F("fEtaPhi_den","#eta x #phi track;#phi;#eta",200,0.,5,50,-1.,1.);
        fEtaPhi_data=new TH2F("fEtaPhi_data","#eta x #phi track;#phi;#eta",200,0.,5,50,-1.,1.);
        
        fOutputList->Add(fEtaPhi_num);
        fOutputList->Add(fEtaPhi_den);
        fOutputList->Add(fEtaPhi_data);
        
        fpt_reco_pt_MC_num=new TH2F("fpt_reco_pt_MC_num","pt reco x pt MC;pt reco; pt MC",300,0.,30,300,0.,30);
        fpt_reco_pt_MC_den=new TH2F("fpt_reco_pt_MC_den","pt reco x pt MC;pt reco; pt MC",300,0.,30,300,0.,30);
        
        fOutputList->Add(fpt_reco_pt_MC_num);
        fOutputList->Add(fpt_reco_pt_MC_den);
        
        
    }
    
    fECluster_pure= new TH1F("fECluster_pure", ";ECluster pure",2000,0,100);
    fOutputList->Add(fECluster_pure);
    
    if(fUseEMCal)
    {
        
        fECluster_not_exotic= new TH1F("fECluster_not_exotic", ";ECluster not exotic - function ",2000,0,100);
        
        fECluster_not_exotic1= new TH1F("fECluster_not_exotic1", ";ECluster not exotic Ncells > E/3+1",2000,0,100);
        
        fECluster_not_exotic2= new TH1F("fECluster_not_exotic2", ";ECluster not exotic 2",2000,0,100);
        fECluster_exotic= new TH1F("fECluster_exotic", ";ECluster exotic",2000,0,100);
        
        fOutputList->Add(fECluster_not_exotic);
        fOutputList->Add(fECluster_not_exotic1);
        fOutputList->Add(fECluster_not_exotic2);
        fOutputList->Add(fECluster_exotic);
    }
    
    //not associated with tracks
    fNCluster_pure= new TH1F("fNCluster_pure", ";Number of clusters - pure",2000,-1,1999);
    fNCluster_pure_aod= new TH1F("fNCluster_pure_aod", ";Number of clusters - pure -aod",2000,-1,1999);
    
    if(fUseEMCal)
    {
        fNCluster_ECluster= new TH2F("fNCluster_ECluster", ";Number of clusters vs. Energy of Cluster",2000,-1,1999, 4000, -1, 1999);
        fNcells_energy= new TH2F("fNcells_energy", "all clusters;Number of cells;Energy of Cluster",100,0,100, 2000, 0, 100);
        fNcells_energy_elec_selected= new TH2F("fNcells_energy_elec_selected", "clusters for electrons on TPC;Number of cells;Energy of Cluster",100,0,100, 2000, 0, 100);
        fNcells_energy_not_exotic= new TH2F("fNcells_energy_not_exotic", "not exotic cluster;Number of cells;Energy of Cluster ",100,0,100, 2000, 0, 100);
        fNcells_Energy_highE0= new TH2F("fNcells_Energy_highE0", "High Energy;Energy of Cluster;Number of cells",2000,0,100, 100, 0, 100);
        fNcells_Energy_highE1= new TH2F("fNcells_Energy_highE1", "High Energy;Energy of Cluster;Number of cells",2000,0,100, 100, 0, 100);
        
        fOutputList->Add(fNCluster_ECluster);
        fOutputList->Add(fNcells_energy);
        fOutputList->Add(fNcells_energy_elec_selected);
        fOutputList->Add(fNcells_energy_not_exotic);
        fOutputList->Add(fNcells_Energy_highE0);
        fOutputList->Add(fNcells_Energy_highE1);
        
    }
    
    
    if(fUseEMCal){
        
        if(!fIsAOD){
            fTime = new TH2D("fTime","Cells Cluster Time; p_{T} (GeV/c); Time (s)",300,0,30,1000,1e-8,1e-5);
            fTime2 = new TH2D("fTime2","Cells Cluster Time;  p_{T} (GeV/c); Time (s)",300,0,30,1000,1e-8,1e-5);
        }
        
        ftimingEle = new TH2D("ftimingEle","Cluster Time;  p_{T} (GeV/c); Time (ns)",300,0,30,2000,-100,100);
        ftimingEle2 = new TH2D("ftimingEle2","Cluster Time;  p_{T} (GeV/c); Time (ns)",300,0,30,2000,-100,100);
        
        fShowerShape_ha = new TH2F("fShowerShape_ha","Shower Shape hadrons;M02;M20",500,0,1.8,500,0,1.8);
        fShowerShape_ele = new TH2F("fShowerShape_ele","Shower Shape electrons;M02;M20",500,0,1.8,500,0,1.8);
        
        fShowerShapeM02_EoverP = new TH2F("fShowerShapeM02_EoverP","Shower Shape;M02;E/p",500,0,1.8,500,0,1.8);
        fShowerShapeM20_EoverP = new TH2F("fShowerShapeM20_EoverP","Shower Shape;M20;E/p",500,0,1.8,500,0,1.8);
        
        if(fIsMC){
            fShowerShapeM02_EoverP_Ecut12_MC = new TH2F("fShowerShapeM02_EoverP_Ecut12_MC","Shower Shape E >12 GeV;M02;E/p",500,0,1.8,500,0,1.8);
            fShowerShapeM20_EoverP_Ecut12_MC = new TH2F("fShowerShapeM20_EoverP_Ecut12_MC","Shower Shape E >12 GeV;M20;E/p",500,0,1.8,500,0,1.8);
            
            fShowerShapeM02_EoverP_Ecut8_MC = new TH2F("fShowerShapeM02_EoverP_Ecut8_MC","Shower Shape E >8 GeV;M02;E/p",500,0,1.8,500,0,1.8);
            fShowerShapeM20_EoverP_Ecut8_MC = new TH2F("fShowerShapeM20_EoverP_Ecut8_MC","Shower Shape E >8 GeV;M20;E/p",500,0,1.8,500,0,1.8);
        }
        
        
        
        
    }
    
    fOutputList->Add(fTOF01);
    fOutputList->Add(fTOF02);
    fOutputList->Add(fTOF03);
    
    
    
    fOutputList->Add(fPtElec_Inc);
    fOutputList->Add(fPtElec_ULS);
    fOutputList->Add(fPtElec_LS);
    
    
    if(fIsMC)
    {
        fOutputList->Add(fPtElec_ULS_NoPid);
        fOutputList->Add(fPtElec_LS_NoPid);
        fOutputList->Add(fPtElec_ULS_MC);
        fOutputList->Add(fPtElec_ULS_MC_weight);
        
        fOutputList->Add(fPtPrim);
        fOutputList->Add(fPtSec);
        fOutputList->Add(fPtPrim2);
        fOutputList->Add(fPtSec2);
        
        fOutputList->Add(fPtElec_ULS_weight);
        fOutputList->Add(fPtElec_LS_weight);
        fOutputList->Add(fPtElec_ULS2);
        fOutputList->Add(fPtElec_LS2);
        fOutputList->Add(fPtElec_ULS2_weight);
        fOutputList->Add(fPtElec_LS2_weight);
        
        fOutputList->Add(fPtElec_ULS_mc_closure);
        fOutputList->Add(fPtElec_LS_mc_closure);
        fOutputList->Add(fPtElec_ULS2_mc_closure);
        fOutputList->Add(fPtElec_LS2_mc_closure);
    }
    
    
    fOutputList->Add(fPtTrigger_Inc);
    fOutputList->Add(fTPCnsigma_pt_2D);
    
    fOutputList->Add(fTPCnsigma_pt_2D0);
    fOutputList->Add(fTPCnsigma_pt_2D1);
    fOutputList->Add(fTPCnsigma_pt_2D2);
    fOutputList->Add(fTPCnsigma_pt_2D3);
    fOutputList->Add(fTPCnsigma_pt_2D4);
    fOutputList->Add(fTPCnsigma_pt_2D5);
    
    fOutputList->Add(fTPCnsigma_p_TPC);
    fOutputList->Add(fTPCnsigma_p_TPC_on_EMCal_acc);
    fOutputList->Add(fTPCnsigma_p_TPC_EoverP_cut);
    
    
    
    
    
    
    fOutputList->Add(fNCluster_pure);
    fOutputList->Add(fNCluster_pure_aod);
    
    
    
    
    
    
    if(fUseEMCal){
        
        if(!fIsAOD){
            fOutputList->Add(fTime);
            fOutputList->Add(fTime2);
            
        }
        
        fOutputList->Add(ftimingEle);
        fOutputList->Add(ftimingEle2);
        
        fOutputList->Add(fShowerShape_ha);
        fOutputList->Add(fShowerShape_ele);
        
        fOutputList->Add(fShowerShapeM02_EoverP);
        fOutputList->Add(fShowerShapeM20_EoverP);
        
        if(fIsMC){
            fOutputList->Add(fShowerShapeM02_EoverP_Ecut12_MC);
            fOutputList->Add(fShowerShapeM20_EoverP_Ecut12_MC);
            
            fOutputList->Add(fShowerShapeM02_EoverP_Ecut8_MC);
            fOutputList->Add(fShowerShapeM20_EoverP_Ecut8_MC);
        }
        
        
        
    }
    
    fVtxZ_new1= new  TH1F("fVtxZ_new1","fVtxZ_new1",4000, -50,50);
    fVtxZ_new2= new  TH1F("fVtxZ_new2","fVtxZ_new2",4000, -50,50);
    fVtxZ_new3= new  TH1F("fVtxZ_new3","fVtxZ_new3",4000, -50,50);
    fVtxZ_new4= new  TH1F("fVtxZ_new4","fVtxZ_new4",4000, -50,50);
    
    fzRes1= new  TH1F("fzRes1","fzRes1",4000, 0,1);
    fzRes2= new  TH1F("fzRes2","fzRes2",4000, 0,1);
    fSPD_track_vtx1= new  TH1F("fSPD_track_vtx1","fSPD_track_vtx1",4000, -5,5);
    fSPD_track_vtx2= new  TH1F("fSPD_track_vtx2","fSPD_track_vtx2",4000, -5,5);
    
    
    //General Histograms
    
    //Steps
    //Step 1: Before Track cuts
    //Step 2: Before PID
    //Step 3: After PID
    
    fEoverP_pt = new TH2F *[3];
    fTPC_p = new TH2F *[3];
    fTPCnsigma_p = new TH2F *[3];
    
    fTPCnsigma_p_eta1 = new TH2F *[3];
    fTPCnsigma_p_eta2 = new TH2F *[3];
    fTPCnsigma_p_eta3 = new TH2F *[3];
    fTPCnsigma_p_eta4 = new TH2F *[3];
    
    fECluster= new TH1F *[3];
    fEtaPhi= new TH2F *[3];
    fVtxZ= new  TH1F *[3];
    fEtad= new  TH1F *[3];
    fNTracks= new  TH1F *[3];
    
    fNTracks_pt= new  TH2F *[3];
    fNTracks_eta= new  TH2F *[3];
    fNTracks_phi= new  TH2F *[3];
    
    fNClusters= new TH1F *[3];
    fTPCNcls_EoverP= new TH2F *[3];
    fTPCNcls_pid=new TH2F *[4];
    
    if(fIsMC)
    {
        fEoverP_pt_true_electrons = new TH2F("fEoverP_pt_true_electrons",";p_{T} (GeV/c);E/p ",1000,0,30,2000,0,2);
        fOutputList->Add(fEoverP_pt_true_electrons);
        
        fEoverP_pt_true_electrons_weight = new TH2F("fEoverP_pt_true_electrons_weight",";p_{T} (GeV/c);E/p ",1000,0,30,2000,0,2);
        fOutputList->Add(fEoverP_pt_true_electrons_weight);
        
        fEoverP_pt_true_HFE = new TH2F("fEoverP_pt_true_HFE",";p_{T} (GeV/c);E/p ",1000,0,30,2000,0,2);
        fOutputList->Add(fEoverP_pt_true_HFE);
        
        fEoverP_pt_not_HFE = new TH2F("fEoverP_pt_not_HFE",";p_{T} (GeV/c);E/p ",1000,0,30,2000,0,2);
        fOutputList->Add(fEoverP_pt_not_HFE);
        
        
    }
    
    fEoverP_ntracks_matched = new TH2F("fEoverP_ntracks_matched",";E/p (GeV/c);NTracksMatched ",200,0,2,10,0,10);
    fOutputList->Add(fEoverP_ntracks_matched);
    
    fEoverP_ncells = new TH2F("fEoverP_ncells",";E/p (GeV/c);ncells ",1000,0,10,30,0,30);
    fOutputList->Add(fEoverP_ncells);
    
    if(fIsMC)
    {
        fEmc_Ereco_gamma0 = new TH2F("fEmc_Ereco_gamma0",";E_{reco} (GeV);E_{true} (GeV)",300,0,30,300,0,30);
        fOutputList->Add(fEmc_Ereco_gamma0);
        fEmc_Ereco_gamma_ratio0 = new TH2F("fEmc_Ereco_gamma_ratio0",";E_{reco} (GeV);E_{true}/E_{reco}",300,0,30,200,0,2);
        fOutputList->Add(fEmc_Ereco_gamma_ratio0);
        
        fEmc_Ereco_gamma_all = new TH2F("fEmc_Ereco_gamma_all",";E_{reco} (GeV);E_{true} (GeV)",300,0,30,300,0,30);
        fOutputList->Add(fEmc_Ereco_gamma_all);
        fEmc_Ereco_gamma_ratio_all = new TH2F("fEmc_Ereco_gamma_ratio_all",";E_{reco} (GeV);E_{true}/E_{reco}",300,0,30,200,0,2);
        fOutputList->Add(fEmc_Ereco_gamma_ratio_all);
        
        
        fEmc_Ereco_gamma = new TH2F("fEmc_Ereco_gamma",";E_{reco} (GeV);E_{true} (GeV)",300,0,30,300,0,30);
        fOutputList->Add(fEmc_Ereco_gamma);
        fEmc_Ereco_gamma_ratio = new TH2F("fEmc_Ereco_gamma_ratio",";E_{reco} (GeV);E_{true}/E_{reco}",300,0,30,200,0,2);
        fOutputList->Add(fEmc_Ereco_gamma_ratio);
    }
    
    fEmc= new TH1F("fEmc",";E_{true}",1000,0,100);
    fOutputList->Add(fEmc);
    
    fEmc_Ereco_ele0 = new TH2F("fEmc_Ereco_ele0",";E_{reco} (GeV);E_{true} (GeV)",300,0,30,300,0,30);
    fOutputList->Add(fEmc_Ereco_ele0);
    fEmc_Ereco_ele_ratio0 = new TH2F("fEmc_Ereco_ele_ratio0",";E_{reco} (GeV);E_{true}/E_{reco}",300,0,30,200,0,2);
    fOutputList->Add(fEmc_Ereco_ele_ratio0);
    
    if(fIsMC)
    {
        fEmc_Ereco_ele_all = new TH2F("fEmc_Ereco_ele_all",";E_{reco} (GeV);E_{true} (GeV)",300,0,30,300,0,30);
        fOutputList->Add(fEmc_Ereco_ele_all);
        fEmc_Ereco_ele_ratio_all = new TH2F("fEmc_Ereco_ele_ratio_all",";E_{reco} (GeV);E_{true}/E_{reco}",300,0,30,200,0,2);
        fOutputList->Add(fEmc_Ereco_ele_ratio_all);
    }
    
    
    fEmc_Ereco_ele = new TH2F("fEmc_Ereco_ele",";E_{reco} (GeV);E_{true} (GeV)",300,0,30,300,0,30);
    fOutputList->Add(fEmc_Ereco_ele);
    fEmc_Ereco_ele_ratio = new TH2F("fEmc_Ereco_ele_ratio",";E_{reco} (GeV);E_{true}/E_{reco}",300,0,30,200,0,2);
    fOutputList->Add(fEmc_Ereco_ele_ratio);
    
    
    fEmc_Ereco_ratio_large_EoverP0 = new TH2F("fEmc_Ereco_ratio_large_EoverP0",";E_{reco} (GeV);E_{true}/E_{reco}",300,0,30,200,0,2);
    fOutputList->Add(fEmc_Ereco_ratio_large_EoverP0);
    
    fEmc_Ereco_ratio_small_EoverP0 = new TH2F("fEmc_Ereco_ratio_small_EoverP0",";E_{reco} (GeV);E_{true}/E_{reco}",300,0,30,200,0,2);
    fOutputList->Add(fEmc_Ereco_ratio_small_EoverP0);
    
    fEmc_Ereco_ratio_large_EoverP = new TH2F("fEmc_Ereco_ratio_large_EoverP",";E_{reco} (GeV);E_{true}/E_{reco}",300,0,30,200,0,2);
    fOutputList->Add(fEmc_Ereco_ratio_large_EoverP);
    
    fEmc_Ereco_ratio_small_EoverP = new TH2F("fEmc_Ereco_ratio_small_EoverP",";E_{reco} (GeV);E_{true}/E_{reco}",300,0,30,200,0,2);
    fOutputList->Add(fEmc_Ereco_ratio_small_EoverP);
    
    
    
    fEoverP_pt_true_electrons0 = new TH2F("fEoverP_pt_true_electrons0",";p_{T} (GeV/c);E/p ",1000,0,30,2000,0,2);
    fOutputList->Add(fEoverP_pt_true_electrons0);
    
    fEoverP_pt_true_hadrons = new TH2F("fEoverP_pt_true_hadrons",";p_{T} (GeV/c);E/p ",1000,0,30,2000,0,2);
    fOutputList->Add(fEoverP_pt_true_hadrons);
    fEoverP_pt_true_hadrons0 = new TH2F("fEoverP_pt_true_hadrons0",";p_{T} (GeV/c);E/p ",1000,0,30,2000,0,2);
    fOutputList->Add(fEoverP_pt_true_hadrons0);
    
    fEtaPhi_large_EoverP= new TH2F("fEtaPhi_large_EoverP","#eta x #phi Clusters;#phi;#eta",200,0.,5,50,-1.,1.);
    fEtaPhi_small_EoverP= new TH2F("fEtaPhi_small_EoverP","#eta x #phi Clusters;#phi;#eta",200,0.,5,50,-1.,1.);
    fOutputList->Add(fEtaPhi_large_EoverP);
    fOutputList->Add(fEtaPhi_small_EoverP);
    
    
    if(fIsMC)
    {
        fEoverP_pt_highE0 = new TH2F("fEoverP_pt_highE0",";p_{t} (GeV/c);E / p ",1000,0,30,2000,0,2);
        fEoverP_pt_highE1 = new TH2F("fEoverP_pt_highE1",";p_{t} (GeV/c);E / p ",1000,0,30,2000,0,2);
        fOutputList->Add(fEoverP_pt_highE0);
        fOutputList->Add(fEoverP_pt_highE1);
    }
    
    fECluster_highpT0= new TH1F("fECluster_highpT0", ";ECluster",2000, 0,100);
    fECluster_highpT1= new TH1F("fECluster_highpT1", ";ECluster",2000, 0,100);
    fOutputList->Add(fECluster_highpT0);
    fOutputList->Add(fECluster_highpT1);
    
    fTPC_momentum = new TH2F("fTPC_momentum",";p (GeV/c);TPC dE/dx (a. u.)",1000,0,30,1000,-20,200);
    fTPC_eta = new TH2F("fTPC_eta",";eta (GeV/c);TPC dE/dx (a. u.)",4000,-2,2,1000,-20,200);
    fOutputList->Add(fTPC_momentum);
    fOutputList->Add(fTPC_eta);
    fTPC_momentum1 = new TH2F("fTPC_momentum1",";p (GeV/c);TPC dE/dx (a. u.)",1000,0,30,1000,-20,200);
    fTPC_eta1 = new TH2F("fTPC_eta1",";eta (GeV/c);TPC dE/dx (a. u.)",4000,-2,2,1000,-20,200);
    fOutputList->Add(fTPC_momentum1);
    fOutputList->Add(fTPC_eta1);
    
    for(Int_t i = 0; i < 3; i++)
    {
        fEoverP_pt[i] = new TH2F(Form("fEoverP_pt%d",i),";p_{t} (GeV/c);E / p ",1000,0,30,2000,0,2);
        fTPC_p[i] = new TH2F(Form("fTPC_p%d",i),";pt (GeV/c);TPC dE/dx (a. u.)",1000,0.3,15,1000,-20,200);
        fTPCnsigma_p[i] = new TH2F(Form("fTPCnsigma_p%d",i),";p (GeV/c);TPC Electron N#sigma",1000,0.3,15,1000,-15,10);
        
        fTPCnsigma_p_eta1[i] = new TH2F(Form("fTPCnsigma_p_eta1%d",i),"-0.6 < #eta < -0.3;p (GeV/c);TPC Electron N#sigma",1000,0.3,15,1000,-15,10);
        fTPCnsigma_p_eta2[i] = new TH2F(Form("fTPCnsigma_p_eta2%d",i),"-0.3 < #eta < 0;p (GeV/c);TPC Electron N#sigma",1000,0.3,15,1000,-15,10);
        fTPCnsigma_p_eta3[i] = new TH2F(Form("fTPCnsigma_p_eta3%d",i),"0 < #eta < 0.3;p (GeV/c);TPC Electron N#sigma",1000,0.3,15,1000,-15,10);
        fTPCnsigma_p_eta4[i] = new TH2F(Form("fTPCnsigma_p_eta4%d",i),"0.3 < #eta < 0.6;p (GeV/c);TPC Electron N#sigma",1000,0.3,15,1000,-15,10);
        
        
        fECluster[i]= new TH1F(Form("fECluster%d",i), ";ECluster",2000, 0,100);
        fEtaPhi[i]= new TH2F(Form("fEtaPhi%d",i),"#eta x #phi Clusters;#phi;#eta",500,0.,5,200,-1.,1.);
        fVtxZ[i]= new  TH1F(Form("fVtxZ%d",i),"VtxZ",1000, -50,50);
        fEtad[i]= new  TH1F(Form("fEtad%d",i),"Eta distribution",200, -1.2,1.2);
        fNTracks[i]= new  TH1F(Form("fNTracks%d",i),"NTracks",1000, 0,5000);
        
        fNTracks_pt[i]= new  TH2F(Form("fNTracks_pt%d",i),"NTracks vs. pt",1000, 0,5000, 1000, 0, 100);
        fNTracks_eta[i]= new  TH2F(Form("fNTracks_eta%d",i),"NTracks vs. pt",1000, 0,5000, 500, -1.0, 1.0);
        fNTracks_phi[i]= new  TH2F(Form("fNTracks_phi%d",i),"NTracks vs. pt",1000, 0,5000, 500, 0, 5.0);
        
        
        
        fNClusters[i]= new TH1F(Form("fNClusters%d",i),"fNClusters0",200, 0,100);
        fTPCNcls_EoverP[i]= new TH2F(Form("fTPCNcls_EoverP%d",i),"TPCNcls_EoverP",1000,0,200,200,0,2);
        
        
        fOutputList->Add(fEoverP_pt[i]);
        fOutputList->Add(fTPC_p[i]);
        fOutputList->Add(fTPCnsigma_p[i]);
        
        fOutputList->Add(fTPCnsigma_p_eta1[i]);
        fOutputList->Add(fTPCnsigma_p_eta2[i]);
        fOutputList->Add(fTPCnsigma_p_eta3[i]);
        fOutputList->Add(fTPCnsigma_p_eta4[i]);
        
        
        fOutputList->Add(fECluster[i]);
        fOutputList->Add(fEtaPhi[i]);
        fOutputList->Add(fVtxZ[i]);
        fOutputList->Add(fEtad[i]);
        fOutputList->Add(fNTracks[i]);
        
        fOutputList->Add(fNTracks_pt[i]);
        fOutputList->Add(fNTracks_eta[i]);
        fOutputList->Add(fNTracks_phi[i]);
        
        fOutputList->Add(fNClusters[i]);
        fOutputList->Add(fTPCNcls_EoverP[i]);
    }
    
    if(fIsMC)
    {
        fEtaPhi_highE0= new TH2F("fEtaPhi_highE0","#eta x #phi Clusters;#phi;#eta",200,0.,5,50,-1.,1.);
        fEtaPhi_highE1= new TH2F("fEtaPhi_highE1","#eta x #phi Clusters;#phi;#eta",200,0.,5,50,-1.,1.);
        fOutputList->Add(fEtaPhi_highE0);
        fOutputList->Add(fEtaPhi_highE1);
        
        
        fEta_highE0=new TH1F("fEta_highE0", ";#eta; counts",100, -0.1,0.1);
        fPhi_highE0=new TH1F("fPhi_highE0", ";#phi; counts", 100, -0.1,0.1);
        fR_highE0=new TH1F("fR_highE0", ";R; counts", 100, -0.1,0.1);
        
        fEta_highE1=new TH1F("fEta_highE1", ";#eta; counts",100, -0.1,0.1);
        fPhi_highE1=new TH1F("fPhi_highE1", ";#phi; counts", 100, -0.1,0.1);
        fR_highE1=new TH1F("fR_highE1", ";R; counts", 100, -0.1,0.1);
        
        fNCluster_highE0= new TH1F("fNCluster_highE0","fNClusters0 in high E",200, 0,100);
        fNCluster_highE1= new TH1F("fNCluster_highE1","fNClusters1 in high E",200, 0,100);
        
        
        
        fOutputList->Add(fEta_highE0);
        fOutputList->Add(fEta_highE1);
        fOutputList->Add(fPhi_highE0);
        fOutputList->Add(fPhi_highE1);
        fOutputList->Add(fR_highE0);
        fOutputList->Add(fR_highE1);
        fOutputList->Add(fNCluster_highE0);
        fOutputList->Add(fNCluster_highE1);
    }
    
    
    fTrack_Multi= new  TH1F("fTrack_Multi","fTrack_Multi",1000, 0,1000);
    
    for(Int_t i = 0; i < 4; i++)
    {
        fTPCNcls_pid[i]= new TH2F(Form("fTPCNcls_pid%d",i),"fTPCNcls_pid;NCls;NCls for PID",200,0,200,200,0,200);
        fOutputList->Add(fTPCNcls_pid[i]);
    }
    
    //pt bin
    
    Double_t fPtBin[7];
    if(fUseAlternativeBinnig)
    {
        fPtBin[0] = 0.5;
        fPtBin[1] = 1.0;
        fPtBin[2] = 1.5;
        fPtBin[3] = 2.0;
        fPtBin[4] = 3.0;
        fPtBin[5] = 4.0;
        fPtBin[6] = 6.0;
    }
    else
    {
        fPtBin[0] = 1.0;
        fPtBin[1] = 2.0;
        fPtBin[2] = 4.0;
        fPtBin[3] = 6.0;
        fPtBin[4] = 8.0;
        fPtBin[5] = 10.0;
        fPtBin[6] = 15.0;
    }
    
    fEoverP_tpc = new TH2F *[6];
    fTPC_pt = new TH1F *[6];
    fTPCnsigma_pt = new TH1F *[6];
    
    fEta=new TH1F *[6];
    fPhi=new TH1F *[6];
    fR=new TH1F *[6];
    fR_EoverP=new TH2F *[6];
    fNcells=new TH1F *[6];
    fNcells_EoverP=new TH2F *[6];
    fM02_EoverP= new TH2F *[6];
    fM20_EoverP= new TH2F *[6];
    fEoverP_ptbins=new TH1F *[6];
    fECluster_ptbins=new TH1F *[6];
    fEoverP_wSSCut=new TH1F *[6];
    fNcells_electrons=new TH1F *[6];
    fNcells_hadrons=new TH1F *[6];
    fTPCnsigma_eta_electrons=new TH2F *[6];
    fTPCnsigma_eta_hadrons=new TH2F *[6];
    
    if(fCorrelationFlag)
    {
        fCEtaPhi_Inc = new TH2F *[6];
        
        fCEtaPhi_Inc_DiHadron = new TH2F *[6];
        
        fCEtaPhi_ULS = new TH2F *[6];
        fCEtaPhi_LS = new TH2F *[6];
        fCEtaPhi_ULS_NoP = new TH2F *[6];
        fCEtaPhi_LS_NoP = new TH2F *[6];
        
        fCEtaPhi_ULS_Weight = new TH2F *[6];
        fCEtaPhi_LS_Weight = new TH2F *[6];
        fCEtaPhi_ULS_NoP_Weight = new TH2F *[6];
        fCEtaPhi_LS_NoP_Weight = new TH2F *[6];
        
        if(fIsMC)
        {
            fCEtaPhi_Inc_PH = new TH2F *[6];
            fCEtaPhi_ULS_Weight_PH= new TH2F *[6];
            fCEtaPhi_LS_Weight_PH = new TH2F *[6];
            fCEtaPhi_ULS_NoP_Weight_PH = new TH2F *[6];
            fCEtaPhi_LS_NoP_Weight_PH = new TH2F *[6];
        }
        
        
        fCEtaPhi_Inc_EM = new TH2F *[6];
        
        fCEtaPhi_ULS_EM = new TH2F *[6];
        fCEtaPhi_LS_EM = new TH2F *[6];
        
        fCEtaPhi_ULS_Weight_EM = new TH2F *[6];
        fCEtaPhi_LS_Weight_EM = new TH2F *[6];
        
        fInvMass = new TH1F("fInvMass","",200,0,0.3);
        fInvMassBack = new TH1F("fInvMassBack","",200,0,0.3);
        fDCA = new TH1F("fDCA","",200,0,1);
        fDCABack = new TH1F("fDCABack","",200,0,1);
        fOpAngle = new TH1F("fOpAngle","",200,0,0.5);
        fOpAngleBack = new TH1F("fOpAngleBack","",200,0,0.5);
        
        //15-09 Add new Histograms for MC clousure test
        if(fIsMC)
        {
            fCetaPhi_MC_with_partner_greater = new TH2F *[6];
            fCetaPhi_MC_with_partner_below = new TH2F *[6];
            fCetaPhi_MC_NHFE_1partner_reco = new TH2F *[6];
            fpT_MC_with_partner_greater = new TH1F("fpT_MC_with_partner_greater","Two reco; p_{T} (GeV/c); Count",300,0,30);
            fpT_MC_with_partner_below = new TH1F("fpT_MC_with_partner_below","One reco; p_{T} (GeV/c); Count",300,0,30);
            fOutputList->Add(fpT_MC_with_partner_greater);
            fOutputList->Add(fpT_MC_with_partner_below);
        }
        
        //end of new Histograms
        
        fOutputList->Add(fInvMass);
        fOutputList->Add(fInvMassBack);
        fOutputList->Add(fDCA);
        fOutputList->Add(fDCABack);
        fOutputList->Add(fOpAngle);
        fOutputList->Add(fOpAngleBack);
    }
    
    if(fFillBackground){
        
        fInvMass = new TH1F("fInvMass","",5000,0,5);
        fInvMassBack = new TH1F("fInvMassBack","",5000,0,5);
        fDCA = new TH1F("fDCA","",200,0,1);
        fDCABack = new TH1F("fDCABack","",200,0,1);
        fOpAngle = new TH1F("fOpAngle","",200,0,0.5);
        fOpAngleBack = new TH1F("fOpAngleBack","",200,0,0.5);
        
        fOutputList->Add(fInvMass);
        fOutputList->Add(fInvMassBack);
        fOutputList->Add(fDCA);
        fOutputList->Add(fDCABack);
        fOutputList->Add(fOpAngle);
        fOutputList->Add(fOpAngleBack);
        
        //histos for TPC-only
        fInvMass2 = new TH1F("fInvMass2","",5000,0,5);
        fInvMassBack2 = new TH1F("fInvMassBack2","",5000,0,5);
        
        fInvMass2_weight = new TH1F("fInvMass2_weight","",5000,0,5);
        fInvMassBack2_weight = new TH1F("fInvMassBack2_weight","",5000,0,5);
        
        fDCA2 = new TH1F("fDCA2","",200,0,1);
        fDCABack2 = new TH1F("fDCABack2","",200,0,1);
        fOpAngle2 = new TH1F("fOpAngle2","",200,0,0.5);
        fOpAngleBack2 = new TH1F("fOpAngleBack2","",200,0,0.5);
        
        fOutputList->Add(fInvMass2);
        fOutputList->Add(fInvMassBack2);
        
        fOutputList->Add(fInvMass2_weight);
        fOutputList->Add(fInvMassBack2_weight);
        
        
        fOutputList->Add(fDCA2);
        fOutputList->Add(fDCABack2);
        fOutputList->Add(fOpAngle2);
        fOutputList->Add(fOpAngleBack2);
        
    }
    
    //new histo for trigger data
    if (fUseTrigger)
    {
        Double_t fPtBin_trigger[11] = {1,2,4,6,8,10,12,14,16,18,20};
        fEoverP_tpc_p_trigger = new TH2F *[10];
        fEoverP_tpc_pt_trigger = new TH2F *[10];
        fInvMass_pT = new TH1F *[10];
        fInvMassBack_pT = new TH1F *[10];
        
        for(Int_t i = 0; i < 10; i++)
        {
            fEoverP_tpc_pt_trigger[i] = new TH2F(Form("fEoverP_tpc_pt_trigger%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;TPC Electron N#sigma; E/p ",fPtBin_trigger[i],fPtBin_trigger[i+1]),1000,-15,15,100,0,2);
            fEoverP_tpc_p_trigger[i] = new TH2F(Form("fEoverP_tpc_p_trigger%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;TPC Electron N#sigma; E/p ",fPtBin_trigger[i],fPtBin_trigger[i+1]),1000,-15,15,100,0,2);
            
            
            fInvMass_pT[i] = new TH1F(Form("fInvMass_pT%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c; Mass ; Counts",fPtBin_trigger[i],fPtBin_trigger[i+1]),5000,0,5);
            fInvMassBack_pT[i] = new TH1F(Form("fInvMassBack_pT%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;Mass;Counts",fPtBin_trigger[i],fPtBin_trigger[i+1]),5000,0,5);
            
            fOutputList->Add(fEoverP_tpc_pt_trigger[i]);
            fOutputList->Add(fEoverP_tpc_p_trigger[i]);
            fOutputList->Add(fInvMass_pT[i]);
            fOutputList->Add(fInvMassBack_pT[i]);
            
        }
    }
    
    
    for(Int_t i = 0; i < 6; i++)
    {
        fEoverP_tpc[i] = new TH2F(Form("fEoverP_tpc%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;TPC Electron N#sigma;E / p ",fPtBin[i],fPtBin[i+1]),1000,-15,15,100,0,2);
        fTPC_pt[i] = new TH1F(Form("fTPC_pt%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;TPC Electron N#sigma;Count",fPtBin[i],fPtBin[i+1]),200,20,200);
        fTPCnsigma_pt[i] = new TH1F(Form("fTPCnsigma_pt%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;TPC Electron N#sigma;Count",fPtBin[i],fPtBin[i+1]),200,-15,10);
        
        fEta[i]=new TH1F(Form("fEta%d",i), Form("%1.1f < p_{t} < %1.1f GeV/c;#eta; counts",fPtBin[i],fPtBin[i+1]),100, -0.1,0.1);
        fPhi[i]=new TH1F(Form("fPhi%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#phi; counts )",fPtBin[i],fPtBin[i+1]), 100, -0.1,0.1);
        fR[i]=new TH1F(Form("fR%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;R;counts )",fPtBin[i],fPtBin[i+1]), 100, -0.1,0.1);
        fR_EoverP[i]=new TH2F(Form("fR_EoverP%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;R;E / p ",fPtBin[i],fPtBin[i+1]),100, 0,0.1,1000,0,10);
        fNcells[i]=new TH1F(Form("fNcells%d",i), Form("%1.1f < p_{t} < %1.1f GeV/c;ncells;counts ",fPtBin[i],fPtBin[i+1]),100, 0, 30);
        fNcells_electrons[i]=new TH1F(Form("fNcells_electrons%d",i), Form("%1.1f < p_{t} < %1.1f GeV/c;ncells;counts ",fPtBin[i],fPtBin[i+1]),100, 0, 30);
        fNcells_hadrons[i]=new TH1F(Form("fNcells_hadrons%d",i), Form("%1.1f < p_{t} < %1.1f GeV/c;ncells;counts ",fPtBin[i],fPtBin[i+1]),100, 0, 30);
        fNcells_EoverP[i]=new TH2F(Form("fNcells_EoverP%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c; E/p; Ncells ",fPtBin[i],fPtBin[i+1]),1000, 0,10,30,0,30);
        fM02_EoverP[i]= new TH2F(Form("fM02_EoverP%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c; M02; E / p ",fPtBin[i],fPtBin[i+1]),1000,0,100,100,0,2);
        fM20_EoverP[i]= new TH2F(Form("fM20_EoverP%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c; M20; E / p ",fPtBin[i],fPtBin[i+1]),1000,0,100,100,0,2);
        fEoverP_ptbins[i] = new TH1F(Form("fEoverP_ptbins%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;E / p ",fPtBin[i],fPtBin[i+1]),500,0,2);
        fECluster_ptbins[i]= new TH1F(Form("fECluster_ptbins%d",i), Form("%1.1f < p_{t} < %1.1f GeV/c;ECluster; Counts ",fPtBin[i],fPtBin[i+1]),2000, 0,100);
        fEoverP_wSSCut[i]=new TH1F(Form("fEoverP_wSSCut%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;E / p ; Counts",fPtBin[i],fPtBin[i+1]),500,0,2);
        fTPCnsigma_eta_electrons[i]=new TH2F(Form("fTPCnsigma_eta_electrons%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;TPC Electron N#sigma;Eta ",fPtBin[i],fPtBin[i+1]),1000,-15,15,100,-1,1);
        fTPCnsigma_eta_hadrons[i]=new TH2F(Form("fTPCnsigma_eta_hadrons%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;TPC Electron N#sigma;Eta ",fPtBin[i],fPtBin[i+1]),1000,-15,15,100,-1,1);
        
        fOutputList->Add(fEoverP_tpc[i]);
        fOutputList->Add(fTPC_pt[i]);
        fOutputList->Add(fTPCnsigma_pt[i]);
        
        fOutputList->Add(fEta[i]);
        fOutputList->Add(fPhi[i]);
        fOutputList->Add(fR[i]);
        fOutputList->Add(fR_EoverP[i]);
        fOutputList->Add(fNcells[i]);
        fOutputList->Add(fNcells_electrons[i]);
        fOutputList->Add(fNcells_hadrons[i]);
        fOutputList->Add(fNcells_EoverP[i]);
        fOutputList->Add(fECluster_ptbins[i]);
        fOutputList->Add(fEoverP_ptbins[i]);
        fOutputList->Add(fEoverP_wSSCut[i]);
        fOutputList->Add(fM02_EoverP[i]);
        fOutputList->Add(fM20_EoverP[i]);
        fOutputList->Add(fTPCnsigma_eta_electrons[i]);
        fOutputList->Add(fTPCnsigma_eta_hadrons[i]);
        
        
        
        if(fCorrelationFlag)
        {
            //new Histograms 15-09
            if (fIsMC)
            {
                fCetaPhi_MC_with_partner_greater[i] = new TH2F(Form("fCetaPhi_MC_with_partner_greater%d",i),Form("MC with part greather than min %1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
                fCetaPhi_MC_with_partner_below[i] = new TH2F(Form("fCetaPhi_MC_with_partner_below%d",i),Form("MC with part below than min %1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
                
                fCetaPhi_MC_NHFE_1partner_reco[i] = new TH2F(Form("fCetaPhi_MC_NHFE_1partner_reco%d",i),Form("MC with part below than min %1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
                
                fOutputList->Add(fCetaPhi_MC_with_partner_greater[i]);
                fOutputList->Add(fCetaPhi_MC_with_partner_below[i]);
                fOutputList->Add(fCetaPhi_MC_NHFE_1partner_reco[i]);
                
                //end new histograms
            }
            fCEtaPhi_Inc[i] = new TH2F(Form("fCEtaPhi_Inc%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
            fCEtaPhi_Inc_DiHadron[i] = new TH2F(Form("fCEtaPhi_Inc_DiHadron%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
            
            fCEtaPhi_ULS[i] = new TH2F(Form("fCEtaPhi_ULS%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
            fCEtaPhi_LS[i] = new TH2F(Form("fCEtaPhi_LS%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
            fCEtaPhi_ULS_NoP[i] = new TH2F(Form("fCEtaPhi_ULS_NoP%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
            fCEtaPhi_LS_NoP[i] = new TH2F(Form("fCEtaPhi_LS_NoP%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
            
            fCEtaPhi_ULS_Weight[i] = new TH2F(Form("fCEtaPhi_ULS_Weight%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
            fCEtaPhi_LS_Weight[i] = new TH2F(Form("fCEtaPhi_LS_Weight%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
            fCEtaPhi_ULS_NoP_Weight[i] = new TH2F(Form("fCEtaPhi_ULS_NoP_Weight%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
            fCEtaPhi_LS_NoP_Weight[i] = new TH2F(Form("fCEtaPhi_LS_NoP_Weight%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
            
            if(fIsMC)
            {
                fCEtaPhi_Inc_PH[i] = new TH2F(Form("fCEtaPhi_Inc_PH%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
                fCEtaPhi_ULS_Weight_PH[i] = new TH2F(Form("fCEtaPhi_ULS_Weight_PH%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
                fCEtaPhi_LS_Weight_PH[i] = new TH2F(Form("fCEtaPhi_LS_Weight_PH%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
                fCEtaPhi_ULS_NoP_Weight_PH[i] = new TH2F(Form("fCEtaPhi_ULS_NoP_Weight_PH%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
                fCEtaPhi_LS_NoP_Weight_PH[i] = new TH2F(Form("fCEtaPhi_LS_NoP_Weight_PH%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
                fOutputList->Add(fCEtaPhi_Inc_PH[i]);
                fOutputList->Add(fCEtaPhi_ULS_Weight_PH[i]);
                fOutputList->Add(fCEtaPhi_LS_Weight_PH[i]);
                fOutputList->Add(fCEtaPhi_ULS_NoP_Weight_PH[i]);
                fOutputList->Add(fCEtaPhi_LS_NoP_Weight_PH[i]);
            }
            
            
            
            fOutputList->Add(fCEtaPhi_Inc[i]);
            fOutputList->Add(fCEtaPhi_Inc_DiHadron[i]);
            
            fOutputList->Add(fCEtaPhi_ULS[i]);
            fOutputList->Add(fCEtaPhi_LS[i]);
            fOutputList->Add(fCEtaPhi_ULS_NoP[i]);
            fOutputList->Add(fCEtaPhi_LS_NoP[i]);
            
            fOutputList->Add(fCEtaPhi_ULS_Weight[i]);
            fOutputList->Add(fCEtaPhi_LS_Weight[i]);
            fOutputList->Add(fCEtaPhi_ULS_NoP_Weight[i]);
            fOutputList->Add(fCEtaPhi_LS_NoP_Weight[i]);
            
            
            
            if(fEventMixingFlag)
            {
                fCEtaPhi_Inc_EM[i] = new TH2F(Form("fCEtaPhi_Inc_EM%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
                
                fCEtaPhi_ULS_EM[i] = new TH2F(Form("fCEtaPhi_ULS_EM%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
                fCEtaPhi_LS_EM[i] = new TH2F(Form("fCEtaPhi_LS_EM%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
                
                fCEtaPhi_ULS_Weight_EM[i] = new TH2F(Form("fCEtaPhi_ULS_Weight_EM%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
                fCEtaPhi_LS_Weight_EM[i] = new TH2F(Form("fCEtaPhi_LS_Weight_EM%d",i),Form("%1.1f < p_{t} < %1.1f GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
                
                fOutputList->Add(fCEtaPhi_Inc_EM[i]);
                
                fOutputList->Add(fCEtaPhi_ULS_EM[i]);
                fOutputList->Add(fCEtaPhi_LS_EM[i]);
                
                fOutputList->Add(fCEtaPhi_ULS_Weight_EM[i]);
                fOutputList->Add(fCEtaPhi_LS_Weight_EM[i]);
            }
        }
    }
    
    //pt integrated
    fTPCnsigma_eta = new TH2F("fTPCnsigma_eta",";Pseudorapidity #eta; TPC signal - <TPC signal>_{elec} (#sigma)",200,-0.9,0.9,200,-15,15);
    fTPCnsigma_phi = new TH2F("fTPCnsigma_phi",";Azimuthal Angle #phi; TPC signal - <TPC signal>_{elec} (#sigma)",200,0,2*TMath::Pi(),200,-15,15);
    
    
    
    fNcells_pt=new TH2F("fNcells_pt","fNcells_pt",1000, 0,20,100,0,30);
    fEoverP_pt_pions= new TH2F("fEoverP_pt_pions","fEoverP_pt_pions",1000,0,30,2000,0,2);
    
    ftpc_p_EoverPcut= new TH2F("ftpc_p_EoverPcut","ftpc_p_EoverPcut",1000,0,30,200,20,200);
    fnsigma_p_EoverPcut= new TH2F("fnsigma_p_EoverPcut","fnsigma_p_EoverPcut",1000,0,30,500,-15,15);
    
    fEoverP_pt_pions2= new TH2F("fEoverP_pt_pions2","fEoverP_pt_pions2",1000,0,30,2000,0,2);
    fEoverP_pt_pions3= new TH2F("fEoverP_pt_pions3","fEoverP_pt_pions3",1000,0,30,2000,0,2);
    
    fEoverP_pt_pions2_highE1= new TH2F("fEoverP_pt_pions2_highE0","fEoverP_pt_pions2",1000,0,30,2000,0,2);
    fEoverP_pt_pions2_highE0= new TH2F("fEoverP_pt_pions2_highE1","fEoverP_pt_pions2",1000,0,30,2000,0,2);
    
    fEoverP_pt_hadrons= new TH2F("fEoverP_pt_hadrons","fEoverP_pt_hadrons",1000,0,30,2000,0,2);
    
    
    fOutputList->Add(fTPCnsigma_eta);
    fOutputList->Add(fTPCnsigma_phi);
    
    fOutputList->Add(fNcells_pt);
    fOutputList->Add(fEoverP_pt_pions);
    
    fOutputList->Add(ftpc_p_EoverPcut);
    fOutputList->Add(fnsigma_p_EoverPcut);
    
    fOutputList->Add(fEoverP_pt_pions2);
    fOutputList->Add(fEoverP_pt_pions3);
    
    fOutputList->Add(fEoverP_pt_pions2_highE0);
    fOutputList->Add(fEoverP_pt_pions2_highE1);
    
    fOutputList->Add(fEoverP_pt_hadrons);
    
    fOutputList->Add(fVtxZ_new1);
    fOutputList->Add(fVtxZ_new2);
    fOutputList->Add(fVtxZ_new3);
    fOutputList->Add(fVtxZ_new4);
    
    fOutputList->Add(fzRes1);
    fOutputList->Add(fzRes2);
    fOutputList->Add(fSPD_track_vtx1);
    fOutputList->Add(fSPD_track_vtx2);
    
    
    
    //__________________________________________________________________
    //Efficiency studies
    if(fIsMC)
    {
        fPtBackgroundBeforeReco = new TH1F("fPtBackgroundBeforeReco",";p_{T} (GeV/c);Count",300,0,30);
        fPtBackgroundBeforeReco_weight = new TH1F("fPtBackgroundBeforeReco_weight",";p_{T} (GeV/c);Count",300,0,30);
        if(fFillBackground)fPtBackgroundBeforeReco2 = new TH1F("fPtBackgroundBeforeReco2",";p_{T} (GeV/c);Count",300,0,30);
        if(fFillBackground)fPtBackgroundBeforeReco2_weight = new TH1F("fPtBackgroundBeforeReco2_weight",";p_{T} (GeV/c);Count",300,0,30);
        fpT_m_electron= new TH2F("fpT_m_electron","fpT_m_electron",1000,0,100,1000,0,100);
        fpT_gm_electron= new TH2F("fpT_gm_electron","fpT_gm_electron",1000,0,100,1000,0,100);
        
        fPtBackgroundAfterReco = new TH1F("fPtBackgroundAfterReco",";p_{T} (GeV/c);Count",300,0,30);
        fPtMCparticleAll = new TH1F("fPtMCparticleAll",";p_{T} (GeV/c);Count",200,0,40);
        fPtMCparticleReco = new TH1F("fPtMCparticleReco",";p_{T} (GeV/c);Count",200,0,40);
        
        fPtMCparticleAll_nonPrimary = new TH1F("fPtMCparticleAll_nonPrimary",";p_{T} (GeV/c);Count",200,0,40);
        fPtMCparticleAlle_nonPrimary = new TH1F("fPtMCparticleAlle_nonPrimary",";p_{T} (GeV/c);Count",200,0,40);
        fPtMCparticleAlle_Primary = new TH1F("fPtMCparticleAlle_Primary",";p_{T} (GeV/c);Count",200,0,40);
        
        fPtMCparticleReco_nonPrimary = new TH1F("fPtMCparticleReco_nonPrimary",";p_{T} (GeV/c);Count",200,0,40);
        
        fPtMCparticleAllHfe1 = new TH1F("fPtMCparticleAllHfe1",";p_{t} (GeV/c);Count",200,0,40);
        fPtMCparticleRecoHfe1 = new TH1F("fPtMCparticleRecoHfe1",";p_{t} (GeV/c);Count",200,0,40);
        fPtMCparticleAllHfe2 = new TH1F("fPtMCparticleAllHfe2",";p_{t} (GeV/c);Count",200,0,40);
        fPtMCparticleRecoHfe2 = new TH1F("fPtMCparticleRecoHfe2",";p_{t} (GeV/c);Count",200,0,40);
        
        fPtMCelectronAfterAll = new TH1F("fPtMCelectronAfterAll",";p_{T} (GeV/c);Count",200,0,40);
        fPtMCelectronAfterAll_unfolding = new TH1F("fPtMCelectronAfterAll_unfolding",";p_{T} (GeV/c);Count",200,0,40);
        fPtMCelectronAfterAll_nonPrimary = new TH1F("fPtMCelectronAfterAll_nonPrimary",";p_{T} (GeV/c);Count",200,0,40);
        fPtMCelectronAfterAll_Primary = new TH1F("fPtMCelectronAfterAll_Primary",";p_{T} (GeV/c);Count",200,0,40);
        
        
        
        fPtMCpi0 = new TH1F("fPtMCpi0",";p_{t} (GeV/c);Count",2000,0,100);
        fPtMCeta = new TH1F("fPtMCeta",";p_{T} (GeV/c);Count",2000,0,100);
        fPtMCpi02 = new TH1F("fPtMCpi02",";p_{t} (GeV/c);Count",2000,0,100);
        fPtMCeta2 = new TH1F("fPtMCeta2",";p_{T} (GeV/c);Count",2000,0,100);
        fPtMCpi03 = new TH1F("fPtMCpi03",";p_{t} (GeV/c);Count",2000,0,100);
        fPtMCeta3 = new TH1F("fPtMCeta3",";p_{T} (GeV/c);Count",2000,0,100);
        
        fPtMC_EMCal_All= new TH1F("fPtMC_EMCal_All",";p_{t} (GeV/c);Count",200,0,40);
        fPtMC_EMCal_Selected= new TH1F("fPtMC_EMCal_Selected",";p_{t} (GeV/c);Count",200,0,40);
        fPtMC_TPC_All= new TH1F("fPtMC_TPC_All",";p_{T} (GeV/c);Count",200,0,40);
        fPtMC_TPC_Selected = new TH1F("fPtMC_TPC_Selected",";p_{T} (GeV/c);Count",200,0,40);
        
        fPt_track_match_den = new TH1F("fPt_track_match_den",";p_{T} (GeV/c);Count",200,0,40);
        fPt_track_match_num = new TH1F("fPt_track_match_num",";p_{T} (GeV/c);Count",200,0,40);
        
        fPtMCWithLabel = new TH1F("fPtMCWithLabel",";p_{t} (GeV/c);Count",200,0,40);
        fPtMCWithoutLabel = new TH1F("fPtMCWithoutLabel",";p_{t} (GeV/c);Count",200,0,40);
        fPtIsPhysicaPrimary = new TH1F("fPtIsPhysicaPrimary",";p_{t} (GeV/c);Count",200,0,40);
        
        fMCEffPID_beforePID = new TH1F("fMCEffPID_beforePID",";p_{t} (GeV/c);Count",2000,0,100);
        fMCEffPID_afterPID = new TH1F("fMCEffPID_afterPID",";p_{t} (GeV/c);Count",2000,0,100);
        fpTMCGeneratedBKG = new TH1F("fpTMCGeneratedBKG",";p_{t} (GeV/c);Count",2000,0,100);
        
        fOutputList->Add(fMCEffPID_beforePID);
        fOutputList->Add(fMCEffPID_afterPID);
        fOutputList->Add(fpTMCGeneratedBKG);
        
        
        fOutputList->Add(fPtBackgroundBeforeReco);
        fOutputList->Add(fPtBackgroundBeforeReco_weight);
        
        if(fFillBackground) fOutputList->Add(fPtBackgroundBeforeReco2);
        if(fFillBackground) fOutputList->Add(fPtBackgroundBeforeReco2_weight);
        
        fOutputList->Add(fpT_m_electron);
        fOutputList->Add(fpT_gm_electron);
        
        fOutputList->Add(fPtBackgroundAfterReco);
        fOutputList->Add(fPtMCparticleAll);
        fOutputList->Add(fPtMCparticleReco);
        
        fOutputList->Add(fPtMCparticleAll_nonPrimary);
        fOutputList->Add(fPtMCparticleAlle_nonPrimary);
        
        fOutputList->Add(fPtMCparticleAlle_Primary);
        fOutputList->Add(fPtMCparticleReco_nonPrimary);
        
        fOutputList->Add(fPtMCparticleAllHfe1);
        fOutputList->Add(fPtMCparticleRecoHfe1);
        fOutputList->Add(fPtMCparticleAllHfe2);
        fOutputList->Add(fPtMCparticleRecoHfe2);
        fOutputList->Add(fPtMCelectronAfterAll);
        fOutputList->Add(fPtMCelectronAfterAll_unfolding);
        
        fOutputList->Add(fPtMCelectronAfterAll_nonPrimary);
        fOutputList->Add(fPtMCelectronAfterAll_Primary);
        
        
        
        fOutputList->Add(fPtMCpi0);
        fOutputList->Add(fPtMCeta);
        fOutputList->Add(fPtMCpi02);
        fOutputList->Add(fPtMCeta2);
        fOutputList->Add(fPtMCpi03);
        fOutputList->Add(fPtMCeta3);
        fOutputList->Add(fPtMC_EMCal_All);
        fOutputList->Add(fPtMC_EMCal_Selected);
        fOutputList->Add(fPtMC_TPC_All);
        fOutputList->Add(fPtMC_TPC_Selected);
        
        fOutputList->Add(fPt_track_match_den);
        fOutputList->Add(fPt_track_match_num);
        
        fOutputList->Add(fPtMCWithLabel);
        fOutputList->Add(fPtMCWithoutLabel);
        fOutputList->Add(fPtIsPhysicaPrimary);
    }
    
    if(!fIspp){
        fCentralityHist = new TH1F("fCentralityHist",";Centrality (%); Count",1000000,0,100);
        fCentralityHistPass = new TH1F("fCentralityHistPass",";Centrality (%); Count",1000000,0,100);
        fOutputList->Add(fCentralityHist);
        fOutputList->Add(fCentralityHistPass);
        
    }
    
    
    //______________________________________________________________________
    //Mixed event analysis -- it was removed from is pp because
    if(fEventMixingFlag)
    {
        fPoolNevents = new TH1F("fPoolNevents","Event Mixing Statistics; Number of events; Count",1000,0,1000);
        fOutputList->Add(fPoolNevents);
        
        Int_t trackDepth = 2000; // number of objects (tracks) kept per event buffer bin. Once the number of stored objects (tracks) is above that limit, the oldest ones are removed.
        Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
        
        Int_t nCentralityBins  = 15;
        Double_t centralityBins[] = { 0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1 };
        
        Int_t nZvtxBins  = 9;
        Double_t vertexBins[] = {-10, -7, -5, -3, -1, 1, 3, 5, 7, 10};
        
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
    
    //______________________________________________________________________
    //Vertex Selection
    if(!fIspp){
        
        fNevent2->Fill(8);
        
        //printf("\n\n =============================\n =============================\n fCalibrateTPC = %d  \n =============================\n=============================\n",fCalibrateTPC);
        //cout << "\n\n =============================\n =============================\n cout fCalibrateTPC = " << fCalibrateTPC << "\n =============================\n=============================\n"<< endl;
        
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
        
        
        
        //new line to apply the events cuts using the HFE package, cuts defined in the Conf. file
        //seems not work here... why?  cuts are done by hand for vertex and pileup
        /*
         fNevent2->Fill(10);
         
         if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fAOD)){
         return;
         }
         //from Andrea
         fCFM->SetRecEventInfo(fAOD);
         */
        
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
    
    //______________________________________________________________________
    //EMCal Trigger Selection (Threshold selection)
    
    TString firedTrigger;
    TString TriggerEG1("EG1"); //takes trigger with name with EG1, ex: CEMC7EG1-B-NOPF-CENTNOTRD
    TString TriggerEG2("EG2");
    //Jan 17, 2014
    //TString TriggerEJE("EJE");
    
    if(fAOD) firedTrigger = fAOD->GetFiredTriggerClasses();
    else if(fESD) firedTrigger = fESD->GetFiredTriggerClasses();
    
    //Bool_t IsEventEMCALL0=kTRUE;
    Bool_t IsEventEMCALL1=kFALSE;
    
    if(firedTrigger.Contains(TriggerEG1)){
        fNevent->Fill(2);
        IsEventEMCALL1=kTRUE;
    }
    if(firedTrigger.Contains(TriggerEG2)){
        fNevent->Fill(3);
        IsEventEMCALL1=kTRUE;
    }
    
    //if the flag is for a given threshold and it was not fired, return.
    
    if(fEMCEG1){
        if(!firedTrigger.Contains(TriggerEG1))return;
        if(firedTrigger.Contains(TriggerEG2)){
            fNevent->Fill(4);
            
        }
        
    }
    
    
    if(fEMCEG2){
        if(!firedTrigger.Contains(TriggerEG2))return;
        if(firedTrigger.Contains(TriggerEG1)){
            fNevent->Fill(5);
        }
        
    }
    
    
    
    //______________________________________________________________________
    //Testing if there is an overlap EGA and EJE
    //none
    /*
     if(!(firedTrigger.Contains(TriggerEG1) && firedTrigger.Contains(TriggerEG2) ) && !firedTrigger.Contains(TriggerEJE))
     {
     fNevent->Fill(6);
     }
     //only GA
     if((firedTrigger.Contains(TriggerEG1) || firedTrigger.Contains(TriggerEG2)) && !firedTrigger.Contains(TriggerEJE))
     {
     fNevent->Fill(7);
     }
     //only JE
     if(!(firedTrigger.Contains(TriggerEG1) && firedTrigger.Contains(TriggerEG2)) && firedTrigger.Contains(TriggerEJE))
     {
     fNevent->Fill(8);
     }
     //both
     if((firedTrigger.Contains(TriggerEG1) || firedTrigger.Contains(TriggerEG2)) && firedTrigger.Contains(TriggerEJE))
     {
     fNevent->Fill(9);
     }
     */
    
    
    
    
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
            Float_t centrality = -1;
            
            if(fIsAOD)
            {
                //fCentrality = fAOD->GetHeader()->GetCentralityP();
                fCentrality = ((AliAODHeader*)fAOD->GetHeader())->GetCentralityP();
            }
            else
            {
                fCentrality = fESD->GetCentrality();
            }
            
            if(fEstimator==1) centrality = fCentrality->GetCentralityPercentile("ZEMvsZDC");
            else centrality = fCentrality->GetCentralityPercentile("V0A");
            
            //printf("Centrality ZDC= %f VOA= %f  ZNA = %f ZNC = %f  ZPA = %f TRK = %f \n", fCentrality->GetCentralityPercentile("ZEMvsZDC") , fCentrality->GetCentralityPercentile("V0A"),fCentrality->GetCentralityPercentile("ZNA"), fCentrality->GetCentralityPercentile("ZNC"), fCentrality->GetCentralityPercentile("ZPA"), fCentrality->GetCentralityPercentile("TRK") );
            
            fCentralityHist->Fill(centrality);
            
            if(centrality<fCentralityMin || centrality>fCentralityMax) return;
            
            fCentralityHistPass->Fill(centrality);
        }
    }
    //______________________________________________________________________
    
    
    fNevent->Fill(17);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    //To use tender
    Int_t ClsNo0 = -999;
    //ClsNo0 = fAOD->GetNumberOfCaloClusters();
    if(fUseTender){
        //TClonesArray  *fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AODFilterTracks"));
        //NTracks = fTracks_tender->GetEntries();
        TClonesArray  *fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));
        ClsNo0 = fCaloClusters_tender->GetEntries();
        
        
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    //______________________________________________________________________
    
    if(fIsMC)
    {
        if(fIsAOD)
        {
            fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
            
            if(!fMCarray)
            {
                AliError("Array of MC particles not found");
                return;
            }
            
            fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
            
            if(!fMCheader)
            {
                AliError("Could not find MC Header in AOD");
                return;
            }
            
            for(Int_t iMC = 0; iMC < fMCarray->GetEntries(); iMC++)
            {
                fMCparticle = (AliAODMCParticle*) fMCarray->At(iMC);
                if(fMCparticle->GetMother()>0) fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                
                Int_t pdg = fMCparticle->GetPdgCode();
                
                //====================================================================
                //To test Non-Lin corrections
                Double_t energyMC_all=fMCparticle->E();
                
                
                
                if(fUseTender){
                    TClonesArray  *fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));
                    ClsNo0 = fCaloClusters_tender->GetEntries();
                    
                }
                
                for (Int_t icl=0; icl< ClsNo0; icl++ ){
                    
                    AliVCluster *clus0 = 0x0;
                    
                    if(fUseTender){
                        TClonesArray  *fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));
                        clus0 = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(icl));
                    }
                    
                    if(!fUseTender){
                        clus0 = (AliVCluster*) fAOD->GetCaloCluster(icl);
                    }
                    
                    
                    if(!clus0->IsEMCAL()) continue;
                    if(clus0->GetLabel() == iMC){
                        
                        if( TMath::Abs(pdg) == 22){
                            fEmc_Ereco_gamma_all->Fill( clus0->E(), energyMC_all);
                            fEmc_Ereco_gamma_ratio_all->Fill( clus0->E(), energyMC_all / clus0->E());
                        }
                        
                        if( TMath::Abs(pdg) == 11){
                            fEmc_Ereco_ele_all->Fill( clus0->E(), energyMC_all);
                            fEmc_Ereco_ele_ratio_all->Fill( clus0->E(), energyMC_all / clus0->E());
                            
                        }
                        
                    }
                }
                //====================================================================
                //End of test of Non-Lin corrections
                
                
                
                //====================================================================
                //trying take pions spectra 27/May/2014
                //IsPrimary only take events from pythia
                //IsPhysicalPrimariee: (all prompt particles, including strong decay products plus weak decay product from heavy-flavor).
                //eta cut same as MinJung
                
                if(fMCparticle->Eta()>=-0.8 && fMCparticle->Eta()<=0.8)
                {
                    if(fMCparticle->IsPrimary()){
                        
                        if(TMath::Abs(pdg)==111) fPtMCpi0->Fill(fMCparticle->Pt());
                        if(TMath::Abs(pdg)==221) fPtMCeta->Fill(fMCparticle->Pt());
                        //eta cut same as MinJung
                    }
                    
                    if(fMCparticle->IsPhysicalPrimary()){
                        if(TMath::Abs(pdg)==111) fPtMCpi02->Fill(fMCparticle->Pt());
                        if(TMath::Abs(pdg)==221) fPtMCeta2->Fill(fMCparticle->Pt());
                        
                    }
                    
                    if(TMath::Abs(pdg)==111) fPtMCpi03->Fill(fMCparticle->Pt());
                    if(TMath::Abs(pdg)==221) fPtMCeta3->Fill(fMCparticle->Pt());
                }
                //====================================================================
                
                double proX = fMCparticle->Xv();
                double proY = fMCparticle->Yv();
                double proR = sqrt(pow(proX,2)+pow(proY,2));
                
                
                if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax && fMCparticle->Charge()!=0)
                {
                    //to correct background
                    if (TMath::Abs(pdg) == 11 && fMCparticle->GetMother()>0){
                        Int_t mpdg = fMCparticleMother->GetPdgCode();
                        
                        if(TMath::Abs(mpdg) == 221 || TMath::Abs(mpdg) == 22 || TMath::Abs(mpdg) == 111){
                            
                            if(proR<7){
                                fPtMCparticleAlle_nonPrimary->Fill(fMCparticle->Pt()); //denominator for total efficiency for all electrons, and not primary
                                
                            }
                            
                        }
                    }
                    
                    if (TMath::Abs(pdg) == 11 && fMCparticle->IsPhysicalPrimary()){
                        fPtMCparticleAlle_Primary->Fill(fMCparticle->Pt()); //denominator for total efficiency for all electrons primary
                    }
                    
                    if( TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 2212 || TMath::Abs(pdg) == 321 || TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13 )
                    {
                        
                        fPtMCparticleAll_nonPrimary->Fill(fMCparticle->Pt()); //denominator for total efficiency for all particles, and not primary
                        if(fMCparticle->IsPhysicalPrimary())
                        {
                            fPtMCparticleAll->Fill(fMCparticle->Pt());
                            
                            Bool_t MotherFound = FindMother(iMC);
                            //Bool_t MotherFound = FindMother(track->GetLabel());
                            if(MotherFound)
                            {
                                if(fIsHFE1){
                                    //denominator for total efficiency and tracking
                                    //unfolding: denominator is pt_MC and numerator is pt_reco
                                    fPtMCparticleAllHfe1->Fill(fMCparticle->Pt());
                                    fEtaPhi_den->Fill(fMCparticle->Phi(),fMCparticle->Eta());
                                    
                                    
                                } //denominator for total efficiency and tracking
                                if(fIsHFE2){
                                    fPtMCparticleAllHfe2->Fill(fMCparticle->Pt());
                                }
                            }
                        }
                    }
                }//eta cut
                
                
                
            }//loop tracks
            
            
            
            //second loop over track, but only the primaries ones
            //only primary pions --> how to take the primaries ones in AOD?
            /*
             for(Int_t iMC = 0; iMC < fMCarray->GetNPrimary(); iMC++){
             fMCparticle = (AliAODMCParticle*) fMCarray->At(iMC);
             pdg = fMCparticle->GetPdgCode();
             
             if(TMath::Abs(pdg)==111) fPtMCpi0->Fill(fMCparticle->Pt());
             if(TMath::Abs(pdg)==221) fPtMCeta->Fill(fMCparticle->Pt());
             
             if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax)
             {
             
             if(TMath::Abs(pdg)==111) fPtMCpi02->Fill(fMCparticle->Pt());
             if(TMath::Abs(pdg)==221) fPtMCeta2->Fill(fMCparticle->Pt());
             
             }
             }
             */
            
            
        }//AOD
        else
        {
            fEventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
            if (!fEventHandler) {
                Printf("ERROR: Could not retrieve MC event handler");
                return;
            }
            
            
            
            fMCevent = fEventHandler->MCEvent();
            if (!fMCevent) {
                Printf("ERROR: Could not retrieve MC event");
                return;
            }
            
            fMCstack = fMCevent->Stack();
            
            //pion and eta spectrum
            //MinJung code
            
            //----------------------------------------------------------------------------------------------------
            AliVParticle *mctrack2 = NULL;
            AliMCParticle *mctrack0 = NULL;
            
            
            for(Int_t imc = 0; imc <fMCEvent->GetNumberOfPrimaries(); imc++){
                if(!(mctrack2 = fMCEvent->GetTrack(imc))) continue;
                TParticle* mcpart0 = fMCEvent->Stack()->Particle(imc);
                if(!mcpart0) continue;
                mctrack0 = dynamic_cast<AliMCParticle *>(mctrack2);
                if(!mctrack0) continue;
                
                if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8){
                    
                    if(TMath::Abs(mctrack0->PdgCode()) == 111) // pi0
                    {
                        fPtMCpi0->Fill(mctrack0->Pt());
                    }
                    
                    if(TMath::Abs(mctrack0->PdgCode()) == 221) // eta
                    {
                        fPtMCeta->Fill(mctrack0->Pt());
                    }
                    
                }
                
            }
            // end of MinJung
            //----------------------------------------------------------------------------------------------------
            
            
            for(Int_t iMC = 0; iMC < fMCstack->GetNtrack(); iMC++)
            {
                
                fMCtrack = fMCstack->Particle(iMC);
                if(fMCtrack->GetFirstMother()>0) fMCtrackMother = fMCstack->Particle(fMCtrack->GetFirstMother());
                TParticle *particle=fMCstack->Particle(iMC);
                
                Int_t pdg = fMCtrack->GetPdgCode();
                
                
                if(TMath::Abs(pdg)==111) fPtMCpi0->Fill(fMCtrack->Pt());
                if(TMath::Abs(pdg)==221) fPtMCeta->Fill(fMCtrack->Pt());
                
                
                if(fMCtrack->Eta()>=fEtaCutMin && fMCtrack->Eta()<=fEtaCutMax)
                {
                    
                    //to correct background
                    if (TMath::Abs(pdg) == 11 && fMCtrack->GetFirstMother()>0){
                        Int_t mpdg = fMCtrackMother->GetPdgCode();
                        if(TMath::Abs(mpdg) == 221 || TMath::Abs(mpdg) == 22 || TMath::Abs(mpdg) == 111){
                            Double_t proR=particle->R();
                            if(proR<7){
                                fPtMCparticleAlle_nonPrimary->Fill(fMCtrack->Pt()); //denominator for total efficiency for all electrons, and not primary
                            }
                        }
                    }
                    
                    if(TMath::Abs(pdg) == 11 && fMCstack->IsPhysicalPrimary(iMC)){
                        
                        fPtMCparticleAlle_Primary->Fill(fMCtrack->Pt());
                    }
                    
                    if( TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 2212 || TMath::Abs(pdg) == 321 || TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13 )
                    {
                        fPtMCparticleAll_nonPrimary->Fill(fMCtrack->Pt());//denominator for total efficiency for all particle, non Primary track
                        
                        if(fMCstack->IsPhysicalPrimary(iMC))
                        {
                            fPtMCparticleAll->Fill(fMCtrack->Pt());
                            
                            Bool_t MotherFound = FindMother(iMC);
                            //Bool_t MotherFound = FindMother(track->GetLabel());
                            if(MotherFound)
                            {
                                if(fIsHFE1){
                                    fPtMCparticleAllHfe1->Fill(fMCtrack->Pt());//denominator for total efficiency and tracking
                                    fEtaPhi_den->Fill(fMCtrack->Phi(),fMCtrack->Eta());
                                }
                                if(fIsHFE2){
                                    fPtMCparticleAllHfe2->Fill(fMCtrack->Pt());
                                }
                            }
                        }//Is Physical primary
                    }
                }//eta cut
            }//loop tracks
        }//ESD
    }//Is MC
    
    //______________________________________________________________________
    //threshold selection was here
    //______________________________________________________________________
    //all events selected
    
    fNevent->Fill(0);
    
    
    //______________________________________________________________________
    //events in the threshold
    
    if(firedTrigger.Contains(TriggerEG1))
    {
        if(fEMCEG1){
            fNevent->Fill(18);
            if(!firedTrigger.Contains(TriggerEG2)) fNevent->Fill(19);
            //if(firedTrigger.Contains(TriggerEG2)) return;
        }
    }
    
    
    //EG2
    if(firedTrigger.Contains(TriggerEG2))
    {
        if(fEMCEG2){
            fNevent->Fill(20);
            if(!firedTrigger.Contains(TriggerEG1)) fNevent->Fill(21);
            //if(firedTrigger.Contains(TriggerEG1)) return;
        }
    }
    
    
    
    //New cluster information
    //after trigger threshold selection
    Int_t ClsNo2 = -999;
    ClsNo2 = fVevent->GetNumberOfCaloClusters();
    fNCluster_pure->Fill(ClsNo2);
    
    
    /*
     if(ClsNo2<=0){
     fNevent->Fill(22); //events with no cluster
     return;
     }
     */
    
    //in order to include time cut
    //fEMCALCells = fAOD->GetEMCALCells();
    //Double_t tof = clus->GetTOF();
    //clus->GetNCells()
    //if ( clus->E() < minE ) continue ;
    
    
    /*
     if(fIsAOD){
     
     //AliAODHeader * aodh = fAOD->GetHeader();
     //Int_t bc= aodh->GetBunchCrossNumber();
     
     
     Int_t ClsNo_aod = -999;
     ClsNo_aod = fAOD->GetNumberOfCaloClusters();
     fNCluster_pure_aod->Fill(ClsNo_aod);
     //Bool_t exotic=kTRUE;
     
     for (Int_t i=0; i< ClsNo_aod; i++ ){
     
     //fClus = fVevent->GetCaloCluster(i);
     //to be compatible with Shingo
     AliVCluster *clust = 0x0;
     clust = (AliVCluster*) fAOD->GetCaloCluster(i);
     
     if(clust && clust->IsEMCAL())
     {
     //pure cluster information
     fECluster_pure->Fill(clust->E());
     
     fNcells_energy->Fill(clust->GetNCells(),clust->E());
     fNCluster_ECluster->Fill(ClsNo_aod,clust->E());
     
     if(clust->E()>1000) fNevent->Fill(23);
     
     //exotic
     
     exotic   = fEMCALRecoUtils->IsExoticCluster(clust, (AliVCaloCells*)fAOD->GetEMCALCells(), bc);
     if(exotic == kFALSE){
					fECluster_not_exotic->Fill(clust->E());
					fNcells_energy_not_exotic->Fill(clust->GetNCells(),clust->E());
     }
     
     
     //approximation to remove exotics
     if(clust->GetNCells()<5 && clust->E()>15.0){
					fECluster_exotic->Fill(clust->E());
     }
     //Marcel cut (another approximation to remove exotics)
     else if((clust->GetNCells())> ((clust->E())/3+1)){
					fECluster_not_exotic1->Fill(clust->E());
     }
     else{
					fECluster_not_exotic2->Fill(clust->E());
     }
     
     
     }
     
     
     }
     }
     */
    
    
    
    fNevent->Fill(24); //events with cluster
    
    
    fVtxZ_new4->Fill(fZvtx);
    
    
    Int_t ClsNo = -999;
    if(!fIsAOD) ClsNo = fESD->GetNumberOfCaloClusters();
    else ClsNo = fAOD->GetNumberOfCaloClusters();
    
    //______________________________________________________________________
    
    ///_____________________________________________________________________
    ///Track loop
    Int_t NTracks=0;
    AliVCluster *clust = 0x0;
    
    if(!fUseTender){
        if(fIsAOD){
            fNCluster_pure_aod->Fill(ClsNo);
            for (Int_t i=0; i< ClsNo; i++ ){
                clust = (AliVCluster*) fAOD->GetCaloCluster(i);
                
                if(clust && clust->IsEMCAL())
                {
                    fECluster_pure->Fill(clust->E());
                }
            }
        }
    }
    
    
    if(!fUseTender) NTracks=fVevent->GetNumberOfTracks();
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    //To use tender
    if(fUseTender){
        TClonesArray  *fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AODFilterTracks"));
        NTracks = fTracks_tender->GetEntries();
        TClonesArray  *fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));
        ClsNo = fCaloClusters_tender->GetEntries();
        
        
        //For cluster information from tender
        for (Int_t i=0; i< ClsNo; i++ ){
            TClonesArray  *fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));
            clust = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(i));
            if (!clust) {
                //printf("ERROR: Could not receive cluster matched calibrated from track %d\n", iTracks);
                continue;
            }
            if(clust && clust->IsEMCAL())
            {
                fECluster_pure->Fill(clust->E());
            }
            
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    //To calculate the efficiency of the Hadrons we look in all the Monte Carlo sample, which of them have a track associeted with them, ad if they have, we check if they meet the requeriments to be reconstructed
    
    if (fIsMC)
    {
        //Loop on the MC tracks to find the reco efficiency
        for(Int_t iTracksMC = 0; iTracksMC < fMCarray->GetEntries(); iTracksMC++)
        {
            //first we check if the MC particle is charged and is at the acceptance
            
            fMCparticle = (AliAODMCParticle*) fMCarray->At(iTracksMC);
            
            if (fMCparticle->Eta()<fEtaCutMin || fMCparticle->Eta()>fEtaCutMax) continue;
            
            if (fMCparticle->Charge() == 0) continue;
            if (!fMCparticle->IsPhysicalPrimary()) continue; //Physical primary
            
            //Save the pT of all Charged hadrons in the acceptance (This is the denominator of the efficiency)
            fpTPhiEffHadronsMC->Fill(fMCparticle->Pt(),fMCparticle->Phi());
            fpTEtaEffHadronsMC->Fill(fMCparticle->Pt(),fMCparticle->Eta());
            
            //Now loop over the tracks
            for(Int_t iTracks = 0; iTracks < NTracks; iTracks++)
            {
                AliVParticle* Vtrack = 0x0;
                
                if(!fUseTender) Vtrack  = fVevent->GetTrack(iTracks);
                
                AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
                
                AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);
                
                //Them we check if it has a matching AliVTrack
                
                if (track->GetLabel() == iTracksMC)
                {
                    //Now we check if the track has requisited bits, the eta range and the minimun number of TPC clusters
                    if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax && atrack->TestFilterMask(AliAODTrack::kTrkTPCOnly) && atrack->GetStatus()&AliESDtrack::kITSrefit && atrack->GetStatus()&AliESDtrack::kTPCrefit && atrack->GetTPCNcls() >= 80)
                    {
                        //Last check: DCA cut
                        //DCA cut for hadrons
                        if(fIsAOD && fUseDCACutforHadrons)
                        {
                            Double_t d0z0[2], cov[3];
                            AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
                            track->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 20., d0z0, cov);
                            Double_t DCAxy = d0z0[0];
                            Double_t DCAz = d0z0[1];
                            if(TMath::Abs(DCAxy) >= fDCAcutrHadron || TMath::Abs(DCAz)>=fDCAcutzHadron) continue;
                        }
                        
                        
                        //If positive, the particle is a reconstructed particle and should be saved in the pT spectrum (This is the numerator of the efficiency)
                        
                        fpTPhiEffHadronsReco->Fill(fMCparticle->Pt(),fMCparticle->Phi());
                        fpTEtaEffHadronsReco->Fill(fMCparticle->Pt(),fMCparticle->Eta());
                        
                        
                    }
                }
            }
            
        }
        
        //Loop on the real tracks to find the contamination from secoundary particles (and check consistency with correlation)
        for(Int_t iTracks = 0; iTracks < NTracks; iTracks++)
        {
            AliVParticle* Vtrack = 0x0;
            Vtrack  = fVevent->GetTrack(iTracks);
            AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
            AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);
            
            if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax && atrack->TestFilterMask(AliAODTrack::kTrkTPCOnly) && atrack->GetStatus()&AliESDtrack::kITSrefit && atrack->GetStatus()&AliESDtrack::kTPCrefit && atrack->GetTPCNcls() >= 80)
            {
                //Last check: DCA cut
                //DCA cut for hadrons
                if(fIsAOD && fUseDCACutforHadrons)
                {
                    Double_t d0z0[2], cov[3];
                    AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
                    track->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 20., d0z0, cov);
                    Double_t DCAxy = d0z0[0];
                    Double_t DCAz = d0z0[1];
                    if(TMath::Abs(DCAxy) >= fDCAcutrHadron || TMath::Abs(DCAz)>=fDCAcutzHadron) continue;
                }
                
                //This track is an associeated particle for the analysis, we do it's pT spectrum
                
                fpTPhiHadron->Fill(track->Pt(),track->Phi());
                fpTEtaHadron->Fill(track->Pt(),track->Eta());
                
                if (track->GetLabel()<0)
                {
                    fpTPhiSecPartNoLabel->Fill(track->Pt(),track->Phi());
                    fpTEtaSecPartNoLabel->Fill(track->Pt(),track->Eta());
                }
                else
                {
                    AliAODMCParticle *hadron_MonteCarlo = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
                    if (hadron_MonteCarlo->IsPhysicalPrimary())
                    {
                        fpTPhiSecPartPhysPri->Fill(track->Pt(),track->Phi());
                        fpTEtaSecPartPhysPri->Fill(track->Pt(),track->Eta());
                    }

                    else
                    {
                        fpTPhiSecPartNonPhysPri->Fill(track->Pt(),track->Phi());
                        fpTEtaSecPartNonPhysPri->Fill(track->Pt(),track->Eta());
                    }
                }
                
                
                
                
            }
            
            
            
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    ///_____________________________________________________________________
    ///Track loop
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    for(Int_t iTracks = 0; iTracks < NTracks; iTracks++)
    {
        //AliVParticle* Vtrack = fVevent->GetTrack(iTracks);
        AliVParticle* Vtrack = 0x0;
        if(!fUseTender) Vtrack  = fVevent->GetTrack(iTracks);
        
        
        //just a test of pdg information
        //Vtrack  = fVevent->GetTrack(iTracks);
        //AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
        //printf("\n\n Track label of track %d is %d  - from Vevent \n\n", iTracks, track->GetLabel());
        
        
        /*
         if(fIsMC  && track->GetLabel()>=0)
         {
         if(fIsAOD)
         {
         fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
         Int_t pdg = fMCparticle->GetPdgCode();
         printf("\n\n pdg code of track %d is %d  - from Vevent \n\n", iTracks, pdg);
         }
         }
         //just a test of pdg information -- remove until here
         */
        
        if(fUseTender){
            TClonesArray  *fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AODFilterTracks"));
            Vtrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(iTracks));
            TClonesArray  *fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));
        }
        
        if (!Vtrack)
        {
            printf("ERROR: Could not receive track %d\n", iTracks);
            continue;
        }
        
        
        AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
        AliESDtrack *etrack = dynamic_cast<AliESDtrack*>(Vtrack);
        AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);
        
        if(fIsMC  && atrack->GetLabel()>=0)
        {
            if(fIsAOD)
            {
                fMCparticle = (AliAODMCParticle*) fMCarray->At(atrack->GetLabel());
                
                
                //checking the same filter bits and cuts we use for hadrons and saving them to separete a possible pT shift. Only Reco hadrons are considered
                if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax && atrack->TestFilterMask(AliAODTrack::kTrkTPCOnly) && atrack->GetStatus()&AliESDtrack::kITSrefit && atrack->GetStatus()&AliESDtrack::kTPCrefit && atrack->GetTPCNcls() >= 80)
                {
                    fpTShifHadronsMC->Fill(fMCparticle->Pt());
                    fpTShiftHadronsReco->Fill(track->Pt());
                    
                }
                
                
                
                Int_t pdg = fMCparticle->GetPdgCode();
                if (TMath::Abs(pdg) == 11 && fMCparticle->Eta()>fEtaCutMin && fMCparticle->Eta()<fEtaCutMax)
                {
                    
                    if (fMCparticle->IsPrimary())
                    {
                        fMCEffPID_beforePID->Fill(fMCparticle->Pt());
                        Int_t Mother_pdg;
                        
                        if (fMCparticle->GetMother()>0)
                        {
                            fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                            Mother_pdg = fMCparticleMother->GetPdgCode();
                            
                        }
                        else
                            Mother_pdg = 0;
                        
                        if( TMath::Abs(Mother_pdg)==22 || TMath::Abs(Mother_pdg)==111 || TMath::Abs(Mother_pdg)==221)
                        {
                            fpTMCGeneratedBKG->Fill(fMCparticle->Pt());
                        }
                        
                        
                    }
                    
                }
                
            }
        }
        
        
        //printf("\n\n Track label of track %d is %d  - after tender \n\n", iTracks, track_tender->GetLabel());
        
        //remove
        /*
         if(fIsMC  && atrack->GetLabel()>=0)
         {
         if(fIsAOD)
         {
         fMCparticle = (AliAODMCParticle*) fMCarray->At(atrack->GetLabel());
         Int_t pdg = fMCparticle->GetPdgCode();
         printf("\n\n pdg code of track %d is %d  - from Tender \n\n", iTracks, pdg);
         }
         }
         //until here
         */
        
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        //To use tender
        /*
         if(fUseTender){
         TClonesArray  *fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AODFilterTracks"));
         track = static_cast<AliAODTrack*>(fTracks_tender->At(iTracks));
         if (!track) {
         printf("ERROR: Could not receive track calibrated %d\n", iTracks);
         continue;
         }
         }
         */
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        Double_t fTPCnSigma0 = -999;
        Double_t fTPCnSigma = -999;
        Double_t fTOFnSigma = -999;
        Double_t fTPCnSigma_pion = -999;
        Double_t fTPCnSigma_proton = -999;
        Double_t fTPCnSigma_kaon = -999;
        Double_t fTPCsignal = -999;
        Double_t fPt = -999;
        Double_t fP = -999;
        
        //December 9th 2013
        //Etacut test on the begging
        fEtad[0]->Fill(track->Eta());
        //if(track->Eta()<fEtaCutMin || track->Eta()>fEtaCutMax) continue;
        
        
        
        
        
        ///_____________________________________________________________________________
        ///Fill QA plots without track selection
        fPt = track->Pt();
        fP = TMath::Sqrt((track->Pt())*(track->Pt()) + (track->Pz())*(track->Pz()));
        
        fTPCsignal = track->GetTPCsignal();
        fTPCnSigma = fPidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        fTOFnSigma = fPidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
        fTPCnSigma_pion = fPidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
        fTPCnSigma_proton = fPidResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        fTPCnSigma_kaon = fPidResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        
        
        fTPC_momentum->Fill(fP,fTPCsignal);
        fTPC_eta->Fill(track->Eta(),fTPCsignal);
        
        
        
        if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax ){
            fTPCnsigma_pt_2D0->Fill(fPt,fTPCnSigma);
        }
        
        //printf("fCalibrateTPC = %d\n",fCalibrateTPC);
        
        if(fCalibrateTPC==kTRUE){
            //printf("Calibration of TPC is turned on \n");
            fTPCnSigma0 =(fTPCnSigma-fCalibrateTPC_mean)/fCalibrateTPC_sigma;
            //fTPCnSigma0 =(fTPCnSigma-3);
            //printf("TPCnsigma is %f, Mean is %f, Sigma is %f, TPCnsigmaCorr is %f\n", fTPCnSigma, fCalibrateTPC_mean,fCalibrateTPC_sigma, fTPCnSigma0);
            
        }
        
        if(fCalibrateTPC==kTRUE && track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax){
            fTPCnSigma = fTPCnSigma0;
            if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax ){
                fTPCnsigma_pt_2D1->Fill(fPt,fTPCnSigma);
                fTPCnsigma_pt_2D2->Fill(fPt,fTPCnSigma0);
            }
        }
        
        
        
        
        
        //================================================================================
        // Checks on SPD hits vs. SPD tracklets
        
        //Int_t ITSNCls = atrack->GetITSNcls();
        
        //printf("Track %d has ITS NCls =%d ", iTracks, ITSNCls);
        
        //if(atrack->HasPointOnITSLayer(0)) //printf("Track %d has point on first layer", iTracks);
        
        //	if(atrack->TestFilterMask(AliAODTrack::kTrkITSsa)){
        
        
        //} //printf("Track %d is in the filter bit kTrkITSsa", iTracks);
        
        
        //Bool_t IsSPDClusterVsTrackletBG(const AliESDEvent* esd, Bool_t fillHists = kFALSE);
        
        //================================================================================
        
        fTPC_p[0]->Fill(fPt,fTPCsignal);
        fTPCnsigma_p[0]->Fill(fP,fTPCnSigma);
        
        
        if(track->Eta()>-0.6 && track->Eta()<-0.3){
            fTPCnsigma_p_eta1[0]->Fill(fP,fTPCnSigma);
        }
        if(track->Eta()>-0.3 && track->Eta()<0){
            fTPCnsigma_p_eta2[0]->Fill(fP,fTPCnSigma);
        }
        if(track->Eta()>0 && track->Eta()<0.3){
            fTPCnsigma_p_eta3[0]->Fill(fP,fTPCnSigma);
        }
        if(track->Eta()>0.3 && track->Eta()<0.6){
            fTPCnsigma_p_eta4[0]->Fill(fP,fTPCnSigma);
        }
        
        
        
        Float_t TPCNcls = track->GetTPCNcls();
        //TPC Ncls for pid
        Float_t TPCNcls_pid = track->GetTPCsignalN();
        
        
        
        Float_t pos[3]={0,0,0};
        
        Double_t fEMCflag = kFALSE;
        if(track->GetEMCALcluster()>0)
        {
            //fClus = fVevent->GetCaloCluster(track->GetEMCALcluster());
            
            if(!fUseTender) fClus = fVevent->GetCaloCluster(track->GetEMCALcluster());
            
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            //to use tender
            if(fUseTender){
                int EMCalIndex = -1;
                EMCalIndex = track->GetEMCALcluster();
                if(EMCalIndex>0){
                    TClonesArray  *fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));
                    fClus = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(EMCalIndex));
                    if (!fClus) {
                        //printf("ERROR: Could not receive cluster matched calibrated from track %d\n", iTracks);
                        continue;
                    }
                }
            }
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            // printf("\n\n !!! test !!! cluster label of track %d is %d  - from tender \n\n", iTracks, fClus->GetLabel());
            
            
            if(fClus->IsEMCAL())
            {
                
                //only for charged tracks
                fECluster[0]->Fill(fClus->E());
                
                
                if(TMath::Abs(fClus->GetTrackDx())<=fdPhiCut && TMath::Abs(fClus->GetTrackDz())<=fdEtaCut)
                {
                    fEMCflag = kTRUE;
                    fEoverP_pt[0]->Fill(fPt,(fClus->E() / fP));
                    
                    
                    Float_t Energy	= fClus->E();
                    Float_t EoverP	= Energy/track->P();
                    //Float_t M02	= fClus->GetM02();
                    //Float_t M20	= fClus->GetM20();
                    
                    /////////////// for Eta Phi distribution
                    fClus->GetPosition(pos);
                    TVector3 vpos(pos[0],pos[1],pos[2]);
                    Double_t cphi = vpos.Phi();
                    Double_t ceta = vpos.Eta();
                    fEtaPhi[0]->Fill(cphi,ceta);
                    
                    
                    fTPCNcls_EoverP[0]->Fill(TPCNcls, EoverP);
                    
                    
                    if(fIsMC && fIsAOD && track->GetLabel()>=0){
                        fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
                        Int_t pdg = fMCparticle->GetPdgCode();
                        
                        
                        //to study the response of detector -> NonLinearity function
                        Double_t energyMC=fMCparticle->E();
                        
                        
                        if( TMath::Abs(pdg) == 22){
                            fEmc_Ereco_gamma0->Fill( fClus->E(), energyMC);
                            fEmc_Ereco_gamma_ratio0->Fill( fClus->E(), energyMC / fClus->E());
                        }
                        
                        if( TMath::Abs(pdg) == 11){
                            fEmc_Ereco_ele0->Fill( fClus->E(), energyMC);
                            fEmc_Ereco_ele_ratio0->Fill( fClus->E(), energyMC / fClus->E());
                            
                        }
                        
                        if(fClus->E() / fP > 1.2)fEmc_Ereco_ratio_large_EoverP0->Fill( fClus->E(), energyMC / fClus->E());
                        if(fClus->E() / fP < 1.2)fEmc_Ereco_ratio_small_EoverP0->Fill( fClus->E(), energyMC / fClus->E());
                        
                    }
                    
                    
                    
                }
            }
        }
        
        //______________________________________________________________
        // Vertex
        
        fVtxZ[0]->Fill(fZvtx);
        if(iTracks == 0)fNTracks[0]->Fill(NTracks);
        fNTracks_pt[0]->Fill(NTracks, fPt);
        fNTracks_eta[0]->Fill(NTracks, track->Eta());
        fNTracks_phi[0]->Fill(NTracks, track->Phi());
        
        
        fNClusters[0]->Fill(ClsNo);
        fTPCNcls_pid[0]->Fill(TPCNcls, TPCNcls_pid);
        //______________________________________________________________
        
        ///Fill QA plots without track selection
        ///_________________________________________________________________________
        //__________________________________________________________________________
        //Track Selection Cuts
        
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
            
            //if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;
            
            
            
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
        
        //RecPrim
        
        //it was not working on aod... testing again
        //July 29th, 2014: aparently working again
        
        ///if(!fIsAOD)
        //{
        
        if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;
        //}
        
        //HFEcuts: ITS layers cuts
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;
        
        //HFE cuts: TPC PID cleanup
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;
        
        //DCA cut done by hand --  has to be in absolute values!!!! Fixed
        //not using because line if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;   already apply DCA cut
        /*
         if(fIsAOD){
         Double_t d0z0[2], cov[3];
         AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
         track->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 20., d0z0, cov);
         Double_t DCAxy = d0z0[0];
         Double_t DCAz = d0z0[1];
         
         if(TMath::Abs(DCAxy) >= fDCAcutr || TMath::Abs(DCAz)>=fDCAcutz) continue;
         }
         */
        
        fEtad[1]->Fill(track->Eta());
        
        
        //______________________________________________________________________________________
        /*
         if(fIsAOD){
         //AOD test -- Francesco suggestion
         //aod test -- Francesco suggestion
         AliAODTrack *aod_track=dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));
         
         Int_t type=aod_track->GetType();
         if(type==AliAODTrack::kPrimary) fPtPrim->Fill(aod_track->Pt());
         if(type==AliAODTrack::kSecondary) fPtSec->Fill(aod_track->Pt());
         
         //Int_t type2=track->GetType();
         //if(type==AliAODTrack::kPrimary) fPtPrim2->Fill(track->Pt());
         //if(type==AliAODTrack::kSecondary) fPtSec2->Fill(track->Pt());
         }
         */
        
        
        ///_____________________________________________________________
        ///QA plots after track selection
        if(fIsMC)
        {
            if(track->GetLabel()>=0) fPtMCWithLabel->Fill(fPt);
            else fPtMCWithoutLabel->Fill(fPt);
        }
        
        if(fIsMC  && track->GetLabel()>=0)
        {
            if(fIsAOD)
            {
                fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
                
                
                if(fMCparticle->IsPhysicalPrimary())
                {
                    fPtIsPhysicaPrimary->Fill(fPt);
                }
                
                Int_t pdg = fMCparticle->GetPdgCode();
                
                //printf("\n\n pdg code of track %d is %d \n\n", iTracks, pdg);
                
                if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax && fMCparticle->Charge()!=0)
                {
                    
                    
                    if( TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 2212 || TMath::Abs(pdg) == 321 || TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13 )
                    {
                        fPtMCparticleReco_nonPrimary->Fill(fMCparticle->Pt()); //not Primary track
                        
                        if(fMCparticle->IsPhysicalPrimary())
                        {
                            fPtMCparticleReco->Fill(fMCparticle->Pt());
                            
                            Bool_t MotherFound = FindMother(track->GetLabel());
                            if(MotherFound)
                            {
                                if(fIsHFE1)
                                {
                                    fPtMCparticleRecoHfe1->Fill(fMCparticle->Pt());//numerator tracking
                                    //unfolding
                                    fpt_reco_pt_MC_den->Fill(track->Pt(),fMCparticle->Pt());
                                    
                                }
                                if(fIsHFE2){
                                    fPtMCparticleRecoHfe2->Fill(fMCparticle->Pt());
                                }
                            }
                        }
                    }
                }
                
                
                
            }//close AOD
            //ESD
            else
            {
                
                
                if(fMCtrack->Eta()>=fEtaCutMin && fMCtrack->Eta()<=fEtaCutMax)
                {
                    
                    fMCtrack = fMCstack->Particle(track->GetLabel());
                    Int_t pdg = fMCtrack->GetPdgCode();
                    
                    if( TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 2212 || TMath::Abs(pdg) == 321 || TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13 )
                    {
                        fPtMCparticleReco_nonPrimary->Fill(fMCtrack->Pt());//not Primary track
                    }
                    
                    
                    if(fMCstack->IsPhysicalPrimary(track->GetLabel()))
                    {
                        fPtIsPhysicaPrimary->Fill(fPt);
                        
                        
                        if( TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 2212 || TMath::Abs(pdg) == 321 || TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13 )
                        {
                            fPtMCparticleReco->Fill(fMCtrack->Pt());
                            
                            Bool_t MotherFound = FindMother(track->GetLabel());
                            if(MotherFound)
                            {
                                if(fIsHFE1) fPtMCparticleRecoHfe1->Fill(fMCtrack->Pt());//numerator tracking
                                if(fIsHFE2) fPtMCparticleRecoHfe2->Fill(fMCtrack->Pt());
                            }
                        }
                    }
                }
            }//close ESD
        }//close IsMC
        
        fTPC_momentum1->Fill(fP,fTPCsignal);
        fTPC_eta1->Fill(track->Eta(),fTPCsignal);
        
        //Calibrated TPCsignal after tracks selection
        if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax ){
            fTPCnsigma_pt_2D3->Fill(fPt,fTPCnSigma);
        }
        
        
        if(track->Eta()>-0.6 && track->Eta()<-0.3){
            fTPCnsigma_p_eta1[1]->Fill(fP,fTPCnSigma);
        }
        if(track->Eta()>-0.3 && track->Eta()<0){
            fTPCnsigma_p_eta2[1]->Fill(fP,fTPCnSigma);
        }
        if(track->Eta()>0 && track->Eta()<0.3){
            fTPCnsigma_p_eta3[1]->Fill(fP,fTPCnSigma);
        }
        if(track->Eta()>0.3 && track->Eta()<0.6){
            fTPCnsigma_p_eta4[1]->Fill(fP,fTPCnSigma);
        }
        
        
        fTPC_p[1]->Fill(fPt,fTPCsignal);
        fTPCnsigma_p[1]->Fill(fP,fTPCnSigma);
        Double_t fPtBin[7];
        if(fUseAlternativeBinnig)
        {
            fPtBin[0] = 0.5;
            fPtBin[1] = 1.0;
            fPtBin[2] = 1.5;
            fPtBin[3] = 2.0;
            fPtBin[4] = 3.0;
            fPtBin[5] = 4.0;
            fPtBin[6] = 6.0;
        }
        else
        {
            fPtBin[0] = 1.0;
            fPtBin[1] = 2.0;
            fPtBin[2] = 4.0;
            fPtBin[3] = 6.0;
            fPtBin[4] = 8.0;
            fPtBin[5] = 10.0;
            fPtBin[6] = 15.0;
        }
        
        Double_t fPtBin_trigger[11] = {1,2,4,6,8,10,12,14,16,18,20};
        
        TPCNcls = track->GetTPCNcls();
        Float_t pos2[3]={0,0,0};
        
        if(track->GetEMCALcluster()>0)
        {
            if(!fUseTender) fClus = fVevent->GetCaloCluster(track->GetEMCALcluster());
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            //to use tender
            if(fUseTender){
                int EMCalIndex = -1;
                EMCalIndex = track->GetEMCALcluster();
                if(EMCalIndex>0){
                    TClonesArray  *fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));
                    fClus = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(EMCalIndex));
                    if (!fClus) {
                        //printf("ERROR: Could not receive cluster matched calibrated from track %d\n", iTracks);
                        continue;
                    }
                }
            }
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            if(fClus->IsEMCAL())
            {
                if(TMath::Abs(fClus->GetTrackDx())<=fdPhiCut && TMath::Abs(fClus->GetTrackDz())<=fdEtaCut)
                {
                    fEoverP_pt[1]->Fill(fPt,(fClus->E() / fP));
                    
                    Float_t Energy	= fClus->E();
                    Float_t EoverP	= Energy/track->P();
                    Float_t M02	= fClus->GetM02();
                    Float_t M20	= fClus->GetM20();
                    Float_t ncells	= fClus->GetNCells();
                    
                    /////////////// for Eta Phi distribution
                    fClus->GetPosition(pos2);
                    TVector3 vpos(pos2[0],pos2[1],pos2[2]);
                    Double_t cphi = vpos.Phi();
                    Double_t ceta = vpos.Eta();
                    fEtaPhi[1]->Fill(cphi,ceta);
                    
                    fECluster[1]->Fill(Energy);
                    fTPCNcls_EoverP[1]->Fill(TPCNcls, EoverP);
                    
                    
                    //for EMCal trigger performance
                    if(EoverP > 0.9){
                        ftpc_p_EoverPcut->Fill(track->P(), fTPCsignal);
                        fnsigma_p_EoverPcut->Fill(track->P(), fTPCnSigma);
                        
                    }
                    
                    
                    //for hadron contamination calculations
                    
                    
                    // EtaCut -> dados
                    if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax ){
                        
                        //-------------------------------------------------------------------
                        //true hadrons and true electrons E/p shape before cuts
                        if(fIsMC && fIsAOD && track->GetLabel()>=0){
                            fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
                            Int_t pdg = fMCparticle->GetPdgCode();
                            
                            
                            //true hadrons
                            if( TMath::Abs(pdg) != 11){
                                fEoverP_pt_true_hadrons0->Fill(fPt,(fClus->E() / fP));
                            }
                            //true electrons
                            if( TMath::Abs(pdg) == 11){
                                fEoverP_pt_true_electrons0->Fill(fPt,(fClus->E() / fP));
                            }
                            
                            
                            //to study the response of detector -> NonLinearity function
                            Double_t energyMC=fMCparticle->E();
                            //printf("\n MC energy is =%f\n", energyMC);
                            fEmc->Fill(energyMC);
                            
                            if( TMath::Abs(pdg) == 22){
                                fEmc_Ereco_gamma->Fill( fClus->E(), energyMC);
                                fEmc_Ereco_gamma_ratio->Fill( fClus->E(), energyMC / fClus->E());
                            }
                            
                            if( TMath::Abs(pdg) == 11){
                                fEmc_Ereco_ele->Fill( fClus->E(), energyMC);
                                fEmc_Ereco_ele_ratio->Fill( fClus->E(), energyMC / fClus->E());
                                
                            }
                            
                            if(fClus->E() / fP > 1.2)fEmc_Ereco_ratio_large_EoverP->Fill( fClus->E(), energyMC / fClus->E());
                            if(fClus->E() / fP < 1.2)fEmc_Ereco_ratio_small_EoverP->Fill( fClus->E(), energyMC / fClus->E());
                            
                        }
                        
                        
                        
                        //-------------------------------------------------------------------
                        
                        
                        //main
                        if(TMath::Abs(fTPCnSigma_pion)<3 || TMath::Abs(fTPCnSigma_proton)<3 || TMath::Abs(fTPCnSigma_kaon)<3 ){
                            
                            if(fTPCnSigma<-3.5){
                                fEoverP_pt_hadrons->Fill(fPt,EoverP);
                                if(fUseEMCal) fShowerShape_ha->Fill(M02,M20);
                                
                                
                                //after cut -> Are they really hadrons? E/p shape
                                if(fIsMC && fIsAOD && track->GetLabel()>=0){
                                    fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
                                    Int_t pdg = fMCparticle->GetPdgCode();
                                    
                                    if( TMath::Abs(pdg) != 11){
                                        fEoverP_pt_true_hadrons->Fill(fPt,(fClus->E() / fP));
                                    }
                                }
                                
                                
                            }
                        }
                        
                        //for systematic studies of hadron contamination
                        if(fTPCnSigma < -3){
                            fEoverP_pt_pions->Fill(fPt, EoverP);
                            
                        }
                        
                        if(fTPCnSigma < -3.5){
                            fEoverP_pt_pions2->Fill(fPt, EoverP);
                            
                            if(fClus->E()>8.0)fEoverP_pt_pions2_highE1->Fill(fPt, EoverP);
                            if(fClus->E()>12.0)fEoverP_pt_pions2_highE0->Fill(fPt, EoverP);
                            
                        }
                        if(fTPCnSigma < -4.0){
                            fEoverP_pt_pions3->Fill(fPt, EoverP);
                            
                            
                        }
                        
                        
                    }
                    
                    
                    
                    
                    for(Int_t i = 0; i < 6; i++)
                    {
                        if(fPt>=fPtBin[i] && fPt<fPtBin[i+1])
                        {
                            
                            if(fTPCnSigma < -5 && fTPCnSigma > -10){
                                fNcells_hadrons[i]->Fill(ncells);
                            }
                            //hadrons selection using E/p
                            if((fClus->E() / fP) >= 0.2 && (fClus->E() / fP) <= 0.7){
                                fTPCnsigma_eta_hadrons[i]->Fill(fTPCnSigma, ceta);
                            }
                            //electrons selection using E/p
                            if((fClus->E() / fP) >= fEoverPCutMin && (fClus->E() / fP) <= fEoverPCutMax) {
                                fTPCnsigma_eta_electrons[i]->Fill(fTPCnSigma, ceta);
                            }
                        }
                    }
                    
                    if((fClus->E() / fP) >= fEoverPCutMin && (fClus->E() / fP) <= fEoverPCutMax)
                    {
                        fTPCnsigma_eta->Fill(track->Eta(),fTPCnSigma);
                        fTPCnsigma_phi->Fill(track->Phi(),fTPCnSigma);
                        
                        if(fUseEMCal)
                        {
                            if(fTPCnSigma < 3.5 && fCorrelationFlag)
                            {
                                fPtTrigger_Inc->Fill(fPt);
                                DiHadronCorrelation(track, iTracks);
                            }
                        }
                    }
                    
                }
            }
        }
        
        //______________________________________________________________
        // Vertex
        
        fVtxZ[1]->Fill(fZvtx);
        if(iTracks == 0)fNTracks[1]->Fill(NTracks);
        fNTracks_pt[1]->Fill(NTracks, fPt);
        fNTracks_eta[1]->Fill(NTracks, track->Eta());
        fNTracks_phi[1]->Fill(NTracks, track->Phi());
        fNClusters[1]->Fill(ClsNo);
        fTPCNcls_pid[1]->Fill(TPCNcls, TPCNcls_pid);
        
        //______________________________________________________________
        
        ///______________________________________________________________________
        ///Histograms for PID Studies
        //Double_t fPtBin[6] = {2,4,6,8,10,15};
        
        for(Int_t i = 0; i < 6; i++)
        {
            if(fPt>=fPtBin[i] && fPt<fPtBin[i+1])
            {
                if(fEMCflag) fEoverP_tpc[i]->Fill(fTPCnSigma,(fClus->E() / fP));
                fTPC_pt[i]->Fill(fTPCsignal);
                fTPCnsigma_pt[i]->Fill(fTPCnSigma);
                
            }
        }
        if(fUseTrigger)
        {
            if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax ){
                
                //new pt bins for trigger data
                
                for(Int_t i = 0; i < 10; i++)
                {
                    if(fP>=fPtBin_trigger[i] && fP<fPtBin_trigger[i+1])
                    {
                        if(fEMCflag)fEoverP_tpc_p_trigger[i]->Fill(fTPCnSigma,(fClus->E() / fP));
                        
                    }
                    
                    if(fPt>=fPtBin_trigger[i] && fPt<fPtBin_trigger[i+1])
                    {
                        if(fEMCflag)fEoverP_tpc_pt_trigger[i]->Fill(fTPCnSigma,(fClus->E() / fP));
                        
                    }
                }
            }
            
            
            //new way to calculate TPCnsigma distribution: TPCnsigma in function of p, with/without E/p cut
            fTPCnsigma_p_TPC->Fill(fP, fTPCnSigma);
            
            if(fEMCflag){
                
                /*
                 if(fEMCEG1)
                 {
                 if(fClus->E()>=12.00){
                 fTPCnsigma_p_TPC_on_EMCal_acc->Fill(fP, fTPCnSigma);
                 }
                 }
                 if(fEMCEG2)
                 {
                 if(fClus->E()>=8.00){
                 fTPCnsigma_p_TPC_on_EMCal_acc->Fill(fP, fTPCnSigma);
                 }
                 }
                 */
                
                //else
                //{
                
                fTPCnsigma_p_TPC_on_EMCal_acc->Fill(fP, fTPCnSigma);
                //}
                
                
                if((fClus->E() / fP) >= fEoverPCutMin && (fClus->E() / fP) <= fEoverPCutMax){
                    /*
                     if(fEMCEG1)
                     {
                     if(fClus->E()>=12.00){
                     fTPCnsigma_p_TPC_EoverP_cut->Fill(fP, fTPCnSigma);
                     
                     }
                     }
                     if(fEMCEG2)
                     {
                     if(fClus->E()>=8.00){
                     fTPCnsigma_p_TPC_EoverP_cut->Fill(fP, fTPCnSigma);
                     
                     }
                     }
                     */
                    
                    //else
                    //{
                    fTPCnsigma_p_TPC_EoverP_cut->Fill(fP, fTPCnSigma);
                    
                    //}
                    
                    
                    
                    
                }
                
            }//close EMCflag
            
        }//close eta cut
        
        ///QA plots after track selection
        ///_____________________________________________________________
        
        //_______________________________________________________
        //Correlation Analysis - DiHadron
        if(!fUseEMCal)
        {
            if(fTPCnSigma < 3.5 && fCorrelationFlag)
            {
                fPtTrigger_Inc->Fill(fPt);
                DiHadronCorrelation(track, iTracks);
            }
        }
        //_______________________________________________________
        
        
        ///////////////////////////////////////////////////////////////////
        ///TPC - efficiency calculation //
        
        /// changing start here
        if(fIsMC && fIsAOD && track->GetLabel()>=0)
        {
            fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
            Int_t pdg = fMCparticle->GetPdgCode();
            
            //
            if(fMCparticle->IsPhysicalPrimary()){
                
                
                if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax ){
                    
                    Bool_t MotherFound = FindMother(track->GetLabel());
                    if(MotherFound){
                        fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                        if( TMath::Abs(pdg) == 11 && fMCparticleMother->GetPdgCode()!=22 ){
                            if(fIsHFE1) fPtMC_TPC_All->Fill(fMCparticle->Pt());
                        }
                    }
                }
            }
        }///until here
        
        else if(fIsMC && track->GetLabel()>=0)//ESD
        {
            
            if(fMCstack->IsPhysicalPrimary(track->GetLabel())){
                fMCtrack = fMCstack->Particle(track->GetLabel());
                
                Int_t pdg = fMCtrack->GetPdgCode();
                if(fMCtrack->Eta()>=fEtaCutMin && fMCtrack->Eta()<=fEtaCutMax ){
                    
                    if(fMCtrack->GetFirstMother()>0){
                        fMCtrackMother = fMCstack->Particle(fMCtrack->GetFirstMother());
                        Bool_t MotherFound = FindMother(track->GetLabel());
                        if(MotherFound){
                            if( TMath::Abs(pdg) == 11 && fMCtrackMother->GetPdgCode()!=22 ){
                                if(fIsHFE1) fPtMC_TPC_All->Fill(fMCtrack->Pt());
                            }
                        }
                    }
                }
            }
        }
        
        
        if(fPt>1 && fPt<2) fTOF01->Fill(fTOFnSigma,fTPCnSigma);
        if(fPt>2 && fPt<4) fTOF02->Fill(fTOFnSigma,fTPCnSigma);
        if(fPt>4 && fPt<6) fTOF03->Fill(fTOFnSigma,fTPCnSigma);
        
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
        ///________________________________________________________________________
        
        //Should be done only for 13d, which has shifted TPCnsigma
        if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax){
            fTPCnsigma_pt_2D4->Fill(fPt,fTPCnSigma);
        }
        
        if(fCalibrateTPC==kTRUE){
            
            if(fTPCnSigma < fTPCcal_CutMin || fTPCnSigma > fTPCcal_CutMax) continue;
            
        }
        
        if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax){
            fTPCnsigma_pt_2D5->Fill(fPt,fTPCnSigma);
        }
        
        if(track->Eta()>-0.6 && track->Eta()<-0.3){
            fTPCnsigma_p_eta1[2]->Fill(fP,fTPCnSigma);
        }
        if(track->Eta()>-0.3 && track->Eta()<0){
            fTPCnsigma_p_eta2[2]->Fill(fP,fTPCnSigma);
        }
        if(track->Eta()>0 && track->Eta()<0.3){
            fTPCnsigma_p_eta3[2]->Fill(fP,fTPCnSigma);
        }
        if(track->Eta()>0.3 && track->Eta()<0.6){
            fTPCnsigma_p_eta4[2]->Fill(fP,fTPCnSigma);
        }
        
        
        ////////////////////////////////////////////////////////////////////
        ///TPC efficiency calculations
        
        
        if(fIsMC && fIsAOD && track->GetLabel()>=0)
        {
            fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
            Int_t pdg = fMCparticle->GetPdgCode();
            
            //
            if(fMCparticle->IsPhysicalPrimary()){
                
                
                if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax ){
                    
                    Bool_t MotherFound = FindMother(track->GetLabel());
                    if(MotherFound){
                        fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                        if( TMath::Abs(pdg) == 11 && fMCparticleMother->GetPdgCode()!=22 ){
                            if(fIsHFE1) fPtMC_TPC_Selected->Fill(fMCparticle->Pt());
                        }
                    }
                }
            }
        }///until here
        
        else if(fIsMC && track->GetLabel()>=0)//ESD
        {
            
            if(fMCstack->IsPhysicalPrimary(track->GetLabel())){
                fMCtrack = fMCstack->Particle(track->GetLabel());
                
                Int_t pdg = fMCtrack->GetPdgCode();
                if(fMCtrack->Eta()>=fEtaCutMin && fMCtrack->Eta()<=fEtaCutMax ){
                    
                    if(fMCtrack->GetFirstMother()>0){
                        fMCtrackMother = fMCstack->Particle(fMCtrack->GetFirstMother());
                        Bool_t MotherFound = FindMother(track->GetLabel());
                        if(MotherFound){
                            if( TMath::Abs(pdg) == 11 && fMCtrackMother->GetPdgCode()!=22 ){
                                if(fIsHFE1) fPtMC_TPC_Selected->Fill(fMCtrack->Pt());
                            }
                        }
                    }
                }
            }
        }
        
        //Eta Cut for TPC only
        if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax ){
            fTPCnsigma_pt_2D->Fill(fPt,fTPCnSigma);
        }
        
        //Background for TPC only
        if(fFillBackground){
            
            //efficiency without SS cut for TPC only
            if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax){
                //with weight for TPConly
                Background(track, iTracks, Vtrack, kTRUE, kTRUE, kFALSE); //IsTPConly=kTRUE, IsWeight=kTRUE
                //pt bins
                Background(track, iTracks, Vtrack, kTRUE, kTRUE, kTRUE); //IsTPConly=kTRUE, IsWeight=kTRUE, MassPtBins=kTRUE
                
                //default
                Background(track, iTracks, Vtrack, kTRUE, kFALSE, kFALSE); //IsTPConly=kTRUE
                
                
            } //Eta cut to be consistent with other efficiency
        }
        
        
        fTPCnsigma_p[2]->Fill(fP,fTPCnSigma);
        fTPC_p[2]->Fill(fP,fTPCsignal);
        TPCNcls = track->GetTPCNcls();
        Float_t pos3[3]={0,0,0};
        
        
        //here denominator for track-matching efficiency
        if(fIsMC && fIsAOD && track->GetLabel()>=0)
        {
            fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
            Int_t pdg = fMCparticle->GetPdgCode();
            
            //
            if(fMCparticle->IsPhysicalPrimary()){
                
                
                if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax ){
                    
                    Bool_t MotherFound = FindMother(track->GetLabel());
                    if(MotherFound){
                        fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                        if( TMath::Abs(pdg) == 11 && fMCparticleMother->GetPdgCode()!=22 ){
                            if(fIsHFE1) fPt_track_match_den->Fill(fMCparticle->Pt());
                        }
                    }
                }
            }
        }///until here
        
        
        if(track->GetEMCALcluster()>0)
        {
            fClus = fVevent->GetCaloCluster(track->GetEMCALcluster());
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            //to use tender
            if(fUseTender){
                int EMCalIndex = -1;
                EMCalIndex = track->GetEMCALcluster();
                if(track->GetEMCALcluster()>0){
                    TClonesArray  *fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));
                    fClus = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(EMCalIndex));
                    if (!fClus) {
                        //printf("ERROR: Could not receive cluster matched calibrated from track %d\n", iTracks);
                        continue;
                    }
                }
            }
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            if(fClus->IsEMCAL())
            {
                
                //________________________________________________________________________
                
                
                //Cluster timing distribution -- for ESD
                if(fUseEMCal && !fIsAOD){
                    
                    AliESDCaloCells &cells_esd=*(fESD->GetEMCALCells());
                    TRefArray* caloClusters_esd = new TRefArray();
                    fESD->GetEMCALClusters(caloClusters_esd);
                    Int_t nclus_esd = caloClusters_esd->GetEntries();
                    
                    
                    for (Int_t icl = 0; icl < nclus_esd; icl++) {
                        
                        AliESDCaloCluster* clus_esd = (AliESDCaloCluster*)caloClusters_esd->At(icl);
                        
                        if(clus_esd->IsEMCAL()){
                            Float_t ncells_esd	= fClus->GetNCells();
                            UShort_t * index_esd = clus_esd->GetCellsAbsId() ;
                            UShort_t * index2_esd = fClus->GetCellsAbsId() ;
                            
                            
                            for(Int_t i = 0; i < ncells_esd ; i++){
                                
                                Int_t absId_esd =   index_esd[i];
                                fTime->Fill(fPt,cells_esd.GetCellTime(absId_esd));
                                
                                Int_t absId2_esd =   index2_esd[i];
                                fTime2->Fill(fPt,cells_esd.GetCellTime(absId2_esd));
                            }
                            
                        }
                    }
                    
                }
                /* not working!
                 //Cluster timing distribution -- for AOD
                 if(fUseEMCal && fIsAOD){
                 
                 AliAODCaloCells &cells_aod=*(fAOD->GetEMCALCells());
                 
                 TRefArray* caloClusters_aod = new TRefArray();
                 fAOD->GetEMCALClusters(caloClusters_aod);
                 
                 Int_t nclus_aod = caloClusters_aod->GetEntries();
                 
                 for (Int_t icl = 0; icl < nclus_aod; icl++) {
                 
                 AliAODCaloCluster* clus_aod = (AliAODCaloCluster*)caloClusters_aod->At(icl);
                 
                 
                 if(clus_aod->IsEMCAL()){
                 Float_t ncells_aod	= fClus->GetNCells();
                 UShort_t * index_aod = clus_aod->GetCellsAbsId() ;
                 UShort_t * index2_aod = fClus->GetCellsAbsId() ;
                 
                 
                 for(Int_t i = 0; i < ncells_aod ; i++){
                 
                 Int_t absId_aod =   index_aod[i];
                 fTime->Fill(fPt,cells_aod.GetCellTime(absId_aod));
                 
                 Int_t absId2_aod =   index2_aod[i];
                 fTime2->Fill(fPt,cells_aod.GetCellTime(absId2_aod));
                 }
                 
                 }
                 }
                 
                 }
                 */
                
                
                if(fUseEMCal){
                    double emctof = fClus->GetTOF();
                    ftimingEle->Fill(fPt,1e9*emctof);
                }
                //________________________________________________________________________
                
                
                
                
                // Residuals
                Double_t Dx = fClus->GetTrackDx();
                Double_t Dz = fClus->GetTrackDz();
                Double_t R=TMath::Sqrt(Dx*Dx+Dz*Dz);
                
                for(Int_t i = 0; i < 6; i++)
                {
                    if(fPt>=fPtBin[i] && fPt<fPtBin[i+1])
                    {
                        
                        fEta[i]->Fill(Dz);
                        fPhi[i]->Fill(Dx);
                        fR[i]->Fill(R);
                    }
                }
                
                if(TMath::Abs(fClus->GetTrackDx())<=fdPhiCut && TMath::Abs(fClus->GetTrackDz())<=fdEtaCut)
                {
                    Float_t Energy	= fClus->E();
                    Float_t EoverP	= Energy/track->P();
                    Float_t M02	= fClus->GetM02();
                    Float_t M20	= fClus->GetM20();
                    Float_t ncells	= fClus->GetNCells();
                    
                    //here numerator for track-matching efficiency
                    if(fIsMC && fIsAOD && track->GetLabel()>=0)
                    {
                        fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
                        Int_t pdg = fMCparticle->GetPdgCode();
                        
                        //
                        if(fMCparticle->IsPhysicalPrimary()){
                            
                            
                            if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax ){
                                
                                Bool_t MotherFound = FindMother(track->GetLabel());
                                if(MotherFound){
                                    fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                                    if( TMath::Abs(pdg) == 11 && fMCparticleMother->GetPdgCode()!=22 ){
                                        if(fIsHFE1) fPt_track_match_num->Fill(fMCparticle->Pt());
                                    }
                                }
                            }
                        }
                    }///until here
                    
                    
                    
                    //----------------------------------------------------------------------------------------
                    //
                    //EtaCut electrons histogram
                    //Shower Shape Cut
                    if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax ){
                        
                        if(fUseShowerShapeCut){
                            if(M02 >= fM02CutMin && M02<=fM02CutMax && M20>=fM20CutMin && M20<=fM20CutMax){
                                fEoverP_pt[2]->Fill(fPt,(fClus->E() / fP));
                                fShowerShapeCut->Fill(M02,M20);
                                
                                if(fClus->E()> 8.0)fEoverP_pt_highE1->Fill(fPt,(fClus->E() / fP));
                                if(fClus->E()> 12.0)fEoverP_pt_highE0->Fill(fPt,(fClus->E() / fP));
                                
                                
                                //in order to check if there are exotic cluster in this selected cluster (27 may 2014)
                                //fNcells_energy_elec_selected->Fill(ncells,Energy);
                                
                                //to check how many tracks matches the cluster
                                Int_t tracks_matched = fClus->GetNTracksMatched();
                                fEoverP_ntracks_matched->Fill(fClus->E()/fP, tracks_matched);
                                fEoverP_ncells->Fill(fClus->E()/fP, ncells);
                                
                                //true electrons E/p shape
                                if(fIsMC && fIsAOD && track->GetLabel()>=0){
                                    fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
                                    Int_t pdg = fMCparticle->GetPdgCode();
                                    
                                    
                                    
                                    if( TMath::Abs(pdg) == 11){
                                        fEoverP_pt_true_electrons->Fill(fPt,(fClus->E() / fP));
                                    }
                                    
                                    //==========================E/p with weight==============================
                                    if( TMath::Abs(pdg) == 11){
                                        if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==22 || TMath::Abs(fMCparticleMother->GetPdgCode())==111 || TMath::Abs(fMCparticleMother->GetPdgCode())==221))
                                        {
                                            Double_t weight=1;
                                            //----------------------------------------------------------------------------
                                            if(TMath::Abs(fMCparticleMother->GetPdgCode())==111 || TMath::Abs(fMCparticleMother->GetPdgCode())==221){
                                                Double_t mPt=fMCparticleMother->Pt();
                                                
                                                
                                                if(TMath::Abs(fMCparticleMother->GetPdgCode())==111){
                                                    Double_t x=mPt;
                                                    weight=CalculateWeight(111, x);
                                                    
                                                }
                                                if(TMath::Abs(fMCparticleMother->GetPdgCode())==221){
                                                    Double_t x=mPt;
                                                    weight=CalculateWeight(221, x);
                                                    
                                                }
                                                fEoverP_pt_true_electrons_weight->Fill((fPt/weight),(fClus->E() / fP));
                                            }//// mother pion or eta
                                            else if(fMCparticleMother->GetMother()>0 && (TMath::Abs(fMCparticleGMother->GetPdgCode())==111 || TMath::Abs(fMCparticleGMother->GetPdgCode())==221 )){
                                                Double_t gmPt=fMCparticleGMother->Pt();
                                                
                                                
                                                if(TMath::Abs(fMCparticleGMother->GetPdgCode())==111){
                                                    Double_t x=gmPt;
                                                    weight=CalculateWeight(111, x);
                                                }
                                                
                                                
                                                if(TMath::Abs(fMCparticleGMother->GetPdgCode())==221){
                                                    Double_t x=gmPt;
                                                    weight=CalculateWeight(221, x);
                                                }
                                                fEoverP_pt_true_electrons_weight->Fill((fPt/weight),(fClus->E() / fP));
                                            }//grandmother pion or eta
                                            else{
                                                fEoverP_pt_true_electrons_weight->Fill((fPt/weight),(fClus->E() / fP));
                                            }
                                        }
                                    }
                                    //==========================E/p with weight==============================
                                    
                                    //==============
                                    //true HFE electrons
                                    Bool_t MotherFound = FindMother(track->GetLabel());
                                    if(MotherFound){
                                        
                                        if(fIsHFE1){
                                            if( TMath::Abs(pdg) == 11){
                                                fEoverP_pt_true_HFE->Fill(fPt,(fClus->E() / fP));
                                            }
                                        }
                                        if(!fIsHFE1){
                                            if( TMath::Abs(pdg) == 11){
                                                fEoverP_pt_not_HFE->Fill(fPt,(fClus->E() / fP));
                                            }
                                        }
                                        
                                    }
                                    //==============
                                    
                                }
                                
                                
                                
                            }
                            
                        }
                        if(!fUseShowerShapeCut){
                            fEoverP_pt[2]->Fill(fPt,(fClus->E() / fP));
                            fShowerShapeCut->Fill(M02,M20);
                            
                            //to check how many tracks matches the cluster
                            Int_t tracks_matched = fClus->GetNTracksMatched();
                            fEoverP_ntracks_matched->Fill(fClus->E()/fP, tracks_matched);
                            
                            fEoverP_ncells->Fill(fClus->E()/fP, ncells);
                            
                            //after cut -> Are they really electrons ? E/p shape
                            if(fIsMC && fIsAOD && track->GetLabel()>=0){
                                fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
                                Int_t pdg = fMCparticle->GetPdgCode();
                                
                                if( TMath::Abs(pdg) == 11){
                                    fEoverP_pt_true_electrons->Fill(fPt,(fClus->E() / fP));
                                }
                                
                                
                                //==============
                                //true HFE electrons
                                Bool_t MotherFound = FindMother(track->GetLabel());
                                if(MotherFound){
                                    
                                    if(fIsHFE1){
                                        if( TMath::Abs(pdg) == 11){
                                            fEoverP_pt_true_HFE->Fill(fPt,(fClus->E() / fP));
                                        }
                                    }
                                    if(!fIsHFE1){
                                        if( TMath::Abs(pdg) == 11){
                                            fEoverP_pt_not_HFE->Fill(fPt,(fClus->E() / fP));
                                        }
                                    }
                                    
                                }
                                //==============
                                
                                
                                
                            }
                            
                            
                        }
                        if(fUseEMCal) fShowerShape_ele->Fill(M02,M20);
                        
                        //for shower shape cut studies - now with TPC PID Cut
                        if(fUseEMCal){
                            if(fEMCEG1){
                                if(fClus->E()>=12.00){
                                    fShowerShapeM02_EoverP->Fill(M02,EoverP);
                                    fShowerShapeM20_EoverP->Fill(M20,EoverP);
                                    fNcells_energy_elec_selected->Fill(ncells,Energy);
                                    
                                }
                                
                            }
                            if(fEMCEG2){
                                if(fClus->E()>=8.00){
                                    fShowerShapeM02_EoverP->Fill(M02,EoverP);
                                    fShowerShapeM20_EoverP->Fill(M20,EoverP);
                                    fNcells_energy_elec_selected->Fill(ncells,Energy);
                                    
                                }
                                
                            }
                            if(!fEMCEG1 && ! fEMCEG2){
                                fShowerShapeM02_EoverP->Fill(M02,EoverP);
                                fShowerShapeM20_EoverP->Fill(M20,EoverP);
                                fNcells_energy_elec_selected->Fill(ncells,Energy);
                                
                                
                            }
                            if(fIsMC){
                                if(fClus->E()>=12.00){
                                    fShowerShapeM02_EoverP_Ecut12_MC->Fill(M02,EoverP);
                                    fShowerShapeM20_EoverP_Ecut12_MC->Fill(M20,EoverP);
                                    
                                }
                                if(fClus->E()>=8.00){
                                    fShowerShapeM02_EoverP_Ecut8_MC->Fill(M02,EoverP);
                                    fShowerShapeM20_EoverP_Ecut8_MC->Fill(M20,EoverP);
                                    
                                }
                                
                                
                            }
                            
                        }
                        
                    }
                    
                    //----------------------------------------------------------------------------------------
                    
                    
                    
                    // for Eta Phi distribution
                    fClus->GetPosition(pos3);
                    TVector3 vpos(pos3[0],pos3[1],pos3[2]);
                    Double_t cphi = vpos.Phi();
                    Double_t ceta = vpos.Eta();
                    fEtaPhi[2]->Fill(cphi,ceta);
                    
                    if(EoverP>1.2)fEtaPhi_large_EoverP->Fill(cphi,ceta);
                    if(EoverP<1.2)fEtaPhi_small_EoverP->Fill(cphi,ceta);
                    
                    fTPCNcls_EoverP[2]->Fill(TPCNcls, EoverP);
                    
                    for(Int_t i = 0; i < 6; i++)
                    {
                        if(fPt>=fPtBin[i] && fPt<fPtBin[i+1])
                        {
                            
                            fR_EoverP[i]->Fill(R, EoverP);
                            fNcells[i]->Fill(ncells);
                            fNcells_EoverP[i]->Fill(EoverP, ncells);
                            fM02_EoverP[i]->Fill(M02,EoverP);
                            fM20_EoverP[i]->Fill(M20,EoverP);
                            fECluster_ptbins[i]->Fill(Energy);
                            fEoverP_ptbins[i]->Fill(EoverP);
                            
                            if((fClus->E() / fP) >= fEoverPCutMin && (fClus->E() / fP) <= fEoverPCutMax) {
                                fNcells_electrons[i]->Fill(ncells);
                            }
                            
                            if(M02<0.5 && M20<0.3) {
                                fEoverP_wSSCut[i]->Fill(EoverP);
                            }
                        }
                    }
                    
                    fNcells_pt->Fill(fPt, ncells);
                    
                    
                    ////////////////////////////////////////////////////////////////////
                    ///EMCal - Efficiency calculations
                    
                    
                    if(fIsMC && fIsAOD && track->GetLabel()>=0)
                    {
                        fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
                        Int_t pdg = fMCparticle->GetPdgCode();
                        
                        //
                        if(fMCparticle->IsPhysicalPrimary()){
                            
                            //Phi cut && fMCparticle->Phi()>=(TMath::Pi()*80/180) && fMCparticle->Phi()<=TMath::Pi()
                            if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax  ){
                                
                                Bool_t MotherFound = FindMother(track->GetLabel());
                                if(MotherFound){
                                    fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                                    if( TMath::Abs(pdg) == 11 && fMCparticleMother->GetPdgCode()!=22 ){
                                        if(fIsHFE1)fPtMC_EMCal_All->Fill(fMCparticle->Pt());
                                    }
                                }
                            }
                        }
                    }
                    
                    else if(fIsMC && track->GetLabel()>=0)//ESD
                    {
                        
                        if(fMCstack->IsPhysicalPrimary(track->GetLabel()))
                        {
                            
                            fMCtrack = fMCstack->Particle(track->GetLabel());
                            
                            Int_t pdg = fMCtrack->GetPdgCode();
                            //Phi cut && fMCtrack->Phi()>=(TMath::Pi()*80/180) && fMCtrack->Phi()<=TMath::Pi()
                            if(fMCtrack->Eta()>=fEtaCutMin && fMCtrack->Eta()<=fEtaCutMax )
                            {
                                Bool_t MotherFound = FindMother(track->GetLabel());
                                //MotherFound included
                                if(MotherFound){
                                    if(fMCtrack->GetFirstMother()>0){
                                        fMCtrackMother = fMCstack->Particle(fMCtrack->GetFirstMother());
                                        if( TMath::Abs(pdg) == 11 && fMCtrackMother->GetPdgCode()!=22 ){
                                            if(fIsHFE1)fPtMC_EMCal_All->Fill(fMCtrack->Pt());
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    //_______________________________________________________
                    //PID using EMCal
                    if(fEoverPnsigma)SetEoverPCutPtDependentMC(track->Pt());
                    
                    if((fClus->E() / fP) >= fEoverPCutMin && (fClus->E() / fP) <= fEoverPCutMax)
                    {
                        
                        
                        fECluster[2]->Fill(Energy);
                        fTPCNcls_pid[3]->Fill(TPCNcls, TPCNcls_pid);
                        
                        if(track->Pt()>=12)fECluster_highpT0->Fill(Energy);
                        if(track->Pt()>=8)fECluster_highpT1->Fill(Energy);
                        
                        
                        
                        
                        if(fUseEMCal)
                        {
                            if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax ){
                                fPtElec_Inc->Fill(fPt);
                                //eta phi distribution for data
                                fEtaPhi_data->Fill(track->Phi(),track->Eta());
                                
                                
                                //studies to see how is the clusters in trigger and MB for a given selected cluster after all cuts
                                //for the 7 GeV threshold
                                if(fClus->E() >= 8.0 && fClus->E() <= 16.0){
                                    fNcells_Energy_highE1->Fill(Energy, ncells);
                                    fEtaPhi_highE1->Fill(cphi,ceta);
                                    fEta_highE1->Fill(Dz);
                                    fPhi_highE1->Fill(Dx);
                                    fR_highE1->Fill(R);
                                    fNCluster_highE1->Fill(ClsNo);
                                }
                                
                                //for the 11 GeV threshold
                                if(fClus->E() >= 12.0 && fClus->E() <= 20.0){
                                    fNcells_Energy_highE0->Fill(Energy, ncells);
                                    fEtaPhi_highE0->Fill(cphi,ceta);
                                    fEta_highE0->Fill(Dz);
                                    fPhi_highE0->Fill(Dx);
                                    fR_highE0->Fill(R);
                                    fNCluster_highE0->Fill(ClsNo);
                                }
                                
                                
                                
                                
                            }
                            
                            //Eta cut for background
                            if(fFillBackground){
                                fEtad[2]->Fill(track->Eta());
                                
                                //background for triggered data: trigger electron must have same cuts on shower shape  06/Jan/2014
                                if(fUseShowerShapeCut){
                                    if(M02 >= fM02CutMin && M02<=fM02CutMax && M20>=fM20CutMin && M20<=fM20CutMax){
                                        if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax){
                                            Background(track, iTracks, Vtrack, kFALSE, kFALSE, kFALSE);
                                        }
                                    }
                                }
                                else{
                                    if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax){
                                        Background(track, iTracks, Vtrack, kFALSE, kFALSE, kFALSE);
                                    }
                                }
                                
                            }
                            
                            double emctof2 = fClus->GetTOF();
                            ftimingEle2->Fill(fPt,1e9*emctof2);
                            //Correlation Analysis
                            if(fCorrelationFlag)
                            {
                                ElectronHadronCorrelation(track, iTracks, Vtrack);
                            }
                        }
                        //_______________________________________________________
                        
                        ////////////////////////////////////////////////////////////////////
                        ///EMCal - efficiency calculations
                        
                        if(track->Charge()<0)  fCharge_n->Fill(fPt);
                        if(track->Charge()>0)  fCharge_p->Fill(fPt);
                        
                        
                        
                        if(fIsMC && fIsAOD && track->GetLabel()>=0)//AOD
                        {
                            if(track->Charge()<0)  fCharge_n->Fill(fPt);
                            if(track->Charge()>0)  fCharge_p->Fill(fPt);
                            
                            fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
                            if(fMCparticle->GetMother()>0) fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                            Int_t pdg = fMCparticle->GetPdgCode();
                            
                            double proX = fMCparticle->Xv();
                            double proY = fMCparticle->Yv();
                            double proR = sqrt(pow(proX,2)+pow(proY,2));
                            
                            
                            if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax && fMCparticle->Phi()>=(TMath::Pi()*80/180) && fMCparticle->Phi()<=TMath::Pi()  ){
                                
                                if( TMath::Abs(pdg) == 11 && fMCparticle->GetMother()>0 ){
                                    Int_t mpdg = fMCparticleMother->GetPdgCode();
                                    if(TMath::Abs(mpdg) == 221 || TMath::Abs(mpdg) == 22 || TMath::Abs(mpdg) == 111){
                                        if(proR<7)fPtMCelectronAfterAll_nonPrimary->Fill(fMCparticle->Pt()); //numerator for the total efficiency, non Primary track
                                    }
                                }
                                if( TMath::Abs(pdg) == 11 && fMCparticle->IsPhysicalPrimary()) fPtMCelectronAfterAll_Primary->Fill(fMCparticle->Pt());
                            }
                            
                            
                            //
                            if(fMCparticle->IsPhysicalPrimary()){
                                
                                
                                if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax ){
                                    
                                    Bool_t MotherFound = FindMother(track->GetLabel());
                                    if(MotherFound){
                                        
                                        if(!fUseShowerShapeCut){
                                            if(fIsHFE1){
                                                //Unfolding   pt_reco/pt_MC  in the efficiency
                                                fPtMCelectronAfterAll->Fill(fMCparticle->Pt());
                                                fPtMCelectronAfterAll_unfolding->Fill(track->Pt());
                                                fEtaPhi_num->Fill(fMCparticle->Phi(),fMCparticle->Eta());
                                                
                                                //new histo to estimate how different is pt reco from pt MC
                                                fpt_reco_pt_MC_num->Fill(track->Pt(),fMCparticle->Pt());
                                            }//numerator for the total efficiency AOD
                                        }
                                        //November 11 for efficiency of triggered data
                                        if(fUseShowerShapeCut){
                                            if(M02 >= fM02CutMin && M02<=fM02CutMax && M20>=fM20CutMin && M20<=fM20CutMax){
                                                if(fIsHFE1){
                                                    //Unfolding   pt_reco/pt_MC  in the efficiency
                                                    fPtMCelectronAfterAll->Fill(fMCparticle->Pt());
                                                    fPtMCelectronAfterAll_unfolding->Fill(track->Pt());
                                                    fEtaPhi_num->Fill(fMCparticle->Phi(),fMCparticle->Eta());
                                                    
                                                    //new histo to estimate how different is pt reco from pt MC
                                                    fpt_reco_pt_MC_num->Fill(track->Pt(),fMCparticle->Pt());
                                                }//numerator for the total efficiency AOD
                                            }
                                        }
                                        
                                        
                                        
                                        fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                                        if( TMath::Abs(pdg) == 11 && fMCparticleMother->GetPdgCode()!=22 ){
                                            if(fIsHFE1)fPtMC_EMCal_Selected->Fill(fMCparticle->Pt());
                                        }
                                    }//if MotherFound
                                }//eta cut
                            }
                        }///close AOD
                        
                        else if(fIsMC && track->GetLabel()>=0)//ESD
                        {
                            if(track->Charge()<0)  fCharge_n->Fill(fPt);
                            if(track->Charge()>0)  fCharge_p->Fill(fPt);
                            
                            fMCtrack = fMCstack->Particle(track->GetLabel());
                            if(fMCtrack->GetFirstMother()>0)
                            {
                                fMCtrackMother = fMCstack->Particle(fMCtrack->GetFirstMother());
                            }
                            TParticle *particle=fMCstack->Particle(track->GetLabel());
                            
                            Int_t pdg = fMCtrack->GetPdgCode();
                            
                            
                            if(fMCtrack->Eta()>=fEtaCutMin && fMCtrack->Eta()<=fEtaCutMax)
                            {
                                if( TMath::Abs(pdg) == 11 && fMCtrack->GetFirstMother()>0 )
                                {
                                    Int_t mpdg = fMCtrackMother->GetPdgCode();
                                    if(TMath::Abs(mpdg) == 221 || TMath::Abs(mpdg) == 22 || TMath::Abs(mpdg) == 111)
                                    {
                                        Double_t proR=particle->R();
                                        if(proR<7)
                                        {
                                            fPtMCelectronAfterAll_nonPrimary->Fill(fMCtrack->Pt()); //numerator for the total efficiency, non Primary track
                                        }
                                    }
                                }
                                if( TMath::Abs(pdg) == 11 && fMCstack->IsPhysicalPrimary(track->GetLabel()))
                                {
                                    fPtMCelectronAfterAll_Primary->Fill(fMCtrack->Pt());
                                }
                            }
                            
                            if(fMCstack->IsPhysicalPrimary(track->GetLabel()))
                            {
                                Bool_t MotherFound = FindMother(track->GetLabel());
                                
                                if(MotherFound)
                                {
                                    if(fMCtrack->Eta()>=fEtaCutMin && fMCtrack->Eta()<=fEtaCutMax)
                                    {
                                        if(!fUseShowerShapeCut){
                                            if(fIsHFE1){
                                                fPtMCelectronAfterAll->Fill(fMCtrack->Pt()); //numerator for the total efficiency ESD
                                                fEtaPhi_num->Fill(fMCtrack->Phi(),fMCtrack->Eta());
                                            }
                                        }
                                        //November 11 for efficiency of triggered data
                                        if(fUseShowerShapeCut){
                                            if(M02 >= fM02CutMin && M02<=fM02CutMax && M20>=fM20CutMin && M20<=fM20CutMax){
                                                if(fIsHFE1){
                                                    fPtMCelectronAfterAll->Fill(fMCtrack->Pt()); //numerator for the total efficiency ESD
                                                    fEtaPhi_num->Fill(fMCtrack->Phi(),fMCtrack->Eta());
                                                }
                                            }
                                        }
                                        
                                        
                                        
                                    }
                                }
                                
                                
                                
                                // Phi cut && fMCtrack->Phi()>=(TMath::Pi()*80/180) && fMCtrack->Phi()<=TMath::Pi()
                                if(fMCtrack->Eta()>=fEtaCutMin && fMCtrack->Eta()<=fEtaCutMax)
                                {
                                    //included MotherFound
                                    
                                    if(MotherFound)
                                    {
                                        if(fMCtrack->GetFirstMother()>0){
                                            fMCtrackMother = fMCstack->Particle(fMCtrack->GetFirstMother());
                                            if( TMath::Abs(pdg) == 11 && fMCtrackMother->GetPdgCode()!=22 ){
                                                
                                                if(fIsHFE1){fPtMC_EMCal_Selected->Fill(fMCtrack->Pt());}
                                            }
                                        }
                                    }
                                }
                            }
                        }//close ESD
                        ///////////////////////////////////////////////////////////////////
                        
                        
                    }
                }
            }
        }
        
        //______________________________________________________________
        // Vertex
        
        fVtxZ[2]->Fill(fZvtx);
        if(iTracks == 0)fNTracks[2]->Fill(NTracks);
        fNTracks_pt[2]->Fill(NTracks, fPt);
        fNTracks_eta[2]->Fill(NTracks, track->Eta());
        fNTracks_phi[2]->Fill(NTracks, track->Phi());
        fNClusters[2]->Fill(ClsNo);
        fTPCNcls_pid[2]->Fill(TPCNcls, TPCNcls_pid);
        
        //______________________________________________________________
        
        //_______________________________________________________
        //Correlation Analysis
        
        if(fIsMC  && atrack->GetLabel()>=0)
        {
            if(fIsAOD)
            {
                fMCparticle = (AliAODMCParticle*) fMCarray->At(atrack->GetLabel());
                Int_t pdg = fMCparticle->GetPdgCode();
                if (TMath::Abs(pdg) == 11)
                {
                    if (fMCparticle->IsPrimary())
                        fMCEffPID_afterPID->Fill(fMCparticle->Pt());
                }
                
            }
        }
        
        
        if(!fUseEMCal)
        {
            
            fPtElec_Inc->Fill(fPt);
            
            if(fCorrelationFlag)
            {
                ElectronHadronCorrelation(track, iTracks, Vtrack);
            }
        }
        //_______________________________________________________
        
        ///________________________________________________________________________
    }
    
    //__________________________________________________________________
    //Event Mixing Analysis
    //Filling pool
    if(fEventMixingFlag)
    {
        if(fIspp)
        {
            fPool = fPoolMgr->GetEventPool(1.5, fZvtx); // Get the buffer associated with the current centrality and z-vtx
            
            //if(!fPool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fCentrality->GetCentralityPercentile("V0A"), fZvtx));
            if(!fPool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f", 1.5, fZvtx));
            
            fPool->UpdatePool(SelectedHadrons());
        }
        else
        {
            fPool = fPoolMgr->GetEventPool(fCentrality->GetCentralityPercentile("V0A"), fZvtx); // Get the buffer associated with the current centrality and z-vtx
            
            if(!fPool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fCentrality->GetCentralityPercentile("V0A"), fZvtx));
            
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


void AliAnalysisTaskHFEpACorrelation::Background(AliVTrack *track, Int_t trackIndex, AliVParticle *vtrack, Bool_t IsTPConly, Bool_t IsWeight, Bool_t MassPtBins)
{
    ///_________________________________________________________________
    ///MC analysis
    //Bool_t IsMCefix=kFALSE; //to make correction on efix, use kTRUE (do not change the efficiency, so I will keep the correction only for d3)
    
    if(fIsMC)
    {
        if(track->GetLabel() < 0)
        {
            AliWarning(Form("The track %d does not have a valid MC label",trackIndex));
            return;
        }
        
        if(fIsAOD)
        {
            fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
            
            if(fMCparticle->GetMother()<0) return;
            
            fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
            if(fMCparticleMother->GetMother()>0)fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
            
            if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==22 || TMath::Abs(fMCparticleMother->GetPdgCode())==111 || TMath::Abs(fMCparticleMother->GetPdgCode())==221))
            {
                //Is Background
                if(!IsTPConly)fPtBackgroundBeforeReco->Fill(track->Pt());
                if(IsTPConly)fPtBackgroundBeforeReco2->Fill(track->Pt());
                
                
                //October 08th weighted histos
                if(TMath::Abs(fMCparticleMother->GetPdgCode())==111 || TMath::Abs(fMCparticleMother->GetPdgCode())==221 ){
                    
                    Double_t mPt=fMCparticleMother->Pt();
                    Double_t mweight=1;
                    
                    
                    //________________________________________________________________
                    //correction for d3 based on data
                    
                    if(TMath::Abs(fMCparticleMother->GetPdgCode())==111){
                        Double_t x=mPt;
                        mweight=CalculateWeight(111, x);
                        
                    }
                    if(TMath::Abs(fMCparticleMother->GetPdgCode())==221){
                        Double_t x=mPt;
                        mweight=CalculateWeight(221, x);
                        
                    }
                    
                    //________________________________________________________________
                    
                    //Histo pT mother versus pT electron
                    fpT_m_electron->Fill(mPt, track->Pt());
                    
                    if(!IsTPConly)fPtBackgroundBeforeReco_weight->Fill(track->Pt(), 1./mweight);
                    if(IsTPConly)fPtBackgroundBeforeReco2_weight->Fill(track->Pt(), 1./mweight);
                }
                else if(fMCparticleMother->GetMother()>0 && (TMath::Abs(fMCparticleGMother->GetPdgCode())==111 || TMath::Abs(fMCparticleGMother->GetPdgCode())==221 )){
                    
                    Double_t gmPt=fMCparticleGMother->Pt();
                    Double_t gmweight=1;
                    
                    //________________________________________________________________
                    //correction for d3 based on data
                    
                    if(TMath::Abs(fMCparticleGMother->GetPdgCode())==111){
                        Double_t x=gmPt;
                        gmweight=CalculateWeight(111, x);
                    }
                    if(TMath::Abs(fMCparticleGMother->GetPdgCode())==221){
                        Double_t x=gmPt;
                        gmweight=CalculateWeight(221, x);
                    }
                    
                    
                    //________________________________________________________________
                    //Histo pT gmother versus pT electron
                    
                    fpT_gm_electron->Fill(gmPt, track->Pt());
                    
                    if(!IsTPConly)fPtBackgroundBeforeReco_weight->Fill(track->Pt(), 1./gmweight);
                    if(IsTPConly)fPtBackgroundBeforeReco2_weight->Fill(track->Pt(), 1./gmweight);
                }
                else{
                    if(!IsTPConly)fPtBackgroundBeforeReco_weight->Fill(track->Pt());
                    if(IsTPConly)fPtBackgroundBeforeReco2_weight->Fill(track->Pt());
                }
            }//particle kind
        }//IsAOD
        //ESD
        else
        {
            fMCtrack = fMCstack->Particle(track->GetLabel());
            
            if(fMCtrack->GetFirstMother()<0) return;
            
            fMCtrackMother = fMCstack->Particle(fMCtrack->GetFirstMother());
            
            if(TMath::Abs(fMCtrack->GetPdgCode())==11 && (TMath::Abs(fMCtrackMother->GetPdgCode())==22 || TMath::Abs(fMCtrackMother->GetPdgCode())==111 || TMath::Abs(fMCtrackMother->GetPdgCode())==221))
            {
                //Is Background
                if(!IsTPConly)fPtBackgroundBeforeReco->Fill(track->Pt());
                if(IsTPConly)fPtBackgroundBeforeReco2->Fill(track->Pt());
            }
        }
    }//IsMC
    
    ///_________________________________________________________________
    
    //________________________________________________
    //Associated particle cut
    fPartnerCuts->SetAcceptKinkDaughters(kFALSE);
    fPartnerCuts->SetRequireITSRefit(kTRUE);
    fPartnerCuts->SetRequireTPCRefit(kTRUE);
    fPartnerCuts->SetEtaRange(-0.9,0.9);
    fPartnerCuts->SetMaxChi2PerClusterTPC(4.0);
    fPartnerCuts->SetMinNClustersTPC(80);
    fPartnerCuts->SetPtRange(0,1e10);
    //fPartnerCuts->SetRequireSigmaToVertex(kTRUE);
    //fPartnerCuts->SetMaxDCAToVertexXY(1);
    //fPartnerCuts->SetMaxDCAToVertexZ(3);
    //_________________________________________________
    
    ///#################################################################
    //Non-HFE reconstruction
    fNonHFE = new AliSelectNonHFE();
    fNonHFE->SetAODanalysis(fIsAOD);
    if(fMassCutFlag) fNonHFE->SetInvariantMassCut(fMassCut);
    if(fAngleCutFlag) fNonHFE->SetOpeningAngleCut(fAngleCut);
    if(fChi2CutFlag) fNonHFE->SetChi2OverNDFCut(fChi2Cut);
    if(fDCAcutFlag) fNonHFE->SetDCACut(fDCAcut);
    fNonHFE->SetAlgorithm("DCA"); //KF
    fNonHFE->SetPIDresponse(fPidResponse);
    fNonHFE->SetTrackCuts(-3.5,3.5,fPartnerCuts);
    fNonHFE->SetAdditionalCuts(fPtMinAsso,fTpcNclsAsso);
    
    //Electron Information
    Double_t fPhiE = -999;
    Double_t fEtaE = -999;
    Double_t fPtE = -999;
    fPhiE = track->Phi();
    fEtaE = track->Eta();
    fPtE = track->Pt();
    
    
    
    
    if(!IsTPConly && !IsWeight){
        fNonHFE->SetHistAngleBack(fOpAngleBack);
        fNonHFE->SetHistAngle(fOpAngle);
        fNonHFE->SetHistDCABack(fDCABack);
        fNonHFE->SetHistDCA(fDCA);
        fNonHFE->SetHistMassBack(fInvMassBack);
        fNonHFE->SetHistMass(fInvMass);
        
    }
    
    
    if(IsTPConly && !IsWeight){
        fNonHFE->SetHistAngleBack(fOpAngleBack2);
        fNonHFE->SetHistAngle(fOpAngle2);
        fNonHFE->SetHistDCABack(fDCABack2);
        fNonHFE->SetHistDCA(fDCA2);
        fNonHFE->SetHistMassBack(fInvMassBack2);
        fNonHFE->SetHistMass(fInvMass2);
        
    }
    /*
     //to look the Invariant mass in pT bins
     if(MassPtBins && !IsWeight){
     Double_t fPtBin_trigger2[11] = {1,2,4,6,8,10,12,14,16,18,20};
     for(Int_t i = 0; i < 10; i++)
     {
     if(fPtE>=fPtBin_trigger2[i] && fPtE<fPtBin_trigger2[i+1])
     {
					fNonHFE->SetHistMassBack(fInvMassBack_pT[i]);
					fNonHFE->SetHistMass(fInvMass_pT[i]);
     }
     }
     }
     //end of Invariant mass in pT bins
     */
    
    
    
    //included to apply weight in the invariant mass distribution  21 May 2015
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    /*
     if(fIsMC)
     {
     if(fIsAOD && IsTPConly && IsWeight)
     {
     if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==22 || TMath::Abs(fMCparticleMother->GetPdgCode())==111 || TMath::Abs(fMCparticleMother->GetPdgCode())==221))
     {
     
					Double_t weight=1;
					//----------------------------------------------------------------------------
					if(TMath::Abs(fMCparticleMother->GetPdgCode())==111 || TMath::Abs(fMCparticleMother->GetPdgCode())==221){
     Double_t mPt=fMCparticleMother->Pt();
     
     
     if(TMath::Abs(fMCparticleMother->GetPdgCode())==111){
     Double_t x=mPt;
     weight=CalculateWeight(111, x);
     
     }
     if(TMath::Abs(fMCparticleMother->GetPdgCode())==221){
     Double_t x=mPt;
     weight=CalculateWeight(221, x);
     
     }
     
     
     fNonHFE->SetHistMassBack(fInvMassBack2_weight,1./weight);
     fNonHFE->SetHistMass(fInvMass2_weight,1./weight);
     
					}// mother pion or eta
					//----------------------------------------------------------------------------
					else if(fMCparticleMother->GetMother()>0 && (TMath::Abs(fMCparticleGMother->GetPdgCode())==111 || TMath::Abs(fMCparticleGMother->GetPdgCode())==221 )){
     Double_t gmPt=fMCparticleGMother->Pt();
     
     
     if(TMath::Abs(fMCparticleGMother->GetPdgCode())==111){
     Double_t x=gmPt;
     weight=CalculateWeight(111, x);
     }
     
     
     if(TMath::Abs(fMCparticleGMother->GetPdgCode())==221){
     Double_t x=gmPt;
     weight=CalculateWeight(221, x);
     }
     
     fNonHFE->SetHistMassBack(fInvMassBack2_weight,1./weight);
     fNonHFE->SetHistMass(fInvMass2_weight,1./weight);
     
     
					}//grandmother pion or eta
					//----------------------------------------------------------------------------
					else{
     fNonHFE->SetHistMassBack(fInvMassBack2_weight,1./weight);
     fNonHFE->SetHistMass(fInvMass2_weight,1./weight);
					}
     
     
     }
     }
     }
     // end of part included to apply weight in the invariant mass distribution  21 May 2015
     //----------------------------------------------------------------------------
     //----------------------------------------------------------------------------
     */
    
    fNonHFE->FindNonHFE(trackIndex,vtrack,fVevent);
    
    //index of track selected as partner
    Int_t *fUlsPartner = fNonHFE->GetPartnersULS();
    
    
    ///_________________________________________________________________
    ///MC analysis
    if(fIsMC)
    {
        if(fIsAOD)
        {
            if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==22 || TMath::Abs(fMCparticleMother->GetPdgCode())==111 || TMath::Abs(fMCparticleMother->GetPdgCode())==221))
            {
                
                Double_t weight=1;
                
                if(!IsTPConly){
                    if(fNonHFE->IsULS()) fPtElec_ULS->Fill(fPtE,fNonHFE->GetNULS());
                    if(fNonHFE->IsLS()) fPtElec_LS->Fill(fPtE,fNonHFE->GetNLS());
                    
                    
                    
                    
                    if(TMath::Abs(fMCparticleMother->GetPdgCode())==111 || TMath::Abs(fMCparticleMother->GetPdgCode())==221){
                        Double_t mPt=fMCparticleMother->Pt();
                        Double_t mweight1=1;
                        Double_t mweight2=1;
                        //Double_t weight=1;
                        
                        
                        //----------------------------------------------------------------------------
                        //correction based on data only
                        if(TMath::Abs(fMCparticleMother->GetPdgCode())==111){
                            Double_t x=mPt;
                            weight=CalculateWeight(111, x);
                        }
                        if(TMath::Abs(fMCparticleMother->GetPdgCode())==221){
                            Double_t x=mPt;
                            weight=CalculateWeight(221, x);
                        }
                        
                        
                        //----------------------------------------------------------------------------
                        
                        //check this
                        if(fNonHFE->IsULS()) mweight1=(fNonHFE->GetNULS())/weight;
                        if(fNonHFE->IsLS())  mweight2=(fNonHFE->GetNLS())/weight;
                        
                        //fill histos
                        if(fNonHFE->IsULS())fPtElec_ULS_weight->Fill(fPtE, mweight1);
                        if(fNonHFE->IsLS())fPtElec_LS_weight->Fill(fPtE, mweight2);
                    }
                    else if(fMCparticleMother->GetMother()>0 && (TMath::Abs(fMCparticleGMother->GetPdgCode())==111 || TMath::Abs(fMCparticleGMother->GetPdgCode())==221 )){
                        Double_t gmPt=fMCparticleGMother->Pt();
                        Double_t gmweight1=1;
                        Double_t gmweight2=1;
                        //Double_t weight=1;
                        
                        //----------------------------------------------------------------------------
                        
                        //----------------------------------------------------------------------------
                        
                        //correction based on data only for pi0
                        
                        if(TMath::Abs(fMCparticleGMother->GetPdgCode())==111){
                            Double_t x=gmPt;
                            weight=CalculateWeight(111, x);
                        }
                        if(TMath::Abs(fMCparticleGMother->GetPdgCode())==221){
                            Double_t x=gmPt;
                            weight=CalculateWeight(221, x);
                        }
                        
                        
                        
                        
                        
                        //check this
                        if(fNonHFE->IsULS()) gmweight1=(fNonHFE->GetNULS())/weight;
                        if(fNonHFE->IsLS())  gmweight2=(fNonHFE->GetNLS())/weight;
                        
                        //fill histos
                        if(fNonHFE->IsULS())fPtElec_ULS_weight->Fill(fPtE, gmweight1);
                        if(fNonHFE->IsLS())fPtElec_LS_weight->Fill(fPtE, gmweight2);
                    }
                    else{
                        if(fNonHFE->IsULS()) fPtElec_ULS_weight->Fill(fPtE,fNonHFE->GetNULS());
                        if(fNonHFE->IsLS()) fPtElec_LS_weight->Fill(fPtE,fNonHFE->GetNLS());
                    }
                    
                    
                }//!IsTPConly
                
                if(IsTPConly){
                    if(fNonHFE->IsULS()) fPtElec_ULS2->Fill(fPtE,fNonHFE->GetNULS());
                    if(fNonHFE->IsLS()) fPtElec_LS2->Fill(fPtE,fNonHFE->GetNLS());
                    
                    //new 08 October	//weighted histograms
                    if(TMath::Abs(fMCparticleMother->GetPdgCode())==111 || TMath::Abs(fMCparticleMother->GetPdgCode())==221){
                        Double_t mPt=fMCparticleMother->Pt();
                        
                        Double_t mweight1=1;
                        Double_t mweight2=1;
                        //Double_t weight=1;
                        
                        
                        //----------------------------------------------------------------------------
                        
                        //correction based on data for d3
                        
                        if(TMath::Abs(fMCparticleMother->GetPdgCode())==111){
                            Double_t x=mPt;
                            weight=CalculateWeight(111, x);
                            
                        }
                        if(TMath::Abs(fMCparticleMother->GetPdgCode())==221){
                            Double_t x=mPt;
                            weight=CalculateWeight(221, x);
                            
                        }
                        
                        
                        
                        
                        //check this
                        if(fNonHFE->IsULS()) mweight1=(fNonHFE->GetNULS())/weight;
                        if(fNonHFE->IsLS())  mweight2=(fNonHFE->GetNLS())/weight;
                        
                        //fill histos
                        if(fNonHFE->IsULS())fPtElec_ULS2_weight->Fill(fPtE, mweight1);
                        if(fNonHFE->IsLS())fPtElec_LS2_weight->Fill(fPtE, mweight2);
                        
                        
                        
                        
                    }
                    else if(fMCparticleMother->GetMother()>0 && (TMath::Abs(fMCparticleGMother->GetPdgCode())==111 || TMath::Abs(fMCparticleGMother->GetPdgCode())==221 )){
                        Double_t gmPt=fMCparticleGMother->Pt();
                        Double_t gmweight1=1;
                        Double_t gmweight2=1;
                        //Double_t weight=1;
                        
                        
                        //----------------------------------------------------------------------------
                        //correction based on data only for pi0
                        
                        if(TMath::Abs(fMCparticleGMother->GetPdgCode())==111){
                            Double_t x=gmPt;
                            weight=CalculateWeight(111, x);
                        }
                        
                        
                        if(TMath::Abs(fMCparticleGMother->GetPdgCode())==221){
                            Double_t x=gmPt;
                            weight=CalculateWeight(221, x);
                        }
                        
                        //----------------------------------------------------------------------------
                        
                        
                        
                        
                        
                        //check this
                        if(fNonHFE->IsULS()) gmweight1=(fNonHFE->GetNULS())/weight;
                        if(fNonHFE->IsLS())  gmweight2=(fNonHFE->GetNLS())/weight;
                        
                        //fill histos
                        if(fNonHFE->IsULS())fPtElec_ULS2_weight->Fill(fPtE, gmweight1);
                        if(fNonHFE->IsLS())fPtElec_LS2_weight->Fill(fPtE, gmweight2);
                        
                        
                        
                    }
                    else{
                        
                        
                        
                        if(fNonHFE->IsULS()) fPtElec_ULS2_weight->Fill(fPtE,fNonHFE->GetNULS());
                        if(fNonHFE->IsLS()) fPtElec_LS2_weight->Fill(fPtE,fNonHFE->GetNLS());
                        
                        
                        
                    }
                    
                    
                    //----------------------------------------------------------------------------
                    //to check other way to calculate efficiency
                    //ULS with no weight from ULS-LS original
                    //we have to know if track2 comes from same mother!!!
                    //----------------------------------------------------------------------------
                    if(fNonHFE->IsULS()){
                        
                        for(Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++)
                        {
                            
                            AliVParticle* Vtrack2 = fVevent->GetTrack(iTracks);
                            if (!Vtrack2)
                            {
                                printf("ERROR: Could not receive track %d\n", iTracks);
                                continue;
                            }
                            AliVTrack *track2 = dynamic_cast<AliVTrack*>(Vtrack2);
                            if(track2->GetLabel()<0) continue;
                            fMCparticle2 = (AliAODMCParticle*) fMCarray->At(track2->GetLabel());
                            if(fMCparticle2->GetMother()<0) continue;
                            
                            for(Int_t i = 0; i < fNonHFE->GetNULS(); i++)
                            {
                                if(fUlsPartner[i]==iTracks){
                                    //only fill if it has same mother
                                    //with weight to take into account the number of partners
                                    if(fMCparticle2->GetMother()==fMCparticle->GetMother()) fPtElec_ULS_MC->Fill(fPtE, fNonHFE->GetNULS());
                                    
                                    //-----------------------------------------------------------------------------------------------------------
                                    //weight for mother
                                    //Double_t weight2=1;
                                    Double_t mPt=fMCparticleMother->Pt();
                                    
                                    
                                    if(TMath::Abs(fMCparticleMother->GetPdgCode())==111){
                                        Double_t x=mPt;
                                        weight=CalculateWeight(111, x);
                                        
                                        
                                    }
                                    if(TMath::Abs(fMCparticleMother->GetPdgCode())==221){
                                        Double_t x=mPt;
                                        weight=CalculateWeight(221, x);
                                        
                                        
                                    }
                                    
                                    //weight for grandmother
                                    
                                    if((fMCparticleMother->GetMother()>0) && TMath::Abs(((fMCparticleGMother->GetPdgCode())==111))){
                                        Double_t gmPt=fMCparticleGMother->Pt();
                                        Double_t x=gmPt;
                                        weight=CalculateWeight(111, x);
                                        
                                    }
                                    if((fMCparticleMother->GetMother()>0) && TMath::Abs(((fMCparticleGMother->GetPdgCode())==221))){
                                        Double_t gmPt=fMCparticleGMother->Pt();
                                        Double_t x=gmPt;
                                        weight=CalculateWeight(221, x);
                                        
                                    }
                                    
                                    
                                    if(fMCparticle2->GetMother()==fMCparticle->GetMother()) fPtElec_ULS_MC_weight->Fill(fPtE, (fNonHFE->GetNULS())*1./weight);
                                    
                                    //-----------------------------------------------------------------------------------------------------------
                                    //end of weight
                                    
                                }//partner found same as track
                            }//loop in all partner
                            
                        }//track
                    }//is ULS
                    //----------------------------------------------------------------------------
                    //end of part to check other way to calculate efficiency
                    //----------------------------------------------------------------------------
                    
                }//IsTPConly
                
            }//particle kind
            
            if(IsTPConly){
                //ULS-LS with no pid AOD
                if(fNonHFE->IsULS()) fPtElec_ULS_NoPid->Fill(fPtE,fNonHFE->GetNULS());
                if(fNonHFE->IsLS()) fPtElec_LS_NoPid->Fill(fPtE,fNonHFE->GetNLS());
            }
            
        }//close IsAOD
        //It is ESD
        else
        {
            if(TMath::Abs(fMCtrack->GetPdgCode())==11 && (TMath::Abs(fMCtrackMother->GetPdgCode())==22 || TMath::Abs(fMCtrackMother->GetPdgCode())==111 || TMath::Abs(fMCtrackMother->GetPdgCode())==221))
            {
                if(!IsTPConly){
                    if(fNonHFE->IsULS()) fPtElec_ULS->Fill(fPtE,fNonHFE->GetNULS());
                    if(fNonHFE->IsLS()) fPtElec_LS->Fill(fPtE,fNonHFE->GetNLS());
                }
                
                if(IsTPConly){
                    if(fNonHFE->IsULS()) fPtElec_ULS2->Fill(fPtE,fNonHFE->GetNULS());
                    if(fNonHFE->IsLS()) fPtElec_LS2->Fill(fPtE,fNonHFE->GetNLS());
                }
            }
            
            
        }
    }//close IsMC
    ///_________________________________________________________________
    //not MC
    else
    {
        if(!IsTPConly){
            if(fNonHFE->IsULS()) fPtElec_ULS->Fill(fPtE,fNonHFE->GetNULS());
            if(fNonHFE->IsLS()) fPtElec_LS->Fill(fPtE,fNonHFE->GetNLS());
        }
        
        if(IsTPConly){
            if(fNonHFE->IsULS()) fPtElec_ULS2->Fill(fPtE,fNonHFE->GetNULS());
            if(fNonHFE->IsLS()) fPtElec_LS2->Fill(fPtE,fNonHFE->GetNLS());
        }
    }
    
    //for MC closure test
    // to be used as data, with MC input
    if(!IsTPConly){
        if(fNonHFE->IsULS()) fPtElec_ULS_mc_closure->Fill(fPtE,fNonHFE->GetNULS());
        if(fNonHFE->IsLS()) fPtElec_LS_mc_closure->Fill(fPtE,fNonHFE->GetNLS());
    }
    
    if(IsTPConly){
        if(fNonHFE->IsULS()) fPtElec_ULS2_mc_closure->Fill(fPtE,fNonHFE->GetNULS());
        if(fNonHFE->IsLS()) fPtElec_LS2_mc_closure->Fill(fPtE,fNonHFE->GetNLS());
    }
    
    
    
}//end of fPtElec_ULS function

//______________________________________________________________________
void AliAnalysisTaskHFEpACorrelation::ElectronHadronCorrelation(AliVTrack *track, Int_t trackIndex, AliVParticle *vtrack)
{
    
    ///_________________________________________________________________
    ///MC analysis
    
    fIsHFE1 = kFALSE;
    
    if(fIsMC)
    {
        if(track->GetLabel() < 0)
        {
            AliWarning(Form("The track %d does not have a valid MC label",trackIndex));
            //return;
        }
        else
        {
            if(fIsAOD)
            {
                fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
                
                if(fMCparticle->GetMother()<0)
                    AliWarning(Form("The track %d does not have a valid Mother MC label",trackIndex));
                else
                {
                    fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                    if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==22 || TMath::Abs(fMCparticleMother->GetPdgCode())==111 || TMath::Abs(fMCparticleMother->GetPdgCode())==221))
                    {
                        //Is Background
                        fPtBackgroundBeforeReco->Fill(track->Pt());
                        //new 17-09
                        for(int hadron_index = 0; hadron_index < fMCarray->GetEntries(); hadron_index++)
                        {
                            if (hadron_index == track->GetLabel() ) continue; //self
                            
                            AliAODMCParticle *hadron_MC = (AliAODMCParticle*) fMCarray->At(hadron_index);
                            
                            if(!hadron_MC->IsPrimary()) continue;
                            
                            if (hadron_MC->Eta()<fEtaCutMin || hadron_MC->Eta()>fEtaCutMax) continue;
                            
                            if (hadron_MC->Charge() == 0) continue;
                            
                            if(hadron_MC->Pt()<fAssHadronPtMin || hadron_MC->Pt()>fAssHadronPtMax) continue;
                            
                            float delta_phi_MC;
                            float delta_eta_MC;
                            
                            double pi =  TMath::Pi();
                            
                            delta_phi_MC = fMCparticle->Phi() - hadron_MC->Phi();
                            
                            delta_eta_MC = fMCparticle->Eta() - hadron_MC->Eta();
                            
                            if (delta_phi_MC> 3*pi/2) delta_phi_MC = delta_phi_MC - 2*pi;
                            if (delta_phi_MC < -pi/2)  delta_phi_MC= delta_phi_MC + 2*pi;
                            
                            
                            Double_t fPtBin[7];
                            if(fUseAlternativeBinnig)
                            {
                                fPtBin[0] = 0.5;
                                fPtBin[1] = 1.0;
                                fPtBin[2] = 1.5;
                                fPtBin[3] = 2.0;
                                fPtBin[4] = 3.0;
                                fPtBin[5] = 4.0;
                                fPtBin[6] = 6.0;
                            }
                            else
                            {
                                fPtBin[0] = 1.0;
                                fPtBin[1] = 2.0;
                                fPtBin[2] = 4.0;
                                fPtBin[3] = 6.0;
                                fPtBin[4] = 8.0;
                                fPtBin[5] = 10.0;
                                fPtBin[6] = 15.0;
                            }
                            
                            
                            
                            for(Int_t i = 0; i < 6; i++)
                            {
                                if(fMCparticle->Pt()>=fPtBin[i] && fMCparticle->Pt()<fPtBin[i+1])
                                {
                                    fCetaPhi_MC_NHFE_1partner_reco[i]->Fill(delta_phi_MC,delta_eta_MC);
                                }
                            }
                            
                        }//end new 17-09
                    }
                }
                
                Bool_t FindMotherHFE = kFALSE;
                
                FindMotherHFE = FindMother(track->GetLabel());
                if (fIsHFE1)
                {
                    fpT_MC_HFE_RECO_MC->Fill(fMCparticle->Pt());
                    fpTReco_vs_MC->Fill(track->Pt(),fMCparticle->Pt());
                    
                    for(int hadron_index = 0; hadron_index < fMCarray->GetEntries(); hadron_index++)
                    {
                        if (hadron_index == track->GetLabel() ) continue;//self
                        
                        AliAODMCParticle *hadron_MC = (AliAODMCParticle*) fMCarray->At(hadron_index);
                        
                        if (hadron_MC->Eta()<fEtaCutMin || hadron_MC->Eta()>fEtaCutMax) continue;
                        
                        if (hadron_MC->Charge() == 0) continue;
                        
                        if(hadron_MC->Pt()<fAssHadronPtMin || hadron_MC->Pt()>fAssHadronPtMax) continue;
                        
                        float delta_phi_MC;
                        float delta_eta_MC;
                        double pi =  TMath::Pi();
                        delta_phi_MC = fMCparticle->Phi() - hadron_MC->Phi();
                        delta_eta_MC = fMCparticle->Eta() - hadron_MC->Eta();
                        
                        if (delta_phi_MC> 3*pi/2) delta_phi_MC = delta_phi_MC - 2*pi;
                        if (delta_phi_MC < -pi/2)  delta_phi_MC= delta_phi_MC + 2*pi;
                        
                        
                        Double_t fPtBin[7];
                        if(fUseAlternativeBinnig)
                        {
                            fPtBin[0] = 0.5;
                            fPtBin[1] = 1.0;
                            fPtBin[2] = 1.5;
                            fPtBin[3] = 2.0;
                            fPtBin[4] = 3.0;
                            fPtBin[5] = 4.0;
                            fPtBin[6] = 6.0;
                        }
                        else
                        {
                            fPtBin[0] = 1.0;
                            fPtBin[1] = 2.0;
                            fPtBin[2] = 4.0;
                            fPtBin[3] = 6.0;
                            fPtBin[4] = 8.0;
                            fPtBin[5] = 10.0;
                            fPtBin[6] = 15.0;
                        }
                        
                        
                        
                        for(Int_t i = 0; i < 6; i++)
                        {
                            if(fMCparticle->Pt()>=fPtBin[i] && fMCparticle->Pt()<fPtBin[i+1])
                            {
                                if(hadron_MC->IsPhysicalPrimary())
                                    fCetaPhi_MC_HFE_RECO_MC_PhyPrimH[i]->Fill(delta_phi_MC,delta_eta_MC);
                            }
                        }
                        
                        for(Int_t i = 0; i < 6; i++)
                        {
                            if(track->Pt()>=fPtBin[i] && track->Pt()<fPtBin[i+1])
                            {
                                if(hadron_MC->IsPhysicalPrimary())
                                    fCetaPhi_MC_HFE_pTofReco[i]->Fill(delta_phi_MC,delta_eta_MC);
                            }
                        }
                    }
                    
                    
                    
                    
                    //now I correlate the Reco HFE with the MC phys primary hadrons
                    
                    
                    for(int hadron_index = 0; hadron_index < fMCarray->GetEntries(); hadron_index++)
                    {
                        if (hadron_index == track->GetLabel() ) continue;//self
                        
                        AliAODMCParticle *hadron_MC = (AliAODMCParticle*) fMCarray->At(hadron_index);
                        
                        if (hadron_MC->Eta()<fEtaCutMin || hadron_MC->Eta()>fEtaCutMax) continue;
                        
                        if (hadron_MC->Charge() == 0) continue;
                        
                        if(hadron_MC->Pt()<fAssHadronPtMin || hadron_MC->Pt()>fAssHadronPtMax) continue;
                        
                        float delta_phi_MC;
                        float delta_eta_MC;
                        double pi =  TMath::Pi();
                        delta_phi_MC = track->Phi() - hadron_MC->Phi();
                        delta_eta_MC = track->Eta() - hadron_MC->Eta();
                        
                        if (delta_phi_MC> 3*pi/2) delta_phi_MC = delta_phi_MC - 2*pi;
                        if (delta_phi_MC < -pi/2)  delta_phi_MC= delta_phi_MC + 2*pi;
                        
                        
                        Double_t fPtBin[7];
                        if(fUseAlternativeBinnig)
                        {
                            fPtBin[0] = 0.5;
                            fPtBin[1] = 1.0;
                            fPtBin[2] = 1.5;
                            fPtBin[3] = 2.0;
                            fPtBin[4] = 3.0;
                            fPtBin[5] = 4.0;
                            fPtBin[6] = 6.0;
                        }
                        else
                        {
                            fPtBin[0] = 1.0;
                            fPtBin[1] = 2.0;
                            fPtBin[2] = 4.0;
                            fPtBin[3] = 6.0;
                            fPtBin[4] = 8.0;
                            fPtBin[5] = 10.0;
                            fPtBin[6] = 15.0;
                        }
                        
                        
                        
                        for(Int_t i = 0; i < 6; i++)
                        {
                            if(track->Pt()>=fPtBin[i] && track->Pt()<fPtBin[i+1])
                            {
                                if(hadron_MC->IsPhysicalPrimary())
                                    fCetaPhi_Reco_with_MChadrons[i]->Fill(delta_phi_MC,delta_eta_MC);
                            }
                        }
                    }
                    
                    
                }
                
            }
            else
            {
                fMCtrack = fMCstack->Particle(track->GetLabel());
                
                if(fMCtrack->GetFirstMother()<0) return;
                
                fMCtrackMother = fMCstack->Particle(fMCtrack->GetFirstMother());
                
                if(TMath::Abs(fMCtrack->GetPdgCode())==11 && (TMath::Abs(fMCtrackMother->GetPdgCode())==22 || TMath::Abs(fMCtrackMother->GetPdgCode())==111 || TMath::Abs(fMCtrackMother->GetPdgCode())==221))
                {
                    //Is Background
                    fPtBackgroundBeforeReco->Fill(track->Pt());
                }
            }
        }
    }
    ///_________________________________________________________________
    
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
    
    ///#################################################################
    //Non-HFE reconstruction
    //fNonHFE = new AliSelectNonHFE(); I do not need to create if again (almost sure)
    fNonHFE->SetAODanalysis(fIsAOD);
    if(fMassCutFlag) fNonHFE->SetInvariantMassCut(fMassCut);
    if(fAngleCutFlag) fNonHFE->SetOpeningAngleCut(fAngleCut);
    if(fChi2CutFlag) fNonHFE->SetChi2OverNDFCut(fChi2Cut);
    if(fDCAcutFlag) fNonHFE->SetDCACut(fDCAcut);
    fNonHFE->SetAlgorithm("DCA"); //KF
    fNonHFE->SetPIDresponse(fPidResponse);
    fNonHFE->SetTrackCuts(-3.5,3.5,fPartnerCuts);
    fNonHFE->SetAdditionalCuts(fPtMinAsso,fTpcNclsAsso);
    
    
    fNonHFE->SetHistAngleBack(fOpAngleBack);
    fNonHFE->SetHistAngle(fOpAngle);
    fNonHFE->SetHistDCABack(fDCABack);
    fNonHFE->SetHistDCA(fDCA);
    fNonHFE->SetHistMassBack(fInvMassBack);
    fNonHFE->SetHistMass(fInvMass);
    
    fNonHFE->FindNonHFE(trackIndex,vtrack,fVevent);
    
    Int_t *fUlsPartner = fNonHFE->GetPartnersULS();
    Int_t *fLsPartner = fNonHFE->GetPartnersLS();
    Bool_t fUlsIsPartner = kFALSE;
    Bool_t fLsIsPartner = kFALSE;
    ///#################################################################
    
    
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
    ///MC analysis
    
    
    ///--------------------------------------------
    //--------new MC analysis AOD in 18-09---------
    
    if(fIsMC && fIsAOD)
    {
        //is and NHFE?
        if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==22 || TMath::Abs(fMCparticleMother->GetPdgCode())==111 || TMath::Abs(fMCparticleMother->GetPdgCode())==221))
        {
            //Do they have the same mother?
            Bool_t IsFromTheSameMother = kFALSE;
            
            int mother_electron_ID = fMCparticleMother->GetLabel();
            
            for(Int_t i = 0; i < fNonHFE->GetNULS(); i++)
            {
                
                AliVParticle* Vtrack_uls_partner = fVevent->GetTrack(fUlsPartner[i]);
                if (!Vtrack_uls_partner)
                {
                    printf("ERROR: Could not receive track %d\n", fUlsPartner[i]);
                    continue;
                }
                AliVTrack *track_partner_ULS = dynamic_cast<AliVTrack*>(Vtrack_uls_partner);
                if(track_partner_ULS->GetLabel()<0)
                {
                    //printf(" \n ULS partner dos not have a valid MC label! \n ");
                    continue;
                }
                
                AliAODMCParticle* fMCparticle_uls_partner = (AliAODMCParticle*) fMCarray->At(track_partner_ULS->GetLabel());
                
                if( mother_electron_ID == fMCparticle_uls_partner->GetMother() )
                {
                    IsFromTheSameMother = kTRUE;
                }
                
            }
            
            if(IsFromTheSameMother)
                fpT_MC_with_partner_greater->Fill(fMCparticle->Pt());
            else
                fpT_MC_with_partner_below->Fill(fMCparticle->Pt());
            
            //correlate with hadrons that are primary
            for(int hadron_index = 0; hadron_index < fMCarray->GetEntries(); hadron_index++)
            {
                if (hadron_index == track->GetLabel() ) continue;//self
                
                AliAODMCParticle *hadron_MC = (AliAODMCParticle*) fMCarray->At(hadron_index);
                
                if(!hadron_MC->IsPrimary()) continue;
                
                if (hadron_MC->Eta()<fEtaCutMin || hadron_MC->Eta()>fEtaCutMax) continue;
                
                if (hadron_MC->Charge() == 0) continue;
                
                if(hadron_MC->Pt()<fAssHadronPtMin || hadron_MC->Pt()>fAssHadronPtMax) continue;
                
                float delta_phi_MC;
                float delta_eta_MC;
                
                double pi =  TMath::Pi();
                
                delta_phi_MC = fMCparticle->Phi() - hadron_MC->Phi();
                
                delta_eta_MC = fMCparticle->Eta() - hadron_MC->Eta();
                
                if (delta_phi_MC> 3*pi/2) delta_phi_MC = delta_phi_MC - 2*pi;
                if (delta_phi_MC < -pi/2)  delta_phi_MC= delta_phi_MC + 2*pi;
                
                
                Double_t fPtBin[7];
                if(fUseAlternativeBinnig)
                {
                    fPtBin[0] = 0.5;
                    fPtBin[1] = 1.0;
                    fPtBin[2] = 1.5;
                    fPtBin[3] = 2.0;
                    fPtBin[4] = 3.0;
                    fPtBin[5] = 4.0;
                    fPtBin[6] = 6.0;
                }
                else
                {
                    fPtBin[0] = 1.0;
                    fPtBin[1] = 2.0;
                    fPtBin[2] = 4.0;
                    fPtBin[3] = 6.0;
                    fPtBin[4] = 8.0;
                    fPtBin[5] = 10.0;
                    fPtBin[6] = 15.0;
                }
                
                
                
                for(Int_t i = 0; i < 6; i++)
                {
                    if(fMCparticle->Pt()>=fPtBin[i] && fMCparticle->Pt()<fPtBin[i+1])
                    {
                        //If they are from the same mother, this means that the partner has a pT greater than the minimum
                        // and had fulfill all the requeriments and is recontructed. This distribution contain only the case with two recontructed
                        if (IsFromTheSameMother)
                            fCetaPhi_MC_with_partner_greater[i]->Fill(delta_phi_MC,delta_eta_MC);
                        else
                            fCetaPhi_MC_with_partner_below[i]->Fill(delta_phi_MC,delta_eta_MC);
                        
                        //If no partner has the same mother as the trigger, than the
                    }
                }
                
            }
            
            
        } //end is NHFE
    }
    
    
    
    
    
    ///-----------end of new MC analysis-----------
    ///--------------------------------------------
    
    
    
    
    
    if(fIsMC)
    {
        if(fIsAOD)
        {
            if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==22 || TMath::Abs(fMCparticleMother->GetPdgCode())==111 || TMath::Abs(fMCparticleMother->GetPdgCode())==221))
            {
                if(fNonHFE->IsULS()) fPtElec_ULS->Fill(fPtE,fNonHFE->GetNULS());
                if(fNonHFE->IsLS()) fPtElec_LS->Fill(fPtE,fNonHFE->GetNLS());
                
                
            }
            
            if(fNonHFE->IsULS()) fPtElec_ULS_NoPid->Fill(fPtE,fNonHFE->GetNULS());
            if(fNonHFE->IsLS()) fPtElec_LS_NoPid->Fill(fPtE,fNonHFE->GetNLS());
        }
        else
        {
            if(TMath::Abs(fMCtrack->GetPdgCode())==11 && (TMath::Abs(fMCtrackMother->GetPdgCode())==22 || TMath::Abs(fMCtrackMother->GetPdgCode())==111 || TMath::Abs(fMCtrackMother->GetPdgCode())==221))
            {
                if(fNonHFE->IsULS()) fPtElec_ULS->Fill(fPtE,fNonHFE->GetNULS());
                if(fNonHFE->IsLS()) fPtElec_LS->Fill(fPtE,fNonHFE->GetNLS());
                
            }
            
            
        }
    }
    ///_________________________________________________________________
    else
    {
        if(fNonHFE->IsULS()) fPtElec_ULS->Fill(fPtE,fNonHFE->GetNULS());
        if(fNonHFE->IsLS()) fPtElec_LS->Fill(fPtE,fNonHFE->GetNLS());
    }
    
    
    
    
    //__________________________________________________________________
    //Event Mixing Analysis - Hadron Loop
    //Retrieve
    
    if (fIsMC)
    {
        if (fIsHFE1)
        {
            fpT_Data_HFE_RECO_Data->Fill(fPtE);
        }
    }
    
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
            fPool = fPoolMgr->GetEventPool(fCentrality->GetCentralityPercentile("V0A"), fZvtx); // Get the buffer associated with the current centrality and z-vtx
            
            if(!fPool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f",fCentrality->GetCentralityPercentile("V0A"), fZvtx));
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
                    
                    Double_t fPtBin[7];
                    if(fUseAlternativeBinnig)
                    {
                        fPtBin[0] = 0.5;
                        fPtBin[1] = 1.0;
                        fPtBin[2] = 1.5;
                        fPtBin[3] = 2.0;
                        fPtBin[4] = 3.0;
                        fPtBin[5] = 4.0;
                        fPtBin[6] = 6.0;
                    }
                    else
                    {
                        fPtBin[0] = 1.0;
                        fPtBin[1] = 2.0;
                        fPtBin[2] = 4.0;
                        fPtBin[3] = 6.0;
                        fPtBin[4] = 8.0;
                        fPtBin[5] = 10.0;
                        fPtBin[6] = 15.0;
                    }
                    
                    
                    
                    for(Int_t i = 0; i < 6; i++)
                    {
                        if(fPtE>=fPtBin[i] && fPtE<fPtBin[i+1])
                        {
                            fCEtaPhi_Inc_EM[i]->Fill(fDphi,fDeta);
                            
                            if(fNonHFE->IsULS()) fCEtaPhi_ULS_EM[i]->Fill(fDphi,fDeta);
                            if(fNonHFE->IsLS()) fCEtaPhi_LS_EM[i]->Fill(fDphi,fDeta);
                            
                            if(fNonHFE->IsULS()) fCEtaPhi_ULS_Weight_EM[i]->Fill(fDphi,fDeta,fNonHFE->GetNULS());
                            if(fNonHFE->IsLS()) fCEtaPhi_LS_Weight_EM[i]->Fill(fDphi,fDeta,fNonHFE->GetNLS());
                            
                            
                        }
                    }
                    
                    
                    
                    
                    // TODO your code: do event mixing with current event and bgTracks
                    // note that usually the content filled now is weighted by 1 / pool->GetCurrentNEvents()
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
        
        if(track2->Eta()<fEtaCutMin || track2->Eta()>fEtaCutMax ) continue;
        
        if(fIsAOD)
        {
            AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(Vtrack2);
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
        
        
        
        //Systematics: test of seconday particle contamination
        Bool_t IsPhyPrimaryHadron;
        
        if(fIsMC)
        {
            if(fIsAOD)
            {
                AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(Vtrack2);
                if (track2->GetLabel()<0)
                    IsPhyPrimaryHadron = kFALSE;
                else
                {
                    AliAODMCParticle *hadron_MonteCarlo = (AliAODMCParticle*) fMCarray->At(track2->GetLabel());
                    if (hadron_MonteCarlo->IsPhysicalPrimary())
                        IsPhyPrimaryHadron = kTRUE;
                    else
                        IsPhyPrimaryHadron = kFALSE;
                }
            }
        }
        
        fPhiH = track2->Phi();
        fEtaH = track2->Eta();
        fPtH = track2->Pt();
        
        if(fPtH<fAssHadronPtMin || fPtH>fAssHadronPtMax) continue;
        
        fDphi = fPhiE - fPhiH;
        
        if (fDphi > 3*pi/2) fDphi = fDphi - 2*pi;
        if (fDphi < -pi/2)  fDphi = fDphi + 2*pi;
        
        fDeta = fEtaE - fEtaH;
        
        Double_t fPtBin[7];
        if(fUseAlternativeBinnig)
        {
            fPtBin[0] = 0.5;
            fPtBin[1] = 1.0;
            fPtBin[2] = 1.5;
            fPtBin[3] = 2.0;
            fPtBin[4] = 3.0;
            fPtBin[5] = 4.0;
            fPtBin[6] = 6.0;
        }
        else
        {
            fPtBin[0] = 1.0;
            fPtBin[1] = 2.0;
            fPtBin[2] = 4.0;
            fPtBin[3] = 6.0;
            fPtBin[4] = 8.0;
            fPtBin[5] = 10.0;
            fPtBin[6] = 15.0;
        }
        
        
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
        
        for(Int_t i = 0; i < 6; i++)
        {
            if(fPtE>=fPtBin[i] && fPtE<fPtBin[i+1])
            {
                Double_t HadronWeight = GetHadronWCorrection(fPtH);
                //Filling histograms
                fCEtaPhi_Inc[i]->Fill(fDphi,fDeta,HadronWeight);
                if(fNonHFE->IsULS()) fCEtaPhi_ULS[i]->Fill(fDphi,fDeta,HadronWeight);
                if(fNonHFE->IsLS()) fCEtaPhi_LS[i]->Fill(fDphi,fDeta,HadronWeight);
                if(fNonHFE->IsULS() && !fUlsIsPartner && !fLsIsPartner) fCEtaPhi_ULS_NoP[i]->Fill(fDphi,fDeta,HadronWeight);
                if(fNonHFE->IsLS() && !fUlsIsPartner && !fLsIsPartner) fCEtaPhi_LS_NoP[i]->Fill(fDphi,fDeta,HadronWeight);
                if(fNonHFE->IsULS()) fCEtaPhi_ULS_Weight[i]->Fill(fDphi,fDeta,fNonHFE->GetNULS()*HadronWeight);
                if(fNonHFE->IsLS()) fCEtaPhi_LS_Weight[i]->Fill(fDphi,fDeta,fNonHFE->GetNLS()*HadronWeight);
                if(fNonHFE->IsULS() && !fUlsIsPartner && !fLsIsPartner) fCEtaPhi_ULS_NoP_Weight[i]->Fill(fDphi,fDeta,fNonHFE->GetNULS()*HadronWeight);
                if(fNonHFE->IsLS() && !fUlsIsPartner && !fLsIsPartner) fCEtaPhi_LS_NoP_Weight[i]->Fill(fDphi,fDeta,fNonHFE->GetNLS()*HadronWeight);
                
                if(fIsMC)
                {
                    if(IsPhyPrimaryHadron)
                    {
                        fCEtaPhi_Inc_PH[i]->Fill(fDphi,fDeta);
                        if(fNonHFE->IsULS()) fCEtaPhi_ULS_Weight_PH[i]->Fill(fDphi,fDeta,fNonHFE->GetNULS());
                        if(fNonHFE->IsLS()) fCEtaPhi_LS_Weight_PH[i]->Fill(fDphi,fDeta,fNonHFE->GetNLS());
                        if(fNonHFE->IsULS() && !fUlsIsPartner && !fLsIsPartner) fCEtaPhi_ULS_NoP_Weight_PH[i]->Fill(fDphi,fDeta,fNonHFE->GetNULS());
                        if(fNonHFE->IsLS() && !fUlsIsPartner && !fLsIsPartner) fCEtaPhi_LS_NoP_Weight_PH[i]->Fill(fDphi,fDeta,fNonHFE->GetNLS());
                        
                    }
                    
                }
                
                if (fIsMC)
                {
                    if (fIsHFE1)
                    {
                        fCetaPhi_Data_Data_RECO_MC_PhyPrimH[i]->Fill(fDphi,fDeta);
                        if(IsPhyPrimaryHadron)
                            fCEtaPhi_DataHFE_with_onlyPhysPriHadron[i]->Fill(fDphi,fDeta);
                        
                        
                    }
                    
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
        
        if(track2->Eta()<fEtaCutMin || track2->Eta()>fEtaCutMax ) continue;
        
        if(fIsAOD)
        {
            AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(Vtrack2);
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
        
        fTracksClone->Add(new AliHFEHCParticle(track2->Eta(), track2->Phi(), track2->Pt()));
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
        if(track2->Eta()<fEtaCutMin || track2->Eta()>fEtaCutMax ) continue;
        
        if(fIsAOD)
        {
            AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(Vtrack2);
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
        
        
        fPhiH = track2->Phi();
        fEtaH = track2->Eta();
        fPtH = track2->Pt();
        
        if(fPtH<fAssHadronPtMin || fPtH>fAssHadronPtMax) continue;
        
        fDphi = fPhiE - fPhiH;
        
        if (fDphi > 3*pi/2) fDphi = fDphi - 2*pi;
        if (fDphi < -pi/2)  fDphi = fDphi + 2*pi;
        
        fDeta = fEtaE - fEtaH;
        
        Double_t fPtBin[7];
        if(fUseAlternativeBinnig)
        {
            fPtBin[0] = 0.5;
            fPtBin[1] = 1.0;
            fPtBin[2] = 1.5;
            fPtBin[3] = 2.0;
            fPtBin[4] = 3.0;
            fPtBin[5] = 4.0;
            fPtBin[6] = 6.0;
        }
        else
        {
            fPtBin[0] = 1.0;
            fPtBin[1] = 2.0;
            fPtBin[2] = 4.0;
            fPtBin[3] = 6.0;
            fPtBin[4] = 8.0;
            fPtBin[5] = 10.0;
            fPtBin[6] = 15.0;
        }
        
        
        
        for(Int_t i = 0; i < 6; i++)
        {
            if(fPtE>=fPtBin[i] && fPtE<fPtBin[i+1])
            {
                fCEtaPhi_Inc_DiHadron[i]->Fill(fDphi,fDeta);
            }
        }
    }
}
//____________________________________________________________________________________________________________

//______________________________________________________________________
Bool_t AliAnalysisTaskHFEpACorrelation::FindMother(Int_t mcIndex)
{
    fIsHFE1 = kFALSE;
    fIsHFE2 = kFALSE;
    fIsNonHFE = kFALSE;
    fIsFromD = kFALSE;
    fIsFromB = kFALSE;
    fIsFromPi0 = kFALSE;
    fIsFromEta = kFALSE;
    fIsFromGamma = kFALSE;
    
    if(mcIndex < 0 || !fIsMC)
    {
        return kFALSE;
    }
    
    Int_t pdg = -99999;
    Int_t mpdg = -99999;
    Int_t gmpdg = -99999;
    Int_t ggmpdg = -99999;
    Int_t gggmpdg = -99999;
    
    if(fIsAOD)
    {
        fMCparticle = (AliAODMCParticle*) fMCarray->At(mcIndex);
        
        pdg = TMath::Abs(fMCparticle->GetPdgCode());
        
        
        if(pdg!=11)
        {
            fIsHFE1 = kFALSE;
            fIsHFE2 = kFALSE;
            fIsNonHFE = kFALSE;
            fIsFromD = kFALSE;
            fIsFromB = kFALSE;
            fIsFromPi0 = kFALSE;
            fIsFromEta = kFALSE;
            fIsFromGamma = kFALSE;
            return kFALSE;
        }
        
        if(fMCparticle->GetMother()<0)
        {
            fIsHFE1 = kFALSE;
            fIsHFE2 = kFALSE;
            fIsNonHFE = kFALSE;
            fIsFromD = kFALSE;
            fIsFromB = kFALSE;
            fIsFromPi0 = kFALSE;
            fIsFromEta = kFALSE;
            fIsFromGamma = kFALSE;
            return kFALSE;
        }
        
        fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
        mpdg = TMath::Abs(fMCparticleMother->GetPdgCode());
        
        if(fMCparticleMother->GetMother()<0)
        {
            gmpdg = 0;
            ggmpdg = 0;
            gggmpdg = 0;
        }
        else
        {
            fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
            gmpdg = TMath::Abs(fMCparticleGMother->GetPdgCode());
            if(fMCparticleGMother->GetMother()<0)
            {
                ggmpdg = 0;
                gggmpdg = 0;
            }
            else
            {
                fMCparticleGGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleGMother->GetMother());
                ggmpdg = TMath::Abs(fMCparticleGGMother->GetPdgCode());
                if(fMCparticleGGMother->GetMother()<0)
                {
                    gggmpdg = 0;
                }
                else
                {
                    fMCparticleGGGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleGGMother->GetMother());
                    gggmpdg = TMath::Abs(fMCparticleGGGMother->GetPdgCode());
                }
            }
        }
    }
    else
    {
        fMCtrack = fMCstack->Particle(mcIndex);
        
        pdg = TMath::Abs(fMCtrack->GetPdgCode());
        
        if(pdg!=11)
        {
            fIsHFE1 = kFALSE;
            fIsHFE2 = kFALSE;
            fIsNonHFE = kFALSE;
            fIsFromD = kFALSE;
            fIsFromB = kFALSE;
            fIsFromPi0 = kFALSE;
            fIsFromEta = kFALSE;
            fIsFromGamma = kFALSE;
            return kFALSE;
        }
        
        if(fMCtrack->GetFirstMother()<0)
        {
            fIsHFE1 = kFALSE;
            fIsHFE2 = kFALSE;
            fIsNonHFE = kFALSE;
            fIsFromD = kFALSE;
            fIsFromB = kFALSE;
            fIsFromPi0 = kFALSE;
            fIsFromEta = kFALSE;
            fIsFromGamma = kFALSE;
            return kFALSE;
        }
        
        fMCtrackMother = fMCstack->Particle(fMCtrack->GetFirstMother());
        mpdg = TMath::Abs(fMCtrackMother->GetPdgCode());
        
        if(fMCtrackMother->GetFirstMother()<0)
        {
            gmpdg = 0;
            ggmpdg = 0;
            gggmpdg = 0;
        }
        else
        {
            fMCtrackGMother = fMCstack->Particle(fMCtrackMother->GetFirstMother());
            gmpdg = TMath::Abs(fMCtrackGMother->GetPdgCode());
            
            if(fMCtrackGMother->GetFirstMother()<0)
            {
                ggmpdg = 0;
                gggmpdg = 0;
            }
            else
            {
                fMCtrackGGMother = fMCstack->Particle(fMCtrackGMother->GetFirstMother());
                ggmpdg = TMath::Abs(fMCtrackGGMother->GetPdgCode());
                
                if(fMCtrackGGMother->GetFirstMother()<0)
                {
                    gggmpdg = 0;
                }
                else
                {
                    fMCtrackGGGMother = fMCstack->Particle(fMCtrackGGMother->GetFirstMother());
                    gggmpdg = TMath::Abs(fMCtrackGGGMother->GetPdgCode());
                }
            }
        }
    }
    
    //Tag Electron Source
    if(mpdg==111 || mpdg==221 || mpdg==22)
    {
        fIsHFE1 = kFALSE;
        fIsHFE2 = kFALSE;
        fIsNonHFE = kTRUE;
        fIsFromD = kFALSE;
        fIsFromB = kFALSE;
        
        fIsFromPi0 = kFALSE;
        fIsFromEta = kFALSE;
        fIsFromGamma = kFALSE;
        
        if(mpdg==111) fIsFromPi0 = kFALSE;
        if(mpdg==221)fIsFromEta = kFALSE;
        if(mpdg==22) fIsFromGamma = kFALSE;
        
        return kTRUE;
    }
    else
    {
        fIsHFE1 = kFALSE;
        fIsHFE2 = kTRUE;
        
        fIsFromPi0 = kFALSE;
        fIsFromEta = kFALSE;
        fIsFromGamma = kFALSE;
        
        fIsNonHFE = kFALSE;
        
        fIsFromD = kFALSE;
        fIsFromB = kFALSE;
        
        if(mpdg>400 && mpdg<500)
        {
            if((gmpdg>500 && gmpdg<600) || (ggmpdg>500 && ggmpdg<600) || (gggmpdg>500 && gggmpdg<600))
            {
                fIsHFE1 = kTRUE;
                fIsFromD = kFALSE;
                fIsFromB = kTRUE;
                return kTRUE;
            }
            else
            {
                fIsHFE1 = kTRUE;
                fIsFromD = kTRUE;
                fIsFromB = kFALSE;
                return kTRUE;
            }
        }
        else if(mpdg>500 && mpdg<600)
        {
            fIsHFE1 = kTRUE;
            fIsFromD = kFALSE;
            fIsFromB = kTRUE;
            return kTRUE;
        }
        else
        {
            fIsHFE1 = kFALSE;
            fIsFromD = kFALSE;
            fIsFromB = kFALSE;
            return kFALSE;
        }
    }
}

Double_t AliAnalysisTaskHFEpACorrelation::CalculateWeight(Int_t pdg_particle, Double_t x)
{
    //weight for d3 based on MinJung parametrization //sent by Jan
    Double_t weight=1;
    //MB
    if(!fUseTrigger){
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
        else weight=1;
        
    }//not trigger
    
    //trigger case LHC14b3a,b,c
    if(fUseTrigger){
        
        
        if(pdg_particle==111){
            if(x>= 0.000000 &&  x < 0.150000 ) weight=5.413029271505721773394270712743;
            if(x>= 0.150000 &&  x < 0.300000 ) weight=1.843465358197068582057909225114;
            if(x>= 0.300000 &&  x < 0.450000 ) weight=1.562291242835710747272059961688;
            if(x>= 0.450000 &&  x < 0.600000 ) weight=1.390412237790664473635615649982;
            if(x>= 0.600000 &&  x < 0.750000 ) weight=1.246791913502574455918647799990;
            if(x>= 0.750000 &&  x < 0.900000 ) weight=1.146357892193592631358001199260;
            if(x>= 0.900000 &&  x < 1.050000 ) weight=1.099457735840871253074624291912;
            if(x>= 1.050000 &&  x < 1.200000 ) weight=1.108489584943816996798204854713;
            if(x>= 1.200000 &&  x < 1.350000 ) weight=1.168142381612309543825745095091;
            if(x>= 1.350000 &&  x < 1.500000 ) weight=1.279815553925982118954607358319;
            if(x>= 1.500000 &&  x < 1.650000 ) weight=1.447372173802167649725447517994;
            if(x>= 1.650000 &&  x < 1.800000 ) weight=1.675200256094552031527200597338;
            if(x>= 1.800000 &&  x < 1.950000 ) weight=1.965825287403190735346925066551;
            if(x>= 1.950000 &&  x < 2.100000 ) weight=2.334029600768103840380263136467;
            if(x>= 2.100000 &&  x < 2.250000 ) weight=2.812404460026804553507417949731;
            if(x>= 2.250000 &&  x < 2.400000 ) weight=3.399170813174843708281969156815;
            if(x>= 2.400000 &&  x < 2.550000 ) weight=4.113418509691398661232142330846;
            if(x>= 2.550000 &&  x < 2.700000 ) weight=5.028714136553443125876583508216;
            if(x>= 2.700000 &&  x < 2.850000 ) weight=6.067468388912953258795823785476;
            if(x>= 2.850000 &&  x < 3.000000 ) weight=7.341577057834249409040694445139;
            if(x>= 3.000000 &&  x < 3.150000 ) weight=8.812073861687228060191046097316;
            if(x>= 3.150000 &&  x < 3.300000 ) weight=10.618539863553172253318734874483;
            if(x>= 3.300000 &&  x < 3.450000 ) weight=12.685034855137985232431674376130;
            if(x>= 3.450000 &&  x < 3.600000 ) weight=15.005011593250202395211090333760;
            if(x>= 3.600000 &&  x < 3.750000 ) weight=17.725824581097811005747644230723;
            if(x>= 3.750000 &&  x < 3.900000 ) weight=20.871579377218040463048964738846;
            if(x>= 3.900000 &&  x < 4.050000 ) weight=24.360848683753186350031683105044;
            if(x>= 4.050000 &&  x < 4.200000 ) weight=28.488188673874866196911170845851;
            if(x>= 4.200000 &&  x < 4.350000 ) weight=32.938704676150628358755056979135;
            if(x>= 4.350000 &&  x < 4.500000 ) weight=38.004919741139133293472696095705;
            if(x>= 4.500000 &&  x < 4.650000 ) weight=43.905847481516879327045899117365;
            if(x>= 4.650000 &&  x < 4.800000 ) weight=50.316565609458564267697511240840;
            if(x>= 4.800000 &&  x < 4.950000 ) weight=57.346497358756437279225792735815;
            if(x>= 4.950000 &&  x < 5.100000 ) weight=65.527849013898290309043659362942;
            if(x>= 5.100000 &&  x < 5.250000 ) weight=74.047274070650502153512206859887;
            if(x>= 5.250000 &&  x < 5.400000 ) weight=83.615531495903567815730639267713;
            if(x>= 5.400000 &&  x < 5.550000 ) weight=94.362690320779350372504268307239;
            if(x>= 5.550000 &&  x < 5.700000 ) weight=106.045286496794076924743421841413;
            if(x>= 5.700000 &&  x < 5.850000 ) weight=119.012985953534311533985601272434;
            if(x>= 5.850000 &&  x < 6.000000 ) weight=132.880542792841708887863205745816;
            if(x>= 6.000000 &&  x < 6.150000 ) weight=147.654408358203852458245819434524;
            if(x>= 6.150000 &&  x < 6.300000 ) weight=164.856176227258970357070211321115;
            if(x>= 6.300000 &&  x < 6.450000 ) weight=182.443341289820949668865068815649;
            if(x>= 6.450000 &&  x < 6.600000 ) weight=201.580388007889069967859541065991;
            if(x>= 6.600000 &&  x < 6.750000 ) weight=222.096149387092708593627321533859;
            if(x>= 6.750000 &&  x < 6.900000 ) weight=245.731662642815820163377793505788;
            if(x>= 6.900000 &&  x < 7.050000 ) weight=269.694558759388030466652708128095;
            if(x>= 7.050000 &&  x < 7.200000 ) weight=296.517907159985099951882148161530;
            if(x>= 7.200000 &&  x < 7.350000 ) weight=323.392152451308618310577003285289;
            if(x>= 7.350000 &&  x < 7.500000 ) weight=352.243120627169730596506269648671;
            if(x>= 7.500000 &&  x < 7.650000 ) weight=385.261158912064445303258253261447;
            if(x>= 7.650000 &&  x < 7.800000 ) weight=418.118715856623225590738002210855;
            if(x>= 7.800000 &&  x < 7.950000 ) weight=455.113425937699730638996697962284;
            if(x>= 7.950000 &&  x < 8.100000 ) weight=495.331199298804051522893132641912;
            if(x>= 8.100000 &&  x < 8.250000 ) weight=535.194415493466749467188492417336;
            if(x>= 8.250000 &&  x < 8.400000 ) weight=580.109285047343064434244297444820;
            if(x>= 8.400000 &&  x < 8.550000 ) weight=626.069818030219266802305355668068;
            if(x>= 8.550000 &&  x < 8.700000 ) weight=676.485514594400115129246842116117;
            if(x>= 8.700000 &&  x < 8.850000 ) weight=727.181488442636236868565902113914;
            if(x>= 8.850000 &&  x < 9.000000 ) weight=780.772897612080214457819238305092;
            if(x>= 9.000000 &&  x < 9.150000 ) weight=843.460956984931044644326902925968;
            if(x>= 9.150000 &&  x < 9.300000 ) weight=903.110782152460615179734304547310;
            if(x>= 9.300000 &&  x < 9.450000 ) weight=966.468528833553591539384797215462;
            if(x>= 9.450000 &&  x < 9.600000 ) weight=1036.956504896397518677986226975918;
            if(x>= 9.600000 &&  x < 9.750000 ) weight=1100.209142769169375242199748754501;
            if(x>= 9.750000 &&  x < 9.900000 ) weight=1180.666364373458009140449576079845;
            if(x>= 9.900000 &&  x < 10.050000 ) weight=1262.838652389473509174422360956669;
            if(x>= 10.050000 &&  x < 10.200000 ) weight=1345.690689151621654673363082110882;
            if(x>= 10.200000 &&  x < 10.350000 ) weight=1441.650245938341868168208748102188;
            if(x>= 10.350000 &&  x < 10.500000 ) weight=1522.357306183298987889429554343224;
            if(x>= 10.500000 &&  x < 10.650000 ) weight=1622.668605949874290672596544027328;
            if(x>= 10.650000 &&  x < 10.800000 ) weight=1717.472191984062419578549452126026;
            if(x>= 10.800000 &&  x < 10.950000 ) weight=1818.572721180298685794696211814880;
            if(x>= 10.950000 &&  x < 11.100000 ) weight=1943.911193132668358884984627366066;
            if(x>= 11.100000 &&  x < 11.250000 ) weight=2054.165600845920835126889869570732;
            if(x>= 11.250000 &&  x < 11.400000 ) weight=2176.172676436633082630578428506851;
            if(x>= 11.400000 &&  x < 11.550000 ) weight=2288.558854875290307973045855760574;
            if(x>= 11.550000 &&  x < 11.700000 ) weight=2413.537715361238042532932013273239;
            if(x>= 11.700000 &&  x < 11.850000 ) weight=2557.627323171029274817556142807007;
            if(x>= 11.850000 &&  x < 12.000000 ) weight=2706.099549780678898969199508428574;
            if(x>= 12.000000 &&  x < 12.150000 ) weight=2851.892789594861369550926610827446;
            if(x>= 12.150000 &&  x < 12.300000 ) weight=2986.992952172975492430850863456726;
            if(x>= 12.300000 &&  x < 12.450000 ) weight=3150.172253940655991755193099379539;
            if(x>= 12.450000 &&  x < 12.600000 ) weight=3356.781323179806804546387866139412;
            if(x>= 12.600000 &&  x < 12.750000 ) weight=3509.653394369353463844163343310356;
            if(x>= 12.750000 &&  x < 12.900000 ) weight=3703.187316206644027261063456535339;
            if(x>= 12.900000 &&  x < 13.050000 ) weight=3874.919703406244934740243479609489;
            if(x>= 13.050000 &&  x < 13.200000 ) weight=4088.973986739532392675755545496941;
            if(x>= 13.200000 &&  x < 13.350000 ) weight=4319.189309648509151884354650974274;
            if(x>= 13.350000 &&  x < 13.500000 ) weight=4474.371121178378416516352444887161;
            if(x>= 13.500000 &&  x < 13.650000 ) weight=4710.164143890304330852814018726349;
            if(x>= 13.650000 &&  x < 13.800000 ) weight=4962.722731031956755032297223806381;
            if(x>= 13.800000 &&  x < 13.950000 ) weight=5212.695189196409955911803990602493;
            if(x>= 13.950000 &&  x < 14.100000 ) weight=5455.828904156822318327613174915314;
            if(x>= 14.100000 &&  x < 14.250000 ) weight=5747.811890884755484876222908496857;
            if(x>= 14.250000 &&  x < 14.400000 ) weight=6027.093776399433409096673130989075;
            if(x>= 14.400000 &&  x < 14.550000 ) weight=6247.732360090687507181428372859955;
            if(x>= 14.550000 &&  x < 14.700000 ) weight=6571.407631286909236223436892032623;
            if(x>= 14.700000 &&  x < 14.850000 ) weight=6886.337364156141120474785566329956;
            if(x>= 14.850000 &&  x < 15.000000 ) weight=7217.840294404419182683341205120087;
            if(x>= 15.000000 &&  x < 15.150000 ) weight=7618.389897352722073264885693788528;
            if(x>= 15.150000 &&  x < 15.300000 ) weight=7868.994076338297418260481208562851;
            if(x>= 15.300000 &&  x < 15.450000 ) weight=8257.157222000565525377169251441956;
            if(x>= 15.450000 &&  x < 15.600000 ) weight=8670.134419143883860670030117034912;
            if(x>= 15.600000 &&  x < 15.750000 ) weight=9084.681221171831566607579588890076;
            if(x>= 15.750000 &&  x < 15.900000 ) weight=9473.392376507903463789261877536774;
            if(x>= 15.900000 &&  x < 16.050000 ) weight=9927.918276660913761588744819164276;
            if(x>= 16.050000 &&  x < 16.200000 ) weight=10377.021055909193819388747215270996;
            if(x>= 16.200000 &&  x < 16.350000 ) weight=10864.833881717548138112761080265045;
            if(x>= 16.350000 &&  x < 16.500000 ) weight=11361.390229444130454794503748416901;
            if(x>= 16.500000 &&  x < 16.650000 ) weight=11813.237809490665313205681741237640;
            if(x>= 16.650000 &&  x < 16.800000 ) weight=12308.067018938214459922164678573608;
            if(x>= 16.800000 &&  x < 16.950000 ) weight=12924.491625300110172247514128684998;
            if(x>= 16.950000 &&  x < 17.100000 ) weight=13449.865155071940534980967640876770;
            if(x>= 17.100000 &&  x < 17.250000 ) weight=14101.318650777528091566637158393860;
            if(x>= 17.250000 &&  x < 17.400000 ) weight=14677.964300364441442070528864860535;
            if(x>= 17.400000 &&  x < 17.550000 ) weight=15385.504008052823337493464350700378;
            if(x>= 17.550000 &&  x < 17.700000 ) weight=15999.518692071294935885816812515259;
            if(x>= 17.700000 &&  x < 17.850000 ) weight=16762.538547729142010211944580078125;
            if(x>= 17.850000 &&  x < 18.000000 ) weight=17442.520843171583692310377955436707;
            if(x>= 18.000000 &&  x < 18.150000 ) weight=18191.438214885354682337492704391479;
            if(x>= 18.150000 &&  x < 18.300000 ) weight=19030.792132568531087599694728851318;
            if(x>= 18.300000 &&  x < 18.450000 ) weight=19851.796984364278614521026611328125;
            if(x>= 18.450000 &&  x < 18.600000 ) weight=20659.185993013343249913305044174194;
            if(x>= 18.600000 &&  x < 18.750000 ) weight=21717.881632693373830989003181457520;
            if(x>= 18.750000 &&  x < 18.900000 ) weight=22720.156667120721976971253752708435;
            if(x>= 18.900000 &&  x < 19.050000 ) weight=23591.610302217453863704577088356018;
            if(x>= 19.050000 &&  x < 19.200000 ) weight=24544.346729397140734363347291946411;
            if(x>= 19.200000 &&  x < 19.350000 ) weight=25610.402231535714236088097095489502;
            if(x>= 19.350000 &&  x < 19.500000 ) weight=26672.140630630929081235080957412720;
            if(x>= 19.500000 &&  x < 19.650000 ) weight=27845.838049204277922399342060089111;
            if(x>= 19.650000 &&  x < 19.800000 ) weight=29256.513710316798096755519509315491;
            if(x>= 19.800000 &&  x < 19.950000 ) weight=30371.685495973357319599017500877380;
            if(x>= 19.950000 &&  x < 20.100000 ) weight=31844.736605136160505935549736022949;
            if(x>= 20.100000 &&  x < 20.250000 ) weight=32877.035440022846159990876913070679;
            if(x>= 20.250000 &&  x < 20.400000 ) weight=34516.841773331318108830600976943970;
            if(x>= 20.400000 &&  x < 20.550000 ) weight=36083.787926580851490143686532974243;
            if(x>= 20.550000 &&  x < 20.700000 ) weight=37415.988561684040178079158067703247;
            if(x>= 20.700000 &&  x < 20.850000 ) weight=39030.587295537974569015204906463623;
            if(x>= 20.850000 &&  x < 21.000000 ) weight=40657.588691992459644097834825515747;
            if(x>= 21.000000 &&  x < 21.150000 ) weight=42783.512493353933678008615970611572;
            if(x>= 21.150000 &&  x < 21.300000 ) weight=44268.586319730238756164908409118652;
            if(x>= 21.300000 &&  x < 21.450000 ) weight=46575.167703819803136866539716720581;
            if(x>= 21.450000 &&  x < 21.600000 ) weight=48189.231645005631435196846723556519;
            if(x>= 21.600000 &&  x < 21.750000 ) weight=50150.034718143819191027432680130005;
            if(x>= 21.750000 &&  x < 21.900000 ) weight=51927.931111845769919455051422119141;
            if(x>= 21.900000 &&  x < 22.050000 ) weight=54636.474645266687730327248573303223;
            if(x>= 22.050000 &&  x < 22.200000 ) weight=56816.398614056255610194057226181030;
            if(x>= 22.200000 &&  x < 22.350000 ) weight=58973.583427455858327448368072509766;
            if(x>= 22.350000 &&  x < 22.500000 ) weight=61755.543408548677689395844936370850;
            if(x>= 22.500000 &&  x < 22.650000 ) weight=64208.755480811916640959680080413818;
            if(x>= 22.650000 &&  x < 22.800000 ) weight=66551.974592001977725885808467864990;
            if(x>= 22.800000 &&  x < 22.950000 ) weight=69338.159464728945749811828136444092;
            if(x>= 22.950000 &&  x < 23.100000 ) weight=71973.643574706831714138388633728027;
            if(x>= 23.100000 &&  x < 23.250000 ) weight=74936.432684315950609743595123291016;
            if(x>= 23.250000 &&  x < 23.400000 ) weight=78037.212588223876082338392734527588;
            if(x>= 23.400000 &&  x < 23.550000 ) weight=80977.782036187811172567307949066162;
            if(x>= 23.550000 &&  x < 23.700000 ) weight=83696.716434489077073521912097930908;
            if(x>= 23.700000 &&  x < 23.850000 ) weight=87407.139563519041985273361206054688;
            if(x>= 23.850000 &&  x < 24.000000 ) weight=90733.192175820004194974899291992188;
            if(x>= 24.000000 &&  x < 24.150000 ) weight=94085.642081807862268760800361633301;
            if(x>= 24.150000 &&  x < 24.300000 ) weight=97672.646482858079252764582633972168;
            if(x>= 24.300000 &&  x < 24.450000 ) weight=101836.834366340306587517261505126953;
            if(x>= 24.450000 &&  x < 24.600000 ) weight=105447.864113618183182552456855773926;
            if(x>= 24.600000 &&  x < 24.750000 ) weight=110145.801423295109998434782028198242;
            if(x>= 24.750000 &&  x < 24.900000 ) weight=113584.949773772576008923351764678955;
            if(x>= 24.900000 &&  x < 25.050000 ) weight=118302.821438762286561541259288787842;
            if(x>= 25.050000 &&  x < 25.200000 ) weight=121273.432330085255671292543411254883;
            if(x>= 25.200000 &&  x < 25.350000 ) weight=125523.107626564087695442140102386475;
            if(x>= 25.350000 &&  x < 25.500000 ) weight=130433.129234553649439476430416107178;
            if(x>= 25.500000 &&  x < 25.650000 ) weight=135284.289055481407558545470237731934;
            if(x>= 25.650000 &&  x < 25.800000 ) weight=139855.187847698485711589455604553223;
            if(x>= 25.800000 &&  x < 25.950000 ) weight=145256.667306432180339470505714416504;
            if(x>= 25.950000 &&  x < 26.100000 ) weight=150857.108721479278756305575370788574;
            if(x>= 26.100000 &&  x < 26.250000 ) weight=155763.772786196263041347265243530273;
            if(x>= 26.250000 &&  x < 26.400000 ) weight=161249.083903852559160441160202026367;
            if(x>= 26.400000 &&  x < 26.550000 ) weight=168117.625096614996436983346939086914;
            if(x>= 26.550000 &&  x < 26.700000 ) weight=172168.221335684007499366998672485352;
            if(x>= 26.700000 &&  x < 26.850000 ) weight=179201.009541476814774796366691589355;
            if(x>= 26.850000 &&  x < 27.000000 ) weight=184691.624567136197583749890327453613;
            if(x>= 27.000000 &&  x < 27.150000 ) weight=190189.644101943325949832797050476074;
            if(x>= 27.150000 &&  x < 27.300000 ) weight=197708.960882001003483310341835021973;
            if(x>= 27.300000 &&  x < 27.450000 ) weight=204437.298149541282327845692634582520;
            if(x>= 27.450000 &&  x < 27.600000 ) weight=210117.450069875980261713266372680664;
            if(x>= 27.600000 &&  x < 27.750000 ) weight=218221.452991388534428551793098449707;
            if(x>= 27.750000 &&  x < 27.900000 ) weight=227323.858492963161552324891090393066;
            if(x>= 27.900000 &&  x < 28.050000 ) weight=232577.545975567889399826526641845703;
            if(x>= 28.050000 &&  x < 28.200000 ) weight=242159.628700424829730764031410217285;
            if(x>= 28.200000 &&  x < 28.350000 ) weight=248470.437760056491242721676826477051;
            if(x>= 28.350000 &&  x < 28.500000 ) weight=255491.113505635614274069666862487793;
            if(x>= 28.500000 &&  x < 28.650000 ) weight=264586.516716481186449527740478515625;
            if(x>= 28.650000 &&  x < 28.800000 ) weight=272049.106330239854287356138229370117;
            if(x>= 28.800000 &&  x < 28.950000 ) weight=280514.584537865419406443834304809570;
            if(x>= 28.950000 &&  x < 29.100000 ) weight=293411.364126523665618151426315307617;
            if(x>= 29.100000 &&  x < 29.250000 ) weight=299323.141511107154656201601028442383;
            if(x>= 29.250000 &&  x < 29.400000 ) weight=309205.560403771232813596725463867188;
            if(x>= 29.400000 &&  x < 29.550000 ) weight=321390.270140645850915461778640747070;
            if(x>= 29.550000 &&  x < 29.700000 ) weight=329480.952456880826503038406372070312;
            if(x>= 29.700000 &&  x < 29.850000 ) weight=338103.102605441992636770009994506836;
            if(x>= 29.850000 &&  x < 30.000000 ) weight=348397.732848623592872172594070434570;
            if(x> 30.000000 ) weight=348397.732848623592872172594070434570;
        }
        
        //eta
        else if(pdg_particle==221)
        {
            if(x>= 0.000000 &&  x < 0.150000 ) weight=2.607013439379831876152593395091;
            if(x>= 0.150000 &&  x < 0.300000 ) weight=1.787246000906128617913282141672;
            if(x>= 0.300000 &&  x < 0.450000 ) weight=1.511981376652269171145803738909;
            if(x>= 0.450000 &&  x < 0.600000 ) weight=1.286308872388486124194173498836;
            if(x>= 0.600000 &&  x < 0.750000 ) weight=1.107481267390605150779947507544;
            if(x>= 0.750000 &&  x < 0.900000 ) weight=0.983018842055391695033961241279;
            if(x>= 0.900000 &&  x < 1.050000 ) weight=0.904161514569871038737858270906;
            if(x>= 1.050000 &&  x < 1.200000 ) weight=0.869832490575683014988328523032;
            if(x>= 1.200000 &&  x < 1.350000 ) weight=0.872583717444833983911678387813;
            if(x>= 1.350000 &&  x < 1.500000 ) weight=0.911565039819943079280051279056;
            if(x>= 1.500000 &&  x < 1.650000 ) weight=0.975231405941574758244883014413;
            if(x>= 1.650000 &&  x < 1.800000 ) weight=1.070131730327214203413177529001;
            if(x>= 1.800000 &&  x < 1.950000 ) weight=1.198776473473240899991765218147;
            if(x>= 1.950000 &&  x < 2.100000 ) weight=1.356288539982308805065258638933;
            if(x>= 2.100000 &&  x < 2.250000 ) weight=1.554186035072071758733613933146;
            if(x>= 2.250000 &&  x < 2.400000 ) weight=1.801105813416342460442365336348;
            if(x>= 2.400000 &&  x < 2.550000 ) weight=2.115664727510948139865831763018;
            if(x>= 2.550000 &&  x < 2.700000 ) weight=2.490247534275529250891167976079;
            if(x>= 2.700000 &&  x < 2.850000 ) weight=2.919512409247132289635828783503;
            if(x>= 2.850000 &&  x < 3.000000 ) weight=3.448718623375256697727309074253;
            if(x>= 3.000000 &&  x < 3.150000 ) weight=4.085821854664533958612082642503;
            if(x>= 3.150000 &&  x < 3.300000 ) weight=4.837239314627054476147804962238;
            if(x>= 3.300000 &&  x < 3.450000 ) weight=5.719576299297343346950128761819;
            if(x>= 3.450000 &&  x < 3.600000 ) weight=6.731485418529241648855077073677;
            if(x>= 3.600000 &&  x < 3.750000 ) weight=7.974702024985921511301967257168;
            if(x>= 3.750000 &&  x < 3.900000 ) weight=9.356885402495796810740102955606;
            if(x>= 3.900000 &&  x < 4.050000 ) weight=10.977954892290387789444139343686;
            if(x>= 4.050000 &&  x < 4.200000 ) weight=12.928445305499170814300669007935;
            if(x>= 4.200000 &&  x < 4.350000 ) weight=15.102790395218876895455650810618;
            if(x>= 4.350000 &&  x < 4.500000 ) weight=17.491979839702356258612780948170;
            if(x>= 4.500000 &&  x < 4.650000 ) weight=20.551324095288894255872946814634;
            if(x>= 4.650000 &&  x < 4.800000 ) weight=23.832444884096748438651047763415;
            if(x>= 4.800000 &&  x < 4.950000 ) weight=27.640620950409040545991956605576;
            if(x>= 4.950000 &&  x < 5.100000 ) weight=31.590453505559374036693043308333;
            if(x>= 5.100000 &&  x < 5.250000 ) weight=36.638488193189267860816471511498;
            if(x>= 5.250000 &&  x < 5.400000 ) weight=41.719615164071718993454851442948;
            if(x>= 5.400000 &&  x < 5.550000 ) weight=47.984489686038052980165957706049;
            if(x>= 5.550000 &&  x < 5.700000 ) weight=54.315247031924201337460544891655;
            if(x>= 5.700000 &&  x < 5.850000 ) weight=62.798772310475520441741537069902;
            if(x>= 5.850000 &&  x < 6.000000 ) weight=71.393158518553590852206980343908;
            if(x>= 6.000000 &&  x < 6.150000 ) weight=80.622124771637729168105579447001;
            if(x>= 6.150000 &&  x < 6.300000 ) weight=91.702053808570084925122500862926;
            if(x>= 6.300000 &&  x < 6.450000 ) weight=103.029975018339229109187726862729;
            if(x>= 6.450000 &&  x < 6.600000 ) weight=116.251935859840514808638545218855;
            if(x>= 6.600000 &&  x < 6.750000 ) weight=130.500478951039184494220535270870;
            if(x>= 6.750000 &&  x < 6.900000 ) weight=145.766943118013614366645924746990;
            if(x>= 6.900000 &&  x < 7.050000 ) weight=165.069622092609876062851981259882;
            if(x>= 7.050000 &&  x < 7.200000 ) weight=183.439451237236966107957414351404;
            if(x>= 7.200000 &&  x < 7.350000 ) weight=203.642690319174761270915041677654;
            if(x>= 7.350000 &&  x < 7.500000 ) weight=227.170526232498076524279895238578;
            if(x>= 7.500000 &&  x < 7.650000 ) weight=255.379316151416901448101270943880;
            if(x>= 7.650000 &&  x < 7.800000 ) weight=280.038994243541083051241002976894;
            if(x>= 7.800000 &&  x < 7.950000 ) weight=311.676204926161574348952854052186;
            if(x>= 7.950000 &&  x < 8.100000 ) weight=345.865987867621299756137887015939;
            if(x>= 8.100000 &&  x < 8.250000 ) weight=382.607356343874243975733406841755;
            if(x>= 8.250000 &&  x < 8.400000 ) weight=421.800556643482991603377740830183;
            if(x>= 8.400000 &&  x < 8.550000 ) weight=462.123578488440898581757210195065;
            if(x>= 8.550000 &&  x < 8.700000 ) weight=511.417289917822529332624981179833;
            if(x>= 8.700000 &&  x < 8.850000 ) weight=560.539627712715969209966715425253;
            if(x>= 8.850000 &&  x < 9.000000 ) weight=614.873142401480890839593484997749;
            if(x>= 9.000000 &&  x < 9.150000 ) weight=675.757278095039055187953636050224;
            if(x>= 9.150000 &&  x < 9.300000 ) weight=735.834297420662664990231860429049;
            if(x>= 9.300000 &&  x < 9.450000 ) weight=809.291165959088971249002497643232;
            if(x>= 9.450000 &&  x < 9.600000 ) weight=891.113928071516170348331797868013;
            if(x>= 9.600000 &&  x < 9.750000 ) weight=962.050308558426650051842443645000;
            if(x>= 9.750000 &&  x < 9.900000 ) weight=1040.726467983726706734159961342812;
            if(x>= 9.900000 &&  x < 10.050000 ) weight=1132.417394300917749205837026238441;
            if(x>= 10.050000 &&  x < 10.200000 ) weight=1237.242134882600794298923574388027;
            if(x>= 10.200000 &&  x < 10.350000 ) weight=1342.500948588073697465006262063980;
            if(x>= 10.350000 &&  x < 10.500000 ) weight=1459.701620782434929424198344349861;
            if(x>= 10.500000 &&  x < 10.650000 ) weight=1580.949867679561975819524377584457;
            if(x>= 10.650000 &&  x < 10.800000 ) weight=1699.710610093904278983245603740215;
            if(x>= 10.800000 &&  x < 10.950000 ) weight=1841.103734446948692493606358766556;
            if(x>= 10.950000 &&  x < 11.100000 ) weight=1988.412856448577258561272174119949;
            if(x>= 11.100000 &&  x < 11.250000 ) weight=2155.130319765818967425730079412460;
            if(x>= 11.250000 &&  x < 11.400000 ) weight=2311.939723316035724565153941512108;
            if(x>= 11.400000 &&  x < 11.550000 ) weight=2509.365523567474610899807885289192;
            if(x>= 11.550000 &&  x < 11.700000 ) weight=2666.786892999638894252711907029152;
            if(x>= 11.700000 &&  x < 11.850000 ) weight=2882.932045434884912538109347224236;
            if(x>= 11.850000 &&  x < 12.000000 ) weight=3098.431345687115026521496474742889;
            if(x>= 12.000000 &&  x < 12.150000 ) weight=3318.461154103130866133142262697220;
            if(x>= 12.150000 &&  x < 12.300000 ) weight=3561.833355517213931307196617126465;
            if(x>= 12.300000 &&  x < 12.450000 ) weight=3824.711382142780621506972238421440;
            if(x>= 12.450000 &&  x < 12.600000 ) weight=4124.000415043714383500628173351288;
            if(x>= 12.600000 &&  x < 12.750000 ) weight=4365.575846899258067423943430185318;
            if(x>= 12.750000 &&  x < 12.900000 ) weight=4719.798416369199003383982926607132;
            if(x>= 12.900000 &&  x < 13.050000 ) weight=5020.193922976839530747383832931519;
            if(x>= 13.050000 &&  x < 13.200000 ) weight=5372.197690875512307684402912855148;
            if(x>= 13.200000 &&  x < 13.350000 ) weight=5699.544998009871960675809532403946;
            if(x>= 13.350000 &&  x < 13.500000 ) weight=6134.392488418478933454025536775589;
            if(x>= 13.500000 &&  x < 13.650000 ) weight=6487.706140926851730910129845142365;
            if(x>= 13.650000 &&  x < 13.800000 ) weight=6903.037574226802462362684309482574;
            if(x>= 13.800000 &&  x < 13.950000 ) weight=7427.179071277792900218628346920013;
            if(x>= 13.950000 &&  x < 14.100000 ) weight=7819.318293623327008390333503484726;
            if(x>= 14.100000 &&  x < 14.250000 ) weight=8438.071455734352639410644769668579;
            if(x>= 14.250000 &&  x < 14.400000 ) weight=8789.399015917477299808524549007416;
            if(x>= 14.400000 &&  x < 14.550000 ) weight=9461.991986986633492051623761653900;
            if(x>= 14.550000 &&  x < 14.700000 ) weight=10051.828561189437095890752971172333;
            if(x>= 14.700000 &&  x < 14.850000 ) weight=10584.915246637565360288135707378387;
            if(x>= 14.850000 &&  x < 15.000000 ) weight=11267.159592693669765139929950237274;
            if(x>= 15.000000 &&  x < 15.150000 ) weight=11962.630896410348213976249098777771;
            if(x>= 15.150000 &&  x < 15.300000 ) weight=12621.696522844902574433945119380951;
            if(x>= 15.300000 &&  x < 15.450000 ) weight=13427.227611217056619352661073207855;
            if(x>= 15.450000 &&  x < 15.600000 ) weight=14131.633453236521745566278696060181;
            if(x>= 15.600000 &&  x < 15.750000 ) weight=15104.205798568980753771029412746429;
            if(x>= 15.750000 &&  x < 15.900000 ) weight=15927.992671955022160545922815799713;
            if(x>= 15.900000 &&  x < 16.050000 ) weight=16835.557825357074761996045708656311;
            if(x>= 16.050000 &&  x < 16.200000 ) weight=17557.688697662833874346688389778137;
            if(x>= 16.200000 &&  x < 16.350000 ) weight=18604.082638004285399802029132843018;
            if(x>= 16.350000 &&  x < 16.500000 ) weight=19551.366895061677496414631605148315;
            if(x>= 16.500000 &&  x < 16.650000 ) weight=20825.192958685602206969633698463440;
            if(x>= 16.650000 &&  x < 16.800000 ) weight=21848.545611173365614376962184906006;
            if(x>= 16.800000 &&  x < 16.950000 ) weight=22946.463071799153112806379795074463;
            if(x>= 16.950000 &&  x < 17.100000 ) weight=24183.065913520775211509317159652710;
            if(x>= 17.100000 &&  x < 17.250000 ) weight=25690.376357720280793728306889533997;
            if(x>= 17.250000 &&  x < 17.400000 ) weight=26979.187151967260433593764901161194;
            if(x>= 17.400000 &&  x < 17.550000 ) weight=28376.784267282273503951728343963623;
            if(x>= 17.550000 &&  x < 17.700000 ) weight=29774.573194361157220555469393730164;
            if(x>= 17.700000 &&  x < 17.850000 ) weight=31316.698672573380463290959596633911;
            if(x>= 17.850000 &&  x < 18.000000 ) weight=32628.808620870371669298037886619568;
            if(x>= 18.000000 &&  x < 18.150000 ) weight=34532.814676859503379091620445251465;
            if(x>= 18.150000 &&  x < 18.300000 ) weight=36639.596010188462969381362199783325;
            if(x>= 18.300000 &&  x < 18.450000 ) weight=38054.399591543384303804486989974976;
            if(x>= 18.450000 &&  x < 18.600000 ) weight=39898.590350935504829976707696914673;
            if(x>= 18.600000 &&  x < 18.750000 ) weight=41850.088052719816914759576320648193;
            if(x>= 18.750000 &&  x < 18.900000 ) weight=43862.632534997981565538793802261353;
            if(x>= 18.900000 &&  x < 19.050000 ) weight=45702.669671746487438213080167770386;
            if(x>= 19.050000 &&  x < 19.200000 ) weight=48289.448494917916832491755485534668;
            if(x>= 19.200000 &&  x < 19.350000 ) weight=50708.267688120380626060068607330322;
            if(x>= 19.350000 &&  x < 19.500000 ) weight=52594.793891625355172436684370040894;
            if(x>= 19.500000 &&  x < 19.650000 ) weight=54868.226217383067705668509006500244;
            if(x>= 19.650000 &&  x < 19.800000 ) weight=58019.566332522525044623762369155884;
            if(x>= 19.800000 &&  x < 19.950000 ) weight=61000.537363612951594404876232147217;
            if(x>= 19.950000 &&  x < 20.100000 ) weight=63139.129938553269312251359224319458;
            if(x>= 20.100000 &&  x < 20.250000 ) weight=66123.099508769795647822320461273193;
            if(x>= 20.250000 &&  x < 20.400000 ) weight=69788.492502370529109612107276916504;
            if(x>= 20.400000 &&  x < 20.550000 ) weight=72151.247086239280179142951965332031;
            if(x>= 20.550000 &&  x < 20.700000 ) weight=74783.421325135976076126098632812500;
            if(x>= 20.700000 &&  x < 20.850000 ) weight=78498.445060129743069410324096679688;
            if(x>= 20.850000 &&  x < 21.000000 ) weight=82994.656149944654316641390323638916;
            if(x>= 21.000000 &&  x < 21.150000 ) weight=85986.635790739441290497779846191406;
            if(x>= 21.150000 &&  x < 21.300000 ) weight=90219.603250956613919697701930999756;
            if(x>= 21.300000 &&  x < 21.450000 ) weight=93191.413883653251104988157749176025;
            if(x>= 21.450000 &&  x < 21.600000 ) weight=97705.202995941304834559559822082520;
            if(x>= 21.600000 &&  x < 21.750000 ) weight=101399.924614956413279287517070770264;
            if(x>= 21.750000 &&  x < 21.900000 ) weight=105515.724870105230365879833698272705;
            if(x>= 21.900000 &&  x < 22.050000 ) weight=110352.370892453691340051591396331787;
            if(x>= 22.050000 &&  x < 22.200000 ) weight=114820.565161417529452592134475708008;
            if(x>= 22.200000 &&  x < 22.350000 ) weight=118395.101410173578187823295593261719;
            if(x>= 22.350000 &&  x < 22.500000 ) weight=123984.958204261231003329157829284668;
            if(x>= 22.500000 &&  x < 22.650000 ) weight=128744.001995247468585148453712463379;
            if(x>= 22.650000 &&  x < 22.800000 ) weight=134532.759974437212804332375526428223;
            if(x>= 22.800000 &&  x < 22.950000 ) weight=140415.005636782327201217412948608398;
            if(x>= 22.950000 &&  x < 23.100000 ) weight=145332.324425637722015380859375000000;
            if(x>= 23.100000 &&  x < 23.250000 ) weight=152464.808183897112030535936355590820;
            if(x>= 23.250000 &&  x < 23.400000 ) weight=156397.875369209592463448643684387207;
            if(x>= 23.400000 &&  x < 23.550000 ) weight=163383.278671121777733787894248962402;
            if(x>= 23.550000 &&  x < 23.700000 ) weight=169408.483720124728279188275337219238;
            if(x>= 23.700000 &&  x < 23.850000 ) weight=176358.580724780389573425054550170898;
            if(x>= 23.850000 &&  x < 24.000000 ) weight=183755.458544989698566496372222900391;
            if(x>= 24.000000 &&  x < 24.150000 ) weight=190759.434472528402693569660186767578;
            if(x>= 24.150000 &&  x < 24.300000 ) weight=198138.960127636528341099619865417480;
            if(x>= 24.300000 &&  x < 24.450000 ) weight=205671.252855026308679953217506408691;
            if(x>= 24.450000 &&  x < 24.600000 ) weight=211954.498173586500342935323715209961;
            if(x>= 24.600000 &&  x < 24.750000 ) weight=221988.260443281324114650487899780273;
            if(x>= 24.750000 &&  x < 24.900000 ) weight=229269.115259110985789448022842407227;
            if(x>= 24.900000 &&  x < 25.050000 ) weight=236272.412393393809907138347625732422;
            if(x>= 25.050000 &&  x < 25.200000 ) weight=247536.035648414283059537410736083984;
            if(x>= 25.200000 &&  x < 25.350000 ) weight=256745.661791664140764623880386352539;
            if(x>= 25.350000 &&  x < 25.500000 ) weight=265959.407035561802331358194351196289;
            if(x>= 25.500000 &&  x < 25.650000 ) weight=272792.205935468897223472595214843750;
            if(x>= 25.650000 &&  x < 25.800000 ) weight=284958.879410992725752294063568115234;
            if(x>= 25.800000 &&  x < 25.950000 ) weight=293135.933526697917841374874114990234;
            if(x>= 25.950000 &&  x < 26.100000 ) weight=306789.669870740908663719892501831055;
            if(x>= 26.100000 &&  x < 26.250000 ) weight=314656.642340345308184623718261718750;
            if(x>= 26.250000 &&  x < 26.400000 ) weight=325491.535578076262027025222778320312;
            if(x>= 26.400000 &&  x < 26.550000 ) weight=335557.709334257640875875949859619141;
            if(x>= 26.550000 &&  x < 26.700000 ) weight=348618.967548974556848406791687011719;
            if(x>= 26.700000 &&  x < 26.850000 ) weight=362826.927372189005836844444274902344;
            if(x>= 26.850000 &&  x < 27.000000 ) weight=372229.699274422135204076766967773438;
            if(x>= 27.000000 &&  x < 27.150000 ) weight=382681.114112682989798486232757568359;
            if(x>= 27.150000 &&  x < 27.300000 ) weight=398779.963024294236674904823303222656;
            if(x>= 27.300000 &&  x < 27.450000 ) weight=411694.789588931133039295673370361328;
            if(x>= 27.450000 &&  x < 27.600000 ) weight=426843.753943774267099797725677490234;
            if(x>= 27.600000 &&  x < 27.750000 ) weight=441188.173102431930601596832275390625;
            if(x>= 27.750000 &&  x < 27.900000 ) weight=457830.728830497420858591794967651367;
            if(x>= 27.900000 &&  x < 28.050000 ) weight=476308.351607583230361342430114746094;
            if(x>= 28.050000 &&  x < 28.200000 ) weight=488209.120811497559770941734313964844;
            if(x>= 28.200000 &&  x < 28.350000 ) weight=503660.894718579307664185762405395508;
            if(x>= 28.350000 &&  x < 28.500000 ) weight=515678.671460700104944407939910888672;
            if(x>= 28.500000 &&  x < 28.650000 ) weight=532079.018234341405332088470458984375;
            if(x>= 28.650000 &&  x < 28.800000 ) weight=552077.714736754074692726135253906250;
            if(x>= 28.800000 &&  x < 28.950000 ) weight=565908.432159291580319404602050781250;
            if(x>= 28.950000 &&  x < 29.100000 ) weight=589745.260512855718843638896942138672;
            if(x>= 29.100000 &&  x < 29.250000 ) weight=606388.982299502473324537277221679688;
            if(x>= 29.250000 &&  x < 29.400000 ) weight=623132.245116776903159916400909423828;
            if(x>= 29.400000 &&  x < 29.550000 ) weight=648147.859134702244773507118225097656;
            if(x>= 29.550000 &&  x < 29.700000 ) weight=665900.618200576282106339931488037109;
            if(x>= 29.700000 &&  x < 29.850000 ) weight=684239.760351214441470801830291748047;
            if(x>= 29.850000 &&  x < 30.000000 ) weight=702028.128093079198151826858520507812;
            if(x>30.000000 ) weight=702028.128093079198151826858520507812;
            
            
        }
        
        else weight=1;
        
    }//close fUseTrigger
    
    return weight/40000;
    
}
Double_t AliAnalysisTaskHFEpACorrelation::SetEoverPCutPtDependentMC(Double_t pt)
{
    
    fEoverPCutMin=0.8;
    fEoverPCutMax=1.2;
    
    /*
     ////======================================================================================
     ////======================================================================================
     ////============================= data 2 sigma ===========================================
     if(!fIsMC){
     //data 2 sigma
     if(pt>= 2.00 &&  pt < 2.50){
     fEoverPCutMin=0.7719;
     fEoverPCutMax=1.1023;
     }
     if(pt>= 2.50 &&  pt < 3.00){
     fEoverPCutMin=0.7966;
     fEoverPCutMax=1.1088;
     }
     if(pt>= 3.00 &&  pt < 4.00){
     fEoverPCutMin=0.8175;
     fEoverPCutMax=1.1101;
     }
     if(pt>= 4.00 &&  pt < 5.00){
     fEoverPCutMin=0.8302;
     fEoverPCutMax=1.1194;
     }
     if(pt>= 5.00 &&  pt < 6.00){
     fEoverPCutMin=0.8517;
     fEoverPCutMax=1.1177;
     }
     if(pt>= 6.00 &&  pt < 8.00){
     fEoverPCutMin=0.8901;
     fEoverPCutMax=1.1139;
     }
     if(pt>= 8.00 &&  pt < 10.00){
     fEoverPCutMin=0.8703;
     fEoverPCutMax=1.1377;
     }
     if(pt>= 10.00 &&  pt < 12.00){
     fEoverPCutMin=0.9043;
     fEoverPCutMax=1.1977;
     }
     }
     
     ////=========================================== MC 2 sigma ===========================================
     
     
     if(fIsMC){
     
     if(pt>= 2.00 &&  pt < 2.50){
     fEoverPCutMin=0.7712;
     fEoverPCutMax=1.0746;
     }
     if(pt>= 2.50 &&  pt < 3.00){
     fEoverPCutMin=0.7946;
     fEoverPCutMax=1.0708;
     }
     if(pt>= 3.00 &&  pt < 4.00){
     fEoverPCutMin=0.8196;
     fEoverPCutMax=1.0678;
     }
     if(pt>= 4.00 &&  pt < 5.00){
     fEoverPCutMin=0.8396;
     fEoverPCutMax=1.0646;
     }
     if(pt>= 5.00 &&  pt < 6.00){
     fEoverPCutMin=0.8527;
     fEoverPCutMax=1.0647;
     }
     if(pt>= 6.00 &&  pt < 8.00){
     fEoverPCutMin=0.8652;
     fEoverPCutMax=1.0670;
     }
     if(pt>= 8.00 &&  pt < 10.00){
     fEoverPCutMin=0.8748;
     fEoverPCutMax=1.0804;
     }
     if(pt>= 10.00 &&  pt < 12.00){
     fEoverPCutMin=0.8814;
     fEoverPCutMax=1.0842;
     }
     
     }
     */
    ////======================================================================================
    ////======================================================================================
    ////=========================== data 1.5 sigma ===========================================
    /*
     if(!fIsMC){
     //data 1.5 sigmas
     if(pt>= 2.00 &&  pt < 2.50){
     fEoverPCutMin=0.8132;
     fEoverPCutMax=1.0610;
     }
     if(pt>= 2.50 &&  pt < 3.00){
     fEoverPCutMin=0.8356;
     fEoverPCutMax=1.0698;
     }
     if(pt>= 3.00 &&  pt < 4.00){
     fEoverPCutMin=0.8541;
     fEoverPCutMax=1.0735;
     }
     if(pt>= 4.00 &&  pt < 5.00){
     fEoverPCutMin=0.8664;
     fEoverPCutMax=1.0832;
     }
     if(pt>= 5.00 &&  pt < 6.00){
     fEoverPCutMin=0.8849;
     fEoverPCutMax=1.0845;
     }
     if(pt>= 6.00 &&  pt < 8.00){
     fEoverPCutMin=0.9180;
     fEoverPCutMax=1.0860;
     }
     if(pt>= 8.00 &&  pt < 10.00){
     fEoverPCutMin=0.9037;
     fEoverPCutMax=1.1043;
     }
     if(pt>= 10.00 &&  pt < 12.00){
     fEoverPCutMin=0.9409;
     fEoverPCutMax=1.1611;
     }
     }
     
     ////======================================================================================
     ////======================================================================================
     ////=========================== MC 1.5 sigma ===========================================
     
     if(fIsMC){
     //MC 1.5 sigmas
     if(pt>= 2.00 &&  pt < 2.50){
     fEoverPCutMin=0.8091;
     fEoverPCutMax=1.0367;
     }
     if(pt>= 2.50 &&  pt < 3.00){
     fEoverPCutMin=0.8292;
     fEoverPCutMax=1.0362;
     }
     if(pt>= 3.00 &&  pt < 4.00){
     fEoverPCutMin=0.8506;
     fEoverPCutMax=1.0368;
     }
     if(pt>= 4.00 &&  pt < 5.00){
     fEoverPCutMin=0.8677;
     fEoverPCutMax=1.0365;
     }
     if(pt>= 5.00 &&  pt < 6.00){
     fEoverPCutMin=0.8792;
     fEoverPCutMax=1.0382;
     }
     if(pt>= 6.00 &&  pt < 8.00){
     fEoverPCutMin=0.8904;
     fEoverPCutMax=1.0418;
     }
     if(pt>= 8.00 &&  pt < 10.00){
     fEoverPCutMin=0.9005;
     fEoverPCutMax=1.0547;
     }
     if(pt>= 10.00 &&  pt < 12.00){
     fEoverPCutMin=0.9067;
     fEoverPCutMax=1.0589;
     }
     
     }
     */
    ////======================================================================================
    ////======================================================================================
    ////=========================== data 2.5 sigma ===========================================
    
    if(!fIsMC){
        //data 2.5 sigmas
        if(pt>= 2.00 &&  pt < 2.50){
            fEoverPCutMin=0.7306;
            fEoverPCutMax=1.1436;
        }
        if(pt>= 2.50 &&  pt < 3.00){
            fEoverPCutMin=0.7575;
            fEoverPCutMax=1.1479;
        }
        if(pt>= 3.00 &&  pt < 4.00){
            fEoverPCutMin=0.7809;
            fEoverPCutMax=1.1467;
        }
        if(pt>= 4.00 &&  pt < 5.00){
            fEoverPCutMin=0.7941;
            fEoverPCutMax=1.1555;
        }
        if(pt>= 5.00 &&  pt < 6.00){
            fEoverPCutMin=0.8184;
            fEoverPCutMax=1.1510;
        }
        if(pt>= 6.00 &&  pt < 8.00){
            fEoverPCutMin=0.8621;
            fEoverPCutMax=1.1419;
        }
        if(pt>= 8.00 &&  pt < 10.00){
            fEoverPCutMin=0.8368;
            fEoverPCutMax=1.1712;
        }
        if(pt>= 10.00 &&  pt < 12.00){
            fEoverPCutMin=0.8676;
            fEoverPCutMax=1.2344;
        }
        
    }
    
    ////======================================================================================
    ////======================================================================================
    ////=========================== MC 2.5 sigma ===========================================
    
    if(fIsMC){
        //MC 2.5 sigmas
        if(pt>= 2.00 &&  pt < 2.50){
            fEoverPCutMin=0.7333;
            fEoverPCutMax=1.1125;
        }
        if(pt>= 2.50 &&  pt < 3.00){
            fEoverPCutMin=0.7601;
            fEoverPCutMax=1.1053;
        }
        if(pt>= 3.00 &&  pt < 4.00){
            fEoverPCutMin=0.7885;
            fEoverPCutMax=1.0989;
        }
        if(pt>= 4.00 &&  pt < 5.00){
            fEoverPCutMin=0.8115;
            fEoverPCutMax=1.0927;
        }
        if(pt>= 5.00 &&  pt < 6.00){
            fEoverPCutMin=0.8262;
            fEoverPCutMax=1.0912;
        }
        if(pt>= 6.00 &&  pt < 8.00){
            fEoverPCutMin=0.8400;
            fEoverPCutMax=1.0922;
        }
        if(pt>= 8.00 &&  pt < 10.00){
            fEoverPCutMin=0.8492;
            fEoverPCutMax=1.1060;
        }
        if(pt>= 10.00 &&  pt < 12.00){
            fEoverPCutMin=0.8560;
            fEoverPCutMax=1.1096;
        }
    }
    
    
    return fEoverPCutMin;
    return fEoverPCutMax;
    
}

Double_t AliAnalysisTaskHFEpACorrelation::GetHadronWCorrection(Double_t pt)
{
    Double_t Bins[] = {0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
    Int_t NumberofBins = 15;
    //Compatibility with old analysis
    Double_t HadronWeight = 1.0;
    

    
    //"new Default" analysis
    if (fDCAcutrHadron == 0.25)
    {
        if (fCentralityMin == 0)
        {
            Double_t HadronEfficiency[] = {0.8202776, 0.8292172, 0.8325606, 0.8336015, 0.8340769, 0.8343634, 0.8341395, 0.8338822, 0.8333004, 0.8319232,0.8291333, 0.8251111, 0.8183794,0.8098981, 0.8009887};
            
            for(Int_t i = 0; i <NumberofBins; i++)
            {
                if ( (pt > Bins[i]) && (pt<Bins[i+1]))
                    HadronWeight = (1.0/HadronEfficiency[i]);
            }
            
        }
        if (fCentralityMin == 20)
        {
            Double_t HadronEfficiency[] = {0.8214875, 0.8302184, 0.8336327, 0.8348229, 0.8352569, 0.8350877, 0.8349991, 0.8346928, 0.8342593, 0.8326481, 0.8300479, 0.8258925, 0.8192205, 0.8104965, 0.8022359};
            
            for(Int_t i = 0; i <NumberofBins; i++)
            {
                if ( (pt > Bins[i]) && (pt<Bins[i+1]))
                    HadronWeight = (1.0/HadronEfficiency[i]);
            }
            
        }
        
        if (fCentralityMin == 60)
        {
            Double_t HadronEfficiency[] = {0.8230401,0.831619,0.8348156, 0.8356443, 0.8358399, 0.8360357, 0.8360404, 0.8354132, 0.8337131,0.8318167, 0.8296427, 0.8254036, 0.8188907, 0.8112107, 0.8027191};
                for(Int_t i = 0; i <NumberofBins; i++)
                {
                    if ( (pt > Bins[i]) && (pt<Bins[i+1]))
                        HadronWeight = (1.0/HadronEfficiency[i]);
                }

            
        }

    }
    //printf("Using pt: %f HW: %f\n",pt, HadronWeight);
    
    return HadronWeight; //other configurations not implemeted yet - need the update in the Config file;
    
    
}

