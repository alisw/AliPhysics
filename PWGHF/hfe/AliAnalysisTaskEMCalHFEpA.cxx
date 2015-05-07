
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
	//		version: April 22, 2015.								      //
	//                                                                    //
	//	    Authors 							                          //
	//		Elienos Pereira de Oliveira Filho (epereira@cern.ch)	      //
	//		Cristiane Jahnke		(cristiane.jahnke@cern.ch)		      //
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
#include "AliAnalysisTaskEMCalHFEpA.h"
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
ClassImp(AliAnalysisTaskEMCalHFEpA)

	//______________________________________________________________________
AliAnalysisTaskEMCalHFEpA::AliAnalysisTaskEMCalHFEpA(const char *name) 
: AliAnalysisTaskSE(name)
,fCorrelationFlag(0)
,fIsMC(0)
,fUseEMCal(kFALSE)
,fUseTrigger(kFALSE)
,fUseTender(kFALSE)
,fUseShowerShapeCut(kFALSE)
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
,fEoverP_pt_true_hadrons(0)
,fEoverP_pt_true_electrons0(0)
,fEoverP_pt_true_hadrons0(0)
,fEoverP_pt(0)
,fEoverP_tpc(0)
,fEoverP_tpc_p_trigger(0)
,fEoverP_tpc_pt_trigger(0)
,fTPC_pt(0)
,fTPC_p(0)
,fTPCnsigma_pt(0)
,fTPCnsigma_p(0)
,fTPCnsigma_p_TPC(0)
,fTPCnsigma_p_TPC_on_EMCal_acc(0)
,fTPCnsigma_p_TPC_EoverP_cut(0)
,fTPCnsigma_pt_2D(0)
,fShowerShapeCut(0)
,fShowerShapeM02_EoverP(0)
,fShowerShapeM20_EoverP(0)
,fShowerShape_ha(0)
,fShowerShape_ele(0)
,fTPCnsigma_eta(0)
,fTPCnsigma_phi(0)
,fECluster(0)
,fECluster_pure(0)
,fECluster_not_exotic(0)
,fECluster_not_exotic1(0)
,fECluster_not_exotic2(0)
,fECluster_exotic(0)
,fNCluster_pure(0)
,fNCluster_pure_aod(0)
,fNCluster_ECluster(0)
,fNcells_energy(0)
,fNcells_energy_elec_selected(0)
,fNcells_energy_not_exotic(0)

,fEtaPhi(0)
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
,fInvMass(0)
,fInvMassBack(0)
,fDCA(0)
,fDCABack(0)
,fOpAngle(0)
,fOpAngleBack(0)
,fInvMass2(0)
,fInvMassBack2(0)
,fDCA2(0)
,fDCABack2(0)
,fOpAngle2(0)
,fOpAngleBack2(0)
,fMassCut(0.1)
,fEtaCutMin(-0.9)
,fEtaCutMax(0.9)
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
//,fEMCALRecoUtils(0)//exotic

{
		//Named constructor 
		// Define input and output slots here
		// Input slot #0 works with a TChain
	
		//exotic
		//fEMCALRecoUtils  = new AliEMCALRecoUtils();
	
	DefineInput(0, TChain::Class());
		// Output slot #0 id reserved by the base class for AOD
		// Output slot #1 writes into a TH1 container
		// DefineOutput(1, TH1I::Class());
	DefineOutput(1, TList::Class());
		//  DefineOutput(3, TTree::Class());
}

	//________________________________________________________________________
AliAnalysisTaskEMCalHFEpA::AliAnalysisTaskEMCalHFEpA() 
: AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTaskEMCalHFEpA")
,fCorrelationFlag(0)
,fIsMC(0)
,fUseEMCal(kFALSE)
,fUseTrigger(kFALSE)
,fUseTender(kFALSE)
,fUseShowerShapeCut(kFALSE)
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
,fEoverP_pt_true_hadrons(0)
,fEoverP_pt_true_electrons0(0)
,fEoverP_pt_true_hadrons0(0)
,fEoverP_pt(0)
,fEoverP_tpc(0)
,fEoverP_tpc_p_trigger(0)
,fEoverP_tpc_pt_trigger(0)
,fTPC_pt(0)
,fTPC_p(0)
,fTPCnsigma_pt(0)
,fTPCnsigma_p(0)
,fTPCnsigma_p_TPC(0)
,fTPCnsigma_p_TPC_on_EMCal_acc(0)
,fTPCnsigma_p_TPC_EoverP_cut(0)
,fTPCnsigma_pt_2D(0)
,fShowerShapeCut(0)
,fShowerShapeM02_EoverP(0)
,fShowerShapeM20_EoverP(0)
,fShowerShape_ha(0)
,fShowerShape_ele(0)
,fTPCnsigma_eta(0)
,fTPCnsigma_phi(0)
,fECluster(0)
,fECluster_pure(0)
,fECluster_not_exotic(0)
,fECluster_not_exotic1(0)
,fECluster_not_exotic2(0)
,fECluster_exotic(0)
,fNCluster_pure(0)
,fNCluster_pure_aod(0)
,fNCluster_ECluster(0)
,fNcells_energy(0)
,fNcells_energy_elec_selected(0)
,fNcells_energy_not_exotic(0)
,fEtaPhi(0)
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
,fInvMass(0)
,fInvMassBack(0)
,fDCA(0)
,fDCABack(0)
,fOpAngle(0)
,fOpAngleBack(0)
,fInvMass2(0)
,fInvMassBack2(0)
,fDCA2(0)
,fDCABack2(0)
,fOpAngle2(0)
,fOpAngleBack2(0)
,fMassCut(0.1)
,fEtaCutMin(-0.9)
,fEtaCutMax(0.9)
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
	//,fEMCALRecoUtils(0)//exotic
{
		// Constructor
		// Define input and output slots here
		// Input slot #0 works with a TChain
	
		//exotic
		// fEMCALRecoUtils  = new AliEMCALRecoUtils();
	
	DefineInput(0, TChain::Class());
		// Output slot #0 id reserved by the base class for AOD
		// Output slot #1 writes into a TH1 container
		// DefineOutput(1, TH1I::Class());
	DefineOutput(1, TList::Class());
		//DefineOutput(3, TTree::Class());
}

	//______________________________________________________________________
AliAnalysisTaskEMCalHFEpA::~AliAnalysisTaskEMCalHFEpA()
{
	//Destructor 
	delete fOutputList;
	delete fPID;
	delete fCFM;
	delete fPIDqa;
	//Lucile
	//delete reader; 
	//if(fEMCALRecoUtils)   delete fEMCALRecoUtils ;
}

	//______________________________________________________________________
	//Create Output Objects
	//Here we can define the histograms and others output files
	//Called once
void AliAnalysisTaskEMCalHFEpA::UserCreateOutputObjects()
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
	
	fCuts->Initialize(fCFM);
		//______________________________________________________________________
	
		///___________________//Lucile
	
		//if(fIsAOD) {
		//	reader = new AliCaloTrackAODReader();
		//}
			//___________________________________________________
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
	
	fPtPrim = new TH1F("fPtPrim","Primary Electrons aod track; p_{T} (GeV/c); Count",300,0,30);
	fPtSec = new TH1F("fPtSec","Secundary Electrons aod track; p_{T} (GeV/c); Count",300,0,30);
	fPtPrim2 = new TH1F("fPtPrim2","Primary Electrons vtrack; p_{T} (GeV/c); Count",300,0,30);
	fPtSec2 = new TH1F("fPtSec2","Secundary Electrons vtrack; p_{T} (GeV/c); Count",300,0,30);

	
	fPtElec_ULS = new TH1F("fPtElec_ULS","ULS; p_{T} (GeV/c); Count",300,0,30);
	fPtElec_LS = new TH1F("fPtElec_LS","LS; p_{T} (GeV/c); Count",300,0,30);
	
	fPtElec_ULS_NoPid = new TH1F("fPtElec_ULS_NoPid","ULS; p_{T} (GeV/c); Count",300,0,30);
	fPtElec_LS_NoPid = new TH1F("fPtElec_LS_NoPid","LS; p_{T} (GeV/c); Count",300,0,30);
	
	fPtElec_ULS_MC = new TH1F("fPtElec_ULS_MC","ULS; p_{T} (GeV/c); Count",300,0,30);
	fPtElec_ULS_MC_weight = new TH1F("fPtElec_ULS_MC_weight","ULS; p_{T} (GeV/c); Count",300,0,30);
	
	fPtElec_ULS_weight = new TH1F("fPtElec_ULS_weight","ULS; p_{T} (GeV/c); Count",300,0,30);
	fPtElec_LS_weight = new TH1F("fPtElec_LS_weight","LS; p_{T} (GeV/c); Count",300,0,30);
	
	fTOF01 = new TH2F("fTOF01","",200,-20,20,200,-20,20);
	fTOF02 = new TH2F("fTOF02","",200,-20,20,200,-20,20);
	fTOF03 = new TH2F("fTOF03","",200,-20,20,200,-20,20);
	
	if(fFillBackground){
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
	
	//new histos for TPC signal -> Can be used for any p range
	fTPCnsigma_p_TPC = new TH2F("fTPCnsigma_p_TPC",";p (GeV/c);TPC Electron N#sigma",3000,0,30,1000,-15,10);
	fTPCnsigma_p_TPC_on_EMCal_acc = new TH2F("fTPCnsigma_p_TPC_on_EMCal_acc",";p (GeV/c);TPC Electron N#sigma",3000,0,30,1000,-15,10);
	fTPCnsigma_p_TPC_EoverP_cut = new TH2F("fTPCnsigma_p_TPC_EoverP_cut",";p (GeV/c);TPC Electron N#sigma",3000,0,30,1000,-15,10);
	
	
	fShowerShapeCut = new TH2F("fShowerShapeCut","Shower Shape;M02;M20",500,0,1.8,500,0,1.8);
	fEtaPhi_num=new TH2F("fEtaPhi_num","#eta x #phi track;#phi;#eta",200,0.,5,50,-1.,1.);
	fEtaPhi_den=new TH2F("fEtaPhi_den","#eta x #phi track;#phi;#eta",200,0.,5,50,-1.,1.);
	fEtaPhi_data=new TH2F("fEtaPhi_data","#eta x #phi track;#phi;#eta",200,0.,5,50,-1.,1.);
		
	fpt_reco_pt_MC_num=new TH2F("fpt_reco_pt_MC_num","pt reco x pt MC;pt reco; pt MC",300,0.,30,300,0.,30);
	fpt_reco_pt_MC_den=new TH2F("fpt_reco_pt_MC_den","pt reco x pt MC;pt reco; pt MC",300,0.,30,300,0.,30);
	
	
	fCharge_n = new TH1F("fCharge_n","Inclusive Electrons (Negative Charge); p_{t} (GeV/c); Count",200,0,30);
	fCharge_p = new TH1F("fCharge_p","Inclusive Positrons (Positive Charge); p_{t} (GeV/c); Count",200,0,30);
	
	fECluster_pure= new TH1F("fECluster_pure", ";ECluster pure",2000,0,100);
	fECluster_not_exotic= new TH1F("fECluster_not_exotic", ";ECluster not exotic - function ",2000,0,100);
	
	fECluster_not_exotic1= new TH1F("fECluster_not_exotic1", ";ECluster not exotic Ncells > E/3+1",2000,0,100);

	fECluster_not_exotic2= new TH1F("fECluster_not_exotic2", ";ECluster not exotic 2",2000,0,100);
	fECluster_exotic= new TH1F("fECluster_exotic", ";ECluster exotic",2000,0,100);
	
	//not associated with tracks
	fNCluster_pure= new TH1F("fNCluster_pure", ";Number of clusters - pure",2000,-1,1999);
	fNCluster_pure_aod= new TH1F("fNCluster_pure_aod", ";Number of clusters - pure -aod",2000,-1,1999);
	fNCluster_ECluster= new TH2F("fNCluster_ECluster", ";Number of clusters vs. Energy of Cluster",2000,-1,1999, 4000, -1, 1999);
	
	fNcells_energy= new TH2F("fNcells_energy", "all clusters;Number of cells;Energy of Cluster",100,0,100, 2000, 0, 100);
	fNcells_energy_elec_selected= new TH2F("fNcells_energy_elec_selected", "clusters for electrons on TPC;Number of cells;Energy of Cluster",100,0,100, 2000, 0, 100);
	fNcells_energy_not_exotic= new TH2F("fNcells_energy_not_exotic", "not exotic cluster;Number of cells;Energy of Cluster ",100,0,100, 2000, 0, 100);
	
	if(fUseEMCal){
		
		if(!fIsAOD){
			fTime = new TH2D("fTime","Cells Cluster Time; p_{T} (GeV/c); Time (s)",300,0,30,1000,1e-8,1e-5);
			fTime2 = new TH2D("fTime2","Cells Cluster Time;  p_{T} (GeV/c); Time (s)",300,0,30,1000,1e-8,1e-5);
		}
		
		ftimingEle = new TH2D("ftimingEle","Cluster Time;  p_{T} (GeV/c); Time (s)",300,0,30,1000,1e-8,1e-5);
		ftimingEle2 = new TH2D("ftimingEle2","Cluster Time;  p_{T} (GeV/c); Time (s)",300,0,30,1000,1e-8,1e-5);
		
		fShowerShape_ha = new TH2F("fShowerShape_ha","Shower Shape hadrons;M02;M20",500,0,1.8,500,0,1.8);
		fShowerShape_ele = new TH2F("fShowerShape_ele","Shower Shape electrons;M02;M20",500,0,1.8,500,0,1.8);
		
		fShowerShapeM02_EoverP = new TH2F("fShowerShapeM02_EoverP","Shower Shape;M02;E/p",500,0,1.8,500,0,1.8);
		fShowerShapeM20_EoverP = new TH2F("fShowerShapeM20_EoverP","Shower Shape;M20;E/p",500,0,1.8,500,0,1.8);
		
	}
	
	fOutputList->Add(fTOF01);
	fOutputList->Add(fTOF02);
	fOutputList->Add(fTOF03);
	
	fOutputList->Add(fEtaPhi_num);
	fOutputList->Add(fEtaPhi_den);
	fOutputList->Add(fEtaPhi_data);
	
	fOutputList->Add(fpt_reco_pt_MC_num);
	fOutputList->Add(fpt_reco_pt_MC_den);
	
	
	fOutputList->Add(fPtElec_Inc);
	fOutputList->Add(fPtElec_ULS);
	fOutputList->Add(fPtElec_LS);
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
	
	if(fFillBackground){
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
	
	fOutputList->Add(fTPCnsigma_p_TPC);
	fOutputList->Add(fTPCnsigma_p_TPC_on_EMCal_acc);
	fOutputList->Add(fTPCnsigma_p_TPC_EoverP_cut);
	
	
	
	fOutputList->Add(fShowerShapeCut);
	
	fOutputList->Add(fCharge_n);
	fOutputList->Add(fCharge_p);
	
	fOutputList->Add(fECluster_pure);
	fOutputList->Add(fECluster_not_exotic);
	fOutputList->Add(fECluster_not_exotic1);
	fOutputList->Add(fECluster_not_exotic2);
	fOutputList->Add(fECluster_exotic);
	
	fOutputList->Add(fNCluster_pure);
	fOutputList->Add(fNCluster_pure_aod);

	fOutputList->Add(fNCluster_ECluster);
	fOutputList->Add(fNcells_energy);
	fOutputList->Add(fNcells_energy_elec_selected);
	fOutputList->Add(fNcells_energy_not_exotic);

	
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
	
	fEoverP_pt_true_electrons = new TH2F("fEoverP_pt_true_electrons",";p_{T} (GeV/c);E/p ",1000,0,30,2000,0,2);
	fOutputList->Add(fEoverP_pt_true_electrons);
	fEoverP_pt_true_electrons0 = new TH2F("fEoverP_pt_true_electrons0",";p_{T} (GeV/c);E/p ",1000,0,30,2000,0,2);
	fOutputList->Add(fEoverP_pt_true_electrons0);

	
	fEoverP_pt_true_hadrons = new TH2F("fEoverP_pt_true_hadrons",";p_{T} (GeV/c);E/p ",1000,0,30,2000,0,2);
	fOutputList->Add(fEoverP_pt_true_hadrons);
	fEoverP_pt_true_hadrons0 = new TH2F("fEoverP_pt_true_hadrons0",";p_{T} (GeV/c);E/p ",1000,0,30,2000,0,2);
	fOutputList->Add(fEoverP_pt_true_hadrons0);



	
	for(Int_t i = 0; i < 3; i++)
	{
		fEoverP_pt[i] = new TH2F(Form("fEoverP_pt%d",i),";p_{t} (GeV/c);E / p ",1000,0,30,2000,0,2);
		fTPC_p[i] = new TH2F(Form("fTPC_p%d",i),";pt (GeV/c);TPC dE/dx (a. u.)",1000,0.3,15,1000,-20,200);
		fTPCnsigma_p[i] = new TH2F(Form("fTPCnsigma_p%d",i),";p (GeV/c);TPC Electron N#sigma",1000,0.3,15,1000,-15,10);
		
		
		fECluster[i]= new TH1F(Form("fECluster%d",i), ";ECluster",2000, 0,100);
		fEtaPhi[i]= new TH2F(Form("fEtaPhi%d",i),"#eta x #phi Clusters;#phi;#eta",200,0.,5,50,-1.,1.);
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
	
		

	
	fTrack_Multi= new  TH1F("fTrack_Multi","fTrack_Multi",1000, 0,1000);
	
	for(Int_t i = 0; i < 4; i++)
	{
		fTPCNcls_pid[i]= new TH2F(Form("fTPCNcls_pid%d",i),"fTPCNcls_pid;NCls;NCls for PID",200,0,200,200,0,200);
		fOutputList->Add(fTPCNcls_pid[i]);
	}
	
		//pt bin
	Int_t fPtBin[7] = {1,2,4,6,8,10,15};
	
	
	Int_t fPtBin_trigger[11] = {1,2,4,6,8,10,12,14,16,18,20};
	fEoverP_tpc_p_trigger = new TH2F *[10];
	fEoverP_tpc_pt_trigger = new TH2F *[10];
	
	
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
		
		fOutputList->Add(fInvMass);
		fOutputList->Add(fInvMassBack);
		fOutputList->Add(fDCA);
		fOutputList->Add(fDCABack);
		fOutputList->Add(fOpAngle);
		fOutputList->Add(fOpAngleBack);
	}
	
	if(fFillBackground){
		
		fInvMass = new TH1F("fInvMass","",200,0,0.3);
		fInvMassBack = new TH1F("fInvMassBack","",200,0,0.3);
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
		fInvMass2 = new TH1F("fInvMass2","",200,0,0.3);
		fInvMassBack2 = new TH1F("fInvMassBack2","",200,0,0.3);
		fDCA2 = new TH1F("fDCA2","",200,0,1);
		fDCABack2 = new TH1F("fDCABack2","",200,0,1);
		fOpAngle2 = new TH1F("fOpAngle2","",200,0,0.5);
		fOpAngleBack2 = new TH1F("fOpAngleBack2","",200,0,0.5);
		
		fOutputList->Add(fInvMass2);
		fOutputList->Add(fInvMassBack2);
		fOutputList->Add(fDCA2);
		fOutputList->Add(fDCABack2);
		fOutputList->Add(fOpAngle2);
		fOutputList->Add(fOpAngleBack2);
		
	}
	
		//new histo for trigger data
	
	for(Int_t i = 0; i < 10; i++)
	{
			fEoverP_tpc_pt_trigger[i] = new TH2F(Form("fEoverP_tpc_pt_trigger%d",i),Form("%d < p_{t} < %d GeV/c;TPC Electron N#sigma; E/p ",fPtBin_trigger[i],fPtBin_trigger[i+1]),1000,-15,15,100,0,2);
			fEoverP_tpc_p_trigger[i] = new TH2F(Form("fEoverP_tpc_p_trigger%d",i),Form("%d < p < %d GeV/c;TPC Electron N#sigma; E/p ",fPtBin_trigger[i],fPtBin_trigger[i+1]),1000,-15,15,100,0,2);
			fOutputList->Add(fEoverP_tpc_pt_trigger[i]);
			fOutputList->Add(fEoverP_tpc_p_trigger[i]);

	}
	
	
	for(Int_t i = 0; i < 6; i++)
	{
		fEoverP_tpc[i] = new TH2F(Form("fEoverP_tpc%d",i),Form("%d < p_{t} < %d GeV/c;TPC Electron N#sigma;E / p ",fPtBin[i],fPtBin[i+1]),1000,-15,15,100,0,2);
		fTPC_pt[i] = new TH1F(Form("fTPC_pt%d",i),Form("%d < p_{t} < %d GeV/c;TPC Electron N#sigma;Count",fPtBin[i],fPtBin[i+1]),200,20,200);
		fTPCnsigma_pt[i] = new TH1F(Form("fTPCnsigma_pt%d",i),Form("%d < p_{t} < %d GeV/c;TPC Electron N#sigma;Count",fPtBin[i],fPtBin[i+1]),200,-15,10);
		
		fEta[i]=new TH1F(Form("fEta%d",i), Form("%d < p_{t} < %d GeV/c;#eta; counts",fPtBin[i],fPtBin[i+1]),100, -0.1,0.1);
		fPhi[i]=new TH1F(Form("fPhi%d",i),Form("%d < p_{t} < %d GeV/c;#phi; counts )",fPtBin[i],fPtBin[i+1]), 100, -0.1,0.1);
		fR[i]=new TH1F(Form("fR%d",i),Form("%d < p_{t} < %d GeV/c;R;counts )",fPtBin[i],fPtBin[i+1]), 100, -0.1,0.1);
		fR_EoverP[i]=new TH2F(Form("fR_EoverP%d",i),Form("%d < p_{t} < %d GeV/c;R;E / p ",fPtBin[i],fPtBin[i+1]),100, 0,0.1,1000,0,10);
		fNcells[i]=new TH1F(Form("fNcells%d",i), Form("%d < p_{t} < %d GeV/c;ncells;counts ",fPtBin[i],fPtBin[i+1]),100, 0, 30);
		fNcells_electrons[i]=new TH1F(Form("fNcells_electrons%d",i), Form("%d < p_{t} < %d GeV/c;ncells;counts ",fPtBin[i],fPtBin[i+1]),100, 0, 30);
		fNcells_hadrons[i]=new TH1F(Form("fNcells_hadrons%d",i), Form("%d < p_{t} < %d GeV/c;ncells;counts ",fPtBin[i],fPtBin[i+1]),100, 0, 30);
		fNcells_EoverP[i]=new TH2F(Form("fNcells_EoverP%d",i),Form("%d < p_{t} < %d GeV/c; Ncells; E / p ",fPtBin[i],fPtBin[i+1]),1000, 0,20,100,0,30);
		fM02_EoverP[i]= new TH2F(Form("fM02_EoverP%d",i),Form("%d < p_{t} < %d GeV/c; M02; E / p ",fPtBin[i],fPtBin[i+1]),1000,0,100,100,0,2);
		fM20_EoverP[i]= new TH2F(Form("fM20_EoverP%d",i),Form("%d < p_{t} < %d GeV/c; M20; E / p ",fPtBin[i],fPtBin[i+1]),1000,0,100,100,0,2);
		fEoverP_ptbins[i] = new TH1F(Form("fEoverP_ptbins%d",i),Form("%d < p_{t} < %d GeV/c;E / p ",fPtBin[i],fPtBin[i+1]),500,0,2);
		fECluster_ptbins[i]= new TH1F(Form("fECluster_ptbins%d",i), Form("%d < p_{t} < %d GeV/c;ECluster; Counts ",fPtBin[i],fPtBin[i+1]),2000, 0,100);
		fEoverP_wSSCut[i]=new TH1F(Form("fEoverP_wSSCut%d",i),Form("%d < p_{t} < %d GeV/c;E / p ; Counts",fPtBin[i],fPtBin[i+1]),500,0,2);
		fTPCnsigma_eta_electrons[i]=new TH2F(Form("fTPCnsigma_eta_electrons%d",i),Form("%d < p_{t} < %d GeV/c;TPC Electron N#sigma;Eta ",fPtBin[i],fPtBin[i+1]),1000,-15,15,100,-1,1);
		fTPCnsigma_eta_hadrons[i]=new TH2F(Form("fTPCnsigma_eta_hadrons%d",i),Form("%d < p_{t} < %d GeV/c;TPC Electron N#sigma;Eta ",fPtBin[i],fPtBin[i+1]),1000,-15,15,100,-1,1);
		
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
			fCEtaPhi_Inc[i] = new TH2F(Form("fCEtaPhi_Inc%d",i),Form("%d < p_{t} < %d GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
			fCEtaPhi_Inc_DiHadron[i] = new TH2F(Form("fCEtaPhi_Inc_DiHadron%d",i),Form("%d < p_{t} < %d GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
			
			fCEtaPhi_ULS[i] = new TH2F(Form("fCEtaPhi_ULS%d",i),Form("%d < p_{t} < %d GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
			fCEtaPhi_LS[i] = new TH2F(Form("fCEtaPhi_LS%d",i),Form("%d < p_{t} < %d GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
			fCEtaPhi_ULS_NoP[i] = new TH2F(Form("fCEtaPhi_ULS_NoP%d",i),Form("%d < p_{t} < %d GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
			fCEtaPhi_LS_NoP[i] = new TH2F(Form("fCEtaPhi_LS_NoP%d",i),Form("%d < p_{t} < %d GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
			
			fCEtaPhi_ULS_Weight[i] = new TH2F(Form("fCEtaPhi_ULS_Weight%d",i),Form("%d < p_{t} < %d GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
			fCEtaPhi_LS_Weight[i] = new TH2F(Form("fCEtaPhi_LS_Weight%d",i),Form("%d < p_{t} < %d GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
			fCEtaPhi_ULS_NoP_Weight[i] = new TH2F(Form("fCEtaPhi_ULS_NoP_Weight%d",i),Form("%d < p_{t} < %d GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
			fCEtaPhi_LS_NoP_Weight[i] = new TH2F(Form("fCEtaPhi_LS_NoP_Weight%d",i),Form("%d < p_{t} < %d GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
			
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
				fCEtaPhi_Inc_EM[i] = new TH2F(Form("fCEtaPhi_Inc_EM%d",i),Form("%d < p_{t} < %d GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
				
				fCEtaPhi_ULS_EM[i] = new TH2F(Form("fCEtaPhi_ULS_EM%d",i),Form("%d < p_{t} < %d GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
				fCEtaPhi_LS_EM[i] = new TH2F(Form("fCEtaPhi_LS_EM%d",i),Form("%d < p_{t} < %d GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
				
				fCEtaPhi_ULS_Weight_EM[i] = new TH2F(Form("fCEtaPhi_ULS_Weight_EM%d",i),Form("%d < p_{t} < %d GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
				fCEtaPhi_LS_Weight_EM[i] = new TH2F(Form("fCEtaPhi_LS_Weight_EM%d",i),Form("%d < p_{t} < %d GeV/c;#DeltaPhi (rad);#Delta#eta",fPtBin[i],fPtBin[i+1]),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-2,2);
				
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
	fEoverP_pt_hadrons= new TH2F("fEoverP_pt_hadrons","fEoverP_pt_hadrons",1000,0,30,2000,0,2);
	
	
	fOutputList->Add(fTPCnsigma_eta);
	fOutputList->Add(fTPCnsigma_phi);
	
	fOutputList->Add(fNcells_pt);
	fOutputList->Add(fEoverP_pt_pions);
	
	fOutputList->Add(ftpc_p_EoverPcut);
	fOutputList->Add(fnsigma_p_EoverPcut);
	
	fOutputList->Add(fEoverP_pt_pions2);
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
		fpT_m_electron= new TH2F("fpT_m_electron","fpT_m_electron",300,0,30,300,0,30);
		fpT_gm_electron= new TH2F("fpT_gm_electron","fpT_gm_electron",300,0,30,300,0,30);
		
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
		

		
		fPtMCpi0 = new TH1F("fPtMCpi0",";p_{t} (GeV/c);Count",200,0,30);
		fPtMCeta = new TH1F("fPtMCeta",";p_{T} (GeV/c);Count",200,0,30);
		fPtMCpi02 = new TH1F("fPtMCpi02",";p_{t} (GeV/c);Count",200,0,30);
		fPtMCeta2 = new TH1F("fPtMCeta2",";p_{T} (GeV/c);Count",200,0,30);
		fPtMCpi03 = new TH1F("fPtMCpi03",";p_{t} (GeV/c);Count",200,0,30);
		fPtMCeta3 = new TH1F("fPtMCeta3",";p_{T} (GeV/c);Count",200,0,30);
		
		fPtMC_EMCal_All= new TH1F("fPtMC_EMCal_All",";p_{t} (GeV/c);Count",200,0,40);
		fPtMC_EMCal_Selected= new TH1F("fPtMC_EMCal_Selected",";p_{t} (GeV/c);Count",200,0,40);
		fPtMC_TPC_All= new TH1F("fPtMC_TPC_All",";p_{T} (GeV/c);Count",200,0,40);
		fPtMC_TPC_Selected = new TH1F("fPtMC_TPC_Selected",";p_{T} (GeV/c);Count",200,0,40);
		
		fPt_track_match_den = new TH1F("fPt_track_match_den",";p_{T} (GeV/c);Count",200,0,40);
		fPt_track_match_num = new TH1F("fPt_track_match_num",";p_{T} (GeV/c);Count",200,0,40);
		
		fPtMCWithLabel = new TH1F("fPtMCWithLabel",";p_{t} (GeV/c);Count",200,0,40);
		fPtMCWithoutLabel = new TH1F("fPtMCWithoutLabel",";p_{t} (GeV/c);Count",200,0,40);
		fPtIsPhysicaPrimary = new TH1F("fPtIsPhysicaPrimary",";p_{t} (GeV/c);Count",200,0,40);
		
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
	
	fCentralityHist = new TH1F("fCentralityHist",";Centrality (%); Count",1000000,0,100);
	fCentralityHistPass = new TH1F("fCentralityHistPass",";Centrality (%); Count",1000000,0,100);
	fOutputList->Add(fCentralityHist);
	fOutputList->Add(fCentralityHistPass);
	
	//______________________________________________________________________
	//Mixed event analysis
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
void AliAnalysisTaskEMCalHFEpA::UserExec(Option_t *) 
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
	if(fIsAOD)
	{
		const AliAODVertex* trkVtx = fAOD->GetPrimaryVertex();
			
		//events without vertex
		if(!trkVtx) fNevent2->Fill(0);
		
		if(!trkVtx || trkVtx->GetNContributors()<=0) return;
		TString vtxTtl = trkVtx->GetTitle();
		if(!vtxTtl.Contains("VertexerTracks")) return;
			//Float_t zvtx = trkVtx->GetZ();
		Float_t zvtx = -100;
		zvtx=trkVtx->GetZ();
		fZvtx = zvtx;
		
		fVtxZ_new1->Fill(fZvtx);
		
		const AliAODVertex* spdVtx = fAOD->GetPrimaryVertexSPD();
		if(spdVtx->GetNContributors()<=0) return;
		TString vtxTyp = spdVtx->GetTitle();
		Double_t cov[6]={0};
		spdVtx->GetCovarianceMatrix(cov);
		Double_t zRes = TMath::Sqrt(cov[5]);
		
		fzRes1->Fill(zRes);
		if(vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return;
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
	//new track loop to select events
	//track pt cut (at least 2)
	/*
	if(fUseTrigger){
		if(fIsAOD){
			double fTrackMulti=0;
			for(Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++) 
			{
				AliVParticle* Vtrack = fVevent->GetTrack(iTracks);
				if (!Vtrack) continue;
				
				AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
					//AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);
		
				if((track->Pt())<0.2 || (track->Pt())>1000.0) continue;
			   	else fTrackMulti=fTrackMulti+1;
		
			}
				//Only take event if track multiplicity is bigger than 2.
			if(fTrackMulti<2) return;
		}
	}
	fNevent->Fill(13);	
     */
	//______________________________________________________________________
	//Using more cuts than I have beeing using
	//eta cut and primary (at least 2)
    /*
	if(fUseTrigger){
		if(fIsAOD){
			double fTrackMulti2=0;
			for(Int_t i = 0; i < fVevent->GetNumberOfTracks(); i++) 
			{
				AliVParticle* Vtrack2 = fVevent->GetTrack(i);
				if (!Vtrack2) continue;
				
				
				AliVTrack *track_new = dynamic_cast<AliVTrack*>(Vtrack2);
				AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(Vtrack2);
				
				
				if(aodtrack)
				{
					
					
				    if(TMath::Abs(track_new->Eta())> 0.9) continue;
					if (aodtrack->GetType()!= AliAODTrack::kPrimary) continue ;
				    else fTrackMulti2=fTrackMulti2+1;
				}
			}
				//Only take event if track multiplicity is bigger than 2.
			if(fTrackMulti2<2) return;

			
		}
	}
	fNevent->Fill(14);	
//______________________________________________________________________
//Using more cuts than I have beeing using
//hybrid (at least2)
	if(fUseTrigger){
		if(fIsAOD){
			double fTrackMulti3=0;
			for(Int_t i = 0; i < fVevent->GetNumberOfTracks(); i++) 
			{
				AliVParticle* Vtrack3 = fVevent->GetTrack(i);
				if (!Vtrack3) continue;
								
					//AliVTrack *track_new = dynamic_cast<AliVTrack*>(Vtrack3);
				AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(Vtrack3);
				
				
				if(aodtrack)
				{
					
					if (!aodtrack->IsHybridGlobalConstrainedGlobal()) continue ;
						//another option if I don't want to use hybrid
						//if ( aodtrack->TestFilterBit(128)==kFALSE) continue ;
				    else fTrackMulti3=fTrackMulti3+1;
				}
			}
			//Only take event if track multiplicity is bigger than 2.
			if(fTrackMulti3<2) return;

		}
	}
	fNevent->Fill(15);	
//______________________________________________________________________
	
	
	if(fUseTrigger){
		if(fIsAOD){
			double fTrackMulti4=0;
			for(Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++) 
			{
				AliVParticle* Vtrack4 = fVevent->GetTrack(iTracks);
				if (!Vtrack4) continue;
				
				
					//AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack4);
				AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack4);
				
				if(!atrack->TestFilterBit(768)) continue;
				if(!atrack->IsHybridGlobalConstrainedGlobal()) continue ;
				
				
				else fTrackMulti4=fTrackMulti4+1;
				
			}
			//Only take event if track multiplicity is bigger than 2.
			if(fTrackMulti4<2) return;
			fTrack_Multi->Fill(fTrackMulti4);
		}
	}
	fNevent->Fill(16);	
     */
//______________________________________________________________________
	
	
//______________________________________________________________________
//Centrality Selection
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
		
		if(fEstimator==1) centrality = fCentrality->GetCentralityPercentile("ZDC");
		else centrality = fCentrality->GetCentralityPercentile("V0A");
		
			//cout << "Centrality = " << centrality << " %" << endl;
		
		fCentralityHist->Fill(centrality);
		
		if(centrality<fCentralityMin || centrality>fCentralityMax) return;
		
		fCentralityHistPass->Fill(centrality);
	}
	//______________________________________________________________________
	
	
	fNevent->Fill(17);
	
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
				/*
				exotic   = fEMCALRecoUtils->IsExoticCluster(clust, (AliVCaloCells*)fAOD->GetEMCALCells(), bc);
				if(exotic == kFALSE){ 
					fECluster_not_exotic->Fill(clust->E());
					fNcells_energy_not_exotic->Fill(clust->GetNCells(),clust->E());
				}
				*/
				
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
			/*
			//______________________________________________________________________
			//Trying to remove events with bad cells and find patches
			//First, I will try to count them
			//______________________________________________________________________

			if(clust && clust->IsEMCAL())
			{
				Bool_t badchannel = ContainsBadChannel("EMCAL", clust->GetCellsAbsId(),clust->GetNCells() );
				printf("Contém bad channel? %d ", badchannel);
				if(badchannel)fNevent2->Fill(0); 
				
				//trying to find patches
				TArrayI patches_found=GetTriggerPatches(IsEventEMCALL0, IsEventEMCALL1);
				printf("N patches %d, first %d, last %d\n",patches_found.GetSize(),  patches_found.At(0), patches_found.At(patches_found.GetSize()-1));

			}
			
			//end of bad cells
			//______________________________________________________________________
			*/
			
		}
	}

	
	
	fNevent->Fill(24); //events with cluster
		
	
	fVtxZ_new4->Fill(fZvtx);
	

	Int_t ClsNo = -999;
	if(!fIsAOD) ClsNo = fESD->GetNumberOfCaloClusters(); 
	else ClsNo = fAOD->GetNumberOfCaloClusters(); 
	
	//______________________________________________________________________
	
	///_____________________________________________________________________
	///Track loop
	Int_t NTracks=0;
	
    if(!fUseTender) NTracks=fVevent->GetNumberOfTracks();
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//To use tender
	if(fUseTender){
		TClonesArray  *fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AODFilterTracks"));
		NTracks = fTracks_tender->GetEntries();
		TClonesArray  *fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));
		ClsNo = fCaloClusters_tender->GetEntries();
		
		
	}
	
	
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
		
		fTPC_p[0]->Fill(fPt,fTPCsignal);
		fTPCnsigma_p[0]->Fill(fP,fTPCnSigma);
		
		
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
            printf("\n\n cluster label of track %d is %d  - from tender \n\n", iTracks, fClus->GetLabel());
            
            
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
		
		fTPC_p[1]->Fill(fPt,fTPCsignal);
		fTPCnsigma_p[1]->Fill(fP,fTPCnSigma);
		Double_t fPtBin[7] = {1,2,4,6,8,10,15};
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
						//true hadrons E/p shape before cuts
						if(fIsMC && fIsAOD && track->GetLabel()>=0){
							fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
							Int_t pdg = fMCparticle->GetPdgCode();
							
							if( TMath::Abs(pdg) != 11){
								fEoverP_pt_true_hadrons0->Fill(fPt,(fClus->E() / fP));
							}
						}
						
						
						//true electrons E/p shape before cuts
						if(fIsMC && fIsAOD && track->GetLabel()>=0){
							fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
							Int_t pdg = fMCparticle->GetPdgCode();
							
							if( TMath::Abs(pdg) == 11){
								fEoverP_pt_true_electrons0->Fill(fPt,(fClus->E() / fP));
							}
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
		
			
			//new way to calculate TPCnsigma distribution: TPCnsigma in function of p, with/without E/p cut 
			fTPCnsigma_p_TPC->Fill(fP, fTPCnSigma);
			if(fEMCflag){
			
				fTPCnsigma_p_TPC_on_EMCal_acc->Fill(fP, fTPCnSigma);
			
				if((fClus->E() / fP) >= fEoverPCutMin && (fClus->E() / fP) <= fEoverPCutMax){
					fTPCnsigma_p_TPC_EoverP_cut->Fill(fP, fTPCnSigma);
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
				Background(track, iTracks, Vtrack, kTRUE); //IsTPConly=kTRUE
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
					ftimingEle->Fill(fPt,emctof);
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
								//in order to check if there are exotic cluster in this selected cluster (27 may 2014)
								fNcells_energy_elec_selected->Fill(ncells,Energy);
								
									//true electrons E/p shape
								if(fIsMC && fIsAOD && track->GetLabel()>=0){
									fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
									Int_t pdg = fMCparticle->GetPdgCode();
									
									if( TMath::Abs(pdg) == 11){
										fEoverP_pt_true_electrons->Fill(fPt,(fClus->E() / fP));
									}
								}
								
								
								
							}
							
						}
						if(!fUseShowerShapeCut){
							fEoverP_pt[2]->Fill(fPt,(fClus->E() / fP));
							fShowerShapeCut->Fill(M02,M20);
							fNcells_energy_elec_selected->Fill(ncells,Energy);
							
							//after cut -> Are they really electrons ? E/p shape
							if(fIsMC && fIsAOD && track->GetLabel()>=0){
								fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
								Int_t pdg = fMCparticle->GetPdgCode();
								
								if( TMath::Abs(pdg) == 11){
									fEoverP_pt_true_electrons->Fill(fPt,(fClus->E() / fP));
								}
							}

							
						}
						if(fUseEMCal) fShowerShape_ele->Fill(M02,M20);
						
						//for shower shape cut studies - now with TPC PID Cut
						if(fUseEMCal){
							fShowerShapeM02_EoverP->Fill(M02,EoverP);
							fShowerShapeM20_EoverP->Fill(M20,EoverP);
						}
						
					}
					
					//----------------------------------------------------------------------------------------
					
					
					
					// for Eta Phi distribution
					fClus->GetPosition(pos3);
					TVector3 vpos(pos3[0],pos3[1],pos3[2]);
					Double_t cphi = vpos.Phi();
					Double_t ceta = vpos.Eta();
					fEtaPhi[2]->Fill(cphi,ceta);
					
					
					
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
						
						
						if(fUseEMCal)
						{
							if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax ){
								fPtElec_Inc->Fill(fPt);
								//eta phi distribution for data
								fEtaPhi_data->Fill(track->Phi(),track->Eta());
							}
							
							//Eta cut for background
							if(fFillBackground){
								fEtad[2]->Fill(track->Eta());
								
								//background for triggered data: trigger electron must have same cuts on shower shape  06/Jan/2014
								if(fUseShowerShapeCut){
									if(M02 >= fM02CutMin && M02<=fM02CutMax && M20>=fM20CutMin && M20<=fM20CutMax){
										if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax){
											Background(track, iTracks, Vtrack, kFALSE);
										}
									}
								}
								else{
									if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax){
										Background(track, iTracks, Vtrack, kFALSE);
									}
								}
								
							}
							
							double emctof2 = fClus->GetTOF();
							ftimingEle2->Fill(fPt,emctof2);
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
		fPool = fPoolMgr->GetEventPool(fCentrality->GetCentralityPercentile("V0A"), fZvtx); // Get the buffer associated with the current centrality and z-vtx
		
		if(!fPool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fCentrality->GetCentralityPercentile("V0A"), fZvtx));
		
		fPool->UpdatePool(SelectedHadrons()); // fill the tracks in the event buffer. The ownership is handed over to the event mixing class. We are not allowed to delete tracksClone anymore!
		
		
	}
	
		//__________________________________________________________________
	
	delete fListOfmotherkink;
	PostData(1, fOutputList);
}      

	//______________________________________________________________________
void AliAnalysisTaskEMCalHFEpA::Terminate(Option_t *) 
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
Bool_t AliAnalysisTaskEMCalHFEpA::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
		//Check single track cuts for a given cut step
		//Note this function is called inside the UserExec function
	const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
	if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
	return kTRUE;
}


	//______________________________________________________________________


void AliAnalysisTaskEMCalHFEpA::Background(AliVTrack *track, Int_t trackIndex, AliVParticle *vtrack, Bool_t IsTPConly)
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
		
		if(!IsTPConly){
			fNonHFE->SetHistAngleBack(fOpAngleBack);
			fNonHFE->SetHistAngle(fOpAngle);
			fNonHFE->SetHistDCABack(fDCABack);
			fNonHFE->SetHistDCA(fDCA);
			fNonHFE->SetHistMassBack(fInvMassBack);
			fNonHFE->SetHistMass(fInvMass);
		}
		if(IsTPConly){
			fNonHFE->SetHistAngleBack(fOpAngleBack2);
			fNonHFE->SetHistAngle(fOpAngle2);
			fNonHFE->SetHistDCABack(fDCABack2);
			fNonHFE->SetHistDCA(fDCA2);
			fNonHFE->SetHistMassBack(fInvMassBack2);
			fNonHFE->SetHistMass(fInvMass2);
		}
		
		fNonHFE->FindNonHFE(trackIndex,vtrack,fVevent);
		
			//index of track selected as partner
		Int_t *fUlsPartner = fNonHFE->GetPartnersULS();
		
			//Electron Information
		Double_t fPhiE = -999;
		Double_t fEtaE = -999;
		Double_t fPtE = -999;
		fPhiE = track->Phi();
		fEtaE = track->Eta();
		fPtE = track->Pt();
		
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
void AliAnalysisTaskEMCalHFEpA::ElectronHadronCorrelation(AliVTrack *track, Int_t trackIndex, AliVParticle *vtrack)
{
	
	///_________________________________________________________________
	///MC analysis
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
	        
	        if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==22 || TMath::Abs(fMCparticleMother->GetPdgCode())==111 || TMath::Abs(fMCparticleMother->GetPdgCode())==221))
	        {
					//Is Background
				fPtBackgroundBeforeReco->Fill(track->Pt());
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
	if(fEventMixingFlag)
	{
		fPool = fPoolMgr->GetEventPool(fCentrality->GetCentralityPercentile("V0A"), fZvtx); // Get the buffer associated with the current centrality and z-vtx
		
		if(!fPool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f",fCentrality->GetCentralityPercentile("V0A"), fZvtx));
		
		if(fPool->GetCurrentNEvents() >= 5) // start mixing when 5 events are in the buffer
		{
			fPoolNevents->Fill(fPool->GetCurrentNEvents());
			
			for (Int_t jMix = 0; jMix < fPool->GetCurrentNEvents(); jMix++)  // mix with each event in the buffer
			{
				TObjArray* bgTracks = fPool->GetEvent(jMix);
				
				for (Int_t kMix = 0; kMix < bgTracks->GetEntriesFast(); kMix++)  // mix with each track in the event
				{
					const AliEHCParticle* MixedTrack(dynamic_cast<AliEHCParticle*>(bgTracks->At(kMix)));
					if (NULL == MixedTrack) continue;
					
					fPhiH = MixedTrack->Phi();
					fEtaH = MixedTrack->Eta();
					fPtH = MixedTrack->Pt();
					
					if(fPtH<fAssHadronPtMin || fPtH>fAssHadronPtMax) continue;
					
					fDphi = fPhiE - fPhiH;
					
					if (fDphi > 3*pi/2) fDphi = fDphi - 2*pi;
					if (fDphi < -pi/2)  fDphi = fDphi + 2*pi;
					
					fDeta = fEtaE - fEtaH;
					
					Double_t fPtBin[7] = {1,2,4,6,8,10,15};
					
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
		
		fPhiH = track2->Phi();
		fEtaH = track2->Eta();
		fPtH = track2->Pt();
		
		if(fPtH<fAssHadronPtMin || fPtH>fAssHadronPtMax) continue;
		
		fDphi = fPhiE - fPhiH;
		
		if (fDphi > 3*pi/2) fDphi = fDphi - 2*pi;
		if (fDphi < -pi/2)  fDphi = fDphi + 2*pi;
		
		fDeta = fEtaE - fEtaH;
		
		Double_t fPtBin[7] = {1,2,4,6,8,10,15};
		
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
				fCEtaPhi_Inc[i]->Fill(fDphi,fDeta);
				
				if(fNonHFE->IsULS()) fCEtaPhi_ULS[i]->Fill(fDphi,fDeta);
				if(fNonHFE->IsLS()) fCEtaPhi_LS[i]->Fill(fDphi,fDeta);
				if(fNonHFE->IsULS() && !fUlsIsPartner) fCEtaPhi_ULS_NoP[i]->Fill(fDphi,fDeta);
				if(fNonHFE->IsLS() && !fLsIsPartner) fCEtaPhi_LS_NoP[i]->Fill(fDphi,fDeta);
				
				if(fNonHFE->IsULS()) fCEtaPhi_ULS_Weight[i]->Fill(fDphi,fDeta,fNonHFE->GetNULS());
				if(fNonHFE->IsLS()) fCEtaPhi_LS_Weight[i]->Fill(fDphi,fDeta,fNonHFE->GetNLS());
				if(fNonHFE->IsULS() && !fUlsIsPartner) fCEtaPhi_ULS_NoP_Weight[i]->Fill(fDphi,fDeta,fNonHFE->GetNULS());
				if(fNonHFE->IsLS() && !fLsIsPartner) fCEtaPhi_LS_NoP_Weight[i]->Fill(fDphi,fDeta,fNonHFE->GetNLS());
		    }
		}
	}
}

	//____________________________________________________________________________________________________________
	//Create a TObjArray with selected hadrons, for the mixed event analysis
TObjArray* AliAnalysisTaskEMCalHFEpA::SelectedHadrons()
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
		
		fTracksClone->Add(new AliEHCParticle(track2->Eta(), track2->Phi(), track2->Pt()));
	}
	return fTracksClone;
}
	//____________________________________________________________________________________________________________

	//______________________________________________________________________
void AliAnalysisTaskEMCalHFEpA::DiHadronCorrelation(AliVTrack *track, Int_t trackIndex)
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
		
		fPhiH = track2->Phi();
		fEtaH = track2->Eta();
		fPtH = track2->Pt();
		
		if(fPtH<fAssHadronPtMin || fPtH>fAssHadronPtMax) continue;
		
		fDphi = fPhiE - fPhiH;
		
		if (fDphi > 3*pi/2) fDphi = fDphi - 2*pi;
		if (fDphi < -pi/2)  fDphi = fDphi + 2*pi;
		
		fDeta = fEtaE - fEtaH;
		
		Double_t fPtBin[7] = {1,2,4,6,8,10,15};
		
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
Bool_t AliAnalysisTaskEMCalHFEpA::FindMother(Int_t mcIndex)
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
/*
Bool_t AliAnalysisTaskEMCalHFEpA::ContainsBadChannel(TString calorimeter,UShort_t* cellList, Int_t nCells)
{
		// Check that in the cluster cells, there is no bad channel of those stored 
		// in fEMCALBadChannelMap 
		// adapted from AliCalorimeterUtils
	
		//if (!fRemoveBadChannels) return kFALSE;
		//printf("fEMCALBadChannelMap %p, fPHOSBadChannelMap %p \n",fEMCALBadChannelMap,fPHOSBadChannelMap);
	if( !fEMCALRecoUtils->GetEMCALChannelStatusMap(0)) return kFALSE;
		
		//Int_t icol = -1;
		//Int_t irow = -1;
		//Int_t imod = -1;
	for(Int_t iCell = 0; iCell<nCells; iCell++){
		
			//Get the column and row
		if(calorimeter == "EMCAL"){
			return fEMCALRecoUtils->ClusterContainsBadChannel((AliEMCALGeometry*)fEMCALGeo,cellList,nCells);
		}
		else return kFALSE;
		
	}// cell cluster loop
	
	return kFALSE;
	
}
*/
/*

//________________________________________________________________________________
TArrayI AliAnalysisTaskEMCalHFEpA::GetTriggerPatches(Bool_t IsEventEMCALL0, Bool_t IsEventEMCALL1)
{
	// Select the patches that triggered
	// Depend on L0 or L1
	
	// some variables
	//Int_t  trigtimes[30], globCol, globRow,ntimes, i;
	Int_t   globCol, globRow;

	Int_t  absId  = -1; //[100];
	Int_t  nPatch = 0;
	
	TArrayI patches(0);
	
		// get object pointer
	AliVCaloTrigger *caloTrigger = InputEvent()->GetCaloTrigger( "EMCAL" );
	
		// class is not empty
	if( caloTrigger->GetEntries() > 0 )
	{
			// must reset before usage, or the class will fail
		caloTrigger->Reset();
		
			// go throuth the trigger channels
		while( caloTrigger->Next() )
		{
				// get position in global 2x2 tower coordinates
			caloTrigger->GetPosition( globCol, globRow );
			
			//L0
			if(IsEventEMCALL0)
			{
			    //not implemented
			}
				
					
			else if(IsEventEMCALL1) // L1
			{
					Int_t bit = 0;
					caloTrigger->GetTriggerBits(bit);
					
					Bool_t isEGA = ((bit >> fBitEGA) & 0x1);
					//Bool_t isEJE = ((bit >> fBitEJE) & 0x1) && IsEventEMCALL1Jet  () ;
					
					if(!isEGA) continue;
					
					Int_t patchsize = -1;
					if(isEGA) patchsize =  2;
					//else if (isEJE) patchsize = 16;
					
					// add 2x2 (EGA) or 16x16 (EJE) patches
					for(Int_t irow=0; irow < patchsize; irow++)
					{
						for(Int_t icol=0; icol < patchsize; icol++)
						{
							GetCaloUtils()->GetEMCALGeometry()->GetAbsFastORIndexFromPositionInEMCAL(globCol+icol,globRow+irow, absId);
							//printf("pass the time cut globCol %d, globRow %d absId %d\n",globCol,globRow, absIDTrig[nTrig]);
							patches.Set(nPatch+1); // Set size of this array to nPatch+1 ints.
							patches.AddAt(absId,nPatch++); //Add Int_t absId at position nPatch++
						}
					}
					
			} // L1
			
		} // trigger iterator
	} // go thorough triggers
	
	printf("N patches %d, test %d,first %d, last %d\n",patches.GetSize(), nPatch, patches.At(0), patches.At(patches.GetSize()-1));
	
	return patches;
}
 */
Double_t AliAnalysisTaskEMCalHFEpA::CalculateWeight(Int_t pdg_particle, Double_t x)
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
			if(x>= 0.000000 &&  x < 0.150000 ) weight=0.000135325731788;
			if(x>= 0.150000 &&  x < 0.300000 ) weight=0.000046086633955;
			if(x>= 0.300000 &&  x < 0.450000 ) weight=0.000039057281071;
			if(x>= 0.450000 &&  x < 0.600000 ) weight=0.000034760305945;
			if(x>= 0.600000 &&  x < 0.750000 ) weight=0.000031169797838;
			if(x>= 0.750000 &&  x < 0.900000 ) weight=0.000028658947305;
			if(x>= 0.900000 &&  x < 1.050000 ) weight=0.000027486443396;
			if(x>= 1.050000 &&  x < 1.200000 ) weight=0.000027712239624;
			if(x>= 1.200000 &&  x < 1.350000 ) weight=0.000029203559540;
			if(x>= 1.350000 &&  x < 1.500000 ) weight=0.000031995388848;
			if(x>= 1.500000 &&  x < 1.650000 ) weight=0.000036184304345;
			if(x>= 1.650000 &&  x < 1.800000 ) weight=0.000041880006402;
			if(x>= 1.800000 &&  x < 1.950000 ) weight=0.000049145632185;
			if(x>= 1.950000 &&  x < 2.100000 ) weight=0.000058350740019;
			if(x>= 2.100000 &&  x < 2.250000 ) weight=0.000070310111501;
			if(x>= 2.250000 &&  x < 2.400000 ) weight=0.000084979270329;
			if(x>= 2.400000 &&  x < 2.550000 ) weight=0.000102835462742;
			if(x>= 2.550000 &&  x < 2.700000 ) weight=0.000125717853414;
			if(x>= 2.700000 &&  x < 2.850000 ) weight=0.000151686709723;
			if(x>= 2.850000 &&  x < 3.000000 ) weight=0.000183539426446;
			if(x>= 3.000000 &&  x < 3.150000 ) weight=0.000220301846542;
			if(x>= 3.150000 &&  x < 3.300000 ) weight=0.000265463496589;
			if(x>= 3.300000 &&  x < 3.450000 ) weight=0.000317125871378;
			if(x>= 3.450000 &&  x < 3.600000 ) weight=0.000375125289831;
			if(x>= 3.600000 &&  x < 3.750000 ) weight=0.000443145614527;
			if(x>= 3.750000 &&  x < 3.900000 ) weight=0.000521789484430;
			if(x>= 3.900000 &&  x < 4.050000 ) weight=0.000609021217094;
			if(x>= 4.050000 &&  x < 4.200000 ) weight=0.000712204716847;
			if(x>= 4.200000 &&  x < 4.350000 ) weight=0.000823467616904;
			if(x>= 4.350000 &&  x < 4.500000 ) weight=0.000950122993528;
			if(x>= 4.500000 &&  x < 4.650000 ) weight=0.001097646187038;
			if(x>= 4.650000 &&  x < 4.800000 ) weight=0.001257914140236;
			if(x>= 4.800000 &&  x < 4.950000 ) weight=0.001433662433969;
			if(x>= 4.950000 &&  x < 5.100000 ) weight=0.001638196225347;
			if(x>= 5.100000 &&  x < 5.250000 ) weight=0.001851181851766;
			if(x>= 5.250000 &&  x < 5.400000 ) weight=0.002090388287398;
			if(x>= 5.400000 &&  x < 5.550000 ) weight=0.002359067258019;
			if(x>= 5.550000 &&  x < 5.700000 ) weight=0.002651132162420;
			if(x>= 5.700000 &&  x < 5.850000 ) weight=0.002975324648838;
			if(x>= 5.850000 &&  x < 6.000000 ) weight=0.003322013569821;
			if(x>= 6.000000 &&  x < 6.150000 ) weight=0.003691360208955;
			if(x>= 6.150000 &&  x < 6.300000 ) weight=0.004121404405681;
			if(x>= 6.300000 &&  x < 6.450000 ) weight=0.004561083532246;
			if(x>= 6.450000 &&  x < 6.600000 ) weight=0.005039509700197;
			if(x>= 6.600000 &&  x < 6.750000 ) weight=0.005552403734677;
			if(x>= 6.750000 &&  x < 6.900000 ) weight=0.006143291566070;
			if(x>= 6.900000 &&  x < 7.050000 ) weight=0.006742363968985;
			if(x>= 7.050000 &&  x < 7.200000 ) weight=0.007412947679000;
			if(x>= 7.200000 &&  x < 7.350000 ) weight=0.008084803811283;
			if(x>= 7.350000 &&  x < 7.500000 ) weight=0.008806078015679;
			if(x>= 7.500000 &&  x < 7.650000 ) weight=0.009631528972802;
			if(x>= 7.650000 &&  x < 7.800000 ) weight=0.010452967896416;
			if(x>= 7.800000 &&  x < 7.950000 ) weight=0.011377835648442;
			if(x>= 7.950000 &&  x < 8.100000 ) weight=0.012383279982470;
			if(x>= 8.100000 &&  x < 8.250000 ) weight=0.013379860387337;
			if(x>= 8.250000 &&  x < 8.400000 ) weight=0.014502732126184;
			if(x>= 8.400000 &&  x < 8.550000 ) weight=0.015651745450755;
			if(x>= 8.550000 &&  x < 8.700000 ) weight=0.016912137864860;
			if(x>= 8.700000 &&  x < 8.850000 ) weight=0.018179537211066;
			if(x>= 8.850000 &&  x < 9.000000 ) weight=0.019519322440302;
			if(x>= 9.000000 &&  x < 9.150000 ) weight=0.021086523924623;
			if(x>= 9.150000 &&  x < 9.300000 ) weight=0.022577769553812;
			if(x>= 9.300000 &&  x < 9.450000 ) weight=0.024161713220839;
			if(x>= 9.450000 &&  x < 9.600000 ) weight=0.025923912622410;
			if(x>= 9.600000 &&  x < 9.750000 ) weight=0.027505228569229;
			if(x>= 9.750000 &&  x < 9.900000 ) weight=0.029516659109336;
			if(x>= 9.900000 &&  x < 10.050000 ) weight=0.031570966309737;
			if(x>= 10.050000 &&  x < 10.200000 ) weight=0.033642267228791;
			if(x>= 10.200000 &&  x < 10.350000 ) weight=0.036041256148459;
			if(x>= 10.350000 &&  x < 10.500000 ) weight=0.038058932654582;
			if(x>= 10.500000 &&  x < 10.650000 ) weight=0.040566715148747;
			if(x>= 10.650000 &&  x < 10.800000 ) weight=0.042936804799602;
			if(x>= 10.800000 &&  x < 10.950000 ) weight=0.045464318029507;
			if(x>= 10.950000 &&  x < 11.100000 ) weight=0.048597779828317;
			if(x>= 11.100000 &&  x < 11.250000 ) weight=0.051354140021148;
			if(x>= 11.250000 &&  x < 11.400000 ) weight=0.054404316910916;
			if(x>= 11.400000 &&  x < 11.550000 ) weight=0.057213971371882;
			if(x>= 11.550000 &&  x < 11.700000 ) weight=0.060338442884031;
			if(x>= 11.700000 &&  x < 11.850000 ) weight=0.063940683079276;
			if(x>= 11.850000 &&  x < 12.000000 ) weight=0.067652488744517;
			if(x>= 12.000000 &&  x < 12.150000 ) weight=0.071297319739872;
			if(x>= 12.150000 &&  x < 12.300000 ) weight=0.074674823804324;
			if(x>= 12.300000 &&  x < 12.450000 ) weight=0.078754306348516;
			if(x>= 12.450000 &&  x < 12.600000 ) weight=0.083919533079495;
			if(x>= 12.600000 &&  x < 12.750000 ) weight=0.087741334859234;
			if(x>= 12.750000 &&  x < 12.900000 ) weight=0.092579682905166;
			if(x>= 12.900000 &&  x < 13.050000 ) weight=0.096872992585156;
			if(x>= 13.050000 &&  x < 13.200000 ) weight=0.102224349668488;
			if(x>= 13.200000 &&  x < 13.350000 ) weight=0.107979732741213;
			if(x>= 13.350000 &&  x < 13.500000 ) weight=0.111859278029459;
			if(x>= 13.500000 &&  x < 13.650000 ) weight=0.117754103597258;
			if(x>= 13.650000 &&  x < 13.800000 ) weight=0.124068068275799;
			if(x>= 13.800000 &&  x < 13.950000 ) weight=0.130317379729910;
			if(x>= 13.950000 &&  x < 14.100000 ) weight=0.136395722603921;
			if(x>= 14.100000 &&  x < 14.250000 ) weight=0.143695297272119;
			if(x>= 14.250000 &&  x < 14.400000 ) weight=0.150677344409986;
			if(x>= 14.400000 &&  x < 14.550000 ) weight=0.156193309002267;
			if(x>= 14.550000 &&  x < 14.700000 ) weight=0.164285190782173;
			if(x>= 14.700000 &&  x < 14.850000 ) weight=0.172158434103904;
			if(x>= 14.850000 &&  x < 15.000000 ) weight=0.180446007360110;
			if(x>= 15.000000 &&  x < 15.150000 ) weight=0.190459747433818;
			if(x>= 15.150000 &&  x < 15.300000 ) weight=0.196724851908457;
			if(x>= 15.300000 &&  x < 15.450000 ) weight=0.206428930550014;
			if(x>= 15.450000 &&  x < 15.600000 ) weight=0.216753360478597;
			if(x>= 15.600000 &&  x < 15.750000 ) weight=0.227117030529296;
			if(x>= 15.750000 &&  x < 15.900000 ) weight=0.236834809412698;
			if(x>= 15.900000 &&  x < 16.050000 ) weight=0.248197956916523;
			if(x>= 16.050000 &&  x < 16.200000 ) weight=0.259425526397730;
			if(x>= 16.200000 &&  x < 16.350000 ) weight=0.271620847042939;
			if(x>= 16.350000 &&  x < 16.500000 ) weight=0.284034755736103;
			if(x>= 16.500000 &&  x < 16.650000 ) weight=0.295330945237267;
			if(x>= 16.650000 &&  x < 16.800000 ) weight=0.307701675473455;
			if(x>= 16.800000 &&  x < 16.950000 ) weight=0.323112290632503;
			if(x>= 16.950000 &&  x < 17.100000 ) weight=0.336246628876799;
			if(x>= 17.100000 &&  x < 17.250000 ) weight=0.352532966269438;
			if(x>= 17.250000 &&  x < 17.400000 ) weight=0.366949107509111;
			if(x>= 17.400000 &&  x < 17.550000 ) weight=0.384637600201321;
			if(x>= 17.550000 &&  x < 17.700000 ) weight=0.399987967301782;
			if(x>= 17.700000 &&  x < 17.850000 ) weight=0.419063463693229;
			if(x>= 17.850000 &&  x < 18.000000 ) weight=0.436063021079290;
			if(x>= 18.000000 &&  x < 18.150000 ) weight=0.454785955372134;
			if(x>= 18.150000 &&  x < 18.300000 ) weight=0.475769803314213;
			if(x>= 18.300000 &&  x < 18.450000 ) weight=0.496294924609107;
			if(x>= 18.450000 &&  x < 18.600000 ) weight=0.516479649825334;
			if(x>= 18.600000 &&  x < 18.750000 ) weight=0.542947040817334;
			if(x>= 18.750000 &&  x < 18.900000 ) weight=0.568003916678018;
			if(x>= 18.900000 &&  x < 19.050000 ) weight=0.589790257555436;
			if(x>= 19.050000 &&  x < 19.200000 ) weight=0.613608668234929;
			if(x>= 19.200000 &&  x < 19.350000 ) weight=0.640260055788393;
			if(x>= 19.350000 &&  x < 19.500000 ) weight=0.666803515765773;
			if(x>= 19.500000 &&  x < 19.650000 ) weight=0.696145951230107;
			if(x>= 19.650000 &&  x < 19.800000 ) weight=0.731412842757920;
			if(x>= 19.800000 &&  x < 19.950000 ) weight=0.759292137399334;
			if(x>= 19.950000 &&  x < 20.100000 ) weight=0.796118415128404;
			if(x>= 20.100000 &&  x < 20.250000 ) weight=0.821925886000571;
			if(x>= 20.250000 &&  x < 20.400000 ) weight=0.862921044333283;
			if(x>= 20.400000 &&  x < 20.550000 ) weight=0.902094698164521;
			if(x>= 20.550000 &&  x < 20.700000 ) weight=0.935399714042101;
			if(x>= 20.700000 &&  x < 20.850000 ) weight=0.975764682388449;
			if(x>= 20.850000 &&  x < 21.000000 ) weight=1.016439717299811;
			if(x>= 21.000000 &&  x < 21.150000 ) weight=1.069587812333848;
			if(x>= 21.150000 &&  x < 21.300000 ) weight=1.106714657993256;
			if(x>= 21.300000 &&  x < 21.450000 ) weight=1.164379192595495;
			if(x>= 21.450000 &&  x < 21.600000 ) weight=1.204730791125141;
			if(x>= 21.600000 &&  x < 21.750000 ) weight=1.253750867953596;
			if(x>= 21.750000 &&  x < 21.900000 ) weight=1.298198277796144;
			if(x>= 21.900000 &&  x < 22.050000 ) weight=1.365911866131667;
			if(x>= 22.050000 &&  x < 22.200000 ) weight=1.420409965351406;
			if(x>= 22.200000 &&  x < 22.350000 ) weight=1.474339585686397;
			if(x>= 22.350000 &&  x < 22.500000 ) weight=1.543888585213717;
			if(x>= 22.500000 &&  x < 22.650000 ) weight=1.605218887020298;
			if(x>= 22.650000 &&  x < 22.800000 ) weight=1.663799364800050;
			if(x>= 22.800000 &&  x < 22.950000 ) weight=1.733453986618224;
			if(x>= 22.950000 &&  x < 23.100000 ) weight=1.799341089367671;
			if(x>= 23.100000 &&  x < 23.250000 ) weight=1.873410817107899;
			if(x>= 23.250000 &&  x < 23.400000 ) weight=1.950930314705597;
			if(x>= 23.400000 &&  x < 23.550000 ) weight=2.024444550904695;
			if(x>= 23.550000 &&  x < 23.700000 ) weight=2.092417910862227;
			if(x>= 23.700000 &&  x < 23.850000 ) weight=2.185178489087976;
			if(x>= 23.850000 &&  x < 24.000000 ) weight=2.268329804395500;
			if(x>= 24.000000 &&  x < 24.150000 ) weight=2.352141052045197;
			if(x>= 24.150000 &&  x < 24.300000 ) weight=2.441816162071452;
			if(x>= 24.300000 &&  x < 24.450000 ) weight=2.545920859158508;
			if(x>= 24.450000 &&  x < 24.600000 ) weight=2.636196602840454;
			if(x>= 24.600000 &&  x < 24.750000 ) weight=2.753645035582378;
			if(x>= 24.750000 &&  x < 24.900000 ) weight=2.839623744344314;
			if(x>= 24.900000 &&  x < 25.050000 ) weight=2.957570535969057;
			if(x>= 25.050000 &&  x < 25.200000 ) weight=3.031835808252132;
			if(x>= 25.200000 &&  x < 25.350000 ) weight=3.138077690664102;
			if(x>= 25.350000 &&  x < 25.500000 ) weight=3.260828230863841;
			if(x>= 25.500000 &&  x < 25.650000 ) weight=3.382107226387035;
			if(x>= 25.650000 &&  x < 25.800000 ) weight=3.496379696192462;
			if(x>= 25.800000 &&  x < 25.950000 ) weight=3.631416682660805;
			if(x>= 25.950000 &&  x < 26.100000 ) weight=3.771427718036982;
			if(x>= 26.100000 &&  x < 26.250000 ) weight=3.894094319654906;
			if(x>= 26.250000 &&  x < 26.400000 ) weight=4.031227097596314;
			if(x>= 26.400000 &&  x < 26.550000 ) weight=4.202940627415375;
			if(x>= 26.550000 &&  x < 26.700000 ) weight=4.304205533392100;
			if(x>= 26.700000 &&  x < 26.850000 ) weight=4.480025238536920;
			if(x>= 26.850000 &&  x < 27.000000 ) weight=4.617290614178405;
			if(x>= 27.000000 &&  x < 27.150000 ) weight=4.754741102548583;
			if(x>= 27.150000 &&  x < 27.300000 ) weight=4.942724022050025;
			if(x>= 27.300000 &&  x < 27.450000 ) weight=5.110932453738532;
			if(x>= 27.450000 &&  x < 27.600000 ) weight=5.252936251746900;
			if(x>= 27.600000 &&  x < 27.750000 ) weight=5.455536324784713;
			if(x>= 27.750000 &&  x < 27.900000 ) weight=5.683096462324079;
			if(x>= 27.900000 &&  x < 28.050000 ) weight=5.814438649389198;
			if(x>= 28.050000 &&  x < 28.200000 ) weight=6.053990717510620;
			if(x>= 28.200000 &&  x < 28.350000 ) weight=6.211760944001412;
			if(x>= 28.350000 &&  x < 28.500000 ) weight=6.387277837640890;
			if(x>= 28.500000 &&  x < 28.650000 ) weight=6.614662917912030;
			if(x>= 28.650000 &&  x < 28.800000 ) weight=6.801227658255996;
			if(x>= 28.800000 &&  x < 28.950000 ) weight=7.012864613446635;
			if(x>= 28.950000 &&  x < 29.100000 ) weight=7.335284103163092;
			if(x>= 29.100000 &&  x < 29.250000 ) weight=7.483078537777679;
			if(x>= 29.250000 &&  x < 29.400000 ) weight=7.730139010094280;
			if(x>= 29.400000 &&  x < 29.550000 ) weight=8.034756753516147;
			if(x>= 29.550000 &&  x < 29.700000 ) weight=8.237023811422020;
			if(x>= 29.700000 &&  x < 29.850000 ) weight=8.452577565136050;
			if(x>= 29.850000 &&  x < 30.000000 ) weight=8.709943321215590;
			if(x> 30.000000 ) weight=8.709943321215590;
		}
		
			//eta
        else if(pdg_particle==221)
        {
			if(x>= 0.000000 &&  x < 0.150000 ) weight=0.000065175335984;
			if(x>= 0.150000 &&  x < 0.300000 ) weight=0.000044681150023;
			if(x>= 0.300000 &&  x < 0.450000 ) weight=0.000037799534416;
			if(x>= 0.450000 &&  x < 0.600000 ) weight=0.000032157721810;
			if(x>= 0.600000 &&  x < 0.750000 ) weight=0.000027687031685;
			if(x>= 0.750000 &&  x < 0.900000 ) weight=0.000024575471051;
			if(x>= 0.900000 &&  x < 1.050000 ) weight=0.000022604037864;
			if(x>= 1.050000 &&  x < 1.200000 ) weight=0.000021745812264;
			if(x>= 1.200000 &&  x < 1.350000 ) weight=0.000021814592936;
			if(x>= 1.350000 &&  x < 1.500000 ) weight=0.000022789125995;
			if(x>= 1.500000 &&  x < 1.650000 ) weight=0.000024380785149;
			if(x>= 1.650000 &&  x < 1.800000 ) weight=0.000026753293258;
			if(x>= 1.800000 &&  x < 1.950000 ) weight=0.000029969411837;
			if(x>= 1.950000 &&  x < 2.100000 ) weight=0.000033907213500;
			if(x>= 2.100000 &&  x < 2.250000 ) weight=0.000038854650877;
			if(x>= 2.250000 &&  x < 2.400000 ) weight=0.000045027645335;
			if(x>= 2.400000 &&  x < 2.550000 ) weight=0.000052891618188;
			if(x>= 2.550000 &&  x < 2.700000 ) weight=0.000062256188357;
			if(x>= 2.700000 &&  x < 2.850000 ) weight=0.000072987810231;
			if(x>= 2.850000 &&  x < 3.000000 ) weight=0.000086217965584;
			if(x>= 3.000000 &&  x < 3.150000 ) weight=0.000102145546367;
			if(x>= 3.150000 &&  x < 3.300000 ) weight=0.000120930982866;
			if(x>= 3.300000 &&  x < 3.450000 ) weight=0.000142989407482;
			if(x>= 3.450000 &&  x < 3.600000 ) weight=0.000168287135463;
			if(x>= 3.600000 &&  x < 3.750000 ) weight=0.000199367550625;
			if(x>= 3.750000 &&  x < 3.900000 ) weight=0.000233922135062;
			if(x>= 3.900000 &&  x < 4.050000 ) weight=0.000274448872307;
			if(x>= 4.050000 &&  x < 4.200000 ) weight=0.000323211132637;
			if(x>= 4.200000 &&  x < 4.350000 ) weight=0.000377569759880;
			if(x>= 4.350000 &&  x < 4.500000 ) weight=0.000437299495993;
			if(x>= 4.500000 &&  x < 4.650000 ) weight=0.000513783102382;
			if(x>= 4.650000 &&  x < 4.800000 ) weight=0.000595811122102;
			if(x>= 4.800000 &&  x < 4.950000 ) weight=0.000691015523760;
			if(x>= 4.950000 &&  x < 5.100000 ) weight=0.000789761337639;
			if(x>= 5.100000 &&  x < 5.250000 ) weight=0.000915962204830;
			if(x>= 5.250000 &&  x < 5.400000 ) weight=0.001042990379102;
			if(x>= 5.400000 &&  x < 5.550000 ) weight=0.001199612242151;
			if(x>= 5.550000 &&  x < 5.700000 ) weight=0.001357881175798;
			if(x>= 5.700000 &&  x < 5.850000 ) weight=0.001569969307762;
			if(x>= 5.850000 &&  x < 6.000000 ) weight=0.001784828962964;
			if(x>= 6.000000 &&  x < 6.150000 ) weight=0.002015553119291;
			if(x>= 6.150000 &&  x < 6.300000 ) weight=0.002292551345214;
			if(x>= 6.300000 &&  x < 6.450000 ) weight=0.002575749375458;
			if(x>= 6.450000 &&  x < 6.600000 ) weight=0.002906298396496;
			if(x>= 6.600000 &&  x < 6.750000 ) weight=0.003262511973776;
			if(x>= 6.750000 &&  x < 6.900000 ) weight=0.003644173577950;
			if(x>= 6.900000 &&  x < 7.050000 ) weight=0.004126740552315;
			if(x>= 7.050000 &&  x < 7.200000 ) weight=0.004585986280931;
			if(x>= 7.200000 &&  x < 7.350000 ) weight=0.005091067257979;
			if(x>= 7.350000 &&  x < 7.500000 ) weight=0.005679263155812;
			if(x>= 7.500000 &&  x < 7.650000 ) weight=0.006384482903785;
			if(x>= 7.650000 &&  x < 7.800000 ) weight=0.007000974856089;
			if(x>= 7.800000 &&  x < 7.950000 ) weight=0.007791905123154;
			if(x>= 7.950000 &&  x < 8.100000 ) weight=0.008646649696691;
			if(x>= 8.100000 &&  x < 8.250000 ) weight=0.009565183908597;
			if(x>= 8.250000 &&  x < 8.400000 ) weight=0.010545013916087;
			if(x>= 8.400000 &&  x < 8.550000 ) weight=0.011553089462211;
			if(x>= 8.550000 &&  x < 8.700000 ) weight=0.012785432247946;
			if(x>= 8.700000 &&  x < 8.850000 ) weight=0.014013490692818;
			if(x>= 8.850000 &&  x < 9.000000 ) weight=0.015371828560037;
			if(x>= 9.000000 &&  x < 9.150000 ) weight=0.016893931952376;
			if(x>= 9.150000 &&  x < 9.300000 ) weight=0.018395857435517;
			if(x>= 9.300000 &&  x < 9.450000 ) weight=0.020232279148977;
			if(x>= 9.450000 &&  x < 9.600000 ) weight=0.022277848201788;
			if(x>= 9.600000 &&  x < 9.750000 ) weight=0.024051257713961;
			if(x>= 9.750000 &&  x < 9.900000 ) weight=0.026018161699593;
			if(x>= 9.900000 &&  x < 10.050000 ) weight=0.028310434857523;
			if(x>= 10.050000 &&  x < 10.200000 ) weight=0.030931053372065;
			if(x>= 10.200000 &&  x < 10.350000 ) weight=0.033562523714702;
			if(x>= 10.350000 &&  x < 10.500000 ) weight=0.036492540519561;
			if(x>= 10.500000 &&  x < 10.650000 ) weight=0.039523746691989;
			if(x>= 10.650000 &&  x < 10.800000 ) weight=0.042492765252348;
			if(x>= 10.800000 &&  x < 10.950000 ) weight=0.046027593361174;
			if(x>= 10.950000 &&  x < 11.100000 ) weight=0.049710321411214;
			if(x>= 11.100000 &&  x < 11.250000 ) weight=0.053878257994145;
			if(x>= 11.250000 &&  x < 11.400000 ) weight=0.057798493082901;
			if(x>= 11.400000 &&  x < 11.550000 ) weight=0.062734138089187;
			if(x>= 11.550000 &&  x < 11.700000 ) weight=0.066669672324991;
			if(x>= 11.700000 &&  x < 11.850000 ) weight=0.072073301135872;
			if(x>= 11.850000 &&  x < 12.000000 ) weight=0.077460783642178;
			if(x>= 12.000000 &&  x < 12.150000 ) weight=0.082961528852578;
			if(x>= 12.150000 &&  x < 12.300000 ) weight=0.089045833887930;
			if(x>= 12.300000 &&  x < 12.450000 ) weight=0.095617784553570;
			if(x>= 12.450000 &&  x < 12.600000 ) weight=0.103100010376093;
			if(x>= 12.600000 &&  x < 12.750000 ) weight=0.109139396172481;
			if(x>= 12.750000 &&  x < 12.900000 ) weight=0.117994960409230;
			if(x>= 12.900000 &&  x < 13.050000 ) weight=0.125504848074421;
			if(x>= 13.050000 &&  x < 13.200000 ) weight=0.134304942271888;
			if(x>= 13.200000 &&  x < 13.350000 ) weight=0.142488624950247;
			if(x>= 13.350000 &&  x < 13.500000 ) weight=0.153359812210462;
			if(x>= 13.500000 &&  x < 13.650000 ) weight=0.162192653523171;
			if(x>= 13.650000 &&  x < 13.800000 ) weight=0.172575939355670;
			if(x>= 13.800000 &&  x < 13.950000 ) weight=0.185679476781945;
			if(x>= 13.950000 &&  x < 14.100000 ) weight=0.195482957340583;
			if(x>= 14.100000 &&  x < 14.250000 ) weight=0.210951786393359;
			if(x>= 14.250000 &&  x < 14.400000 ) weight=0.219734975397937;
			if(x>= 14.400000 &&  x < 14.550000 ) weight=0.236549799674666;
			if(x>= 14.550000 &&  x < 14.700000 ) weight=0.251295714029736;
			if(x>= 14.700000 &&  x < 14.850000 ) weight=0.264622881165939;
			if(x>= 14.850000 &&  x < 15.000000 ) weight=0.281678989817342;
			if(x>= 15.000000 &&  x < 15.150000 ) weight=0.299065772410259;
			if(x>= 15.150000 &&  x < 15.300000 ) weight=0.315542413071123;
			if(x>= 15.300000 &&  x < 15.450000 ) weight=0.335680690280426;
			if(x>= 15.450000 &&  x < 15.600000 ) weight=0.353290836330913;
			if(x>= 15.600000 &&  x < 15.750000 ) weight=0.377605144964225;
			if(x>= 15.750000 &&  x < 15.900000 ) weight=0.398199816798876;
			if(x>= 15.900000 &&  x < 16.050000 ) weight=0.420888945633927;
			if(x>= 16.050000 &&  x < 16.200000 ) weight=0.438942217441571;
			if(x>= 16.200000 &&  x < 16.350000 ) weight=0.465102065950107;
			if(x>= 16.350000 &&  x < 16.500000 ) weight=0.488784172376542;
			if(x>= 16.500000 &&  x < 16.650000 ) weight=0.520629823967140;
			if(x>= 16.650000 &&  x < 16.800000 ) weight=0.546213640279334;
			if(x>= 16.800000 &&  x < 16.950000 ) weight=0.573661576794979;
			if(x>= 16.950000 &&  x < 17.100000 ) weight=0.604576647838019;
			if(x>= 17.100000 &&  x < 17.250000 ) weight=0.642259408943007;
			if(x>= 17.250000 &&  x < 17.400000 ) weight=0.674479678799181;
			if(x>= 17.400000 &&  x < 17.550000 ) weight=0.709419606682057;
			if(x>= 17.550000 &&  x < 17.700000 ) weight=0.744364329859029;
			if(x>= 17.700000 &&  x < 17.850000 ) weight=0.782917466814334;
			if(x>= 17.850000 &&  x < 18.000000 ) weight=0.815720215521759;
			if(x>= 18.000000 &&  x < 18.150000 ) weight=0.863320366921488;
			if(x>= 18.150000 &&  x < 18.300000 ) weight=0.915989900254712;
			if(x>= 18.300000 &&  x < 18.450000 ) weight=0.951359989788585;
			if(x>= 18.450000 &&  x < 18.600000 ) weight=0.997464758773388;
			if(x>= 18.600000 &&  x < 18.750000 ) weight=1.046252201317995;
			if(x>= 18.750000 &&  x < 18.900000 ) weight=1.096565813374950;
			if(x>= 18.900000 &&  x < 19.050000 ) weight=1.142566741793662;
			if(x>= 19.050000 &&  x < 19.200000 ) weight=1.207236212372948;
			if(x>= 19.200000 &&  x < 19.350000 ) weight=1.267706692203010;
			if(x>= 19.350000 &&  x < 19.500000 ) weight=1.314869847290634;
			if(x>= 19.500000 &&  x < 19.650000 ) weight=1.371705655434577;
			if(x>= 19.650000 &&  x < 19.800000 ) weight=1.450489158313063;
			if(x>= 19.800000 &&  x < 19.950000 ) weight=1.525013434090324;
			if(x>= 19.950000 &&  x < 20.100000 ) weight=1.578478248463832;
			if(x>= 20.100000 &&  x < 20.250000 ) weight=1.653077487719245;
			if(x>= 20.250000 &&  x < 20.400000 ) weight=1.744712312559263;
			if(x>= 20.400000 &&  x < 20.550000 ) weight=1.803781177155982;
			if(x>= 20.550000 &&  x < 20.700000 ) weight=1.869585533128399;
			if(x>= 20.700000 &&  x < 20.850000 ) weight=1.962461126503244;
			if(x>= 20.850000 &&  x < 21.000000 ) weight=2.074866403748616;
			if(x>= 21.000000 &&  x < 21.150000 ) weight=2.149665894768486;
			if(x>= 21.150000 &&  x < 21.300000 ) weight=2.255490081273915;
			if(x>= 21.300000 &&  x < 21.450000 ) weight=2.329785347091331;
			if(x>= 21.450000 &&  x < 21.600000 ) weight=2.442630074898533;
			if(x>= 21.600000 &&  x < 21.750000 ) weight=2.534998115373910;
			if(x>= 21.750000 &&  x < 21.900000 ) weight=2.637893121752631;
			if(x>= 21.900000 &&  x < 22.050000 ) weight=2.758809272311342;
			if(x>= 22.050000 &&  x < 22.200000 ) weight=2.870514129035438;
			if(x>= 22.200000 &&  x < 22.350000 ) weight=2.959877535254340;
			if(x>= 22.350000 &&  x < 22.500000 ) weight=3.099623955106531;
			if(x>= 22.500000 &&  x < 22.650000 ) weight=3.218600049881187;
			if(x>= 22.650000 &&  x < 22.800000 ) weight=3.363318999360930;
			if(x>= 22.800000 &&  x < 22.950000 ) weight=3.510375140919558;
			if(x>= 22.950000 &&  x < 23.100000 ) weight=3.633308110640943;
			if(x>= 23.100000 &&  x < 23.250000 ) weight=3.811620204597428;
			if(x>= 23.250000 &&  x < 23.400000 ) weight=3.909946884230240;
			if(x>= 23.400000 &&  x < 23.550000 ) weight=4.084581966778044;
			if(x>= 23.550000 &&  x < 23.700000 ) weight=4.235212093003118;
			if(x>= 23.700000 &&  x < 23.850000 ) weight=4.408964518119510;
			if(x>= 23.850000 &&  x < 24.000000 ) weight=4.593886463624743;
			if(x>= 24.000000 &&  x < 24.150000 ) weight=4.768985861813210;
			if(x>= 24.150000 &&  x < 24.300000 ) weight=4.953474003190913;
			if(x>= 24.300000 &&  x < 24.450000 ) weight=5.141781321375658;
			if(x>= 24.450000 &&  x < 24.600000 ) weight=5.298862454339663;
			if(x>= 24.600000 &&  x < 24.750000 ) weight=5.549706511082033;
			if(x>= 24.750000 &&  x < 24.900000 ) weight=5.731727881477775;
			if(x>= 24.900000 &&  x < 25.050000 ) weight=5.906810309834845;
			if(x>= 25.050000 &&  x < 25.200000 ) weight=6.188400891210357;
			if(x>= 25.200000 &&  x < 25.350000 ) weight=6.418641544791604;
			if(x>= 25.350000 &&  x < 25.500000 ) weight=6.648985175889045;
			if(x>= 25.500000 &&  x < 25.650000 ) weight=6.819805148386722;
			if(x>= 25.650000 &&  x < 25.800000 ) weight=7.123971985274818;
			if(x>= 25.800000 &&  x < 25.950000 ) weight=7.328398338167448;
			if(x>= 25.950000 &&  x < 26.100000 ) weight=7.669741746768523;
			if(x>= 26.100000 &&  x < 26.250000 ) weight=7.866416058508633;
			if(x>= 26.250000 &&  x < 26.400000 ) weight=8.137288389451907;
			if(x>= 26.400000 &&  x < 26.550000 ) weight=8.388942733356441;
			if(x>= 26.550000 &&  x < 26.700000 ) weight=8.715474188724365;
			if(x>= 26.700000 &&  x < 26.850000 ) weight=9.070673184304725;
			if(x>= 26.850000 &&  x < 27.000000 ) weight=9.305742481860554;
			if(x>= 27.000000 &&  x < 27.150000 ) weight=9.567027852817075;
			if(x>= 27.150000 &&  x < 27.300000 ) weight=9.969499075607356;
			if(x>= 27.300000 &&  x < 27.450000 ) weight=10.292369739723279;
			if(x>= 27.450000 &&  x < 27.600000 ) weight=10.671093848594357;
			if(x>= 27.600000 &&  x < 27.750000 ) weight=11.029704327560799;
			if(x>= 27.750000 &&  x < 27.900000 ) weight=11.445768220762435;
			if(x>= 27.900000 &&  x < 28.050000 ) weight=11.907708790189581;
			if(x>= 28.050000 &&  x < 28.200000 ) weight=12.205228020287439;
			if(x>= 28.200000 &&  x < 28.350000 ) weight=12.591522367964483;
			if(x>= 28.350000 &&  x < 28.500000 ) weight=12.891966786517502;
			if(x>= 28.500000 &&  x < 28.650000 ) weight=13.301975455858535;
			if(x>= 28.650000 &&  x < 28.800000 ) weight=13.801942868418852;
			if(x>= 28.800000 &&  x < 28.950000 ) weight=14.147710803982289;
			if(x>= 28.950000 &&  x < 29.100000 ) weight=14.743631512821393;
			if(x>= 29.100000 &&  x < 29.250000 ) weight=15.159724557487563;
			if(x>= 29.250000 &&  x < 29.400000 ) weight=15.578306127919422;
			if(x>= 29.400000 &&  x < 29.550000 ) weight=16.203696478367558;
			if(x>= 29.550000 &&  x < 29.700000 ) weight=16.647515455014407;
			if(x>= 29.700000 &&  x < 29.850000 ) weight=17.105994008780360;
			if(x>= 29.850000 &&  x < 30.000000 ) weight=17.550703202326979;
			if(x>30.000000 ) weight=17.550703202326979;
			
			
		}
		
		/*
		 if(pdg_particle==111){
		 if(x>= 0.000000 &&  x < 0.150000 ) weight=5.413029;
		 if(x>= 0.150000 &&  x < 0.300000 ) weight=1.843465;
		 if(x>= 0.300000 &&  x < 0.450000 ) weight=1.562291;
		 if(x>= 0.450000 &&  x < 0.600000 ) weight=1.390412;
		 if(x>= 0.600000 &&  x < 0.750000 ) weight=1.246792;
		 if(x>= 0.750000 &&  x < 0.900000 ) weight=1.146358;
		 if(x>= 0.900000 &&  x < 1.050000 ) weight=1.099458;
		 if(x>= 1.050000 &&  x < 1.200000 ) weight=1.108490;
		 if(x>= 1.200000 &&  x < 1.350000 ) weight=1.168142;
		 if(x>= 1.350000 &&  x < 1.500000 ) weight=1.279816;
		 if(x>= 1.500000 &&  x < 1.650000 ) weight=1.447372;
		 if(x>= 1.650000 &&  x < 1.800000 ) weight=1.675200;
		 if(x>= 1.800000 &&  x < 1.950000 ) weight=1.965825;
		 if(x>= 1.950000 &&  x < 2.100000 ) weight=2.334030;
		 if(x>= 2.100000 &&  x < 2.250000 ) weight=2.812404;
		 if(x>= 2.250000 &&  x < 2.400000 ) weight=3.399171;
		 if(x>= 2.400000 &&  x < 2.550000 ) weight=4.113419;
		 if(x>= 2.550000 &&  x < 2.700000 ) weight=5.028714;
		 if(x>= 2.700000 &&  x < 2.850000 ) weight=6.067468;
		 if(x>= 2.850000 &&  x < 3.000000 ) weight=7.341577;
		 if(x>= 3.000000 &&  x < 3.150000 ) weight=8.812074;
		 if(x>= 3.150000 &&  x < 3.300000 ) weight=10.618540;
		 if(x>= 3.300000 &&  x < 3.450000 ) weight=12.685035;
		 if(x>= 3.450000 &&  x < 3.600000 ) weight=15.005012;
		 if(x>= 3.600000 &&  x < 3.750000 ) weight=17.725825;
		 if(x>= 3.750000 &&  x < 3.900000 ) weight=20.871579;
		 if(x>= 3.900000 &&  x < 4.050000 ) weight=24.360849;
		 if(x>= 4.050000 &&  x < 4.200000 ) weight=28.488189;
		 if(x>= 4.200000 &&  x < 4.350000 ) weight=32.938705;
		 if(x>= 4.350000 &&  x < 4.500000 ) weight=38.004920;
		 if(x>= 4.500000 &&  x < 4.650000 ) weight=43.905847;
		 if(x>= 4.650000 &&  x < 4.800000 ) weight=50.316566;
		 if(x>= 4.800000 &&  x < 4.950000 ) weight=57.346497;
		 if(x>= 4.950000 &&  x < 5.100000 ) weight=65.527849;
		 if(x>= 5.100000 &&  x < 5.250000 ) weight=74.047274;
		 if(x>= 5.250000 &&  x < 5.400000 ) weight=83.615531;
		 if(x>= 5.400000 &&  x < 5.550000 ) weight=94.362690;
		 if(x>= 5.550000 &&  x < 5.700000 ) weight=106.045286;
		 if(x>= 5.700000 &&  x < 5.850000 ) weight=119.012986;
		 if(x>= 5.850000 &&  x < 6.000000 ) weight=132.880543;
		 if(x>= 6.000000 &&  x < 6.150000 ) weight=147.654408;
		 if(x>= 6.150000 &&  x < 6.300000 ) weight=164.856176;
		 if(x>= 6.300000 &&  x < 6.450000 ) weight=182.443341;
		 if(x>= 6.450000 &&  x < 6.600000 ) weight=201.580388;
		 if(x>= 6.600000 &&  x < 6.750000 ) weight=222.096149;
		 if(x>= 6.750000 &&  x < 6.900000 ) weight=245.731663;
		 if(x>= 6.900000 &&  x < 7.050000 ) weight=269.694559;
		 if(x>= 7.050000 &&  x < 7.200000 ) weight=296.517907;
		 if(x>= 7.200000 &&  x < 7.350000 ) weight=323.392152;
		 if(x>= 7.350000 &&  x < 7.500000 ) weight=352.243121;
		 if(x>= 7.500000 &&  x < 7.650000 ) weight=385.261159;
		 if(x>= 7.650000 &&  x < 7.800000 ) weight=418.118716;
		 if(x>= 7.800000 &&  x < 7.950000 ) weight=455.113426;
		 if(x>= 7.950000 &&  x < 8.100000 ) weight=495.331199;
		 if(x>= 8.100000 &&  x < 8.250000 ) weight=535.194415;
		 if(x>= 8.250000 &&  x < 8.400000 ) weight=580.109285;
		 if(x>= 8.400000 &&  x < 8.550000 ) weight=626.069818;
		 if(x>= 8.550000 &&  x < 8.700000 ) weight=676.485515;
		 if(x>= 8.700000 &&  x < 8.850000 ) weight=727.181488;
		 if(x>= 8.850000 &&  x < 9.000000 ) weight=780.772898;
		 if(x>= 9.000000 &&  x < 9.150000 ) weight=843.460957;
		 if(x>= 9.150000 &&  x < 9.300000 ) weight=903.110782;
		 if(x>= 9.300000 &&  x < 9.450000 ) weight=966.468529;
		 if(x>= 9.450000 &&  x < 9.600000 ) weight=1036.956505;
		 if(x>= 9.600000 &&  x < 9.750000 ) weight=1100.209143;
		 if(x>= 9.750000 &&  x < 9.900000 ) weight=1180.666364;
		 if(x>= 9.900000 &&  x < 10.050000 ) weight=1262.838652;
		 if(x>= 10.050000 &&  x < 10.200000 ) weight=1345.690689;
		 if(x>= 10.200000 &&  x < 10.350000 ) weight=1441.650246;
		 if(x>= 10.350000 &&  x < 10.500000 ) weight=1522.357306;
		 if(x>= 10.500000 &&  x < 10.650000 ) weight=1622.668606;
		 if(x>= 10.650000 &&  x < 10.800000 ) weight=1717.472192;
		 if(x>= 10.800000 &&  x < 10.950000 ) weight=1818.572721;
		 if(x>= 10.950000 &&  x < 11.100000 ) weight=1943.911193;
		 if(x>= 11.100000 &&  x < 11.250000 ) weight=2054.165601;
		 if(x>= 11.250000 &&  x < 11.400000 ) weight=2176.172676;
		 if(x>= 11.400000 &&  x < 11.550000 ) weight=2288.558855;
		 if(x>= 11.550000 &&  x < 11.700000 ) weight=2413.537715;
		 if(x>= 11.700000 &&  x < 11.850000 ) weight=2557.627323;
		 if(x>= 11.850000 &&  x < 12.000000 ) weight=2706.099550;
		 if(x>= 12.000000 &&  x < 12.150000 ) weight=2851.892790;
		 if(x>= 12.150000 &&  x < 12.300000 ) weight=2986.992952;
		 if(x>= 12.300000 &&  x < 12.450000 ) weight=3150.172254;
		 if(x>= 12.450000 &&  x < 12.600000 ) weight=3356.781323;
		 if(x>= 12.600000 &&  x < 12.750000 ) weight=3509.653394;
		 if(x>= 12.750000 &&  x < 12.900000 ) weight=3703.187316;
		 if(x>= 12.900000 &&  x < 13.050000 ) weight=3874.919703;
		 if(x>= 13.050000 &&  x < 13.200000 ) weight=4088.973987;
		 if(x>= 13.200000 &&  x < 13.350000 ) weight=4319.189310;
		 if(x>= 13.350000 &&  x < 13.500000 ) weight=4474.371121;
		 if(x>= 13.500000 &&  x < 13.650000 ) weight=4710.164144;
		 if(x>= 13.650000 &&  x < 13.800000 ) weight=4962.722731;
		 if(x>= 13.800000 &&  x < 13.950000 ) weight=5212.695189;
		 if(x>= 13.950000 &&  x < 14.100000 ) weight=5455.828904;
		 if(x>= 14.100000 &&  x < 14.250000 ) weight=5747.811891;
		 if(x>= 14.250000 &&  x < 14.400000 ) weight=6027.093776;
		 if(x>= 14.400000 &&  x < 14.550000 ) weight=6247.732360;
		 if(x>= 14.550000 &&  x < 14.700000 ) weight=6571.407631;
		 if(x>= 14.700000 &&  x < 14.850000 ) weight=6886.337364;
		 if(x>= 14.850000 &&  x < 15.000000 ) weight=7217.840294;
		 if(x>= 15.000000 &&  x < 15.150000 ) weight=7618.389897;
		 if(x>= 15.150000 &&  x < 15.300000 ) weight=7868.994076;
		 if(x>= 15.300000 &&  x < 15.450000 ) weight=8257.157222;
		 if(x>= 15.450000 &&  x < 15.600000 ) weight=8670.134419;
		 if(x>= 15.600000 &&  x < 15.750000 ) weight=9084.681221;
		 if(x>= 15.750000 &&  x < 15.900000 ) weight=9473.392377;
		 if(x>= 15.900000 &&  x < 16.050000 ) weight=9927.918277;
		 if(x>= 16.050000 &&  x < 16.200000 ) weight=10377.021056;
		 if(x>= 16.200000 &&  x < 16.350000 ) weight=10864.833882;
		 if(x>= 16.350000 &&  x < 16.500000 ) weight=11361.390229;
		 if(x>= 16.500000 &&  x < 16.650000 ) weight=11813.237809;
		 if(x>= 16.650000 &&  x < 16.800000 ) weight=12308.067019;
		 if(x>= 16.800000 &&  x < 16.950000 ) weight=12924.491625;
		 if(x>= 16.950000 &&  x < 17.100000 ) weight=13449.865155;
		 if(x>= 17.100000 &&  x < 17.250000 ) weight=14101.318651;
		 if(x>= 17.250000 &&  x < 17.400000 ) weight=14677.964300;
		 if(x>= 17.400000 &&  x < 17.550000 ) weight=15385.504008;
		 if(x>= 17.550000 &&  x < 17.700000 ) weight=15999.518692;
		 if(x>= 17.700000 &&  x < 17.850000 ) weight=16762.538548;
		 if(x>= 17.850000 &&  x < 18.000000 ) weight=17442.520843;
		 if(x>= 18.000000 &&  x < 18.150000 ) weight=18191.438215;
		 if(x>= 18.150000 &&  x < 18.300000 ) weight=19030.792133;
		 if(x>= 18.300000 &&  x < 18.450000 ) weight=19851.796984;
		 if(x>= 18.450000 &&  x < 18.600000 ) weight=20659.185993;
		 if(x>= 18.600000 &&  x < 18.750000 ) weight=21717.881633;
		 if(x>= 18.750000 &&  x < 18.900000 ) weight=22720.156667;
		 if(x>= 18.900000 &&  x < 19.050000 ) weight=23591.610302;
		 if(x>= 19.050000 &&  x < 19.200000 ) weight=24544.346729;
		 if(x>= 19.200000 &&  x < 19.350000 ) weight=25610.402232;
		 if(x>= 19.350000 &&  x < 19.500000 ) weight=26672.140631;
		 if(x>= 19.500000 &&  x < 19.650000 ) weight=27845.838049;
		 if(x>= 19.650000 &&  x < 19.800000 ) weight=29256.513710;
		 if(x>= 19.800000 &&  x < 19.950000 ) weight=30371.685496;
		 if(x>= 19.950000 &&  x < 20.100000 ) weight=31844.736605;
		 if(x>= 20.100000 &&  x < 20.250000 ) weight=32877.035440;
		 if(x>= 20.250000 &&  x < 20.400000 ) weight=34516.841773;
		 if(x>= 20.400000 &&  x < 20.550000 ) weight=36083.787927;
		 if(x>= 20.550000 &&  x < 20.700000 ) weight=37415.988562;
		 if(x>= 20.700000 &&  x < 20.850000 ) weight=39030.587296;
		 if(x>= 20.850000 &&  x < 21.000000 ) weight=40657.588692;
		 if(x>= 21.000000 &&  x < 21.150000 ) weight=42783.512493;
		 if(x>= 21.150000 &&  x < 21.300000 ) weight=44268.586320;
		 if(x>= 21.300000 &&  x < 21.450000 ) weight=46575.167704;
		 if(x>= 21.450000 &&  x < 21.600000 ) weight=48189.231645;
		 if(x>= 21.600000 &&  x < 21.750000 ) weight=50150.034718;
		 if(x>= 21.750000 &&  x < 21.900000 ) weight=51927.931112;
		 if(x>= 21.900000 &&  x < 22.050000 ) weight=54636.474645;
		 if(x>= 22.050000 &&  x < 22.200000 ) weight=56816.398614;
		 if(x>= 22.200000 &&  x < 22.350000 ) weight=58973.583427;
		 if(x>= 22.350000 &&  x < 22.500000 ) weight=61755.543409;
		 if(x>= 22.500000 &&  x < 22.650000 ) weight=64208.755481;
		 if(x>= 22.650000 &&  x < 22.800000 ) weight=66551.974592;
		 if(x>= 22.800000 &&  x < 22.950000 ) weight=69338.159465;
		 if(x>= 22.950000 &&  x < 23.100000 ) weight=71973.643575;
		 if(x>= 23.100000 &&  x < 23.250000 ) weight=74936.432684;
		 if(x>= 23.250000 &&  x < 23.400000 ) weight=78037.212588;
		 if(x>= 23.400000 &&  x < 23.550000 ) weight=80977.782036;
		 if(x>= 23.550000 &&  x < 23.700000 ) weight=83696.716434;
		 if(x>= 23.700000 &&  x < 23.850000 ) weight=87407.139564;
		 if(x>= 23.850000 &&  x < 24.000000 ) weight=90733.192176;
		 if(x>= 24.000000 &&  x < 24.150000 ) weight=94085.642082;
		 if(x>= 24.150000 &&  x < 24.300000 ) weight=97672.646483;
		 if(x>= 24.300000 &&  x < 24.450000 ) weight=101836.834366;
		 if(x>= 24.450000 &&  x < 24.600000 ) weight=105447.864114;
		 if(x>= 24.600000 &&  x < 24.750000 ) weight=110145.801423;
		 if(x>= 24.750000 &&  x < 24.900000 ) weight=113584.949774;
		 if(x>= 24.900000 &&  x < 25.050000 ) weight=118302.821439;
		 if(x>= 25.050000 &&  x < 25.200000 ) weight=121273.432330;
		 if(x>= 25.200000 &&  x < 25.350000 ) weight=125523.107627;
		 if(x>= 25.350000 &&  x < 25.500000 ) weight=130433.129235;
		 if(x>= 25.500000 &&  x < 25.650000 ) weight=135284.289055;
		 if(x>= 25.650000 &&  x < 25.800000 ) weight=139855.187848;
		 if(x>= 25.800000 &&  x < 25.950000 ) weight=145256.667306;
		 if(x>= 25.950000 &&  x < 26.100000 ) weight=150857.108721;
		 if(x>= 26.100000 &&  x < 26.250000 ) weight=155763.772786;
		 if(x>= 26.250000 &&  x < 26.400000 ) weight=161249.083904;
		 if(x>= 26.400000 &&  x < 26.550000 ) weight=168117.625097;
		 if(x>= 26.550000 &&  x < 26.700000 ) weight=172168.221336;
		 if(x>= 26.700000 &&  x < 26.850000 ) weight=179201.009541;
		 if(x>= 26.850000 &&  x < 27.000000 ) weight=184691.624567;
		 if(x>= 27.000000 &&  x < 27.150000 ) weight=190189.644102;
		 if(x>= 27.150000 &&  x < 27.300000 ) weight=197708.960882;
		 if(x>= 27.300000 &&  x < 27.450000 ) weight=204437.298150;
		 if(x>= 27.450000 &&  x < 27.600000 ) weight=210117.450070;
		 if(x>= 27.600000 &&  x < 27.750000 ) weight=218221.452991;
		 if(x>= 27.750000 &&  x < 27.900000 ) weight=227323.858493;
		 if(x>= 27.900000 &&  x < 28.050000 ) weight=232577.545976;
		 if(x>= 28.050000 &&  x < 28.200000 ) weight=242159.628700;
		 if(x>= 28.200000 &&  x < 28.350000 ) weight=248470.437760;
		 if(x>= 28.350000 &&  x < 28.500000 ) weight=255491.113506;
		 if(x>= 28.500000 &&  x < 28.650000 ) weight=264586.516716;
		 if(x>= 28.650000 &&  x < 28.800000 ) weight=272049.106330;
		 if(x>= 28.800000 &&  x < 28.950000 ) weight=280514.584538;
		 if(x>= 28.950000 &&  x < 29.100000 ) weight=293411.364127;
		 if(x>= 29.100000 &&  x < 29.250000 ) weight=299323.141511;
		 if(x>= 29.250000 &&  x < 29.400000 ) weight=309205.560404;
		 if(x>= 29.400000 &&  x < 29.550000 ) weight=321390.270141;
		 if(x>= 29.550000 &&  x < 29.700000 ) weight=329480.952457;
		 if(x>= 29.700000 &&  x < 29.850000 ) weight=338103.102605;
		 if(x>= 29.850000 &&  x < 30.000000 ) weight=348397.732849;
		 if(x>= 30) weight=348397.732849;
		 
		 }
		 //eta
		 else if(pdg_particle==221)
		 {
		 if(x>= 0.000000 &&  x < 0.150000 ) weight=2.607013;
		 if(x>= 0.150000 &&  x < 0.300000 ) weight=1.787246;
		 if(x>= 0.300000 &&  x < 0.450000 ) weight=1.511981;
		 if(x>= 0.450000 &&  x < 0.600000 ) weight=1.286309;
		 if(x>= 0.600000 &&  x < 0.750000 ) weight=1.107481;
		 if(x>= 0.750000 &&  x < 0.900000 ) weight=0.983019;
		 if(x>= 0.900000 &&  x < 1.050000 ) weight=0.904162;
		 if(x>= 1.050000 &&  x < 1.200000 ) weight=0.869832;
		 if(x>= 1.200000 &&  x < 1.350000 ) weight=0.872584;
		 if(x>= 1.350000 &&  x < 1.500000 ) weight=0.911565;
		 if(x>= 1.500000 &&  x < 1.650000 ) weight=0.975231;
		 if(x>= 1.650000 &&  x < 1.800000 ) weight=1.070132;
		 if(x>= 1.800000 &&  x < 1.950000 ) weight=1.198776;
		 if(x>= 1.950000 &&  x < 2.100000 ) weight=1.356289;
		 if(x>= 2.100000 &&  x < 2.250000 ) weight=1.554186;
		 if(x>= 2.250000 &&  x < 2.400000 ) weight=1.801106;
		 if(x>= 2.400000 &&  x < 2.550000 ) weight=2.115665;
		 if(x>= 2.550000 &&  x < 2.700000 ) weight=2.490248;
		 if(x>= 2.700000 &&  x < 2.850000 ) weight=2.919512;
		 if(x>= 2.850000 &&  x < 3.000000 ) weight=3.448719;
		 if(x>= 3.000000 &&  x < 3.150000 ) weight=4.085822;
		 if(x>= 3.150000 &&  x < 3.300000 ) weight=4.837239;
		 if(x>= 3.300000 &&  x < 3.450000 ) weight=5.719576;
		 if(x>= 3.450000 &&  x < 3.600000 ) weight=6.731485;
		 if(x>= 3.600000 &&  x < 3.750000 ) weight=7.974702;
		 if(x>= 3.750000 &&  x < 3.900000 ) weight=9.356885;
		 if(x>= 3.900000 &&  x < 4.050000 ) weight=10.977955;
		 if(x>= 4.050000 &&  x < 4.200000 ) weight=12.928445;
		 if(x>= 4.200000 &&  x < 4.350000 ) weight=15.102790;
		 if(x>= 4.350000 &&  x < 4.500000 ) weight=17.491980;
		 if(x>= 4.500000 &&  x < 4.650000 ) weight=20.551324;
		 if(x>= 4.650000 &&  x < 4.800000 ) weight=23.832445;
		 if(x>= 4.800000 &&  x < 4.950000 ) weight=27.640621;
		 if(x>= 4.950000 &&  x < 5.100000 ) weight=31.590454;
		 if(x>= 5.100000 &&  x < 5.250000 ) weight=36.638488;
		 if(x>= 5.250000 &&  x < 5.400000 ) weight=41.719615;
		 if(x>= 5.400000 &&  x < 5.550000 ) weight=47.984490;
		 if(x>= 5.550000 &&  x < 5.700000 ) weight=54.315247;
		 if(x>= 5.700000 &&  x < 5.850000 ) weight=62.798772;
		 if(x>= 5.850000 &&  x < 6.000000 ) weight=71.393159;
		 if(x>= 6.000000 &&  x < 6.150000 ) weight=80.622125;
		 if(x>= 6.150000 &&  x < 6.300000 ) weight=91.702054;
		 if(x>= 6.300000 &&  x < 6.450000 ) weight=103.029975;
		 if(x>= 6.450000 &&  x < 6.600000 ) weight=116.251936;
		 if(x>= 6.600000 &&  x < 6.750000 ) weight=130.500479;
		 if(x>= 6.750000 &&  x < 6.900000 ) weight=145.766943;
		 if(x>= 6.900000 &&  x < 7.050000 ) weight=165.069622;
		 if(x>= 7.050000 &&  x < 7.200000 ) weight=183.439451;
		 if(x>= 7.200000 &&  x < 7.350000 ) weight=203.642690;
		 if(x>= 7.350000 &&  x < 7.500000 ) weight=227.170526;
		 if(x>= 7.500000 &&  x < 7.650000 ) weight=255.379316;
		 if(x>= 7.650000 &&  x < 7.800000 ) weight=280.038994;
		 if(x>= 7.800000 &&  x < 7.950000 ) weight=311.676205;
		 if(x>= 7.950000 &&  x < 8.100000 ) weight=345.865988;
		 if(x>= 8.100000 &&  x < 8.250000 ) weight=382.607356;
		 if(x>= 8.250000 &&  x < 8.400000 ) weight=421.800557;
		 if(x>= 8.400000 &&  x < 8.550000 ) weight=462.123578;
		 if(x>= 8.550000 &&  x < 8.700000 ) weight=511.417290;
		 if(x>= 8.700000 &&  x < 8.850000 ) weight=560.539628;
		 if(x>= 8.850000 &&  x < 9.000000 ) weight=614.873142;
		 if(x>= 9.000000 &&  x < 9.150000 ) weight=675.757278;
		 if(x>= 9.150000 &&  x < 9.300000 ) weight=735.834297;
		 if(x>= 9.300000 &&  x < 9.450000 ) weight=809.291166;
		 if(x>= 9.450000 &&  x < 9.600000 ) weight=891.113928;
		 if(x>= 9.600000 &&  x < 9.750000 ) weight=962.050309;
		 if(x>= 9.750000 &&  x < 9.900000 ) weight=1040.726468;
		 if(x>= 9.900000 &&  x < 10.050000 ) weight=1132.417394;
		 if(x>= 10.050000 &&  x < 10.200000 ) weight=1237.242135;
		 if(x>= 10.200000 &&  x < 10.350000 ) weight=1342.500949;
		 if(x>= 10.350000 &&  x < 10.500000 ) weight=1459.701621;
		 if(x>= 10.500000 &&  x < 10.650000 ) weight=1580.949868;
		 if(x>= 10.650000 &&  x < 10.800000 ) weight=1699.710610;
		 if(x>= 10.800000 &&  x < 10.950000 ) weight=1841.103734;
		 if(x>= 10.950000 &&  x < 11.100000 ) weight=1988.412856;
		 if(x>= 11.100000 &&  x < 11.250000 ) weight=2155.130320;
		 if(x>= 11.250000 &&  x < 11.400000 ) weight=2311.939723;
		 if(x>= 11.400000 &&  x < 11.550000 ) weight=2509.365524;
		 if(x>= 11.550000 &&  x < 11.700000 ) weight=2666.786893;
		 if(x>= 11.700000 &&  x < 11.850000 ) weight=2882.932045;
		 if(x>= 11.850000 &&  x < 12.000000 ) weight=3098.431346;
		 if(x>= 12.000000 &&  x < 12.150000 ) weight=3318.461154;
		 if(x>= 12.150000 &&  x < 12.300000 ) weight=3561.833356;
		 if(x>= 12.300000 &&  x < 12.450000 ) weight=3824.711382;
		 if(x>= 12.450000 &&  x < 12.600000 ) weight=4124.000415;
		 if(x>= 12.600000 &&  x < 12.750000 ) weight=4365.575847;
		 if(x>= 12.750000 &&  x < 12.900000 ) weight=4719.798416;
		 if(x>= 12.900000 &&  x < 13.050000 ) weight=5020.193923;
		 if(x>= 13.050000 &&  x < 13.200000 ) weight=5372.197691;
		 if(x>= 13.200000 &&  x < 13.350000 ) weight=5699.544998;
		 if(x>= 13.350000 &&  x < 13.500000 ) weight=6134.392488;
		 if(x>= 13.500000 &&  x < 13.650000 ) weight=6487.706141;
		 if(x>= 13.650000 &&  x < 13.800000 ) weight=6903.037574;
		 if(x>= 13.800000 &&  x < 13.950000 ) weight=7427.179071;
		 if(x>= 13.950000 &&  x < 14.100000 ) weight=7819.318294;
		 if(x>= 14.100000 &&  x < 14.250000 ) weight=8438.071456;
		 if(x>= 14.250000 &&  x < 14.400000 ) weight=8789.399016;
		 if(x>= 14.400000 &&  x < 14.550000 ) weight=9461.991987;
		 if(x>= 14.550000 &&  x < 14.700000 ) weight=10051.828561;
		 if(x>= 14.700000 &&  x < 14.850000 ) weight=10584.915247;
		 if(x>= 14.850000 &&  x < 15.000000 ) weight=11267.159593;
		 if(x>= 15.000000 &&  x < 15.150000 ) weight=11962.630896;
		 if(x>= 15.150000 &&  x < 15.300000 ) weight=12621.696523;
		 if(x>= 15.300000 &&  x < 15.450000 ) weight=13427.227611;
		 if(x>= 15.450000 &&  x < 15.600000 ) weight=14131.633453;
		 if(x>= 15.600000 &&  x < 15.750000 ) weight=15104.205799;
		 if(x>= 15.750000 &&  x < 15.900000 ) weight=15927.992672;
		 if(x>= 15.900000 &&  x < 16.050000 ) weight=16835.557825;
		 if(x>= 16.050000 &&  x < 16.200000 ) weight=17557.688698;
		 if(x>= 16.200000 &&  x < 16.350000 ) weight=18604.082638;
		 if(x>= 16.350000 &&  x < 16.500000 ) weight=19551.366895;
		 if(x>= 16.500000 &&  x < 16.650000 ) weight=20825.192959;
		 if(x>= 16.650000 &&  x < 16.800000 ) weight=21848.545611;
		 if(x>= 16.800000 &&  x < 16.950000 ) weight=22946.463072;
		 if(x>= 16.950000 &&  x < 17.100000 ) weight=24183.065914;
		 if(x>= 17.100000 &&  x < 17.250000 ) weight=25690.376358;
		 if(x>= 17.250000 &&  x < 17.400000 ) weight=26979.187152;
		 if(x>= 17.400000 &&  x < 17.550000 ) weight=28376.784267;
		 if(x>= 17.550000 &&  x < 17.700000 ) weight=29774.573194;
		 if(x>= 17.700000 &&  x < 17.850000 ) weight=31316.698673;
		 if(x>= 17.850000 &&  x < 18.000000 ) weight=32628.808621;
		 if(x>= 18.000000 &&  x < 18.150000 ) weight=34532.814677;
		 if(x>= 18.150000 &&  x < 18.300000 ) weight=36639.596010;
		 if(x>= 18.300000 &&  x < 18.450000 ) weight=38054.399592;
		 if(x>= 18.450000 &&  x < 18.600000 ) weight=39898.590351;
		 if(x>= 18.600000 &&  x < 18.750000 ) weight=41850.088053;
		 if(x>= 18.750000 &&  x < 18.900000 ) weight=43862.632535;
		 if(x>= 18.900000 &&  x < 19.050000 ) weight=45702.669672;
		 if(x>= 19.050000 &&  x < 19.200000 ) weight=48289.448495;
		 if(x>= 19.200000 &&  x < 19.350000 ) weight=50708.267688;
		 if(x>= 19.350000 &&  x < 19.500000 ) weight=52594.793892;
		 if(x>= 19.500000 &&  x < 19.650000 ) weight=54868.226217;
		 if(x>= 19.650000 &&  x < 19.800000 ) weight=58019.566333;
		 if(x>= 19.800000 &&  x < 19.950000 ) weight=61000.537364;
		 if(x>= 19.950000 &&  x < 20.100000 ) weight=63139.129939;
		 if(x>= 20.100000 &&  x < 20.250000 ) weight=66123.099509;
		 if(x>= 20.250000 &&  x < 20.400000 ) weight=69788.492502;
		 if(x>= 20.400000 &&  x < 20.550000 ) weight=72151.247086;
		 if(x>= 20.550000 &&  x < 20.700000 ) weight=74783.421325;
		 if(x>= 20.700000 &&  x < 20.850000 ) weight=78498.445060;
		 if(x>= 20.850000 &&  x < 21.000000 ) weight=82994.656150;
		 if(x>= 21.000000 &&  x < 21.150000 ) weight=85986.635791;
		 if(x>= 21.150000 &&  x < 21.300000 ) weight=90219.603251;
		 if(x>= 21.300000 &&  x < 21.450000 ) weight=93191.413884;
		 if(x>= 21.450000 &&  x < 21.600000 ) weight=97705.202996;
		 if(x>= 21.600000 &&  x < 21.750000 ) weight=101399.924615;
		 if(x>= 21.750000 &&  x < 21.900000 ) weight=105515.724870;
		 if(x>= 21.900000 &&  x < 22.050000 ) weight=110352.370892;
		 if(x>= 22.050000 &&  x < 22.200000 ) weight=114820.565161;
		 if(x>= 22.200000 &&  x < 22.350000 ) weight=118395.101410;
		 if(x>= 22.350000 &&  x < 22.500000 ) weight=123984.958204;
		 if(x>= 22.500000 &&  x < 22.650000 ) weight=128744.001995;
		 if(x>= 22.650000 &&  x < 22.800000 ) weight=134532.759974;
		 if(x>= 22.800000 &&  x < 22.950000 ) weight=140415.005637;
		 if(x>= 22.950000 &&  x < 23.100000 ) weight=145332.324426;
		 if(x>= 23.100000 &&  x < 23.250000 ) weight=152464.808184;
		 if(x>= 23.250000 &&  x < 23.400000 ) weight=156397.875369;
		 if(x>= 23.400000 &&  x < 23.550000 ) weight=163383.278671;
		 if(x>= 23.550000 &&  x < 23.700000 ) weight=169408.483720;
		 if(x>= 23.700000 &&  x < 23.850000 ) weight=176358.580725;
		 if(x>= 23.850000 &&  x < 24.000000 ) weight=183755.458545;
		 if(x>= 24.000000 &&  x < 24.150000 ) weight=190759.434473;
		 if(x>= 24.150000 &&  x < 24.300000 ) weight=198138.960128;
		 if(x>= 24.300000 &&  x < 24.450000 ) weight=205671.252855;
		 if(x>= 24.450000 &&  x < 24.600000 ) weight=211954.498174;
		 if(x>= 24.600000 &&  x < 24.750000 ) weight=221988.260443;
		 if(x>= 24.750000 &&  x < 24.900000 ) weight=229269.115259;
		 if(x>= 24.900000 &&  x < 25.050000 ) weight=236272.412393;
		 if(x>= 25.050000 &&  x < 25.200000 ) weight=247536.035648;
		 if(x>= 25.200000 &&  x < 25.350000 ) weight=256745.661792;
		 if(x>= 25.350000 &&  x < 25.500000 ) weight=265959.407036;
		 if(x>= 25.500000 &&  x < 25.650000 ) weight=272792.205935;
		 if(x>= 25.650000 &&  x < 25.800000 ) weight=284958.879411;
		 if(x>= 25.800000 &&  x < 25.950000 ) weight=293135.933527;
		 if(x>= 25.950000 &&  x < 26.100000 ) weight=306789.669871;
		 if(x>= 26.100000 &&  x < 26.250000 ) weight=314656.642340;
		 if(x>= 26.250000 &&  x < 26.400000 ) weight=325491.535578;
		 if(x>= 26.400000 &&  x < 26.550000 ) weight=335557.709334;
		 if(x>= 26.550000 &&  x < 26.700000 ) weight=348618.967549;
		 if(x>= 26.700000 &&  x < 26.850000 ) weight=362826.927372;
		 if(x>= 26.850000 &&  x < 27.000000 ) weight=372229.699274;
		 if(x>= 27.000000 &&  x < 27.150000 ) weight=382681.114113;
		 if(x>= 27.150000 &&  x < 27.300000 ) weight=398779.963024;
		 if(x>= 27.300000 &&  x < 27.450000 ) weight=411694.789589;
		 if(x>= 27.450000 &&  x < 27.600000 ) weight=426843.753944;
		 if(x>= 27.600000 &&  x < 27.750000 ) weight=441188.173102;
		 if(x>= 27.750000 &&  x < 27.900000 ) weight=457830.728830;
		 if(x>= 27.900000 &&  x < 28.050000 ) weight=476308.351608;
		 if(x>= 28.050000 &&  x < 28.200000 ) weight=488209.120811;
		 if(x>= 28.200000 &&  x < 28.350000 ) weight=503660.894719;
		 if(x>= 28.350000 &&  x < 28.500000 ) weight=515678.671461;
		 if(x>= 28.500000 &&  x < 28.650000 ) weight=532079.018234;
		 if(x>= 28.650000 &&  x < 28.800000 ) weight=552077.714737;
		 if(x>= 28.800000 &&  x < 28.950000 ) weight=565908.432159;
		 if(x>= 28.950000 &&  x < 29.100000 ) weight=589745.260513;
		 if(x>= 29.100000 &&  x < 29.250000 ) weight=606388.982300;
		 if(x>= 29.250000 &&  x < 29.400000 ) weight=623132.245117;
		 if(x>= 29.400000 &&  x < 29.550000 ) weight=648147.859135;
		 if(x>= 29.550000 &&  x < 29.700000 ) weight=665900.618201;
		 if(x>= 29.700000 &&  x < 29.850000 ) weight=684239.760351;
		 if(x>= 29.850000 &&  x < 30.000000 ) weight=702028.128093;
		 if(x>= 30) weight=702028.128093;
		 
		 }
		 
		 */
        
		else weight=1;
        
    }//close fUseTrigger
	
    return weight;
	
}
Double_t AliAnalysisTaskEMCalHFEpA::SetEoverPCutPtDependentMC(Double_t pt)
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
	


