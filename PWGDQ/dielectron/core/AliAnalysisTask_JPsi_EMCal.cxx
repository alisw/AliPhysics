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
//      Task for J/psi analysis using EMCal and                       //
//      EMCal correction framework                                    //
//																	  //
//		v1.0														  //
//                                                                    //
//	    Authors 							                          //
//		                                                              //
//		Cristiane Jahnke		(cristiane.jahnke@cern.ch)            //
//      22 January, 2021 -> TPC calibrations for 2017 and 2018 data   //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"
#include "TVector3.h"
#include <TRandom3.h>
#include "TProfile.h"
#include "TProfile2D.h"

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
#include "AliHFEtools.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliSelectNonHFE.h"
#include "AliHFEpidTPC.h"
#include "AliAnalysisTask_JPsi_EMCal.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"
#include "AliESDHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
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
#include "AliGenEventHeader.h"
#include "AliTrackerBase.h"
#include "AliAODVZERO.h"
#include "AliAODTracklets.h"
#include "AliESDUtils.h"
#include "AliAnalysisUtils.h"

//______________________________________________________________________

//______________________________________________________________________
ClassImp(AliAnalysisTask_JPsi_EMCal)

//______________________________________________________________________
AliAnalysisTask_JPsi_EMCal::AliAnalysisTask_JPsi_EMCal(const char *name)
  : AliAnalysisTaskSE(name)

,fIsMC(0)
,fUseTender(kFALSE)
,fMultiAnalysis(kFALSE)
,fFill_ESparse(kFALSE)
,fFill_ESparseTPC(kFALSE)
,fFill_MSparse(kFALSE)
,fIs_TPC_calibration(kFALSE)

//to select events with high energy cluster (to mimic the trigger)
,fSelect_trigger_events1(kFALSE)
,fSelect_trigger_events2(kFALSE)

//new Tender organization, using global variables
,fTenderClusterName("caloClusters")//TMClustersBranch
,fTenderTrackName("tracks")//usedefault
//,fTenderClusterName("TMClustersBranch")//TMClustersBranch
//,fTenderTrackName("usedefault")//usedefault
,fTracks_tender(0)
,fCaloClusters_tender(0)

//Used in the function FindMother
,fIsHFE1(kFALSE)
,fIsHFE2(kFALSE)
,fIsNonHFE(kFALSE)
,fIsFromD(kFALSE)
,fIsFromB(kFALSE)
,fIsFromPi0(kFALSE)
,fIsFromEta(kFALSE)
,fIsFromGamma(kFALSE)

//General variables
,fESD(0)
,fAOD(0)
,fVevent(0)
,fOutputList(0)
,fListProfiles(0)
,fPidResponse(0)
//,fNonHFE(new AliSelectNonHFE())
,fIsAOD(kFALSE)

,fEMCEG1(kFALSE)
,fEMCEG2(kFALSE)

,fEMCEG1DG1(kFALSE)//to run both EMCal and DCal triggers together
,fEMCEG2DG2(kFALSE)

,fEMCDG1(kFALSE)
,fEMCDG2(kFALSE)

,fIsTrack1Emcal(kFALSE)
,fIsTrack1Dcal(kFALSE)
,fIsTrack2Emcal(kFALSE)
,fIsTrack2Dcal(kFALSE)


,fIsEMCalCls(kFALSE)
,fIsDCalCls(kFALSE)

//analysis cuts
,fVertexCut(10)
//track cuts
,fEtaCutMin(-0.9)
,fEtaCutMax(0.9)
,fPtCutMainEle(1)
,fPtCutPartner(1)
,fRejectKinkMother(kTRUE)
,fAODGlobalTracks(kFALSE)
,fTPCandITSrefit(kTRUE)
,fITSncls(2)
,fITSpixel(1)
,fTPCncls(85)
,fTPCnCrossedRows(70)
,fTPCnclsPID(85)
,fTPCchi2(4)
,fITSchi2(36)
,fDCAxyCut(1)
,fDCAzCut(3)

,fTPCnsigmaCutMin(-1.5)
,fTPCnsigmaCutMax(3)

,fEnergyCut(1)
,fEoverPCutMin(0.8)
,fEoverPCutMax(1.3)
,fMassCutMin(2.92)
,fMassCutMax(3.16)

,fZvtx(0)
//global variable for multiplicity analysis
,fV0Mult(0)
,fSPDMult(0)
,fV0Mult_corr(0)
,fV0Mult_corr2(0)
,fSPDMult_corr(0)
,fRefMult(12.00)
,gRandom(new TRandom3(1607260721))
//,gRandom(new TRandom3(0))

,fRefMult_V0(139.0)
,gRandom_V0(new TRandom3(1607260721))
//,gRandom_V0(new TRandom3(0))

,fClus(0)
,fClus2(0)
,fClusAOD(0)

//Histograms for the analysis
,fNevent(0)
,fNevent2(0)
,fTPC_vs_ITScls(0)
,fPDG_values(0)
,fNevent_SPD_multi(0)
,fNevent_V0_multi(0)
,fEoverP_pt(0)
,fTPC_p(0)
,fTPCnsigma_p(0)
,fTPCnsigma_p_beforeCalibration(0)
,fTPCnsigma_p_afterCalibration(0)
,fTOF_p(0)
,fTOFnsigma_p(0)
,fTPCnsigma_EoverP(0)
,fECluster(0)
,fECluster_emcal(0)
,fECluster_dcal(0)
,fTracksPt(0)
,fTracksQAPt(0)
,fTracksMCPt(0)
,fVtxZ(0)
//histos for SPD and V0 multiplicity
,fVtxZ_V0(0)
,fVtxZ_SPD(0)
,fV0_SPD(0)
,fV0_nch(0)
,fSPD_nch(0)

,fV0(0)
,fV01(0)
,fV02(0)

,fSPD(0)
,fSPD1(0)
,fSPD2(0)



,fNClusters(0)
,fNClusters_pure(0)

,fEoverP_ntracks_matched(0)
,fEoverP_ncells(0)

,fECluster_pure(0)
,fECluster_pure_emcal(0)

,fECluster_pure_emcal_SPD1(0)
,fECluster_pure_emcal_SPD2(0)
,fECluster_pure_emcal_SPD3(0)
,fECluster_pure_emcal_SPD4(0)
,fECluster_pure_emcal_SPD5(0)

,fECluster_pure_emcal_V01(0)
,fECluster_pure_emcal_V02(0)
,fECluster_pure_emcal_V03(0)
,fECluster_pure_emcal_V04(0)
,fECluster_pure_emcal_V05(0)


,fECluster_pure_dcal(0)

,fEtaPhi_both(0)
,fEtaPhi_emcal(0)
,fEtaPhi_dcal(0)

//For the HFE package
//,fCuts(0)
//,fCFM(0)
//,fPID(new AliHFEpid("hfePid"))
//,fPIDqa(0)


,fMCarray(0)
,fMCheader(0)
,fMCparticle(0)
,fMCparticleMother(0)


,fMCparticle2(0)
,fMCparticleMother2(0)


,fMCparticleGMother(0)
,fMCparticleGGMother(0)
,fMCparticleGGGMother(0)
,fEventHandler(0)
,fMCevent(0)

	//JPsi histos
//,fHist_InvMass_pt_ULS(0)
//,fHist_InvMass_pt_LS(0)

	//KF
,fHist_InvMass_pt_ULS_KF(0)
,fHist_InvMass_pt_LS_KF(0)

,fHist_Correlation_leg1_emcal_leg2_not(0)
,fHist_Correlation_leg1_not_leg2_emcal(0)
,fHist_Correlation_leg1_emcal_leg2_emcal(0)


//multiplicity histos
,fHist_InvMass_pt_ULS_KF_SPDmulti_1(0)
,fHist_InvMass_pt_ULS_KF_SPDmulti_2(0)
,fHist_InvMass_pt_ULS_KF_SPDmulti_3(0)
,fHist_InvMass_pt_ULS_KF_SPDmulti_4(0)
,fHist_InvMass_pt_ULS_KF_SPDmulti_5(0)

,fHist_InvMass_pt_ULS_KF_V0multi_1(0)
,fHist_InvMass_pt_ULS_KF_V0multi_2(0)
,fHist_InvMass_pt_ULS_KF_V0multi_3(0)
,fHist_InvMass_pt_ULS_KF_V0multi_4(0)
,fHist_InvMass_pt_ULS_KF_V0multi_5(0)

//with weight
//KF

,fHist_InvMass_pt_ULS_KF_weight(0)


//multiplicity histos
,fHist_InvMass_pt_ULS_KF_SPDmulti_1_weight(0)
,fHist_InvMass_pt_ULS_KF_SPDmulti_2_weight(0)
,fHist_InvMass_pt_ULS_KF_SPDmulti_3_weight(0)
,fHist_InvMass_pt_ULS_KF_SPDmulti_4_weight(0)
,fHist_InvMass_pt_ULS_KF_SPDmulti_5_weight(0)

,fHist_InvMass_pt_ULS_KF_V0multi_1_weight(0)
,fHist_InvMass_pt_ULS_KF_V0multi_2_weight(0)
,fHist_InvMass_pt_ULS_KF_V0multi_3_weight(0)
,fHist_InvMass_pt_ULS_KF_V0multi_4_weight(0)
,fHist_InvMass_pt_ULS_KF_V0multi_5_weight(0)



//generators
//BB
,fHist_InvMass_pt_ULS_KF_BB(0)
,fHist_InvMass_pt_LS_KF_BB(0)
//CC
,fHist_InvMass_pt_ULS_KF_CC(0)
,fHist_InvMass_pt_LS_KF_CC(0)
//B
,fHist_InvMass_pt_ULS_KF_B(0)
,fHist_InvMass_pt_LS_KF_B(0)
//JPsi
,fHist_InvMass_pt_ULS_KF_Jpsi(0)
,fHist_InvMass_pt_LS_KF_Jpsi(0)
//BJpsi
,fHist_InvMass_pt_ULS_KF_BJpsi(0)
,fHist_InvMass_pt_LS_KF_BJpsi(0)



	//leg 1 on EMCal
,fHist_InvMass_pt_ULS1(0)
,fHist_InvMass_pt_LS1(0)
	//leg 2 on EMCal
,fHist_InvMass_pt_ULS2(0)
,fHist_InvMass_pt_LS2(0)
	//both legs on EMCal
,fHist_InvMass_pt_ULSboth(0)
,fHist_InvMass_pt_LSboth(0)

,fHist_InvMass_pt_ULStpc(0)
,fHist_InvMass_pt_LStpc(0)

,fHist_InvMass_pt_ULStpc_wMatching(0)
,fHist_InvMass_pt_LStpc_wMatching(0)

	//new histos
,fdEta_dPhi(0)


,fSparseElectron(0)
,fSparseElectronTPC(0)
,fvalueElectron(0)
,fvalueElectronTPC(0)
,fSparseMulti(0)
,fvalueMulti(0)


	//MC efficiencies
,fPtMCparticleAllHfe1(0)
,fPtMCparticleRecoHfe1(0)
,fPtMCparticleAll_e_from_JPsi(0)
,fPtMCparticleAll_JPsi_pT(0)

,fPtMCparticleAll_e_from_JPsi_electron(0)
,fPtMCparticleAll_JPsi_pT_electron(0)
,fPtMCparticleAll_e_from_JPsi_positron(0)
,fPtMCparticleAll_JPsi_pT_positron(0)

,fPtMCparticleAll_electrons(0)
,fPtMCparticleAll_particles(0)


,fPtMCparticleAll_trueJPsi_pT(0)
,fPtMCparticleAll_trueJPsi_pT_weight(0)
,fPtMCparticleAll_trueJPsi_pT_weight_prompt(0)
,fPtMCparticleReco_e_from_JPsi(0)

//tracking efficiency
,fPtMCparticleReco_electrons(0)
,fPtMCparticleReco_electrons_no_gamma(0)
,fPtMCparticleReco_particles(0)
//TPC PID efficiency
,fPtMCparticle_TPCpid_e_from_JPsi(0)
,fPtMCparticle_TPCpid_electrons(0)
,fPtMCparticle_TPCpid_e_from_JPsi_num(0)
,fPtMCparticle_TPCpid_electrons_num(0)
//EMCal PID efficiency
,fPtMCparticle_EMCalpid_leg1(0)
,fPtMCparticle_EMCalpid_leg2(0)


,fPtMCparticle_EMCal_TM_e_from_JPsi(0)
,fPtMCparticle_EMCal_TM_electrons(0)
,fPtMCparticle_EMCalpid_leg1_e_from_JPsi(0)
,fPtMCparticle_EMCalpid_leg2_e_from_JPsi(0)
,fPtMCparticle_EMCalpid_both_leg1_e_from_JPsi(0)
,fPtMCparticle_EMCalpid_both_leg2_e_from_JPsi(0)
,fPtMCparticle_Total_JPsi_pT(0)


//J/Psi reco
,fPtMCparticle_JPsi(0)
,fPtMCparticle_JPsi_num(0)
//J/Psi mass cut
,fPtMCparticle_JPsi_mass(0)
,fPtMCparticle_JPsi_mass_num(0)


,fPtMCparticle_Total_e_from_JPsi(0)
,fPtMCparticle_Total_e_from_JPsi_sameMother(0)
,fPtMCparticle_TotalplusMass_e_from_JPsi(0)
,fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother(0)

,fPtMCparticle_TotalplusMass_JPsi_pT(0)
,fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother(0)
,fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight(0)
,fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight_prompt(0)

{
  //Named constructor
	
		//fvalueElectron = new Double_t[5];
    
    for(Int_t i=0; i<1; i++) fMultEstimatorAvg[i]=0;
     for(Int_t i=0; i<1; i++) fMultEstimatorV0[i]=0;
    
	
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(1, TH1I::Class());
  DefineOutput(1, TList::Class());
  //  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTask_JPsi_EMCal::AliAnalysisTask_JPsi_EMCal()
  : AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTask_JPsi_EMCal")

,fIsMC(0)
,fUseTender(kFALSE)
,fMultiAnalysis(kFALSE)
,fFill_ESparse(kFALSE)
,fFill_ESparseTPC(kFALSE)
,fFill_MSparse(kFALSE)
,fIs_TPC_calibration(kFALSE)

//to select events with high energy cluster (to mimic the trigger)
,fSelect_trigger_events1(kFALSE)
,fSelect_trigger_events2(kFALSE)

//new Tender organization, uisng global variables
,fTenderClusterName("caloClusters")
,fTenderTrackName("tracks")

//,fTenderClusterName("TMClustersBranch")//TMClustersBranch
//,fTenderTrackName("usedefault")//usedefault

,fTracks_tender(0)
,fCaloClusters_tender(0)


//Used in the function FindMother
,fIsHFE1(kFALSE)
,fIsHFE2(kFALSE)
,fIsNonHFE(kFALSE)
,fIsFromD(kFALSE)
,fIsFromB(kFALSE)
,fIsFromPi0(kFALSE)
,fIsFromEta(kFALSE)
,fIsFromGamma(kFALSE)

//General variables
,fESD(0)
,fAOD(0)
,fVevent(0)
,fOutputList(0)
,fListProfiles(0)
,fPidResponse(0)
//,fNonHFE(new AliSelectNonHFE())
,fIsAOD(kFALSE)

,fEMCEG1(kFALSE)
,fEMCEG2(kFALSE)

,fEMCEG1DG1(kFALSE)//to run both EMCal and DCal triggers together
,fEMCEG2DG2(kFALSE)

,fEMCDG1(kFALSE)
,fEMCDG2(kFALSE)

,fIsTrack1Emcal(kFALSE)
,fIsTrack1Dcal(kFALSE)
,fIsTrack2Emcal(kFALSE)
,fIsTrack2Dcal(kFALSE)

,fIsEMCalCls(kFALSE)
,fIsDCalCls(kFALSE)

//analysis cuts
,fVertexCut(10)
//track cuts
,fEtaCutMin(-0.9)
,fEtaCutMax(0.9)
,fPtCutMainEle(1)
,fPtCutPartner(1)
,fRejectKinkMother(kTRUE)
,fAODGlobalTracks(kFALSE)
,fTPCandITSrefit(kTRUE)
,fITSncls(2)
,fITSpixel(1)
,fTPCncls(85)
,fTPCnCrossedRows(70)
,fTPCnclsPID(85)
,fTPCchi2(4)
,fITSchi2(36)
,fDCAxyCut(1)
,fDCAzCut(3)

,fTPCnsigmaCutMin(-1.5)
,fTPCnsigmaCutMax(3)

,fEnergyCut(1)
,fEoverPCutMin(0.8)
,fEoverPCutMax(1.3)
,fMassCutMin(2.92)
,fMassCutMax(3.16)

,fZvtx(0)
//global variable for multiplicity analysis
,fV0Mult(0)
,fSPDMult(0)
,fV0Mult_corr(0)
,fV0Mult_corr2(0)
,fSPDMult_corr(0)
,fRefMult(12.00)
,gRandom(new TRandom3(1607260721))
//,gRandom(new TRandom3(0))

,fRefMult_V0(139.0)
,gRandom_V0(new TRandom3(1607260721))
//,gRandom_V0(new TRandom3(0))

,fClus(0)
,fClus2(0)
,fClusAOD(0)

//Histograms for the analysis
,fNevent(0)
,fNevent2(0)
,fTPC_vs_ITScls(0)
,fPDG_values(0)
,fNevent_SPD_multi(0)
,fNevent_V0_multi(0)
,fEoverP_pt(0)
,fTPC_p(0)
,fTPCnsigma_p(0)
,fTPCnsigma_p_beforeCalibration(0)
,fTPCnsigma_p_afterCalibration(0)
,fTOF_p(0)
,fTOFnsigma_p(0)
,fTPCnsigma_EoverP(0)
,fECluster(0)
,fECluster_emcal(0)
,fECluster_dcal(0)
,fECluster_pure(0)
,fECluster_pure_emcal(0)

,fECluster_pure_emcal_SPD1(0)
,fECluster_pure_emcal_SPD2(0)
,fECluster_pure_emcal_SPD3(0)
,fECluster_pure_emcal_SPD4(0)
,fECluster_pure_emcal_SPD5(0)

,fECluster_pure_emcal_V01(0)
,fECluster_pure_emcal_V02(0)
,fECluster_pure_emcal_V03(0)
,fECluster_pure_emcal_V04(0)
,fECluster_pure_emcal_V05(0)


,fECluster_pure_dcal(0)

,fEtaPhi_both(0)
,fEtaPhi_emcal(0)
,fEtaPhi_dcal(0)


,fTracksPt(0)
,fTracksQAPt(0)
,fTracksMCPt(0)
,fVtxZ(0)
//histos for SPD and V0 multiplicity
,fVtxZ_V0(0)
,fVtxZ_SPD(0)
,fV0_SPD(0)
,fV0_nch(0)
,fSPD_nch(0)

,fV0(0)
,fV01(0)
,fV02(0)

,fSPD(0)
,fSPD1(0)
,fSPD2(0)

,fNClusters(0)
,fNClusters_pure(0)

,fEoverP_ntracks_matched(0)
,fEoverP_ncells(0)

,fMCarray(0)
,fMCheader(0)
,fMCparticle(0)
,fMCparticleMother(0)


,fMCparticle2(0)
,fMCparticleMother2(0)

,fMCparticleGMother(0)
,fMCparticleGGMother(0)
,fMCparticleGGGMother(0)
,fEventHandler(0)
,fMCevent(0)


	//JPsi histos
//,fHist_InvMass_pt_ULS(0)
//,fHist_InvMass_pt_LS(0)

	//KF
,fHist_InvMass_pt_ULS_KF(0)
,fHist_InvMass_pt_LS_KF(0)

,fHist_Correlation_leg1_emcal_leg2_not(0)
,fHist_Correlation_leg1_not_leg2_emcal(0)
,fHist_Correlation_leg1_emcal_leg2_emcal(0)

//multiplicity histos
,fHist_InvMass_pt_ULS_KF_SPDmulti_1(0)
,fHist_InvMass_pt_ULS_KF_SPDmulti_2(0)
,fHist_InvMass_pt_ULS_KF_SPDmulti_3(0)
,fHist_InvMass_pt_ULS_KF_SPDmulti_4(0)
,fHist_InvMass_pt_ULS_KF_SPDmulti_5(0)

,fHist_InvMass_pt_ULS_KF_V0multi_1(0)
,fHist_InvMass_pt_ULS_KF_V0multi_2(0)
,fHist_InvMass_pt_ULS_KF_V0multi_3(0)
,fHist_InvMass_pt_ULS_KF_V0multi_4(0)
,fHist_InvMass_pt_ULS_KF_V0multi_5(0)


//with weight
//KF
,fHist_InvMass_pt_ULS_KF_weight(0)


//multiplicity histos
,fHist_InvMass_pt_ULS_KF_SPDmulti_1_weight(0)
,fHist_InvMass_pt_ULS_KF_SPDmulti_2_weight(0)
,fHist_InvMass_pt_ULS_KF_SPDmulti_3_weight(0)
,fHist_InvMass_pt_ULS_KF_SPDmulti_4_weight(0)
,fHist_InvMass_pt_ULS_KF_SPDmulti_5_weight(0)

,fHist_InvMass_pt_ULS_KF_V0multi_1_weight(0)
,fHist_InvMass_pt_ULS_KF_V0multi_2_weight(0)
,fHist_InvMass_pt_ULS_KF_V0multi_3_weight(0)
,fHist_InvMass_pt_ULS_KF_V0multi_4_weight(0)
,fHist_InvMass_pt_ULS_KF_V0multi_5_weight(0)


//generators
//BB
,fHist_InvMass_pt_ULS_KF_BB(0)
,fHist_InvMass_pt_LS_KF_BB(0)
//CC
,fHist_InvMass_pt_ULS_KF_CC(0)
,fHist_InvMass_pt_LS_KF_CC(0)
//B
,fHist_InvMass_pt_ULS_KF_B(0)
,fHist_InvMass_pt_LS_KF_B(0)
//JPsi
,fHist_InvMass_pt_ULS_KF_Jpsi(0)
,fHist_InvMass_pt_LS_KF_Jpsi(0)
//BJpsi
,fHist_InvMass_pt_ULS_KF_BJpsi(0)
,fHist_InvMass_pt_LS_KF_BJpsi(0)


	//leg 1 on EMCal
,fHist_InvMass_pt_ULS1(0)
,fHist_InvMass_pt_LS1(0)
	//leg 2 on EMCal
,fHist_InvMass_pt_ULS2(0)
,fHist_InvMass_pt_LS2(0)
	//both legs on EMCal
,fHist_InvMass_pt_ULSboth(0)
,fHist_InvMass_pt_LSboth(0)


,fHist_InvMass_pt_ULStpc(0)
,fHist_InvMass_pt_LStpc(0)

,fHist_InvMass_pt_ULStpc_wMatching(0)
,fHist_InvMass_pt_LStpc_wMatching(0)

	//new histos
,fdEta_dPhi(0)

,fSparseElectron(0)
,fSparseElectronTPC(0)
,fvalueElectron(0)
,fvalueElectronTPC(0)
,fSparseMulti(0)
,fvalueMulti(0)


	//MC efficiencies
,fPtMCparticleAllHfe1(0)
,fPtMCparticleRecoHfe1(0)
,fPtMCparticleAll_e_from_JPsi(0)
,fPtMCparticleAll_JPsi_pT(0)

,fPtMCparticleAll_e_from_JPsi_electron(0)
,fPtMCparticleAll_JPsi_pT_electron(0)
,fPtMCparticleAll_e_from_JPsi_positron(0)
,fPtMCparticleAll_JPsi_pT_positron(0)

,fPtMCparticleAll_electrons(0)
,fPtMCparticleAll_particles(0)

,fPtMCparticleAll_trueJPsi_pT(0)
,fPtMCparticleAll_trueJPsi_pT_weight(0)
,fPtMCparticleAll_trueJPsi_pT_weight_prompt(0)
,fPtMCparticleReco_e_from_JPsi(0)

,fPtMCparticleReco_electrons(0)
,fPtMCparticleReco_electrons_no_gamma(0)
,fPtMCparticleReco_particles(0)

//TPC PID efficiency
,fPtMCparticle_TPCpid_e_from_JPsi(0)
,fPtMCparticle_TPCpid_electrons(0)
,fPtMCparticle_TPCpid_e_from_JPsi_num(0)
,fPtMCparticle_TPCpid_electrons_num(0)
//EMCal PID efficiency
,fPtMCparticle_EMCalpid_leg1(0)
,fPtMCparticle_EMCalpid_leg2(0)

,fPtMCparticle_EMCal_TM_e_from_JPsi(0)
,fPtMCparticle_EMCal_TM_electrons(0)
,fPtMCparticle_EMCalpid_leg1_e_from_JPsi(0)
,fPtMCparticle_EMCalpid_leg2_e_from_JPsi(0)
,fPtMCparticle_EMCalpid_both_leg1_e_from_JPsi(0)
,fPtMCparticle_EMCalpid_both_leg2_e_from_JPsi(0)
,fPtMCparticle_Total_JPsi_pT(0)


//J/Psi reco
,fPtMCparticle_JPsi(0)
,fPtMCparticle_JPsi_num(0)
//J/Psi mass cut
,fPtMCparticle_JPsi_mass(0)
,fPtMCparticle_JPsi_mass_num(0)




,fPtMCparticle_Total_e_from_JPsi(0)
,fPtMCparticle_Total_e_from_JPsi_sameMother(0)
,fPtMCparticle_TotalplusMass_e_from_JPsi(0)
,fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother(0)

,fPtMCparticle_TotalplusMass_JPsi_pT(0)
,fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother(0)
,fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight(0)
,fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight_prompt(0)


{
	// Constructor
	
		//fvalueElectron = new Double_t[5];
	for(Int_t i=0; i<1; i++) fMultEstimatorAvg[i]=0;
     for(Int_t i=0; i<1; i++) fMultEstimatorV0[i]=0;
    
	// Define input and output slots here
	// Input slot #0 works with a TChain
	DefineInput(0, TChain::Class());
	// Output slot #0 id reserved by the base class for AOD
	// Output slot #1 writes into a TH1 container
	// DefineOutput(1, TH1I::Class());
	DefineOutput(1, TList::Class());
	//DefineOutput(3, TTree::Class());
}

//______________________________________________________________________
AliAnalysisTask_JPsi_EMCal::~AliAnalysisTask_JPsi_EMCal()
{
	//Destructor 
	delete fOutputList;
	//delete fPID;
	//delete fCFM;
	//delete fPIDqa;
	
	delete fSparseElectron;
    delete fSparseElectronTPC;
	delete []fvalueElectron;
    delete []fvalueElectronTPC;
    delete fSparseMulti;
    delete []fvalueMulti;
    //new organization of Tender with global variables
    delete fTracks_tender;
    delete fCaloClusters_tender;
    
    for(Int_t i=0; i<1; i++) {
        if (fMultEstimatorAvg[i]) delete fMultEstimatorAvg[i];
    }
    for(Int_t i=0; i<1; i++) {
        if (fMultEstimatorV0[i]) delete fMultEstimatorV0[i];
    }
    delete fListProfiles;
    delete gRandom;
    delete gRandom_V0;
    
}
//_____________________________________________________________________________
void AliAnalysisTask_JPsi_EMCal::Init()
{
    // Initialization of SPD corrections profiles!!!
    
    
    fListProfiles = new TList();
    fListProfiles->SetOwner();
  //  TString period[2];
   // Int_t nProfiles=2;
    TString period[1];
    Int_t nProfiles=1;
    Int_t nProfilesV0=1;
    period[0]="LHC16l";
   // period[1]="LHC16k";
    
    
    for(Int_t i=0; i<nProfiles; i++){
        if(fMultEstimatorAvg[i]){
            TProfile2D* hprof=new TProfile2D(*fMultEstimatorAvg[i]);
            hprof->SetName("ProfileSPD \n");
            fListProfiles->Add(hprof);
        }
    }
    for(Int_t i=0; i<nProfilesV0; i++){
        if(fMultEstimatorV0[i]){
            TProfile2D* hprofV0=new TProfile2D(*fMultEstimatorV0[i]);
            hprofV0->SetName("ProfileV0 \n");
            fListProfiles->Add(hprofV0);
            
        }
    }
    
   // PostData(2,fListProfiles);
    
    
    return;
}

//______________________________________________________________________
//Create Output Objects
//Here we can define the histograms and others output files
//Called once
void AliAnalysisTask_JPsi_EMCal::UserCreateOutputObjects()
{
///______________________________________________________________________
///Output Tlist
//Create TList
	fOutputList = new TList();
	fOutputList->SetOwner();
    
    
    gRandom->SetSeed(1607260721);


//Store the number of events
	//Define the histo
	fNevent = new TH1F("fNevent","Number of Events",30,-0.5,29.5);
    fNevent2 = new TH1F("fNevent2","Number of Events",20,-0.5,19.5);
    fPDG_values = new TH1F("fPDG_values","PDG of generated particles",6000,-3000,3000);
   
	//And then, add to the output list
	fOutputList->Add(fNevent);
    fOutputList->Add(fNevent2);
    fOutputList->Add(fPDG_values);
    
    fNevent_SPD_multi = new TH1F("fNevent_SPD_multi","Number of Events in SPD bins",10,-0.5,9.5);
    fOutputList->Add(fNevent_SPD_multi);
    
    fNevent_V0_multi = new TH1F("fNevent_V0_multi","Number of Events in V0 bins",10,-0.5,9.5);
    fOutputList->Add(fNevent_V0_multi);
    
    //pileup check
    fTPC_vs_ITScls= new TH2F *[4];
	
	//General Histograms
	
	//Steps
	//Step 1: Before Track cuts
	//Step 2: Before PID
	//Step 3: After PID
	
	fEoverP_pt = new TH2F *[3];
	fTPC_p = new TH2F *[3];
	fTPCnsigma_p = new TH2F *[3];
    //TOF
    fTOF_p = new TH2F *[3];
    fTOFnsigma_p = new TH2F *[3];
    
	fTPCnsigma_EoverP = new TH2F *[3];
	fECluster= new TH1F *[4];
   
	
	fECluster_emcal= new TH1F *[3];
	fECluster_dcal= new TH1F *[3];
	
	fVtxZ= new  TH1F *[3];
	fNClusters= new TH1F *[3];

	fdEta_dPhi = new TH2F *[3];
	
	for(Int_t i = 0; i < 3; i++)
	{
	  fEoverP_pt[i] = new TH2F(Form("fEoverP_pt%d",i),";p_{t} (GeV/c);E / p ",400,0,40,500,0,2);
	  fTPC_p[i] = new TH2F(Form("fTPC_p%d",i),";p (GeV/c);TPC dE/dx (a. u.)",400,0,40,1000,-20,200);
	  fTPCnsigma_p[i] = new TH2F(Form("fTPCnsigma_p%d",i),";p (GeV/c);TPC Electron N#sigma",400,0,40,1000,-15,10);
        
        fTOF_p[i] = new TH2F(Form("fTOF_p%d",i),";p (GeV/c);TOF (a. u.)",400,0,40,300,-1,2);
        fTOFnsigma_p[i] = new TH2F(Form("fTOFnsigma_p%d",i),";p (GeV/c);TOF Electron N#sigma",400,0,40,1000,-15,10);
	  
		
		fECluster_emcal[i]= new TH1F(Form("fECluster_emcal%d",i), ";ECluster EMCal",1000, 0,100);
		fECluster_dcal[i]= new TH1F(Form("fECluster_dcal%d",i), ";ECluster DCal",1000, 0,100);
		
	  fVtxZ[i]= new  TH1F(Form("fVtxZ%d",i),"VtxZ",1000, -50,50);
	  fNClusters[i]= new TH1F(Form("fNClusters%d",i),"fNClusters0",100, 0,100);
		
	  fTPCnsigma_EoverP[i] = new TH2F(Form("fTPCnsigma_EoverP%d",i),";TPC Electron N#sigma; E/p",1600,-20,20,200,0,2);
		
			//new histos
	  fdEta_dPhi[i] = new TH2F(Form("fdEta_dPhi%d",i),"Distance of EMCAL cluster to its closest track ;#phi;z",100,-0.3,0.3,100,-0.3,0.3);

	  		
	  fOutputList->Add(fEoverP_pt[i]);
	  fOutputList->Add(fTPC_p[i]);
	  fOutputList->Add(fTPCnsigma_p[i]);
        
      fOutputList->Add(fTOF_p[i]);
      fOutputList->Add(fTOFnsigma_p[i]);
        
      fOutputList->Add(fTPCnsigma_EoverP[i]);
	  
		
	  fOutputList->Add(fECluster_emcal[i]);
	  fOutputList->Add(fECluster_dcal[i]);
		
		
	  fOutputList->Add(fVtxZ[i]);
	  fOutputList->Add(fNClusters[i]);
		
      //new histos
	  fOutputList->Add(fdEta_dPhi[i]);

	}
    
    fTPCnsigma_p_beforeCalibration = new TH2F("fTPCnsigma_p_beforeCalibration",";p (GeV/c);TPC Electron N#sigma (bef. calibration)",400,0,40,1000,-15,10);
    fTPCnsigma_p_afterCalibration = new TH2F("fTPCnsigma_p_afterCalibration", ";p (GeV/c);TPC Electron N#sigma (aft. calibration)",400,0,40,1000,-15,10);
    
    fOutputList->Add(fTPCnsigma_p_beforeCalibration);
    fOutputList->Add(fTPCnsigma_p_afterCalibration);
    
    
    
    for(Int_t i = 0; i < 4; i++)
    {
        fECluster[i]= new TH1F(Form("fECluster%d",i), ";ECluster",2000, 0,100);
        fOutputList->Add(fECluster[i]);
        
        //pileup histos
        fTPC_vs_ITScls[i]= new TH2F(Form("fTPC_vs_ITScls%d",i), ";# TPC clusters; #SSD and SDD clusters",10000, 0,10000, 2000, 0, 2000);
        fOutputList->Add(fTPC_vs_ITScls[i]);
    }
    
    //=================================================================================================================================================================
    // Multiplicity histos
    
    fNClusters_pure= new TH1F("fNClusters_pure","fNClusters_pure",100, 0,100);
    fOutputList->Add(fNClusters_pure);
    
    fEoverP_ntracks_matched = new TH2F("fEoverP_ntracks_matched","fEoverP_ntracks_matched;E/p; N tracks matched to a cluster",200,0,2,20,0,20);
    fOutputList->Add(fEoverP_ntracks_matched);
  
    fEoverP_ncells = new TH2F("fEoverP_ncells","fEoverP_ncells;E/p; N cells",200,0,2,100,0,100);
    fOutputList->Add(fEoverP_ncells);
    
 
    
    
    
    fVtxZ_V0 = new TH2F("fVtxZ_V0","V0 multi vs. VtxZ ;VtxZ; V0 multiplicity",400,-20,20,500,0,1000);
    fOutputList->Add(fVtxZ_V0);
    
    fVtxZ_SPD = new TH2F("fVtxZ_SPD","SPD multi vs. VtxZ ;VtxZ; SPD multiplicity",400,-20,20,500,0,500);
    fOutputList->Add(fVtxZ_SPD);
    
    fV0_SPD = new TH2F("fV0_SPD","SPD multi vs. V0 ;V0; SPD multiplicity",250,0,500,50,0,100);
    fOutputList->Add(fV0_SPD);
    
    fV0_nch = new TH2F("fV0_nch","V0 ;nch;V0 multiplicity",500, 0,1000,500,0,1000);
    fOutputList->Add(fV0_nch);
    
    fSPD_nch = new TH2F("fSPD_nch","SPD ;nch;SPD multiplicity",500, 0,1000,500,0,1000);
    fOutputList->Add(fSPD_nch);
    
    
    //multi histos with pT cut
    fV0 = new TH1F("fV0","V0 multi;V0 multiplicity",100,0,1000);
    fOutputList->Add(fV0);
    fV01 = new TH1F("fV01","V0 multi;V0 multiplicity EG1",100,0,1000);
    fOutputList->Add(fV01);
    fV02 = new TH1F("fV02","V0 multi;V0 multiplicity EG2",100,0,1000);
    fOutputList->Add(fV02);
    
    fSPD = new TH1F("fSPD","SPD multi;SPD multiplicity",100,0,100);
    fOutputList->Add(fSPD);
    fSPD1 = new TH1F("fSPD1","SPD multi;SPD multiplicity EG1",100,0,100);
    fOutputList->Add(fSPD1);
    fSPD2 = new TH1F("fSPD2","SPD multi;SPD multiplicity EG2",100,0,100);
    fOutputList->Add(fSPD2);
    
    
    
    
    
	//=================================================================================================================================================================
    
	fECluster_pure= new TH1F("fECluster_pure", ";ECluster pure",2000,0,100);
	fOutputList->Add(fECluster_pure);
	
    //emcal and dcal separated
	fECluster_pure_emcal= new TH1F("fECluster_pure_emcal", ";ECluster pure EMCal",2000,0,100);
    
    if(fMultiAnalysis){
        fECluster_pure_emcal_SPD1= new TH1F("fECluster_pure_emcal_SPD1", ";ECluster pure EMCal",2000,0,100);
        fECluster_pure_emcal_SPD2= new TH1F("fECluster_pure_emcal_SPD2", ";ECluster pure EMCal",2000,0,100);
        fECluster_pure_emcal_SPD3= new TH1F("fECluster_pure_emcal_SPD3", ";ECluster pure EMCal",2000,0,100);
        fECluster_pure_emcal_SPD4= new TH1F("fECluster_pure_emcal_SPD4", ";ECluster pure EMCal",2000,0,100);
        fECluster_pure_emcal_SPD5= new TH1F("fECluster_pure_emcal_SPD5", ";ECluster pure EMCal",2000,0,100);
    
        fECluster_pure_emcal_V01= new TH1F("fECluster_pure_emcal_V01", ";ECluster pure EMCal",2000,0,100);
        fECluster_pure_emcal_V02= new TH1F("fECluster_pure_emcal_V02", ";ECluster pure EMCal",2000,0,100);
        fECluster_pure_emcal_V03= new TH1F("fECluster_pure_emcal_V03", ";ECluster pure EMCal",2000,0,100);
        fECluster_pure_emcal_V04= new TH1F("fECluster_pure_emcal_V04", ";ECluster pure EMCal",2000,0,100);
        fECluster_pure_emcal_V05= new TH1F("fECluster_pure_emcal_V05", ";ECluster pure EMCal",2000,0,100);
    }
    
	fOutputList->Add(fECluster_pure_emcal);
    
    if(fMultiAnalysis){
        fOutputList->Add(fECluster_pure_emcal_SPD1);
        fOutputList->Add(fECluster_pure_emcal_SPD2);
        fOutputList->Add(fECluster_pure_emcal_SPD3);
        fOutputList->Add(fECluster_pure_emcal_SPD4);
        fOutputList->Add(fECluster_pure_emcal_SPD5);
    
        fOutputList->Add(fECluster_pure_emcal_V01);
        fOutputList->Add(fECluster_pure_emcal_V02);
        fOutputList->Add(fECluster_pure_emcal_V03);
        fOutputList->Add(fECluster_pure_emcal_V04);
        fOutputList->Add(fECluster_pure_emcal_V05);
    }
    
	fECluster_pure_dcal= new TH1F("fECluster_pure_dcal", ";ECluster pure DCal",2000,0,100);
	fOutputList->Add(fECluster_pure_dcal);
	
	
	fEtaPhi_both= new TH2F("fEtaPhi_both","#eta x #phi Clusters;#phi;#eta",300,0.,6,200,-1.,1.);
	fEtaPhi_emcal= new TH2F("fEtaPhi_emcal","#eta x #phi Clusters EMCal;#phi;#eta",300,0.,6,200,-1.,1.);
	fEtaPhi_dcal= new TH2F("fEtaPhi_dcal","#eta x #phi Clusters DCal;#phi;#eta",300,0.,6,200,-1.,1.);
	
	fOutputList->Add(fEtaPhi_both);
	fOutputList->Add(fEtaPhi_emcal);
	fOutputList->Add(fEtaPhi_dcal);
	
	
	fTracksPt=new TH1F *[12];
	for(Int_t i=0; i<12; i++){
		fTracksPt[i]= new TH1F(Form("fTracksPt%d", i), ";p_{T} (GeV/c); Counts ", 300, 0, 30);
		fOutputList->Add(fTracksPt[i]);
	}
	
	fTracksQAPt=new TH1F *[12];
	for(Int_t i=0; i<12; i++){
		fTracksQAPt[i]= new TH1F(Form("fTracksQAPt%d", i), ";p_{T} (GeV/c); Counts ", 300, 0, 30);
		fOutputList->Add(fTracksQAPt[i]);
	}
    
    fTracksMCPt=new TH1F *[12];
    for(Int_t i=0; i<12; i++){
        fTracksMCPt[i]= new TH1F(Form("fTracksMCPt%d", i), ";p_{T} (GeV/c); Counts ", 400, 0, 40);
        fOutputList->Add(fTracksMCPt[i]);
    }
	
	
	
	//JPsi analysis histograms
    /*
	fHist_InvMass_pt_ULS = new TH2F("fHist_InvMass_pt_ULS","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_ULS);
	fHist_InvMass_pt_LS = new TH2F("fHist_InvMass_pt_LS","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_LS);
    */
	
    //KFParticle
	fHist_InvMass_pt_ULS_KF = new TH2F("fHist_InvMass_pt_ULS_KF","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",50,0,50,1000,0,10);
	fOutputList->Add(fHist_InvMass_pt_ULS_KF);
	fHist_InvMass_pt_LS_KF = new TH2F("fHist_InvMass_pt_LS_KF","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",50,0,50,1000,0,10);
	fOutputList->Add(fHist_InvMass_pt_LS_KF);
    
    //Correlation btween leg1 and leg2
    fHist_Correlation_leg1_emcal_leg2_not = new TH2F("fHist_Correlation_leg1_emcal_leg2_not","leg1 vs leg2;p_{T} leg1 (GeV/c);p_{T} leg2 (GeV/c)",50,0,50,50,0,50);
    fHist_Correlation_leg1_not_leg2_emcal = new TH2F("fHist_Correlation_leg1_not_leg2_emcal","leg1 vs leg2;p_{T} leg1 (GeV/c);p_{T} leg2 (GeV/c)",50,0,50,50,0,50);
    fHist_Correlation_leg1_emcal_leg2_emcal = new TH2F("fHist_Correlation_leg1_emcal_leg2_emcal","leg1 vs leg2;p_{T} leg1 (GeV/c);p_{T} leg2 (GeV/c)",50,0,50,50,0,50);
    
    fOutputList->Add(fHist_Correlation_leg1_emcal_leg2_not);
    fOutputList->Add(fHist_Correlation_leg1_not_leg2_emcal);
    fOutputList->Add(fHist_Correlation_leg1_emcal_leg2_emcal);
    
    
    
    //multiplicity histos
    
     if(fMultiAnalysis){
         fHist_InvMass_pt_ULS_KF_SPDmulti_1 = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_1","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fHist_InvMass_pt_ULS_KF_SPDmulti_2 = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_2","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fHist_InvMass_pt_ULS_KF_SPDmulti_3 = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_3","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fHist_InvMass_pt_ULS_KF_SPDmulti_4 = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_4","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fHist_InvMass_pt_ULS_KF_SPDmulti_5 = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_5","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
     
    
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_1);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_2);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_3);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_4);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_5);
    
         fHist_InvMass_pt_ULS_KF_V0multi_1 = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_1","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fHist_InvMass_pt_ULS_KF_V0multi_2 = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_2","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fHist_InvMass_pt_ULS_KF_V0multi_3 = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_3","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fHist_InvMass_pt_ULS_KF_V0multi_4 = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_4","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fHist_InvMass_pt_ULS_KF_V0multi_5 = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_5","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);

         fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_1);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_2);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_3);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_4);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_5);
 
         //multiplicity histos with weight
         //multiplicity histos
    
         //KFParticle
    
         fHist_InvMass_pt_ULS_KF_weight = new TH2F("fHist_InvMass_pt_ULS_KF_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_weight);
    
    
         fHist_InvMass_pt_ULS_KF_SPDmulti_1_weight = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_1_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fHist_InvMass_pt_ULS_KF_SPDmulti_2_weight = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_2_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fHist_InvMass_pt_ULS_KF_SPDmulti_3_weight = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_3_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fHist_InvMass_pt_ULS_KF_SPDmulti_4_weight = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_4_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fHist_InvMass_pt_ULS_KF_SPDmulti_5_weight = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_5_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
    
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_1_weight);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_2_weight);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_3_weight);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_4_weight);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_5_weight);
    
         fHist_InvMass_pt_ULS_KF_V0multi_1_weight = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_1_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fHist_InvMass_pt_ULS_KF_V0multi_2_weight = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_2_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fHist_InvMass_pt_ULS_KF_V0multi_3_weight = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_3_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fHist_InvMass_pt_ULS_KF_V0multi_4_weight = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_4_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
         fHist_InvMass_pt_ULS_KF_V0multi_5_weight = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_5_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
    
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_1_weight);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_2_weight);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_3_weight);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_4_weight);
         fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_5_weight);

    
     }
	
	//=================================================================================================================================================================
	//MC generators
    
    if(fIsMC){
        //BB
        fHist_InvMass_pt_ULS_KF_BB = new TH2F("fHist_InvMass_pt_ULS_KF_BB","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
        fOutputList->Add(fHist_InvMass_pt_ULS_KF_BB);
        fHist_InvMass_pt_LS_KF_BB = new TH2F("fHist_InvMass_pt_LS_KF_BB","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",500,0,50,500,0,5);
        fOutputList->Add(fHist_InvMass_pt_LS_KF_BB);
	
        //CC
        fHist_InvMass_pt_ULS_KF_CC = new TH2F("fHist_InvMass_pt_ULS_KF_CC","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
        fOutputList->Add(fHist_InvMass_pt_ULS_KF_CC);
        fHist_InvMass_pt_LS_KF_CC = new TH2F("fHist_InvMass_pt_LS_KF_CC","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",500,0,50,500,0,5);
        fOutputList->Add(fHist_InvMass_pt_LS_KF_CC);
	
        //B
        fHist_InvMass_pt_ULS_KF_B = new TH2F("fHist_InvMass_pt_ULS_KF_B","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
        fOutputList->Add(fHist_InvMass_pt_ULS_KF_B);
        fHist_InvMass_pt_LS_KF_B = new TH2F("fHist_InvMass_pt_LS_KF_B","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",500,0,50,500,0,5);
        fOutputList->Add(fHist_InvMass_pt_LS_KF_B);
	
        //Jpsi
        fHist_InvMass_pt_ULS_KF_Jpsi = new TH2F("fHist_InvMass_pt_ULS_KF_Jpsi","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
        fOutputList->Add(fHist_InvMass_pt_ULS_KF_Jpsi);
        fHist_InvMass_pt_LS_KF_Jpsi = new TH2F("fHist_InvMass_pt_LS_KF_Jpsi","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",500,0,50,500,0,5);
        fOutputList->Add(fHist_InvMass_pt_LS_KF_Jpsi);
	
        //BJpsi
        fHist_InvMass_pt_ULS_KF_BJpsi = new TH2F("fHist_InvMass_pt_ULS_KF_BJpsi","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",500,0,50,500,0,5);
        fOutputList->Add(fHist_InvMass_pt_ULS_KF_BJpsi);
        fHist_InvMass_pt_LS_KF_BJpsi = new TH2F("fHist_InvMass_pt_LS_KF_BJpsi","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",500,0,50,500,0,5);
        fOutputList->Add(fHist_InvMass_pt_LS_KF_BJpsi);
        
    }
	
	//=================================================================================================================================================================

	
	fHist_InvMass_pt_ULS1 = new TH2F("fHist_InvMass_pt_ULS1","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",50,0,50,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_ULS1);
	fHist_InvMass_pt_LS1 = new TH2F("fHist_InvMass_pt_LS1","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",50,0,50,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_LS1);
	
	fHist_InvMass_pt_ULS2 = new TH2F("fHist_InvMass_pt_ULS2","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",50,0,50,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_ULS2);
	fHist_InvMass_pt_LS2 = new TH2F("fHist_InvMass_pt_LS2","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",50,0,50,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_LS2);
	
	fHist_InvMass_pt_ULSboth = new TH2F("fHist_InvMass_pt_ULSboth","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",50,0,50,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_ULSboth);
	fHist_InvMass_pt_LSboth = new TH2F("fHist_InvMass_pt_LSboth","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",50,0,50,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_LSboth);
	
	
	
	fHist_InvMass_pt_ULStpc = new TH2F("fHist_InvMass_pt_ULStpc","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",50,0,50,100,0,10);
	fOutputList->Add(fHist_InvMass_pt_ULStpc);
	
	fHist_InvMass_pt_LStpc = new TH2F("fHist_InvMass_pt_LStpc","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",50,0,50,100,0,10);
	fOutputList->Add(fHist_InvMass_pt_LStpc);
    
    
    fHist_InvMass_pt_ULStpc_wMatching = new TH2F("fHist_InvMass_pt_ULStpc_wMatching","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",50,0,50,1000,0,10);
    fOutputList->Add(fHist_InvMass_pt_ULStpc_wMatching);
    
    fHist_InvMass_pt_LStpc_wMatching = new TH2F("fHist_InvMass_pt_LStpc_wMatching","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",50,0,50,1000,0,10);
    fOutputList->Add(fHist_InvMass_pt_LStpc_wMatching);
	
    
    
    //MC efficiencies
    if(fIsMC){
        
        fPtMCparticleRecoHfe1 = new TH1F("fPtMCparticleRecoHfe1",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticleAllHfe1 = new TH1F("fPtMCparticleAllHfe1",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticleAll_e_from_JPsi = new TH1F("fPtMCparticleAll_e_from_JPsi",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticleAll_JPsi_pT = new TH1F("fPtMCparticleAll_JPsi_pT",";p_{T} (GeV/c);Count",500,0,50);
        
        //new
        fPtMCparticleAll_e_from_JPsi_electron = new TH1F("fPtMCparticleAll_e_from_JPsi_electron",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticleAll_JPsi_pT_electron = new TH1F("fPtMCparticleAll_JPsi_pT_electron",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticleAll_e_from_JPsi_positron = new TH1F("fPtMCparticleAll_e_from_JPsi_positron",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticleAll_JPsi_pT_positron = new TH1F("fPtMCparticleAll_JPsi_pT_positron",";p_{T} (GeV/c);Count",500,0,50);
        
        //denominator tracking efficiency
        fPtMCparticleAll_electrons = new TH1F("fPtMCparticleAll_electrons",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticleAll_particles = new TH1F("fPtMCparticleAll_particles",";p_{T} (GeV/c);Count",500,0,50);
        
        
        fPtMCparticleAll_trueJPsi_pT = new TH1F("fPtMCparticleAll_trueJPsi_pT",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticleAll_trueJPsi_pT_weight = new TH1F("fPtMCparticleAll_trueJPsi_pT_weight",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticleAll_trueJPsi_pT_weight_prompt = new TH1F("fPtMCparticleAll_trueJPsi_pT_weight_prompt",";p_{T} (GeV/c);Count",500,0,50);
        
        
        fPtMCparticleReco_e_from_JPsi = new TH1F("fPtMCparticleReco_e_from_JPsi",";p_{T} (GeV/c);Count",500,0,50);
        
        fPtMCparticleReco_electrons = new TH1F("fPtMCparticleReco_electrons",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticleReco_electrons_no_gamma = new TH1F("fPtMCparticleReco_electrons_no_gamma",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticleReco_particles = new TH1F("fPtMCparticleReco_particles",";p_{T} (GeV/c);Count",500,0,50);
        
        
        //TPC PID efficiency
        fPtMCparticle_TPCpid_e_from_JPsi = new TH1F("fPtMCparticle_TPCpid_e_from_JPsi",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_TPCpid_electrons = new TH1F("fPtMCparticle_TPCpid_electrons",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_TPCpid_e_from_JPsi_num = new TH1F("fPtMCparticle_TPCpid_e_from_JPsi_num",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_TPCpid_electrons_num = new TH1F("fPtMCparticle_TPCpid_electrons_num",";p_{T} (GeV/c);Count",500,0,50);
        //EMCal PID efficiency
        fPtMCparticle_EMCalpid_leg1 = new TH1F("fPtMCparticle_EMCalpid_leg1",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_EMCalpid_leg2 = new TH1F("fPtMCparticle_EMCalpid_leg2",";p_{T} (GeV/c);Count",500,0,50);
        
        
        fPtMCparticle_EMCal_TM_e_from_JPsi = new TH1F("fPtMCparticle_EMCal_TM_e_from_JPsi",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_EMCal_TM_electrons = new TH1F("fPtMCparticle_EMCal_TM_electrons",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_EMCalpid_leg1_e_from_JPsi = new TH1F("fPtMCparticle_EMCalpid_leg1_e_from_JPsi",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_EMCalpid_leg2_e_from_JPsi = new TH1F("fPtMCparticle_EMCalpid_leg2_e_from_JPsi",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_EMCalpid_both_leg1_e_from_JPsi = new TH1F("fPtMCparticle_EMCalpid_both_leg1_e_from_JPsi",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_EMCalpid_both_leg2_e_from_JPsi = new TH1F("fPtMCparticle_EMCalpid_both_leg2_e_from_JPsi",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_Total_JPsi_pT = new TH1F("fPtMCparticle_Total_JPsi_pT",";p_{T} (GeV/c);Count",500,0,50);
        
        
        
        
        //J/Psi reco
        fPtMCparticle_JPsi = new TH1F("fPtMCparticle_JPsi",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_JPsi_num = new TH1F("fPtMCparticle_JPsi_num",";p_{T} (GeV/c);Count",500,0,50);
        //J/Psi mass cut
        fPtMCparticle_JPsi_mass = new TH1F("fPtMCparticle_JPsi_mass",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_JPsi_mass_num = new TH1F("fPtMCparticle_JPsi_mass_num",";p_{T} (GeV/c);Count",500,0,50);
        
        
	
        fPtMCparticle_Total_e_from_JPsi = new TH1F("fPtMCparticle_Total_e_from_JPsi",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_Total_e_from_JPsi_sameMother = new TH1F("fPtMCparticle_Total_e_from_JPsi_sameMother",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_TotalplusMass_e_from_JPsi = new TH1F("fPtMCparticle_TotalplusMass_e_from_JPsi",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother = new TH1F("fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother",";p_{T} (GeV/c);Count",500,0,50);
   
        fPtMCparticle_TotalplusMass_JPsi_pT = new TH1F("fPtMCparticle_TotalplusMass_JPsi_pT",";p_{T} (GeV/c);Count",500,0,50);
       
         fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother = new TH1F("fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight = new TH1F("fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight",";p_{T} (GeV/c);Count",500,0,50);
        fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight_prompt = new TH1F("fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight_prompt",";p_{T} (GeV/c);Count",500,0,50);
    
	
        fOutputList->Add(fPtMCparticleRecoHfe1);
        fOutputList->Add(fPtMCparticleAllHfe1);
        fOutputList->Add(fPtMCparticleAll_e_from_JPsi);
        fOutputList->Add(fPtMCparticleAll_JPsi_pT);
        
        fOutputList->Add(fPtMCparticleAll_e_from_JPsi_electron);
        fOutputList->Add(fPtMCparticleAll_JPsi_pT_electron);
        
        fOutputList->Add(fPtMCparticleAll_e_from_JPsi_positron);
        fOutputList->Add(fPtMCparticleAll_JPsi_pT_positron);
        
        fOutputList->Add(fPtMCparticleAll_electrons);
        fOutputList->Add(fPtMCparticleAll_particles);
        
        
        fOutputList->Add(fPtMCparticleAll_trueJPsi_pT);
        fOutputList->Add(fPtMCparticleAll_trueJPsi_pT_weight);
        fOutputList->Add(fPtMCparticleAll_trueJPsi_pT_weight_prompt);
        
        fOutputList->Add(fPtMCparticleReco_e_from_JPsi);
        
        fOutputList->Add(fPtMCparticleReco_electrons);
        fOutputList->Add(fPtMCparticleReco_electrons_no_gamma);
        fOutputList->Add(fPtMCparticleReco_particles);
        
        //TPC PID efficiency
        fOutputList->Add(fPtMCparticle_TPCpid_e_from_JPsi);
        fOutputList->Add(fPtMCparticle_TPCpid_electrons);
        fOutputList->Add(fPtMCparticle_TPCpid_e_from_JPsi_num);
        fOutputList->Add(fPtMCparticle_TPCpid_electrons_num);
        //EMCal PID efficiency
        fOutputList->Add(fPtMCparticle_EMCalpid_leg1);
        fOutputList->Add(fPtMCparticle_EMCalpid_leg2);
        
        fOutputList->Add(fPtMCparticle_EMCal_TM_e_from_JPsi);
        fOutputList->Add(fPtMCparticle_EMCal_TM_electrons);
        fOutputList->Add(fPtMCparticle_EMCalpid_leg1_e_from_JPsi);
        fOutputList->Add(fPtMCparticle_EMCalpid_leg2_e_from_JPsi);
        fOutputList->Add(fPtMCparticle_EMCalpid_both_leg1_e_from_JPsi);
        fOutputList->Add(fPtMCparticle_EMCalpid_both_leg2_e_from_JPsi);
        fOutputList->Add(fPtMCparticle_Total_JPsi_pT);
        
        
        //J/Psi reco
        fOutputList->Add(fPtMCparticle_JPsi);
        fOutputList->Add(fPtMCparticle_JPsi_num);
        //J/Psi mass cut
        fOutputList->Add(fPtMCparticle_JPsi_mass);
        fOutputList->Add(fPtMCparticle_JPsi_mass_num);
        
        
        fOutputList->Add(fPtMCparticle_Total_e_from_JPsi);
        fOutputList->Add(fPtMCparticle_Total_e_from_JPsi_sameMother);
        fOutputList->Add(fPtMCparticle_TotalplusMass_e_from_JPsi);
        fOutputList->Add(fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother);
    
        fOutputList->Add(fPtMCparticle_TotalplusMass_JPsi_pT);
        fOutputList->Add(fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother);
        fOutputList->Add(fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight);
        fOutputList->Add(fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight_prompt);
        
        
        
    }
	
	
   
	
	fvalueElectron = new Double_t[8];
    fvalueElectronTPC = new Double_t[6];
    fvalueMulti = new Double_t[6];
	
    //electron Sparse
	Int_t bins[8]={58,40,40,20,20, 40, 6, 40}; // pt, TPCnsig, E/p, M20, M02, E,phi, TPCnsig_old
	Double_t xmin[8]={2,-15,0,0,0,0,0,-15};
	Double_t xmax[8]={60,5,2,2,2,40,6, 5};
	fSparseElectron = new THnSparseD ("Electron","Electron",8,bins,xmin,xmax);
	fOutputList->Add(fSparseElectron);
    
    //electron Sparse TPC
    //electron Sparse
    Int_t binse[6]={40, 40, 40, 18, 6, 40}; // p, pt, TPCnsigma, eta, phi
    Double_t xmine[6]={1, 1,-15,-0.9,0, -15};
    Double_t xmaxe[6]={41,41,5,0.9,6, 5};
    fSparseElectronTPC = new THnSparseD ("Electrons TPC","Electrons TPC",6,binse,xmine,xmaxe);
    fOutputList->Add(fSparseElectronTPC);
    
    
    //multi Sparse
    Int_t binsm[6]    =       {20,50,20,50,20,100};
    Double_t xminm[6]    =    {-10,0,0,0,0,0};
    Double_t xmaxm[6]    =    { 10,500,200,500,200,1000};
    fSparseMulti         = new THnSparseD ("Multiplicity","Multiplicity;zvtx;V0M;SPDTracklets;Corrected_V0M;Corrected_SPDTracklets;ncharge;",6,binsm,xminm,xmaxm);
    fOutputList->Add(fSparseMulti);
			
//______________________________________________________________________
	
	PostData(1, fOutputList);
    //PostData(2,fListProfiles);
	
///______________________________________________________________________
}

//______________________________________________________________________
//Main loop
//Called for each event
void AliAnalysisTask_JPsi_EMCal::UserExec(Option_t *)
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
  
	
//PID response
	fPidResponse = fInputHandler->GetPIDResponse();
		

//Check PID response
	if(!fPidResponse)
	{
		AliDebug(1, "Using default PID Response");
		fPidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class()); 
	}
  

	Double_t *fListOfmotherkink = 0;
	Int_t fNumberOfVertices = 0; 
	Int_t fNumberOfMotherkink = 0;


//Vertex Selection
	
	fNevent->Fill(30);
    if(fIsAOD)
    {
        
           const AliAODVertex* trkVtx = fAOD->GetPrimaryVertex();
        
            Float_t zvtx = trkVtx->GetZ();
            fZvtx = zvtx;
            //printf("vertex =%f, vertex global =%f \n", zvtx, fZvtx);
            //all events with reconstructed vertex:
            fVtxZ[0]->Fill(fZvtx);
        
        
        
            //to cross check also SPD vertex//
            const AliAODVertex* spdVtx = fAOD->GetPrimaryVertexSPD();
            if(!spdVtx || spdVtx->GetNContributors()<=0)
            {
                fNevent->Fill(29);
            }
            if(spdVtx)
            {
                fNevent->Fill(28);
                if((!trkVtx || trkVtx->GetNContributors()<=0) && (spdVtx->GetNContributors()<=0))
                {
                    fNevent->Fill(27);
                }
            }
            //end of SPD cross check
        
        
        
        
           if(!trkVtx || trkVtx->GetNContributors()<=0)
           {//no vertex from tracks
              fNevent->Fill(26);
              return;
           }
           
           
           fNevent->Fill(25);
           //any vertex
           fVtxZ[1]->Fill(fZvtx);
           
           if(TMath::Abs(fZvtx) > fVertexCut) return;
        
           fVtxZ[2]->Fill(fZvtx);
           
           
    }
	else
	{
        const AliESDVertex *trkVtx = fESD->GetPrimaryVertex();
		Float_t zvtx = trkVtx->GetZ();
		fZvtx = zvtx;
		if(TMath::Abs(zvtx) > fVertexCut) return;
	}

	fNevent->Fill(24);
    
//Look for kink mother for AOD
	if(fIsAOD)
	{
        //fNumberOfVertices = 0;
        //fNumberOfMotherkink = 0;
        
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
	
    fNevent->Fill(23);
    
//----------V0M Multiplicity------------------
    AliAODVZERO *vzeroAOD = dynamic_cast<AliAODVZERO *>( dynamic_cast<AliAODEvent *>(fAOD)->GetVZEROData());
    Int_t V0AMult = static_cast<Int_t>(vzeroAOD->GetMTotV0A());
    Int_t V0CMult = static_cast<Int_t>(vzeroAOD->GetMTotV0C());
    Int_t V0Mult=V0AMult+V0CMult;
    
//------------SPDTracklets--------------------
    Int_t nTracklets = 0;
    Int_t nAcc = 0;
    Double_t etaRange = 1.0;
    
    AliAODTracklets *tracklets = static_cast<const AliAODEvent*>(fAOD)->GetTracklets();
    nTracklets = tracklets->GetNumberOfTracklets();
    for (Int_t nn = 0; nn < nTracklets; nn++) {
        Double_t theta = tracklets->GetTheta(nn);
        Double_t eta = -TMath::Log(TMath::Tan(theta/2.0));
        if (TMath::Abs(eta) < etaRange) nAcc++;
    }
    
    Double_t SPDMult = nAcc;
    
    fV0Mult=V0Mult;
    fSPDMult=SPDMult;
    
    //printf("Multiplicity from V0 = %d, multiplicity from SPD =%d, VertexZ = %f\n", fV0Mult, fSPDMult, fZvtx);
    
    
  
    //printf("Seed is %u \n", gRandom->GetSeed());
    //printf("Random is %f \n", gRandom->PoissonD(3));//fixed difference to check random generated
    
    //=======
    //correction for multiplicity
    TProfile2D* estimatorAvg = GetEstimatorHistogram(fAOD);
    //if(!isMC){
        if(estimatorAvg){
            //printf("Estimator SPD exists!\n");
            //correctednAcc=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,nAcc,Zvertex1,fRefMult));
            fSPDMult_corr = AliAnalysisTask_JPsi_EMCal::GetTrackletsMeanCorrection(estimatorAvg,nAcc,fZvtx,fRefMult, fAOD->GetRunNumber());
        
            //countCorr=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,countMult,Zvertex1,fRefMult));
        }
        else{
            fSPDMult_corr=fSPDMult;
            //if we do not load any correction, variable is the same as uncorrected
        }
    //}

 
    //printf("SPD =%f, SPD_corrected =%f\n", SPDMult, fSPDMult_corr);
    

    //V0M Correction
    Int_t vzeroMultACorr=V0AMult, vzeroMultCCorr=V0CMult, vzeroMultCorr=V0Mult;
    vzeroMultACorr = static_cast<Int_t>(AliESDUtils::GetCorrV0A(V0AMult,fZvtx));
    vzeroMultCCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0C(V0CMult,fZvtx));
    fV0Mult_corr = vzeroMultACorr + vzeroMultCCorr; // corrected V0M
    
    //V0 correction applied here
    TProfile2D* estimatorV0 = GetEstimatorHistogram_V0(fAOD);
    if(estimatorV0){
       // printf("Estimator V0 exists!\n");
        fV0Mult_corr2 = AliAnalysisTask_JPsi_EMCal::GetV0MeanCorrection(estimatorV0,fV0Mult_corr,fZvtx,fRefMult_V0, fAOD->GetRunNumber());
    }
    else{
        fV0Mult_corr2 =fV0Mult_corr;
    }

    
    Int_t TPCcls_event = fAOD->GetNumberOfTPCClusters();
    Int_t ITScls_layer1 =  fAOD->GetNumberOfITSClusters(0);//SPD
    Int_t ITScls_layer2 =  fAOD->GetNumberOfITSClusters(1);//SPD
    Int_t ITScls_layer3 =  fAOD->GetNumberOfITSClusters(2);
    Int_t ITScls_layer4 =  fAOD->GetNumberOfITSClusters(3);
    Int_t ITScls_layer5 =  fAOD->GetNumberOfITSClusters(4);
    Int_t ITScls_layer6 =  fAOD->GetNumberOfITSClusters(5);
    
    Int_t SSD_plus_SDD = ITScls_layer3+ITScls_layer4+ITScls_layer5+ITScls_layer6;
    
    //printf("Number of cluster on TPC: %d and SSD_plus_SDD: %d \n",TPCcls_event, SSD_plus_SDD);
    
    //printf("Number of cluster on ITS layers: %d, %d, %d, %d, %d, %d \n",ITScls_layer1, ITScls_layer2, ITScls_layer3, ITScls_layer4, ITScls_layer5, ITScls_layer6);
    
    fTPC_vs_ITScls[0]->Fill(TPCcls_event, SSD_plus_SDD);
    //printf("V0 =%d, V0_corrected =%f,  V0_corrected2 =%f\n", V0Mult, fV0Mult_corr, fV0Mult_corr2);
    
    if(fAOD->IsPileupFromSPDInMultBins()){
        //printf("This event is pileUp from AOD\n");
        fNevent->Fill(22);
        return;
    }
    
    fTPC_vs_ITScls[1]->Fill(TPCcls_event, SSD_plus_SDD);
    
    if(fAOD->IsPileupFromSPD(3.,0.8,3.,2.,5.)){
        //printf("This event is pileUp from AOD\n");
        fNevent->Fill(21);
        return;
    }
    
    fTPC_vs_ITScls[2]->Fill(TPCcls_event, SSD_plus_SDD);
    
    //new pileUp rejection
    Int_t minContributors=5;    //minimum contributors to the pilup vertices, multi-vertex
    Double_t minChi2=5.;
    Double_t minWeiZDiff=15;   //minimum of the sqrt of weighted distance between the primary and the pilup vertex, multi-vertex
    Bool_t checkPlpFromDifferentBC=kFALSE;
    
    AliAnalysisUtils utils;
    utils.SetMinPlpContribMV(minContributors); //Multi Vertex pileup selection
    utils.SetMaxPlpChi2MV(minChi2);   //max value of Chi2perNDF of the pileup vertex, multi-vertex
    utils.SetMinWDistMV(minWeiZDiff);
    utils.SetCheckPlpFromDifferentBCMV(checkPlpFromDifferentBC); //SPD Pileup slection
    Bool_t isPileupFromMV = utils.IsPileUpMV(fAOD);      //check for multi-vertexer pile-up
    
    if(isPileupFromMV){
        fNevent->Fill(20);
        return;
    }
    
    fTPC_vs_ITScls[3]->Fill(TPCcls_event, SSD_plus_SDD);
 
//______________________________________________________________________	
	
//Only events with at least 2 tracks are accepted
	Int_t fNOtrks =  fVevent->GetNumberOfTracks();
    if(fNOtrks<2){
        fNevent->Fill(19);
        return;
    }
	fNevent->Fill(18);
	
	
//______________________________________________________________________
//EMCal Trigger Selection (Threshold selection)
	
	TString firedTrigger;
	TString TriggerEG1("EG1"); //takes trigger with name with EG1, ex: CEMC7EG1-B-NOPF-CENTNOTRD  
	TString TriggerEG2("EG2");
	
    //DCAL
	TString TriggerDG1("DG1"); //takes trigger with name with EG1, ex: CEMC7EG1-B-NOPF-CENTNOTRD  
	TString TriggerDG2("DG2");

	
	if(fAOD) firedTrigger = fAOD->GetFiredTriggerClasses();
	else if(fESD) firedTrigger = fESD->GetFiredTriggerClasses();
	
    //Bool_t IsEventEMCALL0=kTRUE;
	Bool_t IsEventEMCALL1=kFALSE;
	
	if(firedTrigger.Contains(TriggerEG1)){ 
		fNevent->Fill(17);
		IsEventEMCALL1=kTRUE;
	}
	if(firedTrigger.Contains(TriggerEG2)){
		fNevent->Fill(16);
		IsEventEMCALL1=kTRUE;
	}
	
    //if the flag is for a given threshold and it was not fired, return.
    //EMCal trigger word
	if(fEMCEG1){
		if(!firedTrigger.Contains(TriggerEG1))return;
		if(firedTrigger.Contains(TriggerEG2)){
			fNevent->Fill(15);
            //EG2 has to be removed from EG1, because all EG2 events are used.
            return;
			
		}
		
	}
	
	if(fEMCEG2){
		if(!firedTrigger.Contains(TriggerEG2))return;
		if(firedTrigger.Contains(TriggerEG1)){
			fNevent->Fill(14);
		}
		
	}
	//=====================================================
	//DCal trigger word
	if(fEMCDG1){
		if(!firedTrigger.Contains(TriggerDG1))return;
		if(firedTrigger.Contains(TriggerDG2)){
				fNevent->Fill(13);
			
		}
		
	}
	
	if(fEMCDG2){
		if(!firedTrigger.Contains(TriggerDG2))return;
		if(firedTrigger.Contains(TriggerDG1)){
				fNevent->Fill(12);
		}
		
	}
    
    //=====================================================
    fNevent->Fill(11);
    //EMCal + DCal trigger words together
    if(fEMCEG1DG1){
        fNevent->Fill(10);
        if(!firedTrigger.Contains(TriggerDG1) && !firedTrigger.Contains(TriggerEG1)) return;
        
        //to remove double count from EG2 on EG1 (only for EG1 case... for EG2 we should take all events). We remove EG2 from EG1, since it is already used on EG2.
        if(firedTrigger.Contains(TriggerDG2) || firedTrigger.Contains(TriggerEG2)) return;
        
        fNevent->Fill(9);
        if(firedTrigger.Contains(TriggerDG1)){
            fNevent->Fill(8);//if passed, how much is DCal trigger
        }
        if(firedTrigger.Contains(TriggerEG1)){
            fNevent->Fill(7);//if passed, how much is EMCal trigger
        }
        
        
        
    }
    
    if(fEMCEG2DG2){
        fNevent->Fill(6);
        if(!firedTrigger.Contains(TriggerDG2) && !firedTrigger.Contains(TriggerEG2)) return;
        
        //(all EG2 events are used... )
        fNevent->Fill(5);
        if(firedTrigger.Contains(TriggerDG2)){
            fNevent->Fill(4);//if passed, how much is DCal trigger
        }
        if(firedTrigger.Contains(TriggerEG2)){
            fNevent->Fill(3);//if passed, how much is EMCal trigger
        }
    }
    fNevent->Fill(2);


	
	//==============================================
		//For event generator
	    Bool_t IsMB_gen = kFALSE;
		Bool_t IsPythiaCC_gen = kFALSE;
		Bool_t IsPythiaBB_gen = kFALSE;
		Bool_t IsPythiaB_gen = kFALSE;
		Bool_t IsJpsi2ee_gen = kFALSE;
		Bool_t IsB2JPsi2ee_gen = kFALSE;
	//==============================================

	
//======================================================================
	Int_t Nch = 0;
	if(fIsMC){

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
			
			//==========================================================================================================
			//Checking the generator of MC used for pp 13 TeV data (from Ivan)
			/*
			TList* GetCocktailList()
			gives the CocktailHeaders when reading ESDs/AODs (corresponding to fExteral=kFALSE/kTRUE)
			the AODMC header (and the aodmc array) is passed as an instance to MCEvent by the AliAODInputHandler
			*/
			
			//printf("\n Checking generator in MC: \n ");
			TString genname;
			TList *l = (TList*)fMCEvent->GetCocktailList();
			for (Int_t i = l->GetEntries()-1; i >= 0; i--){
				AliGenEventHeader* gh=(AliGenEventHeader*)l->At(i);
				genname=gh->GetName();
			//	printf("Generator Index %d - %s\n", i, genname.Data());
				
				//set flags for each kind of generator:
				if(genname.Contains("MB")){
					IsMB_gen =kTRUE;
					//printf("Contains MB\n");
						//fNevent->Fill(13);
				}
				if(genname.Contains("CC")){
					IsPythiaCC_gen =kTRUE;
                    //printf("Contains CC\n");
						//fNevent->Fill(14);
				}
				if(genname.Contains("BB")){
					IsPythiaBB_gen =kTRUE;
					//printf("Contains BB\n");
					//	//fNevent->Fill(15);
				}
				
				if((genname.Contains("B_1")) && (!genname.Contains("BB_1"))){
					IsPythiaB_gen  =kTRUE;
					//printf("Contains B_1\n");
						//fNevent->Fill(16);
				}
				if(genname.Contains("Jpsi2ee")){
					IsJpsi2ee_gen  =kTRUE;
					//printf("Contains Jpsi2ee\n");
						//fNevent->Fill(17);
				}
				if(genname.Contains("B2JPsi2ee")){
					IsB2JPsi2ee_gen=kTRUE;
					//printf("Contains B2JPsi2ee\n");
						//fNevent->Fill(18);
				}
				
			}
			
						
			/*
				Bool_t IsPythiaCC_gen = kFALSE;
				Bool_t IsPythiaBB_gen = kFALSE;
				Bool_t IsPythiaB_gen = kFALSE;
				Bool_t IsJpsi2ee_gen = kFALSE;
				Bool_t IsB2JPsi2ee_gen = kFALSE;
		    */
						
		
			
			
			if(IsMB_gen)fNevent2->Fill(0);
			if(IsPythiaCC_gen)fNevent2->Fill(1);
			if(IsPythiaBB_gen)fNevent2->Fill(2);
			if(IsPythiaB_gen)fNevent2->Fill(3);
			if(IsJpsi2ee_gen)fNevent2->Fill(4);
			if(IsB2JPsi2ee_gen)fNevent2->Fill(5);
			 
			
			
			
			//==========================================================================================================

			
			for(Int_t iMC = 0; iMC < fMCarray->GetEntries(); iMC++)
			{
				fMCparticle = (AliAODMCParticle*) fMCarray->At(iMC);
				if(fMCparticle->GetMother()>0) fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
				
				
					//=================================================================
					//checking the generator of each particle of event
					//if(fMCparticle->GetGeneratorIndex()==0)printf("This is a particle from MB event!\n");
					//if(fMCparticle->GetGeneratorIndex()==1)printf("This is a particle from OTHER event!\n");

				
				Int_t pdg = fMCparticle->GetPdgCode();
				
				//removed &&fMCparticle->Charge()!=0 requirement
				if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax)
				{
					//check pT distribution of each MC generator
                    if(IsMB_gen)fTracksMCPt[0]->Fill(fMCparticle->Pt());
                    if(IsPythiaCC_gen)fTracksMCPt[1]->Fill(fMCparticle->Pt());
                    if(IsPythiaBB_gen)fTracksMCPt[2]->Fill(fMCparticle->Pt());
                    if(IsPythiaB_gen)fTracksMCPt[3]->Fill(fMCparticle->Pt());
                    if(IsJpsi2ee_gen)fTracksMCPt[4]->Fill(fMCparticle->Pt());
                    if(IsB2JPsi2ee_gen)fTracksMCPt[5]->Fill(fMCparticle->Pt());
                    
                    
                    fPDG_values->Fill(fMCparticle->GetPdgCode());
                    
                    //if(fMCparticle->Charge()==0) printf("pdg code is %d\n",fMCparticle->GetPdgCode());
                    //Take all J/psi generated
                    //if(fMCparticle->IsPhysicalPrimary()){
                        
                        if(TMath::Abs(fMCparticle->GetPdgCode())==443)
                        {
                            fPtMCparticleAll_trueJPsi_pT->Fill(fMCparticle->Pt());
                            
                            Double_t weight = CalculateWeight(fMCparticle->Pt());//weight based on JPsi pT
                            fPtMCparticleAll_trueJPsi_pT_weight->Fill(fMCparticle->Pt(), weight);
                                                        
                            //to check for prompt J/psi
                            if(fMCparticle->GetMother()>0){
                                Int_t mpdg = fMCparticleMother->GetPdgCode();
                                if(mpdg > 500 && mpdg < 600)fPtMCparticleAll_trueJPsi_pT_weight_prompt->Fill(fMCparticle->Pt());//from B wo weight
                                else fPtMCparticleAll_trueJPsi_pT_weight_prompt->Fill(fMCparticle->Pt(), weight);//has mother but it is not B
                            }
                            if(fMCparticle->GetMother()<=0)fPtMCparticleAll_trueJPsi_pT_weight_prompt->Fill(fMCparticle->Pt(), weight);//no mother (prompt)
                            
                            //check J/psi pT distribution of each generator
                            if(IsMB_gen)fTracksMCPt[6]->Fill(fMCparticle->Pt());
                            if(IsPythiaCC_gen)fTracksMCPt[7]->Fill(fMCparticle->Pt());
                            if(IsPythiaBB_gen)fTracksMCPt[8]->Fill(fMCparticle->Pt());
                            if(IsPythiaB_gen)fTracksMCPt[9]->Fill(fMCparticle->Pt());
                            if(IsJpsi2ee_gen)fTracksMCPt[10]->Fill(fMCparticle->Pt());
                            if(IsB2JPsi2ee_gen)fTracksMCPt[11]->Fill(fMCparticle->Pt());
                        }
                   // }
     
					if( TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 2212 || TMath::Abs(pdg) == 321 || TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13 ) 
					{

						if(fMCparticle->IsPhysicalPrimary()) 
						{
							
							Bool_t MotherFound = FindMother(iMC);
							
                            //For JPsi analysis
							if(fMCparticle->GetMother()<0) return;
							
							fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
							if(fMCparticleMother->GetMother()>0)fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
							
							if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==443))
							{ 
								fPtMCparticleAll_e_from_JPsi->Fill(fMCparticle->Pt());
                                fPtMCparticleAll_JPsi_pT->Fill(fMCparticleMother->Pt());
							}
                            
                            //new histos for efficiency checks
                            if(fMCparticle->GetPdgCode()==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==443))
                            {
                                fPtMCparticleAll_e_from_JPsi_electron->Fill(fMCparticle->Pt());
                                fPtMCparticleAll_JPsi_pT_electron->Fill(fMCparticleMother->Pt());
                            }
                            if(fMCparticle->GetPdgCode()==-11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==443))
                            {
                                fPtMCparticleAll_e_from_JPsi_positron->Fill(fMCparticle->Pt());
                                fPtMCparticleAll_JPsi_pT_positron->Fill(fMCparticleMother->Pt());
                            }
                            
                            //denominator for tracking efficiency using all electrons
                            if(TMath::Abs(fMCparticle->GetPdgCode())==11)
                            {
                                fPtMCparticleAll_electrons->Fill(fMCparticle->Pt());
                            }
							//denominator for tracking efficiency using all particles
                            fPtMCparticleAll_particles->Fill(fMCparticle->Pt());
                            
								//
							
							if(MotherFound)
							{
								if(fIsHFE1){
										//denominator for total efficiency and tracking
										//unfolding: denominator is pt_MC and numerator is pt_reco
									fPtMCparticleAllHfe1->Fill(fMCparticle->Pt());
									
										//fEtaPhi_den->Fill(fMCparticle->Phi(),fMCparticle->Eta());

								} //denominator for total efficiency and tracking
							}
						}
					}
				}//eta cut
			}//loop tracks
            
            
            //===========
            //Generated MC charged particle sin eta < 1
            
 
            
            // loop over all tracks
            for (Int_t igen = 0; igen < fMCarray->GetEntriesFast(); igen++){
                AliAODMCParticle *mctrack=(AliAODMCParticle*)fMCarray->UncheckedAt(igen);
                Int_t charge = mctrack->Charge();
                Double_t eta = mctrack->Eta();
                Bool_t isPhysPrim = mctrack->IsPhysicalPrimary();
                if(charge!=0){
                    if(eta > -1.0 && eta < 1.0){
                        if(isPhysPrim){
                            Nch++;
                        }
                    }
                }
            }
            
            
		}//AOD
        
        
	}//close isMC
    
  //only fill after all event selection
  
    //reject event if SPD tracklet is less than 1
    if(fMultiAnalysis){
        if(fSPDMult_corr<=0){
            fNevent->Fill(1);
            return;
        }
    }
    
    
    fvalueMulti[0] = fZvtx;
    fvalueMulti[1] = fV0Mult_corr;
    fvalueMulti[2] = fSPDMult;
    fvalueMulti[3] = fV0Mult_corr2;
    fvalueMulti[4] = fSPDMult_corr;
    fvalueMulti[5] = Nch;
    
    
    if(fFill_MSparse)fSparseMulti->Fill(fvalueMulti);    // multiplicity from tracklets
    
    // Multiplicity histos WITH correctionS
    fVtxZ_V0->Fill(fZvtx, fV0Mult_corr2);
    fVtxZ_SPD->Fill(fZvtx, fSPDMult_corr);
    fV0_SPD->Fill(fV0Mult_corr2,fSPDMult_corr);
    
    
    
    fV0_nch->Fill(Nch, fV0Mult_corr2);
    fSPD_nch->Fill(Nch,fSPDMult_corr);
    
	
//______________________________________________________________________
	
    
    // Number of events
    fNevent->Fill(0);
    
    if(fMultiAnalysis){
        if(fSPDMult_corr>0 && fSPDMult_corr < 10)   fNevent_SPD_multi->Fill(1);
        if(fSPDMult_corr>=10 && fSPDMult_corr < 20) fNevent_SPD_multi->Fill(2);
        if(fSPDMult_corr>=20 && fSPDMult_corr < 30) fNevent_SPD_multi->Fill(3);
        if(fSPDMult_corr>=30 && fSPDMult_corr < 40) fNevent_SPD_multi->Fill(4);
        if(fSPDMult_corr>=40 && fSPDMult_corr < 100) fNevent_SPD_multi->Fill(5);
    
        if(fV0Mult_corr2>0 && fV0Mult_corr2<100)     fNevent_V0_multi->Fill(1);
        if(fV0Mult_corr2>=100 && fV0Mult_corr2<200)  fNevent_V0_multi->Fill(2);
        if(fV0Mult_corr2>=200 && fV0Mult_corr2<300) fNevent_V0_multi->Fill(3);
        if(fV0Mult_corr2>=300 && fV0Mult_corr2<400) fNevent_V0_multi->Fill(4);
        if(fV0Mult_corr2>=400 && fV0Mult_corr2<800) fNevent_V0_multi->Fill(5);
    }
    

	
	Int_t ClsNo = -999;
	if(!fIsAOD) ClsNo = fESD->GetNumberOfCaloClusters(); 
	else ClsNo = fAOD->GetNumberOfCaloClusters(); 
	
   
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
    Bool_t hasCls_aboveEG1=kFALSE;
    Bool_t hasCls_aboveEG2=kFALSE;
	
	AliVCluster *clust = 0x0;
	
	if(!fUseTender){
		if(fIsAOD){
            
             fNClusters_pure->Fill(ClsNo);
            
			for (Int_t i=0; i< ClsNo; i++ ){
				clust = (AliVCluster*) fAOD->GetCaloCluster(i);
				
				if(clust && clust->IsEMCAL())
				{
                    
                   
                    
					fECluster_pure->Fill(clust->E());
                    
                    if((clust->E())>=5.0) hasCls_aboveEG2=kTRUE;
                       
                    if((clust->E())>=10.0) hasCls_aboveEG1=kTRUE;
					
					/////////////// for Eta Phi distribution
					Float_t pos[3]={0,0,0};
					clust->GetPosition(pos);
					TVector3 vpos(pos[0],pos[1],pos[2]);
					Double_t cphi = vpos.Phi();
					Double_t ceta = vpos.Eta();
					
					

					
					///from emcal QA task
					if(cphi < 0) cphi = cphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.
															  // if(cphi > 1.39 && cphi < 3.265) ; //EMCAL : 80 < phi < 187
															  // if(cphi > 4.53 && cphi < 5.708) ; //DCAL  : 260 < phi < 327
                    
                    fEtaPhi_both->Fill(cphi,ceta);

					//emcal
					if(cphi > 1.39 && cphi < 3.265){
						fEtaPhi_emcal->Fill(cphi,ceta);
						fECluster_pure_emcal->Fill(clust->E());
                        
                         if(fMultiAnalysis){
                             if(fSPDMult_corr>0 && fSPDMult_corr < 10)   fECluster_pure_emcal_SPD1->Fill(clust->E());
                             if(fSPDMult_corr>=10 && fSPDMult_corr < 20) fECluster_pure_emcal_SPD2->Fill(clust->E());
                             if(fSPDMult_corr>=20 && fSPDMult_corr < 30) fECluster_pure_emcal_SPD3->Fill(clust->E());
                             if(fSPDMult_corr>=30 && fSPDMult_corr < 40) fECluster_pure_emcal_SPD4->Fill(clust->E());
                             if(fSPDMult_corr>=40 && fSPDMult_corr < 100) fECluster_pure_emcal_SPD5->Fill(clust->E());
                        
                             if(fV0Mult_corr2>0 && fV0Mult_corr2<100)      fECluster_pure_emcal_V01->Fill(clust->E());
                             if(fV0Mult_corr2>=100 && fV0Mult_corr2<200)   fECluster_pure_emcal_V02->Fill(clust->E());
                             if(fV0Mult_corr2>=200 && fV0Mult_corr2<300)  fECluster_pure_emcal_V03->Fill(clust->E());
                             if(fV0Mult_corr2>=300 && fV0Mult_corr2<400)  fECluster_pure_emcal_V04->Fill(clust->E());
                             if(fV0Mult_corr2>=400 && fV0Mult_corr2<800)  fECluster_pure_emcal_V05->Fill(clust->E());
                         }
 
					}
					
					//dcal
					if(cphi > 4.53 && cphi < 5.708){
						fEtaPhi_dcal->Fill(cphi,ceta);
						fECluster_pure_dcal->Fill(clust->E());
					}
	
				}
		
		 }
		}//is AOD
	}//is not Tender

	
    //for the track loop
    Int_t NTracks=0;
    if(!fUseTender) NTracks=fVevent->GetNumberOfTracks();
    
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    //To use tender
	if(fUseTender){
        
        //printf("It is using TENDER! \n\n");
        
        //global variables
        fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTenderTrackName));
        fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTenderClusterName)); //emcal correction
        NTracks = fTracks_tender->GetEntries();
        ClsNo = fCaloClusters_tender->GetEntries();
        
        
		
        //For cluster information from tender
        Int_t ClsNo_emcal=0.;
		for (Int_t i=0; i< ClsNo; i++ ){
			
			clust = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(i));
			if (!clust) {
			continue;
			}
			if(clust && clust->IsEMCAL())
			{
				fECluster_pure->Fill(clust->E());
                
                ClsNo_emcal++;
                
                if((clust->E())>=5.0) hasCls_aboveEG2=kTRUE;
                   
                if((clust->E())>=10.0) hasCls_aboveEG1=kTRUE;
                      
				
				/////////////// for Eta Phi distribution
				Float_t pos[3]={0,0,0};
				clust->GetPosition(pos);
				TVector3 vpos(pos[0],pos[1],pos[2]);
				Double_t cphi = vpos.Phi();
				Double_t ceta = vpos.Eta();
                
                ///from emcal QA task
                if(cphi < 0) cphi = cphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.
                // if(cphi > 1.39 && cphi < 3.265) ; //EMCAL : 80 < phi < 187
                // if(cphi > 4.53 && cphi < 5.708) ; //DCAL  : 260 < phi < 327
				
				fEtaPhi_both->Fill(cphi,ceta);
				
				//emcal
				if(cphi<3.9){
					fEtaPhi_emcal->Fill(cphi,ceta);
                    fECluster_pure_emcal->Fill(clust->E());
                    
				}
				
				//dcal
				if(cphi>=3.9){
					fEtaPhi_dcal->Fill(cphi,ceta);
					fECluster_pure_dcal->Fill(clust->E());
				}
                
                //emcal+dcal
                 if(fMultiAnalysis){
                     if(fSPDMult_corr>0 && fSPDMult_corr < 10)   fECluster_pure_emcal_SPD1->Fill(clust->E());
                     if(fSPDMult_corr>=10 && fSPDMult_corr < 20) fECluster_pure_emcal_SPD2->Fill(clust->E());
                     if(fSPDMult_corr>=20 && fSPDMult_corr < 30) fECluster_pure_emcal_SPD3->Fill(clust->E());
                     if(fSPDMult_corr>=30 && fSPDMult_corr < 40) fECluster_pure_emcal_SPD4->Fill(clust->E());
                     if(fSPDMult_corr>=40 && fSPDMult_corr < 100) fECluster_pure_emcal_SPD5->Fill(clust->E());
                
                     if(fV0Mult_corr2>0 && fV0Mult_corr2<100)      fECluster_pure_emcal_V01->Fill(clust->E());
                     if(fV0Mult_corr2>=100 && fV0Mult_corr2<200)   fECluster_pure_emcal_V02->Fill(clust->E());
                     if(fV0Mult_corr2>=200 && fV0Mult_corr2<300)  fECluster_pure_emcal_V03->Fill(clust->E());
                     if(fV0Mult_corr2>=300 && fV0Mult_corr2<400)  fECluster_pure_emcal_V04->Fill(clust->E());
                     if(fV0Mult_corr2>=400 && fV0Mult_corr2<800)  fECluster_pure_emcal_V05->Fill(clust->E());
                 }

			}
			
		}
        //to count only clusters on EMCal
        fNClusters_pure->Fill(ClsNo_emcal);
	}
   //all cases
    fV0->Fill(fV0Mult_corr2);
    fSPD->Fill(fSPDMult_corr);
    
    //cases with particle above trigger threshold (EG1)
    if(hasCls_aboveEG1){
        fV01->Fill(fV0Mult_corr2);
        fSPD1->Fill(fSPDMult_corr);
    }
    
    //cases with particle above trigger threshold (EG2)
    if(hasCls_aboveEG2){
        fV02->Fill(fV0Mult_corr2);
        fSPD2->Fill(fSPDMult_corr);
    }
    
    
    
    if(fSelect_trigger_events1 && hasCls_aboveEG1==kFALSE){
       // printf("Events rejected since hasCls_aboveEG1 is false \n");
        return;
    }
    if(fSelect_trigger_events2 && hasCls_aboveEG2==kFALSE){
        //printf("Events rejected since hasCls_aboveEG2 is false \n");
        return;
    }
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
    
	
//=======================================================================
//=======================================================================
//=======================================================================
///Track loop
	for(Int_t iTracks = 0; iTracks < NTracks; iTracks++) 
	{
        
		AliVParticle* Vtrack = 0x0;
		if(!fUseTender) Vtrack  = fVevent->GetTrack(iTracks);
		if(fUseTender) Vtrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(iTracks));

		if (!Vtrack) 
		{
			//printf("ERROR: Could not receive track %d\n", iTracks);
			continue;
		}
     
		AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
		AliESDtrack *etrack = dynamic_cast<AliESDtrack*>(Vtrack);
		AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);
		
		
		Double_t eta =0;
		eta = track->Eta();
		if(eta > fEtaCutMax || eta < fEtaCutMin) continue;
		
		
		Double_t fTPCnSigma_old = -999;
		Double_t fTPCnSigma = -999;
		Double_t fTPCnSigma_pion = -999;
		Double_t fTPCnSigma_proton = -999;
		Double_t fTPCnSigma_kaon = -999;
		Double_t fTPCsignal = -999;
     	Double_t fPt = -999;
		Double_t fP = -999;
		Double_t fP2 = -999;
		Double_t fPt2 = -999;
		
		
		//TOF
		Double_t fTOFnSigma = -999;
		Double_t fTOFsignal = -999;
		
		
		Float_t pos0[3]={0,0,0};
		Float_t pos1[3]={0,0,0};
		Float_t pos2[3]={0,0,0};

		///_____________________________________________________________________________
		///Fill QA plots without track selection
		fPt = track->Pt();
		fP = TMath::Sqrt((track->Pt())*(track->Pt()) + (track->Pz())*(track->Pz()));
		
		fTPCsignal = track->GetTPCsignal();
		fTPCnSigma_old = fPidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
		fTPCnSigma_pion = fPidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
		fTPCnSigma_proton = fPidResponse->NumberOfSigmasTPC(track, AliPID::kProton);
		fTPCnSigma_kaon = fPidResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        
        //Apply TPC calibration here!
        fTPCnsigma_p_beforeCalibration->Fill(fP,fTPCnSigma_old);
        
        if(fIs_TPC_calibration) fTPCnSigma = GetTPCCalibration(fAOD->GetRunNumber(), fTPCnSigma_old);
        if(!fIs_TPC_calibration) fTPCnSigma = fTPCnSigma_old;
        
        fTPCnsigma_p_afterCalibration->Fill(fP,fTPCnSigma);
        
        //TOF
        fTOFsignal = track->GetTOFsignal();
        fTOFnSigma = fPidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
		

        fTPC_p[0]->Fill(fP,fTPCsignal);
        fTPCnsigma_p[0]->Fill(fP,fTPCnSigma);
        
        fTOF_p[0]->Fill(fP,fTOFsignal);
        fTOFnsigma_p[0]->Fill(fP,fTOFnSigma);
        
        if(track->GetEMCALcluster()>0)
        {
				
			if(!fUseTender) fClus = fVevent->GetCaloCluster(track->GetEMCALcluster());
			if(fUseTender){
				int EMCalIndex = -1;
				EMCalIndex = track->GetEMCALcluster();
				if(EMCalIndex>0){
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
				
				fdEta_dPhi[0]->Fill(fClus->GetTrackDx(), fClus->GetTrackDz());                
				if(TMath::Abs(fClus->GetTrackDx())<=0.05 && TMath::Abs(fClus->GetTrackDz())<=0.05)
                {
                    fEoverP_pt[0]->Fill(fPt,(fClus->E() / fP));
                    
                    Float_t Energy	= fClus->E();
                    fECluster[0]->Fill(Energy);
					fTPCnsigma_EoverP[0]->Fill(fTPCnSigma, (fClus->E() / fP));
					
					
					//======================================// for Eta Phi distribution
					fClus->GetPosition(pos0);
					TVector3 vpos0(pos0[0],pos0[1],pos0[2]);
					Double_t cphi = vpos0.Phi();
					Double_t ceta = vpos0.Eta();
					
						///from emcal QA task
					if(cphi < 0) cphi = cphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.
						// if(cphi > 1.39 && cphi < 3.265) ; //EMCAL : 80 < phi < 187
						// if(cphi > 4.53 && cphi < 5.708) ; //DCAL  : 260 < phi < 327
					
					
					//emcal
					if(cphi > 1.39 && cphi < 3.265){
						fECluster_emcal[0]->Fill(fClus->E());
					}
					
                    //dcal
					if(cphi > 4.53 && cphi < 5.708){
						fECluster_dcal[0]->Fill(fClus->E());
					}
					//======================================
 
                }
            }
        }
        
       // fVtxZ[0]->Fill(fZvtx);
		fTracksPt[0]->Fill(fPt);
		
		
//=======================================================================
// Track Selection Cuts are applied here
//=======================================================================
        
        

		if(fAOD){
			
			//TPCncls
			//if(atrack->GetTPCNcls() < fTPCncls) continue;
            
            //TPC N crossedRows
            if(atrack->GetTPCCrossedRows() < fTPCnCrossedRows) continue;
            //if(RatioTPCclusters < fRatioCrossedRowOverFindable) return 0;
            //if(nclusN< fTPCNclusPID) return 0 ;
            
            
			fTracksQAPt[1]->Fill(fPt);
            if(fTPCandITSrefit){
                if((!(atrack->GetStatus()&AliESDtrack::kITSrefit))|| (!(atrack->GetStatus()&AliESDtrack::kTPCrefit))) continue;
            }
			fTracksQAPt[2]->Fill(fPt);
            //kAny
            if(fITSpixel==1){
                if(!(atrack->HasPointOnITSLayer(0) || atrack->HasPointOnITSLayer(1))) continue;
            }
            //kBoth
            if(fITSpixel==2){
                if(!(atrack->HasPointOnITSLayer(0) && atrack->HasPointOnITSLayer(1))) continue;
            }
            //kFirst
            if(fITSpixel==3){
                if(!(atrack->HasPointOnITSLayer(0))) continue;
            }
			fTracksQAPt[3]->Fill(fPt);
            
            //ITS Ncls
            if((atrack->GetITSNcls()) < fITSncls) continue;
            
            fTracksQAPt[4]->Fill(fPt);
			
			//DCA cut
			Double_t d0z0[2], cov[3];
			AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
            if(atrack->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0, cov)){
                if(TMath::Abs(d0z0[0]) > fDCAxyCut || TMath::Abs(d0z0[1]) > fDCAzCut) continue;
                
            }
			fTracksQAPt[5]->Fill(fPt);
			
            //Pt cut
			if(fPt<fPtCutMainEle) continue;
			fTracksQAPt[6]->Fill(fPt);
			
			//reject kink mother
            if(fRejectKinkMother){
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
            
			fTracksQAPt[7]->Fill(fPt);
			
            //chi2 per cluster
        
           if(((track->GetTPCchi2())/(atrack->GetTPCNcls())) > fTPCchi2){
                continue;
            }
            fTracksQAPt[8]->Fill(fPt);
            
            //ITS Chi2
            if(((atrack->GetITSchi2())/(atrack->GetITSNcls())) > fITSchi2){
                continue;
            }
				
            fTracksQAPt[9]->Fill(fPt);
            
            if(fAODGlobalTracks){
                if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts
            }
            fTracksQAPt[10]->Fill(fPt);
            
            
            

		}
		
        //if(atrack->GetTPCsignalN() < 80) TPCfor PID
        //if(atrack->GetTPCNclsF() < 0.6) findable
	
		        
//=======================================================================
// QA plots after track selection
//=======================================================================

		fTracksPt[1]->Fill(fPt);
        
		fTPC_p[1]->Fill(fP,fTPCsignal);
		fTPCnsigma_p[1]->Fill(fP,fTPCnSigma);
        
        fTOF_p[1]->Fill(fP,fTOFsignal);
        fTOFnsigma_p[1]->Fill(fP,fTOFnSigma);
        
        
        
       // printf("Track: %d, p: %f, pt: %f, TPCnsigma: %f, eta: %f, phi: %f \n", iTracks, fP,track->Pt(),fTPCnSigma, track->Eta(), track->Phi());
        
        if(track->Pt() >=2){
            
            fvalueElectronTPC[0] = fP;
            fvalueElectronTPC[1] = track->Pt();
            fvalueElectronTPC[2] = fTPCnSigma;
            fvalueElectronTPC[3] = track->Eta();
            fvalueElectronTPC[4] = track->Phi();
            fvalueElectronTPC[5] = fTPCnSigma_old;
            
            if(fFill_ESparseTPC)fSparseElectronTPC->Fill(fvalueElectronTPC);
        }
        
       // printf("SparseElectronTPC was filled \n");
		
			//MC studies
		
   if(fIsMC)
   {
	if(fIsAOD)
	{
		fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
		 
	    Int_t pdg = fMCparticle->GetPdgCode();
		 
	    if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax && fMCparticle->Charge()!=0)
	    {
		 
		    if( TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 2212 || TMath::Abs(pdg) == 321 || TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13 )
		    {
		 
		     if(fMCparticle->IsPhysicalPrimary()){
		 
               if(fMCparticle->GetMother()>0)
               {
		 
                   fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                   if(fMCparticleMother->GetMother()>0)fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
		 
                   if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==443)){
			    	 fPtMCparticleReco_e_from_JPsi->Fill(track->Pt()); //reconstructed pT
                   }
                   
                   //numerator tracking efficiency using all electrons
                   if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())!=22)){
                       fPtMCparticleReco_electrons_no_gamma->Fill(track->Pt()); //reconstructed pT
                   }
                   if(TMath::Abs(fMCparticle->GetPdgCode())==11){
                       fPtMCparticleReco_electrons->Fill(track->Pt()); //reconstructed pT
                   }
               }
               
               //numerator tracking efficiency using all particles
               fPtMCparticleReco_particles->Fill(track->Pt()); //reconstructed pT
               
		   }
		 }
	   }//eta cut
     }//close AOD
   }//close IsMC
		 
		
    //====== end of MC studies
        fIsTrack1Emcal=kFALSE;
        fIsTrack1Dcal=kFALSE;

		if(track->GetEMCALcluster()>0)
		{
				
            
            
			
			if(!fUseTender) fClus = fVevent->GetCaloCluster(track->GetEMCALcluster());
			
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            //to use tender
			if(fUseTender){
				int EMCalIndex = -1;
				EMCalIndex = track->GetEMCALcluster();
				if(EMCalIndex>0){
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
				fdEta_dPhi[1]->Fill(fClus->GetTrackDx(), fClus->GetTrackDz()); 
				//if(TMath::Abs(fClus->GetTrackDx())<=0.05 && TMath::Abs(fClus->GetTrackDz())<=0.05)
			  //{
			      fEoverP_pt[1]->Fill(fPt,(fClus->E() / fP));
                
                
				   
				  Float_t Energy	= fClus->E();
				  fECluster[1]->Fill(Energy);
                  fTracksQAPt[0]->Fill(fPt);
                
                
                //to check how many tracks matches the cluster
                fEoverP_ntracks_matched->Fill(fClus->E()/fP, fClus->GetNTracksMatched());
                fEoverP_ncells->Fill(fClus->E()/fP, fClus->GetNCells());
                
                
                
                //Ecluster for electrons on TPC
                if(fTPCnSigma >= fTPCnsigmaCutMin && fTPCnSigma <= fTPCnsigmaCutMax){
                     fECluster[2]->Fill(fClus->E());
                    //Ecluster for electrons on TPC and on EMCal
                    if((fClus->E() / fP) >= fEoverPCutMin && (fClus->E() / fP) <=fEoverPCutMax){
                         fECluster[3]->Fill(fClus->E());
                    }
                }
                
                //TOF signal for tracks matched to emcal
                fTOF_p[2]->Fill(fP,fTOFsignal);
                fTOFnsigma_p[2]->Fill(fP,fTOFnSigma);
                
				  fTPCnsigma_EoverP[1]->Fill(fTPCnSigma, (fClus->E() / fP));
				  
				  fNClusters[1]->Fill(ClsNo);
				  
                //======================================// for Eta Phi distribution
				  fClus->GetPosition(pos1);
				  TVector3 vpos1(pos1[0],pos1[1],pos1[2]);
				  Double_t cphi = vpos1.Phi();
				  Double_t ceta = vpos1.Eta();
		
				//======================================
				///from emcal QA task
				  if(cphi < 0) cphi = cphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.
															// if(cphi > 1.39 && cphi < 3.265) ; //EMCAL : 80 < phi < 187
															// if(cphi > 4.53 && cphi < 5.708) ; //DCAL  : 260 < phi < 327
				  
               
                //emcal
				  if(cphi > 1.39 && cphi < 3.265){
					  fECluster_emcal[1]->Fill(fClus->E());
                      //fIsTrack1Emcal=kTRUE;
				  }
				  
                //dcal
				  if(cphi > 4.53 && cphi < 5.708){
					  fECluster_dcal[1]->Fill(fClus->E());
                     // fIsTrack1Dcal=kTRUE;
				  }
                
                fIsTrack1Emcal=kTRUE;
			  
			    //}
                //EID THnsparse
				
				fvalueElectron[0] = track->Pt();
				fvalueElectron[1] = fTPCnSigma;
				fvalueElectron[2] = fClus->E() / fP;
				fvalueElectron[3] = fClus->GetM20();
				fvalueElectron[4] = fClus->GetM02();
				fvalueElectron[5] = fClus->E(); // to check rejection factor for electrons
                fvalueElectron[6] = cphi; //to separate emcal and dcal
                fvalueElectron[7] = fTPCnSigma_old; //to separate emcal and dcal
                //fvalueElectron[7] = fV0Mult;//to check RF in bins of multiplicity (bins not exactly same as in the analysis...)
               // fvalueElectron[8] = fSPDMult;//to check RF in bins of multiplicity (bins not exactly same as in the analysis...)
				
				if(fFill_ESparse)fSparseElectron->Fill(fvalueElectron);
				
			}
		}
		
		//fVtxZ[1]->Fill(fZvtx);
		
       // if(fSelect_trigger_events1 || fSelect_trigger_events2){
            //printf("Only Electron sparse is filled... rest of analysis is stopped here \n");
            //for J/psi analysis, set both as kFALSE
           // return;
            
       // }
			
//=======================================================================
//denominator for TPC PID efficiency
//=======================================================================
        
        if(fIsMC)
        {
            if(fIsAOD)
            {
                fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
                
                if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax && fMCparticle->Charge()!=0)
                {
                    
                    if(fMCparticle->GetMother()>0){
                        
                        fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                        if(fMCparticleMother->GetMother()>0){
                            fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
                        }
                        
                        if(fMCparticle->IsPhysicalPrimary())
                        {
                            if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==443)){
                                fPtMCparticle_TPCpid_e_from_JPsi->Fill(track->Pt()); //reconstructed pT
                            }
                            
                            //denominator TPCpid efficiency using all electrons
                            if(TMath::Abs(fMCparticle->GetPdgCode())==11  && (TMath::Abs(fMCparticleMother->GetPdgCode())!=22) ){
                                fPtMCparticle_TPCpid_electrons->Fill(track->Pt()); //reconstructed pT
                            }
                        }
                        
                    }// has mother
                }//eta cut
            }//close AOD
        }//close IsMC
        
//=======================================================================
// Here the PID cut defined in the file "Config.C" is applied
//=======================================================================

		
		if(fTPCnSigma < fTPCnsigmaCutMin || fTPCnSigma > fTPCnsigmaCutMax) continue;
        
        //printf("Main leg: Track1 on Electron band with fPt=%f\n", fPt);

	    fTracksQAPt[11]->Fill(fPt);
        
//=======================================================================
//numerator for TPC pid efficiency
//=======================================================================
       
        if(fIsMC)
        {
            if(fIsAOD)
            {
                fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
                
                if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax && fMCparticle->Charge()!=0)
                {
                    
                    if(fMCparticle->GetMother()>0){
                        
                           fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                           if(fMCparticleMother->GetMother()>0){
                               fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
                           }
                        
                           if(fMCparticle->IsPhysicalPrimary())
                           {
                                if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==443)){
                                    fPtMCparticle_TPCpid_e_from_JPsi_num->Fill(track->Pt()); //reconstructed pT
                                }
                            
                                //denominator TPCpid efficiency using all electrons
                                if(TMath::Abs(fMCparticle->GetPdgCode())==11  && (TMath::Abs(fMCparticleMother->GetPdgCode())!=22) ){
                                    fPtMCparticle_TPCpid_electrons_num->Fill(track->Pt()); //reconstructed pT
                                }
                            }
                        
                    }// has mother
                }//eta cut
            }//close AOD
        }//close IsMC
	
        
//Here I will check how many electrons matches the EMCal (track-matching efficiency for all electrons on TPC)
        if(fIsTrack1Emcal){
            //=======================================================================
            //numerator for EMCal track-matching efficiency --> denominator is num from TPC pid
            //=======================================================================
            
            if(fIsMC)
            {
                if(fIsAOD)
                {
                    fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
                    
                    if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax && fMCparticle->Charge()!=0)
                    {
                        
                        if(fMCparticle->GetMother()>0){
                            
                            fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                            if(fMCparticleMother->GetMother()>0){
                                fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
                            }
                            
                            if(fMCparticle->IsPhysicalPrimary())
                            {
                                if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==443)){
                                    fPtMCparticle_EMCal_TM_e_from_JPsi->Fill(track->Pt()); //reconstructed pT
                                }
                                
                                //denominator TPCpid efficiency using all electrons
                                if(TMath::Abs(fMCparticle->GetPdgCode())==11  && (TMath::Abs(fMCparticleMother->GetPdgCode())!=22) ){
                                    fPtMCparticle_EMCal_TM_electrons->Fill(track->Pt()); //reconstructed pT
                                }
                            }
                            
                        }// has mother
                    }//eta cut
                }//close AOD
            }//close IsMC
            
            
        }//close 'IsTrack1Emcal'
		
        
//=======================================================================
// QA plots after PID selection of first track
//=======================================================================
	  
		fTracksPt[2]->Fill(fPt);
		fTPC_p[2]->Fill(fP,fTPCsignal);
		fTPCnsigma_p[2]->Fill(fP,fTPCnSigma);
		
        ///selecting second track on TPC for the invariant mass
        Float_t charge1		= track->Charge();
        TLorentzVector v1(track->Px(),track->Py(),track->Pz(),track->P());
		
		for (Int_t lTracks = iTracks+1; lTracks < fNOtrks; lTracks++) {
				
				AliVParticle* Vtrack2 = 0x0;
				Vtrack2  = fVevent->GetTrack(lTracks);
				
				
				if (!Vtrack2) {
					//printf("ERROR: Could not receive track %d\n", lTracks);
					continue;
				}
			
				AliVTrack *track2 = dynamic_cast<AliVTrack*>(Vtrack2);
				AliESDtrack *etrack2 = dynamic_cast<AliESDtrack*>(Vtrack2);
				AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(Vtrack2);
			
			   fPt2 = track2->Pt();
			   fP2 = TMath::Sqrt((track2->Pt())*(track2->Pt()) + (track2->Pz())*(track2->Pz()));
			
			
				Double_t eta2=0;
			    eta2 = track2->Eta();
			    if(eta2 > fEtaCutMax || eta2 < fEtaCutMin) continue;

			   fTracksPt[3]->Fill(fPt2);
			
							
				//=======================================================================
				// Track Selection Cuts are applied here
				//=======================================================================
			
        if(fAOD){
                
                //TPCncls
                //if(atrack2->GetTPCNcls() < fTPCncls) continue;
                if(atrack2->GetTPCCrossedRows() < fTPCnCrossedRows) continue;
            
                if(fTPCandITSrefit){
                    if((!(atrack2->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack2->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
                }
                
                //kAny
                if(fITSpixel==1){
                    if(!(atrack2->HasPointOnITSLayer(0) || atrack2->HasPointOnITSLayer(1))) continue;
                }
                //kBoth
                if(fITSpixel==2){
                    if(!(atrack2->HasPointOnITSLayer(0) && atrack2->HasPointOnITSLayer(1))) continue;
                }
                //kFirst
                if(fITSpixel==3){
                    if(!(atrack2->HasPointOnITSLayer(0))) continue;
                }
              

                //ITS Ncls
                if((atrack2->GetITSNcls()) < fITSncls) continue;
                
                
                //DCA cut
                Double_t d0z03[2], cov2[3];
                AliAODVertex *pVtx3 = fAOD->GetPrimaryVertex();
                if(atrack2->PropagateToDCA(pVtx3, fAOD->GetMagneticField(), 20., d0z03, cov2)){
                    if(TMath::Abs(d0z03[0]) > fDCAxyCut || TMath::Abs(d0z03[1]) > fDCAzCut) continue;
                }
               
                //Pt cut
                if(fPt2<fPtCutPartner) continue;
                
                //reject kink mother
                if(fRejectKinkMother){
                    Bool_t kinkmotherpass3 = kTRUE;
                    for(Int_t kinkmother3 = 0; kinkmother3 < fNumberOfMotherkink; kinkmother3++)
                    {
                        if(track2->GetID() == fListOfmotherkink[kinkmother3])
                        {
                            kinkmotherpass3 = kFALSE;
                            continue;
                        }
                    }
                    if(!kinkmotherpass3) continue;
                }

                //chi2 per cluster
                 if(((track2->GetTPCchi2())/(atrack2->GetTPCNcls())) > fTPCchi2) continue;
                
                //ITS Chi2
                if(((atrack2->GetITSchi2())/(atrack2->GetITSNcls())) > fITSchi2){
                    continue;
                }
               

                if(fAODGlobalTracks){
                    if(!atrack2->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts
                }
              
                
            }
			
			fTracksPt[4]->Fill(fPt2);

			 
			 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				Double_t dEdx3 =-999, fTPCnSigma2=-999;
				dEdx3 = track2->GetTPCsignal();
				fTPCnSigma2 = fPidResponse->NumberOfSigmasTPC(track2, AliPID::kElectron);
				if(fTPCnSigma2 >= fTPCnsigmaCutMin && fTPCnSigma2 <= fTPCnsigmaCutMax){
                    
                    //printf("Second leg: Track2 on Electron band with fPt2=%f\n", fPt2);
					
					fTracksPt[5]->Fill(fPt2);
					
					Float_t charge2	= track2->Charge();
					
					TLorentzVector v2(track2->Px(),track2->Py(),track2->Pz(),track2->P());
					Float_t invmass3 = (v1+v2).M(); // Inv mass of pair
					Float_t pt3 = (v1+v2).Pt();  // pt of pair
					
					//=====================================================================================================
					//KFParticle to take the particle pT on the right position to calculate the mass
					//only for AOD
					
					
					Double_t bfield = fAOD->GetMagneticField();
					
					Int_t fPDGtrack1 = 11; 
					Int_t fPDGtrack2 = 11;
					
					if(charge1>0) fPDGtrack1 = -11;
					if(charge2>0) fPDGtrack2 = -11;
					
					AliKFParticle::SetField(bfield);
					AliKFParticle fKFtrack1(*track, fPDGtrack1);
					AliKFParticle fKFtrack2(*track2, fPDGtrack2);
					AliKFParticle fRecoPair(fKFtrack1, fKFtrack2);
					
					//Reconstruction Cuts
					if(fRecoPair.GetNDF()<1) continue;
					// Double_t chi2OverNDF = fRecoPair.GetChi2()/fRecoGamma.GetNDF();
					//if(TMath::Sqrt(TMath::Abs(chi2OverNDF))>fChi2OverNDFCut) continue;
					
					//Invariant Mass
					Double_t imass; 
					Double_t err_imass;//sigma
					fRecoPair.GetMass(imass,err_imass);
					
					//Pt
					Double_t pt_kf; 
					Double_t err_pt_kf;//sigma
					fRecoPair.GetPt(pt_kf,err_pt_kf);
                    
                    //weight to correct pT spectra by the efficiency of J/psi reco and Id
                    
                    /*
                    //if use pT of electrons
                    Double_t a = 0.0000754390;
                    Double_t b = 0.0357673;
                    Double_t c = 6.50168;
                    Double_t d = 1.31532;
                    Double_t weight = 1./(a+(b/(1.0+exp(-(pt_kf-c)/d))));  //   weight = 1/eff
                    */
                    
                    //if use pT of J/psi (correction for EG1 pT > 10 GeV/c)
                   // Double_t a = 11.11;
                    //Double_t b = 10.4;
                    //Double_t c = 0.0333;
                    //if use pT of J/psi (correction for EG1 pT > 5 GeV/c)
                    Double_t a = 5.20;
                    Double_t b = 8.6;
                    Double_t c = 0.0342;
                    Double_t weight = 1./(TMath::Erf((pt_kf-a)/b)*c);  //   weight = 1/eff
                    Double_t weight2 = 1./(TMath::Erf((25.0-a)/b)*c);  //   weight2 to normalize correction and keep right counts and right errors.
                    
                   // printf("weight is %f \n", weight);
                    
					
                    //printf("Mass from TLorentzVector: %f, Mass from KFParticle: %f\n", invmass3, imass);
                    //printf("Momentum from TLorentzVector: %f, Momentum from KFParticle: %f\n", pt3, pt_kf);
					//=====================================================================================================
					
					if(charge1*charge2 <0) fHist_InvMass_pt_ULStpc->Fill(pt_kf,imass);
					if(charge1*charge2 >0) fHist_InvMass_pt_LStpc->Fill(pt_kf,imass);
					  
					
					
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//Using the tracks selected by TPC, now we test if the second one is on EMCal
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					
					
						//second leg
                    fIsTrack2Emcal=kFALSE;
                    fIsTrack2Dcal=kFALSE;
			     	if(track2->GetEMCALcluster()>0){

						if(!fUseTender) fClus2 = fVevent->GetCaloCluster(track2->GetEMCALcluster());
                        
                        if(fUseTender){
                            int EMCalIndex3 = -1;
                            EMCalIndex3 = track2->GetEMCALcluster();
                            if(EMCalIndex3>0){
                                fClus2 = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(EMCalIndex3));
                                if (!fClus2) {
                                    continue;
                                }
                            }
                        }
                       
						if(fClus2->IsEMCAL())
						{
							fTracksPt[7]->Fill(fPt2);
							fClus2->GetPosition(pos2);
							TVector3 vpos2(pos2[0],pos2[1],pos2[2]);
							Double_t cphi = vpos2.Phi();
							Double_t ceta = vpos2.Eta();
							
							if(cphi < 0) cphi = cphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.
																	  // if(cphi > 1.39 && cphi < 3.265) ; //EMCAL : 80 < phi < 187
																	  // if(cphi > 4.53 && cphi < 5.708) ; //DCAL  : 260 < phi < 327
							
                            //emcal
							if(cphi > 1.39 && cphi < 3.265){
								fECluster_emcal[2]->Fill(fClus2->E());
  							}
                            //dcal
							if(cphi > 4.53 && cphi < 5.708){
								fECluster_dcal[2]->Fill(fClus2->E());
                            }
                            fIsTrack2Emcal=kTRUE;
						
						}
					}
					
					//Filling the invariant mass spectrum
					
					if(fIsTrack1Emcal && (!fIsTrack2Emcal)){
                        
						
                        //tpc electrons but at least one leg with track-matching
                        if(charge1*charge2 <0) fHist_InvMass_pt_ULStpc_wMatching->Fill(pt_kf,imass);
                        if(charge1*charge2 >0) fHist_InvMass_pt_LStpc_wMatching->Fill(pt_kf,imass);
                        
						  if((fClus->E() / fP) >=fEoverPCutMin && (fClus->E() / fP) <=fEoverPCutMax && (fClus->E()) >= fEnergyCut){
							
                            //KFParticle
                              if(charge1*charge2 <0){
                                
                                  fHist_InvMass_pt_ULS_KF->Fill(pt_kf,imass);//multi integrated
                               
                                 if(fMultiAnalysis) fHist_InvMass_pt_ULS_KF_weight->Fill(pt_kf,imass, weight/weight2);//multi integrated with weight
                                  
                                 //correlation between leg1 and leg2
                                  fHist_Correlation_leg1_emcal_leg2_not->Fill(fPt, fPt2);
                                  
                              }
							  if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF->Fill(pt_kf,imass);
                              
                            //multiplicity bins histos (only ULS for SPDmulti and V0multi)
                              if(charge1*charge2 <0){
                                  
                                  if(fMultiAnalysis){
                                      if(fSPDMult_corr>0 && fSPDMult_corr < 10)   fHist_InvMass_pt_ULS_KF_SPDmulti_1->Fill(pt_kf,imass);
                                      if(fSPDMult_corr>=10 && fSPDMult_corr < 20) fHist_InvMass_pt_ULS_KF_SPDmulti_2->Fill(pt_kf,imass);
                                      if(fSPDMult_corr>=20 && fSPDMult_corr < 30) fHist_InvMass_pt_ULS_KF_SPDmulti_3->Fill(pt_kf,imass);
                                      if(fSPDMult_corr>=30 && fSPDMult_corr < 40) fHist_InvMass_pt_ULS_KF_SPDmulti_4->Fill(pt_kf,imass);
                                      if(fSPDMult_corr>=40 && fSPDMult_corr < 100) fHist_InvMass_pt_ULS_KF_SPDmulti_5->Fill(pt_kf,imass);
                                  
                                      if(fV0Mult_corr2>0 && fV0Mult_corr2<100)     fHist_InvMass_pt_ULS_KF_V0multi_1->Fill(pt_kf,imass);
                                      if(fV0Mult_corr2>=100 && fV0Mult_corr2<200)  fHist_InvMass_pt_ULS_KF_V0multi_2->Fill(pt_kf,imass);
                                      if(fV0Mult_corr2>=200 && fV0Mult_corr2<300) fHist_InvMass_pt_ULS_KF_V0multi_3->Fill(pt_kf,imass);
                                      if(fV0Mult_corr2>=300 && fV0Mult_corr2<400) fHist_InvMass_pt_ULS_KF_V0multi_4->Fill(pt_kf,imass);
                                      if(fV0Mult_corr2>=400 && fV0Mult_corr2<800) fHist_InvMass_pt_ULS_KF_V0multi_5->Fill(pt_kf,imass);
                                  
                                      //with weight
                                      if(fSPDMult_corr>0 && fSPDMult_corr < 10)   fHist_InvMass_pt_ULS_KF_SPDmulti_1_weight->Fill(pt_kf,imass,weight/weight2);
                                      if(fSPDMult_corr>=10 && fSPDMult_corr < 20) fHist_InvMass_pt_ULS_KF_SPDmulti_2_weight->Fill(pt_kf,imass,weight/weight2);
                                      if(fSPDMult_corr>=20 && fSPDMult_corr < 30) fHist_InvMass_pt_ULS_KF_SPDmulti_3_weight->Fill(pt_kf,imass,weight/weight2);
                                      if(fSPDMult_corr>=30 && fSPDMult_corr < 40) fHist_InvMass_pt_ULS_KF_SPDmulti_4_weight->Fill(pt_kf,imass,weight/weight2);
                                      if(fSPDMult_corr>=40 && fSPDMult_corr < 100) fHist_InvMass_pt_ULS_KF_SPDmulti_5_weight->Fill(pt_kf,imass,weight/weight2);
                                  
                                      if(fV0Mult_corr2>0 && fV0Mult_corr2<100)     fHist_InvMass_pt_ULS_KF_V0multi_1_weight->Fill(pt_kf,imass,weight/weight2);
                                      if(fV0Mult_corr2>=100 && fV0Mult_corr2<200)  fHist_InvMass_pt_ULS_KF_V0multi_2_weight->Fill(pt_kf,imass,weight/weight2);
                                      if(fV0Mult_corr2>=200 && fV0Mult_corr2<300) fHist_InvMass_pt_ULS_KF_V0multi_3_weight->Fill(pt_kf,imass,weight/weight2);
                                      if(fV0Mult_corr2>=300 && fV0Mult_corr2<400) fHist_InvMass_pt_ULS_KF_V0multi_4_weight->Fill(pt_kf,imass,weight/weight2);
                                      if(fV0Mult_corr2>=400 && fV0Mult_corr2<800) fHist_InvMass_pt_ULS_KF_V0multi_5_weight->Fill(pt_kf,imass,weight/weight2);
                                  
                                  }
   
                                }
							
							//leg 1 on emcal
							if(charge1*charge2 <0) fHist_InvMass_pt_ULS1->Fill(pt_kf,imass);
							if(charge1*charge2 >0) fHist_InvMass_pt_LS1->Fill(pt_kf,imass);
							
							fTracksPt[8]->Fill(fPt2);
							
							
							//===============================
							
							if(fIsMC)
							{

								//MC generators
								if(IsPythiaCC_gen){
									if(charge1*charge2 <0) fHist_InvMass_pt_ULS_KF_CC->Fill(pt_kf,imass);
									if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF_CC->Fill(pt_kf,imass);
								}
								if(IsPythiaBB_gen){
									if(charge1*charge2 <0) fHist_InvMass_pt_ULS_KF_BB->Fill(pt_kf,imass);
									if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF_BB->Fill(pt_kf,imass);
								}
								if(IsPythiaB_gen){
									if(charge1*charge2 <0) fHist_InvMass_pt_ULS_KF_B->Fill(pt_kf,imass);
									if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF_B->Fill(pt_kf,imass);
								}
								if(IsJpsi2ee_gen){
									if(charge1*charge2 <0) fHist_InvMass_pt_ULS_KF_Jpsi->Fill(pt_kf,imass);
									if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF_Jpsi->Fill(pt_kf,imass);
								}
								if(IsB2JPsi2ee_gen){
									if(charge1*charge2 <0) fHist_InvMass_pt_ULS_KF_BJpsi->Fill(pt_kf,imass);
									if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF_BJpsi->Fill(pt_kf,imass);
								}
								
								
								if(fIsAOD)
								{
									fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
									fMCparticle2 = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track2->GetLabel()));
									
																		
									Int_t pdg = fMCparticle->GetPdgCode();
									Int_t pdg3 = fMCparticle2->GetPdgCode();
									
									if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax && fMCparticle->Charge()!=0)
									{
									  if(fMCparticle2->Eta()>=fEtaCutMin && fMCparticle2->Eta()<=fEtaCutMax && fMCparticle2->Charge()!=0){
										
										if( TMath::Abs(pdg) == 11 && TMath::Abs(pdg3) == 11 ) 
										{	
											
											if(fMCparticle->IsPhysicalPrimary() && fMCparticle2->IsPhysicalPrimary())
											{
 
												if(fMCparticle->GetMother()<0 || fMCparticle2->GetMother()<0) return;
												
												fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
												fMCparticleMother2 = (AliAODMCParticle*) fMCarray->At(fMCparticle2->GetMother());
												
													
												if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==443))
												{ 
													if(TMath::Abs(fMCparticle2->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother2->GetPdgCode())==443))
													{
														
														
														fPtMCparticle_Total_e_from_JPsi->Fill(track->Pt());
                                                        //leg1 passed emcal cuts. Leg2 is an electron from JPsi and is not on EMCal
                                                        fPtMCparticle_EMCalpid_leg1_e_from_JPsi->Fill(track->Pt());
                                                        
                                                        //checking if they are from same mother
                                                        if((fMCparticleMother->GetLabel())==(fMCparticleMother2->GetLabel())){
                                                            //printf("electrons from same mother\n\n");
                                                            fPtMCparticle_Total_e_from_JPsi_sameMother->Fill(track->Pt()); //reconstructed pT
                                                            fPtMCparticle_Total_JPsi_pT->Fill(pt_kf);
                                               
                                                            if(invmass3>=fMassCutMin && invmass3<=fMassCutMax){
                                                                fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother->Fill(track->Pt()); //reconstructed pT
                                                                fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother->Fill(pt_kf);//spectrum of reconstructed J/Psi
                                                                
                                                                //weights calculated based on J/Psi true MC pT, but applied on e+e- pair pt
                                                                Double_t weight2 = CalculateWeight(fMCparticleMother->Pt());
                                                                fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight->Fill(pt_kf, weight2);
                                                                
                                                                //to check for prompt J/psi -> check JPsi mother
                                                                if(fMCparticleMother->GetMother()>0){
                                                                    fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
                                                                    Int_t gmpdg = fMCparticleGMother->GetPdgCode();
                                                                    //printf("mother of J/psi is %d\n", gmpdg);
                                                                    if(gmpdg > 500 && gmpdg < 600)fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight_prompt->Fill(pt_kf);//from B without weight
                                                                    else fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight_prompt->Fill(pt_kf, weight2);//prompt with weight (has mother but it is not B meson)
                                                                }
                                                                //if there is no mother, has to fill (prompt J/psi)
                                                                if(fMCparticleMother->GetMother()<=0){
                                                                    //printf("case of J/psi without mother\n");
                                                                    fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight_prompt->Fill(pt_kf, weight2);//prompt with weight
                                                                }
                                                                
               
														}//mass cut
                                                      }//same mother
													}//track2 true e from true J/psi
												}//track1 true e from true J/psi
											}//Is Physical primary cut
										  }//etacut second leg
										}
									}//eta cut
								}//close AOD
							}//close IsMC

							
							//===============================

						}
                      //}//close trigger for emcal
						
					}
					
					if((!fIsTrack1Emcal) && fIsTrack2Emcal){
                       // printf("Track2 is on EMCal and track1 is not \n");
                        
                       // if(fEMCEG1 || fEMCEG2){
                            
                         // printf("The cuts are: %f , E/p < %f and E > %f\n",fEoverPCutMin, fEoverPCutMax, fEnergyCut);
                        
                        //tpc electrons but at least one leg with track-matching
                        if(charge1*charge2 <0) fHist_InvMass_pt_ULStpc_wMatching->Fill(pt_kf,imass);
                        if(charge1*charge2 >0) fHist_InvMass_pt_LStpc_wMatching->Fill(pt_kf,imass);
                        
                        
						
						 if((fClus2->E()/fP2) >=fEoverPCutMin && (fClus2->E()/fP2) <=fEoverPCutMax && (fClus2->E()) >= fEnergyCut){
							
                           // printf("Track2 PASSED the cuts \n");
                           // printf("Track2 has pt=%f \n", fPt2);
							
							//if(charge1*charge2 <0) fHist_InvMass_pt_ULS->Fill(pt3,invmass3);
							//if(charge1*charge2 >0) fHist_InvMass_pt_LS->Fill(pt3,invmass3);
							
								//KFParticle
                             if(charge1*charge2 <0){
                                 fHist_InvMass_pt_ULS_KF->Fill(pt_kf,imass);//multi integrated
                                  if(fMultiAnalysis)fHist_InvMass_pt_ULS_KF_weight->Fill(pt_kf,imass,weight/weight2);//multi integrated with weight
                                 
                                 //correlation between leg1 and leg2
                                 fHist_Correlation_leg1_not_leg2_emcal->Fill(fPt, fPt2);
                                 
                                 
                             }
							if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF->Fill(pt_kf,imass);
                             
                             //multiplicity bins histos (only ULS for SPDmulti and V0multi)
                             if(charge1*charge2 <0){
                                 //printf("track2: fSPDMult = %f, fVOMult = %f", fSPDMult,fV0Mult );
                                 
                                 if(fMultiAnalysis){
                                     if(fSPDMult_corr>0 && fSPDMult_corr < 10)   fHist_InvMass_pt_ULS_KF_SPDmulti_1->Fill(pt_kf,imass);
                                     if(fSPDMult_corr>=10 && fSPDMult_corr < 20) fHist_InvMass_pt_ULS_KF_SPDmulti_2->Fill(pt_kf,imass);
                                     if(fSPDMult_corr>=20 && fSPDMult_corr < 30) fHist_InvMass_pt_ULS_KF_SPDmulti_3->Fill(pt_kf,imass);
                                     if(fSPDMult_corr>=30 && fSPDMult_corr < 40) fHist_InvMass_pt_ULS_KF_SPDmulti_4->Fill(pt_kf,imass);
                                     if(fSPDMult_corr>=40 && fSPDMult_corr < 100) fHist_InvMass_pt_ULS_KF_SPDmulti_5->Fill(pt_kf,imass);
                                 
                                     if(fV0Mult_corr2>0 && fV0Mult_corr2<100)     fHist_InvMass_pt_ULS_KF_V0multi_1->Fill(pt_kf,imass);
                                     if(fV0Mult_corr2>=100 && fV0Mult_corr2<200)  fHist_InvMass_pt_ULS_KF_V0multi_2->Fill(pt_kf,imass);
                                     if(fV0Mult_corr2>=200 && fV0Mult_corr2<300) fHist_InvMass_pt_ULS_KF_V0multi_3->Fill(pt_kf,imass);
                                     if(fV0Mult_corr2>=300 && fV0Mult_corr2<400) fHist_InvMass_pt_ULS_KF_V0multi_4->Fill(pt_kf,imass);
                                     if(fV0Mult_corr2>=400 && fV0Mult_corr2<800) fHist_InvMass_pt_ULS_KF_V0multi_5->Fill(pt_kf,imass);
                                 
                                     //with weight
                                     if(fSPDMult_corr>0 && fSPDMult_corr < 10)   fHist_InvMass_pt_ULS_KF_SPDmulti_1_weight->Fill(pt_kf,imass,weight/weight2);
                                     if(fSPDMult_corr>=10 && fSPDMult_corr < 20) fHist_InvMass_pt_ULS_KF_SPDmulti_2_weight->Fill(pt_kf,imass,weight/weight2);
                                     if(fSPDMult_corr>=20 && fSPDMult_corr < 30) fHist_InvMass_pt_ULS_KF_SPDmulti_3_weight->Fill(pt_kf,imass,weight/weight2);
                                     if(fSPDMult_corr>=30 && fSPDMult_corr < 40) fHist_InvMass_pt_ULS_KF_SPDmulti_4_weight->Fill(pt_kf,imass,weight/weight2);
                                     if(fSPDMult_corr>=40 && fSPDMult_corr < 100) fHist_InvMass_pt_ULS_KF_SPDmulti_5_weight->Fill(pt_kf,imass,weight/weight2);
                                 
                                     if(fV0Mult_corr2>0 && fV0Mult_corr2<100)     fHist_InvMass_pt_ULS_KF_V0multi_1_weight->Fill(pt_kf,imass,weight/weight2);
                                     if(fV0Mult_corr2>=100 && fV0Mult_corr2<200)  fHist_InvMass_pt_ULS_KF_V0multi_2_weight->Fill(pt_kf,imass,weight/weight2);
                                     if(fV0Mult_corr2>=200 && fV0Mult_corr2<300) fHist_InvMass_pt_ULS_KF_V0multi_3_weight->Fill(pt_kf,imass,weight/weight2);
                                     if(fV0Mult_corr2>=300 && fV0Mult_corr2<400) fHist_InvMass_pt_ULS_KF_V0multi_4_weight->Fill(pt_kf,imass,weight/weight2);
                                     if(fV0Mult_corr2>=400 && fV0Mult_corr2<800) fHist_InvMass_pt_ULS_KF_V0multi_5_weight->Fill(pt_kf,imass,weight/weight2);
                                 
                                 }
                                 
                             }
                             
                             
							
								//leg 2 on emcal
							if(charge1*charge2 <0) fHist_InvMass_pt_ULS2->Fill(pt_kf,imass);
							if(charge1*charge2 >0) fHist_InvMass_pt_LS2->Fill(pt_kf,imass);
							
							fTracksPt[9]->Fill(fPt2);
							
							
							//===============================
							
							if(fIsMC)
							{
								
								//MC generators
								if(IsPythiaCC_gen){
									if(charge1*charge2 <0) fHist_InvMass_pt_ULS_KF_CC->Fill(pt_kf,imass);
									if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF_CC->Fill(pt_kf,imass);
								}
								if(IsPythiaBB_gen){
									if(charge1*charge2 <0) fHist_InvMass_pt_ULS_KF_BB->Fill(pt_kf,imass);
									if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF_BB->Fill(pt_kf,imass);
								}
								if(IsPythiaB_gen){
									if(charge1*charge2 <0) fHist_InvMass_pt_ULS_KF_B->Fill(pt_kf,imass);
									if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF_B->Fill(pt_kf,imass);
								}
								if(IsJpsi2ee_gen){
									if(charge1*charge2 <0) fHist_InvMass_pt_ULS_KF_Jpsi->Fill(pt_kf,imass);
									if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF_Jpsi->Fill(pt_kf,imass);
								}
								if(IsB2JPsi2ee_gen){
									if(charge1*charge2 <0) fHist_InvMass_pt_ULS_KF_BJpsi->Fill(pt_kf,imass);
									if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF_BJpsi->Fill(pt_kf,imass);
								}

								
								
								
								if(fIsAOD)
								{
									fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
									fMCparticle2 = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track2->GetLabel()));
									
									Int_t pdg = fMCparticle->GetPdgCode();
									Int_t pdg3 = fMCparticle2->GetPdgCode();
									
									if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax && fMCparticle->Charge()!=0)
									{
										if(fMCparticle2->Eta()>=fEtaCutMin && fMCparticle2->Eta()<=fEtaCutMax && fMCparticle2->Charge()!=0){
											
											if( TMath::Abs(pdg) == 11 && TMath::Abs(pdg3) == 11 ) 
											{	
												
												if(fMCparticle->IsPhysicalPrimary() && fMCparticle2->IsPhysicalPrimary())
												{
													
														//For JPsi analysis
													if(fMCparticle->GetMother()<0 || fMCparticle2->GetMother()<0) return;
													
													fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
													fMCparticleMother2 = (AliAODMCParticle*) fMCarray->At(fMCparticle2->GetMother());
													
													
													if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==443))
													{ 
														if(TMath::Abs(fMCparticle2->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother2->GetPdgCode())==443))
                                                        {
                                                            
                                                            
                                                            fPtMCparticle_Total_e_from_JPsi->Fill(track->Pt());
                                                            //leg2 passed emcal cuts. Leg1 is an electron from JPsi and is not on EMCal
                                                            fPtMCparticle_EMCalpid_leg2_e_from_JPsi->Fill(track2->Pt());
                                                            
                                                            //checking if they are from same mother
                                                            if((fMCparticleMother->GetLabel())==(fMCparticleMother2->GetLabel())){
                                                                //printf("electrons from same mother\n\n");
                                                                fPtMCparticle_Total_e_from_JPsi_sameMother->Fill(track->Pt()); //reconstructed pT
                                                                fPtMCparticle_Total_JPsi_pT->Fill(pt_kf);
                                                                
                                                                if(invmass3>=fMassCutMin && invmass3<=fMassCutMax){
                                                                    fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother->Fill(track->Pt()); //reconstructed pT
                                                                    fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother->Fill(pt_kf);//spectrum of reconstructed J/Psi
                                                                    
                                                                    //weights calculated based on J/Psi true MC pT, but applied on e+e- pair pt
                                                                    Double_t weight2 = CalculateWeight(fMCparticleMother->Pt());
                                                                    fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight->Fill(pt_kf, weight2);
                                                                    
                                                                    //to check for prompt J/psi -> check JPsi mother
                                                                    if(fMCparticleMother->GetMother()>0){
                                                                        fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
                                                                        Int_t gmpdg = fMCparticleGMother->GetPdgCode();
                                                                       // printf("mother of J/psi is %d\n", gmpdg);
                                                                        if(gmpdg > 500 && gmpdg < 600)fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight_prompt->Fill(pt_kf);//from B without weight
                                                                        else fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight_prompt->Fill(pt_kf, weight2);//prompt with weight (has mother but it is not B meson)
                                                                    }
                                                                    //if there is no mother, has to fill (prompt J/psi)
                                                                    if(fMCparticleMother->GetMother()<=0){
                                                                        //printf("case of J/psi without mother\n");
                                                                        fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight_prompt->Fill(pt_kf, weight2);//prompt with weight
                                                                    }
                                                                    
                                                                }//mass cut
                                                            }//same mother
                                                        }
													}
												}
											}//etacut second leg
										}
									}//eta cut
								}//close AOD
							}//close IsMC
						
                            //===============================

						 }
                      //}//close emcal trigger condition
						
					}
					
					if((fIsTrack1Emcal) && (fIsTrack2Emcal)){
                       // printf("Both tracks are  on EMCal! \n");
                        
                        //if(fEMCEG1 || fEMCEG2){
                            
                        // printf("The cuts are: %f , E/p < %f and E > %f\n",fEoverPCutMin, fEoverPCutMax, fEnergyCut);
                        
                        //printf("fClus->E() =%f,   fP= %f  \n",fClus->E(), fP);
                       // printf("fClus2->E() =%f,   fP2= %f  \n",fClus2->E(), fP2);
                        
                        //tpc electrons but at least one leg with track-matching
                        if(charge1*charge2 <0) fHist_InvMass_pt_ULStpc_wMatching->Fill(pt_kf,imass);
                        if(charge1*charge2 >0) fHist_InvMass_pt_LStpc_wMatching->Fill(pt_kf,imass);
                        
						
						if(((fClus->E() / fP) >=fEoverPCutMin && (fClus->E() / fP) <=fEoverPCutMax && (fClus->E()) >= fEnergyCut)||((fClus2->E()/fP2) >=fEoverPCutMin && (fClus2->E()/fP2) <=fEoverPCutMax && (fClus2->E()) >= fEnergyCut)){
							
                           // printf("One of them PASSED the cuts: \n");
                            
                            if(((fClus->E() / fP) >=fEoverPCutMin && (fClus->E() / fP) <=fEoverPCutMax && (fClus->E()) >= fEnergyCut)){
                                //printf("Track1 passed the cut and has pt=%f\n", fPt);
                            }
                            if(((fClus2->E()/fP2) >=fEoverPCutMin && (fClus2->E()/fP2) <=fEoverPCutMax && (fClus2->E()) >= fEnergyCut)){
                                // printf("Track2 passed the cut and has pt=%f\n", fPt2);
                            }
                            
						
							//if(charge1*charge2 <0) fHist_InvMass_pt_ULS->Fill(pt3,invmass3);
							//if(charge1*charge2 >0) fHist_InvMass_pt_LS->Fill(pt3,invmass3);
							
								//KFParticle
                            if(charge1*charge2 <0){
                                fHist_InvMass_pt_ULS_KF->Fill(pt_kf,imass);//multi integrated
                                if(fMultiAnalysis) fHist_InvMass_pt_ULS_KF_weight->Fill(pt_kf,imass,weight/weight2);//multi integrated with weight
                           
                                //correlation between leg1 and leg2
                                fHist_Correlation_leg1_emcal_leg2_emcal->Fill(fPt, fPt2);//not both above the threshold, at least one above the threshold
                            
                            }
							if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF->Fill(pt_kf,imass);
                            
                            //multiplicity bins histos (only ULS for SPDmulti and V0multi)
                            if(charge1*charge2 <0){
                                
                                if(fMultiAnalysis){
                                    // printf("both tracks: fSPDMult = %f, fVOMult = %f", fSPDMult,fV0Mult );
                                    if(fSPDMult_corr>0 && fSPDMult_corr < 10)   fHist_InvMass_pt_ULS_KF_SPDmulti_1->Fill(pt_kf,imass);
                                    if(fSPDMult_corr>=10 && fSPDMult_corr < 20) fHist_InvMass_pt_ULS_KF_SPDmulti_2->Fill(pt_kf,imass);
                                    if(fSPDMult_corr>=20 && fSPDMult_corr < 30) fHist_InvMass_pt_ULS_KF_SPDmulti_3->Fill(pt_kf,imass);
                                    if(fSPDMult_corr>=30 && fSPDMult_corr < 40) fHist_InvMass_pt_ULS_KF_SPDmulti_4->Fill(pt_kf,imass);
                                    if(fSPDMult_corr>=40 && fSPDMult_corr < 100) fHist_InvMass_pt_ULS_KF_SPDmulti_5->Fill(pt_kf,imass);
                                
                                    if(fV0Mult_corr2>0 && fV0Mult_corr2<100)     fHist_InvMass_pt_ULS_KF_V0multi_1->Fill(pt_kf,imass);
                                    if(fV0Mult_corr2>=100 && fV0Mult_corr2<200)  fHist_InvMass_pt_ULS_KF_V0multi_2->Fill(pt_kf,imass);
                                    if(fV0Mult_corr2>=200 && fV0Mult_corr2<300) fHist_InvMass_pt_ULS_KF_V0multi_3->Fill(pt_kf,imass);
                                    if(fV0Mult_corr2>=300 && fV0Mult_corr2<400) fHist_InvMass_pt_ULS_KF_V0multi_4->Fill(pt_kf,imass);
                                    if(fV0Mult_corr2>=400 && fV0Mult_corr2<800) fHist_InvMass_pt_ULS_KF_V0multi_5->Fill(pt_kf,imass);
                                
                                
                                    //with weight
                                    if(fSPDMult_corr>0 && fSPDMult_corr < 10)   fHist_InvMass_pt_ULS_KF_SPDmulti_1_weight->Fill(pt_kf,imass,weight/weight2);
                                    if(fSPDMult_corr>=10 && fSPDMult_corr < 20) fHist_InvMass_pt_ULS_KF_SPDmulti_2_weight->Fill(pt_kf,imass,weight/weight2);
                                    if(fSPDMult_corr>=20 && fSPDMult_corr < 30) fHist_InvMass_pt_ULS_KF_SPDmulti_3_weight->Fill(pt_kf,imass,weight/weight2);
                                    if(fSPDMult_corr>=30 && fSPDMult_corr < 40) fHist_InvMass_pt_ULS_KF_SPDmulti_4_weight->Fill(pt_kf,imass,weight/weight2);
                                    if(fSPDMult_corr>=40 && fSPDMult_corr < 100) fHist_InvMass_pt_ULS_KF_SPDmulti_5_weight->Fill(pt_kf,imass,weight/weight2);
                                
                                    if(fV0Mult_corr2>0 && fV0Mult_corr2<100)     fHist_InvMass_pt_ULS_KF_V0multi_1_weight->Fill(pt_kf,imass,weight/weight2);
                                    if(fV0Mult_corr2>=100 && fV0Mult_corr2<200)  fHist_InvMass_pt_ULS_KF_V0multi_2_weight->Fill(pt_kf,imass,weight/weight2);
                                    if(fV0Mult_corr2>=200 && fV0Mult_corr2<300) fHist_InvMass_pt_ULS_KF_V0multi_3_weight->Fill(pt_kf,imass,weight/weight2);
                                    if(fV0Mult_corr2>=300 && fV0Mult_corr2<400) fHist_InvMass_pt_ULS_KF_V0multi_4_weight->Fill(pt_kf,imass,weight/weight2);
                                    if(fV0Mult_corr2>=400 && fV0Mult_corr2<800) fHist_InvMass_pt_ULS_KF_V0multi_5_weight->Fill(pt_kf,imass,weight/weight2);
                                
                                }
                                
                            }
                            
							
								//both legs on emcal
							if(charge1*charge2 <0) fHist_InvMass_pt_ULSboth->Fill(pt_kf,imass);
							if(charge1*charge2 >0) fHist_InvMass_pt_LSboth->Fill(pt_kf,imass);
							
							fTracksPt[10]->Fill(fPt2);
							
								//===============================
							
							if(fIsMC)
							{
								
									//MC generators
								if(IsPythiaCC_gen){
									if(charge1*charge2 <0) fHist_InvMass_pt_ULS_KF_CC->Fill(pt_kf,imass);
									if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF_CC->Fill(pt_kf,imass);
								}
								if(IsPythiaBB_gen){
									if(charge1*charge2 <0) fHist_InvMass_pt_ULS_KF_BB->Fill(pt_kf,imass);
									if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF_BB->Fill(pt_kf,imass);
								}
								if(IsPythiaB_gen){
									if(charge1*charge2 <0) fHist_InvMass_pt_ULS_KF_B->Fill(pt_kf,imass);
									if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF_B->Fill(pt_kf,imass);
								}
								if(IsJpsi2ee_gen){
									if(charge1*charge2 <0) fHist_InvMass_pt_ULS_KF_Jpsi->Fill(pt_kf,imass);
									if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF_Jpsi->Fill(pt_kf,imass);
								}
								if(IsB2JPsi2ee_gen){
									if(charge1*charge2 <0) fHist_InvMass_pt_ULS_KF_BJpsi->Fill(pt_kf,imass);
									if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF_BJpsi->Fill(pt_kf,imass);
								}

								
								
								//Double_t fEtaCutMin=-0.9;
								//Double_t fEtaCutMax=0.9;
								
								//Double_t fMassCutMin=2.92;
								//Double_t fMassCutMax=3.16;
								
								if(fIsAOD)
								{
									fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
									fMCparticle2 = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track2->GetLabel()));
									
									Int_t pdg = fMCparticle->GetPdgCode();
									Int_t pdg3 = fMCparticle2->GetPdgCode();
									
									if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax && fMCparticle->Charge()!=0)
									{
										if(fMCparticle2->Eta()>=fEtaCutMin && fMCparticle2->Eta()<=fEtaCutMax && fMCparticle2->Charge()!=0){
											
											if( TMath::Abs(pdg) == 11 && TMath::Abs(pdg3) == 11 )
											{	
												
												if(fMCparticle->IsPhysicalPrimary() && fMCparticle2->IsPhysicalPrimary())
												{
													
														//For JPsi analysis
													if(fMCparticle->GetMother()<0 || fMCparticle2->GetMother()<0) return;
													
													fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
													fMCparticleMother2 = (AliAODMCParticle*) fMCarray->At(fMCparticle2->GetMother());
													
													
													if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==443))
													{ 
														if(TMath::Abs(fMCparticle2->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother2->GetPdgCode())==443))
                                                        {
                                                            
                                                            
                                                            fPtMCparticle_Total_e_from_JPsi->Fill(track->Pt());
                                                            //leg1 and leg2 on emcal. At least one of them passed the emcal cuts
                                                            fPtMCparticle_EMCalpid_both_leg1_e_from_JPsi->Fill(track->Pt());
                                                            fPtMCparticle_EMCalpid_both_leg2_e_from_JPsi->Fill(track2->Pt());
                                                            
                                                            //checking if they are from same mother
                                                            if((fMCparticleMother->GetLabel())==(fMCparticleMother2->GetLabel())){
                                                                //printf("electrons from same mother\n\n");
                                                                fPtMCparticle_Total_e_from_JPsi_sameMother->Fill(track->Pt()); //reconstructed pT
                                                                fPtMCparticle_Total_JPsi_pT->Fill(pt_kf);
                                                                
                                                                if(invmass3>=fMassCutMin && invmass3<=fMassCutMax){
                                                                    fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother->Fill(track->Pt()); //reconstructed pT
                                                                    fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother->Fill(pt_kf);//spectrum of reconstructed J/Psi
                                                                    
                                                                    //weights calculated based on J/Psi true MC pT, but applied on e+e- pair pt
                                                                    Double_t weight2 = CalculateWeight(fMCparticleMother->Pt());
                                                                    fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight->Fill(pt_kf, weight2);
                                                                    
                                                                    //to check for prompt J/psi -> check JPsi mother
                                                                    if(fMCparticleMother->GetMother()>0){
                                                                        fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
                                                                        Int_t gmpdg = fMCparticleGMother->GetPdgCode();
                                                                        //printf("mother of J/psi is %d\n", gmpdg);
                                                                        if(gmpdg > 500 && gmpdg < 600)fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight_prompt->Fill(pt_kf);//from B without weight
                                                                        else fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight_prompt->Fill(pt_kf, weight2);//prompt with weight (has mother but it is not B meson)
                                                                    }
                                                                    //if there is no mother, has to fill (prompt J/psi)
                                                                    if(fMCparticleMother->GetMother()<=0){
                                                                        //printf("case of J/psi without mother\n");
                                                                        fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight_prompt->Fill(pt_kf, weight2);//prompt with weight
                                                                    }
                                                                    
                                                                    
                                                                }//mass cut
                                                            }//same mother
                                                        }
													}
												}
											}//etacut second leg
										}
									}//eta cut
								}//close AOD
							}//close IsMC
							
							
							//===============================

						}
					 //}//close emcal trigger condition
                    }
					if((!fIsTrack1Emcal) && (!fIsTrack2Emcal)){
						//printf("Track1 and track2 are NOT on EMCal \n");
						
						//printf("In this case the inv mass is not filled \n");
						
						fTracksPt[11]->Fill(fPt2);

					}
                    //at least one leg on EMCal, but not necessarily above the threshold
                    if(fIsTrack1Emcal || fIsTrack2Emcal){
                        
                        if(fIsMC){
                            
                            if(fIsAOD)
                            {
                                fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
                                fMCparticle2 = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track2->GetLabel()));
                                
                                Int_t pdg = fMCparticle->GetPdgCode();
                                Int_t pdg3 = fMCparticle2->GetPdgCode();
                                
                                if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax && fMCparticle->Charge()!=0)
                                {
                                    if(fMCparticle2->Eta()>=fEtaCutMin && fMCparticle2->Eta()<=fEtaCutMax && fMCparticle2->Charge()!=0){
                                        
                                        if( TMath::Abs(pdg) == 11 && TMath::Abs(pdg3) == 11 )
                                        {
                                            
                                            if(fMCparticle->IsPhysicalPrimary() && fMCparticle2->IsPhysicalPrimary()){
                                                
                                                if(fIsTrack1Emcal)fPtMCparticle_EMCalpid_leg1->Fill(track->Pt());
                                                if(fIsTrack2Emcal)fPtMCparticle_EMCalpid_leg2->Fill(track2->Pt());
                                                
                                            }//Is Physical primary cut
                                        }//etacut second leg
                                    }
                                }//eta cut
                            }//close AOD
                        }
                    }//at least one leg on emcal, just to calculate EMCal PID and energy cut efficiency
					
				}
				
			}//close lTrack
		
	//////////////////////////////////////////////////////////////////////////////////////////////////
			
		//fVtxZ[2]->Fill(fZvtx);
		
		
		
        }
	
		
//=======================================================================
	
	delete fListOfmotherkink;
	PostData(1, fOutputList);
    //PostData(2,fListProfiles);
}      

//=======================================================================
void AliAnalysisTask_JPsi_EMCal::Terminate(Option_t *)
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

//=======================================================================
/*
 Bool_t AliAnalysisTask_JPsi_EMCal::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
//Check single track cuts for a given cut step
//Note this function is called inside the UserExec function
	const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
	if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
	return kTRUE;
}*/

//=======================================================================
Bool_t AliAnalysisTask_JPsi_EMCal::FindMother(Int_t mcIndex)
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
	//ESD part never used and then deleted
    
    
}
//____________________________________________________________________________
TProfile2D* AliAnalysisTask_JPsi_EMCal::GetEstimatorHistogram(const AliAODEvent* fAOD)
{
    
   // printf("Inside 'GetEstimatorHistogram \n'");
   // Int_t runNo  = fAOD->GetRunNumber();
    //cout<<"run number"<<runNo<<endl;
   Int_t period = -1;
    //Int_t group_number = -1;
    
    //period = 0;
    /*
    if (runNo>258883 && runNo<260187) period = 0;//16l  259668
    if (runNo>256504 && runNo<258574) period = 1;//16k
    if (period < 0 || period > 1) return 0;
    */
    /*
    //run groups from Steffen (period here means group of runs)
    if (runNo>0 && runNo<271868) group_number = 0;
    if (runNo>=271868 && runNo< 273591) group_number =  1;
    if (runNo>=273591 && runNo< 276644) group_number =  0;
    if (runNo>=276644 && runNo< 285009) group_number =  2;
    if (runNo>=285009 && runNo< 288743) group_number =  3;
    if (runNo>=288743 && runNo< 288861) group_number =  4;
    if (runNo>=288861 && runNo< 288943) group_number =  5;
    if (runNo>=288943 && runNo< 290549) group_number =  4;
    if (runNo>=290549 && runNo< 291285) group_number =  6;
    if (runNo>=291285 && runNo< 291416) group_number =  7;
    if (runNo>=291416 && runNo< 291690) group_number =  8;
    if (runNo>=291690 && runNo< 291942) group_number =  9;
    if (runNo>=291942 && runNo< 292012) group_number =  10;
    if (runNo>=292012 && runNo< 292163) group_number =  11;
    if (runNo>=292163 && runNo< 292265) group_number =  12;
    if (runNo>=292265 && runNo< 292809) group_number =  11;
    if (runNo>=292809 && runNo< 293357) group_number =  13;
    if (runNo>=293357) group_number =  4;
    */
    
    /* 0          --> 0
     271868     --> 1
     273591     --> 0
     276644     --> 2
     285009     --> 3
     288743     --> 4
     288861     --> 5
     288943     --> 4
     290549     --> 6
     291285     --> 7
     291416     --> 8
     291690     --> 9
     291942     --> 10
     292012     --> 11
     292163     --> 12
     292265     --> 11
     292809     --> 13
     293357     --> 4
     */
 
   // cout<<"using period = 0 for all (just a test)"<<period<<endl;
    period = 0;//same file for all periods
    
    return fMultEstimatorAvg[period];
}
//____________________________________________________________________________
TProfile2D* AliAnalysisTask_JPsi_EMCal::GetEstimatorHistogram_V0(const AliAODEvent* fAOD)
{
    
    //printf("Inside 'GetEstimatorHistogram_V0 \n'");
   // Int_t runNo  = fAOD->GetRunNumber();
    //cout<<"run number"<<runNo<<endl;
    Int_t period = -1;
    
    //period = 0;
    /*
    if (runNo>258883 && runNo<260187) period = 0;//16l  259668
    if (runNo>256504 && runNo<258574) period = 1;//16k
    if (period < 0 || period > 1) return 0;
     */
    
    
    
  //  cout<<"using period = 0 for all (just a test -- V0)"<<period<<endl;
    period = 0;
 
    return fMultEstimatorV0[period];
}
//______________________________________________________________________________
Double_t AliAnalysisTask_JPsi_EMCal::GetTrackletsMeanCorrection(TProfile2D* estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult, Int_t runNo)
{
    
   //printf("Inside 'GetTrackletsMeanCorrection' for run number %d and vertex =%f \n ", runNo, vtxZ);
    
    if(TMath::Abs(vtxZ)>10.0){
        //    printf("ERROR: Z vertex out of range for correction of multiplicity\n");
        return uncorrectedNacc;
    }
    
    if(!estimatorAvg){
        printf("ERROR: Missing TProfile for correction of multiplicity\n");
        return uncorrectedNacc;
    }
    Int_t group_number=-1;
    
    //run groups from Steffen
    if (runNo>0 && runNo<271868) group_number = 0;
    if (runNo>=271868 && runNo< 273591) group_number =  1;
    if (runNo>=273591 && runNo< 276644) group_number =  0;
    if (runNo>=276644 && runNo< 285009) group_number =  2;
    if (runNo>=285009 && runNo< 288743) group_number =  3;
    if (runNo>=288743 && runNo< 288861) group_number =  4;
    if (runNo>=288861 && runNo< 288943) group_number =  5;
    if (runNo>=288943 && runNo< 290549) group_number =  4;
    if (runNo>=290549 && runNo< 291285) group_number =  6;
    if (runNo>=291285 && runNo< 291416) group_number =  7;
    if (runNo>=291416 && runNo< 291690) group_number =  8;
    if (runNo>=291690 && runNo< 291942) group_number =  9;
    if (runNo>=291942 && runNo< 292012) group_number =  10;
    if (runNo>=292012 && runNo< 292163) group_number =  11;
    if (runNo>=292163 && runNo< 292265) group_number =  12;
    if (runNo>=292265 && runNo< 292809) group_number =  11;
    if (runNo>=292809 && runNo< 293357) group_number =  13;
    if (runNo>=293357) group_number =  4;
    
   // printf("Group number =%d\n", group_number);
    
    TH1F *hx0= (TH1F*) estimatorAvg->ProfileX();
    TH1F *hy0= (TH1F*) estimatorAvg->ProfileY();
    
    Int_t BinX=hx0->FindBin(group_number);
    Int_t BinY=hy0->FindBin(vtxZ);
    
   // printf("TEST2: binX = %d, binY=%d \n", BinX, BinY);
    
    
    delete hx0;
    delete hy0;
    
    
   
    
    Double_t localAvg = estimatorAvg->GetBinContent(BinX, BinY);
   // printf("LocalAvg = %f  \n", localAvg);
    
    if(localAvg==0){
        
        return uncorrectedNacc;
    }
    //printf("LocalAvg = %f for vertex= %f and group number =%d \n", localAvg, vtxZ, group_number);
    
    Double_t deltaM = 0;
    deltaM = uncorrectedNacc*(refMult/localAvg - 1);
    
    
   
    
    Double_t correctedNacc = uncorrectedNacc + (deltaM>0 ? 1 : -1) * gRandom->PoissonD(TMath::Abs(deltaM));
    
    if(correctedNacc<0) correctedNacc=0;
    
   // printf("Inside 'GetTrackletsMeanCorrection' after 2D correction: old = %f and new = %f \n",uncorrectedNacc, correctedNacc);
    
    return correctedNacc;
    
}

//______________________________________________________________________________
Double_t AliAnalysisTask_JPsi_EMCal::GetV0MeanCorrection(TProfile2D* estimatorV0, Double_t uncorrectedV0, Double_t vtxZ, Double_t refMult_V0, Int_t run_number)
{
    
    //printf("Inside GetV0MeanCorrection \n");
    if(TMath::Abs(vtxZ)>10.0){
        //    printf("ERROR: Z vertex out of range for correction of multiplicity\n");
        return uncorrectedV0;
    }
    
    if(!estimatorV0){
        printf("ERROR: Missing TProfile for correction of multiplicity\n");
        return uncorrectedV0;
    }
  
    //TH1F *hx= (TH1F*) estimatorV0->ProfileX();
    TH1F *hy= (TH1F*) estimatorV0->ProfileY();
    
    //profile as a function of integer from 0 to 1552. Finding relationship with run number:
    Int_t runs[1552]={252235,252248,252271,252310,252317,252319,252322,252325,252330,253437,253478,253481,253482,253488,253517,253529,253530,253563,253589,253591,254128,254147,254148,254149,254174,254175,254178,254193,254199,254204,254205,254293,254302,254303,254304,254330,254331,254332,254418,254419,254422,254604,254606,254608,254629,254630,254632,254640,254644,254646,254648,254649,254651,254652,254653,254654,255079,255082,255085,255086,255091,255111,255154,255159,255162,255167,255171,255173,255174,255176,255177,255240,255242,255247,255248,255249,255251,255252,255255,255256,255275,255276,255280,255283,255350,255351,255352,255398,255402,255407,255415,255418,255419,255420,255421,255440,255442,255447,255463,255465,255466,255467,255539,255540,255541,255542,255543,255577,255582,255583,255591,255614,255615,255616,255617,255618,256204,256207,256210,256212,256213,256215,256219,256222,256223,256227,256228,256231,256281,256282,256283,256284,256287,256289,256290,256292,256295,256297,256298,256299,256302,256307,256309,256311,256356,256357,256361,256362,256363,256364,256365,256366,256368,256371,256372,256373,256415,256417,256418,256420,256941,256942,256944,257011,257012,257021,257026,257028,257077,257080,257082,257083,257084,257086,257092,257095,257100,257136,257137,257139,257140,257141,257142,257144,257145,257204,257206,257209,257224,257260,257318,257320,257322,257330,257358,257364,257433,257457,257468,257474,257487,257488,257490,257491,257492,257530,257531,257537,257539,257540,257541,257560,257561,257562,257566,257587,257588,257590,257592,257594,257595,257601,257604,257605,257606,257630,257632,257635,257636,257642,257644,257682,257684,257685,257687,257688,257689,257691,257692,257694,257697,257724,257725,257727,257733,257734,257735,257737,257754,257757,257765,257773,257797,257798,257799,257800,257803,257804,257850,257851,257853,257855,257893,257936,257937,257939,257957,257958,257963,257979,257986,257989,257992,258003,258008,258012,258014,258017,258019,258039,258041,258042,258045,258048,258049,258053,258059,258060,258062,258063,258107,258108,258109,258113,258114,258117,258178,258197,258198,258202,258203,258204,258256,258257,258258,258270,258271,258273,258274,258278,258299,258301,258302,258303,258306,258307,258332,258336,258359,258387,258391,258393,258426,258452,258454,258456,258477,258499,258537,258962,258964,259088,259090,259091,259096,259099,259117,259118,259162,259164,259204,259257,259261,259263,259264,259269,259270,259271,259272,259273,259274,259302,259303,259305,259307,259334,259336,259339,259340,259341,259342,259378,259382,259388,259389,259394,259395,259396,259473,259477,259747,259748,259750,259751,259752,259756,259781,259788,259789,259822,259841,259842,259860,259866,259867,259868,259888,262424,262425,262426,262428,262705,262706,262708,262713,262717,262719,262723,262725,262727,262760,262768,262776,262777,262778,262841,262842,262844,262847,262849,262853,262855,262858,263487,263490,263496,263497,263529,263647,263652,263653,263654,263657,263662,263663,263682,263689,263690,263691,263737,263738,263739,263741,263743,263744,263784,263785,263786,263787,263790,263792,263793,263803,263810,263861,263863,263866,263905,263916,263917,263920,263923,263977,263978,263981,263984,263985,264033,264035,264076,264078,264082,264085,264086,264109,264110,264129,264137,264138,264139,264164,264168,264188,264190,264194,264197,264198,264232,264233,264235,264238,264259,264260,264261,264262,264264,264265,264266,264267,264273,264277,264279,264281,264305,264306,264312,264336,264341,264345,264346,264347,271868,271870,271871,271873,271874,271880,271881,271886,272018,272020,272036,272038,272039,272040,272041,272042,272076,272100,272101,272123,272151,272152,272153,272154,272155,272156,272194,272335,272340,272359,272360,272388,272389,272394,272395,272399,272400,272411,272413,272461,272462,272463,272466,272468,272521,272574,272575,272577,272585,272607,272608,272610,272620,272690,272691,272692,272712,272746,272747,272749,272760,272763,272764,272782,272783,272784,272828,272829,272833,272834,272836,272870,272871,272873,272880,272903,272905,272932,272933,272934,272935,272939,272947,272949,272976,272983,272985,273009,273010,273077,273099,273100,273103,273591,273592,273593,273653,273654,273687,273689,273690,273695,273709,273711,273719,273824,273825,273885,273886,273887,273889,273918,273942,273943,273946,273985,273986,274058,274064,274092,274094,274125,274147,274148,274174,274212,274232,274258,274259,274263,274264,274266,274268,274269,274270,274271,274276,274278,274280,274281,274283,274329,274351,274352,274360,274363,274364,274385,274386,274387,274388,274389,274390,274442,274593,274594,274595,274596,274601,274653,274657,274667,274669,274671,274690,274708,274801,274802,274803,274806,274807,274811,274815,274817,274821,274822,274877,274878,274882,274886,274978,274979,275067,275068,275073,275075,275076,275150,275151,275173,275174,275177,275180,275184,275188,275283,275314,275322,275324,275326,275328,275332,275333,275360,275361,275369,275372,275394,275395,275401,275404,275406,275443,275448,275452,275453,275456,275457,275459,275467,275471,275472,275515,275558,275559,275612,275617,275621,275622,275623,275624,275647,275648,275650,275657,275661,275664,275847,275924,275925,276012,276013,276017,276019,276020,276040,276041,276045,276097,276098,276099,276102,276104,276135,276140,276145,276166,276169,276170,276177,276178,276205,276230,276257,276259,276290,276291,276292,276294,276297,276302,276307,276312,276348,276351,276435,276437,276438,276439,276462,276506,276507,276508,276551,276552,276556,276557,276608,276644,276670,276671,276672,276674,276675,276762,276916,276917,276920,276967,276969,276970,276971,276972,277015,277016,277017,277037,277073,277076,277079,277082,277087,277091,277117,277121,277155,277180,277182,277183,277184,277188,277189,277193,277194,277196,277197,277256,277257,277262,277293,277310,277312,277314,277360,277383,277384,277385,277386,277389,277416,277417,277418,277472,277473,277476,277477,277478,277479,277530,277531,277534,277536,277537,277574,277575,277576,277577,277721,277722,277723,277725,277745,277746,277747,277749,277794,277795,277799,277800,277801,277802,277805,277834,277836,277841,277842,277845,277847,277848,277870,277876,277897,277898,277899,277900,277901,277903,277904,277907,277930,277952,277987,277989,277991,277996,278121,278122,278123,278126,278127,278158,278163,278164,278166,278167,278189,278191,278215,278216,278914,278915,278936,278939,278941,278959,278960,278963,278964,278999,279000,279005,279007,279008,279035,279036,279041,279043,279044,279068,279069,279073,279074,279075,279106,279107,279117,279118,279122,279123,279130,279155,279157,279199,279201,279207,279208,279232,279234,279235,279238,279242,279264,279265,279267,279268,279270,279273,279274,279309,279310,279312,279342,279344,279348,279349,279354,279355,279391,279410,279435,279439,279441,279483,279487,279488,279491,279550,279559,279630,279632,279641,279642,279676,279677,279679,279682,279683,279684,279687,279688,279689,279715,279718,279719,279747,279749,279773,279826,279827,279830,279853,279854,279879,279880,280051,280052,280066,280107,280108,280111,280114,280118,280126,280131,280134,280135,280140,280282,280283,280284,280285,280286,280290,280310,280312,280348,280349,280350,280351,280374,280375,280403,280406,280412,280413,280415,280419,280443,280445,280446,280447,280448,280490,280499,280518,280519,280546,280547,280550,280551,280574,280575,280576,280581,280583,280613,280634,280636,280637,280639,280645,280647,280671,280679,280681,280705,280706,280729,280753,280754,280755,280756,280757,280761,280762,280763,280764,280765,280766,280767,280768,280786,280787,280792,280793,280842,280844,280845,280847,280848,280849,280854,280856,280880,280881,280890,280897,280936,280940,280943,280947,280990,280994,280996,280997,280998,280999,281032,281033,281035,281036,281060,281061,281062,281080,281081,281179,281180,281181,281189,281190,281191,281212,281213,281240,281241,281242,281243,281244,281271,281273,281275,281277,281301,281321,281415,281441,281443,281444,281446,281449,281450,281475,281477,281509,281511,281557,281562,281563,281568,281569,281574,281583,281592,281633,281892,281893,281894,281895,281915,281916,281918,281920,281928,281931,281932,281939,281940,281953,281956,281961,282528,282544,282545,282546,282573,282575,282579,282580,282606,282607,282608,282609,282618,282620,282622,282629,282651,282653,282666,282667,282668,282670,282671,282673,282676,282677,282700,282702,282703,282704,285978,285979,285980,286014,286025,286064,286124,286127,286129,286130,286159,286198,286199,286201,286202,286203,286229,286230,286231,286254,286255,286257,286258,286261,286263,286282,286284,286287,286288,286289,286308,286309,286310,286311,286312,286314,286336,286337,286340,286341,286345,286348,286349,286350,286380,286426,286427,286428,286454,286455,286482,286501,286502,286508,286509,286511,286566,286567,286568,286569,286591,286592,286594,286633,286653,286661,286695,286731,286799,286801,286805,286809,286810,286846,286848,286850,286852,286874,286876,286877,286907,286908,286910,286911,286930,286931,286932,286933,286936,286937,288861,288862,288863,288864,288868,288902,288903,288908,288909,290323,290324,290327,290350,290374,290375,290376,290399,290401,290404,290411,290412,290423,290425,290426,290427,290428,290456,290458,290459,290467,290469,290499,290500,290501,290549,290550,290553,290588,290590,290627,290632,290645,290658,290660,290665,290687,290689,290692,290696,290699,290721,290764,290766,290769,290774,290776,290787,290790,290841,290843,290846,290848,290860,290862,290886,290887,290888,290892,290894,290895,290932,290935,290941,290943,290944,290948,290974,290975,290976,290979,290980,291002,291003,291004,291005,291006,291035,291037,291041,291065,291066,291069,291093,291100,291101,291110,291111,291116,291143,291188,291209,291240,291257,291263,291265,291266,291282,291283,291284,291285,291286,291360,291361,291362,291363,291375,291377,291397,291399,291400,291402,291416,291417,291419,291420,291424,291446,291447,291451,291453,291456,291457,291481,291482,291484,291485,291590,291614,291615,291618,291622,291624,291626,291690,291692,291694,291697,291698,291702,291706,291729,291755,291756,291760,291762,291766,291768,291769,291795,291796,291803,291942,291944,291945,291946,291948,291953,291976,291977,291982,292012,292040,292060,292061,292062,292067,292075,292077,292080,292081,292106,292107,292108,292109,292114,292115,292140,292160,292161,292162,292163,292164,292166,292167,292168,292192,292218,292240,292241,292242,292265,292269,292270,292273,292274,292298,292397,292398,292405,292406,292428,292429,292430,292432,292434,292456,292457,292460,292461,292495,292496,292497,292500,292521,292523,292524,292526,292553,292554,292557,292559,292560,292563,292584,292586,292693,292695,292696,292698,292701,292704,292737,292739,292744,292747,292748,292750,292752,292754,292758,292803,292804,292809,292810,292811,292831,292832,292834,292836,292839};
    
   
    Int_t BinX=-99;
    for(Int_t i=0;i<1552;i++){
        if(runs[i] == run_number){
           // printf("==============================  (%d)  position for runs[%d]= %d  is %d\n",run_number, i,runs[i], i);
            BinX=i;
        }
    }
   // printf("TEST_V0: binx for run %d is %d \n", run_number, BinX);
    
    //Int_t BinX=hx->FindBin(run_number);
    Int_t BinY=hy->FindBin(vtxZ);
    
    //printf("run= %d, vtx=%f ==>   Bins inside task are binx = %d, biny = %d \n", run_number, vtxZ, BinX, BinY);
    
    //delete hx;
    delete hy;
    
    
    Double_t localV0 = estimatorV0->GetBinContent(BinX, BinY);//first argument is run number/group, second is vertex on Y axis
   // printf("TEST_V0: localV0 = %f\n", localV0);
    
    if(localV0==0){
       // printf("LocalV0 = 0 for vertex = %f and  run =  %d \n", vtxZ, run_number);
        return uncorrectedV0;
    }
    
    Double_t deltaMV0 = 0;
    deltaMV0 = uncorrectedV0*(refMult_V0/localV0 - 1);
    
    //printf("Inside 'GetV0MeanCorrection' (after deltaM) \n");
    
    Double_t correctedV0 = uncorrectedV0 + (deltaMV0>0 ? 1 : -1) * gRandom_V0->PoissonD(TMath::Abs(deltaMV0));
    
    if(correctedV0<0) correctedV0=0;
    
  // printf("Inside 'GetV0MeanCorrection' (end) \n");
    
    return correctedV0;
    
}
Double_t AliAnalysisTask_JPsi_EMCal::CalculateWeight(Double_t x)
{
    Double_t weight=1;
    
    if(x>= 0.000 &&  x < 0.200 ) weight=0.614277;
    if(x>= 0.200 &&  x < 0.400 ) weight=0.613708;
    if(x>= 0.400 &&  x < 0.600 ) weight=0.616810;
    if(x>= 0.600 &&  x < 0.800 ) weight=0.621853;
    if(x>= 0.800 &&  x < 1.000 ) weight=0.630272;
    if(x>= 1.000 &&  x < 1.200 ) weight=0.639764;
    if(x>= 1.200 &&  x < 1.400 ) weight=0.650248;
    if(x>= 1.400 &&  x < 1.600 ) weight=0.662680;
    if(x>= 1.600 &&  x < 1.800 ) weight=0.674690;
    if(x>= 1.800 &&  x < 2.000 ) weight=0.690799;
    if(x>= 2.000 &&  x < 2.200 ) weight=0.706232;
    if(x>= 2.200 &&  x < 2.400 ) weight=0.722133;
    if(x>= 2.400 &&  x < 2.600 ) weight=0.740174;
    if(x>= 2.600 &&  x < 2.800 ) weight=0.759401;
    if(x>= 2.800 &&  x < 3.000 ) weight=0.776361;
    if(x>= 3.000 &&  x < 3.200 ) weight=0.791707;
    if(x>= 3.200 &&  x < 3.400 ) weight=0.812156;
    if(x>= 3.400 &&  x < 3.600 ) weight=0.827268;
    if(x>= 3.600 &&  x < 3.800 ) weight=0.844515;
    if(x>= 3.800 &&  x < 4.000 ) weight=0.861771;
    if(x>= 4.000 &&  x < 4.200 ) weight=0.877448;
    if(x>= 4.200 &&  x < 4.400 ) weight=0.893644;
    if(x>= 4.400 &&  x < 4.600 ) weight=0.911980;
    if(x>= 4.600 &&  x < 4.800 ) weight=0.923118;
    if(x>= 4.800 &&  x < 5.000 ) weight=0.940259;
    if(x>= 5.000 &&  x < 5.200 ) weight=0.952314;
    if(x>= 5.200 &&  x < 5.400 ) weight=0.967387;
    if(x>= 5.400 &&  x < 5.600 ) weight=0.980388;
    if(x>= 5.600 &&  x < 5.800 ) weight=0.994637;
    if(x>= 5.800 &&  x < 6.000 ) weight=1.000000;
    if(x>= 6.000 &&  x < 6.200 ) weight=0.674078;
    if(x>= 6.200 &&  x < 6.400 ) weight=0.648798;
    if(x>= 6.400 &&  x < 6.600 ) weight=0.624233;
    if(x>= 6.600 &&  x < 6.800 ) weight=0.595486;
    if(x>= 6.800 &&  x < 7.000 ) weight=0.563306;
    if(x>= 7.000 &&  x < 7.200 ) weight=0.528029;
    if(x>= 7.200 &&  x < 7.400 ) weight=0.496518;
    if(x>= 7.400 &&  x < 7.600 ) weight=0.469522;
    if(x>= 7.600 &&  x < 7.800 ) weight=0.435593;
    if(x>= 7.800 &&  x < 8.000 ) weight=0.401749;
    if(x>= 8.000 &&  x < 8.200 ) weight=0.373716;
    if(x>= 8.200 &&  x < 8.400 ) weight=0.346930;
    if(x>= 8.400 &&  x < 8.600 ) weight=0.315269;
    if(x>= 8.600 &&  x < 8.800 ) weight=0.292034;
    if(x>= 8.800 &&  x < 9.000 ) weight=0.269475;
    if(x>= 9.000 &&  x < 9.200 ) weight=0.245979;
    if(x>= 9.200 &&  x < 9.400 ) weight=0.223318;
    if(x>= 9.400 &&  x < 9.600 ) weight=0.205727;
    if(x>= 9.600 &&  x < 9.800 ) weight=0.185729;
    if(x>= 9.800 &&  x < 10.000 ) weight=0.171318;
    if(x>= 10.000 &&  x < 10.200 ) weight=0.153972;
    if(x>= 10.200 &&  x < 10.400 ) weight=0.140388;
    if(x>= 10.400 &&  x < 10.600 ) weight=0.128839;
    if(x>= 10.600 &&  x < 10.800 ) weight=0.116666;
    if(x>= 10.800 &&  x < 11.000 ) weight=0.105318;
    if(x>= 11.000 &&  x < 11.200 ) weight=0.097047;
    if(x>= 11.200 &&  x < 11.400 ) weight=0.088199;
    if(x>= 11.400 &&  x < 11.600 ) weight=0.079519;
    if(x>= 11.600 &&  x < 11.800 ) weight=0.074030;
    if(x>= 11.800 &&  x < 12.000 ) weight=0.067971;
    if(x>= 12.000 &&  x < 12.200 ) weight=0.060773;
    if(x>= 12.200 &&  x < 12.400 ) weight=0.054765;
    if(x>= 12.400 &&  x < 12.600 ) weight=0.050744;
    if(x>= 12.600 &&  x < 12.800 ) weight=0.046101;
    if(x>= 12.800 &&  x < 13.000 ) weight=0.042468;
    if(x>= 13.000 &&  x < 13.200 ) weight=0.038402;
    if(x>= 13.200 &&  x < 13.400 ) weight=0.034957;
    if(x>= 13.400 &&  x < 13.600 ) weight=0.031370;
    if(x>= 13.600 &&  x < 13.800 ) weight=0.028945;
    if(x>= 13.800 &&  x < 14.000 ) weight=0.026552;
    if(x>= 14.000 &&  x < 14.200 ) weight=0.024333;
    if(x>= 14.200 &&  x < 14.400 ) weight=0.023281;
    if(x>= 14.400 &&  x < 14.600 ) weight=0.021372;
    if(x>= 14.600 &&  x < 14.800 ) weight=0.019601;
    if(x>= 14.800 &&  x < 15.000 ) weight=0.017743;
    if(x>= 15.000 &&  x < 15.200 ) weight=0.017142;
    if(x>= 15.200 &&  x < 15.400 ) weight=0.015805;
    if(x>= 15.400 &&  x < 15.600 ) weight=0.014111;
    if(x>= 15.600 &&  x < 15.800 ) weight=0.013415;
    if(x>= 15.800 &&  x < 16.000 ) weight=0.011595;
    if(x>= 16.000 &&  x < 16.200 ) weight=0.011130;
    if(x>= 16.200 &&  x < 16.400 ) weight=0.010443;
    if(x>= 16.400 &&  x < 16.600 ) weight=0.010265;
    if(x>= 16.600 &&  x < 16.800 ) weight=0.009009;
    if(x>= 16.800 &&  x < 17.000 ) weight=0.008355;
    if(x>= 17.000 &&  x < 17.200 ) weight=0.007828;
    if(x>= 17.200 &&  x < 17.400 ) weight=0.007548;
    if(x>= 17.400 &&  x < 17.600 ) weight=0.006651;
    if(x>= 17.600 &&  x < 17.800 ) weight=0.006077;
    if(x>= 17.800 &&  x < 18.000 ) weight=0.006180;
    if(x>= 18.000 &&  x < 18.200 ) weight=0.005638;
    if(x>= 18.200 &&  x < 18.400 ) weight=0.005469;
    if(x>= 18.400 &&  x < 18.600 ) weight=0.004448;
    if(x>= 18.600 &&  x < 18.800 ) weight=0.004038;
    if(x>= 18.800 &&  x < 19.000 ) weight=0.003919;
    if(x>= 19.000 &&  x < 19.200 ) weight=0.003747;
    if(x>= 19.200 &&  x < 19.400 ) weight=0.003555;
    if(x>= 19.400 &&  x < 19.600 ) weight=0.003505;
    if(x>= 19.600 &&  x < 19.800 ) weight=0.003087;
    if(x>= 19.800 &&  x < 20.000 ) weight=0.002941;
    if(x>= 20.000 &&  x < 20.200 ) weight=0.002754;
    if(x>= 20.200 &&  x < 20.400 ) weight=0.002387;
    if(x>= 20.400 &&  x < 20.600 ) weight=0.002710;
    if(x>= 20.600 &&  x < 20.800 ) weight=0.002113;
    if(x>= 20.800 &&  x < 21.000 ) weight=0.002486;
    if(x>= 21.000 &&  x < 21.200 ) weight=0.002096;
    if(x>= 21.200 &&  x < 21.400 ) weight=0.002194;
    if(x>= 21.400 &&  x < 21.600 ) weight=0.001721;
    if(x>= 21.600 &&  x < 21.800 ) weight=0.001628;
    if(x>= 21.800 &&  x < 22.000 ) weight=0.001981;
    if(x>= 22.000 &&  x < 22.200 ) weight=0.001892;
    if(x>= 22.200 &&  x < 22.400 ) weight=0.001488;
    if(x>= 22.400 &&  x < 22.600 ) weight=0.001412;
    if(x>= 22.600 &&  x < 22.800 ) weight=0.001390;
    if(x>= 22.800 &&  x < 23.000 ) weight=0.001336;
    if(x>= 23.000 &&  x < 23.200 ) weight=0.001326;
    if(x>= 23.200 &&  x < 23.400 ) weight=0.001355;
    if(x>= 23.400 &&  x < 23.600 ) weight=0.001102;
    if(x>= 23.600 &&  x < 23.800 ) weight=0.001032;
    if(x>= 23.800 &&  x < 24.000 ) weight=0.001129;
    if(x>= 24.000 &&  x < 24.200 ) weight=0.001078;
    if(x>= 24.200 &&  x < 24.400 ) weight=0.000867;
    if(x>= 24.400 &&  x < 24.600 ) weight=0.000866;
    if(x>= 24.600 &&  x < 24.800 ) weight=0.000722;
    if(x>= 24.800 &&  x < 25.000 ) weight=0.000862;
    if(x>= 25.000 &&  x < 25.200 ) weight=0.000818;
    if(x>= 25.200 &&  x < 25.400 ) weight=0.001302;
    if(x>= 25.400 &&  x < 25.600 ) weight=0.000799;
    if(x>= 25.600 &&  x < 25.800 ) weight=0.000603;
    if(x>= 25.800 &&  x < 26.000 ) weight=0.000619;
    if(x>= 26.000 &&  x < 26.200 ) weight=0.000766;
    if(x>= 26.200 &&  x < 26.400 ) weight=0.000535;
    if(x>= 26.400 &&  x < 26.600 ) weight=0.000420;
    if(x>= 26.600 &&  x < 26.800 ) weight=0.000450;
    if(x>= 26.800 &&  x < 27.000 ) weight=0.000365;
    if(x>= 27.000 &&  x < 27.200 ) weight=0.000488;
    if(x>= 27.200 &&  x < 27.400 ) weight=0.000347;
    if(x>= 27.400 &&  x < 27.600 ) weight=0.000386;
    if(x>= 27.600 &&  x < 27.800 ) weight=0.000325;
    if(x>= 27.800 &&  x < 28.000 ) weight=0.000320;
    if(x>= 28.000 &&  x < 28.200 ) weight=0.000260;
    if(x>= 28.200 &&  x < 28.400 ) weight=0.000237;
    if(x>= 28.400 &&  x < 28.600 ) weight=0.000251;
    if(x>= 28.600 &&  x < 28.800 ) weight=0.000382;
    if(x>= 28.800 &&  x < 29.000 ) weight=0.000187;
    if(x>= 29.000 &&  x < 29.200 ) weight=0.000317;
    if(x>= 29.200 &&  x < 29.400 ) weight=0.000168;
    if(x>= 29.400 &&  x < 29.600 ) weight=0.000182;
    if(x>= 29.600 &&  x < 29.800 ) weight=0.000213;
    if(x>= 29.800 &&  x < 30.000 ) weight=0.000332;
    if(x>= 30.000 &&  x < 30.200 ) weight=0.000169;
    if(x>= 30.200 &&  x < 30.400 ) weight=0.000188;
    if(x>= 30.400 &&  x < 30.600 ) weight=0.000135;
    if(x>= 30.600 &&  x < 30.800 ) weight=0.000196;
    if(x>= 30.800 &&  x < 31.000 ) weight=0.000126;
    if(x>= 31.000 &&  x < 31.200 ) weight=0.000142;
    if(x>= 31.200 &&  x < 31.400 ) weight=0.000305;
    if(x>= 31.400 &&  x < 31.600 ) weight=0.000084;
    if(x>= 31.600 &&  x < 31.800 ) weight=0.000139;
    if(x>= 31.800 &&  x < 32.000 ) weight=0.000466;
    if(x>= 32.000 &&  x < 32.200 ) weight=0.000129;
    if(x>= 32.200 &&  x < 32.400 ) weight=0.000132;
    if(x>= 32.400 &&  x < 32.600 ) weight=0.000065;
    if(x>= 32.600 &&  x < 32.800 ) weight=0.000093;
    if(x>= 32.800 &&  x < 33.000 ) weight=0.000133;
    if(x>= 33.000 &&  x < 33.200 ) weight=0.000086;
    if(x>= 33.200 &&  x < 33.400 ) weight=0.000099;
    if(x>= 33.400 &&  x < 33.600 ) weight=0.000067;
    if(x>= 33.600 &&  x < 33.800 ) weight=0.000135;
    if(x>= 33.800 &&  x < 34.000 ) weight=0.000089;
    if(x>= 34.000 &&  x < 34.200 ) weight=0.000073;
    if(x>= 34.200 &&  x < 34.400 ) weight=0.000090;
    if(x>= 34.400 &&  x < 34.600 ) weight=0.000050;
    if(x>= 34.600 &&  x < 34.800 ) weight=0.000076;
    if(x>= 34.800 &&  x < 35.000 ) weight=0.000142;
    if(x>= 35.000 &&  x < 35.200 ) weight=0.000098;
    if(x>= 35.200 &&  x < 35.400 ) weight=0.000186;
    if(x>= 35.400 &&  x < 35.600 ) weight=0.000068;
    if(x>= 35.600 &&  x < 35.800 ) weight=0.000039;
    if(x>= 35.800 &&  x < 36.000 ) weight=0.000097;
    if(x>= 36.000 &&  x < 36.200 ) weight=0.000068;
    if(x>= 36.200 &&  x < 36.400 ) weight=0.000084;
    if(x>= 36.400 &&  x < 36.600 ) weight=0.000039;
    if(x>= 36.600 &&  x < 36.800 ) weight=0.000072;
    if(x>= 36.800 &&  x < 37.000 ) weight=0.000029;
    if(x>= 37.000 &&  x < 37.200 ) weight=0.000095;
    if(x>= 37.200 &&  x < 37.400 ) weight=0.000044;
    if(x>= 37.400 &&  x < 37.600 ) weight=0.000033;
    if(x>= 37.600 &&  x < 37.800 ) weight=0.000038;
    if(x>= 37.800 &&  x < 38.000 ) weight=0.000024;
    if(x>= 38.000 &&  x < 38.200 ) weight=0.000033;
    if(x>= 38.200 &&  x < 38.400 ) weight=0.000024;
    if(x>= 38.400 &&  x < 38.600 ) weight=0.000019;
    if(x>= 38.600 &&  x < 38.800 ) weight=0.000021;
    if(x>= 38.800 &&  x < 39.000 ) weight=0.000014;
    if(x>= 39.000 &&  x < 39.200 ) weight=0.000016;
    if(x>= 39.200 &&  x < 39.400 ) weight=0.000011;
    if(x>= 39.400 &&  x < 39.600 ) weight=0.000014;
    if(x>= 39.600 &&  x < 39.800 ) weight=0.000019;
    if(x>= 39.800 &&  x < 40.000 ) weight=0.000024;
    if(x>= 40.000 &&  x < 40.200 ) weight=0.000017;
    if(x>= 40.200 &&  x < 40.400 ) weight=0.000010;
    if(x>= 40.400 &&  x < 40.600 ) weight=0.000058;
    if(x>= 40.600 &&  x < 40.800 ) weight=0.000011;
    if(x>= 40.800 &&  x < 41.000 ) weight=0.000007;
    if(x>= 41.000 &&  x < 41.200 ) weight=0.000010;
    if(x>= 41.200 &&  x < 41.400 ) weight=0.000011;
    if(x>= 41.400 &&  x < 41.600 ) weight=0.000008;
    if(x>= 41.600 &&  x < 41.800 ) weight=0.000006;
    if(x>= 41.800 &&  x < 42.000 ) weight=0.000006;
    if(x>= 42.000 &&  x < 42.200 ) weight=0.000017;
    if(x>= 42.200 &&  x < 42.400 ) weight=0.000016;
    if(x>= 42.400 &&  x < 42.600 ) weight=0.000004;
    if(x>= 42.600 &&  x < 42.800 ) weight=0.000014;
    if(x>= 42.800 &&  x < 43.000 ) weight=0.000005;
    if(x>= 43.000 &&  x < 43.200 ) weight=0.000013;
    if(x>= 43.200 &&  x < 43.400 ) weight=0.000012;
    if(x>= 43.400 &&  x < 43.600 ) weight=0.000004;
    if(x>= 43.600 &&  x < 43.800 ) weight=0.000003;
    if(x>= 43.800 &&  x < 44.000 ) weight=0.000007;
    if(x>= 44.000 &&  x < 44.200 ) weight=0.000003;
    if(x>= 44.200 &&  x < 44.400 ) weight=0.000006;
    if(x>= 44.400 &&  x < 44.600 ) weight=0.000002;
    if(x>= 44.600 &&  x < 44.800 ) weight=0.000005;
    if(x>= 44.800 &&  x < 45.000 ) weight=0.000002;
    if(x>= 45.000 &&  x < 45.200 ) weight=0.000003;
    if(x>= 45.200 &&  x < 45.400 ) weight=0.000002;
    if(x>= 45.400 &&  x < 45.600 ) weight=0.000002;
    if(x>= 45.600 &&  x < 45.800 ) weight=0.000002;
    if(x>= 45.800 &&  x < 46.000 ) weight=0.000005;
    if(x>= 46.000 &&  x < 46.200 ) weight=0.000003;
    if(x>= 46.200 &&  x < 46.400 ) weight=0.000002;
    if(x>= 46.400 &&  x < 46.600 ) weight=0.000009;
    if(x>= 46.600 &&  x < 46.800 ) weight=0.000002;
    if(x>= 46.800 &&  x < 47.000 ) weight=0.000008;
    if(x>= 47.000 &&  x < 47.200 ) weight=0.000003;
    if(x>= 47.200 &&  x < 47.400 ) weight=0.000002;
    if(x>= 47.400 &&  x < 47.600 ) weight=0.000001;
    if(x>= 47.600 &&  x < 47.800 ) weight=0.000001;
    if(x>= 47.800 &&  x < 48.000 ) weight=0.000001;
    if(x>= 48.000 &&  x < 48.200 ) weight=0.000006;
    if(x>= 48.200 &&  x < 48.400 ) weight=0.000001;
    if(x>= 48.400 &&  x < 48.600 ) weight=0.000001;
    if(x>= 48.600 &&  x < 48.800 ) weight=0.000000;
    if(x>= 48.800 &&  x < 49.000 ) weight=0.000002;
    if(x>= 49.000 &&  x < 49.200 ) weight=0.000001;
    if(x>= 49.200 &&  x < 49.400 ) weight=0.000000;
    if(x>= 49.400 &&  x < 49.600 ) weight=0.000000;
    if(x>= 49.600 &&  x < 49.800 ) weight=0.000002;
    if(x>= 49.800 &&  x < 50.000 ) weight=0.000000;
    
    return weight;
}
//______________________________________________________________________________
Double_t AliAnalysisTask_JPsi_EMCal::GetTPCCalibration(Int_t runNo, Double_t TPCnsigma0)
{
   
    Double_t mean_shift=0.00;
    Double_t sigma_norm=1.00;
    
    
    if(runNo == 258454){//16k test
        mean_shift = 0.2;
        sigma_norm = 1;
     
    }
    
    if (runNo>271839 && runNo<273103){//17h
        mean_shift = 0.08;
        sigma_norm = 1.03;
    }
    if (runNo>=273486 && runNo< 274442){//17i
        mean_shift = 0.15;
        sigma_norm = 1.03;
    }
    if (runNo>=274690 && runNo< 276508){//17k
        mean_shift = 0.15;
        sigma_norm = 1.03;
    }
    if (runNo>=276551 && runNo< 278729){//17l
        mean_shift = 0.15;
        sigma_norm = 1.03;
    }
    if (runNo>=278818 && runNo< 280140){//17m
        mean_shift = 0.15;
        sigma_norm = 1.03;
    }
    if (runNo>=280282 && runNo< 281961){//17o
        mean_shift = 0.25;
        sigma_norm = 1.03;
    }
    if (runNo>=282504 && runNo< 282704){//17r
        mean_shift = 0.15;
        sigma_norm = 1.03;
    }
    
    //2018 periods
    if (runNo>=286982 && runNo< 287977){//18f
        mean_shift = -0.05;
        sigma_norm = 0.8;
    }
    if (runNo>=288804 && runNo< 288806){//18h
        mean_shift = -0.05;
        sigma_norm = 0.8;
    }
    if (runNo>=289165 && runNo< 289201){//18k
        mean_shift = -0.05;
        sigma_norm = 0.8;
    }
    if (runNo>=289240 && runNo< 289971){//18l
        mean_shift = -0.05;
        sigma_norm = 0.8;
    }
    if (runNo>=293368 && runNo< 293898){//18o
        mean_shift = -0.05;
        sigma_norm = 0.8;
    }
    if (runNo>=294009 && runNo< 295232){//18p
        mean_shift = -0.05;
        sigma_norm = 0.8;
    }
    
    
    //corrected TPCnsigma
    Double_t TPCnsigma_corr = (TPCnsigma0 - mean_shift)/sigma_norm;
    
   
    return TPCnsigma_corr;
    
}


/*
Bool_t AliAnalysisTask_JPsi_EMCal::TrackCuts(AliVTrack *track, AliAODTrack *atrack, AliESDtrack *etrack, Bool_t fIsAOD){
	
		//AOD (Test Filter Bit)
	if(fIsAOD)
	{
		
		if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE; 
	}
	
		//RecKine: ITSTPC cuts  
	if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) return kFALSE;
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
					return kFALSE;
				}
			}
			if(!kinkmotherpass) return kFALSE;
		}
		else
		{
			if(etrack->GetKinkIndex(0) != 0) return kFALSE;
		}
	} 
    
		//RecPrim
	
	if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) return kFALSE; //applies DCA cut!!!
	
    
		//HFEcuts: ITS layers cuts
	if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) return kFALSE;
    
		//HFE cuts: TPC PID cleanup
	if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) return kFALSE;
	
	return kTRUE;

}
 */

