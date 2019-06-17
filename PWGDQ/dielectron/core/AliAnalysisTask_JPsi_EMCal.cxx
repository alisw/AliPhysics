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
//		Cristiane Jahnke		(cristiane.jahnke@cern.ch)		      //
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
#include "AliGenEventHeader.h"
#include "AliTrackerBase.h"
#include "AliAODVZERO.h"
#include "AliAODTracklets.h"
#include "AliESDUtils.h"

//______________________________________________________________________

//______________________________________________________________________
ClassImp(AliAnalysisTask_JPsi_EMCal)

//______________________________________________________________________
AliAnalysisTask_JPsi_EMCal::AliAnalysisTask_JPsi_EMCal(const char *name)
  : AliAnalysisTaskSE(name)

,fIsMC(0)
,fUseTender(kFALSE)

//new Tender organization, using global variables
,fTenderClusterName("caloClusters")
,fTenderTrackName("tracks")
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
,fTPCnclsPID(85)
,fTPCchi2(4)
,fDCAxyCut(1)
,fDCAzCut(3)

,fTPCnsigmaCutMin(-2.25)
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
,gRandom(new TRandom3(0))

,fRefMult_V0(86.0)
,gRandom_V0(new TRandom3(0))


,fClus(0)
,fClus2(0)
,fClusAOD(0)

//Histograms for the analysis
,fNevent(0)
,fPDG_values(0)
,fNevent_SPD_multi(0)
,fNevent_V0_multi(0)
,fEoverP_pt(0)
,fTPC_p(0)
,fTPCnsigma_p(0)

,fTOF_p(0)
,fTOFnsigma_p(0)


,fTPCnsigma_EoverP(0)
,fECluster(0)
,fECluster_emcal(0)
,fECluster_dcal(0)
,fTracksPt(0)
,fTracksQAPt(0)
,fVtxZ(0)
//histos for SPD and V0 multiplicity
,fVtxZ_V0(0)
,fVtxZ_SPD(0)
,fV0_SPD(0)
,fV0(0)
,fSPD(0)


,fNClusters(0)
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
,fCuts(0)
//,fCFM(0)
//,fPID(new AliHFEpid("hfePid"))
//,fPIDqa(0)

//For MC
,fMCstack(0)

,fMCtrack(0)
,fMCtrackMother(0)
,fMCtrackGMother(0)
,fMCtrackGGMother(0)
,fMCtrackGGGMother(0)
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

	//new histos
,fdEta_dPhi(0)


,fSparseElectron(0)
,fvalueElectron(0)
,fSparseMulti(0)
,fvalueMulti(0)
,fIspp2011(kFALSE)

	//MC efficiencies
,fPtMCparticleAllHfe1(0)
,fPtMCparticleRecoHfe1(0)
,fPtMCparticleAll_e_from_JPsi(0)
,fPtMCparticleAll_JPsi_pT(0)
,fPtMCparticleAll_trueJPsi_pT(0)
,fPtMCparticleReco_e_from_JPsi(0)
,fPtMCparticle_Total_e_from_JPsi(0)
,fPtMCparticle_Total_e_from_JPsi_sameMother(0)
,fPtMCparticle_TotalplusMass_e_from_JPsi(0)
,fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother(0)

,fPtMCparticle_TotalplusMass_JPsi_pT(0)
,fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother(0)

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

//new Tender organization, uisng global variables
,fTenderClusterName("caloClusters")
,fTenderTrackName("tracks")
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
,fTPCnclsPID(85)
,fTPCchi2(4)
,fDCAxyCut(1)
,fDCAzCut(3)

,fTPCnsigmaCutMin(-2.25)
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
,gRandom(new TRandom3(0))

,fRefMult_V0(86.0)
,gRandom_V0(new TRandom3(0))

,fClus(0)
,fClus2(0)
,fClusAOD(0)

//Histograms for the analysis
,fNevent(0)
,fPDG_values(0)
,fNevent_SPD_multi(0)
,fNevent_V0_multi(0)
,fEoverP_pt(0)
,fTPC_p(0)
,fTPCnsigma_p(0)

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
,fVtxZ(0)
//histos for SPD and V0 multiplicity
,fVtxZ_V0(0)
,fVtxZ_SPD(0)
,fV0_SPD(0)
,fV0(0)
,fSPD(0)

,fNClusters(0)

//For the HFE package
,fCuts(0)
//,fCFM(0)
//,fPID(new AliHFEpid("hfePid"))
//,fPIDqa(0)

//For MC
,fMCstack(0)

,fMCtrack(0)
,fMCtrackMother(0)
,fMCtrackGMother(0)
,fMCtrackGGMother(0)
,fMCtrackGGGMother(0)
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

	//new histos
,fdEta_dPhi(0)

,fSparseElectron(0)
,fvalueElectron(0)
,fSparseMulti(0)
,fvalueMulti(0)
,fIspp2011(kFALSE)

	//MC efficiencies
,fPtMCparticleAllHfe1(0)
,fPtMCparticleRecoHfe1(0)
,fPtMCparticleAll_e_from_JPsi(0)
,fPtMCparticleAll_JPsi_pT(0)
,fPtMCparticleAll_trueJPsi_pT(0)
,fPtMCparticleReco_e_from_JPsi(0)
,fPtMCparticle_Total_e_from_JPsi(0)
,fPtMCparticle_Total_e_from_JPsi_sameMother(0)
,fPtMCparticle_TotalplusMass_e_from_JPsi(0)
,fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother(0)

,fPtMCparticle_TotalplusMass_JPsi_pT(0)
,fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother(0)


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
	delete []fvalueElectron;
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
            TProfile* hprof=new TProfile(*fMultEstimatorAvg[i]);
            hprof->SetName(Form("ProfileTrkVsZvtx%s\n",period[i].Data()));
            fListProfiles->Add(hprof);
        }
    }
    for(Int_t i=0; i<nProfilesV0; i++){
        if(fMultEstimatorV0[i]){
            TProfile2D* hprofV0=new TProfile2D(*fMultEstimatorV0[i]);
            hprofV0->SetName(Form("ProfileV0%s\n",period[i].Data()));
            fListProfiles->Add(hprofV0);
            
        }
    }
    
    PostData(2,fListProfiles);
    
    
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


//Store the number of events
	//Define the histo
	fNevent = new TH1F("fNevent","Number of Events",20,-0.5,19.5);
    fPDG_values = new TH1F("fPDG_values","PDG of generated particles",6000,-3000,3000);
   
	//And then, add to the output list
	fOutputList->Add(fNevent);
    fOutputList->Add(fPDG_values);
    
    fNevent_SPD_multi = new TH1F("fNevent_SPD_multi","Number of Events in SPD bins",10,-0.5,9.5);
    fOutputList->Add(fNevent_SPD_multi);
    
    fNevent_V0_multi = new TH1F("fNevent_V0_multi","Number of Events in V0 bins",10,-0.5,9.5);
    fOutputList->Add(fNevent_V0_multi);
	
	//General Histograms
	
	//Steps
	//Step 1: Before Track cuts
	//Step 2: Before PID
	//Step 3: After PID
	
	fEoverP_pt = new TH2F *[3];
	fTPC_p = new TH2F *[3];
	fTPCnsigma_p = new TH2F *[3];
	fTPCnsigma_EoverP = new TH2F *[3];
	fECluster= new TH1F *[3];
	
	fECluster_emcal= new TH1F *[3];
	fECluster_dcal= new TH1F *[3];
	
	fVtxZ= new  TH1F *[3];
	fNClusters= new TH1F *[3];

	fdEta_dPhi = new TH2F *[3];
	
	for(Int_t i = 0; i < 3; i++)
	{
	  fEoverP_pt[i] = new TH2F(Form("fEoverP_pt%d",i),";p_{t} (GeV/c);E / p ",600,0,30,500,0,2);
	  fTPC_p[i] = new TH2F(Form("fTPC_p%d",i),";p (GeV/c);TPC dE/dx (a. u.)",1000,0,20,1000,-20,200);
	  fTPCnsigma_p[i] = new TH2F(Form("fTPCnsigma_p%d",i),";p (GeV/c);TPC Electron N#sigma",1000,0,20,1000,-15,10);
	  fECluster[i]= new TH1F(Form("fECluster%d",i), ";ECluster",2000, 0,100);
		
		fECluster_emcal[i]= new TH1F(Form("fECluster_emcal%d",i), ";ECluster EMCal",2000, 0,100);
		fECluster_dcal[i]= new TH1F(Form("fECluster_dcal%d",i), ";ECluster DCal",2000, 0,100);
		
	  fVtxZ[i]= new  TH1F(Form("fVtxZ%d",i),"VtxZ",1000, -50,50);
	  fNClusters[i]= new TH1F(Form("fNClusters%d",i),"fNClusters0",100, 0,100);
		
	  fTPCnsigma_EoverP[i] = new TH2F(Form("fTPCnsigma_EoverP%d",i),";TPC Electron N#sigma; E/p",1600,-20,20,200,0,2);
		
			//new histos
	  fdEta_dPhi[i] = new TH2F(Form("fdEta_dPhi%d",i),"Distance of EMCAL cluster to its closest track ;#phi;z",100,-0.3,0.3,100,-0.3,0.3);

	  		
	  fOutputList->Add(fEoverP_pt[i]);
	  fOutputList->Add(fTPC_p[i]);
	  fOutputList->Add(fTPCnsigma_p[i]);
      fOutputList->Add(fTPCnsigma_EoverP[i]);
	  fOutputList->Add(fECluster[i]);
		
	  fOutputList->Add(fECluster_emcal[i]);
	  fOutputList->Add(fECluster_dcal[i]);
		
		
	  fOutputList->Add(fVtxZ[i]);
	  fOutputList->Add(fNClusters[i]);
		
      //new histos
	  fOutputList->Add(fdEta_dPhi[i]);

	}
    
    
    //=================================================================================================================================================================
    // Multiplicity histos
    
    fVtxZ_V0 = new TH2F("fVtxZ_V0","V0 multi vs. VtxZ ;VtxZ; V0 multiplicity",600,-30,30,1500,0,1500);
    fOutputList->Add(fVtxZ_V0);
    
    fVtxZ_SPD = new TH2F("fVtxZ_SPD","SPD multi vs. VtxZ ;VtxZ; SPD multiplicity",600,-30,30,500,0,500);
    fOutputList->Add(fVtxZ_SPD);
    
    fV0_SPD = new TH2F("fV0_SPD","SPD multi vs. V0 ;V0; SPD multiplicity",500,0,500,1500,0,1500);
    fOutputList->Add(fV0_SPD);
    
    fV0 = new TH1F("fV0","V0 ;V0 multiplicity",1500,0,1500);
    fOutputList->Add(fV0);
    
    fSPD = new TH1F("fSPD","SPD ;SPD multiplicity",500,0,500);
    fOutputList->Add(fSPD);
    
    
	//=================================================================================================================================================================
    
	fECluster_pure= new TH1F("fECluster_pure", ";ECluster pure",2000,0,100);
	fOutputList->Add(fECluster_pure);
	
    //emcal and dcal separated
	fECluster_pure_emcal= new TH1F("fECluster_pure_emcal", ";ECluster pure EMCal",2000,0,100);
    
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
    
    
	fOutputList->Add(fECluster_pure_emcal);
    
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
	
	fTracksQAPt=new TH1F *[11];
	for(Int_t i=0; i<11; i++){
		fTracksQAPt[i]= new TH1F(Form("fTracksQAPt%d", i), ";p_{T} (GeV/c); Counts ", 300, 0, 30);
		fOutputList->Add(fTracksQAPt[i]);
	}
	
	
	
	//JPsi analysis histograms
    /*
	fHist_InvMass_pt_ULS = new TH2F("fHist_InvMass_pt_ULS","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_ULS);
	fHist_InvMass_pt_LS = new TH2F("fHist_InvMass_pt_LS","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_LS);
     */
    
	
		//KFParticle
	fHist_InvMass_pt_ULS_KF = new TH2F("fHist_InvMass_pt_ULS_KF","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_ULS_KF);
	fHist_InvMass_pt_LS_KF = new TH2F("fHist_InvMass_pt_LS_KF","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_LS_KF);
    
    //multiplicity histos
    
    
    fHist_InvMass_pt_ULS_KF_SPDmulti_1 = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_1","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fHist_InvMass_pt_ULS_KF_SPDmulti_2 = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_2","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fHist_InvMass_pt_ULS_KF_SPDmulti_3 = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_3","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fHist_InvMass_pt_ULS_KF_SPDmulti_4 = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_4","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fHist_InvMass_pt_ULS_KF_SPDmulti_5 = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_5","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_1);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_2);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_3);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_4);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_5);
    
    fHist_InvMass_pt_ULS_KF_V0multi_1 = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_1","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fHist_InvMass_pt_ULS_KF_V0multi_2 = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_2","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fHist_InvMass_pt_ULS_KF_V0multi_3 = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_3","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fHist_InvMass_pt_ULS_KF_V0multi_4 = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_4","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fHist_InvMass_pt_ULS_KF_V0multi_5 = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_5","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);

    fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_1);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_2);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_3);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_4);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_5);
    
    //multiplicity histos with weight
    //multiplicity histos
    
    //KFParticle
    fHist_InvMass_pt_ULS_KF_weight = new TH2F("fHist_InvMass_pt_ULS_KF_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_weight);
    
    
    fHist_InvMass_pt_ULS_KF_SPDmulti_1_weight = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_1_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fHist_InvMass_pt_ULS_KF_SPDmulti_2_weight = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_2_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fHist_InvMass_pt_ULS_KF_SPDmulti_3_weight = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_3_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fHist_InvMass_pt_ULS_KF_SPDmulti_4_weight = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_4_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fHist_InvMass_pt_ULS_KF_SPDmulti_5_weight = new TH2F("fHist_InvMass_pt_ULS_KF_SPDmulti_5_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_1_weight);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_2_weight);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_3_weight);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_4_weight);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_SPDmulti_5_weight);
    
    fHist_InvMass_pt_ULS_KF_V0multi_1_weight = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_1_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fHist_InvMass_pt_ULS_KF_V0multi_2_weight = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_2_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fHist_InvMass_pt_ULS_KF_V0multi_3_weight = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_3_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fHist_InvMass_pt_ULS_KF_V0multi_4_weight = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_4_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    fHist_InvMass_pt_ULS_KF_V0multi_5_weight = new TH2F("fHist_InvMass_pt_ULS_KF_V0multi_5_weight","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
    
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_1_weight);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_2_weight);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_3_weight);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_4_weight);
    fOutputList->Add(fHist_InvMass_pt_ULS_KF_V0multi_5_weight);

    
	
	//=================================================================================================================================================================
	//MC generatros
	//BB
	fHist_InvMass_pt_ULS_KF_BB = new TH2F("fHist_InvMass_pt_ULS_KF_BB","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_ULS_KF_BB);
	fHist_InvMass_pt_LS_KF_BB = new TH2F("fHist_InvMass_pt_LS_KF_BB","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_LS_KF_BB);
	
	//CC
	fHist_InvMass_pt_ULS_KF_CC = new TH2F("fHist_InvMass_pt_ULS_KF_CC","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_ULS_KF_CC);
	fHist_InvMass_pt_LS_KF_CC = new TH2F("fHist_InvMass_pt_LS_KF_CC","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_LS_KF_CC);
	
	//B
	fHist_InvMass_pt_ULS_KF_B = new TH2F("fHist_InvMass_pt_ULS_KF_B","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_ULS_KF_B);
	fHist_InvMass_pt_LS_KF_B = new TH2F("fHist_InvMass_pt_LS_KF_B","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_LS_KF_B);
	
	//Jpsi
	fHist_InvMass_pt_ULS_KF_Jpsi = new TH2F("fHist_InvMass_pt_ULS_KF_Jpsi","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_ULS_KF_Jpsi);
	fHist_InvMass_pt_LS_KF_Jpsi = new TH2F("fHist_InvMass_pt_LS_KF_Jpsi","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_LS_KF_Jpsi);
	
	//BJpsi
	fHist_InvMass_pt_ULS_KF_BJpsi = new TH2F("fHist_InvMass_pt_ULS_KF_BJpsi","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_ULS_KF_BJpsi);
	fHist_InvMass_pt_LS_KF_BJpsi = new TH2F("fHist_InvMass_pt_LS_KF_BJpsi","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_LS_KF_BJpsi);
	
	//=================================================================================================================================================================

	
	fHist_InvMass_pt_ULS1 = new TH2F("fHist_InvMass_pt_ULS1","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_ULS1);
	fHist_InvMass_pt_LS1 = new TH2F("fHist_InvMass_pt_LS1","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_LS1);
	
	fHist_InvMass_pt_ULS2 = new TH2F("fHist_InvMass_pt_ULS2","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_ULS2);
	fHist_InvMass_pt_LS2 = new TH2F("fHist_InvMass_pt_LS2","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_LS2);
	
	fHist_InvMass_pt_ULSboth = new TH2F("fHist_InvMass_pt_ULSboth","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_ULSboth);
	fHist_InvMass_pt_LSboth = new TH2F("fHist_InvMass_pt_LSboth","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_LSboth);
	
	
	
	fHist_InvMass_pt_ULStpc = new TH2F("fHist_InvMass_pt_ULStpc","Invariant mass e^{-}e^{+} ;p_{T} (GeV/c); M_{e^{-}e^{+}}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_ULStpc);
	
	fHist_InvMass_pt_LStpc = new TH2F("fHist_InvMass_pt_LStpc","Invariant mass ee (like-sign) ;p_{T} (GeV/c); M_{ee}",300,0,30,500,0,5);
	fOutputList->Add(fHist_InvMass_pt_LStpc);
	
		//MC efficiencies
	fPtMCparticleRecoHfe1 = new TH1F("fPtMCparticleRecoHfe1",";p_{T} (GeV/c);Count",200,0,40);	
	fPtMCparticleAllHfe1 = new TH1F("fPtMCparticleAllHfe1",";p_{T} (GeV/c);Count",200,0,40);	
	fPtMCparticleAll_e_from_JPsi = new TH1F("fPtMCparticleAll_e_from_JPsi",";p_{T} (GeV/c);Count",200,0,40);
    fPtMCparticleAll_JPsi_pT = new TH1F("fPtMCparticleAll_JPsi_pT",";p_{T} (GeV/c);Count",200,0,40);
    fPtMCparticleAll_trueJPsi_pT = new TH1F("fPtMCparticleAll_trueJPsi_pT",";p_{T} (GeV/c);Count",200,0,40);
	fPtMCparticleReco_e_from_JPsi = new TH1F("fPtMCparticleReco_e_from_JPsi",";p_{T} (GeV/c);Count",200,0,40);
	
	fPtMCparticle_Total_e_from_JPsi = new TH1F("fPtMCparticle_Total_e_from_JPsi",";p_{T} (GeV/c);Count",200,0,40);
    fPtMCparticle_Total_e_from_JPsi_sameMother = new TH1F("fPtMCparticle_Total_e_from_JPsi_sameMother",";p_{T} (GeV/c);Count",200,0,40);
	fPtMCparticle_TotalplusMass_e_from_JPsi = new TH1F("fPtMCparticle_TotalplusMass_e_from_JPsi",";p_{T} (GeV/c);Count",200,0,40);
    fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother = new TH1F("fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother",";p_{T} (GeV/c);Count",200,0,40);
   
    fPtMCparticle_TotalplusMass_JPsi_pT = new TH1F("fPtMCparticle_TotalplusMass_JPsi_pT",";p_{T} (GeV/c);Count",200,0,40);
    fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother = new TH1F("fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother",";p_{T} (GeV/c);Count",200,0,40);
    
	
	fOutputList->Add(fPtMCparticleRecoHfe1);
	fOutputList->Add(fPtMCparticleAllHfe1);
	fOutputList->Add(fPtMCparticleAll_e_from_JPsi);
    fOutputList->Add(fPtMCparticleAll_JPsi_pT);
    fOutputList->Add(fPtMCparticleAll_trueJPsi_pT);
	fOutputList->Add(fPtMCparticleReco_e_from_JPsi);
	fOutputList->Add(fPtMCparticle_Total_e_from_JPsi);
    fOutputList->Add(fPtMCparticle_Total_e_from_JPsi_sameMother);
	fOutputList->Add(fPtMCparticle_TotalplusMass_e_from_JPsi);
    fOutputList->Add(fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother);
    
    fOutputList->Add(fPtMCparticle_TotalplusMass_JPsi_pT);
    fOutputList->Add(fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother);
	
	
    //TOF
	fTOF_p = new TH2F("fTOF_p",";p (GeV/c);TOF signal (a. u.)",1000,0,20,1000,-20,20);
	fTOFnsigma_p = new TH2F("fTOFnsigma_p",";p (GeV/c);TOF Electron N#sigma",1000,0,20,1000,-15,10);
	fOutputList->Add(fTOF_p);
	fOutputList->Add(fTOFnsigma_p);
	
	fvalueElectron = new Double_t[9];
    fvalueMulti = new Double_t[6];
	
    //electron Sparse
	Int_t bins[9]={600,300,400,200,200, 200, 120, 6,6}; // pt, TPCnsig, E/p, M20, M02, E, V0, SPD
	Double_t xmin[9]={0,-15,0,0,0,0,0,0,0};
	Double_t xmax[9]={30,15,2,2,2,40,6, 450, 90};
	fSparseElectron = new THnSparseD ("Electron","Electron",9,bins,xmin,xmax);
	fOutputList->Add(fSparseElectron);
    
    
    //multi Sparse
    Int_t binsm[6]    =          {200,1000,200,1000,200,1000};
    Double_t xminm[6]    =    {-10,0,0,0,0,0};
    Double_t xmaxm[6]    =    { 10,1000,200,1000,200,1000};
    fSparseMulti         = new THnSparseD ("Multiplicity","Multiplicity;zvtx;V0M;SPDTracklets;Corrected_V0M;Corrected_SPDTracklets;ncharge;",6,binsm,xminm,xmaxm);
    fOutputList->Add(fSparseMulti);
			
//______________________________________________________________________
	
	PostData(1, fOutputList);
    PostData(2,fListProfiles);
	
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
	
	fNevent->Fill(10);
	if(fIsAOD)
	{
		const AliAODVertex* trkVtx = fAOD->GetPrimaryVertex();
		    Float_t zvtx = trkVtx->GetZ();
			fZvtx = zvtx;
			const AliAODVertex* spdVtx = fAOD->GetPrimaryVertexSPD();
			//Any vertex
            if((!trkVtx || trkVtx->GetNContributors()<=0) && (spdVtx->GetNContributors()<=0)) return;
            fNevent->Fill(9);
		    if(TMath::Abs(zvtx) > fVertexCut) return;
	}
	else
	{
        const AliESDVertex *trkVtx = fESD->GetPrimaryVertex();
		Float_t zvtx = trkVtx->GetZ();
		fZvtx = zvtx;
		if(TMath::Abs(zvtx) > fVertexCut) return;
	}

	fNevent->Fill(8);
    
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
	
    fNevent->Fill(7);
    
    
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
    
    // Multiplicity histos without correction
    fVtxZ_V0->Fill(fZvtx, fV0Mult);
    fVtxZ_SPD->Fill(fZvtx, fSPDMult);
    fV0_SPD->Fill(fV0Mult,fSPDMult);
    
    fV0->Fill(fV0Mult);
    fSPD->Fill(fSPDMult);
    
    
    //=======
    //correction for multiplicity
    TProfile* estimatorAvg = GetEstimatorHistogram(fAOD);
    //if(!isMC){
        if(estimatorAvg){
            //printf("Estimator SPD exists!\n");
            //correctednAcc=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,nAcc,Zvertex1,fRefMult));
            fSPDMult_corr = AliAnalysisTask_JPsi_EMCal::GetTrackletsMeanCorrection(estimatorAvg,nAcc,fZvtx,fRefMult);
        
            //countCorr=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,countMult,Zvertex1,fRefMult));
        }
        else{
            fSPDMult_corr=fSPDMult;
            //if we do not load any correction, variable is the same as uncorrected
        }
    //}

 
   // printf("SPD =%f, SPD_corrected =%f\n", SPDMult, fSPDMult_corr);
    

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

    
    
    //printf("V0 =%d, V0_corrected =%f,  V0_corrected2 =%f\n", V0Mult, fV0Mult_corr, fV0Mult_corr2);
   
    
    if(fAOD->IsPileupFromSPDInMultBins()){
        //printf("This event is pileUp from AOD\n");
        fNevent->Fill(6);
        return;
    }
 
//______________________________________________________________________	
	
//Only events with at least 2 tracks are accepted
	Int_t fNOtrks =  fVevent->GetNumberOfTracks();
	if(fNOtrks<2) return;
	fNevent->Fill(5);
	
	
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
		fNevent->Fill(4);
		IsEventEMCALL1=kTRUE;
	}
	if(firedTrigger.Contains(TriggerEG2)){
		fNevent->Fill(3);
		IsEventEMCALL1=kTRUE;
	}
	
		//if the flag is for a given threshold and it was not fired, return.
		//EMCal trigger word
	if(fEMCEG1){
		if(!firedTrigger.Contains(TriggerEG1))return;
		if(firedTrigger.Contains(TriggerEG2)){
			fNevent->Fill(2);
			
		}
		
	}
	
	if(fEMCEG2){
		if(!firedTrigger.Contains(TriggerEG2))return;
		if(firedTrigger.Contains(TriggerEG1)){
			fNevent->Fill(1);
		}
		
	}
	//=====================================================
	//DCal trigger word
	if(fEMCDG1){
		if(!firedTrigger.Contains(TriggerDG1))return;
		if(firedTrigger.Contains(TriggerDG2)){
				//fNevent->Fill(2);
			
		}
		
	}
	
	if(fEMCDG2){
		if(!firedTrigger.Contains(TriggerDG2))return;
		if(firedTrigger.Contains(TriggerDG1)){
				//fNevent->Fill(1);
		}
		
	}

	
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
						
		
			
			
			if(IsMB_gen)fNevent->Fill(13);
			if(IsPythiaCC_gen)fNevent->Fill(14);
			if(IsPythiaBB_gen)fNevent->Fill(15);
			if(IsPythiaB_gen)fNevent->Fill(16);
			if(IsJpsi2ee_gen)fNevent->Fill(17);
			if(IsB2JPsi2ee_gen)fNevent->Fill(18);
			 
			
			
			
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
				
				
				if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax && fMCparticle->Charge()!=0)
				{
					
                    fPDG_values->Fill(fMCparticle->GetPdgCode());
                    //Take all J/psi generated
                    //if(fMCparticle->IsPhysicalPrimary()){
                        
                        if(TMath::Abs(fMCparticle->GetPdgCode())==443)
                        {
                            fPtMCparticleAll_trueJPsi_pT->Fill(fMCparticle->Pt());
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
    
    fvalueMulti[0] = fZvtx;
    fvalueMulti[1] = fV0Mult_corr;
    fvalueMulti[2] = fSPDMult;
    fvalueMulti[3] = fV0Mult_corr2;
    fvalueMulti[4] = fSPDMult_corr;
    fvalueMulti[5] = Nch;
    
    
    fSparseMulti->Fill(fvalueMulti);    // multiplicity from tracklets
    
    
    
	
//______________________________________________________________________
	
    
    // Number of events
    fNevent->Fill(0);
    
    
    if(fSPDMult_corr>0 && fSPDMult_corr < 15)   fNevent_SPD_multi->Fill(1);
    if(fSPDMult_corr>=15 && fSPDMult_corr < 25) fNevent_SPD_multi->Fill(2);
    if(fSPDMult_corr>=25 && fSPDMult_corr < 35) fNevent_SPD_multi->Fill(3);
    if(fSPDMult_corr>=35 && fSPDMult_corr < 45) fNevent_SPD_multi->Fill(4);
    if(fSPDMult_corr>=45 && fSPDMult_corr < 88) fNevent_SPD_multi->Fill(5);
    
    if(fV0Mult_corr2>0 && fV0Mult_corr2<75)     fNevent_V0_multi->Fill(1);
    if(fV0Mult_corr2>=75 && fV0Mult_corr2<150)  fNevent_V0_multi->Fill(2);
    if(fV0Mult_corr2>=150 && fV0Mult_corr2<225) fNevent_V0_multi->Fill(3);
    if(fV0Mult_corr2>=225 && fV0Mult_corr2<300) fNevent_V0_multi->Fill(4);
    if(fV0Mult_corr2>=300 && fV0Mult_corr2<680) fNevent_V0_multi->Fill(5);
    
    

	
	Int_t ClsNo = -999;
	if(!fIsAOD) ClsNo = fESD->GetNumberOfCaloClusters(); 
	else ClsNo = fAOD->GetNumberOfCaloClusters(); 
	
   
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
 
	
	AliVCluster *clust = 0x0;
	
	if(!fUseTender){
		if(fIsAOD){
			for (Int_t i=0; i< ClsNo; i++ ){
				clust = (AliVCluster*) fAOD->GetCaloCluster(i);
				
				if(clust && clust->IsEMCAL())
				{
					fECluster_pure->Fill(clust->E());
					
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
                        
                        if(fSPDMult_corr>0 && fSPDMult_corr < 15)   fECluster_pure_emcal_SPD1->Fill(clust->E());
                        if(fSPDMult_corr>=15 && fSPDMult_corr < 25) fECluster_pure_emcal_SPD2->Fill(clust->E());
                        if(fSPDMult_corr>=25 && fSPDMult_corr < 35) fECluster_pure_emcal_SPD3->Fill(clust->E());
                        if(fSPDMult_corr>=35 && fSPDMult_corr < 45) fECluster_pure_emcal_SPD4->Fill(clust->E());
                        if(fSPDMult_corr>=45 && fSPDMult_corr < 88) fECluster_pure_emcal_SPD5->Fill(clust->E());
                        
                        if(fV0Mult_corr2>0 && fV0Mult_corr2<75)      fECluster_pure_emcal_V01->Fill(clust->E());
                        if(fV0Mult_corr2>=75 && fV0Mult_corr2<150)   fECluster_pure_emcal_V02->Fill(clust->E());
                        if(fV0Mult_corr2>=150 && fV0Mult_corr2<225)  fECluster_pure_emcal_V03->Fill(clust->E());
                        if(fV0Mult_corr2>=225 && fV0Mult_corr2<300)  fECluster_pure_emcal_V04->Fill(clust->E());
                        if(fV0Mult_corr2>=300 && fV0Mult_corr2<680)  fECluster_pure_emcal_V05->Fill(clust->E());
                        
 
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
		for (Int_t i=0; i< ClsNo; i++ ){
			
			clust = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(i));
			if (!clust) {
			continue;
			}
			if(clust && clust->IsEMCAL())
			{
				fECluster_pure->Fill(clust->E());
				
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
                    
                    if(fSPDMult_corr>0 && fSPDMult_corr < 15)   fECluster_pure_emcal_SPD1->Fill(clust->E());
                    if(fSPDMult_corr>=15 && fSPDMult_corr < 25) fECluster_pure_emcal_SPD2->Fill(clust->E());
                    if(fSPDMult_corr>=25 && fSPDMult_corr < 35) fECluster_pure_emcal_SPD3->Fill(clust->E());
                    if(fSPDMult_corr>=35 && fSPDMult_corr < 45) fECluster_pure_emcal_SPD4->Fill(clust->E());
                    if(fSPDMult_corr>=45 && fSPDMult_corr < 88) fECluster_pure_emcal_SPD5->Fill(clust->E());
                    
                    if(fV0Mult_corr2>0 && fV0Mult_corr2<75)      fECluster_pure_emcal_V01->Fill(clust->E());
                    if(fV0Mult_corr2>=75 && fV0Mult_corr2<150)   fECluster_pure_emcal_V02->Fill(clust->E());
                    if(fV0Mult_corr2>=150 && fV0Mult_corr2<225)  fECluster_pure_emcal_V03->Fill(clust->E());
                    if(fV0Mult_corr2>=225 && fV0Mult_corr2<300)  fECluster_pure_emcal_V04->Fill(clust->E());
                    if(fV0Mult_corr2>=300 && fV0Mult_corr2<680)  fECluster_pure_emcal_V05->Fill(clust->E());
                    
				}
				
				//dcal
				if(cphi>=3.9){
					fEtaPhi_dcal->Fill(cphi,ceta);
					fECluster_pure_dcal->Fill(clust->E());
				}

			}
			
		}
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
			printf("ERROR: Could not receive track %d\n", iTracks);
			continue;
		}
     
		AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
		AliESDtrack *etrack = dynamic_cast<AliESDtrack*>(Vtrack);
		AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);
		
		
		Double_t eta =0;
		eta = track->Eta();
		if(eta > fEtaCutMax || eta < fEtaCutMin) continue;
		
		
		if(fIspp2011){
			Double_t phi=0;
			phi = track->Phi();
			if(phi<0 || phi>4 ) continue;
		}
		
		
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
		Double_t fTOFnsigma = -999;
		Double_t fTOFsignal = -999;
		
		
		Float_t pos0[3]={0,0,0};
		Float_t pos1[3]={0,0,0};
		Float_t pos2[3]={0,0,0};

		///_____________________________________________________________________________
		///Fill QA plots without track selection
		fPt = track->Pt();
		fP = TMath::Sqrt((track->Pt())*(track->Pt()) + (track->Pz())*(track->Pz()));
		
		fTPCsignal = track->GetTPCsignal();
		fTPCnSigma = fPidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
		fTPCnSigma_pion = fPidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
		fTPCnSigma_proton = fPidResponse->NumberOfSigmasTPC(track, AliPID::kProton);
		fTPCnSigma_kaon = fPidResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
		
		//TOF
		fTOFsignal = track->GetTOFsignal();
		fTOFnsigma = fPidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
		
		fTOF_p->Fill(fP,fTOFsignal);
        fTOFnsigma_p->Fill(fP,fTOFnsigma);
		
        
        fTPC_p[0]->Fill(fP,fTPCsignal);
        fTPCnsigma_p[0]->Fill(fP,fTPCnSigma);
        
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
					fNClusters[0]->Fill(ClsNo);
					
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
        
        fVtxZ[0]->Fill(fZvtx);
		fTracksPt[0]->Fill(fPt);
		fTracksQAPt[0]->Fill(fPt);
		
//=======================================================================
// Track Selection Cuts are applied here
//=======================================================================
        
        

		if(fAOD){
			
			//TPCncls
			if(atrack->GetTPCNcls() < fTPCncls) continue;
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
           // printf("TPCchi2/Ncls = %f, cut =%f\n",((track->GetTPCchi2())/(atrack->GetTPCNcls())), fTPCchi2);
           if(((track->GetTPCchi2())/(atrack->GetTPCNcls())) > fTPCchi2){
                continue;
            }
            fTracksQAPt[8]->Fill(fPt);
							
            if(fAODGlobalTracks){
                if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts
            }
            fTracksQAPt[9]->Fill(fPt);

		}
		
        //if(atrack->GetTPCsignalN() < 80)
        //if(atrack->GetTPCNclsF() < 0.6)
	
		        
//=======================================================================
// QA plots after track selection
//=======================================================================

		fTracksPt[1]->Fill(fPt);
		fTPC_p[1]->Fill(fP,fTPCsignal);
		fTPCnsigma_p[1]->Fill(fP,fTPCnSigma);
		
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
		 
		   if(fMCparticle->IsPhysicalPrimary()) 
		   {
		 
		 
		     Bool_t MotherFound = FindMother(TMath::Abs(track->GetLabel()));
		 
		     //For JPsi analysis
               if(fMCparticle->GetMother()>0){
		 
                   fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                   if(fMCparticleMother->GetMother()>0)fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
		 
                   if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==443)){
			    	 fPtMCparticleReco_e_from_JPsi->Fill(track->Pt()); //reconstructed pT
                   }
               }
		 
               if(MotherFound){
                    if(fIsHFE1){
                        fPtMCparticleRecoHfe1->Fill(track->Pt());//numerator tracking  reconstructed pT (unfolding)
														 //fpt_reco_pt_MC_den->Fill(track->Pt(),fMCparticle->Pt());
                    }
                }
               
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
                      fIsTrack1Emcal=kTRUE;
				  }
				  
                //dcal
				  if(cphi > 4.53 && cphi < 5.708){
					  fECluster_dcal[1]->Fill(fClus->E());
                      fIsTrack1Dcal=kTRUE;
				  }
			  
			    //}
                //EID THnsparse
				
				fvalueElectron[0] = track->Pt();
				fvalueElectron[1] = fTPCnSigma;
				fvalueElectron[2] = fClus->E() / fP;
				fvalueElectron[3] = fClus->GetM20();
				fvalueElectron[4] = fClus->GetM02();
				fvalueElectron[5] = fClus->E(); // to check rejection factor for electrons
                fvalueElectron[6] = cphi; //to separate emcal and dcal
                fvalueElectron[7] = fV0Mult;//to check RF in bins of multiplicity (bins not exactly same as in the analysis...)
                fvalueElectron[8] = fSPDMult;//to check RF in bins of multiplicity (bins not exactly same as in the analysis...)
				
				fSparseElectron->Fill(fvalueElectron);
				
			}
		}
		
		fVtxZ[1]->Fill(fZvtx);
		

			
//=======================================================================
// Here the PID cut defined in the file "ConfigEMCalHFEpA.C" is applied
//=======================================================================
		
		if(fTPCnSigma < fTPCnsigmaCutMin || fTPCnSigma > fTPCnsigmaCutMax) continue;
        
        //printf("Main leg: Track1 on Electron band with fPt=%f\n", fPt);

	    fTracksQAPt[10]->Fill(fPt);
	
		
        
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
					printf("ERROR: Could not receive track %d\n", lTracks);
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
                if(atrack2->GetTPCNcls() < fTPCncls) continue;
            
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
               

                if(fAODGlobalTracks){
                    if(!atrack2->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts
                }
              
                
            }
			
			fTracksPt[4]->Fill(fPt2);

			 
			 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				Double_t dEdx3 =-999, fTPCnSigma2=-999;
				dEdx3 = track2->GetTPCsignal();
				fTPCnSigma2 = fPidResponse->NumberOfSigmasTPC(track2, AliPID::kElectron);
				if(fTPCnSigma2 > fTPCnsigmaCutMin && fTPCnSigma2 < fTPCnsigmaCutMax){
                    
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
                    
                    //if use pT of J/psi
                    Double_t a = 7.21;
                    Double_t b = 9.1;
                    Double_t c = 0.0335;
                    Double_t weight = 1./(TMath::Erf((pt_kf-a)/b)*c);  //   weight = 1/eff
                    
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
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
                        //to use tender
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
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
                        
						if(fClus2->IsEMCAL())
						{
							
														
								//if(TMath::Abs(fClus2->GetTrackDx())<=0.05 && TMath::Abs(fClus2->GetTrackDz())<=0.05)
								//{
								
								
							 	fECluster[2]->Fill(fClus2->E());
							    fTracksPt[7]->Fill(fPt2);
							
					//======================================// for Eta Phi distribution
							fClus2->GetPosition(pos2);
							TVector3 vpos2(pos2[0],pos2[1],pos2[2]);
							Double_t cphi = vpos2.Phi();
							Double_t ceta = vpos2.Eta();
							
				
							
								///from emcal QA task
							if(cphi < 0) cphi = cphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.
																	  // if(cphi > 1.39 && cphi < 3.265) ; //EMCAL : 80 < phi < 187
																	  // if(cphi > 4.53 && cphi < 5.708) ; //DCAL  : 260 < phi < 327
							
							
                           
                            
                            //emcal
							if(cphi > 1.39 && cphi < 3.265){
								fECluster_emcal[2]->Fill(fClus2->E());
                                fIsTrack2Emcal=kTRUE;
							}
							
								//dcal
							if(cphi > 4.53 && cphi < 5.708){
								fECluster_dcal[2]->Fill(fClus2->E());
                                fIsTrack2Dcal=kTRUE;

							}
							
					//======================================
								//}
						}
					}
					//==================================
					//Filling the invariant mass spectrum
					
					
					if(fIsTrack1Emcal && (!fIsTrack2Emcal)){
						//printf("Track1 is on EMCal and track2 is not \n");
                       // if(fEMCEG1 || fEMCEG2){
                            
                        //printf("The cuts are: %f , E/p < %f and E > %f\n",fEoverPCutMin, fEoverPCutMax, fEnergyCut);
                            
						  if((fClus->E() / fP) >=fEoverPCutMin && (fClus->E() / fP) <=fEoverPCutMax && (fClus->E()) >= fEnergyCut){
							
						//	printf("Track1 PASSED the cuts \n");
                     //   printf("weigh=% f \n", weight);
                       // printf("Track1 has pt=%f \n", fPt);
                        
							
							//sum of all possibilities on emcal
							//if(charge1*charge2 <0) fHist_InvMass_pt_ULS->Fill(pt3,invmass3);
							//if(charge1*charge2 >0) fHist_InvMass_pt_LS->Fill(pt3,invmass3);
							
                            //KFParticle
                              if(charge1*charge2 <0){
                                 // printf("Inside first if \n");
                                  fHist_InvMass_pt_ULS_KF->Fill(pt_kf,imass);//multi integrated
                                //  printf("weigh=% f \n", weight);
                                 // printf("passed first histo \n");
                                  fHist_InvMass_pt_ULS_KF_weight->Fill(pt_kf,imass,weight);//multi integrated with weight
                                  //printf("passed second histo \n");
                              }
							if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF->Fill(pt_kf,imass);
                              
                            //multiplicity bins histos (only ULS for SPDmulti and V0multi)
                              if(charge1*charge2 <0){
                                  
                                //  printf("track1: fSPDMult = %f, fVOMult = %f", fSPDMult,fV0Mult );
                                  if(fSPDMult_corr>0 && fSPDMult_corr < 15)   fHist_InvMass_pt_ULS_KF_SPDmulti_1->Fill(pt_kf,imass);
                                  if(fSPDMult_corr>=15 && fSPDMult_corr < 25) fHist_InvMass_pt_ULS_KF_SPDmulti_2->Fill(pt_kf,imass);
                                  if(fSPDMult_corr>=25 && fSPDMult_corr < 35) fHist_InvMass_pt_ULS_KF_SPDmulti_3->Fill(pt_kf,imass);
                                  if(fSPDMult_corr>=35 && fSPDMult_corr < 45) fHist_InvMass_pt_ULS_KF_SPDmulti_4->Fill(pt_kf,imass);
                                  if(fSPDMult_corr>=45 && fSPDMult_corr < 88) fHist_InvMass_pt_ULS_KF_SPDmulti_5->Fill(pt_kf,imass);
                                  
                                  if(fV0Mult_corr2>0 && fV0Mult_corr2<75)     fHist_InvMass_pt_ULS_KF_V0multi_1->Fill(pt_kf,imass);
                                  if(fV0Mult_corr2>=75 && fV0Mult_corr2<150)  fHist_InvMass_pt_ULS_KF_V0multi_2->Fill(pt_kf,imass);
                                  if(fV0Mult_corr2>=150 && fV0Mult_corr2<225) fHist_InvMass_pt_ULS_KF_V0multi_3->Fill(pt_kf,imass);
                                  if(fV0Mult_corr2>=225 && fV0Mult_corr2<300) fHist_InvMass_pt_ULS_KF_V0multi_4->Fill(pt_kf,imass);
                                  if(fV0Mult_corr2>=300 && fV0Mult_corr2<680) fHist_InvMass_pt_ULS_KF_V0multi_5->Fill(pt_kf,imass);
                                  
                                  //with weight
                                  if(fSPDMult_corr>0 && fSPDMult_corr < 15)   fHist_InvMass_pt_ULS_KF_SPDmulti_1_weight->Fill(pt_kf,imass,weight);
                                  if(fSPDMult_corr>=15 && fSPDMult_corr < 25) fHist_InvMass_pt_ULS_KF_SPDmulti_2_weight->Fill(pt_kf,imass,weight);
                                  if(fSPDMult_corr>=25 && fSPDMult_corr < 35) fHist_InvMass_pt_ULS_KF_SPDmulti_3_weight->Fill(pt_kf,imass,weight);
                                  if(fSPDMult_corr>=35 && fSPDMult_corr < 45) fHist_InvMass_pt_ULS_KF_SPDmulti_4_weight->Fill(pt_kf,imass,weight);
                                  if(fSPDMult_corr>=45 && fSPDMult_corr < 88) fHist_InvMass_pt_ULS_KF_SPDmulti_5_weight->Fill(pt_kf,imass,weight);
                                  
                                  if(fV0Mult_corr2>0 && fV0Mult_corr2<75)     fHist_InvMass_pt_ULS_KF_V0multi_1_weight->Fill(pt_kf,imass,weight);
                                  if(fV0Mult_corr2>=75 && fV0Mult_corr2<150)  fHist_InvMass_pt_ULS_KF_V0multi_2_weight->Fill(pt_kf,imass,weight);
                                  if(fV0Mult_corr2>=150 && fV0Mult_corr2<225) fHist_InvMass_pt_ULS_KF_V0multi_3_weight->Fill(pt_kf,imass,weight);
                                  if(fV0Mult_corr2>=225 && fV0Mult_corr2<300) fHist_InvMass_pt_ULS_KF_V0multi_4_weight->Fill(pt_kf,imass,weight);
                                  if(fV0Mult_corr2>=300 && fV0Mult_corr2<680) fHist_InvMass_pt_ULS_KF_V0multi_5_weight->Fill(pt_kf,imass,weight);
                                  
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
								
								
								
								//printf("It is on MC if for case: Track1 PASSED the cuts \n");
								
							
								
								if(fIsAOD)
								{
									fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
									fMCparticle2 = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track2->GetLabel()));
									
									//=================================================================
									//checking the generator of each particle of event
									//if(fMCparticle->GetGeneratorIndex()==0)printf("This is a particle from MB event!\n");
									
									//if(fMCparticle->GetGeneratorIndex()==1)printf("This is a particle from OTHER event!\n");
									//=================================================================
									
									
									
									
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
														
														//printf("Label leg1 %d, leg2 %d\n", fMCparticle->GetLabel(),fMCparticle2->GetLabel());
                                                        //printf("Label mother leg1 %d, leg2 %d\n", fMCparticleMother->GetLabel(),fMCparticleMother2->GetLabel());
														
														
														fPtMCparticle_Total_e_from_JPsi->Fill(track->Pt()); //reconstructed pT
                                                        //checking if they are from same mother
                                                        if((fMCparticleMother->GetLabel())==(fMCparticleMother2->GetLabel())){
                                                            //printf("electrons from same mother\n\n");
                                                            fPtMCparticle_Total_e_from_JPsi_sameMother->Fill(track->Pt()); //reconstructed pT
                                                        }

														if(invmass3>=fMassCutMin && invmass3<=fMassCutMax){
															fPtMCparticle_TotalplusMass_e_from_JPsi->Fill(track->Pt()); //reconstructed pT from electron
                                                            
                                                            fPtMCparticle_TotalplusMass_JPsi_pT->Fill(pt_kf);//spectrum of reconstructed J/Psi
                                                            
                                                            if((fMCparticleMother->GetLabel())==(fMCparticleMother2->GetLabel())){
                                                               // printf("electrons from same mother + inv mass cut\n\n");
                                                                fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother->Fill(track->Pt()); //reconstructed pT
                                                                fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother->Fill(pt_kf);//spectrum of reconstructed J/Psi
                                                                
                                                            }
                                                            
														}
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
                      //}//close trigger for emcal
						
					}
					
					if((!fIsTrack1Emcal) && fIsTrack2Emcal){
                       // printf("Track2 is on EMCal and track1 is not \n");
                        
                       // if(fEMCEG1 || fEMCEG2){
                            
                         // printf("The cuts are: %f , E/p < %f and E > %f\n",fEoverPCutMin, fEoverPCutMax, fEnergyCut);
						
						 if((fClus2->E()/fP2) >=fEoverPCutMin && (fClus2->E()/fP2) <=fEoverPCutMax && (fClus2->E()) >= fEnergyCut){
							
                           // printf("Track2 PASSED the cuts \n");
                           // printf("Track2 has pt=%f \n", fPt2);
							
							//if(charge1*charge2 <0) fHist_InvMass_pt_ULS->Fill(pt3,invmass3);
							//if(charge1*charge2 >0) fHist_InvMass_pt_LS->Fill(pt3,invmass3);
							
								//KFParticle
                             if(charge1*charge2 <0){
                                 fHist_InvMass_pt_ULS_KF->Fill(pt_kf,imass);//multi integrated
                                 fHist_InvMass_pt_ULS_KF_weight->Fill(pt_kf,imass,weight);//multi integrated with weight
                             }
							if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF->Fill(pt_kf,imass);
                             
                             //multiplicity bins histos (only ULS for SPDmulti and V0multi)
                             if(charge1*charge2 <0){
                                 //printf("track2: fSPDMult = %f, fVOMult = %f", fSPDMult,fV0Mult );
                                 if(fSPDMult_corr>0 && fSPDMult_corr < 15)   fHist_InvMass_pt_ULS_KF_SPDmulti_1->Fill(pt_kf,imass);
                                 if(fSPDMult_corr>=15 && fSPDMult_corr < 25) fHist_InvMass_pt_ULS_KF_SPDmulti_2->Fill(pt_kf,imass);
                                 if(fSPDMult_corr>=25 && fSPDMult_corr < 35) fHist_InvMass_pt_ULS_KF_SPDmulti_3->Fill(pt_kf,imass);
                                 if(fSPDMult_corr>=35 && fSPDMult_corr < 45) fHist_InvMass_pt_ULS_KF_SPDmulti_4->Fill(pt_kf,imass);
                                 if(fSPDMult_corr>=45 && fSPDMult_corr < 88) fHist_InvMass_pt_ULS_KF_SPDmulti_5->Fill(pt_kf,imass);
                                 
                                 if(fV0Mult_corr2>0 && fV0Mult_corr2<75)     fHist_InvMass_pt_ULS_KF_V0multi_1->Fill(pt_kf,imass);
                                 if(fV0Mult_corr2>=75 && fV0Mult_corr2<150)  fHist_InvMass_pt_ULS_KF_V0multi_2->Fill(pt_kf,imass);
                                 if(fV0Mult_corr2>=150 && fV0Mult_corr2<225) fHist_InvMass_pt_ULS_KF_V0multi_3->Fill(pt_kf,imass);
                                 if(fV0Mult_corr2>=225 && fV0Mult_corr2<300) fHist_InvMass_pt_ULS_KF_V0multi_4->Fill(pt_kf,imass);
                                 if(fV0Mult_corr2>=300 && fV0Mult_corr2<680) fHist_InvMass_pt_ULS_KF_V0multi_5->Fill(pt_kf,imass);
                                 
                                 //with weight
                                 if(fSPDMult_corr>0 && fSPDMult_corr < 15)   fHist_InvMass_pt_ULS_KF_SPDmulti_1_weight->Fill(pt_kf,imass,weight);
                                 if(fSPDMult_corr>=15 && fSPDMult_corr < 25) fHist_InvMass_pt_ULS_KF_SPDmulti_2_weight->Fill(pt_kf,imass,weight);
                                 if(fSPDMult_corr>=25 && fSPDMult_corr < 35) fHist_InvMass_pt_ULS_KF_SPDmulti_3_weight->Fill(pt_kf,imass,weight);
                                 if(fSPDMult_corr>=35 && fSPDMult_corr < 45) fHist_InvMass_pt_ULS_KF_SPDmulti_4_weight->Fill(pt_kf,imass,weight);
                                 if(fSPDMult_corr>=45 && fSPDMult_corr < 88) fHist_InvMass_pt_ULS_KF_SPDmulti_5_weight->Fill(pt_kf,imass,weight);
                                 
                                 if(fV0Mult_corr2>0 && fV0Mult_corr2<75)     fHist_InvMass_pt_ULS_KF_V0multi_1_weight->Fill(pt_kf,imass,weight);
                                 if(fV0Mult_corr2>=75 && fV0Mult_corr2<150)  fHist_InvMass_pt_ULS_KF_V0multi_2_weight->Fill(pt_kf,imass,weight);
                                 if(fV0Mult_corr2>=150 && fV0Mult_corr2<225) fHist_InvMass_pt_ULS_KF_V0multi_3_weight->Fill(pt_kf,imass,weight);
                                 if(fV0Mult_corr2>=225 && fV0Mult_corr2<300) fHist_InvMass_pt_ULS_KF_V0multi_4_weight->Fill(pt_kf,imass,weight);
                                 if(fV0Mult_corr2>=300 && fV0Mult_corr2<680) fHist_InvMass_pt_ULS_KF_V0multi_5_weight->Fill(pt_kf,imass,weight);
                                 
                                 
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
															
																//printf("Label mother leg 1 %d, leg %d\n", fMCparticle->GetLabel(),fMCparticle2->GetLabel());
															
																
															fPtMCparticle_Total_e_from_JPsi->Fill(track->Pt()); //reconstructed pT
                                                            //checking if they are from same mother
                                                            if((fMCparticleMother->GetLabel())==(fMCparticleMother2->GetLabel())){
                                                                //printf("electrons from same mother\n\n");
                                                                fPtMCparticle_Total_e_from_JPsi_sameMother->Fill(track->Pt()); //reconstructed pT
                                                            }
															
															if(invmass3>=fMassCutMin && invmass3<=fMassCutMax){
                                                                fPtMCparticle_TotalplusMass_e_from_JPsi->Fill(track->Pt()); //reconstructed pT from electron
                                                                
                                                                fPtMCparticle_TotalplusMass_JPsi_pT->Fill(pt_kf);//spectrum of reconstructed J/Psi
                                                                
                                                                if((fMCparticleMother->GetLabel())==(fMCparticleMother2->GetLabel())){
                                                                    //printf("electrons from same mother + inv mass cut\n\n");
                                                                    fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother->Fill(track->Pt()); //reconstructed pT
                                                                    fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother->Fill(pt_kf);//spectrum of reconstructed J/Psi
                                                                    
                                                                }
                                                                
															}
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
                                fHist_InvMass_pt_ULS_KF_weight->Fill(pt_kf,imass,weight);//multi integrated with weight
                            }
							if(charge1*charge2 >0) fHist_InvMass_pt_LS_KF->Fill(pt_kf,imass);
                            
                            //multiplicity bins histos (only ULS for SPDmulti and V0multi)
                            if(charge1*charge2 <0){
                               // printf("both tracks: fSPDMult = %f, fVOMult = %f", fSPDMult,fV0Mult );
                                if(fSPDMult_corr>0 && fSPDMult_corr < 15)   fHist_InvMass_pt_ULS_KF_SPDmulti_1->Fill(pt_kf,imass);
                                if(fSPDMult_corr>=15 && fSPDMult_corr < 25) fHist_InvMass_pt_ULS_KF_SPDmulti_2->Fill(pt_kf,imass);
                                if(fSPDMult_corr>=25 && fSPDMult_corr < 35) fHist_InvMass_pt_ULS_KF_SPDmulti_3->Fill(pt_kf,imass);
                                if(fSPDMult_corr>=35 && fSPDMult_corr < 45) fHist_InvMass_pt_ULS_KF_SPDmulti_4->Fill(pt_kf,imass);
                                if(fSPDMult_corr>=45 && fSPDMult_corr < 88) fHist_InvMass_pt_ULS_KF_SPDmulti_5->Fill(pt_kf,imass);
                                
                                if(fV0Mult_corr2>0 && fV0Mult_corr2<75)     fHist_InvMass_pt_ULS_KF_V0multi_1->Fill(pt_kf,imass);
                                if(fV0Mult_corr2>=75 && fV0Mult_corr2<150)  fHist_InvMass_pt_ULS_KF_V0multi_2->Fill(pt_kf,imass);
                                if(fV0Mult_corr2>=150 && fV0Mult_corr2<225) fHist_InvMass_pt_ULS_KF_V0multi_3->Fill(pt_kf,imass);
                                if(fV0Mult_corr2>=225 && fV0Mult_corr2<300) fHist_InvMass_pt_ULS_KF_V0multi_4->Fill(pt_kf,imass);
                                if(fV0Mult_corr2>=300 && fV0Mult_corr2<680) fHist_InvMass_pt_ULS_KF_V0multi_5->Fill(pt_kf,imass);
                                
                                //with weight
                                if(fSPDMult_corr>0 && fSPDMult_corr < 15)   fHist_InvMass_pt_ULS_KF_SPDmulti_1_weight->Fill(pt_kf,imass,weight);
                                if(fSPDMult_corr>=15 && fSPDMult_corr < 25) fHist_InvMass_pt_ULS_KF_SPDmulti_2_weight->Fill(pt_kf,imass,weight);
                                if(fSPDMult_corr>=25 && fSPDMult_corr < 35) fHist_InvMass_pt_ULS_KF_SPDmulti_3_weight->Fill(pt_kf,imass,weight);
                                if(fSPDMult_corr>=35 && fSPDMult_corr < 45) fHist_InvMass_pt_ULS_KF_SPDmulti_4_weight->Fill(pt_kf,imass,weight);
                                if(fSPDMult_corr>=45 && fSPDMult_corr < 88) fHist_InvMass_pt_ULS_KF_SPDmulti_5_weight->Fill(pt_kf,imass,weight);
                                
                                if(fV0Mult_corr2>0 && fV0Mult_corr2<75)     fHist_InvMass_pt_ULS_KF_V0multi_1_weight->Fill(pt_kf,imass,weight);
                                if(fV0Mult_corr2>=75 && fV0Mult_corr2<150)  fHist_InvMass_pt_ULS_KF_V0multi_2_weight->Fill(pt_kf,imass,weight);
                                if(fV0Mult_corr2>=150 && fV0Mult_corr2<225) fHist_InvMass_pt_ULS_KF_V0multi_3_weight->Fill(pt_kf,imass,weight);
                                if(fV0Mult_corr2>=225 && fV0Mult_corr2<300) fHist_InvMass_pt_ULS_KF_V0multi_4_weight->Fill(pt_kf,imass,weight);
                                if(fV0Mult_corr2>=300 && fV0Mult_corr2<680) fHist_InvMass_pt_ULS_KF_V0multi_5_weight->Fill(pt_kf,imass,weight);
                                
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
											
											if( TMath::Abs(pdg) == 11 && TMath::Abs(pdg) == 11 ) 
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
															
                                                            //printf("Label mother leg 1 %d, leg %d\n", fMCparticle->GetLabel(),fMCparticle2->GetLabel());
															fPtMCparticle_Total_e_from_JPsi->Fill(track->Pt()); //reconstructed pT
                                                            //checking if they are from same mother
                                                            if((fMCparticleMother->GetLabel())==(fMCparticleMother2->GetLabel())){
                                                               // printf("electrons from same mother\n\n");
                                                                fPtMCparticle_Total_e_from_JPsi_sameMother->Fill(track->Pt()); //reconstructed pT
                                                            }
															
															if(invmass3>=fMassCutMin && invmass3<=fMassCutMax){
                                                                fPtMCparticle_TotalplusMass_e_from_JPsi->Fill(track->Pt()); //reconstructed pT from electron
                                                                
                                                                fPtMCparticle_TotalplusMass_JPsi_pT->Fill(pt_kf);//spectrum of reconstructed J/Psi
                                                                
                                                                if((fMCparticleMother->GetLabel())==(fMCparticleMother2->GetLabel())){
                                                                   // printf("electrons from same mother + inv mass cut\n\n");
                                                                    fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother->Fill(track->Pt()); //reconstructed pT
                                                                    fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother->Fill(pt_kf);//spectrum of reconstructed J/Psi
                                                                    
                                                                }
															}
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
					
				}
				
			}//close lTrack
		
	//////////////////////////////////////////////////////////////////////////////////////////////////
			
		fVtxZ[2]->Fill(fZvtx);
		
		
		
        }
	
		
//=======================================================================
	
	delete fListOfmotherkink;
	PostData(1, fOutputList);
    PostData(2,fListProfiles);
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
//____________________________________________________________________________
TProfile* AliAnalysisTask_JPsi_EMCal::GetEstimatorHistogram(const AliAODEvent* fAOD)
{
    
   // printf("Inside 'GetEstimatorHistogram \n'");
    Int_t runNo  = fAOD->GetRunNumber();
    //cout<<"run number"<<runNo<<endl;
    Int_t period = -1;
    
    //period = 0;
    
    if (runNo>258883 && runNo<260187) period = 0;//16l  259668
    if (runNo>256504 && runNo<258574) period = 1;//16k
    if (period < 0 || period > 1) return 0;
    
   
    
   // cout<<"using period = 0 for all (just a test)"<<period<<endl;
    period = 0;
    
    return fMultEstimatorAvg[period];
}
//____________________________________________________________________________
TProfile2D* AliAnalysisTask_JPsi_EMCal::GetEstimatorHistogram_V0(const AliAODEvent* fAOD)
{
    
    //printf("Inside 'GetEstimatorHistogram_V0 \n'");
    Int_t runNo  = fAOD->GetRunNumber();
    //cout<<"run number"<<runNo<<endl;
    Int_t period = -1;
    
    //period = 0;
    
    if (runNo>258883 && runNo<260187) period = 0;//16l  259668
    if (runNo>256504 && runNo<258574) period = 1;//16k
    if (period < 0 || period > 1) return 0;
    
    
    
  //  cout<<"using period = 0 for all (just a test -- V0)"<<period<<endl;
    period = 0;
 
    return fMultEstimatorV0[period];
}
//______________________________________________________________________________
Double_t AliAnalysisTask_JPsi_EMCal::GetTrackletsMeanCorrection(TProfile* estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult)
{
    
    if(TMath::Abs(vtxZ)>10.0){
        //    printf("ERROR: Z vertex out of range for correction of multiplicity\n");
        return uncorrectedNacc;
    }
    
    if(!estimatorAvg){
        printf("ERROR: Missing TProfile for correction of multiplicity\n");
        return uncorrectedNacc;
    }
   
    
    Double_t localAvg = estimatorAvg->GetBinContent(estimatorAvg->FindBin(vtxZ));
    
    if(localAvg==0){
        //printf("LocalAvg = %f for vertex = %f  \n", localAvg, vtxZ);
        return uncorrectedNacc;
    }
    
    
    Double_t deltaM = 0;
    deltaM = uncorrectedNacc*(refMult/localAvg - 1);
    
   
    
    Double_t correctedNacc = uncorrectedNacc + (deltaM>0 ? 1 : -1) * gRandom->PoissonD(TMath::Abs(deltaM));
    
    if(correctedNacc<0) correctedNacc=0;
    
  //  printf("Inside 'GetTrackletsMeanCorrection' (end) \n");
    
    return correctedNacc;
    
}

//______________________________________________________________________________
Double_t AliAnalysisTask_JPsi_EMCal::GetV0MeanCorrection(TProfile2D* estimatorV0, Double_t uncorrectedV0, Double_t vtxZ, Double_t refMult_V0, Int_t run_number)
{
    
    
    if(TMath::Abs(vtxZ)>10.0){
        //    printf("ERROR: Z vertex out of range for correction of multiplicity\n");
        return uncorrectedV0;
    }
    
    
    if(!estimatorV0){
        printf("ERROR: Missing TProfile for correction of multiplicity\n");
        return uncorrectedV0;
    }
  
    TH1F *hx= (TH1F*) estimatorV0->ProfileX();
    TH1F *hy= (TH1F*) estimatorV0->ProfileY();
    
   // Int_t BinX=(estimatorV0->ProfileX())->FindBin(run_number);
   // Int_t BinY=(estimatorV0->ProfileY())->FindBin(vtxZ);
    
    Int_t BinX=hx->FindBin(run_number);
    Int_t BinY=hy->FindBin(vtxZ);
    
    delete hx;
    delete hy;
    
    
    Double_t localV0 = estimatorV0->GetBinContent(BinX, BinY);//first argument is run number/group, second is vertex on Y axis
   // printf("localV0 = %f\n", localV0);
    
    if(localV0==0){
       // printf("LocalV0 = 0 for vertex = %f and  run =  %d \n", vtxZ, run_number);
        return uncorrectedV0;
    }
    
    Double_t deltaMV0 = 0;
    deltaMV0 = uncorrectedV0*(refMult_V0/localV0 - 1);
    
    //printf("Inside 'GetV0MeanCorrection' (after deltaM) \n");
    
    Double_t correctedV0 = uncorrectedV0 + (deltaMV0>0 ? 1 : -1) * gRandom_V0->PoissonD(TMath::Abs(deltaMV0));
    
    if(correctedV0<0) correctedV0=0;
    
 //   printf("Inside 'GetV0MeanCorrection' (end) \n");
    
    return correctedV0;
    
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

