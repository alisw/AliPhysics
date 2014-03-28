
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
	//		version: March 23, 2014.								      //
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
,fUseShowerShapeCut(kFALSE)
,fFillBackground(kFALSE)
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
,fPtElec_ULS_weight(0)
,fPtElec_LS_weight(0)
,fPtElec_ULS2_weight(0)
,fPtElec_LS2_weight(0)
,fTOF01(0)
,fTOF02(0)
,fTOF03(0)
,fpid(0)
,fEoverP_pt(0)
,fEoverP_tpc(0)
,fTPC_pt(0)
,fTPC_p(0)
,fTPCnsigma_pt(0)
,fTPCnsigma_p(0)
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
,fEtaPhi(0)
,fEtaPhi_num(0)
,fEtaPhi_den(0)
,fpt_reco_pt_MC_num(0)
,fpt_reco_pt_MC_den(0)
,fVtxZ(0)

,fVtxZ_new1(0)
,fVtxZ_new2(0)
,fVtxZ_new3(0)
,fVtxZ_new4(0)

,fEtad(0)
,fNTracks(0)
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
,fDCAcut(999)
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
{
		//Named constructor 
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
AliAnalysisTaskEMCalHFEpA::AliAnalysisTaskEMCalHFEpA() 
: AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTaskEMCalHFEpA")
,fCorrelationFlag(0)
,fIsMC(0)
,fUseEMCal(kFALSE)
,fUseTrigger(kFALSE)
,fUseShowerShapeCut(kFALSE)
,fFillBackground(kFALSE)
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
,fPtElec_ULS_weight(0)
,fPtElec_LS_weight(0)
,fPtElec_ULS2_weight(0)
,fPtElec_LS2_weight(0)
,fTOF01(0)
,fTOF02(0)
,fTOF03(0)
,fpid(0)
,fEoverP_pt(0)
,fEoverP_tpc(0)
,fTPC_pt(0)
,fTPC_p(0)
,fTPCnsigma_pt(0)
,fTPCnsigma_p(0)
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
,fEtaPhi(0)
,fEtaPhi_num(0)
,fEtaPhi_den(0)
,fpt_reco_pt_MC_num(0)
,fpt_reco_pt_MC_den(0)
,fVtxZ(0)

,fVtxZ_new1(0)
,fVtxZ_new2(0)
,fVtxZ_new3(0)
,fVtxZ_new4(0)

,fEtad(0)
,fNTracks(0)
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
,fDCAcut(999)
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
{
		// Constructor
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
AliAnalysisTaskEMCalHFEpA::~AliAnalysisTaskEMCalHFEpA()
{
	//Destructor 
	delete fOutputList;
	delete fPID;
	delete fCFM;
	delete fPIDqa;
	//Lucile
	//delete reader; 
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
	fNevent = new TH1F("fNevent","Number of Events",15,0,15);
		//And then, add to the output list
	fOutputList->Add(fNevent);
	
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

	}
	
	fPtTrigger_Inc = new TH1F("fPtTrigger_Inc","pT dist for Hadron Contamination; p_{t} (GeV/c); Count",300,0,30);
	fTPCnsigma_pt_2D = new TH2F("fTPCnsigma_pt_2D",";pt (GeV/c);TPC Electron N#sigma",1000,0.3,30,1000,-15,10);
	fShowerShapeCut = new TH2F("fShowerShapeCut","Shower Shape;M02;M20",500,0,1.8,500,0,1.8);
	fEtaPhi_num=new TH2F("fEtaPhi_num","#eta x #phi track;#phi;#eta",200,0.,5,50,-1.,1.);
	fEtaPhi_den=new TH2F("fEtaPhi_den","#eta x #phi track;#phi;#eta",200,0.,5,50,-1.,1.);
		
	fpt_reco_pt_MC_num=new TH2F("fpt_reco_pt_MC_num","pt reco x pt MC;pt reco; pt MC",300,0.,30,300,0.,30);
	fpt_reco_pt_MC_den=new TH2F("fpt_reco_pt_MC_den","pt reco x pt MC;pt reco; pt MC",300,0.,30,300,0.,30);
	
	
	fCharge_n = new TH1F("fCharge_n","Inclusive Electrons (Negative Charge); p_{t} (GeV/c); Count",200,0,30);
	fCharge_p = new TH1F("fCharge_p","Inclusive Positrons (Positive Charge); p_{t} (GeV/c); Count",200,0,30);
	
	fECluster_pure= new TH1F("fECluster_pure", ";ECluster pure",2000,0,100);
	
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
	}
	
	
	fOutputList->Add(fPtTrigger_Inc);
	fOutputList->Add(fTPCnsigma_pt_2D);
	fOutputList->Add(fShowerShapeCut);
	
	fOutputList->Add(fCharge_n);
	fOutputList->Add(fCharge_p);
	
	fOutputList->Add(fECluster_pure);
	
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
	
	fVtxZ_new1= new  TH1F("fVtxZ_new1","fVtxZ_new1",1000, -50,50);
	fVtxZ_new2= new  TH1F("fVtxZ_new2","fVtxZ_new2",1000, -50,50);
	fVtxZ_new3= new  TH1F("fVtxZ_new3","fVtxZ_new3",1000, -50,50);
	fVtxZ_new4= new  TH1F("fVtxZ_new4","fVtxZ_new4",1000, -50,50);
	
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
	
	
	
	for(Int_t i = 0; i < 3; i++)
	{
		fEoverP_pt[i] = new TH2F(Form("fEoverP_pt%d",i),";p_{t} (GeV/c);E / p ",1000,0,30,500,0,2);
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
	
	for(Int_t i = 0; i < 4; i++)
	{
		fTPCNcls_pid[i]= new TH2F(Form("fTPCNcls_pid%d",i),"fTPCNcls_pid;NCls;NCls for PID",200,0,200,200,0,200);
		fOutputList->Add(fTPCNcls_pid[i]);
	}
	
		//pt bin
	Int_t fPtBin[7] = {1,2,4,6,8,10,15};
	
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
	fEoverP_pt_pions= new TH2F("fEoverP_pt_pions","fEoverP_pt_pions",1000,0,30,500,0,2);
	
	ftpc_p_EoverPcut= new TH2F("ftpc_p_EoverPcut","ftpc_p_EoverPcut",1000,0,30,200,20,200);
	fnsigma_p_EoverPcut= new TH2F("fnsigma_p_EoverPcut","fnsigma_p_EoverPcut",1000,0,30,500,-15,15);
	
	fEoverP_pt_pions2= new TH2F("fEoverP_pt_pions2","fEoverP_pt_pions2",1000,0,30,500,0,2);
	fEoverP_pt_hadrons= new TH2F("fEoverP_pt_hadrons","fEoverP_pt_hadrons",1000,0,30,500,0,2);
	
	
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
	
		//______________________________________________________________________
		//Vertex Selection
	if(fIsAOD)
	{
			
		const AliAODVertex* trkVtx = fAOD->GetPrimaryVertex();
		if(!trkVtx || trkVtx->GetNContributors()<=0) return;
		TString vtxTtl = trkVtx->GetTitle();
		if(!vtxTtl.Contains("VertexerTracks")) return;
		Float_t zvtx = trkVtx->GetZ();
		fZvtx = zvtx;
		
		fVtxZ_new1->Fill(fZvtx);
		
		const AliAODVertex* spdVtx = fAOD->GetPrimaryVertexSPD();
		if(spdVtx->GetNContributors()<=0) return;
		TString vtxTyp = spdVtx->GetTitle();
		Double_t cov[6]={0};
		spdVtx->GetCovarianceMatrix(cov);
		Double_t zRes = TMath::Sqrt(cov[5]);
		if(vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return;
		if(TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return;
		if(TMath::Abs(zvtx) > 10) return;
		
		fVtxZ_new2->Fill(fZvtx);
		
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
		Float_t zvtx = trkVtx->GetZ();
		
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
	
	
		
		//Only events with at least 2 tracks are accepted
	Int_t fNOtrks =  fVevent->GetNumberOfTracks();
	
	//if(fIsAOD) Int_t fNOtrks =  fAOD->GetNumberOfTracks();
	//if(!fIsAOD) Int_t fNOtrks =  fESD->GetNumberOfTracks();
	
		//commented to test
	if(fNOtrks<2) return;

	
	//new track loop to select events
	
	fNevent->Fill(13);
	/*
	if(fUseTrigger){
		if(fIsAOD){
			double fTrackMulti=0;
			for(Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++) 
			{
				AliVParticle* Vtrack = fVevent->GetTrack(iTracks);
				if (!Vtrack) 
				{
					printf("ERROR: Could not receive track %d\n", iTracks);
					continue;
				}
		
				AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
					//AliESDtrack *etrack = dynamic_cast<AliESDtrack*>(Vtrack);
				AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);
		
				if((track->Pt())<0.2 || (track->Pt())>1000.0) continue;
			   	//if it is not a hybrid track, continue
				if(!atrack->TestFilterBit(768)) continue;
				else fTrackMulti=fTrackMulti+1;
		
			}
				//Only take event if track multiplicity is bigger than 2.
			if(fTrackMulti<2) return;
		}
	}
	 */
	
	fVtxZ_new3->Fill(fZvtx);
	
	fNevent->Fill(14);
	
		
	//trying to use same as Lucile
	/*if(fIsAOD) {
			//reader = new AliCaloTrackAODReader();
		reader->SwitchOnCTS();
		reader->SetCTSPtMin(0.2);
		reader->SetCTSPtMax(1000);	
		reader->SwitchOffRecalculateVertexBC();
		reader->SwitchOffVertexBCEventSelection();
		reader->SwitchOffUseTrackTimeCut();	
		reader->SwitchOffUseTrackDCACut();	
		reader->SwitchOnAODHybridTrackSelection();
		reader->SwitchOnRejectNoTrackEvents();
	}*/
	
	
		//______________________________________________________________________
		//Centrality Selection
	if(fHasCentralitySelection)
	{
		Float_t centrality = -1;
		
		if(fIsAOD) 
		{
			fCentrality = fAOD->GetHeader()->GetCentralityP();
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
					
					if (TMath::Abs(pdg) == 11 && fMCparticle->IsPhysicalPrimary()) fPtMCparticleAlle_Primary->Fill(fMCparticle->Pt()); //denominator for total efficiency for all electrons primary
					
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
								if(fIsHFE2) fPtMCparticleAllHfe2->Fill(fMCparticle->Pt());
							}
						}
					}
				}//eta cut
				
				//only primary pions 
				if(fMCparticle->IsPhysicalPrimary()){
					if(TMath::Abs(pdg)==111) fPtMCpi0->Fill(fMCparticle->Pt());
					if(TMath::Abs(pdg)==221) fPtMCeta->Fill(fMCparticle->Pt());
					
					if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax)
					{
						
						if(TMath::Abs(pdg)==111) fPtMCpi02->Fill(fMCparticle->Pt());
						if(TMath::Abs(pdg)==221) fPtMCeta2->Fill(fMCparticle->Pt());
						
					}
					
				}
				
			}//loop tracks
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
								if(fIsHFE2) fPtMCparticleAllHfe2->Fill(fMCtrack->Pt());
							}
						}//Is Physical primary
					}	
				}//eta cut
	        }//loop tracks
		}//ESD
	}//Is MC
	
	//______________________________________________________________________
	//EMCal Trigger Selection (Threshold selection)
	TString firedTrigger;
	TString TriggerEG1("CEMC7EG1"); //takes trigger with name with EG1, ex: CEMC7EG1-B-NOPF-CENTNOTRD  
	TString TriggerEG2("CEMC7EG2");
	//Jan 17, 2014
	TString TriggerEJE("EJE");
		
	if(fAOD) firedTrigger = fAOD->GetFiredTriggerClasses();
	else if(fESD) firedTrigger = fESD->GetFiredTriggerClasses();
	
	fNevent->Fill(0);
	if(firedTrigger.Contains(TriggerEG1)) fNevent->Fill(1);
	if(firedTrigger.Contains(TriggerEG2)) fNevent->Fill(2);
	
	// Jan 06, 2014: I changed the counters: Only fill with 3 or 4 if we want the trigger threshold selected.
	//EG1
	if(firedTrigger.Contains(TriggerEG1))
	{ 
		if(fEMCEG1) fNevent->Fill(3);
	}
	else 
	{
		if(fEMCEG1) return;
	}
	
	//EG2
	if(firedTrigger.Contains(TriggerEG2))
	{ 
		if(fEMCEG2) fNevent->Fill(4);
	}
	else
	{ 
		if(fEMCEG2) return;
	}
	
	//______________________________________________________________________
	//Testing if there is an overlap EGA and EJE
	//none
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
	
		
	//New cluster information 
	//after trigger threshold selection
	Int_t ClsNo2 = fVevent->GetNumberOfCaloClusters(); 
	
	if(ClsNo2<=0){
		fNevent->Fill(11); //events with no cluster
		return;
	} 
	for (Int_t i=0; i< ClsNo2; i++ ){
		
		fClus = fVevent->GetCaloCluster(i);
		if(fClus->IsEMCAL())
		{
			//pure cluster information
			fECluster_pure->Fill(fClus->E());
		}
	}
	
	fNevent->Fill(12); //events with cluster
	
	
	fVtxZ_new4->Fill(fZvtx);
	
	//__________________________________________________________________
	
	Int_t ClsNo = -999;
	if(!fIsAOD) ClsNo = fESD->GetNumberOfCaloClusters(); 
	else ClsNo = fAOD->GetNumberOfCaloClusters(); 
	
	//______________________________________________________________________
	
	///_____________________________________________________________________
	///Track loop
	for(Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++) 
	{
		AliVParticle* Vtrack = fVevent->GetTrack(iTracks);
		if (!Vtrack) 
		{
			printf("ERROR: Could not receive track %d\n", iTracks);
			continue;
		}
		
		AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
		AliESDtrack *etrack = dynamic_cast<AliESDtrack*>(Vtrack);
		AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);
		
		//aod test -- Francesco suggestion
		AliAODTrack *aod_track=fAOD->GetTrack(iTracks);
		
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
		fEtad[1]->Fill(track->Eta());
		
			
		
		
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
			fClus = fVevent->GetCaloCluster(track->GetEMCALcluster());
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
		if(iTracks == 0)fNTracks[0]->Fill(fNOtrks);
		fNTracks_pt[0]->Fill(fNOtrks, fPt);
		fNTracks_eta[0]->Fill(fNOtrks, track->Eta());
		fNTracks_phi[0]->Fill(fNOtrks, track->Phi());

		
		fNClusters[0]->Fill(ClsNo);
		fTPCNcls_pid[0]->Fill(TPCNcls, TPCNcls_pid);
			//______________________________________________________________
		
///Fill QA plots without track selection
///_____________________________________________________________________________
//______________________________________________________________________________________
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
		if(!fIsAOD)
		{
			if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;
		}
		
//HFEcuts: ITS layers cuts
		if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;
		
//HFE cuts: TPC PID cleanup
		if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;
//______________________________________________________________________________________
		
		if(fIsAOD){	
				//AOD test -- Fancesco suggestion
			Int_t type=aod_track->GetType();
			if(type==AliAODTrack::kPrimary) fPtPrim->Fill(aod_track->Pt());
			if(type==AliAODTrack::kSecondary) fPtSec->Fill(aod_track->Pt());
		
				//Int_t type2=track->GetType();
			if(type==AliAODTrack::kPrimary) fPtPrim2->Fill(track->Pt());
			if(type==AliAODTrack::kSecondary) fPtSec2->Fill(track->Pt());
		}
			
		
///_____________________________________________________________
///QA plots after track selection
		if(fIsMC)
		{
			if(track->GetLabel()>=0) fPtMCWithLabel->Fill(fPt);
			else fPtMCWithoutLabel->Fill(fPt);
		}
		
		if(fIsMC  && track->GetLabel()>=0)
		{
			if(fIsAOD){
				fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
			
											
				if(fMCparticle->IsPhysicalPrimary()) fPtIsPhysicaPrimary->Fill(fPt);
			
				Int_t pdg = fMCparticle->GetPdgCode();
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
								if(fIsHFE1){
									fPtMCparticleRecoHfe1->Fill(fMCparticle->Pt());//numerator tracking
									//unfolding
									fpt_reco_pt_MC_den->Fill(track->Pt(),fMCparticle->Pt());

								}
								if(fIsHFE2) fPtMCparticleRecoHfe2->Fill(fMCparticle->Pt());
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
		
		TPCNcls = track->GetTPCNcls();
		Float_t pos2[3]={0,0,0};
		
		if(track->GetEMCALcluster()>0)
		{
			fClus = fVevent->GetCaloCluster(track->GetEMCALcluster());
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
							//main
						if(TMath::Abs(fTPCnSigma_pion)<3 || TMath::Abs(fTPCnSigma_proton)<3 || TMath::Abs(fTPCnSigma_kaon)<3 ){
							
							if(fTPCnSigma<-3.5){
								fEoverP_pt_hadrons->Fill(fPt,EoverP);
								if(fUseEMCal) fShowerShape_ha->Fill(M02,M20);
							}
						}
							//for systematic studies of hadron contamination
						if(fTPCnSigma < -4){
							fEoverP_pt_pions->Fill(fPt, EoverP);
							
						}
						
						if(fTPCnSigma < -5){
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
		if(iTracks == 0)fNTracks[1]->Fill(fNOtrks);
		fNTracks_pt[1]->Fill(fNOtrks, fPt);
		fNTracks_eta[1]->Fill(fNOtrks, track->Eta());
		fNTracks_phi[1]->Fill(fNOtrks, track->Phi());
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
					
					// EtaCut electrons histogram
					//Shower Shape Cut
					if(track->Eta()>=fEtaCutMin && track->Eta()<=fEtaCutMax ){
						
						if(fUseShowerShapeCut){
							if(M02 >= fM02CutMin && M02<=fM02CutMax && M20>=fM20CutMin && M20<=fM20CutMax){
								fEoverP_pt[2]->Fill(fPt,(fClus->E() / fP));
								fShowerShapeCut->Fill(M02,M20);
								
							}
							
						}
						if(!fUseShowerShapeCut){
							fEoverP_pt[2]->Fill(fPt,(fClus->E() / fP));
							fShowerShapeCut->Fill(M02,M20);
							
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
					if((fClus->E() / fP) >= fEoverPCutMin && (fClus->E() / fP) <= fEoverPCutMax)
					{	
						
					    fECluster[2]->Fill(Energy);
						fTPCNcls_pid[3]->Fill(TPCNcls, TPCNcls_pid);
						
						if(fUseEMCal)
						{
							fPtElec_Inc->Fill(fPt);
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
							if(fMCtrack->GetFirstMother()>0) fMCtrackMother = fMCstack->Particle(fMCtrack->GetFirstMother());
							TParticle *particle=fMCstack->Particle(track->GetLabel());

							Int_t pdg = fMCtrack->GetPdgCode();
							
							
							if(fMCtrack->Eta()>=fEtaCutMin && fMCtrack->Eta()<=fEtaCutMax){
								if( TMath::Abs(pdg) == 11 && fMCtrack->GetFirstMother()>0 ){
									Int_t mpdg = fMCtrackMother->GetPdgCode();
									if(TMath::Abs(mpdg) == 221 || TMath::Abs(mpdg) == 22 || TMath::Abs(mpdg) == 111){
										Double_t proR=particle->R();
										if(proR<7){
										  fPtMCelectronAfterAll_nonPrimary->Fill(fMCtrack->Pt()); //numerator for the total efficiency, non Primary track
										}
									}
								}
								if( TMath::Abs(pdg) == 11 && fMCstack->IsPhysicalPrimary(track->GetLabel())) fPtMCelectronAfterAll_Primary->Fill(fMCtrack->Pt());
							}
							
							if(fMCstack->IsPhysicalPrimary(track->GetLabel()))
							{
								Bool_t MotherFound = FindMother(track->GetLabel());
							    
								if(MotherFound)
								{
									if(fMCtrack->Eta()>=fEtaCutMin && fMCtrack->Eta()<=fEtaCutMax){
										
										
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
											
												if(fIsHFE1)fPtMC_EMCal_Selected->Fill(fMCtrack->Pt());	
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
		if(iTracks == 0)fNTracks[2]->Fill(fNOtrks);
		fNTracks_pt[2]->Fill(fNOtrks, fPt);
		fNTracks_eta[2]->Fill(fNOtrks, track->Eta());
		fNTracks_phi[2]->Fill(fNOtrks, track->Phi());
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
	Bool_t IsMCefix=kFALSE; //to make correction on efix, use kTRUE (do not change the efficiency, so I will keep the correction only for d3)
		
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
							//correction for efix based on data - parametrization from MinJung
						if(IsMCefix){
							if(TMath::Abs(fMCparticleMother->GetPdgCode())==111){
								Double_t x=mPt;
								if(0.100000 <= x < 0.112797 ) mweight=1.030419;
								if(0.112797 <= x < 0.127231 ) mweight=1.044554;
								if(0.127231 <= x < 0.143512 ) mweight=1.062733;
								if(0.143512 <= x < 0.161877 ) mweight=1.085332;
								if(0.161877 <= x < 0.182592 ) mweight=1.115248;
								if(0.182592 <= x < 0.205957 ) mweight=1.153990;
								if(0.205957 <= x < 0.232313 ) mweight=1.201346;
								if(0.232313 <= x < 0.262041 ) mweight=1.257332;
								if(0.262041 <= x < 0.295573 ) mweight=1.315488;
								if(0.295573 <= x < 0.333397 ) mweight=1.369138;
								if(0.333397 <= x < 0.376060 ) mweight=1.407632;
								if(0.376060 <= x < 0.424183 ) mweight=1.422232;
								if(0.424183 <= x < 0.478465 ) mweight=1.406922;
								if(0.478465 <= x < 0.539692 ) mweight=1.360082;
								if(0.539692 <= x < 0.608754 ) mweight=1.284405;
								if(0.608754 <= x < 0.686654 ) mweight=1.182017;
								if(0.686654 <= x < 0.774523 ) mweight=1.062002;
								if(0.774523 <= x < 0.873636 ) mweight=0.935533;
								if(0.873636 <= x < 0.985432 ) mweight=0.816081;
								if(0.985432 <= x < 1.111534 ) mweight=0.717527;
								if(1.111534 <= x < 1.253773 ) mweight=0.647465;
								if(1.253773 <= x < 1.414214 ) mweight=0.607212;
								if(1.414214 <= x < 1.595185 ) mweight=0.589750;
								if(1.595185 <= x < 1.799315 ) mweight=0.587406;
								if(1.799315 <= x < 2.029567 ) mweight=0.592858;
								if(2.029567 <= x < 2.289283 ) mweight=0.601059;
								if(2.289283 <= x < 2.582235 ) mweight=0.608003;
								if(2.582235 <= x < 2.912674 ) mweight=0.611705;
								if(2.912674 <= x < 3.285398 ) mweight=0.610086;
								if(3.285398 <= x < 3.705818 ) mweight=0.605015;
								if(3.705818 <= x < 4.180038 ) mweight=0.596299;
								if(4.180038 <= x < 4.714942 ) mweight=0.590727;
								if(4.714942 <= x < 5.318296 ) mweight=0.585358;
								if(5.318296 <= x < 5.998859 ) mweight=0.585257;
								if(5.998859 <= x < 6.766511 ) mweight=0.580812;
								if(6.766511 <= x < 7.632396 ) mweight=0.576207;
								if(7.632396 <= x < 8.609086 ) mweight=0.575912;
								if(8.609086 <= x < 9.710759 ) mweight=0.558718;
								if(9.710759 <= x < 10.953409 ) mweight=0.555625;
								if(10.953409 <= x < 12.355077 ) mweight=0.558886;
								if(12.355077 <= x < 13.936111 ) mweight=0.545318;
								if(13.936111 <= x < 15.719464 ) mweight=0.517607;
								if(15.719464 <= x < 17.731026 ) mweight=0.512366;
								if(17.731026 <= x < 20.000000 ) mweight=0.497034;
								
								
							}
							
							
						}//end of IsMCefix
						
							//________________________________________________________________
							//correction for d3 based on data //from Jan
						if(!IsMCefix){
							if(TMath::Abs(fMCparticleMother->GetPdgCode())==111){
								Double_t x=mPt;
								if(0.100000 <= x < 0.112797 ) mweight=1.262120;
								if(0.112797 <= x < 0.127231 ) mweight=1.277765;
								if(0.127231 <= x < 0.143512 ) mweight=1.295605;
								if(0.143512 <= x < 0.161877 ) mweight=1.318155;
								if(0.161877 <= x < 0.182592 ) mweight=1.348693;
								if(0.182592 <= x < 0.205957 ) mweight=1.388636;
								if(0.205957 <= x < 0.232313 ) mweight=1.439122;
								if(0.232313 <= x < 0.262041 ) mweight=1.497452;
								if(0.262041 <= x < 0.295573 ) mweight=1.559409;
								if(0.295573 <= x < 0.333397 ) mweight=1.615169;
								if(0.333397 <= x < 0.376060 ) mweight=1.654954;
								if(0.376060 <= x < 0.424183 ) mweight=1.668753;
								if(0.424183 <= x < 0.478465 ) mweight=1.652225;
								if(0.478465 <= x < 0.539692 ) mweight=1.603119;
								if(0.539692 <= x < 0.608754 ) mweight=1.526049;
								if(0.608754 <= x < 0.686654 ) mweight=1.426724;
								if(0.686654 <= x < 0.774523 ) mweight=1.312684;
								if(0.774523 <= x < 0.873636 ) mweight=1.195395;
								if(0.873636 <= x < 0.985432 ) mweight=1.086264;
								if(0.985432 <= x < 1.111534 ) mweight=0.993666;
								if(1.111534 <= x < 1.253773 ) mweight=0.922587;
								if(1.253773 <= x < 1.414214 ) mweight=0.875739;
								if(1.414214 <= x < 1.595185 ) mweight=0.852181;
								if(1.595185 <= x < 1.799315 ) mweight=0.847828;
								if(1.799315 <= x < 2.029567 ) mweight=0.863875;
								if(2.029567 <= x < 2.289283 ) mweight=0.899112;
								if(2.289283 <= x < 2.582235 ) mweight=0.955194;
								if(2.582235 <= x < 2.912674 ) mweight=1.033824;
								if(2.912674 <= x < 3.285398 ) mweight=1.133714;
								if(3.285398 <= x < 3.705818 ) mweight=1.259471;
								if(3.705818 <= x < 4.180038 ) mweight=1.406883;
								if(4.180038 <= x < 4.714942 ) mweight=1.578923;
								if(4.714942 <= x < 5.318296 ) mweight=1.778513;
								if(5.318296 <= x < 5.998859 ) mweight=2.001171;
								if(5.998859 <= x < 6.766511 ) mweight=2.223161;
								if(6.766511 <= x < 7.632396 ) mweight=2.449445;
								if(7.632396 <= x < 8.609086 ) mweight=2.661734;
								if(8.609086 <= x < 9.710759 ) mweight=2.851935;
								if(9.710759 <= x < 10.953409 ) mweight=2.974319;
								if(10.953409 <= x < 12.355077 ) mweight=3.106314;
								if(12.355077 <= x < 13.936111 ) mweight=3.126815;
								if(13.936111 <= x < 15.719464 ) mweight=3.150053;
								if(15.719464 <= x < 17.731026 ) mweight=3.218509;
								if(17.731026 <= x < 20.000000 ) mweight=3.252141;
								
								
							}
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
							//correction for efix based on data - parametrization from MinJung
						if(IsMCefix){
							if(TMath::Abs(fMCparticleGMother->GetPdgCode())==111){
								Double_t x=gmPt;
								if(0.100000 <= x < 0.112797 ) gmweight=1.030419;
								if(0.112797 <= x < 0.127231 ) gmweight=1.044554;
								if(0.127231 <= x < 0.143512 ) gmweight=1.062733;
								if(0.143512 <= x < 0.161877 ) gmweight=1.085332;
								if(0.161877 <= x < 0.182592 ) gmweight=1.115248;
								if(0.182592 <= x < 0.205957 ) gmweight=1.153990;
								if(0.205957 <= x < 0.232313 ) gmweight=1.201346;
								if(0.232313 <= x < 0.262041 ) gmweight=1.257332;
								if(0.262041 <= x < 0.295573 ) gmweight=1.315488;
								if(0.295573 <= x < 0.333397 ) gmweight=1.369138;
								if(0.333397 <= x < 0.376060 ) gmweight=1.407632;
								if(0.376060 <= x < 0.424183 ) gmweight=1.422232;
								if(0.424183 <= x < 0.478465 ) gmweight=1.406922;
								if(0.478465 <= x < 0.539692 ) gmweight=1.360082;
								if(0.539692 <= x < 0.608754 ) gmweight=1.284405;
								if(0.608754 <= x < 0.686654 ) gmweight=1.182017;
								if(0.686654 <= x < 0.774523 ) gmweight=1.062002;
								if(0.774523 <= x < 0.873636 ) gmweight=0.935533;
								if(0.873636 <= x < 0.985432 ) gmweight=0.816081;
								if(0.985432 <= x < 1.111534 ) gmweight=0.717527;
								if(1.111534 <= x < 1.253773 ) gmweight=0.647465;
								if(1.253773 <= x < 1.414214 ) gmweight=0.607212;
								if(1.414214 <= x < 1.595185 ) gmweight=0.589750;
								if(1.595185 <= x < 1.799315 ) gmweight=0.587406;
								if(1.799315 <= x < 2.029567 ) gmweight=0.592858;
								if(2.029567 <= x < 2.289283 ) gmweight=0.601059;
								if(2.289283 <= x < 2.582235 ) gmweight=0.608003;
								if(2.582235 <= x < 2.912674 ) gmweight=0.611705;
								if(2.912674 <= x < 3.285398 ) gmweight=0.610086;
								if(3.285398 <= x < 3.705818 ) gmweight=0.605015;
								if(3.705818 <= x < 4.180038 ) gmweight=0.596299;
								if(4.180038 <= x < 4.714942 ) gmweight=0.590727;
								if(4.714942 <= x < 5.318296 ) gmweight=0.585358;
								if(5.318296 <= x < 5.998859 ) gmweight=0.585257;
								if(5.998859 <= x < 6.766511 ) gmweight=0.580812;
								if(6.766511 <= x < 7.632396 ) gmweight=0.576207;
								if(7.632396 <= x < 8.609086 ) gmweight=0.575912;
								if(8.609086 <= x < 9.710759 ) gmweight=0.558718;
								if(9.710759 <= x < 10.953409 ) gmweight=0.555625;
								if(10.953409 <= x < 12.355077 ) gmweight=0.558886;
								if(12.355077 <= x < 13.936111 ) gmweight=0.545318;
								if(13.936111 <= x < 15.719464 ) gmweight=0.517607;
								if(15.719464 <= x < 17.731026 ) gmweight=0.512366;
								if(17.731026 <= x < 20.000000 ) gmweight=0.497034;
								
								
								
							}
							
							
						}//end of IsMCefix
						
							//________________________________________________________________
							//correction for d3 based on data
						if(!IsMCefix){
							if(TMath::Abs(fMCparticleGMother->GetPdgCode())==111){
								Double_t x=gmPt;
								if(0.100000 <= x < 0.112797 ) gmweight=1.262120;
								if(0.112797 <= x < 0.127231 ) gmweight=1.277765;
								if(0.127231 <= x < 0.143512 ) gmweight=1.295605;
								if(0.143512 <= x < 0.161877 ) gmweight=1.318155;
								if(0.161877 <= x < 0.182592 ) gmweight=1.348693;
								if(0.182592 <= x < 0.205957 ) gmweight=1.388636;
								if(0.205957 <= x < 0.232313 ) gmweight=1.439122;
								if(0.232313 <= x < 0.262041 ) gmweight=1.497452;
								if(0.262041 <= x < 0.295573 ) gmweight=1.559409;
								if(0.295573 <= x < 0.333397 ) gmweight=1.615169;
								if(0.333397 <= x < 0.376060 ) gmweight=1.654954;
								if(0.376060 <= x < 0.424183 ) gmweight=1.668753;
								if(0.424183 <= x < 0.478465 ) gmweight=1.652225;
								if(0.478465 <= x < 0.539692 ) gmweight=1.603119;
								if(0.539692 <= x < 0.608754 ) gmweight=1.526049;
								if(0.608754 <= x < 0.686654 ) gmweight=1.426724;
								if(0.686654 <= x < 0.774523 ) gmweight=1.312684;
								if(0.774523 <= x < 0.873636 ) gmweight=1.195395;
								if(0.873636 <= x < 0.985432 ) gmweight=1.086264;
								if(0.985432 <= x < 1.111534 ) gmweight=0.993666;
								if(1.111534 <= x < 1.253773 ) gmweight=0.922587;
								if(1.253773 <= x < 1.414214 ) gmweight=0.875739;
								if(1.414214 <= x < 1.595185 ) gmweight=0.852181;
								if(1.595185 <= x < 1.799315 ) gmweight=0.847828;
								if(1.799315 <= x < 2.029567 ) gmweight=0.863875;
								if(2.029567 <= x < 2.289283 ) gmweight=0.899112;
								if(2.289283 <= x < 2.582235 ) gmweight=0.955194;
								if(2.582235 <= x < 2.912674 ) gmweight=1.033824;
								if(2.912674 <= x < 3.285398 ) gmweight=1.133714;
								if(3.285398 <= x < 3.705818 ) gmweight=1.259471;
								if(3.705818 <= x < 4.180038 ) gmweight=1.406883;
								if(4.180038 <= x < 4.714942 ) gmweight=1.578923;
								if(4.714942 <= x < 5.318296 ) gmweight=1.778513;
								if(5.318296 <= x < 5.998859 ) gmweight=2.001171;
								if(5.998859 <= x < 6.766511 ) gmweight=2.223161;
								if(6.766511 <= x < 7.632396 ) gmweight=2.449445;
								if(7.632396 <= x < 8.609086 ) gmweight=2.661734;
								if(8.609086 <= x < 9.710759 ) gmweight=2.851935;
								if(9.710759 <= x < 10.953409 ) gmweight=2.974319;
								if(10.953409 <= x < 12.355077 ) gmweight=3.106314;
								if(12.355077 <= x < 13.936111 ) gmweight=3.126815;
								if(13.936111 <= x < 15.719464 ) gmweight=3.150053;
								if(15.719464 <= x < 17.731026 ) gmweight=3.218509;
								if(17.731026 <= x < 20.000000 ) gmweight=3.252141;
								
							}
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
							if(IsMCefix){
								if(TMath::Abs(fMCparticleMother->GetPdgCode())==111){
									Double_t x=mPt;
									if(0.100000 <= x < 0.112797 ) weight=1.030419;
									if(0.112797 <= x < 0.127231 ) weight=1.044554;
									if(0.127231 <= x < 0.143512 ) weight=1.062733;
									if(0.143512 <= x < 0.161877 ) weight=1.085332;
									if(0.161877 <= x < 0.182592 ) weight=1.115248;
									if(0.182592 <= x < 0.205957 ) weight=1.153990;
									if(0.205957 <= x < 0.232313 ) weight=1.201346;
									if(0.232313 <= x < 0.262041 ) weight=1.257332;
									if(0.262041 <= x < 0.295573 ) weight=1.315488;
									if(0.295573 <= x < 0.333397 ) weight=1.369138;
									if(0.333397 <= x < 0.376060 ) weight=1.407632;
									if(0.376060 <= x < 0.424183 ) weight=1.422232;
									if(0.424183 <= x < 0.478465 ) weight=1.406922;
									if(0.478465 <= x < 0.539692 ) weight=1.360082;
									if(0.539692 <= x < 0.608754 ) weight=1.284405;
									if(0.608754 <= x < 0.686654 ) weight=1.182017;
									if(0.686654 <= x < 0.774523 ) weight=1.062002;
									if(0.774523 <= x < 0.873636 ) weight=0.935533;
									if(0.873636 <= x < 0.985432 ) weight=0.816081;
									if(0.985432 <= x < 1.111534 ) weight=0.717527;
									if(1.111534 <= x < 1.253773 ) weight=0.647465;
									if(1.253773 <= x < 1.414214 ) weight=0.607212;
									if(1.414214 <= x < 1.595185 ) weight=0.589750;
									if(1.595185 <= x < 1.799315 ) weight=0.587406;
									if(1.799315 <= x < 2.029567 ) weight=0.592858;
									if(2.029567 <= x < 2.289283 ) weight=0.601059;
									if(2.289283 <= x < 2.582235 ) weight=0.608003;
									if(2.582235 <= x < 2.912674 ) weight=0.611705;
									if(2.912674 <= x < 3.285398 ) weight=0.610086;
									if(3.285398 <= x < 3.705818 ) weight=0.605015;
									if(3.705818 <= x < 4.180038 ) weight=0.596299;
									if(4.180038 <= x < 4.714942 ) weight=0.590727;
									if(4.714942 <= x < 5.318296 ) weight=0.585358;
									if(5.318296 <= x < 5.998859 ) weight=0.585257;
									if(5.998859 <= x < 6.766511 ) weight=0.580812;
									if(6.766511 <= x < 7.632396 ) weight=0.576207;
									if(7.632396 <= x < 8.609086 ) weight=0.575912;
									if(8.609086 <= x < 9.710759 ) weight=0.558718;
									if(9.710759 <= x < 10.953409 ) weight=0.555625;
									if(10.953409 <= x < 12.355077 ) weight=0.558886;
									if(12.355077 <= x < 13.936111 ) weight=0.545318;
									if(13.936111 <= x < 15.719464 ) weight=0.517607;
									if(15.719464 <= x < 17.731026 ) weight=0.512366;
									if(17.731026 <= x < 20.000000 ) weight=0.497034;
									
									
								}
								
								
							}//end of IsMCefix
							 //----------------------------------------------------------------------------
							 //correction based on data only for pi0
							if(!IsMCefix){
								if(TMath::Abs(fMCparticleMother->GetPdgCode())==111){
									Double_t x=mPt;
									if(0.100000 <= x < 0.112797 ) weight=1.262120;
									if(0.112797 <= x < 0.127231 ) weight=1.277765;
									if(0.127231 <= x < 0.143512 ) weight=1.295605;
									if(0.143512 <= x < 0.161877 ) weight=1.318155;
									if(0.161877 <= x < 0.182592 ) weight=1.348693;
									if(0.182592 <= x < 0.205957 ) weight=1.388636;
									if(0.205957 <= x < 0.232313 ) weight=1.439122;
									if(0.232313 <= x < 0.262041 ) weight=1.497452;
									if(0.262041 <= x < 0.295573 ) weight=1.559409;
									if(0.295573 <= x < 0.333397 ) weight=1.615169;
									if(0.333397 <= x < 0.376060 ) weight=1.654954;
									if(0.376060 <= x < 0.424183 ) weight=1.668753;
									if(0.424183 <= x < 0.478465 ) weight=1.652225;
									if(0.478465 <= x < 0.539692 ) weight=1.603119;
									if(0.539692 <= x < 0.608754 ) weight=1.526049;
									if(0.608754 <= x < 0.686654 ) weight=1.426724;
									if(0.686654 <= x < 0.774523 ) weight=1.312684;
									if(0.774523 <= x < 0.873636 ) weight=1.195395;
									if(0.873636 <= x < 0.985432 ) weight=1.086264;
									if(0.985432 <= x < 1.111534 ) weight=0.993666;
									if(1.111534 <= x < 1.253773 ) weight=0.922587;
									if(1.253773 <= x < 1.414214 ) weight=0.875739;
									if(1.414214 <= x < 1.595185 ) weight=0.852181;
									if(1.595185 <= x < 1.799315 ) weight=0.847828;
									if(1.799315 <= x < 2.029567 ) weight=0.863875;
									if(2.029567 <= x < 2.289283 ) weight=0.899112;
									if(2.289283 <= x < 2.582235 ) weight=0.955194;
									if(2.582235 <= x < 2.912674 ) weight=1.033824;
									if(2.912674 <= x < 3.285398 ) weight=1.133714;
									if(3.285398 <= x < 3.705818 ) weight=1.259471;
									if(3.705818 <= x < 4.180038 ) weight=1.406883;
									if(4.180038 <= x < 4.714942 ) weight=1.578923;
									if(4.714942 <= x < 5.318296 ) weight=1.778513;
									if(5.318296 <= x < 5.998859 ) weight=2.001171;
									if(5.998859 <= x < 6.766511 ) weight=2.223161;
									if(6.766511 <= x < 7.632396 ) weight=2.449445;
									if(7.632396 <= x < 8.609086 ) weight=2.661734;
									if(8.609086 <= x < 9.710759 ) weight=2.851935;
									if(9.710759 <= x < 10.953409 ) weight=2.974319;
									if(10.953409 <= x < 12.355077 ) weight=3.106314;
									if(12.355077 <= x < 13.936111 ) weight=3.126815;
									if(13.936111 <= x < 15.719464 ) weight=3.150053;
									if(15.719464 <= x < 17.731026 ) weight=3.218509;
									if(17.731026 <= x < 20.000000 ) weight=3.252141;
									
									
								}
								
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
							if(IsMCefix){
								if(TMath::Abs(fMCparticleGMother->GetPdgCode())==111){
									Double_t x=gmPt;
									if(0.100000 <= x < 0.112797 ) weight=1.030419;
									if(0.112797 <= x < 0.127231 ) weight=1.044554;
									if(0.127231 <= x < 0.143512 ) weight=1.062733;
									if(0.143512 <= x < 0.161877 ) weight=1.085332;
									if(0.161877 <= x < 0.182592 ) weight=1.115248;
									if(0.182592 <= x < 0.205957 ) weight=1.153990;
									if(0.205957 <= x < 0.232313 ) weight=1.201346;
									if(0.232313 <= x < 0.262041 ) weight=1.257332;
									if(0.262041 <= x < 0.295573 ) weight=1.315488;
									if(0.295573 <= x < 0.333397 ) weight=1.369138;
									if(0.333397 <= x < 0.376060 ) weight=1.407632;
									if(0.376060 <= x < 0.424183 ) weight=1.422232;
									if(0.424183 <= x < 0.478465 ) weight=1.406922;
									if(0.478465 <= x < 0.539692 ) weight=1.360082;
									if(0.539692 <= x < 0.608754 ) weight=1.284405;
									if(0.608754 <= x < 0.686654 ) weight=1.182017;
									if(0.686654 <= x < 0.774523 ) weight=1.062002;
									if(0.774523 <= x < 0.873636 ) weight=0.935533;
									if(0.873636 <= x < 0.985432 ) weight=0.816081;
									if(0.985432 <= x < 1.111534 ) weight=0.717527;
									if(1.111534 <= x < 1.253773 ) weight=0.647465;
									if(1.253773 <= x < 1.414214 ) weight=0.607212;
									if(1.414214 <= x < 1.595185 ) weight=0.589750;
									if(1.595185 <= x < 1.799315 ) weight=0.587406;
									if(1.799315 <= x < 2.029567 ) weight=0.592858;
									if(2.029567 <= x < 2.289283 ) weight=0.601059;
									if(2.289283 <= x < 2.582235 ) weight=0.608003;
									if(2.582235 <= x < 2.912674 ) weight=0.611705;
									if(2.912674 <= x < 3.285398 ) weight=0.610086;
									if(3.285398 <= x < 3.705818 ) weight=0.605015;
									if(3.705818 <= x < 4.180038 ) weight=0.596299;
									if(4.180038 <= x < 4.714942 ) weight=0.590727;
									if(4.714942 <= x < 5.318296 ) weight=0.585358;
									if(5.318296 <= x < 5.998859 ) weight=0.585257;
									if(5.998859 <= x < 6.766511 ) weight=0.580812;
									if(6.766511 <= x < 7.632396 ) weight=0.576207;
									if(7.632396 <= x < 8.609086 ) weight=0.575912;
									if(8.609086 <= x < 9.710759 ) weight=0.558718;
									if(9.710759 <= x < 10.953409 ) weight=0.555625;
									if(10.953409 <= x < 12.355077 ) weight=0.558886;
									if(12.355077 <= x < 13.936111 ) weight=0.545318;
									if(13.936111 <= x < 15.719464 ) weight=0.517607;
									if(15.719464 <= x < 17.731026 ) weight=0.512366;
									if(17.731026 <= x < 20.000000 ) weight=0.497034;
									
									
								}
								
								
							}//end of IsMCefix
							 //----------------------------------------------------------------------------
							
								//correction based on data only for pi0
							if(!IsMCefix){
								if(TMath::Abs(fMCparticleGMother->GetPdgCode())==111){
									Double_t x=gmPt;
									if(0.100000 <= x < 0.112797 ) weight=1.262120;
									if(0.112797 <= x < 0.127231 ) weight=1.277765;
									if(0.127231 <= x < 0.143512 ) weight=1.295605;
									if(0.143512 <= x < 0.161877 ) weight=1.318155;
									if(0.161877 <= x < 0.182592 ) weight=1.348693;
									if(0.182592 <= x < 0.205957 ) weight=1.388636;
									if(0.205957 <= x < 0.232313 ) weight=1.439122;
									if(0.232313 <= x < 0.262041 ) weight=1.497452;
									if(0.262041 <= x < 0.295573 ) weight=1.559409;
									if(0.295573 <= x < 0.333397 ) weight=1.615169;
									if(0.333397 <= x < 0.376060 ) weight=1.654954;
									if(0.376060 <= x < 0.424183 ) weight=1.668753;
									if(0.424183 <= x < 0.478465 ) weight=1.652225;
									if(0.478465 <= x < 0.539692 ) weight=1.603119;
									if(0.539692 <= x < 0.608754 ) weight=1.526049;
									if(0.608754 <= x < 0.686654 ) weight=1.426724;
									if(0.686654 <= x < 0.774523 ) weight=1.312684;
									if(0.774523 <= x < 0.873636 ) weight=1.195395;
									if(0.873636 <= x < 0.985432 ) weight=1.086264;
									if(0.985432 <= x < 1.111534 ) weight=0.993666;
									if(1.111534 <= x < 1.253773 ) weight=0.922587;
									if(1.253773 <= x < 1.414214 ) weight=0.875739;
									if(1.414214 <= x < 1.595185 ) weight=0.852181;
									if(1.595185 <= x < 1.799315 ) weight=0.847828;
									if(1.799315 <= x < 2.029567 ) weight=0.863875;
									if(2.029567 <= x < 2.289283 ) weight=0.899112;
									if(2.289283 <= x < 2.582235 ) weight=0.955194;
									if(2.582235 <= x < 2.912674 ) weight=1.033824;
									if(2.912674 <= x < 3.285398 ) weight=1.133714;
									if(3.285398 <= x < 3.705818 ) weight=1.259471;
									if(3.705818 <= x < 4.180038 ) weight=1.406883;
									if(4.180038 <= x < 4.714942 ) weight=1.578923;
									if(4.714942 <= x < 5.318296 ) weight=1.778513;
									if(5.318296 <= x < 5.998859 ) weight=2.001171;
									if(5.998859 <= x < 6.766511 ) weight=2.223161;
									if(6.766511 <= x < 7.632396 ) weight=2.449445;
									if(7.632396 <= x < 8.609086 ) weight=2.661734;
									if(8.609086 <= x < 9.710759 ) weight=2.851935;
									if(9.710759 <= x < 10.953409 ) weight=2.974319;
									if(10.953409 <= x < 12.355077 ) weight=3.106314;
									if(12.355077 <= x < 13.936111 ) weight=3.126815;
									if(13.936111 <= x < 15.719464 ) weight=3.150053;
									if(15.719464 <= x < 17.731026 ) weight=3.218509;
									if(17.731026 <= x < 20.000000 ) weight=3.252141;
									
									
								}
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
							if(IsMCefix){
								if(TMath::Abs(fMCparticleMother->GetPdgCode())==111){
									Double_t x=mPt;
									if(0.100000 <= x < 0.112797 ) weight=1.030419;
									if(0.112797 <= x < 0.127231 ) weight=1.044554;
									if(0.127231 <= x < 0.143512 ) weight=1.062733;
									if(0.143512 <= x < 0.161877 ) weight=1.085332;
									if(0.161877 <= x < 0.182592 ) weight=1.115248;
									if(0.182592 <= x < 0.205957 ) weight=1.153990;
									if(0.205957 <= x < 0.232313 ) weight=1.201346;
									if(0.232313 <= x < 0.262041 ) weight=1.257332;
									if(0.262041 <= x < 0.295573 ) weight=1.315488;
									if(0.295573 <= x < 0.333397 ) weight=1.369138;
									if(0.333397 <= x < 0.376060 ) weight=1.407632;
									if(0.376060 <= x < 0.424183 ) weight=1.422232;
									if(0.424183 <= x < 0.478465 ) weight=1.406922;
									if(0.478465 <= x < 0.539692 ) weight=1.360082;
									if(0.539692 <= x < 0.608754 ) weight=1.284405;
									if(0.608754 <= x < 0.686654 ) weight=1.182017;
									if(0.686654 <= x < 0.774523 ) weight=1.062002;
									if(0.774523 <= x < 0.873636 ) weight=0.935533;
									if(0.873636 <= x < 0.985432 ) weight=0.816081;
									if(0.985432 <= x < 1.111534 ) weight=0.717527;
									if(1.111534 <= x < 1.253773 ) weight=0.647465;
									if(1.253773 <= x < 1.414214 ) weight=0.607212;
									if(1.414214 <= x < 1.595185 ) weight=0.589750;
									if(1.595185 <= x < 1.799315 ) weight=0.587406;
									if(1.799315 <= x < 2.029567 ) weight=0.592858;
									if(2.029567 <= x < 2.289283 ) weight=0.601059;
									if(2.289283 <= x < 2.582235 ) weight=0.608003;
									if(2.582235 <= x < 2.912674 ) weight=0.611705;
									if(2.912674 <= x < 3.285398 ) weight=0.610086;
									if(3.285398 <= x < 3.705818 ) weight=0.605015;
									if(3.705818 <= x < 4.180038 ) weight=0.596299;
									if(4.180038 <= x < 4.714942 ) weight=0.590727;
									if(4.714942 <= x < 5.318296 ) weight=0.585358;
									if(5.318296 <= x < 5.998859 ) weight=0.585257;
									if(5.998859 <= x < 6.766511 ) weight=0.580812;
									if(6.766511 <= x < 7.632396 ) weight=0.576207;
									if(7.632396 <= x < 8.609086 ) weight=0.575912;
									if(8.609086 <= x < 9.710759 ) weight=0.558718;
									if(9.710759 <= x < 10.953409 ) weight=0.555625;
									if(10.953409 <= x < 12.355077 ) weight=0.558886;
									if(12.355077 <= x < 13.936111 ) weight=0.545318;
									if(13.936111 <= x < 15.719464 ) weight=0.517607;
									if(15.719464 <= x < 17.731026 ) weight=0.512366;
									if(17.731026 <= x < 20.000000 ) weight=0.497034;
									
									
								}
								
								
							}//end of IsMCefix
							 //----------------------------------------------------------------------------
							
								//correction based on data only for pi0 for d3
							if(!IsMCefix){
								if(TMath::Abs(fMCparticleMother->GetPdgCode())==111){
									Double_t x=mPt;
									if(0.100000 <= x < 0.112797 ) weight=1.262120;
									if(0.112797 <= x < 0.127231 ) weight=1.277765;
									if(0.127231 <= x < 0.143512 ) weight=1.295605;
									if(0.143512 <= x < 0.161877 ) weight=1.318155;
									if(0.161877 <= x < 0.182592 ) weight=1.348693;
									if(0.182592 <= x < 0.205957 ) weight=1.388636;
									if(0.205957 <= x < 0.232313 ) weight=1.439122;
									if(0.232313 <= x < 0.262041 ) weight=1.497452;
									if(0.262041 <= x < 0.295573 ) weight=1.559409;
									if(0.295573 <= x < 0.333397 ) weight=1.615169;
									if(0.333397 <= x < 0.376060 ) weight=1.654954;
									if(0.376060 <= x < 0.424183 ) weight=1.668753;
									if(0.424183 <= x < 0.478465 ) weight=1.652225;
									if(0.478465 <= x < 0.539692 ) weight=1.603119;
									if(0.539692 <= x < 0.608754 ) weight=1.526049;
									if(0.608754 <= x < 0.686654 ) weight=1.426724;
									if(0.686654 <= x < 0.774523 ) weight=1.312684;
									if(0.774523 <= x < 0.873636 ) weight=1.195395;
									if(0.873636 <= x < 0.985432 ) weight=1.086264;
									if(0.985432 <= x < 1.111534 ) weight=0.993666;
									if(1.111534 <= x < 1.253773 ) weight=0.922587;
									if(1.253773 <= x < 1.414214 ) weight=0.875739;
									if(1.414214 <= x < 1.595185 ) weight=0.852181;
									if(1.595185 <= x < 1.799315 ) weight=0.847828;
									if(1.799315 <= x < 2.029567 ) weight=0.863875;
									if(2.029567 <= x < 2.289283 ) weight=0.899112;
									if(2.289283 <= x < 2.582235 ) weight=0.955194;
									if(2.582235 <= x < 2.912674 ) weight=1.033824;
									if(2.912674 <= x < 3.285398 ) weight=1.133714;
									if(3.285398 <= x < 3.705818 ) weight=1.259471;
									if(3.705818 <= x < 4.180038 ) weight=1.406883;
									if(4.180038 <= x < 4.714942 ) weight=1.578923;
									if(4.714942 <= x < 5.318296 ) weight=1.778513;
									if(5.318296 <= x < 5.998859 ) weight=2.001171;
									if(5.998859 <= x < 6.766511 ) weight=2.223161;
									if(6.766511 <= x < 7.632396 ) weight=2.449445;
									if(7.632396 <= x < 8.609086 ) weight=2.661734;
									if(8.609086 <= x < 9.710759 ) weight=2.851935;
									if(9.710759 <= x < 10.953409 ) weight=2.974319;
									if(10.953409 <= x < 12.355077 ) weight=3.106314;
									if(12.355077 <= x < 13.936111 ) weight=3.126815;
									if(13.936111 <= x < 15.719464 ) weight=3.150053;
									if(15.719464 <= x < 17.731026 ) weight=3.218509;
									if(17.731026 <= x < 20.000000 ) weight=3.252141;
									
									
								}
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
							if(IsMCefix){
								if(TMath::Abs(fMCparticleGMother->GetPdgCode())==111){
									Double_t x=gmPt;
									if(0.100000 <= x < 0.112797 ) weight=1.030419;
									if(0.112797 <= x < 0.127231 ) weight=1.044554;
									if(0.127231 <= x < 0.143512 ) weight=1.062733;
									if(0.143512 <= x < 0.161877 ) weight=1.085332;
									if(0.161877 <= x < 0.182592 ) weight=1.115248;
									if(0.182592 <= x < 0.205957 ) weight=1.153990;
									if(0.205957 <= x < 0.232313 ) weight=1.201346;
									if(0.232313 <= x < 0.262041 ) weight=1.257332;
									if(0.262041 <= x < 0.295573 ) weight=1.315488;
									if(0.295573 <= x < 0.333397 ) weight=1.369138;
									if(0.333397 <= x < 0.376060 ) weight=1.407632;
									if(0.376060 <= x < 0.424183 ) weight=1.422232;
									if(0.424183 <= x < 0.478465 ) weight=1.406922;
									if(0.478465 <= x < 0.539692 ) weight=1.360082;
									if(0.539692 <= x < 0.608754 ) weight=1.284405;
									if(0.608754 <= x < 0.686654 ) weight=1.182017;
									if(0.686654 <= x < 0.774523 ) weight=1.062002;
									if(0.774523 <= x < 0.873636 ) weight=0.935533;
									if(0.873636 <= x < 0.985432 ) weight=0.816081;
									if(0.985432 <= x < 1.111534 ) weight=0.717527;
									if(1.111534 <= x < 1.253773 ) weight=0.647465;
									if(1.253773 <= x < 1.414214 ) weight=0.607212;
									if(1.414214 <= x < 1.595185 ) weight=0.589750;
									if(1.595185 <= x < 1.799315 ) weight=0.587406;
									if(1.799315 <= x < 2.029567 ) weight=0.592858;
									if(2.029567 <= x < 2.289283 ) weight=0.601059;
									if(2.289283 <= x < 2.582235 ) weight=0.608003;
									if(2.582235 <= x < 2.912674 ) weight=0.611705;
									if(2.912674 <= x < 3.285398 ) weight=0.610086;
									if(3.285398 <= x < 3.705818 ) weight=0.605015;
									if(3.705818 <= x < 4.180038 ) weight=0.596299;
									if(4.180038 <= x < 4.714942 ) weight=0.590727;
									if(4.714942 <= x < 5.318296 ) weight=0.585358;
									if(5.318296 <= x < 5.998859 ) weight=0.585257;
									if(5.998859 <= x < 6.766511 ) weight=0.580812;
									if(6.766511 <= x < 7.632396 ) weight=0.576207;
									if(7.632396 <= x < 8.609086 ) weight=0.575912;
									if(8.609086 <= x < 9.710759 ) weight=0.558718;
									if(9.710759 <= x < 10.953409 ) weight=0.555625;
									if(10.953409 <= x < 12.355077 ) weight=0.558886;
									if(12.355077 <= x < 13.936111 ) weight=0.545318;
									if(13.936111 <= x < 15.719464 ) weight=0.517607;
									if(15.719464 <= x < 17.731026 ) weight=0.512366;
									if(17.731026 <= x < 20.000000 ) weight=0.497034;
									
									
								}
								
								
							}//end of IsMCefix
							 //----------------------------------------------------------------------------
							 //correction based on data only for pi0
							if(!IsMCefix){
								if(TMath::Abs(fMCparticleGMother->GetPdgCode())==111){
									Double_t x=gmPt;
									if(0.100000 <= x < 0.112797 ) weight=1.262120;
									if(0.112797 <= x < 0.127231 ) weight=1.277765;
									if(0.127231 <= x < 0.143512 ) weight=1.295605;
									if(0.143512 <= x < 0.161877 ) weight=1.318155;
									if(0.161877 <= x < 0.182592 ) weight=1.348693;
									if(0.182592 <= x < 0.205957 ) weight=1.388636;
									if(0.205957 <= x < 0.232313 ) weight=1.439122;
									if(0.232313 <= x < 0.262041 ) weight=1.497452;
									if(0.262041 <= x < 0.295573 ) weight=1.559409;
									if(0.295573 <= x < 0.333397 ) weight=1.615169;
									if(0.333397 <= x < 0.376060 ) weight=1.654954;
									if(0.376060 <= x < 0.424183 ) weight=1.668753;
									if(0.424183 <= x < 0.478465 ) weight=1.652225;
									if(0.478465 <= x < 0.539692 ) weight=1.603119;
									if(0.539692 <= x < 0.608754 ) weight=1.526049;
									if(0.608754 <= x < 0.686654 ) weight=1.426724;
									if(0.686654 <= x < 0.774523 ) weight=1.312684;
									if(0.774523 <= x < 0.873636 ) weight=1.195395;
									if(0.873636 <= x < 0.985432 ) weight=1.086264;
									if(0.985432 <= x < 1.111534 ) weight=0.993666;
									if(1.111534 <= x < 1.253773 ) weight=0.922587;
									if(1.253773 <= x < 1.414214 ) weight=0.875739;
									if(1.414214 <= x < 1.595185 ) weight=0.852181;
									if(1.595185 <= x < 1.799315 ) weight=0.847828;
									if(1.799315 <= x < 2.029567 ) weight=0.863875;
									if(2.029567 <= x < 2.289283 ) weight=0.899112;
									if(2.289283 <= x < 2.582235 ) weight=0.955194;
									if(2.582235 <= x < 2.912674 ) weight=1.033824;
									if(2.912674 <= x < 3.285398 ) weight=1.133714;
									if(3.285398 <= x < 3.705818 ) weight=1.259471;
									if(3.705818 <= x < 4.180038 ) weight=1.406883;
									if(4.180038 <= x < 4.714942 ) weight=1.578923;
									if(4.714942 <= x < 5.318296 ) weight=1.778513;
									if(5.318296 <= x < 5.998859 ) weight=2.001171;
									if(5.998859 <= x < 6.766511 ) weight=2.223161;
									if(6.766511 <= x < 7.632396 ) weight=2.449445;
									if(7.632396 <= x < 8.609086 ) weight=2.661734;
									if(8.609086 <= x < 9.710759 ) weight=2.851935;
									if(9.710759 <= x < 10.953409 ) weight=2.974319;
									if(10.953409 <= x < 12.355077 ) weight=3.106314;
									if(12.355077 <= x < 13.936111 ) weight=3.126815;
									if(13.936111 <= x < 15.719464 ) weight=3.150053;
									if(15.719464 <= x < 17.731026 ) weight=3.218509;
									if(17.731026 <= x < 20.000000 ) weight=3.252141;
									
									
								}
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
						
							//ULS with no weight from ULS-LS original
							// we have to know if track2 comes from same mother!!!
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
										
										if(fMCparticle2->GetMother()==fMCparticle->GetMother()) fPtElec_ULS_MC->Fill(fPtE);
										
										//-----------------------------------------------------------------------------------------------------------
										//weight for mother
										Double_t weight2=1;
										Double_t mPt=fMCparticleMother->Pt();
										
																				
										if(TMath::Abs(fMCparticleMother->GetPdgCode())==111){
											Double_t x=mPt;
											
											if(!IsMCefix){
												if(0.100000 <= x < 0.112797 ) weight=1.262120;
												if(0.112797 <= x < 0.127231 ) weight=1.277765;
												if(0.127231 <= x < 0.143512 ) weight=1.295605;
												if(0.143512 <= x < 0.161877 ) weight=1.318155;
												if(0.161877 <= x < 0.182592 ) weight=1.348693;
												if(0.182592 <= x < 0.205957 ) weight=1.388636;
												if(0.205957 <= x < 0.232313 ) weight=1.439122;
												if(0.232313 <= x < 0.262041 ) weight=1.497452;
												if(0.262041 <= x < 0.295573 ) weight=1.559409;
												if(0.295573 <= x < 0.333397 ) weight=1.615169;
												if(0.333397 <= x < 0.376060 ) weight=1.654954;
												if(0.376060 <= x < 0.424183 ) weight=1.668753;
												if(0.424183 <= x < 0.478465 ) weight=1.652225;
												if(0.478465 <= x < 0.539692 ) weight=1.603119;
												if(0.539692 <= x < 0.608754 ) weight=1.526049;
												if(0.608754 <= x < 0.686654 ) weight=1.426724;
												if(0.686654 <= x < 0.774523 ) weight=1.312684;
												if(0.774523 <= x < 0.873636 ) weight=1.195395;
												if(0.873636 <= x < 0.985432 ) weight=1.086264;
												if(0.985432 <= x < 1.111534 ) weight=0.993666;
												if(1.111534 <= x < 1.253773 ) weight=0.922587;
												if(1.253773 <= x < 1.414214 ) weight=0.875739;
												if(1.414214 <= x < 1.595185 ) weight=0.852181;
												if(1.595185 <= x < 1.799315 ) weight=0.847828;
												if(1.799315 <= x < 2.029567 ) weight=0.863875;
												if(2.029567 <= x < 2.289283 ) weight=0.899112;
												if(2.289283 <= x < 2.582235 ) weight=0.955194;
												if(2.582235 <= x < 2.912674 ) weight=1.033824;
												if(2.912674 <= x < 3.285398 ) weight=1.133714;
												if(3.285398 <= x < 3.705818 ) weight=1.259471;
												if(3.705818 <= x < 4.180038 ) weight=1.406883;
												if(4.180038 <= x < 4.714942 ) weight=1.578923;
												if(4.714942 <= x < 5.318296 ) weight=1.778513;
												if(5.318296 <= x < 5.998859 ) weight=2.001171;
												if(5.998859 <= x < 6.766511 ) weight=2.223161;
												if(6.766511 <= x < 7.632396 ) weight=2.449445;
												if(7.632396 <= x < 8.609086 ) weight=2.661734;
												if(8.609086 <= x < 9.710759 ) weight=2.851935;
												if(9.710759 <= x < 10.953409 ) weight=2.974319;
												if(10.953409 <= x < 12.355077 ) weight=3.106314;
												if(12.355077 <= x < 13.936111 ) weight=3.126815;
												if(13.936111 <= x < 15.719464 ) weight=3.150053;
												if(15.719464 <= x < 17.731026 ) weight=3.218509;
												if(17.731026 <= x < 20.000000 ) weight=3.252141;
												
											}
											if(IsMCefix){
												
												if(0.100000 <= x < 0.112797 ) weight2=1.030419;
												if(0.112797 <= x < 0.127231 ) weight2=1.044554;
												if(0.127231 <= x < 0.143512 ) weight2=1.062733;
												if(0.143512 <= x < 0.161877 ) weight2=1.085332;
												if(0.161877 <= x < 0.182592 ) weight2=1.115248;
												if(0.182592 <= x < 0.205957 ) weight2=1.153990;
												if(0.205957 <= x < 0.232313 ) weight2=1.201346;
												if(0.232313 <= x < 0.262041 ) weight2=1.257332;
												if(0.262041 <= x < 0.295573 ) weight2=1.315488;
												if(0.295573 <= x < 0.333397 ) weight2=1.369138;
												if(0.333397 <= x < 0.376060 ) weight2=1.407632;
												if(0.376060 <= x < 0.424183 ) weight2=1.422232;
												if(0.424183 <= x < 0.478465 ) weight2=1.406922;
												if(0.478465 <= x < 0.539692 ) weight2=1.360082;
												if(0.539692 <= x < 0.608754 ) weight2=1.284405;
												if(0.608754 <= x < 0.686654 ) weight2=1.182017;
												if(0.686654 <= x < 0.774523 ) weight2=1.062002;
												if(0.774523 <= x < 0.873636 ) weight2=0.935533;
												if(0.873636 <= x < 0.985432 ) weight2=0.816081;
												if(0.985432 <= x < 1.111534 ) weight2=0.717527;
												if(1.111534 <= x < 1.253773 ) weight2=0.647465;
												if(1.253773 <= x < 1.414214 ) weight2=0.607212;
												if(1.414214 <= x < 1.595185 ) weight2=0.589750;
												if(1.595185 <= x < 1.799315 ) weight2=0.587406;
												if(1.799315 <= x < 2.029567 ) weight2=0.592858;
												if(2.029567 <= x < 2.289283 ) weight2=0.601059;
												if(2.289283 <= x < 2.582235 ) weight2=0.608003;
												if(2.582235 <= x < 2.912674 ) weight2=0.611705;
												if(2.912674 <= x < 3.285398 ) weight2=0.610086;
												if(3.285398 <= x < 3.705818 ) weight2=0.605015;
												if(3.705818 <= x < 4.180038 ) weight2=0.596299;
												if(4.180038 <= x < 4.714942 ) weight2=0.590727;
												if(4.714942 <= x < 5.318296 ) weight2=0.585358;
												if(5.318296 <= x < 5.998859 ) weight2=0.585257;
												if(5.998859 <= x < 6.766511 ) weight2=0.580812;
												if(6.766511 <= x < 7.632396 ) weight2=0.576207;
												if(7.632396 <= x < 8.609086 ) weight2=0.575912;
												if(8.609086 <= x < 9.710759 ) weight2=0.558718;
												if(9.710759 <= x < 10.953409 ) weight2=0.555625;
												if(10.953409 <= x < 12.355077 ) weight2=0.558886;
												if(12.355077 <= x < 13.936111 ) weight2=0.545318;
												if(13.936111 <= x < 15.719464 ) weight2=0.517607;
												if(15.719464 <= x < 17.731026 ) weight2=0.512366;
												if(17.731026 <= x < 20.000000 ) weight2=0.497034;
												
												
											}
											
										}
											//weight for grandmother
										Double_t gmPt=fMCparticleGMother->Pt();
										if(TMath::Abs((fMCparticleMother->GetMother()>0) && ((fMCparticleGMother->GetPdgCode())==111))){
											Double_t x=gmPt;
											
											if(!IsMCefix){
												if(0.100000 <= x < 0.112797 ) weight=1.262120;
												if(0.112797 <= x < 0.127231 ) weight=1.277765;
												if(0.127231 <= x < 0.143512 ) weight=1.295605;
												if(0.143512 <= x < 0.161877 ) weight=1.318155;
												if(0.161877 <= x < 0.182592 ) weight=1.348693;
												if(0.182592 <= x < 0.205957 ) weight=1.388636;
												if(0.205957 <= x < 0.232313 ) weight=1.439122;
												if(0.232313 <= x < 0.262041 ) weight=1.497452;
												if(0.262041 <= x < 0.295573 ) weight=1.559409;
												if(0.295573 <= x < 0.333397 ) weight=1.615169;
												if(0.333397 <= x < 0.376060 ) weight=1.654954;
												if(0.376060 <= x < 0.424183 ) weight=1.668753;
												if(0.424183 <= x < 0.478465 ) weight=1.652225;
												if(0.478465 <= x < 0.539692 ) weight=1.603119;
												if(0.539692 <= x < 0.608754 ) weight=1.526049;
												if(0.608754 <= x < 0.686654 ) weight=1.426724;
												if(0.686654 <= x < 0.774523 ) weight=1.312684;
												if(0.774523 <= x < 0.873636 ) weight=1.195395;
												if(0.873636 <= x < 0.985432 ) weight=1.086264;
												if(0.985432 <= x < 1.111534 ) weight=0.993666;
												if(1.111534 <= x < 1.253773 ) weight=0.922587;
												if(1.253773 <= x < 1.414214 ) weight=0.875739;
												if(1.414214 <= x < 1.595185 ) weight=0.852181;
												if(1.595185 <= x < 1.799315 ) weight=0.847828;
												if(1.799315 <= x < 2.029567 ) weight=0.863875;
												if(2.029567 <= x < 2.289283 ) weight=0.899112;
												if(2.289283 <= x < 2.582235 ) weight=0.955194;
												if(2.582235 <= x < 2.912674 ) weight=1.033824;
												if(2.912674 <= x < 3.285398 ) weight=1.133714;
												if(3.285398 <= x < 3.705818 ) weight=1.259471;
												if(3.705818 <= x < 4.180038 ) weight=1.406883;
												if(4.180038 <= x < 4.714942 ) weight=1.578923;
												if(4.714942 <= x < 5.318296 ) weight=1.778513;
												if(5.318296 <= x < 5.998859 ) weight=2.001171;
												if(5.998859 <= x < 6.766511 ) weight=2.223161;
												if(6.766511 <= x < 7.632396 ) weight=2.449445;
												if(7.632396 <= x < 8.609086 ) weight=2.661734;
												if(8.609086 <= x < 9.710759 ) weight=2.851935;
												if(9.710759 <= x < 10.953409 ) weight=2.974319;
												if(10.953409 <= x < 12.355077 ) weight=3.106314;
												if(12.355077 <= x < 13.936111 ) weight=3.126815;
												if(13.936111 <= x < 15.719464 ) weight=3.150053;
												if(15.719464 <= x < 17.731026 ) weight=3.218509;
												if(17.731026 <= x < 20.000000 ) weight=3.252141;
												
											}
											if(IsMCefix){
												
												if(0.100000 <= x < 0.112797 ) weight2=1.030419;
												if(0.112797 <= x < 0.127231 ) weight2=1.044554;
												if(0.127231 <= x < 0.143512 ) weight2=1.062733;
												if(0.143512 <= x < 0.161877 ) weight2=1.085332;
												if(0.161877 <= x < 0.182592 ) weight2=1.115248;
												if(0.182592 <= x < 0.205957 ) weight2=1.153990;
												if(0.205957 <= x < 0.232313 ) weight2=1.201346;
												if(0.232313 <= x < 0.262041 ) weight2=1.257332;
												if(0.262041 <= x < 0.295573 ) weight2=1.315488;
												if(0.295573 <= x < 0.333397 ) weight2=1.369138;
												if(0.333397 <= x < 0.376060 ) weight2=1.407632;
												if(0.376060 <= x < 0.424183 ) weight2=1.422232;
												if(0.424183 <= x < 0.478465 ) weight2=1.406922;
												if(0.478465 <= x < 0.539692 ) weight2=1.360082;
												if(0.539692 <= x < 0.608754 ) weight2=1.284405;
												if(0.608754 <= x < 0.686654 ) weight2=1.182017;
												if(0.686654 <= x < 0.774523 ) weight2=1.062002;
												if(0.774523 <= x < 0.873636 ) weight2=0.935533;
												if(0.873636 <= x < 0.985432 ) weight2=0.816081;
												if(0.985432 <= x < 1.111534 ) weight2=0.717527;
												if(1.111534 <= x < 1.253773 ) weight2=0.647465;
												if(1.253773 <= x < 1.414214 ) weight2=0.607212;
												if(1.414214 <= x < 1.595185 ) weight2=0.589750;
												if(1.595185 <= x < 1.799315 ) weight2=0.587406;
												if(1.799315 <= x < 2.029567 ) weight2=0.592858;
												if(2.029567 <= x < 2.289283 ) weight2=0.601059;
												if(2.289283 <= x < 2.582235 ) weight2=0.608003;
												if(2.582235 <= x < 2.912674 ) weight2=0.611705;
												if(2.912674 <= x < 3.285398 ) weight2=0.610086;
												if(3.285398 <= x < 3.705818 ) weight2=0.605015;
												if(3.705818 <= x < 4.180038 ) weight2=0.596299;
												if(4.180038 <= x < 4.714942 ) weight2=0.590727;
												if(4.714942 <= x < 5.318296 ) weight2=0.585358;
												if(5.318296 <= x < 5.998859 ) weight2=0.585257;
												if(5.998859 <= x < 6.766511 ) weight2=0.580812;
												if(6.766511 <= x < 7.632396 ) weight2=0.576207;
												if(7.632396 <= x < 8.609086 ) weight2=0.575912;
												if(8.609086 <= x < 9.710759 ) weight2=0.558718;
												if(9.710759 <= x < 10.953409 ) weight2=0.555625;
												if(10.953409 <= x < 12.355077 ) weight2=0.558886;
												if(12.355077 <= x < 13.936111 ) weight2=0.545318;
												if(13.936111 <= x < 15.719464 ) weight2=0.517607;
												if(15.719464 <= x < 17.731026 ) weight2=0.512366;
												if(17.731026 <= x < 20.000000 ) weight2=0.497034;
												
											}
										}
										
										if(fMCparticle2->GetMother()==fMCparticle->GetMother()) fPtElec_ULS_MC_weight->Fill(fPtE, 1./weight2);
										
											//-----------------------------------------------------------------------------------------------------------
											//end of weight
										
									}//partner found same as track
								}//loop in all partner
								
							}//track
						}//is ULS
						
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
		
		
		
}//end of Background function

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
