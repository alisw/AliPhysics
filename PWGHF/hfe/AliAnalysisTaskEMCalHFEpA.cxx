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
	//		version: September 12th, 2013.								  //
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
	//______________________________________________________________________

	//______________________________________________________________________
ClassImp(AliAnalysisTaskEMCalHFEpA)

	//______________________________________________________________________
AliAnalysisTaskEMCalHFEpA::AliAnalysisTaskEMCalHFEpA(const char *name) 
: AliAnalysisTaskSE(name)
,fCorrelationFlag(0)
,fIsMC(0)
,fUseEMCal(kFALSE)
,fUseShowerShapeCut(kFALSE)
,fFillBackground(kFALSE)


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

,fCharge_n(0)
,fCharge_p(0)

,fTime(0)
,fTime2(0)
,ftimingEle(0)
,ftimingEle2(0)


,fPtElec_ULS(0)
,fPtElec_LS(0)
,fPtElec_ULS2(0)
,fPtElec_LS2(0)
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
,fEtaPhi(0)
,fVtxZ(0)
,fNTracks(0)
,fNClusters(0)
,fTPCNcls_EoverP(0)
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
,fPtBackgroundBeforeReco(0)
,fPtBackgroundBeforeReco2(0)
,fPtBackgroundAfterReco(0)

,fPtMinAsso(0.3)
,fTpcNclsAsso(80)

,fPtMCparticleAll(0)
,fPtMCparticleAll_nonPrimary(0)
,fPtMCparticleAlle_nonPrimary(0)
,fPtMCparticleReco(0)
,fPtMCparticleReco_nonPrimary(0)
,fPtMCparticleAllHfe1(0)
,fPtMCparticleRecoHfe1(0)
,fPtMCparticleAllHfe2(0)
,fPtMCparticleRecoHfe2(0)
,fPtMCelectronAfterAll(0)
,fPtMCelectronAfterAll_nonPrimary(0)
,fPtMCpi0(0)
,fPtMCeta(0)
,fPtMC_EMCal_All(0)
,fPtMC_EMCal_Selected(0)
,fPtMC_TPC_All(0)
,fPtMC_TPC_Selected(0)
,fPtMCWithLabel(0)
,fPtMCWithoutLabel(0)
,fPtIsPhysicaPrimary(0)
,fCuts(0)
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
,fUseShowerShapeCut(kFALSE)
,fFillBackground(kFALSE)

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

,fCharge_n(0)
,fCharge_p(0)

,fTime(0)
,fTime2(0)
,ftimingEle(0)
,ftimingEle2(0)



,fPtElec_ULS(0)
,fPtElec_LS(0)
,fPtElec_ULS2(0)
,fPtElec_LS2(0)
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
,fEtaPhi(0)
,fVtxZ(0)
,fNTracks(0)
,fNClusters(0)
,fTPCNcls_EoverP(0)
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
,fPtBackgroundBeforeReco(0)
,fPtBackgroundBeforeReco2(0)
,fPtBackgroundAfterReco(0)

,fPtMinAsso(0.3)
,fTpcNclsAsso(80)


,fPtMCparticleAll(0)
,fPtMCparticleAll_nonPrimary(0)
,fPtMCparticleAlle_nonPrimary(0)
,fPtMCparticleReco(0)
,fPtMCparticleReco_nonPrimary(0)

,fPtMCparticleAllHfe1(0)
,fPtMCparticleRecoHfe1(0)
,fPtMCparticleAllHfe2(0)
,fPtMCparticleRecoHfe2(0)
,fPtMCelectronAfterAll(0)
,fPtMCelectronAfterAll_nonPrimary(0)

,fPtMCpi0(0)
,fPtMCeta(0)
,fPtMC_EMCal_All(0)
,fPtMC_EMCal_Selected(0)
,fPtMC_TPC_All(0)
,fPtMC_TPC_Selected(0)
,fPtMCWithLabel(0)
,fPtMCWithoutLabel(0)
,fPtIsPhysicaPrimary(0)
,fCuts(0)
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
	
		///______________________________________________________________________
		///Output Tlist
		//Create TList
	fOutputList = new TList();
	fOutputList->SetOwner();	
	
		//PIDqa
	fOutputList->Add(fPIDqa->MakeList("PIDQA"));
	
		//Store the number of events
		//Define the histo
	fNevent = new TH1F("fNevent","Number of Events",5,-0.5,4.5);
		//And then, add to the output list
	fOutputList->Add(fNevent);
	
	fpid = new TH1F("fpid","PID flag",5,0,5);
	fOutputList->Add(fpid);
	
		//pt Distribution
	fPtElec_Inc = new TH1F("fPtElec_Inc","Inclusive Electrons; p_{T} (GeV/c); Count",300,0,30);
	
	fPtElec_ULS = new TH1F("fPtElec_ULS","Inclusive Electrons; p_{T} (GeV/c); Count",300,0,30);
	fPtElec_LS = new TH1F("fPtElec_LS","Inclusive Electrons; p_{T} (GeV/c); Count",300,0,30);
	
	if(fFillBackground){
		fPtElec_ULS2 = new TH1F("fPtElec_ULS2","Inclusive Electrons; p_{T} (GeV/c); Count",300,0,30);
		fPtElec_LS2 = new TH1F("fPtElec_LS2","Inclusive Electrons; p_{T} (GeV/c); Count",300,0,30);
	}
	
	fPtTrigger_Inc = new TH1F("fPtTrigger_Inc","pT dist for Hadron Contamination; p_{t} (GeV/c); Count",300,0,30);
	fTPCnsigma_pt_2D = new TH2F("fTPCnsigma_pt_2D",";pt (GeV/c);TPC Electron N#sigma",1000,0.3,30,1000,-15,10);
	fShowerShapeCut = new TH2F("fShowerShapeCut","Shower Shape;M02;M20",500,0,1.8,500,0,1.8);
	
	
	
	fCharge_n = new TH1F("fCharge_n","Inclusive Electrons (Negative Charge); p_{t} (GeV/c); Count",200,0,30);
	fCharge_p = new TH1F("fCharge_p","Inclusive Positrons (Positive Charge); p_{t} (GeV/c); Count",200,0,30);
	
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
	
	fOutputList->Add(fPtElec_Inc);
	fOutputList->Add(fPtElec_ULS);
	fOutputList->Add(fPtElec_LS);
	
	if(fFillBackground){
		fOutputList->Add(fPtElec_ULS2);
		fOutputList->Add(fPtElec_LS2);
	}
	
	
	fOutputList->Add(fPtTrigger_Inc);
	fOutputList->Add(fTPCnsigma_pt_2D);
	fOutputList->Add(fShowerShapeCut);
	
	fOutputList->Add(fCharge_n);
	fOutputList->Add(fCharge_p);
	
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
	fNTracks= new  TH1F *[3];
	fNClusters= new TH1F *[3];
	fTPCNcls_EoverP= new TH2F *[3];	
	
	for(Int_t i = 0; i < 3; i++)
	{
		fEoverP_pt[i] = new TH2F(Form("fEoverP_pt%d",i),";p_{t} (GeV/c);E / p ",1000,0,30,500,0,2);
		fTPC_p[i] = new TH2F(Form("fTPC_p%d",i),";pt (GeV/c);TPC dE/dx (a. u.)",1000,0.3,15,1000,-20,200);
		fTPCnsigma_p[i] = new TH2F(Form("fTPCnsigma_p%d",i),";p (GeV/c);TPC Electron N#sigma",1000,0.3,15,1000,-15,10);
		
		
		fECluster[i]= new TH1F(Form("fECluster%d",i), ";ECluster",2000, 0,100);
		fEtaPhi[i]= new TH2F(Form("fEtaPhi%d",i),"#eta x #phi Clusters;#phi;#eta",200,0.,5,50,-1.,1.);
		fVtxZ[i]= new  TH1F(Form("fVtxZ%d",i),"VtxZ",1000, -50,50);
		fNTracks[i]= new  TH1F(Form("fNTracks%d",i),"NTracks",1000, 0,1000);
		fNClusters[i]= new TH1F(Form("fNClusters%d",i),"fNClusters0",200, 0,100);
		fTPCNcls_EoverP[i]= new TH2F(Form("fTPCNcls_EoverP%d",i),"TPCNcls_EoverP",1000,0,200,200,0,2);	
		
		
		fOutputList->Add(fEoverP_pt[i]);
		fOutputList->Add(fTPC_p[i]);
		fOutputList->Add(fTPCnsigma_p[i]);
		
		
		fOutputList->Add(fECluster[i]);
		fOutputList->Add(fEtaPhi[i]);
		fOutputList->Add(fVtxZ[i]);
		fOutputList->Add(fNTracks[i]);
		fOutputList->Add(fNClusters[i]);
		fOutputList->Add(fTPCNcls_EoverP[i]);
	}
	
		//pt bin
	Int_t fPtBin[6] = {2,4,6,8,10,15};
	
	fEoverP_tpc = new TH2F *[5];
	fTPC_pt = new TH1F *[5];
	fTPCnsigma_pt = new TH1F *[5];
	
	fEta=new TH1F *[5];
	fPhi=new TH1F *[5];
	fR=new TH1F *[5];
	fR_EoverP=new TH2F *[5];
	fNcells=new TH1F *[5];
	fNcells_EoverP=new TH2F *[5];
	fM02_EoverP= new TH2F *[5];
	fM20_EoverP= new TH2F *[5];
	fEoverP_ptbins=new TH1F *[5];
	fECluster_ptbins=new TH1F *[5];
	fEoverP_wSSCut=new TH1F *[5];
	fNcells_electrons=new TH1F *[5];
	fNcells_hadrons=new TH1F *[5];
	fTPCnsigma_eta_electrons=new TH2F *[5];
	fTPCnsigma_eta_hadrons=new TH2F *[5];
	
	if(fCorrelationFlag)
	{
		fCEtaPhi_Inc = new TH2F *[5];
		fCEtaPhi_Inc_DiHadron = new TH2F *[5];
		
		fCEtaPhi_ULS = new TH2F *[5];
		fCEtaPhi_LS = new TH2F *[5];
		fCEtaPhi_ULS_NoP = new TH2F *[5];
		fCEtaPhi_LS_NoP = new TH2F *[5];
		
		fCEtaPhi_ULS_Weight = new TH2F *[5];
		fCEtaPhi_LS_Weight = new TH2F *[5];
		fCEtaPhi_ULS_NoP_Weight = new TH2F *[5];
		fCEtaPhi_LS_NoP_Weight = new TH2F *[5];
		
		fCEtaPhi_Inc_EM = new TH2F *[5];
		
		fCEtaPhi_ULS_EM = new TH2F *[5];
		fCEtaPhi_LS_EM = new TH2F *[5];
		
		fCEtaPhi_ULS_Weight_EM = new TH2F *[5];
		fCEtaPhi_LS_Weight_EM = new TH2F *[5];
		
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
	
	for(Int_t i = 0; i < 5; i++)
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
	
		//__________________________________________________________________
		//Efficiency studies
	if(fIsMC)
	{
		fPtBackgroundBeforeReco = new TH1F("fPtBackgroundBeforeReco",";p_{T} (GeV/c);Count",300,0,30);
		if(fFillBackground)fPtBackgroundBeforeReco2 = new TH1F("fPtBackgroundBeforeReco2",";p_{T} (GeV/c);Count",300,0,30);
		fPtBackgroundAfterReco = new TH1F("fPtBackgroundAfterReco",";p_{T} (GeV/c);Count",300,0,30);	
		fPtMCparticleAll = new TH1F("fPtMCparticleAll",";p_{T} (GeV/c);Count",200,0,40);	
		fPtMCparticleReco = new TH1F("fPtMCparticleReco",";p_{T} (GeV/c);Count",200,0,40);
		
		fPtMCparticleAll_nonPrimary = new TH1F("fPtMCparticleAll_nonPrimary",";p_{T} (GeV/c);Count",200,0,40);	
		fPtMCparticleAlle_nonPrimary = new TH1F("fPtMCparticleAlle_nonPrimary",";p_{T} (GeV/c);Count",200,0,40);	
		fPtMCparticleReco_nonPrimary = new TH1F("fPtMCparticleReco_nonPrimary",";p_{T} (GeV/c);Count",200,0,40);
		
		fPtMCparticleAllHfe1 = new TH1F("fPtMCparticleAllHfe1",";p_{t} (GeV/c);Count",200,0,40);
		fPtMCparticleRecoHfe1 = new TH1F("fPtMCparticleRecoHfe1",";p_{t} (GeV/c);Count",200,0,40);
		fPtMCparticleAllHfe2 = new TH1F("fPtMCparticleAllHfe2",";p_{t} (GeV/c);Count",200,0,40);
		fPtMCparticleRecoHfe2 = new TH1F("fPtMCparticleRecoHfe2",";p_{t} (GeV/c);Count",200,0,40);
		
		fPtMCelectronAfterAll = new TH1F("fPtMCelectronAfterAll",";p_{T} (GeV/c);Count",200,0,40);
		fPtMCelectronAfterAll_nonPrimary = new TH1F("fPtMCelectronAfterAll_nonPrimary",";p_{T} (GeV/c);Count",200,0,40);
	

		
		fPtMCpi0 = new TH1F("fPtMCpi0",";p_{t} (GeV/c);Count",200,0,30);
		fPtMCeta = new TH1F("fPtMCeta",";p_{T} (GeV/c);Count",200,0,30);
		fPtMC_EMCal_All= new TH1F("fPtMC_EMCal_All",";p_{t} (GeV/c);Count",200,0,40);
		fPtMC_EMCal_Selected= new TH1F("fPtMC_EMCal_Selected",";p_{t} (GeV/c);Count",200,0,40);
		fPtMC_TPC_All= new TH1F("fPtMC_TPC_All",";p_{t} (GeV/c);Count",200,0,40);
		fPtMC_TPC_Selected = new TH1F("fPtMC_TPC_Selected",";p_{t} (GeV/c);Count",200,0,40);
		fPtMCWithLabel = new TH1F("fPtMCWithLabel",";p_{t} (GeV/c);Count",200,0,40);
		fPtMCWithoutLabel = new TH1F("fPtMCWithoutLabel",";p_{t} (GeV/c);Count",200,0,40);
		fPtIsPhysicaPrimary = new TH1F("fPtIsPhysicaPrimary",";p_{t} (GeV/c);Count",200,0,40);
		
		fOutputList->Add(fPtBackgroundBeforeReco);
		if(fFillBackground) fOutputList->Add(fPtBackgroundBeforeReco2);
		fOutputList->Add(fPtBackgroundAfterReco);
		fOutputList->Add(fPtMCparticleAll);
		fOutputList->Add(fPtMCparticleReco);
		
		fOutputList->Add(fPtMCparticleAll_nonPrimary);
		fOutputList->Add(fPtMCparticleAlle_nonPrimary);
		fOutputList->Add(fPtMCparticleReco_nonPrimary);
		
		fOutputList->Add(fPtMCparticleAllHfe1);
		fOutputList->Add(fPtMCparticleRecoHfe1);
		fOutputList->Add(fPtMCparticleAllHfe2);
		fOutputList->Add(fPtMCparticleRecoHfe2);
		fOutputList->Add(fPtMCelectronAfterAll);
		
		fOutputList->Add(fPtMCelectronAfterAll_nonPrimary);
		
		
		fOutputList->Add(fPtMCpi0);
		fOutputList->Add(fPtMCeta);
		fOutputList->Add(fPtMC_EMCal_All);
		fOutputList->Add(fPtMC_EMCal_Selected);
		fOutputList->Add(fPtMC_TPC_All);
		fOutputList->Add(fPtMC_TPC_Selected);
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
		const AliAODVertex* spdVtx = fAOD->GetPrimaryVertexSPD();
		if(spdVtx->GetNContributors()<=0) return;
		TString vtxTyp = spdVtx->GetTitle();
		Double_t cov[6]={0};
		spdVtx->GetCovarianceMatrix(cov);
		Double_t zRes = TMath::Sqrt(cov[5]);
		if(vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return;
		if(TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return;
		if(TMath::Abs(zvtx) > 10) return;
		
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
	if(fNOtrks<2) return;
	
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
				
				Int_t pdg = fMCparticle->GetPdgCode();
				if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax && fMCparticle->Charge()!=0)
				{
					if (TMath::Abs(pdg) == 11) fPtMCparticleAlle_nonPrimary->Fill(fMCparticle->Pt()); //denominator for total efficiency for all electrons, and not primary
					
					if( TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 2212 || TMath::Abs(pdg) == 321 || TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13 ) 
					{
						
						fPtMCparticleAll_nonPrimary->Fill(fMCparticle->Pt()); //denominator for total efficiency for all particles, and not primary
						if(fMCparticle->IsPhysicalPrimary()) 
						{
							fPtMCparticleAll->Fill(fMCparticle->Pt());
							
							Bool_t MotherFound = FindMother(iMC);
							if(MotherFound)
							{
								if(fIsHFE1) fPtMCparticleAllHfe1->Fill(fMCparticle->Pt()); //denominator for total efficiency
								if(fIsHFE2) fPtMCparticleAllHfe2->Fill(fMCparticle->Pt());
							}
						}
					}
				}
				if(TMath::Abs(pdg)==111) fPtMCpi0->Fill(fMCparticle->Pt());
				if(TMath::Abs(pdg)==211) fPtMCeta->Fill(fMCparticle->Pt());
			}
		}
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
				Int_t pdg = fMCtrack->GetPdgCode();
				
				if(TMath::Abs(pdg)==111) fPtMCpi0->Fill(fMCtrack->Pt());
				if(TMath::Abs(pdg)==211) fPtMCeta->Fill(fMCtrack->Pt());

				
								
				
				if(fMCtrack->Eta()>=fEtaCutMin && fMCtrack->Eta()<=fEtaCutMax)
				{
					
					if (TMath::Abs(pdg) == 11)  fPtMCparticleAlle_nonPrimary->Fill(fMCtrack->Pt());//denominator for total efficiency for all electrons, and not primary

					
					if( TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 2212 || TMath::Abs(pdg) == 321 || TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13 )
					{
						fPtMCparticleAll_nonPrimary->Fill(fMCtrack->Pt());//denominator for total efficiency for all particle, non Primary track
						
						if(!fMCstack->IsPhysicalPrimary(iMC)) continue;
						fPtMCparticleAll->Fill(fMCtrack->Pt());
						
						Bool_t MotherFound = FindMother(iMC);
						if(MotherFound)
						{
							if(fIsHFE1) fPtMCparticleAllHfe1->Fill(fMCtrack->Pt());//denominator for total efficiency
							if(fIsHFE2) fPtMCparticleAllHfe2->Fill(fMCtrack->Pt());
						}
					}	
				}
	        }
		}
	}
	
		//______________________________________________________________________
		//EMCal Trigger Selection (Threshould selection)
	TString firedTrigger;
	TString TriggerEG1("EG1");
	TString TriggerEG2("EG2");
	
	if(fAOD) firedTrigger = fAOD->GetFiredTriggerClasses();
	else if(fESD) firedTrigger = fESD->GetFiredTriggerClasses();
	
	fNevent->Fill(0);
	if(firedTrigger.Contains(TriggerEG1)) fNevent->Fill(1);
	if(firedTrigger.Contains(TriggerEG2)) fNevent->Fill(2);
	
		//EG1
	if(firedTrigger.Contains(TriggerEG1))
	{ 
		fNevent->Fill(3);
	}
	else 
	{
		if(fEMCEG1) return;
	}
	
		//EG2
	if(firedTrigger.Contains(TriggerEG2))
	{ 
		fNevent->Fill(4);
	}
	else
	{ 
		if(fEMCEG2) return;
	}
	
		//______________________________________________________________________
	
	Int_t ClsNo = -999;
	if(!fIsAOD) ClsNo = fESD->GetNumberOfCaloClusters(); 
	else ClsNo = fAOD->GetNumberOfCaloClusters(); 
	
		//______________________________________________________________________
	
		///______________________________________________________________________
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
		
		Double_t fTPCnSigma = -999;
		Double_t fTPCnSigma_pion = -999;
		Double_t fTPCnSigma_proton = -999;
		Double_t fTPCnSigma_kaon = -999;
		Double_t fTPCsignal = -999;
		Double_t fPt = -999;
		Double_t fP = -999;
		
			///_____________________________________________________________________________
			///Fill QA plots without track selection
		fPt = track->Pt();
		fP = TMath::Sqrt((track->Pt())*(track->Pt()) + (track->Pz())*(track->Pz()));
		
		fTPCsignal = track->GetTPCsignal();
		fTPCnSigma = fPidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
		fTPCnSigma_pion = fPidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
		fTPCnSigma_proton = fPidResponse->NumberOfSigmasTPC(track, AliPID::kProton);
		fTPCnSigma_kaon = fPidResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
		
		fTPC_p[0]->Fill(fPt,fTPCsignal);
		fTPCnsigma_p[0]->Fill(fP,fTPCnSigma);
		
		
		Float_t TPCNcls = track->GetTPCNcls();
		Float_t pos[3]={0,0,0};
		
		Double_t fEMCflag = kFALSE;
		if(track->GetEMCALcluster()>0)
		{
			fClus = fVevent->GetCaloCluster(track->GetEMCALcluster());
			if(fClus->IsEMCAL())
			{
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
					
					fECluster[0]->Fill(Energy);
					fTPCNcls_EoverP[0]->Fill(TPCNcls, EoverP);
				}
			}
		}
		
			//______________________________________________________________
			// Vertex
		
		fVtxZ[0]->Fill(fZvtx);
		fNTracks[0]->Fill(fNOtrks);
		fNClusters[0]->Fill(ClsNo);
		
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
		
			///_____________________________________________________________
			///QA plots after track selection
		if(fIsMC)
		{
			if(track->GetLabel()>=0) fPtMCWithLabel->Fill(fPt);
			else fPtMCWithoutLabel->Fill(fPt);
		}
		
		if(fIsMC && fIsAOD && track->GetLabel()>=0)
		{
			fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
			
			if(fMCparticle->IsPhysicalPrimary()) fPtIsPhysicaPrimary->Fill(fPt);
			
			Int_t pdg = fMCparticle->GetPdgCode();
			if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax && fMCparticle->Charge()!=0)
			{
				
				fPtMCparticleReco_nonPrimary->Fill(fMCparticle->Pt()); //not Primary track
				
				if( TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 2212 || TMath::Abs(pdg) == 321 || TMath::Abs(pdg) == 11 || TMath::Abs(pdg) == 13 ) 
				{	
					fPtMCparticleReco_nonPrimary->Fill(fMCparticle->Pt()); //not Primary track
					
					if(fMCparticle->IsPhysicalPrimary()) 
					{
						fPtMCparticleReco->Fill(fMCparticle->Pt());
						
						Bool_t MotherFound = FindMother(track->GetLabel());
						if(MotherFound)
						{
							if(fIsHFE1) fPtMCparticleRecoHfe1->Fill(fMCparticle->Pt());
							if(fIsHFE2) fPtMCparticleRecoHfe2->Fill(fMCparticle->Pt());
						}
					}
				}
			}
		}
		else if(fIsMC && track->GetLabel()>=0)
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
								if(fIsHFE1) fPtMCparticleRecoHfe1->Fill(fMCtrack->Pt());
								if(fIsHFE2) fPtMCparticleRecoHfe2->Fill(fMCtrack->Pt());
							}
						}
				}
			}
		}
		
		fTPC_p[1]->Fill(fPt,fTPCsignal);
		fTPCnsigma_p[1]->Fill(fP,fTPCnSigma);
		Double_t fPtBin[6] = {2,4,6,8,10,15};
		
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
						
						if(fTPCnSigma < -3.5){
							fEoverP_pt_pions->Fill(fPt, EoverP);
							
						}
						
						if(fTPCnSigma < -3.5 && fTPCnSigma > -10){
							fEoverP_pt_pions2->Fill(fPt, EoverP);
							
						}
						
						
					}
					
					
					
					
					for(Int_t i = 0; i < 5; i++)
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
		fNTracks[1]->Fill(fNOtrks);
		fNClusters[1]->Fill(ClsNo);
			//______________________________________________________________
		
			///______________________________________________________________________
			///Histograms for PID Studies
			//Double_t fPtBin[6] = {2,4,6,8,10,15};
		
		for(Int_t i = 0; i < 5; i++)
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
				        if( TMath::Abs(pdg) == 11 && fMCtrackMother->GetPdgCode()!=22 ){
							if(fIsHFE1) fPtMC_TPC_All->Fill(fMCtrack->Pt());	
				        }
					}
				}
			}
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
			///________________________________________________________________________		
		
		
			////////////////////////////////////////////////////////////////////
			///TPC efficiency calculations 
		
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
				        if( TMath::Abs(pdg) == 11 && fMCtrackMother->GetPdgCode()!=22 ){
							if(fIsHFE1) fPtMC_TPC_Selected->Fill(fMCtrack->Pt());	
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
			Background(track, iTracks, Vtrack, kTRUE); //IsTPConly=kTRUE
		}
		
		
		fTPCnsigma_p[2]->Fill(fP,fTPCnSigma);
		fTPC_p[2]->Fill(fP,fTPCsignal);
		TPCNcls = track->GetTPCNcls();
		Float_t pos3[3]={0,0,0};
		
		if(track->GetEMCALcluster()>0)
		{
			fClus = fVevent->GetCaloCluster(track->GetEMCALcluster());
			if(fClus->IsEMCAL())
			{
				
					//________________________________________________________________________		
				
				
					//Cluster timing distribution
				if(fUseEMCal && !fIsAOD ){
					AliESDCaloCells &cells=*(fESD->GetEMCALCells());
						//	Int_t nTotalCells = cells.GetNumberOfCells() ;  
						//Int_t type        = cells.GetType();
						//for (Int_t icell=  0; icell <  nTotalCells; icell++) {
						//fTime->Fill(cells.GetTime(icell));
						//}
					
					TRefArray* caloClusters = new TRefArray();
					fESD->GetEMCALClusters(caloClusters);
					
					Int_t nclus = caloClusters->GetEntries();
					
						//fClusESD = fESD->GetCaloCluster(track->GetEMCALcluster());
					
					for (Int_t icl = 0; icl < nclus; icl++) {
						
						AliESDCaloCluster* clus = (AliESDCaloCluster*)caloClusters->At(icl);
						
						if(fClus->IsEMCAL()){
							Float_t ncells	= fClus->GetNCells();
							UShort_t * index = clus->GetCellsAbsId() ;
							UShort_t * index2 = fClus->GetCellsAbsId() ;
							
							
							for(Int_t i = 0; i < ncells ; i++){
								
								Int_t absId =   index[i];
								fTime->Fill(fPt,cells.GetCellTime(absId));
								
								Int_t absId2 =   index2[i];
								fTime2->Fill(fPt,cells.GetCellTime(absId2));
							}
							
						}
					}
					
					
					
					
				}
				
				if(fUseEMCal){
					double emctof = fClus->GetTOF();
					ftimingEle->Fill(fPt,emctof);
				}
					//________________________________________________________________________		
				
				
				
				
					/////////////// Residuals
				Double_t Dx = fClus->GetTrackDx();
				Double_t Dz = fClus->GetTrackDz();
				Double_t R=TMath::Sqrt(Dx*Dx+Dz*Dz);
				for(Int_t i = 0; i < 5; i++)
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
					
					
					
						/////////////// for Eta Phi distribution
					fClus->GetPosition(pos3);
					TVector3 vpos(pos3[0],pos3[1],pos3[2]);
					Double_t cphi = vpos.Phi();
					Double_t ceta = vpos.Eta();
					fEtaPhi[2]->Fill(cphi,ceta);
					
					
					
					fTPCNcls_EoverP[2]->Fill(TPCNcls, EoverP);
					
					for(Int_t i = 0; i < 5; i++)
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
										fPtMC_EMCal_All->Fill(fMCparticle->Pt());	
									}
								}
							}
						}
					}///until here
					
					else if(fIsMC && track->GetLabel()>=0)//ESD
					{
						
						if(fMCstack->IsPhysicalPrimary(track->GetLabel()))
						{
							
							fMCtrack = fMCstack->Particle(track->GetLabel());
							
							Int_t pdg = fMCtrack->GetPdgCode();
							if(fMCtrack->Eta()>=fEtaCutMin && fMCtrack->Eta()<=fEtaCutMax && fMCtrack->Phi()>=(TMath::Pi()*80/180) && fMCtrack->Phi()<=TMath::Pi())
							{
								if(fMCtrack->GetFirstMother()>0){
									fMCtrackMother = fMCstack->Particle(fMCtrack->GetFirstMother());
									if( TMath::Abs(pdg) == 11 && fMCtrackMother->GetPdgCode()!=22 ){
										
										fPtMC_EMCal_All->Fill(fMCtrack->Pt());	
									}
								}
							}
						}
					}
					
					
					if((fClus->E() / fP) >= fEoverPCutMin && (fClus->E() / fP) <= fEoverPCutMax)
					{	
						
					    fECluster[2]->Fill(Energy);
							//_______________________________________________________
							//Correlation Analysis
						if(fUseEMCal)
						{
							fPtElec_Inc->Fill(fPt);
								//new function to fill non-HFE histos
							if(fFillBackground){
								Background(track, iTracks, Vtrack, kFALSE);
							}
							
							double emctof2 = fClus->GetTOF();
							ftimingEle2->Fill(fPt,emctof2);
							
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
						
						
							/// changing start here
						if(fIsMC && fIsAOD && track->GetLabel()>=0)
						{
							if(track->Charge()<0)  fCharge_n->Fill(fPt);
							if(track->Charge()>0)  fCharge_p->Fill(fPt);
							
							fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
							Int_t pdg = fMCparticle->GetPdgCode();
							
							if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax ){
									
								if( TMath::Abs(pdg) == 11) fPtMCelectronAfterAll_nonPrimary->Fill(fMCparticle->Pt()); //numerator for the total efficiency, non Primary track
							}	
								
							
							//
							if(fMCparticle->IsPhysicalPrimary()){
								
								
								if(fMCparticle->Eta()>=fEtaCutMin && fMCparticle->Eta()<=fEtaCutMax ){
										//check
									if(fIsHFE1) fPtMCelectronAfterAll->Fill(fMCparticle->Pt()); //numerator for the total efficiency
									
									
									Bool_t MotherFound = FindMother(track->GetLabel());
									if(MotherFound){
										fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
										if( TMath::Abs(pdg) == 11 && fMCparticleMother->GetPdgCode()!=22 ){
											fPtMC_EMCal_Selected->Fill(fMCparticle->Pt());	
										}
									}
								}
							}
						}///until here
						
						else if(fIsMC && track->GetLabel()>=0)//ESD
						{
							if(track->Charge()<0)  fCharge_n->Fill(fPt);
							if(track->Charge()>0)  fCharge_p->Fill(fPt);
							
							fMCtrack = fMCstack->Particle(track->GetLabel());
							Int_t pdg = fMCtrack->GetPdgCode();
							
							if(fMCtrack->Eta()>=fEtaCutMin && fMCtrack->Eta()<=fEtaCutMax){
								if( TMath::Abs(pdg) == 11) fPtMCelectronAfterAll_nonPrimary->Fill(fMCtrack->Pt()); //numerator for the total efficiency, non Primary track
							}
							
							if(fMCstack->IsPhysicalPrimary(track->GetLabel()))
							{
								Bool_t MotherFound = FindMother(track->GetLabel());
							    
								if(MotherFound)
								{
									if(fMCtrack->Eta()>=fEtaCutMin && fMCtrack->Eta()<=fEtaCutMax){
										if(fIsHFE1) fPtMCelectronAfterAll->Fill(fMCtrack->Pt()); //numerator for the total efficiency
									}
								}
								
								
								
								
								if(fMCtrack->Eta()>=fEtaCutMin && fMCtrack->Eta()<=fEtaCutMax && fMCtrack->Phi()>=(TMath::Pi()*80/180) && fMCtrack->Phi()<=TMath::Pi())
								{
									if(fMCtrack->GetFirstMother()>0){
										fMCtrackMother = fMCstack->Particle(fMCtrack->GetFirstMother());
										if( TMath::Abs(pdg) == 11 && fMCtrackMother->GetPdgCode()!=22 ){
											
											fPtMC_EMCal_Selected->Fill(fMCtrack->Pt());	
										}
									}
								}
							}
						}
							///////////////////////////////////////////////////////////////////
						
						
					}
				}
			}
		}
		
			//______________________________________________________________
			// Vertex
		
		fVtxZ[2]->Fill(fZvtx);
		fNTracks[2]->Fill(fNOtrks);
		fNClusters[2]->Fill(ClsNo);
		
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
				if(!IsTPConly)fPtBackgroundBeforeReco->Fill(track->Pt());
				if(IsTPConly)fPtBackgroundBeforeReco2->Fill(track->Pt());
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
				if(!IsTPConly)fPtBackgroundBeforeReco->Fill(track->Pt());
				if(IsTPConly)fPtBackgroundBeforeReco2->Fill(track->Pt());
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
	}
		///_________________________________________________________________
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
	
	
	
}


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
					
					if(fPtH>fPtE || fPtH<1) continue;
					
					fDphi = fPhiE - fPhiH;
					
					if (fDphi > 3*pi/2) fDphi = fDphi - 2*pi;
					if (fDphi < -pi/2)  fDphi = fDphi + 2*pi;
					
					fDeta = fEtaE - fEtaH;
					
					Double_t fPtBin[6] = {2,4,6,8,10,15};
					
					for(Int_t i = 0; i < 5; i++)
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
		
		if(fIsAOD) 
		{
			AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(Vtrack2);
			if(!atrack2->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
			if((!(atrack2->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack2->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
			if(atrack2->GetTPCNcls() < 80) continue; 
		}
		else
		{   
			AliESDtrack *etrack2 = dynamic_cast<AliESDtrack*>(Vtrack2); 
			if(!fPartnerCuts->AcceptTrack(etrack2)) continue; 
		}
		
		fPhiH = track2->Phi();
		fEtaH = track2->Eta();
		fPtH = track2->Pt();
		
		if(fPtH>fPtE || fPtH<1) continue;
		
		fDphi = fPhiE - fPhiH;
		
		if (fDphi > 3*pi/2) fDphi = fDphi - 2*pi;
		if (fDphi < -pi/2)  fDphi = fDphi + 2*pi;
		
		fDeta = fEtaE - fEtaH;
		
		Double_t fPtBin[6] = {2,4,6,8,10,15};
		
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
		
		for(Int_t i = 0; i < 5; i++)
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
		
		if(fIsAOD) 
		{
			AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(Vtrack2);
			if(!atrack2->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
			if((!(atrack2->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack2->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
			if(atrack2->GetTPCNcls() < 80) continue; 
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
		
		if(fIsAOD) 
		{
			AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(Vtrack2);
			if(!atrack2->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
			if((!(atrack2->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack2->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
			if(atrack2->GetTPCNcls() < 80) continue; 
		}
		else
		{   
			AliESDtrack *etrack2 = dynamic_cast<AliESDtrack*>(Vtrack2); 
			if(!fPartnerCuts->AcceptTrack(etrack2)) continue; 
		}
		
		fPhiH = track2->Phi();
		fEtaH = track2->Eta();
		fPtH = track2->Pt();
		
		if(fPtH>fPtE || fPtH<1) continue;
		
		fDphi = fPhiE - fPhiH;
		
		if (fDphi > 3*pi/2) fDphi = fDphi - 2*pi;
		if (fDphi < -pi/2)  fDphi = fDphi + 2*pi;
		
		fDeta = fEtaE - fEtaH;
		
		Double_t fPtBin[6] = {2,4,6,8,10,15};
		
		for(Int_t i = 0; i < 5; i++)
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
