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
//      Task for Heavy-flavour electron analysis in pp collisions     //
//      (and also Pb-Pb)             								  //
//																	  //
//		v1.0														  //
//                                                                    //
//	    Authors 							                          //
//		Camila de Conti (camila.de.conti@cern.ch)				      //
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
#include "AliAnalysisHFETPCTOF.h"
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
#include "AliMultSelection.h"
//______________________________________________________________________

//______________________________________________________________________
ClassImp(AliAnalysisHFETPCTOF)

//______________________________________________________________________
AliAnalysisHFETPCTOF::AliAnalysisHFETPCTOF(const char *name) 
  : AliAnalysisTaskSE(name)
,fIsMC(0)


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
,fPidResponse(0)
,fNonHFE(new AliSelectNonHFE())
,fIsAOD(kFALSE)
,fIsPP(kFALSE)
,fZvtx(0)
,fClus(0)

//Elienos class for nonHFE reconstruction by IM

,fMassCutFlag(kTRUE)
,fAngleCutFlag(kFALSE)
,fMassCut(0.1)
,fAngleCut(999)
,fChi2CutFlag(kFALSE)
,fChi2Cut(3.5)
,fDCAcutFlag(kFALSE)
,fDCAcut(999)//dca between two tracks
,fPartnerCuts(new AliESDtrackCuts())
,fPtMinAsso(0.3)
,fTpcNclsAsso(80)
,fInvMass(0)
,fInvMassBack(0)
,fInvMass1(0)
,fInvMassBack1(0)
,fInvMass2(0)
,fInvMassBack2(0)
,fInvMass3(0)
,fInvMassBack3(0)
,fInvMass4(0)
,fInvMassBack4(0)
,fDCA(0)
,fDCABack(0)
,fAngle(0)
,fAngleBack(0)
,fPtElec_ULS(0)
,fPtElec_LS(0)
,fTPCnSigmaPElec_ULS(0)
,fTPCnSigmaPElec_LS(0)
,fPtElec_ULS_MC(0)
,fPtElec_LS_MC(0)
,fPtElec_ULS_MC_from_pi0(0)
,fPtElec_LS_MC_from_pi0(0)
,fPtElec_ULS_MC_from_eta(0)
,fPtElec_LS_MC_from_eta(0)
,fPtElec_ULS_MC_from_gamma(0)
,fPtElec_LS_MC_from_gamma(0)
,fPtElec(0)


//Histograms for the analysis
,fNevent(0)
,fNevent2(0)
,fCent(0)
,fTPC_p1(0)
,fTPC_p2(0)
,fTPC_p3(0)
,fTPCnsigma_p1(0)
,fTPCnsigma_p2(0)
,fTPCnsigma_p3(0)
,fTPCnsigma_pt1(0)
,fTPCnsigma_pt2(0)
,fTPCnsigma_pt3(0)
,fTPCnsigma_p_after_tof(0)
,fTPCnsigma_pt_after_tof(0)
,fTPCnsigma_p_after_tof_its(0)
,fTPCnsigma_pt_after_tof_its(0)
,fTPCnsigma_p_after_tof_inverse(0)
,fTPCnsigma_pt_after_tof_inverse(0)
,fTOFnsigma_p1(0)
,fTOFnsigma_p2(0)
,fTOFnsigma_p3(0)
,fTOFnsigma_pt1(0)
,fTOFnsigma_pt2(0)
,fTOFnsigma_pt3(0)
,fITSnsigma_p1(0)
,fITSnsigma_p2(0)
,fITSnsigma_p3(0)
,fITSnsigma_pt1(0)
,fITSnsigma_pt2(0)
,fITSnsigma_pt3(0)
//,fTPCnsigma_EoverP(0)
,fTPCnsigma_TOFnsigma1(0)
,fTPCnsigma_TOFnsigma2(0)
,fTPCnsigma_TOFnsigma3(0)
//,fVtxZ(0)
//,fNClusters(0)
//,fEtaPhiTracks(0)
//,fEtaPhiClusters(0)
,fDCAxy_pt(0)
,fDCAz_pt(0)
,fDCAxy_pt_gamma(0)
,fDCAz_pt_gamma(0)
,fDCAxy_pt_light(0)
,fDCAz_pt_light(0)
,fDCAxy_pt_B(0)
,fDCAz_pt_B(0)
,fDCAxy_pt_C(0)
,fDCAz_pt_C(0)
,fPt_elec_phot(0)
,fPt_elec_from_pi0(0)
,fPt_elec_from_eta(0)
,fPt_elec_from_gamma(0)
,fPt_HFelec_passTOFandTPC_2(0)
,fPt_HFelec_passTPC_2(0)
,fPt_HFelec_pass_track(0)
,fPt_HFelec_pass_trackDCA(0)
,fPt_HFelec_pass_trackITSlayer(0)
,fPt_HFelec_pass_trackFB(0)
,fPt_HFelec_pass_trackTPCITS(0)
,fPt_HFelec_pass_trackKink(0)
,fPt_HFelec_pass_trackandTOF(0)


//For the HFE package
,fCuts(0)
,fCFM(0)
,fPID(new AliHFEpid("hfePid"))
,fPIDqa(0)

//For MC
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
,fPtMCpi0(0)
,fPtMCpi0_feeddown(0)
,fPtMCpi0_nomother(0)
,fPtMCeta(0)
,fPtMCeta_feeddown(0)
,fPtMCeta_nomother(0)
,fPtMCpi0_afterNorm(0)
,fPtMCeta_afterNorm(0)
,fPtMCpi02(0)
,fPtMCeta2(0)
,fPtMCpi03(0)
,fPtMCeta3(0)
,fPtHFEMC_2(0)

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
AliAnalysisHFETPCTOF::AliAnalysisHFETPCTOF() 
  : AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisHFETPCTOF")

,fIsMC(0)

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
,fPidResponse(0)
,fNonHFE(new AliSelectNonHFE())
,fIsAOD(kFALSE)
,fIsPP(kFALSE)
,fZvtx(0)
,fClus(0)

//Elienos class for nonHFE reconstruction by IM
,fMassCutFlag(kTRUE)
,fAngleCutFlag(kFALSE)
,fMassCut(0.1)
,fAngleCut(999)
,fChi2CutFlag(kFALSE)
,fChi2Cut(3.5)
,fDCAcutFlag(kFALSE)
,fDCAcut(999)//dca between two tracks
,fPartnerCuts(new AliESDtrackCuts())
,fPtMinAsso(0.3)
,fTpcNclsAsso(80)
,fInvMass(0)
,fInvMassBack(0)
,fInvMass1(0)
,fInvMassBack1(0)
,fInvMass2(0)
,fInvMassBack2(0)
,fInvMass3(0)
,fInvMassBack3(0)
,fInvMass4(0)
,fInvMassBack4(0)
,fDCA(0)
,fDCABack(0)
,fAngle(0)
,fAngleBack(0)
,fPtElec_ULS(0)
,fPtElec_LS(0)
,fTPCnSigmaPElec_ULS(0)
,fTPCnSigmaPElec_LS(0)
,fPtElec_ULS_MC(0)
,fPtElec_LS_MC(0)
,fPtElec_ULS_MC_from_pi0(0)
,fPtElec_LS_MC_from_pi0(0)
,fPtElec_ULS_MC_from_eta(0)
,fPtElec_LS_MC_from_eta(0)
,fPtElec_ULS_MC_from_gamma(0)
,fPtElec_LS_MC_from_gamma(0)
,fPtElec(0)


//Histograms for the analysis
,fNevent(0)
,fNevent2(0)
,fCent(0)
,fTPC_p1(0)
,fTPC_p2(0)
,fTPC_p3(0)
,fTPCnsigma_p1(0)
,fTPCnsigma_p2(0)
,fTPCnsigma_p3(0)
,fTPCnsigma_pt1(0)
,fTPCnsigma_pt2(0)
,fTPCnsigma_pt3(0)
,fTPCnsigma_p_after_tof(0)
,fTPCnsigma_pt_after_tof(0)
,fTPCnsigma_p_after_tof_its(0)
,fTPCnsigma_pt_after_tof_its(0)
,fTPCnsigma_p_after_tof_inverse(0)
,fTPCnsigma_pt_after_tof_inverse(0)
,fTOFnsigma_p1(0)
,fTOFnsigma_p2(0)
,fTOFnsigma_p3(0)
,fTOFnsigma_pt1(0)
,fTOFnsigma_pt2(0)
,fTOFnsigma_pt3(0)
,fITSnsigma_p1(0)
,fITSnsigma_p2(0)
,fITSnsigma_p3(0)
,fITSnsigma_pt1(0)
,fITSnsigma_pt2(0)
,fITSnsigma_pt3(0)
//,fTPCnsigma_EoverP(0)
,fTPCnsigma_TOFnsigma1(0)
,fTPCnsigma_TOFnsigma2(0)
,fTPCnsigma_TOFnsigma3(0)
//,fVtxZ(0)
//,fNClusters(0)
//,fEtaPhiTracks(0)
//,fEtaPhiClusters(0)
,fDCAxy_pt(0)
,fDCAz_pt(0)
,fDCAxy_pt_gamma(0)
,fDCAz_pt_gamma(0)
,fDCAxy_pt_light(0)
,fDCAz_pt_light(0)
,fDCAxy_pt_B(0)
,fDCAz_pt_B(0)
,fDCAxy_pt_C(0)
,fDCAz_pt_C(0)
,fPt_elec_phot(0)
,fPt_elec_from_pi0(0)
,fPt_elec_from_eta(0)
,fPt_elec_from_gamma(0)
,fPt_HFelec_passTOFandTPC_2(0)
,fPt_HFelec_passTPC_2(0)
,fPt_HFelec_pass_track(0)
,fPt_HFelec_pass_trackDCA(0)
,fPt_HFelec_pass_trackITSlayer(0)
,fPt_HFelec_pass_trackFB(0)
,fPt_HFelec_pass_trackTPCITS(0)
,fPt_HFelec_pass_trackKink(0)
,fPt_HFelec_pass_trackandTOF(0)

//For the HFE package
,fCuts(0)
,fCFM(0)
,fPID(new AliHFEpid("hfePid"))
,fPIDqa(0)

//For MC
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
,fPtMCpi0(0)
,fPtMCpi0_feeddown(0)
,fPtMCpi0_nomother(0)
,fPtMCeta(0)
,fPtMCeta_feeddown(0)
,fPtMCeta_nomother(0)
,fPtMCpi0_afterNorm(0)
,fPtMCeta_afterNorm(0)
,fPtMCpi02(0)
,fPtMCeta2(0)
,fPtMCpi03(0)
,fPtMCeta3(0)
,fPtHFEMC_2(0)


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
AliAnalysisHFETPCTOF::~AliAnalysisHFETPCTOF()
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
void AliAnalysisHFETPCTOF::UserCreateOutputObjects()
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
	fNevent2 = new TH1F("fNevent2","Number of Events 2",5,-0.5,4.5);
	fCent = new TH1F("fCent","Centrality",408,-1,101);
	//And then, add to the output list
	fOutputList->Add(fNevent);
	fOutputList->Add(fNevent2);
	fOutputList->Add(fCent);
	
	//General Histograms
	
	//Steps
	//Step 1: Before Track cuts
	//Step 2: Before PID
	//Step 3: After PID
	
	/*
	fTPC_p = new TH2F *[3];
	fTPCnsigma_p = new TH2F *[3];
	fTPCnsigma_pt = new TH2F *[3];
	fTOFnsigma_p = new TH2F *[3];
	fTOFnsigma_pt = new TH2F *[3];
	fITSnsigma_p = new TH2F *[3];
	fITSnsigma_pt = new TH2F *[3];
	//fTPCnsigma_EoverP = new TH2F *[3];
	fTPCnsigma_TOFnsigma = new TH2F *[3];
	//fVtxZ= new  TH1F *[3];
	//fNClusters= new TH1F *[3];
	//fEtaPhiTracks= new TH2F *[3];
	//fEtaPhiClusters= new TH2F *[3];
	
	
	///Initialize pointers to histograms
	
	for(Int_t i = 0; i < 3; i++)
	{   
		fTPC_p[i] = new TH2F(Form("fTPC_p%d",i),";p (GeV/c);TPC dE/dx (a. u.)",300,0,15,400,-20,200);
		fTPCnsigma_p[i] = new TH2F(Form("fTPCnsigma_p%d",i),";p (GeV/c);TPC Electron N#sigma",300,0,15,200,-15,10);
		fTPCnsigma_pt[i] = new TH2F(Form("fTPCnsigma_pt%d",i),";pt (GeV/c);TPC Electron N#sigma",300,0,15,200,-15,10);
		fTOFnsigma_p[i] = new TH2F(Form("fTOFnsigma_p%d",i),";p (GeV/c);TOF Electron N#sigma",300,0,15,200,-10,30);
		fTOFnsigma_pt[i] = new TH2F(Form("fTOFnsigma_pt%d",i),";pt (GeV/c);TOF Electron N#sigma",300,0,15,200,-10,30);
		fITSnsigma_p[i] = new TH2F(Form("fITSnsigma_p%d",i),";p (GeV/c);ITS Electron N#sigma",300,0,15,200,-10,30);
		fITSnsigma_pt[i] = new TH2F(Form("fITSnsigma_pt%d",i),";pt (GeV/c);ITS Electron N#sigma",300,0,15,200,-10,30);
		//fTPCnsigma_EoverP[i] = new TH2F(Form("fTPCnsigma_EoverP%d",i),";E/p;TPC Electron N#sigma",200,0,2,200,-15,10);
		fTPCnsigma_TOFnsigma[i] = new TH2F(Form("fTPCnsigma_TOFnsigma%d",i),";TOF Electron N#sigma;TPC Electron N#sigma",200,-10,30,200,-15,10);
		//fVtxZ[i]= new  TH1F(Form("fVtxZ%d",i),"VtxZ",500, -50,50);
		//fNClusters[i]= new TH1F(Form("fNClusters%d",i),"fNClusters0",200, 0,100);
		//fEtaPhiTracks[i]= new TH2F(Form("fEtaPhiTracks%d",i),"fEtaPhiTracks0",100,-3,3,100,0,7);
		//fEtaPhiClusters[i]= new TH2F(Form("fEtaPhiClusters%d",i),"fEtaPhiClusters0",100,-3,3,100,0,7);
	  		
		fOutputList->Add(fTPC_p[i]);
		fOutputList->Add(fTPCnsigma_p[i]);
		fOutputList->Add(fTPCnsigma_pt[i]);
		fOutputList->Add(fTOFnsigma_p[i]);
		fOutputList->Add(fTOFnsigma_pt[i]);
		fOutputList->Add(fITSnsigma_p[i]);
		fOutputList->Add(fITSnsigma_pt[i]);
		//fOutputList->Add(fTPCnsigma_EoverP[i]);
		fOutputList->Add(fTPCnsigma_TOFnsigma[i]);
		//fOutputList->Add(fVtxZ[i]);
		//fOutputList->Add(fNClusters[i]);
		//fOutputList->Add(fEtaPhiTracks[i]);
		//fOutputList->Add(fEtaPhiClusters[i]);

	}
	
	*/
	
	///*****************************************************************
	///Substituting the pointers to normal histograms

	fTPC_p1 = new TH2F("fTPC_p1","p (GeV/c);TPC dE/dx (a. u.)",300,0,15,400,-20,200);
	fOutputList->Add(fTPC_p1);
	
	fTPC_p2 = new TH2F("fTPC_p2","p (GeV/c);TPC dE/dx (a. u.)",300,0,15,400,-20,200);
	fOutputList->Add(fTPC_p2);
	
	fTPC_p3 = new TH2F("fTPC_p3","p (GeV/c);TPC dE/dx (a. u.)",300,0,15,400,-20,200);
	fOutputList->Add(fTPC_p3);
	
	fTPCnsigma_p1 = new TH2F("fTPCnsigma_p1","p (GeV/c);TPC Electron N#sigma",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnsigma_p1);
	
	fTPCnsigma_p2 = new TH2F("fTPCnsigma_p2","p (GeV/c);TPC Electron N#sigma",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnsigma_p2);
	
	fTPCnsigma_p3 = new TH2F("fTPCnsigma_p3","p (GeV/c);TPC Electron N#sigma",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnsigma_p3);
	
	fTPCnsigma_pt1 = new TH2F("fTPCnsigma_pt1","pt (GeV/c);TPC Electron N#sigma",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnsigma_pt1);
	
	fTPCnsigma_pt2 = new TH2F("fTPCnsigma_pt2","pt (GeV/c);TPC Electron N#sigma",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnsigma_pt2);
	
	fTPCnsigma_pt3 = new TH2F("fTPCnsigma_pt3","pt (GeV/c);TPC Electron N#sigma",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnsigma_pt3);
	
	fTOFnsigma_p1 = new TH2F("fTOFnsigma_p1","p (GeV/c);TOF Electron N#sigma",300,0,15,200,-10,30);
	fOutputList->Add(fTOFnsigma_p1);
	
	fTOFnsigma_p2 = new TH2F("fTOFnsigma_p2","p (GeV/c);TOF Electron N#sigma",300,0,15,200,-10,30);
	fOutputList->Add(fTOFnsigma_p2);
	
	fTOFnsigma_p3 = new TH2F("fTOFnsigma_p3","p (GeV/c);TOF Electron N#sigma",300,0,15,200,-10,30);
	fOutputList->Add(fTOFnsigma_p3);
	
	fTOFnsigma_pt1 = new TH2F("fTOFnsigma_pt1","pt (GeV/c);TOF Electron N#sigma",300,0,15,200,-10,30);
	fOutputList->Add(fTOFnsigma_pt1);
	
	fTOFnsigma_pt2 = new TH2F("fTOFnsigma_pt2","pt (GeV/c);TOF Electron N#sigma",300,0,15,200,-10,30);
	fOutputList->Add(fTOFnsigma_pt2);
	
	fTOFnsigma_pt3 = new TH2F("fTOFnsigma_pt3","pt (GeV/c);TOF Electron N#sigma",300,0,15,200,-10,30);
	fOutputList->Add(fTOFnsigma_pt3);
	
	fITSnsigma_p1 = new TH2F("fITSnsigma_p1","p (GeV/c);ITS Electron N#sigma",300,0,15,200,-10,30);
	fOutputList->Add(fITSnsigma_p1);
	
	fITSnsigma_p2 = new TH2F("fITSnsigma_p2","p (GeV/c);ITS Electron N#sigma",300,0,15,200,-10,30);
	fOutputList->Add(fITSnsigma_p2);
	
	fITSnsigma_p3 = new TH2F("fITSnsigma_p3","p (GeV/c);ITS Electron N#sigma",300,0,15,200,-10,30);
	fOutputList->Add(fITSnsigma_p3);
	
	fITSnsigma_pt1 = new TH2F("fITSnsigma_pt1","pt (GeV/c);ITS Electron N#sigma",300,0,15,200,-10,30);
	fOutputList->Add(fITSnsigma_pt1);
	
	fITSnsigma_pt2 = new TH2F("fITSnsigma_pt2","pt (GeV/c);ITS Electron N#sigma",300,0,15,200,-10,30);
	fOutputList->Add(fITSnsigma_pt2);
	
	fITSnsigma_pt3 = new TH2F("fITSnsigma_pt3","pt (GeV/c);ITS Electron N#sigma",300,0,15,200,-10,30);
	fOutputList->Add(fITSnsigma_pt3);
	
	fTPCnsigma_TOFnsigma1 = new TH2F("fTPCnsigma_TOFnsigma1","TOF Electron N#sigma;TPC Electron N#sigma",200,-10,30,200,-15,10);
	fOutputList->Add(fTPCnsigma_TOFnsigma1);
	
	fTPCnsigma_TOFnsigma2 = new TH2F("fTPCnsigma_TOFnsigma2","TOF Electron N#sigma;TPC Electron N#sigma",200,-10,30,200,-15,10);
	fOutputList->Add(fTPCnsigma_TOFnsigma2);
	
	fTPCnsigma_TOFnsigma3 = new TH2F("fTPCnsigma_TOFnsigma3","TOF Electron N#sigma;TPC Electron N#sigma",200,-10,30,200,-15,10);
	fOutputList->Add(fTPCnsigma_TOFnsigma3);
	
	
	///*****************************************************************

	fInvMass = new TH1F("fInvMass","",200,0,0.3);
	fOutputList->Add(fInvMass);
	
	fInvMass1 = new TH1F("fInvMass1","",200,0,0.3);
	fOutputList->Add(fInvMass1);
	
	fInvMass2 = new TH1F("fInvMass2","",200,0,0.3);
	fOutputList->Add(fInvMass2);
	
	fInvMass3 = new TH1F("fInvMass3","",200,0,0.3);
	fOutputList->Add(fInvMass3);
	
	fInvMass4 = new TH1F("fInvMass4","",200,0,0.3);
	fOutputList->Add(fInvMass4);
	
	fInvMassBack = new TH1F("fInvMassBack","",200,0,0.3);
	fOutputList->Add(fInvMassBack);		
	
	fInvMassBack1 = new TH1F("fInvMassBack1","",200,0,0.3);
	fOutputList->Add(fInvMassBack1);
	
	fInvMassBack2 = new TH1F("fInvMassBack2","",200,0,0.3);
	fOutputList->Add(fInvMassBack2);
	
	fInvMassBack3 = new TH1F("fInvMassBack3","",200,0,0.3);
	fOutputList->Add(fInvMassBack3);
	
	fInvMassBack4 = new TH1F("fInvMassBack4","",200,0,0.3);
	fOutputList->Add(fInvMassBack4);
	
	fDCA = new TH1F("fDCA","",5000,0,10);
	fOutputList->Add(fDCA);
	
    fDCABack = new TH1F("fDCABack","",5000,0,10);
	fOutputList->Add(fDCABack);
	
	fAngle = new TH1F("fAngle","",5000,0,5);
	fOutputList->Add(fAngle);
	
    fAngleBack = new TH1F("fAngleBack","",5000,0,5);
	fOutputList->Add(fAngleBack);
	
	fPtElec_ULS = new TH1F("fPtElec_ULS","ULS; p_{T} [GeV/c]; Count",300,0.1,30.1);
	fOutputList->Add(fPtElec_ULS);
	
	fPtElec_LS = new TH1F("fPtElec_LS","LS; p_{T} [GeV/c]; Count",300,0.1,30.1);
    fOutputList->Add(fPtElec_LS);
    
    fTPCnSigmaPElec_ULS = new TH2F("fTPCnSigmaPElec_ULS","p (GeV/c);TPC Electron N#sigma",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnSigmaPElec_ULS);   
	
	fTPCnSigmaPElec_LS = new TH2F("fTPCnSigmaPElec_LS","p (GeV/c);TPC Electron N#sigma",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnSigmaPElec_LS);   
    
    fPtElec = new TH1F("fPtElec","; p_{T} [GeV/c]; Count",300,0,30);
    fOutputList->Add(fPtElec);
    
    fDCAxy_pt = new TH2F("fDCAxy_pt",";p_{t} (GeV/c);DCAxy ",300,0,30,200,-0.2,0.2);
    fOutputList->Add(fDCAxy_pt);
    
    fDCAz_pt = new TH2F("fDCAz_pt",";p_{t} (GeV/c);DCAz ",300,0,30,200,-0.2,0.2);
    fOutputList->Add(fDCAz_pt);
	
	fDCAxy_pt_gamma = new TH2F("fDCAxy_pt_gamma",";p_{t} (GeV/c);DCAxy ",300,0,30,200,-0.2,0.2);
    fOutputList->Add(fDCAxy_pt_gamma);
    
    fDCAz_pt_gamma = new TH2F("fDCAz_pt_gamma",";p_{t} (GeV/c);DCAz ",300,0,30,200,-0.2,0.2);
    fOutputList->Add(fDCAz_pt_gamma);
	
	fDCAxy_pt_light = new TH2F("fDCAxy_pt_light",";p_{t} (GeV/c);DCAxy ",300,0,30,200,-0.2,0.2);
    fOutputList->Add(fDCAxy_pt_light);
    
    fDCAz_pt_light = new TH2F("fDCAz_pt_light",";p_{t} (GeV/c);DCAz ",300,0,30,200,-0.2,0.2);
    fOutputList->Add(fDCAz_pt_light);
    
    fDCAxy_pt_B = new TH2F("fDCAxy_pt_B",";p_{t} (GeV/c);DCAxy ",300,0,30,200,-0.2,0.2);
    fOutputList->Add(fDCAxy_pt_B);
    
    fDCAz_pt_B = new TH2F("fDCAz_pt_B",";p_{t} (GeV/c);DCAz ",300,0,30,200,-0.2,0.2);
    fOutputList->Add(fDCAz_pt_B);
    
    fDCAxy_pt_C = new TH2F("fDCAxy_pt_C",";p_{t} (GeV/c);DCAxy ",300,0,30,200,-0.2,0.2);
    fOutputList->Add(fDCAxy_pt_C);
    
    fDCAz_pt_C = new TH2F("fDCAz_pt_C",";p_{t} (GeV/c);DCAz ",300,0,30,200,-0.2,0.2);
    fOutputList->Add(fDCAz_pt_C);
    
    fPt_elec_phot = new TH1F("fPt_elec_phot","",300,0.1,30.1);
	fOutputList->Add(fPt_elec_phot);
	
	fPt_elec_from_pi0 = new TH1F("fPt_elec_from_pi0","",300,0.1,30.1);
	fOutputList->Add(fPt_elec_from_pi0);
	
	fPt_elec_from_eta = new TH1F("fPt_elec_from_eta","",300,0.1,30.1);
	fOutputList->Add(fPt_elec_from_eta);
	
	fPt_elec_from_gamma = new TH1F("fPt_elec_from_gamma","",300,0.1,30.1);
	fOutputList->Add(fPt_elec_from_gamma);
	
	fPt_HFelec_passTOFandTPC_2 = new TH1F("fPt_HFelec_passTOFandTPC_2","",300,0.1,30.1);
	fOutputList->Add(fPt_HFelec_passTOFandTPC_2);
	
	fPt_HFelec_passTPC_2 = new TH1F("fPt_HFelec_passTPC_2","",300,0.1,30.1);
	fOutputList->Add(fPt_HFelec_passTPC_2);
	
	fPt_HFelec_pass_track = new TH1F("fPt_HFelec_pass_track","",300,0.1,30.1);
	fOutputList->Add(fPt_HFelec_pass_track);
	
	fPt_HFelec_pass_trackDCA = new TH1F("fPt_HFelec_pass_trackDCA","",300,0.1,30.1);
	fOutputList->Add(fPt_HFelec_pass_trackDCA);
	
	fPt_HFelec_pass_trackITSlayer = new TH1F("fPt_HFelec_pass_trackITSlayer","",300,0.1,30.1);
	fOutputList->Add(fPt_HFelec_pass_trackITSlayer);
	
	fPt_HFelec_pass_trackFB = new TH1F("fPt_HFelec_pass_trackFB","",300,0.1,30.1);
	fOutputList->Add(fPt_HFelec_pass_trackFB);
	
	fPt_HFelec_pass_trackTPCITS = new TH1F("fPt_HFelec_pass_trackTPCITS","",300,0.1,30.1);
	fOutputList->Add(fPt_HFelec_pass_trackTPCITS);
	
	fPt_HFelec_pass_trackKink = new TH1F("fPt_HFelec_pass_trackKink","",300,0.1,30.1);
	fOutputList->Add(fPt_HFelec_pass_trackKink);
	
	
	fPt_HFelec_pass_trackandTOF = new TH1F("fPt_HFelec_pass_trackandTOF","",300,0.1,30.1);
	fOutputList->Add(fPt_HFelec_pass_trackandTOF);
	
	
	fPtMCpi0 = new TH1F("fPtMCpi0",";p_{t} (GeV/c)",2000,0,100);
	fOutputList->Add(fPtMCpi0);
	
	fPtMCpi0_feeddown = new TH1F("fPtMCpi0_feeddown",";p_{t} (GeV/c)",2000,0,100);
	fOutputList->Add(fPtMCpi0_feeddown);
	
	fPtMCpi0_nomother = new TH1F("fPtMCpi0_nomother",";p_{t} (GeV/c)",2000,0,100);
	fOutputList->Add(fPtMCpi0_nomother);
	
	
	fPtMCeta = new TH1F("fPtMCeta",";p_{t} (GeV/c)",2000,0,100);
	fOutputList->Add(fPtMCeta);
	
	fPtMCeta_feeddown = new TH1F("fPtMCeta_feeddown",";p_{t} (GeV/c)",2000,0,100);
	fOutputList->Add(fPtMCeta_feeddown);
	
	fPtMCeta_nomother = new TH1F("fPtMCeta_nomother",";p_{t} (GeV/c)",2000,0,100);
	fOutputList->Add(fPtMCeta_nomother);
	
	fPtMCpi0_afterNorm = new TH1F("fPtMCpi0_afterNorm",";p_{t} (GeV/c)",2000,0,100);
	fOutputList->Add(fPtMCpi0_afterNorm);
	
	fPtMCeta_afterNorm = new TH1F("fPtMCeta_afterNorm",";p_{t} (GeV/c)",2000,0,100);
	fOutputList->Add(fPtMCeta_afterNorm);
	
	fPtMCpi02 = new TH1F("fPtMCpi02",";p_{t} (GeV/c)",2000,0,100);
	fOutputList->Add(fPtMCpi02);
	
	fPtMCeta2 = new TH1F("fPtMCeta2",";p_{t} (GeV/c)",2000,0,100);
	fOutputList->Add(fPtMCeta2);
	
	fPtMCpi03 = new TH1F("fPtMCpi03",";p_{t} (GeV/c)",2000,0,100);
	fOutputList->Add(fPtMCpi03);
	
	fPtMCeta3 = new TH1F("fPtMCeta3",";p_{t} (GeV/c)",2000,0,100);
	fOutputList->Add(fPtMCeta3);
	
	fPtElec_ULS_MC = new TH1F("fPtElec_ULS_MC","p_{T} [GeV/c]",300,0.1,30.1);
	fOutputList->Add(fPtElec_ULS_MC);
	
	fPtElec_LS_MC = new TH1F("fPtElec_LS_MC","p_{T} [GeV/c]",300,0.1,30.1);
    fOutputList->Add(fPtElec_LS_MC);
	
	fPtElec_ULS_MC_from_pi0 = new TH1F("fPtElec_ULS_MC_from_pi0","p_{T} [GeV/c]",300,0.1,30.1);
	fOutputList->Add(fPtElec_ULS_MC_from_pi0);
	
	fPtElec_LS_MC_from_pi0 = new TH1F("fPtElec_LS_MC_from_pi0","p_{T} [GeV/c]",300,0.1,30.1);
    fOutputList->Add(fPtElec_LS_MC_from_pi0);
	
	fPtElec_ULS_MC_from_eta = new TH1F("fPtElec_ULS_MC_from_eta","p_{T} [GeV/c]",300,0.1,30.1);
	fOutputList->Add(fPtElec_ULS_MC_from_eta);
	
	fPtElec_LS_MC_from_eta = new TH1F("fPtElec_LS_MC_from_eta","p_{T} [GeV/c]",300,0.1,30.1);
    fOutputList->Add(fPtElec_LS_MC_from_eta);
	
	fPtElec_ULS_MC_from_gamma = new TH1F("fPtElec_ULS_MC_from_gamma","p_{T} [GeV/c]",300,0.1,30.1);
	fOutputList->Add(fPtElec_ULS_MC_from_gamma);
	
	fPtElec_LS_MC_from_gamma = new TH1F("fPtElec_LS_MC_from_gamma","p_{T} [GeV/c]",300,0.1,30.1);
    fOutputList->Add(fPtElec_LS_MC_from_gamma);
	
	
	fPtHFEMC_2 = new TH1F("fPtHFEMC_2",";p_{t} (GeV/c)",300,0.1,30.1);
	fOutputList->Add(fPtHFEMC_2);
	
	fTPCnsigma_p_after_tof = new TH2F("fTPCnsigma_p_after_tof","p (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnsigma_p_after_tof);
	
	fTPCnsigma_p_after_tof_its = new TH2F("fTPCnsigma_p_after_tof_its","p (GeV/c);TPC Electron N#sigma after TOF and ITS cuts",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnsigma_p_after_tof_its);
	
	fTPCnsigma_p_after_tof_inverse = new TH2F("fTPCnsigma_p_after_tof_inverse","p (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnsigma_p_after_tof_inverse);
	
	fTPCnsigma_pt_after_tof = new TH2F("fTPCnsigma_pt_after_tof","pt (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnsigma_pt_after_tof);
    
    fTPCnsigma_pt_after_tof_its = new TH2F("fTPCnsigma_pt_after_tof_its","pt (GeV/c);TPC Electron N#sigma after TOF and ITS cuts",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnsigma_pt_after_tof_its);
    
    fTPCnsigma_pt_after_tof_inverse = new TH2F("fTPCnsigma_pt_after_tof_inverse","pt (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnsigma_pt_after_tof_inverse);
	
//______________________________________________________________________
	
	PostData(1, fOutputList);
	
///______________________________________________________________________
}

	Int_t pdg = -99999;
	Int_t pdg_mother = -99999;
	Int_t pdg_gmother = -99999;
	Double_t pt_gmother = -99999;
	Double_t pt_mother = -99999;
	Double_t weight = -99999;
	
	Double_t fTPCnSigma = -999;
	Double_t fTOFnSigma = -999;
	Double_t fITSnSigma = -999;
	Double_t fTPCnSigma_pion = -999;
	Double_t fTPCnSigma_proton = -999;
	Double_t fTPCnSigma_kaon = -999;
	Double_t fTPCsignal = -999;
	Double_t fPt = -999;
	Double_t fEta = -999;
	Double_t fPhi = -999;
	Double_t fP = -999;


//______________________________________________________________________
//Main loop
//Called for each event
void AliAnalysisHFETPCTOF::UserExec(Option_t *) 
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
	
	
//Select events which the reconstruction DOES include SDD
	fNevent->Fill(0);
	if(fIsPP){
		if(!fIsMC){ 
		
			TString firedTriggerClasses = static_cast<const AliAODEvent*>(InputEvent())->GetFiredTriggerClasses();
			if (!(firedTriggerClasses.Contains("ALLNOTRD"))){
				fNevent2->Fill(0);
				return; // events without SDD
			}
			//firedTriggerClasses.Contains("ALLNOTRD")
		
		}
	}
	//if(fIsPP){
	//	cout<<"IS PP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	//}
		
	//TString trcl = fVevent->GetFiredTriggerClasses();	
	//TString trcl = esdEvent->GetFiredTriggerClasses();
	//if (!(trcl->Contains("ALLNOTRD"))) return; // events without SDD

//______________________________________________________________________
//Vertex Selection for pPb
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
				
		const AliESDVertex *trkVtx = fESD->GetPrimaryVertex();
		Float_t zvtx = trkVtx->GetZ();
		fZvtx = zvtx;
		if(TMath::Abs(zvtx) > 10) return;
		
		
	}

//______________________________________________________________________	
	
//Only events with at least 2 tracks are accepted
	Int_t fNOtrks =  fVevent->GetNumberOfTracks();
	if(fNOtrks<2) return;
//______________________________________________________________________
	
//Centrality selection
	if(!fIsPP){
		if(!fIsMC){     
			Float_t lPercentile = 300; 
			AliMultSelection *MultSelection = 0x0; 
			MultSelection = (AliMultSelection * ) fVevent->FindListObject("MultSelection");
			if( !MultSelection) {
				//If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
				AliWarning("AliMultSelection object not found!");
			}else{
				lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
			}  
    
			Double_t cent = MultSelection->GetMultiplicityPercentile("V0M", true); // returns centrality only for events used in calibration
			//printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!cent = %f\n", cent);
			//if(cent>=0.0 && cent<=10.0){
			//    printf("******************************************************************************cent = %f\n", cent);
			//}	
    
			//Centrality
			fCent->Fill(cent);
    
			if(cent<0.5 || cent>10.5) return;
		}
	}
//______________________________________________________________________
    
    

	
	Int_t ClsNo = -999;
	if(!fIsAOD) ClsNo = fESD->GetNumberOfCaloClusters(); 
	else ClsNo = fAOD->GetNumberOfCaloClusters(); 
	
//=======================================================================
//=======================================================================
///Initialization for MC analysis:

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
			 
			 //pt spectra for pi0 and eta
			 //*********************************************************
			 
			 weight = 1;
			 for(Int_t iMC = 0; iMC < fMCarray->GetEntries(); iMC++)
			 {
				 fMCparticle = (AliAODMCParticle*) fMCarray->At(iMC);
				 pdg = fMCparticle->GetPdgCode();
				 
				     if(fMCparticle->Eta()>=-0.8 && fMCparticle->Eta()<=0.8)
				     {
						 
						// if(fMCparticle->IsPrimary()){
							 if(TMath::Abs(pdg) == 11){
								Bool_t MotherFound = FindMother(iMC);
								if(MotherFound){
									if(fIsHFE1){
										//cout<<"00000000000000000000000000000000000000000000pdg_mother = "<<pdg_mother<<endl;
										fPtHFEMC_2->Fill(fMCparticle->Pt());
									}
								}
					       
							}
						// }
						 
									 

			            if(fMCparticle->IsPrimary()){ ///Does not include particles from weak decays or created in an interaction with the material 

							if(TMath::Abs(pdg)==111) fPtMCpi03->Fill(fMCparticle->Pt());
							if(TMath::Abs(pdg)==221) fPtMCeta3->Fill(fMCparticle->Pt());
							
							
							if(TMath::Abs(pdg)==111){ 
								if(fMCparticle->GetMother()>0){ 
									fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
									pdg_mother = TMath::Abs(fMCparticleMother->GetPdgCode());
									///if mother is a pi0, eta, omega, phi, etaprim, or rho
									if(pdg_mother == 111 || pdg_mother == 221 || pdg_mother == 223 || pdg_mother == 333 || pdg_mother == 331 || pdg_mother == 113){
										fPtMCpi0_feeddown->Fill(fMCparticle->Pt());
									}
									else{
										fPtMCpi0->Fill(fMCparticle->Pt());
									}
									//if(TMath::Abs(pdg)==221) fPtMCeta->Fill(fMCparticle->Pt());
								}
								else{
									///enhancement
									fPtMCpi0_nomother->Fill(fMCparticle->Pt());

									weight = CalculateWeight(111,fMCparticle->Pt());
									fPtMCpi0_afterNorm->Fill(fMCparticle->Pt(),weight);
 				        
									
								}
							}
								
							if(TMath::Abs(pdg)==221){ 
								if(fMCparticle->GetMother()>0){ 
									fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
									pdg_mother = TMath::Abs(fMCparticleMother->GetPdgCode());
									///if mother is a pi0, eta, omega, phi, etaprim, or rho
									if(pdg_mother == 111 || pdg_mother == 221 || pdg_mother == 223 || pdg_mother == 333 || pdg_mother == 331 || pdg_mother == 113){
										fPtMCeta_feeddown->Fill(fMCparticle->Pt());
									}
									else{
										fPtMCeta->Fill(fMCparticle->Pt());
									}
									//if(TMath::Abs(pdg)==221) fPtMCeta->Fill(fMCparticle->Pt());
								}
								else{
									///enhancement
									fPtMCeta_nomother->Fill(fMCparticle->Pt());
									
									weight = CalculateWeight(221,fMCparticle->Pt());
									fPtMCeta_afterNorm->Fill(fMCparticle->Pt(),weight);
								}
							}	
								
								
						} //End of IsPrimary
						
					    if(fMCparticle->IsPhysicalPrimary()){
						    if(TMath::Abs(pdg)==111) fPtMCpi02->Fill(fMCparticle->Pt());
						    if(TMath::Abs(pdg)==221) fPtMCeta2->Fill(fMCparticle->Pt());
						
					    }
								
					   
					
					   

			       }//end of eta cut
			  // }
		   }
		   
			 //*********************************************************
			 
			 
         }    
    }
    
//=======================================================================
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
		
		
		
		if((track->Eta() > -0.8) && (track->Eta() < 0.8)){///********************************************************************************************       


		
		///_____________________________________________________________________________
		///Fill QA plots without track selection
		fPt = track->Pt();
		fEta = track->Eta();
		fPhi = track->Phi();
		fP = TMath::Sqrt((track->Pt())*(track->Pt()) + (track->Pz())*(track->Pz()));
		
		fTPCsignal = track->GetTPCsignal();
		fTPCnSigma = fPidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
		fTOFnSigma = fPidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
		fITSnSigma = fPidResponse->NumberOfSigmasITS(track, AliPID::kElectron);
		fTPCnSigma_pion = fPidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
		fTPCnSigma_proton = fPidResponse->NumberOfSigmasTPC(track, AliPID::kProton);
		fTPCnSigma_kaon = fPidResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        
        fTPC_p1->Fill(fP,fTPCsignal);
        fTPCnsigma_p1->Fill(fP,fTPCnSigma);
        fTPCnsigma_pt1->Fill(fPt,fTPCnSigma);
        Float_t pos2[3]={0,0,0};
        
        fTPCnsigma_TOFnsigma1->Fill(fTOFnSigma,fTPCnSigma);
        fTOFnsigma_p1->Fill(fP,fTOFnSigma);
        fTOFnsigma_pt1->Fill(fPt,fTOFnSigma);
        
        fITSnsigma_p1->Fill(fP,fITSnSigma);
        fITSnsigma_pt1->Fill(fPt,fITSnSigma);
        
        //fVtxZ[0]->Fill(fZvtx);
        //fNClusters[0]->Fill(ClsNo);
        //fEtaPhiTracks[0]->Fill(fEta,fPhi);
        
        
        
        
        
		
//=======================================================================
// Track Selection Cuts are applied here
//=======================================================================


//AOD (Test Filter Bit)---------------------------------------------------------------------
		if(fIsAOD)
		{
						
			if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; 
		}
		
		///****How many HF electrons pass the cuts (eta and filterbit selection)*******
		if(fIsMC && fIsAOD && track->GetLabel()>=0){           	            
			fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
			//if(fMCparticle->IsPrimary()){
				pdg = fMCparticle->GetPdgCode();
				///Is electron:
				if(TMath::Abs(pdg) == 11){
					//printf("pdg = %d\n !!!",pdg);
				
					Bool_t MotherFound = FindMother(track->GetLabel());
					if(MotherFound){
						if(fIsHFE1){
							fPt_HFelec_pass_trackFB->Fill(fPt);
						}
					}
						
					
				}
			//}
		}
		///*********************************************************
//-----------------------------------------------------------------------------------
		
//RecKine: ITSTPC cuts --------------------------------------------------------------
		if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
		
		///****How many HF electrons pass the cuts (eta and TPCITS selection)*******
		if(fIsMC && fIsAOD && track->GetLabel()>=0){           	            
			fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
			//if(fMCparticle->IsPrimary()){
				pdg = fMCparticle->GetPdgCode();
				///Is electron:
				if(TMath::Abs(pdg) == 11){
					//printf("pdg = %d\n !!!",pdg);
				
					Bool_t MotherFound = FindMother(track->GetLabel());
					if(MotherFound){
						if(fIsHFE1){
							fPt_HFelec_pass_trackTPCITS->Fill(fPt);
						}
					}
						
					
				}
			//}
		}
		///*********************************************************
//--------------------------------------------------------------------------------	
		
		
//RecKink -------------------------------------------------------------------------
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
		
		
		///****How many HF electrons pass the cuts (eta and TPCITS selection)*******
		if(fIsMC && fIsAOD && track->GetLabel()>=0){           	            
			fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
			//if(fMCparticle->IsPrimary()){
				pdg = fMCparticle->GetPdgCode();
				///Is electron:
				if(TMath::Abs(pdg) == 11){
					//printf("pdg = %d\n !!!",pdg);
				
					Bool_t MotherFound = FindMother(track->GetLabel());
					if(MotherFound){
						if(fIsHFE1){
							fPt_HFelec_pass_trackKink->Fill(fPt);
						}
					}
						
					
				}
			//}
		}
		///*********************************************************
//---------------------------------------------------------------------------------    
    
//RecPrim -------------------------------------------------------------------------
		//if(!fIsAOD)
		//{
			if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;				///ProcessCutStep(Int_t cutStep, AliVParticle *track)
		//}
		
		///****How many HF electrons pass the cuts (eta and track DCA selection)*******
		if(fIsMC && fIsAOD && track->GetLabel()>=0){           	            
			fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
			//if(fMCparticle->IsPrimary()){
				pdg = fMCparticle->GetPdgCode();
				///Is electron:
				if(TMath::Abs(pdg) == 11){
					//printf("pdg = %d\n !!!",pdg);
				
					Bool_t MotherFound = FindMother(track->GetLabel());
					if(MotherFound){
						if(fIsHFE1){
							fPt_HFelec_pass_trackDCA->Fill(fPt);
						}
					}
						
					
				}
			//}
		}
		///*************************************************************
//-------------------------------------------------------------------------------------		
    
//HFEcuts: ITS layers cuts -------------------------------------------------------------
		if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;
		
		///****How many HF electrons pass the cuts (eta and track ITS layer selection)*******
		if(fIsMC && fIsAOD && track->GetLabel()>=0){           	            
			fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
			//if(fMCparticle->IsPrimary()){
				pdg = fMCparticle->GetPdgCode();
				///Is electron:
				if(TMath::Abs(pdg) == 11){
					//printf("pdg = %d\n !!!",pdg);
				
					Bool_t MotherFound = FindMother(track->GetLabel());
					if(MotherFound){
						if(fIsHFE1){
							fPt_HFelec_pass_trackITSlayer->Fill(fPt);
						}
					}
						
					
				}
			//}
		}
		///*********************************************************
//--------------------------------------------------------------------------------------		
		
		
    
//HFE cuts: TPC PID cleanup
		if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;
        
//=======================================================================
// QA plots after track selection
//=======================================================================

		//if((track->Eta() > -0.8) && (track->Eta() < 0.8)){///********************************************************************************************
		
		
		///****How many HF electrons pass the cuts (eta and track selection)*******
			
		if(fIsMC && fIsAOD && track->GetLabel()>=0){           	            
			fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
			//if(fMCparticle->IsPrimary()){
				pdg = fMCparticle->GetPdgCode();
				///Is electron:
				if(TMath::Abs(pdg) == 11){
					//printf("pdg = %d\n !!!",pdg);
				
					Bool_t MotherFound = FindMother(track->GetLabel());
					if(MotherFound){
						if(fIsHFE1){
							fPt_HFelec_pass_track->Fill(fPt);
						}
					}
						
					
				}
			//}
		}
			
			
			
			///*********************************************************
		
		
		
		
		fTPC_p2->Fill(fP,fTPCsignal);
		
		fTPCnsigma_TOFnsigma2->Fill(fTOFnSigma,fTPCnSigma);
		fTOFnsigma_p2->Fill(fP,fTOFnSigma);
		fTOFnsigma_pt2->Fill(fPt,fTOFnSigma);
		fTPCnsigma_p2->Fill(fP,fTPCnSigma);
		fTPCnsigma_pt2->Fill(fPt,fTPCnSigma);
		fITSnsigma_p2->Fill(fP,fITSnSigma);
		fITSnsigma_pt2->Fill(fPt,fITSnSigma);
			
		
		if(fTOFnSigma >= -3 && fTOFnSigma <= 3){
			fTPCnsigma_p_after_tof->Fill(fP,fTPCnSigma);
			fTPCnsigma_pt_after_tof->Fill(fPt,fTPCnSigma);
			
			if(fITSnSigma >= -2 && fITSnSigma <= 2){
				fTPCnsigma_p_after_tof_its->Fill(fP,fTPCnSigma);
				fTPCnsigma_pt_after_tof_its->Fill(fPt,fTPCnSigma);
			}
			
			///****How many HF electrons pass the cuts (TOF and track selection)*******
			
			if(fIsMC && fIsAOD && track->GetLabel()>=0){           	            
				fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
				//if(fMCparticle->IsPrimary()){
				pdg = fMCparticle->GetPdgCode();
				///Is electron:
				if(TMath::Abs(pdg) == 11){
					//printf("pdg = %d\n !!!",pdg);
					
					Bool_t MotherFound = FindMother(track->GetLabel());
					if(MotherFound){
						if(fIsHFE1){
							fPt_HFelec_pass_trackandTOF->Fill(fPt);
						}
					}
						
					
				}
			//}
			}
		///*********************************************************
			
			
		}
		if(fTOFnSigma <= -3 || fTOFnSigma >= 3){
			fTPCnsigma_p_after_tof_inverse->Fill(fP,fTPCnSigma);
			fTPCnsigma_pt_after_tof_inverse->Fill(fPt,fTPCnSigma);
		}
		
				
		//fVtxZ[1]->Fill(fZvtx);
		//fNClusters[1]->Fill(ClsNo);
		//fEtaPhiTracks[1]->Fill(fEta,fPhi);
		
		
		
		
		

			
//=======================================================================
// Here the PID cuts defined in the file "ConfigEMCalHFEpA.C" is applied
//=======================================================================
	    Int_t pidpassed = 1;
	    AliHFEpidObject hfetrack;
	    hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
	    hfetrack.SetRecTrack(track);
	    hfetrack.SetPP();	//proton-proton analysis
	    if(!fPID->IsSelected(&hfetrack, NULL, "", fPIDqa)) pidpassed = 0;
	   		
	    if(pidpassed==0) continue;
        
//=======================================================================
// QA plots after PID selection
//=======================================================================
		
		
	///****How many HF electrons pass the cuts (TPC and track selection)*******
			
		if(fIsMC && fIsAOD && track->GetLabel()>=0){           	            
			fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
			//if(fMCparticle->IsPrimary()){
				pdg = fMCparticle->GetPdgCode();
				///Is electron:
				if(TMath::Abs(pdg) == 11){
					//printf("pdg = %d\n !!!",pdg);
					
					Bool_t MotherFound = FindMother(track->GetLabel());
					if(MotherFound){
						if(fIsHFE1){
							fPt_HFelec_passTPC_2->Fill(fPt);
						}
					}
						
					
				}
			//}
		}
	///*********************************************************
	
	
	    if(fTOFnSigma >= -3 && fTOFnSigma <= 3){
					
			fTPC_p3->Fill(fP,fTPCsignal);
			fTPCnsigma_p3->Fill(fP,fTPCnSigma);
			fTPCnsigma_pt3->Fill(fPt,fTPCnSigma);
		
			fTPCnsigma_TOFnsigma3->Fill(fTOFnSigma,fTPCnSigma);
			fTOFnsigma_p3->Fill(fP,fTOFnSigma);
			fTOFnsigma_pt3->Fill(fPt,fTOFnSigma);
			
			fITSnsigma_p3->Fill(fP,fITSnSigma);
			fITSnsigma_pt3->Fill(fPt,fITSnSigma);
			
			fPtElec->Fill(fPt);
			
			///****How many HF electrons pass the cuts (TOF, TPC and track selection)*******
			
			if(fIsMC && fIsAOD && track->GetLabel()>=0){           	            
				fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
				//if(fMCparticle->IsPrimary()){
					pdg = fMCparticle->GetPdgCode();
					///Is electron:
					if(TMath::Abs(pdg) == 11){
						//printf("pdg = %d\n !!!",pdg);
					
						Bool_t MotherFound = FindMother(track->GetLabel());
						if(MotherFound){
							if(fIsHFE1){
								fPt_HFelec_passTOFandTPC_2->Fill(fPt);
							}
						}
						
					
					}
				//}
			}
			
			
			
			///*********************************************************
			
			
        
			///InvMass calculation (to find NonHFE):
			///**********************************************************
			//Associated particle cut
		            
			fNonHFE = new AliSelectNonHFE();
		            
			fPartnerCuts->SetAcceptKinkDaughters(kFALSE);
			fPartnerCuts->SetRequireITSRefit(kTRUE);
			fPartnerCuts->SetRequireTPCRefit(kTRUE);
			fPartnerCuts->SetEtaRange(-0.9,0.9);
			fPartnerCuts->SetMaxChi2PerClusterTPC(4.0);
			fPartnerCuts->SetMinNClustersTPC(80);
			fPartnerCuts->SetPtRange(0.3,1e10);

			fNonHFE->SetAODanalysis(fIsAOD);

			if(fMassCutFlag) fNonHFE->SetInvariantMassCut(fMassCut);
			if(fAngleCutFlag) fNonHFE->SetOpeningAngleCut(fAngleCut);
			if(fChi2CutFlag) fNonHFE->SetChi2OverNDFCut(fChi2Cut);
			if(fDCAcutFlag) fNonHFE->SetDCACut(fDCAcut);
			fNonHFE->SetAlgorithm("DCA"); //KF
			fNonHFE->SetPIDresponse(fPidResponse);
			fNonHFE->SetTrackCuts(-3.5,3.5,fPartnerCuts);
			fNonHFE->SetAdditionalCuts(fPtMinAsso,fTpcNclsAsso);
										
										
			
			fNonHFE->SetHistDCABack(fDCABack);
			fNonHFE->SetHistDCA(fDCA);
			fNonHFE->SetHistAngleBack(fAngleBack);
			fNonHFE->SetHistAngle(fAngle);
			           
			if(fPt>0 && fPt<=1){
				fNonHFE->SetHistMassBack(fInvMassBack1);
				fNonHFE->SetHistMass(fInvMass1);
			}           
			if(fPt>1 && fPt<=2){
				fNonHFE->SetHistMassBack(fInvMassBack2);
				fNonHFE->SetHistMass(fInvMass2);
			}
			if(fPt>2 && fPt<=3){
				fNonHFE->SetHistMassBack(fInvMassBack3);
				fNonHFE->SetHistMass(fInvMass3);
			}
			if(fPt>3 && fPt<=4){
				fNonHFE->SetHistMassBack(fInvMassBack4);
				fNonHFE->SetHistMass(fInvMass4);
			}
			if(fPt>4){
				fNonHFE->SetHistMassBack(fInvMassBack);
				fNonHFE->SetHistMass(fInvMass);
			}
			       
			           
			fNonHFE->FindNonHFE(iTracks,Vtrack,fVevent); //Calling the function that will fill the histograms above
			//(Int_t iTrack1, AliVParticle *Vtrack1, AliVEvent *fVevent)
			        
			
			           
			if(fNonHFE->IsULS()) fPtElec_ULS->Fill(fPt,fNonHFE->GetNULS());
			if(fNonHFE->IsLS()) fPtElec_LS->Fill(fPt,fNonHFE->GetNLS());
					
			if(fNonHFE->IsULS()) fTPCnSigmaPElec_ULS->Fill(fP,fTPCnSigma*fNonHFE->GetNULS());
			if(fNonHFE->IsLS()) fTPCnSigmaPElec_LS->Fill(fP,fTPCnSigma*fNonHFE->GetNLS());
					
			///*************************************************************
        
            ///*************Invariant mass efficiency (tagging efficiency)********************
            
            if(fIsMC && fIsAOD && track->GetLabel()>=0){      
					            
				fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
				//if(fMCparticle->IsPrimary()){
					pdg = fMCparticle->GetPdgCode();
					///Is electron:
					if(TMath::Abs(pdg) == 11){
						//printf("pdg = %d\n !!!",pdg);
						if(fMCparticle->GetMother()>0){
							fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
							pdg_mother = fMCparticleMother->GetPdgCode();
							pt_mother = fMCparticleMother->Pt();
							
							
							
							//printf("pdg_mother = %d\n ",pdg_mother);
							weight=1;
						
							
							///Photonic Electrons:
							if(TMath::Abs(pdg_mother) == 22 || TMath::Abs(pdg_mother) == 111 || TMath::Abs(pdg_mother) == 221){
							
								///Electrons from pi0 (Dalitz):
								if(TMath::Abs(pdg_mother) == 111){
									///The ones used to calculate the weight
									if(fMCparticleMother->IsPrimary()){
										
											if(fMCparticleMother->GetMother() <= 0){
							   
												weight = CalculateWeight(pdg_mother,pt_mother);
									   
												fPt_elec_phot->Fill(fPt,weight);  
												if(fNonHFE->IsULS()) fPtElec_ULS_MC->Fill(fPt,fNonHFE->GetNULS()*weight);
												if(fNonHFE->IsLS()) fPtElec_LS_MC->Fill(fPt,fNonHFE->GetNLS()*weight);
			
												fPt_elec_from_pi0->Fill(fPt,weight);
												if(fNonHFE->IsULS()) fPtElec_ULS_MC_from_pi0->Fill(fPt,fNonHFE->GetNULS()*weight);
												if(fNonHFE->IsLS()) fPtElec_LS_MC_from_pi0->Fill(fPt,fNonHFE->GetNLS()*weight);
										
										
										
											}
										
										
										
									}
								}
								///Electrons from eta (Dalitz)
								if(TMath::Abs(pdg_mother) == 221){
									///The ones used to calculate the weight
									if(fMCparticleMother->IsPrimary()){
										
										
												if(fMCparticleMother->GetMother() <= 0){   
					
												weight = CalculateWeight(pdg_mother,pt_mother);
									   
												fPt_elec_phot->Fill(fPt,weight); 
												if(fNonHFE->IsULS()) fPtElec_ULS_MC->Fill(fPt,fNonHFE->GetNULS()*weight);
												if(fNonHFE->IsLS()) fPtElec_LS_MC->Fill(fPt,fNonHFE->GetNLS()*weight);
								
												fPt_elec_from_eta->Fill(fPt,weight);
												if(fNonHFE->IsULS()) fPtElec_ULS_MC_from_eta->Fill(fPt,fNonHFE->GetNULS()*weight);
												if(fNonHFE->IsLS()) fPtElec_LS_MC_from_eta->Fill(fPt,fNonHFE->GetNLS()*weight);
											}
										
																				
									}
			
								}
								///Electron from gamma
								if(TMath::Abs(pdg_mother) == 22){
									///when gamma has a mother:
									if(fMCparticleMother->GetMother()>0){
										fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
										pdg_gmother = fMCparticleGMother->GetPdgCode();
										pt_gmother = fMCparticleGMother->Pt();
								
										///And the mother is a pi0 or an eta
										if(TMath::Abs(pdg_gmother) == 111 || TMath::Abs(pdg_gmother) == 221){
											
											if(fMCparticleGMother->IsPrimary()){
												//if(fIsPP){
													if(fMCparticleGMother->GetMother() <= 0){ 

														weight = CalculateWeight(pdg_gmother,pt_gmother);
	
														fPt_elec_phot->Fill(fPt,weight);
														if(fNonHFE->IsULS()) fPtElec_ULS_MC->Fill(fPt,fNonHFE->GetNULS()*weight);
														if(fNonHFE->IsLS()) fPtElec_LS_MC->Fill(fPt,fNonHFE->GetNLS()*weight);
									
														fPt_elec_from_gamma->Fill(fPt,weight);
														if(fNonHFE->IsULS()) fPtElec_ULS_MC_from_gamma->Fill(fPt,fNonHFE->GetNULS()*weight);
														if(fNonHFE->IsLS()) fPtElec_LS_MC_from_gamma->Fill(fPt,fNonHFE->GetNLS()*weight);
													}
												//}
												
												
												
											}
	
										}
										///And the mother is not a pi0 or an eta:
										/*
										else{		   
											fPt_elec_phot->Fill(fPt);
											if(fNonHFE->IsULS()) fPtElec_ULS_MC->Fill(fPt,fNonHFE->GetNULS());
											if(fNonHFE->IsLS()) fPtElec_LS_MC->Fill(fPt,fNonHFE->GetNLS());
										
											fPt_elec_from_gamma->Fill(fPt);
											if(fNonHFE->IsULS()) fPtElec_ULS_MC_from_gamma->Fill(fPt,fNonHFE->GetNULS());
											if(fNonHFE->IsLS()) fPtElec_LS_MC_from_gamma->Fill(fPt,fNonHFE->GetNLS());
										}
										*/
									}
									///electron from gamma, when the gamma does not have a mother (direct gamma):	  
									/*
									else{		   
										fPt_elec_phot->Fill(fPt);
										if(fNonHFE->IsULS()) fPtElec_ULS_MC->Fill(fPt,fNonHFE->GetNULS());
										if(fNonHFE->IsLS()) fPtElec_LS_MC->Fill(fPt,fNonHFE->GetNLS());
									
										fPt_elec_from_gamma->Fill(fPt);
										if(fNonHFE->IsULS()) fPtElec_ULS_MC_from_gamma->Fill(fPt,fNonHFE->GetNULS());
										if(fNonHFE->IsLS()) fPtElec_LS_MC_from_gamma->Fill(fPt,fNonHFE->GetNLS());
									}
									*/
								
								}
														
							      
							} ///End of photonic electrons
		                
					    
					    
					    
		                
						
						}
					}
				//}
				
					
			}
            
            ///**********************************************************
            
        
		
			//fVtxZ[2]->Fill(fZvtx);
			//fNClusters[2]->Fill(ClsNo);
			//fEtaPhiTracks[2]->Fill(fEta,fPhi);
		
	    }
		
		
					/*
					Double_t d0z0[2], cov[3];
		            AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
		            track->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 20., d0z0, cov);
		            Double_t DCAxy = d0z0[0];
		            Double_t DCAz = d0z0[1];
		            
		            fDCAxy_pt->Fill(fPt,DCAxy);
		            fDCAz_pt->Fill(fPt,DCAz);
		            */
		            
							            
					
        
		} //End of eta selection
        
		
        }//End of track loop
	
//=======================================================================
	
	delete fListOfmotherkink;
	PostData(1, fOutputList);
}      

//=======================================================================
void AliAnalysisHFETPCTOF::Terminate(Option_t *) 
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
Double_t AliAnalysisHFETPCTOF::CalculateWeight(Int_t pdg_particle, Double_t x)
{
     Double_t weight=1;
     
     if(fIsPP){
		if(pdg_particle==111){ //pi0


///Weight for only enhanced spectra (choosing IsPrimary in both MC)

if(x>=0 && x<0.1) weight = 482.74;
if(x>=0.1 && x<0.2) weight = 774.008;
if(x>=0.2 && x<0.3) weight = 772.56;
if(x>=0.3 && x<0.4) weight = 641.061;
if(x>=0.4 && x<0.5) weight = 498.833;
if(x>=0.5 && x<0.6) weight = 383.013;
if(x>=0.6 && x<0.7) weight = 286.518;
if(x>=0.7 && x<0.8) weight = 217.158;
if(x>=0.8 && x<0.9) weight = 165.97;
if(x>=0.9 && x<1) weight = 128.828;
if(x>=1 && x<1.1) weight = 99.3552;
if(x>=1.1 && x<1.2) weight = 78.2052;
if(x>=1.2 && x<1.3) weight = 62.2902;
if(x>=1.3 && x<1.4) weight = 49.6331;
if(x>=1.4 && x<1.5) weight = 40.2988;
if(x>=1.5 && x<1.6) weight = 32.6421;
if(x>=1.6 && x<1.7) weight = 26.6991;
if(x>=1.7 && x<1.8) weight = 22.3921;
if(x>=1.8 && x<1.9) weight = 18.5378;
if(x>=1.9 && x<2) weight = 15.3481;
if(x>=2 && x<2.1) weight = 12.9108;
if(x>=2.1 && x<2.2) weight = 10.8542;
if(x>=2.2 && x<2.3) weight = 9.27531;
if(x>=2.3 && x<2.4) weight = 7.88082;
if(x>=2.4 && x<2.5) weight = 6.75634;
if(x>=2.5 && x<2.6) weight = 5.77543;
if(x>=2.6 && x<2.7) weight = 4.99225;
if(x>=2.7 && x<2.8) weight = 4.33398;
if(x>=2.8 && x<2.9) weight = 3.75435;
if(x>=2.9 && x<3) weight = 3.29254;
if(x>=3 && x<3.1) weight = 2.86516;
if(x>=3.1 && x<3.2) weight = 2.51818;
if(x>=3.2 && x<3.3) weight = 2.22884;
if(x>=3.3 && x<3.4) weight = 1.96287;
if(x>=3.4 && x<3.5) weight = 1.74215;
if(x>=3.5 && x<3.6) weight = 1.52082;
if(x>=3.6 && x<3.7) weight = 1.37077;
if(x>=3.7 && x<3.8) weight = 1.21582;
if(x>=3.8 && x<3.9) weight = 1.10746;
if(x>=3.9 && x<4) weight = 0.971212;
if(x>=4 && x<4.1) weight = 0.884194;
if(x>=4.1 && x<4.2) weight = 0.792149;
if(x>=4.2 && x<4.3) weight = 0.713165;
if(x>=4.3 && x<4.4) weight = 0.646046;
if(x>=4.4 && x<4.5) weight = 0.583767;
if(x>=4.5 && x<4.6) weight = 0.525869;
if(x>=4.6 && x<4.7) weight = 0.483417;
if(x>=4.7 && x<4.8) weight = 0.439922;
if(x>=4.8 && x<4.9) weight = 0.395916;
if(x>=4.9 && x<5) weight = 0.368723;
if(x>=5 && x<5.1) weight = 0.334619;
if(x>=5.1 && x<5.2) weight = 0.304714;
if(x>=5.2 && x<5.3) weight = 0.277995;
if(x>=5.3 && x<5.4) weight = 0.255232;
if(x>=5.4 && x<5.5) weight = 0.236205;
if(x>=5.5 && x<5.6) weight = 0.220869;
if(x>=5.6 && x<5.7) weight = 0.203044;
if(x>=5.7 && x<5.8) weight = 0.18389;
if(x>=5.8 && x<5.9) weight = 0.170851;
if(x>=5.9 && x<6) weight = 0.158248;
if(x>=6 && x<6.1) weight = 0.143708;
if(x>=6.1 && x<6.2) weight = 0.134629;
if(x>=6.2 && x<6.3) weight = 0.123699;
if(x>=6.3 && x<6.4) weight = 0.115016;
if(x>=6.4 && x<6.5) weight = 0.105829;
if(x>=6.5 && x<6.6) weight = 0.0983275;
if(x>=6.6 && x<6.7) weight = 0.0950082;
if(x>=6.7 && x<6.8) weight = 0.0853564;
if(x>=6.8 && x<6.9) weight = 0.0809471;
if(x>=6.9 && x<7) weight = 0.0756388;
if(x>=7 && x<7.1) weight = 0.0713881;
if(x>=7.1 && x<7.2) weight = 0.0661315;
if(x>=7.2 && x<7.3) weight = 0.0633135;
if(x>=7.3 && x<7.4) weight = 0.0575921;
if(x>=7.4 && x<7.5) weight = 0.0542962;
if(x>=7.5 && x<7.6) weight = 0.0505808;
if(x>=7.6 && x<7.7) weight = 0.0467498;
if(x>=7.7 && x<7.79999) weight = 0.0439359;
if(x>=7.79999 && x<7.89999) weight = 0.0415339;
if(x>=7.89999 && x<7.99999) weight = 0.0384786;
if(x>=7.99999 && x<8.09999) weight = 0.0362496;
if(x>=8.09999 && x<8.2) weight = 0.0348247;
if(x>=8.2 && x<8.3) weight = 0.0323179;
if(x>=8.3 && x<8.4) weight = 0.0299472;
if(x>=8.4 && x<8.5) weight = 0.0286392;
if(x>=8.5 && x<8.6) weight = 0.0265999;
if(x>=8.6 && x<8.7) weight = 0.0251215;
if(x>=8.7 && x<8.8) weight = 0.023278;
if(x>=8.8 && x<8.9) weight = 0.0230172;
if(x>=8.9 && x<9) weight = 0.0228469;
if(x>=9 && x<9.1) weight = 0.0220172;
if(x>=9.1 && x<9.2) weight = 0.020236;
if(x>=9.2 && x<9.3) weight = 0.0184216;
if(x>=9.3 && x<9.4) weight = 0.0166398;
if(x>=9.4 && x<9.5) weight = 0.0159598;
if(x>=9.5 && x<9.6) weight = 0.0154042;
if(x>=9.6 && x<9.7) weight = 0.0149367;
if(x>=9.7 && x<9.8) weight = 0.0141611;
if(x>=9.8 && x<9.9) weight = 0.0138181;
if(x>=9.9 && x<10) weight = 0.0122691;
if(x>=10 && x<10.1) weight = 0.012681;
if(x>=10.1 && x<10.2) weight = 0.0130359;
if(x>=10.2 && x<10.3) weight = 0.0108548;
if(x>=10.3 && x<10.4) weight = 0.0106733;
if(x>=10.4 && x<10.5) weight = 0.0103862;
if(x>=10.5 && x<10.6) weight = 0.00957966;
if(x>=10.6 && x<10.7) weight = 0.00906199;
if(x>=10.7 && x<10.8) weight = 0.0087744;
if(x>=10.8 && x<10.9) weight = 0.00769096;
if(x>=10.9 && x<11) weight = 0.00753246;
if(x>=11 && x<11.1) weight = 0.00790024;
if(x>=11.1 && x<11.2) weight = 0.00741829;
if(x>=11.2 && x<11.3) weight = 0.00694702;
if(x>=11.3 && x<11.4) weight = 0.00642725;
if(x>=11.4 && x<11.5) weight = 0.00603252;
if(x>=11.5 && x<11.6) weight = 0.00578299;
if(x>=11.6 && x<11.7) weight = 0.00556848;
if(x>=11.7 && x<11.8) weight = 0.00554197;
if(x>=11.8 && x<11.9) weight = 0.00500483;
if(x>=11.9 && x<12) weight = 0.00499248;
if(x>=12 && x<12.1) weight = 0.00455349;
if(x>=12.1 && x<12.2) weight = 0.00493805;
if(x>=12.2 && x<12.3) weight = 0.00434806;
if(x>=12.3 && x<12.4) weight = 0.00458255;
if(x>=12.4 && x<12.5) weight = 0.00434505;
if(x>=12.5 && x<12.6) weight = 0.00414693;
if(x>=12.6 && x<12.7) weight = 0.00405504;
if(x>=12.7 && x<12.8) weight = 0.00373716;
if(x>=12.8 && x<12.9) weight = 0.00373466;
if(x>=12.9 && x<13) weight = 0.00333056;
if(x>=13 && x<13.1) weight = 0.00332144;
if(x>=13.1 && x<13.2) weight = 0.00294159;
if(x>=13.2 && x<13.3) weight = 0.00247827;
if(x>=13.3 && x<13.4) weight = 0.00313991;
if(x>=13.4 && x<13.5) weight = 0.0026059;
if(x>=13.5 && x<13.6) weight = 0.00269088;
if(x>=13.6 && x<13.7) weight = 0.00237075;
if(x>=13.7 && x<13.8) weight = 0.00218079;
if(x>=13.8 && x<13.9) weight = 0.00246232;
if(x>=13.9 && x<14) weight = 0.00197782;
if(x>=14 && x<14.1) weight = 0.00240045;
if(x>=14.1 && x<14.2) weight = 0.00233604;
if(x>=14.2 && x<14.3) weight = 0.00215253;
if(x>=14.3 && x<14.4) weight = 0.00201289;
if(x>=14.4 && x<14.5) weight = 0.00167549;
if(x>=14.5 && x<14.6) weight = 0.00154853;
if(x>=14.6 && x<14.7) weight = 0.00160888;
if(x>=14.7 && x<14.8) weight = 0.00178278;
if(x>=14.8 && x<14.9) weight = 0.00168676;
if(x>=14.9 && x<15) weight = 0.00144662;
if(x>=15 && x<15.1) weight = 0.00158789;
if(x>=15.1 && x<15.2) weight = 0.0019857;
if(x>=15.2 && x<15.3) weight = 0.00150416;
if(x>=15.3 && x<15.4) weight = 0.00154834;
if(x>=15.4 && x<15.5) weight = 0.00144237;
if(x>=15.5 && x<15.6) weight = 0.0016616;
if(x>=15.6 && x<15.7) weight = 0.00125402;
if(x>=15.7 && x<15.8) weight = 0.00135707;
if(x>=15.8 && x<15.9) weight = 0.00125916;
if(x>=15.9 && x<16) weight = 0.00121498;
if(x>=16 && x<16.1) weight = 0.00109933;
if(x>=16.1 && x<16.2) weight = 0.000924099;
if(x>=16.2 && x<16.3) weight = 0.00120679;
if(x>=16.3 && x<16.4) weight = 0.000971354;
if(x>=16.4 && x<16.5) weight = 0.000941452;
if(x>=16.5 && x<16.6) weight = 0.000653479;
if(x>=16.6 && x<16.7) weight = 0.000886997;
if(x>=16.7 && x<16.8) weight = 0.000711136;
if(x>=16.8 && x<16.9) weight = 0.000914559;
if(x>=16.9 && x<17) weight = 0.00082989;
if(x>=17 && x<17.1) weight = 0.00091638;
if(x>=17.1 && x<17.2) weight = 0.000866092;
if(x>=17.2 && x<17.3) weight = 0.000736713;
if(x>=17.3 && x<17.4) weight = 0.000894596;
if(x>=17.4 && x<17.5) weight = 0.000618506;
if(x>=17.5 && x<17.6) weight = 0.000707564;
if(x>=17.6 && x<17.7) weight = 0.000510159;
if(x>=17.7 && x<17.8) weight = 0.000637552;
if(x>=17.8 && x<17.9) weight = 0.000692545;
if(x>=17.9 && x<18) weight = 0.000690938;
if(x>=18 && x<18.1) weight = 0.000511599;
if(x>=18.1 && x<18.2) weight = 0.000728229;
if(x>=18.2 && x<18.3) weight = 0.000707076;
if(x>=18.3 && x<18.4) weight = 0.000566121;
if(x>=18.4 && x<18.5) weight = 0.000459486;
if(x>=18.5 && x<18.6) weight = 0.000458101;
if(x>=18.6 && x<18.7) weight = 0.000582154;
if(x>=18.7 && x<18.8) weight = 0.000492732;
if(x>=18.8 && x<18.9) weight = 0.000424651;
if(x>=18.9 && x<19) weight = 0.000527398;
if(x>=19 && x<19.1) weight = 0.000459729;
if(x>=19.1 && x<19.2) weight = 0.000373486;
if(x>=19.2 && x<19.3) weight = 0.000533618;
if(x>=19.3 && x<19.4) weight = 0.00042597;
if(x>=19.4 && x<19.5) weight = 0.000424148;
if(x>=19.5 && x<19.6) weight = 0.0006327;
if(x>=19.6 && x<19.7) weight = 0.000390057;
if(x>=19.7 && x<19.8) weight = 0.000354032;
if(x>=19.8 && x<19.9) weight = 0.000404253;
if(x>=19.9 && x<20) weight = 0.00042523;
if(x>=20 && x<20.1) weight = 0.000353001;
if(x>=20.1 && x<20.2) weight = 0.000441119;
if(x>=20.2 && x<20.3) weight = 0.000211595;
if(x>=20.3 && x<20.4) weight = 0.000426386;
if(x>=20.4 && x<20.5) weight = 0.0002817;
if(x>=20.5 && x<20.6) weight = 0.00026401;
if(x>=20.6 && x<20.7) weight = 0.000351983;
if(x>=20.7 && x<20.8) weight = 0.000299686;
if(x>=20.8 && x<20.9) weight = 0.000141363;
if(x>=20.9 && x<21) weight = 0.000191674;
if(x>=21 && x<21.1) weight = 0.000317707;
if(x>=21.1 && x<21.2) weight = 0.000444808;
if(x>=21.2 && x<21.3) weight = 0.000246561;
if(x>=21.3 && x<21.4) weight = 0.000264816;
if(x>=21.4 && x<21.5) weight = 0.000266747;
if(x>=21.5 && x<21.6) weight = 0.000140553;
if(x>=21.6 && x<21.7) weight = 0.000265778;
if(x>=21.7 && x<21.8) weight = 0.000193441;
if(x>=21.8 && x<21.9) weight = 0.000247845;
if(x>=21.9 && x<22) weight = 0.000158543;
if(x>=22 && x<22.1) weight = 0.000229609;
if(x>=22.1 && x<22.2) weight = 0.000244837;
if(x>=22.2 && x<22.3) weight = 0.000160034;
if(x>=22.3 && x<22.4) weight = 0.000141121;
if(x>=22.4 && x<22.5) weight = 0.000229682;
if(x>=22.5 && x<22.6) weight = 0.00015943;
if(x>=22.6 && x<22.7001) weight = 0.000193839;
if(x>=22.7001 && x<22.8001) weight = 0.000195091;
if(x>=22.8001 && x<22.9001) weight = 7.07952e-05;
if(x>=22.9001 && x<23.0001) weight = 0.000123697;
if(x>=23.0001 && x<23.1001) weight = 0.000247972;
if(x>=23.1001 && x<23.2001) weight = 0.000212288;
if(x>=23.2001 && x<23.3001) weight = 7.10076e-05;
if(x>=23.3001 && x<23.4001) weight = 7.07551e-05;
if(x>=23.4001 && x<23.5001) weight = 0.000212502;
if(x>=23.5001 && x<23.6001) weight = 0.000106163;
if(x>=23.6001 && x<23.7001) weight = 0.00014039;
if(x>=23.7001 && x<23.8001) weight = 0.000159929;
if(x>=23.8001 && x<23.9001) weight = 0.000176028;
if(x>=23.9001 && x<24.0001) weight = 0.0001053;
if(x>=24.0001 && x<24.1001) weight = 0.000141525;
if(x>=24.1001 && x<24.2001) weight = 7.09534e-05;
if(x>=24.2001 && x<24.3001) weight = 0.000265825;
if(x>=24.3001 && x<24.4001) weight = 0.000141638;
if(x>=24.4001 && x<24.5001) weight = 0.000141288;
if(x>=24.5001 && x<24.6001) weight = 0.000141213;
if(x>=24.6001 && x<24.7001) weight = 0.000176697;
if(x>=24.7001 && x<24.8001) weight = 5.34207e-05;
if(x>=24.8001 && x<24.9001) weight = 0.000124054;
if(x>=24.9001 && x<25.0001) weight = 0.000141221;
if(x>=25.0001 && x<25.1001) weight = 7.07626e-05;
if(x>=25.1001 && x<25.2001) weight = 0.000159515;
if(x>=25.2001 && x<25.3001) weight = 8.81181e-05;
if(x>=25.3001 && x<25.4001) weight = 0.000106589;
if(x>=25.4001 && x<25.5001) weight = 0;
if(x>=25.5001 && x<25.6001) weight = 5.28728e-05;
if(x>=25.6001 && x<25.7001) weight = 0.000123311;
if(x>=25.7001 && x<25.8001) weight = 0.000106234;
if(x>=25.8001 && x<25.9001) weight = 0.00015905;
if(x>=25.9001 && x<26.0001) weight = 0.000159932;
if(x>=26.0001 && x<26.1001) weight = 7.07864e-05;
if(x>=26.1001 && x<26.2001) weight = 0.000176616;
if(x>=26.2001 && x<26.3001) weight = 8.82675e-05;
if(x>=26.3001 && x<26.4001) weight = 0.00014081;
if(x>=26.4001 && x<26.5001) weight = 8.82846e-05;
if(x>=26.5001 && x<26.6001) weight = 0.000106387;
if(x>=26.6001 && x<26.7001) weight = 8.82021e-05;
if(x>=26.7001 && x<26.8001) weight = 0.000106629;
if(x>=26.8001 && x<26.9001) weight = 8.84643e-05;
if(x>=26.9001 && x<27.0001) weight = 5.28095e-05;
if(x>=27.0001 && x<27.1001) weight = 5.30448e-05;
if(x>=27.1001 && x<27.2001) weight = 6.98544e-05;
if(x>=27.2001 && x<27.3001) weight = 0.000105917;
if(x>=27.3001 && x<27.4001) weight = 8.83439e-05;
if(x>=27.4001 && x<27.5001) weight = 5.29362e-05;
if(x>=27.5001 && x<27.6001) weight = 8.87012e-05;
if(x>=27.6001 && x<27.7001) weight = 5.29932e-05;
if(x>=27.7001 && x<27.8001) weight = 7.02593e-05;
if(x>=27.8001 && x<27.9001) weight = 0.00012357;
if(x>=27.9001 && x<28.0001) weight = 1.77349e-05;
if(x>=28.0001 && x<28.1001) weight = 0.000106667;
if(x>=28.1001 && x<28.2001) weight = 5.30195e-05;
if(x>=28.2001 && x<28.3001) weight = 7.08391e-05;
if(x>=28.3001 && x<28.4001) weight = 8.7747e-05;
if(x>=28.4001 && x<28.5001) weight = 0.000105947;
if(x>=28.5001 && x<28.6001) weight = 1.75165e-05;
if(x>=28.6001 && x<28.7001) weight = 5.32926e-05;
if(x>=28.7001 && x<28.8001) weight = 1.76697e-05;
if(x>=28.8001 && x<28.9001) weight = 5.29502e-05;
if(x>=28.9001 && x<29.0001) weight = 5.26556e-05;
if(x>=29.0001 && x<29.1001) weight = 1.75639e-05;
if(x>=29.1001 && x<29.2001) weight = 0;
if(x>=29.2001 && x<29.3001) weight = 5.2803e-05;
if(x>=29.3001 && x<29.4001) weight = 0;
if(x>=29.4001 && x<29.5001) weight = 3.51661e-05;
if(x>=29.5001 && x<29.6001) weight = 1.76991e-05;
if(x>=29.6001 && x<29.7001) weight = 7.07589e-05;
if(x>=29.7001 && x<29.8001) weight = 1.78079e-05;
if(x>=29.8001 && x<29.9001) weight = 0.000106542;
if(x>=29.9001 && x<30.0001) weight = 5.3072e-05;
if(x>=30.0001 && x<30.1001) weight = 3.54597e-05;


		}
		else if(pdg_particle==221){ //eta
		

///Weight for only enhanced spectra (choosing IsPrimary in both MC)

if(x>=0 && x<0.1) weight = 16.6482;
if(x>=0.1 && x<0.2) weight = 35.3778;
if(x>=0.2 && x<0.3) weight = 49.3503;
if(x>=0.3 && x<0.4) weight = 55.1089;
if(x>=0.4 && x<0.5) weight = 54.5517;
if(x>=0.5 && x<0.6) weight = 50.0072;
if(x>=0.6 && x<0.7) weight = 43.8912;
if(x>=0.7 && x<0.8) weight = 36.9331;
if(x>=0.8 && x<0.9) weight = 30.8918;
if(x>=0.9 && x<1) weight = 25.4834;
if(x>=1 && x<1.1) weight = 21.0065;
if(x>=1.1 && x<1.2) weight = 17.4224;
if(x>=1.2 && x<1.3) weight = 14.6096;
if(x>=1.3 && x<1.4) weight = 12.1366;
if(x>=1.4 && x<1.5) weight = 10.1297;
if(x>=1.5 && x<1.6) weight = 8.56067;
if(x>=1.6 && x<1.7) weight = 7.18474;
if(x>=1.7 && x<1.8) weight = 6.06375;
if(x>=1.8 && x<1.9) weight = 5.15489;
if(x>=1.9 && x<2) weight = 4.3911;
if(x>=2 && x<2.1) weight = 3.78573;
if(x>=2.1 && x<2.2) weight = 3.20524;
if(x>=2.2 && x<2.3) weight = 2.76667;
if(x>=2.3 && x<2.4) weight = 2.358;
if(x>=2.4 && x<2.5) weight = 2.05094;
if(x>=2.5 && x<2.6) weight = 1.79703;
if(x>=2.6 && x<2.7) weight = 1.54752;
if(x>=2.7 && x<2.8) weight = 1.35948;
if(x>=2.8 && x<2.9) weight = 1.19015;
if(x>=2.9 && x<3) weight = 1.03754;
if(x>=3 && x<3.1) weight = 0.911192;
if(x>=3.1 && x<3.2) weight = 0.799433;
if(x>=3.2 && x<3.3) weight = 0.719964;
if(x>=3.3 && x<3.4) weight = 0.642558;
if(x>=3.4 && x<3.5) weight = 0.567598;
if(x>=3.5 && x<3.6) weight = 0.502407;
if(x>=3.6 && x<3.7) weight = 0.449138;
if(x>=3.7 && x<3.8) weight = 0.409392;
if(x>=3.8 && x<3.9) weight = 0.360166;
if(x>=3.9 && x<4) weight = 0.328277;
if(x>=4 && x<4.1) weight = 0.29631;
if(x>=4.1 && x<4.2) weight = 0.271259;
if(x>=4.2 && x<4.3) weight = 0.239338;
if(x>=4.3 && x<4.4) weight = 0.215612;
if(x>=4.4 && x<4.5) weight = 0.197521;
if(x>=4.5 && x<4.6) weight = 0.17469;
if(x>=4.6 && x<4.7) weight = 0.162464;
if(x>=4.7 && x<4.8) weight = 0.15101;
if(x>=4.8 && x<4.9) weight = 0.134913;
if(x>=4.9 && x<5) weight = 0.124206;
if(x>=5 && x<5.1) weight = 0.114006;
if(x>=5.1 && x<5.2) weight = 0.104045;
if(x>=5.2 && x<5.3) weight = 0.0947138;
if(x>=5.3 && x<5.4) weight = 0.0872681;
if(x>=5.4 && x<5.5) weight = 0.0826883;
if(x>=5.5 && x<5.6) weight = 0.0763556;
if(x>=5.6 && x<5.7) weight = 0.0698036;
if(x>=5.7 && x<5.8) weight = 0.0602925;
if(x>=5.8 && x<5.9) weight = 0.0602648;
if(x>=5.9 && x<6) weight = 0.0554156;
if(x>=6 && x<6.1) weight = 0.0512432;
if(x>=6.1 && x<6.2) weight = 0.0465133;
if(x>=6.2 && x<6.3) weight = 0.0423586;
if(x>=6.3 && x<6.4) weight = 0.0406576;
if(x>=6.4 && x<6.5) weight = 0.0372029;
if(x>=6.5 && x<6.6) weight = 0.0344116;
if(x>=6.6 && x<6.7) weight = 0.033169;
if(x>=6.7 && x<6.8) weight = 0.0302002;
if(x>=6.8 && x<6.9) weight = 0.0282926;
if(x>=6.9 && x<7) weight = 0.0267922;
if(x>=7 && x<7.1) weight = 0.0251031;
if(x>=7.1 && x<7.2) weight = 0.0228793;
if(x>=7.2 && x<7.3) weight = 0.02128;
if(x>=7.3 && x<7.4) weight = 0.0193766;
if(x>=7.4 && x<7.5) weight = 0.0193382;
if(x>=7.5 && x<7.6) weight = 0.0176085;
if(x>=7.6 && x<7.7) weight = 0.017612;
if(x>=7.7 && x<7.79999) weight = 0.0156457;
if(x>=7.79999 && x<7.89999) weight = 0.0154013;
if(x>=7.89999 && x<7.99999) weight = 0.0136434;
if(x>=7.99999 && x<8.09999) weight = 0.0134926;
if(x>=8.09999 && x<8.2) weight = 0.0121375;
if(x>=8.2 && x<8.3) weight = 0.0119262;
if(x>=8.3 && x<8.4) weight = 0.0106799;
if(x>=8.4 && x<8.5) weight = 0.010546;
if(x>=8.5 && x<8.6) weight = 0.0095009;
if(x>=8.6 && x<8.7) weight = 0.00924352;
if(x>=8.7 && x<8.8) weight = 0.00927782;
if(x>=8.8 && x<8.9) weight = 0.00867205;
if(x>=8.9 && x<9) weight = 0.00769503;
if(x>=9 && x<9.1) weight = 0.00690412;
if(x>=9.1 && x<9.2) weight = 0.00707209;
if(x>=9.2 && x<9.3) weight = 0.00621283;
if(x>=9.3 && x<9.4) weight = 0.00680911;
if(x>=9.4 && x<9.5) weight = 0.00654103;
if(x>=9.5 && x<9.6) weight = 0.00578403;
if(x>=9.6 && x<9.7) weight = 0.00569871;
if(x>=9.7 && x<9.8) weight = 0.00515655;
if(x>=9.8 && x<9.9) weight = 0.00482338;
if(x>=9.9 && x<10) weight = 0.00450786;
if(x>=10 && x<10.1) weight = 0.00471897;
if(x>=10.1 && x<10.2) weight = 0.00451627;
if(x>=10.2 && x<10.3) weight = 0.00392212;
if(x>=10.3 && x<10.4) weight = 0.00387302;
if(x>=10.4 && x<10.5) weight = 0.0033188;
if(x>=10.5 && x<10.6) weight = 0.00322433;
if(x>=10.6 && x<10.7) weight = 0.00342829;
if(x>=10.7 && x<10.8) weight = 0.00291041;
if(x>=10.8 && x<10.9) weight = 0.0031883;
if(x>=10.9 && x<11) weight = 0.00296307;
if(x>=11 && x<11.1) weight = 0.00242486;
if(x>=11.1 && x<11.2) weight = 0.00215176;
if(x>=11.2 && x<11.3) weight = 0.00271103;
if(x>=11.3 && x<11.4) weight = 0.00248752;
if(x>=11.4 && x<11.5) weight = 0.00260421;
if(x>=11.5 && x<11.6) weight = 0.0022926;
if(x>=11.6 && x<11.7) weight = 0.00196885;
if(x>=11.7 && x<11.8) weight = 0.00192537;
if(x>=11.8 && x<11.9) weight = 0.00188217;
if(x>=11.9 && x<12) weight = 0.00159044;
if(x>=12 && x<12.1) weight = 0.00179997;
if(x>=12.1 && x<12.2) weight = 0.00184061;
if(x>=12.2 && x<12.3) weight = 0.00172453;
if(x>=12.3 && x<12.4) weight = 0.00143019;
if(x>=12.4 && x<12.5) weight = 0.00114075;
if(x>=12.5 && x<12.6) weight = 0.00123192;
if(x>=12.6 && x<12.7) weight = 0.00132884;
if(x>=12.7 && x<12.8) weight = 0.00123684;
if(x>=12.8 && x<12.9) weight = 0.00111245;
if(x>=12.9 && x<13) weight = 0.00106313;
if(x>=13 && x<13.1) weight = 0.00104033;
if(x>=13.1 && x<13.2) weight = 0.000883221;
if(x>=13.2 && x<13.3) weight = 0.00107319;
if(x>=13.3 && x<13.4) weight = 0.000918923;
if(x>=13.4 && x<13.5) weight = 0.000876117;
if(x>=13.5 && x<13.6) weight = 0.00105402;
if(x>=13.6 && x<13.7) weight = 0.000850461;
if(x>=13.7 && x<13.8) weight = 0.000817066;
if(x>=13.8 && x<13.9) weight = 0.00088881;
if(x>=13.9 && x<14) weight = 0.000760981;
if(x>=14 && x<14.1) weight = 0.000515987;
if(x>=14.1 && x<14.2) weight = 0.000707889;
if(x>=14.2 && x<14.3) weight = 0.000793427;
if(x>=14.3 && x<14.4) weight = 0.000632889;
if(x>=14.4 && x<14.5) weight = 0.000583297;
if(x>=14.5 && x<14.6) weight = 0.000547094;
if(x>=14.6 && x<14.7) weight = 0.000727144;
if(x>=14.7 && x<14.8) weight = 0.000514403;
if(x>=14.8 && x<14.9) weight = 0.000705978;
if(x>=14.9 && x<15) weight = 0.00059817;
if(x>=15 && x<15.1) weight = 0.000476073;
if(x>=15.1 && x<15.2) weight = 0.000632644;
if(x>=15.2 && x<15.3) weight = 0.000527389;
if(x>=15.3 && x<15.4) weight = 0.000478172;
if(x>=15.4 && x<15.5) weight = 0.000369731;
if(x>=15.5 && x<15.6) weight = 0.000459664;
if(x>=15.6 && x<15.7) weight = 0.000580403;
if(x>=15.7 && x<15.8) weight = 0.000527361;
if(x>=15.8 && x<15.9) weight = 0.000548906;
if(x>=15.9 && x<16) weight = 0.0004772;
if(x>=16 && x<16.1) weight = 0.000301312;
if(x>=16.1 && x<16.2) weight = 0.000446476;
if(x>=16.2 && x<16.3) weight = 0.000493531;
if(x>=16.3 && x<16.4) weight = 0.000353382;
if(x>=16.4 && x<16.5) weight = 0.000246574;
if(x>=16.5 && x<16.6) weight = 0.000247394;
if(x>=16.6 && x<16.7) weight = 0.000281279;
if(x>=16.7 && x<16.8) weight = 0.00031778;
if(x>=16.8 && x<16.9) weight = 0.000371754;
if(x>=16.9 && x<17) weight = 0.000387072;
if(x>=17 && x<17.1) weight = 0.000282965;
if(x>=17.1 && x<17.2) weight = 0.000283166;
if(x>=17.2 && x<17.3) weight = 0.000141405;
if(x>=17.3 && x<17.4) weight = 0.000351871;
if(x>=17.4 && x<17.5) weight = 0.000245761;
if(x>=17.5 && x<17.6) weight = 0.00022982;
if(x>=17.6 && x<17.7) weight = 0.000106538;
if(x>=17.7 && x<17.8) weight = 0.000229682;
if(x>=17.8 && x<17.9) weight = 0.000176491;
if(x>=17.9 && x<18) weight = 0.000141039;
if(x>=18 && x<18.1) weight = 0.000176454;
if(x>=18.1 && x<18.2) weight = 0.000177478;
if(x>=18.2 && x<18.3) weight = 0.000176177;
if(x>=18.3 && x<18.4) weight = 0.000209644;
if(x>=18.4 && x<18.5) weight = 0.000176922;
if(x>=18.5 && x<18.6) weight = 0.000282795;
if(x>=18.6 && x<18.7) weight = 0.000194449;
if(x>=18.7 && x<18.8) weight = 0.000122181;
if(x>=18.8 && x<18.9) weight = 0.000265623;
if(x>=18.9 && x<19) weight = 0.000141138;
if(x>=19 && x<19.1) weight = 0.000212228;
if(x>=19.1 && x<19.2) weight = 7.08793e-05;
if(x>=19.2 && x<19.3) weight = 0.000210534;
if(x>=19.3 && x<19.4) weight = 8.77285e-05;
if(x>=19.4 && x<19.5) weight = 0.000158744;
if(x>=19.5 && x<19.6) weight = 0.00015898;
if(x>=19.6 && x<19.7) weight = 0.000123359;
if(x>=19.7 && x<19.8) weight = 0.000175346;
if(x>=19.8 && x<19.9) weight = 0.000106185;
if(x>=19.9 && x<20) weight = 0.000159351;
if(x>=20 && x<20.1) weight = 0.000106204;
if(x>=20.1 && x<20.2) weight = 5.28793e-05;
if(x>=20.2 && x<20.3) weight = 0.000105435;
if(x>=20.3 && x<20.4) weight = 0.000229674;
if(x>=20.4 && x<20.5) weight = 7.08215e-05;
if(x>=20.5 && x<20.6) weight = 8.7991e-05;
if(x>=20.6 && x<20.7) weight = 0.000105326;
if(x>=20.7 && x<20.8) weight = 0.000106236;
if(x>=20.8 && x<20.9) weight = 0.000105865;
if(x>=20.9 && x<21) weight = 8.83689e-05;
if(x>=21 && x<21.1) weight = 8.85881e-05;
if(x>=21.1 && x<21.2) weight = 5.31378e-05;
if(x>=21.2 && x<21.3) weight = 8.87406e-05;
if(x>=21.3 && x<21.4) weight = 7.04362e-05;
if(x>=21.4 && x<21.5) weight = 0.000159269;
if(x>=21.5 && x<21.6) weight = 7.10846e-05;
if(x>=21.6 && x<21.7) weight = 5.30701e-05;
if(x>=21.7 && x<21.8) weight = 5.30579e-05;
if(x>=21.8 && x<21.9) weight = 7.04411e-05;
if(x>=21.9 && x<22) weight = 7.0878e-05;
if(x>=22 && x<22.1) weight = 5.28243e-05;
if(x>=22.1 && x<22.2) weight = 8.87894e-05;
if(x>=22.2 && x<22.3) weight = 0.000106099;
if(x>=22.3 && x<22.4) weight = 1.76059e-05;
if(x>=22.4 && x<22.5) weight = 5.35313e-05;
if(x>=22.5 && x<22.6) weight = 3.52479e-05;
if(x>=22.6 && x<22.7001) weight = 3.55372e-05;
if(x>=22.7001 && x<22.8001) weight = 1.755e-05;
if(x>=22.8001 && x<22.9001) weight = 7.09409e-05;
if(x>=22.9001 && x<23.0001) weight = 5.30504e-05;
if(x>=23.0001 && x<23.1001) weight = 0.000106005;
if(x>=23.1001 && x<23.2001) weight = 3.54365e-05;
if(x>=23.2001 && x<23.3001) weight = 0.000122684;
if(x>=23.3001 && x<23.4001) weight = 5.30485e-05;
if(x>=23.4001 && x<23.5001) weight = 0.000141864;
if(x>=23.5001 && x<23.6001) weight = 0;
if(x>=23.6001 && x<23.7001) weight = 7.00341e-05;
if(x>=23.7001 && x<23.8001) weight = 3.52839e-05;
if(x>=23.8001 && x<23.9001) weight = 0.000124098;
if(x>=23.9001 && x<24.0001) weight = 0.000106519;
if(x>=24.0001 && x<24.1001) weight = 0.000123925;
if(x>=24.1001 && x<24.2001) weight = 3.52709e-05;
if(x>=24.2001 && x<24.3001) weight = 3.53788e-05;
if(x>=24.3001 && x<24.4001) weight = 8.85724e-05;
if(x>=24.4001 && x<24.5001) weight = 1.76732e-05;
if(x>=24.5001 && x<24.6001) weight = 3.54321e-05;
if(x>=24.6001 && x<24.7001) weight = 0;
if(x>=24.7001 && x<24.8001) weight = 1.75389e-05;
if(x>=24.8001 && x<24.9001) weight = 3.55082e-05;
if(x>=24.9001 && x<25.0001) weight = 5.29241e-05;
if(x>=25.0001 && x<25.1001) weight = 5.28123e-05;
if(x>=25.1001 && x<25.2001) weight = 1.76735e-05;
if(x>=25.2001 && x<25.3001) weight = 7.11718e-05;
if(x>=25.3001 && x<25.4001) weight = 1.76725e-05;
if(x>=25.4001 && x<25.5001) weight = 5.3434e-05;
if(x>=25.5001 && x<25.6001) weight = 0;
if(x>=25.6001 && x<25.7001) weight = 7.0263e-05;
if(x>=25.7001 && x<25.8001) weight = 3.54943e-05;
if(x>=25.8001 && x<25.9001) weight = 0;
if(x>=25.9001 && x<26.0001) weight = 3.52802e-05;
if(x>=26.0001 && x<26.1001) weight = 1.76757e-05;
if(x>=26.1001 && x<26.2001) weight = 1.77368e-05;
if(x>=26.2001 && x<26.3001) weight = 3.52361e-05;
if(x>=26.3001 && x<26.4001) weight = 1.75193e-05;
if(x>=26.4001 && x<26.5001) weight = 1.76813e-05;
if(x>=26.5001 && x<26.6001) weight = 1.7546e-05;
if(x>=26.6001 && x<26.7001) weight = 3.53957e-05;
if(x>=26.7001 && x<26.8001) weight = 1.76171e-05;
if(x>=26.8001 && x<26.9001) weight = 3.54868e-05;
if(x>=26.9001 && x<27.0001) weight = 1.76146e-05;
if(x>=27.0001 && x<27.1001) weight = 0;
if(x>=27.1001 && x<27.2001) weight = 0;
if(x>=27.2001 && x<27.3001) weight = 5.32756e-05;
if(x>=27.3001 && x<27.4001) weight = 0;
if(x>=27.4001 && x<27.5001) weight = 1.77365e-05;
if(x>=27.5001 && x<27.6001) weight = 3.52914e-05;
if(x>=27.6001 && x<27.7001) weight = 0;
if(x>=27.7001 && x<27.8001) weight = 1.77519e-05;
if(x>=27.8001 && x<27.9001) weight = 0;
if(x>=27.9001 && x<28.0001) weight = 1.75405e-05;
if(x>=28.0001 && x<28.1001) weight = 3.53682e-05;
if(x>=28.1001 && x<28.2001) weight = 0;
if(x>=28.2001 && x<28.3001) weight = 3.54736e-05;
if(x>=28.3001 && x<28.4001) weight = 5.2967e-05;
if(x>=28.4001 && x<28.5001) weight = 3.54679e-05;
if(x>=28.5001 && x<28.6001) weight = 0;
if(x>=28.6001 && x<28.7001) weight = 5.31953e-05;
if(x>=28.7001 && x<28.8001) weight = 1.76056e-05;
if(x>=28.8001 && x<28.9001) weight = 1.77418e-05;
if(x>=28.9001 && x<29.0001) weight = 1.7744e-05;
if(x>=29.0001 && x<29.1001) weight = 1.76373e-05;
if(x>=29.1001 && x<29.2001) weight = 5.3167e-05;
if(x>=29.2001 && x<29.3001) weight = 0;
if(x>=29.3001 && x<29.4001) weight = 0;
if(x>=29.4001 && x<29.5001) weight = 3.52877e-05;
if(x>=29.5001 && x<29.6001) weight = 1.76641e-05;
if(x>=29.6001 && x<29.7001) weight = 0;
if(x>=29.7001 && x<29.8001) weight = 0;
if(x>=29.8001 && x<29.9001) weight = 5.2843e-05;
if(x>=29.9001 && x<30.0001) weight = 1.76479e-05;
if(x>=30.0001 && x<30.1001) weight = 0;

		}

	return weight;
	}
	
	///Weight for Pb-Pb
	if(!fIsPP){
		
		if(pdg_particle==111){ //pi0
			
if(x>=0 && x<0.1) weight = 3.32419;
if(x>=0.1 && x<0.2) weight = 2.54442;
if(x>=0.2 && x<0.3) weight = 1.82559;
if(x>=0.3 && x<0.4) weight = 1.3684;
if(x>=0.4 && x<0.5) weight = 1.08552;
if(x>=0.5 && x<0.6) weight = 0.915234;
if(x>=0.6 && x<0.7) weight = 0.81781;
if(x>=0.7 && x<0.8) weight = 0.753406;
if(x>=0.8 && x<0.9) weight = 0.72134;
if(x>=0.9 && x<1) weight = 0.700434;
if(x>=1 && x<1.1) weight = 0.690641;
if(x>=1.1 && x<1.2) weight = 0.681571;
if(x>=1.2 && x<1.3) weight = 0.672113;
if(x>=1.3 && x<1.4) weight = 0.661982;
if(x>=1.4 && x<1.5) weight = 0.654034;
if(x>=1.5 && x<1.6) weight = 0.641112;
if(x>=1.6 && x<1.7) weight = 0.630549;
if(x>=1.7 && x<1.8) weight = 0.613403;
if(x>=1.8 && x<1.9) weight = 0.595767;
if(x>=1.9 && x<2) weight = 0.586549;
if(x>=2 && x<2.1) weight = 0.566889;
if(x>=2.1 && x<2.2) weight = 0.559352;
if(x>=2.2 && x<2.3) weight = 0.555745;
if(x>=2.3 && x<2.4) weight = 0.528827;
if(x>=2.4 && x<2.5) weight = 0.532158;
if(x>=2.5 && x<2.6) weight = 0.499868;
if(x>=2.6 && x<2.7) weight = 0.505048;
if(x>=2.7 && x<2.8) weight = 0.472686;
if(x>=2.8 && x<2.9) weight = 0.457234;
if(x>=2.9 && x<3) weight = 0.442031;
if(x>=3 && x<3.1) weight = 0.433296;
if(x>=3.1 && x<3.2) weight = 0.408119;
if(x>=3.2 && x<3.3) weight = 0.397404;
if(x>=3.3 && x<3.4) weight = 0.380257;
if(x>=3.4 && x<3.5) weight = 0.374139;
if(x>=3.5 && x<3.6) weight = 0.343083;
if(x>=3.6 && x<3.7) weight = 0.329866;
if(x>=3.7 && x<3.8) weight = 0.341286;
if(x>=3.8 && x<3.9) weight = 0.283544;
if(x>=3.9 && x<4) weight = 0.290994;
if(x>=4 && x<4.1) weight = 0.280361;
if(x>=4.1 && x<4.2) weight = 0.265524;
if(x>=4.2 && x<4.3) weight = 0.261704;
if(x>=4.3 && x<4.4) weight = 0.227433;
if(x>=4.4 && x<4.5) weight = 0.23766;
if(x>=4.5 && x<4.6) weight = 0.219045;
if(x>=4.6 && x<4.7) weight = 0.211493;
if(x>=4.7 && x<4.8) weight = 0.188017;
if(x>=4.8 && x<4.9) weight = 0.191124;
if(x>=4.9 && x<5) weight = 0.173571;
if(x>=5 && x<5.1) weight = 0.184366;
if(x>=5.1 && x<5.2) weight = 0.146323;
if(x>=5.2 && x<5.3) weight = 0.165052;
if(x>=5.3 && x<5.4) weight = 0.145557;
if(x>=5.4 && x<5.5) weight = 0.129918;
if(x>=5.5 && x<5.6) weight = 0.120794;
if(x>=5.6 && x<5.7) weight = 0.123121;
if(x>=5.7 && x<5.8) weight = 0.114557;
if(x>=5.8 && x<5.9) weight = 0.0999146;
if(x>=5.9 && x<6) weight = 0.0884551;
if(x>=6 && x<6.1) weight = 0.0768881;
if(x>=6.1 && x<6.2) weight = 0.091316;
if(x>=6.2 && x<6.3) weight = 0.0807397;
if(x>=6.3 && x<6.4) weight = 0.0838189;
if(x>=6.4 && x<6.5) weight = 0.0922858;
if(x>=6.5 && x<6.6) weight = 0.0758401;
if(x>=6.6 && x<6.7) weight = 0.0607455;
if(x>=6.7 && x<6.8) weight = 0.0555556;
if(x>=6.8 && x<6.9) weight = 0.0562771;
if(x>=6.9 && x<7) weight = 0.0599613;
if(x>=7 && x<7.1) weight = 0.0525813;
if(x>=7.1 && x<7.2) weight = 0.044316;
if(x>=7.2 && x<7.3) weight = 0.0477612;
if(x>=7.3 && x<7.4) weight = 0.0403789;
if(x>=7.4 && x<7.5) weight = 0.0400792;
if(x>=7.5 && x<7.6) weight = 0.0425321;
if(x>=7.6 && x<7.7) weight = 0.0401587;
if(x>=7.7 && x<7.79999) weight = 0.0321429;
if(x>=7.79999 && x<7.89999) weight = 0.0386399;
if(x>=7.89999 && x<7.99999) weight = 0.0284858;
if(x>=7.99999 && x<8.09999) weight = 0.0332311;
if(x>=8.09999 && x<8.2) weight = 0.022761;
if(x>=8.2 && x<8.3) weight = 0.0233051;
if(x>=8.3 && x<8.4) weight = 0.020979;
if(x>=8.4 && x<8.5) weight = 0.0222683;
if(x>=8.5 && x<8.6) weight = 0.024666;
if(x>=8.6 && x<8.7) weight = 0.0135065;
if(x>=8.7 && x<8.8) weight = 0.0167539;
if(x>=8.8 && x<8.9) weight = 0.0166493;
if(x>=8.9 && x<9) weight = 0.0129611;
if(x>=9 && x<9.1) weight = 0.018797;
if(x>=9.1 && x<9.2) weight = 0.0155233;
if(x>=9.2 && x<9.3) weight = 0.0119792;
if(x>=9.3 && x<9.4) weight = 0.0133476;
if(x>=9.4 && x<9.5) weight = 0.0106952;
if(x>=9.5 && x<9.6) weight = 0.0106045;
if(x>=9.6 && x<9.7) weight = 0.0135883;
if(x>=9.7 && x<9.8) weight = 0.0103734;
if(x>=9.8 && x<9.9) weight = 0.00968784;
if(x>=9.9 && x<10) weight = 0.0119177;
if(x>=10 && x<10.1) weight = 0.0115789;
if(x>=10.1 && x<10.2) weight = 0.00997899;
if(x>=10.2 && x<10.3) weight = 0.0060208;
if(x>=10.3 && x<10.4) weight = 0.00959488;
if(x>=10.4 && x<10.5) weight = 0.00730689;
if(x>=10.5 && x<10.6) weight = 0.00498615;
if(x>=10.6 && x<10.7) weight = 0.00757986;
if(x>=10.7 && x<10.8) weight = 0.0102987;
if(x>=10.8 && x<10.9) weight = 0.00477707;
if(x>=10.9 && x<11) weight = 0.00841662;
if(x>=11 && x<11.1) weight = 0.00793651;
if(x>=11.1 && x<11.2) weight = 0.00624675;
if(x>=11.2 && x<11.3) weight = 0.00513347;
if(x>=11.3 && x<11.4) weight = 0.00431267;
if(x>=11.4 && x<11.5) weight = 0.00595883;
if(x>=11.5 && x<11.6) weight = 0.00485699;
if(x>=11.6 && x<11.7) weight = 0.00529661;
if(x>=11.7 && x<11.8) weight = 0.00593952;
if(x>=11.8 && x<11.9) weight = 0.00546747;
if(x>=11.9 && x<12) weight = 0.00562114;
if(x>=12 && x<12.1) weight = 0.00304414;
if(x>=12.1 && x<12.2) weight = 0.00370566;
if(x>=12.2 && x<12.3) weight = 0.00271739;
if(x>=12.3 && x<12.4) weight = 0.00474934;
if(x>=12.4 && x<12.5) weight = 0.00478723;
if(x>=12.5 && x<12.6) weight = 0.0032967;
if(x>=12.6 && x<12.7) weight = 0.00377358;
if(x>=12.7 && x<12.8) weight = 0.00215633;
if(x>=12.8 && x<12.9) weight = 0.00309119;
if(x>=12.9 && x<13) weight = 0.000529661;
if(x>=13 && x<13.1) weight = 0;
if(x>=13.1 && x<13.2) weight = 0.00164384;
if(x>=13.2 && x<13.3) weight = 0;
if(x>=13.3 && x<13.4) weight = 0.00425306;
if(x>=13.4 && x<13.5) weight = 0.00101112;
if(x>=13.5 && x<13.6) weight = 0.00210416;
if(x>=13.6 && x<13.7) weight = 0.00112549;
if(x>=13.7 && x<13.8) weight = 0.00165654;
if(x>=13.8 && x<13.9) weight = 0.000547945;
if(x>=13.9 && x<14) weight = 0.000541418;
if(x>=14 && x<14.1) weight = 0.00162866;
if(x>=14.1 && x<14.2) weight = 0.00108578;
if(x>=14.2 && x<14.3) weight = 0.00172117;
if(x>=14.3 && x<14.4) weight = 0;
if(x>=14.4 && x<14.5) weight = 0.00163488;
if(x>=14.5 && x<14.6) weight = 0.00107759;
if(x>=14.6 && x<14.7) weight = 0.00214362;
if(x>=14.7 && x<14.8) weight = 0.000517063;
if(x>=14.8 && x<14.9) weight = 0.00209754;
if(x>=14.9 && x<15) weight = 0.000523013;
if(x>=15) weight = 0;
		
		}
		
		else if(pdg_particle==221){ //eta
			
if(x>=0 && x<0.1) weight = 0.498069;
if(x>=0.1 && x<0.2) weight = 0.474634;
if(x>=0.2 && x<0.3) weight = 0.435564;
if(x>=0.3 && x<0.4) weight = 0.399732;
if(x>=0.4 && x<0.5) weight = 0.372838;
if(x>=0.5 && x<0.6) weight = 0.345814;
if(x>=0.6 && x<0.7) weight = 0.328843;
if(x>=0.7 && x<0.8) weight = 0.319475;
if(x>=0.8 && x<0.9) weight = 0.315264;
if(x>=0.9 && x<1) weight = 0.317843;
if(x>=1 && x<1.1) weight = 0.312453;
if(x>=1.1 && x<1.2) weight = 0.311473;
if(x>=1.2 && x<1.3) weight = 0.313118;
if(x>=1.3 && x<1.4) weight = 0.310571;
if(x>=1.4 && x<1.5) weight = 0.303369;
if(x>=1.5 && x<1.6) weight = 0.305545;
if(x>=1.6 && x<1.7) weight = 0.29942;
if(x>=1.7 && x<1.8) weight = 0.292572;
if(x>=1.8 && x<1.9) weight = 0.284717;
if(x>=1.9 && x<2) weight = 0.281138;
if(x>=2 && x<2.1) weight = 0.280071;
if(x>=2.1 && x<2.2) weight = 0.278327;
if(x>=2.2 && x<2.3) weight = 0.258907;
if(x>=2.3 && x<2.4) weight = 0.264442;
if(x>=2.4 && x<2.5) weight = 0.250617;
if(x>=2.5 && x<2.6) weight = 0.24021;
if(x>=2.6 && x<2.7) weight = 0.231397;
if(x>=2.7 && x<2.8) weight = 0.232913;
if(x>=2.8 && x<2.9) weight = 0.221626;
if(x>=2.9 && x<3) weight = 0.22948;
if(x>=3 && x<3.1) weight = 0.201028;
if(x>=3.1 && x<3.2) weight = 0.192935;
if(x>=3.2 && x<3.3) weight = 0.192692;
if(x>=3.3 && x<3.4) weight = 0.201646;
if(x>=3.4 && x<3.5) weight = 0.182323;
if(x>=3.5 && x<3.6) weight = 0.180809;
if(x>=3.6 && x<3.7) weight = 0.163413;
if(x>=3.7 && x<3.8) weight = 0.148758;
if(x>=3.8 && x<3.9) weight = 0.135414;
if(x>=3.9 && x<4) weight = 0.142185;
if(x>=4 && x<4.1) weight = 0.128976;
if(x>=4.1 && x<4.2) weight = 0.134124;
if(x>=4.2 && x<4.3) weight = 0.118199;
if(x>=4.3 && x<4.4) weight = 0.112536;
if(x>=4.4 && x<4.5) weight = 0.108796;
if(x>=4.5 && x<4.6) weight = 0.101893;
if(x>=4.6 && x<4.7) weight = 0.104934;
if(x>=4.7 && x<4.8) weight = 0.0956955;
if(x>=4.8 && x<4.9) weight = 0.0751717;
if(x>=4.9 && x<5) weight = 0.084803;
if(x>=5 && x<5.1) weight = 0.0730558;
if(x>=5.1 && x<5.2) weight = 0.0685805;
if(x>=5.2 && x<5.3) weight = 0.0729927;
if(x>=5.3 && x<5.4) weight = 0.0604805;
if(x>=5.4 && x<5.5) weight = 0.0645564;
if(x>=5.5 && x<5.6) weight = 0.0614574;
if(x>=5.6 && x<5.7) weight = 0.0449343;
if(x>=5.7 && x<5.8) weight = 0.0472944;
if(x>=5.8 && x<5.9) weight = 0.0436474;
if(x>=5.9 && x<6) weight = 0.0449537;
if(x>=6 && x<6.1) weight = 0.0556561;
if(x>=6.1 && x<6.2) weight = 0.0457457;
if(x>=6.2 && x<6.3) weight = 0.0402778;
if(x>=6.3 && x<6.4) weight = 0.0355598;
if(x>=6.4 && x<6.5) weight = 0.0341394;
if(x>=6.5 && x<6.6) weight = 0.0342169;
if(x>=6.6 && x<6.7) weight = 0.030784;
if(x>=6.7 && x<6.8) weight = 0.031;
if(x>=6.8 && x<6.9) weight = 0.030086;
if(x>=6.9 && x<7) weight = 0.0221675;
if(x>=7 && x<7.1) weight = 0.0309735;
if(x>=7.1 && x<7.2) weight = 0.0273119;
if(x>=7.2 && x<7.3) weight = 0.0240286;
if(x>=7.3 && x<7.4) weight = 0.0218718;
if(x>=7.4 && x<7.5) weight = 0.0142717;
if(x>=7.5 && x<7.6) weight = 0.0181543;
if(x>=7.6 && x<7.7) weight = 0.0139032;
if(x>=7.7 && x<7.79999) weight = 0.0173096;
if(x>=7.79999 && x<7.89999) weight = 0.0123775;
if(x>=7.89999 && x<7.99999) weight = 0.0175077;
if(x>=7.99999 && x<8.09999) weight = 0.0126449;
if(x>=8.09999 && x<8.2) weight = 0.0163016;
if(x>=8.2 && x<8.3) weight = 0.0142204;
if(x>=8.3 && x<8.4) weight = 0.0105152;
if(x>=8.4 && x<8.5) weight = 0.0102197;
if(x>=8.5 && x<8.6) weight = 0.0100636;
if(x>=8.6 && x<8.7) weight = 0.00938478;
if(x>=8.7 && x<8.8) weight = 0.0120421;
if(x>=8.8 && x<8.9) weight = 0.00873138;
if(x>=8.9 && x<9) weight = 0.0123203;
if(x>=9 && x<9.1) weight = 0.00819672;
if(x>=9.1 && x<9.2) weight = 0.0100582;
if(x>=9.2 && x<9.3) weight = 0.00680985;
if(x>=9.3 && x<9.4) weight = 0.00843882;
if(x>=9.4 && x<9.5) weight = 0.00854701;
if(x>=9.5 && x<9.6) weight = 0.00582011;
if(x>=9.6 && x<9.7) weight = 0.00918919;
if(x>=9.7 && x<9.8) weight = 0.00591398;
if(x>=9.8 && x<9.9) weight = 0.00886802;
if(x>=9.9 && x<10) weight = 0.00475436;
if(x>=10 && x<10.1) weight = 0.0042735;
if(x>=10.1 && x<10.2) weight = 0.00632244;
if(x>=10.2 && x<10.3) weight = 0.00537346;
if(x>=10.3 && x<10.4) weight = 0.00316289;
if(x>=10.4 && x<10.5) weight = 0.00656455;
if(x>=10.5 && x<10.6) weight = 0.00428266;
if(x>=10.6 && x<10.7) weight = 0.00485699;
if(x>=10.7 && x<10.8) weight = 0.00428266;
if(x>=10.8 && x<10.9) weight = 0.00374532;
if(x>=10.9 && x<11) weight = 0.00221117;
if(x>=11 && x<11.1) weight = 0.00272331;
if(x>=11.1 && x<11.2) weight = 0.00279642;
if(x>=11.2 && x<11.3) weight = 0.00271592;
if(x>=11.3 && x<11.4) weight = 0.00208877;
if(x>=11.4 && x<11.5) weight = 0.000517866;
if(x>=11.5 && x<11.6) weight = 0.00326619;
if(x>=11.6 && x<11.7) weight = 0.0015748;
if(x>=11.7 && x<11.8) weight = 0.00170358;
if(x>=11.8 && x<11.9) weight = 0.00164204;
if(x>=11.9 && x<12) weight = 0.00286041;
if(x>=12 && x<12.1) weight = 0.0010929;
if(x>=12.1 && x<12.2) weight = 0.0021645;
if(x>=12.2 && x<12.3) weight = 0.000534474;
if(x>=12.3 && x<12.4) weight = 0.00105932;
if(x>=12.4 && x<12.5) weight = 0.00105097;
if(x>=12.5 && x<12.6) weight = 0.000534188;
if(x>=12.6 && x<12.7) weight = 0.00168919;
if(x>=12.7 && x<12.8) weight = 0.00108814;
if(x>=12.8 && x<12.9) weight = 0.00156576;
if(x>=12.9 && x<13) weight = 0;
if(x>=13 && x<13.1) weight = 0.00214362;
if(x>=13.1 && x<13.2) weight = 0.00160858;
if(x>=13.2 && x<13.3) weight = 0.00161031;
if(x>=13.3 && x<13.4) weight = 0.00107817;
if(x>=13.4 && x<13.5) weight = 0.00107066;
if(x>=13.5 && x<13.6) weight = 0.00168919;
if(x>=13.6 && x<13.7) weight = 0.00109111;
if(x>=13.7 && x<13.8) weight = 0.00107411;
if(x>=13.8 && x<13.9) weight = 0.00165198;
if(x>=13.9 && x<14) weight = 0.00109589;
if(x>=14 && x<14.1) weight = 0.000531915;
if(x>=14.1 && x<14.2) weight = 0.00108578;
if(x>=14.2 && x<14.3) weight = 0.00221117;
if(x>=14.3 && x<14.4) weight = 0.000558347;
if(x>=14.4 && x<14.5) weight = 0.00103306;
if(x>=14.5 && x<14.6) weight = 0.00108696;
if(x>=14.6 && x<14.7) weight = 0;
if(x>=14.7 && x<14.8) weight = 0;
if(x>=14.8 && x<14.9) weight = 0.00110011;
if(x>=14.9 && x<15) weight = 0;
if(x>=15) weight = 0;
		
		}
		
		return weight;
	
	}
}

//=======================================================================
Bool_t AliAnalysisHFETPCTOF::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
//Check single track cuts for a given cut step
//Note this function is called inside the UserExec function
	const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
	if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
	return kTRUE;
}

//=======================================================================
Bool_t AliAnalysisHFETPCTOF::FindMother(Int_t mcIndex)
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
		
		if(mpdg>400 && mpdg<500) //charmed mesons (many kind of mesons D)
		{
			if((gmpdg>500 && gmpdg<600) || (ggmpdg>500 && ggmpdg<600) || (gggmpdg>500 && gggmpdg<600)) //when the D comes from B
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
		else if(mpdg>500 && mpdg<600) //bottom mesons (many kind of mesons B)
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
