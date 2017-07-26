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

#include "AliAnalysisUtils.h"
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
#include "AliAnalysisHFETPCTOF_2.h"
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
ClassImp(AliAnalysisHFETPCTOF_2)

//______________________________________________________________________
AliAnalysisHFETPCTOF_2::AliAnalysisHFETPCTOF_2(const char *name) 
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
,fHadrons(0)
,fPi0w(0)
,fPi0w2(0)
,fPi0w3(0)
,fEtaw(0)
,fEtaw2(0)
,fEtaw3(0)
,fMass(0)
,fMinPt(0)
,fTpcNclusAsso(0)
,ftpcPIDmincut(0)
,ftpcPIDmaxcut(0)
,ftofPIDmincut(0)
,ftofPIDmaxcut(0)
,fEtaMin(0)
,fEtaMax(0)

//Elienos class for nonHFE reconstruction by IM

,fMassCutFlag(kTRUE)
,fAngleCutFlag(kFALSE)
,fAngleCut(999)
,fChi2CutFlag(kTRUE)
,fChi2Cut(4.0)
,fDCAcutFlag(kFALSE)
,fDCAcut(999)//dca between two tracks
,fPartnerCuts(new AliESDtrackCuts())
,fTpcNclsAsso(60)
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
,fPtElec_ULS_MC2(0)
,fPtElec_LS_MC2(0)
,fPtElec_ULS_MC_from_pi0(0)
,fPtElec_LS_MC_from_pi0(0)

,fPtElec_ULS_MC_from_pi02(0)
,fPtElec_LS_MC_from_pi02(0)

,fPtElec_ULS_MC_from_eta2(0)
,fPtElec_LS_MC_from_eta2(0)

,fPtElec_ULS_MC_from_gamma2(0)
,fPtElec_LS_MC_from_gamma2(0)

,fPtElec(0)
,fPElec(0)
,fPHad_f(0)
,fPtHad_f(0)

//Histograms for the analysis
,fVertex1(0)
,fNevent(0)
,fNeventT0(0)
,fNevent_3(0)
,fNevent_T0b(0)
,fNevent_corrcut(0)
,fNevent_no_vertex(0)
,fNevent_no_vertex_2(0)
//,fNevent2(0)
,fCent(0)
,fCent2(0)
,fTPC_p1(0)
,fTPC_p2(0)
,fTPC_p3(0)
,fPt_1(0)
,fPt_2(0)
,fITSnClus_1(0)
,fTPCnClus_1(0)
,fITSnClus_2(0)
,fTPCnClus_2(0)
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
,fPt_elec_MB(0)
,fPt_elec_MB_ULS(0)
,fPt_elec_MB_LS(0)

,fPt_elec_pi0_MB2(0)
,fPt_elec_pi0_MB_ULS2(0)
,fPt_elec_pi0_MB_LS2(0)

,fPt_elec_eta_MB2(0)
,fPt_elec_eta_MB_ULS2(0)
,fPt_elec_eta_MB_LS2(0)

,fPt_elec_gamma_MB(0)
,fPt_elec_gamma_MB_ULS(0)
,fPt_elec_gamma_MB_LS(0)

,fPt_elec_gamma_MB2(0)
,fPt_elec_gamma_MB_ULS2(0)
,fPt_elec_gamma_MB_LS2(0)

,fPt_elec_phot(0)
,fPt_elec_phot2(0)
,fPt_elec_from_pi02(0)
,fPt_elec_from_eta2(0)
,fPt_elec_from_gamma2(0)
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
,fPtMCeta(0)
,fPtMCeta_MB_afterNorm(0)
,fPtMCpi03(0)
,fPtMCeta3(0)
,fPtHFEMC_2(0)
,fPi0EtaSpectra(0)


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
AliAnalysisHFETPCTOF_2::AliAnalysisHFETPCTOF_2() 
  : AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisHFETPCTOF_2")

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
,fHadrons(0)
,fPi0w(0)
,fPi0w2(0)
,fPi0w3(0)
,fEtaw(0)
,fEtaw2(0)
,fEtaw3(0)
,fMass(0)
,fMinPt(0)
,fTpcNclusAsso(0)
,ftpcPIDmincut(0)
,ftpcPIDmaxcut(0)
,ftofPIDmincut(0)
,ftofPIDmaxcut(0)
,fEtaMin(0)
,fEtaMax(0)

//Elienos class for nonHFE reconstruction by IM
,fMassCutFlag(kTRUE)
,fAngleCutFlag(kFALSE)
,fAngleCut(999)
,fChi2CutFlag(kTRUE)
,fChi2Cut(4.0)
,fDCAcutFlag(kFALSE)
,fDCAcut(999)//dca between two tracks
,fPartnerCuts(new AliESDtrackCuts())
,fTpcNclsAsso(60)
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
,fPtElec_ULS_MC2(0)
,fPtElec_LS_MC2(0)
,fPtElec_ULS_MC_from_pi0(0)
,fPtElec_LS_MC_from_pi0(0)

,fPtElec_ULS_MC_from_pi02(0)
,fPtElec_LS_MC_from_pi02(0)

,fPtElec_ULS_MC_from_eta2(0)
,fPtElec_LS_MC_from_eta2(0)

,fPtElec_ULS_MC_from_gamma2(0)
,fPtElec_LS_MC_from_gamma2(0)

,fPtElec(0)
,fPElec(0)
,fPHad_f(0)
,fPtHad_f(0)

//Histograms for the analysis
,fVertex1(0)
,fNevent(0)
,fNeventT0(0)
,fNevent_3(0)
,fNevent_T0b(0)
,fNevent_corrcut(0)
,fNevent_no_vertex(0)
,fNevent_no_vertex_2(0)
,fCent(0)
,fCent2(0)
,fTPC_p1(0)
,fTPC_p2(0)
,fTPC_p3(0)
,fPt_1(0)
,fPt_2(0)
,fITSnClus_1(0)
,fTPCnClus_1(0)
,fITSnClus_2(0)
,fTPCnClus_2(0)
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
,fPt_elec_MB(0)
,fPt_elec_MB_ULS(0)
,fPt_elec_MB_LS(0)

,fPt_elec_pi0_MB2(0)
,fPt_elec_pi0_MB_ULS2(0)
,fPt_elec_pi0_MB_LS2(0)

,fPt_elec_eta_MB2(0)
,fPt_elec_eta_MB_ULS2(0)
,fPt_elec_eta_MB_LS2(0)

,fPt_elec_gamma_MB(0)
,fPt_elec_gamma_MB_ULS(0)
,fPt_elec_gamma_MB_LS(0)

,fPt_elec_gamma_MB2(0)
,fPt_elec_gamma_MB_ULS2(0)
,fPt_elec_gamma_MB_LS2(0)

,fPt_elec_phot(0)
,fPt_elec_phot2(0)
,fPt_elec_from_pi02(0)
,fPt_elec_from_eta2(0)
,fPt_elec_from_gamma2(0)
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
,fPtMCeta(0)
,fPtMCpi03(0)
,fPtMCeta3(0)
,fPtHFEMC_2(0)
,fPi0EtaSpectra(0)


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
AliAnalysisHFETPCTOF_2::~AliAnalysisHFETPCTOF_2()
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
void AliAnalysisHFETPCTOF_2::UserCreateOutputObjects()
{
//______________________________________________________________________
//Initialize PID
	if(!fPID->GetNumberOfPIDdetectors()) 
    {
		fPID->AddDetector("TPC", 0);
		fPID->AddDetector("TOF", 1);
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
	fNeventT0 = new TH1F("fNeventT0","Number of Events",5,-0.5,4.5);
	fVertex1 = new TH1F("fVertex1","vertex",500,-50,50);
	fNevent_3 = new TH1F("fNevent_3","Number of Events",5,-0.5,4.5);
	fNevent_T0b = new TH1F("fNevent_T0b","Number of Events",5,-0.5,4.5);
	fNevent_corrcut = new TH1F("fNevent_corrcut","Number of Events",5,-0.5,4.5);
	fNevent_no_vertex = new TH1F("fNevent_no_vertex","Number of Events",5,-0.5,4.5);
	fNevent_no_vertex_2 = new TH1F("fNevent_no_vertex_2","Number of Events",5,-0.5,4.5);
	fCent = new TH1F("fCent","Centrality",408,-1,101);
	fCent2 = new TH1F("fCent2","Centrality",408,-1,101);
	//And then, add to the output list
	fOutputList->Add(fNevent);
	fOutputList->Add(fNeventT0);
	fOutputList->Add(fVertex1);
	fOutputList->Add(fNevent_3);
	fOutputList->Add(fNevent_T0b);
	fOutputList->Add(fNevent_corrcut);
	fOutputList->Add(fNevent_no_vertex);
	fOutputList->Add(fNevent_no_vertex_2);
	fOutputList->Add(fCent);
	fOutputList->Add(fCent2);
	
	//General Histograms
	
	//Steps
	//Step 1: Before Track cuts
	//Step 2: Before PID
	//Step 3: After PID

	
	Double_t ptbinning[25] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5};

	fTPC_p1 = new TH2F("fTPC_p1","p (GeV/c);TPC dE/dx (a. u.)",300,0,15,400,-20,200);
	fOutputList->Add(fTPC_p1);
	
	fTPC_p2 = new TH2F("fTPC_p2","p (GeV/c);TPC dE/dx (a. u.)",300,0,15,400,-20,200);
	fOutputList->Add(fTPC_p2);
	
	fTPC_p3 = new TH2F("fTPC_p3","p (GeV/c);TPC dE/dx (a. u.)",300,0,15,400,-20,200);
	fOutputList->Add(fTPC_p3);
	
	fPt_1 = new TH1F("fPt_1","pt (GeV/c)",1000,0,100);
	fOutputList->Add(fPt_1);
	
	fPt_2 = new TH1F("fPt_2","pt (GeV/c)",1000,0,100);
	fOutputList->Add(fPt_2);
	
	fITSnClus_1 = new TH1F("fITSnClus_1","fITSnClus_1",1000,0,80);
	fOutputList->Add(fITSnClus_1);
	
	fITSnClus_2 = new TH1F("fITSnClus_2","fITSnClus_2",1000,0,80);
	fOutputList->Add(fITSnClus_2);
	
	fTPCnClus_1 = new TH1F("fTPCnClus_1","fTPCnClus_1",1000,0,400);
	fOutputList->Add(fTPCnClus_1);
	
	fTPCnClus_2 = new TH1F("fTPCnClus_2","fTPCnClus_2",1000,0,400);
	fOutputList->Add(fTPCnClus_2);	
		
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
	
	fPtElec_ULS = new TH1F("fPtElec_ULS","ULS; p_{T} [GeV/c]; Count",24,ptbinning);
	fOutputList->Add(fPtElec_ULS);
	
	//fPtElec_LS = new TH1F("fPtElec_LS","LS; p_{T} [GeV/c]; Count",24,ptbinning);
	fPtElec_LS = new TH1F("fPtElec_LS","LS; p_{T} [GeV/c]; Count",24,ptbinning);
    fOutputList->Add(fPtElec_LS);
    
    fTPCnSigmaPElec_ULS = new TH2F("fTPCnSigmaPElec_ULS","p (GeV/c);TPC Electron N#sigma",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnSigmaPElec_ULS);   
	
	fTPCnSigmaPElec_LS = new TH2F("fTPCnSigmaPElec_LS","p (GeV/c);TPC Electron N#sigma",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnSigmaPElec_LS);   
    
    fPtElec = new TH1F("fPtElec","; p_{T} [GeV/c]; Count",24,ptbinning);
    fOutputList->Add(fPtElec);
    
    fPElec = new TH1F("fPElec","; p [GeV/c]; Count",24,ptbinning);
    fOutputList->Add(fPElec);
    
  
     
    fPtHad_f = new TH1F("fPtHad_f","; p_{T} [GeV/c]; Count",24,ptbinning);
    fOutputList->Add(fPtHad_f);
    fPtHad_f->Sumw2();
    
    fPHad_f = new TH1F("fPHad_f","; p [GeV/c]; Count",24,ptbinning);
    fOutputList->Add(fPHad_f);
    fPHad_f->Sumw2();
    
   
    fDCAxy_pt = new TH2F("fDCAxy_pt",";p_{t} (GeV/c);DCAxy ",300,0,30,200,-3,3);
    fOutputList->Add(fDCAxy_pt);
    
    fDCAz_pt = new TH2F("fDCAz_pt",";p_{t} (GeV/c);DCAz ",300,0,30,200,-3,3);
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
    
    fPt_elec_phot = new TH1F("fPt_elec_phot","",24,ptbinning);
	fOutputList->Add(fPt_elec_phot);
	fPt_elec_phot->Sumw2();
	
	fPt_elec_phot2 = new TH1F("fPt_elec_phot2","",24,ptbinning);
	fOutputList->Add(fPt_elec_phot2);
	fPt_elec_phot2->Sumw2();
	
	
	fPt_elec_pi0_MB2 = new TH1F("fPt_elec_pi0_MB2","",24,ptbinning);
	fOutputList->Add(fPt_elec_pi0_MB2);
	fPt_elec_pi0_MB2->Sumw2();
	
	fPt_elec_pi0_MB_ULS2 = new TH1F("fPt_elec_pi0_MB_ULS2","",24,ptbinning);
	fOutputList->Add(fPt_elec_pi0_MB_ULS2);
	fPt_elec_pi0_MB_ULS2->Sumw2();
	
	fPt_elec_pi0_MB_LS2 = new TH1F("fPt_elec_pi0_MB_LS2","",24,ptbinning);
	fOutputList->Add(fPt_elec_pi0_MB_LS2);
	fPt_elec_pi0_MB_LS2->Sumw2();
	
	
	fPt_elec_MB = new TH1F("fPt_elec_MB","",24,ptbinning);
	fOutputList->Add(fPt_elec_MB);
	fPt_elec_MB->Sumw2();
	
	fPt_elec_MB_ULS = new TH1F("fPt_elec_MB_ULS","",24,ptbinning);
	fOutputList->Add(fPt_elec_MB_ULS);
	fPt_elec_MB_ULS->Sumw2();
	
	fPt_elec_MB_LS = new TH1F("fPt_elec_MB_LS","",24,ptbinning);
	fOutputList->Add(fPt_elec_MB_LS);
	fPt_elec_MB_LS->Sumw2();
	
	fPt_elec_eta_MB2 = new TH1F("fPt_elec_eta_MB2","",24,ptbinning);
	fOutputList->Add(fPt_elec_eta_MB2);
	fPt_elec_eta_MB2->Sumw2();
	
	fPt_elec_eta_MB_ULS2 = new TH1F("fPt_elec_eta_MB_ULS2","",24,ptbinning);
	fOutputList->Add(fPt_elec_eta_MB_ULS2);
	fPt_elec_eta_MB_ULS2->Sumw2();
	
	fPt_elec_eta_MB_LS2 = new TH1F("fPt_elec_eta_MB_LS2","",24,ptbinning);
	fOutputList->Add(fPt_elec_eta_MB_LS2);
	fPt_elec_eta_MB_LS2->Sumw2();
	
	
	
	fPt_elec_gamma_MB = new TH1F("fPt_elec_gamma_MB","",24,ptbinning);
	fOutputList->Add(fPt_elec_gamma_MB);
	fPt_elec_gamma_MB->Sumw2();
	
	fPt_elec_gamma_MB_ULS = new TH1F("fPt_elec_gamma_MB_ULS","",24,ptbinning);
	fOutputList->Add(fPt_elec_gamma_MB_ULS);
	fPt_elec_gamma_MB_ULS->Sumw2();
	
	fPt_elec_gamma_MB_LS = new TH1F("fPt_elec_gamma_MB_LS","",24,ptbinning);
	fOutputList->Add(fPt_elec_gamma_MB_LS);
	fPt_elec_gamma_MB_LS->Sumw2();
	
	fPt_elec_gamma_MB2 = new TH1F("fPt_elec_gamma_MB2","",24,ptbinning);
	fOutputList->Add(fPt_elec_gamma_MB2);
	fPt_elec_gamma_MB2->Sumw2();
	
	fPt_elec_gamma_MB_ULS2 = new TH1F("fPt_elec_gamma_MB_ULS2","",24,ptbinning);
	fOutputList->Add(fPt_elec_gamma_MB_ULS2);
	fPt_elec_gamma_MB_ULS2->Sumw2();
	
	fPt_elec_gamma_MB_LS2 = new TH1F("fPt_elec_gamma_MB_LS2","",24,ptbinning);
	fOutputList->Add(fPt_elec_gamma_MB_LS2);
	fPt_elec_gamma_MB_LS2->Sumw2();
	
	
	fPt_elec_from_pi02 = new TH1F("fPt_elec_from_pi02","",24,ptbinning);
	fOutputList->Add(fPt_elec_from_pi02);
	fPt_elec_from_pi02->Sumw2();
	
	fPt_elec_from_eta2 = new TH1F("fPt_elec_from_eta2","",24,ptbinning);
	fOutputList->Add(fPt_elec_from_eta2);
	fPt_elec_from_eta2->Sumw2();
	

	fPt_elec_from_gamma2 = new TH1F("fPt_elec_from_gamma2","",24,ptbinning);
	fOutputList->Add(fPt_elec_from_gamma2);
	fPt_elec_from_gamma2->Sumw2();
	
	fPt_HFelec_passTOFandTPC_2 = new TH1F("fPt_HFelec_passTOFandTPC_2","",24,ptbinning);
	fOutputList->Add(fPt_HFelec_passTOFandTPC_2);
	
	fPt_HFelec_passTPC_2 = new TH1F("fPt_HFelec_passTPC_2","",24,ptbinning);
	fOutputList->Add(fPt_HFelec_passTPC_2);
	
	fPt_HFelec_pass_track = new TH1F("fPt_HFelec_pass_track","",24,ptbinning);
	fOutputList->Add(fPt_HFelec_pass_track);
	
	fPt_HFelec_pass_trackDCA = new TH1F("fPt_HFelec_pass_trackDCA","",24,ptbinning);
	fOutputList->Add(fPt_HFelec_pass_trackDCA);
	
	fPt_HFelec_pass_trackITSlayer = new TH1F("fPt_HFelec_pass_trackITSlayer","",24,ptbinning);
	fOutputList->Add(fPt_HFelec_pass_trackITSlayer);
	
	fPt_HFelec_pass_trackFB = new TH1F("fPt_HFelec_pass_trackFB","",24,ptbinning);
	fOutputList->Add(fPt_HFelec_pass_trackFB);
	
	fPt_HFelec_pass_trackTPCITS = new TH1F("fPt_HFelec_pass_trackTPCITS","",24,ptbinning);
	fOutputList->Add(fPt_HFelec_pass_trackTPCITS);
	
	fPt_HFelec_pass_trackKink = new TH1F("fPt_HFelec_pass_trackKink","",24,ptbinning);
	fOutputList->Add(fPt_HFelec_pass_trackKink);
	
	
	fPt_HFelec_pass_trackandTOF = new TH1F("fPt_HFelec_pass_trackandTOF","",24,ptbinning);
	fOutputList->Add(fPt_HFelec_pass_trackandTOF);
		
			
	fPtMCeta = new TH1F("fPtMCeta",";p_{t} (GeV/c)",2000,0,100);
	fOutputList->Add(fPtMCeta);
	
	fPtMCeta_MB_afterNorm = new TH1F("fPtMCeta_MB_afterNorm",";p_{t} (GeV/c)",2000,0,100);
	fOutputList->Add(fPtMCeta_MB_afterNorm);
	fPtMCeta_MB_afterNorm->Sumw2();
	
			
	fPtMCpi03 = new TH1F("fPtMCpi03",";p_{t} (GeV/c)",2000,0,100);
	fOutputList->Add(fPtMCpi03);
	
	fPtMCeta3 = new TH1F("fPtMCeta3",";p_{t} (GeV/c)",2000,0,100);
	fOutputList->Add(fPtMCeta3);
	
	fPtElec_ULS_MC = new TH1F("fPtElec_ULS_MC","p_{T} [GeV/c]",24,ptbinning);
	fOutputList->Add(fPtElec_ULS_MC);
	fPtElec_ULS_MC->Sumw2();
	
	fPtElec_ULS_MC2 = new TH1F("fPtElec_ULS_MC2","p_{T} [GeV/c]",24,ptbinning);
	fOutputList->Add(fPtElec_ULS_MC2);
	fPtElec_ULS_MC2->Sumw2();
	
		
	fPtElec_LS_MC = new TH1F("fPtElec_LS_MC","p_{T} [GeV/c]",24,ptbinning);
    fOutputList->Add(fPtElec_LS_MC);
	fPtElec_LS_MC->Sumw2();
	
	fPtElec_LS_MC2 = new TH1F("fPtElec_LS_MC2","p_{T} [GeV/c]",24,ptbinning);
    fOutputList->Add(fPtElec_LS_MC2);
	fPtElec_LS_MC2->Sumw2();

	
	fPtElec_ULS_MC_from_pi0 = new TH1F("fPtElec_ULS_MC_from_pi0","p_{T} [GeV/c]",24,ptbinning);
	fOutputList->Add(fPtElec_ULS_MC_from_pi0);
	fPtElec_ULS_MC_from_pi0->Sumw2();
	
	fPtElec_LS_MC_from_pi0 = new TH1F("fPtElec_LS_MC_from_pi0","p_{T} [GeV/c]",24,ptbinning);
    fOutputList->Add(fPtElec_LS_MC_from_pi0);
	fPtElec_LS_MC_from_pi0->Sumw2();
	
	fPtElec_ULS_MC_from_pi02 = new TH1F("fPtElec_ULS_MC_from_pi02","p_{T} [GeV/c]",24,ptbinning);
	fOutputList->Add(fPtElec_ULS_MC_from_pi02);
	fPtElec_ULS_MC_from_pi02->Sumw2();
	
	fPtElec_LS_MC_from_pi02 = new TH1F("fPtElec_LS_MC_from_pi02","p_{T} [GeV/c]",24,ptbinning);
    fOutputList->Add(fPtElec_LS_MC_from_pi02);
	fPtElec_LS_MC_from_pi02->Sumw2();
	
	
	fPtElec_ULS_MC_from_eta2 = new TH1F("fPtElec_ULS_MC_from_eta2","p_{T} [GeV/c]",24,ptbinning);
	fOutputList->Add(fPtElec_ULS_MC_from_eta2);
	fPtElec_ULS_MC_from_eta2->Sumw2();
	
	fPtElec_LS_MC_from_eta2 = new TH1F("fPtElec_LS_MC_from_eta2","p_{T} [GeV/c]",24,ptbinning);
    fOutputList->Add(fPtElec_LS_MC_from_eta2);
	fPtElec_LS_MC_from_eta2->Sumw2();
		
	fPtElec_ULS_MC_from_gamma2 = new TH1F("fPtElec_ULS_MC_from_gamma2","p_{T} [GeV/c]",24,ptbinning);
	fOutputList->Add(fPtElec_ULS_MC_from_gamma2);
	fPtElec_ULS_MC_from_gamma2->Sumw2();
	
	fPtElec_LS_MC_from_gamma2 = new TH1F("fPtElec_LS_MC_from_gamma2","p_{T} [GeV/c]",24,ptbinning);
    fOutputList->Add(fPtElec_LS_MC_from_gamma2);
	fPtElec_LS_MC_from_gamma2->Sumw2();
	
	fPtHFEMC_2 = new TH1F("fPtHFEMC_2",";p_{t} (GeV/c)",24,ptbinning);
	fOutputList->Add(fPtHFEMC_2);
	
	fTPCnsigma_p_after_tof = new TH2F("fTPCnsigma_p_after_tof","p (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnsigma_p_after_tof);
	
	fTPCnsigma_p_after_tof_its = new TH2F("fTPCnsigma_p_after_tof_its","p (GeV/c);TPC Electron N#sigma after TOF and ITS cuts",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnsigma_p_after_tof_its);
	
	
	fTPCnsigma_pt_after_tof = new TH2F("fTPCnsigma_pt_after_tof","pt (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnsigma_pt_after_tof);
    
    fTPCnsigma_pt_after_tof_its = new TH2F("fTPCnsigma_pt_after_tof_its","pt (GeV/c);TPC Electron N#sigma after TOF and ITS cuts",300,0,15,200,-15,10);
	fOutputList->Add(fTPCnsigma_pt_after_tof_its);
	
	
	// ----- weights for tagging efficiency -----
    // QA weights
    Int_t nBinspdg = 3;
    Double_t minpdg = 0.;
    Double_t maxpdg = 3.;
    Double_t binLimpdg[nBinspdg+1];
    for(Int_t i=0; i<=nBinspdg; i++) binLimpdg[i]=(Double_t)minpdg + (maxpdg-minpdg)/nBinspdg*(Double_t)i ;
    
    Int_t nBinsg = 2;
    Double_t ming = 0.;
    Double_t maxg = 2.;
    Double_t binLimg[nBinsg+1];
    for(Int_t i=0; i<=nBinsg; i++) binLimg[i]=(Double_t)ming + (maxg-ming)/nBinsg*(Double_t)i ;
    
    Int_t nBinstype = 7;
    Double_t mintype = -1.;
    Double_t maxtype = 6.;
    Double_t binLimtype[nBinstype+1];
    for(Int_t i=0; i<=nBinstype; i++) binLimtype[i]=(Double_t)mintype + (maxtype-mintype)/nBinstype*(Double_t)i ;
    
    Int_t nBinspt = 44;
    Double_t binLimpt[45]= {0.1,0.112797,0.127231,0.143512,0.161877,0.182592,0.205957,0.232313,0.262041,0.295573,0.333397,
0.37606,0.424183,0.478465,0.539692,0.608754,0.686654,0.774523,0.873636,0.985432,1.11153,1.25377,1.41421,1.59519,1.79932,2.02957,
2.28928,2.58223,2.91267,3.2854,3.70582,4.18004,4.71494,5.3183,5.99886,6.76651,7.6324,8.60909,9.71076,10.9534,12.3551,13.9361,15.7195,17.731,20};//bin limits from the measured pi0 spectrum
    //for(Int_t i=0; i<=nBinspt; i++) binLimpt[i]=(Double_t)minpt + (maxpt-minpt)/nBinspt*(Double_t)i ;
    
    const Int_t nDima=4;
    Int_t nBina[nDima] = {44,nBinspdg,nBinsg,nBinstype};
    fPi0EtaSpectra = new THnSparseF("Pi0EtaSpectra","Pi0EtaSpectra",nDima,nBina);
    fPi0EtaSpectra->SetBinEdges(0,binLimpt);
    fPi0EtaSpectra->SetBinEdges(1,binLimpdg);
    fPi0EtaSpectra->SetBinEdges(2,binLimg);
    fPi0EtaSpectra->SetBinEdges(3,binLimtype);
    fPi0EtaSpectra->Sumw2();
    //fOutputList->Add(fPi0EtaSpectra); ---> put at the end of the list
    // ------------------------------------------
    fOutputList->Add(fPi0EtaSpectra);
	
	
    
    
//______________________________________________________________________
	
	PostData(1, fOutputList);
	
///______________________________________________________________________
}




//______________________________________________________________________
//Main loop
//Called for each event
void AliAnalysisHFETPCTOF_2::UserExec(Option_t *) 
{


	Int_t pdg = -99999;
	Int_t pdg_mother = -99999;
	Double_t weight = -99999;
	Double_t weight_had_f = -99999;
	Double_t weight_pi0 = -99999;
	Double_t weight_eta = -99999;
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
	Int_t fITSnClus = 99999;
	Int_t fTPCnClus = 99999; 
	Double_t qaweights[5];

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
				
			
		///Checking Nv0 and NT0 before the vertex selection--------------
		fNevent->Fill(0);///triggered by v0
		
		///also triggered by T0
		if(fIsPP){
			if(!fIsMC){ 
		
				TString firedTriggerClasses = static_cast<const AliAODEvent*>(InputEvent())->GetFiredTriggerClasses();
				if ((firedTriggerClasses.Contains("C0TVX-B-NOPF-CENT"))){
					fNevent_T0b->Fill(0);
					
				}
		
			}
		}
		///-------------------------------------------------------------

		AliAODVertex* vtTrc = fAOD->GetPrimaryVertex();
				
		///Events with no vertex by tracks-------------------------------		
		if(!vtTrc){
			 fNevent_no_vertex->Fill(0);
			 return;	
		}
		TString vtxTtl = vtTrc->GetTitle();
		if(!vtxTtl.Contains("VertexerTracks")) fNevent_no_vertex_2->Fill(0); 
		///-------------------------------------------------------------
		
		///checking the minimum number of contributors (this need to be done in pp, Pb-Pb, p-Pb...)
		if(vtTrc->GetNContributors()<2) return;
		
		
		Bool_t isPileupfromSPDmulbins=fAOD->IsPileupFromSPDInMultBins();
		
		///minContributors=5; minChi2=5.; minWeiZDiff=15; checkPlpFromDifferentBC=kFALSE;
		//AliAnalysisUtils *fUtils;
		AliAnalysisUtils utils;
		utils.SetMinPlpContribMV(5);
		utils.SetMaxPlpChi2MV(5.);
		utils.SetMinWDistMV(15);
		utils.SetCheckPlpFromDifferentBCMV(kFALSE);
		Bool_t isPileupFromMV = utils.IsPileUpMV(fAOD);
		
		
		///Number of events rejected by the correlation cuts between SPD and track vertexes
		fNevent_corrcut->Fill(0);
		
		///Track vertex cut:
		Float_t zvtx = vtTrc->GetZ();
		if(TMath::Abs(zvtx) > 10) return;
		fVertex1->Fill(zvtx);
		
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
		
		
		///-------------------------------------------------------------
		///Checking Nv0 and NT0 after the vertex selection
		fNevent_3->Fill(0); ///Number of analysed events
		
		///Number of events in kint7 that is also triggered by T0
		if(fIsPP){
			if(!fIsMC){ 
		
				TString firedTriggerClasses = static_cast<const AliAODEvent*>(InputEvent())->GetFiredTriggerClasses();
				if ((firedTriggerClasses.Contains("C0TVX-B-NOPF-CENT"))){
					fNeventT0->Fill(0);
					
				}
		
			}
		}
		///-------------------------------------------------------------
		
		
		
	}
	
	
//Centrality selection
	if(!fIsPP){

//______________________________________________________________________	
	
//Only events with at least 2 tracks are accepted
	Int_t fNOtrks =  fVevent->GetNumberOfTracks();
	if(fNOtrks<2) return;
//______________________________________________________________________


		//if(!fIsMC){     
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
    
			//if(cent<0.5 || cent>10.5) return;
			if(cent<30.5 || cent>50.5) return; //Comparing with Denise
			
			fCent2->Fill(cent);
		//}
	}
//______________________________________________________________________

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
			 if(fMCarray->GetEntries() < 1) return;

			 
			 weight = 1;
			 for(Int_t iMC = 0; iMC < fMCarray->GetEntries(); iMC++)
			 {
				 fMCparticle = (AliAODMCParticle*) fMCarray->At(iMC);
				 pdg = fMCparticle->GetPdgCode();
				 
                 
                 ///For the reconstruction efficiency of HFE:-----------
				 if(fMCparticle->Eta()>=fEtaMin && fMCparticle->Eta()<=fEtaMax)
				 {					
					if(TMath::Abs(pdg) == 11){
						
						Bool_t MotherFound = FindMother(iMC);
							if(MotherFound){
								if(fIsHFE1){
									fPtHFEMC_2->Fill(fMCparticle->Pt());
								}
							}
					}
				 }
				 ///----------------------------------------------------
				     
				///pt spectra for pi0 and eta
				if(fMCparticle->Eta()<-1.2 || fMCparticle->Eta()>1.2) continue;
				    
				///Using thnSparse--------------------------------------
				//Pt
				qaweights[0] = fMCparticle->Pt();
						
				// What pdg
				qaweights[1]=-1.;
				if (TMath::Abs(fMCparticle->GetPdgCode())==111) qaweights[1]=0.2;  // pdg=111 pi0
				if (TMath::Abs(fMCparticle->GetPdgCode())==221) qaweights[1]=1.2;  // pdg=221 eta
				if (TMath::Abs(fMCparticle->GetPdgCode())==22) qaweights[1]=2.2;   // pdg=22  gamma
						
				// What type
				Int_t type = GetPi0EtaType(fMCparticle,fMCarray);
				qaweights[3]=type;
					
				// Fill
				if(qaweights[1]>0.) fPi0EtaSpectra->Fill(qaweights);
				///-----------------------------------------------------
		
		       			
					

			}///end of loop in MC array
		   
			 
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
		
		///Since I just look at the spectra up to 4.3 GeV
		if(track->Pt()>5) continue;

		///_____________________________________________________________________________
		///Fill QA plots without track selection
		fPt = track->Pt();
		fEta = track->Eta();
		fPhi = track->Phi();
		fP = TMath::Sqrt((track->Pt())*(track->Pt()) + (track->Pz())*(track->Pz()));
		
		fITSnClus =  track->GetNumberOfITSClusters();
		fTPCnClus =  track->GetNumberOfTPCClusters(); 
		
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
        
        fPt_1->Fill(fPt);
        
        fITSnClus_1->Fill(fITSnClus);
		fTPCnClus_1->Fill(fTPCnClus);
        
        if((track->Eta() < fEtaMin) || (track->Eta() > fEtaMax)) continue;
           
		
//=======================================================================
// Track Selection Cuts are applied here
//=======================================================================


//AOD (Test Filter Bit)---------------------------------------------------------------------
		if(fIsAOD)
		{
						
			if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; 
		}
		
		///How many HF electrons pass this cut
		if(fIsMC && fIsAOD){
			Bool_t IsHFEMC = IsHFelectronsMC(track);
			if(IsHFEMC) fPt_HFelec_pass_trackFB->Fill(fPt);
		}
		///-------------------------------------------------------------
		
		
//-----------------------------------------------------------------------------------
		
//RecKine: ITSTPC cuts --------------------------------------------------------------
		if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
		
		///How many HF electrons pass this cut
		if(fIsMC && fIsAOD){
			Bool_t IsHFEMC = IsHFelectronsMC(track);
			if(IsHFEMC) fPt_HFelec_pass_trackTPCITS->Fill(fPt);
		}
		///-------------------------------------------------------------------
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
		
		
		///How many HF electrons pass this cut
		if(fIsMC && fIsAOD){
			Bool_t IsHFEMC = IsHFelectronsMC(track);
			if(IsHFEMC) fPt_HFelec_pass_trackKink->Fill(fPt);
		}
		///--------------------------------------------------------------------
//---------------------------------------------------------------------------------    
    
//RecPrim -------------------------------------------------------------------------
		if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;				///ProcessCutStep(Int_t cutStep, AliVParticle *track)
				
		///How many HF electrons pass this cut
		if(fIsMC && fIsAOD){
			Bool_t IsHFEMC = IsHFelectronsMC(track);
			if(IsHFEMC) fPt_HFelec_pass_trackDCA->Fill(fPt);
		}
		///-------------------------------------------------------------
//-------------------------------------------------------------------------------------		
    
//HFEcuts: ITS layers cuts -------------------------------------------------------------
		if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;
		
		///How many HF electrons pass this cut		
		if(fIsMC && fIsAOD){
			Bool_t IsHFEMC = IsHFelectronsMC(track);
			if(IsHFEMC) fPt_HFelec_pass_trackITSlayer->Fill(fPt);
		}
		///-------------------------------------------------------------
//--------------------------------------------------------------------------------------		
		
//HFE cuts: TPC PID cleanup
		if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;
        
//=======================================================================
// QA plots after track selection
//=======================================================================

		///How many HF electrons pass the cuts (eta and track selection)
		if(fIsMC && fIsAOD){
			Bool_t IsHFEMC = IsHFelectronsMC(track);
			if(IsHFEMC) fPt_HFelec_pass_track->Fill(fPt);
		}
		///------------------------------------------------------------
		
		fTPC_p2->Fill(fP,fTPCsignal);
		
		fTPCnsigma_TOFnsigma2->Fill(fTOFnSigma,fTPCnSigma);
		fTOFnsigma_p2->Fill(fP,fTOFnSigma);
		fTOFnsigma_pt2->Fill(fPt,fTOFnSigma);
		fTPCnsigma_p2->Fill(fP,fTPCnSigma);
		fTPCnsigma_pt2->Fill(fPt,fTPCnSigma);
		fITSnsigma_p2->Fill(fP,fITSnSigma);
		fITSnsigma_pt2->Fill(fPt,fITSnSigma);
		
		fPt_2->Fill(fPt);
		
		fITSnClus_2->Fill(fITSnClus);
		fTPCnClus_2->Fill(fTPCnClus);
		
		///Checking DCA---------	
		Double_t d0z0[2], cov[3];
		AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
		track->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 20., d0z0, cov);
		Double_t DCAxy = d0z0[0];
		Double_t DCAz = d0z0[1];
		            
		fDCAxy_pt->Fill(fPt,DCAxy);
		fDCAz_pt->Fill(fPt,DCAz);
		///---------------------	
			
		
		if(fTOFnSigma >= ftofPIDmincut && fTOFnSigma <= ftofPIDmaxcut){
			fTPCnsigma_p_after_tof->Fill(fP,fTPCnSigma);
			fTPCnsigma_pt_after_tof->Fill(fPt,fTPCnSigma);
			
			if(fITSnSigma >= -2 && fITSnSigma <= 2){
				fTPCnsigma_p_after_tof_its->Fill(fP,fTPCnSigma);
				fTPCnsigma_pt_after_tof_its->Fill(fPt,fTPCnSigma);
			}
			
			///How many HF electrons pass the cuts (TOF and track selection)		
			if(fIsMC && fIsAOD){
				Bool_t IsHFEMC = IsHFelectronsMC(track);
				if(IsHFEMC) fPt_HFelec_pass_trackandTOF->Fill(fPt);
			}
			///----------------------------------------------------------			
		}
						
				
		
		if(fTPCnSigma >= ftpcPIDmincut && fTPCnSigma <= ftpcPIDmaxcut){
			
			///How many HF electrons pass the cuts (TPC and track selection)
			if(fIsMC && fIsAOD){
				Bool_t IsHFEMC = IsHFelectronsMC(track);
				if(IsHFEMC) fPt_HFelec_passTPC_2->Fill(fPt);
			}
			///---------------------------------------------------------
		
			if(fTOFnSigma >= ftofPIDmincut && fTOFnSigma <= ftofPIDmaxcut){
				
				///How many HF electrons pass the cuts (TOF, TPC and track selection)		
				if(fIsMC && fIsAOD){
					Bool_t IsHFEMC = IsHFelectronsMC(track);
					if(IsHFEMC) fPt_HFelec_passTOFandTPC_2->Fill(fPt);
				}
				///-----------------------------------------------------
				
			}	
		}
			
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

	    				
			fTPC_p3->Fill(fP,fTPCsignal);
			fTPCnsigma_p3->Fill(fP,fTPCnSigma);
			fTPCnsigma_pt3->Fill(fPt,fTPCnSigma);
		
			fTPCnsigma_TOFnsigma3->Fill(fTOFnSigma,fTPCnSigma);
			fTOFnsigma_p3->Fill(fP,fTOFnSigma);
			fTOFnsigma_pt3->Fill(fPt,fTOFnSigma);
			
			fITSnsigma_p3->Fill(fP,fITSnSigma);
			fITSnsigma_pt3->Fill(fPt,fITSnSigma);
			
			fPtElec->Fill(fPt);
			fPElec->Fill(fP);
			
			///Hadron Contamination Subtraction-------------------------
			///Calculating the weight via the fitted function:
			weight_had_f = 0;
			weight_had_f = fHadrons->Eval(fP);
			
			fPtHad_f->Fill(fPt,weight_had_f);
			fPHad_f->Fill(fP,weight_had_f);
			///---------------------------------------------------------
		
        
			///InvMass calculation (to find NonHFE):
			///---------------------------------------------------------
			//Associated particle cut
		            
			fNonHFE = new AliSelectNonHFE();
		      
			fNonHFE->SetAODanalysis(fIsAOD);

			if(fMassCutFlag) fNonHFE->SetInvariantMassCut(fMass);
			if(fAngleCutFlag) fNonHFE->SetOpeningAngleCut(fAngleCut);
			if(fChi2CutFlag) fNonHFE->SetChi2OverNDFCut(fChi2Cut);
			if(fDCAcutFlag) fNonHFE->SetDCACut(fDCAcut);
			fNonHFE->SetAlgorithm("DCA"); //KF
			fNonHFE->SetPIDresponse(fPidResponse);
			fNonHFE->SetTrackCuts(-3.5,3.5); 
			fNonHFE->SetAdditionalCuts(fMinPt,fTpcNclusAsso);											
			
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
					
			///---------------------------------------------------------
        
            ///----------------- Tagging efficiency --------------------
           
            if(fIsMC && fIsAOD){      
					            
				fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
					pdg = fMCparticle->GetPdgCode();
					///Is electron:
					if(TMath::Abs(pdg) != 11) continue;
						if(fMCparticle->GetMother()<=0) continue;
							fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
							pdg_mother = fMCparticleMother->GetPdgCode();

							weight=1;
						
							
							///Photonic Electrons:
							if(TMath::Abs(pdg_mother) == 22 || TMath::Abs(pdg_mother) == 111 || TMath::Abs(pdg_mother) == 221){
									
							Double_t ptmotherw = -1.;
							Int_t electronsource = GetElecSourceType(fMCparticle,fMCarray,ptmotherw);		
							if((electronsource==kPi0NoFeedDown)){
								fPt_elec_pi0_MB2->Fill(fPt);  
								if(fNonHFE->IsULS()) fPt_elec_pi0_MB_ULS2->Fill(fPt,fNonHFE->GetNULS());
								if(fNonHFE->IsLS()) fPt_elec_pi0_MB_LS2->Fill(fPt,fNonHFE->GetNLS());
							
								weight_pi0=0;
								if(ptmotherw >= 0.5 && ptmotherw < 2) weight_pi0 = fPi0w->Eval(ptmotherw);
								if(ptmotherw >= 2 && ptmotherw < 6) weight_pi0 = fPi0w2->Eval(ptmotherw);
								if(ptmotherw >= 6 && ptmotherw < 20) weight_pi0 = fPi0w3->Eval(ptmotherw);
								
								fPt_elec_from_pi02->Fill(fPt,weight_pi0);
								if(fNonHFE->IsULS()) fPtElec_ULS_MC_from_pi02->Fill(fPt,fNonHFE->GetNULS()*weight_pi0);
								if(fNonHFE->IsLS()) fPtElec_LS_MC_from_pi02->Fill(fPt,fNonHFE->GetNLS()*weight_pi0);
								
								fPt_elec_phot2->Fill(fPt,weight_pi0);  
								if(fNonHFE->IsULS()) fPtElec_ULS_MC2->Fill(fPt,fNonHFE->GetNULS()*weight_pi0);
								if(fNonHFE->IsLS()) fPtElec_LS_MC2->Fill(fPt,fNonHFE->GetNLS()*weight_pi0);
							
							}	
							
								
							if((electronsource==kEtaNoFeedDown)) {
								fPt_elec_eta_MB2->Fill(fPt);  
								if(fNonHFE->IsULS()) fPt_elec_eta_MB_ULS2->Fill(fPt,fNonHFE->GetNULS());
								if(fNonHFE->IsLS()) fPt_elec_eta_MB_LS2->Fill(fPt,fNonHFE->GetNLS());
								
								weight_eta=0;
								if(ptmotherw >= 0.5 && ptmotherw < 2) weight_eta = fEtaw->Eval(ptmotherw);
								if(ptmotherw >= 2 && ptmotherw < 6) weight_eta = fEtaw2->Eval(ptmotherw);
								if(ptmotherw >= 6 && ptmotherw < 20) weight_eta = fEtaw3->Eval(ptmotherw);
								//cout<<"ptmotherw = "<<ptmotherw<<" weight_eta = "<<weight_eta<<endl;
								fPt_elec_from_eta2->Fill(fPt,weight_eta);
								if(fNonHFE->IsULS()) fPtElec_ULS_MC_from_eta2->Fill(fPt,fNonHFE->GetNULS()*weight_eta);
								if(fNonHFE->IsLS()) fPtElec_LS_MC_from_eta2->Fill(fPt,fNonHFE->GetNLS()*weight_eta);
								
								fPt_elec_phot2->Fill(fPt,weight_eta);  
								if(fNonHFE->IsULS()) fPtElec_ULS_MC2->Fill(fPt,fNonHFE->GetNULS()*weight_eta);
								if(fNonHFE->IsLS()) fPtElec_LS_MC2->Fill(fPt,fNonHFE->GetNLS()*weight_eta);
								
							}		
							
							
							if((electronsource==kGPi0NoFeedDown)) {
								fPt_elec_gamma_MB2->Fill(fPt);  
								if(fNonHFE->IsULS()) fPt_elec_gamma_MB_ULS2->Fill(fPt,fNonHFE->GetNULS());
								if(fNonHFE->IsLS()) fPt_elec_gamma_MB_LS2->Fill(fPt,fNonHFE->GetNLS());
								
								weight_pi0=0;
								if(ptmotherw >= 0.5 && ptmotherw < 2) weight_pi0 = fPi0w->Eval(ptmotherw);
								if(ptmotherw >= 2 && ptmotherw < 6) weight_pi0 = fPi0w2->Eval(ptmotherw);
								if(ptmotherw >= 6 && ptmotherw < 20) weight_pi0 = fPi0w3->Eval(ptmotherw);
								fPt_elec_from_gamma2->Fill(fPt,weight_pi0);
								if(fNonHFE->IsULS()) fPtElec_ULS_MC_from_gamma2->Fill(fPt,fNonHFE->GetNULS()*weight_pi0);
								if(fNonHFE->IsLS()) fPtElec_LS_MC_from_gamma2->Fill(fPt,fNonHFE->GetNLS()*weight_pi0);
								
								fPt_elec_phot2->Fill(fPt,weight_pi0);  
								if(fNonHFE->IsULS()) fPtElec_ULS_MC2->Fill(fPt,fNonHFE->GetNULS()*weight_pi0);
								if(fNonHFE->IsLS()) fPtElec_LS_MC2->Fill(fPt,fNonHFE->GetNLS()*weight_pi0);
								
								
							}			
							
							if((electronsource==kGEtaNoFeedDown)) {
								fPt_elec_gamma_MB2->Fill(fPt);  
								if(fNonHFE->IsULS()) fPt_elec_gamma_MB_ULS2->Fill(fPt,fNonHFE->GetNULS());
								if(fNonHFE->IsLS()) fPt_elec_gamma_MB_LS2->Fill(fPt,fNonHFE->GetNLS());
								
								weight_eta=0;
								if(ptmotherw >= 0.5 && ptmotherw < 2) weight_eta = fEtaw->Eval(ptmotherw);
								if(ptmotherw >= 2 && ptmotherw < 6) weight_eta = fEtaw2->Eval(ptmotherw);
								if(ptmotherw >= 6 && ptmotherw < 20) weight_eta = fEtaw3->Eval(ptmotherw);
								fPt_elec_from_gamma2->Fill(fPt,weight_eta);
								if(fNonHFE->IsULS()) fPtElec_ULS_MC_from_gamma2->Fill(fPt,fNonHFE->GetNULS()*weight_eta);
								if(fNonHFE->IsLS()) fPtElec_LS_MC_from_gamma2->Fill(fPt,fNonHFE->GetNLS()*weight_eta);
								
								fPt_elec_phot2->Fill(fPt,weight_eta);  
								if(fNonHFE->IsULS()) fPtElec_ULS_MC2->Fill(fPt,fNonHFE->GetNULS()*weight_eta);
								if(fNonHFE->IsLS()) fPtElec_LS_MC2->Fill(fPt,fNonHFE->GetNLS()*weight_eta);
								
							}
							
														
							      
					} ///End of photonic electrons
		                

			}
            ///---------------------------------------------------------

		        
		
     }//End of track loop
	
//=======================================================================
	
	delete fListOfmotherkink;
	PostData(1, fOutputList);
}      

//=======================================================================
void AliAnalysisHFETPCTOF_2::Terminate(Option_t *) 
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

//=======================================================================
Bool_t AliAnalysisHFETPCTOF_2::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
//Check single track cuts for a given cut step
//Note this function is called inside the UserExec function
	const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
	if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
	return kTRUE;
}

//Setter for the Partner cuts
void AliAnalysisHFETPCTOF_2::SetPartnerCuts(Float_t Mass, Float_t MinPt, Float_t TpcNclus) {
	fMass = Mass; 
	fMinPt = MinPt; 
	fTpcNclusAsso = TpcNclus;
}

//Setter for the PID cuts (TOF and TPC)
void AliAnalysisHFETPCTOF_2::SetPIDCuts(Float_t tpcPIDmincut, Float_t tpcPIDmaxcut, Float_t tofPIDmincut, Float_t tofPIDmaxcut) {
	ftpcPIDmincut = tpcPIDmincut; 
	ftpcPIDmaxcut = tpcPIDmaxcut;
	ftofPIDmincut = tofPIDmincut;
	ftofPIDmaxcut = tofPIDmaxcut;
}
//Setter for the Eta cut
void AliAnalysisHFETPCTOF_2::SetEtaCut(Float_t EtaMin, Float_t EtaMax){
	fEtaMin = EtaMin; 
	fEtaMax = EtaMax;
}


//=======================================================================
Bool_t AliAnalysisHFETPCTOF_2::FindMother(Int_t mcIndex)
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
//======================================================================

Bool_t AliAnalysisHFETPCTOF_2::IsHFelectronsMC(AliVTrack *track){

          	            
			fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
			float pdg = fMCparticle->GetPdgCode();
			///Is electron:
			if(TMath::Abs(pdg) == 11){
				Bool_t MotherFound = FindMother(TMath::Abs(track->GetLabel()));
				if(MotherFound){
					if(fIsHFE1){
						return kTRUE;
					}
					else{
						return kFALSE;
					}
				}
				else{
					return kFALSE;
				}		
			}
			else{
				return kFALSE;
			}
}

//_________________________________________ for taggign efficiency weights ________________________________________
Int_t AliAnalysisHFETPCTOF_2::GetPi0EtaType(AliAODMCParticle *pi0eta, TClonesArray *mcArray){
    //
    // Return what type of pi0, eta it is
    //
   
    // IsPrimary
    Bool_t primMC = pi0eta->IsPrimary();
    if(!primMC) return kNoIsPrimary;
    
    // Mother
    Int_t motherlabel = pi0eta->GetMother();
    if(motherlabel<0) return kNoMother;
    else {
    
        AliAODMCParticle *mother = (AliAODMCParticle*)fMCarray->At(motherlabel);
        Int_t motherpdg = TMath::Abs(mother->GetPdgCode());
        if(motherpdg == 111 || motherpdg == 221 || motherpdg == 223 || motherpdg == 333 || motherpdg == 331 || motherpdg == 113 || motherpdg == 213 || motherpdg == 313 || motherpdg == 323) return kLightMesons;
        
        if ( (int(TMath::Abs(motherpdg)/100.)%10) == 5 || (int(TMath::Abs(motherpdg)/1000.)%10) == 5 ) return kBeauty;
        if ( (int(TMath::Abs(motherpdg)/100.)%10) == 4 || (int(TMath::Abs(motherpdg)/1000.)%10) == 4 ) return kCharm;
        return kNoFeedDown;    
        
    }
}


//___________________________________________________________________________________________________________

Int_t AliAnalysisHFETPCTOF_2::GetElecSourceType(AliAODMCParticle *electron, TClonesArray *mcArray,Double_t &ptm){
    //
    // Return what type of gammax it is
    //
    
    
    
    // Mother
    Int_t motherlabel = electron->GetMother();
    if(motherlabel<0) return kNoMotherE;
    else {
        
        AliAODMCParticle *mother = (AliAODMCParticle*)mcArray->At(motherlabel);
        Int_t motherpdg = TMath::Abs(mother->GetPdgCode());
        ptm=mother->Pt();
        if(motherpdg == 111){
            Int_t typepi0eta = GetPi0EtaType(mother,mcArray);
            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) 	return kPi0NoFeedDown;
        }
        if(motherpdg == 221){
            Int_t typepi0eta = GetPi0EtaType(mother,mcArray);
            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kEtaNoFeedDown;
        }
        if(motherpdg == 22) {
            Int_t gmotherlabel = mother->GetMother();
            if(gmotherlabel<0) return kDirectGamma;
            else {
                AliAODMCParticle *gmother = (AliAODMCParticle*)mcArray->At(gmotherlabel);
                ptm=gmother->Pt();
                Int_t gmotherpdg = TMath::Abs(gmother->GetPdgCode());
                if(gmotherpdg == 111) {
                    Int_t typepi0eta = GetPi0EtaType(gmother,mcArray);
                    if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGPi0NoFeedDown;
                }
                if(gmotherpdg == 221) {
                    Int_t typepi0eta = GetPi0EtaType(gmother,mcArray);
                    if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGEtaNoFeedDown;
                }
                if(gmotherpdg == 22) {
                    Int_t ggmotherlabel = gmother->GetMother();
                    if(ggmotherlabel<0) return kDirectGamma;
                    else {
                        AliAODMCParticle *ggmother = (AliAODMCParticle*)mcArray->At(ggmotherlabel);
                        ptm=ggmother->Pt();
                        Int_t ggmotherpdg = TMath::Abs(ggmother->GetPdgCode());
                        if(ggmotherpdg == 111) {
                            Int_t typepi0eta = GetPi0EtaType(ggmother,mcArray);
                            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGPi0NoFeedDown;
                        }
                        if(ggmotherpdg == 221) {
                            Int_t typepi0eta = GetPi0EtaType(ggmother,mcArray);
                            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGEtaNoFeedDown;
                        }
                    }
                }
            }
        }
    }
    
    return kOthersE;
    
}





