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
//      Task for Beauty analysis in Pb-Pb collisions   				  //
//      															  //
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
#include "AliAnalysisHFETPCTOFBeauty.h"
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
#include "AliGenEventHeader.h"
#include "TRandom3.h"
//______________________________________________________________________

//______________________________________________________________________
ClassImp(AliAnalysisHFETPCTOFBeauty)

//______________________________________________________________________
AliAnalysisHFETPCTOFBeauty::AliAnalysisHFETPCTOFBeauty(const char *name)
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
,ftpcPIDmincut(0)
,ftpcPIDmaxcut(0)
,ftofPIDmincut(0)
,ftofPIDmaxcut(0)
,fEtaMin(0)
,fEtaMax(0)

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
,fTPCnsigma_p_after_tof_p(0)
,fTPCnsigma_p_after_tof_pion(0)
,fTPCnsigma_p_after_tof_k(0)
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
,fTPCnsigma_TOFnsigma1(0)
,fTPCnsigma_TOFnsigma2(0)
,fTPCnsigma_TOFnsigma3(0)
,fDCAxy_pt_had(0)
,fDCAz_pt_had(0)
,fDCAxy_pt_ele(0)
,fDCAz_pt_ele(0)
,hCharmMotherPt(0)
,hCharmMotherPt_vsElecPt(0)
,hElecPt_vsCharmMotherPt(0)
,hCharmMotherPt_corr(0)
,hCharmMotherPt_corr2(0)
,hCharmMotherPt_corr3(0)
,hCharmMotherPt_corr4(0)
,hCharmMotherPt1_corr(0)
,hCharmMotherPt2_corr(0)
,hCharmMotherPt3_corr(0)
,hCharmMotherPt4_corr(0)
,hCharmMotherPt5_corr(0)
,hCharmMotherPt6_corr(0)
,hCharmMotherPt7_corr(0)
,hCharmMotherPt8_corr(0)
,hCharmMotherPt9_corr(0)
,hCharmMotherPt10_corr(0)
,hCharmMotherPt11_corr(0)
,hCharmMotherPt12_corr(0)
,hCharmMotherPt13_corr(0)
,hCharmMotherPt14_corr(0)
,hCharmMotherPt15_corr(0)
,hBeautyMotherPt(0)
,fPtBeautyGenerated(0)
,fPtBeautyReconstructedTracks(0)
,fPtBeautyReconstructedTracksPID(0)

//For the HFE package
,fCuts(0)
,fCFM(0)
,fPID(new AliHFEpid("hfePid"))
,fPIDqa(0)

//For MC
,fMCstack(0)
,fRejectKinkMother(kTRUE)
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
,fD0(0)
,fD0Data(0)
,fNpureMC(0)
,fNembMCpi0(0)
,fNembMCeta(0)
,fNTotMCpart(0)

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
AliAnalysisHFETPCTOFBeauty::AliAnalysisHFETPCTOFBeauty()
: AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisHFETPCTOFBeauty")

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
,ftpcPIDmincut(0)
,ftpcPIDmaxcut(0)
,ftofPIDmincut(0)
,ftofPIDmaxcut(0)
,fEtaMin(0)
,fEtaMax(0)

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
,fTPCnsigma_p_after_tof_p(0)
,fTPCnsigma_p_after_tof_pion(0)
,fTPCnsigma_p_after_tof_k(0)
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
,fTPCnsigma_TOFnsigma1(0)
,fTPCnsigma_TOFnsigma2(0)
,fTPCnsigma_TOFnsigma3(0)
,fDCAxy_pt_had(0)
,fDCAz_pt_had(0)
,fDCAxy_pt_ele(0)
,fDCAz_pt_ele(0)
,hCharmMotherPt(0)
,hCharmMotherPt_vsElecPt(0)
,hElecPt_vsCharmMotherPt(0)
,hCharmMotherPt_corr(0)
,hCharmMotherPt_corr2(0)
,hCharmMotherPt_corr3(0)
,hCharmMotherPt_corr4(0)
,hCharmMotherPt1_corr(0)
,hCharmMotherPt2_corr(0)
,hCharmMotherPt3_corr(0)
,hCharmMotherPt4_corr(0)
,hCharmMotherPt5_corr(0)
,hCharmMotherPt6_corr(0)
,hCharmMotherPt7_corr(0)
,hCharmMotherPt8_corr(0)
,hCharmMotherPt9_corr(0)
,hCharmMotherPt10_corr(0)
,hCharmMotherPt11_corr(0)
,hCharmMotherPt12_corr(0)
,hCharmMotherPt13_corr(0)
,hCharmMotherPt14_corr(0)
,hCharmMotherPt15_corr(0)
,hBeautyMotherPt(0)
,fPtBeautyGenerated(0)
,fPtBeautyReconstructedTracks(0)
,fPtBeautyReconstructedTracksPID(0)

//For the HFE package
,fCuts(0)
,fCFM(0)
,fPID(new AliHFEpid("hfePid"))
,fPIDqa(0)

//For MC
,fMCstack(0)
,fRejectKinkMother(kTRUE)
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
,fD0(0)
,fD0Data(0)
,fNpureMC(0)
,fNembMCpi0(0)
,fNembMCeta(0)
,fNTotMCpart(0)


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
AliAnalysisHFETPCTOFBeauty::~AliAnalysisHFETPCTOFBeauty()
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
void AliAnalysisHFETPCTOFBeauty::UserCreateOutputObjects()
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
    fCent = new TH1F("fCent","Centrality",100,0,100);
    fCent2 = new TH1F("fCent2","Centrality",100,0,100);
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
    
    
    Double_t ptbinning[33] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10};
    
    Double_t ptbinningHF[14] = {1,2,3,4,5,6,7,8,10,12,16,24,36,50};
    
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
    
    hCharmMotherPt = new TH1F("hCharmMotherPt","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt);
    
    hCharmMotherPt_vsElecPt = new TH2F("hCharmMotherPt_vsElecPt","; p_{T} [GeV/c]; Count",100,0,20,100,0,20);
    fOutputList->Add(hCharmMotherPt_vsElecPt);
    
    hElecPt_vsCharmMotherPt = new TH2F("hElecPt_vsCharmMotherPt","; p_{T} [GeV/c]; Count",100,0,20,100,0,20);
    fOutputList->Add(hElecPt_vsCharmMotherPt);
    
    hCharmMotherPt_corr = new TH1F("hCharmMotherPt_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt_corr);
    
    hCharmMotherPt_corr2 = new TH1F("hCharmMotherPt_corr2","; p_{T} [GeV/c]; Count",100,0,50);
    fOutputList->Add(hCharmMotherPt_corr2);

	hCharmMotherPt_corr3 = new TH1F("hCharmMotherPt_corr3","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt_corr3);
    
    hCharmMotherPt_corr4 = new TH1F("hCharmMotherPt_corr4","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt_corr4);
    
    hCharmMotherPt1_corr = new TH1F("hCharmMotherPt1_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt1_corr);
    
	hCharmMotherPt2_corr = new TH1F("hCharmMotherPt2_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt2_corr);
    
	hCharmMotherPt3_corr = new TH1F("hCharmMotherPt3_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt3_corr);
    
	hCharmMotherPt4_corr = new TH1F("hCharmMotherPt4_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt4_corr);
    
	hCharmMotherPt5_corr = new TH1F("hCharmMotherPt5_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt5_corr);
    
	hCharmMotherPt6_corr = new TH1F("hCharmMotherPt6_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt6_corr);
    
	hCharmMotherPt7_corr = new TH1F("hCharmMotherPt7_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt7_corr);
    
	hCharmMotherPt8_corr = new TH1F("hCharmMotherPt8_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt8_corr); 
    
	hCharmMotherPt9_corr = new TH1F("hCharmMotherPt9_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt9_corr); 
    
	hCharmMotherPt10_corr = new TH1F("hCharmMotherPt10_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt10_corr);
    
	hCharmMotherPt11_corr = new TH1F("hCharmMotherPt11_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt11_corr);
    
	hCharmMotherPt12_corr = new TH1F("hCharmMotherPt12_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt12_corr);
    
	hCharmMotherPt13_corr = new TH1F("hCharmMotherPt13_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt13_corr);
    
	hCharmMotherPt14_corr = new TH1F("hCharmMotherPt14_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt14_corr);
    
	hCharmMotherPt15_corr = new TH1F("hCharmMotherPt15_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt15_corr);

	hBeautyMotherPt = new TH1F("hBeautyMotherPt","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hBeautyMotherPt);

    fPtElec = new TH1F("fPtElec","; p_{T} [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPtElec);
    
    hPtD0 = new TH1F("hPtD0","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hPtD0);
    
    hPtLambdaC = new TH1F("hPtLambdaC","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hPtLambdaC);
     
    fPElec = new TH1F("fPElec","; p [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPElec);

    fPtHad_f = new TH1F("fPtHad_f","; p_{T} [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPtHad_f);
    fPtHad_f->Sumw2();
    
    fPHad_f = new TH1F("fPHad_f","; p [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPHad_f);
    fPHad_f->Sumw2();
       
    fDCAxy_pt_had = new TH2F("fDCAxy_pt_had",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_pt_had);
    
    fDCAz_pt_had = new TH2F("fDCAz_pt_had",";p_{t} (GeV/c);DCAz hadrons",300,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAz_pt_had);
    
    fDCAxy_pt_ele = new TH2F("fDCAxy_pt_ele",";p_{t} (GeV/c);DCAxy ",300,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_pt_ele);
    
    fDCAz_pt_ele = new TH2F("fDCAz_pt_ele",";p_{t} (GeV/c);DCAz ",300,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAz_pt_ele); 
    
    fPtMCeta = new TH1F("fPtMCeta",";p_{t} (GeV/c)",2000,0,100);
    fOutputList->Add(fPtMCeta);
    
    fTPCnsigma_p_after_tof = new TH2F("fTPCnsigma_p_after_tof","p (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_p_after_tof);
    
    fTPCnsigma_p_after_tof_p = new TH2F("fTPCnsigma_p_after_tof_p","p (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_p_after_tof_p);
    
    fTPCnsigma_p_after_tof_pion = new TH2F("fTPCnsigma_p_after_tof_pion","p (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_p_after_tof_pion);
    
    fTPCnsigma_p_after_tof_k = new TH2F("fTPCnsigma_p_after_tof_k","p (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_p_after_tof_k);
       
    fTPCnsigma_p_after_tof_its = new TH2F("fTPCnsigma_p_after_tof_its","p (GeV/c);TPC Electron N#sigma after TOF and ITS cuts",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_p_after_tof_its);
    
    
    fTPCnsigma_pt_after_tof = new TH2F("fTPCnsigma_pt_after_tof","pt (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_pt_after_tof);
    
    fTPCnsigma_pt_after_tof_its = new TH2F("fTPCnsigma_pt_after_tof_its","pt (GeV/c);TPC Electron N#sigma after TOF and ITS cuts",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_pt_after_tof_its);
    
    fPtBeautyGenerated = new TH1F("fPtBeautyGenerated","; p_{T} [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPtBeautyGenerated);
    
    fPtBeautyReconstructedTracks = new TH1F("fPtBeautyReconstructedTracks","; p_{T} [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPtBeautyReconstructedTracks);
    
    fPtBeautyReconstructedTracksPID = new TH1F("fPtBeautyReconstructedTracksPID","; p_{T} [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPtBeautyReconstructedTracksPID);
      
    ///THnSparse to store DCA of different particle species-------------
    
    Int_t nBinspdg2 = 12;
    Double_t minpdg2 = 0.;
    Double_t maxpdg2 = 12.;
    Double_t binLimpdg2[nBinspdg2+1];
    for(Int_t i=0; i<=nBinspdg2; i++) binLimpdg2[i]=(Double_t)minpdg2 + (maxpdg2-minpdg2)/nBinspdg2*(Double_t)i ;
    
    Int_t nBinsdcaxy = 800;
    Double_t mindcaxy = -0.2;
    Double_t maxdcaxy = 0.2;
    Double_t binLimdcaxy[nBinsdcaxy+1];
    for(Int_t i=0; i<=nBinsdcaxy; i++) binLimdcaxy[i]=(Double_t)mindcaxy + (maxdcaxy-mindcaxy)/nBinsdcaxy*(Double_t)i ;
    
    Int_t nBinsg = 4;
    Double_t ming = 0.;
    Double_t maxg = 4.;
    Double_t binLimg[nBinsg+1];
    for(Int_t i=0; i<=nBinsg; i++) binLimg[i]=(Double_t)ming + (maxg-ming)/nBinsg*(Double_t)i ;
    
    Int_t nBinsR = 1800;
    Double_t minR = 0;
    Double_t maxR = 60;
    Double_t binLimR[nBinsR+1];
    for(Int_t i=0; i<=nBinsR; i++) binLimR[i]=(Double_t)minR + (maxR-minR)/nBinsR*(Double_t)i ;
    
    Int_t nBinsITSchi2 = 400;
    Double_t minITSchi2 = 0.;
    Double_t maxITSchi2 = 100;
    Double_t binLimITSchi2[nBinsITSchi2+1];
    for(Int_t i=0; i<=nBinsITSchi2; i++) binLimITSchi2[i]=(Double_t)minITSchi2 + (maxITSchi2-minITSchi2)/nBinsITSchi2*(Double_t)i ;
    
    Int_t nBinsITSsha = 100;
    Double_t minITSsha = 0.;
    Double_t maxITSsha = 1;
    Double_t binLimITSsha[nBinsITSsha+1];
    for(Int_t i=0; i<=nBinsITSsha; i++) binLimITSsha[i]=(Double_t)minITSsha + (maxITSsha-minITSsha)/nBinsITSsha*(Double_t)i ;
    
    Int_t nBinstype = 7;
    Double_t mintype = -1.;
    Double_t maxtype = 6.;
    Double_t binLimtype[nBinstype+1];
    for(Int_t i=0; i<=nBinstype; i++) binLimtype[i]=(Double_t)mintype + (maxtype-mintype)/nBinstype*(Double_t)i ;
    
    /*
     Double_t binLimpt2[45]= {0.1,0.112797,0.127231,0.143512,0.161877,0.182592,0.205957,0.232313,0.262041,0.295573,0.333397,
     0.37606,0.424183,0.478465,0.539692,0.608754,0.686654,0.774523,0.873636,0.985432,1.11153,1.25377,1.41421,1.59519,1.79932,2.02957,
     2.28928,2.58223,2.91267,3.2854,3.70582,4.18004,4.71494,5.3183,5.99886,6.76651,7.6324,8.60909,9.71076,10.9534,12.3551,13.9361,15.7195,17.731,20};//bin limits from the measured pi0 spectrum
     */
    
    const Int_t nDima2=8;
    Int_t nBina2[nDima2] = {32,nBinspdg2,nBinsdcaxy,nBinsg,nBinsR,nBinsITSchi2,nBinsITSsha,nBinstype};
    fD0 = new THnSparseF("fD0","fD0",nDima2,nBina2);
    fD0->SetBinEdges(0,ptbinning); ///pt spectra -> same binning as other histograms
    fD0->SetBinEdges(1,binLimpdg2); /// electrons from D,charm baryons, B, beauty baryons, gamma, pi0, eta, Dcorrected, Dcorrected by weight, protons, kaons
    fD0->SetBinEdges(2,binLimdcaxy); ///dca distribution
    fD0->SetBinEdges(3,binLimg);  ///From which generator (Hijing, else, pi0, eta)
    fD0->SetBinEdges(4,binLimR); ///Position where the electron is created
    fD0->SetBinEdges(5,binLimITSchi2); ///ITS chi2 
    fD0->SetBinEdges(6,binLimITSsha); ///fraction ITS shared clusters 
    fD0->SetBinEdges(7,binLimtype); ///pi0 and eta type  ///kNoMother, kNoFeedDown, kNoIsPrimary, kLightMesons, kBeauty, kCharm
    fD0->Sumw2();
    fOutputList->Add(fD0);
    
    ///-----------------------------------------------------------------
    
    
    
    ///THnSparse to store DCA in Data
    const Int_t nDima3=5;
    Int_t nBina3[nDima3] = {32,nBinsdcaxy,nBinsITSchi2,nBinsITSsha,nBinspdg2};
    fD0Data = new THnSparseF("fD0Data","fD0Data",nDima3,nBina3);
    fD0Data->SetBinEdges(0,ptbinning); ///pt spectra -> same binning as other histograms
    fD0Data->SetBinEdges(1,binLimdcaxy); ///dca distribution
    fD0Data->SetBinEdges(2,binLimITSchi2); ///ITS chi2 
    fD0Data->SetBinEdges(3,binLimITSsha); ///fraction ITS shared clusters 
    fD0Data->SetBinEdges(4,binLimpdg2); /// electrons and pions
    fD0Data->Sumw2();
    fOutputList->Add(fD0Data);
    
    
    PostData(1, fOutputList);
    
    ///______________________________________________________________________
}

	


//______________________________________________________________________
//Main loop
//Called for each event
void AliAnalysisHFETPCTOFBeauty::UserExec(Option_t *)
{
	gRandom->SetSeed(0);
	//gRandom = new TRandom3(0);
    Int_t pdg = -99999;
    Int_t pdg_mother = -99999;
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
    Int_t fITSnClus = 99999;
    Int_t fTPCnClus = 99999;
    Double_t qadca[7];
    Double_t qadcaData[7];
    
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
 
 	///Checking Nv0 and NT0 before event and vertex selection--------------
    fNevent->Fill(0);///triggered by v0
        
    ///Number of events in kint7 that is also triggered by T0
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
    //if(!vtxTtl.Contains("VertexerTracks") || !(fAOD->GetPrimaryVertexSPD())) fNevent_no_vertex_2->Fill(0);
    ///-------------------------------------------------------------
        
    ////////////////////
	//Correlation cuts///
	////////////////////
	if(!PassCorrCuts(fAOD)) return;
        
    ///Number of events rejected by the correlation cuts between SPD and track vertexes
    fNevent_corrcut->Fill(0);
        
    ////////////////////
	//Track vertex cut//
	////////////////////
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

       
     ///Checking Nv0 and NT0 after the vertex selection----------------
     fNevent_3->Fill(0); ///Number of analysed events ///triggered by v0
        
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
        
   
    /////////////////////////
	//Centrality Selection///
	/////////////////////////
	if(!fIsPP){  

        Float_t lPercentile = 300;
        AliMultSelection *MultSelection = 0x0;
        MultSelection = (AliMultSelection * ) fVevent->FindListObject("MultSelection");
        if(!MultSelection) {
            //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
            AliWarning("AliMultSelection object not found!");
        }
        else{
            lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
        }
        
        Double_t cent = MultSelection->GetMultiplicityPercentile("V0M", true); // returns centrality only for events used in calibration
        
        //Centrality
        fCent->Fill(cent);
        
        if(cent<0 || cent>10) return;
        //if(cent<30.5 || cent>50.5) return; //Comparing with Denise
        
        fCent2->Fill(cent);
    }
    
   
    ///Sign of the magnetic field
	//cout<<"fAOD->GetMagneticField() = "<<fAOD->GetMagneticField()<<endl;
    int signB = 0;
    if(fAOD->GetMagneticField() < 0) signB = -1;
    if(fAOD->GetMagneticField() > 0) signB = 1;
    
    ///Function from which I get a random number between 0 and 1
    TF1 *fconstant = new TF1("fconstant","x",0,1);
    
    ///Functions used to correct D meson spectrum from MC
    TF1 *fDmesonShape1 = new TF1("fDmesonShape1","(1/[0])*(CrossOverLc(2.75563e+00,2.34650e+00,x[0])*(-4.22294e-01)/TMath::Power((-1.49187e+00)+x[0]/(4.14660e-07),1.71933e-01)+CrossOverRc(2.75563e+00,2.34650e+00,x[0])*1.64461e+00/TMath::Power((-2.58810e-01)+x/2.21576e+00,4.29324e+00)*CrossOverLc(5.67592e+00,(7.62946e-01)+CrossOverRc(5.67592e+00,7.62946e-01,x[0])*(5.49622e-01)/TMath::Power(x[0],(-2.15261e+00)),x[0]))/(CrossOverLc(3.38860e+00,4.05992e-01,x[0])*(-7.12382e-02)/TMath::Power((4.57260e-01)+x[0]/(5.36517e+02),4.07632e+00)+CrossOverRc(3.38860e+00,4.05992e-01,x[0])*(3.47577e-03)/TMath::Power(8.63986e-02+x/(2.06117e+01),2.58547e+00)*CrossOverLc(1.00000e+00,(9.99900e+03)+CrossOverRc(1.00000e+00,9.99900e+03,x[0])*(0.00000e+00)/TMath::Power(x[0],0.00000e+00),x[0]))",2,30);
    fDmesonShape1->SetParameter(0,8.50802);
    
    TF1 *fDmesonShape2 = new TF1("fDmesonShape2","(1/[0])*(CrossOverLc(2.75563e+00,2.34650e+00,x[0])*(-4.22294e-01)/TMath::Power((-1.49187e+00)+x[0]/(4.14660e-07),1.71933e-01)+CrossOverRc(2.75563e+00,2.34650e+00,x[0])*1.64461e+00/TMath::Power((-2.58810e-01)+x/2.21576e+00,4.29324e+00)*CrossOverLc(5.67592e+00,(7.62946e-01)+CrossOverRc(5.67592e+00,7.62946e-01,x[0])*(5.49622e-01)/TMath::Power(x[0],(-2.15261e+00)),x[0]))/(CrossOverLc(3.38860e+00,4.05992e-01,x[0])*(-7.12382e-02)/TMath::Power((4.57260e-01)+x[0]/(5.36517e+02),4.07632e+00)+CrossOverRc(3.38860e+00,4.05992e-01,x[0])*(3.47577e-03)/TMath::Power(8.63986e-02+x/(2.06117e+01),2.58547e+00)*CrossOverLc(1.00000e+00,(9.99900e+03)+CrossOverRc(1.00000e+00,9.99900e+03,x[0])*(0.00000e+00)/TMath::Power(x[0],0.00000e+00),x[0]))",2.25,30);   
    fDmesonShape2->SetParameter(0,8.31432);
    
    TF1 *fDmesonShape3 = new TF1("fDmesonShape3","(1/[0])*(CrossOverLc(2.75563e+00,2.34650e+00,x[0])*(-4.22294e-01)/TMath::Power((-1.49187e+00)+x[0]/(4.14660e-07),1.71933e-01)+CrossOverRc(2.75563e+00,2.34650e+00,x[0])*1.64461e+00/TMath::Power((-2.58810e-01)+x/2.21576e+00,4.29324e+00)*CrossOverLc(5.67592e+00,(7.62946e-01)+CrossOverRc(5.67592e+00,7.62946e-01,x[0])*(5.49622e-01)/TMath::Power(x[0],(-2.15261e+00)),x[0]))/(CrossOverLc(3.38860e+00,4.05992e-01,x[0])*(-7.12382e-02)/TMath::Power((4.57260e-01)+x[0]/(5.36517e+02),4.07632e+00)+CrossOverRc(3.38860e+00,4.05992e-01,x[0])*(3.47577e-03)/TMath::Power(8.63986e-02+x/(2.06117e+01),2.58547e+00)*CrossOverLc(1.00000e+00,(9.99900e+03)+CrossOverRc(1.00000e+00,9.99900e+03,x[0])*(0.00000e+00)/TMath::Power(x[0],0.00000e+00),x[0]))",2.5,30);   
    fDmesonShape3->SetParameter(0,8.01983);
    
    TF1 *fDmesonShape4 = new TF1("fDmesonShape4","(1/[0])*(CrossOverLc(2.75563e+00,2.34650e+00,x[0])*(-4.22294e-01)/TMath::Power((-1.49187e+00)+x[0]/(4.14660e-07),1.71933e-01)+CrossOverRc(2.75563e+00,2.34650e+00,x[0])*1.64461e+00/TMath::Power((-2.58810e-01)+x/2.21576e+00,4.29324e+00)*CrossOverLc(5.67592e+00,(7.62946e-01)+CrossOverRc(5.67592e+00,7.62946e-01,x[0])*(5.49622e-01)/TMath::Power(x[0],(-2.15261e+00)),x[0]))/(CrossOverLc(3.38860e+00,4.05992e-01,x[0])*(-7.12382e-02)/TMath::Power((4.57260e-01)+x[0]/(5.36517e+02),4.07632e+00)+CrossOverRc(3.38860e+00,4.05992e-01,x[0])*(3.47577e-03)/TMath::Power(8.63986e-02+x/(2.06117e+01),2.58547e+00)*CrossOverLc(1.00000e+00,(9.99900e+03)+CrossOverRc(1.00000e+00,9.99900e+03,x[0])*(0.00000e+00)/TMath::Power(x[0],0.00000e+00),x[0]))",2.75,30);   
    fDmesonShape4->SetParameter(0,7.26397);
    
    TF1 *fDmesonShape5 = new TF1("fDmesonShape5","(1/[0])*(CrossOverLc(2.75563e+00,2.34650e+00,x[0])*(-4.22294e-01)/TMath::Power((-1.49187e+00)+x[0]/(4.14660e-07),1.71933e-01)+CrossOverRc(2.75563e+00,2.34650e+00,x[0])*1.64461e+00/TMath::Power((-2.58810e-01)+x/2.21576e+00,4.29324e+00)*CrossOverLc(5.67592e+00,(7.62946e-01)+CrossOverRc(5.67592e+00,7.62946e-01,x[0])*(5.49622e-01)/TMath::Power(x[0],(-2.15261e+00)),x[0]))/(CrossOverLc(3.38860e+00,4.05992e-01,x[0])*(-7.12382e-02)/TMath::Power((4.57260e-01)+x[0]/(5.36517e+02),4.07632e+00)+CrossOverRc(3.38860e+00,4.05992e-01,x[0])*(3.47577e-03)/TMath::Power(8.63986e-02+x/(2.06117e+01),2.58547e+00)*CrossOverLc(1.00000e+00,(9.99900e+03)+CrossOverRc(1.00000e+00,9.99900e+03,x[0])*(0.00000e+00)/TMath::Power(x[0],0.00000e+00),x[0]))",3,30);   
    fDmesonShape5->SetParameter(0,6.11613);
    
    TF1 *fDmesonShape6 = new TF1("fDmesonShape6","(1/[0])*(CrossOverLc(2.75563e+00,2.34650e+00,x[0])*(-4.22294e-01)/TMath::Power((-1.49187e+00)+x[0]/(4.14660e-07),1.71933e-01)+CrossOverRc(2.75563e+00,2.34650e+00,x[0])*1.64461e+00/TMath::Power((-2.58810e-01)+x/2.21576e+00,4.29324e+00)*CrossOverLc(5.67592e+00,(7.62946e-01)+CrossOverRc(5.67592e+00,7.62946e-01,x[0])*(5.49622e-01)/TMath::Power(x[0],(-2.15261e+00)),x[0]))/(CrossOverLc(3.38860e+00,4.05992e-01,x[0])*(-7.12382e-02)/TMath::Power((4.57260e-01)+x[0]/(5.36517e+02),4.07632e+00)+CrossOverRc(3.38860e+00,4.05992e-01,x[0])*(3.47577e-03)/TMath::Power(8.63986e-02+x/(2.06117e+01),2.58547e+00)*CrossOverLc(1.00000e+00,(9.99900e+03)+CrossOverRc(1.00000e+00,9.99900e+03,x[0])*(0.00000e+00)/TMath::Power(x[0],0.00000e+00),x[0]))",3.5,30);   
    fDmesonShape6->SetParameter(0,4.16545);
    
    TF1 *fDmesonShape7 = new TF1("fDmesonShape7","(1/2.72057)*((CrossOverLc(2.75563e+00,2.34650e+00,x[0])*(-4.22294e-01)/TMath::Power((-1.49187e+00)+x[0]/(4.14660e-07),1.71933e-01)+CrossOverRc(2.75563e+00,2.34650e+00,x[0])*1.64461e+00/TMath::Power((-2.58810e-01)+x/2.21576e+00,4.29324e+00)*CrossOverLc(5.67592e+00,(7.62946e-01)+CrossOverRc(5.67592e+00,7.62946e-01,x[0])*(5.49622e-01)/TMath::Power(x[0],(-2.15261e+00)),x[0]))/((CrossOverLc(3.38860e+00,4.05992e-01,x[0])*(-7.12382e-02)/TMath::Power((4.57260e-01)+x[0]/(5.36517e+02),4.07632e+00)+CrossOverRc(3.38860e+00,4.05992e-01,x[0])*(3.47577e-03)/TMath::Power(8.63986e-02+x/(2.06117e+01),2.58547e+00)*CrossOverLc(1.00000e+00,(9.99900e+03)+CrossOverRc(1.00000e+00,9.99900e+03,x[0])*(0.00000e+00)/TMath::Power(x[0],0.00000e+00),x[0]))))",4,30);   
    //TF1 *fDmesons_bin7 = new TF1("fDmesons_bin7","10*x[0]",8,10);   
    //fDmesons_bin7->FixParameter(0,0.382748);
    
    TF1 *fDmesonShape8 = new TF1("fDmesonShape8","(1/[0])*(CrossOverLc(2.75563e+00,2.34650e+00,x[0])*(-4.22294e-01)/TMath::Power((-1.49187e+00)+x[0]/(4.14660e-07),1.71933e-01)+CrossOverRc(2.75563e+00,2.34650e+00,x[0])*1.64461e+00/TMath::Power((-2.58810e-01)+x/2.21576e+00,4.29324e+00)*CrossOverLc(5.67592e+00,(7.62946e-01)+CrossOverRc(5.67592e+00,7.62946e-01,x[0])*(5.49622e-01)/TMath::Power(x[0],(-2.15261e+00)),x[0]))/(CrossOverLc(3.38860e+00,4.05992e-01,x[0])*(-7.12382e-02)/TMath::Power((4.57260e-01)+x[0]/(5.36517e+02),4.07632e+00)+CrossOverRc(3.38860e+00,4.05992e-01,x[0])*(3.47577e-03)/TMath::Power(8.63986e-02+x/(2.06117e+01),2.58547e+00)*CrossOverLc(1.00000e+00,(9.99900e+03)+CrossOverRc(1.00000e+00,9.99900e+03,x[0])*(0.00000e+00)/TMath::Power(x[0],0.00000e+00),x[0]))",4.5,30);   
    fDmesonShape8->SetParameter(0,1.87864);
    
    TF1 *fDmesonShape9 = new TF1("fDmesonShape9","(1/[0])*(CrossOverLc(2.75563e+00,2.34650e+00,x[0])*(-4.22294e-01)/TMath::Power((-1.49187e+00)+x[0]/(4.14660e-07),1.71933e-01)+CrossOverRc(2.75563e+00,2.34650e+00,x[0])*1.64461e+00/TMath::Power((-2.58810e-01)+x/2.21576e+00,4.29324e+00)*CrossOverLc(5.67592e+00,(7.62946e-01)+CrossOverRc(5.67592e+00,7.62946e-01,x[0])*(5.49622e-01)/TMath::Power(x[0],(-2.15261e+00)),x[0]))/(CrossOverLc(3.38860e+00,4.05992e-01,x[0])*(-7.12382e-02)/TMath::Power((4.57260e-01)+x[0]/(5.36517e+02),4.07632e+00)+CrossOverRc(3.38860e+00,4.05992e-01,x[0])*(3.47577e-03)/TMath::Power(8.63986e-02+x/(2.06117e+01),2.58547e+00)*CrossOverLc(1.00000e+00,(9.99900e+03)+CrossOverRc(1.00000e+00,9.99900e+03,x[0])*(0.00000e+00)/TMath::Power(x[0],0.00000e+00),x[0]))",5,30);   
    fDmesonShape9->SetParameter(0,1.36889);
    
    TF1 *fDmesonShape10 = new TF1("fDmesonShape10","(1/[0])*(CrossOverLc(2.75563e+00,2.34650e+00,x[0])*(-4.22294e-01)/TMath::Power((-1.49187e+00)+x[0]/(4.14660e-07),1.71933e-01)+CrossOverRc(2.75563e+00,2.34650e+00,x[0])*1.64461e+00/TMath::Power((-2.58810e-01)+x/2.21576e+00,4.29324e+00)*CrossOverLc(5.67592e+00,(7.62946e-01)+CrossOverRc(5.67592e+00,7.62946e-01,x[0])*(5.49622e-01)/TMath::Power(x[0],(-2.15261e+00)),x[0]))/(CrossOverLc(3.38860e+00,4.05992e-01,x[0])*(-7.12382e-02)/TMath::Power((4.57260e-01)+x[0]/(5.36517e+02),4.07632e+00)+CrossOverRc(3.38860e+00,4.05992e-01,x[0])*(3.47577e-03)/TMath::Power(8.63986e-02+x/(2.06117e+01),2.58547e+00)*CrossOverLc(1.00000e+00,(9.99900e+03)+CrossOverRc(1.00000e+00,9.99900e+03,x[0])*(0.00000e+00)/TMath::Power(x[0],0.00000e+00),x[0]))",5.5,30);   
    fDmesonShape10->SetParameter(0,1.03924);
    
    TF1 *fDmesonShape11 = new TF1("fDmesonShape11","(1/[0])*(CrossOverLc(2.75563e+00,2.34650e+00,x[0])*(-4.22294e-01)/TMath::Power((-1.49187e+00)+x[0]/(4.14660e-07),1.71933e-01)+CrossOverRc(2.75563e+00,2.34650e+00,x[0])*1.64461e+00/TMath::Power((-2.58810e-01)+x/2.21576e+00,4.29324e+00)*CrossOverLc(5.67592e+00,(7.62946e-01)+CrossOverRc(5.67592e+00,7.62946e-01,x[0])*(5.49622e-01)/TMath::Power(x[0],(-2.15261e+00)),x[0]))/(CrossOverLc(3.38860e+00,4.05992e-01,x[0])*(-7.12382e-02)/TMath::Power((4.57260e-01)+x[0]/(5.36517e+02),4.07632e+00)+CrossOverRc(3.38860e+00,4.05992e-01,x[0])*(3.47577e-03)/TMath::Power(8.63986e-02+x/(2.06117e+01),2.58547e+00)*CrossOverLc(1.00000e+00,(9.99900e+03)+CrossOverRc(1.00000e+00,9.99900e+03,x[0])*(0.00000e+00)/TMath::Power(x[0],0.00000e+00),x[0]))",6,30);   
    fDmesonShape11->SetParameter(0,0.814814);
    
    TF1 *fDmesonShape12 = new TF1("fDmesonShape12","(1/[0])*(CrossOverLc(2.75563e+00,2.34650e+00,x[0])*(-4.22294e-01)/TMath::Power((-1.49187e+00)+x[0]/(4.14660e-07),1.71933e-01)+CrossOverRc(2.75563e+00,2.34650e+00,x[0])*1.64461e+00/TMath::Power((-2.58810e-01)+x/2.21576e+00,4.29324e+00)*CrossOverLc(5.67592e+00,(7.62946e-01)+CrossOverRc(5.67592e+00,7.62946e-01,x[0])*(5.49622e-01)/TMath::Power(x[0],(-2.15261e+00)),x[0]))/(CrossOverLc(3.38860e+00,4.05992e-01,x[0])*(-7.12382e-02)/TMath::Power((4.57260e-01)+x[0]/(5.36517e+02),4.07632e+00)+CrossOverRc(3.38860e+00,4.05992e-01,x[0])*(3.47577e-03)/TMath::Power(8.63986e-02+x/(2.06117e+01),2.58547e+00)*CrossOverLc(1.00000e+00,(9.99900e+03)+CrossOverRc(1.00000e+00,9.99900e+03,x[0])*(0.00000e+00)/TMath::Power(x[0],0.00000e+00),x[0]))",6.5,30);   
    fDmesonShape12->SetParameter(0,0.655621);
    
    TF1 *fDmesonShape13 = new TF1("fDmesonShape13","(1/[0])*(CrossOverLc(2.75563e+00,2.34650e+00,x[0])*(-4.22294e-01)/TMath::Power((-1.49187e+00)+x[0]/(4.14660e-07),1.71933e-01)+CrossOverRc(2.75563e+00,2.34650e+00,x[0])*1.64461e+00/TMath::Power((-2.58810e-01)+x/2.21576e+00,4.29324e+00)*CrossOverLc(5.67592e+00,(7.62946e-01)+CrossOverRc(5.67592e+00,7.62946e-01,x[0])*(5.49622e-01)/TMath::Power(x[0],(-2.15261e+00)),x[0]))/(CrossOverLc(3.38860e+00,4.05992e-01,x[0])*(-7.12382e-02)/TMath::Power((4.57260e-01)+x[0]/(5.36517e+02),4.07632e+00)+CrossOverRc(3.38860e+00,4.05992e-01,x[0])*(3.47577e-03)/TMath::Power(8.63986e-02+x/(2.06117e+01),2.58547e+00)*CrossOverLc(1.00000e+00,(9.99900e+03)+CrossOverRc(1.00000e+00,9.99900e+03,x[0])*(0.00000e+00)/TMath::Power(x[0],0.00000e+00),x[0]))",7,30);   
    fDmesonShape13->SetParameter(0,0.538845);
    
    TF1 *fDmesonShape14 = new TF1("fDmesonShape14","(1/[0])*(CrossOverLc(2.75563e+00,2.34650e+00,x[0])*(-4.22294e-01)/TMath::Power((-1.49187e+00)+x[0]/(4.14660e-07),1.71933e-01)+CrossOverRc(2.75563e+00,2.34650e+00,x[0])*1.64461e+00/TMath::Power((-2.58810e-01)+x/2.21576e+00,4.29324e+00)*CrossOverLc(5.67592e+00,(7.62946e-01)+CrossOverRc(5.67592e+00,7.62946e-01,x[0])*(5.49622e-01)/TMath::Power(x[0],(-2.15261e+00)),x[0]))/(CrossOverLc(3.38860e+00,4.05992e-01,x[0])*(-7.12382e-02)/TMath::Power((4.57260e-01)+x[0]/(5.36517e+02),4.07632e+00)+CrossOverRc(3.38860e+00,4.05992e-01,x[0])*(3.47577e-03)/TMath::Power(8.63986e-02+x/(2.06117e+01),2.58547e+00)*CrossOverLc(1.00000e+00,(9.99900e+03)+CrossOverRc(1.00000e+00,9.99900e+03,x[0])*(0.00000e+00)/TMath::Power(x[0],0.00000e+00),x[0]))",8,30);   
    fDmesonShape14->SetParameter(0,0.382748);
    
    TF1 *fDmesonShape15 = new TF1("fDmesonShape15","(1/[0])*(CrossOverLc(2.75563e+00,2.34650e+00,x[0])*(-4.22294e-01)/TMath::Power((-1.49187e+00)+x[0]/(4.14660e-07),1.71933e-01)+CrossOverRc(2.75563e+00,2.34650e+00,x[0])*1.64461e+00/TMath::Power((-2.58810e-01)+x/2.21576e+00,4.29324e+00)*CrossOverLc(5.67592e+00,(7.62946e-01)+CrossOverRc(5.67592e+00,7.62946e-01,x[0])*(5.49622e-01)/TMath::Power(x[0],(-2.15261e+00)),x[0]))/(CrossOverLc(3.38860e+00,4.05992e-01,x[0])*(-7.12382e-02)/TMath::Power((4.57260e-01)+x[0]/(5.36517e+02),4.07632e+00)+CrossOverRc(3.38860e+00,4.05992e-01,x[0])*(3.47577e-03)/TMath::Power(8.63986e-02+x/(2.06117e+01),2.58547e+00)*CrossOverLc(1.00000e+00,(9.99900e+03)+CrossOverRc(1.00000e+00,9.99900e+03,x[0])*(0.00000e+00)/TMath::Power(x[0],0.00000e+00),x[0]))",9,30);   
    fDmesonShape15->SetParameter(0,0.286206);
       
    TF1 *fDmesonShape_norm = new TF1("fDmesonShape_norm","(1/[0])*(CrossOverLc(2.75563e+00,2.34650e+00,x[0])*(-4.22294e-01)/TMath::Power((-1.49187e+00)+x[0]/(4.14660e-07),1.71933e-01)+CrossOverRc(2.75563e+00,2.34650e+00,x[0])*1.64461e+00/TMath::Power((-2.58810e-01)+x/2.21576e+00,4.29324e+00)*CrossOverLc(5.67592e+00,(7.62946e-01)+CrossOverRc(5.67592e+00,7.62946e-01,x[0])*(5.49622e-01)/TMath::Power(x[0],(-2.15261e+00)),x[0]))/(CrossOverLc(3.38860e+00,4.05992e-01,x[0])*(-7.12382e-02)/TMath::Power((4.57260e-01)+x[0]/(5.36517e+02),4.07632e+00)+CrossOverRc(3.38860e+00,4.05992e-01,x[0])*(3.47577e-03)/TMath::Power(8.63986e-02+x/(2.06117e+01),2.58547e+00)*CrossOverLc(1.00000e+00,(9.99900e+03)+CrossOverRc(1.00000e+00,9.99900e+03,x[0])*(0.00000e+00)/TMath::Power(x[0],0.00000e+00),x[0]))",2,50);   
    fDmesonShape_norm->SetParameter(0,8.50802);
    
    //=======================================================================
    ///Initialization for MC analysis:
    
    if(fIsMC)
    {
        if(fIsAOD)
        {
            fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
            if(!fMCarray){
                AliError("Array of MC particles not found");
                return;
            }
            if(fMCarray->GetEntries() < 1) return;
            
            fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
            if (!fMCheader) {
                AliError("Could not find MC Header in AOD");
                return;
            }
            
            Bool_t test = GetNMCPartProduced();
            //cout << "=====================" << endl;
            weight = 1;
            
             
            ////--------------------------------------------------------
            
            ///Get proportion between generated D0 and LambdaC
            for(Int_t iMC = 0; iMC < fNTotMCpart; iMC++)
            {
				fMCparticle = (AliAODMCParticle*) fMCarray->At(iMC);
				
				if((fMCparticle->Eta() < fEtaMin) || (fMCparticle->Eta() > fEtaMax)) continue;
                     				
				Int_t TrackPDG = TMath::Abs(fMCparticle->GetPdgCode());
                if(TrackPDG == 421) hPtD0->Fill(fMCparticle->Pt());///D0
                if(TrackPDG == 4122) hPtLambdaC->Fill(fMCparticle->Pt()); ///LambdaC
                
                
                ///For the beauty reconstruction efficiency-----------
				if(TrackPDG == 11){	
					Bool_t MotherFound = FindMother(iMC);
					if(fIsFromB){
                        fPtBeautyGenerated->Fill(fMCparticle->Pt());
					}
				}
				///----------------------------------------------------
                
                
				  
            }
            

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


        if(track->Pt()>30) continue;
        
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
        
        
        //RecKine: ITSTPC cuts --------------------------------------------------------------
        if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
        
       
        
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
        
              
        //RecPrim -------------------------------------------------------------------------
        if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;				///ProcessCutStep(Int_t cutStep, AliVParticle *track)
        
       
        //HFEcuts: ITS layers cuts -------------------------------------------------------------
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;
        
        
        //HFE cuts: TPC PID cleanup
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;
        
        //=======================================================================
        // QA plots after track selection
        //=======================================================================
        
        
        ///For the beauty reconstruction efficiency-----------
        if(fIsMC && fIsAOD){
			Bool_t IsHFEMC = IsHFelectronsMC(track);
			if(!IsHFEMC){
				if(fIsFromB){
					fPtBeautyReconstructedTracks->Fill(fPt);
				}
			}
		}
		
		///----------------------------------------------------
        
             
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
        
        ///Calculating DCA
        Double_t d0z0[2], cov[3];
        AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
        if(!(track->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 100., d0z0, cov))) continue;
        Double_t DCAxy = d0z0[0];
        Double_t DCAz = d0z0[1];
        
        //cout<<"fAOD->GetMagneticField() = "<<fAOD->GetMagneticField()<<endl;
        
               
        
        ///Checking nsigmaTPC after PID cuts 
        if(fTOFnSigma >= ftofPIDmincut && fTOFnSigma <= ftofPIDmaxcut){
            fTPCnsigma_p_after_tof->Fill(fP,fTPCnSigma);
            fTPCnsigma_pt_after_tof->Fill(fPt,fTPCnSigma);
            
             ///-------------------------------
			///Cheking the hadron nsigmaTPC after TOF cut
			if(fIsMC){
				fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
				if(TMath::Abs(fMCparticle->GetPdgCode()) == 2212)  fTPCnsigma_p_after_tof_p->Fill(fP,fTPCnSigma);
				if(TMath::Abs(fMCparticle->GetPdgCode()) == 130 || TMath::Abs(fMCparticle->GetPdgCode()) == 310 || TMath::Abs(fMCparticle->GetPdgCode()) == 311 || TMath::Abs(fMCparticle->GetPdgCode()) == 321)  fTPCnsigma_p_after_tof_k->Fill(fP,fTPCnSigma);
				if(TMath::Abs(fMCparticle->GetPdgCode()) == 211)  fTPCnsigma_p_after_tof_pion->Fill(fP,fTPCnSigma);
			}
			///-------------------------------        
            
            if(fITSnSigma >= -2 && fITSnSigma <= 2){
                fTPCnsigma_p_after_tof_its->Fill(fP,fTPCnSigma);
                fTPCnsigma_pt_after_tof_its->Fill(fPt,fTPCnSigma);
            }
        }
        
              
        if(fTPCnSigma >= -5 && fTPCnSigma <= -3){
			fDCAxy_pt_had->Fill(fPt,DCAxy*track->Charge()*signB);
			fDCAz_pt_had->Fill(fPt,DCAz);
        }
        
		///Using thnSparse to store the DCA in Data-----------------------
         if(!fIsMC){
			 qadcaData[0] = fPt;
         
			 qadcaData[1] = DCAxy*track->Charge()*signB;
			
			 qadcaData[4] = -1.;
			
			 ///Charged pions
			 if(fTPCnSigma >= -5 && fTPCnSigma <= -3){
				  qadcaData[4] = 0.5;
			 }
			 ///Electron candidates
			 if(fTPCnSigma >= ftpcPIDmincut && fTPCnSigma <= ftpcPIDmaxcut){
				if(fTOFnSigma >= ftofPIDmincut && fTOFnSigma <= ftofPIDmaxcut){
					qadcaData[4] = 1.5;					
				}
			 }
        
			Double_t ITSNcls = atrack->GetITSNcls();
			//cout<<"atrack->GetITSNcls() = "<<atrack->GetITSNcls()<<endl;
            
			///ITS Chi2            
			qadcaData[2] = atrack->GetITSchi2()/ITSNcls; 
			//cout<<"track->GetITSchi2() = "<<track->GetITSchi2()<<endl;
            
			///Fraction of shared clusters in the ITS
			Bool_t HasSharedCls = kFALSE;
			Double_t ITSNSharedcls = 0;
			for(int itsL = 0; itsL < 6; itsL++){
				HasSharedCls = atrack->HasSharedPointOnITSLayer(itsL);
				if(HasSharedCls) ITSNSharedcls++;
			}
			//cout<<"ITSNSharedcls = "<<ITSNSharedcls<<endl;
            
                        
			Double_t fsharedclsITS = ITSNSharedcls/ITSNcls;
			//cout<<"fsharedclsITS = "<<fsharedclsITS<<endl;
            
			qadcaData[3] = fsharedclsITS; 
         
			if(qadcaData[4]>0.) fD0Data->Fill(qadcaData);
        }
        ///-------------------------------------------------------------          
        
               
        //=======================================================================
        // Here the PID cuts defined in the file "Config.C" is applied
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
        
        
        ///For the beauty reconstruction efficiency-----------
         if(fIsMC && fIsAOD){
			Bool_t IsHFEMC = IsHFelectronsMC(track);
			if(!IsHFEMC){
				if(fIsFromB){
					fPtBeautyReconstructedTracksPID->Fill(fPt);
				}
			}
		}
		///----------------------------------------------------
        
        
        ///After TOF and TPC cuts
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
 
        ///DCA electron candidates
        fDCAxy_pt_ele->Fill(fPt,DCAxy*track->Charge()*signB);
        fDCAz_pt_ele->Fill(fPt,DCAz);
              
        
        ///Using thnSparse to store the DCA in MC-----------------------
        if(fIsMC && fIsAOD){
            
            Int_t trkLabel = TMath::Abs(track->GetLabel());
            Int_t labelm = GetPrimary(trkLabel,fMCarray);///gives the label of first mother
            AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCarray->At(labelm);
            Int_t trkIndexPrimHFE= AODMCtrack->GetLabel();///gives index of the particle in original MCparticle array (labelm and trkIndexPrimHFE are the same, so I don't understand why this is done)
            
            ///Electron generator:
            Int_t MChijingHFE=kHijing;
            if(trkIndexPrimHFE >= fNpureMC) MChijingHFE = kPhytia;///check if the particle comes from hijing or from enhanc
            //cout<<"trkIndexPrimHFE = "<<trkIndexPrimHFE<<endl;
                        
            qadca[3]=MChijingHFE;
			
			///Position where the electron is created:
			fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
			Double_t fVx = fMCparticle->Xv();
			Double_t fVy = fMCparticle->Yv();
			Double_t Rconv = TMath::Sqrt(fVx*fVx+fVy*fVy);
			qadca[4]=Rconv; 
			
			///Pt
            qadca[0] = fPt;
            
            ///Selecting electron source
            qadca[1]=-1.;
            
            ///------------
            if(fMCparticle->GetPdgCode() == 2212) qadca[1]=10.5; ///to check DCA of protons
            if(fMCparticle->GetPdgCode() == 321) qadca[1]=11.5; ///to check DCA of kaons
            ///------------
            
            
            Bool_t IsHFEMC = IsHFelectronsMC(track);
            
            if(IsHFEMC){
                if(fIsFromD){
                    fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
                    fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                    pdg_mother = fMCparticleMother->GetPdgCode();
                    if(TMath::Abs(pdg_mother)>400 && TMath::Abs(pdg_mother)<500){///charmed meson 
						 qadca[1]=0.5;
						 hCharmMotherPt->Fill(fMCparticleMother->Pt());
						 
						 hCharmMotherPt_vsElecPt->Fill(fPt,fMCparticleMother->Pt());
						 hElecPt_vsCharmMotherPt->Fill(fMCparticleMother->Pt(),fPt);
						 
						 ///(Weight for each pT bin starting at one)
						 float probAcceptD = -999;
						 if(fMCparticleMother->Pt() > 30) continue;
						 if(fPt >=2 && fPt < 2.25) probAcceptD = fDmesonShape1->Eval(fMCparticleMother->Pt());
						 if(fPt >=2.25 && fPt < 2.5) probAcceptD = fDmesonShape2->Eval(fMCparticleMother->Pt());
						 if(fPt >=2.5 && fPt < 2.75) probAcceptD = fDmesonShape3->Eval(fMCparticleMother->Pt());
						 if(fPt >=2.75 && fPt < 3) probAcceptD = fDmesonShape4->Eval(fMCparticleMother->Pt());
						 if(fPt >=3 && fPt < 3.5) probAcceptD = fDmesonShape5->Eval(fMCparticleMother->Pt());
						 if(fPt >=3.5 && fPt < 4) probAcceptD = fDmesonShape6->Eval(fMCparticleMother->Pt());
						 if(fPt >=4 && fPt < 4.5) probAcceptD = fDmesonShape7->Eval(fMCparticleMother->Pt());
						 if(fPt >=4.5 && fPt < 5) probAcceptD = fDmesonShape8->Eval(fMCparticleMother->Pt());
						 if(fPt >=5 && fPt < 5.5) probAcceptD = fDmesonShape9->Eval(fMCparticleMother->Pt());
						 if(fPt >=5.5 && fPt < 6) probAcceptD = fDmesonShape10->Eval(fMCparticleMother->Pt());
						 if(fPt >=6 && fPt < 6.5) probAcceptD = fDmesonShape11->Eval(fMCparticleMother->Pt());
						 if(fPt >=6.5 && fPt < 7) probAcceptD = fDmesonShape12->Eval(fMCparticleMother->Pt());
						 if(fPt >=7 && fPt < 8) probAcceptD = fDmesonShape13->Eval(fMCparticleMother->Pt());
						 if(fPt >=8 && fPt < 9) probAcceptD = fDmesonShape14->Eval(fMCparticleMother->Pt());
						 if(fPt >=9 && fPt < 10) probAcceptD = fDmesonShape15->Eval(fMCparticleMother->Pt());
						 
											 
						 ///(A single weight that starts at one in the first pT bin)
						 float probAcceptD3 = -999;
						 if(fPt >=2 && fPt < 10) probAcceptD3 = fDmesonShape_norm->Eval(fMCparticleMother->Pt());
						 
						 float a = gRandom->Uniform(0,1);
						 //float a = fconstant->GetRandom();
						 //cout<<"fMCparticleMother->Pt() = "<<fMCparticleMother->Pt()<<endl;
						 //cout<<"probAcceptD = "<<probAcceptD<<endl;
						 //cout<<"a = "<<a<<endl;
						 if(a < probAcceptD){
							 qadca[1]=7.5;
							 hCharmMotherPt_corr->Fill(fMCparticleMother->Pt());
							 hCharmMotherPt_corr2->Fill(fMCparticleMother->Pt());
							 
							if(fPt >=2 && fPt < 2.25) hCharmMotherPt1_corr->Fill(fMCparticleMother->Pt());
							if(fPt >=2.25 && fPt < 2.5) hCharmMotherPt2_corr->Fill(fMCparticleMother->Pt());
							if(fPt >=2.5 && fPt < 2.75) hCharmMotherPt3_corr->Fill(fMCparticleMother->Pt());
							if(fPt >=2.75 && fPt < 3) hCharmMotherPt4_corr->Fill(fMCparticleMother->Pt());
							if(fPt >=3 && fPt < 3.5) hCharmMotherPt5_corr->Fill(fMCparticleMother->Pt());
							if(fPt >=3.5 && fPt < 4) hCharmMotherPt6_corr->Fill(fMCparticleMother->Pt());
							if(fPt >=4 && fPt < 4.5) hCharmMotherPt7_corr->Fill(fMCparticleMother->Pt());
							if(fPt >=4.5 && fPt < 5) hCharmMotherPt8_corr->Fill(fMCparticleMother->Pt());
							if(fPt >=5 && fPt < 5.5) hCharmMotherPt9_corr->Fill(fMCparticleMother->Pt());
							if(fPt >=5.5 && fPt < 6) hCharmMotherPt10_corr->Fill(fMCparticleMother->Pt());
							if(fPt >=6 && fPt < 6.5) hCharmMotherPt11_corr->Fill(fMCparticleMother->Pt());
							if(fPt >=6.5 && fPt < 7) hCharmMotherPt12_corr->Fill(fMCparticleMother->Pt());
							if(fPt >=7 && fPt < 8) hCharmMotherPt13_corr->Fill(fMCparticleMother->Pt());
							if(fPt >=8 && fPt < 9) hCharmMotherPt14_corr->Fill(fMCparticleMother->Pt());
							if(fPt >=9 && fPt < 10) hCharmMotherPt15_corr->Fill(fMCparticleMother->Pt());
							 
							 
							 
						 }
						 
						 if(a < probAcceptD3){
							 hCharmMotherPt_corr3->Fill(fMCparticleMother->Pt());
							 hCharmMotherPt_corr4->Fill(fMCparticleMother->Pt());
							 qadca[1]=8.5;	
						 }
						 
					}
					
					if(TMath::Abs(pdg_mother)>4000 && TMath::Abs(pdg_mother)<5000){///charmed baryon
						 qadca[1]=1.5;
					}
                }
                
                if(fIsFromB){
                    fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
                    fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                    pdg_mother = fMCparticleMother->GetPdgCode();
                    if(TMath::Abs(pdg_mother)>500 && TMath::Abs(pdg_mother)<600){///beauty meson 
						qadca[1]=2.5;
						hBeautyMotherPt->Fill(fMCparticleMother->Pt());
					}
					if(TMath::Abs(pdg_mother)>5000 && TMath::Abs(pdg_mother)<6000){///beauty baryon
						qadca[1]=3.5;
					}
                }
            }
           
           if(!IsHFEMC){
		
                fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
                pdg = fMCparticle->GetPdgCode();
                ///Is electron:
                if(TMath::Abs(pdg) != 11) continue;
                if(fMCparticle->GetMother()<=0) continue;
                fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                pdg_mother = fMCparticleMother->GetPdgCode();
                
                ///Photonic Electrons:
                if(TMath::Abs(pdg_mother) == 111 || TMath::Abs(pdg_mother) == 221){
					
					if(TMath::Abs(pdg_mother) == 111) qadca[1]=4.5; 
					if(TMath::Abs(pdg_mother) == 221) qadca[1]=5.5;
					Int_t fType = GetPi0EtaType(fMCparticleMother,fMCarray);
					qadca[7]=fType;
                }
                
                if(TMath::Abs(pdg_mother) == 22){
                    Bool_t primMC = fMCparticleMother->IsPrimary();
                    if(primMC) qadca[1]=6.5;
                }
            }
            
                       
            ///DCAxy
            qadca[2]=DCAxy*track->Charge()*signB;
            
            
            Double_t ITSNcls = atrack->GetITSNcls();
            //cout<<"atrack->GetITSNcls() = "<<atrack->GetITSNcls()<<endl;
            
            ///ITS Chi2            
            qadca[5] = atrack->GetITSchi2()/ITSNcls; 
            //cout<<"track->GetITSchi2() = "<<track->GetITSchi2()<<endl;
            
            ///Fraction of shared clusters in the ITS
            Bool_t HasSharedCls = kFALSE;
            Double_t ITSNSharedcls = 0;
            for(int itsL = 0; itsL < 6; itsL++){
				HasSharedCls = atrack->HasSharedPointOnITSLayer(itsL);
				if(HasSharedCls) ITSNSharedcls++;
			}
			//cout<<"ITSNSharedcls = "<<ITSNSharedcls<<endl;
            
                        
            Double_t fsharedclsITS = ITSNSharedcls/ITSNcls;
            //cout<<"fsharedclsITS = "<<fsharedclsITS<<endl;
            
            qadca[6] = fsharedclsITS; 
                        
            ///Fill
            if(qadca[1]>0.) fD0->Fill(qadca);
            
            
        }
        ///-------------------------------------------------------------
       

    }//End of track loop
    
    //=======================================================================
    
    delete fListOfmotherkink;
    PostData(1, fOutputList);
}

//=======================================================================
void AliAnalysisHFETPCTOFBeauty::Terminate(Option_t *)
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
Bool_t AliAnalysisHFETPCTOFBeauty::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
    //Check single track cuts for a given cut step
    //Note this function is called inside the UserExec function
    const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
    if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
    return kTRUE;
}

//Setter for the PID cuts (TOF and TPC)
void AliAnalysisHFETPCTOFBeauty::SetPIDCuts(Float_t tpcPIDmincut, Float_t tpcPIDmaxcut, Float_t tofPIDmincut, Float_t tofPIDmaxcut) {
    ftpcPIDmincut = tpcPIDmincut;
    ftpcPIDmaxcut = tpcPIDmaxcut;
    ftofPIDmincut = tofPIDmincut;
    ftofPIDmaxcut = tofPIDmaxcut;
}
//Setter for the Eta cut
void AliAnalysisHFETPCTOFBeauty::SetEtaCut(Float_t EtaMin, Float_t EtaMax){
    fEtaMin = EtaMin;
    fEtaMax = EtaMax;
}


//=======================================================================
Bool_t AliAnalysisHFETPCTOFBeauty::FindMother(Int_t mcIndex)
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
        
        if((mpdg>400 && mpdg<500) || (mpdg>4000 && mpdg<5000)) //charmed mesons and baryons
        {
            if((gmpdg>500 && gmpdg<600) || (ggmpdg>500 && ggmpdg<600) || (gggmpdg>500 && gggmpdg<600) || (gmpdg>5000 && gmpdg<6000) || (ggmpdg>5000 && ggmpdg<6000) || (gggmpdg>5000 && gggmpdg<6000)) //when the charm comes from beauty
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
        else if((mpdg>500 && mpdg<600) || (mpdg>5000 && mpdg<6000)) //beauty mesons and baryons
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

Bool_t AliAnalysisHFETPCTOFBeauty::IsHFelectronsMC(AliVTrack *track){
    
    
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
Int_t AliAnalysisHFETPCTOFBeauty::GetPi0EtaType(AliAODMCParticle *pi0eta, TClonesArray *mcArray){
    //
    // Return what type of pi0, eta it is
    //
    
    // IsPrimary
    Bool_t primMC = pi0eta->IsPrimary(); ///Does not include particles from weak decays or created in an interaction with the material
    if(!primMC) return kNoIsPrimary;
    
    // Mother
    Int_t motherlabel = pi0eta->GetMother();
    if(motherlabel<0) return kNoMother;
    else {
        
        AliAODMCParticle *mother = (AliAODMCParticle*)fMCarray->At(motherlabel);
        Int_t motherpdg = TMath::Abs(mother->GetPdgCode());
        ///pi0, eta, omega, phi, eta',rho0, rho+,k*0,K*+,lambda(strangeness)
        if(motherpdg == 111 || motherpdg == 221 || motherpdg == 223 || motherpdg == 333 || motherpdg == 331 || motherpdg == 113 || motherpdg == 213 || motherpdg == 313 || motherpdg == 323 || motherpdg == 3122) return kLightMesons;
        
        if ( (int(TMath::Abs(motherpdg)/100.)%10) == 5 || (int(TMath::Abs(motherpdg)/1000.)%10) == 5 ) return kBeauty; ///beauty mesons and barions
        if ( (int(TMath::Abs(motherpdg)/100.)%10) == 4 || (int(TMath::Abs(motherpdg)/1000.)%10) == 4 ) return kCharm; ///charmed mesons and barions
        return kNoFeedDown;    
        
    }
}


//___________________________________________________________________________________________________________

Int_t AliAnalysisHFETPCTOFBeauty::GetElecSourceType(AliAODMCParticle *electron, TClonesArray *mcArray,Double_t &ptm){
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
            if(typepi0eta==kNoMother){ 	
            
            //Int_t motherlabelpi0 = mother->GetMother();
            //cout<<"motherlabelpi0 = "<<motherlabelpi0<<endl;
            
			return kPi0NoFeedDown;
			}
        }
        if(motherpdg == 221){
            Int_t typepi0eta = GetPi0EtaType(mother,mcArray);
            if(typepi0eta==kNoMother) return kEtaNoFeedDown;
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
                    if(typepi0eta==kNoMother) return kGPi0NoFeedDown;
                }
                if(gmotherpdg == 221) {
                    Int_t typepi0eta = GetPi0EtaType(gmother,mcArray);
                    if(typepi0eta==kNoMother) return kGEtaNoFeedDown;
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
                            if(typepi0eta==kNoMother) return kGPi0NoFeedDown;
                        }
                        if(ggmotherpdg == 221) {
                            Int_t typepi0eta = GetPi0EtaType(ggmother,mcArray);
                            if(typepi0eta==kNoMother) return kGEtaNoFeedDown;
                        }
                    }
                }
            }
        }
    }
    
    return kOthersE;
    
}

//_______________________________________________________________________________________
Bool_t AliAnalysisHFETPCTOFBeauty::GetNMCPartProduced()
{
    //Get number of MC particles produced by generators.
    
    TList *lh = fMCheader->GetCocktailHeaders();
    fNTotMCpart = 0;
    fNembMCpi0 = 0;
    fNembMCeta = 0;
    fNpureMC = 0;
    
    TString MCgen;
    TString embpi0("pi");
    TString embeta("eta");
    
    if(!lh){
        AliError("no MC header");
        return (0);
    }
    
    
    for(int igene=0; igene<lh->GetEntries(); igene++)
    {
        AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(igene);
        if(!gh) continue;
        
        MCgen =  gh->GetName();
        //cout << "Gen name, N produced = " << gh->GetName() << ", " << gh->NProduced() << endl;
        
        if(igene==0) fNpureMC = gh->NProduced();  // generated by HIJING
        //if(MCgen.Contains(embpi0))cout << MCgen << endl;
        //if(MCgen.Contains(embeta))cout << MCgen << endl;
        
        if(MCgen.Contains(embpi0))fNembMCpi0 = fNTotMCpart;
        if(MCgen.Contains(embeta))fNembMCeta = fNTotMCpart;
        
        fNTotMCpart += gh->NProduced();
    }
    //cout << "fNpureMC, fNembMCpi0, fNembMCeta, fNTotMCpart : " <<fNpureMC << ", " << fNembMCpi0 << ", " << fNembMCeta << ", " << fNTotMCpart << endl;
    
    return kTRUE;
}

//_________________________________________

Int_t AliAnalysisHFETPCTOFBeauty::GetPrimary(Int_t id, TClonesArray *mcArray){
    
    // Return number of primary that has generated track
    int current, parent;
    parent=id;
    while (1) {
        current=parent;
        AliAODMCParticle *Part = (AliAODMCParticle*)mcArray->At(current);
        parent=Part->GetMother();
        //  cout << "GetPartArr momid :"  << parent << endl;
        if(parent<0) return current;
    }
}
//_________________________________________

//_________________________________________
Bool_t AliAnalysisHFETPCTOFBeauty::PassCorrCuts(AliAODEvent *fAOD)
{
  //event selection cuts
	if(!fIsPP){
		Int_t ntracks = -999;
		ntracks = fAOD->GetNumberOfTracks();
		if(ntracks < 1) return kFALSE;

		AliAODVertex* vtTrc = fAOD->GetPrimaryVertex();
		Double_t NcontV = vtTrc->GetNContributors();
		AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();
		Double_t NcontSPD = vtSPD->GetNContributors();
		if(NcontV<2 || NcontSPD<1)return kFALSE;

		Double_t covTrc[6],covSPD[6];
		vtTrc->GetCovarianceMatrix(covTrc);
		vtSPD->GetCovarianceMatrix(covSPD);
		Double_t dz = vtTrc->GetZ() - vtSPD->GetZ();
		Double_t errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
		Double_t errTrc = TMath::Sqrt(covTrc[5]);
		Double_t nsigTot = TMath::Abs(dz)/errTot;
		Double_t nsigTrc = TMath::Abs(dz)/errTrc;
		if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20) return kFALSE;// bad vertexing
		return kTRUE;
	}
	
	if(fIsPP){
		AliAODVertex* vtTrc = fAOD->GetPrimaryVertex();
		if(vtTrc->GetNContributors()<2) return kFALSE;
        
        Bool_t isPileupfromSPDmulbins=fAOD->IsPileupFromSPDInMultBins();
        if(isPileupfromSPDmulbins) return kFALSE;
        
        ///minContributors=5; minChi2=5.; minWeiZDiff=15; checkPlpFromDifferentBC=kFALSE;
        AliAnalysisUtils utils;
        utils.SetMinPlpContribMV(5);
        utils.SetMaxPlpChi2MV(5.);
        utils.SetMinWDistMV(15);
        utils.SetCheckPlpFromDifferentBCMV(kFALSE);
        Bool_t isPileupFromMV = utils.IsPileUpMV(fAOD);
		if(isPileupFromMV) return kFALSE;
		return kTRUE;
	}
  
}

Double_t CrossOverLc(const double a, const double b, const double x){
    if(x<b-a/2) return 1.0;
    else if(x>b+a/2) return 0.0;
    else return cos(((x-b)/a+0.5)*TMath::Pi())/2+0.5;
}
Double_t CrossOverRc(const double a, const double b, const double x){
    return 1-CrossOverLc(a,b,x);
}
//_________________________________________





