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
//      Task for Beauty analysis in Pb-Pb central collisions   				  //
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
#include "AliHFEextraCuts.h"
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
#include "AliAnalysisTaskSEImproveITS.h"
#include "AliHFEV0taginfo.h"
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
,fIsFromBarionB(kFALSE)
,fIsFromMesonB(kFALSE)
,fIsFromBarionBD(kFALSE)
,fIsFromMesonBD(kFALSE)
,fIsFromPi0(kFALSE)
,fIsFromEta(kFALSE)
,fIsFromGamma(kFALSE)

//General variables
,fESD(0)
,fAOD(0)
,fVevent(0)
,fOutputList(0)
,fPidResponse(0)
,fExtraCuts(NULL)
,fV0Tagger(NULL)
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
,fNevent_passvertex(0)
,fNevent_corrcut(0)
,fNevent_no_vertex(0)
,fNevent_no_vertex_2(0)
,fNeventAnalized(0)
//,fNevent2(0)
,fCent(0)
,fCent2(0)
,fTPC_p1(0)
,fTPC_p2(0)
,fTPC_p3(0)
,fPt_1(0)
,fPt_2(0)
,fNAnalizedTracks(0)
,fNAnalizedTracksHijing(0)
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
,fTPCnsigma_p_after_V0selection(0)
,fTPCnsigma_p_after_tof(0)
,fTPCnsigma_p_after_tof_v2(0)
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
,fDCAxy_pt_charmbef(0)
,fDCAxy_pt_charmaft(0)
,fDCAxy_pt_beautybef(0)
,fDCAxy_pt_beautyaft(0)
,fDCAxy_pt_had_onlyDCA(0)
,fDCAxy_pt_had_onlyDCA_Hijing(0)
,fDCAxy_pt_had_onlyDCA_Phytia(0)
,fDCAz_pt_had(0)
,fDCAxy_pt_ele(0)
,fDCAz_pt_ele(0)
,hCharmMotherPt(0)
,hCharmMotherPt_vsElecPt(0)
,hElecPt_vsCharmMotherPt(0)
,hCharmMotherPt_vsElecPt_corr(0)
,hElecPt_vsCharmMotherPt_corr(0)
,hCharmMotherPt_corr(0)
,hCharmMotherPt_corr2(0)
,hBeautyMotherPt(0)
,fPtBeautyGenerated(0)
,fPtGeneratedBmesons(0)
,fPtBeautyReconstructedAll(0)
,fPtBeautyReconstructedTracks(0)
,fPtBeautyReconstructedTracksPID(0)
,fPtBeautyReconstructedTracksPIDTPC(0)
,fPtBeautyReconstructedTracksPIDTOF(0)
,fPtBeautyReconstructedTracksPIDITS(0)
,hTOFEffDen(0)
,hTOFEffNum(0)
,fPtBeautyPtrecVsPtparticle(0)
,hPtLambdaC(0)
,hPtD0(0)

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
,fBcorr(0)
,fDcorr1(0)
,fDcorr2(0)
,fDcorr3(0)
,fDcorr4(0)
,fDcorr5(0)
,fDcorr6(0)
,fDcorr7(0)
,fDcorr8(0)
,fDcorr9(0)
,fDcorr10(0)
,fDcorr11(0)
,fDcorr12(0)
,fDcorr13(0)
,fDcorr14(0)
,fDcorr15(0)
,fDcorr16(0)
,fDcorr17(0)
,fDcorr18(0)
,fDcorr19(0)
,fDcorr20(0)
,fDcorr21(0)
,fDcorr22(0)

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
,fIsFromBarionB(kFALSE)
,fIsFromMesonB(kFALSE)
,fIsFromBarionBD(kFALSE)
,fIsFromMesonBD(kFALSE)
,fIsFromPi0(kFALSE)
,fIsFromEta(kFALSE)
,fIsFromGamma(kFALSE)

//General variables
,fESD(0)
,fAOD(0)
,fVevent(0)
,fOutputList(0)
,fPidResponse(0)
,fExtraCuts(NULL)
,fV0Tagger(NULL)
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
,fNevent_passvertex(0)
,fNevent_corrcut(0)
,fNevent_no_vertex(0)
,fNevent_no_vertex_2(0)
,fNeventAnalized(0)
,fCent(0)
,fCent2(0)
,fTPC_p1(0)
,fTPC_p2(0)
,fTPC_p3(0)
,fPt_1(0)
,fPt_2(0)
,fNAnalizedTracks(0)
,fNAnalizedTracksHijing(0)
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
,fTPCnsigma_p_after_V0selection(0)
,fTPCnsigma_p_after_tof(0)
,fTPCnsigma_p_after_tof_v2(0)
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
,fDCAxy_pt_charmbef(0)
,fDCAxy_pt_charmaft(0)
,fDCAxy_pt_beautybef(0)
,fDCAxy_pt_beautyaft(0)
,fDCAxy_pt_had_onlyDCA(0)
,fDCAxy_pt_had_onlyDCA_Hijing(0)
,fDCAxy_pt_had_onlyDCA_Phytia(0)
,fDCAz_pt_had(0)
,fDCAxy_pt_ele(0)
,fDCAz_pt_ele(0)
,hCharmMotherPt(0)
,hCharmMotherPt_vsElecPt(0)
,hElecPt_vsCharmMotherPt(0)
,hCharmMotherPt_vsElecPt_corr(0)
,hElecPt_vsCharmMotherPt_corr(0)
,hCharmMotherPt_corr(0)
,hCharmMotherPt_corr2(0)
,hBeautyMotherPt(0)
,fPtBeautyGenerated(0)
,fPtGeneratedBmesons(0)
,fPtBeautyReconstructedAll(0)
,fPtBeautyReconstructedTracks(0)
,fPtBeautyReconstructedTracksPID(0)
,fPtBeautyReconstructedTracksPIDTPC(0)
,fPtBeautyReconstructedTracksPIDTOF(0)
,fPtBeautyReconstructedTracksPIDITS(0)
,hTOFEffDen(0)
,hTOFEffNum(0)
,fPtBeautyPtrecVsPtparticle(0)
,hPtLambdaC(0)
,hPtD0(0)

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
,fBcorr(0)
,fDcorr1(0)
,fDcorr2(0)
,fDcorr3(0)
,fDcorr4(0)
,fDcorr5(0)
,fDcorr6(0)
,fDcorr7(0)
,fDcorr8(0)
,fDcorr9(0)
,fDcorr10(0)
,fDcorr11(0)
,fDcorr12(0)
,fDcorr13(0)
,fDcorr14(0)
,fDcorr15(0)
,fDcorr16(0)
,fDcorr17(0)
,fDcorr18(0)
,fDcorr19(0)
,fDcorr20(0)
,fDcorr21(0)
,fDcorr22(0)


{
    // Constructor
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    // DefineOutput(1, TH1I::Class());
    DefineOutput(1, TList::Class());
    fV0Tagger = new AliHFEV0taginfo("Tagger");
}

//______________________________________________________________________
AliAnalysisHFETPCTOFBeauty::~AliAnalysisHFETPCTOFBeauty()
{
    //Destructor
    delete fOutputList;
    delete fPID;
    delete fCFM;
    delete fPIDqa;
    delete fV0Tagger;
}

//______________________________________________________________________
//Create Output Objects
//Here we can define the histograms and others output files
//Called once
void AliAnalysisHFETPCTOFBeauty::UserCreateOutputObjects()
{	
    ///Initialize PID
    if(!fPID->GetNumberOfPIDdetectors())
    {
        fPID->AddDetector("TPC", 0);
        fPID->AddDetector("TOF", 1);
    }
    
    fPID->SortDetectors();
    
    fPIDqa = new AliHFEpidQAmanager();
    fPIDqa->Initialize(fPID);
        
    ///Initialize correction Framework and Cuts
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
    
    
    
    ///Output Tlist
    //Create TList
    fOutputList = new TList();
    fOutputList->SetOwner();
    
    //PIDqa
    fOutputList->Add(fPIDqa->MakeList("PIDQA"));
    
    //Define the histo
    fNevent = new TH1F("fNevent","Number of Events",5,-0.5,4.5);
    fVertex1 = new TH1F("fVertex1","vertex",500,-50,50);
    fNevent_passvertex = new TH1F("fNevent_passvertex","Number of Events",5,-0.5,4.5);
    fNevent_corrcut = new TH1F("fNevent_corrcut","Number of Events",5,-0.5,4.5);
    fNevent_no_vertex = new TH1F("fNevent_no_vertex","Number of Events",5,-0.5,4.5);
    fNevent_no_vertex_2 = new TH1F("fNevent_no_vertex_2","Number of Events",5,-0.5,4.5);
    fNeventAnalized = new TH1F("fNeventAnalized","Number of Events",5,-0.5,4.5);
    fCent = new TH1F("fCent","Centrality",100,0,100);
    fCent2 = new TH1F("fCent2","Centrality",100,0,100);
    //And then, add to the output list
    fOutputList->Add(fNevent);
    fOutputList->Add(fVertex1);
    fOutputList->Add(fNevent_passvertex);
    fOutputList->Add(fNevent_corrcut);
    fOutputList->Add(fNevent_no_vertex);
    fOutputList->Add(fNevent_no_vertex_2);
    fOutputList->Add(fNeventAnalized);
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
    
    fNAnalizedTracks = new TH1F("fNAnalizedTracks","pt (GeV/c)",5000,0,5000);
    fOutputList->Add(fNAnalizedTracks);
    
    fNAnalizedTracksHijing = new TH1F("fNAnalizedTracksHijing","pt (GeV/c)",5000,0,5000);
    fOutputList->Add(fNAnalizedTracksHijing);
    
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
    
    hCharmMotherPt_vsElecPt = new TH2F("hCharmMotherPt_vsElecPt","; p_{T} [GeV/c]; Count",1000,0,50,1000,0,50);
    fOutputList->Add(hCharmMotherPt_vsElecPt);
    
    hElecPt_vsCharmMotherPt = new TH2F("hElecPt_vsCharmMotherPt","; p_{T} [GeV/c]; Count",1000,0,50,1000,0,50);
    fOutputList->Add(hElecPt_vsCharmMotherPt);
    
    hCharmMotherPt_vsElecPt_corr = new TH2F("hCharmMotherPt_vsElecPt_corr","; p_{T} [GeV/c]; Count",1000,0,50,1000,0,50);
    fOutputList->Add(hCharmMotherPt_vsElecPt_corr);
    
    hElecPt_vsCharmMotherPt_corr = new TH2F("hElecPt_vsCharmMotherPt_corr","; p_{T} [GeV/c]; Count",1000,0,50,1000,0,50);
    fOutputList->Add(hElecPt_vsCharmMotherPt_corr);
    
    hCharmMotherPt_corr = new TH1F("hCharmMotherPt_corr","; p_{T} [GeV/c]; Count",13,ptbinningHF);
    fOutputList->Add(hCharmMotherPt_corr);
    
    hCharmMotherPt_corr2 = new TH1F("hCharmMotherPt_corr2","; p_{T} [GeV/c]; Count",100,0,50);
    fOutputList->Add(hCharmMotherPt_corr2);
    
	hBeautyMotherPt = new TH2F("hBeautyMotherPt","; p_{T} [GeV/c]; Count",1000,0,50,1000,0,50);
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
    
    fDCAxy_pt_charmbef = new TH2F("fDCAxy_pt_charmbef",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_pt_charmbef);
    
    fDCAxy_pt_charmaft = new TH2F("fDCAxy_pt_charmaft",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_pt_charmaft);
    
    fDCAxy_pt_beautybef = new TH2F("fDCAxy_pt_beautybef",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_pt_beautybef);
    
    fDCAxy_pt_beautyaft = new TH2F("fDCAxy_pt_beautyaft",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_pt_beautyaft); 
    
    fDCAxy_pt_had_onlyDCA = new TH2F("fDCAxy_pt_had_onlyDCA",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_pt_had_onlyDCA);
    
    fDCAxy_pt_had_onlyDCA_Hijing = new TH2F("fDCAxy_pt_had_onlyDCA_Hijing",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_pt_had_onlyDCA_Hijing);
       
    fDCAxy_pt_had_onlyDCA_Phytia = new TH2F("fDCAxy_pt_had_onlyDCA_Phytia",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_pt_had_onlyDCA_Phytia);
       
    fDCAz_pt_had = new TH2F("fDCAz_pt_had",";p_{t} (GeV/c);DCAz hadrons",300,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAz_pt_had);
    
    fDCAxy_pt_ele = new TH2F("fDCAxy_pt_ele",";p_{t} (GeV/c);DCAxy ",300,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAxy_pt_ele);
    
    fDCAz_pt_ele = new TH2F("fDCAz_pt_ele",";p_{t} (GeV/c);DCAz ",300,0,30,800,-0.2,0.2);
    fOutputList->Add(fDCAz_pt_ele); 
    
    fPtMCeta = new TH1F("fPtMCeta",";p_{t} (GeV/c)",2000,0,100);
    fOutputList->Add(fPtMCeta);
       
    fTPCnsigma_p_after_V0selection = new TH2F("fTPCnsigma_p_after_V0selection","p (GeV/c);TPC Electron N#sigma after V0 selection",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_p_after_V0selection);
    
    fTPCnsigma_p_after_tof = new TH2F("fTPCnsigma_p_after_tof","p (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_p_after_tof);
    
    fTPCnsigma_p_after_tof_v2 = new TH2F("fTPCnsigma_p_after_tof_v2","p (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_p_after_tof_v2);
    
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
    
    fPtGeneratedBmesons = new TH1F("fPtGeneratedBmesons","; p_{T} [GeV/c]; Count",1000,0,200);
    fOutputList->Add(fPtGeneratedBmesons);
   
    fPtBeautyReconstructedAll = new TH1F("fPtBeautyReconstructedAll","; p_{T} [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPtBeautyReconstructedAll);
    
    fPtBeautyReconstructedTracks = new TH1F("fPtBeautyReconstructedTracks","; p_{T} [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPtBeautyReconstructedTracks);
    
    fPtBeautyReconstructedTracksPID = new TH1F("fPtBeautyReconstructedTracksPID","; p_{T} [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPtBeautyReconstructedTracksPID);
    
    fPtBeautyReconstructedTracksPIDTPC = new TH1F("fPtBeautyReconstructedTracksPIDTPC","; p_{T} [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPtBeautyReconstructedTracksPIDTPC);
    
    fPtBeautyReconstructedTracksPIDTOF = new TH1F("fPtBeautyReconstructedTracksPIDTOF","; p_{T} [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPtBeautyReconstructedTracksPIDTOF);
    
    fPtBeautyReconstructedTracksPIDITS = new TH1F("fPtBeautyReconstructedTracksPIDITS","; p_{T} [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPtBeautyReconstructedTracksPIDITS);
    
    hTOFEffDen = new TH1F("hTOFEffDen","; p_{T} [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(hTOFEffDen);
    
    hTOFEffNum = new TH1F("hTOFEffNum","; p_{T} [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(hTOFEffNum);

    fPtBeautyPtrecVsPtparticle = new TH2F("fPtBeautyPtrecVsPtparticle","; p_{T} [GeV/c]; Count",32,ptbinning,32,ptbinning);
    fOutputList->Add(fPtBeautyPtrecVsPtparticle);
    
    
      
    ///THnSparse to store DCA of different particle species in MC-------------
    Int_t nBinspdg2 = 21;
    Double_t minpdg2 = 0.;
    Double_t maxpdg2 = 21.;
    Double_t binLimpdg2[nBinspdg2+1];
    for(Int_t i=0; i<=nBinspdg2; i++) binLimpdg2[i]=(Double_t)minpdg2 + (maxpdg2-minpdg2)/nBinspdg2*(Double_t)i ;
    
    Int_t nBinsdcaxy = 800; ///bin size 0.0005 cm
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
    /*
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
    */
    Int_t nBinstype = 9;
    Double_t mintype = -1.;
    Double_t maxtype = 8.;
    Double_t binLimtype[nBinstype+1];
    for(Int_t i=0; i<=nBinstype; i++) binLimtype[i]=(Double_t)mintype + (maxtype-mintype)/nBinstype*(Double_t)i ;
    
    Int_t nBinscharge = 4;
    Double_t mincharge = -2.;
    Double_t maxcharge = 2.;
    Double_t binLimcharge[nBinscharge+1];
    for(Int_t i=0; i<=nBinscharge; i++) binLimcharge[i]=(Double_t)mincharge + (maxcharge-mincharge)/nBinscharge*(Double_t)i ;
    
    /*
     Double_t binLimpt2[45]= {0.1,0.112797,0.127231,0.143512,0.161877,0.182592,0.205957,0.232313,0.262041,0.295573,0.333397,
     0.37606,0.424183,0.478465,0.539692,0.608754,0.686654,0.774523,0.873636,0.985432,1.11153,1.25377,1.41421,1.59519,1.79932,2.02957,
     2.28928,2.58223,2.91267,3.2854,3.70582,4.18004,4.71494,5.3183,5.99886,6.76651,7.6324,8.60909,9.71076,10.9534,12.3551,13.9361,15.7195,17.731,20};//bin limits from the measured pi0 spectrum
     */
    
    const Int_t nDima2=8;
    Int_t nBina2[nDima2] = {32,nBinspdg2,nBinsdcaxy,nBinsg,nBinsR,nBinstype,nBinsdcaxy,nBinscharge};
    fD0 = new THnSparseF("fD0","fD0",nDima2,nBina2);
    fD0->SetBinEdges(0,ptbinning); ///pt spectra -> same binning as other histograms
    fD0->SetBinEdges(1,binLimpdg2); /// electrons from D,charm baryons, B, beauty baryons, gamma, pi0, eta, Dcorrected, Dcorrected by weight, protons, kaons, D0_corr, D+-_corr,Ds_corr,Lc_corr, D0, D+-,Ds,Lc
    fD0->SetBinEdges(2,binLimdcaxy); ///dca distribution x charge x B
    fD0->SetBinEdges(3,binLimg);  ///From which generator (Hijing, else, pi0, eta)
    fD0->SetBinEdges(4,binLimR); ///Position where the electron is created
    //fD0->SetBinEdges(5,binLimITSchi2); ///ITS chi2 
    //fD0->SetBinEdges(6,binLimITSsha); ///fraction ITS shared clusters 
    fD0->SetBinEdges(5,binLimtype); ///pi0 and eta type  ///kNoMother, kNoFeedDown, kNoIsPrimary, kLightMesons, kBeauty, kCharm, kKaonFromHF, kKaonFromNonHF
    fD0->SetBinEdges(6,binLimdcaxy); /// dca distribution
    fD0->SetBinEdges(7,binLimcharge); /// track charge
    fD0->Sumw2();
    fOutputList->Add(fD0);
    ///-----------------------------------------------------------------
    
    ///THnSparse to store DCA in Data
    const Int_t nDima3=3;
    //Int_t nBina3[nDima3] = {32,nBinsdcaxy,nBinsITSchi2,nBinsITSsha,nBinspdg2};
    Int_t nBina3[nDima3] = {32,nBinsdcaxy,nBinspdg2};
    fD0Data = new THnSparseF("fD0Data","fD0Data",nDima3,nBina3);
    fD0Data->SetBinEdges(0,ptbinning); ///pt spectra -> same binning as other histograms
    fD0Data->SetBinEdges(1,binLimdcaxy); ///dca distribution x charge x B
    //fD0Data->SetBinEdges(2,binLimITSchi2); ///ITS chi2 
    //fD0Data->SetBinEdges(3,binLimITSsha); ///fraction ITS shared clusters 
    fD0Data->SetBinEdges(2,binLimpdg2); /// electrons and pions
    fD0Data->Sumw2();
    fOutputList->Add(fD0Data);
    
    PostData(1, fOutputList);
    
    
}

//////////////
//Main loop///
//////////////

///Called for each event
void AliAnalysisHFETPCTOFBeauty::UserExec(Option_t *)
{
	gRandom->SetSeed(0);
	//gRandom = new TRandom3(0);
    Int_t pdg = -99999;
    Int_t pdg_mother = -99999;
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
    Double_t qadca[10];
	Double_t qadcaData[10];
    for(int j=0;j<10;j++){
		qadca[j] = 0;
		qadcaData[j] = 0;
	}
    
    ///Check Event
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
    
    ///Check Cuts
    if(!fCuts)
    {
        AliError("HFE cuts not available");
        return;
    }
    ///Check PID
    if(!fPID->IsInitialized())
    {
        /// Initialize PID with the given run number
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
    
    ///PID response
    fPidResponse = fInputHandler->GetPIDResponse();
    
    
    ///Check PID response
    if(!fPidResponse)
    {
        AliDebug(1, "Using default PID Response");
        fPidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class());
    }
    
    fPID->SetPIDResponse(fPidResponse);
    
    fCFM->SetRecEventInfo(fVevent);
    
    ///Initialize V0 electron tagger
    if(fV0Tagger){
		fV0Tagger->Reset();
		fV0Tagger->TagV0Tracks(fAOD);
	}
	
	///Initialize ExtraCuts: for DCA calculation
	if(!fExtraCuts){
			fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
		}
	fExtraCuts->SetRecEventInfo(fAOD);
    
    
    Double_t *fListOfmotherkink = 0;
    Int_t fNumberOfVertices = 0;
    Int_t fNumberOfMotherkink = 0;
    
    
    
    ////////////////////
	//Vertex Selection//
	////////////////////
 
 	///Checking Nv0 before vertex selection (events triggered by v0)
    fNevent->Fill(0);

    ///Getting primary vertex    
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
	//Correlation cuts//
	////////////////////
	if(!PassCorrCuts(fAOD)) return;
        
    ///Number of events that passes the correlation cuts between SPD and track vertexes
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

       
     ///Checking Nv0 after the vertex selection (triggered by v0)
     fNevent_passvertex->Fill(0); 
        
        
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
        
        fCent2->Fill(cent);
    }
    
    ///Number of analysed events
    fNeventAnalized->Fill(0);
    
    ///Sign of the magnetic field
	//cout<<"fAOD->GetMagneticField() = "<<fAOD->GetMagneticField()<<endl;
    int signB = 0;
    if(fAOD->GetMagneticField() < 0) signB = -1;
    if(fAOD->GetMagneticField() > 0) signB = 1;
    
    int NAnalizedTracks = 0;   
    int NAnalizedTracksHijing = 0;         
        
        
    //=======================================================================    
    ///////////////////////////////////
	//Initialization for MC analysis///
	///////////////////////////////////
    if(fIsMC && fIsAOD)
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
            
        Bool_t test = GetNMCPartProduced(); ///Getting number of particles produced by the MC generator
           
        for(Int_t iMC = 0; iMC < fNTotMCpart; iMC++)
        {
			fMCparticle = (AliAODMCParticle*) fMCarray->At(iMC);
			///Pseudo-rapidity cut
			if((fMCparticle->Eta() < fEtaMin) || (fMCparticle->Eta() > fEtaMax)) continue;     				
			
            ///Get proportion between generated D0 and LambdaC:
            Int_t TrackPDG = TMath::Abs(fMCparticle->GetPdgCode());
            if(TrackPDG == 421) hPtD0->Fill(fMCparticle->Pt());///D0
            if(TrackPDG == 4122) hPtLambdaC->Fill(fMCparticle->Pt()); ///LambdaC
            if((TrackPDG>500 && TrackPDG<600)) fPtGeneratedBmesons->Fill(fMCparticle->Pt()); 
                    
            ///Beauty reconstruction efficiency block-----------
			if(TrackPDG == 11){	
				Bool_t MotherFound = FindMother(iMC);
				if(fIsFromMesonB || fIsFromBarionB || fIsFromBarionBD || fIsFromMesonBD){
					fPtBeautyGenerated->Fill(fMCparticle->Pt());
                    //cout<<iMC<<endl;
				}
			}
			///----------------------------------------------------  
        }
        //cout<<"fNTotMCpart = "<<fNTotMCpart<<endl; 
        //cout<<"fMCarray->GetEntries() = "<<fMCarray->GetEntries()<<endl; 
    }
    //cout<<"fVevent->GetNumberOfTracks() = "<<fVevent->GetNumberOfTracks()<<endl;
    //=======================================================================
    
    
    ///////////////
	//Track loop///
	///////////////
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
        
        ///Test Filter Bit
        if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; ///(Apply filter bit here, otherwise you count the track many times!!)

        ///////////////////////////
		//BEFORE TRACK SELECTION///
		///////////////////////////
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
        
        fTPCnsigma_TOFnsigma1->Fill(fTOFnSigma,fTPCnSigma);
        fTOFnsigma_p1->Fill(fP,fTOFnSigma);
        fTOFnsigma_pt1->Fill(fPt,fTOFnSigma);
        
        fITSnsigma_p1->Fill(fP,fITSnSigma);
        fITSnsigma_pt1->Fill(fPt,fITSnSigma);
        
        fPt_1->Fill(fPt);
        
        fITSnClus_1->Fill(fITSnClus);
        fTPCnClus_1->Fill(fTPCnClus);

        ///Pseudo-rapidity cut
        if((track->Eta() < fEtaMin) || (track->Eta() > fEtaMax)) continue;
        
        ///Beauty reconstruction efficiency block-----------
        if(fIsMC && fIsAOD){
			Bool_t IsHFEMC = IsHFelectronsMC(track);
			if(IsHFEMC){
				if(fIsFromMesonB || fIsFromBarionB || fIsFromBarionBD || fIsFromMesonBD){
					fPtBeautyReconstructedAll->Fill(fPt);
					//cout<<"reconstructed by track cut"<<endl;
				}
			}
		}
		///----------------------------------------------------
        
        
        //=======================================================================
        // Track Selection Cuts are applied here
        //=======================================================================        
                
        ///RecKine: ITS TPC cuts 
        if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
        
        ///RecKink 
        if(fRejectKinkMother)
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
            
        ///RecPrim
        if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;				///ProcessCutStep(Int_t cutStep, AliVParticle *track)
        
        ///HFEcuts: ITS layers cuts 
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;
        
        ///HFE cuts: TPC PID cleanup
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;
        
        ///ITS Chi2 
        Double_t ITSNcls = atrack->GetITSNcls();
			//cout<<"atrack->GetITSNcls() = "<<atrack->GetITSNcls()<<endl;
        if((atrack->GetITSchi2()/ITSNcls) > 5) continue; 
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
        if(fsharedclsITS > 0.3) continue;    
		
        
        
        
        ///////////////////////////
		//AFTER TRACK SELECTION////
		///////////////////////////
		
        ///Beauty reconstruction efficiency block-----------
        if(fIsMC && fIsAOD){
			Bool_t IsHFEMC = IsHFelectronsMC(track);
			if(IsHFEMC){
				if(fIsFromMesonB || fIsFromBarionB || fIsFromBarionBD || fIsFromMesonBD){
					fPtBeautyReconstructedTracks->Fill(fPt);
					//cout<<"reconstructed by track cut"<<endl;
				}
			}
		}
		///----------------------------------------------------
        
        ///Number of analised tracks
		NAnalizedTracks = NAnalizedTracks+1;
        
        ///Number of analysed tracks from Hijing:
        if(fIsMC){
			Int_t trkLabel = TMath::Abs(track->GetLabel());
			Int_t labelm = GetPrimary(trkLabel,fMCarray);///gives the label of first mother
			AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCarray->At(labelm);
			Int_t trkIndexPrimHFE= AODMCtrack->GetLabel();///gives index of the particle in original MCparticle array (labelm and trkIndexPrimHFE are the same, so I don't understand why this is done)
   
			if(trkIndexPrimHFE < fNpureMC){
				NAnalizedTracksHijing = NAnalizedTracksHijing+1;
			}
		}

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
        
        
        ////////////////////
		//Calculating DCA///
		////////////////////
		
		///The way it was done before
		//Double_t d0z0[2], cov[3];
		//AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
        //if(!(track->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 100., d0z0, cov))) continue;
        //Double_t DCAxy = d0z0[0];
        //Double_t DCAz = d0z0[1];
        //cout<<"---------------------- fPt = "<<fPt<<endl;
        //cout<<"---------------------- DCAxy = "<<DCAxy<<endl;
        
        
        ///Using AliHFEextraCuts method, that includes primary vertex recalculation if the number of contributors is smaller than 30       
        Double_t hfeDCA[2];
        hfeDCA[0] = -999.;
        hfeDCA[1] = -999.;
        Double_t hfeDCAResol[3];
        hfeDCAResol[0] = -999.;
        hfeDCAResol[1] = -999.;
        hfeDCAResol[2] = -999.;
		fExtraCuts->GetHFEImpactParameters((AliVTrack *)track, hfeDCA, hfeDCAResol);
		//Float_t hfeDCA = -999.;
        //Float_t hfeDCAResol = -999.;
		//fExtraCuts->GetImpactParameters((AliVTrack *)track, hfeDCA, hfeDCAResol);
        //cout<<"---------------------- hfeDCA = "<<hfeDCA[0]<<endl;
        //printf("----------------------- hfeDCA = %0.4f",hfeDCA);
        //printf("----------------------- hfeDCAResol = %0.4f",hfeDCAResol);
        Double_t DCAxy = hfeDCA[0];
        Double_t DCAz = hfeDCA[1];
  
  
        ///Checking nsigmaTPC after PID cuts in tof and its
        /*
        if(fTOFnSigma >= ftofPIDmincut && fTOFnSigma <= ftofPIDmaxcut){
            fTPCnsigma_p_after_tof->Fill(fP,fTPCnSigma);
            fTPCnsigma_pt_after_tof->Fill(fPt,fTPCnSigma);
            
			///Cheking the hadron nsigmaTPC after TOF cut
			if(fIsMC){
				fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
				if(TMath::Abs(fMCparticle->GetPdgCode()) == 2212)  fTPCnsigma_p_after_tof_p->Fill(fP,fTPCnSigma);
				if(TMath::Abs(fMCparticle->GetPdgCode()) == 130 || TMath::Abs(fMCparticle->GetPdgCode()) == 310 || TMath::Abs(fMCparticle->GetPdgCode()) == 311 || TMath::Abs(fMCparticle->GetPdgCode()) == 321)  fTPCnsigma_p_after_tof_k->Fill(fP,fTPCnSigma);
				if(TMath::Abs(fMCparticle->GetPdgCode()) == 211)  fTPCnsigma_p_after_tof_pion->Fill(fP,fTPCnSigma);
			}
			
            if(fITSnSigma >= -2 && fITSnSigma <= 2){
                fTPCnsigma_p_after_tof_its->Fill(fP,fTPCnSigma);
                fTPCnsigma_pt_after_tof_its->Fill(fPt,fTPCnSigma);
            }
        }
        */
        ///New PID cuts
        if(fP <= 2){
			if(fTOFnSigma >= -4 && fTOFnSigma <= 2){
				if(fITSnSigma >= -3 && fITSnSigma <= 3){
					fTPCnsigma_p_after_tof_v2->Fill(fP,fTPCnSigma);
				}
			}
		}
		if(fP > 2){
			if(fTOFnSigma >= ftofPIDmincut && fTOFnSigma <= ftofPIDmaxcut){
				fTPCnsigma_p_after_tof_v2->Fill(fP,fTPCnSigma);
			}
		}
        
        
        //////////////////////////////////////////////
		//Hadron DCA  --- for DCA resolution studies//
		//////////////////////////////////////////////
        if(fTPCnSigma >= -5 && fTPCnSigma <= -3){
			fDCAxy_pt_had_onlyDCA->Fill(fPt,DCAxy);
			fDCAxy_pt_had->Fill(fPt,DCAxy*track->Charge()*signB);
			fDCAz_pt_had->Fill(fPt,DCAz);	
			
			///Checking the effect of the improver in the resolution for hijing events separetely
			if(fIsMC){
				Int_t trkLabel = TMath::Abs(track->GetLabel());
				Int_t labelm = GetPrimary(trkLabel,fMCarray);///gives the label of first mother
				AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCarray->At(labelm);
				Int_t trkIndexPrimHFE= AODMCtrack->GetLabel();///gives index of the particle in original MCparticle array (labelm and trkIndexPrimHFE are the same, so I don't understand why this is done)
            
				///Particle generator:
				if(trkIndexPrimHFE < fNpureMC){
					fDCAxy_pt_had_onlyDCA_Hijing->Fill(fPt,DCAxy);
				}
				if(trkIndexPrimHFE >= fNpureMC){
					fDCAxy_pt_had_onlyDCA_Phytia->Fill(fPt,DCAxy);
				}			
			}
        }
        
        /////////////////////////////////////////////////////////////////////////////////////////////
		//THnSparse to store the DCA information of Data (electron candidates and charged hadrons)///
		/////////////////////////////////////////////////////////////////////////////////////////////
		qadcaData[0] = fPt;
         
		qadcaData[1] = DCAxy*track->Charge()*signB;
			
		qadcaData[2] = -1.;
			
		///Charged pions
		if(fTPCnSigma >= -5 && fTPCnSigma <= -3){
		  qadcaData[2] = 0.5;
		}
		///Electron candidates
		if(fTPCnSigma >= ftpcPIDmincut && fTPCnSigma <= ftpcPIDmaxcut){
			if(fP <= 2){
				if(fTOFnSigma >= -4 && fTOFnSigma <= 2){
					if(fITSnSigma >= -3 && fITSnSigma <= 3){
						qadcaData[2] = 1.5;	
					}
				}
			}
			if(fP > 2){
				if(fTOFnSigma >= ftofPIDmincut && fTOFnSigma <= ftofPIDmaxcut){
					qadcaData[2] = 1.5;					
				}
			}
		}

		if(qadcaData[2]>0.) fD0Data->Fill(qadcaData);
       
        
        /////////////////////////////////////////////////////////
		//Efficiencies of each PID cut for the beauty electrons//
		////////////////////////////////////////////////////////    
               
		///Checking PID cuts separately (for the efficiency)
		
		///TPC	
        if(fTPCnSigma >= ftpcPIDmincut && fTPCnSigma <= ftpcPIDmaxcut){
			///For the beauty reconstruction efficiency-----------
			if(fIsMC && fIsAOD){
				Bool_t IsHFEMC = IsHFelectronsMC(track);
				if(IsHFEMC){
					if(fIsFromMesonB || fIsFromBarionB || fIsFromBarionBD || fIsFromMesonBD){
						fPtBeautyReconstructedTracksPIDTPC->Fill(fPt);
					}
				}
			}
			///----------------------------------------------------
		}
        
        ///TOF low pt
        if(fP <= 2){
			if(fTOFnSigma >= -4 && fTOFnSigma <= 2){
				///For the beauty reconstruction efficiency-----------
				if(fIsMC && fIsAOD){
					Bool_t IsHFEMC = IsHFelectronsMC(track);
					if(IsHFEMC){
						if(fIsFromMesonB || fIsFromBarionB || fIsFromBarionBD || fIsFromMesonBD){
							fPtBeautyReconstructedTracksPIDTOF->Fill(fPt);
						}
					}
				}
			///----------------------------------------------------
			}
		}
		///TOF high pt
        if(fP > 2){
			if(fTOFnSigma >= ftofPIDmincut && fTOFnSigma <= ftofPIDmaxcut){
				///For the beauty reconstruction efficiency-----------
				if(fIsMC && fIsAOD){
					Bool_t IsHFEMC = IsHFelectronsMC(track);
					if(IsHFEMC){
						if(fIsFromMesonB || fIsFromBarionB || fIsFromBarionBD || fIsFromMesonBD){
							fPtBeautyReconstructedTracksPIDTOF->Fill(fPt);
						}
					}
				}
			///----------------------------------------------------
			}
		} 
		///ITS low pt
		if(fP <= 2){
			if(fITSnSigma >= -3 && fITSnSigma <= 3){
				///For the beauty reconstruction efficiency-----------
				if(fIsMC && fIsAOD){
					Bool_t IsHFEMC = IsHFelectronsMC(track);
					if(IsHFEMC){
						if(fIsFromMesonB || fIsFromBarionB || fIsFromBarionBD || fIsFromMesonBD){
							fPtBeautyReconstructedTracksPIDITS->Fill(fPt);
						}
					}
				}
			///----------------------------------------------------
			}
		}
		///ITS high pt
		if(fP > 2){
				///For the beauty reconstruction efficiency-----------
				if(fIsMC && fIsAOD){
					Bool_t IsHFEMC = IsHFelectronsMC(track);
					if(IsHFEMC){
						if(fIsFromMesonB || fIsFromBarionB || fIsFromBarionBD || fIsFromMesonBD){
							fPtBeautyReconstructedTracksPIDITS->Fill(fPt);
						}
					}
				}
			///----------------------------------------------------
		}         
        
        ///////////////////////////////////////////////////////
		//V0 electrons from systematic studies of TOF PID cut//
		/////////////////////////////////////////////////////// 
		
		AliPID::EParticleType myv0pid = fV0Tagger->GetV0Info(track->GetID()); /// enum EParticleType: kElectron = 0, kMuon = 1, kPion = 2, etc
		if(myv0pid == AliPID::kElectron){
			//cout<<"----------------------V0 electrons --------------------------"<<endl;
			fTPCnsigma_p_after_V0selection->Fill(fP,fTPCnSigma); ///Cross checking electron signal
			if(fTPCnSigma >= ftpcPIDmincut && fTPCnSigma <= ftpcPIDmaxcut){
				hTOFEffDen->Fill(fPt);
				if(fP > 2){
					if(fTOFnSigma >= ftofPIDmincut && fTOFnSigma <= ftofPIDmaxcut){
						hTOFEffNum->Fill(fPt);
					}
				}
				if(fP <= 2){
					if(fTOFnSigma >= -4 && fTOFnSigma <= 2){
						hTOFEffNum->Fill(fPt);
					}
				}
			}
		}
		
 
        //=======================================================================
        // Here the PID cuts are applied
        //=======================================================================        
        
        if(fTPCnSigma <= ftpcPIDmincut || fTPCnSigma >= ftpcPIDmaxcut) continue;
		if(fP <= 2){
			if(fTOFnSigma <= -4 || fTOFnSigma >= 2) continue;
			if(fITSnSigma <= -3 || fITSnSigma >= 3) continue;
		}
		if(fP > 2){
			if(fTOFnSigma <= ftofPIDmincut || fTOFnSigma >= ftofPIDmaxcut) continue;
		}
		//cout<< "fTOFnSigma = " <<fTOFnSigma<< "fP = "<< fP <<endl;

        /////////////////////////
		//AFTER PID SELECTION////
		/////////////////////////
        ///Beauty reconstruction efficiency block-----------
         if(fIsMC && fIsAOD){
			Bool_t IsHFEMC = IsHFelectronsMC(track);
			if(IsHFEMC){
				if(fIsFromMesonB || fIsFromBarionB || fIsFromBarionBD || fIsFromMesonBD){
					fPtBeautyReconstructedTracksPID->Fill(fPt);
					
					///Matrix for the unfolding
					fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
					fPtBeautyPtrecVsPtparticle->Fill(fMCparticle->Pt(),fPt);
					
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
              
        
        /////////////////////////////////////////////////
		//THnSparse to store the DCA information of MC///
		/////////////////////////////////////////////////
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
            
            ///Selecting particle
            qadca[1]=-1.; 
            
            ///------------
            if(TMath::Abs(fMCparticle->GetPdgCode()) == 2212) qadca[1]=10.5; ///to check DCA of protons
            if(TMath::Abs(fMCparticle->GetPdgCode()) == 321) qadca[1]=11.5; ///to check DCA of kaons
            ///------------

            /////////////////////////
			//Electrons from charm///
			/////////////////////////
            Bool_t IsHFEMC = IsHFelectronsMC(track);
            if(IsHFEMC){
                if(fIsFromD){
                    fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
                    fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                    pdg_mother = fMCparticleMother->GetPdgCode();
                    
                    if(TMath::Abs(pdg_mother)>400 && TMath::Abs(pdg_mother)<500){///charmed meson 
						qadca[1]=0.5; 
						if(TMath::Abs(pdg_mother) == 421) qadca[1]=16.5; ///to check DCA D0
						if(TMath::Abs(pdg_mother) == 411) qadca[1]=17.5; ///to check DCA D+-
						if(TMath::Abs(pdg_mother) == 433) qadca[1]=18.5; ///to check DCA Ds 
											 
						 hCharmMotherPt->Fill(fMCparticleMother->Pt());
						 hCharmMotherPt_vsElecPt->Fill(fPt,fMCparticleMother->Pt());
						 hElecPt_vsCharmMotherPt->Fill(fMCparticleMother->Pt(),fPt);
						 fDCAxy_pt_charmbef->Fill(fPt,DCAxy*track->Charge()*signB);
						 
						 ///Correcting pT spectrum
						 ///(Weight for each electron pT bin starting at one - charm correction)
						 float probAcceptD = -999;
						 if(fMCparticleMother->Pt() > 30) continue;
						 
						 if(fMCparticleMother->Pt() < fPt){ ///Accept all D mesons with pt smaller than the electrons pt
							 probAcceptD = 1;
						 }
						 else{
							if(fPt >=1.0 && fPt < 1.1) probAcceptD = fDcorr22->Eval(fMCparticleMother->Pt());
							if(fPt >=1.1 && fPt < 1.2) probAcceptD = fDcorr21->Eval(fMCparticleMother->Pt());
							if(fPt >=1.2 && fPt < 1.3) probAcceptD = fDcorr20->Eval(fMCparticleMother->Pt());
							if(fPt >=1.3 && fPt < 1.4) probAcceptD = fDcorr19->Eval(fMCparticleMother->Pt());
							if(fPt >=1.4 && fPt < 1.5) probAcceptD = fDcorr18->Eval(fMCparticleMother->Pt());
							if(fPt >=1.5 && fPt < 1.75) probAcceptD = fDcorr17->Eval(fMCparticleMother->Pt());
							if(fPt >=1.75 && fPt < 2) probAcceptD = fDcorr16->Eval(fMCparticleMother->Pt());
							if(fPt >=2 && fPt < 2.25) probAcceptD = fDcorr1->Eval(fMCparticleMother->Pt());
							if(fPt >=2.25 && fPt < 2.5) probAcceptD = fDcorr2->Eval(fMCparticleMother->Pt());
							if(fPt >=2.5 && fPt < 2.75) probAcceptD = fDcorr3->Eval(fMCparticleMother->Pt());
							if(fPt >=2.75 && fPt < 3) probAcceptD = fDcorr4->Eval(fMCparticleMother->Pt());
							if(fPt >=3 && fPt < 3.5) probAcceptD = fDcorr5->Eval(fMCparticleMother->Pt());
							if(fPt >=3.5 && fPt < 4) probAcceptD = fDcorr6->Eval(fMCparticleMother->Pt());
							if(fPt >=4 && fPt < 4.5) probAcceptD = fDcorr7->Eval(fMCparticleMother->Pt());
							if(fPt >=4.5 && fPt < 5) probAcceptD = fDcorr8->Eval(fMCparticleMother->Pt());
							if(fPt >=5 && fPt < 5.5) probAcceptD = fDcorr9->Eval(fMCparticleMother->Pt());
							if(fPt >=5.5 && fPt < 6) probAcceptD = fDcorr10->Eval(fMCparticleMother->Pt());
							if(fPt >=6 && fPt < 6.5) probAcceptD = fDcorr11->Eval(fMCparticleMother->Pt());
							if(fPt >=6.5 && fPt < 7) probAcceptD = fDcorr12->Eval(fMCparticleMother->Pt());
							if(fPt >=7 && fPt < 8) probAcceptD = fDcorr13->Eval(fMCparticleMother->Pt());
							if(fPt >=8 && fPt < 9) probAcceptD = fDcorr14->Eval(fMCparticleMother->Pt());
							if(fPt >=9 && fPt < 10) probAcceptD = fDcorr15->Eval(fMCparticleMother->Pt());
						}
						 
						 float a = gRandom->Uniform(0,1); ///Random number between 0 and 1
						 //cout<<"fPt = "<<fPt<<endl;
						 //cout<<"fMCparticleMother->Pt() = "<<fMCparticleMother->Pt()<<endl;
						 //cout<<"probAcceptD = "<<probAcceptD<<endl;
						 //cout<<"a = "<<a<<endl;
						 if(a < probAcceptD){
							qadca[1]=7.5;
							if(TMath::Abs(pdg_mother) == 421) qadca[1]=12.5; ///to check DCA D0 - corrected
							if(TMath::Abs(pdg_mother) == 411) qadca[1]=13.5; ///to check DCA D+- - corrected
							if(TMath::Abs(pdg_mother) == 433) qadca[1]=14.5; ///to check DCA Ds  - corrected	
							
							hCharmMotherPt_corr->Fill(fMCparticleMother->Pt());
							hCharmMotherPt_corr2->Fill(fMCparticleMother->Pt());
							hCharmMotherPt_vsElecPt_corr->Fill(fPt,fMCparticleMother->Pt());
							hElecPt_vsCharmMotherPt_corr->Fill(fMCparticleMother->Pt(),fPt);
							fDCAxy_pt_charmaft->Fill(fPt,DCAxy*track->Charge()*signB);						 
						 }												 
					}
					
					if(TMath::Abs(pdg_mother)>4000 && TMath::Abs(pdg_mother)<5000){///charmed baryon
						 qadca[1]=1.5;
						 if(TMath::Abs(pdg_mother) == 4122) qadca[1]=15.5; ///to check DCA Lc
					}
                }//end of electrons from charm
                
                //////////////////////////
				//Electrons from beauty///
				//////////////////////////
                if(fIsFromMesonB || fIsFromBarionB || fIsFromBarionBD || fIsFromMesonBD){
                    fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
                    fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                    pdg_mother = fMCparticleMother->GetPdgCode();
                    
                    if(fIsFromMesonB || fIsFromMesonBD){///beauty meson 
						qadca[1]=2.5;
						hBeautyMotherPt->Fill(fMCparticleMother->Pt(),fPt);
						fDCAxy_pt_beautybef->Fill(fPt,DCAxy*track->Charge()*signB);
						
						///Correcting the pT spectrum
						float probAcceptB = -999;
						if(fIsFromMesonB && (fMCparticleMother->Pt() < 30)){
							 probAcceptB = fBcorr->Eval(fMCparticleMother->Pt());///always evaluating the function in the pT of the B meson
							 //cout<<"pdg_mother = "<<pdg_mother<<endl;	
						}
						if(fIsFromMesonBD){
							if(fMCparticleMother->GetMother() > 0){
								fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
								float pdg_gmother = fMCparticleGMother->GetPdgCode();
								
								if(TMath::Abs(pdg_gmother)>500 && TMath::Abs(pdg_gmother)<600 && (fMCparticleGMother->Pt()<30)){
									 probAcceptB = fBcorr->Eval(fMCparticleGMother->Pt());
									 /*
									cout<<"---------------"<<endl;
									cout<<"pdg_mother BD = "<<pdg_mother<<endl;
									cout<<"pdg_gmother = "<<pdg_gmother<<endl;	
									cout<<"fMCparticleGMother->Pt() = "<<fMCparticleGMother->Pt()<<endl;		
									cout<<"probAcceptB = "<<probAcceptB<<endl;		
									cout<<"---------------"<<endl;
									*/
								}
								else if(fMCparticleGMother->GetMother() > 0){
									fMCparticleGGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleGMother->GetMother());
									float pdg_ggmother = fMCparticleGGMother->GetPdgCode();
									if(TMath::Abs(pdg_ggmother)>500 && TMath::Abs(pdg_ggmother)<600 && (fMCparticleGGMother->Pt()<30)) probAcceptB = fBcorr->Eval(fMCparticleGGMother->Pt());
									else if(fMCparticleGGMother->GetMother() > 0){
										fMCparticleGGGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleGGMother->GetMother());
										float pdg_gggmother = fMCparticleGGGMother->GetPdgCode();
										if (TMath::Abs(pdg_gggmother)>500 && TMath::Abs(pdg_gggmother)<600 && (fMCparticleGGGMother->Pt()<30)) probAcceptB = fBcorr->Eval(fMCparticleGGGMother->Pt());
									}
								}
							}
						}
						
						float b = gRandom->Uniform(0,1);
						//cout<<"fPt = "<<fPt<<endl;
						//cout<<"fMCparticleMother->Pt() = "<<fMCparticleMother->Pt()<<endl;
						//cout<<"probAcceptB = "<<probAcceptB<<endl;
						//cout<<"b = "<<b<<endl;
						
						if(b < probAcceptB){
							fDCAxy_pt_beautyaft->Fill(fPt,DCAxy*track->Charge()*signB);
							qadca[1]=19.5;
						}
					}
					
					if(fIsFromBarionB || fIsFromBarionBD){///beauty baryon
						qadca[1]=3.5;
					}
                }//end of electrons from beauty
            }//end of hfe
            
            
           ////////////////////////////////
		   //Non Heavy-Flavour electrons///
		   ////////////////////////////////
           if(!IsHFEMC){
                fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
                pdg = fMCparticle->GetPdgCode();
                ///Is electron:
                if(TMath::Abs(pdg) != 11) continue;
                if(fMCparticle->GetMother()<=0) continue;
                fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                pdg_mother = fMCparticleMother->GetPdgCode();
                
                ///Photonic Electrons:
                if(TMath::Abs(pdg_mother) == 111 || TMath::Abs(pdg_mother) == 221 || TMath::Abs(pdg_mother) == 22){
					if(TMath::Abs(pdg_mother) == 111) qadca[1]=4.5; 
					if(TMath::Abs(pdg_mother) == 221) qadca[1]=5.5;
					if(TMath::Abs(pdg_mother) == 22) qadca[1]=6.5;
					Int_t fType = GetPi0EtaType(fMCparticleMother,fMCarray);
					qadca[5]=fType;
                }
            }

     
            ///DCAxy x charge x B sign
            qadca[2]=DCAxy*track->Charge()*signB;
            
            ///DCAxy
            qadca[6]=DCAxy;
            
            ///track charge
            qadca[7]=track->Charge();
            
            /*
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
            */
                        
            ///Fill
            if(qadca[1]>0.) fD0->Fill(qadca);
        }
        
       

    }//End of track loop
    
    fNAnalizedTracks->Fill(NAnalizedTracks);
	fNAnalizedTracksHijing->Fill(NAnalizedTracksHijing);
	
       
    delete fListOfmotherkink;
    PostData(1, fOutputList);
    
}//end of main loop





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
//=======================================================================


//=======================================================================
//Setter for the PID cuts (TOF and TPC)
void AliAnalysisHFETPCTOFBeauty::SetPIDCuts(Float_t tpcPIDmincut, Float_t tpcPIDmaxcut, Float_t tofPIDmincut, Float_t tofPIDmaxcut) {
    ftpcPIDmincut = tpcPIDmincut;
    ftpcPIDmaxcut = tpcPIDmaxcut;
    ftofPIDmincut = tofPIDmincut;
    ftofPIDmaxcut = tofPIDmaxcut;
}
//=======================================================================


//=======================================================================
//Setter for the Eta cut
void AliAnalysisHFETPCTOFBeauty::SetEtaCut(Float_t EtaMin, Float_t EtaMax){
    fEtaMin = EtaMin;
    fEtaMax = EtaMax;
}
//=======================================================================


//=======================================================================
Bool_t AliAnalysisHFETPCTOFBeauty::FindMother(Int_t mcIndex)
{
    fIsHFE1 = kFALSE;
    fIsHFE2 = kFALSE;
    fIsNonHFE = kFALSE;
    fIsFromD = kFALSE;
    fIsFromBarionB = kFALSE;
	fIsFromMesonB = kFALSE;
    fIsFromBarionBD =kFALSE;
    fIsFromMesonBD = kFALSE;
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
            fIsFromBarionB = kFALSE;
			fIsFromMesonB = kFALSE;
			fIsFromBarionBD =kFALSE;
			fIsFromMesonBD = kFALSE;
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
            fIsFromBarionB = kFALSE;
			fIsFromMesonB = kFALSE;
			fIsFromBarionBD =kFALSE;
			fIsFromMesonBD = kFALSE;
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
            fIsFromBarionB = kFALSE;
			fIsFromMesonB = kFALSE;
			fIsFromBarionBD =kFALSE;
			fIsFromMesonBD = kFALSE;
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
            fIsFromBarionB = kFALSE;
			fIsFromMesonB = kFALSE;
			fIsFromBarionBD =kFALSE;
			fIsFromMesonBD = kFALSE;
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
    
    ///Tag Electron Source
    if(mpdg==111 || mpdg==221 || mpdg==22)
    {
        fIsHFE1 = kFALSE;
        fIsHFE2 = kFALSE;
        fIsNonHFE = kTRUE;
        fIsFromD = kFALSE;
        fIsFromBarionB = kFALSE;
		fIsFromMesonB = kFALSE;
        fIsFromBarionBD =kFALSE;
        fIsFromMesonBD = kFALSE;
        
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
        fIsHFE1 = kTRUE;
        
        fIsFromPi0 = kFALSE;
        fIsFromEta = kFALSE;
        fIsFromGamma = kFALSE;
        
        fIsNonHFE = kFALSE;
        
        fIsFromD = kFALSE;
        fIsFromBarionB = kFALSE;
        fIsFromMesonB = kFALSE;
        fIsFromBarionBD =kFALSE;
        fIsFromMesonBD = kFALSE;
        
        if((mpdg>400 && mpdg<500) || (mpdg>4000 && mpdg<5000)) //charmed mesons and baryons
        {
            if((gmpdg>500 && gmpdg<600) || (ggmpdg>500 && ggmpdg<600) || (gggmpdg>500 && gggmpdg<600)) //when the charm comes from beauty meson
            {
                fIsHFE1 = kTRUE;
                fIsFromD = kFALSE;
                fIsFromBarionB = kFALSE;
				fIsFromMesonB = kFALSE;
                fIsFromBarionBD = kFALSE;
                fIsFromMesonBD = kTRUE;
                return kTRUE;
            }
            else if((gmpdg>5000 && gmpdg<6000) || (ggmpdg>5000 && ggmpdg<6000) || (gggmpdg>5000 && gggmpdg<6000)) //when the charm comes from beauty barion
            {
                fIsHFE1 = kTRUE;
                fIsFromD = kFALSE;
                fIsFromBarionB = kFALSE;
				fIsFromMesonB = kFALSE;
                fIsFromBarionBD = kTRUE;
                fIsFromMesonBD = kFALSE;
                return kTRUE;
            }
            
            else
            {
                fIsHFE1 = kTRUE;
                fIsFromD = kTRUE;
                fIsFromBarionB = kFALSE;
				fIsFromMesonB = kFALSE;
                fIsFromBarionBD =kFALSE;
                fIsFromMesonBD = kFALSE;
                return kTRUE;
            }
        }
        else if((mpdg>500 && mpdg<600)) //beauty mesons 
        {
            fIsHFE1 = kTRUE;
            fIsFromD = kFALSE;
            fIsFromBarionB = kFALSE;
			fIsFromMesonB = kTRUE;
            fIsFromBarionBD =kFALSE;
            fIsFromMesonBD = kFALSE;
            return kTRUE;
        }
        else if((mpdg>5000 && mpdg<6000)) //beauty baryons
        {
            fIsHFE1 = kTRUE;
            fIsFromD = kFALSE;
            fIsFromBarionB = kTRUE;
			fIsFromMesonB = kFALSE;
            fIsFromBarionBD =kFALSE;
            fIsFromMesonBD = kFALSE;
            return kTRUE;
        }
        
        else
        {
            fIsHFE1 = kFALSE;
            fIsFromD = kFALSE;
            fIsFromBarionB = kFALSE;
			fIsFromMesonB = kFALSE;
            fIsFromBarionBD =kFALSE;
            fIsFromMesonBD = kFALSE;
            return kFALSE;
        }
    }
}
//======================================================================


//=======================================================================

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
//=======================================================================


//=======================================================================
Int_t AliAnalysisHFETPCTOFBeauty::GetPi0EtaType(AliAODMCParticle *pi0eta, TClonesArray *mcArray){
    //
    // Return what type of pi0, eta it is
    //
    
    //IsPrimary
    Bool_t primMC = pi0eta->IsPrimary(); ///Does not include particles from weak decays or created in an interaction with the material
    if(!primMC) return kNoIsPrimary;
    
    // Mother
    Int_t motherlabel = pi0eta->GetMother();
    if(motherlabel<0) return kNoMother;
    else {
        
        AliAODMCParticle *mother = (AliAODMCParticle*)fMCarray->At(motherlabel);
        Int_t motherpdg = TMath::Abs(mother->GetPdgCode());
        ///pi0, eta, omega, phi, eta',rho0, rho+
        if(motherpdg == 111 || motherpdg == 221 || motherpdg == 223 || motherpdg == 333 || motherpdg == 331 || motherpdg == 113 || motherpdg == 213) return kLightMesons;
        
        ///If the mother is kaons from heavy-flavour decay (K0L,K0S,K0,K+,,k*0,,K*+)
        if(motherpdg == 130 || motherpdg == 310 || motherpdg == 311 || motherpdg == 321 || motherpdg == 313 || motherpdg == 323){
			Int_t kaonmotherlabel = mother->GetMother();
			if(kaonmotherlabel>0){
				AliAODMCParticle *kaonmother = (AliAODMCParticle*)fMCarray->At(kaonmotherlabel);
				Int_t kaonmotherpdg = TMath::Abs(kaonmother->GetPdgCode());
				if ( (int(TMath::Abs(kaonmotherpdg)/100.)%10) == 5 || (int(TMath::Abs(kaonmotherpdg)/1000.)%10) == 5 || (int(TMath::Abs(kaonmotherpdg)/100.)%10) == 4 || (int(TMath::Abs(kaonmotherpdg)/1000.)%10) == 4 ) return kKaonFromHF;
				else return kKaonFromNonHF;
			}
			else return kKaonFromNonHF;
        }
        
        if ( (int(TMath::Abs(motherpdg)/100.)%10) == 5 || (int(TMath::Abs(motherpdg)/1000.)%10) == 5 ) return kBeauty; ///beauty mesons and barions
        if ( (int(TMath::Abs(motherpdg)/100.)%10) == 4 || (int(TMath::Abs(motherpdg)/1000.)%10) == 4 ) return kCharm; ///charmed mesons and barions
        return kNoFeedDown;    
        
    }
}
//=======================================================================


//=======================================================================

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
//=======================================================================


//=======================================================================
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
//=======================================================================


//=======================================================================
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
//=======================================================================


//=======================================================================
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
//=======================================================================





