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
//      Task for Beauty analysis in p-p collisions   				  //
//      															  //
//																	  //
//		v1.0														  //
//                                                                    //
//	    Authors 							                          //
//		Sudhir Pandurang Rode (sudhir.pandurang.rode@cern.ch)				      //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "AliAnalysisUtils.h"
#include "TChain.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
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
#include "AliHFEextraCuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliSelectNonHFE.h"
#include "AliHFEpidTPC.h"
#include "AliAnalysisHFEppTPCTOFBeauty.h"
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
#include "AliHFEmcQA.h"
#include "AliHFEsignalCuts.h"
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
//______________________________________________________________________

//______________________________________________________________________
ClassImp(AliAnalysisHFEppTPCTOFBeauty)

//______________________________________________________________________
AliAnalysisHFEppTPCTOFBeauty::AliAnalysisHFEppTPCTOFBeauty(const char *name)
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
,fSignalCuts(NULL)
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
//,fNAnalizedTracksHijing(0)
//,fITSnClus_1(0)
,fTPCnClus_1(0)
//,fITSnClus_2(0)
,fTPCnClus_2(0)
,fTPCnsigma_p1(0)
,fTPCnsigma_p2(0)
,fTPCnsigma_p3(0)
,fTPCnsigma_pt1(0)
,fTPCnsigma_pt2(0)
,fTPCnsigma_pt3(0)
,fTPCnsigma_p_after_tof(0)
,fTPCnsigma_proton_p_after_tof(0)
,fTPCnsigma_p_after_tof_p(0)
,fTPCnsigma_p_after_tof_pion(0)
,fTPCnsigma_p_after_tof_k(0)
,fTPCnsigma_pt_after_tof(0)
//,fTPCnsigma_p_after_tof_its(0)
//,fTPCnsigma_pt_after_tof_its(0)
,fTOFnsigma_p1(0)
,fTOFnsigma_p2(0)
,fTOFnsigma_p3(0)
,fTOFnsigma_pt1(0)
,fTOFnsigma_pt2(0)
,fTOFnsigma_pt3(0)
//,fITSnsigma_p1(0)
//,fITSnsigma_p2(0)
//,fITSnsigma_p3(0)
//,fITSnsigma_pt1(0)
//,fITSnsigma_pt2(0)
//,fITSnsigma_pt3(0)
,fTPCnsigma_TOFnsigma1(0)
,fTPCnsigma_TOFnsigma2(0)
,fTPCnsigma_TOFnsigma3(0)
,fDCAxy_pt_had(0)
,fDCAxy_pt_had_WoPID(0)
,fDCAxy_pt_charmbef(0)
,fDCAxy_pt_charmaft(0)
,fDCAxy_pt_beautybef(0)
,fDCAxy_pt_beautyaft(0)
,fDCAxy_pt_MesonB_beautyaft(0)
,fDCAxy_pt_MesonB_beautybef(0)
,fDCAxy_pt_MesonBD_beautyaft(0)
,fDCAxy_pt_MesonBD_beautybef(0)
,fDCAxy_pt_BaryonB_beautybef(0)
,fDCAxy_pt_BaryonBD_beautybef(0)
,fDCAxy_pt_had_onlyDCA(0)
,fDCAxy_pt_had_onlyDCA_WoPID(0)
//,fDCAxy_pt_had_onlyDCA_Hijing(0)
,fDCAxy_pt_had_onlyDCA_Phytia(0)
,fDCAz_pt_had(0)
,fDCAz_pt_had_WoPID(0)
,fDCAxy_pt_had_onlyDCA_phi1(0)
,fDCAxy_pt_had_onlyDCA_phi2(0)
,fDCAxy_pt_had_onlyDCA_phi3(0)
,fDCAxy_pt_had_onlyDCA_phi4(0)
,fDCAxy_pt_had_phi1_ChB(0)
,fDCAxy_pt_had_phi1_B(0)
,fDCAxy_pt_had_phi2_ChB(0)
,fDCAxy_pt_had_phi2_B(0)
,fDCAxy_pt_had_phi3_ChB(0)
,fDCAxy_pt_had_phi3_B(0)
,fDCAxy_pt_had_phi4_ChB(0)
,fDCAxy_pt_had_phi4_B(0)
,fDCAxy_pt_had_ResCorr_phi1(0)
,fDCAxy_pt_had_ResCorr_phi2(0)
,fDCAxy_pt_had_ResCorr_phi3(0)
,fDCAxy_pt_had_ResCorr_phi4(0)
,fResGausCorr_phi1(0)
,fResGausCorr_phi2(0)
,fResGausCorr_phi3(0)
,fResGausCorr_phi4(0)
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
/*,hDCAPtProtons(0)
,hDCAPtProtons2(0)
,hDCAPtProtons3(0)*/
,hBeautyMotherPtbef(0)
,hBeautyMotherPtaft(0)
,hBeautyMotherPt2Daft(0)
,fPtBeautyGenerated(0)
,fPtGeneratedBmesons(0)
,fPtBeautyReconstructedAll(0)
,fPtBeautyReconstructedTracks(0)
,fPtBeautyReconstructedTracksPID(0)
,fPtBeautyReconstructedTracksPIDTPC(0)
,fPtBeautyReconstructedTracksPIDTOF(0)
,hPtLambdaC(0)
,hPtD0(0)
,hPtDstar(0)
,hPtDp(0)
,hPtDs(0)
,fLcD01(0)
,hEleDstar_WCorr_WFCorr(0)
,hEleDp_WCorr_WFCorr(0)
,hEleDs_WCorr_WFCorr(0)
,hEleLc_WCorr_WFCorr(0)
,hEleD0_WCorr(0)
,HistWeightLcD0(0)
,HistWeightDstarD0(0)
,HistWeightDsD0(0)
,HistWeightDpD0(0)
//For the HFE package
,fCuts(0)
,fCFM(0)
,fPID(new AliHFEpid("hfePid"))
,fPIDqa(0)

//For MC
,fRejectKinkMother(kTRUE)
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
//,fD0HC(0)
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
,probAcceptB(0)
,probAcceptD(0)
,fPionsPt(0)
,fKaonsPt(0)
//,fHC(0)
/*,fDCAxy_pt_Dalitz(0)
,fDCAxy_pt_DalitzFromFeeddown(0)
,fDCAxy_pt_Conversions(0)
,fDCAxy_pt_ConversionsFromFeeddown(0)
,fDCAxy_pt_ConversionsFromStrangeFeeddown(0)*/
,fDCAxy_pt_Dalitz2(0)
//,fDCAxy_pt_DalitzFromFeeddown2(0)
,fDCAxy_pt_Conversions2(0)
,fDCAxy_pt_Beauty2(0)
,fDCAxy_pt_Charm2(0)
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
AliAnalysisHFEppTPCTOFBeauty::AliAnalysisHFEppTPCTOFBeauty()
: AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisHFEppTPCTOFBeauty")

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
,fSignalCuts(NULL)
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
//,fNAnalizedTracksHijing(0)
//,fITSnClus_1(0)
,fTPCnClus_1(0)
//,fITSnClus_2(0)
,fTPCnClus_2(0)
,fTPCnsigma_p1(0)
,fTPCnsigma_p2(0)
,fTPCnsigma_p3(0)
,fTPCnsigma_pt1(0)
,fTPCnsigma_pt2(0)
,fTPCnsigma_pt3(0)
,fTPCnsigma_p_after_tof(0)
,fTPCnsigma_proton_p_after_tof(0)
,fTPCnsigma_p_after_tof_p(0)
,fTPCnsigma_p_after_tof_pion(0)
,fTPCnsigma_p_after_tof_k(0)
,fTPCnsigma_pt_after_tof(0)
//,fTPCnsigma_p_after_tof_its(0)
//,fTPCnsigma_pt_after_tof_its(0)
,fTOFnsigma_p1(0)
,fTOFnsigma_p2(0)
,fTOFnsigma_p3(0)
,fTOFnsigma_pt1(0)
,fTOFnsigma_pt2(0)
,fTOFnsigma_pt3(0)
//,fITSnsigma_p1(0)
//,fITSnsigma_p2(0)
//,fITSnsigma_p3(0)
//,fITSnsigma_pt1(0)
//,fITSnsigma_pt2(0)
//,fITSnsigma_pt3(0)
,fTPCnsigma_TOFnsigma1(0)
,fTPCnsigma_TOFnsigma2(0)
,fTPCnsigma_TOFnsigma3(0)
,fDCAxy_pt_had(0)
,fDCAxy_pt_had_WoPID(0)
,fDCAxy_pt_charmbef(0)
,fDCAxy_pt_charmaft(0)
,fDCAxy_pt_beautybef(0)
,fDCAxy_pt_beautyaft(0)
,fDCAxy_pt_MesonB_beautyaft(0)
,fDCAxy_pt_MesonB_beautybef(0)
,fDCAxy_pt_MesonBD_beautyaft(0)
,fDCAxy_pt_MesonBD_beautybef(0)
,fDCAxy_pt_BaryonB_beautybef(0)
,fDCAxy_pt_BaryonBD_beautybef(0)
,fDCAxy_pt_had_onlyDCA(0)
,fDCAxy_pt_had_onlyDCA_WoPID(0)
//,fDCAxy_pt_had_onlyDCA_Hijing(0)
,fDCAxy_pt_had_onlyDCA_Phytia(0)
,fDCAz_pt_had(0)
,fDCAz_pt_had_WoPID(0)
,fDCAxy_pt_had_onlyDCA_phi1(0)
,fDCAxy_pt_had_onlyDCA_phi2(0)
,fDCAxy_pt_had_onlyDCA_phi3(0)
,fDCAxy_pt_had_onlyDCA_phi4(0)
,fDCAxy_pt_had_phi1_ChB(0)
,fDCAxy_pt_had_phi1_B(0)
,fDCAxy_pt_had_phi2_ChB(0)
,fDCAxy_pt_had_phi2_B(0)
,fDCAxy_pt_had_phi3_ChB(0)
,fDCAxy_pt_had_phi3_B(0)
,fDCAxy_pt_had_phi4_ChB(0)
,fDCAxy_pt_had_phi4_B(0)
,fDCAxy_pt_had_ResCorr_phi1(0)
,fDCAxy_pt_had_ResCorr_phi2(0)
,fDCAxy_pt_had_ResCorr_phi3(0)
,fDCAxy_pt_had_ResCorr_phi4(0)
,fResGausCorr_phi1(0)
,fResGausCorr_phi2(0)
,fResGausCorr_phi3(0)
,fResGausCorr_phi4(0)
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
/*,hDCAPtProtons(0)
,hDCAPtProtons2(0)
,hDCAPtProtons3(0)*/
,hBeautyMotherPtbef(0)
,hBeautyMotherPtaft(0)
,hBeautyMotherPt2Daft(0)
,fPtBeautyGenerated(0)
,fPtGeneratedBmesons(0)
,fPtBeautyReconstructedAll(0)
,fPtBeautyReconstructedTracks(0)
,fPtBeautyReconstructedTracksPID(0)
,fPtBeautyReconstructedTracksPIDTPC(0)
,fPtBeautyReconstructedTracksPIDTOF(0)
,hPtLambdaC(0)
,hPtD0(0)
,hPtDstar(0)
,hPtDp(0)
,hPtDs(0)
,fLcD01(0)
,hEleDstar_WCorr_WFCorr(0)
,hEleDp_WCorr_WFCorr(0)
,hEleDs_WCorr_WFCorr(0)
,hEleLc_WCorr_WFCorr(0)
,hEleD0_WCorr(0)
,HistWeightLcD0(0)
,HistWeightDstarD0(0)
,HistWeightDsD0(0)
,HistWeightDpD0(0)
//For the HFE package
,fCuts(0)
,fCFM(0)
,fPID(new AliHFEpid("hfePid"))
,fPIDqa(0)

//For MC
,fRejectKinkMother(kTRUE)
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
//,fD0HC(0)
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
,probAcceptB(0)
,probAcceptD(0)
,fPionsPt(0)
,fKaonsPt(0)
//,fHC(0)
/*,fDCAxy_pt_Dalitz(0)
,fDCAxy_pt_DalitzFromFeeddown(0)
,fDCAxy_pt_Conversions(0)
,fDCAxy_pt_ConversionsFromFeeddown(0)
,fDCAxy_pt_ConversionsFromStrangeFeeddown(0)*/
,fDCAxy_pt_Dalitz2(0)
//,fDCAxy_pt_DalitzFromFeeddown2(0)
,fDCAxy_pt_Conversions2(0)
,fDCAxy_pt_Beauty2(0)
,fDCAxy_pt_Charm2(0)
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
AliAnalysisHFEppTPCTOFBeauty::~AliAnalysisHFEppTPCTOFBeauty()
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
void AliAnalysisHFEppTPCTOFBeauty::UserCreateOutputObjects()
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
    
    Double_t ptbinningHF[15] = {0,1,2,3,4,5,6,7,8,10,12,16,24,36,50};
    
    Double_t ptbinningHF2[15] = {0,1,2,3,4,5,6,7,8,10,12,16,24,36,50};
    
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
    
    fPionsPt = new TH1F("fPionsPt","pt (GeV/c)",1000,0,100);
    fOutputList->Add(fPionsPt);
        
    fKaonsPt = new TH1F("fKaonsPt","pt (GeV/c)",1000,0,100);
    fOutputList->Add(fKaonsPt);


    //fNAnalizedTracksHijing = new TH1F("fNAnalizedTracksHijing","pt (GeV/c)",5000,0,5000);
    //fOutputList->Add(fNAnalizedTracksHijing);
    /*
    fITSnClus_1 = new TH1F("fITSnClus_1","fITSnClus_1",1000,0,80);
    fOutputList->Add(fITSnClus_1);
    
    fITSnClus_2 = new TH1F("fITSnClus_2","fITSnClus_2",1000,0,80);
    fOutputList->Add(fITSnClus_2);
    */
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
    
   /* fITSnsigma_p1 = new TH2F("fITSnsigma_p1","p (GeV/c);ITS Electron N#sigma",300,0,15,200,-10,30);
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
    */
    fTPCnsigma_TOFnsigma1 = new TH2F("fTPCnsigma_TOFnsigma1","TOF Electron N#sigma;TPC Electron N#sigma",200,-10,30,200,-15,10);
    fOutputList->Add(fTPCnsigma_TOFnsigma1);
    
    fTPCnsigma_TOFnsigma2 = new TH2F("fTPCnsigma_TOFnsigma2","TOF Electron N#sigma;TPC Electron N#sigma",200,-10,30,200,-15,10);
    fOutputList->Add(fTPCnsigma_TOFnsigma2);
    
    fTPCnsigma_TOFnsigma3 = new TH2F("fTPCnsigma_TOFnsigma3","TOF Electron N#sigma;TPC Electron N#sigma",200,-10,30,200,-15,10);
    fOutputList->Add(fTPCnsigma_TOFnsigma3);
    
    hCharmMotherPt = new TH1F("hCharmMotherPt","; p_{T} [GeV/c]; Count",14,ptbinningHF);
    fOutputList->Add(hCharmMotherPt);
    
    hCharmMotherPt_vsElecPt = new TH2F("hCharmMotherPt_vsElecPt","; p_{T} [GeV/c]; Count",1000,0,50,1000,0,50);
    fOutputList->Add(hCharmMotherPt_vsElecPt);
    
    hElecPt_vsCharmMotherPt = new TH2F("hElecPt_vsCharmMotherPt","; p_{T} [GeV/c]; Count",1000,0,50,1000,0,50);
    fOutputList->Add(hElecPt_vsCharmMotherPt);
    
    hCharmMotherPt_vsElecPt_corr = new TH2F("hCharmMotherPt_vsElecPt_corr","; p_{T} [GeV/c]; Count",1000,0,50,1000,0,50);
    fOutputList->Add(hCharmMotherPt_vsElecPt_corr);
    
    hElecPt_vsCharmMotherPt_corr = new TH2F("hElecPt_vsCharmMotherPt_corr","; p_{T} [GeV/c]; Count",1000,0,50,1000,0,50);
    fOutputList->Add(hElecPt_vsCharmMotherPt_corr);
    
    hCharmMotherPt_corr = new TH1F("hCharmMotherPt_corr","; p_{T} [GeV/c]; Count",14,ptbinningHF);
    fOutputList->Add(hCharmMotherPt_corr);
    
    hCharmMotherPt_corr2 = new TH1F("hCharmMotherPt_corr2","; p_{T} [GeV/c]; Count",100,0,50);
    fOutputList->Add(hCharmMotherPt_corr2);
    
    hBeautyMotherPtbef = new TH1F("hBeautyMotherPtbef","; p_{T} [GeV/c]; Count",14,ptbinningHF2);
    fOutputList->Add(hBeautyMotherPtbef);
    
    hBeautyMotherPtaft = new TH1F("hBeautyMotherPtaft","; p_{T} [GeV/c]; Count",14,ptbinningHF2);
    fOutputList->Add(hBeautyMotherPtaft);
    
    hBeautyMotherPt = new TH2F("hBeautyMotherPt","; p_{T} [GeV/c]; Count",1000,0,50,1000,0,50);
    fOutputList->Add(hBeautyMotherPt);
    
    /*hDCAPtProtons = new TH2F("hDCAPtProtons","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
    fOutputList->Add(hDCAPtProtons);
    
    hDCAPtProtons2 = new TH2F("hDCAPtProtons2","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
    fOutputList->Add(hDCAPtProtons2);
    
    hDCAPtProtons3 = new TH2F("hDCAPtProtons3","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
    fOutputList->Add(hDCAPtProtons3);
    */
    hBeautyMotherPt2Daft = new TH2F("hBeautyMotherPt2Daft","; p_{T} [GeV/c]; Count",1000,0,50,1000,0,50);
    fOutputList->Add(hBeautyMotherPt2Daft);

    fPtElec = new TH1F("fPtElec","; p_{T} [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPtElec);
    
    hPtD0 = new TH1F("hPtD0","; p_{T} [GeV/c]; Count",1000,0, 50);
    fOutputList->Add(hPtD0);
    
    hPtLambdaC = new TH1F("hPtLambdaC","; p_{T} [GeV/c]; Count",1000,0, 50);
    fOutputList->Add(hPtLambdaC);
    
    hPtDp = new TH1F("hPtDp","; p_{T} [GeV/c]; Count",1000,0, 50);
    fOutputList->Add(hPtDp);
    
    hPtDstar = new TH1F("hPtDstar","; p_{T} [GeV/c]; Count",1000,0, 50);
    fOutputList->Add(hPtDstar);
    
    hPtDs = new TH1F("hPtDs","; p_{T} [GeV/c]; Count",1000,0, 50);
    fOutputList->Add(hPtDs);
    
        
    fPElec = new TH1F("fPElec","; p [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPElec);

    fPtHad_f = new TH1F("fPtHad_f","; p_{T} [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPtHad_f);
    fPtHad_f->Sumw2();
    
    fPHad_f = new TH1F("fPHad_f","; p [GeV/c]; Count",32,ptbinning);
    fOutputList->Add(fPHad_f);
    fPHad_f->Sumw2();
       
    fDCAxy_pt_had = new TH2F("fDCAxy_pt_had",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had);
    
    fDCAxy_pt_had_WoPID = new TH2F("fDCAxy_pt_had_WoPID",";p_{t} (GeV/c);DCAxy hadrons_WoPID",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_WoPID);
    
    fDCAxy_pt_charmbef = new TH2F("fDCAxy_pt_charmbef",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_charmbef);
    
    fDCAxy_pt_charmaft = new TH2F("fDCAxy_pt_charmaft",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_charmaft);
    
    fDCAxy_pt_beautybef = new TH2F("fDCAxy_pt_beautybef",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_beautybef);
    
    fDCAxy_pt_beautyaft = new TH2F("fDCAxy_pt_beautyaft",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_beautyaft); 
    
    fDCAxy_pt_MesonB_beautybef = new TH2F("fDCAxy_pt_MesonB_beautybef",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_MesonB_beautybef);
    
    fDCAxy_pt_MesonB_beautyaft = new TH2F("fDCAxy_pt_MesonB_beautyaft",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_MesonB_beautyaft);
    
    fDCAxy_pt_MesonBD_beautybef = new TH2F("fDCAxy_pt_MesonBD_beautybef",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_MesonBD_beautybef);
    
    fDCAxy_pt_MesonBD_beautyaft = new TH2F("fDCAxy_pt_MesonBD_beautyaft",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_MesonBD_beautyaft);
    
    fDCAxy_pt_BaryonB_beautybef = new TH2F("fDCAxy_pt_BaryonB_beautybef",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_BaryonB_beautybef);
     
    fDCAxy_pt_BaryonBD_beautybef = new TH2F("fDCAxy_pt_BaryonBD_beautybef",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_BaryonBD_beautybef);
       
    fDCAxy_pt_had_onlyDCA_WoPID = new TH2F("fDCAxy_pt_had_onlyDCA_WoPID",";p_{t} (GeV/c);DCAxy hadrons_WoPID",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_onlyDCA_WoPID);
    
    fDCAxy_pt_had_onlyDCA = new TH2F("fDCAxy_pt_had_onlyDCA",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_onlyDCA);
    /*
    fDCAxy_pt_had_onlyDCA_Hijing = new TH2F("fDCAxy_pt_had_onlyDCA_Hijing",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_onlyDCA_Hijing);
       */
    fDCAxy_pt_had_onlyDCA_Phytia = new TH2F("fDCAxy_pt_had_onlyDCA_Phytia",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_onlyDCA_Phytia);
       
    fDCAz_pt_had = new TH2F("fDCAz_pt_had",";p_{t} (GeV/c);DCAz hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAz_pt_had);
    
    fDCAxy_pt_had_onlyDCA_phi1 = new TH2F("fDCAxy_pt_had_onlyDCA_phi1",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_onlyDCA_phi1);
    
    fDCAxy_pt_had_phi1_ChB = new TH2F("fDCAxy_pt_had_phi1_ChB",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_phi1_ChB);
    
    fDCAxy_pt_had_phi1_B = new TH2F("fDCAxy_pt_had_phi1_B",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_phi1_B);
    
    fDCAxy_pt_had_onlyDCA_phi2 = new TH2F("fDCAxy_pt_had_onlyDCA_phi2",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_onlyDCA_phi2);
    
    fDCAxy_pt_had_phi2_ChB = new TH2F("fDCAxy_pt_had_phi2_ChB",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_phi2_ChB);
    
    fDCAxy_pt_had_phi2_B = new TH2F("fDCAxy_pt_had_phi2_B",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_phi2_B);
    
    fDCAxy_pt_had_onlyDCA_phi3 = new TH2F("fDCAxy_pt_had_onlyDCA_phi3",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_onlyDCA_phi3);
    
    fDCAxy_pt_had_phi3_ChB = new TH2F("fDCAxy_pt_had_phi3_ChB",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_phi3_ChB);
    
    fDCAxy_pt_had_phi3_B = new TH2F("fDCAxy_pt_had_phi3_B",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_phi3_B);
    
    fDCAxy_pt_had_onlyDCA_phi4 = new TH2F("fDCAxy_pt_had_onlyDCA_phi4",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_onlyDCA_phi4);
    
    fDCAxy_pt_had_phi4_ChB = new TH2F("fDCAxy_pt_had_phi4_ChB",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_phi4_ChB);
    
    fDCAxy_pt_had_phi4_B = new TH2F("fDCAxy_pt_had_phi4_B",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_phi4_B);
    
       
    fDCAxy_pt_had_ResCorr_phi1 = new TH2F("fDCAxy_pt_had_ResCorr_phi1",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_ResCorr_phi1);
    
    fDCAxy_pt_had_ResCorr_phi2 = new TH2F("fDCAxy_pt_had_ResCorr_phi2",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_ResCorr_phi2);
    
    fDCAxy_pt_had_ResCorr_phi3 = new TH2F("fDCAxy_pt_had_ResCorr_phi3",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_ResCorr_phi3);
    
    fDCAxy_pt_had_ResCorr_phi4 = new TH2F("fDCAxy_pt_had_ResCorr_phi4",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_ResCorr_phi4);
    
    
    fDCAz_pt_had_WoPID = new TH2F("fDCAz_pt_had_WoPID",";p_{t} (GeV/c);DCAz hadrons_WoPID",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAz_pt_had_WoPID);
    
    fDCAxy_pt_ele = new TH2F("fDCAxy_pt_ele",";p_{t} (GeV/c);DCAxy ",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele);
    
    fDCAz_pt_ele = new TH2F("fDCAz_pt_ele",";p_{t} (GeV/c);DCAz ",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAz_pt_ele); 
    
    fPtMCeta = new TH1F("fPtMCeta",";p_{t} (GeV/c)",2000,0,100);
    fOutputList->Add(fPtMCeta);
    
    fTPCnsigma_p_after_tof = new TH2F("fTPCnsigma_p_after_tof","p (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_p_after_tof);
    
    fTPCnsigma_p_after_tof_p = new TH2F("fTPCnsigma_p_after_tof_p","p (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_p_after_tof_p);
    
    fTPCnsigma_proton_p_after_tof = new TH2F("fTPCnsigma_proton_p_after_tof","p (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_proton_p_after_tof);
    
    fTPCnsigma_p_after_tof_pion = new TH2F("fTPCnsigma_p_after_tof_pion","p (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_p_after_tof_pion);
    
    fTPCnsigma_p_after_tof_k = new TH2F("fTPCnsigma_p_after_tof_k","p (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_p_after_tof_k);
       
   /* fTPCnsigma_p_after_tof_its = new TH2F("fTPCnsigma_p_after_tof_its","p (GeV/c);TPC Electron N#sigma after TOF and ITS cuts",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_p_after_tof_its);
    */
    
    fTPCnsigma_pt_after_tof = new TH2F("fTPCnsigma_pt_after_tof","pt (GeV/c);TPC Electron N#sigma after TOF cut",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_pt_after_tof);
    
   /* fTPCnsigma_pt_after_tof_its = new TH2F("fTPCnsigma_pt_after_tof_its","pt (GeV/c);TPC Electron N#sigma after TOF and ITS cuts",300,0,15,200,-15,10);
    fOutputList->Add(fTPCnsigma_pt_after_tof_its);
    */
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
    
  /*fDCAxy_pt_Dalitz = new TH2F("fDCAxy_pt_Dalitz","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_Dalitz->Sumw2();
  fOutputList->Add(fDCAxy_pt_Dalitz);
  
  fDCAxy_pt_DalitzFromFeeddown = new TH2F("fDCAxy_pt_DalitzFromFeeddown","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_DalitzFromFeeddown->Sumw2();
  fOutputList->Add(fDCAxy_pt_DalitzFromFeeddown);
  
  fDCAxy_pt_Conversions = new TH2F("fDCAxy_pt_Conversions","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_Conversions->Sumw2();
  fOutputList->Add(fDCAxy_pt_Conversions);
  
  fDCAxy_pt_ConversionsFromFeeddown = new TH2F("fDCAxy_pt_ConversionsFromFeeddown","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_ConversionsFromFeeddown->Sumw2();
  fOutputList->Add(fDCAxy_pt_ConversionsFromFeeddown);
  
  fDCAxy_pt_ConversionsFromStrangeFeeddown = new TH2F("fDCAxy_pt_ConversionsFromStrangeFeeddown","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_ConversionsFromStrangeFeeddown->Sumw2();
  fOutputList->Add(fDCAxy_pt_ConversionsFromStrangeFeeddown);
  */
  
  fDCAxy_pt_Dalitz2 = new TH2F("fDCAxy_pt_Dalitz2","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_Dalitz2->Sumw2();
  fOutputList->Add(fDCAxy_pt_Dalitz2);
  
 /* fDCAxy_pt_DalitzFromFeeddown2 = new TH2F("fDCAxy_pt_DalitzFromFeeddown2","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_DalitzFromFeeddown2->Sumw2();
  fOutputList->Add(fDCAxy_pt_DalitzFromFeeddown2);
  */
  fDCAxy_pt_Conversions2 = new TH2F("fDCAxy_pt_Conversions2","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_Conversions2->Sumw2();
  fOutputList->Add(fDCAxy_pt_Conversions2);
  
  fDCAxy_pt_Beauty2 = new TH2F("fDCAxy_pt_Beauty2","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_Beauty2->Sumw2();
  fOutputList->Add(fDCAxy_pt_Beauty2);
  
  fDCAxy_pt_Charm2 = new TH2F("fDCAxy_pt_Charm2","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_Charm2->Sumw2();
  fOutputList->Add(fDCAxy_pt_Charm2);
  
  
  hEleLc_WCorr_WFCorr = new TH2F("hEleLc_WCorr_WFCorr","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
  hEleLc_WCorr_WFCorr->Sumw2();
  fOutputList->Add(hEleLc_WCorr_WFCorr);
  
  hEleDstar_WCorr_WFCorr = new TH2F("hEleDstar_WCorr_WFCorr","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
  hEleDstar_WCorr_WFCorr->Sumw2();
  fOutputList->Add(hEleDstar_WCorr_WFCorr);
  
  hEleDs_WCorr_WFCorr = new TH2F("hEleDs_WCorr_WFCorr","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
  hEleDs_WCorr_WFCorr->Sumw2();
  fOutputList->Add(hEleDs_WCorr_WFCorr);
  
  hEleDp_WCorr_WFCorr = new TH2F("hEleDp_WCorr_WFCorr","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
  hEleDp_WCorr_WFCorr->Sumw2();
  fOutputList->Add(hEleDp_WCorr_WFCorr);
  
  hEleD0_WCorr = new TH2F("hEleD0_WCorr","; p_{T} [GeV/c]; Count",32,ptbinning,2000,-0.5,0.5);
  hEleD0_WCorr->Sumw2();
  fOutputList->Add(hEleD0_WCorr);
      
      
    ///THnSparse to store DCA of different particle species in MC-------------
    Int_t nBinspdg2 = 30;
    Double_t minpdg2 = 0.;
    Double_t maxpdg2 = 30.;
    Double_t binLimpdg2[nBinspdg2+1];
    for(Int_t i=0; i<=nBinspdg2; i++) binLimpdg2[i]=(Double_t)minpdg2 + (maxpdg2-minpdg2)/nBinspdg2*(Double_t)i ;
    
  //  Int_t nBinsdcaxy = 3200; //0.000125 cm
    Int_t nBinsdcaxy = 2000; //0.000125 cm
   // Double_t mindcaxy = -0.2;
   // Double_t maxdcaxy = 0.2;
    Double_t mindcaxy = -0.5;
    Double_t maxdcaxy = 0.5;
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
    
    Int_t nBinstype = 9;
    Double_t mintype = -1.;
    Double_t maxtype = 8.;
    Double_t binLimtype[nBinstype+1];
    for(Int_t i=0; i<=nBinstype; i++) binLimtype[i]=(Double_t)mintype + (maxtype-mintype)/nBinstype*(Double_t)i ;
    
    /*
     Double_t binLimpt2[45]= {0.1,0.112797,0.127231,0.143512,0.161877,0.182592,0.205957,0.232313,0.262041,0.295573,0.333397,
     0.37606,0.424183,0.478465,0.539692,0.608754,0.686654,0.774523,0.873636,0.985432,1.11153,1.25377,1.41421,1.59519,1.79932,2.02957,
     2.28928,2.58223,2.91267,3.2854,3.70582,4.18004,4.71494,5.3183,5.99886,6.76651,7.6324,8.60909,9.71076,10.9534,12.3551,13.9361,15.7195,17.731,20};//bin limits from the measured pi0 spectrum
     */
    if(fIsMC){
    const Int_t nDima2=11;
    Int_t nBina2[nDima2] = {32,nBinspdg2,nBinsdcaxy,nBinsg,nBinsR,nBinsITSchi2,nBinsITSsha,nBinstype,nBinspdg2,32,nBinsdcaxy};
    fD0 = new THnSparseF("fD0","fD0",nDima2,nBina2);
    fD0->SetBinEdges(0,ptbinning); ///pt spectra -> same binning as other histograms
    fD0->SetBinEdges(1,binLimpdg2); /// storing particles (charm and beauty) before correction:
    fD0->SetBinEdges(2,binLimdcaxy); ///dca distribution without Manual Mean and Sigma correction
    fD0->SetBinEdges(3,binLimg);  ///From which generator (Hijing, else, pi0, eta)
    fD0->SetBinEdges(4,binLimR); ///Position where the electron is created
    fD0->SetBinEdges(5,binLimITSchi2); ///ITS chi2 
    fD0->SetBinEdges(6,binLimITSsha); ///fraction ITS shared clusters 
    fD0->SetBinEdges(7,binLimtype); ///pi0 and eta type  ///kNoMother, kNoFeedDown, kNoIsPrimary, kLightMesons, kBeauty, kCharm, kKaonFromHF, kKaonFromNonHF
    fD0->SetBinEdges(8,binLimpdg2);  /// electrons from D,charm baryons, B, beauty baryons, gamma, pi0, eta, Dcorrected, Dcorrected by weight, protons, kaons, D0_corr, D+-_corr,Ds_corr,Lc_corr, D0, D+-,Ds,Lc
    fD0->SetBinEdges(9,ptbinning); ///dca distribution with Manual Mean and Sigma correction
    fD0->SetBinEdges(10,binLimdcaxy);
    fD0->Sumw2();
    fOutputList->Add(fD0);
    }
    ///-----------------------------------------------------------------
    
    if(!fIsMC){
    ///THnSparse to store DCA in Data
    const Int_t nDima3=6;
    Int_t nBina3[nDima3] = {32,nBinsdcaxy,nBinsITSchi2,nBinsITSsha,nBinspdg2,nBinspdg2};
    fD0Data = new THnSparseF("fD0Data","fD0Data",nDima3,nBina3);
    fD0Data->SetBinEdges(0,ptbinning); ///pt spectra -> same binning as other histograms
    fD0Data->SetBinEdges(1,binLimdcaxy); ///dca distribution
    fD0Data->SetBinEdges(2,binLimITSchi2); ///ITS chi2 
    fD0Data->SetBinEdges(3,binLimITSsha); ///fraction ITS shared clusters 
    fD0Data->SetBinEdges(4,binLimpdg2); /// electrons and pions
    fD0Data->SetBinEdges(5,binLimpdg2); /// electrons and pions
    fD0Data->Sumw2();
    fOutputList->Add(fD0Data);
   } 
    /*
    const Int_t nDima4=5;
    Int_t nBina4[nDima4] = {32,nBinsdcaxy,nBinsITSchi2,nBinsITSsha,nBinspdg2};
    fD0HC = new THnSparseF("fD0HC","fD0HC",nDima4,nBina4);
    fD0HC->SetBinEdges(0,ptbinning); ///pt spectra -> same binning as other histograms
    fD0HC->SetBinEdges(1,binLimdcaxy); ///dca distribution
    fD0HC->SetBinEdges(2,binLimITSchi2); ///ITS chi2 
    fD0HC->SetBinEdges(3,binLimITSsha); ///fraction ITS shared clusters 
    fD0HC->SetBinEdges(4,binLimpdg2); /// electrons and pions
    fD0HC->Sumw2();
    fOutputList->Add(fD0HC);
    */
    
    PostData(1, fOutputList);
    
    
}

//////////////
//Main loop///
//////////////

///Called for each event
void AliAnalysisHFEppTPCTOFBeauty::UserExec(Option_t *)
{

//cout<<"=======================================INSIDE THE UserExec======================================="<<endl;
	gRandom->SetSeed(0);
	//gRandom = new TRandom3(0);
    Int_t pdg = -99999;
    Int_t pdg_mother = -99999;
    Double_t fTPCnSigma = -999;
    Double_t fTOFnSigma = -999;
   // Double_t fITSnSigma = -999;
    Double_t fTPCnSigma_pion = -999;
    Double_t fTPCnSigma_proton = -999;
    Double_t fTOFnSigma_proton = -999;
    Double_t fTPCnSigma_kaon = -999;
    Double_t fTPCsignal = -999;
    Double_t fPt = -999;
    Double_t fEta = -999;
    Double_t fPhi = -999;
    Double_t fP = -999;
   // Int_t fITSnClus = 99999;
    Int_t fTPCnClus = 99999;
    Double_t qadca[10];
	Double_t qadcaData[10];
	Double_t qadcaHC[10];
    for(int j=0;j<10;j++){
		qadca[j] = 0;
		qadcaData[j] = 0;
		qadcaHC[j] = 0;
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
           
        for(Int_t iMC = 0; iMC < fMCarray->GetEntries(); iMC++)
        {
			fMCparticle = (AliAODMCParticle*) fMCarray->At(iMC);
			///Pseudo-rapidity cut
			if((fMCparticle->Eta() < fEtaMin) || (fMCparticle->Eta() > fEtaMax)) continue;     				
			
            ///Get proportion between generated D0 and LambdaC:
            Int_t TrackPDG = TMath::Abs(fMCparticle->GetPdgCode());
            if(TrackPDG == 421) hPtD0->Fill(fMCparticle->Pt());///D0
            if(TrackPDG == 4122) hPtLambdaC->Fill(fMCparticle->Pt()); ///LambdaC
            if(TrackPDG == 411) hPtDp->Fill(fMCparticle->Pt()); ///D+
            if(TrackPDG == 413) hPtDstar->Fill(fMCparticle->Pt()); ///D*
            if(TrackPDG == 431) hPtDs->Fill(fMCparticle->Pt()); ///Ds
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

        ///////////////////////////
		//BEFORE TRACK SELECTION///
		///////////////////////////
		fPt = track->Pt();
        fEta = track->Eta();
        fPhi = track->Phi();
        fP = TMath::Sqrt((track->Pt())*(track->Pt()) + (track->Pz())*(track->Pz()));
      /*  
        if(fP <= 2.0 ){
        ftofPIDmincut = -2;
        ftofPIDmaxcut = 2;
        }
        else{
        ftofPIDmincut = -3;
        ftofPIDmaxcut = 3;
        }
        */
        //fITSnClus =  track->GetNumberOfITSClusters();
        fTPCnClus =  track->GetNumberOfTPCClusters();
        
        fTPCsignal = track->GetTPCsignal();
        fTPCnSigma = fPidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        fTOFnSigma = fPidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
       // fITSnSigma = fPidResponse->NumberOfSigmasITS(track, AliPID::kElectron);
        fTPCnSigma_pion = fPidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
        fTPCnSigma_proton = fPidResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        fTOFnSigma_proton = fPidResponse->NumberOfSigmasTOF(track, AliPID::kProton);
        fTPCnSigma_kaon = fPidResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        
        fTPC_p1->Fill(fP,fTPCsignal);
        fTPCnsigma_p1->Fill(fP,fTPCnSigma);
        fTPCnsigma_pt1->Fill(fPt,fTPCnSigma);
        
        fTPCnsigma_TOFnsigma1->Fill(fTOFnSigma,fTPCnSigma);
        fTOFnsigma_p1->Fill(fP,fTOFnSigma);
        fTOFnsigma_pt1->Fill(fPt,fTOFnSigma);
        
        //fITSnsigma_p1->Fill(fP,fITSnSigma);
        //fITSnsigma_pt1->Fill(fPt,fITSnSigma);
        
        fPt_1->Fill(fPt);
        
       // fITSnClus_1->Fill(fITSnClus);
        fTPCnClus_1->Fill(fTPCnClus);

        ///Pseudo-rapidity cut
        if((track->Eta() < fEtaMin) || (track->Eta() > fEtaMax)) continue;
        
        ///Beauty reconstruction efficiency block-----------
        if(fIsMC && fIsAOD){
			Bool_t IsHFEMC = IsHFelectronsMC(track);
			if(IsHFEMC){
				if(fIsFromMesonB || fIsFromBarionB || fIsFromBarionBD || fIsFromMesonBD){
					fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
					fPtBeautyReconstructedAll->Fill(fMCparticle->Pt());
					//cout<<"reconstructed by track cut"<<endl;
				}
			}
		}
		///----------------------------------------------------
        
        
        //=======================================================================
        // Track Selection Cuts are applied here
        //=======================================================================        
        ///Test Filter Bit
         if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;
        
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
        

         ////////////////////
         //Calculating DCA///
         ////////////////////
         //
         if(!fExtraCuts){
                fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
          }
         fExtraCuts->SetRecEventInfo(fAOD);
         Double_t d0z0[2]={-999,-999}, cov[3]={999,999,999};
         AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
          if(!(track->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 3., d0z0, cov))) continue; 
         //cout<<d0z0[0]<<"    "<<d0z0[1]<<endl;
        // fExtraCuts->GetHFEImpactParameters(track, d0z0, cov); // recalculation of vertex is done here, earlier was not done in the task and was the reason for "shoulder" shape in the DCA templates. Also, this is not giving effect for PbPb but only pp. ====> Sudhir 19 January, 2019 ///Solved
                                                                                                                                                         		Double_t DCAxy = d0z0[0];
                                                                                                                                                         		Double_t DCAz = d0z0[1];
	
	if(TMath::Abs(DCAxy) > 1.0 || TMath::Abs(DCAz) > 2.0) continue;
         

        
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
   	
   	/*Int_t trkLabel = TMath::Abs(track->GetLabel());
		Int_t labelm = GetPrimary(trkLabel,fMCarray);///gives the label of first mother
		AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCarray->At(labelm);
		Int_t trkIndexPrimHFE= AODMCtrack->GetLabel();///gives index of the particle in original MCparticle array (labelm and trkIndexPrimHFE are the same, so I don't understand why this is done)
   
		if(trkIndexPrimHFE < fNpureMC){
			NAnalizedTracksHijing = NAnalizedTracksHijing+1;     // commented by Sudhir 14 june 2018 since code was not working with trigger kINT7 so it was commented out and this piece was giving an error === only for the check....
		}*/

        fTPC_p2->Fill(fP,fTPCsignal);
        fTPCnsigma_TOFnsigma2->Fill(fTOFnSigma,fTPCnSigma);
        fTOFnsigma_p2->Fill(fP,fTOFnSigma);
        fTOFnsigma_pt2->Fill(fPt,fTOFnSigma);
        fTPCnsigma_p2->Fill(fP,fTPCnSigma);
        fTPCnsigma_pt2->Fill(fPt,fTPCnSigma);
       // fITSnsigma_p2->Fill(fP,fITSnSigma);
       // fITSnsigma_pt2->Fill(fPt,fITSnSigma);
        fPt_2->Fill(fPt);
       // fITSnClus_2->Fill(fITSnClus);
        fTPCnClus_2->Fill(fTPCnClus);
        
         
  	if(fTOFnSigma_proton >= -3 && fTOFnSigma_proton <= 3){
  	fTPCnsigma_proton_p_after_tof->Fill(fP,fTPCnSigma_proton);
  	}
        ///Checking nsigmaTPC after PID cuts in tof and its
        if(fTOFnSigma >= ftofPIDmincut && fTOFnSigma <= ftofPIDmaxcut){
            fTPCnsigma_p_after_tof->Fill(fP,fTPCnSigma);
            
            fTPCnsigma_pt_after_tof->Fill(fPt,fTPCnSigma);
            
			///Cheking the hadron nsigmaTPC after TOF cut
			if(fIsMC){
				fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
				if(TMath::Abs(fMCparticle->GetPdgCode()) == 2212)  fTPCnsigma_p_after_tof_p->Fill(fP,fTPCnSigma);
				if(TMath::Abs(fMCparticle->GetPdgCode()) == 130 || TMath::Abs(fMCparticle->GetPdgCode()) == 310 || TMath::Abs(fMCparticle->GetPdgCode()) == 311 || TMath::Abs(fMCparticle->GetPdgCode()) == 321)  fTPCnsigma_p_after_tof_k->Fill(fP,fTPCnSigma);
				if(TMath::Abs(fMCparticle->GetPdgCode()) == 211)  fTPCnsigma_p_after_tof_pion->Fill(fP,fTPCnSigma);
				if(TMath::Abs(fMCparticle->GetPdgCode()) == 211) fPionsPt->Fill(fPt);
				if(TMath::Abs(fMCparticle->GetPdgCode()) == 321) fKaonsPt->Fill(fPt);
				}
			
            /*if(fITSnSigma >= -2 && fITSnSigma <= 2){
                fTPCnsigma_p_after_tof_its->Fill(fP,fTPCnSigma);
                fTPCnsigma_pt_after_tof_its->Fill(fPt,fTPCnSigma);
            }*/
        }
        
        
        //////////////////////////////////////////////
		//Hadron DCA  --- for DCA resolution studies//
		//////////////////////////////////////////////
		
		
		 ///////////////////
		// Without PID cuts //     // By Sudhir on Monday Oct 8, 2018 For Improver check
	       ///////////////////
		fDCAxy_pt_had_onlyDCA_WoPID->Fill(fPt,DCAxy);
			fDCAxy_pt_had_WoPID->Fill(fPt,DCAxy*track->Charge()*signB);
			fDCAz_pt_had_WoPID->Fill(fPt,TMath::Sqrt(cov[0]));	

		 ///////////////////
		// With PID cuts //
	       ///////////////////
	       
	       Double_t phi_d0 = (track->Phi()*180.0/TMath::Pi());
        if(fTPCnSigma >= -5 && fTPCnSigma <= -3){
			fDCAxy_pt_had_onlyDCA->Fill(fPt,DCAxy);
			fDCAxy_pt_had->Fill(fPt,DCAxy*track->Charge()*signB);
			fDCAz_pt_had->Fill(fPt,TMath::Sqrt(cov[0]));	
			
			if(phi_d0 > 315.0 || phi_d0 < 45.0){
			fDCAxy_pt_had_onlyDCA_phi1->Fill(fPt,DCAxy);
			fDCAxy_pt_had_phi1_ChB->Fill(fPt,DCAxy*track->Charge()*signB);
			fDCAxy_pt_had_phi1_B->Fill(fPt,TMath::Sqrt(cov[0]));
			
			float DCAMCRes_phi1 = GetDCAResolMC_phi1(fPt); ///resolution of the MC
			float DCAMCMean_phi1 = GetDCAMeanMC_phi1(fPt); ///mean of the MC
			float correction_phi1 = gRandom->Gaus(DCAMCMean_phi1,DCAMCRes_phi1); 
			float DCAResCorr_phi1 =  DCAxy + correction_phi1;
			
			//fResGausCorr_phi1->Fill(correction_phi1);
			fDCAxy_pt_had_ResCorr_phi1->Fill(fPt,DCAResCorr_phi1);
			//cout<<"Phi1 value:===   "<<phi_d0<<endl;
			}
			
			if(phi_d0 > 45.0 && phi_d0 < 135.0){
			fDCAxy_pt_had_onlyDCA_phi2->Fill(fPt,DCAxy);
			fDCAxy_pt_had_phi2_ChB->Fill(fPt,DCAxy*track->Charge()*signB);
			fDCAxy_pt_had_phi2_B->Fill(fPt,TMath::Sqrt(cov[0]));
			float DCAMCRes_phi2 = GetDCAResolMC_phi2(fPt); ///resolution of the MC
			float DCAMCMean_phi2 = GetDCAMeanMC_phi2(fPt); ///mean of the MC
			float correction_phi2 = gRandom->Gaus(DCAMCMean_phi2,DCAMCRes_phi2); 
			float DCAResCorr_phi2 =  DCAxy + correction_phi2;
			
			//fResGausCorr_phi2->Fill(correction_phi2);
			fDCAxy_pt_had_ResCorr_phi2->Fill(fPt,DCAResCorr_phi2);
			//cout<<"Phi2 value:===   "<<phi_d0<<endl;
			}
			
			if(phi_d0 > 135.0 && phi_d0 < 225.0){
			fDCAxy_pt_had_onlyDCA_phi3->Fill(fPt,DCAxy);
			fDCAxy_pt_had_phi3_ChB->Fill(fPt,DCAxy*track->Charge()*signB);
			fDCAxy_pt_had_phi3_B->Fill(fPt,TMath::Sqrt(cov[0]));
			float DCAMCRes_phi3 = GetDCAResolMC_phi3(fPt); ///resolution of the MC
			float DCAMCMean_phi3 = GetDCAMeanMC_phi3(fPt); ///mean of the MC
			float correction_phi3 = gRandom->Gaus(DCAMCMean_phi3,DCAMCRes_phi3); 
			float DCAResCorr_phi3 =  DCAxy + correction_phi3;
			
			//fResGausCorr_phi3->Fill(correction_phi3);
			fDCAxy_pt_had_ResCorr_phi3->Fill(fPt,DCAResCorr_phi3);
			//cout<<"Phi3 value:===   "<<phi_d0<<endl;
			}
			
			if(phi_d0 > 225.0 && phi_d0 < 315.0){
			fDCAxy_pt_had_onlyDCA_phi4->Fill(fPt,DCAxy);
			fDCAxy_pt_had_phi4_ChB->Fill(fPt,DCAxy*track->Charge()*signB);
			fDCAxy_pt_had_phi4_B->Fill(fPt,TMath::Sqrt(cov[0]));
			float DCAMCRes_phi4 = GetDCAResolMC_phi4(fPt); ///resolution of the MC
			float DCAMCMean_phi4 = GetDCAMeanMC_phi4(fPt); ///mean of the MC
			float correction_phi4 = gRandom->Gaus(DCAMCMean_phi4,DCAMCRes_phi4); 
			float DCAResCorr_phi4 =  DCAxy + correction_phi4;
			
			//fResGausCorr_phi4->Fill(correction_phi4);
			fDCAxy_pt_had_ResCorr_phi4->Fill(fPt,DCAResCorr_phi4);
			//cout<<"Phi4 value:===   "<<phi_d0<<endl;
			}
			
			
			///Checking the effect of the improver in the resolution for hijing events separetely
			if(fIsMC){
				Int_t trkLabel = TMath::Abs(track->GetLabel());
				Int_t labelm = GetPrimary(trkLabel,fMCarray);///gives the label of first mother
				AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCarray->At(labelm);
				Int_t trkIndexPrimHFE= AODMCtrack->GetLabel();///gives index of the particle in original MCparticle array (labelm and trkIndexPrimHFE are the same, so I don't understand why this is done)
            
				///Particle generator:
				if(trkIndexPrimHFE < fNpureMC){
					//fDCAxy_pt_had_onlyDCA_Hijing->Fill(fPt,DCAxy);
				}
				if(trkIndexPrimHFE >= fNpureMC){
					fDCAxy_pt_had_onlyDCA_Phytia->Fill(fPt,DCAxy);
				}			
			}
        }
        
        ///////////////////////////////////////////////////
		//THnSparse to store the DCA information of Data///
		///////////////////////////////////////////////////
         if(!fIsMC){
			 qadcaData[0] = fPt;
         
			 qadcaData[1] = DCAxy*track->Charge()*signB;
			
			 qadcaData[4] = -1.;
			
			 qadcaData[5] = 10.;
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
			 
			  ///Proton candidates
			 if(fTPCnSigma_proton >= -3 && fTPCnSigma_proton <= 8){
				if(fTOFnSigma_proton >= -3 && fTOFnSigma_proton <= 3){
				//hDCAPtProtons->Fill(fPt,DCAxy*track->Charge()*signB);
					//qadcaData[4] = 3.5;
										
				}
			 }
			 
			   ///Proton candidates
			 if(fTPCnSigma_proton >= -3 && fTPCnSigma_proton <= 3){
				if(fTOFnSigma_proton >= -3 && fTOFnSigma_proton <= 3){
				//hDCAPtProtons3->Fill(fPt,DCAxy*track->Charge()*signB);
					qadcaData[5] = 3.5;
										
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
         //cout<<fPt<<endl;
			if(qadcaData[4]>0.) fD0Data->Fill(qadcaData);
        }
        
       
       
       
        ///////////////////////////////////////////////////
		//THnSparse to store the DCA information of Hadron Contamination///
		///////////////////////////////////////////////////
         //if(!fIsMC){
		/*	 qadcaHC[0] = fPt;
         
			 qadcaHC[1] = DCAxy*track->Charge()*signB;
			
			 qadcaHC[4] = -1.;
			
			 ///Charged pions
			 if(fTPCnSigma >= -5 && fTPCnSigma <= -3){
				  qadcaHC[4] = 0.5;
			 }
			 ///Electron candidates
			 if(fTPCnSigma >= ftpcPIDmincut && fTPCnSigma <= ftpcPIDmaxcut){
				if(fTOFnSigma >= ftofPIDmincut && fTOFnSigma <= ftofPIDmaxcut){
					qadcaHC[4] = 1.5;					
				}
			 }
        
			Double_t ITSNcls2 = atrack->GetITSNcls();
			//cout<<"atrack->GetITSNcls() = "<<atrack->GetITSNcls()<<endl;
            
			///ITS Chi2            
			qadcaHC[2] = atrack->GetITSchi2()/ITSNcls2; 
			//cout<<"track->GetITSchi2() = "<<track->GetITSchi2()<<endl;
            
			///Fraction of shared clusters in the ITS
			Bool_t HasSharedCls2 = kFALSE;
			Double_t ITSNSharedcls2 = 0;
			for(int itsL = 0; itsL < 6; itsL++){
				HasSharedCls2 = atrack->HasSharedPointOnITSLayer(itsL);
				if(HasSharedCls2) ITSNSharedcls2++;
			}
			//cout<<"ITSNSharedcls = "<<ITSNSharedcls<<endl;
         
			Double_t fsharedclsITS2 = ITSNSharedcls2/ITSNcls2;
			//cout<<"fsharedclsITS = "<<fsharedclsITS<<endl;
            
			qadcaHC[3] = fsharedclsITS2; 
         		Double_t WeightHC = fHC->Eval(fP);
         		//cout<<WeightHC<<endl;
			if(qadcaHC[4]>0.) fD0HC->Fill(qadcaHC, WeightHC);
        //}
        */
       
        
               
		///Checking PID cuts separately (for the efficiency)	
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
   //     if(fTOFnSigma >= ftofPIDmincut && fTOFnSigma <= ftofPIDmaxcut){
  /* if(fIsMC && fIsAOD){
   
   	 if(TMath::Abs(fMCparticle->GetPdgCode()) == 2212){ 
            hDCAPtProtons2->Fill(fPt,DCAxy*track->Charge()*signB);
           // cout<<"Hello:  "<<TMath::Abs(fMCparticle->GetPdgCode())<<endl;
            }
   
         if(fTPCnSigma >= ftpcPIDmincut && fTPCnSigma <= ftpcPIDmaxcut){
        if(TMath::Abs(fMCparticle->GetPdgCode()) == 2212){ 
         //   hDCAPtProtons->Fill(fPt,DCAxy*track->Charge()*signB);
         //   cout<<"Hello:  "<<TMath::Abs(fMCparticle->GetPdgCode())<<endl;
            }

        }
        
        if(fTPCnSigma >= -5 && fTPCnSigma <= 5){
        if(TMath::Abs(fMCparticle->GetPdgCode()) == 2212){ 
         //   hDCAPtProtons3->Fill(fPt,DCAxy*track->Charge()*signB);
         //   cout<<"Hello:  "<<TMath::Abs(fMCparticle->GetPdgCode())<<endl;
            }

        }
      }*/
        //=======================================================================
        // Here the PID cuts defined in the file "Config.C" is applied
        //=======================================================================
        Int_t pidpassed = 1;
        AliHFEpidObject hfetrack;
        hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
        hfetrack.SetRecTrack(track);
        hfetrack.SetPP();	//proton-proton analysis
        if(!fPID->IsSelected(&hfetrack, NULL, "", fPIDqa)) pidpassed = 0;
        //cout<<"Before the pidpassed = = = = = = = = = = = ======================================"<<endl;
        if(pidpassed==0) continue;
        //cout<<"After  the pidpassed = = = = = = = = = = = ======================================"<<endl;
        
        //if(fTPCnSigma < ftpcPIDmincut || fTPCnSigma > ftpcPIDmaxcut && fTOFnSigma < ftofPIDmincut || fTOFnSigma > ftofPIDmaxcut) continue;   // Applying simultaneous pid cuts manually by sudhir, since above cut function does not work in grid...
        //if(fTOFnSigma == -999 || fTPCnSigma == -999)continue;
        //cout<<ftpcPIDmincut<<"    "<<fTPCnSigma<<"     "<<ftpcPIDmaxcut<<"       "<<ftofPIDmincut<<"     "<<fTOFnSigma<<"     "<<ftofPIDmaxcut<<endl;
        //if(fTPCnSigma >= ftpcPIDmincut && fTPCnSigma <= ftpcPIDmaxcut && fTOFnSigma >= ftofPIDmincut && fTOFnSigma <= ftofPIDmaxcut){
        // cout<<fTPCnSigma<<"       "<<fTOFnSigma<<endl;
        
       // cout<<fP<<"      "<<ftofPIDmincut<<"      "<<ftofPIDmaxcut<<endl;
                /////////////////////////
		//AFTER PID SELECTION////
		/////////////////////////
        ///Beauty reconstruction efficiency block-----------
         if(fIsMC && fIsAOD){
			Bool_t IsHFEMC = IsHFelectronsMC(track);
			if(IsHFEMC){
				if(fIsFromMesonB || fIsFromBarionB || fIsFromBarionBD || fIsFromMesonBD){
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
        
        //fITSnsigma_p3->Fill(fP,fITSnSigma);
        //fITSnsigma_pt3->Fill(fPt,fITSnSigma);
        
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
            qadca[9] = fMCparticle->Pt();
            ///Selecting particle
            qadca[1]=-1.; 
             //cout<<"Hello:"<<endl;           
            qadca[8]=29.5; //if noone passes the correction HFE then 29.5 is filled
            ///------------
            if(TMath::Abs(fMCparticle->GetPdgCode()) == 2212){ 
            qadca[1]=10.5; ///to check DCA of protons
            //hDCAPtProtons->Fill(fPt,DCAxy*track->Charge()*signB);
            //cout<<"Hello:"<<endl;
            }
            if(TMath::Abs(fMCparticle->GetPdgCode()) == 321) qadca[1]=11.5; ///to check DCA of kaons
            ///------------

            /////////////////////////
			//Electrons from charm///
			/////////////////////////
            Int_t PhotonicType2 = 999;
            /*
           // if(!fIsMC){
             AliHFEsignalCuts *fSignalCuts = new AliHFEsignalCuts("HFEsignalCuts", "HFE MC Signal definition");
      	     if(fMCarray){
      	     fSignalCuts->SetMCAODInfo(fMCarray);
             }
      	    PhotonicType2 = fSignalCuts->GetSignalSource(track);*/
         //   }	   
	/*
            fSignalCuts = new AliHFEsignalCuts("HFEsignalCuts", "HFE MC Signal definition");
             if(fMCarray){
             fSignalCuts->SetMCAODInfo(fMCarray);
             }
            PhotonicType2 = fSignalCuts->GetSignalSource(track);
 */
 
 	    ///DCAxy
            	float DCAResCorr = 999;
            	
           	if(phi_d0 > 315.0 || phi_d0 < 45.0){
		float DCAMCRes_phi1 = GetDCAResolMC_phi1(fPt); ///resolution of the MC
		float DCAMCMean_phi1 = GetDCAMeanMC_phi1(fPt); ///mean of the MC
		float correction_phi1 = gRandom->Gaus(DCAMCMean_phi1,DCAMCRes_phi1); 
		DCAResCorr =  (DCAxy + correction_phi1)*track->Charge()*signB;
		}
			
		if(phi_d0 > 45.0 && phi_d0 < 135.0){
		float DCAMCRes_phi2 = GetDCAResolMC_phi2(fPt); ///resolution of the MC
		float DCAMCMean_phi2 = GetDCAMeanMC_phi2(fPt); ///mean of the MC
		float correction_phi2 = gRandom->Gaus(DCAMCMean_phi2,DCAMCRes_phi2); 
		DCAResCorr =  (DCAxy + correction_phi2)*track->Charge()*signB;
		}
			
		if(phi_d0 > 135.0 && phi_d0 < 225.0){
		float DCAMCRes_phi3 = GetDCAResolMC_phi3(fPt); ///resolution of the MC
		float DCAMCMean_phi3 = GetDCAMeanMC_phi3(fPt); ///mean of the MC
		float correction_phi3 = gRandom->Gaus(DCAMCMean_phi3,DCAMCRes_phi3); 
		DCAResCorr =  (DCAxy + correction_phi3)*track->Charge()*signB;
		}
			
		if(phi_d0 > 225.0 && phi_d0 < 315.0){
		float DCAMCRes_phi4 = GetDCAResolMC_phi4(fPt); ///resolution of the MC
		float DCAMCMean_phi4 = GetDCAMeanMC_phi4(fPt); ///mean of the MC
		float correction_phi4 = gRandom->Gaus(DCAMCMean_phi4,DCAMCRes_phi4);  
		DCAResCorr =  (DCAxy + correction_phi4)*track->Charge()*signB;
		}
 
 
	    Bool_t IsHFEMC = IsHFelectronsMC(track);
            if(IsHFEMC){
            
               Double_t fWeight_Dp = 999;
               Double_t fWeight_Ds = 999;
               Double_t fWeight_Dstar = 999;
               Double_t fWeight_Lc = 999;
     // AliHFEmcQA *hfemcqa = new AliHFEmcQA();
     // if(hfemcqa) hfemcqa->SetMCArray(fMCarray);
            
                if(fIsFromD){
                    fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
                    fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                    fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
                    fMCparticleGGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleGMother->GetMother());
                    pdg_mother = fMCparticleMother->GetPdgCode();
                    float pdg_gmother = fMCparticleGMother->GetPdgCode();
                    float pdg_ggmother = fMCparticleGGMother->GetPdgCode();
                    
                    if(TMath::Abs(pdg_mother)>400 && TMath::Abs(pdg_mother)<500){///charmed meson 
						qadca[1]=7.5;
						if(TMath::Abs(pdg_mother) == 421) qadca[1]=20.5; ///to check DCA  prompt D0 bef corr
						if(TMath::Abs(pdg_mother) == 411) qadca[1]=23.5; ///to check DCA D+ bef corr
						
						if(TMath::Abs(pdg_mother) == 431) qadca[1]=25.5; ///to check DCA Ds+ bef corr
						
						if(TMath::Abs(pdg_gmother) == 413) qadca[1]=26.5; ///to check DCA Dstar+ bef corr
						if(TMath::Abs(pdg_gmother) == 413){ 
						if(TMath::Abs(pdg_mother) == 421){ 
						qadca[1]=21.5; ///to check DCA Dstar+ bef corr
						}
						}
						 
						if(TMath::Abs(pdg_gmother) == 413){ 
						if(TMath::Abs(pdg_mother) == 411){ 
						qadca[1]=22.5; ///to check DCA Dstar+ bef corr
						}
						}
						
						if(TMath::Abs(pdg_gmother) == 413){ 
						if(TMath::Abs(pdg_mother) == 431){ 
						qadca[1]=24.5; ///to check DCA Dstar+ bef corr
						}
						}
						 hCharmMotherPt->Fill(fMCparticleMother->Pt());
						 hCharmMotherPt_vsElecPt->Fill(fPt,fMCparticleMother->Pt());
						 hElecPt_vsCharmMotherPt->Fill(fMCparticleMother->Pt(),fPt);
						 fDCAxy_pt_charmbef->Fill(fPt,DCAxy*track->Charge()*signB);
						 
						 ///Correcting pT spectrum
						 ///(Weight for each electron pT bin starting at one - charm correction)
						 probAcceptD = -999;
						 if(fMCparticleMother->Pt() > 36) continue;
						 
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
						qadca[8]=7.5;
						if(TMath::Abs(pdg_mother) == 421){ 
						qadca[8]=20.5; ///to check DCA  prompt D0 aft corr
						hEleD0_WCorr->Fill(fPt,DCAResCorr);
						}
						if(TMath::Abs(pdg_mother) == 411){ 
						qadca[8]=23.5; ///to check DCA D+ aft corr
						if(fMCparticleMother->Pt() > 30) continue;
						fWeight_Dp = HistWeightDpD0->Eval(fMCparticleMother->Pt());
						hEleDp_WCorr_WFCorr->Fill(fPt, DCAResCorr, fWeight_Dp);
						}
						
						if(TMath::Abs(pdg_mother) == 431){ 
						qadca[8]=25.5; ///to check DCA Ds+ aft corr
						if(fMCparticleMother->Pt() > 30) continue;
						fWeight_Ds = HistWeightDsD0->Eval(fMCparticleMother->Pt());
						hEleDs_WCorr_WFCorr->Fill(fPt, DCAResCorr, fWeight_Ds);
						}
						
						if(TMath::Abs(pdg_gmother) == 413){ 
						qadca[8]=26.5; ///to check DCA Dstar+ aft corr
						if(fMCparticleGMother->Pt() > 30) continue;
						fWeight_Dstar = HistWeightDstarD0->Eval(fMCparticleGMother->Pt());
						hEleDstar_WCorr_WFCorr->Fill(fPt, DCAResCorr, fWeight_Dstar);
						}
						if(TMath::Abs(pdg_gmother) == 413){ 
						if(TMath::Abs(pdg_mother) == 421){ 
						qadca[8]=21.5; ///to check DCA DstarD0 aft corr
						if(fMCparticleGMother->Pt() > 30) continue;
						fWeight_Dstar = HistWeightDstarD0->Eval(fMCparticleGMother->Pt());
						hEleDstar_WCorr_WFCorr->Fill(fPt, DCAResCorr, fWeight_Dstar);
						}
						}
						 
						if(TMath::Abs(pdg_gmother) == 413){ 
						if(TMath::Abs(pdg_mother) == 411){ 
						qadca[8]=22.5; ///to check DCA DstarD+ aft corr
						if(fMCparticleMother->Pt() > 30) continue;
						fWeight_Dp = HistWeightDpD0->Eval(fMCparticleMother->Pt());
						hEleDp_WCorr_WFCorr->Fill(fPt, DCAResCorr, fWeight_Dp);
						}
						}
						
						if(TMath::Abs(pdg_gmother) == 413){ 
						if(TMath::Abs(pdg_mother) == 431){ 
						qadca[8]=24.5; ///to check DCA Dstar+ aft corr
						if(fMCparticleMother->Pt() > 30) continue;
						fWeight_Ds = HistWeightDsD0->Eval(fMCparticleMother->Pt());
						hEleDs_WCorr_WFCorr->Fill(fPt, DCAResCorr, fWeight_Ds);
						}
						}
							hCharmMotherPt_corr->Fill(fMCparticleMother->Pt());
							hCharmMotherPt_corr2->Fill(fMCparticleMother->Pt());
							hCharmMotherPt_vsElecPt_corr->Fill(fPt,fMCparticleMother->Pt());
							hElecPt_vsCharmMotherPt_corr->Fill(fMCparticleMother->Pt(),fPt);
							fDCAxy_pt_charmaft->Fill(fPt,DCAxy*track->Charge()*signB);						 
						 }												 
					}
					
							
					
					if(TMath::Abs(pdg_mother)>4000 && TMath::Abs(pdg_mother)<5000){///charmed baryon
						 qadca[1]=1.5;
						 if(TMath::Abs(pdg_mother) == 4122){ 
						 qadca[1]=27.5; ///to check DCA prompt Lc+
						 if(fMCparticleMother->Pt() > 30) continue;
						 fWeight_Lc = HistWeightLcD0->Eval(fMCparticleMother->Pt());
						 hEleLc_WCorr_WFCorr->Fill(fPt, DCAResCorr, fWeight_Lc);
						 }

						
					}
                }//end of electrons from charm
                
                //////////////////////////
				//Electrons from beauty///
				//////////////////////////
                if(fIsFromMesonB || fIsFromBarionB || fIsFromBarionBD || fIsFromMesonBD){
                    fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
                    fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                    pdg_mother = fMCparticleMother->GetPdgCode();
                    
                    if(fIsFromMesonB){
                    fDCAxy_pt_MesonB_beautybef->Fill(fPt,DCAxy*track->Charge()*signB);
                    }
                    if(fIsFromMesonBD){
                    
                    fDCAxy_pt_MesonBD_beautybef->Fill(fPt,DCAxy*track->Charge()*signB);
                    }
                    // check by Sudhir on 27 Oct 2018
                    if(fIsFromBarionB){
                    fDCAxy_pt_BaryonB_beautybef->Fill(fPt,DCAxy*track->Charge()*signB);
                    }
                    if(fIsFromBarionBD){
                    
                    fDCAxy_pt_BaryonBD_beautybef->Fill(fPt,DCAxy*track->Charge()*signB);
                    }
                    // check by Sudhir on 27 Oct 2018
                    if(fIsFromMesonB || fIsFromMesonBD){///beauty meson 
						qadca[1]=2.5;
						hBeautyMotherPt->Fill(fMCparticleMother->Pt(),fPt);
						hBeautyMotherPtbef->Fill(fMCparticleMother->Pt());
						fDCAxy_pt_beautybef->Fill(fPt,DCAxy*track->Charge()*signB);
						
						///Correcting the pT spectrum
						probAcceptB = -999;
						if(fIsFromMesonB && (fMCparticleMother->Pt() < 50)){
							 probAcceptB = fBcorr->Eval(fMCparticleMother->Pt());///always evaluating the function in the pT of the B meson
							 //cout<<"pdg_mother = "<<pdg_mother<<endl;	
						}
						if(fIsFromMesonBD){
							if(fMCparticleMother->GetMother() > 0){
								fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
								float pdg_gmother = fMCparticleGMother->GetPdgCode();
								
								if(TMath::Abs(pdg_gmother)>500 && TMath::Abs(pdg_gmother)<600 && (fMCparticleGMother->Pt()<50)){
									 probAcceptB = fBcorr->Eval(fMCparticleGMother->Pt());
									
								}
								else if(fMCparticleGMother->GetMother() > 0){
									fMCparticleGGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleGMother->GetMother());
									float pdg_ggmother = fMCparticleGGMother->GetPdgCode();
									if(TMath::Abs(pdg_ggmother)>500 && TMath::Abs(pdg_ggmother)<600 && (fMCparticleGGMother->Pt()<50)) probAcceptB = fBcorr->Eval(fMCparticleGGMother->Pt());
									else if(fMCparticleGGMother->GetMother() > 0){
										fMCparticleGGGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleGGMother->GetMother());
										float pdg_gggmother = fMCparticleGGGMother->GetPdgCode();
										if (TMath::Abs(pdg_gggmother)>500 && TMath::Abs(pdg_gggmother)<600 && (fMCparticleGGGMother->Pt()<50)) probAcceptB = fBcorr->Eval(fMCparticleGGGMother->Pt());
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
							hBeautyMotherPtaft->Fill(fMCparticleMother->Pt());
							hBeautyMotherPt2Daft->Fill(fMCparticleMother->Pt(),fPt);
							qadca[8]=2.5;
							
						if(fIsFromMesonB){
                    					fDCAxy_pt_MesonB_beautyaft->Fill(fPt,DCAxy*track->Charge()*signB);
                    				}
                    				
                    				if(fIsFromMesonBD){
                        	                	fDCAxy_pt_MesonBD_beautyaft->Fill(fPt,DCAxy*track->Charge()*signB);
                    				}
						// check by Sudhir on 27 Oct 2018	
						}
					}
					
					if(fIsFromBarionB || fIsFromBarionBD){///beauty baryon
						qadca[1]=3.5;
					}
                }//end of electrons from beauty
                 Double_t DCANew = DCAxy*track->Charge()*signB;
        	//if(fIsMC) GetHFElectronTemplates2(fMCarray, track, fPt, DCANew, PhotonicType2);
        	if(PhotonicType2 == 0){ // electrons from photon conversions from light mesons
      		fDCAxy_pt_Charm2->Fill(fMCparticle->Pt(), DCANew);
      		}
      
      		if(PhotonicType2 == 1){ // electrons from photon conversions from strange
      		fDCAxy_pt_Beauty2->Fill(fMCparticle->Pt(), DCANew);
      		}
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
                
                Double_t DCANew = DCAxy*track->Charge()*signB;
                
        	//if(fIsMC) GetGammaAndDalitzElectronTemplates(fMCarray, track, fPt, DCANew);
        	//if(fIsMC) GetGammaAndDalitzElectronTemplates2(fMCarray, track, fPt, DCANew, PhotonicType2);
                
                if(PhotonicType2 == 3){ // electrons from Dalitz 
      		fDCAxy_pt_Dalitz2->Fill(fPt, DCAResCorr);
		//cout<<"Pdg pTy 3:   "<<pdg_mother<<endl;
      		}
          
      		if(PhotonicType2 == 2){ // electrons from photon conversions
      		fDCAxy_pt_Conversions2->Fill(fPt, DCAResCorr);
		//cout<<"Pdg pTy 2:   "<<pdg_mother<<endl;
      		}
      
                
                ///Photonic Electrons:
                if(TMath::Abs(pdg_mother) == 111 || TMath::Abs(pdg_mother) == 221 || TMath::Abs(pdg_mother) == 22){
					if(TMath::Abs(pdg_mother) == 111) qadca[1]=28.5;
						//cout<<"Pdg pTy old dal1:   "<<pdg_mother<<endl;} 
					if(TMath::Abs(pdg_mother) == 221) qadca[1]=29.5;
						//cout<<"Pdg pTy old dal2:   "<<pdg_mother<<endl;}
					if(TMath::Abs(pdg_mother) == 22) qadca[1]=4.5;
						//cout<<"Pdg pTy old gamma:   "<<pdg_mother<<endl;}
					/*if(PhotonicType2 == 2){ 
						qadca[1] = 5.5; // gamma
						cout<<"Pdg pTy 2:   "<<pdg_mother<<endl;}
					if(PhotonicType2 == 3){ 
						qadca[1] = 6.5; // dalitz
						cout<<"Pdg pTy 3:   "<<pdg_mother<<endl;}*/
					Int_t fType = GetPi0EtaType(fMCparticleMother,fMCarray);
					//if(fType == 0)cout<<"Feed down from Light mesons:   "<<fType<<"    "<<DCAResCorr<<endl;
					qadca[7]=fType;
                }
            }

     
             
		
	    qadca[2]=DCAxy*track->Charge()*signB;	
            
            qadca[10]=DCAResCorr;
            
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
                //    cout<<"Filling DCA MC================================================================="<<endl;    
            ///Fill
            if(qadca[1]>0.) fD0->Fill(qadca);
        }
        
    //   } // manual PID selection

    }//End of track loop
    
    fNAnalizedTracks->Fill(NAnalizedTracks);
	//fNAnalizedTracksHijing->Fill(NAnalizedTracksHijing);
	
       
    delete fListOfmotherkink;
    PostData(1, fOutputList);
    
}//end of main loop





//=======================================================================
void AliAnalysisHFEppTPCTOFBeauty::Terminate(Option_t *)
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
Bool_t AliAnalysisHFEppTPCTOFBeauty::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
    //Check single track cuts for a given cut step
    //Note this function is called inside the UserExec function
    const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
    if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
    return kTRUE;
}
//=======================================================================
/*
void AliAnalysisHFEppTPCTOFBeauty::GetGammaAndDalitzElectronTemplates(TClonesArray *fMCarray, AliVTrack *track, Double_t fpt, Double_t NewDCA){

       //Float_t dcaxy = -999., dcaz = -999.;
       //fExtraCuts->GetImpactParameters(track, dcaxy, dcaz);

      const AliAODMCParticle  *fMCparticle = (const AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
      //AliHFEmcQA *aodmcqa = dynamic_cast< AliAODmcQA *>(mctrack);
     
      AliHFEmcQA *hfemcqa = new AliHFEmcQA();
      if(hfemcqa) hfemcqa->SetMCArray(fMCarray);
      //fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
      //Int_t PhotonicType = fSignalCuts->GetMCQAObject()->GetElecSource(fMCparticle,kTRUE);
      Int_t PhotonicType = hfemcqa->GetElecSource(fMCparticle, kTRUE);
      if(PhotonicType == 5 || PhotonicType == 8){ // electrons from Dalitz 
      fDCAxy_pt_Dalitz->Fill(fpt, NewDCA);
      }
     
      if(PhotonicType == 31 || PhotonicType == 32 || PhotonicType == 33 || PhotonicType == 34 || PhotonicType == 39){ // electrons from Dalitz from feeddown strange
      fDCAxy_pt_DalitzFromFeeddown->Fill(fpt, NewDCA);
      }
      
      if(PhotonicType == 13 || PhotonicType == 14){ // electrons from photon conversions
      fDCAxy_pt_Conversions->Fill(fpt, NewDCA);
      }
      
      if(PhotonicType == 15 || PhotonicType == 16 || PhotonicType == 17 || PhotonicType == 18){ // electrons from photon conversions from light mesons
      fDCAxy_pt_ConversionsFromFeeddown->Fill(fpt, NewDCA);
      }
      
      if(PhotonicType == 24){ // electrons from photon conversions from light mesons
      fDCAxy_pt_DalitzFromFeeddown2->Fill(fpt, NewDCA);
      }
      
      if(PhotonicType == 35 || PhotonicType == 36 || PhotonicType == 37 || PhotonicType == 38 || PhotonicType == 40){ // electrons from photon conversions from strange
      fDCAxy_pt_ConversionsFromStrangeFeeddown->Fill(fpt, NewDCA);
      }

      
    }*/
/*
void AliAnalysisHFEppTPCTOFBeauty::GetGammaAndDalitzElectronTemplates2(TClonesArray *fMCarray2, AliVTrack *track2, Double_t fpt2, Double_t NewDCA2, Int_t PhotonicType){

       //Float_t dcaxy = -999., dcaz = -999.;
       //fExtraCuts->GetImpactParameters(track, dcaxy, dcaz);

      const AliAODMCParticle  *fMCparticle = (const AliAODMCParticle*) fMCarray2->At(TMath::Abs(track2->GetLabel()));
      //AliHFEmcQA *aodmcqa = dynamic_cast< AliAODmcQA *>(mctrack);
      //AliHFEsignalCuts *fSignalCuts = new AliHFEsignalCuts("HFEsignalCuts", "HFE MC Signal definition");
      
      //fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
      Double_t mcmpt = -999;
      Int_t mcmpdg = 0;
//      Int_t PhotonicType = fSignalCuts->GetSignalSource(track2);
      if(PhotonicType == 3){ // electrons from Dalitz 
      fDCAxy_pt_Dalitz2->Fill(fpt2, NewDCA2);
      }
          
      if(PhotonicType == 2){ // electrons from photon conversions
      fDCAxy_pt_Conversions2->Fill(fpt2, NewDCA2);
      }
      
  
      
    }


void AliAnalysisHFEppTPCTOFBeauty::GetHFElectronTemplates2(TClonesArray *fMCarray2, AliVTrack *track2, Double_t fpt2, Double_t NewDCA2, Int_t PhotonicType){

       //Float_t dcaxy = -999., dcaz = -999.;
       //fExtraCuts->GetImpactParameters(track, dcaxy, dcaz);

      const AliAODMCParticle  *fMCparticle = (const AliAODMCParticle*) fMCarray2->At(TMath::Abs(track2->GetLabel()));
      //AliHFEmcQA *aodmcqa = dynamic_cast< AliAODmcQA *>(mctrack);
      //AliHFEsignalCuts *fSignalCuts = new AliHFEsignalCuts("HFEsignalCuts", "HFE MC Signal definition");
      /* AliHFEsignalCuts *fSignalCuts = new AliHFEsignalCuts("HFEsignalCuts", "HFE MC Signal definition");
      if(fMCarray){
      fSignalCuts->SetMCAODInfo(fMCarray);
      }
      
      AliHFEmcQA *hfemcqa = new AliHFEmcQA();
      if(hfemcqa) hfemcqa->SetMCArray(fMCarray2);*/
      //fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
     /* Double_t mcmpt = -999;
      Int_t mcmpdg = 0;
      //Int_t PhotonicType = fSignalCuts->GetSignalSource(track2);
      
      if(PhotonicType == 0){ // electrons from photon conversions from light mesons
      fDCAxy_pt_Charm2->Fill(fpt2, NewDCA2);
      }
      
      if(PhotonicType == 1){ // electrons from photon conversions from strange
      fDCAxy_pt_Beauty2->Fill(fpt2, NewDCA2);
      }

      
    }

*/

//=======================================================================
//Setter for the PID cuts (TOF and TPC)
void AliAnalysisHFEppTPCTOFBeauty::SetPIDCuts(Float_t tpcPIDmincut, Float_t tpcPIDmaxcut, Float_t tofPIDmincut, Float_t tofPIDmaxcut) {
    ftpcPIDmincut = tpcPIDmincut;
    ftpcPIDmaxcut = tpcPIDmaxcut;
    ftofPIDmincut = tofPIDmincut;
    ftofPIDmaxcut = tofPIDmaxcut;
}
//=======================================================================


//=======================================================================
//Setter for the Eta cut
void AliAnalysisHFEppTPCTOFBeauty::SetEtaCut(Float_t EtaMin, Float_t EtaMax){
    fEtaMin = EtaMin;
    fEtaMax = EtaMax;
}
//=======================================================================


//=======================================================================
Bool_t AliAnalysisHFEppTPCTOFBeauty::FindMother(Int_t mcIndex)
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
            //    cout<<fMCparticle->GetMother()<<"   "<<mpdg<<"    "<<gmpdg<<"    "<<ggmpdg<<"    "<<gggmpdg<<endl;
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

Bool_t AliAnalysisHFEppTPCTOFBeauty::IsHFelectronsMC(AliVTrack *track){
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
Int_t AliAnalysisHFEppTPCTOFBeauty::GetPi0EtaType(AliAODMCParticle *pi0eta, TClonesArray *mcArray){
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
        if(motherpdg == 111 || motherpdg == 221 || motherpdg == 223 || motherpdg == 333 || motherpdg == 331 || motherpdg == 113 || motherpdg == 213){ 

	Int_t pionmotherlabel = mother->GetMother();
	if(pionmotherlabel > 0){
	AliAODMCParticle *pionmother = (AliAODMCParticle*)fMCarray->At(pionmotherlabel);
        Int_t pionmotherpdg = TMath::Abs(pionmother->GetPdgCode());
	//cout<<"PionMotherPDG:     "<<pionmotherpdg<<endl;
	if ( (int(TMath::Abs(pionmotherpdg)/100.)%10) == 5 || (int(TMath::Abs(pionmotherpdg)/1000.)%10) == 5 ) return kBeauty; ///beauty mesons and barions
        if ( (int(TMath::Abs(pionmotherpdg)/100.)%10) == 4 || (int(TMath::Abs(pionmotherpdg)/1000.)%10) == 4 ) return kCharm; ///charmed mesons and barions

	if (pionmotherpdg == 111 || pionmotherpdg == 221 || pionmotherpdg == 223 || pionmotherpdg == 333 || pionmotherpdg == 331 || pionmotherpdg == 113 || pionmotherpdg == 213){ 
	Int_t gpionmotherlabel = pionmother->GetMother();
        if(gpionmotherlabel > 0){
        AliAODMCParticle *gpionmother = (AliAODMCParticle*)fMCarray->At(gpionmotherlabel);
        Int_t gpionmotherpdg = TMath::Abs(gpionmother->GetPdgCode());
	//cout<<"PionGMotherPDG:     "<<gpionmotherpdg<<endl;
	if ( (int(TMath::Abs(gpionmotherpdg)/100.)%10) == 5 || (int(TMath::Abs(gpionmotherpdg)/1000.)%10) == 5 ) return kBeauty; ///beauty mesons and barions
        if ( (int(TMath::Abs(gpionmotherpdg)/100.)%10) == 4 || (int(TMath::Abs(gpionmotherpdg)/1000.)%10) == 4 ) return kCharm; ///charmed mesons and barions
	 if (gpionmotherpdg == 111 || gpionmotherpdg == 221 || gpionmotherpdg == 223 || gpionmotherpdg == 333 || gpionmotherpdg == 331 || gpionmotherpdg == 113 || gpionmotherpdg == 213){
	Int_t ggpionmotherlabel = gpionmother->GetMother();
        if(ggpionmotherlabel > 0){
        AliAODMCParticle *ggpionmother = (AliAODMCParticle*)fMCarray->At(ggpionmotherlabel);
        Int_t ggpionmotherpdg = TMath::Abs(ggpionmother->GetPdgCode());
        //cout<<"PionGGMotherPDG:     "<<ggpionmotherpdg<<endl;
        if ( (int(TMath::Abs(ggpionmotherpdg)/100.)%10) == 5 || (int(TMath::Abs(ggpionmotherpdg)/1000.)%10) == 5 ) return kBeauty; ///beauty mesons and barions
        if ( (int(TMath::Abs(ggpionmotherpdg)/100.)%10) == 4 || (int(TMath::Abs(ggpionmotherpdg)/1000.)%10) == 4 ) return kCharm; ///charmed mesons and barions
	return kNoMother;
	}else return kNoMother;
	}
	return kNoMother;
	}
	else return kNoMother; ///kGammaM2M
	}
	}
	else{
	return kLightMesons;
	}
	return kLightMesons;
}
        
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

Int_t AliAnalysisHFEppTPCTOFBeauty::GetElecSourceType(AliAODMCParticle *electron, TClonesArray *mcArray,Double_t &ptm){
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
Bool_t AliAnalysisHFEppTPCTOFBeauty::GetNMCPartProduced()
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
Int_t AliAnalysisHFEppTPCTOFBeauty::GetPrimary(Int_t id, TClonesArray *mcArray){
    
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
Bool_t AliAnalysisHFEppTPCTOFBeauty::PassCorrCuts(AliAODEvent *fAOD)
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

/*
Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAMeanMC_phi1(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float MeanGaus = 0;

if (x >= 0.90 && x < 1.10) MeanGaus = 0.000109; 
if (x >= 1.10 && x < 1.30) MeanGaus = 0.000086; 
if (x >= 1.30 && x < 1.50) MeanGaus = 0.000126; 
if (x >= 1.50 && x < 1.70) MeanGaus = 0.000155; 
if (x >= 1.70 && x < 1.90) MeanGaus = 0.000173; 
if (x >= 1.90 && x < 2.10) MeanGaus = 0.000231; 
if (x >= 2.10 && x < 2.30) MeanGaus = 0.000254; 
if (x >= 2.30 && x < 2.50) MeanGaus = 0.000261; 
if (x >= 2.50 && x < 2.70) MeanGaus = 0.000310; 
if (x >= 2.70 && x < 2.90) MeanGaus = 0.000329; 
if (x >= 2.90 && x < 3.10) MeanGaus = 0.000383; 
if (x >= 3.10 && x < 3.30) MeanGaus = 0.000344; 
if (x >= 3.30 && x < 3.50) MeanGaus = 0.000367; 
if (x >= 3.50 && x < 3.70) MeanGaus = 0.000382; 
if (x >= 3.70 && x < 3.90) MeanGaus = 0.000365; 
if (x >= 3.90 && x < 4.10) MeanGaus = 0.000374; 
if (x >= 4.10 && x < 4.30) MeanGaus = 0.000394; 
if (x >= 4.30 && x < 4.50) MeanGaus = 0.000420; 
if (x >= 4.50 && x < 4.70) MeanGaus = 0.000411; 
if (x >= 4.70 && x < 5.00) MeanGaus = 0.000426; 
if (x >= 5.00 && x < 5.50) MeanGaus = 0.000371; 
if (x >= 5.50 && x < 6.00) MeanGaus = 0.000495; 
if (x >= 6.00 && x < 6.50) MeanGaus = 0.000400; 
if (x >= 6.50 && x < 7.00) MeanGaus = 0.000422; 
if (x >= 7.00 && x < 8.00) MeanGaus = 0.000379; 

return MeanGaus;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAMeanMC_phi2(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float MeanGaus = 0;

if (x >= 0.90 && x < 1.10) MeanGaus = -0.001688; 
if (x >= 1.10 && x < 1.30) MeanGaus = -0.001688; 
if (x >= 1.30 && x < 1.50) MeanGaus = -0.001554; 
if (x >= 1.50 && x < 1.70) MeanGaus = -0.001371; 
if (x >= 1.70 && x < 1.90) MeanGaus = -0.001220; 
if (x >= 1.90 && x < 2.10) MeanGaus = -0.001140; 
if (x >= 2.10 && x < 2.30) MeanGaus = -0.001087; 
if (x >= 2.30 && x < 2.50) MeanGaus = -0.001022; 
if (x >= 2.50 && x < 2.70) MeanGaus = -0.000985; 
if (x >= 2.70 && x < 2.90) MeanGaus = -0.000883; 
if (x >= 2.90 && x < 3.10) MeanGaus = -0.000879; 
if (x >= 3.10 && x < 3.30) MeanGaus = -0.000844; 
if (x >= 3.30 && x < 3.50) MeanGaus = -0.000769; 
if (x >= 3.50 && x < 3.70) MeanGaus = -0.000690; 
if (x >= 3.70 && x < 3.90) MeanGaus = -0.000732; 
if (x >= 3.90 && x < 4.10) MeanGaus = -0.000695; 
if (x >= 4.10 && x < 4.30) MeanGaus = -0.000652; 
if (x >= 4.30 && x < 4.50) MeanGaus = -0.000640; 
if (x >= 4.50 && x < 4.70) MeanGaus = -0.000600; 
if (x >= 4.70 && x < 5.00) MeanGaus = -0.000500; 
if (x >= 5.00 && x < 5.50) MeanGaus = -0.000483; 
if (x >= 5.50 && x < 6.00) MeanGaus = -0.000488; 
if (x >= 6.00 && x < 6.50) MeanGaus = -0.000383; 
if (x >= 6.50 && x < 7.00) MeanGaus = -0.000447; 
if (x >= 7.00 && x < 8.00) MeanGaus = -0.000396;   

return MeanGaus;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAMeanMC_phi3(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float MeanGaus = 0;

if (x >= 0.90 && x < 1.10) MeanGaus = -0.000696; 
if (x >= 1.10 && x < 1.30) MeanGaus = -0.000704; 
if (x >= 1.30 && x < 1.50) MeanGaus = -0.000601; 
if (x >= 1.50 && x < 1.70) MeanGaus = -0.000553; 
if (x >= 1.70 && x < 1.90) MeanGaus = -0.000481; 
if (x >= 1.90 && x < 2.10) MeanGaus = -0.000453; 
if (x >= 2.10 && x < 2.30) MeanGaus = -0.000442; 
if (x >= 2.30 && x < 2.50) MeanGaus = -0.000422; 
if (x >= 2.50 && x < 2.70) MeanGaus = -0.000403; 
if (x >= 2.70 && x < 2.90) MeanGaus = -0.000317; 
if (x >= 2.90 && x < 3.10) MeanGaus = -0.000334; 
if (x >= 3.10 && x < 3.30) MeanGaus = -0.000247; 
if (x >= 3.30 && x < 3.50) MeanGaus = -0.000291; 
if (x >= 3.50 && x < 3.70) MeanGaus = -0.000297; 
if (x >= 3.70 && x < 3.90) MeanGaus = -0.000216; 
if (x >= 3.90 && x < 4.10) MeanGaus = -0.000255; 
if (x >= 4.10 && x < 4.30) MeanGaus = -0.000140; 
if (x >= 4.30 && x < 4.50) MeanGaus = -0.000185; 
if (x >= 4.50 && x < 4.70) MeanGaus = -0.000150; 
if (x >= 4.70 && x < 5.00) MeanGaus = -0.000205; 
if (x >= 5.00 && x < 5.50) MeanGaus = -0.000106; 
if (x >= 5.50 && x < 6.00) MeanGaus = -0.000172; 
if (x >= 6.00 && x < 6.50) MeanGaus = -0.000188; 
if (x >= 6.50 && x < 7.00) MeanGaus = -0.000045; 
if (x >= 7.00 && x < 8.00) MeanGaus = -0.000096; 
return MeanGaus;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAMeanMC_phi4(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float MeanGaus = 0;

if (x >= 0.90 && x < 1.10) MeanGaus = -0.001975; 
if (x >= 1.10 && x < 1.30) MeanGaus = -0.001952; 
if (x >= 1.30 && x < 1.50) MeanGaus = -0.001801; 
if (x >= 1.50 && x < 1.70) MeanGaus = -0.001590; 
if (x >= 1.70 && x < 1.90) MeanGaus = -0.001448; 
if (x >= 1.90 && x < 2.10) MeanGaus = -0.001278; 
if (x >= 2.10 && x < 2.30) MeanGaus = -0.001204; 
if (x >= 2.30 && x < 2.50) MeanGaus = -0.001104; 
if (x >= 2.50 && x < 2.70) MeanGaus = -0.001036; 
if (x >= 2.70 && x < 2.90) MeanGaus = -0.000929; 
if (x >= 2.90 && x < 3.10) MeanGaus = -0.000852; 
if (x >= 3.10 && x < 3.30) MeanGaus = -0.000819; 
if (x >= 3.30 && x < 3.50) MeanGaus = -0.000760; 
if (x >= 3.50 && x < 3.70) MeanGaus = -0.000689; 
if (x >= 3.70 && x < 3.90) MeanGaus = -0.000655; 
if (x >= 3.90 && x < 4.10) MeanGaus = -0.000534; 
if (x >= 4.10 && x < 4.30) MeanGaus = -0.000568; 
if (x >= 4.30 && x < 4.50) MeanGaus = -0.000500; 
if (x >= 4.50 && x < 4.70) MeanGaus = -0.000443; 
if (x >= 4.70 && x < 5.00) MeanGaus = -0.000412; 
if (x >= 5.00 && x < 5.50) MeanGaus = -0.000367; 
if (x >= 5.50 && x < 6.00) MeanGaus = -0.000361; 
if (x >= 6.00 && x < 6.50) MeanGaus = -0.000272; 
if (x >= 6.50 && x < 7.00) MeanGaus = -0.000374; 
if (x >= 7.00 && x < 8.00) MeanGaus = -0.000201; 

return MeanGaus;   
}


Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAResolMC_phi1(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float sigmaG = 0;

if (x >= 0.90 && x < 1.10) sigmaG = 0.003523; 
if (x >= 1.10 && x < 1.30) sigmaG = 0.003597; 
if (x >= 1.30 && x < 1.50) sigmaG = 0.003269; 
if (x >= 1.50 && x < 1.70) sigmaG = 0.003005; 
if (x >= 1.70 && x < 1.90) sigmaG = 0.002800; 
if (x >= 1.90 && x < 2.10) sigmaG = 0.002715; 
if (x >= 2.10 && x < 2.30) sigmaG = 0.002586; 
if (x >= 2.30 && x < 2.50) sigmaG = 0.002534; 
if (x >= 2.50 && x < 2.70) sigmaG = 0.002456; 
if (x >= 2.70 && x < 2.90) sigmaG = 0.002425; 
if (x >= 2.90 && x < 3.10) sigmaG = 0.002392; 
if (x >= 3.10 && x < 3.30) sigmaG = 0.002368; 
if (x >= 3.30 && x < 3.50) sigmaG = 0.002203; 
if (x >= 3.50 && x < 3.70) sigmaG = 0.002209; 
if (x >= 3.70 && x < 3.90) sigmaG = 0.002102; 
if (x >= 3.90 && x < 4.10) sigmaG = 0.002033; 
if (x >= 4.10 && x < 4.30) sigmaG = 0.002168; 
if (x >= 4.30 && x < 4.50) sigmaG = 0.002091; 
if (x >= 4.50 && x < 4.70) sigmaG = 0.001998; 
if (x >= 4.70 && x < 5.00) sigmaG = 0.001942; 
if (x >= 5.00 && x < 5.50) sigmaG = 0.001864; 
if (x >= 5.50 && x < 6.00) sigmaG = 0.001783; 
if (x >= 6.00 && x < 6.50) sigmaG = 0.001491; 
if (x >= 6.50 && x < 7.00) sigmaG = 0.001569; 
if (x >= 7.00 && x < 8.00) sigmaG = 0.001591;

return sigmaG;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAResolMC_phi2(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float sigmaG = 0;

if (x >= 0.90 && x < 1.10) sigmaG = 0.004456; 
if (x >= 1.10 && x < 1.30) sigmaG = 0.004643; 
if (x >= 1.30 && x < 1.50) sigmaG = 0.004339; 
if (x >= 1.50 && x < 1.70) sigmaG = 0.004005; 
if (x >= 1.70 && x < 1.90) sigmaG = 0.003767; 
if (x >= 1.90 && x < 2.10) sigmaG = 0.003611; 
if (x >= 2.10 && x < 2.30) sigmaG = 0.003501; 
if (x >= 2.30 && x < 2.50) sigmaG = 0.003461; 
if (x >= 2.50 && x < 2.70) sigmaG = 0.003427; 
if (x >= 2.70 && x < 2.90) sigmaG = 0.003316; 
if (x >= 2.90 && x < 3.10) sigmaG = 0.003208; 
if (x >= 3.10 && x < 3.30) sigmaG = 0.003076; 
if (x >= 3.30 && x < 3.50) sigmaG = 0.003071; 
if (x >= 3.50 && x < 3.70) sigmaG = 0.003018; 
if (x >= 3.70 && x < 3.90) sigmaG = 0.002872; 
if (x >= 3.90 && x < 4.10) sigmaG = 0.002829; 
if (x >= 4.10 && x < 4.30) sigmaG = 0.002774; 
if (x >= 4.30 && x < 4.50) sigmaG = 0.002707; 
if (x >= 4.50 && x < 4.70) sigmaG = 0.002721; 
if (x >= 4.70 && x < 5.00) sigmaG = 0.002718; 
if (x >= 5.00 && x < 5.50) sigmaG = 0.002518; 
if (x >= 5.50 && x < 6.00) sigmaG = 0.002492; 
if (x >= 6.00 && x < 6.50) sigmaG = 0.002546; 
if (x >= 6.50 && x < 7.00) sigmaG = 0.002241; 
if (x >= 7.00 && x < 8.00) sigmaG = 0.001996; 

return sigmaG;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAResolMC_phi3(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float sigmaG = 0;

if (x >= 0.90 && x < 1.10) sigmaG = 0.004687; 
if (x >= 1.10 && x < 1.30) sigmaG = 0.004552; 
if (x >= 1.30 && x < 1.50) sigmaG = 0.004081; 
if (x >= 1.50 && x < 1.70) sigmaG = 0.003625; 
if (x >= 1.70 && x < 1.90) sigmaG = 0.003370; 
if (x >= 1.90 && x < 2.10) sigmaG = 0.003223; 
if (x >= 2.10 && x < 2.30) sigmaG = 0.003101; 
if (x >= 2.30 && x < 2.50) sigmaG = 0.003038; 
if (x >= 2.50 && x < 2.70) sigmaG = 0.002978; 
if (x >= 2.70 && x < 2.90) sigmaG = 0.002852; 
if (x >= 2.90 && x < 3.10) sigmaG = 0.002712; 
if (x >= 3.10 && x < 3.30) sigmaG = 0.002661; 
if (x >= 3.30 && x < 3.50) sigmaG = 0.002737; 
if (x >= 3.50 && x < 3.70) sigmaG = 0.002538; 
if (x >= 3.70 && x < 3.90) sigmaG = 0.002447; 
if (x >= 3.90 && x < 4.10) sigmaG = 0.002486; 
if (x >= 4.10 && x < 4.30) sigmaG = 0.002483; 
if (x >= 4.30 && x < 4.50) sigmaG = 0.002566; 
if (x >= 4.50 && x < 4.70) sigmaG = 0.002467; 
if (x >= 4.70 && x < 5.00) sigmaG = 0.002358; 
if (x >= 5.00 && x < 5.50) sigmaG = 0.002283; 
if (x >= 5.50 && x < 6.00) sigmaG = 0.002238; 
if (x >= 6.00 && x < 6.50) sigmaG = 0.002261; 
if (x >= 6.50 && x < 7.00) sigmaG = 0.002097; 
if (x >= 7.00 && x < 8.00) sigmaG = 0.001937;   
return sigmaG;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAResolMC_phi4(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float sigmaG = 0;

if (x >= 0.90 && x < 1.10) sigmaG = 0.004601; 
if (x >= 1.10 && x < 1.30) sigmaG = 0.004646; 
if (x >= 1.30 && x < 1.50) sigmaG = 0.004214; 
if (x >= 1.50 && x < 1.70) sigmaG = 0.003825; 
if (x >= 1.70 && x < 1.90) sigmaG = 0.003556; 
if (x >= 1.90 && x < 2.10) sigmaG = 0.003325; 
if (x >= 2.10 && x < 2.30) sigmaG = 0.003158; 
if (x >= 2.30 && x < 2.50) sigmaG = 0.003140; 
if (x >= 2.50 && x < 2.70) sigmaG = 0.003017; 
if (x >= 2.70 && x < 2.90) sigmaG = 0.002917; 
if (x >= 2.90 && x < 3.10) sigmaG = 0.002934; 
if (x >= 3.10 && x < 3.30) sigmaG = 0.002827; 
if (x >= 3.30 && x < 3.50) sigmaG = 0.002844; 
if (x >= 3.50 && x < 3.70) sigmaG = 0.002636; 
if (x >= 3.70 && x < 3.90) sigmaG = 0.002595; 
if (x >= 3.90 && x < 4.10) sigmaG = 0.002581; 
if (x >= 4.10 && x < 4.30) sigmaG = 0.002691; 
if (x >= 4.30 && x < 4.50) sigmaG = 0.002451; 
if (x >= 4.50 && x < 4.70) sigmaG = 0.002302; 
if (x >= 4.70 && x < 5.00) sigmaG = 0.002417; 
if (x >= 5.00 && x < 5.50) sigmaG = 0.002371; 
if (x >= 5.50 && x < 6.00) sigmaG = 0.002210; 
if (x >= 6.00 && x < 6.50) sigmaG = 0.002085; 
if (x >= 6.50 && x < 7.00) sigmaG = 0.002134; 
if (x >= 7.00 && x < 8.00) sigmaG = 0.002009; 

return sigmaG;   
}*/

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAMeanMC_phi1(Float_t x){

// Return the DCA resolution of the track (in MC) accordingly to its pT

float MeanGaus = 0;

if (x >= 0.90 && x < 1.10) MeanGaus = 0.000086; 
if (x >= 1.10 && x < 1.30) MeanGaus = 0.000070; 
if (x >= 1.30 && x < 1.50) MeanGaus = 0.000090; 
if (x >= 1.50 && x < 1.70) MeanGaus = 0.000101; 
if (x >= 1.70 && x < 1.90) MeanGaus = 0.000118; 
if (x >= 1.90 && x < 2.10) MeanGaus = 0.000134; 
if (x >= 2.10 && x < 2.30) MeanGaus = 0.000142; 
if (x >= 2.30 && x < 2.50) MeanGaus = 0.000154; 
if (x >= 2.50 && x < 2.70) MeanGaus = 0.000172; 
if (x >= 2.70 && x < 2.90) MeanGaus = 0.000170; 
if (x >= 2.90 && x < 3.10) MeanGaus = 0.000180; 
if (x >= 3.10 && x < 3.30) MeanGaus = 0.000178; 
if (x >= 3.30 && x < 3.50) MeanGaus = 0.000183; 
if (x >= 3.50 && x < 3.70) MeanGaus = 0.000182; 
if (x >= 3.70 && x < 3.90) MeanGaus = 0.000173; 
if (x >= 3.90 && x < 4.10) MeanGaus = 0.000166; 
if (x >= 4.10 && x < 4.30) MeanGaus = 0.000178; 
if (x >= 4.30 && x < 4.50) MeanGaus = 0.000175; 
if (x >= 4.50 && x < 4.70) MeanGaus = 0.000166; 
if (x >= 4.70 && x < 5.00) MeanGaus = 0.000179; 
if (x >= 5.00 && x < 5.50) MeanGaus = 0.000155; 
if (x >= 5.50 && x < 6.00) MeanGaus = 0.000170; 
if (x >= 6.00 && x < 6.50) MeanGaus = 0.000148; 
if (x >= 6.50 && x < 7.00) MeanGaus = 0.000143; 
if (x >= 7.00 && x < 8.00) MeanGaus = 0.000137; 

return MeanGaus;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAMeanMC_phi2(Float_t x){

// Return the DCA resolution of the track (in MC) accordingly to its pT

float MeanGaus = 0;

if (x >= 0.90 && x < 1.10) MeanGaus = -0.001080; 
if (x >= 1.10 && x < 1.30) MeanGaus = -0.000994; 
if (x >= 1.30 && x < 1.50) MeanGaus = -0.000854; 
if (x >= 1.50 && x < 1.70) MeanGaus = -0.000714; 
if (x >= 1.70 && x < 1.90) MeanGaus = -0.000604; 
if (x >= 1.90 && x < 2.10) MeanGaus = -0.000551; 
if (x >= 2.10 && x < 2.30) MeanGaus = -0.000515; 
if (x >= 2.30 && x < 2.50) MeanGaus = -0.000482; 
if (x >= 2.50 && x < 2.70) MeanGaus = -0.000459; 
if (x >= 2.70 && x < 2.90) MeanGaus = -0.000401; 
if (x >= 2.90 && x < 3.10) MeanGaus = -0.000377; 
if (x >= 3.10 && x < 3.30) MeanGaus = -0.000347; 
if (x >= 3.30 && x < 3.50) MeanGaus = -0.000326; 
if (x >= 3.50 && x < 3.70) MeanGaus = -0.000284; 
if (x >= 3.70 && x < 3.90) MeanGaus = -0.000289; 
if (x >= 3.90 && x < 4.10) MeanGaus = -0.000258; 
if (x >= 4.10 && x < 4.30) MeanGaus = -0.000234; 
if (x >= 4.30 && x < 4.50) MeanGaus = -0.000224; 
if (x >= 4.50 && x < 4.70) MeanGaus = -0.000201; 
if (x >= 4.70 && x < 5.00) MeanGaus = -0.000170; 
if (x >= 5.00 && x < 5.50) MeanGaus = -0.000169; 
if (x >= 5.50 && x < 6.00) MeanGaus = -0.000141; 
if (x >= 6.00 && x < 6.50) MeanGaus = -0.000133; 
if (x >= 6.50 && x < 7.00) MeanGaus = -0.000111; 
if (x >= 7.00 && x < 8.00) MeanGaus = -0.000124;  

return MeanGaus;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAMeanMC_phi3(Float_t x){

// Return the DCA resolution of the track (in MC) accordingly to its pT

float MeanGaus = 0;

if (x >= 0.90 && x < 1.10) MeanGaus = -0.000504; 
if (x >= 1.10 && x < 1.30) MeanGaus = -0.000454; 
if (x >= 1.30 && x < 1.50) MeanGaus = -0.000375; 
if (x >= 1.50 && x < 1.70) MeanGaus = -0.000307; 
if (x >= 1.70 && x < 1.90) MeanGaus = -0.000273; 
if (x >= 1.90 && x < 2.10) MeanGaus = -0.000263; 
if (x >= 2.10 && x < 2.30) MeanGaus = -0.000240; 
if (x >= 2.30 && x < 2.50) MeanGaus = -0.000202; 
if (x >= 2.50 && x < 2.70) MeanGaus = -0.000195; 
if (x >= 2.70 && x < 2.90) MeanGaus = -0.000170; 
if (x >= 2.90 && x < 3.10) MeanGaus = -0.000147; 
if (x >= 3.10 && x < 3.30) MeanGaus = -0.000092; 
if (x >= 3.30 && x < 3.50) MeanGaus = -0.000106; 
if (x >= 3.50 && x < 3.70) MeanGaus = -0.000088; 
if (x >= 3.70 && x < 3.90) MeanGaus = -0.000082; 
if (x >= 3.90 && x < 4.10) MeanGaus = -0.000075; 
if (x >= 4.10 && x < 4.30) MeanGaus = -0.000032; 
if (x >= 4.30 && x < 4.50) MeanGaus = -0.000046; 
if (x >= 4.50 && x < 4.70) MeanGaus = -0.000039; 
if (x >= 4.70 && x < 5.00) MeanGaus = -0.000050; 
if (x >= 5.00 && x < 5.50) MeanGaus = -0.000015; 
if (x >= 5.50 && x < 6.00) MeanGaus = -0.000029; 
if (x >= 6.00 && x < 6.50) MeanGaus = -0.000001; 
if (x >= 6.50 && x < 7.00) MeanGaus = -0.000016; 
if (x >= 7.00 && x < 8.00) MeanGaus = -0.000032; 
return MeanGaus;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAMeanMC_phi4(Float_t x){

// Return the DCA resolution of the track (in MC) accordingly to its pT

float MeanGaus = 0;

if (x >= 0.90 && x < 1.10) MeanGaus = -0.001437; 
if (x >= 1.10 && x < 1.30) MeanGaus = -0.001372; 
if (x >= 1.30 && x < 1.50) MeanGaus = -0.001207; 
if (x >= 1.50 && x < 1.70) MeanGaus = -0.000999; 
if (x >= 1.70 && x < 1.90) MeanGaus = -0.000844; 
if (x >= 1.90 && x < 2.10) MeanGaus = -0.000724; 
if (x >= 2.10 && x < 2.30) MeanGaus = -0.000658; 
if (x >= 2.30 && x < 2.50) MeanGaus = -0.000586; 
if (x >= 2.50 && x < 2.70) MeanGaus = -0.000513; 
if (x >= 2.70 && x < 2.90) MeanGaus = -0.000466; 
if (x >= 2.90 && x < 3.10) MeanGaus = -0.000410; 
if (x >= 3.10 && x < 3.30) MeanGaus = -0.000359; 
if (x >= 3.30 && x < 3.50) MeanGaus = -0.000338; 
if (x >= 3.50 && x < 3.70) MeanGaus = -0.000298; 
if (x >= 3.70 && x < 3.90) MeanGaus = -0.000259; 
if (x >= 3.90 && x < 4.10) MeanGaus = -0.000212; 
if (x >= 4.10 && x < 4.30) MeanGaus = -0.000206; 
if (x >= 4.30 && x < 4.50) MeanGaus = -0.000182; 
if (x >= 4.50 && x < 4.70) MeanGaus = -0.000171; 
if (x >= 4.70 && x < 5.00) MeanGaus = -0.000147; 
if (x >= 5.00 && x < 5.50) MeanGaus = -0.000115; 
if (x >= 5.50 && x < 6.00) MeanGaus = -0.000088; 
if (x >= 6.00 && x < 6.50) MeanGaus = -0.000064; 
if (x >= 6.50 && x < 7.00) MeanGaus = -0.000067; 
if (x >= 7.00 && x < 8.00) MeanGaus = -0.000058;  

return MeanGaus;   
}


Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAResolMC_phi1(Float_t x){

// Return the DCA resolution of the track (in MC) accordingly to its pT

float sigmaG = 0;

if (x >= 0.90 && x < 1.10) sigmaG = 0.002169; 
if (x >= 1.10 && x < 1.30) sigmaG = 0.000993; 
if (x >= 1.30 && x < 1.50) sigmaG = 0.000909; 
if (x >= 1.50 && x < 1.70) sigmaG = 0.000698; 
if (x >= 1.70 && x < 1.90) sigmaG = 0.000479; 
if (x >= 1.90 && x < 2.10) sigmaG = 0.000614; 
if (x >= 2.10 && x < 2.30) sigmaG = 0.000672; 
if (x >= 2.30 && x < 2.50) sigmaG = 0.000599; 
if (x >= 2.50 && x < 2.70) sigmaG = 0.000521; 
if (x >= 2.70 && x < 2.90) sigmaG = 0.000383; 
if (x >= 2.90 && x < 3.10) sigmaG = 0.000306; 
if (x >= 3.10 && x < 3.30) sigmaG = 0.000214; 
if (x >= 3.30 && x < 3.50) sigmaG = 0.000330; 
if (x >= 3.50 && x < 3.70) sigmaG = 0.000387; 
if (x >= 3.70 && x < 3.90) sigmaG = 0.000202; 
if (x >= 3.90 && x < 4.10) sigmaG = 0.000241; 
if (x >= 4.10 && x < 4.30) sigmaG = 0.000062; 
if (x >= 4.30 && x < 4.50) sigmaG = 0.000206; 
if (x >= 4.50 && x < 4.70) sigmaG = 0.000454; 
if (x >= 4.70 && x < 5.00) sigmaG = 0.000643; 
if (x >= 5.00 && x < 5.50) sigmaG = 0.000320; 
if (x >= 5.50 && x < 6.00) sigmaG = 0.000525; 
if (x >= 6.00 && x < 6.50) sigmaG = 0.000697; 
if (x >= 6.50 && x < 7.00) sigmaG = 0.000151; 
if (x >= 7.00 && x < 8.00) sigmaG = 0.000526; 

return sigmaG;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAResolMC_phi2(Float_t x){

// Return the DCA resolution of the track (in MC) accordingly to its pT

float sigmaG = 0;

if (x >= 0.90 && x < 1.10) sigmaG = 0.001767; 
if (x >= 1.10 && x < 1.30) sigmaG = 0.001370; 
if (x >= 1.30 && x < 1.50) sigmaG = 0.001260; 
if (x >= 1.50 && x < 1.70) sigmaG = 0.001041; 
if (x >= 1.70 && x < 1.90) sigmaG = 0.000906; 
if (x >= 1.90 && x < 2.10) sigmaG = 0.001000; 
if (x >= 2.10 && x < 2.30) sigmaG = 0.000963; 
if (x >= 2.30 && x < 2.50) sigmaG = 0.000955; 
if (x >= 2.50 && x < 2.70) sigmaG = 0.000847; 
if (x >= 2.70 && x < 2.90) sigmaG = 0.000688; 
if (x >= 2.90 && x < 3.10) sigmaG = 0.000564; 
if (x >= 3.10 && x < 3.30) sigmaG = 0.000448; 
if (x >= 3.30 && x < 3.50) sigmaG = 0.000489; 
if (x >= 3.50 && x < 3.70) sigmaG = 0.000597; 
if (x >= 3.70 && x < 3.90) sigmaG = 0.000390; 
if (x >= 3.90 && x < 4.10) sigmaG = 0.000454; 
if (x >= 4.10 && x < 4.30) sigmaG = 0.000349; 
if (x >= 4.30 && x < 4.50) sigmaG = 0.000417; 
if (x >= 4.50 && x < 4.70) sigmaG = 0.000318; 
if (x >= 4.70 && x < 5.00) sigmaG = 0.000365; 
if (x >= 5.00 && x < 5.50) sigmaG = 0.000345; 
if (x >= 5.50 && x < 6.00) sigmaG = 0.000101; 
if (x >= 6.00 && x < 6.50) sigmaG = 0.000156; 
if (x >= 6.50 && x < 7.00) sigmaG = 0.000423; 
if (x >= 7.00 && x < 8.00) sigmaG = 0.000571; 

return sigmaG;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAResolMC_phi3(Float_t x){

// Return the DCA resolution of the track (in MC) accordingly to its pT

float sigmaG = 0;

if (x >= 0.90 && x < 1.10) sigmaG = 0.002756; 
if (x >= 1.10 && x < 1.30) sigmaG = 0.001274; 
if (x >= 1.30 && x < 1.50) sigmaG = 0.001667; 
if (x >= 1.50 && x < 1.70) sigmaG = 0.001323; 
if (x >= 1.70 && x < 1.90) sigmaG = 0.001157; 
if (x >= 1.90 && x < 2.10) sigmaG = 0.001088; 
if (x >= 2.10 && x < 2.30) sigmaG = 0.001124; 
if (x >= 2.30 && x < 2.50) sigmaG = 0.001120; 
if (x >= 2.50 && x < 2.70) sigmaG = 0.001058; 
if (x >= 2.70 && x < 2.90) sigmaG = 0.000790; 
if (x >= 2.90 && x < 3.10) sigmaG = 0.000799; 
if (x >= 3.10 && x < 3.30) sigmaG = 0.000727; 
if (x >= 3.30 && x < 3.50) sigmaG = 0.000719; 
if (x >= 3.50 && x < 3.70) sigmaG = 0.000575; 
if (x >= 3.70 && x < 3.90) sigmaG = 0.000660; 
if (x >= 3.90 && x < 4.10) sigmaG = 0.000380; 
if (x >= 4.10 && x < 4.30) sigmaG = 0.000623; 
if (x >= 4.30 && x < 4.50) sigmaG = 0.000462; 
if (x >= 4.50 && x < 4.70) sigmaG = 0.000401; 
if (x >= 4.70 && x < 5.00) sigmaG = 0.000581; 
if (x >= 5.00 && x < 5.50) sigmaG = 0.000410; 
if (x >= 5.50 && x < 6.00) sigmaG = 0.000476; 
if (x >= 6.00 && x < 6.50) sigmaG = 0.000257; 
if (x >= 6.50 && x < 7.00) sigmaG = 0.000387; 
if (x >= 7.00 && x < 8.00) sigmaG = 0.000279;  
return sigmaG;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAResolMC_phi4(Float_t x){

// Return the DCA resolution of the track (in MC) accordingly to its pT

float sigmaG = 0;

if (x >= 0.90 && x < 1.10) sigmaG = 0.002412; 
if (x >= 1.10 && x < 1.30) sigmaG = 0.001436; 
if (x >= 1.30 && x < 1.50) sigmaG = 0.001695; 
if (x >= 1.50 && x < 1.70) sigmaG = 0.001356; 
if (x >= 1.70 && x < 1.90) sigmaG = 0.001193; 
if (x >= 1.90 && x < 2.10) sigmaG = 0.001136; 
if (x >= 2.10 && x < 2.30) sigmaG = 0.001075; 
if (x >= 2.30 && x < 2.50) sigmaG = 0.001062; 
if (x >= 2.50 && x < 2.70) sigmaG = 0.000939; 
if (x >= 2.70 && x < 2.90) sigmaG = 0.000751; 
if (x >= 2.90 && x < 3.10) sigmaG = 0.000605; 
if (x >= 3.10 && x < 3.30) sigmaG = 0.000590; 
if (x >= 3.30 && x < 3.50) sigmaG = 0.000704; 
if (x >= 3.50 && x < 3.70) sigmaG = 0.000419; 
if (x >= 3.70 && x < 3.90) sigmaG = 0.000517; 
if (x >= 3.90 && x < 4.10) sigmaG = 0.000598; 
if (x >= 4.10 && x < 4.30) sigmaG = 0.000677; 
if (x >= 4.30 && x < 4.50) sigmaG = 0.000674; 
if (x >= 4.50 && x < 4.70) sigmaG = 0.000508; 
if (x >= 4.70 && x < 5.00) sigmaG = 0.000537; 
if (x >= 5.00 && x < 5.50) sigmaG = 0.000612; 
if (x >= 5.50 && x < 6.00) sigmaG = 0.000579; 
if (x >= 6.00 && x < 6.50) sigmaG = 0.000588; 
if (x >= 6.50 && x < 7.00) sigmaG = 0.000406; 
if (x >= 7.00 && x < 8.00) sigmaG = 0.000467;  

return sigmaG;   
}



//=======================================================================





