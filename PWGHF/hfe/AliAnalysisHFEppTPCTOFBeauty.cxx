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
,fBHadpT(0)
,fBMesonpT(0)
,fBMesonpT_Corr(0)
//,fBMesonpTG(0)
//,fBMesonpTGG(0)
,fDHadpT(0)
,fDMesonpT(0)
,fBDHadpT(0)
,fD0pT(0)
,fLambdaCpT(0)
,fLambdaCpT_AFCorr(0)
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
 ,fDCAxy_pt_beautybaryons(0)
 ,fDCAxy_pt_beautybaryons_corr(0)
,fDCAxy_pt_DstarDplusbef(0)
,fDCAxy_pt_Dplusbef(0)
,fDCAxy_pt_DstarDzerobef(0)
,fDCAxy_pt_Dzerobef(0)
,fDCAxy_pt_DstarDsbef(0)
,fDCAxy_pt_Dsbef(0)
,fDCAxy_pt_charmmesonsbef(0)
,fDCAxy_pt_DstarDplusAft(0)
,fDCAxy_pt_DplusAft(0)
,fDCAxy_pt_DstarDzeroAft(0)
,fDCAxy_pt_DzeroAft(0)
,fDCAxy_pt_DstarDsAft(0)
,fDCAxy_pt_DsAft(0)
,fDCAxy_pt_charmmesonsAft(0)
,fDCAxy_pt_Lc(0)
,fDCAxy_pt_charmbaryons(0)
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
,fDCAxy_pt_ele_onlyDCA(0)
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
,fDCAxy_pt_ele_onlyDCA_phi1(0)
,fDCAxy_pt_ele_onlyDCA_phi2(0)
,fDCAxy_pt_ele_onlyDCA_phi3(0)
,fDCAxy_pt_ele_onlyDCA_phi4(0)
,fDCAxy_pt_ele_phi1_ChB(0)
,fDCAxy_pt_ele_phi1_B(0)
,fDCAxy_pt_ele_phi2_ChB(0)
,fDCAxy_pt_ele_phi2_B(0)
,fDCAxy_pt_ele_phi3_ChB(0)
,fDCAxy_pt_ele_phi3_B(0)
,fDCAxy_pt_ele_phi4_ChB(0)
,fDCAxy_pt_ele_phi4_B(0)
,fDCAxy_pt_ele_ResCorr_phi1(0)
,fDCAxy_pt_ele_ResCorr_phi2(0)
,fDCAxy_pt_ele_ResCorr_phi3(0)
,fDCAxy_pt_ele_ResCorr_phi4(0)
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
,hPtDstarD0(0)
,hPtDstarDp(0)
,hPtDp(0)
,hPtDs(0)
,hPtD0_AFCorr(0)
,hPtDstarD0_AFCorr(0)
,hPtDp_AFCorr(0)
,hPtDs_AFCorr(0)
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
,HistWeightD0(0)
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
,fKPicorr(0)
,fNpureMC(0)
,fNembMCpi0(0)
,fNembMCeta(0)
,fNTotMCpart(0)
,fBcorr(0)
,fDcorr(0)
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
,fCalculateBDMesonpTWeights(0)
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
,fBHadpT(0)
,fBMesonpT(0)
,fBMesonpT_Corr(0)
//,fBMesonpTG(0)
//,fBMesonpTGG(0)
,fDHadpT(0)
,fDMesonpT(0)
,fBDHadpT(0)
,fD0pT(0)
,fLambdaCpT(0)
,fLambdaCpT_AFCorr(0)
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
 ,fDCAxy_pt_beautybaryons(0)
 ,fDCAxy_pt_beautybaryons_corr(0)
,fDCAxy_pt_DstarDplusbef(0)
,fDCAxy_pt_Dplusbef(0)
,fDCAxy_pt_DstarDzerobef(0)
,fDCAxy_pt_Dzerobef(0)
,fDCAxy_pt_DstarDsbef(0)
,fDCAxy_pt_Dsbef(0)
,fDCAxy_pt_charmmesonsbef(0)
,fDCAxy_pt_DstarDplusAft(0)
,fDCAxy_pt_DplusAft(0)
,fDCAxy_pt_DstarDzeroAft(0)
,fDCAxy_pt_DzeroAft(0)
,fDCAxy_pt_DstarDsAft(0)
,fDCAxy_pt_DsAft(0)
,fDCAxy_pt_charmmesonsAft(0)
,fDCAxy_pt_Lc(0)
,fDCAxy_pt_charmbaryons(0)
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
,fDCAxy_pt_ele_onlyDCA(0)
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
,fDCAxy_pt_ele_onlyDCA_phi1(0)
,fDCAxy_pt_ele_onlyDCA_phi2(0)
,fDCAxy_pt_ele_onlyDCA_phi3(0)
,fDCAxy_pt_ele_onlyDCA_phi4(0)
,fDCAxy_pt_ele_phi1_ChB(0)
,fDCAxy_pt_ele_phi1_B(0)
,fDCAxy_pt_ele_phi2_ChB(0)
,fDCAxy_pt_ele_phi2_B(0)
,fDCAxy_pt_ele_phi3_ChB(0)
,fDCAxy_pt_ele_phi3_B(0)
,fDCAxy_pt_ele_phi4_ChB(0)
,fDCAxy_pt_ele_phi4_B(0)
,fDCAxy_pt_ele_ResCorr_phi1(0)
,fDCAxy_pt_ele_ResCorr_phi2(0)
,fDCAxy_pt_ele_ResCorr_phi3(0)
,fDCAxy_pt_ele_ResCorr_phi4(0)
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
,hPtDstarD0(0)
,hPtDstarDp(0)
,hPtDp(0)
,hPtDs(0)
,fKPicorr(0)
,hPtD0_AFCorr(0)
,hPtDstarD0_AFCorr(0)
,hPtDp_AFCorr(0)
,hPtDs_AFCorr(0)
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
,HistWeightD0(0)
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
,fDcorr(0)
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
,fCalculateBDMesonpTWeights(0)
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
    
    if(fIsMC && fCalculateBDMesonpTWeights){
    
    	
    
        fBHadpT = new TH1F("fBHadpT","B hadron pT;p_{T} (GeV/c);counts",250,0,50);
        fBHadpT->Sumw2();
        fOutputList->Add(fBHadpT);
        
        fBMesonpT = new TH1F("fBMesonpT","B meson pT;p_{T} (GeV/c);counts",250,0,50);
        fBMesonpT->Sumw2();
        fOutputList->Add(fBMesonpT);
        
        fBMesonpT_Corr = new TH1F("fBMesonpT_Corr","B meson pT;p_{T} (GeV/c);counts",250,0,50);
        fBMesonpT_Corr->Sumw2();
        fOutputList->Add(fBMesonpT_Corr);
        
     /*   fBMesonpTG = new TH1F("fBMesonpTG","B meson pT;p_{T} (GeV/c);counts",250,0,50);
        fBMesonpTG->Sumw2();
        fOutputList->Add(fBMesonpTG);
        
        fBMesonpTGG = new TH1F("fBMesonpTGG","B meson pT;p_{T} (GeV/c);counts",250,0,50);
        fBMesonpTGG->Sumw2();
        fOutputList->Add(fBMesonpTGG);
        */
        fBDHadpT = new TH1F("fBDHadpT","D (<- B) hadron pT;p_{T} (GeV/c);counts",250,0,50);
        fBDHadpT->Sumw2();
        fOutputList->Add(fBDHadpT);
        
        fDHadpT = new TH1F("fDHadpT","Prompt D hadron pT;p_{T} (GeV/c);counts",250,0,50);
        fDHadpT->Sumw2();
        fOutputList->Add(fDHadpT);
        
        fDMesonpT = new TH1F("fDMesonpT","Prompt D meson pT;p_{T} (GeV/c);counts",250,0,50);
        fDMesonpT->Sumw2();
        fOutputList->Add(fDMesonpT);

        fD0pT = new TH1F("fD0pT","Prompt D0 meson pT;p_{T} (GeV/c);counts",250,0,50);
        fD0pT->Sumw2();
        fOutputList->Add(fD0pT);
        
        fLambdaCpT = new TH1F("fLambdaCpT","Prompt Lammda_c pT;p_{T} (GeV/c);counts",250,0,50);
        fLambdaCpT->Sumw2();
        fOutputList->Add(fLambdaCpT);
        
        fLambdaCpT_AFCorr = new TH1F("fLambdaCpT_AFCorr","Prompt Lammda_c _AFCorrpT;p_{T} (GeV/c);counts",250,0,50);
        fLambdaCpT_AFCorr->Sumw2();
        fOutputList->Add(fLambdaCpT_AFCorr);

  }
    
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
    hPtD0->Sumw2();
    fOutputList->Add(hPtD0);
    
    hPtLambdaC = new TH1F("hPtLambdaC","; p_{T} [GeV/c]; Count",1000,0, 50);
    hPtLambdaC->Sumw2();
    fOutputList->Add(hPtLambdaC);
    
    hPtDp = new TH1F("hPtDp","; p_{T} [GeV/c]; Count",1000,0, 50);
    hPtDp->Sumw2();
    fOutputList->Add(hPtDp);
    
    hPtDstar = new TH1F("hPtDstar","; p_{T} [GeV/c]; Count",1000,0, 50);
    hPtDstar->Sumw2();
    fOutputList->Add(hPtDstar);
    
    hPtDstarD0 = new TH1F("hPtDstarD0","; p_{T} [GeV/c]; Count",1000,0, 50);
    hPtDstarD0->Sumw2();
    fOutputList->Add(hPtDstarD0);
    
    hPtDstarDp = new TH1F("hPtDstarDp","; p_{T} [GeV/c]; Count",1000,0, 50);
    hPtDstarDp->Sumw2();
    fOutputList->Add(hPtDstarDp);
    
    hPtDs = new TH1F("hPtDs","; p_{T} [GeV/c]; Count",1000,0, 50);
    hPtDs->Sumw2();
    fOutputList->Add(hPtDs);
    
    hPtD0_AFCorr = new TH1F("hPtD0_AFCorr","; p_{T} [GeV/c]; Count",1000,0, 50);
    hPtD0_AFCorr->Sumw2();
    fOutputList->Add(hPtD0_AFCorr);
    
    hPtDp_AFCorr = new TH1F("hPtDp_AFCorr","; p_{T} [GeV/c]; Count",1000,0, 50);
    hPtDp_AFCorr->Sumw2();
    fOutputList->Add(hPtDp_AFCorr);
    
    hPtDstarD0_AFCorr = new TH1F("hPtDstarD0_AFCorr","; p_{T} [GeV/c]; Count",1000,0, 50);
    hPtDstarD0_AFCorr->Sumw2();
    fOutputList->Add(hPtDstarD0_AFCorr);
       
    hPtDs_AFCorr = new TH1F("hPtDs_AFCorr","; p_{T} [GeV/c]; Count",1000,0, 50);
    hPtDs_AFCorr->Sumw2();
    fOutputList->Add(hPtDs_AFCorr);
    
        
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
    
    if(fIsMC){

  fDCAxy_pt_charmbef = new TH2F("fDCAxy_pt_charmbef",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  //fDCAxy_pt_charmbef->Sumw2();
  fOutputList->Add(fDCAxy_pt_charmbef);
  
  
  fDCAxy_pt_DstarDplusbef = new TH2F("fDCAxy_pt_DstarDplusbef",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  //fDCAxy_pt_DstarDplusbef->Sumw2();
  fOutputList->Add(fDCAxy_pt_DstarDplusbef);
  
  fDCAxy_pt_DstarDplusAft = new TH2F("fDCAxy_pt_DstarDplusAft",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_DstarDplusAft->Sumw2();
  fOutputList->Add(fDCAxy_pt_DstarDplusAft);
  
  fDCAxy_pt_DstarDzerobef = new TH2F("fDCAxy_pt_DstarDzerobef",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  //fDCAxy_pt_DstarDzerobef->Sumw2();
  fOutputList->Add(fDCAxy_pt_DstarDzerobef);
  
  fDCAxy_pt_DstarDzeroAft = new TH2F("fDCAxy_pt_DstarDzeroAft",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_DstarDzeroAft->Sumw2();
  fOutputList->Add(fDCAxy_pt_DstarDzeroAft);
  
  fDCAxy_pt_DstarDsbef = new TH2F("fDCAxy_pt_DstarDsbef",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  //fDCAxy_pt_DstarDsbef->Sumw2();
  fOutputList->Add(fDCAxy_pt_DstarDsbef);
  
  fDCAxy_pt_DstarDsAft = new TH2F("fDCAxy_pt_DstarDsAft",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_DstarDsAft->Sumw2();
  fOutputList->Add(fDCAxy_pt_DstarDsAft);
  
  fDCAxy_pt_Dplusbef = new TH2F("fDCAxy_pt_Dplusbef",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  //fDCAxy_pt_Dplusbef->Sumw2();
  fOutputList->Add(fDCAxy_pt_Dplusbef);
  
  fDCAxy_pt_Dzerobef = new TH2F("fDCAxy_pt_Dzerobef",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  //fDCAxy_pt_Dzerobef->Sumw2();
  fOutputList->Add(fDCAxy_pt_Dzerobef);
  
  fDCAxy_pt_Dsbef = new TH2F("fDCAxy_pt_Dsbef",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  //fDCAxy_pt_Dsbef->Sumw2();
  fOutputList->Add(fDCAxy_pt_Dsbef);
  
  fDCAxy_pt_DplusAft = new TH2F("fDCAxy_pt_DplusAft",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_DplusAft->Sumw2();
  fOutputList->Add(fDCAxy_pt_DplusAft);
  
  fDCAxy_pt_DzeroAft = new TH2F("fDCAxy_pt_DzeroAft",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_DzeroAft->Sumw2();
  fOutputList->Add(fDCAxy_pt_DzeroAft);
  
  fDCAxy_pt_DsAft = new TH2F("fDCAxy_pt_DsAft",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_DsAft->Sumw2();
  fOutputList->Add(fDCAxy_pt_DsAft);
  
  fDCAxy_pt_charmmesonsbef = new TH2F("fDCAxy_pt_charmmesonsbef",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  //fDCAxy_pt_charmmesonsbef->Sumw2();
  fOutputList->Add(fDCAxy_pt_charmmesonsbef);
  
  fDCAxy_pt_charmmesonsAft = new TH2F("fDCAxy_pt_charmmesonsAft",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_charmmesonsAft->Sumw2();
  fOutputList->Add(fDCAxy_pt_charmmesonsAft);
  
  fDCAxy_pt_Lc = new TH2F("fDCAxy_pt_Lc",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  //fDCAxy_pt_Lc->Sumw2();
  fOutputList->Add(fDCAxy_pt_Lc);
  
  fDCAxy_pt_charmbaryons = new TH2F("fDCAxy_pt_charmbaryons",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  //fDCAxy_pt_charmbaryons->Sumw2();
  fOutputList->Add(fDCAxy_pt_charmbaryons);
  
  fDCAxy_pt_beautybaryons = new TH2F("fDCAxy_pt_beautybaryons",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  //fDCAxy_pt_beautybaryons->Sumw2();
  fOutputList->Add(fDCAxy_pt_beautybaryons);
  
  
  /* fDCAxy_pt_Lc_corr = new TH2F("fDCAxy_pt_Lc_corr",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_Lc_corr->Sumw2();
  fOutputList->Add(fDCAxy_pt_Lc_corr);
  
  fDCAxy_pt_charmbaryons_corr = new TH2F("fDCAxy_pt_charmbaryons_corr",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_charmbaryons_corr->Sumw2();
  fOutputList->Add(fDCAxy_pt_charmbaryons_corr);
  */
  fDCAxy_pt_beautybaryons_corr = new TH2F("fDCAxy_pt_beautybaryons_corr",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_beautybaryons_corr->Sumw2();
  fOutputList->Add(fDCAxy_pt_beautybaryons_corr);

  fDCAxy_pt_charmaft = new TH2F("fDCAxy_pt_charmaft",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_charmaft->Sumw2();
  fOutputList->Add(fDCAxy_pt_charmaft);

  fDCAxy_pt_beautybef = new TH2F("fDCAxy_pt_beautybef",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  //fDCAxy_pt_beautybef->Sumw2();
  fOutputList->Add(fDCAxy_pt_beautybef);

  fDCAxy_pt_beautyaft = new TH2F("fDCAxy_pt_beautyaft",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_beautyaft->Sumw2();
  fOutputList->Add(fDCAxy_pt_beautyaft); 

  fDCAxy_pt_MesonB_beautybef = new TH2F("fDCAxy_pt_MesonB_beautybef",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  //fDCAxy_pt_MesonB_beautybef->Sumw2();
  fOutputList->Add(fDCAxy_pt_MesonB_beautybef);

  fDCAxy_pt_MesonB_beautyaft = new TH2F("fDCAxy_pt_MesonB_beautyaft",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_MesonB_beautyaft->Sumw2();
  fOutputList->Add(fDCAxy_pt_MesonB_beautyaft);

  fDCAxy_pt_MesonBD_beautybef = new TH2F("fDCAxy_pt_MesonBD_beautybef",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  //fDCAxy_pt_MesonBD_beautybef->Sumw2();
  fOutputList->Add(fDCAxy_pt_MesonBD_beautybef);

  fDCAxy_pt_MesonBD_beautyaft = new TH2F("fDCAxy_pt_MesonBD_beautyaft",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  fDCAxy_pt_MesonBD_beautyaft->Sumw2();
  fOutputList->Add(fDCAxy_pt_MesonBD_beautyaft);

  fDCAxy_pt_BaryonB_beautybef = new TH2F("fDCAxy_pt_BaryonB_beautybef",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
 // fDCAxy_pt_BaryonB_beautybef->Sumw2();
  fOutputList->Add(fDCAxy_pt_BaryonB_beautybef);

  fDCAxy_pt_BaryonBD_beautybef = new TH2F("fDCAxy_pt_BaryonBD_beautybef",";p_{t} (GeV/c);DCAxy hadrons",32,ptbinning,2000,-0.5,0.5);
  //fDCAxy_pt_BaryonBD_beautybef->Sumw2();
  fOutputList->Add(fDCAxy_pt_BaryonBD_beautybef);
  }
       
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
    

    fDCAxy_pt_ele_onlyDCA = new TH2F("fDCAxy_pt_ele_onlyDCA",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_onlyDCA);


    fDCAxy_pt_ele_onlyDCA_phi1 = new TH2F("fDCAxy_pt_ele_onlyDCA_phi1",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_onlyDCA_phi1);

    fDCAxy_pt_ele_phi1_ChB = new TH2F("fDCAxy_pt_ele_phi1_ChB",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_phi1_ChB);

    fDCAxy_pt_ele_phi1_B = new TH2F("fDCAxy_pt_ele_phi1_B",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_phi1_B);

    fDCAxy_pt_ele_onlyDCA_phi2 = new TH2F("fDCAxy_pt_ele_onlyDCA_phi2",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_onlyDCA_phi2);

    fDCAxy_pt_ele_phi2_ChB = new TH2F("fDCAxy_pt_ele_phi2_ChB",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_phi2_ChB);

    fDCAxy_pt_ele_phi2_B = new TH2F("fDCAxy_pt_ele_phi2_B",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_phi2_B);

    fDCAxy_pt_ele_onlyDCA_phi3 = new TH2F("fDCAxy_pt_ele_onlyDCA_phi3",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_onlyDCA_phi3);

    fDCAxy_pt_ele_phi3_ChB = new TH2F("fDCAxy_pt_ele_phi3_ChB",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_phi3_ChB);

    fDCAxy_pt_ele_phi3_B = new TH2F("fDCAxy_pt_ele_phi3_B",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_phi3_B);

    fDCAxy_pt_ele_onlyDCA_phi4 = new TH2F("fDCAxy_pt_ele_onlyDCA_phi4",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_onlyDCA_phi4);

    fDCAxy_pt_ele_phi4_ChB = new TH2F("fDCAxy_pt_ele_phi4_ChB",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_phi4_ChB);

    fDCAxy_pt_ele_phi4_B = new TH2F("fDCAxy_pt_ele_phi4_B",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_phi4_B);
    
       
    fDCAxy_pt_had_ResCorr_phi1 = new TH2F("fDCAxy_pt_had_ResCorr_phi1",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_ResCorr_phi1);
    
    fDCAxy_pt_had_ResCorr_phi2 = new TH2F("fDCAxy_pt_had_ResCorr_phi2",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_ResCorr_phi2);
    
    fDCAxy_pt_had_ResCorr_phi3 = new TH2F("fDCAxy_pt_had_ResCorr_phi3",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_ResCorr_phi3);
    
    fDCAxy_pt_had_ResCorr_phi4 = new TH2F("fDCAxy_pt_had_ResCorr_phi4",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_had_ResCorr_phi4);
    

    fDCAxy_pt_ele_ResCorr_phi1 = new TH2F("fDCAxy_pt_ele_ResCorr_phi1",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_ResCorr_phi1);

    fDCAxy_pt_ele_ResCorr_phi2 = new TH2F("fDCAxy_pt_ele_ResCorr_phi2",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_ResCorr_phi2);

    fDCAxy_pt_ele_ResCorr_phi3 = new TH2F("fDCAxy_pt_ele_ResCorr_phi3",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_ResCorr_phi3);

    fDCAxy_pt_ele_ResCorr_phi4 = new TH2F("fDCAxy_pt_ele_ResCorr_phi4",";p_{t} (GeV/c);DCAxy hadrons",300,0,30,2000,-0.5,0.5);
    fOutputList->Add(fDCAxy_pt_ele_ResCorr_phi4);

    
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
    Int_t nBina3[nDima3] = {32,nBinsdcaxy,nBinsITSchi2,nBinsITSsha,nBinspdg2,nBinsdcaxy};
    fD0Data = new THnSparseF("fD0Data","fD0Data",nDima3,nBina3);
    fD0Data->SetBinEdges(0,ptbinning); ///pt spectra -> same binning as other histograms
    fD0Data->SetBinEdges(1,binLimdcaxy); ///dca distribution
    fD0Data->SetBinEdges(2,binLimITSchi2); ///ITS chi2 
    fD0Data->SetBinEdges(3,binLimITSsha); ///fraction ITS shared clusters 
    fD0Data->SetBinEdges(4,binLimpdg2); /// electrons and pions
    fD0Data->SetBinEdges(5,binLimdcaxy); /// electrons and pions
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
    Double_t fTOFnSigma_kaon = -999;
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
        
        if(fCalculateBDMesonpTWeights) GetMCTemplateWeight(fMCarray); // For B meson weight 
        
        for(Int_t iMC = 0; iMC < fMCarray->GetEntries(); iMC++)
        {
			fMCparticle = (AliAODMCParticle*) fMCarray->At(iMC);
			///Pseudo-rapidity cut
			if((fMCparticle->Eta() < fEtaMin) || (fMCparticle->Eta() > fEtaMax)) continue;     				
			
            ///Get proportion between generated D0 and LambdaC:
            Int_t TrackPDG = TMath::Abs(fMCparticle->GetPdgCode());
            //if(TrackPDG == 421) hPtD0->Fill(fMCparticle->Pt());///D0
            if(TrackPDG == 4122) hPtLambdaC->Fill(fMCparticle->Pt()); ///LambdaC
           // if(TrackPDG == 411) hPtDp->Fill(fMCparticle->Pt()); ///D+
           // if(TrackPDG == 413) hPtDstar->Fill(fMCparticle->Pt()); ///D*
           // if(TrackPDG == 431) hPtDs->Fill(fMCparticle->Pt()); ///Ds
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
	fTOFnSigma_kaon = fPidResponse->NumberOfSigmasTOF(track, AliPID::kKaon);        

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
         Double_t d0z0[2] = {999, 999}, cov[3] = {999, 999, 999};
         //AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
         // if(!(track->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 3., d0z0, cov))) continue; 
         //cout<<d0z0[0]<<"    "<<d0z0[1]<<endl;
         fExtraCuts->GetHFEImpactParameters(track, d0z0, cov); // recalculation of vertex is done here, earlier was not done in the task and was the reason for "shoulder" shape in the DCA templates. Also, this is not giving effect for PbPb but only pp. ====> Sudhir 19 January, 2019 ///Solved
                                                                                                                                                          Double_t DCAxy = d0z0[0];                                                                                                                       Double_t DCAz = d0z0[1];
                                                                                                                                                         //cout<<"After   "<<DCAxy<<"         "<<DCAz<<endl;
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
			//	if(TMath::Abs(fMCparticle->GetPdgCode()) == 321) fKaonsPt->Fill(fPt);
				}
			
            /*if(fITSnSigma >= -2 && fITSnSigma <= 2){
                fTPCnsigma_p_after_tof_its->Fill(fP,fTPCnSigma);
                fTPCnsigma_pt_after_tof_its->Fill(fPt,fTPCnSigma);
            }*/
        }
        if(fIsMC){
	if(fTOFnSigma_kaon >= -3 && fTOFnSigma_kaon <= 3){
	   if(fTPCnSigma_kaon >= -3 && fTPCnSigma_kaon <= 3){
                    if(TMath::Abs(fMCparticle->GetPdgCode()) == 321) fKaonsPt->Fill(fPt);
	  }
	}
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
	       float DCAMCMean_phi1 = 999;
            	float DCAMCMean_phi2 = 999;
            	float DCAMCMean_phi3 = 999;
            	float DCAMCMean_phi4 = 999;
            	
            	float DCAMCRes_phi1 = 999;
            	float DCAMCRes_phi2 = 999;
            	float DCAMCRes_phi3 = 999;
            	float DCAMCRes_phi4 = 999;
            	
            	float correction_phi1 = 999;
            	float correction_phi2 = 999;
            	float correction_phi3 = 999;
            	float correction_phi4 = 999;
	       Double_t phi_d0 = (track->Phi()*180.0/TMath::Pi());
        if(fTPCnSigma >= -5 && fTPCnSigma <= -3){
			fDCAxy_pt_had_onlyDCA->Fill(fPt,DCAxy);
			fDCAxy_pt_had->Fill(fPt,DCAxy*track->Charge()*signB);
			fDCAz_pt_had->Fill(fPt,TMath::Sqrt(cov[0]));	
			
			if(phi_d0 > 315.0 || phi_d0 < 45.0){
			fDCAxy_pt_had_onlyDCA_phi1->Fill(fPt,DCAxy);
			fDCAxy_pt_had_phi1_ChB->Fill(fPt,DCAxy*track->Charge()*signB);
			fDCAxy_pt_had_phi1_B->Fill(fPt,TMath::Sqrt(cov[0]));
			
			if(fIsMC){
			DCAMCRes_phi1 = GetDCAResolMC_phi1(fPt); ///resolution of the MC
			}else{
			DCAMCRes_phi1 = 0.0;
			}
			if(!fIsMC){
			DCAMCMean_phi1 = GetDCAMeanData_phi1(fPt); ///mean of the MC
			}
			else{
			DCAMCMean_phi1 = GetDCAMeanMC_phi1(fPt);
			}
			correction_phi1 = gRandom->Gaus(DCAMCMean_phi1,DCAMCRes_phi1); 
			float DCAResCorr_phi1 =  DCAxy + correction_phi1;
			
			//fResGausCorr_phi1->Fill(correction_phi1);
			fDCAxy_pt_had_ResCorr_phi1->Fill(fPt,DCAResCorr_phi1);
			//cout<<"Phi1 value:===   "<<phi_d0<<endl;
			}
			
			if(phi_d0 > 45.0 && phi_d0 < 135.0){
			fDCAxy_pt_had_onlyDCA_phi2->Fill(fPt,DCAxy);
			fDCAxy_pt_had_phi2_ChB->Fill(fPt,DCAxy*track->Charge()*signB);
			fDCAxy_pt_had_phi2_B->Fill(fPt,TMath::Sqrt(cov[0]));
			if(fIsMC){
			DCAMCRes_phi2 = GetDCAResolMC_phi2(fPt); ///resolution of the MC
			}else{
			DCAMCRes_phi2 = 0.0;
			}
			if(!fIsMC){
			DCAMCMean_phi2 = GetDCAMeanData_phi2(fPt); ///mean of the MC
			}
			else{
			DCAMCMean_phi2 = GetDCAMeanMC_phi2(fPt);
			}
			correction_phi2 = gRandom->Gaus(DCAMCMean_phi2,DCAMCRes_phi2); 
			float DCAResCorr_phi2 =  DCAxy + correction_phi2;
			
			//fResGausCorr_phi2->Fill(correction_phi2);
			fDCAxy_pt_had_ResCorr_phi2->Fill(fPt,DCAResCorr_phi2);
			//cout<<"Phi2 value:===   "<<phi_d0<<endl;
			}
			
			if(phi_d0 > 135.0 && phi_d0 < 225.0){
			fDCAxy_pt_had_onlyDCA_phi3->Fill(fPt,DCAxy);
			fDCAxy_pt_had_phi3_ChB->Fill(fPt,DCAxy*track->Charge()*signB);
			fDCAxy_pt_had_phi3_B->Fill(fPt,TMath::Sqrt(cov[0]));
			if(fIsMC){
			DCAMCRes_phi3 = GetDCAResolMC_phi3(fPt); ///resolution of the MC
			}else{
			DCAMCRes_phi3 = 0.0;
			}
			if(!fIsMC){
			DCAMCMean_phi3 = GetDCAMeanData_phi3(fPt); ///mean of the MC
			}
			else{
			DCAMCMean_phi3 = GetDCAMeanMC_phi3(fPt);
			}
			correction_phi3 = gRandom->Gaus(DCAMCMean_phi3,DCAMCRes_phi3); 
			float DCAResCorr_phi3 =  DCAxy + correction_phi3;
			
			//fResGausCorr_phi3->Fill(correction_phi3);
			fDCAxy_pt_had_ResCorr_phi3->Fill(fPt,DCAResCorr_phi3);
			//cout<<"Phi3 value:===   "<<phi_d0<<endl;
			}
			
			if(phi_d0 > 225.0 && phi_d0 < 315.0){
			fDCAxy_pt_had_onlyDCA_phi4->Fill(fPt,DCAxy);
			fDCAxy_pt_had_phi4_ChB->Fill(fPt,DCAxy*track->Charge()*signB);
			fDCAxy_pt_had_phi4_B->Fill(fPt,TMath::Sqrt(cov[0]));
			if(fIsMC){
			DCAMCRes_phi4 = GetDCAResolMC_phi4(fPt); ///resolution of the MC
			}else{
			DCAMCRes_phi4 = 0.0;
			}
			if(!fIsMC){
			DCAMCMean_phi4 = GetDCAMeanData_phi4(fPt); ///mean of the MC
			}
			else{
			DCAMCMean_phi4 = GetDCAMeanMC_phi4(fPt);
			}
			correction_phi4 = gRandom->Gaus(DCAMCMean_phi4,DCAMCRes_phi4); 
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
		float DataDCAResCorr = 999;
            	float DCAMeanCorr_phi1 = 999;
            	float DCAMeanCorr_phi2 = 999;
            	float DCAMeanCorr_phi3 = 999;
            	float DCAMeanCorr_phi4 = 999;
            	
		if(!fIsMC){
		if(phi_d0 > 315.0 || phi_d0 < 45.0){
		DCAMCMean_phi1 = GetDCAMeanData_phi1(fPt); ///mean of the MC
		float Meancorrection_phi1 = gRandom->Gaus(DCAMCMean_phi1,0.0); 
		DataDCAResCorr =  (DCAxy + Meancorrection_phi1)*track->Charge()*signB;
		}
			
		if(phi_d0 > 45.0 && phi_d0 < 135.0){
		DCAMCMean_phi2 = GetDCAMeanData_phi2(fPt); ///mean of the MC
		float Meancorrection_phi2 = gRandom->Gaus(DCAMCMean_phi2,0.0); 
		DataDCAResCorr =  (DCAxy + Meancorrection_phi2)*track->Charge()*signB;
		}
			
		if(phi_d0 > 135.0 && phi_d0 < 225.0){
		DCAMCMean_phi3 = GetDCAMeanData_phi3(fPt); ///mean of the MC
		float Meancorrection_phi3 = gRandom->Gaus(DCAMCMean_phi3,0.0); 
		DataDCAResCorr =  (DCAxy + Meancorrection_phi3)*track->Charge()*signB;
		}
			
		if(phi_d0 > 225.0 && phi_d0 < 315.0){
		DCAMCMean_phi4 = GetDCAMeanData_phi4(fPt); ///mean of the MC
		float Meancorrection_phi4 = gRandom->Gaus(DCAMCMean_phi4,0.0);  
		DataDCAResCorr =  (DCAxy + Meancorrection_phi4)*track->Charge()*signB;
		}
 		}
		
         if(!fIsMC){
			 qadcaData[0] = fPt;
         
			 qadcaData[5] = DCAxy*track->Charge()*signB;
			
			 qadcaData[4] = -1.;
			
			 qadcaData[1] = DataDCAResCorr;
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
        //if(fTPCnSigma >= ftpcPIDmincut && fTPCnSigma <= ftpcPIDmaxcut) cout<<"PIDssssss  "<<fTPCnSigma<<endl;
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
        fDCAxy_pt_ele_onlyDCA->Fill(fPt,DCAxy);
        fDCAxy_pt_ele->Fill(fPt,DCAxy*track->Charge()*signB);
        fDCAz_pt_ele->Fill(fPt,TMath::Sqrt(cov[0]));
        
	// Manual correction to the mean and resolution of the elcttron dca dated 
		/*	if(phi_d0 > 315.0 || phi_d0 < 45.0){
			fDCAxy_pt_ele_onlyDCA_phi1->Fill(fPt,DCAxy);
			fDCAxy_pt_ele_phi1_ChB->Fill(fPt,DCAxy*track->Charge()*signB);
			fDCAxy_pt_ele_phi1_B->Fill(fPt,TMath::Sqrt(cov[0]));
			float DCAMCRes_phi1 = GetDCAResolMC_phi1(fPt); ///resolution of the MC
			float DCAMCMean_phi1 = GetDCAMeanMC_phi1(fPt); ///mean of the MC
			float correction_phi1 = gRandom->Gaus(DCAMCMean_phi1,DCAMCRes_phi1); 
			float DCAResCorr_phi1 =  DCAxy + correction_phi1;
          		fDCAxy_pt_ele_ResCorr_phi1->Fill(fPt,DCAResCorr_phi1);
			}	
	
			if(phi_d0 > 45.0 && phi_d0 < 135.0){
                        fDCAxy_pt_ele_onlyDCA_phi2->Fill(fPt,DCAxy);
                        fDCAxy_pt_ele_phi2_ChB->Fill(fPt,DCAxy*track->Charge()*signB);
                        fDCAxy_pt_ele_phi2_B->Fill(fPt,TMath::Sqrt(cov[0]));
                        float DCAMCRes_phi2 = GetDCAResolMC_phi2(fPt); ///resolution of the MC
                        float DCAMCMean_phi2 = GetDCAMeanMC_phi2(fPt); ///mean of the MC
                        float correction_phi2 = gRandom->Gaus(DCAMCMean_phi2,DCAMCRes_phi2);
                        float DCAResCorr_phi2 =  DCAxy + correction_phi2;
                        fDCAxy_pt_ele_ResCorr_phi2->Fill(fPt,DCAResCorr_phi2);
                        }

			if(phi_d0 > 135.0 && phi_d0 < 225.0){
                        fDCAxy_pt_ele_onlyDCA_phi3->Fill(fPt,DCAxy);
                        fDCAxy_pt_ele_phi3_ChB->Fill(fPt,DCAxy*track->Charge()*signB);
                        fDCAxy_pt_ele_phi3_B->Fill(fPt,TMath::Sqrt(cov[0]));
                        float DCAMCRes_phi3 = GetDCAResolMC_phi3(fPt); ///resolution of the MC
                        float DCAMCMean_phi3 = GetDCAMeanMC_phi3(fPt); ///mean of the MC
                        float correction_phi3 = gRandom->Gaus(DCAMCMean_phi3,DCAMCRes_phi3);
                        float DCAResCorr_phi3 =  DCAxy + correction_phi3;
                        fDCAxy_pt_ele_ResCorr_phi3->Fill(fPt,DCAResCorr_phi3);
                        }

			if(phi_d0 > 225.0 && phi_d0 < 315.0){
                        fDCAxy_pt_ele_onlyDCA_phi4->Fill(fPt,DCAxy);
                        fDCAxy_pt_ele_phi4_ChB->Fill(fPt,DCAxy*track->Charge()*signB);
                        fDCAxy_pt_ele_phi4_B->Fill(fPt,TMath::Sqrt(cov[0]));
                        float DCAMCRes_phi4 = GetDCAResolMC_phi4(fPt); ///resolution of the MC
                        float DCAMCMean_phi4 = GetDCAMeanMC_phi4(fPt); ///mean of the MC
                        float correction_phi4 = gRandom->Gaus(DCAMCMean_phi4,DCAMCRes_phi4);
                        float DCAResCorr_phi4 =  DCAxy + correction_phi4;
                        fDCAxy_pt_ele_ResCorr_phi4->Fill(fPt,DCAResCorr_phi4);
                        }*/



        
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
            	//float DCAMeanCorr = 999;
            	
            	float DCAResCorr_phi1 = 999;
            	float DCAResCorr_phi2 = 999;
            	float DCAResCorr_phi3 = 999;
            	float DCAResCorr_phi4 = 999;
            	
            	
            	
            	
            	if(fIsMC){
           	if(phi_d0 > 315.0 || phi_d0 < 45.0){
           	DCAMCMean_phi1 = GetDCAMeanMC_phi1(fPt);
		DCAMCRes_phi1 = GetDCAResolMC_phi1(fPt); ///resolution of the MC
		float Rescorrection_phi1 = gRandom->Gaus(DCAMCMean_phi1,DCAMCRes_phi1); 
		DCAResCorr =  (DCAxy + Rescorrection_phi1)*track->Charge()*signB;
		}
			
		if(phi_d0 > 45.0 && phi_d0 < 135.0){
		DCAMCMean_phi2 = GetDCAMeanMC_phi2(fPt);
		DCAMCRes_phi2 = GetDCAResolMC_phi2(fPt); ///resolution of the MC
		float Rescorrection_phi2 = gRandom->Gaus(DCAMCMean_phi2,DCAMCRes_phi2); 
		DCAResCorr =  (DCAxy + Rescorrection_phi2)*track->Charge()*signB;
		}
			
		if(phi_d0 > 135.0 && phi_d0 < 225.0){
		DCAMCMean_phi3 = GetDCAMeanMC_phi3(fPt);
		DCAMCRes_phi3 = GetDCAResolMC_phi3(fPt); ///resolution of the MC
		float Rescorrection_phi3 = gRandom->Gaus(DCAMCMean_phi3,DCAMCRes_phi3); 
		DCAResCorr =  (DCAxy + Rescorrection_phi3)*track->Charge()*signB;
		}
			
		if(phi_d0 > 225.0 && phi_d0 < 315.0){
		DCAMCMean_phi4 = GetDCAMeanMC_phi4(fPt);
		DCAMCRes_phi4 = GetDCAResolMC_phi4(fPt); ///resolution of the MC
		float Rescorrection_phi4 = gRandom->Gaus(DCAMCMean_phi4,DCAMCRes_phi4); 
		DCAResCorr =  (DCAxy + Rescorrection_phi4)*track->Charge()*signB;
		}
		}
		
		
		
 
	    Bool_t IsHFEMC = IsHFelectronsMC(track);
            if(IsHFEMC){
            
               Double_t fWeight_Dp = 1.0;
               Double_t fWeight_D0 = 1.0;
               Double_t fWeight_Ds = 1.0;
               Double_t fWeight_Dstar = 1.0;
               Double_t fWeight_Lc = 1.0;
               Double_t fDWeight = 1.0;
               Double_t fBWeight = 1.0;
               
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
						 if(fMCparticleMother->Pt() > 36) continue;
						 
						
	    fDWeight = fDcorr->Eval(fMCparticleMother->Pt());
            if(TMath::Abs(pdg_mother) == 411 || TMath::Abs(pdg_mother) == 421 || TMath::Abs(pdg_mother) == 431){
                if(TMath::Abs(pdg_mother) == 411){
                  if(TMath::Abs(pdg_gmother) == 413){
                  if(fMCparticleMother->Pt() < 30) fWeight_Dp = HistWeightDpD0->Eval(fMCparticleMother->Pt());
                  hEleDp_WCorr_WFCorr->Fill(fPt, DCAResCorr, fDWeight*fWeight_Dp);
		  hPtDstarDp->Fill(fMCparticleGMother->Pt(), fDWeight);///DstarDp
		  hPtDp->Fill(fMCparticleMother->Pt(), fDWeight);///Dp
		  hPtDp_AFCorr->Fill(fMCparticleMother->Pt(), fDWeight*fWeight_Dp);///Dp
                  fDCAxy_pt_DstarDplusAft->Fill(fPt,DCAResCorr, fDWeight);
                  hCharmMotherPt_corr->Fill(fMCparticleMother->Pt(),fDWeight);
                  hCharmMotherPt_corr2->Fill(fMCparticleMother->Pt(),fDWeight);
                  hElecPt_vsCharmMotherPt_corr->Fill(fMCparticleMother->Pt(),fPt,fDWeight);
                  fDCAxy_pt_charmaft->Fill(fPt,DCAResCorr,fDWeight);
                  }else{
                  if(fMCparticleMother->Pt() < 30) fWeight_Dp = HistWeightDpD0->Eval(fMCparticleMother->Pt());
		  hEleDp_WCorr_WFCorr->Fill(fPt, DCAResCorr, fDWeight*fWeight_Dp);
		  hPtDp->Fill(fMCparticleMother->Pt(), fDWeight);///Dp
		  hPtDp_AFCorr->Fill(fMCparticleMother->Pt(), fDWeight*fWeight_Dp);///Dp
                  fDCAxy_pt_DplusAft->Fill(fPt,DCAResCorr, fDWeight);
                  hCharmMotherPt_corr->Fill(fMCparticleMother->Pt(),fDWeight);
                  hCharmMotherPt_corr2->Fill(fMCparticleMother->Pt(),fDWeight);
                  hElecPt_vsCharmMotherPt_corr->Fill(fMCparticleMother->Pt(),fPt,fDWeight);
                  fDCAxy_pt_charmaft->Fill(fPt,DCAResCorr,fDWeight);
                  }
                  }
                if(TMath::Abs(pdg_mother) == 421){
                  if(TMath::Abs(pdg_gmother) == 413){
                  if(fMCparticleMother->Pt() < 30){   
                  fWeight_D0 = HistWeightD0->Eval(fMCparticleMother->Pt());            
                  hPtD0_AFCorr->Fill(fMCparticleMother->Pt(), fDWeight*fWeight_D0);///D0 
                  hPtD0->Fill(fMCparticleMother->Pt(), fDWeight);///D0
                   }

		  if(fMCparticleGMother->Pt() < 30) fWeight_Dstar = HistWeightDstarD0->Eval(fMCparticleGMother->Pt());
		  hEleDstar_WCorr_WFCorr->Fill(fPt, DCAResCorr, fDWeight*fWeight_Dstar);
		  hPtDstarD0->Fill(fMCparticleGMother->Pt(), fDWeight);///DstarD0
		  hPtDstarD0_AFCorr->Fill(fMCparticleGMother->Pt(), fDWeight*fWeight_Dstar);///Dp
		  fDCAxy_pt_DstarDzeroAft->Fill(fPt,DCAResCorr, fDWeight);
                  hCharmMotherPt_corr->Fill(fMCparticleMother->Pt(),fDWeight);
                  hCharmMotherPt_corr2->Fill(fMCparticleMother->Pt(),fDWeight);
                  hElecPt_vsCharmMotherPt_corr->Fill(fMCparticleMother->Pt(),fPt,fDWeight);
                  fDCAxy_pt_charmaft->Fill(fPt,DCAResCorr,fDWeight);
 
                  }else{
                  if(fMCparticleMother->Pt() < 30) fWeight_D0 = HistWeightD0->Eval(fMCparticleMother->Pt());
                  hPtD0->Fill(fMCparticleMother->Pt(), fDWeight);///D0
                  hPtD0_AFCorr->Fill(fMCparticleMother->Pt(), fDWeight*fWeight_D0);///D0
                  hEleD0_WCorr->Fill(fPt, DCAResCorr, fDWeight*fWeight_D0);
                  fDCAxy_pt_DzeroAft->Fill(fPt,DCAResCorr, fDWeight);// create weight for D0 to regain the actual fraction
                  hCharmMotherPt_corr->Fill(fMCparticleMother->Pt(),fDWeight);
                  hCharmMotherPt_corr2->Fill(fMCparticleMother->Pt(),fDWeight);
                  hElecPt_vsCharmMotherPt_corr->Fill(fMCparticleMother->Pt(),fPt,fDWeight);
                  fDCAxy_pt_charmaft->Fill(fPt,DCAResCorr,fDWeight);
                  }
                  }
                if(TMath::Abs(pdg_mother) == 431){
                  if(TMath::Abs(pdg_gmother) == 413){
                  if(fMCparticleMother->Pt() < 30) fWeight_Ds = HistWeightDsD0->Eval(fMCparticleMother->Pt());
		  hEleDs_WCorr_WFCorr->Fill(fPt, DCAResCorr, fDWeight*fWeight_Ds);
		  hPtDs->Fill(fMCparticleMother->Pt(), fDWeight);///Ds
		  hPtDs_AFCorr->Fill(fMCparticleMother->Pt(), fDWeight*fWeight_Ds);///Ds
                  fDCAxy_pt_DstarDsAft->Fill(fPt,DCAResCorr, fDWeight);
                  hCharmMotherPt_corr->Fill(fMCparticleMother->Pt(),fDWeight);
                  hCharmMotherPt_corr2->Fill(fMCparticleMother->Pt(),fDWeight);
                  hElecPt_vsCharmMotherPt_corr->Fill(fMCparticleMother->Pt(),fPt,fDWeight);
                  fDCAxy_pt_charmaft->Fill(fPt,DCAResCorr,fDWeight);
                  }
                  else{
                  if(fMCparticleMother->Pt() < 30) fWeight_Ds = HistWeightDsD0->Eval(fMCparticleMother->Pt());
		  hEleDs_WCorr_WFCorr->Fill(fPt, DCAResCorr, fDWeight*fWeight_Ds);
		  hPtDs->Fill(fMCparticleMother->Pt(), fDWeight);///Ds
		  hPtDs_AFCorr->Fill(fMCparticleMother->Pt(), fDWeight*fWeight_Ds);///Ds
                  fDCAxy_pt_DsAft->Fill(fPt,DCAResCorr, fDWeight);
                  hCharmMotherPt_corr->Fill(fMCparticleMother->Pt(),fDWeight);
                  hCharmMotherPt_corr2->Fill(fMCparticleMother->Pt(),fDWeight);
                  hElecPt_vsCharmMotherPt_corr->Fill(fMCparticleMother->Pt(),fPt,fDWeight);
                  fDCAxy_pt_charmaft->Fill(fPt,DCAResCorr,fDWeight);
                  }
                  }
                  }
                  else{
                  fDCAxy_pt_charmmesonsAft->Fill(fPt,DCAResCorr, fDWeight);
                  hCharmMotherPt_corr->Fill(fMCparticleMother->Pt(),fDWeight);
                  hCharmMotherPt_corr2->Fill(fMCparticleMother->Pt(),fDWeight);
                  hElecPt_vsCharmMotherPt_corr->Fill(fMCparticleMother->Pt(),fPt,fDWeight);
                  fDCAxy_pt_charmaft->Fill(fPt,DCAResCorr,fDWeight);
                  }
	      
                  hCharmMotherPt_vsElecPt_corr->Fill(fPt,fMCparticleMother->Pt());
                  	
                  }
					
							
					
		 if(TMath::Abs(pdg_mother)>4000 && TMath::Abs(pdg_mother)<5000){///charmed baryon
			 qadca[1]=1.5;
				 if(TMath::Abs(pdg_mother) == 4122){ 
					 qadca[1]=27.5; ///to check DCA prompt Lc+
					 if(fMCparticleMother->Pt() < 30) fWeight_Lc = HistWeightLcD0->Eval(fMCparticleMother->Pt());
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
                    fDCAxy_pt_MesonB_beautybef->Fill(fPt,DCAResCorr);
                    }
                    if(fIsFromMesonBD){
                    
                    fDCAxy_pt_MesonBD_beautybef->Fill(fPt,DCAResCorr);
                    }
                    // check by Sudhir on 27 Oct 2018
                    if(fIsFromBarionB){
                    fDCAxy_pt_BaryonB_beautybef->Fill(fPt,DCAResCorr);
                    }
                    if(fIsFromBarionBD){
                    
                    fDCAxy_pt_BaryonBD_beautybef->Fill(fPt,DCAResCorr);
                    }
                    // check by Sudhir on 27 Oct 2018
                    if(fIsFromMesonB || fIsFromMesonBD){///beauty meson 
						qadca[1]=2.5;
						fDCAxy_pt_beautybef->Fill(fPt,DCAResCorr);
						            ///Correcting the pT spectrum
            
            if(fIsFromMesonB && (fMCparticleMother->Pt() < 50)){
              fBWeight = fBcorr->Eval(fMCparticleMother->Pt());///always evaluating the function in the pT of the B meson
              //cout<<"pdg_mother = "<<pdg_mother<<endl;
              hBeautyMotherPt->Fill(fMCparticleMother->Pt(),fPt);
             // hBeautyMotherPtbef->Fill(fMCparticleMother->Pt());
              hBeautyMotherPtaft->Fill(fMCparticleMother->Pt(),fBWeight);
              hBeautyMotherPt2Daft->Fill(fMCparticleMother->Pt(),fPt,fBWeight);
              fDCAxy_pt_MesonB_beautyaft->Fill(fPt,DCAResCorr, fBWeight);
              fDCAxy_pt_beautyaft->Fill(fPt,DCAResCorr, fBWeight);	
            }
            if(fIsFromMesonBD){
              if(fMCparticleMother->GetMother() > 0){
                fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
                float pdg_gmother = fMCparticleGMother->GetPdgCode();

                if(TMath::Abs(pdg_gmother)>500 && TMath::Abs(pdg_gmother)<600 && (fMCparticleGMother->Pt()<50)){
                  fBWeight = fBcorr->Eval(fMCparticleGMother->Pt());
                  hBeautyMotherPt->Fill(fMCparticleGMother->Pt(),fPt);
                 // hBeautyMotherPtbef->Fill(fMCparticleGMother->Pt());
                  hBeautyMotherPtaft->Fill(fMCparticleGMother->Pt(),fBWeight);
                  hBeautyMotherPt2Daft->Fill(fMCparticleGMother->Pt(),fPt,fBWeight);
                  fDCAxy_pt_MesonBD_beautyaft->Fill(fPt,DCAResCorr, fBWeight);
                  fDCAxy_pt_beautyaft->Fill(fPt,DCAResCorr, fBWeight);
                 
                }
                else{
                if(fMCparticleGMother->GetMother() > 0){
                  fMCparticleGGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleGMother->GetMother());
                  float pdg_ggmother = fMCparticleGGMother->GetPdgCode();
                  if(TMath::Abs(pdg_ggmother)>500 && TMath::Abs(pdg_ggmother)<600 && (fMCparticleGGMother->Pt()<50)){ 
                  fBWeight = fBcorr->Eval(fMCparticleGGMother->Pt());
                  hBeautyMotherPt->Fill(fMCparticleGGMother->Pt(),fPt);
                 // hBeautyMotherPtbef->Fill(fMCparticleGGMother->Pt());
                  hBeautyMotherPtaft->Fill(fMCparticleGGMother->Pt(),fBWeight);
                  hBeautyMotherPt2Daft->Fill(fMCparticleGGMother->Pt(),fPt,fBWeight);
                  fDCAxy_pt_MesonBD_beautyaft->Fill(fPt,DCAResCorr, fBWeight);
                  fDCAxy_pt_beautyaft->Fill(fPt,DCAResCorr, fBWeight);
                  }else{ 
                  if(fMCparticleGGMother->GetMother() > 0){
                    fMCparticleGGGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleGGMother->GetMother());
                    float pdg_gggmother = fMCparticleGGGMother->GetPdgCode();
                    if (TMath::Abs(pdg_gggmother)>500 && TMath::Abs(pdg_gggmother)<600 && (fMCparticleGGGMother->Pt()<50)){             
                    fBWeight = fBcorr->Eval(fMCparticleGGGMother->Pt());
                    hBeautyMotherPt->Fill(fMCparticleGGGMother->Pt(),fPt);
            	  //  hBeautyMotherPtbef->Fill(fMCparticleGGGMother->Pt());
                    hBeautyMotherPtaft->Fill(fMCparticleGGGMother->Pt(),fBWeight);
                    hBeautyMotherPt2Daft->Fill(fMCparticleGGGMother->Pt(),fPt,fBWeight);
                    fDCAxy_pt_MesonBD_beautyaft->Fill(fPt,DCAResCorr, fBWeight);
                    fDCAxy_pt_beautyaft->Fill(fPt,DCAResCorr, fBWeight);
                    }
                  }
                  }
                }
                }
              }
            }	
						
					}
					
	if(fIsFromBarionB || fIsFromBarionBD){///beauty baryon
          if(fIsFromBarionB && (fMCparticleMother->Pt()<50)){
          fBWeight = fBcorr->Eval(fMCparticleMother->Pt());
          fDCAxy_pt_beautybaryons_corr->Fill(fPt,DCAResCorr, fBWeight);
          }
          if(fIsFromBarionBD){
          if(fMCparticleMother->GetMother() > 0){
          fMCparticleGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleMother->GetMother());
          float pdg_gmother = fMCparticleGMother->GetPdgCode();
          if((TMath::Abs(pdg_gmother)>5000 && TMath::Abs(pdg_gmother)<6000) && (fMCparticleGMother->Pt()<50)){
          fBWeight = fBcorr->Eval(fMCparticleGMother->Pt());
          fDCAxy_pt_beautybaryons_corr->Fill(fPt,DCAResCorr, fBWeight);
          }
          else{
          if(fMCparticleGMother->GetMother() > 0){
          fMCparticleGGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleGMother->GetMother());
          float pdg_ggmother = fMCparticleGGMother->GetPdgCode();
          if((TMath::Abs(pdg_ggmother)>5000 && TMath::Abs(pdg_ggmother)<6000) && (fMCparticleGGMother->Pt()<50)){
          fBWeight = fBcorr->Eval(fMCparticleGGMother->Pt());
          fDCAxy_pt_beautybaryons_corr->Fill(fPt,DCAResCorr, fBWeight);
          }
          else{
          if(fMCparticleGGMother->GetMother() > 0){
          fMCparticleGGGMother = (AliAODMCParticle*) fMCarray->At(fMCparticleGGMother->GetMother());
          float pdg_gggmother = 999;
          if(fMCparticleGGGMother) float pdg_gggmother = fMCparticleGGGMother->GetPdgCode();
          if((TMath::Abs(pdg_gggmother)>5000 && TMath::Abs(pdg_gggmother)<6000) && (fMCparticleGGGMother->Pt()<50)){
          fBWeight = fBcorr->Eval(fMCparticleGGGMother->Pt());
          fDCAxy_pt_beautybaryons_corr->Fill(fPt,DCAResCorr, fBWeight);
          }
          }
          }
          }
          }
          }
          }
          qadca[1]=3.5;
          fDCAxy_pt_beautybaryons->Fill(fPt,DCAResCorr);
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
                //Int_t motherlabel = pi0eta->GetMother();
                Int_t gmoth = fMCparticleMother->GetMother();
		Int_t motherpdg = 999;
		if(gmoth > 0){
		AliAODMCParticle *mother = (AliAODMCParticle*)fMCarray->At(gmoth);
                motherpdg = TMath::Abs(mother->GetPdgCode());
                }
        	//if(fIsMC) GetGammaAndDalitzElectronTemplates(fMCarray, track, fPt, DCANew);
        	//if(fIsMC) GetGammaAndDalitzElectronTemplates2(fMCarray, track, fPt, DCANew, PhotonicType2);
                
                
      
                
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
					if((TMath::Abs(pdg_mother) == 111 || TMath::Abs(pdg_mother) == 221)){
					if(TMath::Abs(motherpdg) == 321 && fType == 4){ // electrons from Dalitz 
      					fDCAxy_pt_Dalitz2->Fill(fPt, DCAResCorr, fKPicorr->Eval(fMCparticleMother->Pt()));
      					}
					if((TMath::Abs(motherpdg) != 321 && fType == 4) || fType == 0 || fType == 1 || fType == 2 || fType == 3){
      					fDCAxy_pt_Dalitz2->Fill(fPt, DCAResCorr);
      					}
          				}
					if(TMath::Abs(pdg_mother) == 22){
      					if(TMath::Abs(motherpdg) == 321 && fType == 4){ // electrons from photon conversions
      					fDCAxy_pt_Conversions2->Fill(fPt, DCAResCorr, fKPicorr->Eval(fMCparticleMother->Pt()));
      					}
					if((TMath::Abs(motherpdg) != 321 && fType == 4) || fType == 0 || fType == 1 || fType == 2 || fType == 3){
      					fDCAxy_pt_Conversions2->Fill(fPt, DCAResCorr);
      					}
					}
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

void AliAnalysisHFEppTPCTOFBeauty::GetMCTemplateWeight(TClonesArray *fMCArray)
{
    //Get histograms for D,B and Lamdac weight calculation
    
    AliAODMCParticle *MCPart;
    AliAODMCParticle *MCPartMom;
    AliAODMCParticle *MCPartGMom;
    AliAODMCParticle *MCPartGGMom;
    
    Double_t PartPt = -999, fBWeight = -999;
    Int_t iMCmom = -999, iMCgmom = -999;
    Int_t MomPDG = -999, GMomPDG=-999;
    
    for(Int_t imcArrayL=0; imcArrayL< fMCArray->GetEntries(); imcArrayL++){
        MCPart = (AliAODMCParticle*)fMCArray->At(imcArrayL);
        Int_t PDGcode = TMath::Abs(MCPart->GetPdgCode());
        
        iMCmom = -999, iMCgmom = -999;
        MomPDG = -999, GMomPDG=-999;
        PartPt = -999;
        
        Bool_t IsMCHF = kFALSE, IsMCD = kFALSE, IsMCB = kFALSE, IsMCBD = kFALSE;
        
        if(TMath::Abs(MCPart->Eta()) > 0.8) continue;
        
        PartPt = MCPart->Pt();
        
        if((PDGcode>400 && PDGcode<600) || (PDGcode>4000 && PDGcode<6000)){
            IsMCHF = kTRUE;
            
            if((PDGcode>500 && PDGcode<600) || (PDGcode>5000 && PDGcode<6000)){
                IsMCB = kTRUE;
                fBHadpT->Fill(PartPt);
                
                if(PDGcode>500 && PDGcode<600){ 
                fBMesonpT->Fill(PartPt);
                fBWeight = fBcorr->Eval(PartPt);
                fBMesonpT_Corr->Fill(PartPt,fBWeight);                
                }
            }
            else{
                iMCmom = MCPart->GetMother();
                if(iMCmom > 0){
                    MCPartMom = (AliAODMCParticle*)fMCArray->At(iMCmom);
                    MomPDG = TMath::Abs(MCPartMom->GetPdgCode());
                    
                    if((MomPDG>500 && MomPDG<600) || (MomPDG>5000 && MomPDG<6000)){
                        IsMCB = kTRUE;
                        IsMCBD = kTRUE;
                        fBDHadpT->Fill(MCPartMom->Pt());
                     //   if(MomPDG>500 && MomPDG<600) fBMesonpTG->Fill(MCPartMom->Pt());
                    }
                    else{
                        iMCgmom = MCPartMom->GetMother();
                        if(iMCgmom > 0){
                            MCPartGMom = (AliAODMCParticle*)fMCArray->At(iMCgmom);
                            GMomPDG = TMath::Abs(MCPartGMom->GetPdgCode());
                            
                            if((GMomPDG>500 && GMomPDG<600) || (GMomPDG>5000 && GMomPDG<6000)){
                                IsMCB = kTRUE;
                                IsMCBD = kTRUE;
                                fBDHadpT->Fill(MCPartGMom->Pt());
                              //  if(GMomPDG>500 && GMomPDG<600) fBMesonpTGG->Fill(MCPartGMom->Pt());
                            }
                        }
                    }
                }
            }
            
            if(!IsMCB) {
                if((PDGcode>400 && PDGcode<500) || (PDGcode>4000 && PDGcode<5000)) fDHadpT->Fill(PartPt);
                if(PDGcode > 400 && PDGcode < 500) fDMesonpT->Fill(PartPt);
                if(PDGcode == 421) fD0pT->Fill(PartPt);
                if(PDGcode == 4122) fLambdaCpT->Fill(PartPt);
                if(PDGcode == 4122){
                if(PartPt > 30) continue;
		fLambdaCpT_AFCorr->Fill(PartPt, HistWeightLcD0->Eval(PartPt));
                }
            }
        }
    }
}


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


Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAMeanMC_phi1(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float MeanMCGaus = 0;

if (x >= 0.90 && x < 1.10) MeanMCGaus = 0.000287; 
if (x >= 1.10 && x < 1.30) MeanMCGaus = 0.000247; 
if (x >= 1.30 && x < 1.50) MeanMCGaus = 0.000264; 
if (x >= 1.50 && x < 1.70) MeanMCGaus = 0.000273; 
if (x >= 1.70 && x < 1.90) MeanMCGaus = 0.000286; 
if (x >= 1.90 && x < 2.10) MeanMCGaus = 0.000329; 
if (x >= 2.10 && x < 2.30) MeanMCGaus = 0.000344; 
if (x >= 2.30 && x < 2.50) MeanMCGaus = 0.000350; 
if (x >= 2.50 && x < 2.70) MeanMCGaus = 0.000374; 
if (x >= 2.70 && x < 2.90) MeanMCGaus = 0.000383; 
if (x >= 2.90 && x < 3.10) MeanMCGaus = 0.000392; 
if (x >= 3.10 && x < 3.30) MeanMCGaus = 0.000365; 
if (x >= 3.30 && x < 3.50) MeanMCGaus = 0.000405; 
if (x >= 3.50 && x < 3.70) MeanMCGaus = 0.000384; 
if (x >= 3.70 && x < 3.90) MeanMCGaus = 0.000386; 
if (x >= 3.90 && x < 4.10) MeanMCGaus = 0.000391; 
if (x >= 4.10 && x < 4.30) MeanMCGaus = 0.000385; 
if (x >= 4.30 && x < 4.50) MeanMCGaus = 0.000455; 
if (x >= 4.50 && x < 4.70) MeanMCGaus = 0.000401; 
if (x >= 4.70 && x < 5.00) MeanMCGaus = 0.000431; 
if (x >= 5.00 && x < 5.50) MeanMCGaus = 0.000345; 
if (x >= 5.50 && x < 6.00) MeanMCGaus = 0.000381; 
if (x >= 6.00 && x < 6.50) MeanMCGaus = 0.000322; 
if (x >= 6.50 && x < 7.00) MeanMCGaus = 0.000414; 
if (x >= 7.00 && x < 8.00) MeanMCGaus = 0.000333; 




return MeanMCGaus;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAMeanMC_phi2(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float MeanMCGaus = 0;

if (x >= 0.90 && x < 1.10) MeanMCGaus = -0.000481; 
if (x >= 1.10 && x < 1.30) MeanMCGaus = -0.000431; 
if (x >= 1.30 && x < 1.50) MeanMCGaus = -0.000380; 
if (x >= 1.50 && x < 1.70) MeanMCGaus = -0.000332; 
if (x >= 1.70 && x < 1.90) MeanMCGaus = -0.000292; 
if (x >= 1.90 && x < 2.10) MeanMCGaus = -0.000272; 
if (x >= 2.10 && x < 2.30) MeanMCGaus = -0.000251; 
if (x >= 2.30 && x < 2.50) MeanMCGaus = -0.000223; 
if (x >= 2.50 && x < 2.70) MeanMCGaus = -0.000216; 
if (x >= 2.70 && x < 2.90) MeanMCGaus = -0.000166; 
if (x >= 2.90 && x < 3.10) MeanMCGaus = -0.000175; 
if (x >= 3.10 && x < 3.30) MeanMCGaus = -0.000141; 
if (x >= 3.30 && x < 3.50) MeanMCGaus = -0.000147; 
if (x >= 3.50 && x < 3.70) MeanMCGaus = -0.000111; 
if (x >= 3.70 && x < 3.90) MeanMCGaus = -0.000148; 
if (x >= 3.90 && x < 4.10) MeanMCGaus = -0.000124; 
if (x >= 4.10 && x < 4.30) MeanMCGaus = -0.000110; 
if (x >= 4.30 && x < 4.50) MeanMCGaus = -0.000094; 
if (x >= 4.50 && x < 4.70) MeanMCGaus = -0.000103; 
if (x >= 4.70 && x < 5.00) MeanMCGaus = -0.000039; 
if (x >= 5.00 && x < 5.50) MeanMCGaus = -0.000058; 
if (x >= 5.50 && x < 6.00) MeanMCGaus = -0.000028; 
if (x >= 6.00 && x < 6.50) MeanMCGaus = -0.000105; 
if (x >= 6.50 && x < 7.00) MeanMCGaus = -0.000033; 
if (x >= 7.00 && x < 8.00) MeanMCGaus = -0.000032; 



return MeanMCGaus;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAMeanMC_phi3(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float MeanMCGaus = 0;

if (x >= 0.90 && x < 1.10) MeanMCGaus = 0.000319; 
if (x >= 1.10 && x < 1.30) MeanMCGaus = 0.000265; 
if (x >= 1.30 && x < 1.50) MeanMCGaus = 0.000269; 
if (x >= 1.50 && x < 1.70) MeanMCGaus = 0.000268; 
if (x >= 1.70 && x < 1.90) MeanMCGaus = 0.000273; 
if (x >= 1.90 && x < 2.10) MeanMCGaus = 0.000288; 
if (x >= 2.10 && x < 2.30) MeanMCGaus = 0.000279; 
if (x >= 2.30 && x < 2.50) MeanMCGaus = 0.000261; 
if (x >= 2.50 && x < 2.70) MeanMCGaus = 0.000270; 
if (x >= 2.70 && x < 2.90) MeanMCGaus = 0.000264; 
if (x >= 2.90 && x < 3.10) MeanMCGaus = 0.000240; 
if (x >= 3.10 && x < 3.30) MeanMCGaus = 0.000281; 
if (x >= 3.30 && x < 3.50) MeanMCGaus = 0.000248; 
if (x >= 3.50 && x < 3.70) MeanMCGaus = 0.000231; 
if (x >= 3.70 && x < 3.90) MeanMCGaus = 0.000265; 
if (x >= 3.90 && x < 4.10) MeanMCGaus = 0.000244; 
if (x >= 4.10 && x < 4.30) MeanMCGaus = 0.000284; 
if (x >= 4.30 && x < 4.50) MeanMCGaus = 0.000215; 
if (x >= 4.50 && x < 4.70) MeanMCGaus = 0.000193; 
if (x >= 4.70 && x < 5.00) MeanMCGaus = 0.000207; 
if (x >= 5.00 && x < 5.50) MeanMCGaus = 0.000260; 
if (x >= 5.50 && x < 6.00) MeanMCGaus = 0.000213; 
if (x >= 6.00 && x < 6.50) MeanMCGaus = 0.000275; 
if (x >= 6.50 && x < 7.00) MeanMCGaus = 0.000161; 
if (x >= 7.00 && x < 8.00) MeanMCGaus = 0.000179; 

return MeanMCGaus;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAMeanMC_phi4(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float MeanMCGaus = 0;

if (x >= 0.90 && x < 1.10) MeanMCGaus = -0.000715; 
if (x >= 1.10 && x < 1.30) MeanMCGaus = -0.000714; 
if (x >= 1.30 && x < 1.50) MeanMCGaus = -0.000609; 
if (x >= 1.50 && x < 1.70) MeanMCGaus = -0.000480; 
if (x >= 1.70 && x < 1.90) MeanMCGaus = -0.000378; 
if (x >= 1.90 && x < 2.10) MeanMCGaus = -0.000264; 
if (x >= 2.10 && x < 2.30) MeanMCGaus = -0.000195; 
if (x >= 2.30 && x < 2.50) MeanMCGaus = -0.000122; 
if (x >= 2.50 && x < 2.70) MeanMCGaus = -0.000088; 
if (x >= 2.70 && x < 2.90) MeanMCGaus = -0.000031; 
if (x >= 2.90 && x < 3.10) MeanMCGaus = -0.000017; 
if (x >= 3.10 && x < 3.30) MeanMCGaus = 0.000041; 
if (x >= 3.30 && x < 3.50) MeanMCGaus = 0.000067; 
if (x >= 3.50 && x < 3.70) MeanMCGaus = 0.000112; 
if (x >= 3.70 && x < 3.90) MeanMCGaus = 0.000139; 
if (x >= 3.90 && x < 4.10) MeanMCGaus = 0.000177; 
if (x >= 4.10 && x < 4.30) MeanMCGaus = 0.000169; 
if (x >= 4.30 && x < 4.50) MeanMCGaus = 0.000262; 
if (x >= 4.50 && x < 4.70) MeanMCGaus = 0.000213; 
if (x >= 4.70 && x < 5.00) MeanMCGaus = 0.000265; 
if (x >= 5.00 && x < 5.50) MeanMCGaus = 0.000276; 
if (x >= 5.50 && x < 6.00) MeanMCGaus = 0.000286; 
if (x >= 6.00 && x < 6.50) MeanMCGaus = 0.000278; 
if (x >= 6.50 && x < 7.00) MeanMCGaus = 0.000298; 
if (x >= 7.00 && x < 8.00) MeanMCGaus = 0.000324; 




return MeanMCGaus;   
}


Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAMeanData_phi1(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float MeanDataGaus = 0;

if (x >= 0.90 && x < 1.10) MeanDataGaus = 0.000240; 
if (x >= 1.10 && x < 1.30) MeanDataGaus = 0.000222; 
if (x >= 1.30 && x < 1.50) MeanDataGaus = 0.000208; 
if (x >= 1.50 && x < 1.70) MeanDataGaus = 0.000191; 
if (x >= 1.70 && x < 1.90) MeanDataGaus = 0.000184; 
if (x >= 1.90 && x < 2.10) MeanDataGaus = 0.000162; 
if (x >= 2.10 && x < 2.30) MeanDataGaus = 0.000144; 
if (x >= 2.30 && x < 2.50) MeanDataGaus = 0.000130; 
if (x >= 2.50 && x < 2.70) MeanDataGaus = 0.000117; 
if (x >= 2.70 && x < 2.90) MeanDataGaus = 0.000093; 
if (x >= 2.90 && x < 3.10) MeanDataGaus = 0.000077; 
if (x >= 3.10 && x < 3.30) MeanDataGaus = 0.000068; 
if (x >= 3.30 && x < 3.50) MeanDataGaus = 0.000089; 
if (x >= 3.50 && x < 3.70) MeanDataGaus = 0.000087; 
if (x >= 3.70 && x < 3.90) MeanDataGaus = 0.000067; 
if (x >= 3.90 && x < 4.10) MeanDataGaus = 0.000057; 
if (x >= 4.10 && x < 4.30) MeanDataGaus = 0.000036; 
if (x >= 4.30 && x < 4.50) MeanDataGaus = 0.000049; 
if (x >= 4.50 && x < 4.70) MeanDataGaus = 0.000036; 
if (x >= 4.70 && x < 5.00) MeanDataGaus = 0.000046; 
if (x >= 5.00 && x < 5.50) MeanDataGaus = 0.000052; 
if (x >= 5.50 && x < 6.00) MeanDataGaus = 0.000016; 
if (x >= 6.00 && x < 6.50) MeanDataGaus = 0.000032; 
if (x >= 6.50 && x < 7.00) MeanDataGaus = 0.000025; 
if (x >= 7.00 && x < 8.00) MeanDataGaus = 0.000026; 




return MeanDataGaus;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAMeanData_phi2(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float MeanDataGaus = 0;

if (x >= 0.90 && x < 1.10) MeanDataGaus = 0.001212; 
if (x >= 1.10 && x < 1.30) MeanDataGaus = 0.001207; 
if (x >= 1.30 && x < 1.50) MeanDataGaus = 0.001090; 
if (x >= 1.50 && x < 1.70) MeanDataGaus = 0.000980; 
if (x >= 1.70 && x < 1.90) MeanDataGaus = 0.000864; 
if (x >= 1.90 && x < 2.10) MeanDataGaus = 0.000815; 
if (x >= 2.10 && x < 2.30) MeanDataGaus = 0.000778; 
if (x >= 2.30 && x < 2.50) MeanDataGaus = 0.000747; 
if (x >= 2.50 && x < 2.70) MeanDataGaus = 0.000719; 
if (x >= 2.70 && x < 2.90) MeanDataGaus = 0.000650; 
if (x >= 2.90 && x < 3.10) MeanDataGaus = 0.000628; 
if (x >= 3.10 && x < 3.30) MeanDataGaus = 0.000607; 
if (x >= 3.30 && x < 3.50) MeanDataGaus = 0.000550; 
if (x >= 3.50 && x < 3.70) MeanDataGaus = 0.000512; 
if (x >= 3.70 && x < 3.90) MeanDataGaus = 0.000521; 
if (x >= 3.90 && x < 4.10) MeanDataGaus = 0.000521; 
if (x >= 4.10 && x < 4.30) MeanDataGaus = 0.000456; 
if (x >= 4.30 && x < 4.50) MeanDataGaus = 0.000466; 
if (x >= 4.50 && x < 4.70) MeanDataGaus = 0.000437; 
if (x >= 4.70 && x < 5.00) MeanDataGaus = 0.000367; 
if (x >= 5.00 && x < 5.50) MeanDataGaus = 0.000348; 
if (x >= 5.50 && x < 6.00) MeanDataGaus = 0.000328; 
if (x >= 6.00 && x < 6.50) MeanDataGaus = 0.000233; 
if (x >= 6.50 && x < 7.00) MeanDataGaus = 0.000287; 
if (x >= 7.00 && x < 8.00) MeanDataGaus = 0.000331; 



return MeanDataGaus;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAMeanData_phi3(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float MeanDataGaus = 0;

if (x >= 0.90 && x < 1.10) MeanDataGaus = 0.000752; 
if (x >= 1.10 && x < 1.30) MeanDataGaus = 0.000691; 
if (x >= 1.30 && x < 1.50) MeanDataGaus = 0.000595; 
if (x >= 1.50 && x < 1.70) MeanDataGaus = 0.000563; 
if (x >= 1.70 && x < 1.90) MeanDataGaus = 0.000523; 
if (x >= 1.90 && x < 2.10) MeanDataGaus = 0.000527; 
if (x >= 2.10 && x < 2.30) MeanDataGaus = 0.000492; 
if (x >= 2.30 && x < 2.50) MeanDataGaus = 0.000483; 
if (x >= 2.50 && x < 2.70) MeanDataGaus = 0.000455; 
if (x >= 2.70 && x < 2.90) MeanDataGaus = 0.000400; 
if (x >= 2.90 && x < 3.10) MeanDataGaus = 0.000386; 
if (x >= 3.10 && x < 3.30) MeanDataGaus = 0.000362; 
if (x >= 3.30 && x < 3.50) MeanDataGaus = 0.000362; 
if (x >= 3.50 && x < 3.70) MeanDataGaus = 0.000364; 
if (x >= 3.70 && x < 3.90) MeanDataGaus = 0.000353; 
if (x >= 3.90 && x < 4.10) MeanDataGaus = 0.000326; 
if (x >= 4.10 && x < 4.30) MeanDataGaus = 0.000275; 
if (x >= 4.30 && x < 4.50) MeanDataGaus = 0.000241; 
if (x >= 4.50 && x < 4.70) MeanDataGaus = 0.000286; 
if (x >= 4.70 && x < 5.00) MeanDataGaus = 0.000296; 
if (x >= 5.00 && x < 5.50) MeanDataGaus = 0.000259; 
if (x >= 5.50 && x < 6.00) MeanDataGaus = 0.000256; 
if (x >= 6.00 && x < 6.50) MeanDataGaus = 0.000259; 
if (x >= 6.50 && x < 7.00) MeanDataGaus = 0.000189; 
if (x >= 7.00 && x < 8.00) MeanDataGaus = 0.000193; 


return MeanDataGaus;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAMeanData_phi4(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float MeanDataGaus = 0;

if (x >= 0.90 && x < 1.10) MeanDataGaus = 0.001308; 
if (x >= 1.10 && x < 1.30) MeanDataGaus = 0.001265; 
if (x >= 1.30 && x < 1.50) MeanDataGaus = 0.001213; 
if (x >= 1.50 && x < 1.70) MeanDataGaus = 0.001140; 
if (x >= 1.70 && x < 1.90) MeanDataGaus = 0.001086; 
if (x >= 1.90 && x < 2.10) MeanDataGaus = 0.001029; 
if (x >= 2.10 && x < 2.30) MeanDataGaus = 0.000997; 
if (x >= 2.30 && x < 2.50) MeanDataGaus = 0.000952; 
if (x >= 2.50 && x < 2.70) MeanDataGaus = 0.000910; 
if (x >= 2.70 && x < 2.90) MeanDataGaus = 0.000881; 
if (x >= 2.90 && x < 3.10) MeanDataGaus = 0.000829; 
if (x >= 3.10 && x < 3.30) MeanDataGaus = 0.000823; 
if (x >= 3.30 && x < 3.50) MeanDataGaus = 0.000795; 
if (x >= 3.50 && x < 3.70) MeanDataGaus = 0.000776; 
if (x >= 3.70 && x < 3.90) MeanDataGaus = 0.000748; 
if (x >= 3.90 && x < 4.10) MeanDataGaus = 0.000699; 
if (x >= 4.10 && x < 4.30) MeanDataGaus = 0.000682; 
if (x >= 4.30 && x < 4.50) MeanDataGaus = 0.000697; 
if (x >= 4.50 && x < 4.70) MeanDataGaus = 0.000629; 
if (x >= 4.70 && x < 5.00) MeanDataGaus = 0.000628; 
if (x >= 5.00 && x < 5.50) MeanDataGaus = 0.000593; 
if (x >= 5.50 && x < 6.00) MeanDataGaus = 0.000575; 
if (x >= 6.00 && x < 6.50) MeanDataGaus = 0.000549; 
if (x >= 6.50 && x < 7.00) MeanDataGaus = 0.000514; 
if (x >= 7.00 && x < 8.00) MeanDataGaus = 0.000497; 




return MeanDataGaus;   
}


Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAResolMC_phi1(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float sigmaG = 0;

if (x >= 0.90 && x < 1.10) sigmaG = 0.003428; 
if (x >= 1.10 && x < 1.30) sigmaG = 0.003423; 
if (x >= 1.30 && x < 1.50) sigmaG = 0.003152; 
if (x >= 1.50 && x < 1.70) sigmaG = 0.002919; 
if (x >= 1.70 && x < 1.90) sigmaG = 0.002768; 
if (x >= 1.90 && x < 2.10) sigmaG = 0.002610; 
if (x >= 2.10 && x < 2.30) sigmaG = 0.002536; 
if (x >= 2.30 && x < 2.50) sigmaG = 0.002457; 
if (x >= 2.50 && x < 2.70) sigmaG = 0.002381; 
if (x >= 2.70 && x < 2.90) sigmaG = 0.002357; 
if (x >= 2.90 && x < 3.10) sigmaG = 0.002301; 
if (x >= 3.10 && x < 3.30) sigmaG = 0.002264; 
if (x >= 3.30 && x < 3.50) sigmaG = 0.002099; 
if (x >= 3.50 && x < 3.70) sigmaG = 0.002012; 
if (x >= 3.70 && x < 3.90) sigmaG = 0.002033; 
if (x >= 3.90 && x < 4.10) sigmaG = 0.001888; 
if (x >= 4.10 && x < 4.30) sigmaG = 0.002147; 
if (x >= 4.30 && x < 4.50) sigmaG = 0.001799; 
if (x >= 4.50 && x < 4.70) sigmaG = 0.001770; 
if (x >= 4.70 && x < 5.00) sigmaG = 0.001844; 
if (x >= 5.00 && x < 5.50) sigmaG = 0.001677; 
if (x >= 5.50 && x < 6.00) sigmaG = 0.001887; 
if (x >= 6.00 && x < 6.50) sigmaG = 0.000853; 
if (x >= 6.50 && x < 7.00) sigmaG = 0.000596; 
if (x >= 7.00 && x < 8.00) sigmaG = 0.001381; 


return sigmaG;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAResolMC_phi2(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float sigmaG = 0;

if (x >= 0.90 && x < 1.10) sigmaG = 0.004476; 
if (x >= 1.10 && x < 1.30) sigmaG = 0.004611; 
if (x >= 1.30 && x < 1.50) sigmaG = 0.004288; 
if (x >= 1.50 && x < 1.70) sigmaG = 0.004009; 
if (x >= 1.70 && x < 1.90) sigmaG = 0.003741; 
if (x >= 1.90 && x < 2.10) sigmaG = 0.003590; 
if (x >= 2.10 && x < 2.30) sigmaG = 0.003467; 
if (x >= 2.30 && x < 2.50) sigmaG = 0.003400; 
if (x >= 2.50 && x < 2.70) sigmaG = 0.003280; 
if (x >= 2.70 && x < 2.90) sigmaG = 0.003161; 
if (x >= 2.90 && x < 3.10) sigmaG = 0.002931; 
if (x >= 3.10 && x < 3.30) sigmaG = 0.002851; 
if (x >= 3.30 && x < 3.50) sigmaG = 0.002780; 
if (x >= 3.50 && x < 3.70) sigmaG = 0.002647; 
if (x >= 3.70 && x < 3.90) sigmaG = 0.002662; 
if (x >= 3.90 && x < 4.10) sigmaG = 0.002580; 
if (x >= 4.10 && x < 4.30) sigmaG = 0.002494; 
if (x >= 4.30 && x < 4.50) sigmaG = 0.002365; 
if (x >= 4.50 && x < 4.70) sigmaG = 0.002301; 
if (x >= 4.70 && x < 5.00) sigmaG = 0.002402; 
if (x >= 5.00 && x < 5.50) sigmaG = 0.001829; 
if (x >= 5.50 && x < 6.00) sigmaG = 0.002080; 
if (x >= 6.00 && x < 6.50) sigmaG = 0.001396; 
if (x >= 6.50 && x < 7.00) sigmaG = 0.001628; 
if (x >= 7.00 && x < 8.00) sigmaG = 0.001690; 


return sigmaG;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAResolMC_phi3(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float sigmaG = 0;

if (x >= 0.90 && x < 1.10) sigmaG = 0.004735; 
if (x >= 1.10 && x < 1.30) sigmaG = 0.004451; 
if (x >= 1.30 && x < 1.50) sigmaG = 0.003968; 
if (x >= 1.50 && x < 1.70) sigmaG = 0.003580; 
if (x >= 1.70 && x < 1.90) sigmaG = 0.003393; 
if (x >= 1.90 && x < 2.10) sigmaG = 0.003200; 
if (x >= 2.10 && x < 2.30) sigmaG = 0.003074; 
if (x >= 2.30 && x < 2.50) sigmaG = 0.002994; 
if (x >= 2.50 && x < 2.70) sigmaG = 0.002927; 
if (x >= 2.70 && x < 2.90) sigmaG = 0.002727; 
if (x >= 2.90 && x < 3.10) sigmaG = 0.002650; 
if (x >= 3.10 && x < 3.30) sigmaG = 0.002579; 
if (x >= 3.30 && x < 3.50) sigmaG = 0.002454; 
if (x >= 3.50 && x < 3.70) sigmaG = 0.002488; 
if (x >= 3.70 && x < 3.90) sigmaG = 0.002288; 
if (x >= 3.90 && x < 4.10) sigmaG = 0.002302; 
if (x >= 4.10 && x < 4.30) sigmaG = 0.002333; 
if (x >= 4.30 && x < 4.50) sigmaG = 0.002051; 
if (x >= 4.50 && x < 4.70) sigmaG = 0.001959; 
if (x >= 4.70 && x < 5.00) sigmaG = 0.002156; 
if (x >= 5.00 && x < 5.50) sigmaG = 0.001944; 
if (x >= 5.50 && x < 6.00) sigmaG = 0.001678; 
if (x >= 6.00 && x < 6.50) sigmaG = 0.001600; 
if (x >= 6.50 && x < 7.00) sigmaG = 0.001869; 
if (x >= 7.00 && x < 8.00) sigmaG = 0.001882; 


return sigmaG;   
}

Float_t AliAnalysisHFEppTPCTOFBeauty::GetDCAResolMC_phi4(Float_t x){
    
// Return the DCA resolution of the track (in MC) accordingly to its pT

float sigmaG = 0;

if (x >= 0.90 && x < 1.10) sigmaG = 0.004705; 
if (x >= 1.10 && x < 1.30) sigmaG = 0.004542; 
if (x >= 1.30 && x < 1.50) sigmaG = 0.004177; 
if (x >= 1.50 && x < 1.70) sigmaG = 0.003836; 
if (x >= 1.70 && x < 1.90) sigmaG = 0.003587; 
if (x >= 1.90 && x < 2.10) sigmaG = 0.003383; 
if (x >= 2.10 && x < 2.30) sigmaG = 0.003227; 
if (x >= 2.30 && x < 2.50) sigmaG = 0.003058; 
if (x >= 2.50 && x < 2.70) sigmaG = 0.003007; 
if (x >= 2.70 && x < 2.90) sigmaG = 0.002960; 
if (x >= 2.90 && x < 3.10) sigmaG = 0.002778; 
if (x >= 3.10 && x < 3.30) sigmaG = 0.002682; 
if (x >= 3.30 && x < 3.50) sigmaG = 0.002704; 
if (x >= 3.50 && x < 3.70) sigmaG = 0.002549; 
if (x >= 3.70 && x < 3.90) sigmaG = 0.002497; 
if (x >= 3.90 && x < 4.10) sigmaG = 0.002218; 
if (x >= 4.10 && x < 4.30) sigmaG = 0.002384; 
if (x >= 4.30 && x < 4.50) sigmaG = 0.002116; 
if (x >= 4.50 && x < 4.70) sigmaG = 0.002154; 
if (x >= 4.70 && x < 5.00) sigmaG = 0.002214; 
if (x >= 5.00 && x < 5.50) sigmaG = 0.001752; 
if (x >= 5.50 && x < 6.00) sigmaG = 0.001890; 
if (x >= 6.00 && x < 6.50) sigmaG = 0.002123; 
if (x >= 6.50 && x < 7.00) sigmaG = 0.001453; 
if (x >= 7.00 && x < 8.00) sigmaG = 0.001961; 

return sigmaG;   
}

//=======================================================================





