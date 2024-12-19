/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, proviyaded that the above copyright notice appears in all *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is         *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                        //
//                        Analysis for Jet Hadrochemistry                                 //
//                                                                                        //
//    This analysis extracts pT-spectra of charged kaons, protons, and pions              //
//                      for the inclusive event and in jets.                              //
//   It is based on particles identification via the dE/dx signal of the TPC              //
//                   and the time of flight nsigma from the TOF.                          //
//                              This is the AOD version.                                  //
//                                                                                        //
// Author: Sierra Cantway (Weyhmiller) <sierra.lisa.weyhmiller@cern.ch>, Yale University  //
//      Author: Mesut Arslandok <mesut.arslandok@cern.ch>, Yale University                //
//                                                                                        //
////////////////////////////////////////////////////////////////////////////////////////////

//MONTE CARLO TEST VERSION
#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TROOT.h"
#include "TGrid.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#include "THnSparse.h"
#include "TList.h"
#include "TMath.h"
#include "TVectorF.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliPDG.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliLumiTools.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenHepMCEventHeader.h"
#include "AliGenEposEventHeader.h"
#include "AliPID.h"
#include "AliAODVertex.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include <AliAODTrack.h>
#include <AliVTrack.h>
#include "AliVParticle.h"
#include "AliJetContainer.h"
#include "AliTrackContainer.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "AliPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliLog.h"
#include "AliAnalysisTaskJetHadroAOD.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"
#include "AliEventCuts.h"
#include "AliAnalysisUtils.h"
#include "AliFJWrapper.h"
#include "AliEmcalJet.h"
#include "AliEmcalJetTask.h"
#include "AliRhoParameter.h"
#include "AliEmcalContainerUtils.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <bitset>
#include <TObjString.h>
using namespace std;
using std::cout;
using std::setw;

ClassImp(AliAnalysisTaskJetHadroAOD)

#define USE_STREAMER 1


// -----------------------------------------------------------------------
//                            Constructor and Destructor
// -----------------------------------------------------------------------
//________________________________________________________________________
AliAnalysisTaskJetHadroAOD::AliAnalysisTaskJetHadroAOD()
: AliAnalysisTaskEmcalJet("JetHadro"), fEventCuts(0), fPIDResponse(0), fAOD(0), fMCEvent(0), fListHist(0),
fAOD_FilterBits(0),
fPIDCombined(0x0),
fMCStack(0x0),
fVertex_AOD(0x0),
fTreeSRedirector(0x0),
fTreejetsEMCconst(0x0),
fTreeMC(0x0),
fTreeCuts(0x0),
fTreejetsEMCBGconst(0x0),
fTreejetsFJ(0x0),
fTreejetsFJBG(0x0),
fTreejetsFJconst(0x0),
fTreejetsFJBGconst(0x0),
fTreejetEvents(0x0),
fTreejetsEMC(0x0),
fTreejetsEMCBG(0x0),
fRandom(0x0),
fPeriodName(""),
fYear(0),
fPassIndex(0),
fPercentageOfEvents(0),
fSmallOut(kFALSE),
fMCtrue(kFALSE),
fIncludeITS(kTRUE),
fFilljetsFJBGTree(kTRUE),
fFillJetsFJBGConst(kTRUE),
fDoIncTracks(kTRUE),
fDoPerpCone(kTRUE),
fDoRandCone(kTRUE),
fDoFastJet(kTRUE),
fDoEMCJet(kTRUE),
fDoEventTree(kFALSE),
fDoEventHistos(kTRUE),
fFillFastJet(kTRUE),
fFillJetsFJConst(kTRUE),
fFillEMCJet(kTRUE),
fFillIncTracks(kTRUE),
fFill_TPC(kTRUE),
fFillpTPC_pT(kTRUE),
fFillp_pT(kTRUE),
fFillpTPC_p(kTRUE),
fFill_TOF(kTRUE),
fFill_TOF_expecs(kTRUE),
fFill_TPC_expecs(kTRUE),
fFillJetsEMCConst(kTRUE),
fFillJetsEMCBG(kFALSE),
fFillJetsEMCBGConst(kFALSE),
fcent_min(0.0),
fcent_max(100.0),
fDoRapCut(kFALSE),
fEtaCut(0.9),
fYCut(100.0),
fjetMinPtSub(-1000.0),
fjetMaxPtSub(1000.0),
fjetMinArea(-1000.0),
fFillTreeMC(kFALSE),
fUseCouts(kFALSE),
fSetTPCmom(0),
fSetTOFmom(0),
fSetBetamom(0),
fSetEta(0),
fNTPCMom_Bins(2000),
fTPCMom_Bins(0),
fTOFMom_NBins(1000),
fTOFMom_Bins(0),
fNEta_Bins(9),
fEta_Bins(0),
fNdEdxBins(2000),
fdEdxBins(0),
fNBetaBins(1000),
fBetaBins(0),
fNTOFNSigmaBins(1000),
fTOFNSigmaBins(0),
fNSigmasElTOF(0),
fNSigmasMuTOF(0),
fNSigmasPiTOF(0),
fNSigmasKaTOF(0),
fNSigmasPrTOF(0),
fNSigmasDeTOF(0),
fDEdxEl(0),
fDEdxKa(0),
fDEdxPi(0),
fDEdxPr(0),
fDEdxDe(0),
fSigmaEl(0),
fSigmaKa(0),
fSigmaPi(0),
fSigmaPr(0),
fSigmaDe(0),
fNSigmasElTPC(0),
fNSigmasPiTPC(0),
fNSigmasKaTPC(0),
fNSigmasPrTPC(0),
fNSigmasDeTPC(0),
fTPCSignalMC(0),
fPtotMC(0),
fPtMC(0),
fEtaMC(0),
fSignMC(0),
fMCImpactParameter(0),
fNHardScatters(0),
fNProjectileParticipants(0),
fNTargetParticipants(0),
fNNColl(0),
fNNwColl(0),
fNwNColl(0),
fNwNwColl(0),
fPtot(0),
fPVertex(0),
fPt(0),
fY(0),
fCentrality(0),
fCentImpBin(0),
fVz(0),
fEventCountInFile(0),
fTPCSignal(0),
fEta(0),
fNContributors(0),
fPhi(0),
fSign(0),
fJetContainer(0),
fbgJetContainer(0),
fRecoJetContainer(0),
fEmbJetContainer(0),
fGenJetContainer(0),
fJetPt(0),
fJetEta(0),
fJetPhi(0),
fjetRhoVal(0),
fjetRecoRhoVal(0),
fjetEmbRhoVal(0),
fjetGenRhoVal(0),
frhoFJ(0),
fisGoodIncEvent(0),
fFilledUECone_Emb(0),
fFilledUECone_Gen(0),
fFilledUECone_Rec(0),
fhasAcceptedFJjet(0),
fhasAcceptedEMCjet(0),
fhasRealFJjet(0),
fhasRealEMCjet(0),
fNumRealFJJets(0),
fNumRealEMCJets(0),
fIsEmbeddedEvent(kFALSE),
fDoPartLevelMatching(kFALSE),
fDoDetLevelMatching(kFALSE),
fMCParticleArrayName("mcparticles"),
fMCParticleArray(0),
fJetMatchingRadius(0.3),
fTruthMinLabel(0),
fTruthMaxLabel(100000),
fSaveMCInformation(0),
fJetMatchingSharedPtFraction(0.5),
fTrackProbElTPC(0),
fTrackProbPiTPC(0),
fTrackProbKaTPC(0),
fTrackProbPrTPC(0),
fTrackProbDeTPC(0),
fTrackProbElTOF(0),
fTrackProbPiTOF(0),
fTrackProbKaTOF(0),
fTrackProbPrTOF(0),
fTrackProbDeTOF(0),
fHistCentrality(0),
fHistImpParam(0),
fHistVertex(0),
fHist_Is_Good_Inc_Event(0),
fHist_Has_Acc_Real_FJ_Jet(0),
fHist_Has_Acc_Real_EMC_Jet(0),
fHist_Num_Real_FJ_Jets(0),
fHist_Num_Real_EMC_Jets(0),
fHist_RhoFJ(0),
fHist_RhoEMC(0),
fHist_RhoEMCReco(0),
fHist_RhoEMCGen(0),
fHist_RhoEMCEmb(0),
fHist_Filled_UEC_Emb(0),
fHist_Filled_UEC_Gen(0),
fHist_Filled_UEC_Rec(0),
fHistIncTracks_dEdx(0),
fHistIncTracks_moms(0),
fHistIncTracks_moms_p(0),
fHistIncTracks_moms_pTPC_p(0),
fHistIncTracks_kin(0),
fHistIncTracks_beta(0),
fHistIncTracks_t0(0),
fHistIncTracks_TOFpi_nsigma(0),
fHistIncTracks_TOFka_nsigma(0),
fHistIncTracks_TOFpr_nsigma(0),
fHistJetTracks_dEdx(0),
fHistJetTracks_moms(0),
fHistJetTracks_moms_p(0),
fHistJetTracks_moms_pTPC_p(0),
fHistJetTracks_kin(0),
fHistJetTracks_beta(0),
fHistJetTracks_TOFpi_nsigma(0),
fHistJetTracks_TOFka_nsigma(0),
fHistJetTracks_TOFpr_nsigma(0),
fHistIncTracks_mpi(0),
fHistBetaExpec_pi(0),
fHistBetaExpec_ka(0),
fHistBetaExpec_pr(0),
fHistjet_BetaExpec_pi(0),
fHistjet_BetaExpec_ka(0),
fHistjet_BetaExpec_pr(0),
fHist_pi_mismatch(0),
fHist_ka_mismatch(0),
fHist_pr_mismatch(0),
fHist_jet_pi_mismatch(0),
fHist_jet_ka_mismatch(0),
fHist_jet_pr_mismatch(0),
fHist_kaExpec_pihyp(0),
fHist_prExpec_pihyp(0),
fHist_piExpec_kahyp(0),
fHist_prExpec_kahyp(0),
fHist_piExpec_prhyp(0),
fHist_kaExpec_prhyp(0),
fHist_deExpec_prhyp(0),
fHist_jet_elExpec_pihyp(0),
fHist_jet_muExpec_pihyp(0),
fHist_jet_kaExpec_pihyp(0),
fHist_jet_prExpec_pihyp(0),
fHist_jet_piExpec_kahyp(0),
fHist_jet_prExpec_kahyp(0),
fHist_jet_piExpec_prhyp(0),
fHist_jet_kaExpec_prhyp(0),
fHist_jet_deExpec_prhyp(0),
fHistTOFSigmaExpec_pi(0),
fHistTOFSigmaExpec_ka(0),
fHistTOFSigmaExpec_pr(0),
fHistjet_TOFSigmaExpec_pi(0),
fHistjet_TOFSigmaExpec_ka(0),
fHistjet_TOFSigmaExpec_pr(0),
fHistIncTracks_spi(0),
fHistIncTracks_mel(0),
fHistIncTracks_sel(0),
fHistIncTracks_mka(0),
fHistIncTracks_ska(0),
fHistIncTracks_mpr(0),
fHistIncTracks_spr(0),
fHistJet_ptsub_v_area(0),
fHistJet_kin(0),
fHistJet_moms(0),
fHist_pi_DCAxy(0),
fHist_pr_DCAxy(0),
fHistJet_pi_DCAxy(0),
fHistJet_pr_DCAxy(0),
fHistMCTruth_TrackEff_Den_pi(0),
fHistMCTruth_TrackEff_Den_ka(0),
fHistMCTruth_TrackEff_Den_pr(0),
fHistMCReco_TrackEff_Num_pi(0),
fHistMCReco_TrackEff_Num_ka(0),
fHistMCReco_TrackEff_Num_pr(0),
fHistMCReco_MatchEff_Num_pi(0),
fHistMCReco_MatchEff_Num_ka(0),
fHistMCReco_MatchEff_Num_pr(0),
fHistMCReco_pi_prim_DCAxy(0),
fHistMCReco_pi_scdweak_DCAxy(0),
fHistMCReco_pr_prim_DCAxy(0),
fHistMCReco_pr_scdweak_DCAxy(0),
fHistMCReco_pr_scdmat_DCAxy(0),
fHist_Jet_MCTruth_TrackEff_Den_pi(0),
fHist_Jet_MCTruth_TrackEff_Den_ka(0),
fHist_Jet_MCTruth_TrackEff_Den_pr(0),
fHist_Jet_MCReco_TrackEff_Num_pi(0),
fHist_Jet_MCReco_TrackEff_Num_ka(0),
fHist_Jet_MCReco_TrackEff_Num_pr(0),
fHist_Jet_MCReco_MatchEff_Num_pi(0),
fHist_Jet_MCReco_MatchEff_Num_ka(0),
fHist_Jet_MCReco_MatchEff_Num_pr(0),
fHist_JetMatch_MCTruth_TrackEff_Den_pi(0),
fHist_JetMatch_MCTruth_TrackEff_Den_ka(0),
fHist_JetMatch_MCTruth_TrackEff_Den_pr(0),
fHist_JetMatch_MCReco_TrackEff_Num_pi(0),
fHist_JetMatch_MCReco_TrackEff_Num_ka(0),
fHist_JetMatch_MCReco_TrackEff_Num_pr(0),
fHist_JetMatch_MCReco_MatchEff_Num_pi(0),
fHist_JetMatch_MCReco_MatchEff_Num_ka(0),
fHist_JetMatch_MCReco_MatchEff_Num_pr(0),
fHistMCT_Jet_ptsub_v_area(0),
fHistMCR_Jet_ptsub_v_area(0),
fall_reco_jets_w_multiple_matches(0),
fall_reco_jets_w_matches(0),
fHistMC_Jet_shared_pt_frac(0),
fHistMC_Jet_deltaR(0),
fFastJetWrapper(0x0),
fFastJetWrapperBG(0x0),
fFastJetWrapper_Rec(0x0),
fFastJetWrapperBG_Rec(0x0),
fFastJetWrapper_Gen(0x0),
fFastJetWrapperBG_Gen(0x0),
fHist_TrueJetPtFraction(0),
fHist_MatchJetPts(0),
fHist_MatchJetEtas(0),
fHist_MatchJetDeltaRs(0),
fHist_JERS_PbPbunid_tru(0),
fHist_JERS_PbPbunid_det(0),
fHist_JERS_PbPbunid_truhyb(0),
fHist_JERS_PbPbunid_dethyb(0),
fHist_JERS_truPythia(0),
fHist_JERS_truPythia_pi(0),
fHist_JERS_truPythia_ka(0),
fHist_JERS_truPythia_pr(0),
fHist_JERS_detPythia(0),
fHist_JERS_detPythia_pi(0),
fHist_JERS_detPythia_ka(0),
fHist_JERS_detPythia_pr(0),
fHist_JERS_truhybPythia(0),
fHist_JERS_truhybPythia_pi(0),
fHist_JERS_truhybPythia_ka(0),
fHist_JERS_truhybPythia_pr(0),
fHist_JERS_dethybPythia(0),
fHist_JERS_dethybPythia_pi(0),
fHist_JERS_dethybPythia_ka(0),
fHist_JERS_dethybPythia_pr(0),
fHist_PC_spectra_Pythia(0),
fHist_PC_spectra_Pythia_pi(0),
fHist_PC_spectra_Pythia_ka(0),
fHist_PC_spectra_Pythia_pr(0),
fHist_det_spectra_Pythia(0),
fHist_det_spectra_Pythia_pi(0),
fHist_det_spectra_Pythia_ka(0),
fHist_det_spectra_Pythia_pr(0),
fHist_part_spectra_Pythia(0),
fHist_part_spectra_Pythia_pi(0),
fHist_part_spectra_Pythia_ka(0),
fHist_part_spectra_Pythia_pr(0),
fDoEmbedding(kFALSE),
fMultTrueFlag(kFALSE)
{
  // default Constructor
}

//________________________________________________________________________
AliAnalysisTaskJetHadroAOD::AliAnalysisTaskJetHadroAOD(const char *name)
: AliAnalysisTaskEmcalJet(name), fEventCuts(0), fPIDResponse(0), fAOD(0), fMCEvent(0), fListHist(0),
fAOD_FilterBits(0),
fPIDCombined(0x0),
fMCStack(0x0),
fVertex_AOD(0x0),
fTreeSRedirector(0x0),
fTreejetsEMCconst(0x0),
fTreejetsEMCBGconst(0x0),
fTreeMC(0x0),
fTreeCuts(0x0),
fTreejetsFJ(0x0),
fTreejetsFJBG(0x0),
fTreejetsFJconst(0x0),
fTreejetsFJBGconst(0x0),
fTreejetEvents(0x0),
fTreejetsEMC(0x0),
fTreejetsEMCBG(0x0),
fRandom(0x0),
fPeriodName(""),
fYear(0),
fPassIndex(0),
fPercentageOfEvents(0),
fSmallOut(kFALSE),
fMCtrue(kFALSE),
fIncludeITS(kTRUE),
fFilljetsFJBGTree(kTRUE),
fFillJetsFJBGConst(kTRUE),
fDoIncTracks(kTRUE),
fDoPerpCone(kTRUE),
fDoRandCone(kTRUE),
fDoFastJet(kTRUE),
fDoEMCJet(kTRUE),
fDoEventTree(kFALSE),
fDoEventHistos(kTRUE),
fFillFastJet(kTRUE),
fFillJetsFJConst(kTRUE),
fFillEMCJet(kTRUE),
fFillJetsEMCConst(kTRUE),
fFillJetsEMCBG(kFALSE),
fFillJetsEMCBGConst(kFALSE),
fFillIncTracks(kTRUE),
fFill_TPC(kTRUE),
fFillpTPC_pT(kTRUE),
fFillp_pT(kTRUE),
fFillpTPC_p(kTRUE),
fFill_TOF(kTRUE),
fFill_TOF_expecs(kTRUE),
fFill_TPC_expecs(kTRUE),
fcent_min(0.0),
fcent_max(100.0),
fDoRapCut(kFALSE),
fEtaCut(0.9),
fYCut(100.0),
fjetMinPtSub(-1000.0),
fjetMaxPtSub(1000.0),
fjetMinArea(-1000.0),
fFillTreeMC(kFALSE),
fUseCouts(kFALSE),
fSetTPCmom(0),
fSetTOFmom(0),
fSetBetamom(0),
fSetEta(0),
fNTPCMom_Bins(2000),
fTPCMom_Bins(0),
fTOFMom_NBins(1000),
fTOFMom_Bins(0),
fNEta_Bins(9),
fEta_Bins(0),
fNdEdxBins(2000),
fdEdxBins(0),
fNBetaBins(1000),
fBetaBins(0),
fNTOFNSigmaBins(1000),
fTOFNSigmaBins(0),
fNSigmasElTOF(0),
fNSigmasMuTOF(0),
fNSigmasPiTOF(0),
fNSigmasKaTOF(0),
fNSigmasPrTOF(0),
fNSigmasDeTOF(0),
fDEdxEl(0),
fDEdxKa(0),
fDEdxPi(0),
fDEdxPr(0),
fDEdxDe(0),
fSigmaEl(0),
fSigmaKa(0),
fSigmaPi(0),
fSigmaPr(0),
fSigmaDe(0),
fNSigmasElTPC(0),
fNSigmasPiTPC(0),
fNSigmasKaTPC(0),
fNSigmasPrTPC(0),
fNSigmasDeTPC(0),
fTPCSignalMC(0),
fPtotMC(0),
fPtMC(0),
fEtaMC(0),
fSignMC(0),
fMCImpactParameter(0),
fNHardScatters(0),
fNProjectileParticipants(0),
fNTargetParticipants(0),
fNNColl(0),
fNNwColl(0),
fNwNColl(0),
fNwNwColl(0),
fPtot(0),
fPVertex(0),
fPt(0),
fY(0),
fCentrality(0),
fCentImpBin(0),
fVz(0),
fEventCountInFile(0),
fTPCSignal(0),
fEta(0),
fNContributors(0),
fPhi(0),
fSign(0),
fJetContainer(0),
fbgJetContainer(0),
fRecoJetContainer(0),
fEmbJetContainer(0),
fGenJetContainer(0),
fJetPt(0),
fJetEta(0),
fJetPhi(0),
fjetRhoVal(0),
fjetRecoRhoVal(0),
fjetEmbRhoVal(0),
fjetGenRhoVal(0),
frhoFJ(0),
fisGoodIncEvent(0),
fFilledUECone_Emb(0),
fFilledUECone_Gen(0),
fFilledUECone_Rec(0),
fhasAcceptedFJjet(0),
fhasAcceptedEMCjet(0),
fhasRealFJjet(0),
fhasRealEMCjet(0),
fNumRealFJJets(0),
fNumRealEMCJets(0),
fIsEmbeddedEvent(kFALSE),
fDoPartLevelMatching(kFALSE),
fDoDetLevelMatching(kFALSE),
fMCParticleArrayName("mcparticles"),
fMCParticleArray(0),
fJetMatchingRadius(0.3),
fTruthMinLabel(0),
fTruthMaxLabel(100000),
fSaveMCInformation(0),
fJetMatchingSharedPtFraction(0.5),
fTrackProbElTPC(0),
fTrackProbPiTPC(0),
fTrackProbKaTPC(0),
fTrackProbPrTPC(0),
fTrackProbDeTPC(0),
fTrackProbElTOF(0),
fTrackProbPiTOF(0),
fTrackProbKaTOF(0),
fTrackProbPrTOF(0),
fTrackProbDeTOF(0),
fHistCentrality(0),
fHistImpParam(0),
fHistVertex(0),
fHist_Is_Good_Inc_Event(0),
fHist_Has_Acc_Real_FJ_Jet(0),
fHist_Has_Acc_Real_EMC_Jet(0),
fHist_Num_Real_FJ_Jets(0),
fHist_Num_Real_EMC_Jets(0),
fHist_RhoFJ(0),
fHist_RhoEMC(0),
fHist_RhoEMCReco(0),
fHist_RhoEMCGen(0),
fHist_RhoEMCEmb(0),
fHist_Filled_UEC_Emb(0),
fHist_Filled_UEC_Gen(0),
fHist_Filled_UEC_Rec(0),
fHistIncTracks_dEdx(0),
fHistIncTracks_moms(0),
fHistIncTracks_moms_p(0),
fHistIncTracks_moms_pTPC_p(0),
fHistIncTracks_kin(0),
fHistIncTracks_beta(0),
fHistIncTracks_t0(0),
fHistIncTracks_TOFpi_nsigma(0),
fHistIncTracks_TOFka_nsigma(0),
fHistIncTracks_TOFpr_nsigma(0),
fHistJetTracks_dEdx(0),
fHistJetTracks_moms(0),
fHistJetTracks_moms_p(0),
fHistJetTracks_moms_pTPC_p(0),
fHistJetTracks_kin(0),
fHistJetTracks_beta(0),
fHistJetTracks_TOFpi_nsigma(0),
fHistJetTracks_TOFka_nsigma(0),
fHistJetTracks_TOFpr_nsigma(0),
fHistIncTracks_mpi(0),
fHistBetaExpec_pi(0),
fHistBetaExpec_ka(0),
fHistBetaExpec_pr(0),
fHistjet_BetaExpec_pi(0),
fHistjet_BetaExpec_ka(0),
fHistjet_BetaExpec_pr(0),
fHist_pi_mismatch(0),
fHist_ka_mismatch(0),
fHist_pr_mismatch(0),
fHist_jet_pi_mismatch(0),
fHist_jet_ka_mismatch(0),
fHist_jet_pr_mismatch(0),
fHist_kaExpec_pihyp(0),
fHist_prExpec_pihyp(0),
fHist_piExpec_kahyp(0),
fHist_prExpec_kahyp(0),
fHist_piExpec_prhyp(0),
fHist_kaExpec_prhyp(0),
fHist_deExpec_prhyp(0),
fHist_jet_elExpec_pihyp(0),
fHist_jet_muExpec_pihyp(0),
fHist_jet_kaExpec_pihyp(0),
fHist_jet_prExpec_pihyp(0),
fHist_jet_piExpec_kahyp(0),
fHist_jet_prExpec_kahyp(0),
fHist_jet_piExpec_prhyp(0),
fHist_jet_kaExpec_prhyp(0),
fHist_jet_deExpec_prhyp(0),
fHistTOFSigmaExpec_pi(0),
fHistTOFSigmaExpec_ka(0),
fHistTOFSigmaExpec_pr(0),
fHistjet_TOFSigmaExpec_pi(0),
fHistjet_TOFSigmaExpec_ka(0),
fHistjet_TOFSigmaExpec_pr(0),
fHistIncTracks_spi(0),
fHistIncTracks_mel(0),
fHistIncTracks_sel(0),
fHistIncTracks_mka(0),
fHistIncTracks_ska(0),
fHistIncTracks_mpr(0),
fHistIncTracks_spr(0),
fHistJet_ptsub_v_area(0),
fHistJet_kin(0),
fHistJet_moms(0),
fHist_pi_DCAxy(0),
fHist_pr_DCAxy(0),
fHistJet_pi_DCAxy(0),
fHistJet_pr_DCAxy(0),
fHistMCTruth_TrackEff_Den_pi(0),
fHistMCTruth_TrackEff_Den_ka(0),
fHistMCTruth_TrackEff_Den_pr(0),
fHistMCReco_TrackEff_Num_pi(0),
fHistMCReco_TrackEff_Num_ka(0),
fHistMCReco_TrackEff_Num_pr(0),
fHistMCReco_MatchEff_Num_pi(0),
fHistMCReco_MatchEff_Num_ka(0),
fHistMCReco_MatchEff_Num_pr(0),
fHistMCReco_pi_prim_DCAxy(0),
fHistMCReco_pi_scdweak_DCAxy(0),
fHistMCReco_pr_prim_DCAxy(0),
fHistMCReco_pr_scdweak_DCAxy(0),
fHistMCReco_pr_scdmat_DCAxy(0),
fHist_Jet_MCTruth_TrackEff_Den_pi(0),
fHist_Jet_MCTruth_TrackEff_Den_ka(0),
fHist_Jet_MCTruth_TrackEff_Den_pr(0),
fHist_Jet_MCReco_TrackEff_Num_pi(0),
fHist_Jet_MCReco_TrackEff_Num_ka(0),
fHist_Jet_MCReco_TrackEff_Num_pr(0),
fHist_Jet_MCReco_MatchEff_Num_pi(0),
fHist_Jet_MCReco_MatchEff_Num_ka(0),
fHist_Jet_MCReco_MatchEff_Num_pr(0),
fHist_JetMatch_MCTruth_TrackEff_Den_pi(0),
fHist_JetMatch_MCTruth_TrackEff_Den_ka(0),
fHist_JetMatch_MCTruth_TrackEff_Den_pr(0),
fHist_JetMatch_MCReco_TrackEff_Num_pi(0),
fHist_JetMatch_MCReco_TrackEff_Num_ka(0),
fHist_JetMatch_MCReco_TrackEff_Num_pr(0),
fHist_JetMatch_MCReco_MatchEff_Num_pi(0),
fHist_JetMatch_MCReco_MatchEff_Num_ka(0),
fHist_JetMatch_MCReco_MatchEff_Num_pr(0),
fHistMCT_Jet_ptsub_v_area(0),
fHistMCR_Jet_ptsub_v_area(0),
fall_reco_jets_w_multiple_matches(0),
fall_reco_jets_w_matches(0),
fHistMC_Jet_shared_pt_frac(0),
fHistMC_Jet_deltaR(0),
fFastJetWrapper(0x0),
fFastJetWrapperBG(0x0),
fFastJetWrapper_Rec(0x0),
fFastJetWrapperBG_Rec(0x0),
fFastJetWrapper_Gen(0x0),
fFastJetWrapperBG_Gen(0x0),
fHist_TrueJetPtFraction(0),
fHist_MatchJetPts(0),
fHist_MatchJetEtas(0),
fHist_MatchJetDeltaRs(0),
fHist_JERS_PbPbunid_tru(0),
fHist_JERS_PbPbunid_det(0),
fHist_JERS_PbPbunid_truhyb(0),
fHist_JERS_PbPbunid_dethyb(0),
fHist_JERS_truPythia(0),
fHist_JERS_truPythia_pi(0),
fHist_JERS_truPythia_ka(0),
fHist_JERS_truPythia_pr(0),
fHist_JERS_detPythia(0),
fHist_JERS_detPythia_pi(0),
fHist_JERS_detPythia_ka(0),
fHist_JERS_detPythia_pr(0),
fHist_JERS_truhybPythia(0),
fHist_JERS_truhybPythia_pi(0),
fHist_JERS_truhybPythia_ka(0),
fHist_JERS_truhybPythia_pr(0),
fHist_JERS_dethybPythia(0),
fHist_JERS_dethybPythia_pi(0),
fHist_JERS_dethybPythia_ka(0),
fHist_JERS_dethybPythia_pr(0),
fHist_PC_spectra_Pythia(0),
fHist_PC_spectra_Pythia_pi(0),
fHist_PC_spectra_Pythia_ka(0),
fHist_PC_spectra_Pythia_pr(0),
fHist_det_spectra_Pythia(0),
fHist_det_spectra_Pythia_pi(0),
fHist_det_spectra_Pythia_ka(0),
fHist_det_spectra_Pythia_pr(0),
fHist_part_spectra_Pythia(0),
fHist_part_spectra_Pythia_pi(0),
fHist_part_spectra_Pythia_ka(0),
fHist_part_spectra_Pythia_pr(0),
fDoEmbedding(kFALSE),
fMultTrueFlag(kFALSE)
{
  //
  //         standard constructur which should be used
  //
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  std::cout << " Info::siweyhmi:***************** CONSTRUCTOR CALLED: AliAnalysisTaskJetHadroAOD  *****************"<< std::endl;
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  std::cout << " Info::siweyhmi:===================================================================================="<< std::endl;
  // ==========================================
  //
  // ==========================================
  //
  // Define outputs
  DefineOutput(1, TList::Class());

  if (fDoEventTree){
    DefineOutput(2, TTree::Class());
    DefineOutput(3, TTree::Class());
    DefineOutput(4, TTree::Class());
    DefineOutput(5, TTree::Class());
    DefineOutput(6, TTree::Class());
    DefineOutput(7, TTree::Class());
    DefineOutput(8, TTree::Class());
    DefineOutput(9, TTree::Class());
    DefineOutput(10, TTree::Class());
    DefineOutput(11, TTree::Class());
    DefineOutput(12, TTree::Class());
  }

  // ==========================================

}
//________________________________________________________________________
AliAnalysisTaskJetHadroAOD::~AliAnalysisTaskJetHadroAOD()
{

  //
  // Destructor
  //
  std::cout << " Info::siweyhmi: ===== In the Destructor ===== " << std::endl;
  if (fHistCentrality)      delete fHistCentrality;
  if (fHistImpParam)        delete fHistImpParam;
  if (fHistVertex)          delete fHistVertex;
  if (fHist_Is_Good_Inc_Event)          delete fHist_Is_Good_Inc_Event;
  if (fHist_Has_Acc_Real_FJ_Jet)          delete fHist_Has_Acc_Real_FJ_Jet;
  if (fHist_Has_Acc_Real_EMC_Jet)          delete fHist_Has_Acc_Real_EMC_Jet;
  if (fHist_Num_Real_FJ_Jets)          delete fHist_Num_Real_FJ_Jets;
  if (fHist_Num_Real_EMC_Jets)          delete fHist_Num_Real_EMC_Jets;
  if (fHist_RhoFJ)          delete fHist_RhoFJ;
  if (fHist_RhoEMC)          delete fHist_RhoEMC;
  if (fHist_RhoEMCReco)          delete fHist_RhoEMCReco;
  if (fHist_RhoEMCGen)          delete fHist_RhoEMCGen;
  if (fHist_RhoEMCEmb)          delete fHist_RhoEMCEmb;
  if (fHist_Filled_UEC_Emb)          delete fHist_Filled_UEC_Emb;
  if (fHist_Filled_UEC_Gen)          delete fHist_Filled_UEC_Gen;
  if (fHist_Filled_UEC_Rec)          delete fHist_Filled_UEC_Rec;

  if (fHistIncTracks_dEdx)          delete fHistIncTracks_dEdx;
  if (fHistIncTracks_moms)          delete fHistIncTracks_moms;
  if (fHistIncTracks_moms_p)          delete fHistIncTracks_moms_p;
  if (fHistIncTracks_moms_pTPC_p)          delete fHistIncTracks_moms_pTPC_p;
  if (fHistIncTracks_kin)          delete fHistIncTracks_kin;
  if (fHistIncTracks_beta)          delete fHistIncTracks_beta;
  if (fHistJetTracks_beta)          delete fHistJetTracks_beta;
  if (fHistIncTracks_TOFpi_nsigma)          delete fHistIncTracks_TOFpi_nsigma;
  if (fHistIncTracks_TOFka_nsigma)          delete fHistIncTracks_TOFka_nsigma;
  if (fHistIncTracks_TOFpr_nsigma)          delete fHistIncTracks_TOFpr_nsigma;
  if (fHistJetTracks_dEdx)          delete fHistJetTracks_dEdx;
  if (fHistJetTracks_moms)          delete fHistJetTracks_moms;
  if (fHistJetTracks_moms_p)          delete fHistJetTracks_moms_p;
  if (fHistJetTracks_moms_pTPC_p)          delete fHistJetTracks_moms_pTPC_p;
  if (fHistJetTracks_kin)          delete fHistJetTracks_kin;
  if (fHistIncTracks_t0)          delete fHistIncTracks_t0;
  if (fHistJetTracks_TOFpi_nsigma)          delete fHistJetTracks_TOFpi_nsigma;
  if (fHistJetTracks_TOFka_nsigma)          delete fHistJetTracks_TOFka_nsigma;
  if (fHistJetTracks_TOFpr_nsigma)          delete fHistJetTracks_TOFpr_nsigma;
  if (fHistBetaExpec_pi) delete fHistBetaExpec_pi;
  if (fHistBetaExpec_ka) delete fHistBetaExpec_ka;
  if (fHistBetaExpec_pr) delete fHistBetaExpec_pr;
  if (fHistjet_BetaExpec_pi) delete fHistjet_BetaExpec_pi;
  if (fHistjet_BetaExpec_ka) delete fHistjet_BetaExpec_ka;
  if (fHistjet_BetaExpec_pr) delete fHistjet_BetaExpec_pr;
  if (fHist_pi_mismatch) delete fHist_pi_mismatch;
  if (fHist_ka_mismatch) delete fHist_ka_mismatch;
  if (fHist_pr_mismatch) delete fHist_pr_mismatch;
  if (fHist_jet_pi_mismatch) delete fHist_jet_pi_mismatch;
  if (fHist_jet_ka_mismatch) delete fHist_jet_ka_mismatch;
  if (fHist_jet_pr_mismatch) delete fHist_jet_pr_mismatch;
  if (fHist_elExpec_pihyp) delete fHist_elExpec_pihyp;
  if (fHist_muExpec_pihyp) delete fHist_muExpec_pihyp;
  if (fHist_kaExpec_pihyp) delete fHist_kaExpec_pihyp;
  if (fHist_prExpec_pihyp) delete fHist_prExpec_pihyp;
  if (fHist_piExpec_kahyp) delete fHist_piExpec_kahyp;
  if (fHist_prExpec_kahyp) delete fHist_prExpec_kahyp;
  if (fHist_piExpec_prhyp) delete fHist_piExpec_prhyp;
  if (fHist_kaExpec_prhyp) delete fHist_kaExpec_prhyp;
  if (fHist_deExpec_prhyp) delete fHist_deExpec_prhyp;
  if (fHist_jet_elExpec_pihyp) delete fHist_jet_elExpec_pihyp;
  if (fHist_jet_muExpec_pihyp) delete fHist_jet_muExpec_pihyp;
  if (fHist_jet_kaExpec_pihyp) delete fHist_jet_kaExpec_pihyp;
  if (fHist_jet_prExpec_pihyp) delete fHist_jet_prExpec_pihyp;
  if (fHist_jet_piExpec_kahyp) delete fHist_jet_piExpec_kahyp;
  if (fHist_jet_prExpec_kahyp) delete fHist_jet_prExpec_kahyp;
  if (fHist_jet_piExpec_prhyp) delete fHist_jet_piExpec_prhyp;
  if (fHist_jet_kaExpec_prhyp) delete fHist_jet_kaExpec_prhyp;
  if (fHist_jet_deExpec_prhyp) delete fHist_jet_deExpec_prhyp;
  if (fHistTOFSigmaExpec_pi) delete fHistTOFSigmaExpec_pi;
  if (fHistTOFSigmaExpec_ka) delete fHistTOFSigmaExpec_ka;
  if (fHistTOFSigmaExpec_pr) delete fHistTOFSigmaExpec_pr;
  if (fHistjet_TOFSigmaExpec_pi) delete fHistjet_TOFSigmaExpec_pi;
  if (fHistjet_TOFSigmaExpec_ka) delete fHistjet_TOFSigmaExpec_ka;
  if (fHistjet_TOFSigmaExpec_pr) delete fHistjet_TOFSigmaExpec_pr;
  if (fHistIncTracks_mpi)          delete fHistIncTracks_mpi;
  if (fHistIncTracks_spi)          delete fHistIncTracks_spi;
  if (fHistIncTracks_mel)          delete fHistIncTracks_mel;
  if (fHistIncTracks_sel)          delete fHistIncTracks_sel;
  if (fHistIncTracks_mka)          delete fHistIncTracks_mka;
  if (fHistIncTracks_ska)          delete fHistIncTracks_ska;
  if (fHistIncTracks_mpr)          delete fHistIncTracks_mpr;
  if (fHistIncTracks_spr)          delete fHistIncTracks_spr;
  if (fHistJet_ptsub_v_area)          delete fHistJet_ptsub_v_area;
  if (fHistJet_kin)          delete fHistJet_kin;
  if (fHistJet_moms)          delete fHistJet_moms;
  if (fHist_pi_DCAxy)          delete fHist_pi_DCAxy;
  if (fHist_pr_DCAxy)          delete fHist_pr_DCAxy;
  if (fHistJet_pi_DCAxy)          delete fHistJet_pi_DCAxy;
  if (fHistJet_pr_DCAxy)          delete fHistJet_pr_DCAxy;
  if (fHistMCTruth_TrackEff_Den_pi)          delete fHistMCTruth_TrackEff_Den_pi;
  if (fHistMCTruth_TrackEff_Den_ka)          delete fHistMCTruth_TrackEff_Den_ka;
  if (fHistMCTruth_TrackEff_Den_pr)          delete fHistMCTruth_TrackEff_Den_pr;
  if (fHistMCReco_TrackEff_Num_pi)          delete fHistMCReco_TrackEff_Num_pi;
  if (fHistMCReco_TrackEff_Num_ka)          delete fHistMCReco_TrackEff_Num_ka;
  if (fHistMCReco_TrackEff_Num_pr)          delete fHistMCReco_TrackEff_Num_pr;
  if (fHistMCReco_MatchEff_Num_pi)          delete fHistMCReco_MatchEff_Num_pi;
  if (fHistMCReco_MatchEff_Num_ka)          delete fHistMCReco_MatchEff_Num_ka;
  if (fHistMCReco_MatchEff_Num_pr)          delete fHistMCReco_MatchEff_Num_pr;
  if (fHistMCReco_pi_prim_DCAxy)          delete fHistMCReco_pi_prim_DCAxy;
  if (fHistMCReco_pi_scdweak_DCAxy)          delete fHistMCReco_pi_scdweak_DCAxy;
  if (fHistMCReco_pr_prim_DCAxy)          delete fHistMCReco_pr_prim_DCAxy;
  if (fHistMCReco_pr_scdweak_DCAxy)          delete fHistMCReco_pr_scdweak_DCAxy;
  if (fHistMCReco_pr_scdmat_DCAxy)          delete fHistMCReco_pr_scdmat_DCAxy;
  if (fHist_Jet_MCTruth_TrackEff_Den_pi)          delete fHist_Jet_MCTruth_TrackEff_Den_pi;
  if (fHist_Jet_MCTruth_TrackEff_Den_ka)          delete fHist_Jet_MCTruth_TrackEff_Den_ka;
  if (fHist_Jet_MCTruth_TrackEff_Den_pr)          delete fHist_Jet_MCTruth_TrackEff_Den_pr;
  if (fHist_Jet_MCReco_TrackEff_Num_pi)          delete fHist_Jet_MCReco_TrackEff_Num_pi;
  if (fHist_Jet_MCReco_TrackEff_Num_ka)          delete fHist_Jet_MCReco_TrackEff_Num_ka;
  if (fHist_Jet_MCReco_TrackEff_Num_pr)          delete fHist_Jet_MCReco_TrackEff_Num_pr;
  if (fHist_Jet_MCReco_MatchEff_Num_pi)          delete fHist_Jet_MCReco_MatchEff_Num_pi;
  if (fHist_Jet_MCReco_MatchEff_Num_ka)          delete fHist_Jet_MCReco_MatchEff_Num_ka;
  if (fHist_Jet_MCReco_MatchEff_Num_pr)          delete fHist_Jet_MCReco_MatchEff_Num_pr;
  if (fHist_JetMatch_MCTruth_TrackEff_Den_pi)          delete fHist_JetMatch_MCTruth_TrackEff_Den_pi;
  if (fHist_JetMatch_MCTruth_TrackEff_Den_ka)          delete fHist_JetMatch_MCTruth_TrackEff_Den_ka;
  if (fHist_JetMatch_MCTruth_TrackEff_Den_pr)          delete fHist_JetMatch_MCTruth_TrackEff_Den_pr;
  if (fHist_JetMatch_MCReco_TrackEff_Num_pi)          delete fHist_JetMatch_MCReco_TrackEff_Num_pi;
  if (fHist_JetMatch_MCReco_TrackEff_Num_ka)          delete fHist_JetMatch_MCReco_TrackEff_Num_ka;
  if (fHist_JetMatch_MCReco_TrackEff_Num_pr)          delete fHist_JetMatch_MCReco_TrackEff_Num_pr;
  if (fHist_JetMatch_MCReco_MatchEff_Num_pi)          delete fHist_JetMatch_MCReco_MatchEff_Num_pi;
  if (fHist_JetMatch_MCReco_MatchEff_Num_ka)          delete fHist_JetMatch_MCReco_MatchEff_Num_ka;
  if (fHist_JetMatch_MCReco_MatchEff_Num_pr)          delete fHist_JetMatch_MCReco_MatchEff_Num_pr;
  if (fHistMCT_Jet_ptsub_v_area)          delete fHistMCT_Jet_ptsub_v_area;
  if (fHistMCR_Jet_ptsub_v_area)          delete fHistMCR_Jet_ptsub_v_area;
  if (fHistMC_Jet_shared_pt_frac)          delete fHistMC_Jet_shared_pt_frac;
  if (fHistMC_Jet_deltaR)          delete fHistMC_Jet_deltaR;

  if (fHist_TrueJetPtFraction)          delete fHist_TrueJetPtFraction;
  if (fHist_MatchJetPts)          delete fHist_MatchJetPts;
  if (fHist_MatchJetEtas)          delete fHist_MatchJetEtas;
  if (fHist_MatchJetDeltaRs)          delete fHist_MatchJetDeltaRs;

  if (fHist_JERS_PbPbunid_tru)          delete fHist_JERS_PbPbunid_tru;
  if (fHist_JERS_PbPbunid_det)          delete fHist_JERS_PbPbunid_det;
  if (fHist_JERS_PbPbunid_truhyb)          delete fHist_JERS_PbPbunid_truhyb;
  if (fHist_JERS_PbPbunid_dethyb)          delete fHist_JERS_PbPbunid_dethyb;

  if (fHist_JERS_truPythia)          delete fHist_JERS_truPythia;
  if (fHist_JERS_truPythia_pi)          delete fHist_JERS_truPythia_pi;
  if (fHist_JERS_truPythia_ka)          delete fHist_JERS_truPythia_ka;
  if (fHist_JERS_truPythia_pr)          delete fHist_JERS_truPythia_pr;

  if (fHist_JERS_detPythia)          delete fHist_JERS_detPythia;
  if (fHist_JERS_detPythia_pi)          delete fHist_JERS_detPythia_pi;
  if (fHist_JERS_detPythia_ka)          delete fHist_JERS_detPythia_ka;
  if (fHist_JERS_detPythia_pr)          delete fHist_JERS_detPythia_pr;

  if (fHist_JERS_truhybPythia)          delete fHist_JERS_truhybPythia;
  if (fHist_JERS_truhybPythia_pi)          delete fHist_JERS_truhybPythia_pi;
  if (fHist_JERS_truhybPythia_ka)          delete fHist_JERS_truhybPythia_ka;
  if (fHist_JERS_truhybPythia_pr)          delete fHist_JERS_truhybPythia_pr;

  if (fHist_JERS_dethybPythia)          delete fHist_JERS_dethybPythia;
  if (fHist_JERS_dethybPythia_pi)          delete fHist_JERS_dethybPythia_pi;
  if (fHist_JERS_dethybPythia_ka)          delete fHist_JERS_dethybPythia_ka;
  if (fHist_JERS_dethybPythia_pr)          delete fHist_JERS_dethybPythia_pr;

  if (fHist_PC_spectra_Pythia)          delete fHist_PC_spectra_Pythia;
  if (fHist_PC_spectra_Pythia_pi)          delete fHist_PC_spectra_Pythia_pi;
  if (fHist_PC_spectra_Pythia_ka)          delete fHist_PC_spectra_Pythia_ka;
  if (fHist_PC_spectra_Pythia_pr)          delete fHist_PC_spectra_Pythia_pr;

  if (fHist_det_spectra_Pythia)          delete fHist_det_spectra_Pythia;
  if (fHist_det_spectra_Pythia_pi)          delete fHist_det_spectra_Pythia_pi;
  if (fHist_det_spectra_Pythia_ka)          delete fHist_det_spectra_Pythia_ka;
  if (fHist_det_spectra_Pythia_pr)          delete fHist_det_spectra_Pythia_pr;

  if (fHist_part_spectra_Pythia)          delete fHist_part_spectra_Pythia;
  if (fHist_part_spectra_Pythia_pi)          delete fHist_part_spectra_Pythia_pi;
  if (fHist_part_spectra_Pythia_ka)          delete fHist_part_spectra_Pythia_ka;
  if (fHist_part_spectra_Pythia_pr)          delete fHist_part_spectra_Pythia_pr;

  if (fPIDCombined) delete fPIDCombined;
  if (fTreeSRedirector)       delete fTreeSRedirector;

  if (fFastJetWrapper) delete fFastJetWrapper;
  if (fFastJetWrapperBG) delete fFastJetWrapperBG;
  if (fFastJetWrapper_Rec) delete fFastJetWrapper_Rec;
  if (fFastJetWrapperBG_Rec) delete fFastJetWrapperBG_Rec;
  if (fFastJetWrapper_Gen) delete fFastJetWrapper_Gen;
  if (fFastJetWrapperBG_Gen) delete fFastJetWrapperBG_Gen;
  if (fRandom) delete fRandom;

}
//
// ---------------------------------------------------------------------------------
//                                     Functions
// ---------------------------------------------------------------------------------
//
void AliAnalysisTaskJetHadroAOD::Initialize()
{
  //
  std::cout << " Info::siweyhmi: ===== In the Initialize ===== " << std::endl;
  //
  // ------------------------------------------------
  //
  //

  //
  //
  std::cout << " Info::siweyhmi: ===================================================== " << std::endl;
  std::cout << " Info::siweyhmi: =============== Summary of Track Cuts =============== " << std::endl;
  std::cout << " Info::siweyhmi: ===================================================== " << std::endl;
}
//________________________________________________________________________
void AliAnalysisTaskJetHadroAOD::UserCreateOutputObjects()
{
  //
  // Create output histograms, trees of the analysis (called once)
  //
  Initialize();
  std::cout << " Info::siweyhmi: ===== In the UserCreateOutputObjects ===== " << std::endl;

  fFastJetWrapper = new AliFJWrapper("fFastJetWrapper","fFastJetWrapper");
  fFastJetWrapper->Clear();

  fFastJetWrapperBG = new AliFJWrapper("fFastJetWrapperBG","fFastJetWrapperBG");
  fFastJetWrapperBG->Clear();

  fFastJetWrapper_Rec = new AliFJWrapper("fFastJetWrapper_Rec","fFastJetWrapper_Rec");
  fFastJetWrapper_Rec->Clear();

  fFastJetWrapperBG_Rec = new AliFJWrapper("fFastJetWrapperBG_Rec","fFastJetWrapperBG_Rec");
  fFastJetWrapperBG_Rec->Clear();

  fFastJetWrapper_Gen = new AliFJWrapper("fFastJetWrapper_Gen","fFastJetWrapper_Gen");
  fFastJetWrapper_Gen->Clear();

  fFastJetWrapperBG_Gen = new AliFJWrapper("fFastJetWrapperBG_Gen","fFastJetWrapperBG_Gen");
  fFastJetWrapperBG_Gen->Clear();

  fRandom = new TRandom3(0);

  for(Int_t iCont=0; iCont<fParticleCollArray.GetEntriesFast(); iCont++){
    if (GetParticleContainer(iCont)->GetIsEmbedding()){
      fIsEmbeddedEvent = kTRUE;
    }
  }

  if (fDoEmbedding)
  {
    if(!fIsEmbeddedEvent){
      fMCParticleArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fMCParticleArrayName.Data()));
    }
    else
    {
      // In case of embedding, the MC particle array needs to be fetched differently
      AliVEvent* event = AliEmcalContainerUtils::GetEvent(InputEvent(), kTRUE);
      fMCParticleArray = dynamic_cast<TClonesArray*>(event->FindListObject(fMCParticleArrayName.Data()));
    }
  }


  // ------------  setup PIDCombined  ---------------
  fPIDCombined = new AliPIDCombined("pidCombined","");
  //
  // **********************   Input handler to get the PID object *********************
  AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if (!inputHandler)
  AliFatal("Input handler needed");
  else {
    fPIDResponse = inputHandler->GetPIDResponse();       // PID response object
    if (!fPIDResponse) std::cout << " Info::siweyhmi: ======= PIDResponse object was not created ====== " << std::endl;
  }
  //
  // ************************************************************************
  //   OpenFile output --> one can open several files
  // ************************************************************************
  //
  OpenFile(1);
  fTreeSRedirector = new TTreeSRedirector();
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);
  //
  fEventCuts.AddQAplotsToList(fListHist);

  Float_t tpc_mom_bins[58] = {};
  for (int i=0; i < (57+1); i++){
    if (i < 20) tpc_mom_bins[i] = 0.0 + i*0.05;
    else if (i < 30) tpc_mom_bins[i] = 1.0 + (i-20)*0.1;
    else if (i < 40) tpc_mom_bins[i] = 2.0 + (i-30)*0.2;
    else if (i < 46) tpc_mom_bins[i] = 4.0 + (i-40)*0.5;
    else if (i < 55) tpc_mom_bins[i] = 7.0 + (i-46)*1.0;
    else tpc_mom_bins[i] = 16.0 + (i-55)*2.0;
  }
  Int_t n_tpc_mom_bins = sizeof(tpc_mom_bins)/sizeof(Float_t) - 1;

  Float_t TOF_mom_bins[48] = {};
  for (int i=0; i < (47+1); i++){
    if (i < 20) TOF_mom_bins[i] = 0.0 + i*0.05;
    else if (i < 30) TOF_mom_bins[i] = 1.0 + (i-20)*0.1;
    else if (i < 40) TOF_mom_bins[i] = 2.0 + (i-30)*0.2;
    else if (i < 46) TOF_mom_bins[i] = 4.0 + (i-40)*0.5;
    else TOF_mom_bins[i] = 7.0 + (i-46)*1.0;
  }
  Int_t n_TOF_mom_bins = sizeof(TOF_mom_bins)/sizeof(Float_t) - 1;

  Float_t dedx_bins[2001] = {};
  for (int i=0; i < (2000+1); i++){
    dedx_bins[i] = 0.0 + i*1.0;
  }
  Int_t n_dedx_bins = sizeof(dedx_bins)/sizeof(Float_t) - 1;

  Float_t beta_bins[1001] = {};
  for (int i=0; i < (1000+1); i++){
    beta_bins[i] = 0.0 + i*0.002;
  }
  Int_t n_beta_bins = sizeof(beta_bins)/sizeof(Float_t) - 1;

  Float_t tof_n_sigma_bins[2001] = {};
  for (int i=0; i < (2000+1); i++){
    tof_n_sigma_bins[i] = -100.0 + i*0.1;
  }
  Int_t n_tof_n_sigma_bins = sizeof(tof_n_sigma_bins)/sizeof(Float_t) - 1;

  Float_t eta_bins[5] = {};
  for (int i=0; i < (4+1); i++){
    if (i < 3) eta_bins[i] = 0.0 + i*0.2;
    else eta_bins[i] = 0.6 + (i-3)*0.3;
  }
  Int_t n_eta_bins = sizeof(eta_bins)/sizeof(Float_t) - 1;

  Float_t dcaxy_bins[2401] = {};
  for (int i=0; i < (2400+1); i++){
    dcaxy_bins[i] = -3.0 + i*0.0025; //600 from -3 to 3
  }
  Int_t n_dcaxy_bins = sizeof(dcaxy_bins)/sizeof(Float_t) - 1;

  //
  // ************************************************************************
  //   Event histograms
  // ************************************************************************
  //
  fHistCentrality        = new TH1F("hCentrality",           "control histogram for centrality"           , 11,  -10., 100.);
  fHistVertex            = new TH1F("hVertex",               "control histogram for vertex Z position"    , 200, -20., 20.);

  fHist_Is_Good_Inc_Event            = new TH1F("hIs_Good_Inc_Event", "Good inclusive event or not"    , 2, 0., 2.);
  fHist_Has_Acc_Real_FJ_Jet            = new TH1F("hHas_Acc_Real_FJ_Jet", "Good jet event: 0 no accepted jets, 1 accepted jet but not real, 2 has real jet" , 3, 0., 3.);
  fHist_Has_Acc_Real_EMC_Jet            = new TH1F("hHas_Acc_Real_EMC_Jet", "Good jet event: 0 no accepted jets, 1 accepted jet but not real, 2 has real jet" , 3, 0., 3.);
  fHist_Num_Real_FJ_Jets            = new TH1F("hNum_Real_FJ_Jets", "Number of Real FJ Jets", 10, 0., 10.);
  fHist_Num_Real_EMC_Jets            = new TH1F("hNum_Real_EMC_Jets", "Number of Real FJ Jets", 10, 0., 10.);
  fHist_RhoFJ            = new TH1F("hRhoFJ", "FJ Rho Value"    , 600, 0., 600.);
  fHist_RhoEMC            = new TH1F("hRhoEMC", "EMC Rho Value"    , 600, 0., 600.);
  fHist_RhoEMCReco            = new TH1F("hRhoEMCReco", "EMC Reco Rho Value", 600, 0., 600.);
  fHist_RhoEMCGen            = new TH1F("hRhoEMCGen", "EMC Gen Rho Value", 600, 0., 600.);
  fHist_RhoEMCEmb            = new TH1F("hRhoEMCEmb", "EMC Emb Rho Value", 600, 0., 600.);

  fHist_Filled_UEC_Emb            = new TH1F("hFilled_UEC_Emb", "Filled Emb UE Cone or not", 2, 0., 2.);
  fHist_Filled_UEC_Gen            = new TH1F("hFilled_UEC_Gen", "Filled Gen UE Cone or not", 2, 0., 2.);
  fHist_Filled_UEC_Rec            = new TH1F("hFilled_UEC_Rec", "Filled Rec UE Cone or not", 2, 0., 2.);

  if (!fMCtrue || fDoEmbedding){

    if (fFill_TPC) {
      fHistIncTracks_dEdx    = new TH3F("fHistIncTracks_dEdx",   "dEdx histogram for inclusive tracks"        , n_tpc_mom_bins,tpc_mom_bins,n_dedx_bins,dedx_bins, n_eta_bins,eta_bins);
      fHistJetTracks_dEdx    = new TH3F("fHistJetTracks_dEdx",   "dEdx histogram for Jet tracks"            , n_tpc_mom_bins,tpc_mom_bins,n_dedx_bins,dedx_bins, n_eta_bins,eta_bins);
    }

    if (fFillpTPC_pT){
      fHistIncTracks_moms    = new TH2F("fHistIncTracks_moms",   "All mom types for inclusive tracks"         , n_tpc_mom_bins,tpc_mom_bins, n_tpc_mom_bins,tpc_mom_bins);
      fHistJetTracks_moms    = new TH2F("fHistJetTracks_moms",   "All mom types for Jet tracks"         , n_tpc_mom_bins,tpc_mom_bins, n_tpc_mom_bins,tpc_mom_bins);
    }
    if (fFillp_pT) {
      fHistIncTracks_moms_p    = new TH2F("fHistIncTracks_moms_p",   "All mom types for inclusive tracks"         , n_tpc_mom_bins,tpc_mom_bins, n_tpc_mom_bins,tpc_mom_bins);
      fHistJetTracks_moms_p    = new TH2F("fHistJetTracks_moms_p",   "All mom types for Jet tracks"         , n_tpc_mom_bins,tpc_mom_bins, n_tpc_mom_bins,tpc_mom_bins);
    }
    if (fFillpTPC_p) {
      fHistIncTracks_moms_pTPC_p    = new TH2F("fHistIncTracks_moms_pTPC_p",   "All mom types for inclusive tracks"         , n_tpc_mom_bins,tpc_mom_bins, n_tpc_mom_bins,tpc_mom_bins);
      fHistJetTracks_moms_pTPC_p    = new TH2F("fHistJetTracks_moms_pTPC_p",   "All mom types for Jet tracks"         , n_tpc_mom_bins,tpc_mom_bins, n_tpc_mom_bins,tpc_mom_bins);
    }

    fHistIncTracks_kin    = new TH3F("fHistIncTracks_kin",     "Kinematics histogram for inclusive tracks"  , 200, 0., 20., 9, 0.0, 0.9, 64, -3.2, 3.2);

    if (fFill_TOF){
      fHistIncTracks_beta    = new TH2F("fHistIncTracks_beta",   "Beta histogram for inclusive tracks"        , n_TOF_mom_bins,TOF_mom_bins,   n_beta_bins,beta_bins);
      fHistIncTracks_TOFpi_nsigma    = new TH3F("fHistIncTracks_TOFpi_nsigma",   "TOF Nsigma histogram for inclusive tracks under the pion hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHistIncTracks_TOFka_nsigma    = new TH3F("fHistIncTracks_TOFka_nsigma",   "TOF Nsigma histogram for inclusive tracks under the kaon hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHistIncTracks_TOFpr_nsigma    = new TH3F("fHistIncTracks_TOFpr_nsigma",   "TOF Nsigma histogram for inclusive tracks under the proton hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);

      fHistJetTracks_beta    = new TH2F("fHistJetTracks_beta",   "Beta histogram for jet tracks"        , n_TOF_mom_bins,TOF_mom_bins,   n_beta_bins,beta_bins);
      fHistJetTracks_TOFpi_nsigma    = new TH3F("fHistJetTracks_TOFpi_nsigma",   "TOF Nsigma histogram for jet tracks under the pion hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHistJetTracks_TOFka_nsigma    = new TH3F("fHistJetTracks_TOFka_nsigma",   "TOF Nsigma histogram for jet tracks under the kaon hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHistJetTracks_TOFpr_nsigma    = new TH3F("fHistJetTracks_TOFpr_nsigma",   "TOF Nsigma histogram for jet tracks under the proton hypothesis"        , n_TOF_mom_bins,TOF_mom_bins,   n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);

      Float_t t0_bins[141];
      for (int i=0; i < (140+1); i++){
        t0_bins[i] = -700.0 + i*10.0;
      }
      fHistIncTracks_t0   = new TH2F("fHistIncTracks_t0",   "TO histogram for inc tracks", n_TOF_mom_bins,TOF_mom_bins,   140,t0_bins);
    }

    fHistJetTracks_kin    = new TH3F("fHistJetTracks_kin",     "Kinematics histogram for Jet tracks"  , 200, 0., 20., 9, 0., 0.9, 64, -3.2, 3.2);

    if (fFill_TOF_expecs){
      fHistBetaExpec_pi    = new TH2F("fHistBetaExpec_pi",   "Expected Beta histogram for pions", n_TOF_mom_bins,TOF_mom_bins, n_beta_bins,beta_bins);
      fHistBetaExpec_ka    = new TH2F("fHistBetaExpec_ka",   "Expected Beta histogram for kaons", n_TOF_mom_bins,TOF_mom_bins,  n_beta_bins,beta_bins);
      fHistBetaExpec_pr    = new TH2F("fHistBetaExpec_pr",   "Expected Beta histogram for protons", n_TOF_mom_bins,TOF_mom_bins,  n_beta_bins,beta_bins);

      fHistjet_BetaExpec_pi    = new TH2F("fHistjet_BetaExpec_pi",   "Expected Beta histogram for pions", n_TOF_mom_bins,TOF_mom_bins, n_beta_bins,beta_bins);
      fHistjet_BetaExpec_ka    = new TH2F("fHistjet_BetaExpec_ka",   "Expected Beta histogram for kaons", n_TOF_mom_bins,TOF_mom_bins,  n_beta_bins,beta_bins);
      fHistjet_BetaExpec_pr    = new TH2F("fHistjet_BetaExpec_pr",   "Expected Beta histogram for protons", n_TOF_mom_bins,TOF_mom_bins, n_beta_bins,beta_bins);

      fHist_pi_mismatch    = new TH3F("fHist_pi_mismatch",   "Pion mismatch", n_TOF_mom_bins,TOF_mom_bins,  n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_ka_mismatch    = new TH3F("fHist_ka_mismatch",   "Kaon mismatch", n_TOF_mom_bins,TOF_mom_bins,  n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_pr_mismatch    = new TH3F("fHist_pr_mismatch",   "Proton mismatch", n_TOF_mom_bins,TOF_mom_bins,  n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);

      fHist_jet_pi_mismatch    = new TH3F("fHist_jet_pi_mismatch",   "Jet Pion mismatch", n_TOF_mom_bins,TOF_mom_bins,  n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_jet_ka_mismatch    = new TH3F("fHist_jet_ka_mismatch",   "Jet Kaon mismatch", n_TOF_mom_bins,TOF_mom_bins,  n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_jet_pr_mismatch    = new TH3F("fHist_jet_pr_mismatch",   "Jet Proton mismatch", n_TOF_mom_bins,TOF_mom_bins,  n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);

      fHist_elExpec_pihyp    = new TH3F("fHist_elExpec_pihyp",   "Expected Nsigma histogram for electrons under the pion hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_muExpec_pihyp    = new TH3F("fHist_muExpec_pihyp",   "Expected Nsigma histogram for muons under the pion hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_kaExpec_pihyp    = new TH3F("fHist_kaExpec_pihyp",   "Expected Nsigma histogram for kaons under the pion hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_prExpec_pihyp    = new TH3F("fHist_prExpec_pihyp",   "Expected Nsigma histogram for protons under the pion hypothesis", n_TOF_mom_bins,TOF_mom_bins,  n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_piExpec_kahyp    = new TH3F("fHist_piExpec_kahyp",   "Expected Nsigma histogram for pions under the kaon hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_prExpec_kahyp    = new TH3F("fHist_prExpec_kahyp",   "Expected Nsigma histogram for protons under the kaon hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins,n_eta_bins,eta_bins);
      fHist_piExpec_prhyp    = new TH3F("fHist_piExpec_prhyp",   "Expected Nsigma histogram for pions under the proton hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_kaExpec_prhyp    = new TH3F("fHist_kaExpec_prhyp",   "Expected Nsigma histogram for kaons under the proton hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_deExpec_prhyp    = new TH3F("fHist_deExpec_prhyp",   "Expected Nsigma histogram for deuterons under the proton hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);

      fHist_jet_elExpec_pihyp    = new TH3F("fHist_jet_elExpec_pihyp",   "Expected Nsigma histogram for electrons under the pion hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_jet_muExpec_pihyp    = new TH3F("fHist_jet_muExpec_pihyp",   "Expected Nsigma histogram for muons under the pion hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_jet_kaExpec_pihyp    = new TH3F("fHist_jet_kaExpec_pihyp",   "Expected Nsigma histogram for kaons under the pion hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_jet_prExpec_pihyp    = new TH3F("fHist_jet_prExpec_pihyp",   "Expected Nsigma histogram for protons under the pion hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_jet_piExpec_kahyp    = new TH3F("fHist_jet_piExpec_kahyp",   "Expected Nsigma histogram for pions under the kaon hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_jet_prExpec_kahyp    = new TH3F("fHist_jet_prExpec_kahyp",   "Expected Nsigma histogram for protons under the kaon hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_jet_piExpec_prhyp    = new TH3F("fHist_jet_piExpec_prhyp",   "Expected Nsigma histogram for pions under the proton hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_jet_kaExpec_prhyp    = new TH3F("fHist_jet_kaExpec_prhyp",   "Expected Nsigma histogram for kaons under the proton hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);
      fHist_jet_deExpec_prhyp    = new TH3F("fHist_jet_deExpec_prhyp",   "Expected Nsigma histogram for deuterons under the proton hypothesis", n_TOF_mom_bins,TOF_mom_bins, n_tof_n_sigma_bins,tof_n_sigma_bins, n_eta_bins,eta_bins);

      Float_t pi_expec_bins[401];
      for (int i=0; i < (400+1); i++){
        pi_expec_bins[i] = 60.0 + i*0.1;
      }
      Float_t ka_expec_bins[5401];
      for (int i=0; i < (5400+1); i++){
        ka_expec_bins[i] = 60.0 + i*0.1;
      }
      Float_t pr_expec_bins[5401];
      for (int i=0; i < (5400+1); i++){
        pr_expec_bins[i] = 60.0 + i*0.1;
      }

      fHistTOFSigmaExpec_pi    = new TH3F("fHistTOFSigmaExpec_pi",   "Expected TOF Sigma histogram for pions", n_TOF_mom_bins,TOF_mom_bins,   400,pi_expec_bins, n_eta_bins,eta_bins);
      fHistTOFSigmaExpec_ka    = new TH3F("fHistTOFSigmaExpec_ka",   "Expected TOF Sigma histogram for kaons", n_TOF_mom_bins,TOF_mom_bins,   5400,ka_expec_bins, n_eta_bins,eta_bins);
      fHistTOFSigmaExpec_pr    = new TH3F("fHistTOFSigmaExpec_pr",   "Expected TOF Sigma histogram for protons", n_TOF_mom_bins,TOF_mom_bins,   5400,pr_expec_bins, n_eta_bins,eta_bins);

      fHistjet_TOFSigmaExpec_pi    = new TH3F("fHistjet_TOFSigmaExpec_pi",   "Expected TOF Sigma histogram for pions", n_TOF_mom_bins,TOF_mom_bins,  400,pi_expec_bins, n_eta_bins,eta_bins);
      fHistjet_TOFSigmaExpec_ka    = new TH3F("fHistjet_TOFSigmaExpec_ka",   "Expected TOF Sigma histogram for kaons", n_TOF_mom_bins,TOF_mom_bins,   5400,ka_expec_bins, n_eta_bins,eta_bins);
      fHistjet_TOFSigmaExpec_pr    = new TH3F("fHistjet_TOFSigmaExpec_pr",   "Expected TOF Sigma histogram for protons", n_TOF_mom_bins,TOF_mom_bins,   5400,pr_expec_bins, n_eta_bins,eta_bins);
    }

    if (fFill_TPC_expecs){
      Float_t mean_expec_bins[361];
      for (int i=0; i < (360+1); i++){
        mean_expec_bins[i] = 40.0 + i*1.0;
      }
      Float_t sigma_expec_bins[26];
      for (int i=0; i < (25+1); i++){
        sigma_expec_bins[i] = 0.0 + i*1.0;
      }

      fHistIncTracks_mpi  = new TH3F("fHistIncTracks_mpi",     "Expected mean pion histogram for inclusive tracks"  , n_tpc_mom_bins,tpc_mom_bins, 360, mean_expec_bins, n_eta_bins,eta_bins);
      fHistIncTracks_spi  = new TH3F("fHistIncTracks_spi",     "Expected sigma pion histogram for inclusive tracks"  , n_tpc_mom_bins,tpc_mom_bins, 25, sigma_expec_bins, n_eta_bins,eta_bins);

      fHistIncTracks_mel  = new TH3F("fHistIncTracks_mel",     "Expected mean electron histogram for inclusive tracks"  , n_tpc_mom_bins,tpc_mom_bins, 360, mean_expec_bins, n_eta_bins,eta_bins);
      fHistIncTracks_sel  = new TH3F("fHistIncTracks_sel",     "Expected sigma electron histogram for inclusive tracks"  , n_tpc_mom_bins,tpc_mom_bins, 25, sigma_expec_bins, n_eta_bins,eta_bins);

      fHistIncTracks_mka  = new TH3F("fHistIncTracks_mka",     "Expected mean kaon histogram for inclusive tracks"  , n_tpc_mom_bins,tpc_mom_bins, 360, mean_expec_bins, n_eta_bins,eta_bins);
      fHistIncTracks_ska  = new TH3F("fHistIncTracks_ska",     "Expected sigma kaon histogram for inclusive tracks"  , n_tpc_mom_bins,tpc_mom_bins, 25, sigma_expec_bins, n_eta_bins,eta_bins);

      fHistIncTracks_mpr  = new TH3F("fHistIncTracks_mpr",     "Expected mean proton histogram for inclusive tracks"  , n_tpc_mom_bins,tpc_mom_bins, 360, mean_expec_bins, n_eta_bins,eta_bins);
      fHistIncTracks_spr  = new TH3F("fHistIncTracks_spr",     "Expected sigma proton histogram for inclusive tracks"  , n_tpc_mom_bins,tpc_mom_bins, 25, sigma_expec_bins, n_eta_bins,eta_bins);
    }

    fHistJet_ptsub_v_area  = new TH2F("fHistJet_ptsub_v_area", "Before cuts, Jet pt after subtraction vs jet area"  , 100, 0., 1., 300, 0., 300.);
    fHistJet_kin  = new TH3F("fHistJet_kin", "Kinematics histogram for Jets"  , 300, 0., 300., 48, -0.6, 0.6, 130, 0., 6.5);
    fHistJet_moms  = new TH2F("fHistJet_moms", "All mom types for jets"  , 300, 0., 300., 600, 0., 600.);

    fHist_pi_DCAxy  = new TH2F("fHist_pi_DCAxy", "Data DCAxy distribution pion", n_tpc_mom_bins,tpc_mom_bins, n_dcaxy_bins,dcaxy_bins);
    fHist_pr_DCAxy  = new TH2F("fHist_pr_DCAxy", "Data DCAxy distribution proton", n_tpc_mom_bins,tpc_mom_bins, n_dcaxy_bins,dcaxy_bins);

    fHistJet_pi_DCAxy  = new TH2F("fHistJet_pi_DCAxy", "Data DCAxy distribution pion", n_tpc_mom_bins,tpc_mom_bins, n_dcaxy_bins,dcaxy_bins);
    fHistJet_pr_DCAxy  = new TH2F("fHistJet_pr_DCAxy", "Data DCAxy distribution proton", n_tpc_mom_bins,tpc_mom_bins, n_dcaxy_bins,dcaxy_bins);
  }

  Float_t jet_pT_bins[49] = {};
  for (int i=0; i < (48+1); i++){
    if (i< 15) jet_pT_bins[i] = -300.0 + i*20.0;
    else if (i< 39) jet_pT_bins[i] = 0.0 + (i-15)*5.0;
    else jet_pT_bins[i] = 120.0 + (i-39)*20.0; 
    //size 20 GeV bins from -300 to 0, size 5 GeV bins from 0 to 120, size 20 GeV bins from 120 to 300
  }
  Int_t n_jet_pT_bins = sizeof(jet_pT_bins)/sizeof(Float_t) - 1;

  if (fMCtrue && !fDoEmbedding){

    fHistMCTruth_TrackEff_Den_pi  = new TH2F("fHistMCTruth_TrackEff_Den_pi", "MC Truth Level Tracking Efficiency Denominator pion", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins); //ensure the bins also work for TOF
    fHistMCTruth_TrackEff_Den_ka  = new TH2F("fHistMCTruth_TrackEff_Den_ka", "MC Truth Level Tracking Efficiency Denominator kaon", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins); //ensure the bins also work for TOF
    fHistMCTruth_TrackEff_Den_pr  = new TH2F("fHistMCTruth_TrackEff_Den_pr", "MC Truth Level Tracking Efficiency Denominator proton", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins); //ensure the bins also work for TOF
    fHistMCReco_TrackEff_Num_pi  = new TH2F("fHistMCReco_TrackEff_Num_pi", "MC Reconstructed Level Tracking Efficiency Numerator pion", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins);
    fHistMCReco_TrackEff_Num_ka  = new TH2F("fHistMCReco_TrackEff_Num_ka", "MC Reconstructed Level Tracking Efficiency Numerator kaon", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins);
    fHistMCReco_TrackEff_Num_pr  = new TH2F("fHistMCReco_TrackEff_Num_pr", "MC Reconstructed Level Tracking Efficiency Numerator proton", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins);
    fHistMCReco_MatchEff_Num_pi  = new TH2F("fHistMCReco_MatchEff_Num_pi", "MC Reconstructed Level TOF Matching Efficiency Numerator pion", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins);
    fHistMCReco_MatchEff_Num_ka  = new TH2F("fHistMCReco_MatchEff_Num_ka", "MC Reconstructed Level TOF Matching Efficiency Numerator kaon", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins);
    fHistMCReco_MatchEff_Num_pr  = new TH2F("fHistMCReco_MatchEff_Num_pr", "MC Reconstructed Level TOF Matching Efficiency Numerator proton", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins);

    fHistMCReco_pi_prim_DCAxy  = new TH2F("fHistMCReco_pi_prim_DCAxy", "MC Reconstructed Level Primary DCAxy dist pion", n_tpc_mom_bins,tpc_mom_bins, n_dcaxy_bins,dcaxy_bins);
    fHistMCReco_pi_scdweak_DCAxy  = new TH2F("fHistMCReco_pi_scdweak_DCAxy", "MC Reconstructed Level Secondary from Weak Decay DCAxy dist pion", n_tpc_mom_bins,tpc_mom_bins, n_dcaxy_bins,dcaxy_bins);
    fHistMCReco_pr_prim_DCAxy  = new TH2F("fHistMCReco_pr_prim_DCAxy", "MC Reconstructed Level Primary DCAxy dist pion", n_tpc_mom_bins,tpc_mom_bins, n_dcaxy_bins,dcaxy_bins);
    fHistMCReco_pr_scdweak_DCAxy  = new TH2F("fHistMCReco_pr_scdweak_DCAxy", "MC Reconstructed Level Secondary from Weak Decay DCAxy dist proton", n_tpc_mom_bins,tpc_mom_bins, n_dcaxy_bins,dcaxy_bins);
    fHistMCReco_pr_scdmat_DCAxy  = new TH2F("fHistMCReco_pr_scdmat_DCAxy", "MC Reconstructed Level Secondary from Material DCAxy dist proton", n_tpc_mom_bins,tpc_mom_bins, n_dcaxy_bins,dcaxy_bins);

    fHist_Jet_MCTruth_TrackEff_Den_pi  = new TH2F("fHist_Jet_MCTruth_TrackEff_Den_pi", "MC Truth Level Tracking Efficiency Denominator pion", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins); //ensure the bins also work for TOF
    fHist_Jet_MCTruth_TrackEff_Den_ka  = new TH2F("fHist_Jet_MCTruth_TrackEff_Den_ka", "MC Truth Level Tracking Efficiency Denominator kaon", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins); //ensure the bins also work for TOF
    fHist_Jet_MCTruth_TrackEff_Den_pr  = new TH2F("fHist_Jet_MCTruth_TrackEff_Den_pr", "MC Truth Level Tracking Efficiency Denominator proton", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins); //ensure the bins also work for TOF
    fHist_Jet_MCReco_TrackEff_Num_pi  = new TH2F("fHist_Jet_MCReco_TrackEff_Num_pi", "MC Reconstructed Level Tracking Efficiency Numerator pion", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins);
    fHist_Jet_MCReco_TrackEff_Num_ka  = new TH2F("fHist_Jet_MCReco_TrackEff_Num_ka", "MC Reconstructed Level Tracking Efficiency Numerator kaon", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins);
    fHist_Jet_MCReco_TrackEff_Num_pr  = new TH2F("fHist_Jet_MCReco_TrackEff_Num_pr", "MC Reconstructed Level Tracking Efficiency Numerator proton", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins);
    fHist_Jet_MCReco_MatchEff_Num_pi  = new TH2F("fHist_Jet_MCReco_MatchEff_Num_pi", "MC Reconstructed Level TOF Matching Efficiency Numerator pion", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins);
    fHist_Jet_MCReco_MatchEff_Num_ka  = new TH2F("fHist_Jet_MCReco_MatchEff_Num_ka", "MC Reconstructed Level TOF Matching Efficiency Numerator kaon", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins);
    fHist_Jet_MCReco_MatchEff_Num_pr  = new TH2F("fHist_Jet_MCReco_MatchEff_Num_pr", "MC Reconstructed Level TOF Matching Efficiency Numerator proton", n_tpc_mom_bins,tpc_mom_bins, n_eta_bins,eta_bins);

    fHistMCT_Jet_ptsub_v_area  = new TH2F("fHistMCT_Jet_ptsub_v_area", "Before cuts, Jet pt after subtraction vs jet area"  , 100, 0., 1., 300, 0., 300.);
    fHistMCR_Jet_ptsub_v_area  = new TH2F("fHistMCR_Jet_ptsub_v_area", "Before cuts, Jet pt after subtraction vs jet area"  , 100, 0., 1., 300, 0., 300.);

    fHistMC_Jet_shared_pt_frac  = new TH2F("fHistMC_Jet_shared_pt_frac", "fHistMC_Jet_shared_pt_frac"  , 300, 0., 300., 222, -10., 1.1);
    fHistMC_Jet_deltaR  = new TH2F("fHistMC_Jet_deltaR", "fHistMC_Jet_deltaR"  , 300, 0., 300., 70, 0., 3.5);

    fHist_JetMatch_MCTruth_TrackEff_Den_pi = new TH3F("fHist_JetMatch_MCTruth_TrackEff_Den_pi", "MC Truth Level Tracking Efficiency Denominator pion for matched jets only", n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins,jet_pT_bins,n_jet_pT_bins,jet_pT_bins);
    fHist_JetMatch_MCTruth_TrackEff_Den_ka = new TH3F("fHist_JetMatch_MCTruth_TrackEff_Den_ka", "MC Truth Level Tracking Efficiency Denominator kaon for matched jets only", n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins,jet_pT_bins,n_jet_pT_bins,jet_pT_bins);
    fHist_JetMatch_MCTruth_TrackEff_Den_pr = new TH3F("fHist_JetMatch_MCTruth_TrackEff_Den_pr", "MC Truth Level Tracking Efficiency Denominator proton for matched jets only", n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins,jet_pT_bins,n_jet_pT_bins,jet_pT_bins);

    fHist_JetMatch_MCReco_TrackEff_Num_pi  = new TH3F("fHist_JetMatch_MCReco_TrackEff_Num_pi", "MC Reconstructed Level Tracking Efficiency Numerator pion for matched jets only",n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins,jet_pT_bins,n_jet_pT_bins,jet_pT_bins);
    fHist_JetMatch_MCReco_TrackEff_Num_ka  = new TH3F("fHist_JetMatch_MCReco_TrackEff_Num_ka", "MC Reconstructed Level Tracking Efficiency Numerator kaon for matched jets only",n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins,jet_pT_bins,n_jet_pT_bins,jet_pT_bins);
    fHist_JetMatch_MCReco_TrackEff_Num_pr  = new TH3F("fHist_JetMatch_MCReco_TrackEff_Num_pr", "MC Reconstructed Level Tracking Efficiency Numerator proton for matched jets only",n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins,jet_pT_bins,n_jet_pT_bins,jet_pT_bins);

    fHist_JetMatch_MCReco_MatchEff_Num_pi  = new TH3F("fHist_JetMatch_MCReco_MatchEff_Num_pi", "MC Reconstructed Level Matching Efficiency Numerator pion for matched jets only", n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins,jet_pT_bins,n_jet_pT_bins,jet_pT_bins);
    fHist_JetMatch_MCReco_MatchEff_Num_ka  = new TH3F("fHist_JetMatch_MCReco_MatchEff_Num_ka", "MC Reconstructed Level Matching Efficiency Numerator kaon for matched jets only", n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins,jet_pT_bins,n_jet_pT_bins,jet_pT_bins);
    fHist_JetMatch_MCReco_MatchEff_Num_pr  = new TH3F("fHist_JetMatch_MCReco_MatchEff_Num_pr", "MC Reconstructed Level Matching Efficiency Numerator proton for matched jets only",n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins,jet_pT_bins,n_jet_pT_bins,jet_pT_bins);
  }

  Float_t true_frac_bins[101] = {};
  for (int i=0; i < (100+1); i++){
    true_frac_bins[i] = 0.0 + i*0.01;
    //size 20 GeV bins from -300 to 0, size 5 GeV bins from 0 to 120, size 20 GeV bins from 120 to 300
  }

  Int_t n_true_frac_bins = sizeof(true_frac_bins)/sizeof(Float_t) - 1;

  Float_t jet_pT_bins2[25] = {};
  for (int i=0; i < (24+1); i++){
    jet_pT_bins2[i] = 0.0 + 5.0*i;
    //size 20 GeV bins from -300 to 0, size 5 GeV bins from 0 to 120, size 20 GeV bins from 120 to 300
  }
  Int_t n_jet_pT_bins2 = sizeof(jet_pT_bins2)/sizeof(Float_t) - 1;


  Float_t deltaR_bins[31] = {};
  for (int i=0; i < (30+1); i++){
    deltaR_bins[i] = 0.0 + i*0.01;
    //size 20 GeV bins from -300 to 0, size 5 GeV bins from 0 to 120, size 20 GeV bins from 120 to 300
  }

  Int_t n_deltaR_bins = sizeof(deltaR_bins)/sizeof(Float_t) - 1;


  Float_t JERS_bins[401] = {};
  for (int i=0; i < (400+1); i++){
    JERS_bins[i] = -10.0 + i*0.05;
    //size 20 GeV bins from -300 to 0, size 5 GeV bins from 0 to 120, size 20 GeV bins from 120 to 300
  }

  Int_t n_JERS_bins = sizeof(JERS_bins)/sizeof(Float_t) - 1;

  if (fDoEmbedding){

    fHistMCT_Jet_ptsub_v_area  = new TH2F("fHistMCT_Jet_ptsub_v_area", "Before cuts, Jet pt after subtraction vs jet area"  , 100, 0., 1., 300, 0., 300.);
    fHistMCR_Jet_ptsub_v_area  = new TH2F("fHistMCR_Jet_ptsub_v_area", "Before cuts, Jet pt after subtraction vs jet area"  , 100, 0., 1., 300, 0., 300.);

    fHist_TrueJetPtFraction = new TH3F("fHist_TrueJetPtFraction", "Jet pt fraction from truth particles", n_jet_pT_bins,jet_pT_bins,n_true_frac_bins,true_frac_bins,n_true_frac_bins,true_frac_bins);

    fHist_MatchJetPts = new TH3F("fHist_MatchJetPts", "Jet pts of matched jets", n_jet_pT_bins,jet_pT_bins,n_jet_pT_bins2,jet_pT_bins2,n_jet_pT_bins2,jet_pT_bins2);
    fHist_MatchJetEtas = new TH3F("fHist_MatchJetEtas", "Jet etas of matched jets", n_eta_bins,eta_bins,n_eta_bins,eta_bins,n_eta_bins,eta_bins);
    fHist_MatchJetDeltaRs = new TH3F("fHist_MatchJetDeltaRs", "Jet distances of matched jets", n_jet_pT_bins,jet_pT_bins,n_deltaR_bins,deltaR_bins,n_deltaR_bins,deltaR_bins);


    fHist_JERS_PbPbunid_tru = new TH3F("fHist_JERS_PbPbunid_tru", "part JERS for unid as fxn of part pT", n_JERS_bins, JERS_bins, n_jet_pT_bins2,jet_pT_bins2,n_tpc_mom_bins,tpc_mom_bins);

    fHist_JERS_PbPbunid_det = new TH3F("fHist_JERS_PbPbunid_det", "det JERS for unid as fxn of det pT", n_JERS_bins, JERS_bins, n_jet_pT_bins2,jet_pT_bins2,n_tpc_mom_bins,tpc_mom_bins);

    fHist_JERS_PbPbunid_truhyb = new TH3F("fHist_JERS_PbPbunid_truhyb", "part JERS for unid as fxn of hyb pT", n_JERS_bins, JERS_bins, n_jet_pT_bins,jet_pT_bins,n_tpc_mom_bins,tpc_mom_bins);

    fHist_JERS_PbPbunid_dethyb = new TH3F("fHist_JERS_PbPbunid_dethyb", "det JERS for unid as fxn of hyb pT", n_JERS_bins, JERS_bins, n_jet_pT_bins,jet_pT_bins,n_tpc_mom_bins,tpc_mom_bins);


    fHist_JERS_truPythia = new TH3F("fHist_JERS_truPythia", "part JERS for Pythia unid as fxn of part pT", n_JERS_bins, JERS_bins, n_jet_pT_bins2,jet_pT_bins2,n_tpc_mom_bins,tpc_mom_bins);

    fHist_JERS_truPythia_pi = new TH3F("fHist_JERS_truPythia_pi", "part JERS for Pythia pi as fxn of part pT", n_JERS_bins, JERS_bins, n_jet_pT_bins2,jet_pT_bins2,n_tpc_mom_bins,tpc_mom_bins);
    fHist_JERS_truPythia_ka = new TH3F("fHist_JERS_truPythia_ka", "part JERS for Pythia ka as fxn of part pT", n_JERS_bins, JERS_bins, n_jet_pT_bins2,jet_pT_bins2,n_tpc_mom_bins,tpc_mom_bins);
    fHist_JERS_truPythia_pr = new TH3F("fHist_JERS_truPythia_pr", "part JERS for Pythia pr as fxn of part pT", n_JERS_bins, JERS_bins, n_jet_pT_bins2,jet_pT_bins2,n_tpc_mom_bins,tpc_mom_bins);


    fHist_JERS_detPythia = new TH3F("fHist_JERS_detPythia", "det JERS for Pythia unid as fxn of det pT", n_JERS_bins, JERS_bins, n_jet_pT_bins2,jet_pT_bins2,n_tpc_mom_bins,tpc_mom_bins);

    fHist_JERS_detPythia_pi = new TH3F("fHist_JERS_detPythia_pi", "det JERS for Pythia pi as fxn of det pT", n_JERS_bins, JERS_bins, n_jet_pT_bins2,jet_pT_bins2,n_tpc_mom_bins,tpc_mom_bins);
    fHist_JERS_detPythia_ka = new TH3F("fHist_JERS_detPythia_ka", "det JERS for Pythia ka as fxn of det pT", n_JERS_bins, JERS_bins, n_jet_pT_bins2,jet_pT_bins2,n_tpc_mom_bins,tpc_mom_bins);
    fHist_JERS_detPythia_pr = new TH3F("fHist_JERS_detPythia_pr", "det JERS for Pythia pr as fxn of det pT", n_JERS_bins, JERS_bins, n_jet_pT_bins2,jet_pT_bins2,n_tpc_mom_bins,tpc_mom_bins);


    fHist_JERS_truhybPythia = new TH3F("fHist_JERS_truhybPythia", "part JERS for Pythia unid as fxn of hyb pT", n_JERS_bins, JERS_bins, n_jet_pT_bins,jet_pT_bins,n_tpc_mom_bins,tpc_mom_bins);

    fHist_JERS_truhybPythia_pi = new TH3F("fHist_JERS_truhybPythia_pi", "part JERS for Pythia pi as fxn of hyb pT", n_JERS_bins, JERS_bins, n_jet_pT_bins,jet_pT_bins,n_tpc_mom_bins,tpc_mom_bins);
    fHist_JERS_truhybPythia_ka = new TH3F("fHist_JERS_truhybPythia_ka", "part JERS for Pythia ka as fxn of hyb pT", n_JERS_bins, JERS_bins, n_jet_pT_bins,jet_pT_bins,n_tpc_mom_bins,tpc_mom_bins);
    fHist_JERS_truhybPythia_pr = new TH3F("fHist_JERS_truhybPythia_pr", "part JERS for Pythia pr as fxn of hyb pT", n_JERS_bins, JERS_bins, n_jet_pT_bins,jet_pT_bins,n_tpc_mom_bins,tpc_mom_bins);


    fHist_JERS_dethybPythia = new TH3F("fHist_JERS_dethybPythia", "det JERS for Pythia unid as fxn of hyb pT", n_JERS_bins, JERS_bins, n_jet_pT_bins,jet_pT_bins,n_tpc_mom_bins,tpc_mom_bins);

    fHist_JERS_dethybPythia_pi = new TH3F("fHist_JERS_dethybPythia_pi", "det JERS for Pythia pi as fxn of hyb pT", n_JERS_bins, JERS_bins, n_jet_pT_bins,jet_pT_bins,n_tpc_mom_bins,tpc_mom_bins);
    fHist_JERS_dethybPythia_ka = new TH3F("fHist_JERS_dethybPythia_ka", "det JERS for Pythia ka as fxn of hyb pT", n_JERS_bins, JERS_bins, n_jet_pT_bins,jet_pT_bins,n_tpc_mom_bins,tpc_mom_bins);
    fHist_JERS_dethybPythia_pr = new TH3F("fHist_JERS_dethybPythia_pr", "det JERS for Pythia pr as fxn of hyb pT", n_JERS_bins, JERS_bins, n_jet_pT_bins,jet_pT_bins,n_tpc_mom_bins,tpc_mom_bins);

    fHist_PC_spectra_Pythia = new TH1F("fHist_PC_spectra_Pythia", "Pythia parts in PC unid spectra", n_tpc_mom_bins,tpc_mom_bins);
    fHist_PC_spectra_Pythia_pi = new TH1F("fHist_PC_spectra_Pythia_pi", "Pythia parts in PC pi spectra", n_tpc_mom_bins,tpc_mom_bins);
    fHist_PC_spectra_Pythia_ka = new TH1F("fHist_PC_spectra_Pythia_ka", "Pythia parts in PC ka spectra", n_tpc_mom_bins,tpc_mom_bins);
    fHist_PC_spectra_Pythia_pr = new TH1F("fHist_PC_spectra_Pythia_pr", "Pythia parts in PC pr spectra", n_tpc_mom_bins,tpc_mom_bins);

    fHist_det_spectra_Pythia = new TH2F("fHist_det_spectra_Pythia", "Pythia parts in det unid spectra", n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins2,jet_pT_bins2);
    fHist_det_spectra_Pythia_pi = new TH2F("fHist_det_spectra_Pythia_pi", "Pythia parts in det pi spectra", n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins2,jet_pT_bins2);
    fHist_det_spectra_Pythia_ka = new TH2F("fHist_det_spectra_Pythia_ka", "Pythia parts in det ka spectra", n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins2,jet_pT_bins2);
    fHist_det_spectra_Pythia_pr = new TH2F("fHist_det_spectra_Pythia_pr", "Pythia parts in det pr spectra", n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins2,jet_pT_bins2);

    fHist_part_spectra_Pythia = new TH2F("fHist_part_spectra_Pythia", "Pythia parts in part unid spectra", n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins2,jet_pT_bins2);
    fHist_part_spectra_Pythia_pi = new TH2F("fHist_part_spectra_Pythia_pi", "Pythia parts in part pi spectra", n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins2,jet_pT_bins2);
    fHist_part_spectra_Pythia_ka = new TH2F("fHist_part_spectra_Pythia_ka", "Pythia parts in part ka spectra", n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins2,jet_pT_bins2);
    fHist_part_spectra_Pythia_pr = new TH2F("fHist_part_spectra_Pythia_pr", "Pythia parts in part pr spectra", n_tpc_mom_bins,tpc_mom_bins,n_jet_pT_bins2,jet_pT_bins2);

  }

  if (fDoEventHistos){
    fListHist->Add(fHistCentrality);
    fListHist->Add(fHist_Is_Good_Inc_Event);
    fListHist->Add(fHist_Has_Acc_Real_FJ_Jet);
    fListHist->Add(fHist_Has_Acc_Real_EMC_Jet);
    fListHist->Add(fHist_Num_Real_FJ_Jets);
    fListHist->Add(fHist_Num_Real_EMC_Jets);
    fListHist->Add(fHist_RhoFJ);
    fListHist->Add(fHist_RhoEMC);
    fListHist->Add(fHist_Filled_UEC_Rec);
    if (fMCtrue || fDoEmbedding){
      fListHist->Add(fHist_RhoEMCReco);
      fListHist->Add(fHist_RhoEMCGen);
      fListHist->Add(fHist_RhoEMCEmb);
      fListHist->Add(fHist_Filled_UEC_Emb);
      fListHist->Add(fHist_Filled_UEC_Gen);
    }
  }


  if (!fMCtrue || fDoEmbedding){

    if (fFill_TPC) {
      fListHist->Add(fHistIncTracks_dEdx);
      fListHist->Add(fHistJetTracks_dEdx);
    }
    if (fFillpTPC_pT){
      fListHist->Add(fHistIncTracks_moms);
      fListHist->Add(fHistJetTracks_moms);
    }

    if (fFillp_pT){
      fListHist->Add(fHistIncTracks_moms_p);
      fListHist->Add(fHistJetTracks_moms_p);
    }

    if (fFillpTPC_p){
      fListHist->Add(fHistIncTracks_moms_pTPC_p);
      fListHist->Add(fHistJetTracks_moms_pTPC_p);
    }

      fListHist->Add(fHistIncTracks_kin);

    if (fFill_TOF){
      fListHist->Add(fHistIncTracks_beta);
      fListHist->Add(fHistIncTracks_TOFpi_nsigma);
      fListHist->Add(fHistIncTracks_TOFka_nsigma);
      fListHist->Add(fHistIncTracks_TOFpr_nsigma);

      fListHist->Add(fHistJetTracks_beta);
      fListHist->Add(fHistJetTracks_TOFpi_nsigma);
      fListHist->Add(fHistJetTracks_TOFka_nsigma);
      fListHist->Add(fHistJetTracks_TOFpr_nsigma);

      fListHist->Add(fHistIncTracks_t0);
    }

    fListHist->Add(fHistJetTracks_kin);

    if (fFill_TOF_expecs){
      fListHist->Add(fHistBetaExpec_pi);
      fListHist->Add(fHistBetaExpec_ka);
      fListHist->Add(fHistBetaExpec_pr);
      fListHist->Add(fHistjet_BetaExpec_pi);
      fListHist->Add(fHistjet_BetaExpec_ka);
      fListHist->Add(fHistjet_BetaExpec_pr);
      fListHist->Add(fHist_pi_mismatch);
      fListHist->Add(fHist_ka_mismatch);
      fListHist->Add(fHist_pr_mismatch);
      fListHist->Add(fHist_jet_pi_mismatch);
      fListHist->Add(fHist_jet_ka_mismatch);
      fListHist->Add(fHist_jet_pr_mismatch);
      fListHist->Add(fHist_elExpec_pihyp);
      fListHist->Add(fHist_muExpec_pihyp);
      fListHist->Add(fHist_kaExpec_pihyp);
      fListHist->Add(fHist_prExpec_pihyp);
      fListHist->Add(fHist_piExpec_kahyp);
      fListHist->Add(fHist_prExpec_kahyp);
      fListHist->Add(fHist_piExpec_prhyp);
      fListHist->Add(fHist_kaExpec_prhyp);
      fListHist->Add(fHist_deExpec_prhyp);

      fListHist->Add(fHist_jet_elExpec_pihyp);
      fListHist->Add(fHist_jet_muExpec_pihyp);
      fListHist->Add(fHist_jet_kaExpec_pihyp);
      fListHist->Add(fHist_jet_prExpec_pihyp);
      fListHist->Add(fHist_jet_piExpec_kahyp);
      fListHist->Add(fHist_jet_prExpec_kahyp);
      fListHist->Add(fHist_jet_piExpec_prhyp);
      fListHist->Add(fHist_jet_kaExpec_prhyp);
      fListHist->Add(fHist_jet_deExpec_prhyp);

      fListHist->Add(fHistTOFSigmaExpec_pi);
      fListHist->Add(fHistTOFSigmaExpec_ka);
      fListHist->Add(fHistTOFSigmaExpec_pr);
      fListHist->Add(fHistjet_TOFSigmaExpec_pi);
      fListHist->Add(fHistjet_TOFSigmaExpec_ka);
      fListHist->Add(fHistjet_TOFSigmaExpec_pr);
    }

    fListHist->Add(fHistJet_ptsub_v_area);
    fListHist->Add(fHistJet_kin);
    fListHist->Add(fHistJet_moms);

    fListHist->Add(fHist_pi_DCAxy);
    fListHist->Add(fHist_pr_DCAxy);

    fListHist->Add(fHistJet_pi_DCAxy);
    fListHist->Add(fHistJet_pr_DCAxy);

    if (fFill_TPC_expecs){
      fListHist->Add(fHistIncTracks_mpi);
      fListHist->Add(fHistIncTracks_spi);

      fListHist->Add(fHistIncTracks_mel);
      fListHist->Add(fHistIncTracks_sel);

      fListHist->Add(fHistIncTracks_mka);
      fListHist->Add(fHistIncTracks_ska);

      fListHist->Add(fHistIncTracks_mpr);
      fListHist->Add(fHistIncTracks_spr);
    }
  }

  if (fMCtrue && !fDoEmbedding){
    fListHist->Add(fHistMCTruth_TrackEff_Den_pi);
    fListHist->Add(fHistMCTruth_TrackEff_Den_ka);
    fListHist->Add(fHistMCTruth_TrackEff_Den_pr);
    fListHist->Add(fHistMCReco_TrackEff_Num_pi);
    fListHist->Add(fHistMCReco_TrackEff_Num_ka);
    fListHist->Add(fHistMCReco_TrackEff_Num_pr);
    fListHist->Add(fHistMCReco_MatchEff_Num_pi);
    fListHist->Add(fHistMCReco_MatchEff_Num_ka);
    fListHist->Add(fHistMCReco_MatchEff_Num_pr);

    fListHist->Add(fHistMCReco_pi_prim_DCAxy);
    fListHist->Add(fHistMCReco_pi_scdweak_DCAxy);
    fListHist->Add(fHistMCReco_pr_prim_DCAxy);
    fListHist->Add(fHistMCReco_pr_scdweak_DCAxy);
    fListHist->Add(fHistMCReco_pr_scdmat_DCAxy);

    fListHist->Add(fHist_Jet_MCTruth_TrackEff_Den_pi);
    fListHist->Add(fHist_Jet_MCTruth_TrackEff_Den_ka);
    fListHist->Add(fHist_Jet_MCTruth_TrackEff_Den_pr);
    fListHist->Add(fHist_Jet_MCReco_TrackEff_Num_pi);
    fListHist->Add(fHist_Jet_MCReco_TrackEff_Num_ka);
    fListHist->Add(fHist_Jet_MCReco_TrackEff_Num_pr);
    fListHist->Add(fHist_Jet_MCReco_MatchEff_Num_pi);
    fListHist->Add(fHist_Jet_MCReco_MatchEff_Num_ka);
    fListHist->Add(fHist_Jet_MCReco_MatchEff_Num_pr);

    fListHist->Add(fHist_JetMatch_MCTruth_TrackEff_Den_pi);
    fListHist->Add(fHist_JetMatch_MCTruth_TrackEff_Den_ka);
    fListHist->Add(fHist_JetMatch_MCTruth_TrackEff_Den_pr);
    fListHist->Add(fHist_JetMatch_MCReco_TrackEff_Num_pi);
    fListHist->Add(fHist_JetMatch_MCReco_TrackEff_Num_ka);
    fListHist->Add(fHist_JetMatch_MCReco_TrackEff_Num_pr);
    fListHist->Add(fHist_JetMatch_MCReco_MatchEff_Num_pi);
    fListHist->Add(fHist_JetMatch_MCReco_MatchEff_Num_ka);
    fListHist->Add(fHist_JetMatch_MCReco_MatchEff_Num_pr);

    fListHist->Add(fHistMC_Jet_shared_pt_frac);
    fListHist->Add(fHistMC_Jet_deltaR);

    fListHist->Add(fHistMCT_Jet_ptsub_v_area);
    fListHist->Add(fHistMCR_Jet_ptsub_v_area);
  }

  if (fDoEmbedding){
    
    fListHist->Add(fHistMCT_Jet_ptsub_v_area);
    fListHist->Add(fHistMCR_Jet_ptsub_v_area);

    fListHist->Add(fHist_TrueJetPtFraction);

    fListHist->Add(fHist_MatchJetPts);
    fListHist->Add(fHist_MatchJetEtas);
    fListHist->Add(fHist_MatchJetDeltaRs);

    fListHist->Add(fHist_JERS_PbPbunid_tru);
    fListHist->Add(fHist_JERS_PbPbunid_det);
    fListHist->Add(fHist_JERS_PbPbunid_truhyb);
    fListHist->Add(fHist_JERS_PbPbunid_dethyb);

    fListHist->Add(fHist_JERS_truPythia);
    fListHist->Add(fHist_JERS_truPythia_pi);
    fListHist->Add(fHist_JERS_truPythia_ka);
    fListHist->Add(fHist_JERS_truPythia_pr);

    fListHist->Add(fHist_JERS_detPythia);
    fListHist->Add(fHist_JERS_detPythia_pi);
    fListHist->Add(fHist_JERS_detPythia_ka);
    fListHist->Add(fHist_JERS_detPythia_pr);

    fListHist->Add(fHist_JERS_truhybPythia);
    fListHist->Add(fHist_JERS_truhybPythia_pi);
    fListHist->Add(fHist_JERS_truhybPythia_ka);
    fListHist->Add(fHist_JERS_truhybPythia_pr);

    fListHist->Add(fHist_JERS_dethybPythia);
    fListHist->Add(fHist_JERS_dethybPythia_pi);
    fListHist->Add(fHist_JERS_dethybPythia_ka);
    fListHist->Add(fHist_JERS_dethybPythia_pr);

    fListHist->Add(fHist_PC_spectra_Pythia);
    fListHist->Add(fHist_PC_spectra_Pythia_pi);
    fListHist->Add(fHist_PC_spectra_Pythia_ka);
    fListHist->Add(fHist_PC_spectra_Pythia_pr);

    fListHist->Add(fHist_det_spectra_Pythia);
    fListHist->Add(fHist_det_spectra_Pythia_pi);
    fListHist->Add(fHist_det_spectra_Pythia_ka);
    fListHist->Add(fHist_det_spectra_Pythia_pr);

    fListHist->Add(fHist_part_spectra_Pythia);
    fListHist->Add(fHist_part_spectra_Pythia_pi);
    fListHist->Add(fHist_part_spectra_Pythia_ka);
    fListHist->Add(fHist_part_spectra_Pythia_pr);

  }

  //
  // ************************************************************************
  //   Trees
  // ************************************************************************
  //
  fTreejetsEMCconst  = ((*fTreeSRedirector)<<"jetsEMCconst").GetTree();
  fTreejetsEMCBGconst  = ((*fTreeSRedirector)<<"jetsEMCBGconst").GetTree();
  fTreejetsFJ        = ((*fTreeSRedirector)<<"jetsFJ").GetTree();
  fTreejetsFJBG      = ((*fTreeSRedirector)<<"jetsFJBG").GetTree();
  fTreejetsEMC       = ((*fTreeSRedirector)<<"jetsEMC").GetTree();
  fTreejetsEMCBG       = ((*fTreeSRedirector)<<"jetsEMCBG").GetTree();
  fTreejetsFJconst   = ((*fTreeSRedirector)<<"jetsFJconst").GetTree();
  fTreejetsFJBGconst   = ((*fTreeSRedirector)<<"jetsFJBGconst").GetTree();
  fTreejetEvents     = ((*fTreeSRedirector)<<"jeteventInfo").GetTree();
  fTreeMC        = ((*fTreeSRedirector)<<"fTreeMC").GetTree();
  fTreeCuts      = ((*fTreeSRedirector)<<"tracks").GetTree();
  //
  // ************************************************************************
  //   Send output objects to container
  // ************************************************************************
  //
  PostData(1, fListHist);

  if (fDoEventTree){
    PostData(2, fTreejetsEMCconst);
    PostData(3, fTreejetsEMCBGconst);
    PostData(4, fTreejetsFJ);
    PostData(5, fTreejetsFJBG);
    PostData(6, fTreejetsFJconst);
    PostData(7, fTreejetsFJBGconst);
    PostData(8, fTreejetEvents);
    PostData(9, fTreejetsEMC);
    PostData(10, fTreejetsEMCBG);
    PostData(11, fTreeMC);
    PostData(12, fTreeCuts);
  } 

  std::cout << " Info::siweyhmi: ===== Out of UserCreateOutputObjects ===== " << std::endl;

}
//________________________________________________________________________
Bool_t AliAnalysisTaskJetHadroAOD::Run()
{
  //
  // main event loop
  //
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the UserExec ===== " << std::endl;
  //
  // Check Monte Carlo information and other access first:
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) {
    fMCtrue = kFALSE;
    if (fUseCouts) std::cout << " Info::siweyhmi: ===== There is no MCEventHandler ===== " << std::endl;
  }
  fEventCountInFile++;
  //
  //  get the filename
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { Printf(" Error::siweyhmi: Could not receive input chain"); return kFALSE;}
  TObjString fileName(chain->GetCurrentFile()->GetName());
  //
  // ======================================================================
  // ========================== See if MC or Real =========================
  // ======================================================================
  //
  if (eventHandler) fMCEvent = MCEvent(); //MCEvent(); //eventHandler->MCEvent();
  if (!fMCEvent) if (fUseCouts) std::cout << " Info::siweyhmi: ===== There is no MCEvent ===== " << std::endl;

  //
  // If AODs exist get some event variables
  //
  fCentrality = -5;
  fCentImpBin =-10.;

  fAOD = dynamic_cast<AliAODEvent*>( InputEvent() );
  if (fAOD)
  {
    //
    // Init magnetic filed for golden chi2 cut
    fAOD->InitMagneticField();
    //
    // event selection
      if ( (fPassIndex==3 || fPassIndex==2) && fYear>2013){
        //
        fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,1); // standard
        if (!fEventCuts.AcceptEvent(fAOD)) {cout<< "pileup event " << endl; return kFALSE;}
      }
    //
    //

    if (fYear!=2017){
      AliCentrality    *Centrality = 0x0;
      AliMultSelection *MultSelection = 0x0;
      AliMultSelectionTask *MultSelectionTask = 0x0;
      Centrality = fAOD->GetCentrality();
      MultSelection = (AliMultSelection*) fAOD-> FindListObject("MultSelection");

      if (MultSelection) {
        if (fMultTrueFlag) fCentrality = MultSelection->GetMultiplicityPercentile("V0M",kTRUE);
        else fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
      } else if (Centrality) {
        fCentrality = Centrality->GetCentralityPercentile("V0M");
      } else {
        std::cout << " Info::siweyhmi: Error: There is no cent info " << std::endl;
      }

    }

  }

  //
  // Get rid of "E-AliESDpid::GetTPCsignalTunedOnData: Tune On Data requested, but MC event not set. Call SetCurrentMCEvent before!" errors
  //if (!fPIDResponse) fPIDResponse = ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetESDpid();
  if (fMCEvent && fPIDResponse) fPIDResponse->SetCurrentMCEvent(fMCEvent);
  //
  if(fMCtrue)
  {
    //
    // Add different generator particles to PDG Data Base to avoid problems when reading MC generator particles
    AliPDG::AddParticlesToPdgDataBase();
  }
  
  //
  // ========================== Real =========================
  //
  if (!fAOD)          { Printf(" Error::siweyhmi: fAOD not available"); return kFALSE; }
  //
  // ------------------------------------------------
  // ------- monitor vertex position =---------------
  // ------------------------------------------------
  //
  Bool_t isVertexOk = kTRUE;
  fVertex_AOD = fAOD->GetPrimaryVertexTracks();
  if( fVertex_AOD->GetNContributors()<1) isVertexOk = kFALSE;
  if( fVertex_AOD->GetNContributors()>1) {
    fVz    = fVertex_AOD->GetZ();
  }
  fNContributors   = fVertex_AOD->GetNContributors();
  //
  // ------------------------------------------------
  // ------- event vertex cut along Z ---------------
  // ------------------------------------------------
  //

  if (TMath::Abs(fVz)>10) return kFALSE; //Data/ Reco level
  //
  if (fVertex_AOD && isVertexOk) fHistVertex->Fill(fVz);
  else return kFALSE;
  //
  // ------------------------------------------------
  // ---------- Centrality definition ---------------
  // ------------------------------------------------
  //


  //
  if (fUseCouts) {
    std::cout << " Info::siweyhmi: =============================================================================================== " << std::endl;
    std::cout << " Info::siweyhmi: Event counter = " << fEventCountInFile << " - cent =  " << fCentrality << std::endl;
    std::cout << " Info::siweyhmi: =============================================================================================== " << std::endl;
  }
  fHistCentrality->Fill(fCentrality);  // count events after physics and vertex selection
  //
  // in case small stat is enough
  if (fPercentageOfEvents>0 && (fEventCountInFile%fPercentageOfEvents)!=0) return kFALSE;

  // if centrality estimation failed
  if (fCentrality > 100.) return kFALSE;
  //
  // ======================================================================
  //   Filling part
  // ======================================================================
  //
  // Real Data Analysis
  //
  if (!fMCtrue && fAOD && fCentrality>=fcent_min && fCentrality<fcent_max){
    fisGoodIncEvent = 0;
    fFilledUECone_Emb = 0;
    fFilledUECone_Gen = 0;
    fFilledUECone_Rec = 0;
    fhasAcceptedFJjet = 0;
    fhasRealFJjet = 0;
    fhasAcceptedEMCjet = 0;
    fhasRealEMCjet = 0;
    frhoFJ = 0.0;
    fjetRhoVal = 0.0;
    fNumRealFJJets = 0;
    fNumRealEMCJets = 0;

    if (fDoIncTracks){
      FillIncTracksReal();
    }
    if (fDoEMCJet){
      FindJetsEMC();
    }
    if (fDoFastJet){
      FindJetsFJ();
    }
    if (fDoEventTree){
      FillEventTree();
    }
    if (fDoEventHistos){
      FillEventHistos();
    }
    if (fUseCouts)  std::cout << " Info::siweyhmi: (Real Data Analysis) End of Filling part = " << fEventCountInFile << std::endl;
    return kTRUE;
  }
  //
  // full MC analysis
  //

  if (fMCtrue && fAOD  && fCentrality>=fcent_min && fCentrality<fcent_max){
    fisGoodIncEvent = 0;
    fFilledUECone_Emb = 0;
    fFilledUECone_Gen = 0;
    fFilledUECone_Rec = 0;
    fhasAcceptedFJjet = 0;
    fhasRealFJjet = 0;
    fhasAcceptedEMCjet = 0;
    fhasRealEMCjet = 0;
    frhoFJ = 0.0;
    fjetRhoVal = 0.0;
    fjetGenRhoVal = 0.0;
    fjetRecoRhoVal = 0.0;
    fjetEmbRhoVal = 0.0;
    fNumRealFJJets = 0;
    fNumRealEMCJets = 0;
    fall_reco_jets_w_multiple_matches = 0;
    fall_reco_jets_w_matches = 0;

    if (fDoIncTracks){
      FillTreeMC();
    }
    if (fDoFastJet){
      FillMCJets();
    }

    if (fDoEmbedding){
      FillEmbJets();
    }
    if (fDoEventTree){
      FillEventTree();
    }
    if (fDoEventHistos){
      FillEventHistos();
    }
    if (fUseCouts)  std::cout << " Info::siweyhmi: (full MC analysis) End of Filling part = " << fEventCountInFile << std::endl;

    if (fUseCouts) cout << "num_reco_jets_w_multiple_matches over all events is " <<  fall_reco_jets_w_multiple_matches << endl;
    if (fUseCouts) cout << "num_reco_jets_w_matches over all events is " <<  fall_reco_jets_w_matches << endl;

    return kTRUE;
  }
  return kTRUE;

}
//________________________________________________________________________
void AliAnalysisTaskJetHadroAOD::FindJetsEMC()
{

  //
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FindJetsEMC ===== " << std::endl;
  AliTOFPIDResponse fTOFPIDResponse = fPIDResponse->GetTOFResponse();

  //
  Float_t pT_sub_min = fjetMinPtSub;
  if (fUseCouts) cout << "Minimum jet pt after subtraction is " << fjetMinPtSub << endl;
  if (fUseCouts) cout << "Maximum jet pt after subtraction is " << fjetMaxPtSub << endl;
  // Get the jet container
  fJetContainer = this->GetJetContainer("detJets");
  if (fUseCouts) cout << "fYear is " << fYear << endl;
  if (fYear!=2017){
    TString fRhoName = fJetContainer->GetRhoName();
    if (fUseCouts) cout << "Rho Name is " << fRhoName << endl;

    if (fJetContainer->GetRhoParameter()) fjetRhoVal = fJetContainer->GetRhoVal();
  }
  else {
    fjetRhoVal=0.0;
  }
  if (fUseCouts) cout << "In FindJetsEMC Rho value is " << fjetRhoVal << endl;

  //
  // Get some jet container properties
  Int_t Njets         = fJetContainer->GetNJets();
  Int_t NAcceptedjets = fJetContainer->GetNAcceptedJets();
  Float_t jetRadius  = fJetContainer->GetJetRadius();
  Float_t jetEtaMin = fJetContainer->GetJetEtaMin();
  Float_t jetEtaMax = fJetContainer->GetJetEtaMax();
  UInt_t jetAcceptanceType = fJetContainer->GetAcceptanceType();
  //
  // loop over jets
  Float_t leadJetPhi = 999;
  Float_t leadJetEta = 999;
  Float_t leadJetPtSub = -999;

  for(auto jet : fJetContainer->accepted())
  {
    //if (TMath::Abs(jet->Eta()) >= 0.9-jetRadius) continue; //this is not needed when using jetcontainer->accepted
    fhasAcceptedEMCjet = 1;
    fJetPt = jet->Pt();
    fJetEta = jet->Eta();
    fJetPhi = jet->Phi();
    Float_t JetM = jet->M();
    Int_t JetLabel = jet->GetLabel();
    Float_t JetAreaPt = jet->AreaPt(); //jet transverse area
    Int_t JetNumberOfConstituents = jet->GetNumberOfParticleConstituents(); //tracks + clusters = particles

    Float_t particleMaxChargedPt = jet->MaxChargedPt(); //max pT of charged particle in jet
    Float_t jetptsub = jet->PtSub(fjetRhoVal, kFALSE);

    if (jetptsub > pT_sub_min && jetptsub < fjetMaxPtSub && JetAreaPt > fjetMinArea)
    {
      fhasRealEMCjet = 1;
      fNumRealEMCJets +=1;
    }

    if(!fTreeSRedirector) return;
    if (fFillEMCJet){
      (*fTreeSRedirector)<<"jetsEMC"<<
      "NAcceptedjets=" << NAcceptedjets << // Number of accepted jets in event
      "jetRhoVal=" << fjetRhoVal <<
      "jetpt="     << fJetPt    << // jetPt
      "jeteta="    << fJetEta   << // jetEta
      "jetphi="    << fJetPhi   << // jetPhi
      "jetArea="   << JetAreaPt  <<
      "jetNumberOfConstituents="    << JetNumberOfConstituents  <<
      "particleMaxChargedPt="  << particleMaxChargedPt <<
      "jetptsub="  << jetptsub <<
      "cent="      << fCentrality  <<  //  centrality
      "\n";
    }

    fHistJet_ptsub_v_area->Fill(JetAreaPt, jetptsub);

    if (JetAreaPt <= fjetMinArea) continue;

    fHistJet_kin->Fill(jetptsub, fJetEta, fJetPhi);
    fHistJet_moms->Fill(jetptsub, fJetPt);

    if (jetptsub <= pT_sub_min || jetptsub >= fjetMaxPtSub) continue;

    if (jetptsub > leadJetPtSub){
      leadJetPtSub = jetptsub;
      leadJetEta = fJetEta;
      leadJetPhi = fJetPhi;
    }


    for(Int_t i = 0; i < jet->GetNumberOfParticleConstituents(); i++)
    {
      const AliVParticle* particle = jet->Track(i);
      AliAODTrack* AODtrack = (AliAODTrack*)(particle);
      if (!AODtrack) continue;

      //Track cuts start
      fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
      fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
      //
      //if (!fAODtrackCuts->AcceptTrack(AODtrack))  continue;    // default cuts - redundant since these track cuts are passed to jet finder
      //
      // Get the track variables
      GetExpecteds(AODtrack);
      SetCutBitsAndSomeTrackVariables(AODtrack);
      Float_t length    = AODtrack->GetIntegratedLength();
      Float_t tofSignal = AODtrack->GetTOFsignal();
      Float_t t0 = fTOFPIDResponse.GetStartTime(fPVertex);
      Float_t beta = -.05;
      //int nTOFClusters = AODtrack->GetNTOFclusters(); //All matchable clusters AOD doesn't have this
      if((length > 0) && (tofSignal > 0)) beta = length / (2.99792458e-2*(tofSignal - t0));

      Float_t fTPCmom_choice = fPtot;
      if (fSetTPCmom == 1) fTPCmom_choice = fPVertex;
      if (fSetTPCmom == 2) fTPCmom_choice = fPt;

      Float_t fTOFmom_choice = fPt;
      if (fSetTOFmom == 1) fTPCmom_choice = fPVertex;
      if (fSetTOFmom == 2) fTPCmom_choice = fPtot;

      Float_t fBetamom_choice = fPVertex;
      if (fSetBetamom == 1) fTPCmom_choice = fPt;
      if (fSetBetamom == 2) fTPCmom_choice = fPtot;

      Float_t fEta_choice = fEta;
      if (fSetEta == 1) fEta_choice = fY;

      Float_t DCAxy = -999.;
      Float_t DCAz = -999.;
      AODtrack->GetImpactParameters(DCAxy, DCAz);
      Float_t comb_sig_pi = TMath::Sqrt(fNSigmasPiTPC * fNSigmasPiTPC + fNSigmasPiTOF * fNSigmasPiTOF);
      Float_t comb_sig_pr = TMath::Sqrt(fNSigmasPrTPC * fNSigmasPrTPC + fNSigmasPrTOF * fNSigmasPrTOF);

      if (TMath::Abs(fEta) < fEtaCut){
        if (fFill_TPC) fHistJetTracks_dEdx->Fill(fTPCmom_choice,fTPCSignal,TMath::Abs(fEta_choice));
        if (fFillpTPC_pT) fHistJetTracks_moms->Fill(fPt,fPtot);
        if (fFillp_pT) fHistJetTracks_moms_p->Fill(fPt,fPVertex);
        if (fFillpTPC_p) fHistJetTracks_moms_pTPC_p->Fill(fPtot,fPVertex);
        fHistJetTracks_kin->Fill(fPt,TMath::Abs(fEta),fPhi);
        if (comb_sig_pi < 2.0) fHistJet_pi_DCAxy->Fill(fPt, DCAxy);
        if (comb_sig_pr < 2.0) fHistJet_pr_DCAxy->Fill(fPt, DCAxy);
      }

      //TPC-TOF Matching conditions: Use standard TPC tracks, then require kTIME and kTOFout
      Bool_t fTOFout = kFALSE;
      if ((AODtrack->GetStatus() & AliAODTrack::kTOFout) != 0) { //Track has the kTOFout flag
      fTOFout = kTRUE;
      }

      //
      //kTIME flag
      Bool_t fTime = kFALSE;
      if ((AODtrack->GetStatus() & AliAODTrack::kTIME) != 0) { //Track has the kTIME flag
      fTime = kTRUE;
      }

      Double_t inttime[5]; //6 is needed to account for earlier species - 0 = electron, 1 = muon, 2 = pion, 3 = kaon, 4 = proton, ESD only: 5 = deuteron
      AODtrack->GetIntegratedTimes(inttime, 5); // Returns the array with integrated times for each particle hypothesis
      Float_t fTOFMismatchTime = -20000000.;
      fTOFMismatchTime = AliTOFPIDResponse::GetMismatchRandomValue(AODtrack->Eta());
      if (fTOFMismatchTime <= 0){
        fTOFMismatchTime = -20000000.;
      }

      for (Int_t i = 0; i < 3; i++) {
        const Float_t beta_expec = length / (2.99792458e-2 * (inttime[i+2] - t0));
        if (i==0){
          fHistjet_BetaExpec_pi->Fill(fTPCmom_choice, beta_expec);
          fHist_jet_pi_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          if (TMath::Abs(fNSigmasPiTOF) < 2.0){
            fHistjet_TOFSigmaExpec_pi->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasElTOF) < 2.0){
            fHist_jet_elExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasMuTOF) < 2.0){
            fHist_jet_muExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasKaTOF) < 2.0){
            fHist_jet_kaExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasPrTOF) < 2.0){
            fHist_jet_prExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
        }
        if (i==1){
          fHistjet_BetaExpec_ka->Fill(fTPCmom_choice, beta_expec);
          fHist_jet_ka_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
          if (TMath::Abs(fNSigmasKaTOF) < 2.0){
            fHistjet_TOFSigmaExpec_ka->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasPiTOF) < 2.0){
            fHist_jet_piExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasPrTOF) < 2.0){
            fHist_jet_prExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
          }
        }
        if (i==2){
          fHistjet_BetaExpec_pr->Fill(fTPCmom_choice, beta_expec);
          fHist_jet_pr_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
          if (TMath::Abs(fNSigmasPrTOF) < 2.0){
            fHistjet_TOFSigmaExpec_pr->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasPiTOF) < 2.0){
            fHist_jet_piExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasKaTOF) < 2.0){
            fHist_jet_kaExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasDeTOF) < 2.0){
            fHist_jet_deExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
          }
        }
      }

      if (TMath::Abs(fEta) < fEtaCut  && fTOFout && fTime && fFill_TOF){
        fHistJetTracks_beta->Fill(fBetamom_choice,beta);
        fHistJetTracks_TOFpi_nsigma->Fill(fTOFmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
        fHistJetTracks_TOFka_nsigma->Fill(fTOFmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
        fHistJetTracks_TOFpr_nsigma->Fill(fTOFmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
        /*if (nTOFClusters < 2) { //ESD only
          fHistJetTracks_TOFpi_nsigma_1cls->Fill(fTOFmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
          fHistJetTracks_TOFka_nsigma_1cls->Fill(fTOFmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
          fHistJetTracks_TOFpr_nsigma_1cls->Fill(fTOFmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
        }*/
      }

      if (fFillJetsEMCConst && fTOFout && fTime){

        (*fTreeSRedirector)<<"jetsEMCconst"<<
        //
        "jetRadius=" << jetRadius << // jet Radius
        "Njets="     << Njets << // Number of jets in event
        "NAcceptedjets=" << NAcceptedjets << // Number of accepted jets in event
        "jetRhoVal=" << fjetRhoVal <<
        "jetpt="     << fJetPt    << // jetPt
        "jeteta="    << fJetEta   << // jetEta
        "jetphi="    << fJetPhi   << // jetEta
        "jetArea="   << JetAreaPt  <<
        "jetNumberOfConstituents="    << JetNumberOfConstituents  <<
        "particleMaxChargedPt="  << particleMaxChargedPt <<
        "jetptsub="  << jetptsub <<
        "dEdx="      << fTPCSignal            <<  //  dEdx of the track
        "sign="      << fSign                 <<  //  charge
        "ptot="      << fPtot                 <<  //  TPC momentum
        "p="         << fPVertex              <<  //  momentum at vertex
        "pT="        << fPt                   <<  // transverse momentum
        "eta="       << fEta                  <<  //  eta
        "phi="       << fPhi                  <<  //  phi
        "dEdxMeanEl="  << fDEdxEl              << //mean dEdx for electrons
        "dEdxSigmaEl=" << fSigmaEl            << //sigma dEdx for electrons
        "dEdxMeanPi="  << fDEdxPi              <<
        "dEdxSigmaPi=" << fSigmaPi            <<
        "dEdxMeanKa="  << fDEdxKa              <<
        "dEdxSigmaKa=" << fSigmaKa            <<
        "dEdxMeanPr="  << fDEdxPr              <<
        "dEdxSigmaPr=" << fSigmaPr            <<
        "cent="      << fCentrality;
        if (!fSmallOut){
          (*fTreeSRedirector)<<"jetsEMCconst"<<
          "primMult="  << fNContributors <<          //  #prim tracks
          "eltpcpid="  << fNSigmasElTPC         <<  // nsigma TPC for electrons
          "pitpcpid="  << fNSigmasPiTPC         <<
          "katpcpid="  << fNSigmasKaTPC         <<
          "prtpcpid="  << fNSigmasPrTPC         <<
          "tofSignal=" << tofSignal         <<
          "prtofpid="  << fNSigmasPrTOF<<
          "dEdxMeanDe="  << fDEdxDe              <<
          "dEdxSigmaDe=" << fSigmaDe            <<
          "beta="      << beta;
        }
        (*fTreeSRedirector)<<"jetsEMCconst"<<"\n";
      }
    }
    
  }

  if ( (fDoPerpCone || fDoRandCone) && (leadJetPtSub > pT_sub_min) ) {

    //now we find the eta,phi perp (and R.C.) to it
    //From the AliAnalysisTaskRhoPerpCone Task
    Double_t dPhi1 = 999.;
    Double_t dPhi2 = 999.;
    Double_t dEta = 999.;
    Double_t Axis1 = 999, Axis2 = 999;

    if (fDoRandCone){
      double OffsetRandom = (1 + fRandom->Rndm()) * TMath::Pi()/3; // creating random Phi in the range of PI/3 and 2PI/3
      Axis1 = leadJetPhi + OffsetRandom;  // adding the random Phi to the Phi of leading jet, to get the Phi of random cone
      Axis2 = leadJetPhi - OffsetRandom; // mirroring it, so you get Phi of the 2nd random cone
      if(Axis1 > TMath::TwoPi()) Axis1 -= TMath::TwoPi(); // checking if it's larger than 2PI (leading jet phi is distributed in range of 0 and 2PI
      if(Axis2 < 0) Axis2 += TMath::TwoPi(); // checking if 2nd cone is smaller than 0
    }

    if (fDoPerpCone){
      Axis1 = ((leadJetPhi + (TMath::Pi() / 2.)) > TMath::TwoPi()) ? leadJetPhi - ((3. / 2.) * TMath::Pi()) : leadJetPhi + (TMath::Pi() / 2.);
      Axis2 = ((leadJetPhi - (TMath::Pi() / 2.)) < 0) ? leadJetPhi + ((3. / 2.) * TMath::Pi()) : leadJetPhi - (TMath::Pi() / 2.);
    }

    //loop over all inclusive tracks
    AliVEvent *event=InputEvent();
    for (Int_t itrack=0;itrack<event->GetNumberOfTracks();++itrack) {   // Track loop
      //
      AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(itrack));
      if (!track) continue;
      //
      // --------------------------------------------------------------
      //      Get relevant track info and set cut bits
      // --------------------------------------------------------------
      //

      Bool_t ifDefaultCuts = track->TestFilterBit(fAOD_FilterBits);
      if (ifDefaultCuts != 1) continue;
      //
      SetCutBitsAndSomeTrackVariables(track);
      if (TMath::Abs(fEta) > fEtaCut) continue;
              
      Float_t mod_track_phi = track->Phi() + TMath::Pi();
      //Check if the track is within the R=0.4 cone in eta, phi
      dPhi1 = TMath::Abs(mod_track_phi - Axis1);
      dPhi1 = (dPhi1 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi1 : dPhi1;
      dPhi2 = TMath::Abs(mod_track_phi - Axis2);
      dPhi2 = (dPhi2 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi2 : dPhi2;
      dEta = leadJetEta - track->Eta();

      if ((TMath::Sqrt(dPhi1 * dPhi1 + dEta * dEta) > 0.4) && (TMath::Sqrt(dPhi2 * dPhi2 + dEta * dEta) > 0.4)) continue; //scale the yields by 1/(2Ncones*piR^2) offline

      //
      // --------------------------------------------------------------
      //   Fill the trees
      // --------------------------------------------------------------
      //

      fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
      fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
      // Get the track variables
      GetExpecteds(track);

      Float_t length = track->GetIntegratedLength();
      Float_t tofSignal = track->GetTOFsignal();
      Float_t t0 = fTOFPIDResponse.GetStartTime(fPVertex);
      Float_t beta = -.05;
      //int nTOFClusters = track->GetNTOFclusters(); //All matchable clusters ESD only
      if((length > 0) && (tofSignal > 0)) beta = length / (2.99792458e-2*(tofSignal - t0));

      if (fPt>100.0) continue; //So we can match the jets that we throw out w/ max track pT>100

      Float_t fTPCmom_choice = fPtot;
      if (fSetTPCmom == 1) fTPCmom_choice = fPVertex;
      if (fSetTPCmom == 2) fTPCmom_choice = fPt;

      Float_t fTOFmom_choice = fPt;
      if (fSetTOFmom == 1) fTPCmom_choice = fPVertex;
      if (fSetTOFmom == 2) fTPCmom_choice = fPtot;

      Float_t fBetamom_choice = fPVertex;
      if (fSetBetamom == 1) fTPCmom_choice = fPt;
      if (fSetBetamom == 2) fTPCmom_choice = fPtot;

      Float_t fEta_choice = fEta;
      if (fSetEta == 1) fEta_choice = fY;

      Float_t DCAxy = -999.;
      Float_t DCAz = -999.;
      track->GetImpactParameters(DCAxy, DCAz);
      Float_t comb_sig_pi = TMath::Sqrt(fNSigmasPiTPC * fNSigmasPiTPC + fNSigmasPiTOF * fNSigmasPiTOF);
      Float_t comb_sig_pr = TMath::Sqrt(fNSigmasPrTPC * fNSigmasPrTPC + fNSigmasPrTOF * fNSigmasPrTOF);

      if (ifDefaultCuts == 1 && TMath::Abs(fEta) < fEtaCut){
        if (fFill_TPC) fHistIncTracks_dEdx->Fill(fTPCmom_choice,fTPCSignal,TMath::Abs(fEta_choice));
        if (fFillpTPC_pT) fHistIncTracks_moms->Fill(fPt,fPtot);
        if (fFillp_pT) fHistIncTracks_moms_p->Fill(fPt,fPVertex);
        if (fFillpTPC_p) fHistIncTracks_moms_pTPC_p->Fill(fPtot,fPVertex);
        fHistIncTracks_kin->Fill(fPt,TMath::Abs(fEta),fPhi);
        if (comb_sig_pi < 2.0) fHist_pi_DCAxy->Fill(fPt, DCAxy);
        if (comb_sig_pr < 2.0) fHist_pr_DCAxy->Fill(fPt, DCAxy);

        if (fFill_TPC_expecs){
          fHistIncTracks_mpi->Fill(fTPCmom_choice,fDEdxPi,TMath::Abs(fEta_choice));
          fHistIncTracks_spi->Fill(fTPCmom_choice,fSigmaPi,TMath::Abs(fEta_choice));

          fHistIncTracks_mel->Fill(fTPCmom_choice,fDEdxEl,TMath::Abs(fEta_choice));
          fHistIncTracks_sel->Fill(fTPCmom_choice,fSigmaEl,TMath::Abs(fEta_choice));

          fHistIncTracks_mka->Fill(fTPCmom_choice,fDEdxKa,TMath::Abs(fEta_choice));
          fHistIncTracks_ska->Fill(fTPCmom_choice,fSigmaKa,TMath::Abs(fEta_choice));

          fHistIncTracks_mpr->Fill(fTPCmom_choice,fDEdxPr,TMath::Abs(fEta_choice));
          fHistIncTracks_spr->Fill(fTPCmom_choice,fSigmaPr,TMath::Abs(fEta_choice));
        }
      }
      //
      //TPC-TOF Matching conditions: Use standard TPC tracks, then require kTIME and kTOFout
      Bool_t fTOFout = kFALSE;
      if ((track->GetStatus() & AliAODTrack::kTOFout) != 0) { //Track has the kTOFout flag
      fTOFout = kTRUE;
      }

      //kTIME flag
      Bool_t fTime = kFALSE;
      if ((track->GetStatus() & AliAODTrack::kTIME) != 0) { //Track has the kTIME flag
      fTime = kTRUE;
      }

      if (ifDefaultCuts == 1 && TMath::Abs(fEta) < fEtaCut && fTOFout && fTime && fFill_TOF_expecs){
        Double_t inttime[5]; //6 is needed to account for earlier species - 0 = electron, 1 = muon, 2 = pion, 3 = kaon, 4 = proton, ESD only 5 = deuteron
        track->GetIntegratedTimes(inttime, 5); // Returns the array with integrated times for each particle hypothesis

        Float_t fTOFMismatchTime = -20000000.;
        fTOFMismatchTime = AliTOFPIDResponse::GetMismatchRandomValue(track->Eta());
        if (fTOFMismatchTime <= 0){
          fTOFMismatchTime = -20000000.;
        }

        for (Int_t i = 0; i < 3; i++) {
          const Double_t beta_expec = length / (2.99792458e-2 * (inttime[i+2] - t0));
          if (i==0){
            fHistBetaExpec_pi->Fill(fTPCmom_choice, beta_expec);
            fHist_pi_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion), TMath::Abs(fEta_choice));
            if (TMath::Abs(fNSigmasPiTOF) < 2.0){
              fHistTOFSigmaExpec_pi->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion), TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasElTOF) < 2.0){
              fHist_elExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasMuTOF) < 2.0){
              fHist_muExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasKaTOF) < 2.0){
              fHist_kaExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasPrTOF) < 2.0){
              fHist_prExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
            }
          }
          if (i==1){
            fHistBetaExpec_ka->Fill(fTPCmom_choice, beta_expec);
            fHist_ka_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
            if (TMath::Abs(fNSigmasKaTOF) < 2.0){
              fHistTOFSigmaExpec_ka->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasPiTOF) < 2.0){
              fHist_piExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasPrTOF) < 2.0){
              fHist_prExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
            }
          }
          if (i==2){
            fHistBetaExpec_pr->Fill(fTPCmom_choice, beta_expec);
            fHist_pr_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
            if (TMath::Abs(fNSigmasPrTOF) < 2.0){
              fHistTOFSigmaExpec_pr->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasPiTOF) < 2.0){
              fHist_piExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasKaTOF) < 2.0){
              fHist_kaExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasDeTOF) < 2.0){
              fHist_deExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
            }
          }
        }
      }

      if (ifDefaultCuts == 1 && TMath::Abs(fEta) < fEtaCut && fTOFout && fTime && fFill_TOF){
        fHistIncTracks_beta->Fill(fBetamom_choice,beta);
        fHistIncTracks_t0->Fill(fBetamom_choice,t0);
        fHistIncTracks_TOFpi_nsigma->Fill(fTPCmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
        fHistIncTracks_TOFka_nsigma->Fill(fTPCmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
        fHistIncTracks_TOFpr_nsigma->Fill(fTPCmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
        /*if (nTOFClusters < 2) { //ESD only
          fHistIncTracks_TOFpi_nsigma_1cls->Fill(fTPCmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
          fHistIncTracks_TOFka_nsigma_1cls->Fill(fTPCmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
          fHistIncTracks_TOFpr_nsigma_1cls->Fill(fTPCmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
        }*/
      }
    } //end inc track loop for perp cones
    fFilledUECone_Rec = kTRUE;
  }//end perp cone if

  if (fFillJetsEMCBG || fFillJetsEMCBGConst){
  fbgJetContainer = this->GetJetContainer("bgJets");
  //
  // Get some jet container properties
  Int_t Njetsbg         = fbgJetContainer->GetNJets();
  Int_t NAcceptedjetsbg = fbgJetContainer->GetNAcceptedJets();
  Float_t jetRadiusbg  = fbgJetContainer->GetJetRadius();
  Float_t jetEtaMinbg = fbgJetContainer->GetJetEtaMin();
  Float_t jetEtaMaxbg = fbgJetContainer->GetJetEtaMax();
  //
  // loop over jets
  for(auto jet : fbgJetContainer->accepted())
  {
    fJetPt = jet->Pt();
    fJetEta = jet->Eta();
    fJetPhi = jet->Phi();
    Float_t JetM = jet->M();
    Int_t JetLabel = jet->GetLabel();
    Float_t JetAreaPt = jet->AreaPt(); //jet transverse area
    Int_t JetNumberOfConstituents = jet->GetNumberOfParticleConstituents(); //tracks + clusters
    Bool_t IsJetMc = jet->IsMC(); //using >0.0 of jet area inside emcal
    Float_t particleMaxChargedPt = jet->MaxChargedPt(); //max pT of charged particle in jet
    Float_t jetMCPt = jet->MCPt();
    Float_t jetptsub = jet->PtSub(fjetRhoVal, kFALSE); //bg sub pt

    if (fFillJetsEMCBG){
    if(!fTreeSRedirector) return;
    (*fTreeSRedirector)<<"jetsEMCBG"<<
    "jetRadius=" << jetRadiusbg << // jet Radius
    "Njets="     << Njetsbg << // Number of jets in event
    "NAcceptedjets=" << NAcceptedjetsbg << // Number of accepted jets in event
    "jetEtaMinbg=" << jetEtaMinbg <<
    "jetEtaMaxbg=" << jetEtaMaxbg <<
    "jetRhoVal=" << fjetRhoVal <<
    "jetpt="     << fJetPt    << // jetPt
    "jeteta="    << fJetEta   << // jetEta
    "jetphi="    << fJetPhi   << // jetPhi
    "jetArea="   << JetAreaPt  <<
    "jetNumberOfConstituents="    << JetNumberOfConstituents  <<
    "particleMaxChargedPt="  << particleMaxChargedPt <<
    "jetptsub="  << jetptsub <<
    "cent="      << fCentrality  <<  //  centrality
    "\n";
    }

    for(Int_t i = 0; i < jet->GetNumberOfParticleConstituents(); i++)
    {
      const AliVParticle* particle = jet->Track(i);
      AliAODTrack* AODtrack = (AliAODTrack*)(particle);
      if (!AODtrack) continue;

      //Track cuts start
      fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
      fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
      //
      //if (!fAODtrackCuts->AcceptTrack(AODtrack))  continue;    // default cuts - redundant since these track cuts are passed to jet finder
      //
      // Get the track variables
      GetExpecteds(AODtrack);
      SetCutBitsAndSomeTrackVariables(AODtrack);
      Double_t length    = AODtrack->GetIntegratedLength();
      Double_t tofSignal = AODtrack->GetTOFsignal();
      Double_t t0 = fTOFPIDResponse.GetStartTime(fPVertex);
      Double_t beta = -.05;
      //int nTOFClusters = AODtrack->GetNTOFclusters(); //All matchable clusters
      if((length > 0) && (tofSignal > 0)) beta = length / (2.99792458e-2*(tofSignal - t0));

      // Fill track constituent information
      if (fFillJetsEMCBGConst){

        (*fTreeSRedirector)<<"jetsEMCBGconst"<<
        //
        "jetRadius=" << jetRadiusbg << // jet Radius
        "Njets="     << Njetsbg << // Number of jets in event
        "NAcceptedjets=" << NAcceptedjetsbg << // Number of accepted jets in event
        "jetRhoVal=" << fjetRhoVal <<
        "jetpt="     << fJetPt    << // jetPt
        "jeteta="    << fJetEta   << // jetEta
        "jetphi="    << fJetPhi   << // jetPhi
        "jetArea="   << JetAreaPt  <<
        "jetNumberOfConstituents="    << JetNumberOfConstituents  <<
        "particleMaxChargedPt="  << particleMaxChargedPt <<
        "jetptsub="  << jetptsub <<
        "dEdx="      << fTPCSignal            <<  //  dEdx of the track
        "sign="      << fSign                 <<  //  charge
        "ptot="      << fPtot                 <<  //  TPC momentum
        "p="         << fPVertex              <<  //  momentum at vertex
        "pT="        << fPt                   <<  // transverse momentum
        "eta="       << fEta                  <<  //  eta
        "phi="       << fPhi                  <<  //  phi
        "dEdxMeanEl="  << fDEdxEl              << //mean dEdx for electrons
        "dEdxSigmaEl=" << fSigmaEl            << //sigma dEdx for electrons
        "dEdxMeanPi="  << fDEdxPi              <<
        "dEdxSigmaPi=" << fSigmaPi            <<
        "dEdxMeanKa="  << fDEdxKa              <<
        "dEdxSigmaKa=" << fSigmaKa            <<
        "dEdxMeanPr="  << fDEdxPr              <<
        "dEdxSigmaPr=" << fSigmaPr            <<
        "cent="      << fCentrality;
        if (!fSmallOut){
          (*fTreeSRedirector)<<"jetsEMCBGconst"<<
          "primMult="  << fNContributors <<          //  #prim tracks
          "eltpcpid="  << fNSigmasElTPC         <<  // nsigma TPC for electrons
          "pitpcpid="  << fNSigmasPiTPC         <<
          "katpcpid="  << fNSigmasKaTPC         <<
          "prtpcpid="  << fNSigmasPrTPC         <<
          "tofSignal=" << tofSignal         <<
          "prtofpid="  << fNSigmasPrTOF<<
          "dEdxMeanDe="  << fDEdxDe              <<
          "dEdxSigmaDe=" << fSigmaDe            <<
          "beta="      << beta;
        }
        (*fTreeSRedirector)<<"jetsEMCBGconst"<<"\n";
      }
    }

  }
  }

}
//________________________________________________________________________
void AliAnalysisTaskJetHadroAOD::FindJetsFJ()
{
  //
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FindJetsFJ ===== " << std::endl;
  AliTOFPIDResponse fTOFPIDResponse = fPIDResponse->GetTOFResponse();
  //
  // Create jetwrapper with the same settings used in FindJetsEMC
  int nJetRadiusBins = 1;
  float fTrackPt = 0.15;
  float fTrackEta = fEtaCut;
  float fGhostArea = 0.005;
  float bgJetAbsEtaCut = 0.7;           // fixed
  float bgJetRadius = 0.2;              // fixed
  //
  float pT_sub_min = fjetMinPtSub;    // can be 40, 60, 80
  std::vector<float>  fJetRadius;
  std::vector<float>  fPtSubMin;
  fJetRadius.resize(nJetRadiusBins);
  fJetRadius[0] = 0.4;
  fJetRadius[1] = 0.2;
  fJetRadius[2] = 0.6;
  fPtSubMin.resize(nJetRadiusBins);

    for (int iJetRadius=0; iJetRadius<nJetRadiusBins; iJetRadius++){
      Float_t jetAbsEtaCut = fTrackEta-fJetRadius[iJetRadius];   // fixed

      //
      //SOME CODE FROM NIMA
      fFastJetWrapper->Clear();
      fFastJetWrapper->SetR(fJetRadius[iJetRadius]);
      fFastJetWrapper->SetAlgorithm(fastjet::JetAlgorithm::antikt_algorithm);
      fFastJetWrapper->SetRecombScheme(fastjet::RecombinationScheme::pt_scheme);
      //fFastJetWrapper->SetStrategy(fastjet::Strategy::Best); //this doesn't need to be here - in FJ constructor already
      fFastJetWrapper->SetGhostArea(fGhostArea);
      fFastJetWrapper->SetAreaType(fastjet::AreaType::active_area);
      fFastJetWrapper->SetMaxRap(1.0);

      fFastJetWrapperBG->Clear();
      fFastJetWrapperBG->SetR(bgJetRadius);
      fFastJetWrapperBG->SetAlgorithm(fastjet::JetAlgorithm::kt_algorithm);
      fFastJetWrapperBG->SetRecombScheme(fastjet::RecombinationScheme::pt_scheme);
      //fFastJetWrapperBG->SetStrategy(fastjet::Strategy::Best); //this doesn't need to be here - in FJ constructor already
      fFastJetWrapperBG->SetGhostArea(fGhostArea);
      fFastJetWrapperBG->SetAreaType(fastjet::AreaType::active_area_explicit_ghosts);
      fFastJetWrapperBG->SetMaxRap(1.0);

      //std::vector<fastjet::PseudoJet> particlesEmbeddedSubtracted; //will be filled with your subtracted event
      std::vector<fastjet::PseudoJet> particlesEmbedded; //fill this with your event
      float particleEtaCut = fEtaCut;
      //
      // loop over AOD tracks and add their four vector to wrapper --> identical to track container in EMC jet
      for (Int_t iTrack = 0; iTrack < fAOD->GetNumberOfTracks(); iTrack++) {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
        if(!track || !track->TestFilterBit(fAOD_FilterBits)) continue; //hybrid track cuts filterbit
        if (track->Pt() < fTrackPt || TMath::Abs(track->Eta()) >= particleEtaCut) continue;
        fFastJetWrapper->AddInputVector(track->Px(), track->Py(), track->Pz(), track->E(), iTrack);//TMath::Sqrt(track->P()*track->P()+0.13957*0.13957),iTrack);
        fFastJetWrapperBG->AddInputVector(track->Px(), track->Py(), track->Pz(), track->E(), iTrack);//TMath::Sqrt(track->P()*track->P()+0.13957*0.13957),iTrack);
        particlesEmbedded.push_back(fastjet::PseudoJet(track->Px(), track->Py(), track->Pz(), track->E()));// TMath::Sqrt(track->P()*track->P()+0.13957*0.13957) ) );
      }
      //
      // background jet definitions
      fastjet::JetMedianBackgroundEstimator bgE;
      fastjet::Selector selectorBG = !fastjet::SelectorNHardest(2); //set the max eta cut on the estimator, then get rid of 2 highest pt jets
      bgE.set_selector(selectorBG);

      fastjet::JetDefinition jetDefBG(fastjet::kt_algorithm, bgJetRadius, fastjet::pt_scheme, fastjet::Best); //define the kT jet finding which will do the average background estimation
      fastjet::GhostedAreaSpec ghostSpecBG(particleEtaCut, 1, fGhostArea); //this ghost area might be too small and increase processing time too much
      fastjet::AreaDefinition areaDefBG(fastjet::active_area_explicit_ghosts, ghostSpecBG);
      fastjet::ClusterSequenceArea cluster_seq_BG(particlesEmbedded, jetDefBG, areaDefBG);

      fastjet::Selector selectorBGjets = !fastjet::SelectorIsPureGhost() * fastjet::SelectorAbsEtaMax(bgJetAbsEtaCut) * fastjet::SelectorPtRange(fTrackPt, 1000.0);
      fFastJetWrapperBG->Run();
      std::vector<fastjet::PseudoJet> jetsBG = sorted_by_pt(selectorBGjets(fFastJetWrapperBG->GetInclusiveJets()));

      if (jetsBG.size() > 0) {
        bgE.set_jets(jetsBG);  // give the kT jets to the background estimator
        frhoFJ = bgE.rho();
      }
      if (fUseCouts) std::cout << "frhoFJ is " << frhoFJ << std::endl;
      //
      // start of background jet loop
      if (fFilljetsFJBGTree || fFillJetsFJBGConst){
        for (Int_t ijet=0; ijet<Int_t(jetsBG.size()); ijet++) {
          fastjet::PseudoJet jet = jetsBG[ijet];
          //if (jet.pt() < fTrackPt || jet.perp() > 1000.0 || TMath::Abs(jet.eta()) >= bgJetAbsEtaCut) continue; //redundant because of the selector
          Float_t jetpt = jet.pt();
          Float_t jetphi = jet.phi();
          Float_t jeteta = jet.eta();
          Float_t jetArea = jet.area();
          Float_t jetptsub = jetpt - frhoFJ*jetArea;
          Int_t jetNum = jetsBG.size();
          Int_t Njets = jetsBG.size();
          std::vector<fastjet::PseudoJet> constituents = sorted_by_pt(jet.constituents());
          Int_t nConstituents = constituents.size();
          fastjet::PseudoJet highestpt_const = constituents[0];
          if (highestpt_const.pt() > 100.0) continue;

          if (fFilljetsFJBGTree)
          {
            (*fTreeSRedirector)<<"jetsFJBG"<<
            "bjJetRadius="    << bgJetRadius << // jet Radius
            "bgJetAbsEtaCut=" << bgJetAbsEtaCut << //abs eta cut for jet
            "jetNum="         << jetNum <<    //  number of jets
            "jetpt="          << jetpt <<
            "jetphi="         << jetphi <<
            "jeteta="         << jeteta <<
            "jetptsub="       << jetptsub << //bg sub jet pt (pt - rho*Area)
            "nConstituents="  << nConstituents << //num consts
            "rhoFJ="          << frhoFJ << //event rho
            "jetArea="        << jetArea << //jet area
            "cent="           << fCentrality  <<  //  centrality
            "\n";
          }

          for(Int_t i = 0; i < nConstituents; i++)
          {
            fastjet::PseudoJet &constituent = constituents[i];
            Float_t pt = constituent.pt();
            if (pt<1.e-10) continue;
            Float_t phi = constituent.phi();
            Float_t eta = constituent.eta();

            if (fFillJetsFJBGConst){
              (*fTreeSRedirector)<<"jetsFJBGconst"<<
              "pT="        << pt                   <<  // transverse momentum
              "eta="       << eta                  <<  //  eta
              "phi="       << phi                  <<  //  phi
              "\n";
            }
          }
        } // end of background jet loop
      }
      //
      // background subtraction on the constituent level TODO
      // fastjet::contrib::ConstituentSubtractor subtractorConstituent(&bgE); //add the background estimator to the correct subtractor
      // subtractorConstituent.set_common_bge_for_rho_and_rhom(true); //CHECK : should not be the case since particles have mass
      // subtractorConstituent.set_max_standardDeltaR(0.25); // set the max event wise subtraction distance
      // particlesEmbeddedSubtracted = subtractorConstituent.subtract_event(particlesEmbedded, particleEtaCut); //perform subtraction and fill the subtracted event container
      //

      // run jet finder using wrapper
      fFastJetWrapper->Run();
      std::vector<fastjet::PseudoJet> jets = sorted_by_pt(fFastJetWrapper->GetInclusiveJets());
      Float_t leadJetPhi = 999;
      Float_t leadJetEta = 999;
      Float_t leadJetPtSub = -999;
      // start of jet loop
      for (Int_t ijet=0; ijet<Int_t(jets.size()); ijet++){
        fastjet::PseudoJet jet = jets[ijet];
        if (jet.pt() < fTrackPt || jet.perp() > 1000.0 || TMath::Abs(jet.eta()) >= jetAbsEtaCut) continue;
        Float_t jetpt = jet.pt();
        Float_t jetphi = jet.phi();
        Float_t jeteta = jet.eta();
        Float_t jetArea = jet.area();
        Float_t jetptsub = jetpt - frhoFJ*jetArea;

        Int_t jetNum = jets.size();
        std::vector<fastjet::PseudoJet> constituents = sorted_by_pt(jets[ijet].constituents());

        Int_t nConstituents = constituents.size();
        if (nConstituents<1) continue;
        fastjet::PseudoJet highestpt_const = constituents[0];
        if (highestpt_const.pt() > 100.0) continue;
        fhasAcceptedFJjet = 1;

        if (jetptsub > pT_sub_min && jetptsub < fjetMaxPtSub  && jetArea > fjetMinArea)
        {
          fhasRealFJjet = 1;
          fNumRealFJJets +=1;
        }

        if (fFillFastJet){
          (*fTreeSRedirector)<<"jetsFJ"<<
          "jetNum="       << jetNum <<    //  number of jets
          "jetpt="        << jetpt <<
          "jetphi="       << jetphi <<
          "jeteta="       << jeteta <<
          "nConst="       << nConstituents <<    //  global event ID //this includes ghosts! Be careful
          "cent="         << fCentrality           <<  //  centrality
          "jetptsub="     << jetptsub << //bg sub jet pt (pt - rho*Area)
          "rhoFJ="        << frhoFJ << //event rho
          "jetArea="      << jetArea << //jet area
          "\n";
        }

        fHistJet_ptsub_v_area->Fill(jetArea, jetptsub);

        if (jetArea <= fjetMinArea) continue;
        fHistJet_kin->Fill(jetptsub, jeteta, jetphi);
        fHistJet_moms->Fill(jetptsub, jetpt);

        if (jetptsub <= pT_sub_min || jetptsub >= fjetMaxPtSub) continue;

        if (jetptsub > leadJetPtSub){
          leadJetPtSub = jetptsub;
          leadJetEta = jeteta;
          leadJetPhi = jetphi;
        }

        for(Int_t i = 0; i < nConstituents; i++)
        {
          fastjet::PseudoJet &constituent = constituents[i];
          Float_t pt = constituent.pt();
          if (pt<1.e-10) continue;

          Int_t trackIndex = constituent.user_index();
          AliAODTrack* trackConst = static_cast<AliAODTrack*>(fAOD->GetTrack(trackIndex));

          //Track cuts start
          fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
          fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
          //
          // --------------------------------------------------------------
          //      Get relevant track info and set cut bits
          // --------------------------------------------------------------
          //
          if(!trackConst || !trackConst->TestFilterBit(fAOD_FilterBits)) continue;
          Bool_t ifDefaultCuts = 1.0;
          //
          // Get the track variables
          GetExpecteds(trackConst);
          SetCutBitsAndSomeTrackVariables(trackConst);

          Float_t length    = trackConst->GetIntegratedLength();
          Float_t tofSignal = trackConst->GetTOFsignal();
          Float_t t0 = fTOFPIDResponse.GetStartTime(fPVertex);
          Float_t beta = -.05;
          //int nTOFClusters = trackConst->GetNTOFclusters(); //All matchable clusters ESD only
          if((length > 0) && (tofSignal > 0)) beta = length / (2.99792458e-2*(tofSignal - t0));

          Float_t fTPCmom_choice = fPtot;
          if (fSetTPCmom == 1) fTPCmom_choice = fPVertex;
          if (fSetTPCmom == 2) fTPCmom_choice = fPt;

          Float_t fTOFmom_choice = fPt;
          if (fSetTOFmom == 1) fTPCmom_choice = fPVertex;
          if (fSetTOFmom == 2) fTPCmom_choice = fPtot;

          Float_t fBetamom_choice = fPVertex;
          if (fSetBetamom == 1) fTPCmom_choice = fPt;
          if (fSetBetamom == 2) fTPCmom_choice = fPtot;

          Float_t fEta_choice = fEta;
          if (fSetEta == 1) fEta_choice = fY;

          Float_t DCAxy = -999.;
          Float_t DCAz = -999.;
          trackConst->GetImpactParameters(DCAxy, DCAz);
          Float_t comb_sig_pi = TMath::Sqrt(fNSigmasPiTPC * fNSigmasPiTPC + fNSigmasPiTOF * fNSigmasPiTOF);
          Float_t comb_sig_pr = TMath::Sqrt(fNSigmasPrTPC * fNSigmasPrTPC + fNSigmasPrTOF * fNSigmasPrTOF);

          if (ifDefaultCuts == 1 && TMath::Abs(fEta) < fEtaCut){
            if (fFill_TPC) fHistJetTracks_dEdx->Fill(fTPCmom_choice,fTPCSignal,TMath::Abs(fEta_choice));
            if (fFillpTPC_pT) fHistJetTracks_moms->Fill(fPt,fPtot);
            if (fFillp_pT) fHistJetTracks_moms_p->Fill(fPt,fPVertex);
            if (fFillpTPC_p) fHistJetTracks_moms_pTPC_p->Fill(fPtot,fPVertex);
            fHistJetTracks_kin->Fill(fPt,TMath::Abs(fEta),fPhi);
            if (comb_sig_pi < 2.0) fHistJet_pi_DCAxy->Fill(fPt, DCAxy);
            if (comb_sig_pr < 2.0) fHistJet_pr_DCAxy->Fill(fPt, DCAxy);
          }

          //TPC-TOF Matching conditions: Use standard TPC tracks, then require kTIME and kTOFout
          Bool_t fTOFout = kFALSE;
          if ((trackConst->GetStatus() & AliAODTrack::kTOFout) != 0) { //Track has the kTOFout flag
          fTOFout = kTRUE;
          }

          //
          //kTIME flag
          Bool_t fTime = kFALSE;
          if ((trackConst->GetStatus() & AliAODTrack::kTIME) != 0) { //Track has the kTIME flag
          fTime = kTRUE;
          }

          Double_t inttime[5]; //6 is needed to account for earlier species - 0 = electron, 1 = muon, 2 = pion, 3 = kaon, 4 = proton, ESD only: 5 = deuteron
          trackConst->GetIntegratedTimes(inttime, 5); // Returns the array with integrated times for each particle hypothesis

          Float_t fTOFMismatchTime = -20000000.;
          fTOFMismatchTime = AliTOFPIDResponse::GetMismatchRandomValue(trackConst->Eta());
          if (fTOFMismatchTime <= 0){
            fTOFMismatchTime = -20000000.;
          }

          for (Int_t i = 0; i < 3; i++) {
            const Float_t beta_expec = length / (2.99792458e-2 * (inttime[i+2] - t0));
            if (i==0){
              fHistjet_BetaExpec_pi->Fill(fTPCmom_choice, beta_expec);
              fHist_jet_pi_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
              if (TMath::Abs(fNSigmasPiTOF) < 2.0){
                fHistjet_TOFSigmaExpec_pi->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasElTOF) < 2.0){
                fHist_jet_elExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasMuTOF) < 2.0){
                fHist_jet_muExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasKaTOF) < 2.0){
                fHist_jet_kaExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasPrTOF) < 2.0){
                fHist_jet_prExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
              }
            }
            if (i==1){
              fHistjet_BetaExpec_ka->Fill(fTPCmom_choice, beta_expec);
              fHist_jet_ka_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
              if (TMath::Abs(fNSigmasKaTOF) < 2.0){
                fHistjet_TOFSigmaExpec_ka->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasPiTOF) < 2.0){
                fHist_jet_piExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasPrTOF) < 2.0){
                fHist_jet_prExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
              }
            }
            if (i==2){
              fHistjet_BetaExpec_pr->Fill(fTPCmom_choice, beta_expec);
              fHist_jet_pr_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
              if (TMath::Abs(fNSigmasPrTOF) < 2.0){
                fHistjet_TOFSigmaExpec_pr->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasPiTOF) < 2.0){
                fHist_jet_piExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasKaTOF) < 2.0){
                fHist_jet_kaExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasDeTOF) < 2.0){
                fHist_jet_deExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
              }
            }
          }

          if (ifDefaultCuts == 1 && TMath::Abs(fEta) < fEtaCut  && fTOFout && fTime && fFill_TOF){
            fHistJetTracks_beta->Fill(fBetamom_choice,beta);
            fHistJetTracks_TOFpi_nsigma->Fill(fTOFmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
            fHistJetTracks_TOFka_nsigma->Fill(fTOFmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
            fHistJetTracks_TOFpr_nsigma->Fill(fTOFmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
            /*if (nTOFClusters < 2) { //ESD only
              fHistJetTracks_TOFpi_nsigma_1cls->Fill(fTOFmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
              fHistJetTracks_TOFka_nsigma_1cls->Fill(fTOFmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
              fHistJetTracks_TOFpr_nsigma_1cls->Fill(fTOFmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
            }*/
          }

          if (fFillJetsFJConst && fTOFout && fTime && ifDefaultCuts)
          {

          (*fTreeSRedirector)<<"jetsFJconst"<<
          "jetRadius=" << fJetRadius[iJetRadius] << // jet Radius
          "jetNum="    << jetNum <<    //  number of jets
          "jetpt="     << jetpt <<     //  global event ID
          "jetphi="    << jetphi <<    //  global event ID
          "jeteta="    << jeteta <<    //  global event ID
          "nConst="    << nConstituents <<    //  global event ID

          "jetptsub="  << jetptsub << //bg sub jet pt (pt - rho*Area)
          "rhoFJ="     << frhoFJ << //event rho
          "jetArea="   << jetArea << //jet area

          //"defCut="    << ifDefaultCuts <<  // default cuts tuned by hand

          "dEdx="      << fTPCSignal            <<  //  dEdx of the track
          "sign="      << fSign                 <<  //  charge
          "ptot="      << fPtot                 <<  //  TPC momentum
          "p="         << fPVertex              <<  //  momentum at vertex
          "pT="        << fPt                   <<  // transverse momentum
          "eta="       << fEta                  <<  //  eta
          "cent="      << fCentrality           <<  //  centrality
          "phi="       << fPhi                  <<  //  phi

          "dEdxMeanEl=" << fDEdxEl              << //mean dEdx for electrons
          "dEdxSigmaEl=" << fSigmaEl            << //sigma dEdx for electrons
          "dEdxMeanPi=" << fDEdxPi              <<
          "dEdxSigmaPi=" << fSigmaPi            <<
          "dEdxMeanKa=" << fDEdxKa              <<
          "dEdxSigmaKa=" << fSigmaKa            <<
          "dEdxMeanPr=" << fDEdxPr              <<
          "dEdxSigmaPr=" << fSigmaPr            <<
          "dEdxMeanDe=" << fDEdxDe              <<
          "dEdxSigmaDe=" << fSigmaDe            <<

          "beta=" << beta         <<

          "eltofpid="  << fNSigmasElTOF         <<  // nsigma TPC for electrons
          "pitofpid="  << fNSigmasPiTOF         <<
          "katofpid="  << fNSigmasKaTOF         <<
          "prtofpid="  << fNSigmasPrTOF         <<
          "detofpid="  << fNSigmasDeTOF;

          if (!fSmallOut){
            (*fTreeSRedirector)<<"jetsFJconst"<<
            //"bit128="    << fBit128 <<        // TPC only tracks cuts
            //"bit768="    << fBit768 <<        // Hybrid track cuts
            "primMult="  << fNContributors <<          //  #prim tracks
            "eltpcpid="  << fNSigmasElTPC         <<  // nsigma TPC for electrons
            "pitpcpid="  << fNSigmasPiTPC         <<
            "katpcpid="  << fNSigmasKaTPC         <<
            "prtpcpid="  << fNSigmasPrTPC         <<
            "detpcpid="  << fNSigmasDeTPC         <<

            "tofSignal=" << tofSignal;
          }
        (*fTreeSRedirector)<<"jetsFJconst"<<"\n";
        }
        }
      } // end of jet loop

      if ( (fDoPerpCone || fDoRandCone) && (leadJetPtSub > pT_sub_min) ) {
        //now we find the eta,phi perp (and R.C.) to it
        //From the AliAnalysisTaskRhoPerpCone Task
        Double_t dPhi1 = 999.;
        Double_t dPhi2 = 999.;
        Double_t dEta = 999.;
        Double_t Axis1 = 999, Axis2 = 999;

        if (fDoRandCone){
          double OffsetRandom = (1 + fRandom->Rndm()) * TMath::Pi()/3; // creating random Phi in the range of PI/3 and 2PI/3
          Axis1 = leadJetPhi + OffsetRandom;  // adding the random Phi to the Phi of leading jet, to get the Phi of random cone
          Axis2 = leadJetPhi - OffsetRandom; // mirroring it, so you get Phi of the 2nd random cone
          if(Axis1 > TMath::TwoPi()) Axis1 -= TMath::TwoPi(); // checking if it's larger than 2PI (leading jet phi is distributed in range of 0 and 2PI
          if(Axis2 < 0) Axis2 += TMath::TwoPi(); // checking if 2nd cone is smaller than 0
        }

        if (fDoPerpCone){
          Axis1 = ((leadJetPhi + (TMath::Pi() / 2.)) > TMath::TwoPi()) ? leadJetPhi - ((3. / 2.) * TMath::Pi()) : leadJetPhi + (TMath::Pi() / 2.);
          Axis2 = ((leadJetPhi - (TMath::Pi() / 2.)) < 0) ? leadJetPhi + ((3. / 2.) * TMath::Pi()) : leadJetPhi - (TMath::Pi() / 2.);
        }

        //loop over all inclusive tracks
        AliVEvent *event=InputEvent();
        for (Int_t itrack=0;itrack<event->GetNumberOfTracks();++itrack) {   // Track loop
          //
          AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(itrack));
          if (!track) continue;
          //
          // --------------------------------------------------------------
          //      Get relevant track info and set cut bits
          // --------------------------------------------------------------
          //

          Bool_t ifDefaultCuts = track->TestFilterBit(fAOD_FilterBits);
          if (ifDefaultCuts != 1) continue;
          //
          SetCutBitsAndSomeTrackVariables(track);
          if (TMath::Abs(fEta) > fEtaCut) continue;
                  
          Float_t mod_track_phi = track->Phi() + TMath::Pi();
          //Check if the track is within the R=0.4 cone in eta, phi
          dPhi1 = TMath::Abs(mod_track_phi - Axis1);
          dPhi1 = (dPhi1 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi1 : dPhi1;
          dPhi2 = TMath::Abs(mod_track_phi - Axis2);
          dPhi2 = (dPhi2 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi2 : dPhi2;
          dEta = leadJetEta - track->Eta();

          if ((TMath::Sqrt(dPhi1 * dPhi1 + dEta * dEta) > 0.4) && (TMath::Sqrt(dPhi2 * dPhi2 + dEta * dEta) > 0.4)) continue; //scale the yields by 1/(2Ncones*piR^2) offline

          //
          // --------------------------------------------------------------
          //   Fill the trees
          // --------------------------------------------------------------
          //

          fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
          fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
          // Get the track variables
          GetExpecteds(track);

          Float_t length = track->GetIntegratedLength();
          Float_t tofSignal = track->GetTOFsignal();
          Float_t t0 = fTOFPIDResponse.GetStartTime(fPVertex);
          Float_t beta = -.05;
          //int nTOFClusters = track->GetNTOFclusters(); //All matchable clusters ESD only
          if((length > 0) && (tofSignal > 0)) beta = length / (2.99792458e-2*(tofSignal - t0));

          if (fPt>100.0) continue; //So we can match the jets that we throw out w/ max track pT>100

          Float_t fTPCmom_choice = fPtot;
          if (fSetTPCmom == 1) fTPCmom_choice = fPVertex;
          if (fSetTPCmom == 2) fTPCmom_choice = fPt;

          Float_t fTOFmom_choice = fPt;
          if (fSetTOFmom == 1) fTPCmom_choice = fPVertex;
          if (fSetTOFmom == 2) fTPCmom_choice = fPtot;

          Float_t fBetamom_choice = fPVertex;
          if (fSetBetamom == 1) fTPCmom_choice = fPt;
          if (fSetBetamom == 2) fTPCmom_choice = fPtot;

          Float_t fEta_choice = fEta;
          if (fSetEta == 1) fEta_choice = fY;

          Float_t DCAxy = -999.;
          Float_t DCAz = -999.;
          track->GetImpactParameters(DCAxy, DCAz);
          Float_t comb_sig_pi = TMath::Sqrt(fNSigmasPiTPC * fNSigmasPiTPC + fNSigmasPiTOF * fNSigmasPiTOF);
          Float_t comb_sig_pr = TMath::Sqrt(fNSigmasPrTPC * fNSigmasPrTPC + fNSigmasPrTOF * fNSigmasPrTOF);

          if (ifDefaultCuts == 1 && TMath::Abs(fEta) < fEtaCut){
            if (fFill_TPC) fHistIncTracks_dEdx->Fill(fTPCmom_choice,fTPCSignal,TMath::Abs(fEta_choice));
            if (fFillpTPC_pT) fHistIncTracks_moms->Fill(fPt,fPtot);
            if (fFillp_pT) fHistIncTracks_moms_p->Fill(fPt,fPVertex);
            if (fFillpTPC_p) fHistIncTracks_moms_pTPC_p->Fill(fPtot,fPVertex);
            fHistIncTracks_kin->Fill(fPt,TMath::Abs(fEta),fPhi);
            if (comb_sig_pi < 2.0) fHist_pi_DCAxy->Fill(fPt, DCAxy);
            if (comb_sig_pr < 2.0) fHist_pr_DCAxy->Fill(fPt, DCAxy);

            if (fFill_TPC_expecs){
              fHistIncTracks_mpi->Fill(fTPCmom_choice,fDEdxPi,TMath::Abs(fEta_choice));
              fHistIncTracks_spi->Fill(fTPCmom_choice,fSigmaPi,TMath::Abs(fEta_choice));

              fHistIncTracks_mel->Fill(fTPCmom_choice,fDEdxEl,TMath::Abs(fEta_choice));
              fHistIncTracks_sel->Fill(fTPCmom_choice,fSigmaEl,TMath::Abs(fEta_choice));

              fHistIncTracks_mka->Fill(fTPCmom_choice,fDEdxKa,TMath::Abs(fEta_choice));
              fHistIncTracks_ska->Fill(fTPCmom_choice,fSigmaKa,TMath::Abs(fEta_choice));

              fHistIncTracks_mpr->Fill(fTPCmom_choice,fDEdxPr,TMath::Abs(fEta_choice));
              fHistIncTracks_spr->Fill(fTPCmom_choice,fSigmaPr,TMath::Abs(fEta_choice));
            }
          }
          //
          //TPC-TOF Matching conditions: Use standard TPC tracks, then require kTIME and kTOFout
          Bool_t fTOFout = kFALSE;
          if ((track->GetStatus() & AliAODTrack::kTOFout) != 0) { //Track has the kTOFout flag
          fTOFout = kTRUE;
          }

          //kTIME flag
          Bool_t fTime = kFALSE;
          if ((track->GetStatus() & AliAODTrack::kTIME) != 0) { //Track has the kTIME flag
          fTime = kTRUE;
          }

          if (ifDefaultCuts == 1 && TMath::Abs(fEta) < fEtaCut && fTOFout && fTime && fFill_TOF_expecs){
            Double_t inttime[5]; //6 is needed to account for earlier species - 0 = electron, 1 = muon, 2 = pion, 3 = kaon, 4 = proton, ESD only 5 = deuteron
            track->GetIntegratedTimes(inttime, 5); // Returns the array with integrated times for each particle hypothesis

            Float_t fTOFMismatchTime = -20000000.;
            fTOFMismatchTime = AliTOFPIDResponse::GetMismatchRandomValue(track->Eta());
            if (fTOFMismatchTime <= 0){
              fTOFMismatchTime = -20000000.;
            }

            for (Int_t i = 0; i < 3; i++) {
              const Double_t beta_expec = length / (2.99792458e-2 * (inttime[i+2] - t0));
              if (i==0){
                fHistBetaExpec_pi->Fill(fTPCmom_choice, beta_expec);
                fHist_pi_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion), TMath::Abs(fEta_choice));
                if (TMath::Abs(fNSigmasPiTOF) < 2.0){
                  fHistTOFSigmaExpec_pi->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion), TMath::Abs(fEta_choice));
                }
                if (TMath::Abs(fNSigmasElTOF) < 2.0){
                  fHist_elExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
                }
                if (TMath::Abs(fNSigmasMuTOF) < 2.0){
                  fHist_muExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
                }
                if (TMath::Abs(fNSigmasKaTOF) < 2.0){
                  fHist_kaExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
                }
                if (TMath::Abs(fNSigmasPrTOF) < 2.0){
                  fHist_prExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
                }
              }
              if (i==1){
                fHistBetaExpec_ka->Fill(fTPCmom_choice, beta_expec);
                fHist_ka_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
                if (TMath::Abs(fNSigmasKaTOF) < 2.0){
                  fHistTOFSigmaExpec_ka->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
                }
                if (TMath::Abs(fNSigmasPiTOF) < 2.0){
                  fHist_piExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
                }
                if (TMath::Abs(fNSigmasPrTOF) < 2.0){
                  fHist_prExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
                }
              }
              if (i==2){
                fHistBetaExpec_pr->Fill(fTPCmom_choice, beta_expec);
                fHist_pr_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
                if (TMath::Abs(fNSigmasPrTOF) < 2.0){
                  fHistTOFSigmaExpec_pr->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
                }
                if (TMath::Abs(fNSigmasPiTOF) < 2.0){
                  fHist_piExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
                }
                if (TMath::Abs(fNSigmasKaTOF) < 2.0){
                  fHist_kaExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
                }
                if (TMath::Abs(fNSigmasDeTOF) < 2.0){
                  fHist_deExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
                }
              }
            }
          }

          if (ifDefaultCuts == 1 && TMath::Abs(fEta) < fEtaCut && fTOFout && fTime && fFill_TOF){
            fHistIncTracks_beta->Fill(fBetamom_choice,beta);
            fHistIncTracks_t0->Fill(fBetamom_choice,t0);
            fHistIncTracks_TOFpi_nsigma->Fill(fTPCmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
            fHistIncTracks_TOFka_nsigma->Fill(fTPCmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
            fHistIncTracks_TOFpr_nsigma->Fill(fTPCmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
            /*if (nTOFClusters < 2) { //ESD only
              fHistIncTracks_TOFpi_nsigma_1cls->Fill(fTPCmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
              fHistIncTracks_TOFka_nsigma_1cls->Fill(fTPCmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
              fHistIncTracks_TOFpr_nsigma_1cls->Fill(fTPCmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
            }*/
          }
        } //end inc track loop for perp cones
        fFilledUECone_Rec = kTRUE;
      }//end perp cone if

  } //end jet radius loop
}
//________________________________________________________________________
void AliAnalysisTaskJetHadroAOD::FillIncTracksReal()
{
  //
  // Fill dEdx information for the TPC and also clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FillIncTracksReal ===== " << std::endl;
  AliTOFPIDResponse fTOFPIDResponse = fPIDResponse->GetTOFResponse();
  // --------------------------------------------------------------
  // Get the event
  AliVEvent *event=InputEvent();
  if (CountEmptyEvents()) return;
  fisGoodIncEvent = 1;
  //
  // --------------------------------------------------------------
  //  Main track loop
  // --------------------------------------------------------------
  //
  //
  for (Int_t itrack=0;itrack<event->GetNumberOfTracks();++itrack) {   // Track loop

    fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
    fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
    //
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(itrack));
    if (!track) continue;
    //
    // --------------------------------------------------------------
    //      Get relevant track info and set cut bits
    // --------------------------------------------------------------
    //

    Bool_t ifDefaultCuts = track->TestFilterBit(fAOD_FilterBits);
    //
    // Get the track variables
    GetExpecteds(track);
    SetCutBitsAndSomeTrackVariables(track);
    //
    // --------------------------------------------------------------
    //   Fill the trees
    // --------------------------------------------------------------
    //

    Float_t length = track->GetIntegratedLength();
    Float_t tofSignal = track->GetTOFsignal();
    Float_t t0 = fTOFPIDResponse.GetStartTime(fPVertex);
    Float_t beta = -.05;
    //int nTOFClusters = track->GetNTOFclusters(); //All matchable clusters ESD only
    if((length > 0) && (tofSignal > 0)) beta = length / (2.99792458e-2*(tofSignal - t0));

    if (fPt>100.0) continue; //So we can match the jets that we throw out w/ max track pT>100

    Float_t fTPCmom_choice = fPtot;
    if (fSetTPCmom == 1) fTPCmom_choice = fPVertex;
    if (fSetTPCmom == 2) fTPCmom_choice = fPt;

    Float_t fTOFmom_choice = fPt;
    if (fSetTOFmom == 1) fTPCmom_choice = fPVertex;
    if (fSetTOFmom == 2) fTPCmom_choice = fPtot;

    Float_t fBetamom_choice = fPVertex;
    if (fSetBetamom == 1) fTPCmom_choice = fPt;
    if (fSetBetamom == 2) fTPCmom_choice = fPtot;

    Float_t fEta_choice = fEta;
    if (fSetEta == 1) fEta_choice = fY;

    Float_t DCAxy = -999.;
    Float_t DCAz = -999.;
    track->GetImpactParameters(DCAxy, DCAz);
    Float_t comb_sig_pi = TMath::Sqrt(fNSigmasPiTPC * fNSigmasPiTPC + fNSigmasPiTOF * fNSigmasPiTOF);
    Float_t comb_sig_pr = TMath::Sqrt(fNSigmasPrTPC * fNSigmasPrTPC + fNSigmasPrTOF * fNSigmasPrTOF);

    if (ifDefaultCuts == 1 && TMath::Abs(fEta) < fEtaCut && TMath::Abs(fY) < fYCut){
      if (fFill_TPC) fHistIncTracks_dEdx->Fill(fTPCmom_choice,fTPCSignal,TMath::Abs(fEta_choice));
      if (fFillpTPC_pT) fHistIncTracks_moms->Fill(fPt,fPtot);
      if (fFillp_pT) fHistIncTracks_moms_p->Fill(fPt,fPVertex);
      if (fFillpTPC_p) fHistIncTracks_moms_pTPC_p->Fill(fPtot,fPVertex);
      fHistIncTracks_kin->Fill(fPt,TMath::Abs(fEta),fPhi);
      if (comb_sig_pi < 2.0) fHist_pi_DCAxy->Fill(fPt, DCAxy);
      if (comb_sig_pr < 2.0) fHist_pr_DCAxy->Fill(fPt, DCAxy);

      if (fFill_TPC_expecs){
        fHistIncTracks_mpi->Fill(fTPCmom_choice,fDEdxPi,TMath::Abs(fEta_choice));
        fHistIncTracks_spi->Fill(fTPCmom_choice,fSigmaPi,TMath::Abs(fEta_choice));

        fHistIncTracks_mel->Fill(fTPCmom_choice,fDEdxEl,TMath::Abs(fEta_choice));
        fHistIncTracks_sel->Fill(fTPCmom_choice,fSigmaEl,TMath::Abs(fEta_choice));

        fHistIncTracks_mka->Fill(fTPCmom_choice,fDEdxKa,TMath::Abs(fEta_choice));
        fHistIncTracks_ska->Fill(fTPCmom_choice,fSigmaKa,TMath::Abs(fEta_choice));

        fHistIncTracks_mpr->Fill(fTPCmom_choice,fDEdxPr,TMath::Abs(fEta_choice));
        fHistIncTracks_spr->Fill(fTPCmom_choice,fSigmaPr,TMath::Abs(fEta_choice));
      }
    }
    //

    //TPC-TOF Matching conditions: Use standard TPC tracks, then require kTIME and kTOFout
    Bool_t fTOFout = kFALSE;
    if ((track->GetStatus() & AliAODTrack::kTOFout) != 0) { //Track has the kTOFout flag
    fTOFout = kTRUE;
    }

    //kTIME flag
    Bool_t fTime = kFALSE;
    if ((track->GetStatus() & AliAODTrack::kTIME) != 0) { //Track has the kTIME flag
    fTime = kTRUE;
    }

    if (ifDefaultCuts == 1 && TMath::Abs(fEta) < fEtaCut && TMath::Abs(fY) < fYCut && fTOFout && fTime && fFill_TOF_expecs){
      Double_t inttime[5]; //6 is needed to account for earlier species - 0 = electron, 1 = muon, 2 = pion, 3 = kaon, 4 = proton, ESD only 5 = deuteron
      track->GetIntegratedTimes(inttime, 5); // Returns the array with integrated times for each particle hypothesis

      Float_t fTOFMismatchTime = -20000000.;
      fTOFMismatchTime = AliTOFPIDResponse::GetMismatchRandomValue(track->Eta());
      if (fTOFMismatchTime <= 0){
        fTOFMismatchTime = -20000000.;
      }

      for (Int_t i = 0; i < 3; i++) {
        const Double_t beta_expec = length / (2.99792458e-2 * (inttime[i+2] - t0));
        if (i==0){
          fHistBetaExpec_pi->Fill(fTPCmom_choice, beta_expec);
          fHist_pi_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion), TMath::Abs(fEta_choice));
          if (TMath::Abs(fNSigmasPiTOF) < 2.0){
            fHistTOFSigmaExpec_pi->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion), TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasElTOF) < 2.0){
            fHist_elExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasMuTOF) < 2.0){
            fHist_muExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasKaTOF) < 2.0){
            fHist_kaExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasPrTOF) < 2.0){
            fHist_prExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
          }
        }
        if (i==1){
          fHistBetaExpec_ka->Fill(fTPCmom_choice, beta_expec);
          fHist_ka_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
          if (TMath::Abs(fNSigmasKaTOF) < 2.0){
            fHistTOFSigmaExpec_ka->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasPiTOF) < 2.0){
            fHist_piExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasPrTOF) < 2.0){
            fHist_prExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
          }
        }
        if (i==2){
          fHistBetaExpec_pr->Fill(fTPCmom_choice, beta_expec);
          fHist_pr_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
          if (TMath::Abs(fNSigmasPrTOF) < 2.0){
            fHistTOFSigmaExpec_pr->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasPiTOF) < 2.0){
            fHist_piExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasKaTOF) < 2.0){
            fHist_kaExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
          }
          if (TMath::Abs(fNSigmasDeTOF) < 2.0){
            fHist_deExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
          }
        }
      }
    }

    if (ifDefaultCuts == 1 && TMath::Abs(fEta) < fEtaCut && TMath::Abs(fY) < fYCut && fTOFout && fTime && fFill_TOF){
      fHistIncTracks_beta->Fill(fBetamom_choice,beta);
      fHistIncTracks_t0->Fill(fBetamom_choice,t0);
      fHistIncTracks_TOFpi_nsigma->Fill(fTPCmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
      fHistIncTracks_TOFka_nsigma->Fill(fTPCmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
      fHistIncTracks_TOFpr_nsigma->Fill(fTPCmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
      /*if (nTOFClusters < 2) { //ESD only
        fHistIncTracks_TOFpi_nsigma_1cls->Fill(fTPCmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
        fHistIncTracks_TOFka_nsigma_1cls->Fill(fTPCmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
        fHistIncTracks_TOFpr_nsigma_1cls->Fill(fTPCmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
      }*/
    }

	  if (fFillIncTracks && fTOFout && fTime && ifDefaultCuts)
    {
      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"tracks"<<
      //
      "defCut="    << ifDefaultCuts <<  // default cuts tuned by hand
      "dEdx="      << fTPCSignal            <<  //  dEdx of the track
      "sign="      << fSign                 <<  //  charge
      "ptot="      << fPtot                 <<  //  TPC momentum
      "p="         << fPVertex              <<  //  momentum at vertex
      "pT="        << fPt                   <<  // transverse momentum
      "eta="       << fEta                  <<  //  eta
      "phi="       << fPhi                  <<  //  phi
      "dEdxMeanEl=" << fDEdxEl              << //mean dEdx for electrons
      "dEdxSigmaEl=" << fSigmaEl            << //sigma dEdx for electrons
      "dEdxMeanPi=" << fDEdxPi              <<
      "dEdxSigmaPi=" << fSigmaPi            <<
      "dEdxMeanKa=" << fDEdxKa              <<
      "dEdxSigmaKa=" << fSigmaKa            <<
      "dEdxMeanPr=" << fDEdxPr              <<
      "dEdxSigmaPr=" << fSigmaPr            <<
      "cent="      << fCentrality           <<
      "eltofpid="  << fNSigmasElTOF         <<  // nsigma TPC for electrons
      "pitofpid="  << fNSigmasPiTOF         <<
      "katofpid="  << fNSigmasKaTOF         <<
      "prtofpid="  << fNSigmasPrTOF         <<
      "detofpid="  << fNSigmasDeTOF         <<
      "beta="      << beta;

      if (!fSmallOut){
        (*fTreeSRedirector)<<"tracks"<<
        "primMult="  << fNContributors <<          //  #prim tracks
        "eltpcpid="  << fNSigmasElTPC         <<  // nsigma TPC for electrons
        "pitpcpid="  << fNSigmasPiTPC         <<
        "katpcpid="  << fNSigmasKaTPC         <<
        "prtpcpid="  << fNSigmasPrTPC         <<
        "tofSignal=" << tofSignal         <<
        "prtofpid="  << fNSigmasPrTOF<<
        "dEdxMeanDe=" << fDEdxDe              <<
        "dEdxSigmaDe=" << fSigmaDe;

      }
      (*fTreeSRedirector)<<"tracks"<<"\n";
    }

  }// end of track loop

}
//________________________________________________________________________
void AliAnalysisTaskJetHadroAOD::FillTreeMC()
{
  //Some parts from AliAnalysisTaskTOFSpectra
  //Track eff = Num1/Den
  //TOF Match Eff = Num2/Num1

  //Den = # MC Truth particles that pass: 1. Is Primary 2. eta cut 3. pT cut 4. rap cut (if applicable) 5. Event selections (same as reco) 6. MC Truth PDG code
  //Num1 = # MC Reco particles that pass: 1. Is Primary 2. All data track selection cuts 3. eta cut 3. pT cut 4. rap cut (if applicable) 5. Event selections 6. MC Truth PDG code
  //Num2 = # MC Reco particles that pass: all of Num1 + 2. TOFMatch==0

  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FillTreeMC ===== " << std::endl;

  AliTOFPIDResponse fTOFPIDResponse = fPIDResponse->GetTOFResponse();

  //Int_t numMCTracks = fMCEvent->GetNumberOfTracks();      //ESD Number of gen particles - used for looping over tracks
  
  TClonesArray* AODMCTrackArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray == NULL) return;
  if (fUseCouts) std::cout << "AODMCTrackArray->GetEntriesFast() is " << AODMCTrackArray->GetEntriesFast() << std::endl; //AOD Number of gen particles - used for looping over tracks

    // Loop over all primary MC particle
    for(Long_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {

      AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
      if (!particle) {
        if (fUseCouts) std::cout << "couldn't find particle" << std::endl;
        continue;
      }

      if(!particle->IsPhysicalPrimary()) {
        if (fUseCouts) std::cout << "particle wasnt primary" << std::endl;
        continue;
      }

      Float_t fPMC_Truth = particle->P();
      Float_t fPtMC_Truth = particle->Pt();
      Float_t fEtaMC_Truth = particle->Eta();
      Float_t fPhiMC_Truth = particle->Phi();
      
      if (TMath::Abs(particle->Y()) > fYCut && fDoRapCut) {
        if (fUseCouts) std::cout << "particle didnt pass y cut" << std::endl;
        continue;
      }

      if (TMath::Abs(particle->Eta()) > fEtaCut) {
        if (fUseCouts) std::cout << "particle didn't pass eta cut" << std::endl;
        continue;
      }

      if (particle->Pt() < 0.15 || particle->Pt() > 100.0) {
        if (fUseCouts) std::cout << "particle didn't pass pT cut" << std::endl;
        continue;
      }

      Int_t TfPdgcode = particle->GetPdgCode();
      int fPdgIndex = -1;
      if ((TMath::Abs(TfPdgcode) == 211))
        fPdgIndex = 0; //Particle is a Pion
      else if ((TMath::Abs(TfPdgcode) == 321))
        fPdgIndex = 1; //Particle is a Kaon
      else if ((TMath::Abs(TfPdgcode) == 2212))
        fPdgIndex = 2; //Particle is a Proton

      Bool_t fSignMC = kFALSE;
      if (TfPdgcode > 0)
        fSignMC = kFALSE; //Particle is positive
      else
        fSignMC = kTRUE; //Particle is negative

      int iMCPartLabel = TMath::Abs(particle->GetLabel());
      bool particleIsPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iMCPartLabel,MCEvent());
      if (particleIsPileup) {
        if (fUseCouts) std::cout << "This MC particle was from OOB pileup" << std::endl;
        continue;
      }

      //Save truth level by PID
      if (fPdgIndex == 0) fHistMCTruth_TrackEff_Den_pi->Fill(fPtMC_Truth, TMath::Abs(fEtaMC_Truth));
      if (fPdgIndex == 1) fHistMCTruth_TrackEff_Den_ka->Fill(fPtMC_Truth, TMath::Abs(fEtaMC_Truth));
      if (fPdgIndex == 2) fHistMCTruth_TrackEff_Den_pr->Fill(fPtMC_Truth, TMath::Abs(fEtaMC_Truth));
    }

  //
  // ======================================================================
  // ------   reconstructed MC particles   -------------
  // ======================================================================
  //

  Int_t trackOrigin = -10;

  for(Int_t irectrack = 0; irectrack < fAOD->GetNumberOfTracks(); irectrack++)
  {
    AliAODTrack* trackReal = static_cast<AliAODTrack*>(fAOD->GetTrack(irectrack));
    if (trackReal==NULL) {
      if (fUseCouts) std::cout << "Didn't have a reco track" << std::endl;
      continue;
    }

    if(!trackReal || !trackReal->TestFilterBit(fAOD_FilterBits)) {
      continue;
    }

    //
    // Get generated track info
    Int_t lab = TMath::Abs(trackReal->GetLabel());
    Int_t imc=trackReal->GetLabel();

    const AliAODMCParticle *mc=static_cast<AliAODMCParticle*>(AODMCTrackArray->At(lab));
    Int_t pdg = mc->GetPdgCode();

    if(mc->IsPhysicalPrimary()) {
      trackOrigin = 0;
    }
    if(mc->IsSecondaryFromWeakDecay()) {
      trackOrigin = 1;
    }
    if(mc->IsSecondaryFromMaterial()) {
      trackOrigin = 2;
    }

    Float_t DCAxy = -999.;
    Float_t DCAz = -999.;

    trackReal->GetImpactParameters(DCAxy, DCAz);

    if (fUseCouts) cout << "DCAxy for trackOrigin= " << trackOrigin << " is " << DCAxy << endl;

    Int_t TOFTrkLabel[3] = { -1 }; //This can contain three particles wich occupy the same cluster

    Int_t fMCTOFMatch = -2;
    trackReal->GetTOFLabel(TOFTrkLabel);

    if (TOFTrkLabel[0] == -1) {// Track was not matched to any TOF hit.
      fMCTOFMatch = -1;
    }
    else if (lab == TOFTrkLabel[0]) {// Track was correctly matched to a TOF hit.
      fMCTOFMatch = 0;
    }
    else {// Track was matched to a TOF hit but comes from mismatch!
      fMCTOFMatch = 1;
    }
    
    // match the track with mc track
    Int_t iPart = -1;
    if (TMath::Abs(pdg) == kPDGel) { iPart = 0; } // select el
    if (TMath::Abs(pdg) == kPDGpi) { iPart = 1; } // select pi
    if (TMath::Abs(pdg) == kPDGka) { iPart = 2; } // select ka
    if (TMath::Abs(pdg) == kPDGpr) { iPart = 3; } // select pr
    if (TMath::Abs(pdg) == kPDGde) { iPart = 4; } // select de
    if (TMath::Abs(pdg) == kPDGmu) { iPart = 5; } // select mu

    //
    GetExpecteds(trackReal);
    SetCutBitsAndSomeTrackVariables(trackReal);
    //
    //if (trackReal-> GetInnerParam()){ AODs don't have this info
      fPtotMC       = trackReal->P(); //trackReal->GetInnerParam()->GetP();
      fTPCSignalMC  = trackReal->GetTPCsignal();
    //}
    fEtaMC        = trackReal->Eta();
    Float_t fYMC          = trackReal->Y();
    fPtMC         = trackReal->Pt();
    fSignMC       = trackReal->GetSign();
    Float_t pMC   = trackReal->P();
    Float_t fPhiMC= trackReal->Phi();

    Float_t fRapidityMC = TMath::ASinH(fPt / TMath::Sqrt(AliPID::ParticleMass(iPart + 1) * AliPID::ParticleMass(iPart + 1) + fPt * fPt) * TMath::SinH(fEta));

    if (TMath::Abs(fRapidityMC) > fYCut && fDoRapCut) {
      if (fUseCouts) std::cout << "particle didn't pass rapidity cut with MC PID" << std::endl;
      continue;
    }

    if (fPtMC<0.15 || fPtMC>100.0) continue; //So we can match the jets that we throw out w/ max track pT>100

    if (TMath::Abs(fEtaMC) > fEtaCut) continue;

    if (iPart == 1 && trackOrigin==0) fHistMCReco_pi_prim_DCAxy->Fill(fPtMC, DCAxy);
    if (iPart == 1 && trackOrigin==1) fHistMCReco_pi_scdweak_DCAxy->Fill(fPtMC, DCAxy);
    if (iPart == 3 && trackOrigin==0) fHistMCReco_pr_prim_DCAxy->Fill(fPtMC, DCAxy);
    if (iPart == 3 && trackOrigin==1) fHistMCReco_pr_scdweak_DCAxy->Fill(fPtMC, DCAxy);
    if (iPart == 3 && trackOrigin==2) fHistMCReco_pr_scdmat_DCAxy->Fill(fPtMC, DCAxy);

    if(!mc->IsPhysicalPrimary()) { //CHeck that this is the same as !fMCEvent->IsPhysicalPrimary(lab)
      if (fUseCouts) std::cout << "particle wasnt primary" << std::endl; //most aren't primaries
      continue;
    }

    if (iPart == 1) fHistMCReco_TrackEff_Num_pi->Fill(fPtMC, TMath::Abs(fEtaMC));
    if (iPart == 2) fHistMCReco_TrackEff_Num_ka->Fill(fPtMC, TMath::Abs(fEtaMC));
    if (iPart == 3) fHistMCReco_TrackEff_Num_pr->Fill(fPtMC, TMath::Abs(fEtaMC));

    if (iPart == 1 && fMCTOFMatch == 0) fHistMCReco_MatchEff_Num_pi->Fill(fPtMC, TMath::Abs(fEtaMC));
    if (iPart == 2 && fMCTOFMatch == 0) fHistMCReco_MatchEff_Num_ka->Fill(fPtMC, TMath::Abs(fEtaMC));
    if (iPart == 3 && fMCTOFMatch == 0) fHistMCReco_MatchEff_Num_pr->Fill(fPtMC, TMath::Abs(fEtaMC));

  } // ======= end of track loop for MC dEdx filling =======

}
//________________________________________________________________________
void AliAnalysisTaskJetHadroAOD::FillMCJets()
{
  //Track eff = Num1/Den
  //TOF Match Eff = Num2/Num1

  //Den = # MC Truth particles that pass: 1. Is Primary 2. eta cut 3. pT cut 4. rap cut (if applicable) 5. Event selections (same as reco) 6. MC Truth PDG code
  //Num1 = # MC Reco particles that pass: 1. Is Primary 2. All data track selection cuts 3. eta cut 3. pT cut 4. rap cut (if applicable) 5. Event selections 6. MC Truth PDG code
  //Num2 = # MC Reco particles that pass: all of Num1 + 2. TOFMatch==0

  //To implement jet matching, we want: 1. To see if any MC jet matches w/in some delta R and delta pT threshold of the given reco jet. 2. To find the closest (in R, R and pT?) MC jet to the reco jet
  //3. Only count an MC jet if it is in the list of closest jets to reco jets 4. Only count a reco jet if it has a closest MC jet (aka is in the list).
  //
  AliTOFPIDResponse fTOFPIDResponse = fPIDResponse->GetTOFResponse();

  //Make jet wrappers
  float fTrackPt = 0.15;
  float fTrackEta = fEtaCut;
  float fGhostArea = 0.005;
  float bgJetAbsEtaCut = 0.7;
  float bgJetRadius = 0.2;

  Float_t jetAbsEtaCut = fTrackEta-0.4;

  //SOME CODE FROM NIMA
  fFastJetWrapper_Gen->Clear();
  fFastJetWrapper_Gen->SetR(0.4);
  fFastJetWrapper_Gen->SetAlgorithm(fastjet::JetAlgorithm::antikt_algorithm);
  fFastJetWrapper_Gen->SetRecombScheme(fastjet::RecombinationScheme::pt_scheme);
  //fFastJetWrapper_Gen->SetStrategy(fastjet::Strategy::Best); //this doesn't need to be here - in FJ constructor already
  fFastJetWrapper_Gen->SetGhostArea(fGhostArea);
  fFastJetWrapper_Gen->SetAreaType(fastjet::AreaType::active_area);
  fFastJetWrapper_Gen->SetMaxRap(1.0);

  fFastJetWrapperBG_Gen->Clear();
  fFastJetWrapperBG_Gen->SetR(bgJetRadius);
  fFastJetWrapperBG_Gen->SetAlgorithm(fastjet::JetAlgorithm::kt_algorithm);
  fFastJetWrapperBG_Gen->SetRecombScheme(fastjet::RecombinationScheme::pt_scheme);
  //fFastJetWrapperBG_Gen->SetStrategy(fastjet::Strategy::Best); //this doesn't need to be here - in FJ constructor already
  fFastJetWrapperBG_Gen->SetGhostArea(fGhostArea);
  fFastJetWrapperBG_Gen->SetAreaType(fastjet::AreaType::active_area_explicit_ghosts);
  fFastJetWrapperBG_Gen->SetMaxRap(1.0);

  fFastJetWrapper_Rec->Clear();
  fFastJetWrapper_Rec->SetR(0.4);
  fFastJetWrapper_Rec->SetAlgorithm(fastjet::JetAlgorithm::antikt_algorithm);
  fFastJetWrapper_Rec->SetRecombScheme(fastjet::RecombinationScheme::pt_scheme);
  //fFastJetWrapper->SetStrategy(fastjet::Strategy::Best); //this doesn't need to be here - in FJ constructor already
  fFastJetWrapper_Rec->SetGhostArea(fGhostArea);
  fFastJetWrapper_Rec->SetAreaType(fastjet::AreaType::active_area);
  fFastJetWrapper_Rec->SetMaxRap(1.0);

  fFastJetWrapperBG_Rec->Clear();
  fFastJetWrapperBG_Rec->SetR(bgJetRadius);
  fFastJetWrapperBG_Rec->SetAlgorithm(fastjet::JetAlgorithm::kt_algorithm);
  fFastJetWrapperBG_Rec->SetRecombScheme(fastjet::RecombinationScheme::pt_scheme);
  //fFastJetWrapperBG->SetStrategy(fastjet::Strategy::Best); //this doesn't need to be here - in FJ constructor already
  fFastJetWrapperBG_Rec->SetGhostArea(fGhostArea);
  fFastJetWrapperBG_Rec->SetAreaType(fastjet::AreaType::active_area_explicit_ghosts);
  fFastJetWrapperBG_Rec->SetMaxRap(1.0);

  float particleEtaCut = fEtaCut;

  //Int_t numMCTracks = fMCEvent->GetNumberOfTracks();      //ESD Number of gen particles - used for looping over tracks
  
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FillMCJets ===== " << std::endl;

  //Get truth level particles
  TClonesArray* AODMCTrackArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray == NULL) return;
  if (fUseCouts) std::cout << "AODMCTrackArray->GetEntriesFast() is " << AODMCTrackArray->GetEntriesFast() << std::endl; //AOD Number of gen particles - used for looping over tracks

  // Loop over all primary MC particle
  for(Long_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {

    AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
    if (!particle) {
      if (fUseCouts) std::cout << "couldn't find particle" << std::endl;
      continue;
    }

    if(!particle->IsPhysicalPrimary()) {
      if (fUseCouts) std::cout << "particle wasnt primary" << std::endl;
      continue;
    }

    Float_t fPMC_Truth = particle->P();
    Float_t fPtMC_Truth = particle->Pt();
    Float_t fEtaMC_Truth = particle->Eta();
    Float_t fPhiMC_Truth = particle->Phi();
    
    if (TMath::Abs(particle->Y()) > fYCut && fDoRapCut) {
      if (fUseCouts) std::cout << "particle didnt pass y cut" << std::endl;
      continue;
    }

    if (TMath::Abs(particle->Eta()) > fEtaCut) {
      if (fUseCouts) std::cout << "particle didn't pass eta cut" << std::endl;
      continue;
    }

    if (particle->Pt() < 0.15) { //|| particle->Pt() > 100.0 needs to be after jet finding
      if (fUseCouts) std::cout << "particle didn't pass pT cut" << std::endl;
      continue;
    }

    int iMCPartLabel = TMath::Abs(particle->GetLabel());
    bool particleIsPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iMCPartLabel,MCEvent());
    if (particleIsPileup) {
      if (fUseCouts) std::cout << "This jet MC particle was from OOB pileup" << std::endl;
      continue;
    }

    //Make truth level jets
    // loop over AOD tracks and add their four vector to wrapper --> identical to track container in EMC jet
    fFastJetWrapper_Gen->AddInputVector(particle->Px(), particle->Py(), particle->Pz(), particle->E(), i);//TMath::Sqrt(track->P()*track->P()+0.13957*0.13957),iTrack);
    fFastJetWrapperBG_Gen->AddInputVector(particle->Px(), particle->Py(), particle->Pz(), particle->E(), i);//TMath::Sqrt(track->P()*track->P()+0.13957*0.13957),iTrack);
  }

  //
  // background jet definitions
  fastjet::JetMedianBackgroundEstimator bgE;
  fastjet::Selector selectorBG = !fastjet::SelectorNHardest(2); //set the max eta cut on the estimator, then get rid of 2 highest pt jets
  bgE.set_selector(selectorBG);

  fastjet::Selector selectorBGjets = !fastjet::SelectorIsPureGhost() * fastjet::SelectorAbsEtaMax(bgJetAbsEtaCut) * fastjet::SelectorPtRange(fTrackPt, 1000.0);
  fFastJetWrapperBG_Gen->Run();
  std::vector<fastjet::PseudoJet> jetsBG = sorted_by_pt(selectorBGjets(fFastJetWrapperBG_Gen->GetInclusiveJets()));

  Float_t rho_MC = 0.0;
  if (jetsBG.size() > 0) {
    bgE.set_jets(jetsBG);  // give the kT jets to the background estimator
    rho_MC = bgE.rho();
  }
  if (fUseCouts) std::cout << "rho_MC is " << rho_MC << std::endl;

  // run jet finder using wrapper
  fFastJetWrapper_Gen->Run();
  std::vector<fastjet::PseudoJet> jets_truth = sorted_by_pt(fFastJetWrapper_Gen->GetInclusiveJets());


  //Get Reco level particles
  Int_t trackOrigin = -10;

  int num_reco_jets_w_multiple_matches = 0;
  int num_reco_jets_w_matches = 0;
  for (Int_t iTrack = 0; iTrack < fAOD->GetNumberOfTracks(); iTrack++) {
    AliAODTrack* trackReal = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if (trackReal==NULL) {
      if (fUseCouts) std::cout << "Didn't have a reco track" << std::endl;
      continue;
    }

    Int_t lab = TMath::Abs(trackReal->GetLabel());
    const AliAODMCParticle *mc=static_cast<AliAODMCParticle*>(AODMCTrackArray->At(lab));

    if(!mc->IsPhysicalPrimary()) {
      if (fUseCouts) std::cout << "particle wasnt primary" << std::endl;
      continue;
    }

    if(!trackReal->TestFilterBit(fAOD_FilterBits)) continue; //hybrid track cuts filterbit
    if (trackReal->Pt() < fTrackPt || TMath::Abs(trackReal->Eta()) > particleEtaCut) continue;
    if (TMath::Abs(trackReal->Y()) > fYCut && fDoRapCut) continue;
    fFastJetWrapper_Rec->AddInputVector(trackReal->Px(), trackReal->Py(), trackReal->Pz(), trackReal->E(), iTrack);//TMath::Sqrt(track->P()*track->P()+0.13957*0.13957),iTrack);
    fFastJetWrapperBG_Rec->AddInputVector(trackReal->Px(), trackReal->Py(), trackReal->Pz(), trackReal->E(), iTrack);//TMath::Sqrt(track->P()*track->P()+0.13957*0.13957),iTrack);
  }
  //
  // background jet definitions
  fastjet::JetMedianBackgroundEstimator bgERec;
  fastjet::Selector selectorBGRec = !fastjet::SelectorNHardest(2); //set the max eta cut on the estimator, then get rid of 2 highest pt jets
  bgERec.set_selector(selectorBGRec);

  fastjet::Selector selectorBGjetsRec = !fastjet::SelectorIsPureGhost() * fastjet::SelectorAbsEtaMax(bgJetAbsEtaCut) * fastjet::SelectorPtRange(fTrackPt, 1000.0);
  fFastJetWrapperBG_Rec->Run();
  std::vector<fastjet::PseudoJet> jetsBGRec = sorted_by_pt(selectorBGjetsRec(fFastJetWrapperBG_Rec->GetInclusiveJets()));

  Float_t rho_MC_Rec = 0.0;
  if (jetsBGRec.size() > 0) {
    bgERec.set_jets(jetsBGRec);  // give the kT jets to the background estimator
    rho_MC_Rec = bgERec.rho();
  }
  if (fUseCouts) std::cout << "rho_MC_Rec is " << rho_MC_Rec << std::endl;

  fFastJetWrapper_Rec->Run();
  std::vector<fastjet::PseudoJet> jets_rec = sorted_by_pt(fFastJetWrapper_Rec->GetInclusiveJets());

  vector<pair<int, int>> matched_MC_jets_indexes;

  Float_t leadJetPhi_Rec = 999;
  Float_t leadJetEta_Rec = 999;
  Float_t leadJetPtSub_Rec = -999;

  Float_t leadJetPhi_Gen = 999;
  Float_t leadJetEta_Gen = 999;
  Float_t leadJetPtSub_Gen = -999;

  // start of reco jet loop
  for (Int_t ijet_reco=0; ijet_reco<Int_t(jets_rec.size()); ijet_reco++){
    fastjet::PseudoJet jet_reco = jets_rec[ijet_reco];
    if (jet_reco.pt() < fTrackPt || jet_reco.perp() > 1000.0 || TMath::Abs(jet_reco.eta()) >= jetAbsEtaCut) continue;
    Float_t jet_recopt = jet_reco.pt();
    Float_t jet_recophi = jet_reco.phi();
    Float_t jet_recoeta = jet_reco.eta();
    Float_t jet_recoArea = jet_reco.area();
    Float_t jet_recoptsub = jet_recopt - rho_MC_Rec*jet_recoArea;
    Int_t jet_recoNum = jets_rec.size();
    std::vector<fastjet::PseudoJet> constituents = sorted_by_pt(jets_rec[ijet_reco].constituents());

    Int_t nConstituents = constituents.size();
    if (nConstituents<1) continue;
    fastjet::PseudoJet highestpt_const = constituents[0];
    if (highestpt_const.pt() > 100.0) continue;

    fHistMCR_Jet_ptsub_v_area->Fill(jet_recoArea, jet_recoptsub);

    if (jet_recoArea <= fjetMinArea) continue;

    //could add jet pT > 30 cut here if needed to reduce processing time (min for share pT > 0.5 for 60 GeV jets)

    int i_best_match = -1;
    double min_share_pT = 0.5;
    bool found_match = false;
    double matched_gen_pTsub = -1000.;

    //match to truth jets
    for (Int_t ijet_truth=0; ijet_truth<Int_t(jets_truth.size()); ijet_truth++){
      fastjet::PseudoJet jet_MC = jets_truth[ijet_truth];
      if (jet_MC.pt() < fTrackPt || jet_MC.perp() > 1000.0 || TMath::Abs(jet_MC.eta()) >= jetAbsEtaCut) continue;
      std::vector<fastjet::PseudoJet> constituents_MC = sorted_by_pt(jets_truth[ijet_truth].constituents());
      Int_t nConstituents_MC = constituents_MC.size();
      if (nConstituents_MC<1) continue;
      fastjet::PseudoJet highestpt_const_MC = constituents_MC[0];
      if (highestpt_const_MC.pt() > 100.0) continue;
      Float_t jetArea_MC = jet_MC.area();
      Float_t jetpt_MC = jet_MC.pt();
      Float_t jetptsub_MC = jetpt_MC - rho_MC*jetArea_MC;
      if (jetArea_MC <= fjetMinArea) continue;

      //could add jet pT > 30 cut here if needed to reduce processing time (min for share pT > 0.5 for 60 GeV jets)
        
      double dR = jet_reco.delta_R(jet_MC);
      double pT_share_frac = (jetptsub_MC - jet_recoptsub)/jetptsub_MC;
      fHistMC_Jet_shared_pt_frac->Fill(jetptsub_MC,pT_share_frac,jet_recoptsub);
      fHistMC_Jet_deltaR->Fill(jetptsub_MC,dR,jet_recoptsub);

      if (dR > 0.3) continue;
      if (TMath::Abs(pT_share_frac) < 0.5 && min_share_pT!=0.5) num_reco_jets_w_multiple_matches++; //we have more than one match candidate
      if (TMath::Abs(pT_share_frac) < min_share_pT) {
        if (min_share_pT==0.5) num_reco_jets_w_matches++;
        min_share_pT = TMath::Abs(pT_share_frac);
        i_best_match = ijet_truth;
        found_match = true;
        matched_gen_pTsub = jetptsub_MC;
      }
    }

    if (found_match) {
        matched_MC_jets_indexes.push_back(std::make_pair(i_best_match,ijet_reco));
    }

    if (jet_recoptsub > leadJetPtSub_Rec){
      leadJetPtSub_Rec = jet_recoptsub;
      leadJetEta_Rec = jet_recoeta;
      leadJetPhi_Rec = jet_recophi;
    }


    for(Int_t i = 0; i < nConstituents; i++)
    {
      fastjet::PseudoJet &constituent = constituents[i];
      Float_t pt = constituent.pt();
      if (pt<0.15) continue;  

      Int_t trackIndex = constituent.user_index();
      AliAODTrack* trackReal = static_cast<AliAODTrack*>(fAOD->GetTrack(trackIndex));

      if (trackReal==NULL) {
        if (fUseCouts) std::cout << "Didn't have a reco track" << std::endl;
        continue;
      }
      //
      // Get generated track info
      Int_t lab = TMath::Abs(trackReal->GetLabel());
      Int_t imc=trackReal->GetLabel();

      const AliAODMCParticle *mc=static_cast<AliAODMCParticle*>(AODMCTrackArray->At(lab));
      Int_t pdg = mc->GetPdgCode();

      Int_t TOFTrkLabel[3] = { -1 };

      Int_t fMCTOFMatch = -2;
      trackReal->GetTOFLabel(TOFTrkLabel);
      if (TOFTrkLabel[0] == -1) {// Track was not matched to any TOF hit.
        fMCTOFMatch = -1;
      }
      else if (lab == TOFTrkLabel[0]) {// Track was correctly matched to a TOF hit.
        fMCTOFMatch = 0;
      }
      else {// Track was matched to a TOF hit but comes from mismatch!
        fMCTOFMatch = 1;
      }
       
      // check the origin of the track
      trackOrigin = -10;

      if(!trackReal || !trackReal->TestFilterBit(fAOD_FilterBits)) {
        continue;
      }
      //
      //
      // match the track with mc track
      Int_t iPart = -1; //-10;
      if (TMath::Abs(pdg) == kPDGel) { iPart = 0; } // select el
      if (TMath::Abs(pdg) == kPDGpi) { iPart = 1; } // select pi
      if (TMath::Abs(pdg) == kPDGka) { iPart = 2; } // select ka
      if (TMath::Abs(pdg) == kPDGpr) { iPart = 3; } // select pr
      if (TMath::Abs(pdg) == kPDGde) { iPart = 4; } // select de
      if (TMath::Abs(pdg) == kPDGmu) { iPart = 5; } // select mu
      //
      GetExpecteds(trackReal);
      SetCutBitsAndSomeTrackVariables(trackReal);
      //
      //if (trackReal-> GetInnerParam()){ AODs don't have this info
        fPtotMC       = trackReal->P(); //trackReal->GetInnerParam()->GetP();
        fTPCSignalMC  = trackReal->GetTPCsignal();
      //}
      fEtaMC        = trackReal->Eta();
      Float_t fYMC          = trackReal->Y();
      fPtMC         = trackReal->Pt();
      fSignMC       = trackReal->GetSign();
      Float_t pMC   = trackReal->P();
      Float_t fPhiMC= trackReal->Phi();

      Float_t fRapidityMC = TMath::ASinH(fPt / TMath::Sqrt(AliPID::ParticleMass(iPart + 1) * AliPID::ParticleMass(iPart + 1) + fPt * fPt) * TMath::SinH(fEta));

      if (TMath::Abs(fRapidityMC) > fYCut && fDoRapCut) {
        if (fUseCouts) std::cout << "didn't pass rap cut" << std::endl;
        continue;
      }
      //
      //
      //

      if (fPtMC<0.15 || fPtMC>100.0) continue; //So we can match the jets that we throw out w/ max track pT>100

      if (TMath::Abs(fEtaMC) > fEtaCut) continue;

      if (found_match==true && iPart == 1) fHist_JetMatch_MCReco_TrackEff_Num_pi->Fill(fPtMC, jet_recoptsub, matched_gen_pTsub);
      if (found_match==true && iPart == 2) fHist_JetMatch_MCReco_TrackEff_Num_ka->Fill(fPtMC, jet_recoptsub, matched_gen_pTsub);
      if (found_match==true && iPart == 3) fHist_JetMatch_MCReco_TrackEff_Num_pr->Fill(fPtMC, jet_recoptsub, matched_gen_pTsub);

      if (found_match==true && iPart == 1 && fMCTOFMatch == 0) fHist_JetMatch_MCReco_MatchEff_Num_pi->Fill(fPtMC, jet_recoptsub, matched_gen_pTsub);
      if (found_match==true && iPart == 2 && fMCTOFMatch == 0) fHist_JetMatch_MCReco_MatchEff_Num_ka->Fill(fPtMC, jet_recoptsub, matched_gen_pTsub);
      if (found_match==true && iPart == 3 && fMCTOFMatch == 0) fHist_JetMatch_MCReco_MatchEff_Num_pr->Fill(fPtMC, jet_recoptsub, matched_gen_pTsub);

      if (jet_recoptsub <= fjetMinPtSub || jet_recoptsub >= fjetMaxPtSub) continue;

      if (iPart == 1) fHist_Jet_MCReco_TrackEff_Num_pi->Fill(fPtMC, TMath::Abs(fEtaMC));
      if (iPart == 2) fHist_Jet_MCReco_TrackEff_Num_ka->Fill(fPtMC, TMath::Abs(fEtaMC));
      if (iPart == 3) fHist_Jet_MCReco_TrackEff_Num_pr->Fill(fPtMC, TMath::Abs(fEtaMC));

      if (iPart == 1 && fMCTOFMatch == 0) fHist_Jet_MCReco_MatchEff_Num_pi->Fill(fPtMC, TMath::Abs(fEtaMC));
      if (iPart == 2 && fMCTOFMatch == 0) fHist_Jet_MCReco_MatchEff_Num_ka->Fill(fPtMC, TMath::Abs(fEtaMC));
      if (iPart == 3 && fMCTOFMatch == 0) fHist_Jet_MCReco_MatchEff_Num_pr->Fill(fPtMC, TMath::Abs(fEtaMC));


    } // ======= end of track loop for MC jet reco filling =======

  }

  if ( (fDoPerpCone || fDoRandCone) && (leadJetPtSub_Rec > fjetMinPtSub) ) {
    //now we find the eta,phi perp (and R.C.) to it
    //From the AliAnalysisTaskRhoPerpCone Task
    Double_t dPhi1 = 999.;
    Double_t dPhi2 = 999.;
    Double_t dEta = 999.;
    Double_t Axis1 = 999, Axis2 = 999;

    if (fDoRandCone){
      double OffsetRandom = (1 + fRandom->Rndm()) * TMath::Pi()/3; // creating random Phi in the range of PI/3 and 2PI/3
      Axis1 = leadJetPhi_Rec + OffsetRandom;  // adding the random Phi to the Phi of leading jet, to get the Phi of random cone
      Axis2 = leadJetPhi_Rec - OffsetRandom; // mirroring it, so you get Phi of the 2nd random cone
      if(Axis1 > TMath::TwoPi()) Axis1 -= TMath::TwoPi(); // checking if it's larger than 2PI (leading jet phi is distributed in range of 0 and 2PI
      if(Axis2 < 0) Axis2 += TMath::TwoPi(); // checking if 2nd cone is smaller than 0
    }

    if (fDoPerpCone){
      Axis1 = ((leadJetPhi_Rec + (TMath::Pi() / 2.)) > TMath::TwoPi()) ? leadJetPhi_Rec - ((3. / 2.) * TMath::Pi()) : leadJetPhi_Rec + (TMath::Pi() / 2.);
      Axis2 = ((leadJetPhi_Rec - (TMath::Pi() / 2.)) < 0) ? leadJetPhi_Rec + ((3. / 2.) * TMath::Pi()) : leadJetPhi_Rec - (TMath::Pi() / 2.);
    }

    Int_t trackOrigin = -10;

    for(Int_t irectrack = 0; irectrack < fAOD->GetNumberOfTracks(); irectrack++)
    {
      AliAODTrack* trackReal = static_cast<AliAODTrack*>(fAOD->GetTrack(irectrack));
      if (trackReal==NULL) {
        if (fUseCouts) std::cout << "Didn't have a reco track" << std::endl;
        continue;
      }

      if(!trackReal || !trackReal->TestFilterBit(fAOD_FilterBits)) {
        continue;
      }

      fEtaMC        = trackReal->Eta();
      if (TMath::Abs(fEta) > fEtaCut) continue;
                    
      Float_t mod_track_phi = trackReal->Phi() + TMath::Pi();
      //Check if the track is within the R=0.4 cone in eta, phi
      dPhi1 = TMath::Abs(mod_track_phi - Axis1);
      dPhi1 = (dPhi1 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi1 : dPhi1;
      dPhi2 = TMath::Abs(mod_track_phi - Axis2);
      dPhi2 = (dPhi2 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi2 : dPhi2;
      dEta = leadJetEta_Rec - trackReal->Eta();

      if ((TMath::Sqrt(dPhi1 * dPhi1 + dEta * dEta) > 0.4) && (TMath::Sqrt(dPhi2 * dPhi2 + dEta * dEta) > 0.4)) continue; //scale the yields by 1/(2Ncones*piR^2) offline

      //
      // Get generated track info
      Int_t lab = TMath::Abs(trackReal->GetLabel());
      Int_t imc=trackReal->GetLabel();

      const AliAODMCParticle *mc=static_cast<AliAODMCParticle*>(AODMCTrackArray->At(lab));
      Int_t pdg = mc->GetPdgCode();

      if(mc->IsPhysicalPrimary()) {
        trackOrigin = 0;
      }
      if(mc->IsSecondaryFromWeakDecay()) {
        trackOrigin = 1;
      }
      if(mc->IsSecondaryFromMaterial()) {
        trackOrigin = 2;
      }

      Float_t DCAxy = -999.;
      Float_t DCAz = -999.;

      trackReal->GetImpactParameters(DCAxy, DCAz);

      if (fUseCouts) cout << "DCAxy for trackOrigin= " << trackOrigin << " is " << DCAxy << endl;

      Int_t TOFTrkLabel[3] = { -1 }; //This can contain three particles wich occupy the same cluster

      Int_t fMCTOFMatch = -2;
      trackReal->GetTOFLabel(TOFTrkLabel);

      if (TOFTrkLabel[0] == -1) {// Track was not matched to any TOF hit.
        fMCTOFMatch = -1;
      }
      else if (lab == TOFTrkLabel[0]) {// Track was correctly matched to a TOF hit.
        fMCTOFMatch = 0;
      }
      else {// Track was matched to a TOF hit but comes from mismatch!
        fMCTOFMatch = 1;
      }
      
      // match the track with mc track
      Int_t iPart = -1;
      if (TMath::Abs(pdg) == kPDGel) { iPart = 0; } // select el
      if (TMath::Abs(pdg) == kPDGpi) { iPart = 1; } // select pi
      if (TMath::Abs(pdg) == kPDGka) { iPart = 2; } // select ka
      if (TMath::Abs(pdg) == kPDGpr) { iPart = 3; } // select pr
      if (TMath::Abs(pdg) == kPDGde) { iPart = 4; } // select de
      if (TMath::Abs(pdg) == kPDGmu) { iPart = 5; } // select mu

      //
      GetExpecteds(trackReal);
      SetCutBitsAndSomeTrackVariables(trackReal);
      //
      //if (trackReal-> GetInnerParam()){ AODs don't have this info
        fPtotMC       = trackReal->P(); //trackReal->GetInnerParam()->GetP();
        fTPCSignalMC  = trackReal->GetTPCsignal();
      //}
      fEtaMC        = trackReal->Eta();
      Float_t fYMC          = trackReal->Y();
      fPtMC         = trackReal->Pt();
      fSignMC       = trackReal->GetSign();
      Float_t pMC   = trackReal->P();
      Float_t fPhiMC= trackReal->Phi();

      Float_t fRapidityMC = TMath::ASinH(fPt / TMath::Sqrt(AliPID::ParticleMass(iPart + 1) * AliPID::ParticleMass(iPart + 1) + fPt * fPt) * TMath::SinH(fEta));

      if (TMath::Abs(fRapidityMC) > fYCut && fDoRapCut) {
        if (fUseCouts) std::cout << "particle didn't pass rapidity cut with MC PID" << std::endl;
        continue;
      }

      if (fPtMC<0.15 || fPtMC>100.0) continue; //So we can match the jets that we throw out w/ max track pT>100

      if (TMath::Abs(fEtaMC) > fEtaCut) continue;

      if (iPart == 1 && trackOrigin==0) fHistMCReco_pi_prim_DCAxy->Fill(fPtMC, DCAxy);
      if (iPart == 1 && trackOrigin==1) fHistMCReco_pi_scdweak_DCAxy->Fill(fPtMC, DCAxy);
      if (iPart == 3 && trackOrigin==0) fHistMCReco_pr_prim_DCAxy->Fill(fPtMC, DCAxy);
      if (iPart == 3 && trackOrigin==1) fHistMCReco_pr_scdweak_DCAxy->Fill(fPtMC, DCAxy);
      if (iPart == 3 && trackOrigin==2) fHistMCReco_pr_scdmat_DCAxy->Fill(fPtMC, DCAxy);

      if(!mc->IsPhysicalPrimary()) { //CHeck that this is the same as !fMCEvent->IsPhysicalPrimary(lab)
        if (fUseCouts) std::cout << "particle wasnt primary" << std::endl; //most aren't primaries
        continue;
      }

      if (iPart == 1) fHistMCReco_TrackEff_Num_pi->Fill(fPtMC, TMath::Abs(fEtaMC));
      if (iPart == 2) fHistMCReco_TrackEff_Num_ka->Fill(fPtMC, TMath::Abs(fEtaMC));
      if (iPart == 3) fHistMCReco_TrackEff_Num_pr->Fill(fPtMC, TMath::Abs(fEtaMC));

      if (iPart == 1 && fMCTOFMatch == 0) fHistMCReco_MatchEff_Num_pi->Fill(fPtMC, TMath::Abs(fEtaMC));
      if (iPart == 2 && fMCTOFMatch == 0) fHistMCReco_MatchEff_Num_ka->Fill(fPtMC, TMath::Abs(fEtaMC));
      if (iPart == 3 && fMCTOFMatch == 0) fHistMCReco_MatchEff_Num_pr->Fill(fPtMC, TMath::Abs(fEtaMC));

    } // end inc tracks loop for cone
    fFilledUECone_Rec = kTRUE;
  } //end cone loop


  // start of truth jet loop
  for (Int_t ijet_truth=0; ijet_truth<Int_t(jets_truth.size()); ijet_truth++){
    fastjet::PseudoJet jet = jets_truth[ijet_truth];
    if (jet.pt() < fTrackPt || jet.perp() > 1000.0 || TMath::Abs(jet.eta()) >= jetAbsEtaCut) continue;
    Float_t jetpt = jet.pt();
    Float_t jetphi = jet.phi();
    Float_t jeteta = jet.eta();
    Float_t jetArea = jet.area();
    Float_t jetptsub = jetpt - rho_MC*jetArea;
    Int_t jetNum = jets_truth.size();
    std::vector<fastjet::PseudoJet> constituents = sorted_by_pt(jets_truth[ijet_truth].constituents());

    Int_t nConstituents = constituents.size();
    if (nConstituents<1) continue;
    fastjet::PseudoJet highestpt_const = constituents[0];
    if (highestpt_const.pt() > 100.0) continue;

    fHistMCT_Jet_ptsub_v_area->Fill(jetArea, jetptsub);

    if (jetArea <= fjetMinArea) continue;

    //could add jet pT > 30 cut here if needed to reduce processing time (min for share pT > 0.5 for 60 GeV jets)

    bool found = false;
    double matched_reco_jet_i = -1;

    //Check if this is a matched truth jet
    for (const auto& pair : matched_MC_jets_indexes) {
        if (pair.first == ijet_truth) {
            found = true;
            matched_reco_jet_i = pair.second;
            break; // Exit the loop early if the match is found
        }
    }


    if (jetptsub > leadJetPtSub_Gen){
      leadJetPtSub_Gen = jetptsub;
      leadJetEta_Gen = jeteta;
      leadJetPhi_Gen = jetphi;
    }

    for(Int_t i = 0; i < nConstituents; i++)
    {
      fastjet::PseudoJet &constituent = constituents[i];
      Float_t pt = constituent.pt();
      if (pt<0.15) continue;
      
      Int_t trackIndex = constituent.user_index();

      AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(trackIndex));

      Int_t TfPdgcode = particle->GetPdgCode();
      int fPdgIndex = -1;
      if ((TMath::Abs(TfPdgcode) == 211))
        fPdgIndex = 0; //Particle is a Pion
      else if ((TMath::Abs(TfPdgcode) == 321))
        fPdgIndex = 1; //Particle is a Kaon
      else if ((TMath::Abs(TfPdgcode) == 2212))
        fPdgIndex = 2; //Particle is a Proton

      Bool_t fSignMC = kFALSE;
      if (TfPdgcode > 0)
        fSignMC = kFALSE; //Particle is positive
      else
        fSignMC = kTRUE; //Particle is negative

      Float_t fPMC_Truth = particle->P();
      Float_t fPtMC_Truth = particle->Pt();
      Float_t fEtaMC_Truth = particle->Eta();
      Float_t fPhiMC_Truth = particle->Phi();

      if (found == true) {
        fastjet::PseudoJet jet_reco = jets_rec[matched_reco_jet_i];
        double matched_rec_pTsub = jet_reco.pt()-rho_MC_Rec*jet_reco.area();
        if (fPdgIndex == 0) fHist_JetMatch_MCTruth_TrackEff_Den_pi->Fill(fPtMC_Truth, matched_rec_pTsub, jetptsub);
        if (fPdgIndex == 1) fHist_JetMatch_MCTruth_TrackEff_Den_ka->Fill(fPtMC_Truth, matched_rec_pTsub, jetptsub);
        if (fPdgIndex == 2) fHist_JetMatch_MCTruth_TrackEff_Den_pr->Fill(fPtMC_Truth, matched_rec_pTsub, jetptsub);
      }

      if (jetptsub <= fjetMinPtSub || jetptsub >= fjetMaxPtSub) continue;

      if (fPdgIndex == 0) fHist_Jet_MCTruth_TrackEff_Den_pi->Fill(fPtMC_Truth, TMath::Abs(fEtaMC_Truth));
      if (fPdgIndex == 1) fHist_Jet_MCTruth_TrackEff_Den_ka->Fill(fPtMC_Truth, TMath::Abs(fEtaMC_Truth));
      if (fPdgIndex == 2) fHist_Jet_MCTruth_TrackEff_Den_pr->Fill(fPtMC_Truth, TMath::Abs(fEtaMC_Truth));

    }
  }

  if ( (fDoPerpCone || fDoRandCone) && (leadJetPtSub_Gen > fjetMinPtSub) ) {
    //now we find the eta,phi perp (and R.C.) to it
    //From the AliAnalysisTaskRhoPerpCone Task
    Double_t dPhi1 = 999.;
    Double_t dPhi2 = 999.;
    Double_t dEta = 999.;
    Double_t Axis1 = 999, Axis2 = 999;

    if (fDoRandCone){
      double OffsetRandom = (1 + fRandom->Rndm()) * TMath::Pi()/3; // creating random Phi in the range of PI/3 and 2PI/3
      Axis1 = leadJetPhi_Gen + OffsetRandom;  // adding the random Phi to the Phi of leading jet, to get the Phi of random cone
      Axis2 = leadJetPhi_Gen - OffsetRandom; // mirroring it, so you get Phi of the 2nd random cone
      if(Axis1 > TMath::TwoPi()) Axis1 -= TMath::TwoPi(); // checking if it's larger than 2PI (leading jet phi is distributed in range of 0 and 2PI
      if(Axis2 < 0) Axis2 += TMath::TwoPi(); // checking if 2nd cone is smaller than 0
    }

    if (fDoPerpCone){
      Axis1 = ((leadJetPhi_Gen + (TMath::Pi() / 2.)) > TMath::TwoPi()) ? leadJetPhi_Gen - ((3. / 2.) * TMath::Pi()) : leadJetPhi_Gen + (TMath::Pi() / 2.);
      Axis2 = ((leadJetPhi_Gen - (TMath::Pi() / 2.)) < 0) ? leadJetPhi_Gen + ((3. / 2.) * TMath::Pi()) : leadJetPhi_Gen - (TMath::Pi() / 2.);
    }

    TClonesArray* AODMCTrackArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if (AODMCTrackArray == NULL) return;
    if (fUseCouts) std::cout << "AODMCTrackArray->GetEntriesFast() is " << AODMCTrackArray->GetEntriesFast() << std::endl; //AOD Number of gen particles - used for looping over tracks

    // Loop over all primary MC particle
    for(Long_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {

      AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
      if (!particle) {
        if (fUseCouts) std::cout << "couldn't find particle" << std::endl;
        continue;
      }

      Float_t fEtaMC_Truth = particle->Eta();
      if (TMath::Abs(fEtaMC_Truth) > fEtaCut) continue;
                    
      Float_t mod_track_phi = particle->Phi() + TMath::Pi();
      //Check if the track is within the R=0.4 cone in eta, phi
      dPhi1 = TMath::Abs(mod_track_phi - Axis1);
      dPhi1 = (dPhi1 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi1 : dPhi1;
      dPhi2 = TMath::Abs(mod_track_phi - Axis2);
      dPhi2 = (dPhi2 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi2 : dPhi2;
      dEta = leadJetEta_Gen - particle->Eta();

      if ((TMath::Sqrt(dPhi1 * dPhi1 + dEta * dEta) > 0.4) && (TMath::Sqrt(dPhi2 * dPhi2 + dEta * dEta) > 0.4)) continue; //scale the yields by 1/(2Ncones*piR^2) offline

      if(!particle->IsPhysicalPrimary()) {
        if (fUseCouts) std::cout << "particle wasnt primary" << std::endl;
        continue;
      }

      Float_t fPMC_Truth = particle->P();
      Float_t fPtMC_Truth = particle->Pt();
      Float_t fPhiMC_Truth = particle->Phi();
      
      if (TMath::Abs(particle->Y()) > fYCut && fDoRapCut) {
        if (fUseCouts) std::cout << "particle didnt pass y cut" << std::endl;
        continue;
      }

      if (TMath::Abs(particle->Eta()) > fEtaCut) {
        if (fUseCouts) std::cout << "particle didn't pass eta cut" << std::endl;
        continue;
      }

      if (particle->Pt() < 0.15 || particle->Pt() > 100.0) {
        if (fUseCouts) std::cout << "particle didn't pass pT cut" << std::endl;
        continue;
      }

      Int_t TfPdgcode = particle->GetPdgCode();
      int fPdgIndex = -1;
      if ((TMath::Abs(TfPdgcode) == 211))
        fPdgIndex = 0; //Particle is a Pion
      else if ((TMath::Abs(TfPdgcode) == 321))
        fPdgIndex = 1; //Particle is a Kaon
      else if ((TMath::Abs(TfPdgcode) == 2212))
        fPdgIndex = 2; //Particle is a Proton

      Bool_t fSignMC = kFALSE;
      if (TfPdgcode > 0)
        fSignMC = kFALSE; //Particle is positive
      else
        fSignMC = kTRUE; //Particle is negative

      int iMCPartLabel = TMath::Abs(particle->GetLabel());
      bool particleIsPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iMCPartLabel,MCEvent());
      if (particleIsPileup) {
        if (fUseCouts) std::cout << "This MC particle was from OOB pileup" << std::endl;
        continue;
      }

      //Save truth level by PID
      if (fPdgIndex == 0) fHistMCTruth_TrackEff_Den_pi->Fill(fPtMC_Truth, TMath::Abs(fEtaMC_Truth));
      if (fPdgIndex == 1) fHistMCTruth_TrackEff_Den_ka->Fill(fPtMC_Truth, TMath::Abs(fEtaMC_Truth));
      if (fPdgIndex == 2) fHistMCTruth_TrackEff_Den_pr->Fill(fPtMC_Truth, TMath::Abs(fEtaMC_Truth));
    } //end gen particle loop for cones
    fFilledUECone_Gen = kTRUE;
  } //end gen cone loop  

  if (fUseCouts) cout << "num_reco_jets_w_multiple_matches is " << num_reco_jets_w_multiple_matches << endl;
  if (fUseCouts) cout << "num_reco_jets_w_matches is " << num_reco_jets_w_matches << endl;
  fall_reco_jets_w_multiple_matches += num_reco_jets_w_multiple_matches;
  fall_reco_jets_w_matches += num_reco_jets_w_matches;

}
//________________________________________________________________________

//________________________________________________________________________
void AliAnalysisTaskJetHadroAOD::FillEmbJets()
{
  //Parts from jet extractor task
  //
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FillEmbJets ===== " << std::endl;
  AliTOFPIDResponse fTOFPIDResponse = fPIDResponse->GetTOFResponse();

  //
  Float_t pT_sub_min = fjetMinPtSub;
  if (fUseCouts) cout << "Minimum jet pt after subtraction is " << fjetMinPtSub << endl;
  if (fUseCouts) cout << "Maximum jet pt after subtraction is " << fjetMaxPtSub << endl;

  Int_t trackOrigin = -10;

  Float_t leadJetPhi_Emb = 999;
  Float_t leadJetEta_Emb = 999;
  Float_t leadJetPtSub_Emb = -999;

  Float_t leadJetPhi_Rec = 999;
  Float_t leadJetEta_Rec = 999;
  Float_t leadJetPtSub_Rec = -999;

  Float_t leadJetPhi_Gen = 999;
  Float_t leadJetEta_Gen = 999;
  Float_t leadJetPtSub_Gen = -999;


  //Gen Signal JETS
  // Get the jet container
  fGenJetContainer = this->GetJetContainer("GenJets"); //2 = Pythia particle truth jets
  if (this->GetJetContainer("GenJets")) fDoPartLevelMatching = kTRUE;
  TString fGenRhoName = fGenJetContainer->GetRhoName();
  if (fUseCouts) cout << "Gen Rho Name is " << fGenRhoName << endl;

  if (fGenJetContainer->GetRhoParameter()) fjetGenRhoVal = fGenJetContainer->GetRhoVal();
  if (fUseCouts) cout << "In the FillEmbJets Gen Rho value is " << fjetGenRhoVal << endl;

  //RECO Signal JETS
  // Get the jet container
  fRecoJetContainer = this->GetJetContainer("RecoJets"); //1 = Pythia det level jets
  if (this->GetJetContainer("RecoJets")) fDoDetLevelMatching = kTRUE; 
  TString fRecoRhoName = fRecoJetContainer->GetRhoName();
  if (fUseCouts) cout << "Reco Rho Name is " << fRecoRhoName << endl;

  if (fRecoJetContainer->GetRhoParameter()) fjetRecoRhoVal = fRecoJetContainer->GetRhoVal();
  if (fUseCouts) cout << "In the FillEmbJets Reco Rho value is " << fjetRecoRhoVal << endl;


  //EMB Signal JETS
  // Get the jet container
  fEmbJetContainer = this->GetJetContainer("EmbJets"); //0 = embedded = Data + Pythia det level jets
  TString fEmbRhoName = fEmbJetContainer->GetRhoName();
  if (fUseCouts) cout << "Emb Rho Name is " << fEmbRhoName << endl;
  AliParticleContainer* embParticles = fEmbJetContainer->GetParticleContainer();
  if (embParticles) cout << "embParticles name is " << embParticles->GetName() << endl;


  if (fEmbJetContainer->GetRhoParameter()) fjetEmbRhoVal = fEmbJetContainer->GetRhoVal();
  if (fUseCouts) cout << "In the FillEmbJets Emb Rho value is " << fjetEmbRhoVal << endl;

  if(fDoDetLevelMatching) DoJetMatching(); //this matches by R

  fEmbJetContainer->ResetCurrentID();

  for(auto jet_emb : fEmbJetContainer->accepted())
  {

    fhasAcceptedEMCjet = 1;
    //if (jet_emb->Pt() < fTrackPt || jet_emb->Pt() > 1000.0 || TMath::Abs(jet_emb->Eta()) >= jetAbsEtaCut) continue; //this is not needed when using jetcontainer->accepted
    Float_t jet_embpt = jet_emb->Pt();
    Float_t jet_embphi = jet_emb->Phi();
    Float_t jet_embeta = jet_emb->Eta();
    Float_t jet_embArea = jet_emb->AreaPt();
    Float_t jet_embptsub = jet_emb->PtSub(fjetEmbRhoVal, kFALSE);

    Double_t matchedJetDistance_Det = -0.1;
    Double_t matchedJetPt_Det = -0.1;
    Double_t matchedJetEta_Det = -1.1;
    Double_t matchedJetPhi_Det = -0.1;
    Double_t matchedJetDistance_Part = -0.1;
    Double_t matchedJetPt_Part = -0.1;
    Double_t matchedJetEta_Part = -1.1;
    Double_t matchedJetPhi_Part = -0.1;
    Double_t truePtFraction = 0;
    Double_t truePtFraction_PartLevel = 0;

    // Get true estimators: for pt
    GetTrueJetPtFraction(jet_emb, truePtFraction, truePtFraction_PartLevel);
    
    fHist_TrueJetPtFraction->Fill(jet_embptsub, truePtFraction, truePtFraction_PartLevel);

    //HERE is where the pT matching condition comes in
    GetMatchedJetObservables(jet_emb, matchedJetPt_Det, matchedJetPt_Part, matchedJetPhi_Det, matchedJetEta_Det, matchedJetPhi_Part, matchedJetEta_Part, matchedJetDistance_Det, matchedJetDistance_Part);

    fHist_MatchJetPts->Fill(jet_embptsub, matchedJetPt_Det, matchedJetPt_Part);
    fHist_MatchJetEtas->Fill(jet_embeta, matchedJetEta_Det, matchedJetEta_Part);
    fHist_MatchJetDeltaRs->Fill(jet_embptsub, matchedJetDistance_Det, matchedJetDistance_Part);


    if (matchedJetPt_Det==-0.1 || matchedJetPt_Part==-0.1) {
      if (fUseCouts) cout << "THIS JET WASNT MATCHED" << endl;
      continue;
    }


    if (jet_embptsub > pT_sub_min && jet_embptsub < fjetMaxPtSub && jet_embArea > fjetMinArea)
    {
      fhasRealEMCjet = 1;
      fNumRealEMCJets +=1;
    }

    fHistJet_ptsub_v_area->Fill(jet_embArea, jet_embptsub);

    if (jet_embArea <= fjetMinArea) continue;

    fHistJet_kin->Fill(jet_embptsub, jet_embeta, jet_embphi);
    fHistJet_moms->Fill(jet_embptsub, jet_embpt);


    if (jet_embptsub > leadJetPtSub_Emb){
      leadJetPtSub_Emb = jet_embptsub;
      leadJetEta_Emb = jet_embeta;
      leadJetPhi_Emb = jet_embphi;
    }

    for(Int_t i = 0; i < jet_emb->GetNumberOfParticleConstituents(); i++)
    {
      Bool_t isPbPbPart = kFALSE;
      const AliVParticle* particle = jet_emb->Track(i);
      AliAODTrack* trackReal = (AliAODTrack*)(particle);
      if (trackReal==NULL) {
        if (fUseCouts) std::cout << "Didn't have a emb track" << std::endl;
        continue;
      }

      Int_t lab = TMath::Abs(trackReal->GetLabel());
      Int_t imc=trackReal->GetLabel();
      
      if (imc==-1) isPbPbPart = kTRUE; //A Pb-pb particle

      if(!trackReal || !trackReal->TestFilterBit(fAOD_FilterBits)) {
        continue;
      }

      //Track cuts start
      fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
      fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
      //
      //if (!fAODtrackCuts->AcceptTrack(trackReal))  continue;    // default cuts - redundant since these track cuts are passed to jet finder
      //
      // Get the track variables

      //
      //if (trackReal-> GetInnerParam()){ AODs don't have this info
        fPtotMC       = trackReal->P(); //trackReal->GetInnerParam()->GetP();
        fTPCSignalMC  = trackReal->GetTPCsignal();
      //}
      fEtaMC        = trackReal->Eta();
      Float_t fYMC          = trackReal->Y();
      fPtMC         = trackReal->Pt();
      fSignMC       = trackReal->GetSign();
      Float_t pMC   = trackReal->P();
      Float_t fPhiMC= trackReal->Phi();

      if (fPtMC<0.15 || fPtMC>100.0) continue; //So we can match the jets that we throw out w/ max track pT>100

      if (TMath::Abs(fEtaMC) > fEtaCut) continue;

      fHist_JERS_PbPbunid_tru->Fill( (jet_embptsub - matchedJetPt_Part)/matchedJetPt_Part, matchedJetPt_Part, fPtMC);
      fHist_JERS_PbPbunid_det->Fill( (jet_embptsub - matchedJetPt_Det)/matchedJetPt_Det, matchedJetPt_Det, fPtMC);

      fHist_JERS_PbPbunid_truhyb->Fill( (jet_embptsub - matchedJetPt_Part)/matchedJetPt_Part, jet_embptsub, fPtMC);
      fHist_JERS_PbPbunid_dethyb->Fill( (jet_embptsub - matchedJetPt_Det)/matchedJetPt_Det, jet_embptsub, fPtMC);


      if (isPbPbPart == kFALSE){

        if (fMCParticleArray == NULL) continue;
        const AliAODMCParticle *mc=static_cast<AliAODMCParticle*>(fMCParticleArray->At(lab));
        Int_t pdg = mc->GetPdgCode();

        if(!mc->IsPhysicalPrimary()) {
          if (fUseCouts) std::cout << "particle wasnt primary" << std::endl;
          continue;
        }

        Int_t TOFTrkLabel[3] = { -1 };

        Int_t fMCTOFMatch = -2;
        trackReal->GetTOFLabel(TOFTrkLabel);
        if (TOFTrkLabel[0] == -1) {// Track was not matched to any TOF hit.
          fMCTOFMatch = -1;
        }
        else if (lab == TOFTrkLabel[0]) {// Track was correctly matched to a TOF hit.
          fMCTOFMatch = 0;
        }
        else {// Track was matched to a TOF hit but comes from mismatch!
          fMCTOFMatch = 1;
        }
        
        // check the origin of the track
        trackOrigin = -10;

        //
        //
        // match the track with mc track
        Int_t iPart = -1; //-10;
        if (TMath::Abs(pdg) == kPDGel) { iPart = 0; } // select el
        if (TMath::Abs(pdg) == kPDGpi) { iPart = 1; } // select pi
        if (TMath::Abs(pdg) == kPDGka) { iPart = 2; } // select ka
        if (TMath::Abs(pdg) == kPDGpr) { iPart = 3; } // select pr
        if (TMath::Abs(pdg) == kPDGde) { iPart = 4; } // select de
        if (TMath::Abs(pdg) == kPDGmu) { iPart = 5; } // select mu

        Float_t fRapidityMC = TMath::ASinH(fPt / TMath::Sqrt(AliPID::ParticleMass(iPart + 1) * AliPID::ParticleMass(iPart + 1) + fPt * fPt) * TMath::SinH(fEta));

        if (TMath::Abs(fRapidityMC) > fYCut && fDoRapCut) {
          if (fUseCouts) std::cout << "didn't pass rap cut" << std::endl;
          continue;
        }

        fHist_JERS_truPythia->Fill( (jet_embptsub - matchedJetPt_Part)/matchedJetPt_Part, matchedJetPt_Part, fPtMC);
        fHist_JERS_detPythia->Fill( (jet_embptsub - matchedJetPt_Det)/matchedJetPt_Det, matchedJetPt_Det, fPtMC);

        fHist_JERS_truhybPythia->Fill( (jet_embptsub - matchedJetPt_Part)/matchedJetPt_Part, jet_embptsub, fPtMC);
        fHist_JERS_dethybPythia->Fill( (jet_embptsub - matchedJetPt_Det)/matchedJetPt_Det, jet_embptsub, fPtMC);


        if (iPart==1) fHist_JERS_truPythia_pi->Fill( (jet_embptsub - matchedJetPt_Part)/matchedJetPt_Part, matchedJetPt_Part, fPtMC);
        if (iPart==1) fHist_JERS_detPythia_pi->Fill( (jet_embptsub - matchedJetPt_Det)/matchedJetPt_Det, matchedJetPt_Det, fPtMC);

        if (iPart==1) fHist_JERS_truhybPythia_pi->Fill( (jet_embptsub - matchedJetPt_Part)/matchedJetPt_Part, jet_embptsub, fPtMC);
        if (iPart==1) fHist_JERS_dethybPythia_pi->Fill( (jet_embptsub - matchedJetPt_Det)/matchedJetPt_Det, jet_embptsub, fPtMC);

        if (iPart==2) fHist_JERS_truPythia_ka->Fill( (jet_embptsub - matchedJetPt_Part)/matchedJetPt_Part, matchedJetPt_Part, fPtMC);
        if (iPart==2) fHist_JERS_detPythia_ka->Fill( (jet_embptsub - matchedJetPt_Det)/matchedJetPt_Det, matchedJetPt_Det, fPtMC);

        if (iPart==2) fHist_JERS_truhybPythia_ka->Fill( (jet_embptsub - matchedJetPt_Part)/matchedJetPt_Part, jet_embptsub, fPtMC);
        if (iPart==2) fHist_JERS_dethybPythia_ka->Fill( (jet_embptsub - matchedJetPt_Det)/matchedJetPt_Det, jet_embptsub, fPtMC);

        if (iPart==3) fHist_JERS_truPythia_pr->Fill( (jet_embptsub - matchedJetPt_Part)/matchedJetPt_Part, matchedJetPt_Part, fPtMC);
        if (iPart==3) fHist_JERS_detPythia_pr->Fill( (jet_embptsub - matchedJetPt_Det)/matchedJetPt_Det, matchedJetPt_Det, fPtMC);

        if (iPart==3) fHist_JERS_truhybPythia_pr->Fill( (jet_embptsub - matchedJetPt_Part)/matchedJetPt_Part, jet_embptsub, fPtMC);
        if (iPart==3) fHist_JERS_dethybPythia_pr->Fill( (jet_embptsub - matchedJetPt_Det)/matchedJetPt_Det, jet_embptsub, fPtMC);
      }


      if (jet_embptsub <= fjetMinPtSub || jet_embptsub >= fjetMaxPtSub) continue; 

      if (isPbPbPart == kTRUE){

        GetExpecteds(trackReal);
        SetCutBitsAndSomeTrackVariables(trackReal);
        Float_t length    = trackReal->GetIntegratedLength();
        Float_t tofSignal = trackReal->GetTOFsignal();
        Float_t t0 = fTOFPIDResponse.GetStartTime(fPVertex);
        Float_t beta = -.05;
        //int nTOFClusters = trackReal->GetNTOFclusters(); //All matchable clusters AOD doesn't have this
        if((length > 0) && (tofSignal > 0)) beta = length / (2.99792458e-2*(tofSignal - t0));

        Float_t fTPCmom_choice = fPtot;
        if (fSetTPCmom == 1) fTPCmom_choice = fPVertex;
        if (fSetTPCmom == 2) fTPCmom_choice = fPt;

        Float_t fTOFmom_choice = fPt;
        if (fSetTOFmom == 1) fTPCmom_choice = fPVertex;
        if (fSetTOFmom == 2) fTPCmom_choice = fPtot;

        Float_t fBetamom_choice = fPVertex;
        if (fSetBetamom == 1) fTPCmom_choice = fPt;
        if (fSetBetamom == 2) fTPCmom_choice = fPtot;

        Float_t fEta_choice = fEta;
        if (fSetEta == 1) fEta_choice = fY;

        Float_t DCAxy = -999.;
        Float_t DCAz = -999.;
        trackReal->GetImpactParameters(DCAxy, DCAz);
        Float_t comb_sig_pi = TMath::Sqrt(fNSigmasPiTPC * fNSigmasPiTPC + fNSigmasPiTOF * fNSigmasPiTOF);
        Float_t comb_sig_pr = TMath::Sqrt(fNSigmasPrTPC * fNSigmasPrTPC + fNSigmasPrTOF * fNSigmasPrTOF);

        if (TMath::Abs(fEta) < fEtaCut){
          if (fFill_TPC) fHistJetTracks_dEdx->Fill(fTPCmom_choice,fTPCSignal,TMath::Abs(fEta_choice));
          if (fFillpTPC_pT) fHistJetTracks_moms->Fill(fPt,fPtot);
          if (fFillp_pT) fHistJetTracks_moms_p->Fill(fPt,fPVertex);
          if (fFillpTPC_p) fHistJetTracks_moms_pTPC_p->Fill(fPtot,fPVertex);
          fHistJetTracks_kin->Fill(fPt,TMath::Abs(fEta),fPhi);
          if (comb_sig_pi < 2.0) fHistJet_pi_DCAxy->Fill(fPt, DCAxy);
          if (comb_sig_pr < 2.0) fHistJet_pr_DCAxy->Fill(fPt, DCAxy);
        }

        //TPC-TOF Matching conditions: Use standard TPC tracks, then require kTIME and kTOFout
        Bool_t fTOFout = kFALSE;
        if ((trackReal->GetStatus() & AliAODTrack::kTOFout) != 0) { //Track has the kTOFout flag
        fTOFout = kTRUE;
        }

        //
        //kTIME flag
        Bool_t fTime = kFALSE;
        if ((trackReal->GetStatus() & AliAODTrack::kTIME) != 0) { //Track has the kTIME flag
        fTime = kTRUE;
        }

        Double_t inttime[5]; //6 is needed to account for earlier species - 0 = electron, 1 = muon, 2 = pion, 3 = kaon, 4 = proton, ESD only: 5 = deuteron
        trackReal->GetIntegratedTimes(inttime, 5); // Returns the array with integrated times for each particle hypothesis
        Float_t fTOFMismatchTime = -20000000.;
        fTOFMismatchTime = AliTOFPIDResponse::GetMismatchRandomValue(trackReal->Eta());
        if (fTOFMismatchTime <= 0){
          fTOFMismatchTime = -20000000.;
        }

        for (Int_t i = 0; i < 3; i++) {
          const Float_t beta_expec = length / (2.99792458e-2 * (inttime[i+2] - t0));
          if (i==0){
            fHistjet_BetaExpec_pi->Fill(fTPCmom_choice, beta_expec);
            fHist_jet_pi_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
            if (TMath::Abs(fNSigmasPiTOF) < 2.0){
              fHistjet_TOFSigmaExpec_pi->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasElTOF) < 2.0){
              fHist_jet_elExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasMuTOF) < 2.0){
              fHist_jet_muExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasKaTOF) < 2.0){
              fHist_jet_kaExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasPrTOF) < 2.0){
              fHist_jet_prExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
            }
          }
          if (i==1){
            fHistjet_BetaExpec_ka->Fill(fTPCmom_choice, beta_expec);
            fHist_jet_ka_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
            if (TMath::Abs(fNSigmasKaTOF) < 2.0){
              fHistjet_TOFSigmaExpec_ka->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasPiTOF) < 2.0){
              fHist_jet_piExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasPrTOF) < 2.0){
              fHist_jet_prExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
            }
          }
          if (i==2){
            fHistjet_BetaExpec_pr->Fill(fTPCmom_choice, beta_expec);
            fHist_jet_pr_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
            if (TMath::Abs(fNSigmasPrTOF) < 2.0){
              fHistjet_TOFSigmaExpec_pr->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasPiTOF) < 2.0){
              fHist_jet_piExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasKaTOF) < 2.0){
              fHist_jet_kaExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
            }
            if (TMath::Abs(fNSigmasDeTOF) < 2.0){
              fHist_jet_deExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
            }
          }
        }

        if (TMath::Abs(fEta) < fEtaCut  && fTOFout && fTime && fFill_TOF){
          fHistJetTracks_beta->Fill(fBetamom_choice,beta);
          fHistJetTracks_TOFpi_nsigma->Fill(fTOFmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
          fHistJetTracks_TOFka_nsigma->Fill(fTOFmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
          fHistJetTracks_TOFpr_nsigma->Fill(fTOFmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
          /*if (nTOFClusters < 2) { //ESD only
            fHistJetTracks_TOFpi_nsigma_1cls->Fill(fTOFmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
            fHistJetTracks_TOFka_nsigma_1cls->Fill(fTOFmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
            fHistJetTracks_TOFpr_nsigma_1cls->Fill(fTOFmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
          }*/
        }
      } //end isPbPbCondition
    } // ======= end of track loop for MC jet emb filling =======

  }//end emb jet loop

  if ( (fDoPerpCone || fDoRandCone) && (leadJetPtSub_Emb > fjetMinPtSub) ) {

    Bool_t isPbPbPart = kFALSE;

    //now we find the eta,phi perp (and R.C.) to it
    //From the AliAnalysisTaskRhoPerpCone Task
    Double_t dPhi1 = 999.;
    Double_t dPhi2 = 999.;
    Double_t dEta = 999.;
    Double_t Axis1 = 999, Axis2 = 999;

    if (fDoRandCone){
      double OffsetRandom = (1 + fRandom->Rndm()) * TMath::Pi()/3; // creating random Phi in the range of PI/3 and 2PI/3
      Axis1 = leadJetPhi_Emb + OffsetRandom;  // adding the random Phi to the Phi of leading jet, to get the Phi of random cone
      Axis2 = leadJetPhi_Emb - OffsetRandom; // mirroring it, so you get Phi of the 2nd random cone
      if(Axis1 > TMath::TwoPi()) Axis1 -= TMath::TwoPi(); // checking if it's larger than 2PI (leading jet phi is distributed in range of 0 and 2PI
      if(Axis2 < 0) Axis2 += TMath::TwoPi(); // checking if 2nd cone is smaller than 0
    }

    if (fDoPerpCone){
      Axis1 = ((leadJetPhi_Emb + (TMath::Pi() / 2.)) > TMath::TwoPi()) ? leadJetPhi_Emb - ((3. / 2.) * TMath::Pi()) : leadJetPhi_Emb + (TMath::Pi() / 2.);
      Axis2 = ((leadJetPhi_Emb - (TMath::Pi() / 2.)) < 0) ? leadJetPhi_Emb + ((3. / 2.) * TMath::Pi()) : leadJetPhi_Emb - (TMath::Pi() / 2.);
    }

    Int_t trackOrigin = -10; 

    AliParticleContainer * partCont = 0;
    TIter nextPartCont(&fParticleCollArray);
    while ((partCont = static_cast<AliParticleContainer*>(nextPartCont()))) {
      AliParticleIterableMomentumContainer itcont = partCont->accepted_momentum();
      for (AliParticleIterableMomentumContainer::iterator it = itcont.begin(); it != itcont.end(); it++) {
        AliVTrack * particle = static_cast<AliVTrack*>(it->second); //partIter.second);
      
        AliAODTrack* trackReal = (AliAODTrack*)(particle);
        
        if (trackReal==NULL) {
          if (fUseCouts) std::cout << "Didn't have a reco track" << std::endl;
          continue;
        }

        if(!trackReal || !trackReal->TestFilterBit(fAOD_FilterBits)) {
          continue;
        }

        Bool_t ifDefaultCuts = trackReal->TestFilterBit(fAOD_FilterBits);

        //
        SetCutBitsAndSomeTrackVariables(trackReal);
        fEtaMC        = trackReal->Eta();
        if (TMath::Abs(fEta) > fEtaCut) {
          continue;
        }
                      
        Float_t mod_track_phi = trackReal->Phi() + TMath::Pi();
        //Check if the track is within the R=0.4 cone in eta, phi
        dPhi1 = TMath::Abs(mod_track_phi - Axis1);
        dPhi1 = (dPhi1 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi1 : dPhi1;
        dPhi2 = TMath::Abs(mod_track_phi - Axis2);
        dPhi2 = (dPhi2 > TMath::Pi()) ? 2 * TMath::Pi() - dPhi2 : dPhi2;
        dEta = leadJetEta_Emb - trackReal->Eta();

        if ((TMath::Sqrt(dPhi1 * dPhi1 + dEta * dEta) > 0.4) && (TMath::Sqrt(dPhi2 * dPhi2 + dEta * dEta) > 0.4)) continue; //scale the yields by 1/(2Ncones*piR^2) offline
        
        //
        // --------------------------------------------------------------
        //   Fill the trees
        // --------------------------------------------------------------
        //

        fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
        fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
        // Get the track variables
        GetExpecteds(trackReal);

        Float_t length = trackReal->GetIntegratedLength();
        Float_t tofSignal = trackReal->GetTOFsignal();
        Float_t t0 = fTOFPIDResponse.GetStartTime(fPVertex);
        Float_t beta = -.05;
        //int nTOFClusters = trackReal->GetNTOFclusters(); //All matchable clusters ESD only
        if((length > 0) && (tofSignal > 0)) beta = length / (2.99792458e-2*(tofSignal - t0));

        if (fPt>100.0) continue; //So we can match the jets that we throw out w/ max track pT>100

        Float_t fTPCmom_choice = fPtot;
        if (fSetTPCmom == 1) fTPCmom_choice = fPVertex;
        if (fSetTPCmom == 2) fTPCmom_choice = fPt;

        Float_t fTOFmom_choice = fPt;
        if (fSetTOFmom == 1) fTPCmom_choice = fPVertex;
        if (fSetTOFmom == 2) fTPCmom_choice = fPtot;

        Float_t fBetamom_choice = fPVertex;
        if (fSetBetamom == 1) fTPCmom_choice = fPt;
        if (fSetBetamom == 2) fTPCmom_choice = fPtot;

        Float_t fEta_choice = fEta;
        if (fSetEta == 1) fEta_choice = fY;

        Float_t DCAxy = -999.;
        Float_t DCAz = -999.;
        trackReal->GetImpactParameters(DCAxy, DCAz);
        Float_t comb_sig_pi = TMath::Sqrt(fNSigmasPiTPC * fNSigmasPiTPC + fNSigmasPiTOF * fNSigmasPiTOF);
        Float_t comb_sig_pr = TMath::Sqrt(fNSigmasPrTPC * fNSigmasPrTPC + fNSigmasPrTOF * fNSigmasPrTOF);

        if (ifDefaultCuts == 1 && TMath::Abs(fEta) < fEtaCut){
          if (fFill_TPC) fHistIncTracks_dEdx->Fill(fTPCmom_choice,fTPCSignal,TMath::Abs(fEta_choice));
          if (fFillpTPC_pT) fHistIncTracks_moms->Fill(fPt,fPtot);
          if (fFillp_pT) fHistIncTracks_moms_p->Fill(fPt,fPVertex);
          if (fFillpTPC_p) fHistIncTracks_moms_pTPC_p->Fill(fPtot,fPVertex);
          fHistIncTracks_kin->Fill(fPt,TMath::Abs(fEta),fPhi);
          if (comb_sig_pi < 2.0) fHist_pi_DCAxy->Fill(fPt, DCAxy);
          if (comb_sig_pr < 2.0) fHist_pr_DCAxy->Fill(fPt, DCAxy);

          if (fFill_TPC_expecs){
            fHistIncTracks_mpi->Fill(fTPCmom_choice,fDEdxPi,TMath::Abs(fEta_choice));
            fHistIncTracks_spi->Fill(fTPCmom_choice,fSigmaPi,TMath::Abs(fEta_choice));

            fHistIncTracks_mel->Fill(fTPCmom_choice,fDEdxEl,TMath::Abs(fEta_choice));
            fHistIncTracks_sel->Fill(fTPCmom_choice,fSigmaEl,TMath::Abs(fEta_choice));

            fHistIncTracks_mka->Fill(fTPCmom_choice,fDEdxKa,TMath::Abs(fEta_choice));
            fHistIncTracks_ska->Fill(fTPCmom_choice,fSigmaKa,TMath::Abs(fEta_choice));

            fHistIncTracks_mpr->Fill(fTPCmom_choice,fDEdxPr,TMath::Abs(fEta_choice));
            fHistIncTracks_spr->Fill(fTPCmom_choice,fSigmaPr,TMath::Abs(fEta_choice));
          }
        }
        //
        //TPC-TOF Matching conditions: Use standard TPC tracks, then require kTIME and kTOFout
        Bool_t fTOFout = kFALSE;
        if ((trackReal->GetStatus() & AliAODTrack::kTOFout) != 0) { //Track has the kTOFout flag
        fTOFout = kTRUE;
        }

        //kTIME flag
        Bool_t fTime = kFALSE;
        if ((trackReal->GetStatus() & AliAODTrack::kTIME) != 0) { //Track has the kTIME flag
        fTime = kTRUE;
        }

        if (ifDefaultCuts == 1 && TMath::Abs(fEta) < fEtaCut && fTOFout && fTime && fFill_TOF_expecs){
          Double_t inttime[5]; //6 is needed to account for earlier species - 0 = electron, 1 = muon, 2 = pion, 3 = kaon, 4 = proton, ESD only 5 = deuteron
          trackReal->GetIntegratedTimes(inttime, 5); // Returns the array with integrated times for each particle hypothesis

          Float_t fTOFMismatchTime = -20000000.;
          fTOFMismatchTime = AliTOFPIDResponse::GetMismatchRandomValue(trackReal->Eta());
          if (fTOFMismatchTime <= 0){
            fTOFMismatchTime = -20000000.;
          }

          for (Int_t i = 0; i < 3; i++) {
            const Double_t beta_expec = length / (2.99792458e-2 * (inttime[i+2] - t0));
            if (i==0){
              fHistBetaExpec_pi->Fill(fTPCmom_choice, beta_expec);
              fHist_pi_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion), TMath::Abs(fEta_choice));
              if (TMath::Abs(fNSigmasPiTOF) < 2.0){
                fHistTOFSigmaExpec_pi->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion), TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasElTOF) < 2.0){
                fHist_elExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasMuTOF) < 2.0){
                fHist_muExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasKaTOF) < 2.0){
                fHist_kaExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasPrTOF) < 2.0){
                fHist_prExpec_pihyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kPion),TMath::Abs(fEta_choice));
              }
            }
            if (i==1){
              fHistBetaExpec_ka->Fill(fTPCmom_choice, beta_expec);
              fHist_ka_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
              if (TMath::Abs(fNSigmasKaTOF) < 2.0){
                fHistTOFSigmaExpec_ka->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasPiTOF) < 2.0){
                fHist_piExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasPrTOF) < 2.0){
                fHist_prExpec_kahyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kKaon),TMath::Abs(fEta_choice));
              }
            }
            if (i==2){
              fHistBetaExpec_pr->Fill(fTPCmom_choice, beta_expec);
              fHist_pr_mismatch->Fill(fTPCmom_choice, (fTOFMismatchTime - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
              if (TMath::Abs(fNSigmasPrTOF) < 2.0){
                fHistTOFSigmaExpec_pr->Fill(fTPCmom_choice, fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::ParticleMass(i+2)),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasPiTOF) < 2.0){
                fHist_piExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasKaTOF) < 2.0){
                fHist_kaExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
              }
              if (TMath::Abs(fNSigmasDeTOF) < 2.0){
                fHist_deExpec_prhyp->Fill(fTPCmom_choice, (tofSignal - fPIDResponse->GetTOFResponse().GetStartTime(fPVertex) - inttime[i+2])/fTOFPIDResponse.GetExpectedSigma(fPVertex, inttime[i+2], AliPID::kProton),TMath::Abs(fEta_choice));
              }
            }
          }
        }

        if (ifDefaultCuts == 1 && TMath::Abs(fEta) < fEtaCut && fTOFout && fTime && fFill_TOF){
          fHistIncTracks_beta->Fill(fBetamom_choice,beta);
          fHistIncTracks_t0->Fill(fBetamom_choice,t0);
          fHistIncTracks_TOFpi_nsigma->Fill(fTPCmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
          fHistIncTracks_TOFka_nsigma->Fill(fTPCmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
          fHistIncTracks_TOFpr_nsigma->Fill(fTPCmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
          /*if (nTOFClusters < 2) { //ESD only
            fHistIncTracks_TOFpi_nsigma_1cls->Fill(fTPCmom_choice,fNSigmasPiTOF,TMath::Abs(fEta_choice));
            fHistIncTracks_TOFka_nsigma_1cls->Fill(fTPCmom_choice,fNSigmasKaTOF,TMath::Abs(fEta_choice));
            fHistIncTracks_TOFpr_nsigma_1cls->Fill(fTPCmom_choice,fNSigmasPrTOF,TMath::Abs(fEta_choice));
          }*/
        }

        //
        // Get generated track info
        Int_t lab = TMath::Abs(trackReal->GetLabel());
        Int_t imc=trackReal->GetLabel();
        if (imc==-1) isPbPbPart = kTRUE;


        if (isPbPbPart==kFALSE){

          if (fMCParticleArray == NULL) continue;
  
          const AliAODMCParticle *mc=static_cast<AliAODMCParticle*>(fMCParticleArray->At(lab));
          Int_t pdg = mc->GetPdgCode();

          if(mc->IsPhysicalPrimary()) {
            trackOrigin = 0;
          }
          if(mc->IsSecondaryFromWeakDecay()) {
            trackOrigin = 1;
          }
          if(mc->IsSecondaryFromMaterial()) {
            trackOrigin = 2;
          }

          Float_t DCAxy = -999.;
          Float_t DCAz = -999.;

          trackReal->GetImpactParameters(DCAxy, DCAz);

          if (fUseCouts) cout << "DCAxy for trackOrigin= " << trackOrigin << " is " << DCAxy << endl;

          Int_t TOFTrkLabel[3] = { -1 }; //This can contain three particles wich occupy the same cluster

          Int_t fMCTOFMatch = -2;
          trackReal->GetTOFLabel(TOFTrkLabel);

          if (TOFTrkLabel[0] == -1) {// Track was not matched to any TOF hit.
            fMCTOFMatch = -1;
          }
          else if (lab == TOFTrkLabel[0]) {// Track was correctly matched to a TOF hit.
            fMCTOFMatch = 0;
          }
          else {// Track was matched to a TOF hit but comes from mismatch!
            fMCTOFMatch = 1;
          }
          
          // match the track with mc track
          Int_t iPart = -1;
          if (TMath::Abs(pdg) == kPDGel) { iPart = 0; } // select el
          if (TMath::Abs(pdg) == kPDGpi) { iPart = 1; } // select pi
          if (TMath::Abs(pdg) == kPDGka) { iPart = 2; } // select ka
          if (TMath::Abs(pdg) == kPDGpr) { iPart = 3; } // select pr
          if (TMath::Abs(pdg) == kPDGde) { iPart = 4; } // select de
          if (TMath::Abs(pdg) == kPDGmu) { iPart = 5; } // select mu

          //
          GetExpecteds(trackReal);
          SetCutBitsAndSomeTrackVariables(trackReal);
          //
          //if (trackReal-> GetInnerParam()){ AODs don't have this info
            fPtotMC       = trackReal->P(); //trackReal->GetInnerParam()->GetP();
            fTPCSignalMC  = trackReal->GetTPCsignal();
          //}
          fEtaMC        = trackReal->Eta();
          Float_t fYMC          = trackReal->Y();
          fPtMC         = trackReal->Pt();
          fSignMC       = trackReal->GetSign();
          Float_t pMC   = trackReal->P();
          Float_t fPhiMC= trackReal->Phi();

          Float_t fRapidityMC = TMath::ASinH(fPt / TMath::Sqrt(AliPID::ParticleMass(iPart + 1) * AliPID::ParticleMass(iPart + 1) + fPt * fPt) * TMath::SinH(fEta));

          if (TMath::Abs(fRapidityMC) > fYCut && fDoRapCut) {
            if (fUseCouts) std::cout << "particle didn't pass rapidity cut with MC PID" << std::endl;
            continue;
          }

          if (fPtMC<0.15 || fPtMC>100.0) continue; //So we can match the jets that we throw out w/ max track pT>100

          if (TMath::Abs(fEtaMC) > fEtaCut) continue;

          if(!mc->IsPhysicalPrimary()) {
            if (fUseCouts) std::cout << "particle wasnt primary" << std::endl; //most aren't primaries
            continue;
          }

          fHist_PC_spectra_Pythia->Fill(fPt);

          if (iPart==1) fHist_PC_spectra_Pythia_pi->Fill(fPt);
          if (iPart==2) fHist_PC_spectra_Pythia_ka->Fill(fPt);
          if (iPart==3) fHist_PC_spectra_Pythia_pr->Fill(fPt);

        }//end is PbPbCondition
      } // end inc tracks loop for cone
    }
    fFilledUECone_Emb = kTRUE;
  } //end cone loop


  for(auto jet_reco : fRecoJetContainer->accepted())
  {
    //if (jet_reco->Pt() < fTrackPt || jet_reco->Pt() > 1000.0 || TMath::Abs(jet_reco->Eta()) >= jetAbsEtaCut) continue; //this is not needed when using jetcontainer->accepted
    Float_t jet_recopt = jet_reco->Pt();
    Float_t jet_recophi = jet_reco->Phi();
    Float_t jet_recoeta = jet_reco->Eta();
    Float_t jet_recoArea = jet_reco->AreaPt();

    Float_t jet_recoptsub = jet_reco->PtSub(fjetRecoRhoVal, kFALSE);

    fHistMCR_Jet_ptsub_v_area->Fill(jet_recoArea, jet_recoptsub);

    if (jet_recoArea <= fjetMinArea) continue;

    for(Int_t i = 0; i < jet_reco->GetNumberOfParticleConstituents(); i++)
    {
      const AliVParticle* particle = jet_reco->Track(i);
      AliAODTrack* trackReal = (AliAODTrack*)(particle);
      if (trackReal==NULL) {
        if (fUseCouts) std::cout << "Didn't have a reco track" << std::endl;
        continue;
      }

      // Get generated track info
      Int_t lab = TMath::Abs(trackReal->GetLabel());
      Int_t imc=trackReal->GetLabel();

      if (fMCParticleArray == NULL) continue;
      const AliAODMCParticle *mc=static_cast<AliAODMCParticle*>(fMCParticleArray->At(lab));
      Int_t pdg = mc->GetPdgCode();

      //EMCal jet definition does not have this. Will in principle affect the jets selected
      if(!mc->IsPhysicalPrimary()) {
        if (fUseCouts) std::cout << "particle wasnt primary" << std::endl;
        continue;
      }

      Int_t TOFTrkLabel[3] = { -1 };

      Int_t fMCTOFMatch = -2;
      trackReal->GetTOFLabel(TOFTrkLabel);
      if (TOFTrkLabel[0] == -1) {// Track was not matched to any TOF hit.
        fMCTOFMatch = -1;
      }
      else if (lab == TOFTrkLabel[0]) {// Track was correctly matched to a TOF hit.
        fMCTOFMatch = 0;
      }
      else {// Track was matched to a TOF hit but comes from mismatch!
        fMCTOFMatch = 1;
      }
       
      // check the origin of the track
      trackOrigin = -10;

      if(!trackReal || !trackReal->TestFilterBit(fAOD_FilterBits)) {
        continue;
      }
      //
      //
      // match the track with mc track
      Int_t iPart = -1; //-10;
      if (TMath::Abs(pdg) == kPDGel) { iPart = 0; } // select el
      if (TMath::Abs(pdg) == kPDGpi) { iPart = 1; } // select pi
      if (TMath::Abs(pdg) == kPDGka) { iPart = 2; } // select ka
      if (TMath::Abs(pdg) == kPDGpr) { iPart = 3; } // select pr
      if (TMath::Abs(pdg) == kPDGde) { iPart = 4; } // select de
      if (TMath::Abs(pdg) == kPDGmu) { iPart = 5; } // select mu
      //
      GetExpecteds(trackReal);
      SetCutBitsAndSomeTrackVariables(trackReal);
      //
      //if (trackReal-> GetInnerParam()){ AODs don't have this info
        fPtotMC       = trackReal->P(); //trackReal->GetInnerParam()->GetP();
        fTPCSignalMC  = trackReal->GetTPCsignal();
      //}
      fEtaMC        = trackReal->Eta();
      Float_t fYMC          = trackReal->Y();
      fPtMC         = trackReal->Pt();
      fSignMC       = trackReal->GetSign();
      Float_t pMC   = trackReal->P();
      Float_t fPhiMC= trackReal->Phi();

      Float_t fRapidityMC = TMath::ASinH(fPt / TMath::Sqrt(AliPID::ParticleMass(iPart + 1) * AliPID::ParticleMass(iPart + 1) + fPt * fPt) * TMath::SinH(fEta));

      if (TMath::Abs(fRapidityMC) > fYCut && fDoRapCut) {
        if (fUseCouts) std::cout << "didn't pass rap cut" << std::endl;
        continue;
      }
      //
      //
      //

      if (fPtMC<0.15 || fPtMC>100.0) continue; //So we can match the jets that we throw out w/ max track pT>100

      if (TMath::Abs(fEtaMC) > fEtaCut) continue;

      fHist_det_spectra_Pythia->Fill(fPtMC, jet_recopt);
      if (iPart==1) fHist_det_spectra_Pythia_pi->Fill(fPtMC, jet_recopt);
      if (iPart==2) fHist_det_spectra_Pythia_ka->Fill(fPtMC, jet_recopt);
      if (iPart==3) fHist_det_spectra_Pythia_pr->Fill(fPtMC, jet_recopt);

    } // ======= end of track loop for MC jet reco filling =======

  }//end reco jet loop

  //Gen jets
  for(auto jet : fGenJetContainer->accepted())
  {
    //if (jet->Pt() < fTrackPt || jet->Pt() > 1000.0 || TMath::Abs(jet->Eta()) >= jetAbsEtaCut) continue; //this is not needed when using jetcontainer->accepted
    Float_t jetpt = jet->Pt();

    Float_t jetphi = jet->Phi();
    Float_t jeteta = jet->Eta();
    Float_t jetArea = jet->AreaPt();

    Float_t jetptsub = jet->PtSub(fjetGenRhoVal, kFALSE);

    fHistMCT_Jet_ptsub_v_area->Fill(jetArea, jetptsub);

    if (jetArea <= fjetMinArea) continue;

    for(Int_t i = 0; i < jet->GetNumberOfParticleConstituents(); i++)
    {
      const PWG::JETFW::AliEmcalParticleJetConstituent* prepart1 = jet->ParticleConstituentAt(i); 
      const AliVParticle* prepart2 = prepart1->GetParticle();
      AliAODMCParticle* particle = (AliAODMCParticle*)(prepart2);
      if (particle==NULL) {
        if (fUseCouts) std::cout << "Didn't have a gen track" << std::endl;
        continue;
      }

      int iMCPartLabel = TMath::Abs(particle->GetLabel());
      /*bool particleIsPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(iMCPartLabel,MCEvent());
      if (particleIsPileup) {
        if (fUseCouts) std::cout << "This jet MC particle was from OOB pileup" << std::endl;
        continue;
      }*/ //jet MC doesn't have pileup

      if (TMath::Abs(particle->Y()) > fYCut && fDoRapCut) {
        if (fUseCouts) std::cout << "particle didnt pass y cut" << std::endl;
        continue;
      }

      Int_t TfPdgcode = particle->GetPdgCode();
      int fPdgIndex = -1;
      if ((TMath::Abs(TfPdgcode) == 211))
        fPdgIndex = 0; //Particle is a Pion
      else if ((TMath::Abs(TfPdgcode) == 321))
        fPdgIndex = 1; //Particle is a Kaon
      else if ((TMath::Abs(TfPdgcode) == 2212))
        fPdgIndex = 2; //Particle is a Proton

      Bool_t fSignMC = kFALSE;
      if (TfPdgcode > 0)
        fSignMC = kFALSE; //Particle is positive
      else
        fSignMC = kTRUE; //Particle is negative

      Float_t fPMC_Truth = particle->P();
      Float_t fPtMC_Truth = particle->Pt();
      Float_t fEtaMC_Truth = particle->Eta();
      Float_t fPhiMC_Truth = particle->Phi();


      fHist_part_spectra_Pythia->Fill(fPMC_Truth, jetpt);
      if (fPdgIndex==0) fHist_part_spectra_Pythia_pi->Fill(fPMC_Truth, jetpt);
      if (fPdgIndex==1) fHist_part_spectra_Pythia_ka->Fill(fPMC_Truth, jetpt);
      if (fPdgIndex==2) fHist_part_spectra_Pythia_pr->Fill(fPMC_Truth, jetpt);

    }
  }//end gen jets loop

}
//________________________________________________________________________
void AliAnalysisTaskJetHadroAOD::FillEventTree()
{
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FillEventTree ===== " << std::endl;
  (*fTreeSRedirector)<<"jeteventInfo"<<
  "rhoFJ="      << frhoFJ << //event rho
  "rhoEMC="  << fjetRhoVal <<
  "rhoEMCReco="  << fjetRecoRhoVal <<
  "rhoEMCGen="  << fjetGenRhoVal <<
  "rhoEMCEmb=" << fjetEmbRhoVal <<
  "cent="      << fCentrality  <<  //  centrality
  "isGoodIncEvent="   << fisGoodIncEvent <<
  "filledUECone_Emb="   << fFilledUECone_Emb <<
  "filledUECone_Gen="   << fFilledUECone_Gen <<
  "filledUECone_Rec="   << fFilledUECone_Rec <<
  "hasAcceptedFJjet="   << fhasAcceptedFJjet <<
  "hasRealFJjet="   << fhasRealFJjet <<
  "hasAcceptedEMCjet="   << fhasAcceptedEMCjet <<
  "hasRealEMCjet="   << fhasRealEMCjet <<
  "NumRealFJJets="   << fNumRealFJJets <<
  "NumRealEMCJets="   << fNumRealEMCJets <<
  "fall_reco_jets_w_multiple_matches="   << fall_reco_jets_w_multiple_matches <<
  "fall_reco_jets_w_matches="   << fall_reco_jets_w_matches <<
  "\n";
}
//________________________________________________________________________
void AliAnalysisTaskJetHadroAOD::FillEventHistos()
{
  if (fUseCouts) std::cout << " Info::siweyhmi: ===== In the FillEventHistos ===== " << std::endl;

  //Good event info (inclusive and jets)
  fHist_Is_Good_Inc_Event->Fill(fisGoodIncEvent+0.5);
  if (fhasAcceptedFJjet){
    if (fhasRealFJjet){
      fHist_Has_Acc_Real_FJ_Jet->Fill(2.5);
    }
    else{
      fHist_Has_Acc_Real_FJ_Jet->Fill(1.5);
    }
  }
  else {
    fHist_Has_Acc_Real_FJ_Jet->Fill(0.5);
  }
  if (fhasAcceptedEMCjet){
    if (fhasRealEMCjet){
      fHist_Has_Acc_Real_EMC_Jet->Fill(2.5);
    }
    else{
      fHist_Has_Acc_Real_EMC_Jet->Fill(1.5);
    }
  }
  else {
    fHist_Has_Acc_Real_EMC_Jet->Fill(0.5);
  }

  //Number of Jets Info
  fHist_Num_Real_FJ_Jets->Fill(fNumRealFJJets+0.5);
  fHist_Num_Real_EMC_Jets->Fill(fNumRealEMCJets+0.5);

  //Rho Info
  fHist_RhoFJ->Fill(frhoFJ);
  fHist_RhoEMC->Fill(fjetRhoVal);
  fHist_RhoEMCReco->Fill(fjetRecoRhoVal);
  fHist_RhoEMCGen->Fill(fjetGenRhoVal);
  fHist_RhoEMCEmb->Fill(fjetEmbRhoVal);

  //UE Cone Info
  fHist_Filled_UEC_Emb->Fill(fFilledUECone_Emb+0.5);
  fHist_Filled_UEC_Gen->Fill(fFilledUECone_Gen+0.5);
  fHist_Filled_UEC_Rec->Fill(fFilledUECone_Rec+0.5);
}
//________________________________________________________________________
void AliAnalysisTaskJetHadroAOD::GetExpecteds(AliAODTrack* track)
{

  //
  // bettaGamma is not well deifned below bg=0.01 --> below 200MeV protons and deuterons
  Float_t ptotForBetaGamma = track->P(); //track->GetInnerParam()->GetP(); AODs don't have this info
  Float_t ptotForBetaGammaThr = 0.2;
  // if (ptotForBetaGamma<ptotForBetaGammaThr) return;
  //
  // --------------------------------------------------------------
  //  Calculates expected sigma and dEdx for a given track and returns colesest expected particle and its index
  // --------------------------------------------------------------
  //
  fNSigmasElTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
  fNSigmasPiTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  fNSigmasKaTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
  //
  fNSigmasElTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron, fPIDResponse->GetTOFResponse().GetTimeZero());
  fNSigmasMuTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kMuon, fPIDResponse->GetTOFResponse().GetTimeZero());
  fNSigmasPiTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion,     fPIDResponse->GetTOFResponse().GetTimeZero());
  fNSigmasKaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon,     fPIDResponse->GetTOFResponse().GetTimeZero());
  //
  if (ptotForBetaGamma>ptotForBetaGammaThr){
    fNSigmasPrTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    fNSigmasDeTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
    fNSigmasPrTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton,   fPIDResponse->GetTOFResponse().GetTimeZero());
    fNSigmasDeTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron, fPIDResponse->GetTOFResponse().GetTimeZero());
  }

  //
  //
  //
  // Electron Expected mean and sigma within 3nsigmaTPC
  if (TMath::Abs(fNSigmasElTPC)<3 && ptotForBetaGamma>ptotForBetaGammaThr) {
    fDEdxEl  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kElectron);
    fSigmaEl = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kElectron);
  }
  //
  // Pion Expected mean and sigma within 3nsigmaTPC
  if (TMath::Abs(fNSigmasPiTPC)<3 && ptotForBetaGamma>ptotForBetaGammaThr) {
    fDEdxPi  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kPion);
    fSigmaPi = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kPion);
  }
  //
  // Kaon Expected mean and sigma within 3nsigmaTPC
  if (TMath::Abs(fNSigmasKaTPC)<3 && ptotForBetaGamma>ptotForBetaGammaThr) {
    fDEdxKa  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kKaon);
    fSigmaKa = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kKaon);
  }
  //
  // Proton Expected mean and sigma within 3nsigmaTPC
  if (TMath::Abs(fNSigmasPrTPC)<3 && ptotForBetaGamma>ptotForBetaGammaThr) {
    fDEdxPr  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kProton);
    fSigmaPr = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kProton);
  }
  //
  // Deuteron Expected mean and sigma within 3nsigmaTPC
  if (TMath::Abs(fNSigmasDeTPC)<3 && ptotForBetaGamma>ptotForBetaGammaThr) {
    fDEdxDe  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kDeuteron);
    fSigmaDe = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kDeuteron);
  }

}
//________________________________________________________________________
void AliAnalysisTaskJetHadroAOD::SetCutBitsAndSomeTrackVariables(AliAODTrack* track)
{
  //
  // Set some track variables
  //
  fPVertex = track->P();
  fSign= track->GetSign();
  fPt  =track->Pt();
  fY   =track->Y();
  fPhi =track->Phi()-TMath::Pi();
  fEta = track->Eta();

  //
  // TPC related quantities
  //if (track->GetInnerParam()){ AODs don't have this info
    fPtot      = track->P(); //track->GetInnerParam()->GetP();
    fTPCSignal = track->GetTPCsignal();
    //
    Float_t ptotForBetaGamma = track->P(); //track->GetInnerParam()->GetP(); AODs don't have this info
    Float_t ptotForBetaGammaThr = 0.2;
    //
    // --------------------------------------------------------------
    //      Bayesian PID part
    // --------------------------------------------------------------
    //
    /*if (ptotForBetaGamma>ptotForBetaGammaThr) {
      fPIDCombined->SetDefaultTPCPriors();
      Double_t probTPC[AliPID::kSPECIES]={0.};
      Double_t probTOF[AliPID::kSPECIES]={0.};
      // Get TPC probabilities
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
      fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPC);
      fTrackProbPiTPC = probTPC[AliPID::kPion];
      fTrackProbKaTPC = probTPC[AliPID::kKaon];
      fTrackProbPrTPC = probTPC[AliPID::kProton];
      //fTrackProbDeTPC = probTPC[AliPID::kDeuteron];
      // Get TOF probabilities
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
      fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTOF);
      fTrackProbPiTOF = probTOF[AliPID::kPion];
      fTrackProbKaTOF = probTOF[AliPID::kKaon];
      fTrackProbPrTOF = probTOF[AliPID::kProton];
    }*/
  //}
  //

}
//________________________________________________________________________
Bool_t AliAnalysisTaskJetHadroAOD::CountEmptyEvents()
{
  //
  // count Empty Events
  //
  Bool_t emptyEvent= kTRUE;
  /*Int_t iTracks(fAOD->GetNumberOfTracks());           // see how many tracks there are in the event
  for(Int_t i(0); i < iTracks; i++) {                 // loop ove rall these tracks
      AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
      if(!track || !track->TestFilterBit(1)) continue;                            // if we failed, skip this track
      cout << "GOOD!" << endl;                  // plot the pt value of the track in a histogram
  }*/

  for (Int_t itrack=0;itrack<fAOD->GetNumberOfTracks();itrack++) {   // Track loop
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(itrack));
    if(!track || !track->TestFilterBit(fAOD_FilterBits)) continue;
    if (track->Pt()<100.0) { //aod no track cuts, need filterbit
      emptyEvent= kFALSE;
      break;
    }
  }
  //
  // check if the event is empty
  if (emptyEvent) {
    std::cout << " Info::siweyhmi: Empty event " << std::endl;
  }
  else {
    if (fUseCouts) std::cout << " Info::siweyhmi: ====== EVENT IS COOL GO AHEAD ======= " << std::endl;
  }
  return emptyEvent;

}
//________________________________________________________________________
void AliAnalysisTaskJetHadroAOD::DoJetMatching(){
  //From jet extractor task

   // Perform the matching before the main jet loop (Only if we use the tagger for matching)                                              
  AliJetContainer * jetsHybrid = this->GetJetContainer("EmbJets");
  AliJetContainer * jetsDetLevel = this->GetJetContainer("RecoJets");
  AliJetContainer * jetsPartLevel = this->GetJetContainer("GenJets");

  
  // Now, begin the actual matching.
  // Hybrid <-> det first
  AliDebugStream(2) << "Matching hybrid to detector level jets.\n";
  // First, we reset the tagging
    for(auto j : jetsHybrid->all()){
    j->ResetMatching();
  }
  for(auto j : jetsDetLevel->all()){
    j->ResetMatching();
  }
  // Next, we perform the matching
  PerformGeometricalJetMatching(*jetsHybrid, *jetsDetLevel, fJetMatchingRadius);
  // Now, begin the next matching stage
  // det <-> particle
  AliDebugStream(2) << "Matching detector level to particle level jets.\n";
  // First, we reset the tagging. We need to reset the det matching again to ensure
  // that it doesn't accidentally keep some latent matches to the hybrid jets.
  for(auto j : jetsDetLevel->all()){
    j->ResetMatching();
  }
  // if we do not need to do particle level matching return
  if(!fDoPartLevelMatching) return;
  
  for(auto j : jetsPartLevel->all()){
    j->ResetMatching();
  }
  // Next, we perform the matching
  PerformGeometricalJetMatching(*jetsHybrid, *jetsDetLevel, fJetMatchingRadius);
  // Now, begin the next matching stage 
  // det <-> particle
  AliDebugStream(2) << "Matching detector level to particle level jets.\n";
  // First, we reset the tagging. We need to reset the det matching again to ensure
  // that it doesn't accidentally keep some latent matches to the hybrid jets.
  for(auto j : jetsDetLevel->all()){
    j->ResetMatching();
  }
  for(auto j : jetsPartLevel->all()){
    j->ResetMatching();
  }
  // Next, we perform the matching
  PerformGeometricalJetMatching(*jetsDetLevel, *jetsPartLevel, fJetMatchingRadius);
}
//________________________________________________________________________
bool AliAnalysisTaskJetHadroAOD::PerformGeometricalJetMatching(AliJetContainer& contBase,
                                    AliJetContainer& contTag, double maxDist) 
{

  //from jet extractor task
  // Note that this function is also utilized in /PWGJE/EMCALJetTasks/UserTasks/AliAnalysisTaskEmcalJetHPerformance.cxx. For more details, see this file.
  // Setup
  const Int_t kNacceptedBase = contBase.GetNAcceptedJets(), kNacceptedTag = contTag.GetNAcceptedJets();
  if (!(kNacceptedBase && kNacceptedTag)) {
    return false;
  }

  // Build up vectors of jet pointers to use when assigning the closest jets.
  // The storages are needed later for applying the tagging, in order to avoid multiple occurrence of jet selection
  std::vector<AliEmcalJet*> jetsBase(kNacceptedBase), jetsTag(kNacceptedTag);

  int countBase(0), countTag(0);
  for (auto jb : contBase.accepted()) {
    jetsBase[countBase] = jb;
    countBase++;
  }
  for (auto jt : contTag.accepted()) {
    jetsTag[countTag] = jt;
    countTag++;
  }

  TArrayI faMatchIndexTag(kNacceptedBase), faMatchIndexBase(kNacceptedTag);
  faMatchIndexBase.Reset(-1);
  faMatchIndexTag.Reset(-1);

  // find the closest distance to the base jet
  countBase = 0;
  for (auto jet1 : contBase.accepted()) {
    double distance = maxDist;

    // Loop over all accepted jets and brute force search for the closest jet.
    // NOTE: current_index() returns the jet index in the underlying array, not
    //       the index within the accepted jets that are returned.
    int contTagAcceptedIndex = 0;
    for (auto jet2 : contTag.accepted()) {
      double dR = jet1->DeltaR(jet2);
      if (dR < distance && dR < maxDist) {
        faMatchIndexTag[countBase] = contTagAcceptedIndex;
        distance = dR;
      }
      contTagAcceptedIndex++;
    }


    countBase++;
  }

  // other way around
  countTag = 0;
  for (auto jet1 : contTag.accepted()) {
    double distance = maxDist;

    // Loop over all accepted jets and brute force search for the closest jet.
    // NOTE: current_index() returns the jet index in the underlying array, not
    //       the index within the accepted jets that are returned.
    int contBaseAcceptedIndex = 0;
    for (auto jet2 : contBase.accepted()) {
      double dR = jet1->DeltaR(jet2);
      if (dR < distance && dR < maxDist) {
        faMatchIndexBase[countTag] = contBaseAcceptedIndex;
        distance = dR;
      }
      contBaseAcceptedIndex++;
    }
    countTag++;
  }

  // check for "true" correlations
  // these are pairs where the base jet is the closest to the tag jet and vice versa
  // As the lists are linear a loop over the outer base jet is sufficient.
  AliDebugStream(1) << "Starting true jet loop: nbase(" << kNacceptedBase << "), ntag(" << kNacceptedTag << ")\n";
  for (int ibase = 0; ibase < kNacceptedBase; ibase++) {
    AliDebugStream(2) << "base jet " << ibase << ": match index in tag jet container " << faMatchIndexTag[ibase]
             << "\n";
    if (faMatchIndexTag[ibase] > -1) {
      AliDebugStream(2) << "tag jet " << faMatchIndexTag[ibase] << ": matched base jet " << faMatchIndexBase[faMatchIndexTag[ibase]] << "\n";
    }
    // We have a true correlation where each jet points to the other.
    if (faMatchIndexTag[ibase] > -1 && faMatchIndexBase[faMatchIndexTag[ibase]] == ibase) {
      AliDebugStream(2) << "found a true match \n";
      AliEmcalJet *jetBase = jetsBase[ibase], *jetTag = jetsTag[faMatchIndexTag[ibase]];
      // We have a valid pair of matched jets, so set the closest jet properties.
      if (jetBase && jetTag) {
        Double_t dR = jetBase->DeltaR(jetTag);
        jetBase->SetClosestJet(jetTag, dR);
        jetTag->SetClosestJet(jetBase, dR);
      }
    }
  }
  return true;
}
//________________________________________________________________________
void AliAnalysisTaskJetHadroAOD::GetTrueJetPtFraction(AliEmcalJet* jet, Double_t& truePtFraction, Double_t& truePtFraction_mcparticles)
{
  //from jet extractor task
  // #################################################################################
  // ##### FRACTION OF TRUE PT IN JET: Defined as "not from toy"

  Double_t pt_truth = 0.;
  Double_t pt_truth_mcparticles = 0.;
  Double_t pt_all   = 0.;
  truePtFraction = 0;
  truePtFraction_mcparticles = 0;

  // ### Loop over all tracks constituents
  for(Int_t iConst = 0; iConst < jet->GetNumberOfParticleConstituents(); iConst++)
  {
    const AliVParticle* particle = jet->GetParticleConstituents()[iConst].GetParticle();
    if(!particle) continue;
    if(particle->Pt() < 1e-6) continue;

    // Particles marked w/ labels within label range OR explicitly set as embedded tracks are considered to be from truth
    if (  (fIsEmbeddedEvent && jet->GetParticleConstituents()[iConst].IsFromEmbeddedEvent()) ||
          (!fIsEmbeddedEvent && ((particle->GetLabel() >= fTruthMinLabel) && (particle->GetLabel() < fTruthMaxLabel)))  ){
            pt_truth += particle->Pt();
    }
    pt_all += particle->Pt();
  }

  // ### Loop over all primary (charged) MC particles and check if they have a corresponding track
  //     Correspondence is checked geometrically, sum of matched particles pT is truth
  Double_t jetRadius = 0.4;
  if(fMCParticleArray)
    for(Int_t iPart=0; iPart<fMCParticleArray->GetEntriesFast();iPart++)
    {
      AliAODMCParticle* part = (AliAODMCParticle*)fMCParticleArray->At(iPart);
      if(!part) continue;
      if(!part->IsPhysicalPrimary()) continue;
      if(!part->Charge()) continue;
      if(part->Pt() < 1e-6) continue;

      Double_t deltaR = GetDistance(part->Eta(), jet->Eta(), part->Phi(), jet->Phi());
      if(deltaR <= jetRadius){
        pt_truth_mcparticles += part->Pt();
      }
    }

  if(pt_all)
  {
    truePtFraction = (pt_truth/pt_all);
    truePtFraction_mcparticles = (pt_truth_mcparticles/pt_all);
  }
}
//________________________________________________________________________
void AliAnalysisTaskJetHadroAOD::GetMatchedJetObservables(AliEmcalJet* jet, Double_t& detJetPt, Double_t& partJetPt, Double_t& detJetPhi, Double_t& detJetEta, Double_t& partJetPhi, Double_t& partJetEta, Double_t& detJetDistance, Double_t& partJetDistance)
{
  //from jet extractor task
  // Get the Matched Observables                                                                                                                   
  AliJetContainer * hybridJetCont = this->GetJetContainer("EmbJets");
  AliJetContainer * detJetCont    = this->GetJetContainer("RecoJets");
  AliJetContainer * partJetCont   = this->GetJetContainer("GenJets");
  
  // hybrid to detector level matching 
  AliEmcalJet * jet2 = jet->ClosestJet();
  if (!jet2) { return;}// if there is no match return.
  Double_t ptJet2 = jet2->Pt() - detJetCont->GetRhoVal() * jet2->Area();
  // This will retrieve the fraction of jet2's momentum in jet1.
  Double_t fraction = hybridJetCont->GetFractionSharedPt(jet); //jet is emb jet
  if (fraction < fJetMatchingSharedPtFraction) { return; }

  // if we are not doing the particle level match, fill all observables here and return
  if(!fDoPartLevelMatching){
    detJetPt          = ptJet2;
    detJetDistance    = jet->DeltaR(jet2);
    return; 
  }
  
  // detector to particle matching
  AliEmcalJet * jet3 = jet2->ClosestJet();
  if(!jet3){return;}
  Double_t ptJet3 = jet3->Pt() - partJetCont->GetRhoVal() * jet3->Area(); 
  // make no fraction cut on the detector to particle

  // if we are doing particle level matching, require that there is both a particle and detector level match
  detJetPt  = ptJet2;
  partJetPt = ptJet3;
  detJetPhi = jet2->Phi();
  partJetPhi = jet3->Phi();
  detJetEta = jet2->Eta();
  partJetEta = jet3->Eta();
  detJetDistance  = jet->DeltaR(jet2);
  partJetDistance = jet2->DeltaR(jet3);

  return;
}
//________________________________________________________________________
void AliAnalysisTaskJetHadroAOD::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  std::cout << " Info::siweyhmi: ===== In the Terminate ===== " << std::endl;

  /*
  AliTPCPIDResponse fTPCresponse = fPIDResponse->GetTPCResponse();
  Double_t temp_mip = fTPCresponse.GetMIP();
  std::cout << "temp_mip is " << temp_mip << std::endl;
  for (int j=0; j<8; j++){
    AliPID::EParticleType species = AliPID::kPion;
    if (j==1) {
      species = AliPID::kElectron;
      std::cout << "el" << std::endl;
    }
    if (j==2) {
      species = AliPID::kKaon;
      std::cout << "ka" << std::endl;
    }
    if (j==3) {
      species = AliPID::kProton;
      std::cout << "pr" << std::endl;
    }
    if (j==4) {
      species = AliPID::kDeuteron;
      std::cout << "d" << std::endl;
    }
    if (j==5) {
      species = AliPID::kTriton;
      std::cout << "t" << std::endl;
    }
    if (j==6) {
      species = AliPID::kHe3;
      std::cout << "He3" << std::endl;
    }
    if (j==7) {
      species = AliPID::kAlpha;
      std::cout << "He4" << std::endl;
    }
    TSpline3* resp_fxn = fTPCresponse.GetResponseFunction(species, AliTPCPIDResponse::ETPCgainScenario::kDefault);
    Double_t mass=AliPID::ParticleMassZ(species);
    std::cout << "The mass is " << AliPID::ParticleMassZ(species) << std::endl;
    //charge factor. BB goes with z^2, however in reality it is slightly larger (calibration, threshold effects, ...)
    // !!! Splines for light nuclei need to be normalised to this factor !!!
    const Double_t chargeFactor = TMath::Power(AliPID::ParticleCharge(species),2.3);
    std::cout << "The particle charge is " << AliPID::ParticleCharge(species) << std::endl;

    stringstream name;
    name << "resp_fxn_" << j;
    //gSystem->ChangeDirectory("/Users/sierraweyhmiller/Sierra/Caines_Research_Yr2/JetHadro/most_recent_code");
    //resp_fxn->Write();
    resp_fxn->SaveAs(name.str().c_str());
  }
  */

}
